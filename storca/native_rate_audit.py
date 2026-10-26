"""Reference-free native evidence checks for isolated known-reaction rates.

This is an internal-consistency audit, not a mechanism, uncertainty or physical
accuracy certificate. No experimental answer, extra temperature or fitted
multiplier is read. Electronic-structure programs are never dispatched here.
The native replay must run in the explicitly selected RMG interpreter.
"""
from __future__ import annotations

import argparse
from collections import Counter
import importlib
import json
import math
from pathlib import Path
import re
import sys
import time

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from storca.composite_rate_recipe import (
    NUMBER, _input, _native, bind_artifact, reopen, ArtifactBinding,
    reopen_composite_recipe_assessment,
)
from storca.arkane_rotor import artifact, read_document, validate_contract
from src.orca_runner import run_bounded_kinetics_command

AVOGADRO = 6.02214076e23
HARTREE_TO_J_MOL = 2625499.6394799
_NUCLEAR_CHARGE = {"H": 1, "C": 6, "O": 8}
_CHEMICAL_CORE = {"H": 0, "C": 2, "O": 2}
_QZVPP_FUNCTIONS = {"H": 30, "C": 57, "O": 57}


class UnsupportedRateModelError(ValueError):
    """The requested physical model is outside this audit's admitted scope."""


class PhysicalRateValidationError(ValueError):
    """Loaded evidence violates a required physical/numerical rate identity."""


def _finite(value, name, *, positive=False):
    if isinstance(value, bool):
        raise ValueError("Boolean is not a physical " + name)
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError("Invalid " + name)
    return value


def bimolecular_rate_cm3_molecule_s(value):
    """Convert m³ mol⁻¹ s⁻¹ to cm³ molecule⁻¹ s⁻¹ without changing kinetics."""
    return _finite(value, "bimolecular rate", positive=True) * 1e6 / AVOGADRO


def read_electronic_single_point(input_path, output_path, *, method, charge, multiplicity):
    """Reopen exact ORCA 6.1 H/C/O def2-QZVPP native SP evidence.

    Counts are derived from nuclei, never OH/methane labels. The only atomic
    exception is an explicitly requested one-electron HF calculation in the
    same basis; unsupported correlated-atom output is not guessed to be valid.
    """
    if (isinstance(charge, bool) or not isinstance(charge, int)
            or isinstance(multiplicity, bool) or not isinstance(multiplicity, int) or multiplicity < 1):
        raise ValueError("Explicit integer charge and positive spin multiplicity required")
    inp, out = bind_artifact(Path(input_path)), bind_artifact(Path(output_path))
    reference = _input(reopen(inp), method)
    if Path(output_path).resolve().with_suffix(".inp") != Path(input_path).resolve():
        raise ValueError("Electronic input/output basenames differ")
    if (reference["charge"], reference["multiplicity"]) != (charge, multiplicity):
        raise ValueError("Electronic declaration/native state mismatch")
    if {"freq", "anfreq", "numfreq"} & reference["keywords"]:
        raise ValueError("Separate electronic single point required")
    if method.split("/")[-1].casefold() != "def2-qzvpp" or set(reference["labels"]) - set(_NUCLEAR_CHARGE):
        raise ValueError("Electronic inventory policy currently covers H/C/O def2-QZVPP only")
    text, energy = _native(reopen(out), reference)
    formula = Counter(reference["labels"])
    electrons = sum(_NUCLEAR_CHARGE[s] * n for s, n in formula.items()) - charge
    basis = sum(_QZVPP_FUNCTIONS[s] * n for s, n in formula.items())
    native_basis = re.findall(r"Number of basis functions\s+\.{2,}\s+(\d+)\b", text)
    native_electrons = re.findall(r"Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)\b", text)
    if not native_basis or set(map(int, native_basis)) != {basis} \
            or not native_electrons or set(map(int, native_electrons)) != {electrons}:
        raise ValueError("Native basis/electron inventory differs from declared nuclei")
    atomic = len(reference["labels"]) == 1
    if atomic:
        if formula != {"H": 1} or (charge, multiplicity, electrons) != (0, 2, 1) \
                or method.split("/")[0].casefold() != "hf":
            raise ValueError("Unsupported atomic electronic model: explicit one-electron HF evidence required")
        result = {"energy_policy": "one_electron_hydrogen_native_hf_same_basis_v1",
                  "reference": "UHF", "frozen_core_electrons": 0, "correlated_electrons": 0,
                  "CC_residual": None}
    else:
        if method.split("/")[0].casefold() != "ccsd(t)":
            raise ValueError("Molecular composite electronic audit requires declared canonical CCSD(T)")
        frozen = sum(_CHEMICAL_CORE[s] * n for s, n in formula.items())
        native_core = re.findall(r"Frozen core treatment\s+\.{2,}\s+chemical core \((\d+) el\)", text)
        # ORCA can explicitly report no frozen core for all-hydrogen systems.
        no_core = bool(re.search(r"Frozen core treatment\s+\.{2,}\s+(?:none|no frozen core)\b", text, re.I))
        if not ((native_core and set(map(int, native_core)) == {frozen}) or (frozen == 0 and no_core)):
            raise ValueError("Native chemical frozen-core inventory differs from nuclei")
        expected_ref = "RHF" if multiplicity == 1 else "UHF"
        refs = re.findall(r"Reference Wavefunction\s+\.{2,}\s+(\w+)", text)
        correlated = re.findall(r"Number of correlated electrons\s+\.{2,}\s+(\d+)", text)
        totals = re.findall(rf"E\(CCSD\(T\)\)\s+\.{{2,}}\s+({NUMBER})", text)
        table = re.search(r"Iter\s+E\(tot\).*?\n(.*?)--- The Coupled-Cluster iterations have converged", text, re.S)
        rows = [] if table is None else re.findall(
            rf"^\s*\d+\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})", table[1], re.M)
        residual = None if not rows else float(rows[-1][-1].replace("D", "E"))
        if refs != [expected_ref] or correlated != [str(electrons - frozen)] or "QR-MO GENERATION" in text:
            raise ValueError("Native reference/correlated-electron convention mismatch")
        if len(totals) != 1 or not math.isclose(float(totals[0].replace("D", "E")), energy, rel_tol=0, abs_tol=1e-8) \
                or residual is None or not math.isfinite(residual) or not 0 <= residual <= 1e-6:
            raise ValueError("Native canonical CC completion/energy/residual evidence unresolved")
        result = {"energy_policy": "canonical_ccsd_t_chemical_core_v1", "reference": expected_ref,
                  "frozen_core_electrons": frozen, "correlated_electrons": electrons - frozen,
                  "CC_residual": residual}
    reopen(inp); reopen(out)
    return {**result, "energy_hartree": energy, "formula": dict(sorted(formula.items())),
            "charge": charge, "multiplicity": multiplicity, "electrons": electrons,
            "basis_function_count": basis, "input": inp.record(), "output": out.record(),
            "coordinates_angstrom": reference["xyz"].tolist(), "nucleus_order": list(reference["labels"]),
            "T1_diagnostic": [float(value) for value in re.findall(rf"T1 diagnostic\s+\.{{2,}}\s+({NUMBER})", text)],
            "SCF_S_squared": [float(value) for value in re.findall(rf"Expectation value of <S\*\*2>\s*:\s*({NUMBER})", text)],
            "physical_accuracy_claim_allowed": False}


def _h_ch4_classical_rotor_flags(conformer, atom_count):
    """Check this transfer's stated classical model, not a general selector."""
    rotors = [mode for mode in conformer.modes
              if type(mode).__name__ in {"LinearRotor", "NonlinearRotor"}]
    if len(rotors) != int(atom_count > 1):
        raise UnsupportedRateModelError("H + CH4 requires its complete declared classical rotor inventory")
    flags = [getattr(mode, "quantum", None) for mode in rotors]
    if any(flag is not False for flag in flags):
        raise UnsupportedRateModelError("H + CH4 frozen high-temperature model requires classical rotors, not quantum or missing controls")
    return flags


def audit_native_rate(input_path, contract_path):
    """Load the exact staged native model and audit its selected 2→2 rate.

    Native imports and statmech work occur only when called. The caller selects
    the interpreter/environment; this function does not discover or install it.
    """
    from storca.arkane_launcher import register_adapter
    from storca.arkane_rotor import reopen as reopen_rotor
    contract = validate_contract(contract_path)
    if Path(input_path).resolve() != reopen_rotor(contract["input"]):
        raise ValueError("Native audit input differs from admitted rotor contract")
    provenance = read_document(reopen_rotor(contract["provenance"]))
    recipe_source = reopen_rotor(provenance["composite_recipe"])
    _, assessment = reopen_composite_recipe_assessment(recipe_source)
    # The wire representation uses arrays, independently of dataclass tuples.
    recipe = json.loads(json.dumps(assessment["recipe"]))
    by_label = {point["label"]: point for point in assessment["points"]}
    if set(by_label) != {point["label"] for point in contract["stationary_points"]}:
        raise ValueError("Native/composite stationary-point inventories differ")
    staged_folder = Path(input_path).resolve().parent
    def inventory():
        return {str(path): (artifact(path)["sha256"], path.stat().st_mtime_ns)
                for path in staged_folder.rglob("*") if path.is_file()}
    before = inventory()
    adapter = register_adapter(contract)
    from arkane.input import load_input_file
    from arkane.statmech import StatMechJob
    from arkane.kinetics import KineticsJob
    from rmgpy import constants
    from rmgpy.reaction import Reaction
    jobs, reactions, species, transition_states, networks, level = load_input_file(str(Path(input_path).resolve()))
    statjobs = [job for job in jobs if isinstance(job, StatMechJob)]
    kineticjobs = [job for job in jobs if isinstance(job, KineticsJob)]
    if len(statjobs) != len(by_label) or len(kineticjobs) != 1 or len(jobs) != len(by_label) + 1 \
            or networks or len(reactions) != 1:
        raise UnsupportedRateModelError("Native audit supports one isolated declared elementary reaction only")
    points = {}
    for job in statjobs:
        row = by_label.get(job.species.label)
        if row is None or job.includeHinderedRotors or job.applyAtomEnergyCorrections \
                or job.applyBondEnergyCorrections or job.useIsodesmicReactions:
            raise ValueError("Native point/uncorrected RRHO assumptions differ from recipe")
        if float(job.frequencyScaleFactor) != recipe["frequency_scale_factor"]:
            raise ValueError("Native frequency scale differs from frozen recipe")
        job.execute(output_directory=None, plot=False, pdep=False)
        conformer = job.species.conformer
        if recipe["label"] == "transfer-h_ch4":
            _h_ch4_classical_rotor_flags(conformer, row["atom_count"])
        audit = adapter.storca_audits.get(job.species.label)
        if not audit or audit["status"] != "passed" or not all(audit["preserved"].values()):
            raise ValueError("Native mode/rotor admission is incomplete")
        expected_E0 = assessment["gauged_electronic_energies_J_mol"][job.species.label] \
            + row["raw_harmonic_ZPE_J_mol"] * recipe["frequency_scale_factor"] / 1.014
        observed_E0 = _finite(conformer.E0.value_si, "native E0")
        if not math.isclose(observed_E0, expected_E0, rel_tol=0, abs_tol=0.05):
            raise ValueError("Native E0 differs from declared electronic gauge/raw-ZPE convention")
        if conformer.spin_multiplicity != row["multiplicity"] or conformer.optical_isomers != row["optical_isomers"]:
            raise ValueError("Native spin/optical declaration changed")
        points[job.species.label] = {"E0_J_mol": observed_E0, "expected_E0_J_mol": expected_E0,
            "mode_types": [type(mode).__name__ for mode in conformer.modes],
            "spin_multiplicity": conformer.spin_multiplicity, "optical_isomers": conformer.optical_isomers,
            "rotor_admission": audit}
    if set(points) != set(by_label):
        raise ValueError("Native point evidence is incomplete or duplicated")
    job = kineticjobs[0]
    reaction, ts = job.reaction, job.reaction.transition_state
    if [item.label for item in reaction.reactants] != recipe["reactant_labels"] \
            or [item.label for item in reaction.products] != recipe["product_labels"] \
            or ts.label != recipe["transition_state_label"] or not reaction.is_balanced() \
            or len(reaction.reactants) != 2 or len(reaction.products) != 2 or reaction.degeneracy != 1:
        raise ValueError("Native reaction differs from exact declared two-to-two recipe")
    grid = [_finite(value, "temperature", positive=True) for value in job.Tlist.value_si]
    if grid != recipe["temperatures_K"]:
        raise ValueError("Native temperature grid differs from frozen prediction inputs")
    E_reac = sum(item.conformer.E0.value_si for item in reaction.reactants)
    E_prod = sum(item.conformer.E0.value_si for item in reaction.products)
    E_ts = ts.conformer.E0.value_si
    selected = "none" if ts.tunneling is None else type(ts.tunneling).__name__.casefold()
    if selected != recipe["tunneling"].casefold() or selected not in {"none", "wigner", "eckart"}:
        raise ValueError("Native tunneling model differs from explicit recipe")
    if selected == "eckart" and (E_ts <= E_reac or E_ts <= E_prod):
        raise PhysicalRateValidationError("Eckart requires strictly positive forward and reverse E0 barriers")
    if reaction.kinetics is not None:
        raise ValueError("Existing native kinetics would bypass TST")
    job.generate_kinetics()
    if not job.usedTST:
        raise ValueError("Native selected-model rate did not use TST")
    reverse = Reaction(reactants=reaction.products, products=reaction.reactants, transition_state=ts)
    rows = []
    for temperature in grid:
        kf = _finite(reaction.calculate_tst_rate_coefficient(temperature), "forward rate", positive=True)
        kr = _finite(reverse.calculate_tst_rate_coefficient(temperature), "reverse rate", positive=True)
        kappa = 1. if ts.tunneling is None else _finite(ts.tunneling.calculate_tunneling_factor(temperature), "tunneling factor", positive=True)
        standard_volume = constants.R * temperature / 101325.
        qts = ts.get_partition_function(temperature) / standard_volume
        qr = math.prod(item.get_partition_function(temperature) / standard_volume for item in reaction.reactants)
        qp = math.prod(item.get_partition_function(temperature) / standard_volume for item in reaction.products)
        manual = constants.kB * temperature / constants.h * qts / qr \
            * math.exp(-(E_ts - E_reac) / constants.R / temperature) * kappa
        equilibrium = qp / qr * math.exp(-(E_prod - E_reac) / constants.R / temperature)
        native_eq = reaction.get_equilibrium_constant(temperature, type="Kc")
        errors = {"native_TST_vs_explicit_partition": abs(kf / manual - 1),
                  "forward_reverse_vs_partition_equilibrium": abs((kf / kr) / equilibrium - 1),
                  "forward_reverse_vs_native_equilibrium": abs((kf / kr) / native_eq - 1)}
        if not all(math.isfinite(value) for value in errors.values()) or max(errors.values()) > 1e-8:
            raise PhysicalRateValidationError("Native TST/detailed-balance identity failed")
        fitted = _finite(reaction.kinetics.get_rate_coefficient(temperature), "fitted rate", positive=True)
        rows.append({"temperature_K": temperature, "forward_direct_m3_mol_s": kf,
                     "reverse_direct_m3_mol_s": kr,
                     "forward_direct_cm3_molecule_s": bimolecular_rate_cm3_molecule_s(kf),
                     "forward_fitted_m3_mol_s": fitted, "log_fit_over_direct": math.log(fitted / kf),
                     "tunneling_factor": kappa, "identity_relative_errors": errors})
    for binding in (contract["input"], contract["provenance"], provenance["composite_recipe"]):
        reopen_rotor(binding)
    validate_contract(contract)
    if before != inventory():
        raise ValueError("Native rate replay changed admitted or previous staged artifacts")
    return {"schema_version": 1, "status": "conditional_native_rate_consistency_passed",
            "reaction_label": reaction.label, "selected_tunneling": selected,
            "temperatures_K": grid, "stationary_points": points, "rate_rows": rows,
            "barriers_J_mol": {"forward_E0": E_ts - E_reac, "reverse_E0": E_ts - E_prod},
            "declared_reaction_degeneracy": float(reaction.degeneracy),
            "input": contract["input"], "contract": artifact(contract_path),
            "composite_recipe": provenance["composite_recipe"],
            "native_interpreter": artifact(sys.executable),
            "native_modules": {name: artifact(importlib.import_module(name).__file__) for name in (
                "arkane.input", "arkane.statmech", "arkane.kinetics", "rmgpy.reaction", "rmgpy.constants")},
            "experimental_answers_read": False, "electronic_structure_jobs_launched": 0,
            "all_admitted_and_previous_staged_artifacts_unchanged": True,
            "physical_rate_accuracy_claim_allowed": False, "calibrated_uncertainty_claim_allowed": False,
            "lifetime_claim_allowed": False, "network_reinjection_allowed": False,
            "production_defaults_changed": False,
            "limitations": ["Internal consistency does not validate physical RRHO/TST/tunneling assumptions.",
                            "Counting and electronic-state assumptions are explicit declarations, not calibrated facts.",
                            "No experimental comparison, pressure dependence or mechanism completeness claim."]}


def _write_new_document(path, document):
    with Path(path).open("x") as stream:
        json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def _validate_audit_document(folder, audit):
    """Reopen portable native sources and check the declared numerical claims."""
    from storca.arkane_rotor import reopen as reopen_rotor
    if audit.get("status") != "conditional_native_rate_consistency_passed":
        raise ValueError("Native rate audit did not complete conditional consistency checks")
    for field in ("physical_rate_accuracy_claim_allowed", "calibrated_uncertainty_claim_allowed",
                  "lifetime_claim_allowed", "network_reinjection_allowed", "production_defaults_changed",
                  "experimental_answers_read"):
        if audit.get(field) is not False:
            raise ValueError("Native rate audit contains an inadmissible scientific promotion")
    if audit.get("electronic_structure_jobs_launched") != 0 \
            or audit.get("all_admitted_and_previous_staged_artifacts_unchanged") is not True:
        raise ValueError("Native replay preservation/job scope evidence is absent")
    expected = {"input": folder / "input.py", "contract": folder / "arkane-rotor-contract.json",
                "composite_recipe": folder / "composite-rate-recipe.json"}
    for name, path in expected.items():
        if reopen_rotor(audit[name]) != path:
            raise ValueError("Native rate audit source belongs to a different staged model")
    validate_contract(folder / "arkane-rotor-contract.json")
    recipe, assessed = reopen_composite_recipe_assessment(folder / "composite-rate-recipe.json")
    if audit.get("reaction_label") != recipe.label or audit.get("selected_tunneling") != recipe.tunneling.casefold() \
            or audit.get("temperatures_K") != list(recipe.temperatures_K):
        raise ValueError("Native rate audit differs from the exact recipe/grid")
    modules = audit.get("native_modules")
    required_modules = {"arkane.input", "arkane.statmech", "arkane.kinetics", "rmgpy.reaction", "rmgpy.constants"}
    if not isinstance(modules, dict) or set(modules) != required_modules:
        raise ValueError("Native rate audit module inventory is incomplete")
    for binding in [*modules.values(), audit["native_interpreter"]]:
        reopen_rotor(binding)
    points = audit.get("stationary_points")
    if not isinstance(points, dict) or set(points) != {point["label"] for point in assessed["points"]}:
        raise ValueError("Native rate audit stationary inventory is incomplete")
    for row in assessed["points"]:
        point = points[row["label"]]
        expected_e0 = assessed["gauged_electronic_energies_J_mol"][row["label"]] \
            + row["raw_harmonic_ZPE_J_mol"] * recipe.frequency_scale_factor / 1.014
        if not math.isclose(_finite(point["E0_J_mol"], "retained native E0"), expected_e0, rel_tol=0, abs_tol=.05) \
                or point.get("spin_multiplicity") != row["multiplicity"] \
                or point.get("optical_isomers") != row["optical_isomers"]:
            raise ValueError("Retained native atomic/molecular state/E0 differs from source recipe")
    if audit.get("declared_reaction_degeneracy") != 1.0:
        raise ValueError("Retained native rate contains an undeclared additional counting factor")
    barriers = audit.get("barriers_J_mol")
    expected_barriers = {direction + "_E0": points[recipe.transition_state_label]["E0_J_mol"]
                         - sum(points[label]["E0_J_mol"] for label in labels)
                         for direction, labels in (("forward", recipe.reactant_labels), ("reverse", recipe.product_labels))}
    if not isinstance(barriers, dict) or set(barriers) != set(expected_barriers) \
            or any(not math.isclose(_finite(barriers[name], "retained E0 barrier"), value,
                                    rel_tol=0, abs_tol=1e-6) for name, value in expected_barriers.items()):
        raise ValueError("Retained native barrier inventory differs from its actual endpoint E0 values")
    if recipe.tunneling.casefold() == "eckart" and any(value <= 0 for value in expected_barriers.values()):
        raise ValueError("Retained Eckart model lacks positive forward/reverse E0 barriers")
    rows = audit.get("rate_rows")
    if not isinstance(rows, list) or [row.get("temperature_K") for row in rows] != list(recipe.temperatures_K):
        raise ValueError("Retained native rate grid is incomplete or altered")
    rate_points = []
    for row in rows:
        kf = _finite(row["forward_direct_m3_mol_s"], "retained forward rate", positive=True)
        _finite(row["reverse_direct_m3_mol_s"], "retained reverse rate", positive=True)
        converted = bimolecular_rate_cm3_molecule_s(kf)
        if not math.isclose(_finite(row["forward_direct_cm3_molecule_s"], "retained converted rate", positive=True), converted, rel_tol=1e-12):
            raise ValueError("Retained native bimolecular rate units changed")
        errors = row.get("identity_relative_errors")
        required_errors = {"native_TST_vs_explicit_partition", "forward_reverse_vs_partition_equilibrium",
                           "forward_reverse_vs_native_equilibrium"}
        if not isinstance(errors, dict) or set(errors) != required_errors \
                or any(not 0 <= _finite(value, "retained identity residual") <= 1e-8 for value in errors.values()):
            raise ValueError("Retained native TST/detailed-balance audit failed")
        rate_points.append({"temperature_K": row["temperature_K"], "rate_constant": converted,
                            "rate_units": "cm^3 molecule^-1 s^-1"})
    return rate_points


def run_native_rate_audit(staged_folder, *, rmg_env):
    """One bounded, source-bound selected-environment replay, with checked reuse.

    A charged/ambiguous prior invocation is never silently rerun. This function
    dispatches only the existing native kinetics command supervisor, not ORCA.
    """
    if not isinstance(rmg_env, str) or not rmg_env.strip():
        raise ValueError("Select the RMG/Arkane environment explicitly with --rmg-env")
    folder = Path(staged_folder).resolve(strict=True)
    validate_contract(folder / "arkane-rotor-contract.json")
    reopen_composite_recipe_assessment(folder / "composite-rate-recipe.json")
    started_path = folder / "native-rate-audit-started.json"
    execution_path = folder / "native-rate-audit-execution.json"
    output_path = folder / "native-rate-audit.json"
    helper = Path(__file__).resolve()
    expected = {"schema_version": 1, "kind": "storca_known_reaction_native_rate_replay",
                "rmg_env": rmg_env, "helper": artifact(helper), "input": artifact(folder / "input.py"),
                "contract": artifact(folder / "arkane-rotor-contract.json"),
                "composite_recipe": artifact(folder / "composite-rate-recipe.json"),
                "timeout_seconds": 21600.0, "electronic_structure_jobs_launched": 0}
    if started_path.exists() or execution_path.exists() or output_path.exists():
        if not all(path.is_file() and not path.is_symlink() for path in (started_path, execution_path, output_path)):
            raise ValueError("Prior native replay is active, failed or ambiguous; preserved evidence forbids redispatch")
        if read_document(started_path) != expected:
            raise ValueError("Native replay source/environment contract changed")
        execution = read_document(execution_path)
        if execution.get("status") != "completed" or execution.get("returncode") != 0 \
                or execution.get("started") != artifact(started_path) or execution.get("audit") != artifact(output_path):
            raise ValueError("Prior native replay did not complete a source-bound successful audit")
        audit = read_document(output_path)
        return {**audit, "rate_points": _validate_audit_document(folder, audit),
                "execution_receipt": artifact(execution_path), "reused": True}
    _write_new_document(started_path, expected)  # Exclusive durable reservation before any subprocess.
    command = ["conda", "run", "--no-capture-output", "-n", rmg_env, "python", "-B", str(helper),
               str(folder / "input.py"), str(folder / "arkane-rotor-contract.json"), str(output_path)]
    started = time.monotonic()
    try:
        completed = run_bounded_kinetics_command(command, cwd=folder, timeout_seconds=21600.0)
        receipt = {"schema_version": 1, "status": "completed" if completed.returncode == 0 else "failed_execution",
                   "started": artifact(started_path), "command": command, "returncode": completed.returncode,
                   "stdout": completed.stdout, "stderr": completed.stderr,
                   "elapsed_seconds": time.monotonic() - started,
                   "audit": artifact(output_path) if output_path.is_file() else None}
        if completed.returncode != 0 or not output_path.is_file():
            receipt["status"] = "failed_execution"
            if output_path.is_file():
                failure = read_document(output_path)
                if failure.get("status") == "failed_native_rate_admission":
                    receipt.update({"status": "failed_" + failure["failure_category"],
                                    "failure_reason": failure["failure_reason"]})
            _write_new_document(execution_path, receipt)
            raise ValueError("Native rate replay failed; execution evidence retained, no automatic retry")
        audit = read_document(output_path)
        try:
            rates = _validate_audit_document(folder, audit)
        except Exception as error:
            receipt.update({"status": "failed_reader_or_validation", "failure_reason": str(error)})
            _write_new_document(execution_path, receipt)
            raise
        _write_new_document(execution_path, receipt)
        return {**audit, "rate_points": rates, "execution_receipt": artifact(execution_path), "reused": False}
    except BaseException as error:
        if not execution_path.exists():
            _write_new_document(execution_path, {"schema_version": 1, "status": "failed_execution",
                "started": artifact(started_path), "command": command, "returncode": None,
                "failure_type": type(error).__name__, "failure_reason": str(error),
                "stdout": getattr(error, "stdout", getattr(error, "output", "")),
                "stderr": getattr(error, "stderr", ""), "elapsed_seconds": time.monotonic() - started})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("contract", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink() or not args.output.parent.is_dir():
        raise ValueError("Native rate audit requires a fresh output file in an existing directory")
    try:
        result = audit_native_rate(args.input, args.contract)
    except Exception as error:
        category = ("unsupported_model" if isinstance(error, UnsupportedRateModelError) else
                    "physical_validation" if isinstance(error, PhysicalRateValidationError) else "reader_validation")
        _write_new_document(args.output, {"schema_version": 1, "status": "failed_native_rate_admission",
            "failure_category": category, "failure_type": type(error).__name__, "failure_reason": str(error),
            "input": artifact(args.input), "contract": artifact(args.contract), "helper": artifact(__file__),
            "electronic_structure_jobs_launched": 0, "physical_rate_accuracy_claim_allowed": False})
        raise
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
