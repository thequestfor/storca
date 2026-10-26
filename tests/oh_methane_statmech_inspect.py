"""Read-only native Arkane replay for an isolated five-point OH/CH4 rate audit.

Run with the explicitly selected RMG interpreter and ``-B``. This helper does
not run Arkane's main entry point, write Arkane outputs, or dispatch electronic
structure jobs. Its sole persistent write is a fresh requested JSON diagnostic.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
TEMPERATURE_SOURCE = ROOT / "benchmarks/stability/reference_sources/2026-10-07-berg2024-oh-methane/source-audit.json"


def artifact(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError("Expected a regular source file: " + str(path))
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def snapshot(paths):
    return {str(path): {**artifact(path), "mtime_ns": path.stat().st_mtime_ns}
            for path in sorted(set(paths))}


def finite(value, name, positive=False):
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError("Nonfinite or invalid " + name)
    return value


def plain(value):
    """Serialize native arrays without silently stringifying physical numbers."""
    if hasattr(value, "tolist"):
        return plain(value.tolist())
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if isinstance(value, dict):
        return {str(key): plain(item) for key, item in value.items()}
    if value is None or isinstance(value, (str, bool, int)):
        return value
    return finite(value, "native scalar")


def quantity(value):
    if value is None:
        return None
    return {"value": plain(value.value), "units": value.units,
            "value_si": plain(value.value_si)}


def mode_record(mode):
    row = {"class": type(mode).__name__, "native_repr": repr(mode)}
    for name in ("mass", "inertia", "frequencies"):
        if hasattr(mode, name):
            row[name] = quantity(getattr(mode, name))
    for name in ("symmetry", "quantum"):
        if hasattr(mode, name):
            row[name] = plain(getattr(mode, name))
    return row


def inspect(input_path, contract_path):
    from storca.arkane_rotor import reopen, validate_contract
    from storca.arkane_launcher import register_adapter

    contract = validate_contract(contract_path)
    if reopen(contract["input"]) != input_path:
        raise ValueError("Input is not the exact source admitted by the rotor contract")
    points = contract["stationary_points"]
    if len(points) != 5 or Counter(point["kind"] for point in points) != {
            "species": 4, "transition_state": 1}:
        raise ValueError("This diagnostic requires four fragments and one transition state")
    adapter = register_adapter(contract)
    from arkane.input import load_input_file
    from arkane.statmech import StatMechJob
    from arkane.kinetics import KineticsJob
    from rmgpy import constants
    from rmgpy.kinetics.tunneling import Eckart, Wigner
    from rmgpy.reaction import Reaction

    jobs, reactions, species, transition_states, networks, level = load_input_file(str(input_path))
    statjobs = [job for job in jobs if isinstance(job, StatMechJob)]
    kineticjobs = [job for job in jobs if isinstance(job, KineticsJob)]
    if len(statjobs) != 5 or len(kineticjobs) != 1 or len(jobs) != 6 or networks \
            or len(reactions) != 1 or len(species) != 4 or len(transition_states) != 1:
        raise ValueError("Unexpected native jobs/network; only the scoped five-point replay is supported")
    if {job.species.label for job in statjobs} != {point["label"] for point in points}:
        raise ValueError("Native stationary-point labels differ from the contract")
    if any(job.includeHinderedRotors or job.applyAtomEnergyCorrections
           or job.applyBondEnergyCorrections or job.useIsodesmicReactions for job in statjobs):
        raise ValueError("This diagnostic admits only the declared uncorrected RRHO audit")
    if len({float(job.frequencyScaleFactor) for job in statjobs}) != 1:
        raise ValueError("Mixed frequency scaling across the five-point audit")

    point_rows, warnings = {}, []
    point_by_label = {point["label"]: point for point in points}
    for job in statjobs:
        if Path(job.path).resolve(strict=True) != reopen(point_by_label[job.species.label]["statmech_file"]):
            raise ValueError("Statmech job path differs from its admitted source")
        job.execute(output_directory=None, plot=False, pdep=False)
        point = point_by_label[job.species.label]
        conformer = job.species.conformer
        native_log = adapter(point["sources"]["staged_output"]["path"])
        auto_optical, auto_symmetry, auto_group = native_log.get_symmetry_properties()
        audit = adapter.storca_audits.get(job.species.label)
        if not audit or audit["status"] != "passed" or not all(audit["preserved"].values()):
            raise ValueError("Missing successful process-local rotor admission")
        info = job.supporting_info
        if len(info) != 14 or info[0] != job.species.label:
            raise ValueError("Installed native supporting-info schema differs from audited interface")
        if not math.isclose(conformer.E0.value_si, info[9], rel_tol=1e-13, abs_tol=1e-7):
            raise ValueError("Loaded native E0 differs from its supporting evidence")
        if conformer.spin_multiplicity != point["multiplicity"] \
                or conformer.optical_isomers != point["optical_isomers"]:
            raise ValueError("Loaded spin/optical declaration differs from the contract")
        row = {
            "label": job.species.label, "kind": point["kind"], "formula": point["formula"],
            "sources": point["sources"], "statmech_file": point["statmech_file"],
            "E0_J_mol": finite(conformer.E0.value_si, "E0"),
            "electronic_energy_J_mol": None if info[7] is None else finite(info[7], "electronic energy"),
            "electronic_plus_scaled_ZPE_J_mol": None if info[8] is None else finite(info[8], "electronic+ZPE"),
            "scaled_ZPE_J_mol": None if info[7] is None or info[8] is None else finite(info[8] - info[7], "ZPE"),
            "frequency_scale_factor": float(job.frequencyScaleFactor),
            "native_ZPE_scale_factor": float(job.frequencyScaleFactor) / 1.014,
            "unscaled_positive_frequencies_cm_inverse": plain(info[5]),
            "unscaled_negative_frequency_cm_inverse": plain(info[6]),
            "scaled_TS_frequency": quantity(getattr(job.species, "frequency", None)),
            "modes": [mode_record(mode) for mode in conformer.modes],
            "spin_multiplicity": int(conformer.spin_multiplicity),
            "imposed_optical_isomers": int(point["optical_isomers"]),
            "loaded_optical_isomers": int(conformer.optical_isomers),
            "loaded_external_symmetry": int(info[1]),
            "automatic_symmetry": {"point_group": str(auto_group), "symmetry_number": int(auto_symmetry),
                                   "optical_isomers": int(auto_optical)},
            "symbols": plain(info[10]), "coordinates_angstrom": plain(info[11]),
            "rotor_admission": audit,
        }
        if auto_optical != conformer.optical_isomers:
            warnings.append({"label": job.species.label, "kind": "automatic_optical_vs_imposed",
                             "automatic": int(auto_optical), "imposed": int(conformer.optical_isomers),
                             "meaning": "Diagnostic only; no silent chirality/degeneracy correction is applied."})
        point_rows[job.species.label] = row

    job = kineticjobs[0]
    reaction = job.reaction
    ts = reaction.transition_state
    if len(reaction.reactants) != 2 or len(reaction.products) != 2 or reaction.degeneracy != 1:
        raise ValueError("The scoped audit requires a two-to-two reaction with no extra declared degeneracy")
    def formulas(items):
        return sorted(tuple(sorted(point_by_label[item.label]["formula"].items())) for item in items)
    if formulas(reaction.reactants) != sorted([(("H", 1), ("O", 1)), (("C", 1), ("H", 4))]) \
            or formulas(reaction.products) != sorted([(("H", 2), ("O", 1)), (("C", 1), ("H", 3))]) \
            or ts.label not in transition_states or not reaction.is_balanced():
        raise ValueError("Input is not the admitted OH+CH4 to H2O+CH3 reaction")
    grid = [finite(value, "temperature", positive=True) for value in job.Tlist.value_si]
    if len(grid) < 3 or grid != sorted(set(grid)):
        raise ValueError("A unique increasing native fit grid with at least three temperatures is required")
    observation_set = json.loads(TEMPERATURE_SOURCE.read_text())["audited_observation_sets"][0]
    if not observation_set["transcription_frozen"] or len(observation_set["observations"]) != 14:
        raise ValueError("Exact frozen fourteen-temperature evaluation source required")
    evaluation_points = [{"point_id": point["point_id"],
                          "temperature_K": finite(point["temperature_K"], "measurement temperature", positive=True)}
                         for point in observation_set["observations"]]
    evaluation_grid = sorted(set(grid) | {point["temperature_K"] for point in evaluation_points})
    if reaction.kinetics is not None:
        raise ValueError("Pre-existing kinetics would bypass native TST generation")
    selected = "classical" if ts.tunneling is None else type(ts.tunneling).__name__
    if selected not in ("classical", "Wigner", "Eckart"):
        raise ValueError("Unsupported selected tunneling model")
    job.generate_kinetics()
    if not job.usedTST:
        raise ValueError("Native generation did not use TST")
    selected_fit = reaction.kinetics
    original_tunneling = ts.tunneling
    E_reac = sum(item.conformer.E0.value_si for item in reaction.reactants)
    E_prod = sum(item.conformer.E0.value_si for item in reaction.products)
    E_ts = ts.conformer.E0.value_si
    freq = finite(ts.frequency.value_si, "TS frequency")
    if freq >= 0:
        raise ValueError("Transition-state frequency is not negative")
    models = {
        "classical": None,
        "Wigner": Wigner(frequency=(freq, "cm^-1")),
        "Eckart": Eckart(frequency=(freq, "cm^-1"), E0_reac=(E_reac, "J/mol"),
                         E0_TS=(E_ts, "J/mol"), E0_prod=(E_prod, "J/mol")),
    }
    rates = {}
    try:
        for name, model in models.items():
            ts.tunneling = model
            forward = Reaction(reactants=reaction.reactants, products=reaction.products, transition_state=ts)
            reverse = Reaction(reactants=reaction.products, products=reaction.reactants, transition_state=ts)
            fit_job = KineticsJob(forward, Tlist=(grid, "K"), three_params=job.three_params)
            fit_job.generate_kinetics()
            rows = []
            for T in evaluation_grid:
                kf = finite(forward.calculate_tst_rate_coefficient(T), "forward rate", positive=True)
                kr = finite(reverse.calculate_tst_rate_coefficient(T), "reverse rate", positive=True)
                kappa = 1.0 if model is None else finite(model.calculate_tunneling_factor(T), "tunneling", positive=True)
                # Reproduce the installed native standard-concentration convention, not a guessed unit conversion.
                standard_volume = constants.R * T / 101325.0
                Qts = ts.get_partition_function(T) / standard_volume
                Qreac = math.prod(item.get_partition_function(T) / standard_volume for item in reaction.reactants)
                Qprod = math.prod(item.get_partition_function(T) / standard_volume for item in reaction.products)
                manual = constants.kB * T / constants.h * Qts / Qreac * math.exp(-(E_ts - E_reac) / constants.R / T) * kappa
                Keq_Q = Qprod / Qreac * math.exp(-(E_prod - E_reac) / constants.R / T)
                Keq_native = finite(reaction.get_equilibrium_constant(T, type="Kc"), "native equilibrium", positive=True)
                consistency = {"native_TST_vs_explicit_Q": abs(kf / manual - 1.0),
                               "forward_reverse_vs_partition_equilibrium": abs((kf / kr) / Keq_Q - 1.0),
                               "forward_reverse_vs_native_equilibrium": abs((kf / kr) / Keq_native - 1.0)}
                if max(consistency.values()) > 1e-8:
                    raise ValueError("Native TST/equilibrium identity mismatch at " + str(T))
                fit = finite(forward.kinetics.get_rate_coefficient(T), "fitted rate", positive=True)
                row = {"temperature_K": T, "in_original_fit_grid": T in grid,
                       "forward_direct_m3_mol_s": kf, "reverse_direct_m3_mol_s": kr,
                       "tunneling_factor": kappa, "explicit_partition_TST_m3_mol_s": manual,
                       "equilibrium_Kc_dimensionless_native": Keq_native,
                       "equilibrium_Kc_dimensionless_partition": Keq_Q,
                       "forward_fitted_m3_mol_s": fit, "fit_over_direct": fit / kf,
                       "log_fit_over_direct": math.log(fit / kf), "identity_relative_errors": consistency}
                if name == selected:
                    actual_fit = finite(selected_fit.get_rate_coefficient(T), "selected input fitted rate", positive=True)
                    if not math.isclose(actual_fit, fit, rel_tol=1e-10, abs_tol=0):
                        raise ValueError("Regenerated selected-model fit differs from actual input fit")
                    row["actual_input_fitted_m3_mol_s"] = actual_fit
                rows.append(row)
            rates[name] = {"native_tunneling_repr": repr(model), "native_fit_repr": repr(forward.kinetics),
                           "fit_parameters": {key: quantity(getattr(forward.kinetics, key))
                                              for key in ("A", "n", "Ea", "T0", "Tmin", "Tmax")},
                           "maximum_absolute_log_fit_error_original_grid": max(abs(row["log_fit_over_direct"])
                                                                                 for row in rows if row["in_original_fit_grid"]),
                           "maximum_absolute_log_fit_error_evaluation_grid": max(abs(row["log_fit_over_direct"]) for row in rows),
                           "grid": rows}
    finally:
        ts.tunneling = original_tunneling
        reaction.kinetics = selected_fit
    if set(adapter.storca_audits) != set(point_rows):
        raise ValueError("Incomplete adapter point evidence")
    native_modules = ("arkane.input", "arkane.statmech", "arkane.kinetics", "arkane.common",
                      "arkane.ess.orca", "arkane.ess.factory", "rmgpy.reaction", "rmgpy.species",
                      "rmgpy.statmech.conformer", "rmgpy.statmech.rotation", "rmgpy.statmech.vibration",
                      "rmgpy.statmech.translation", "rmgpy.kinetics.tunneling", "rmgpy.kinetics.arrhenius",
                      "rmgpy.quantity", "rmgpy.constants")
    import rmgpy
    return {
        "status": "passed", "schema_version": 1, "reaction_label": reaction.label,
        "reactants": [item.label for item in reaction.reactants], "products": [item.label for item in reaction.products],
        "transition_state": ts.label, "level_of_theory": repr(level), "selected_tunneling": selected,
        "declared_reaction_degeneracy": float(reaction.degeneracy), "additional_degeneracy_multiplier": 1.0,
        "original_fit_temperature_grid_K": grid, "evaluation_temperature_grid_K": evaluation_grid,
        "extra_evaluation_temperature_source": artifact(TEMPERATURE_SOURCE),
        "measurement_temperature_points": evaluation_points,
        "experimental_rates_read_into_physical_model": False,
        "extra_temperatures_used_for_fitting": False,
        "stationary_points": point_rows, "rates": rates,
        "barriers_J_mol": {"forward_E0": E_ts - E_reac, "reverse_E0": E_ts - E_prod,
                            "reaction_E0": E_prod - E_reac},
        "actual_input_fit_repr": repr(selected_fit), "warnings": warnings,
        "native": {"rmg_version": str(rmgpy.__version__), "python": sys.version,
                   "interpreter": artifact(sys.executable),
                   "symmetry_executable": artifact(Path(sys.executable).resolve().parent / "symmetry"),
                   "modules": {name: artifact(importlib.import_module(name).__file__) for name in native_modules}},
        "limitations": ["Internal statmech/TST/fit diagnostic, not experimental calibration or a lifetime claim.",
                        "No additional reaction-path degeneracy or chirality correction has been applied.",
                        "RRHO and native one-dimensional Wigner/Eckart assumptions remain unvalidated physically.",
                        "Two-to-two equilibrium cancels the native TST/thermo standard-pressure convention difference.",
                        "A fit residual is reported, not hidden or treated as evidence of rate accuracy."],
        "electronic_structure_jobs_launched": 0, "production_defaults_changed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("contract", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    input_path = args.input.resolve(strict=True)
    contract_path = args.contract.resolve(strict=True)
    output = args.output.absolute()
    if output.exists() or output.is_symlink() or not output.parent.is_dir():
        raise ValueError("Requested diagnostic must be a fresh file in an existing directory")
    output = output.parent.resolve(strict=True) / output.name
    workdir = input_path.parent
    from storca.arkane_rotor import validate_contract
    contract = validate_contract(contract_path)
    sources = {input_path, contract_path, Path(__file__).resolve(), TEMPERATURE_SOURCE.resolve(strict=True)}
    sources.update(Path(record["path"]).resolve(strict=True) for record in contract["adapter_sources"].values())
    sources.add(Path(contract["provenance"]["path"]).resolve(strict=True))
    for point in contract["stationary_points"]:
        sources.add(Path(point["statmech_file"]["path"]).resolve(strict=True))
        sources.update(Path(record["path"]).resolve(strict=True) for record in point["sources"].values())
    staged_files = {path.resolve() for path in workdir.rglob("*") if path.is_file()}
    if output in sources or output in staged_files:
        raise ValueError("Diagnostic output overlaps admitted or previous evidence")
    before = snapshot(sources | staged_files)
    original_cwd = Path.cwd()
    original_env = {name: os.environ.get(name) for name in ("MPLCONFIGDIR", "XDG_CACHE_HOME", "PATH")}
    document, error = None, None
    try:
        with tempfile.TemporaryDirectory(prefix="storca-statmech-inspect-") as temporary:
            os.environ["MPLCONFIGDIR"] = temporary
            os.environ["XDG_CACHE_HOME"] = temporary
            # The caller selects this environment by interpreter, without an
            # implicit shell activation. Bind native symmetry to that same env.
            symmetry = Path(sys.executable).resolve().parent / "symmetry"
            if not symmetry.is_file() or not os.access(symmetry, os.X_OK):
                raise ValueError("Selected RMG interpreter environment lacks native symmetry")
            os.environ["PATH"] = str(symmetry.parent) + os.pathsep + os.environ.get("PATH", "")
            document = inspect(input_path, contract_path)
    except Exception as caught:
        error = caught
        document = {"schema_version": 1, "status": "failed", "failure_type": type(caught).__name__,
                    "failure_reason": str(caught), "electronic_structure_jobs_launched": 0}
    finally:
        os.chdir(original_cwd)
        for name, value in original_env.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
    after_files = {path.resolve() for path in workdir.rglob("*") if path.is_file()}
    try:
        unchanged = before == snapshot(sources | staged_files) and staged_files == after_files
    except (OSError, ValueError):
        unchanged = False
    document.update({"input": artifact(input_path), "contract": artifact(contract_path),
                     "helper": artifact(__file__), "prior_source_inventory": before,
                     "all_admitted_and_previous_artifacts_unchanged": unchanged,
                     "only_persistent_write": str(output)})
    if not unchanged:
        document.update({"status": "failed", "failure_reason": "Admitted/previous artifact preservation failed"})
        error = error or ValueError(document["failure_reason"])
    with output.open("x") as stream:
        json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    if error:
        raise error
    print("NATIVE_STATMECH_INSPECTION_PASSED " + str(output), flush=True)


if __name__ == "__main__":
    main()
