"""Append-only, zero-ORCA electronic-state replay of the retained OH/CH4 audit.

This is an experiment-specific consistency diagnostic, not a production rate
adapter. Native RMG is used only for typed modal thermochemistry/TST and freshly
constructed Wigner/Eckart factors. Retained Arkane input is never executed.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
POINTS = ("OH", "CH4", "H2O", "CH3", "TS")
SPINS = {"OH": 2, "CH4": 1, "H2O": 1, "CH3": 2, "TS": 2}
FORMULAS = {"OH": {"O": 1, "H": 1}, "CH4": {"C": 1, "H": 4},
            "H2O": {"O": 1, "H": 2}, "CH3": {"C": 1, "H": 3},
            "TS": {"C": 1, "O": 1, "H": 5}}
MODE_COUNTS = {"OH": 1, "CH4": 9, "H2O": 3, "CH3": 6, "TS": 14}
REFERENCE = ROOT / "benchmarks/stability/reference_sources/2026-10-07-berg2024-oh-methane/source-audit.json"
LEVEL_SOURCE = "https://truhlar.chem.umn.edu/sites/truhlar.chem.umn.edu/files/2022-05/final704.pdf#equation-29"
CENTROID_SOURCE = "https://www2.sci.u-szeged.hu/czako/papers/PCCP_22_14560_2020.pdf"
DELTA_J_MOL = 6.62607015e-34 * 299792458.0 * 100.0 * 6.02214076e23 * 140.0
MODELS = ("classical", "Wigner", "Eckart")
NATIVE_MODULES = ("rmgpy.constants", "rmgpy.quantity", "rmgpy.species", "rmgpy.reaction",
                  "rmgpy.statmech.conformer", "rmgpy.statmech.rotation",
                  "rmgpy.statmech.translation", "rmgpy.statmech.vibration",
                  "rmgpy.kinetics.tunneling")


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate JSON key: " + key)
        result[key] = value
    return result


def read_json(path):
    return json.loads(Path(path).read_text(), object_pairs_hook=unique_object,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError("Nonfinite JSON")))


def artifact(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError("Expected a regular source file")
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def reopen(record):
    if artifact(record["path"]) != {key: record[key] for key in ("path", "sha256")}:
        raise ValueError("Source hash changed: " + record["path"])
    path = Path(record["path"])
    if "mtime_ns" in record and path.stat().st_mtime_ns != record["mtime_ns"]:
        raise ValueError("Source modification time changed: " + str(path))
    return path


def write_new(path, data):
    with Path(path).open("x") as stream:
        json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def inventory(paths):
    return {str(Path(path).resolve(strict=True)): {**artifact(path), "mtime_ns": Path(path).stat().st_mtime_ns}
            for path in sorted(set(paths))}


def reopen_inventory(records):
    for path, record in records.items():
        if path != record["path"]:
            raise ValueError("Inventory key differs from source path")
        reopen(record)


def bound_paths(value):
    """Declared records only; never recursively inventory the user's runs tree."""
    result = set()
    if isinstance(value, dict):
        if "path" in value and "sha256" in value:
            result.add(Path(value["path"]))
        for item in value.values():
            result.update(bound_paths(item))
    elif isinstance(value, list):
        for item in value:
            result.update(bound_paths(item))
    return result


def state_models():
    from storca.electronic_states import ElectronicLevel, ElectronicStateModel
    baseline = {}
    for name in POINTS:
        baseline[name] = ElectronicStateModel(
            label=name, energy_reference="retained native spin-only E0",
            reference_offset_J_mol=0.0, native_spin_multiplicity=SPINS[name],
            levels=(ElectronicLevel("native_spin_level", 0.0, SPINS[name], "retained_native_recipe"),),
            assumptions=("Retained single-surface spin-only model; not a state-resolved accuracy claim.",),
            source_identifiers=("retained_native_recipe",))
    assumptions = ("Separable two-level OH electronic and unchanged rotational/vibrational partitions.",
                   "Retained spin-free OH energy is the equal-degeneracy centroid.",
                   "TS spin-orbit quenching is assumed, not computed on this geometry.",
                   "One reactive doublet surface and thermal electronic populations are assumed.")
    centroid = {**baseline, "OH": ElectronicStateModel(
        label="OH", energy_reference="spin-free equal-degeneracy centroid",
        reference_offset_J_mol=0.0, native_spin_multiplicity=2,
        levels=(ElectronicLevel("2Pi_3/2", -DELTA_J_MOL / 2, 2, LEVEL_SOURCE),
                ElectronicLevel("2Pi_1/2", DELTA_J_MOL / 2, 2, LEVEL_SOURCE)),
        assumptions=assumptions, source_identifiers=(LEVEL_SOURCE, CENTROID_SOURCE))}
    ground = {**baseline, "OH": ElectronicStateModel(
        label="OH", energy_reference="2Pi_3/2 ground level",
        reference_offset_J_mol=-DELTA_J_MOL / 2, native_spin_multiplicity=2,
        levels=(ElectronicLevel("2Pi_3/2", 0.0, 2, LEVEL_SOURCE),
                ElectronicLevel("2Pi_1/2", DELTA_J_MOL, 2, LEVEL_SOURCE)),
        assumptions=assumptions, source_identifiers=(LEVEL_SOURCE, CENTROID_SOURCE))}
    return {key: {name: model.to_dict() for name, model in models.items()}
            for key, models in (("spin_only", baseline), ("OH_two_level_centroid", centroid),
                                ("OH_two_level_ground_reference", ground))}


def branch_specs():
    return [{"branch_id": state + "_TS_count" + str(count), "state_model": state,
             "TS_optical_isomers": count,
             "counting_assumption": ("Retained single harmonic TS basin; physical basin coverage unresolved."
                 if count == 1 else "Two independent mirror harmonic TS basins assumed; not established by geometry alone."),
             "angular_model": "retained harmonic soft mode, no rotor correction",
             "additional_channel_multiplier": 1}
            for state in ("spin_only", "OH_two_level_centroid") for count in (1, 2)]


def check_audit(audit):
    if (audit.get("status") != "passed" or audit.get("all_admitted_and_previous_artifacts_unchanged") is not True
            or audit.get("reactants") != ["OH", "CH4"] or audit.get("products") != ["H2O", "CH3"]
            or audit.get("transition_state") != "TS" or set(audit.get("stationary_points", {})) != set(POINTS)
            or audit.get("declared_reaction_degeneracy") != 1
            or audit.get("additional_degeneracy_multiplier") != 1):
        raise ValueError("Exact retained five-point, two-to-two audit required")
    grid = audit["evaluation_temperature_grid_K"]
    if grid != sorted(set(grid)) or any(type(T) not in (int, float) or not math.isfinite(T) or T <= 0 for T in grid):
        raise ValueError("Invalid evaluation grid")
    if len(audit["measurement_temperature_points"]) != 14:
        raise ValueError("All fourteen retained temperature points required")
    if any(point["temperature_K"] not in grid for point in audit["measurement_temperature_points"]):
        raise ValueError("Measurement temperature missing from evaluation grid")
    for name in POINTS:
        point = audit["stationary_points"][name]
        if (point["label"] != name or point["formula"] != FORMULAS[name]
                or point["spin_multiplicity"] != SPINS[name] or point["loaded_optical_isomers"] != 1
                or point["imposed_optical_isomers"] != 1):
            raise ValueError("Retained state/composition/counting mismatch: " + name)
    return audit


def load_reporter(work):
    report = read_json(work / "rate-comparison.json")
    path = reopen(report["report_driver"])
    if path != Path(__file__).with_name("report_oh_methane_rate_audit.py").resolve():
        raise ValueError("Unexpected retained reporter")
    spec = importlib.util.spec_from_file_location("retained_oh_rate_reporter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    closure = module.closeout(work)  # Read-only; NEVER call its file-writing report().
    if closure != report["evidence_closeout"] or report["reference"] != artifact(REFERENCE):
        raise ValueError("Retained evidence closeout/reference no longer reproduces")
    observations = read_json(REFERENCE)["audited_observation_sets"][0]["observations"]
    audits = {}
    for stage in ("baseline", "composite"):
        audit = check_audit(read_json(work / stage / "native-statmech-audit.json"))
        if module.score(audit, observations) != report["stages"][stage]["comparison"]:
            raise ValueError("Retained comparison no longer reproduces")
        audits[stage] = audit
    return report, closure, audits


def freeze_sources(work, folder, report, closure, audits):
    declared = {Path(__file__).resolve(), ROOT / "storca/electronic_states.py", REFERENCE,
                work / "rate-comparison.json", work / "contract.json", work / "energy-contract.json",
                work / "energy-result.json", work / "builder-revision-admission.json"}
    for stage in audits:
        paths = (work / (stage + "-result.json"), work / stage / "composite-rate-recipe.json",
                 work / stage / "provenance.json", work / stage / "native-statmech-audit.json")
        declared.update(paths)
        for path in paths:
            declared.update(bound_paths(read_json(path)))
    declared.update(bound_paths(report))
    declared.update(bound_paths(closure))
    declared.update(bound_paths(read_json(work / "contract.json")))
    # Installed RMG constants may be a compiled extension. Bind its readable
    # sibling too, then require actual native values to match those literals.
    declared.add(native_constants_source(audits["composite"]["native"]["modules"]["rmgpy.constants"]))
    records = inventory(declared)
    snapshots = []
    destination = folder / "source-snapshots"
    destination.mkdir()
    for index, record in enumerate(records.values()):
        original = reopen(record)
        if original.suffix != ".py":
            continue  # Native binaries are hash/version bound, not redundantly copied.
        saved = destination / (str(index) + "-" + original.name)
        with saved.open("xb") as stream:
            stream.write(original.read_bytes())
        if artifact(saved)["sha256"] != record["sha256"]:
            raise ValueError("Source snapshot differs")
        snapshots.append({"live_source": record, "snapshot": artifact(saved)})
    revision = read_json(work / "builder-revision-admission.json")
    current = reopen(revision["admitted_builder"]).read_bytes()
    delta = revision["exact_delta"]
    if current.count(delta["new_line"].encode()) != 1:
        raise ValueError("Original generator reconstruction is not uniquely defined")
    original = current.replace(delta["new_line"].encode(), delta["old_line"].encode(), 1)
    if hashlib.sha256(original).hexdigest() != revision["original_builder"]["sha256"]:
        raise ValueError("Original generator version cannot be reconstructed exactly")
    saved = destination / "original-composite-generator.py"
    with saved.open("xb") as stream:
        stream.write(original)
    historical = {"original_binding": revision["original_builder"], "snapshot": artifact(saved),
                  "admission_receipt": artifact(work / "builder-revision-admission.json")}
    reopen_inventory(records)
    return records, snapshots, historical


def build_contract(work, folder):
    report, closure, audits = load_reporter(work)
    sources, snapshots, historical = freeze_sources(work, folder, report, closure, audits)
    return {"schema_version": 1, "kind": "conditional_oh_methane_electronic_state_replay",
            "source_workdir": str(work), "read_only_evidence_closeout": closure,
            "retained_audits": {stage: artifact(work / stage / "native-statmech-audit.json") for stage in audits},
            "retained_source_inventory": sources, "source_snapshots": snapshots,
            "historic_generator_snapshot": historical, "models": state_models(), "branches": branch_specs(),
            "native_expected": {**audits["composite"]["native"],
                                "constants": read_native_constants(artifact(native_constants_source(
                                    audits["composite"]["native"]["modules"]["rmgpy.constants"])))},
            "OH_splitting_cm_inverse": 140.0, "OH_splitting_J_mol": DELTA_J_MOL,
            "conversion": "Exact-SI h, c and NA convert 140 cm^-1 to J/mol; native R is explicitly used for native thermodynamics.",
            "literature_PDF_bytes_frozen": False,
            "literature_scope": "Explicit source pointers support a declared conditional model, not computed SOC or source-audited calibration.",
            "native_replay_timeout_seconds": 1800, "experimental_reference": artifact(REFERENCE),
            "environment": "isolated ideal gas, single adiabatic doublet surface",
            "verification_status": "diagnostic_not_admitted", "network_reinjection_allowed": False,
            "physical_rate_accuracy_claim_allowed": False, "calibrated_uncertainty_claim_allowed": False,
            "lifetime_claim_allowed": False, "electronic_structure_jobs_launched": 0,
            "production_defaults_changed": False}


def native_constants_source(record):
    path = reopen(record)
    return path if path.suffix == ".py" else path.with_name("constants.py").resolve(strict=True)


def read_native_constants(record):
    """Read literal constants from the already hash-bound native source, no exec."""
    selected = {"R", "Na", "h", "kB", "c"}
    values = {}
    tree = ast.parse(reopen(record).read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in selected:
                if name in values:
                    raise ValueError("Duplicate native constant assignment")
                values[name] = finite_positive(ast.literal_eval(node.value), "native constant")
    if set(values) != selected:
        raise ValueError("Bound native constants are not the admitted literal schema")
    return values


def finite_positive(value, description):
    if type(value) not in (float, int) or not math.isfinite(value) or value <= 0:
        raise ValueError("Invalid positive finite " + description)
    return value


def identity(left, right, description, tolerance=1e-8):
    error = abs(finite_positive(left, description) / finite_positive(right, description) - 1)
    if error > tolerance:
        raise ValueError("Identity mismatch: " + description)
    return error


def reconstruct_points(audit, models, optical_count):
    from rmgpy.statmech import Conformer, IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator
    from rmgpy.species import Species, TransitionState
    from storca.electronic_states import evaluate_electronic_partition
    from rmgpy import constants

    objects = {}
    for name, row in audit["stationary_points"].items():
        modes = []
        for mode in row["modes"]:
            kind = mode["class"]
            if type(mode["quantum"]) is not bool:
                raise ValueError("Invalid mode quantum flag")
            def quantity(key, units):
                record = mode[key]
                if record["units"] != units:
                    raise ValueError("Unexpected mode units")
                values = record["value"] if isinstance(record["value"], list) else [record["value"]]
                for value in values:
                    finite_positive(value, key)
                return (record["value"], units)
            if kind == "IdealGasTranslation":
                modes.append(IdealGasTranslation(mass=quantity("mass", "amu"), quantum=mode["quantum"]))
            elif kind in ("LinearRotor", "NonlinearRotor"):
                if type(mode["symmetry"]) is not int or mode["symmetry"] < 1:
                    raise ValueError("Invalid rotor symmetry")
                cls = LinearRotor if kind == "LinearRotor" else NonlinearRotor
                modes.append(cls(inertia=quantity("inertia", "amu*angstrom^2"),
                                 symmetry=mode["symmetry"], quantum=mode["quantum"]))
            elif kind == "HarmonicOscillator":
                if len(mode["frequencies"]["value"]) != MODE_COUNTS[name]:
                    raise ValueError("Physical vibration inventory differs")
                modes.append(HarmonicOscillator(frequencies=quantity("frequencies", "cm^-1"), quantum=mode["quantum"]))
            else:
                raise ValueError("Unsupported native mode class; repr is never evaluated")
        expected_rotor = "LinearRotor" if name == "OH" else "NonlinearRotor"
        if [type(mode).__name__ for mode in modes] != ["IdealGasTranslation", expected_rotor, "HarmonicOscillator"]:
            raise ValueError("Unexpected native modal inventory")
        state = evaluate_electronic_partition(models[name], 298.15, gas_constant_J_mol_K=constants.R)
        conformer = Conformer(E0=(row["E0_J_mol"] + state["ground_offset_J_mol"], "J/mol"),
                              modes=modes, spin_multiplicity=SPINS[name],
                              optical_isomers=optical_count if name == "TS" else 1)
        if name == "TS":
            frequency = row["scaled_TS_frequency"]
            if frequency["units"] != "cm^-1" or frequency["value"] >= 0:
                raise ValueError("Invalid TS imaginary frequency")
            objects[name] = TransitionState(label=name, conformer=conformer,
                                            frequency=(frequency["value"], "cm^-1"))
        else:
            objects[name] = Species(label=name, conformer=conformer)
    return objects


def native_branch(audit, models, baseline_models, optical_count):
    from rmgpy import constants
    from rmgpy.reaction import Reaction
    from rmgpy.kinetics.tunneling import Wigner, Eckart
    from storca.electronic_states import evaluate_electronic_partition, evaluate_reaction_electronic_correction

    points = reconstruct_points(audit, models, optical_count)
    baseline = reconstruct_points(audit, baseline_models, 1)
    reactants = [points[name] for name in audit["reactants"]]
    products = [points[name] for name in audit["products"]]
    ts = points["TS"]
    E_r = sum(item.conformer.E0.value_si for item in reactants)
    E_p = sum(item.conformer.E0.value_si for item in products)
    E_ts = ts.conformer.E0.value_si
    if E_ts <= max(E_r, E_p):
        raise ValueError("Single-saddle Eckart requires two positive barriers")
    barriers = {"forward_E0": E_ts - E_r, "reverse_E0": E_ts - E_p, "reaction_E0": E_p - E_r}
    tunneling = {"classical": None, "Wigner": Wigner(frequency=(ts.frequency.value_si, "cm^-1")),
                "Eckart": Eckart(frequency=(ts.frequency.value_si, "cm^-1"), E0_reac=(E_r, "J/mol"),
                                  E0_TS=(E_ts, "J/mol"), E0_prod=(E_p, "J/mol"))}
    rates = {}
    for name, model in tunneling.items():
        ts.tunneling = model
        forward = Reaction(reactants=reactants, products=products, transition_state=ts)
        reverse = Reaction(reactants=products, products=reactants, transition_state=ts)
        rows = []
        for T in audit["evaluation_temperature_grid_K"]:
            states = {label: evaluate_electronic_partition(state, T, gas_constant_J_mol_K=constants.R)
                      for label, state in models.items()}
            correction = evaluate_reaction_electronic_correction(
                tuple(models[label] for label in audit["reactants"]), tuple(models[label] for label in audit["products"]),
                models["TS"], T, gas_constant_J_mol_K=constants.R)
            kappa = 1.0 if model is None else float(model.calculate_tunneling_factor(T))
            # Native rate has corrected ground energies/new Eckart but still
            # spin-only Q. Replace that electronic factor exactly once.
            kf = forward.calculate_tst_rate_coefficient(T) * correction["forward"]["population_factor"]
            kr = reverse.calculate_tst_rate_coefficient(T) * correction["reverse"]["population_factor"]
            volume = constants.R * T / 101325.0
            q = {label: item.get_partition_function(T) / volume * states[label]["ground_partition_relative_native_spin"]
                 for label, item in points.items()}
            Qr = math.prod(q[label] for label in audit["reactants"])
            Qp = math.prod(q[label] for label in audit["products"])
            manual_f = constants.kB * T / constants.h * q["TS"] / Qr * math.exp(-barriers["forward_E0"] / constants.R / T) * kappa
            manual_r = constants.kB * T / constants.h * q["TS"] / Qp * math.exp(-barriers["reverse_E0"] / constants.R / T) * kappa
            Kq = Qp / Qr * math.exp(-barriers["reaction_E0"] / constants.R / T)
            # Native Gibbs is spin-only. Add state ΔG relative to the ORIGINAL
            # baseline, not to a conformer whose E0 has already been shifted.
            G = {label: baseline[label].get_free_energy(T)
                 + state["thermodynamics_relative_native_baseline"]["delta_G_J_mol"]
                 for label, state in states.items() if label != "TS"}
            dG = sum(G[label] for label in audit["products"]) - sum(G[label] for label in audit["reactants"])
            Kg = math.exp(-dG / constants.R / T)
            errors = {"forward_explicit_partition_TST": identity(kf, manual_f, "forward TST"),
                      "reverse_explicit_partition_TST": identity(kr, manual_r, "reverse TST"),
                      "forward_reverse_partition_equilibrium": identity(kf / kr, Kq, "partition detailed balance"),
                      "partition_corrected_Gibbs_equilibrium": identity(Kq, Kg, "state-aware thermo equilibrium")}
            rows.append({"temperature_K": T, "forward_direct_m3_mol_s": float(kf),
                         "reverse_direct_m3_mol_s": float(kr), "tunneling_factor": kappa,
                         "corrected_partition_Kc": Kq, "corrected_thermo_Kc": Kg,
                         "ground_partition_per_standard_volume": q,
                         "native_baseline_point_G_J_mol": {label: baseline[label].get_free_energy(T)
                                                            for label in audit["reactants"] + audit["products"]},
                         "native_spin_only_thermo_Kc_after_ground_E0_shift": forward.get_equilibrium_constant(T, type="Kc"),
                         "electronic_correction": correction, "point_electronic_partitions": states,
                         "identity_relative_errors": errors})
        rates[name] = {"native_tunneling_repr": repr(model), "grid": rows}
    return {"barriers_J_mol": barriers, "rates": rates}


def validate_branches(audit, branches, reference_equivalence):
    baseline = branches["spin_only_TS_count1"]
    for name in MODELS:
        for actual, old in zip(baseline["rates"][name]["grid"], audit["rates"][name]["grid"], strict=True):
            if actual["temperature_K"] != old["temperature_K"]:
                raise ValueError("Baseline temperature order differs")
            for key in ("forward_direct_m3_mol_s", "reverse_direct_m3_mol_s", "tunneling_factor"):
                identity(actual[key], old[key], "retained baseline " + key)
    for state in ("spin_only", "OH_two_level_centroid"):
        one = branches[state + "_TS_count1"]
        two = branches[state + "_TS_count2"]
        if one["barriers_J_mol"] != two["barriers_J_mol"]:
            raise ValueError("TS counting branch changed barrier")
        for name in MODELS:
            for a, b in zip(one["rates"][name]["grid"], two["rates"][name]["grid"], strict=True):
                identity(b["forward_direct_m3_mol_s"], 2 * a["forward_direct_m3_mol_s"], "joint forward basin counting")
                identity(b["reverse_direct_m3_mol_s"], 2 * a["reverse_direct_m3_mol_s"], "joint reverse basin counting")
                identity(b["corrected_partition_Kc"], a["corrected_partition_Kc"], "counting-invariant equilibrium")
                identity(b["tunneling_factor"], a["tunneling_factor"], "counting-invariant tunneling")
                for label in POINTS:
                    identity(b["ground_partition_per_standard_volume"][label],
                             a["ground_partition_per_standard_volume"][label] * (2 if label == "TS" else 1),
                             "joint partition basin counting " + label)
    centroid = branches["OH_two_level_centroid_TS_count1"]
    if centroid["barriers_J_mol"] != reference_equivalence["barriers_J_mol"]:
        raise ValueError("Electronic energy-reference transformation changed barrier")
    for name in MODELS:
        for a, b in zip(centroid["rates"][name]["grid"], reference_equivalence["rates"][name]["grid"], strict=True):
            for key in ("forward_direct_m3_mol_s", "reverse_direct_m3_mol_s", "corrected_thermo_Kc", "tunneling_factor"):
                identity(a[key], b[key], "ground/centroid reference equivalence " + key)


def native_replay(contract_path, output):
    contract_binding = artifact(contract_path)
    contract = read_json(contract_path)
    reopen_inventory(contract["retained_source_inventory"])
    for snapshot in contract["source_snapshots"]:
        reopen(snapshot["snapshot"])
    reopen(contract["historic_generator_snapshot"]["snapshot"])
    expected = contract["native_expected"]
    if artifact(sys.executable) != expected["interpreter"]:
        raise ValueError("Selected RMG interpreter differs from retained runtime")
    native_modules = {name: artifact(importlib.import_module(name).__file__) for name in NATIVE_MODULES}
    if any(record != expected["modules"][name] for name, record in native_modules.items()):
        raise ValueError("Selected RMG native module differs from retained runtime")
    from rmgpy import constants
    from storca.electronic_states import ElectronicStateModel
    models = {key: {label: ElectronicStateModel.from_dict(value) for label, value in row.items()}
              for key, row in contract["models"].items()}
    results = {}
    for stage, record in contract["retained_audits"].items():
        audit = check_audit(read_json(reopen(record)))
        branches = {spec["branch_id"]: native_branch(audit, models[spec["state_model"]], models["spin_only"],
                                                    spec["TS_optical_isomers"])
                    for spec in contract["branches"]}
        equivalence = native_branch(audit, models["OH_two_level_ground_reference"], models["spin_only"], 1)
        validate_branches(audit, branches, equivalence)
        results[stage] = {"branches": branches, "ground_centroid_reference_equivalence": "passed",
                          "spin_only_baseline_reproduction": "passed",
                          "joint_TS_counting_invariants": "passed"}
    reopen(contract_binding)
    reopen_inventory(contract["retained_source_inventory"])
    document = {"schema_version": 1, "kind": "native_oh_methane_state_replay", "status": "passed",
                "contract": contract_binding, "stages": results,
                "native": {"interpreter": artifact(sys.executable), "modules": native_modules,
                           "constants": {name: float(getattr(constants, name)) for name in ("R", "Na", "h", "kB", "c")}},
                "electronic_structure_jobs_launched": 0, "retained_artifacts_unchanged": True,
                "experimental_rates_read_into_physical_model": False,
                "network_reinjection_allowed": False, "production_defaults_changed": False,
                "verification_status": "diagnostic_not_admitted"}
    write_new(output, document)
    print("NATIVE_ELECTRONIC_STATE_REPLAY_PASSED", flush=True)


def validate_native_receipt(contract, native):
    """Recompute semantic witnesses; a success flag or rehashed k is not enough."""
    from storca.electronic_states import ElectronicStateModel, evaluate_electronic_partition, evaluate_reaction_electronic_correction
    if (native.get("status") != "passed" or native.get("verification_status") != "diagnostic_not_admitted"
            or native.get("retained_artifacts_unchanged") is not True
            or native.get("electronic_structure_jobs_launched") != 0
            or native.get("network_reinjection_allowed") is not False
            or native.get("production_defaults_changed") is not False
            or set(native.get("stages", {})) != set(contract["retained_audits"])):
        raise ValueError("Invalid diagnostic-only native disposition")
    runtime = native["native"]
    expected_runtime = contract["native_expected"]
    if (runtime["interpreter"] != expected_runtime["interpreter"]
            or runtime["modules"] != {name: expected_runtime["modules"][name] for name in NATIVE_MODULES}
            or runtime["constants"] != expected_runtime["constants"]):
        raise ValueError("Native receipt runtime/constants differ from source-bound contract")
    reopen(runtime["interpreter"])
    for record in runtime["modules"].values():
        reopen(record)
    constants = runtime["constants"]
    R, kB, h = (finite_positive(constants[key], "native constant") for key in ("R", "kB", "h"))
    models = {key: {label: ElectronicStateModel.from_dict(value) for label, value in values.items()}
              for key, values in contract["models"].items()}
    for stage, result in native["stages"].items():
        audit = check_audit(read_json(reopen(contract["retained_audits"][stage])))
        specs = {item["branch_id"]: item for item in contract["branches"]}
        if set(result["branches"]) != set(specs):
            raise ValueError("Missing or extra counting/state branch")
        baseline = result["branches"]["spin_only_TS_count1"]
        for branch_id, branch in result["branches"].items():
            spec = specs[branch_id]
            states = models[spec["state_model"]]
            for name in MODELS:
                rows = branch["rates"][name]["grid"]
                old_rows = audit["rates"][name]["grid"]
                base_rows = baseline["rates"][name]["grid"]
                if [row["temperature_K"] for row in rows] != audit["evaluation_temperature_grid_K"]:
                    raise ValueError("Native replay omitted/duplicated/reordered evaluation temperatures")
                for row, base, old in zip(rows, base_rows, old_rows, strict=True):
                    T = row["temperature_K"]
                    evaluated = {label: evaluate_electronic_partition(model, T, gas_constant_J_mol_K=R)
                                 for label, model in states.items()}
                    expected = evaluate_reaction_electronic_correction(
                        tuple(states[label] for label in audit["reactants"]),
                        tuple(states[label] for label in audit["products"]), states["TS"], T,
                        gas_constant_J_mol_K=R)
                    if evaluated != row["point_electronic_partitions"] or expected != row["electronic_correction"]:
                        raise ValueError("Electronic partition/correction witness changed")
                    bars = branch["barriers_J_mol"]
                    for direction in ("forward", "reverse"):
                        target = audit["barriers_J_mol"][direction + "_E0"] + expected[direction]["barrier_shift_J_mol"]
                        if abs(bars[direction + "_E0"] - target) > 1e-6:
                            raise ValueError("State-aware barrier witness changed")
                    if abs(bars["forward_E0"] - bars["reverse_E0"] - bars["reaction_E0"]) > 1e-6:
                        raise ValueError("Reaction-energy witness changed")
                    q = row["ground_partition_per_standard_volume"]
                    for label in POINTS:
                        factor = evaluated[label]["ground_partition_relative_native_spin"]
                        if label == "TS":
                            factor *= spec["TS_optical_isomers"]
                        identity(q[label], base["ground_partition_per_standard_volume"][label] * factor,
                                 "state-aware partition witness " + label)
                    Qr = math.prod(q[label] for label in audit["reactants"])
                    Qp = math.prod(q[label] for label in audit["products"])
                    for direction, endpoint_q in (("forward", Qr), ("reverse", Qp)):
                        explicit = kB * T / h * q["TS"] / endpoint_q * math.exp(-bars[direction + "_E0"] / R / T) * row["tunneling_factor"]
                        identity(row[direction + "_direct_m3_mol_s"], explicit, "reported rate witness")
                        if branch_id == "spin_only_TS_count1":
                            identity(row[direction + "_direct_m3_mol_s"], old[direction + "_direct_m3_mol_s"], "retained rate witness")
                    Kq = Qp / Qr * math.exp(-bars["reaction_E0"] / R / T)
                    identity(row["corrected_partition_Kc"], Kq, "reported partition equilibrium witness")
                    G = {label: base["native_baseline_point_G_J_mol"][label]
                         + evaluated[label]["thermodynamics_relative_native_baseline"]["delta_G_J_mol"]
                         for label in audit["reactants"] + audit["products"]}
                    Kg = math.exp(-(sum(G[label] for label in audit["products"]) - sum(G[label] for label in audit["reactants"])) / R / T)
                    identity(row["corrected_thermo_Kc"], Kg, "reported thermo witness")
                    identity(Kq, Kg, "corrected equilibrium witness")
                    identity(row["forward_direct_m3_mol_s"] / row["reverse_direct_m3_mol_s"], Kq, "reported detailed balance")
                    expected_rate_ratio = expected["forward"]["classical_rate_factor"] * spec["TS_optical_isomers"] * row["tunneling_factor"] / base["tunneling_factor"]
                    identity(row["forward_direct_m3_mol_s"], base["forward_direct_m3_mol_s"] * expected_rate_ratio,
                             "joint state/Eckart rate ratio")
                    if name == "classical":
                        identity(row["tunneling_factor"], 1.0, "classical transmission")
                    elif name == "Wigner" or branch_id == "spin_only_TS_count1":
                        identity(row["tunneling_factor"], old["tunneling_factor"], "baseline/frequency-only transmission")


def compare_measurements(native):
    frozen = read_json(REFERENCE)["audited_observation_sets"][0]
    if not frozen["transcription_frozen"] or len(frozen["observations"]) != 14:
        raise ValueError("Exact frozen fourteen-row reference required")
    comparison = {}
    for stage, result in native["stages"].items():
        comparison[stage] = {}
        for branch, values in result["branches"].items():
            models = {}
            for name in MODELS:
                by_T = {row["temperature_K"]: row for row in values["rates"][name]["grid"]}
                rows = []
                for point in frozen["observations"]:
                    if point["rate_units"] != "cm^3 molecule^-1 s^-1":
                        raise ValueError("Unexpected experimental units")
                    predicted = finite_positive(by_T[point["temperature_K"]]["forward_direct_m3_mol_s"], "rate")
                    measured = finite_positive(point["rate"] * 6.02214076e23 * 1e-6, "reference")
                    rows.append({"point_id": point["point_id"], "temperature_K": point["temperature_K"],
                                 "forward_direct_m3_mol_s": predicted, "measured_m3_mol_s": measured,
                                 "prediction_over_measurement": predicted / measured})
                models[name] = {"rows": rows, "ratio_min": min(row["prediction_over_measurement"] for row in rows),
                                "ratio_max": max(row["prediction_over_measurement"] for row in rows)}
            comparison[stage][branch] = models
    return comparison


def run(work, rmg_env):
    from src.orca_runner import run_bounded_kinetics_command
    work = Path(work).resolve(strict=True)
    if ROOT / "runs" not in work.parents:
        raise ValueError("Retained source must be an owned repository runs child")
    # Preflight before creating an experiment or invoking the native runtime.
    load_reporter(work)
    folder = Path(tempfile.mkdtemp(prefix="oh-methane-state-replay-20261011-", dir=ROOT / "runs"))
    print("Fresh diagnostic: " + str(folder), flush=True)
    started = time.monotonic()
    try:
        contract = build_contract(work, folder)
        write_new(folder / "contract.json", contract)
        binding = artifact(folder / "contract.json")
        command = ["conda", "run", "--no-capture-output", "-n", rmg_env, "python", "-B", str(Path(__file__).resolve()),
                   "--native", str(folder / "contract.json"), str(folder / "native-replay.json")]
        completed = run_bounded_kinetics_command(command, cwd=folder, timeout_seconds=contract["native_replay_timeout_seconds"])
        write_new(folder / "execution.json", {"command": command, "returncode": completed.returncode,
                  "stdout": completed.stdout, "stderr": completed.stderr, "elapsed_seconds": time.monotonic() - started})
        if completed.returncode != 0:
            raise ValueError("Native electronic-state replay failed; execution streams retained")
        native = read_json(folder / "native-replay.json")
        if native["status"] != "passed" or native["contract"] != binding:
            raise ValueError("Native receipt differs from frozen contract")
        validate_native_receipt(contract, native)
        reopen(binding)
        reopen_inventory(contract["retained_source_inventory"])
        load_reporter(work)
        write_new(folder / "result.json", {"schema_version": 1, "kind": contract["kind"],
                  "status": "completed_conditional_state_diagnostic", "verification_status": "diagnostic_not_admitted",
                  "contract": binding, "native_replay": artifact(folder / "native-replay.json"),
                  "execution": artifact(folder / "execution.json"), "comparison": compare_measurements(native),
                  "retained_artifacts_unchanged": True, "experimental_studies": 1, "measurement_rows": 14,
                  "electronic_structure_jobs_launched": 0, "network_reinjection_allowed": False,
                  "production_defaults_changed": False, "physical_rate_accuracy_claim_allowed": False,
                  "calibrated_uncertainty_claim_allowed": False, "lifetime_claim_allowed": False,
                  "limitations": ["Electronic-level model and TS SOC quenching are assumptions, not computed SOC.",
                    "TS-count1/count2 are joint harmonic-basin alternatives; neither is selected as the true rate.",
                    "Soft-mode harmonic validity, tunneling dynamics and finite-basis error remain unresolved.",
                    "Humid-air experiment versus isolated-gas applicability remains conditional.",
                    "No mechanism-completeness, lifetime, environment-support or calibrated-confidence claim."]})
    except Exception as error:
        write_new(folder / "failure.json", {"status": "failed", "failure_type": type(error).__name__,
                  "failure_reason": str(error), "elapsed_seconds": time.monotonic() - started,
                  "stdout": getattr(error, "stdout", getattr(error, "output", None)),
                  "stderr": getattr(error, "stderr", None), "electronic_structure_jobs_launched": 0,
                  "network_reinjection_allowed": False})
        raise
    print("CONDITIONAL_STATE_DIAGNOSTIC_COMPLETED " + str(folder / "result.json"), flush=True)
    return folder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("native_output", type=Path, nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("--rmg-env", default="rmg4_env")
    parser.add_argument("--native", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.native:
        if args.native_output is None or args.native_output.exists() or args.native_output.is_symlink():
            raise ValueError("Fresh native diagnostic output required")
        if args.native_output.parent.resolve() != args.source.parent.resolve():
            raise ValueError("Native output must be in the fresh contract directory")
        native_replay(args.source, args.native_output)
    else:
        if args.native_output is not None:
            parser.error("Only the retained source directory is a public positional argument")
        run(args.source, args.rmg_env)


if __name__ == "__main__":
    main()
