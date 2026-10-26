"""Portable, source-bound admission for STORCA's ORCA rotor adapter.

This module deliberately imports no chemistry packages.  A receipt proves the
declared stationary points passed one particular statmech loading policy; it
does not establish rate accuracy or pressure/conformer applicability.
"""

from __future__ import annotations

import ast
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import uuid

POLICY_ID = "storca-orca-diatomic-rotor-v1"
ATOMIC_ENERGY_POLICY = "one_electron_hydrogen_native_hf_same_basis_v1"
CONTRACT_NAME = "arkane-rotor-contract.json"
RECEIPT_NAME = "arkane-statmech-validation.json"
FAILURE_NAME = "arkane-statmech-validation.failed.json"
# A refusal boundary, never permission to flatten a polyatomic geometry.
NEAR_LINEAR_RATIO = 1.0e-6
COLLINEAR_ROUNDOFF_RATIO = 1.0e-12
SYMBOLS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
           "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr "
           "Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr "
           "Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl "
           "Mc Lv Ts Og").split()


def artifact(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError("Arkane rotor artifact is not a file")
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def reopen(record):
    if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
        raise ValueError("Arkane rotor artifact record is invalid")
    try:
        actual = artifact(record["path"])
    except (OSError, TypeError) as error:
        raise ValueError("Arkane rotor artifact is absent") from error
    if actual != record:
        raise ValueError("Arkane rotor artifact hash mismatch")
    return Path(record["path"])


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def _unique(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate Arkane rotor JSON key")
        result[key] = value
    return result


def read_document(path):
    try:
        value = json.loads(Path(path).read_text(), object_pairs_hook=_unique)
    except (OSError, TypeError, json.JSONDecodeError) as error:
        raise ValueError("Arkane rotor document is absent or invalid") from error
    if not isinstance(value, dict):
        raise ValueError("Arkane rotor document must be an object")
    return value


def _write_new(path, value):
    with Path(path).open("x") as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _positive(value, field):
    if isinstance(value, bool) or not isinstance(value, (int, float)) \
            or not math.isfinite(value) or value <= 0:
        raise ValueError("Invalid " + field)
    return float(value)


def _geometry(coords, numbers, masses):
    try:
        if any(isinstance(value, bool) for value in masses):
            raise ValueError("Invalid native atom mass")
        xyz = [[float(item) for item in row] for row in coords]
        nums = list(numbers)
        mass = [_positive(float(value), "native atom mass") for value in masses]
    except (TypeError, OverflowError) as error:
        raise ValueError("Malformed native geometry") from error
    if not xyz or len(xyz) != len(nums) or len(xyz) != len(mass) \
            or any(len(row) != 3 or not all(math.isfinite(v) for v in row) for row in xyz):
        raise ValueError("Malformed or nonfinite native geometry")
    for value in nums:
        try:
            valid = not isinstance(value, bool) and int(value) == value and 1 <= value < len(SYMBOLS)
        except (TypeError, ValueError, OverflowError):
            valid = False
        if not valid:
            raise ValueError("Invalid native atomic number")
    return xyz, [int(value) for value in nums], mass


def diatomic_inertia(coords, masses):
    """Analytic transverse inertia using native amu masses and angstroms."""
    xyz, _, mass = _geometry(coords, [1, 1], masses)
    if len(xyz) != 2:
        raise ValueError("Diatomic rotor requires exactly two atoms")
    # Explicit summation is portable between the selected RMG interpreter and
    # the application interpreter (built-in sum changed in Python 3.12).
    distance2 = math.fsum((b - a) ** 2 for a, b in zip(xyz[0], xyz[1]))
    reduced_mass = mass[0] / (mass[0] + mass[1]) * mass[1]
    return _positive(reduced_mass * distance2, "diatomic transverse inertia")


def assess_geometry(coords, numbers, masses, *, allow_linear=False):
    """Admit atoms/diatomics and clearly bent polyatomics; refuse ambiguity."""
    xyz, nums, mass = _geometry(coords, numbers, masses)
    count = len(xyz)
    result = {"atom_count": count, "coordinates_angstrom": xyz,
              "atomic_numbers": nums, "atomic_masses_amu": mass}
    if count == 1:
        return {**result, "classification": "atom"}
    pairs = [(math.fsum((b - a) ** 2 for a, b in zip(xyz[i], xyz[j])), i, j)
             for i in range(count) for j in range(i + 1, count)]
    if any(not math.isfinite(value) or value <= 0 for value, _, _ in pairs):
        raise ValueError("Coincident or invalid native atoms")
    span2, first, last = max(pairs)
    if count == 2:
        return {**result, "classification": "diatomic", "bond_length_squared_angstrom2": span2,
                "analytic_inertia_amu_angstrom2": diatomic_inertia(xyz, mass)}
    direction = [b - a for a, b in zip(xyz[first], xyz[last])]
    transverse2 = []
    for point in xyz:
        displacement = [b - a for a, b in zip(xyz[first], point)]
        projection = math.fsum(a * b for a, b in zip(displacement, direction)) / span2
        transverse2.append(math.fsum((a - projection * b) ** 2 for a, b in zip(displacement, direction)))
    ratio = math.sqrt(max(transverse2) / span2)
    if allow_linear and math.isfinite(ratio) and ratio <= COLLINEAR_ROUNDOFF_RATIO:
        return {**result, "classification": "linear_polyatomic", "linearity_ratio": ratio}
    if not math.isfinite(ratio) or ratio <= NEAR_LINEAR_RATIO:
        raise ValueError("Polyatomic external rotor is linear or near-linear; classification unresolved")
    return {**result, "classification": "nonlinear_polyatomic", "linearity_ratio": ratio}


def _assignments(path):
    tree = ast.parse(Path(path).read_text())
    terms = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1 \
                or not isinstance(node.targets[0], ast.Name) or node.targets[0].id in terms:
            raise ValueError("Unexpected Arkane species declaration")
        terms[node.targets[0].id] = node.value
    return tree, terms


def _stationary_points(input_path, provenance):
    if provenance.get("kind") != "verified_arkane_route" \
            or provenance.get("arkane_statmech_policy") != POLICY_ID:
        raise ValueError("Arkane rotor policy is absent from verified provenance")
    records = provenance.get("orca_artifacts")
    if not isinstance(records, list) or not records:
        raise ValueError("Arkane rotor stationary-point inventory is absent")
    tree = ast.parse(input_path.read_text())
    declarations = {}
    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call) \
                and isinstance(node.value.func, ast.Name) \
                and node.value.func.id in {"species", "transitionState"} and len(node.value.args) >= 2:
            call = node.value
            if len(call.args) != 2:
                raise ValueError("Unexpected stationary-point declaration")
            label = ast.literal_eval(call.args[0])
            if not isinstance(label, str) or not label or label in declarations:
                raise ValueError("Duplicate or invalid stationary-point label")
            declarations[label] = (call.func.id, (input_path.parent / ast.literal_eval(call.args[1])).resolve())
    allowed_logs = {str(Path(item["staged_output"]).resolve()) for item in records}
    result = []
    for record in records:
        label, kind = record.get("label"), record.get("kind")
        if kind not in {"species", "transition_state"} or label not in declarations:
            raise ValueError("Stationary point is absent from the Arkane input")
        statmech = Path(record["statmech_file"]).resolve(strict=True)
        if declarations[label] != ("species" if kind == "species" else "transitionState", statmech):
            raise ValueError("Stationary-point input/provenance mismatch")
        formula = record.get("formula")
        if not isinstance(formula, dict) or not formula \
                or any(symbol not in SYMBOLS[1:] or isinstance(count, bool)
                       or not isinstance(count, int) or count <= 0 for symbol, count in formula.items()):
            raise ValueError("Stationary-point atom composition is invalid")
        if isinstance(record.get("multiplicity"), bool) or not isinstance(record.get("multiplicity"), int) \
                or record["multiplicity"] <= 0:
            raise ValueError("Stationary-point spin multiplicity is invalid")
        atomic = record.get("evidence_kind") == "monatomic_single_point"
        if atomic and (formula != {"H": 1} or kind != "species" or record.get("charge") != 0
                       or record["multiplicity"] != 2 or record.get("isotope") not in (0, 1)
                       or record.get("atomic_energy_policy") != ATOMIC_ENERGY_POLICY
                       or not math.isclose(_positive(record.get("atomic_mass_amu"), "atomic mass"),
                                           1.00782503224, rel_tol=0, abs_tol=1e-8)
                       or "hessian" in record or "staged_hessian" in record):
            raise ValueError("Unsupported explicit atomic statmech evidence")
        source_names = ("output", "staged_output") if atomic else ("output", "hessian", "staged_output", "staged_hessian")
        sources = {name: artifact(record[name]) for name in source_names}
        source_checks = (("output", "output"), ("staged_output", "output")) if atomic else (
            ("output", "output"), ("staged_output", "output"), ("hessian", "hessian"), ("staged_hessian", "hessian"))
        for name, key in source_checks:
            if sources[name]["sha256"] != (record.get("checksums") or {}).get(key):
                raise ValueError("Stationary-point source hash mismatch")
        if not atomic and Path(sources["staged_output"]["path"]).with_suffix(".hess") != Path(sources["staged_hessian"]["path"]):
            raise ValueError("Arkane Hessian is not colocated with its admitted log")
        species_tree, terms = _assignments(statmech)
        if not {"energy", "geometry", "frequencies", "spin_multiplicity", "opticalIsomers"} <= terms.keys():
            raise ValueError("Incomplete Arkane species declarations")
        if ast.literal_eval(terms["spin_multiplicity"]) != record["multiplicity"]:
            raise ValueError("Arkane species spin/provenance mismatch")
        if isinstance(ast.literal_eval(terms["opticalIsomers"]), bool) \
                or not isinstance(ast.literal_eval(terms["opticalIsomers"]), int) \
                or ast.literal_eval(terms["opticalIsomers"]) < 1:
            raise ValueError("Arkane optical-isomer declaration is invalid")
        for field in ("geometry", "frequencies"):
            value = terms[field]
            if not isinstance(value, ast.Call) or not isinstance(value.func, ast.Name) \
                    or value.func.id != "Log" or len(value.args) != 1 or value.keywords \
                    or str(Path(ast.literal_eval(value.args[0])).resolve()) != sources["staged_output"]["path"]:
                raise ValueError("Arkane geometry/frequency source is not admitted")
        for call in (node for node in ast.walk(species_tree) if isinstance(node, ast.Call)):
            if isinstance(call.func, ast.Name) and call.func.id == "Log" and len(call.args) == 1 \
                    and not call.keywords and str(Path(ast.literal_eval(call.args[0])).resolve()) in allowed_logs:
                continue
            if isinstance(call.func, ast.Attribute) and call.func.attr == "load_energy" \
                    and isinstance(call.func.value, ast.Call) and isinstance(call.func.value.func, ast.Name) \
                    and call.func.value.func.id == "Log" and not call.args and not call.keywords:
                continue
            raise ValueError("Arkane species file uses an unknown log or unsupported callable")
        item = {"label": label, "kind": kind, "atom_count": sum(formula.values()),
                       "formula": formula, "charge": record.get("charge"),
                       "multiplicity": record["multiplicity"], "optical_isomers": ast.literal_eval(terms["opticalIsomers"]),
                       "statmech_file": artifact(statmech),
                       "sources": sources}
        if atomic:
            if ast.literal_eval(terms.get("externalSymmetry")) != 1 or item["optical_isomers"] != 1:
                raise ValueError("Atomic symmetry and optical counts must both be one")
            item.update({"evidence_kind": "monatomic_single_point", "atomic_mass_amu": record["atomic_mass_amu"],
                         "isotope": record["isotope"], "atomic_energy_policy": record["atomic_energy_policy"]})
        result.append(item)
    if len(result) != len(declarations) or len({item["label"] for item in result}) != len(result):
        raise ValueError("Arkane stationary-point inventory is incomplete or duplicated")
    return result


def create_rotor_contract(workdir: Path) -> Path:
    workdir = Path(workdir).resolve()
    for name in (CONTRACT_NAME, RECEIPT_NAME, FAILURE_NAME, "output.py"):
        if (workdir / name).exists() or (workdir / name).is_symlink():
            raise ValueError("Arkane rotor execution requires a fresh workspace: " + name)
    input_path = workdir / "input.py"
    provenance_path = workdir / "provenance.json"
    contract = {"schema_version": 1, "kind": "storca_arkane_rotor_contract", "policy_id": POLICY_ID,
                "run_id": str(uuid.uuid4()), "input": artifact(input_path),
                "provenance": artifact(provenance_path), "expected_output": str(workdir / "output.py"),
                "adapter_sources": {name: artifact(Path(__file__).with_name(name)) for name in
                                    ("arkane_rotor.py", "arkane_launcher.py")},
                "stationary_points": _stationary_points(input_path, read_document(provenance_path))}
    contract["route_id"] = (read_document(provenance_path).get("route") or {}).get("label")
    contract["scientific_hash"] = _hash(contract)
    path = workdir / CONTRACT_NAME
    _write_new(path, contract)
    return path


def validate_contract(contract):
    contract = read_document(contract) if not isinstance(contract, dict) else contract
    plain = {key: value for key, value in contract.items() if key != "scientific_hash"}
    if contract.get("schema_version") != 1 or contract.get("kind") != "storca_arkane_rotor_contract" \
            or contract.get("policy_id") != POLICY_ID or not isinstance(contract.get("run_id"), str) \
            or not contract["run_id"] or contract.get("scientific_hash") != _hash(plain):
        raise ValueError("Arkane rotor contract is invalid")
    input_path, provenance_path = reopen(contract["input"]), reopen(contract["provenance"])
    if input_path.parent != provenance_path.parent \
            or Path(contract["expected_output"]) != input_path.parent / "output.py":
        raise ValueError("Arkane rotor workspace binding is invalid")
    if contract.get("stationary_points") != _stationary_points(input_path, read_document(provenance_path)):
        raise ValueError("Arkane rotor stationary-point contract changed")
    if contract.get("route_id") != (read_document(provenance_path).get("route") or {}).get("label"):
        raise ValueError("Arkane rotor contract reaction identity changed")
    sources = contract.get("adapter_sources")
    if not isinstance(sources, dict) or set(sources) != {"arkane_rotor.py", "arkane_launcher.py"}:
        raise ValueError("Arkane rotor adapter source inventory is incomplete")
    for source in sources.values():
        reopen(source)
    return contract


def _mode_snapshot(call):
    if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Name) or call.args:
        raise ValueError("Unexpected retained statmech mode")
    kwargs = {item.arg: item.value for item in call.keywords}
    if None in kwargs or len(kwargs) != len(call.keywords):
        raise ValueError("Unexpected retained statmech mode keywords")
    name = call.func.id
    if name not in {"LinearRotor", "NonlinearRotor"}:
        return name, None
    inertia, units = ast.literal_eval(kwargs["inertia"])
    if units != "amu*angstrom^2":
        raise ValueError("Retained rotor inertia units are invalid")
    values = [inertia] if isinstance(inertia, (float, int)) else list(inertia)
    return name, {"class": name, "inertia_amu_angstrom2": [float(value) for value in values],
                  "symmetry": ast.literal_eval(kwargs["symmetry"]),
                  "quantum": ast.literal_eval(kwargs["quantum"]) if "quantum" in kwargs else False}


def _validate_output(path, points, audits):
    calls = [node.value for node in ast.parse(path.read_text()).body if isinstance(node, ast.Expr)
             and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name)
             and node.value.func.id == "conformer"]
    output = {}
    for call in calls:
        kwargs = {item.arg: item.value for item in call.keywords}
        label = ast.literal_eval(kwargs["label"])
        if label in output or not isinstance(kwargs.get("modes"), (ast.List, ast.Tuple)):
            raise ValueError("Duplicate or incomplete retained conformer")
        modes = [_mode_snapshot(mode) for mode in kwargs["modes"].elts]
        output[label] = {"mode_types": [item[0] for item in modes],
                         "rotors": [item[1] for item in modes if item[1] is not None],
                         "spin_multiplicity": ast.literal_eval(kwargs["spin_multiplicity"]),
                         "optical_isomers": ast.literal_eval(kwargs["optical_isomers"])}
    expected = {point["label"] for point in points}
    if set(output) != expected or not isinstance(audits, list) \
            or len(audits) != len(expected) or {row.get("label") for row in audits} != expected:
        raise ValueError("Arkane rotor receipt/output stationary-point inventory is incomplete")
    by_label = {row["label"]: row for row in audits}
    for point in points:
        row, saved = by_label[point["label"]], output[point["label"]]
        if row.get("status") != "passed" or row.get("source_path") != point["sources"]["staged_output"]["path"] \
                or (row.get("geometry") or {}).get("atom_count") != point["atom_count"] \
                or saved["spin_multiplicity"] != point["multiplicity"] \
                or saved["optical_isomers"] != point["optical_isomers"] \
                or saved["optical_isomers"] != row.get("optical_isomers") \
                or set(row.get("preserved") or {}) != {"nonrotor_modes", "energy", "unscaled_frequencies",
                                                       "spin_multiplicity", "optical_isomers"} \
                or not all(value is True for value in row["preserved"].values()):
            raise ValueError("Arkane rotor audit has inconsistent state or preservation evidence")
        after = row.get("after") or {}
        geometry = row["geometry"]
        recomputed = assess_geometry(geometry.get("coordinates_angstrom"), geometry.get("atomic_numbers"),
                                     geometry.get("atomic_masses_amu"), allow_linear=True)
        if recomputed != geometry or dict(Counter(SYMBOLS[number] for number in geometry["atomic_numbers"])) != point["formula"]:
            raise ValueError("Arkane rotor audit geometry/composition is inconsistent")
        if after.get("mode_types") != saved["mode_types"]:
            raise ValueError("Arkane output mode inventory differs from rotor audit")
        before_modes = (row.get("before") or {}).get("mode_types")
        if not isinstance(before_modes, list) or [mode for mode in before_modes
                if mode not in {"LinearRotor", "NonlinearRotor"}] != [mode for mode in saved["mode_types"]
                if mode not in {"LinearRotor", "NonlinearRotor"}]:
            raise ValueError("Arkane rotor audit changed its nonrotational mode inventory")
        if point["atom_count"] == 1:
            if saved["rotors"] or after.get("rotor") is not None or row.get("action") != "unchanged_atom" \
                    or row.get("before") != after or saved["mode_types"] != ["IdealGasTranslation"]:
                raise ValueError("An atom must not have an external rotor")
            continue
        if len(saved["rotors"]) != 1 or not isinstance(after.get("rotor"), dict):
            raise ValueError("Arkane external rotor inventory is invalid")
        native, retained = after["rotor"], saved["rotors"][0]
        if native.get("class") != retained["class"] or native.get("symmetry") != retained["symmetry"] \
                or native.get("quantum") is not retained["quantum"] \
                or not isinstance(native.get("symmetry"), int) or isinstance(native["symmetry"], bool) \
                or native["symmetry"] < 1 or not isinstance(native.get("quantum"), bool):
            raise ValueError("Arkane retained rotor metadata differs from the native audit")
        expected_class = "LinearRotor" if geometry["classification"] in {"diatomic", "linear_polyatomic"} else "NonlinearRotor"
        if native["class"] != expected_class:
            raise ValueError("Arkane retained external rotor classification is invalid")
        values = native.get("inertia_amu_angstrom2")
        size = 1 if expected_class == "LinearRotor" else 3
        if not isinstance(values, list) or len(values) != size or len(retained["inertia_amu_angstrom2"]) != size:
            raise ValueError("Arkane rotor inertia dimensions are invalid")
        serialized = native.get("serialized_inertia_amu_angstrom2")
        if not isinstance(serialized, list) or len(serialized) != size:
            raise ValueError("Arkane native rotor serialization evidence is absent")
        for left, right, printed in zip(values, retained["inertia_amu_angstrom2"], serialized):
            _positive(left, "audited rotor inertia")
            _positive(right, "retained rotor inertia")
            # Installed RMG's quantity repr uses {:g}, i.e. six significant
            # digits. Compare that exact serialization, not a model tolerance;
            # full native precision remains authoritative in the audit.
            if right != printed or printed not in (left, float(format(left, "g"))):
                raise ValueError("Arkane retained rotor inertia differs from the native audit")
        if point["atom_count"] == 2:
            analytic = _positive(geometry.get("analytic_inertia_amu_angstrom2"), "analytic diatomic inertia")
            previous_rotor = (row.get("before") or {}).get("rotor") or {}
            if geometry.get("classification") != "diatomic" or row.get("action") != "diatomic_linear_rotor" \
                    or previous_rotor.get("symmetry") != native["symmetry"] \
                    or previous_rotor.get("quantum") is not native["quantum"] \
                    or previous_rotor.get("class") not in {"LinearRotor", "NonlinearRotor"} \
                    or not math.isclose(values[0], analytic, rel_tol=1e-12, abs_tol=0):
                raise ValueError("Arkane diatomic inertia is not its admitted analytic inertia")
        elif geometry["classification"] == "linear_polyatomic":
            if geometry["linearity_ratio"] > COLLINEAR_ROUNDOFF_RATIO \
                    or row.get("action") != "unchanged_linear_polyatomic" or row.get("before") != after:
                raise ValueError("Arkane linear polyatomic rotor was not preserved")
            moments = row.get("native_principal_inertia_amu_angstrom2")
            if not isinstance(moments, list) or len(moments) != 3 \
                    or not all(isinstance(value, (float, int)) and math.isfinite(value) for value in moments):
                raise ValueError("Arkane linear polyatomic inertia evidence is incomplete")
            smallest, first, second = sorted(moments)
            if abs(smallest) > COLLINEAR_ROUNDOFF_RATIO * _positive(first, "polyatomic transverse inertia") \
                    or not math.isclose(first, _positive(second, "polyatomic transverse inertia"), rel_tol=1e-12, abs_tol=0) \
                    or not math.isclose(values[0], first, rel_tol=1e-12, abs_tol=0):
                raise ValueError("Arkane linear polyatomic transverse moments are inconsistent")
        elif geometry.get("classification") != "nonlinear_polyatomic" \
                or _positive(geometry.get("linearity_ratio"), "polyatomic bend ratio") <= NEAR_LINEAR_RATIO \
                or row.get("action") != "unchanged_nonlinear_polyatomic" or row.get("before") != after:
            raise ValueError("Arkane polyatomic rotor classification is unresolved")


def mint_rotor_receipt(workdir, contract, audits, native_metadata, *, execution_scope="native_arkane"):
    """Finalize a successful native run only after its output is complete."""
    workdir = Path(workdir).resolve()
    contract = validate_contract(contract)
    if Path(contract["expected_output"]).parent != workdir:
        raise ValueError("Arkane rotor receipt workspace differs from its contract")
    _validate_output(workdir / "output.py", contract["stationary_points"], audits)
    receipt = {"schema_version": 1, "kind": "storca_arkane_statmech_validation", "policy_id": POLICY_ID,
               "status": "passed", "execution_scope": execution_scope,
               "route_id": contract["route_id"],
               "contract": artifact(workdir / CONTRACT_NAME), "input": contract["input"],
               "provenance": contract["provenance"], "output": artifact(workdir / "output.py"),
               "adapter_sources": contract["adapter_sources"], "native": native_metadata,
               "stationary_points": audits}
    receipt["scientific_hash"] = _hash(receipt)
    _validate_receipt_document(receipt)
    _write_new(workdir / RECEIPT_NAME, receipt)
    return validate_rotor_receipt(workdir)


def _validate_receipt_document(receipt):
    plain = {key: value for key, value in receipt.items() if key != "scientific_hash"}
    if receipt.get("schema_version") != 1 or receipt.get("kind") != "storca_arkane_statmech_validation" \
            or receipt.get("policy_id") != POLICY_ID or receipt.get("status") != "passed" \
            or receipt.get("execution_scope") != "native_arkane" or receipt.get("scientific_hash") != _hash(plain):
        raise ValueError("Arkane statmech receipt is invalid or not a completed native run")
    contract = validate_contract(reopen(receipt["contract"]))
    if receipt.get("route_id") != contract.get("route_id"):
        raise ValueError("Arkane statmech receipt reaction identity mismatch")
    for name in ("input", "provenance"):
        if receipt.get(name) != contract[name]:
            raise ValueError("Arkane rotor receipt/contract binding mismatch")
        reopen(receipt[name])
    output = reopen(receipt["output"])
    if str(output) != contract["expected_output"] or receipt.get("adapter_sources") != contract["adapter_sources"]:
        raise ValueError("Arkane rotor receipt output/adapter binding mismatch")
    native = receipt.get("native")
    if not isinstance(native, dict) or not isinstance(native.get("rmg_version"), str) or not native["rmg_version"] \
            or not isinstance(native.get("source_artifacts"), dict) or not native["source_artifacts"]:
        raise ValueError("Arkane native source/version evidence is incomplete")
    for source in native["source_artifacts"].values():
        reopen(source)
    _validate_output(output, contract["stationary_points"], receipt.get("stationary_points"))
    return receipt


def validate_rotor_receipt(workdir: Path) -> dict:
    path = Path(workdir).resolve() / RECEIPT_NAME
    try:
        evidence = {"policy_id": POLICY_ID, "status": "passed", "receipt": artifact(path)}
        receipt = validate_rotor_evidence(evidence)
        if Path(receipt["output"]["path"]) != Path(workdir).resolve() / "output.py":
            raise ValueError("Arkane rotor receipt output belongs to a different workspace")
        return evidence
    except (OSError, TypeError) as error:
        raise ValueError("Arkane statmech receipt is absent or invalid") from error


def validate_rotor_evidence(evidence: dict) -> dict:
    """Reopen all hashes and validate a manifest's portable admission receipt."""
    try:
        if not isinstance(evidence, dict) or set(evidence) != {"policy_id", "status", "receipt"} \
                or evidence.get("policy_id") != POLICY_ID or evidence.get("status") != "passed":
            raise ValueError("Arkane statmech admission evidence is absent or invalid")
        path = reopen(evidence["receipt"])
        receipt = _validate_receipt_document(read_document(path))
        if path.name != RECEIPT_NAME or Path(receipt["output"]["path"]).parent != path.parent:
            raise ValueError("Arkane statmech receipt belongs to a different workspace")
        return receipt
    except (OSError, KeyError, TypeError, SyntaxError, OverflowError) as error:
        raise ValueError("Arkane statmech evidence is incomplete or invalid") from error
