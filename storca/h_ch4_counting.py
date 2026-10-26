"""Small, opt-in geometric admission of the declared H + CH4 counting model.

This does not find a saddle, verify an IRC, modify coordinates, diagnose the
electronic wavefunction or validate RRHO. Numerical symmetry tolerances are
engineering declarations, not bounds on a rate prediction's physical error.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import re

import numpy as np


POLICY_ID = "h_ch4_native_geometric_counting_v1"
LABELS = ("r1", "r2", "p1", "p2", "TS")
RMSD_TOLERANCE_ANGSTROM = 1e-3
MAX_DISPLACEMENT_TOLERANCE_ANGSTROM = 2e-3
_MAPS = {"r1": {1: "H"}, "r2": {2: "C", 3: "H", 4: "H", 5: "H", 6: "H"},
         "p1": {1: "H", 6: "H"}, "p2": {2: "C", 3: "H", 4: "H", 5: "H"}}
_EDGES = {"r1": set(), "r2": {(2, n) for n in (3, 4, 5, 6)},
          "p1": {(1, 6)}, "p2": {(2, n) for n in (3, 4, 5)}}


def expected_counts():
    """Fresh declarations: site equivalence is already in rotor symmetries."""
    return {label: {"external_symmetry": count, "optical_isomers": 1}
            for label, count in zip(LABELS, (1, 12, 2, 6, 3))}


def _binding(path):
    if not isinstance(path, (str, Path)) or isinstance(path, bool):
        raise ValueError("native XYZ path is missing or malformed")
    target = Path(path).resolve(strict=True)
    if not target.is_file():
        raise ValueError("native XYZ source is not a file")
    return {"path": str(target), "sha256": hashlib.sha256(target.read_bytes()).hexdigest()}


def _reopen(binding):
    if (not isinstance(binding, dict) or set(binding) != {"path", "sha256"}
            or not isinstance(binding["path"], str) or not Path(binding["path"]).is_absolute()
            or not isinstance(binding["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", binding["sha256"]) is None
            or _binding(binding["path"]) != binding):
        raise ValueError("native XYZ source hash changed or binding is malformed")
    return Path(binding["path"])


def _xyz(binding):
    lines = _reopen(binding).read_text().splitlines()
    if len(lines) < 3 or re.fullmatch(r"[1-6]", lines[0].strip()) is None:
        raise ValueError("native XYZ requires one to six explicit nuclei")
    count = int(lines[0])
    rows = [line.split() for line in lines[2:]]
    while rows and not rows[-1]:
        rows.pop()
    if len(rows) != count or any(len(row) != 4 for row in rows):
        raise ValueError("native XYZ is incomplete or contains extra records")
    elements = [row[0] for row in rows]
    xyz = np.asarray([[float(item) for item in row[1:]] for row in rows], dtype=float)
    if not np.isfinite(xyz).all():
        raise ValueError("native XYZ contains nonfinite coordinates")
    if any(np.linalg.norm(a - b) <= 1e-8 for a, b in itertools.combinations(xyz, 2)):
        raise ValueError("native XYZ contains coincident nuclei")
    return elements, xyz


def _protium_identities(elements, isotopes):
    if not isinstance(isotopes, (list, tuple)) or len(isotopes) != len(elements):
        raise ValueError("per-nucleus isotope declaration is malformed")
    identities = []
    for element, isotope in zip(elements, isotopes):
        if type(isotope) is not int:
            raise ValueError("isotope numbers must be integers, never booleans")
        if element == "H" and isotope in (0, 1):
            identities.append(("H", 1))
        elif element == "C" and isotope in (0, 12):
            identities.append(("C", 12))
        else:
            raise ValueError("this H/CH4 counting policy supports protium and carbon-12 only; no D/T")
    return identities


def _species_identity(label, smiles, elements):
    from rdkit import Chem
    if not isinstance(smiles, str) or not smiles:
        raise ValueError("complete mapped SMILES is required for every species")
    parameters = Chem.SmilesParserParams()
    parameters.removeHs = False
    mol = Chem.MolFromSmiles(smiles, parameters)
    if mol is None:
        raise ValueError("mapped species SMILES is malformed")
    mol = Chem.AddHs(mol)
    atoms = list(mol.GetAtoms())
    maps = [atom.GetAtomMapNum() for atom in atoms]
    if not maps or any(n <= 0 for n in maps) or len(set(maps)) != len(maps):
        raise ValueError("complete unique nuclear maps are required, including every hydrogen")
    if {atom.GetAtomMapNum(): atom.GetSymbol() for atom in atoms} != _MAPS[label]:
        raise ValueError("species nuclear maps differ from the declared H/CH4 reaction")
    if [atom.GetSymbol() for atom in atoms] != elements:
        raise ValueError("native XYZ nuclear order differs from mapped SMILES order")
    edges = {tuple(sorted((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())))
             for bond in mol.GetBonds()}
    if (edges != _EDGES[label] or any(bond.GetBondTypeAsDouble() != 1 for bond in mol.GetBonds())
            or any(atom.GetFormalCharge() for atom in atoms)
            or sum(atom.GetNumRadicalElectrons() for atom in atoms) != (label in {"r1", "p2"})):
        raise ValueError("mapped species graph/state differs from the declared neutral H/CH4 reaction")
    isotopes = [int(atom.GetIsotope()) for atom in atoms]
    return isotopes, maps


def _permutations(identities):
    groups = {}
    for index, identity in enumerate(identities):
        groups.setdefault(identity, []).append(index)
    count = 1
    choices = []
    for indices in groups.values():
        options = list(itertools.permutations(indices))
        count *= len(options)
        if count > 120:
            raise ValueError("H/CH4 nuclear-permutation bound exceeded")
        choices.append((indices, options))
    for selection in itertools.product(*(options for _, options in choices)):
        permutation = list(range(len(identities)))
        for (indices, _), targets in zip(choices, selection):
            for original, target in zip(indices, targets):
                permutation[original] = target
        yield permutation


def _proper_fit(source, target):
    x, y = source - source.mean(axis=0), target - target.mean(axis=0)
    left, _, right = np.linalg.svd(x.T @ y)
    parity = 1. if np.linalg.det(left @ right) >= 0 else -1.
    rotation = left @ np.diag((1., 1., parity)) @ right
    distances = np.linalg.norm(x @ rotation - y, axis=1)
    return {"rmsd_angstrom": float(np.sqrt(np.mean(distances**2))),
            "max_displacement_angstrom": float(np.max(distances))}


def _admitted(fit):
    return (fit["rmsd_angstrom"] <= RMSD_TOLERANCE_ANGSTROM
            and fit["max_displacement_angstrom"] <= MAX_DISPLACEMENT_TOLERANCE_ANGSTROM)


def _geometry(xyz, identities):
    proper, mirror = [], []
    reflected = (xyz - xyz.mean(axis=0)) @ np.diag((-1., 1., 1.))
    examined = 0
    for permutation in _permutations(identities):
        examined += 1
        fit = _proper_fit(xyz, xyz[permutation])
        if _admitted(fit):
            proper.append({"permutation": permutation, **fit})
        reflected_fit = _proper_fit(reflected, xyz[permutation])
        mirror.append({"permutation": permutation, **reflected_fit})
    operations = {tuple(item["permutation"]) for item in proper}
    closed = all(tuple(left[right[n]] for n in range(len(xyz))) in operations
                 for left in operations for right in operations)
    best = min(mirror, key=lambda fit: (max(fit["rmsd_angstrom"] / RMSD_TOLERANCE_ANGSTROM,
                                          fit["max_displacement_angstrom"] / MAX_DISPLACEMENT_TOLERANCE_ANGSTROM),
                                      fit["permutation"]))
    return {"permutations_examined": examined, "observed_external_symmetry": len(proper),
            "proper_symmetry_operations": proper, "accepted_operations_closed": closed,
            "achiral_within_declared_tolerances": _admitted(best),
            "best_reflected_superposition": best}


def _resolved_shape(label, xyz):
    """Reject collapsed/linear imitations before counting approximate fits."""
    values = np.linalg.svd(xyz - xyz.mean(axis=0), compute_uv=False)
    required_rank = {"r1": 0, "r2": 3, "p1": 1, "p2": 2, "TS": 3}[label]
    resolved_rank = int(sum(value > MAX_DISPLACEMENT_TOLERANCE_ANGSTROM for value in values))
    if resolved_rank < required_rank:
        raise ValueError("native nuclear geometry is collapsed, linear or rank-unresolved for this counting model")
    return {"centered_coordinate_singular_values_angstrom": [float(value) for value in values],
            "minimum_required_coordinate_rank": required_rank,
            "rank_resolved_above_max_displacement_tolerance": resolved_rank}


def assess_h_ch4_counting(points, counting):
    """Return a source-bound pass/refusal receipt; never launch or alter work."""
    reasons, observed, inputs = [], {}, {}
    try:
        declared = json.loads(json.dumps(counting, allow_nan=False))
    except (TypeError, ValueError):
        declared = None
    expected = expected_counts()
    if (not isinstance(counting, dict) or set(counting) != set(LABELS)
            or any(not isinstance(row, dict) or set(row) != {"external_symmetry", "optical_isomers"}
                   or any(type(n) is not int for n in row.values())
                   for row in counting.values()) or counting != expected):
        reasons.append("declared_counting_differs_from_frozen_H_CH4_model")
    if not isinstance(points, dict) or set(points) != set(LABELS):
        reasons.append("exactly_five_H_CH4_stationary_points_required")
    else:
        for label in LABELS:
            try:
                point = points[label]
                if not isinstance(point, dict):
                    raise ValueError("stationary-point record is malformed")
                source = _binding(point.get("optimized_xyz"))
                elements, xyz = _xyz(source)
                if label == "TS":
                    if point.get("smiles") is not None or Counter(elements) != Counter({"C": 1, "H": 5}):
                        raise ValueError("TS requires the fixed C1H5 nucleus inventory, not a species SMILES")
                    isotopes, maps = [0] * len(elements), None
                else:
                    isotopes, maps = _species_identity(label, point.get("smiles"), elements)
                identities = _protium_identities(elements, isotopes)
                if "isotopes" in point:
                    if _protium_identities(elements, point["isotopes"]) != identities:
                        raise ValueError("native/mapped isotope declarations disagree")
                if "isotope" in point:
                    if len(elements) != 1 or _protium_identities(elements, [point["isotope"]]) != identities:
                        raise ValueError("atomic isotope declaration disagrees")
                inputs[label] = {"optimized_xyz": source, "smiles": point.get("smiles"),
                                 **({"isotopes": list(point["isotopes"])} if "isotopes" in point else {}),
                                 **({"isotope": point["isotope"]} if "isotope" in point else {})}
                shape = _resolved_shape(label, xyz)
                result = _geometry(xyz, identities)
                observed[label] = {"native_xyz": source, "elements": elements, "mapped_atom_numbers": maps,
                                   "isotope_identities": [list(identity) for identity in identities], **shape, **result}
                if result["observed_external_symmetry"] != expected[label]["external_symmetry"]:
                    reasons.append(label + ":native_geometric_symmetry_differs_from_declaration")
                if not result["accepted_operations_closed"]:
                    reasons.append(label + ":symmetry_ambiguous_at_declared_numerical_tolerances")
                if not result["achiral_within_declared_tolerances"]:
                    reasons.append(label + ":optical_one_not_supported_by_reflected_superposition")
                _reopen(source)
            except (ValueError, TypeError, KeyError, OSError, np.linalg.LinAlgError) as error:
                reasons.append(label + ":" + str(error))
    return {"schema_version": 1, "kind": POLICY_ID, "status": "failed" if reasons else "passed",
            "expected_counting": expected, "declared_counting": declared,
            "inputs": inputs, "points": observed, "unresolved_reasons": reasons,
            "policy": {"rmsd_tolerance_angstrom": RMSD_TOLERANCE_ANGSTROM,
                       "max_displacement_tolerance_angstrom": MAX_DISPLACEMENT_TOLERANCE_ANGSTROM,
                       "tolerance_scope": "Engineering numerical symmetry, not physical rate accuracy.",
                       "maximum_permutations_per_point": 120,
                       "isotope_scope": "Protium and carbon-12 only; unspecified isotopes mean these declared defaults.",
                       "reaction_degeneracy": 1, "additional_site_multiplier": 1},
            "scientific_basis": {"url": "https://d-nb.info/1239993528/34", "equation": "6",
                "interpretation": "For H + CH4 through C3v TS, rotor symmetry ratio 1*12/3 already supplies four equivalent sites; no additional factor four."},
            "limitations": ["Geometric symmetry does not verify an electronic wavefunction, saddle, reaction mode or IRC.",
                            "No coordinates are symmetrized; no harmonic, tunneling or kinetic-model accuracy is certified.",
                            "Mapped nuclear identities are not artificial physical distinguishability of identical isotopes."],
            "native_jobs_launched": 0, "physical_rate_accuracy_claim_allowed": False,
            "route_verification_claim_allowed": False, "network_reinjection_allowed": False}


def reopen_h_ch4_counting(receipt):
    """Reopen bound XYZ bytes and reproduce a passed receipt without mutation."""
    if isinstance(receipt, (str, Path)):
        receipt = json.loads(Path(receipt).read_text())
    if (not isinstance(receipt, dict) or receipt.get("kind") != POLICY_ID
            or receipt.get("status") != "passed" or type(receipt.get("schema_version")) is not int
            or receipt["schema_version"] != 1 or set(receipt.get("inputs", {})) != set(LABELS)):
        raise ValueError("A passed H/CH4 native counting receipt is required")
    points = {}
    for label, record in receipt["inputs"].items():
        if not isinstance(record, dict) or not {"optimized_xyz", "smiles"} <= set(record) \
                or set(record) - {"optimized_xyz", "smiles", "isotopes", "isotope"}:
            raise ValueError("Malformed source-bound counting inputs")
        points[label] = {**record, "optimized_xyz": _reopen(record["optimized_xyz"])}
    repeated = assess_h_ch4_counting(points, receipt.get("declared_counting"))
    if repeated != receipt:
        raise ValueError("Native counting receipt does not reproduce its geometric/source claims")
    return repeated
