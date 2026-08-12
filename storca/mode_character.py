"""Internal-coordinate fingerprints for stable vibrational band identity."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Mapping

import numpy as np

from .ir_modes import NormalModeSet, _hungarian_minimize


BOHR_PER_ANGSTROM = 1.889726125
_COVALENT_RADII_ANGSTROM = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}


@dataclass(frozen=True)
class InternalCoordinate:
    kind: str
    atoms: tuple[int, ...]
    label: str
    key: str


def infer_bonds(
    labels: Iterable[str], coordinates_bohr: np.ndarray, *, atom_indices: Iterable[int] | None = None,
) -> list[tuple[int, int]]:
    """Infer a conservative covalent graph from element radii and geometry."""
    labels = tuple(labels)
    coordinates = np.asarray(coordinates_bohr, dtype=float) / BOHR_PER_ANGSTROM
    selected = sorted(set(range(len(labels))) if atom_indices is None else {int(i) for i in atom_indices})
    bonds: list[tuple[int, int]] = []
    for offset, left in enumerate(selected):
        for right in selected[offset + 1:]:
            radius = _COVALENT_RADII_ANGSTROM.get(labels[left], 0.77) + _COVALENT_RADII_ANGSTROM.get(labels[right], 0.77)
            distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
            if 0.35 < distance <= 1.25 * radius + 0.15:
                bonds.append((left, right))
    return bonds


def _coordinate_label(kind: str, atoms: tuple[int, ...], labels: tuple[str, ...], local: dict[int, int]) -> str:
    elements = "-".join(labels[index] for index in atoms)
    local_atoms = "-".join(str(local[index]) for index in atoms)
    if kind == "stretch" and set(labels[index] for index in atoms) == {"C", "O"}:
        # Geometry supplies a useful non-authoritative carbonyl distinction.
        name = "C-O"
    else:
        name = elements
    return f"{kind}:{name}:{local_atoms}"


def _stretch_label(
    atoms: tuple[int, int], mode_set: NormalModeSet, local: dict[int, int],
) -> str:
    left, right = atoms
    elements = tuple(sorted((mode_set.atom_labels[left], mode_set.atom_labels[right])))
    distance = float(np.linalg.norm(
        mode_set.coordinates_bohr[left] - mode_set.coordinates_bohr[right]
    ) / BOHR_PER_ANGSTROM)
    chemical = {
        ("H", "O"): "O-H stretch",
        ("H", "N"): "N-H stretch",
        ("C", "H"): "C-H stretch",
        ("C", "N"): "C-N stretch",
    }.get(elements)
    if elements == ("C", "O"):
        chemical = "C=O stretch" if distance <= 1.30 else "C-O stretch"
    if chemical is None:
        chemical = f"{mode_set.atom_labels[left]}-{mode_set.atom_labels[right]} stretch"
    return f"{chemical}:{local[left]}-{local[right]}"


def _cycle_basis(adjacency: dict[int, set[int]], maximum_size: int = 8) -> list[tuple[int, ...]]:
    cycles: set[tuple[int, ...]] = set()
    nodes = sorted(adjacency)
    for start in nodes:
        stack: list[tuple[int, tuple[int, ...]]] = [(start, (start,))]
        while stack:
            current, path = stack.pop()
            if len(path) > maximum_size:
                continue
            for neighbor in adjacency[current]:
                if neighbor == start and len(path) >= 3:
                    rotations = [path[i:] + path[:i] for i in range(len(path))]
                    reverse = tuple(reversed(path))
                    rotations += [reverse[i:] + reverse[:i] for i in range(len(reverse))]
                    cycles.add(min(rotations))
                elif neighbor > start and neighbor not in path:
                    stack.append((neighbor, path + (neighbor,)))
    # Retain chordless rings so fused systems do not explode into redundant cycles.
    chordless = []
    for cycle in sorted(cycles, key=lambda item: (len(item), item)):
        edges = {
            tuple(sorted((cycle[i], cycle[(i + 1) % len(cycle)]))) for i in range(len(cycle))
        }
        if any(
            tuple(sorted((cycle[i], cycle[j]))) not in edges
            and cycle[j] in adjacency[cycle[i]]
            for i in range(len(cycle)) for j in range(i + 1, len(cycle))
        ):
            continue
        chordless.append(cycle)
    return chordless


def build_internal_coordinates(
    mode_set: NormalModeSet, *, atom_indices: Iterable[int] | None = None,
    bonds: Iterable[tuple[int, int]] | None = None,
    atom_mapping: Mapping[int, int] | None = None,
) -> list[InternalCoordinate]:
    """Build stretches, bends, torsions, and ring deformations for target atoms."""
    selected = sorted(set(range(len(mode_set.atom_labels))) if atom_indices is None else {int(i) for i in atom_indices})
    if not selected or min(selected) < 0 or max(selected) >= len(mode_set.atom_labels):
        raise ValueError("Mode-character target atoms are outside the Hessian")
    local = {
        atom: int(atom_mapping[atom]) if atom_mapping and atom in atom_mapping else position
        for position, atom in enumerate(selected)
    }
    if len(set(local.values())) != len(local):
        raise ValueError("Mode-character atom mapping must be one-to-one")
    labels = mode_set.atom_labels
    resolved_bonds = sorted({tuple(sorted((int(a), int(b)))) for a, b in (
        infer_bonds(labels, mode_set.coordinates_bohr, atom_indices=selected) if bonds is None else bonds
    ) if a in local and b in local and a != b})
    adjacency = {atom: set() for atom in selected}
    coordinates: list[InternalCoordinate] = []
    for left, right in resolved_bonds:
        adjacency[left].add(right)
        adjacency[right].add(left)
        atoms = (left, right)
        coordinates.append(InternalCoordinate(
            "stretch", atoms, _stretch_label(atoms, mode_set, local),
            "stretch:" + "-".join(str(index) for index in sorted((local[left], local[right]))),
        ))
    for center in selected:
        neighbors = sorted(adjacency[center])
        for i, left in enumerate(neighbors):
            for right in neighbors[i + 1:]:
                atoms = (left, center, right)
                outer = sorted((local[left], local[right]))
                coordinates.append(InternalCoordinate(
                    "bend", atoms, _coordinate_label("bend", atoms, labels, local),
                    f"bend:{outer[0]}-{local[center]}-{outer[1]}",
                ))
    seen_torsions: set[tuple[int, ...]] = set()
    for center_left, center_right in resolved_bonds:
        for outer_left in adjacency[center_left] - {center_right}:
            for outer_right in adjacency[center_right] - {center_left}:
                atoms = (outer_left, center_left, center_right, outer_right)
                canonical = min(atoms, tuple(reversed(atoms)))
                canonical_local = min(
                    tuple(local[index] for index in atoms),
                    tuple(local[index] for index in reversed(atoms)),
                )
                if canonical in seen_torsions:
                    continue
                seen_torsions.add(canonical)
                coordinates.append(InternalCoordinate(
                    "torsion", atoms, _coordinate_label("torsion", atoms, labels, local),
                    "torsion:" + "-".join(str(index) for index in canonical_local),
                ))
    for cycle_index, cycle in enumerate(_cycle_basis(adjacency)):
        mapped_cycle = tuple(local[index] for index in cycle)
        reverse = tuple(reversed(mapped_cycle))
        canonical_cycle = min(
            *(mapped_cycle[i:] + mapped_cycle[:i] for i in range(len(mapped_cycle))),
            *(reverse[i:] + reverse[:i] for i in range(len(reverse))),
        )
        coordinates.append(InternalCoordinate(
            "ring_deformation", cycle, f"ring deformation:{len(cycle)}-membered",
            "ring_deformation:" + "-".join(str(index) for index in canonical_cycle),
        ))
    return coordinates


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    left, right = a - b, c - b
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator <= 1e-14:
        return 0.0
    return math.acos(float(np.clip(np.dot(left, right) / denominator, -1.0, 1.0)))


def _dihedral(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    axis = c - b
    norm = float(np.linalg.norm(axis))
    if norm <= 1e-14:
        return 0.0
    axis /= norm
    left = (a - b) - np.dot(a - b, axis) * axis
    right = (d - c) - np.dot(d - c, axis) * axis
    if np.linalg.norm(left) <= 1e-14 or np.linalg.norm(right) <= 1e-14:
        return 0.0
    return math.atan2(float(np.dot(axis, np.cross(left, right))), float(np.dot(left, right)))


def _coordinate_value(coordinate: InternalCoordinate, geometry: np.ndarray) -> float:
    atoms = coordinate.atoms
    if coordinate.kind == "stretch":
        return float(np.linalg.norm(geometry[atoms[0]] - geometry[atoms[1]]))
    if coordinate.kind == "bend":
        return _angle(geometry[atoms[0]], geometry[atoms[1]], geometry[atoms[2]])
    if coordinate.kind == "torsion":
        return _dihedral(*(geometry[index] for index in atoms))
    if coordinate.kind == "ring_deformation":
        points = geometry[list(atoms)]
        centered = points - points.mean(axis=0)
        _, _, right_t = np.linalg.svd(centered, full_matrices=False)
        normal = right_t[-1]
        return float(np.sqrt(np.mean(np.square(centered @ normal))))
    raise ValueError(f"Unsupported internal coordinate kind '{coordinate.kind}'")


def _periodic_delta(left: float, right: float) -> float:
    return (left - right + math.pi) % (2.0 * math.pi) - math.pi


def mode_character_fingerprints(
    mode_set: NormalModeSet, *, atom_indices: Iterable[int] | None = None,
    bonds: Iterable[tuple[int, int]] | None = None, maximum_components: int = 8,
    atom_mapping: Mapping[int, int] | None = None,
) -> list[dict]:
    """Project every Cartesian normal mode onto a target internal-coordinate basis."""
    coordinates = build_internal_coordinates(
        mode_set, atom_indices=atom_indices, bonds=bonds, atom_mapping=atom_mapping,
    )
    target = sorted(set(range(len(mode_set.atom_labels))) if atom_indices is None else {int(i) for i in atom_indices})
    fingerprints: list[dict] = []
    epsilon = 1e-4
    for mode_index, vector in enumerate(mode_set.vectors):
        target_norm = float(np.linalg.norm(vector[target]))
        full_norm = float(np.linalg.norm(vector))
        direction = vector / full_norm if full_norm > 0 else vector
        plus = mode_set.coordinates_bohr + epsilon * direction
        minus = mode_set.coordinates_bohr - epsilon * direction
        raw = []
        for coordinate in coordinates:
            plus_value = _coordinate_value(coordinate, plus)
            minus_value = _coordinate_value(coordinate, minus)
            delta = (
                _periodic_delta(plus_value, minus_value)
                if coordinate.kind == "torsion" else plus_value - minus_value
            )
            raw.append((coordinate, (delta / (2.0 * epsilon)) ** 2))
        total = sum(value for _, value in raw)
        components = [
            {"key": coordinate.key, "kind": coordinate.kind, "label": coordinate.label,
             "atoms": list(coordinate.atoms), "participation": float(value / total)}
            for coordinate, value in sorted(raw, key=lambda item: item[1], reverse=True)
            if total > 0 and value / total >= 1e-6
        ][:maximum_components]
        kinds: dict[str, float] = {}
        for component in components:
            kinds[component["kind"]] = kinds.get(component["kind"], 0.0) + component["participation"]
        fingerprints.append({
            "schema_version": 1,
            "mode": mode_index,
            "frequency_cm-1": float(mode_set.frequencies_cm_1[mode_index]),
            "target_localization": 0.0 if full_norm <= 0 else min(1.0, target_norm ** 2 / full_norm ** 2),
            "components": components,
            "coordinate_kind_participation": kinds,
            "dominant_coordinate": components[0]["label"] if components else None,
            "character_status": "resolved" if components else "unresolved",
        })
    return fingerprints


def fingerprint_similarity(left: dict, right: dict) -> float:
    left_values = {item["key"]: float(item["participation"]) for item in left.get("components", [])}
    right_values = {item["key"]: float(item["participation"]) for item in right.get("components", [])}
    keys = sorted(set(left_values) | set(right_values))
    if not keys:
        return 0.0
    a = np.asarray([left_values.get(key, 0.0) for key in keys])
    b = np.asarray([right_values.get(key, 0.0) for key in keys])
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return 0.0 if denominator <= 0 else float(np.dot(a, b) / denominator)


def match_mode_fingerprints(reference: list[dict], candidate: list[dict]) -> list[dict]:
    """Match equal target-coordinate spaces using chemical character first."""
    if len(reference) > len(candidate):
        raise ValueError("Candidate mode set is smaller than the reference mode set")
    similarity = np.asarray([
        [fingerprint_similarity(left, right) for right in candidate] for left in reference
    ], dtype=float)
    gaps = np.asarray([
        [abs(float(left["frequency_cm-1"]) - float(right["frequency_cm-1"])) for right in candidate]
        for left in reference
    ])
    localization = np.asarray([
        [min(float(left.get("target_localization", 0.0)), float(right.get("target_localization", 0.0)))
         for right in candidate] for left in reference
    ])
    combined = 0.75 * similarity + 0.15 * np.exp(-np.square(gaps / 300.0)) + 0.10 * localization
    assignment = _hungarian_minimize(1.0 - combined)
    matches = []
    for reference_index, candidate_index in enumerate(assignment):
        score = float(similarity[reference_index, candidate_index])
        confidence = "high" if score >= 0.80 else "medium" if score >= 0.50 else "low"
        matches.append({
            "reference_mode": int(reference[reference_index]["mode"]),
            "candidate_mode": int(candidate[candidate_index]["mode"]),
            "mode_character_similarity": score,
            "displacement_overlap": score,
            "frequency_gap_cm-1": float(gaps[reference_index, candidate_index]),
            "confidence": confidence,
            "matching_basis": "internal_coordinate_fingerprint",
            "mode_character": candidate[candidate_index],
        })
    return matches
