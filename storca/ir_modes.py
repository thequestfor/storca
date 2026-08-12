"""ORCA Hessian normal-mode parsing and geometry-aligned mode matching."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class NormalModeSet:
    frequencies_cm_1: np.ndarray
    vectors: np.ndarray
    atom_labels: tuple[str, ...]
    coordinates_bohr: np.ndarray


def _section(lines: list[str], name: str) -> list[str]:
    marker = f"${name}"
    try:
        start = next(index for index, line in enumerate(lines) if line.strip() == marker) + 1
    except StopIteration as error:
        raise ValueError(f"ORCA Hessian has no {marker} section") from error
    end = next((index for index in range(start, len(lines)) if lines[index].lstrip().startswith("$")), len(lines))
    return lines[start:end]


def parse_orca_normal_modes(path: Path) -> NormalModeSet:
    """Parse frequencies, Cartesian modes, and atoms from an ORCA ``.hess`` file."""
    lines = Path(path).read_text(errors="replace").splitlines()

    frequency_lines = [line.split() for line in _section(lines, "vibrational_frequencies") if line.strip()]
    if not frequency_lines:
        raise ValueError("Empty ORCA vibrational-frequency section")
    count = int(frequency_lines[0][0])
    frequencies = np.full(count, np.nan, dtype=float)
    for fields in frequency_lines[1:]:
        if len(fields) >= 2:
            frequencies[int(fields[0])] = float(fields[1].replace("D", "E").replace("d", "e"))
    if not np.all(np.isfinite(frequencies)):
        raise ValueError("Incomplete ORCA vibrational-frequency section")

    mode_lines = _section(lines, "normal_modes")
    position = next((index for index, line in enumerate(mode_lines) if line.strip()), None)
    if position is None:
        raise ValueError("Empty ORCA normal-mode section")
    dimensions = mode_lines[position].split()
    rows, columns = int(dimensions[0]), int(dimensions[1])
    matrix = np.full((rows, columns), np.nan, dtype=float)
    position += 1
    while position < len(mode_lines):
        if np.all(np.isfinite(matrix)):
            break
        if not mode_lines[position].strip():
            position += 1
            continue
        if mode_lines[position].lstrip().startswith("#"):
            break
        column_ids = [int(value) for value in mode_lines[position].split()]
        position += 1
        populated_rows = 0
        while position < len(mode_lines) and populated_rows < rows:
            fields = mode_lines[position].split()
            position += 1
            if not fields:
                continue
            if len(fields) != len(column_ids) + 1:
                raise ValueError("Malformed ORCA normal-mode matrix row")
            row = int(fields[0])
            matrix[row, column_ids] = [float(value.replace("D", "E").replace("d", "e")) for value in fields[1:]]
            populated_rows += 1
    if rows != columns or columns != count or not np.all(np.isfinite(matrix)) or rows % 3:
        raise ValueError("Inconsistent ORCA normal-mode matrix dimensions")

    atom_lines = [line.split() for line in _section(lines, "atoms") if line.strip()]
    atom_count = int(atom_lines[0][0])
    if len(atom_lines) - 1 < atom_count or atom_count * 3 != rows:
        raise ValueError("ORCA Hessian atom count does not match normal modes")
    labels = tuple(fields[0] for fields in atom_lines[1:atom_count + 1])
    coordinates = np.asarray([[float(value) for value in fields[2:5]] for fields in atom_lines[1:atom_count + 1]])
    # ORCA stores Cartesian components in rows and modes in columns.
    vectors = matrix.T.reshape(columns, atom_count, 3)
    return NormalModeSet(frequencies, vectors, labels, coordinates)


def _alignment_rotation(moving: np.ndarray, reference: np.ndarray) -> np.ndarray:
    moving_centered = moving - moving.mean(axis=0)
    reference_centered = reference - reference.mean(axis=0)
    left, _, right_t = np.linalg.svd(moving_centered.T @ reference_centered)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    return rotation


def _hungarian_minimize(cost: np.ndarray) -> np.ndarray:
    """Return one minimum-cost column per row using a rectangular Hungarian method."""
    cost = np.asarray(cost, dtype=float)
    rows, columns = cost.shape
    if rows > columns:
        raise ValueError("Hungarian implementation requires rows <= columns")
    u = np.zeros(rows + 1)
    v = np.zeros(columns + 1)
    p = np.zeros(columns + 1, dtype=int)
    way = np.zeros(columns + 1, dtype=int)
    for row in range(1, rows + 1):
        p[0] = row
        column0 = 0
        minimum = np.full(columns + 1, np.inf)
        used = np.zeros(columns + 1, dtype=bool)
        while True:
            used[column0] = True
            row0 = p[column0]
            delta = np.inf
            column1 = 0
            for column in range(1, columns + 1):
                if used[column]:
                    continue
                current = cost[row0 - 1, column - 1] - u[row0] - v[column]
                if current < minimum[column]:
                    minimum[column] = current
                    way[column] = column0
                if minimum[column] < delta:
                    delta = minimum[column]
                    column1 = column
            for column in range(columns + 1):
                if used[column]:
                    u[p[column]] += delta
                    v[column] -= delta
                else:
                    minimum[column] -= delta
            column0 = column1
            if p[column0] == 0:
                break
        while True:
            column1 = way[column0]
            p[column0] = p[column1]
            column0 = column1
            if column0 == 0:
                break
    assignment = np.full(rows, -1, dtype=int)
    for column in range(1, columns + 1):
        if p[column]:
            assignment[p[column] - 1] = column - 1
    return assignment


def match_normal_modes(reference: NormalModeSet, candidate: NormalModeSet) -> list[dict]:
    """Match candidate modes to a reference using aligned displacement overlap."""
    if reference.atom_labels != candidate.atom_labels:
        raise ValueError("Normal-mode matching requires identical atom ordering")
    if reference.vectors.shape != candidate.vectors.shape:
        raise ValueError("Normal-mode matching requires equal mode dimensions")
    rotation = _alignment_rotation(candidate.coordinates_bohr, reference.coordinates_bohr)
    reference_vectors = reference.vectors.reshape(len(reference.vectors), -1)
    candidate_vectors = (candidate.vectors @ rotation).reshape(len(candidate.vectors), -1)
    reference_norms = np.linalg.norm(reference_vectors, axis=1)
    candidate_norms = np.linalg.norm(candidate_vectors, axis=1)
    denominator = np.outer(reference_norms, candidate_norms)
    overlap = np.divide(
        np.abs(reference_vectors @ candidate_vectors.T), denominator,
        out=np.zeros_like(denominator), where=denominator > 0,
    )
    gaps = np.abs(reference.frequencies_cm_1[:, None] - candidate.frequencies_cm_1[None, :])
    # Frequency is a soft guard against assigning unrelated modes that happen
    # to have similar Cartesian character.  Displacement overlap remains the
    # dominant term for crossings and reordered modes.
    frequency_score = np.exp(-np.square(gaps / 250.0))
    combined = 0.85 * overlap + 0.15 * frequency_score
    assignment = _hungarian_minimize(1.0 - combined)
    matches = []
    for reference_index, candidate_index in enumerate(assignment):
        value = float(overlap[reference_index, candidate_index])
        gap = float(gaps[reference_index, candidate_index])
        confidence = "high" if value >= 0.75 else "medium" if value >= 0.45 else "low"
        matches.append({
            "reference_mode": reference_index,
            "candidate_mode": int(candidate_index),
            "displacement_overlap": value,
            "frequency_gap_cm-1": gap,
            "confidence": confidence,
        })
    return matches


def mode_localization_fraction(mode_set: NormalModeSet, mode: int, atom_count: int) -> float:
    """Return the Cartesian displacement fraction on the first ``atom_count`` atoms."""
    if not 0 < atom_count <= len(mode_set.atom_labels) or not 0 <= mode < len(mode_set.vectors):
        raise ValueError("Invalid localization atom or mode range")
    squared = np.square(mode_set.vectors[mode]).sum(axis=1)
    total = float(squared.sum())
    return 0.0 if total == 0 else float(squared[:atom_count].sum() / total)


def local_stretch_mode_assignments(
    mode_set: NormalModeSet, bonds: list[dict], *, minimum_frequency_cm_1: float = 1000.0,
) -> list[dict]:
    """Assign one normal mode to each local X-H coordinate across any cluster size."""
    assignments = []
    used_modes: set[int] = set()
    mode_norms = np.linalg.norm(mode_set.vectors.reshape(len(mode_set.vectors), -1), axis=1)
    for bond in bonds:
        heavy = int(bond["heavy_atom"])
        hydrogen = int(bond["hydrogen_atom"])
        if min(heavy, hydrogen) < 0 or max(heavy, hydrogen) >= len(mode_set.atom_labels):
            continue
        axis = mode_set.coordinates_bohr[hydrogen] - mode_set.coordinates_bohr[heavy]
        norm = float(np.linalg.norm(axis))
        if norm <= 0:
            continue
        axis /= norm
        relative = mode_set.vectors[:, hydrogen, :] - mode_set.vectors[:, heavy, :]
        projections = np.abs(relative @ axis)
        scores = np.divide(projections, mode_norms, out=np.zeros_like(projections), where=mode_norms > 0)
        eligible = [
            index for index, frequency in enumerate(mode_set.frequencies_cm_1)
            if frequency >= minimum_frequency_cm_1 and index not in used_modes
        ]
        if not eligible:
            continue
        selected = max(eligible, key=lambda index: (float(scores[index]), float(mode_set.frequencies_cm_1[index])))
        used_modes.add(selected)
        value = min(1.0, float(scores[selected]))
        assignments.append({
            **bond,
            "mode": selected,
            "local_stretch_projection": value,
            "confidence": "high" if value >= 0.45 else "medium" if value >= 0.25 else "low",
        })
    return assignments
