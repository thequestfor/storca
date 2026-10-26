"""Recover an evaluated ORCA optimization frame, not a relabelled EnGrad header.

This pure, mapped seven-atom OH/CH4 reader joins property-file geometry, total
energy and gradient by GeometryIndex, then checks the native output witnesses.
It does not read files, dispatch jobs, establish process ownership, check a
saddle index, or admit a rate. Callers bind all source bytes independently.
"""
from __future__ import annotations

import json
import math
import re

import numpy as np


LABELS = ("O", "H", "C", "H", "H", "H", "H")
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
_KNOWN = {"Geometry", "Single_Point_Data", "SCF_Nuc_Gradient"}
_ENERGY_TOL = 1e-12
_GRADIENT_PRINT_TOL = 5.1e-10
_FRAME_PRINT_TOL_A = 6e-7
_FINAL_FILE_TOL_A = 1e-9


def _finite(token):
    if not isinstance(token, str) or re.fullmatch(_NUMBER, token) is None:
        raise ValueError("Finite native numeric field required")
    value = float(token.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("Nonfinite native numeric field")
    return value


def _one(pattern, text, label):
    matches = list(re.finditer(pattern, text, re.M))
    if len(matches) != 1:
        raise ValueError("Unique " + label + " required")
    return matches[0]


def _field_count(body, field):
    if len(re.findall(r"(?m)^\s*&" + re.escape(field) + r"\b", body)) != 1:
        raise ValueError("Unique property field " + field + " required")


def _integer(body, field, expected=None):
    _field_count(body, field)
    match = _one(r"^\s*&" + re.escape(field)
                 + r'\s*\[\s*&Type\s*"Integer"\s*\]\s+(\d+)\s*$', body, field)
    value = int(match[1])
    if expected is not None and value != expected:
        raise ValueError("Unexpected property " + field)
    return value


def _string(body, field, expected):
    _field_count(body, field)
    match = _one(r"^\s*&" + re.escape(field)
                 + r'\s*\[\s*&Type\s*"String"\s*\]\s+"([^"\n]*)"\s*$', body, field)
    if match[1] != expected:
        raise ValueError("Unexpected property " + field)


def _sections(text):
    sections = {}
    matches = list(re.finditer(r"^\$(\w+)[ \t]*\n(.*?)^\$End[ \t]*$", text, re.M | re.S))
    headers = re.findall(r"(?m)^\$(Geometry|Single_Point_Data|SCF_Nuc_Gradient)[ \t]*$", text)
    relevant = [match for match in matches if match[1] in _KNOWN]
    if len(relevant) != len(headers):
        raise ValueError("Truncated or nested required property section")
    for match in relevant:
        body = match[2]
        if re.search(r"(?m)^\$", body):
            raise ValueError("Nested required property section")
        _field_count(body, "GeometryIndex")
        index_match = _one(r"^\s*&GeometryIndex\s+(\d+)\s*$", body, "GeometryIndex")
        index = int(index_match[1])
        key = (match[1], index)
        if key in sections:
            raise ValueError("Duplicate property section/GeometryIndex: " + repr(key))
        sections[key] = body
    if not sections:
        raise ValueError("Required evaluated-gradient property sections missing")
    return sections


def _geometry(body):
    _integer(body, "NAtoms", 7)
    _field_count(body, "CartesianCoordinates")
    header = _one(
        r'^\s*&CartesianCoordinates\s*\[\s*&Type\s*"Coordinates"\s*,\s*'
        r'&Dim\s*\(\s*7\s*,\s*4\s*\)\s*,\s*&Units\s*"Bohr"\s*\][ \t]*$',
        body, "mapped Bohr coordinate array")
    rows = [line.split() for line in body[header.end():].splitlines() if line.strip()]
    if len(rows) != 7 or any(len(row) != 4 for row in rows) or tuple(row[0] for row in rows) != LABELS:
        raise ValueError("Seven mapped O/H/C/H/H/H/H property coordinates required")
    return np.asarray([[_finite(value) for value in row[1:]] for row in rows])


def _energy(body):
    _field_count(body, "FinalEnergy")
    match = _one(r'^\s*&FinalEnergy\s*\[\s*&Type\s*"Double"\s*\]\s+('
                 + _NUMBER + r')(?:[ \t]+"[^"\n]*")?[ \t]*$', body, "total FinalEnergy")
    _field_count(body, "Converged")
    converged = _one(r'^\s*&Converged\s*\[\s*&Type\s*"Boolean"\s*\]\s+(true|false)[ \t]*$',
                     body, "SCF convergence flag")
    if converged[1] != "true":
        raise ValueError("Selected/final total energy is not converged")
    return _finite(match[1])


def _gradient(body):
    _integer(body, "NAtoms", 7)
    _integer(body, "Mult", 2)
    _integer(body, "State", 0)
    _integer(body, "Irrep", 0)
    _string(body, "Method", "SCF")
    _string(body, "Level", "Relaxed density")
    _field_count(body, "grad")
    header = _one(r'^\s*&grad\s*\[\s*&Type\s*"ArrayOfDoubles"\s*,\s*'
                  r'&Dim\s*\(\s*21\s*,\s*1\s*\)\s*\][ \t]*$', body, "21-component gradient")
    tail = body[header.end():]
    end = re.search(r"(?m)^\s*&", tail)
    data = tail[:end.start()] if end else tail
    rows = [line.split() for line in data.splitlines() if line.strip()]
    if (len(rows) != 22 or rows[0] != ["0"] or any(len(row) != 2 for row in rows[1:])
            or [row[0] for row in rows[1:]] != [str(index) for index in range(21)]):
        raise ValueError("Complete uniquely indexed 21-component gradient required")
    gradient = np.asarray([_finite(row[1]) for row in rows[1:]])
    return gradient


def _native_frames(text):
    frames = []
    markers = list(re.finditer(r"(?m)^\s*CARTESIAN COORDINATES \(ANGSTROEM\)[ \t]*$", text))
    for marker in markers:
        # ORCA prints these fixed seven-row tables in well under 4 KiB. Do not
        # repeatedly copy/split a many-megabyte MPI output tail for each frame.
        remainder = text[marker.end():marker.end() + 4096]
        lines = remainder.splitlines()
        position = 0
        while position < len(lines) and (not lines[position].strip() or set(lines[position].strip()) == {"-"}):
            position += 1
        rows = [line.split() for line in lines[position:position + 7]]
        if (len(rows) != 7 or any(len(row) != 4 for row in rows)
                or tuple(row[0] for row in rows) != LABELS):
            raise ValueError("Mapped native Cartesian frame missing or malformed")
        following = lines[position + 7].split() if len(lines) > position + 7 else []
        if len(following) == 4 and re.fullmatch(r"[A-Z][a-z]?", following[0]):
            raise ValueError("Native Cartesian frame has more than seven atoms")
        frame = np.asarray([[_finite(token) for token in row[1:]] for row in rows])
        frames.append((marker.start(), frame))
    if not frames:
        raise ValueError("Native mapped Cartesian frames missing")
    return frames


def _native_last_gradient(text):
    markers = list(re.finditer(r"(?m)^\s*CARTESIAN GRADIENT[ \t]*$", text))
    if not markers:
        raise ValueError("Native full Cartesian gradient missing")
    marker = markers[-1]
    tail = text[marker.end():marker.end() + 4096]
    lines = tail.splitlines()
    position = 0
    while position < len(lines) and (not lines[position].strip() or set(lines[position].strip()) == {"-"}):
        position += 1
    rows = [line.split() for line in lines[position:position + 7]]
    if (len(rows) != 7 or any(len(row) != 6 or row[2] != ":" for row in rows)
            or [row[0] for row in rows] != [str(index) for index in range(1, 8)]
            or tuple(row[1] for row in rows) != LABELS):
        raise ValueError("Complete mapped native 7-by-3 gradient required")
    following = lines[position + 7].split() if len(lines) > position + 7 else []
    if len(following) == 6 and following[0].isdigit() and following[2] == ":":
        raise ValueError("Native gradient has more than seven atoms")
    return marker.start(), np.asarray([_finite(value) for row in rows for value in row[3:]])


def evaluated_gradient_receipt(property_text, native_text, final_xyz_A, bohr_A):
    """Join exact evaluated E/g/geometry by index and check native witnesses.

    Select the last gradient-bearing property geometry at/before the final
    geometry. A final-energy pass may update EnGrad coordinates without a fresh
    gradient; its header is deliberately not an input to this API.
    """
    if not isinstance(property_text, str) or not isinstance(native_text, str):
        raise ValueError("Raw property/native output text required")
    if type(bohr_A) not in (float, int) or not math.isfinite(bohr_A) or bohr_A <= 0:
        raise ValueError("Explicit finite positive native Bohr-to-angstrom conversion required")
    final_supplied = np.asarray(final_xyz_A, dtype=float)
    if final_supplied.shape != (7, 3) or not np.isfinite(final_supplied).all():
        raise ValueError("Finite mapped final seven-atom XYZ required")
    sections = _sections(property_text)
    geometries = {index: _geometry(body) for (kind, index), body in sections.items() if kind == "Geometry"}
    energy_indices = {index for kind, index in sections if kind == "Single_Point_Data"}
    gradients = {index: _gradient(body) for (kind, index), body in sections.items() if kind == "SCF_Nuc_Gradient"}
    if not geometries or not gradients or not set(gradients) <= geometries.keys():
        raise ValueError("Each gradient must bind an existing property geometry")
    final_index, selected_index = max(geometries), max(gradients)
    if selected_index > final_index or not {selected_index, final_index} <= energy_indices:
        raise ValueError("Selected/final geometry lacks its same-index total energy")
    selected_xyz = geometries[selected_index]
    selected_A, final_A = selected_xyz * bohr_A, geometries[final_index] * bohr_A
    energy = _energy(sections[("Single_Point_Data", selected_index)])
    final_energy = _energy(sections[("Single_Point_Data", final_index)])
    gradient = gradients[selected_index]
    file_gap = float(np.max(abs(final_A - final_supplied)))
    if file_gap > _FINAL_FILE_TOL_A:
        raise ValueError("Final property geometry differs from final XYZ")

    marker = "THE OPTIMIZATION HAS CONVERGED"
    if (native_text.count(marker) != 1 or native_text.count("ORCA TERMINATED NORMALLY") != 1
            or native_text.index(marker) >= native_text.index("ORCA TERMINATED NORMALLY")
            or re.search(r"ORCA.*(?:error termination|TERMINATED ABNORMALLY)", native_text, re.I)):
        raise ValueError("Unique converged normal native optimization required")
    frames = _native_frames(native_text)
    native_energies = list(re.finditer(r"(?m)^\s*FINAL SINGLE POINT ENERGY\s+(" + _NUMBER + r")[ \t]*$", native_text))
    gradient_position, native_gradient = _native_last_gradient(native_text)
    earlier_frames = [item for item in frames if item[0] < gradient_position]
    earlier_energies = [item for item in native_energies if item.start() < gradient_position]
    if not earlier_frames or not earlier_energies or not native_energies:
        raise ValueError("Native evaluated gradient lacks its preceding frame/energy")
    if earlier_energies[-1].start() < earlier_frames[-1][0]:
        raise ValueError("Native evaluated energy does not belong to the gradient's preceding frame")
    cycle_markers = list(re.finditer(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)\b", native_text))
    cycle_numbers = [int(item[1]) for item in cycle_markers]
    if not cycle_numbers or cycle_numbers != sorted(set(cycle_numbers)):
        raise ValueError("Unique ordered native optimization cycles required")
    before_cycles = [item for item in cycle_markers if item.start() < gradient_position]
    if gradient_position < native_text.index(marker):
        if not before_cycles or int(before_cycles[-1][1]) != selected_index:
            raise ValueError("Property gradient index differs from its native evaluated cycle")
        cycle = before_cycles[-1]
        if (earlier_frames[-1][0] < cycle.end() or earlier_energies[-1].start() < earlier_frames[-1][0]
                or selected_index != cycle_numbers[-1]):
            raise ValueError("Last native gradient is not the final evaluated optimization cycle")
        phase = "last_evaluated_optimization_cycle"
        selected_native_cycle = selected_index
    else:
        if selected_index != final_index or earlier_frames[-1][0] < native_text.index(marker):
            raise ValueError("Post-optimization gradient lacks its final-geometry witness")
        phase = "post_optimization_final_evaluation"
        selected_native_cycle = None
    final_evaluations = re.findall(r"FINAL ENERGY EVALUATION AT THE STATIONARY POINT.*?\(AFTER\s+(\d+)\s+CYCLES\)",
                                   native_text, re.S)
    if final_evaluations != [str(cycle_numbers[-1])]:
        raise ValueError("Final native evaluation/cycle count missing or inconsistent")
    if frames[-1][0] < native_text.index(marker) or native_energies[-1].start() < frames[-1][0]:
        raise ValueError("Final native geometry/energy evaluation missing")
    selected_frame_gap = float(np.max(abs(earlier_frames[-1][1] - selected_A)))
    final_frame_gap = float(np.max(abs(frames[-1][1] - final_A)))
    gradient_gap = float(np.max(abs(native_gradient - gradient)))
    selected_energy_gap = abs(_finite(earlier_energies[-1][1]) - energy)
    final_energy_gap = abs(_finite(native_energies[-1][1]) - final_energy)
    if selected_frame_gap > _FRAME_PRINT_TOL_A or final_frame_gap > _FRAME_PRINT_TOL_A:
        raise ValueError("Property/native evaluated or final Cartesian frame differs")
    if gradient_gap > _GRADIENT_PRINT_TOL:
        raise ValueError("Property gradient differs from final native full gradient")
    if selected_energy_gap > _ENERGY_TOL or final_energy_gap > _ENERGY_TOL:
        raise ValueError("Same-index total energy differs from native evaluated/final energy")
    displacement = final_A - selected_A
    norms = np.linalg.norm(displacement, axis=1)
    receipt = {
        "kind": "ORCA_mapped_evaluated_gradient_frame", "schema_version": 1,
        "selected_geometry_index": selected_index, "final_geometry_index": final_index,
        "selected_frame_is_native_final": selected_index == final_index,
        "mapped_labels": list(LABELS), "xyz_bohr": selected_xyz.tolist(), "xyz_A": selected_A.tolist(),
        "energy_Eh": energy, "gradient_Eh_bohr": gradient.tolist(), "native_bohr_to_angstrom": float(bohr_A),
        "energy_scope": "Same-index Single_Point_Data.FinalEnergy, including dispersion; not bare SCF energy",
        "native_gradient_witness": {"phase": phase, "cycle": selected_native_cycle,
            "property_vs_native_frame_max_abs_A": selected_frame_gap,
            "property_vs_native_gradient_max_abs_Eh_bohr": gradient_gap,
            "property_vs_native_energy_abs_Eh": selected_energy_gap},
        "final_xyz_witness": {"xyz_bohr": geometries[final_index].tolist(), "xyz_A": final_A.tolist(),
            "supplied_final_xyz_A": final_supplied.tolist(), "energy_Eh": final_energy,
            "property_vs_file_max_abs_A": file_gap, "property_vs_native_frame_max_abs_A": final_frame_gap,
            "property_vs_native_energy_abs_Eh": final_energy_gap},
        "evaluated_vs_final_displacements": {"vectors_A": displacement.tolist(),
            "per_atom_A": norms.tolist(), "RMS_A": float(np.sqrt(np.mean(norms**2))),
            "maximum_atomic_A": float(norms.max()), "maximum_component_A": float(np.max(abs(displacement))),
            "alignment_or_permutation_applied": False},
        "scope": "Selected evaluated same-index geometry/E/g; do not attribute its gradient to a different final XYZ. Candidate stationarity, state, ownership and Hessian gates belong to caller.",
        "engrad_header_used_as_evaluation_identity": False, "saddle_character_verified": False,
        "rate_claim_allowed": False, "lifetime_claim_allowed": False, "native_jobs_dispatched": 0,
    }
    json.dumps(receipt, allow_nan=False)
    return receipt
