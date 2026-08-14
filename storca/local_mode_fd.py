"""Projected ORCA finite differences for local bond-stretch frequencies."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from typing import Iterable

import numpy as np

from src.inputgen import create_orca_gradient_input
from src.orca_runner import find_orca, run_orca, validate_orca_output
from .environment_refinement import parse_orca_engrad


LOCAL_MODE_SCHEMA_VERSION = 1
HARTREE_VIBRATIONAL_WAVENUMBER_CM_1 = 219474.6313632
ATOMIC_MASS_UNIT_IN_ELECTRON_MASSES = 1822.888486209
DEBYE_PER_ATOMIC_DIPOLE = 2.541746473
ANGSTROM_PER_BOHR = 0.529177210903
IR_INTENSITY_CONVERSION_KM_MOL = 42.255

# Standard atomic weights are sufficient for the intended local-mode mass
# model. Isotopologues need an explicit future mass contract.
ATOMIC_MASSES_U = {
    "H": 1.00794, "D": 2.014101778, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998403163, "Si": 28.085,
    "P": 30.973761998, "S": 32.06, "Cl": 35.45, "Br": 79.904,
    "I": 126.90447,
}


@dataclass(frozen=True)
class LocalModeFiniteDifferenceConfig:
    schema_version: int = LOCAL_MODE_SCHEMA_VERSION
    displacement_steps_angstrom: tuple[float, float] = (0.005, 0.010)
    maximum_frequency_step_disagreement_cm_1: float = 20.0
    maximum_negative_curvature_hartree_per_bohr2: float = -1.0e-8


@dataclass(frozen=True)
class LocalModeHessianValidationConfig:
    maximum_subspace_center_error_cm_1: float = 30.0
    maximum_total_intensity_relative_error: float = 0.35
    minimum_modes: int = 1


class OrcaInvocationLedger:
    """Persistent hard cap over actual ORCA process invocations in one stage."""

    def __init__(self, path: Path, hard_cap: int):
        self.path = Path(path)
        if hard_cap < 1:
            raise ValueError("ORCA invocation hard cap must be positive")
        if self.path.is_file():
            self.payload = json.loads(self.path.read_text())
            if int(self.payload["hard_cap"]) != int(hard_cap):
                raise RuntimeError("Retained ORCA ledger has a different hard cap")
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.payload = {
                "schema_version": 1,
                "kind": "orca_process_invocation_ledger",
                "hard_cap": int(hard_cap),
                "invocations": [],
            }
            self._write()

    def _write(self) -> None:
        self.payload["invocations_used"] = len(self.payload["invocations"])
        self.payload["invocations_remaining"] = (
            int(self.payload["hard_cap"]) - len(self.payload["invocations"])
        )
        self.path.write_text(json.dumps(self.payload, indent=2, sort_keys=True) + "\n")

    def is_completed(self, logical_job_id: str, input_sha256: str) -> bool:
        return any(
            item["logical_job_id"] == logical_job_id
            and item["input_sha256"] == input_sha256
            and item["status"] == "completed"
            for item in self.payload["invocations"]
        )

    def require_capacity(self, invocations: int) -> None:
        if invocations < 0:
            raise ValueError("Required ORCA invocation count cannot be negative")
        remaining = int(self.payload["hard_cap"]) - len(self.payload["invocations"])
        if invocations > remaining:
            raise RuntimeError(
                f"ORCA invocation budget exhausted: need {invocations}, have {remaining}"
            )

    def begin(self, logical_job_id: str, input_sha256: str) -> str | None:
        if self.is_completed(logical_job_id, input_sha256):
            return None
        self.require_capacity(1)
        attempt = 1 + sum(
            item["logical_job_id"] == logical_job_id
            for item in self.payload["invocations"]
        )
        invocation_id = f"{logical_job_id}:attempt-{attempt}"
        self.payload["invocations"].append({
            "invocation_id": invocation_id,
            "logical_job_id": logical_job_id,
            "input_sha256": input_sha256,
            "status": "running",
            "started_at": datetime.now(timezone.utc).isoformat(),
            "finished_at": None,
            "error": None,
        })
        self._write()
        return invocation_id

    def finish(self, invocation_id: str, *, error: str | None = None) -> None:
        matches = [
            item for item in self.payload["invocations"]
            if item["invocation_id"] == invocation_id
        ]
        if len(matches) != 1:
            raise RuntimeError(f"Unknown ORCA invocation: {invocation_id}")
        matches[0]["status"] = "failed" if error else "completed"
        matches[0]["finished_at"] = datetime.now(timezone.utc).isoformat()
        matches[0]["error"] = error
        self._write()


def _read_xyz(path: Path) -> tuple[list[str], np.ndarray, str]:
    lines = Path(path).read_text().splitlines()
    try:
        atom_count = int(lines[0])
    except (IndexError, ValueError) as error:
        raise ValueError(f"Invalid local-mode XYZ: {path}") from error
    rows = [line.split() for line in lines[2:2 + atom_count]]
    if len(rows) != atom_count or any(len(row) < 4 for row in rows):
        raise ValueError(f"Incomplete local-mode XYZ: {path}")
    symbols = [row[0] for row in rows]
    try:
        coordinates = np.asarray([[float(value) for value in row[1:4]] for row in rows])
    except ValueError as error:
        raise ValueError(f"Non-numeric local-mode XYZ: {path}") from error
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"Non-finite local-mode XYZ: {path}")
    comment = lines[1] if len(lines) > 1 else ""
    return symbols, coordinates, comment


def _write_xyz(path: Path, symbols: list[str], coordinates: np.ndarray, comment: str) -> None:
    rows = [str(len(symbols)), comment]
    rows.extend(
        f"{symbol:<2} {point[0]: .12f} {point[1]: .12f} {point[2]: .12f}"
        for symbol, point in zip(symbols, coordinates)
    )
    Path(path).write_text("\n".join(rows) + "\n")


def local_stretch_displacement(
    symbols: list[str], coordinates: np.ndarray, heavy_atom: int,
    hydrogen_atom: int, displacement_angstrom: float,
) -> tuple[np.ndarray, dict]:
    """Displace one bond length while preserving the pair center of mass."""
    if heavy_atom == hydrogen_atom or min(heavy_atom, hydrogen_atom) < 0:
        raise ValueError("Local stretch needs two distinct nonnegative atom indices")
    if max(heavy_atom, hydrogen_atom) >= len(symbols):
        raise ValueError("Local stretch atom index is outside the geometry")
    if symbols[hydrogen_atom] not in {"H", "D"}:
        raise ValueError("The local stretch hydrogen atom must be H or D")
    try:
        heavy_mass = ATOMIC_MASSES_U[symbols[heavy_atom]]
        hydrogen_mass = ATOMIC_MASSES_U[symbols[hydrogen_atom]]
    except KeyError as error:
        raise ValueError(f"No atomic mass available for {error.args[0]}") from error
    vector = coordinates[hydrogen_atom] - coordinates[heavy_atom]
    length = float(np.linalg.norm(vector))
    if length <= 1.0e-8:
        raise ValueError("Cannot displace a zero-length local bond")
    unit = vector / length
    total_mass = heavy_mass + hydrogen_mass
    displaced = np.asarray(coordinates, dtype=float).copy()
    displaced[hydrogen_atom] += displacement_angstrom * heavy_mass / total_mass * unit
    displaced[heavy_atom] -= displacement_angstrom * hydrogen_mass / total_mass * unit
    return displaced, {
        "unit_vector": unit.tolist(),
        "heavy_mass_u": heavy_mass,
        "hydrogen_mass_u": hydrogen_mass,
        "reduced_mass_u": heavy_mass * hydrogen_mass / total_mass,
        "initial_bond_length_angstrom": length,
    }


def parse_orca_dipole(path: Path) -> np.ndarray:
    """Return the last total dipole vector printed by ORCA, in atomic units."""
    pattern = re.compile(
        r"Total\s+Dipole\s+Moment\s*:\s*"
        r"([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)",
        re.IGNORECASE,
    )
    matches = pattern.findall(Path(path).read_text(errors="replace"))
    if not matches:
        raise ValueError(f"ORCA output contains no total dipole vector: {path}")
    vector = np.asarray([
        float(value.replace("D", "E").replace("d", "e"))
        for value in matches[-1]
    ])
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"ORCA output contains a non-finite dipole vector: {path}")
    return vector


def projected_local_mode_from_differences(
    plus_gradient: np.ndarray, minus_gradient: np.ndarray,
    plus_dipole_au: np.ndarray, minus_dipole_au: np.ndarray, *,
    unit_vector: np.ndarray, heavy_atom: int, hydrogen_atom: int,
    heavy_mass_u: float, hydrogen_mass_u: float, step_angstrom: float,
) -> dict:
    """Calculate projected curvature, frequency, and double-harmonic intensity."""
    gradients = [np.asarray(value, dtype=float) for value in (plus_gradient, minus_gradient)]
    if any(value.ndim == 1 for value in gradients):
        gradients = [value.reshape((-1, 3)) for value in gradients]
    if gradients[0].shape != gradients[1].shape or gradients[0].shape[1:] != (3,):
        raise ValueError("Projected local mode needs equal Cartesian gradient arrays")
    if max(heavy_atom, hydrogen_atom) >= len(gradients[0]) or step_angstrom <= 0:
        raise ValueError("Invalid local-mode atom index or displacement step")
    unit = np.asarray(unit_vector, dtype=float)
    if unit.shape != (3,) or not np.isclose(np.linalg.norm(unit), 1.0, atol=1e-8):
        raise ValueError("Local-mode direction must be a normalized Cartesian vector")
    total_mass = heavy_mass_u + hydrogen_mass_u
    heavy_fraction = hydrogen_mass_u / total_mass
    hydrogen_fraction = heavy_mass_u / total_mass

    def generalized_force(gradient: np.ndarray) -> float:
        return float(
            hydrogen_fraction * np.dot(gradient[hydrogen_atom], unit)
            - heavy_fraction * np.dot(gradient[heavy_atom], unit)
        )

    step_bohr = step_angstrom / ANGSTROM_PER_BOHR
    curvature = (generalized_force(gradients[0]) - generalized_force(gradients[1])) / (2 * step_bohr)
    reduced_mass_u = heavy_mass_u * hydrogen_mass_u / total_mass
    frequency = (
        HARTREE_VIBRATIONAL_WAVENUMBER_CM_1
        * math.sqrt(curvature / (reduced_mass_u * ATOMIC_MASS_UNIT_IN_ELECTRON_MASSES))
        if curvature > 0 else 0.0
    )
    dipole_derivative_au_per_bohr = (
        np.asarray(plus_dipole_au, dtype=float) - np.asarray(minus_dipole_au, dtype=float)
    ) / (2 * step_bohr)
    dipole_derivative_debye_per_angstrom = (
        dipole_derivative_au_per_bohr * DEBYE_PER_ATOMIC_DIPOLE / ANGSTROM_PER_BOHR
    )
    intensity = (
        IR_INTENSITY_CONVERSION_KM_MOL
        * float(np.dot(dipole_derivative_debye_per_angstrom, dipole_derivative_debye_per_angstrom))
        / reduced_mass_u
    )
    return {
        "step_angstrom": step_angstrom,
        "curvature_hartree_per_bohr2": curvature,
        "reduced_mass_u": reduced_mass_u,
        "frequency_cm-1": frequency,
        "dipole_derivative_au_per_bohr": dipole_derivative_au_per_bohr.tolist(),
        "intensity_km_mol": intensity,
    }


def _richardson_result(coarse: dict, fine: dict, maximum_disagreement: float) -> dict:
    if fine["step_angstrom"] >= coarse["step_angstrom"]:
        raise ValueError("Richardson local-mode steps must be supplied coarse then fine")
    extrapolated_curvature = (4.0 * fine["curvature_hartree_per_bohr2"] - coarse["curvature_hartree_per_bohr2"]) / 3.0
    reduced_mass = float(fine["reduced_mass_u"])
    frequency = (
        HARTREE_VIBRATIONAL_WAVENUMBER_CM_1
        * math.sqrt(extrapolated_curvature / (reduced_mass * ATOMIC_MASS_UNIT_IN_ELECTRON_MASSES))
        if extrapolated_curvature > 0 else 0.0
    )
    disagreement = abs(float(fine["frequency_cm-1"]) - float(coarse["frequency_cm-1"]))
    dipole = (
        4.0 * np.asarray(fine["dipole_derivative_au_per_bohr"])
        - np.asarray(coarse["dipole_derivative_au_per_bohr"])
    ) / 3.0
    dipole_debye_per_angstrom = dipole * DEBYE_PER_ATOMIC_DIPOLE / ANGSTROM_PER_BOHR
    intensity = IR_INTENSITY_CONVERSION_KM_MOL * float(np.dot(dipole_debye_per_angstrom, dipole_debye_per_angstrom)) / reduced_mass
    return {
        "frequency_cm-1": frequency,
        "intensity_km_mol": intensity,
        "curvature_hartree_per_bohr2": extrapolated_curvature,
        "frequency_step_disagreement_cm-1": disagreement,
        "displacement_stability_status": (
            "passed" if disagreement <= maximum_disagreement else "failed"
        ),
        "step_results": [coarse, fine],
    }


def validate_local_modes_against_harmonic_subspace(
    local_modes: Iterable[dict], harmonic_modes: Iterable[dict], *,
    config: LocalModeHessianValidationConfig | None = None,
) -> dict:
    """Validate localized oscillators against the corresponding coupled subspace.

    Individual localized bond stretches are not expected to equal individual
    symmetric/antisymmetric normal modes.  Their mean frequency and summed IR
    intensity should reproduce the same stretch subspace.
    """
    resolved = config or LocalModeHessianValidationConfig()
    local = [dict(item) for item in local_modes if item.get("status") == "validated"]
    harmonic = [dict(item) for item in harmonic_modes]
    failures = []
    if len(local) < resolved.minimum_modes or len(harmonic) < resolved.minimum_modes:
        failures.append("insufficient_modes")
    if len(local) != len(harmonic):
        failures.append("subspace_dimension_mismatch")
    local_center = float(np.mean([float(item["frequency_cm-1"]) for item in local])) if local else None
    harmonic_center = float(np.mean([float(item["frequency_cm-1"]) for item in harmonic])) if harmonic else None
    center_error = (
        abs(local_center - harmonic_center)
        if local_center is not None and harmonic_center is not None else None
    )
    local_intensity = sum(max(0.0, float(item["intensity_km_mol"])) for item in local)
    harmonic_intensity = sum(max(0.0, float(item["intensity_km_mol"])) for item in harmonic)
    intensity_error = (
        abs(local_intensity - harmonic_intensity) / harmonic_intensity
        if harmonic_intensity > 0 else None
    )
    if center_error is None or center_error > resolved.maximum_subspace_center_error_cm_1:
        failures.append("subspace_center_error_too_large")
    if intensity_error is None or intensity_error > resolved.maximum_total_intensity_relative_error:
        failures.append("subspace_intensity_error_too_large")
    return {
        "schema_version": LOCAL_MODE_SCHEMA_VERSION,
        "kind": "local_mode_stationary_hessian_subspace_validation",
        "status": "validated" if not failures else "failed_validation",
        "local_mode_count": len(local),
        "harmonic_mode_count": len(harmonic),
        "local_center_cm-1": local_center,
        "harmonic_subspace_center_cm-1": harmonic_center,
        "center_error_cm-1": center_error,
        "local_total_intensity_km_mol": local_intensity,
        "harmonic_subspace_total_intensity_km_mol": harmonic_intensity,
        "total_intensity_relative_error": intensity_error,
        "thresholds": asdict(resolved),
        "failures": failures,
    }


def calculate_orca_local_modes(
    xyz_path: Path, bonds: Iterable[dict], job_dir: Path, *, hard_orca_invocation_cap: int,
    charge: int = 0, multiplicity: int = 1, ncores: int = 1,
    method_keywords: list[str] | None = None,
    config: LocalModeFiniteDifferenceConfig | None = None,
    orca_ledger_path: Path | None = None,
    ledger_namespace: str = "",
    point_charge_file: Path | None = None,
) -> dict:
    """Run cached ±gradient calculations and return validated local modes."""
    resolved = config or LocalModeFiniteDifferenceConfig()
    if ncores > 1 and shutil.which("mpirun") is None:
        raise RuntimeError(
            "Parallel ORCA local modes require mpirun on PATH; refusing to consume "
            "an ORCA invocation before the MPI runtime is configured"
        )
    steps = tuple(sorted({float(value) for value in resolved.displacement_steps_angstrom}))
    if len(steps) != 2 or steps[0] <= 0:
        raise ValueError("Local-mode validation requires two positive displacement sizes")
    xyz_path, job_dir = Path(xyz_path), Path(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)
    source = xyz_path.read_bytes()
    point_charge_source = Path(point_charge_file).read_bytes() if point_charge_file is not None else None
    input_xyz = job_dir / "input.xyz"
    if input_xyz.is_file() and input_xyz.read_bytes() != source:
        raise RuntimeError("Refusing to replace a different retained local-mode input")
    input_xyz.write_bytes(source)
    symbols, coordinates, _ = _read_xyz(input_xyz)
    bonds = [dict(item) for item in bonds]
    contract = {
        "schema_version": LOCAL_MODE_SCHEMA_VERSION,
        "kind": "orca_projected_local_mode_finite_differences",
        "input_xyz_sha256": hashlib.sha256(source).hexdigest(),
        "orca_executable": str(Path(find_orca()).resolve()),
        "charge": int(charge), "multiplicity": int(multiplicity), "ncores": int(ncores),
        "method_keywords": list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]),
        "bonds": bonds,
        "configuration": {
            **asdict(resolved),
            "displacement_steps_angstrom": list(resolved.displacement_steps_angstrom),
        },
        "orca_invocations_required": 4 * len(bonds),
        "electrostatic_embedding": (
            {
                "point_charge_sha256": hashlib.sha256(point_charge_source).hexdigest(),
                "point_charge_count": int(point_charge_source.splitlines()[0]),
            } if point_charge_source is not None else None
        ),
    }
    contract_path, result_path = job_dir / "calculation-contract.json", job_dir / "local-modes.json"
    if contract_path.is_file() and json.loads(contract_path.read_text()) != contract:
        raise RuntimeError("Retained local-mode calculation has a different contract")
    if contract_path.is_file() and result_path.is_file():
        return json.loads(result_path.read_text())
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    ledger = OrcaInvocationLedger(
        orca_ledger_path or (job_dir / "orca-invocations.json"),
        hard_orca_invocation_cap,
    )
    logical_jobs = []
    for bond_index, bond in enumerate(bonds):
        for step in steps:
            for sign in (-1, 1):
                logical_jobs.append((bond_index, bond, step, sign))
    unfinished = 0
    for bond_index, _, step, sign in logical_jobs:
        prefix = f"{ledger_namespace}:" if ledger_namespace else ""
        logical_id = f"{prefix}bond-{bond_index:03d}:step-{step:.6f}:{'plus' if sign > 0 else 'minus'}"
        existing_input = job_dir / f"bond-{bond_index:03d}" / f"step-{step:.6f}-{'plus' if sign > 0 else 'minus'}" / "gradient.inp"
        input_hash = hashlib.sha256(existing_input.read_bytes()).hexdigest() if existing_input.is_file() else "not_generated"
        if not ledger.is_completed(logical_id, input_hash):
            unfinished += 1
    ledger.require_capacity(unfinished)

    collected: dict[tuple[int, float, int], tuple[np.ndarray, np.ndarray]] = {}
    for bond_index, bond, step, sign in logical_jobs:
        heavy, hydrogen = int(bond["heavy_atom"]), int(bond["hydrogen_atom"])
        displaced, displacement = local_stretch_displacement(
            symbols, coordinates, heavy, hydrogen, sign * step,
        )
        sign_label = "plus" if sign > 0 else "minus"
        calculation_dir = job_dir / f"bond-{bond_index:03d}" / f"step-{step:.6f}-{sign_label}"
        calculation_dir.mkdir(parents=True, exist_ok=True)
        geometry = calculation_dir / "geometry.xyz"
        _write_xyz(geometry, symbols, displaced, f"local mode bond {heavy}-{hydrogen} {sign * step:+.6f} A")
        calculation_point_charges = None
        if point_charge_source is not None:
            calculation_point_charges = calculation_dir / "embedding.pc"
            if calculation_point_charges.is_file() and calculation_point_charges.read_bytes() != point_charge_source:
                raise RuntimeError("Refusing to replace different retained embedding point charges")
            calculation_point_charges.write_bytes(point_charge_source)
        gradient_input = create_orca_gradient_input(
            geometry, charge=charge, multiplicity=multiplicity, label="gradient",
            ncores=ncores, method_keywords=method_keywords,
            point_charge_file=calculation_point_charges,
        )
        input_hash = hashlib.sha256(gradient_input.read_bytes()).hexdigest()
        prefix = f"{ledger_namespace}:" if ledger_namespace else ""
        logical_id = f"{prefix}bond-{bond_index:03d}:step-{step:.6f}:{sign_label}"
        output, engrad = gradient_input.with_suffix(".out"), gradient_input.with_suffix(".engrad")
        invocation_id = ledger.begin(logical_id, input_hash)
        if invocation_id is not None:
            try:
                run_orca(gradient_input)
                ledger.finish(invocation_id)
            except Exception as error:
                ledger.finish(invocation_id, error=str(error))
                raise
        if not engrad.is_file() or not validate_orca_output(output)["normal_termination"]:
            raise RuntimeError(f"Incomplete retained local-mode gradient: {calculation_dir}")
        gradient = np.asarray(parse_orca_engrad(engrad)["gradient_hartree_per_bohr"])
        collected[(bond_index, step, sign)] = (gradient, parse_orca_dipole(output))

    mode_results = []
    for bond_index, bond in enumerate(bonds):
        heavy, hydrogen = int(bond["heavy_atom"]), int(bond["hydrogen_atom"])
        _, displacement = local_stretch_displacement(symbols, coordinates, heavy, hydrogen, steps[0])
        per_step = []
        for step in reversed(steps):
            plus_gradient, plus_dipole = collected[(bond_index, step, 1)]
            minus_gradient, minus_dipole = collected[(bond_index, step, -1)]
            per_step.append(projected_local_mode_from_differences(
                plus_gradient, minus_gradient, plus_dipole, minus_dipole,
                unit_vector=np.asarray(displacement["unit_vector"]),
                heavy_atom=heavy, hydrogen_atom=hydrogen,
                heavy_mass_u=displacement["heavy_mass_u"],
                hydrogen_mass_u=displacement["hydrogen_mass_u"],
                step_angstrom=step,
            ))
        result = _richardson_result(
            per_step[0], per_step[1],
            resolved.maximum_frequency_step_disagreement_cm_1,
        )
        curvature_ok = result["curvature_hartree_per_bohr2"] > resolved.maximum_negative_curvature_hartree_per_bohr2
        result.update({
            "bond": bond,
            "heavy_atom": heavy,
            "hydrogen_atom": hydrogen,
            "status": (
                "validated" if curvature_ok and result["frequency_cm-1"] > 0
                and result["displacement_stability_status"] == "passed"
                else "failed_validation"
            ),
        })
        mode_results.append(result)
    payload = {
        "schema_version": LOCAL_MODE_SCHEMA_VERSION,
        "kind": "orca_projected_local_modes",
        "status": "completed" if mode_results and all(item["status"] == "validated" for item in mode_results) else "failed_validation",
        "population_interpretation": "representative_local_properties_not_populations",
        "modes": mode_results,
        "calculation_contract": str(contract_path),
        "orca_invocation_ledger": str(ledger.path),
    }
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload
