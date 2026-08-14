"""Environment-preserving DFT refinement and unrestrained gradient gating."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np

from src.inputgen import (create_orca_environment_refinement_input,
                          create_orca_gradient_input)
from src.orca_runner import find_orca, run_orca


REFINEMENT_SCHEMA_VERSION = 2


@dataclass(frozen=True)
class EnvironmentRefinementConfig:
    schema_version: int = REFINEMENT_SCHEMA_VERSION
    maximum_iterations: int = 120
    constrain_angles: bool = True
    maximum_gradient_rms_hartree_per_bohr: float = 3.0e-4
    maximum_gradient_component_hartree_per_bohr: float = 1.0e-3


def parse_orca_engrad(path: Path) -> dict:
    """Parse ORCA's machine-readable ``.engrad`` gradient file."""
    values = [
        line.strip() for line in Path(path).read_text(errors="replace").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    try:
        atom_count = int(values[0])
        energy = float(values[1].replace("D", "E").replace("d", "e"))
        gradient = np.asarray([
            float(value.replace("D", "E").replace("d", "e"))
            for value in values[2:2 + 3 * atom_count]
        ])
    except (IndexError, ValueError) as error:
        raise ValueError(f"Malformed ORCA gradient file: {path}") from error
    if atom_count < 1 or len(gradient) != 3 * atom_count or not np.all(np.isfinite(gradient)):
        raise ValueError(f"Incomplete or non-finite ORCA gradient file: {path}")
    return {
        "atom_count": atom_count,
        "energy_hartree": energy,
        "gradient_hartree_per_bohr": gradient.reshape((atom_count, 3)).tolist(),
        "gradient_rms_hartree_per_bohr": float(np.sqrt(np.mean(np.square(gradient)))),
        "gradient_maximum_component_hartree_per_bohr": float(np.max(np.abs(gradient))),
    }


def _read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = Path(path).read_text().splitlines()
    count = int(lines[0])
    rows = [line.split() for line in lines[2:2 + count]]
    if len(rows) != count or any(len(row) < 4 for row in rows):
        raise ValueError(f"Invalid refinement XYZ: {path}")
    return [row[0] for row in rows], np.asarray([
        [float(value) for value in row[1:4]] for row in rows
    ])


def _angle(first: np.ndarray, vertex: np.ndarray, third: np.ndarray) -> float:
    left, right = first - vertex, third - vertex
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator <= 1e-12:
        raise ValueError("Cannot constrain an interaction containing coincident atoms")
    return math.degrees(math.acos(float(np.clip(np.dot(left, right) / denominator, -1.0, 1.0))))


def _resolved_interactions(xyz_path: Path, interactions: list[dict]) -> list[dict]:
    _, coordinates = _read_xyz(xyz_path)
    output = []
    for item in interactions:
        donor, hydrogen, acceptor = (
            int(item["donor_atom"]), int(item["donor_hydrogen"]), int(item["acceptor_atom"])
        )
        if min(donor, hydrogen, acceptor) < 0 or max(donor, hydrogen, acceptor) >= len(coordinates):
            raise ValueError("Environment interaction atom is outside the refinement geometry")
        output.append({
            "donor_atom": donor, "donor_hydrogen": hydrogen, "acceptor_atom": acceptor,
            "h_bond_distance_angstrom": float(np.linalg.norm(coordinates[hydrogen] - coordinates[acceptor])),
            "donor_h_acceptor_angle_degrees": _angle(
                coordinates[donor], coordinates[hydrogen], coordinates[acceptor],
            ),
        })
    if not output:
        raise ValueError("Environment refinement requires at least one H-bond interaction")
    return output


def refine_orca_environment(
    xyz_path: Path, interactions: list[dict], job_dir: Path, *, charge: int = 0,
    multiplicity: int = 1, ncores: int = 1,
    method_keywords: list[str] | None = None,
    config: EnvironmentRefinementConfig | None = None,
) -> dict:
    """Refine one environment, then gate full-Hessian use by an unrestrained gradient."""
    resolved = config or EnvironmentRefinementConfig()
    xyz_path, job_dir = Path(xyz_path), Path(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)
    input_path = job_dir / "input.xyz"
    source = xyz_path.read_bytes()
    if input_path.is_file() and input_path.read_bytes() != source:
        raise RuntimeError(f"Refusing to replace a different environment refinement input: {input_path}")
    input_path.write_bytes(source)
    resolved_interactions = _resolved_interactions(input_path, interactions)
    contract = {
        "schema_version": REFINEMENT_SCHEMA_VERSION,
        "input_xyz_sha256": hashlib.sha256(source).hexdigest(),
        "orca_executable": str(Path(find_orca()).resolve()),
        "method_keywords": list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]),
        "charge": int(charge), "multiplicity": int(multiplicity), "ncores": int(ncores),
        "optimization_constraints": resolved_interactions,
        "frequency_constraints": None,
        "gradient_constraints": None,
        "configuration": asdict(resolved),
    }
    contract_path, record_path = job_dir / "calculation-contract.json", job_dir / "refinement-record.json"
    retained_contract = json.loads(contract_path.read_text()) if contract_path.is_file() else None
    retained_record = json.loads(record_path.read_text()) if record_path.is_file() else None
    if retained_contract is not None and retained_contract != contract:
        raise RuntimeError(f"Retained environment refinement has a different contract: {job_dir}")
    refined_path = job_dir / "environment-refine.xyz"
    gradient_path = job_dir / "unrestrained-gradient.engrad"
    if retained_contract == contract and retained_record and refined_path.is_file() and gradient_path.is_file():
        return retained_record
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    try:
        refinement_input = create_orca_environment_refinement_input(
            input_path, interactions=resolved_interactions, charge=charge,
            multiplicity=multiplicity, ncores=ncores,
            max_iterations=resolved.maximum_iterations,
            constrain_angles=resolved.constrain_angles,
            method_keywords=method_keywords,
        )
        refinement = run_orca(refinement_input)
        if not Path(refinement["xyz"]).is_file():
            raise RuntimeError("Environment refinement completed without a refined XYZ")
        gradient_input = create_orca_gradient_input(
            Path(refinement["xyz"]), charge=charge, multiplicity=multiplicity,
            label="unrestrained-gradient", ncores=ncores, method_keywords=method_keywords,
        )
        run_orca(gradient_input)
        gradient = parse_orca_engrad(gradient_path)
        full_hessian_permitted = bool(
            gradient["gradient_rms_hartree_per_bohr"] <= resolved.maximum_gradient_rms_hartree_per_bohr
            and gradient["gradient_maximum_component_hartree_per_bohr"]
            <= resolved.maximum_gradient_component_hartree_per_bohr
        )
        record = {
            "schema_version": REFINEMENT_SCHEMA_VERSION,
            "status": "completed",
            "refined_xyz": str(refined_path),
            "gradient": gradient,
            "stationarity_status": "usable" if full_hessian_permitted else "poor",
            "full_hessian_use": "permitted" if full_hessian_permitted else "diagnostic_only",
            "local_mode_policy": (
                "full_hessian_allowed" if full_hessian_permitted
                else "use_local_mode_finite_differences_not_full_hessian"
            ),
            "error": None,
        }
    except Exception as error:
        record = {
            "schema_version": REFINEMENT_SCHEMA_VERSION, "status": "failed",
            "refined_xyz": None, "gradient": None,
            "stationarity_status": "unknown", "full_hessian_use": "not_permitted",
            "local_mode_policy": "not_available", "error": str(error),
        }
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def representative_vibrational_route(refinement_record: dict) -> dict:
    """Resolve the only permitted vibrational treatment after refinement."""
    if (
        refinement_record.get("status") == "completed"
        and refinement_record.get("full_hessian_use") == "permitted"
        and refinement_record.get("stationarity_status") == "usable"
    ):
        return {
            "route": "full_orca_hessian",
            "full_hessian_use": "permitted",
            "reason": "unrestrained_gradient_gate_passed",
        }
    return {
        "route": "orca_projected_local_mode_finite_differences",
        "full_hessian_use": "not_permitted",
        "reason": "unrestrained_gradient_gate_not_passed",
    }


def refine_selected_orca_environments(
    run_dir: Path, *, charge: int = 0, multiplicity: int = 1, ncores: int = 1,
    method_keywords: list[str] | None = None, progress=None,
) -> tuple[list[Path], Path]:
    """Refine selected environments and update their execution manifest in place."""
    run_dir = Path(run_dir)
    manifest_path = run_dir / "clusters" / "selected-conformers.json"
    if not manifest_path.is_file():
        raise RuntimeError("Environment refinement requires selected representatives")
    manifest = json.loads(manifest_path.read_text())
    entries = sorted(
        manifest.get("conformers") or [], key=lambda item: int(item["selected_position"]),
    )
    if not entries:
        raise RuntimeError("Environment refinement manifest has no representatives")
    records = []
    refined_paths = []
    refinement_root = run_dir / "clusters" / "environment-refinements"
    for position, entry in enumerate(entries, start=1):
        candidate_id = str(entry["source_xtb_candidate_id"])
        if progress:
            progress(
                f"Environment DFT refinement {position}/{len(entries)}: {candidate_id}"
            )
        source_path = Path(entry.get("source_xtb_xyz") or entry["xyz"])
        record = refine_orca_environment(
            source_path, entry.get("hydrogen_bond_interactions") or [],
            refinement_root / candidate_id, charge=charge,
            multiplicity=multiplicity, ncores=ncores,
            method_keywords=method_keywords,
        )
        route = representative_vibrational_route(record)
        refined_path = Path(record["refined_xyz"]) if record.get("refined_xyz") else source_path
        entry.update({
            "xyz": str(refined_path),
            "pre_refinement_xyz": str(source_path),
            "environment_refinement": record,
            "environment_refinement_artifact": str(
                refinement_root / candidate_id / "refinement-record.json"
            ),
            "vibrational_route": route,
            "gradient_rms_hartree_per_bohr": (
                (record.get("gradient") or {}).get("gradient_rms_hartree_per_bohr")
            ),
            "gradient_maximum_component_hartree_per_bohr": (
                (record.get("gradient") or {}).get(
                    "gradient_maximum_component_hartree_per_bohr"
                )
            ),
        })
        records.append({
            "source_xtb_candidate_id": candidate_id,
            "input_xyz": str(source_path),
            "execution_xyz": str(refined_path),
            "refinement": record,
            "vibrational_route": route,
        })
        refined_paths.append(refined_path)
    manifest["schema_version"] = max(3, int(manifest.get("schema_version", 0)))
    manifest["environment_refinement"] = {
        "kind": "constrained_dft_refinement_then_unrestrained_gradient_gate",
        "requested": len(entries),
        "completed": sum(item["refinement"].get("status") == "completed" for item in records),
        "full_hessian_permitted": sum(
            item["vibrational_route"].get("full_hessian_use") == "permitted"
            for item in records
        ),
        "records": records,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    report_path = run_dir / "environment-refinement.json"
    report_path.write_text(json.dumps({
        "schema_version": REFINEMENT_SCHEMA_VERSION,
        "kind": "selected_environment_dft_refinement",
        "manifest": str(manifest_path),
        "records": records,
    }, indent=2, sort_keys=True) + "\n")
    return refined_paths, report_path
