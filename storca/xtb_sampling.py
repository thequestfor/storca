"""Resumable, restraint-preserving GFN2-xTB environment sampling."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Callable

import numpy as np

from src.orca_runner import find_xtb

from .cluster_ir import (classify_local_stretch_bonds, generate_stratified_dimer_poses,
                         generate_stratified_trimer_poses)


HARTREE_TO_KCAL_MOL = 627.509474
XTB_SAMPLING_SCHEMA_VERSION = 1
DEFAULT_RESTRAINT_FORCE_CONSTANT = 0.05


def xtb_sampling_defaults(fidelity: str) -> dict:
    fidelity = fidelity.strip().lower()
    candidate_counts = {"fast": 0, "auto": 40, "balanced": 75, "accurate": 100}
    if fidelity not in candidate_counts:
        raise ValueError("Unsupported xTB sampling fidelity")
    return {
        "schema_version": XTB_SAMPLING_SCHEMA_VERSION,
        "fidelity": fidelity,
        "candidate_count": candidate_counts[fidelity],
        "method": "GFN2-xTB",
        "optimization_level": "loose",
        "restraint_force_constant": DEFAULT_RESTRAINT_FORCE_CONSTANT,
        "population_model": "not_assigned",
        "population_warning": (
            "Stratified restrained geometries provide coverage, not condensed-phase populations."
        ),
    }


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _json_write(path: Path, value: dict) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return path


def _write_exact(path: Path, text: str) -> Path:
    path = Path(path)
    if path.is_file() and path.read_text() != text:
        raise RuntimeError(f"Refusing to overwrite a different retained sampling input: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _xyz_text(symbols: list[str], coordinates: np.ndarray, comment: str) -> str:
    return "\n".join([
        str(len(symbols)),
        comment,
        *(
            f"{symbol} {position[0]:.10f} {position[1]:.10f} {position[2]:.10f}"
            for symbol, position in zip(symbols, coordinates)
        ),
        "",
    ])


def _read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = Path(path).read_text().splitlines()
    try:
        atom_count = int(lines[0])
    except (IndexError, ValueError) as error:
        raise ValueError(f"Invalid XYZ file: {path}") from error
    rows = [line.split() for line in lines[2:2 + atom_count]]
    if len(rows) != atom_count or any(len(row) < 4 for row in rows):
        raise ValueError(f"Incomplete XYZ file: {path}")
    symbols = [row[0] for row in rows]
    try:
        coordinates = np.asarray([[float(value) for value in row[1:4]] for row in rows], dtype=float)
    except ValueError as error:
        raise ValueError(f"Non-numeric XYZ coordinates: {path}") from error
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"Non-finite XYZ coordinates: {path}")
    return symbols, coordinates


def _angle_degrees(first: np.ndarray, vertex: np.ndarray, third: np.ndarray) -> float:
    left = first - vertex
    right = third - vertex
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator < 1e-12:
        raise ValueError("Cannot calculate an angle containing coincident atoms")
    cosine = float(np.dot(left, right) / denominator)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _dihedral_degrees(coordinates: np.ndarray, atoms: tuple[int, int, int, int]) -> float | None:
    first, second, third, fourth = (coordinates[index] for index in atoms)
    bond1 = second - first
    bond2 = third - second
    bond3 = fourth - third
    normal1 = np.cross(bond1, bond2)
    normal2 = np.cross(bond2, bond3)
    if min(np.linalg.norm(normal1), np.linalg.norm(normal2), np.linalg.norm(bond2)) < 1e-10:
        return None
    normal1 /= np.linalg.norm(normal1)
    normal2 /= np.linalg.norm(normal2)
    direction = bond2 / np.linalg.norm(bond2)
    value = math.degrees(math.atan2(np.dot(np.cross(normal1, normal2), direction), np.dot(normal1, normal2)))
    return value % 360.0


def _aligned_heavy_atom_rmsd(
    initial: np.ndarray, optimized: np.ndarray, symbols: list[str],
) -> float:
    indices = [index for index, symbol in enumerate(symbols) if symbol != "H"]
    if not indices:
        indices = list(range(len(symbols)))
    reference = initial[indices] - initial[indices].mean(axis=0)
    candidate = optimized[indices] - optimized[indices].mean(axis=0)
    covariance = candidate.T @ reference
    left, _, right = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right))
    rotation = left @ correction @ right
    aligned = candidate @ rotation
    return float(np.sqrt(np.mean(np.square(aligned - reference))))


def _xtb_version(executable: str) -> str:
    result = subprocess.run(
        [executable, "--version"], capture_output=True, text=True, timeout=30,
    )
    text = f"{result.stdout}\n{result.stderr}"
    match = re.search(r"xtb version\s+([^\s]+)", text, flags=re.IGNORECASE)
    if result.returncode != 0 or not match:
        raise RuntimeError(f"Could not determine xTB version from {executable}")
    return match.group(1)


def _xcontrol_text(pose: dict, force_constant: float) -> str:
    if not math.isfinite(force_constant) or force_constant <= 0:
        raise ValueError("xTB restraint force constant must be finite and positive")
    lines = [
        "$constrain",
        f"  force constant={force_constant:.8f}",
    ]
    interactions = pose.get("interactions") or [{
        "donor_atom": pose["donor_atom"], "donor_hydrogen": pose["donor_hydrogen"],
        "acceptor_atom": pose["acceptor_atom"], "orientation_atoms": pose.get("orientation_atoms"),
        "target_geometry": pose["target_geometry"],
    }]
    for interaction in interactions:
        donor = int(interaction["donor_atom"]) + 1
        hydrogen = int(interaction["donor_hydrogen"]) + 1
        acceptor = int(interaction["acceptor_atom"]) + 1
        target = interaction["target_geometry"]
        lines.append(f"  distance: {hydrogen},{acceptor},{float(target['h_bond_distance_angstrom']):.8f}")
        lines.append(f"  angle: {donor},{hydrogen},{acceptor},{float(target['donor_h_acceptor_angle_degrees']):.8f}")
        orientation_atoms = interaction.get("orientation_atoms")
        if orientation_atoms:
            target_dihedral = _dihedral_degrees(
                np.asarray(pose["coordinates"], dtype=float), tuple(orientation_atoms),
            )
            if target_dihedral is not None:
                atoms = ",".join(str(int(index) + 1) for index in orientation_atoms)
                lines.append(f"  dihedral: {atoms},{target_dihedral:.8f}")
    lines.extend(["$end", ""])
    return "\n".join(lines)


def _extract_json_energy(value) -> float | None:
    if isinstance(value, dict):
        normalized = {str(key).strip().lower().replace("_", " "): item for key, item in value.items()}
        for key in ("total energy", "totalenergy", "energy"):
            candidate = normalized.get(key)
            if isinstance(candidate, (int, float)) and math.isfinite(candidate):
                return float(candidate)
        for item in value.values():
            candidate = _extract_json_energy(item)
            if candidate is not None:
                return candidate
    return None


def _parse_xtb_energy(json_path: Path, output_text: str) -> float:
    if json_path.is_file():
        try:
            energy = _extract_json_energy(json.loads(json_path.read_text()))
        except (OSError, ValueError, TypeError):
            energy = None
        if energy is not None:
            return energy
    matches = re.findall(
        r"TOTAL ENERGY\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)", output_text,
        flags=re.IGNORECASE,
    )
    if not matches:
        raise ValueError("xTB output did not contain a finite total energy")
    energy = float(matches[-1])
    if not math.isfinite(energy):
        raise ValueError("xTB output contained a non-finite total energy")
    return energy


def _validate_optimized_pose(pose: dict, optimized_path: Path) -> dict:
    symbols, optimized = _read_xyz(optimized_path)
    initial = np.asarray(pose["coordinates"], dtype=float)
    if symbols != list(pose["symbols"]) or optimized.shape != initial.shape:
        raise ValueError("Optimized xTB geometry changed atom identities or ordering")
    broken_bonds = []
    for left, right in pose["bonds"]:
        initial_length = float(np.linalg.norm(initial[left] - initial[right]))
        optimized_length = float(np.linalg.norm(optimized[left] - optimized[right]))
        if optimized_length < 0.6 * initial_length or optimized_length > 1.8 * initial_length:
            broken_bonds.append([left, right])
    atom_count = int(pose["monomer_atom_count"])
    cluster_size = int(pose.get("cluster_size", len(optimized) // atom_count))
    inter_distances = []
    for left in range(cluster_size):
        for right in range(left + 1, cluster_size):
            inter_distances.append(np.linalg.norm(
                optimized[left * atom_count:(left + 1) * atom_count, None, :]
                - optimized[None, right * atom_count:(right + 1) * atom_count, :], axis=2,
            ))
    minimum_inter_distance = min(float(values.min()) for values in inter_distances)
    interactions = pose.get("interactions") or [{
        "donor_atom": pose["donor_atom"], "donor_hydrogen": pose["donor_hydrogen"],
        "acceptor_atom": pose["acceptor_atom"], "orientation_atoms": pose.get("orientation_atoms"),
        "target_geometry": pose["target_geometry"],
    }]
    interaction_results = []
    for interaction in interactions:
        donor = int(interaction["donor_atom"])
        hydrogen = int(interaction["donor_hydrogen"])
        acceptor = int(interaction["acceptor_atom"])
        target = interaction["target_geometry"]
        distance = float(np.linalg.norm(optimized[hydrogen] - optimized[acceptor]))
        angle = _angle_degrees(optimized[donor], optimized[hydrogen], optimized[acceptor])
        interaction_results.append({
            "h_bond_distance_angstrom": distance,
            "donor_h_acceptor_angle_degrees": angle,
            "donor_acceptor_distance_angstrom": float(np.linalg.norm(optimized[donor] - optimized[acceptor])),
            "intermolecular_orientation_degrees": (
                _dihedral_degrees(optimized, tuple(interaction["orientation_atoms"]))
                if interaction.get("orientation_atoms") else None
            ),
            "distance_deviation_angstrom": abs(distance - float(target["h_bond_distance_angstrom"])),
            "angle_deviation_degrees": abs(angle - float(target["donor_h_acceptor_angle_degrees"])),
        })
    failures = []
    if broken_bonds:
        failures.append("intramolecular_connectivity_changed")
    if minimum_inter_distance < 0.65:
        failures.append("severe_intermolecular_collision")
    if any(not 1.2 <= item["h_bond_distance_angstrom"] <= 3.2 for item in interaction_results):
        failures.append("hydrogen_bond_contact_dissociated_or_covalent")
    if any(item["distance_deviation_angstrom"] > 0.35 for item in interaction_results):
        failures.append("distance_left_sampled_stratum")
    if any(item["angle_deviation_degrees"] > 25.0 for item in interaction_results):
        failures.append("angle_left_sampled_stratum")
    distances = [item["h_bond_distance_angstrom"] for item in interaction_results]
    angles = [item["donor_h_acceptor_angle_degrees"] for item in interaction_results]
    donor_acceptor_distances = [item["donor_acceptor_distance_angstrom"] for item in interaction_results]
    orientations = [item["intermolecular_orientation_degrees"] for item in interaction_results if item["intermolecular_orientation_degrees"] is not None]
    features = {
        "geometry_cluster_id": None,
        "cluster_size": cluster_size,
        "topology": pose.get("topology", "dimer"),
        "h_bond_distance_angstrom": float(np.mean(distances)),
        "minimum_h_bond_distance_angstrom": min(distances),
        "maximum_h_bond_distance_angstrom": max(distances),
        "h_bond_distance_span_angstrom": max(distances) - min(distances),
        "donor_h_acceptor_angle_degrees": float(np.mean(angles)),
        "minimum_h_bond_angle_degrees": min(angles),
        "maximum_h_bond_angle_degrees": max(angles),
        "donor_acceptor_distance_angstrom": float(np.mean(donor_acceptor_distances)),
        "intermolecular_orientation_degrees": orientations[0] if orientations else None,
        "heavy_atom_rmsd_angstrom": _aligned_heavy_atom_rmsd(initial, optimized, symbols),
        "heavy_atom_rmsd_basis": "initial_to_xtb_relaxed_after_global_alignment",
        "local_electric_field": None,
        "estimated_local_frequency_cm-1": None,
    }
    return {
        "valid": not failures,
        "failures": failures,
        "broken_bonds": broken_bonds,
        "minimum_inter_distance_angstrom": minimum_inter_distance,
        "distance_deviation_angstrom": max(item["distance_deviation_angstrom"] for item in interaction_results),
        "angle_deviation_degrees": max(item["angle_deviation_degrees"] for item in interaction_results),
        "interactions": interaction_results,
        "features": features,
    }


def run_xtb_candidate(
    pose: dict,
    candidate_dir: Path,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    ncores: int = 1,
    force_constant: float = DEFAULT_RESTRAINT_FORCE_CONSTANT,
    executable: str | None = None,
    executable_version: str | None = None,
    timeout_seconds: float = 600.0,
    retry_failed: bool = False,
) -> dict:
    """Run or reuse one isolated restrained xTB optimization."""
    if multiplicity < 1 or ncores < 1:
        raise ValueError("Multiplicity and xTB core count must be positive")
    candidate_dir = Path(candidate_dir)
    candidate_dir.mkdir(parents=True, exist_ok=True)
    initial_text = _xyz_text(
        list(pose["symbols"]), np.asarray(pose["coordinates"], dtype=float),
        f"{pose['candidate_id']} stratified hydrogen-bond environment",
    )
    input_path = _write_exact(candidate_dir / "initial.xyz", initial_text)
    xcontrol_text = _xcontrol_text(pose, force_constant)
    xcontrol_path = _write_exact(candidate_dir / "xcontrol.inp", xcontrol_text)
    executable = executable or find_xtb()
    executable = str(Path(executable).resolve())
    executable_version = executable_version or _xtb_version(executable)
    contract = {
        "schema_version": XTB_SAMPLING_SCHEMA_VERSION,
        "candidate_id": pose["candidate_id"],
        "input_xyz_sha256": _sha256_text(initial_text),
        "xcontrol_sha256": _sha256_text(xcontrol_text),
        "xtb_executable": executable,
        "xtb_version": executable_version,
        "method": "GFN2-xTB",
        "optimization_level": "loose",
        "charge": int(charge),
        "unpaired_electrons": int(multiplicity - 1),
        "ncores": int(ncores),
    }
    contract_path = candidate_dir / "calculation-contract.json"
    record_path = candidate_dir / "sampling-record.json"
    retained_contract = json.loads(contract_path.read_text()) if contract_path.is_file() else None
    retained_record = json.loads(record_path.read_text()) if record_path.is_file() else None
    if retained_contract is not None and retained_contract != contract:
        raise RuntimeError(f"Retained xTB candidate has a different calculation contract: {candidate_dir}")
    optimized_path = candidate_dir / "optimized.xyz"
    if retained_contract == contract and retained_record:
        if retained_record.get("optimization_status") == "completed" and optimized_path.is_file():
            return retained_record
        if retained_record.get("optimization_status") == "failed" and not retry_failed:
            return retained_record
    _json_write(contract_path, contract)
    output_path = candidate_dir / "xtb.out"
    runtime_environment = os.environ.copy()
    runtime_environment.update({
        "OMP_NUM_THREADS": str(ncores),
        "MKL_NUM_THREADS": str(ncores),
        "OMP_STACKSIZE": runtime_environment.get("OMP_STACKSIZE", "1G"),
    })
    command = [
        executable, input_path.name, "--gfn", "2", "--opt", "loose",
        "--input", xcontrol_path.name, "--chrg", str(charge),
        "--uhf", str(multiplicity - 1), "--json", "--norestart",
    ]
    base_record = {
        "schema_version": XTB_SAMPLING_SCHEMA_VERSION,
        "candidate_id": pose["candidate_id"],
        "cluster_size": int(pose.get("cluster_size", 2)),
        "topology": pose.get("topology", "dimer"),
        "molecule_atom_ranges": pose.get("molecule_atom_ranges"),
        "local_stretch_bonds": classify_local_stretch_bonds(
            pose.get("local_stretch_bonds", []), pose.get("interactions") or [pose],
            int(pose["monomer_atom_count"]),
        ),
        "hydrogen_bond_interactions": [
            {key: item[key] for key in ("donor_atom", "donor_hydrogen", "acceptor_atom")}
            for item in (pose.get("interactions") or [pose])
        ],
        "target_geometry": pose["target_geometry"],
        "site_identity": pose["site_identity"],
        "restraints": {
            "force_constant": force_constant,
            "distance": True,
            "angle": True,
            "dihedral": any(item.get("orientation_atoms") for item in (pose.get("interactions") or [pose])),
            "interaction_count": len(pose.get("interactions") or [pose]),
        },
        "population_weight": None,
        "population_model": "not_assigned",
        "acquisition_round": int(pose.get("acquisition_round", 0)),
    }
    try:
        with output_path.open("w") as output_handle:
            result = subprocess.run(
                command, cwd=candidate_dir, env=runtime_environment,
                stdout=output_handle, stderr=subprocess.STDOUT,
                timeout=timeout_seconds,
            )
        output_text = output_path.read_text(errors="replace")
        if result.returncode != 0 or "normal termination of xtb" not in output_text.lower():
            raise RuntimeError(f"xTB did not terminate normally (return code {result.returncode})")
        native_optimized = candidate_dir / "xtbopt.xyz"
        if not native_optimized.is_file():
            raise RuntimeError("xTB completed without xtbopt.xyz")
        shutil.copy2(native_optimized, optimized_path)
        validation = _validate_optimized_pose(pose, optimized_path)
        energy = _parse_xtb_energy(candidate_dir / "xtbout.json", output_text)
        record = {
            **base_record,
            "optimization_status": "completed",
            "sampling_status": "retained" if validation["valid"] else "rejected_geometry",
            "xtb_energy_hartree": energy,
            "relative_xtb_energy_kcal_mol": None,
            "optimized_xyz": str(optimized_path),
            "geometry_validation": validation,
            "environment_features": validation["features"],
            "error": None,
        }
    except Exception as error:
        record = {
            **base_record,
            "optimization_status": "failed",
            "sampling_status": "failed",
            "xtb_energy_hartree": None,
            "relative_xtb_energy_kcal_mol": None,
            "optimized_xyz": None,
            "geometry_validation": None,
            "environment_features": None,
            "error": str(error),
        }
    _json_write(record_path, record)
    return record


def _merge_sampling_report(path: Path, sampling: dict) -> Path:
    retained = json.loads(path.read_text()) if path.is_file() else {}
    retained["schema_version"] = max(2, int(retained.get("schema_version", 0)))
    retained["sampling"] = sampling
    return _json_write(path, retained)


def sample_xtb_dimer_environments(
    smiles: str,
    run_dir: Path,
    *,
    fidelity: str,
    monomer_xyz: Path,
    charge: int = 0,
    multiplicity: int = 1,
    ncores: int = 1,
    seed: int = 42,
    candidate_count: int | None = None,
    force_constant: float = DEFAULT_RESTRAINT_FORCE_CONSTANT,
    executable: str | None = None,
    progress: Callable[[str], None] | None = None,
    include_trimers: bool = False,
    candidate_id_prefix: str = "",
    acquisition_round: int = 0,
    merge_existing: bool = False,
) -> tuple[list[dict], Path]:
    """Generate, optimize, validate, and record a dimer/network ensemble."""
    settings = xtb_sampling_defaults(fidelity)
    requested = settings["candidate_count"] if candidate_count is None else int(candidate_count)
    if requested < 1:
        raise ValueError("This fidelity tier does not request xTB environment sampling")
    trimer_count = int(round(requested * 0.55)) if include_trimers else 0
    dimer_count = requested - trimer_count
    poses = generate_stratified_dimer_poses(
        smiles, candidate_count=dimer_count, seed=seed, monomer_xyz=monomer_xyz,
    )
    if trimer_count:
        poses.extend(generate_stratified_trimer_poses(
            smiles, candidate_count=trimer_count, seed=seed + 1, monomer_xyz=monomer_xyz,
        ))
    if candidate_id_prefix:
        for position, pose in enumerate(poses, start=1):
            pose["candidate_id"] = f"{candidate_id_prefix}-{position:04d}"
    for pose in poses:
        pose["acquisition_round"] = int(acquisition_round)
    executable = executable or find_xtb()
    executable_version = _xtb_version(str(Path(executable).resolve()))
    geometry_dir = Path(run_dir) / "environment-geometries"
    records: list[dict] = []
    for position, pose in enumerate(poses, start=1):
        if progress:
            progress(f"xTB environment {position}/{len(poses)}: {pose['candidate_id']}")
        records.append(run_xtb_candidate(
            pose, geometry_dir / pose["candidate_id"], charge=charge,
            multiplicity=multiplicity, ncores=ncores, force_constant=force_constant,
            executable=executable, executable_version=executable_version,
        ))
        records[-1]["acquisition_round"] = int(acquisition_round)
        _json_write(
            geometry_dir / records[-1]["candidate_id"] / "sampling-record.json",
            records[-1],
        )
    if merge_existing:
        sampling_path = Path(run_dir) / "environment-sampling.json"
        retained_payload = json.loads(sampling_path.read_text()) if sampling_path.is_file() else {}
        existing = (retained_payload.get("sampling") or {}).get("candidates") or []
        by_id = {record["candidate_id"]: record for record in existing}
        for record in records:
            if record["candidate_id"] in by_id and by_id[record["candidate_id"]] != record:
                raise RuntimeError(f"Extension candidate collides with retained record: {record['candidate_id']}")
            by_id[record["candidate_id"]] = record
        records = list(existing) + [
            record for record in records if record["candidate_id"] not in {
                item["candidate_id"] for item in existing
            }
        ]
    energy_minima = {}
    for record in records:
        energy = record.get("xtb_energy_hartree")
        if record.get("sampling_status") == "retained" and isinstance(energy, (int, float)):
            stratum = (record.get("cluster_size", 2), record.get("topology", "dimer"))
            energy_minima[stratum] = min(float(energy), energy_minima.get(stratum, math.inf))
    if energy_minima:
        for record in records:
            energy = record.get("xtb_energy_hartree")
            stratum = (record.get("cluster_size", 2), record.get("topology", "dimer"))
            if isinstance(energy, (int, float)) and stratum in energy_minima:
                relative = (float(energy) - energy_minima[stratum]) * HARTREE_TO_KCAL_MOL
                record["relative_xtb_energy_kcal_mol"] = relative
                features = record.get("environment_features")
                if isinstance(features, dict):
                    features["relative_xtb_energy_kcal_mol"] = relative
                _json_write(
                    geometry_dir / record["candidate_id"] / "sampling-record.json", record,
                )
    summary = {
        "configuration": {
            **settings,
            "candidate_count": requested,
            "seed": seed,
            "restraint_force_constant": force_constant,
            "xtb_executable": str(Path(executable).resolve()),
            "xtb_version": executable_version,
            "monomer_xyz": str(Path(monomer_xyz)),
            "cluster_model": "dimers_plus_trimers" if include_trimers else "dimers_only",
            "dimer_candidate_count": dimer_count,
            "trimer_candidate_count": trimer_count,
            "acquisition_round": int(acquisition_round),
            "extension_merge": bool(merge_existing),
            "cumulative_candidate_count": len(records),
        },
        "candidate_summary": {
            "requested": len(records),
            "generated": len(records),
            "optimization_completed": sum(record["optimization_status"] == "completed" for record in records),
            "retained": sum(record["sampling_status"] == "retained" for record in records),
            "rejected_geometry": sum(record["sampling_status"] == "rejected_geometry" for record in records),
            "failed": sum(record["sampling_status"] == "failed" for record in records),
            "by_topology": {
                topology: sum(record.get("topology", "dimer") == topology and record.get("sampling_status") == "retained" for record in records)
                for topology in sorted({record.get("topology", "dimer") for record in records})
            },
        },
        "population_model": {
            "kind": "not_assigned",
            "warning": settings["population_warning"],
            "vacuum_xtb_energies_used_as_liquid_populations": False,
        },
        "candidates": records,
    }
    manifest = _merge_sampling_report(Path(run_dir) / "environment-sampling.json", summary)
    return records, manifest
