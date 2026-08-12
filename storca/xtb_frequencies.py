"""Resumable unrestrained xTB Hessians for retained environment snapshots."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Callable, Iterable

import numpy as np

from src.orca_runner import find_xtb

from .cluster_ir import classify_local_stretch_bonds, parse_xcontrol_interactions
from .ir_modes import NormalModeSet, local_stretch_mode_assignments
from .mode_character import mode_character_fingerprints


BOHR_PER_ANGSTROM = 1.889726125
XTB_FREQUENCY_SCHEMA_VERSION = 1

_ELEMENTS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu "
    "Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe"
).split()


@dataclass(frozen=True)
class XTBFrequencyResult:
    mode_set: NormalModeSet
    intensities_km_mol: np.ndarray


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _xtb_version(executable: str) -> str:
    result = subprocess.run(
        [executable, "--version"], capture_output=True, text=True, timeout=30,
    )
    match = re.search(
        r"xtb version\s+([^\s]+)", f"{result.stdout}\n{result.stderr}", re.IGNORECASE,
    )
    if result.returncode or not match:
        raise RuntimeError(f"Could not determine xTB version from {executable}")
    return match.group(1)


def parse_xtb_g98(path: Path) -> XTBFrequencyResult:
    """Parse xTB's Gaussian-style coordinates, modes, and IR intensities."""
    lines = Path(path).read_text(errors="replace").splitlines()
    orientation = next(
        (index for index, line in enumerate(lines) if "Standard orientation:" in line), None,
    )
    if orientation is None:
        raise ValueError("xTB g98.out has no standard orientation")
    position = orientation + 1
    while position < len(lines) and "----" not in lines[position]:
        position += 1
    position += 1
    while position < len(lines) and "----" not in lines[position]:
        position += 1
    position += 1
    labels: list[str] = []
    coordinates: list[list[float]] = []
    while position < len(lines) and "----" not in lines[position]:
        fields = lines[position].split()
        if len(fields) >= 6:
            atomic_number = int(fields[1])
            if not 0 < atomic_number < len(_ELEMENTS):
                raise ValueError("xTB g98.out contains an unsupported atomic number")
            labels.append(_ELEMENTS[atomic_number])
            coordinates.append([float(value) * BOHR_PER_ANGSTROM for value in fields[3:6]])
        position += 1
    if not labels:
        raise ValueError("xTB g98.out contains no atoms")

    frequencies: list[float] = []
    intensities: list[float] = []
    vectors: list[np.ndarray] = []
    position = 0
    while position < len(lines):
        if "Frequencies --" not in lines[position]:
            position += 1
            continue
        block_frequencies = [
            float(value.replace("D", "E"))
            for value in lines[position].split("--", 1)[1].split()
        ]
        count = len(block_frequencies)
        intensity_position = next(
            (index for index in range(position + 1, min(position + 10, len(lines)))
             if "IR Inten" in lines[index]), None,
        )
        atom_header = next(
            (index for index in range(position + 1, min(position + 14, len(lines)))
             if lines[index].strip().startswith("Atom AN")), None,
        )
        if intensity_position is None or atom_header is None:
            raise ValueError("Malformed xTB g98.out frequency block")
        block_intensities = [
            float(value.replace("D", "E"))
            for value in lines[intensity_position].split("--", 1)[1].split()
        ]
        if len(block_intensities) != count:
            raise ValueError("xTB g98.out frequency/intensity count differs")
        block_vectors = np.zeros((count, len(labels), 3), dtype=float)
        for atom in range(len(labels)):
            fields = lines[atom_header + 1 + atom].split()
            if len(fields) != 2 + 3 * count or int(fields[0]) != atom + 1:
                raise ValueError("Malformed xTB g98.out normal-coordinate row")
            values = [float(value.replace("D", "E")) for value in fields[2:]]
            for mode in range(count):
                block_vectors[mode, atom] = values[3 * mode:3 * mode + 3]
        frequencies.extend(block_frequencies)
        intensities.extend(block_intensities)
        vectors.extend(block_vectors)
        position = atom_header + len(labels) + 1
    arrays = (np.asarray(frequencies), np.asarray(intensities), np.asarray(vectors))
    if not frequencies or any(not np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("xTB g98.out contains incomplete or non-finite modes")
    expected = 3 * len(labels) - (5 if len(labels) == 2 else 6)
    if len(frequencies) != expected:
        raise ValueError(
            f"xTB g98.out contains {len(frequencies)} vibrational modes; expected {expected}"
        )
    return XTBFrequencyResult(
        NormalModeSet(
            arrays[0], arrays[2], tuple(labels), np.asarray(coordinates, dtype=float),
        ),
        arrays[1],
    )


def _record_from_result(result: XTBFrequencyResult, local_bonds: list[dict]) -> dict:
    mode_set = result.mode_set
    fingerprints = mode_character_fingerprints(mode_set)
    local_modes = local_stretch_mode_assignments(mode_set, local_bonds)
    modes = [
        {
            "mode": index,
            "freq": float(frequency),
            "intensity": float(result.intensities_km_mol[index]),
        }
        for index, frequency in enumerate(mode_set.frequencies_cm_1)
    ]
    for assignment in local_modes:
        mode = modes[int(assignment["mode"])]
        mode.update({key: value for key, value in assignment.items() if key != "mode"})
    imaginary = [float(value) for value in mode_set.frequencies_cm_1 if value < 0]
    return {
        "mode_count": len(modes),
        "modes": modes,
        "local_stretch_modes": local_modes,
        "mode_character_fingerprints": fingerprints,
        "snapshot_reliability": {
            "imaginary_mode_count": len(imaginary),
            "lowest_frequency_cm-1": min(map(float, mode_set.frequencies_cm_1)),
            "material_imaginary_mode_count": sum(value < -100.0 for value in imaginary),
            "status": "usable_snapshot" if not any(value < -250.0 for value in imaginary)
                      else "large_imaginary_curvature",
        },
    }


def run_xtb_snapshot_frequency(
    sampling_record: dict,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    ncores: int = 1,
    executable: str | None = None,
    executable_version: str | None = None,
    timeout_seconds: float = 600.0,
    retry_failed: bool = False,
) -> dict:
    """Run or reuse an unrestrained numerical Hessian at one retained snapshot."""
    if multiplicity < 1 or ncores < 1:
        raise ValueError("Multiplicity and xTB core count must be positive")
    source = Path(str(sampling_record.get("optimized_xyz") or ""))
    if sampling_record.get("sampling_status") != "retained" or not source.is_file():
        raise ValueError("xTB snapshot frequencies require a retained optimized candidate")
    job_dir = source.parent / "xtb-frequency"
    job_dir.mkdir(parents=True, exist_ok=True)
    input_path = job_dir / "input.xyz"
    source_text = source.read_text()
    if input_path.is_file() and input_path.read_text() != source_text:
        raise RuntimeError(f"Refusing to overwrite a different retained xTB Hessian input: {input_path}")
    input_path.write_text(source_text)
    executable = str(Path(executable or find_xtb()).resolve())
    version = executable_version or _xtb_version(executable)
    contract = {
        "schema_version": XTB_FREQUENCY_SCHEMA_VERSION,
        "candidate_id": sampling_record["candidate_id"],
        "input_xyz_sha256": _sha256_file(input_path),
        "xtb_executable": executable,
        "xtb_version": version,
        "method": "GFN2-xTB",
        "calculation": "unrestrained_snapshot_hessian",
        "charge": int(charge),
        "unpaired_electrons": int(multiplicity - 1),
        "ncores": int(ncores),
    }
    contract_path = job_dir / "calculation-contract.json"
    record_path = job_dir / "frequency-record.json"
    retained_contract = json.loads(contract_path.read_text()) if contract_path.is_file() else None
    retained_record = json.loads(record_path.read_text()) if record_path.is_file() else None
    local_bonds = sampling_record.get("local_stretch_bonds") or []
    legacy_local_bonds = bool(
        local_bonds and any("spectral_band_class" not in bond for bond in local_bonds)
    )
    if legacy_local_bonds:
        ranges = sampling_record.get("molecule_atom_ranges") or []
        if ranges:
            monomer_atom_count = int(ranges[0][1]) - int(ranges[0][0])
        else:
            atom_count = int(source.read_text().splitlines()[0])
            monomer_atom_count = atom_count // max(1, int(sampling_record.get("cluster_size") or 1))
        interactions = sampling_record.get("hydrogen_bond_interactions") or []
        if not interactions and (source.parent / "xcontrol.inp").is_file():
            interactions = parse_xcontrol_interactions(source.parent / "xcontrol.inp")
        local_bonds = classify_local_stretch_bonds(local_bonds, interactions, monomer_atom_count)
    if retained_contract is not None and retained_contract != contract:
        raise RuntimeError(f"Retained xTB Hessian has a different calculation contract: {job_dir}")
    if retained_contract == contract and retained_record:
        if retained_record.get("frequency_status") == "completed":
            changed = False
            if (
                (legacy_local_bonds or any(
                    "spectral_band_class" not in item
                    for item in retained_record.get("local_stretch_modes") or []
                ))
                and (job_dir / "g98.out").is_file()
            ):
                retained_record.update(_record_from_result(
                    parse_xtb_g98(job_dir / "g98.out"), local_bonds,
                ))
                changed = True
            additions = {
                "cluster_size": int(sampling_record.get("cluster_size") or 1),
                "topology": sampling_record.get("topology"),
                "environment_features": sampling_record.get("environment_features"),
                "acquisition_round": int(sampling_record.get("acquisition_round", 0)),
            }
            for key, value in additions.items():
                if retained_record.get(key) != value:
                    retained_record[key] = value
                    changed = True
            if changed:
                record_path.write_text(
                    json.dumps(retained_record, indent=2, sort_keys=True) + "\n"
                )
            return retained_record
        if retained_record.get("frequency_status") == "failed" and not retry_failed:
            return retained_record
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    output_path = job_dir / "xtb.out"
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": str(ncores), "MKL_NUM_THREADS": str(ncores),
        "OMP_STACKSIZE": environment.get("OMP_STACKSIZE", "1G"),
    })
    command = [
        executable, input_path.name, "--gfn", "2", "--hess", "--chrg", str(charge),
        "--uhf", str(multiplicity - 1), "--norestart",
    ]
    base = {
        "schema_version": XTB_FREQUENCY_SCHEMA_VERSION,
        "candidate_id": sampling_record["candidate_id"],
        "calculation": "unrestrained_snapshot_hessian",
        "restraints_in_hessian": False,
        "population_weight": sampling_record.get("population_weight"),
        "cluster_size": int(sampling_record.get("cluster_size") or 1),
        "topology": sampling_record.get("topology"),
        "environment_features": sampling_record.get("environment_features"),
        "acquisition_round": int(sampling_record.get("acquisition_round", 0)),
    }
    try:
        with output_path.open("w") as handle:
            completed = subprocess.run(
                command, cwd=job_dir, env=environment, stdout=handle,
                stderr=subprocess.STDOUT, timeout=timeout_seconds,
            )
        output_text = output_path.read_text(errors="replace")
        g98_path = job_dir / "g98.out"
        if completed.returncode or "normal termination of xtb" not in output_text.lower():
            raise RuntimeError(f"xTB Hessian did not terminate normally (return code {completed.returncode})")
        if not g98_path.is_file():
            raise RuntimeError("xTB Hessian completed without g98.out")
        parsed = _record_from_result(
            parse_xtb_g98(g98_path), local_bonds,
        )
        record = {**base, "frequency_status": "completed", **parsed, "error": None}
    except Exception as error:
        record = {**base, "frequency_status": "failed", "error": str(error)}
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def sample_xtb_snapshot_frequencies(
    records: Iterable[dict], run_dir: Path, *, charge: int = 0, multiplicity: int = 1,
    ncores: int = 1, executable: str | None = None,
    progress: Callable[[str], None] | None = None,
) -> tuple[list[dict], Path]:
    """Evaluate retained candidates and write the ensemble frequency manifest."""
    record_list = list(records)
    retained = [record for record in record_list if record.get("sampling_status") == "retained"]
    executable = str(Path(executable or find_xtb()).resolve())
    version = _xtb_version(executable)
    results = []
    for position, record in enumerate(retained, start=1):
        if progress:
            progress(f"xTB snapshot frequency {position}/{len(retained)}: {record['candidate_id']}")
        result = run_xtb_snapshot_frequency(
            record, charge=charge, multiplicity=multiplicity, ncores=ncores,
            executable=executable, executable_version=version,
        )
        results.append(result)
        local = [
            float(mode["freq"]) for mode in result.get("modes", [])
            if mode.get("spectral_band_class") in {"hydrogen_bonded_oh", "non_donating_oh"}
        ]
        if local and isinstance(record.get("environment_features"), dict):
            record["environment_features"]["estimated_local_frequency_cm-1"] = float(np.mean(local))
            result["environment_features"] = record["environment_features"]
            frequency_record_path = (
                Path(str(record["optimized_xyz"])).parent
                / "xtb-frequency" / "frequency-record.json"
            )
            frequency_record_path.write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n"
            )
            sampling_path = Path(str(record["optimized_xyz"])).parent / "sampling-record.json"
            if sampling_path.is_file():
                sampling_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    manifest = Path(run_dir) / "xtb-snapshot-frequencies.json"
    payload = {
        "schema_version": XTB_FREQUENCY_SCHEMA_VERSION,
        "kind": "unrestrained_xtb_snapshot_hessian_ensemble",
        "method": "GFN2-xTB",
        "restraints_in_hessian": False,
        "requested": len(retained),
        "completed": sum(item.get("frequency_status") == "completed" for item in results),
        "failed": sum(item.get("frequency_status") == "failed" for item in results),
        "candidates": results,
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    sampling_manifest = Path(run_dir) / "environment-sampling.json"
    if sampling_manifest.is_file():
        retained_payload = json.loads(sampling_manifest.read_text())
        sampling = retained_payload.get("sampling") or {}
        by_id = {record["candidate_id"]: record for record in record_list}
        sampling["candidates"] = [
            by_id.get(record.get("candidate_id"), record)
            for record in sampling.get("candidates") or []
        ]
        retained_payload["sampling"] = sampling
        retained_payload["xtb_snapshot_frequencies"] = {
            "manifest": str(manifest), "completed": payload["completed"],
            "failed": payload["failed"], "restraints_in_hessian": False,
        }
        sampling_manifest.write_text(
            json.dumps(retained_payload, indent=2, sort_keys=True) + "\n"
        )
    return results, manifest
