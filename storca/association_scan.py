"""Resumable, endpoint-seeded ORCA radical-association scans."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.inputgen import create_orca_constrained_opt_input
from src.orca_runner import run_orca
from src.parser import parse_orca_energy

from .route_geometry import read_xyz, write_xyz


def _normal(path: Path) -> bool:
    return Path(path).is_file() and "ORCA TERMINATED NORMALLY" in Path(path).read_text(errors="replace")


def _endpoint_job(repair: dict, smiles: str) -> dict:
    return next(job for job in repair["endpoint_jobs"] if job["endpoint"]["smiles"] == smiles)


def _seed(oxygen_xyz: Path, radical_xyz: Path, output: Path, *, orientation: int) -> tuple[Path, tuple[int, int]]:
    """Combine optimized endpoints in a reproducible, H-avoiding orientation."""
    oxygen_elements, oxygen, _ = read_xyz(oxygen_xyz)
    radical_elements, radical, _ = read_xyz(radical_xyz)
    # RMG's [NH]NO retains its radical centre as atom 0; this only accepts the
    # already resolved NNO collision-repair route and refuses silent guessing.
    if oxygen_elements != ["O", "O"] or not radical_elements or radical_elements[0] != "N":
        raise ValueError("Endpoint atom ordering is not supported for the association scan")
    oxygen -= oxygen[0]
    offset = np.array([5.5, 0.9 if orientation == 1 else -0.9, 0.8 if orientation == 3 else 0.0])
    oxygen += radical[0] + offset
    write_xyz(output, oxygen_elements + radical_elements, np.vstack((oxygen, radical)),
              "validated-endpoint O2/radical association seed")
    return output, (0, len(oxygen_elements))


def _distance(xyz: Path, pair: tuple[int, int]) -> float:
    _, coordinates, _ = read_xyz(xyz)
    return float(np.linalg.norm(coordinates[pair[0]] - coordinates[pair[1]]))


def _validate_geometry(xyz: Path, pair: tuple[int, int], target: float) -> dict:
    elements, coordinates, _ = read_xyz(xyz)
    actual = _distance(xyz, pair)
    # The constrained bond must be retained; an H---O contact shorter than
    # 0.75 A is a clear malformed/rearranged geometry, not curve evidence.
    contacts = []
    for i, element in enumerate(elements):
        if element != "H":
            continue
        for j, other in enumerate(elements):
            if other != "O" or i == j:
                continue
            value = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            if value < 0.75:
                contacts.append({"hydrogen_index": i, "oxygen_index": j, "distance_angstrom": value})
    valid = abs(actual - target) <= 0.08 and not contacts
    return {"valid": valid, "constrained_distance_angstrom": actual,
            "target_distance_angstrom": target, "unintended_short_h_o_contacts": contacts,
            "status": "validated" if valid else "coordinate_inadequate_or_competing_rearrangement"}


def _point(folder: Path, seed: Path, pair: tuple[int, int], distance: float, *, ncores: int,
           method_keywords: list[str] | None, max_iterations: int, retry: bool) -> dict:
    folder.mkdir(parents=True, exist_ok=True)
    local_seed = folder / "seed.xyz"
    if not local_seed.is_file():
        local_seed.write_text(Path(seed).read_text())
    label = "opt-retry" if retry else "opt"
    inp = create_orca_constrained_opt_input(local_seed, bond_atom_indices=pair, distance_angstrom=distance,
        charge=0, multiplicity=2, label=label, ncores=ncores, max_iterations=max_iterations,
        fresh_hessian=retry or distance >= 5.49, method_keywords=method_keywords)
    out = inp.with_suffix(".out")
    try:
        print(f"[STORCA] ORCA association point {distance:.2f} Å ({'retry' if retry else 'primary'}) started", flush=True)
        if not _normal(out):
            artifacts = run_orca(inp)
            out = Path(artifacts["out"])
        xyz = inp.with_suffix(".xyz")
        if not xyz.is_file():
            raise RuntimeError("ORCA did not retain an optimized point geometry")
        validation = _validate_geometry(xyz, pair, distance)
        status = "completed" if validation["valid"] else validation["status"]
        print(f"[STORCA] ORCA association point {distance:.2f} Å {status}", flush=True)
        return {"status": status, "input": str(inp),
                "output": str(out), "optimized_xyz": str(xyz), "energy_hartree": parse_orca_energy(out),
                "validation": validation, "retry": retry}
    except Exception as error:
        print(f"[STORCA] ORCA association point {distance:.2f} Å failed: {error}", flush=True)
        return {"status": "failed", "input": str(inp), "output": str(out), "failure_reason": str(error), "retry": retry}


def run_resumable_association_scans(run_dir: Path, *, execute: bool = False, ncores: int = 1,
                                    points: int = 16, start_distance: float = 5.5,
                                    end_distance: float = 2.0, method_keywords: list[str] | None = None) -> dict:
    """Prepare or run three independent doublet O2/radical approach curves."""
    if points < 5 or start_distance <= end_distance:
        raise ValueError("Association scan needs at least five descending points")
    run_dir = Path(run_dir)
    repair_path = run_dir / "collision-repair" / "collision-repair.json"
    repair = json.loads(repair_path.read_text())["repairs"][0]
    oxygen = _endpoint_job(repair, "[O][O]")
    radical = next(job for job in repair["endpoint_jobs"] if job["endpoint"].get("multiplicity") == 2 and job["name"].startswith("reactant"))
    if execute and not all(job.get("status") == "completed" and (job.get("frequency_check") or {}).get("IsMinimum") for job in repair["endpoint_jobs"]):
        raise ValueError("All ORCA endpoint jobs must be completed local minima before executing association scans")
    distances = [float(value) for value in np.linspace(start_distance, end_distance, points)]
    orientations = []
    for orientation in (1, 2, 3):
        print(f"[STORCA] Preparing association orientation {orientation}/3 ({points} points)", flush=True)
        base = run_dir / "collision-repair" / "collision-route-1" / "association-scan" / f"orientation-{orientation}"
        initial, pair = _seed(Path(oxygen.get("optimized_xyz") or oxygen["input_xyz"]),
                              Path(radical.get("optimized_xyz") or radical["input_xyz"]), base / "initial.xyz", orientation=orientation)
        records, previous = [], initial
        for number, distance in enumerate(distances):
            folder = base / f"point-{number:02d}-{distance:.2f}A"
            if execute:
                record = _point(folder, previous, pair, distance, ncores=ncores, method_keywords=method_keywords,
                                max_iterations=300, retry=False)
                if record["status"] == "failed":
                    record["retry_record"] = _point(folder, previous, pair, distance, ncores=ncores,
                        method_keywords=method_keywords, max_iterations=500, retry=True)
                    if record["retry_record"]["status"] == "completed":
                        record = record["retry_record"]
                if record["status"] == "completed":
                    previous = Path(record["optimized_xyz"])
            else:
                seed = folder / "seed.xyz"; folder.mkdir(parents=True, exist_ok=True); seed.write_text(Path(previous).read_text())
                inp = create_orca_constrained_opt_input(seed, bond_atom_indices=pair, distance_angstrom=distance,
                    charge=0, multiplicity=2, label="opt", ncores=ncores, max_iterations=300,
                    fresh_hessian=number == 0, method_keywords=method_keywords)
                record = {"status": "prepared", "input": str(inp), "seed_xyz": str(seed), "retry": False}
            record.update(point_index=number, target_distance_angstrom=distance)
            records.append(record)
        valid = [point for point in records if point["status"] == "completed"]
        orientations.append({"orientation": orientation, "coordinate_atom_indices": list(pair), "points": records,
                             "valid_point_count": len(valid), "coverage_angstrom": (max((x["target_distance_angstrom"] for x in valid), default=0) - min((x["target_distance_angstrom"] for x in valid), default=0))})
    result = {"schema_version": 2, "kind": "resumable_orca_association_scans", "status": "computed" if execute else "prepared",
              "source_repair": str(repair_path), "distances_angstrom": distances, "orientations": orientations,
              "limitations": ["Each curve is a constrained doublet ground-state scan, not a rate or a transition-state verification."]}
    path = run_dir / "collision-repair" / "association-scans.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {**result, "path": path}
