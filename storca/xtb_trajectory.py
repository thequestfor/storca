"""Restrained finite-temperature xTB trajectories and decorrelated snapshots."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Callable

import numpy as np

from src.orca_runner import find_xtb

from .xtb_sampling import _angle_degrees, _read_xyz, _xtb_version, _xyz_text


TRAJECTORY_SCHEMA_VERSION = 1
TRAJECTORY_POPULATION_MODEL = "restrained_xtb_trajectory_occupancy"
TRAJECTORY_POPULATION_WARNING = (
    "Weights are finite-sample occupancies within the restrained seeded GFN2-xTB "
    "trajectory protocol; they are not unbiased bulk-liquid equilibrium populations."
)


def xtb_trajectory_defaults(fidelity: str) -> dict:
    """Return the bounded trajectory contract for an experimental fidelity tier."""
    fidelity = fidelity.strip().lower()
    settings = {
        "fast": (0, 0.0, 0),
        "auto": (4, 2.0, 40),
        "balanced": (6, 4.0, 75),
        "accurate": (8, 5.0, 100),
    }
    if fidelity not in settings:
        raise ValueError("Unsupported xTB trajectory fidelity")
    seed_count, production_time_ps, target_snapshots = settings[fidelity]
    return {
        "schema_version": TRAJECTORY_SCHEMA_VERSION,
        "fidelity": fidelity,
        "method": "GFN2-xTB",
        "temperature_K": 298.15,
        "trajectory_seed_count": seed_count,
        "production_time_ps": production_time_ps,
        "time_step_fs": 0.5,
        "dump_interval_fs": 10.0,
        "hydrogen_mass_amu": 4,
        "equilibration_fraction": 0.25,
        "target_snapshot_count": target_snapshots,
        "decorrelation_method": "initial_positive_sequence_statistical_inefficiency",
        "population_model": TRAJECTORY_POPULATION_MODEL,
        "population_warning": TRAJECTORY_POPULATION_WARNING,
    }


def read_xyz_trajectory(path: Path) -> list[dict]:
    """Parse an xTB multi-XYZ trajectory without assuming a comment format."""
    lines = Path(path).read_text(errors="replace").splitlines()
    frames: list[dict] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError as error:
            raise ValueError(f"Invalid XYZ trajectory atom count at line {cursor + 1}: {path}") from error
        if atom_count < 1 or cursor + atom_count + 1 >= len(lines):
            raise ValueError(f"Incomplete XYZ trajectory frame at line {cursor + 1}: {path}")
        comment = lines[cursor + 1]
        rows = [line.split() for line in lines[cursor + 2:cursor + 2 + atom_count]]
        if len(rows) != atom_count or any(len(row) < 4 for row in rows):
            raise ValueError(f"Incomplete XYZ trajectory coordinates at line {cursor + 1}: {path}")
        try:
            coordinates = np.asarray([[float(value) for value in row[1:4]] for row in rows])
        except ValueError as error:
            raise ValueError(f"Non-numeric XYZ trajectory coordinates: {path}") from error
        if not np.all(np.isfinite(coordinates)):
            raise ValueError(f"Non-finite XYZ trajectory coordinates: {path}")
        frames.append({
            "frame_index": len(frames),
            "comment": comment,
            "symbols": [row[0] for row in rows],
            "coordinates": coordinates,
        })
        cursor += atom_count + 2
    if not frames:
        raise ValueError(f"xTB trajectory contains no frames: {path}")
    return frames


def statistical_inefficiency(values: list[float] | np.ndarray) -> float:
    """Estimate g = 1 + 2 sum rho(t) using the initial positive sequence."""
    series = np.asarray(values, dtype=float)
    if series.ndim != 1 or len(series) < 2 or not np.all(np.isfinite(series)):
        return 1.0
    centered = series - series.mean()
    variance = float(np.dot(centered, centered) / len(centered))
    if variance <= 1e-20:
        return 1.0
    correlation_sum = 0.0
    # Stop well before the noisy tail and at the first non-positive estimate.
    for lag in range(1, min(len(series) // 2 + 1, 1000)):
        correlation = float(np.dot(centered[:-lag], centered[lag:]) / ((len(series) - lag) * variance))
        if not math.isfinite(correlation) or correlation <= 0.0:
            break
        correlation_sum += correlation
    return max(1.0, min(float(len(series)), 1.0 + 2.0 * correlation_sum))


def _interaction_observable(frame: dict, interactions: list[dict]) -> float:
    coordinates = frame["coordinates"]
    if interactions:
        distances = [
            float(np.linalg.norm(
                coordinates[int(item["donor_hydrogen"])] - coordinates[int(item["acceptor_atom"])]
            ))
            for item in interactions
        ]
        return float(np.mean(distances))
    centered = coordinates - coordinates.mean(axis=0)
    return float(np.sqrt(np.mean(np.square(centered))))


def decorrelated_frame_indices(
    frames: list[dict], *, interactions: list[dict] | None = None,
    equilibration_fraction: float = 0.25, maximum_snapshots: int | None = None,
) -> tuple[list[int], dict]:
    """Discard equilibration and retain frames at no less than one estimated g."""
    if not frames or not 0.0 <= equilibration_fraction < 1.0:
        raise ValueError("Trajectory frames and a valid equilibration fraction are required")
    first = min(len(frames) - 1, int(math.ceil(len(frames) * equilibration_fraction)))
    production = frames[first:]
    values = [_interaction_observable(frame, interactions or []) for frame in production]
    inefficiency = statistical_inefficiency(values)
    stride = max(1, int(math.ceil(inefficiency)))
    indices = list(range(first, len(frames), stride))
    if maximum_snapshots is not None and maximum_snapshots > 0 and len(indices) > maximum_snapshots:
        positions = np.linspace(0, len(indices) - 1, maximum_snapshots)
        indices = [indices[int(round(position))] for position in positions]
        indices = list(dict.fromkeys(indices))
    return indices, {
        "raw_frame_count": len(frames),
        "equilibration_frames_discarded": first,
        "production_frame_count": len(production),
        "statistical_inefficiency_frames": inefficiency,
        "decorrelation_stride_frames": stride,
        "retained_frame_count": len(indices),
        "observable": (
            "mean_restrained_hydrogen_acceptor_distance_angstrom"
            if interactions else "root_mean_square_radius_angstrom"
        ),
    }


def _md_input(restraint_text: str, settings: dict) -> str:
    return restraint_text.rstrip() + "\n\n" + "\n".join([
        "$md",
        f"  temp={float(settings['temperature_K']):.8f}",
        f"  time={float(settings['production_time_ps']):.8f}",
        f"  dump={float(settings['dump_interval_fs']):.8f}",
        f"  step={float(settings['time_step_fs']):.8f}",
        "  velo=false",
        "  nvt=true",
        f"  hmass={int(settings['hydrogen_mass_amu'])}",
        "  shake=0",
        "  sccacc=1.0",
        "$end",
        "",
    ])


def _frame_features(seed: dict, symbols: list[str], coordinates: np.ndarray) -> tuple[dict, list[str]]:
    interactions = seed.get("hydrogen_bond_interactions") or []
    measurements = []
    failures: list[str] = []
    for item in interactions:
        donor = int(item["donor_atom"])
        hydrogen = int(item["donor_hydrogen"])
        acceptor = int(item["acceptor_atom"])
        distance = float(np.linalg.norm(coordinates[hydrogen] - coordinates[acceptor]))
        angle = _angle_degrees(coordinates[donor], coordinates[hydrogen], coordinates[acceptor])
        measurements.append((distance, angle, float(np.linalg.norm(coordinates[donor] - coordinates[acceptor]))))
    if measurements and any(not 1.1 <= item[0] <= 3.5 for item in measurements):
        failures.append("restrained_contact_dissociated_or_covalent")
    ranges = seed.get("molecule_atom_ranges") or []
    minimum_inter_distance = math.inf
    for left in range(len(ranges)):
        left_start, left_end = map(int, ranges[left])
        for right in range(left + 1, len(ranges)):
            right_start, right_end = map(int, ranges[right])
            distances = np.linalg.norm(
                coordinates[left_start:left_end, None, :] - coordinates[None, right_start:right_end, :], axis=2,
            )
            minimum_inter_distance = min(minimum_inter_distance, float(distances.min()))
    if minimum_inter_distance < 0.65:
        failures.append("severe_intermolecular_collision")
    seed_symbols, seed_coordinates = _read_xyz(Path(seed["optimized_xyz"]))
    if seed_symbols != symbols or seed_coordinates.shape != coordinates.shape:
        raise ValueError("Trajectory frame changed atom identities or ordering")
    indices = [index for index, symbol in enumerate(symbols) if symbol != "H"] or list(range(len(symbols)))
    reference = seed_coordinates[indices] - seed_coordinates[indices].mean(axis=0)
    candidate = coordinates[indices] - coordinates[indices].mean(axis=0)
    u_matrix, _, v_matrix = np.linalg.svd(candidate.T @ reference)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u_matrix @ v_matrix))
    rmsd = float(np.sqrt(np.mean(np.square(candidate @ (u_matrix @ correction @ v_matrix) - reference))))
    base = dict(seed.get("environment_features") or {})
    if measurements:
        base.update({
            "h_bond_distance_angstrom": float(np.mean([item[0] for item in measurements])),
            "minimum_h_bond_distance_angstrom": min(item[0] for item in measurements),
            "maximum_h_bond_distance_angstrom": max(item[0] for item in measurements),
            "h_bond_distance_span_angstrom": max(item[0] for item in measurements) - min(item[0] for item in measurements),
            "donor_h_acceptor_angle_degrees": float(np.mean([item[1] for item in measurements])),
            "donor_acceptor_distance_angstrom": float(np.mean([item[2] for item in measurements])),
        })
    base.update({
        "geometry_cluster_id": None,
        "heavy_atom_rmsd_angstrom": rmsd,
        "heavy_atom_rmsd_basis": "trajectory_seed_to_snapshot_after_global_alignment",
    })
    return base, failures


def _sha256(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def run_restrained_xtb_trajectory(
    seed: dict, trajectory_dir: Path, *, settings: dict, charge: int = 0,
    multiplicity: int = 1, ncores: int = 1, executable: str | None = None,
    executable_version: str | None = None, maximum_snapshots: int | None = None,
    timeout_seconds: float = 1800.0,
) -> tuple[list[dict], dict]:
    """Run or reuse one restrained NVT trajectory and extract decorrelated frames."""
    if multiplicity < 1 or ncores < 1:
        raise ValueError("Multiplicity and xTB core count must be positive")
    trajectory_dir = Path(trajectory_dir)
    trajectory_dir.mkdir(parents=True, exist_ok=True)
    seed_text = Path(seed["optimized_xyz"]).read_text()
    source_restraint = Path(seed["optimized_xyz"]).parent / "xcontrol.inp"
    if not source_restraint.is_file():
        raise RuntimeError(f"Trajectory seed has no retained restraint input: {source_restraint}")
    md_text = _md_input(source_restraint.read_text(), settings)
    input_path = trajectory_dir / "initial.xyz"
    control_path = trajectory_dir / "md.inp"
    for path, content in ((input_path, seed_text), (control_path, md_text)):
        if path.is_file() and path.read_text() != content:
            raise RuntimeError(f"Refusing to overwrite a different trajectory input: {path}")
        path.write_text(content)
    executable = str(Path(executable or find_xtb()).resolve())
    executable_version = executable_version or _xtb_version(executable)
    contract = {
        "schema_version": TRAJECTORY_SCHEMA_VERSION,
        "seed_candidate_id": seed["candidate_id"],
        "seed_xyz_sha256": _sha256(seed_text),
        "md_input_sha256": _sha256(md_text),
        "xtb_executable": executable,
        "xtb_version": executable_version,
        "charge": int(charge),
        "unpaired_electrons": int(multiplicity - 1),
        "ncores": int(ncores),
        "settings": settings,
    }
    contract_path = trajectory_dir / "calculation-contract.json"
    report_path = trajectory_dir / "trajectory-record.json"
    retained_contract = json.loads(contract_path.read_text()) if contract_path.is_file() else None
    if retained_contract is not None and retained_contract != contract:
        raise RuntimeError(f"Retained xTB trajectory has a different calculation contract: {trajectory_dir}")
    if retained_contract == contract and report_path.is_file():
        report = json.loads(report_path.read_text())
        snapshots = report.get("snapshots") or []
        if report.get("trajectory_status") == "completed" and all(Path(item["optimized_xyz"]).is_file() for item in snapshots):
            return snapshots, report
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    runtime_environment = os.environ.copy()
    runtime_environment.update({
        "OMP_NUM_THREADS": str(ncores), "MKL_NUM_THREADS": str(ncores),
        "OMP_STACKSIZE": runtime_environment.get("OMP_STACKSIZE", "1G"),
    })
    output_path = trajectory_dir / "xtb.out"
    command = [
        executable, input_path.name, "--gfn", "2", "--md", "--input", control_path.name,
        "--chrg", str(charge), "--uhf", str(multiplicity - 1), "--norestart",
    ]
    try:
        with output_path.open("w") as output_handle:
            result = subprocess.run(command, cwd=trajectory_dir, env=runtime_environment,
                                    stdout=output_handle, stderr=subprocess.STDOUT,
                                    timeout=timeout_seconds)
        output_text = output_path.read_text(errors="replace")
        if result.returncode != 0 or not (trajectory_dir / "xtbmdok").is_file():
            raise RuntimeError(f"xTB trajectory did not complete normally (return code {result.returncode})")
        if "md is unstable" in output_text.lower() or "emergency exit" in output_text.lower():
            raise RuntimeError("xTB reported an unstable MD emergency exit")
        trajectory_path = trajectory_dir / "xtb.trj"
        frames = read_xyz_trajectory(trajectory_path)
        expected_frames = max(1, int(round(
            float(settings["production_time_ps"]) * 1000.0
            / float(settings["dump_interval_fs"])
        )))
        minimum_frames = max(1, int(math.floor(0.95 * expected_frames)))
        if len(frames) < minimum_frames:
            raise RuntimeError(
                f"xTB trajectory retained only {len(frames)}/{expected_frames} expected frames"
            )
        indices, decorrelation = decorrelated_frame_indices(
            frames, interactions=seed.get("hydrogen_bond_interactions") or [],
            equilibration_fraction=float(settings["equilibration_fraction"]),
            maximum_snapshots=maximum_snapshots,
        )
        snapshot_dir = trajectory_dir / "snapshots"
        snapshot_dir.mkdir(exist_ok=True)
        snapshots = []
        for frame_index in indices:
            frame = frames[frame_index]
            features, failures = _frame_features(seed, frame["symbols"], frame["coordinates"])
            snapshot_id = f"traj-{seed['candidate_id']}-f{frame_index:06d}"
            frame_dir = snapshot_dir / f"frame-{frame_index:06d}"
            frame_dir.mkdir(exist_ok=True)
            snapshot_path = frame_dir / "optimized.xyz"
            snapshot_path.write_text(_xyz_text(frame["symbols"], frame["coordinates"], snapshot_id))
            snapshots.append({
                **{key: seed.get(key) for key in (
                    "cluster_size", "topology", "molecule_atom_ranges", "local_stretch_bonds",
                    "hydrogen_bond_interactions", "site_identity", "target_geometry",
                )},
                "schema_version": TRAJECTORY_SCHEMA_VERSION,
                "candidate_id": snapshot_id,
                "seed_candidate_id": seed["candidate_id"],
                "trajectory_frame_index": frame_index,
                "trajectory_time_fs": frame_index * float(settings["dump_interval_fs"]),
                "sampling_status": "retained" if not failures else "rejected_geometry",
                "optimization_status": "trajectory_snapshot",
                "optimized_xyz": str(snapshot_path),
                "environment_features": features,
                "geometry_validation": {"valid": not failures, "failures": failures},
                "xtb_energy_hartree": None,
                "relative_xtb_energy_kcal_mol": None,
                "population_weight": None,
                "population_model": TRAJECTORY_POPULATION_MODEL,
                "population_warning": TRAJECTORY_POPULATION_WARNING,
                "trajectory_decorrelation": decorrelation,
                "error": None,
            })
        report = {
            "schema_version": TRAJECTORY_SCHEMA_VERSION,
            "trajectory_status": "completed",
            "seed_candidate_id": seed["candidate_id"],
            "trajectory_xyz": str(trajectory_path),
            "decorrelation": decorrelation,
            "snapshots": snapshots,
            "error": None,
        }
    except Exception as error:
        report = {
            "schema_version": TRAJECTORY_SCHEMA_VERSION,
            "trajectory_status": "failed",
            "seed_candidate_id": seed["candidate_id"],
            "snapshots": [],
            "error": str(error),
        }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return list(report["snapshots"]), report


def _select_trajectory_seeds(records: list[dict], count: int) -> list[dict]:
    eligible = sorted((item for item in records if item.get("sampling_status") == "retained"),
                      key=lambda item: item["candidate_id"])
    if len(eligible) <= count:
        return eligible
    from .environment_selection import deterministic_k_medoids, pairwise_environment_distances
    matrix, _ = pairwise_environment_distances(eligible)
    medoids, _ = deterministic_k_medoids(eligible, matrix, count)
    return [eligible[index] for index in medoids]


def sample_restrained_xtb_trajectories(
    records: list[dict], run_dir: Path, *, fidelity: str, charge: int = 0,
    multiplicity: int = 1, ncores: int = 1, temperature_K: float = 298.15,
    executable: str | None = None, progress: Callable[[str], None] | None = None,
) -> tuple[list[dict], Path]:
    """Run diverse restrained trajectories and assign decorrelated frame occupancies."""
    settings = {**xtb_trajectory_defaults(fidelity), "temperature_K": float(temperature_K)}
    if not math.isfinite(float(temperature_K)) or float(temperature_K) <= 0:
        raise ValueError("Trajectory temperature must be finite and positive")
    if int(settings["trajectory_seed_count"]) < 1 or int(settings["target_snapshot_count"]) < 1:
        raise ValueError("This fidelity tier does not request restrained xTB trajectories")
    seeds = _select_trajectory_seeds(records, int(settings["trajectory_seed_count"]))
    if not seeds:
        raise RuntimeError("No retained xTB environments are available as trajectory seeds")
    executable = str(Path(executable or find_xtb()).resolve())
    executable_version = _xtb_version(executable)
    target_snapshots = int(settings["target_snapshot_count"])
    base_limit, extra_limits = divmod(target_snapshots, len(seeds))
    snapshot_limits = [max(1, base_limit + (position < extra_limits)) for position in range(len(seeds))]
    all_snapshots: list[dict] = []
    trajectory_reports = []
    root = Path(run_dir) / "environment-trajectories"
    for position, seed in enumerate(seeds, start=1):
        if progress:
            progress(
                f"Restrained {float(settings['temperature_K']):g} K xTB trajectory "
                f"{position}/{len(seeds)}: {seed['candidate_id']}"
            )
        snapshots, report = run_restrained_xtb_trajectory(
            seed, root / seed["candidate_id"], settings=settings, charge=charge,
            multiplicity=multiplicity, ncores=ncores, executable=executable,
            executable_version=executable_version,
            maximum_snapshots=snapshot_limits[position - 1],
        )
        trajectory_reports.append(report)
        all_snapshots.extend(item for item in snapshots if item.get("sampling_status") == "retained")
    if not all_snapshots:
        raise RuntimeError("No valid decorrelated xTB trajectory snapshots were retained")
    completed_trajectories = sum(
        report.get("trajectory_status") == "completed" for report in trajectory_reports
    )
    if completed_trajectories != len(seeds):
        raise RuntimeError(
            f"Only {completed_trajectories}/{len(seeds)} restrained xTB trajectories completed; "
            "refusing biased partial-seed occupancy weights"
        )
    snapshots_by_seed = {
        seed["candidate_id"]: [
            snapshot for snapshot in all_snapshots
            if snapshot["seed_candidate_id"] == seed["candidate_id"]
        ]
        for seed in seeds
    }
    populated_seeds = [items for items in snapshots_by_seed.values() if items]
    for seed_snapshots in populated_seeds:
        weight = 1.0 / len(populated_seeds) / len(seed_snapshots)
        for snapshot in seed_snapshots:
            snapshot["population_weight"] = weight
            # The aggregate manifest is authoritative; seed reports retain raw extractions.
            snapshot["population_weight_basis"] = (
                "equal_seed_mass_then_equal_time_decorrelated_frame_occupancy"
            )
    report = {
        "schema_version": TRAJECTORY_SCHEMA_VERSION,
        "kind": "restrained_finite_temperature_xtb_environment_ensemble",
        "configuration": settings,
        "seed_candidate_ids": [seed["candidate_id"] for seed in seeds],
        "trajectory_summary": {
            "requested": len(seeds),
            "completed": sum(item.get("trajectory_status") == "completed" for item in trajectory_reports),
            "failed": sum(item.get("trajectory_status") == "failed" for item in trajectory_reports),
            "retained_decorrelated_snapshots": len(all_snapshots),
        },
        "population_model": {
            "kind": TRAJECTORY_POPULATION_MODEL,
            "warning": TRAJECTORY_POPULATION_WARNING,
            "weight_basis": "equal_seed_mass_then_equal_time_decorrelated_frame_occupancy",
            "weights_normalized": True,
            "vacuum_xtb_energies_used_as_liquid_populations": False,
        },
        "trajectories": trajectory_reports,
        "candidates": all_snapshots,
    }
    manifest = Path(run_dir) / "xtb-trajectories.json"
    manifest.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    sampling_path = Path(run_dir) / "environment-sampling.json"
    sampling = json.loads(sampling_path.read_text()) if sampling_path.is_file() else {}
    sampling["trajectory_sampling"] = {
        "manifest": str(manifest),
        "configuration": settings,
        "candidate_summary": report["trajectory_summary"],
        "population_model": report["population_model"],
    }
    sampling_path.write_text(json.dumps(sampling, indent=2, sort_keys=True) + "\n")
    return all_snapshots, manifest
