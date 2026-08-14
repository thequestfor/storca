"""Geometry/occupancy-driven adaptive acquisition of direct local DFT modes."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from .direct_local_dft import coordination_class, spectrum_band_metrics
from .spectrum import write_spectrum_csv


@dataclass(frozen=True)
class AdaptiveLocalModeConfig:
    schema_version: int = 1
    initial_representatives: int = 6
    batch_size: int = 3
    initial_target_representatives: int = 12
    maximum_representatives: int = 18
    invocations_per_oh_bond: int = 4
    ncores_per_invocation: int = 8
    required_consecutive_passes: int = 2
    center_tolerance_cm_1: float = 5.0
    fwhm_absolute_tolerance_cm_1: float = 10.0
    fwhm_relative_tolerance: float = 0.10
    minimum_envelope_overlap: float = 0.98
    maximum_center_ci_width_cm_1: float = 15.0
    minimum_effective_sample_size_per_important_class: float = 5.0
    minimum_important_class_occupancy_coverage: float = 0.80
    important_class_minimum_occupancy: float = 0.05
    bootstrap_samples: int = 500
    bootstrap_seed: int = 29815


def _feature(record: dict, key: str, default: float = 0.0) -> float:
    value = (record.get("environment_features") or {}).get(key, default)
    try:
        value = float(value)
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def _classes(record: dict) -> set[str]:
    return {coordination_class(item) for item in record.get("local_stretch_bonds") or []}


def select_adaptive_local_mode_representatives(
    records: list[dict], distance_matrix: np.ndarray, representative_count: int, *,
    already_selected_ids: Iterable[str] = (),
) -> tuple[list[int], dict]:
    """Select using occupancy and geometry only; xTB frequencies are never accepted."""
    matrix = np.asarray(distance_matrix, dtype=float)
    if matrix.shape != (len(records), len(records)) or not np.all(np.isfinite(matrix)):
        raise ValueError("Adaptive local-mode selection requires a finite square distance matrix")
    if not 0 <= representative_count <= len(records):
        raise ValueError("Representative count is outside the candidate pool")
    by_id = {str(item["candidate_id"]): index for index, item in enumerate(records)}
    selected = list(dict.fromkeys(
        by_id[str(item)] for item in already_selected_ids if str(item) in by_id
    ))
    class_occupancy: dict[str, float] = {}
    for record in records:
        occupancy = max(0.0, float(record.get("population_weight", 0.0)))
        for key in _classes(record):
            class_occupancy[key] = class_occupancy.get(key, 0.0) + occupancy
    decisions = []
    while len(selected) < representative_count:
        candidates = [index for index in range(len(records)) if index not in selected]
        selected_classes = {key for index in selected for key in _classes(records[index])}
        selected_topologies = {str(records[index].get("topology", "unknown")) for index in selected}
        selected_distances = [_feature(records[index], "h_bond_distance_angstrom") for index in selected]
        selected_angles = [_feature(records[index], "donor_h_acceptor_angle_degrees") for index in selected]

        def score(index: int) -> tuple[float, float, str]:
            record = records[index]
            occupancy = max(0.0, float(record.get("population_weight", 0.0)))
            uncovered = sum(class_occupancy.get(key, 0.0) for key in _classes(record) - selected_classes)
            topology = 1.0 if str(record.get("topology", "unknown")) not in selected_topologies else 0.0
            distance_novelty = min(
                (abs(_feature(record, "h_bond_distance_angstrom") - value) for value in selected_distances),
                default=0.0,
            )
            angle_novelty = min(
                (abs(_feature(record, "donor_h_acceptor_angle_degrees") - value) / 180.0 for value in selected_angles),
                default=0.0,
            )
            separation = min((float(matrix[index, old]) for old in selected), default=0.0)
            total = 4.0 * occupancy + 3.0 * uncovered + topology + distance_novelty + angle_novelty + separation
            return total, separation, str(record["candidate_id"])

        chosen = max(candidates, key=score)
        record = records[chosen]
        selected.append(chosen)
        decisions.append({
            "position": len(selected), "candidate_id": record["candidate_id"],
            "trajectory_occupancy": record.get("population_weight", 0.0),
            "coordination_classes": sorted(_classes(record)),
            "h_bond_distance_angstrom": _feature(record, "h_bond_distance_angstrom"),
            "h_bond_angle_degrees": _feature(record, "donor_h_acceptor_angle_degrees"),
            "cluster_topology": record.get("topology"),
            "distance_from_calculated_environments": score(chosen)[1],
            "selection_score": score(chosen)[0],
        })
    return selected, {
        "kind": "geometry_occupancy_direct_local_dft_acquisition",
        "selection_inputs": [
            "trajectory_occupancy", "donor_acceptor_coordination", "h_bond_distance",
            "h_bond_angle", "cluster_topology", "geometric_diversity",
            "distance_from_already_calculated_environments",
        ],
        "prohibited_selection_input": "xtb_frequency",
        "xtb_frequency_used": False,
        "selected_candidate_ids": [records[index]["candidate_id"] for index in selected],
        "decisions": decisions,
    }


def _cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left, right) / denominator) if denominator else 0.0


def cluster_bootstrap_center_ci(
    observations: list[dict], *, samples: int = 500, seed: int = 29815,
) -> dict:
    """Resample trajectory clusters as blocks, never individual snapshots."""
    clusters: dict[str, list[dict]] = {}
    for item in observations:
        cluster = str(item.get("trajectory_cluster_id") or item.get("independent_environment_id"))
        clusters.setdefault(cluster, []).append(item)
    if not clusters:
        return {"status": "insufficient_clusters", "ci95_cm-1": [None, None], "width_cm-1": None}
    keys = sorted(clusters)
    rng = np.random.default_rng(seed)
    centers = []
    for _ in range(samples):
        chosen = rng.choice(keys, size=len(keys), replace=True)
        rows = [item for key in chosen for item in clusters[str(key)]]
        weights = np.asarray([max(0.0, float(item.get("weight", 0.0))) for item in rows])
        frequencies = np.asarray([float(item["frequency_cm-1"]) for item in rows])
        if weights.sum() <= 0:
            weights[:] = 1.0
        centers.append(float(np.average(frequencies, weights=weights)))
    low, high = np.quantile(centers, [0.025, 0.975])
    return {
        "status": "completed", "resampling_unit": "trajectory_cluster",
        "clusters": len(keys), "samples": samples,
        "ci95_cm-1": [float(low), float(high)], "width_cm-1": float(high - low),
    }


def evaluate_adaptive_batches(
    batches: list[dict], *, config: AdaptiveLocalModeConfig | None = None,
) -> dict:
    """Apply all convergence gates for two consecutive three-representative batches."""
    resolved = config or AdaptiveLocalModeConfig()
    reports, consecutive, converged = [], 0, False
    previous = None
    for batch in batches:
        grid = np.asarray(batch["grid"], dtype=float)
        spectrum = np.asarray(batch["spectrum"], dtype=float)
        metrics = spectrum_band_metrics(grid, spectrum, 2800.0, 3900.0)
        observations = list(batch.get("observations") or [])
        bootstrap = cluster_bootstrap_center_ci(
            observations, samples=resolved.bootstrap_samples, seed=resolved.bootstrap_seed + len(reports),
        )
        comparison = None
        if previous is not None:
            center_change = abs(float(metrics["center_cm-1"]) - float(previous["metrics"]["center_cm-1"]))
            width_change = abs(float(metrics["fwhm_cm-1"]) - float(previous["metrics"]["fwhm_cm-1"]))
            width_relative = width_change / max(float(previous["metrics"]["fwhm_cm-1"]), 1.0)
            overlap = _cosine(previous["spectrum"], spectrum)
            class_reports = batch.get("coordination_classes") or {}
            important = [item for item in class_reports.values() if float(item.get("occupancy_fraction", 0.0)) >= resolved.important_class_minimum_occupancy]
            gates = {
                "center": center_change <= resolved.center_tolerance_cm_1,
                "fwhm": width_change <= resolved.fwhm_absolute_tolerance_cm_1 and width_relative <= resolved.fwhm_relative_tolerance,
                "envelope_overlap": overlap >= resolved.minimum_envelope_overlap,
                "bootstrap_center_ci": bootstrap.get("width_cm-1") is not None and bootstrap["width_cm-1"] <= resolved.maximum_center_ci_width_cm_1,
                "effective_sample_size": bool(important) and all(float(item.get("effective_sample_size", 0.0)) >= resolved.minimum_effective_sample_size_per_important_class for item in important),
                "occupancy_coverage": bool(important) and all(float(item.get("occupancy_coverage_fraction", 0.0)) >= resolved.minimum_important_class_occupancy_coverage for item in important),
            }
            passed = all(gates.values())
            comparison = {
                "center_change_cm-1": center_change, "fwhm_change_cm-1": width_change,
                "fwhm_relative_change": width_relative, "whole_oh_envelope_overlap": overlap,
                "gates": gates, "passed": passed,
            }
            consecutive = consecutive + 1 if passed else 0
        report = {
            "representatives": int(batch["representatives"]), "metrics": metrics,
            "bootstrap": bootstrap, "comparison_to_previous": comparison,
            "consecutive_passing_batches": consecutive,
        }
        reports.append(report)
        previous = {"metrics": metrics, "spectrum": spectrum}
        if consecutive >= resolved.required_consecutive_passes:
            converged = True
            break
    evaluated = reports[-1]["representatives"] if reports else 0
    return {
        "schema_version": 1, "kind": "adaptive_direct_local_dft_convergence",
        "configuration": asdict(resolved), "batches": reports,
        "summary": {
            "status": "converged" if converged else "hard_budget_reached_unconverged" if evaluated >= resolved.maximum_representatives else "additional_batch_required",
            "converged": converged, "representatives_evaluated": evaluated,
            "promotion_to_bulk_embedding": not converged and evaluated >= resolved.maximum_representatives,
        },
    }


def run_adaptive_local_mode_acquisition(
    run_dir: Path, selected_records: list[dict], *,
    acquire_batch: Callable[[list[dict], int], dict],
    config: AdaptiveLocalModeConfig | None = None,
) -> dict:
    """Run resumable batches of three through a caller-supplied ORCA acquisition callback."""
    resolved = config or AdaptiveLocalModeConfig()
    run_dir = Path(run_dir)
    artifact = run_dir / "local-mode-acquisition.json"
    retained = json.loads(artifact.read_text()) if artifact.is_file() else {"completed_batches": []}
    batches = list(retained.get("completed_batches") or [])
    start = max(resolved.initial_representatives, batches[-1]["representatives"] if batches else 0)
    payload = {
        "schema_version": 1, "kind": "adaptive_local_mode_acquisition",
        "configuration": asdict(resolved), "frequency_selection_source": "geometry_and_occupancy_only",
        "completed_batches": batches,
        "convergence": evaluate_adaptive_batches(batches, config=resolved),
        "final_representative_manifest": [item.get("candidate_id") for item in selected_records[:start]],
    }
    for endpoint in range(start + resolved.batch_size, resolved.maximum_representatives + 1, resolved.batch_size):
        result = acquire_batch(selected_records[:endpoint], endpoint)
        batches = [item for item in batches if int(item["representatives"]) != endpoint]
        batches.append(result)
        batches.sort(key=lambda item: int(item["representatives"]))
        convergence = evaluate_adaptive_batches(batches, config=resolved)
        payload = {
            "schema_version": 1, "kind": "adaptive_local_mode_acquisition",
            "configuration": asdict(resolved), "frequency_selection_source": "geometry_and_occupancy_only",
            "completed_batches": batches, "convergence": convergence,
            "final_representative_manifest": [item.get("candidate_id") for item in selected_records[:endpoint]],
        }
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        batch_path = run_dir / f"local-mode-convergence-batch-{endpoint:02d}.json"
        batch_path.write_text(json.dumps(convergence["batches"][-1], indent=2, sort_keys=True) + "\n")
        uncertainty_path = run_dir / "local-mode-uncertainty-bands.json"
        uncertainty_path.write_text(json.dumps({
            "schema_version": 1, "resampling_unit": "trajectory_cluster",
            "batches": [
                {"representatives": item["representatives"], "band_center_bootstrap": item["bootstrap"]}
                for item in convergence["batches"]
            ],
        }, indent=2, sort_keys=True) + "\n")
        manifest_path = run_dir / "final-local-mode-representative-manifest.json"
        manifest_path.write_text(json.dumps({
            "schema_version": 1, "representatives": selected_records[:endpoint],
            "selection_frequency_source": "none_geometry_and_occupancy_only",
        }, indent=2, sort_keys=True) + "\n")
        payload["artifacts"] = {
            "uncertainty_bands": str(uncertainty_path),
            "final_representative_manifest": str(manifest_path),
            "per_batch_convergence_report": str(batch_path),
        }
        if convergence["summary"]["converged"]:
            final_batch = batches[-1]
            spectrum_path = write_spectrum_csv(
                run_dir / "spectrum_direct_local_dft_converged.csv",
                np.asarray(final_batch["grid"], dtype=float),
                np.asarray(final_batch["spectrum"], dtype=float), "absorption_strength",
            )
            payload["artifacts"]["converged_direct_local_dft_spectrum"] = str(spectrum_path)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        if convergence["summary"]["converged"]:
            break
    payload["artifact"] = str(artifact)
    return payload
