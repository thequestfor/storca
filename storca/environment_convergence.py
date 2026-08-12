"""Adaptive convergence for diversity-selected environment Hessians."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class EnvironmentConvergenceConfig:
    schema_version: int = 1
    initial_batch_size: int = 2
    required_consecutive_passes: int = 2
    center_tolerance_cm_1: float = 5.0
    width_absolute_tolerance_cm_1: float = 10.0
    width_relative_tolerance: float = 0.10
    minimum_spectrum_overlap: float = 0.98
    important_band_relative_intensity: float = 0.05
    comparison_fwhm_cm_1: float = 2.0
    band_matching_window_cm_1: float = 75.0
    minimum_mode_overlap: float = 0.45


def convergence_schedule(fidelity: str, representative_count: int) -> list[int]:
    """Return cumulative, budget-capped representative batch endpoints."""
    if representative_count < 1:
        return []
    fidelity = fidelity.strip().lower()
    templates = {
        "auto": [2, 3, 4],
        "balanced": [2, 4, 6],
        "accurate": [3, 6, 9, 12],
    }
    if fidelity not in templates:
        return [representative_count]
    endpoints = [min(value, representative_count) for value in templates[fidelity]]
    endpoints.append(representative_count)
    return sorted(set(endpoints))


def spectrum_cosine_overlap(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape or left.ndim != 1 or not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise ValueError("Spectrum overlap requires equal finite one-dimensional arrays")
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left, right) / denominator) if denominator > 0 else 0.0


def _important_bands(bands: list[dict], threshold: float) -> list[dict]:
    maximum = max((float(band.get("population_weighted_intensity", 0.0)) for band in bands), default=0.0)
    if maximum <= 0:
        return []
    return [
        band for band in bands
        if float(band.get("population_weighted_intensity", 0.0)) >= threshold * maximum
    ]


def _match_bands(previous: list[dict], current: list[dict], window: float) -> tuple[list[tuple[dict, dict]], list[dict]]:
    """Greedily match important bands by mode identity, then frequency proximity."""
    available = set(range(len(previous)))
    pairs: list[tuple[dict, dict]] = []
    unmatched: list[dict] = []
    for band in sorted(current, key=lambda item: -float(item["population_weighted_intensity"])):
        same_mode = [index for index in available if previous[index].get("mode") == band.get("mode")]
        candidates = same_mode or list(available)
        if not candidates:
            unmatched.append(band)
            continue
        index = min(
            candidates,
            key=lambda candidate: abs(
                float(previous[candidate]["mean_frequency_cm-1"])
                - float(band["mean_frequency_cm-1"])
            ),
        )
        gap = abs(float(previous[index]["mean_frequency_cm-1"]) - float(band["mean_frequency_cm-1"]))
        if gap > window:
            unmatched.append(band)
            continue
        available.remove(index)
        pairs.append((previous[index], band))
    return pairs, unmatched


def compare_environment_batches(
    previous: dict,
    current: dict,
    *,
    config: EnvironmentConvergenceConfig | None = None,
) -> dict:
    resolved = config or EnvironmentConvergenceConfig()
    previous_important = _important_bands(previous["bands"], resolved.important_band_relative_intensity)
    current_important = _important_bands(current["bands"], resolved.important_band_relative_intensity)
    pairs, new_bands = _match_bands(
        previous_important, current_important, resolved.band_matching_window_cm_1,
    )
    band_changes = []
    center_pass = bool(pairs) and not new_bands
    width_pass = bool(pairs) and not new_bands
    mode_overlap_pass = bool(current_important)
    for old, new in pairs:
        center_change = abs(float(new["mean_frequency_cm-1"]) - float(old["mean_frequency_cm-1"]))
        old_width = float(old["ensemble_fwhm_equivalent_cm-1"])
        new_width = float(new["ensemble_fwhm_equivalent_cm-1"])
        width_change = abs(new_width - old_width)
        relative_change = width_change / max(abs(old_width), 1.0)
        center_ok = center_change < resolved.center_tolerance_cm_1
        width_ok = (
            width_change < resolved.width_absolute_tolerance_cm_1
            or relative_change < resolved.width_relative_tolerance
        )
        overlap = new.get("minimum_displacement_overlap")
        overlap_ok = isinstance(overlap, (int, float)) and float(overlap) >= resolved.minimum_mode_overlap
        center_pass = center_pass and center_ok
        width_pass = width_pass and width_ok
        mode_overlap_pass = mode_overlap_pass and overlap_ok
        band_changes.append({
            "mode": new.get("mode"),
            "previous_center_cm-1": old["mean_frequency_cm-1"],
            "current_center_cm-1": new["mean_frequency_cm-1"],
            "center_change_cm-1": center_change,
            "previous_fwhm_cm-1": old_width,
            "current_fwhm_cm-1": new_width,
            "width_change_cm-1": width_change,
            "width_relative_change": relative_change,
            "center_pass": center_ok,
            "width_pass": width_ok,
            "mode_overlap_pass": overlap_ok,
        })
    spectrum_overlap = spectrum_cosine_overlap(previous["spectrum"], current["spectrum"])
    spectrum_pass = spectrum_overlap > resolved.minimum_spectrum_overlap
    no_new_band_class = not new_bands
    integrated_intensity = float(current["integrated_intensity"])
    previous_integrated = float(previous["integrated_intensity"])
    intensity_change = abs(integrated_intensity - previous_integrated)
    passed = all((center_pass, width_pass, spectrum_pass, no_new_band_class, mode_overlap_pass))
    return {
        "passed": passed,
        "center_pass": center_pass,
        "width_pass": width_pass,
        "spectrum_pass": spectrum_pass,
        "mode_overlap_pass": mode_overlap_pass,
        "no_new_significant_band_class": no_new_band_class,
        "whole_spectrum_overlap": spectrum_overlap,
        "integrated_intensity_change": intensity_change,
        "integrated_intensity_relative_change": intensity_change / max(abs(previous_integrated), 1.0e-12),
        "important_bands_previous": len(previous_important),
        "important_bands_current": len(current_important),
        "new_significant_bands": [band["mean_frequency_cm-1"] for band in new_bands],
        "band_changes": band_changes,
    }


def _batch_snapshot(records: list[dict], scale_factor: float, config: EnvironmentConvergenceConfig) -> dict:
    from .experimental_ir import ensemble_band_statistics
    from .spectrum import build_ir_spectrum

    grid, spectrum, successful = build_ir_spectrum(
        records, scale_factor=scale_factor, fwhm=config.comparison_fwhm_cm_1,
    )
    bands = ensemble_band_statistics(
        successful, scale_factor=scale_factor, environment_sampling=True,
    )
    return {
        "records": records,
        "successful_environments": len(successful),
        "grid": grid,
        "spectrum": spectrum,
        "integrated_intensity": float(np.trapezoid(spectrum, grid)),
        "bands": bands,
    }


def _serializable_batch(snapshot: dict) -> dict:
    return {
        "successful_environments": snapshot["successful_environments"],
        "integrated_intensity": snapshot["integrated_intensity"],
        "bands": snapshot["bands"],
    }


def _manifest_sha256(run_dir: Path) -> str | None:
    path = Path(run_dir) / "clusters" / "selected-conformers.json"
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def _write_report(run_dir: Path, report: dict) -> Path:
    path = Path(run_dir) / "environment-convergence.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    sampling_path = Path(run_dir) / "environment-sampling.json"
    sampling = json.loads(sampling_path.read_text()) if sampling_path.is_file() else {}
    sampling["schema_version"] = max(3, int(sampling.get("schema_version", 0)))
    sampling["adaptive_convergence"] = {
        "artifact": str(path),
        "status": report["summary"]["status"],
        "converged": report["summary"]["converged"],
        "stop_reason": report["summary"]["stop_reason"],
        "environments_evaluated": report["summary"]["environments_evaluated"],
    }
    sampling_path.write_text(json.dumps(sampling, indent=2, sort_keys=True) + "\n")
    return path


def evaluate_completed_environment_batches(
    record_batches: list[tuple[int, list[dict]]],
    run_dir: Path,
    *,
    fidelity: str,
    scale_factor: float,
    representative_budget: int,
    budget_exhausted: bool,
    config: EnvironmentConvergenceConfig | None = None,
) -> dict:
    """Evaluate cumulative completed batches and persist an auditable report."""
    resolved = config or EnvironmentConvergenceConfig()
    report_batches: list[dict] = []
    previous = None
    consecutive = 0
    converged = False
    evaluated = 0
    for requested, records in record_batches:
        try:
            snapshot = _batch_snapshot(records, scale_factor, resolved)
            evaluated = max(evaluated, requested)
            comparison = (
                compare_environment_batches(previous, snapshot, config=resolved)
                if previous is not None else None
            )
            if comparison and comparison["passed"] and snapshot["successful_environments"] >= 3:
                consecutive += 1
            elif comparison:
                consecutive = 0
            report_batches.append({
                "cumulative_environments_requested": requested,
                **_serializable_batch(snapshot),
                "comparison_to_previous": comparison,
                "consecutive_passing_comparisons": consecutive,
            })
            previous = snapshot
            if consecutive >= resolved.required_consecutive_passes:
                converged = True
                break
        except Exception as error:
            consecutive = 0
            report_batches.append({
                "cumulative_environments_requested": requested,
                "successful_environments": 0,
                "error": str(error),
                "comparison_to_previous": None,
                "consecutive_passing_comparisons": consecutive,
            })
    final_successful = report_batches[-1].get("successful_environments", 0) if report_batches else 0
    if converged:
        status, stop_reason = "converged", "two_consecutive_batches_passed"
    elif final_successful < 3:
        status, stop_reason = "insufficient_environment_sampling", "fewer_than_three_successful_environments"
    elif budget_exhausted:
        status, stop_reason = "environment_width_unconverged", "hard_environment_orca_budget_reached"
    else:
        status, stop_reason = "adaptive_convergence_incomplete", "not_all_planned_batches_evaluated"
    report = {
        "schema_version": 1,
        "kind": "adaptive_environment_convergence",
        "configuration": asdict(resolved),
        "fidelity": fidelity,
        "representative_manifest_sha256": _manifest_sha256(run_dir),
        "representative_budget": representative_budget,
        "batches": report_batches,
        "summary": {
            "status": status,
            "converged": converged,
            "stop_reason": stop_reason,
            "environments_evaluated": evaluated,
            "successful_environments": final_successful,
            "consecutive_passing_comparisons": consecutive,
            "budget_exhausted": bool(budget_exhausted),
            "converged_important_modes": (
                [
                    band.get("mode") for band in _important_bands(
                        previous["bands"], resolved.important_band_relative_intensity,
                    )
                ] if converged and previous is not None else []
            ),
        },
    }
    _write_report(run_dir, report)
    return report


def run_adaptive_environment_convergence(
    cluster_paths: list[Path],
    cluster_dir: Path,
    run_dir: Path,
    *,
    fidelity: str,
    scale_factor: float,
    run_batch: Callable[[list[Path]], dict],
    progress: Callable[[str], None] | None = None,
    after_batch: Callable[[int, dict], list[Path] | None] | None = None,
    config: EnvironmentConvergenceConfig | None = None,
) -> tuple[dict, dict]:
    """Submit cumulative representative prefixes until stable or budget-limited."""
    resolved = config or EnvironmentConvergenceConfig()
    endpoints = convergence_schedule(fidelity, len(cluster_paths))
    record_batches: list[tuple[int, list[dict]]] = []
    last_result: dict | None = None
    for endpoint in endpoints:
        if progress:
            progress(f"Environment convergence batch: {endpoint}/{len(cluster_paths)} representatives")
        try:
            last_result = run_batch(cluster_paths[:endpoint])
            record_batches.append((endpoint, last_result["conformers"]))
            if after_batch is not None:
                reordered = after_batch(endpoint, last_result)
                if reordered:
                    if len(reordered) != len(cluster_paths):
                        raise RuntimeError("Adaptive acquisition changed the hard representative budget")
                    cluster_paths[:] = reordered
        except Exception as error:
            record_batches.append((endpoint, []))
            if progress:
                progress(f"Environment convergence batch {endpoint} produced no usable spectrum ({error})")
        provisional = evaluate_completed_environment_batches(
            record_batches, run_dir, fidelity=fidelity, scale_factor=scale_factor,
            representative_budget=len(cluster_paths), budget_exhausted=endpoint == len(cluster_paths),
            config=resolved,
        )
        if provisional["summary"]["converged"]:
            break
    if last_result is None:
        raise RuntimeError("No environment convergence batch produced a usable ORCA spectrum")
    return last_result, provisional
