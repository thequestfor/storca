"""Transparent, held-out calibration of harmonic IR frequency scaling."""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import numpy as np

from .analysis import analyze_ir_spectra
from .benchmark import compare_spectra
from .spectrum import build_ir_spectrum, format_spectrum, load_completed_ir_records, write_spectrum_csv


def scale_grid(start: float, stop: float, step: float) -> list[float]:
    """Create a stable inclusive scale-factor grid."""
    if (not all(math.isfinite(value) for value in (start, stop, step))
            or not 0.1 <= start <= stop <= 1.2 or step < 1e-6):
        raise ValueError("Scale range must be within 0.1..1.2 and step must be at least 1e-6")
    values = [start]
    while values[-1] + step < stop - 1e-12:
        values.append(values[-1] + step)
    values.append(stop)
    rounded = [round(float(value), 9) for value in values]
    unique = list(dict.fromkeys(rounded))
    if len(unique) < 2 and start != stop:
        raise ValueError("Scale step is too small to produce distinct retained factors")
    return unique


def _manifest_entries(manifest_path: Path, profile: str) -> list[dict]:
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported IR benchmark manifest schema")
    entries = []
    seen_ids: set[str] = set()
    for entry in manifest.get("entries", []):
        entry_id = entry.get("id")
        if not entry_id or entry_id in seen_ids:
            raise ValueError("IR benchmark manifest IDs must be present and unique")
        seen_ids.add(entry_id)
        relative_prediction = entry.get("baseline_spectra", {}).get(profile)
        if not relative_prediction:
            continue
        prediction = (manifest_path.parent / relative_prediction).resolve()
        reference = (manifest_path.parent / entry["reference_spectrum"]).resolve()
        if not prediction.is_file() or not reference.is_file():
            continue
        entries.append({"id": entry["id"], "run_dir": prediction.parent, "reference": reference})
    return entries


def _score_scale(entries: list[dict], scale_factor: float, fwhm: float, temperature: float, scratch: Path) -> dict:
    results = []
    for entry in entries:
        records = entry["records"]
        grid, intensity, _ = build_ir_spectrum(records, scale_factor=scale_factor, fwhm=fwhm)
        display, column = format_spectrum(intensity, style="transmittance")
        predicted = scratch / f"{entry['id']}-{scale_factor:.6f}.csv"
        write_spectrum_csv(predicted, grid, display, column)
        trace = compare_spectra(predicted, entry["reference"])
        peaks = analyze_ir_spectra(predicted, entry["reference"])
        results.append({
            "id": entry["id"],
            "pearson_correlation": trace["pearson_correlation"],
            "rmse": trace["rmse"],
            "matched_bands": len(peaks["matches"]),
            "reference_bands": len(peaks["reference_peaks"]),
            "unmatched_reference_bands": len(peaks["unmatched_reference"]),
            "mean_absolute_position_error_cm-1": peaks["mean_absolute_position_error_cm-1"],
        })
    correlations = [row["pearson_correlation"] for row in results if row["pearson_correlation"] is not None]
    position_errors = [row["mean_absolute_position_error_cm-1"] for row in results if row["mean_absolute_position_error_cm-1"] is not None]
    total_matched = sum(row["matched_bands"] for row in results)
    total_unmatched_reference = sum(row["unmatched_reference_bands"] for row in results)
    matched_error_sum = sum(
        (row["mean_absolute_position_error_cm-1"] or 0.0) * row["matched_bands"] for row in results
    )
    unmatched_penalty_cm_1 = 40.0
    scored_reference_bands = total_matched + total_unmatched_reference
    penalized_error = (
        (matched_error_sum + unmatched_penalty_cm_1 * total_unmatched_reference) / scored_reference_bands
        if scored_reference_bands else None
    )
    return {
        "scale_factor": scale_factor,
        "entries": results,
        "aggregate": {
            "mean_pearson_correlation": float(np.mean(correlations)) if correlations else None,
            "mean_rmse": float(np.mean([row["rmse"] for row in results])) if results else None,
            "mean_absolute_position_error_cm-1": float(np.mean(position_errors)) if position_errors else None,
            "coverage_penalized_position_error_cm-1": penalized_error,
            "unmatched_reference_penalty_cm-1": unmatched_penalty_cm_1,
            "matched_bands": int(total_matched),
            "unmatched_reference_bands": int(total_unmatched_reference),
        },
    }


def calibrate_harmonic_scale(
    manifest_path: Path,
    *,
    profile: str,
    training_ids: list[str],
    holdout_ids: list[str],
    scale_start: float = 0.94,
    scale_stop: float = 1.00,
    scale_step: float = 0.002,
    fwhm: float = 15.0,
    temperature: float = 298.15,
) -> dict:
    """Choose a global scale on training molecules and report untouched holdouts.

    No ORCA job is executed and no existing run artifact is overwritten. The
    selection criterion is mean prominent-band position error, then correlation
    as a deterministic tiebreaker. It is deliberately not a learned
    SMILES-to-spectrum model.
    """
    if not training_ids or not holdout_ids:
        raise ValueError("Provide at least one training ID and one holdout ID")
    if len(set(training_ids)) != len(training_ids) or len(set(holdout_ids)) != len(holdout_ids):
        raise ValueError("Training and holdout IDs must not contain duplicates")
    if set(training_ids) & set(holdout_ids):
        raise ValueError("Training and holdout IDs must not overlap")
    if not all(math.isfinite(value) for value in (fwhm, temperature)) or fwhm <= 0 or temperature <= 0:
        raise ValueError("FWHM and temperature must be positive")
    entries = _manifest_entries(Path(manifest_path), profile)
    by_id = {entry["id"]: entry for entry in entries}
    unknown = sorted((set(training_ids) | set(holdout_ids)) - set(by_id))
    if unknown:
        raise ValueError(f"No usable '{profile}' spectrum is registered for: {', '.join(unknown)}")
    requested_ids = set(training_ids) | set(holdout_ids)
    for entry_id in requested_ids:
        entry = by_id[entry_id]
        entry["records"] = load_completed_ir_records(entry["run_dir"], temperature=temperature)
    training = [by_id[item] for item in training_ids]
    holdout = [by_id[item] for item in holdout_ids]
    with tempfile.TemporaryDirectory(prefix="storca-scale-") as temporary:
        scratch = Path(temporary)
        candidates = [_score_scale(training, value, fwhm, temperature, scratch) for value in scale_grid(scale_start, scale_stop, scale_step)]
        def rank(candidate: dict) -> tuple[float, int, float]:
            aggregate = candidate["aggregate"]
            error = aggregate["coverage_penalized_position_error_cm-1"]
            correlation = aggregate["mean_pearson_correlation"]
            return (float("inf") if error is None else error, -aggregate["matched_bands"],
                    -(float("-inf") if correlation is None else correlation))
        selected = min(candidates, key=rank)
        if selected["aggregate"]["coverage_penalized_position_error_cm-1"] is None:
            raise ValueError("Training references contained no detectable bands for calibration")
        holdout_report = _score_scale(holdout, selected["scale_factor"], fwhm, temperature, scratch)
    return {
        "schema_version": 1,
        "kind": "transparent_harmonic_scale_calibration",
        "profile": profile,
        "selection": {
            "criterion": "minimum coverage-penalized prominent-band position error; matched-band count and correlation are tiebreakers",
            "selected_scale_factor": selected["scale_factor"],
            "training_metrics": selected["aggregate"],
        },
        "settings": {"fwhm_cm-1": fwhm, "temperature_K": temperature, "scale_grid": scale_grid(scale_start, scale_stop, scale_step)},
        "training_ids": training_ids,
        "holdout_ids": holdout_ids,
        "candidates": candidates,
        "holdout": holdout_report,
        "limitations": "A global harmonic scale does not model solvent, phase, hydrogen bonding, anharmonicity, or instrument line shape. References must be compatible with the intended use case before accepting this calibration.",
    }


def write_calibration(path: Path, report: dict) -> Path:
    path = Path(path)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return path
