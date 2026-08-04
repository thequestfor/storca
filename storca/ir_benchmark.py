"""Aggregate a fixed IR benchmark manifest without rerunning calculations."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .benchmark import compare_spectra


def evaluate_ir_manifest(manifest_path: Path, profile: str) -> dict:
    """Evaluate one profile's existing spectra listed in a manifest.

    The manifest deliberately preserves heterogeneous measurement conditions.
    Its aggregate is a development score, never a claim of instrument-specific
    prediction accuracy.
    """
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported IR benchmark manifest schema")
    results, missing = [], []
    seen_ids: set[str] = set()
    for entry in manifest.get("entries", []):
        entry_id = entry.get("id")
        if not entry_id or entry_id in seen_ids:
            raise ValueError("IR benchmark manifest IDs must be present and unique")
        seen_ids.add(entry_id)
        spectrum_rel = entry.get("baseline_spectra", {}).get(profile)
        if not spectrum_rel:
            missing.append({"id": entry.get("id"), "reason": f"No spectrum registered for profile '{profile}'"})
            continue
        prediction = (manifest_path.parent / spectrum_rel).resolve()
        reference = (manifest_path.parent / entry["reference_spectrum"]).resolve()
        if not prediction.is_file() or not reference.is_file():
            missing.append({"id": entry.get("id"), "reason": "Prediction or reference file is missing", "prediction": str(prediction), "reference": str(reference)})
            continue
        metrics = compare_spectra(prediction, reference)
        results.append({"id": entry["id"], "smiles": entry["smiles"],
                        "reference_conditions": entry.get("reference_conditions"), "metrics": metrics})
    if not results:
        raise ValueError(f"No usable spectra are registered for profile '{profile}'")
    correlations = [item["metrics"]["pearson_correlation"] for item in results if item["metrics"]["pearson_correlation"] is not None]
    rmses = [item["metrics"]["rmse"] for item in results]
    return {
        "manifest": str(manifest_path),
        "profile": profile,
        "evaluated_entries": len(results),
        "missing_entries": missing,
        "aggregate": {
            "mean_pearson_correlation": float(np.mean(correlations)) if correlations else None,
            "median_pearson_correlation": float(np.median(correlations)) if correlations else None,
            "mean_rmse": float(np.mean(rmses)) if rmses else None,
        },
        "entries": results,
        "limitations": "Aggregate scores span heterogeneous experimental references and are for method-ranking development only.",
    }


def write_ir_benchmark(path: Path, report: dict) -> Path:
    Path(path).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return Path(path)
