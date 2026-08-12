"""Aggregate a fixed IR benchmark manifest without rerunning calculations."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .benchmark import compare_spectra
from .ir_contracts import (assert_partition_separation, condition_compatibility,
                           validate_reference_entry)


def load_ir_manifest(manifest_path: Path) -> tuple[dict, list[dict]]:
    """Load a legacy development manifest or a strict schema-v2 manifest."""
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    version = manifest.get("schema_version")
    if version not in {1, 2}:
        raise ValueError("Unsupported IR benchmark manifest schema")
    entries = []
    for source in manifest.get("entries", []):
        if version == 2:
            entry = validate_reference_entry(source, manifest_path.parent)
        else:
            entry = {
                **source,
                "partition": "development_unassigned",
                "condition": source.get("reference_conditions") or {
                    "phase": "unspecified", "measurement": {"geometry": "unspecified"},
                },
                "contract_status": "legacy_development_only",
            }
        entries.append(entry)
    assert_partition_separation(entries)
    return manifest, entries


def evaluate_ir_manifest(manifest_path: Path, profile: str) -> dict:
    """Evaluate one profile's existing spectra listed in a manifest.

    The manifest deliberately preserves heterogeneous measurement conditions.
    Its aggregate is a development score, never a claim of instrument-specific
    prediction accuracy.
    """
    manifest_path = Path(manifest_path)
    manifest, entries = load_ir_manifest(manifest_path)
    results, missing = [], []
    seen_ids: set[str] = set()
    partitions: dict[str, int] = {}
    for entry in entries:
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
        partition = entry.get("partition", "development_unassigned")
        partitions[partition] = partitions.get(partition, 0) + 1
        predicted_condition_path = prediction.parent / "experimental-condition.json"
        compatibility = {
            "status": "prediction_condition_unavailable",
            "quantitative_comparison_allowed": False,
            "mismatches": [],
            "unknown_or_unspecified": ["predicted_condition"],
        }
        if predicted_condition_path.is_file():
            compatibility = condition_compatibility(
                json.loads(predicted_condition_path.read_text()), entry.get("condition") or {},
            )
        results.append({
            "id": entry["id"], "smiles": entry["smiles"],
            "partition": partition,
            "reference_conditions": entry.get("condition"),
            "condition_compatibility": compatibility,
            "metrics": metrics,
        })
    if not results:
        raise ValueError(f"No usable spectra are registered for profile '{profile}'")
    correlations = [item["metrics"]["pearson_correlation"] for item in results if item["metrics"]["pearson_correlation"] is not None]
    rmses = [item["metrics"]["rmse"] for item in results]
    return {
        "manifest": str(manifest_path),
        "profile": profile,
        "evaluated_entries": len(results),
        "missing_entries": missing,
        "partitions": partitions,
        "condition_matched_entries": sum(
            item["condition_compatibility"]["quantitative_comparison_allowed"] for item in results
        ),
        "aggregate": {
            "mean_pearson_correlation": float(np.mean(correlations)) if correlations else None,
            "median_pearson_correlation": float(np.median(correlations)) if correlations else None,
            "mean_rmse": float(np.mean(rmses)) if rmses else None,
        },
        "entries": results,
        "limitations": (
            "The aggregate includes trace-shape development metrics. Only entries with "
            "condition_compatibility.quantitative_comparison_allowed=true support condition-matched claims."
        ),
    }


def write_ir_benchmark(path: Path, report: dict) -> Path:
    Path(path).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return Path(path)
