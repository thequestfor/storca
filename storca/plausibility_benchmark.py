"""Validation harness for condition-specific molecular persistence assessments."""

from __future__ import annotations

import json
from pathlib import Path


PERSISTENCE_CATEGORIES = {
    "ordinary_condition_persistent",
    "reactive_but_isolable",
    "transient",
    "special_conditions_only",
    "not_supported",
    "insufficient_evidence",
}
NONPERSISTENT_CATEGORIES = {
    "reactive_but_isolable",
    "transient",
    "special_conditions_only",
    "not_supported",
}


def load_plausibility_manifest(path: Path) -> dict:
    """Read and validate a benchmark manifest without accepting draft labels."""
    path = Path(path)
    manifest = json.loads(path.read_text())
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported plausibility benchmark manifest schema")
    identifiers = set()
    for entry in manifest.get("entries", []):
        identifier = entry.get("id")
        if not identifier or identifier in identifiers:
            raise ValueError("Every plausibility benchmark entry needs a unique id")
        identifiers.add(identifier)
        status = entry.get("evidence_status")
        category = entry.get("expected_category")
        if status == "accepted":
            if category not in PERSISTENCE_CATEGORIES - {"insufficient_evidence"}:
                raise ValueError(f"Accepted entry '{identifier}' needs a concrete expected_category")
            if not entry.get("evidence_note") or not entry.get("evidence_source"):
                raise ValueError(f"Accepted entry '{identifier}' needs an evidence note and source")
        elif status != "draft":
            raise ValueError(f"Entry '{identifier}' must have evidence_status 'draft' or 'accepted'")
    return manifest


def _result_category(result: dict) -> str | None:
    assessment = result.get("assessment", {})
    category = assessment.get("persistence_category")
    return category if category in PERSISTENCE_CATEGORIES else None


def evaluate_plausibility_manifest(manifest_path: Path, results_dir: Path) -> dict:
    """Score completed dossiers against accepted, condition-matched references.

    A missing or malformed dossier is reported as missing evidence. It is never
    converted into a persistent prediction. False reassurance is reported
    separately because it is the most consequential benchmark error.
    """
    manifest = load_plausibility_manifest(manifest_path)
    results_dir = Path(results_dir)
    entries, missing = [], []
    for reference in manifest.get("entries", []):
        if reference["evidence_status"] != "accepted":
            continue
        result_path = results_dir / f"{reference['id']}.json"
        if not result_path.is_file():
            missing.append({"id": reference["id"], "reason": "No plausibility dossier"})
            continue
        try:
            prediction = _result_category(json.loads(result_path.read_text()))
        except json.JSONDecodeError:
            prediction = None
        if prediction is None:
            missing.append({"id": reference["id"], "reason": "Dossier has no valid persistence category"})
            continue
        expected = reference["expected_category"]
        entries.append({
            "id": reference["id"],
            "expected_category": expected,
            "predicted_category": prediction,
            "exact_match": prediction == expected,
            "false_reassurance": expected in NONPERSISTENT_CATEGORIES and prediction == "ordinary_condition_persistent",
        })
    return {
        "schema_version": 1,
        "kind": "plausibility_benchmark_report",
        "manifest": str(manifest_path),
        "accepted_entries": sum(item.get("evidence_status") == "accepted" for item in manifest.get("entries", [])),
        "evaluated_entries": len(entries),
        "missing_entries": missing,
        "entries": entries,
        "aggregate": {
            "exact_category_accuracy": sum(item["exact_match"] for item in entries) / len(entries) if entries else None,
            "false_reassurance_count": sum(item["false_reassurance"] for item in entries),
            "false_reassurance_rate": sum(item["false_reassurance"] for item in entries) / len(entries) if entries else None,
        },
        "limitations": "Only accepted entries with explicitly matched conditions are scored. Draft entries cannot calibrate or validate a persistence threshold.",
    }


def write_plausibility_benchmark(path: Path, report: dict) -> Path:
    path = Path(path)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return path
