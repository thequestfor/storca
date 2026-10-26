"""Scientific benchmark and conformance gate for specialist kinetic providers.

The benchmark deliberately separates provider conformance from scientific
accuracy.  A provider can produce well-formed evidence and still fail the rate
benchmark; a development reference set can be evaluated without authorizing a
supported-lane claim.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
from statistics import median
from typing import Any, Mapping


SPECIALIST_RATE_BENCHMARK_SCHEMA_VERSION = 1
PREDICTION_SCHEMA_VERSION = 1


def _canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _finite_positive(value: object, *, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be a finite positive number") from error
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{label} must be a finite positive number")
    return number


def _validate_source(source: Mapping[str, Any], *, case_id: str) -> None:
    required = ("citation", "url", "source_kind", "reference_scope")
    missing = [key for key in required if not str(source.get(key) or "").strip()]
    if missing:
        raise ValueError(
            f"Benchmark case {case_id} has incomplete reference source: {', '.join(missing)}"
        )


def load_specialist_rate_benchmark_manifest(path: Path) -> dict[str, Any]:
    """Load a condition-complete, provenance-backed specialist-rate benchmark."""
    payload = json.loads(Path(path).read_text())
    if payload.get("schema_version") != SPECIALIST_RATE_BENCHMARK_SCHEMA_VERSION:
        raise ValueError("Unsupported specialist-rate benchmark manifest schema")
    if payload.get("kind") != "storca_specialist_rate_benchmark_manifest":
        raise ValueError("Invalid specialist-rate benchmark manifest kind")
    gate = payload.get("acceptance_gate") or {}
    for field in (
        "maximum_median_absolute_log10_error",
        "maximum_absolute_log10_error",
        "minimum_uncertainty_coverage_fraction",
        "maximum_false_applicability_count",
    ):
        if field not in gate:
            raise ValueError(f"Specialist-rate benchmark acceptance gate needs {field}")
    cases = payload.get("cases") or ()
    if not cases:
        raise ValueError("Specialist-rate benchmark manifest has no cases")
    identifiers: set[str] = set()
    partitions: set[str] = set()
    for case in cases:
        case_id = str(case.get("case_id") or "").strip()
        if not case_id or case_id in identifiers:
            raise ValueError("Specialist-rate benchmark case IDs must be present and unique")
        identifiers.add(case_id)
        partition = str(case.get("partition") or "")
        if partition not in {"development", "holdout", "challenge"}:
            raise ValueError(f"Benchmark case {case_id} has invalid partition")
        partitions.add(partition)
        if not str(case.get("reaction_equation") or "").strip():
            raise ValueError(f"Benchmark case {case_id} needs a reaction equation")
        expected = str(case.get("expected_provider_disposition") or "")
        if expected not in {"verified", "unresolved"}:
            raise ValueError(f"Benchmark case {case_id} has invalid expected disposition")
        observable_kind = str(case.get("observable_kind") or "rate_points")
        if observable_kind not in {"rate_points", "falloff_curve"}:
            raise ValueError(f"Benchmark case {case_id} has invalid observable kind")
        sources = case.get("reference_sources") or ()
        if not sources:
            raise ValueError(f"Benchmark case {case_id} needs reference provenance")
        for source in sources:
            _validate_source(source, case_id=case_id)
        points = case.get("reference_points") or ()
        if expected == "verified" and not points:
            raise ValueError(f"Verified benchmark case {case_id} needs reference points")
        point_ids: set[str] = set()
        for point in points:
            point_id = str(point.get("point_id") or "").strip()
            if not point_id or point_id in point_ids:
                raise ValueError(f"Benchmark case {case_id} point IDs must be present and unique")
            point_ids.add(point_id)
            _finite_positive(point.get("temperature_K"), label=f"{case_id}/{point_id} temperature")
            _finite_positive(point.get("reference_rate"), label=f"{case_id}/{point_id} reference rate")
            _finite_positive(
                point.get("reference_uncertainty_factor"),
                label=f"{case_id}/{point_id} reference uncertainty factor",
            )
            if float(point["reference_uncertainty_factor"]) < 1.0:
                raise ValueError(f"{case_id}/{point_id} uncertainty factor must be at least one")
            if not str(point.get("rate_units") or "").strip():
                raise ValueError(f"{case_id}/{point_id} needs explicit rate units")
            if str(point.get("pressure_regime") or "") not in {
                "low_pressure_limit", "high_pressure_limit", "finite_pressure",
            }:
                raise ValueError(f"{case_id}/{point_id} has invalid pressure regime")
            if point.get("pressure_regime") == "finite_pressure":
                _finite_positive(
                    point.get("pressure_bar"), label=f"{case_id}/{point_id} pressure",
                )
        if observable_kind == "falloff_curve" and expected == "verified":
            if any(point.get("pressure_regime") != "finite_pressure" for point in points):
                raise ValueError(f"Falloff case {case_id} may contain only finite-pressure points")
            pressures = sorted({float(point["pressure_bar"]) for point in points})
            if len(pressures) < 3:
                raise ValueError(f"Falloff case {case_id} needs at least three unique pressures")

    pending_ids: set[str] = set()
    for pending in payload.get("pending_cases") or ():
        case_id = str(pending.get("case_id") or "").strip()
        if not case_id or case_id in identifiers or case_id in pending_ids:
            raise ValueError("Pending benchmark case IDs must be present and globally unique")
        pending_ids.add(case_id)
        if pending.get("curation_status") != "reference_data_pending":
            raise ValueError(f"Pending benchmark case {case_id} needs an explicit pending status")
        if not pending.get("missing_reference_data"):
            raise ValueError(f"Pending benchmark case {case_id} must name its missing data")
        sources = pending.get("reference_sources") or ()
        if not sources:
            raise ValueError(f"Pending benchmark case {case_id} needs candidate sources")
        for source in sources:
            _validate_source(source, case_id=case_id)

    if payload.get("status") == "locked" or payload.get("holdout_locked"):
        holdout_cases = [case for case in cases if case.get("partition") == "holdout"]
        lock = payload.get("holdout_lock") or {}
        for field in (
            "locked_at_utc", "protocol_sha256", "holdout_content_sha256", "custodian",
        ):
            if not str(lock.get(field) or "").strip():
                raise ValueError(f"Locked benchmark needs holdout_lock.{field}")
        for field in ("protocol_sha256", "holdout_content_sha256"):
            if not re.fullmatch(r"sha256:[0-9a-f]{64}", str(lock.get(field) or "")):
                raise ValueError(f"Locked benchmark has invalid holdout_lock.{field}")
        if lock.get("provider_access") != "withheld_until_predictions_frozen":
            raise ValueError("Locked holdout must remain withheld until predictions are frozen")
        if lock.get("holdout_content_sha256") != _canonical_hash(holdout_cases):
            raise ValueError("Locked holdout content hash does not match the manifest")
    payload["partitions_present"] = sorted(partitions)
    return payload


def _artifact_conformance(result: Mapping[str, Any], root: Path) -> list[str]:
    blockers: list[str] = []
    artifacts = result.get("evidence_artifacts") or ()
    if not artifacts:
        return ["prediction_evidence_artifacts_missing"]
    for artifact in artifacts:
        raw_path = Path(str((artifact or {}).get("path") or ""))
        path = raw_path if raw_path.is_absolute() else root / raw_path
        expected = str((artifact or {}).get("sha256") or "")
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        if not actual or actual != expected:
            blockers.append("prediction_evidence_artifact_missing_or_mutated")
    return blockers


def _score_case(
    case: Mapping[str, Any], result: Mapping[str, Any], root: Path, *, holdout_locked: bool,
) -> dict[str, Any]:
    case_id = str(case["case_id"])
    blockers: list[str] = []
    if result.get("schema_version") != PREDICTION_SCHEMA_VERSION:
        blockers.append("prediction_schema_unsupported")
    if result.get("kind") != "storca_specialist_rate_benchmark_prediction":
        blockers.append("prediction_kind_invalid")
    if result.get("case_id") != case_id:
        blockers.append("prediction_case_mismatch")
    provider = result.get("provider") or {}
    if not all(str(provider.get(key) or "").strip() for key in ("name", "version")):
        blockers.append("prediction_provider_identity_incomplete")

    expected = str(case["expected_provider_disposition"])
    disposition = str(result.get("disposition") or "")
    false_applicability = expected == "unresolved" and disposition == "verified"
    missed_applicability = expected == "verified" and disposition != "verified"
    if disposition not in {"verified", "unresolved", "outside_applicability"}:
        blockers.append("prediction_disposition_invalid")
    if false_applicability:
        blockers.append("false_method_applicability")
    if missed_applicability:
        blockers.append("required_case_not_verified")
    if holdout_locked and case.get("partition") == "holdout":
        blinding = result.get("holdout_blinding") or {}
        artifact_sha = str(blinding.get("provider_artifact_sha256") or "")
        if not (
            blinding.get("reference_access_before_freeze") is False
            and blinding.get("prediction_frozen_before_unblinding") is True
            and len(artifact_sha) == 64
            and all(character in "0123456789abcdef" for character in artifact_sha)
        ):
            blockers.append("holdout_blinding_evidence_missing_or_invalid")

    applicability = result.get("applicability") or {}
    if not applicability.get("diagnostics"):
        blockers.append("prediction_applicability_diagnostics_missing")
    if disposition == "verified" and applicability.get("status") != "applicable":
        blockers.append("prediction_applicability_not_established")
    if disposition != "verified" and not str(result.get("unresolved_reason") or "").strip():
        blockers.append("prediction_unresolved_reason_missing")
    required_reason = str(case.get("required_unresolved_reason") or "").strip()
    if required_reason and result.get("unresolved_reason") != required_reason:
        blockers.append("prediction_unresolved_reason_mismatch")

    point_scores: list[dict[str, Any]] = []
    if disposition == "verified":
        blockers.extend(_artifact_conformance(result, root))
        predicted = {
            str(item.get("point_id") or ""): item
            for item in result.get("predictions") or ()
        }
        expected_ids = {str(item["point_id"]) for item in case.get("reference_points") or ()}
        if set(predicted) != expected_ids:
            blockers.append("prediction_grid_mismatch")
        for reference in case.get("reference_points") or ():
            point_id = str(reference["point_id"])
            item = predicted.get(point_id)
            if item is None:
                continue
            if item.get("rate_units") != reference.get("rate_units"):
                blockers.append("prediction_rate_units_mismatch")
                continue
            try:
                predicted_temperature = float(item.get("temperature_K"))
            except (TypeError, ValueError):
                predicted_temperature = math.nan
            if not math.isfinite(predicted_temperature) or not math.isclose(
                predicted_temperature, float(reference["temperature_K"]),
                rel_tol=0.0, abs_tol=1.0e-9,
            ):
                blockers.append("prediction_temperature_mismatch")
                continue
            if item.get("pressure_regime") != reference.get("pressure_regime"):
                blockers.append("prediction_pressure_regime_mismatch")
                continue
            if reference.get("pressure_regime") == "finite_pressure":
                try:
                    predicted_pressure = float(item.get("pressure_bar"))
                except (TypeError, ValueError):
                    predicted_pressure = math.nan
                if not math.isfinite(predicted_pressure) or not math.isclose(
                    predicted_pressure, float(reference["pressure_bar"]),
                    rel_tol=0.0, abs_tol=1.0e-12,
                ):
                    blockers.append("prediction_pressure_mismatch")
                    continue
            try:
                predicted_rate = _finite_positive(
                    item.get("predicted_rate"), label=f"{case_id}/{point_id} prediction"
                )
                lower = _finite_positive(
                    item.get("lower_rate"), label=f"{case_id}/{point_id} lower rate"
                )
                upper = _finite_positive(
                    item.get("upper_rate"), label=f"{case_id}/{point_id} upper rate"
                )
            except ValueError:
                blockers.append("prediction_rate_or_interval_invalid")
                continue
            if not lower <= predicted_rate <= upper:
                blockers.append("prediction_interval_does_not_bound_central_rate")
                continue
            reference_rate = float(reference["reference_rate"])
            factor = float(reference["reference_uncertainty_factor"])
            absolute_log_error = abs(math.log10(predicted_rate / reference_rate))
            point_scores.append({
                "point_id": point_id,
                "temperature_K": float(reference["temperature_K"]),
                "pressure_bar": reference.get("pressure_bar"),
                "absolute_log10_error": absolute_log_error,
                "within_reference_uncertainty": (
                    reference_rate / factor <= predicted_rate <= reference_rate * factor
                ),
                "provider_interval_covers_reference": lower <= reference_rate <= upper,
            })

    curve_metrics = None
    if str(case.get("observable_kind") or "rate_points") == "falloff_curve" and point_scores:
        errors = [float(point["absolute_log10_error"]) for point in point_scores]
        pressures = [float(point["pressure_bar"]) for point in point_scores]
        curve_metrics = {
            "point_count": len(point_scores),
            "pressure_span_decades": math.log10(max(pressures) / min(pressures)),
            "rmse_log10_rate": math.sqrt(sum(error * error for error in errors) / len(errors)),
            "maximum_absolute_log10_error": max(errors),
        }
    blockers = list(dict.fromkeys(blockers))
    return {
        "case_id": case_id,
        "partition": case["partition"],
        "expected_provider_disposition": expected,
        "provider_disposition": disposition or None,
        "conformance_status": "passed" if not blockers else "failed",
        "blocking_reasons": blockers,
        "false_applicability": false_applicability,
        "missed_applicability": missed_applicability,
        "point_scores": point_scores,
        "curve_metrics": curve_metrics,
    }


def evaluate_specialist_rate_benchmark(
    manifest_path: Path,
    results_dir: Path,
) -> dict[str, Any]:
    """Evaluate retained provider predictions without silently granting support."""
    manifest = load_specialist_rate_benchmark_manifest(manifest_path)
    root = Path(results_dir)
    scores: list[dict[str, Any]] = []
    for case in manifest["cases"]:
        case_id = str(case["case_id"])
        result_path = root / f"{case_id}.json"
        if not result_path.is_file():
            scores.append({
                "case_id": case_id,
                "partition": case["partition"],
                "expected_provider_disposition": case["expected_provider_disposition"],
                "provider_disposition": None,
                "conformance_status": "failed",
                "blocking_reasons": ["benchmark_prediction_missing"],
                "false_applicability": False,
                "missed_applicability": case["expected_provider_disposition"] == "verified",
                "point_scores": [],
            })
            continue
        try:
            result = json.loads(result_path.read_text())
        except (OSError, json.JSONDecodeError):
            result = {}
        scores.append(_score_case(
            case, result, result_path.parent,
            holdout_locked=bool(manifest.get("holdout_locked")),
        ))

    points = [point for case in scores for point in case["point_scores"]]
    absolute_errors = [point["absolute_log10_error"] for point in points]
    gate = manifest.get("acceptance_gate") or {}
    complete = not any(case["blocking_reasons"] for case in scores)
    numerical_pass = bool(points) and (
        median(absolute_errors) <= float(gate.get("maximum_median_absolute_log10_error", 0.0))
        and max(absolute_errors) <= float(gate.get("maximum_absolute_log10_error", 0.0))
        and sum(point["provider_interval_covers_reference"] for point in points) / len(points)
        >= float(gate.get("minimum_uncertainty_coverage_fraction", 1.0))
    )
    curves = [case["curve_metrics"] for case in scores if case.get("curve_metrics")]
    minimum_curve_count = int(gate.get("minimum_falloff_curve_count", 0))
    minimum_pressure_span = float(gate.get("minimum_falloff_pressure_span_decades", 0.0))
    maximum_curve_rmse = float(gate.get("maximum_falloff_curve_rmse_log10", math.inf))
    curve_pass = (
        len(curves) >= minimum_curve_count
        and all(curve["pressure_span_decades"] >= minimum_pressure_span for curve in curves)
        and all(curve["rmse_log10_rate"] <= maximum_curve_rmse for curve in curves)
    )
    locked = manifest.get("status") == "locked" and bool(manifest.get("holdout_locked"))
    required_partitions = set(manifest.get("required_partitions") or ())
    partitions_present = {str(case["partition"]) for case in manifest["cases"]}
    partition_complete = required_partitions <= partitions_present
    allowed_false_applicability = int(gate.get("maximum_false_applicability_count", 0))
    false_applicability_count = sum(case["false_applicability"] for case in scores)
    curation_complete = not bool(manifest.get("pending_cases"))
    support_allowed = (
        complete and numerical_pass and curve_pass and locked and partition_complete
        and curation_complete
        and false_applicability_count <= allowed_false_applicability
    )
    if not locked:
        status = "scientific_curation_incomplete"
    elif not partition_complete:
        status = "required_partitions_missing"
    elif not curation_complete:
        status = "reference_data_curation_incomplete"
    elif support_allowed:
        status = "acceptance_gate_passed"
    else:
        status = "acceptance_gate_failed"
    return {
        "schema_version": SPECIALIST_RATE_BENCHMARK_SCHEMA_VERSION,
        "kind": "storca_specialist_rate_benchmark_report",
        "status": status,
        "support_claim_allowed": support_allowed,
        "manifest_status": manifest.get("status"),
        "case_count": len(scores),
        "point_count": len(points),
        "conforming_case_count": sum(case["conformance_status"] == "passed" for case in scores),
        "false_applicability_count": false_applicability_count,
        "missed_applicability_count": sum(case["missed_applicability"] for case in scores),
        "median_absolute_log10_error": median(absolute_errors) if absolute_errors else None,
        "maximum_absolute_log10_error": max(absolute_errors) if absolute_errors else None,
        "uncertainty_coverage_fraction": (
            sum(point["provider_interval_covers_reference"] for point in points) / len(points)
            if points else None
        ),
        "falloff_curve_count": len(curves),
        "falloff_curve_gate_passed": curve_pass,
        "pending_case_count": len(manifest.get("pending_cases") or ()),
        "reference_data_curation_complete": curation_complete,
        "cases": scores,
        "interpretation": (
            "Conformance is not scientific accuracy. A support claim additionally requires "
            "a preregistered locked holdout, complete provenance, numerical accuracy, calibrated "
            "uncertainty coverage, and zero false applicability in the declared challenge set."
        ),
    }


def write_specialist_rate_benchmark(path: Path, report: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(report), indent=2, sort_keys=True) + "\n")
    return path


__all__ = [
    "PREDICTION_SCHEMA_VERSION",
    "SPECIALIST_RATE_BENCHMARK_SCHEMA_VERSION",
    "evaluate_specialist_rate_benchmark",
    "load_specialist_rate_benchmark_manifest",
    "write_specialist_rate_benchmark",
]
