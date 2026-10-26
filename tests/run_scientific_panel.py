"""Score stored panel predictions and export reference-free prepared requests.

This bounded offline driver launches no ORCA/RMG/Arkane jobs. Reuse the P1
preregistration and retained outcome with --preregistration and --p1-result;
that adapter admits only the existing verified lifetime field and observed
driver wall time. It never derives a verified rate from a screening lifetime.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from run_p1_preflight import THREAD_LIMITS, digest, save


def _named_paths(values):
    paths = {}
    for value in values:
        name, delimiter, path = value.partition("=")
        if not delimiter or not name.strip() or not path or name in paths:
            raise ValueError("Expected unique METHOD=PATH arguments")
        paths[name] = Path(path)
    return paths


def run(manifest_path, output, *, result_directories, preregistration_paths=(), p1_results=None,
        rate_preregistration_paths=()):
    import storca.scientific_panel as evaluator
    from storca.scientific_panel import (
        canonical_hash, evaluate_scientific_panel, load_scientific_panel_manifest,
        prediction_from_p1_result, prediction_request,
    )

    os.environ.update(THREAD_LIMITS)
    if not result_directories:
        raise ValueError("At least one stored prediction method is required")
    p1_results = p1_results or {}
    if set(p1_results) - set(result_directories):
        raise ValueError("P1 outcomes must name a declared prediction method")
    if p1_results and len(preregistration_paths) != 1:
        raise ValueError("P1 outcome adaptation needs exactly one P1 preregistration")
    manifest_path, output = Path(manifest_path).resolve(), Path(output).resolve()
    panel = load_scientific_panel_manifest(manifest_path, preregistration_paths=preregistration_paths,
                                         rate_preregistration_paths=rate_preregistration_paths)
    requests = [prediction_request(case) for case in panel["cases"] if case.get("scoring") or case.get("rate_only_scoring")]
    reports = {}
    for name, directory in result_directories.items():
        retained = {}
        if name in p1_results:
            prediction = prediction_from_p1_result(preregistration_paths[0], p1_results[name])
            retained[prediction["case_id"]] = prediction
        reports[name] = evaluate_scientific_panel(
            manifest_path, directory, preregistration_paths=preregistration_paths,
            retained_predictions=retained,
            rate_preregistration_paths=rate_preregistration_paths,
        )
    comparisons = []
    for request in requests:
        case_id = request["case_id"]
        method_scores = {name: next(score for score in report["cases"] + report["rate_only"]["cases"] if score["case_id"] == case_id)
                         for name, report in reports.items()}
        fallback_scope = "rate_only" if request["kind"] == "storca_scientific_panel_rate_request" else "full_case"
        scopes = []
        for score in method_scores.values():
            scope = score.get("observable_scope", fallback_scope)
            if not isinstance(scope, str) or not scope.strip():
                raise ValueError(f"Invalid observable scope for {case_id}")
            scopes.append(scope)
        if len(set(scopes)) != 1:
            raise ValueError(f"Conflicting method observable scopes for {case_id}")
        budgets = [score["execution_budget"] for score in method_scores.values()]
        matched = (len(budgets) >= 2 and all(isinstance(item, dict) and item for item in budgets)
                   and all(isinstance(item.get("overall_timeout_seconds"), (int, float))
                           and not isinstance(item["overall_timeout_seconds"], bool)
                           and 0 < item["overall_timeout_seconds"] < float("inf")
                           and isinstance(item.get("maximum_cores"), int)
                           and not isinstance(item["maximum_cores"], bool) and item["maximum_cores"] > 0
                           for item in budgets)
                   and len({canonical_hash(item) for item in budgets}) == 1)
        comparisons.append({"case_id": case_id, "scope": scopes[0],
                            "matched_resource_ceilings": bool(matched),
                            "method_outcomes": {name: score["outcome"] for name, score in method_scores.items()}})
    output.mkdir(parents=True, exist_ok=False)
    # Preserve hashes, not reference-only source material, in calculation exports.
    save(output, "prediction-requests.json", {
        "schema_version": 1, "kind": "storca_scientific_panel_prepared_requests",
        "execution_admitted": False, "requests": requests,
        "request_hashes": {item["case_id"]: canonical_hash(item) for item in requests},
        "missing_request_case_count": len(panel["cases"]) - len(requests),
    })
    report = {
        "schema_version": 1, "kind": "storca_scientific_panel_batch_report",
        "support_claim_allowed": False, "calibration_claim_allowed": False,
        "chemistry_jobs_executed_by_evaluator": 0, "methods": reports,
        "comparisons": comparisons,
        "provenance": {"panel_sha256": digest(manifest_path.read_bytes()),
                       "driver_sha256": digest(Path(__file__).read_bytes()),
                       "evaluator_sha256": digest(Path(evaluator.__file__).read_bytes()),
                       "result_directories": {name: str(Path(path).resolve()) for name, path in result_directories.items()},
                       "p1_outcomes": {name: str(Path(path).resolve()) for name, path in p1_results.items()},
                       "preregistration_sha256": [digest(Path(path).read_bytes()) for path in preregistration_paths]},
    }
    report["provenance"]["rate_preregistration_sha256"] = [digest(Path(path).read_bytes()) for path in rate_preregistration_paths]
    save(output, "report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--results", action="append", required=True, metavar="METHOD=DIR")
    parser.add_argument("--preregistration", action="append", type=Path, default=[])
    parser.add_argument("--rate-preregistration", action="append", type=Path, default=[],
                        help="Internal offline rate-only development protocol; grants no chemistry execution")
    parser.add_argument("--p1-result", action="append", default=[], metavar="METHOD=RESULT_JSON")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.panel, args.output_dir, result_directories=_named_paths(args.results),
                 preregistration_paths=args.preregistration, p1_results=_named_paths(args.p1_result),
                 rate_preregistration_paths=args.rate_preregistration)
    print(json.dumps({"output_dir": str(args.output_dir.resolve()),
                      "method_count": len(report["methods"]), "support_claim_allowed": False}))


if __name__ == "__main__":
    main()
