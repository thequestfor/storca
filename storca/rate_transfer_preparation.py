"""Preparation-only, known-reaction transfer protocol; no native dispatcher.

Unlike an autonomous fate request, this protocol intentionally supplies products.
Its separate scoring dossier must never travel to a prediction backend. The
nominal temperatures are curation targets, not manufactured measurement rows.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile

from .scientific_panel import canonical_hash


CASE_IDS = ("h_ch4", "oh_h2")
KIND = "storca_rate_transfer_preparation_request"
FLAGS = ("rate_claim_allowed", "lifetime_claim_allowed", "network_reinjection_allowed",
         "production_defaults_changed", "calibrated_uncertainty_claim_allowed",
         "angular_partition_admitted")
PROPOSED_BUDGET = {
    "authorization": "proposed_only_not_authorized", "maximum_cores": 4,
    "timeout_seconds_per_job": 21600, "maxcore_MB_per_rank": 1000,
    "minimum_free_disk_bytes": 4 * 1024**3,
    "per_case_attempt_caps": {
        case: {"optimization_or_path_jobs": 6, "NumFreq": 4, "IRC": 2, "single_point": 5}
        for case in CASE_IDS
    },
    "maximum_total_native_wall_seconds_per_case": 86400,
    "maximum_total_native_core_seconds_per_case": 345600,
    "work_unit_policy": "Every top-level engine dispatch consumes a category attempt; NEB/path/recovery uses optimization_or_path_jobs; separately dispatched image screening uses single_point. Internal workers/displacements consume the owning job's wall/core budget. Four NumFreq caps cover three molecular species plus TS; the monatomic H asymptote has no vibration.",
    "failed_attempts_count": True, "reused_identical_jobs_are_not_redispatched": True,
}
SCORING_RULES = {
    "maximum_geometric_mean_multiplicative_error_per_reaction": 3.0,
    "maximum_pointwise_multiplicative_error": 10.0,
    "errors": ["signed_log10", "absolute_log10", "temperature_dependence_log10_error_change"],
    "aggregation": "one_summary_per_reaction_all_frozen_cases_in_denominator",
    "meaning": "engineering_usefulness_targets_not_calibrated_error_or_support",
    "author_fit_values_are_measured_rows": False,
}


def binding(path: Path) -> dict:
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _reopen(record: dict) -> Path:
    if (not isinstance(record, dict) or set(record) != {"path", "sha256"}
            or not isinstance(record["path"], str) or not Path(record["path"]).is_absolute()
            or not re.fullmatch(r"[0-9a-f]{64}", str(record["sha256"]))):
        raise ValueError("Exact absolute artifact binding required")
    path = Path(record["path"])
    if binding(path) != record:
        raise ValueError("Bound transfer source changed: " + str(path))
    return path


def _json(path: Path) -> dict:
    def reject_constant(value):
        raise ValueError("Nonfinite JSON constant: " + value)
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate JSON key: " + key)
            result[key] = value
        return result
    value = json.loads(path.read_text(), parse_constant=reject_constant,
                       object_pairs_hook=unique_object)
    if not isinstance(value, dict):
        raise ValueError("JSON object required")
    return value


def _candidate(record: dict) -> dict:
    candidate = _json(_reopen(record))
    if (type(candidate.get("schema_version")) is not int or candidate["schema_version"] != 1
            or candidate.get("kind") != "OH_CH4_frozen_conditional_rate_candidate"
            or candidate.get("status") != "diagnostic_not_admitted"
            or any(candidate.get(flag) is not False for flag in FLAGS)
            or candidate.get("physical_rate_accuracy_claim_allowed", False) is not False
            or candidate.get("recipe_candidate_has_unresolved_assumptions") is not True):
        raise ValueError("An explicitly frozen conditional, non-admitted candidate is required")
    for key in ("campaign", "recipe", "electronic_state_replay"):
        _reopen(candidate[key])
    if not candidate.get("source_bindings") or not candidate.get("source_snapshots"):
        raise ValueError("Frozen candidate native evidence and implementation snapshots required")
    for source in candidate["source_bindings"]:
        _reopen(source)
    for pair in candidate["source_snapshots"]:
        _reopen(pair["snapshot"])
        if pair["original"]["sha256"] != pair["snapshot"]["sha256"]:
            raise ValueError("Candidate implementation snapshot differs from original revision")
    return candidate


def prediction_request(case_id: str, candidate: dict) -> dict:
    """Construct an exact, answer-free request, with no opaque metadata fields."""
    _candidate(candidate)
    return _prediction_inputs(case_id, candidate)


def _prediction_inputs(case_id: str, candidate: dict) -> dict:
    """Build fixed inputs after the caller has validated the shared candidate."""
    if case_id == "h_ch4":
        species = [("[H:1]", 2), ("[C:2]([H:3])([H:4])([H:5])[H:6]", 1),
                   ("[H:1][H:6]", 1), ("[C:2]([H:3])([H:4])[H:5]", 2)]
        temperatures = [1000.0, 1200.0, 1500.0]
    elif case_id == "oh_h2":
        species = [("[O:1][H:2]", 2), ("[H:3][H:4]", 1),
                   ("[O:1]([H:2])[H:3]", 1), ("[H:4]", 2)]
        temperatures = [298.0, 375.0, 450.0]
    else:
        raise ValueError("Unsupported known-reaction transfer case")
    records = [{"mapped_smiles": smiles, "charge": 0, "multiplicity": spin,
                "stoichiometric_coefficient": 1} for smiles, spin in species]
    return {
        "schema_version": 1, "kind": "storca_known_reaction_rate_request",
        "protocol": "known_reaction_elementary_rate", "case_id": case_id,
        "reactants": records[:2], "products": records[2:],
        "total_charge": 0, "surface_multiplicity": 2,
        "environment": {"phase": "gas", "surface_model": "single_adiabatic_surface",
                        "pressure_bar": None, "pressure_regime": "unresolved",
                        "bath_composition": None, "light_condition": "thermal_dark"},
        "temperature_grid_K": temperatures,
        "temperature_grid_status": "proposed_selection_targets_not_measurements",
        "requested_observable": "forward_elementary_bimolecular_rate_coefficient",
        "rate_units": "cm^3 molecule^-1 s^-1",
        "method_recipe": {"candidate": candidate, "use": "conditional_diagnostic_only",
                          "transfer_policy": "Transfer the procedure, recompute new species and TS energies, states and channel/symmetry counts; never copy OH/CH4 numeric counts."},
        "execution_admitted": False,
    }


def _inputs(request_path: Path) -> tuple[dict, dict, dict]:
    request = _json(request_path)
    if (set(request) != {"schema_version", "kind", "candidate_recipe", "reference_curation"}
            or type(request["schema_version"]) is not int or request["schema_version"] != 1
            or request["kind"] != KIND):
        raise ValueError("Unsupported preparation protocol or non-allowlisted request fields")
    candidate = _candidate(request["candidate_recipe"])
    curation = _json(_reopen(request["reference_curation"]))
    if (type(curation.get("schema_version")) is not int or curation["schema_version"] != 1
            or curation.get("protocol") != "known_reaction_elementary_rate"
            or [case.get("case_id") for case in curation.get("cases", [])] != list(CASE_IDS)):
        raise ValueError("Ordered known-reaction reference dossier required")
    return request, candidate, curation


def _bundle(request_path: Path) -> tuple[dict, dict, dict]:
    request, _, curation = _inputs(request_path)
    predictions = {"schema_version": 1, "kind": "storca_known_reaction_transfer_inputs",
                   "cases": [_prediction_inputs(case, request["candidate_recipe"]) for case in CASE_IDS]}
    scoring = {"schema_version": 1, "kind": "storca_rate_transfer_scoring_dossier",
               "reference_curation": request["reference_curation"],
               "cases": curation["cases"], "rules": SCORING_RULES,
               "reference_data_is_prediction_input": False}
    gaps = ["separate_compute_budget_authorization_required",
            "known_reaction_execution_adapter_not_wired_in_this_increment",
            "conditional_candidate_physical_assumptions_not_admitted",
            "primary_audited_exact_temperature_and_environment_lock_required",
            "h_ch4:reaction_degeneracy_TS_symmetry_and_channel_counting_unresolved",
            "oh_h2:OH_state_TS_counting_and_H2_ortho_para_rotational_population_unresolved",
            "transfer_new_species_TS_counts_must_not_copy_OH_CH4_values"]
    for case in curation["cases"]:
        for key in ("measured_rows_audited", "conditions_audited", "uncertainty_audited",
                    "corrections_audited", "ready_for_scoring"):
            if case.get(key) is not True:
                gaps.append(case["case_id"] + ":" + key + "_pending")
        if not case.get("measured_rows"):
            gaps.append(case["case_id"] + ":actual_measured_rows_missing")
    manifest = {
        "schema_version": 1, "kind": "storca_rate_transfer_preparation",
        "status": "prepared_not_authorized", "execution_ready": False,
        "execution_admitted": False, "native_jobs_launched": 0,
        "scientific_pass": False, "support_claim_allowed": False,
        "calibrated_uncertainty_claim_allowed": False,
        "partition": "exposed_development_transfer_not_blind_holdout",
        "request": binding(request_path), "implementation": binding(Path(__file__)),
        "candidate_recipe": request["candidate_recipe"],
        "reference_curation": request["reference_curation"],
        "prediction_inputs_sha256": canonical_hash(predictions),
        "scoring_dossier_sha256": canonical_hash(scoring),
        "proposed_budget": PROPOSED_BUDGET, "readiness_gaps": gaps,
        "next_step": "Audit primary rows/corrections/conditions and freeze revised condition lock. Before separately authorized native execution, wire this known-reaction request to existing route verification and source-bound composite/state evaluators; this preparation CLI cannot launch a benchmark.",
    }
    return manifest, predictions, scoring


def _publish(path: Path, value: dict) -> None:
    content = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
    if path.exists():
        if path.read_bytes() != content:
            raise ValueError("Refusing to overwrite changed frozen transfer artifact: " + str(path))
        return
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".transfer-", delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
            try:
                os.link(temporary, path)
            except FileExistsError:
                if path.read_bytes() != content:
                    raise ValueError("Concurrent transfer artifact differs")
        finally:
            temporary.unlink()


def prepare_transfer_benchmark(request_path: Path, output_dir: Path) -> dict:
    """Freeze three immutable preparation files; identical replay launches zero jobs."""
    manifest, predictions, scoring = _bundle(Path(request_path))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, payload in (("prediction-requests.json", predictions), ("scoring-dossier.json", scoring),
                          ("manifest.json", manifest)):
        _publish(output_dir / name, payload)
    return validate_preparation(output_dir)


def validate_preparation(output_dir: Path) -> dict:
    output_dir = Path(output_dir)
    supplied = _json(output_dir / "manifest.json")
    request_path = _reopen(supplied["request"])
    expected = _bundle(request_path)
    for name, payload in zip(("manifest.json", "prediction-requests.json", "scoring-dossier.json"), expected):
        if canonical_hash(_json(output_dir / name)) != canonical_hash(payload):
            raise ValueError("Frozen transfer preparation changed: " + name)
    return supplied


def score_transfer_outcomes(scoring: dict, outcomes: dict) -> dict:
    """Descriptive preregistered arithmetic, never a physical admission decision.

    Every frozen case stays in the denominator, including absent or unusable
    results. Audited measurement rows, not author-fit evaluations, are required.
    """
    if (scoring.get("kind") != "storca_rate_transfer_scoring_dossier"
            or scoring.get("rules") != SCORING_RULES
            or [case.get("case_id") for case in scoring.get("cases", [])] != list(CASE_IDS)
            or set(outcomes) - set(CASE_IDS)):
        raise ValueError("Unsupported frozen transfer scoring protocol")
    reports = []
    for case in scoring["cases"]:
        result = outcomes.get(case["case_id"], {})
        disposition = result.get("disposition", "missing")
        if disposition not in {"prediction", "missing", "unresolved", "crashed", "budget_stopped",
                               "outside_applicability"}:
            disposition = "invalid"
        report = {"case_id": case["case_id"], "outcome": disposition,
                  "engineering_target_met": False, "scored_points": 0}
        references = case.get("measured_rows", [])
        reference_ready = bool(references) and all(case.get(field) is True for field in (
            "measured_rows_audited", "conditions_audited", "uncertainty_audited",
            "corrections_audited", "ready_for_scoring"))
        if disposition == "prediction" and not reference_ready:
            report["outcome"] = "reference_pending"
        elif disposition == "prediction":
            try:
                rows = result["rate_points"]
                if len(rows) != len(references):
                    raise ValueError("Incomplete prediction grid")
                errors = []
                for ref, row in zip(references, rows):
                    temperatures = (row["temperature_K"], ref["temperature_K"])
                    if (any(type(t) not in (int, float) or not math.isfinite(t) or t <= 0 for t in temperatures)
                            or row["temperature_K"] != ref["temperature_K"]
                            or row["rate_units"] != ref["rate_units"]
                            or row["rate_units"] != "cm^3 molecule^-1 s^-1"):
                        raise ValueError("Rate conditions or units differ")
                    conditions = ref["conditions"]
                    if (not isinstance(conditions, dict) or conditions.get("phase") != "gas"
                            or type(conditions.get("pressure_bar")) not in (int, float)
                            or not math.isfinite(conditions["pressure_bar"]) or conditions["pressure_bar"] <= 0
                            or not conditions.get("bath_composition") or not conditions.get("state_preparation")
                            or canonical_hash(conditions) != canonical_hash(row["conditions"])):
                        raise ValueError("Full frozen pressure/bath/state preparation identity required")
                    values = (row["rate_constant"], ref["rate_constant"])
                    if any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0 for v in values):
                        raise ValueError("Positive finite rate constants required")
                    errors.append(math.log10(values[0]) - math.log10(values[1]))
                geometric = 10 ** (sum(map(abs, errors)) / len(errors))
                worst = 10 ** max(map(abs, errors))
                report.update(scored_points=len(errors), signed_log10_errors=errors,
                              absolute_log10_errors=list(map(abs, errors)),
                              geometric_mean_multiplicative_error=geometric,
                              worst_multiplicative_error=worst,
                              temperature_dependence_log10_error_change=(errors[-1] - errors[0] if len(errors) > 1 else None),
                              engineering_target_met=geometric <= 3 and worst <= 10)
            except (KeyError, TypeError, ValueError, OverflowError):
                report["outcome"] = "invalid"
        reports.append(report)
    return {"cases": reports, "declared_case_count": len(CASE_IDS),
            "quantitatively_scored_case_count": sum(row["scored_points"] > 0 for row in reports),
            "engineering_target_met": all(row["engineering_target_met"] for row in reports),
            "scientific_pass": False, "support_claim_allowed": False,
            "calibrated_uncertainty_claim_allowed": False}
