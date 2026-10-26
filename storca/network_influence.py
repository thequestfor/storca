"""Source-bound declared-model scenarios, never physical uncertainty bounds.

The native environment supplies operative rates and stiff integration. This
module owns scenario meaning, evidence retention and next-evidence requests;
it does not launch discovery, electronic structure, or rate construction.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path


POLICY_ID = "retained-network-declared-scenarios-v1"
REPORT_KIND = "expanded_retained_network_assessment_report"
_CLAIMS = {
    "physical_rate_bounds_established": False, "uncertainty_calibrated": False,
    "terminal_robustness_eligible": False, "lifetime_claim_allowed": False,
    "conclusive_observable_allowed": False, "physical_irrelevance_certified": False,
}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _artifact(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        handle.write(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
    return _artifact(path)


def _implementation_sources():
    directory = Path(__file__).parent
    names = ("network_influence.py", "retained_linear_network.py", "bounded_handoff.py",
             "method_applicability.py", "rmg_bridge.py", "native_linear_diagnostics.py",
             "reverse_relevance.py", "reaction_correspondence.py", "rmg_identity_bridge.py",
             "flux_verification.py", "native_rate_sources.py", "retained_rate_audit.py",
             "retained_rate_recipe.py")
    return [_artifact(directory / name) for name in names if (directory / name).is_file()]


def network_influence_policy():
    return {"policy_id": POLICY_ID, "max_species": 32, "max_directed_channels": 128,
            "max_scenarios": 16, "retention_reporting_resolution": 1e-4,
            "numerical_parent_absolute_tolerance": 1e-6,
            "scenario_scope": "declared_counterfactual_models_not_rate_uncertainty"}


def build_declared_network_scenarios(network, state_audit):
    """Challenge identifiable state-changing entry/return/escape subgraphs.

    No uncertain rate multiplier is inferred. Removing a directed return is
    explicitly a counterfactual that does not preserve the original ratio.
    """
    parent, product = network["parent_id"], network["product_id"]
    audits = {item["channel_id"]: item for item in state_audit.get("channels", [])}
    flagged = []
    for channel in network["channels"]:
        record = audits.get(channel["channel_id"], {})
        detail = record.get("endpoint_audit", record.get("audit", record))
        if detail.get("status") == "spin_changing_rate_model_unresolved":
            flagged.append(channel)
    scenarios = [{"scenario_id": "baseline", "kind": "retained_full_model",
                  "multipliers": {}, "detailed_balance_status": "unchanged_declared_model"}]
    if flagged:
        scenarios.append({"scenario_id": "state_changing_subnetwork_removed",
                          "kind": "model_comparison", "multipliers": {
                              item["channel_id"]: 0.0 for item in flagged},
                          "detailed_balance_status": "flagged_channels_removed_not_physical_bounds"})
    intermediates = sorted({next(iter(item["products"])) for item in flagged
                            if item["reactants"] == {parent: 1}
                            and item["products"] != {product: 1}})
    for intermediate in intermediates:
        returns = [item["channel_id"] for item in flagged
                   if item["reactants"] == {intermediate: 1} and item["products"] == {parent: 1}]
        escapes = [item["channel_id"] for item in flagged
                   if item["reactants"] == {intermediate: 1} and item["products"] == {product: 1}]
        for name, channels in (("return", returns), ("escape", escapes)):
            if channels:
                scenarios.append({"scenario_id": f"{intermediate}:{name}_removed",
                                  "kind": "directed_channel_counterfactual",
                                  "intermediate_id": intermediate,
                                  "multipliers": {item: 0.0 for item in channels},
                                  "detailed_balance_status": "original_ratio_not_preserved_not_physical_bounds"})
    if len(scenarios) > network_influence_policy()["max_scenarios"]:
        raise ValueError("declared_network_scenario_limit_exceeded")
    return scenarios


def build_network_scenario_response(network, state_audit, native_result, scenarios):
    """Interpret only complete, numerically checked retained-grid responses."""
    result = {"schema_version": 1, "kind": "declared_network_scenario_response",
              "status": "numerical_unresolved", "policy": network_influence_policy(),
              "network_scientific_hash": network["scientific_hash"],
              "scenario_definitions": copy.deepcopy(scenarios),
              "electronic_state_audit": state_audit, "native_diagnostics": native_result,
              "scenarios": [], "potentially_influential_under_declared_scenarios": None,
              "limitations": ["Counterfactual scenarios are not a physical uncertainty distribution.",
                              "Insensitivity in these scenarios cannot certify an unknown branch irrelevant.",
                              "No completeness, pressure adequacy or electronic accuracy is established."],
              **_CLAIMS}
    if native_result.get("status") != "completed":
        result["blocking_reasons"] = native_result.get("blocking_reasons") or [
            native_result.get("failure_reason") or "native_numerical_checks_unresolved"]
        result["scientific_hash"] = _digest(result)
        return result
    checks = native_result.get("numerical_checks") or {}
    if (native_result.get("kind") != "native_linear_model_diagnostics"
            or native_result.get("condition") != network["condition"]
            or native_result.get("source_artifacts") != network["source_artifacts"]
            or checks.get("status") != "passed" or checks.get("native_baseline_checked") is not True
            or checks.get("source_artifacts_unchanged") is not True):
        raise ValueError("declared_network_native_numerical_or_source_evidence_missing")
    rows = native_result.get("scenarios") or []
    expected_ids = [item["scenario_id"] for item in scenarios]
    ids = [item.get("scenario_id") for item in rows]
    if len(ids) != len(set(ids)) or set(ids) != set(expected_ids):
        raise ValueError("declared_network_scenario_inventory_mismatch")
    by_id = {item["scenario_id"]: item for item in rows}
    baseline = by_id["baseline"]
    base_profile = baseline["target_profile"]
    horizon = network["condition_contract"]["target_duration_seconds"]
    base_times = [item["time_seconds"] for item in base_profile]
    if not base_times or base_times[-1] != horizon or base_times[0] != 0.0 \
            or any(b <= a for a, b in zip(base_times, base_times[1:])):
        raise ValueError("declared_network_exact_horizon_or_time_grid_missing")
    for definition in scenarios:
        row = by_id[definition["scenario_id"]]
        if row.get("multipliers") != definition["multipliers"] \
                or (row.get("numerical_checks") or {}).get("status") != "passed":
            raise ValueError("declared_network_scenario_definition_or_numerical_check_mismatch")
        profile = row["target_profile"]
        if [item["time_seconds"] for item in profile] != base_times:
            raise ValueError("declared_network_scenario_time_grid_mismatch")
        values = [row["horizon_parent_fraction"],
                  *(item["target_fraction_remaining"] for item in profile)]
        if any(isinstance(value, bool) or not isinstance(value, (int, float))
               or not math.isfinite(value) or value < -1e-12 or value > 1 + 1e-12 for value in values):
            raise ValueError("declared_network_parent_fraction_invalid")
        if values[0] != profile[-1]["target_fraction_remaining"]:
            raise ValueError("declared_network_horizon_value_mismatch")
        maximum = max(abs(a["target_fraction_remaining"] - b["target_fraction_remaining"])
                      for a, b in zip(profile, base_profile))
        result["scenarios"].append({**definition, **row,
                                    "change_from_baseline_at_horizon": (
                                        row["horizon_parent_fraction"] - baseline["horizon_parent_fraction"]),
                                    "maximum_retained_grid_change_from_baseline": maximum,
                                    "potentially_influential_under_this_scenario": (
                                        maximum >= result["policy"]["retention_reporting_resolution"])})
    result.update(status="completed_declared_scenarios", blocking_reasons=[],
                  potentially_influential_under_declared_scenarios=any(
                      item["potentially_influential_under_this_scenario"] for item in result["scenarios"]))
    result["scientific_hash"] = _digest(result)
    return result


def next_network_evidence_requests(network, state_audit, response, rate_audit=None):
    flagged = [item for item in state_audit.get("reactions", [])
               if not item.get("ordinary_path_eligible", False)]
    requests = []
    if flagged:
        rate_rows = [] if rate_audit is None else [item for item in rate_audit["channels"]
                                                  if not item["ordinary_path_eligible"]]
        requests.append({"action": "resolve_spin_changing_branch_model_or_supported_relevance_bound",
                         "routes": flagged, "reason": "Retained endpoint-state concern; scenario influence is diagnostic only.",
                         "retained_rate_source_decisions": rate_rows,
                         "potentially_influential_under_declared_scenarios": response.get(
                             "potentially_influential_under_declared_scenarios"),
                         "evidence_required": ["supported electronic-state/channel interpretation",
                                               "justified return/escape or thermochemistry constraints",
                                               "nonadiabatic rate evidence if the branch remains influential"],
                         "executable_with_current_backend": False, "estimated_cost": None})
    requests.extend([
        {"action": "establish_requested_pressure_applicability_and_rate_uncertainty",
         "reason": "High-pressure central fit is not a requested-pressure or accuracy certificate.",
         "evidence_required": ["condition/composition-matched pressure adequacy", "applicable rate uncertainty"],
         "executable_with_current_backend": False, "estimated_cost": None},
        {"action": "resolve_conformer_and_statistical_counting_basis",
         "reason": "Reached conformer and statistical/path counting are not made complete by network assessment.",
         "evidence_required": ["source-bound conformer/statistical recipe assessment"],
         "executable_with_current_backend": False, "estimated_cost": None},
        {"action": "establish_prediction_conditioned_mechanism_coverage",
         "reason": "The retained model is the discovered inventory, not proof no influential chemistry is missing.",
         "evidence_required": ["observable-conditioned discovery/coverage evidence"],
         "executable_with_current_backend": False, "estimated_cost": None},
    ])
    for item in requests:
        item.update(network_scientific_hash=network["scientific_hash"],
                    condition_contract_id=network["condition_contract"]["contract_id"],
                    information_gain_score=None, physical_prediction_claim_allowed=False)
    if rate_audit is not None:
        for item in requests:
            item["rate_applicability_audit_scientific_hash"] = rate_audit["scientific_hash"]
    return requests


def validate_expanded_network_report(report, condition_contract):
    """Reopen immutable evidence and reconstruct scenario interpretation."""
    from .bounded_handoff import validate_expanded_forward_only_model_assessment
    from .retained_linear_network import validate_retained_linear_network
    from .method_applicability import audit_retained_network_states
    from .retained_rate_audit import validate_retained_network_rate_audit
    if report.get("kind") != REPORT_KIND or report.get("condition_contract") != condition_contract \
            or any(report.get(key) is not value for key, value in _CLAIMS.items()):
        raise ValueError("expanded_network_report_condition_or_claim_mismatch")
    if report.get("implementation_sources") != _implementation_sources():
        raise ValueError("expanded_network_report_implementation_changed")
    for record in report.get("artifacts", {}).values():
        if _artifact(record["path"]) != record:
            raise ValueError("expanded_network_report_artifact_changed")
    for field, key in (("retained_network", "retained_network"),
                       ("expanded_model_assessment", "expanded_model_assessment"),
                       ("scenario_response", "scenario_response"),
                       ("electronic_state_audit", "electronic_state_audit"),
                       ("rate_applicability_audit", "rate_applicability_audit"),
                       ("next_evidence_requests", "next_evidence_requests")):
        record = report["artifacts"][key]
        if json.loads(Path(record["path"]).read_text()) != report[field]:
            raise ValueError("expanded_network_report_sidecar_mismatch")
    network = validate_retained_linear_network(report["retained_network"], condition_contract)
    validate_expanded_forward_only_model_assessment(report["expanded_model_assessment"], condition_contract)
    audit = audit_retained_network_states(network)
    rate_audit = validate_retained_network_rate_audit(report["rate_applicability_audit"], network, audit)
    definitions = build_declared_network_scenarios(network, audit)
    response = build_network_scenario_response(network, audit, report["scenario_response"]["native_diagnostics"], definitions)
    if response != report["scenario_response"] or audit != report["electronic_state_audit"] \
            or next_network_evidence_requests(network, audit, response, rate_audit) != report["next_evidence_requests"]:
        raise ValueError("expanded_network_report_diagnostic_changed")
    allowed = report["expanded_model_assessment"].get("bounded_approximation_allowed") is True
    expected_status = "passed" if allowed and response["status"] == "completed_declared_scenarios" else "unresolved"
    if report.get("status") != expected_status or report.get("policy") != network_influence_policy():
        raise ValueError("expanded_network_report_status_or_policy_mismatch")
    expected = dict(report)
    expected.pop("scientific_hash", None)
    if report.get("scientific_hash") != _digest(expected):
        raise ValueError("expanded_network_report_hash_mismatch")
    return report


def assess_retained_expanded_network(model_file, replacement, condition_contract, workdir, *,
                                     rmg_env, bridge=None):
    """Assess a retained full model, without discovery or electronic jobs."""
    from .bounded_handoff import assess_expanded_forward_only_model
    from .retained_linear_network import build_retained_linear_network
    from .reaction_correspondence import resolve_replacement_correspondence
    from .flux_verification import verify_applied_replacement
    from .method_applicability import audit_retained_network_states
    from .retained_rate_audit import audit_retained_network_rates
    from .rmg_bridge_client import run_rmg_bridge

    bridge = run_rmg_bridge if bridge is None else bridge
    model_file, root = Path(model_file).resolve(strict=True), Path(workdir).resolve()
    document = json.loads(model_file.read_text())
    if document.get("condition_contract") != condition_contract:
        raise ValueError("expanded_network_condition_mismatch")
    evidence = document["rmg_evidence"]
    correspondence = resolve_replacement_correspondence(evidence, replacement)
    applied = verify_applied_replacement(evidence, replacement, condition_contract)
    if correspondence.get("status") != "matched" or applied.get("status") != "passed":
        raise ValueError("expanded_network_native_insertion_or_identity_unresolved")
    report_path = root / "assessment-report.json"
    if report_path.exists():
        saved = json.loads(report_path.read_text())
        # Bind cache lookup to the actual call, not merely the cached conditions.
        if saved.get("input_model") != _artifact(model_file) or saved.get("replacement_digest") != _digest(replacement):
            raise ValueError("expanded_network_cache_input_changed")
        return validate_expanded_network_report(saved, condition_contract)
    artifacts = evidence["artifacts"]
    for source in (model_file, Path(replacement["library"]), *(Path(artifacts[key]) for key in
                    ("chemkin", "species_dictionary") if artifacts.get(key))):
        if root == source or root.is_relative_to(source):
            raise ValueError("expanded_network_workspace_must_not_modify_source")
    mapping = correspondence["correspondence"]["species_correspondence"]
    parent = mapping["reactants"][0]["model_label"]
    species_ids = [item["chemkin_identifier"] for item in evidence["mechanism_inspection"]["species"]]
    native_payload = {"chemkin": artifacts["chemkin"], "dictionary": artifacts["species_dictionary"],
                      "transport": artifacts.get("transport"),
                      "temperature_K": condition_contract["temperature_K"],
                      "pressure_bar": condition_contract["pressure_bar"],
                      "condition_contract_id": condition_contract["contract_id"],
                      "initial_fractions": {item: float(item == parent) for item in species_ids},
                      "parent_id": parent, "target_duration_seconds": condition_contract["target_duration_seconds"]}
    snapshot = bridge({**native_payload, "operation": "linear_model_snapshot",
                       "include_rate_source_snapshot": True}, rmg_env=rmg_env)
    snapshot_record = _save(root / "native-rate-snapshot.json", snapshot)
    network = build_retained_linear_network(model_file, replacement, condition_contract, snapshot_record["path"])
    records = {"native_rate_snapshot": snapshot_record,
               "retained_network": _save(root / "retained-network.json", network)}
    retained_rate_file = replacement["forward_only_approximation_evidence"]["retained_rate_evidence"]["path"]
    assessment = assess_expanded_forward_only_model(model_file, retained_rate_file, replacement,
        condition_contract, snapshot_record["path"])
    records["expanded_model_assessment"] = _save(root / "expanded-model-assessment.json", assessment)
    audit = audit_retained_network_states(network)
    rate_audit = audit_retained_network_rates(network, audit)
    records["rate_applicability_audit"] = _save(root / "rate-applicability-audit.json", rate_audit)
    definitions = build_declared_network_scenarios(network, audit)
    native_result = bridge({**native_payload, "operation": "linear_model_diagnostics",
                            "inventory": network, "source_artifacts": network["source_artifacts"],
                            "scenarios": definitions}, rmg_env=rmg_env)
    response = build_network_scenario_response(network, audit, native_result, definitions)
    records["scenario_response"] = _save(root / "scenario-response.json", response)
    records["electronic_state_audit"] = _save(root / "electronic-state-audit.json", audit)
    requests = next_network_evidence_requests(network, audit, response, rate_audit)
    records["next_evidence_requests"] = _save(root / "next-evidence-requests.json", requests)
    allowed = assessment.get("bounded_approximation_allowed") is True
    report = {"schema_version": 2, "kind": REPORT_KIND,
              "status": "passed" if allowed and response["status"] == "completed_declared_scenarios" else "unresolved",
              "condition_contract": copy.deepcopy(condition_contract), "input_model": _artifact(model_file),
              "replacement_digest": _digest(replacement), "policy": network_influence_policy(),
              "implementation_sources": _implementation_sources(), "retained_network": network,
              "expanded_model_assessment": assessment, "electronic_state_audit": audit,
              "rate_applicability_audit": rate_audit,
              "scenario_response": response, "next_evidence_requests": requests,
              "artifacts": records, "chemistry_jobs_executed": 0, "arkane_jobs_executed": 0,
              "rmg_discovery_jobs_executed": 0, **_CLAIMS}
    report["scientific_hash"] = _digest(report)
    _save(report_path, report)
    return report
