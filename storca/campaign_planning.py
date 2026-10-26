"""Preparation-only planning for immutable multi-lane stability campaigns."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping

from .lanes import (
    IMPLEMENTED_GAS_CAPABILITIES,
    SOLUTION_MASS_ACTION_CAPABILITIES,
    SOLUTION_UNIMOLECULAR_CAPABILITIES,
    admit_lane,
    build_stability_lane_registry,
)
from .stability_campaign import CampaignLaneRequest, StabilityCampaignRequest
from .lanes.gas import GAS_LANE_SCENARIOS


CAMPAIGN_PLAN_SCHEMA_VERSION = 1
CAMPAIGN_PLAN_STATES = {
    "ready_to_prepare", "partially_ready", "inputs_missing",
    "capability_unavailable", "observable_unsupported",
    "outside_method_applicability", "declared_only", "invalid_request",
}


def _write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)
    return path


def _available_capabilities(lane_id: str) -> frozenset[str]:
    if lane_id.startswith("gas-"):
        return IMPLEMENTED_GAS_CAPABILITIES
    if lane_id == "solution-neutral-implicit-v1":
        return SOLUTION_UNIMOLECULAR_CAPABILITIES
    if lane_id == "solution-neutral-implicit-mass-action-v1":
        return SOLUTION_MASS_ACTION_CAPABILITIES
    return frozenset()


def _positive_parameter(
    parameters: Mapping[str, Any], names: tuple[str, ...], *, label: str,
) -> str | None:
    value = next((parameters[name] for name in names if parameters.get(name) is not None), None)
    if value is None:
        return f"missing_{label}"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return f"invalid_{label}"
    if not math.isfinite(number) or number <= 0.0:
        return f"invalid_{label}"
    return None


def _environment_preflight(lane: CampaignLaneRequest) -> list[dict[str, Any]]:
    parameters = lane.environment.parameters
    findings: list[dict[str, Any]] = []
    if lane.lane_id.startswith("gas-"):
        checks = (
            (("temperature_K",), "temperature_K"),
            (("pressure_bar",), "pressure_bar"),
            (("target_duration_seconds", "maximum_duration_seconds"), "duration_seconds"),
        )
        for names, label in checks:
            reason = _positive_parameter(parameters, names, label=label)
            if reason:
                findings.append({"kind": "invalid_environment", "reason": reason})
        retention = parameters.get("retention_fraction")
        try:
            valid_retention = retention is not None and 0.0 < float(retention) < 1.0
        except (TypeError, ValueError):
            valid_retention = False
        if not valid_retention:
            findings.append({
                "kind": "invalid_environment", "reason": "invalid_retention_fraction",
            })
        if lane.environment.model_name != "storca_condition_spec_gas_v1":
            findings.append({
                "kind": "invalid_environment",
                "reason": "unsupported_gas_environment_model",
            })
        expected_scenario = GAS_LANE_SCENARIOS[lane.lane_id]
        if parameters.get("scenario") != expected_scenario:
            findings.append({
                "kind": "invalid_environment",
                "reason": f"scenario_must_be_{expected_scenario}",
            })
        expected_light = (
            "sunlight" if lane.lane_id.startswith("gas-sunlight-") else "dark"
        )
        if parameters.get("light_condition") != expected_light:
            findings.append({
                "kind": "invalid_environment",
                "reason": f"light_condition_must_be_{expected_light}",
            })
        if lane.lane_id in {"gas-humid-air-v1", "gas-sunlight-humid-air-v1"}:
            humidity = parameters.get("relative_humidity")
            try:
                valid_humidity = humidity is not None and 0.0 <= float(humidity) <= 1.0
            except (TypeError, ValueError):
                valid_humidity = False
            if not valid_humidity:
                findings.append({
                    "kind": "invalid_environment",
                    "reason": "invalid_relative_humidity",
                })
        if lane.lane_id.startswith("gas-sunlight-"):
            if not parameters.get("light_model"):
                findings.append({
                    "kind": "missing_input", "reason": "light_model_missing",
                })
            if "photolysis_evidence" not in {
                item.input_kind for item in lane.evidence_inputs
            }:
                findings.append({
                    "kind": "missing_input", "reason": "photolysis_evidence_missing",
                })
    elif lane.lane_id in {
        "solution-neutral-implicit-v1", "solution-neutral-implicit-mass-action-v1",
    }:
        try:
            if lane.lane_id == "solution-neutral-implicit-v1":
                from .solution_kinetics import ImplicitSolutionEnvironment

                solution_environment = ImplicitSolutionEnvironment.from_dict(parameters)
            else:
                from .solution_mass_action import MassActionSolutionEnvironment

                solution_environment = MassActionSolutionEnvironment.from_dict(parameters)
            if (
                lane.method_profile is not None
                and solution_environment.electronic_structure_profile != lane.method_profile
            ):
                findings.append({
                    "kind": "method_environment_mismatch",
                    "reason": "method_profile_does_not_match_solution_environment",
                })
        except (TypeError, ValueError, RuntimeError) as error:
            findings.append({
                "kind": "invalid_environment",
                "reason": f"{type(error).__name__}: {error}",
            })
    return findings


def _status_for_admission(
    *, support_level: str, admission: Mapping[str, Any], findings: list[dict[str, Any]],
    method_profile: str | None,
) -> str:
    if support_level == "declared_only":
        return "declared_only"
    if any(item["kind"] in {"invalid_environment", "method_environment_mismatch"} for item in findings):
        return "outside_method_applicability"
    if any(item["kind"] == "unresolved_user_applicability_flag" for item in findings):
        return "outside_method_applicability"
    if method_profile is None or any(item["kind"] == "missing_input" for item in findings):
        return "inputs_missing"
    status = admission["status"]
    reason = str(admission.get("reason") or "")
    if status == "supported":
        return "ready_to_prepare"
    if status == "inputs_missing":
        return "inputs_missing"
    if status == "capability_unavailable":
        return "capability_unavailable"
    if "does not support requested observable" in reason:
        return "observable_unsupported"
    return "outside_method_applicability"


def _observable_environment_findings(
    lane: CampaignLaneRequest, observable: Any,
) -> list[dict[str, Any]]:
    parameters = lane.environment.parameters
    findings: list[dict[str, Any]] = []
    environment_retention = parameters.get("retention_fraction")
    if observable.retention_fraction is not None and environment_retention is not None:
        try:
            matches = math.isclose(
                float(observable.retention_fraction), float(environment_retention),
                rel_tol=0.0, abs_tol=1.0e-12,
            )
        except (TypeError, ValueError):
            matches = False
        if not matches:
            findings.append({
                "kind": "observable_environment_mismatch",
                "reason": "retention_fraction_does_not_match_lane_environment",
            })
    duration = next((
        parameters[name]
        for name in ("target_duration_seconds", "maximum_duration_seconds")
        if parameters.get(name) is not None
    ), None)
    if observable.time_horizon_seconds is not None and duration is not None:
        try:
            matches = math.isclose(
                float(observable.time_horizon_seconds), float(duration),
                rel_tol=1.0e-12, abs_tol=1.0e-12,
            )
        except (TypeError, ValueError):
            matches = False
        if not matches:
            findings.append({
                "kind": "observable_environment_mismatch",
                "reason": "time_horizon_does_not_match_lane_environment",
            })
    return findings


def _plan_lane(
    request: StabilityCampaignRequest,
    lane: CampaignLaneRequest,
    observables_by_id: Mapping[str, Any],
) -> dict[str, Any]:
    definition = build_stability_lane_registry().get(lane.lane_id)
    capabilities = _available_capabilities(lane.lane_id)
    available_inputs = tuple(item.input_kind for item in lane.evidence_inputs)
    findings = _environment_preflight(lane)
    if (
        len(request.chemical_system.components) != 1
        and lane.lane_id != "solution-neutral-implicit-mass-action-v1"
    ):
        findings.append({
            "kind": "invalid_environment",
            "reason": "multi_component_campaign_execution_not_supported",
        })
    if lane.lane_id == "solution-neutral-implicit-v1":
        target = request.chemical_system.target
        if target.charge != 0 or target.multiplicity != 1:
            findings.append({
                "kind": "invalid_environment",
                "reason": "solution_lane_requires_neutral_closed_shell_target",
            })
        input_kinds = {item.input_kind for item in lane.evidence_inputs}
        direct = {"solution_environment", "verified_solution_network"}
        session = "solution_verification_session_request"
        if session in input_kinds and input_kinds & direct:
            findings.append({
                "kind": "invalid_environment",
                "reason": "ambiguous_solution_evidence_modes",
            })
        supported = {*direct, session}
        if input_kinds - supported:
            findings.append({
                "kind": "invalid_environment",
                "reason": "unsupported_solution_evidence_kind",
            })
        if session in input_kinds:
            try:
                temperature = float(lane.environment.parameters["temperature_K"])
                session_temperature_valid = math.isclose(
                    temperature, 298.15, rel_tol=0.0, abs_tol=0.05,
                )
            except (KeyError, TypeError, ValueError):
                session_temperature_valid = False
            if not session_temperature_valid:
                findings.append({
                    "kind": "invalid_environment",
                    "reason": "solution_verification_session_requires_298.15_K",
                })
    elif lane.lane_id == "solution-neutral-implicit-mass-action-v1":
        target = request.chemical_system.target
        if target.charge != 0 or target.multiplicity != 1:
            findings.append({
                "kind": "invalid_environment",
                "reason": "mass_action_solution_lane_requires_neutral_closed_shell_target",
            })
        input_kinds = {item.input_kind for item in lane.evidence_inputs}
        supported = {"solution_environment", "verified_solution_mass_action_network"}
        if input_kinds - supported:
            findings.append({
                "kind": "invalid_environment",
                "reason": "unsupported_mass_action_solution_evidence_kind",
            })
    for name, active in request.chemical_system.user_declared_applicability_flags.items():
        if active:
            findings.append({
                "kind": "unresolved_user_applicability_flag",
                "reason": name,
            })
    if lane.environment.environment_kind != definition.environment_kind:
        findings.append({
            "kind": "invalid_environment",
            "reason": (
                f"lane_requires_{definition.environment_kind}_not_"
                f"{lane.environment.environment_kind}"
            ),
        })
    if lane.method_profile is None and definition.support_level != "declared_only":
        findings.append({"kind": "missing_input", "reason": "method_profile_missing"})

    observable_plans = []
    for observable_id in lane.observable_ids:
        observable = observables_by_id[observable_id]
        observable_findings = _observable_environment_findings(lane, observable)
        admission = admit_lane(
            definition,
            lane.environment,
            requested_observable=observable.kind,
            available_capabilities=capabilities,
            available_inputs=available_inputs,
        ).as_dict()
        status = _status_for_admission(
            support_level=definition.support_level,
            admission=admission,
            findings=[*findings, *observable_findings],
            method_profile=lane.method_profile,
        )
        if observable_findings and status == "ready_to_prepare":
            status = "outside_method_applicability"
        observable_plans.append({
            "observable_id": observable_id,
            "observable_kind": observable.kind,
            "status": status,
            "admission": admission,
            "preflight_findings": observable_findings,
            "ready_for_execution_adapter": status == "ready_to_prepare",
        })

    states = {item["status"] for item in observable_plans}
    if states == {"ready_to_prepare"}:
        lane_status = "ready_to_prepare"
    elif "ready_to_prepare" in states:
        lane_status = "partially_ready"
    else:
        priority = (
            "declared_only", "outside_method_applicability", "inputs_missing",
            "capability_unavailable", "observable_unsupported", "invalid_request",
        )
        lane_status = next((item for item in priority if item in states), sorted(states)[0])
    if lane_status not in CAMPAIGN_PLAN_STATES:
        raise AssertionError(f"Unknown campaign lane planning state: {lane_status}")

    return {
        "schema_version": CAMPAIGN_PLAN_SCHEMA_VERSION,
        "kind": "storca_stability_campaign_lane_plan",
        "lane_request_id": lane.lane_request_id,
        "lane_id": lane.lane_id,
        "status": lane_status,
        "lane_definition": definition.as_dict(),
        "environment_contract": lane.environment.as_dict(),
        "method_profile": lane.method_profile,
        "compute_budget": dict(lane.compute_budget),
        "evidence_input_kinds": list(available_inputs),
        "preflight_findings": findings,
        "observable_plans": observable_plans,
        "planned_scientific_stages": {
            "discovery_backends": list(definition.discovery_backends),
            "path_verifier": definition.path_verifier,
            "rate_model": definition.rate_model,
            "propagator": definition.propagator,
            "convergence_policy": definition.convergence_policy,
        },
        "cost_estimate": {
            "status": "not_estimated",
            "reason": (
                "No calculation count is fabricated before a lane adapter resolves "
                "candidate and evidence inventories. Immutable configured ceilings are "
                "retained in compute_budget."
            ),
        },
        "external_programs_launched": [],
    }


def plan_stability_campaign(request: StabilityCampaignRequest) -> dict[str, Any]:
    """Return a deterministic scientific plan without launching calculations."""
    observables = {item.observable_id: item for item in request.observables}
    lane_plans = [
        _plan_lane(request, lane, observables) for lane in request.lane_requests
    ]
    counts: dict[str, int] = {}
    for lane_plan in lane_plans:
        counts[lane_plan["status"]] = counts.get(lane_plan["status"], 0) + 1
    return {
        "schema_version": CAMPAIGN_PLAN_SCHEMA_VERSION,
        "kind": "storca_stability_campaign_plan",
        "status": "prepared",
        "execution_policy": request.campaign_policy["execution_policy"],
        "campaign_contract_id": request.contract_id,
        "chemical_system_contract_id": request.chemical_system.as_dict()["contract_id"],
        "global_compute_budget": dict(request.global_compute_budget),
        "summary": {
            "lane_request_count": len(lane_plans),
            "status_counts": dict(sorted(counts.items())),
            "ready_lane_request_ids": [
                item["lane_request_id"]
                for item in lane_plans if item["status"] == "ready_to_prepare"
            ],
            "blocked_lane_request_ids": [
                item["lane_request_id"]
                for item in lane_plans if item["status"] != "ready_to_prepare"
            ],
            "scientific_interpretation": (
                "Lane plans are independent conditional models. This preparation does not "
                "combine their future observable values or establish reaction coverage."
            ),
        },
        "lane_plans": lane_plans,
        "external_programs_launched": [],
    }


def prepare_stability_campaign(
    request: StabilityCampaignRequest,
    campaign_dir: Path,
) -> dict[str, Any]:
    """Retain or reproduce a preparation-only campaign under an exact contract."""
    campaign_dir = Path(campaign_dir)
    campaign_dir.mkdir(parents=True, exist_ok=True)
    request_path = campaign_dir / "campaign-request.json"
    if request_path.is_file():
        retained = StabilityCampaignRequest.from_dict(json.loads(request_path.read_text()))
        if retained.contract_id != request.contract_id:
            raise ValueError(
                "Campaign directory contains a different immutable scientific request"
            )
    else:
        _write_json(request_path, request.as_dict())

    metadata_path = campaign_dir / "campaign-metadata.json"
    metadata_history: list[dict[str, Any]] = []
    if metadata_path.is_file():
        metadata_history = list(json.loads(metadata_path.read_text()).get("observed_metadata") or [])
    metadata = dict(request.metadata)
    if metadata not in metadata_history:
        metadata_history.append(metadata)
    _write_json(metadata_path, {
        "schema_version": 1,
        "kind": "storca_stability_campaign_metadata_history",
        "campaign_contract_id": request.contract_id,
        "observed_metadata": metadata_history,
        "note": "Metadata is retained but excluded from the immutable scientific contract hash.",
    })

    plan = plan_stability_campaign(request)
    plan_path = _write_json(campaign_dir / "campaign-plan.json", plan)
    lane_dir = campaign_dir / "lane-plans"
    lane_artifacts = {}
    for lane_plan in plan["lane_plans"]:
        path = _write_json(
            lane_dir / f"{lane_plan['lane_request_id']}.json", lane_plan,
        )
        lane_artifacts[lane_plan["lane_request_id"]] = str(path)
    summary = {
        "schema_version": 1,
        "kind": "storca_stability_campaign_summary",
        "status": "prepared",
        "campaign_contract_id": request.contract_id,
        **plan["summary"],
        "observable_claims_emitted": [],
        "reason": "Preparation-only schema v1 launched no chemistry calculations.",
    }
    summary_path = _write_json(campaign_dir / "campaign-summary.json", summary)
    return {
        "status": "prepared",
        "campaign_contract_id": request.contract_id,
        "campaign_dir": str(campaign_dir),
        "request_json": str(request_path),
        "plan_json": str(plan_path),
        "summary_json": str(summary_path),
        "lane_plan_artifacts": lane_artifacts,
        "summary": summary,
        "external_programs_launched": [],
    }


def load_and_prepare_stability_campaign(request_path: Path, campaign_dir: Path) -> dict[str, Any]:
    return prepare_stability_campaign(
        StabilityCampaignRequest.from_json(request_path), campaign_dir,
    )


__all__ = [
    "CAMPAIGN_PLAN_SCHEMA_VERSION", "CAMPAIGN_PLAN_STATES",
    "load_and_prepare_stability_campaign", "plan_stability_campaign",
    "prepare_stability_campaign",
]
