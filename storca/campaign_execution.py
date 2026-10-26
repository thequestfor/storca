"""Execution orchestration for stability-campaign lane adapters.

Schema v1 has live adapters for the existing homogeneous-gas ladder and the
restricted neutral implicit-solution model. Every other lane remains an
explicit, retained non-execution result.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .campaign_planning import prepare_stability_campaign
from .gas_campaign import GasCampaignExecutor
from .solution_campaign import SolutionCampaignExecutor
from .stability_campaign import CampaignLaneRequest, StabilityCampaignRequest


CAMPAIGN_EXECUTION_SCHEMA_VERSION = 1
SUPPORTED_GLOBAL_CAMPAIGN_BUDGETS = {
    "maximum_gas_stages", "maximum_solution_lanes",
    "maximum_solution_estimated_core_hours",
}


def _write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


def _gas_observable_claims(
    request: StabilityCampaignRequest,
    gas_result: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Emit only condition-scoped observables that crossed existing verification gates."""
    if gas_result.get("status") != "completed":
        return []
    report = gas_result.get("report") or {}
    observables = {item.observable_id: item for item in request.observables}
    lanes = {item.lane_request_id: item for item in request.lane_requests}
    claims: list[dict[str, Any]] = []
    for stage in report.get("stages") or []:
        lane_request_id = str(stage.get("id") or "")
        lane = lanes.get(lane_request_id)
        lane_result = stage.get("lane_result") or {}
        if lane is None or lane_result.get("observable_status") != "verified":
            continue
        for observable_id in lane.observable_ids:
            observable = observables[observable_id]
            if observable.kind != "t95":
                continue
            claims.append({
                "schema_version": CAMPAIGN_EXECUTION_SCHEMA_VERSION,
                "kind": "storca_condition_scoped_observable_claim",
                "claim_id": f"claim:{request.contract_id[7:23]}:{lane_request_id}:{observable_id}",
                "campaign_contract_id": request.contract_id,
                "lane_request_id": lane_request_id,
                "lane_id": lane.lane_id,
                "environment_contract_id": lane.environment.as_dict()["contract_id"],
                "observable_id": observable_id,
                "observable_kind": observable.kind,
                "value": lane_result["observable_value"],
                "units": lane_result["observable_units"],
                "status": "verified",
                "evidence": {
                    "gas_execution_contract_id": (
                        gas_result.get("execution_contract") or {}
                    ).get("contract_id"),
                    "gas_result_json": gas_result.get("result_json"),
                    "stage_id": lane_request_id,
                },
                "scope": (
                    "Conditional on this retained environment, discovered network, "
                    "verified rates, propagation, and convergence evidence."
                ),
                "not_a_claim": "No cross-environment or universal lifetime is inferred.",
            })
    return claims


def execute_stability_campaign(
    request: StabilityCampaignRequest,
    campaign_dir: Path,
    *,
    request_base_dir: Path,
    ncores: int = 1,
    rmg_env: str = "rmg_env",
    graph_xtb_executable: str | None = None,
    screen_tier: str = "quick-screen",
    rmg_start_method: str = "auto",
    rmg_database_libraries: tuple[str, ...] = (),
    auto_reference_libraries: bool = True,
    gas_executor: GasCampaignExecutor | None = None,
    solution_executor: SolutionCampaignExecutor | None = None,
) -> dict[str, Any]:
    """Prepare, select, and execute admitted gas and restricted solution lanes.

    Execution requires the request's immutable policy and this function call;
    either control alone is insufficient. Unsupported lanes and observables are
    retained without launching substitute physics.
    """
    if request.campaign_policy["execution_policy"] != "explicit_execute":
        raise ValueError(
            "Campaign request is prepare_only; set execution_policy='explicit_execute' "
            "and use the explicit execution command"
        )
    unsupported_global = sorted(
        set(request.global_compute_budget) - SUPPORTED_GLOBAL_CAMPAIGN_BUDGETS
    )
    if unsupported_global:
        raise ValueError(
            "Campaign execution cannot claim to enforce global budgets: "
            f"{unsupported_global}"
        )
    prepared = prepare_stability_campaign(request, campaign_dir)
    plan = json.loads(Path(prepared["plan_json"]).read_text())
    plan_by_id = {item["lane_request_id"]: item for item in plan["lane_plans"]}
    blocked = [
        item for item in plan["lane_plans"] if item["status"] != "ready_to_prepare"
    ]
    if blocked and not request.campaign_policy["continue_after_blocked_lane"]:
        raise ValueError(
            "Campaign has blocked lanes and continue_after_blocked_lane is false: "
            + ", ".join(item["lane_request_id"] for item in blocked)
        )

    gas_lanes: list[CampaignLaneRequest] = []
    solution_lanes: list[CampaignLaneRequest] = []
    lane_execution: list[dict[str, Any]] = []
    for lane in request.lane_requests:
        lane_plan = plan_by_id[lane.lane_request_id]
        if lane_plan["status"] != "ready_to_prepare":
            lane_execution.append({
                "lane_request_id": lane.lane_request_id,
                "lane_id": lane.lane_id,
                "status": "not_executed",
                "reason": f"planning_status:{lane_plan['status']}",
            })
            continue
        if lane.lane_id in {
            "solution-neutral-implicit-v1", "solution-neutral-implicit-mass-action-v1",
        }:
            requested_kinds = {
                item["observable_kind"] for item in lane_plan["observable_plans"]
                if item["status"] == "ready_to_prepare"
            }
            if (
                requested_kinds == {"t95"}
                or (
                    lane.lane_id == "solution-neutral-implicit-mass-action-v1"
                    and requested_kinds
                    <= {
                        "t95", "retention_at_duration", "concentration_at_time",
                        "product_yield", "dominant_product", "branching_fraction",
                        "pathway_flux",
                    }
                )
            ):
                solution_lanes.append(lane)
            else:
                lane_execution.append({
                    "lane_request_id": lane.lane_request_id,
                    "lane_id": lane.lane_id,
                    "status": "adapter_not_implemented",
                    "reason": "The live solution campaign adapter currently executes t95 only.",
                })
            continue
        if not lane.lane_id.startswith("gas-"):
            lane_execution.append({
                "lane_request_id": lane.lane_request_id,
                "lane_id": lane.lane_id,
                "status": "adapter_not_implemented",
                "reason": (
                    "The lane is scientifically admitted for preparation, but campaign "
                    "schema v1 has no live execution adapter for this environment."
                ),
            })
            continue
        requested_kinds = {
            item["observable_kind"] for item in lane_plan["observable_plans"]
            if item["status"] == "ready_to_prepare"
        }
        if requested_kinds != {"t95"}:
            lane_execution.append({
                "lane_request_id": lane.lane_request_id,
                "lane_id": lane.lane_id,
                "status": "adapter_not_implemented",
                "reason": "The live gas campaign adapter currently executes t95 only.",
            })
            continue
        gas_lanes.append(lane)

    if not gas_lanes:
        gas_result: dict[str, Any] = {
            "status": "not_requested",
            "reason": "No admitted t95 gas lane was requested",
            "execution_contract": None,
            "result_json": None,
        }
    else:
        executor = gas_executor or GasCampaignExecutor()
        try:
            gas_result = executor.execute(
                request,
                tuple(gas_lanes),
                Path(campaign_dir) / "gas-execution",
                request_base_dir=Path(request_base_dir),
                ncores=ncores,
                rmg_env=rmg_env,
                graph_xtb_executable=graph_xtb_executable,
                screen_tier=screen_tier,
                rmg_start_method=rmg_start_method,
                rmg_database_libraries=rmg_database_libraries,
                auto_reference_libraries=auto_reference_libraries,
            )
        except Exception as error:
            gas_result = {
                "status": "failed_closed",
                "reason": f"{type(error).__name__}: {error}",
                "execution_contract": None,
                "result_json": None,
            }
    if not solution_lanes:
        solution_result: dict[str, Any] = {
            "status": "not_requested",
            "reason": "No admitted restricted implicit-solution t95 lane was requested",
            "lane_executions": [],
            "observable_claims": [],
        }
    else:
        executor = solution_executor or SolutionCampaignExecutor()
        try:
            solution_result = executor.execute(
                request,
                tuple(solution_lanes),
                Path(campaign_dir) / "solution-execution",
                request_base_dir=Path(request_base_dir),
                ncores=ncores,
            )
        except Exception as error:
            solution_result = {
                "status": "failed_closed",
                "reason": f"{type(error).__name__}: {error}",
                "lane_executions": [],
                "observable_claims": [],
            }
    report = gas_result.get("report") or {}
    stages = {str(item.get("id")): item for item in report.get("stages") or []}
    for lane in gas_lanes:
        stage = stages.get(lane.lane_request_id) or {}
        lane_execution.append({
            "lane_request_id": lane.lane_request_id,
            "lane_id": lane.lane_id,
            "status": (
                str((stage.get("lane_result") or {}).get("execution_status"))
                if stage else str(gas_result.get("status"))
            ),
            "observable_status": (stage.get("lane_result") or {}).get("observable_status"),
            "lane_result": stage.get("lane_result"),
            "evidence_artifact": gas_result.get("result_json"),
        })
    solution_by_id = {
        item["lane_request_id"]: item
        for item in solution_result.get("lane_executions") or []
    }
    for lane in solution_lanes:
        retained = solution_by_id.get(lane.lane_request_id) or {}
        retained_lane_result = retained.get("lane_result") or {}
        retained_lane_results = retained.get("lane_results") or (
            [retained_lane_result] if retained_lane_result else []
        )
        individual_statuses = [
            item.get("observable_status") for item in retained_lane_results
        ]
        if len(individual_statuses) > 1:
            aggregate_observable_status = (
                "all_requested_observables_resolved"
                if all(
                    item in {"verified", "not_reached_within_horizon"}
                    for item in individual_statuses
                )
                else "one_or_more_requested_observables_unresolved"
            )
        else:
            aggregate_observable_status = retained_lane_result.get("observable_status")
        lane_execution.append({
            "lane_request_id": lane.lane_request_id,
            "lane_id": lane.lane_id,
            "status": str(retained.get("status") or solution_result.get("status")),
            "native_status": retained.get("native_status"),
            "observable_status": aggregate_observable_status,
            "individual_observable_statuses": individual_statuses,
            "lane_result": retained.get("lane_result"),
            "lane_results": retained_lane_results,
            "evidence_artifact": retained.get("result_json"),
            "execution_mode": retained.get("execution_mode"),
            "resumed": retained.get("resumed", False),
            "reason": retained.get("reason"),
        })
    order = {lane.lane_request_id: index for index, lane in enumerate(request.lane_requests)}
    lane_execution.sort(key=lambda item: order[item["lane_request_id"]])
    claims = [
        *_gas_observable_claims(request, gas_result),
        *(solution_result.get("observable_claims") or []),
    ]
    resolved_observable_statuses = {
        "verified", "not_reached_within_horizon", "all_requested_observables_resolved",
    }
    unresolved = [
        item["lane_request_id"] for item in lane_execution
        if item.get("observable_status") not in resolved_observable_statuses
    ]
    adapter_statuses = [
        item.get("status") for item in (gas_result, solution_result)
        if item.get("status") != "not_requested"
    ]
    any_adapter_completed = any(status == "completed" for status in adapter_statuses)
    status = (
        "failed_closed" if not any_adapter_completed
        else "completed_with_unresolved_lanes" if unresolved
        else "completed"
    )
    execution = {
        "schema_version": CAMPAIGN_EXECUTION_SCHEMA_VERSION,
        "kind": "storca_stability_campaign_execution",
        "status": status,
        "campaign_contract_id": request.contract_id,
        "gas_execution": {
            "status": gas_result.get("status"),
            "result_json": gas_result.get("result_json"),
            "execution_contract_id": (
                gas_result.get("execution_contract") or {}
            ).get("contract_id"),
            "resumed": gas_result.get("resumed", False),
            "reason": gas_result.get("reason"),
        },
        "solution_execution": {
            "status": solution_result.get("status"),
            "reason": solution_result.get("reason"),
            "failed_lane_request_ids": solution_result.get(
                "failed_lane_request_ids", [],
            ),
        },
        "lane_executions": lane_execution,
        "observable_claims": claims,
        "unresolved_lane_request_ids": unresolved,
        "resolved_lane_request_ids": [
            item["lane_request_id"] for item in lane_execution
            if item.get("observable_status") in resolved_observable_statuses
        ],
        "cross_environment_aggregation": {
            "status": "not_performed",
            "reason": (
                "Lane observables are conditional on distinct physical models; schema v1 "
                "does not fabricate a minimum or universal lifetime."
            ),
        },
    }
    execution_path = _write_json(
        Path(campaign_dir) / "campaign-execution.json", execution,
    )
    summary = {
        "schema_version": CAMPAIGN_EXECUTION_SCHEMA_VERSION,
        "kind": "storca_stability_campaign_summary",
        "status": status,
        "campaign_contract_id": request.contract_id,
        **plan["summary"],
        "observable_claims_emitted": claims,
        "unresolved_lane_request_ids": unresolved,
        "resolved_lane_request_ids": [
            item["lane_request_id"] for item in lane_execution
            if item.get("observable_status") in resolved_observable_statuses
        ],
        "reason": (
            "Implemented gas and restricted implicit-solution adapters were executed "
            "independently; unsupported lanes and unresolved evidence remain explicit."
        ),
        "execution_json": str(execution_path),
    }
    summary_path = _write_json(Path(campaign_dir) / "campaign-summary.json", summary)
    return {
        "status": status,
        "campaign_contract_id": request.contract_id,
        "campaign_dir": str(campaign_dir),
        "plan_json": prepared["plan_json"],
        "execution_json": str(execution_path),
        "summary_json": str(summary_path),
        "summary": summary,
        "gas_execution": gas_result,
        "solution_execution": solution_result,
    }


def load_and_execute_stability_campaign(
    request_path: Path,
    campaign_dir: Path,
    **kwargs: Any,
) -> dict[str, Any]:
    request_path = Path(request_path)
    return execute_stability_campaign(
        StabilityCampaignRequest.from_json(request_path),
        campaign_dir,
        request_base_dir=request_path.parent,
        **kwargs,
    )


__all__ = [
    "CAMPAIGN_EXECUTION_SCHEMA_VERSION", "execute_stability_campaign",
    "load_and_execute_stability_campaign",
]
