"""Campaign adapter for the existing bounded gas stability machinery."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Mapping

from .condition_ladder import humid_air_scenario, run_gas_stage_sequence
from .campaign_evidence import evidence_contract_records, resolve_lane_evidence
from .lanes import IMPLEMENTED_GAS_CAPABILITIES, admit_lane, build_stability_lane_registry
from .lanes.gas import GAS_LANE_SCENARIOS
from .methods import harmonic_method_profile
from .stability_campaign import CampaignLaneRequest, StabilityCampaignRequest


GAS_EXECUTION_SCHEMA_VERSION = 1
ENFORCEABLE_GAS_BUDGETS = {
    "rmg_attempt_timeout_seconds",
    "rmg_max_processes",
    "rmg_initial_max_iterations",
    "rmg_initial_max_edge_species",
    "rmg_max_objects_per_iteration",
    "rmg_maximum_heavy_atoms",
    "rmg_attempts",
    "rmg_heartbeat_timeout_seconds",
    "verification_max_iterations",
    "verification_path_timeout_seconds",
    "verification_orientations",
    "verification_neb_images",
}
ENFORCEABLE_GLOBAL_GAS_BUDGETS = {"maximum_gas_stages"}
KNOWN_NON_GAS_GLOBAL_BUDGETS = {
    "maximum_solution_lanes", "maximum_solution_estimated_core_hours",
}


def _write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


def _canonical_hash(value: Mapping[str, Any]) -> str:
    import hashlib

    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False, default=str,
    )
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _duration_seconds(parameters: Mapping[str, Any]) -> float:
    value = parameters.get("maximum_duration_seconds")
    if value is None:
        value = parameters.get("target_duration_seconds")
    duration = float(value)
    if duration <= 0.0:
        raise ValueError("Gas lane duration must be positive")
    return duration


def gas_stage_from_campaign_lane(lane: CampaignLaneRequest) -> dict[str, Any]:
    """Translate one admitted campaign lane into an existing gas-stage record."""
    definition = build_stability_lane_registry().get(lane.lane_id)
    if not lane.lane_id.startswith("gas-"):
        raise ValueError(f"Not a gas campaign lane: {lane.lane_id}")
    parameters = lane.environment.parameters
    temperature = float(parameters["temperature_K"])
    pressure = float(parameters["pressure_bar"])
    retention = float(parameters["retention_fraction"])
    duration = _duration_seconds(parameters)
    legacy_id = str(definition.metadata.get("legacy_stage_id") or lane.lane_id)
    expected_scenario = GAS_LANE_SCENARIOS[lane.lane_id]
    if parameters.get("scenario") != expected_scenario:
        raise ValueError(
            f"Gas lane {lane.lane_request_id} requires scenario={expected_scenario!r}"
        )
    expected_light = "sunlight" if lane.lane_id.startswith("gas-sunlight-") else "dark"
    if parameters.get("light_condition") != expected_light:
        raise ValueError(
            f"Gas lane {lane.lane_request_id} requires light_condition={expected_light!r}"
        )
    stage: dict[str, Any] = {
        "id": lane.lane_request_id,
        "campaign_lane_request_id": lane.lane_request_id,
        "campaign_compute_budget": dict(lane.compute_budget),
        "legacy_stage_id": legacy_id,
        "kind": "photolysis" if lane.lane_id.startswith("gas-sunlight-") else "rmg",
        "temperature_K": temperature,
        "pressure_bar": pressure,
        "retention_fraction": retention,
        "maximum_duration_seconds": duration,
        "lane_id": lane.lane_id,
        "lane_contract": definition.as_dict(),
        "environment_contract": lane.environment.as_dict(),
    }
    if lane.lane_id in {"gas-humid-air-v1", "gas-sunlight-humid-air-v1"}:
        if parameters.get("relative_humidity") is None:
            raise ValueError(f"Humid gas lane {lane.lane_request_id} needs relative_humidity")
        humidity = float(parameters["relative_humidity"])
        stage["scenario_config"] = humid_air_scenario(
            temperature_K=temperature,
            pressure_bar=pressure,
            relative_humidity=humidity,
        )
        stage["relative_humidity"] = humidity
    else:
        stage["scenario"] = expected_scenario
    if stage["kind"] == "photolysis":
        stage["light_model"] = parameters.get("light_model")
        if not stage["light_model"]:
            raise ValueError(f"Sunlight gas lane {lane.lane_request_id} needs light_model")
    admission = admit_lane(
        definition,
        lane.environment,
        requested_observable="t95",
        available_capabilities=IMPLEMENTED_GAS_CAPABILITIES,
    )
    if not admission.permits_observable:
        raise ValueError(
            f"Gas lane {lane.lane_request_id} is not admitted: {admission.reason}"
        )
    stage["lane_admission"] = admission.as_dict()
    return stage


class GasCampaignExecutor:
    """Execute selected gas lanes through existing RMG/ORCA/Arkane propagation."""

    def __init__(
        self,
        *,
        stage_runner: Callable[..., dict] | None = None,
        graph_discoverer: Callable[..., Mapping[str, Any]] | None = None,
    ) -> None:
        self._injected_stage_runner = stage_runner
        self._graph_discoverer = graph_discoverer

    @staticmethod
    def _validate_budgets(
        request: StabilityCampaignRequest, lanes: tuple[CampaignLaneRequest, ...],
    ) -> None:
        unsupported_global = sorted(
            set(request.global_compute_budget)
            - ENFORCEABLE_GLOBAL_GAS_BUDGETS
            - KNOWN_NON_GAS_GLOBAL_BUDGETS
        )
        if unsupported_global:
            raise ValueError(
                "Gas execution cannot claim to enforce global budgets: "
                f"{unsupported_global}. Use preparation-only mode until a global ledger exists."
            )
        maximum_stages = request.global_compute_budget.get("maximum_gas_stages")
        if maximum_stages is not None:
            if int(maximum_stages) != maximum_stages or int(maximum_stages) < 1:
                raise ValueError("maximum_gas_stages must be a positive integer")
            if len(lanes) > int(maximum_stages):
                raise ValueError("Requested gas lanes exceed maximum_gas_stages")
        integer_budgets = {
            "rmg_max_processes", "rmg_initial_max_iterations",
            "rmg_initial_max_edge_species",
            "rmg_max_objects_per_iteration", "rmg_maximum_heavy_atoms",
            "rmg_attempts", "verification_max_iterations",
            "verification_orientations", "verification_neb_images",
        }
        positive_time_budgets = {
            "rmg_attempt_timeout_seconds", "rmg_heartbeat_timeout_seconds",
            "verification_path_timeout_seconds",
        }
        for lane in lanes:
            unsupported = sorted(set(lane.compute_budget) - ENFORCEABLE_GAS_BUDGETS)
            if unsupported:
                raise ValueError(
                    f"Gas lane {lane.lane_request_id} contains unenforceable budgets: "
                    f"{unsupported}"
                )
            for name in integer_budgets & set(lane.compute_budget):
                value = lane.compute_budget[name]
                minimum = 1
                if int(value) != value or int(value) < minimum:
                    raise ValueError(
                        f"Gas budget {name} must be an integer >= {minimum}"
                    )
            for name in positive_time_budgets & set(lane.compute_budget):
                if float(lane.compute_budget[name]) <= 0.0:
                    raise ValueError(f"Gas budget {name} must be positive")

    @staticmethod
    def _append_event(root: Path, payload: Mapping[str, Any]) -> Path:
        event_dir = Path(root) / "execution-events"
        event_dir.mkdir(parents=True, exist_ok=True)
        index = len(list(event_dir.glob("*.json")))
        return _write_json(event_dir / f"{index:06d}.json", {
            "schema_version": GAS_EXECUTION_SCHEMA_VERSION,
            "kind": "storca_gas_campaign_execution_event",
            **dict(payload),
        })

    def execute(
        self,
        request: StabilityCampaignRequest,
        gas_lanes: tuple[CampaignLaneRequest, ...],
        workdir: Path,
        *,
        request_base_dir: Path,
        ncores: int = 1,
        rmg_env: str = "rmg_env",
        graph_xtb_executable: str | None = None,
        screen_tier: str = "quick-screen",
        rmg_start_method: str = "auto",
        rmg_database_libraries: tuple[str, ...] = (),
        auto_reference_libraries: bool = True,
    ) -> dict[str, Any]:
        if not gas_lanes:
            raise ValueError("GasCampaignExecutor needs at least one gas lane")
        if ncores < 1:
            raise ValueError("ncores must be positive")
        self._validate_budgets(request, gas_lanes)
        methods = {lane.method_profile for lane in gas_lanes}
        if None in methods or len(methods) != 1:
            raise ValueError("One gas execution sequence requires one shared method_profile")
        method_name = str(next(iter(methods)))
        method_profile = harmonic_method_profile(method_name)
        discovery_policies = {
            json.dumps(lane.discovery_policy, sort_keys=True) for lane in gas_lanes
        }
        if len(discovery_policies) != 1:
            raise ValueError("Shared gas execution requires identical discovery_policy values")
        discovery_policy = dict(gas_lanes[0].discovery_policy)
        resolved = {
            lane.lane_request_id: resolve_lane_evidence(
                lane, base_dir=request_base_dir,
            )
            for lane in gas_lanes
        }
        photo_records = [
            item for evidence in resolved.values()
            for item in evidence.get("photolysis_evidence", [])
        ]
        if any(lane.lane_id.startswith("gas-sunlight-") for lane in gas_lanes):
            if not photo_records:
                raise ValueError("Sunlight gas execution requires photolysis_evidence")
            if len({item["content_sha256"] for item in photo_records}) != 1:
                raise ValueError(
                    "One gas sequence currently requires shared photolysis evidence"
                )
        graph_policy = discovery_policy["graph_xtb"]
        if self._injected_stage_runner is None:
            from src.orca_runner import find_orca

            orca_resolution = find_orca()
        else:
            orca_resolution = "injected_stage_runner_for_test_or_embedding"
        precomputed_graph_records = {
            lane_id: evidence.get("precomputed_graph_xtb_discovery", [])
            for lane_id, evidence in resolved.items()
        }
        for lane_id, records in precomputed_graph_records.items():
            if len(records) > 1:
                raise ValueError(
                    f"Gas lane {lane_id} supplies more than one precomputed graph/xTB artifact"
                )
        needs_live_graph = graph_policy == "bounded" and any(
            not records for records in precomputed_graph_records.values()
        )
        resolved_graph_executable = graph_xtb_executable
        if needs_live_graph and self._graph_discoverer is None:
            from src.orca_runner import find_xtb

            resolved_graph_executable = graph_xtb_executable or find_xtb()
        target = request.chemical_system.target
        execution_payload = {
            "schema_version": GAS_EXECUTION_SCHEMA_VERSION,
            "kind": "storca_gas_campaign_execution_contract",
            "campaign_contract_id": request.contract_id,
            "lane_request_ids": [lane.lane_request_id for lane in gas_lanes],
            "method_profile": method_profile,
            "ncores": ncores,
            "rmg_env": rmg_env,
            "screen_tier": screen_tier,
            "rmg_start_method": rmg_start_method,
            "rmg_database_libraries": list(rmg_database_libraries),
            "auto_reference_libraries": bool(auto_reference_libraries),
            "graph_xtb_executable": resolved_graph_executable,
            "discovery_policy": discovery_policy,
            "lane_compute_budgets": {
                lane.lane_request_id: dict(lane.compute_budget) for lane in gas_lanes
            },
            "evidence": {
                lane_id: evidence_contract_records(evidence)
                for lane_id, evidence in resolved.items()
            },
            "software_resolution": {
                "orca_executable": orca_resolution,
                "rmg_environment": rmg_env,
                "graph_xtb_executable": resolved_graph_executable,
                "note": (
                    "Exact version outputs remain retained by their owning stage artifacts; "
                    "this contract fixes preflight-resolvable paths and selected environments."
                ),
            },
        }
        execution_contract_id = _canonical_hash(execution_payload)
        execution_contract = {**execution_payload, "contract_id": execution_contract_id}
        root = Path(workdir)
        root.mkdir(parents=True, exist_ok=True)
        contract_path = root / "gas-execution-contract.json"
        if contract_path.is_file():
            retained_contract = json.loads(contract_path.read_text())
            if retained_contract.get("contract_id") != execution_contract_id:
                raise ValueError("Gas execution directory contains a different execution contract")
        else:
            _write_json(contract_path, execution_contract)
        state_path = root / "gas-execution-state.json"
        if state_path.is_file():
            state = json.loads(state_path.read_text())
            if state.get("execution_contract_id") != execution_contract_id:
                raise ValueError("Retained gas state does not match the execution contract")
            if state.get("status") == "completed":
                report_path = Path(state["result_json"])
                if report_path.is_file():
                    report = json.loads(report_path.read_text())
                    return {
                        "status": "completed", "resumed": True,
                        "execution_contract": execution_contract,
                        "report": report, "result_json": str(report_path),
                    }
        else:
            state = {
                "schema_version": GAS_EXECUTION_SCHEMA_VERSION,
                "kind": "storca_gas_campaign_execution_state",
                "status": "pending",
                "campaign_contract_id": request.contract_id,
                "execution_contract_id": execution_contract_id,
                "stage_results": {},
            }
        state["status"] = "running"
        _write_json(state_path, state)
        self._append_event(root, {
            "event": "execution_started",
            "execution_contract_id": execution_contract_id,
        })

        def fail_closed(error: Exception) -> dict[str, Any]:
            state["status"] = "failed_closed"
            state["reason"] = f"{type(error).__name__}: {error}"
            _write_json(state_path, state)
            self._append_event(root, {
                "event": "execution_failed_closed", "reason": state["reason"],
            })
            return {
                "status": "failed_closed",
                "reason": state["reason"],
                "execution_contract": execution_contract,
                "result_json": None,
            }

        try:
            graph_by_lane: dict[str, Mapping[str, Any]] = {}
            from .reaction_discovery import DiscoveryResult

            for lane_id, records in precomputed_graph_records.items():
                if not records:
                    continue
                artifact = json.loads(Path(records[0]["path"]).read_text())
                graph_by_lane[lane_id] = DiscoveryResult.from_dict(artifact).as_dict()
            stages = [gas_stage_from_campaign_lane(lane) for lane in gas_lanes]
        except Exception as error:
            return fail_closed(error)
        stage_by_id = {lane.lane_request_id: lane for lane in gas_lanes}
        shared_parent_orca: dict[str, Any] | None = next((
            dict(result["orca_evidence"])
            for result in (state.get("stage_results") or {}).values()
            if isinstance(result, Mapping)
            and isinstance(result.get("orca_evidence"), Mapping)
            and result["orca_evidence"].get("status") == "completed"
        ), None)

        if self._injected_stage_runner is not None:
            injected_stage_runner = self._injected_stage_runner

            def stage_runner(stage_smiles: str, stage_dir: Path, **stage_kwargs) -> dict:
                lane_id = str(stage_kwargs["campaign_lane_request_id"])
                graph_artifact = graph_by_lane.get(lane_id)
                if graph_policy == "bounded" and graph_artifact is None:
                    if self._graph_discoverer is not None:
                        scenario_config = stage_kwargs.get("scenario_config")
                        if scenario_config is None:
                            from .stability import resolve_stability_configuration

                            scenario_config, _ = resolve_stability_configuration(
                                str(stage_kwargs.get("scenario")), screen_tier,
                            )
                        reactive_partners = [
                            {
                                key: item[key]
                                for key in (
                                    "label", "smiles", "charge", "multiplicity",
                                    "reactive",
                                )
                                if key in item
                            }
                            for item in scenario_config.get("additional_species") or []
                            if item.get("reactive") is True
                        ]
                        graph_artifact = dict(self._graph_discoverer(
                            target=target,
                            policy=discovery_policy,
                            system={
                                "smiles": target.smiles,
                                "charge": target.charge,
                                "multiplicity": target.multiplicity,
                                "target_label": "stability",
                                "reactive_partners": reactive_partners,
                            },
                            condition={
                                "lane_request_id": lane_id,
                                "scenario": stage_kwargs.get("scenario"),
                                "scenario_config": scenario_config,
                                "temperature_K": stage_kwargs.get("temperature"),
                                "pressure_bar": stage_kwargs.get("pressure"),
                            },
                            workdir=stage_dir / "reaction-discovery" / "graph-xtb",
                        ))
                    else:
                        stage_kwargs["graph_xtb_discovery"] = True
                        stage_kwargs["graph_xtb_executable"] = resolved_graph_executable
                if graph_artifact is not None:
                    stage_kwargs["precomputed_graph_xtb_discovery"] = graph_artifact
                return injected_stage_runner(stage_smiles, stage_dir, **stage_kwargs)
        else:
            from .kinetic_uncertainty import load_rate_uncertainty_evidence
            from .stability import run_stability_screen

            def stage_runner(stage_smiles: str, stage_dir: Path, **stage_kwargs) -> dict:
                nonlocal shared_parent_orca
                lane_id = str(stage_kwargs.pop("campaign_lane_request_id"))
                lane = stage_by_id[lane_id]
                budget = dict(stage_kwargs.pop("campaign_compute_budget", {}))
                evidence = resolved[lane_id]
                libraries = [
                    Path(item["path"]) for item in evidence.get("reaction_library", [])
                ]
                libraries.extend(stage_kwargs.pop("reaction_libraries", []))
                uncertainty = tuple(
                    record
                    for item in evidence.get("rate_uncertainty_evidence", [])
                    for record in load_rate_uncertainty_evidence(Path(item["path"]))
                )
                graph_artifact = graph_by_lane.get(lane_id)
                result = run_stability_screen(
                    stage_smiles,
                    stage_dir,
                    scenario=stage_kwargs.pop("scenario"),
                    scenario_config=stage_kwargs.pop("scenario_config"),
                    target_duration_hours=stage_kwargs.pop("target_duration_hours"),
                    retention_fraction=stage_kwargs.pop("retention_fraction"),
                    light_condition=stage_kwargs.pop("light_condition", "dark"),
                    light_model=stage_kwargs.pop("light_model", None),
                    temperature=stage_kwargs.pop("temperature"),
                    pressure=stage_kwargs.pop("pressure"),
                    rmg_env=rmg_env,
                    screen_tier=screen_tier,
                    charge=target.charge,
                    multiplicity=target.multiplicity,
                    ncores=ncores,
                    method_profile=method_profile,
                    precomputed_orca_evidence=shared_parent_orca,
                    reaction_libraries=libraries,
                    database_reaction_libraries=list(rmg_database_libraries),
                    auto_reference_libraries=auto_reference_libraries,
                    auto_verify_routes=True,
                    rmg_hard_timeout_seconds=budget.get("rmg_attempt_timeout_seconds"),
                    rmg_max_processes=budget.get("rmg_max_processes"),
                    rmg_max_iterations=budget.get("rmg_initial_max_iterations"),
                    rmg_max_edge_species=budget.get("rmg_initial_max_edge_species"),
                    rmg_max_objects_per_iteration=budget.get("rmg_max_objects_per_iteration"),
                    rmg_maximum_heavy_atoms=budget.get("rmg_maximum_heavy_atoms"),
                    rmg_retry_attempts=int(budget.get("rmg_attempts", 2)),
                    rmg_heartbeat_timeout_seconds=budget.get("rmg_heartbeat_timeout_seconds"),
                    rmg_start_method=rmg_start_method,
                    verification_max_iterations=int(budget.get("verification_max_iterations", 12)),
                    verification_timeout_seconds=budget.get(
                        "verification_path_timeout_seconds", 14400.0,
                    ),
                    verification_orientations=int(budget.get("verification_orientations", 3)),
                    verification_neb_images=int(budget.get("verification_neb_images", 8)),
                    graph_xtb_budget={
                        "max_candidates_for_verification": discovery_policy[
                            "graph_xtb_max_candidates"
                        ],
                        "energy_window_kcal_mol": discovery_policy[
                            "graph_xtb_energy_window_kcal_mol"
                        ],
                    },
                    graph_xtb_discovery=(
                        graph_policy == "bounded" and graph_artifact is None
                    ),
                    graph_xtb_executable=resolved_graph_executable,
                    precomputed_graph_xtb_discovery=graph_artifact,
                    rate_uncertainty_evidence=uncertainty,
                )
                if (
                    shared_parent_orca is None
                    and result.get("orca_evidence", {}).get("status") == "completed"
                ):
                    shared_parent_orca = result["orca_evidence"]
                return result

        photo_path: Path | None = None
        if any(stage["kind"] == "photolysis" for stage in stages):
            photo_path = Path(photo_records[0]["path"])

        def checkpoint(stage: dict, result: dict) -> None:
            state["stage_results"][str(stage["id"])] = result
            state["completed_stage_ids"] = list(state["stage_results"])
            _write_json(state_path, state)
            self._append_event(root, {
                "event": "stage_retained",
                "lane_request_id": stage["id"],
                "stage_status": result.get("status"),
            })

        try:
            report = run_gas_stage_sequence(
                target.smiles,
                root,
                stages=stages,
                stage_runner=stage_runner,
                photolysis_evidence=photo_path,
                computational_light_spectrum=None,
                campaign_mode=str(request.campaign_policy["mode"]),
                retained_stage_results=state.get("stage_results") or {},
                stage_checkpoint=checkpoint,
                rmg_env=rmg_env,
                ncores=ncores,
                charge=target.charge,
                multiplicity=target.multiplicity,
                method_keywords=method_profile["orca_keywords"],
            )
        except Exception as error:
            return fail_closed(error)
        state["status"] = "completed"
        state["result_json"] = str(report["result_json"])
        _write_json(state_path, state)
        self._append_event(root, {
            "event": "execution_completed", "result_json": str(report["result_json"]),
        })
        return {
            "status": "completed",
            "resumed": False,
            "execution_contract": execution_contract,
            "report": report,
            "result_json": str(report["result_json"]),
        }


__all__ = [
    "ENFORCEABLE_GAS_BUDGETS", "ENFORCEABLE_GLOBAL_GAS_BUDGETS",
    "GAS_EXECUTION_SCHEMA_VERSION", "GasCampaignExecutor",
    "gas_stage_from_campaign_lane",
]
