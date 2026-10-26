"""Campaign adapter for restricted homogeneous implicit-solution stability."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Mapping

from .campaign_evidence import evidence_contract_records, resolve_lane_evidence
from .lanes import build_stability_lane_registry
from .lanes.base import LaneResult
from .solution_kinetics import ImplicitSolutionEnvironment, VerifiedSolutionNetwork
from .stability_campaign import CampaignLaneRequest, StabilityCampaignRequest


SOLUTION_EXECUTION_SCHEMA_VERSION = 1
DIRECT_INPUT_KINDS = frozenset({"solution_environment", "verified_solution_network"})
MASS_ACTION_DIRECT_INPUT_KINDS = frozenset({
    "solution_environment", "verified_solution_mass_action_network",
})
SESSION_INPUT_KIND = "solution_verification_session_request"
DIRECT_SOLUTION_BUDGETS = frozenset({"profile_points"})
MASS_ACTION_DIRECT_SOLUTION_BUDGETS = frozenset({
    "profile_points", "relative_tolerance", "absolute_tolerance", "maximum_steps",
})
SESSION_SOLUTION_BUDGETS = frozenset({
    "maximum_orca_paths", "maximum_estimated_core_hours",
})


def _write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


def _canonical_hash(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False, default=str,
    )
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _load_single(
    evidence: Mapping[str, list[dict[str, Any]]], input_kind: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    records = evidence.get(input_kind) or []
    if len(records) != 1:
        raise ValueError(f"Solution execution requires exactly one {input_kind} input")
    path = Path(records[0]["path"])
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise TypeError(f"Solution evidence {input_kind} must contain a JSON object")
    return payload, records[0]


def _evidence_mode(
    evidence: Mapping[str, list[dict[str, Any]]], *, mass_action: bool = False,
) -> str:
    kinds = set(evidence)
    direct_kinds = MASS_ACTION_DIRECT_INPUT_KINDS if mass_action else DIRECT_INPUT_KINDS
    direct_present = bool(kinds & direct_kinds)
    session_present = SESSION_INPUT_KIND in kinds
    if direct_present and session_present:
        raise ValueError("Solution evidence ambiguously supplies both execution modes")
    if session_present:
        if kinds != {SESSION_INPUT_KIND}:
            raise ValueError("Verification-session mode contains unsupported evidence kinds")
        return "verification_session"
    if kinds == set(direct_kinds):
        return "verified_input"
    raise ValueError(
        "Solution execution requires either solution_environment plus "
        "verified_solution_network, or solution_verification_session_request"
    )


def _validate_target(target_smiles: str, request: StabilityCampaignRequest) -> None:
    from .reaction_candidates import canonical_smiles

    if canonical_smiles(target_smiles) != canonical_smiles(
        request.chemical_system.target.smiles
    ):
        raise ValueError("Solution evidence target does not match the campaign target")


def _validate_environment(
    environment: Any,
    lane: CampaignLaneRequest,
) -> None:
    retained = environment.as_environment_spec().as_dict()
    expected = lane.environment.as_dict()
    if retained != expected:
        raise ValueError("Solution evidence environment does not match the campaign lane")
    if environment.electronic_structure_profile != lane.method_profile:
        raise ValueError("Solution evidence method does not match the campaign lane")


def _validated_direct_budget(lane: CampaignLaneRequest, *, mass_action: bool = False) -> dict[str, Any]:
    allowed = MASS_ACTION_DIRECT_SOLUTION_BUDGETS if mass_action else DIRECT_SOLUTION_BUDGETS
    unknown = sorted(set(lane.compute_budget) - allowed)
    if unknown:
        raise ValueError(
            f"Verified-input solution lane contains unenforceable budgets: {unknown}"
        )
    raw = lane.compute_budget.get("profile_points", 180 if mass_action else 160)
    if isinstance(raw, bool) or int(raw) != raw or int(raw) < 20:
        raise ValueError("profile_points must be an integer of at least 20")
    result: dict[str, Any] = {"profile_points": int(raw)}
    if mass_action:
        result.update({
            "relative_tolerance": float(lane.compute_budget.get("relative_tolerance", 1.0e-7)),
            "absolute_tolerance": float(lane.compute_budget.get("absolute_tolerance", 1.0e-12)),
            "maximum_steps": int(lane.compute_budget.get("maximum_steps", 200_000)),
        })
        if not 0.0 < result["relative_tolerance"] <= 1.0e-3:
            raise ValueError("relative_tolerance is outside the audited range")
        if not 0.0 < result["absolute_tolerance"] <= 1.0e-6:
            raise ValueError("absolute_tolerance is outside the audited range")
        if result["maximum_steps"] < 100:
            raise ValueError("maximum_steps must be at least 100")
    return result


def _validated_session_budget(
    lane: CampaignLaneRequest,
    session_request: Mapping[str, Any],
) -> None:
    unknown = sorted(set(lane.compute_budget) - SESSION_SOLUTION_BUDGETS)
    if unknown:
        raise ValueError(
            f"Solution verification-session lane contains unenforceable budgets: {unknown}"
        )
    if set(lane.compute_budget) != set(SESSION_SOLUTION_BUDGETS):
        raise ValueError(
            "Solution verification-session mode requires maximum_orca_paths and "
            "maximum_estimated_core_hours"
        )
    retained = dict(session_request.get("compute_budget") or {})
    if retained != dict(lane.compute_budget):
        raise ValueError(
            "Solution session compute_budget must exactly match the campaign lane budget"
        )


def _lane_result_from_session(
    lane: CampaignLaneRequest,
    result: Mapping[str, Any],
) -> dict[str, Any]:
    definition_contract_id = (
        build_stability_lane_registry().get(lane.lane_id).as_dict()["contract_id"]
    )
    claim = result.get("scientific_claim") or {}
    central = claim.get("central_t95_seconds")
    if result.get("status") == "prediction_converged" and central is not None:
        return LaneResult(
            lane_id=lane.lane_id,
            lane_contract_id=definition_contract_id,
            admission_status="supported",
            support_level="observable_supported",
            execution_status="completed_verified_t95",
            observable_status="verified",
            requested_observable="t95",
            observable_value=float(central),
            observable_units="s",
            reason="Solution verification session satisfied its t95 convergence contract.",
        ).as_dict()
    if result.get("status") == "prediction_converged":
        observable_status = "not_reached_within_horizon"
        reason = "Verified solution rate bounds did not cross t95 within the horizon."
    else:
        observable_status = "unresolved"
        reason = str(result.get("reason") or "Solution verification remains unresolved.")
    return LaneResult(
        lane_id=lane.lane_id,
        lane_contract_id=definition_contract_id,
        admission_status="supported",
        support_level="observable_supported",
        execution_status=str(result.get("status") or "invalid_or_incomplete_evidence"),
        observable_status=observable_status,
        requested_observable="t95",
        reason=reason,
    ).as_dict()


def _condition_scoped_claim(
    request: StabilityCampaignRequest,
    lane: CampaignLaneRequest,
    lane_result: Mapping[str, Any],
    scientific_claim: Mapping[str, Any],
    *,
    mode: str,
    result_json: str | None,
    execution_contract_id: str,
    observable_id: str | None = None,
) -> dict[str, Any] | None:
    if lane_result.get("observable_status") != "verified":
        return None
    observable_id = observable_id or lane.observable_ids[0]
    observable = next(item for item in request.observables if item.observable_id == observable_id)
    observable_kind = observable.kind
    observable_data = dict(lane_result.get("observable_data") or {})
    if observable_kind == "t95":
        bounds = {
            "lower_seconds": scientific_claim.get("lower_t95_seconds")
            or (scientific_claim.get("t95_bounds_seconds") or [None, None])[0],
            "upper_seconds": scientific_claim.get("upper_t95_seconds")
            or (scientific_claim.get("t95_bounds_seconds") or [None, None])[1],
        }
    else:
        bounds = {
            "scenario_min": observable_data.get("rate_bound_scenario_min"),
            "scenario_max": observable_data.get("rate_bound_scenario_max"),
        }
    return {
        "schema_version": SOLUTION_EXECUTION_SCHEMA_VERSION,
        "kind": "storca_condition_scoped_observable_claim",
        "claim_id": f"claim:{request.contract_id[7:23]}:{lane.lane_request_id}:{observable_id}",
        "campaign_contract_id": request.contract_id,
        "lane_request_id": lane.lane_request_id,
        "lane_id": lane.lane_id,
        "environment_contract_id": lane.environment.as_dict()["contract_id"],
        "observable_id": observable_id,
        "observable_kind": observable_kind,
        "target_component_id": observable.target_component_id,
        "value": lane_result["observable_value"],
        "units": lane_result["observable_units"],
        "observable_data": observable_data,
        "bounds": bounds,
        "status": "verified",
        "evidence": {
            "solution_execution_contract_id": execution_contract_id,
            "solution_result_json": result_json,
            "source_claim_id": scientific_claim.get("claim_id"),
            "network_coverage_claim_id": (
                scientific_claim.get("network_coverage_claim_id")
                or scientific_claim.get("conditional_on_coverage_claim_id")
                or scientific_claim.get("conditional_on_network_coverage_claim_id")
            ),
            "rate_claim_ids": list(scientific_claim.get("rate_claim_ids") or []),
            "path_claim_ids": list(scientific_claim.get("path_claim_ids") or []),
            "execution_mode": mode,
        },
        "scope": (
            "Conditional on the retained dark neutral homogeneous implicit-solution "
            + (
                "environment, ideal-activity mass-action coverage, dimensional rate evidence, "
                "diffusion diagnostics and uncertainty."
                if lane.lane_id == "solution-neutral-implicit-mass-action-v1"
                else "environment, restricted unimolecular coverage, rate evidence and uncertainty."
            )
        ),
        "not_a_claim": "No gas-phase, cross-environment, or universal lifetime is inferred.",
    }


class SolutionCampaignExecutor:
    """Execute independent restricted implicit-solution campaign lanes."""

    def __init__(
        self,
        *,
        bundle_runner: Callable[..., dict[str, Any]] | None = None,
        session_runner: Callable[..., dict[str, Any]] | None = None,
    ) -> None:
        self._bundle_runner = bundle_runner
        self._session_runner = session_runner

    def _execute_lane(
        self,
        request: StabilityCampaignRequest,
        lane: CampaignLaneRequest,
        workdir: Path,
        *,
        request_base_dir: Path,
        ncores: int,
    ) -> dict[str, Any]:
        mass_action = lane.lane_id == "solution-neutral-implicit-mass-action-v1"
        if len(request.chemical_system.components) != 1 and not mass_action:
            raise ValueError(
                "Restricted solution execution does not support additional mixture components"
            )
        if (
            request.chemical_system.target.charge != 0
            or request.chemical_system.target.multiplicity != 1
        ):
            raise ValueError("Restricted solution execution requires a neutral closed-shell target")
        evidence = resolve_lane_evidence(lane, base_dir=request_base_dir)
        mode = _evidence_mode(evidence, mass_action=mass_action)
        observables_by_id = {item.observable_id: item for item in request.observables}
        observables = [observables_by_id[item] for item in lane.observable_ids]
        if not mass_action and len(observables) != 1:
            raise ValueError("The linear solution lane requires one t95 observable per request")
        if not mass_action and observables[0].kind != "t95":
            raise ValueError("The linear solution campaign lane supports t95 only")
        observable_parameters_by_id: dict[str, dict[str, Any]] = {}

        session_request: dict[str, Any] | None = None
        direct_bundle: dict[str, Any] | None = None
        direct_budget: dict[str, Any] | None = None
        if mode == "verified_input":
            environment_payload, _ = _load_single(evidence, "solution_environment")
            network_kind = (
                "verified_solution_mass_action_network" if mass_action
                else "verified_solution_network"
            )
            network_payload, network_reference = _load_single(evidence, network_kind)
            if mass_action:
                from .solution_mass_action import (
                    MassActionSolutionEnvironment, VerifiedMassActionNetwork,
                )

                environment = MassActionSolutionEnvironment.from_dict(environment_payload)
                network = VerifiedMassActionNetwork.from_dict(network_payload)
            else:
                environment = ImplicitSolutionEnvironment.from_dict(environment_payload)
                network = VerifiedSolutionNetwork.from_dict(network_payload)
            _validate_environment(environment, lane)
            target_species = next(
                item for item in network.species
                if item.species_id == network.target_species_id
            )
            _validate_target(target_species.smiles, request)
            if mass_action:
                from .reaction_candidates import canonical_smiles

                declared_components = {
                    canonical_smiles(item.smiles) for item in request.chemical_system.components
                }
                positive_initial_species = {
                    canonical_smiles(item.smiles) for item in network.species
                    if item.initial_concentration_mol_L > 0.0
                }
                if not positive_initial_species <= declared_components:
                    raise ValueError(
                        "Positive-concentration mass-action species are absent from the campaign chemical system"
                    )
                component_by_id = {
                    item.component_id: item for item in request.chemical_system.components
                }
                species_by_identity = {
                    canonical_smiles(item.smiles): item.species_id for item in network.species
                }

                def species_for_component(component_id: str) -> str:
                    component = component_by_id[component_id]
                    try:
                        return species_by_identity[canonical_smiles(component.smiles)]
                    except KeyError as error:
                        raise ValueError(
                            f"Observable component {component_id} is absent from the mass-action network"
                        ) from error

                for observable in observables:
                    if observable.kind in {"t95", "retention_at_duration"}:
                        if observable.target_component_id != request.chemical_system.target_component_id:
                            raise ValueError(f"{observable.kind} must target the campaign stability target")
                        parameters = {"species_id": network.target_species_id}
                    elif observable.kind in {"concentration_at_time", "product_yield"}:
                        parameters = {
                            "species_id": species_for_component(observable.target_component_id),
                            "component_id": observable.target_component_id,
                        }
                        if (
                            observable.kind == "product_yield"
                            and observable.target_component_id == request.chemical_system.target_component_id
                        ):
                            raise ValueError("product_yield must target a declared non-target product")
                    elif observable.kind == "dominant_product":
                        if observable.target_component_id != request.chemical_system.target_component_id:
                            raise ValueError("dominant_product must target the campaign stability target")
                        requested_components = tuple(
                            observable.parameters.get("candidate_component_ids") or ()
                        )
                        if not requested_components:
                            requested_components = tuple(
                                item.component_id for item in request.chemical_system.components
                                if item.role == "product"
                            )
                        if not requested_components:
                            raise ValueError(
                                "dominant_product requires declared product components or candidate_component_ids"
                            )
                        unknown = set(requested_components) - set(component_by_id)
                        if unknown:
                            raise ValueError("dominant_product references unknown candidate components")
                        species_to_component = {
                            species_for_component(component_id): component_id
                            for component_id in requested_components
                        }
                        parameters = {
                            "candidate_species_ids": list(species_to_component),
                            "species_to_component_id": species_to_component,
                        }
                    elif observable.kind in {"branching_fraction", "pathway_flux"}:
                        reaction_id = str(observable.parameters.get("reaction_id") or "").strip()
                        if not reaction_id:
                            raise ValueError(f"{observable.kind} requires parameters.reaction_id")
                        parameters = {"reaction_id": reaction_id}
                        if observable.kind == "pathway_flux":
                            parameters["metric"] = str(
                                observable.parameters.get("metric", "integrated_extent")
                            )
                    else:
                        raise ValueError(
                            f"Unsupported mass-action campaign observable: {observable.kind}"
                        )
                    observable_parameters_by_id[observable.observable_id] = parameters
            if (
                network_reference.get("claim_id") is not None
                and network_reference["claim_id"] != network.coverage_claim_id
            ):
                raise ValueError(
                    "Verified solution network evidence claim_id does not match its "
                    "coverage_claim_id"
                )
            direct_budget = _validated_direct_budget(lane, mass_action=mass_action)
            direct_bundle = {
                "schema_version": 1,
                "kind": (
                    "verified_implicit_solution_mass_action_stability_bundle"
                    if mass_action else "verified_implicit_solution_stability_bundle"
                ),
                "target_smiles": request.chemical_system.target.smiles,
                "environment": environment_payload,
                "network": network_payload,
                "producer": {
                    "kind": "storca_stability_campaign_adapter",
                    "campaign_contract_id": request.contract_id,
                },
            }
        else:
            session_request, _ = _load_single(evidence, SESSION_INPUT_KIND)
            environment = ImplicitSolutionEnvironment.from_dict(
                session_request.get("environment") or {},
            )
            _validate_environment(environment, lane)
            _validate_target(str(session_request.get("target_smiles") or ""), request)
            if session_request.get("observable", "t95") != "t95":
                raise ValueError("Solution session observable must be t95")
            if session_request.get("method_profile") != lane.method_profile:
                raise ValueError("Solution session method does not match the campaign lane")
            _validated_session_budget(lane, session_request)

        if mode == "verification_session" and self._session_runner is None:
            from src.orca_runner import find_orca

            software_resolution = {"orca_executable": find_orca()}
        elif mode == "verification_session":
            software_resolution = {
                "orca_executable": "injected_session_runner_for_test_or_embedding",
            }
        else:
            software_resolution = {"orca_executable": "not_required_for_verified_input"}

        execution_payload = {
            "schema_version": SOLUTION_EXECUTION_SCHEMA_VERSION,
            "kind": "storca_solution_campaign_execution_contract",
            "campaign_contract_id": request.contract_id,
            "lane_request_id": lane.lane_request_id,
            "lane_id": lane.lane_id,
            "environment_contract_id": lane.environment.as_dict()["contract_id"],
            "method_profile": lane.method_profile,
            "observable_ids": list(lane.observable_ids),
            "observable_id": lane.observable_ids[0] if len(lane.observable_ids) == 1 else None,
            "execution_mode": mode,
            "compute_budget": dict(lane.compute_budget),
            "ncores": ncores if mode == "verification_session" else None,
            "software_resolution": software_resolution,
            "evidence": evidence_contract_records(evidence),
        }
        contract_id = _canonical_hash(execution_payload)
        contract = {**execution_payload, "contract_id": contract_id}
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        contract_path = workdir / "solution-execution-contract.json"
        result_path = workdir / "solution-campaign-result.json"
        if contract_path.is_file():
            retained = json.loads(contract_path.read_text())
            if retained.get("contract_id") != contract_id:
                raise ValueError(
                    "Solution execution directory contains a different execution contract"
                )
        else:
            _write_json(contract_path, contract)
        if result_path.is_file():
            retained_result = json.loads(result_path.read_text())
            if (
                retained_result.get("execution_contract_id") == contract_id
                and retained_result.get("terminal") is True
                and retained_result.get("status") == "completed"
            ):
                return {**retained_result, "resumed": True}

        try:
            if mode == "verified_input":
                if self._bundle_runner is None:
                    if mass_action:
                        from .solution_mass_action_stability import (
                            run_implicit_solution_mass_action_bundle,
                        )

                        runner = run_implicit_solution_mass_action_bundle
                    else:
                        from .solution_stability import run_implicit_solution_stability_bundle

                        runner = run_implicit_solution_stability_bundle
                else:
                    runner = self._bundle_runner
                native_results = []
                lane_results = []
                scientific_claims = []
                for observable in observables:
                    observable_dir = (
                        workdir / "verified-input"
                        if len(observables) == 1
                        else workdir / "verified-input" / observable.observable_id
                    )
                    native_result = runner(
                        direct_bundle, observable_dir, **(direct_budget or {}),
                        **(
                            {
                                "requested_observable": observable.kind,
                                "observable_parameters": observable_parameters_by_id[
                                    observable.observable_id
                                ],
                            }
                            if mass_action else {}
                        ),
                    )
                    lane_result = (
                        (native_result.get("lane_run") or {}).get("lane_result") or {}
                    )
                    scientific_claim = (
                        ((native_result.get("lane_run") or {}).get("evidence") or {})
                        .get("scientific_claim") or {}
                    )
                    if not lane_result:
                        lane_result = LaneResult(
                            lane_id=lane.lane_id,
                            lane_contract_id=(
                                build_stability_lane_registry().get(lane.lane_id)
                                .as_dict()["contract_id"]
                            ),
                            admission_status="supported",
                            support_level="observable_supported",
                            execution_status=str(
                                native_result.get("status")
                                or "invalid_or_incomplete_evidence"
                            ),
                            observable_status="unresolved",
                            requested_observable=observable.kind,
                            reason=str(
                                native_result.get("reason")
                                or "Verified-input solution propagation did not return a lane result."
                            ),
                        ).as_dict()
                    native_results.append(native_result)
                    lane_results.append(lane_result)
                    scientific_claims.append(scientific_claim)
            else:
                if self._session_runner is None:
                    from .solution_verification_session import (
                        run_solution_verification_session,
                    )

                    runner = run_solution_verification_session
                else:
                    runner = self._session_runner
                native_result = runner(
                    session_request, workdir / "verification-session",
                    execute=True, ncores=ncores,
                )
                lane_result = _lane_result_from_session(lane, native_result)
                scientific_claim = native_result.get("scientific_claim") or {}
                native_results = [native_result]
                lane_results = [lane_result]
                scientific_claims = [scientific_claim]
            observable_claims = []
            for observable, native_result, lane_result, scientific_claim in zip(
                observables, native_results, lane_results, scientific_claims,
            ):
                claim = _condition_scoped_claim(
                    request, lane, lane_result, scientific_claim,
                    mode=mode, result_json=native_result.get("result_json"),
                    execution_contract_id=contract_id,
                    observable_id=observable.observable_id,
                )
                if claim is not None:
                    observable_claims.append(claim)
            native_result_json = native_results[0].get("result_json")
            terminal = True
            result = {
                "schema_version": SOLUTION_EXECUTION_SCHEMA_VERSION,
                "kind": "storca_solution_campaign_lane_execution",
                "status": "completed",
                "terminal": terminal,
                "campaign_contract_id": request.contract_id,
                "execution_contract_id": contract_id,
                "lane_request_id": lane.lane_request_id,
                "lane_id": lane.lane_id,
                "execution_mode": mode,
                "lane_result": lane_results[0],
                "lane_results": lane_results,
                "scientific_claim": scientific_claims[0],
                "scientific_claims": scientific_claims,
                "observable_claim": observable_claims[0] if len(observable_claims) == 1 else None,
                "observable_claims": observable_claims,
                "native_status": native_results[0].get("status"),
                "native_statuses": [item.get("status") for item in native_results],
                "native_result_json": native_result_json,
                "result_json": str(result_path),
                "resumed": False,
            }
        except Exception as error:
            result = {
                "schema_version": SOLUTION_EXECUTION_SCHEMA_VERSION,
                "kind": "storca_solution_campaign_lane_execution",
                "status": "failed_closed",
                "terminal": True,
                "campaign_contract_id": request.contract_id,
                "execution_contract_id": contract_id,
                "lane_request_id": lane.lane_request_id,
                "lane_id": lane.lane_id,
                "execution_mode": mode,
                "lane_result": None,
                "lane_results": [],
                "observable_claim": None,
                "observable_claims": [],
                "reason": f"{type(error).__name__}: {error}",
                "result_json": str(result_path),
                "resumed": False,
            }
        _write_json(result_path, result)
        return result

    def execute(
        self,
        request: StabilityCampaignRequest,
        solution_lanes: tuple[CampaignLaneRequest, ...],
        workdir: Path,
        *,
        request_base_dir: Path,
        ncores: int = 1,
    ) -> dict[str, Any]:
        if ncores < 1:
            raise ValueError("ncores must be positive")
        if not solution_lanes:
            raise ValueError("SolutionCampaignExecutor needs at least one solution lane")
        maximum_lanes = request.global_compute_budget.get("maximum_solution_lanes")
        if maximum_lanes is not None:
            if int(maximum_lanes) != maximum_lanes or int(maximum_lanes) < 1:
                raise ValueError("maximum_solution_lanes must be a positive integer")
            if len(solution_lanes) > int(maximum_lanes):
                raise ValueError("Requested solution lanes exceed maximum_solution_lanes")
        maximum_core_hours = request.global_compute_budget.get(
            "maximum_solution_estimated_core_hours"
        )
        requested_core_hours = sum(
            float(lane.compute_budget.get("maximum_estimated_core_hours", 0.0))
            for lane in solution_lanes
        )
        if (
            maximum_core_hours is not None
            and requested_core_hours > float(maximum_core_hours)
        ):
            raise ValueError(
                "Requested solution lane budgets exceed "
                "maximum_solution_estimated_core_hours"
            )
        lane_results = []
        for lane in solution_lanes:
            lane_dir = Path(workdir) / lane.lane_request_id
            try:
                lane_result = self._execute_lane(
                    request, lane, lane_dir,
                    request_base_dir=request_base_dir, ncores=ncores,
                )
            except Exception as error:
                lane_result = {
                    "schema_version": SOLUTION_EXECUTION_SCHEMA_VERSION,
                    "kind": "storca_solution_campaign_lane_execution",
                    "status": "failed_closed",
                    "terminal": True,
                    "campaign_contract_id": request.contract_id,
                    "execution_contract_id": None,
                    "lane_request_id": lane.lane_request_id,
                    "lane_id": lane.lane_id,
                    "execution_mode": None,
                    "lane_result": None,
                    "lane_results": [],
                    "observable_claim": None,
                    "observable_claims": [],
                    "reason": f"{type(error).__name__}: {error}",
                    "result_json": str(lane_dir / "solution-campaign-result.json"),
                    "resumed": False,
                }
                _write_json(Path(lane_result["result_json"]), lane_result)
            lane_results.append(lane_result)
        failed = [item for item in lane_results if item["status"] == "failed_closed"]
        return {
            "schema_version": SOLUTION_EXECUTION_SCHEMA_VERSION,
            "kind": "storca_solution_campaign_execution",
            "status": "failed_closed" if len(failed) == len(lane_results) else "completed",
            "lane_executions": lane_results,
            "observable_claims": [
                claim for item in lane_results
                for claim in (
                    item.get("observable_claims")
                    or ([item["observable_claim"]] if item.get("observable_claim") else [])
                )
            ],
            "failed_lane_request_ids": [item["lane_request_id"] for item in failed],
        }


__all__ = [
    "DIRECT_INPUT_KINDS", "MASS_ACTION_DIRECT_INPUT_KINDS", "SESSION_INPUT_KIND",
    "SolutionCampaignExecutor",
]
