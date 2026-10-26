"""Family-neutral iterative route verification state machine.

The engine owns evidence gates, iteration order, persistence, and terminal
claims.  Chemistry-specific path finding, rate construction, and the full RMG
rerun are injected callbacks.  A callback failure never becomes a lifetime;
the engine retains the failed iteration and emits ``verification_incomplete``.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Callable, Iterable, Mapping

from .flux_verification import (
    DEFAULT_FLUX_COVERAGE_TARGET,
    DEFAULT_MAX_FLUX_CLOSURE_ERROR,
    DEFAULT_T95_RELATIVE_TOLERANCE,
    assess_repaired_mechanism,
    assess_t95_robustness,
    assess_target_loss_commitment,
    build_flux_ranked_verification_plan,
    prepare_combined_prediction_rate_bounds,
    prepare_combined_unresolved_rate_bounds,
    reaction_signature,
    run_collision_limit_rate_envelope,
    run_combined_prediction_rate_sensitivity,
    run_combined_unresolved_sensitivity,
    run_t95_sensitivity_attribution,
    verify_applied_replacement,
    _retained_network_routes,
    _solver_gates,
)
from .generated_kinetics import validate_rate_replacement_evidence
from .kinetic_uncertainty import (
    RateUncertaintyRecord,
    coerce_rate_uncertainty_evidence,
    rate_uncertainty_artifact,
)
from .surface_resolution import SURFACE_TERMINAL_STATES


VerifyRoute = Callable[[dict, Path], dict]
BuildRate = Callable[[dict, dict, Path], dict]
RerunModel = Callable[[list[Path], Path], dict]
BuildSpecialistRate = Callable[[dict, Path], dict]

_PATH_VERIFIED_STATUSES = {
    "verified", "orca_verified", "irc_endpoint_verified", "barrierless_capture_verified",
    "barriered_ts_verified", "verified_barrierless_capture", "verified_irc_endpoints",
    "verified_transition_state_path", "verified_multistep_path",
}
_RATE_VERIFIED_STATUSES = {
    "validated", "verified", "rate_verified", "arkane_verified", "collision_bounded_verified",
    "verified_rate", "verified_arkane", "verified_collision_bounded",
}
_DEFERABLE_MAPPING_STATUSES = {
    # These are explicit fail-closed outcomes from resolve_route_atom_mapping.
    # They invalidate one candidate route, not the verification engine itself;
    # the engine must retain the evidence and continue to the next ranked route.
    "all_atom_element_mismatch",
    "source_atom_mapping_evidence_invalid",
    "ambiguous_atom_mapping",
    "ambiguous_hydrogen_transfer_mapping",
    "atom_mapping_unresolved",
    "declared_multiplicity_not_common",
    "declared_total_charge_mismatch",
    "fragment_charge_contract_incomplete",
    "fragment_spin_contract_incomplete",
    "incomplete_all_atom_mapping",
    "invalid_endpoint",
    "mapping_search_exhausted",
    "mapping_symmetry_class_exhausted",
    "multiple_common_spin_surfaces_require_selection",
    "no_common_adiabatic_spin_surface",
    "unbalanced_elements",
    "unbalanced_explicit_hydrogens",
    "unbalanced_hydrogen_isotopes",
    "unbalanced_smiles_atom_count",
    "unbalanced_total_charge",
}
_DEFERABLE_PATH_STATUSES = {
    *_DEFERABLE_MAPPING_STATUSES,
    *SURFACE_TERMINAL_STATES,
    "barrierless_rate_model_required",
    "endpoint_validation_failed",
    "execution_infrastructure_failure",
    "intermediate_detected_requires_segmented_verification",
    "irc_endpoint_validation_failed",
    "mapping_unresolved",
    "path_execution_failed",
    "stationary_point_validation_failed",
    "surface_unresolved",
    "thermochemistry_bounded_below_loss_threshold",
    "ts_frequency_validation_failed",
}


@dataclass(frozen=True)
class VerificationEngineConfig:
    """Explicit non-chemical stopping and numerical policy."""

    flux_coverage_target: float = DEFAULT_FLUX_COVERAGE_TARGET
    max_flux_closure_error: float = DEFAULT_MAX_FLUX_CLOSURE_ERROR
    t95_relative_tolerance: float = DEFAULT_T95_RELATIVE_TOLERANCE
    max_iterations: int = 12

    def __post_init__(self) -> None:
        if not 0.0 < self.flux_coverage_target <= 1.0:
            raise ValueError("Flux coverage target must be in (0, 1]")
        if not 0.0 <= self.max_flux_closure_error < 1.0:
            raise ValueError("Maximum flux closure error must be in [0, 1)")
        if self.t95_relative_tolerance < 0.0:
            raise ValueError("t95 relative tolerance must be nonnegative")
        if self.max_iterations < 1:
            raise ValueError("Verification engine needs at least one iteration")


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _write_json(path: Path, payload: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n")
    return str(path)


def _condition_hash(condition_contract: dict) -> str:
    encoded = json.dumps(_jsonable(condition_contract), sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _explicit_gas_environment(environment: Mapping[str, Any]) -> bool:
    """Reject every conflicting/unsupported declaration, not just the first."""
    allowed = {
        "gas", "gas phase", "homogeneous gas", "homogeneous gas phase",
        "homogeneous gas phase surrogate", "homogeneous gas phase low pressure surrogate",
        "homogeneous gas phase pure target model",
    }
    declarations = [
        " ".join(str(environment[key]).strip().casefold().replace("_", " ").replace("-", " ").split())
        for key in ("phase", "phase_model", "model", "environment_kind") if key in environment
    ]
    return bool(declarations) and all(value in allowed for value in declarations) and not any(
        environment.get(key) not in (None, "") for key in ("solvent", "solvent_model")
    )


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _resolve_retained_path(raw_path: object, *, state_path: Path) -> Path:
    """Resolve an artifact path without guessing outside the retained run."""
    path = Path(str(raw_path))
    if path.is_absolute():
        return path
    if path.exists():
        return path.resolve()
    run_dir = state_path.resolve().parent.parent
    suffix = None
    if run_dir.name in path.parts:
        marker = len(path.parts) - 1 - tuple(reversed(path.parts)).index(run_dir.name)
        suffix = Path(*path.parts[marker + 1:])
    candidates = tuple(item for item in (
        run_dir / path,
        run_dir.parent / path,
        run_dir / suffix if suffix is not None else None,
    ) if item is not None)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return path


def _unique_session(root: Path) -> Path:
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    base = root / "verification-engine"
    if not base.exists():
        base.mkdir()
        return base
    index = 2
    while True:
        candidate = root / f"verification-engine-{index:03d}"
        if not candidate.exists():
            candidate.mkdir()
            return candidate
        index += 1


def _finite(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _extract_model_result(result: dict) -> tuple[dict, dict | None]:
    if not isinstance(result, dict):
        raise TypeError("Model callback must return a dictionary")
    for key in ("rmg_evidence", "repaired_rmg_evidence", "model_evidence", "evidence"):
        if isinstance(result.get(key), dict):
            return result[key], result.get("condition_contract") or result[key].get("condition_contract")
    return result, result.get("condition_contract")


def _sensitivity_envelope(result: dict, evidence: dict) -> dict | None:
    for source in (result, evidence):
        for key in ("unresolved_sensitivity_envelope", "combined_unresolved_sensitivity",
                    "combined_unresolved_route_sensitivity",
                    "combined_prediction_rate_sensitivity"):
            if isinstance(source.get(key), dict):
                return source[key]
    return None


def _nested_dict(result: dict, keys: Iterable[str]) -> dict | None:
    for key in keys:
        value = result.get(key)
        if isinstance(value, dict):
            return value
    return None


def _verified_status(result: dict, accepted: set[str], *, nested_keys: Iterable[str],
                     allow_verified_prefix: bool = False) -> str | None:
    blocks = [result]
    nested = _nested_dict(result, nested_keys)
    if nested is not None:
        blocks.insert(0, nested)
    for block in blocks:
        for key in ("route_verification_status", "rate_verification_status", "verification_status", "status"):
            status = str(block.get(key) or "")
            if status in accepted or (allow_verified_prefix and status.startswith("verified_")):
                return status
    return None


def _unverified_path_detail(result: dict) -> str:
    """Expose the chemistry-level stop instead of only the adapter contract."""
    block = _nested_dict(result, (
        "path_result", "route_verification", "verification", "orca_verification",
    )) or result
    status = str(block.get("status") or "status_missing")
    reasons = []
    if block.get("failure_reason"):
        reasons.append(str(block["failure_reason"]))
    orientations = result.get("orientations") or block.get("orientations") or []
    for item in orientations:
        if item.get("failure_reason"):
            detail = str(item["failure_reason"])
            if detail not in reasons:
                reasons.append(detail)
    suffix = "; ".join(reasons[:3])
    return f"route_status={status}" + (f"; {suffix}" if suffix else "")


def _unverified_path_status(result: dict) -> str:
    block = _nested_dict(result, (
        "path_result", "route_verification", "verification", "orca_verification",
    )) or result
    return str(block.get("status") or "status_missing")


def _route_for_plan(evidence: dict, plan: dict, condition_contract: dict | None = None) -> dict:
    source = plan.get("selected_route_source") or "candidate_routes"
    selected = plan.get("selected_route")
    if isinstance(selected, dict) and source != "retained_partial_routes":
        route = dict(selected)
        if not route.get("route_id"):
            raise ValueError("Supplemental verification plan selected a route without an ID")
        if plan.get("selected_route_id") and str(route["route_id"]) != str(plan["selected_route_id"]):
            raise ValueError("Supplemental verification plan route and route ID disagree")
        from .rmg_mapping import bind_rmg_mapping_source
        return bind_rmg_mapping_source(route, evidence, condition_contract)
    if source == "retained_partial_routes":
        from .rmg_partial_candidates import validated_retained_partial_routes
        retained = validated_retained_partial_routes(evidence, condition_contract or {})
        if retained["blocking_reasons"]:
            raise ValueError("Retained attempt evidence rejected: " + ", ".join(retained["blocking_reasons"]))
        routes = retained["routes"]
    elif source == "candidate_routes":
        routes = evidence.get("candidate_routes") or []
    elif source == "network_routes":
        routes = _retained_network_routes(evidence)
    else:
        raise ValueError(f"Unknown route source: {source}")
    index = plan.get("selected_route_index")
    if not isinstance(index, int) or index < 0 or index >= len(routes):
        raise ValueError("Flux plan selected no resolvable candidate route")
    route = dict(routes[index])
    route.setdefault("route_id", plan.get("selected_route_id"))
    if str(route.get("route_id")) != str(plan.get("selected_route_id")):
        raise ValueError("Flux plan route index and route ID disagree")
    from .rmg_mapping import bind_rmg_mapping_source
    return bind_rmg_mapping_source(route, evidence, condition_contract)


def _normalize_replacement(rate_result: dict, route: dict) -> tuple[dict, list[Path]]:
    if not isinstance(rate_result, dict):
        raise TypeError("Rate callback must return a dictionary")
    block = _nested_dict(rate_result, (
        "replacement", "rate_evidence", "verified_rate", "generated_kinetics_library", "kinetics",
    )) or rate_result
    replacement = dict(block)
    replacement.setdefault("route_id", route.get("route_id"))
    if route.get("native_candidate_id"):
        replacement.setdefault("native_candidate_id", route.get("native_candidate_id"))
        replacement.setdefault(
            "discovery_sources", list(route.get("discovery_sources") or ()),
        )
    replacement.setdefault("reaction_equation", route.get("reaction_equation") or route.get("printed_reaction_equation"))
    source = str(replacement.get("source") or "")
    if source.startswith("arkane_") or rate_result.get("kind") == "storca_specialist_rate_result":
        replacement["rate_uncertainty_propagation_required"] = True
    if source == "arkane_pdep":
        replacement["required_uncertainty_components"] = [
            "electronic_structure", "conformer_sampling",
            "harmonic_statistical_mechanics", "tunneling",
            "energy_transfer", "transport", "pressure_dependence_model",
        ]
    elif source == "mess_pdep":
        replacement["required_uncertainty_components"] = [
            "electronic_structure", "conformer_sampling", "statistical_mechanics",
            "capture_or_transition_state", "energy_transfer",
            "collision_transport", "master_equation_numerics",
        ]
    elif source == "arkane_tst":
        replacement["required_uncertainty_components"] = [
            "electronic_structure", "conformer_sampling",
            "harmonic_statistical_mechanics", "tunneling",
        ]
    explicit_status = _verified_status(rate_result, _RATE_VERIFIED_STATUSES, nested_keys=(
        "replacement", "rate_evidence", "verified_rate", "generated_kinetics_library", "kinetics",
    ), allow_verified_prefix=True)
    if not replacement.get("verification_status") and explicit_status:
        replacement["verification_status"] = (
            explicit_status if explicit_status.startswith("verified") else f"verified_{explicit_status}"
        )
    if not replacement.get("condition_rate"):
        nested_rate = _nested_dict(rate_result, ("condition_rate", "rate_constant"))
        if nested_rate:
            replacement["condition_rate"] = nested_rate

    raw_paths: list[Any] = []
    for source in (rate_result, block):
        for key in ("library_path", "library"):
            value = source.get(key)
            if value:
                raw_paths.append(value)
        value = source.get("library_paths")
        if isinstance(value, (list, tuple)):
            raw_paths.extend(value)
    paths = []
    for value in raw_paths:
        path = Path(value)
        if path not in paths:
            paths.append(path)
    if not paths:
        raise ValueError("Verified rate callback did not return a replacement library path")
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Replacement library path does not exist: " + ", ".join(missing))
    return replacement, paths


def _planning_replacements(
    replacements: list[dict], evidence: dict, condition_contract: dict,
) -> list[dict]:
    """Rebind only rates proved present in this model, without changing IDs.

    The source replacements remain immutable scientific records. These are
    planning-only views in the current model's label space, bound to the same
    structural correspondence and numerical checks used for repair acceptance.
    """
    candidate_routes = list(evidence.get("candidate_routes") or [])
    routes = list(candidate_routes)
    candidate_signatures = {
        reaction_signature(
            route.get("reaction_equation") or route.get("printed_reaction_equation") or "",
            normalized=True,
        )
        for route in routes
    }
    routes.extend(
        route for route in _retained_network_routes(evidence)
        if reaction_signature(
            route.get("reaction_equation") or route.get("equation") or "", normalized=True,
        )
        not in candidate_signatures
    )
    rebound = []
    for replacement in replacements:
        item = dict(replacement)
        if "forward_only_approximation_evidence" in item:
            # Network insertion under a finite-horizon omission certificate
            # does not establish physical rate uncertainty or applicability.
            # Never allow this central-model approximation to earn the same
            # observable-convergence credit as a fully admitted verified rate.
            item["verification_status"] = "conditional_model_rate_not_observable_verified"
            item["planning_applied_rate_validation"] = {
                "status": "failed",
                "blocking_reasons": ["conditional_model_rate_not_observable_verified"],
            }
            rebound.append(item)
            continue
        applied = verify_applied_replacement(evidence, replacement, condition_contract)
        item["planning_applied_rate_validation"] = applied
        identity = applied.get("identity_resolution") or {}
        item["identity_source_route_id"] = replacement.get("route_id")
        item["source_reaction_equation"] = (
            replacement.get("reaction_equation") or replacement.get("replaces_reaction")
        )
        if applied.get("status") != "passed":
            rebound.append(item)
            continue
        if identity.get("status") == "matched":
            # Graph matching has already been proved by the shared resolver.
            # Here we only locate that exact reaction in the model's route
            # views; normalized formula-like labels are not identity evidence.
            model_signature = reaction_signature(identity.get("source_oriented_equation") or "")
            reaction_index = applied.get("matched_reaction_index")

            def matches_model_route(route: dict) -> bool:
                equation = route.get("reaction_equation") or route.get("printed_reaction_equation") or route.get("equation") or ""
                if model_signature is None or reaction_signature(equation) != model_signature:
                    return False
                retained = route.get("retained_rate_evidence") or {}
                index = retained.get("reaction_index", route.get("reaction_index"))
                return index is None or index == reaction_index

            matches = [route for route in candidate_routes if matches_model_route(route)]
            if not matches:
                matches = [route for route in _retained_network_routes(evidence) if matches_model_route(route)]
            item["reaction_equation"] = identity.get("source_oriented_equation")
        else:
            # Historical-only compatibility: preserve the original normalized
            # label match, after proving the applied rate. Failed structural
            # evidence cannot reach this branch.
            signature = reaction_signature(item.get("reaction_equation") or item.get("replaces_reaction") or "", normalized=True)
            matches = [route for route in routes if reaction_signature(
                route.get("reaction_equation") or route.get("printed_reaction_equation") or "", normalized=True,
            ) == signature]
        if len(matches) == 1 and matches[0].get("route_id"):
            item["route_id"] = matches[0].get("route_id")
        else:
            item["planning_applied_rate_validation"] = {
                **applied,
                "status": "failed",
                "blocking_reasons": [
                    *(applied.get("blocking_reasons") or []),
                    "planning_model_route_missing_or_ambiguous",
                ],
            }
        rebound.append(item)
    return rebound


def _propagation_loss(propagation: dict) -> float | None:
    value = _finite(propagation.get("target_loss_fraction"))
    if value is not None:
        return value
    profile = propagation.get("target_profile") or []
    final = _finite(profile[-1].get("target_fraction_remaining")) if profile else None
    return max(0.0, 1.0 - final) if final is not None else None


def _model_loss(evidence: dict) -> float | None:
    return _propagation_loss(evidence.get("independent_cantera_propagation") or {})


def _unresolved_loss_upper_bound(evidence: dict, envelope: dict | None) -> float:
    baseline = _model_loss(evidence) or 0.0
    if not envelope:
        return 0.0
    maximum = _finite(envelope.get("maximum_loss_upper_bound_fraction"))
    if maximum is not None:
        return max(0.0, maximum - baseline)
    losses = [baseline]
    for key in ("lower_bound_propagation", "upper_bound_propagation"):
        value = _propagation_loss(envelope.get(key) or {})
        if value is not None:
            losses.append(value)
    return max(0.0, max(losses) - baseline)


class VerificationEngine:
    """Execute one verified route/repaired-model iteration at a time."""

    def __init__(
        self,
        initial_rmg_evidence: dict,
        condition_contract: dict,
        workdir: Path,
        *,
        verify_route: VerifyRoute,
        build_rate: BuildRate,
        rerun_model: RerunModel,
        config: VerificationEngineConfig | None = None,
        supplemental_candidate_routes: Iterable[dict] | None = None,
        build_specialist_rate: BuildSpecialistRate | None = None,
        rate_uncertainties: Iterable[dict | RateUncertaintyRecord] = (),
        rmg_env: str | None = None,
        initial_replacements: Iterable[dict] = (),
        initial_library_paths: Iterable[Path] = (),
        initial_deferred_routes: dict[str, dict] | None = None,
        calculation_workspaces: dict[str, Path] | None = None,
        path_verifier_policy_id: str = "generic-orca-path-v1",
        route_policy_attempts: dict[str, Iterable[str]] | None = None,
        same_policy_execution_recovery_route_ids: Iterable[str] = (),
        retained_specialist_path_results: Mapping[str, Mapping[str, Any]] | None = None,
        continuation_evidence: dict | None = None,
    ) -> None:
        if not isinstance(condition_contract, dict) or not condition_contract:
            raise ValueError("Verification requires a complete immutable condition contract")
        retention = _finite(condition_contract.get("retention_fraction"))
        if retention is None or not 0.0 < retention < 1.0:
            raise ValueError("Condition contract retention_fraction must be in (0, 1)")
        self.initial_result = initial_rmg_evidence
        self.condition_contract = _jsonable(condition_contract)
        self.condition_sha256 = _condition_hash(self.condition_contract)
        self.rate_uncertainties = coerce_rate_uncertainty_evidence(rate_uncertainties)
        self.rmg_env = rmg_env
        self.session = _unique_session(Path(workdir))
        self.verify_route_callback = verify_route
        self.build_rate_callback = build_rate
        self.rerun_model_callback = rerun_model
        self.build_specialist_rate_callback = build_specialist_rate
        self.config = config or VerificationEngineConfig()
        self.replacements = [dict(item) for item in initial_replacements]
        self.library_paths = [Path(item) for item in initial_library_paths]
        self.iterations: list[dict] = []
        # Calculation evidence is not the same claim as admission to a repaired
        # kinetic model.  Keep rejected handoffs outside replacements/libraries
        # so they cannot receive verified flux credit or trigger propagation.
        self.conditional_rate_records: list[dict] = []
        self.deferred_routes = {
            str(route_id): dict(item)
            for route_id, item in (initial_deferred_routes or {}).items()
        }
        self.calculation_workspaces = {
            str(route_id): Path(path)
            for route_id, path in (calculation_workspaces or {}).items()
        }
        self.path_verifier_policy_id = str(path_verifier_policy_id).strip()
        if not self.path_verifier_policy_id:
            raise ValueError("Path-verifier policy ID cannot be empty")
        self.route_policy_attempts = {
            str(route_id): [str(policy) for policy in policies]
            for route_id, policies in (route_policy_attempts or {}).items()
        }
        self.same_policy_execution_recovery_route_ids = {
            str(route_id) for route_id in same_policy_execution_recovery_route_ids
        }
        self.retained_specialist_path_results = {
            str(route_id): dict(record)
            for route_id, record in (retained_specialist_path_results or {}).items()
        }
        self.continuation_evidence = _jsonable(continuation_evidence)
        # These routes are deliberately not inserted into the RMG flux table:
        # before a verified rate is reinjected their kinetic relevance is
        # unknown.  They are a bounded supplemental queue, not hidden flux.
        self.supplemental_candidate_routes = tuple(
            dict(route) for route in (supplemental_candidate_routes or ())
            if isinstance(route, dict)
        )

    def _supplemental_plan(self, *, require_direct_reactants: bool = False,
                           require_independent_gas: bool = False) -> dict:
        """Select one retained discovery hypothesis without fabricating flux."""
        from .reaction_discovery import graph_candidate_verification_plan

        verified = {str(item.get("route_id")) for item in self.replacements if item.get("route_id")}
        ready = [
            route for route in self.supplemental_candidate_routes
            if route.get("route_id")
            and str(route["route_id"]) not in self.deferred_routes
            and str(route["route_id"]) not in verified
        ]
        if require_independent_gas:
            # Missing RMG solver evidence removes flux ranking, not an
            # independently screened structural hypothesis. This existing
            # gas engine must not adopt a solvent/surface candidate or treat
            # RMG-only discovery as independent evidence.
            condition_id = self.condition_contract.get("contract_id")
            independent = []
            for route in ready if condition_id and _explicit_gas_environment(self.condition_contract) else ():
                sources = set(route.get("discovery_sources") or str(route.get("source") or "").split("+"))
                environment = route.get("environment") or {}
                if (
                    route.get("native_candidate_id") and "graph_xtb" in sources and "rmg" not in sources
                    and route.get("verification_condition_contract_id") == condition_id
                    and _explicit_gas_environment(environment)
                ):
                    independent.append(route)
            ready = independent
        if require_direct_reactants:
            from .flux_verification import normalize_rmg_label

            target = normalize_rmg_label(str(self.condition_contract.get("target_label") or "stability"))
            initial = {
                normalize_rmg_label(str(label))
                for label, amount in (self.condition_contract.get("composition") or {}).items()
                if (_finite(amount) or 0.0) > 0.0
            }
            direct = []
            for route in ready:
                endpoints = route.get("resolved_endpoints") or {}
                reactants = endpoints.get("reactants") or []
                products = endpoints.get("products") or []
                labels = [normalize_rmg_label(str(item.get("label") or "")) for item in reactants]
                route_condition = route.get("verification_condition_contract_id")
                if (
                    labels and products and target in labels
                    and set(labels) <= initial
                    and (route.get("atom_balance") or {}).get("status") == "passed"
                    and route.get("screening_status") == "candidate_unverified"
                    and (not route_condition or route_condition == self.condition_contract.get("contract_id"))
                ):
                    direct.append(route)
            ready = direct
        return graph_candidate_verification_plan(ready)

    def _retain_unexpected_candidates(self, parent: dict, path_result: dict) -> dict:
        """Queue bounded, trace-backed hypotheses; never admit a path or rate."""
        from collections import Counter
        from .reaction_candidates import ReactionCandidate, canonical_smiles

        report = {"kind": "unexpected_candidate_queue_handoff", "schema_version": 1,
                  "condition_sha256": self.condition_sha256, "records": [],
                  "claims_prohibited": ["verified_path", "rate", "lifetime"]}
        pending = [path_result]
        proposals = []
        visited = 0
        while pending:
            node = pending.pop()
            visited += 1
            if visited > 10000:
                raise ValueError("Unexpected-candidate evidence exceeds bounded traversal")
            if isinstance(node, dict):
                proposal = node.get("unexpected_endpoint_discovery")
                if isinstance(proposal, dict):
                    proposals.append(proposal)
                pending.extend(value for key, value in node.items()
                               if key != "unexpected_endpoint_discovery" and isinstance(value, (dict, list)))
            elif isinstance(node, list):
                pending.extend(value for value in node if isinstance(value, (dict, list)))
        known = {str(route.get("route_id")) for route in self.supplemental_candidate_routes}
        known.update(str(route.get("route_id")) for route in self.replacements)
        parent_id = str(parent.get("route_id") or "")
        for proposal in proposals:
            for raw in proposal.get("native_candidates") or []:
                record = {"parent_route_id": parent_id, "status": "not_queued"}
                report["records"].append(record)
                try:
                    candidate = ReactionCandidate.from_dict(raw)
                    record["candidate_id"] = candidate.candidate_id
                    if candidate.candidate_id in known or candidate.candidate_id == parent_id:
                        record["status"] = "already_retained"
                        continue
                    if int(parent.get("unexpected_discovery_depth") or 0) >= 1:
                        raise ValueError("Unexpected discovery depth bound reached")
                    if len(self.supplemental_candidate_routes) >= self.config.max_iterations:
                        raise ValueError("Supplemental candidate budget reached")
                    if proposal.get("status") != "candidate_proposed" or not (
                        proposal.get("trace_atom_mapping") or {}).get("valid"):
                        raise ValueError("Completed trace atom correspondence is unresolved")
                    if candidate.environment != (parent.get("environment") or {}):
                        raise ValueError("Candidate environment differs from parent route")
                    parent_reactants = parent.get("reactant_smiles") or [
                        item.get("smiles") for item in
                        (parent.get("resolved_endpoints") or {}).get("reactants", [])
                    ]
                    if (not parent_reactants or not all(parent_reactants) or
                            Counter(canonical_smiles(value) for value in parent_reactants)
                            != Counter(canonical_smiles(value) for value in candidate.reactants)):
                        raise ValueError("Candidate reactants differ from verified parent anchor")
                    condition_id = parent.get("verification_condition_contract_id")
                    if condition_id and condition_id != self.condition_contract.get("contract_id"):
                        raise ValueError("Parent route condition differs from engine")
                    provenance = candidate.provenance
                    if (provenance.get("parent_route_id") != parent_id or
                            provenance.get("source") != "unexpected_completed_irc_endpoint" or
                            provenance.get("parent_rate_reuse_allowed") is not False or
                            candidate.verification_status != "unexpected_channel_unverified" or
                            candidate.screening_status != "candidate_unverified" or
                            candidate.estimated_barrier_kcal_mol is not None):
                        raise ValueError("Candidate lacks unverified parent-trace contract")
                    for field in ("source_xyz", "source_trajectory"):
                        path = Path(str(provenance.get(field) or ""))
                        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != provenance.get(field + "_sha256"):
                            raise ValueError("Retained discovery geometry hash mismatch")
                    parent_smiles = parent.get("parent_smiles") or (path_result.get("route") or {}).get("parent_smiles")
                    if not parent_smiles:
                        raise ValueError("Parent chemical identity missing")
                    route = candidate.as_route_dict(parent_smiles=parent_smiles)
                    route.update(
                        unexpected_discovery_depth=1,
                        verification_condition_contract_id=self.condition_contract.get("contract_id"),
                        unexpected_discovery_evidence={
                            "parent_route_id": parent_id, "condition_sha256": self.condition_sha256,
                            "proposal": proposal, "retained_ts_reuse_admitted": False,
                        },
                    )
                    self.supplemental_candidate_routes += (route,)
                    known.add(candidate.candidate_id)
                    record["status"] = "queued_unverified"
                except (ValueError, TypeError, OSError) as error:
                    record["reason"] = str(error)
        return report

    def _deferred_supplemental_route_ids(self) -> list[str]:
        retained = {str(route.get("route_id")) for route in self.supplemental_candidate_routes}
        return sorted(route_id for route_id in self.deferred_routes if route_id in retained)

    def _persist_state(self, status: str, **extra: Any) -> None:
        _write_json(self.session / "engine-state.json", {
            "schema_version": 1,
            "kind": "generic_verification_engine_state",
            "status": status,
            "condition_contract": self.condition_contract,
            "condition_contract_sha256": self.condition_sha256,
            "verified_replacements": self.replacements,
            "conditional_rate_records": self.conditional_rate_records,
            "library_paths": self.library_paths,
            "iterations": self.iterations,
            "deferred_unverified_routes": self.deferred_routes,
            "supplemental_candidate_routes": self.supplemental_candidate_routes,
            "specialist_rate_provider_configured": bool(self.build_specialist_rate_callback),
            "rate_uncertainty_evidence": rate_uncertainty_artifact(self.rate_uncertainties),
            "rmg_environment": self.rmg_env,
            "engine_config": _jsonable(self.config.__dict__),
            "path_verifier_policy_id": self.path_verifier_policy_id,
            "route_policy_attempts": self.route_policy_attempts,
            "same_policy_execution_recovery_route_ids": sorted(
                self.same_policy_execution_recovery_route_ids
            ),
            "retained_specialist_path_route_ids": sorted(
                self.retained_specialist_path_results
            ),
            "continuation_evidence": self.continuation_evidence,
            **extra,
        })

    def _retain_blocked_rate_handoff(
        self, *, route_id: str, path_status: str, rate_result: dict,
        replacement: dict, validation: dict, library_paths: list[Path],
        iteration: dict, network_validation: dict | None = None,
        repropagated: bool = False,
    ) -> dict:
        """Retain a calculated elementary rate without admitting a model repair.

        The conditional forward check only separates the conformer/network
        compatibility gate from the pre-existing rate checks.  It is not a
        pressure-applicability or uncertainty certification and cannot alter
        the authoritative replacement validation supplied to this method.
        """
        reasons = list((network_validation or validation).get("blocking_reasons") or ())
        forward_blockers = [reason for reason in validation.get("blocking_reasons") or ()
                            if reason != "conformer_thermo_basis_unresolved"]
        producer_forward = (rate_result.get("conditional_forward_rate_evidence")
                            or replacement.get("conditional_forward_rate_evidence"))
        record = {
            "schema_version": 1,
            "kind": "conditional_elementary_rate_handoff",
            "route_id": replacement.get("route_id") or route_id,
            "parent_route_id": route_id,
            "condition_contract_sha256": self.condition_sha256,
            "path_verification_status": path_status,
            "stationary_point_verification_status": iteration.get(
                "stationary_point_verification_status", "not_reported"),
            "status": ("conditional_model_rate_calculated" if not forward_blockers
                       else "calculated_rate_not_validated"),
            "condition_rate": replacement.get("condition_rate"),
            "source": replacement.get("source"),
            "verification_status": replacement.get("verification_status"),
            "forward_rate_validation": {
                "status": "passed" if not forward_blockers else "failed",
                "blocking_reasons": forward_blockers,
                "scope": "existing_elementary_rate_gates_excluding_conformer_network_compatibility",
                "pressure_applicability_established": False,
                "calibrated_rate_uncertainty_established": False,
                "interpretation": (
                    "A calculated rate within its declared physical model is not a "
                    "condition-specific network rate or a lifetime. Pressure applicability "
                    "and uncertainty require their own evidence."
                ),
            },
            "replacement_validation": validation,
            "network_admission_validation": network_validation or validation,
            "network_admission": {
                "status": "blocked", "blocking_reasons": reasons,
                "network_reinjection_allowed": False,
                "verified_flux_credit": False, "repropagated": repropagated,
                "observable_claim_allowed": False,
            },
            "artifacts": {
                "path_result": iteration["artifacts"].get("path_result"),
                "rate_result": iteration["artifacts"].get("rate_result"),
                "replacement_validation": iteration["artifacts"].get("replacement_validation"),
                "unadmitted_library_paths": [str(path) for path in library_paths],
            },
        }
        for key in ("conformer_rate_basis", "reverse_relevance_assessment",
                    "forward_only_approximation_evidence", "method_applicability_evidence"):
            value = replacement.get(key, rate_result.get(key, iteration.get(key)))
            if value is not None:
                record[key] = value
        if producer_forward is not None:
            record["conditional_forward_rate_evidence"] = producer_forward
        self.conditional_rate_records.append(record)
        iteration.setdefault("conditional_rate_route_ids", []).append(record["route_id"])
        iteration["network_admission"] = record["network_admission"]
        iteration["rate_verification_status"] = record["status"]
        iteration["artifacts"]["conditional_rate_handoffs"] = _write_json(
            Path(iteration["artifacts"]["directory"]) / "conditional-rate-handoffs.json",
            [item for item in self.conditional_rate_records
             if item["parent_route_id"] == route_id],
        )
        return record

    def _try_bounded_forward_only_handoff(
        self, *, replacement: dict, library_paths: list[Path], route: dict,
        current_evidence: dict, iteration: dict,
    ) -> tuple[dict, list[Path]]:
        """Try a scoped, evidence-only repair; unsupported chemistry stays blocked.

        This is only used for single ordinary rates. The producer, not these
        eligibility checks, establishes the condition-bound omission claim.
        Original rate results, their library, and the unresolved conformer
        basis are never edited or retrospectively relabelled.
        """
        basis = replacement.get("conformer_rate_basis")
        if (not isinstance(basis, dict) or basis.get("network_reinjection_allowed") is not False
                or "forward_only_approximation_evidence" in replacement):
            return replacement, library_paths
        handoff_dir = Path(iteration["artifacts"]["directory"]) / "model-handoff"
        diagnostic = {
            "kind": "bounded_forward_only_handoff_attempt",
            "route_id": replacement.get("route_id"),
            "condition_contract_sha256": self.condition_sha256,
            "status": "not_attempted", "blocking_reasons": [],
            "original_rate_result": iteration["artifacts"].get("rate_result"),
            "original_library_paths": [str(path) for path in library_paths],
            "new_quantum_calculations": 0,
            "conclusive_observable_allowed": False, "lifetime_claim_allowed": False,
        }
        forward = replacement.get("conditional_forward_rate_evidence")
        reverse = replacement.get("reverse_rate_evidence")
        if not isinstance(forward, dict) or not isinstance(reverse, dict):
            diagnostic["blocking_reasons"].append("retained_structured_forward_reverse_evidence_missing")
        reactants, products = replacement.get("reactants"), replacement.get("products")
        labels = lambda values: [re.sub(r"\(\d+\)$", "", str(value).strip()) for value in values or ()]
        direct_channels = [
            item for item in (current_evidence.get("mechanism_inspection") or {}).get("reactions") or []
            if ((labels(item.get("reactant_labels")) == labels(reactants)
                 and labels(item.get("product_labels")) == labels(products))
                or (labels(item.get("reactant_labels")) == labels(products)
                    and labels(item.get("product_labels")) == labels(reactants)))
            and (_finite(item.get("evaluated_equilibrium_constant_Kc")) or 0.0) > 0.0
        ]
        if not reactants or not products or len(direct_channels) != 1:
            diagnostic["blocking_reasons"].append("direct_native_equilibrium_constant_missing_or_ambiguous")
        if len(library_paths) != 1:
            diagnostic["blocking_reasons"].append("single_original_rate_library_required")
        if not diagnostic["blocking_reasons"]:
            try:
                from .bounded_handoff import write_forward_only_approximation_library

                model_file = _write_json(handoff_dir / "source-model.json", {
                    "condition_contract": self.condition_contract, "rmg_evidence": current_evidence,
                })
                rate_file = _write_json(handoff_dir / "retained-rate-evidence.json", {
                    "conditional_forward_rate_evidence": forward, "reverse_rate_evidence": reverse,
                })
                diagnostic["source_model"] = model_file
                diagnostic["retained_rate_evidence"] = rate_file
                derived = write_forward_only_approximation_library(
                    library_paths[0], model_file, rate_file, self.condition_contract,
                    handoff_dir / "derived-library",
                )
                derived_replacement, derived_paths = _normalize_replacement(
                    {"replacement": derived, "library": derived.get("library")}, route,
                )
                for key in ("method_applicability", "method_applicability_evidence",
                            "rate_uncertainty_propagation_required", "required_uncertainty_components"):
                    if key in replacement:
                        derived_replacement.setdefault(key, replacement[key])
                diagnostic.update(status="derived_conditional_library_ready_for_validation",
                                  derived_library_paths=[str(path) for path in derived_paths],
                                  reverse_relevance_assessment=derived_replacement.get("reverse_relevance_assessment"))
                iteration["artifacts"]["derived_rate_result"] = _write_json(
                    handoff_dir / "derived-rate-result.json", derived_replacement,
                )
                replacement, library_paths = derived_replacement, derived_paths
            except (ValueError, KeyError, TypeError, OSError, SyntaxError) as error:
                # Unsupported scope/insufficient bound is a scientific refusal,
                # not a failed Arkane calculation. The original authoritative
                # conformer gate remains in effect below.
                diagnostic.update(status="bounded_handoff_unresolved",
                                  blocking_reasons=[str(error)], error_type=type(error).__name__)
                assessment = getattr(error, "assessment", None)
                if isinstance(assessment, dict):
                    diagnostic["reverse_relevance_assessment"] = assessment
                    iteration["reverse_relevance_assessment"] = assessment
        iteration["artifacts"]["bounded_handoff_diagnostic"] = _write_json(
            handoff_dir / "bounded-handoff-diagnostic.json", diagnostic,
        )
        return replacement, library_paths

    def _summary(self, status: str, *, reason: str, evidence: dict | None = None,
                 plan: dict | None = None, robustness: dict | None = None,
                 envelope: dict | None = None,
                 conditional_prediction: dict | None = None) -> dict:
        propagation = (evidence or {}).get("independent_cantera_propagation") or {}
        flux = propagation.get("reaction_flux_attribution") or {}
        commitment = (plan or {}).get("target_loss_commitment") or assess_target_loss_commitment(
            flux, str(self.condition_contract.get("target_label") or "stability"),
        )
        t95 = _finite((robustness or {}).get("t95_seconds"))
        target_loss = _model_loss(evidence or {})
        verified_coverage = _finite((plan or {}).get("verified_flux_fraction"))
        if (plan or {}).get("status") == "no_target_destruction_flux" and status == "verified_below_loss_threshold":
            verified_coverage = 1.0
        unresolved_upper = _unresolved_loss_upper_bound(evidence or {}, envelope)
        collision_bounded = (
            status == "verified_below_loss_threshold"
            and (envelope or {}).get("kind") == "collision_limit_rate_envelope"
            and bool((envelope or {}).get("observable_classification_robust"))
        )
        model_only_no_loss = status == "no_target_loss_in_completed_rmg_model"
        structural_only_handoffs = [
            {
                "route_id": item["route_id"],
                "path_verification_status": item.get("path_verification_status"),
                "path_result": item["artifacts"].get("path_result"),
                "rmg_solver_blocking_reasons": item.get("rmg_solver_blocking_reasons", []),
                "rate_claim_allowed": False, "observable_claim_allowed": False,
            }
            for item in self.iterations
            if item.get("status") == "structural_evidence_retained_without_kinetic_model"
        ]
        structural_only = bool(structural_only_handoffs) and not self.replacements
        verified_path_route_ids = [
            item["route_id"] for item in self.iterations
            if item.get("path_verification_status") in _PATH_VERIFIED_STATUSES
            or item.get("path_verification_status") == "verified_specialist_rate_model"
        ]
        blocked_rates = list(self.conditional_rate_records)
        conditional_model = status == "conditional_model_prediction"
        conditional_replacements = [
            item for item in self.replacements
            if "forward_only_approximation_evidence" in item
        ]
        observable_replacements = [
            item for item in self.replacements
            if "forward_only_approximation_evidence" not in item
        ]
        if conditional_model:
            # The repaired fixed model is useful evidence, but its central
            # rates and omission bound do not certify a physical observable.
            target_loss = unresolved_upper = None
            verified_coverage = 0.0
            commitment = {"status": "not_certified_for_conditional_prediction",
                          "verified_flux_credit": False}
        elif blocked_rates and conditional_replacements:
            # The full repaired model remains retained in iteration artifacts,
            # but its refused conditional rates earn no physical loss credit.
            target_loss = unresolved_upper = None
            verified_coverage = 0.0
        elif structural_only:
            target_loss = verified_coverage = unresolved_upper = None
        elif model_only_no_loss:
            verified_coverage = None
        if conditional_model:
            route_status = "verified"
            rate_status = "conditional_model_rates_inserted"
            propagation_status = "completed_conditional_model_not_verified_observable"
        elif structural_only:
            route_status = "independent_structural_evidence_only"
            rate_status = "unavailable_without_completed_kinetic_model"
            propagation_status = "unavailable_for_independent_structural_handoff"
        elif model_only_no_loss:
            route_status = rate_status = "not_required_for_retained_model"
            propagation_status = "completed_model"
        elif collision_bounded:
            route_status = "not_required_for_collision_bounded_observable"
            rate_status = "collision_bounded_for_observable_not_rate_verified"
            propagation_status = "verified_converged_under_collision_rate_envelope"
        elif blocked_rates:
            route_status = "verified" if verified_path_route_ids else "unresolved"
            rate_status = (
                "conditional_model_rate_calculated_network_admission_blocked"
                if all(item["status"] == "conditional_model_rate_calculated"
                       for item in blocked_rates) else
                "calculated_rate_not_validated_network_admission_blocked"
            )
            propagation_status = (
                "repropagated_conditional_model_scope_rejected"
                if any(item["network_admission"]["repropagated"] for item in blocked_rates)
                else "not_repropagated_after_blocked_rate_admission"
            )
        elif self.deferred_routes:
            route_status = "deferred_unverified"
            rate_status = "unavailable_for_deferred_routes"
            propagation_status = "screening_model_only"
        else:
            route_status = "verified" if self.replacements or verified_path_route_ids else "not_required"
            rate_status = ("conditional_model_rate_calculated" if conditional_replacements else
                           "validated" if self.replacements else
                           "unavailable_after_verified_path" if verified_path_route_ids else "not_required")
            propagation_status = ("verified_converged" if status in {
                "verified_t95_converged", "verified_below_loss_threshold",
            } else "incomplete")
        if conditional_model:
            claim_scope = "conditional_fixed_kinetic_model_prediction_no_physical_lifetime"
        elif model_only_no_loss:
            claim_scope = "completed_retained_rmg_model_only"
        elif collision_bounded:
            claim_scope = "condition_specific_retained_network_collision_bounded_observable"
        elif status in {"verified_t95_converged", "verified_below_loss_threshold"}:
            claim_scope = "verified_condition_specific_kinetic_model"
        else:
            claim_scope = "verification_incomplete_no_condition_lifetime"
        summary = {
            "schema_version": 1,
            "kind": "verification_summary",
            "status": status,
            "reason": reason,
            "condition_contract": self.condition_contract,
            "condition_contract_sha256": self.condition_sha256,
            "condition_contract_match": True,
            "estimated_time_to_retention_seconds": t95 if status == "verified_t95_converged" else None,
            "t95_seconds": t95 if status == "verified_t95_converged" else None,
            "target_loss_fraction": target_loss,
            "target_loss_commitment": commitment,
            "committed_target_loss_kmol": _finite(commitment.get("committed_target_loss_kmol")),
            "parent_reformation_credit_kmol": _finite(commitment.get("integrated_parent_reformation_kmol")),
            "verified_flux_coverage": verified_coverage,
            "unresolved_loss_upper_bound": unresolved_upper,
            "route_verification_status": route_status,
            "rate_verification_status": rate_status,
            "propagation_status": propagation_status,
            "verified_route_ids": [item.get("route_id") for item in observable_replacements],
            "conditional_model_route_ids": [item.get("route_id") for item in conditional_replacements],
            "verified_path_route_ids": verified_path_route_ids,
            "conditional_rate_records": blocked_rates,
            "conditional_rate_route_ids": [item["route_id"] for item in blocked_rates],
            "network_admission": {
                "status": ("blocked" if blocked_rates else
                           "admitted_for_conditional_model_only" if conditional_model else
                           "provisional_conditional_model_repair_unresolved" if conditional_replacements else
                           "admitted" if self.replacements else "not_attempted"),
                "blocking_reasons": sorted({
                    reason for item in blocked_rates
                    for reason in item["network_admission"]["blocking_reasons"]
                }),
                "blocked_route_ids": [item["route_id"] for item in blocked_rates],
                "conditional_rates_credited_as_verified_flux": False,
            },
            "stationary_point_evidence": [
                {"route_id": item["route_id"],
                 "status": item["stationary_point_verification_status"],
                 "path_result": item["artifacts"].get("path_result")}
                for item in self.iterations
                if item.get("stationary_point_verification_status")
            ],
            "native_candidate_handoffs": [
                {
                    "route_id": item.get("route_id"),
                    "native_candidate_id": item.get("native_candidate_id"),
                    "discovery_sources": list(item.get("discovery_sources") or ()),
                    "status": ("conditional_model_rate_inserted_and_repropagated"
                               if "forward_only_approximation_evidence" in item else
                               "verified_rate_inserted_and_repropagated"),
                }
                for item in self.replacements
                if item.get("native_candidate_id")
            ],
            "structural_only_handoffs": structural_only_handoffs,
            "deferred_unverified_routes": list(self.deferred_routes.values()),
            "iterations_completed": len(self.iterations),
            "robustness": robustness,
            "collision_limit_rate_envelope": envelope if collision_bounded else None,
            "conditional_model_prediction": conditional_prediction,
            "next_evidence_requests": ((conditional_prediction or {}).get("next_evidence_requests") or [
                request for item in self.iterations
                for request in item.get("expanded_network_next_evidence_requests", [])
            ]),
            "retained_expanded_network_assessments": [
                item["artifacts"]["expanded_network_assessment"] for item in self.iterations
                if item.get("artifacts", {}).get("expanded_network_assessment")
            ],
            "allocation_policy": (plan or {}).get("allocation_policy"),
            "artifacts": {
                "session_directory": str(self.session),
                "engine_state": str(self.session / "engine-state.json"),
                "verification_summary": str(self.session / "verification-summary.json"),
            },
            "screening_estimate_prohibited": status in {
                "verification_incomplete", "no_target_loss_in_completed_rmg_model",
                "conditional_model_prediction",
            },
            "claim_scope": claim_scope,
            "model_coverage": None if structural_only else propagation.get("coverage"),
        }
        _write_json(self.session / "verification-summary.json", summary)
        self._persist_state("completed" if status != "verification_incomplete" else "failed_closed",
                            verification_summary=summary)
        return summary

    def _incomplete(self, reason: str, *, evidence: dict | None = None,
                    plan: dict | None = None, error: Exception | None = None) -> dict:
        detail = reason if error is None else f"{reason}: {type(error).__name__}: {error}"
        return self._summary("verification_incomplete", reason=detail, evidence=evidence, plan=plan)

    def _terminal_if_robust(
        self,
        evidence: dict,
        plan: dict,
        envelope: dict | None,
        *,
        rate_uncertainty_required_route_ids: Iterable[str] = (),
    ) -> dict | None:
        if any("forward_only_approximation_evidence" in item for item in self.replacements):
            return None
        kinetics = evidence.get("kinetics_validation") or {}
        kinetics_passed = kinetics.get("status") == "passed" and not kinetics.get("violation_count", 0)
        if not kinetics_passed:
            return None
        if plan.get("status") == "no_target_destruction_flux" and not self.replacements:
            bound_status = str((envelope or {}).get("verification_status") or "")
            if not bound_status.startswith("verified"):
                return self._summary(
                    "no_target_loss_in_completed_rmg_model",
                    reason="No target-loss flux occurred in the completed retained RMG model; this is not universal stability evidence.",
                    evidence=evidence, plan=plan,
                )
        robustness = assess_t95_robustness(
            evidence, plan,
            retention_fraction=float(self.condition_contract["retention_fraction"]),
            unresolved_sensitivity_envelope=envelope,
            rate_uncertainty_required_route_ids=rate_uncertainty_required_route_ids,
            relative_tolerance=self.config.t95_relative_tolerance,
        )
        if robustness.get("status") == "robust_verified_t95":
            return self._summary("verified_t95_converged", reason=robustness.get("reason", "robust_t95"),
                                 evidence=evidence, plan=plan, robustness=robustness, envelope=envelope)
        if robustness.get("status") == "robust_no_t95_within_horizon":
            loss = _model_loss(evidence)
            unresolved_upper = _unresolved_loss_upper_bound(evidence, envelope)
            threshold = 1.0 - float(self.condition_contract["retention_fraction"])
            if loss is not None and loss + unresolved_upper < threshold:
                return self._summary(
                    "verified_below_loss_threshold", reason=robustness.get("reason", "robust_below_threshold"),
                    evidence=evidence, plan=plan, robustness=robustness, envelope=envelope,
                )
        return None

    def run(self) -> dict:
        self._persist_state("started")
        try:
            return self._run_started()
        except BaseException as error:
            # Normal fail-closed returns persist their own terminal summary.
            # This guard covers an exception or user interrupt which escapes
            # the chemistry callbacks, preventing a dead process from leaving
            # an apparently live ``running`` engine state behind.
            self._persist_state(
                "failed_closed",
                terminal_error={
                    "error_type": type(error).__name__,
                    "message": str(error),
                },
            )
            raise

    def _run_started(self) -> dict:
        try:
            current_evidence, _ = _extract_model_result(self.initial_result)
        except Exception as error:
            return self._incomplete("initial_model_evidence_invalid", error=error)
        current_result = self.initial_result

        # The extra pass performs the terminal robustness check after the last
        # permitted chemistry iteration; it may not start another verification.
        for iteration_number in range(1, self.config.max_iterations + 2):
            propagation = current_evidence.get("independent_cantera_propagation") or {}
            coverage = propagation.get("coverage") or {}
            screening_t95 = _finite(propagation.get("estimated_time_to_retention_seconds"))
            screening_loss = _finite(propagation.get("target_loss_fraction"))
            kinetics = current_evidence.get("kinetics_validation") or {}
            threshold = 1.0 - float(self.condition_contract["retention_fraction"])
            if (
                not self.replacements
                and not self.supplemental_candidate_routes
                and not _solver_gates(
                    current_evidence, max_closure_error=self.config.max_flux_closure_error,
                    condition_contract=self.condition_contract,
                )
                and propagation.get("status") == "completed"
                and coverage.get("complete")
                and screening_t95 is None
                and screening_loss is not None
                and screening_loss < threshold
                and kinetics.get("status") == "passed"
                and not kinetics.get("violation_count", 0)
            ):
                return self._summary(
                    "no_target_loss_in_completed_rmg_model",
                    reason=("Completed retained RMG propagation stayed strictly below the immutable loss threshold; "
                            "no ORCA route verification was started and this is not universal stability evidence."),
                    evidence=current_evidence,
                )
            planning_replacements = _planning_replacements(
                self.replacements, current_evidence, self.condition_contract,
            )
            rate_uncertainty_required_route_ids = {
                str(item.get("route_id")) for item in planning_replacements
                if item.get("route_id") and item.get("rate_uncertainty_propagation_required")
            }
            required_uncertainty_components = {
                str(item["route_id"]): tuple(item.get("required_uncertainty_components") or ())
                for item in planning_replacements
                if item.get("route_id") and item.get("required_uncertainty_components")
            }
            rmg_plan = build_flux_ranked_verification_plan(
                current_evidence, self.condition_contract,
                verified_replacements=planning_replacements,
                deferred_route_ids=self.deferred_routes,
                rate_uncertainties=self.rate_uncertainties,
                flux_coverage_target=self.config.flux_coverage_target,
                max_closure_error=self.config.max_flux_closure_error,
            )
            allocation = rmg_plan.get("allocation_policy") or {}
            missing_sensitivity = allocation.get("missing_t95_sensitivity_route_ids") or []
            if (
                self.rate_uncertainties
                and not (allocation.get("missing_rate_uncertainty_route_ids") or [])
                and missing_sensitivity
                and rmg_plan.get("status") not in {
                    "insufficient_solver_evidence", "insufficient_route_attribution",
                    "controlling_route_has_unresolved_initiation", "partial_enumeration_no_direct_route",
                }
            ):
                try:
                    sensitivity = run_t95_sensitivity_attribution(
                        current_evidence,
                        rmg_plan,
                        self.condition_contract,
                        route_ids=missing_sensitivity,
                        rmg_env=self.rmg_env,
                    )
                except Exception as error:
                    sensitivity = {
                        "schema_version": 1,
                        "kind": "reaction_t95_sensitivity_attribution",
                        "status": "incomplete",
                        "reason": f"{type(error).__name__}: {error}",
                        "requested_route_ids": list(missing_sensitivity),
                        "reactions": [],
                        "interpretation": (
                            "The optional allocation diagnostic failed; route ordering fell back "
                            "without changing any scientific convergence gate."
                        ),
                    }
                propagation["reaction_sensitivity_attribution"] = sensitivity
                _write_json(
                    self.session / f"t95-sensitivity-{iteration_number:03d}.json",
                    sensitivity,
                )
                rmg_plan = build_flux_ranked_verification_plan(
                    current_evidence, self.condition_contract,
                    verified_replacements=planning_replacements,
                    deferred_route_ids=self.deferred_routes,
                    rate_uncertainties=self.rate_uncertainties,
                    flux_coverage_target=self.config.flux_coverage_target,
                    max_closure_error=self.config.max_flux_closure_error,
                )
            _write_json(self.session / f"plan-{iteration_number:03d}.json", rmg_plan)
            partial_fallback = (
                rmg_plan.get("partial_enumeration") is True
                and rmg_plan.get("status") == "partial_enumeration_no_direct_route"
            )
            independent_structural_fallback = (
                rmg_plan.get("status") == "insufficient_solver_evidence"
                and not self.replacements
            )
            supplemental_plan = self._supplemental_plan(
                require_direct_reactants=partial_fallback or independent_structural_fallback,
                require_independent_gas=independent_structural_fallback,
            )
            use_partial_fallback = (
                partial_fallback and supplemental_plan.get("status") == "verification_required"
            )
            use_structural_fallback = (
                independent_structural_fallback
                and supplemental_plan.get("status") == "verification_required"
            )
            if rmg_plan.get("status") in {
                "insufficient_solver_evidence", "insufficient_route_attribution",
                "controlling_route_has_unresolved_initiation", "partial_enumeration_no_direct_route",
            } and not (use_partial_fallback or use_structural_fallback):
                return self._incomplete(f"flux_plan_{rmg_plan.get('status')}", evidence=current_evidence, plan=rmg_plan)

            use_supplemental = (
                supplemental_plan.get("status") == "verification_required"
                and (use_partial_fallback or use_structural_fallback
                     or rmg_plan.get("status") in {"no_target_destruction_flux", "flux_coverage_target_met"})
            )
            if use_partial_fallback:
                supplemental_plan["partial_enumeration"] = True
                supplemental_plan["rmg_planning_status"] = rmg_plan["status"]
                supplemental_plan["observable_claim_allowed"] = False
            if use_structural_fallback:
                supplemental_plan.update(
                    structural_only=True, rate_claim_allowed=False, observable_claim_allowed=False,
                    rmg_planning_status=rmg_plan["status"],
                    rmg_solver_blocking_reasons=list(rmg_plan.get("blocking_reasons") or ()),
                )
            plan = supplemental_plan if use_supplemental else rmg_plan
            if use_supplemental:
                _write_json(self.session / f"supplemental-plan-{iteration_number:03d}.json", plan)

            envelope = _sensitivity_envelope(current_result, current_evidence)
            if rate_uncertainty_required_route_ids and envelope is not None \
                    and envelope.get("uncertain_route_ids") is None:
                # An unresolved-only envelope cannot establish robustness for
                # uncertainty in newly calculated central rates.
                envelope = None
            if self.replacements and envelope is None:
                prediction_rate_envelope = bool(rate_uncertainty_required_route_ids)
                bound_preparation = (
                    prepare_combined_prediction_rate_bounds(
                        rmg_plan,
                        rate_uncertainty_required_route_ids=rate_uncertainty_required_route_ids,
                        required_uncertainty_components_by_route=required_uncertainty_components,
                    )
                    if prediction_rate_envelope else
                    prepare_combined_unresolved_rate_bounds(current_evidence, rmg_plan)
                )
                _write_json(
                    self.session / (
                        f"combined-prediction-rate-bound-preparation-{iteration_number:03d}.json"
                        if prediction_rate_envelope else
                        f"combined-rate-bound-preparation-{iteration_number:03d}.json"
                    ),
                    bound_preparation,
                )
                if bound_preparation.get("status") == "completed":
                    try:
                        envelope = (
                            run_combined_prediction_rate_sensitivity(
                                current_evidence,
                                self.condition_contract,
                                prepared_bounds=bound_preparation,
                                rmg_env=self.rmg_env,
                            )
                            if prediction_rate_envelope else
                            run_combined_unresolved_sensitivity(
                                current_evidence,
                                rmg_plan,
                                self.condition_contract,
                                rate_bounds=bound_preparation["rate_bounds"],
                                rmg_env=self.rmg_env,
                            )
                        )
                    except Exception as error:
                        envelope = {
                            "schema_version": 1,
                            "kind": (
                                "combined_prediction_rate_sensitivity" if prediction_rate_envelope
                                else "combined_unresolved_route_sensitivity"
                            ),
                            "status": "incomplete",
                            "reason": f"{type(error).__name__}: {error}",
                            "combined_perturbation": True,
                            (
                                "uncertain_route_ids" if prediction_rate_envelope
                                else "unresolved_route_ids"
                            ): bound_preparation[
                                "uncertain_route_ids" if prediction_rate_envelope
                                else "unresolved_route_ids"
                            ],
                        }
                    _write_json(
                        self.session / (
                            f"combined-prediction-rate-sensitivity-{iteration_number:03d}.json"
                            if prediction_rate_envelope else
                            f"combined-unresolved-sensitivity-{iteration_number:03d}.json"
                        ),
                        envelope,
                    )
                    if isinstance(current_result, dict):
                        current_result[
                            "combined_prediction_rate_sensitivity" if prediction_rate_envelope
                            else "combined_unresolved_route_sensitivity"
                        ] = envelope
            # A direct unimolecular RMG route is deliberately queued after an
            # invalid kinetic screen even when its estimated flux is zero.  It
            # is not eligible for the usual no-loss robustness shortcut until
            # ORCA has accepted or deferred that actual RMG hypothesis.
            deferred_supplemental = self._deferred_supplemental_route_ids()
            collision_envelope = None
            collision_envelope_eligible = (
                str(kinetics.get("status") or "") == "kinetics_unreliable"
                and int(kinetics.get("violation_count") or 0) > 0
                and not _solver_gates(
                    current_evidence, max_closure_error=self.config.max_flux_closure_error,
                    condition_contract=self.condition_contract,
                )
                and propagation.get("status") == "completed"
                and coverage.get("complete")
                and screening_t95 is None
                and screening_loss is not None
                and screening_loss < threshold
                and rmg_plan.get("status") in {
                    "no_target_destruction_flux", "flux_coverage_target_met",
                }
                and not rmg_plan.get("intrinsic_failover_required")
                and supplemental_plan.get("status") == "no_graph_candidate_ready"
                and not deferred_supplemental
                and not self.deferred_routes
            )
            if collision_envelope_eligible:
                try:
                    collision_envelope = run_collision_limit_rate_envelope(
                        current_evidence,
                        self.condition_contract,
                        rmg_env=self.rmg_env,
                    )
                except Exception as error:
                    collision_envelope = {
                        "schema_version": 1,
                        "kind": "collision_limit_rate_envelope",
                        "status": "incomplete",
                        "reason": f"{type(error).__name__}: {error}",
                        "observable_classification_robust": False,
                    }
                _write_json(
                    self.session / f"collision-limit-rate-envelope-{iteration_number:03d}.json",
                    collision_envelope,
                )
                maximum_loss = _finite(
                    collision_envelope.get("maximum_loss_upper_bound_fraction")
                )
                if (
                    collision_envelope.get("status") == "completed"
                    and collision_envelope.get("observable_classification_robust") is True
                    and maximum_loss is not None
                    and maximum_loss < threshold
                ):
                    robustness = {
                        "status": "robust_no_t95_within_horizon",
                        "reason": "collision_bounded_retained_network_remains_below_loss_threshold",
                        "t95_seconds": None,
                        "loss_upper_bound_fraction": maximum_loss,
                        "loss_threshold_fraction": threshold,
                        "scope": "retained_network_collision_rate_envelope",
                    }
                    return self._summary(
                        "verified_below_loss_threshold",
                        reason=robustness["reason"],
                        evidence=current_evidence,
                        plan=rmg_plan,
                        robustness=robustness,
                        envelope=collision_envelope,
                    )
            terminal = None if (
                use_supplemental
                or deferred_supplemental
                or plan.get("intrinsic_failover_required")
            ) else self._terminal_if_robust(
                current_evidence, plan, envelope,
                rate_uncertainty_required_route_ids=rate_uncertainty_required_route_ids,
            )
            if terminal is not None:
                return terminal

            required_calculated_routes = [
                item for item in plan.get("ranked_routes") or ()
                if item.get("verified_replacement")
                and str(item.get("route_id")) in rate_uncertainty_required_route_ids
                and (_finite(item.get("gross_target_destruction_fraction")) or 0.0) > 0.0
            ]
            missing_calculated_uncertainty = sorted(
                str(item.get("route_id")) for item in required_calculated_routes
                if not item.get("rate_uncertainty_evidence")
            )
            incomplete_component_coverage = {
                str(item.get("route_id")): sorted(
                    set(required_uncertainty_components.get(str(item.get("route_id")), ()))
                    - set(
                        (item.get("rate_uncertainty_evidence") or {}).get(
                            "uncertainty_components"
                        ) or ()
                    )
                )
                for item in required_calculated_routes
                if "total_rate_model" not in set(
                    (item.get("rate_uncertainty_evidence") or {}).get(
                        "uncertainty_components"
                    ) or ()
                )
            }
            incomplete_component_coverage = {
                route_id: missing for route_id, missing in incomplete_component_coverage.items()
                if missing
            }
            if missing_calculated_uncertainty and plan.get("status") in {
                "flux_coverage_target_met", "no_target_destruction_flux",
            }:
                return self._incomplete(
                    "calculated_rate_uncertainty_missing_for_"
                    + ",".join(missing_calculated_uncertainty),
                    evidence=current_evidence, plan=plan,
                )
            if incomplete_component_coverage and plan.get("status") in {
                "flux_coverage_target_met", "no_target_destruction_flux",
            }:
                detail = ";".join(
                    f"{route_id}:{','.join(components)}"
                    for route_id, components in sorted(incomplete_component_coverage.items())
                )
                return self._incomplete(
                    "calculated_rate_uncertainty_components_missing_for_" + detail,
                    evidence=current_evidence, plan=plan,
                )
            if required_calculated_routes \
                    and (envelope or {}).get("status") != "completed" \
                    and plan.get("status") in {
                        "flux_coverage_target_met", "no_target_destruction_flux",
                    }:
                return self._incomplete(
                    "calculated_rate_uncertainty_propagation_incomplete",
                    evidence=current_evidence, plan=plan,
                )

            # A deferred graph candidate cannot be discarded just because the
            # retained RMG model has no loss flux.  It is independent,
            # unresolved mechanism evidence, so keep the scientific outcome
            # fail-closed rather than restoring the RMG-only no-loss claim.
            if (
                not use_supplemental
                and rmg_plan.get("status") in {"no_target_destruction_flux", "flux_coverage_target_met"}
                and deferred_supplemental
            ):
                return self._incomplete(
                    "supplemental_graph_candidates_deferred_without_verified_rate",
                    evidence=current_evidence,
                    plan=rmg_plan,
                )

            if iteration_number > self.config.max_iterations:
                return self._incomplete("maximum_verification_iterations_reached",
                                        evidence=current_evidence, plan=plan)

            # Reaching the configured coverage is not terminal by itself.  If
            # no physical sensitivity envelope proves robustness, continue to
            # the next nonzero-flux route until coverage is complete.
            selection_plan = plan
            if not use_supplemental and plan.get("status") in {"flux_coverage_target_met", "no_target_destruction_flux"}:
                selection_plan = build_flux_ranked_verification_plan(
                    current_evidence, self.condition_contract,
                    verified_replacements=planning_replacements,
                    deferred_route_ids=self.deferred_routes,
                    rate_uncertainties=self.rate_uncertainties,
                    flux_coverage_target=1.0,
                    max_closure_error=self.config.max_flux_closure_error,
                )
                _write_json(self.session / f"full-coverage-plan-{iteration_number:03d}.json", selection_plan)
            if selection_plan.get("status") != "verification_required":
                return self._incomplete(
                    f"robustness_not_established_after_{selection_plan.get('status')}",
                    evidence=current_evidence, plan=selection_plan,
                )

            try:
                multistep_stage: str | None = None
                route = _route_for_plan(current_evidence, selection_plan, self.condition_contract)
            except Exception as error:
                return self._incomplete("selected_route_invalid", evidence=current_evidence,
                                        plan=selection_plan, error=error)
            route_id = str(route.get("route_id"))
            folder_name = re.sub(r"[^A-Za-z0-9_.-]+", "-", route_id).strip("-") or "route"
            iteration_dir = self.session / f"iteration-{iteration_number:03d}-{folder_name}"
            iteration_dir.mkdir(parents=True, exist_ok=False)
            iteration = {
                "iteration": iteration_number,
                "route_id": route_id,
                "status": "started",
                "artifacts": {
                    "directory": str(iteration_dir),
                    "route": _write_json(iteration_dir / "route.json", route),
                    "selection_plan": _write_json(iteration_dir / "selection-plan.json", selection_plan),
                    "condition_contract": _write_json(iteration_dir / "condition-contract.json", self.condition_contract),
                },
            }
            self.iterations.append(iteration)
            self._persist_state("running")

            try:
                attempts = self.route_policy_attempts.setdefault(route_id, [])
                specialist_rate_result: dict | None = None
                path_result: dict | None = None
                same_policy_recovery = route_id in self.same_policy_execution_recovery_route_ids
                retained_path = self.retained_specialist_path_results.get(route_id)
                if self.path_verifier_policy_id in attempts \
                        and not same_policy_recovery and retained_path is None:
                    raise ValueError(
                        f"Route {route_id} was already attempted under path-verifier policy "
                        f"{self.path_verifier_policy_id}; an unchanged policy cannot retry it"
                    )
                if self.path_verifier_policy_id not in attempts:
                    attempts.append(self.path_verifier_policy_id)
                if same_policy_recovery:
                    # This is another execution segment of the same scientific
                    # attempt, not permission to spend a second policy retry.
                    self.same_policy_execution_recovery_route_ids.remove(route_id)
                    iteration["same_policy_execution_recovery"] = True
                calculation_dir = self.calculation_workspaces.get(route_id, iteration_dir)
                calculation_dir.mkdir(parents=True, exist_ok=True)
                iteration["path_verifier_policy_id"] = self.path_verifier_policy_id
                iteration["artifacts"]["calculation_workspace"] = str(calculation_dir)
                iteration["resumed_calculation_workspace"] = calculation_dir != iteration_dir
                self._persist_state("running")
                if retained_path is not None:
                    path_result = dict(retained_path["path_result"])
                    retained_status = _unverified_path_status(path_result)
                    if retained_status not in {
                        "barrierless_rate_model_required",
                        "dissociation_rate_model_required",
                    }:
                        raise ValueError(
                            "Retained path is not eligible for rate-only specialist continuation: "
                            + retained_status
                        )
                    iteration["retained_path_reused"] = True
                    iteration["artifacts"]["retained_path_reuse"] = _write_json(
                        iteration_dir / "retained-path-reuse.json", {
                            "kind": "retained_specialist_path_reuse",
                            "route_id": route_id,
                            "source_path_result": retained_path.get("source_path"),
                            "source_path_result_sha256": retained_path.get("source_sha256"),
                            "path_status": retained_status,
                            "path_verifier_reexecuted": False,
                        },
                    )
                else:
                    path_result = self.verify_route_callback(route, calculation_dir / "path")
                iteration["artifacts"]["path_result"] = _write_json(iteration_dir / "path-result.json", path_result)
                discovery_handoff = self._retain_unexpected_candidates(route, path_result)
                if discovery_handoff["records"]:
                    iteration["artifacts"]["unexpected_candidate_handoff"] = _write_json(
                        iteration_dir / "unexpected-candidate-handoff.json", discovery_handoff,
                    )
                    self._persist_state("running")
                path_status = _verified_status(path_result, _PATH_VERIFIED_STATUSES, nested_keys=(
                    "path_result", "route_verification", "verification", "orca_verification",
                ))
                if path_status is not None:
                    iteration["path_verification_status"] = path_status
                    path_block = _nested_dict(path_result, (
                        "path_result", "route_verification", "verification", "orca_verification",
                    )) or path_result
                    stationary = path_block.get("stationary_points")
                    if isinstance(stationary, dict):
                        iteration["stationary_point_verification_status"] = (
                            "validated" if stationary.get("valid") is True else "unresolved"
                        )
                if use_structural_fallback:
                    # Preserve the actual path outcome, including a verified
                    # saddle or a specialist-required refusal, but do not
                    # manufacture a rate/repair/observable when the quantitative
                    # RMG model is unavailable. The ordinary bounded queue
                    # continues to the next independent structural candidate.
                    retained_status = path_status or _unverified_path_status(path_result)
                    blockers = list(rmg_plan.get("blocking_reasons") or ())
                    iteration.update(
                        status="structural_evidence_retained_without_kinetic_model",
                        path_verification_status=retained_status,
                        rmg_solver_blocking_reasons=blockers,
                    )
                    self.deferred_routes[route_id] = {
                        "route_id": route_id,
                        "status": "independent_structural_path_retained_kinetic_model_unavailable",
                        "path_verification_status": retained_status,
                        "reason": "RMG quantitative model unavailable; no rate or observable admission",
                        "rmg_solver_blocking_reasons": blockers,
                        "path_result": iteration["artifacts"]["path_result"],
                        "rate_claim_allowed": False, "observable_claim_allowed": False,
                    }
                    self._persist_state("running_with_deferred_routes")
                    continue
                if path_status is None:
                    unverified_status = _unverified_path_status(path_result)
                    if (
                        unverified_status in {
                            "barrierless_rate_model_required",
                            "dissociation_rate_model_required",
                        }
                        and self.build_specialist_rate_callback is not None
                    ):
                        from .specialist_rates import (
                            build_specialist_rate_request,
                            specialist_rate_uncertainty_record,
                            validate_specialist_rate_result,
                        )

                        multistep_stage = "specialist_rate_generation"
                        request = build_specialist_rate_request(
                            route, path_result, self.condition_contract,
                        )
                        iteration["artifacts"]["specialist_rate_request"] = _write_json(
                            iteration_dir / "specialist-rate-request.json", request,
                        )
                        specialist_rate_result = self.build_specialist_rate_callback(
                            request, iteration_dir / "specialist-rate",
                        )
                        iteration["artifacts"]["specialist_rate_result"] = _write_json(
                            iteration_dir / "specialist-rate-result.json",
                            specialist_rate_result,
                        )
                        specialist_validation = validate_specialist_rate_result(
                            specialist_rate_result, request,
                        )
                        iteration["artifacts"]["specialist_rate_validation"] = _write_json(
                            iteration_dir / "specialist-rate-validation.json",
                            specialist_validation,
                        )
                        if specialist_validation.get("status") != "passed":
                            raise ValueError(
                                "Specialist rate failed validation: "
                                + ", ".join(
                                    specialist_validation.get("blocking_reasons") or ()
                                )
                            )
                        specialist_uncertainty = specialist_rate_uncertainty_record(
                            specialist_rate_result, request,
                        )
                        self.rate_uncertainties = tuple(
                            record for record in self.rate_uncertainties
                            if record.route_id != specialist_uncertainty.route_id
                        ) + (specialist_uncertainty,)
                        path_status = "verified_specialist_rate_model"
                    if unverified_status in _DEFERABLE_PATH_STATUSES:
                        if specialist_rate_result is None:
                            detail = _unverified_path_detail(path_result)
                            from .kinetic_uncertainty import (
                                uncertainty_from_forward_rate_upper_bound,
                            )
                            physical_uncertainty = uncertainty_from_forward_rate_upper_bound(
                                route, path_result, self.condition_contract,
                            )
                            if physical_uncertainty is not None:
                                self.rate_uncertainties = tuple(
                                    record for record in self.rate_uncertainties
                                    if record.route_id != physical_uncertainty.route_id
                                ) + (physical_uncertainty,)
                            iteration.update(
                                status="deferred_unverified",
                                path_verification_status=unverified_status,
                                failure_reason=detail,
                            )
                            self.deferred_routes[route_id] = {
                                "route_id": route_id,
                                "status": unverified_status,
                                "reason": detail,
                                "path_result": iteration["artifacts"]["path_result"],
                                "rate_uncertainty": (
                                    physical_uncertainty.as_dict()
                                    if physical_uncertainty is not None else None
                                ),
                            }
                            self._persist_state("running_with_deferred_routes")
                            continue
                    elif specialist_rate_result is None:
                        raise ValueError(
                            "Route callback did not return an explicitly verified path status: "
                            + _unverified_path_detail(path_result)
                        )

                if path_status == "verified_multistep_path":
                    multistep_stage = "segment_rate_generation"
                rate_result = (
                    specialist_rate_result
                    if specialist_rate_result is not None else
                    self.build_rate_callback(route, path_result, calculation_dir / "rate")
                )
                iteration["artifacts"]["rate_result"] = _write_json(iteration_dir / "rate-result.json", rate_result)
                if rate_result.get("kind") == "storca_verified_multistep_rate_bundle":
                    if rate_result.get("status") != "completed_multistep":
                        raise ValueError("Multistep rate bundle is not complete")
                    segment_rates = list(rate_result.get("segments") or [])
                    if len(segment_rates) < 2:
                        raise ValueError("Multistep rate bundle needs at least two elementary rates")
                    normalized: list[tuple[dict, list[Path], dict]] = []
                    seen_segment_ids: set[str] = set()
                    for segment in segment_rates:
                        segment_route = segment.get("route") or {}
                        segment_id = str(segment_route.get("route_id") or "")
                        if not segment_id or segment_id in seen_segment_ids:
                            raise ValueError("Multistep rate bundle has missing or duplicate segment identities")
                        seen_segment_ids.add(segment_id)
                        replacement, new_paths = _normalize_replacement(
                            segment.get("rate_result") or {}, segment_route,
                        )
                        validation = validate_rate_replacement_evidence(
                            replacement,
                            temperature_K=float(self.condition_contract["temperature_K"]),
                            pressure_bar=float(self.condition_contract["pressure_bar"]),
                            condition_contract=self.condition_contract,
                        )
                        normalized.append((replacement, new_paths, validation))
                    validations = {
                        "schema_version": 1,
                        "kind": "multistep_replacement_validation",
                        "status": (
                            "passed" if all(item[2].get("status") == "passed" for item in normalized)
                            else "failed"
                        ),
                        "parent_route_id": route_id,
                        "segments": [item[2] for item in normalized],
                    }
                    iteration["artifacts"]["replacement_validation"] = _write_json(
                        iteration_dir / "replacement-validation.json", validations,
                    )
                    if validations["status"] != "passed":
                        reasons = [
                            reason
                            for _, _, validation in normalized
                            for reason in validation.get("blocking_reasons") or []
                        ]
                        for replacement, new_paths, validation in normalized:
                            if validation.get("status") != "passed":
                                self._retain_blocked_rate_handoff(
                                    route_id=route_id, path_status=path_status,
                                    rate_result=rate_result, replacement=replacement,
                                    validation=validation, library_paths=new_paths,
                                    iteration=iteration,
                                )
                        if reasons and set(reasons) == {"conformer_thermo_basis_unresolved"}:
                            iteration.update(
                                status="network_admission_blocked",
                                scientific_failure_state="network_admission_blocked",
                                failure_reason=", ".join(sorted(set(reasons))),
                            )
                            return self._incomplete(
                                "network_admission_blocked: " + iteration["failure_reason"],
                                evidence=current_evidence, plan=selection_plan,
                            )
                        raise ValueError(
                            "Multistep replacement rates failed validation: " + ", ".join(reasons)
                        )
                    for replacement, new_paths, _ in normalized:
                        replacement["parent_multistep_route_id"] = route_id
                        self.replacements.append(replacement)
                        for path in new_paths:
                            if path not in self.library_paths:
                                self.library_paths.append(path)
                    replacement_route_ids = [item[0].get("route_id") for item in normalized]
                    replacement = normalized[-1][0]
                    multistep_stage = "network_insertion_and_repropagation"
                else:
                    replacement, new_paths = _normalize_replacement(rate_result, route)
                    replacement, new_paths = self._try_bounded_forward_only_handoff(
                        replacement=replacement, library_paths=new_paths, route=route,
                        current_evidence=current_evidence, iteration=iteration,
                    )
                    validation = validate_rate_replacement_evidence(
                        replacement,
                        temperature_K=float(self.condition_contract["temperature_K"]),
                        pressure_bar=float(self.condition_contract["pressure_bar"]),
                        condition_contract=self.condition_contract,
                    )
                    iteration["artifacts"]["replacement_validation"] = _write_json(
                        iteration_dir / "replacement-validation.json", validation,
                    )
                    if validation.get("status") != "passed":
                        self._retain_blocked_rate_handoff(
                            route_id=route_id, path_status=path_status,
                            rate_result=rate_result, replacement=replacement,
                            validation=validation, library_paths=new_paths,
                            iteration=iteration,
                        )
                        reasons = list(validation.get("blocking_reasons") or ())
                        if reasons and set(reasons) == {"conformer_thermo_basis_unresolved"}:
                            iteration.update(
                                status="network_admission_blocked",
                                scientific_failure_state="network_admission_blocked",
                                failure_reason=", ".join(reasons),
                            )
                            return self._incomplete(
                                "network_admission_blocked: " + iteration["failure_reason"],
                                evidence=current_evidence, plan=selection_plan,
                            )
                        raise ValueError("Replacement rate failed validation: " + ", ".join(validation.get("blocking_reasons") or []))
                    self.replacements.append(replacement)
                    for path in new_paths:
                        if path not in self.library_paths:
                            self.library_paths.append(path)
                    replacement_route_ids = [replacement.get("route_id")]

                if route.get("native_candidate_id"):
                    handoff = {
                        "schema_version": 1,
                        "kind": "storca_native_candidate_kinetic_handoff",
                        "status": (
                            "conditional_model_rate_ready_for_scoped_network_insertion"
                            if any("forward_only_approximation_evidence" in item
                                   for item in self.replacements
                                   if item.get("route_id") in replacement_route_ids)
                            else "verified_rate_ready_for_network_insertion"
                        ),
                        "route_id": route_id,
                        "native_candidate_id": route.get("native_candidate_id"),
                        "discovery_sources": list(route.get("discovery_sources") or ()),
                        "path_verification_status": path_status,
                        "replacement_route_ids": replacement_route_ids,
                        "replacement_library_paths": [str(path) for path in self.library_paths],
                        "condition_contract_sha256": self.condition_sha256,
                        "scientific_scope": (
                            "The discovery hypothesis has a verified path and calculated rate model, but it "
                            "cannot affect an observable until the repaired network accepts the "
                            "reaction and is repropagated under the identical condition contract."
                        ),
                    }
                    iteration["artifacts"]["native_candidate_kinetic_handoff"] = _write_json(
                        iteration_dir / "native-candidate-kinetic-handoff.json", handoff,
                    )

                rerun_result = self.rerun_model_callback(list(self.library_paths), calculation_dir / "rerun")
                iteration["artifacts"]["rerun_result"] = _write_json(iteration_dir / "rerun-result.json", rerun_result)
                repaired_evidence, retained_contract = _extract_model_result(rerun_result)
                if not isinstance(retained_contract, dict) or _condition_hash(retained_contract) != self.condition_sha256:
                    raise ValueError("Full model rerun did not retain the exact immutable condition contract")
                repair = assess_repaired_mechanism(repaired_evidence, self.replacements, self.condition_contract)
                iteration["artifacts"]["repaired_model_validation"] = _write_json(
                    iteration_dir / "repaired-model-validation.json", repair,
                )
                if repair.get("status") not in {
                    "accepted_for_flux_reranking",
                    "accepted_for_flux_reranking_with_unrepaired_collision_violations",
                }:
                    raise ValueError("Repaired model failed validation: " + ", ".join(repair.get("blocking_reasons") or []))
                conditional_inserted = [
                    item for item in self.replacements
                    if "forward_only_approximation_evidence" in item
                ]
                conditional_certificates = []
                expanded_network_reports = []
                if conditional_inserted:
                    from .bounded_handoff import validate_forward_only_repaired_model

                    for item in conditional_inserted:
                        certificate = validate_forward_only_repaired_model(
                            repaired_evidence, item, self.condition_contract,
                        )
                        # Never widen a historical certificate. A current,
                        # structurally admitted expanded linear model receives
                        # a separate source-bound assessment and diagnostics.
                        if (len(conditional_inserted) == 1
                                and certificate.get("status") != "passed"
                                and certificate.get("blocking_reasons") == [
                                    "bounded_handoff_repaired_inventory_changed_new_species_or_channels"]):
                            try:
                                from .network_influence import assess_retained_expanded_network

                                expanded = assess_retained_expanded_network(
                                    iteration["artifacts"]["rerun_result"], item,
                                    self.condition_contract, iteration_dir / "expanded-network",
                                    rmg_env=self.rmg_env,
                                )
                                expanded_network_reports.append(expanded)
                                iteration["expanded_network_next_evidence_requests"] = expanded["next_evidence_requests"]
                                iteration["artifacts"]["expanded_network_assessment"] = str(
                                    iteration_dir / "expanded-network" / "assessment-report.json")
                                old_certificate = certificate
                                certificate = {
                                    "kind": "expanded_conditional_model_inventory_validation",
                                    "status": "passed" if expanded["status"] == "passed" else "failed",
                                    "blocking_reasons": ([] if expanded["status"] == "passed" else
                                                         ["expanded_model_assessment_or_numerics_unresolved"]),
                                    "prior_scope_validation": old_certificate,
                                    "expanded_model_assessment": expanded["expanded_model_assessment"],
                                    "maximum_abs_parent_retention_error": expanded[
                                        "expanded_model_assessment"].get("maximum_abs_parent_retention_error"),
                                    "conclusive_observable_allowed": False,
                                }
                            except (ValueError, KeyError, TypeError, OSError, RuntimeError) as error:
                                certificate = {**certificate,
                                               "expanded_assessment_failure": f"{type(error).__name__}: {error}"}
                        conditional_certificates.append(certificate)
                        if certificate.get("status") != "passed":
                            self._retain_blocked_rate_handoff(
                                route_id=route_id, path_status=path_status,
                                rate_result=rate_result, replacement=item,
                                validation=validate_rate_replacement_evidence(
                                    item, temperature_K=float(self.condition_contract["temperature_K"]),
                                    pressure_bar=float(self.condition_contract["pressure_bar"]),
                                    condition_contract=self.condition_contract,
                                ),
                                library_paths=[Path(item["library"])] if item.get("library") else [],
                                iteration=iteration, network_validation=certificate,
                                repropagated=True,
                            )
                    iteration["artifacts"]["conditional_model_inventory_validation"] = _write_json(
                        iteration_dir / "conditional-model-inventory-validation.json",
                        {"kind": "conditional_repaired_model_inventory_validation",
                         "status": ("passed" if all(item.get("status") == "passed"
                                                     for item in conditional_certificates) else "failed"),
                         "certificates": conditional_certificates},
                    )
                    if any(item.get("status") != "passed" for item in conditional_certificates):
                        reasons = sorted({reason for item in conditional_certificates
                                          for reason in item.get("blocking_reasons") or ()})
                        iteration.update(
                            status="network_admission_blocked",
                            scientific_failure_state="conditional_repaired_model_scope_unresolved",
                            failure_reason=", ".join(reasons),
                        )
                        return self._incomplete(
                            "conditional_repaired_model_scope_unresolved: " + iteration["failure_reason"],
                            evidence=repaired_evidence, plan=selection_plan,
                        )
                    iteration["network_admission"] = {
                        "status": "admitted_for_conditional_model_only",
                        "blocking_reasons": [], "repropagated": True,
                        "verified_flux_credit": False, "observable_claim_allowed": False,
                        "lifetime_claim_allowed": False,
                    }
                if route.get("native_candidate_id"):
                    handoff.update(
                        status=("conditional_model_rate_inserted_and_repropagated"
                                if conditional_inserted else "verified_rate_inserted_and_repropagated"),
                        repaired_model_validation_status=repair.get("status"),
                        repaired_model_artifact=iteration["artifacts"]["repaired_model_validation"],
                    )
                    iteration["artifacts"]["native_candidate_kinetic_handoff"] = _write_json(
                        iteration_dir / "native-candidate-kinetic-handoff.json", handoff,
                    )
                    iteration["native_candidate_handoff_status"] = handoff["status"]
                iteration.update(status="completed", path_verification_status=path_status,
                                 replacement_route_id=replacement.get("route_id"),
                                 replacement_route_ids=replacement_route_ids,
                                 parent_multistep_route_resolved=(
                                     route_id if len(replacement_route_ids) > 1 else None
                                 ))
                current_evidence = repaired_evidence
                current_result = rerun_result
                self._persist_state("running")
                if conditional_inserted:
                    propagation = repaired_evidence.get("independent_cantera_propagation") or {}
                    expanded_baseline = (
                        next((row for row in expanded_network_reports[0]["scenario_response"]["scenarios"]
                              if row["scenario_id"] == "baseline"), None)
                        if len(expanded_network_reports) == 1 else None
                    )
                    conditional_prediction = {
                        "schema_version": 1,
                        "kind": "conditional_repaired_kinetic_model_prediction",
                        "model_scope": "declared_fixed_kinetic_model_only",
                        "approximation_error_scope": "omitted_reverse_channels_only",
                        "reverse_omission_bound_not_total_prediction_error": True,
                        "total_numerical_error_bound_established": False,
                        "condition_contract_sha256": self.condition_sha256,
                        "target_duration_seconds": self.condition_contract["target_duration_seconds"],
                        "target_profile": (expanded_baseline["target_profile"] if expanded_baseline
                                           else propagation.get("target_profile")),
                        "target_loss_fraction_under_declared_model": (
                            1.0 - expanded_baseline["horizon_parent_fraction"] if expanded_baseline
                            else _model_loss(repaired_evidence)),
                        "absolute_retention_approximation_error_bound": (
                            conditional_certificates[0].get("maximum_abs_parent_retention_error")
                            if len(conditional_inserted) == 1 else None
                        ),
                        "conditional_model_route_ids": [item.get("route_id") for item in conditional_inserted],
                        "omission_bounds": [
                            {"route_id": item.get("route_id"),
                             "absolute_retention_approximation_error_bound": certificate.get(
                                 "maximum_abs_parent_retention_error"),
                             "assessment": (certificate.get("expanded_model_assessment") or
                                            item["forward_only_approximation_evidence"].get("assessment"))}
                            for item, certificate in zip(conditional_inserted, conditional_certificates)
                        ],
                        "elementary_rate_evidence": [
                            {"route_id": item.get("route_id"),
                             "condition_rate": item.get("condition_rate"),
                             "conformer_rate_basis": item.get("conformer_rate_basis"),
                             "conditional_forward_rate_evidence": item.get("conditional_forward_rate_evidence"),
                             "forward_only_approximation_evidence": item["forward_only_approximation_evidence"]}
                            for item in conditional_inserted
                        ],
                        "repaired_model_inventory_certificates": conditional_certificates,
                        "expanded_network_assessments": expanded_network_reports,
                        "next_evidence_requests": [request for report in expanded_network_reports
                                                   for request in report["next_evidence_requests"]],
                        "pressure_applicability_verified": False,
                        "rate_uncertainty_resolved": False,
                        "conformer_completeness_verified": False,
                        "network_completeness_verified": False,
                        "conclusive_observable_allowed": False,
                        "lifetime_claim_allowed": False,
                        "artifacts": {
                            "rerun_result": iteration["artifacts"]["rerun_result"],
                            "repaired_model_validation": iteration["artifacts"]["repaired_model_validation"],
                            "conditional_model_inventory_validation": iteration["artifacts"]["conditional_model_inventory_validation"],
                        },
                    }
                    iteration.update(status="completed_conditional_model",
                                     observable_claim_allowed=False,
                                     rate_verification_status="conditional_model_rate_calculated")
                    iteration["artifacts"]["conditional_model_prediction"] = _write_json(
                        iteration_dir / "conditional-model-prediction.json", conditional_prediction,
                    )
                    return self._summary(
                        "conditional_model_prediction",
                        reason=("Forward-only repaired model was propagated within its certified "
                                "finite-horizon omission scope; physical applicability and rate "
                                "uncertainty remain unresolved, so no lifetime is claimed."),
                        evidence=repaired_evidence, plan=selection_plan,
                        conditional_prediction=conditional_prediction,
                    )
            except Exception as error:
                if multistep_stage == "specialist_rate_generation" and path_result is not None:
                    retained_status = _unverified_path_status(path_result)
                    self.deferred_routes[route_id] = {
                        "route_id": route_id,
                        "status": retained_status,
                        "reason": _unverified_path_detail(path_result),
                        "path_result": iteration["artifacts"].get("path_result"),
                        "specialist_rate_retry_eligible": True,
                    }
                precise_state = (
                    "verified_multistep_path_network_insertion_pending"
                    if multistep_stage == "network_insertion_and_repropagation" else
                    "segment_rate_model_unavailable"
                    if multistep_stage == "segment_rate_generation" else
                    "specialist_rate_model_unavailable"
                    if multistep_stage == "specialist_rate_generation" else
                    "verification_iteration_failed"
                )
                iteration.update(
                    status="failed",
                    scientific_failure_state=precise_state,
                    failure_reason=f"{type(error).__name__}: {error}",
                )
                iteration["artifacts"]["failure"] = _write_json(iteration_dir / "failure.json", {
                    "stage": "route_path_rate_or_model_rerun",
                    "scientific_failure_state": precise_state,
                    "error_type": type(error).__name__,
                    "message": str(error),
                })
                self._persist_state("failed_closed")
                return self._incomplete(precise_state, evidence=current_evidence,
                                        plan=selection_plan, error=error)

        return self._incomplete("verification_engine_reached_unexpected_terminal_state", evidence=current_evidence)


def run_verification_engine(
    initial_rmg_evidence: dict,
    condition_contract: dict,
    workdir: Path,
    *,
    verify_route: VerifyRoute,
    build_rate: BuildRate,
    rerun_model: RerunModel,
    config: VerificationEngineConfig | None = None,
    supplemental_candidate_routes: Iterable[dict] | None = None,
    build_specialist_rate: BuildSpecialistRate | None = None,
    rate_uncertainties: Iterable[dict | RateUncertaintyRecord] = (),
    rmg_env: str | None = None,
    path_verifier_policy_id: str = "generic-orca-path-v1",
) -> dict:
    """Convenience API for a fresh retained verification-engine session."""
    return VerificationEngine(
        initial_rmg_evidence, condition_contract, workdir,
        verify_route=verify_route,
        build_rate=build_rate,
        rerun_model=rerun_model,
        config=config,
        supplemental_candidate_routes=supplemental_candidate_routes,
        build_specialist_rate=build_specialist_rate,
        rate_uncertainties=rate_uncertainties,
        rmg_env=rmg_env,
        path_verifier_policy_id=path_verifier_policy_id,
    ).run()


def continue_verification_engine(
    initial_rmg_evidence: dict,
    condition_contract: dict,
    workdir: Path,
    predecessor_state: Path,
    *,
    verify_route: VerifyRoute,
    build_rate: BuildRate,
    rerun_model: RerunModel,
    path_verifier_policy_id: str,
    retry_deferred_route_ids: Iterable[str] | None = None,
    recover_in_progress_route_ids: Iterable[str] | None = None,
    execution_recovery_evidence: dict | None = None,
    scientific_policy_amendment: dict | None = None,
    config: VerificationEngineConfig | None = None,
    build_specialist_rate: BuildSpecialistRate | None = None,
    retained_specialist_path_results: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict:
    """Continue a failed-closed session under the exact retained contract.

    Engine records remain append-only.  Calculation workspaces may be resumed
    in place so their program-level checkpoints and validated endpoint caches
    remain usable; the new session records every result and the predecessor
    state hash.  A scientifically deferred route can be retried only once per
    explicit path-verifier policy ID.
    """
    state_path = Path(predecessor_state)
    if not state_path.is_file():
        raise FileNotFoundError(f"Verification engine state not found: {state_path}")
    state = json.loads(state_path.read_text())
    if state.get("kind") != "generic_verification_engine_state":
        raise ValueError("Predecessor is not a generic verification engine state")
    retained_contract = state.get("condition_contract")
    retained_hash = state.get("condition_contract_sha256")
    if not isinstance(retained_contract, dict) or retained_hash != _condition_hash(retained_contract):
        raise ValueError("Predecessor condition contract is missing or has been mutated")
    if _condition_hash(condition_contract) != retained_hash:
        raise ValueError("Continuation condition contract does not exactly match the predecessor")

    if retry_deferred_route_ids is not None and recover_in_progress_route_ids is not None:
        raise ValueError("A continuation cannot be both a policy retry and an execution recovery")

    deferred = {
        str(route_id): dict(item)
        for route_id, item in (state.get("deferred_unverified_routes") or {}).items()
    }
    execution_recovery = recover_in_progress_route_ids is not None
    requested = tuple(str(item) for item in (
        recover_in_progress_route_ids if execution_recovery
        else retry_deferred_route_ids if retry_deferred_route_ids is not None
        else deferred
    ))
    if not requested:
        raise ValueError("Predecessor has no deferred route eligible for continuation")
    if execution_recovery:
        recovery_kind = (execution_recovery_evidence or {}).get("kind")
        if recovery_kind in {"orca-6.1-neb-preopt-parser-repair-v1", "orca-neb-time-limit-checkpoint-continuation-v1",
                             "orca-6.1-neb-restart-preopt-parser-repair-v1"}:
            # Reopen the exact failed job proof, rather than trusting a caller's
            # recovery label as permission to repeat a chemical failure.
            from .verification_resume import (_recoverable_neb_preopt_parser_failure, _recoverable_neb_time_limit_checkpoint,
                                              _recoverable_neb_restart_preopt_parser_failure)
            inspector = {"orca-6.1-neb-preopt-parser-repair-v1": _recoverable_neb_preopt_parser_failure,
                         "orca-neb-time-limit-checkpoint-continuation-v1": _recoverable_neb_time_limit_checkpoint,
                         "orca-6.1-neb-restart-preopt-parser-repair-v1": _recoverable_neb_restart_preopt_parser_failure}[recovery_kind]
            actual_recovery = inspector(predecessor_state, state)
            if (not actual_recovery
                    or not set(requested).issubset(actual_recovery["route_ids"])
                    or any((execution_recovery_evidence or {}).get(key) != value
                           for key, value in actual_recovery.items())):
                raise ValueError("NEB execution recovery lacks the exact current source-bound failure proof")
        latest_by_route = {}
        for iteration in state.get("iterations") or []:
            route_id = str(iteration.get("route_id") or "")
            if route_id:
                latest_by_route[route_id] = iteration
        ineligible = sorted(
            route_id for route_id in requested
            if route_id not in latest_by_route
            or (
                latest_by_route[route_id].get("status") != "started"
                and not (
                    (execution_recovery_evidence or {}).get("kind") in {
                        "orca_startup_infrastructure_execution_failure",
                        "surface_budget_policy_serialization_execution_failure",
                        "surface_context_policy_amendment_execution_failure",
                        "surface_attempt_policy_context_identity_execution_failure",
                        "orca-6.1-neb-preopt-parser-repair-v1",
                        "orca-neb-time-limit-checkpoint-continuation-v1",
                        "orca-6.1-neb-restart-preopt-parser-repair-v1",
                    }
                    and latest_by_route[route_id].get("status") in {
                        "deferred_unverified", "failed",
                    }
                )
            )
        )
        if ineligible:
            raise ValueError(
                "Execution recovery is restricted to retained started routes: "
                + ", ".join(ineligible)
            )
    else:
        missing = sorted(set(requested) - set(deferred))
        if missing:
            raise ValueError(
                "Requested retry routes are not deferred in the predecessor: "
                + ", ".join(missing)
            )

    prior_policy = str(state.get("path_verifier_policy_id") or "legacy-unversioned-path-policy")
    attempts = {
        str(route_id): [str(policy) for policy in policies]
        for route_id, policies in (state.get("route_policy_attempts") or {}).items()
    }
    for iteration in state.get("iterations") or []:
        route_id = str(iteration.get("route_id") or "")
        if route_id:
            policy = str(iteration.get("path_verifier_policy_id") or prior_policy)
            if policy not in attempts.setdefault(route_id, []):
                attempts[route_id].append(policy)
    current_policy = str(path_verifier_policy_id).strip()
    if not current_policy:
        raise ValueError("Continuation needs an explicit path-verifier policy ID")
    specialist_continuation = bool(build_specialist_rate) and not execution_recovery and all(
        str((deferred.get(route_id) or {}).get("status")) in {
            "barrierless_rate_model_required", "dissociation_rate_model_required",
        }
        for route_id in requested
    )
    repeated = sorted(route_id for route_id in requested if current_policy in attempts.get(route_id, ()))
    if repeated and not execution_recovery and not specialist_continuation:
        raise ValueError(
            "Deferred routes were already attempted under the current path-verifier policy: "
            + ", ".join(repeated)
        )

    workspaces: dict[str, Path] = {}
    retained_paths: dict[str, dict[str, Any]] = {
        str(route_id): dict(record)
        for route_id, record in (retained_specialist_path_results or {}).items()
    }
    for route_id in requested:
        matching = [
            item for item in (state.get("iterations") or [])
            if str(item.get("route_id")) == route_id
        ]
        if not matching:
            raise ValueError(f"Deferred route {route_id} has no retained iteration workspace")
        raw_directory = (matching[-1].get("artifacts") or {}).get("calculation_workspace") \
            or (matching[-1].get("artifacts") or {}).get("directory")
        workspace = _resolve_retained_path(raw_directory, state_path=state_path)
        if not workspace.is_dir():
            raise FileNotFoundError(f"Retained calculation workspace is missing: {workspace}")
        workspaces[route_id] = workspace
        if specialist_continuation and route_id not in retained_paths:
            raw_path_result = (deferred[route_id].get("path_result")
                               or (matching[-1].get("artifacts") or {}).get("path_result"))
            path_result_path = _resolve_retained_path(raw_path_result, state_path=state_path)
            if not path_result_path.is_file():
                raise FileNotFoundError(
                    f"Retained specialist path result is missing: {path_result_path}"
                )
            retained_paths[route_id] = {
                "path_result": json.loads(path_result_path.read_text()),
                "source_path": str(path_result_path),
                "source_sha256": _file_sha256(path_result_path),
            }
        deferred.pop(route_id, None)

    replacements = [dict(item) for item in (state.get("verified_replacements") or [])]
    library_paths = [
        _resolve_retained_path(item, state_path=state_path)
        for item in (state.get("library_paths") or [])
    ]
    missing_libraries = [str(path) for path in library_paths if not path.exists()]
    if missing_libraries:
        raise FileNotFoundError("Retained verified kinetics libraries are missing: " + ", ".join(missing_libraries))

    model_result: dict = initial_rmg_evidence
    if replacements:
        completed = [
            item for item in (state.get("iterations") or [])
            if item.get("status") == "completed" and (item.get("artifacts") or {}).get("rerun_result")
        ]
        if not completed:
            raise ValueError("Predecessor has replacements but no retained repaired-model result")
        rerun_path = _resolve_retained_path(
            completed[-1]["artifacts"]["rerun_result"], state_path=state_path,
        )
        if not rerun_path.is_file():
            raise FileNotFoundError(f"Retained repaired-model result is missing: {rerun_path}")
        model_result = json.loads(rerun_path.read_text())
        _, model_contract = _extract_model_result(model_result)
        if not isinstance(model_contract, dict) or _condition_hash(model_contract) != retained_hash:
            raise ValueError("Retained repaired-model result does not match the immutable condition")

    uncertainty = (state.get("rate_uncertainty_evidence") or {}).get("records") or ()
    continuation = {
        "kind": (
            "verification_engine_stale_execution_recovery" if execution_recovery
            else "verification_engine_specialist_rate_continuation" if specialist_continuation
            else "verification_engine_continuation"
        ),
        "predecessor_state": str(state_path),
        "predecessor_state_sha256": _file_sha256(state_path),
        "predecessor_status": state.get("status"),
        "prior_path_verifier_policy_id": prior_policy,
        "path_verifier_policy_id": current_policy,
        "retried_route_ids": [] if execution_recovery else list(requested),
        "recovered_in_progress_route_ids": list(requested) if execution_recovery else [],
        "scientific_policy_attempt_reused": execution_recovery or specialist_continuation,
        "path_verifier_reexecuted": not specialist_continuation,
        "scientific_policy_amendment": _jsonable(scientific_policy_amendment),
        "execution_recovery_evidence": (
            _jsonable(execution_recovery_evidence) if execution_recovery else None
        ),
    }
    return VerificationEngine(
        model_result, condition_contract, workdir,
        verify_route=verify_route,
        build_rate=build_rate,
        rerun_model=rerun_model,
        config=config,
        supplemental_candidate_routes=state.get("supplemental_candidate_routes") or (),
        build_specialist_rate=build_specialist_rate,
        rate_uncertainties=uncertainty,
        rmg_env=state.get("rmg_environment"),
        initial_replacements=replacements,
        initial_library_paths=library_paths,
        initial_deferred_routes=deferred,
        calculation_workspaces=workspaces,
        path_verifier_policy_id=current_policy,
        route_policy_attempts=attempts,
        same_policy_execution_recovery_route_ids=requested if execution_recovery else (),
        retained_specialist_path_results=retained_paths,
        continuation_evidence=continuation,
    ).run()


__all__ = [
    "BuildRate", "BuildSpecialistRate", "RerunModel", "VerifyRoute", "VerificationEngine", "VerificationEngineConfig",
    "continue_verification_engine", "run_verification_engine",
]
