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
from typing import Any, Callable, Iterable

from .flux_verification import (
    DEFAULT_FLUX_COVERAGE_TARGET,
    DEFAULT_MAX_FLUX_CLOSURE_ERROR,
    DEFAULT_T95_RELATIVE_TOLERANCE,
    assess_repaired_mechanism,
    assess_t95_robustness,
    assess_target_loss_commitment,
    build_flux_ranked_verification_plan,
    reaction_signature,
    _retained_network_routes,
    _solver_gates,
)
from .generated_kinetics import validate_rate_replacement_evidence


VerifyRoute = Callable[[dict, Path], dict]
BuildRate = Callable[[dict, dict, Path], dict]
RerunModel = Callable[[list[Path], Path], dict]

_PATH_VERIFIED_STATUSES = {
    "verified", "orca_verified", "irc_endpoint_verified", "barrierless_capture_verified",
    "barriered_ts_verified", "verified_barrierless_capture", "verified_irc_endpoints",
    "verified_transition_state_path",
}
_RATE_VERIFIED_STATUSES = {
    "validated", "verified", "rate_verified", "arkane_verified", "collision_bounded_verified",
    "verified_rate", "verified_arkane", "verified_collision_bounded",
}
_DEFERABLE_PATH_STATUSES = {
    "ambiguous_atom_mapping",
    "barrierless_rate_model_required",
    "endpoint_validation_failed",
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
                    "combined_unresolved_route_sensitivity"):
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


def _route_for_plan(evidence: dict, plan: dict) -> dict:
    source = plan.get("selected_route_source") or "candidate_routes"
    routes = ((evidence.get("candidate_routes") or []) if source == "candidate_routes"
              else _retained_network_routes(evidence))
    index = plan.get("selected_route_index")
    if not isinstance(index, int) or index < 0 or index >= len(routes):
        raise ValueError("Flux plan selected no resolvable candidate route")
    route = dict(routes[index])
    route.setdefault("route_id", plan.get("selected_route_id"))
    if str(route.get("route_id")) != str(plan.get("selected_route_id")):
        raise ValueError("Flux plan route index and route ID disagree")
    return route


def _normalize_replacement(rate_result: dict, route: dict) -> tuple[dict, list[Path]]:
    if not isinstance(rate_result, dict):
        raise TypeError("Rate callback must return a dictionary")
    block = _nested_dict(rate_result, (
        "replacement", "rate_evidence", "verified_rate", "generated_kinetics_library", "kinetics",
    )) or rate_result
    replacement = dict(block)
    replacement.setdefault("route_id", route.get("route_id"))
    replacement.setdefault("reaction_equation", route.get("reaction_equation") or route.get("printed_reaction_equation"))
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


def _planning_replacements(replacements: list[dict], evidence: dict) -> list[dict]:
    """Rebind verified reaction signatures if RMG reordered candidate IDs."""
    routes = list(evidence.get("candidate_routes") or [])
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
        signature = reaction_signature(
            replacement.get("reaction_equation") or replacement.get("replaces_reaction") or "", normalized=True,
        )
        matches = [route for route in routes if reaction_signature(
            route.get("reaction_equation") or route.get("printed_reaction_equation") or "", normalized=True,
        ) == signature]
        item = dict(replacement)
        if len(matches) == 1:
            item["route_id"] = matches[0].get("route_id")
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
    ) -> None:
        if not isinstance(condition_contract, dict) or not condition_contract:
            raise ValueError("Verification requires a complete immutable condition contract")
        retention = _finite(condition_contract.get("retention_fraction"))
        if retention is None or not 0.0 < retention < 1.0:
            raise ValueError("Condition contract retention_fraction must be in (0, 1)")
        self.initial_result = initial_rmg_evidence
        self.condition_contract = _jsonable(condition_contract)
        self.condition_sha256 = _condition_hash(self.condition_contract)
        self.session = _unique_session(Path(workdir))
        self.verify_route_callback = verify_route
        self.build_rate_callback = build_rate
        self.rerun_model_callback = rerun_model
        self.config = config or VerificationEngineConfig()
        self.replacements: list[dict] = []
        self.library_paths: list[Path] = []
        self.iterations: list[dict] = []
        self.deferred_routes: dict[str, dict] = {}

    def _persist_state(self, status: str, **extra: Any) -> None:
        _write_json(self.session / "engine-state.json", {
            "schema_version": 1,
            "kind": "generic_verification_engine_state",
            "status": status,
            "condition_contract": self.condition_contract,
            "condition_contract_sha256": self.condition_sha256,
            "verified_replacements": self.replacements,
            "library_paths": self.library_paths,
            "iterations": self.iterations,
            "deferred_unverified_routes": self.deferred_routes,
            **extra,
        })

    def _summary(self, status: str, *, reason: str, evidence: dict | None = None,
                 plan: dict | None = None, robustness: dict | None = None,
                 envelope: dict | None = None) -> dict:
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
        model_only_no_loss = status == "no_target_loss_in_completed_rmg_model"
        if model_only_no_loss:
            verified_coverage = None
        if model_only_no_loss:
            route_status = rate_status = "not_required_for_retained_model"
            propagation_status = "completed_model"
        elif self.deferred_routes:
            route_status = "deferred_unverified"
            rate_status = "unavailable_for_deferred_routes"
            propagation_status = "screening_model_only"
        else:
            route_status = "verified" if self.replacements else "not_required"
            rate_status = "validated" if self.replacements else "not_required"
            propagation_status = ("verified_converged" if status in {
                "verified_t95_converged", "verified_below_loss_threshold",
            } else "incomplete")
        if model_only_no_loss:
            claim_scope = "completed_retained_rmg_model_only"
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
            "verified_route_ids": [item.get("route_id") for item in self.replacements],
            "deferred_unverified_routes": list(self.deferred_routes.values()),
            "iterations_completed": len(self.iterations),
            "robustness": robustness,
            "artifacts": {
                "session_directory": str(self.session),
                "engine_state": str(self.session / "engine-state.json"),
                "verification_summary": str(self.session / "verification-summary.json"),
            },
            "screening_estimate_prohibited": status in {
                "verification_incomplete", "no_target_loss_in_completed_rmg_model",
            },
            "claim_scope": claim_scope,
            "model_coverage": propagation.get("coverage"),
        }
        _write_json(self.session / "verification-summary.json", summary)
        self._persist_state("completed" if status != "verification_incomplete" else "failed_closed",
                            verification_summary=summary)
        return summary

    def _incomplete(self, reason: str, *, evidence: dict | None = None,
                    plan: dict | None = None, error: Exception | None = None) -> dict:
        detail = reason if error is None else f"{reason}: {type(error).__name__}: {error}"
        return self._summary("verification_incomplete", reason=detail, evidence=evidence, plan=plan)

    def _terminal_if_robust(self, evidence: dict, plan: dict, envelope: dict | None) -> dict | None:
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
                and not _solver_gates(
                    current_evidence, max_closure_error=self.config.max_flux_closure_error,
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
            planning_replacements = _planning_replacements(self.replacements, current_evidence)
            plan = build_flux_ranked_verification_plan(
                current_evidence, self.condition_contract,
                verified_replacements=planning_replacements,
                deferred_route_ids=self.deferred_routes,
                flux_coverage_target=self.config.flux_coverage_target,
                max_closure_error=self.config.max_flux_closure_error,
            )
            _write_json(self.session / f"plan-{iteration_number:03d}.json", plan)
            if plan.get("status") in {
                "insufficient_solver_evidence", "insufficient_route_attribution",
                "controlling_route_has_unresolved_initiation", "partial_enumeration_no_direct_route",
            }:
                return self._incomplete(f"flux_plan_{plan.get('status')}", evidence=current_evidence, plan=plan)

            envelope = _sensitivity_envelope(current_result, current_evidence)
            # A direct unimolecular RMG route is deliberately queued after an
            # invalid kinetic screen even when its estimated flux is zero.  It
            # is not eligible for the usual no-loss robustness shortcut until
            # ORCA has accepted or deferred that actual RMG hypothesis.
            terminal = None if plan.get("intrinsic_failover_required") else self._terminal_if_robust(
                current_evidence, plan, envelope,
            )
            if terminal is not None:
                return terminal

            if iteration_number > self.config.max_iterations:
                return self._incomplete("maximum_verification_iterations_reached",
                                        evidence=current_evidence, plan=plan)

            # Reaching the configured coverage is not terminal by itself.  If
            # no physical sensitivity envelope proves robustness, continue to
            # the next nonzero-flux route until coverage is complete.
            selection_plan = plan
            if plan.get("status") in {"flux_coverage_target_met", "no_target_destruction_flux"}:
                selection_plan = build_flux_ranked_verification_plan(
                    current_evidence, self.condition_contract,
                    verified_replacements=planning_replacements,
                    deferred_route_ids=self.deferred_routes,
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
                route = _route_for_plan(current_evidence, selection_plan)
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
                path_result = self.verify_route_callback(route, iteration_dir / "path")
                iteration["artifacts"]["path_result"] = _write_json(iteration_dir / "path-result.json", path_result)
                path_status = _verified_status(path_result, _PATH_VERIFIED_STATUSES, nested_keys=(
                    "path_result", "route_verification", "verification", "orca_verification",
                ))
                if path_status is None:
                    unverified_status = _unverified_path_status(path_result)
                    if unverified_status in _DEFERABLE_PATH_STATUSES:
                        detail = _unverified_path_detail(path_result)
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
                        }
                        self._persist_state("running_with_deferred_routes")
                        continue
                    raise ValueError(
                        "Route callback did not return an explicitly verified path status: "
                        + _unverified_path_detail(path_result)
                    )

                rate_result = self.build_rate_callback(route, path_result, iteration_dir / "rate")
                iteration["artifacts"]["rate_result"] = _write_json(iteration_dir / "rate-result.json", rate_result)
                replacement, new_paths = _normalize_replacement(rate_result, route)
                validation = validate_rate_replacement_evidence(
                    replacement,
                    temperature_K=float(self.condition_contract["temperature_K"]),
                    pressure_bar=float(self.condition_contract["pressure_bar"]),
                )
                iteration["artifacts"]["replacement_validation"] = _write_json(
                    iteration_dir / "replacement-validation.json", validation,
                )
                if validation.get("status") != "passed":
                    raise ValueError("Replacement rate failed validation: " + ", ".join(validation.get("blocking_reasons") or []))
                self.replacements.append(replacement)
                for path in new_paths:
                    if path not in self.library_paths:
                        self.library_paths.append(path)

                rerun_result = self.rerun_model_callback(list(self.library_paths), iteration_dir / "rerun")
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
                iteration.update(status="completed", path_verification_status=path_status,
                                 replacement_route_id=replacement.get("route_id"))
                current_evidence = repaired_evidence
                current_result = rerun_result
                self._persist_state("running")
            except Exception as error:
                iteration.update(status="failed", failure_reason=f"{type(error).__name__}: {error}")
                iteration["artifacts"]["failure"] = _write_json(iteration_dir / "failure.json", {
                    "stage": "route_path_rate_or_model_rerun",
                    "error_type": type(error).__name__,
                    "message": str(error),
                })
                self._persist_state("failed_closed")
                return self._incomplete("verification_iteration_failed", evidence=current_evidence,
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
) -> dict:
    """Convenience API for a fresh retained verification-engine session."""
    return VerificationEngine(
        initial_rmg_evidence, condition_contract, workdir,
        verify_route=verify_route, build_rate=build_rate, rerun_model=rerun_model, config=config,
    ).run()


__all__ = [
    "BuildRate", "RerunModel", "VerifyRoute", "VerificationEngine", "VerificationEngineConfig",
    "run_verification_engine",
]
