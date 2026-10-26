"""Bounded, observable-conditioned verification sessions for solution t95."""

from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping

from .solution_allocation import allocate_next_solution_verification
from .solution_candidate_adapters import candidate_network_from_session_request
from .solution_candidates import SolutionCandidateNetwork
from .solution_kinetics import VerifiedSolutionNetwork, VerifiedSolutionReaction
from .solution_observables import (
    assess_solution_session_convergence, evaluate_solution_t95_bounds,
)
from .solution_path_verification import run_implicit_solution_path_verification
from .solution_path_verification import THERMOCHEMISTRY_TEMPERATURE_CONTRACT
from .orca_thermochemistry_recipe import THERMOCHEMISTRY_RECIPE_CONTRACT
from .conformer_rates import CONFORMER_RATE_CONTRACT
from .solution_stability import run_implicit_solution_stability_bundle


PathVerifier = Callable[..., dict[str, Any]]


def _write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


def _canonical_hash(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False, default=str)
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _environment_dict(network: SolutionCandidateNetwork) -> dict[str, Any]:
    environment = network.environment
    return {
        "solvent": environment.solvent, "solvent_model": environment.solvent_model,
        "temperature_K": environment.temperature_K,
        "target_duration_seconds": environment.target_duration_seconds,
        "retention_fraction": environment.retention_fraction,
        "standard_state_mol_L": environment.standard_state_mol_L,
        "volume_L": environment.volume_L, "light_condition": environment.light_condition,
        "electronic_structure_profile": environment.electronic_structure_profile,
    }


def _validate_session_policy(request: Mapping[str, Any], network: SolutionCandidateNetwork) -> dict[str, Any]:
    if str(request.get("observable") or "t95") != "t95":
        raise ValueError("The first solution verification session supports observable='t95' only")
    convergence = dict(request.get("convergence") or {})
    budget = dict(request.get("compute_budget") or {})
    try:
        ratio = float(convergence["maximum_t95_ratio"])
        max_iterations = int(convergence["maximum_iterations"])
        max_paths = int(budget["maximum_orca_paths"])
        max_core_hours = float(budget["maximum_estimated_core_hours"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "Session requires maximum_t95_ratio, maximum_iterations, maximum_orca_paths, "
            "and maximum_estimated_core_hours"
        ) from error
    if not math.isfinite(ratio) or ratio < 1.0 or max_iterations < 1 or max_paths < 0:
        raise ValueError("Invalid convergence or ORCA-path budget")
    if not math.isfinite(max_core_hours) or max_core_hours < 0.0:
        raise ValueError("maximum_estimated_core_hours must be finite and nonnegative")
    profile = str(request.get("method_profile") or "").strip()
    if not profile or profile != network.environment.electronic_structure_profile:
        raise ValueError("method_profile must match environment.electronic_structure_profile exactly")
    # The path producer now explicitly requests and proves single-T thermal
    # corrections. This is not calibration of named solvent parameters at T.
    from .methods import harmonic_method_profile
    harmonic_method_profile(profile)
    return {
        "maximum_t95_ratio": ratio,
        "maximum_iterations": max_iterations,
        "maximum_orca_paths": max_paths,
        "maximum_estimated_core_hours": max_core_hours,
        "require_same_horizon_classification": bool(
            convergence.get("require_same_horizon_classification", True)
        ),
        "method_profile": profile,
        "thermochemistry_temperature_K": network.environment.temperature_K,
        "thermochemistry_temperature_contract": THERMOCHEMISTRY_TEMPERATURE_CONTRACT,
        "thermochemistry_recipe_contract": THERMOCHEMISTRY_RECIPE_CONTRACT,
        "conformer_rate_contract": CONFORMER_RATE_CONTRACT,
        "solvent_temperature_treatment": "fixed_ORCA_named_solvent_parameters",
    }


def _cache_key(candidate, network: SolutionCandidateNetwork, method_profile: str) -> str:
    return _canonical_hash({
        "candidate_id": candidate.candidate_id,
        "reactant_id": candidate.reactant_id,
        "products": candidate.products,
        "route": candidate.route,
        "environment_contract_id": network.environment.as_environment_spec().as_dict()["contract_id"],
        "method_profile": method_profile,
        "thermochemistry_temperature_contract": THERMOCHEMISTRY_TEMPERATURE_CONTRACT,
        "thermochemistry_recipe_contract": THERMOCHEMISTRY_RECIPE_CONTRACT,
        "conformer_rate_contract": CONFORMER_RATE_CONTRACT,
        "uncertainty": candidate.uncertainty,
        "schema_version": 1,
    })


def _scientific_production_policy(policy: Mapping[str, Any], network: SolutionCandidateNetwork) -> dict[str, Any]:
    """Credit native session-produced rates under an explicit producer contract.

    This does not validate underlying calculations or recredit imported rates.
    It prevents a newer producer from silently trusting its older session state.
    """
    return {
        "kind": "solution_session_scientific_production_policy", "schema_version": 1,
        **{key: policy[key] for key in (
            "method_profile", "thermochemistry_temperature_K",
            "thermochemistry_temperature_contract", "thermochemistry_recipe_contract",
            "conformer_rate_contract", "solvent_temperature_treatment",
        )},
        "environment_contract_id": network.environment.as_environment_spec().as_dict()["contract_id"],
    }


def _rate_evidence_origins(
    network: SolutionCandidateNetwork, initial_network: SolutionCandidateNetwork,
) -> tuple[dict[str, str], list[str]]:
    # Import identity comes from the immutable original request, not an editable
    # state annotation or a reaction source string that happens to say STORCA.
    imported = {
        item.reaction_id: _canonical_hash(item.as_dict())
        for item in initial_network.verified_reactions
    }
    origins, session_tst_ids = {}, []
    for reaction in network.verified_reactions:
        if imported.get(reaction.reaction_id) == _canonical_hash(reaction.as_dict()):
            origins[reaction.reaction_id] = "user_declared_imported"
        elif reaction.verification_status == "verified_solution_tst" or reaction.rate_theory == "solution_eyring_tst":
            origins[reaction.reaction_id] = "session_produced_tst"
            session_tst_ids.append(reaction.reaction_id)
        else:
            origins[reaction.reaction_id] = "session_produced_other_rate"
    return origins, sorted(session_tst_ids)


def _path_request(network: SolutionCandidateNetwork, candidate) -> dict[str, Any]:
    unresolved_for_path = sorted({
        *network.unresolved_discovery_ids,
        *(item.candidate_id for item in network.candidates if item.candidate_id != candidate.candidate_id),
    })
    return {
        "schema_version": 1,
        "kind": "implicit_solution_path_verification_request",
        "target_smiles": network.target_smiles,
        "environment": _environment_dict(network),
        "network": {
            "target_species_id": network.target_species_id,
            "coverage_status": (
                "bounded_search_incomplete" if unresolved_for_path else network.coverage_status
            ),
            "coverage_claim_id": network.coverage_claim_id,
            "discovery_sources": list(network.discovery_sources),
            "unresolved_candidate_ids": unresolved_for_path,
            "species": [item.as_dict() for item in network.species],
            "reactions": [item.as_dict() for item in network.verified_reactions],
        },
        "candidate_reaction": {
            "reaction_id": candidate.candidate_id,
            "reactant_id": candidate.reactant_id,
            "products": dict(candidate.products), "molecularity": 1,
            "applicability_flags": dict(candidate.applicability_flags),
            "uncertainty": dict(candidate.uncertainty),
        },
        "route": dict(candidate.route),
    }


def _terminal_status(network, convergence, allocation, used_paths, policy) -> tuple[str, str]:
    target_escalation = [
        item.candidate_id for item in network.candidates
        if item.reactant_id == network.target_species_id
        and (item.active_applicability_flags or item.verification_state == "method_escalation_required")
    ]
    target_unresolved = [
        item.candidate_id for item in network.candidates
        if item.reactant_id == network.target_species_id and item.verification_state == "path_unresolved"
    ]
    if target_escalation:
        return "method_escalation_required", "Influential target-loss chemistry requires a specialized method."
    if target_unresolved:
        return "path_verification_unresolved", "At least one influential target-loss path remains unresolved."
    if network.coverage_status != "verified_complete_within_declared_unimolecular_scope" or network.unresolved_discovery_ids:
        return "network_coverage_unresolved", "Reaction discovery does not support the restricted coverage claim."
    if used_paths >= policy["maximum_orca_paths"]:
        return "verification_budget_exhausted", "The bounded ORCA-path budget was exhausted."
    if allocation.get("selected_candidate_id") is None:
        return "rate_uncertainty_unresolved", "No eligible ordinary path can resolve the remaining t95 uncertainty."
    return "not_converged", "The requested convergence contract has not been satisfied."


def _verified_bundle(
    network: SolutionCandidateNetwork,
    claim: Mapping[str, Any],
    workdir: Path,
    retained_path_claims: list[dict[str, Any]],
) -> tuple[str | None, dict[str, Any] | None, str | None]:
    if network.candidates:
        return None, None, (
            "A t95-only claim converged, but non-target candidate reactions remain; STORCA will not "
            "emit a general concentration/product bundle that omits their kinetics."
        )
    verified = VerifiedSolutionNetwork(
        target_species_id=network.target_species_id, species=network.species,
        reactions=network.verified_reactions, coverage_status=network.coverage_status,
        coverage_claim_id=network.coverage_claim_id,
        discovery_sources=network.discovery_sources,
        unresolved_candidate_ids=network.unresolved_discovery_ids,
    )
    bundle = {
        "schema_version": 1, "kind": "verified_implicit_solution_stability_bundle",
        "target_smiles": network.target_smiles, "environment": _environment_dict(network),
        "network": verified.as_dict(),
        "scientific_claims": [*retained_path_claims, dict(claim)],
        "producer": {"kind": "storca_solution_verification_session", "schema_version": 1},
    }
    bundle_path = Path(workdir) / "verified-solution-bundle.json"
    _write_json(bundle_path, bundle)
    propagation = run_implicit_solution_stability_bundle(bundle, Path(workdir) / "propagation")
    return str(bundle_path), propagation, None


def _retained_path_claims(workdir: Path) -> list[dict[str, Any]]:
    claims: list[dict[str, Any]] = []
    seen: set[str] = set()
    for path in sorted(Path(workdir).glob("iteration-*/path-result.json")):
        try:
            result = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        for key in ("path_claim", "rate_claim"):
            claim = result.get(key)
            if not isinstance(claim, dict):
                continue
            identity = str(claim.get("claim_id") or _canonical_hash(claim))
            if identity not in seen:
                claims.append(claim)
                seen.add(identity)
    return claims


def run_solution_verification_session(
    request: Mapping[str, Any],
    workdir: Path,
    *,
    execute: bool = False,
    ncores: int = 1,
    path_verifier: PathVerifier = run_implicit_solution_path_verification,
) -> dict[str, Any]:
    """Run or prepare one append-only t95-conditioned verification session."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    request_path = workdir / "session-request.json"
    state_path = workdir / "session-state.json"
    result_path = workdir / "solution-verification-session.json"
    try:
        if ncores < 1:
            raise ValueError("ncores must be at least one")
        request_hash = _canonical_hash(request)
        if request_path.is_file():
            retained = json.loads(request_path.read_text())
            if _canonical_hash(retained) != request_hash:
                return {
                    "schema_version": 1, "kind": "solution_verification_session",
                    "status": "invalid_or_incomplete_evidence",
                    "reason": "Session directory contains a different immutable request contract",
                    "verified_bundle": None, "result_json": None,
                    "retained_result_json": str(result_path) if result_path.is_file() else None,
                    "requested_request_sha256": request_hash,
                    "retained_request_sha256": _canonical_hash(retained),
                }
        else:
            _write_json(request_path, dict(request))
        initial_network = candidate_network_from_session_request(request)
        retained_state = None
        if state_path.is_file():
            state = json.loads(state_path.read_text())
            retained_state = state
            network = SolutionCandidateNetwork.from_dict(state["candidate_network"])
            used_paths = int(state.get("used_orca_paths", 0))
            used_core_hours = float(state.get("used_estimated_core_hours", 0.0))
            completed_iterations = int(state.get("completed_iterations", 0))
            artifact_index = int(state.get("next_artifact_index", len(list(workdir.glob("iteration-*")))))
        else:
            network = initial_network
            used_paths = 0
            used_core_hours = 0.0
            completed_iterations = 0
            artifact_index = 0
        policy = _validate_session_policy(request, network)
        scientific_production_policy = _scientific_production_policy(policy, network)
        rate_origins, session_tst_ids = _rate_evidence_origins(network, initial_network)
        if retained_state is not None and session_tst_ids:
            retained_policy = retained_state.get("scientific_production_policy")
            if (not isinstance(retained_policy, Mapping)
                    or _canonical_hash(retained_policy) != _canonical_hash(scientific_production_policy)):
                # Do not replace a successful old report/bundle/state with a
                # refusal, nor create a new iteration that grants rate credit.
                return {
                    "schema_version": 1, "kind": "solution_verification_session",
                    "status": "invalid_or_incomplete_evidence",
                    "reason": (
                        "Retained session-produced TST rates have a missing or changed scientific production policy; "
                        "their temperature, thermochemistry-recipe and conformer admission must be revalidated "
                        "before prediction convergence can be claimed. Retained artifacts were preserved."
                    ),
                    "verified_bundle": None, "result_json": None,
                    "retained_result_json": str(result_path) if result_path.is_file() else None,
                    "state_json": str(state_path),
                    "session_produced_tst_reaction_ids": session_tst_ids,
                    "rate_evidence_origins": rate_origins,
                    "required_scientific_production_policy": scientific_production_policy,
                    "retained_scientific_production_policy": retained_policy,
                }
        # Active specialized-physics findings are state, not failed ORCA jobs.
        network = replace(network, candidates=tuple(
            replace(item, verification_state="method_escalation_required")
            if item.active_applicability_flags and item.verification_state == "pending" else item
            for item in network.candidates
        ))
        latest_assessment = None
        latest_convergence = None
        latest_allocation = None
        final_status = "not_converged"
        final_reason = "Session did not execute an iteration."
        selected_id = None

        while completed_iterations < policy["maximum_iterations"]:
            iteration = artifact_index
            iteration_dir = workdir / f"iteration-{iteration:03d}"
            iteration_dir.mkdir(parents=True, exist_ok=True)
            latest_assessment = evaluate_solution_t95_bounds(network)
            latest_convergence = assess_solution_session_convergence(
                network, latest_assessment,
                maximum_t95_ratio=policy["maximum_t95_ratio"],
                require_same_horizon_classification=policy["require_same_horizon_classification"],
            )
            latest_allocation = allocate_next_solution_verification(network, latest_assessment)
            _write_json(iteration_dir / "observable-bounds.json", latest_assessment)
            _write_json(iteration_dir / "convergence.json", latest_convergence)
            _write_json(iteration_dir / "allocation.json", latest_allocation)

            if latest_convergence["status"] == "prediction_converged":
                final_status = "prediction_converged"
                final_reason = "Verified rate bounds satisfy the requested t95 convergence contract."
                _write_json(iteration_dir / "decision.json", {
                    "status": final_status, "reason": final_reason, "selected_candidate_id": None,
                })
                completed_iterations += 1
                artifact_index += 1
                break
            selected_id = latest_allocation.get("selected_candidate_id")
            if not execute:
                final_status = "prepared"
                final_reason = (
                    "The next observable-conditioned verification decision is prepared; no ORCA jobs were run."
                )
                decision = {
                    "status": final_status, "reason": final_reason,
                    "selected_candidate_id": selected_id,
                }
                if selected_id:
                    candidate = next(item for item in network.candidates if item.candidate_id == selected_id)
                    decision["cache_key"] = _cache_key(candidate, network, policy["method_profile"])
                    _write_json(iteration_dir / "path-request.json", _path_request(network, candidate))
                _write_json(iteration_dir / "decision.json", decision)
                artifact_index += 1
                break
            if not selected_id:
                final_status, final_reason = _terminal_status(
                    network, latest_convergence, latest_allocation, used_paths, policy,
                )
                _write_json(iteration_dir / "decision.json", {
                    "status": final_status, "reason": final_reason, "selected_candidate_id": None,
                })
                completed_iterations += 1
                artifact_index += 1
                break
            candidate = next(item for item in network.candidates if item.candidate_id == selected_id)
            if (
                used_paths >= policy["maximum_orca_paths"]
                or used_core_hours + candidate.estimated_cost_core_hours > policy["maximum_estimated_core_hours"]
            ):
                final_status = "verification_budget_exhausted"
                final_reason = "The next influential path would exceed the immutable compute budget."
                _write_json(iteration_dir / "decision.json", {
                    "status": final_status, "reason": final_reason,
                    "selected_candidate_id": selected_id,
                })
                completed_iterations += 1
                artifact_index += 1
                break
            path_request = _path_request(network, candidate)
            _write_json(iteration_dir / "path-request.json", path_request)
            path_dir = iteration_dir / "path-verification"
            path_result = path_verifier(
                path_request, path_dir, execute=True, ncores=ncores,
                method_profile_name=policy["method_profile"],
                require_complete_network_coverage=False,
            )
            _write_json(iteration_dir / "path-result.json", path_result)
            used_paths += 1
            used_core_hours += candidate.estimated_cost_core_hours
            if path_result.get("verified_reaction"):
                reaction = VerifiedSolutionReaction.from_dict(path_result["verified_reaction"])
                network = network.replace_candidate_with_verified(selected_id, reaction)
                decision_status = "candidate_verified_and_network_repaired"
            else:
                state = (
                    "method_escalation_required"
                    if path_result.get("status") == "method_escalation_required"
                    else "path_unresolved"
                )
                network = network.replace_candidate_state(
                    selected_id, state, str(path_result.get("result_json") or iteration_dir / "path-result.json"),
                )
                decision_status = state
            _write_json(iteration_dir / "decision.json", {
                "status": decision_status, "selected_candidate_id": selected_id,
                "cache_key": _cache_key(candidate, network, policy["method_profile"]),
                "used_orca_paths": used_paths,
                "used_estimated_core_hours": used_core_hours,
            })
            completed_iterations += 1
            artifact_index += 1
            rate_origins, session_tst_ids = _rate_evidence_origins(network, initial_network)
            _write_json(state_path, {
                "schema_version": 1, "kind": "solution_verification_session_state",
                "request_sha256": request_hash, "candidate_network": network.as_dict(),
                "used_orca_paths": used_paths,
                "used_estimated_core_hours": used_core_hours,
                "completed_iterations": completed_iterations,
                "next_artifact_index": artifact_index,
                "scientific_production_policy": scientific_production_policy,
                "rate_evidence_origins": rate_origins,
                "session_produced_tst_reaction_ids": session_tst_ids,
            })
        else:
            final_status = "verification_budget_exhausted"
            final_reason = "The maximum verification-session iteration count was reached."

        # Re-evaluate after the final repair even when the iteration budget was
        # consumed by that repair.  Convergence checks do not allocate another
        # expensive calculation and must not be lost to an off-by-one policy.
        latest_assessment = evaluate_solution_t95_bounds(network)
        latest_convergence = assess_solution_session_convergence(
            network, latest_assessment, maximum_t95_ratio=policy["maximum_t95_ratio"],
            require_same_horizon_classification=policy["require_same_horizon_classification"],
        )
        latest_allocation = allocate_next_solution_verification(network, latest_assessment)
        if latest_convergence["status"] == "prediction_converged":
            final_status = "prediction_converged"
            final_reason = "Verified rate bounds satisfy the requested t95 convergence contract."
        elif execute and final_status in {"not_converged", "verification_budget_exhausted"}:
            terminal_status, terminal_reason = _terminal_status(
                network, latest_convergence, latest_allocation, used_paths, policy,
            )
            if terminal_status != "not_converged":
                final_status, final_reason = terminal_status, terminal_reason
        claim = {
            "schema_version": 1, "kind": "scientific_claim", "claim_type": "solution_t95",
            "status": "verified_within_declared_model" if final_status == "prediction_converged" else "unresolved",
            "claim": (
                "Solution t95 is bounded by verified first-order target-loss rates under the declared environment."
                if final_status == "prediction_converged" else
                "No defensible solution t95 is available from the retained verification session."
            ),
            "central_t95_seconds": (
                latest_assessment.get("central_t95_seconds")
                if final_status == "prediction_converged" else None
            ),
            "lower_t95_seconds": (
                latest_assessment.get("fast_t95_seconds")
                if final_status == "prediction_converged" else None
            ),
            "upper_t95_seconds": (
                latest_assessment.get("slow_t95_seconds")
                if final_status == "prediction_converged" else None
            ),
            "screening_assessment_retained_separately": final_status != "prediction_converged",
            "conditional_on_coverage_claim_id": network.coverage_claim_id,
            "convergence": latest_convergence,
        }
        claim["claim_id"] = "claim:solution-t95:" + _canonical_hash(claim).split(":", 1)[1]
        bundle_path = None
        propagation = None
        bundle_limitation = None
        if final_status == "prediction_converged":
            bundle_path, propagation, bundle_limitation = _verified_bundle(
                network, claim, workdir, _retained_path_claims(workdir),
            )
        rate_origins, session_tst_ids = _rate_evidence_origins(network, initial_network)
        state = {
            "schema_version": 1, "kind": "solution_verification_session_state",
            "request_sha256": request_hash, "candidate_network": network.as_dict(),
            "used_orca_paths": used_paths, "used_estimated_core_hours": used_core_hours,
            "completed_iterations": completed_iterations,
            "next_artifact_index": artifact_index,
            "scientific_production_policy": scientific_production_policy,
            "rate_evidence_origins": rate_origins,
            "session_produced_tst_reaction_ids": session_tst_ids,
        }
        _write_json(state_path, state)
        result = {
            "schema_version": 1, "kind": "solution_verification_session",
            "status": final_status, "reason": final_reason,
            "request_sha256": request_hash, "policy": policy,
            "scientific_production_policy": scientific_production_policy,
            "rate_evidence_origins": rate_origins,
            "session_produced_tst_reaction_ids": session_tst_ids,
            "used_orca_paths": used_paths,
            "used_estimated_core_hours": used_core_hours,
            "completed_iterations": completed_iterations,
            "selected_candidate_id": selected_id,
            "candidate_network": network.as_dict(),
            "observable_assessment": latest_assessment,
            "convergence": latest_convergence,
            "allocation": latest_allocation,
            "scientific_claim": claim, "verified_bundle": bundle_path,
            "bundle_limitation": bundle_limitation,
            "propagation_result": propagation,
            "result_json": str(result_path), "state_json": str(state_path),
        }
        _write_json(result_path, result)
        return result
    except Exception as error:
        result = {
            "schema_version": 1, "kind": "solution_verification_session",
            "status": "invalid_or_incomplete_evidence",
            "reason": f"{type(error).__name__}: {error}",
            "verified_bundle": None, "result_json": str(result_path),
        }
        _write_json(result_path, result)
        return result


def load_and_run_solution_verification_session(
    request_path: Path, workdir: Path, **kwargs: Any,
) -> dict[str, Any]:
    request = json.loads(Path(request_path).read_text())
    return run_solution_verification_session(request, workdir, **kwargs)


__all__ = ["load_and_run_solution_verification_session", "run_solution_verification_session"]
