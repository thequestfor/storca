"""Observable bounds for the restricted monotonic first-order solution lane."""

from __future__ import annotations

import math
from typing import Any

from .solution_candidates import SolutionCandidateNetwork


def evaluate_solution_t95_bounds(network: SolutionCandidateNetwork) -> dict[str, Any]:
    """Evaluate exact direct-loss t95 bounds without promoting screening rates.

    Under the lane's no-target-reformation contract, downstream reactions do
    not affect target disappearance.  The sum of direct target-consuming
    first-order rates therefore determines t95 exactly for each rate vector.
    """
    target = network.target_species_id
    contributions = []
    for reaction in network.verified_reactions:
        if reaction.reactant_id != target:
            continue
        lower = reaction.rate_constant_s * float(reaction.lower_rate_multiplier or 1.0)
        upper = reaction.rate_constant_s * float(reaction.upper_rate_multiplier or 1.0)
        contributions.append({
            "reaction_id": reaction.reaction_id, "record_kind": "verified_reaction",
            "verification_state": "verified", "central_rate_s": reaction.rate_constant_s,
            "lower_rate_s": lower, "upper_rate_s": upper,
            "prediction_bound_admissible": True, "source": reaction.rate_claim_id,
        })
    for candidate in network.candidates:
        if candidate.reactant_id != target:
            continue
        rate = candidate.approximate_rate
        contributions.append({
            "reaction_id": candidate.candidate_id, "record_kind": "candidate_reaction",
            "verification_state": candidate.verification_state,
            "central_rate_s": rate.central_rate_s,
            "lower_rate_s": rate.lower_rate_s,
            "upper_rate_s": rate.upper_rate_s,
            "prediction_bound_admissible": rate.use_scope == "prediction_bound",
            "source": rate.source, "rate_status": rate.status,
        })
    finite_central = [float(item["central_rate_s"]) for item in contributions if item["central_rate_s"] is not None]
    total_central = sum(finite_central)
    total_lower = sum(float(item["lower_rate_s"] or 0.0) for item in contributions)
    upper_unbounded = any(item["upper_rate_s"] is None for item in contributions)
    total_upper = None if upper_unbounded else sum(float(item["upper_rate_s"]) for item in contributions)
    crossing = -math.log(network.environment.retention_fraction)

    def t95(rate):
        if rate is None:
            return 0.0
        return crossing / rate if rate > 0.0 else None

    central_t95 = t95(total_central)
    fast_t95 = t95(total_upper)
    slow_t95 = t95(total_lower)
    horizon = network.environment.target_duration_seconds
    if slow_t95 is not None and slow_t95 <= horizon:
        horizon_status = "crossing_within_horizon_for_all_bounds"
    elif fast_t95 is None or fast_t95 > horizon:
        horizon_status = "no_crossing_within_horizon_for_any_bound"
    else:
        horizon_status = "crossing_depends_on_rate_bounds"
    interval_ratio = (
        slow_t95 / fast_t95
        if slow_t95 is not None and fast_t95 is not None and fast_t95 > 0.0
        else None
    )
    for item in contributions:
        central = item["central_rate_s"]
        item["absolute_log_t95_sensitivity"] = (
            float(central) / total_central if central is not None and total_central > 0.0 else None
        )
        upper = item["upper_rate_s"]
        item["maximum_upper_rate_fraction"] = (
            float(upper) / total_upper
            if upper is not None and total_upper is not None and total_upper > 0.0 else None
        )
        item["affects_t95"] = True
    downstream = [
        {
            "reaction_id": item.reaction_id, "record_kind": "verified_reaction",
            "affects_t95": False, "reason": "reactant_is_not_target_under_monotonic_loss_contract",
        }
        for item in network.verified_reactions if item.reactant_id != target
    ] + [
        {
            "reaction_id": item.candidate_id, "record_kind": "candidate_reaction",
            "affects_t95": False, "reason": "reactant_is_not_target_under_monotonic_loss_contract",
        }
        for item in network.candidates if item.reactant_id != target
    ]
    pending_target = [
        item.candidate_id for item in network.candidates
        if item.reactant_id == target and item.verification_state == "pending"
    ]
    unresolved_target = [
        item.candidate_id for item in network.candidates
        if item.reactant_id == target and item.verification_state in {
            "path_unresolved", "method_escalation_required",
        }
    ]
    inadmissible_bounds = [
        item["reaction_id"] for item in contributions
        if item["record_kind"] == "candidate_reaction" and not item["prediction_bound_admissible"]
    ]
    return {
        "schema_version": 1, "kind": "solution_t95_bound_assessment",
        "status": "completed", "observable": "t95",
        "retention_fraction": network.environment.retention_fraction,
        "horizon_seconds": horizon, "central_t95_seconds": central_t95,
        "fast_t95_seconds": fast_t95, "slow_t95_seconds": slow_t95,
        "t95_interval_ratio": interval_ratio, "horizon_status": horizon_status,
        "total_central_target_loss_rate_s": total_central,
        "total_lower_target_loss_rate_s": total_lower,
        "total_upper_target_loss_rate_s": total_upper,
        "upper_rate_unbounded": upper_unbounded,
        "contributions": contributions, "downstream_irrelevant_reactions": downstream,
        "pending_target_candidate_ids": pending_target,
        "unresolved_target_candidate_ids": unresolved_target,
        "candidate_bounds_not_admissible_for_prediction": inadmissible_bounds,
        "interpretation": (
            "Candidate intervals are used to allocate verification. Only verified reaction bounds "
            "are admissible for a final solution-lane prediction."
        ),
    }


def assess_solution_session_convergence(
    network: SolutionCandidateNetwork,
    assessment: dict[str, Any],
    *,
    maximum_t95_ratio: float,
    require_same_horizon_classification: bool = True,
) -> dict[str, Any]:
    if not math.isfinite(maximum_t95_ratio) or maximum_t95_ratio < 1.0:
        raise ValueError("maximum_t95_ratio must be finite and at least one")
    blockers = []
    if network.coverage_status != "verified_complete_within_declared_unimolecular_scope":
        blockers.append("network_coverage_unresolved")
    if not network.coverage_claim_id or not network.discovery_sources:
        blockers.append("network_coverage_provenance_missing")
    if network.unresolved_discovery_ids:
        blockers.append("unresolved_discovery_candidates")
    if assessment["pending_target_candidate_ids"]:
        blockers.append("target_rate_verification_pending")
    if assessment["unresolved_target_candidate_ids"]:
        blockers.append("target_path_or_method_unresolved")
    if assessment["candidate_bounds_not_admissible_for_prediction"]:
        blockers.append("screening_bounds_not_prediction_evidence")
    if require_same_horizon_classification and assessment["horizon_status"] == "crossing_depends_on_rate_bounds":
        blockers.append("horizon_classification_not_robust")
    ratio = assessment["t95_interval_ratio"]
    if ratio is None or ratio > maximum_t95_ratio:
        blockers.append("t95_interval_too_wide")
    if not network.verified_reactions:
        blockers.append("no_verified_solution_rates")
    return {
        "schema_version": 1, "kind": "solution_session_convergence_assessment",
        "status": "prediction_converged" if not blockers else "not_converged",
        "blockers": sorted(set(blockers)),
        "maximum_t95_ratio": maximum_t95_ratio,
        "require_same_horizon_classification": require_same_horizon_classification,
        "observed_t95_ratio": ratio,
        "horizon_status": assessment["horizon_status"],
    }


__all__ = ["assess_solution_session_convergence", "evaluate_solution_t95_bounds"]
