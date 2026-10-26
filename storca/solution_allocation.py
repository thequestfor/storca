"""Auditable next-reaction allocation for solution t95 verification."""

from __future__ import annotations

import math
from typing import Any

from .solution_candidates import SolutionCandidateNetwork


def allocate_next_solution_verification(
    network: SolutionCandidateNetwork,
    assessment: dict[str, Any],
) -> dict[str, Any]:
    contribution = {
        item["reaction_id"]: item for item in assessment.get("contributions") or []
        if item.get("record_kind") == "candidate_reaction"
    }
    total_upper = assessment.get("total_upper_target_loss_rate_s")
    crossing_rate = -math.log(network.environment.retention_fraction) / network.environment.target_duration_seconds
    entries = []
    for candidate in network.candidates:
        item = contribution.get(candidate.candidate_id)
        affects_t95 = candidate.reactant_id == network.target_species_id
        active_flags = candidate.active_applicability_flags
        eligible = candidate.verification_state == "pending" and affects_t95 and not active_flags
        if not affects_t95:
            priority_class = 6
            reason = "irrelevant_to_t95_under_monotonic_target_loss"
        elif candidate.verification_state != "pending":
            priority_class = 7
            reason = f"verification_state_{candidate.verification_state}"
        elif active_flags:
            priority_class = 0
            reason = "specialized_method_escalation_required"
            eligible = False
        elif candidate.approximate_rate.status in {"unbounded", "unavailable"} or candidate.approximate_rate.upper_rate_s is None:
            priority_class = 1
            reason = "influential_rate_unbounded"
        else:
            upper = float(candidate.approximate_rate.upper_rate_s)
            upper_without = None if total_upper is None else max(0.0, float(total_upper) - upper)
            if total_upper is None or (float(total_upper) >= crossing_rate and upper_without < crossing_rate):
                priority_class = 2
                reason = "candidate_can_change_horizon_classification"
            elif candidate.approximate_rate.use_scope != "prediction_bound":
                priority_class = 3
                reason = "screening_interval_not_admissible_for_prediction"
            else:
                priority_class = 4
                reason = "bounded_residual_t95_uncertainty"
        lower = float(candidate.approximate_rate.lower_rate_s or 0.0)
        upper = candidate.approximate_rate.upper_rate_s
        span = None if upper is None else max(0.0, float(upper) - lower)
        influence = (
            item.get("maximum_upper_rate_fraction") if item else 0.0
        )
        if influence is None:
            influence = 1.0 if affects_t95 else 0.0
        information_term = (
            50.0 if span is None else math.log1p(span / max(crossing_rate, 1e-300))
        )
        source_factor = 1.0 + 0.1 * min(5, len(candidate.discovery_sources))
        score = float(influence) * information_term * source_factor / candidate.estimated_cost_core_hours
        entries.append({
            "candidate_id": candidate.candidate_id,
            "eligible_for_current_verifier": eligible,
            "affects_t95": affects_t95,
            "priority_class": priority_class, "priority_reason": reason,
            "observable_influence": influence,
            "rate_interval_width_s": span,
            "information_term": information_term,
            "discovery_source_factor": source_factor,
            "estimated_cost_core_hours": candidate.estimated_cost_core_hours,
            "score": score, "active_applicability_flags": list(active_flags),
            "score_formula": "observable_influence * information_term * discovery_source_factor / estimated_cost_core_hours",
        })
    eligible = [item for item in entries if item["eligible_for_current_verifier"]]
    selected = min(
        eligible,
        key=lambda item: (item["priority_class"], -item["score"], item["candidate_id"]),
        default=None,
    )
    return {
        "schema_version": 1, "kind": "solution_verification_allocation",
        "status": "candidate_selected" if selected else "no_eligible_candidate",
        "selected_candidate_id": selected["candidate_id"] if selected else None,
        "selected_reason": selected["priority_reason"] if selected else None,
        "entries": sorted(entries, key=lambda item: (item["priority_class"], -item["score"], item["candidate_id"])),
        "policy": {
            "observable": "t95", "first_apply_priority_classes": True,
            "screening_rates_are_allocation_only": True,
            "downstream_reactions_receive_zero_t95_priority": True,
        },
    }


__all__ = ["allocate_next_solution_verification"]
