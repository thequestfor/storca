"""Flux-guided, family-agnostic kinetic verification policy.

This module deliberately separates *selection* from quantum-chemistry path
generation.  It consumes retained solver/Cantera evidence, chooses one
controlling target-loss route at a time, and says when a repaired lifetime is
insensitive to the routes that remain unresolved.  It never estimates a flux,
rate bound, or lifetime that is absent from an artifact.
"""

from __future__ import annotations

import math
import re
from typing import Iterable

from .generated_kinetics import validate_rate_replacement_evidence
from .rmg_evidence import normalize_rmg_label


DEFAULT_FLUX_COVERAGE_TARGET = 0.95
DEFAULT_T95_RELATIVE_TOLERANCE = 0.20
DEFAULT_MAX_FLUX_CLOSURE_ERROR = 0.20
# A closure failure is only waived when all relevant amounts are smaller than
# this fraction of the initial target inventory.  This is a numerical
# classification, never a kinetic/lifetime criterion.
DEFAULT_NUMERICAL_MATERIALITY_RELATIVE_FLOOR = 1e-10
DEFAULT_NUMERICAL_MATERIALITY_ABSOLUTE_FLOOR_KMOL = 1e-18
# The retained RMG graph can contain hundreds of mutually dependent radical
# routes.  Route selection is a bounded triage operation, never an exhaustive
# graph-enumeration job.
MAX_DEPENDENCY_SEARCH_STATES = 256
MAX_PRODUCERS_PER_INTERMEDIATE = 12
_ARROW = re.compile(r"<=>|=>|=")
_AVOGADRO = 6.02214076e23


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _clean_equation(equation: str) -> str:
    """Discard retained comments/rate columns while preserving the equation."""
    lines = [line.strip() for line in str(equation).splitlines() if line.strip() and not line.lstrip().startswith("!")]
    text = lines[-1] if lines else ""
    match = re.search(r"(.+?(?:<=>|=>|=).+?)(?:\s{2,}[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+)?$", text)
    return (match.group(1) if match else text).strip()


def _side_signature(side: str, *, normalized: bool) -> tuple[tuple[str, float], ...]:
    counts: dict[str, float] = {}
    for token in side.split("+"):
        token = token.strip()
        if not token:
            continue
        match = re.match(r"^(?:(\d+(?:\.\d+)?)\s+)?(.+?)$", token)
        coefficient = float(match.group(1) or 1.0)
        label = match.group(2).strip()
        if normalized:
            label = normalize_rmg_label(label)
        counts[label] = counts.get(label, 0.0) + coefficient
    return tuple(sorted(counts.items()))


def reaction_signature(equation: str, *, normalized: bool = False) -> tuple | None:
    """Return an orientation-independent stoichiometric reaction signature."""
    cleaned = _clean_equation(equation)
    parts = _ARROW.split(cleaned, maxsplit=1)
    if len(parts) != 2:
        return None
    sides = sorted((_side_signature(parts[0], normalized=normalized),
                    _side_signature(parts[1], normalized=normalized)))
    return tuple(sides)


def _oriented_signature(equation: str, *, normalized: bool = False) -> tuple | None:
    parts = _ARROW.split(_clean_equation(equation), maxsplit=1)
    if len(parts) != 2:
        return None
    return (_side_signature(parts[0], normalized=normalized),
            _side_signature(parts[1], normalized=normalized))


def _rate_to_si(value: float, units: str, molecularity: int) -> float | None:
    compact = units.lower().replace("^", "").replace("(", "").replace(")", "")
    compact = re.sub(r"\s+", "/", compact.replace("*", "/")).strip("/")
    compact = re.sub(r"/+", "/", compact)
    aliases = {
        (1, "s-1"): 1.0,
        (2, "m3/mol/s"): 1.0,
        (2, "m3mol-1s-1"): 1.0,
        (2, "cm3/mol/s"): 1e-6,
        (2, "cm3mol-1s-1"): 1e-6,
        (2, "cm3/molecule/s"): 1e-6 * _AVOGADRO,
        (2, "cm3molecule-1s-1"): 1e-6 * _AVOGADRO,
        (3, "m6/mol2/s"): 1.0,
        (3, "m6mol-2s-1"): 1.0,
        (3, "cm6/mol2/s"): 1e-12,
        (3, "cm6mol-2s-1"): 1e-12,
        (3, "cm6/molecule2/s"): 1e-12 * _AVOGADRO ** 2,
        (3, "cm6molecule-2s-1"): 1e-12 * _AVOGADRO ** 2,
    }
    factor = aliases.get((molecularity, compact))
    if factor is None:
        factor = aliases.get((molecularity, compact.replace("/", "")))
    return value * factor if factor is not None else None


def verify_applied_replacement(
    repaired_rmg_evidence: dict,
    replacement: dict,
    condition_contract: dict,
    *,
    relative_tolerance: float = 1e-4,
) -> dict:
    """Match a retained mechanism reaction to its verified replacement rate."""
    blockers = []
    expected_equation = replacement.get("reaction_equation") or replacement.get("replaces_reaction") or ""
    expected_signature = reaction_signature(expected_equation, normalized=True)
    expected_oriented = _oriented_signature(expected_equation, normalized=True)
    inspection = repaired_rmg_evidence.get("mechanism_inspection") or {}
    candidates = [item for item in inspection.get("reactions") or []
                  if reaction_signature(item.get("equation", ""), normalized=True) == expected_signature]
    if not expected_signature:
        blockers.append("replacement_reaction_signature_missing")
    if len(candidates) != 1:
        blockers.append("applied_reaction_missing_or_ambiguous")
        candidate = {}
    else:
        candidate = candidates[0]
    applied_oriented = (_oriented_signature(candidate.get("equation", ""), normalized=True)
                        if candidate else None)
    reversed_orientation = bool(
        candidate and expected_oriented is not None and applied_oriented == tuple(reversed(expected_oriented))
    )
    if candidate and applied_oriented != expected_oriented and not reversed_orientation:
        blockers.append("applied_reaction_direction_mismatch")
    molecularity = int(replacement.get("molecularity") or 0)
    applied_molecularity = int(candidate.get(
        "reverse_molecularity" if reversed_orientation else "molecularity",
    ) or 0) if candidate else 0
    if candidate and applied_molecularity != molecularity:
        blockers.append("applied_reaction_molecularity_mismatch")
    expected = replacement.get("condition_rate") or replacement.get("rate_constant") or {}
    if isinstance(expected, (int, float)):
        expected = {"value": expected, "units": replacement.get("rate_units")}
    expected_value = _finite(expected.get("value"))
    expected_si = (_rate_to_si(expected_value, str(expected.get("units") or ""), molecularity)
                   if expected_value is not None else None)
    applied_value = _finite(candidate.get(
        "evaluated_reverse_rate_coefficient_si" if reversed_orientation
        else "evaluated_forward_rate_coefficient_si"
    )) if candidate else None
    evaluated_at = candidate.get("evaluated_at") or {}
    if expected_si is None:
        blockers.append("replacement_rate_units_not_convertible_to_si")
    if applied_value is None:
        blockers.append("applied_condition_rate_missing")
    if (_finite(evaluated_at.get("temperature_K")) != _finite(condition_contract.get("temperature_K")) or
            _finite(evaluated_at.get("pressure_bar")) != _finite(condition_contract.get("pressure_bar"))):
        blockers.append("applied_rate_condition_mismatch")
    relative_error = None
    if expected_si is not None and applied_value is not None:
        scale = max(abs(expected_si), abs(applied_value), 1e-300)
        relative_error = abs(expected_si - applied_value) / scale
        if relative_error > relative_tolerance:
            blockers.append("applied_rate_does_not_match_replacement")
    provenance_values = [str(candidate.get("source_library") or ""), str(candidate.get("kinetics_comment") or "")]
    expected_source = str(replacement.get("source_library") or replacement.get("library_name") or "")
    provenance_available = any(provenance_values)
    if expected_source and provenance_available and not any(expected_source in value for value in provenance_values):
        blockers.append("applied_rate_library_provenance_mismatch")
    return {
        "status": "passed" if not blockers else "failed",
        "route_id": replacement.get("route_id"),
        "matched_reaction_index": candidate.get("reaction_index") if candidate else None,
        "applied_direction": "reverse" if reversed_orientation else "forward",
        "expected_rate_si": expected_si,
        "applied_rate_si": applied_value,
        "relative_error": relative_error,
        "relative_tolerance": relative_tolerance,
        "library_provenance_checked": bool(expected_source and provenance_available),
        "blocking_reasons": blockers,
    }


def _route_equations(route: dict) -> list[str]:
    values = [route.get("reaction_equation"), route.get("printed_reaction_equation")]
    endpoints = route.get("resolved_endpoints") or {}
    reactants = endpoints.get("reactants") or []
    products = endpoints.get("products") or []
    if reactants and products:
        values.append(
            " + ".join(str(item.get("label")) for item in reactants) + " <=> " +
            " + ".join(str(item.get("label")) for item in products)
        )
    return [str(value) for value in values if value]


def _invalid_screening_kinetics(rmg_evidence: dict) -> bool:
    """True when retained RMG kinetics cannot safely rank loss routes."""
    kinetics = rmg_evidence.get("kinetics_validation") or {}
    return (
        str(kinetics.get("status") or "") in {"kinetics_unreliable", "failed"}
        or int(kinetics.get("violation_count") or 0) > 0
    )


def _direct_unimolecular_rmg_candidates(routes: list[dict], target_label: str) -> list[tuple[int, dict]]:
    """Select *only RMG-generated* direct target decomposition hypotheses.

    This deliberately does not enumerate bonds, rearrangements, or products.
    It is the failover queue used when an RMG rate which would otherwise rank
    the network is physically invalid (for example, above a collision bound).
    """
    target = normalize_rmg_label(str(target_label))
    selected: list[tuple[int, dict]] = []
    for index, raw_route in enumerate(routes):
        endpoints = raw_route.get("resolved_endpoints") or {}
        reactants = list(endpoints.get("reactants") or [])
        products = list(endpoints.get("products") or [])
        if len(reactants) != 1 or not products:
            continue
        label = normalize_rmg_label(str((reactants[0] or {}).get("label") or ""))
        if label != target:
            continue
        # A direct unimolecular route has no external initiation requirement.
        # Endpoint/surface issues are retained for ORCA to resolve/defer; they
        # are not silently discarded by selection.
        if raw_route.get("initiation_status") not in {
            None, "directly_present", "thermally_reachable",
        }:
            continue
        selected.append((index, raw_route))
    return selected


def _direct_partial_enumeration_candidates(routes: list[dict]) -> list[tuple[int, dict]]:
    """Keep only direct, already-enumerated routes after an RMG solver stop.

    Partial Chemkin output is never kinetic evidence.  It can nevertheless
    provide a bounded ORCA queue if a route consumes the target and every
    non-target reactant was declared in the condition contract.  Generated
    radicals/intermediates are intentionally excluded until their producer is
    validated in a later iteration.
    """
    selected: list[tuple[int, dict]] = []
    for index, raw_route in enumerate(routes):
        endpoints = raw_route.get("resolved_endpoints") or {}
        if not endpoints.get("reactants") or not endpoints.get("products"):
            continue
        if (raw_route.get("atom_balance") or {}).get("status") != "passed":
            continue
        if raw_route.get("initiation_status") not in {None, "directly_present"}:
            continue
        co_reactants = list(raw_route.get("co_reactants") or [])
        if any(item.get("origin") not in {None, "initial_scenario_species"} for item in co_reactants):
            continue
        if float(raw_route.get("target_stoichiometry") or 0.0) <= 0.0:
            continue
        selected.append((index, raw_route))
    return selected


def _partial_enumeration_plan(
    rmg_evidence: dict,
    condition_contract: dict,
    *,
    deferred_ids: set[str],
    verified_ids: set[str],
) -> dict | None:
    """Build an ORCA-only queue from a failed RMG propagation artifact."""
    if rmg_evidence.get("enumeration_status") != "partial_enumeration":
        return None
    candidates = _direct_partial_enumeration_candidates(list(rmg_evidence.get("candidate_routes") or []))
    selectable = [
        item for item in candidates
        if str(item[1].get("route_id") or f"rmg-route-{item[0] + 1}") not in deferred_ids
        and str(item[1].get("route_id") or f"rmg-route-{item[0] + 1}") not in verified_ids
    ]
    if not selectable:
        return {
            "schema_version": 1,
            "kind": "flux_ranked_route_verification_plan",
            "status": "partial_enumeration_no_direct_route",
            "selected_route_id": None,
            "ranked_routes": [],
            "blocking_reasons": ["rmg_propagation_incomplete"],
            "rule": "Partial RMG enumeration contained no direct, declared-reactant target-loss route eligible for ORCA.",
        }
    selected_index, selected = selectable[0]
    route_id = str(selected.get("route_id") or f"rmg-route-{selected_index + 1}")
    return {
        "schema_version": 1,
        "kind": "flux_ranked_route_verification_plan",
        "status": "verification_required",
        "ranking_basis": "direct_reachable_route_from_partial_rmg_enumeration",
        "selection_reason": "RMG propagation failed, but this direct atom-resolved route was already enumerated and must be checked by ORCA rather than discarded.",
        "partial_enumeration": True,
        "selected_route_id": route_id,
        "selected_route_index": selected_index,
        "selected_route_source": "candidate_routes",
        "ranked_routes": [{
            "rank": rank,
            "route_index": index,
            "route_source": "candidate_routes",
            "route_id": str(route.get("route_id") or f"rmg-route-{index + 1}"),
            "reaction_equation": route.get("reaction_equation") or route.get("printed_reaction_equation"),
            "initiation_status": route.get("initiation_status") or "directly_present",
            "verification_queue_class": "direct_partial_rmg_enumeration_candidate",
            "verified_replacement": False,
        } for rank, (index, route) in enumerate(selectable, 1)],
        "target_loss_commitment": {"status": "unavailable_partial_rmg_enumeration"},
        "blocking_reasons": [],
        "requires_repropagation_after_each_verification": True,
        "rule": "No RMG rate, flux, or t95 is used from partial enumeration; ORCA must validate every queued route.",
    }


def _numerically_negligible_closure(flux: dict) -> bool:
    """Whether an otherwise poor relative closure is immaterial numerical dust.

    New bridge artifacts provide an inventory and an absolute discrepancy.
    Legacy artifacts deliberately fail this test so they retain the old
    fail-closed behavior.
    """
    initial = _finite(flux.get("initial_target_amount_kmol"))
    absolute = _finite(flux.get("numerical_closure_absolute_error_kmol"))
    profile_loss = _finite(flux.get("amount_profile_target_loss_kmol"))
    integrated_loss = _finite(flux.get("integrated_net_target_loss_kmol"))
    if None in {initial, absolute, profile_loss, integrated_loss} or initial is None or initial <= 0.0:
        return False
    floor = max(
        DEFAULT_NUMERICAL_MATERIALITY_ABSOLUTE_FLOOR_KMOL,
        initial * DEFAULT_NUMERICAL_MATERIALITY_RELATIVE_FLOOR,
    )
    return max(abs(absolute), abs(profile_loss), abs(integrated_loss)) <= floor


def _solver_gates(rmg_evidence: dict, *, max_closure_error: float) -> list[str]:
    reasons = []
    if rmg_evidence.get("status") != "completed":
        reasons.append("rmg_model_not_completed")
    execution = rmg_evidence.get("execution_assessment") or {}
    if execution and execution.get("status") != "completed":
        reasons.append("rmg_execution_not_completed")
    if not (rmg_evidence.get("time_coverage") or {}).get("complete"):
        reasons.append("requested_time_coverage_incomplete")
    propagation = rmg_evidence.get("independent_cantera_propagation") or {}
    if propagation.get("status") != "completed":
        reasons.append("independent_propagation_missing_or_incomplete")
    if not (propagation.get("coverage") or {}).get("complete"):
        reasons.append("independent_propagation_time_coverage_incomplete")
    flux = propagation.get("reaction_flux_attribution") or {}
    if flux.get("status") != "completed":
        reasons.append("quantitative_reaction_flux_missing_or_invalid")
    closure = _finite(flux.get("numerical_closure_relative_error"))
    if closure is None:
        reasons.append("reaction_flux_numerical_closure_missing")
    elif closure > max_closure_error and not _numerically_negligible_closure(flux):
        reasons.append("reaction_flux_numerical_closure_failed")
    return reasons


def _flux_matches(routes: list[dict], reactions: list[dict]) -> tuple[dict[int, list[int]], list[int]]:
    """Map flux records to routes, preferring exact indexed labels.

    Normalized-label fallback is admitted only when it identifies exactly one
    route.  This prevents two RMG isomers sharing a formula-like compatibility
    label from silently receiving one another's flux.
    """
    exact: dict[tuple, list[int]] = {}
    normalized: dict[tuple, list[int]] = {}
    for route_index, route in enumerate(routes):
        for equation in _route_equations(route):
            signature = reaction_signature(equation)
            if signature is not None:
                exact.setdefault(signature, []).append(route_index)
            signature = reaction_signature(equation, normalized=True)
            if signature is not None:
                normalized.setdefault(signature, []).append(route_index)
    matches: dict[int, list[int]] = {}
    unmatched = []
    for flux_index, reaction in enumerate(reactions):
        equation = reaction.get("reaction_equation", "")
        candidates = sorted(set(exact.get(reaction_signature(equation), [])))
        if not candidates:
            candidates = sorted(set(normalized.get(reaction_signature(equation, normalized=True), [])))
            if len(candidates) != 1:
                candidates = []
        if candidates:
            matches[flux_index] = candidates
        else:
            unmatched.append(flux_index)
    return matches, unmatched


def _sensitivity_by_signature(propagation: dict) -> dict[tuple, dict]:
    evidence = propagation.get("reaction_sensitivity_attribution") or {}
    if evidence.get("status") != "completed":
        return {}
    result = {}
    for item in evidence.get("reactions") or []:
        signature = reaction_signature(item.get("reaction_equation", ""), normalized=True)
        sensitivity = _finite(item.get("absolute_log_t95_sensitivity"))
        if signature is not None and sensitivity is not None:
            result[signature] = item
    return result


def assess_target_loss_commitment(flux: dict, target_label: str) -> dict:
    """Attribute horizon-net parent loss after cancelling return pathways.

    Every retained reaction is oriented by its integrated net extent. Positive
    target-loss reactions create one or more non-parent species; negative
    target-loss reactions consume non-parent species while reforming the
    parent. A recovery reaction is credited only to loss reactions whose
    product-side species can reach that recovery precursor through retained,
    nonzero net flux. Unattributable recovery fails closed instead of being
    used to reassure.
    """
    target = normalize_rmg_label(str(target_label))
    records = [dict(item) for item in flux.get("reactions") or []]
    required = ("integrated_net_extent_kmol", "integrated_net_target_loss_kmol")
    if not records or any(_finite(item.get(key)) is None for item in records for key in required):
        return {
            "schema_version": 1,
            "kind": "target_loss_commitment_attribution",
            "status": "unavailable_legacy_flux",
            "reason": "Directional net extents and per-reaction net target loss are required.",
            "reaction_commitment": [],
        }

    directed: list[dict] = []
    adjacency: dict[str, set[str]] = {}
    for position, item in enumerate(records):
        oriented = _oriented_signature(item.get("reaction_equation", ""), normalized=True)
        if oriented is None:
            return {
                "schema_version": 1,
                "kind": "target_loss_commitment_attribution",
                "status": "unresolved_reaction_equation",
                "reason": f"Could not orient retained reaction {item.get('reaction_index', position)}.",
                "reaction_commitment": [],
            }
        extent = float(item["integrated_net_extent_kmol"])
        inputs, outputs = oriented if extent >= 0.0 else tuple(reversed(oriented))
        input_labels = {label for label, coefficient in inputs if coefficient > 0.0 and label != target}
        output_labels = {label for label, coefficient in outputs if coefficient > 0.0 and label != target}
        if abs(extent) > 0.0:
            for source in input_labels:
                adjacency.setdefault(source, set()).update(output_labels)
        directed.append({
            "position": position,
            "reaction_index": item.get("reaction_index", position),
            "reaction_equation": item.get("reaction_equation"),
            "input_species": sorted(input_labels),
            "output_species": sorted(output_labels),
            "net_target_loss_kmol": float(item["integrated_net_target_loss_kmol"]),
        })

    def reachable(starts: set[str], goals: set[str]) -> bool:
        if starts & goals:
            return True
        pending = list(starts)
        visited = set(starts)
        while pending:
            node = pending.pop()
            for neighbor in adjacency.get(node, ()):
                if neighbor in goals:
                    return True
                if neighbor not in visited:
                    visited.add(neighbor)
                    pending.append(neighbor)
        return False

    losses = [item for item in directed if item["net_target_loss_kmol"] > 0.0]
    recoveries = [item for item in directed if item["net_target_loss_kmol"] < 0.0]
    remaining = {item["position"]: item["net_target_loss_kmol"] for item in losses}
    credits = {item["position"]: 0.0 for item in losses}
    recovery_records = []
    unallocated_recovery = 0.0
    for recovery in sorted(recoveries, key=lambda item: abs(item["net_target_loss_kmol"]), reverse=True):
        amount = -recovery["net_target_loss_kmol"]
        goals = set(recovery["input_species"])
        eligible = [
            loss for loss in losses
            if remaining[loss["position"]] > 0.0
            and goals
            and all(reachable(set(loss["output_species"]), {goal}) for goal in goals)
        ]
        capacity = sum(remaining[item["position"]] for item in eligible)
        allocated = min(amount, capacity)
        if allocated > 0.0 and capacity > 0.0:
            for loss in eligible:
                position = loss["position"]
                share = allocated * remaining[position] / capacity
                credits[position] += share
                remaining[position] = max(0.0, remaining[position] - share)
        unallocated = max(0.0, amount - allocated)
        unallocated_recovery += unallocated
        recovery_records.append({
            "reaction_index": recovery["reaction_index"],
            "reaction_equation": recovery["reaction_equation"],
            "parent_reformation_kmol": amount,
            "allocated_parent_return_credit_kmol": allocated,
            "unallocated_parent_return_kmol": unallocated,
            "credited_loss_reaction_indices": [item["reaction_index"] for item in eligible],
        })

    reaction_commitment = []
    for loss in losses:
        position = loss["position"]
        gross_net = loss["net_target_loss_kmol"]
        committed = remaining[position]
        reaction_commitment.append({
            "reaction_index": loss["reaction_index"],
            "reaction_equation": loss["reaction_equation"],
            "positive_net_target_loss_kmol": gross_net,
            "allocated_parent_return_credit_kmol": credits[position],
            "committed_target_loss_kmol": committed,
            "commitment_fraction": committed / gross_net if gross_net > 0.0 else None,
            "target_derived_product_species": loss["output_species"],
        })
    positive_loss = sum(item["net_target_loss_kmol"] for item in losses)
    total_recovery = sum(-item["net_target_loss_kmol"] for item in recoveries)
    committed_total = sum(item["committed_target_loss_kmol"] for item in reaction_commitment)
    retained_net = _finite(flux.get("integrated_net_target_loss_kmol"))
    expected_net = max(0.0, positive_loss - total_recovery)
    scale = max(positive_loss, total_recovery, abs(retained_net or 0.0), 1e-300)
    attribution_error = abs(committed_total - expected_net) / scale
    retained_error = (
        abs((retained_net or 0.0) - (positive_loss - total_recovery)) / scale
        if retained_net is not None else None
    )
    complete = unallocated_recovery / scale <= 1e-8 and attribution_error <= 1e-8
    status = (
        "no_committed_target_loss" if complete and committed_total <= scale * 1e-15
        else "completed" if complete
        else "unresolved_parent_return_attribution"
    )
    return {
        "schema_version": 1,
        "kind": "target_loss_commitment_attribution",
        "status": status,
        "basis": "integrated_net_flux_with_flux_supported_parent_return_reachability",
        "target_label": target,
        "positive_net_target_loss_kmol": positive_loss,
        "integrated_parent_reformation_kmol": total_recovery,
        "committed_target_loss_kmol": committed_total,
        "unallocated_parent_reformation_kmol": unallocated_recovery,
        "commitment_fraction": committed_total / positive_loss if positive_loss > 0.0 else 0.0,
        "attribution_closure_relative_error": attribution_error,
        "retained_net_loss_closure_relative_error": retained_error,
        "reaction_commitment": reaction_commitment,
        "parent_reformation_credit": recovery_records,
        "interpretation": (
            "Committed loss is parent consumption not cancelled by a retained net-flux path that reforms "
            "the parent within the requested horizon. It is a network-flow attribution, not a molecular "
            "trajectory or solvent-cage model."
        ),
    }


def _validated_replacement_ids(replacements: Iterable[dict], condition_contract: dict) -> tuple[set[str], list[dict]]:
    verified = set()
    validation = []
    for replacement in replacements:
        report = validate_rate_replacement_evidence(
            replacement,
            temperature_K=_finite(condition_contract.get("temperature_K")),
            pressure_bar=_finite(condition_contract.get("pressure_bar")),
        )
        validation.append(report)
        route_id = replacement.get("route_id")
        if report["status"] == "passed" and route_id:
            verified.add(str(route_id))
    return verified, validation


def _endpoint_labels(route: dict, side: str) -> set[str]:
    endpoints = route.get("resolved_endpoints") or {}
    return {normalize_rmg_label(str(item.get("label"))) for item in endpoints.get(side) or []
            if item.get("label")}


def _generated_reactants(route: dict, initial_labels: set[str]) -> set[str]:
    declared = {
        normalize_rmg_label(str(item.get("label")))
        for item in route.get("co_reactants") or []
        if item.get("label") and item.get("origin") == "generated_intermediate"
    }
    if declared:
        return declared
    # Endpoint fallback is used only when the route has already been classified
    # as dependent on modeled/generated chemistry.  Direct routes often carry
    # no co-reactant records and must not acquire invented dependencies.
    if (route.get("_route_source") != "network_routes" and
            route.get("initiation_status") not in {"thermally_reachable", "requires_generated_intermediate"}):
        return set()
    return _endpoint_labels(route, "reactants") - initial_labels


def _retained_network_routes(rmg_evidence: dict) -> list[dict]:
    raw = rmg_evidence.get("network_routes")
    if raw is None:
        raw = (rmg_evidence.get("mechanism_inspection") or {}).get("reactions") or []
    result = []
    for index, item in enumerate(raw or []):
        if not isinstance(item, dict):
            continue
        route = dict(item)
        route.setdefault("route_id", f"network-{index + 1}")
        route.setdefault("reaction_equation", route.get("equation"))
        if not route.get("resolved_endpoints"):
            reactants = [{"label": label} for label in route.get("reactant_labels") or []]
            products = [{"label": label} for label in route.get("product_labels") or []]
            if reactants and products:
                route["resolved_endpoints"] = {"reactants": reactants, "products": products}
        result.append(route)
    return result


def _producer_extent(route: dict, flux_reactions: list[dict]) -> float:
    """Return integrated extent in the route's declared forward direction."""
    equation = route.get("reaction_equation") or route.get("equation") or ""
    oriented = _oriented_signature(equation)
    normalized = _oriented_signature(equation, normalized=True)
    extent = 0.0
    for record in flux_reactions:
        record_equation = record.get("reaction_equation", "")
        record_oriented = _oriented_signature(record_equation)
        record_normalized = _oriented_signature(record_equation, normalized=True)
        same = oriented is not None and record_oriented == oriented
        reverse = oriented is not None and record_oriented == tuple(reversed(oriented))
        if not same and not reverse:
            same = normalized is not None and record_normalized == normalized
            reverse = normalized is not None and record_normalized == tuple(reversed(normalized))
        if same:
            value = _finite(record.get("integrated_forward_extent_kmol"))
        elif reverse:
            value = _finite(record.get("integrated_reverse_extent_kmol"))
        else:
            continue
        # Older retained target-loss records lack directional extents.  Their
        # explicit target destruction is still valid evidence for that same
        # candidate direction, but cannot rescue a non-target network route.
        if value is None and route.get("_route_source") == "candidate_routes":
            value = _finite(record.get("integrated_target_destruction_kmol"))
        if value is not None and value > 0.0:
            extent += value
    return extent


def _dependency_order(
    controlling_index: int,
    routes: list[dict],
    producer_relevance: list[float],
    initial_labels: set[str],
    verified_ids: set[str],
    deferred_ids: set[str] | None = None,
) -> tuple[list[int], list[str], list[dict]]:
    """Backtrack through viable retained producers for generated reactants.

    A failed ORCA branch is excluded only from the current verification
    session.  The search then tries the next flux-supported producer instead
    of allowing one failed prerequisite to block the whole downstream route.
    """
    deferred_ids = deferred_ids or set()
    diagnostics: list[dict] = []
    state_count = 0
    budget_exhausted = False
    resolved_cache: dict[int, tuple[int, ...]] = {}

    def route_id(route_index: int) -> str:
        return str(routes[route_index].get("route_id") or f"rmg-route-{route_index + 1}")

    producers_by_label: dict[str, list[tuple[bool, float, int]]] = {}
    for candidate_index, candidate in enumerate(routes):
        relevance = producer_relevance[candidate_index]
        if relevance <= 0.0:
            continue
        directly_reachable = not _generated_reactants(candidate, initial_labels)
        for product_label in _endpoint_labels(candidate, "products"):
            producers_by_label.setdefault(product_label, []).append(
                (directly_reachable, relevance, candidate_index),
            )
    for producers in producers_by_label.values():
        producers.sort(key=lambda item: (item[0], item[1], -item[2]), reverse=True)

    def search(route_index: int, visiting: tuple[int, ...]) -> tuple[bool, list[int], set[str]]:
        nonlocal state_count, budget_exhausted
        state_count += 1
        if state_count > MAX_DEPENDENCY_SEARCH_STATES:
            if not budget_exhausted:
                diagnostics.append({
                    "status": "dependency_search_budget_exhausted",
                    "max_states": MAX_DEPENDENCY_SEARCH_STATES,
                })
            budget_exhausted = True
            return False, [], {"dependency_search_budget_exhausted"}
        identity = route_id(route_index)
        if identity in deferred_ids and identity not in verified_ids:
            diagnostics.append({
                "route_id": identity,
                "status": "producer_deferred_after_unverified_path",
            })
            return False, [], set()
        if route_index in visiting:
            diagnostics.append({
                "route_id": identity,
                "status": "dependency_cycle_rejected",
                "cycle_route_ids": [route_id(index) for index in (*visiting, route_index)],
            })
            return False, [], {f"cycle_at_route:{identity}"}
        cached = resolved_cache.get(route_index)
        if cached is not None and not any(index in visiting for index in cached):
            return True, list(cached), set()

        branch_order: list[int] = []
        next_visiting = (*visiting, route_index)
        for generated_label in sorted(_generated_reactants(routes[route_index], initial_labels)):
            producers = []
            rejected_deferred = []
            for directly_reachable, relevance, candidate_index in producers_by_label.get(generated_label, []):
                if candidate_index == route_index:
                    continue
                candidate_id = route_id(candidate_index)
                if candidate_id in deferred_ids and candidate_id not in verified_ids:
                    rejected_deferred.append(candidate_id)
                    continue
                producers.append((directly_reachable, relevance, candidate_index))
            omitted = max(0, len(producers) - MAX_PRODUCERS_PER_INTERMEDIATE)
            producers = producers[:MAX_PRODUCERS_PER_INTERMEDIATE]
            attempted = []
            selected = None
            collected_unresolved: set[str] = set()
            for _, relevance, producer_index in producers:
                producer_id = route_id(producer_index)
                success, producer_order, unresolved = search(producer_index, next_visiting)
                attempted.append({
                    "route_id": producer_id,
                    "producer_integrated_extent_kmol": relevance,
                    "status": "selected" if success else "branch_unresolved",
                    "unresolved_dependencies": sorted(unresolved),
                })
                if success:
                    selected = (producer_index, producer_order)
                    break
                collected_unresolved.update(unresolved)
            diagnostics.append({
                "consumer_route_id": identity,
                "generated_species": generated_label,
                "status": "producer_selected" if selected is not None else "no_viable_producer",
                "selected_producer_route_id": (
                    route_id(selected[0]) if selected is not None else None
                ),
                "deferred_producer_route_ids": sorted(rejected_deferred),
                "producer_candidates_omitted_by_cap": omitted,
                "attempted_producers": attempted,
            })
            if selected is None:
                return False, [], collected_unresolved | {generated_label}
            for index in selected[1]:
                if index not in branch_order:
                    branch_order.append(index)

        if identity not in verified_ids and route_index not in branch_order:
            branch_order.append(route_index)
        resolved_cache[route_index] = tuple(branch_order)
        return True, branch_order, set()

    success, ordered, unresolved = search(controlling_index, ())
    diagnostics.append({
        "status": "dependency_search_completed",
        "states_explored": state_count,
        "max_states": MAX_DEPENDENCY_SEARCH_STATES,
        "budget_exhausted": budget_exhausted,
    })
    return (ordered if success else []), sorted(unresolved), diagnostics


def build_flux_ranked_verification_plan(
    rmg_evidence: dict,
    condition_contract: dict,
    *,
    verified_replacements: Iterable[dict] = (),
    deferred_route_ids: Iterable[str] = (),
    flux_coverage_target: float = DEFAULT_FLUX_COVERAGE_TARGET,
    max_closure_error: float = DEFAULT_MAX_FLUX_CLOSURE_ERROR,
) -> dict:
    """Rank retained target-loss routes and select the next controlling route.

    One route is selected per iteration.  After its verified replacement is
    installed, the mechanism must be propagated and ranked again because flux
    can move to another pathway.
    """
    if not 0.0 < flux_coverage_target <= 1.0:
        raise ValueError("Flux coverage target must be in (0, 1]")
    if not 0.0 <= max_closure_error < 1.0:
        raise ValueError("Maximum flux closure error must be in [0, 1)")
    deferred_ids = {str(route_id) for route_id in deferred_route_ids}
    verified_ids, _ = _validated_replacement_ids(verified_replacements, condition_contract)
    partial_plan = _partial_enumeration_plan(
        rmg_evidence, condition_contract,
        deferred_ids=deferred_ids, verified_ids=verified_ids,
    )
    if partial_plan is not None:
        return partial_plan
    blockers = _solver_gates(rmg_evidence, max_closure_error=max_closure_error)
    if blockers:
        return {
            "schema_version": 1,
            "kind": "flux_ranked_route_verification_plan",
            "status": "insufficient_solver_evidence",
            "blocking_reasons": blockers,
            "selected_route_id": None,
            "ranked_routes": [],
            "rule": "No route is prioritized without completed, closed quantitative target-loss flux evidence.",
        }

    candidate_routes = []
    for index, raw_route in enumerate(rmg_evidence.get("candidate_routes") or []):
        route = dict(raw_route)
        route["_route_source"] = "candidate_routes"
        route["_source_route_index"] = index
        candidate_routes.append(route)
    candidate_signatures = {
        reaction_signature(equation, normalized=True)
        for route in candidate_routes for equation in _route_equations(route)
        if reaction_signature(equation, normalized=True) is not None
    }
    network_routes = []
    retained_network_signatures = set()
    for index, raw_route in enumerate(_retained_network_routes(rmg_evidence)):
        physical_signature = reaction_signature(
            raw_route.get("reaction_equation") or raw_route.get("equation") or "", normalized=True,
        )
        oriented_signature = _oriented_signature(
            raw_route.get("reaction_equation") or raw_route.get("equation") or "", normalized=True,
        )
        # The target candidate carries richer direction/reachability evidence;
        # neither canonical Chemkin orientation of that same reversible core
        # reaction may re-enter as a support route.  Distinct forward/reverse
        # support directions are retained for reactions which are not already
        # target candidates because either direction may produce an initiator.
        if (physical_signature is None or oriented_signature is None
                or physical_signature in candidate_signatures
                or oriented_signature in retained_network_signatures):
            continue
        retained_network_signatures.add(oriented_signature)
        route = dict(raw_route)
        route["_route_source"] = "network_routes"
        route["_source_route_index"] = index
        network_routes.append(route)
    routes = [*candidate_routes, *network_routes]
    propagation = rmg_evidence["independent_cantera_propagation"]
    flux = propagation["reaction_flux_attribution"]
    commitment = assess_target_loss_commitment(
        flux, str(condition_contract.get("target_label") or "stability"),
    )
    commitment_by_index = {
        item.get("reaction_index"): item
        for item in commitment.get("reaction_commitment") or []
    }
    commitment_available = commitment.get("status") in {"completed", "no_committed_target_loss"}
    if commitment_available:
        total = _finite(commitment.get("committed_target_loss_kmol"))
        ranking_basis = "flux_committed_net_target_loss"
    else:
        total = _finite(flux.get("total_integrated_target_destruction_kmol"))
        ranking_basis = "gross_target_destruction_fallback_commitment_unresolved"
    # If RMG's rate screen contains a collision-limit violation, its flux
    # ranking is not safe.  Do not invent chemistry to recover: promote only
    # direct unimolecular reactions already emitted by RMG into an explicit
    # ORCA verification queue.  This is especially important when the invalid
    # route masks an intrinsic target decomposition with low estimated flux.
    intrinsic_failover = (
        _direct_unimolecular_rmg_candidates(
            candidate_routes, str(condition_contract.get("target_label") or "stability"),
        ) if _invalid_screening_kinetics(rmg_evidence) else []
    )
    if total is None or total <= 0.0:
        if intrinsic_failover:
            selectable_failover = [
                item for item in intrinsic_failover
                if str(item[1].get("route_id") or f"rmg-route-{item[0] + 1}") not in deferred_ids
                and str(item[1].get("route_id") or f"rmg-route-{item[0] + 1}") not in verified_ids
            ]
            if not selectable_failover:
                return {
                    "schema_version": 1,
                    "kind": "flux_ranked_route_verification_plan",
                    "status": "all_unverified_routes_deferred",
                    "ranking_basis": "rmg_direct_unimolecular_failover_after_invalid_screening_kinetics",
                    "intrinsic_failover_required": True,
                    "selected_route_id": None,
                    "selected_route_index": None,
                    "selected_route_source": None,
                    "ranked_routes": [],
                    "target_loss_commitment": commitment,
                    "blocking_reasons": [],
                }
            selected_index, selected_route = selectable_failover[0]
            route_id = str(selected_route.get("route_id") or f"rmg-route-{selected_index + 1}")
            return {
                "schema_version": 1,
                "kind": "flux_ranked_route_verification_plan",
                "status": "verification_required",
                "ranking_basis": "rmg_direct_unimolecular_failover_after_invalid_screening_kinetics",
                "selection_reason": "invalid_screening_kinetics_requires_orca_check_of_rmg_direct_unimolecular_routes",
                "intrinsic_failover_required": True,
                "flux_coverage_target": flux_coverage_target,
                "max_flux_closure_error": max_closure_error,
                "selected_route_id": route_id,
                "selected_route_index": selected_index,
                "selected_route_source": "candidate_routes",
                "ranked_routes": [{
                    "rank": position,
                    "route_index": route_index,
                    "route_source": "candidate_routes",
                    "route_id": str(route.get("route_id") or f"rmg-route-{route_index + 1}"),
                    "reaction_equation": route.get("reaction_equation") or route.get("printed_reaction_equation"),
                    "initiation_status": route.get("initiation_status") or "directly_present",
                    "verification_queue_class": "direct_unimolecular_rmg_candidate",
                    "integrated_committed_target_loss_kmol": 0.0,
                    "committed_target_loss_fraction": 0.0,
                    "verified_replacement": False,
                    "deferred_unverified": str(route.get("route_id") or f"rmg-route-{route_index + 1}") in deferred_ids,
                } for position, (route_index, route) in enumerate(intrinsic_failover, 1)],
                "target_loss_commitment": commitment,
                "blocking_reasons": [],
                "requires_repropagation_after_each_verification": True,
                "rule": "RMG supplied these direct unimolecular routes; ORCA must test them because an invalid RMG rate cannot rank the network.",
            }
        return {
            "schema_version": 1, "kind": "flux_ranked_route_verification_plan",
            "status": "no_target_destruction_flux", "blocking_reasons": [],
            "selected_route_id": None, "ranked_routes": [],
            "flux_coverage_target": flux_coverage_target,
            "target_loss_commitment": commitment,
        }

    flux_reactions = []
    for raw_item in flux.get("reactions") or []:
        item = dict(raw_item)
        committed = commitment_by_index.get(item.get("reaction_index")) or {}
        metric = (
            _finite(committed.get("committed_target_loss_kmol"))
            if commitment_available else _finite(item.get("integrated_target_destruction_kmol"))
        )
        item["verification_ranking_target_loss_kmol"] = metric or 0.0
        item["commitment_evidence"] = committed or None
        if (metric or 0.0) > 0.0:
            flux_reactions.append(item)
    matches, unmatched_indices = _flux_matches(candidate_routes, flux_reactions)
    sensitivities = _sensitivity_by_signature(propagation)
    route_flux = [0.0] * len(candidate_routes)
    route_gross_flux = [0.0] * len(candidate_routes)
    route_positive_net_flux = [0.0] * len(candidate_routes)
    route_evidence: list[list[dict]] = [[] for _ in candidate_routes]
    ambiguous_flux = 0.0
    for flux_index, route_indices in matches.items():
        record = flux_reactions[flux_index]
        amount = float(record["verification_ranking_target_loss_kmol"])
        if len(route_indices) != 1:
            ambiguous_flux += amount
            continue
        route_index = route_indices[0]
        route_flux[route_index] += amount
        route_gross_flux[route_index] += _finite(record.get("integrated_target_destruction_kmol")) or 0.0
        route_positive_net_flux[route_index] += max(
            0.0, _finite(record.get("integrated_net_target_loss_kmol")) or 0.0,
        )
        route_evidence[route_index].append(record)
    unmatched_flux = sum(float(flux_reactions[index]["verification_ranking_target_loss_kmol"])
                         for index in unmatched_indices) + ambiguous_flux

    verified_ids, replacement_validation = _validated_replacement_ids(verified_replacements, condition_contract)
    ranked = []
    total_gross = _finite(flux.get("total_integrated_target_destruction_kmol")) or 0.0
    for route_index, (route, amount) in enumerate(zip(candidate_routes, route_flux)):
        if amount <= 0.0:
            continue
        equation = route.get("reaction_equation") or route.get("printed_reaction_equation") or ""
        sensitivity = sensitivities.get(reaction_signature(equation, normalized=True), {})
        sensitivity_value = _finite(sensitivity.get("absolute_log_t95_sensitivity"))
        ranked.append({
            "route_index": route_index,
            "route_source": "candidate_routes",
            "route_id": route.get("route_id") or f"rmg-route-{route_index + 1}",
            "reaction_equation": equation,
            "initiation_status": route.get("initiation_status"),
            "integrated_target_destruction_kmol": route_gross_flux[route_index],
            "integrated_positive_net_target_loss_kmol": route_positive_net_flux[route_index],
            "integrated_committed_target_loss_kmol": amount,
            "gross_target_destruction_fraction": (
                route_gross_flux[route_index] / total_gross if total_gross > 0.0 else None
            ),
            "committed_target_loss_fraction": amount / total,
            "absolute_log_t95_sensitivity": sensitivity_value,
            "sensitivity_evidence": sensitivity or None,
            "flux_reaction_indices": [item.get("reaction_index") for item in route_evidence[route_index]],
            "verified_replacement": str(route.get("route_id") or f"rmg-route-{route_index + 1}") in verified_ids,
            "deferred_unverified": str(route.get("route_id") or f"rmg-route-{route_index + 1}") in deferred_ids,
        })
    # Keep direct RMG unimolecular candidates visible even if the invalid
    # screen assigned them zero/tiny flux.  The RMG-derived route, not a
    # topology-generated substitute, is what ORCA will test.
    present_indices = {item["route_index"] for item in ranked}
    for route_index, route in intrinsic_failover:
        if route_index in present_indices:
            continue
        route_id = str(route.get("route_id") or f"rmg-route-{route_index + 1}")
        ranked.append({
            "route_index": route_index,
            "route_source": "candidate_routes",
            "route_id": route_id,
            "reaction_equation": route.get("reaction_equation") or route.get("printed_reaction_equation") or "",
            "initiation_status": route.get("initiation_status") or "directly_present",
            "verification_queue_class": "direct_unimolecular_rmg_candidate",
            "integrated_target_destruction_kmol": 0.0,
            "integrated_positive_net_target_loss_kmol": 0.0,
            "integrated_committed_target_loss_kmol": 0.0,
            "gross_target_destruction_fraction": 0.0,
            "committed_target_loss_fraction": 0.0,
            "absolute_log_t95_sensitivity": None,
            "sensitivity_evidence": None,
            "flux_reaction_indices": [],
            "verified_replacement": route_id in verified_ids,
            "deferred_unverified": route_id in deferred_ids,
        })
    ranked.sort(key=lambda item: (
        item.get("verification_queue_class") == "direct_unimolecular_rmg_candidate",
        item["absolute_log_t95_sensitivity"] is not None,
        item["absolute_log_t95_sensitivity"] or 0.0,
        item["committed_target_loss_fraction"],
    ), reverse=True)

    cumulative = 0.0
    verified_coverage = 0.0
    for rank, item in enumerate(ranked, 1):
        cumulative += item["committed_target_loss_fraction"]
        item["rank"] = rank
        item["cumulative_mapped_flux_fraction"] = cumulative
        if item["verified_replacement"]:
            verified_coverage += item["committed_target_loss_fraction"]
    attributable_fraction = max(0.0, min(1.0, 1.0 - unmatched_flux / total))
    unresolved = [item for item in ranked if not item["verified_replacement"]]
    selectable_unresolved = [item for item in unresolved if not item["deferred_unverified"]]
    highest_flux_unverified_route = selectable_unresolved[0] if selectable_unresolved else None
    controlling_route = None
    next_route = None
    dependency_order: list[int] = []
    unresolved_dependencies: list[str] = []
    dependency_search: list[dict] = []
    selected_dependency_diagnostics: list[dict] = []
    producer_relevance = [
        _producer_extent(route, flux.get("reactions") or []) for route in routes
    ]
    if verified_coverage < flux_coverage_target:
        initial_labels = {normalize_rmg_label(str(label))
                          for label in (condition_contract.get("composition") or {})}
        all_unresolved_dependencies: set[str] = set()
        for candidate in selectable_unresolved:
            candidate_order, candidate_unresolved, diagnostics = _dependency_order(
                candidate["route_index"], routes, producer_relevance, initial_labels,
                verified_ids, deferred_ids,
            )
            eligible_order = [
                index for index in candidate_order
                if str(routes[index].get("route_id") or f"rmg-route-{index + 1}") not in deferred_ids
                and str(routes[index].get("route_id") or f"rmg-route-{index + 1}") not in verified_ids
            ]
            search_record = {
                "controlling_route_id": candidate["route_id"],
                "status": "viable" if eligible_order and not candidate_unresolved else "branch_unresolved",
                "ordered_route_ids": [
                    str(routes[index].get("route_id") or f"rmg-route-{index + 1}")
                    for index in candidate_order
                ],
                "unresolved_initiation_species": candidate_unresolved,
                "producer_search": diagnostics,
            }
            dependency_search.append(search_record)
            if candidate_unresolved or not eligible_order:
                all_unresolved_dependencies.update(candidate_unresolved)
                continue
            controlling_route = candidate
            dependency_order = candidate_order
            selected_dependency_diagnostics = diagnostics
            selected_index = eligible_order[0]
            selected = routes[selected_index]
            if selected_index < len(candidate_routes):
                next_route = next(
                    (item for item in ranked if item["route_index"] == selected_index),
                    {
                        "route_index": selected["_source_route_index"],
                        "route_source": "candidate_routes",
                        "route_id": selected.get("route_id"),
                        "reaction_equation": selected.get("reaction_equation"),
                        "initiation_status": selected.get("initiation_status"),
                        "integrated_target_destruction_kmol": 0.0,
                        "gross_target_destruction_fraction": 0.0,
                        "committed_target_loss_fraction": 0.0,
                        "producer_integrated_extent_kmol": producer_relevance[selected_index],
                        "verified_replacement": False,
                        "deferred_unverified": False,
                        "prerequisite_for_route_id": candidate["route_id"],
                    },
                )
            else:
                next_route = {
                    "route_index": selected["_source_route_index"],
                    "route_source": "network_routes",
                    "route_id": selected.get("route_id"),
                    "reaction_equation": selected.get("reaction_equation") or selected.get("equation"),
                    "initiation_status": selected.get("initiation_status"),
                    "integrated_target_destruction_kmol": 0.0,
                    "gross_target_destruction_fraction": 0.0,
                    "committed_target_loss_fraction": 0.0,
                    "producer_integrated_extent_kmol": producer_relevance[selected_index],
                    "verified_replacement": False,
                    "deferred_unverified": False,
                    "prerequisite_for_route_id": candidate["route_id"],
                }
            break
        if controlling_route is None:
            unresolved_dependencies = sorted(all_unresolved_dependencies)
    if attributable_fraction + 1e-12 < flux_coverage_target:
        status = "insufficient_route_attribution"
        next_route = None
    elif verified_coverage + 1e-12 >= flux_coverage_target:
        status = "flux_coverage_target_met"
    elif unresolved_dependencies:
        status = "controlling_route_has_unresolved_initiation"
    elif next_route:
        status = "verification_required"
    elif unresolved and all(item["deferred_unverified"] for item in unresolved):
        status = "all_unverified_routes_deferred"
    else:
        status = "insufficient_route_attribution"
    return {
        "schema_version": 1,
        "kind": "flux_ranked_route_verification_plan",
        "status": status,
        "ranking_basis": (
            f"retained_t95_sensitivity_then_{ranking_basis}" if sensitivities else ranking_basis
        ),
        "intrinsic_failover_required": bool(intrinsic_failover),
        "selection_reason": (
            "invalid_screening_kinetics_requires_orca_check_of_rmg_direct_unimolecular_routes"
            if intrinsic_failover else "flux_and_reachability_ranked_route_selection"
        ),
        "flux_coverage_target": flux_coverage_target,
        "max_flux_closure_error": max_closure_error,
        "attributable_flux_fraction": attributable_fraction,
        "unattributed_flux_fraction": max(0.0, 1.0 - attributable_fraction),
        "verified_flux_fraction": verified_coverage,
        "attempted_unverified_route_ids": sorted(deferred_ids),
        "deferred_unverified_flux_fraction": sum(
            item["committed_target_loss_fraction"] for item in unresolved
            if item["deferred_unverified"]
        ),
        "target_loss_commitment": commitment,
        "highest_flux_unverified_route_id": (
            highest_flux_unverified_route["route_id"] if highest_flux_unverified_route else None
        ),
        "controlling_route_id": controlling_route["route_id"] if controlling_route else None,
        "selected_route_id": next_route["route_id"] if next_route else None,
        "selected_route_index": next_route["route_index"] if next_route else None,
        "selected_route_source": next_route.get("route_source") if next_route else None,
        "ordered_dependency_route_ids": [
            str(routes[index].get("route_id") or f"rmg-route-{index + 1}") for index in dependency_order
        ],
        "ordered_dependency_routes": [
            {
                "route_id": str(routes[index].get("route_id") or f"rmg-route-{index + 1}"),
                "route_source": routes[index].get("_route_source"),
                "route_index": routes[index].get("_source_route_index"),
            }
            for index in dependency_order
        ],
        "unresolved_initiation_species": unresolved_dependencies,
        "dependency_branch_search": dependency_search,
        "selected_dependency_search": selected_dependency_diagnostics,
        "ranked_routes": ranked,
        "deferred_route_ids": list(dict.fromkeys([
            route_id for route_id in [
                *[str(routes[index].get("route_id") or f"rmg-route-{index + 1}")
                  for index in dependency_order[1:]],
                *[item["route_id"] for item in unresolved],
            ] if route_id != (next_route["route_id"] if next_route else None)
        ])) if next_route else [],
        "replacement_validation": replacement_validation,
        "requires_repropagation_after_each_verification": True,
        "rule": "Verify one controlling route, install its validated rate, repropagate, then recompute flux before choosing another.",
    }


def _complete_propagation(propagation: dict) -> bool:
    return (propagation.get("status") == "completed" and
            bool((propagation.get("coverage") or {}).get("complete")) and
            bool(propagation.get("target_profile")))


def _final_fraction(propagation: dict) -> float | None:
    profile = propagation.get("target_profile") or []
    return _finite(profile[-1].get("target_fraction_remaining")) if profile else None


def assess_t95_robustness(
    repaired_rmg_evidence: dict,
    verification_plan: dict,
    *,
    retention_fraction: float,
    unresolved_sensitivity_envelope: dict | None = None,
    relative_tolerance: float = DEFAULT_T95_RELATIVE_TOLERANCE,
) -> dict:
    """Decide whether a repaired t95 is insensitive to unresolved routes.

    If any nonzero-flux route remains unresolved, a combined lower/upper-bound
    propagation artifact is mandatory.  Individual one-at-a-time sensitivity
    coefficients are useful for ranking but cannot prove a robust lifetime.
    """
    if not 0.0 < retention_fraction < 1.0:
        raise ValueError("Retention fraction must be in (0, 1)")
    if relative_tolerance < 0.0:
        raise ValueError("Relative t95 tolerance must be nonnegative")
    propagation = repaired_rmg_evidence.get("independent_cantera_propagation") or {}
    if not _complete_propagation(propagation):
        return {"status": "not_robust", "reason": "repaired_propagation_incomplete", "t95_seconds": None}
    if verification_plan.get("status") == "insufficient_solver_evidence":
        return {"status": "not_robust", "reason": "flux_plan_has_insufficient_solver_evidence", "t95_seconds": None}
    ranked = verification_plan.get("ranked_routes") or []
    unresolved_ids = sorted(str(item["route_id"]) for item in ranked
                            if not item.get("verified_replacement") and
                            (_finite(item.get("gross_target_destruction_fraction")) or 0.0) > 0.0)
    if (_finite(verification_plan.get("unattributed_flux_fraction")) or 0.0) > 1e-12:
        unresolved_ids.append("__unattributed_target_loss_flux__")
    baseline_t95 = _finite(propagation.get("estimated_time_to_retention_seconds"))
    if not unresolved_ids:
        return {
            "status": "robust_verified_t95" if baseline_t95 is not None else "robust_no_t95_within_horizon",
            "reason": "all_nonzero_target_loss_flux_routes_have_validated_replacements",
            "t95_seconds": baseline_t95,
            "unresolved_route_ids": [],
        }

    envelope = unresolved_sensitivity_envelope or {}
    if envelope.get("status") != "completed" or not envelope.get("combined_perturbation"):
        return {"status": "not_robust", "reason": "combined_unresolved_route_sensitivity_missing",
                "t95_seconds": None, "unresolved_route_ids": unresolved_ids}
    if sorted(str(value) for value in envelope.get("unresolved_route_ids") or []) != sorted(unresolved_ids):
        return {"status": "not_robust", "reason": "sensitivity_envelope_route_set_mismatch",
                "t95_seconds": None, "unresolved_route_ids": unresolved_ids}
    bounds = envelope.get("rate_bounds") or []
    bounded = {str(item.get("route_id")) for item in bounds
               if item.get("bound_source") and _finite(item.get("lower_rate")) is not None
               and _finite(item.get("upper_rate")) is not None}
    if bounded != set(unresolved_ids):
        return {"status": "not_robust", "reason": "unresolved_rate_bounds_missing_provenance",
                "t95_seconds": None, "unresolved_route_ids": unresolved_ids}
    lower = envelope.get("lower_bound_propagation") or {}
    upper = envelope.get("upper_bound_propagation") or {}
    if not _complete_propagation(lower) or not _complete_propagation(upper):
        return {"status": "not_robust", "reason": "sensitivity_bound_propagation_incomplete",
                "t95_seconds": None, "unresolved_route_ids": unresolved_ids}
    bound_t95 = [_finite(lower.get("estimated_time_to_retention_seconds")),
                 _finite(upper.get("estimated_time_to_retention_seconds"))]
    if baseline_t95 is not None:
        if any(value is None for value in bound_t95):
            return {"status": "not_robust", "reason": "t95_crossing_not_preserved_across_bounds",
                    "t95_seconds": None, "unresolved_route_ids": unresolved_ids}
        values = [baseline_t95, *bound_t95]
        relative_span = (max(values) - min(values)) / baseline_t95
        if relative_span > relative_tolerance:
            return {"status": "not_robust", "reason": "t95_sensitive_to_unresolved_routes",
                    "t95_seconds": None, "relative_span": relative_span,
                    "unresolved_route_ids": unresolved_ids}
        return {"status": "robust_verified_t95", "reason": "combined_rate_bounds_preserve_t95",
                "t95_seconds": baseline_t95, "t95_bounds_seconds": [min(values), max(values)],
                "relative_span": relative_span, "relative_tolerance": relative_tolerance,
                "unresolved_route_ids": unresolved_ids}

    final_fractions = [_final_fraction(item) for item in (propagation, lower, upper)]
    if all(value is not None and value >= retention_fraction for value in final_fractions):
        return {"status": "robust_no_t95_within_horizon",
                "reason": "all_combined_rate_bounds_remain_above_retention_at_complete_horizon",
                "t95_seconds": None, "unresolved_route_ids": unresolved_ids}
    return {"status": "not_robust", "reason": "retention_conclusion_changes_across_unresolved_rate_bounds",
            "t95_seconds": None, "unresolved_route_ids": unresolved_ids}


def _propagation_payload(rmg_evidence: dict, condition_contract: dict) -> dict:
    artifacts = rmg_evidence.get("artifacts") or {}
    required = {
        "chemkin": artifacts.get("chemkin"),
        "dictionary": artifacts.get("species_dictionary"),
    }
    missing = [name for name, path in required.items() if not path]
    if missing:
        raise ValueError("Retained repaired mechanism artifacts are missing: " + ", ".join(missing))
    composition = condition_contract.get("composition") or {}
    required_conditions = ("temperature_K", "pressure_bar", "target_duration_seconds", "retention_fraction")
    absent = [name for name in required_conditions if _finite(condition_contract.get(name)) is None]
    if absent or not composition:
        raise ValueError("Immutable condition contract is incomplete for propagation: " + ", ".join(absent))
    return {
        "operation": "propagate",
        "chemkin": required["chemkin"],
        "dictionary": required["dictionary"],
        "initial_mole_fractions": composition,
        "temperature_K": condition_contract["temperature_K"],
        "pressure_bar": condition_contract["pressure_bar"],
        "target_label": condition_contract.get("target_label", "stability"),
        "target_duration_seconds": condition_contract["target_duration_seconds"],
        "retention_fraction": condition_contract["retention_fraction"],
    }


def run_combined_unresolved_sensitivity(
    repaired_rmg_evidence: dict,
    verification_plan: dict,
    condition_contract: dict,
    *,
    rate_bounds: Iterable[dict],
    rmg_env: str | None,
    timeout_seconds: float = 300.0,
) -> dict:
    """Propagate all unresolved routes at retained lower and upper bounds.

    Bounds and their provenance are caller-supplied evidence.  STORCA merely
    converts each absolute bound to a multiplier of the retained mechanism's
    explicitly supplied baseline rate.  It never invents an uncertainty
    factor.  Two combined runs capture route interactions more conservatively
    than one-at-a-time sensitivities.
    """
    ranked = verification_plan.get("ranked_routes") or []
    unresolved = [item for item in ranked if not item.get("verified_replacement") and
                  (_finite(item.get("gross_target_destruction_fraction")) or 0.0) > 0.0]
    if (_finite(verification_plan.get("unattributed_flux_fraction")) or 0.0) > 1e-12:
        raise ValueError("Cannot perturb unattributed target-loss flux")
    unresolved_by_id = {str(item["route_id"]): item for item in unresolved}
    supplied = {str(item.get("route_id")): item for item in rate_bounds if item.get("route_id")}
    if set(supplied) != set(unresolved_by_id):
        raise ValueError("Rate bounds must cover exactly the unresolved nonzero-flux routes")
    lower_multipliers: dict[str, float] = {}
    upper_multipliers: dict[str, float] = {}
    retained_bounds = []
    for route_id, route in unresolved_by_id.items():
        bound = supplied[route_id]
        baseline = _finite(bound.get("baseline_rate"))
        lower = _finite(bound.get("lower_rate"))
        upper = _finite(bound.get("upper_rate"))
        units = str(bound.get("units") or "").strip()
        source = bound.get("bound_source")
        if baseline is None or baseline <= 0.0 or lower is None or upper is None or lower < 0.0 or upper < lower:
            raise ValueError(f"Invalid absolute rate bounds for unresolved route {route_id}")
        if not units or not source:
            raise ValueError(f"Rate bounds lack units/provenance for unresolved route {route_id}")
        reaction_indices = route.get("flux_reaction_indices") or []
        if not reaction_indices or any(index is None for index in reaction_indices):
            raise ValueError(f"No retained Cantera reaction index for unresolved route {route_id}")
        for reaction_index in reaction_indices:
            key = str(int(reaction_index))
            if key in lower_multipliers:
                raise ValueError(f"Cantera reaction index {key} maps to multiple unresolved routes")
            lower_multipliers[key] = lower / baseline
            upper_multipliers[key] = upper / baseline
        retained_bounds.append({
            "route_id": route_id, "baseline_rate": baseline, "lower_rate": lower, "upper_rate": upper,
            "units": units, "bound_source": source, "reaction_indices": [int(value) for value in reaction_indices],
        })
    payload = _propagation_payload(repaired_rmg_evidence, condition_contract)
    from .rmg_bridge_client import run_rmg_bridge

    lower = run_rmg_bridge({**payload, "reaction_multipliers": lower_multipliers},
                           rmg_env=rmg_env, timeout_seconds=timeout_seconds)
    upper = run_rmg_bridge({**payload, "reaction_multipliers": upper_multipliers},
                           rmg_env=rmg_env, timeout_seconds=timeout_seconds)
    status = "completed" if _complete_propagation(lower) and _complete_propagation(upper) else "incomplete"
    return {
        "schema_version": 1,
        "kind": "combined_unresolved_route_sensitivity",
        "status": status,
        "combined_perturbation": True,
        "unresolved_route_ids": sorted(unresolved_by_id),
        "rate_bounds": retained_bounds,
        "lower_bound_propagation": lower,
        "upper_bound_propagation": upper,
        "interpretation": "These are model responses to retained physical rate bounds, not newly verified reaction rates.",
    }


def assess_repaired_mechanism(
    repaired_rmg_evidence: dict,
    replacements: Iterable[dict],
    condition_contract: dict,
) -> dict:
    """Gate a repropagated mechanism after verified rate replacement."""
    replacements = list(replacements)
    verified_ids, validation = _validated_replacement_ids(replacements, condition_contract)
    applied_rate_validation = [
        verify_applied_replacement(repaired_rmg_evidence, replacement, condition_contract)
        for replacement in replacements
    ]
    expected_ids = {str(item.get("route_id")) for item in replacements if item.get("route_id")}
    applied_ids = {str(item.get("route_id")) for item in
                   (repaired_rmg_evidence.get("generated_kinetics_libraries") or []) if item.get("route_id")}
    hard_blockers = []
    if verified_ids != expected_ids:
        hard_blockers.append("one_or_more_replacement_rates_failed_validation")
    if not expected_ids <= applied_ids:
        hard_blockers.append("validated_replacement_library_not_present_in_repaired_model")
    if any(item["status"] != "passed" for item in applied_rate_validation):
        hard_blockers.append("applied_replacement_not_rate_matched")
    kinetics = repaired_rmg_evidence.get("kinetics_validation") or {}
    collision_violation_remains = bool(
        kinetics.get("status") != "passed" or kinetics.get("violation_count", 0)
    )
    hard_blockers.extend(_solver_gates(repaired_rmg_evidence, max_closure_error=DEFAULT_MAX_FLUX_CLOSURE_ERROR))
    if hard_blockers:
        status = "rejected"
    elif collision_violation_remains:
        # A downstream invalid rate cannot support a terminal lifetime, but it
        # must not prevent reranking after a valid upstream initiation repair.
        status = "accepted_for_flux_reranking_with_unrepaired_collision_violations"
    else:
        status = "accepted_for_flux_reranking"
    blockers = [*hard_blockers]
    if collision_violation_remains:
        blockers.append("collision_limit_violation_remains")
    propagation = repaired_rmg_evidence.get("independent_cantera_propagation") or {}
    return {
        "schema_version": 1,
        "kind": "repaired_kinetic_model_validation",
        "status": status,
        "blocking_reasons": list(dict.fromkeys(blockers)),
        "hard_blocking_reasons": list(dict.fromkeys(hard_blockers)),
        "terminal_kinetics_eligible": not hard_blockers and not collision_violation_remains,
        "replacement_validation": validation,
        "applied_rate_validation": applied_rate_validation,
        "expected_replacement_route_ids": sorted(expected_ids),
        "applied_replacement_route_ids": sorted(applied_ids & expected_ids),
        "estimated_time_to_retention_seconds": (
            propagation.get("estimated_time_to_retention_seconds")
            if not hard_blockers and not collision_violation_remains else None
        ),
        "rule": "A repaired lifetime is unavailable until every replacement is validated, loaded, collision-safe, and repropagated with complete flux evidence.",
    }
