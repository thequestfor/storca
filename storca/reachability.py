"""Triage RMG routes by declared initiation conditions and modeled relevance."""

from __future__ import annotations

from .conditions import ConditionSpec
from .rmg_evidence import normalize_rmg_label, time_to_retention_seconds


RADICAL_PRESENCE_FLOOR = 1e-12


def _profile_record(profile: dict | None, normalized_label: str) -> dict | None:
    if not profile:
        return None
    for label, record in profile["species"].items():
        if normalize_rmg_label(label) == normalized_label:
            return record
    return None


def enrich_candidate_route(
    route: dict,
    *,
    species_dictionary: dict[str, dict],
    profile: dict | None,
    conditions: ConditionSpec,
) -> dict:
    """Resolve route species and determine whether its initiation is credible."""
    initial_labels = set(dict(conditions.composition))
    co_reactants = []
    generated_radical = False
    for reactant in route.get("co_reactants", []):
        normalized = normalize_rmg_label(reactant["label"])
        structure = species_dictionary.get(reactant["label"]) or species_dictionary.get(normalized)
        profile_record = _profile_record(profile, normalized)
        origin = "initial_scenario_species" if normalized in initial_labels else "generated_intermediate"
        maximum = profile_record["maximum_mole_fraction"] if profile_record else None
        is_radical = bool(structure and structure["is_radical"])
        if origin == "generated_intermediate" and is_radical:
            generated_radical = True
        co_reactants.append({
            **reactant,
            "normalized_label": normalized,
            "origin": origin,
            "species": structure,
            "maximum_modeled_mole_fraction": maximum,
            "modeled_above_presence_floor": maximum is not None and maximum >= RADICAL_PRESENCE_FLOOR,
        })

    if not co_reactants:
        initiation_status = "directly_present"
    elif all(item["origin"] == "initial_scenario_species" for item in co_reactants):
        initiation_status = "directly_present"
    elif generated_radical and conditions.radical_sources == ("none",):
        if any(item["modeled_above_presence_floor"] for item in co_reactants if item["origin"] == "generated_intermediate"):
            initiation_status = "thermally_reachable"
        else:
            initiation_status = "requires_generated_intermediate"
    elif generated_radical:
        initiation_status = "requires_declared_initiator"
    else:
        initiation_status = "unresolved"

    retention_loss_limit = 1.0 - conditions.retention_fraction
    system_loss = profile["target_loss_fraction"] if profile else None
    if initiation_status in {"requires_generated_intermediate", "requires_declared_initiator", "unresolved"}:
        relevance = "no_credible_initiation_demonstrated"
    elif system_loss is not None and system_loss < retention_loss_limit:
        relevance = "reachable_but_below_loss_threshold"
    else:
        relevance = "kinetically_relevant_candidate"
    return {
        **route,
        "co_reactants": co_reactants,
        "initiation_status": initiation_status,
        "kinetic_relevance": relevance,
        "radical_presence_floor_mole_fraction": RADICAL_PRESENCE_FLOOR,
    }


def assess_kinetic_relevance(routes: list[dict], profile: dict | None, conditions: ConditionSpec) -> dict:
    """Summarize whether any candidate can support modeled loss above the threshold."""
    relevant = [route for route in routes if route["kinetic_relevance"] == "kinetically_relevant_candidate"]
    t_retention = time_to_retention_seconds(profile, conditions.retention_fraction)
    if relevant:
        return {
            "status": "kinetically_relevant_candidate",
            "relevant_route_count": len(relevant),
            "estimated_time_to_retention_seconds": t_retention,
        }
    reachable = [route for route in routes if route["initiation_status"] in {"directly_present", "thermally_reachable"}]
    if reachable:
        return {
            "status": "reachable_but_below_loss_threshold",
            "relevant_route_count": 0,
            "estimated_time_to_retention_seconds": t_retention,
        }
    if routes:
        return {"status": "no_credible_initiation_demonstrated", "relevant_route_count": 0,
                "estimated_time_to_retention_seconds": t_retention}
    return {"status": "no_candidate_route_found", "relevant_route_count": 0,
            "estimated_time_to_retention_seconds": t_retention}
