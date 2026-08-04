"""Dependency-ordered selection of RMG routes for ORCA verification."""

from __future__ import annotations

import re

from .rmg_evidence import normalize_rmg_label


def _labels(side: list[dict]) -> set[str]:
    return {str(item.get("label")) for item in side if item.get("label")}


def _collision_reactants(equation: str) -> set[str]:
    arrow = "<=>" if "<=>" in equation else "=>"
    left = equation.split(arrow, 1)[0]
    return {re.sub(r"^\s*\d+(?:\.\d+)?\s*", "", token).strip() for token in left.split("+")}


def build_verification_dependency_plan(rmg_evidence: dict, condition_contract: dict) -> dict:
    """Select upstream generation routes before downstream collision violators."""
    routes = rmg_evidence.get("candidate_routes", [])
    initial = set((condition_contract.get("composition") or {}).keys())
    collision_routes = []
    prerequisites = []
    seen = set()
    for violator in ((rmg_evidence.get("kinetics_validation") or {}).get("violators") or
                     rmg_evidence.get("collision_rate_violators") or []):
        equation = violator.get("reaction_equation", "")
        generated = {label for label in _collision_reactants(equation)
                     if normalize_rmg_label(label) not in initial}
        blocked = {"reaction_equation": equation, "generated_reactants": sorted(generated)}
        collision_routes.append(blocked)
        for generated_label in generated:
            for index, route in enumerate(routes):
                endpoints = route.get("resolved_endpoints") or {}
                products = _labels(endpoints.get("products", []))
                reactants = endpoints.get("reactants", [])
                normalized_reactants = {normalize_rmg_label(item.get("label", "")) for item in reactants}
                directly_reachable = normalized_reactants and normalized_reactants <= initial
                consumes_target = "stability" in normalized_reactants
                if generated_label in products and directly_reachable and consumes_target and index not in seen:
                    seen.add(index)
                    prerequisites.append({"route_index": index, "route_id": route.get("route_id"),
                                          "reaction_equation": route.get("reaction_equation"),
                                          "produces": generated_label,
                                          "initiation_status": "directly_present"})
    if prerequisites:
        return {"status": "upstream_initiation_verification_required",
                "selected_route_index": prerequisites[0]["route_index"],
                "ordered_prerequisites": prerequisites,
                "blocked_downstream_routes": collision_routes,
                "rule": "Verify the earliest route from declared initial species before downstream generated-intermediate chemistry."}
    return {"status": "collision_route_has_no_resolved_upstream_prerequisite",
            "selected_route_index": None, "ordered_prerequisites": [],
            "blocked_downstream_routes": collision_routes}
