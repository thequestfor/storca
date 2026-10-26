"""Source-neutral discovery adapters and conservative verification hand-off."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping, Protocol, Sequence

from .reaction_candidates import ReactionCandidate, candidate_from_rmg_route, merge_reaction_candidates


DISCOVERY_RESULT_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class DiscoveryResult:
    """Retained output of one bounded discovery backend.

    ``status`` describes the search, not reaction existence.  In particular,
    ``discovery_truncated`` is positive evidence that the configured search
    did not cover its entire nominal edit space.
    """

    backend: str
    status: str
    candidates: tuple[ReactionCandidate, ...] = ()
    records: tuple[dict[str, Any], ...] = ()
    provenance: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": DISCOVERY_RESULT_SCHEMA_VERSION,
            "kind": "storca_reaction_discovery_result",
            "backend": self.backend,
            "status": self.status,
            "candidates": [candidate.as_dict() for candidate in self.candidates],
            "records": [dict(record) for record in self.records],
            "provenance": dict(self.provenance),
            "warnings": list(self.warnings),
            "raw_candidate_count": len(self.records),
            "candidate_count": len(self.candidates),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DiscoveryResult":
        """Restore a deterministic discovery artifact for reuse in another lane."""
        if not isinstance(value, Mapping):
            raise TypeError("Reaction-discovery artifact must be a mapping")
        kind = value.get("kind")
        if kind not in {None, "storca_reaction_discovery_result"}:
            raise ValueError(f"Unsupported reaction-discovery artifact kind: {kind}")
        backend = str(value.get("backend") or "")
        status = str(value.get("status") or "")
        if not backend or not status:
            raise ValueError("Reaction-discovery artifact needs backend and status")
        raw_candidates = value.get("candidates") or ()
        raw_records = value.get("records") or ()
        if any(not isinstance(candidate, Mapping) for candidate in raw_candidates):
            raise ValueError("Reaction-discovery artifact contains a malformed candidate record")
        if any(not isinstance(record, Mapping) for record in raw_records):
            raise ValueError("Reaction-discovery artifact contains a malformed search record")
        return cls(
            backend=backend,
            status=status,
            candidates=tuple(
                ReactionCandidate.from_dict(dict(candidate))
                for candidate in raw_candidates
            ),
            records=tuple(
                dict(record) for record in raw_records
            ),
            provenance=dict(value.get("provenance") or {}),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


class ReactionDiscoveryBackend(Protocol):
    """Small common interface for bounded candidate discovery."""

    name: str

    def discover(
        self,
        system: str | Mapping[str, Any],
        conditions: Mapping[str, Any] | None,
        budget: Any,
        *,
        workdir: Path | None = None,
    ) -> DiscoveryResult:
        ...


def _system_smiles(system: str | Mapping[str, Any]) -> str:
    if isinstance(system, str):
        return system
    for key in ("smiles", "target_smiles", "parent_smiles"):
        value = system.get(key)
        if value:
            return str(value)
    raise ValueError("Discovery system needs a target SMILES string")


class RmgReactionDiscoveryBackend:
    """Adapter over retained RMG candidate-route evidence.

    It intentionally does not run or reconfigure RMG.  The existing bounded
    RMG workflow remains responsible for generation, execution, and all
    kinetic assumptions.
    """

    name = "rmg"

    def __init__(self, rmg_evidence: Mapping[str, Any] | None = None) -> None:
        self.rmg_evidence = dict(rmg_evidence or {})

    def discover(
        self,
        system: str | Mapping[str, Any],
        conditions: Mapping[str, Any] | None = None,
        budget: Any = None,
        *,
        workdir: Path | None = None,
    ) -> DiscoveryResult:
        parent_smiles = _system_smiles(system)
        evidence = dict(self.rmg_evidence or (system.get("rmg_evidence") if isinstance(system, Mapping) else {}) or {})
        candidates: list[ReactionCandidate] = []
        records: list[dict[str, Any]] = []
        from .rmg_partial_candidates import validated_retained_partial_routes
        retained = validated_retained_partial_routes(evidence, dict(conditions or {}))
        routes = [("candidate_routes", index, route)
                  for index, route in enumerate(evidence.get("candidate_routes") or [])]
        routes.extend(("retained_partial_routes", index, route)
                      for index, route in enumerate(retained["routes"]))
        for source, index, route in routes:
            try:
                route_record = dict(route)
                candidate = candidate_from_rmg_route(parent_smiles, route_record)
            except Exception as error:
                route_record = dict(route) if isinstance(route, Mapping) else {}
                candidate = None
                adapter_failure = f"{type(error).__name__}: {error}"
            else:
                adapter_failure = None
            if candidate is None:
                records.append({
                    "route_index": index,
                    "route_source": source,
                    "status": "not_generated",
                    "reason": (
                        "rmg_route_adapter_failed" if adapter_failure else
                        "rmg_route_lacks_resolved_endpoint_smiles"
                    ),
                    "adapter_failure": adapter_failure,
                    "route_id": route_record.get("route_id"),
                })
            else:
                if source == "retained_partial_routes":
                    provenance = {
                        "source_attempt": route_record["source_attempt"],
                        "source_route_index": route_record["source_route_index"],
                        "discovery_only": True,
                    }
                    candidate = replace(
                        candidate,
                        discovery_assertions=tuple({**assertion, **provenance}
                                                   for assertion in candidate.discovery_assertions),
                        provenance={**candidate.provenance, **provenance},
                    )
                candidates.append(candidate)
                records.append({
                    "route_index": index,
                    "route_source": source,
                    "status": "candidate_unverified",
                    "route_id": route_record.get("route_id"),
                    "candidate_id": candidate.candidate_id,
                })
        search_status = "completed" if evidence.get("status") == "completed" else "incomplete"
        return DiscoveryResult(
            backend=self.name,
            status=search_status,
            candidates=tuple(candidates),
            records=tuple(records),
            provenance={
                "rmg_status": evidence.get("status"),
                "retained_attempt_blocking_reasons": retained["blocking_reasons"],
                "retained_attempt_findings": retained["findings"],
                "controller": evidence.get("controller"),
                "condition_contract_id": (conditions or {}).get("contract_id"),
                "adapter": "retained_rmg_candidate_routes",
                "search_scope": "condition_specific_rmg_family_and_library_generation",
                "reaction_classes": ["rmg_retained_bounded_generation"],
                "scope_interpretation": (
                    "RMG searched enabled families and libraries under its retained model-"
                    "generation bounds; missing family coverage remains possible."
                ),
            },
            warnings=(
                "RMG adapter only represents routes retained by the bounded RMG search; absence is not chemical absence.",
            ),
        )


def merge_discovery_results(results: Sequence[DiscoveryResult]) -> DiscoveryResult:
    """Merge candidate assertions while preserving every backend result."""
    candidates = merge_reaction_candidates(
        candidate for result in results for candidate in result.candidates
    )
    statuses = {result.status for result in results}
    status = "discovery_truncated" if "discovery_truncated" in statuses else (
        "completed" if statuses <= {"completed"} else "incomplete"
    )
    return DiscoveryResult(
        backend="merged",
        status=status,
        candidates=tuple(candidates),
        records=tuple(
            {
                "backend": result.backend,
                "status": result.status,
                "candidate_count": len(result.candidates),
                "raw_candidate_count": len(result.records),
            }
            for result in results
        ),
        provenance={
            "backends": [result.backend for result in results],
            "backend_results": [result.as_dict() for result in results],
        },
        warnings=tuple(sorted({warning for result in results for warning in result.warnings})),
    )


def graph_only_candidates_for_verification(
    candidates: Sequence[ReactionCandidate],
    *,
    parent_smiles: str,
    target_label: str = "stability",
    deferred_route_ids: Sequence[str] = (),
    verified_route_ids: Sequence[str] = (),
    maximum_candidates: int | None = None,
) -> list[dict[str, Any]]:
    """Return a bounded, explicitly non-kinetic graph-only verification queue.

    Candidates jointly found by RMG are intentionally left to the existing
    flux/dependency planner. Graph-only candidates have unknown kinetic
    relevance; ordering is solely by surviving xTB reaction energy then a
    deterministic ID, not a fabricated rate or barrier score.
    """
    deferred = {str(value) for value in deferred_route_ids}
    verified = {str(value) for value in verified_route_ids}
    eligible = [
        candidate for candidate in candidates
        if "graph_xtb" in candidate.discovery_sources
        and "rmg" not in candidate.discovery_sources
        and candidate.screening_status == "candidate_unverified"
        and candidate.verification_status == "candidate_unverified"
        and candidate.candidate_id not in deferred
        and candidate.candidate_id not in verified
    ]
    eligible.sort(key=lambda candidate: (
        candidate.estimated_reaction_energy_kcal_mol is None,
        candidate.estimated_reaction_energy_kcal_mol
        if candidate.estimated_reaction_energy_kcal_mol is not None else float("inf"),
        candidate.candidate_id,
    ))
    if maximum_candidates is not None:
        eligible = eligible[:max(0, int(maximum_candidates))]
    routes = []
    for rank, candidate in enumerate(eligible, 1):
        route = candidate.as_route_dict(parent_smiles=parent_smiles, target_label=target_label)
        route.update({
            "graph_only_rank": rank,
            "verification_queue_class": "graph_xtb_screened_kinetic_relevance_unknown",
            "selection_reason": (
                "Graph candidate survived bounded structural and xTB reaction-energy screening; "
                "kinetic relevance is unknown and no xTB rate/barrier was assigned."
            ),
        })
        routes.append(route)
    return routes
def graph_candidate_verification_plan(routes: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Make one supplemental route consumable by :class:`VerificationEngine`."""
    if not routes:
        return {
            "schema_version": 1,
            "kind": "graph_xtb_candidate_verification_plan",
            "status": "no_graph_candidate_ready",
            "selected_route_id": None,
            "ranked_routes": [],
            "rule": "No graph-only route is selected without surviving bounded xTB screening.",
        }
    selected = dict(routes[0])
    unexpected = bool(selected.get("unexpected_discovery_evidence"))
    earlier_rmg = selected.get("verification_queue_class") == "earlier_rmg_attempt_discovery_kinetic_relevance_unknown"
    ranked = [
        {
            "rank": int(route.get("graph_only_rank") or index),
            "route_id": route.get("route_id"),
            "reaction_equation": route.get("reaction_equation"),
            "estimated_reaction_energy_kcal_mol": route.get("estimated_reaction_energy_kcal_mol"),
            "screening_status": route.get("screening_status"),
            "kinetic_relevance": "unknown",
            "verification_queue_class": (
                "earlier_rmg_attempt_discovery_kinetic_relevance_unknown"
                if route.get("verification_queue_class") == "earlier_rmg_attempt_discovery_kinetic_relevance_unknown" else
                "unexpected_irc_channel_kinetic_relevance_unknown"
                if route.get("unexpected_discovery_evidence") else "graph_xtb_screened_kinetic_relevance_unknown"
            ),
        }
        for index, route in enumerate(routes, 1)
    ]
    return {
        "schema_version": 1,
        "kind": "earlier_rmg_candidate_verification_plan" if earlier_rmg else "unexpected_irc_candidate_verification_plan" if unexpected else "graph_xtb_candidate_verification_plan",
        "status": "verification_required",
        "ranking_basis": ("bounded_earlier_rmg_structural_discovery_not_kinetics" if earlier_rmg else "bounded_unexpected_irc_discovery_not_kinetics" if unexpected
                          else "bounded_graph_xtb_reaction_energy_screen_not_kinetics"),
        "selected_route_id": selected["route_id"],
        "selected_route_source": "retained_partial_routes" if earlier_rmg else "supplemental_unexpected_irc_candidates" if unexpected else "supplemental_graph_candidates",
        **({"selected_route_index": selected["retained_partial_route_index"], "observable_claim_allowed": False}
           if earlier_rmg else {}),
        "selected_route": selected,
        "ranked_routes": ranked,
        "requires_repropagation_after_each_verification": True,
        "rule": (
            "A discovery survivor is a candidate hypothesis. It has unknown kinetic relevance and must pass the "
            "existing ORCA TS/frequency/IRC and Arkane/repropagation gates before it can affect a lifetime."
        ),
    }


__all__ = [
    "DISCOVERY_RESULT_SCHEMA_VERSION", "DiscoveryResult", "ReactionDiscoveryBackend",
    "RmgReactionDiscoveryBackend", "graph_candidate_verification_plan",
    "graph_only_candidates_for_verification", "merge_discovery_results",
]
