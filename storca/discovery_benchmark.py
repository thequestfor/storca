"""Small, evidence-oriented scaffolding for reaction-discovery benchmarks.

This module intentionally scores candidate coverage, not kinetic accuracy.
It is a harness for comparing RMG-only, graph/xTB-only, and merged candidate
pools against explicitly curated net-reaction fixtures; it does not assert
that a candidate is a validated elementary mechanism.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from .reaction_candidates import ReactionCandidate, reaction_net_identity
from .reaction_discovery import DiscoveryResult, merge_discovery_results


DISCOVERY_BENCHMARK_SCHEMA_VERSION = 1


def _identity_key(identity: tuple) -> str:
    return json.dumps(identity, separators=(",", ":"), sort_keys=True)


@dataclass(frozen=True)
class DiscoveryBenchmarkCase:
    """A curated candidate-recall fixture, never a kinetic ground truth."""

    identifier: str
    system: dict[str, Any]
    expected_reactions: tuple[dict[str, Any], ...]
    notes: str = ""

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DiscoveryBenchmarkCase":
        return cls(
            identifier=str(value["id"]),
            system=dict(value.get("system") or {}),
            expected_reactions=tuple(dict(item) for item in value.get("expected_reactions") or ()),
            notes=str(value.get("notes") or ""),
        )

    @property
    def expected_identities(self) -> set[str]:
        return {
            _identity_key(reaction_net_identity(item["reactants"], item["products"]))
            for item in self.expected_reactions
        }


def load_benchmark_manifest(path: Path) -> tuple[DiscoveryBenchmarkCase, ...]:
    payload = json.loads(Path(path).read_text())
    if int(payload.get("schema_version", 0)) != DISCOVERY_BENCHMARK_SCHEMA_VERSION:
        raise ValueError("Unsupported reaction-discovery benchmark manifest schema")
    return tuple(DiscoveryBenchmarkCase.from_dict(item) for item in payload.get("cases") or ())


def _score_pool(case: DiscoveryBenchmarkCase, candidates: Sequence[ReactionCandidate]) -> dict[str, Any]:
    expected = case.expected_identities
    found = {_identity_key(candidate.net_identity) for candidate in candidates}
    true_positive = sorted(expected & found)
    false_positive = sorted(found - expected)
    missing = sorted(expected - found)
    return {
        "expected_reaction_count": len(expected),
        "candidate_count": len(candidates),
        "mechanism_recall": len(true_positive) / len(expected) if expected else None,
        "candidate_precision": len(true_positive) / len(found) if found else None,
        "matched_expected_reactions": true_positive,
        "missed_expected_reactions": missing,
        "unmatched_candidate_reactions": false_positive,
        "candidate_sources": {
            candidate.candidate_id: list(candidate.discovery_sources)
            for candidate in candidates
        },
    }


def compare_discovery_results(
    case: DiscoveryBenchmarkCase,
    rmg_result: DiscoveryResult,
    graph_xtb_result: DiscoveryResult,
) -> dict[str, Any]:
    """Summarize A/RMG, B/graph-xTB, and C/merged candidate coverage."""
    merged = merge_discovery_results((rmg_result, graph_xtb_result))
    rmg_score = _score_pool(case, rmg_result.candidates)
    graph_score = _score_pool(case, graph_xtb_result.candidates)
    merged_score = _score_pool(case, merged.candidates)
    graph_statuses = {}
    for candidate in graph_xtb_result.candidates:
        graph_statuses[candidate.screening_status] = graph_statuses.get(candidate.screening_status, 0) + 1
    return {
        "schema_version": DISCOVERY_BENCHMARK_SCHEMA_VERSION,
        "kind": "storca_reaction_discovery_benchmark_result",
        "case_id": case.identifier,
        "case_notes": case.notes,
        "A_rmg_only": rmg_score,
        "B_graph_xtb_only": {
            **graph_score,
            "search_status": graph_xtb_result.status,
            "raw_candidates": len(graph_xtb_result.records),
            "structurally_valid_candidates": sum(
                1 for record in graph_xtb_result.records if record.get("status") == "structurally_valid"
            ),
            "xTB_endpoint_optimizations": (graph_xtb_result.provenance or {}).get("xtb_endpoint_optimizations"),
            "xTB_screening_status_counts": graph_statuses,
        },
        "C_merged": {
            **merged_score,
            "source_agreement_count": sum(
                1 for candidate in merged.candidates
                if {"rmg", "graph_xtb"} <= set(candidate.discovery_sources)
            ),
            "rmg_only_count": sum(
                1 for candidate in merged.candidates
                if "rmg" in candidate.discovery_sources
                and "graph_xtb" not in candidate.discovery_sources
            ),
            "graph_xtb_only_count": sum(
                1 for candidate in merged.candidates
                if "graph_xtb" in candidate.discovery_sources
                and "rmg" not in candidate.discovery_sources
            ),
            "graph_only_candidates_eligible_for_high_level_verification": sum(
                1 for candidate in merged.candidates
                if "graph_xtb" in candidate.discovery_sources
                and "rmg" not in candidate.discovery_sources
                and candidate.screening_status == "candidate_unverified"
            ),
            "candidate_inflation_over_rmg": (
                len(merged.candidates) - len(rmg_result.candidates)
            ),
            "recall_gain_over_rmg": (
                merged_score["mechanism_recall"] - rmg_score["mechanism_recall"]
                if merged_score["mechanism_recall"] is not None
                and rmg_score["mechanism_recall"] is not None else None
            ),
        },
        "interpretation": (
            "Candidate recall is not mechanism validation. A match only records that a backend proposed the curated "
            "net connectivity; TS, rate, and network evidence are evaluated separately."
        ),
    }


__all__ = [
    "DISCOVERY_BENCHMARK_SCHEMA_VERSION", "DiscoveryBenchmarkCase", "compare_discovery_results",
    "load_benchmark_manifest",
]
