"""Observable-conditioned mechanism-coverage evidence.

This module evaluates what bounded discovery and kinetic verification evidence
can support.  It deliberately does not estimate a percentage probability that
the true mechanism is complete.  A completed search is only complete within
its retained backend, reaction-class, budget, environment, and condition
contract.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping


DISCOVERY_COVERAGE_SCHEMA_VERSION = 3
_OBSERVABLE_TERMINAL_STATUSES = {
    "verified_t95_converged",
    "verified_below_loss_threshold",
    "no_target_loss_in_completed_rmg_model",
}


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            str(key): _jsonable(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def _content_id(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        _jsonable(value), sort_keys=True, separators=(",", ":"),
        allow_nan=False, default=str,
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _assembly_dict(value: Any) -> dict[str, Any]:
    if isinstance(value, Mapping):
        return dict(value)
    as_dict = getattr(value, "as_dict", None)
    if callable(as_dict):
        restored = as_dict()
        if isinstance(restored, Mapping):
            return dict(restored)
    raise TypeError("Coverage assessment needs a discovery assembly outcome or mapping")


@dataclass(frozen=True)
class CoveragePolicy:
    """Auditable decision policy for one observable and environment lane."""

    policy_name: str
    observable_kind: str
    environment_kind: str
    required_backends: tuple[str, ...]
    required_reaction_classes: tuple[str, ...] = ()
    require_independent_discovery: bool = True
    minimum_independent_backend_count: int = 2
    block_empty_shared_search_without_calibration: bool = True
    require_search_expansion: bool = False
    allow_observable_claim_within_declared_scope: bool = True
    benchmark_profile: str | None = None
    calibration_status: str = "not_blind_spot_calibrated"
    declared_scope_notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.policy_name.strip():
            raise ValueError("CoveragePolicy needs policy_name")
        if not self.observable_kind.strip() or not self.environment_kind.strip():
            raise ValueError("CoveragePolicy needs observable_kind and environment_kind")
        backends = tuple(sorted({str(item) for item in self.required_backends}))
        if not backends:
            raise ValueError("CoveragePolicy needs at least one required backend")
        if self.minimum_independent_backend_count < 1:
            raise ValueError("minimum_independent_backend_count must be positive")
        if self.minimum_independent_backend_count > len(backends):
            raise ValueError(
                "minimum_independent_backend_count cannot exceed required_backends"
            )
        object.__setattr__(self, "required_backends", backends)
        object.__setattr__(
            self, "required_reaction_classes",
            tuple(sorted({str(item) for item in self.required_reaction_classes})),
        )
        object.__setattr__(
            self, "declared_scope_notes",
            tuple(str(item) for item in self.declared_scope_notes),
        )

    @classmethod
    def for_assembly(
        cls,
        assembly: Any,
        *,
        observable_kind: str,
        environment_kind: str,
        benchmark_profile: str | None = None,
    ) -> "CoveragePolicy":
        retained = _assembly_dict(assembly)
        contract = retained.get("contract") or {}
        return cls(
            policy_name="storca_discovery_disagreement_coverage_risk_v2",
            observable_kind=observable_kind,
            environment_kind=environment_kind,
            required_backends=tuple(sorted({
                *(contract.get("required_backends") or ("rmg",)),
                "rmg", "graph_xtb",
            })),
            benchmark_profile=benchmark_profile,
            declared_scope_notes=(
                "A positive decision is conditional on the declared bounded searches and is not a universal mechanism-completeness claim.",
                "Independent discovery agreement is supporting evidence, not proof that omitted chemistry does not exist.",
            ),
        )

    def as_dict(self) -> dict[str, Any]:
        payload = {
            "schema_version": DISCOVERY_COVERAGE_SCHEMA_VERSION,
            "kind": "storca_discovery_coverage_policy",
            "policy_name": self.policy_name,
            "observable_kind": self.observable_kind,
            "environment_kind": self.environment_kind,
            "required_backends": list(self.required_backends),
            "required_reaction_classes": list(self.required_reaction_classes),
            "require_independent_discovery": self.require_independent_discovery,
            "minimum_independent_backend_count": self.minimum_independent_backend_count,
            "block_empty_shared_search_without_calibration": (
                self.block_empty_shared_search_without_calibration
            ),
            "require_search_expansion": self.require_search_expansion,
            "allow_observable_claim_within_declared_scope": (
                self.allow_observable_claim_within_declared_scope
            ),
            "benchmark_profile": self.benchmark_profile,
            "calibration_status": self.calibration_status,
            "declared_scope_notes": list(self.declared_scope_notes),
        }
        return {**payload, "policy_id": _content_id(payload)}


def compare_discovery_expansion(base: Any, expanded: Any) -> dict[str, Any]:
    """Compare two bounded assemblies without treating saturation as completeness."""
    base_value = _assembly_dict(base)
    expanded_value = _assembly_dict(expanded)
    base_candidates = {
        str(item.get("candidate_id"))
        for item in (base_value.get("merged_result") or {}).get("candidates") or []
        if item.get("candidate_id")
    }
    expanded_candidates = {
        str(item.get("candidate_id"))
        for item in (expanded_value.get("merged_result") or {}).get("candidates") or []
        if item.get("candidate_id")
    }
    new_ids = sorted(expanded_candidates - base_candidates)
    missing_ids = sorted(base_candidates - expanded_candidates)
    comparable = (
        base_value.get("status") == "completed_declared_search"
        and expanded_value.get("status") == "completed_declared_search"
    )
    if not comparable:
        status = "not_comparable_incomplete_search"
    elif new_ids or missing_ids:
        status = "candidate_set_changed_under_expansion"
    else:
        status = "stable_within_declared_expansion"
    payload = {
        "schema_version": DISCOVERY_COVERAGE_SCHEMA_VERSION,
        "kind": "storca_discovery_search_expansion_evidence",
        "status": status,
        "base_contract_id": (base_value.get("contract") or {}).get("contract_id"),
        "expanded_contract_id": (expanded_value.get("contract") or {}).get("contract_id"),
        "base_candidate_count": len(base_candidates),
        "expanded_candidate_count": len(expanded_candidates),
        "new_candidate_ids": new_ids,
        "missing_candidate_ids": missing_ids,
        "mechanism_completeness_claimed": False,
        "interpretation": (
            "Candidate-set stability applies only to this declared budget expansion; "
            "it does not establish reaction-network completeness."
        ),
    }
    return {**payload, "evidence_id": _content_id(payload)}


def _candidate_reaction_classes(candidate: Mapping[str, Any]) -> set[str]:
    classes: set[str] = set()
    for assertion in candidate.get("discovery_assertions") or []:
        if not isinstance(assertion, Mapping):
            continue
        edit = assertion.get("edit") or {}
        if isinstance(edit, Mapping) and edit.get("kind"):
            classes.add(str(edit["kind"]))
        if assertion.get("reaction_family"):
            classes.add("rmg_family:" + str(assertion["reaction_family"]))
        if assertion.get("source") == "rmg":
            classes.add("rmg_retained_bounded_generation")
    return classes or {"unclassified_connectivity_hypothesis"}


def _reaction_class_matrix(assembly: Mapping[str, Any]) -> list[dict[str, Any]]:
    searched_by: dict[str, set[str]] = {}
    found: dict[str, set[str]] = {}
    agreement: dict[str, set[str]] = {}
    for result in assembly.get("backend_results") or []:
        backend = str(result.get("backend") or "unknown")
        provenance = result.get("provenance") or {}
        classes = list(provenance.get("reaction_classes") or [])
        if backend == "rmg":
            classes.append("rmg_retained_bounded_generation")
        for reaction_class in classes:
            searched_by.setdefault(str(reaction_class), set()).add(backend)
    candidates = (assembly.get("merged_result") or {}).get("candidates") or []
    for candidate in candidates:
        candidate_id = str(candidate.get("candidate_id") or "")
        sources = {str(item) for item in candidate.get("discovery_sources") or []}
        for reaction_class in _candidate_reaction_classes(candidate):
            found.setdefault(reaction_class, set()).add(candidate_id)
            if {"rmg", "graph_xtb"} <= sources:
                agreement.setdefault(reaction_class, set()).add(candidate_id)
    rows = []
    for reaction_class in sorted(set(searched_by) | set(found)):
        rows.append({
            "reaction_class": reaction_class,
            "searched_by": sorted(searched_by.get(reaction_class, set())),
            "candidate_ids": sorted(found.get(reaction_class, set())),
            "source_agreement_candidate_ids": sorted(
                agreement.get(reaction_class, set())
            ),
            "interpretation": (
                "searched_by reports explicit backend scope only; an empty candidate "
                "set means not found within that search, not chemically absent."
            ),
        })
    return rows


def _coarse_candidate_classes(candidate: Mapping[str, Any]) -> set[str]:
    """Infer only topology classes needed to decide backend search opportunity."""
    classes = set(_candidate_reaction_classes(candidate))
    reactants = tuple(candidate.get("reactants") or ())
    products = tuple(candidate.get("products") or ())
    if len(reactants) == 1 and len(products) > 1:
        classes.add("single_bond_break")
        classes.add("coarse_topology:fragmentation")
    elif len(reactants) > 1 and len(products) == 1:
        classes.add("single_bond_formation")
        classes.add("coarse_topology:association")
    elif len(reactants) == 1 and len(products) == 1:
        classes.add("coarse_topology:unimolecular_rearrangement")
    else:
        classes.add("coarse_topology:multifragment_rearrangement")
    return classes


def analyze_discovery_disagreement(assembly: Any) -> dict[str, Any]:
    """Retain symmetric source agreement and disagreement without scoring truth.

    A backend absence is called a disagreement only when that backend completed
    and the candidate lies inside its explicitly retained scope.  This avoids
    treating "not searched" as "searched and not found".
    """
    retained = _assembly_dict(assembly)
    backend_assessments = retained.get("backend_assessments") or {}
    backend_scopes: dict[str, dict[str, Any]] = {}
    for result in retained.get("backend_results") or []:
        backend = str(result.get("backend") or "unknown")
        provenance = result.get("provenance") or {}
        backend_scopes[backend] = {
            "backend": backend,
            "assessment_status": (backend_assessments.get(backend) or {}).get("status"),
            "search_scope": provenance.get("search_scope"),
            "reaction_classes": sorted({
                str(item) for item in provenance.get("reaction_classes") or ()
            }),
            "condition_contract_id": provenance.get("condition_contract_id"),
        }
    completed = {
        backend for backend, scope in backend_scopes.items()
        if scope.get("assessment_status") in {
            "completed_declared_search", "completed_no_candidates_in_declared_search",
        }
    }

    rows: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    for candidate in (retained.get("merged_result") or {}).get("candidates") or []:
        candidate_id = str(candidate.get("candidate_id") or "")
        sources = sorted({str(item) for item in candidate.get("discovery_sources") or ()})
        classes = _coarse_candidate_classes(candidate)
        eligible: set[str] = set(sources)
        for backend, scope in backend_scopes.items():
            declared = set(scope.get("reaction_classes") or ())
            if declared & classes:
                eligible.add(backend)
            elif backend == "rmg" and backend in completed:
                # RMG is configured as a broad family/library generator.  This
                # means graph-only chemistry is genuine source disagreement,
                # not proof that RMG was required to contain a matching family.
                eligible.add(backend)
        missing = sorted((eligible & completed) - set(sources))
        if len(sources) >= 2:
            status = "source_agreement"
        elif missing:
            status = "source_disagreement"
        else:
            status = "single_source_outside_other_declared_scope"
        counts[status] = counts.get(status, 0) + 1
        rows.append({
            "candidate_id": candidate_id,
            "status": status,
            "found_by": sources,
            "eligible_completed_backends": sorted(eligible & completed),
            "missing_eligible_backends": missing,
            "reaction_classes": sorted(classes),
            "screening_status": candidate.get("screening_status"),
            "interpretation": (
                "Source disagreement is a mechanism-coverage risk signal. It does not "
                "identify which backend is chemically correct."
            ),
        })

    merged_candidates = (retained.get("merged_result") or {}).get("candidates") or []
    empty_shared_search = not merged_candidates and len(completed) >= 2
    independence = {
        "status": (
            "methodologically_distinct_with_shared_representation_dependencies"
            if {"rmg", "graph_xtb"} <= set(backend_scopes)
            else "insufficient_independent_backends"
        ),
        "backend_pairs": ([{
            "backends": ["graph_xtb", "rmg"],
            "distinct_discovery_logic": True,
            "shared_dependencies": [
                "molecular_graph_representation",
                "endpoint_identity_and_state_adaptation",
            ],
            "not_independent_ground_truth": True,
        }] if {"rmg", "graph_xtb"} <= set(backend_scopes) else []),
    }
    payload = {
        "schema_version": DISCOVERY_COVERAGE_SCHEMA_VERSION,
        "kind": "storca_discovery_disagreement_evidence",
        "assembly_contract_id": (retained.get("contract") or {}).get("contract_id"),
        "backend_scopes": [backend_scopes[key] for key in sorted(backend_scopes)],
        "completed_backend_count": len(completed),
        "completed_backends": sorted(completed),
        "candidate_count": len(merged_candidates),
        "candidate_disagreements": rows,
        "status_counts": {key: counts[key] for key in sorted(counts)},
        "empty_shared_search": empty_shared_search,
        "independence_assessment": independence,
        "mechanism_completeness_claimed": False,
        "interpretation": (
            "Agreement and disagreement are conditioned on explicit backend scopes. "
            "Agreement cannot exclude a shared blind spot, and an empty shared search "
            "is not evidence that no chemistry exists."
        ),
    }
    return {**payload, "evidence_id": _content_id(payload)}


def _independent_scope_gaps(assembly: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Identify high-level chemistry no independent non-RMG backend searched."""
    contract = assembly.get("contract") or {}
    system = contract.get("system") or {}
    rmg_contract = contract.get("rmg_contract") or {}
    scenario = rmg_contract.get("scenario") or {}
    additional_species = scenario.get("additional_species") or []
    reactive_partners = sorted({
        str(item.get("label") or item.get("smiles") or "unknown")
        for item in additional_species
        if isinstance(item, Mapping) and item.get("reactive") is True
    })
    non_rmg_scopes = [
        str((result.get("provenance") or {}).get("search_scope") or "")
        for result in assembly.get("backend_results") or []
        if result.get("backend") != "rmg"
    ]
    gaps: list[dict[str, Any]] = []
    if reactive_partners and not any("bimolecular" in scope for scope in non_rmg_scopes):
        gaps.append({
            "code": "independent_bimolecular_discovery_scope_missing",
            "reactive_partners": reactive_partners,
            "reason": (
                "The condition contains reactive co-species, but the independent backend "
                "searched only target-unimolecular graph edits."
            ),
            "required_next_evidence": (
                "an independent condition-aware bimolecular discovery source or a "
                "validated restricted-scope justification"
            ),
        })
    if int(system.get("multiplicity", 1)) > 1 and not any(
        "bimolecular" in scope or "radical" in scope for scope in non_rmg_scopes
    ):
        gaps.append({
            "code": "independent_radical_reaction_scope_missing",
            "reason": (
                "The target is open-shell, but no independent backend searched radical "
                "encounters or self-reactions under the retained conditions."
            ),
            "required_next_evidence": "independent radical/bimolecular discovery",
        })
    for result in assembly.get("backend_results") or []:
        if result.get("backend") == "rmg":
            continue
        provenance = result.get("provenance") or {}
        for assessment in provenance.get("reactive_partner_assessments") or []:
            if not isinstance(assessment, Mapping):
                continue
            if assessment.get("status") == "completed_declared_search":
                continue
            gaps.append({
                "code": "independent_reactive_partner_scope_unresolved",
                "reactive_partner": assessment.get("label"),
                "partner_status": assessment.get("status"),
                "reason": str(
                    assessment.get("reason")
                    or "The independent backend did not complete its declared reactive-partner search."
                ),
                "required_next_evidence": (
                    "a declared partner electronic state and completed bounded pairwise search"
                ),
            })
    return gaps


def _verification_ids(summary: Mapping[str, Any] | None) -> tuple[set[str], set[str]]:
    if not summary:
        return set(), set()
    verified = {
        str(item) for item in summary.get("verified_route_ids") or [] if item
    }
    deferred = {
        str(item.get("route_id"))
        for item in summary.get("deferred_unverified_routes") or []
        if isinstance(item, Mapping) and item.get("route_id")
    }
    return verified, deferred


def _candidate_influence(
    assembly: Mapping[str, Any],
    verification_summary: Mapping[str, Any] | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    verified, deferred = _verification_ids(verification_summary)
    verification_status = str((verification_summary or {}).get("status") or "not_run")
    queue_ids = {
        str(item.get("route_id"))
        for item in assembly.get("graph_only_verification_queue") or []
        if item.get("route_id")
    }
    assessments: list[dict[str, Any]] = []
    blockers: list[dict[str, Any]] = []
    candidates = (assembly.get("merged_result") or {}).get("candidates") or []
    for candidate in candidates:
        candidate_id = str(candidate.get("candidate_id") or "")
        sources = sorted({str(item) for item in candidate.get("discovery_sources") or []})
        screening = str(candidate.get("screening_status") or "candidate_unverified")
        graph_only = "graph_xtb" in sources and "rmg" not in sources
        if screening == "xtb_screened_out":
            influence = "excluded_by_declared_endpoint_energy_screen"
            reason = "Candidate exceeded the declared xTB endpoint-energy window; no barrier or absence claim is implied."
        elif graph_only and candidate_id in deferred:
            influence = "unresolved_candidate_kinetics"
            reason = "Independent candidate could not obtain verified path/rate evidence."
        elif graph_only and candidate_id in verified and verification_status in _OBSERVABLE_TERMINAL_STATUSES:
            influence = "verified_and_propagated"
            reason = "Candidate received verified rate evidence and was included in observable-robust repropagation."
        elif graph_only and candidate_id in verified:
            influence = "verified_rate_network_influence_unresolved"
            reason = "A rate replacement exists, but observable convergence has not been established."
        elif graph_only:
            influence = "potentially_influential_unknown_rate"
            reason = (
                "Independent target-connectivity candidate has no verified rate; xTB "
                "reaction energy cannot establish negligible observable influence."
            )
        else:
            influence = "represented_in_retained_kinetic_model"
            reason = (
                "Candidate is represented by RMG kinetic evidence; its observable influence "
                "is governed by condition-matched flux, sensitivity, and verification gates."
            )
        record = {
            "candidate_id": candidate_id,
            "discovery_sources": sources,
            "reaction_classes": sorted(_candidate_reaction_classes(candidate)),
            "screening_status": screening,
            "influence_state": influence,
            "queued_for_verification": candidate_id in queue_ids,
            "verified_route_evidence_retained": candidate_id in verified,
            "reason": reason,
        }
        assessments.append(record)
        if influence in {
            "unresolved_candidate_kinetics",
            "verified_rate_network_influence_unresolved",
            "potentially_influential_unknown_rate",
        }:
            blockers.append({
                "code": "unresolved_candidate_kinetics",
                "candidate_id": candidate_id,
                "reason": reason,
                "required_next_evidence": (
                    "condition-applicable path verification, rate evidence, repaired-network "
                    "repropagation, and observable robustness"
                ),
            })
    return assessments, blockers


def assess_discovery_coverage(
    assembly: Any,
    policy: CoveragePolicy,
    *,
    observable_context: Mapping[str, Any] | None = None,
    verification_summary: Mapping[str, Any] | None = None,
    applicability_gaps: Iterable[Mapping[str, Any]] = (),
    search_expansion_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a deterministic, fail-closed mechanism-coverage assessment."""
    retained = _assembly_dict(assembly)
    policy_artifact = policy.as_dict()
    blockers: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    backend_assessments = retained.get("backend_assessments") or {}
    retained_backends = set(backend_assessments)
    for backend in policy.required_backends:
        if backend not in retained_backends:
            blockers.append({
                "code": "required_backend_missing", "backend": backend,
                "reason": "The coverage policy requires a backend with no retained assessment.",
            })
            continue
        backend_status = str((backend_assessments.get(backend) or {}).get("status") or "")
        if backend_status in {"backend_failed", "backend_unavailable"}:
            blockers.append({
                "code": "unresolved_backend_failure", "backend": backend,
                "backend_status": backend_status,
                "reason": "A required discovery backend failed or was unavailable.",
            })
        elif backend_status == "discovery_truncated":
            blockers.append({
                "code": "unresolved_search_truncation", "backend": backend,
                "backend_status": backend_status,
                "reason": "A required bounded discovery search exhausted its declared budget.",
            })
        elif backend_status in {"incomplete", "completed_with_unresolved_candidates"}:
            blockers.append({
                "code": "unresolved_backend_evidence", "backend": backend,
                "backend_status": backend_status,
                "reason": "A required backend retained unresolved discovery evidence.",
            })

    disagreement = analyze_discovery_disagreement(retained)
    completed_required = {
        backend for backend in policy.required_backends
        if str((backend_assessments.get(backend) or {}).get("status") or "") in {
            "completed_declared_search", "completed_no_candidates_in_declared_search",
        }
    }
    if (
        policy.require_independent_discovery
        and len(completed_required) < policy.minimum_independent_backend_count
    ):
        blockers.append({
            "code": "insufficient_independent_discovery",
            "completed_required_backends": sorted(completed_required),
            "required_backend_count": policy.minimum_independent_backend_count,
            "reason": (
                "The observable coverage policy requires multiple completed, "
                "methodologically distinct discovery sources."
            ),
        })
    if (
        policy.block_empty_shared_search_without_calibration
        and disagreement.get("empty_shared_search")
        and policy.calibration_status != "blind_spot_calibrated"
    ):
        blockers.append({
            "code": "empty_search_shared_blind_spot_unresolved",
            "completed_backends": disagreement.get("completed_backends") or [],
            "reason": (
                "Multiple bounded backends returned no candidate, but the coverage policy "
                "is not calibrated to rule out their shared blind spots."
            ),
            "required_next_evidence": (
                "bounded search expansion, a third complementary discovery source, or "
                "preregistered blind-spot calibration applicable to this chemical scope"
            ),
        })
    for row in disagreement.get("candidate_disagreements") or []:
        if row.get("status") == "source_disagreement":
            warnings.append({
                "code": "discovery_source_disagreement",
                "candidate_id": row.get("candidate_id"),
                "found_by": row.get("found_by") or [],
                "missing_eligible_backends": row.get("missing_eligible_backends") or [],
                "reason": (
                    "A completed eligible backend did not propose this candidate. The "
                    "candidate must be resolved on its evidence; no backend is presumed correct."
                ),
            })

    independent_scope_gaps = _independent_scope_gaps(retained)
    for gap in independent_scope_gaps:
        blockers.append(dict(gap))

    accounting = retained.get("candidate_accounting") or {}
    for identifier in accounting.get("unresolved_candidate_or_assertion_ids") or []:
        blockers.append({
            "code": "unresolved_candidate_or_assertion",
            "candidate_or_assertion_id": str(identifier),
            "reason": "Discovery evidence could not be adapted or screened to a resolved candidate state.",
        })
    independently_retained_ids = {
        str(identifier)
        for identifier in accounting.get("source_agreement_candidate_ids") or []
    }
    for identifier in accounting.get("excluded_by_declared_screen_ids") or []:
        identifier = str(identifier)
        independently_retained = identifier in independently_retained_ids
        warnings.append({
            "code": "candidate_excluded_by_declared_screen",
            "candidate_id": identifier,
            "independently_retained": independently_retained,
            "reason": (
                "The xTB exclusion is conditional on the declared endpoint screen and "
                "is not proof of reaction absence."
                + (
                    " The same native reaction was retained by an independent backend, "
                    "so this disagreement did not remove it from the candidate network."
                    if independently_retained else ""
                )
            ),
        })
        if not independently_retained:
            blockers.append({
                "code": "uncalibrated_endpoint_screen_exclusion",
                "candidate_id": identifier,
                "reason": (
                    "The candidate was excluded by an xTB endpoint reaction-energy "
                    "threshold. Endpoint energy is not a barrier and cannot establish "
                    "negligible network influence. STORCA has no validated, hash-bound "
                    "calibration contract that can admit this exclusion."
                ),
                "required_next_evidence": (
                    "retain and verify the candidate or expand/remove the endpoint-energy "
                    "screen; a future bypass requires a retained benchmark calibration "
                    "bound to the exact screen and decision scope"
                ),
            })

    class_matrix = _reaction_class_matrix(retained)
    searched_classes = {
        str(row["reaction_class"]) for row in class_matrix if row.get("searched_by")
    }
    for reaction_class in policy.required_reaction_classes:
        if reaction_class not in searched_classes:
            blockers.append({
                "code": "required_reaction_class_unsearched",
                "reaction_class": reaction_class,
                "reason": "The lane policy requires a reaction class that no retained backend explicitly searched.",
            })

    normalized_gaps = []
    for index, item in enumerate(applicability_gaps):
        gap = _jsonable(dict(item))
        gap.setdefault("gap_id", f"applicability-gap-{index + 1}")
        gap.setdefault("blocking", True)
        normalized_gaps.append(gap)
        if gap["blocking"]:
            blockers.append({
                "code": "outside_method_applicability",
                "gap_id": gap["gap_id"],
                "reason": str(gap.get("reason") or "A declared applicability gap is unresolved."),
                "required_method": gap.get("required_method"),
            })

    candidate_influence, candidate_blockers = _candidate_influence(
        retained, verification_summary,
    )
    blockers.extend(candidate_blockers)

    expansion = _jsonable(search_expansion_evidence) if search_expansion_evidence else None
    if policy.require_search_expansion and expansion is None:
        blockers.append({
            "code": "search_expansion_evidence_missing",
            "reason": "The coverage policy requires a bounded search-expansion comparison.",
        })
    elif expansion is not None and expansion.get("status") != "stable_within_declared_expansion":
        blockers.append({
            "code": "unresolved_budget_expansion",
            "reason": "The candidate set changed or could not be compared under the declared expansion.",
            "expansion_status": expansion.get("status"),
            "new_candidate_ids": expansion.get("new_candidate_ids") or [],
        })

    verification_status = str((verification_summary or {}).get("status") or "not_run")
    assembly_condition = (retained.get("contract") or {}).get("condition_contract") or {}
    verification_condition = (verification_summary or {}).get("condition_contract") or {}
    if verification_summary and (
        assembly_condition.get("contract_id")
        and verification_condition.get("contract_id")
        and assembly_condition.get("contract_id") != verification_condition.get("contract_id")
    ):
        blockers.append({
            "code": "condition_contract_mismatch",
            "reason": "Verification and discovery coverage evidence use different condition contracts.",
        })
    if verification_summary and verification_status not in _OBSERVABLE_TERMINAL_STATUSES:
        blockers.append({
            "code": "observable_convergence_unresolved",
            "verification_status": verification_status,
            "reason": "Path/rate/network evidence has not established a robust terminal observable.",
        })

    # Deduplicate deterministic blocker records without hiding distinct targets.
    unique_blockers: list[dict[str, Any]] = []
    seen: set[str] = set()
    for blocker in blockers:
        key = json.dumps(_jsonable(blocker), sort_keys=True, separators=(",", ":"))
        if key not in seen:
            unique_blockers.append(blocker)
            seen.add(key)
    blockers = unique_blockers

    codes = {str(item.get("code")) for item in blockers}
    if "outside_method_applicability" in codes:
        status = "outside_method_applicability"
    elif "unresolved_backend_failure" in codes:
        status = "unresolved_backend_failure"
    elif "unresolved_search_truncation" in codes:
        status = "unresolved_search_truncation"
    elif codes & {
        "insufficient_independent_discovery",
        "required_backend_missing",
        "empty_search_shared_blind_spot_unresolved",
        "independent_bimolecular_discovery_scope_missing",
        "independent_radical_reaction_scope_missing",
        "independent_reactive_partner_scope_unresolved",
        "uncalibrated_endpoint_screen_exclusion",
    }:
        status = "unresolved_mechanism_coverage_risk"
    elif codes & {
        "unresolved_candidate_kinetics", "unresolved_candidate_or_assertion",
        "unresolved_backend_evidence",
    }:
        status = "unresolved_candidate_kinetics"
    elif codes & {"unresolved_budget_expansion", "search_expansion_evidence_missing"}:
        status = "unresolved_budget_expansion"
    elif blockers:
        status = "unresolved_observable_convergence"
    elif (
        verification_status in _OBSERVABLE_TERMINAL_STATUSES
        and policy.allow_observable_claim_within_declared_scope
    ):
        status = "adequate_within_declared_scope"
    else:
        status = "consistent_within_declared_searches"

    supports_observable = status == "adequate_within_declared_scope"
    if supports_observable:
        claim = (
            "The requested observable is supported under the exact retained conditions "
            "within the declared discovery and method scope."
        )
    elif blockers:
        claim = "No condition-specific observable claim is defensible until the retained blockers are resolved."
    else:
        claim = (
            "The declared searches are internally consistent, but terminal observable "
            "verification has not been supplied."
        )
    if status in {
        "outside_method_applicability", "unresolved_backend_failure",
        "unresolved_search_truncation", "unresolved_mechanism_coverage_risk",
    }:
        risk_level = "critical"
    elif blockers:
        risk_level = "high"
    elif warnings:
        risk_level = "guarded"
    else:
        risk_level = "controlled_within_declared_scope"
    uncertainty_register = {
        "mechanism_coverage": {
            "status": "controlled_within_declared_scope" if supports_observable else "unresolved",
            "risk_level": risk_level,
            "evidence_id": disagreement.get("evidence_id"),
            "blocking_codes": sorted(codes),
        },
        **{
            category: {
                "status": "not_assessed_by_discovery_coverage",
                "reason": "This uncertainty domain requires its own retained evidence and policy.",
            }
            for category in (
                "electronic_structure", "conformer_sampling", "transition_state_search",
                "environment_representation", "kinetic_theory", "tunneling",
                "pressure_dependence", "electronic_state", "network_propagation",
                "numerical_error",
            )
        },
    }
    payload = {
        "schema_version": DISCOVERY_COVERAGE_SCHEMA_VERSION,
        "kind": "storca_discovery_coverage_assessment",
        "status": status,
        "supports_observable_claim": supports_observable,
        "observable_context": _jsonable(dict(observable_context or {})),
        "policy": policy_artifact,
        "assembly_contract_id": (retained.get("contract") or {}).get("contract_id"),
        "condition_contract_id": assembly_condition.get("contract_id"),
        "assembly_status": retained.get("status"),
        "verification_status": verification_status,
        "reaction_class_coverage": class_matrix,
        "discovery_disagreement": disagreement,
        "independent_scope_gaps": independent_scope_gaps,
        "candidate_influence": candidate_influence,
        "applicability_gaps": normalized_gaps,
        "search_expansion_evidence": expansion,
        "blocking_findings": blockers,
        "warnings": warnings,
        "scientific_claim": claim,
        "mechanism_completeness_claimed": False,
        "calibration_status": policy.calibration_status,
        "endpoint_screen_calibration_status": "not_calibrated_no_bypass_supported",
        "risk_summary": {
            "level": risk_level,
            "blocking_finding_count": len(blockers),
            "warning_count": len(warnings),
            "unqualified_mechanism_completeness_claim_allowed": False,
            "observable_claim_qualification": (
                "conditional_on_declared_search_and_method_scope"
                if supports_observable else "observable_withheld"
            ),
        },
        "uncertainty_register": uncertainty_register,
        "interpretation": (
            "This categorical decision is observable-, condition-, environment-, and "
            "declared-search-scoped. It is not a probability that the true mechanism is complete."
        ),
    }
    return {**payload, "assessment_id": _content_id(payload)}


def write_discovery_coverage(path: Path, assessment: Mapping[str, Any]) -> Path:
    """Retain an immutable content-addressed revision plus the latest pointer."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    retained = _jsonable(dict(assessment))
    assessment_id = str(retained.get("assessment_id") or "")
    unhashed = {
        key: value for key, value in retained.items()
        if key not in {"assessment_id", "artifact", "revision_artifact"}
    }
    expected_id = _content_id(unhashed)
    if assessment_id != expected_id:
        raise ValueError("Discovery-coverage assessment_id does not match its content")
    revision_dir = path.parent / "coverage-assessments"
    revision_dir.mkdir(parents=True, exist_ok=True)
    revision_path = revision_dir / f"{assessment_id.split(':', 1)[-1]}.json"
    if revision_path.is_file():
        if json.loads(revision_path.read_text()) != retained:
            raise ValueError("Retained discovery-coverage revision was mutated")
    else:
        revision_temporary = revision_path.with_suffix(".json.tmp")
        revision_temporary.write_text(
            json.dumps(retained, indent=2, sort_keys=True, default=str) + "\n"
        )
        revision_temporary.replace(revision_path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(retained, indent=2, sort_keys=True, default=str) + "\n"
    )
    temporary.replace(path)
    return path


def score_coverage_benchmark(cases: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    """Score held-out coverage decisions with false reassurance as primary risk.

    Each case must declare whether an influential reference pathway is missing
    and supply the corresponding assessment.  This is a harness, not a bundled
    claim that the repository currently contains a scientifically curated
    blind-spot dataset.
    """
    records = []
    for item in cases:
        case_id = str(item.get("case_id") or "")
        if not case_id:
            raise ValueError("Coverage benchmark case needs case_id")
        if "influential_reference_path_missing" not in item:
            raise ValueError(
                f"Coverage benchmark case {case_id} needs influential_reference_path_missing"
            )
        assessment = item.get("assessment") or {}
        supports = bool(assessment.get("supports_observable_claim"))
        missing = bool(item["influential_reference_path_missing"])
        records.append({
            "case_id": case_id,
            "partition": str(item.get("partition") or "unspecified"),
            "stratum": str(item.get("stratum") or "unspecified"),
            "influential_reference_path_missing": missing,
            "supports_observable_claim": supports,
            "false_reassurance": missing and supports,
            "correct_abstention": missing and not supports,
            "false_abstention": not missing and not supports,
            "correct_support": not missing and supports,
        })
    missing_cases = [item for item in records if item["influential_reference_path_missing"]]
    covered_cases = [item for item in records if not item["influential_reference_path_missing"]]
    false_reassurance = sum(item["false_reassurance"] for item in records)
    correct_abstention = sum(item["correct_abstention"] for item in records)
    false_abstention = sum(item["false_abstention"] for item in records)
    correct_support = sum(item["correct_support"] for item in records)
    strata = {}
    for stratum in sorted({item["stratum"] for item in records}):
        subset = [item for item in records if item["stratum"] == stratum]
        subset_missing = [
            item for item in subset if item["influential_reference_path_missing"]
        ]
        subset_covered = [
            item for item in subset if not item["influential_reference_path_missing"]
        ]
        strata[stratum] = {
            "case_count": len(subset),
            "missing_influential_path_case_count": len(subset_missing),
            "false_reassurance_count": sum(
                item["false_reassurance"] for item in subset
            ),
            "false_reassurance_rate": (
                sum(item["false_reassurance"] for item in subset_missing)
                / len(subset_missing) if subset_missing else None
            ),
            "false_abstention_count": sum(item["false_abstention"] for item in subset),
            "false_abstention_rate": (
                sum(item["false_abstention"] for item in subset_covered)
                / len(subset_covered) if subset_covered else None
            ),
        }
    payload = {
        "schema_version": DISCOVERY_COVERAGE_SCHEMA_VERSION,
        "kind": "storca_discovery_coverage_benchmark_score",
        "case_count": len(records),
        "missing_influential_path_case_count": len(missing_cases),
        "false_reassurance_count": false_reassurance,
        "false_reassurance_rate": (
            false_reassurance / len(missing_cases) if missing_cases else None
        ),
        "correct_abstention_count": correct_abstention,
        "correct_abstention_rate": (
            correct_abstention / len(missing_cases) if missing_cases else None
        ),
        "reference_covered_case_count": len(covered_cases),
        "false_abstention_count": false_abstention,
        "false_abstention_rate": (
            false_abstention / len(covered_cases) if covered_cases else None
        ),
        "correct_support_count": correct_support,
        "correct_support_rate": (
            correct_support / len(covered_cases) if covered_cases else None
        ),
        "strata": strata,
        "cases": records,
        "interpretation": (
            "Scientific use requires preregistered, provenance-backed held-out cases; "
            "synthetic unit fixtures validate software behavior only."
        ),
    }
    return {**payload, "score_id": _content_id(payload)}


__all__ = [
    "DISCOVERY_COVERAGE_SCHEMA_VERSION", "CoveragePolicy",
    "analyze_discovery_disagreement", "assess_discovery_coverage",
    "compare_discovery_expansion",
    "score_coverage_benchmark", "write_discovery_coverage",
]
