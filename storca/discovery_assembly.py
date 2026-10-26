"""Immutable orchestration for bounded multi-source reaction discovery.

The assembly layer coordinates discovery evidence; it does not verify a path,
construct a rate, or propagate a kinetic model.  RMG remains the condition-
specific kinetic-model producer, while graph/xTB remains a reusable bounded
connectivity and endpoint-screening backend.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Mapping

from .reaction_candidates import ReactionCandidate, canonical_smiles
from .reaction_discovery import (
    DiscoveryResult,
    ReactionDiscoveryBackend,
    RmgReactionDiscoveryBackend,
    graph_only_candidates_for_verification,
    merge_discovery_results,
)


DISCOVERY_ASSEMBLY_SCHEMA_VERSION = 1
RmgEvidenceProducer = Callable[[], Mapping[str, Any]]


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


def _canonical_hash(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        _jsonable(value), sort_keys=True, separators=(",", ":"),
        allow_nan=False, default=str,
    )
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


def _append_event(root: Path, value: Mapping[str, Any]) -> Path:
    event_dir = Path(root) / "assembly-events"
    event_dir.mkdir(parents=True, exist_ok=True)
    index = len(list(event_dir.glob("*.json")))
    return _write_json(event_dir / f"{index:06d}.json", {
        "schema_version": DISCOVERY_ASSEMBLY_SCHEMA_VERSION,
        "kind": "storca_discovery_assembly_event",
        **dict(value),
    })


@dataclass(frozen=True)
class DiscoveryAssemblyContract:
    """Path-independent identity for one condition-bound discovery assembly."""

    system: Mapping[str, Any]
    condition_contract: Mapping[str, Any]
    required_backends: tuple[str, ...]
    rmg_contract: Mapping[str, Any]
    graph_xtb_budget: Mapping[str, Any] | None
    precomputed_graph_xtb_sha256: str | None

    def __post_init__(self) -> None:
        required = tuple(sorted({str(item) for item in self.required_backends}))
        if not required or not set(required) <= {"rmg", "graph_xtb"}:
            raise ValueError("Discovery assembly requires rmg and/or graph_xtb backends")
        if "rmg" not in required:
            raise ValueError("STORCA discovery assembly currently requires RMG")
        target = str(self.system.get("smiles") or self.system.get("target_smiles") or "")
        if not target:
            raise ValueError("Discovery assembly system needs a target SMILES")
        if not self.condition_contract:
            raise ValueError("Discovery assembly needs a complete condition contract")
        if "graph_xtb" in required and self.graph_xtb_budget is None:
            raise ValueError("Requested graph_xtb discovery needs a resolved budget")
        object.__setattr__(self, "required_backends", required)

    def scientific_dict(self) -> dict[str, Any]:
        target = str(self.system.get("smiles") or self.system.get("target_smiles"))
        reactive_partners = []
        for item in self.system.get("reactive_partners") or ():
            if not isinstance(item, Mapping):
                raise ValueError("Discovery assembly reactive partners must be mappings")
            partner_smiles = str(item.get("smiles") or "")
            if not partner_smiles:
                raise ValueError("Discovery assembly reactive partner needs SMILES")
            reactive_partners.append({
                "label": str(item.get("label") or canonical_smiles(partner_smiles)),
                "canonical_smiles": canonical_smiles(partner_smiles),
                "charge": (
                    int(item["charge"]) if item.get("charge") is not None else None
                ),
                "multiplicity": (
                    int(item["multiplicity"])
                    if item.get("multiplicity") is not None else None
                ),
                "reactive": bool(item.get("reactive", True)),
            })
        reactive_partners.sort(key=lambda value: (
            value["label"], value["canonical_smiles"],
            value["charge"] if value["charge"] is not None else 0,
            value["multiplicity"] if value["multiplicity"] is not None else 0,
        ))
        return {
            "schema_version": DISCOVERY_ASSEMBLY_SCHEMA_VERSION,
            "kind": "storca_discovery_assembly_contract",
            "system": {
                "canonical_smiles": canonical_smiles(target),
                "charge": int(self.system.get("charge", 0)),
                "multiplicity": int(self.system.get("multiplicity", 1)),
                "reactive_partners": reactive_partners,
                "environment": _jsonable(self.system.get("environment") or {}),
            },
            "condition_contract": _jsonable(self.condition_contract),
            "required_backends": list(self.required_backends),
            "rmg_contract": _jsonable(self.rmg_contract),
            "graph_xtb_budget": _jsonable(self.graph_xtb_budget),
            "precomputed_graph_xtb_sha256": self.precomputed_graph_xtb_sha256,
            "reuse_policy": {
                "rmg": "condition_specific_not_reusable_across_conditions",
                "graph_xtb": "bounded_screen_reusable_after_exact_target_state_environment_and_budget_match",
            },
        }

    def as_dict(self) -> dict[str, Any]:
        value = self.scientific_dict()
        return {**value, "contract_id": _canonical_hash(value)}


def _restore_discovery(value: DiscoveryResult | Mapping[str, Any]) -> DiscoveryResult:
    if isinstance(value, DiscoveryResult):
        return value
    if isinstance(value, Mapping):
        return DiscoveryResult.from_dict(value)
    raise TypeError("Precomputed discovery evidence must be a DiscoveryResult or mapping")


def _validate_precomputed_graph(
    result: DiscoveryResult,
    *,
    system: Mapping[str, Any],
    expected_budget: Mapping[str, Any],
) -> None:
    if result.backend != "graph_xtb":
        raise ValueError("Precomputed discovery artifact must have backend='graph_xtb'")
    provenance = result.provenance
    target = str(system.get("smiles") or system.get("target_smiles") or "")
    if canonical_smiles(str(provenance.get("target_smiles") or "")) != canonical_smiles(target):
        raise ValueError("Precomputed graph/xTB discovery target does not match")
    if int(provenance.get("target_charge", system.get("charge", 0))) != int(
        system.get("charge", 0)
    ):
        raise ValueError("Precomputed graph/xTB discovery charge does not match")
    if int(provenance.get("target_multiplicity", system.get("multiplicity", 1))) != int(
        system.get("multiplicity", 1)
    ):
        raise ValueError("Precomputed graph/xTB discovery multiplicity does not match")
    if provenance.get("search_scope") != "single_target_unimolecular_graph_edit":
        raise ValueError("Precomputed graph/xTB discovery search scope is not reusable")
    if _jsonable(provenance.get("budget") or {}) != _jsonable(expected_budget):
        raise ValueError("Precomputed graph/xTB discovery budget does not match")
    from .graph_xtb_discovery import _screening_environment
    expected_environment = _screening_environment(system.get("environment") or {})
    observed_environment = provenance.get("screening_environment") or {"phase": "gas"}
    if observed_environment != expected_environment:
        raise ValueError("Precomputed graph/xTB screening environment does not match")


def _backend_assessment(result: DiscoveryResult) -> dict[str, Any]:
    unresolved_ids: list[str] = []
    excluded_ids: list[str] = []
    if result.backend == "rmg":
        unresolved_records = [
            record for record in result.records
            if record.get("status") == "not_generated"
        ]
        unresolved_ids = [
            str(record.get("route_id") or f"rmg-route-{record.get('route_index')}")
            for record in unresolved_records
        ]
    elif result.backend == "graph_xtb":
        unresolved_ids = [
            item.candidate_id for item in result.candidates
            if item.screening_status in {"xtb_failed", "search_truncated"}
        ]
        excluded_ids = [
            item.candidate_id for item in result.candidates
            if item.screening_status == "xtb_screened_out"
        ]

    if result.status == "discovery_truncated" or any(
        item.screening_status == "search_truncated" for item in result.candidates
    ):
        status = "discovery_truncated"
    elif result.status == "backend_unavailable":
        status = "backend_unavailable"
    elif result.status == "backend_failed":
        status = "backend_failed"
    elif result.status != "completed":
        status = "incomplete"
    elif unresolved_ids:
        status = "completed_with_unresolved_candidates"
    elif not result.candidates:
        status = "completed_no_candidates_in_declared_search"
    else:
        status = "completed_declared_search"
    return {
        "backend": result.backend,
        "status": status,
        "backend_reported_status": result.status,
        "candidate_count": len(result.candidates),
        "unresolved_candidate_or_assertion_ids": sorted(unresolved_ids),
        "excluded_by_declared_screen_ids": sorted(excluded_ids),
        "interpretation": (
            "Completion describes only the retained bounded backend search; it is not proof "
            "that no omitted chemistry exists."
        ),
    }


def _assembly_status(assessments: Mapping[str, Mapping[str, Any]]) -> str:
    statuses = {str(item.get("status")) for item in assessments.values()}
    if "backend_failed" in statuses:
        return "failed_closed"
    if "backend_unavailable" in statuses:
        return "backend_unavailable"
    if "discovery_truncated" in statuses:
        return "discovery_truncated"
    if statuses & {"incomplete", "completed_with_unresolved_candidates"}:
        return "completed_with_unresolved_candidates"
    return "completed_declared_search"


def _candidate_accounting(
    merged: DiscoveryResult,
    assessments: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    agreement = []
    rmg_only = []
    graph_only = []
    for candidate in merged.candidates:
        sources = set(candidate.discovery_sources)
        if {"rmg", "graph_xtb"} <= sources:
            agreement.append(candidate.candidate_id)
        elif "rmg" in sources:
            rmg_only.append(candidate.candidate_id)
        elif "graph_xtb" in sources:
            graph_only.append(candidate.candidate_id)
    unresolved = sorted({
        identifier
        for assessment in assessments.values()
        for identifier in assessment.get("unresolved_candidate_or_assertion_ids") or []
    })
    excluded = sorted({
        identifier
        for assessment in assessments.values()
        for identifier in assessment.get("excluded_by_declared_screen_ids") or []
    })
    return {
        "merged_candidate_count": len(merged.candidates),
        "source_agreement_candidate_ids": sorted(agreement),
        "rmg_only_candidate_ids": sorted(rmg_only),
        "graph_xtb_only_candidate_ids": sorted(graph_only),
        "unresolved_candidate_or_assertion_ids": unresolved,
        "excluded_by_declared_screen_ids": excluded,
        "mechanism_completeness_claimed": False,
    }


@dataclass(frozen=True)
class DiscoveryAssemblyOutcome:
    status: str
    contract: Mapping[str, Any]
    rmg_evidence: Mapping[str, Any]
    backend_results: tuple[DiscoveryResult, ...]
    merged_result: DiscoveryResult
    backend_assessments: Mapping[str, Mapping[str, Any]]
    candidate_accounting: Mapping[str, Any]
    graph_only_verification_queue: tuple[Mapping[str, Any], ...]
    reuse: Mapping[str, Any] | None
    artifacts: Mapping[str, Any]
    resumed: bool = False
    earlier_rmg_verification_queue: tuple[Mapping[str, Any], ...] = ()
    earlier_rmg_queue_assessment: Mapping[str, Any] | None = None

    @property
    def supplemental_verification_queue(self) -> tuple[Mapping[str, Any], ...]:
        """Share the existing bounded engine queue across discovery sources."""
        return (*self.graph_only_verification_queue, *self.earlier_rmg_verification_queue)

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": DISCOVERY_ASSEMBLY_SCHEMA_VERSION,
            "kind": "storca_discovery_assembly_outcome",
            "status": self.status,
            "terminal": True,
            "contract": dict(self.contract),
            "rmg_evidence": _jsonable(self.rmg_evidence),
            "backend_results": [item.as_dict() for item in self.backend_results],
            "merged_result": self.merged_result.as_dict(),
            "backend_assessments": _jsonable(self.backend_assessments),
            "candidate_accounting": _jsonable(self.candidate_accounting),
            "graph_only_verification_queue": [
                _jsonable(item) for item in self.graph_only_verification_queue
            ],
            "reuse": _jsonable(self.reuse),
            "artifacts": _jsonable(self.artifacts),
            "resumed": self.resumed,
            "earlier_rmg_verification_queue": _jsonable(self.earlier_rmg_verification_queue),
            "earlier_rmg_queue_assessment": _jsonable(self.earlier_rmg_queue_assessment),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DiscoveryAssemblyOutcome":
        return cls(
            status=str(value["status"]),
            contract=dict(value["contract"]),
            rmg_evidence=dict(value.get("rmg_evidence") or {}),
            backend_results=tuple(
                DiscoveryResult.from_dict(item)
                for item in value.get("backend_results") or ()
            ),
            merged_result=DiscoveryResult.from_dict(value["merged_result"]),
            backend_assessments=dict(value.get("backend_assessments") or {}),
            candidate_accounting=dict(value.get("candidate_accounting") or {}),
            graph_only_verification_queue=tuple(
                dict(item) for item in value.get("graph_only_verification_queue") or ()
            ),
            reuse=dict(value["reuse"]) if value.get("reuse") is not None else None,
            artifacts=dict(value.get("artifacts") or {}),
            resumed=bool(value.get("resumed", False)),
            earlier_rmg_verification_queue=tuple(dict(item) for item in value.get("earlier_rmg_verification_queue") or ()),
            earlier_rmg_queue_assessment=dict(value.get("earlier_rmg_queue_assessment") or {}),
        )

    def reaction_discovery_artifact(self) -> dict[str, Any]:
        """Return the backward-compatible merged artifact plus assembly evidence."""
        return {
            **self.merged_result.as_dict(),
            "assembly_status": self.status,
            "assembly_contract_id": self.contract.get("contract_id"),
            "backend_assessments": _jsonable(self.backend_assessments),
            "candidate_accounting": _jsonable(self.candidate_accounting),
            "graph_only_verification_queue": [
                _jsonable(item) for item in self.graph_only_verification_queue
            ],
            "reuse": _jsonable(self.reuse),
            "artifacts": _jsonable(self.artifacts),
            "resumed": self.resumed,
            "earlier_rmg_verification_queue": _jsonable(self.earlier_rmg_verification_queue),
            "earlier_rmg_queue_assessment": _jsonable(self.earlier_rmg_queue_assessment),
        }


class DiscoveryAssemblyExecutor:
    """Run or reuse requested discovery backends under one immutable contract."""

    def execute(
        self,
        *,
        system: Mapping[str, Any],
        condition_contract: Mapping[str, Any],
        workdir: Path,
        rmg_producer: RmgEvidenceProducer,
        rmg_contract: Mapping[str, Any],
        graph_xtb_budget: Mapping[str, Any] | None = None,
        graph_backend: ReactionDiscoveryBackend | None = None,
        precomputed_graph_xtb: DiscoveryResult | Mapping[str, Any] | None = None,
        maximum_graph_candidates_for_verification: int = 0,
        maximum_earlier_rmg_candidates_for_verification: int = 0,
    ) -> DiscoveryAssemblyOutcome:
        from .rmg_partial_candidates import retained_earlier_rmg_verification_queue

        if type(maximum_earlier_rmg_candidates_for_verification) is not int or maximum_earlier_rmg_candidates_for_verification < 0:
            raise ValueError("Earlier-RMG verification ceiling must be a nonnegative integer")
        if maximum_earlier_rmg_candidates_for_verification:
            rmg_contract = {**rmg_contract, "earlier_rmg_supplemental_queue": {
                "policy_id": "source_bound_completed_final_earlier_only_v1",
                "maximum_candidates": maximum_earlier_rmg_candidates_for_verification,
                "historical_kinetic_credit_allowed": False,
            }}
        if graph_backend is not None and precomputed_graph_xtb is not None:
            raise ValueError("Choose graph backend execution or precomputed evidence, not both")
        graph_requested = graph_backend is not None or precomputed_graph_xtb is not None
        precomputed_result = (
            _restore_discovery(precomputed_graph_xtb)
            if precomputed_graph_xtb is not None else None
        )
        if precomputed_result is not None:
            _validate_precomputed_graph(
                precomputed_result, system=system,
                expected_budget=graph_xtb_budget or {},
            )
        precomputed_hash = (
            _canonical_hash(precomputed_result.as_dict())
            if precomputed_result is not None else None
        )
        contract = DiscoveryAssemblyContract(
            system=system,
            condition_contract=condition_contract,
            required_backends=("rmg", "graph_xtb") if graph_requested else ("rmg",),
            rmg_contract=rmg_contract,
            graph_xtb_budget=graph_xtb_budget if graph_requested else None,
            precomputed_graph_xtb_sha256=precomputed_hash,
        ).as_dict()
        root = Path(workdir)
        root.mkdir(parents=True, exist_ok=True)
        contract_path = root / "discovery-assembly-contract.json"
        result_path = root / "discovery-assembly-result.json"
        merged_path = root / "reaction-discovery.json"
        if contract_path.is_file():
            retained_contract = json.loads(contract_path.read_text())
            if retained_contract.get("contract_id") != contract["contract_id"]:
                raise ValueError(
                    "Discovery assembly directory contains a different immutable contract"
                )
        else:
            _write_json(contract_path, contract)
        if result_path.is_file():
            retained = json.loads(result_path.read_text())
            if (
                retained.get("terminal") is True
                and (retained.get("contract") or {}).get("contract_id")
                == contract["contract_id"]
            ):
                outcome = DiscoveryAssemblyOutcome.from_dict(retained)
                if maximum_earlier_rmg_candidates_for_verification:
                    current = retained_earlier_rmg_verification_queue(
                        dict(outcome.rmg_evidence), dict(condition_contract),
                        maximum_candidates=maximum_earlier_rmg_candidates_for_verification,
                    )
                    if current["blocking_reasons"] or current["routes"] != list(outcome.earlier_rmg_verification_queue):
                        raise ValueError("Retained discovery assembly earlier-RMG source/view validation failed")
                return replace(outcome, resumed=True)

        _append_event(root, {
            "event": "assembly_started",
            "contract_id": contract["contract_id"],
            "required_backends": contract["required_backends"],
        })
        try:
            rmg_evidence = dict(rmg_producer())
        except Exception as error:
            rmg_evidence = {
                "status": "incomplete",
                "candidate_routes": [],
                "failure_reason": f"{type(error).__name__}: {error}",
            }
            rmg_result = DiscoveryResult(
                backend="rmg", status="backend_failed",
                provenance={
                    "failure_reason": rmg_evidence["failure_reason"],
                    "condition_contract_id": condition_contract.get("contract_id"),
                },
                warnings=(
                    "RMG discovery failed; no mechanism or observable completeness is implied.",
                ),
            )
        else:
            try:
                rmg_result = RmgReactionDiscoveryBackend(rmg_evidence).discover(
                    system, condition_contract, rmg_contract,
                    workdir=root / "rmg-adapter",
                )
            except Exception as error:
                rmg_result = DiscoveryResult(
                    backend="rmg", status="backend_failed",
                    provenance={
                        "failure_reason": f"{type(error).__name__}: {error}",
                        "condition_contract_id": condition_contract.get("contract_id"),
                        "raw_rmg_evidence_retained": True,
                    },
                    warnings=(
                        "RMG candidate adaptation failed; raw RMG evidence remains retained "
                        "but mechanism coverage is unresolved.",
                    ),
                )
        backend_results = [rmg_result]
        reuse: dict[str, Any] | None = None
        if graph_requested:
            if precomputed_result is not None:
                graph_result = precomputed_result
                reuse = {
                    "status": "reused_precomputed_target_only_discovery",
                    "source_condition_contract_id": graph_result.provenance.get(
                        "condition_contract_id"
                    ),
                    "reuse_condition_contract_id": condition_contract.get("contract_id"),
                    "search_scope": graph_result.provenance.get("search_scope"),
                    "interpretation": (
                        "Only target-connectivity and endpoint-screening evidence was reused; "
                        "RMG, verification, rates, and propagation remain condition-specific."
                    ),
                }
            else:
                try:
                    graph_result = graph_backend.discover(
                        system, condition_contract, graph_xtb_budget,
                        workdir=root / "graph-xtb",
                    )
                    if graph_result.backend != "graph_xtb":
                        raise ValueError(
                            "Configured graph/xTB backend returned a different backend identity"
                        )
                except Exception as error:
                    graph_result = DiscoveryResult(
                        backend="graph_xtb", status="backend_failed",
                        provenance={
                            "failure_reason": f"{type(error).__name__}: {error}",
                            "condition_contract_id": condition_contract.get(
                                "contract_id"
                            ),
                        },
                        warnings=(
                            "Graph/xTB discovery failed; no mechanism or observable "
                            "completeness is implied.",
                        ),
                    )
            backend_results.append(graph_result)

        assessments = {
            result.backend: _backend_assessment(result) for result in backend_results
        }
        earlier_queue = retained_earlier_rmg_verification_queue(
            rmg_evidence, dict(condition_contract),
            maximum_candidates=maximum_earlier_rmg_candidates_for_verification,
        )
        if earlier_queue["blocking_reasons"] or earlier_queue["queue_truncated"]:
            assessment = assessments["rmg"]
            assessment["status"] = "completed_with_unresolved_candidates"
            assessment["unresolved_candidate_or_assertion_ids"] = sorted(set([
                *assessment["unresolved_candidate_or_assertion_ids"],
                *earlier_queue["blocking_reasons"],
                *("earlier-rmg-unqueued:" + identifier for identifier in earlier_queue["unqueued_route_ids"]),
            ]))
        assembly_status = _assembly_status(assessments)
        merged = merge_discovery_results(tuple(backend_results))
        merged_status = (
            "completed" if assembly_status == "completed_declared_search"
            else "discovery_truncated" if assembly_status == "discovery_truncated"
            else "incomplete"
        )
        if merged.status != merged_status:
            merged = replace(merged, status=merged_status)
        accounting = _candidate_accounting(merged, assessments)
        target = str(system.get("smiles") or system.get("target_smiles"))
        queue = graph_only_candidates_for_verification(
            merged.candidates,
            parent_smiles=target,
            target_label=str(system.get("target_label") or "stability"),
            maximum_candidates=maximum_graph_candidates_for_verification,
        )
        for route in queue:
            route["discovery_condition_scope"] = (
                "target_only_connectivity_and_endpoint_screening"
            )
            route["verification_condition_contract_id"] = condition_contract.get(
                "contract_id"
            )
        _write_json(merged_path, merged.as_dict())
        outcome = DiscoveryAssemblyOutcome(
            status=assembly_status,
            contract=contract,
            rmg_evidence=rmg_evidence,
            backend_results=tuple(backend_results),
            merged_result=merged,
            backend_assessments=assessments,
            candidate_accounting=accounting,
            graph_only_verification_queue=tuple(queue),
            earlier_rmg_verification_queue=tuple(earlier_queue["routes"]),
            earlier_rmg_queue_assessment={key: value for key, value in earlier_queue.items() if key != "routes"},
            reuse=reuse,
            artifacts={
                "assembly_contract": str(contract_path),
                "assembly_result": str(result_path),
                "merged_discovery": str(merged_path),
                "graph_xtb_discovery": (
                    str(root / "graph-xtb" / "graph-xtb-discovery.json")
                    if graph_requested and precomputed_result is None else None
                ),
            },
        )
        _write_json(result_path, outcome.as_dict())
        _append_event(root, {
            "event": "assembly_retained",
            "contract_id": contract["contract_id"],
            "status": assembly_status,
            "candidate_accounting": accounting,
        })
        return outcome


__all__ = [
    "DISCOVERY_ASSEMBLY_SCHEMA_VERSION", "DiscoveryAssemblyContract",
    "DiscoveryAssemblyExecutor", "DiscoveryAssemblyOutcome",
]
