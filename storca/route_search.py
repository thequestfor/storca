"""Deterministic, persistent planning for an explicitly requested verified route.

This module neither runs chemistry nor accounts for native calculations.  The
caller supplies the existing execution, budget and final scientific validators.
The append-only search journal records *decisions*, not charged attempts.  A
missing outcome is never permission to repeat a potentially live calculation.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Mapping


class RouteSearchError(RuntimeError):
    """Retained search evidence or an adapter contract is inconsistent."""


DEFAULT_STRATEGIES = (
    "reanalyze_paid_evidence", "validate_reached_minima", "continue_verification", "continue_irc",
    "verify_ts", "verify_ts_seed", "resume_neb", "resume_band", "guided_scan", "neb",
    "validate_intermediate", "verify_segments", "prepare_seeds", "prepare_endpoints",
)
_PATH_KEYS = frozenset({"path", "directory", "workdir", "workspace", "folder", "basename"})
def _bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def _copy(value):
    try:
        return json.loads(_bytes(value))
    except (ValueError, TypeError) as error:
        raise RouteSearchError("Search evidence must be finite JSON data") from error


def _hash(value):
    return hashlib.sha256(_bytes(value)).hexdigest()


def _scientific(value):
    """Location-independent identity; every omitted path needs a content hash."""
    if isinstance(value, dict):
        locations = _PATH_KEYS.intersection(value)
        if locations and not any(key in value for key in ("sha256", "content_sha256", "input_sha256")):
            raise RouteSearchError("Scientific input paths require an adjacent content hash")
        return {key: _scientific(item) for key, item in value.items() if key not in _PATH_KEYS}
    if isinstance(value, list):
        return [_scientific(item) for item in value]
    if isinstance(value, str) and Path(value).is_absolute():
        raise RouteSearchError("Absolute artifact paths belong in content bindings or artifacts, not scientific identity")
    return value


def _cost(value):
    value = _copy(value)
    if not isinstance(value, dict) or any(type(item) is not int or item < 0 for item in value.values()):
        raise RouteSearchError("Category forecasts must be nonnegative integers")
    return value


@dataclass(frozen=True)
class RouteSearchPolicy:
    objective: str = "route_required"
    strategy_order: tuple = DEFAULT_STRATEGIES
    default_completion_cost: Mapping = field(default_factory=dict)
    extra_completion_cost: Mapping = field(default_factory=dict)
    validation_policy: Mapping = field(default_factory=dict)
    source_inventory: Mapping = field(default_factory=dict)
    version: str = "adaptive_verified_route_v1"

    def __post_init__(self):
        if self.objective != "route_required" or not self.version:
            raise RouteSearchError("Only an explicit route_required objective is supported")
        order = tuple(self.strategy_order)
        if not order or len(set(order)) != len(order) or any(not isinstance(item, str) or not item for item in order):
            raise RouteSearchError("An ordered, unique supported strategy set is required")
        object.__setattr__(self, "strategy_order", order)
        object.__setattr__(self, "default_completion_cost", _cost(self.default_completion_cost))
        object.__setattr__(self, "extra_completion_cost", _cost(self.extra_completion_cost))
        object.__setattr__(self, "validation_policy", _copy(self.validation_policy))
        object.__setattr__(self, "source_inventory", _copy(self.source_inventory))

    def as_dict(self):
        return _copy(self.__dict__)

    @classmethod
    def from_dict(cls, value):
        return cls(**value)


@dataclass(frozen=True)
class SearchAction:
    kind: str
    inputs: Mapping
    context: Mapping = field(default_factory=dict)
    parameters: Mapping = field(default_factory=dict)
    parent_evidence: tuple = ()
    cost_envelope: Mapping = field(default_factory=dict)
    completion_cost: Mapping | None = None
    preconditions: Mapping = field(default_factory=dict)
    rank: tuple = ()
    reason: str = ""
    artifacts: Mapping = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.kind, str) or not self.kind:
            raise RouteSearchError("Action kind is required")
        for name in ("inputs", "context", "parameters", "preconditions", "artifacts"):
            value = _copy(getattr(self, name))
            if not isinstance(value, dict):
                raise RouteSearchError(f"Action {name} must be a mapping")
            object.__setattr__(self, name, value)
        if not self.inputs or not _scientific(self.inputs):
            raise RouteSearchError("Action needs scientific inputs, not a workspace name")
        _scientific(self.context)
        _scientific(self.parameters)
        if any(type(value) is not bool for value in self.preconditions.values()):
            raise RouteSearchError("Action preconditions must be explicit booleans")
        if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) for value in self.rank):
            raise RouteSearchError("Action rank must contain finite numeric criteria")
        object.__setattr__(self, "rank", tuple(self.rank))
        object.__setattr__(self, "parent_evidence", tuple(_copy(self.parent_evidence)))
        object.__setattr__(self, "cost_envelope", _cost(self.cost_envelope))
        if self.completion_cost is not None:
            object.__setattr__(self, "completion_cost", _cost(self.completion_cost))

    def fingerprint(self, validation_policy=None, search_context=None):
        # Labels, parent explanations, paths and budget forecasts do not make a
        # failed physical experiment a new experiment.
        return _hash({"kind": self.kind, "inputs": _scientific(self.inputs),
                      "context": _scientific(self.context), "parameters": _scientific(self.parameters),
                      "search_context": _scientific(_copy(search_context or {})),
                      "validation_policy": _copy(validation_policy or {})})

    def as_dict(self):
        return _copy(self.__dict__)

    @classmethod
    def from_dict(cls, value):
        return value if isinstance(value, cls) else cls(**value)


@dataclass(frozen=True)
class ActionOutcome:
    execution: str = "completed"
    reader: str = "passed"
    physical: str = "unresolved"
    applicability: str = "supported"
    evidence: Mapping = field(default_factory=dict)
    native_reservation_ids: tuple = ()
    reason: str = ""

    def __post_init__(self):
        if self.execution not in {"completed", "not_required", "not_started", "failed", "blocked_resources", "ownership_unresolved", "cancelled"}:
            raise RouteSearchError("Unknown execution outcome")
        if self.reader not in {"passed", "failed", "not_run"}:
            raise RouteSearchError("Unknown reader outcome")
        if self.applicability not in {"supported", "unsupported", "unresolved"}:
            raise RouteSearchError("Unknown physical applicability outcome")
        if not isinstance(self.physical, str) or not self.physical:
            raise RouteSearchError("Physical outcome must remain explicit")
        object.__setattr__(self, "evidence", _copy(self.evidence))
        object.__setattr__(self, "native_reservation_ids", tuple(_copy(self.native_reservation_ids)))
        if self.execution == "not_started" and self.native_reservation_ids:
            raise RouteSearchError("A charged/reserved action cannot be labelled not_started")

    def as_dict(self):
        return _copy(self.__dict__)

    @classmethod
    def from_dict(cls, value):
        return value if isinstance(value, cls) else cls(**value)


@dataclass(frozen=True)
class RouteValidation:
    status: str = "unresolved"
    reached_minima_validated: bool = False
    all_required_segments_verified: bool = False
    evidence: Mapping = field(default_factory=dict)
    reason: str = ""

    def __post_init__(self):
        if self.status not in {"unresolved", "verified_route", "verified_multistep_route", "unsupported_physical_model"}:
            raise RouteSearchError("Unknown canonical route validation status")
        if type(self.reached_minima_validated) is not bool or type(self.all_required_segments_verified) is not bool:
            raise RouteSearchError("Final validation flags must be explicit booleans")
        object.__setattr__(self, "evidence", _copy(self.evidence))
        if self.status.startswith("verified_") and not (self.reached_minima_validated and self.all_required_segments_verified and self.evidence):
            raise RouteSearchError("Route success requires canonical segment and reached-minimum evidence")

    def as_dict(self):
        return _copy(self.__dict__)

    @classmethod
    def from_dict(cls, value):
        return value if isinstance(value, cls) else cls(**value)


@dataclass(frozen=True)
class SearchSnapshot:
    policy: RouteSearchPolicy
    context: Mapping
    completed: tuple
    pending: tuple
    frontier: tuple
    budget: Mapping
    validation: RouteValidation = field(default_factory=RouteValidation)
    status: str = "searching"
    exclusions: tuple = ()

    @property
    def outcomes(self):
        return self.completed

    @property
    def evidence(self):
        return tuple(item["outcome"]["evidence"] for item in self.completed)

    def as_dict(self):
        return {"status": self.status, "completed": _copy(self.completed), "pending": _copy(self.pending),
                "frontier": [item.as_dict() for item in self.frontier], "budget": _copy(self.budget),
                "validation": self.validation.as_dict(), "exclusions": _copy(self.exclusions),
                "rate_claim_allowed": False, "lifetime_claim_allowed": False}


def _budget_status(budget):
    if budget.get("pending") or budget.get("pending_reservation_ids"):
        return "ownership_unresolved"
    status = budget.get("status", "ready")
    if status in {"blocked_resources", "ownership_unresolved", "cancelled", "unsupported_physical_model"}:
        return status
    return None


def _eligibility(action, state):
    if action.kind not in state.policy.strategy_order:
        return "unsupported_strategy"
    if not all(action.preconditions.values()):
        return "unmet_preconditions"
    fingerprint = action.fingerprint(state.policy.validation_policy, state.context)
    if any(item["fingerprint"] == fingerprint and not (
        item["outcome"]["execution"] in {"not_started", "blocked_resources"}
        and not item["outcome"]["native_reservation_ids"]
    ) for item in state.completed):
        return "already_attempted_scientific_input"
    if any(item["fingerprint"] == fingerprint for item in state.pending):
        return "pending_ownership"
    reserve = action.completion_cost if action.completion_cost is not None else state.policy.default_completion_cost
    remaining = state.budget.get("remaining", {})
    extra = state.policy.extra_completion_cost
    for category in set(action.cost_envelope) | set(reserve) | set(extra):
        needed = action.cost_envelope.get(category, 0) + reserve.get(category, 0) + extra.get(category, 0)
        available = remaining.get(category, 0)
        if type(available) is not int or available < needed:
            return "insufficient_completion_capacity"
    if any(action.cost_envelope.values()):
        for key in ("wall_seconds_remaining", "core_seconds_remaining"):
            if key in state.budget and state.budget[key] <= 0:
                return "insufficient_completion_capacity"
    return None


def reduce_search_state(policy, decisions=(), evidence=(), native_budget=None, *, context=None, frontier=(), validation=None):
    """Reconstruct planner state; native_budget is supplied by its sole owner.

    decisions have ``action`` and ``fingerprint`` fields. Evidence records have
    the same fingerprint and an ``outcome``. The persistent controller provides
    these records; callers can also use this reducer directly in pure tests.
    """
    policy = policy if isinstance(policy, RouteSearchPolicy) else RouteSearchPolicy.from_dict(policy)
    outcomes = {item["decision_id"]: item for item in evidence}
    completed, pending = [], []
    for decision in decisions:
        action = SearchAction.from_dict(decision["action"])
        if decision["fingerprint"] != action.fingerprint(policy.validation_policy, context or {}):
            raise RouteSearchError("Decision fingerprint no longer matches scientific inputs")
        if decision["decision_id"] in outcomes:
            outcome = ActionOutcome.from_dict(outcomes[decision["decision_id"]]["outcome"])
            completed.append({**_copy(decision), "outcome": outcome.as_dict()})
        else:
            pending.append(_copy(decision))
    final = validation if isinstance(validation, RouteValidation) else RouteValidation.from_dict(validation or {})
    state = SearchSnapshot(policy, _copy(context or {}), tuple(completed), tuple(pending),
                           tuple(SearchAction.from_dict(item) for item in frontier), _copy(native_budget or {}), final)
    exclusions = tuple({"fingerprint": action.fingerprint(policy.validation_policy, state.context), "reason": reason}
                       for action in state.frontier if (reason := _eligibility(action, state)))
    status = final.status if final.status != "unresolved" else "searching"
    if status == "searching" and any(item["outcome"]["applicability"] == "unsupported" for item in completed):
        status = "unsupported_physical_model"
    if status == "searching":
        status = "ownership_unresolved" if pending else (_budget_status(state.budget) or "searching")
    if status == "searching" and not eligible_actions(state):
        status = "budget_exhausted" if any(item["reason"] == "insufficient_completion_capacity" for item in exclusions) else "search_space_exhausted_under_policy"
    return SearchSnapshot(policy, state.context, state.completed, state.pending, state.frontier, state.budget, final, status, exclusions)


def eligible_actions(state):
    """A finite distinct frontier, sorted without calibrated success guesses."""
    unique = {}
    for action in state.frontier:
        if _eligibility(action, state) is None:
            fingerprint = action.fingerprint(state.policy.validation_policy, state.context)
            unique.setdefault(fingerprint, action)
    order = {kind: index for index, kind in enumerate(state.policy.strategy_order)}
    return tuple(sorted(unique.values(), key=lambda item: (order[item.kind], item.rank,
                      item.fingerprint(state.policy.validation_policy, state.context))))


def choose_next_action(state):
    if state.status != "searching" or state.pending or _budget_status(state.budget):
        return None
    choices = eligible_actions(state)
    return choices[0] if choices else None


def _read(path):
    if path.is_symlink() or not path.is_file():
        raise RouteSearchError(f"Missing regular search evidence: {path}")
    data = path.read_bytes()
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise RouteSearchError("Duplicate search evidence key")
            result[key] = value
        return result
    try:
        value = json.loads(data, object_pairs_hook=pairs)
        if _bytes(value) != data:
            raise RouteSearchError("Noncanonical search evidence")
        return value
    except (ValueError, UnicodeError) as error:
        raise RouteSearchError("Invalid search evidence") from error


def _atomic(path, value, replace=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".search-tmp-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(_bytes(value))
            stream.flush()
            os.fsync(stream.fileno())
        if replace:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)  # atomic publication, never overwrite evidence
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class RouteSearchController:
    """Append decisions around injected existing stage adapters, never dispatch.

    ``recover_pending`` must reopen authoritative owned evidence without
    launching a second calculation. Return None for live/ambiguous work. An
    explicit ``not_started`` outcome may permit the original action again only
    when no reservation was ever made. Execution-only retries after reservation
    require a different explicit recovery action bound to its stopped receipt.
    """

    def __init__(self, root, policy, context=None):
        self.root = Path(root).absolute()
        self.policy = policy if isinstance(policy, RouteSearchPolicy) else RouteSearchPolicy.from_dict(policy)
        self.context = _copy(context or {})
        self.contract = {"schema_version": 1, "kind": "storca_adaptive_route_search",
                         "policy": self.policy.as_dict(), "context": self.context}

    @contextmanager
    def _ownership(self):
        if self.contract != {"schema_version": 1, "kind": "storca_adaptive_route_search",
                             "policy": self.policy.as_dict(), "context": _copy(self.context)}:
            raise RouteSearchError("Route policy/context was mutated after controller creation")
        if self.root.is_symlink():
            raise RouteSearchError("Search root cannot be a symlink")
        self.root.mkdir(parents=True, exist_ok=True)
        lock_path = self.root / ".search.lock"
        if lock_path.is_symlink():
            raise RouteSearchError("Search ownership file cannot be a symlink")
        with lock_path.open("a+b") as stream:
            try:
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RouteSearchError("Another route-search planner owns this workspace") from error
            try:
                contract_path = self.root / "contract.json"
                if contract_path.exists():
                    if _read(contract_path) != self.contract:
                        raise RouteSearchError("Immutable route-search contract changed")
                else:
                    if (self.root / "events").exists() and any((self.root / "events").iterdir()):
                        raise RouteSearchError("Missing search contract cannot reset decision history")
                    _atomic(contract_path, self.contract)
                events = self.root / "events"
                if events.is_symlink():
                    raise RouteSearchError("Search events cannot be a symlink")
                events.mkdir(exist_ok=True)
                yield
            finally:
                fcntl.flock(stream, fcntl.LOCK_UN)

    def _events(self):
        events, previous = [], _hash(self.contract)
        decisions, completed = {}, set()
        for path in sorted((self.root / "events").iterdir()):
            if path.name.startswith(".search-tmp-"):
                continue
            value = _read(path)
            expected_sequence = len(events) + 1
            if path.name != f"{expected_sequence:06d}-{_hash(value)}.json" or value.get("previous_sha256") != previous:
                raise RouteSearchError("Search decision journal changed or has a gap")
            if value.get("event") == "decision":
                identifier = value.get("decision_id")
                if identifier != expected_sequence or identifier in decisions:
                    raise RouteSearchError("Invalid search decision identity")
                decisions[identifier] = value
            elif value.get("event") == "outcome":
                identifier = value.get("decision_id")
                if identifier not in decisions or identifier in completed:
                    raise RouteSearchError("Orphan or repeated action outcome")
                ActionOutcome.from_dict(value["outcome"])
                completed.add(identifier)
            else:
                raise RouteSearchError("Unknown search journal event")
            events.append(value)
            previous = _hash(value)
        return events

    def _append(self, fields):
        events = self._events()
        value = {**fields, "previous_sha256": _hash(events[-1]) if events else _hash(self.contract)}
        _atomic(self.root / "events" / f"{len(events)+1:06d}-{_hash(value)}.json", value)
        return value

    def _snapshot(self, budget, frontier=(), validation=None):
        events = self._events()
        return reduce_search_state(self.policy, [item for item in events if item["event"] == "decision"],
            [item for item in events if item["event"] == "outcome"], budget, context=self.context,
            frontier=frontier, validation=validation)

    def _checkpoint(self, state):
        result = state.as_dict()
        _atomic(self.root / "checkpoint.json", result, replace=True)
        return result

    def run(self, actions, execute, budget, validate_route, recover_pending=None, cancelled=None):
        """Run one evidence-driven loop; return its reconstructible checkpoint.

        Unexpected callback exceptions propagate with the immutable decision
        left pending. They are not silently converted into physical failures.
        A canonical validator is called after each outcome, including the last
        available job; exhausting a budget does not hide already earned success.
        """
        with self._ownership():
            while True:
                state = self._snapshot(budget())
                if state.pending:
                    for decision in state.pending:
                        outcome = recover_pending(_copy(decision)) if recover_pending else None
                        if outcome is None:
                            return self._checkpoint(state)
                        outcome = ActionOutcome.from_dict(outcome)
                        self._append({"event": "outcome", "decision_id": decision["decision_id"],
                                      "outcome": outcome.as_dict(), "recovered": True})
                    continue
                final = RouteValidation.from_dict(validate_route(state))
                frontier = tuple(SearchAction.from_dict(item) for item in actions(state))
                state = self._snapshot(budget(), frontier, final)
                if cancelled and cancelled():
                    state = SearchSnapshot(state.policy, state.context, state.completed, state.pending,
                                           state.frontier, state.budget, state.validation, "cancelled", state.exclusions)
                if state.status != "searching":
                    return self._checkpoint(state)
                action = choose_next_action(state)
                identifier = len(self._events()) + 1
                decision = self._append({"event": "decision", "decision_id": identifier,
                    "fingerprint": action.fingerprint(self.policy.validation_policy, self.context),
                    "action": action.as_dict(), "budget_snapshot": _copy(state.budget)})
                outcome = ActionOutcome.from_dict(execute(action))
                self._append({"event": "outcome", "decision_id": decision["decision_id"], "outcome": outcome.as_dict(), "recovered": False})
                if outcome.execution in {"blocked_resources", "ownership_unresolved", "cancelled"} or outcome.applicability == "unsupported":
                    state = self._snapshot(budget(), frontier)
                    status = "unsupported_physical_model" if outcome.applicability == "unsupported" else outcome.execution
                    state = SearchSnapshot(state.policy, state.context, state.completed, state.pending,
                                           state.frontier, state.budget, state.validation, status, state.exclusions)
                    return self._checkpoint(state)
