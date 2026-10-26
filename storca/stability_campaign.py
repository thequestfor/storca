"""Immutable, multi-environment request contracts for stability campaigns.

This module contains no calculation execution.  It normalizes the scientific
request that later planners and lane adapters consume while leaving the
existing gas ``ConditionSpec`` and solution environment contracts intact.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

from .chemical_system import ChemicalSystem
from .lanes.base import CampaignSpec, EnvironmentSpec


STABILITY_CAMPAIGN_SCHEMA_VERSION = 1
STABILITY_CAMPAIGN_KIND = "storca_stability_campaign_request"
CAMPAIGN_MODES = {"environment_map", "worst_case_search"}
OBSERVABLE_KINDS = {
    "t95", "half_life", "retention_at_duration", "concentration_at_time",
    "product_yield", "dominant_product", "selectivity", "branching_fraction",
    "induction_time", "pathway_flux",
}
_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
_SHA256 = re.compile(r"sha256:[0-9a-f]{64}")


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            str(key): _jsonable(item)
            for key, item in sorted(value.items(), key=lambda entry: str(entry[0]))
        }
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    return value


def _contract_id(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        _jsonable(payload), sort_keys=True, separators=(",", ":"), allow_nan=False,
    )
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _text(value: Any, *, name: str, required: bool = True) -> str | None:
    if value is None and not required:
        return None
    text = str(value or "").strip()
    if required and not text:
        raise ValueError(f"{name} is required")
    return text or None


def _identifier(value: Any, *, name: str) -> str:
    identifier = str(_text(value, name=name))
    if not _IDENTIFIER.fullmatch(identifier):
        raise ValueError(f"{name} must be a stable letters/numbers identifier")
    return identifier


def _optional_positive(value: Any, *, name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite positive number")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be a finite positive number")
    return number


def _validated_budget(value: Mapping[str, Any] | None, *, name: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")

    def validate(item: Any, path: str) -> Any:
        if isinstance(item, Mapping):
            return {
                str(key): validate(entry, f"{path}.{key}")
                for key, entry in sorted(item.items(), key=lambda pair: str(pair[0]))
            }
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise ValueError(f"{path} must be a nonnegative finite number")
        number = float(item)
        if not math.isfinite(number) or number < 0.0:
            raise ValueError(f"{path} must be a nonnegative finite number")
        return int(item) if isinstance(item, int) else number

    return validate(value, name)


@dataclass(frozen=True)
class ObservableRequest:
    observable_id: str
    kind: str
    target_component_id: str
    retention_fraction: float | None = None
    time_horizon_seconds: float | None = None
    parameters: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        observable_id = _identifier(self.observable_id, name="observable_id")
        kind = str(_text(self.kind, name=f"{observable_id}.kind")).lower()
        if kind not in OBSERVABLE_KINDS:
            raise ValueError(f"Unsupported observable kind: {kind}")
        target = _identifier(
            self.target_component_id, name=f"{observable_id}.target_component_id",
        )
        retention = self.retention_fraction
        if retention is not None:
            retention = float(retention)
            if not math.isfinite(retention) or not 0.0 < retention < 1.0:
                raise ValueError(f"{observable_id}.retention_fraction must be in (0, 1)")
        horizon = _optional_positive(
            self.time_horizon_seconds, name=f"{observable_id}.time_horizon_seconds",
        )
        if kind == "t95":
            if retention is None:
                retention = 0.95
            if not math.isclose(retention, 0.95, rel_tol=0.0, abs_tol=1.0e-12):
                raise ValueError("The t95 observable requires retention_fraction=0.95")
        if kind in {"retention_at_duration", "concentration_at_time", "product_yield"}:
            if horizon is None:
                raise ValueError(f"Observable {kind} requires time_horizon_seconds")
        object.__setattr__(self, "observable_id", observable_id)
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "target_component_id", target)
        object.__setattr__(self, "retention_fraction", retention)
        object.__setattr__(self, "time_horizon_seconds", horizon)
        object.__setattr__(self, "parameters", _jsonable(dict(self.parameters)))
        _contract_id(self._payload())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ObservableRequest":
        if not isinstance(value, Mapping):
            raise TypeError("observable request must be a mapping")
        observable = cls(
            observable_id=value.get("observable_id"),
            kind=value.get("observable_kind") or value.get("kind"),
            target_component_id=value.get("target_component_id"),
            retention_fraction=value.get("retention_fraction"),
            time_horizon_seconds=value.get("time_horizon_seconds"),
            parameters=value.get("parameters") or {},
        )
        retained = value.get("contract_id")
        if retained is not None and retained != observable.as_dict()["contract_id"]:
            raise ValueError(f"Observable {observable.observable_id} contract_id mismatch")
        return observable

    def _payload(self) -> dict[str, Any]:
        return {
            "schema_version": STABILITY_CAMPAIGN_SCHEMA_VERSION,
            "kind": "storca_stability_observable_request",
            "observable_id": self.observable_id,
            "observable_kind": self.kind,
            "target_component_id": self.target_component_id,
            "retention_fraction": self.retention_fraction,
            "time_horizon_seconds": self.time_horizon_seconds,
            "parameters": _jsonable(self.parameters),
        }

    def as_dict(self) -> dict[str, Any]:
        payload = self._payload()
        return {**payload, "contract_id": _contract_id(payload)}


@dataclass(frozen=True)
class EvidenceInputReference:
    input_kind: str
    content_sha256: str
    input_id: str | None = None
    path: str | None = None
    claim_id: str | None = None

    def __post_init__(self) -> None:
        input_kind = _identifier(self.input_kind, name="evidence input_kind")
        input_id = _identifier(
            self.input_id or input_kind, name=f"{input_kind}.input_id",
        )
        digest = str(_text(self.content_sha256, name=f"{input_kind}.content_sha256")).lower()
        if not _SHA256.fullmatch(digest):
            raise ValueError(
                f"{input_kind}.content_sha256 must be 'sha256:' plus 64 hexadecimal digits"
            )
        path = _text(self.path, name=f"{input_kind}.path", required=False)
        claim_id = _text(self.claim_id, name=f"{input_kind}.claim_id", required=False)
        if path is None and claim_id is None:
            raise ValueError(f"Evidence input {input_kind} needs a path or claim_id")
        object.__setattr__(self, "input_kind", input_kind)
        object.__setattr__(self, "input_id", input_id)
        object.__setattr__(self, "content_sha256", digest)
        object.__setattr__(self, "path", path)
        object.__setattr__(self, "claim_id", claim_id)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "EvidenceInputReference":
        if not isinstance(value, Mapping):
            raise TypeError("evidence input reference must be a mapping")
        return cls(
            input_kind=value.get("input_kind"),
            content_sha256=value.get("content_sha256"),
            input_id=value.get("input_id"),
            path=value.get("path"),
            claim_id=value.get("claim_id"),
        )

    def scientific_dict(self) -> dict[str, Any]:
        return {
            "input_id": self.input_id,
            "input_kind": self.input_kind,
            "content_sha256": self.content_sha256,
            "claim_id": self.claim_id,
        }

    def as_dict(self) -> dict[str, Any]:
        return {**self.scientific_dict(), "path": self.path}


@dataclass(frozen=True)
class CampaignLaneRequest:
    lane_request_id: str
    lane_id: str
    environment: EnvironmentSpec
    observable_ids: tuple[str, ...]
    execution_order: int | None = None
    method_profile: str | None = None
    compute_budget: Mapping[str, Any] = field(default_factory=dict)
    discovery_policy: Mapping[str, Any] = field(default_factory=dict)
    evidence_inputs: tuple[EvidenceInputReference, ...] = ()
    user_assumptions: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        request_id = _identifier(self.lane_request_id, name="lane_request_id")
        lane_id = _identifier(self.lane_id, name=f"{request_id}.lane_id")
        observables = tuple(str(item) for item in self.observable_ids)
        if not observables or any(not item.strip() for item in observables):
            raise ValueError(f"Lane request {request_id} needs at least one observable_id")
        if len(observables) != len(set(observables)):
            raise ValueError(f"Lane request {request_id} has duplicate observable_ids")
        evidence = tuple(sorted(self.evidence_inputs, key=lambda item: str(item.input_id)))
        input_ids = [str(item.input_id) for item in evidence]
        if len(input_ids) != len(set(input_ids)):
            raise ValueError(f"Lane request {request_id} has duplicate evidence input_id values")
        execution_order = self.execution_order
        if execution_order is not None:
            if isinstance(execution_order, bool) or int(execution_order) != execution_order:
                raise ValueError(f"Lane request {request_id} execution_order must be an integer")
            execution_order = int(execution_order)
            if execution_order < 0:
                raise ValueError(f"Lane request {request_id} execution_order must be nonnegative")
        discovery = dict(self.discovery_policy)
        if lane_id.startswith("gas-"):
            allowed_discovery = {
                "rmg", "graph_xtb", "graph_xtb_max_candidates",
                "graph_xtb_energy_window_kcal_mol",
            }
            unknown_discovery = sorted(set(discovery) - allowed_discovery)
            if unknown_discovery:
                raise ValueError(
                    f"Lane request {request_id} has unknown discovery policy fields: "
                    f"{unknown_discovery}"
                )
            rmg_policy = str(discovery.get("rmg", "enabled"))
            graph_policy = str(discovery.get("graph_xtb", "disabled"))
            if rmg_policy != "enabled":
                raise ValueError("Gas campaign schema v1 requires rmg='enabled'")
            if graph_policy not in {"disabled", "bounded"}:
                raise ValueError("graph_xtb discovery policy must be 'disabled' or 'bounded'")
            raw_max_graph = discovery.get("graph_xtb_max_candidates", 4)
            if isinstance(raw_max_graph, bool) or int(raw_max_graph) != raw_max_graph:
                raise ValueError("graph_xtb_max_candidates must be an integer")
            max_graph = int(raw_max_graph)
            graph_window = float(discovery.get("graph_xtb_energy_window_kcal_mol", 80.0))
            if max_graph < 0 or not math.isfinite(graph_window):
                raise ValueError("Graph/xTB discovery bounds must be finite and nonnegative")
            discovery = {
                "rmg": rmg_policy,
                "graph_xtb": graph_policy,
                "graph_xtb_max_candidates": max_graph,
                "graph_xtb_energy_window_kcal_mol": graph_window,
            }
        elif discovery:
            raise ValueError(
                f"Lane request {request_id} cannot use the gas discovery_policy contract"
            )
        object.__setattr__(self, "lane_request_id", request_id)
        object.__setattr__(self, "lane_id", lane_id)
        object.__setattr__(self, "observable_ids", observables)
        object.__setattr__(self, "execution_order", execution_order)
        object.__setattr__(
            self, "method_profile", _text(
                self.method_profile, name=f"{request_id}.method_profile", required=False,
            ),
        )
        object.__setattr__(
            self, "compute_budget", _validated_budget(
                self.compute_budget, name=f"{request_id}.compute_budget",
            ),
        )
        object.__setattr__(self, "discovery_policy", discovery)
        object.__setattr__(self, "evidence_inputs", evidence)
        object.__setattr__(
            self, "user_assumptions", _jsonable(dict(self.user_assumptions)),
        )
        _contract_id(self.scientific_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "CampaignLaneRequest":
        if not isinstance(value, Mapping):
            raise TypeError("campaign lane request must be a mapping")
        environment = value.get("environment")
        if not isinstance(environment, Mapping):
            raise TypeError("lane environment must be a mapping")
        return cls(
            lane_request_id=value.get("lane_request_id"),
            lane_id=value.get("lane_id"),
            environment=EnvironmentSpec(
                environment_kind=environment.get("environment_kind"),
                model_name=environment.get("model_name"),
                parameters=environment.get("parameters") or {},
            ),
            observable_ids=tuple(value.get("observable_ids") or ()),
            execution_order=value.get("execution_order"),
            method_profile=value.get("method_profile"),
            compute_budget=value.get("compute_budget") or {},
            discovery_policy=value.get("discovery_policy") or {},
            evidence_inputs=tuple(
                EvidenceInputReference.from_dict(item)
                for item in value.get("evidence_inputs") or []
            ),
            user_assumptions=value.get("user_assumptions") or {},
        )

    def scientific_dict(self) -> dict[str, Any]:
        return {
            "lane_request_id": self.lane_request_id,
            "lane_id": self.lane_id,
            "environment": self.environment.as_dict(),
            "observable_ids": list(self.observable_ids),
            "execution_order": self.execution_order,
            "method_profile": self.method_profile,
            "compute_budget": _jsonable(self.compute_budget),
            "discovery_policy": _jsonable(self.discovery_policy),
            "evidence_inputs": [item.scientific_dict() for item in self.evidence_inputs],
            "user_assumptions": _jsonable(self.user_assumptions),
        }

    def as_dict(self) -> dict[str, Any]:
        value = self.scientific_dict()
        value["evidence_inputs"] = [item.as_dict() for item in self.evidence_inputs]
        return value


@dataclass(frozen=True)
class StabilityCampaignRequest:
    chemical_system: ChemicalSystem
    observables: tuple[ObservableRequest, ...]
    lane_requests: tuple[CampaignLaneRequest, ...]
    global_compute_budget: Mapping[str, Any] = field(default_factory=dict)
    campaign_policy: Mapping[str, Any] = field(default_factory=lambda: {
        "mode": "environment_map", "execution_policy": "prepare_only",
    })
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        observables = tuple(sorted(self.observables, key=lambda item: item.observable_id))
        lanes = tuple(sorted(
            self.lane_requests,
            key=lambda item: (
                item.execution_order is None,
                item.execution_order if item.execution_order is not None else 0,
                item.lane_request_id,
            ),
        ))
        if not observables or not lanes:
            raise ValueError("A stability campaign needs observables and lane requests")
        observable_ids = [item.observable_id for item in observables]
        lane_request_ids = [item.lane_request_id for item in lanes]
        if len(observable_ids) != len(set(observable_ids)):
            raise ValueError("Campaign observable_id values must be unique")
        if len(lane_request_ids) != len(set(lane_request_ids)):
            raise ValueError("Campaign lane_request_id values must be unique")
        component_ids = {item.component_id for item in self.chemical_system.components}
        if any(item.target_component_id not in component_ids for item in observables):
            raise ValueError("An observable target_component_id is absent from chemical_system")
        known_observables = set(observable_ids)
        if any(set(item.observable_ids) - known_observables for item in lanes):
            raise ValueError("A lane request references an unknown observable_id")

        from .lanes import build_stability_lane_registry

        registered = {item.lane_id for item in build_stability_lane_registry().definitions()}
        unknown_lanes = sorted({item.lane_id for item in lanes} - registered)
        if unknown_lanes:
            raise ValueError(f"Unknown stability lane IDs: {unknown_lanes}")

        policy = dict(self.campaign_policy)
        allowed_policy = {"mode", "execution_policy", "continue_after_blocked_lane"}
        unknown_policy = sorted(set(policy) - allowed_policy)
        if unknown_policy:
            raise ValueError(f"Unknown campaign policy fields: {unknown_policy}")
        mode = str(policy.get("mode", "environment_map"))
        if mode not in CAMPAIGN_MODES:
            raise ValueError(f"Unsupported campaign mode: {mode}")
        execution_policy = str(policy.get("execution_policy", "prepare_only"))
        if execution_policy not in {"prepare_only", "explicit_execute"}:
            raise ValueError(
                "Campaign execution_policy must be 'prepare_only' or 'explicit_execute'"
            )
        continue_after = policy.get("continue_after_blocked_lane", True)
        if not isinstance(continue_after, bool):
            raise ValueError("continue_after_blocked_lane must be boolean")
        policy = {
            "mode": mode,
            "execution_policy": execution_policy,
            "continue_after_blocked_lane": continue_after,
        }
        if mode == "worst_case_search" and len(lanes) > 1:
            if any(not item.lane_id.startswith("gas-") for item in lanes):
                raise ValueError(
                    "worst_case_search across multiple lane requests is currently restricted "
                    "to the comparable homogeneous-gas ladder; use environment_map for "
                    "cross-environment campaigns"
                )
            orders = [item.execution_order for item in lanes]
            if any(item is None for item in orders) or len(orders) != len(set(orders)):
                raise ValueError(
                    "worst_case_search with multiple gas lanes requires unique explicit "
                    "execution_order values"
                )
        object.__setattr__(self, "observables", observables)
        object.__setattr__(self, "lane_requests", lanes)
        object.__setattr__(
            self, "global_compute_budget", _validated_budget(
                self.global_compute_budget, name="global_compute_budget",
            ),
        )
        object.__setattr__(self, "campaign_policy", policy)
        object.__setattr__(self, "metadata", _jsonable(dict(self.metadata)))
        json.dumps(self.metadata, sort_keys=True, allow_nan=False)
        _contract_id(self.scientific_payload())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "StabilityCampaignRequest":
        if not isinstance(value, Mapping):
            raise TypeError("stability campaign request must be a mapping")
        version = int(value.get("schema_version", STABILITY_CAMPAIGN_SCHEMA_VERSION))
        if version != STABILITY_CAMPAIGN_SCHEMA_VERSION:
            raise ValueError(f"Unsupported stability campaign schema_version: {version}")
        kind = value.get("kind", STABILITY_CAMPAIGN_KIND)
        if kind != STABILITY_CAMPAIGN_KIND:
            raise ValueError(f"Unsupported stability campaign kind: {kind}")
        request = cls(
            chemical_system=ChemicalSystem.from_dict(value.get("chemical_system") or {}),
            observables=tuple(
                ObservableRequest.from_dict(item) for item in value.get("observables") or []
            ),
            lane_requests=tuple(
                CampaignLaneRequest.from_dict(item) for item in value.get("lane_requests") or []
            ),
            global_compute_budget=value.get("global_compute_budget") or {},
            campaign_policy=value.get("campaign_policy") or {},
            metadata=value.get("metadata") or {},
        )
        retained = value.get("contract_id")
        if retained is not None and retained != request.contract_id:
            raise ValueError("Stability campaign contract_id does not match its content")
        return request

    @classmethod
    def from_json(cls, path: Path) -> "StabilityCampaignRequest":
        path = Path(path)
        return cls.from_dict(json.loads(path.read_text()))

    def scientific_payload(self) -> dict[str, Any]:
        return {
            "schema_version": STABILITY_CAMPAIGN_SCHEMA_VERSION,
            "kind": STABILITY_CAMPAIGN_KIND,
            "chemical_system": self.chemical_system.as_dict(),
            "observables": [item.as_dict() for item in self.observables],
            "lane_requests": [item.scientific_dict() for item in self.lane_requests],
            "global_compute_budget": _jsonable(self.global_compute_budget),
            "campaign_policy": _jsonable(self.campaign_policy),
        }

    @property
    def contract_id(self) -> str:
        return _contract_id(self.scientific_payload())

    def as_dict(self) -> dict[str, Any]:
        payload = self.scientific_payload()
        payload["lane_requests"] = [item.as_dict() for item in self.lane_requests]
        return {**payload, "metadata": _jsonable(self.metadata), "contract_id": self.contract_id}

    def to_legacy_campaign_spec(self) -> CampaignSpec:
        """Return the old one-observable summary without changing its schema."""
        kinds = {
            observable.kind
            for observable in self.observables
            if any(observable.observable_id in lane.observable_ids for lane in self.lane_requests)
        }
        if len(kinds) != 1:
            raise ValueError("Legacy CampaignSpec can represent exactly one observable kind")
        return CampaignSpec(
            mode=str(self.campaign_policy["mode"]),
            requested_observable=next(iter(kinds)),
            lane_ids=tuple(dict.fromkeys(item.lane_id for item in self.lane_requests)),
            compute_budget=self.global_compute_budget,
        )


__all__ = [
    "CAMPAIGN_MODES", "OBSERVABLE_KINDS", "STABILITY_CAMPAIGN_KIND",
    "STABILITY_CAMPAIGN_SCHEMA_VERSION", "CampaignLaneRequest",
    "EvidenceInputReference", "ObservableRequest", "StabilityCampaignRequest",
]
