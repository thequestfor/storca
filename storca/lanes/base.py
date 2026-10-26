"""Immutable scientific contracts for stability-campaign lanes.

A lane is an auditable composition of existing discovery, verification, rate,
and propagation capabilities.  Registration alone never authorizes an
observable claim; :func:`admit_lane` performs a separate runtime applicability
decision and :class:`LaneResult` prevents an unadmitted lane from retaining a
verified lifetime.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
import re
from typing import Any, Iterable, Mapping


LANE_CONTRACT_SCHEMA_VERSION = 1
CAMPAIGN_CONTRACT_SCHEMA_VERSION = 1
LANE_SUPPORT_LEVELS = {
    "declared_only", "screening_supported", "observable_supported",
}
LANE_ADMISSION_STATUSES = {
    "supported", "screening_only", "inputs_missing", "capability_unavailable",
    "outside_method_applicability", "unresolved",
}


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
        _jsonable(payload), sort_keys=True, separators=(",", ":"),
        allow_nan=False,
    )
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class EnvironmentSpec:
    """General environment identity without weakening gas-only ConditionSpec.

    Lane-specific physical contracts remain responsible for validating their
    own parameters.  For the current gas lanes, ``condition_contract`` points
    to the existing strict :class:`storca.conditions.ConditionSpec` artifact.
    """

    environment_kind: str
    model_name: str
    parameters: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.environment_kind.strip() or not self.model_name.strip():
            raise ValueError("EnvironmentSpec needs environment_kind and model_name")
        object.__setattr__(self, "parameters", _jsonable(dict(self.parameters)))

    def as_dict(self) -> dict[str, Any]:
        payload = {
            "schema_version": LANE_CONTRACT_SCHEMA_VERSION,
            "kind": "storca_environment_spec",
            "environment_kind": self.environment_kind,
            "model_name": self.model_name,
            "parameters": _jsonable(self.parameters),
        }
        return {**payload, "contract_id": _contract_id(payload)}


@dataclass(frozen=True)
class LaneDefinition:
    """Static scientific capability declaration for one environment lane."""

    lane_id: str
    environment_kind: str
    support_level: str
    stage_kind: str
    applicability: str
    discovery_backends: tuple[str, ...]
    path_verifier: str
    rate_model: str
    propagator: str
    convergence_policy: str
    supported_observables: tuple[str, ...]
    required_capabilities: tuple[str, ...] = ()
    required_inputs: tuple[str, ...] = ()
    benchmark_profile: str | None = None
    optional_capabilities: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not re.fullmatch(r"[a-z0-9][a-z0-9_.-]*", self.lane_id):
            raise ValueError("lane_id must be a stable lowercase identifier")
        if self.support_level not in LANE_SUPPORT_LEVELS:
            raise ValueError(f"Unsupported lane support level: {self.support_level}")
        for name in (
            "environment_kind", "stage_kind", "applicability", "path_verifier",
            "rate_model", "propagator", "convergence_policy",
        ):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"LaneDefinition needs {name}")
        for name in (
            "discovery_backends", "supported_observables", "required_capabilities",
            "required_inputs", "optional_capabilities",
        ):
            object.__setattr__(self, name, tuple(str(item) for item in getattr(self, name)))
        object.__setattr__(self, "metadata", _jsonable(dict(self.metadata)))

    def as_dict(self) -> dict[str, Any]:
        payload = {
            "schema_version": LANE_CONTRACT_SCHEMA_VERSION,
            "kind": "storca_stability_lane_definition",
            "lane_id": self.lane_id,
            "environment_kind": self.environment_kind,
            "support_level": self.support_level,
            "stage_kind": self.stage_kind,
            "applicability": self.applicability,
            "discovery_backends": list(self.discovery_backends),
            "path_verifier": self.path_verifier,
            "rate_model": self.rate_model,
            "propagator": self.propagator,
            "convergence_policy": self.convergence_policy,
            "supported_observables": list(self.supported_observables),
            "required_capabilities": list(self.required_capabilities),
            "required_inputs": list(self.required_inputs),
            "optional_capabilities": list(self.optional_capabilities),
            "benchmark_profile": self.benchmark_profile,
            "metadata": _jsonable(self.metadata),
        }
        return {**payload, "contract_id": _contract_id(payload)}


@dataclass(frozen=True)
class LaneAdmission:
    """Runtime decision about whether a lane may attempt an observable."""

    lane_id: str
    status: str
    support_level: str
    requested_observable: str
    reason: str
    missing_capabilities: tuple[str, ...] = ()
    missing_inputs: tuple[str, ...] = ()
    environment_contract_id: str | None = None

    def __post_init__(self) -> None:
        if self.status not in LANE_ADMISSION_STATUSES:
            raise ValueError(f"Unsupported lane admission status: {self.status}")

    @property
    def permits_observable(self) -> bool:
        return self.status == "supported" and self.support_level == "observable_supported"

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": LANE_CONTRACT_SCHEMA_VERSION,
            "kind": "storca_lane_admission",
            "lane_id": self.lane_id,
            "status": self.status,
            "support_level": self.support_level,
            "requested_observable": self.requested_observable,
            "permits_observable": self.permits_observable,
            "reason": self.reason,
            "missing_capabilities": list(self.missing_capabilities),
            "missing_inputs": list(self.missing_inputs),
            "environment_contract_id": self.environment_contract_id,
        }


def admit_lane(
    definition: LaneDefinition,
    environment: EnvironmentSpec,
    *,
    requested_observable: str,
    available_capabilities: Iterable[str],
    available_inputs: Iterable[str] = (),
) -> LaneAdmission:
    """Fail closed unless the declared environment and capability set agree."""
    environment_contract = environment.as_dict()
    common = {
        "lane_id": definition.lane_id,
        "support_level": definition.support_level,
        "requested_observable": requested_observable,
        "environment_contract_id": environment_contract["contract_id"],
    }
    if environment.environment_kind != definition.environment_kind:
        return LaneAdmission(
            **common, status="outside_method_applicability",
            reason=(
                f"Lane expects environment kind '{definition.environment_kind}', not "
                f"'{environment.environment_kind}'."
            ),
        )
    if definition.support_level == "declared_only":
        return LaneAdmission(
            **common, status="outside_method_applicability",
            reason="The lane is registered for planning but has no supported scientific execution model.",
        )
    capabilities = {str(item) for item in available_capabilities}
    missing_capabilities = tuple(sorted(set(definition.required_capabilities) - capabilities))
    if missing_capabilities:
        return LaneAdmission(
            **common, status="capability_unavailable",
            reason="One or more required scientific capabilities are unavailable.",
            missing_capabilities=missing_capabilities,
        )
    inputs = {str(item) for item in available_inputs}
    missing_inputs = tuple(sorted(set(definition.required_inputs) - inputs))
    alternatives = tuple(
        tuple(str(item) for item in alternative)
        for alternative in definition.metadata.get("required_input_alternatives", ())
    )
    alternative_satisfied = any(set(alternative) <= inputs for alternative in alternatives)
    if missing_inputs and not alternative_satisfied:
        described_missing = set(missing_inputs)
        for alternative in alternatives:
            described_missing.update(set(alternative) - inputs)
        return LaneAdmission(
            **common, status="inputs_missing",
            reason="One or more required environment inputs are missing.",
            missing_inputs=tuple(sorted(described_missing)),
        )
    if requested_observable not in definition.supported_observables:
        return LaneAdmission(
            **common, status="outside_method_applicability",
            reason=f"Lane does not support requested observable '{requested_observable}'.",
        )
    if definition.support_level == "screening_supported":
        return LaneAdmission(
            **common, status="screening_only",
            reason="Lane can produce screening evidence but cannot issue this observable as a verified prediction.",
        )
    return LaneAdmission(
        **common, status="supported",
        reason="Environment, observable, and declared scientific capabilities satisfy the lane admission contract.",
    )


@dataclass(frozen=True)
class LaneResult:
    """Small standardized summary layered over an existing stage artifact."""

    lane_id: str
    lane_contract_id: str
    admission_status: str
    support_level: str
    execution_status: str
    observable_status: str
    requested_observable: str
    observable_value: float | None = None
    observable_units: str | None = None
    observable_data: Mapping[str, Any] = field(default_factory=dict)
    reason: str | None = None

    def __post_init__(self) -> None:
        data = _jsonable(dict(self.observable_data))
        if (
            self.observable_status == "verified"
            and self.observable_value is None and not data
        ):
            raise ValueError("A verified lane observable needs retained scalar or structured data")
        if self.observable_value is not None and self.observable_status != "verified":
            raise ValueError("A lane observable value is only valid with observable_status='verified'")
        if data and self.observable_status != "verified":
            raise ValueError("Structured lane observable data is only valid when verified")
        if self.observable_value is not None:
            value = float(self.observable_value)
            positive_only = self.requested_observable in {
                "t95", "half_life", "induction_time",
            }
            if not math.isfinite(value) or value < 0.0 or (positive_only and value <= 0.0):
                qualifier = "positive" if positive_only else "nonnegative"
                raise ValueError(f"Lane observable value must be {qualifier} and finite")
            if self.admission_status != "supported" or self.support_level != "observable_supported":
                raise ValueError("An unadmitted lane cannot retain a verified observable")
        if data and (
            self.admission_status != "supported" or self.support_level != "observable_supported"
        ):
            raise ValueError("An unadmitted lane cannot retain verified observable data")
        object.__setattr__(self, "observable_data", data)

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": LANE_CONTRACT_SCHEMA_VERSION,
            "kind": "storca_stability_lane_result",
            "lane_id": self.lane_id,
            "lane_contract_id": self.lane_contract_id,
            "admission_status": self.admission_status,
            "support_level": self.support_level,
            "execution_status": self.execution_status,
            "observable_status": self.observable_status,
            "requested_observable": self.requested_observable,
            "observable_value": self.observable_value,
            "observable_units": self.observable_units,
            "observable_data": self.observable_data,
            "reason": self.reason,
        }


@dataclass(frozen=True)
class CampaignSpec:
    """Immutable request shared by every lane in one stability campaign."""

    mode: str
    requested_observable: str
    lane_ids: tuple[str, ...]
    compute_budget: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.mode not in {"worst_case_search", "environment_map"}:
            raise ValueError("Unsupported campaign mode")
        if not self.requested_observable.strip() or not self.lane_ids:
            raise ValueError("CampaignSpec needs an observable and at least one lane")
        object.__setattr__(self, "lane_ids", tuple(str(item) for item in self.lane_ids))
        object.__setattr__(self, "compute_budget", _jsonable(dict(self.compute_budget)))

    def as_dict(self) -> dict[str, Any]:
        payload = {
            "schema_version": CAMPAIGN_CONTRACT_SCHEMA_VERSION,
            "kind": "storca_stability_campaign_spec",
            "mode": self.mode,
            "requested_observable": self.requested_observable,
            "lane_ids": list(self.lane_ids),
            "compute_budget": _jsonable(self.compute_budget),
        }
        return {**payload, "contract_id": _contract_id(payload)}


__all__ = [
    "CAMPAIGN_CONTRACT_SCHEMA_VERSION", "LANE_ADMISSION_STATUSES",
    "LANE_CONTRACT_SCHEMA_VERSION", "LANE_SUPPORT_LEVELS", "CampaignSpec",
    "EnvironmentSpec", "LaneAdmission", "LaneDefinition", "LaneResult",
    "admit_lane",
]
