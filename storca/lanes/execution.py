"""Executable stability-lane boundary.

Lane registration describes scientific intent.  This module separately binds
an admitted definition to code and validates the standardized result before an
observable can leave the lane.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Protocol, runtime_checkable

from .base import EnvironmentSpec, LaneDefinition, LaneResult, admit_lane


@dataclass(frozen=True)
class LaneExecutionRequest:
    environment: EnvironmentSpec
    requested_observable: str
    inputs: Mapping[str, Any]
    compute_budget: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class LaneExecutionOutcome:
    lane_result: LaneResult
    evidence: Mapping[str, Any]
    artifacts: Mapping[str, Any] = field(default_factory=dict)


@runtime_checkable
class LaneExecutor(Protocol):
    """Minimal implementation interface for one scientific lane."""

    @property
    def definition(self) -> LaneDefinition: ...

    @property
    def capabilities(self) -> frozenset[str]: ...

    def execute(self, request: LaneExecutionRequest, workdir: Path) -> LaneExecutionOutcome: ...


def execute_lane(
    executor: LaneExecutor,
    request: LaneExecutionRequest,
    workdir: Path,
) -> dict[str, Any]:
    """Admit, execute, and validate a lane without fabricating fallback data."""
    definition = executor.definition
    available_inputs = tuple(str(key) for key, value in request.inputs.items() if value is not None)
    admission = admit_lane(
        definition,
        request.environment,
        requested_observable=request.requested_observable,
        available_capabilities=executor.capabilities,
        available_inputs=available_inputs,
    )
    contract = definition.as_dict()
    common = {
        "schema_version": 1,
        "kind": "storca_executable_lane_run",
        "lane_definition": contract,
        "environment_contract": request.environment.as_dict(),
        "admission": admission.as_dict(),
        "requested_observable": request.requested_observable,
    }
    if not admission.permits_observable:
        result = LaneResult(
            lane_id=definition.lane_id,
            lane_contract_id=contract["contract_id"],
            admission_status=admission.status,
            support_level=definition.support_level,
            execution_status="lane_not_admitted",
            observable_status="unresolved",
            requested_observable=request.requested_observable,
            reason=admission.reason,
        )
        return {**common, "status": "lane_not_admitted", "lane_result": result.as_dict(),
                "evidence": {}, "artifacts": {}}

    outcome = executor.execute(request, Path(workdir))
    result = outcome.lane_result
    if result.lane_id != definition.lane_id:
        raise ValueError("Lane executor returned a result for a different lane")
    if result.lane_contract_id != contract["contract_id"]:
        raise ValueError("Lane executor returned a result for a different lane contract")
    if result.admission_status != admission.status:
        raise ValueError("Lane executor result does not retain the runtime admission status")
    if result.requested_observable != request.requested_observable:
        raise ValueError("Lane executor changed the requested observable")
    return {
        **common,
        "status": result.execution_status,
        "lane_result": result.as_dict(),
        "evidence": dict(outcome.evidence),
        "artifacts": dict(outcome.artifacts),
    }


__all__ = [
    "LaneExecutionOutcome", "LaneExecutionRequest", "LaneExecutor", "execute_lane",
]
