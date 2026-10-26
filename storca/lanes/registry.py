"""Registry and default scientific declarations for stability lanes."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from .base import LaneDefinition


@dataclass
class LaneRegistry:
    """Deterministic lane registry; duplicate IDs are configuration errors."""

    _definitions: dict[str, LaneDefinition] = field(default_factory=dict)

    def register(self, definition: LaneDefinition) -> None:
        if definition.lane_id in self._definitions:
            raise ValueError(f"Duplicate stability lane: {definition.lane_id}")
        self._definitions[definition.lane_id] = definition

    def register_many(self, definitions: Iterable[LaneDefinition]) -> None:
        for definition in definitions:
            self.register(definition)

    def get(self, lane_id: str) -> LaneDefinition:
        try:
            return self._definitions[str(lane_id)]
        except KeyError as error:
            raise KeyError(f"Unknown stability lane: {lane_id}") from error

    def definitions(self) -> tuple[LaneDefinition, ...]:
        return tuple(self._definitions[key] for key in sorted(self._definitions))

    def as_dict(self) -> dict:
        definitions = self.definitions()
        return {
            "schema_version": 1,
            "kind": "storca_stability_lane_registry",
            "lane_count": len(definitions),
            "lanes": [definition.as_dict() for definition in definitions],
        }


def future_lane_declarations() -> tuple[LaneDefinition, ...]:
    """Recognized future environments that deliberately remain unsupported."""
    common = {
        "support_level": "declared_only",
        "stage_kind": "external_specialist",
        "discovery_backends": (),
        "path_verifier": "not_implemented",
        "rate_model": "not_implemented",
        "propagator": "not_implemented",
        "convergence_policy": "no_observable_claim",
        "supported_observables": (),
        "benchmark_profile": None,
    }
    return (
        LaneDefinition(
            lane_id="solution-assisted-microsolvation-v1",
            environment_kind="microsolvated_solution",
            applicability="Declared future explicit solvent-assisted path lane; no rate or propagation model is implemented.",
            **common,
        ),
        LaneDefinition(
            lane_id="surface-periodic-v1",
            environment_kind="periodic_surface",
            applicability="Declared future periodic surface lane; no adsorption, surface-rate, or microkinetic model is implemented.",
            **common,
        ),
        LaneDefinition(
            lane_id="electrochemical-interface-v1",
            environment_kind="electrochemical_interface",
            applicability="Declared future potential-aware interface lane; no constant-potential rate model is implemented.",
            **common,
        ),
    )


def build_stability_lane_registry() -> LaneRegistry:
    from .gas import gas_lane_definitions
    from .solution import implicit_solution_lane_definition
    from .solution_mass_action import implicit_solution_mass_action_lane_definition

    registry = LaneRegistry()
    registry.register_many(gas_lane_definitions())
    registry.register(implicit_solution_lane_definition())
    registry.register(implicit_solution_mass_action_lane_definition())
    registry.register_many(future_lane_declarations())
    return registry


__all__ = [
    "LaneRegistry", "build_stability_lane_registry", "future_lane_declarations",
]
