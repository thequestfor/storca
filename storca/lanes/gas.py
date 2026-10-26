"""Lane declarations adapting the existing stability gas machinery."""

from __future__ import annotations

from typing import Any, Mapping

from .base import EnvironmentSpec, LaneDefinition


THERMAL_GAS_CAPABILITIES = (
    "rmg_bounded_gas_discovery",
    "orca_ground_state_path_verification",
    "arkane_thermal_rate",
    "verified_rmg_library_reinjection",
    "cantera_gas_network_propagation",
    "observable_convergence_gate",
)
PHOTOCHEMICAL_GAS_CAPABILITIES = (
    *THERMAL_GAS_CAPABILITIES,
    "condition_bound_photolysis_evidence",
)
IMPLEMENTED_GAS_CAPABILITIES = frozenset(PHOTOCHEMICAL_GAS_CAPABILITIES)
GAS_LANE_SCENARIOS = {
    "gas-low-pressure-thermal-v1": "low-pressure-intrinsic-gas-screen",
    "gas-dry-inert-v1": "dry-inert-gas-screen",
    "gas-dry-air-v1": "ambient-air-gas-screen",
    "gas-humid-air-v1": "humid-air-gas-screen",
    "gas-sunlight-dry-air-v1": "ambient-air-gas-screen",
    "gas-sunlight-humid-air-v1": "humid-air-gas-screen",
}


def _definition(
    lane_id: str,
    legacy_stage_id: str,
    *,
    illuminated: bool = False,
    applicability: str,
) -> LaneDefinition:
    return LaneDefinition(
        lane_id=lane_id,
        environment_kind="illuminated_homogeneous_gas" if illuminated else "homogeneous_gas",
        support_level="observable_supported",
        stage_kind="photolysis_rmg" if illuminated else "thermal_rmg",
        applicability=applicability,
        discovery_backends=("rmg", "graph_xtb_optional"),
        path_verifier="orca_endpoint_neb_ts_frequency_irc",
        rate_model="photolysis_plus_arkane_thermal" if illuminated else "arkane_thermal",
        propagator="repaired_rmg_model_plus_independent_cantera",
        convergence_policy="verified_rate_reinjection_flux_and_t95_robustness",
        supported_observables=("t95", "retention_at_duration"),
        required_capabilities=(
            PHOTOCHEMICAL_GAS_CAPABILITIES if illuminated else THERMAL_GAS_CAPABILITIES
        ),
        optional_capabilities=("graph_xtb_target_only_discovery",),
        benchmark_profile="current_storca_gas_contract_regression_v1",
        metadata={
            "legacy_stage_id": legacy_stage_id,
            "scientific_support_scope": (
                "Observable support is conditional on all existing route, rate, condition, "
                "propagation, and convergence evidence gates."
            ),
        },
    )


def gas_lane_definitions() -> tuple[LaneDefinition, ...]:
    return (
        _definition(
            "gas-low-pressure-thermal-v1", "dark-low-pressure",
            applicability="Dilute target in nitrogen at low pressure; homogeneous gas surrogate, not literal vacuum.",
        ),
        _definition(
            "gas-dry-inert-v1", "dark-dry-inert",
            applicability="Dilute target in dry nitrogen under the retained homogeneous gas contract.",
        ),
        _definition(
            "gas-dry-air-v1", "dark-dry-air",
            applicability="Dilute target in dry oxygen/nitrogen under the retained homogeneous gas contract.",
        ),
        _definition(
            "gas-humid-air-v1", "dark-humid-air",
            applicability="Dilute target in humid gas; water vapour only, never a liquid or droplet model.",
        ),
        _definition(
            "gas-sunlight-dry-air-v1", "sunlight-dry-air", illuminated=True,
            applicability="Dry homogeneous air with condition-bound photolysis evidence and ground-state thermal chemistry.",
        ),
        _definition(
            "gas-sunlight-humid-air-v1", "sunlight-humid-air", illuminated=True,
            applicability="Humid homogeneous gas with condition-bound photolysis evidence; no liquid-water chemistry.",
        ),
    )


def gas_lane_by_stage_id(stage_id: str) -> LaneDefinition:
    for definition in gas_lane_definitions():
        if definition.metadata.get("legacy_stage_id") == stage_id:
            return definition
    raise KeyError(f"No gas stability lane for stage: {stage_id}")


def environment_for_stage(stage: Mapping[str, Any]) -> EnvironmentSpec:
    definition = gas_lane_by_stage_id(str(stage.get("id") or ""))
    parameters = {
        key: stage.get(key)
        for key in (
            "scenario", "temperature_K", "pressure_bar", "relative_humidity",
            "retention_fraction", "maximum_duration_seconds", "light_model",
        )
        if stage.get(key) is not None
    }
    scenario = stage.get("scenario_config")
    if isinstance(scenario, Mapping):
        parameters["scenario"] = scenario.get("name")
        parameters["phase_model"] = scenario.get("phase")
    return EnvironmentSpec(
        environment_kind=definition.environment_kind,
        model_name="storca_condition_spec_gas_v1",
        parameters=parameters,
    )


__all__ = [
    "GAS_LANE_SCENARIOS", "IMPLEMENTED_GAS_CAPABILITIES",
    "PHOTOCHEMICAL_GAS_CAPABILITIES",
    "THERMAL_GAS_CAPABILITIES", "environment_for_stage", "gas_lane_by_stage_id",
    "gas_lane_definitions",
]
