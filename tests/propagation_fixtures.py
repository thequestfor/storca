"""Explicit synthetic reactor evidence for non-native orchestration tests."""


def isothermal_model(temperature=298.15, pressure=1.0):
    return {
        "policy_id": "isothermal-ideal-gas-constant-pressure-v1",
        "temperature_K": temperature,
        "pressure_bar": pressure,
        "energy_enabled": False,
        "condition_validation": "passed",
        "observed_temperature_range_K": [temperature, temperature],
        "observed_pressure_range_bar": [pressure, pressure],
    }
