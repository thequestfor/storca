"""Bounded, versioned computational light-screen assumptions.

This model deliberately provides a low/central/high *projection*, not a
measured quantum yield.  It is designed to make sunlight stages comparable in
cost and interpretation across molecules.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math


@dataclass(frozen=True)
class ComputationalLightModel:
    name: str = "computational-light-screen-v2"
    maximum_routes: int = 3
    energy_softening_eV: float = 0.35
    minimum_oscillator_strength: float = 0.01
    minimum_route_accessibility: float = 1e-6
    # The actual admitted window is intersected with the immutable source
    # spectrum.  Starting at 180 nm permits UV photolysis screens (for example
    # a declared 254-nm lamp) without pretending that terrestrial AM1.5G has
    # photons where its supplied spectrum has none.
    solar_window_nm: tuple[float, float] = (180.0, 800.0)
    low_reactive_fraction: float = 0.001
    central_reactive_fraction: float = 0.01
    high_reactive_fraction: float = 0.10

    def __post_init__(self) -> None:
        if (self.maximum_routes < 1 or not math.isfinite(self.energy_softening_eV)
                or self.energy_softening_eV <= 0):
            raise ValueError("Light-model route count and energy softening must be positive")
        if (not math.isfinite(self.minimum_oscillator_strength) or self.minimum_oscillator_strength < 0
                or not math.isfinite(self.minimum_route_accessibility)
                or not 0 <= self.minimum_route_accessibility < 1):
            raise ValueError("Minimum route accessibility must be from zero up to one")
        if (len(self.solar_window_nm) != 2 or not all(math.isfinite(value) for value in self.solar_window_nm)
                or not 0 < self.solar_window_nm[0] < self.solar_window_nm[1]):
            raise ValueError("Solar window must be an increasing pair of positive finite wavelengths")
        if not 0 < self.low_reactive_fraction <= self.central_reactive_fraction <= self.high_reactive_fraction <= 1:
            raise ValueError("Reactive fractions must be ordered probabilities from zero through one")

    def as_dict(self) -> dict:
        return asdict(self)

    def profiles(self) -> dict[str, float]:
        return {"low": self.low_reactive_fraction, "central": self.central_reactive_fraction,
                "high": self.high_reactive_fraction}


def energetic_accessibility(excited_state_eV: float, route_energy_eV: float, *, softening_eV: float = 0.35) -> float:
    """Return a smooth 0–1 availability score rather than a discontinuous cutoff."""
    if (not all(math.isfinite(value) for value in (excited_state_eV, route_energy_eV, softening_eV))
            or softening_eV <= 0):
        raise ValueError("softening_eV must be positive")
    exponent = max(-700.0, min(700.0, (excited_state_eV - route_energy_eV) / softening_eV))
    return 1.0 / (1.0 + math.exp(-exponent))


def distribute_reactive_branch(routes: list[dict], reactive_fraction: float) -> list[dict]:
    """Allocate one profile's reactive photon budget across accessible routes."""
    if not 0 <= reactive_fraction <= 1:
        raise ValueError("reactive_fraction must be from zero through one")
    scores = [max(0.0, float(route.get("accessibility", 0.0))) for route in routes]
    if not all(math.isfinite(score) for score in scores):
        raise ValueError("Route accessibility values must be finite")
    total = sum(scores)
    if total == 0:
        return [{**route, "quantum_yield": 0.0} for route in routes]
    # Accessibility must not be normalized away: if all routes are only 1%
    # accessible, at most 1%-of-a-route worth of the profile prior is admitted.
    # Once the summed score exceeds one, cap the total branch at the declared
    # reactive fraction and distribute it proportionally.
    denominator = max(1.0, total)
    return [{**route, "quantum_yield": reactive_fraction * score / denominator}
            for route, score in zip(routes, scores)]
