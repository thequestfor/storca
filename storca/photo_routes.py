"""Generic, bounded photo-route discovery and ranking."""

from __future__ import annotations

from .intrinsic_initiation import enumerate_homolytic_cleavages
from .light_model import ComputationalLightModel, energetic_accessibility


def generic_photo_candidates(smiles: str, *, maximum_routes: int = 3) -> list[dict]:
    """Use the uniform homolysis enumerator; never use functional-group alerts."""
    candidates = enumerate_homolytic_cleavages(smiles)["candidates"]
    # Product energetics, not a hand-written alert list, decides final rank.
    # Canonical atom indices give deterministic bounded work before energies are available.
    ordered = sorted(candidates, key=lambda item: tuple(item["bond_atom_indices"]))[:maximum_routes]
    return [{"route_id": f"homolysis-{index + 1}", "kind": "homolysis", **candidate}
            for index, candidate in enumerate(ordered)]


def rank_photo_routes(candidates: list[dict], transitions: list[dict], route_energies_eV: dict[str, float],
                      model: ComputationalLightModel = ComputationalLightModel(),
                      source_window_nm: tuple[float, float] | None = None) -> list[dict]:
    """Rank computed routes by bright-state energetic accessibility."""
    window = source_window_nm or model.solar_window_nm
    illuminated_states = [state for state in transitions
                          if window[0] <= state["wavelength_nm"] <= window[1]
                          and state["oscillator_strength"] > 0]
    ranked = []
    for candidate in candidates:
        energy = route_energies_eV.get(candidate["route_id"])
        if energy is None:
            ranked.append({**candidate, "status": "energetics_incomplete", "route_energy_eV": None,
                           "accessibility": 0.0})
            continue
        availability = [
            (state, energetic_accessibility(state["energy_eV"], energy,
                                             softening_eV=model.energy_softening_eV))
            for state in illuminated_states
        ]
        oscillator_sum = sum(state["oscillator_strength"] for state, _ in availability)
        accessibility = (
            sum(state["oscillator_strength"] * value for state, value in availability) / oscillator_sum
            if oscillator_sum > 0 else 0.0
        )
        # Keep a dominant state for explanation, but use the smooth weighted
        # score above for kinetics so small oscillator-strength changes cannot
        # switch a route discontinuously between two raw accessibility values.
        dominant_state, _ = max(availability, key=lambda item: item[0]["oscillator_strength"] * item[1],
                                default=(None, 0.0))
        ranked.append({**candidate, "status": "computed", "route_energy_eV": energy,
                       "accessibility": accessibility, "dominant_excited_state": dominant_state,
                       "accessibility_model": "oscillator_strength_weighted_state_average",
                       "state_accessibility": [
                           {"state": state.get("state"), "wavelength_nm": state["wavelength_nm"],
                            "oscillator_strength": state["oscillator_strength"], "accessibility": value}
                           for state, value in availability
                       ]})
    return sorted(ranked, key=lambda route: route["accessibility"], reverse=True)
