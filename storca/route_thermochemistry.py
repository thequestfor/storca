"""Validated endpoint thermochemistry and conservative no-saddle rate bounds."""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any

from .rmg_evidence import normalize_rmg_label
from .route_verify import RouteSpec


HARTREE_TO_J_MOL = 2625499.6394799
GAS_CONSTANT_J_MOL_K = 8.31446261815324
AVOGADRO = 6.02214076e23


def _last_float(text: str, pattern: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.IGNORECASE | re.MULTILINE)
    if not matches:
        return None
    try:
        value = float(matches[-1])
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def parse_orca_thermochemistry(output: Path) -> dict[str, Any]:
    """Parse the condition-bearing totals from a normally terminated ORCA frequency job."""
    output = Path(output)
    if not output.is_file():
        return {"status": "missing_output", "valid": False, "source_output": str(output)}
    text = output.read_text(errors="replace")
    values = {
        "electronic_energy_hartree": _last_float(
            text, r"^\s*FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)\s*$",
        ),
        "enthalpy_hartree": _last_float(text, r"^\s*Total Enthalpy\s+\.\.\.\s+([-+0-9.Ee]+)\s+Eh"),
        "gibbs_free_energy_hartree": _last_float(
            text, r"^\s*Final Gibbs free energy\s+\.\.\.\s+([-+0-9.Ee]+)\s+Eh",
        ),
        "zero_point_energy_hartree": _last_float(
            text, r"^\s*Zero point energy\s+\.\.\.\s+([-+0-9.Ee]+)\s+Eh",
        ),
        "temperature_K": _last_float(text, r"^\s*Temperature\s+\.\.\.\s+([-+0-9.Ee]+)\s+K"),
    }
    normal = "ORCA TERMINATED NORMALLY" in text
    required = ("electronic_energy_hartree", "enthalpy_hartree", "gibbs_free_energy_hartree", "temperature_K")
    valid = normal and all(values[key] is not None for key in required)
    return {
        "schema_version": 1,
        "kind": "orca_frequency_thermochemistry",
        "status": "validated" if valid else "incomplete",
        "valid": valid,
        "normal_termination": normal,
        **values,
        "source_output": str(output),
    }


def assemble_route_thermochemistry(
    route: RouteSpec,
    stationary_points: dict,
    condition_contract: dict | None,
) -> dict[str, Any]:
    """Sum validated species totals into one direction-specific reaction thermochemistry record."""
    if not condition_contract:
        return {"status": "condition_contract_unavailable", "valid": False}
    temperature = float(condition_contract["temperature_K"])
    sides: dict[str, list[dict]] = {}
    for side in ("reactants", "products"):
        records = []
        for endpoint in stationary_points.get(side) or []:
            frequency_output = endpoint.get("frequency_output")
            parsed = parse_orca_thermochemistry(Path(frequency_output)) if frequency_output else {
                "status": "frequency_output_unavailable", "valid": False,
            }
            records.append({
                "label": endpoint.get("label"),
                "smiles": endpoint.get("smiles"),
                "charge": endpoint.get("charge"),
                "multiplicity": endpoint.get("multiplicity"),
                "thermochemistry": parsed,
            })
        sides[side] = records
    all_records = [record for records in sides.values() for record in records]
    valid = bool(all_records) and all(
        record["thermochemistry"].get("valid")
        and math.isclose(
            float(record["thermochemistry"]["temperature_K"]), temperature,
            rel_tol=0.0, abs_tol=0.05,
        )
        for record in all_records
    )
    result: dict[str, Any] = {
        "schema_version": 1,
        "kind": "validated_route_endpoint_thermochemistry",
        "status": "validated" if valid else "incomplete_or_condition_mismatch",
        "valid": valid,
        "temperature_K": temperature,
        "reactant_molecularity": len(route.reactant_smiles),
        "product_molecularity": len(route.product_smiles),
        "reactants": sides["reactants"],
        "products": sides["products"],
    }
    if not valid:
        return result
    totals = {}
    for field in ("electronic_energy_hartree", "enthalpy_hartree", "gibbs_free_energy_hartree"):
        reactant = sum(float(record["thermochemistry"][field]) for record in sides["reactants"])
        product = sum(float(record["thermochemistry"][field]) for record in sides["products"])
        totals[field] = {"reactants": reactant, "products": product, "reaction": product - reactant}
    result["totals_hartree"] = totals
    result["reaction"] = {
        "delta_e_hartree": totals["electronic_energy_hartree"]["reaction"],
        "delta_h_hartree": totals["enthalpy_hartree"]["reaction"],
        "delta_g_hartree": totals["gibbs_free_energy_hartree"]["reaction"],
        "delta_g_kcal_mol": totals["gibbs_free_energy_hartree"]["reaction"] * 627.5094740631,
    }
    return result


def _collision_ceiling_si(smiles_pair: tuple[str, str], temperature_K: float,
                          collision_diameter_angstrom: float,
                          collision_safety_factor: float) -> dict[str, Any]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_pair]
    if any(molecule is None for molecule in molecules):
        raise ValueError("Collision ceiling requires two RDKit-readable species")
    masses_g_mol = [float(Descriptors.MolWt(molecule)) for molecule in molecules]
    masses_kg = [mass / 1000.0 / AVOGADRO for mass in masses_g_mol]
    reduced_mass = masses_kg[0] * masses_kg[1] / sum(masses_kg)
    diameter_m = collision_diameter_angstrom * 1.0e-10
    relative_speed = math.sqrt(
        8.0 * 1.380649e-23 * temperature_K / (math.pi * reduced_mass)
    )
    molecular_rate_m3_s = math.pi * diameter_m ** 2 * relative_speed
    return {
        "value_si_m3_mol_s": molecular_rate_m3_s * AVOGADRO * collision_safety_factor,
        "unscaled_hard_sphere_value_si_m3_mol_s": molecular_rate_m3_s * AVOGADRO,
        "molecular_weights_g_mol": masses_g_mol,
        "collision_diameter_angstrom": collision_diameter_angstrom,
        "collision_safety_factor": collision_safety_factor,
        "model": "safety_factored_hard_sphere_gas_collision_ceiling",
    }


def condition_specific_forward_loss_bound(
    route: RouteSpec,
    route_thermochemistry: dict,
    condition_contract: dict | None,
    *,
    collision_diameter_angstrom: float = 5.0,
    collision_safety_factor: float = 10.0,
    thermochemistry_uncertainty_kcal_mol: float = 10.0,
) -> dict[str, Any]:
    """Bound forward target loss for a validated gas-phase 2-to-2 route.

    For equal molecularity, the concentration-standard-state factor cancels.
    Detailed balance gives ``kf/kr = exp(-deltaG/RT)``; bounding both forward
    and reverse bimolecular rates by their hard-sphere collision ceilings then
    gives a conservative forward ceiling.  No other molecularity is inferred.
    """
    if not condition_contract or not route_thermochemistry.get("valid"):
        return {"status": "thermochemistry_or_condition_unavailable", "applicable": False}
    if len(route.reactant_smiles) != 2 or len(route.product_smiles) != 2:
        return {
            "status": "unsupported_molecularity",
            "applicable": False,
            "supported_molecularity": "gas_phase_2_to_2_only",
        }
    if str(condition_contract.get("phase_model")) != "homogeneous gas-phase surrogate":
        return {"status": "unsupported_phase", "applicable": False}
    if len(route.reactant_charges) != 2 or len(route.product_charges) != 2:
        return {
            "status": "fragment_charge_contract_unavailable", "applicable": False,
            "reason": "A neutral-capture bound requires explicit charge state for every fragment.",
        }
    charges = [*route.reactant_charges, *route.product_charges]
    if any(int(charge) != 0 for charge in charges):
        return {
            "status": "ionic_capture_bound_unsupported", "applicable": False,
            "reason": "Long-range ionic capture is not safely bounded by the neutral hard-sphere model.",
        }
    if collision_safety_factor < 1.0 or thermochemistry_uncertainty_kcal_mol < 0.0:
        raise ValueError("Rate-bound safety margins must be nonnegative and collision factor at least one")
    temperature = float(condition_contract["temperature_K"])
    pressure_bar = float(condition_contract["pressure_bar"])
    duration = float(condition_contract["target_duration_seconds"])
    retention = float(condition_contract["retention_fraction"])
    delta_g_hartree = float(route_thermochemistry["reaction"]["delta_g_hartree"])
    # Lower delta-G by an explicit uncertainty allowance before bounding the
    # forward rate.  This makes screening harder, preventing modest electronic-
    # structure errors from becoming false reassurance through exponentiation.
    effective_delta_g_j_mol = (
        delta_g_hartree * HARTREE_TO_J_MOL
        - thermochemistry_uncertainty_kcal_mol * 4184.0
    )
    exponent = -effective_delta_g_j_mol / (GAS_CONSTANT_J_MOL_K * temperature)
    equilibrium_ratio = math.exp(max(-745.0, min(709.0, exponent)))
    forward_collision = _collision_ceiling_si(
        tuple(route.reactant_smiles), temperature, collision_diameter_angstrom,
        collision_safety_factor,
    )
    reverse_collision = _collision_ceiling_si(
        tuple(route.product_smiles), temperature, collision_diameter_angstrom,
        collision_safety_factor,
    )
    forward_upper = min(
        float(forward_collision["value_si_m3_mol_s"]),
        float(reverse_collision["value_si_m3_mol_s"]) * equilibrium_ratio,
    )

    target = normalize_rmg_label(str(condition_contract.get("target_label") or "stability"))
    labels = list(route.reactant_labels) if len(route.reactant_labels) == 2 else list(route.reactant_smiles)
    target_indices = [index for index, label in enumerate(labels) if normalize_rmg_label(str(label)) == target]
    if len(target_indices) != 1:
        return {
            "status": "target_stoichiometry_unsupported", "applicable": False,
            "target_reactant_occurrences": len(target_indices),
        }
    co_index = 1 - target_indices[0]
    composition = {
        normalize_rmg_label(str(label)): float(value)
        for label, value in (condition_contract.get("composition") or {}).items()
    }
    co_label = normalize_rmg_label(str(labels[co_index]))
    mole_fraction = composition.get(co_label)
    concentration_source = "declared_condition_composition"
    if mole_fraction is None:
        mole_fraction = 1.0
        concentration_source = "conservative_total_gas_concentration_ceiling"
    total_concentration_mol_m3 = pressure_bar * 1.0e5 / (GAS_CONSTANT_J_MOL_K * temperature)
    co_concentration = max(0.0, min(1.0, mole_fraction)) * total_concentration_mol_m3
    pseudo_first_order_upper = forward_upper * co_concentration
    loss_upper = -math.expm1(-min(745.0, pseudo_first_order_upper * duration))
    t_retention_lower = (
        -math.log(retention) / pseudo_first_order_upper if pseudo_first_order_upper > 0.0 else math.inf
    )
    threshold = 1.0 - retention
    status = (
        "forward_loss_below_retention_threshold_upper_bound"
        if loss_upper < threshold else "forward_loss_upper_bound_crosses_retention_threshold"
    )
    return {
        "schema_version": 1,
        "kind": "detailed_balance_forward_loss_upper_bound",
        "status": status,
        "applicable": True,
        "temperature_K": temperature,
        "pressure_bar": pressure_bar,
        "target_duration_seconds": duration,
        "retention_fraction": retention,
        "delta_g_hartree": delta_g_hartree,
        "thermochemistry_uncertainty_kcal_mol": thermochemistry_uncertainty_kcal_mol,
        "effective_conservative_delta_g_kcal_mol": effective_delta_g_j_mol / 4184.0,
        "equilibrium_forward_reverse_rate_ratio": equilibrium_ratio,
        "forward_collision_ceiling": forward_collision,
        "reverse_collision_ceiling": reverse_collision,
        "forward_rate_upper_bound_si_m3_mol_s": forward_upper,
        "co_reactant_label": labels[co_index],
        "co_reactant_mole_fraction_upper_bound": mole_fraction,
        "co_reactant_concentration_mol_m3": co_concentration,
        "co_reactant_concentration_source": concentration_source,
        "pseudo_first_order_rate_upper_bound_s": pseudo_first_order_upper,
        "target_loss_fraction_upper_bound": loss_upper,
        "t_retention_seconds_lower_bound": t_retention_lower if math.isfinite(t_retention_lower) else None,
        "retention_threshold_crossed_by_upper_bound": loss_upper >= threshold,
        "assumptions": [
            "Validated ideal-gas ORCA endpoint Gibbs energies at the condition temperature.",
            "Equal 2-to-2 molecularity, so the concentration standard-state factor cancels.",
            "Reverse and forward elementary rates cannot exceed retained hard-sphere collision ceilings.",
            "Collision ceilings include an explicit multiplicative safety factor for neutral long-range attraction.",
            "Reaction Gibbs energy is lowered by the stated uncertainty allowance before exponentiation.",
            "The bound is an upper limit, not a fitted or verified reaction rate.",
        ],
    }


__all__ = [
    "assemble_route_thermochemistry",
    "condition_specific_forward_loss_bound",
    "parse_orca_thermochemistry",
]
