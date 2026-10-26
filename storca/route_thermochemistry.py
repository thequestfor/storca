"""Validated endpoint thermochemistry and conservative no-saddle rate bounds."""

from __future__ import annotations

import math
import hashlib
from pathlib import Path
import re
from typing import Any

from .rmg_evidence import normalize_rmg_label
from .route_verify import RouteSpec


HARTREE_TO_J_MOL = 2625499.6394799
GAS_CONSTANT_J_MOL_K = 8.31446261815324
AVOGADRO = 6.02214076e23


def atomic_hydrogen_identity(smiles: str, charge: int, multiplicity: int) -> dict | None:
    """Identify the supported atom without pretending it has a molecular Hessian."""
    from rdkit import Chem
    parameters = Chem.SmilesParserParams()
    parameters.removeHs = False
    molecule = Chem.MolFromSmiles(smiles, parameters)
    if molecule is None:
        raise ValueError("Invalid atomic endpoint SMILES")
    molecule = Chem.AddHs(molecule)
    if molecule.GetNumAtoms() != 1:
        return None
    atom = molecule.GetAtomWithIdx(0)
    if (atom.GetAtomicNum() != 1 or charge != 0 or multiplicity != 2
            or atom.GetFormalCharge() != 0 or atom.GetIsotope() not in (0, 1, 2, 3)):
        raise ValueError("Monatomic endpoint modeling currently supports neutral ground-state H doublet only")
    isotope = int(atom.GetIsotope())
    table = Chem.GetPeriodicTable()
    mass = (table.GetMassForIsotope(1, isotope) if isotope
            else table.GetMostCommonIsotopeMass(1))
    return {"element": "H", "isotope": isotope, "mass_amu": float(mass), "atomic_mass_amu": float(mass),
            "mass_policy": "declared_isotope" if isotope else "unspecified_H_assumed_protium",
            "electronic_state_model": "neutral_H_ground_doublet", "electronic_degeneracy": 2}


def read_atomic_hydrogen_energy(input_xyz: Path, energy_input: Path, energy_output: Path,
                               *, smiles: str, charge: int, multiplicity: int) -> dict:
    """Reopen native atom SP evidence; no synthetic Opt/Freq/mode evidence."""
    import numpy as np
    from .route_geometry import read_xyz

    identity = atomic_hydrogen_identity(smiles, charge, multiplicity)
    if identity is None:
        raise ValueError("Atomic evidence requires exactly one hydrogen nucleus")
    input_xyz, energy_input, energy_output = map(Path, (input_xyz, energy_input, energy_output))
    elements, coordinates, _ = read_xyz(input_xyz)
    if elements != ["H"] or not np.isfinite(coordinates).all():
        raise ValueError("Atomic input must contain exactly one finite H nucleus")
    native_input = energy_input.read_text()
    keywords = " ".join(line[1:] for line in native_input.splitlines() if line.lstrip().startswith("!"))
    if re.search(r"\b(?:Opt\w*|Freq|NumFreq|NEB\S*|IRC|GOAT)\b", keywords, re.I):
        raise ValueError("Atomic evidence requires an energy-only job, not a geometry/frequency job")
    states = re.findall(r"(?mi)^\s*\*\s*xyzfile\s+([-+]?\d+)\s+(\d+)\s+(\S+)\s*$", native_input)
    inline = re.findall(r"(?mis)^\s*\*\s*xyz\s+([-+]?\d+)\s+(\d+)\s*\n(.*?)^\s*\*\s*$", native_input)
    if len(states) + len(inline) != 1:
        raise ValueError("Atomic input state contract is missing or inconsistent")
    selected = (states or inline)[0]
    if tuple(map(int, selected[:2])) != (charge, multiplicity):
        raise ValueError("Atomic input state contract is missing or inconsistent")
    if states:
        referenced = (energy_input.parent / selected[2].strip('"')).resolve(strict=True)
        if referenced != input_xyz.resolve(strict=True):
            raise ValueError("Atomic energy input does not reference the bound nucleus geometry")
    else:
        rows = [line.split() for line in selected[2].splitlines() if line.strip()]
        if (len(rows) != 1 or len(rows[0]) != 4 or rows[0][0] != "H"
                or not np.allclose([[float(value) for value in rows[0][1:]]], coordinates,
                                   atol=1e-9, rtol=0)):
            raise ValueError("Atomic native input coordinates do not match the bound nucleus geometry")
    text = energy_output.read_text(errors="replace")
    if (text.count("ORCA TERMINATED NORMALLY") != 1 or "SCF CONVERGED AFTER" not in text
            or re.search(r"ORCA.*(?:error termination|TERMINATED ABNORMALLY)", text, re.I)):
        raise ValueError("Atomic energy requires unique normal completion and converged SCF")
    charges = re.findall(r"(?m)^\s*Total Charge\s+Charge\s+\.{2,}\s+([-+]?\d+)\s*$", text)
    spins = re.findall(r"(?m)^\s*Multiplicity\s+Mult\s+\.{2,}\s+(\d+)\s*$", text)
    if not charges or not spins or set(map(int, charges)) != {charge} or set(map(int, spins)) != {multiplicity}:
        raise ValueError("Atomic native output does not establish the declared electronic state")
    geometry_count = 0
    for section in text.split("CARTESIAN COORDINATES (ANGSTROEM)")[1:]:
        rows = []
        for line in section.splitlines():
            fields = line.split()
            if len(fields) == 4 and re.fullmatch(r"[A-Z][a-z]?", fields[0]):
                rows.append((fields[0], [float(value) for value in fields[1:]]))
            elif rows:
                break
        if rows:
            if [row[0] for row in rows] != elements or not np.allclose(
                    [row[1] for row in rows], coordinates, atol=1e-5, rtol=0):
                raise ValueError("Atomic native output nucleus geometry differs from its input")
            geometry_count += 1
    if not geometry_count:
        raise ValueError("Atomic native output nucleus geometry is absent")
    energies = re.findall(r"(?m)^\s*FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)\s*$", text)
    if len(energies) != 1 or not math.isfinite(float(energies[0])):
        raise ValueError("A unique finite native atomic electronic energy is required")
    return {"kind": "monatomic_stationary_evidence", "status": "validated_atom", **identity,
            "smiles": smiles, "charge": charge, "multiplicity": multiplicity,
            "input_xyz": str(input_xyz), "optimized_xyz": str(input_xyz),
            "energy_input": str(energy_input), "energy_output": str(energy_output),
            "electronic_energy_hartree": float(energies[0]),
            "source_evidence_hashes": {name: hashlib.sha256(path.read_bytes()).hexdigest()
                for name, path in (("input_xyz", input_xyz), ("energy_input", energy_input),
                                   ("energy_output", energy_output))},
            "optimization_required": False, "vibrational_hessian_required": False,
            "scope": "isolated_neutral_ground_state_atom_not_a_molecular_minimum"}


def atomic_hydrogen_thermochemistry(endpoint: dict, temperature_K: float,
                                   pressure_pa: float = 101325.0) -> dict:
    """Source-bound ideal-gas translation plus declared electronic degeneracy."""
    reopened = read_atomic_hydrogen_energy(
        Path(endpoint["input_xyz"]), Path(endpoint["energy_input"]), Path(endpoint["energy_output"]),
        smiles=endpoint["smiles"], charge=endpoint["charge"], multiplicity=endpoint["multiplicity"],
    )
    for key in ("source_evidence_hashes", "electronic_energy_hartree", "mass_amu", "atomic_mass_amu", "isotope"):
        if endpoint.get(key) != reopened[key]:
            raise ValueError("Atomic stationary evidence changed or declared energy/mass was altered")
    temperature, pressure = float(temperature_K), float(pressure_pa)
    if not math.isfinite(temperature) or temperature <= 0 or not math.isfinite(pressure) or pressure <= 0:
        raise ValueError("Atomic ideal-gas thermochemistry requires finite positive T and standard pressure")
    mass_kg = reopened["mass_amu"] * 1.66053906660e-27
    q_translation = (2.0 * math.pi * mass_kg * 1.380649e-23 * temperature / 6.62607015e-34 ** 2) ** 1.5
    q_translation *= 1.380649e-23 * temperature / pressure
    entropy = GAS_CONSTANT_J_MOL_K * (math.log(q_translation) + 2.5 + math.log(2.0))
    thermal_enthalpy = 2.5 * GAS_CONSTANT_J_MOL_K * temperature
    electronic = reopened["electronic_energy_hartree"]
    return {"kind": "ideal_gas_atomic_thermochemistry", "status": "validated", "valid": True,
            "temperature_K": temperature, "standard_pressure_pa": pressure,
            "electronic_energy_hartree": electronic, "zero_point_energy_hartree": 0.0,
            "enthalpy_hartree": electronic + thermal_enthalpy / HARTREE_TO_J_MOL,
            "gibbs_free_energy_hartree": electronic + (thermal_enthalpy - temperature * entropy) / HARTREE_TO_J_MOL,
            "entropy_J_mol_K": entropy, "mass_amu": reopened["mass_amu"],
            "translational_degrees": 3, "rotational_degrees": 0, "vibrational_degrees": 0,
            "electronic_degeneracy": 2, "source_evidence_hashes": reopened["source_evidence_hashes"],
            "assumptions": ["isolated_ideal_gas_atom", "neutral_H_ground_doublet", reopened["mass_policy"]]}


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
            if endpoint.get("kind") == "monatomic_stationary_evidence":
                try:
                    parsed = atomic_hydrogen_thermochemistry(endpoint, temperature)
                except (OSError, ValueError, KeyError, TypeError) as error:
                    parsed = {"status": "atomic_evidence_invalid", "valid": False, "reason": str(error)}
                records.append({"label": endpoint.get("label"), "smiles": endpoint.get("smiles"),
                                "charge": endpoint.get("charge"), "multiplicity": endpoint.get("multiplicity"),
                                "thermochemistry": parsed})
                continue
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
    thermochemistry_standard_pressure_pa: float = 101325.0,
) -> dict[str, Any]:
    """Bound forward target loss for selected neutral gas-phase routes.

    For 2-to-2 exchange, the concentration-standard-state factor cancels.
    For 1-to-2 dissociation, detailed balance gives
    ``kdiss/kassoc = exp(-deltaG/RT) * c_standard``.  Bounding the reverse
    bimolecular association by a deliberately inflated hard-sphere collision
    ceiling therefore bounds the forward unimolecular loss rate without
    pretending to provide VTST or RRKM kinetics.  Other molecularities are not
    inferred.
    """
    if not condition_contract or not route_thermochemistry.get("valid"):
        return {"status": "thermochemistry_or_condition_unavailable", "applicable": False}
    molecularity = (len(route.reactant_smiles), len(route.product_smiles))
    if molecularity not in {(2, 2), (1, 2)}:
        return {
            "status": "unsupported_molecularity",
            "applicable": False,
            "supported_molecularities": ["gas_phase_2_to_2", "gas_phase_1_to_2_dissociation"],
        }
    if condition_contract.get("phase_model") not in (
        "homogeneous gas-phase surrogate",
        "homogeneous gas-phase low-pressure surrogate",
    ):
        return {"status": "unsupported_phase", "applicable": False}
    if (len(route.reactant_charges), len(route.product_charges)) != molecularity:
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
    if (collision_safety_factor < 1.0 or thermochemistry_uncertainty_kcal_mol < 0.0
            or thermochemistry_standard_pressure_pa <= 0.0):
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
    reverse_collision = _collision_ceiling_si(
        tuple(route.product_smiles), temperature, collision_diameter_angstrom,
        collision_safety_factor,
    )

    target = normalize_rmg_label(str(condition_contract.get("target_label") or "stability"))
    labels = (
        list(route.reactant_labels)
        if len(route.reactant_labels) == len(route.reactant_smiles)
        else list(route.reactant_smiles)
    )
    target_indices = [index for index, label in enumerate(labels) if normalize_rmg_label(str(label)) == target]
    if len(target_indices) != 1:
        return {
            "status": "target_stoichiometry_unsupported", "applicable": False,
            "target_reactant_occurrences": len(target_indices),
        }
    if molecularity == (2, 2):
        forward_collision = _collision_ceiling_si(
            tuple(route.reactant_smiles), temperature, collision_diameter_angstrom,
            collision_safety_factor,
        )
        forward_upper = min(
            float(forward_collision["value_si_m3_mol_s"]),
            float(reverse_collision["value_si_m3_mol_s"]) * equilibrium_ratio,
        )
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
        bound_basis = "equal_molecularity_detailed_balance_and_bidirectional_collision_ceilings"
        rate_units = "m^3/(mol*s)"
        standard_concentration = None
    else:
        if target_indices != [0]:
            return {
                "status": "target_stoichiometry_unsupported", "applicable": False,
                "reason": "The supported 1-to-2 bound requires the target to be the sole reactant.",
            }
        # ORCA ideal-gas thermochemistry uses a pressure standard state.  Its
        # value is explicit here because it does not cancel for 1 -> 2.
        standard_concentration = (
            thermochemistry_standard_pressure_pa / (GAS_CONSTANT_J_MOL_K * temperature)
        )
        forward_upper = (
            float(reverse_collision["value_si_m3_mol_s"])
            * equilibrium_ratio * standard_concentration
        )
        forward_collision = None
        mole_fraction = None
        concentration_source = None
        co_concentration = None
        pseudo_first_order_upper = forward_upper
        bound_basis = "unimolecular_dissociation_detailed_balance_from_reverse_collision_ceiling"
        rate_units = "s^-1"
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
        "phase_model": condition_contract["phase_model"],
        "temperature_K": temperature,
        "pressure_bar": pressure_bar,
        "target_duration_seconds": duration,
        "retention_fraction": retention,
        "delta_g_hartree": delta_g_hartree,
        "thermochemistry_uncertainty_kcal_mol": thermochemistry_uncertainty_kcal_mol,
        "effective_conservative_delta_g_kcal_mol": effective_delta_g_j_mol / 4184.0,
        "equilibrium_forward_reverse_rate_ratio": equilibrium_ratio,
        "molecularity": {"reactants": molecularity[0], "products": molecularity[1]},
        "bound_basis": bound_basis,
        "forward_collision_ceiling": forward_collision,
        "reverse_collision_ceiling": reverse_collision,
        "forward_rate_upper_bound_si": forward_upper,
        "forward_rate_upper_bound_units": rate_units,
        # Retain the historical explicit field for 2 -> 2 consumers.
        "forward_rate_upper_bound_si_m3_mol_s": (
            forward_upper if molecularity == (2, 2) else None
        ),
        "thermochemistry_standard_pressure_pa": (
            thermochemistry_standard_pressure_pa if molecularity == (1, 2) else None
        ),
        "standard_concentration_mol_m3": standard_concentration,
        "co_reactant_label": labels[co_index] if molecularity == (2, 2) else None,
        "co_reactant_mole_fraction_upper_bound": mole_fraction,
        "co_reactant_concentration_mol_m3": co_concentration,
        "co_reactant_concentration_source": concentration_source,
        "pseudo_first_order_rate_upper_bound_s": pseudo_first_order_upper,
        "target_loss_fraction_upper_bound": loss_upper,
        "t_retention_seconds_lower_bound": t_retention_lower if math.isfinite(t_retention_lower) else None,
        "retention_threshold_crossed_by_upper_bound": loss_upper >= threshold,
        "assumptions": [
            "Validated ideal-gas ORCA endpoint Gibbs energies at the condition temperature.",
            (
                "Equal 2-to-2 molecularity, so the concentration standard-state factor cancels."
                if molecularity == (2, 2) else
                "The explicit ideal-gas thermochemistry pressure standard state is converted to concentration for 1-to-2 detailed balance."
            ),
            (
                "Reverse and forward elementary rates cannot exceed retained hard-sphere collision ceilings."
                if molecularity == (2, 2) else
                "The reverse bimolecular association rate cannot exceed the retained hard-sphere collision ceiling."
            ),
            "Collision ceilings include an explicit multiplicative safety factor for neutral long-range attraction.",
            "Reaction Gibbs energy is lowered by the stated uncertainty allowance before exponentiation.",
            "The bound is an upper limit, not a fitted or verified reaction rate.",
            "The bound does not replace VTST, RRKM/master-equation falloff, capture probability, or recrossing calculations when it crosses the requested loss threshold.",
        ],
    }


__all__ = [
    "assemble_route_thermochemistry",
    "condition_specific_forward_loss_bound",
    "parse_orca_thermochemistry",
]
