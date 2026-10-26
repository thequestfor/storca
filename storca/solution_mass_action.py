"""Verified-input nonlinear mass-action kinetics for homogeneous solution.

This module is deliberately separate from :mod:`storca.solution_kinetics`.
The original lane is a linear, monotonic, first-order model; accepting a
second reactant there would silently change both its mathematics and its
stored-evidence contract.  This sibling model admits only explicitly verified
first- or second-order solution rates and never converts gas-phase kinetics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import itertools
import math
from typing import Any, Mapping

import numpy as np

from .solution_kinetics import (
    ImplicitSolutionEnvironment,
    PROHIBITED_MECHANISM_FLAGS,
    SolutionSpecies,
    SUPPORTED_RATE_STATUSES,
    SUPPORTED_RATE_THEORIES,
    _finite,
    _text,
)


MASS_ACTION_NETWORK_SCHEMA_VERSION = 1
MASS_ACTION_COVERAGE_STATUS = "verified_complete_within_declared_mass_action_scope"
FIRST_ORDER_UNITS = {"s^-1", "s-1", "1/s"}
SECOND_ORDER_UNITS = {
    "l*mol^-1*s^-1", "lmol^-1s^-1", "lmol-1s-1", "l/(mol*s)",
    "m^-1*s^-1", "m^-1s^-1", "m-1s-1",
}
MASS_ACTION_RATE_THEORIES = SUPPORTED_RATE_THEORIES | {"solution_diffusion_corrected_rate"}
MASS_ACTION_RATE_STATUSES = SUPPORTED_RATE_STATUSES | {"verified_solution_diffusion_corrected_rate"}


def _normalized_units(value: Any) -> str:
    return str(value or "").strip().casefold().replace(" ", "")


@dataclass(frozen=True)
class MassActionSolutionEnvironment:
    """Dark, neutral, ideal-activity, constant-volume implicit solution."""

    solvent: str
    solvent_model: str
    temperature_K: float
    target_duration_seconds: float
    retention_fraction: float = 0.95
    standard_state_mol_L: float = 1.0
    volume_L: float = 1.0
    light_condition: str = "dark"
    electronic_structure_profile: str = "user_verified_solution_profile"
    activity_model: str = "ideal_molar_concentration"

    def __post_init__(self) -> None:
        # Reuse the original lane's strict solvent and scalar validation.
        validated = ImplicitSolutionEnvironment(
            solvent=self.solvent,
            solvent_model=self.solvent_model,
            temperature_K=self.temperature_K,
            target_duration_seconds=self.target_duration_seconds,
            retention_fraction=self.retention_fraction,
            standard_state_mol_L=self.standard_state_mol_L,
            volume_L=self.volume_L,
            light_condition=self.light_condition,
            electronic_structure_profile=self.electronic_structure_profile,
        )
        if str(self.activity_model).strip() != "ideal_molar_concentration":
            raise ValueError("The mass-action solution lane supports ideal molar activities only")
        for name in (
            "solvent", "solvent_model", "temperature_K", "target_duration_seconds",
            "retention_fraction", "standard_state_mol_L", "volume_L", "light_condition",
            "electronic_structure_profile",
        ):
            object.__setattr__(self, name, getattr(validated, name))
        object.__setattr__(self, "activity_model", "ideal_molar_concentration")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "MassActionSolutionEnvironment":
        if not isinstance(value, Mapping):
            raise TypeError("mass-action solution environment must be a mapping")
        return cls(
            solvent=value.get("solvent"),
            solvent_model=value.get("solvent_model"),
            temperature_K=value.get("temperature_K"),
            target_duration_seconds=value.get("target_duration_seconds"),
            retention_fraction=value.get("retention_fraction", 0.95),
            standard_state_mol_L=value.get("standard_state_mol_L", 1.0),
            volume_L=value.get("volume_L", 1.0),
            light_condition=value.get("light_condition", "dark"),
            electronic_structure_profile=value.get(
                "electronic_structure_profile", "user_verified_solution_profile",
            ),
            activity_model=value.get("activity_model", "ideal_molar_concentration"),
        )

    def as_environment_spec(self):
        from .lanes.base import EnvironmentSpec

        return EnvironmentSpec(
            environment_kind="homogeneous_solution",
            model_name="implicit_solvent_constant_volume_mass_action_v1",
            parameters={
                "solvent": self.solvent,
                "solvent_model": self.solvent_model,
                "temperature_K": self.temperature_K,
                "target_duration_seconds": self.target_duration_seconds,
                "retention_fraction": self.retention_fraction,
                "standard_state_mol_L": self.standard_state_mol_L,
                "volume_L": self.volume_L,
                "light_condition": self.light_condition,
                "electronic_structure_profile": self.electronic_structure_profile,
                "activity_model": self.activity_model,
                "phase_model": "homogeneous implicit solution",
            },
        )


@dataclass(frozen=True)
class VerifiedMassActionReaction:
    reaction_id: str
    reactants: Mapping[str, int]
    products: Mapping[str, int]
    rate_constant: float
    rate_units: str
    temperature_K: float
    solvent: str
    solvent_model: str
    standard_state_mol_L: float
    verification_status: str
    rate_theory: str
    source: str
    rate_claim_id: str
    path_claim_id: str | None = None
    standard_state_correction: str | None = None
    calculation_hashes: tuple[str, ...] = ()
    evidence_hashes: tuple[str, ...] = ()
    lower_rate_multiplier: float | None = None
    upper_rate_multiplier: float | None = None
    uncertainty_source: str | None = None
    uncertainty_components: tuple[str, ...] = ()
    diffusion_limit_L_mol_s: float | None = None
    diffusion_limit_source: str | None = None
    reversible_partner_id: str | None = None
    applicability_flags: Mapping[str, bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        reaction_id = _text(self.reaction_id, name="reaction_id")

        def validated_stoichiometry(value: Mapping[str, int], *, side: str) -> dict[str, int]:
            result: dict[str, int] = {}
            for key, coefficient in dict(value).items():
                species_id = str(key).strip()
                if (
                    not species_id or isinstance(coefficient, bool)
                    or int(coefficient) != coefficient or int(coefficient) <= 0
                ):
                    raise ValueError(
                        f"Reaction {reaction_id} {side} stoichiometry must use positive integers"
                    )
                result[species_id] = int(coefficient)
            return result

        reactants = validated_stoichiometry(self.reactants, side="reactant")
        products = validated_stoichiometry(self.products, side="product")
        if not reactants or not products:
            raise ValueError(f"Reaction {reaction_id} needs reactants and products")
        order = sum(reactants.values())
        if order not in {1, 2}:
            raise ValueError("Mass-action solution reactions must be first or second order")
        rate = _finite(self.rate_constant, name=f"{reaction_id}.rate_constant")
        if rate <= 0.0:
            raise ValueError("Verified solution rates must be positive")
        units = _normalized_units(self.rate_units)
        allowed_units = FIRST_ORDER_UNITS if order == 1 else SECOND_ORDER_UNITS
        if units not in allowed_units:
            expected = "s^-1" if order == 1 else "L mol^-1 s^-1"
            raise ValueError(f"Reaction {reaction_id} order {order} requires rate units {expected}")
        status = _text(self.verification_status, name="verification_status")
        theory = _text(self.rate_theory, name="rate_theory")
        if status not in MASS_ACTION_RATE_STATUSES or theory not in MASS_ACTION_RATE_THEORIES:
            raise ValueError(f"Unsupported verified solution-rate contract for {reaction_id}")
        path_claim = str(self.path_claim_id).strip() if self.path_claim_id is not None else None
        correction = (
            str(self.standard_state_correction).strip()
            if self.standard_state_correction is not None else None
        )
        hashes = tuple(_text(item, name="calculation_hash") for item in self.calculation_hashes)
        evidence_hashes = tuple(_text(item, name="evidence_hash") for item in self.evidence_hashes)
        if not evidence_hashes or any(
            len(item) != 71 or not item.startswith("sha256:")
            or any(character not in "0123456789abcdef" for character in item[7:].lower())
            for item in evidence_hashes
        ):
            raise ValueError("Every verified mass-action rate requires sha256 evidence hashes")
        if status == "verified_solution_tst":
            if not path_claim or correction != "one_molar_solution" or not hashes:
                raise ValueError(
                    "Verified solution-TST rates require path evidence, one-molar correction, and hashes"
                )
        lower = self.lower_rate_multiplier
        upper = self.upper_rate_multiplier
        uncertainty_source = (
            str(self.uncertainty_source).strip() if self.uncertainty_source is not None else None
        )
        if lower is None or upper is None or not uncertainty_source:
            raise ValueError("Every mass-action rate requires explicit lower/upper uncertainty evidence")
        lower = _finite(lower, name="lower_rate_multiplier")
        upper = _finite(upper, name="upper_rate_multiplier")
        if lower < 0.0 or lower > 1.0 or upper < 1.0:
            raise ValueError("Rate uncertainty must satisfy 0 <= lower <= 1 <= upper")
        uncertainty_components = tuple(
            _text(item, name="uncertainty_component") for item in self.uncertainty_components
        )
        if not uncertainty_components:
            raise ValueError("Every mass-action rate requires named uncertainty components")
        if order == 2 and not any(
            "diffusion" in item.casefold() or "encounter" in item.casefold()
            for item in uncertainty_components
        ):
            raise ValueError(
                "Second-order rate uncertainty must include diffusion or encounter effects"
            )
        diffusion_limit = self.diffusion_limit_L_mol_s
        diffusion_source = (
            str(self.diffusion_limit_source).strip()
            if self.diffusion_limit_source is not None else None
        )
        if order == 2:
            if diffusion_limit is None or not diffusion_source:
                raise ValueError("Second-order rates require an evidenced diffusion ceiling")
            diffusion_limit = _finite(diffusion_limit, name="diffusion_limit_L_mol_s")
            if diffusion_limit <= 0.0 or rate * upper > diffusion_limit * (1.0 + 1.0e-12):
                raise ValueError("Second-order rate uncertainty exceeds its diffusion ceiling")
        elif diffusion_limit is not None or diffusion_source is not None:
            raise ValueError("Diffusion ceilings apply only to second-order reactions")
        flags = {str(key): bool(value) for key, value in dict(self.applicability_flags).items()}
        unknown = set(flags) - PROHIBITED_MECHANISM_FLAGS
        if unknown:
            raise ValueError("Unknown solution applicability flags: " + ", ".join(sorted(unknown)))
        partner = str(self.reversible_partner_id).strip() if self.reversible_partner_id else None
        object.__setattr__(self, "reaction_id", reaction_id)
        object.__setattr__(self, "reactants", dict(sorted(reactants.items())))
        object.__setattr__(self, "products", dict(sorted(products.items())))
        object.__setattr__(self, "rate_constant", rate)
        object.__setattr__(self, "rate_units", "s^-1" if order == 1 else "L mol^-1 s^-1")
        object.__setattr__(self, "temperature_K", _finite(self.temperature_K, name="temperature_K"))
        object.__setattr__(self, "solvent", _text(self.solvent, name="solvent"))
        object.__setattr__(self, "solvent_model", _text(self.solvent_model, name="solvent_model").upper())
        object.__setattr__(self, "standard_state_mol_L", _finite(self.standard_state_mol_L, name="standard_state_mol_L"))
        object.__setattr__(self, "verification_status", status)
        object.__setattr__(self, "rate_theory", theory)
        object.__setattr__(self, "source", _text(self.source, name="source"))
        object.__setattr__(self, "rate_claim_id", _text(self.rate_claim_id, name="rate_claim_id"))
        object.__setattr__(self, "path_claim_id", path_claim)
        object.__setattr__(self, "standard_state_correction", correction)
        object.__setattr__(self, "calculation_hashes", hashes)
        object.__setattr__(self, "evidence_hashes", evidence_hashes)
        object.__setattr__(self, "lower_rate_multiplier", lower)
        object.__setattr__(self, "upper_rate_multiplier", upper)
        object.__setattr__(self, "uncertainty_source", uncertainty_source)
        object.__setattr__(self, "uncertainty_components", uncertainty_components)
        object.__setattr__(self, "diffusion_limit_L_mol_s", diffusion_limit)
        object.__setattr__(self, "diffusion_limit_source", diffusion_source)
        object.__setattr__(self, "reversible_partner_id", partner)
        object.__setattr__(self, "applicability_flags", dict(sorted(flags.items())))

    @property
    def reaction_order(self) -> int:
        return sum(self.reactants.values())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "VerifiedMassActionReaction":
        if not isinstance(value, Mapping):
            raise TypeError("mass-action reaction must be a mapping")
        rate = value.get("rate_constant")
        if rate is None:
            rate = value.get("rate_constant_s")
        return cls(
            reaction_id=value.get("reaction_id"), reactants=value.get("reactants") or {},
            products=value.get("products") or {}, rate_constant=rate,
            rate_units=value.get("rate_units"), temperature_K=value.get("temperature_K"),
            solvent=value.get("solvent"), solvent_model=value.get("solvent_model"),
            standard_state_mol_L=value.get("standard_state_mol_L", 1.0),
            verification_status=value.get("verification_status"), rate_theory=value.get("rate_theory"),
            source=value.get("source"), rate_claim_id=value.get("rate_claim_id"),
            path_claim_id=value.get("path_claim_id"),
            standard_state_correction=value.get("standard_state_correction"),
            calculation_hashes=tuple(value.get("calculation_hashes") or ()),
            evidence_hashes=tuple(value.get("evidence_hashes") or ()),
            lower_rate_multiplier=value.get("lower_rate_multiplier"),
            upper_rate_multiplier=value.get("upper_rate_multiplier"),
            uncertainty_source=value.get("uncertainty_source"),
            uncertainty_components=tuple(value.get("uncertainty_components") or ()),
            diffusion_limit_L_mol_s=value.get("diffusion_limit_L_mol_s"),
            diffusion_limit_source=value.get("diffusion_limit_source"),
            reversible_partner_id=value.get("reversible_partner_id"),
            applicability_flags=value.get("applicability_flags") or {},
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "reaction_id": self.reaction_id, "reactants": dict(self.reactants),
            "products": dict(self.products), "reaction_order": self.reaction_order,
            "rate_constant": self.rate_constant, "rate_units": self.rate_units,
            "temperature_K": self.temperature_K, "solvent": self.solvent,
            "solvent_model": self.solvent_model, "standard_state_mol_L": self.standard_state_mol_L,
            "verification_status": self.verification_status, "rate_theory": self.rate_theory,
            "source": self.source, "rate_claim_id": self.rate_claim_id,
            "path_claim_id": self.path_claim_id,
            "standard_state_correction": self.standard_state_correction,
            "calculation_hashes": list(self.calculation_hashes),
            "evidence_hashes": list(self.evidence_hashes),
            "lower_rate_multiplier": self.lower_rate_multiplier,
            "upper_rate_multiplier": self.upper_rate_multiplier,
            "uncertainty_source": self.uncertainty_source,
            "uncertainty_components": list(self.uncertainty_components),
            "diffusion_limit_L_mol_s": self.diffusion_limit_L_mol_s,
            "diffusion_limit_source": self.diffusion_limit_source,
            "reversible_partner_id": self.reversible_partner_id,
            "applicability_flags": dict(self.applicability_flags),
        }


@dataclass(frozen=True)
class VerifiedMassActionNetwork:
    target_species_id: str
    species: tuple[SolutionSpecies, ...]
    reactions: tuple[VerifiedMassActionReaction, ...]
    coverage_status: str
    coverage_claim_id: str | None = None
    discovery_sources: tuple[str, ...] = ()
    unresolved_candidate_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        target = _text(self.target_species_id, name="target_species_id")
        if not self.species or not self.reactions:
            raise ValueError("A mass-action network needs species and reactions")
        species_ids = [item.species_id for item in self.species]
        reaction_ids = [item.reaction_id for item in self.reactions]
        if len(species_ids) != len(set(species_ids)) or len(reaction_ids) != len(set(reaction_ids)):
            raise ValueError("Mass-action species and reaction IDs must be unique")
        if len({item.smiles for item in self.species}) != len(self.species):
            raise ValueError("Mass-action species must have unique canonical identities")
        if target not in species_ids:
            raise ValueError("Target species is absent from the mass-action network")
        by_id = {item.species_id: item for item in self.species}
        if by_id[target].initial_concentration_mol_L <= 0.0:
            raise ValueError("Target initial concentration must be positive")
        available = set(species_ids)
        for reaction in self.reactions:
            referenced = set(reaction.reactants) | set(reaction.products)
            if not referenced <= available:
                raise ValueError(f"Reaction {reaction.reaction_id} references unknown species")
            self._validate_balance(reaction, by_id)
        partners = {item.reaction_id: item for item in self.reactions}
        for reaction in self.reactions:
            if reaction.reversible_partner_id:
                partner = partners.get(reaction.reversible_partner_id)
                if (
                    partner is None or partner.reversible_partner_id != reaction.reaction_id
                    or dict(partner.reactants) != dict(reaction.products)
                    or dict(partner.products) != dict(reaction.reactants)
                ):
                    raise ValueError(f"Reaction {reaction.reaction_id} has an invalid reversible partner")
        coverage = _text(self.coverage_status, name="coverage_status")
        claim = str(self.coverage_claim_id).strip() if self.coverage_claim_id else None
        sources = tuple(_text(item, name="discovery_source") for item in self.discovery_sources)
        unresolved = tuple(_text(item, name="unresolved_candidate_id") for item in self.unresolved_candidate_ids)
        if coverage == MASS_ACTION_COVERAGE_STATUS and (not claim or not sources or unresolved):
            raise ValueError("Verified mass-action coverage requires a claim, sources, and no unresolved candidates")
        if coverage not in {MASS_ACTION_COVERAGE_STATUS, "bounded_search_incomplete", "unresolved"}:
            raise ValueError("Unsupported mass-action network coverage status")
        object.__setattr__(self, "target_species_id", target)
        object.__setattr__(self, "coverage_status", coverage)
        object.__setattr__(self, "coverage_claim_id", claim)
        object.__setattr__(self, "discovery_sources", sources)
        object.__setattr__(self, "unresolved_candidate_ids", unresolved)

    @staticmethod
    def _validate_balance(reaction: VerifiedMassActionReaction, by_id: Mapping[str, SolutionSpecies]) -> None:
        elements = {element for species in by_id.values() for element in species.element_counts}
        mismatches = {}
        for element in sorted(elements):
            left = sum(coef * by_id[item].element_counts.get(element, 0) for item, coef in reaction.reactants.items())
            right = sum(coef * by_id[item].element_counts.get(element, 0) for item, coef in reaction.products.items())
            if left != right:
                mismatches[element] = {"reactants": left, "products": right}
        left_charge = sum(coef * by_id[item].charge for item, coef in reaction.reactants.items())
        right_charge = sum(coef * by_id[item].charge for item, coef in reaction.products.items())
        if mismatches or left_charge != right_charge:
            raise ValueError(
                f"Reaction {reaction.reaction_id} is not element/charge balanced: "
                f"elements={mismatches}, charge={left_charge}->{right_charge}"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "VerifiedMassActionNetwork":
        if not isinstance(value, Mapping):
            raise TypeError("mass-action network must be a mapping")
        return cls(
            target_species_id=value.get("target_species_id"),
            species=tuple(SolutionSpecies.from_dict(item) for item in value.get("species") or ()),
            reactions=tuple(VerifiedMassActionReaction.from_dict(item) for item in value.get("reactions") or ()),
            coverage_status=value.get("coverage_status", "unresolved"),
            coverage_claim_id=value.get("coverage_claim_id"),
            discovery_sources=tuple(value.get("discovery_sources") or ()),
            unresolved_candidate_ids=tuple(value.get("unresolved_candidate_ids") or ()),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": MASS_ACTION_NETWORK_SCHEMA_VERSION,
            "kind": "verified_implicit_solution_mass_action_network",
            "target_species_id": self.target_species_id,
            "coverage_status": self.coverage_status,
            "coverage_claim_id": self.coverage_claim_id,
            "discovery_sources": list(self.discovery_sources),
            "unresolved_candidate_ids": list(self.unresolved_candidate_ids),
            "species": [item.as_dict() for item in self.species],
            "reactions": [item.as_dict() for item in self.reactions],
        }


def assess_mass_action_applicability(
    environment: MassActionSolutionEnvironment,
    network: VerifiedMassActionNetwork,
) -> dict[str, Any]:
    blockers: list[dict[str, Any]] = []
    if network.coverage_status != MASS_ACTION_COVERAGE_STATUS:
        blockers.append({
            "reason": "reaction_network_coverage_unresolved",
            "coverage_status": network.coverage_status,
            "unresolved_candidate_ids": list(network.unresolved_candidate_ids),
        })
    if not math.isclose(environment.standard_state_mol_L, 1.0, rel_tol=0.0, abs_tol=1e-12):
        blockers.append({"reason": "unsupported_standard_state", "required_mol_L": 1.0})
    ionic = [item.species_id for item in network.species if item.charge != 0]
    if ionic:
        blockers.append({"reason": "ionic_species_outside_neutral_lane", "species_ids": ionic})
    for reaction in network.reactions:
        reasons = []
        if not math.isclose(reaction.temperature_K, environment.temperature_K, rel_tol=0.0, abs_tol=1e-8):
            reasons.append("rate_temperature_mismatch")
        if reaction.solvent.casefold() != environment.solvent.casefold():
            reasons.append("rate_solvent_mismatch")
        if reaction.solvent_model != environment.solvent_model:
            reasons.append("rate_solvent_model_mismatch")
        if not math.isclose(reaction.standard_state_mol_L, environment.standard_state_mol_L, rel_tol=0.0, abs_tol=1e-12):
            reasons.append("rate_standard_state_mismatch")
        reasons.extend(sorted(key for key, active in reaction.applicability_flags.items() if active))
        if reasons:
            blockers.append({"reaction_id": reaction.reaction_id, "reasons": reasons})
    return {
        "schema_version": 1,
        "kind": "implicit_solution_mass_action_applicability_assessment",
        "status": "supported" if not blockers else "outside_method_applicability",
        "blocking_findings": blockers,
        "scope": (
            "Dark neutral homogeneous constant-volume implicit solution; ideal molar activities; "
            "verified first/second-order rates; evidenced diffusion ceilings."
        ),
    }


def _system_matrices(network: VerifiedMassActionNetwork) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    index = {item.species_id: position for position, item in enumerate(network.species)}
    stoichiometry = np.zeros((len(network.species), len(network.reactions)))
    reactant_orders = np.zeros((len(network.reactions), len(network.species)))
    for column, reaction in enumerate(network.reactions):
        for species_id, coefficient in reaction.reactants.items():
            stoichiometry[index[species_id], column] -= coefficient
            reactant_orders[column, index[species_id]] = coefficient
        for species_id, coefficient in reaction.products.items():
            stoichiometry[index[species_id], column] += coefficient
    return stoichiometry, reactant_orders, np.asarray(
        [item.initial_concentration_mol_L for item in network.species], dtype=float,
    )


def _integrate(
    network: VerifiedMassActionNetwork,
    rate_constants: np.ndarray,
    times: np.ndarray,
    *,
    target_index: int,
    threshold: float,
    relative_tolerance: float,
    absolute_tolerance: float,
    maximum_steps: int,
) -> dict[str, Any]:
    stoich, orders, initial_concentrations = _system_matrices(network)
    species_count = len(network.species)
    y = np.concatenate((initial_concentrations, np.zeros(len(network.reactions))))

    def rates_and_derivatives(state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        concentrations = np.maximum(state[:species_count], 0.0)
        rates = rate_constants.copy()
        for reaction_index in range(len(network.reactions)):
            for species_index, order in enumerate(orders[reaction_index]):
                if order:
                    rates[reaction_index] *= concentrations[species_index] ** order
        rate_jacobian = np.zeros((len(network.reactions), species_count))
        for reaction_index in range(len(network.reactions)):
            for differentiated_index, differentiated_order in enumerate(orders[reaction_index]):
                if not differentiated_order:
                    continue
                value = rate_constants[reaction_index] * differentiated_order
                for species_index, order in enumerate(orders[reaction_index]):
                    exponent = order - (1 if species_index == differentiated_index else 0)
                    if exponent:
                        value *= concentrations[species_index] ** exponent
                rate_jacobian[reaction_index, differentiated_index] = value
        return rates, rate_jacobian

    def derivative(state: np.ndarray) -> np.ndarray:
        rates, _ = rates_and_derivatives(state)
        return np.concatenate((stoich @ rates, rates))

    def jacobian(state: np.ndarray) -> np.ndarray:
        _, rate_jacobian = rates_and_derivatives(state)
        result = np.zeros((len(state), len(state)))
        result[:species_count, :species_count] = stoich @ rate_jacobian
        result[species_count:, :species_count] = rate_jacobian
        return result

    def rk4(state: np.ndarray, step: float) -> np.ndarray:
        k1 = derivative(state)
        k2 = derivative(state + 0.5 * step * k1)
        k3 = derivative(state + 0.5 * step * k2)
        k4 = derivative(state + step * k3)
        return state + step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    maximum_linear_condition = 1.0

    def linearly_implicit_euler(state: np.ndarray, step: float) -> np.ndarray:
        nonlocal maximum_linear_condition
        matrix = np.eye(len(state)) - step * jacobian(state)
        condition = float(np.linalg.cond(matrix))
        maximum_linear_condition = max(maximum_linear_condition, condition)
        if not math.isfinite(condition) or condition > 1.0e14:
            raise ValueError("Stiff mass-action linear solve is ill-conditioned")
        return state + np.linalg.solve(matrix, step * derivative(state))

    initial_jacobian = jacobian(y)[:species_count, :species_count]
    eigenvalues = np.linalg.eigvals(initial_jacobian)
    spectral_rate = float(max((abs(value) for value in eigenvalues), default=0.0))
    stiffness_indicator = spectral_rate * float(times[-1])
    use_stiff_method = stiffness_indicator > 50.0
    advance = linearly_implicit_euler if use_stiff_method else rk4
    method_order = 1 if use_stiff_method else 4

    def crossing_in_step(start: np.ndarray, step: float) -> float:
        lower, upper = 0.0, step
        for _ in range(48):
            midpoint = 0.5 * (lower + upper)
            if advance(start, midpoint)[target_index] <= threshold:
                upper = midpoint
            else:
                lower = midpoint
        return upper

    output_states = [y.copy()]
    current_time = 0.0
    next_output = 1
    horizon = float(times[-1])
    step = max(horizon * 1.0e-8, min(horizon / 100.0, 1.0))
    accepted = rejected = 0
    minimum_step = math.inf
    maximum_step = 0.0
    most_negative = 0.0
    first_crossing = 0.0 if y[target_index] <= threshold else None
    while current_time < horizon:
        if accepted + rejected >= maximum_steps:
            raise ValueError("Nonlinear solution propagation exceeded maximum_steps")
        boundary = float(times[next_output]) if next_output < len(times) else horizon
        step = min(step, boundary - current_time, horizon - current_time)
        if step <= max(1.0, horizon) * 1.0e-15:
            raise ValueError("Nonlinear solution propagation step underflow")
        full = advance(y, step)
        half = advance(advance(y, 0.5 * step), 0.5 * step)
        scale = absolute_tolerance + relative_tolerance * np.maximum(np.abs(y), np.abs(half))
        error = float(np.max(np.abs(half - full) / np.maximum(scale, 1.0e-300)))
        concentration_floor = float(np.min(half[:species_count]))
        negative_tolerance = max(absolute_tolerance * 10.0, float(np.max(initial_concentrations)) * 1.0e-10)
        if error > 1.0 or concentration_floor < -negative_tolerance:
            rejected += 1
            step *= max(0.1, 0.9 * max(error, 1.0) ** (-1.0 / (method_order + 1)))
            continue
        previous = y
        previous_time = current_time
        y = half
        most_negative = min(most_negative, concentration_floor)
        y[:species_count][np.abs(y[:species_count]) < negative_tolerance] = 0.0
        current_time += step
        accepted += 1
        minimum_step = min(minimum_step, step)
        maximum_step = max(maximum_step, step)
        if first_crossing is None and previous[target_index] > threshold and y[target_index] <= threshold:
            first_crossing = previous_time + crossing_in_step(previous, step)
        if next_output < len(times) and math.isclose(current_time, boundary, rel_tol=0.0, abs_tol=max(1.0, horizon) * 1.0e-13):
            output_states.append(y.copy())
            next_output += 1
        factor = (
            2.0 if error == 0.0
            else min(2.5, max(0.2, 0.9 * error ** (-1.0 / (method_order + 1))))
        )
        step *= factor
    if len(output_states) != len(times):
        raise ValueError("Nonlinear integrator failed to retain every requested profile time")
    return {
        "states": output_states,
        "first_crossing_seconds": first_crossing,
        "accepted_steps": accepted,
        "rejected_steps": rejected,
        "minimum_step_seconds": minimum_step if accepted else None,
        "maximum_step_seconds": maximum_step,
        "most_negative_raw_concentration_mol_L": most_negative,
        "method": (
            "adaptive_linearly_implicit_euler_step_doubling"
            if use_stiff_method else "adaptive_rk4_step_doubling"
        ),
        "initial_spectral_stiffness_indicator": stiffness_indicator,
        "maximum_linear_solve_condition_number": maximum_linear_condition,
    }


def _element_inventory(network: VerifiedMassActionNetwork, concentrations: np.ndarray) -> dict[str, float]:
    elements = sorted({element for item in network.species for element in item.element_counts})
    return {
        element: float(sum(
            concentrations[index] * species.element_counts.get(element, 0)
            for index, species in enumerate(network.species)
        ))
        for element in elements
    }


def propagate_verified_mass_action_network(
    environment: MassActionSolutionEnvironment,
    network: VerifiedMassActionNetwork,
    *,
    profile_points: int = 180,
    relative_tolerance: float = 1.0e-7,
    absolute_tolerance: float = 1.0e-12,
    maximum_steps: int = 200_000,
) -> dict[str, Any]:
    applicability = assess_mass_action_applicability(environment, network)
    if applicability["status"] != "supported":
        return {
            "schema_version": 1, "kind": "verified_solution_mass_action_propagation",
            "status": "outside_method_applicability", "applicability": applicability,
            "estimated_time_to_retention_seconds": None,
        }
    if profile_points < 20 or profile_points > 1000:
        raise ValueError("profile_points must be in [20, 1000]")
    if (
        isinstance(maximum_steps, bool)
        or int(maximum_steps) != maximum_steps
        or maximum_steps < 100
    ):
        raise ValueError("maximum_steps must be an integer of at least 100")
    if not (0.0 < relative_tolerance <= 1.0e-3 and 0.0 < absolute_tolerance <= 1.0e-6):
        raise ValueError("Integrator tolerances are outside the audited range")
    horizon = environment.target_duration_seconds
    earliest = max(horizon * 1.0e-10, min(horizon, 1.0e-6))
    times = (
        np.asarray([0.0, horizon]) if earliest >= horizon
        else np.concatenate(([0.0], np.geomspace(earliest, horizon, profile_points - 1)))
    )
    species_index = {item.species_id: index for index, item in enumerate(network.species)}
    target_index = species_index[network.target_species_id]
    initial_target = network.species[target_index].initial_concentration_mol_L
    threshold = initial_target * environment.retention_fraction
    central_rates = np.asarray([item.rate_constant for item in network.reactions], dtype=float)
    central = _integrate(
        network, central_rates, times, target_index=target_index, threshold=threshold,
        relative_tolerance=relative_tolerance, absolute_tolerance=absolute_tolerance,
        maximum_steps=maximum_steps,
    )
    final_state = central["states"][-1]
    species_count = len(network.species)
    initial_concentrations = np.asarray(
        [item.initial_concentration_mol_L for item in network.species], dtype=float,
    )
    final_concentrations = final_state[:species_count]
    extents = final_state[species_count:]

    initial_inventory = _element_inventory(network, initial_concentrations)
    final_inventory = _element_inventory(network, final_concentrations)
    inventory_errors = {
        element: abs(final_inventory[element] - initial_inventory[element]) / max(abs(initial_inventory[element]), 1.0e-300)
        for element in initial_inventory
    }
    maximum_inventory_error = max(inventory_errors.values(), default=0.0)
    if maximum_inventory_error > 1.0e-7:
        raise ValueError("Mass-action propagation failed elemental inventory closure")
    initial_charge_inventory = float(sum(
        concentration * species.charge
        for concentration, species in zip(initial_concentrations, network.species)
    ))
    final_charge_inventory = float(sum(
        concentration * species.charge
        for concentration, species in zip(final_concentrations, network.species)
    ))
    charge_inventory_error = abs(final_charge_inventory - initial_charge_inventory)
    if charge_inventory_error > max(1.0, abs(initial_charge_inventory)) * 1.0e-10:
        raise ValueError("Mass-action propagation failed charge inventory closure")
    stoich, _, _ = _system_matrices(network)
    extent_delta = stoich @ extents
    closure = float(np.max(np.abs((final_concentrations - initial_concentrations) - extent_delta)))
    closure_scale = max(1.0, float(np.max(initial_concentrations)))
    if closure / closure_scale > 1.0e-8:
        raise ValueError("Mass-action propagation failed reaction-extent closure")

    def scenario_observables(state: np.ndarray) -> dict[str, Any]:
        concentrations = state[:species_count]
        scenario_extents = state[species_count:]
        gains = {
            item.species_id: max(
                0.0, float(concentrations[index] - initial_concentrations[index]),
            )
            for index, item in enumerate(network.species)
            if item.species_id != network.target_species_id
        }
        dominant = max(gains, key=gains.get) if any(gains.values()) else None
        consumptions = [
            max(0.0, -float(stoich[target_index, index])) * float(scenario_extents[index])
            for index in range(len(network.reactions))
        ]
        total_consumption = sum(consumptions)
        return {
            "final_concentrations_mol_L": {
                item.species_id: float(concentrations[index])
                for index, item in enumerate(network.species)
            },
            "product_yields_relative_to_initial_target": {
                key: value / initial_target for key, value in gains.items()
            },
            "dominant_product_species_id": dominant,
            "integrated_reaction_extents_mol_L": {
                reaction.reaction_id: float(scenario_extents[index])
                for index, reaction in enumerate(network.reactions)
            },
            "gross_target_consumption_fractions": {
                reaction.reaction_id: (
                    consumptions[index] / total_consumption if total_consumption > 0.0 else 0.0
                )
                for index, reaction in enumerate(network.reactions)
            },
        }

    # Deterministic uncertainty sensitivity set.  For small networks retain all
    # corners; for larger ones retain joint and one-at-a-time extremes and say
    # explicitly that this is not a rigorous continuous interval enclosure.
    multipliers = [
        (float(item.lower_rate_multiplier), float(item.upper_rate_multiplier))
        for item in network.reactions
    ]
    if len(network.reactions) <= 10:
        scenario_vectors = list(itertools.product((0, 1), repeat=len(network.reactions)))
        scenario_kind = "all_rate_bound_corners"
    else:
        all_low = tuple(0 for _ in network.reactions)
        all_high = tuple(1 for _ in network.reactions)
        scenario_vectors = [all_low, all_high]
        for index in range(len(network.reactions)):
            low = [None] * len(network.reactions)
            high = [None] * len(network.reactions)
            low[index] = 0
            high[index] = 1
            scenario_vectors.extend((tuple(low), tuple(high)))
        scenario_kind = "joint_and_one_at_a_time_rate_bounds"
    uncertainty_scenarios = []
    for scenario_index, vector in enumerate(scenario_vectors):
        rates = central_rates.copy()
        labels = []
        for index, selector in enumerate(vector):
            if selector is None:
                labels.append("central")
            else:
                rates[index] *= multipliers[index][selector]
                labels.append("lower" if selector == 0 else "upper")
        scenario = _integrate(
            network, rates, np.asarray([0.0, horizon]), target_index=target_index,
            threshold=threshold, relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance, maximum_steps=maximum_steps,
        )
        scenario_metrics = scenario_observables(scenario["states"][-1])
        uncertainty_scenarios.append({
            "scenario_id": f"rate-bounds-{scenario_index:04d}",
            "rate_settings": labels,
            "first_crossing_seconds": scenario["first_crossing_seconds"],
            "target_fraction_remaining": float(scenario["states"][-1][target_index] / initial_target),
            **scenario_metrics,
        })
    crossings = [item["first_crossing_seconds"] for item in uncertainty_scenarios]
    all_cross = all(value is not None for value in crossings)
    none_cross = all(value is None for value in crossings)
    uncertainty_status = (
        "crossing_within_horizon_for_all_retained_scenarios" if all_cross
        else "no_crossing_within_horizon_for_any_retained_scenario" if none_cross
        else "crossing_depends_on_rate_uncertainty"
    )

    profile = []
    for time, state in zip(times, central["states"]):
        concentrations = state[:species_count]
        profile.append({
            "time_seconds": float(time),
            "target_fraction_remaining": float(concentrations[target_index] / initial_target),
            "concentrations_mol_L": {
                item.species_id: float(concentrations[index])
                for index, item in enumerate(network.species)
            },
        })
    net_target_contributions = []
    gross_target_consumption = 0.0
    for index, reaction in enumerate(network.reactions):
        target_coefficient = float(stoich[target_index, index])
        target_consumed = max(0.0, -target_coefficient) * float(extents[index])
        target_formed = max(0.0, target_coefficient) * float(extents[index])
        gross_target_consumption += target_consumed
        net_target_contributions.append({
            "reaction_id": reaction.reaction_id,
            "integrated_extent_mol_L": float(extents[index]),
            "net_target_contribution_mol_L": float(target_coefficient * extents[index]),
            "target_stoichiometric_coefficient": target_coefficient,
            "integrated_target_consumption_mol_L": target_consumed,
            "integrated_target_formation_mol_L": target_formed,
            "rate_claim_id": reaction.rate_claim_id,
            "path_claim_id": reaction.path_claim_id,
        })
    for item in net_target_contributions:
        item["gross_target_consumption_fraction"] = (
            item["integrated_target_consumption_mol_L"] / gross_target_consumption
            if gross_target_consumption > 0.0 else 0.0
        )
    dominant_target_consumption = max(
        net_target_contributions,
        key=lambda item: item["integrated_target_consumption_mol_L"],
    )
    product_gains = {
        item.species_id: max(0.0, float(final_concentrations[index] - initial_concentrations[index]))
        for index, item in enumerate(network.species)
        if item.species_id != network.target_species_id
    }
    dominant_product = max(product_gains, key=product_gains.get) if any(product_gains.values()) else None
    t_values = [float(value) for value in crossings if value is not None]
    return {
        "schema_version": 1,
        "kind": "verified_solution_mass_action_propagation",
        "status": "completed",
        "applicability": applicability,
        "environment_contract": environment.as_environment_spec().as_dict(),
        "network": network.as_dict(),
        "estimated_time_to_retention_seconds": central["first_crossing_seconds"],
        "retention_fraction": environment.retention_fraction,
        "target_duration_seconds": horizon,
        "target_species_id": network.target_species_id,
        "target_fraction_remaining": float(final_concentrations[target_index] / initial_target),
        "final_concentrations_mol_L": {
            item.species_id: float(final_concentrations[index])
            for index, item in enumerate(network.species)
        },
        "product_yields_relative_to_initial_target": {
            key: value / initial_target for key, value in product_gains.items()
        },
        "dominant_product_species_id": dominant_product,
        "reaction_flux_attribution": {
            "status": "completed", "reactions": net_target_contributions,
            "gross_target_consumption_mol_L": gross_target_consumption,
            "dominant_target_consumption_reaction_id": (
                dominant_target_consumption["reaction_id"]
                if gross_target_consumption > 0.0 else None
            ),
            "target_net_change_closure_relative_error": abs(
                sum(item["net_target_contribution_mol_L"] for item in net_target_contributions)
                - (final_concentrations[target_index] - initial_target)
            ) / max(initial_target, 1.0e-300),
        },
        "rate_uncertainty_sensitivity": {
            "status": uncertainty_status,
            "scenario_kind": scenario_kind,
            "scenario_count": len(uncertainty_scenarios),
            "lower_crossing_seconds": min(t_values) if all_cross else None,
            "upper_crossing_seconds": max(t_values) if all_cross else None,
            "scenarios": uncertainty_scenarios,
            "interpretation": (
                "Deterministic propagation of retained rate-bound scenarios; not a statistical confidence "
                "interval and, for more than ten reactions, not a rigorous enclosure of every joint rate combination."
            ),
        },
        "numerical_validation": {
            "status": "passed",
            "method": central["method"],
            "method_selection": {
                "policy": "initial_jacobian_spectral_rate_times_horizon",
                "stiff_threshold": 50.0,
                "indicator": central["initial_spectral_stiffness_indicator"],
            },
            "relative_tolerance": relative_tolerance,
            "absolute_tolerance_mol_L": absolute_tolerance,
            "accepted_steps": central["accepted_steps"],
            "rejected_steps": central["rejected_steps"],
            "minimum_step_seconds": central["minimum_step_seconds"],
            "maximum_step_seconds": central["maximum_step_seconds"],
            "most_negative_raw_concentration_mol_L": central["most_negative_raw_concentration_mol_L"],
            "maximum_element_inventory_relative_error": maximum_inventory_error,
            "element_inventory_relative_errors": inventory_errors,
            "initial_charge_inventory_mol_L": initial_charge_inventory,
            "final_charge_inventory_mol_L": final_charge_inventory,
            "charge_inventory_absolute_error_mol_L": charge_inventory_error,
            "reaction_extent_closure_relative_error": closure / closure_scale,
            "maximum_linear_solve_condition_number": central["maximum_linear_solve_condition_number"],
        },
        "target_profile": profile,
        "claim_scope": (
            "Conditional on the supplied verified solution rates, ideal-activity mass action, and a "
            "complete network claim within the declared neutral homogeneous scope."
        ),
    }


__all__ = [
    "MASS_ACTION_COVERAGE_STATUS", "MassActionSolutionEnvironment",
    "VerifiedMassActionNetwork", "VerifiedMassActionReaction",
    "assess_mass_action_applicability", "propagate_verified_mass_action_network",
]
