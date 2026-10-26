"""Verified-input concentration kinetics for the first solution stability lane.

The implemented domain is intentionally narrow: dark, neutral, homogeneous,
constant-volume solution chemistry represented entirely by verified first-order
elementary rates.  The module does not convert gas rates, infer solvent effects,
or approximate bimolecular encounter kinetics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

import numpy as np


SOLUTION_NETWORK_SCHEMA_VERSION = 1
SUPPORTED_SOLVENT_MODELS = {"CPCM", "SMD"}
SUPPORTED_RATE_THEORIES = {"experimental_solution_rate", "solution_eyring_tst"}
SUPPORTED_RATE_STATUSES = {"verified_experimental_solution_rate", "verified_solution_tst"}
PROHIBITED_MECHANISM_FLAGS = {
    "bimolecular_encounter",
    "diffusion_limited",
    "explicit_solvent_participation",
    "acid_base_speciation_required",
    "ion_pairing_required",
    "photochemical",
    "spin_crossing",
    "surface_or_interface",
    "multireference_unresolved",
}


def _finite(value: Any, *, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite number")
    return number


def _text(value: Any, *, name: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"{name} is required")
    return text


@dataclass(frozen=True)
class ImplicitSolutionEnvironment:
    solvent: str
    solvent_model: str
    temperature_K: float
    target_duration_seconds: float
    retention_fraction: float = 0.95
    standard_state_mol_L: float = 1.0
    volume_L: float = 1.0
    light_condition: str = "dark"
    electronic_structure_profile: str = "user_verified_solution_profile"

    def __post_init__(self) -> None:
        solvent = _text(self.solvent, name="solvent")
        model = _text(self.solvent_model, name="solvent_model").upper()
        if model not in SUPPORTED_SOLVENT_MODELS:
            raise ValueError("solvent_model must be CPCM or SMD")
        temperature = _finite(self.temperature_K, name="temperature_K")
        duration = _finite(self.target_duration_seconds, name="target_duration_seconds")
        retention = _finite(self.retention_fraction, name="retention_fraction")
        standard_state = _finite(self.standard_state_mol_L, name="standard_state_mol_L")
        volume = _finite(self.volume_L, name="volume_L")
        if temperature <= 0.0 or duration <= 0.0 or standard_state <= 0.0 or volume <= 0.0:
            raise ValueError("Solution temperature, duration, standard state, and volume must be positive")
        if not 0.0 < retention < 1.0:
            raise ValueError("retention_fraction must be in (0, 1)")
        if str(self.light_condition).strip().lower() != "dark":
            raise ValueError("The first implicit-solution lane supports dark conditions only")
        object.__setattr__(self, "solvent", solvent)
        object.__setattr__(self, "solvent_model", model)
        object.__setattr__(self, "temperature_K", temperature)
        object.__setattr__(self, "target_duration_seconds", duration)
        object.__setattr__(self, "retention_fraction", retention)
        object.__setattr__(self, "standard_state_mol_L", standard_state)
        object.__setattr__(self, "volume_L", volume)
        object.__setattr__(self, "light_condition", "dark")
        object.__setattr__(
            self, "electronic_structure_profile",
            _text(self.electronic_structure_profile, name="electronic_structure_profile"),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ImplicitSolutionEnvironment":
        if not isinstance(value, Mapping):
            raise TypeError("solution environment must be a mapping")
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
        )

    def as_environment_spec(self):
        from .lanes.base import EnvironmentSpec

        return EnvironmentSpec(
            environment_kind="homogeneous_solution",
            model_name="implicit_solvent_constant_volume_unimolecular_v1",
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
                "phase_model": "homogeneous implicit solution",
            },
        )


@dataclass(frozen=True)
class SolutionSpecies:
    species_id: str
    smiles: str
    initial_concentration_mol_L: float
    element_counts: Mapping[str, int]
    charge: int = 0

    def __post_init__(self) -> None:
        species_id = _text(self.species_id, name="species_id")
        smiles = _text(self.smiles, name=f"{species_id}.smiles")
        concentration = _finite(
            self.initial_concentration_mol_L, name=f"{species_id}.initial_concentration_mol_L",
        )
        if concentration < 0.0:
            raise ValueError("Initial concentrations must be nonnegative")
        elements = {str(key): int(value) for key, value in dict(self.element_counts).items()}
        if not elements or any(not key or value <= 0 for key, value in elements.items()):
            raise ValueError(f"Species {species_id} needs positive elemental counts")
        try:
            from rdkit import Chem
        except ImportError as error:
            raise RuntimeError("RDKit is required to validate solution species identity") from error
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Species {species_id} has invalid SMILES")
        explicit = Chem.AddHs(molecule)
        derived_elements: dict[str, int] = {}
        for atom in explicit.GetAtoms():
            symbol = atom.GetSymbol()
            derived_elements[symbol] = derived_elements.get(symbol, 0) + 1
        derived_charge = int(Chem.GetFormalCharge(molecule))
        if dict(sorted(derived_elements.items())) != dict(sorted(elements.items())):
            raise ValueError(
                f"Species {species_id} elemental counts do not match its SMILES: "
                f"declared={elements}, derived={derived_elements}"
            )
        if derived_charge != int(self.charge):
            raise ValueError(
                f"Species {species_id} charge does not match its SMILES: "
                f"declared={self.charge}, derived={derived_charge}"
            )
        object.__setattr__(self, "species_id", species_id)
        object.__setattr__(self, "smiles", Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True))
        object.__setattr__(self, "initial_concentration_mol_L", concentration)
        object.__setattr__(self, "element_counts", dict(sorted(elements.items())))
        object.__setattr__(self, "charge", int(self.charge))

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SolutionSpecies":
        return cls(
            species_id=value.get("species_id"),
            smiles=value.get("smiles"),
            initial_concentration_mol_L=value.get("initial_concentration_mol_L"),
            element_counts=value.get("element_counts") or {},
            charge=value.get("charge", 0),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "species_id": self.species_id,
            "smiles": self.smiles,
            "initial_concentration_mol_L": self.initial_concentration_mol_L,
            "element_counts": dict(self.element_counts),
            "charge": self.charge,
        }


@dataclass(frozen=True)
class VerifiedSolutionReaction:
    reaction_id: str
    reactant_id: str
    products: Mapping[str, int]
    rate_constant_s: float
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
    lower_rate_multiplier: float | None = None
    upper_rate_multiplier: float | None = None
    uncertainty_source: str | None = None
    uncertainty_components: tuple[str, ...] = ()
    applicability_flags: Mapping[str, bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        reaction_id = _text(self.reaction_id, name="reaction_id")
        reactant = _text(self.reactant_id, name=f"{reaction_id}.reactant_id")
        products = {str(key): int(value) for key, value in dict(self.products).items()}
        if not products or any(not key or value <= 0 for key, value in products.items()):
            raise ValueError(f"Reaction {reaction_id} needs positive integer product stoichiometry")
        rate = _finite(self.rate_constant_s, name=f"{reaction_id}.rate_constant_s")
        if rate <= 0.0:
            raise ValueError("Verified solution rates must be positive")
        status = _text(self.verification_status, name="verification_status")
        theory = _text(self.rate_theory, name="rate_theory")
        if status not in SUPPORTED_RATE_STATUSES:
            raise ValueError(f"Unsupported solution rate verification status: {status}")
        if theory not in SUPPORTED_RATE_THEORIES:
            raise ValueError(f"Unsupported solution rate theory: {theory}")
        path_claim = str(self.path_claim_id).strip() if self.path_claim_id is not None else None
        standard_state_correction = (
            str(self.standard_state_correction).strip()
            if self.standard_state_correction is not None else None
        )
        hashes = tuple(_text(item, name="calculation_hash") for item in self.calculation_hashes)
        if status == "verified_solution_tst":
            if not path_claim:
                raise ValueError("A verified solution-TST rate requires path_claim_id")
            if standard_state_correction != "one_molar_solution":
                raise ValueError(
                    "A verified solution-TST rate requires standard_state_correction='one_molar_solution'"
                )
            if not hashes:
                raise ValueError("A verified solution-TST rate requires calculation_hashes")
        lower_multiplier = self.lower_rate_multiplier
        upper_multiplier = self.upper_rate_multiplier
        uncertainty_source = (
            str(self.uncertainty_source).strip() if self.uncertainty_source is not None else None
        )
        if lower_multiplier is not None or upper_multiplier is not None or uncertainty_source is not None:
            if lower_multiplier is None or upper_multiplier is None or not uncertainty_source:
                raise ValueError(
                    "Rate uncertainty requires lower/upper multipliers and uncertainty_source"
                )
            lower_multiplier = _finite(lower_multiplier, name="lower_rate_multiplier")
            upper_multiplier = _finite(upper_multiplier, name="upper_rate_multiplier")
            if lower_multiplier < 0.0 or lower_multiplier > 1.0 or upper_multiplier < 1.0:
                raise ValueError("Rate uncertainty must satisfy 0 <= lower <= 1 <= upper")
        uncertainty_components = tuple(
            _text(item, name="uncertainty_component") for item in self.uncertainty_components
        )
        flags = {str(key): bool(value) for key, value in dict(self.applicability_flags).items()}
        unknown_flags = set(flags) - PROHIBITED_MECHANISM_FLAGS
        if unknown_flags:
            raise ValueError("Unknown solution mechanism applicability flags: " + ", ".join(sorted(unknown_flags)))
        object.__setattr__(self, "reaction_id", reaction_id)
        object.__setattr__(self, "reactant_id", reactant)
        object.__setattr__(self, "products", dict(sorted(products.items())))
        object.__setattr__(self, "rate_constant_s", rate)
        object.__setattr__(self, "temperature_K", _finite(self.temperature_K, name="temperature_K"))
        object.__setattr__(self, "solvent", _text(self.solvent, name="solvent"))
        object.__setattr__(self, "solvent_model", _text(self.solvent_model, name="solvent_model").upper())
        object.__setattr__(
            self, "standard_state_mol_L",
            _finite(self.standard_state_mol_L, name="standard_state_mol_L"),
        )
        object.__setattr__(self, "verification_status", status)
        object.__setattr__(self, "rate_theory", theory)
        object.__setattr__(self, "source", _text(self.source, name="source"))
        object.__setattr__(self, "rate_claim_id", _text(self.rate_claim_id, name="rate_claim_id"))
        object.__setattr__(self, "path_claim_id", path_claim)
        object.__setattr__(self, "standard_state_correction", standard_state_correction)
        object.__setattr__(self, "calculation_hashes", hashes)
        object.__setattr__(self, "lower_rate_multiplier", lower_multiplier)
        object.__setattr__(self, "upper_rate_multiplier", upper_multiplier)
        object.__setattr__(self, "uncertainty_source", uncertainty_source)
        object.__setattr__(self, "uncertainty_components", uncertainty_components)
        object.__setattr__(self, "applicability_flags", dict(sorted(flags.items())))

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "VerifiedSolutionReaction":
        molecularity = value.get("molecularity", 1)
        if isinstance(molecularity, bool) or int(molecularity) != 1:
            raise ValueError("The first solution lane accepts molecularity=1 only")
        units = str(value.get("rate_units", "s^-1")).strip().casefold().replace(" ", "")
        if units not in {"s^-1", "s-1", "1/s"}:
            raise ValueError("First-order solution rate units must be s^-1")
        if value.get("reactants") is not None:
            reactants = value.get("reactants")
            if not isinstance(reactants, Mapping) or len(reactants) != 1:
                raise ValueError("The first solution lane accepts exactly one reactant species")
            reactant_id, coefficient = next(iter(reactants.items()))
            if int(coefficient) != 1 or str(reactant_id) != str(value.get("reactant_id")):
                raise ValueError("reactants must agree with the unit-stoichiometry reactant_id")
        return cls(
            reaction_id=value.get("reaction_id"), reactant_id=value.get("reactant_id"),
            products=value.get("products") or {}, rate_constant_s=value.get("rate_constant_s"),
            temperature_K=value.get("temperature_K"), solvent=value.get("solvent"),
            solvent_model=value.get("solvent_model"),
            standard_state_mol_L=value.get("standard_state_mol_L", 1.0),
            verification_status=value.get("verification_status"),
            rate_theory=value.get("rate_theory"), source=value.get("source"),
            rate_claim_id=value.get("rate_claim_id"), path_claim_id=value.get("path_claim_id"),
            standard_state_correction=value.get("standard_state_correction"),
            calculation_hashes=tuple(value.get("calculation_hashes") or ()),
            lower_rate_multiplier=value.get("lower_rate_multiplier"),
            upper_rate_multiplier=value.get("upper_rate_multiplier"),
            uncertainty_source=value.get("uncertainty_source"),
            uncertainty_components=tuple(value.get("uncertainty_components") or ()),
            applicability_flags=value.get("applicability_flags") or {},
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "reaction_id": self.reaction_id, "reactant_id": self.reactant_id,
            "products": dict(self.products), "molecularity": 1,
            "rate_constant_s": self.rate_constant_s, "rate_units": "s^-1",
            "temperature_K": self.temperature_K, "solvent": self.solvent,
            "solvent_model": self.solvent_model,
            "standard_state_mol_L": self.standard_state_mol_L,
            "verification_status": self.verification_status, "rate_theory": self.rate_theory,
            "source": self.source, "rate_claim_id": self.rate_claim_id,
            "path_claim_id": self.path_claim_id,
            "standard_state_correction": self.standard_state_correction,
            "calculation_hashes": list(self.calculation_hashes),
            "lower_rate_multiplier": self.lower_rate_multiplier,
            "upper_rate_multiplier": self.upper_rate_multiplier,
            "uncertainty_source": self.uncertainty_source,
            "uncertainty_components": list(self.uncertainty_components),
            "applicability_flags": dict(self.applicability_flags),
        }


@dataclass(frozen=True)
class VerifiedSolutionNetwork:
    target_species_id: str
    species: tuple[SolutionSpecies, ...]
    reactions: tuple[VerifiedSolutionReaction, ...]
    coverage_status: str
    coverage_claim_id: str | None = None
    discovery_sources: tuple[str, ...] = ()
    unresolved_candidate_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        target = _text(self.target_species_id, name="target_species_id")
        if not self.species or not self.reactions:
            raise ValueError("A verified solution network needs species and reactions")
        species_ids = [item.species_id for item in self.species]
        reaction_ids = [item.reaction_id for item in self.reactions]
        if len(set(species_ids)) != len(species_ids):
            raise ValueError("Solution species IDs must be unique")
        smiles_identities = [item.smiles for item in self.species]
        if len(set(smiles_identities)) != len(smiles_identities):
            raise ValueError("Solution species must have unique canonical molecular identities")
        if len(set(reaction_ids)) != len(reaction_ids):
            raise ValueError("Solution reaction IDs must be unique")
        if target not in species_ids:
            raise ValueError("Target species is absent from the solution network")
        target_species = next(item for item in self.species if item.species_id == target)
        if target_species.initial_concentration_mol_L <= 0.0:
            raise ValueError("Target initial concentration must be positive")
        available = set(species_ids)
        for reaction in self.reactions:
            referenced = {reaction.reactant_id, *reaction.products}
            if not referenced <= available:
                raise ValueError(
                    f"Reaction {reaction.reaction_id} references unknown species: "
                    + ", ".join(sorted(referenced - available))
                )
        if any(target in reaction.products for reaction in self.reactions):
            raise ValueError(
                "The first solution lane requires monotonic target loss; target-producing reactions are unsupported"
            )
        coverage_status = _text(self.coverage_status, name="coverage_status")
        if coverage_status not in {
            "verified_complete_within_declared_unimolecular_scope",
            "bounded_search_incomplete",
            "unresolved",
        }:
            raise ValueError("Unsupported solution network coverage status")
        coverage_claim = (
            str(self.coverage_claim_id).strip() if self.coverage_claim_id is not None else None
        )
        sources = tuple(_text(item, name="discovery_source") for item in self.discovery_sources)
        unresolved = tuple(_text(item, name="unresolved_candidate_id") for item in self.unresolved_candidate_ids)
        if coverage_status == "verified_complete_within_declared_unimolecular_scope":
            if not coverage_claim or not sources or unresolved:
                raise ValueError(
                    "Verified restricted network coverage requires a claim ID, discovery sources, "
                    "and no unresolved candidates"
                )
        object.__setattr__(self, "target_species_id", target)
        object.__setattr__(self, "coverage_status", coverage_status)
        object.__setattr__(self, "coverage_claim_id", coverage_claim)
        object.__setattr__(self, "discovery_sources", sources)
        object.__setattr__(self, "unresolved_candidate_ids", unresolved)
        self._validate_balances()

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "VerifiedSolutionNetwork":
        if not isinstance(value, Mapping):
            raise TypeError("solution network must be a mapping")
        return cls(
            target_species_id=value.get("target_species_id"),
            species=tuple(SolutionSpecies.from_dict(item) for item in value.get("species") or ()),
            reactions=tuple(
                VerifiedSolutionReaction.from_dict(item) for item in value.get("reactions") or ()
            ),
            coverage_status=value.get("coverage_status", "unresolved"),
            coverage_claim_id=value.get("coverage_claim_id"),
            discovery_sources=tuple(value.get("discovery_sources") or ()),
            unresolved_candidate_ids=tuple(value.get("unresolved_candidate_ids") or ()),
        )

    def _validate_balances(self) -> None:
        by_id = {item.species_id: item for item in self.species}
        for reaction in self.reactions:
            reactant = by_id[reaction.reactant_id]
            elements = set(reactant.element_counts)
            for product_id in reaction.products:
                elements.update(by_id[product_id].element_counts)
            mismatches = {}
            for element in sorted(elements):
                left = reactant.element_counts.get(element, 0)
                right = sum(
                    coefficient * by_id[product_id].element_counts.get(element, 0)
                    for product_id, coefficient in reaction.products.items()
                )
                if left != right:
                    mismatches[element] = {"reactants": left, "products": right}
            product_charge = sum(
                coefficient * by_id[product_id].charge
                for product_id, coefficient in reaction.products.items()
            )
            if mismatches or reactant.charge != product_charge:
                raise ValueError(
                    f"Reaction {reaction.reaction_id} is not element/charge balanced: "
                    f"elements={mismatches}, charge={reactant.charge}->{product_charge}"
                )

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SOLUTION_NETWORK_SCHEMA_VERSION,
            "kind": "verified_implicit_solution_unimolecular_network",
            "target_species_id": self.target_species_id,
            "coverage_status": self.coverage_status,
            "coverage_claim_id": self.coverage_claim_id,
            "discovery_sources": list(self.discovery_sources),
            "unresolved_candidate_ids": list(self.unresolved_candidate_ids),
            "species": [item.as_dict() for item in self.species],
            "reactions": [item.as_dict() for item in self.reactions],
        }


def assess_solution_network_applicability(
    environment: ImplicitSolutionEnvironment,
    network: VerifiedSolutionNetwork,
) -> dict[str, Any]:
    blockers: list[dict[str, Any]] = []
    if network.coverage_status != "verified_complete_within_declared_unimolecular_scope":
        blockers.append({
            "reason": "reaction_network_coverage_unresolved",
            "coverage_status": network.coverage_status,
            "coverage_claim_id": network.coverage_claim_id,
            "unresolved_candidate_ids": list(network.unresolved_candidate_ids),
        })
    if not math.isclose(environment.standard_state_mol_L, 1.0, rel_tol=0.0, abs_tol=1e-12):
        blockers.append({"reason": "unsupported_standard_state", "required_mol_L": 1.0})
    ionic = [item.species_id for item in network.species if item.charge != 0]
    if ionic:
        blockers.append({"reason": "ionic_species_outside_neutral_lane", "species_ids": ionic})
    for reaction in network.reactions:
        reaction_blockers = []
        if (
            reaction.reactant_id == network.target_species_id
            and (
                reaction.lower_rate_multiplier is None
                or reaction.upper_rate_multiplier is None
                or not reaction.uncertainty_source
            )
        ):
            reaction_blockers.append("target_loss_rate_uncertainty_missing")
        if not math.isclose(reaction.temperature_K, environment.temperature_K, rel_tol=0.0, abs_tol=1e-8):
            reaction_blockers.append("rate_temperature_mismatch")
        if reaction.solvent.casefold() != environment.solvent.casefold():
            reaction_blockers.append("rate_solvent_mismatch")
        if reaction.solvent_model != environment.solvent_model:
            reaction_blockers.append("rate_solvent_model_mismatch")
        if not math.isclose(
            reaction.standard_state_mol_L, environment.standard_state_mol_L,
            rel_tol=0.0, abs_tol=1e-12,
        ):
            reaction_blockers.append("rate_standard_state_mismatch")
        active_flags = sorted(key for key, active in reaction.applicability_flags.items() if active)
        if active_flags:
            reaction_blockers.extend(active_flags)
        if reaction_blockers:
            blockers.append({"reaction_id": reaction.reaction_id, "reasons": reaction_blockers})
    return {
        "schema_version": 1,
        "kind": "implicit_solution_lane_applicability_assessment",
        "status": "supported" if not blockers else "outside_method_applicability",
        "blocking_findings": blockers,
        "scope": (
            "Dark neutral homogeneous implicit solution; verified elementary first-order rates; "
            "constant volume; monotonic target loss."
        ),
    }


_PADE13_COEFFICIENTS = (
    64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
    1187353796428800.0, 129060195264000.0, 10559470521600.0,
    670442572800.0, 33522128640.0, 1323241920.0, 40840800.0,
    960960.0, 16380.0, 182.0, 1.0,
)


def _matrix_exponential(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    """Scaling-and-squaring Padé(13), retaining the linear-solve condition."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not np.all(np.isfinite(matrix)):
        raise ValueError("Matrix exponential requires a finite square matrix")
    size = matrix.shape[0]
    identity = np.eye(size)
    norm = float(np.linalg.norm(matrix, 1))
    if norm == 0.0:
        return identity, 1.0
    scale = max(0, int(math.ceil(math.log2(norm / 5.371920351148152))))
    reduced = matrix / (2.0 ** scale)
    a2 = reduced @ reduced
    a4 = a2 @ a2
    a6 = a4 @ a2
    b = _PADE13_COEFFICIENTS
    u = reduced @ (
        a6 @ (b[13] * a6 + b[11] * a4 + b[9] * a2)
        + b[7] * a6 + b[5] * a4 + b[3] * a2 + b[1] * identity
    )
    v = (
        a6 @ (b[12] * a6 + b[10] * a4 + b[8] * a2)
        + b[6] * a6 + b[4] * a4 + b[2] * a2 + b[0] * identity
    )
    denominator = v - u
    condition = float(np.linalg.cond(denominator))
    if not math.isfinite(condition) or condition > 1.0e14:
        raise ValueError("Matrix exponential linear solve is numerically ill-conditioned")
    result = np.linalg.solve(denominator, v + u)
    for _ in range(scale):
        result = result @ result
    if not np.all(np.isfinite(result)):
        raise ValueError("Matrix exponential produced non-finite values")
    return result, condition


def _kinetic_matrix(network: VerifiedSolutionNetwork) -> tuple[np.ndarray, np.ndarray]:
    species_index = {item.species_id: index for index, item in enumerate(network.species)}
    species_count = len(network.species)
    reaction_count = len(network.reactions)
    kinetic = np.zeros((species_count, species_count), dtype=float)
    extent = np.zeros((reaction_count, species_count), dtype=float)
    for reaction_index, reaction in enumerate(network.reactions):
        reactant_index = species_index[reaction.reactant_id]
        kinetic[reactant_index, reactant_index] -= reaction.rate_constant_s
        extent[reaction_index, reactant_index] = reaction.rate_constant_s
        for product_id, coefficient in reaction.products.items():
            kinetic[species_index[product_id], reactant_index] += coefficient * reaction.rate_constant_s
    return kinetic, extent


def _element_inventory(network: VerifiedSolutionNetwork, concentrations: np.ndarray) -> dict[str, float]:
    elements = sorted({element for item in network.species for element in item.element_counts})
    return {
        element: float(sum(
            concentrations[index] * item.element_counts.get(element, 0)
            for index, item in enumerate(network.species)
        ))
        for element in elements
    }


def propagate_verified_solution_network(
    environment: ImplicitSolutionEnvironment,
    network: VerifiedSolutionNetwork,
    *,
    profile_points: int = 160,
) -> dict[str, Any]:
    applicability = assess_solution_network_applicability(environment, network)
    if applicability["status"] != "supported":
        return {
            "schema_version": 1,
            "kind": "verified_solution_concentration_propagation",
            "status": "outside_method_applicability",
            "applicability": applicability,
            "estimated_time_to_retention_seconds": None,
        }
    if profile_points < 20 or profile_points > 1000:
        raise ValueError("profile_points must be in [20, 1000]")
    kinetic, extent_matrix = _kinetic_matrix(network)
    species_count = len(network.species)
    reaction_count = len(network.reactions)
    augmented = np.zeros((species_count + reaction_count, species_count + reaction_count))
    augmented[:species_count, :species_count] = kinetic
    augmented[species_count:, :species_count] = extent_matrix
    initial_species = np.asarray(
        [item.initial_concentration_mol_L for item in network.species], dtype=float,
    )
    initial = np.concatenate((initial_species, np.zeros(reaction_count)))
    target_index = next(
        index for index, item in enumerate(network.species)
        if item.species_id == network.target_species_id
    )
    initial_target = float(initial_species[target_index])
    horizon = environment.target_duration_seconds
    max_rate = max(reaction.rate_constant_s for reaction in network.reactions)
    earliest = min(horizon, max(horizon * 1e-12, min(1.0, 1.0e-6 / max_rate)))
    if earliest >= horizon:
        times = np.asarray([0.0, horizon])
    else:
        times = np.concatenate(([0.0], np.geomspace(earliest, horizon, profile_points - 1)))
    profile = []
    maximum_condition = 1.0
    most_negative = 0.0

    def state_at(time_seconds: float) -> np.ndarray:
        nonlocal maximum_condition, most_negative
        exponential, condition = _matrix_exponential(augmented * float(time_seconds))
        maximum_condition = max(maximum_condition, condition)
        state = exponential @ initial
        most_negative = min(most_negative, float(np.min(state)))
        tolerance = max(1.0, float(np.max(initial_species))) * 1.0e-10
        if np.min(state) < -tolerance:
            raise ValueError("Concentration propagation produced a material negative value")
        state[np.abs(state) < tolerance] = 0.0
        return state

    final_state = None
    for time_seconds in times:
        state = state_at(float(time_seconds))
        final_state = state
        concentrations = state[:species_count]
        profile.append({
            "time_seconds": float(time_seconds),
            "target_fraction_remaining": float(concentrations[target_index] / initial_target),
            "concentrations_mol_L": {
                item.species_id: float(concentrations[index])
                for index, item in enumerate(network.species)
            },
        })
    assert final_state is not None
    threshold = environment.retention_fraction * initial_target
    final_target = float(final_state[target_index])
    t_retention = None
    if final_target <= threshold:
        lower, upper = 0.0, horizon
        for _ in range(80):
            midpoint = 0.5 * (lower + upper)
            if state_at(midpoint)[target_index] <= threshold:
                upper = midpoint
            else:
                lower = midpoint
        t_retention = upper

    target_loss_reactions = [
        reaction for reaction in network.reactions
        if reaction.reactant_id == network.target_species_id
    ]
    central_target_rate = sum(reaction.rate_constant_s for reaction in target_loss_reactions)
    lower_target_rate = sum(
        reaction.rate_constant_s * float(reaction.lower_rate_multiplier)
        for reaction in target_loss_reactions
    )
    upper_target_rate = sum(
        reaction.rate_constant_s * float(reaction.upper_rate_multiplier)
        for reaction in target_loss_reactions
    )
    crossing_log = -math.log(environment.retention_fraction)
    t95_fast = crossing_log / upper_target_rate if upper_target_rate > 0.0 else None
    t95_slow = crossing_log / lower_target_rate if lower_target_rate > 0.0 else None
    uncertainty_crossing_status = (
        "crossing_within_horizon_for_all_bounds"
        if t95_slow is not None and t95_slow <= horizon
        else "no_crossing_within_horizon_for_any_bound"
        if t95_fast is None or t95_fast > horizon
        else "crossing_depends_on_rate_uncertainty"
    )

    final_concentrations = final_state[:species_count]
    extents = final_state[species_count:]
    initial_inventory = _element_inventory(network, initial_species)
    final_inventory = _element_inventory(network, final_concentrations)
    inventory_errors = {
        element: abs(final_inventory[element] - initial_inventory[element])
        / max(abs(initial_inventory[element]), 1.0e-300)
        for element in initial_inventory
    }
    maximum_inventory_error = max(inventory_errors.values(), default=0.0)
    if maximum_inventory_error > 1.0e-8:
        raise ValueError("Solution propagation failed the elemental inventory closure gate")
    target_loss = max(0.0, initial_target - final_target)
    target_loss_extents = [
        float(extents[index]) if reaction.reactant_id == network.target_species_id else 0.0
        for index, reaction in enumerate(network.reactions)
    ]
    attributed_target_loss = sum(target_loss_extents)
    target_flux_closure = abs(attributed_target_loss - target_loss) / max(target_loss, 1.0e-300)
    if target_loss > initial_target * 1.0e-12 and target_flux_closure > 1.0e-8:
        raise ValueError("Solution target-loss flux attribution failed numerical closure")
    reaction_flux = []
    for index, reaction in enumerate(network.reactions):
        amount = target_loss_extents[index]
        reaction_flux.append({
            "reaction_id": reaction.reaction_id,
            "reactant_id": reaction.reactant_id,
            "products": dict(reaction.products),
            "integrated_extent_mol_L": float(extents[index]),
            "integrated_target_loss_mol_L": amount,
            "target_loss_fraction": amount / target_loss if target_loss > 0.0 else 0.0,
            "rate_claim_id": reaction.rate_claim_id,
            "path_claim_id": reaction.path_claim_id,
        })
    dominant = max(reaction_flux, key=lambda item: item["integrated_target_loss_mol_L"])
    return {
        "schema_version": 1,
        "kind": "verified_solution_concentration_propagation",
        "status": "completed",
        "applicability": applicability,
        "environment_contract": environment.as_environment_spec().as_dict(),
        "network": network.as_dict(),
        "estimated_time_to_retention_seconds": t_retention,
        "t95_uncertainty": {
            "status": uncertainty_crossing_status,
            "lower_seconds": t95_fast,
            "upper_seconds": t95_slow,
            "central_total_target_loss_rate_s": central_target_rate,
            "lower_total_target_loss_rate_s": lower_target_rate,
            "upper_total_target_loss_rate_s": upper_target_rate,
            "bound_sources": sorted({
                reaction.uncertainty_source for reaction in target_loss_reactions
                if reaction.uncertainty_source
            }),
            "interpretation": (
                "Combined explicit multiplicative bounds for all target-consuming first-order rates; "
                "not a statistical confidence interval."
            ),
        },
        "retention_fraction": environment.retention_fraction,
        "target_duration_seconds": horizon,
        "target_species_id": network.target_species_id,
        "target_fraction_remaining": final_target / initial_target,
        "target_loss_fraction": target_loss / initial_target,
        "final_concentrations_mol_L": {
            item.species_id: float(final_concentrations[index])
            for index, item in enumerate(network.species)
        },
        "reaction_flux_attribution": {
            "status": "completed",
            "reactions": reaction_flux,
            "dominant_target_loss_reaction_id": dominant["reaction_id"] if target_loss > 0.0 else None,
            "numerical_closure_relative_error": target_flux_closure,
        },
        "numerical_validation": {
            "status": "passed",
            "method": "linear_first_order_augmented_matrix_exponential_pade13",
            "maximum_element_inventory_relative_error": maximum_inventory_error,
            "element_inventory_relative_errors": inventory_errors,
            "maximum_pade_linear_solve_condition_number": maximum_condition,
            "most_negative_raw_state_value": most_negative,
            "nonnegative_concentrations": most_negative >= -max(1.0, float(np.max(initial_species))) * 1.0e-10,
        },
        "target_profile": profile,
        "claim_scope": (
            "Conditional on the supplied verified first-order solution rates and complete network within "
            "the lane's restricted applicability domain."
        ),
    }


__all__ = [
    "ImplicitSolutionEnvironment", "PROHIBITED_MECHANISM_FLAGS", "SolutionSpecies",
    "VerifiedSolutionNetwork", "VerifiedSolutionReaction",
    "assess_solution_network_applicability", "propagate_verified_solution_network",
]
