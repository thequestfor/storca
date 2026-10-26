"""Speculative reaction-network contracts for solution verification allocation.

Nothing in this module is automatically admissible kinetic evidence.  It
exists beside, not inside, :mod:`storca.solution_kinetics` so screening priors
cannot silently weaken the verified solution-network contract.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from typing import Any, Mapping

from .solution_kinetics import (
    ImplicitSolutionEnvironment,
    PROHIBITED_MECHANISM_FLAGS,
    SolutionSpecies,
    VerifiedSolutionReaction,
)


RATE_INTERVAL_STATUSES = {
    "verified", "estimated_bounded", "screening_prior_only", "unbounded", "unavailable",
}
RATE_INTERVAL_USE_SCOPES = {"prediction_bound", "allocation_only", "none"}
VERIFICATION_STATES = {
    "pending", "verified", "path_unresolved", "method_escalation_required",
    "ruled_out_with_evidence",
}


def _finite(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _text(value: Any, name: str) -> str:
    result = str(value or "").strip()
    if not result:
        raise ValueError(f"{name} is required")
    return result


@dataclass(frozen=True)
class ApproximateRateInterval:
    """Condition-specific rate information with an explicit permitted use."""

    status: str
    lower_rate_s: float | None
    upper_rate_s: float | None
    central_rate_s: float | None = None
    source: str = ""
    source_kind: str = "unknown"
    use_scope: str = "allocation_only"
    temperature_K: float = 298.15
    solvent: str = "unknown"
    solvent_model: str = "unknown"
    uncertainty_components: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.status not in RATE_INTERVAL_STATUSES:
            raise ValueError(f"Unsupported approximate-rate status: {self.status}")
        if self.use_scope not in RATE_INTERVAL_USE_SCOPES:
            raise ValueError(f"Unsupported approximate-rate use scope: {self.use_scope}")
        lower = None if self.lower_rate_s is None else _finite(self.lower_rate_s, "lower_rate_s")
        upper = None if self.upper_rate_s is None else _finite(self.upper_rate_s, "upper_rate_s")
        central = None if self.central_rate_s is None else _finite(self.central_rate_s, "central_rate_s")
        if lower is not None and lower < 0.0:
            raise ValueError("Rate lower bound must be nonnegative")
        if upper is not None and upper <= 0.0:
            raise ValueError("Rate upper bound must be positive")
        if lower is not None and upper is not None and lower > upper:
            raise ValueError("Rate bounds must satisfy lower <= upper")
        if central is not None and central <= 0.0:
            raise ValueError("Central rate must be positive")
        if central is not None and lower is not None and central < lower:
            raise ValueError("Central rate lies below its lower bound")
        if central is not None and upper is not None and central > upper:
            raise ValueError("Central rate lies above its upper bound")
        if self.status in {"verified", "estimated_bounded"} and (lower is None or upper is None):
            raise ValueError(f"{self.status} rate evidence requires finite lower and upper bounds")
        if self.status == "verified" and self.use_scope != "prediction_bound":
            raise ValueError("Verified rate evidence must be usable as a prediction bound")
        if self.use_scope == "prediction_bound" and (lower is None or upper is None):
            raise ValueError("Prediction bounds require finite lower and upper rates")
        object.__setattr__(self, "lower_rate_s", lower)
        object.__setattr__(self, "upper_rate_s", upper)
        object.__setattr__(self, "central_rate_s", central)
        object.__setattr__(self, "source", _text(self.source, "rate interval source"))
        object.__setattr__(self, "source_kind", _text(self.source_kind, "source_kind"))
        object.__setattr__(self, "temperature_K", _finite(self.temperature_K, "temperature_K"))
        object.__setattr__(self, "solvent", _text(self.solvent, "solvent"))
        object.__setattr__(self, "solvent_model", _text(self.solvent_model, "solvent_model").upper())
        object.__setattr__(self, "uncertainty_components", tuple(str(item) for item in self.uncertainty_components))
        object.__setattr__(self, "limitations", tuple(str(item) for item in self.limitations))

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ApproximateRateInterval":
        if not isinstance(value, Mapping):
            raise TypeError("approximate_rate must be a mapping")
        return cls(
            status=str(value.get("status") or "unavailable"),
            lower_rate_s=value.get("lower_rate_s"), upper_rate_s=value.get("upper_rate_s"),
            central_rate_s=value.get("central_rate_s"), source=value.get("source") or "unspecified",
            source_kind=value.get("source_kind") or "unknown",
            use_scope=value.get("use_scope") or "allocation_only",
            temperature_K=value.get("temperature_K", 298.15),
            solvent=value.get("solvent") or "unknown",
            solvent_model=value.get("solvent_model") or "unknown",
            uncertainty_components=tuple(value.get("uncertainty_components") or ()),
            limitations=tuple(value.get("limitations") or ()),
        )

    @classmethod
    def from_verified_reaction(cls, reaction: VerifiedSolutionReaction) -> "ApproximateRateInterval":
        lower = reaction.rate_constant_s * float(reaction.lower_rate_multiplier or 1.0)
        upper = reaction.rate_constant_s * float(reaction.upper_rate_multiplier or 1.0)
        return cls(
            status="verified", lower_rate_s=lower, upper_rate_s=upper,
            central_rate_s=reaction.rate_constant_s, source=reaction.rate_claim_id,
            source_kind="verified_solution_rate_claim", use_scope="prediction_bound",
            temperature_K=reaction.temperature_K, solvent=reaction.solvent,
            solvent_model=reaction.solvent_model,
            uncertainty_components=reaction.uncertainty_components,
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "status": self.status, "lower_rate_s": self.lower_rate_s,
            "upper_rate_s": self.upper_rate_s, "central_rate_s": self.central_rate_s,
            "source": self.source, "source_kind": self.source_kind,
            "use_scope": self.use_scope, "temperature_K": self.temperature_K,
            "solvent": self.solvent, "solvent_model": self.solvent_model,
            "uncertainty_components": list(self.uncertainty_components),
            "limitations": list(self.limitations),
        }


@dataclass(frozen=True)
class SolutionCandidateReaction:
    candidate_id: str
    reactant_id: str
    products: Mapping[str, int]
    route: Mapping[str, Any]
    approximate_rate: ApproximateRateInterval
    discovery_sources: tuple[str, ...]
    discovery_evidence: tuple[Mapping[str, Any], ...] = ()
    uncertainty: Mapping[str, Any] = field(default_factory=dict)
    applicability_flags: Mapping[str, bool] = field(default_factory=dict)
    verification_state: str = "pending"
    estimated_cost_core_hours: float = 1.0
    verification_attempts: int = 0
    retained_result: str | None = None

    def __post_init__(self) -> None:
        candidate_id = _text(self.candidate_id, "candidate_id")
        reactant_id = _text(self.reactant_id, "reactant_id")
        products = {str(key): int(value) for key, value in dict(self.products).items()}
        if not products or any(not key or value <= 0 for key, value in products.items()):
            raise ValueError("Solution candidate needs positive integer product stoichiometry")
        if not isinstance(self.route, Mapping) or not self.route:
            raise ValueError("Solution candidate needs an executable or auditable route record")
        sources = tuple(sorted({_text(item, "discovery_source") for item in self.discovery_sources}))
        if not sources:
            raise ValueError("Solution candidate needs at least one discovery source")
        flags = {str(key): bool(value) for key, value in dict(self.applicability_flags).items()}
        unknown = set(flags) - PROHIBITED_MECHANISM_FLAGS
        if unknown:
            raise ValueError("Unknown solution applicability flags: " + ", ".join(sorted(unknown)))
        state = str(self.verification_state)
        if state not in VERIFICATION_STATES:
            raise ValueError(f"Unsupported solution verification state: {state}")
        cost = _finite(self.estimated_cost_core_hours, "estimated_cost_core_hours")
        if cost <= 0.0:
            raise ValueError("Estimated verification cost must be positive")
        if int(self.verification_attempts) < 0:
            raise ValueError("verification_attempts must be nonnegative")
        object.__setattr__(self, "candidate_id", candidate_id)
        object.__setattr__(self, "reactant_id", reactant_id)
        object.__setattr__(self, "products", dict(sorted(products.items())))
        object.__setattr__(self, "route", dict(self.route))
        object.__setattr__(self, "discovery_sources", sources)
        object.__setattr__(self, "discovery_evidence", tuple(dict(item) for item in self.discovery_evidence))
        object.__setattr__(self, "uncertainty", dict(self.uncertainty))
        object.__setattr__(self, "applicability_flags", dict(sorted(flags.items())))
        object.__setattr__(self, "verification_state", state)
        object.__setattr__(self, "estimated_cost_core_hours", cost)
        object.__setattr__(self, "verification_attempts", int(self.verification_attempts))

    @property
    def active_applicability_flags(self) -> tuple[str, ...]:
        return tuple(sorted(key for key, active in self.applicability_flags.items() if active))

    def with_verification(self, state: str, result: str | None) -> "SolutionCandidateReaction":
        return replace(
            self, verification_state=state, retained_result=result,
            verification_attempts=self.verification_attempts + 1,
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SolutionCandidateReaction":
        if not isinstance(value, Mapping):
            raise TypeError("candidate reaction must be a mapping")
        return cls(
            candidate_id=value.get("candidate_id") or value.get("reaction_id"),
            reactant_id=value.get("reactant_id"), products=value.get("products") or {},
            route=value.get("route") or {},
            approximate_rate=ApproximateRateInterval.from_dict(value.get("approximate_rate") or {}),
            discovery_sources=tuple(value.get("discovery_sources") or ()),
            discovery_evidence=tuple(value.get("discovery_evidence") or ()),
            uncertainty=value.get("uncertainty") or {},
            applicability_flags=value.get("applicability_flags") or {},
            verification_state=value.get("verification_state", "pending"),
            estimated_cost_core_hours=value.get("estimated_cost_core_hours", 1.0),
            verification_attempts=value.get("verification_attempts", 0),
            retained_result=value.get("retained_result"),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": 1, "kind": "storca_solution_candidate_reaction",
            "candidate_id": self.candidate_id, "reactant_id": self.reactant_id,
            "products": dict(self.products), "route": dict(self.route),
            "approximate_rate": self.approximate_rate.as_dict(),
            "discovery_sources": list(self.discovery_sources),
            "discovery_evidence": [dict(item) for item in self.discovery_evidence],
            "uncertainty": dict(self.uncertainty),
            "applicability_flags": dict(self.applicability_flags),
            "active_applicability_flags": list(self.active_applicability_flags),
            "verification_state": self.verification_state,
            "estimated_cost_core_hours": self.estimated_cost_core_hours,
            "verification_attempts": self.verification_attempts,
            "retained_result": self.retained_result,
        }


@dataclass(frozen=True)
class SolutionCandidateNetwork:
    environment: ImplicitSolutionEnvironment
    target_smiles: str
    target_species_id: str
    species: tuple[SolutionSpecies, ...]
    verified_reactions: tuple[VerifiedSolutionReaction, ...]
    candidates: tuple[SolutionCandidateReaction, ...]
    coverage_status: str
    coverage_claim_id: str | None
    discovery_sources: tuple[str, ...]
    unresolved_discovery_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        species_ids = [item.species_id for item in self.species]
        if len(species_ids) != len(set(species_ids)) or self.target_species_id not in species_ids:
            raise ValueError("Candidate network species IDs must be unique and contain the target")
        if any(item.charge != 0 for item in self.species):
            raise ValueError("The first solution verification session accepts neutral species only")
        candidate_ids = [item.candidate_id for item in self.candidates]
        verified_ids = [item.reaction_id for item in self.verified_reactions]
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValueError("Candidate IDs must be unique")
        if set(candidate_ids) & set(verified_ids):
            raise ValueError("A reaction cannot be both speculative and verified")
        by_id = {item.species_id: item for item in self.species}
        from .reaction_candidates import canonical_smiles
        if canonical_smiles(self.target_smiles) != by_id[self.target_species_id].smiles:
            raise ValueError("Candidate network target_smiles does not match its target species")
        if self.coverage_status not in {
            "verified_complete_within_declared_unimolecular_scope",
            "bounded_search_incomplete", "unresolved",
        }:
            raise ValueError("Unsupported candidate-network coverage status")
        for candidate in self.candidates:
            if candidate.reactant_id not in by_id or not set(candidate.products) <= set(by_id):
                raise ValueError(f"Candidate {candidate.candidate_id} references unknown species")
            self._validate_balance(candidate.reactant_id, candidate.products, by_id, candidate.candidate_id)
        for reaction in self.verified_reactions:
            if reaction.reactant_id not in by_id or not set(reaction.products) <= set(by_id):
                raise ValueError(f"Verified reaction {reaction.reaction_id} references unknown species")
            self._validate_balance(reaction.reactant_id, reaction.products, by_id, reaction.reaction_id)
            if (
                not math.isclose(reaction.temperature_K, self.environment.temperature_K, rel_tol=0.0, abs_tol=1e-8)
                or reaction.solvent.casefold() != self.environment.solvent.casefold()
                or reaction.solvent_model != self.environment.solvent_model
            ):
                raise ValueError(f"Verified reaction {reaction.reaction_id} does not match the session environment")
        if any(
            self.target_species_id in reaction.products for reaction in self.verified_reactions
        ) or any(self.target_species_id in candidate.products for candidate in self.candidates):
            raise ValueError("The t95 session requires monotonic target loss without target production")
        object.__setattr__(self, "discovery_sources", tuple(sorted(set(self.discovery_sources))))
        object.__setattr__(self, "unresolved_discovery_ids", tuple(sorted(set(self.unresolved_discovery_ids))))

    @staticmethod
    def _validate_balance(reactant_id, products, by_id, reaction_id) -> None:
        reactant = by_id[reactant_id]
        elements = set(reactant.element_counts)
        for product_id in products:
            elements.update(by_id[product_id].element_counts)
        mismatch = {
            element: {
                "reactants": reactant.element_counts.get(element, 0),
                "products": sum(
                    int(coefficient) * by_id[product_id].element_counts.get(element, 0)
                    for product_id, coefficient in products.items()
                ),
            }
            for element in elements
            if reactant.element_counts.get(element, 0) != sum(
                int(coefficient) * by_id[product_id].element_counts.get(element, 0)
                for product_id, coefficient in products.items()
            )
        }
        product_charge = sum(int(coefficient) * by_id[product_id].charge for product_id, coefficient in products.items())
        if mismatch or product_charge != reactant.charge:
            raise ValueError(f"Reaction {reaction_id} is not element/charge balanced")

    def replace_candidate_with_verified(
        self, candidate_id: str, reaction: VerifiedSolutionReaction,
    ) -> "SolutionCandidateNetwork":
        matches = [item for item in self.candidates if item.candidate_id == candidate_id]
        if len(matches) != 1 or reaction.reaction_id != candidate_id:
            raise ValueError("Verified reaction does not match exactly one candidate")
        return replace(
            self,
            candidates=tuple(item for item in self.candidates if item.candidate_id != candidate_id),
            verified_reactions=(*self.verified_reactions, reaction),
        )

    def replace_candidate_state(
        self, candidate_id: str, state: str, result: str | None,
    ) -> "SolutionCandidateNetwork":
        if sum(item.candidate_id == candidate_id for item in self.candidates) != 1:
            raise ValueError("Candidate state update must match exactly one candidate")
        return replace(
            self,
            candidates=tuple(
                item.with_verification(state, result) if item.candidate_id == candidate_id else item
                for item in self.candidates
            ),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": 1, "kind": "storca_solution_candidate_network",
            "target_smiles": self.target_smiles,
            "environment": {
                "solvent": self.environment.solvent,
                "solvent_model": self.environment.solvent_model,
                "temperature_K": self.environment.temperature_K,
                "target_duration_seconds": self.environment.target_duration_seconds,
                "retention_fraction": self.environment.retention_fraction,
                "standard_state_mol_L": self.environment.standard_state_mol_L,
                "volume_L": self.environment.volume_L,
                "light_condition": self.environment.light_condition,
                "electronic_structure_profile": self.environment.electronic_structure_profile,
            },
            "target_species_id": self.target_species_id,
            "coverage_status": self.coverage_status,
            "coverage_claim_id": self.coverage_claim_id,
            "discovery_sources": list(self.discovery_sources),
            "unresolved_discovery_ids": list(self.unresolved_discovery_ids),
            "species": [item.as_dict() for item in self.species],
            "verified_reactions": [item.as_dict() for item in self.verified_reactions],
            "candidates": [item.as_dict() for item in self.candidates],
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SolutionCandidateNetwork":
        return cls(
            environment=ImplicitSolutionEnvironment.from_dict(value.get("environment") or {}),
            target_smiles=value.get("target_smiles"),
            target_species_id=value.get("target_species_id"),
            species=tuple(SolutionSpecies.from_dict(item) for item in value.get("species") or ()),
            verified_reactions=tuple(
                VerifiedSolutionReaction.from_dict(item) for item in value.get("verified_reactions") or ()
            ),
            candidates=tuple(
                SolutionCandidateReaction.from_dict(item) for item in value.get("candidates") or ()
            ),
            coverage_status=value.get("coverage_status", "unresolved"),
            coverage_claim_id=value.get("coverage_claim_id"),
            discovery_sources=tuple(value.get("discovery_sources") or ()),
            unresolved_discovery_ids=tuple(value.get("unresolved_discovery_ids") or ()),
        )


__all__ = [
    "ApproximateRateInterval", "RATE_INTERVAL_STATUSES", "RATE_INTERVAL_USE_SCOPES",
    "SolutionCandidateNetwork", "SolutionCandidateReaction", "VERIFICATION_STATES",
]
