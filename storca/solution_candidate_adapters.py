"""Adapters from retained discovery records to solution-session candidates."""

from __future__ import annotations

from collections import Counter
from dataclasses import fields, replace
import math
from typing import Any, Mapping

from .reaction_candidates import ReactionCandidate, candidate_from_rmg_route, canonical_smiles
from .route_verify import RouteSpec
from .solution_candidates import (
    ApproximateRateInterval, SolutionCandidateNetwork, SolutionCandidateReaction,
)
from .solution_kinetics import ImplicitSolutionEnvironment, SolutionSpecies, VerifiedSolutionReaction
from .solution_path_verification import REQUIRED_UNCERTAINTY_COMPONENTS


BOLTZMANN_J_K = 1.380649e-23
PLANCK_J_S = 6.62607015e-34


def _clean_route(value: Mapping[str, Any], *, parent_smiles: str) -> dict[str, Any]:
    allowed = {item.name for item in fields(RouteSpec)}
    route = {key: item for key, item in value.items() if key in allowed}
    route.setdefault("parent_smiles", parent_smiles)
    # RouteSpec validates required fields and catches undeclared/zero spin at
    # the path layer; this construction establishes a stable serializable shape.
    spec = RouteSpec(**route)
    return spec.as_dict()


def _species_identity_index(species: tuple[SolutionSpecies, ...]) -> dict[str, str]:
    grouped: dict[str, list[str]] = {}
    for item in species:
        grouped.setdefault(canonical_smiles(item.smiles), []).append(item.species_id)
    ambiguous = {key: ids for key, ids in grouped.items() if len(ids) != 1}
    if ambiguous:
        raise ValueError("Session species contain duplicate canonical molecular identities")
    return {key: ids[0] for key, ids in grouped.items()}


def run_solution_candidate_rate_coefficients(
    value: Mapping[str, Any], environment: Mapping[str, Any], workdir,
    **kwargs: Any,
) -> dict[str, Any]:
    """Native solvent-bound discovery → existing path/conditional coefficients.

    No preassembled species table, concentration history, supplied kinetics or
    barrier-error bound is required. This is an internal candidate adapter, not
    an autonomous discovery engine or admission to verified mass-action rates.
    Gas/unspecified screening evidence is not promoted to solvent-bound evidence.
    """
    native = ReactionCandidate.from_dict(dict(value))
    if not native.discovery_sources:
        raise ValueError("Solution coefficient candidate needs retained discovery evidence")
    if (native.environment.get("solvent_model") != environment.get("solvent_model")
            or native.environment.get("solvent") != environment.get("solvent")):
        raise ValueError("Native coefficient candidate is not bound to the requested solution environment")
    source_validation = None
    if "graph_xtb" in native.discovery_sources:
        from .graph_xtb_discovery import validate_graph_xtb_candidate_sources
        source_validation = validate_graph_xtb_candidate_sources(native)
    route = native.as_route_dict(parent_smiles=native.reactants[0])
    # Identical stoichiometric copies are one chemical species, even when
    # the native route uses per-fragment display labels for atom mapping.
    by_identity = {}
    for side, smiles_values in (("reactants", native.reactants), ("products", native.products)):
        labels = route["reactant_labels" if side == "reactants" else "product_labels"]
        normalized = [by_identity.setdefault(canonical_smiles(smiles), label)
                      for smiles, label in zip(smiles_values, labels)]
        route["reactant_labels" if side == "reactants" else "product_labels"] = normalized
    request = {
        "schema_version": 1, "kind": "implicit_solution_rate_coefficient_request",
        "environment": dict(environment), "route": _clean_route(route, parent_smiles=native.reactants[0]),
        "discovery_candidate": native.as_dict(),
        "applicability_flags": dict(value.get("applicability_flags") or {}),
    }
    if source_validation is not None:
        request["discovery_source_validation"] = source_validation
    from .solution_path_verification import run_implicit_solution_path_verification
    result = run_implicit_solution_path_verification(request, workdir, **kwargs)
    if source_validation is not None and validate_graph_xtb_candidate_sources(native) != source_validation:
        raise ValueError("Native discovery sources changed during solution path production")
    return result


def _native_candidate(
    value: Mapping[str, Any], environment: ImplicitSolutionEnvironment,
    species: tuple[SolutionSpecies, ...], target_smiles: str,
) -> SolutionCandidateReaction:
    native = ReactionCandidate.from_dict(dict(value))
    if len(native.reactants) != 1:
        raise ValueError("The first solution verification session accepts unimolecular candidates only")
    if native.total_charge not in {None, 0}:
        raise ValueError("The first solution verification session accepts neutral candidates only")
    identity = _species_identity_index(species)
    reactant_key = canonical_smiles(native.reactants[0])
    if reactant_key not in identity:
        raise ValueError("Native candidate reactant is not present in the solution species table")
    product_counts = Counter(canonical_smiles(item) for item in native.products)
    if not set(product_counts) <= set(identity):
        raise ValueError("Native candidate products are not all present in the solution species table")
    products = {identity[key]: count for key, count in product_counts.items()}
    source_blind = not native.environment or not native.environment.get("solvent_model")
    if not source_blind:
        model = str(native.environment.get("solvent_model") or "").upper()
        solvent = str(native.environment.get("solvent") or "")
        if model != environment.solvent_model or solvent.casefold() != environment.solvent.casefold():
            raise ValueError("Discovery candidate declares a different solution environment")
    prefactor = BOLTZMANN_J_K * environment.temperature_K / PLANCK_J_S
    prior = dict(value.get("solution_rate_prior") or {})
    central = prior.get("central_rate_s")
    if central is not None:
        central = float(central)
        if not math.isfinite(central) or central <= 0.0 or central > prefactor:
            raise ValueError("Solution scheduling prior must be positive and below kBT/h")
    interval = ApproximateRateInterval(
        status="screening_prior_only" if central is not None else "unbounded",
        lower_rate_s=0.0, upper_rate_s=prefactor, central_rate_s=central,
        source=str(prior.get("source") or "+".join(native.discovery_sources) or "native discovery candidate"),
        source_kind=("solution_screening_prior" if central is not None else "nonnegative_barrier_tst_ceiling"),
        use_scope="allocation_only", temperature_K=environment.temperature_K,
        solvent=environment.solvent, solvent_model=environment.solvent_model,
        uncertainty_components=("mechanism_candidate", "rate_unresolved"),
        limitations=(
            "Discovery evidence is not a verified solution rate.",
            "The kBT/h ceiling is used for allocation only and assumes a nonnegative conventional barrier.",
        ),
    )
    route = _clean_route(
        native.as_route_dict(parent_smiles=native.reactants[0], target_label=identity[reactant_key]),
        parent_smiles=native.reactants[0],
    )
    flags = dict(value.get("applicability_flags") or {})
    return SolutionCandidateReaction(
        candidate_id=native.candidate_id, reactant_id=identity[reactant_key], products=products,
        route=route, approximate_rate=interval,
        discovery_sources=native.discovery_sources or ("unknown_native_discovery",),
        discovery_evidence=native.discovery_assertions,
        uncertainty=value.get("solution_verification_uncertainty") or {},
        applicability_flags=flags,
        estimated_cost_core_hours=value.get("estimated_cost_core_hours", 100.0),
    )


def adapt_solution_candidate(
    value: Mapping[str, Any], environment: ImplicitSolutionEnvironment,
    species: tuple[SolutionSpecies, ...], target_smiles: str,
) -> SolutionCandidateReaction:
    kind = value.get("kind")
    if kind == "storca_reaction_candidate":
        candidate = _native_candidate(value, environment, species, target_smiles)
        _validate_verification_uncertainty(candidate)
        return _validate_candidate_route(candidate, species)
    if kind == "rmg_solution_candidate_adapter":
        route = value.get("rmg_route") or {}
        native = candidate_from_rmg_route(target_smiles, route)
        if native is None:
            raise ValueError("RMG solution candidate lacks resolved endpoint identities")
        native_value = native.as_dict()
        native_value.update({
            "solution_rate_prior": value.get("solution_rate_prior") or {},
            "solution_verification_uncertainty": value.get("solution_verification_uncertainty") or {},
            "estimated_cost_core_hours": value.get("estimated_cost_core_hours", 100.0),
            "applicability_flags": value.get("applicability_flags") or {},
        })
        candidate = _native_candidate(native_value, environment, species, target_smiles)
        _validate_verification_uncertainty(candidate)
        return _validate_candidate_route(candidate, species)
    candidate = SolutionCandidateReaction.from_dict(value)
    interval = candidate.approximate_rate
    if interval.use_scope == "prediction_bound" and (
        not math.isclose(interval.temperature_K, environment.temperature_K, rel_tol=0.0, abs_tol=1e-8)
        or interval.solvent.casefold() != environment.solvent.casefold()
        or interval.solvent_model != environment.solvent_model
    ):
        raise ValueError("Candidate prediction interval does not match the session environment")
    _validate_verification_uncertainty(candidate)
    return _validate_candidate_route(candidate, species)


def _validate_verification_uncertainty(candidate: SolutionCandidateReaction) -> None:
    if candidate.active_applicability_flags:
        return
    uncertainty = candidate.uncertainty
    try:
        barrier = float(uncertainty["barrier_kcal_mol"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"Candidate {candidate.candidate_id} needs a composite barrier uncertainty before allocation"
        ) from error
    components = {str(item) for item in uncertainty.get("components") or ()}
    if (
        not math.isfinite(barrier) or barrier <= 0.0
        or not str(uncertainty.get("source") or "").strip()
        or not REQUIRED_UNCERTAINTY_COMPONENTS <= components
    ):
        raise ValueError(
            f"Candidate {candidate.candidate_id} has incomplete solution-verification uncertainty"
        )


def _validate_candidate_route(
    candidate: SolutionCandidateReaction,
    species: tuple[SolutionSpecies, ...],
) -> SolutionCandidateReaction:
    by_id = {item.species_id: item for item in species}
    if candidate.reactant_id not in by_id or not set(candidate.products) <= set(by_id):
        raise ValueError(f"Candidate {candidate.candidate_id} route references unknown species")
    clean = _clean_route(candidate.route, parent_smiles=by_id[candidate.reactant_id].smiles)
    spec = RouteSpec(**clean)
    if len(spec.reactant_smiles) != 1 or canonical_smiles(spec.reactant_smiles[0]) != by_id[candidate.reactant_id].smiles:
        raise ValueError(f"Candidate {candidate.candidate_id} route reactant identity is inconsistent")
    expected_products = []
    for product_id, coefficient in candidate.products.items():
        expected_products.extend([by_id[product_id].smiles] * int(coefficient))
    if sorted(canonical_smiles(item) for item in spec.product_smiles) != sorted(expected_products):
        raise ValueError(f"Candidate {candidate.candidate_id} route product stoichiometry is inconsistent")
    if spec.charge != 0 or spec.multiplicity < 1:
        raise ValueError(f"Candidate {candidate.candidate_id} route lacks a neutral declared spin surface")
    return replace(candidate, route=clean)


def candidate_network_from_session_request(request: Mapping[str, Any]) -> SolutionCandidateNetwork:
    if request.get("kind") != "implicit_solution_verification_session_request":
        raise ValueError("Request kind must be 'implicit_solution_verification_session_request'")
    environment = ImplicitSolutionEnvironment.from_dict(request.get("environment") or {})
    context = request.get("network_context") or {}
    species = tuple(SolutionSpecies.from_dict(item) for item in context.get("species") or ())
    target_smiles = str(request.get("target_smiles") or "").strip()
    target_id = str(context.get("target_species_id") or "").strip()
    by_id = {item.species_id: item for item in species}
    if target_id not in by_id or canonical_smiles(target_smiles) != by_id[target_id].smiles:
        raise ValueError("Session target_smiles does not match network_context.target_species_id")
    verified = tuple(
        VerifiedSolutionReaction.from_dict(item) for item in context.get("verified_reactions") or ()
    )
    candidates = tuple(
        adapt_solution_candidate(item, environment, species, target_smiles)
        for item in request.get("candidate_reactions") or ()
    )
    if not verified and not candidates:
        raise ValueError("Solution verification session has no reactions to assess")
    return SolutionCandidateNetwork(
        environment=environment, target_smiles=target_smiles,
        target_species_id=target_id, species=species,
        verified_reactions=verified, candidates=candidates,
        coverage_status=str(context.get("coverage_status") or "unresolved"),
        coverage_claim_id=context.get("coverage_claim_id"),
        discovery_sources=tuple(context.get("discovery_sources") or ()),
        unresolved_discovery_ids=tuple(context.get("unresolved_discovery_ids") or ()),
    )


__all__ = ["adapt_solution_candidate", "candidate_network_from_session_request", "run_solution_candidate_rate_coefficients"]
