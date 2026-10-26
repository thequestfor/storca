"""Solvent-consistent ORCA path evidence for the first solution lane.

The verified stability wrapper remains restricted to neutral first-order
networks. A separate coefficient-only request can calculate both directions
of a neutral one-/two-molecule candidate without concentrations or an invented
error bound. Both reuse STORCA's endpoint/NEB/TS/frequency/IRC engine and bind
every ORCA input to one immutable CPCM or SMD environment.
Only a verified saddle with condition- and recipe-matched thermochemistry can
be converted to a conventional solution-Eyring rate. Network completeness is a separate
claim supplied by reaction discovery; path verification never invents it.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

from .methods import harmonic_method_profile
from .conformer_rates import (
    CONFORMER_RATE_CONTRACT, run_reactant_conformer_rate_diagnostic, validate_conformer_production_scope,
)
from .orca_thermochemistry_recipe import (
    THERMOCHEMISTRY_RECIPE_CONTRACT, observed_orca_thermochemistry_recipe,
    matched_orca_thermochemistry_recipes,
)
from .reaction_candidates import canonical_smiles
from .reaction_path import run_generic_reaction_path, validated_orca_temperature_thermochemistry
from .route_thermochemistry import HARTREE_TO_J_MOL
from .route_verify import RouteSpec
from .solution_kinetics import (
    ImplicitSolutionEnvironment,
    PROHIBITED_MECHANISM_FLAGS,
    SolutionSpecies,
    VerifiedSolutionNetwork,
    VerifiedSolutionReaction,
)
from .solution_stability import run_implicit_solution_stability_bundle


BOLTZMANN_J_K = 1.380649e-23
PLANCK_J_S = 6.62607015e-34
GAS_CONSTANT_J_MOL_K = 8.31446261815324
KCAL_TO_J = 4184.0
# Legacy reference value only; production/session admission now requires an
# explicit, condition-matched single-temperature input/output evidence contract.
SUPPORTED_THERMOCHEMISTRY_TEMPERATURE_K = 298.15
THERMOCHEMISTRY_TEMPERATURE_CONTRACT = "explicit_single_temperature_v1"
REQUIRED_UNCERTAINTY_COMPONENTS = frozenset({
    "conformer_sampling",
    "electronic_structure_barrier",
    "harmonic_tst",
    "implicit_solvent_model",
})
_SOLVENT_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9,+-]*$")
_SOLVENT_TOKEN = re.compile(r"\b(?:CPCM|SMD)\s*\([^\n)]*\)", re.IGNORECASE)


@dataclass(frozen=True)
class _CoefficientConditions:
    """Physical rate conditions, deliberately without a fabricated trajectory."""

    solvent: str
    solvent_model: str
    temperature_K: float
    electronic_structure_profile: str
    standard_state_mol_L: float = 1.0
    light_condition: str = "dark"

    def as_environment_spec(self):
        from .lanes.base import EnvironmentSpec
        return EnvironmentSpec(
            "homogeneous_solution", "implicit_solution_selected_path_coefficients_v1",
            {name: getattr(self, name) for name in self.__dataclass_fields__},
        )


def _validate_coefficient_request(request, method_profile_name):
    """Internal candidate-path producer; not autonomous network discovery."""
    if isinstance(request.get("schema_version"), bool) or request.get("schema_version") != 1:
        raise ValueError("Coefficient request requires schema_version=1")
    value = dict(request.get("environment") or {})
    if any(key in request for key in ("network", "uncertainty", "candidate_reaction")):
        raise ValueError("Coefficient production uses a native candidate route and physical conditions, not a supplied network or error bound")
    allowed = set(_CoefficientConditions.__dataclass_fields__)
    if set(value) - allowed:
        raise ValueError("Coefficient conditions contain unsupported fields: " + ", ".join(sorted(set(value) - allowed)))
    temperature = value.get("temperature_K")
    if isinstance(temperature, bool) or not math.isfinite(float(temperature)) or float(temperature) <= 0:
        raise ValueError("Coefficient temperature must be finite and positive")
    if value.get("electronic_structure_profile") != method_profile_name:
        raise ValueError("Coefficient electronic_structure_profile must match the requested method")
    if value.get("solvent_model") not in {"CPCM", "SMD"} or not isinstance(value.get("solvent"), str) or not _SOLVENT_NAME.fullmatch(value["solvent"]):
        raise ValueError("Coefficient production requires one declared CPCM/SMD named solvent")
    standard = value.get("standard_state_mol_L", 1.0)
    if isinstance(standard, bool) or float(standard) != 1.0:
        raise ValueError("Coefficient production currently supports the one-molar standard state only")
    if value.get("light_condition", "dark") != "dark":
        raise ValueError("Conventional solution coefficient production requires dark conditions")
    environment = _CoefficientConditions(
        value["solvent"], value["solvent_model"], float(temperature), method_profile_name,
    )
    route = _route_spec(request.get("route") or {})
    if (not isinstance(route.charge, int) or isinstance(route.charge, bool) or route.charge != 0
            or not isinstance(route.multiplicity, int) or isinstance(route.multiplicity, bool) or route.multiplicity != 1):
        raise ValueError("Coefficient production requires a neutral singlet surface")
    if not str(route.route_id).strip():
        raise ValueError("Coefficient route needs an immutable reaction identity")
    if any(len(side) not in {1, 2} for side in (route.reactant_smiles, route.product_smiles)):
        raise ValueError("Coefficient production supports one or two separated molecules on each side")
    flags = dict(request.get("applicability_flags") or {})
    if any(not isinstance(flag, bool) for flag in flags.values()) or set(flags) - PROHIBITED_MECHANISM_FLAGS:
        raise ValueError("Coefficient applicability flags must be known boolean diagnostics")
    if any(active for name, active in flags.items() if name != "bimolecular_encounter"):
        raise ValueError("Coefficient candidate requires a specialized method lane")
    from rdkit import Chem
    identities = {}
    side_labels = []
    side_counts = []
    for smiles_values, labels, charges, spins in (
        (route.reactant_smiles, route.reactant_labels, route.reactant_charges, route.reactant_multiplicities),
        (route.product_smiles, route.product_labels, route.product_charges, route.product_multiplicities),
    ):
        if charges and (len(charges) != len(smiles_values) or any(isinstance(q, bool) or q != 0 for q in charges)):
            raise ValueError("Coefficient fragments must have neutral declared charges")
        if spins and (len(spins) != len(smiles_values) or any(isinstance(s, bool) or s != 1 for s in spins)):
            raise ValueError("Coefficient fragments must be declared singlets")
        if labels and len(labels) != len(smiles_values):
            raise ValueError("Coefficient fragment labels do not match stoichiometry")
        selected_labels = tuple(labels or tuple(canonical_smiles(s) for s in smiles_values))
        if any(not isinstance(label, str) or not label or label != label.strip() for label in selected_labels):
            raise ValueError("Coefficient species labels must be nonempty canonical strings without surrounding whitespace")
        counts = {}
        for smiles, label in zip(smiles_values, selected_labels):
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
                raise ValueError("Each coefficient fragment must identify one molecule")
            if any(atom.GetFormalCharge() or atom.GetNumRadicalElectrons() or atom.GetIsotope() for atom in molecule.GetAtoms()):
                raise ValueError("Coefficient fragments require neutral closed-shell ordinary-isotope identities")
            label = str(label).strip()
            identity = canonical_smiles(smiles)
            if not label or label in identities and identities[label] != identity:
                raise ValueError("Coefficient species label has conflicting molecular identities")
            identities[label] = identity
            for atom in Chem.AddHs(molecule).GetAtoms():
                element = atom.GetSymbol()
                counts[element] = counts.get(element, 0) + 1
        side_labels.append(selected_labels)
        side_counts.append(counts)
    if side_counts[0] != side_counts[1]:
        raise ValueError("Coefficient route is not element balanced")
    if route.surface_multiplicities and tuple(route.surface_multiplicities) != (1,):
        raise ValueError("Coefficient production does not resolve multiple electronic surfaces")
    # Reuse the path engine's actual physical-model policy before allocation;
    # a gas-discovery route is not silently relabeled as solution evidence.
    from .reaction_path import _surface_environment_applicability
    from .surface_resolution import SurfaceResolutionPolicy
    applicability = _surface_environment_applicability(
        route, {"phase_model": "homogeneous implicit solution", "solvent": environment.solvent,
                "solvent_model": environment.solvent_model}, SurfaceResolutionPolicy(),
        solution_orca_keywords(environment, harmonic_method_profile(method_profile_name)),
    )
    if applicability["status"] != "applicable":
        raise ValueError("Coefficient route environment conflicts with the physical solution request: " + ", ".join(applicability["blockers"]))
    candidate = {"reaction_id": route.route_id, "applicability_flags": flags,
                 "side_labels": side_labels, "species_identities": identities}
    return environment, {}, candidate, route, {}


def solution_eyring_coefficient(
    barrier_hartree: float, temperature_K: float, *, molecularity: int = 1,
    thermochemistry_pressure_atm: float = 1.0, standard_state_mol_L: float = 1.0,
) -> dict[str, Any]:
    """Dimensional selected-saddle TST arithmetic, not path or error evidence.

    Per-species conversion is RT log(C° RT/p). The activation correction is
    (1-n) times that conversion, not +1.89 kcal/mol for every reaction or T.
    Concentrations use mol/L; rate means elementary reaction events, with no
    extra identical-reactant depletion/degeneracy factor silently inserted.
    """
    if not isinstance(molecularity, int) or isinstance(molecularity, bool) or molecularity not in {1, 2}:
        raise ValueError("Solution coefficient molecularity must be one or two")
    values = (barrier_hartree, temperature_K, thermochemistry_pressure_atm, standard_state_mol_L)
    if any(isinstance(v, bool) or not math.isfinite(float(v)) for v in values):
        raise ValueError("Solution coefficient parameters must be finite numbers")
    barrier, temperature, pressure, standard = map(float, values)
    if temperature <= 0 or pressure <= 0 or standard <= 0:
        raise ValueError("Solution coefficient temperature, reference pressure and concentration must be positive")
    RT = GAS_CONSTANT_J_MOL_K * temperature
    per_species = RT * math.log(1000.0 * standard * RT / (pressure * 101325.0))
    correction = (1 - molecularity) * per_species
    barrier_J_mol = barrier * HARTREE_TO_J_MOL + correction
    if barrier_J_mol <= 0:
        raise ValueError("Nonpositive concentration-standard-state barrier requires a specialized rate model")
    logarithm = math.log(BOLTZMANN_J_K * temperature / PLANCK_J_S) + (1 - molecularity) * math.log(standard) - barrier_J_mol / RT
    if not -744.0 < logarithm < 709.0:
        raise ValueError("The Eyring coefficient is outside finite positive numerical range")
    rate = math.exp(logarithm)
    return {
        "rate_constant": rate, "rate_units": "s^-1" if molecularity == 1 else "L mol^-1 s^-1",
        **({"rate_constant_s": rate} if molecularity == 1 else {}),
        "molecularity": molecularity, "temperature_K": temperature,
        "uncorrected_activation_free_energy_hartree": barrier,
        "activation_free_energy_hartree": barrier_J_mol / HARTREE_TO_J_MOL,
        "activation_free_energy_kcal_mol": barrier_J_mol / KCAL_TO_J,
        "transmission_coefficient": 1.0,
        "rate_convention": "elementary_reaction_event_mass_action",
        "standard_state": {
            "name": "one_molar_solution" if standard == 1.0 else "explicit_concentration_standard_state",
            "standard_state_mol_L": standard, "thermochemistry_pressure_atm": pressure,
            "per_species_correction_J_mol": per_species, "correction_J_mol": correction,
            "correction_kcal_mol": correction / KCAL_TO_J,
            "equation": "DeltaG_dagger_C = G_TS - sum(G_i) + (1-n)*R*T*ln(1000*C_standard*R*T/(p_atm*101325))",
        },
    }


def _write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


def _canonical_hash(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False, default=str,
    )
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _file_hash(path: Path) -> str:
    return "sha256:" + hashlib.sha256(Path(path).read_bytes()).hexdigest()


def solution_orca_keywords(environment: ImplicitSolutionEnvironment, profile: Mapping[str, Any]) -> list[str]:
    """Return a safe, single-solvent ORCA keyword contract."""
    solvent = environment.solvent
    if not _SOLVENT_NAME.fullmatch(solvent):
        raise ValueError(
            "The autonomous solution path producer currently accepts ORCA solvent names containing "
            "letters, numbers, commas, plus signs, or hyphens only"
        )
    keywords = [str(item).strip() for item in profile.get("orca_keywords") or ()]
    if not keywords or any(not item for item in keywords):
        raise ValueError("The selected ORCA method profile has no valid keyword contract")
    joined = " ".join(keywords)
    if _SOLVENT_TOKEN.search(joined) or "%cpcm" in joined.casefold():
        raise ValueError("The base method profile must not contain a second solvent model")
    return [*keywords, f"{environment.solvent_model}({solvent})"]


def _route_spec(value: Mapping[str, Any]) -> RouteSpec:
    if not isinstance(value, Mapping):
        raise TypeError("route must be a mapping")
    allowed = {item.name for item in fields(RouteSpec)}
    unknown = set(value) - allowed
    if unknown:
        raise ValueError("Unknown route fields: " + ", ".join(sorted(unknown)))
    tuple_fields = {
        "reactant_smiles", "product_smiles", "reactant_labels", "product_labels",
        "reactant_charges", "product_charges", "reactant_multiplicities",
        "product_multiplicities", "surface_multiplicities", "broken_bonds",
        "formed_bonds", "atom_mapping",
    }
    normalized = dict(value)
    for name in tuple_fields:
        if name in normalized:
            normalized[name] = tuple(
                tuple(item) if name in {"broken_bonds", "formed_bonds", "atom_mapping"} else item
                for item in (normalized[name] or ())
            )
    return RouteSpec(**normalized)


def _expanded_products(products: Mapping[str, Any], species: Mapping[str, SolutionSpecies]) -> list[str]:
    expanded: list[str] = []
    for species_id, coefficient_value in products.items():
        coefficient = int(coefficient_value)
        if coefficient <= 0 or species_id not in species:
            raise ValueError("Candidate products need positive coefficients and known species IDs")
        expanded.extend([species[species_id].smiles] * coefficient)
    if not expanded:
        raise ValueError("Candidate reaction has no products")
    return expanded


def _validate_request(
    request: Mapping[str, Any], method_profile_name: str,
    *, require_complete_network_coverage: bool = True,
) -> tuple[ImplicitSolutionEnvironment, dict[str, Any], dict[str, Any], RouteSpec, dict[str, Any]]:
    if not isinstance(request, Mapping):
        raise TypeError("Solution path request must be a mapping")
    if request.get("kind") != "implicit_solution_path_verification_request":
        raise ValueError("Request kind must be 'implicit_solution_path_verification_request'")
    environment = ImplicitSolutionEnvironment.from_dict(request.get("environment") or {})
    if environment.electronic_structure_profile != method_profile_name:
        raise ValueError(
            "environment.electronic_structure_profile must exactly match --method-profile"
        )
    if not math.isclose(environment.standard_state_mol_L, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("The solution path producer supports the one-molar solution standard state only")

    network = dict(request.get("network") or {})
    species_records = network.get("species") or []
    species_values = [SolutionSpecies.from_dict(item) for item in species_records]
    species = {item.species_id: item for item in species_values}
    if not species_values or len(species) != len(species_values):
        raise ValueError("Network species IDs must be nonempty and unique")
    if any(item.charge != 0 for item in species_values):
        raise ValueError("The neutral implicit-solution path producer rejects ionic species")
    target_id = str(network.get("target_species_id") or "").strip()
    if target_id not in species:
        raise ValueError("network.target_species_id is missing from the species records")
    target_smiles = str(request.get("target_smiles") or "").strip()
    if not target_smiles or canonical_smiles(target_smiles) != species[target_id].smiles:
        raise ValueError("target_smiles does not match the network target species")

    candidate = dict(request.get("candidate_reaction") or {})
    reaction_id = str(candidate.get("reaction_id") or "").strip()
    reactant_id = str(candidate.get("reactant_id") or "").strip()
    if not reaction_id or reactant_id not in species:
        raise ValueError("Candidate reaction needs an ID and one known reactant species")
    molecularity = candidate.get("molecularity", 1)
    if isinstance(molecularity, bool) or int(molecularity) != 1:
        raise ValueError("The first solution path producer supports molecularity=1 only")
    products = dict(candidate.get("products") or {})
    product_smiles = _expanded_products(products, species)
    flags = {str(key): bool(value) for key, value in dict(candidate.get("applicability_flags") or {}).items()}
    unknown_flags = set(flags) - PROHIBITED_MECHANISM_FLAGS
    if unknown_flags:
        raise ValueError("Unknown applicability flags: " + ", ".join(sorted(unknown_flags)))
    active_flags = sorted(key for key, active in flags.items() if active)
    if active_flags:
        raise ValueError(
            "Candidate requires a specialized method lane: " + ", ".join(active_flags)
        )

    uncertainty = dict(candidate.get("uncertainty") or {})
    try:
        barrier_uncertainty = float(uncertainty["barrier_kcal_mol"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("candidate_reaction.uncertainty.barrier_kcal_mol is required") from error
    if not math.isfinite(barrier_uncertainty) or barrier_uncertainty <= 0.0:
        raise ValueError("Composite barrier uncertainty must be finite and positive")
    uncertainty_source = str(uncertainty.get("source") or "").strip()
    components = tuple(sorted({str(item).strip() for item in uncertainty.get("components") or () if str(item).strip()}))
    missing_components = REQUIRED_UNCERTAINTY_COMPONENTS - set(components)
    if not uncertainty_source or missing_components:
        raise ValueError(
            "Rate uncertainty needs a source and components covering: "
            + ", ".join(sorted(REQUIRED_UNCERTAINTY_COMPONENTS))
        )
    candidate["reaction_id"] = reaction_id
    candidate["reactant_id"] = reactant_id
    candidate["products"] = products
    candidate["applicability_flags"] = flags
    candidate["uncertainty"] = {
        "barrier_kcal_mol": barrier_uncertainty,
        "source": uncertainty_source,
        "components": list(components),
    }

    existing_reactions = network.get("reactions") or []
    if any(str(item.get("reaction_id")) == reaction_id for item in existing_reactions):
        raise ValueError("candidate_reaction must not duplicate a retained verified reaction ID")
    verified_existing = [VerifiedSolutionReaction.from_dict(item) for item in existing_reactions]
    coverage_complete = (
        str(network.get("coverage_status"))
        == "verified_complete_within_declared_unimolecular_scope"
    )
    if require_complete_network_coverage and not coverage_complete:
        raise ValueError(
            "Path verification cannot create mechanism completeness; an independent verified restricted "
            "network-coverage claim is required"
        )
    if require_complete_network_coverage and (
        not network.get("coverage_claim_id") or not network.get("discovery_sources")
    ):
        raise ValueError("Verified restricted coverage needs a claim ID and discovery sources")
    if require_complete_network_coverage and network.get("unresolved_candidate_ids"):
        raise ValueError("Unresolved discovery candidates block emission of an admissible solution bundle")

    route = _route_spec(request.get("route") or {})
    if len(route.reactant_smiles) != 1:
        raise ValueError("The solution path route must have exactly one molecular reactant")
    if canonical_smiles(route.reactant_smiles[0]) != species[reactant_id].smiles:
        raise ValueError("Route reactant identity does not match candidate_reaction.reactant_id")
    if canonical_smiles(route.parent_smiles) != species[reactant_id].smiles:
        raise ValueError("Route parent_smiles must identify the unimolecular reactant")
    if sorted(canonical_smiles(item) for item in route.product_smiles) != sorted(product_smiles):
        raise ValueError("Route product identities/stoichiometry do not match candidate_reaction.products")
    if route.charge != 0 or route.multiplicity < 1:
        raise ValueError("The neutral solution path route needs charge=0 and a declared multiplicity")
    validate_conformer_production_scope(route)

    # Build a provisional structurally valid reaction only to invoke the native
    # balance and monotonic-target-loss gates before any ORCA allocation.  Its
    # rate is explicitly a validation sentinel and is never retained.
    sentinel = VerifiedSolutionReaction(
        reaction_id=reaction_id, reactant_id=reactant_id, products=products,
        rate_constant_s=1.0, temperature_K=environment.temperature_K,
        solvent=environment.solvent, solvent_model=environment.solvent_model,
        standard_state_mol_L=1.0,
        verification_status="verified_experimental_solution_rate",
        rate_theory="experimental_solution_rate", source="preflight balance sentinel; not evidence",
        rate_claim_id="preflight-sentinel", lower_rate_multiplier=1.0,
        upper_rate_multiplier=1.0, uncertainty_source="preflight sentinel",
        uncertainty_components=("preflight_only",), applicability_flags=flags,
    )
    provisional_network = {
        **network,
        "species": [item.as_dict() for item in species_values],
        "reactions": [item.as_dict() for item in verified_existing] + [sentinel.as_dict()],
    }
    VerifiedSolutionNetwork.from_dict(provisional_network)
    return environment, network, candidate, route, {
        "species": species,
        "verified_existing": verified_existing,
    }


def _audit_orca_inputs(workdir: Path, expected_keyword: str) -> dict[str, Any]:
    records = []
    for path in sorted(Path(workdir).rglob("*.inp")):
        text = path.read_text(errors="replace")
        tokens = _SOLVENT_TOKEN.findall(text)
        if len(tokens) != 1 or tokens[0].casefold() != expected_keyword.casefold():
            raise ValueError(
                f"Mixed or missing implicit-solvent contract in ORCA input: {path}"
            )
        if "%cpcm" in text.casefold():
            raise ValueError(f"An undeclared %cpcm override was found in ORCA input: {path}")
        records.append({
            "path": str(path), "sha256": _file_hash(path),
            "solvent_keyword": tokens[0],
        })
    if not records:
        raise ValueError("Verified path result retained no auditable ORCA input files")
    return {
        "status": "passed", "expected_solvent_keyword": expected_keyword,
        "input_count": len(records), "inputs": records,
    }


def _condition_matched_thermochemistry(
    path_result: Mapping[str, Any], environment: ImplicitSolutionEnvironment,
) -> tuple[dict[str, Any], dict[str, Any], float]:
    verification = path_result.get("route_verification") or {}
    stationary = verification.get("stationary_points") or {}
    reactants = stationary.get("reactants") or []
    if len(reactants) != 1:
        raise ValueError("Solution Eyring rate generation requires one validated reactant thermochemistry record")
    # For a bound isomerization, reference the exact optimized reactant end of
    # the selected IRC/NEB orientation.  The generic engine also retains a
    # separately generated stationary-species record for Arkane; using that
    # record first could silently mix conformers in the activation barrier.
    endpoint_reactant = (verification.get("endpoint_validation") or {}).get("reactant") or {}
    reactant_output = endpoint_reactant.get("frequency_output")
    reactant_reference = "selected_path_endpoint"
    if not reactant_output:
        # Unimolecular dissociation uses assembled separated-fragment channel
        # endpoints, so its independently validated reactant species is the
        # physical thermochemical reference.
        reactant_output = reactants[0].get("frequency_output")
        reactant_reference = "validated_stationary_reactant"
    ts_output = (verification.get("transition_state") or {}).get("frequency_output")
    if not reactant_output or not ts_output:
        raise ValueError("Reactant and transition-state frequency outputs are required")
    reactant = _condition_matched_frequency_record(reactant_output, environment.temperature_K)
    transition_state = _condition_matched_frequency_record(ts_output, environment.temperature_K)
    matched_orca_thermochemistry_recipes(
        reactant["thermochemistry_recipe"], transition_state["thermochemistry_recipe"],
    )
    reactant["activation_reference_role"] = reactant_reference
    for label, record in (("reactant", reactant), ("transition_state", transition_state)):
        if not record.get("valid"):
            raise ValueError(f"{label} ORCA thermochemistry is incomplete")
    barrier_hartree = (
        float(transition_state["gibbs_free_energy_hartree"])
        - float(reactant["gibbs_free_energy_hartree"])
    )
    if not math.isfinite(barrier_hartree) or barrier_hartree <= 0.0:
        raise ValueError("Conventional TST requires a positive condition- and recipe-matched free-energy barrier")
    return reactant, transition_state, barrier_hartree


def _condition_matched_frequency_record(output, temperature):
    record = validated_orca_temperature_thermochemistry(Path(output), temperature)
    recipe = observed_orca_thermochemistry_recipe(Path(output))
    if (recipe["source_output"]["sha256"] != record["source_output_sha256"]
            or recipe["source_input"]["sha256"] != record["temperature_input_evidence"]["source_input_sha256"]):
        raise ValueError("ORCA thermochemistry changed between temperature and recipe validation")
    return {**record, "thermochemistry_recipe": recipe}


def _coefficient_path_evidence(base_result, route, candidate, environment, method_keywords, workdir):
    """Calculate both directions of the native producer's verified path.

    Called only after run_generic_reaction_path, never as a public admission
    function for user-authored preverified JSON. Path/normal-mode/IRC validity
    remains the existing path producer's responsibility. All frequency source
    files are reopened here; numerical rates carry no network admission credit.
    """
    verification = base_result["path_result"]["route_verification"]
    stationary = verification.get("stationary_points") or {}
    imaginary = ((verification.get("transition_state") or {}).get("frequency_check") or {}).get("NumImaginary")
    if ((verification.get("endpoint_validation") or {}).get("valid") is not True
            or ((verification.get("irc") or {}).get("endpoint_validation") or {}).get("valid") is not True
            or not isinstance(imaginary, int) or isinstance(imaginary, bool) or imaginary != 1
            or stationary.get("valid") is not True or stationary.get("status") != "validated"):
        raise ValueError("Coefficient producer lacks validated native endpoints, minima, saddle or IRC diagnostics")
    from .reaction_path import _validate_retained_orca_job_input
    native_contract_count = 0
    for audited in base_result["orca_input_audit"]["inputs"]:
        input_path = Path(audited["path"])
        if input_path.with_suffix(".contract.json").is_file():
            _validate_retained_orca_job_input(input_path)
            native_contract_count += 1
    artifacts = _artifact_hashes(workdir)
    sources = {str(Path(item["path"]).resolve()): item["sha256"] for item in artifacts}
    for audited in base_result["orca_input_audit"]["inputs"]:
        if sources.get(str(Path(audited["path"]).resolve())) != audited["sha256"]:
            raise ValueError("Already-audited solution input changed before the coefficient source snapshot")
    def bind_native_references(value):
        if isinstance(value, Mapping):
            for key, item in value.items():
                if key == "source_references" and isinstance(item, Mapping):
                    for reference in item.values():
                        if not isinstance(reference, Mapping) or not reference.get("path") or not reference.get("sha256"):
                            raise ValueError("Native path source reference is incomplete")
                        path = str(Path(reference["path"]).resolve())
                        digest = "sha256:" + str(reference["sha256"]).removeprefix("sha256:")
                        if _file_hash(Path(path)) != digest or path in sources and sources[path] != digest:
                            raise ValueError("Native path source reference changed before coefficient admission")
                        sources[path] = digest
                else:
                    bind_native_references(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                bind_native_references(item)
    bind_native_references(verification)
    expected_tokens = sorted(token.casefold() for keyword in method_keywords for token in keyword.split())
    records = []
    def read_frequency(output, role):
        record = _condition_matched_frequency_record(output, environment.temperature_K)
        recipe = record["thermochemistry_recipe"]
        input_path = Path(recipe["source_input"]["path"])
        text = re.sub(r"#[^\n]*", "", input_path.read_text())
        if set(name.casefold() for name in re.findall(r"%([A-Za-z]\w*)", text)) - {"pal", "freq", "maxcore"}:
            raise ValueError("Coefficient frequency contains non-native blocks that can override the electronic-structure model")
        frequency_blocks = re.findall(r"%freq\b(.*?)\bend\b", text, flags=re.IGNORECASE | re.DOTALL)
        allowed_frequency_controls = {"temp", "quasirrho", "qrrhoreffreq", "cutofffreq", "scalfreq", "pressure", "projecttr"}
        for block in frequency_blocks:
            for line in block.splitlines():
                if line.strip() and line.split()[0].casefold() not in allowed_frequency_controls:
                    raise ValueError("Coefficient frequency contains an undeclared thermochemistry or isotope override")
        tokens = " ".join(re.findall(r"(?m)^\s*!([^\n]*)", text)).casefold().split()
        frequency_tokens = {"freq", "anfreq", "numfreq"}
        if len([token for token in tokens if token in frequency_tokens]) != 1 or sorted(
            token for token in tokens if token not in frequency_tokens | {"tightscf"}
        ) != sorted(token for token in expected_tokens if token != "tightscf"):
            raise ValueError("Coefficient frequency electronic-structure method differs from the requested profile")
        states = re.findall(r"(?mi)^\s*\*\s+xyz(?:file)?\s+(-?\d+)\s+(\d+)\b", text)
        if states != [("0", "1")]:
            raise ValueError("Coefficient frequency input charge or multiplicity differs from the neutral singlet contract")
        if "%cpcm" in text.casefold():
            raise ValueError("Coefficient frequency contains an undeclared solvent override")
        for key in ("source_input", "source_output"):
            source = recipe[key]
            path = str(Path(source["path"]).resolve())
            digest = "sha256:" + source["sha256"]
            if path in sources and sources[path] != digest:
                raise ValueError("Coefficient source changed during thermochemistry parsing")
            sources[path] = digest
        record["activation_reference_role"] = role
        records.append(record)
        return record
    sides = []
    for side, smiles_values in (("reactants", route.reactant_smiles), ("products", route.product_smiles)):
        native_records = stationary.get(side) or []
        if len(native_records) != len(smiles_values):
            raise ValueError("Coefficient stationary-species records do not match separated stoichiometry")
        for native, smiles in zip(native_records, smiles_values):
            if (native.get("status") != "validated_minimum"
                    or (native.get("frequency_check") or {}).get("IsMinimum") is not True
                    or canonical_smiles(str(native.get("smiles") or "")) != canonical_smiles(smiles)
                    or isinstance(native.get("charge"), bool) or native.get("charge") != 0
                    or isinstance(native.get("multiplicity"), bool) or native.get("multiplicity") != 1):
                raise ValueError("Coefficient stationary species identity/state/minimum evidence is unresolved")
        endpoint = (verification["endpoint_validation"].get(side[:-1]) or {})
        endpoint_output = endpoint.get("frequency_output") if len(smiles_values) == 1 else None
        # Selected bound endpoints avoid conformer substitution. Two-fragment
        # sides instead reference separated minima, never an encounter complex.
        side_records = [read_frequency(
            endpoint_output or native["frequency_output"],
            "selected_path_endpoint" if endpoint_output else "validated_stationary_separated_species",
        ) for native in native_records]
        sides.append(side_records)
    ts_record = read_frequency(verification["transition_state"]["frequency_output"], "shared_verified_transition_state")
    for record in records:
        matched_orca_thermochemistry_recipes(record["thermochemistry_recipe"], ts_record["thermochemistry_recipe"])
    pressures = [record["thermochemistry_recipe"]["explicit_input_overrides"].get(
        "pressure", record["thermochemistry_recipe"]["recipe"]["thermochemistry_pressure_atm"],
    ) for record in records]
    if any(not math.isclose(value, pressures[0], rel_tol=0.0, abs_tol=1e-12) for value in pressures):
        raise ValueError("Exact frequency thermochemistry pressure contracts differ despite printed rounding")
    pressure = pressures[0]
    method_contract = base_result["method_contract"]
    method_contract["thermochemistry_recipe_id"] = ts_record["thermochemistry_recipe"]["recipe_id"]
    method_contract["observed_thermochemistry_recipe"] = ts_record["thermochemistry_recipe"]["recipe"]
    for path, digest in sorted(sources.items()):
        if _file_hash(Path(path)) != digest:
            raise ValueError("Coefficient calculation source changed before claim admission")
    artifacts = [{"path": path, "sha256": digest} for path, digest in sorted(sources.items())]
    calculation_hashes = sorted({item["sha256"] for item in artifacts
                               if Path(item["path"]).suffix in {".inp", ".out", ".xyz", ".hess"}})
    path_claim = {
        "schema_version": 1, "kind": "scientific_claim", "claim_type": "verified_solution_first_order_saddle",
        "status": "verified_within_declared_model", "reaction_id": candidate["reaction_id"], "route_id": route.route_id,
        "claim": "The native path producer verified a first-order saddle connecting the declared channels on one implicit-solution surface.",
        "method_contract": method_contract,
        "evidence": {key: verification[key] for key in ("endpoint_validation", "transition_state", "irc")},
        "limitations": ["Path verification does not establish rate accuracy, conformer coverage, encounter equilibrium or network completeness."],
    }
    path_claim["evidence"].update(orca_input_audit=base_result["orca_input_audit"], calculation_artifacts=artifacts)
    path_claim["evidence"]["native_input_contract_rechecks"] = {
        "checked_contract_count": native_contract_count,
        "legacy_inputs_without_contract_count": base_result["orca_input_audit"]["input_count"] - native_contract_count,
        "interpretation": "Existing input/reference contracts rechecked read-only; absent legacy contracts receive no invented historical proof. Native path diagnostics remain upstream.",
    }
    path_claim["claim_id"] = "claim:path:" + _canonical_hash(path_claim).split(":", 1)[1]
    packets = []
    sums = [math.fsum(record["gibbs_free_energy_hartree"] for record in side) for side in sides]
    ts_gibbs = ts_record["gibbs_free_energy_hartree"]
    for index, direction in enumerate(("forward", "reverse")):
        coefficient = solution_eyring_coefficient(
            ts_gibbs - sums[index], environment.temperature_K, molecularity=len(sides[index]),
            thermochemistry_pressure_atm=pressure, standard_state_mol_L=environment.standard_state_mol_L,
        )
        rate_claim = {
            "schema_version": 1, "kind": "scientific_claim", "claim_type": "solution_eyring_rate",
            "status": "estimated_with_unresolved_uncertainty", "reaction_id": candidate["reaction_id"] + ":" + direction,
            "path_claim_id": path_claim["claim_id"], "direction": direction,
            "claim": "A conditional selected-saddle conventional solution-TST coefficient, not an ensemble, transport-corrected or calibrated rate.",
            **coefficient,
            "thermochemistry": {"reactants": sides[index], "products": sides[1 - index], "transition_state": ts_record},
            "thermochemistry_recipe": ts_record["thermochemistry_recipe"],
            "uncertainty": {"status": "calibration_required", "lower_rate_multiplier": None, "upper_rate_multiplier": None,
                            "unresolved_components": sorted(REQUIRED_UNCERTAINTY_COMPONENTS | {"transition_state_search", "tunneling", "path_degeneracy", "diffusion_and_encounter", "solution_activities"})},
            "conformer_assessment": {"status": "single_selected_path_conditional_scenario", "network_rate_admission_allowed": False,
                                     "reason": "Endpoint, encounter, interconversion and TS conformer ensembles are not established by this coefficient producer."},
            "limitations": [
                "Conventional TST with transmission coefficient one, no tunneling or recrossing correction.",
                "Implicit-solvent energies with the observed ORCA ideal-gas/QRRHO recipe and an explicit concentration-standard-state conversion.",
                "Named solvent parameters are fixed; temperature-dependent solvent properties, isotope effects and activities are not validated.",
                "No diffusion ceiling/correction or summed path degeneracy; the coefficient is for this selected saddle only.",
                "No network completeness, lifetime, experimental calibration or environment-support promotion.",
            ],
        }
        rate_claim["claim_id"] = "claim:rate:" + _canonical_hash(rate_claim).split(":", 1)[1]
        labels = candidate["side_labels"]
        reactants = {label: labels[index].count(label) for label in set(labels[index])}
        products = {label: labels[1 - index].count(label) for label in set(labels[1 - index])}
        record = {
            "reaction_id": rate_claim["reaction_id"], "reactants": reactants, "products": products,
            "rate_constant": coefficient["rate_constant"], "rate_units": coefficient["rate_units"],
            **({"rate_constant_s": coefficient["rate_constant"]} if len(sides[index]) == 1 else {}),
            "temperature_K": environment.temperature_K, "solvent": environment.solvent, "solvent_model": environment.solvent_model,
            "standard_state_mol_L": 1.0, "standard_state_correction": "one_molar_solution",
            "verification_status": "verified_solution_tst", "rate_theory": "solution_eyring_tst",
            "rate_claim_id": rate_claim["claim_id"], "path_claim_id": path_claim["claim_id"],
            "calculation_hashes": calculation_hashes, "lower_rate_multiplier": None, "upper_rate_multiplier": None,
            "uncertainty_source": None, "network_rate_admission_allowed": False,
            "source": "STORCA native verified path; conditional selected-saddle coefficient with unresolved error",
        }
        packets.append({"schema_version": 1, "kind": "implicit_solution_rate_coefficient_production", "path_claim": path_claim,
                        "rate_claim": rate_claim, "verified_reaction": record, "network_rate_admission_allowed": False})
    # Both directions share the same TS and activity standard. Preserve units
    # in the log ratio; a second-order/first-order ratio is not dimensionless.
    per_species = packets[0]["rate_claim"]["standard_state"]["per_species_correction_J_mol"]
    reaction_J = (sums[1] - sums[0]) * HARTREE_TO_J_MOL + (len(sides[1]) - len(sides[0])) * per_species
    expected_log = -reaction_J / (GAS_CONSTANT_J_MOL_K * environment.temperature_K)
    actual_log = math.log(packets[0]["rate_claim"]["rate_constant"]) - math.log(packets[1]["rate_claim"]["rate_constant"])
    if not math.isclose(actual_log, expected_log, rel_tol=0, abs_tol=1e-9):
        raise ValueError("Forward/reverse coefficients violate their shared thermochemistry detailed balance")
    for path, digest in sources.items():
        if _file_hash(Path(path)) != digest:
            raise ValueError("Coefficient source changed during rate/claim assembly")
    base_result.update(
        status="conditional_rate_coefficients_emitted", rate_coefficients=packets, calculation_artifacts=artifacts,
        detailed_balance={"status": "consistent_with_shared_standard_state", "reaction_free_energy_J_mol": reaction_J,
                          "log_forward_reverse_ratio": actual_log,
                          "rate_ratio_units": "dimensionless" if len(sides[0]) == len(sides[1]) else "L mol^-1" if len(sides[0]) == 2 else "mol L^-1"},
        reason="Conditional forward/reverse selected-path coefficients retained; calibration, ensembles, transport and network coverage remain unresolved. No stability bundle was emitted.",
    )
    return base_result


def _artifact_hashes(workdir: Path, *, exclude: set[Path] | None = None) -> list[dict[str, str]]:
    excluded = {path.resolve() for path in (exclude or set())}
    records = []
    for path in sorted(item for item in Path(workdir).rglob("*") if item.is_file()):
        if path.resolve() in excluded or path.suffix == ".tmp":
            continue
        records.append({"path": str(path), "sha256": _file_hash(path)})
    return records


def _environment_bundle_dict(environment: ImplicitSolutionEnvironment) -> dict[str, Any]:
    return {
        "solvent": environment.solvent,
        "solvent_model": environment.solvent_model,
        "temperature_K": environment.temperature_K,
        "target_duration_seconds": environment.target_duration_seconds,
        "retention_fraction": environment.retention_fraction,
        "standard_state_mol_L": environment.standard_state_mol_L,
        "volume_L": environment.volume_L,
        "light_condition": environment.light_condition,
        "electronic_structure_profile": environment.electronic_structure_profile,
    }


def run_implicit_solution_path_verification(
    request: Mapping[str, Any],
    workdir: Path,
    *,
    execute: bool = False,
    ncores: int = 1,
    method_profile_name: str = "b3lyp-def2-svp",
    timeout_seconds: float | None = 14400.0,
    orientations: int = 3,
    neb_images: int = 8,
    require_complete_network_coverage: bool = True,
) -> dict[str, Any]:
    """Prepare/execute native solution paths, coefficients or a strict bundle.

    Coefficient requests contain physical conditions and an internally selected
    candidate route, but no network, concentrations, duration or invented error
    bound. They never authorize the legacy verified stability bundle.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    retained_request = workdir / "solution-path-request.json"
    result_path = workdir / "solution-path-verification.json"
    bundle_path = workdir / "verified-solution-bundle.json"
    coefficient_only = request.get("kind") == "implicit_solution_rate_coefficient_request" if isinstance(request, Mapping) else False
    try:
        requested_text = json.dumps(dict(request), indent=2, sort_keys=True, allow_nan=False) + "\n"
        if retained_request.is_file() and retained_request.read_text() != requested_text:
            # The returned diagnostic belongs to the attempted new request.
            # Do not replace the old request, successful report, bundle or jobs.
            return {
                "schema_version": 1, "kind": "implicit_solution_path_verification",
                "status": "invalid_or_incomplete_evidence",
                "reason": "Solution directory contains a different immutable request contract; use a new workspace",
                "verified_bundle": None, "result_json": None,
                "retained_result_json": str(result_path) if result_path.is_file() else None,
                "retained_request_sha256": _file_hash(retained_request),
                "requested_request_sha256": "sha256:" + hashlib.sha256(requested_text.encode()).hexdigest(),
            }
        if coefficient_only:
            environment, network, candidate, route, context = _validate_coefficient_request(request, method_profile_name)
        else:
            environment, network, candidate, route, context = _validate_request(
                request, method_profile_name,
                require_complete_network_coverage=require_complete_network_coverage,
            )
        profile = harmonic_method_profile(method_profile_name)
        method_keywords = solution_orca_keywords(environment, profile)
        solvent_keyword = method_keywords[-1]
        if not retained_request.is_file():
            with retained_request.open("x") as stream:
                stream.write(requested_text)
        method_contract = {
            "profile_name": method_profile_name,
            "method": profile["method"], "basis": profile["basis"],
            "orca_keywords": method_keywords,
            "solvent": environment.solvent,
            "solvent_model": environment.solvent_model,
            "thermochemistry_temperature_K": environment.temperature_K,
            "thermochemistry_temperature_control": "%freq Temp",
            "thermochemistry_temperature_contract": THERMOCHEMISTRY_TEMPERATURE_CONTRACT,
            "thermochemistry_recipe_contract": THERMOCHEMISTRY_RECIPE_CONTRACT,
            "conformer_rate_contract": "single_selected_path_conditional_v1" if coefficient_only else CONFORMER_RATE_CONTRACT,
            "solvent_temperature_treatment": "fixed_ORCA_named_solvent_parameters",
            "standard_state_mol_L": 1.0,
            "tunneling_correction": "none",
            "rate_theory": "solution_eyring_tst",
        }
        path_result = run_generic_reaction_path(
            route, workdir / "orca-path", prepare_only=not execute,
            ncores=ncores, method_keywords=method_keywords,
            timeout_seconds=timeout_seconds, orientations=orientations,
            neb_images=neb_images,
            thermochemistry_temperature_K=environment.temperature_K,
            condition_contract={
                "temperature_K": environment.temperature_K,
                **({"target_duration_seconds": environment.target_duration_seconds,
                    "retention_fraction": environment.retention_fraction} if not coefficient_only else {}),
                "phase_model": "homogeneous implicit solution",
                "solvent": environment.solvent,
                "solvent_model": environment.solvent_model,
                "standard_state_mol_L": 1.0,
            },
        )
        verification = path_result.get("route_verification") or {}
        base_result: dict[str, Any] = {
            "schema_version": 1,
            "kind": "implicit_solution_path_verification",
            "status": "prepared" if not execute else "path_unresolved",
            "request_sha256": _file_hash(retained_request),
            "environment": environment.as_environment_spec().as_dict(),
            "method_contract": method_contract,
            "route_id": route.route_id,
            "reaction_id": candidate["reaction_id"],
            "path_result": path_result,
            "verified_bundle": None,
            "result_json": str(result_path),
        }
        if not execute:
            base_result["reason"] = (
                "Endpoint structures were prepared; no ORCA path or rate claim was executed."
            )
            _write_json(result_path, base_result)
            return base_result
        if verification.get("status") != "verified_transition_state_path" or verification.get("path_classification") != "verified_barriered_path":
            base_result["reason"] = (
                "The solvent-consistent path did not establish an ordinary verified barriered saddle; "
                "no conventional solution-TST rate or bundle was emitted."
            )
            _write_json(result_path, base_result)
            return base_result

        input_audit = _audit_orca_inputs(workdir / "orca-path", solvent_keyword)
        base_result["orca_input_audit"] = input_audit
        if coefficient_only:
            base_result = _coefficient_path_evidence(
                base_result, route, candidate, environment, method_keywords, workdir / "orca-path",
            )
            for packet in base_result["rate_coefficients"]:
                claim_id = packet["rate_claim"]["claim_id"]
                packet_path = workdir / ("rate-" + packet["rate_claim"]["direction"] + "-" + claim_id.rsplit(":", 1)[-1][:16] + ".json")
                if packet_path.is_file():
                    if json.loads(packet_path.read_text()) != packet:
                        raise ValueError("Retained coefficient packet differs from its content-addressed claim")
                else:
                    _write_json(packet_path, packet)
                packet["production_json"] = str(packet_path)
                packet["production_sha256"] = _file_hash(packet_path)
            _write_json(result_path, base_result)
            return base_result
        reactant_thermo, ts_thermo, barrier_hartree = _condition_matched_thermochemistry(
            path_result, environment,
        )
        thermochemistry_recipe = matched_orca_thermochemistry_recipes(
            reactant_thermo["thermochemistry_recipe"], ts_thermo["thermochemistry_recipe"],
        )
        method_contract["thermochemistry_recipe_id"] = thermochemistry_recipe["recipe_id"]
        method_contract["observed_thermochemistry_recipe"] = thermochemistry_recipe["recipe"]
        conformer_assessment = run_reactant_conformer_rate_diagnostic(
            route, path_result, workdir / "orca-path" / "reactant-conformers",
            reactant_thermo=reactant_thermo, ts_thermo=ts_thermo, ncores=ncores,
            method_keywords=method_keywords, timeout_seconds=timeout_seconds,
        )
        base_result["conformer_assessment"] = conformer_assessment
        if not conformer_assessment.get("network_rate_admission_allowed"):
            base_result.update(status="conformer_rate_model_unresolved", reason=conformer_assessment["reason"])
            _write_json(result_path, base_result)
            return base_result
        # Conformer refinement creates inputs after the first path audit.
        # Audit the complete retained solution calculation set before credit.
        input_audit = _audit_orca_inputs(workdir / "orca-path", solvent_keyword)
        rate_s = solution_eyring_coefficient(barrier_hartree, environment.temperature_K)["rate_constant_s"]
        uncertainty = candidate["uncertainty"]
        uncertainty_exponent = (
            float(uncertainty["barrier_kcal_mol"]) * KCAL_TO_J
            / (GAS_CONSTANT_J_MOL_K * environment.temperature_K)
        )
        if uncertainty_exponent > 700.0:
            raise ValueError("The declared barrier uncertainty produces unbounded numerical rate multipliers")
        lower_multiplier = math.exp(-uncertainty_exponent)
        upper_multiplier = math.exp(uncertainty_exponent)

        artifacts = _artifact_hashes(workdir / "orca-path")
        for recipe_record in (thermochemistry_recipe["reactant"], thermochemistry_recipe["transition_state"]):
            for key in ("source_input", "source_output"):
                source = recipe_record[key]
                if _file_hash(Path(source["path"])) != "sha256:" + source["sha256"]:
                    raise ValueError("ORCA thermochemistry recipe source changed before rate admission")
        calculation_hashes = sorted({
            item["sha256"] for item in artifacts
            if Path(item["path"]).suffix.lower() in {".inp", ".out", ".xyz", ".hess"}
        })
        if not calculation_hashes:
            raise ValueError("Verified path retained no hashable ORCA calculation artifacts")
        path_claim = {
            "schema_version": 1,
            "kind": "scientific_claim",
            "claim_type": "verified_solution_first_order_saddle",
            "status": "verified_within_declared_model",
            "reaction_id": candidate["reaction_id"],
            "route_id": route.route_id,
            "claim": "A first-order saddle connects the declared reactant and products on one implicit-solution surface.",
            "evidence": {
                "endpoint_validation": verification.get("endpoint_validation"),
                "transition_state": verification.get("transition_state"),
                "irc": verification.get("irc"),
                "orca_input_audit": input_audit,
                "calculation_artifacts": artifacts,
            },
            "method_contract": method_contract,
            "limitations": [
                "One implicit-solvent electronic surface and one selected endpoint/path conformer are represented.",
                "A verified path does not establish reaction-network completeness.",
                "No claim is made for tunneling, explicit solvent participation, spin crossing, or multireference character.",
                "Named implicit-solvent parameters are held fixed; requesting thermal corrections at T does not validate solvent temperature dependence or phase stability.",
            ],
        }
        path_claim_id = "claim:path:" + _canonical_hash({
            "claim_type": path_claim["claim_type"],
            "reaction_id": candidate["reaction_id"],
            "route_id": route.route_id,
            "method_contract": method_contract,
            "calculation_hashes": calculation_hashes,
            "endpoint_valid": bool((verification.get("endpoint_validation") or {}).get("valid")),
            "imaginary_mode_count": (
                (verification.get("transition_state") or {}).get("frequency_check") or {}
            ).get("NumImaginary"),
            "irc_valid": bool(
                ((verification.get("irc") or {}).get("endpoint_validation") or {}).get("valid")
            ),
        }).split(":", 1)[1]
        path_claim["claim_id"] = path_claim_id
        rate_claim = {
            "schema_version": 1,
            "kind": "scientific_claim",
            "claim_type": "solution_eyring_rate",
            "status": "estimated_with_explicit_bounds",
            "reaction_id": candidate["reaction_id"],
            "path_claim_id": path_claim_id,
            "claim": "The reaction has a condition-specific conventional solution-Eyring first-order rate using the observed ORCA free-energy recipe within the declared model.",
            "thermochemistry_recipe": thermochemistry_recipe,
            "conformer_assessment": conformer_assessment,
            "temperature_K": environment.temperature_K,
            "rate_constant_s": rate_s,
            "activation_free_energy_hartree": barrier_hartree,
            "activation_free_energy_kcal_mol": barrier_hartree * 627.5094740631,
            "transmission_coefficient": 1.0,
            "standard_state": {
                "name": "one_molar_solution",
                "correction_kcal_mol": 0.0,
                "reason": "The unimolecular reactant and transition state each have molecularity one, so the concentration standard-state correction cancels.",
            },
            "thermochemistry": {
                "reactant": reactant_thermo,
                "transition_state": ts_thermo,
            },
            "uncertainty": {
                **uncertainty,
                "lower_rate_multiplier": lower_multiplier,
                "upper_rate_multiplier": upper_multiplier,
                "interpretation": "Deterministic multiplicative bounds derived from the declared composite free-energy-barrier bound; not a confidence interval.",
            },
            "limitations": [
                "Conventional TST with transmission coefficient one and the separately recorded observed ORCA entropy treatment.",
                "Quasi-RRHO entropy is not an explicit hindered-rotor or anharmonic treatment.",
                "No tunneling, variational, explicit anharmonic, pressure-dependent, diffusion, or explicit-solvent correction.",
                "The composite barrier bound must be calibrated externally for the intended chemistry class.",
                "Thermal corrections match the requested temperature, but named implicit-solvent parameter temperature dependence is not modeled.",
            ],
        }
        rate_claim_id = "claim:rate:" + _canonical_hash({
            "claim_type": rate_claim["claim_type"],
            "reaction_id": candidate["reaction_id"],
            "path_claim_id": path_claim_id,
            "temperature_K": environment.temperature_K,
            "rate_constant_s": rate_s,
            "activation_free_energy_hartree": barrier_hartree,
            "uncertainty": rate_claim["uncertainty"],
            "standard_state": rate_claim["standard_state"],
            "thermochemistry_recipe_id": thermochemistry_recipe["recipe_id"],
        }).split(":", 1)[1]
        rate_claim["claim_id"] = rate_claim_id

        verified_reaction = VerifiedSolutionReaction(
            reaction_id=candidate["reaction_id"], reactant_id=candidate["reactant_id"],
            products=candidate["products"], rate_constant_s=rate_s,
            temperature_K=environment.temperature_K, solvent=environment.solvent,
            solvent_model=environment.solvent_model, standard_state_mol_L=1.0,
            verification_status="verified_solution_tst", rate_theory="solution_eyring_tst",
            source="STORCA solvent-consistent ORCA endpoint/NEB-TS/frequency/IRC evidence",
            rate_claim_id=rate_claim_id, path_claim_id=path_claim_id,
            standard_state_correction="one_molar_solution",
            calculation_hashes=tuple(calculation_hashes),
            lower_rate_multiplier=lower_multiplier, upper_rate_multiplier=upper_multiplier,
            uncertainty_source=str(uncertainty["source"]),
            uncertainty_components=tuple(uncertainty["components"]),
            applicability_flags=candidate["applicability_flags"],
        )
        final_network_dict = {
            **network,
            "reactions": [item.as_dict() for item in context["verified_existing"]]
            + [verified_reaction.as_dict()],
        }
        final_network = VerifiedSolutionNetwork.from_dict(final_network_dict)
        coverage_complete = (
            final_network.coverage_status
            == "verified_complete_within_declared_unimolecular_scope"
            and bool(final_network.coverage_claim_id)
            and bool(final_network.discovery_sources)
            and not final_network.unresolved_candidate_ids
        )
        if not coverage_complete:
            base_result.update(
                status="verified_rate_evidence_emitted_network_unresolved",
                path_claim=path_claim,
                rate_claim=rate_claim,
                verified_reaction=verified_reaction.as_dict(),
                orca_input_audit=input_audit,
                reason=(
                    "The solvent-consistent path and bounded solution-Eyring rate were retained, "
                    "but no stability bundle was emitted because reaction-network coverage remains unresolved."
                ),
            )
            _write_json(result_path, base_result)
            return base_result
        bundle = {
            "schema_version": 1,
            "kind": "verified_implicit_solution_stability_bundle",
            "target_smiles": request["target_smiles"],
            "environment": _environment_bundle_dict(environment),
            "network": final_network.as_dict(),
            "scientific_claims": [path_claim, rate_claim],
            "producer": {
                "kind": "storca_implicit_solution_path_verifier",
                "schema_version": 1,
                "request_sha256": base_result["request_sha256"],
                "method_contract": method_contract,
            },
        }
        _write_json(bundle_path, bundle)
        propagation = run_implicit_solution_stability_bundle(
            bundle, workdir / "propagation",
        )
        base_result.update(
            status="verified_bundle_emitted",
            path_claim=path_claim,
            rate_claim=rate_claim,
            verified_reaction=verified_reaction.as_dict(),
            orca_input_audit=input_audit,
            verified_bundle=str(bundle_path),
            propagation_result=propagation,
            reason=(
                "A solvent-consistent verified path and bounded solution-Eyring rate were retained; "
                "the network was repropagated through the existing solution lane."
            ),
        )
        _write_json(result_path, base_result)
        return base_result
    except Exception as error:
        result = {
            "schema_version": 1,
            "kind": "implicit_solution_path_verification",
            "status": "invalid_or_incomplete_evidence",
            "reason": f"{type(error).__name__}: {error}",
            "verified_bundle": None,
            "result_json": str(result_path),
        }
        if result_path.is_file():
            result["result_json"] = None
            result["retained_result_json"] = str(result_path)
        else:
            _write_json(result_path, result)
        return result


def load_and_run_implicit_solution_path_verification(
    request_path: Path,
    workdir: Path,
    **kwargs: Any,
) -> dict[str, Any]:
    request = json.loads(Path(request_path).read_text())
    return run_implicit_solution_path_verification(request, workdir, **kwargs)


__all__ = [
    "REQUIRED_UNCERTAINTY_COMPONENTS",
    "load_and_run_implicit_solution_path_verification",
    "run_implicit_solution_path_verification",
    "solution_eyring_coefficient",
    "solution_orca_keywords",
]
