"""Bounded reactant sampling/refinement in the existing solution rate producer.

UFF energies only rank geometric seeds. Only condition/recipe-matched ORCA
minima enter the conditional free-energy calculation. A single verified TS
does not establish fast conformer equilibration or complete TS coverage.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from .orca_thermochemistry_recipe import observed_orca_thermochemistry_recipe
from .reaction_path import (
    _active_surface_budget, _optimize_endpoint, _surface_context, _validate_retained_orca_job_input,
    validated_orca_temperature_thermochemistry,
)
from .route_geometry import _build_endpoint_side, read_xyz, validate_declared_connectivity, write_xyz
from .route_thermochemistry import HARTREE_TO_J_MOL, GAS_CONSTANT_J_MOL_K
from .surface_resolution import SurfaceResolutionPolicy


CONFORMER_RATE_CONTRACT = "bounded_rdkit_orca_reactant_conformers_v1"


@dataclass(frozen=True)
class ConformerRatePolicy:
    policy_id: str = CONFORMER_RATE_CONTRACT
    random_seeds: tuple[int, ...] = (42, 1729)
    seeds_per_random_seed: int = 12
    maximum_refinements: int = 4
    maximum_uff_iterations: int = 200
    maximum_embed_iterations: int = 200
    maximum_symmetry_matches: int = 1024
    maximum_heavy_atoms: int = 20
    seed_rmsd_angstrom: float = 0.25
    refined_rmsd_angstrom: float = 0.15
    duplicate_free_energy_tolerance_hartree: float = 1e-4
    maximum_seed_wall_seconds: float = 30.0

    def __post_init__(self):
        if self.policy_id != CONFORMER_RATE_CONTRACT:
            raise ValueError("Unknown conformer policy contract")
        for name, maximum in (("seeds_per_random_seed", 128), ("maximum_refinements", 16),
                              ("maximum_uff_iterations", 2000), ("maximum_embed_iterations", 2000),
                              ("maximum_symmetry_matches", 4096), ("maximum_heavy_atoms", 40)):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
                raise ValueError(f"Invalid bounded conformer policy: {name}")
        if not self.random_seeds or len(self.random_seeds) > 4 or len(set(self.random_seeds)) != len(self.random_seeds):
            raise ValueError("Conformer policy requires one to four distinct random seeds")
        if any(isinstance(value, bool) or not isinstance(value, int) or not 0 <= value < 2**31 for value in self.random_seeds):
            raise ValueError("Invalid conformer random seed")
        for name in ("seed_rmsd_angstrom", "refined_rmsd_angstrom", "duplicate_free_energy_tolerance_hartree", "maximum_seed_wall_seconds"):
            if isinstance(getattr(self, name), bool) or not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"Invalid conformer policy: {name}")


def _hash(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _identity(value) -> str:
    return hashlib.sha256(json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def _immutable_json(path: Path, value):
    text = json.dumps(_jsonable(value), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if path.is_file():
        if path.read_text() != text:
            raise ValueError("Conformer workspace contains a different immutable contract")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        stream.write(text)


def _validate_refinement_sources(path: Path, contract_sha256: str):
    """Completed refinement evidence is immutable, not fresh cache credit."""
    record = json.loads(path.read_text())
    body = {key: value for key, value in record.items() if key != "evidence_sha256"}
    if record.get("evidence_sha256") != _identity(body) or record.get("contract_sha256") != contract_sha256:
        raise ValueError("Retained conformer refinement evidence contract/hash changed")
    for source in record["sources"]:
        if _hash(Path(source["path"])) != source["sha256"]:
            raise ValueError("Retained conformer refinement source changed")


def _molecule(smiles, policy):
    from rdkit import Chem
    molecule = _build_endpoint_side((smiles,)).explicit_molecules[0]
    if (len(Chem.GetMolFrags(molecule)) != 1 or molecule.GetNumHeavyAtoms() > policy.maximum_heavy_atoms
            or any(atom.GetFormalCharge() or atom.GetNumRadicalElectrons() or atom.GetIsotope() for atom in molecule.GetAtoms())):
        raise ValueError("Conformer producer supports one small neutral closed-shell molecule with verified ordinary isotope treatment only")
    return molecule


def validate_conformer_production_scope(route, *, policy: ConformerRatePolicy | None = None):
    """Deterministic eligibility, shared by preparation and rate production."""
    policy = policy or ConformerRatePolicy()
    if len(route.reactant_smiles) != 1 or route.charge != 0 or route.multiplicity != 1:
        raise ValueError("Conformer producer requires one neutral singlet reactant")
    return _molecule(route.reactant_smiles[0], policy)


def conformer_sampling_decision(smiles: str) -> dict:
    """A trigger, not a proof that untriggered systems have one conformer."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid reactant identity for conformer diagnosis")
    rotors = int(rdMolDescriptors.CalcNumRotatableBonds(molecule))
    flexible_rings = [list(ring) for ring in molecule.GetRingInfo().AtomRings()
                      if len(ring) >= 4 and not all(molecule.GetAtomWithIdx(index).GetIsAromatic() for index in ring)]
    return {
        "policy_id": CONFORMER_RATE_CONTRACT,
        "sampling_required": bool(rotors or flexible_rings),
        "rotatable_bond_count": rotors, "nonaromatic_ring_count_at_least_four": len(flexible_rings),
        "limitations": ["Rotor/ring triggers are not a completeness diagnostic for OH/NH rotation, inversion, tautomers, or other missed conformers."],
    }


def _best_rmsd(molecule, left, right, policy):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    matches = molecule.GetSubstructMatches(molecule, uniquify=False, useChirality=True,
                                          maxMatches=policy.maximum_symmetry_matches + 1)
    if len(matches) > policy.maximum_symmetry_matches:
        raise ValueError("Conformer symmetry-match budget exhausted; population deduplication unresolved")
    if not matches:
        raise ValueError("Conformer identity has no chirality-preserving atom correspondence")
    copies = []
    for coordinates in (left, right):
        copy = Chem.Mol(molecule)
        copy.RemoveAllConformers()
        conformer = Chem.Conformer(copy.GetNumAtoms())
        for index, xyz in enumerate(coordinates):
            conformer.SetAtomPosition(index, tuple(float(value) for value in xyz))
        copy.AddConformer(conformer)
        copies.append(copy)
    return float(rdMolAlign.GetBestRMS(
        copies[0], copies[1], map=[list(enumerate(match)) for match in matches],
        maxMatches=policy.maximum_symmetry_matches, symmetrizeConjugatedTerminalGroups=False, numThreads=1,
    ))


def acquire_reactant_conformer_seeds(smiles: str, workdir: Path, *, policy: ConformerRatePolicy | None = None) -> dict:
    """ETKDGv3/UFF seeds, immutable artifacts, no populations or absence claim."""
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    policy = policy or ConformerRatePolicy()
    molecule = _molecule(smiles, policy)
    workdir = Path(workdir).resolve()
    contract = {
        "kind": "storca_reactant_conformer_seed_contract", "schema_version": 1,
        "ordered_input_smiles": smiles, "rdkit_version": rdBase.rdkitVersion,
        "policy": asdict(policy), "backend": "RDKit_ETKDGv3_UFF", "threads": 1,
        "energy_role": "seed_ranking_only", "population_model": "not_assigned",
    }
    contract_path, manifest_path = workdir / "sampling-contract.json", workdir / "sampling.json"
    _immutable_json(contract_path, contract)
    if manifest_path.is_file():
        record = json.loads(manifest_path.read_text())
        if record.get("contract_sha256") != _hash(contract_path):
            raise ValueError("Conformer sampling manifest contract hash changed")
        if record.get("sampling_evidence_sha256") != _identity({key: value for key, value in record.items() if key != "sampling_evidence_sha256"}):
            raise ValueError("Conformer sampling manifest evidence hash changed")
        for seed in record["selected"]:
            path = Path(seed["xyz"])
            if path.resolve().parent != workdir or not path.is_file() or _hash(path) != seed["xyz_sha256"]:
                raise ValueError("Conformer seed hash changed")
        return record
    started = time.monotonic()
    attempts, embedding_attempts, candidates, generated = [], [], [], 0
    time_exhausted = False
    for seed_number in policy.random_seeds:
        if time.monotonic() - started >= policy.maximum_seed_wall_seconds:
            time_exhausted = True
            break
        trial = Chem.Mol(molecule)
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = seed_number
        parameters.numThreads = 1
        parameters.maxIterations = policy.maximum_embed_iterations
        parameters.pruneRmsThresh = -1.0
        try:
            identifiers = list(AllChem.EmbedMultipleConfs(trial, numConfs=policy.seeds_per_random_seed, params=parameters))
            embedding_attempts.append({
                "random_seed": seed_number, "requested_count": policy.seeds_per_random_seed,
                "embedded_count": len(identifiers),
                "status": "completed" if len(identifiers) == policy.seeds_per_random_seed else "embedding_incomplete",
            })
        except Exception as error:
            embedding_attempts.append({
                "random_seed": seed_number, "requested_count": policy.seeds_per_random_seed,
                "embedded_count": 0, "status": "embedding_failed", "reason": f"{type(error).__name__}: {error}",
            })
            continue
        generated += len(identifiers)
        for identifier in identifiers:
            if time.monotonic() - started >= policy.maximum_seed_wall_seconds:
                time_exhausted = True
                break
            record = {"random_seed": seed_number, "conformer_id": int(identifier)}
            try:
                if not AllChem.UFFHasAllMoleculeParams(trial):
                    raise ValueError("UFF lacks parameters; no force-field fallback is assumed")
                converged = AllChem.UFFOptimizeMolecule(trial, confId=identifier, maxIters=policy.maximum_uff_iterations) == 0
                energy = float(AllChem.UFFGetMoleculeForceField(trial, confId=identifier).CalcEnergy())
                coordinates = np.asarray(trial.GetConformer(identifier).GetPositions())
                if not converged or not math.isfinite(energy) or not np.isfinite(coordinates).all():
                    raise ValueError("UFF seed relaxation incomplete/nonfinite")
                record.update(status="seed_unverified", approximate_energy_kcal_mol=energy)
                candidates.append((energy, seed_number, int(identifier), coordinates.copy()))
            except Exception as error:
                record.update(status="seed_failed", reason=f"{type(error).__name__}: {error}")
            attempts.append(record)
        if time_exhausted:
            break
    unique = []
    for item in sorted(candidates, key=lambda item: item[:3]):
        if time.monotonic() - started >= policy.maximum_seed_wall_seconds:
            time_exhausted = True
            break
        if all(_best_rmsd(molecule, item[3], other[3], policy) > policy.seed_rmsd_angstrom for other in unique):
            unique.append(item)
    selected = []
    elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    for index, (energy, seed_number, identifier, coordinates) in enumerate(unique[:policy.maximum_refinements], 1):
        path = workdir / f"seed-{index:03d}.xyz"
        if path.exists():
            raise ValueError("Partial conformer acquisition has a retained seed without a manifest; use a new workspace")
        write_xyz(path, elements, coordinates, f"RDKit ETKDGv3/UFF seed {seed_number}:{identifier}; not a population or quantum minimum")
        connectivity = validate_declared_connectivity(path, (smiles,))
        if not connectivity.get("valid"):
            raise ValueError("Generated conformer changed declared identity/stereochemistry")
        selected.append({
            "selected_position": index, "xyz": str(path), "xyz_sha256": _hash(path),
            "approximate_energy_kcal_mol": energy, "approximate_energy_role": "seed_ranking_only",
            "random_seed": seed_number, "conformer_id": identifier, "connectivity": connectivity,
        })
    record = {
        "kind": "storca_reactant_conformer_seed_acquisition", "schema_version": 1,
        "contract_sha256": _hash(contract_path), "generated_count": generated,
        "unique_seed_count": len(unique), "selected": selected, "attempts": attempts,
        "embedding_attempts": embedding_attempts,
        "known_seed_failures": any(item["status"] != "completed" for item in embedding_attempts)
            or any(item["status"] == "seed_failed" for item in attempts),
        "selection_truncated": len(unique) > len(selected), "seed_wall_budget_exhausted": time_exhausted,
        "status": "bounded_seeds_acquired" if selected else "no_valid_seed",
        "coverage_complete": False, "population_model": "not_assigned",
        "measured_seed_wall_seconds": time.monotonic() - started,
        "limitations": [
            "UFF ranks seeds only; these are not solution free energies, rates, or populations.",
            "RMSD seed clustering and finite sampling have uncalibrated false-negative risk.",
            "Cheap gas-like seed generation does not resolve solvent-dependent conformer coverage.",
        ],
    }
    record["sampling_evidence_sha256"] = _identity(record)
    _immutable_json(manifest_path, record)
    return record


def conditional_single_path_equilibrium_scenario(free_energies_hartree, selected_free_energy_hartree,
                                                ts_free_energy_hartree, temperature_K) -> dict:
    """Restricted flux through one verified TS, conditional on rapid equilibrium.

    This is NOT an aggregate rate admission. Other TS conformers/channels and
    interconversion may change it. Molecular degeneracies are not inferred.
    """
    original_values = list(free_energies_hartree)
    if any(isinstance(value, bool) for value in [*original_values, selected_free_energy_hartree, ts_free_energy_hartree, temperature_K]):
        raise ValueError("Conditional conformer scenario requires numeric, not boolean, evidence")
    values = [float(value) for value in original_values]
    selected, transition_state, temperature = map(float, (selected_free_energy_hartree, ts_free_energy_hartree, temperature_K))
    if (not values or not all(math.isfinite(value) for value in [*values, selected, transition_state, temperature])
            or temperature <= 0 or selected not in values or transition_state <= selected):
        raise ValueError("Conditional conformer scenario requires finite matched free energies, selected state and positive barriers")
    minimum = min(values)
    beta = HARTREE_TO_J_MOL / (GAS_CONSTANT_J_MOL_K * temperature)
    factors = [math.exp(-(value - minimum) * beta) for value in values]
    partition = math.fsum(factors)
    fraction = math.exp(-(selected - minimum) * beta) / partition
    log_rate = math.log(1.380649e-23 * temperature / 6.62607015e-34) - (transition_state - minimum) * beta - math.log(partition)
    if not -745 < log_rate < 709:
        raise ValueError("Conditional conformer rate is numerically unresolved")
    return {
        "status": "conditional_scenario_only", "network_rate_admission_allowed": False,
        "rate_constant_s": math.exp(log_rate), "temperature_K": temperature,
        "selected_path_reactant_fraction": fraction,
        "assumptions": ["Rapid equilibrium among the explicitly refined reactant minima is assumed, not demonstrated.",
                        "Only the selected verified TS is included; TS-conformer coverage is unresolved.",
                        "Each represented geometry has unit additional degeneracy; symmetry already present in its ORCA thermochemistry is not added twice."],
        "interpretation": "Conditional flux scenario through the selected path, not a confidence bound or total reaction rate.",
    }


def run_reactant_conformer_rate_diagnostic(route, path_result, workdir: Path, *, reactant_thermo: dict,
                                          ts_thermo: dict, ncores: int, method_keywords: list[str],
                                          timeout_seconds: float | None,
                                          policy: ConformerRatePolicy | None = None) -> dict:
    """Refine bounded seeds using the existing endpoint executor and route budget."""
    policy = policy or ConformerRatePolicy()
    molecule = validate_conformer_production_scope(route, policy=policy)
    decision = conformer_sampling_decision(route.reactant_smiles[0])
    if not decision["sampling_required"]:
        return {"status": "sampling_not_triggered", "network_rate_admission_allowed": True,
                "decision": decision, "coverage_complete": False,
                "reason": "Legacy single-path estimate remains conditional on its separately supplied conformer/model bound; this trigger does not prove conformer completeness."}
    if isinstance(ncores, bool) or not isinstance(ncores, int) or not 1 <= ncores <= 4:
        raise ValueError("Conformer refinement requires one to four aggregate sequential cores")
    verification = path_result["route_verification"]
    reactant_endpoint = (verification.get("endpoint_validation") or {}).get("reactant") or {}
    xyz = reactant_endpoint.get("endpoint_xyz")
    if not xyz:
        records = (verification.get("stationary_points") or {}).get("reactants") or []
        xyz = (records[0].get("optimized_xyz") or records[0].get("xyz")) if len(records) == 1 else None
    if not xyz or not Path(xyz).is_file():
        raise ValueError("Conformer refinement lacks the selected validated reactant geometry")
    selected_xyz = Path(xyz)
    identity = validate_declared_connectivity(selected_xyz, route.reactant_smiles)
    if not identity.get("valid"):
        raise ValueError("Selected reactant conformer identity/stereochemistry is unresolved")
    temperature = reactant_thermo["requested_temperature_K"]
    if ts_thermo["requested_temperature_K"] != temperature:
        raise ValueError("Conformer production cannot mix thermochemistry temperatures")
    recipe = reactant_thermo["thermochemistry_recipe"]
    if recipe["recipe_id"] != ts_thermo["thermochemistry_recipe"]["recipe_id"]:
        raise ValueError("Conformer production cannot mix thermochemistry recipes")
    context = _surface_context(Path(workdir))
    electronic = (context or {}).get("electronic_structure") or {}
    if (not context or _identity(context["route"]) != _identity(route.as_dict())
            or electronic.get("method_keywords") != method_keywords
            or electronic.get("charge") != route.charge or electronic.get("multiplicity") != route.multiplicity
            or electronic.get("thermochemistry_temperature_K") != temperature
            or (context.get("environment_applicability") or {}).get("canonical_model") != "implicit_solution"):
        raise ValueError("Conformer refinement requires the matching source-bound solution route/method/temperature context")
    contract = {
        "kind": "storca_reactant_conformer_refinement_contract", "schema_version": 1,
        "policy": asdict(policy), "route": route.as_dict(), "method_keywords": method_keywords,
        "ncores": ncores, "temperature_K": temperature, "thermochemistry_recipe_id": recipe["recipe_id"],
        "selected_geometry": {"path": str(selected_xyz), "sha256": _hash(selected_xyz)},
        "selected_thermochemistry": recipe["source_output"], "selected_input": recipe["source_input"],
        "ts_thermochemistry": ts_thermo["thermochemistry_recipe"]["source_output"],
        "ts_input": ts_thermo["thermochemistry_recipe"]["source_input"],
        "scientific_context_sha256": context["scientific_context_sha256"],
        "population_model": "unresolved_without_interconversion_evidence",
    }
    workdir = Path(workdir) / _identity(contract)
    _immutable_json(workdir / "refinement-contract.json", contract)
    active_budget = _active_surface_budget(workdir)
    if active_budget is None:
        raise ValueError("Conformer refinements require the existing shared route budget; no independent budget is created")
    inherited_values = {key: value for key, value in active_budget[1]["policy"].items() if key not in {"kind", "scientific_hash"}}
    inherited_policy = SurfaceResolutionPolicy(**inherited_values)
    seeds = acquire_reactant_conformer_seeds(route.reactant_smiles[0], workdir / "seeds", policy=policy)
    baseline = {
        "role": "selected_verified_path_reactant", "xyz": str(selected_xyz), "xyz_sha256": _hash(selected_xyz),
        "gibbs_free_energy_hartree": reactant_thermo["gibbs_free_energy_hartree"],
        "thermochemistry": reactant_thermo, "cluster_representative": 0,
    }
    minima, records, unresolved = [baseline], [], []
    selected_elements, selected_coordinates, _ = read_xyz(selected_xyz)
    coordinates = [selected_coordinates]
    for seed in seeds["selected"]:
        record = {"seed": seed}
        endpoint_directory = workdir / f"minimum-{seed['selected_position']:03d}"
        retained_sources = endpoint_directory / "validated-sources.json"
        try:
            contract_sha256 = _hash(workdir / "refinement-contract.json")
            if retained_sources.is_file():
                _validate_refinement_sources(retained_sources, contract_sha256)
            # The native endpoint executor checks the inherited budget only
            # when dispatching. Completed outputs remain revalidatable after
            # exhaustion; calling the gate here would wrongly reject reuse.
            result = _optimize_endpoint(Path(seed["xyz"]), endpoint_directory,
                                        charge=route.charge, multiplicity=route.multiplicity, ncores=ncores,
                                        method_keywords=method_keywords, timeout_seconds=timeout_seconds,
                                        surface_resolution_policy=inherited_policy)
            record["endpoint_workflow"] = _jsonable(result)
            if not (result.get("frequency_check") or {}).get("IsMinimum"):
                raise ValueError("Conformer refinement is not a frequency-validated minimum")
            output = Path(result["frequency"]["out"])
            optimization_output = Path(result["optimization"]["out"])
            for executed_input in (output.with_suffix(".inp"), optimization_output.with_suffix(".inp")):
                if not executed_input.with_suffix(".contract.json").is_file():
                    raise ValueError("New conformer refinement requires exact executed native ORCA job contracts; missing ownership is not legacy credit")
                _validate_retained_orca_job_input(executed_input)
            refined_xyz = Path(result["optimized_xyz"])
            elements, positions, _ = read_xyz(refined_xyz)
            if elements != selected_elements or not validate_declared_connectivity(refined_xyz, route.reactant_smiles).get("valid"):
                raise ValueError("Refined conformer changed composition, identity or declared stereochemistry")
            frequencies = (result.get("frequency_check") or {}).get("AllFrequencies") or []
            if len(frequencies) != 3 * len(elements) or not all(math.isfinite(value) for value in frequencies):
                raise ValueError("Refined conformer requires its complete finite 3N ORCA frequency list, including external modes")
            thermo = validated_orca_temperature_thermochemistry(output, temperature)
            refined_recipe = observed_orca_thermochemistry_recipe(output)
            if (refined_recipe["recipe_id"] != recipe["recipe_id"]
                    or refined_recipe["source_output"]["sha256"] != thermo["source_output_sha256"]
                    or refined_recipe["source_input"]["sha256"] != thermo["temperature_input_evidence"]["source_input_sha256"]):
                raise ValueError("Refined conformer thermochemistry recipe/source differs")
            distances = [_best_rmsd(molecule, positions, previous, policy) for previous in coordinates]
            duplicate = next((index for index, distance in enumerate(distances) if distance <= policy.refined_rmsd_angstrom), None)
            entry = {"xyz": str(refined_xyz), "xyz_sha256": _hash(refined_xyz),
                     "gibbs_free_energy_hartree": thermo["gibbs_free_energy_hartree"],
                     "thermochemistry": {**thermo, "thermochemistry_recipe": refined_recipe}}
            if duplicate is not None:
                if abs(entry["gibbs_free_energy_hartree"] - minima[duplicate]["gibbs_free_energy_hartree"]) > policy.duplicate_free_energy_tolerance_hartree:
                    raise ValueError("Near-duplicate conformer geometry has inconsistent thermochemistry")
                record.update(status="refined_duplicate_minimum", cluster_representative=duplicate, evidence=entry)
            else:
                entry["cluster_representative"] = len(minima)
                minima.append(entry)
                coordinates.append(positions)
                record.update(status="distinct_refined_minimum", evidence=entry)
            source_paths = {refined_xyz, output, output.with_suffix(".inp"),
                            optimization_output, optimization_output.with_suffix(".inp")}
            for source_input in (optimization_output.with_suffix(".inp"), output.with_suffix(".inp")):
                if source_input.with_suffix(".contract.json").is_file():
                    source_paths.add(source_input.with_suffix(".contract.json"))
            source_record = {
                "kind": "storca_validated_conformer_refinement_sources", "schema_version": 1,
                "contract_sha256": contract_sha256,
                "sources": [{"path": str(path), "sha256": _hash(path)} for path in sorted(source_paths)],
            }
            source_record["evidence_sha256"] = _identity(source_record)
            _immutable_json(retained_sources, source_record)
            record["validated_sources"] = {"path": str(retained_sources), "sha256": _hash(retained_sources)}
        except Exception as error:
            record.update(status="conformer_refinement_unresolved", reason=f"{type(error).__name__}: {error}")
            unresolved.append(seed["selected_position"])
        records.append(record)
    for record in records:
        entry = record.get("evidence")
        if not entry:
            continue
        retained_recipe = entry["thermochemistry"]["thermochemistry_recipe"]
        sources = [{"path": entry["xyz"], "sha256": entry["xyz_sha256"]},
                   retained_recipe["source_input"], retained_recipe["source_output"]]
        if record.get("validated_sources"):
            sources.append(record["validated_sources"])
        if any(_hash(Path(source["path"])) != source["sha256"] for source in sources):
            raise ValueError("Refined conformer source changed before rate assessment")
    for source in (contract["selected_geometry"], recipe["source_input"], recipe["source_output"], contract["ts_thermochemistry"], contract["ts_input"]):
        if _hash(Path(source["path"])) != source["sha256"]:
            raise ValueError("Selected path evidence changed during conformer refinement")
    allowed = bool(seeds["selected"] and not unresolved and not seeds["selection_truncated"] and not seeds["known_seed_failures"]
                   and not seeds["seed_wall_budget_exhausted"] and len(minima) == 1)
    status = ("single_refined_conformer_cluster" if allowed else
              "multiple_refined_minima_interconversion_unresolved" if len(minima) > 1 else "conformer_sampling_or_refinement_unresolved")
    assessment = {
        "kind": "storca_reactant_conformer_rate_assessment", "schema_version": 1,
        "status": status, "network_rate_admission_allowed": allowed, "decision": decision,
        "contract": str(workdir / "refinement-contract.json"), "contract_sha256": _hash(workdir / "refinement-contract.json"),
        "seeds": seeds, "refinements": records, "minima": minima,
        "unresolved_refinement_positions": unresolved, "coverage_complete": False,
        "interconversion_evidence": "unavailable", "ts_conformer_coverage": "one_selected_verified_path",
        "reason": ("Bounded refinements returned one geometric minimum cluster; the original single-path estimate still needs its supplied sampling/model bound." if allowed else
                   "Distinct minima, omitted seeds or failed refinements prevent a defensible aggregate rate without interconversion and TS-conformer evidence."),
        "limitations": ["This finite sampler does not prove conformer/network completeness or calibrate unsampled risk.",
                        "Geometric clustering thresholds are protocol settings, not physical error bounds.",
                        "No degeneracy is inferred from the number of generated/duplicate seeds."],
    }
    if len(minima) > 1:
        try:
            assessment["conditional_single_path_equilibrium_scenario"] = conditional_single_path_equilibrium_scenario(
                [entry["gibbs_free_energy_hartree"] for entry in minima], baseline["gibbs_free_energy_hartree"],
                ts_thermo["gibbs_free_energy_hartree"], temperature,
            )
        except ValueError as error:
            assessment["conditional_single_path_equilibrium_scenario"] = {
                "status": "numerically_unresolved", "network_rate_admission_allowed": False,
                "reason": str(error),
            }
    # Reports are content-addressed revisions. Resume revalidates every job;
    # old assessments are not overwritten or used as unconditional rate credit.
    _immutable_json(workdir / f"assessment-{_identity(assessment)}.json", assessment)
    return assessment
