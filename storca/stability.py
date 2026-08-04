"""Combined ORCA local-minimum and RMG kinetic stability evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
from copy import deepcopy

from src.inputgen import create_rmg_input
from src.orca_runner import run_rmg_supervised
from src.parser import parse_chemkin_annotated

from .conditions import ConditionSpec, build_condition_spec
from .reachability import assess_kinetic_relevance, enrich_candidate_route
from .rmg_evidence import parse_annotated_reaction_provenance, parse_final_solver_profile, parse_species_dictionary
from .rmg_execution import (
    assess_rmg_execution,
    merge_collision_rate_validation,
    parse_collision_rate_violators,
    requested_time_coverage,
)
from .runs import write_metadata
from .workflow import run_optimization_and_frequency


STABILITY_SCENARIOS = {
    "low-pressure-intrinsic-gas-screen": {
        "phase": "homogeneous gas-phase low-pressure surrogate",
        "atmosphere": "dilute target in nitrogen at low total pressure",
        "model_applicability": (
            "A dilute target in a low-pressure nitrogen bath used to approach isolated-molecule "
            "thermal behavior while retaining a defined collision partner. It is not a literal vacuum, condensed phase, "
            "surface, container, or excited-state model."
        ),
        "additional_species": [{"label": "nitrogen", "smiles": "N#N", "reactive": False}],
        "initial_mole_fractions": {"stability": 0.01, "nitrogen": 0.99},
    },
    "dry-inert-gas-screen": {
        "phase": "homogeneous gas-phase surrogate",
        "atmosphere": "dry nitrogen",
        "model_applicability": "A dilute target molecule in an inert nitrogen bath; not a condensed-phase, humidity, light-exposure, or container-compatibility model.",
        "additional_species": [{"label": "nitrogen", "smiles": "N#N", "reactive": False}],
        "initial_mole_fractions": {"stability": 0.01, "nitrogen": 0.99},
    },
    "ambient-air-gas-screen": {
        "phase": "homogeneous gas-phase surrogate",
        "atmosphere": "dry air (oxygen/nitrogen)",
        "model_applicability": "A dilute target molecule in dry air. This does not model liquid, solid, dissolved, humid, illuminated, or surface-contact storage.",
        "additional_species": [
            {"label": "oxygen", "smiles": "[O][O]", "reactive": True},
            {"label": "nitrogen", "smiles": "N#N", "reactive": False},
        ],
        "initial_mole_fractions": {"stability": 0.01, "oxygen": 0.2079, "nitrogen": 0.7821},
    },
    "humid-air-gas-screen": {
        "phase": "homogeneous gas-phase surrogate",
        "atmosphere": "humid air (water vapour/oxygen/nitrogen)",
        "model_applicability": "A dilute target in homogeneous humid gas. This does not model droplets, dissolved ions, surfaces, or container compatibility.",
        "additional_species": [
            {"label": "water", "smiles": "O", "reactive": True},
            {"label": "oxygen", "smiles": "[O][O]", "reactive": True},
            {"label": "nitrogen", "smiles": "N#N", "reactive": False},
        ],
        # This named standalone scenario is deliberately dry by default.  The
        # ladder supplies a humidity-derived composition for humid air.
        "initial_mole_fractions": {"stability": 0.01, "water": 0.0, "oxygen": 0.2079, "nitrogen": 0.7821},
    },
}

RMG_RESOURCE_TIERS = {
    "quick-screen": {"walltime": "00:00:10:00", "max_processes": 1, "max_iterations": 100, "max_edge_species": 250, "max_objects_per_iteration": 1},
    "review-screen": {"walltime": "00:00:30:00", "max_processes": 1, "max_iterations": 300, "max_edge_species": 1000, "max_objects_per_iteration": 1},
}


def _reaction_space_budget(smiles: str, scenario: dict, *, intrinsic_scope: bool) -> dict:
    """Return a small, explicit RMG generation budget for this condition.

    The limit is based only on the supplied target and declared *reactive*
    environment species.  It therefore preserves a one-collision association
    product (target + O2, target + H2O, or target + target), while preventing
    an open-ended chain of larger oligomers from consuming the screen.
    """
    try:
        from rdkit import Chem

        target = Chem.MolFromSmiles(smiles)
        target_heavy = target.GetNumHeavyAtoms() if target is not None else 0
        target_radicals = sum(atom.GetNumRadicalElectrons() for atom in target.GetAtoms()) if target else 0
        partner_heavy = 0
        partner_radicals = 0
        partner_labels: list[str] = []
        for raw in scenario.get("additional_species") or []:
            if not raw.get("reactive", True):
                continue
            partner = Chem.MolFromSmiles(str(raw.get("smiles") or ""))
            if partner is None:
                continue
            heavy = partner.GetNumHeavyAtoms()
            radicals = sum(atom.GetNumRadicalElectrons() for atom in partner.GetAtoms())
            if heavy >= partner_heavy:
                partner_heavy = heavy
                partner_radicals = radicals
            partner_labels.append(str(raw.get("label")))
        if target_heavy <= 0:
            raise ValueError("target has no heavy atoms")
    except Exception:
        # A missing cheminformatics dependency must not quietly impose a
        # possibly wrong structural cap.  The retained evidence says so.
        return {
            "maximum_heavy_atoms": None,
            "maximum_radical_electrons": None,
            "filter_reactions": True,
            "scope": "unbounded_due_to_structure_parse_failure",
            "reactive_partner_labels": [],
        }

    # A target-target collision is still physically possible at the stated
    # concentration.  In reactive conditions, retain whichever direct
    # collision creates the larger one-step species.
    maximum_heavy = target_heavy if intrinsic_scope else target_heavy + max(target_heavy, partner_heavy)
    maximum_radicals = max(2, target_radicals + partner_radicals)
    return {
        "maximum_heavy_atoms": maximum_heavy,
        "maximum_radical_electrons": maximum_radicals,
        "filter_reactions": True,
        "scope": (
            "direct_intrinsic_target_decomposition" if intrinsic_scope
            else "one_declared_collision_product_network"
        ),
        "target_heavy_atoms": target_heavy,
        "largest_reactive_partner_heavy_atoms": partner_heavy,
        "reactive_partner_labels": partner_labels,
    }


def resolve_stability_configuration(
    scenario: str,
    screen_tier: str,
    *,
    rmg_walltime: str | None = None,
    rmg_max_processes: int | None = None,
    rmg_max_iterations: int | None = None,
    rmg_max_edge_species: int | None = None,
    rmg_max_objects_per_iteration: int | None = None,
) -> tuple[dict, dict]:
    """Resolve named scenario/tier defaults while retaining explicit overrides."""
    try:
        scenario_config = STABILITY_SCENARIOS[scenario]
    except KeyError as error:
        choices = ", ".join(sorted(STABILITY_SCENARIOS))
        raise ValueError(f"Unsupported stability scenario '{scenario}'. Available: {choices}") from error
    try:
        tier_defaults = RMG_RESOURCE_TIERS[screen_tier]
    except KeyError as error:
        choices = ", ".join(sorted(RMG_RESOURCE_TIERS))
        raise ValueError(f"Unsupported RMG screen tier '{screen_tier}'. Available: {choices}") from error
    resources = {
        "walltime": rmg_walltime or tier_defaults["walltime"],
        "max_processes": rmg_max_processes or tier_defaults["max_processes"],
        "max_iterations": rmg_max_iterations or tier_defaults["max_iterations"],
        "max_edge_species": rmg_max_edge_species or tier_defaults["max_edge_species"],
        "max_objects_per_iteration": rmg_max_objects_per_iteration or tier_defaults["max_objects_per_iteration"],
    }
    return {"name": scenario, **scenario_config}, {"name": screen_tier, **resources}


def _xyz_element_counts(path: Path) -> dict[str, int]:
    """Return element counts from an XYZ file without inferring bonding."""
    lines = Path(path).read_text().splitlines()
    try:
        atom_count = int(lines[0].strip())
        coordinates = lines[2:2 + atom_count]
    except (IndexError, ValueError):
        coordinates = lines
    counts: dict[str, int] = {}
    for line in coordinates:
        fields = line.split()
        if fields:
            counts[fields[0]] = counts.get(fields[0], 0) + 1
    if not counts:
        raise ValueError(f"No atomic coordinates found in {path}")
    return dict(sorted(counts.items()))


def collect_orca_evidence(
    smiles: str,
    run_dir: Path,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    ncores: int = 1,
    method_profile: dict | None = None,
) -> dict:
    """Run ORCA optimization/frequencies and retain local-minimum evidence.

    This is intentionally a structure check, not a kinetic-stability verdict.
    """
    run_dir = Path(run_dir)
    input_xyz = run_dir / "input.xyz"
    artifacts = {
        "input_xyz": str(input_xyz),
        "optimization_input": str(run_dir / "opt.inp"),
        "optimization_output": str(run_dir / "opt.out"),
        "optimized_xyz": str(run_dir / "opt.xyz"),
        "frequency_input": str(run_dir / "freq.inp"),
        "frequency_output": str(run_dir / "freq.out"),
    }
    evidence = {
        "kind": "orca_geometry_and_frequency_check",
        "status": "failed",
        "charge": charge,
        "multiplicity": multiplicity,
        "ncores": ncores,
        "method_profile": method_profile,
        "artifacts": artifacts,
        "interpretation": "An ORCA harmonic frequency calculation tests for a local minimum; it does not establish kinetic persistence.",
    }
    try:
        from src.molecule_tools import sanitize_smiles, smiles_to_xyz

        canonical_smiles = sanitize_smiles(smiles)
        smiles_to_xyz(canonical_smiles, input_xyz)
        input_elements = _xyz_element_counts(input_xyz)
        result = run_optimization_and_frequency(
            input_xyz,
            run_dir,
            charge=charge,
            multiplicity=multiplicity,
            ncores=ncores,
            method_keywords=(method_profile or {}).get("orca_keywords"),
        )
        optimized_elements = _xyz_element_counts(result["optimized_xyz"])
        frequency_check = result["frequency_check"]
        minimum_modes = max(1, 3 * sum(input_elements.values()) - 6)
        frequency_mode_coverage = {
            "observed_modes": len(frequency_check["AllFrequencies"]),
            "minimum_expected_for_nonlinear_structure": minimum_modes,
            "adequate": len(frequency_check["AllFrequencies"]) >= minimum_modes,
        }
        evidence.update(
            status="completed",
            canonical_smiles=canonical_smiles,
            local_minimum=frequency_check["IsMinimum"],
            frequency_check=frequency_check,
            input_element_counts=input_elements,
            optimized_element_counts=optimized_elements,
            element_inventory_retained=input_elements == optimized_elements,
            frequency_mode_coverage=frequency_mode_coverage,
        )
        if input_elements != optimized_elements:
            evidence["status"] = "incomplete"
            evidence["interpretation"] = "ORCA completed, but optimized XYZ does not retain the input element inventory. Review artifacts before using this evidence."
        elif not frequency_mode_coverage["adequate"]:
            evidence["status"] = "incomplete"
            evidence["interpretation"] = "ORCA output lacks enough harmonic modes for a minimum check; review artifacts before using this evidence."
    except Exception as error:
        evidence["failure_reason"] = str(error)
    return evidence


def _reuse_orca_evidence(
    evidence: dict,
    smiles: str,
    *,
    charge: int,
    multiplicity: int,
    method_profile: dict,
) -> dict:
    """Validate and annotate condition-independent parent ORCA evidence.

    The ladder changes the reactor composition and pressure, not the isolated
    parent electronic-structure calculation.  Reuse is therefore safe only
    when molecular identity, charge, spin, method, completed artifacts, and
    the local-minimum validation contract match exactly.
    """
    from src.molecule_tools import sanitize_smiles

    if not isinstance(evidence, dict) or evidence.get("status") != "completed":
        raise ValueError("Reusable ORCA evidence must be a completed evidence record")
    expected_smiles = sanitize_smiles(smiles)
    checks = {
        "canonical_smiles": evidence.get("canonical_smiles") == expected_smiles,
        "charge": evidence.get("charge") == charge,
        "multiplicity": evidence.get("multiplicity") == multiplicity,
        "method_profile": evidence.get("method_profile") == method_profile,
        "element_inventory_retained": evidence.get("element_inventory_retained") is True,
        "frequency_mode_coverage": (evidence.get("frequency_mode_coverage") or {}).get("adequate") is True,
    }
    artifacts = evidence.get("artifacts") or {}
    required_artifacts = ("optimized_xyz", "frequency_output")
    checks["artifacts_present"] = all(
        artifacts.get(name) and Path(artifacts[name]).is_file()
        for name in required_artifacts
    )
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError("Reusable ORCA evidence is incompatible: " + ", ".join(failed))
    reused = deepcopy(evidence)
    reused["reuse"] = {
        "status": "validated_reuse",
        "reason": "The ladder condition changed without changing the isolated parent electronic-state contract.",
        "source_artifacts": {name: artifacts[name] for name in required_artifacts},
        "validation_checks": checks,
    }
    return reused


def _collect_rmg_evidence_once(
    smiles: str,
    run_dir: Path,
    *,
    rmg_env: str | None = None,
    rmg_command: str = "rmg.py",
    barrier_threshold: float = 50.0,
    temperature: float = 298.0,
    pressure: float = 1.0,
    rmg_walltime: str | None = "00:00:10:00",
    rmg_max_processes: int | None = 1,
    rmg_max_iterations: int | None = 100,
    rmg_max_edge_species: int = 250,
    rmg_max_objects_per_iteration: int = 1,
    restart_from_seed: Path | None = None,
    rmg_hard_timeout_seconds: float | None = None,
    rmg_heartbeat_timeout_seconds: float | None = None,
    rmg_start_method: str = "auto",
    requested_phase: str | None = None,
    scenario: dict | None = None,
    conditions: ConditionSpec | None = None,
    reaction_libraries: list[Path] | None = None,
) -> dict:
    """Run a bounded RMG screen and return evidence without a stability verdict."""
    if barrier_threshold <= 0:
        raise ValueError("Barrier threshold must be greater than zero")
    run_dir = Path(run_dir)
    rmg_dir = run_dir / "rmg"
    # Every stage receives an explicit finite reaction-space budget.  The old
    # air-stage exception allowed a target/O2 network to grow indefinitely;
    # that is not a meaningful condition model and prevents ORCA from ever
    # seeing the direct routes RMG did enumerate.
    intrinsic_scope = bool(conditions and (
        "low-pressure" in conditions.scenario or "inert" in conditions.scenario
    ))
    reaction_space_budget = _reaction_space_budget(
        smiles, scenario or {}, intrinsic_scope=intrinsic_scope,
    )
    from .generated_kinetics import validate_generated_library
    reaction_libraries = reaction_libraries or []
    library_manifests = [validate_generated_library(path, temperature_K=temperature, pressure_bar=pressure)
                         for path in reaction_libraries]
    input_file = create_rmg_input(
        "stability", smiles, rmg_dir, temperature=temperature, pressure=pressure,
        max_edge_species=rmg_max_edge_species,
        termination_time=(conditions.target_duration_seconds if conditions else 1e7),
        additional_species=(scenario or {}).get("additional_species"),
        initial_mole_fractions=(scenario or {}).get("initial_mole_fractions"),
        reaction_libraries=reaction_libraries,
        # Association is retained up to the explicit one-collision budget.
        cap_generated_carbon_at_target=False,
        maximum_heavy_atoms=reaction_space_budget["maximum_heavy_atoms"],
        maximum_radical_electrons=reaction_space_budget["maximum_radical_electrons"],
        filter_reactions=reaction_space_budget["filter_reactions"],
        max_objects_per_iteration=rmg_max_objects_per_iteration,
        restart_from_seed=restart_from_seed,
    )
    artifacts = {
        "input": str(input_file), "stdout": str(rmg_dir / "rmg.stdout.log"),
        "stderr": str(rmg_dir / "rmg.stderr.log"), "log": str(rmg_dir / "RMG.log"),
        "chemkin": str(rmg_dir / "chemkin" / "chem_annotated.inp"),
        "species_dictionary": str(rmg_dir / "chemkin" / "species_dictionary.txt"),
        "transport": str(rmg_dir / "chemkin" / "tran.dat"),
        "solver": str(rmg_dir / "solver"),
        "seed": str(rmg_dir / "seed"),
        "collision_rate_violators": str(rmg_dir / "collision_rate_violators.log"),
    }
    evidence = {
        "kind": "rmg_kinetic_network_screen",
        "status": "failed",
        "model_conditions": {
            "reactor": "simpleReactor", "temperature_K": temperature,
            "pressure_bar": pressure, "scenario": (scenario or {}).get("name"),
        },
        "phase_match_to_requested_conditions": "not_established" if requested_phase else "not_applicable",
        "resource_limits": {"walltime": rmg_walltime, "max_processes": rmg_max_processes,
                            "max_iterations": rmg_max_iterations, "max_edge_species": rmg_max_edge_species},
        "model_enlargement": {"max_objects_per_iteration": rmg_max_objects_per_iteration},
        "generated_species_constraints": {
            "maximum_carbon_atoms": None,
            **reaction_space_budget,
            "reason": "RMG may retain only direct one-collision products of the declared target/environment system.",
        },
        "reaction_scope": (
            "intrinsic_target_decomposition_with_declared_bath" if intrinsic_scope
            else "bounded_direct_reactive_network_with_declared_environment"
        ),
        "candidate_routes": [], "artifacts": artifacts,
        "generated_kinetics_libraries": library_manifests,
        "interpretation": "RMG is a bounded homogeneous simple-reactor model screen. Missing routes do not establish persistence.",
    }
    try:
        execution = run_rmg_supervised(input_file, rmg_env=rmg_env, rmg_command=rmg_command,
                                       walltime=rmg_walltime, max_processes=rmg_max_processes,
                                       max_iterations=rmg_max_iterations,
                                       hard_timeout_seconds=rmg_hard_timeout_seconds,
                                       heartbeat_timeout_seconds=rmg_heartbeat_timeout_seconds,
                                       start_method=rmg_start_method)
        Path(artifacts["stdout"]).write_text(execution["stdout"])
        Path(artifacts["stderr"]).write_text(execution["stderr"])
        execution_assessment = assess_rmg_execution(Path(artifacts["log"]), execution)
        if execution["supervisor"]["status"] != "completed":
            execution_assessment.update(status="incomplete", termination_reason=execution["supervisor"]["stop_reason"])
        evidence["execution_assessment"] = execution_assessment
        process_failure_reason = None
        if execution["returncode"] != 0:
            stderr = execution["stderr"]
            if "Could not get database with name: kinetics" in stderr:
                process_failure_reason = (
                    "RMG multiprocessing workers could not access the initialized kinetics database. "
                    "Use one RMG worker in this environment; ORCA cores are independent."
                )
            else:
                process_failure_reason = f"RMG process exited with code {execution['returncode']}; inspect retained stderr."
            # RMG writes Chemkin incrementally.  A solver crash invalidates
            # its propagation, but it must not erase direct, atom-balanced
            # routes which can still be independently checked by ORCA.
            evidence["failure_reason"] = process_failure_reason
        chemkin_file = Path(artifacts["chemkin"])
        if not chemkin_file.is_file():
            evidence.update(status="failed" if process_failure_reason else "incomplete",
                            interpretation="RMG exited without annotated Chemkin output; search coverage is incomplete.")
            return evidence
        parsed = parse_chemkin_annotated(chemkin_file, barrier_threshold=barrier_threshold, label="stability")
        raw_routes = [
            {**route, "screening_threshold_kcal_mol": barrier_threshold}
            for route in parsed["candidate_routes"]
        ]
        dictionary_file = Path(artifacts["species_dictionary"])
        species_dictionary = parse_species_dictionary(dictionary_file) if dictionary_file.is_file() else {}
        if rmg_env:
            # Validate the retained mechanism in the same RMG/Cantera runtime
            # which generated it.  This is deliberately isolated from the
            # app's Python 3.12 environment.
            from .rmg_bridge_client import run_rmg_bridge
            evidence["mechanism_inspection"] = run_rmg_bridge(
                {"operation": "inspect", "chemkin": str(chemkin_file),
                 "dictionary": str(dictionary_file) if dictionary_file.is_file() else None,
                 "transport": artifacts["transport"] if Path(artifacts["transport"]).is_file() else None,
                 "temperature_K": temperature, "pressure_bar": pressure},
                rmg_env=rmg_env,
            )
        provenance = parse_annotated_reaction_provenance(chemkin_file)
        profile = parse_final_solver_profile(Path(artifacts["solver"]), "stability")
        evidence["time_coverage"] = requested_time_coverage(
            profile, conditions.target_duration_seconds if conditions else 0.0,
        )
        evidence["kinetics_validation"] = merge_collision_rate_validation(
            parse_collision_rate_violators(Path(artifacts["collision_rate_violators"])),
            evidence.get("mechanism_inspection"),
            temperature_K=temperature,
            pressure_bar=pressure,
        )
        evidence["candidate_routes"] = [
            enrich_candidate_route(
                {**provenance.get(route.get("printed_reaction_equation", route["reaction_equation"]), {}), **route},
                species_dictionary=species_dictionary,
                profile=profile,
                conditions=conditions,
            ) if conditions else route
            for route in raw_routes
        ]
        evidence["candidate_routes"] = [
            _resolve_rmg_route_endpoints(route, evidence.get("mechanism_inspection"), index)
            for index, route in enumerate(evidence["candidate_routes"])
        ]
        # Retain the full core graph as possible initiation evidence.  Routes
        # that consume the target remain in ``candidate_routes``; formation of
        # a generated radical/intermediate may proceed through a different core
        # reaction and must not disappear merely because its target
        # stoichiometry is zero.
        from .flux_verification import reaction_signature
        network_routes: list[dict] = []
        seen_route_ids = {route.get("route_id") for route in evidence["candidate_routes"]}
        candidate_reaction_signatures = {
            reaction_signature(
                route.get("reaction_equation") or route.get("printed_reaction_equation") or "",
                normalized=True,
            )
            for route in evidence["candidate_routes"]
        }
        for reaction in (evidence.get("mechanism_inspection") or {}).get("reactions", []):
            reactant_labels = list(reaction.get("reactant_labels") or [])
            product_labels = list(reaction.get("product_labels") or [])
            if not reactant_labels or not product_labels:
                continue
            arrow = "<=>" if reaction.get("reversible", True) else "=>"
            directions = [
                (reactant_labels, product_labels, "chemkin_forward"),
                *(([(product_labels, reactant_labels, "reverse_of_chemkin_direction")])
                  if reaction.get("reversible", True) else []),
            ]
            for left_labels, right_labels, direction in directions:
                network_route = _resolve_rmg_route_endpoints(
                    {
                        "reaction": reaction.get("equation"),
                        "reaction_equation": "+".join(left_labels) + arrow + "+".join(right_labels),
                        "printed_reaction_equation": reaction.get("equation"),
                        "direction": direction,
                        "network_reaction_index": reaction.get("reaction_index"),
                        "network_support_route": True,
                        "target_stoichiometry": 0.0,
                        "initiation_status": "network_reachability_pending",
                        "kinetics_comment": reaction.get("kinetics_comment"),
                        "source_library": reaction.get("source_library"),
                    },
                    evidence.get("mechanism_inspection"),
                    len(network_routes),
                )
                physical_signature = reaction_signature(
                    network_route.get("reaction_equation") or "", normalized=True,
                )
                if (physical_signature not in candidate_reaction_signatures
                        and network_route.get("route_id") not in seen_route_ids):
                    network_routes.append(network_route)
                    seen_route_ids.add(network_route.get("route_id"))
        evidence["network_routes"] = network_routes
        if evidence["kinetics_validation"].get("status") == "kinetics_unreliable" and conditions:
            from .verification_plan import build_verification_dependency_plan
            evidence["verification_dependency_plan"] = build_verification_dependency_plan(
                evidence, conditions.as_dict(),
            )
        evidence["species_dictionary"] = species_dictionary
        evidence["solver_profile"] = profile
        relevance_profile = profile
        if rmg_env and execution_assessment["status"] == "completed":
            from .rmg_bridge_client import run_rmg_bridge
            propagation = run_rmg_bridge(
                {
                    "operation": "propagate", "chemkin": str(chemkin_file),
                    "dictionary": str(dictionary_file) if dictionary_file.is_file() else None,
                    "initial_mole_fractions": dict(conditions.composition) if conditions else {"stability": 1.0},
                    "temperature_K": temperature, "pressure_bar": pressure,
                    "target_label": "stability",
                    "target_duration_seconds": conditions.target_duration_seconds if conditions else 1e7,
                    "retention_fraction": conditions.retention_fraction if conditions else 0.95,
                },
                rmg_env=rmg_env,
            )
            evidence["independent_cantera_propagation"] = propagation
            # Keep RMG's species profile for route initiation, but use an
            # amount-based integration for target retention and t95.
            relevance_profile = {**(profile or {}),
                "target_loss_fraction": propagation["target_loss_fraction"],
                "target_time_series": [
                    {"time_seconds": point["time_seconds"],
                     "fraction_remaining": point["target_fraction_remaining"]}
                    for point in propagation["target_profile"]
                ],
                "end_time_seconds": propagation["coverage"]["simulated_seconds"],
            }
        evidence["kinetic_relevance"] = (
            assess_kinetic_relevance(evidence["candidate_routes"], relevance_profile, conditions)
            if conditions else None
        )
        evidence["activation_energy_source_unit"] = parsed["activation_energy_source_unit"]
        evidence["activation_energy_unit"] = parsed["activation_energy_unit"]
        # Chemkin can be saved after a resource stop.  It is useful diagnostic
        # evidence, but never sufficient for a persistence conclusion.
        fully_completed = execution_assessment["status"] == "completed" and evidence["time_coverage"]["complete"]
        evidence["enumeration_status"] = "completed" if fully_completed else (
            "partial_enumeration" if evidence["candidate_routes"] else "unavailable"
        )
        evidence["status"] = "completed" if fully_completed else "incomplete"
        evidence["search_outcome"] = ("candidate_decomposition_pathway_identified" if evidence["candidate_routes"]
                                      else "no_route_found_within_bounded_rmg_model")
        if evidence["status"] != "completed":
            evidence["interpretation"] = (
                "RMG retained diagnostic artifacts but did not complete the requested model and time horizon; "
                "it cannot support a stability or lifetime conclusion. Direct reachable routes in its retained "
                "Chemkin output may be submitted to ORCA as partial enumeration evidence."
            )
        if evidence["candidate_routes"]:
            evidence["orca_verification"] = {
                "status": "required_before_lifetime",
                "eligible_route_indices": [index for index, route in enumerate(evidence["candidate_routes"])
                                           if route.get("kinetic_relevance") == "kinetically_relevant_candidate"],
                "requirements": ["resolve atom-mapped endpoints", "optimize endpoints and TS", "require one imaginary frequency", "verify IRC endpoints", "calculate a verified condition-specific rate"],
            }
    except Exception as error:
        stdout, stderr = getattr(error, "stdout", None), getattr(error, "stderr", None)
        if stdout is not None:
            Path(artifacts["stdout"]).write_text(stdout)
        if stderr is not None:
            Path(artifacts["stderr"]).write_text(stderr)
        evidence["failure_reason"] = str(error)
    return evidence


def collect_rmg_evidence(
    smiles: str,
    run_dir: Path,
    *,
    rmg_retry_attempts: int = 2,
    rmg_hard_timeout_seconds: float | None = None,
    rmg_heartbeat_timeout_seconds: float | None = None,
    rmg_start_method: str = "auto",
    **kwargs,
) -> dict:
    """Run a bounded RMG search and, only if needed, one seed restart.

    Every attempt is retained.  The final record is the last attempt, with a
    controller dossier that makes budget exhaustion and attempt-to-attempt
    changes explicit to the assessment layer.
    """
    from .rmg_controller import convergence_comparison, restart_is_justified, staged_budgets

    if rmg_retry_attempts < 1:
        raise ValueError("rmg_retry_attempts must be at least one")
    resources = {
        "walltime": kwargs.get("rmg_walltime") or "00:00:10:00",
        "max_processes": kwargs.get("rmg_max_processes") or 1,
        "max_iterations": kwargs.get("rmg_max_iterations") or 100,
        "max_edge_species": kwargs.get("rmg_max_edge_species") or 250,
    }
    budgets = staged_budgets(resources, rmg_retry_attempts)
    attempts: list[dict] = []
    restart_seed: Path | None = None
    attempt_kwargs = dict(kwargs)
    for name in ("rmg_walltime", "rmg_max_processes", "rmg_max_iterations", "rmg_max_edge_species",
                 "restart_from_seed", "rmg_hard_timeout_seconds", "rmg_heartbeat_timeout_seconds"):
        attempt_kwargs.pop(name, None)
    for budget in budgets:
        attempt_dir = Path(run_dir) if budget.attempt == 1 else Path(run_dir) / "rmg-restarts" / f"attempt-{budget.attempt:02d}"
        evidence = _collect_rmg_evidence_once(
            smiles, attempt_dir, **attempt_kwargs,
            rmg_walltime=budget.walltime, rmg_max_processes=budget.max_processes,
            rmg_max_iterations=budget.max_iterations, rmg_max_edge_species=budget.max_edge_species,
            restart_from_seed=restart_seed, rmg_hard_timeout_seconds=rmg_hard_timeout_seconds,
            rmg_heartbeat_timeout_seconds=rmg_heartbeat_timeout_seconds,
            rmg_start_method=rmg_start_method,
        )
        evidence["attempt_budget"] = budget.as_dict()
        attempts.append(evidence)
        if evidence.get("status") == "completed" or not restart_is_justified(evidence):
            break
        seed = Path(evidence["artifacts"]["seed"])
        if not seed.is_dir():
            break
        restart_seed = seed
    final = attempts[-1]
    final["controller"] = {
        "status": "completed_converged" if final.get("status") == "completed" else "mechanism_not_converged",
        "attempt_count": len(attempts), "attempt_budgets": [budget.as_dict() for budget in budgets[:len(attempts)]],
        "convergence": convergence_comparison(attempts),
        "attempt_artifacts": [attempt.get("artifacts") for attempt in attempts],
        "interpretation": "A restart is attempted only after genuine RMG resource exhaustion with a retained seed. Partial attempts cannot support persistence.",
    }
    return final


def run_stability_screen(
    smiles: str,
    run_dir: Path,
    *,
    rmg_env: str | None = None,
    rmg_command: str = "rmg.py",
    barrier_threshold: float = 50.0,
    temperature: float = 298.0,
    pressure: float = 1.0,
    rmg_walltime: str | None = None,
    rmg_max_processes: int | None = None,
    rmg_max_iterations: int | None = None,
    rmg_max_edge_species: int | None = None,
    rmg_max_objects_per_iteration: int | None = None,
    rmg_retry_attempts: int = 2,
    rmg_hard_timeout_seconds: float | None = None,
    rmg_heartbeat_timeout_seconds: float | None = None,
    rmg_start_method: str = "auto",
    scenario: str | None = "dry-inert-gas-screen",
    scenario_config: dict | None = None,
    screen_tier: str = "quick-screen",
    target_duration_hours: float = 24.0,
    retention_fraction: float = 0.95,
    light_condition: str = "dark",
    light_model: dict | None = None,
    charge: int = 0,
    multiplicity: int = 1,
    ncores: int = 1,
    method_profile: dict | None = None,
    precomputed_orca_evidence: dict | None = None,
    reaction_libraries: list[Path] | None = None,
    auto_verify_routes: bool = True,
    auto_verify_collision_routes: bool | None = None,
    verification_max_iterations: int = 12,
    verification_timeout_seconds: float | None = 14400.0,
    verification_orientations: int = 3,
    verification_neb_images: int = 8,
) -> dict:
    """Collect ORCA and RMG stability evidence in one retained run directory."""
    # Backward-compatible API alias.  The old flag only handled collision
    # warnings; it now enables/disables the complete generic verifier.
    if auto_verify_collision_routes is not None:
        auto_verify_routes = bool(auto_verify_collision_routes)
    run_dir = Path(run_dir)
    if method_profile is None:
        from .methods import harmonic_method_profile
        method_profile = harmonic_method_profile()
    if scenario_config is None:
        if scenario is None:
            raise ValueError("Either scenario or scenario_config is required")
        scenario_config, resource_config = resolve_stability_configuration(
            scenario, screen_tier, rmg_walltime=rmg_walltime,
            rmg_max_processes=rmg_max_processes, rmg_max_iterations=rmg_max_iterations,
            rmg_max_edge_species=rmg_max_edge_species,
            rmg_max_objects_per_iteration=rmg_max_objects_per_iteration,
        )
    else:
        # Preserve the named resource-tier policy even when the ladder supplies
        # a calculated humid-air composition.
        _, resource_config = resolve_stability_configuration(
            "dry-inert-gas-screen", screen_tier, rmg_walltime=rmg_walltime,
            rmg_max_processes=rmg_max_processes, rmg_max_iterations=rmg_max_iterations,
            rmg_max_edge_species=rmg_max_edge_species,
            rmg_max_objects_per_iteration=rmg_max_objects_per_iteration,
        )
    from .conditions import normalize_target_environment_species
    scenario_config = normalize_target_environment_species(smiles, scenario_config)
    conditions = build_condition_spec(
        scenario_config, temperature_K=temperature, pressure_bar=pressure,
        target_duration_hours=target_duration_hours, retention_fraction=retention_fraction,
        light_condition=light_condition,
        light_model=light_model,
        relative_humidity=scenario_config.get("relative_humidity"),
    )
    condition_contract = conditions.as_dict()
    write_metadata(
        run_dir,
        workflow="orca_rmg_stability_screen",
        smiles=smiles,
        rmg_env=rmg_env,
        rmg_command=rmg_command,
        temperature=temperature,
        pressure=pressure,
        route_energy_ceiling_kcal_mol=barrier_threshold,
        stability_scenario=scenario_config,
        condition_contract=condition_contract,
        rmg_resource_tier=resource_config,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        harmonic_method_profile=method_profile,
        status="running",
        interpretation="Combined ORCA local-minimum and condition-specific RMG kinetic evidence; not a safety, shelf-life, or compatibility assessment.",
    )
    if precomputed_orca_evidence is not None:
        print("[STORCA] Reusing the validated parent ORCA minimum for this ladder condition", flush=True)
        orca = _reuse_orca_evidence(
            precomputed_orca_evidence,
            smiles,
            charge=charge,
            multiplicity=multiplicity,
            method_profile=method_profile,
        )
    else:
        orca = collect_orca_evidence(
            smiles, run_dir, charge=charge, multiplicity=multiplicity,
            ncores=ncores, method_profile=method_profile,
        )
    rmg = collect_rmg_evidence(smiles, run_dir, rmg_env=rmg_env, rmg_command=rmg_command,
                               barrier_threshold=barrier_threshold, temperature=temperature, pressure=pressure,
                               rmg_walltime=resource_config["walltime"], rmg_max_processes=resource_config["max_processes"],
                               rmg_max_iterations=resource_config["max_iterations"], rmg_max_edge_species=resource_config["max_edge_species"],
                               rmg_max_objects_per_iteration=resource_config["max_objects_per_iteration"],
                               scenario=scenario_config, conditions=conditions, reaction_libraries=reaction_libraries,
                               rmg_retry_attempts=rmg_retry_attempts,
                               rmg_hard_timeout_seconds=rmg_hard_timeout_seconds,
                               rmg_heartbeat_timeout_seconds=rmg_heartbeat_timeout_seconds,
                               rmg_start_method=rmg_start_method)
    risk_flags: list[str] = []
    if orca["status"] != "completed":
        risk_flags.append("orca_evidence_incomplete")
    elif not orca["local_minimum"]:
        risk_flags.append("orca_not_a_local_minimum")
    if rmg["status"] != "completed":
        risk_flags.append("rmg_evidence_incomplete")
    elif rmg["candidate_routes"]:
        route_contexts = {route["route_context"] for route in rmg["candidate_routes"]}
        if "unimolecular_decomposition" in route_contexts:
            risk_flags.append("candidate_unimolecular_decomposition")
        if "self_reaction" in route_contexts:
            risk_flags.append("candidate_self_reaction")
        if "radical_or_impurity_initiated" in route_contexts:
            risk_flags.append("candidate_radical_or_impurity_initiated_pathway")
        if "co_reactant_dependent" in route_contexts:
            risk_flags.append("candidate_co_reactant_dependent_pathway")
    from .flux_verification import assess_target_loss_commitment
    retained_flux = (rmg.get("independent_cantera_propagation") or {}).get("reaction_flux_attribution") or {}
    target_loss_commitment = assess_target_loss_commitment(retained_flux, conditions.target_label)
    rmg["target_loss_commitment"] = target_loss_commitment
    if (target_loss_commitment.get("committed_target_loss_kmol") or 0.0) > 0.0:
        risk_flags.append("modeled_committed_target_loss")
    kinetic_status = (rmg.get("kinetic_relevance") or {}).get("status")
    advisories: list[str] = []
    # Completion and coverage are gates, not advisories.  In particular, RMG
    # can save a long solver profile while model generation was interrupted;
    # route relevance from that partial mechanism cannot establish persistence.
    if orca["status"] != "completed" or rmg["status"] != "completed":
        assessment_status = "incomplete_evidence"
        assessment_reason = "Required ORCA or RMG evidence did not complete under the declared condition contract; no stability conclusion is available."
    elif orca.get("local_minimum") is not True:
        assessment_status = "incomplete_evidence"
        assessment_reason = (
            "The parent structure was not validated as a local minimum by the retained ORCA frequency calculation. "
            "A model no-loss result cannot override that failed structural gate."
        )
    elif (rmg.get("kinetics_validation") or {}).get("status") not in {"passed", "kinetics_unreliable"}:
        assessment_status = "incomplete_evidence"
        assessment_reason = (
            "The retained mechanism could not validate every applicable bimolecular rate against a finite "
            "condition-specific collision ceiling. No stability or lifetime conclusion is available."
        )
    elif (rmg.get("kinetics_validation") or {}).get("status") == "kinetics_unreliable":
        if kinetic_status == "kinetically_relevant_candidate":
            assessment_status = "likely_air_reactive_lifetime_unverified"
            assessment_reason = "RMG found a reachable, target-consuming route under the declared conditions, but a controlling rate violates the gas collision limit. Treat this as a qualitative reactivity warning; ORCA repair is required before any lifetime is reported."
        else:
            assessment_status = "kinetics_unreliable"
            assessment_reason = "RMG reported collision-rate-limit violations in the retained core mechanism. Its modeled loss and t95 require ORCA route verification and a repaired kinetic model."
    elif kinetic_status in {"no_credible_initiation_demonstrated", "reachable_but_below_loss_threshold"}:
        advisories = list(risk_flags)
        risk_flags = [flag for flag in risk_flags if not flag.startswith("candidate_")]
        assessment_status = "no_loss_detected_in_rmg_model"
        assessment_reason = "No credible initiation and target loss above the retention threshold were demonstrated in the completed retained RMG model. This is model-scoped, not a universal stability claim."
    elif kinetic_status == "kinetically_relevant_candidate":
        assessment_status = "orca_verification_required"
        assessment_reason = "RMG identified a kinetically relevant candidate route, but this is not an instability verdict or lifetime. An ORCA TS/frequency/IRC route verification is required first."
    elif not risk_flags:
        assessment_status = "no_loss_detected_in_rmg_model"
        assessment_reason = "ORCA found a local minimum and the completed bounded RMG model found no target-loss route crossing the declared retention threshold. This is model-scoped."
    elif "candidate_unimolecular_decomposition" in risk_flags:
        assessment_status = "candidate_intrinsic_instability"
        assessment_reason = "RMG found a candidate low-barrier unimolecular route. Verify the route and barrier with higher-level calculations before drawing a condition-scoped conclusion."
    elif "candidate_self_reaction" in risk_flags:
        assessment_status = "candidate_self_reaction_risk"
        assessment_reason = "RMG found a candidate route requiring multiple target molecules. Review concentration and phase assumptions, then verify material routes with higher-level calculations."
    elif "candidate_radical_or_impurity_initiated_pathway" in risk_flags:
        assessment_status = "condition_dependent_risk"
        assessment_reason = "RMG found a route requiring a radical or impurity co-reactant. It is not evidence of spontaneous decomposition; assess whether the co-reactant is credible under the stated conditions."
    elif "candidate_co_reactant_dependent_pathway" in risk_flags:
        assessment_status = "condition_dependent_risk"
        assessment_reason = "RMG found a route requiring a co-reactant. It is not evidence of spontaneous decomposition; assess whether that co-reactant is credible under the stated conditions."
    else:
        assessment_status = "incomplete_evidence"
        assessment_reason = "At least one required ORCA or RMG evidence source is incomplete or does not support a local-minimum structure."
    result = {
        "schema_version": 2,
        "kind": "combined_orca_rmg_stability_screen",
        "assessment": {"status": assessment_status, "reason": assessment_reason},
        "risk_flags": risk_flags,
        "advisories": advisories,
        "orca_evidence": orca,
        "rmg_evidence": rmg,
        "target_loss_commitment": target_loss_commitment,
        "kinetic_lifetime": {
            "retention_fraction": retention_fraction,
            "estimated_time_to_retention_seconds": None,
            "screening_estimated_time_to_retention_seconds": (rmg.get("kinetic_relevance") or {}).get("estimated_time_to_retention_seconds"),
            "interpretation": "The RMG screen estimate is retained for triage only. This field remains null until ORCA verifies the selected route by TS, one imaginary mode, IRC endpoints, and condition-specific kinetics.",
        },
        "condition_contract": condition_contract,
        "scenario": {key: value for key, value in scenario_config.items() if key != "additional_species"},
        "rmg_resource_tier": resource_config,
        "decomposition_reactions": [route["reaction"] for route in rmg["candidate_routes"]],
        "reaction_barriers_kcal_mol": [route["activation_energy_kcal_mol"] for route in rmg["candidate_routes"]],
        "provenance": {
        "source": "ORCA geometry/frequency workflow plus RMG generated model and annotated Chemkin parser",
        "network_accessed": False,
        "conditions": {"temperature_K": temperature, "pressure_bar": pressure},
        "route_energy_ceiling_kcal_mol": barrier_threshold,
        },
    }
    result["limitations"] = [
        "This is a condition-specific model screen, not a laboratory safety assessment.",
        "No detected route is not evidence of shelf-life, compatibility, or absence of hazards.",
        scenario_config["model_applicability"],
    ]
    result_path = run_dir / "stability.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    rmg_has_orca_queue = (
        rmg.get("status") == "completed"
        or rmg.get("enumeration_status") == "partial_enumeration"
    )
    if (auto_verify_routes and orca.get("status") == "completed" and orca.get("local_minimum") is True
            and rmg_has_orca_queue):
        print("[STORCA] Generic flux-guided ORCA/Arkane verification planning started", flush=True)
        try:
            from .verification_adapters import make_arkane_rate_builder, make_orca_route_verifier
            from .verification_engine import VerificationEngineConfig, run_verification_engine

            base_libraries = list(reaction_libraries or [])
            verify_route = make_orca_route_verifier(
                smiles,
                ncores=ncores,
                method_profile=method_profile,
                condition_contract=condition_contract,
                timeout_seconds=verification_timeout_seconds,
                orientations=verification_orientations,
                neb_images=verification_neb_images,
            )
            build_rate = make_arkane_rate_builder(
                condition_contract,
                method_profile=method_profile,
                rmg_env=rmg_env,
            )

            def rerun_model(library_paths: list[Path], verification_dir: Path) -> dict:
                libraries: list[Path] = []
                # RMG keeps the first already-registered reaction when later
                # libraries/families contain the same elementary step.  Put
                # verified replacements first so they actually replace a base
                # estimate; exact post-rerun rate matching still fails closed.
                for library in [*library_paths, *base_libraries]:
                    resolved_library = Path(library)
                    if resolved_library not in libraries:
                        libraries.append(resolved_library)
                print(
                    f"[STORCA] Full RMG repair rerun started with {len(library_paths)} verified replacement(s)",
                    flush=True,
                )
                repaired = collect_rmg_evidence(
                    smiles,
                    verification_dir,
                    rmg_env=rmg_env,
                    rmg_command=rmg_command,
                    barrier_threshold=barrier_threshold,
                    temperature=temperature,
                    pressure=pressure,
                    rmg_walltime=resource_config["walltime"],
                    rmg_max_processes=resource_config["max_processes"],
                    rmg_max_iterations=resource_config["max_iterations"],
                    rmg_max_edge_species=resource_config["max_edge_species"],
                    rmg_max_objects_per_iteration=resource_config["max_objects_per_iteration"],
                    scenario=scenario_config,
                    conditions=conditions,
                    reaction_libraries=libraries,
                    rmg_retry_attempts=rmg_retry_attempts,
                    rmg_hard_timeout_seconds=rmg_hard_timeout_seconds,
                    rmg_heartbeat_timeout_seconds=rmg_heartbeat_timeout_seconds,
                    rmg_start_method=rmg_start_method,
                )
                print(f"[STORCA] Full RMG repair rerun finished: {repaired.get('status')}", flush=True)
                return {"rmg_evidence": repaired, "condition_contract": condition_contract}

            verification = run_verification_engine(
                rmg,
                condition_contract,
                run_dir,
                verify_route=verify_route,
                build_rate=build_rate,
                rerun_model=rerun_model,
                config=VerificationEngineConfig(max_iterations=verification_max_iterations),
            )
        except Exception as error:
            verification = {
                "schema_version": 1,
                "kind": "verification_summary",
                "status": "verification_incomplete",
                "reason": f"Generic verification could not start or complete: {type(error).__name__}: {error}",
                "condition_contract": condition_contract,
                "condition_contract_match": True,
                "estimated_time_to_retention_seconds": None,
            }
        result["verification_summary"] = verification
        verification_status = verification.get("status")
        if verification_status == "verified_t95_converged":
            lifetime = verification.get("estimated_time_to_retention_seconds")
            result["assessment"] = {
                "status": "orca_verified_condition_dependent_instability",
                "reason": "Flux-controlling ORCA/Arkane rates were exactly reinjected and the repaired condition model converged to the retention threshold.",
            }
            result["kinetic_lifetime"].update(
                estimated_time_to_retention_seconds=lifetime,
                interpretation="Condition-specific t95 from the converged repaired mechanism after ORCA/Arkane verification.",
            )
        elif verification_status == "verified_below_loss_threshold":
            result["assessment"] = {
                "status": "verified_route_below_loss_threshold",
                "reason": "Verified controlling-route rates and conservative unresolved-route bounds remain below the retention threshold in this condition model.",
            }
        elif verification_status == "no_target_loss_in_completed_rmg_model":
            result["assessment"] = {
                "status": "no_loss_detected_in_rmg_model",
                "reason": verification.get("reason"),
            }
        else:
            result["assessment"] = {
                "status": "orca_verification_incomplete",
                "reason": verification.get("reason") or "Generic route/rate verification did not converge.",
            }
        print(f"[STORCA] Generic verification finished: {verification_status}", flush=True)
        result_path.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
        assessment_status = result["assessment"]["status"]
    elif not auto_verify_routes:
        result["verification_summary"] = {
            "schema_version": 1,
            "kind": "verification_summary",
            "status": "verification_disabled",
            "reason": "Automatic generic ORCA/Arkane route verification was disabled by the caller.",
            "condition_contract": condition_contract,
            "condition_contract_match": True,
            "estimated_time_to_retention_seconds": None,
        }
        result_path.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
    else:
        failed_gates = []
        if orca.get("status") != "completed":
            failed_gates.append("orca_evidence_not_completed")
        elif orca.get("local_minimum") is not True:
            failed_gates.append("orca_parent_not_a_validated_local_minimum")
        if not rmg_has_orca_queue:
            failed_gates.append("rmg_evidence_not_completed")
        result["verification_summary"] = {
            "schema_version": 1,
            "kind": "verification_summary",
            "status": "verification_incomplete",
            "reason": "Automatic route verification was not started: " + ", ".join(failed_gates),
            "condition_contract": condition_contract,
            "condition_contract_match": True,
            "estimated_time_to_retention_seconds": None,
            "claim_scope": "verification_incomplete_no_condition_lifetime",
        }
        result_path.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
    terminal_statuses = {
        "no_loss_detected_in_rmg_model", "verified_route_below_loss_threshold",
        "orca_verified_condition_dependent_instability",
    }
    metadata_status = "completed" if assessment_status in terminal_statuses else "completed_with_risk_or_incomplete_evidence"
    write_metadata(run_dir, result_json=str(result_path), status=metadata_status, assessment_status=assessment_status)
    return {**result, "result_json": result_path, "chemkin_file": Path(rmg["artifacts"]["chemkin"])}
def _resolve_rmg_route_endpoints(route: dict, mechanism_inspection: dict | None, index: int) -> dict:
    """Attach exact retained species and an explicit electronic-state contract.

    RMG does not retain a reaction atom map or a unique total-spin surface.
    Those are resolved later by the bounded generic path engine, but endpoint
    identities, stoichiometry, charge, and all common spin surfaces are fixed
    here so no downstream ORCA job can silently assume neutral singlets.
    """
    equation = str(route.get("reaction_equation") or route.get("reaction") or "")
    # Candidate-list positions can change after a repaired RMG rerun.  Bind
    # replacement evidence to a deterministic, direction-aware normalized
    # reaction identity instead of an ephemeral ``rmg-4``-style index.
    from .rmg_evidence import normalize_rmg_label
    arrow_for_id = "<=>" if "<=>" in equation else ("=>" if "=>" in equation else None)
    if arrow_for_id:
        identifier_sides = []
        for raw_side in equation.split(arrow_for_id, 1):
            counts: dict[str, int] = {}
            for token in raw_side.split("+"):
                match = re.match(r"^\s*(\d+)\s+(.+?)\s*$", token)
                coefficient = int(match.group(1)) if match else 1
                label = normalize_rmg_label((match.group(2) if match else token).strip())
                counts[label] = counts.get(label, 0) + coefficient
            identifier_sides.append("+".join(
                f"{coefficient if coefficient != 1 else ''}{label}"
                for label, coefficient in sorted(counts.items())
            ))
        route_identity = arrow_for_id.join(identifier_sides)
    else:
        route_identity = equation.strip()
    stable_route_id = "rmg-" + hashlib.sha256(route_identity.encode("utf-8")).hexdigest()[:12]
    resolved = {
        **route,
        "route_id": str(route.get("route_id") or stable_route_id),
        "stable_route_identity": route_identity,
        "source_route_index": index,
    }
    records = (mechanism_inspection or {}).get("species", [])
    by_name = {}
    for record in records:
        for name in (record.get("label"), record.get("chemkin_identifier")):
            if name:
                by_name[str(name)] = record
    arrow = "<=>" if "<=>" in equation else ("=>" if "=>" in equation else None)
    if not arrow:
        resolved["orca_verification"] = {"status": "unresolved_endpoints", "reason": "RMG reaction equation could not be parsed."}
        return resolved
    def side(value: str) -> list[dict]:
        entries = []
        for token in value.split("+"):
            coefficient_match = re.match(r"^\s*(\d+)\s+(.+?)\s*$", token)
            coefficient = int(coefficient_match.group(1)) if coefficient_match else 1
            label = (coefficient_match.group(2) if coefficient_match else token).strip()
            item = by_name.get(label)
            for copy_index in range(coefficient):
                smiles = item.get("smiles") if item else None
                charge = item.get("net_charge") if item else None
                if charge is None and smiles:
                    try:
                        from rdkit import Chem
                        molecule = Chem.MolFromSmiles(smiles)
                        charge = Chem.GetFormalCharge(molecule) if molecule is not None else None
                    except Exception:
                        charge = None
                entries.append({
                    "label": label,
                    "stoichiometric_copy": copy_index + 1 if coefficient > 1 else None,
                    "smiles": smiles,
                    "charge": charge,
                    "multiplicity": item.get("multiplicity") if item else None,
                })
        return entries
    left, right = equation.split(arrow, 1)
    reactants, products = side(left), side(right)
    try:
        from .flux_verification import reaction_signature
        signature = reaction_signature(equation, normalized=True)
        reaction_records = (mechanism_inspection or {}).get("reactions", [])
        rate_matches = [
            item for item in reaction_records
            if reaction_signature(item.get("equation", ""), normalized=True) == signature
        ]
    except Exception:
        rate_matches = []
    if len(rate_matches) == 1:
        retained = rate_matches[0]
        reverse_direction = route.get("direction") == "reverse_of_chemkin_direction"
        rate_key = "evaluated_reverse_rate_coefficient_si" if reverse_direction else "evaluated_forward_rate_coefficient_si"
        units_key = "evaluated_reverse_rate_units" if reverse_direction else "evaluated_rate_units"
        collision_key = "reverse_collision_limit_si" if reverse_direction else "forward_collision_limit_si"
        collision_units_key = "reverse_collision_limit_units" if reverse_direction else "forward_collision_limit_units"
        resolved["retained_rate_evidence"] = {
            "reaction_index": retained.get("reaction_index"),
            "direction": "reverse" if reverse_direction else "forward",
            "molecularity": len(reactants),
            "rate_coefficient_si": retained.get(rate_key),
            "rate_units": retained.get(units_key),
            "collision_limit_si": retained.get(collision_key),
            "collision_limit_units": retained.get(collision_units_key),
            "collision_limit_source": retained.get("collision_limit_source"),
            "evaluated_at": retained.get("evaluated_at"),
            "source_library": retained.get("source_library"),
            "kinetics_comment": retained.get("kinetics_comment"),
        }
    complete = all(
        item["smiles"] and item["multiplicity"] and item["charge"] is not None
        for item in reactants + products
    )
    resolved["resolved_endpoints"] = {"reactants": reactants, "products": products}
    if complete:
        try:
            from rdkit import Chem

            def inventory(items: list[dict]) -> dict[int, int]:
                counts: dict[int, int] = {}
                for item in items:
                    molecule = Chem.MolFromSmiles(str(item["smiles"]))
                    if molecule is None:
                        raise ValueError(f"invalid endpoint SMILES: {item['smiles']}")
                    # RMG/RDKit endpoint SMILES commonly carry hydrogens
                    # implicitly. Formula balance must include them: a
                    # heavy-atom-only check can accept an impossible H-loss
                    # reaction before ORCA has a chance to reject it.
                    for atom in Chem.AddHs(molecule).GetAtoms():
                        atomic_number = atom.GetAtomicNum()
                        counts[atomic_number] = counts.get(atomic_number, 0) + 1
                return counts

            reactant_inventory = inventory(reactants)
            product_inventory = inventory(products)
            resolved["atom_balance"] = {
                "status": "passed" if reactant_inventory == product_inventory else "failed",
                "reactants": reactant_inventory,
                "products": product_inventory,
            }
        except Exception as error:
            resolved["atom_balance"] = {"status": "unresolved", "reason": str(error)}
    else:
        resolved["atom_balance"] = {"status": "unresolved", "reason": "endpoint structures are incomplete"}
    if complete:
        def allowed_multiplicities(values: list[int]) -> set[int]:
            # Couple twice-spin integers exactly, avoiding half-integer floats.
            coupled = {0}
            for multiplicity in values:
                twice_spin = int(multiplicity) - 1
                coupled = {
                    total + 1
                    for prior in coupled
                    for total in range(abs(prior - twice_spin), prior + twice_spin + 1, 2)
                }
                # ``coupled`` above stores multiplicities; convert back to 2S
                # before the next fragment is coupled.
                coupled = {multiplicity_value - 1 for multiplicity_value in coupled}
            return {twice_spin + 1 for twice_spin in coupled}

        resolved["reactant_smiles"] = [item["smiles"] for item in reactants]
        resolved["product_smiles"] = [item["smiles"] for item in products]
        resolved["reactant_charges"] = [int(item["charge"]) for item in reactants]
        resolved["product_charges"] = [int(item["charge"]) for item in products]
        resolved["reactant_multiplicities"] = [item["multiplicity"] for item in reactants]
        resolved["product_multiplicities"] = [item["multiplicity"] for item in products]
        reactant_charge = sum(resolved["reactant_charges"])
        product_charge = sum(resolved["product_charges"])
        common_surfaces = sorted(
            allowed_multiplicities(resolved["reactant_multiplicities"])
            & allowed_multiplicities(resolved["product_multiplicities"])
        )
        resolved["charge"] = reactant_charge
        resolved["multiplicity"] = common_surfaces[0] if common_surfaces else 1
        resolved["surface_multiplicities"] = common_surfaces
        resolved["electronic_state_validation"] = {
            "status": "resolved" if reactant_charge == product_charge and common_surfaces else "inconsistent",
            "reactant_total_charge": reactant_charge,
            "product_total_charge": product_charge,
            "common_surface_multiplicities": common_surfaces,
        }
        complete = reactant_charge == product_charge and bool(common_surfaces)
    resolved["orca_verification"] = {
        "status": "ready_for_bounded_generic_mapping" if complete else "unresolved_endpoints_or_surface",
        "endpoint_structures_resolved": complete,
        "atom_mapping_resolved": False,
        "can_execute": complete,
        "reason": (
            "Endpoint identity, stoichiometry, charge, and common spin surfaces are resolved; "
            "the generic path engine must still produce an unambiguous atom map and computed path."
            if complete else
            "One or more endpoint identities, charges, multiplicities, or common spin surfaces could not be resolved."
        ),
    }
    return resolved
