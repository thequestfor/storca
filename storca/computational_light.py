"""Bounded ORCA-to-RMG computational sunlight model."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from src.molecule_tools import smiles_to_xyz
from src.parser import parse_orca_energy

from .generated_kinetics import write_photolysis_library
from .conditions import normalize_target_environment_species
from .decomposition_visuals import render_candidate_storyboard
from .light_model import ComputationalLightModel, energetic_accessibility
from .photo_routes import generic_photo_candidates, rank_photo_routes
from .photolysis import integrate_photolysis_rate
from .photophysics import EV_NM, oscillator_strength_cross_sections, run_photostability_screen
from .runs import write_metadata
from .workflow import run_optimization_and_frequency


HARTREE_TO_EV = 27.211386245988


def _route_mechanism_explanation(route: dict) -> dict:
    """Return an honest, user-facing explanation of a photon-route hypothesis."""
    state = route.get("dominant_excited_state") or {}
    products = route.get("fragment_radical_smiles") or []
    return {
        "status": "candidate_mechanism_only",
        "summary": (
            "A bright ORCA TD-DFT state has enough vertical energy to make this RMG-independent "
            "homolysis candidate energetically accessible. This is not yet proof that the excited "
            "state dissociates rather than fluorescing, internally converting, or crossing surfaces."
        ),
        "proposed_elementary_change": {
            "kind": route.get("kind"),
            "bond_atom_indices": route.get("bond_atom_indices"),
            "products": products,
        },
        "vertical_excitation": {
            "state": state.get("state"), "wavelength_nm": state.get("wavelength_nm"),
            "energy_eV": state.get("energy_eV"), "oscillator_strength": state.get("oscillator_strength"),
        },
        "ground_state_fragment_energy_eV": route.get("route_energy_eV"),
        "required_next_evidence": [
            "ORCA excited-state scan along the declared dissociation coordinate with state tracking",
            "validated excited-state product connection or dissociative profile",
            "competing radiative/non-radiative channel evidence before assigning a quantum yield",
        ],
    }


def _multiplicity(smiles: str) -> int:
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid route-product SMILES: {smiles}")
    return max(1, sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()) + 1)


def _route_energies(candidates: list[dict], parent_energy_hartree: float, workdir: Path,
                    *, ncores: int, method_keywords: list[str] | None) -> tuple[dict[str, float], dict[str, dict]]:
    """Calculate capped radical-product energies for generic homolysis candidates."""
    energies, artifacts = {}, {}
    for candidate in candidates:
        product_energies = []
        route_dir = workdir / "routes" / candidate["route_id"]
        try:
            for index, product_smiles in enumerate(candidate["fragment_radical_smiles"], start=1):
                product_dir = route_dir / f"product-{index}"
                product_dir.mkdir(parents=True, exist_ok=True)
                xyz = product_dir / "input.xyz"
                smiles_to_xyz(product_smiles, xyz)
                result = run_optimization_and_frequency(xyz, product_dir, multiplicity=_multiplicity(product_smiles),
                                                       ncores=ncores, run_frequency=False, method_keywords=method_keywords)
                product_energies.append(parse_orca_energy(result["optimization"]["out"]))
            energies[candidate["route_id"]] = (sum(product_energies) - parent_energy_hartree) * HARTREE_TO_EV
            artifacts[candidate["route_id"]] = {"status": "completed", "directory": str(route_dir),
                                                 "product_energies_hartree": product_energies}
        except Exception as error:
            artifacts[candidate["route_id"]] = {"status": "failed", "directory": str(route_dir), "failure_reason": str(error)}
    return energies, artifacts


def _run_computational_light_model(smiles: str, run_dir: Path, *, sunlight_spectrum: Path | None,
                                  rmg_env: str | None = "rmg_env", charge: int = 0, multiplicity: int = 1,
                                  ncores: int = 1, nroots: int = 16,
                                  model: ComputationalLightModel = ComputationalLightModel(),
                                  method_keywords: list[str] | None = None) -> dict:
    """Calculate bounded photo-route projections and optional RMG libraries.

    If no source spectrum is supplied this still discovers/calculates routes,
    but deliberately cannot create a photon rate or RMG library.
    """
    run_dir = Path(run_dir)
    screen = run_photostability_screen(smiles, run_dir / "tddft", charge=charge, multiplicity=multiplicity,
                                       ncores=ncores, nroots=nroots, method_keywords=method_keywords)
    result = {
        "schema_version": 2, "kind": "computational_light_model", "smiles": smiles, "status": "failed",
        "model_contract": model.as_dict(), "source_spectrum": str(sunlight_spectrum) if sunlight_spectrum else None,
        "tddft_screen": str(screen["result_json"]), "profiles": {},
        "limitations": [
            "Reactive photon fractions are low/central/high model priors, not calculated quantum yields.",
            "Candidate routes are generic single-bond homolyses only in v2.",
            "Generated rates are computational screen projections, not experimental photolysis rates.",
            "No admitted channel is missing evidence, not proof of photostability.",
        ],
    }
    if screen["status"] != "completed":
        result["failure_reason"] = screen.get("failure_reason", "TD-DFT screen did not complete")
    else:
        candidates = generic_photo_candidates(smiles, maximum_routes=model.maximum_routes)
        parent_energy = parse_orca_energy(Path(screen["artifacts"]["tddft_output"]).with_name("opt.out"))
        energies, route_artifacts = _route_energies(candidates, parent_energy, run_dir, ncores=ncores,
                                                    method_keywords=method_keywords)
        source_window = None
        if sunlight_spectrum:
            with Path(sunlight_spectrum).open(newline="") as handle:
                source_grid = [float(row["wavelength_nm"]) for row in csv.DictReader(handle)]
            if source_grid:
                source_window = (max(min(source_grid), model.solar_window_nm[0]),
                                 min(max(source_grid), model.solar_window_nm[1]))
                if source_window[0] >= source_window[1]:
                    raise ValueError("Declared light source does not overlap the configured photochemical window")
        ranked = rank_photo_routes(candidates, screen["transitions"], energies, model, source_window_nm=source_window)
        illuminated = [
            state for state in screen["transitions"]
            if source_window and source_window[0] <= float(state["wavelength_nm"]) <= source_window[1]
        ]
        bright_outside = [
            state for state in screen["transitions"]
            if source_window and not (source_window[0] <= float(state["wavelength_nm"]) <= source_window[1])
            and float(state.get("oscillator_strength") or 0.0) > 0.0
        ]
        result["source_coverage_assessment"] = {
            "source_window_nm": list(source_window) if source_window else None,
            "transition_count_inside_source_window": len(illuminated),
            "bright_transition_count_outside_source_window": len(bright_outside),
            "nearest_bright_transition_nm": (
                min((float(item["wavelength_nm"]) for item in bright_outside), default=None)
            ),
            "status": (
                "source_does_not_cover_computed_bright_transition" if source_window and not illuminated and bright_outside
                else "source_and_transition_window_overlap" if illuminated else "source_window_not_declared"
            ),
            "interpretation": (
                "No computed transition lies within the supplied source window; this is a source/vertical-spectrum "
                "mismatch, not evidence that the molecule is photostable."
                if source_window and not illuminated and bright_outside else
                "The source overlaps computed vertical transitions; an excited-state dissociation path is still required."
                if illuminated else "No source spectrum was supplied."
            ),
        }
        for route in ranked:
            route["mechanism_explanation"] = _route_mechanism_explanation(route)
            if route["status"] == "computed":
                route["candidate_storyboard"] = str(render_candidate_storyboard(
                    smiles, route["fragment_radical_smiles"], run_dir / "visuals" / f"{route['route_id']}.png",
                    title=f"{route['route_id']} ({route['kind']})"))
        result.update(status="completed", candidate_routes=ranked, route_calculations=route_artifacts,
                      absorption_cross_sections=str(run_dir / "absorption-cross-sections.csv"))
        spectrum_wavelengths = source_grid if sunlight_spectrum else None
        cross_sections = oscillator_strength_cross_sections(screen["transitions"], wavelengths_nm=spectrum_wavelengths)
        cross_section_csv = Path(result["absorption_cross_sections"])
        cross_section_csv.parent.mkdir(parents=True, exist_ok=True)
        with cross_section_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["wavelength_nm", "absorption_cross_section_cm2_molecule"])
            writer.writeheader(); writer.writerows(cross_sections)
        if sunlight_spectrum:
            if not Path(sunlight_spectrum).is_file():
                raise FileNotFoundError(f"Sunlight spectrum not found: {sunlight_spectrum}")
            computed = [route for route in ranked if route["status"] == "computed"
                        and route["accessibility"] >= model.minimum_route_accessibility]
            for profile, fraction in model.profiles().items():
                libraries, channels = [], []
                spectral_scores: list[list[float]] = []
                for point in cross_sections:
                    photon_energy_eV = EV_NM / float(point["wavelength_nm"])
                    scores = [energetic_accessibility(
                        photon_energy_eV, float(route["route_energy_eV"]),
                        softening_eV=model.energy_softening_eV,
                    ) for route in computed]
                    denominator = max(1.0, sum(scores))
                    spectral_scores.append([fraction * score / denominator for score in scores])
                for route_index, route in enumerate(computed):
                    route_spectrum = run_dir / "spectral-evidence" / profile / f"{route['route_id']}.csv"
                    route_spectrum.parent.mkdir(parents=True, exist_ok=True)
                    route_rows = [
                        {**point, "quantum_yield": spectral_scores[point_index][route_index]}
                        for point_index, point in enumerate(cross_sections)
                    ]
                    with route_spectrum.open("w", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=[
                            "wavelength_nm", "absorption_cross_section_cm2_molecule", "quantum_yield",
                        ])
                        writer.writeheader(); writer.writerows(route_rows)
                    integral = integrate_photolysis_rate(Path(sunlight_spectrum), route_spectrum)
                    if integral["photolysis_rate_constant_s-1"] <= 0:
                        channels.append({"route_id": route["route_id"],
                                         "quantum_yield_range": [min(row["quantum_yield"] for row in route_rows),
                                                                  max(row["quantum_yield"] for row in route_rows)],
                                         "photolysis_rate_constant_s-1": 0.0, "status": "no_source_overlap"})
                        continue
                    products = [(f"{route['route_id']}-p{index}", value)
                                for index, value in enumerate(route["fragment_radical_smiles"], start=1)]
                    library = write_photolysis_library(
                        run_dir / "libraries" / profile / route["route_id"], route_id=f"{profile}-{route['route_id']}",
                        reactant_label="stability", reactant_smiles=smiles, products=products,
                        photolysis_rate_constant_s_1=integral["photolysis_rate_constant_s-1"],
                        photolysis_evidence=route_spectrum, rmg_env=rmg_env,
                    )
                    libraries.append(library["library"])
                    channels.append({"route_id": route["route_id"],
                                     "quantum_yield_model": "wavelength_resolved_energy_gating",
                                     "quantum_yield_range": [min(row["quantum_yield"] for row in route_rows),
                                                              max(row["quantum_yield"] for row in route_rows)],
                                     "photolysis_rate_constant_s-1": integral["photolysis_rate_constant_s-1"],
                                     "spectral_evidence": str(route_spectrum), "library": library["library"]})
                result["profiles"][profile] = {"reactive_photon_fraction": fraction, "channels": channels,
                                                "reaction_libraries": libraries}
        else:
            result["status"] = "source_spectrum_required"
    output = run_dir / "computational-light.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
    write_metadata(run_dir, workflow="computational_light_model", result_json=str(output), light_model=result["model_contract"])
    return {**result, "result_json": output}


def run_computational_light_model(smiles: str, run_dir: Path, *, sunlight_spectrum: Path | None,
                                  rmg_env: str | None = "rmg_env", charge: int = 0, multiplicity: int = 1,
                                  ncores: int = 1, nroots: int = 16,
                                  model: ComputationalLightModel = ComputationalLightModel(),
                                  method_keywords: list[str] | None = None) -> dict:
    """Run the bounded model and retain a diagnostic JSON even after late failure."""
    try:
        return _run_computational_light_model(
            smiles, run_dir, sunlight_spectrum=sunlight_spectrum, rmg_env=rmg_env, charge=charge,
            multiplicity=multiplicity, ncores=ncores, nroots=nroots, model=model, method_keywords=method_keywords,
        )
    except Exception as error:
        run_dir = Path(run_dir); run_dir.mkdir(parents=True, exist_ok=True)
        result = {"schema_version": 2, "kind": "computational_light_model", "smiles": smiles,
                  "status": "failed", "failure_reason": str(error), "model_contract": model.as_dict()}
        output = run_dir / "computational-light.json"
        output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        write_metadata(run_dir, workflow="computational_light_model", result_json=str(output), status="failed")
        return {**result, "result_json": output}


def simulate_computational_light_profiles(light_result: dict, run_dir: Path, *, rmg_env: str | None = "rmg_env",
                                          scenario: str = "ambient-air-gas-screen", screen_tier: str = "quick-screen",
                                          scenario_config: dict | None = None,
                                          temperature: float = 298.15, pressure: float = 1.0,
                                          target_duration_hours: float = 24.0, retention_fraction: float = 0.95) -> dict:
    """Propagate each profile's generated photon reactions through RMG once."""
    from .conditions import build_condition_spec
    from .stability import collect_rmg_evidence, resolve_stability_configuration

    if light_result.get("status") != "completed":
        raise ValueError("Computational light model must complete before kinetic propagation")
    run_dir = Path(run_dir)
    # A molecule can have candidate excited-state routes yet no admitted
    # reactive photon branch in any profile.  In that valid no-channel case no
    # per-profile RMG directory is created, so create the report directory
    # explicitly before retaining the all-no-channel result.
    run_dir.mkdir(parents=True, exist_ok=True)
    resolved_scenario, resources = resolve_stability_configuration(scenario, screen_tier)
    scenario_config = normalize_target_environment_species(
        light_result["smiles"], scenario_config or resolved_scenario,
    )
    condition = build_condition_spec(scenario_config, temperature_K=temperature, pressure_bar=pressure,
                                     target_duration_hours=target_duration_hours, retention_fraction=retention_fraction,
                                     light_condition="sunlight", light_model=light_result["model_contract"])
    profiles = {}
    for name, profile in light_result.get("profiles", {}).items():
        libraries = [Path(value) for value in profile["reaction_libraries"]]
        if not libraries:
            profiles[name] = {"status": "no_reactive_photon_channel", "estimated_time_to_retention_seconds": None}
            continue
        evidence = collect_rmg_evidence(
            light_result["smiles"], run_dir / name, rmg_env=rmg_env, temperature=temperature, pressure=pressure,
            rmg_walltime=resources["walltime"], rmg_max_processes=resources["max_processes"],
            rmg_max_iterations=resources["max_iterations"], rmg_max_edge_species=resources["max_edge_species"],
            scenario=scenario_config, conditions=condition, reaction_libraries=libraries,
        )
        profiles[name] = {"status": evidence["status"],
                          "estimated_time_to_retention_seconds": (evidence.get("kinetic_relevance") or {}).get("estimated_time_to_retention_seconds"),
                          "target_loss_fraction": (evidence.get("solver_profile") or {}).get("target_loss_fraction"),
                          "rmg_evidence": evidence}
    profile_statuses = {str(profile.get("status") or "") for profile in profiles.values()}
    complete_statuses = {"completed", "no_reactive_photon_channel"}
    if profile_statuses and profile_statuses == {"no_reactive_photon_channel"}:
        status = "completed_no_reactive_photon_channel"
        outcome = "no_modeled_sunlight_loss"
    elif profile_statuses and profile_statuses <= complete_statuses:
        status = "completed"
        outcome = "photon_profile_kinetics_completed"
    else:
        status = "completed_with_incomplete_evidence"
        outcome = "photon_profile_kinetics_incomplete"
    result = {"schema_version": 1, "kind": "computational_light_profile_kinetics", "condition_contract": condition.as_dict(),
              "status": status, "outcome": outcome, "profiles": profiles,
                                                "interpretation": "Low/central/high t95 values are consequences of versioned photon-branch priors and RMG propagation, not measured quantum-yield lifetimes. A no-channel result is not photostability evidence."}
    output = run_dir / "computational-light-kinetics.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {**result, "result_json": output}
