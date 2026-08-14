"""Command-line interface for STORCA."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Annotated, Optional

import typer

from .runs import create_run_directory, write_metadata
from .workflow import run_optimization_and_frequency
from src.orca_runner import _is_gnome_orca, find_xtb

app = typer.Typer(help="Reproducible ORCA calculation workflows.", no_args_is_help=True)


@app.command(name="direct-local-dft")
def direct_local_dft(
    run_dir: Annotated[Path, typer.Option(help="Existing six-representative spectrum run.")],
    nist_reference: Annotated[Optional[Path], typer.Option(help="Condition-compatible NIST/reference spectrum CSV.")] = None,
    baseline_spectrum: Annotated[Optional[Path], typer.Option(help="Existing isolated-cluster baseline CSV.")] = None,
    maximum_invocations: Annotated[int, typer.Option(min=1, help="Hard cap for additional ORCA gradient processes.")] = 24,
    cores: Annotated[int, typer.Option(min=1, help="Cores per ORCA gradient invocation.")] = 8,
    execute: Annotated[bool, typer.Option(help="Execute approved ORCA jobs; otherwise write only the fail-closed plan.")] = False,
) -> None:
    """Plan or execute the six-representative direct local-DFT gate."""
    from .direct_local_dft import (
        plan_six_representative_direct_local_dft,
        run_six_representative_direct_local_dft,
    )
    conformers_path = run_dir / "clusters" / "conformers.json"
    if not conformers_path.is_file():
        raise typer.BadParameter(f"Environment conformers not found: {conformers_path}")
    if execute:
        result = run_six_representative_direct_local_dft(
            run_dir, maximum_additional_orca_invocations=maximum_invocations,
            ncores=cores, nist_reference_csv=nist_reference,
            baseline_spectrum_csv=baseline_spectrum,
        )
    else:
        result = plan_six_representative_direct_local_dft(
            json.loads(conformers_path.read_text()),
            maximum_additional_orca_invocations=maximum_invocations,
        )
        path = run_dir / "direct-local-dft-plan.json"
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        result["artifact"] = str(path)
    typer.echo(f"Direct local-DFT status: {result['status']}")
    typer.echo(f"Artifact: {result.get('artifact', '')}")


@app.command(name="spectroscopy-gate-status")
def spectroscopy_gate_status(
    run_dir: Annotated[Path, typer.Option(help="Spectrum run containing retained gate artifacts.")],
) -> None:
    """Report which, if any, expensive spectroscopy stage is authorized next."""
    from .spectroscopy_gates import evaluate_spectroscopy_gate_sequence
    result = evaluate_spectroscopy_gate_sequence(run_dir)
    typer.echo(f"Gate state: {result['state']}")
    typer.echo(f"Authorized next stage: {result['authorized_next_expensive_stage'] or 'none'}")
    typer.echo(f"Artifact: {result['artifact']}")


@app.command(name="bulk-methanol-box")
def bulk_methanol_box(
    output_dir: Annotated[Path, typer.Option(help="Directory for the periodic-box artifacts.")],
    molecules: Annotated[int, typer.Option(min=1)] = 216,
    density: Annotated[float, typer.Option(min=0.001)] = 0.7866,
    seed: Annotated[int, typer.Option()] = 1,
) -> None:
    """Build a seeded periodic methanol starting box at the declared density."""
    from .bulk_embedding import BulkEmbeddingConfig, build_periodic_methanol_box
    result = build_periodic_methanol_box(
        output_dir,
        config=BulkEmbeddingConfig(molecule_count=molecules, density_g_cm3=density),
        seed=seed,
    )
    typer.echo(f"Periodic box: {result['artifact']}")


@app.command(name="sunlight-spectrum-fetch")
def sunlight_spectrum_fetch(
    output: Annotated[Path, typer.Option(help="Destination for normalized ASTM G173 AM1.5 global-tilt CSV.")] = Path("references/astm-g173-am15g.csv"),
) -> None:
    """Fetch and retain STORCA's default terrestrial AM1.5 reference spectrum."""
    from .sunlight import fetch_default_am15g
    try:
        path = fetch_default_am15g(output)
    except Exception as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Fetched AM1.5 global-tilt spectrum: {path}")


def _input_xyz(
    xyz_file: Optional[Path], smiles: Optional[str], run_dir: Path
) -> Path:
    if bool(xyz_file) == bool(smiles):
        raise typer.BadParameter("Provide exactly one of an XYZ file or --smiles.")
    if xyz_file:
        if not xyz_file.is_file():
            raise typer.BadParameter(f"XYZ file not found: {xyz_file}")
        target = run_dir / "input.xyz"
        shutil.copy2(xyz_file, target)
        return target

    try:
        from src.molecule_tools import sanitize_smiles, smiles_to_xyz
    except ImportError as error:
        raise typer.BadParameter("SMILES input requires RDKit. Install storca[chemistry].") from error
    canonical_smiles = sanitize_smiles(smiles or "")
    target = run_dir / "input.xyz"
    smiles_to_xyz(canonical_smiles, target)
    return target


@app.command()
def run(
    xyz_file: Annotated[Optional[Path], typer.Argument(help="Input XYZ geometry.")] = None,
    smiles: Annotated[Optional[str], typer.Option(help="SMILES input.")] = None,
    output_dir: Annotated[Path, typer.Option(help="Parent directory for calculation runs.")] = Path("runs"),
    name: Annotated[Optional[str], typer.Option(help="Name used for the run directory.")] = None,
    charge: Annotated[int, typer.Option(help="Molecular charge.")] = 0,
    multiplicity: Annotated[int, typer.Option(help="Spin multiplicity.")] = 1,
    cores: Annotated[int, typer.Option(min=1, help="ORCA CPU cores.")] = 1,
    skip_frequency: Annotated[bool, typer.Option(help="Run optimization only.")] = False,
) -> None:
    """Optimize a molecule, then calculate vibrational frequencies."""
    if multiplicity < 1:
        raise typer.BadParameter("Multiplicity must be at least 1.")
    source_name = name or (xyz_file.stem if xyz_file else "smiles")
    run_dir = create_run_directory(output_dir, source_name)
    try:
        from .methods import harmonic_method_profile

        input_xyz = _input_xyz(xyz_file, smiles, run_dir)
        write_metadata(
            run_dir,
            source_xyz=str(xyz_file) if xyz_file else None,
            smiles=smiles,
            harmonic_method_profile=harmonic_method_profile(),
            vibrational_treatment="harmonic" if not skip_frequency else "not_requested",
        )
        result = run_optimization_and_frequency(
            input_xyz,
            run_dir,
            charge=charge,
            multiplicity=multiplicity,
            ncores=cores,
            run_frequency=not skip_frequency,
        )
    except Exception:
        typer.echo(f"Run directory retained for diagnosis: {run_dir}", err=True)
        raise

    typer.echo(f"Run completed: {run_dir}")
    typer.echo(f"Optimized geometry: {result['optimized_xyz']}")
    if "frequency_check" in result:
        check = result["frequency_check"]
        typer.echo(f"True minimum: {check['IsMinimum']}")
        typer.echo(f"Imaginary frequencies: {check['NumImaginary']}")


@app.command()
def doctor() -> None:
    """Report availability of external programs and optional Python packages."""
    orca = os.environ.get("STORCA_ORCA_BIN") or shutil.which("orca")
    if orca and _is_gnome_orca(Path(orca)):
        typer.echo(f"WRONG   ORCA ({orca}; this is GNOME Orca, not quantum-chemistry ORCA)")
        orca = None
    else:
        typer.echo(f"{'OK' if orca else 'MISSING':7} ORCA{f' ({orca})' if orca else ''}")
    try:
        xtb = find_xtb()
    except RuntimeError:
        xtb = None
    checks = {
        "xTB": xtb,
        "Open Babel": shutil.which("obabel"),
        "RMG": shutil.which("rmg.py"),
    }
    for name, executable in checks.items():
        typer.echo(f"{'OK' if executable else 'MISSING':7} {name}{f' ({executable})' if executable else ''}")
    for name, module in {"RDKit": "rdkit", "Open Babel Python": "openbabel"}.items():
        found = importlib.util.find_spec(module) is not None
        typer.echo(f"{'OK' if found else 'MISSING':7} {name}")


@app.command()
def describe(
    smiles: Annotated[str, typer.Option(help="SMILES input.")],
    json_output: Annotated[Optional[Path], typer.Option(help="Optional machine-readable description JSON.")] = None,
) -> None:
    """Describe a molecule locally from RDKit; no network lookup is performed."""
    try:
        from .describe import describe_smiles, description_lines
        result = describe_smiles(smiles)
    except ImportError as error:
        raise typer.BadParameter("Molecule descriptions require RDKit. Install storca[chemistry].") from error
    except ValueError as error:
        raise typer.BadParameter(str(error)) from error
    for line in description_lines(result):
        typer.echo(line)
    typer.echo(f"Description: {result['description']}")
    if json_output:
        json_output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        typer.echo(f"JSON: {json_output}")


@app.command()
def structure(
    smiles: Annotated[str, typer.Option(help="SMILES input.")],
    output_dir: Annotated[Path, typer.Option(help="Parent directory for molecule records.")] = Path("runs"),
    name: Annotated[Optional[str], typer.Option(help="Name used for the record directory.")] = None,
) -> None:
    """Create local 2D/3D structure artifacts and a provenance-labelled record."""
    from .structure import build_structure_artifacts

    run_dir = create_run_directory(output_dir, name or "molecule")
    try:
        result = build_structure_artifacts(smiles, run_dir)
        write_metadata(
            run_dir, workflow="local_structure_record",
            smiles=result["description"]["canonical_smiles"],
            network_accessed=False, status="completed",
            molecule_record=result["record_json"],
        )
    except (ImportError, ValueError, RuntimeError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Molecule record: {result['record_json']}")
    typer.echo(f"2D depiction: {result['artifacts']['structure_png']}")
    typer.echo(f"Initial 3D geometry: {result['artifacts']['initial_xyz']}")


@app.command()
def enrich(
    smiles: Annotated[str, typer.Option(help="SMILES input.")],
    source: Annotated[str, typer.Option(help="External source; currently only pubchem.")] = "pubchem",
    json_output: Annotated[Path, typer.Option(help="Source-labelled enrichment JSON.")] = Path("enrichment.json"),
    synonym_limit: Annotated[int, typer.Option(min=0, max=100)] = 10,
    timeout: Annotated[float, typer.Option(min=1.0, max=60.0)] = 10.0,
) -> None:
    """Explicitly look up compound identifiers/properties from an external source."""
    if source.lower() != "pubchem":
        raise typer.BadParameter("Only source='pubchem' is currently supported")
    from .enrich import enrich_pubchem, write_enrichment

    try:
        result = enrich_pubchem(smiles, timeout_seconds=timeout, synonym_limit=synonym_limit)
    except (ImportError, ValueError, RuntimeError) as error:
        raise typer.BadParameter(str(error)) from error
    write_enrichment(json_output, result)
    data = result["pubchem"]
    typer.echo(f"PubChem CID: {data['cid']}")
    typer.echo(f"Name: {data['iupac_name'] or data['title'] or 'not supplied'}")
    typer.echo("Network source: PubChem PUG REST (not a safety assessment)")
    typer.echo(f"JSON: {json_output}")


@app.command()
def plausibility(
    smiles: Annotated[str, typer.Option(help="Proposed molecular SMILES.")],
    output_dir: Annotated[Path, typer.Option(help="Parent directory for the persistence dossier.")] = Path("runs"),
    name: Annotated[Optional[str], typer.Option(help="Name used for the dossier directory.")] = None,
    charge: Annotated[int, typer.Option()] = 0,
    multiplicity: Annotated[int, typer.Option(min=1)] = 1,
    temperature: Annotated[float, typer.Option(min=1.0)] = 298.15,
    pressure: Annotated[float, typer.Option(min=0.001)] = 1.0,
    atmosphere: Annotated[str, typer.Option()] = "air",
    phase: Annotated[str, typer.Option()] = "ordinary condensed-phase handling",
    persistence_target: Annotated[str, typer.Option()] = "at least 24 hours",
    with_rmg: Annotated[bool, typer.Option(help="Run a bounded RMG pathway screen and attach evidence.")] = False,
    rmg_env: Annotated[Optional[str], typer.Option(help="Conda RMG environment, e.g. rmg_env.")] = None,
    rmg_walltime: Annotated[str, typer.Option(help="RMG wall time as DD:HH:MM:SS.")] = "00:00:10:00",
    rmg_max_processes: Annotated[int, typer.Option(min=1)] = 1,
    rmg_max_iterations: Annotated[int, typer.Option(min=1)] = 100,
    rmg_max_edge_species: Annotated[int, typer.Option(min=1)] = 250,
) -> None:
    """Create an evidence-first persistence dossier; no optimizer-only verdict."""
    from .plausibility import attach_rmg_evidence, create_plausibility_dossier, write_plausibility_dossier

    run_dir = create_run_directory(output_dir, name or "plausibility")
    try:
        dossier = create_plausibility_dossier(
            smiles, charge=charge, multiplicity=multiplicity,
            conditions={"temperature_K": temperature, "pressure_bar": pressure,
                        "atmosphere": atmosphere, "phase": phase,
                        "persistence_target": persistence_target},
        )
        if with_rmg:
            from .stability import collect_rmg_evidence
            evidence = collect_rmg_evidence(
                dossier["species"]["canonical_smiles"], run_dir, rmg_env=rmg_env,
                temperature=temperature, pressure=pressure, rmg_walltime=rmg_walltime,
                rmg_max_processes=rmg_max_processes, rmg_max_iterations=rmg_max_iterations,
                rmg_max_edge_species=rmg_max_edge_species, requested_phase=phase,
            )
            attach_rmg_evidence(dossier, evidence)
    except (ImportError, ValueError) as error:
        raise typer.BadParameter(str(error)) from error
    report = write_plausibility_dossier(run_dir / "plausibility.json", dossier)
    write_metadata(run_dir, workflow="molecular_persistence_dossier", smiles=dossier["species"]["canonical_smiles"],
                   conditions=dossier["conditions"], persistence_category="insufficient_evidence", status="completed")
    typer.echo("Persistence category: insufficient_evidence")
    typer.echo("No ORCA evidence has been run yet; this dossier cannot establish persistence.")
    typer.echo(f"Dossier: {report}")


@app.command()
def spectrum(
    smiles: Annotated[str, typer.Option(help="SMILES input.")],
    output_dir: Annotated[Path, typer.Option(help="Parent directory for spectrum runs.")] = Path("runs"),
    name: Annotated[Optional[str], typer.Option(help="Name used for the run directory.")] = None,
    charge: Annotated[int, typer.Option(help="Molecular charge.")] = 0,
    multiplicity: Annotated[int, typer.Option(help="Spin multiplicity.")] = 1,
    cores: Annotated[int, typer.Option(min=1, help="ORCA CPU cores.")] = 1,
    temperature: Annotated[float, typer.Option(min=1.0, help="Boltzmann temperature in K.")] = 298.15,
    pressure: Annotated[Optional[float], typer.Option(min=0.000001, help="Sample pressure in bar when known.")] = None,
    conformer_engine: Annotated[str, typer.Option(help="Conformer engine: goat (recommended) or rdkit.")] = "goat",
    initial_conformers: Annotated[int, typer.Option(min=1, help="RDKit conformers to generate in fallback mode.")] = 20,
    max_conformers: Annotated[int, typer.Option(min=1, help="Maximum conformers sent to ORCA frequencies.")] = 10,
    goat_population: Annotated[float, typer.Option(min=0.01, max=1.0, help="Cumulative GOAT population to retain.")] = 0.95,
    scale_factor: Annotated[Optional[float], typer.Option(min=0.1, max=1.2, help="Harmonic scale-factor override; otherwise use --method-profile.")] = None,
    method_profile: Annotated[str, typer.Option(help="Named harmonic method/scale-factor profile.")] = "b3lyp-def2-svp",
    fwhm: Annotated[float, typer.Option(min=1.0, help="Gaussian band width (FWHM) in cm^-1.")] = 15.0,
    spectrum_model: Annotated[str, typer.Option(help="Spectrum model: raw, practical, or experimental FTIR.")] = "raw",
    spectrum_style: Annotated[str, typer.Option(help="Display: transmittance, absorbance, or relative.")] = "transmittance",
    max_absorbance: Annotated[float, typer.Option(min=0.01, help="Deepest plotted absorbance for transmittance display.")] = 1.0,
    phase: Annotated[str, typer.Option(help="Sample phase for experimental model: gas, solution, liquid, or solid.")] = "liquid",
    measurement: Annotated[str, typer.Option(help="Measurement for experimental model: auto, ATR, transmission, or gas-cell.")] = "auto",
    solvent: Annotated[Optional[str], typer.Option(help="Solvent identity, or 'neat' for a pure liquid.")] = None,
    composition: Annotated[Optional[str], typer.Option(help="JSON mapping of component identities to nonnegative fractions.")] = None,
    concentration: Annotated[Optional[float], typer.Option(min=0.000001, help="Analyte concentration in mol/L when known.")] = None,
    path_length: Annotated[Optional[float], typer.Option(min=0.000001, help="Transmission or gas-cell path length in mm.")] = None,
    atr_crystal: Annotated[Optional[str], typer.Option(help="ATR crystal material when known.")] = None,
    atr_incidence_angle: Annotated[Optional[float], typer.Option(min=0.1, max=89.9, help="ATR incidence angle in degrees.")] = None,
    sample_refractive_index: Annotated[Optional[float], typer.Option(min=0.000001, help="Sample refractive index for ATR rendering when known.")] = None,
    instrument_resolution: Annotated[float, typer.Option(min=0.1, help="Nominal FTIR resolution in cm^-1.")] = 4.0,
    apodization: Annotated[str, typer.Option(help="Instrument line shape: gaussian, triangular, or happ-genzel.")] = "happ-genzel",
    residual_fwhm: Annotated[Optional[float], typer.Option(min=0.1, help="Residual physical linewidth; defaults by phase.")] = None,
    fidelity: Annotated[str, typer.Option(help="Experimental fidelity: fast, auto, or balanced.")] = "auto",
    max_clusters: Annotated[int, typer.Option(min=1, help="Maximum sampled environments sent to ORCA.")] = 7,
    cluster_energy_window: Annotated[float, typer.Option(min=0.1, help="Deprecated legacy force-field screening setting; retained for CLI compatibility.")] = 5.0,
    max_orca_jobs: Annotated[int, typer.Option(min=1, help="Hard cap across monomer and cluster frequency jobs.")] = 12,
    local_mode_orca_invocations: Annotated[int, typer.Option(min=0, help="Additional explicit hard cap for ORCA local-mode gradient process invocations; zero disables fallback execution.")] = 0,
    dry_run: Annotated[bool, typer.Option(help="Write and print the adaptive calculation plan without running ORCA.")] = False,
) -> None:
    """Predict a conformer-weighted IR spectrum from a SMILES string."""
    if multiplicity < 1:
        raise typer.BadParameter("Multiplicity must be at least 1.")
    try:
        from src.molecule_tools import sanitize_smiles, smiles_to_conformers, smiles_to_xyz
        from .methods import harmonic_method_profile, resolve_scale_factor
        from .spectrum import run_goat_search, run_ir_spectrum
    except ImportError as error:
        raise typer.BadParameter("Spectrum calculations require RDKit. Install storca[chemistry].") from error

    run_dir = create_run_directory(output_dir, name or "spectrum")
    try:
        profile = harmonic_method_profile(method_profile)
        resolved_scale_factor, scale_factor_source = resolve_scale_factor(method_profile, scale_factor)
        canonical_smiles = sanitize_smiles(smiles)
        try:
            composition_map = json.loads(composition) if composition else {}
        except json.JSONDecodeError as error:
            raise typer.BadParameter("--composition must be a JSON object") from error
        if not isinstance(composition_map, dict):
            raise typer.BadParameter("--composition must be a JSON object")
        if not composition_map and phase.lower() == "liquid":
            composition_map = {canonical_smiles: 1.0}
        resolved_solvent = solvent or ("neat" if phase.lower() == "liquid" else None)
        experimental_condition_details = {
            "pressure_bar": pressure,
            "composition": composition_map,
            "solvent": resolved_solvent,
            "concentration_mol_L": concentration,
            "path_length_mm": path_length,
            "atr_crystal": atr_crystal,
            "atr_incidence_angle_degrees": atr_incidence_angle,
            "sample_refractive_index": sample_refractive_index,
        }
        engine = conformer_engine.lower()
        if engine not in {"goat", "rdkit"}:
            raise typer.BadParameter("Conformer engine must be 'goat' or 'rdkit'.")
        spectrum_model = spectrum_model.lower()
        if spectrum_model not in {"raw", "practical", "experimental"}:
            raise typer.BadParameter("spectrum_model must be 'raw', 'practical', or 'experimental'.")
        if spectrum_model == "experimental":
            from .experimental_ir import resolve_experimental_profile
            from .ir_contracts import build_experimental_condition
            try:
                measurement_profile = resolve_experimental_profile(
                    phase=phase, measurement=measurement,
                    instrument_resolution_cm_1=instrument_resolution,
                    apodization=apodization, residual_fwhm_cm_1=residual_fwhm,
                )
                build_experimental_condition(
                    phase=measurement_profile.phase,
                    measurement=measurement_profile.measurement,
                    temperature_K=temperature,
                    resolution_cm_1=measurement_profile.instrument_resolution_cm_1,
                    apodization=measurement_profile.apodization,
                    **experimental_condition_details,
                )
            except ValueError as error:
                raise typer.BadParameter(str(error)) from error
        fidelity = fidelity.lower()
        if fidelity not in {"fast", "auto", "balanced"}:
            raise typer.BadParameter("fidelity must be 'fast', 'auto', or 'balanced'; VPT2 and dynamics are not implemented yet.")
        from .cluster_ir import dimer_sampling_plan
        from .environment_selection import allocate_orca_budget, select_xtb_environment_representatives
        from .environment_convergence import EnvironmentConvergenceConfig, convergence_schedule
        from .xtb_sampling import sample_xtb_dimer_environments, xtb_sampling_defaults
        from .xtb_trajectory import xtb_trajectory_defaults
        cluster_plan = dimer_sampling_plan(
            canonical_smiles, phase=phase.lower(), charge=charge, multiplicity=multiplicity,
        )
        requested_clusters = (
            spectrum_model == "experimental" and fidelity in {"auto", "balanced"}
            and cluster_plan["eligible"]
        )
        orca_allocation = allocate_orca_budget(
            fidelity=fidelity, max_conformers=max_conformers,
            max_environments=max_clusters, max_orca_jobs=max_orca_jobs,
            environment_eligible=requested_clusters,
        )
        effective_max_conformers = orca_allocation["maximum_monomer_jobs"]
        xtb_environment_plan = xtb_sampling_defaults(fidelity)
        plan = {
            "schema_version": 1,
            "smiles": canonical_smiles,
            "spectrum_model": spectrum_model,
            "fidelity": fidelity,
            "experimental_condition_details": experimental_condition_details,
            "maximum_monomer_jobs": effective_max_conformers,
            "maximum_cluster_jobs": orca_allocation["reserved_environment_jobs"],
            "maximum_total_orca_jobs": max_orca_jobs,
            "maximum_local_mode_orca_invocations": local_mode_orca_invocations,
            "orca_allocation": orca_allocation,
            "dimer_sampling": cluster_plan,
            "dimer_sampling_requested": requested_clusters,
            "xtb_environment_sampling": {
                **xtb_environment_plan,
                "requested": requested_clusters,
                "xtb_sampling_orca_jobs_consumed": 0,
                "representative_orca_jobs_reserved": orca_allocation["reserved_environment_jobs"],
                "representative_selection": "mode_class_coverage_then_frequency_geometry_diversity",
                "include_trimers": bool(cluster_plan.get("trimer_eligible")),
            },
            "xtb_trajectory_sampling": {
                **xtb_trajectory_defaults(fidelity),
                "temperature_K": temperature,
                "requested": requested_clusters,
                "seed_source": "diverse_restrained_xtb_optimized_strata",
                "snapshot_use": "decorrelated_occupancy_weighted_environment_ensemble",
            },
            "adaptive_environment_convergence": {
                "requested": requested_clusters,
                "configuration": asdict(EnvironmentConvergenceConfig()),
                "cumulative_batch_schedule": convergence_schedule(
                    fidelity, orca_allocation["reserved_environment_jobs"],
                ),
            },
            "unimplemented_tiers": ["microsolvation", "VPT2", "unrestrained_bulk_dynamics", "AIMD"],
        }
        plan_path = run_dir / "spectrum-plan.json"
        plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
        if dry_run:
            write_metadata(run_dir, workflow="ir_spectrum", smiles=canonical_smiles,
                           calculation_plan=str(plan_path), status="planned_not_run")
            typer.echo(json.dumps(plan, indent=2, sort_keys=True))
            typer.echo(f"Plan: {plan_path}")
            return
        seed_xyz = run_dir / "seed.xyz"
        smiles_to_xyz(canonical_smiles, seed_xyz)
        if engine == "goat":
            xyz_paths = run_goat_search(
                seed_xyz, run_dir, charge=charge, multiplicity=multiplicity,
                ncores=cores, population_cutoff=goat_population,
                max_conformers=effective_max_conformers, progress=typer.echo,
            )
        else:
            initial_dir = run_dir / "initial-conformers"
            xyz_paths = smiles_to_conformers(canonical_smiles, initial_dir, n_confs=initial_conformers, max_confs=effective_max_conformers)
            write_metadata(run_dir, conformer_engine="rdkit", smiles=canonical_smiles, initial_conformers=len(xyz_paths))
        cluster_paths: list[Path] = []
        xtb_sampling_records: list[dict] = []
        xtb_frequency_records: list[dict] = []
        xtb_sampling_manifest: Path | None = None
        environment_selection_error: str | None = None
        if requested_clusters:
            try:
                xtb_sampling_records, xtb_sampling_manifest = sample_xtb_dimer_environments(
                    canonical_smiles, run_dir, fidelity=fidelity,
                    monomer_xyz=xyz_paths[0], charge=charge, multiplicity=multiplicity,
                    ncores=cores, progress=typer.echo,
                    include_trimers=bool(cluster_plan.get("trimer_eligible")),
                )
                retained_xtb = sum(
                    record.get("sampling_status") == "retained"
                    for record in xtb_sampling_records
                )
                typer.echo(
                    f"Retained {retained_xtb}/{len(xtb_sampling_records)} restrained xTB environments: "
                    f"{xtb_sampling_manifest}"
                )
                try:
                    from .xtb_trajectory import sample_restrained_xtb_trajectories
                    trajectory_records, trajectory_manifest = sample_restrained_xtb_trajectories(
                        xtb_sampling_records, run_dir, fidelity=fidelity, charge=charge,
                        multiplicity=multiplicity, ncores=cores, temperature_K=temperature,
                        progress=typer.echo,
                    )
                    xtb_sampling_records = trajectory_records
                    typer.echo(
                        f"Retained {len(trajectory_records)} decorrelated restrained-xTB "
                        f"trajectory snapshots: {trajectory_manifest}"
                    )
                except Exception as trajectory_error:
                    write_metadata(
                        run_dir, xtb_trajectory_sampling_status="failed_static_fallback",
                        xtb_trajectory_sampling_error=str(trajectory_error),
                    )
                    typer.echo(
                        "Restrained xTB trajectory sampling failed; using the retained static "
                        f"coverage ensemble ({trajectory_error})", err=True,
                    )
                from .xtb_frequencies import sample_xtb_snapshot_frequencies
                xtb_frequency_records, frequency_manifest = sample_xtb_snapshot_frequencies(
                    xtb_sampling_records, run_dir, charge=charge,
                    multiplicity=multiplicity, ncores=cores, progress=typer.echo,
                )
                typer.echo(
                    f"Completed unrestrained xTB snapshot frequencies for "
                    f"{sum(item.get('frequency_status') == 'completed' for item in xtb_frequency_records)}"
                    f"/{len(xtb_frequency_records)} retained environments: {frequency_manifest}"
                )
            except Exception as error:
                write_metadata(
                    run_dir, xtb_environment_sampling_status="failed",
                    xtb_environment_sampling_error=str(error),
                )
                typer.echo(
                    f"Restrained xTB environment sampling failed; retaining the monomer spectrum ({error})",
                    err=True,
                )
            if xtb_sampling_records:
                try:
                    cluster_paths, cluster_manifest = select_xtb_environment_representatives(
                        xtb_sampling_records, run_dir,
                        representative_count=orca_allocation["reserved_environment_jobs"],
                        frequency_records=xtb_frequency_records,
                    )
                    typer.echo(
                        f"Selected {len(cluster_paths)} diverse xTB environment representatives: "
                        f"{cluster_manifest}"
                    )
                    from .environment_refinement import refine_selected_orca_environments
                    cluster_paths, refinement_manifest = refine_selected_orca_environments(
                        run_dir, charge=0, multiplicity=1, ncores=cores,
                        method_keywords=profile["orca_keywords"], progress=typer.echo,
                    )
                    typer.echo(
                        f"Completed environment-preserving DFT refinement and gradient gates: "
                        f"{refinement_manifest}"
                    )
                except Exception as error:
                    environment_selection_error = str(error)
                    write_metadata(
                        run_dir, environment_selection_status="failed",
                        environment_selection_error=environment_selection_error,
                    )
                    typer.echo(
                        f"xTB diversity selection failed; retaining the monomer spectrum ({error})",
                        err=True,
                    )
            if xtb_sampling_manifest and xtb_sampling_manifest.is_file():
                environment_report = json.loads(xtb_sampling_manifest.read_text())
                environment_report["orca_budget"] = orca_allocation
                xtb_sampling_manifest.write_text(
                    json.dumps(environment_report, indent=2, sort_keys=True) + "\n"
                )
        write_metadata(
            run_dir, workflow="ir_spectrum", smiles=canonical_smiles, charge=charge,
            multiplicity=multiplicity, ncores=cores, temperature=temperature,
            scale_factor=resolved_scale_factor, scale_factor_source=scale_factor_source,
            harmonic_method_profile=profile, vibrational_treatment="harmonic",
            spectrum_model=spectrum_model,
            spectrum_model_description=(
                "Calculated conformer ensemble plus condition-aware FTIR measurement transfer"
                if spectrum_model == "experimental"
                else "Boltzmann-weighted conformer ensemble with Gaussian broadening"
            ),
            fwhm_cm_1=fwhm, spectrum_style=spectrum_style,
            max_absorbance=max_absorbance, phase=phase, measurement=measurement,
            instrument_resolution_cm_1=instrument_resolution, apodization=apodization,
            residual_fwhm_cm_1=residual_fwhm, status="running",
            experimental_condition_details=experimental_condition_details,
            fidelity=fidelity, max_orca_jobs=max_orca_jobs,
            local_mode_orca_invocations=local_mode_orca_invocations,
            orca_allocation=orca_allocation,
            calculation_plan=str(plan_path), cluster_candidates=len(cluster_paths),
            xtb_environment_sampling_status=(
                "completed" if xtb_sampling_manifest is not None
                else "not_requested_or_failed"
            ),
            xtb_environment_candidates=len(xtb_sampling_records),
            xtb_environment_candidates_retained=sum(
                record.get("sampling_status") == "retained"
                for record in xtb_sampling_records
            ),
            xtb_environment_sampling_manifest=str(xtb_sampling_manifest or ""),
            environment_selection_status=(
                "completed" if cluster_paths else
                "failed" if environment_selection_error else
                "not_requested_or_unavailable"
            ),
            environment_selection_error=environment_selection_error,
        )
        result = run_ir_spectrum(
            xyz_paths,
            run_dir,
            charge=charge,
            multiplicity=multiplicity,
            ncores=cores,
            temperature=temperature,
            scale_factor=resolved_scale_factor,
            fwhm=fwhm,
            spectrum_style=spectrum_style,
            max_absorbance=max_absorbance,
            spectrum_model=spectrum_model,
            practical_smiles=canonical_smiles,
            phase=phase,
            measurement=measurement,
            instrument_resolution=instrument_resolution,
            apodization=apodization,
            residual_fwhm=residual_fwhm,
            experimental_condition_details=experimental_condition_details,
            method_keywords=profile["orca_keywords"],
            progress=typer.echo,
        )
        if cluster_paths:
            from .spectrum import assemble_self_dimer_environment
            from .environment_convergence import run_adaptive_environment_convergence
            cluster_dir = run_dir / "clusters"
            try:
                typer.echo(f"Adaptively calculating up to {len(cluster_paths)} dimer/trimer environments")
                def record_transfer_batch(endpoint: int, _batch_result: dict) -> list[Path]:
                    from .environment_acquisition import (record_acquisition_batch,
                                                           reprioritize_pending_representatives)
                    from .frequency_transfer import build_frequency_transfer_artifacts
                    transfer = build_frequency_transfer_artifacts(run_dir)
                    record_acquisition_batch(run_dir, endpoint, transfer)
                    return reprioritize_pending_representatives(run_dir, endpoint)

                cluster_result, convergence_report = run_adaptive_environment_convergence(
                    cluster_paths, cluster_dir, run_dir, fidelity=fidelity,
                    scale_factor=resolved_scale_factor,
                    run_batch=lambda paths: run_ir_spectrum(
                        paths, cluster_dir, charge=0, multiplicity=1,
                        ncores=cores, temperature=temperature,
                        scale_factor=resolved_scale_factor, fwhm=fwhm,
                        spectrum_style=spectrum_style, max_absorbance=max_absorbance,
                        spectrum_model="raw", method_keywords=profile["orca_keywords"],
                        geometry_role="environment_snapshot", progress=typer.echo,
                    ),
                    progress=typer.echo,
                    after_batch=record_transfer_batch,
                )
                result = assemble_self_dimer_environment(
                    result, cluster_result["conformers"], run_dir,
                    scale_factor=resolved_scale_factor, fwhm=fwhm,
                    spectrum_style=spectrum_style, max_absorbance=max_absorbance,
                    phase=phase, measurement=measurement,
                    instrument_resolution=instrument_resolution,
                    apodization=apodization, residual_fwhm=residual_fwhm,
                    environment_convergence=convergence_report,
                    experimental_condition_details=experimental_condition_details,
                )
                if local_mode_orca_invocations > 0:
                    from .environment_local_modes import run_environment_local_mode_fallbacks
                    local_modes = run_environment_local_mode_fallbacks(
                        run_dir,
                        maximum_orca_invocations=local_mode_orca_invocations,
                        ncores=cores, method_keywords=profile["orca_keywords"],
                        progress=typer.echo,
                    )
                    result["environment_local_modes"] = local_modes.get("artifact")
                    write_metadata(
                        run_dir,
                        environment_local_mode_status=local_modes.get("status"),
                        environment_local_mode_artifact=local_modes.get("artifact"),
                    )
                from .frequency_transfer import build_frequency_transfer_artifacts
                transfer_artifacts = build_frequency_transfer_artifacts(run_dir)
                result.update(transfer_artifacts)
                write_metadata(
                    run_dir,
                    frequency_transfer_status=transfer_artifacts.get("status"),
                    frequency_transfer_use_status=transfer_artifacts.get(
                        "frequency_transfer_use_status"
                    ),
                    frequency_transfer_display_basis=transfer_artifacts.get("display_basis"),
                    frequency_transfer_artifact=str(
                        transfer_artifacts.get("frequency_transfer") or ""
                    ),
                    frequency_transfer_validation_artifact=str(
                        transfer_artifacts.get("frequency_transfer_validation") or ""
                    ),
                )
            except Exception as error:
                write_metadata(run_dir, cluster_calculation_status="failed_fallback_to_monomer",
                               cluster_calculation_error=str(error), status="completed")
                typer.echo(f"Dimer calculations incomplete; retaining monomer spectrum ({error})", err=True)
    except Exception:
        typer.echo(f"Run directory retained for diagnosis: {run_dir}", err=True)
        raise
    typer.echo(f"Spectrum run completed: {run_dir}")
    typer.echo(f"Spectrum CSV: {result['spectrum_csv']}")
    typer.echo(f"Spectrum PNG: {result['spectrum_png']}")
    if "raw_spectrum_csv" in result:
        typer.echo(f"Raw calculated spectrum CSV: {result['raw_spectrum_csv']}")
    if "ensemble_spectrum_csv" in result:
        typer.echo(f"Calculated ensemble spectrum CSV: {result['ensemble_spectrum_csv']}")
    if "intrinsic_spectrum_csv" in result:
        typer.echo(f"Intrinsic calculated spectrum CSV: {result['intrinsic_spectrum_csv']}")


@app.command(name="spectrum-extend-environments")
def spectrum_extend_environments(
    run_dir: Annotated[Path, typer.Argument(help="Existing experimental spectrum run.")],
    maximum_candidates: Annotated[int, typer.Option(min=41, max=500, help="Hard cumulative xTB candidate cap.")] = 100,
    batch_size: Annotated[int, typer.Option(min=5, max=100, help="New candidates per balanced acquisition round.")] = 20,
    cores: Annotated[int, typer.Option(min=1, help="xTB CPU cores.")] = 1,
) -> None:
    """Extend a retained xTB environment ensemble in convergence-controlled rounds."""
    from .adaptive_xtb_extension import extend_xtb_environment_ensemble

    run_dir = Path(run_dir)
    metadata_path = run_dir / "metadata.json"
    sampling_path = run_dir / "environment-sampling.json"
    if not metadata_path.is_file() or not sampling_path.is_file():
        raise typer.BadParameter("Run lacks retained experimental environment artifacts")
    metadata = json.loads(metadata_path.read_text())
    sampling = json.loads(sampling_path.read_text()).get("sampling") or {}
    configuration = sampling.get("configuration") or {}
    monomer_xyz = Path(str(configuration.get("monomer_xyz") or ""))
    smiles = metadata.get("smiles")
    if not smiles or not monomer_xyz.is_file():
        raise typer.BadParameter("Run lacks its retained SMILES or monomer sampling geometry")
    result = extend_xtb_environment_ensemble(
        run_dir, smiles=str(smiles), monomer_xyz=monomer_xyz,
        maximum_candidates=maximum_candidates, batch_size=batch_size,
        charge=int(metadata.get("charge", 0)),
        multiplicity=int(metadata.get("multiplicity", 1)),
        ncores=cores,
        include_trimers=bool(configuration.get("cluster_model") == "dimers_plus_trimers"),
        progress=typer.echo,
    )
    write_metadata(
        run_dir, xtb_extension_status=result["status"],
        xtb_extension_artifact=result["artifact"],
        xtb_environment_candidates=result["final_candidates"],
        xtb_extension_numerically_stable=result.get(
            "numerically_stable_across_acquisition_rounds"
        ),
        xtb_extension_bootstrap_precision_pass=result.get(
            "bootstrap_precision_pass"
        ),
    )
    typer.echo(f"xTB extension status: {result['status']}")
    typer.echo(f"Candidates evaluated: {result['final_candidates']}")
    typer.echo(f"Artifact: {result['artifact']}")


@app.command(name="spectrum-analyze")
def spectrum_analyze(
    calculated: Annotated[Path, typer.Argument(help="Calculated spectrum CSV.")],
    reference: Annotated[Path, typer.Argument(help="Experimental/reference spectrum CSV.")],
    output: Annotated[Path, typer.Option(help="Peak-level JSON report.")] = Path("spectrum-analysis.json"),
    overlay: Annotated[Optional[Path], typer.Option(help="Optional annotated PNG overlay.")] = None,
    min_prominence: Annotated[float, typer.Option(min=0.001, max=0.99, help="Minimum peak strength as a fraction of each trace's maximum.")] = 0.05,
    minimum_separation: Annotated[float, typer.Option(min=1.0, help="Minimum separation between reported peaks in cm^-1.")] = 12.0,
    match_window: Annotated[float, typer.Option(min=1.0, help="Maximum calculated/reference peak difference in cm^-1.")] = 40.0,
) -> None:
    """Match prominent IR bands without fitting or shifting either spectrum."""
    from .analysis import analyze_ir_spectra, write_peak_analysis, write_peak_analysis_plot

    report = analyze_ir_spectra(
        calculated, reference, min_prominence=min_prominence,
        minimum_separation_cm_1=minimum_separation, match_window_cm_1=match_window,
    )
    write_peak_analysis(output, report)
    if overlay:
        write_peak_analysis_plot(overlay, calculated=calculated, reference=reference, report=report)
        typer.echo(f"Overlay: {overlay}")
    typer.echo(f"Matched bands: {len(report['matches'])}")
    if report["mean_absolute_position_error_cm-1"] is not None:
        typer.echo(f"Mean absolute position error: {report['mean_absolute_position_error_cm-1']:.2f} cm^-1")
    typer.echo(f"Unmatched reference bands: {len(report['unmatched_reference'])}")
    typer.echo(f"Unmatched calculated bands: {len(report['unmatched_calculated'])}")
    typer.echo(f"Report: {output}")


@app.command(name="spectrum-resume")
def spectrum_resume(
    run_dir: Annotated[Path, typer.Argument(help="Interrupted spectrum run directory.")],
    charge: Annotated[Optional[int], typer.Option(help="Override retained charge.")] = None,
    multiplicity: Annotated[Optional[int], typer.Option(min=1, help="Override retained multiplicity.")] = None,
    cores: Annotated[Optional[int], typer.Option(min=1, help="Override retained ORCA core count.")] = None,
    temperature: Annotated[Optional[float], typer.Option(min=1.0, help="Override retained temperature.")] = None,
    scale_factor: Annotated[Optional[float], typer.Option(min=0.1, max=1.2, help="Override retained scale factor.")] = None,
    fwhm: Annotated[Optional[float], typer.Option(min=1.0, help="Override retained FWHM.")] = None,
    spectrum_style: Annotated[Optional[str], typer.Option(help="Override retained display style.")] = None,
    max_absorbance: Annotated[Optional[float], typer.Option(min=0.01, help="Override retained display absorbance.")] = None,
    phase: Annotated[Optional[str], typer.Option(help="Override retained sample phase.")] = None,
    measurement: Annotated[Optional[str], typer.Option(help="Override retained measurement geometry.")] = None,
    instrument_resolution: Annotated[Optional[float], typer.Option(min=0.1, help="Override retained FTIR resolution.")] = None,
    apodization: Annotated[Optional[str], typer.Option(help="Override retained apodization.")] = None,
    residual_fwhm: Annotated[Optional[float], typer.Option(min=0.1, help="Override retained residual linewidth.")] = None,
    local_mode_orca_invocations: Annotated[Optional[int], typer.Option(min=0, help="Explicit ORCA process-invocation allowance for local-mode fallbacks.")] = None,
) -> None:
    """Continue incomplete ORCA jobs using the original run settings by default."""
    from .spectrum import resume_ir_spectrum
    result = resume_ir_spectrum(
        run_dir, charge=charge, multiplicity=multiplicity, ncores=cores,
        temperature=temperature, scale_factor=scale_factor, fwhm=fwhm,
        spectrum_style=spectrum_style, max_absorbance=max_absorbance,
        phase=phase, measurement=measurement, instrument_resolution=instrument_resolution,
        apodization=apodization, residual_fwhm=residual_fwhm,
        local_mode_orca_invocations=local_mode_orca_invocations,
        progress=typer.echo,
    )
    typer.echo(f"Spectrum CSV: {result['spectrum_csv']}")
    typer.echo(f"Spectrum PNG: {result['spectrum_png']}")


@app.command()
def benchmark(
    predicted: Annotated[Path, typer.Argument(help="Predicted spectrum CSV.")],
    reference: Annotated[Path, typer.Argument(help="Reference spectrum CSV with matching intensity convention.")],
    output: Annotated[Path, typer.Option(help="JSON comparison report.")] = Path("benchmark.json"),
    phase: Annotated[str, typer.Option(help="Reference sample phase.")] = "unspecified",
    solvent: Annotated[str, typer.Option(help="Reference solvent or neat sample.")] = "unspecified",
    temperature: Annotated[Optional[float], typer.Option(help="Reference temperature in K.")] = None,
    measurement: Annotated[str, typer.Option(help="ATR, transmission, gas-phase, etc.")] = "unspecified",
    shift_window: Annotated[float, typer.Option(min=0.0, help="Search +/- this global band shift in cm^-1 (diagnostic only). ")] = 0.0,
    shift_step: Annotated[float, typer.Option(min=0.1, help="Global-shift search increment in cm^-1.")] = 1.0,
    overlay: Annotated[Optional[Path], typer.Option(help="Optional PNG overlay; calculated trace uses the diagnostic best shift.")] = None,
) -> None:
    """Compare a predicted spectrum CSV against a numerical reference spectrum."""
    from .benchmark import compare_spectra, write_benchmark_result, write_comparison_plot
    metrics = compare_spectra(predicted, reference, shift_window=shift_window, shift_step=shift_step)
    report = write_benchmark_result(output, predicted=predicted, reference=reference, metrics=metrics, conditions={"phase": phase, "solvent": solvent, "temperature_K": temperature, "measurement": measurement})
    correlation = metrics["pearson_correlation"]
    typer.echo(f"Correlation: {correlation:.4f}" if correlation is not None else "Correlation: undefined (flat trace)")
    typer.echo(f"RMSE: {metrics['rmse']:.4f}")
    if shift_window:
        aligned = metrics["shift_aligned"]
        typer.echo(f"Best diagnostic shift: {metrics['best_global_shift_cm-1']:+.1f} cm^-1")
        aligned_correlation = aligned["pearson_correlation"]
        typer.echo(f"Shift-aligned correlation: {aligned_correlation:.4f}"
                   if aligned_correlation is not None else "Shift-aligned correlation: undefined (flat trace)")
    if overlay:
        plot = write_comparison_plot(overlay, predicted=predicted, reference=reference, best_shift_cm_1=metrics["best_global_shift_cm-1"])
        typer.echo(f"Overlay: {plot}")
    typer.echo(f"Report: {report}")


@app.command(name="ir-benchmark")
def ir_benchmark(
    manifest: Annotated[Path, typer.Argument(help="IR benchmark manifest JSON.")],
    profile: Annotated[str, typer.Option(help="Method profile to evaluate.")] = "b3lyp-def2-svp",
    output: Annotated[Path, typer.Option(help="Aggregate JSON report.")] = Path("ir-benchmark-report.json"),
) -> None:
    """Aggregate existing raw spectra for one profile; no ORCA jobs are run."""
    from .ir_benchmark import evaluate_ir_manifest, write_ir_benchmark

    try:
        report = evaluate_ir_manifest(manifest, profile)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise typer.BadParameter(str(error)) from error
    write_ir_benchmark(output, report)
    aggregate = report["aggregate"]
    typer.echo(f"Evaluated entries: {report['evaluated_entries']}")
    typer.echo(f"Missing entries: {len(report['missing_entries'])}")
    if aggregate["mean_pearson_correlation"] is not None:
        typer.echo(f"Mean correlation: {aggregate['mean_pearson_correlation']:.4f}")
    if aggregate["mean_rmse"] is not None:
        typer.echo(f"Mean RMSE: {aggregate['mean_rmse']:.4f}")
    typer.echo(f"Report: {output}")


@app.command(name="plausibility-benchmark")
def plausibility_benchmark(
    manifest: Annotated[Path, typer.Argument(help="Plausibility benchmark manifest JSON.")],
    results_dir: Annotated[Path, typer.Argument(help="Directory containing one <benchmark-id>.json dossier per accepted case.")],
    output: Annotated[Path, typer.Option(help="Benchmark report JSON.")] = Path("plausibility-benchmark-report.json"),
) -> None:
    """Score persistence dossiers against accepted condition-matched references."""
    from .plausibility_benchmark import evaluate_plausibility_manifest, write_plausibility_benchmark

    try:
        report = evaluate_plausibility_manifest(manifest, results_dir)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise typer.BadParameter(str(error)) from error
    write_plausibility_benchmark(output, report)
    aggregate = report["aggregate"]
    typer.echo(f"Accepted entries: {report['accepted_entries']}")
    typer.echo(f"Evaluated entries: {report['evaluated_entries']}")
    typer.echo(f"False reassurance count: {aggregate['false_reassurance_count']}")
    typer.echo(f"Report: {output}")


@app.command(name="calibrate-ir-scale")
def calibrate_ir_scale(
    manifest: Annotated[Path, typer.Argument(help="IR benchmark manifest JSON.")],
    training_ids: Annotated[str, typer.Option(help="Comma-separated benchmark IDs used to choose the scale.")],
    holdout_ids: Annotated[str, typer.Option(help="Comma-separated benchmark IDs reserved for evaluation.")],
    profile: Annotated[str, typer.Option(help="Existing ORCA method profile to recalibrate.")] = "b3lyp-def2-svp",
    scale_start: Annotated[float, typer.Option(min=0.1, max=1.2)] = 0.94,
    scale_stop: Annotated[float, typer.Option(min=0.1, max=1.2)] = 1.00,
    scale_step: Annotated[float, typer.Option(min=0.0001, max=0.1)] = 0.002,
    fwhm: Annotated[float, typer.Option(min=1.0, help="Fixed FWHM during scale selection in cm^-1.")] = 15.0,
    temperature: Annotated[float, typer.Option(min=1.0)] = 298.15,
    output: Annotated[Path, typer.Option(help="Calibration report JSON.")] = Path("ir-scale-calibration.json"),
) -> None:
    """Select a global harmonic scale on training IDs and score untouched holdouts."""
    from .calibration import calibrate_harmonic_scale, write_calibration

    train = [item.strip() for item in training_ids.split(",") if item.strip()]
    holdout = [item.strip() for item in holdout_ids.split(",") if item.strip()]
    try:
        report = calibrate_harmonic_scale(
            manifest, profile=profile, training_ids=train, holdout_ids=holdout,
            scale_start=scale_start, scale_stop=scale_stop, scale_step=scale_step,
            fwhm=fwhm, temperature=temperature,
        )
    except (OSError, ValueError, json.JSONDecodeError, RuntimeError) as error:
        raise typer.BadParameter(str(error)) from error
    write_calibration(output, report)
    selection = report["selection"]
    holdout = report["holdout"]["aggregate"]
    typer.echo(f"Selected scale factor: {selection['selected_scale_factor']:.6f}")
    typer.echo(f"Training mean position error: {selection['training_metrics']['mean_absolute_position_error_cm-1']:.2f} cm^-1")
    if holdout["mean_absolute_position_error_cm-1"] is not None:
        typer.echo(f"Holdout mean position error: {holdout['mean_absolute_position_error_cm-1']:.2f} cm^-1")
    typer.echo("No ORCA jobs were run and no existing spectra were overwritten.")
    typer.echo(f"Report: {output}")


@app.command(name="digitize-ir")
def digitize_ir(
    image: Annotated[Path, typer.Argument(help="IR plot image to trace.")],
    output: Annotated[Path, typer.Option(help="Output two-column CSV.")],
    overlay: Annotated[Optional[Path], typer.Option(help="Optional QA overlay PNG/JPG.")] = None,
    left: Annotated[int, typer.Option(help="Plot-area left pixel.")] = 29,
    right: Annotated[int, typer.Option(help="Plot-area right pixel.")] = 713,
    top: Annotated[int, typer.Option(help="Plot-area top pixel.")] = 96,
    bottom: Annotated[int, typer.Option(help="Plot-area bottom pixel.")] = 417,
    x_left: Annotated[float, typer.Option(help="Wavenumber at left axis edge.")] = 4000.0,
    x_right: Annotated[float, typer.Option(help="Wavenumber at right axis edge.")] = 400.0,
    axis: Annotated[str, typer.Option(help="Axis calibration: sdbs (non-linear printed ticks) or linear.")] = "sdbs",
) -> None:
    """Digitize a clean monochrome transmittance plot with explicit calibration."""
    from .digitize import PlotCalibration, digitize_transmittance, write_digitized_csv, write_trace_overlay
    if axis.lower() == "sdbs":
        from .digitize import SDBS_IR_CALIBRATION
        calibration = SDBS_IR_CALIBRATION
    elif axis.lower() == "linear":
        calibration = PlotCalibration(left, right, top, bottom, x_left, x_right)
    else:
        raise typer.BadParameter("axis must be 'sdbs' or 'linear'.")
    wavenumbers, transmittance, pixels = digitize_transmittance(image, calibration)
    write_digitized_csv(output, wavenumbers, transmittance)
    if overlay:
        write_trace_overlay(image, overlay, calibration, pixels)
        typer.echo(f"Overlay: {overlay}")
    typer.echo(f"Digitized {len(wavenumbers)} points: {output}")


@app.command(name="spectrum-finalize")
def spectrum_finalize(
    run_dir: Annotated[Path, typer.Argument(help="Existing spectrum run directory.")],
    temperature: Annotated[Optional[float], typer.Option(min=1.0, help="Override retained Boltzmann temperature.")] = None,
    scale_factor: Annotated[Optional[float], typer.Option(min=0.1, max=1.2, help="Override retained harmonic scale factor.")] = None,
    fwhm: Annotated[Optional[float], typer.Option(min=1.0, help="Override retained Gaussian FWHM.")] = None,
    spectrum_model: Annotated[Optional[str], typer.Option(help="Override retained model: raw, practical, or experimental.")] = None,
    spectrum_style: Annotated[Optional[str], typer.Option(help="Override retained display style.")] = None,
    max_absorbance: Annotated[Optional[float], typer.Option(min=0.01, help="Override retained display absorbance.")] = None,
    phase: Annotated[Optional[str], typer.Option(help="Override retained sample phase.")] = None,
    measurement: Annotated[Optional[str], typer.Option(help="Override retained measurement geometry.")] = None,
    instrument_resolution: Annotated[Optional[float], typer.Option(min=0.1, help="Override retained FTIR resolution.")] = None,
    apodization: Annotated[Optional[str], typer.Option(help="Override retained apodization.")] = None,
    residual_fwhm: Annotated[Optional[float], typer.Option(min=0.1, help="Override retained residual linewidth.")] = None,
) -> None:
    """Build artifacts from completed conformers using retained settings by default."""
    from .spectrum import finalize_ir_spectrum

    try:
        result = finalize_ir_spectrum(
            run_dir, temperature=temperature, scale_factor=scale_factor,
            fwhm=fwhm, spectrum_style=spectrum_style, max_absorbance=max_absorbance,
            spectrum_model=spectrum_model,
            phase=phase, measurement=measurement, instrument_resolution=instrument_resolution,
            apodization=apodization, residual_fwhm=residual_fwhm,
        )
    except Exception as error:
        raise typer.BadParameter(f"Could not finalize spectrum: {error}") from error
    typer.echo(f"Spectrum CSV: {result['spectrum_csv']}")
    typer.echo(f"Spectrum PNG: {result['spectrum_png']}")


@app.command()
def orbitals(
    output_file: Annotated[Path, typer.Argument(help="ORCA output file to parse.")],
    json_output: Annotated[Optional[Path], typer.Option(help="Optional path for machine-readable results.")] = None,
) -> None:
    """Report HOMO and LUMO energies from an ORCA output file."""
    if not output_file.is_file():
        raise typer.BadParameter(f"ORCA output file not found: {output_file}")
    from src.parser import parse_orca_orbitals

    result = parse_orca_orbitals(output_file)
    if result["homo_number"] is None:
        raise typer.BadParameter("No occupied orbitals were found in this ORCA output.")
    typer.echo(f"HOMO {result['homo_number']}: {result['homo_energy']:.6f} eV")
    if result["lumo_number"] is None:
        typer.echo("LUMO: not found")
    else:
        typer.echo(f"LUMO {result['lumo_number']}: {result['lumo_energy']:.6f} eV")
        typer.echo(f"Gap: {result['lumo_energy'] - result['homo_energy']:.6f} eV")
    if json_output:
        json_output = Path(json_output)
        json_output.write_text(json.dumps(result, indent=2) + "\n")
        typer.echo(f"JSON: {json_output}")


@app.command()
def reactivity(
    output_file: Annotated[Path, typer.Argument(help="ORCA output containing an orbital-energy table.")],
    json_output: Annotated[Path, typer.Option(help="Machine-readable qualitative frontier summary.")] = Path("reactivity.json"),
) -> None:
    """Report qualitative molecule-level frontier-orbital reactivity proxies."""
    from .reactivity import frontier_reactivity_summary, write_frontier_reactivity

    if not output_file.is_file():
        raise typer.BadParameter(f"ORCA output file not found: {output_file}")
    try:
        result = frontier_reactivity_summary(output_file)
    except ValueError as error:
        raise typer.BadParameter(str(error)) from error
    write_frontier_reactivity(json_output, result)
    orbitals = result["frontier_orbitals"]
    proxies = result["derived_frontier_proxies"]
    typer.echo(f"HOMO/LUMO gap: {orbitals['gap_eV']:.4f} eV")
    typer.echo(f"Hardness proxy: {proxies['hardness_proxy_eV']:.4f} eV")
    typer.echo(f"Qualitative report: {json_output}")


@app.command()
def stability(
    smiles: Annotated[str, typer.Option(help="SMILES input to screen.")],
    scenario: Annotated[str, typer.Option(help="Declared homogeneous-gas scenario: low-pressure-intrinsic, dry inert, dry air, or humid air.")],
    output_dir: Annotated[Path, typer.Option(help="Parent directory for stability runs.")] = Path("runs"),
    name: Annotated[Optional[str], typer.Option(help="Name used for the run directory.")] = None,
    charge: Annotated[int, typer.Option(help="Molecular charge for the ORCA calculation.")] = 0,
    multiplicity: Annotated[int, typer.Option(min=1, help="Spin multiplicity for the ORCA calculation.")] = 1,
    cores: Annotated[int, typer.Option(min=1, help="ORCA CPU cores.")] = 1,
    method_profile: Annotated[str, typer.Option(help="Named ORCA method profile used for optimization and frequencies.")] = "b3lyp-def2-svp",
    rmg_env: Annotated[Optional[str], typer.Option(help="Optional Conda RMG environment; default uses PATH.")] = None,
    rmg_command: Annotated[str, typer.Option(help="RMG executable or command.")] = "rmg.py",
    temperature: Annotated[float, typer.Option(min=1.0, help="Reactor temperature in K.")] = 298.0,
    pressure: Annotated[float, typer.Option(min=1e-9, help="Reactor pressure in bar.")] = 1.0,
    barrier_threshold: Annotated[float, typer.Option(min=0.1, help="RMG route-energy discovery ceiling in kcal/mol; not a persistence threshold.")] = 50.0,
    screen_tier: Annotated[str, typer.Option(help="RMG resource tier: quick-screen or review-screen.")] = "quick-screen",
    target_duration_hours: Annotated[float, typer.Option(min=0.001, help="Declared condition duration in hours; used as the RMG simulation horizon.")] = 24.0,
    retention_fraction: Annotated[float, typer.Option(min=0.01, max=0.9999, help="Retention criterion recorded in the condition contract.")] = 0.95,
    light_condition: Annotated[str, typer.Option(help="Declared light condition; only 'dark' is currently supported.")] = "dark",
    rmg_walltime: Annotated[Optional[str], typer.Option(help="Optional RMG wall-time override as DD:HH:MM:SS.")] = None,
    rmg_max_processes: Annotated[Optional[int], typer.Option(min=1, help="Optional maximum RMG-process override.")] = None,
    rmg_max_iterations: Annotated[Optional[int], typer.Option(min=1, help="Optional RMG iteration-limit override.")] = None,
    rmg_max_edge_species: Annotated[Optional[int], typer.Option(min=1, help="Optional RMG edge-species-limit override.")] = None,
    rmg_maximum_heavy_atoms: Annotated[Optional[int], typer.Option(min=1, help="Override the generated-species heavy-atom cap; default admits one target-target collision product.")] = None,
    rmg_database_library: Annotated[list[str], typer.Option(help="Named kinetics library from the selected RMG database; repeatable.")] = [],
    auto_reference_libraries: Annotated[bool, typer.Option(help="Automatically enable curated RMG libraries matching the target's elemental scope.")] = True,
    reaction_library: Annotated[list[Path], typer.Option(help="Validated local generated RMG kinetics-library directory; repeatable.")] = [],
    auto_verify_routes: Annotated[bool, typer.Option(help="Automatically verify modeled threshold-crossing routes with generic ORCA paths, Arkane rates, and repaired RMG propagation.")] = True,
    verification_max_iterations: Annotated[int, typer.Option(min=1, help="Maximum flux-guided route verification/reinjection iterations.")] = 12,
    verification_timeout_seconds: Annotated[float, typer.Option(min=1, help="Timeout for each individual ORCA route job; default is four hours.")] = 14400.0,
    verification_orientations: Annotated[int, typer.Option(min=2, help="Maximum endpoint orientations; two agreeing paths are required for barrierless classification.")] = 3,
    verification_neb_images: Annotated[int, typer.Option(min=5, help="Movable ORCA NEB-TS images per route orientation.")] = 8,
) -> None:
    """Collect ORCA local-minimum and RMG pathway evidence under stated conditions."""
    from .stability import run_stability_screen
    from .methods import harmonic_method_profile

    run_dir = create_run_directory(output_dir, name or "stability")
    try:
        profile = harmonic_method_profile(method_profile)
        result = run_stability_screen(
            smiles,
            run_dir,
            charge=charge,
            multiplicity=multiplicity,
            ncores=cores,
            method_profile=profile,
            scenario=scenario,
            screen_tier=screen_tier,
            target_duration_hours=target_duration_hours,
            retention_fraction=retention_fraction,
            light_condition=light_condition,
            rmg_env=rmg_env,
            rmg_command=rmg_command,
            temperature=temperature,
            pressure=pressure,
            barrier_threshold=barrier_threshold,
            rmg_walltime=rmg_walltime,
            rmg_max_processes=rmg_max_processes,
            rmg_max_iterations=rmg_max_iterations,
            rmg_max_edge_species=rmg_max_edge_species,
            rmg_maximum_heavy_atoms=rmg_maximum_heavy_atoms,
            database_reaction_libraries=rmg_database_library,
            auto_reference_libraries=auto_reference_libraries,
            reaction_libraries=reaction_library,
            auto_verify_routes=auto_verify_routes,
            verification_max_iterations=verification_max_iterations,
            verification_timeout_seconds=verification_timeout_seconds,
            verification_orientations=verification_orientations,
            verification_neb_images=verification_neb_images,
        )
    except Exception:
        typer.echo(f"Run directory retained for diagnosis: {run_dir}", err=True)
        raise
    typer.echo(f"Stability screen completed: {run_dir}")
    typer.echo(f"Assessment: {result['assessment']['status']}")
    typer.echo(f"ORCA evidence: {result['orca_evidence']['status']}")
    typer.echo(f"RMG evidence: {result['rmg_evidence']['status']}")
    if result["risk_flags"]:
        typer.echo(f"Risk flags: {', '.join(result['risk_flags'])}")
    typer.echo(f"Result JSON: {result['result_json']}")


@app.command(name="stability-verify")
def stability_verify(
    run_dir: Annotated[Path, typer.Argument(help="Completed stability-screen directory.")],
    route: Annotated[int, typer.Option(min=0, help="Zero-based candidate route index to prepare for ORCA verification.")],
) -> None:
    """Prepare ORCA route verification after kinetic relevance has been established."""
    from .route_verify import prepare_route_verification

    try:
        result = prepare_route_verification(run_dir, route)
    except (FileNotFoundError, ValueError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Route-verification dossier: {result['path']}")


@app.command(name="stability-repair-collision")
def stability_repair_collision(
    run_dir: Annotated[Path, typer.Argument(help="Completed stability-screen directory containing stability.json.")],
    execute_endpoints: Annotated[bool, typer.Option(help="Run ORCA optimization/frequency jobs for all resolved endpoint species.")] = False,
    cores: Annotated[int, typer.Option(min=1, help="ORCA cores for each endpoint job.")] = 1,
) -> None:
    """Prepare ORCA repair work for RMG collision-limit kinetics violations."""
    from .collision_repair import prepare_collision_rate_repairs
    try:
        result = prepare_collision_rate_repairs(run_dir, execute_endpoints=execute_endpoints, ncores=cores)
    except (FileNotFoundError, ValueError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Collision-rate repair dossier: {result['path']}")


@app.command(name="stability-scan-association")
def stability_scan_association(
    run_dir: Annotated[Path, typer.Argument(help="Stability-screen directory after stability-repair-collision.")],
    execute: Annotated[bool, typer.Option(help="Execute the prepared ORCA doublet approach scans.")] = False,
    cores: Annotated[int, typer.Option(min=1, help="ORCA cores per scan.")] = 1,
    steps: Annotated[int, typer.Option(min=5, help="Relaxed-scan points per orientation.")] = 16,
) -> None:
    """Prepare or run spin-consistent scans for barrierless association candidates."""
    from .collision_repair import prepare_association_scans
    try:
        result = prepare_association_scans(run_dir, execute=execute, ncores=cores, steps=steps)
    except (FileNotFoundError, ValueError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Association-scan dossier: {result['path']}")


@app.command(name="stability-finalize-association")
def stability_finalize_association(
    run_dir: Annotated[Path, typer.Argument(help="Stability-screen directory after completed association scans.")],
) -> None:
    """Record whether completed ORCA scans support a barrierless or TS route."""
    from .collision_repair import finalize_association_verification
    try:
        result = finalize_association_verification(run_dir)
    except (FileNotFoundError, ValueError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Association verification: {result['status']}")
    typer.echo(f"Association verification dossier: {result['path']}")


@app.command(name="stability-repair-run")
def stability_repair_run(
    run_dir: Annotated[Path, typer.Argument(help="Collision-limited stability-screen directory.")],
    cores: Annotated[int, typer.Option(min=1, help="ORCA cores per job.")] = 1,
    scan_steps: Annotated[int, typer.Option(min=5, help="Relaxed-scan points per orientation.")] = 16,
) -> None:
    """Run endpoint and association evidence, then issue a fail-closed route result."""
    from .collision_repair import run_collision_repair_pipeline
    try:
        result = run_collision_repair_pipeline(run_dir, ncores=cores, scan_steps=scan_steps)
    except (FileNotFoundError, ValueError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Collision-repair assessment: {result['status']}")
    typer.echo(f"Collision-repair pipeline: {result['path']}")


@app.command(name="stability-explain")
def stability_explain(
    report: Annotated[Path, typer.Argument(help="Stability run directory or retained result JSON.")],
    route: Annotated[Optional[str], typer.Option(help="Route ID or zero-based route index; default selects the leading executable route.")] = None,
    output_dir: Annotated[Optional[Path], typer.Option(help="Explanation directory; defaults to <run>/explanation.")] = None,
    prepare_only: Annotated[bool, typer.Option(help="Build endpoints and a clearly labelled preview without running ORCA paths.")] = False,
    charge: Annotated[Optional[int], typer.Option(help="Override the route's total charge.")] = None,
    multiplicity: Annotated[Optional[int], typer.Option(min=1, help="Override the route's total spin multiplicity.")] = None,
    cores: Annotated[int, typer.Option(min=1, help="ORCA CPU cores.")] = 1,
    method_profile: Annotated[str, typer.Option(help="Named ORCA method profile used consistently across the path.")] = "b3lyp-def2-svp",
    timeout_seconds: Annotated[float, typer.Option(min=1, help="Timeout for each individual ORCA job; default is four hours.")] = 14400.0,
    scan_steps: Annotated[int, typer.Option(min=3, help="Relaxed dissociation-scan points.")] = 20,
    neb_images: Annotated[int, typer.Option(min=3, help="Movable images in an activated-route NEB-TS calculation.")] = 8,
) -> None:
    """Explain a selected instability with a computed path and labelled animation."""
    from .methods import harmonic_method_profile
    from .reaction_path import run_decomposition_explanation

    try:
        profile = harmonic_method_profile(method_profile)
        result = run_decomposition_explanation(
            report,
            output_dir=output_dir,
            route_id=route,
            prepare_only=prepare_only,
            charge=charge,
            multiplicity=multiplicity,
            ncores=cores,
            method_keywords=profile["orca_keywords"],
            timeout_seconds=timeout_seconds,
            scan_steps=scan_steps,
            neb_images=neb_images,
        )
    except (FileNotFoundError, ValueError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Decomposition explanation: {result['result_json']}")
    typer.echo(f"Status: {result['status']}")
    if result.get("visuals", {}).get("animation"):
        typer.echo(f"Animation: {result['visuals']['animation']}")
    if result.get("failure_reason"):
        typer.echo(f"Reason: {result['failure_reason']}", err=True)


@app.command(name="stability-ladder")
def stability_ladder(
    smiles: Annotated[str, typer.Option(help="SMILES input to screen through the ordered condition ladder.")],
    output_dir: Annotated[Path, typer.Option(help="Parent directory for ladder runs.")] = Path("runs"),
    name: Annotated[Optional[str], typer.Option(help="Name used for the run directory.")] = None,
    temperature: Annotated[float, typer.Option(min=273.15, max=323.15, help="Stage temperature in K.")] = 298.15,
    pressure: Annotated[float, typer.Option(min=0.01, help="Ambient-stage pressure in bar.")] = 1.0,
    relative_humidity: Annotated[float, typer.Option(min=0.0, max=1.0, help="Humid-air relative humidity as a fraction.")] = 0.5,
    retention_fraction: Annotated[float, typer.Option(min=0.01, max=0.9999)] = 0.95,
    maximum_duration_days: Annotated[float, typer.Option(min=0.01, help="Simulation horizon for each executable dark-gas stage.")] = 365.0,
    photolysis_evidence: Annotated[Optional[Path], typer.Option(help="JSON with cross section, quantum yield, and calculated photolysis rate.")] = None,
    computational_light_spectrum: Annotated[Optional[Path], typer.Option(help="Bounded TD-DFT light source CSV; defaults to bundled ASTM G173-03 AM1.5.")] = Path("references/astm-g173-am15g.csv"),
    rmg_env: Annotated[Optional[str], typer.Option(help="Conda RMG environment, e.g. rmg_env.")] = None,
    screen_tier: Annotated[str, typer.Option(help="RMG resource tier.")] = "quick-screen",
    charge: Annotated[int, typer.Option()] = 0,
    multiplicity: Annotated[int, typer.Option(min=1)] = 1,
    cores: Annotated[int, typer.Option(min=1)] = 1,
    light_nroots: Annotated[int, typer.Option(min=1, help="Full LR-TDDFT roots for each sunlight screen.")] = 20,
    rmg_workers: Annotated[Optional[int], typer.Option(min=1, help="RMG model-generation worker processes; defaults to the conservative tier setting (1).")]= None,
    rmg_objects_per_iteration: Annotated[int, typer.Option(min=1, help="RMG objects admitted per enlargement cycle; keep 1 unless an A/B mechanism comparison justifies more.")] = 1,
    rmg_maximum_heavy_atoms: Annotated[Optional[int], typer.Option(min=1, help="Override the generated-species heavy-atom cap for every ladder stage.")] = None,
    rmg_database_library: Annotated[list[str], typer.Option(help="Named kinetics library from the selected RMG database; repeatable.")] = [],
    auto_reference_libraries: Annotated[bool, typer.Option(help="Automatically enable curated RMG libraries matching the target's elemental scope.")] = True,
    rmg_retry_attempts: Annotated[int, typer.Option(min=1, max=3, help="Bounded seed-restart attempts after genuine RMG resource exhaustion.")] = 2,
    rmg_heartbeat_timeout_minutes: Annotated[Optional[float], typer.Option(min=1, help="Optional no-log-progress timeout for an RMG attempt.")] = None,
    rmg_start_method: Annotated[str, typer.Option(help="RMG multiprocessing start method: auto uses fork on macOS.")] = "auto",
    method_profile: Annotated[str, typer.Option()] = "b3lyp-def2-svp",
    reaction_library: Annotated[list[Path], typer.Option(help="Validated local generated RMG kinetics-library directory; repeatable.")] = [],
    auto_verify_routes: Annotated[bool, typer.Option(help="Automatically run flux-guided generic ORCA/Arkane verification and repaired RMG propagation.")] = True,
    verification_max_iterations: Annotated[int, typer.Option(min=1, help="Maximum route verification/reinjection iterations per ladder stage.")] = 12,
    verification_timeout_seconds: Annotated[float, typer.Option(min=1, help="Timeout for each individual ORCA route job; default is four hours.")] = 14400.0,
    verification_orientations: Annotated[int, typer.Option(min=2, help="Maximum endpoint orientations per route.")] = 3,
    verification_neb_images: Annotated[int, typer.Option(min=5, help="Movable ORCA NEB-TS images per orientation.")] = 8,
) -> None:
    """Screen an ordered low-pressure/inert/air/humidity/sunlight gas ladder."""
    from .condition_ladder import run_condition_ladder
    from .methods import harmonic_method_profile
    from .stability import run_stability_screen

    run_dir = create_run_directory(output_dir, name or "stability-ladder")
    profile = harmonic_method_profile(method_profile)
    shared_parent_orca: dict | None = None

    def dark_gas_stage(stage_smiles: str, stage_dir: Path, **stage_kwargs) -> dict:
        nonlocal shared_parent_orca
        stage_libraries = [*reaction_library, *stage_kwargs.pop("reaction_libraries", [])]
        stage_temperature = stage_kwargs.pop("temperature", temperature)
        stage_pressure = stage_kwargs.pop("pressure", pressure)
        stage_result = run_stability_screen(
            stage_smiles, stage_dir, scenario=stage_kwargs.pop("scenario"),
            scenario_config=stage_kwargs.pop("scenario_config"),
            target_duration_hours=stage_kwargs.pop("target_duration_hours"),
            retention_fraction=stage_kwargs.pop("retention_fraction"),
            light_condition=stage_kwargs.pop("light_condition", "dark"),
            light_model=stage_kwargs.pop("light_model", None),
            rmg_env=rmg_env, screen_tier=screen_tier, temperature=stage_temperature,
            pressure=stage_pressure, charge=charge, multiplicity=multiplicity, ncores=cores,
            # RMG 4 multiprocessing is not portable on macOS: spawned workers
            # can lose the initialized kinetics database.  ORCA cores and RMG
            # workers are therefore deliberately independent.
            rmg_max_processes=rmg_workers,
            rmg_max_objects_per_iteration=rmg_objects_per_iteration,
            rmg_maximum_heavy_atoms=rmg_maximum_heavy_atoms,
            rmg_retry_attempts=rmg_retry_attempts,
            rmg_heartbeat_timeout_seconds=(rmg_heartbeat_timeout_minutes * 60 if rmg_heartbeat_timeout_minutes else None),
            rmg_start_method=rmg_start_method,
            method_profile=profile,
            precomputed_orca_evidence=shared_parent_orca,
            reaction_libraries=stage_libraries,
            database_reaction_libraries=rmg_database_library,
            auto_reference_libraries=auto_reference_libraries,
            auto_verify_routes=auto_verify_routes,
            verification_max_iterations=verification_max_iterations,
            verification_timeout_seconds=verification_timeout_seconds,
            verification_orientations=verification_orientations,
            verification_neb_images=verification_neb_images,
        )
        if shared_parent_orca is None and stage_result.get("orca_evidence", {}).get("status") == "completed":
            shared_parent_orca = stage_result["orca_evidence"]
        return stage_result

    try:
        result = run_condition_ladder(
            smiles, run_dir, stage_runner=dark_gas_stage, temperature_K=temperature,
            pressure_bar=pressure, relative_humidity=relative_humidity,
            retention_fraction=retention_fraction, maximum_duration_days=maximum_duration_days,
            photolysis_evidence=photolysis_evidence,
            computational_light_spectrum=computational_light_spectrum,
            rmg_env=rmg_env,
            ncores=cores,
            charge=charge,
            multiplicity=multiplicity,
            method_keywords=profile["orca_keywords"],
            light_nroots=light_nroots,
        )
    except Exception:
        typer.echo(f"Run directory retained for diagnosis: {run_dir}", err=True)
        raise
    typer.echo(f"Condition ladder completed: {run_dir}")
    typer.echo(f"Verdict: {result['verdict']}")
    if result["first_t95_stage"]:
        typer.echo(f"First verified t95 stage: {result['first_t95_stage']}")
    typer.echo(f"Result JSON: {result['result_json']}")


@app.command(name="rmg-compare")
def rmg_compare(
    smiles: Annotated[str, typer.Option(help="SMILES input run through both RMG installations.")],
    output_dir: Annotated[Path, typer.Option(help="Parent directory for comparison runs.")] = Path("runs"),
    name: Annotated[Optional[str], typer.Option(help="Name used for the comparison directory.")] = None,
    rmg3_env: Annotated[str, typer.Option(help="Conda environment containing the established RMG baseline.")] = "rmg_env",
    rmg4_env: Annotated[str, typer.Option(help="Conda environment containing RMG 4.")] = "rmg4_env",
    scenario: Annotated[str, typer.Option(help="Identical declared scenario for both runs.")] = "ambient-air-gas-screen",
    screen_tier: Annotated[str, typer.Option(help="Identical RMG resource tier for both runs.")] = "quick-screen",
    temperature: Annotated[float, typer.Option(min=1.0)] = 298.15,
    pressure: Annotated[float, typer.Option(min=0.001)] = 1.0,
    target_duration_hours: Annotated[float, typer.Option(min=0.001)] = 24.0,
    retention_fraction: Annotated[float, typer.Option(min=0.01, max=0.9999)] = 0.95,
) -> None:
    """Compare RMG 3 and RMG 4 under one retained, matched contract."""
    from .rmg_compare import compare_rmg_versions

    run_dir = create_run_directory(output_dir, name or "rmg-compare")
    result = compare_rmg_versions(
        smiles, run_dir, rmg3_env=rmg3_env, rmg4_env=rmg4_env, scenario=scenario,
        screen_tier=screen_tier, temperature=temperature, pressure=pressure,
        target_duration_hours=target_duration_hours, retention_fraction=retention_fraction,
    )
    typer.echo(f"RMG comparison completed: {run_dir}")
    for item in result["runs"]:
        typer.echo(f"{item['label']}: RMG {item['version']['version'] or 'unknown'}, {item['status']}, {item['wall_seconds']:.1f}s")
    typer.echo(f"Result JSON: {result['result_json']}")


@app.command(name="photostability-screen")
def photostability_screen(
    smiles: Annotated[str, typer.Option(help="SMILES input for a TD-DFT sunlight-absorption screen.")],
    output_dir: Annotated[Path, typer.Option(help="Parent directory for retained photophysics runs.")] = Path("runs"),
    name: Annotated[Optional[str], typer.Option(help="Name used for the run directory.")] = None,
    charge: Annotated[int, typer.Option()] = 0,
    multiplicity: Annotated[int, typer.Option(min=1)] = 1,
    cores: Annotated[int, typer.Option(min=1)] = 1,
    nroots: Annotated[int, typer.Option(min=1, help="Number of TD-DFT roots to calculate.")] = 20,
    sunlight_spectrum: Annotated[Optional[Path], typer.Option(help="Declared spectrum CSV: wavelength_nm,irradiance_W_m2_nm.")] = None,
) -> None:
    """Screen whether TD-DFT predicts bright transitions in the solar window."""
    from .photophysics import run_photostability_screen

    if sunlight_spectrum and not sunlight_spectrum.is_file():
        raise typer.BadParameter(f"Sunlight spectrum not found: {sunlight_spectrum}")
    run_dir = create_run_directory(output_dir, name or "photostability")
    result = run_photostability_screen(smiles, run_dir, charge=charge, multiplicity=multiplicity,
                                       ncores=cores, nroots=nroots, sunlight_csv=sunlight_spectrum)
    typer.echo(f"Photostability screen: {run_dir}")
    typer.echo(f"Assessment: {result['assessment']}")
    typer.echo(f"TD-DFT status: {result['status']}")
    typer.echo(f"Result JSON: {result['result_json']}")


@app.command(name="computational-light")
def computational_light(
    smiles: Annotated[str, typer.Option(help="SMILES input for the bounded computational sunlight model.")],
    sunlight_spectrum: Annotated[Path, typer.Option(help="Source spectrum CSV: wavelength_nm,irradiance_W_m2_nm.")] = Path("references/astm-g173-am15g.csv"),
    output_dir: Annotated[Path, typer.Option()] = Path("runs"),
    name: Annotated[Optional[str], typer.Option()] = None,
    rmg_env: Annotated[Optional[str], typer.Option(help="RMG environment used to write compatible local libraries.")] = "rmg_env",
    charge: Annotated[int, typer.Option()] = 0,
    multiplicity: Annotated[int, typer.Option(min=1)] = 1,
    cores: Annotated[int, typer.Option(min=1)] = 1,
    nroots: Annotated[int, typer.Option(min=1)] = 16,
    simulate: Annotated[bool, typer.Option(help="Propagate each generated photon profile through bounded RMG.")] = True,
    scenario: Annotated[str, typer.Option(help="RMG condition used when --simulate is enabled.")] = "ambient-air-gas-screen",
    target_duration_hours: Annotated[float, typer.Option(min=0.001)] = 24.0,
) -> None:
    """Generate low/central/high photon reactions from bounded ORCA calculations."""
    from .computational_light import run_computational_light_model, simulate_computational_light_profiles

    if not sunlight_spectrum.is_file():
        raise typer.BadParameter(f"Sunlight spectrum not found: {sunlight_spectrum}")
    run_dir = create_run_directory(output_dir, name or "computational-light")
    result = run_computational_light_model(smiles, run_dir, sunlight_spectrum=sunlight_spectrum,
                                           rmg_env=rmg_env, charge=charge, multiplicity=multiplicity,
                                           ncores=cores, nroots=nroots)
    typer.echo(f"Computational light model: {run_dir}")
    typer.echo(f"Status: {result['status']}")
    for profile, values in result.get("profiles", {}).items():
        typer.echo(f"{profile}: {len(values['channels'])} generated photon channel(s)")
    if simulate and result["status"] == "completed":
        kinetics = simulate_computational_light_profiles(result, run_dir / "rmg-profiles", rmg_env=rmg_env,
                                                          scenario=scenario, target_duration_hours=target_duration_hours)
        for profile, values in kinetics["profiles"].items():
            typer.echo(f"{profile} t95: {values['estimated_time_to_retention_seconds']}")
        typer.echo(f"Kinetic profiles: {kinetics['result_json']}")
    typer.echo(f"Result JSON: {result['result_json']}")


@app.command(name="arkane-route")
def arkane_route(
    label: Annotated[str, typer.Option(help="Verified elementary-route label.")],
    reactant_label: Annotated[str, typer.Option()],
    reactant_smiles: Annotated[str, typer.Option()],
    reactant_orca: Annotated[Path, typer.Option(help="Completed reactant ORCA output.")],
    product: Annotated[list[str], typer.Option(help="Repeat LABEL:SMILES:ORCA_OUTPUT for each product.")],
    transition_state_orca: Annotated[Path, typer.Option(help="Completed TS ORCA output with one imaginary mode.")],
    output_dir: Annotated[Path, typer.Option()] = Path("runs"),
    name: Annotated[Optional[str], typer.Option()] = None,
    temperatures: Annotated[str, typer.Option(help="Comma-separated temperature grid in K.")] = "298.15",
    pressures: Annotated[str, typer.Option(help="Comma-separated pressure grid in bar.")] = "1e-6,1.0",
    bath_gas: Annotated[str, typer.Option(help="Comma-separated NAME:FRACTION entries.")] = "nitrogen:1.0",
    rmg_env: Annotated[Optional[str], typer.Option()] = "rmg_env",
) -> None:
    """Calculate pressure-dependent kinetics from verified ORCA route artifacts."""
    from .arkane_runner import ArkaneRouteSpec, run_arkane_route

    def parse_product(value: str) -> tuple[str, str, Path]:
        fields = value.split(":", 2)
        if len(fields) != 3:
            raise typer.BadParameter("Each --product must be LABEL:SMILES:ORCA_OUTPUT")
        return fields[0], fields[1], Path(fields[2])
    def parse_grid(value: str) -> tuple[float, ...]:
        try:
            values = tuple(float(item) for item in value.split(",") if item)
        except ValueError as error:
            raise typer.BadParameter("Grid values must be comma-separated numbers") from error
        if not values or any(item <= 0 for item in values):
            raise typer.BadParameter("Grid values must be positive")
        return values
    try:
        products = [parse_product(value) for value in product]
        gases = tuple((item.split(":", 1)[0], float(item.split(":", 1)[1])) for item in bath_gas.split(","))
        spec = ArkaneRouteSpec(
            label=label, reactant_label=reactant_label, reactant_smiles=reactant_smiles,
            reactant_orca_output=reactant_orca, product_labels=tuple(item[0] for item in products),
            product_smiles=tuple(item[1] for item in products), product_orca_outputs=tuple(item[2] for item in products),
            transition_state_orca_output=transition_state_orca, temperatures_K=parse_grid(temperatures),
            pressures_bar=parse_grid(pressures), bath_gas=gases,
        )
    except (ValueError, IndexError) as error:
        raise typer.BadParameter(str(error)) from error
    run_dir = create_run_directory(output_dir, name or "arkane-route")
    result = run_arkane_route(spec, run_dir, rmg_env=rmg_env)
    typer.echo(f"Arkane route status: {result['status']}")
    typer.echo(f"Artifacts: {run_dir}")


if __name__ == "__main__":
    app()
