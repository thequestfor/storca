"""Command-line interface for STORCA."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
from pathlib import Path
from typing import Annotated, Optional

import typer

from .runs import create_run_directory, write_metadata
from .workflow import run_optimization_and_frequency

app = typer.Typer(help="Reproducible ORCA calculation workflows.", no_args_is_help=True)


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
    checks = {
        "ORCA": os.environ.get("STORCA_ORCA_BIN") or shutil.which("orca"),
        "xTB": shutil.which("xtb"),
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
    conformer_engine: Annotated[str, typer.Option(help="Conformer engine: goat (recommended) or rdkit.")] = "goat",
    initial_conformers: Annotated[int, typer.Option(min=1, help="RDKit conformers to generate in fallback mode.")] = 20,
    max_conformers: Annotated[int, typer.Option(min=1, help="Maximum conformers sent to ORCA frequencies.")] = 10,
    goat_population: Annotated[float, typer.Option(min=0.01, max=1.0, help="Cumulative GOAT population to retain.")] = 0.95,
    scale_factor: Annotated[Optional[float], typer.Option(min=0.1, max=1.2, help="Harmonic scale-factor override; otherwise use --method-profile.")] = None,
    method_profile: Annotated[str, typer.Option(help="Named harmonic method/scale-factor profile.")] = "b3lyp-def2-svp",
    fwhm: Annotated[float, typer.Option(min=1.0, help="Gaussian band width (FWHM) in cm^-1.")] = 15.0,
    spectrum_model: Annotated[str, typer.Option(help="Spectrum model: raw ORCA bands or practical rule-calibrated display.")] = "raw",
    spectrum_style: Annotated[str, typer.Option(help="Display: transmittance, absorbance, or relative.")] = "transmittance",
    max_absorbance: Annotated[float, typer.Option(min=0.01, help="Deepest plotted absorbance for transmittance display.")] = 1.0,
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
        engine = conformer_engine.lower()
        if engine not in {"goat", "rdkit"}:
            raise typer.BadParameter("Conformer engine must be 'goat' or 'rdkit'.")
        spectrum_model = spectrum_model.lower()
        if spectrum_model not in {"raw", "practical"}:
            raise typer.BadParameter("spectrum_model must be 'raw' or 'practical'.")
        seed_xyz = run_dir / "seed.xyz"
        smiles_to_xyz(canonical_smiles, seed_xyz)
        if engine == "goat":
            xyz_paths = run_goat_search(
                seed_xyz, run_dir, charge=charge, multiplicity=multiplicity,
                ncores=cores, population_cutoff=goat_population,
                max_conformers=max_conformers, progress=typer.echo,
            )
        else:
            initial_dir = run_dir / "initial-conformers"
            xyz_paths = smiles_to_conformers(canonical_smiles, initial_dir, n_confs=initial_conformers, max_confs=max_conformers)
            write_metadata(run_dir, conformer_engine="rdkit", smiles=canonical_smiles, initial_conformers=len(xyz_paths))
        write_metadata(
            run_dir, workflow="ir_spectrum", smiles=canonical_smiles, charge=charge,
            multiplicity=multiplicity, ncores=cores, temperature=temperature,
            scale_factor=resolved_scale_factor, scale_factor_source=scale_factor_source,
            harmonic_method_profile=profile, vibrational_treatment="harmonic",
            spectrum_model=spectrum_model,
            spectrum_model_description="Boltzmann-weighted conformer ensemble with Gaussian broadening",
            fwhm_cm_1=fwhm, spectrum_style=spectrum_style,
            max_absorbance=max_absorbance, status="running",
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
            method_keywords=profile["orca_keywords"],
            progress=typer.echo,
        )
    except Exception:
        typer.echo(f"Run directory retained for diagnosis: {run_dir}", err=True)
        raise
    typer.echo(f"Spectrum run completed: {run_dir}")
    typer.echo(f"Spectrum CSV: {result['spectrum_csv']}")
    typer.echo(f"Spectrum PNG: {result['spectrum_png']}")
    if "raw_spectrum_csv" in result:
        typer.echo(f"Raw calculated spectrum CSV: {result['raw_spectrum_csv']}")


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
) -> None:
    """Continue incomplete ORCA jobs using the original run settings by default."""
    from .spectrum import resume_ir_spectrum
    result = resume_ir_spectrum(
        run_dir, charge=charge, multiplicity=multiplicity, ncores=cores,
        temperature=temperature, scale_factor=scale_factor, fwhm=fwhm,
        spectrum_style=spectrum_style, max_absorbance=max_absorbance,
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
    spectrum_model: Annotated[Optional[str], typer.Option(help="Override retained model: raw or practical.")] = None,
    spectrum_style: Annotated[Optional[str], typer.Option(help="Override retained display style.")] = None,
    max_absorbance: Annotated[Optional[float], typer.Option(min=0.01, help="Override retained display absorbance.")] = None,
) -> None:
    """Build artifacts from completed conformers using retained settings by default."""
    from .spectrum import finalize_ir_spectrum

    try:
        result = finalize_ir_spectrum(
            run_dir, temperature=temperature, scale_factor=scale_factor,
            fwhm=fwhm, spectrum_style=spectrum_style, max_absorbance=max_absorbance,
            spectrum_model=spectrum_model,
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
            rmg_retry_attempts=rmg_retry_attempts,
            rmg_heartbeat_timeout_seconds=(rmg_heartbeat_timeout_minutes * 60 if rmg_heartbeat_timeout_minutes else None),
            rmg_start_method=rmg_start_method,
            method_profile=profile,
            precomputed_orca_evidence=shared_parent_orca,
            reaction_libraries=stage_libraries,
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
