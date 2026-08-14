"""Conformer-weighted IR spectrum calculations without Multiwfn."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from src.inputgen import create_orca_input
from src.orca_runner import find_xtb, run_orca
from src.parser import attach_xyz_to_conformers, parse_goat_out, parse_orca_energy, parse_orca_ir, parse_xyz_ensemble
from src.stability.freq_check import frequency_stability_check

from .runs import write_metadata

HARTREE_BOLTZMANN_PER_K = 3.166811563e-6


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _calculation_contract(source_xyz: str, *, charge: int, multiplicity: int,
                          method_keywords: list[str] | None,
                          geometry_role: str = "optimized_minimum") -> dict:
    return {
        "schema_version": 1,
        "source_xyz_sha256": _sha256_text(source_xyz),
        "charge": int(charge),
        "multiplicity": int(multiplicity),
        "method_keywords": list(method_keywords) if method_keywords is not None else None,
        "geometry_role": geometry_role,
    }


def _legacy_contract_matches(metadata: dict, previous_xyz: str | None, source_xyz: str,
                             contract: dict) -> bool:
    """Conservatively recognize pre-signature jobs from retained run metadata."""
    if previous_xyz != source_xyz:
        return False
    retained_keywords = (metadata.get("harmonic_method_profile") or {}).get("orca_keywords")
    if retained_keywords is None:
        retained_keywords = metadata.get("orca_method_keywords")
    return (
        contract.get("geometry_role") == "optimized_minimum"
        and
        int(metadata.get("charge", contract["charge"])) == contract["charge"]
        and int(metadata.get("multiplicity", contract["multiplicity"])) == contract["multiplicity"]
        and (list(retained_keywords) if retained_keywords is not None else None) == contract["method_keywords"]
    )


def select_goat_conformers(conformers: Iterable[dict], *, population_cutoff: float = 0.95, max_conformers: int = 10) -> list[dict]:
    """Select population-ranked GOAT conformers for higher-level calculations."""
    if not 0 < population_cutoff <= 1 or max_conformers < 1:
        raise ValueError("Invalid GOAT population cutoff or conformer limit")
    selected: list[dict] = []
    for conformer in conformers:
        if len(selected) >= max_conformers:
            break
        selected.append(conformer)
        if conformer["cum_weight"] >= population_cutoff:
            break
    if not selected:
        raise RuntimeError("GOAT did not report any conformers")
    return selected


def run_goat_search(seed_xyz: Path, run_dir: Path, *, charge: int = 0, multiplicity: int = 1, ncores: int = 1, population_cutoff: float = 0.95, max_conformers: int = 10, progress: Callable[[str], None] | None = None) -> list[Path]:
    """Run ORCA GOAT/xTB and write population-selected XYZ conformers."""
    run_dir = Path(run_dir)
    goat_dir = run_dir / "goat"
    goat_dir.mkdir(parents=True, exist_ok=True)
    input_xyz = goat_dir / "seed.xyz"
    input_xyz.write_text(Path(seed_xyz).read_text())
    if progress:
        progress("GOAT/xTB conformer search: running")
    goat_input = create_orca_input(input_xyz, charge=charge, multiplicity=multiplicity, use_goat=True, label="goat", ncores=ncores)
    xtb_executable = find_xtb()
    goat_outputs = run_orca(goat_input, extra_env={"XTBEXE": xtb_executable})
    final_ensemble = goat_dir / "goat.finalensemble.xyz"
    if not final_ensemble.is_file():
        raise RuntimeError(f"GOAT completed without a final ensemble: {final_ensemble}")
    conformers, _ = parse_goat_out(goat_outputs["out"])
    ensemble = attach_xyz_to_conformers(conformers, parse_xyz_ensemble(final_ensemble))
    selected = select_goat_conformers(ensemble, population_cutoff=population_cutoff, max_conformers=max_conformers)
    selected_dir = run_dir / "selected-conformers"
    selected_dir.mkdir(exist_ok=True)
    paths: list[Path] = []
    for position, conformer in enumerate(selected, start=1):
        path = selected_dir / f"conf_{position:03d}.xyz"
        coordinates = conformer["xyz"]["xyz_lines"]
        path.write_text(f"{len(coordinates)}\nGOAT conformer {conformer['idx']}; population {conformer['weight']:.4f}\n" + "\n".join(coordinates) + "\n")
        paths.append(path)
    selection_manifest = run_dir / "selected-conformers.json"
    selection_manifest.write_text(json.dumps({
        "schema_version": 1,
        "requested_population_cutoff": population_cutoff,
        "requested_max_conformers": max_conformers,
        "retained_cumulative_population": selected[-1]["cum_weight"],
        "cutoff_reached": selected[-1]["cum_weight"] >= population_cutoff,
        "truncated_by_max_conformers": len(selected) >= max_conformers and selected[-1]["cum_weight"] < population_cutoff,
        "conformers": [
            {
                "selected_position": position,
                "goat_index": conformer["idx"],
                "goat_relative_energy": conformer["energy"],
                "goat_degeneracy": conformer["degeneracy"],
                "goat_population": conformer["weight"],
                "goat_cumulative_population": conformer["cum_weight"],
                "xyz": str(paths[position - 1]),
            }
            for position, conformer in enumerate(selected, start=1)
        ],
    }, indent=2, sort_keys=True) + "\n")
    write_metadata(run_dir, conformer_engine="goat", goat_population_cutoff=population_cutoff,
                   goat_retained_cumulative_population=selected[-1]["cum_weight"],
                   goat_population_cutoff_reached=selected[-1]["cum_weight"] >= population_cutoff,
                   goat_conformers_reported=len(conformers), goat_conformers_selected=len(paths),
                   goat_selection_manifest=str(selection_manifest), goat_output=str(goat_outputs["out"]))
    if progress:
        progress(f"GOAT selected {len(paths)} conformers covering {selected[-1]['cum_weight']:.1%} of its ensemble")
    return paths


def boltzmann_weights(energies: Iterable[float], temperature: float,
                      degeneracies: Iterable[float] | None = None) -> list[float]:
    """Return normalized Boltzmann weights for energies in Hartree."""
    values = list(energies)
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("Temperature must be greater than zero")
    if not values or not all(math.isfinite(energy) for energy in values):
        raise ValueError("At least one finite conformer energy is required")
    degeneracy_values = list(degeneracies) if degeneracies is not None else [1.0] * len(values)
    if (len(degeneracy_values) != len(values)
            or not all(math.isfinite(value) and value > 0 for value in degeneracy_values)):
        raise ValueError("Conformer degeneracies must be finite, positive, and match the energies")
    minimum = min(values)
    factors = [degeneracy * math.exp(-(energy - minimum) / (HARTREE_BOLTZMANN_PER_K * temperature))
               for energy, degeneracy in zip(values, degeneracy_values)]
    total = sum(factors)
    return [factor / total for factor in factors]


def build_ir_spectrum(
    conformers: Iterable[dict],
    *,
    frequency_min: float = 400.0,
    frequency_max: float = 4000.0,
    resolution: float = 1.0,
    fwhm: float = 15.0,
    scale_factor: float = 0.97,
    band_transform: Callable[[float, float], tuple[float, float]] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """Build a relative, Gaussian-broadened IR spectrum.

    Only conformers carrying finite energies and at least one positive IR mode
    are included. Their Boltzmann weights are recalculated so failed frequency
    jobs cannot dilute the resulting spectrum.
    """
    if (not all(math.isfinite(value) for value in (frequency_min, frequency_max, resolution, fwhm, scale_factor))
            or frequency_min >= frequency_max or resolution <= 0 or fwhm <= 0 or scale_factor <= 0):
        raise ValueError("Invalid spectrum range, resolution, or FWHM")
    valid: list[dict] = []
    for conformer in conformers:
        energy = conformer.get("energy")
        modes = conformer.get("ir_modes", [])
        if not isinstance(energy, (int, float)) or not math.isfinite(energy):
            continue
        frequency_check = conformer.get("frequency_check")
        if (frequency_check is not None and frequency_check.get("IsMinimum") is not True
                and conformer.get("geometry_role") != "environment_snapshot"):
            continue
        usable_modes = []
        for mode in modes:
            try:
                frequency = float(mode["freq"])
                intensity = float(mode["intensity"])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(frequency) and frequency > 0 and math.isfinite(intensity) and intensity >= 0:
                usable_modes.append({**mode, "freq": frequency, "intensity": intensity})
        if usable_modes:
            valid.append({**conformer, "ir_modes": usable_modes})
    if not valid:
        raise RuntimeError("No conformers completed both energy and IR parsing successfully")

    explicit_weights = [conformer.get("population_weight") for conformer in valid]
    if any(weight is not None for weight in explicit_weights):
        if not all(
            isinstance(weight, (int, float)) and math.isfinite(weight) and weight > 0
            for weight in explicit_weights
        ):
            raise ValueError("Explicit population weights must be finite, positive, and present for every record")
        total_weight = sum(float(weight) for weight in explicit_weights)
        weights = [float(weight) / total_weight for weight in explicit_weights]
    else:
        weights = boltzmann_weights(
            [conformer["energy"] for conformer in valid],
            temperature=conformers_temperature(valid),
            degeneracies=[float(conformer.get("degeneracy", 1.0)) for conformer in valid],
        )
    for conformer, weight in zip(valid, weights):
        conformer["weight"] = weight

    grid = np.arange(frequency_min, frequency_max + resolution * 0.5, resolution)
    spectrum = np.zeros_like(grid, dtype=float)
    for conformer in valid:
        for mode in conformer["ir_modes"]:
            frequency = mode["freq"] * scale_factor
            band_fwhm = fwhm
            if band_transform:
                frequency, band_fwhm = band_transform(frequency, fwhm)
            if not math.isfinite(frequency) or not math.isfinite(band_fwhm) or band_fwhm <= 0:
                raise ValueError("Band transform returned an invalid frequency or FWHM")
            band_sigma = band_fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
            # ORCA IR intensities are integrated band strengths.  Normalize
            # each Gaussian so changing a display width cannot create or
            # destroy integrated intensity.
            gaussian = np.exp(-0.5 * ((grid - frequency) / band_sigma) ** 2) / (band_sigma * math.sqrt(2.0 * math.pi))
            spectrum += conformer["weight"] * mode["intensity"] * gaussian
    return grid, spectrum, valid


def conformers_temperature(conformers: list[dict]) -> float:
    """Read the common workflow temperature stored on conformer records."""
    temperatures = {conformer.get("temperature") for conformer in conformers}
    if len(temperatures) != 1 or None in temperatures:
        raise ValueError("Conformers must share one explicit temperature")
    temperature = float(temperatures.pop())
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("Conformer temperature must be finite and greater than zero")
    return temperature


def format_spectrum(intensity: np.ndarray, *, style: str, max_absorbance: float = 1.0) -> tuple[np.ndarray, str]:
    """Convert relative line intensities into a requested display convention."""
    intensity = np.asarray(intensity, dtype=float)
    if intensity.size == 0 or not np.all(np.isfinite(intensity)) or np.any(intensity < 0):
        raise ValueError("Spectrum intensity must contain finite, nonnegative values")
    if style == "relative":
        peak = float(intensity.max())
        return (np.zeros_like(intensity) if peak == 0 else intensity / peak), "relative_intensity"
    if style not in {"absorbance", "transmittance"}:
        raise ValueError("Spectrum style must be relative, absorbance, or transmittance")
    if not math.isfinite(max_absorbance) or max_absorbance <= 0:
        raise ValueError("Maximum absorbance must be greater than zero")
    peak = float(intensity.max())
    absorbance = np.zeros_like(intensity) if peak == 0 else intensity / peak * max_absorbance
    if style == "absorbance":
        return absorbance, "absorbance"
    return 100.0 * np.power(10.0, -absorbance), "transmittance_percent"


def write_spectrum_csv(path: Path, grid: np.ndarray, intensity: np.ndarray, column_name: str = "relative_intensity") -> Path:
    path = Path(path)
    grid = np.asarray(grid, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    if (grid.ndim != 1 or intensity.ndim != 1 or len(grid) == 0 or len(grid) != len(intensity)
            or not np.all(np.isfinite(grid)) or not np.all(np.isfinite(intensity))):
        raise ValueError("Spectrum grid and intensity must be equal-length finite one-dimensional arrays")
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber_cm-1", column_name])
        writer.writerows(zip(grid, intensity))
    return path


def write_spectrum_plot(path: Path, grid: np.ndarray, intensity: np.ndarray, y_label: str, *, title: str = "Boltzmann-weighted harmonic IR spectrum") -> Path:
    """Write a PNG of the relative IR spectrum with the conventional x-axis."""
    import matplotlib.pyplot as plt

    path = Path(path)
    figure, axis = plt.subplots(figsize=(9, 4.5))
    axis.plot(grid, intensity, color="black", linewidth=1.2)
    axis.invert_xaxis()
    axis.set_xlabel("Wavenumber (cm$^{-1}$)")
    axis.set_ylabel(y_label)
    axis.set_title(title)
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def _write_conformer_summary(path: Path, conformers: list[dict]) -> Path:
    serializable = []
    for conformer in conformers:
        serializable.append({
            "index": conformer["index"],
            "energy_hartree": conformer.get("energy"),
            "weight": conformer.get("weight"),
            "population_weight": conformer.get("population_weight"),
            "population_model": conformer.get("population_model"),
            "geometry_role": conformer.get("geometry_role"),
            "geometry_cluster_id": conformer.get("geometry_cluster_id"),
            "independent_environment_id": conformer.get("independent_environment_id"),
            "environment_features": conformer.get("environment_features"),
            "cluster_size": conformer.get("cluster_size"),
            "topology": conformer.get("topology"),
            "molecule_atom_ranges": conformer.get("molecule_atom_ranges"),
            "local_stretch_bonds": conformer.get("local_stretch_bonds"),
            "hydrogen_bond_interactions": conformer.get("hydrogen_bond_interactions"),
            "degeneracy": conformer.get("degeneracy", 1),
            "source_conformer": conformer.get("source_conformer"),
            "optimized_xyz": str(conformer.get("optimized_xyz", "")),
            "frequency_output": str(conformer.get("frequency_output", "")),
            "frequency_check": conformer.get("frequency_check"),
            "snapshot_hessian_reliability": conformer.get("snapshot_hessian_reliability"),
            "error": conformer.get("error"),
        })
    path.write_text(json.dumps(serializable, indent=2) + "\n")
    return path


def _completed_orca_output(path: Path) -> bool:
    return path.is_file() and "ORCA TERMINATED NORMALLY" in path.read_text(errors="replace")


def _completed_optimization(path: Path) -> bool:
    if not _completed_orca_output(path):
        return False
    text = path.read_text(errors="replace")
    return "THE OPTIMIZATION HAS CONVERGED" in text


def _backfill_local_environment_roles(record: dict) -> None:
    bonds = record.get("local_stretch_bonds") or []
    if not bonds or all(item.get("spectral_band_class") for item in bonds):
        return
    source = record.get("source_conformer") or {}
    source_xyz = source.get("source_xtb_xyz")
    ranges = record.get("molecule_atom_ranges") or []
    if not source_xyz or not ranges:
        return
    xcontrol = Path(source_xyz).parent / "xcontrol.inp"
    if not xcontrol.is_file():
        return
    from .cluster_ir import classify_local_stretch_bonds, parse_xcontrol_interactions
    interactions = parse_xcontrol_interactions(xcontrol)
    atom_count = int(ranges[0][1]) - int(ranges[0][0])
    record["hydrogen_bond_interactions"] = interactions
    record["local_stretch_bonds"] = classify_local_stretch_bonds(bonds, interactions, atom_count)


def _snapshot_hessian_reliability(frequency_check: dict | None) -> dict:
    check = frequency_check or {}
    imaginary = [float(value) for value in check.get("ImaginaryFrequencies", [])]
    material = [value for value in imaginary if value < -50.0]
    return {
        "stationarity_status": "poor" if material else "acceptable_snapshot",
        "full_hessian_use": "diagnostic_only" if material else "usable_with_snapshot_warning",
        "material_imaginary_modes": len(material),
        "lowest_frequency_cm-1": min(imaginary) if imaginary else None,
        "gradient_rms_hartree_per_bohr": None,
        "gradient_status": "not_available_from_retained_frequency_contract",
    }


def load_completed_ir_records(run_dir: Path, *, temperature: float = 298.15) -> list[dict]:
    """Load completed ORCA conformer results without writing or rerunning anything."""
    run_dir = Path(run_dir)
    conformer_dirs = sorted((run_dir / "conformers").glob("conf-*"))
    if not conformer_dirs:
        raise FileNotFoundError(f"No conformer directories found in {run_dir}")
    selection_path = run_dir / "selected-conformers.json"
    selection = json.loads(selection_path.read_text()) if selection_path.is_file() else {}
    selection_by_position = {item["selected_position"]: item for item in selection.get("conformers", [])}
    records: list[dict] = []
    for index, job_dir in enumerate(conformer_dirs, start=1):
        opt_out = job_dir / "opt.out"
        freq_out = job_dir / "freq.out"
        source = selection_by_position.get(index, {})
        geometry_role = source.get("geometry_role", "optimized_minimum")
        record = {"index": index, "temperature": temperature,
                  "degeneracy": source.get("goat_degeneracy", 1),
                  "source_conformer": source or None,
                  "geometry_role": geometry_role}
        for key in (
            "population_weight", "population_model", "population_warning",
            "geometry_cluster_id", "independent_environment_id", "environment_features",
            "sampling_support_count", "source_xtb_candidate_id",
            "cluster_size", "topology", "molecule_atom_ranges", "local_stretch_bonds",
            "hydrogen_bond_interactions",
            "environment_refinement", "environment_refinement_artifact",
            "vibrational_route", "gradient_rms_hartree_per_bohr",
            "gradient_maximum_component_hartree_per_bohr",
        ):
            if key in source:
                record[key] = source[key]
        _backfill_local_environment_roles(record)
        if geometry_role != "environment_snapshot" and not _completed_optimization(opt_out):
            record["error"] = "Optimization is missing, unconverged, or did not terminate normally"
        elif not _completed_orca_output(freq_out):
            record["error"] = "Frequency calculation is missing or did not terminate normally"
        else:
            try:
                record["optimized_xyz"] = (
                    job_dir / "input.xyz" if geometry_role == "environment_snapshot"
                    else job_dir / "opt.xyz"
                )
                record["energy"] = parse_orca_energy(
                    freq_out if geometry_role == "environment_snapshot" else opt_out
                )
                record["frequency_output"] = freq_out
                record["frequency_check"] = frequency_stability_check(freq_out)
                if geometry_role == "environment_snapshot":
                    record["snapshot_hessian_reliability"] = _snapshot_hessian_reliability(record["frequency_check"])
                record["ir_modes"] = parse_orca_ir(freq_out)
                if (not record["frequency_check"].get("IsMinimum")
                        and geometry_role != "environment_snapshot"):
                    record["error"] = "Frequency calculation does not describe a local minimum"
            except Exception as error:
                record["error"] = str(error)
        records.append(record)
    return records


def _write_ir_artifacts(
    records: list[dict],
    run_dir: Path,
    *,
    scale_factor: float,
    fwhm: float,
    spectrum_style: str,
    max_absorbance: float,
    spectrum_model: str,
    practical_smiles: str | None,
    phase: str = "liquid",
    measurement: str = "auto",
    instrument_resolution: float = 4.0,
    apodization: str = "happ-genzel",
    residual_fwhm: float | None = None,
    environment_sampling: bool = False,
    environment_convergence: dict | None = None,
    temperature_K: float | None = None,
    experimental_condition_details: dict | None = None,
) -> dict:
    """Build raw artifacts and an optional practical or measured display."""
    if spectrum_model not in {"raw", "practical", "experimental"}:
        raise ValueError("spectrum_model must be 'raw', 'practical', or 'experimental'")
    raw_grid, raw_intensity, successful = build_ir_spectrum(
        records, scale_factor=scale_factor, fwhm=fwhm,
    )
    weights_by_index = {item["index"]: item["weight"] for item in successful}
    for record in records:
        record["weight"] = weights_by_index.get(record["index"])
    raw_display, column_name = format_spectrum(raw_intensity, style=spectrum_style, max_absorbance=max_absorbance)
    if spectrum_model == "raw":
        csv_path = write_spectrum_csv(run_dir / "spectrum.csv", raw_grid, raw_display, column_name)
        png_path = write_spectrum_plot(run_dir / "spectrum.png", raw_grid, raw_display, column_name.replace("_", " "))
        return {"spectrum_csv": csv_path, "spectrum_png": png_path, "successful": successful,
                "practical_profile": None, "experimental_profile": None}

    if spectrum_model == "experimental":
        from .experimental_ir import (apply_environment_convergence_status, apply_measurement_response,
                                      collapse_insufficient_environment_bands,
                                      ensemble_band_statistics, environment_sufficiency_metadata,
                                      experimental_profile_metadata, resolve_experimental_profile)

        profile = resolve_experimental_profile(
            phase=phase, measurement=measurement,
            instrument_resolution_cm_1=instrument_resolution,
            apodization=apodization, residual_fwhm_cm_1=residual_fwhm,
        )
        from .ir_contracts import build_experimental_condition, write_experimental_condition
        if temperature_K is None:
            temperature_K = next((
                float(record["temperature"]) for record in records
                if isinstance(record.get("temperature"), (int, float))
            ), None)
        details = dict(experimental_condition_details or {})
        condition = build_experimental_condition(
            phase=profile.phase, measurement=profile.measurement,
            temperature_K=temperature_K,
            pressure_bar=details.get("pressure_bar"),
            composition=details.get("composition"),
            solvent=details.get("solvent"),
            concentration_mol_L=details.get("concentration_mol_L"),
            resolution_cm_1=profile.instrument_resolution_cm_1,
            apodization=profile.apodization,
            path_length_mm=details.get("path_length_mm"),
            atr_crystal=details.get("atr_crystal"),
            atr_incidence_angle_degrees=details.get("atr_incidence_angle_degrees"),
            sample_refractive_index=details.get("sample_refractive_index"),
        )
        condition_path = run_dir / "experimental-condition.json"
        if condition_path.is_file():
            retained_condition = json.loads(condition_path.read_text())
            if retained_condition != condition.as_dict():
                raise RuntimeError(
                    "The retained experimental condition contract differs from the requested "
                    "sample or measurement settings; start a new spectrum run"
                )
        else:
            write_experimental_condition(condition_path, condition)
        sampled_grid, sampled_intensity, sampled_successful = build_ir_spectrum(
            records, scale_factor=scale_factor, fwhm=profile.residual_fwhm_cm_1,
        )
        diagnostic_csv = None
        diagnostic_png = None
        if environment_sampling:
            display_records, band_statistics = collapse_insufficient_environment_bands(
                sampled_successful, scale_factor=scale_factor,
            )
            band_statistics = apply_environment_convergence_status(
                band_statistics, environment_convergence,
            )
            ensemble_grid, ensemble_intensity, ensemble_successful = build_ir_spectrum(
                display_records, scale_factor=scale_factor, fwhm=profile.residual_fwhm_cm_1,
            )
            diagnostic_display, diagnostic_column = format_spectrum(
                sampled_intensity, style=spectrum_style, max_absorbance=max_absorbance,
            )
            diagnostic_csv = write_spectrum_csv(
                run_dir / "spectrum_environment_sampled.csv", sampled_grid,
                diagnostic_display, diagnostic_column,
            )
            diagnostic_png = write_spectrum_plot(
                run_dir / "spectrum_environment_sampled.png", sampled_grid,
                diagnostic_display, diagnostic_column.replace("_", " "),
                title="Sampled environment distribution (diagnostic)",
            )
        else:
            ensemble_grid, ensemble_intensity, ensemble_successful = (
                sampled_grid, sampled_intensity, sampled_successful
            )
            band_statistics = ensemble_band_statistics(
                ensemble_successful, scale_factor=scale_factor,
            )
        measured_intensity = apply_measurement_response(
            ensemble_grid, ensemble_intensity, profile=profile,
        )
        raw_csv = write_spectrum_csv(run_dir / "spectrum_raw.csv", raw_grid, raw_display, column_name)
        raw_png = write_spectrum_plot(
            run_dir / "spectrum_raw.png", raw_grid, raw_display, column_name.replace("_", " "),
            title="Harmonic conformer spectrum",
        )
        ensemble_display, ensemble_column = format_spectrum(
            ensemble_intensity, style=spectrum_style, max_absorbance=max_absorbance,
        )
        intrinsic_csv = write_spectrum_csv(
            run_dir / "spectrum_intrinsic.csv", ensemble_grid, ensemble_intensity,
            "absorption_strength",
        )
        intrinsic_png = write_spectrum_plot(
            run_dir / "spectrum_intrinsic.png", ensemble_grid,
            ensemble_intensity / max(float(np.max(ensemble_intensity)), 1e-300),
            "relative absorption strength", title="Intrinsic calculated IR spectrum",
        )
        ensemble_csv = write_spectrum_csv(
            run_dir / "spectrum_ensemble.csv", ensemble_grid, ensemble_display, ensemble_column,
        )
        ensemble_png = write_spectrum_plot(
            run_dir / "spectrum_ensemble.png", ensemble_grid, ensemble_display,
            ensemble_column.replace("_", " "), title="Calculated conformer-ensemble IR spectrum",
        )
        measured_display, measured_column = format_spectrum(
            measured_intensity, style=spectrum_style, max_absorbance=max_absorbance,
        )
        csv_path = write_spectrum_csv(run_dir / "spectrum.csv", ensemble_grid, measured_display, measured_column)
        png_path = write_spectrum_plot(
            run_dir / "spectrum.png", ensemble_grid, measured_display,
            measured_column.replace("_", " "), title="Simulated FTIR measurement",
        )
        for band in band_statistics:
            band["residual_fwhm_cm-1"] = profile.residual_fwhm_cm_1
            band["instrument_resolution_cm-1"] = profile.instrument_resolution_cm_1
        statistics_path = run_dir / "ensemble-bands.json"
        statistics_path.write_text(json.dumps({
            "schema_version": 3 if environment_sampling else 2,
            "matching_basis": sorted({band["matching_basis"] for band in band_statistics}),
            "environment_sufficiency_applied": environment_sampling,
            "bands": band_statistics,
        }, indent=2, sort_keys=True) + "\n")
        mode_character_path = run_dir / "mode-character.json"
        mode_character_path.write_text(json.dumps({
            "schema_version": 1,
            "kind": "internal_coordinate_band_character",
            "bands": [
                {
                    "mode": band["mode"],
                    "matching_basis": band["matching_basis"],
                    "matching_confidence": band["matching_confidence"],
                    "mode_character": band.get("mode_character"),
                }
                for band in band_statistics
            ],
        }, indent=2, sort_keys=True) + "\n")
        environment_sampling_path = None
        environment_summary = None
        if environment_sampling:
            statuses = sorted({band["width_status"] for band in band_statistics})
            environment_summary = {
                "width_status": (
                    statuses[0] if len(statuses) == 1 else "mixed_band_sufficiency"
                ),
                "calculated_environment_fwhm_cm-1": (
                    0.0 if statuses == ["insufficient_environment_sampling"] else None
                ),
                "display_width_source": (
                    "residual_plus_instrument"
                    if statuses == ["insufficient_environment_sampling"]
                    else "per_band_sufficiency_decision"
                ),
                "bands_evaluated": len(band_statistics),
                "bands_with_calculated_environment_width": sum(
                    band["width_status"] == "calculated_environment_width"
                    for band in band_statistics
                ),
                "bands_with_unconverged_environment_width": sum(
                    band["width_status"] == "environment_width_unconverged"
                    for band in band_statistics
                ),
                "bands_with_insufficient_environment_sampling": sum(
                    band["width_status"] == "insufficient_environment_sampling"
                    for band in band_statistics
                ),
                "convergence_status": (
                    (environment_convergence or {}).get("summary", {}).get(
                        "status", "adaptive_convergence_not_evaluated",
                    )
                ),
                "converged": bool(
                    (environment_convergence or {}).get("summary", {}).get("converged", False)
                ),
            }
            environment_sampling_path = run_dir / "environment-sampling.json"
            environment_report = (
                json.loads(environment_sampling_path.read_text())
                if environment_sampling_path.is_file() else {}
            )
            width_sufficiency = {
                "thresholds": environment_sufficiency_metadata(),
                "summary": environment_summary,
                "bands": band_statistics,
                "diagnostic_sampled_spectrum": str(diagnostic_csv),
                "display_policy": (
                    "Bands failing sufficiency are collapsed to their weighted center before "
                    "residual and instrument response are applied."
                ),
            }
            environment_report.update({
                "schema_version": max(2, int(environment_report.get("schema_version", 0))),
                "kind": "environment_sampling_and_width_sufficiency",
                # Retain these aliases for readers of the first gate schema.
                "thresholds": width_sufficiency["thresholds"],
                "summary": width_sufficiency["summary"],
                "bands": width_sufficiency["bands"],
                "width_sufficiency": width_sufficiency,
            })
            environment_sampling_path.write_text(
                json.dumps(environment_report, indent=2, sort_keys=True) + "\n"
            )
        return {
            "spectrum_csv": csv_path,
            "spectrum_png": png_path,
            "raw_spectrum_csv": raw_csv,
            "raw_spectrum_png": raw_png,
            "ensemble_spectrum_csv": ensemble_csv,
            "ensemble_spectrum_png": ensemble_png,
            "intrinsic_spectrum_csv": intrinsic_csv,
            "intrinsic_spectrum_png": intrinsic_png,
            "experimental_condition": condition_path,
            "experimental_condition_contract": condition.as_dict(),
            "ensemble_band_statistics": statistics_path,
            "mode_character": mode_character_path,
            "environment_sampling": environment_sampling_path,
            "environment_sampling_summary": environment_summary,
            "sampled_environment_spectrum_csv": diagnostic_csv,
            "sampled_environment_spectrum_png": diagnostic_png,
            "successful": sampled_successful,
            "practical_profile": None,
            "experimental_profile": experimental_profile_metadata(profile),
        }

    if not practical_smiles:
        raise ValueError("A SMILES string is required for the practical spectrum model")
    from .practical_ir import describe_molecule, practical_band_transform, practical_profile_metadata

    features = describe_molecule(practical_smiles)
    def transform(frequency: float, base_width: float) -> tuple[float, float]:
        adjusted_frequency, adjusted_width, _ = practical_band_transform(
            frequency, baseline_fwhm=base_width, features=features,
        )
        return adjusted_frequency, adjusted_width

    practical_grid, practical_intensity, _ = build_ir_spectrum(
        records, scale_factor=scale_factor, fwhm=fwhm, band_transform=transform,
    )
    practical_display, practical_column = format_spectrum(practical_intensity, style=spectrum_style, max_absorbance=max_absorbance)
    raw_csv = write_spectrum_csv(run_dir / "spectrum_raw.csv", raw_grid, raw_display, column_name)
    raw_png = write_spectrum_plot(run_dir / "spectrum_raw.png", raw_grid, raw_display, column_name.replace("_", " "))
    csv_path = write_spectrum_csv(run_dir / "spectrum.csv", practical_grid, practical_display, practical_column)
    png_path = write_spectrum_plot(
        run_dir / "spectrum.png", practical_grid, practical_display,
        practical_column.replace("_", " "), title="Rule-calibrated practical IR spectrum",
    )
    return {
        "spectrum_csv": csv_path,
        "spectrum_png": png_path,
        "raw_spectrum_csv": raw_csv,
        "raw_spectrum_png": raw_png,
        "successful": successful,
        "practical_profile": practical_profile_metadata(features),
        "experimental_profile": None,
    }


def assemble_self_dimer_environment(
    monomer_artifacts: dict,
    cluster_records: list[dict],
    run_dir: Path,
    *,
    scale_factor: float,
    fwhm: float,
    spectrum_style: str,
    max_absorbance: float,
    phase: str,
    measurement: str,
    instrument_resolution: float,
    apodization: str,
    residual_fwhm: float | None,
    environment_convergence: dict | None = None,
    experimental_condition_details: dict | None = None,
) -> dict:
    """Replace the displayed neat-phase model with per-molecule dimer spectra."""
    from .cluster_ir import normalize_environment_cluster_records

    run_dir = Path(run_dir)
    retained = {
        "spectrum_csv": "spectrum_monomer.csv",
        "spectrum_png": "spectrum_monomer.png",
        "raw_spectrum_csv": "spectrum_monomer_raw.csv",
        "raw_spectrum_png": "spectrum_monomer_raw.png",
        "ensemble_spectrum_csv": "spectrum_monomer_ensemble.csv",
        "ensemble_spectrum_png": "spectrum_monomer_ensemble.png",
        "ensemble_band_statistics": "monomer-ensemble-bands.json",
        "intrinsic_spectrum_csv": "spectrum_monomer_intrinsic.csv",
        "intrinsic_spectrum_png": "spectrum_monomer_intrinsic.png",
    }
    monomer_retained: dict[str, Path] = {}
    for key, filename in retained.items():
        source = monomer_artifacts.get(key)
        if source and Path(source).is_file():
            destination = run_dir / filename
            shutil.copy2(source, destination)
            monomer_retained[key] = destination

    reliability_records = [
        {
            "independent_environment_id": record.get("independent_environment_id"),
            "topology": record.get("topology"),
            **(record.get("snapshot_hessian_reliability") or {}),
        }
        for record in cluster_records
    ]
    reliability_summary = {
        "representatives_evaluated": len(reliability_records),
        "poor_stationarity_representatives": sum(
            item.get("stationarity_status") == "poor" for item in reliability_records
        ),
        "full_hessian_policy": "diagnostic_only_if_material_imaginary_modes_are_present",
        "representatives": reliability_records,
    }
    sampling_path = run_dir / "environment-sampling.json"
    if sampling_path.is_file():
        sampling_report = json.loads(sampling_path.read_text())
        sampling_report["snapshot_hessian_reliability"] = reliability_summary
        sampling_path.write_text(json.dumps(sampling_report, indent=2, sort_keys=True) + "\n")
    normalized = normalize_environment_cluster_records(cluster_records)
    artifacts = _write_ir_artifacts(
        normalized, run_dir, scale_factor=scale_factor, fwhm=fwhm,
        spectrum_style=spectrum_style, max_absorbance=max_absorbance,
        spectrum_model="experimental", practical_smiles=None,
        phase=phase, measurement=measurement,
        instrument_resolution=instrument_resolution, apodization=apodization,
        residual_fwhm=residual_fwhm,
        environment_sampling=True,
        environment_convergence=environment_convergence,
        experimental_condition_details=experimental_condition_details,
    )
    dimer_csv = run_dir / "spectrum_dimer.csv"
    network_csv = run_dir / "spectrum_network.csv"
    dimer_records = [record for record in normalized if int(record.get("cluster_size") or 2) == 2]
    if dimer_records:
        from .experimental_ir import resolve_experimental_profile
        dimer_profile = resolve_experimental_profile(
            phase=phase, measurement=measurement,
            instrument_resolution_cm_1=instrument_resolution, apodization=apodization,
            residual_fwhm_cm_1=residual_fwhm,
        )
        dimer_grid, dimer_intensity, _ = build_ir_spectrum(
            dimer_records, scale_factor=scale_factor, fwhm=dimer_profile.residual_fwhm_cm_1,
        )
        dimer_display, dimer_column = format_spectrum(
            dimer_intensity, style=spectrum_style, max_absorbance=max_absorbance,
        )
        write_spectrum_csv(dimer_csv, dimer_grid, dimer_display, dimer_column)
    else:
        shutil.copy2(artifacts["ensemble_spectrum_csv"], dimer_csv)
    shutil.copy2(artifacts["ensemble_spectrum_csv"], network_csv)
    artifacts["dimer_spectrum_csv"] = dimer_csv
    artifacts["network_spectrum_csv"] = network_csv
    profile = dict(artifacts["experimental_profile"] or {})
    sampling_summary = artifacts.get("environment_sampling_summary") or {}
    convergence_path = run_dir / "environment-convergence.json"
    convergence_summary = (environment_convergence or {}).get("summary") or {}
    artifacts["environment_convergence"] = convergence_path if convergence_path.is_file() else None
    cluster_population_models = sorted({
        str(record["population_model"]) for record in normalized
        if record.get("population_model")
    })
    cluster_population_warnings = sorted({
        str(record["population_warning"]) for record in normalized
        if record.get("population_warning")
    })
    profile.update({
        "environment_source": (
            "sampled_neutral_dimer_trimer_network_configurations"
            if any(int(record.get("cluster_size") or 2) == 3 for record in normalized)
            else "sampled_neutral_self_dimer_configurations"
        ),
        "environment_width_status": sampling_summary.get("width_status"),
        "environment_sampling_artifact": str(artifacts.get("environment_sampling") or ""),
        "environment_convergence_status": convergence_summary.get("status"),
        "environment_convergence_artifact": str(artifacts.get("environment_convergence") or ""),
        "intensity_normalization": "cluster_intensity_divided_by_target_molecule_count",
        "cluster_population_model": (
            cluster_population_models[0]
            if len(cluster_population_models) == 1 else cluster_population_models
        ),
        "population_warning": (
            cluster_population_warnings[0]
            if len(cluster_population_warnings) == 1 else cluster_population_warnings
        ),
        "network_composition": (
            "dimers_plus_trimers" if any(int(record.get("cluster_size") or 2) == 3 for record in normalized)
            else "dimer_only"
        ),
        "snapshot_hessian_reliability": reliability_summary,
        "limitations": (
            "Neutral dimer/trimer snapshot Hessians use the recorded conditional sampling weights, not exact "
            "bulk populations. Explicit solvent, larger aggregates, anharmonic lifetimes, crystal packing, "
            "and vendor-specific instrument response are not represented."
        ),
    })
    artifacts["experimental_profile"] = profile
    write_metadata(
        run_dir,
        experimental_ir_profile=profile,
        environment_model=(
            "neutral_dimer_trimer_network_ensemble"
            if any(int(record.get("cluster_size") or 2) == 3 for record in normalized)
            else "neutral_self_dimer_ensemble"
        ),
        environment_width_status=sampling_summary.get("width_status"),
        environment_sampling_artifact=str(artifacts.get("environment_sampling") or ""),
        environment_convergence_status=convergence_summary.get("status"),
        environment_converged=bool(convergence_summary.get("converged", False)),
        environment_convergence_artifact=str(artifacts.get("environment_convergence") or ""),
        snapshot_hessian_reliability=reliability_summary,
        cluster_conformers_successful=len(artifacts["successful"]),
        monomer_spectrum_csv=str(monomer_retained.get("spectrum_csv", "")),
        spectrum_csv=str(artifacts["spectrum_csv"]),
        spectrum_png=str(artifacts["spectrum_png"]),
        spectrum_intrinsic_csv=str(artifacts.get("intrinsic_spectrum_csv") or ""),
        spectrum_intrinsic_png=str(artifacts.get("intrinsic_spectrum_png") or ""),
        mode_character_artifact=str(artifacts.get("mode_character") or ""),
        spectrum_dimer_csv=str(dimer_csv),
        spectrum_network_csv=str(network_csv),
        status="completed",
    )
    return {
        **{key: value for key, value in artifacts.items() if key not in {"successful", "practical_profile"}},
        "monomer_artifacts": monomer_retained,
        "conformers": normalized,
    }


def finalize_ir_spectrum(
    run_dir: Path,
    *,
    temperature: float | None = None,
    scale_factor: float | None = None,
    fwhm: float | None = None,
    spectrum_style: str | None = None,
    max_absorbance: float | None = None,
    spectrum_model: str | None = None,
    practical_smiles: str | None = None,
    phase: str | None = None,
    measurement: str | None = None,
    instrument_resolution: float | None = None,
    apodization: str | None = None,
    residual_fwhm: float | None = None,
) -> dict:
    """Build spectrum artifacts from completed conformer directories.

    This is intended for recovery after an interrupted workflow. Incomplete
    ORCA calculations are recorded but never used in the weighted spectrum.
    """
    run_dir = Path(run_dir)
    metadata_path = run_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    temperature = float(metadata.get("temperature", 298.15) if temperature is None else temperature)
    scale_factor = float(metadata.get("scale_factor", 0.97) if scale_factor is None else scale_factor)
    fwhm = float(metadata.get("fwhm_cm_1", 15.0) if fwhm is None else fwhm)
    spectrum_style = str(metadata.get("spectrum_style", "transmittance") if spectrum_style is None else spectrum_style)
    max_absorbance = float(metadata.get("max_absorbance", 1.0) if max_absorbance is None else max_absorbance)
    spectrum_model = str(metadata.get("spectrum_model", "raw") if spectrum_model is None else spectrum_model)
    practical_smiles = practical_smiles or metadata.get("smiles")
    phase = str(metadata.get("phase", "liquid") if phase is None else phase)
    measurement = str(metadata.get("measurement", "auto") if measurement is None else measurement)
    instrument_resolution = float(metadata.get("instrument_resolution_cm_1", 4.0)
                                  if instrument_resolution is None else instrument_resolution)
    apodization = str(metadata.get("apodization", "happ-genzel") if apodization is None else apodization)
    if residual_fwhm is None:
        residual_fwhm = metadata.get("residual_fwhm_cm_1")
    experimental_condition_details = dict(metadata.get("experimental_condition_details") or {})
    records = load_completed_ir_records(run_dir, temperature=temperature)

    summary_path = run_dir / "conformers.json"
    explicit_population_models = sorted({
        str(record["population_model"]) for record in records
        if record.get("population_weight") is not None and record.get("population_model")
    })
    weighting_basis = (
        "explicit_population_weights" if explicit_population_models
        else "final_optimized_electronic_energy_hartree"
    )
    try:
        if spectrum_model == "practical" and not practical_smiles:
            raise ValueError("The retained run has no SMILES needed by the practical spectrum model")
        artifacts = _write_ir_artifacts(
            records, run_dir, scale_factor=scale_factor, fwhm=fwhm,
            spectrum_style=spectrum_style, max_absorbance=max_absorbance,
            spectrum_model=spectrum_model, practical_smiles=practical_smiles,
            phase=phase, measurement=measurement, instrument_resolution=instrument_resolution,
            apodization=apodization, residual_fwhm=residual_fwhm,
            temperature_K=temperature,
            experimental_condition_details=experimental_condition_details,
        )
    except Exception:
        _write_conformer_summary(summary_path, records)
        failure_metadata = {
            "last_operation": "ir_spectrum_finalize",
            "status": "failed",
            "conformers_discovered": len(records),
            "conformers_successful": 0,
        }
        if "workflow" not in metadata:
            failure_metadata["workflow"] = "ir_spectrum"
        write_metadata(run_dir, **failure_metadata)
        raise
    _write_conformer_summary(summary_path, records)
    completion_metadata = dict(
        last_operation="ir_spectrum_finalize",
        temperature=temperature,
        conformer_weighting_energy_basis=weighting_basis,
        explicit_population_models=explicit_population_models,
        conformer_degeneracy_model="goat_reported_degeneracy_or_one",
        line_shape_model="unit_area_gaussian_v2",
        scale_factor=scale_factor,
        fwhm_cm_1=fwhm,
        spectrum_style=spectrum_style,
        max_absorbance=max_absorbance,
        spectrum_model=spectrum_model,
        practical_rule_profile=artifacts["practical_profile"],
        experimental_ir_profile=artifacts["experimental_profile"],
        experimental_condition_contract=artifacts.get("experimental_condition_contract"),
        experimental_condition_artifact=str(artifacts.get("experimental_condition") or ""),
        phase=phase,
        measurement=measurement,
        instrument_resolution_cm_1=instrument_resolution,
        apodization=apodization,
        residual_fwhm_cm_1=(artifacts["experimental_profile"] or {}).get("residual_fwhm_cm_1"),
        conformers_discovered=len(records),
        conformers_successful=len(artifacts["successful"]),
        spectrum_csv=str(artifacts["spectrum_csv"]),
        spectrum_png=str(artifacts["spectrum_png"]),
        spectrum_intrinsic_csv=str(artifacts.get("intrinsic_spectrum_csv") or ""),
        spectrum_intrinsic_png=str(artifacts.get("intrinsic_spectrum_png") or ""),
        mode_character_artifact=str(artifacts.get("mode_character") or ""),
        status="completed",
    )
    if "workflow" not in metadata:
        completion_metadata["workflow"] = "ir_spectrum"
    write_metadata(run_dir, **completion_metadata)
    result = {
        **{key: value for key, value in artifacts.items() if key != "successful" and key != "practical_profile"},
        "conformer_summary": summary_path,
        "conformers": records,
    }
    cluster_dir = run_dir / "clusters"
    if spectrum_model == "experimental" and (cluster_dir / "selected-conformers").is_dir():
        cluster_records = load_completed_ir_records(cluster_dir, temperature=temperature)
        _write_conformer_summary(cluster_dir / "conformers.json", cluster_records)
        from .environment_convergence import convergence_schedule, evaluate_completed_environment_batches
        cluster_manifest_path = cluster_dir / "selected-conformers.json"
        cluster_manifest = (
            json.loads(cluster_manifest_path.read_text()) if cluster_manifest_path.is_file() else {}
        )
        representative_budget = len(cluster_manifest.get("conformers", [])) or len(cluster_records)
        fidelity = str(metadata.get("fidelity", "auto"))
        attempted = len(cluster_records)
        endpoints = [
            endpoint for endpoint in convergence_schedule(fidelity, representative_budget)
            if endpoint <= attempted
        ]
        if attempted and attempted not in endpoints:
            endpoints.append(attempted)
        convergence_report = evaluate_completed_environment_batches(
            [(endpoint, cluster_records[:endpoint]) for endpoint in sorted(set(endpoints))],
            run_dir, fidelity=fidelity, scale_factor=scale_factor,
            representative_budget=representative_budget,
            budget_exhausted=attempted >= representative_budget,
        )
        result = assemble_self_dimer_environment(
            result, cluster_records, run_dir,
            scale_factor=scale_factor, fwhm=fwhm,
            spectrum_style=spectrum_style, max_absorbance=max_absorbance,
            phase=phase, measurement=measurement,
            instrument_resolution=instrument_resolution, apodization=apodization,
            residual_fwhm=residual_fwhm,
            environment_convergence=convergence_report,
            experimental_condition_details=experimental_condition_details,
        )
    return result


def run_ir_spectrum(
    conformer_xyzs: Iterable[Path],
    run_dir: Path,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    ncores: int = 1,
    temperature: float = 298.15,
    scale_factor: float = 0.97,
    fwhm: float = 15.0,
    spectrum_style: str = "transmittance",
    max_absorbance: float = 1.0,
    spectrum_model: str = "raw",
    practical_smiles: str | None = None,
    phase: str = "liquid",
    measurement: str = "auto",
    instrument_resolution: float = 4.0,
    apodization: str = "happ-genzel",
    residual_fwhm: float | None = None,
    experimental_condition_details: dict | None = None,
    method_keywords: list[str] | None = None,
    geometry_role: str = "optimized_minimum",
    progress: Callable[[str], None] | None = None,
) -> dict:
    """Optimize conformers, calculate frequencies, and write an IR spectrum."""
    if geometry_role not in {"optimized_minimum", "environment_snapshot"}:
        raise ValueError("Unsupported IR geometry role")
    run_dir = Path(run_dir)
    metadata_path = run_dir / "metadata.json"
    retained_metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    selection_path = run_dir / "selected-conformers.json"
    selection = json.loads(selection_path.read_text()) if selection_path.is_file() else {}
    selection_by_position = {item["selected_position"]: item for item in selection.get("conformers", [])}
    xyz_paths = list(conformer_xyzs)
    records: list[dict] = []
    for index, xyz_path in enumerate(xyz_paths, start=1):
        job_dir = run_dir / "conformers" / f"conf-{index:03d}"
        job_dir.mkdir(parents=True, exist_ok=True)
        input_xyz = job_dir / "input.xyz"
        source_xyz = Path(xyz_path).read_text()
        previous_xyz = input_xyz.read_text() if input_xyz.is_file() else None
        requested_contract = _calculation_contract(
            source_xyz, charge=charge, multiplicity=multiplicity, method_keywords=method_keywords,
            geometry_role=geometry_role,
        )
        contract_path = job_dir / "calculation-contract.json"
        retained_contract = json.loads(contract_path.read_text()) if contract_path.is_file() else None
        base_contract_matches = bool(
            retained_contract
            and all(
                retained_contract.get(key, "optimized_minimum" if key == "geometry_role" else None) == value
                for key, value in requested_contract.items()
            )
        )
        if retained_contract is None:
            base_contract_matches = _legacy_contract_matches(
                retained_metadata, previous_xyz, source_xyz, requested_contract,
            )
        opt_out = job_dir / "opt.out"
        freq_out = job_dir / "freq.out"
        if (_completed_orca_output(opt_out) or _completed_orca_output(freq_out)) and not base_contract_matches:
            raise RuntimeError(
                f"Conformer {index} has completed ORCA artifacts for a different or unverified geometry/electronic contract; "
                "start a new spectrum run rather than mixing calculations"
            )
        input_xyz.write_text(source_xyz)
        source = selection_by_position.get(index, {})
        record = {"index": index, "temperature": temperature,
                  "degeneracy": source.get("goat_degeneracy", 1),
                  "source_conformer": source or None,
                  "geometry_role": geometry_role}
        for key in (
            "population_weight", "population_model", "population_warning",
            "geometry_cluster_id", "independent_environment_id", "environment_features",
            "sampling_support_count", "source_xtb_candidate_id",
            "cluster_size", "topology", "molecule_atom_ranges", "local_stretch_bonds",
            "hydrogen_bond_interactions",
            "environment_refinement", "environment_refinement_artifact",
            "vibrational_route", "gradient_rms_hartree_per_bohr",
            "gradient_maximum_component_hartree_per_bohr",
        ):
            if key in source:
                record[key] = source[key]
        _backfill_local_environment_roles(record)
        try:
            retained_opt_hash = retained_contract.get("optimized_xyz_sha256") if retained_contract else None
            if geometry_role == "environment_snapshot":
                optimized_xyz = input_xyz
                optimized_xyz_hash = _sha256_text(source_xyz)
                optimization_reused = False
                record["optimized_xyz"] = optimized_xyz
                record["stationarity_model"] = "snapshot_hessian_without_unconstrained_dft_optimization"
            else:
                optimized_xyz = job_dir / "opt.xyz"
                current_opt_hash = _sha256_text(optimized_xyz.read_text()) if optimized_xyz.is_file() else None
                geometry_contract_matches = base_contract_matches and (
                    retained_opt_hash == current_opt_hash if retained_opt_hash is not None else current_opt_hash is not None
                )
                if _completed_optimization(opt_out) and optimized_xyz.is_file() and geometry_contract_matches:
                    optimization_reused = True
                    if progress:
                        progress(f"Conformer {index}/{len(xyz_paths)}: reusing completed optimization")
                else:
                    optimization_reused = False
                    if progress:
                        progress(f"Conformer {index}/{len(xyz_paths)}: optimizing")
                    opt_input = create_orca_input(input_xyz, charge, multiplicity, opt=True, label="opt", ncores=ncores, method_keywords=method_keywords)
                    opt = run_orca(opt_input)
                    optimized_xyz = opt["xyz"]
                if not optimized_xyz.exists() or not _completed_optimization(opt_out):
                    raise RuntimeError("ORCA did not produce a normally terminated, converged optimized geometry")
                record["optimized_xyz"] = optimized_xyz
                record["energy"] = parse_orca_energy(opt_out)
                optimized_xyz_hash = _sha256_text(optimized_xyz.read_text())
            frequency_contract_matches = base_contract_matches and (
                retained_opt_hash == optimized_xyz_hash
                if retained_opt_hash is not None else optimization_reused
            )
            if frequency_contract_matches and _completed_orca_output(freq_out):
                if progress:
                    progress(f"Conformer {index}/{len(xyz_paths)}: reusing completed frequencies")
            else:
                if progress:
                    progress(f"Conformer {index}/{len(xyz_paths)}: calculating frequencies")
                freq_input = create_orca_input(optimized_xyz, charge, multiplicity, freq=True, label="freq", ncores=ncores, method_keywords=method_keywords)
                run_orca(freq_input)
            if geometry_role == "environment_snapshot":
                record["energy"] = parse_orca_energy(freq_out)
            record["frequency_output"] = freq_out
            record["frequency_check"] = frequency_stability_check(freq_out)
            if geometry_role == "environment_snapshot":
                record["snapshot_hessian_reliability"] = _snapshot_hessian_reliability(record["frequency_check"])
                refinement = record.get("environment_refinement") or {}
                gradient = refinement.get("gradient") or {}
                if refinement:
                    frequency_reliability = record["snapshot_hessian_reliability"]
                    gradient_permitted = refinement.get("full_hessian_use") == "permitted"
                    frequency_permitted = frequency_reliability.get("stationarity_status") != "poor"
                    record["snapshot_hessian_reliability"] = {
                        **frequency_reliability,
                        "gradient_rms_hartree_per_bohr": gradient.get(
                            "gradient_rms_hartree_per_bohr"
                        ),
                        "gradient_maximum_component_hartree_per_bohr": gradient.get(
                            "gradient_maximum_component_hartree_per_bohr"
                        ),
                        "gradient_status": refinement.get("stationarity_status"),
                        "full_hessian_use": (
                            "permitted" if gradient_permitted and frequency_permitted
                            else "diagnostic_only"
                        ),
                        "stationarity_status": (
                            "usable" if gradient_permitted and frequency_permitted else "poor"
                        ),
                    }
            record["ir_modes"] = parse_orca_ir(freq_out)
            if (not record["frequency_check"].get("IsMinimum")
                    and geometry_role != "environment_snapshot"):
                raise RuntimeError("Frequency calculation does not describe a local minimum")
            if geometry_role == "environment_snapshot":
                record["snapshot_frequency_warning"] = (
                    "The Hessian was evaluated at a restrained xTB environment snapshot; "
                    "it is not asserted to be an unconstrained DFT stationary point."
                )
            contract_path.write_text(json.dumps({
                **requested_contract,
                "optimized_xyz_sha256": optimized_xyz_hash,
            }, indent=2, sort_keys=True) + "\n")
        except Exception as error:  # Retain failed job artifacts and report them in the manifest.
            record["error"] = str(error)
            if progress:
                progress(f"Conformer {index}/{len(xyz_paths)}: failed ({error})")
        records.append(record)

    summary_path = run_dir / "conformers.json"
    explicit_population_models = sorted({
        str(record["population_model"]) for record in records
        if record.get("population_weight") is not None and record.get("population_model")
    })
    weighting_basis = (
        "explicit_population_weights" if explicit_population_models
        else "final_optimized_electronic_energy_hartree"
    )
    try:
        artifacts = _write_ir_artifacts(
            records, run_dir, scale_factor=scale_factor, fwhm=fwhm,
            spectrum_style=spectrum_style, max_absorbance=max_absorbance,
            spectrum_model=spectrum_model, practical_smiles=practical_smiles,
            phase=phase, measurement=measurement, instrument_resolution=instrument_resolution,
            apodization=apodization, residual_fwhm=residual_fwhm,
            temperature_K=temperature,
            experimental_condition_details=experimental_condition_details,
        )
    except Exception:
        # The per-conformer record is still the primary diagnostic artifact
        # when all calculations fail or no usable IR modes were parsed.
        _write_conformer_summary(summary_path, records)
        write_metadata(
            run_dir,
            workflow="ir_spectrum",
            charge=charge,
            multiplicity=multiplicity,
            ncores=ncores,
            temperature=temperature,
            conformer_weighting_energy_basis=weighting_basis,
            explicit_population_models=explicit_population_models,
            conformer_degeneracy_model="goat_reported_degeneracy_or_one",
            line_shape_model="unit_area_gaussian_v2",
            scale_factor=scale_factor,
            fwhm_cm_1=fwhm,
            spectrum_style=spectrum_style,
            max_absorbance=max_absorbance,
            spectrum_model=spectrum_model,
            conformers_attempted=len(records),
            conformers_successful=0,
            status="failed",
        )
        raise
    _write_conformer_summary(summary_path, records)
    write_metadata(
        run_dir,
        workflow="ir_spectrum",
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        temperature=temperature,
        conformer_weighting_energy_basis=weighting_basis,
        explicit_population_models=explicit_population_models,
        conformer_degeneracy_model="goat_reported_degeneracy_or_one",
        line_shape_model="unit_area_gaussian_v2",
        scale_factor=scale_factor,
        fwhm_cm_1=fwhm,
        spectrum_style=spectrum_style,
        max_absorbance=max_absorbance,
        spectrum_model=spectrum_model,
        practical_rule_profile=artifacts["practical_profile"],
        experimental_ir_profile=artifacts["experimental_profile"],
        experimental_condition_contract=artifacts.get("experimental_condition_contract"),
        experimental_condition_artifact=str(artifacts.get("experimental_condition") or ""),
        experimental_condition_details=experimental_condition_details or {},
        phase=phase,
        measurement=measurement,
        instrument_resolution_cm_1=instrument_resolution,
        apodization=apodization,
        residual_fwhm_cm_1=(artifacts["experimental_profile"] or {}).get("residual_fwhm_cm_1"),
        conformers_attempted=len(records),
        conformers_successful=len(artifacts["successful"]),
        spectrum_csv=str(artifacts["spectrum_csv"]),
        spectrum_png=str(artifacts["spectrum_png"]),
        spectrum_intrinsic_csv=str(artifacts.get("intrinsic_spectrum_csv") or ""),
        spectrum_intrinsic_png=str(artifacts.get("intrinsic_spectrum_png") or ""),
        mode_character_artifact=str(artifacts.get("mode_character") or ""),
        status="completed",
    )
    return {
        **{key: value for key, value in artifacts.items() if key != "successful" and key != "practical_profile"},
        "conformer_summary": summary_path,
        "conformers": records,
    }


def resume_ir_spectrum(run_dir: Path, **settings) -> dict:
    """Continue an interrupted spectrum using its retained run contract."""
    run_dir = Path(run_dir)
    xyz_paths = sorted((run_dir / "selected-conformers").glob("*.xyz"))
    if not xyz_paths:
        xyz_paths = sorted((run_dir / "initial-conformers").glob("*.xyz"))
    if not xyz_paths:
        raise FileNotFoundError("No selected or initial conformer geometries found for resume")
    metadata_path = run_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    plan_path = run_dir / "spectrum-plan.json"
    calculation_plan = json.loads(plan_path.read_text()) if plan_path.is_file() else {}
    defaults = {
        "charge": metadata.get("charge", 0),
        "multiplicity": metadata.get("multiplicity", 1),
        "ncores": metadata.get("ncores", 1),
        "temperature": metadata.get("temperature", 298.15),
        "scale_factor": metadata.get("scale_factor", 0.97),
        "fwhm": metadata.get("fwhm_cm_1", 15.0),
        "spectrum_style": metadata.get("spectrum_style", "transmittance"),
        "max_absorbance": metadata.get("max_absorbance", 1.0),
        "spectrum_model": metadata.get("spectrum_model", "raw"),
        "practical_smiles": metadata.get("smiles"),
        "phase": metadata.get("phase", "liquid"),
        "measurement": metadata.get("measurement", "auto"),
        "instrument_resolution": metadata.get("instrument_resolution_cm_1", 4.0),
        "apodization": metadata.get("apodization", "happ-genzel"),
        "residual_fwhm": metadata.get("residual_fwhm_cm_1"),
        "experimental_condition_details": metadata.get("experimental_condition_details") or {},
        "method_keywords": (metadata.get("harmonic_method_profile") or {}).get("orca_keywords"),
        "local_mode_orca_invocations": metadata.get("local_mode_orca_invocations", 0),
    }
    resolved = {key: value for key, value in defaults.items() if value is not None}
    resolved.update({key: value for key, value in settings.items() if value is not None})
    sampling_plan = calculation_plan.get("xtb_environment_sampling") or {}
    sampling_smiles = metadata.get("smiles") or calculation_plan.get("smiles")
    if (
        resolved.get("spectrum_model") == "experimental"
        and sampling_plan.get("requested") is True
        and sampling_smiles
    ):
        from .xtb_sampling import sample_xtb_dimer_environments
        try:
            sampled_records, _ = sample_xtb_dimer_environments(
                sampling_smiles, run_dir,
                fidelity=str(sampling_plan.get("fidelity", metadata.get("fidelity", "auto"))),
                monomer_xyz=xyz_paths[0],
                charge=int(resolved.get("charge", 0)),
                multiplicity=int(resolved.get("multiplicity", 1)),
                ncores=int(resolved.get("ncores", 1)),
                candidate_count=int(sampling_plan.get("candidate_count", 40)),
                include_trimers=bool(sampling_plan.get("include_trimers", False)),
                progress=resolved.get("progress"),
            )
            # The restraints used to preserve diverse sampling strata must not
            # enter the force constants.  Evaluate a separate, unrestrained
            # numerical xTB Hessian at every retained snapshot.
            from .xtb_frequencies import sample_xtb_snapshot_frequencies
            sample_xtb_snapshot_frequencies(
                sampled_records, run_dir,
                charge=int(resolved.get("charge", 0)),
                multiplicity=int(resolved.get("multiplicity", 1)),
                ncores=int(resolved.get("ncores", 1)),
                progress=resolved.get("progress"),
            )
            retained_cluster_paths = sorted((run_dir / "clusters" / "selected-conformers").glob("*.xyz"))
            if not retained_cluster_paths:
                from .environment_selection import select_xtb_environment_representatives
                allocation = calculation_plan.get("orca_allocation") or {}
                representative_count = int(
                    sampling_plan.get(
                        "representative_orca_jobs_reserved",
                        allocation.get("reserved_environment_jobs", 0),
                    )
                )
                if representative_count > 0:
                    select_xtb_environment_representatives(
                        sampled_records, run_dir, representative_count=representative_count,
                    )
        except Exception as error:
            write_metadata(
                run_dir, xtb_environment_sampling_status="failed_during_resume",
                xtb_environment_sampling_error=str(error),
            )
    local_mode_orca_invocations = int(resolved.pop("local_mode_orca_invocations", 0))
    result = run_ir_spectrum(xyz_paths, run_dir, **resolved)
    cluster_dir = run_dir / "clusters"
    cluster_manifest_path = cluster_dir / "selected-conformers.json"
    if cluster_manifest_path.is_file():
        cluster_manifest = json.loads(cluster_manifest_path.read_text())
        cluster_paths = [
            Path(entry["xyz"]) for entry in sorted(
                cluster_manifest.get("conformers", []),
                key=lambda entry: int(entry["selected_position"]),
            )
        ]
    else:
        cluster_paths = sorted((cluster_dir / "selected-conformers").glob("*.xyz"))
    if resolved.get("spectrum_model") == "experimental" and cluster_paths:
        from .environment_convergence import run_adaptive_environment_convergence
        def _record_transfer_batch(endpoint: int, _batch_result: dict) -> list[Path]:
            from .environment_acquisition import (record_acquisition_batch,
                                                   reprioritize_pending_representatives)
            from .frequency_transfer import build_frequency_transfer_artifacts
            transfer = build_frequency_transfer_artifacts(run_dir)
            record_acquisition_batch(run_dir, endpoint, transfer)
            return reprioritize_pending_representatives(run_dir, endpoint)

        fidelity = str(sampling_plan.get("fidelity", metadata.get("fidelity", "auto")))
        cluster_result, convergence_report = run_adaptive_environment_convergence(
            cluster_paths, cluster_dir, run_dir, fidelity=fidelity,
            scale_factor=resolved.get("scale_factor", 0.97),
            run_batch=lambda paths: run_ir_spectrum(
                paths, cluster_dir, charge=0, multiplicity=1,
                ncores=resolved.get("ncores", 1), temperature=resolved.get("temperature", 298.15),
                scale_factor=resolved.get("scale_factor", 0.97), fwhm=resolved.get("fwhm", 15.0),
                spectrum_style=resolved.get("spectrum_style", "transmittance"),
                max_absorbance=resolved.get("max_absorbance", 1.0), spectrum_model="raw",
                method_keywords=resolved.get("method_keywords"), geometry_role="environment_snapshot",
                progress=resolved.get("progress"),
            ),
            progress=resolved.get("progress"),
            after_batch=_record_transfer_batch,
        )
        result = assemble_self_dimer_environment(
            result, cluster_result["conformers"], run_dir,
            scale_factor=resolved.get("scale_factor", 0.97), fwhm=resolved.get("fwhm", 15.0),
            spectrum_style=resolved.get("spectrum_style", "transmittance"),
            max_absorbance=resolved.get("max_absorbance", 1.0),
            phase=resolved.get("phase", "liquid"), measurement=resolved.get("measurement", "auto"),
            instrument_resolution=resolved.get("instrument_resolution", 4.0),
            apodization=resolved.get("apodization", "happ-genzel"),
            residual_fwhm=resolved.get("residual_fwhm"),
            environment_convergence=convergence_report,
            experimental_condition_details=resolved.get("experimental_condition_details"),
        )
        if local_mode_orca_invocations > 0:
            from .environment_local_modes import run_environment_local_mode_fallbacks
            local_modes = run_environment_local_mode_fallbacks(
                run_dir,
                maximum_orca_invocations=local_mode_orca_invocations,
                ncores=int(resolved.get("ncores", 1)),
                method_keywords=resolved.get("method_keywords"),
                progress=resolved.get("progress"),
            )
            result["environment_local_modes"] = local_modes.get("artifact")
            write_metadata(
                run_dir,
                environment_local_mode_status=local_modes.get("status"),
                environment_local_mode_artifact=local_modes.get("artifact"),
            )
        # Frequency transfer is deliberately fail-closed.  Until leave-one-
        # representative-out validation passes, the representative ORCA
        # ensemble remains the displayed spectrum and the larger xTB ensemble
        # is emitted only as a diagnostic.
        try:
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
            write_metadata(
                run_dir, frequency_transfer_status="failed",
                frequency_transfer_error=str(error),
            )
    return result
