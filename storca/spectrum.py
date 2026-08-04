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
from src.orca_runner import run_orca
from src.parser import attach_xyz_to_conformers, parse_goat_out, parse_orca_energy, parse_orca_ir, parse_xyz_ensemble
from src.stability.freq_check import frequency_stability_check

from .runs import write_metadata

HARTREE_BOLTZMANN_PER_K = 3.166811563e-6


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _calculation_contract(source_xyz: str, *, charge: int, multiplicity: int,
                          method_keywords: list[str] | None) -> dict:
    return {
        "schema_version": 1,
        "source_xyz_sha256": _sha256_text(source_xyz),
        "charge": int(charge),
        "multiplicity": int(multiplicity),
        "method_keywords": list(method_keywords) if method_keywords is not None else None,
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
    xtb_executable = shutil.which("xtb")
    if not xtb_executable:
        raise RuntimeError("xTB executable not found in PATH; activate the configured STORCA environment")
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
        if frequency_check is not None and frequency_check.get("IsMinimum") is not True:
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
            "degeneracy": conformer.get("degeneracy", 1),
            "source_conformer": conformer.get("source_conformer"),
            "optimized_xyz": str(conformer.get("optimized_xyz", "")),
            "frequency_output": str(conformer.get("frequency_output", "")),
            "frequency_check": conformer.get("frequency_check"),
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
        record = {"index": index, "temperature": temperature,
                  "degeneracy": source.get("goat_degeneracy", 1),
                  "source_conformer": source or None}
        if not _completed_optimization(opt_out):
            record["error"] = "Optimization is missing, unconverged, or did not terminate normally"
        elif not _completed_orca_output(freq_out):
            record["error"] = "Frequency calculation is missing or did not terminate normally"
        else:
            try:
                record["optimized_xyz"] = job_dir / "opt.xyz"
                record["energy"] = parse_orca_energy(opt_out)
                record["frequency_output"] = freq_out
                record["frequency_check"] = frequency_stability_check(freq_out)
                record["ir_modes"] = parse_orca_ir(freq_out)
                if not record["frequency_check"].get("IsMinimum"):
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
) -> dict:
    """Build raw artifacts and, when requested, an explicit practical display."""
    if spectrum_model not in {"raw", "practical"}:
        raise ValueError("spectrum_model must be 'raw' or 'practical'")
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
        return {"spectrum_csv": csv_path, "spectrum_png": png_path, "successful": successful, "practical_profile": None}

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
    records = load_completed_ir_records(run_dir, temperature=temperature)

    summary_path = run_dir / "conformers.json"
    try:
        if spectrum_model == "practical" and not practical_smiles:
            raise ValueError("The retained run has no SMILES needed by the practical spectrum model")
        artifacts = _write_ir_artifacts(
            records, run_dir, scale_factor=scale_factor, fwhm=fwhm,
            spectrum_style=spectrum_style, max_absorbance=max_absorbance,
            spectrum_model=spectrum_model, practical_smiles=practical_smiles,
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
        conformer_weighting_energy_basis="final_optimized_electronic_energy_hartree",
        conformer_degeneracy_model="goat_reported_degeneracy_or_one",
        line_shape_model="unit_area_gaussian_v2",
        scale_factor=scale_factor,
        fwhm_cm_1=fwhm,
        spectrum_style=spectrum_style,
        max_absorbance=max_absorbance,
        spectrum_model=spectrum_model,
        practical_rule_profile=artifacts["practical_profile"],
        conformers_discovered=len(records),
        conformers_successful=len(artifacts["successful"]),
        spectrum_csv=str(artifacts["spectrum_csv"]),
        spectrum_png=str(artifacts["spectrum_png"]),
        status="completed",
    )
    if "workflow" not in metadata:
        completion_metadata["workflow"] = "ir_spectrum"
    write_metadata(run_dir, **completion_metadata)
    return {
        **{key: value for key, value in artifacts.items() if key != "successful" and key != "practical_profile"},
        "conformer_summary": summary_path,
        "conformers": records,
    }


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
    method_keywords: list[str] | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict:
    """Optimize conformers, calculate frequencies, and write an IR spectrum."""
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
        )
        contract_path = job_dir / "calculation-contract.json"
        retained_contract = json.loads(contract_path.read_text()) if contract_path.is_file() else None
        base_contract_matches = bool(
            retained_contract
            and all(retained_contract.get(key) == value for key, value in requested_contract.items())
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
                  "source_conformer": source or None}
        try:
            optimized_xyz = job_dir / "opt.xyz"
            retained_opt_hash = retained_contract.get("optimized_xyz_sha256") if retained_contract else None
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
                retained_opt_hash == optimized_xyz_hash if retained_opt_hash is not None else optimization_reused
            )
            if optimization_reused and frequency_contract_matches and _completed_orca_output(freq_out):
                if progress:
                    progress(f"Conformer {index}/{len(xyz_paths)}: reusing completed frequencies")
            else:
                if progress:
                    progress(f"Conformer {index}/{len(xyz_paths)}: calculating frequencies")
                freq_input = create_orca_input(optimized_xyz, charge, multiplicity, freq=True, label="freq", ncores=ncores, method_keywords=method_keywords)
                run_orca(freq_input)
            record["frequency_output"] = freq_out
            record["frequency_check"] = frequency_stability_check(freq_out)
            record["ir_modes"] = parse_orca_ir(freq_out)
            if not record["frequency_check"].get("IsMinimum"):
                raise RuntimeError("Frequency calculation does not describe a local minimum")
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
    try:
        artifacts = _write_ir_artifacts(
            records, run_dir, scale_factor=scale_factor, fwhm=fwhm,
            spectrum_style=spectrum_style, max_absorbance=max_absorbance,
            spectrum_model=spectrum_model, practical_smiles=practical_smiles,
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
            conformer_weighting_energy_basis="final_optimized_electronic_energy_hartree",
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
        conformer_weighting_energy_basis="final_optimized_electronic_energy_hartree",
        conformer_degeneracy_model="goat_reported_degeneracy_or_one",
        line_shape_model="unit_area_gaussian_v2",
        scale_factor=scale_factor,
        fwhm_cm_1=fwhm,
        spectrum_style=spectrum_style,
        max_absorbance=max_absorbance,
        spectrum_model=spectrum_model,
        practical_rule_profile=artifacts["practical_profile"],
        conformers_attempted=len(records),
        conformers_successful=len(artifacts["successful"]),
        spectrum_csv=str(artifacts["spectrum_csv"]),
        spectrum_png=str(artifacts["spectrum_png"]),
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
        "method_keywords": (metadata.get("harmonic_method_profile") or {}).get("orca_keywords"),
    }
    resolved = {key: value for key, value in defaults.items() if value is not None}
    resolved.update({key: value for key, value in settings.items() if value is not None})
    return run_ir_spectrum(xyz_paths, run_dir, **resolved)
