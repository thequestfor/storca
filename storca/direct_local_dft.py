"""Direct projected-local-DFT O--H ensemble assembly.

This module is intentionally independent of :mod:`frequency_transfer`.  Local
frequencies and intensities enter only from retained ORCA finite differences;
the remainder of the spectrum comes from retained harmonic ORCA output.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from .experimental_ir import apply_measurement_response, resolve_experimental_profile
from .local_mode_fd import calculate_orca_local_modes
from .spectrum import format_spectrum, write_spectrum_csv


DIRECT_LOCAL_DFT_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class DirectLocalDFTConfig:
    schema_version: int = DIRECT_LOCAL_DFT_SCHEMA_VERSION
    oh_region_min_cm_1: float = 2800.0
    oh_region_max_cm_1: float = 3900.0
    grid_min_cm_1: float = 400.0
    grid_max_cm_1: float = 4000.0
    grid_step_cm_1: float = 0.5
    residual_fwhm_cm_1: float = 6.0
    expected_valid_local_modes: int | None = None
    harmonic_scale_factor: float = 0.97
    minimum_valid_environments_per_class: int = 3
    material_fingerprint_degradation: float = 0.02


def plan_six_representative_direct_local_dft(
    conformers: list[dict], *, maximum_additional_orca_invocations: int = 24,
) -> dict:
    """Plan representatives 4--6 without silently exceeding the stated budget."""
    jobs = []
    for position in range(4, min(6, len(conformers)) + 1):
        conformer = conformers[position - 1]
        bonds = [
            dict(item) for item in (conformer.get("local_stretch_bonds") or [])
            if str(item.get("bond_class", "")).upper().replace("–", "-").startswith("O-H")
        ]
        frequency_output = str(conformer.get("frequency_output") or "")
        retained_path = (
            Path(frequency_output).parent / "local-mode-finite-differences" / "local-modes.json"
            if frequency_output else Path("")
        )
        retained = json.loads(retained_path.read_text()) if frequency_output and retained_path.is_file() else {}
        retained_modes = retained.get("modes") or []
        cached = len(retained_modes) == len(bonds) and all(
            item.get("status") == "validated" for item in retained_modes
        )
        jobs.append({
            "representative_position": position,
            "xyz": str(conformer.get("optimized_xyz") or ""),
            "frequency_output": frequency_output,
            "bonds": bonds,
            "cached_valid_modes": len(retained_modes) if cached else 0,
            "orca_invocations": 0 if cached else 4 * len(bonds),
        })
    required = sum(item["orca_invocations"] for item in jobs)
    existing = collect_direct_local_modes(conformers[:3])
    existing_valid = sum(item["validation_status"] == "validated" for item in existing)
    requested_total_modes = len(existing) + sum(len(item["bonds"]) for item in jobs)
    return {
        "schema_version": DIRECT_LOCAL_DFT_SCHEMA_VERSION,
        "kind": "six_representative_direct_local_dft_plan",
        "ncores_per_orca_invocation": 8,
        "maximum_additional_orca_invocations": maximum_additional_orca_invocations,
        "required_additional_orca_invocations": required,
        "existing_valid_modes_representatives_1_to_3": existing_valid,
        "requested_total_modes_representatives_1_to_6": requested_total_modes,
        "acceptance_expected_total_modes": requested_total_modes,
        "acceptance_mode_count_basis": "derived_from_frozen_six_representative_manifest",
        "jobs": jobs,
        "status": "planned" if required <= maximum_additional_orca_invocations and len(jobs) == 3 else "failed_closed_budget_or_coverage",
        "budget_consistency": (
            "consistent" if required <= maximum_additional_orca_invocations else
            "specified_bonds_require_more_than_24_gradient_invocations"
        ),
        "mode_count_consistency": "manifest_derived",
    }


def run_six_representative_direct_local_dft(
    run_dir: Path, *, maximum_additional_orca_invocations: int = 24,
    ncores: int = 8, method_keywords: list[str] | None = None,
    nist_reference_csv: Path | None = None, baseline_spectrum_csv: Path | None = None,
    progress=None,
) -> dict:
    """Acquire representatives 4--6, then assemble all six; fail before ORCA on over-budget plans."""
    run_dir = Path(run_dir)
    conformers_path = run_dir / "clusters" / "conformers.json"
    if not conformers_path.is_file():
        return {"status": "failed_closed_missing_conformers"}
    conformers = json.loads(conformers_path.read_text())
    plan = plan_six_representative_direct_local_dft(
        conformers, maximum_additional_orca_invocations=maximum_additional_orca_invocations,
    )
    plan_path = run_dir / "direct-local-dft-plan.json"
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    if plan["status"] != "planned":
        return {**plan, "artifact": str(plan_path)}
    ledger_path = run_dir / "direct-local-dft-orca-invocations.json"
    results = []
    failed_closed = False
    for job in plan["jobs"]:
        position = int(job["representative_position"])
        xyz = Path(job["xyz"])
        frequency_output = Path(job["frequency_output"])
        if not xyz.is_file() or not frequency_output.parent.is_dir():
            results.append({**job, "status": "missing_representative_artifact"})
            break
        if progress:
            progress(f"Direct local DFT representative {position}/6 on {ncores} cores")
        try:
            root = frequency_output.parent / "local-mode-finite-differences"
            completed_modes = []
            bond_results = []
            for bond_index, bond in enumerate(job["bonds"]):
                result = calculate_orca_local_modes(
                    xyz, [bond], root / f"direct-bond-{bond_index:03d}",
                    hard_orca_invocation_cap=maximum_additional_orca_invocations,
                    ncores=ncores, method_keywords=method_keywords,
                    orca_ledger_path=ledger_path,
                    ledger_namespace=f"direct-representative-{position:03d}-bond-{bond_index:03d}",
                )
                bond_results.append({
                    "bond_index": bond_index, "status": result["status"],
                    "artifact": str(root / f"direct-bond-{bond_index:03d}" / "local-modes.json"),
                })
                completed_modes.extend(result.get("modes") or [])
                if result["status"] != "completed":
                    failed_closed = True
                    break
            combined = {
                "schema_version": DIRECT_LOCAL_DFT_SCHEMA_VERSION,
                "kind": "orca_projected_local_modes_direct_ensemble_input",
                "status": "completed" if len(completed_modes) == len(job["bonds"]) and all(item.get("status") == "validated" for item in completed_modes) else "failed_validation",
                "modes": completed_modes, "bond_results": bond_results,
                "orca_invocation_ledger": str(ledger_path),
            }
            root.mkdir(parents=True, exist_ok=True)
            combined_path = root / "local-modes.json"
            combined_path.write_text(json.dumps(combined, indent=2, sort_keys=True) + "\n")
            conformers[position - 1]["local_mode_output"] = str(combined_path)
            results.append({**job, "status": combined["status"], "bond_results": bond_results})
            if failed_closed:
                break
        except Exception as error:
            results.append({**job, "status": "failed", "error": str(error)})
            failed_closed = True
            break
    conformers_path.write_text(json.dumps(conformers, indent=2, sort_keys=True) + "\n")
    report = write_direct_local_dft_artifacts(
        run_dir, conformers[:6],
        config=DirectLocalDFTConfig(expected_valid_local_modes=plan["acceptance_expected_total_modes"]),
        nist_reference_csv=nist_reference_csv,
        baseline_spectrum_csv=baseline_spectrum_csv,
    )
    report["acquisition"] = {
        "plan": str(plan_path), "orca_invocation_ledger": str(ledger_path), "results": results,
    }
    Path(report["artifact"]).write_text(json.dumps(
        {key: value for key, value in report.items() if key != "artifact"},
        indent=2, sort_keys=True,
    ) + "\n")
    return report


def _finite(value, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def coordination_class(bond: dict) -> str:
    return ":".join((
        str(bond.get("bond_class", "O-H")),
        str(bond.get("spectral_band_class", "oh_stretch")),
        str(bond.get("coordination_class", "unclassified")),
    ))


def _geometry_provenance(path: Path) -> dict:
    path = Path(path)
    return {
        "path": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None,
        "status": "retained" if path.is_file() else "missing",
    }


def collect_direct_local_modes(conformers: Iterable[dict]) -> list[dict]:
    """Collect an auditable mode row for every requested local oscillator."""
    rows: list[dict] = []
    for position, conformer in enumerate(conformers, start=1):
        geometry = Path(str(conformer.get("optimized_xyz") or ""))
        local_path_value = conformer.get("local_mode_output")
        if not local_path_value and conformer.get("frequency_output"):
            local_path_value = str(
                Path(conformer["frequency_output"]).parent
                / "local-mode-finite-differences" / "local-modes.json"
            )
        local_path = Path(str(local_path_value or ""))
        payload = json.loads(local_path.read_text()) if local_path.is_file() else {}
        modes = payload.get("modes") or []
        expected_bonds = conformer.get("local_stretch_bonds") or []
        count = max(len(expected_bonds), len(modes))
        raw_weight = _finite(
            conformer.get("population_weight", conformer.get("weight", 0.0)), 0.0,
        )
        molecule_count = max(1, int(conformer.get("cluster_size") or 1))
        for bond_index in range(count):
            mode = modes[bond_index] if bond_index < len(modes) else {}
            bond = dict(
                mode.get("bond")
                or (expected_bonds[bond_index] if bond_index < len(expected_bonds) else {})
            )
            status = str(mode.get("status") or "missing")
            source = "orca_projected_local_mode_finite_difference" if mode else "missing"
            rows.append({
                "representative_position": position,
                "independent_environment_id": str(
                    conformer.get("independent_environment_id") or f"representative-{position}"
                ),
                "geometry": _geometry_provenance(geometry),
                "local_mode_artifact": str(local_path),
                "bond": bond,
                "coordination_class": coordination_class(bond),
                "trajectory_occupancy_weight": raw_weight,
                "molecule_normalization_divisor": molecule_count,
                "frequency_cm-1": mode.get("frequency_cm-1"),
                "intensity_km_mol": mode.get("intensity_km_mol"),
                "frequency_step_disagreement_cm-1": mode.get(
                    "frequency_step_disagreement_cm-1"
                ),
                "validation_status": status,
                "frequency_source": source,
                "intensity_source": source,
            })
    return rows


def coordination_class_coverage(mode_rows: Iterable[dict]) -> dict[str, dict]:
    grouped: dict[str, list[dict]] = {}
    for row in mode_rows:
        if "coordination_class" not in row:
            continue
        grouped.setdefault(str(row["coordination_class"]), []).append(row)
    report = {}
    for key, members in sorted(grouped.items()):
        valid = [item for item in members if item["validation_status"] == "validated"]
        environment_weights: dict[str, float] = {}
        all_environment_weights: dict[str, float] = {}
        for item in members:
            env = str(item["independent_environment_id"])
            weight = max(0.0, _finite(item["trajectory_occupancy_weight"]))
            all_environment_weights[env] = max(all_environment_weights.get(env, 0.0), weight)
        for item in valid:
            env = str(item["independent_environment_id"])
            weight = max(0.0, _finite(item["trajectory_occupancy_weight"]))
            environment_weights[env] = max(environment_weights.get(env, 0.0), weight)
        total = sum(all_environment_weights.values())
        covered = sum(environment_weights.values())
        normalized = [value / covered for value in environment_weights.values()] if covered else []
        ess = 1.0 / sum(value * value for value in normalized) if normalized else 0.0
        report[key] = {
            "requested_modes": len(members),
            "valid_modes": len(valid),
            "valid_independent_environments": len(environment_weights),
            "trajectory_occupancy_total": total,
            "trajectory_occupancy_covered": covered,
            "occupancy_coverage_fraction": covered / total if total else 0.0,
            "effective_sample_size": ess,
        }
    return report


def _gaussian(grid: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    return np.exp(-0.5 * np.square((grid - center) / sigma)) / (
        sigma * math.sqrt(2.0 * math.pi)
    )


def assemble_direct_local_dft(
    conformers: list[dict], *, config: DirectLocalDFTConfig | None = None,
) -> tuple[np.ndarray, np.ndarray, list[dict], dict]:
    """Replace only harmonic O--H-region modes with validated direct modes."""
    resolved = config or DirectLocalDFTConfig()
    grid = np.arange(
        resolved.grid_min_cm_1,
        resolved.grid_max_cm_1 + resolved.grid_step_cm_1 / 2,
        resolved.grid_step_cm_1,
    )
    signal = np.zeros_like(grid)
    rows = collect_direct_local_modes(conformers)
    weights = np.asarray([
        max(0.0, _finite(item.get("population_weight", item.get("weight", 0.0))))
        for item in conformers
    ])
    if weights.sum() <= 0 and len(weights):
        weights[:] = 1.0
    if weights.sum() > 0:
        weights /= weights.sum()

    harmonic_rows = []
    for position, conformer in enumerate(conformers, start=1):
        divisor = max(1, int(conformer.get("cluster_size") or 1))
        harmonic_modes = conformer.get("ir_modes") or []
        if not harmonic_modes and conformer.get("frequency_output"):
            try:
                from src.parser import parse_orca_ir
                harmonic_modes = parse_orca_ir(Path(conformer["frequency_output"]))
            except (OSError, ValueError):
                harmonic_modes = []
        for mode in harmonic_modes:
            raw_frequency = _finite(mode.get("freq"), -1.0)
            intensity = max(0.0, _finite(mode.get("intensity")))
            if raw_frequency <= 0 or resolved.oh_region_min_cm_1 <= raw_frequency <= resolved.oh_region_max_cm_1:
                continue
            frequency = raw_frequency * resolved.harmonic_scale_factor
            strength = (weights[position - 1] if len(weights) else 0.0) * intensity / divisor
            signal += strength * _gaussian(grid, frequency, resolved.residual_fwhm_cm_1)
            harmonic_rows.append({
                "representative_position": position,
                "geometry": _geometry_provenance(Path(str(conformer.get("optimized_xyz") or ""))),
                "mode": mode.get("mode"),
                "raw_frequency_cm-1": raw_frequency,
                "frequency_cm-1": frequency,
                "intensity_km_mol": intensity,
                "ensemble_weight": float(weights[position - 1]) if len(weights) else 0.0,
                "molecule_normalization_divisor": divisor,
                "frequency_source": "orca_harmonic",
                "intensity_source": "orca_harmonic",
                "validation_status": "preserved_outside_replacement_region",
                "region": "fingerprint_or_non_oh",
            })

    valid_rows = [item for item in rows if item["validation_status"] == "validated"]
    position_weights = {
        position: float(weights[position - 1])
        for position in range(1, len(weights) + 1)
    }
    for row in valid_rows:
        # Occupancy belongs to an environment, not to each oscillator.  Each
        # oscillator receives that environment weight and is then normalized
        # per methanol molecule below.
        ensemble_weight = position_weights.get(int(row["representative_position"]), 0.0)
        row["ensemble_weight"] = ensemble_weight
        strength = (
            ensemble_weight * max(0.0, _finite(row["intensity_km_mol"]))
            / int(row["molecule_normalization_divisor"])
        )
        row["raw_frequency_cm-1"] = _finite(row["frequency_cm-1"])
        row["frequency_cm-1"] = row["raw_frequency_cm-1"] * resolved.harmonic_scale_factor
        signal += strength * _gaussian(
            grid, _finite(row["frequency_cm-1"]), resolved.residual_fwhm_cm_1,
        )
    provenance = {
        "oh_frequency_source": "orca_projected_local_mode_finite_difference",
        "oh_intensity_source": "orca_projected_local_mode_finite_difference",
        "non_oh_source": "orca_harmonic",
        "replacement_region_cm-1": [
            resolved.oh_region_min_cm_1, resolved.oh_region_max_cm_1,
        ],
        "trajectory_weighting": "retained_trajectory_occupancy_normalized_over_valid_modes",
        "cluster_normalization": "each_contribution_divided_by_methanol_molecule_count",
        "harmonic_scale_factor": resolved.harmonic_scale_factor,
        "scale_factor_scope": "same_method_level_factor_applied_to_orca_harmonic_and_projected_local_modes",
        "xtb_frequency_correction_used_in_oh_region": False,
        "harmonic_modes_preserved": len(harmonic_rows),
    }
    return grid, signal, rows + harmonic_rows, provenance


def spectrum_band_metrics(grid: np.ndarray, signal: np.ndarray, low: float, high: float) -> dict:
    mask = (grid >= low) & (grid <= high)
    x, y = grid[mask], signal[mask]
    area = float(np.trapezoid(y, x)) if len(x) else 0.0
    center = float(np.trapezoid(x * y, x) / area) if area > 0 else None
    if not len(y) or float(y.max()) <= 0:
        return {"center_cm-1": center, "fwhm_cm-1": None, "integrated_intensity": area}
    above = x[y >= 0.5 * float(y.max())]
    width = float(above[-1] - above[0]) if len(above) else None
    return {"center_cm-1": center, "fwhm_cm-1": width, "integrated_intensity": area}


def _read_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "wavenumber_cm-1" not in rows[0]:
        raise ValueError(f"Invalid spectrum CSV: {path}")
    value_key = next(key for key in rows[0] if key != "wavenumber_cm-1")
    grid = np.asarray([float(row["wavenumber_cm-1"]) for row in rows])
    values = np.asarray([float(row[value_key]) for row in rows])
    order = np.argsort(grid)
    return grid[order], values[order]


def _normalized_correlation(
    grid: np.ndarray, signal: np.ndarray, reference_path: Path, low: float, high: float,
) -> float | None:
    ref_x, ref_y = _read_spectrum(reference_path)
    mask = (grid >= low) & (grid <= high)
    if mask.sum() < 2:
        return None
    interpolated = np.interp(grid[mask], ref_x, ref_y)
    left, right = signal[mask], interpolated
    if np.std(left) <= 0 or np.std(right) <= 0:
        return None
    return float(np.corrcoef(left, right)[0, 1])


def _absorbance_from_transmittance_percent(values: np.ndarray) -> np.ndarray:
    transmission = np.asarray(values, dtype=float) / 100.0
    if np.any(~np.isfinite(transmission)) or np.any(transmission <= 0):
        raise ValueError("Transmittance-to-absorbance conversion requires positive finite percentages")
    return -np.log10(transmission)


def _regional_shape_metrics(
    grid: np.ndarray, predicted: np.ndarray, reference: np.ndarray, low: float, high: float,
) -> dict:
    mask = (grid >= low) & (grid <= high)
    left, right = predicted[mask], reference[mask]
    if len(left) < 2:
        return {"correlation": None, "normalized_rmse": None}
    left_scale, right_scale = max(float(left.max()), 1.0e-300), max(float(right.max()), 1.0e-300)
    left, right = left / left_scale, right / right_scale
    correlation = float(np.corrcoef(left, right)[0, 1]) if np.std(left) > 0 and np.std(right) > 0 else None
    return {"correlation": correlation, "normalized_rmse": float(np.sqrt(np.mean(np.square(left - right))))}


def write_direct_local_dft_artifacts(
    run_dir: Path, conformers: list[dict], *,
    config: DirectLocalDFTConfig | None = None,
    nist_reference_csv: Path | None = None,
    baseline_spectrum_csv: Path | None = None,
) -> dict:
    """Assemble spectra, metrics, and a fail-closed promotion decision."""
    resolved = config or DirectLocalDFTConfig()
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    grid, intrinsic, ledger, provenance = assemble_direct_local_dft(
        conformers, config=resolved,
    )
    intrinsic_path = write_spectrum_csv(
        run_dir / "spectrum_direct_local_dft_intrinsic.csv", grid, intrinsic,
        "absorption_strength",
    )
    metadata_path = run_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    profile = resolve_experimental_profile(
        phase=str(metadata.get("phase", "liquid")), measurement="transmission",
        instrument_resolution_cm_1=_finite(metadata.get("instrument_resolution_cm_1"), 4.0),
        apodization=str(metadata.get("apodization", "happ-genzel")),
        residual_fwhm_cm_1=_finite(metadata.get("residual_fwhm_cm_1"), 6.0),
    )
    measured = apply_measurement_response(grid, intrinsic, profile=profile)
    transmission, column = format_spectrum(
        measured, style="transmittance", max_absorbance=_finite(metadata.get("max_absorbance"), 1.0),
    )
    transmission_path = write_spectrum_csv(
        run_dir / "spectrum_direct_local_dft_transmission.csv", grid, transmission, column,
    )
    coverage = coordination_class_coverage(
        item for item in ledger if item.get("frequency_source") != "orca_harmonic"
    )
    valid = [
        item for item in ledger
        if item.get("frequency_source") != "orca_harmonic"
        and item.get("validation_status") == "validated"
    ]
    expected_modes = (
        int(resolved.expected_valid_local_modes)
        if resolved.expected_valid_local_modes is not None else
        sum(len(item.get("local_stretch_bonds") or []) for item in conformers)
    )
    gates = {
        "all_expected_local_modes_step_validated": {
            "passed": len(valid) == expected_modes,
            "valid": len(valid), "expected": expected_modes,
        },
        "no_xtb_frequency_correction_in_oh_region": {
            "passed": not any("xtb" in str(item.get("frequency_source", "")).lower() for item in valid),
        },
        "minimum_three_valid_environments_per_coordination_class": {
            "passed": bool(coverage) and all(
                item["valid_independent_environments"] >= resolved.minimum_valid_environments_per_class
                for item in coverage.values()
            ),
            "minimum": resolved.minimum_valid_environments_per_class,
        },
    }
    metrics = {
        "direct_local_dft": {
            "oh_region": spectrum_band_metrics(
                grid, intrinsic, resolved.oh_region_min_cm_1, resolved.oh_region_max_cm_1,
            ),
            "fingerprint_region": spectrum_band_metrics(grid, intrinsic, 400.0, 1800.0),
        }
    }
    comparison_rows = {
        "wavenumber_cm-1": grid, "direct_local_dft_intrinsic": intrinsic,
        "direct_local_dft_transmittance_percent": transmission,
    }
    if baseline_spectrum_csv and Path(baseline_spectrum_csv).is_file():
        base_x, base_y = _read_spectrum(Path(baseline_spectrum_csv))
        baseline = np.interp(grid, base_x, base_y)
        comparison_rows["baseline_input"] = baseline
        metrics["baseline"] = {
            "oh_region": spectrum_band_metrics(grid, baseline, resolved.oh_region_min_cm_1, resolved.oh_region_max_cm_1),
            "fingerprint_region": spectrum_band_metrics(grid, baseline, 400.0, 1800.0),
        }
    if nist_reference_csv and Path(nist_reference_csv).is_file():
        ref_x, ref_transmission = _read_spectrum(Path(nist_reference_csv))
        reference_transmission = np.interp(grid, ref_x, ref_transmission)
        reference_absorbance = _absorbance_from_transmittance_percent(reference_transmission)
        direct_absorbance = _absorbance_from_transmittance_percent(transmission)
        comparison_rows["nist_transmittance_percent"] = reference_transmission
        comparison_rows["nist_absorbance"] = reference_absorbance
        comparison_rows["direct_local_dft_absorbance"] = direct_absorbance
        nist_oh = spectrum_band_metrics(grid, reference_absorbance, resolved.oh_region_min_cm_1, resolved.oh_region_max_cm_1)
        direct_oh = spectrum_band_metrics(grid, direct_absorbance, resolved.oh_region_min_cm_1, resolved.oh_region_max_cm_1)
        direct_center = direct_oh["center_cm-1"]
        direct_error = abs(direct_center - nist_oh["center_cm-1"]) if direct_center and nist_oh["center_cm-1"] else None
        baseline_error = None
        baseline_shape = None
        if "baseline" in metrics:
            baseline_absorbance = _absorbance_from_transmittance_percent(baseline)
            comparison_rows["baseline_absorbance"] = baseline_absorbance
            baseline_center = spectrum_band_metrics(grid, baseline_absorbance, resolved.oh_region_min_cm_1, resolved.oh_region_max_cm_1)["center_cm-1"]
            baseline_error = abs(baseline_center - nist_oh["center_cm-1"]) if baseline_center and nist_oh["center_cm-1"] else None
            baseline_shape = _regional_shape_metrics(grid, baseline_absorbance, reference_absorbance, 400.0, 1800.0)
        direct_shape = _regional_shape_metrics(grid, direct_absorbance, reference_absorbance, 400.0, 1800.0)
        metrics["nist_comparison"] = {
            "reference_oh": nist_oh,
            "direct_oh": direct_oh,
            "direct_oh_center_error_cm-1": direct_error,
            "baseline_oh_center_error_cm-1": baseline_error,
            "direct_fingerprint": direct_shape,
            "baseline_fingerprint": baseline_shape,
            "comparison_basis": "absorbance_derived_from_condition_matched_percent_transmittance",
        }
        gates["nist_oh_center_improves"] = {
            "passed": direct_error is not None and baseline_error is not None and direct_error < baseline_error,
        }
        gates["fingerprint_not_materially_degraded"] = {
            "passed": direct_shape["correlation"] is not None and baseline_shape is not None
            and baseline_shape["correlation"] is not None
            and direct_shape["correlation"] >= baseline_shape["correlation"] - resolved.material_fingerprint_degradation
            and direct_shape["normalized_rmse"] <= baseline_shape["normalized_rmse"] * 1.05,
            "maximum_allowed_correlation_loss": resolved.material_fingerprint_degradation,
            "maximum_allowed_normalized_rmse_increase_fraction": 0.05,
        }
    else:
        gates["nist_oh_center_improves"] = {"passed": False, "status": "not_evaluated_missing_reference"}
        gates["fingerprint_not_materially_degraded"] = {"passed": False, "status": "not_evaluated_missing_reference"}

    comparison_path = run_dir / "spectrum_direct_local_dft_comparison.csv"
    with comparison_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        keys = list(comparison_rows)
        writer.writerow(keys)
        writer.writerows(zip(*(comparison_rows[key] for key in keys)))
    metrics_path = run_dir / "direct-local-dft-nist-regional-metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    passed = bool(gates) and all(item["passed"] for item in gates.values())
    report = {
        "schema_version": DIRECT_LOCAL_DFT_SCHEMA_VERSION,
        "kind": "direct_local_dft_ensemble",
        "status": "accepted" if passed else "failed_closed",
        "promotion_to_adaptive_acquisition": passed,
        "configuration": asdict(resolved),
        "provenance": provenance,
        "coordination_classes": coverage,
        "modes": ledger,
        "acceptance_gates": gates,
        "artifacts": {
            "intrinsic_spectrum": str(intrinsic_path),
            "transmission_spectrum": str(transmission_path),
            "comparison_spectrum": str(comparison_path),
            "nist_regional_metrics": str(metrics_path),
        },
    }
    artifact = run_dir / "direct-local-dft-ensemble.json"
    artifact.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(artifact)
    return report
