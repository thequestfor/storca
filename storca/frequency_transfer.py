"""Validated additive DFT corrections for inexpensive environment frequencies."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from .ir_modes import local_stretch_mode_assignments, parse_orca_normal_modes
from .mode_character import fingerprint_similarity, mode_character_fingerprints


TRANSFER_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class TransferValidationConfig:
    minimum_representative_pairs_per_class: int = 3
    minimum_mode_character_similarity: float = 0.50
    minimum_raw_error_improvement_fraction: float = 0.20
    maximum_high_frequency_loo_mae_cm_1: float = 30.0
    maximum_other_frequency_loo_mae_cm_1: float = 20.0
    maximum_withheld_error_cm_1: float = 75.0
    maximum_applicability_distance: float = 3.0
    minimum_ensemble_covered_fraction: float = 0.80


def _finite(value, default: float = 0.0) -> float:
    return float(value) if isinstance(value, (int, float)) and math.isfinite(value) else default


def _environment_distance(left: dict, right: dict) -> float:
    """Dimensionless distance using features available for every cluster size."""
    scales = {
        "h_bond_distance_angstrom": 0.20,
        "donor_h_acceptor_angle_degrees": 20.0,
        "donor_acceptor_distance_angstrom": 0.25,
        "heavy_atom_rmsd_angstrom": 0.35,
        "relative_xtb_energy_kcal_mol": 4.0,
    }
    terms = []
    for key, scale in scales.items():
        lvalue, rvalue = left.get(key), right.get(key)
        if all(isinstance(value, (int, float)) and math.isfinite(value) for value in (lvalue, rvalue)):
            terms.append(((float(lvalue) - float(rvalue)) / scale) ** 2)
    l_angle, r_angle = left.get("intermolecular_orientation_degrees"), right.get("intermolecular_orientation_degrees")
    if all(isinstance(value, (int, float)) and math.isfinite(value) for value in (l_angle, r_angle)):
        gap = abs(float(l_angle) - float(r_angle)) % 360.0
        terms.append((min(gap, 360.0 - gap) / 60.0) ** 2)
    return math.sqrt(sum(terms) / len(terms)) if terms else math.inf


def _class_key(item: dict) -> str:
    return ":".join((
        str(item.get("bond_class", "X-H")),
        str(item.get("spectral_band_class", "xh_stretch")),
        str(item.get("coordination_class", "unclassified")),
    ))


def pair_representative_local_modes(
    representative: dict, xtb_record: dict, orca_hessian: Path,
    orca_ir_modes: list[dict] | None = None,
) -> list[dict]:
    """Pair xTB and ORCA X-H modes at an identical representative geometry."""
    if xtb_record.get("frequency_status") != "completed":
        return []
    orca = parse_orca_normal_modes(orca_hessian)
    bonds = representative.get("local_stretch_bonds") or []
    orca_assignments = local_stretch_mode_assignments(orca, bonds)
    xtb_assignments = xtb_record.get("local_stretch_modes") or []
    orca_fingerprints = {
        int(item["mode"]): item for item in mode_character_fingerprints(orca)
    }
    xtb_fingerprints = {
        int(item["mode"]): item for item in xtb_record.get("mode_character_fingerprints") or []
    }
    xtb_modes = {int(item["mode"]): item for item in xtb_record.get("modes") or []}
    orca_intensities = {
        int(item["mode"]): float(item["intensity"]) for item in (orca_ir_modes or [])
        if isinstance(item.get("intensity"), (int, float))
    }
    identity = lambda item: (
        int(item.get("molecule_index", -1)), int(item["heavy_atom"]), int(item["hydrogen_atom"])
    )
    xtb_by_identity = {identity(item): item for item in xtb_assignments}
    pairs = []
    for orca_item in orca_assignments:
        xtb_item = xtb_by_identity.get(identity(orca_item))
        if not xtb_item:
            continue
        orca_mode, xtb_mode = int(orca_item["mode"]), int(xtb_item["mode"])
        if orca_mode not in orca_fingerprints or xtb_mode not in xtb_fingerprints:
            continue
        similarity = fingerprint_similarity(
            orca_fingerprints[orca_mode], xtb_fingerprints[xtb_mode],
        )
        xtb_frequency = float(xtb_modes[xtb_mode]["freq"])
        dft_frequency = float(orca.frequencies_cm_1[orca_mode])
        if not all(math.isfinite(value) and value > 0 for value in (xtb_frequency, dft_frequency)):
            continue
        pairs.append({
            "representative_id": str(representative.get("source_xtb_candidate_id")),
            "environment_id": representative.get("independent_environment_id"),
            "mode_class": _class_key(orca_item),
            "bond": {key: orca_item.get(key) for key in (
                "bond_class", "molecule_index", "heavy_atom", "hydrogen_atom",
                "spectral_band_class", "coordination_class",
            )},
            "xtb_mode": xtb_mode,
            "orca_mode": orca_mode,
            "xtb_frequency_cm-1": xtb_frequency,
            "dft_frequency_cm-1": dft_frequency,
            "dft_intensity_km_mol": orca_intensities.get(orca_mode),
            "additive_correction_cm-1": dft_frequency - xtb_frequency,
            "mode_character_similarity": similarity,
            "xtb_local_projection": float(xtb_item["local_stretch_projection"]),
            "orca_local_projection": float(orca_item["local_stretch_projection"]),
            "environment_features": representative.get("environment_features") or {},
        })
    return pairs


def apply_finite_difference_pair_overrides(
    pairs: Iterable[dict], payload: dict, *, exclude_unvalidated_snapshot_pairs: bool = False,
) -> list[dict]:
    """Prefer validated bond-local DFT finite differences over snapshot modes."""
    local_modes = payload.get("modes") or []
    validated = {}
    for mode in local_modes:
        bond = mode.get("bond") or {}
        if mode.get("status") != "validated":
            continue
        identity = (
            int(bond.get("molecule_index", -1)),
            int(mode.get("heavy_atom", bond.get("heavy_atom", -1))),
            int(mode.get("hydrogen_atom", bond.get("hydrogen_atom", -1))),
        )
        validated[identity] = mode
    output = []
    for pair in pairs:
        bond = pair.get("bond") or {}
        identity = (
            int(bond.get("molecule_index", -1)),
            int(bond.get("heavy_atom", -1)),
            int(bond.get("hydrogen_atom", -1)),
        )
        local = validated.get(identity)
        if local is None:
            if not exclude_unvalidated_snapshot_pairs:
                output.append(dict(pair))
            continue
        dft_frequency = float(local["frequency_cm-1"])
        updated = {
            **pair,
            "dft_frequency_cm-1": dft_frequency,
            "dft_intensity_km_mol": float(local["intensity_km_mol"]),
            "additive_correction_cm-1": dft_frequency - float(pair["xtb_frequency_cm-1"]),
            "dft_frequency_source": "orca_projected_local_mode_finite_difference",
            "dft_intensity_source": "orca_local_dipole_finite_difference",
            "finite_difference_validation": {
                "status": local.get("status"),
                "displacement_stability_status": local.get("displacement_stability_status"),
                "frequency_step_disagreement_cm-1": local.get(
                    "frequency_step_disagreement_cm-1"
                ),
            },
        }
        output.append(updated)
    return output


def validate_frequency_transfer(
    pairs: Iterable[dict], *, config: TransferValidationConfig | None = None,
) -> dict:
    """Leave one representative out and fail closed for every mode class."""
    resolved = config or TransferValidationConfig()
    grouped: dict[str, list[dict]] = {}
    for pair in pairs:
        grouped.setdefault(str(pair["mode_class"]), []).append(pair)
    classes = {}
    for mode_class, observations in sorted(grouped.items()):
        distinct = len({item["representative_id"] for item in observations})
        predictions = []
        for held in observations:
            training = [
                item for item in observations
                if item["representative_id"] != held["representative_id"]
            ]
            if not training:
                continue
            nearest = min(
                training,
                key=lambda item: _environment_distance(
                    held["environment_features"], item["environment_features"],
                ),
            )
            estimated = held["xtb_frequency_cm-1"] + nearest["additive_correction_cm-1"]
            predictions.append({
                "held_out_representative_id": held["representative_id"],
                "training_representative_id": nearest["representative_id"],
                "raw_error_cm-1": abs(held["dft_frequency_cm-1"] - held["xtb_frequency_cm-1"]),
                "corrected_error_cm-1": abs(held["dft_frequency_cm-1"] - estimated),
                "applicability_distance": _environment_distance(
                    held["environment_features"], nearest["environment_features"],
                ),
            })
        raw_mae = float(np.mean([item["raw_error_cm-1"] for item in predictions])) if predictions else math.inf
        corrected_mae = float(np.mean([item["corrected_error_cm-1"] for item in predictions])) if predictions else math.inf
        improvement = (raw_mae - corrected_mae) / raw_mae if math.isfinite(raw_mae) and raw_mae > 0 else -math.inf
        mean_frequency = float(np.mean([item["dft_frequency_cm-1"] for item in observations]))
        maximum_mae = (
            resolved.maximum_high_frequency_loo_mae_cm_1 if mean_frequency >= 2500.0
            else resolved.maximum_other_frequency_loo_mae_cm_1
        )
        failures = []
        if distinct < resolved.minimum_representative_pairs_per_class:
            failures.append("fewer_than_three_independent_representatives")
        if min(item["mode_character_similarity"] for item in observations) < resolved.minimum_mode_character_similarity:
            failures.append("low_mode_character_similarity")
        if improvement < resolved.minimum_raw_error_improvement_fraction:
            failures.append("leave_one_out_does_not_improve_raw_xtb")
        if corrected_mae > maximum_mae:
            failures.append("leave_one_out_mae_too_large")
        if predictions and max(item["corrected_error_cm-1"] for item in predictions) > resolved.maximum_withheld_error_cm_1:
            failures.append("leave_one_out_maximum_error_too_large")
        classes[mode_class] = {
            "status": "validated" if not failures else "insufficient_dft_transfer_validation",
            "failures": failures,
            "representative_pairs": len(observations),
            "independent_representatives": distinct,
            "raw_xtb_loo_mae_cm-1": raw_mae if math.isfinite(raw_mae) else None,
            "corrected_loo_mae_cm-1": corrected_mae if math.isfinite(corrected_mae) else None,
            "raw_error_improvement_fraction": improvement if math.isfinite(improvement) else None,
            "maximum_corrected_loo_error_cm-1": (
                max(item["corrected_error_cm-1"] for item in predictions) if predictions else None
            ),
            "correction_stddev_cm-1": float(np.std([
                item["additive_correction_cm-1"] for item in observations
            ])),
            "leave_one_representative_out": predictions,
        }
    return {
        "schema_version": TRANSFER_SCHEMA_VERSION,
        "kind": "leave_one_representative_out_dft_transfer_validation",
        "thresholds": asdict(resolved),
        "status": (
            "validated" if classes and all(item["status"] == "validated" for item in classes.values())
            else "insufficient_dft_transfer_validation"
        ),
        "classes": classes,
    }


def transfer_snapshot_modes(
    xtb_records: Iterable[dict], pairs: Iterable[dict], validation: dict,
    *, config: TransferValidationConfig | None = None,
) -> list[dict]:
    """Apply the nearest validated additive correction inside its coverage."""
    resolved = config or TransferValidationConfig()
    pair_list = list(pairs)
    output = []
    for record in xtb_records:
        if record.get("frequency_status") != "completed":
            continue
        features = record.get("environment_features") or {}
        transferred_modes = []
        for mode in record.get("modes") or []:
            if not mode.get("bond_class"):
                continue
            mode_class = _class_key(mode)
            class_validation = validation.get("classes", {}).get(mode_class, {})
            training = [item for item in pair_list if item["mode_class"] == mode_class]
            if class_validation.get("status") != "validated" or not training:
                continue
            nearest = min(
                training,
                key=lambda item: _environment_distance(features, item["environment_features"]),
            )
            distance = _environment_distance(features, nearest["environment_features"])
            if distance > resolved.maximum_applicability_distance:
                continue
            correction_spread = _finite(class_validation.get("correction_stddev_cm-1"))
            loo_mae = _finite(class_validation.get("corrected_loo_mae_cm-1"))
            uncertainty = math.sqrt(loo_mae ** 2 + correction_spread ** 2 + (10.0 * distance) ** 2)
            transferred_modes.append({
                **mode,
                "raw_xtb_frequency_cm-1": float(mode["freq"]),
                "freq": float(mode["freq"]) + nearest["additive_correction_cm-1"],
                "dft_correction_cm-1": nearest["additive_correction_cm-1"],
                "transfer_representative_id": nearest["representative_id"],
                "mode_class": mode_class,
                "applicability_distance": distance,
                "estimated_uncertainty_cm-1": uncertainty,
                "intensity_model": "raw_xtb_diagnostic_not_dft_transferred",
            })
        output.append({
            "candidate_id": record["candidate_id"],
            "environment_features": features,
            "population_weight": record.get("population_weight"),
            "cluster_size": int(record.get("cluster_size") or 1),
            "topology": record.get("topology", "dimer"),
            "acquisition_round": int(record.get("acquisition_round", 0)),
            "modes": transferred_modes,
        })
    return output


def apply_ensemble_transfer_coverage(
    validation: dict, xtb_records: Iterable[dict], pairs: Iterable[dict],
    *, config: TransferValidationConfig | None = None,
) -> dict:
    """Require validated representatives to cover the sampled mode classes."""
    resolved = config or TransferValidationConfig()
    pair_list = list(pairs)
    totals: dict[str, int] = {}
    covered: dict[str, int] = {}
    for record in xtb_records:
        if record.get("frequency_status") != "completed":
            continue
        features = record.get("environment_features") or {}
        for mode in record.get("modes") or []:
            if not mode.get("bond_class"):
                continue
            mode_class = _class_key(mode)
            totals[mode_class] = totals.get(mode_class, 0) + 1
            training = [item for item in pair_list if item["mode_class"] == mode_class]
            class_result = validation.setdefault("classes", {}).setdefault(mode_class, {
                "status": "insufficient_dft_transfer_validation",
                "failures": ["no_representative_pairs_for_sampled_class"],
                "representative_pairs": 0,
                "independent_representatives": 0,
            })
            if class_result.get("status") != "validated" or not training:
                continue
            distance = min(
                _environment_distance(features, item["environment_features"])
                for item in training
            )
            if distance <= resolved.maximum_applicability_distance:
                covered[mode_class] = covered.get(mode_class, 0) + 1
    for mode_class, total in totals.items():
        class_result = validation["classes"][mode_class]
        fraction = covered.get(mode_class, 0) / total
        class_result["sampled_observations"] = total
        class_result["covered_sampled_observations"] = covered.get(mode_class, 0)
        class_result["ensemble_covered_fraction"] = fraction
        if fraction < resolved.minimum_ensemble_covered_fraction:
            failures = class_result.setdefault("failures", [])
            if "insufficient_sampled_ensemble_coverage" not in failures:
                failures.append("insufficient_sampled_ensemble_coverage")
            class_result["status"] = "insufficient_dft_transfer_validation"
    validation["status"] = (
        "validated" if validation.get("classes")
        and all(item.get("status") == "validated" for item in validation["classes"].values())
        else "insufficient_dft_transfer_validation"
    )
    return validation


def _write_spectrum(path: Path, records: Iterable[dict], *, raw: bool) -> Path:
    grid = np.arange(400.0, 4000.0 + 0.5, 0.5)
    intensity = np.zeros_like(grid)
    records = list(records)
    weight = 1.0 / len(records) if records else 0.0
    sigma = 6.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    for record in records:
        for mode in record.get("modes") or []:
            frequency = float(mode.get("raw_xtb_frequency_cm-1") if raw else mode["freq"])
            strength = max(0.0, _finite(mode.get("intensity"))) / max(
                1, int(record.get("cluster_size") or 1)
            )
            intensity += weight * strength * np.exp(-0.5 * np.square((grid - frequency) / sigma))
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber_cm-1", "absorption_strength"])
        writer.writerows(zip(grid, intensity))
    return path


def _hybrid_multifidelity_intrinsic(
    run_dir: Path, transferred: list[dict], pairs: list[dict],
    selected: list[dict], conformers: list[dict], *, scale_factor: float,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Combine nonlocal ORCA modes with transferred frequencies and DFT class intensities."""
    from src.parser import parse_orca_ir

    grid = np.arange(400.0, 4000.0 + 0.5, 0.5)
    signal = np.zeros_like(grid)
    sigma = 6.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    normalization = sigma * math.sqrt(2.0 * math.pi)
    local_by_representative: dict[str, set[int]] = {}
    class_intensities: dict[str, list[float]] = {}
    for pair in pairs:
        local_by_representative.setdefault(str(pair["representative_id"]), set()).add(int(pair["orca_mode"]))
        intensity = pair.get("dft_intensity_km_mol")
        if isinstance(intensity, (int, float)) and math.isfinite(intensity) and intensity >= 0:
            class_intensities.setdefault(str(pair["mode_class"]), []).append(float(intensity))
    dft_class_intensity = {
        key: float(np.mean(values)) for key, values in class_intensities.items()
    }
    representative_weight = 1.0 / len(conformers) if conformers else 0.0
    for representative, conformer in zip(selected, conformers):
        output = conformer.get("frequency_output")
        if not output:
            continue
        local_modes = local_by_representative.get(str(representative.get("source_xtb_candidate_id")), set())
        divisor = max(1, int(conformer.get("cluster_size") or representative.get("cluster_size") or 1))
        # Representative environments are deliberately stratified, not a
        # vacuum-energy population.  Give each selected stratum equal weight.
        weight = representative_weight
        for mode in parse_orca_ir(output):
            if int(mode["mode"]) in local_modes or float(mode["freq"]) <= 0:
                continue
            frequency = float(mode["freq"]) * scale_factor
            strength = weight * max(0.0, float(mode["intensity"])) / divisor
            signal += strength * np.exp(-0.5 * np.square((grid - frequency) / sigma)) / normalization
    candidate_weight = 1.0 / len(transferred) if transferred else 0.0
    local_modes_added = 0
    local_modes_expected = 0
    for record in transferred:
        divisor = max(1, int(record.get("cluster_size") or 1))
        for mode in record.get("modes") or []:
            local_modes_expected += 1
            intensity = dft_class_intensity.get(str(mode["mode_class"]))
            if intensity is None:
                continue
            frequency = float(mode["freq"]) * scale_factor
            strength = candidate_weight * intensity / divisor
            signal += strength * np.exp(-0.5 * np.square((grid - frequency) / sigma)) / normalization
            local_modes_added += 1
    return grid, signal, {
        "frequency_source": "validated_dft_corrected_xtb_snapshot_frequencies",
        "local_intensity_source": "representative_orca_mode_class_mean",
        "nonlocal_mode_source": "representative_orca_harmonic_modes",
        "population_source": "stratified_equal_candidate_weights",
        "cluster_normalization": "intensity_divided_by_target_molecule_count",
        "scale_factor": scale_factor,
        "local_modes_added": local_modes_added,
        "local_modes_expected": local_modes_expected,
        "local_intensity_coverage_fraction": (
            local_modes_added / local_modes_expected if local_modes_expected else 0.0
        ),
        "dft_class_intensities_km_mol": dft_class_intensity,
    }


def _write_or_activate_hybrid_spectrum(
    run_dir: Path, grid: np.ndarray, intrinsic: np.ndarray, *, activate: bool,
) -> dict:
    """Write a hybrid diagnostic and activate it only after every gate passes."""
    from .spectrum import (format_spectrum, write_spectrum_csv,
                           write_spectrum_plot)

    diagnostic = write_spectrum_csv(
        run_dir / "spectrum_hybrid_multifidelity_intrinsic.csv", grid,
        intrinsic, "absorption_strength",
    )
    result = {
        "hybrid_intrinsic_spectrum": str(diagnostic),
        "hybrid_spectrum_role": "displayed" if activate else "diagnostic_not_displayed",
    }
    if not activate:
        return result
    from .experimental_ir import apply_measurement_response, resolve_experimental_profile
    metadata_path = run_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    profile = resolve_experimental_profile(
        phase=str(metadata.get("phase", "liquid")),
        measurement=str(metadata.get("measurement", "auto")),
        instrument_resolution_cm_1=float(metadata.get("instrument_resolution_cm_1", 4.0)),
        apodization=str(metadata.get("apodization", "happ-genzel")),
        residual_fwhm_cm_1=metadata.get("residual_fwhm_cm_1"),
    )
    write_spectrum_csv(run_dir / "spectrum_intrinsic.csv", grid, intrinsic, "absorption_strength")
    write_spectrum_plot(
        run_dir / "spectrum_intrinsic.png", grid,
        intrinsic / max(float(np.max(intrinsic)), 1e-300),
        "relative absorption strength", title="Hybrid multi-fidelity intrinsic IR spectrum",
    )
    measured = apply_measurement_response(grid, intrinsic, profile=profile)
    display, column = format_spectrum(
        measured, style=str(metadata.get("spectrum_style", "transmittance")),
        max_absorbance=float(metadata.get("max_absorbance", 1.0)),
    )
    write_spectrum_csv(run_dir / "spectrum.csv", grid, display, column)
    write_spectrum_plot(
        run_dir / "spectrum.png", grid, display, column.replace("_", " "),
        title="Hybrid multi-fidelity simulated FTIR measurement",
    )
    return result


def _stratified_transfer_order(records: list[dict], variant: int = 0) -> list[dict]:
    """Interleave topology/coordination strata and span each stratum's frequencies."""
    strata: dict[tuple, list[dict]] = {}
    for record in records:
        classes = tuple(sorted({str(mode["mode_class"]) for mode in record.get("modes") or []}))
        features = record.get("environment_features") or {}
        distance = _finite(features.get("h_bond_distance_angstrom"), math.inf)
        association = "strong" if distance < 1.85 else "intermediate" if distance <= 2.20 else "weak"
        key = (str(record.get("topology", "dimer")), association, classes)
        strata.setdefault(key, []).append(record)
    ordered_groups = []
    for key in sorted(strata):
        group = sorted(
            strata[key],
            key=lambda record: (
                float(np.mean([float(mode["freq"]) for mode in record.get("modes") or []])),
                str(record["candidate_id"]),
            ),
        )
        # Alternating low/high frequencies prevents one tail entering only in a late batch.
        spanning = []
        while group:
            spanning.append(group.pop(0))
            if group:
                spanning.append(group.pop(-1))
        if spanning:
            shift = variant % len(spanning)
            spanning = spanning[shift:] + spanning[:shift]
        ordered_groups.append((key, spanning))
    if ordered_groups:
        shift = variant % len(ordered_groups)
        ordered_groups = ordered_groups[shift:] + ordered_groups[:shift]
    output = []
    while any(group for _, group in ordered_groups):
        for _, group in ordered_groups:
            if group:
                output.append(group.pop(0))
    return output


def _bootstrap_band_statistics(
    observations: list[tuple[float, float]], *, seed: int, samples: int = 300,
) -> dict:
    frequencies = np.asarray([item[0] for item in observations], dtype=float)
    if len(frequencies) < 3:
        return {"samples": 0, "center_95_ci_cm-1": None, "fwhm_95_ci_cm-1": None}
    generator = np.random.default_rng(seed)
    indices = generator.integers(0, len(frequencies), size=(samples, len(frequencies)))
    draws = frequencies[indices]
    centers = np.mean(draws, axis=1)
    widths = 2.354820045 * np.std(draws, axis=1)
    center_ci = np.percentile(centers, [2.5, 97.5])
    width_ci = np.percentile(widths, [2.5, 97.5])
    return {
        "samples": samples,
        "center_95_ci_cm-1": [float(value) for value in center_ci],
        "center_95_ci_halfwidth_cm-1": float((center_ci[1] - center_ci[0]) / 2.0),
        "fwhm_95_ci_cm-1": [float(value) for value in width_ci],
        "fwhm_95_ci_halfwidth_cm-1": float((width_ci[1] - width_ci[0]) / 2.0),
    }


def evaluate_transferred_ensemble_convergence(records: Iterable[dict], *, batch_size: int = 10) -> dict:
    """Evaluate candidate-independent distributions under several stratified orders."""
    records = list(records)
    if batch_size < 1:
        raise ValueError("Transfer convergence batch size must be positive")

    def summarize(batch: list[dict], endpoint: int) -> tuple[dict, np.ndarray]:
        grouped: dict[str, list[tuple[float, float]]] = {}
        grid = np.arange(400.0, 4000.0 + 1.0, 1.0)
        spectrum = np.zeros_like(grid)
        candidate_weight = 1.0 / len(batch) if batch else 0.0
        sigma = 6.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        for record in batch:
            divisor = max(1, int(record.get("cluster_size") or 1))
            candidate_modes: dict[str, list[dict]] = {}
            for mode in record.get("modes") or []:
                candidate_modes.setdefault(str(mode["mode_class"]), []).append(mode)
            for mode_class, modes in candidate_modes.items():
                frequency = float(np.mean([float(mode["freq"]) for mode in modes]))
                intensity = sum(max(0.0, _finite(mode.get("intensity"))) for mode in modes) / divisor
                grouped.setdefault(mode_class, []).append((frequency, intensity))
                spectrum += candidate_weight * intensity * np.exp(
                    -0.5 * np.square((grid - frequency) / sigma)
                )
        summary = {}
        for mode_class, observations in grouped.items():
            frequencies = np.asarray([item[0] for item in observations])
            intensities = np.asarray([item[1] for item in observations])
            seed = int(hashlib.sha256(f"{mode_class}:{endpoint}".encode()).hexdigest()[:8], 16)
            summary[mode_class] = {
                "center_cm-1": float(np.mean(frequencies)),
                "frequency_stddev_cm-1": float(np.std(frequencies)),
                "equivalent_fwhm_cm-1": float(2.354820045 * np.std(frequencies)),
                "integrated_intensity": float(np.sum(intensities) / max(1, len(batch))),
                "independent_environments": len(observations),
                **_bootstrap_band_statistics(observations, seed=seed),
            }
        return summary, spectrum

    def compare(previous_summary, previous_spectrum, summary, spectrum, consecutive):
        common = sorted(set(previous_summary) & set(summary))
        new_classes = sorted(set(summary) - set(previous_summary))
        center_changes = [abs(summary[key]["center_cm-1"] - previous_summary[key]["center_cm-1"]) for key in common]
        width_changes = [abs(summary[key]["equivalent_fwhm_cm-1"] - previous_summary[key]["equivalent_fwhm_cm-1"]) for key in common]
        relative = [change / max(previous_summary[key]["equivalent_fwhm_cm-1"], 1e-12) for key, change in zip(common, width_changes)]
        denominator = float(np.linalg.norm(spectrum) * np.linalg.norm(previous_spectrum))
        overlap = float(np.dot(spectrum, previous_spectrum) / denominator) if denominator else 0.0
        passed = bool(common and not new_classes and max(center_changes, default=math.inf) < 5.0
                      and all(absolute < 10.0 or fraction < 0.10 for absolute, fraction in zip(width_changes, relative))
                      and overlap > 0.98)
        consecutive = consecutive + 1 if passed else 0
        return {
            "maximum_center_change_cm-1": max(center_changes, default=None),
            "maximum_width_change_cm-1": max(width_changes, default=None),
            "whole_spectrum_overlap": overlap,
            "new_band_classes": new_classes,
            "passed": passed,
            "consecutive_passes": consecutive,
        }, consecutive

    endpoints = list(range(min(batch_size, len(records)), len(records) + 1, batch_size))
    if records and (not endpoints or endpoints[-1] != len(records)):
        endpoints.append(len(records))
    ordering_reports = []
    for variant in range(min(3, max(1, len(records)))):
        ordered = _stratified_transfer_order(records, variant)
        batches, previous_summary, previous_spectrum, consecutive = [], None, None, 0
        for endpoint in endpoints:
            summary, spectrum = summarize(ordered[:endpoint], endpoint)
            comparison = None
            if previous_summary is not None:
                comparison, consecutive = compare(
                    previous_summary, previous_spectrum, summary, spectrum, consecutive,
                )
            batches.append({"environments": endpoint, "bands": summary, "comparison": comparison})
            previous_summary, previous_spectrum = summary, spectrum
        ordering_reports.append({
            "variant": variant,
            "candidate_order": [record["candidate_id"] for record in ordered],
            "converged": consecutive >= 2,
            "consecutive_passing_comparisons": consecutive,
            "batches": batches,
        })
    ordering_robust = bool(
        ordering_reports and all(report["converged"] for report in ordering_reports)
    )
    acquisition_groups: dict[int, list[dict]] = {}
    for record in records:
        acquisition_groups.setdefault(int(record.get("acquisition_round", 0)), []).append(record)
    acquisition_batches = []
    previous_summary, previous_spectrum, acquisition_consecutive = None, None, 0
    cumulative = []
    for acquisition_round in sorted(acquisition_groups):
        cumulative.extend(sorted(
            acquisition_groups[acquisition_round], key=lambda item: str(item["candidate_id"]),
        ))
        summary, spectrum = summarize(cumulative, len(cumulative))
        comparison = None
        if previous_summary is not None:
            comparison, acquisition_consecutive = compare(
                previous_summary, previous_spectrum, summary, spectrum,
                acquisition_consecutive,
            )
        acquisition_batches.append({
            "acquisition_round": acquisition_round,
            "environments": len(cumulative),
            "bands": summary,
            "comparison": comparison,
        })
        previous_summary, previous_spectrum = summary, spectrum
    acquisition_converged = bool(
        len(acquisition_batches) >= 3 and acquisition_consecutive >= 2
    )
    final_acquisition_bands = (
        acquisition_batches[-1]["bands"] if acquisition_batches else {}
    )
    bootstrap_precision_pass = bool(final_acquisition_bands) and all(
        band.get("center_95_ci_halfwidth_cm-1") is not None
        and band.get("fwhm_95_ci_halfwidth_cm-1") is not None
        and float(band["center_95_ci_halfwidth_cm-1"]) <= 15.0
        and float(band["fwhm_95_ci_halfwidth_cm-1"]) <= 25.0
        for band in final_acquisition_bands.values()
    )
    converged = acquisition_converged and bootstrap_precision_pass
    return {
        "schema_version": TRANSFER_SCHEMA_VERSION,
        "kind": "adaptive_dft_transferred_environment_convergence",
        "ordering_model": "three_deterministic_topology_coordination_frequency_stratified_orders_as_sensitivity_diagnostic",
        "statistical_unit": "one_candidate_environment_per_mode_class",
        "bootstrap_model": "deterministic_candidate_resampling_300_draws",
        "thresholds": {
            "center_change_cm-1": 5.0, "width_change_cm-1": 10.0,
            "relative_width_change": 0.10, "whole_spectrum_overlap": 0.98,
            "required_consecutive_batches": 2, "required_orderings_passed": len(ordering_reports),
            "maximum_center_95_ci_halfwidth_cm-1": 15.0,
            "maximum_fwhm_95_ci_halfwidth_cm-1": 25.0,
        },
        "status": "converged" if converged else "sample_limit_reached_unconverged",
        "converged": converged,
        "environments_evaluated": len(records),
        "batches": ordering_reports[0]["batches"] if ordering_reports else [],
        "ordering_diagnostics": ordering_reports,
        "ordering_sensitivity_converged": ordering_robust,
        "bootstrap_precision_pass": bootstrap_precision_pass,
        "acquisition_round_diagnostics": {
            "model": "real_cumulative_balanced_acquisition_rounds",
            "required_consecutive_comparisons": 2,
            "converged": acquisition_converged,
            "consecutive_passing_comparisons": acquisition_consecutive,
            "batches": acquisition_batches,
        },
    }


def build_frequency_transfer_artifacts(run_dir: Path) -> dict:
    """Build transfer validation and diagnostic spectra from completed cached jobs."""
    run_dir = Path(run_dir)
    xtb_path = run_dir / "xtb-snapshot-frequencies.json"
    selected_path = run_dir / "clusters" / "selected-conformers.json"
    conformers_path = run_dir / "clusters" / "conformers.json"
    if not all(path.is_file() for path in (xtb_path, selected_path, conformers_path)):
        return {"status": "missing_frequency_transfer_inputs"}
    xtb_payload = json.loads(xtb_path.read_text())
    xtb_records = xtb_payload.get("candidates") or []
    xtb_by_id = {item["candidate_id"]: item for item in xtb_records}
    selected = sorted(
        json.loads(selected_path.read_text()).get("conformers") or [],
        key=lambda item: int(item["selected_position"]),
    )
    conformers = json.loads(conformers_path.read_text())
    pairs = []
    for representative, conformer in zip(selected, conformers):
        candidate_id = representative.get("source_xtb_candidate_id")
        frequency_output = conformer.get("frequency_output")
        hessian = Path(frequency_output).with_suffix(".hess") if frequency_output else None
        if candidate_id in xtb_by_id and hessian and hessian.is_file():
            from src.parser import parse_orca_ir
            paired_representative = {
                **representative,
                # Older selection manifests predate local coordination labels;
                # the completed conformer manifest backfills them from the
                # retained H-bond interaction identities.
                "local_stretch_bonds": (
                    conformer.get("local_stretch_bonds")
                    or representative.get("local_stretch_bonds") or []
                ),
                "environment_features": (
                    conformer.get("environment_features")
                    or representative.get("environment_features") or {}
                ),
            }
            representative_pairs = pair_representative_local_modes(
                paired_representative, xtb_by_id[candidate_id], hessian,
                parse_orca_ir(frequency_output),
            )
            local_mode_path = conformer.get("local_mode_output")
            if not local_mode_path:
                local_mode_path = (
                    Path(frequency_output).parent
                    / "local-mode-finite-differences" / "local-modes.json"
                )
            local_mode_path = Path(local_mode_path)
            if local_mode_path.is_file():
                representative_pairs = apply_finite_difference_pair_overrides(
                    representative_pairs, json.loads(local_mode_path.read_text()),
                    exclude_unvalidated_snapshot_pairs=(
                        (conformer.get("snapshot_hessian_reliability") or {}).get(
                            "stationarity_status"
                        ) == "poor"
                    ),
                )
            pairs.extend(representative_pairs)
    validation = apply_ensemble_transfer_coverage(
        validate_frequency_transfer(pairs), xtb_records, pairs,
    )
    acquisition_path = run_dir / "environment-acquisition.json"
    if acquisition_path.is_file():
        acquisition = json.loads(acquisition_path.read_text())
        validation["environment_acquisition"] = {
            "artifact": str(acquisition_path),
            "status": acquisition.get("status"),
            "representative_budget": acquisition.get("representative_budget"),
            "remaining_class_deficits": acquisition.get("remaining_class_deficits"),
        }
    validation_path = run_dir / "frequency-transfer-validation.json"
    validation_path.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    transferred = transfer_snapshot_modes(xtb_records, pairs, validation)
    transferred_by_id = {item["candidate_id"]: item for item in transferred}
    validated_classes = {
        key for key, item in validation.get("classes", {}).items()
        if item.get("status") == "validated"
    }
    out_of_domain_candidates = []
    for record in xtb_records:
        expected = {
            _class_key(mode) for mode in record.get("modes") or []
            if mode.get("bond_class") and _class_key(mode) in validated_classes
        }
        covered = {
            str(mode["mode_class"])
            for mode in (transferred_by_id.get(record["candidate_id"], {}).get("modes") or [])
        }
        if expected - covered:
            out_of_domain_candidates.append({
                "candidate_id": record["candidate_id"],
                "missing_mode_classes": sorted(expected - covered),
                "acquisition_round": int(record.get("acquisition_round", 0)),
            })
    transfer_convergence = (
        evaluate_transferred_ensemble_convergence(transferred)
        if validation["status"] == "validated"
        else {
            "schema_version": TRANSFER_SCHEMA_VERSION,
            "kind": "adaptive_dft_transferred_environment_convergence",
            "status": "insufficient_dft_transfer_validation",
            "converged": False,
            "environments_evaluated": len(transferred),
            "batches": [],
        }
    )
    convergence_path = run_dir / "frequency-transfer-convergence.json"
    transfer_convergence["out_of_domain_candidates"] = out_of_domain_candidates
    transfer_convergence["additional_orca_representatives_requested"] = bool(
        out_of_domain_candidates
    )
    convergence_path.write_text(
        json.dumps(transfer_convergence, indent=2, sort_keys=True) + "\n"
    )
    environment_convergence_path = run_dir / "environment-convergence.json"
    if environment_convergence_path.is_file():
        environment_convergence = json.loads(environment_convergence_path.read_text())
        environment_convergence["dft_transferred_ensemble"] = transfer_convergence
        environment_convergence_path.write_text(
            json.dumps(environment_convergence, indent=2, sort_keys=True) + "\n"
        )
    metadata_path = run_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    scale_factor = _finite(metadata.get("scale_factor"), 0.97)
    hybrid_grid, hybrid_intrinsic, hybrid_provenance = _hybrid_multifidelity_intrinsic(
        run_dir, transferred, pairs, selected, conformers, scale_factor=scale_factor,
    )
    activate_hybrid = bool(
        validation["status"] == "validated"
        and transfer_convergence.get("converged") is True
        and hybrid_provenance["local_modes_expected"] > 0
        and hybrid_provenance["local_intensity_coverage_fraction"] >= 1.0
    )
    hybrid_artifacts = _write_or_activate_hybrid_spectrum(
        run_dir, hybrid_grid, hybrid_intrinsic, activate=activate_hybrid,
    )
    if activate_hybrid:
        display_basis = "hybrid_multifidelity_ensemble"
        use_status = "validated_converged_hybrid_displayed"
        display_frequency_source = hybrid_provenance["frequency_source"]
        display_intensity_source = (
            "representative_orca_nonlocal_and_mode_class_local_intensities"
        )
        display_population_source = hybrid_provenance["population_source"]
        display_width_source = "converged_dft_corrected_xtb_environment_distribution"
        display_switch_gate = "passed"
    else:
        display_basis = "representative_orca_ensemble"
        use_status = (
            "validated_frequency_diagnostic_not_displayed"
            if validation["status"] == "validated"
            else "insufficient_dft_transfer_validation"
        )
        display_frequency_source = "representative_orca_harmonic_frequencies"
        display_intensity_source = "representative_orca_ir_intensities"
        display_population_source = "stratified_equal_representative_weights"
        display_width_source = "representative_orca_environment_distribution"
        display_switch_gate = (
            "blocked_until_validation_convergence_and_local_intensity_coverage_pass"
        )
    transfer = {
        "schema_version": TRANSFER_SCHEMA_VERSION,
        "kind": "representative_additive_dft_frequency_transfer",
        "formula": "nu_estimated = nu_xtb_snapshot + (nu_dft_representative - nu_xtb_representative)",
        "status": validation["status"],
        "display_basis": display_basis,
        "transfer_use_status": use_status,
        "frequency_source": "dft_corrected_xtb_snapshot_frequencies",
        "intensity_source": "raw_xtb_diagnostic_only",
        "display_frequency_source": display_frequency_source,
        "display_intensity_source": display_intensity_source,
        "display_population_source": display_population_source,
        "display_width_source": display_width_source,
        "display_switch_gate": display_switch_gate,
        "intensity_transfer_status": (
            "representative_mode_class_intensities_applied"
            if activate_hybrid else "diagnostic_pending_transfer_convergence"
        ),
        "population_model": "stratified_equal_candidate_weights_not_liquid_populations",
        "representative_pairs": pairs,
        "transferred_candidates": transferred,
        "validation_artifact": str(validation_path),
        "convergence_artifact": str(convergence_path),
        "hybrid_multifidelity": {
            **hybrid_provenance,
            **hybrid_artifacts,
            "activation_gate_passed": activate_hybrid,
        },
    }
    transfer_path = run_dir / "frequency-transfer.json"
    transfer_path.write_text(json.dumps(transfer, indent=2, sort_keys=True) + "\n")
    raw_records = [
        {
            "cluster_size": int(item.get("cluster_size") or 1),
            "modes": [mode for mode in item.get("modes") or [] if mode.get("bond_class")],
        }
        for item in xtb_records if item.get("frequency_status") == "completed"
    ]
    raw_path = _write_spectrum(run_dir / "spectrum_xtb_environment.csv", raw_records, raw=False)
    corrected_path = _write_spectrum(
        run_dir / "spectrum_dft_transferred_intrinsic.csv", transferred, raw=False,
    )
    return {
        "status": validation["status"], "frequency_transfer": transfer_path,
        "frequency_transfer_validation": validation_path,
        "frequency_transfer_convergence": convergence_path,
        "xtb_environment_spectrum": raw_path,
        "dft_transferred_intrinsic_spectrum": corrected_path,
        "dft_transferred_intrinsic_spectrum_role": "diagnostic_not_displayed",
        **hybrid_artifacts,
        "frequency_transfer_use_status": use_status,
        "display_basis": display_basis,
        "transferred_candidates": len(transferred),
    }
