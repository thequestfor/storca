"""Condition-aware conversion of molecular IR absorption into an FTIR-like trace.

This module deliberately keeps the molecular ensemble and measurement transfer
functions separate.  It does not invent solvent shifts, aggregation, or
vibrational lifetimes that were not present in the calculated ensemble.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from pathlib import Path
from typing import Iterable

import numpy as np


SUPPORTED_PHASES = {"gas", "solution", "liquid", "solid"}
SUPPORTED_MEASUREMENTS = {"auto", "atr", "transmission", "gas-cell"}
SUPPORTED_APODIZATIONS = {"gaussian", "triangular", "happ-genzel"}


@dataclass(frozen=True)
class ExperimentalIRProfile:
    phase: str
    measurement: str
    instrument_resolution_cm_1: float
    apodization: str
    residual_fwhm_cm_1: float


@dataclass(frozen=True)
class EnvironmentSufficiencyConfig:
    """Versioned minimum evidence required to expose an environment width."""

    schema_version: int = 1
    minimum_independent_environments: int = 3
    minimum_effective_sample_size: float = 2.5
    minimum_frequency_variance_cm_2: float = 1.0
    minimum_mode_overlap: float = 0.45
    minimum_overlap_covered_population: float = 0.8
    minimum_h_bond_distance_span_angstrom: float = 0.15
    minimum_h_bond_angle_span_degrees: float = 15.0
    minimum_orientation_span_degrees: float = 30.0
    minimum_heavy_atom_rmsd_span_angstrom: float = 0.25


def environment_sufficiency_metadata(
    config: EnvironmentSufficiencyConfig | None = None,
) -> dict:
    return asdict(config or EnvironmentSufficiencyConfig())


def resolve_experimental_profile(
    *,
    phase: str,
    measurement: str,
    instrument_resolution_cm_1: float,
    apodization: str,
    residual_fwhm_cm_1: float | None,
) -> ExperimentalIRProfile:
    """Validate and resolve a small, explicit measurement contract."""
    phase = phase.strip().lower()
    measurement = measurement.strip().lower()
    apodization = apodization.strip().lower()
    if phase not in SUPPORTED_PHASES:
        raise ValueError(f"Unsupported phase '{phase}'; choose from {', '.join(sorted(SUPPORTED_PHASES))}")
    if measurement not in SUPPORTED_MEASUREMENTS:
        raise ValueError(
            f"Unsupported measurement '{measurement}'; choose from {', '.join(sorted(SUPPORTED_MEASUREMENTS))}"
        )
    if apodization not in SUPPORTED_APODIZATIONS:
        raise ValueError(
            f"Unsupported apodization '{apodization}'; choose from {', '.join(sorted(SUPPORTED_APODIZATIONS))}"
        )
    if not math.isfinite(instrument_resolution_cm_1) or instrument_resolution_cm_1 <= 0:
        raise ValueError("Instrument resolution must be finite and greater than zero")
    if measurement == "auto":
        measurement = {"gas": "gas-cell", "solution": "transmission", "liquid": "atr", "solid": "atr"}[phase]
    if measurement == "gas-cell" and phase != "gas":
        raise ValueError("gas-cell measurement requires --phase gas")
    if phase == "gas" and measurement == "atr":
        raise ValueError("ATR measurement is incompatible with --phase gas")

    # This is deliberately a small residual width.  Conformer/environmental
    # frequency distributions should provide the chemical broadening.
    default_widths = {"gas": 1.0, "solution": 4.0, "liquid": 6.0, "solid": 8.0}
    residual = default_widths[phase] if residual_fwhm_cm_1 is None else float(residual_fwhm_cm_1)
    if not math.isfinite(residual) or residual <= 0:
        raise ValueError("Residual FWHM must be finite and greater than zero")
    return ExperimentalIRProfile(
        phase=phase,
        measurement=measurement,
        instrument_resolution_cm_1=float(instrument_resolution_cm_1),
        apodization=apodization,
        residual_fwhm_cm_1=residual,
    )


def _instrument_kernel(step: float, resolution: float, apodization: str, size: int) -> np.ndarray:
    half_width = min(max(4, int(math.ceil(8.0 * resolution / step))), max(1, (size - 1) // 2))
    offsets = np.arange(-half_width, half_width + 1, dtype=float) * step
    if apodization == "gaussian":
        sigma = resolution / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        kernel = np.exp(-0.5 * np.square(offsets / sigma))
    else:
        # A triangularly apodized interferogram has a sinc-squared instrument
        # line shape.  The zero spacing below makes its central-lobe FWHM equal
        # to the requested nominal resolution.
        zero_spacing = resolution / 0.886
        triangular = np.square(np.sinc(offsets / zero_spacing))
        if apodization == "triangular":
            kernel = triangular
        else:
            # Happ-Genzel is represented as a documented positive central-lobe
            # approximation.  It suppresses sinc sidelobes more strongly while
            # avoiding a false claim of emulating a particular vendor firmware.
            sigma = resolution / (2.0 * math.sqrt(2.0 * math.log(2.0)))
            gaussian = np.exp(-0.5 * np.square(offsets / sigma))
            kernel = 0.65 * gaussian + 0.35 * triangular
    total = float(kernel.sum())
    if total <= 0 or not math.isfinite(total):
        raise RuntimeError("Could not construct a finite instrument line-shape kernel")
    return kernel / total


def apply_measurement_response(
    grid: np.ndarray,
    absorbance_strength: np.ndarray,
    *,
    profile: ExperimentalIRProfile,
) -> np.ndarray:
    """Apply measurement geometry and finite FTIR resolution to absorption."""
    grid = np.asarray(grid, dtype=float)
    signal = np.asarray(absorbance_strength, dtype=float)
    if (grid.ndim != 1 or signal.ndim != 1 or len(grid) != len(signal) or len(grid) < 3
            or not np.all(np.isfinite(grid)) or not np.all(np.isfinite(signal))
            or np.any(signal < 0) or np.any(np.diff(grid) <= 0)):
        raise ValueError("Measurement response requires equal, finite, increasing spectrum arrays")
    steps = np.diff(grid)
    step = float(np.median(steps))
    if not np.allclose(steps, step, rtol=1e-5, atol=1e-8):
        raise ValueError("Measurement response requires a uniformly spaced wavenumber grid")

    if profile.measurement == "atr":
        # First-order ATR penetration-depth dependence.  Refractive indices and
        # incidence angle are not inferred; a vendor-specific ATR correction is
        # therefore intentionally outside this profile.
        reference = 2000.0
        signal = signal * reference / grid

    kernel = _instrument_kernel(
        step, profile.instrument_resolution_cm_1, profile.apodization, len(grid),
    )
    measured = np.convolve(signal, kernel, mode="same")
    return np.maximum(measured, 0.0)


def _feature_value(features: dict, key: str) -> float | None:
    value = features.get(key)
    if isinstance(value, (int, float)) and math.isfinite(value):
        return float(value)
    return None


def _feature_span(members: list[dict], key: str) -> float | None:
    values = [
        value for member in members
        if (value := _feature_value(member.get("environment_features") or {}, key)) is not None
    ]
    return max(values) - min(values) if len(values) == len(members) and values else None


def _circular_span_degrees(members: list[dict], key: str) -> float | None:
    values = [
        value % 360.0 for member in members
        if (value := _feature_value(member.get("environment_features") or {}, key)) is not None
    ]
    if len(values) != len(members) or not values:
        return None
    return max(
        min(abs(left - right), 360.0 - abs(left - right))
        for left in values for right in values
    )


def evaluate_environment_sufficiency(
    members: Iterable[dict],
    *,
    config: EnvironmentSufficiencyConfig | None = None,
) -> dict:
    """Fail closed unless a matched band has independently diverse evidence."""
    resolved = config or EnvironmentSufficiencyConfig()
    member_list = list(members)
    population = sum(float(member["weight"]) for member in member_list)
    normalized_weights = [float(member["weight"]) / population for member in member_list] if population > 0 else []
    environment_weight_map: dict[str, float] = {}
    geometry_by_environment: dict[str, dict] = {}
    for member in member_list:
        environment_id = str(member.get("independent_environment_id", member.get("position")))
        environment_weight_map[environment_id] = environment_weight_map.get(environment_id, 0.0) + float(member["weight"])
        geometry_by_environment.setdefault(environment_id, member)
    environment_weight_total = sum(environment_weight_map.values())
    environment_weights = [
        weight / environment_weight_total for weight in environment_weight_map.values()
    ] if environment_weight_total > 0 else []
    effective_sample_size = (
        1.0 / sum(weight * weight for weight in environment_weights)
        if environment_weights else 0.0
    )
    environment_ids = set(environment_weight_map)
    frequencies = [float(member["frequency"]) for member in member_list]
    mean = (
        sum(frequency * weight for frequency, weight in zip(frequencies, normalized_weights))
        if normalized_weights else 0.0
    )
    variance = (
        sum(weight * (frequency - mean) ** 2 for frequency, weight in zip(frequencies, normalized_weights))
        if normalized_weights else 0.0
    )
    overlap_population = sum(
        float(member["weight"]) for member in member_list
        if member.get("overlap") is not None
    )
    overlap_coverage = overlap_population / population if population > 0 else 0.0
    overlaps = [float(member["overlap"]) for member in member_list if member.get("overlap") is not None]
    minimum_overlap = min(overlaps) if overlaps else None

    geometry_members = list(geometry_by_environment.values())
    distance_span = _feature_span(geometry_members, "h_bond_distance_angstrom")
    angle_span = _feature_span(geometry_members, "donor_h_acceptor_angle_degrees")
    orientation_span = _circular_span_degrees(geometry_members, "intermolecular_orientation_degrees")
    rmsd_span = _feature_span(geometry_members, "heavy_atom_rmsd_angstrom")
    cluster_ids = {
        str((member.get("environment_features") or {}).get("geometry_cluster_id"))
        for member in geometry_members
        if (member.get("environment_features") or {}).get("geometry_cluster_id") is not None
    }
    distinct_geometry_clusters = len(cluster_ids) if len(cluster_ids) == len(geometry_members) else 0
    secondary_geometry_diversity = any((
        angle_span is not None and angle_span >= resolved.minimum_h_bond_angle_span_degrees,
        orientation_span is not None and orientation_span >= resolved.minimum_orientation_span_degrees,
        rmsd_span is not None and rmsd_span >= resolved.minimum_heavy_atom_rmsd_span_angstrom,
    ))
    geometric_diversity = bool(
        distinct_geometry_clusters >= resolved.minimum_independent_environments
        and distance_span is not None
        and distance_span >= resolved.minimum_h_bond_distance_span_angstrom
        and secondary_geometry_diversity
    )

    failures: list[str] = []
    if len(environment_ids) < resolved.minimum_independent_environments:
        failures.append("fewer_than_minimum_independent_environments")
    if variance < resolved.minimum_frequency_variance_cm_2:
        failures.append("frequency_variance_below_minimum")
    if effective_sample_size < resolved.minimum_effective_sample_size:
        failures.append("effective_sample_size_below_minimum")
    if not geometric_diversity:
        failures.append("insufficient_or_unavailable_geometric_diversity")
    if overlap_coverage < resolved.minimum_overlap_covered_population:
        failures.append("mode_overlap_population_below_minimum")
    if minimum_overlap is None or minimum_overlap < resolved.minimum_mode_overlap:
        failures.append("mode_matching_confidence_below_minimum")

    sufficient = not failures
    sampled_fwhm = 2.0 * math.sqrt(2.0 * math.log(2.0)) * math.sqrt(max(variance, 0.0))
    return {
        "width_status": "calculated_environment_width" if sufficient else "insufficient_environment_sampling",
        "calculated_environment_fwhm_cm-1": sampled_fwhm if sufficient else 0.0,
        "display_width_source": (
            "environment_distribution_plus_residual_plus_instrument"
            if sufficient else "residual_plus_instrument"
        ),
        "environments": len(environment_ids),
        "effective_sample_size": effective_sample_size,
        "effective_sample_size_basis": "independent_environment_aggregated_weights",
        "frequency_variance_cm-2": variance,
        "sampled_environment_fwhm_cm-1": sampled_fwhm,
        "mode_overlap_covered_population": overlap_coverage,
        "minimum_displacement_overlap": minimum_overlap,
        "geometric_diversity_sufficient": geometric_diversity,
        "geometric_diversity": {
            "distinct_geometry_clusters": distinct_geometry_clusters,
            "h_bond_distance_span_angstrom": distance_span,
            "donor_h_acceptor_angle_span_degrees": angle_span,
            "intermolecular_orientation_span_degrees": orientation_span,
            "heavy_atom_rmsd_span_angstrom": rmsd_span,
        },
        "sufficiency_failures": failures,
    }


def _ensemble_band_analysis(
    conformers: Iterable[dict],
    *,
    scale_factor: float,
    environment_sampling: bool = False,
    sufficiency_config: EnvironmentSufficiencyConfig | None = None,
) -> tuple[list[dict], dict[tuple[int, int], int]]:
    """Return matched band statistics and each source mode's band assignment."""
    conformer_list = list(conformers)
    normal_mode_sets: dict[int, object] = {}
    for position, conformer in enumerate(conformer_list):
        output = conformer.get("frequency_output")
        if not output:
            continue
        try:
            from .ir_modes import parse_orca_normal_modes
            hessian = Path(output).with_suffix(".hess")
            if hessian.is_file():
                normal_mode_sets[position] = parse_orca_normal_modes(hessian)
        except (OSError, ValueError):
            continue
    mode_characters: dict[int, dict[int, dict]] = {}
    if normal_mode_sets:
        from .mode_character import mode_character_fingerprints
        for position, mode_set in normal_mode_sets.items():
            try:
                fingerprints = mode_character_fingerprints(mode_set)
            except (ValueError, np.linalg.LinAlgError):
                continue
            mode_characters[position] = {
                int(fingerprint["mode"]): fingerprint for fingerprint in fingerprints
            }
    reference_position = max(
        normal_mode_sets,
        key=lambda position: float(conformer_list[position].get("weight", 0.0)),
        default=None,
    )
    correspondence: dict[int, dict[int, dict]] = {}
    local_weight_divisors: dict[tuple[int, int], int] = {}
    mixed_cluster_sizes = len({
        int(conformer.get("cluster_size") or 2) for conformer in conformer_list
    }) > 1
    local_matching = bool(mixed_cluster_sizes and normal_mode_sets)
    if local_matching:
        from .ir_modes import local_stretch_mode_assignments
        for position, mode_set in normal_mode_sets.items():
            bonds = conformer_list[position].get("local_stretch_bonds") or []
            assignments_for_position = local_stretch_mode_assignments(mode_set, bonds)
            counts: dict[str, int] = {}
            for item in assignments_for_position:
                band_key = ":".join((
                    str(item["bond_class"]),
                    str(item.get("spectral_band_class", "unclassified")),
                    str(item.get("coordination_class", "unclassified")),
                ))
                counts[band_key] = counts.get(band_key, 0) + 1
            correspondence[position] = {}
            for item in assignments_for_position:
                mode = int(item["mode"])
                band_key = ":".join((
                    str(item["bond_class"]),
                    str(item.get("spectral_band_class", "unclassified")),
                    str(item.get("coordination_class", "unclassified")),
                ))
                correspondence[position][mode] = {
                    "reference_mode": f"local:{band_key}",
                    "displacement_overlap": item["local_stretch_projection"],
                    "confidence": item["confidence"],
                    "matching_basis": "target_local_xh_projection",
                    "spectral_band_class": item.get("spectral_band_class", "unclassified"),
                    "coordination_class": item.get("coordination_class", "unclassified"),
                    "mode_character": mode_characters.get(position, {}).get(mode),
                }
                local_weight_divisors[(position, mode)] = counts[band_key]
    elif reference_position is not None:
        from .ir_modes import match_normal_modes
        from .mode_character import match_mode_fingerprints
        reference = normal_mode_sets[reference_position]
        for position, mode_set in normal_mode_sets.items():
            try:
                reference_fingerprints = list(mode_characters[reference_position].values())
                candidate_fingerprints = list(mode_characters[position].values())
                matches = match_mode_fingerprints(reference_fingerprints, candidate_fingerprints)
            except (KeyError, ValueError):
                try:
                    matches = match_normal_modes(reference, mode_set)
                except ValueError:
                    continue
            correspondence[position] = {
                item["candidate_mode"]: item for item in matches
            }

    grouped: dict[int, list[dict]] = {}
    for position, conformer in enumerate(conformer_list):
        weight = conformer.get("weight")
        if not isinstance(weight, (int, float)) or not math.isfinite(weight) or weight <= 0:
            continue
        for mode in conformer.get("ir_modes", []):
            try:
                candidate_mode = int(mode["mode"])
                frequency = float(mode["freq"]) * scale_factor
                intensity = float(mode["intensity"])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(frequency) and math.isfinite(intensity) and frequency > 0 and intensity >= 0:
                match = correspondence.get(position, {}).get(candidate_mode)
                if local_matching and match is None:
                    continue
                mode_index = match["reference_mode"] if match else candidate_mode
                if not isinstance(mode_index, str):
                    mode_index = int(mode_index)
                member_weight = float(weight) / local_weight_divisors.get((position, candidate_mode), 1)
                grouped.setdefault(mode_index, []).append({
                    "position": position,
                    "candidate_mode": candidate_mode,
                    "frequency": frequency,
                    "intensity": intensity,
                    "weight": member_weight,
                    "independent_environment_id": conformer.get("independent_environment_id", position),
                    "environment_features": conformer.get("environment_features"),
                    "matching_basis": (
                        match.get("matching_basis", "normal_mode_overlap") if match else "orca_mode_index"
                    ),
                    "overlap": match.get("displacement_overlap") if match else None,
                    "confidence": match.get("confidence") if match else "fallback",
                    "topology": conformer.get("topology", "dimer"),
                    "cluster_size": int(conformer.get("cluster_size") or 2),
                    "spectral_band_class": match.get("spectral_band_class") if match else None,
                    "coordination_class": match.get("coordination_class") if match else None,
                    "mode_character": (
                        match.get("mode_character") if match
                        else mode_characters.get(position, {}).get(candidate_mode)
                    ),
                })

    statistics: list[dict] = []
    assignments: dict[tuple[int, int], int] = {}
    for mode_index, members in sorted(grouped.items(), key=lambda item: str(item[0])):
        population = sum(member["weight"] for member in members)
        if population <= 0:
            continue
        mean = sum(member["frequency"] * member["weight"] for member in members) / population
        variance = sum(member["weight"] * (member["frequency"] - mean) ** 2 for member in members) / population
        intensity = sum(member["weight"] * member["intensity"] for member in members) / population
        overlaps = [member["overlap"] for member in members if member["overlap"] is not None]
        bases = {member["matching_basis"] for member in members}
        for member in members:
            assignments[(member["position"], member["candidate_mode"])] = mode_index
        statistic = {
            "mode": mode_index,
            "conformers": len(members),
            "covered_population": population,
            "mean_frequency_cm-1": mean,
            "frequency_stddev_cm-1": math.sqrt(max(variance, 0.0)),
            "ensemble_fwhm_equivalent_cm-1": 2.0 * math.sqrt(2.0 * math.log(2.0)) * math.sqrt(max(variance, 0.0)),
            "population_weighted_intensity": intensity,
            "matching_basis": (
                next(iter(bases)) if len(bases) == 1
                else "mixed_or_mode_index_fallback"
            ),
            "minimum_displacement_overlap": min(overlaps) if overlaps else None,
            "matching_confidence": sorted({member["confidence"] for member in members}),
            "contributing_topologies": {
                topology: len({member["independent_environment_id"] for member in members if member["topology"] == topology})
                for topology in sorted({member["topology"] for member in members})
            },
            "contributing_cluster_sizes": sorted({member["cluster_size"] for member in members}),
            "spectral_band_class": next(iter({member["spectral_band_class"] for member in members}), None),
            "coordination_classes": sorted({member["coordination_class"] for member in members if member["coordination_class"]}),
        }
        representative_member = max(members, key=lambda member: member["weight"])
        if representative_member.get("mode_character"):
            statistic["mode_character"] = representative_member["mode_character"]
        if statistic["matching_basis"] == "target_local_xh_projection":
            statistic["local_mode_reliability"] = (
                "usable" if overlaps and min(overlaps) >= (sufficiency_config or EnvironmentSufficiencyConfig()).minimum_mode_overlap
                else "low_confidence"
            )
        if environment_sampling:
            decision = evaluate_environment_sufficiency(members, config=sufficiency_config)
            statistic.update(decision)
            statistic.update({
                "center_cm-1": mean,
                "environment_fwhm_cm-1": decision["calculated_environment_fwhm_cm-1"],
                "converged": False,
                "convergence_status": "adaptive_convergence_not_yet_evaluated",
            })
        statistics.append(statistic)
    return statistics, assignments


def ensemble_band_statistics(
    conformers: Iterable[dict],
    *,
    scale_factor: float,
    environment_sampling: bool = False,
    sufficiency_config: EnvironmentSufficiencyConfig | None = None,
) -> list[dict]:
    """Report matched frequency distributions and optional environment gates."""
    statistics, _ = _ensemble_band_analysis(
        conformers,
        scale_factor=scale_factor,
        environment_sampling=environment_sampling,
        sufficiency_config=sufficiency_config,
    )
    return statistics


def apply_environment_convergence_status(
    bands: Iterable[dict], convergence_report: dict | None,
) -> list[dict]:
    """Combine minimum-width sufficiency with the adaptive stopping decision."""
    copied = [{**band} for band in bands]
    summary = (convergence_report or {}).get("summary") or {}
    convergence_status = str(summary.get("status", "adaptive_convergence_not_evaluated"))
    converged = bool(summary.get("converged", False))
    converged_modes = set(summary.get("converged_important_modes") or [])
    for band in copied:
        if band.get("width_status") == "insufficient_environment_sampling":
            band["converged"] = False
            band["convergence_status"] = "insufficient_environment_sampling"
        elif converged and band.get("mode") in converged_modes:
            band["converged"] = True
            band["convergence_status"] = "adaptive_environment_converged"
        else:
            band["width_status"] = "environment_width_unconverged"
            band["converged"] = False
            band["convergence_status"] = (
                convergence_status if not converged else "not_evaluated_as_important_band"
            )
            band["display_width_source"] = (
                "provisional_environment_distribution_plus_residual_plus_instrument"
            )
    return copied


def collapse_insufficient_environment_bands(
    conformers: Iterable[dict],
    *,
    scale_factor: float,
    sufficiency_config: EnvironmentSufficiencyConfig | None = None,
) -> tuple[list[dict], list[dict]]:
    """Center failed distributions so only residual/instrument width is shown."""
    if not math.isfinite(scale_factor) or scale_factor <= 0:
        raise ValueError("Scale factor must be finite and greater than zero")
    conformer_list = list(conformers)
    statistics, assignments = _ensemble_band_analysis(
        conformer_list,
        scale_factor=scale_factor,
        environment_sampling=True,
        sufficiency_config=sufficiency_config,
    )
    bands = {band["mode"]: band for band in statistics}
    collapsed: list[dict] = []
    for position, conformer in enumerate(conformer_list):
        copied_modes = []
        for mode in conformer.get("ir_modes", []):
            copied_mode = {**mode}
            try:
                key = (position, int(mode["mode"]))
            except (KeyError, TypeError, ValueError):
                copied_modes.append(copied_mode)
                continue
            band = bands.get(assignments.get(key))
            if band and band.get("width_status") == "insufficient_environment_sampling":
                copied_mode["freq"] = float(band["mean_frequency_cm-1"]) / scale_factor
            copied_modes.append(copied_mode)
        collapsed.append({**conformer, "ir_modes": copied_modes})
    return collapsed, statistics


def experimental_profile_metadata(profile: ExperimentalIRProfile) -> dict:
    return {
        "name": "calculation-first-ftir-v1",
        "kind": "ensemble_plus_measurement_transfer",
        **asdict(profile),
        "atr_model": "first_order_inverse_wavenumber_penetration_depth" if profile.measurement == "atr" else None,
        "noise_or_baseline_added": False,
        "limitations": (
            "Conformer broadening is calculated, but aggregation, solvent configurations, anharmonic lifetimes, "
            "and vendor-specific instrument response are not inferred. Happ-Genzel is a central-lobe approximation."
        ),
    }
