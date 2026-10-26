"""Auditable, bounded diagnosis for unresolved reaction surfaces.

This module deliberately contains no electronic-structure execution.  It
turns retained path evidence into a deterministic *diagnosis* and a bounded
strategy contract.  A diagnosis is allocation evidence only: it never proves
that a reaction, transition state, intermediate, rate, or lifetime exists.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


SURFACE_RESOLUTION_SCHEMA_VERSION = 1
SURFACE_RESOLUTION_POLICY_ID = "generic-orca-path-v10-reached-endpoint-minima"


@dataclass(frozen=True)
class EndpointReferenceAssessment:
    """The meaning of an endpoint, independently of its coordinate integrity.

    Fragment minima do not make their finite-separation assembly an asymptote.
    A reference receipt permits comparisons, never a mechanism or rate claim.
    """

    kind: str
    sources: tuple[str, ...] = ()
    reason: str = ""

    def __post_init__(self):
        if self.kind not in {"finite_encounter", "bound_minimum",
                             "qualified_fragment_reference", "unresolved"}:
            raise ValueError("Unknown endpoint reference kind")
        if self.kind in {"bound_minimum", "qualified_fragment_reference"} and not self.sources:
            raise ValueError("Qualified endpoint reference requires native source evidence")

    def as_dict(self):
        return {"kind": self.kind, "sources": list(self.sources), "reason": self.reason,
                "energetic_classification_permitted": self.kind in {
                    "bound_minimum", "qualified_fragment_reference"},
                "mechanism_claim_allowed": False, "rate_claim_allowed": False}

SURFACE_DIAGNOSES = frozenset({
    "single_saddle_candidate",
    "intermediate_candidate",
    "validated_intermediate",
    "intermediate_rejected",
    "endpoint_orientation_artifact",
    "flat_dissociation_candidate",
    "barrierless_candidate",
    "electronic_state_unresolved",
    "surface_resolution_budget_exhausted",
    "surface_unresolved",
})

SURFACE_TERMINAL_STATES = frozenset({
    "irc_path_incomplete",
    "irc_endpoint_identity_unresolved",
    "intermediate_optimization_unconverged",
    "intermediate_nonminimum",
    "intermediate_identity_unresolved",
    "coordinate_scan_no_ts_seed",
    "ts_mode_mismatch",
    "segment_path_unresolved",
    "segment_rate_model_unavailable",
    "verified_multistep_path_network_insertion_pending",
    "barrierless_rate_model_required",
    "dissociation_rate_model_required",
    "electronic_state_escalation_required",
    "environment_model_escalation_required",
    "surface_resolution_budget_exhausted",
    "surface_unresolved",
})


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in sorted(value.items(), key=lambda item: str(item[0]))}
    if isinstance(value, (tuple, list, set, frozenset)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def _scientific_hash(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class SurfaceResolutionPolicy:
    """Immutable bounds and thresholds for one path-resolution attempt."""

    schema_version: int = SURFACE_RESOLUTION_SCHEMA_VERSION
    policy_id: str = SURFACE_RESOLUTION_POLICY_ID
    strategies: tuple[str, ...] = (
        "standard_neb_ts",
        "quantitative_surface_diagnosis",
        "stationary_point_recovery",
        "bond_change_guided_scan",
        "segmented_path_verification",
        "irc_branch_continuation",
        "barrierless_method_escalation",
    )
    maximum_orientations: int = 3
    maximum_stationary_point_recoveries: int = 3
    maximum_scan_points: int = 18
    maximum_scan_coordinates: int = 2
    maximum_intermediates: int = 2
    maximum_segmentation_depth: int = 2
    intermediate_identity_policy_id: str = "rdkit-connectivity-first-charge-spin-map-v3"
    irc_continuation_policy_id: str = "irc-downhill-single-extension-v1"
    maximum_calculations: int = 48
    maximum_core_hours: float = 400.0
    maximum_wall_seconds: float = 7.0 * 24.0 * 3600.0
    meaningful_imaginary_frequency_cm_1: float = -20.0
    energy_extremum_tolerance_hartree: float = 1.0e-4
    endpoint_rmsd_tolerance_angstrom: float = 0.75
    independent_orientation_agreement: int = 2
    allowed_environment_models: tuple[str, ...] = (
        "homogeneous_gas",
        "homogeneous gas-phase low-pressure surrogate",
        "implicit_solution",
    )
    terminal_abstention_states: tuple[str, ...] = tuple(sorted(SURFACE_TERMINAL_STATES))
    endpoint_connectivity_policy_id: str = "distance-gate-exact-closed-shell-ring-diagonal-v1"
    endpoint_geometry_integrity_policy_id: str = "finite-coordinates-and-integer-atom-order-v1"
    declared_stereochemistry_policy_id: str = "declared-graph-rdkit-3d-tetrahedral-ez-v1"
    irc_endpoint_identity_policy_id: str = "mapped-chemistry-conformer-diagnostic-and-reached-minima-v2"

    def __post_init__(self) -> None:
        integer_bounds = (
            self.maximum_orientations,
            self.maximum_stationary_point_recoveries,
            self.maximum_scan_points,
            self.maximum_scan_coordinates,
            self.maximum_intermediates,
            self.maximum_segmentation_depth,
            self.maximum_calculations,
            self.independent_orientation_agreement,
        )
        if any(int(value) < 1 for value in integer_bounds):
            raise ValueError("Surface-resolution integer bounds must be positive")
        if self.maximum_scan_coordinates > 2:
            raise ValueError("The initial surface resolver supports at most two scan coordinates")
        if self.independent_orientation_agreement > self.maximum_orientations:
            raise ValueError("Orientation agreement cannot exceed the orientation budget")
        if not math.isfinite(self.maximum_core_hours) or self.maximum_core_hours <= 0.0:
            raise ValueError("maximum_core_hours must be positive and finite")
        if not math.isfinite(self.maximum_wall_seconds) or self.maximum_wall_seconds <= 0.0:
            raise ValueError("maximum_wall_seconds must be positive and finite")
        if self.meaningful_imaginary_frequency_cm_1 >= 0.0:
            raise ValueError("The meaningful imaginary-frequency threshold must be negative")
        if self.energy_extremum_tolerance_hartree <= 0.0:
            raise ValueError("The energy-extremum tolerance must be positive")
        if not self.strategies or len(set(self.strategies)) != len(self.strategies):
            raise ValueError("Surface-resolution strategies must be nonempty and unique")
        if set(self.terminal_abstention_states) != set(SURFACE_TERMINAL_STATES):
            raise ValueError("Surface-resolution policy must retain every fail-closed terminal state")

    def as_dict(self) -> dict[str, Any]:
        value = _jsonable({
            "kind": "storca_surface_resolution_policy",
            **asdict(self),
        })
        value["scientific_hash"] = _scientific_hash(value)
        return value


@dataclass(frozen=True)
class SurfaceResolutionAttempt:
    """One append-only strategy decision or calculation attempt."""

    attempt_id: str
    strategy: str
    status: str
    input_evidence_hashes: tuple[str, ...] = ()
    scientific_context: Mapping[str, Any] = field(default_factory=dict)
    calculation_parameters: Mapping[str, Any] = field(default_factory=dict)
    budget_consumed: Mapping[str, Any] = field(default_factory=dict)
    process_leases: tuple[Mapping[str, Any], ...] = ()
    artifacts: Mapping[str, Any] = field(default_factory=dict)
    failure_reason: str | None = None
    claims_permitted: tuple[str, ...] = ()
    claims_prohibited: tuple[str, ...] = ("rate", "lifetime")

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SURFACE_RESOLUTION_SCHEMA_VERSION,
            "kind": "storca_surface_resolution_attempt",
            **_jsonable(asdict(self)),
        }


@dataclass(frozen=True)
class SurfaceDiagnosis:
    """A deterministic interpretation of retained path evidence."""

    status: str
    reasons: tuple[str, ...]
    evidence: Mapping[str, Any]
    corroborated: bool = False
    recommended_strategy: str | None = None
    claims_permitted: tuple[str, ...] = ()
    claims_prohibited: tuple[str, ...] = ("mechanism", "rate", "lifetime")

    def __post_init__(self) -> None:
        if self.status not in SURFACE_DIAGNOSES:
            raise ValueError(f"Unsupported surface diagnosis: {self.status}")

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SURFACE_RESOLUTION_SCHEMA_VERSION,
            "kind": "storca_surface_diagnosis",
            **_jsonable(asdict(self)),
        }


@dataclass(frozen=True)
class SurfaceResolutionResult:
    """Terminal or resumable output from the resolution policy."""

    status: str
    diagnosis: SurfaceDiagnosis
    attempts: tuple[SurfaceResolutionAttempt, ...] = ()
    next_strategy: str | None = None
    rate_claim_allowed: bool = False
    lifetime_claim_allowed: bool = False
    evidence: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if (self.rate_claim_allowed or self.lifetime_claim_allowed) and self.status in SURFACE_TERMINAL_STATES:
            raise ValueError("An unresolved surface state cannot permit rate or lifetime claims")

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SURFACE_RESOLUTION_SCHEMA_VERSION,
            "kind": "storca_surface_resolution_result",
            "status": self.status,
            "diagnosis": self.diagnosis.as_dict(),
            "attempts": [item.as_dict() for item in self.attempts],
            "next_strategy": self.next_strategy,
            "rate_claim_allowed": self.rate_claim_allowed,
            "lifetime_claim_allowed": self.lifetime_claim_allowed,
            "evidence": _jsonable(self.evidence),
        }


def _finite_energies(values: Sequence[float] | None) -> list[float] | None:
    if values is None:
        return None
    try:
        energies = [float(value) for value in values]
    except (TypeError, ValueError):
        return None
    if len(energies) < 3 or any(not math.isfinite(value) for value in energies):
        return None
    return energies


def _interior_extrema(values: Sequence[float], tolerance: float) -> tuple[list[int], list[int]]:
    minima: list[int] = []
    maxima: list[int] = []
    for index in range(1, len(values) - 1):
        left = values[index - 1]
        middle = values[index]
        right = values[index + 1]
        if middle <= left - tolerance and middle <= right - tolerance:
            minima.append(index)
        if middle >= left + tolerance and middle >= right + tolerance:
            maxima.append(index)
    return minima, maxima


def _path_geometry_metrics(coordinates: Sequence[Any] | None) -> dict[str, Any]:
    if coordinates is None:
        return {"status": "not_evaluated"}
    try:
        frames = [np.asarray(frame, dtype=float) for frame in coordinates]
    except (TypeError, ValueError):
        return {"status": "invalid_coordinates"}
    if len(frames) < 2 or any(frame.ndim != 2 or frame.shape[1] != 3 for frame in frames):
        return {"status": "invalid_coordinates"}
    if len({frame.shape for frame in frames}) != 1 or any(not np.isfinite(frame).all() for frame in frames):
        return {"status": "invalid_coordinates"}
    step_rmsd = [
        float(np.sqrt(np.mean(np.sum((right - left) ** 2, axis=1))))
        for left, right in zip(frames, frames[1:])
    ]
    return {
        "status": "completed",
        "frame_count": len(frames),
        "path_arc_rms_angstrom": float(sum(step_rmsd)),
        "step_rmsd_angstrom": step_rmsd,
        "maximum_step_rmsd_angstrom": max(step_rmsd),
    }


def _bond_coordinate_metrics(
    coordinates: Sequence[Any] | None,
    bond_changes: Mapping[str, Any] | None,
) -> dict[str, Any]:
    pairs: list[tuple[str, int, int]] = []
    for key, label in (("broken_bonds", "broken"), ("formed_bonds", "formed")):
        for raw in (bond_changes or {}).get(key) or ():
            if len(raw) >= 2:
                pairs.append((label, int(raw[0]), int(raw[1])))
    for raw in (bond_changes or {}).get("bond_order_changes") or ():
        indices = raw.get("atom_indices") if isinstance(raw, Mapping) else None
        if indices and len(indices) >= 2:
            pairs.append(("bond_order_changed", int(indices[0]), int(indices[1])))
    if not pairs or coordinates is None:
        return {"status": "not_evaluated", "coordinates": []}
    try:
        frames = [np.asarray(frame, dtype=float) for frame in coordinates]
        records = []
        for label, left, right in pairs:
            distances = [float(np.linalg.norm(frame[left] - frame[right])) for frame in frames]
            records.append({
                "kind": label,
                "atom_indices": [left, right],
                "distances_angstrom": distances,
                "net_change_angstrom": distances[-1] - distances[0],
                "range_angstrom": max(distances) - min(distances),
            })
        return {"status": "completed", "coordinates": records}
    except (IndexError, TypeError, ValueError):
        return {"status": "invalid_coordinates", "coordinates": []}


def _fragment_motion_metrics(
    coordinates: Sequence[Any] | None,
    fragment_atom_indices: Sequence[Sequence[int]] | None,
) -> dict[str, Any]:
    if coordinates is None or not fragment_atom_indices or len(fragment_atom_indices) < 2:
        return {"status": "not_evaluated"}
    try:
        frames = [np.asarray(frame, dtype=float) for frame in coordinates]
        fragments = [np.asarray(tuple(int(item) for item in fragment), dtype=int) for fragment in fragment_atom_indices]
        if any(len(fragment) == 0 for fragment in fragments):
            raise ValueError("empty fragment")
        first, second = fragments[:2]
        separations = [
            float(np.linalg.norm(frame[first].mean(axis=0) - frame[second].mean(axis=0)))
            for frame in frames
        ]
        rotation = None
        if len(first) >= 2:
            initial_axis = frames[0][first[1]] - frames[0][first[0]]
            final_axis = frames[-1][first[1]] - frames[-1][first[0]]
            denominator = float(np.linalg.norm(initial_axis) * np.linalg.norm(final_axis))
            if denominator > 0.0:
                cosine = float(np.clip(np.dot(initial_axis, final_axis) / denominator, -1.0, 1.0))
                rotation = float(np.degrees(np.arccos(abs(cosine))))
        return {
            "status": "completed",
            "centroid_separation_angstrom": separations,
            "net_separation_change_angstrom": separations[-1] - separations[0],
            "separation_range_angstrom": max(separations) - min(separations),
            "first_fragment_axis_rotation_degrees": rotation,
        }
    except (IndexError, TypeError, ValueError):
        return {"status": "invalid_fragment_partition"}


_COVALENT_RADII_ANGSTROM = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66,
    "F": 0.57, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Br": 1.20, "I": 1.39,
}


def _connectivity_change_metrics(
    atom_labels: Sequence[str] | None,
    coordinates: Sequence[Any] | None,
) -> dict[str, Any]:
    if not atom_labels or coordinates is None:
        return {"status": "not_evaluated"}
    labels = tuple(str(item) for item in atom_labels)
    if any(label not in _COVALENT_RADII_ANGSTROM for label in labels):
        return {"status": "unsupported_element", "atom_labels": list(labels)}
    try:
        frames = [np.asarray(frame, dtype=float) for frame in coordinates]
        if any(frame.shape != (len(labels), 3) or not np.isfinite(frame).all() for frame in frames):
            raise ValueError("invalid frame")
        signatures: list[tuple[tuple[int, int], ...]] = []
        for frame in frames:
            bonds = []
            for left in range(len(labels)):
                for right in range(left + 1, len(labels)):
                    threshold = 1.25 * (
                        _COVALENT_RADII_ANGSTROM[labels[left]]
                        + _COVALENT_RADII_ANGSTROM[labels[right]]
                    )
                    distance = float(np.linalg.norm(frame[left] - frame[right]))
                    if 0.1 < distance <= threshold:
                        bonds.append((left, right))
            signatures.append(tuple(bonds))
        changes = []
        for index, (left, right) in enumerate(zip(signatures, signatures[1:]), start=1):
            broken = sorted(set(left) - set(right))
            formed = sorted(set(right) - set(left))
            if broken or formed:
                changes.append({
                    "between_image_indices": [index - 1, index],
                    "broken_bonds": [list(pair) for pair in broken],
                    "formed_bonds": [list(pair) for pair in formed],
                })
        return {
            "status": "completed",
            "image_count": len(signatures),
            "unique_connectivity_count": len(set(signatures)),
            "connectivity_changed": bool(changes),
            "changes": changes,
            "bond_signatures": [[list(pair) for pair in signature] for signature in signatures],
        }
    except (IndexError, TypeError, ValueError):
        return {"status": "invalid_coordinates"}


def climbing_image_history(output_text: str | None) -> dict[str, Any]:
    """Parse only the image identity retained in ORCA LBFGS records."""
    indices: list[int] = []
    for line in str(output_text or "").splitlines():
        match = re.match(r"^\s*LBFGS\s+\d+\s+(\d+)\b", line, re.IGNORECASE)
        if match:
            indices.append(int(match.group(1)))
    changes = sum(left != right for left, right in zip(indices, indices[1:]))
    return {
        "status": "completed" if indices else "not_available",
        "image_indices": indices,
        "distinct_image_indices": sorted(set(indices)),
        "migration_count": changes,
    }


def diagnose_surface(
    *,
    energies_hartree: Sequence[float] | None,
    coordinates: Sequence[Any] | None = None,
    atom_labels: Sequence[str] | None = None,
    bond_changes: Mapping[str, Any] | None = None,
    fragment_atom_indices: Sequence[Sequence[int]] | None = None,
    explicit_intermediate_indices: Iterable[int] = (),
    output_text: str | None = None,
    evidence_complete: bool = True,
    endpoints_valid: bool = True,
    electronic_state_status: str = "applicable_single_surface",
    reactant_fragment_count: int = 1,
    product_fragment_count: int = 1,
    policy: SurfaceResolutionPolicy | None = None,
) -> SurfaceDiagnosis:
    """Diagnose retained path evidence without promoting it to a claim."""
    policy = policy or SurfaceResolutionPolicy()
    energies = _finite_energies(energies_hartree)
    intermediate_warnings = sorted(set(int(item) for item in explicit_intermediate_indices))
    geometry = _path_geometry_metrics(coordinates)
    bond_metrics = _bond_coordinate_metrics(coordinates, bond_changes)
    fragment_motion = _fragment_motion_metrics(coordinates, fragment_atom_indices)
    climbing = climbing_image_history(output_text)
    connectivity = _connectivity_change_metrics(atom_labels, coordinates)
    evidence = {
        "energies_hartree": energies,
        "evidence_complete": bool(evidence_complete),
        "endpoints_valid": bool(endpoints_valid),
        "electronic_state_status": electronic_state_status,
        "explicit_intermediate_indices": intermediate_warnings,
        "geometry": geometry,
        "bond_coordinates": bond_metrics,
        "fragment_motion": fragment_motion,
        "climbing_image_history": climbing,
        "connectivity_changes": connectivity,
        "reactant_fragment_count": int(reactant_fragment_count),
        "product_fragment_count": int(product_fragment_count),
    }
    if electronic_state_status != "applicable_single_surface":
        return SurfaceDiagnosis(
            status="electronic_state_unresolved",
            reasons=("declared_electronic_surface_is_not_applicable_or_resolved",),
            evidence=evidence,
            recommended_strategy=None,
        )
    if not endpoints_valid or not evidence_complete or energies is None:
        reasons = []
        if not endpoints_valid:
            reasons.append("endpoint_validation_failed")
        if not evidence_complete:
            reasons.append("path_evidence_incomplete")
        if energies is None:
            reasons.append("finite_energy_profile_unavailable")
        return SurfaceDiagnosis(
            status="surface_unresolved",
            reasons=tuple(reasons),
            evidence=evidence,
        )

    minima, maxima = _interior_extrema(energies, policy.energy_extremum_tolerance_hartree)
    differences = np.diff(np.asarray(energies, dtype=float))
    monotonic_up = bool(np.all(differences >= -policy.energy_extremum_tolerance_hartree))
    monotonic_down = bool(np.all(differences <= policy.energy_extremum_tolerance_hartree))
    evidence.update({
        "interior_minima": minima,
        "interior_maxima": maxima,
        "monotonic_up": monotonic_up,
        "monotonic_down": monotonic_down,
        "energy_range_hartree": max(energies) - min(energies),
    })

    warning_corroborated = bool(intermediate_warnings and (
        minima or any(abs(index - candidate) <= 1 for index in minima for candidate in intermediate_warnings)
    ))
    if minima or intermediate_warnings:
        reasons = []
        if minima:
            reasons.append("interior_energy_minimum_detected")
        if intermediate_warnings:
            reasons.append("orca_reported_possible_intermediate")
        if intermediate_warnings and not warning_corroborated:
            reasons.append("orca_warning_not_yet_corroborated_by_energy_extremum")
        return SurfaceDiagnosis(
            status="intermediate_candidate",
            reasons=tuple(reasons),
            evidence=evidence,
            corroborated=warning_corroborated or bool(minima),
            recommended_strategy="stationary_point_recovery",
        )

    separated_channel = reactant_fragment_count != product_fragment_count
    fragment_dominated = (
        fragment_motion.get("status") == "completed"
        and abs(float(fragment_motion.get("net_separation_change_angstrom", 0.0))) >= 1.0
        and not maxima
    )
    if fragment_dominated and not (monotonic_up or monotonic_down):
        return SurfaceDiagnosis(
            status="endpoint_orientation_artifact",
            reasons=("fragment_motion_dominates_without_a_resolved_interior_extremum",),
            evidence=evidence,
            recommended_strategy="bond_change_guided_scan",
        )
    if separated_channel and (monotonic_up or monotonic_down):
        return SurfaceDiagnosis(
            status="flat_dissociation_candidate",
            reasons=("complete_monotonic_profile_changes_fragment_count",),
            evidence=evidence,
            corroborated=False,
            recommended_strategy="barrierless_method_escalation",
        )
    if monotonic_up or monotonic_down:
        return SurfaceDiagnosis(
            status="barrierless_candidate",
            reasons=("complete_monotonic_profile_has_no_interior_extremum",),
            evidence=evidence,
            corroborated=False,
            recommended_strategy="barrierless_method_escalation",
        )
    if len(maxima) == 1:
        return SurfaceDiagnosis(
            status="single_saddle_candidate",
            reasons=("one_interior_energy_maximum_and_no_interior_minimum",),
            evidence=evidence,
            corroborated=True,
            recommended_strategy="bond_change_guided_scan",
        )
    return SurfaceDiagnosis(
        status="surface_unresolved",
        reasons=("energy_profile_does_not_support_a_unique_surface_diagnosis",),
        evidence=evidence,
    )


def aggregate_surface_diagnoses(
    diagnoses: Sequence[SurfaceDiagnosis],
    *,
    policy: SurfaceResolutionPolicy | None = None,
) -> SurfaceDiagnosis:
    """Require independent agreement before a cross-orientation diagnosis."""
    policy = policy or SurfaceResolutionPolicy()
    if not diagnoses:
        return SurfaceDiagnosis(
            status="surface_unresolved",
            reasons=("no_orientation_diagnoses",),
            evidence={"orientation_diagnoses": []},
        )
    counts = {status: sum(item.status == status for item in diagnoses) for status in SURFACE_DIAGNOSES}
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    status, count = ranked[0]
    evidence = {
        "required_agreement": policy.independent_orientation_agreement,
        "status_counts": {key: value for key, value in ranked if value},
        "orientation_diagnoses": [item.as_dict() for item in diagnoses],
    }
    if count < policy.independent_orientation_agreement:
        return SurfaceDiagnosis(
            status="surface_unresolved",
            reasons=("independent_orientations_do_not_agree",),
            evidence=evidence,
        )
    matching = [item for item in diagnoses if item.status == status]
    strategy_counts: dict[str, int] = {}
    for item in matching:
        if item.recommended_strategy:
            strategy_counts[item.recommended_strategy] = strategy_counts.get(item.recommended_strategy, 0) + 1
    strategy = sorted(strategy_counts, key=lambda key: (-strategy_counts[key], key))[0] if strategy_counts else None
    return SurfaceDiagnosis(
        status=status,
        reasons=("independent_orientation_agreement_established",),
        evidence=evidence,
        corroborated=all(item.corroborated for item in matching),
        recommended_strategy=strategy,
    )


def budget_status(
    policy: SurfaceResolutionPolicy,
    *,
    calculation_count: int,
    core_hours: float,
    wall_seconds: float,
) -> dict[str, Any]:
    """Return an auditable budget gate without silently clipping requested work."""
    consumed = {
        "calculation_count": int(calculation_count),
        "core_hours": float(core_hours),
        "wall_seconds": float(wall_seconds),
    }
    exceeded = []
    if calculation_count >= policy.maximum_calculations:
        exceeded.append("maximum_calculations")
    if core_hours >= policy.maximum_core_hours:
        exceeded.append("maximum_core_hours")
    if wall_seconds >= policy.maximum_wall_seconds:
        exceeded.append("maximum_wall_seconds")
    return {
        "schema_version": SURFACE_RESOLUTION_SCHEMA_VERSION,
        "kind": "storca_surface_resolution_budget_status",
        "status": "surface_resolution_budget_exhausted" if exceeded else "within_budget",
        "consumed": consumed,
        "limits": {
            "maximum_calculations": policy.maximum_calculations,
            "maximum_core_hours": policy.maximum_core_hours,
            "maximum_wall_seconds": policy.maximum_wall_seconds,
        },
        "exceeded": exceeded,
        "rate_claim_allowed": False,
        "lifetime_claim_allowed": False,
    }


def build_bond_change_scan_plan(
    reactant_coordinates: Any,
    product_coordinates: Any,
    bond_changes: Mapping[str, Any],
    *,
    policy: SurfaceResolutionPolicy | None = None,
    requested_points: int = 12,
) -> dict[str, Any]:
    """Build a bounded one- or two-coordinate simultaneous relaxed scan.

    Coordinates come only from declared atom-mapped endpoint geometry.  The
    resulting plan is a TS-seed hypothesis, never path or rate evidence.
    """
    policy = policy or SurfaceResolutionPolicy()
    if requested_points < 3:
        raise ValueError("A guided relaxed scan requires at least three points")
    if requested_points > policy.maximum_scan_points:
        raise ValueError("Requested scan points exceed the immutable surface-resolution policy")
    reactant = np.asarray(reactant_coordinates, dtype=float)
    product = np.asarray(product_coordinates, dtype=float)
    if (
        reactant.ndim != 2 or reactant.shape[1] != 3 or product.shape != reactant.shape
        or not np.isfinite(reactant).all() or not np.isfinite(product).all()
    ):
        raise ValueError("Guided scan endpoints need matching finite Cartesian coordinates")

    proposals: list[dict[str, Any]] = []

    def append(kind: str, indices: Sequence[int], priority: int, metadata: Mapping[str, Any] | None = None) -> None:
        if len(indices) < 2:
            return
        left, right = int(indices[0]), int(indices[1])
        if left == right or min(left, right) < 0 or max(left, right) >= len(reactant):
            raise ValueError("Bond-change scan atom index is outside the endpoint geometry")
        start = float(np.linalg.norm(reactant[left] - reactant[right]))
        end = float(np.linalg.norm(product[left] - product[right]))
        if not math.isfinite(start) or not math.isfinite(end) or min(start, end) <= 0.0:
            raise ValueError("Bond-change scan coordinate has an invalid endpoint distance")
        if abs(end - start) < 0.05:
            return
        proposals.append({
            "kind": kind,
            "atom_indices": [left, right],
            "start_distance_angstrom": start,
            "end_distance_angstrom": end,
            "net_change_angstrom": end - start,
            "priority": priority,
            "metadata": _jsonable(metadata or {}),
        })

    for indices in bond_changes.get("broken_bonds") or ():
        append("breaking_bond", indices, 0)
    for item in bond_changes.get("bond_order_changes") or ():
        if isinstance(item, Mapping):
            append("bond_order_change", item.get("atom_indices") or (), 1, item)
    for indices in bond_changes.get("formed_bonds") or ():
        append("forming_bond", indices, 2)

    unique: dict[tuple[int, int], dict[str, Any]] = {}
    for item in sorted(proposals, key=lambda value: (value["priority"], value["atom_indices"])):
        pair = tuple(sorted(item["atom_indices"]))
        unique.setdefault(pair, item)
    selected = list(unique.values())[:policy.maximum_scan_coordinates]
    omitted = list(unique.values())[policy.maximum_scan_coordinates:]
    status = "prepared" if selected else "coordinate_scan_no_ts_seed"
    contract = {
        "schema_version": SURFACE_RESOLUTION_SCHEMA_VERSION,
        "kind": "storca_bond_change_guided_scan_plan",
        "status": status,
        "coordinates": selected,
        "omitted_coordinates": omitted,
        "coordinate_search_truncated": bool(omitted),
        "simultaneous_scan": len(selected) > 1,
        "scan_points": int(requested_points),
        "policy_id": policy.policy_id,
        "claims_permitted": [],
        "claims_prohibited": ["transition_state", "mechanism", "rate", "lifetime"],
    }
    contract["scientific_hash"] = _scientific_hash(contract)
    return contract


def select_interior_scan_maximum(
    energies_hartree: Sequence[float] | None,
    *,
    tolerance_hartree: float = 1.0e-4,
) -> dict[str, Any]:
    """Select a unique interior maximum without treating it as a TS."""
    energies = _finite_energies(energies_hartree)
    if energies is None:
        return {"status": "coordinate_scan_no_ts_seed", "reason": "complete_energy_profile_unavailable"}
    if tolerance_hartree <= 0.0:
        raise ValueError("Scan maximum tolerance must be positive")
    maximum = max(energies)
    candidates = [index for index, value in enumerate(energies) if maximum - value <= tolerance_hartree]
    interior = [index for index in candidates if 0 < index < len(energies) - 1]
    if len(interior) != 1:
        return {
            "status": "coordinate_scan_no_ts_seed",
            "reason": "maximum_is_endpoint_or_not_unique",
            "maximum_indices": candidates,
            "energies_hartree": energies,
        }
    index = interior[0]
    if not (
        energies[index] > energies[index - 1] + tolerance_hartree
        and energies[index] > energies[index + 1] + tolerance_hartree
    ):
        return {
            "status": "coordinate_scan_no_ts_seed",
            "reason": "interior_maximum_is_not_bracketed",
            "maximum_indices": candidates,
            "energies_hartree": energies,
        }
    return {
        "status": "interior_ts_seed_candidate",
        "index": index,
        "energy_hartree": energies[index],
        "left_energy_hartree": energies[index - 1],
        "right_energy_hartree": energies[index + 1],
        "claims_permitted": [],
        "claims_prohibited": ["transition_state", "mechanism", "rate", "lifetime"],
    }


def bond_change_mode_overlap(
    coordinates: Any,
    mode_vector: Any,
    scan_coordinates: Sequence[Mapping[str, Any]],
    *,
    minimum_overlap: float = 0.10,
) -> dict[str, Any]:
    """Check whether a TS mode moves at least one intended bond coordinate."""
    geometry = np.asarray(coordinates, dtype=float)
    vector = np.asarray(mode_vector, dtype=float)
    if geometry.shape != vector.shape or geometry.ndim != 2 or geometry.shape[1] != 3:
        raise ValueError("Mode and geometry must have matching N x 3 coordinates")
    if not np.isfinite(geometry).all() or not np.isfinite(vector).all():
        raise ValueError("Mode-overlap inputs must be finite")
    if not 0.0 < minimum_overlap <= 1.0:
        raise ValueError("minimum_overlap must be in (0, 1]")
    mode_rms = float(np.sqrt(np.mean(np.sum(vector * vector, axis=1))))
    if mode_rms <= 0.0:
        raise ValueError("Transition-state mode has no displacement")
    records = []
    for coordinate in scan_coordinates:
        indices = coordinate.get("atom_indices") or ()
        if len(indices) != 2:
            continue
        left, right = int(indices[0]), int(indices[1])
        bond = geometry[right] - geometry[left]
        distance = float(np.linalg.norm(bond))
        if distance <= 0.0:
            raise ValueError("Mode-overlap coordinate has zero bond distance")
        stretch = float(np.dot(vector[right] - vector[left], bond / distance))
        normalized = abs(stretch) / (math.sqrt(2.0) * mode_rms)
        records.append({
            "atom_indices": [left, right],
            "kind": coordinate.get("kind"),
            "signed_stretch_component": stretch,
            "normalized_overlap": normalized,
        })
    maximum = max((item["normalized_overlap"] for item in records), default=0.0)
    return {
        "schema_version": SURFACE_RESOLUTION_SCHEMA_VERSION,
        "kind": "storca_ts_mode_bond_change_overlap",
        "status": "matched" if maximum >= minimum_overlap else "ts_mode_mismatch",
        "maximum_normalized_overlap": maximum,
        "minimum_required_overlap": minimum_overlap,
        "coordinates": records,
        "rate_claim_allowed": False,
        "lifetime_claim_allowed": False,
    }


__all__ = [
    "SURFACE_DIAGNOSES",
    "SURFACE_RESOLUTION_POLICY_ID",
    "SURFACE_TERMINAL_STATES",
    "SurfaceDiagnosis",
    "SurfaceResolutionAttempt",
    "SurfaceResolutionPolicy",
    "SurfaceResolutionResult",
    "aggregate_surface_diagnoses",
    "bond_change_mode_overlap",
    "build_bond_change_scan_plan",
    "budget_status",
    "climbing_image_history",
    "diagnose_surface",
    "select_interior_scan_maximum",
]
