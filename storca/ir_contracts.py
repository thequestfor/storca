"""Versioned sample, measurement, and reference contracts for FTIR work."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping


IR_CONDITION_SCHEMA_VERSION = 1
IR_REFERENCE_SCHEMA_VERSION = 2
SUPPORTED_PHASES = {"gas", "solution", "liquid", "solid", "unspecified"}
SUPPORTED_MEASUREMENTS = {"atr", "transmission", "gas-cell", "unspecified"}
SUPPORTED_PARTITIONS = {"training", "validation", "holdout", "development_unassigned"}


def _optional_positive(name: str, value: float | None) -> float | None:
    if value is None:
        return None
    resolved = float(value)
    if not math.isfinite(resolved) or resolved <= 0:
        raise ValueError(f"{name} must be finite and greater than zero")
    return resolved


def _optional_text(value: str | None) -> str | None:
    if value is None:
        return None
    resolved = str(value).strip()
    return resolved or None


@dataclass(frozen=True)
class MeasurementCondition:
    geometry: str
    resolution_cm_1: float | None = None
    apodization: str | None = None
    path_length_mm: float | None = None
    atr_crystal: str | None = None
    atr_incidence_angle_degrees: float | None = None
    sample_refractive_index: float | None = None

    def __post_init__(self) -> None:
        geometry = self.geometry.strip().lower()
        if geometry not in SUPPORTED_MEASUREMENTS:
            raise ValueError(f"Unsupported FTIR measurement geometry '{self.geometry}'")
        object.__setattr__(self, "geometry", geometry)
        object.__setattr__(self, "resolution_cm_1", _optional_positive("Resolution", self.resolution_cm_1))
        object.__setattr__(self, "path_length_mm", _optional_positive("Path length", self.path_length_mm))
        object.__setattr__(self, "sample_refractive_index", _optional_positive(
            "Sample refractive index", self.sample_refractive_index,
        ))
        normalized_apodization = _optional_text(self.apodization)
        object.__setattr__(self, "apodization", normalized_apodization.lower() if normalized_apodization else None)
        object.__setattr__(self, "atr_crystal", _optional_text(self.atr_crystal))
        if self.atr_incidence_angle_degrees is not None:
            angle = float(self.atr_incidence_angle_degrees)
            if not math.isfinite(angle) or not 0 < angle < 90:
                raise ValueError("ATR incidence angle must be between zero and 90 degrees")
            object.__setattr__(self, "atr_incidence_angle_degrees", angle)
        atr_fields = (self.atr_crystal, self.atr_incidence_angle_degrees, self.sample_refractive_index)
        if geometry != "atr" and any(value is not None for value in atr_fields):
            raise ValueError("ATR optical properties require an ATR measurement geometry")


@dataclass(frozen=True)
class ExperimentalCondition:
    """Immutable condition contract for one intrinsic or measured IR spectrum."""

    phase: str
    measurement: MeasurementCondition
    temperature_K: float | None = None
    pressure_bar: float | None = None
    composition: Mapping[str, float] = field(default_factory=dict)
    solvent: str | None = None
    concentration_mol_L: float | None = None
    schema_version: int = IR_CONDITION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != IR_CONDITION_SCHEMA_VERSION:
            raise ValueError("Unsupported experimental-condition schema")
        phase = self.phase.strip().lower()
        if phase not in SUPPORTED_PHASES:
            raise ValueError(f"Unsupported FTIR phase '{self.phase}'")
        if phase == "gas" and self.measurement.geometry == "atr":
            raise ValueError("ATR measurement is incompatible with a gas sample")
        if self.measurement.geometry == "gas-cell" and phase not in {"gas", "unspecified"}:
            raise ValueError("Gas-cell measurement requires a gas or unspecified phase")
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "temperature_K", _optional_positive("Temperature", self.temperature_K))
        object.__setattr__(self, "pressure_bar", _optional_positive("Pressure", self.pressure_bar))
        object.__setattr__(self, "concentration_mol_L", _optional_positive(
            "Concentration", self.concentration_mol_L,
        ))
        object.__setattr__(self, "solvent", _optional_text(self.solvent))
        composition: dict[str, float] = {}
        for identity, amount in sorted(dict(self.composition).items()):
            label = str(identity).strip()
            value = float(amount)
            if not label or not math.isfinite(value) or value < 0:
                raise ValueError("Composition requires named, finite, nonnegative fractions")
            composition[label] = value
        if composition and sum(composition.values()) <= 0:
            raise ValueError("Composition must contain a positive amount")
        object.__setattr__(self, "composition", composition)

    def as_dict(self) -> dict:
        return asdict(self)


def build_experimental_condition(
    *, phase: str, measurement: str, temperature_K: float | None = None,
    pressure_bar: float | None = None, composition: Mapping[str, float] | None = None,
    solvent: str | None = None, concentration_mol_L: float | None = None,
    resolution_cm_1: float | None = None, apodization: str | None = None,
    path_length_mm: float | None = None, atr_crystal: str | None = None,
    atr_incidence_angle_degrees: float | None = None,
    sample_refractive_index: float | None = None,
) -> ExperimentalCondition:
    return ExperimentalCondition(
        phase=phase,
        temperature_K=temperature_K,
        pressure_bar=pressure_bar,
        composition=composition or {},
        solvent=solvent,
        concentration_mol_L=concentration_mol_L,
        measurement=MeasurementCondition(
            geometry=measurement,
            resolution_cm_1=resolution_cm_1,
            apodization=apodization,
            path_length_mm=path_length_mm,
            atr_crystal=atr_crystal,
            atr_incidence_angle_degrees=atr_incidence_angle_degrees,
            sample_refractive_index=sample_refractive_index,
        ),
    )


def write_experimental_condition(path: Path, condition: ExperimentalCondition) -> Path:
    path = Path(path)
    path.write_text(json.dumps(condition.as_dict(), indent=2, sort_keys=True) + "\n")
    return path


def condition_compatibility(predicted: Mapping, reference: Mapping) -> dict:
    """Return a fail-closed condition comparison without inventing unknown data."""
    left, right = dict(predicted), dict(reference)
    mismatches: list[str] = []
    unknown: list[str] = []
    for key in ("phase", "temperature_K", "pressure_bar"):
        left_value, right_value = left.get(key), right.get(key)
        if left_value in {None, "", "unspecified"} or right_value in {None, "", "unspecified"}:
            unknown.append(key)
            continue
        if key == "temperature_K":
            if abs(float(left_value) - float(right_value)) > 2.0:
                mismatches.append(key)
        elif key == "pressure_bar":
            if not math.isclose(float(left_value), float(right_value), rel_tol=0.05, abs_tol=0.05):
                mismatches.append(key)
        elif str(left_value).lower() != str(right_value).lower():
            mismatches.append(key)
    phases = {str(left.get("phase", "")).lower(), str(right.get("phase", "")).lower()}
    if phases & {"solution", "liquid", "solid"}:
        left_solvent, right_solvent = left.get("solvent"), right.get("solvent")
        if left_solvent in {None, "", "unspecified"} or right_solvent in {None, "", "unspecified"}:
            unknown.append("solvent")
        elif str(left_solvent).lower() != str(right_solvent).lower():
            mismatches.append("solvent")
    if phases == {"solution"}:
        left_concentration, right_concentration = left.get("concentration_mol_L"), right.get("concentration_mol_L")
        if left_concentration is None or right_concentration is None:
            unknown.append("concentration_mol_L")
        elif not math.isclose(float(left_concentration), float(right_concentration), rel_tol=0.02):
            mismatches.append("concentration_mol_L")
    left_measurement = left.get("measurement") or {}
    right_measurement = right.get("measurement") or {}
    if not isinstance(left_measurement, Mapping):
        left_measurement = {"geometry": left_measurement}
    if not isinstance(right_measurement, Mapping):
        right_measurement = {"geometry": right_measurement}
    left_geometry = left_measurement.get("geometry")
    right_geometry = right_measurement.get("geometry")
    if left_geometry in {None, "", "unspecified"} or right_geometry in {None, "", "unspecified"}:
        unknown.append("measurement.geometry")
    elif str(left_geometry).lower() != str(right_geometry).lower():
        mismatches.append("measurement.geometry")
    left_composition, right_composition = left.get("composition"), right.get("composition")
    if not isinstance(left_composition, Mapping) or not left_composition or not isinstance(right_composition, Mapping) or not right_composition:
        unknown.append("composition")
    else:
        left_total, right_total = sum(map(float, left_composition.values())), sum(map(float, right_composition.values()))
        left_normalized = {str(key): float(value) / left_total for key, value in left_composition.items()}
        right_normalized = {str(key): float(value) / right_total for key, value in right_composition.items()}
        if set(left_normalized) != set(right_normalized) or any(
            not math.isclose(left_normalized[key], right_normalized[key], abs_tol=0.02)
            for key in set(left_normalized) & set(right_normalized)
        ):
            mismatches.append("composition")
    for key, tolerance in (("resolution_cm_1", 0.1),):
        left_value, right_value = left_measurement.get(key), right_measurement.get(key)
        if left_value is None or right_value is None:
            unknown.append(f"measurement.{key}")
        elif not math.isclose(float(left_value), float(right_value), rel_tol=tolerance, abs_tol=tolerance):
            mismatches.append(f"measurement.{key}")
    for key in ("apodization",):
        left_value, right_value = left_measurement.get(key), right_measurement.get(key)
        if left_value in {None, "", "unspecified"} or right_value in {None, "", "unspecified"}:
            unknown.append(f"measurement.{key}")
        elif str(left_value).lower() != str(right_value).lower():
            mismatches.append(f"measurement.{key}")
    if str(left_geometry).lower() == str(right_geometry).lower() == "atr":
        for key in ("atr_crystal", "atr_incidence_angle_degrees", "sample_refractive_index"):
            left_value, right_value = left_measurement.get(key), right_measurement.get(key)
            if left_value is None or right_value is None:
                unknown.append(f"measurement.{key}")
            elif isinstance(left_value, (int, float)) and isinstance(right_value, (int, float)):
                if not math.isclose(float(left_value), float(right_value), rel_tol=0.02, abs_tol=0.1):
                    mismatches.append(f"measurement.{key}")
            elif str(left_value).lower() != str(right_value).lower():
                mismatches.append(f"measurement.{key}")
    if str(left_geometry).lower() == str(right_geometry).lower() in {"transmission", "gas-cell"}:
        left_path, right_path = left_measurement.get("path_length_mm"), right_measurement.get("path_length_mm")
        if left_path is None or right_path is None:
            unknown.append("measurement.path_length_mm")
        elif not math.isclose(float(left_path), float(right_path), rel_tol=0.02, abs_tol=1e-6):
            mismatches.append("measurement.path_length_mm")
    status = "incompatible" if mismatches else "insufficient_metadata" if unknown else "compatible"
    return {
        "status": status,
        "quantitative_comparison_allowed": status == "compatible",
        "mismatches": mismatches,
        "unknown_or_unspecified": unknown,
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_reference_entry(entry: Mapping, root: Path) -> dict:
    """Validate and normalize one schema-v2 benchmark entry."""
    normalized = dict(entry)
    if not normalized.get("id") or not normalized.get("smiles"):
        raise ValueError("Reference entries require id and smiles")
    partition = normalized.get("partition", "development_unassigned")
    if partition not in SUPPORTED_PARTITIONS:
        raise ValueError(f"Unsupported benchmark partition '{partition}'")
    normalized["partition"] = partition
    reference = Path(root) / str(normalized.get("reference_spectrum", ""))
    if not reference.is_file():
        raise ValueError(f"Reference spectrum is missing for {normalized['id']}")
    actual_hash = sha256_file(reference)
    declared_hash = normalized.get("reference_sha256")
    if declared_hash and declared_hash != actual_hash:
        raise ValueError(f"Reference content hash changed for {normalized['id']}")
    normalized["reference_sha256"] = actual_hash
    if not normalized.get("provenance"):
        raise ValueError(f"Reference provenance is required for {normalized['id']}")
    if not isinstance(normalized.get("condition"), Mapping):
        raise ValueError(f"Reference condition contract is required for {normalized['id']}")
    return normalized


def assert_partition_separation(entries: list[Mapping]) -> None:
    ids: set[str] = set()
    identities: dict[str, str] = {}
    for entry in entries:
        entry_id = str(entry["id"])
        if entry_id in ids:
            raise ValueError("IR benchmark manifest IDs must be present and unique")
        ids.add(entry_id)
        partition = str(entry.get("partition", "development_unassigned"))
        identity = str(entry.get("canonical_identity") or entry.get("smiles"))
        previous = identities.get(identity)
        if previous and previous != partition and "development_unassigned" not in {previous, partition}:
            raise ValueError(f"Chemical identity {identity} occurs in multiple benchmark partitions")
        identities[identity] = partition
