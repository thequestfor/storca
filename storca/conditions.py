"""Immutable, explicitly scoped conditions for stability screening."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Any


def _freeze(value: Any) -> Any:
    """Recursively make a condition payload safe to retain as a contract."""
    if isinstance(value, dict):
        return tuple(sorted((str(key), _freeze(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, set):
        return tuple(sorted(_freeze(item) for item in value))
    return value


def _thaw(value: Any) -> Any:
    if isinstance(value, tuple):
        if all(isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], str) for item in value):
            return {key: _thaw(item) for key, item in value}
        return [_thaw(item) for item in value]
    return value


def _canonical_smiles(smiles: str) -> str:
    """Canonicalize structure identity for scenario species de-duplication."""
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid scenario species SMILES: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True)


def normalize_target_environment_species(smiles: str, scenario: dict) -> dict:
    """Merge environmental species structurally identical to the target.

    RMG cannot declare two labels for the same molecular graph. When the
    target itself is water, oxygen, or a bath-gas component, its environmental
    fraction is merged into the stability pool. The retained contract records
    this because the model cannot isotope-tag an identical molecular pool.
    """
    normalized = dict(scenario)
    target = _canonical_smiles(smiles)
    composition = dict(normalized.get("initial_mole_fractions") or {})
    additional = []
    merges = []
    for raw in normalized.get("additional_species") or []:
        species = dict(raw)
        canonical = _canonical_smiles(str(species["smiles"]))
        if canonical != target:
            additional.append(species)
            continue
        label = str(species["label"])
        amount = float(composition.pop(label, 0.0))
        composition["stability"] = float(composition.get("stability", 0.0)) + amount
        merges.append({
            "environment_label": label,
            "canonical_smiles": canonical,
            "merged_mole_fraction": amount,
        })
    if not merges:
        return normalized
    normalized["additional_species"] = additional
    normalized["initial_mole_fractions"] = composition
    normalized["target_identity_merges"] = merges
    normalized["model_applicability"] = (
        f"{normalized['model_applicability']} The target is structurally identical to "
        f"{', '.join(item['environment_label'] for item in merges)}; those fractions are merged into "
        "one chemically indistinguishable target pool."
    )
    return normalized


@dataclass(frozen=True)
class ConditionSpec:
    """The physical/model assumptions under which a screen is interpreted."""

    scenario: str
    phase_model: str
    temperature_K: float
    pressure_bar: float
    composition: tuple[tuple[str, float], ...]
    target_label: str
    target_duration_seconds: float
    retention_fraction: float
    light_condition: str
    radical_sources: tuple[str, ...]
    model_applicability: str
    relative_humidity: float | None = None
    light_model: tuple[tuple[str, Any], ...] = ()
    target_identity_merges: tuple[tuple[str, str, float], ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "composition", tuple(sorted((str(key), float(value)) for key, value in self.composition)))
        object.__setattr__(self, "radical_sources", tuple(self.radical_sources))
        object.__setattr__(self, "light_model", _freeze(dict(self.light_model)))
        object.__setattr__(
            self,
            "target_identity_merges",
            tuple((str(label), str(smiles), float(fraction))
                  for label, smiles, fraction in self.target_identity_merges),
        )
        if self.temperature_K <= 0 or self.pressure_bar <= 0:
            raise ValueError("Temperature and pressure must be positive")
        if self.target_duration_seconds <= 0:
            raise ValueError("Target duration must be positive")
        if not 0 < self.retention_fraction < 1:
            raise ValueError("Retention fraction must be between zero and one")
        if self.light_condition not in {"dark", "sunlight"}:
            raise ValueError("Light condition must be 'dark' or 'sunlight'")
        if self.light_condition == "sunlight" and not self.light_model:
            raise ValueError("Sunlight conditions require an immutable light-model contract")
        if self.radical_sources != ("none",):
            raise ValueError("Unsupported conditions: explicit radical-source models are not implemented yet")
        if not self.phase_model.startswith("homogeneous gas-phase"):
            raise ValueError("Unsupported conditions: only homogeneous gas-phase scenarios are currently modelled")
        composition = dict(self.composition)
        if self.target_label not in composition or abs(sum(composition.values()) - 1.0) > 1e-9:
            raise ValueError("Condition composition must contain the target and sum to one")
        if self.relative_humidity is not None and not 0.0 <= self.relative_humidity <= 1.0:
            raise ValueError("Relative humidity must be a fraction from zero through one")

    @property
    def target_mole_fraction(self) -> float:
        return dict(self.composition)[self.target_label]

    def as_dict(self) -> dict:
        payload = {
            "scenario": self.scenario,
            "phase_model": self.phase_model,
            "temperature_K": self.temperature_K,
            "pressure_bar": self.pressure_bar,
            "composition": dict(self.composition),
            "target_label": self.target_label,
            "target_mole_fraction": self.target_mole_fraction,
            "target_duration_seconds": self.target_duration_seconds,
            "retention_fraction": self.retention_fraction,
            "light_condition": self.light_condition,
            "radical_sources": list(self.radical_sources),
            "model_applicability": self.model_applicability,
            "relative_humidity": self.relative_humidity,
            "light_model": _thaw(self.light_model) or None,
            "target_identity_merges": [
                {
                    "environment_label": label,
                    "canonical_smiles": smiles,
                    "merged_mole_fraction": fraction,
                }
                for label, smiles, fraction in self.target_identity_merges
            ],
        }
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
        return {
            "schema_version": 1,
            "contract_id": "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
            **payload,
        }


def build_condition_spec(
    scenario: dict,
    *,
    temperature_K: float,
    pressure_bar: float,
    target_duration_hours: float = 24.0,
    retention_fraction: float = 0.95,
    light_condition: str = "dark",
    radical_sources: tuple[str, ...] = ("none",),
    relative_humidity: float | None = None,
    light_model: dict[str, Any] | None = None,
) -> ConditionSpec:
    """Build a validated condition contract from a named scenario definition."""
    return ConditionSpec(
        scenario=scenario["name"],
        phase_model=scenario["phase"],
        temperature_K=temperature_K,
        pressure_bar=pressure_bar,
        composition=tuple(sorted(scenario["initial_mole_fractions"].items())),
        target_label="stability",
        target_duration_seconds=target_duration_hours * 3600.0,
        retention_fraction=retention_fraction,
        light_condition=light_condition,
        radical_sources=radical_sources,
        model_applicability=scenario["model_applicability"],
        relative_humidity=relative_humidity,
        light_model=_freeze(light_model or {}),
        target_identity_merges=tuple(
            (
                str(item["environment_label"]),
                str(item["canonical_smiles"]),
                float(item["merged_mole_fraction"]),
            )
            for item in scenario.get("target_identity_merges") or []
        ),
    )
