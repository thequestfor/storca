"""Declared electronic levels relative to a native spin-only E0 baseline.

This pure module performs no jobs, source-file ingestion or applicability
inference. Source identifiers and assumptions are declarations, not proof that
the levels are physically appropriate. Degeneracies include spin already;
the native spin multiplicity is divided out, never multiplied in again.
Thermodynamic differences assume temperature-independent declared levels and
a fixed native E0; all non-electronic modal contributions remain unchanged.

For level energies e_i relative to a declared reference whose offset from the
native E0 is d, the new/native partition ratio is
    exp[-(d + min(e_i))/(R T)] * sum_i g_i exp[-(e_i-min(e_i))/(R T)] / spin.
The ground offset changes asymptotic E0/barriers. It must not be used merely
as a population multiplier when replaying energy-dependent tunneling models.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math


CONTRACT_ID = "declared_electronic_levels_v1"
DEGENERACY_CONVENTION = "total_electronic_including_spin"
GAS_CONSTANT_J_MOL_K = 8.31446261815324


def _text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(label + " must be a nonempty explicit string")
    return value


def _finite(value, label, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(label + " must be a finite number")
    try:
        value = float(value)
    except OverflowError as error:
        raise ValueError(label + " must be a finite number") from error
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(label + " must be finite" + (" and positive" if positive else ""))
    return value


def _positive_integer(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(label + " must be a positive integer")
    return value


def _sum(values, label):
    try:
        return _finite(math.fsum(values), label)
    except OverflowError as error:
        raise ValueError(label + " is unrepresentable") from error


def _strings(value, label):
    if not isinstance(value, tuple) or not value:
        raise ValueError(label + " must be a nonempty immutable tuple")
    for item in value:
        _text(item, label)
    if len(set(value)) != len(value):
        raise ValueError(label + " must not contain duplicate entries")


def _keys(record, expected, label):
    if not isinstance(record, dict) or set(record) != set(expected):
        raise ValueError(label + " has missing or unrecognized fields")


@dataclass(frozen=True)
class ElectronicLevel:
    label: str
    energy_J_mol: float
    degeneracy: int
    source_identifier: str

    def __post_init__(self):
        _text(self.label, "level label")
        object.__setattr__(self, "energy_J_mol", _finite(self.energy_J_mol, "level energy_J_mol"))
        _positive_integer(self.degeneracy, "level degeneracy")
        _text(self.source_identifier, "level source_identifier")

    def to_dict(self) -> dict:
        return {"label": self.label, "energy_J_mol": self.energy_J_mol,
                "degeneracy": self.degeneracy, "source_identifier": self.source_identifier}

    @classmethod
    def from_dict(cls, record: dict) -> ElectronicLevel:
        _keys(record, ("label", "energy_J_mol", "degeneracy", "source_identifier"), "electronic level")
        return cls(**record)


@dataclass(frozen=True)
class ElectronicStateModel:
    label: str
    energy_reference: str
    reference_offset_J_mol: float
    native_spin_multiplicity: int
    levels: tuple[ElectronicLevel, ...]
    assumptions: tuple[str, ...]
    source_identifiers: tuple[str, ...]
    degeneracy_convention: str = DEGENERACY_CONVENTION

    def __post_init__(self):
        _text(self.label, "model label")
        _text(self.energy_reference, "energy_reference")
        object.__setattr__(self, "reference_offset_J_mol",
                           _finite(self.reference_offset_J_mol, "reference_offset_J_mol"))
        _positive_integer(self.native_spin_multiplicity, "native_spin_multiplicity")
        if not isinstance(self.levels, tuple) or not self.levels \
                or any(not isinstance(level, ElectronicLevel) for level in self.levels):
            raise ValueError("levels must be a nonempty immutable tuple of ElectronicLevel")
        if len({level.label for level in self.levels}) != len(self.levels):
            raise ValueError("electronic level labels must be unique")
        _strings(self.assumptions, "assumptions")
        _strings(self.source_identifiers, "source_identifiers")
        if any(level.source_identifier not in self.source_identifiers for level in self.levels):
            raise ValueError("every level source_identifier must be declared by the model")
        if self.degeneracy_convention != DEGENERACY_CONVENTION:
            raise ValueError("degeneracies must explicitly include spin; no native spin double count")

    def to_dict(self) -> dict:
        return {"contract_id": CONTRACT_ID, "label": self.label,
                "energy_reference": self.energy_reference,
                "reference_offset_J_mol": self.reference_offset_J_mol,
                "native_spin_multiplicity": self.native_spin_multiplicity,
                "levels": [level.to_dict() for level in self.levels],
                "assumptions": list(self.assumptions), "source_identifiers": list(self.source_identifiers),
                "degeneracy_convention": self.degeneracy_convention}

    @classmethod
    def from_dict(cls, record: dict) -> ElectronicStateModel:
        _keys(record, ("contract_id", "label", "energy_reference", "reference_offset_J_mol",
                       "native_spin_multiplicity", "levels", "assumptions", "source_identifiers",
                       "degeneracy_convention"), "electronic state model")
        if record["contract_id"] != CONTRACT_ID:
            raise ValueError("unrecognized electronic state contract_id")
        for name in ("levels", "assumptions", "source_identifiers"):
            if not isinstance(record[name], list):
                raise ValueError("serialized " + name + " must be a JSON array")
        return cls(record["label"], record["energy_reference"], record["reference_offset_J_mol"],
                   record["native_spin_multiplicity"],
                   tuple(ElectronicLevel.from_dict(level) for level in record["levels"]),
                   tuple(record["assumptions"]), tuple(record["source_identifiers"]),
                   record["degeneracy_convention"])

    @property
    def sha256(self) -> str:
        encoded = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"),
                             allow_nan=False).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()


def _logsumexp(terms):
    maximum = max(terms)
    return _finite(maximum + math.log(math.fsum(math.exp(term - maximum) for term in terms)),
                   "derived log partition")


def _linear(log_value):
    """None explicitly denotes an unrepresentable linear value; log is retained."""
    try:
        value = math.exp(log_value)
    except OverflowError:
        return None
    return value if math.isfinite(value) and value > 0 else None


def _thermal(temperature_K, gas_constant_J_mol_K):
    temperature = _finite(temperature_K, "temperature_K", positive=True)
    gas_constant = _finite(gas_constant_J_mol_K, "gas_constant_J_mol_K", positive=True)
    return temperature, gas_constant, _finite(temperature * gas_constant, "R*T", positive=True)


def evaluate_electronic_partition(
    model: ElectronicStateModel, temperature_K: float, *,
    gas_constant_J_mol_K: float = GAS_CONSTANT_J_MOL_K,
) -> dict:
    """Evaluate declared levels; the correction is Q_new/Q_native, not its inverse.

Q relative to the ground contains total electronic degeneracies. Q relative
to the native baseline divides by the previously counted native spin and
includes the ground energy offset. Linear overflow/underflow is represented
as None with an explicit flag; finite log values remain authoritative.
"""
    if not isinstance(model, ElectronicStateModel):
        raise ValueError("model must be an ElectronicStateModel")
    temperature, gas_constant, thermal_energy = _thermal(temperature_K, gas_constant_J_mol_K)
    ground = min(level.energy_J_mol for level in model.levels)
    offset = _finite(model.reference_offset_J_mol + ground, "derived ground offset")
    terms, gaps = [], []
    for level in model.levels:
        gap = _finite(level.energy_J_mol - ground, "derived excitation energy")
        term = _finite(math.log(level.degeneracy) - gap / thermal_energy, "derived level log weight")
        terms.append(term)
        gaps.append(gap)
    log_ground = _logsumexp(terms)
    log_ground_spin = _finite(log_ground - math.log(model.native_spin_multiplicity),
                              "derived ground partition relative native spin")
    log_native = _finite(log_ground_spin - offset / thermal_energy, "derived native partition correction")
    log_populations = [term - log_ground for term in terms]
    populations = [math.exp(value) for value in log_populations]
    mean = _sum((gap * population for gap, population in zip(gaps, populations)),
                "derived mean excitation energy")
    try:
        variance = _sum(((math.sqrt(population) * (gap - mean)) ** 2
                         for gap, population in zip(gaps, populations)), "derived excitation variance")
    except OverflowError as error:
        raise ValueError("derived excitation variance is unrepresentable") from error
    thermo = {
        "delta_G_J_mol": _finite(offset - thermal_energy * log_ground_spin, "derived delta_G"),
        "delta_H_J_mol": _finite(offset + mean, "derived delta_H"),
        "delta_S_J_mol_K": _finite(gas_constant * log_ground_spin + mean / temperature, "derived delta_S"),
        "delta_Cp_J_mol_K": _finite(variance / thermal_energy / temperature, "derived delta_Cp"),
    }
    linear = {"partition_relative_ground": _linear(log_ground),
              "ground_partition_relative_native_spin": _linear(log_ground_spin),
              "partition_relative_native_baseline": _linear(log_native)}
    return {"contract_id": CONTRACT_ID, "model_sha256": model.sha256,
            "model": model.to_dict(), "temperature_K": temperature,
            "gas_constant_J_mol_K": gas_constant, "energy_unit": "J/mol",
            "ground_level_energy_relative_declared_reference_J_mol": ground,
            "ground_offset_J_mol": offset, **linear,
            "log_partition_relative_ground": log_ground,
            "log_ground_partition_relative_native_spin": log_ground_spin,
            "log_partition_relative_native_baseline": log_native,
            "log_partition_correction": log_native,
            "level_log_populations": dict(zip((level.label for level in model.levels), log_populations)),
            "mean_excitation_energy_J_mol": mean, "excitation_variance_J2_mol2": variance,
            "thermodynamics_relative_native_baseline": thermo,
            "thermodynamic_scope": "temperature_independent_levels_relative_fixed_native_E0",
            "linear_values_representable": all(value is not None for value in linear.values()),
            "source_evidence_verified": False, "physical_applicability_certified": False,
            "network_reinjection_allowed": False, "native_jobs_dispatched": 0}


def evaluate_reaction_electronic_correction(
    reactants: tuple[ElectronicStateModel, ...], products: tuple[ElectronicStateModel, ...],
    transition_state: ElectronicStateModel, temperature_K: float, *,
    gas_constant_J_mol_K: float = GAS_CONSTANT_J_MOL_K,
) -> dict:
    """Classical state corrections and ground-offset barrier shifts for both sides.

Repeat a model in a side tuple to represent integer stoichiometry. This is
partition algebra, not a mass/state/path conservation certificate. No symmetry,
channel count or tunneling correction is inferred. Replay Eckart with *both*
shifted native E0 barriers; multiplying an old Eckart rate is not supported.
"""
    for name, side in (("reactants", reactants), ("products", products)):
        if not isinstance(side, tuple) or not side \
                or any(not isinstance(model, ElectronicStateModel) for model in side):
            raise ValueError(name + " must be a nonempty immutable tuple of ElectronicStateModel")
    temperature, gas_constant, thermal_energy = _thermal(temperature_K, gas_constant_J_mol_K)
    evaluate = lambda model: evaluate_electronic_partition(
        model, temperature, gas_constant_J_mol_K=gas_constant)
    reactant_states, product_states = tuple(map(evaluate, reactants)), tuple(map(evaluate, products))
    ts = evaluate(transition_state)

    def direction(side):
        shift = _sum((ts["ground_offset_J_mol"],
                      *(-item["ground_offset_J_mol"] for item in side)), "derived barrier shift")
        log_population = _sum((ts["log_ground_partition_relative_native_spin"],
                               *(-item["log_ground_partition_relative_native_spin"] for item in side)),
                              "derived log population factor")
        log_factor = _finite(log_population - shift / thermal_energy, "derived log classical rate factor")
        population, factor = _linear(log_population), _linear(log_factor)
        return {"barrier_shift_J_mol": shift, "population_factor": population,
                "log_population_factor": log_population, "classical_rate_factor": factor,
                "log_classical_rate_factor": log_factor,
                "linear_values_representable": population is not None and factor is not None}

    forward, reverse = direction(reactant_states), direction(product_states)
    log_equilibrium = _sum((*(item["log_partition_correction"] for item in product_states),
                            *(-item["log_partition_correction"] for item in reactant_states)),
                           "derived equilibrium correction")
    log_ratio = _finite(forward["log_classical_rate_factor"] - reverse["log_classical_rate_factor"],
                        "derived forward/reverse correction ratio")
    residual = _finite(log_ratio - log_equilibrium, "derived detailed balance residual")
    equilibrium = _linear(log_equilibrium)
    return {"contract_id": CONTRACT_ID, "temperature_K": temperature,
            "gas_constant_J_mol_K": gas_constant, "energy_unit": "J/mol",
            "reactants": list(reactant_states), "products": list(product_states), "transition_state": ts,
            "forward": forward, "reverse": reverse,
            "equilibrium_factor": equilibrium, "log_equilibrium_factor": log_equilibrium,
            "log_forward_reverse_factor_ratio": log_ratio,
            "detailed_balance_log_residual": residual,
            "linear_values_representable": (forward["linear_values_representable"]
                                             and reverse["linear_values_representable"] and equilibrium is not None),
            "tunneling_recalculation_required_if_barriers_shift": bool(
                forward["barrier_shift_J_mol"] or reverse["barrier_shift_J_mol"]),
            "rate_theory_scope": "classical_electronic_partition_correction_only",
            "source_evidence_verified": False, "physical_applicability_certified": False,
            "network_reinjection_allowed": False, "native_jobs_dispatched": 0}
