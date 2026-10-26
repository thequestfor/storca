"""Generate, execute, and retain Arkane pressure-dependent route jobs."""

from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from typing import Any, Mapping

from src.orca_runner import run_arkane, run_bounded_kinetics_command
from .arkane_rotor import POLICY_ID as ARKANE_STATMECH_POLICY, create_rotor_contract, validate_rotor_receipt


@dataclass(frozen=True)
class VerifiedArkaneRouteSpec:
    """Family-neutral Arkane contract for an ORCA-verified elementary step.

    Every entry points at a completed ORCA frequency calculation.  Reactant
    and product cardinality determine whether Arkane may construct a
    pressure-dependent unimolecular network; bimolecular and exchange
    reactions receive high-pressure-limit TST kinetics only.
    """

    label: str
    reactant_labels: tuple[str, ...]
    reactant_smiles: tuple[str, ...]
    reactant_orca_outputs: tuple[Path, ...]
    reactant_multiplicities: tuple[int, ...]
    product_labels: tuple[str, ...]
    product_smiles: tuple[str, ...]
    product_orca_outputs: tuple[Path, ...]
    product_multiplicities: tuple[int, ...]
    transition_state_orca_output: Path
    transition_state_multiplicity: int
    model_chemistry: str
    reaction_equation: str | None = None
    temperatures_K: tuple[float, ...] = (250.0, 273.15, 298.15, 323.15, 400.0, 500.0, 750.0, 1000.0)
    pressures_bar: tuple[float, ...] = (1.0e-6, 1.0e-4, 1.0e-2, 1.0, 10.0)
    bath_gas: tuple[tuple[str, float], ...] = (("nitrogen", 1.0),)
    pressure_dependent: bool | None = None
    method: str = "modified strong collision"
    interpolation_model: tuple[object, ...] = ("Chebyshev", 6, 4)
    tunneling: str = "Eckart"
    frequency_scale_factor: float = 1.0
    maximum_grain_size_kcal_mol: float = 0.5
    minimum_grain_count: int = 250
    condition_temperature_K: float = 298.15
    condition_pressure_bar: float = 1.0
    collision_limit_m3_mol_s: float | None = None
    collision_limit_source: str | None = None
    explicit_pressure_dependence_model: Mapping[str, Any] | None = None
    conformer_rate_basis: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        if self.conformer_rate_basis is not None:
            basis = self.conformer_rate_basis
            if (not isinstance(basis, Mapping)
                    or basis.get("kind") != "reached_endpoint_conformer_rate_basis"
                    or not isinstance(basis.get("network_reinjection_allowed"), bool)):
                raise ValueError("Invalid reached-endpoint conformer rate basis")
        reactant_lengths = {
            len(self.reactant_labels), len(self.reactant_smiles),
            len(self.reactant_orca_outputs), len(self.reactant_multiplicities),
        }
        product_lengths = {
            len(self.product_labels), len(self.product_smiles),
            len(self.product_orca_outputs), len(self.product_multiplicities),
        }
        if reactant_lengths != {len(self.reactant_labels)} or not self.reactant_labels:
            raise ValueError("Every reactant needs a label, SMILES, ORCA output, and multiplicity")
        if product_lengths != {len(self.product_labels)} or not self.product_labels:
            raise ValueError("Every product needs a label, SMILES, ORCA output, and multiplicity")
        if any(value < 1 for value in (*self.reactant_multiplicities, *self.product_multiplicities,
                                       self.transition_state_multiplicity)):
            raise ValueError("Arkane spin multiplicities must be positive")
        # Retained/imported paths still have to obey the geometry engine's
        # adiabatic spin contract. A positive TS multiplicity alone does not
        # establish that either endpoint can exist on the selected surface.
        from .route_geometry import _coupled_multiplicities

        common_surfaces = set(_coupled_multiplicities(self.reactant_multiplicities)) & set(
            _coupled_multiplicities(self.product_multiplicities)
        )
        if self.transition_state_multiplicity not in common_surfaces:
            raise ValueError(
                "Arkane transition-state multiplicity must belong to a common endpoint spin surface; "
                f"selected={self.transition_state_multiplicity}, common={sorted(common_surfaces)}"
            )
        if any(not str(label).strip() for label in (*self.reactant_labels, *self.product_labels)):
            raise ValueError("Arkane species labels cannot be empty")
        if not self.model_chemistry.strip():
            raise ValueError("A validated Arkane model-chemistry string is required")
        if self.reaction_equation is not None and not re.search(r"<=>|=>|=", self.reaction_equation):
            raise ValueError("Retained Arkane replacement equation lacks a reaction arrow")
        if len(self.temperatures_K) < 3 or any(not math.isfinite(value) or value <= 0 for value in self.temperatures_K):
            raise ValueError("Arkane three-parameter TST fitting needs at least three positive finite temperatures")
        if tuple(sorted(self.temperatures_K)) != self.temperatures_K or len(set(self.temperatures_K)) != len(self.temperatures_K):
            raise ValueError("Arkane temperatures must be unique and strictly increasing")
        if not self.pressures_bar or any(not math.isfinite(value) or value <= 0 for value in self.pressures_bar):
            raise ValueError("Arkane pressures must be positive and finite")
        if tuple(sorted(self.pressures_bar)) != self.pressures_bar or len(set(self.pressures_bar)) != len(self.pressures_bar):
            raise ValueError("Arkane pressures must be unique and strictly increasing")
        if self.use_pressure_dependence:
            if not self.bath_gas or len({name for name, _ in self.bath_gas}) != len(self.bath_gas):
                raise ValueError("Arkane bath gas must contain unique component labels")
            if any(not math.isfinite(value) or value <= 0 for _, value in self.bath_gas):
                raise ValueError("Arkane bath-gas mole fractions must be positive and finite")
            if abs(sum(value for _, value in self.bath_gas) - 1.0) > 1e-9:
                raise ValueError("Arkane bath-gas mole fractions must sum to one")
        if not math.isfinite(self.frequency_scale_factor) or self.frequency_scale_factor <= 0:
            raise ValueError("Arkane frequency scale factor must be positive and finite")
        if not math.isfinite(self.condition_temperature_K) or self.condition_temperature_K <= 0:
            raise ValueError("Arkane condition temperature must be positive and finite")
        if not math.isfinite(self.condition_pressure_bar) or self.condition_pressure_bar <= 0:
            raise ValueError("Arkane condition pressure must be positive and finite")
        if not min(self.temperatures_K) <= self.condition_temperature_K <= max(self.temperatures_K):
            raise ValueError("Condition temperature lies outside the Arkane kinetics grid")
        if len(self.reactant_labels) > 2 or len(self.product_labels) > 2:
            raise ValueError("Arkane verification currently supports only one- and two-species endpoint channels")
        if len(self.reactant_labels) == 2:
            if self.collision_limit_m3_mol_s is None or not math.isfinite(self.collision_limit_m3_mol_s) \
                    or self.collision_limit_m3_mol_s <= 0:
                raise ValueError("Bimolecular Arkane verification requires a positive finite collision limit")
            if not str(self.collision_limit_source or "").strip():
                raise ValueError("Bimolecular Arkane verification requires collision-limit provenance")
        if self.minimum_grain_count < 1 or not math.isfinite(self.maximum_grain_size_kcal_mol) \
                or self.maximum_grain_size_kcal_mol <= 0:
            raise ValueError("Arkane pressure dependence needs a positive grain size and grain count")
        if str(self.tunneling).casefold() not in {"", "none", "wigner", "eckart"}:
            raise ValueError("Arkane tunneling must be none, Wigner, or Eckart")

        model = tuple(self.interpolation_model)
        if not model:
            raise ValueError("Arkane interpolation model cannot be empty")
        if str(model[0]).casefold() == "chebyshev":
            if len(model) != 3 or not all(isinstance(value, int) and value > 0 for value in model[1:]):
                raise ValueError("Chebyshev interpolation requires positive temperature and pressure degrees")
            if model[1] > len(self.temperatures_K) or model[2] > len(self.pressures_bar):
                raise ValueError("Chebyshev polynomial degree cannot exceed its Arkane sampling grid")
        elif str(model[0]).casefold() != "pdeparrhenius" or len(model) != 1:
            raise ValueError("Arkane interpolation must be ('Chebyshev', Tdegree, Pdegree) or ('PDepArrhenius',)")

    @property
    def use_pressure_dependence(self) -> bool:
        # A master equation is needed only when this elementary route connects
        # a bimolecular channel to a single energized well.  Association is
        # oriented around the product well but its reported rate remains in the
        # original 2 -> 1 direction.  Ordinary 1 <-> 1 isomerization and 2 <-> 2
        # exchange receive high-pressure-limit TST.
        requires_pdep = {len(self.reactant_labels), len(self.product_labels)} == {1, 2}
        requested = requires_pdep if self.pressure_dependent is None else self.pressure_dependent
        if bool(requested) != requires_pdep:
            raise ValueError(
                "Verified Arkane branching requires pressure dependence for 1 <-> 2 well/channel routes "
                "and high-pressure-limit TST for 1 <-> 1 or 2 <-> 2 routes"
            )
        if requested:
            if len(self.pressures_bar) < 2:
                raise ValueError("Arkane pressure dependence requires at least two pressures")
            if not min(self.pressures_bar) <= self.condition_pressure_bar <= max(self.pressures_bar):
                raise ValueError("Condition pressure lies outside the Arkane pressure grid")
        return bool(requested)


@dataclass(frozen=True)
class ArkaneRouteSpec:
    """Verified elementary route inputs required by Arkane.

    The ORCA files must be completed endpoint/TS calculations.  This class does
    not infer a TS, spin state, or collisional model from a structural alert.
    """
    label: str
    reactant_label: str
    reactant_smiles: str
    reactant_orca_output: Path
    product_labels: tuple[str, ...]
    product_smiles: tuple[str, ...]
    product_orca_outputs: tuple[Path, ...]
    transition_state_orca_output: Path
    bath_gas: tuple[tuple[str, float], ...] = (("nitrogen", 1.0),)
    temperatures_K: tuple[float, ...] = (298.15,)
    pressures_bar: tuple[float, ...] = (1e-6, 1.0)
    method: str = "modified strong collision"
    interpolation_model: tuple[object, ...] = ("Chebyshev", 6, 4)
    refinement_profile: str | None = None

    def __post_init__(self) -> None:
        if len(self.product_labels) != len(self.product_smiles) or len(self.product_labels) != len(self.product_orca_outputs):
            raise ValueError("Every product needs a label, SMILES, and ORCA output")
        if not self.temperatures_K or not self.pressures_bar:
            raise ValueError("Arkane route needs at least one temperature and pressure")
        if abs(sum(value for _, value in self.bath_gas) - 1.0) > 1e-9:
            raise ValueError("Arkane bath-gas mole fractions must sum to one")


def _py(value: object) -> str:
    return repr(value)


def _safe_label(label: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "-", label).strip("-.")
    return value or "species"


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _finite_positive(value: object, field: str, *, allow_zero: bool = False) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be numeric")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{field} must be numeric") from error
    valid = math.isfinite(number) and (number >= 0.0 if allow_zero else number > 0.0)
    if not valid:
        raise ValueError(f"{field} must be {'nonnegative' if allow_zero else 'positive'} and finite")
    return number


def _validated_explicit_pressure_model(spec: VerifiedArkaneRouteSpec) -> dict | None:
    """Validate an optional route-specific collisional model without inference."""
    raw = spec.explicit_pressure_dependence_model
    if raw is None:
        return None
    if not isinstance(raw, Mapping):
        raise TypeError("Explicit pressure-dependence model must be a mapping")
    model = dict(raw)
    if model.get("schema_version") != 1 \
            or model.get("kind") != "storca_explicit_gas_pressure_dependence_model":
        raise ValueError("Unsupported explicit pressure-dependence model schema")
    if model.get("status") != "applicable":
        raise ValueError("Explicit pressure-dependence model applicability is not established")
    temperature_range = model.get("temperature_range_K") or ()
    pressure_range = model.get("pressure_range_bar") or ()
    if len(temperature_range) != 2 or not (
        float(temperature_range[0]) <= spec.condition_temperature_K <= float(temperature_range[1])
    ):
        raise ValueError("Explicit pressure-dependence model does not cover the condition temperature")
    if len(pressure_range) != 2 or not (
        float(pressure_range[0]) <= spec.condition_pressure_bar <= float(pressure_range[1])
    ):
        raise ValueError("Explicit pressure-dependence model does not cover the condition pressure")
    well_label = (
        spec.reactant_labels[0] if len(spec.reactant_labels) == 1
        else spec.product_labels[0]
    )
    energy = dict(model.get("energy_transfer_model") or {})
    if energy.get("model") != "SingleExponentialDown" \
            or str(energy.get("well_label")) != well_label:
        raise ValueError("Explicit energy-transfer model must be SingleExponentialDown for the exact well")
    normalized_energy = {
        "model": "SingleExponentialDown",
        "well_label": well_label,
        "alpha0_kJ_mol": _finite_positive(energy.get("alpha0_kJ_mol"), "alpha0_kJ_mol"),
        "T0_K": _finite_positive(energy.get("T0_K"), "energy-transfer T0_K"),
        "n": _finite_positive(energy.get("n"), "energy-transfer exponent", allow_zero=True),
        "source": str(energy.get("source") or "").strip(),
        "applicability": str(energy.get("applicability") or "").strip(),
    }
    if not normalized_energy["source"] or not normalized_energy["applicability"]:
        raise ValueError("Explicit energy-transfer model needs source and applicability evidence")

    transport = model.get("transport_models") or {}
    if not isinstance(transport, Mapping):
        raise ValueError("Explicit transport_models must map species labels to models")
    required_labels = {well_label, *(name for name, _ in spec.bath_gas)}
    normalized_transport = {}
    for label in sorted(required_labels):
        item = dict(transport.get(label) or {})
        if item.get("model") != "LennardJones":
            raise ValueError(f"Explicit Lennard-Jones transport model missing for {label}")
        normalized_transport[label] = {
            "model": "LennardJones",
            "sigma_angstrom": _finite_positive(item.get("sigma_angstrom"), f"{label} sigma"),
            "epsilon_K": _finite_positive(item.get("epsilon_K"), f"{label} epsilon"),
            "shape_index": int(item.get("shape_index")),
            "source": str(item.get("source") or "").strip(),
            "applicability": str(item.get("applicability") or "").strip(),
        }
        if normalized_transport[label]["shape_index"] not in {0, 1, 2}:
            raise ValueError(f"Explicit transport shape index is invalid for {label}")
        if not normalized_transport[label]["source"] \
                or not normalized_transport[label]["applicability"]:
            raise ValueError(f"Explicit transport model needs source and applicability for {label}")
    artifact_checks = []
    artifacts = model.get("evidence_artifacts") or ()
    if not artifacts:
        raise ValueError("Explicit pressure-dependence model needs hashed evidence artifacts")
    for raw_artifact in artifacts:
        artifact = Path(str((raw_artifact or {}).get("path") or ""))
        expected = str((raw_artifact or {}).get("sha256") or "")
        actual = _sha256(artifact) if artifact.is_file() else None
        if actual is None or actual != expected:
            raise ValueError("Explicit pressure-dependence evidence artifact is missing or mutated")
        artifact_checks.append({"path": str(artifact.resolve()), "sha256": actual})
    normalized = {
        "schema_version": 1,
        "kind": "storca_explicit_gas_pressure_dependence_model",
        "status": "applicable",
        "temperature_range_K": [float(value) for value in temperature_range],
        "pressure_range_bar": [float(value) for value in pressure_range],
        "energy_transfer_model": normalized_energy,
        "transport_models": normalized_transport,
        "evidence_artifacts": artifact_checks,
    }
    encoded = json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
    return {**normalized, "model_sha256": "sha256:" + hashlib.sha256(encoded).hexdigest()}


def _smiles_identity(smiles: str) -> tuple[dict[str, int], int]:
    """Return the explicit-H elemental formula and formal charge for a SMILES."""
    try:
        from rdkit import Chem
    except ImportError as error:
        raise RuntimeError("RDKit is required to validate Arkane endpoint composition") from error
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Could not parse Arkane endpoint SMILES: {smiles}")
    molecule = Chem.AddHs(molecule)
    formula = Counter(atom.GetSymbol() for atom in molecule.GetAtoms())
    charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    return dict(sorted(formula.items())), int(charge)


def _canonical_smiles(smiles: str) -> str:
    try:
        from rdkit import Chem
    except ImportError as error:
        raise RuntimeError("RDKit is required to validate Arkane endpoint structures") from error
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Could not parse Arkane SMILES: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True)


def _orca_identity(text: str) -> tuple[dict[str, int], int, int]:
    charge_matches = re.findall(r"^\s*Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)\s*$", text, re.MULTILINE)
    multiplicity_matches = re.findall(r"^\s*Multiplicity\s+Mult\s+\.{2,}\s+(\d+)\s*$", text, re.MULTILINE)
    if not charge_matches or not multiplicity_matches:
        raise ValueError("ORCA output does not expose its total charge and multiplicity")
    lines = text.splitlines()
    starts = [index for index, line in enumerate(lines) if "CARTESIAN COORDINATES (ANGSTROEM)" in line]
    if not starts:
        raise ValueError("ORCA output does not contain a Cartesian geometry")
    atoms: list[str] = []
    started = False
    for line in lines[starts[-1] + 1:]:
        stripped = line.strip()
        if not stripped or set(stripped) == {"-"}:
            if started:
                break
            continue
        fields = stripped.split()
        if len(fields) >= 4 and re.fullmatch(r"[A-Z][a-z]?", fields[0]):
            try:
                tuple(float(value.replace("D", "E")) for value in fields[1:4])
            except ValueError:
                if started:
                    break
                continue
            atoms.append(fields[0])
            started = True
        elif started:
            break
    if not atoms:
        raise ValueError("ORCA output has an empty or unparseable final Cartesian geometry")
    return dict(sorted(Counter(atoms).items())), int(charge_matches[-1]), int(multiplicity_matches[-1])


def _combined_identity(smiles_values: tuple[str, ...]) -> tuple[dict[str, int], int]:
    formula: Counter[str] = Counter()
    charge = 0
    for smiles in smiles_values:
        item_formula, item_charge = _smiles_identity(smiles)
        formula.update(item_formula)
        charge += item_charge
    return dict(sorted(formula.items())), charge


def _validate_orca_frequency_artifact(
    path: Path,
    *,
    transition_state: bool,
    expected_formula: dict[str, int],
    expected_charge: int,
    expected_multiplicity: int,
) -> dict:
    """Validate an ORCA statmech artifact before exposing it to Arkane."""
    from src.stability.freq_check import frequency_stability_check

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Required ORCA output not found: {path}")
    text = path.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in text:
        raise ValueError(f"ORCA output did not terminate normally: {path}")
    formula, charge, multiplicity = _orca_identity(text)
    if formula != expected_formula:
        raise ValueError(f"ORCA artifact composition {formula} does not match expected {expected_formula}: {path}")
    if charge != expected_charge:
        raise ValueError(f"ORCA artifact charge {charge} does not match expected {expected_charge}: {path}")
    if multiplicity != expected_multiplicity:
        raise ValueError(
            f"ORCA artifact multiplicity {multiplicity} does not match expected {expected_multiplicity}: {path}"
        )
    hessian = path.with_suffix(".hess")
    if not hessian.is_file():
        raise ValueError(f"Arkane requires the ORCA Hessian beside its frequency output: {hessian}")
    check = frequency_stability_check(path)
    expected = 1 if transition_state else 0
    if check["NumImaginary"] != expected:
        kind = "transition state" if transition_state else "endpoint"
        raise ValueError(
            f"Arkane {kind} artifact has {check['NumImaginary']} significant imaginary modes; expected {expected}: {path}"
        )
    if transition_state:
        negative_modes = [value for value in check["AllFrequencies"] if value < 0.0]
        significant_mode = check["ImaginaryFrequencies"][0]
        # Arkane's ORCA adapter selects the first printed negative mode.  Make
        # sure that is the validated reaction coordinate, even if harmless
        # near-zero translational modes are also printed later.
        if not negative_modes or not math.isclose(negative_modes[0], significant_mode, abs_tol=1e-9):
            raise ValueError(
                "Arkane would not select the validated significant imaginary frequency from this ORCA TS output"
            )
    return {
        "output": str(path.resolve()), "hessian": str(hessian.resolve()), "frequency_check": check,
        "formula": formula, "charge": charge, "multiplicity": multiplicity,
        "checksums": {"output": _sha256(path), "hessian": _sha256(hessian)},
    }


def create_verified_arkane_input(spec: VerifiedArkaneRouteSpec, workdir: Path) -> Path:
    """Create a runnable Arkane job from validated ORCA frequency outputs.

    Arkane's ``species`` and ``transitionState`` positional arguments refer to
    small statmech Python files, not directly to electronic-structure logs.
    Each generated file uses Arkane's ORCA ``Log`` adapter for energy,
    geometry, frequencies, and the colocated Hessian. Electronic energies
    receive a balanced per-atom reference shift before statmech processing.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    # Trigger branching validation even when a caller only creates the input.
    use_pdep = spec.use_pressure_dependence
    bath_gas = spec.bath_gas if use_pdep else ()
    explicit_pressure_model = (
        _validated_explicit_pressure_model(spec) if use_pdep else None
    )
    reactant_formula, reactant_charge = _combined_identity(spec.reactant_smiles)
    product_formula, product_charge = _combined_identity(spec.product_smiles)
    if reactant_formula != product_formula or reactant_charge != product_charge:
        raise ValueError(
            "Arkane reaction endpoints are not atom/charge balanced: "
            f"reactants={reactant_formula}, charge={reactant_charge}; "
            f"products={product_formula}, charge={product_charge}"
        )

    records = []
    species_dir = workdir / "species"
    species_dir.mkdir(exist_ok=True)

    label_contracts: dict[str, tuple[str, Path, int]] = {}
    for label, smiles, output, multiplicity in zip(
        (*spec.reactant_labels, *spec.product_labels),
        (*spec.reactant_smiles, *spec.product_smiles),
        (*spec.reactant_orca_outputs, *spec.product_orca_outputs),
        (*spec.reactant_multiplicities, *spec.product_multiplicities),
    ):
        contract = (smiles, Path(output).resolve(), multiplicity)
        if label in label_contracts and label_contracts[label] != contract:
            raise ValueError(f"Repeated Arkane label {label!r} has inconsistent structure, output, or spin")
        label_contracts[label] = contract
    used_files: set[str] = set()

    def statmech_file(
        label: str,
        output: Path,
        multiplicity: int,
        *,
        formula: dict[str, int],
        charge: int,
        is_ts: bool = False,
    ) -> Path:
        validation = _validate_orca_frequency_artifact(
            output, transition_state=is_ts, expected_formula=formula,
            expected_charge=charge, expected_multiplicity=multiplicity,
        )
        stem = _safe_label(label)
        candidate = stem
        serial = 2
        while candidate in used_files:
            candidate = f"{stem}-{serial}"
            serial += 1
        used_files.add(candidate)
        artifact_dir = species_dir / f"{candidate}-artifacts"
        artifact_dir.mkdir(exist_ok=True)
        staged_output = artifact_dir / "calculation.out"
        staged_hessian = artifact_dir / "calculation.hess"
        shutil.copy2(validation["output"], staged_output)
        shutil.copy2(validation["hessian"], staged_hessian)
        file = species_dir / f"{candidate}.py"
        resolved = str(staged_output.resolve())
        file.write_text(
            f"energy = Log({_py(resolved)})\n"
            f"geometry = Log({_py(resolved)})\n"
            f"frequencies = Log({_py(resolved)})\n"
            f"spin_multiplicity = {multiplicity}\n"
            "opticalIsomers = 1\n"
        )
        records.append({"label": label, "kind": "transition_state" if is_ts else "species",
                        "statmech_file": str(file.resolve()), "staged_output": str(staged_output.resolve()),
                        "staged_hessian": str(staged_hessian.resolve()), **validation})
        return file

    species_files: dict[str, Path] = {}
    for label, (smiles, output, multiplicity) in label_contracts.items():
        formula, charge = _smiles_identity(smiles)
        species_files[label] = statmech_file(
            label, output, multiplicity, formula=formula, charge=charge,
        )
    ts_label = "__storca_ts__"
    while ts_label in label_contracts or any(ts_label == name for name, _ in bath_gas):
        ts_label += "_"
    ts_file = statmech_file(
        ts_label, spec.transition_state_orca_output, spec.transition_state_multiplicity,
        formula=reactant_formula, charge=reactant_charge, is_ts=True,
    )

    # RMG exponentiates absolute energies in several pressure-network routines.
    # Atom corrections are disabled, so raw electronic energies can overflow.
    # A common per-atom gauge preserves every balanced channel energy and TS
    # barrier, including 1 <-> 2 channels where a per-species shift would not.
    # Keep the electronic-energy tag: a two-item tuple would incorrectly tell
    # Arkane that ZPE was already included.
    atom_count = sum(reactant_formula.values())
    by_label = {record["label"]: record for record in records}
    reference_terms = [
        f"Log({_py(by_label[label]['staged_output'])}).load_energy()"
        for label in spec.reactant_labels
    ]
    reference_expression = "(" + " + ".join(reference_terms) + f") / {atom_count}"
    for record in records:
        count = sum(record["formula"].values())
        file = Path(record["statmech_file"])
        original = file.read_text().splitlines(keepends=True)
        original[0] = (
            f"energy = (Log({_py(record['staged_output'])}).load_energy() - "
            f"({reference_expression}) * {count}, 'J/mol', 'e_electronic')\n"
        )
        file.write_text("".join(original))
        record["energy_reference"] = {
            "policy_id": "balanced-per-atom-electronic-gauge-v1",
            "reference_reactant_labels": list(spec.reactant_labels),
            "reference_atom_count": atom_count,
            "species_atom_count": count,
            "zero_point_energy_included": False,
            "absolute_formation_thermochemistry_claim_allowed": False,
        }

    bath_smiles = {
        "nitrogen": "N#N", "N2": "N#N",
        "oxygen": "[O][O]", "O2": "[O][O]",
        "water": "O", "H2O": "O",
        "argon": "[Ar]", "Ar": "[Ar]",
        "helium": "[He]", "He": "[He]",
    }
    unknown = sorted(name for name, _ in bath_gas if name not in bath_smiles)
    if unknown:
        raise ValueError(f"Unsupported Arkane bath-gas labels: {', '.join(unknown)}")
    for name, _ in bath_gas:
        if name in label_contracts and _canonical_smiles(label_contracts[name][0]) != _canonical_smiles(bath_smiles[name]):
            raise ValueError(f"Bath-gas label {name!r} collides with a different reacting-species structure")
    lines = [
        "# Generated by STORCA from validated ORCA frequency/Hessian artifacts.",
        f"modelChemistry = {_py(spec.model_chemistry)}",
        f"frequencyScaleFactor = {spec.frequency_scale_factor!r}",
        "useHinderedRotors = False",
        "useAtomCorrections = False",
        "useBondCorrections = False",
    ]
    for label, (smiles, _, _) in label_contracts.items():
        file = species_files[label]
        collision_arguments = ""
        if explicit_pressure_model is not None \
                and label in explicit_pressure_model["transport_models"]:
            transport = explicit_pressure_model["transport_models"][label]
            collision_arguments += (
                f", collisionModel=TransportData(shapeIndex={transport['shape_index']!r}, "
                f"epsilon=({transport['epsilon_K']!r}, 'K'), "
                f"sigma=({transport['sigma_angstrom']!r}, 'angstrom'))"
            )
        if explicit_pressure_model is not None \
                and label == explicit_pressure_model["energy_transfer_model"]["well_label"]:
            energy = explicit_pressure_model["energy_transfer_model"]
            collision_arguments += (
                f", energyTransferModel=SingleExponentialDown("
                f"alpha0=({energy['alpha0_kJ_mol']!r}, 'kJ/mol'), "
                f"T0=({energy['T0_K']!r}, 'K'), n={energy['n']!r})"
            )
        lines.append(
            f"species({_py(label)}, {_py(str(file.relative_to(workdir)))}, "
            f"structure=SMILES({_py(smiles)}){collision_arguments})"
        )
    for gas, _ in bath_gas:
        if gas not in label_contracts:
            collision_arguments = ""
            if explicit_pressure_model is not None:
                transport = explicit_pressure_model["transport_models"][gas]
                collision_arguments = (
                    f", collisionModel=TransportData(shapeIndex={transport['shape_index']!r}, "
                    f"epsilon=({transport['epsilon_K']!r}, 'K'), "
                    f"sigma=({transport['sigma_angstrom']!r}, 'angstrom'))"
                )
            lines.append(
                f"species({_py(gas)}, structure=SMILES({_py(bath_smiles[gas])}), "
                f"reactive=False{collision_arguments})"
            )
    lines.extend([
        f"transitionState({_py(ts_label)}, {_py(str(ts_file.relative_to(workdir)))})",
        f"reaction({_py(spec.label)}, reactants={list(spec.reactant_labels)!r}, "
        f"products={list(spec.product_labels)!r}, transitionState={_py(ts_label)}, tunneling={_py(spec.tunneling)})",
        f"kinetics({_py(spec.label)}, Tlist=({list(spec.temperatures_K)!r}, 'K'))",
    ])
    if use_pdep:
        if len(spec.reactant_labels) == 1:
            isomers = list(spec.reactant_labels)
            reactant_channels: list[list[str]] = []
            product_channels = [list(spec.product_labels)]
        else:
            isomers = list(spec.product_labels)
            reactant_channels = [list(spec.reactant_labels)]
            product_channels = []
        lines.extend([
            f"network({_py(spec.label + '_network')}, isomers={isomers!r}, "
            f"reactants={reactant_channels!r}, products={product_channels!r}, pathReactions={[spec.label]!r}, "
            f"bathGas={dict(spec.bath_gas)!r})",
            f"pressureDependence({_py(spec.label + '_network')}, "
            f"Tlist=({list(spec.temperatures_K)!r}, 'K'), "
            f"Plist=({list(spec.pressures_bar)!r}, 'bar'), method={_py(spec.method)}, "
            f"interpolationModel={spec.interpolation_model!r}, "
            f"maximumGrainSize=({spec.maximum_grain_size_kcal_mol!r}, 'kcal/mol'), "
            f"minimumGrainCount={spec.minimum_grain_count})",
        ])
    path = workdir / "input.py"
    path.write_text("\n".join(lines) + "\n")
    (workdir / "provenance.json").write_text(json.dumps({
        "schema_version": 2,
        "kind": "verified_arkane_route",
        "arkane_statmech_policy": ARKANE_STATMECH_POLICY,
        "route": {key: value for key, value in spec.__dict__.items()},
        "pressure_dependent": use_pdep,
        "bath_gas_applicable": use_pdep,
        "explicit_pressure_dependence_model": explicit_pressure_model,
        "reaction_identity": {"formula": reactant_formula, "charge": reactant_charge},
        "orca_artifacts": records,
        "assumptions": [
            "Rigid-rotor/harmonic-oscillator statmech; hindered rotors are not inferred.",
            "Atom corrections cancel for the atom-balanced route; bond corrections are disabled rather than inferred.",
            "A common per-atom electronic-energy reference preserves balanced reaction energies and barriers; absolute formation thermochemistry is not supplied. ZPE is still added by Arkane.",
            (
                "No bath-gas or collisional model is used for high-pressure-limit TST."
                if not use_pdep else
                "Route-specific hashed transport and single-exponential-down evidence is used."
                if explicit_pressure_model is not None else
                "Arkane/RMG estimates gas transport and uses its generic single-exponential-down energy-transfer model for pressure dependence."
            ),
        ],
    }, indent=2, sort_keys=True, default=str) + "\n")
    return path


def create_arkane_input(spec: ArkaneRouteSpec, workdir: Path) -> Path:
    """Write an Arkane pdep job only for a validated refinement profile.

    Raw ORCA ``.out`` files are not Arkane species files.  The old scaffold
    generated inputs Arkane could not interpret, so it now fails before
    presenting uncalibrated screen outputs as kinetics.
    """
    if not spec.refinement_profile:
        raise ValueError(
            "Arkane refinement is not configured: provide a validated profile with generated Arkane species/TS files; "
            "screen-level ORCA outputs cannot be used directly."
        )
    if len(spec.temperatures_K) < 2 or len(spec.pressures_bar) < 2:
        raise ValueError("Arkane pressure-dependent kinetics requires at least two temperatures and two pressures")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    for path in (spec.reactant_orca_output, *spec.product_orca_outputs, spec.transition_state_orca_output):
        if not Path(path).is_file():
            raise FileNotFoundError(f"Required ORCA output not found: {path}")
    bath_smiles = {"nitrogen": "N#N", "oxygen": "[O][O]", "argon": "[Ar]"}
    unknown = [name for name, _ in spec.bath_gas if name not in bath_smiles]
    if unknown:
        raise ValueError(f"Unsupported bath-gas labels: {', '.join(unknown)}")
    lines = [
        "# Generated by STORCA. All energies/frequencies come from retained ORCA files.",
        f"modelChemistry = 'ORCA route artifacts; see provenance.json'",
        f"species({_py(spec.reactant_label)}, {_py(str(Path(spec.reactant_orca_output).resolve()))})",
    ]
    for label, output in zip(spec.product_labels, spec.product_orca_outputs):
        lines.append(f"species({_py(label)}, {_py(str(Path(output).resolve()))})")
    for gas, _ in spec.bath_gas:
        lines.append(f"species({_py(gas)}, structure=SMILES({_py(bath_smiles[gas])}), reactive=False)")
    lines.extend([
        f"transitionState('TS', {_py(str(Path(spec.transition_state_orca_output).resolve()))})",
        f"reaction({_py(spec.label)}, reactants=[{_py(spec.reactant_label)}], products={list(spec.product_labels)!r}, transitionState='TS')",
        f"network({_py(spec.label + '_network')}, isomers=[{_py(spec.reactant_label)}], products={[list(spec.product_labels)]!r}, pathReactions=[{_py(spec.label)}], bathGas={dict(spec.bath_gas)!r})",
        f"pressureDependence({_py(spec.label + '_network')}, Tlist={list(spec.temperatures_K)!r}, Plist={[ (pressure, 'bar') for pressure in spec.pressures_bar ]!r}, method={_py(spec.method)}, interpolationModel={spec.interpolation_model!r})",
    ])
    path = workdir / "input.py"
    path.write_text("\n".join(lines) + "\n")
    (workdir / "provenance.json").write_text(json.dumps({
        "kind": "arkane_pressure_dependent_route",
        "route": {key: (str(value) if isinstance(value, Path) else value) for key, value in spec.__dict__.items()},
        "interpretation": "Requires a verified TS/IRC and consistent ORCA endpoint artifacts.",
    }, indent=2, default=str, sort_keys=True) + "\n")
    return path


def parse_arkane_pdep_output(path: Path) -> dict:
    """Return the retained Arkane output and whether it contains pdep kinetics."""
    text = Path(path).read_text()
    reactions = re.findall(r"pdepreaction\(", text)
    return {"path": str(path), "pressure_dependent_reaction_count": len(reactions),
            "status": "completed" if reactions else "completed_without_pdep_reaction"}


def _rmg_python_command(rmg_env: str | None, script: str, *arguments: str) -> list[str]:
    return (["conda", "run", "-n", rmg_env, "python", "-c", script, *arguments]
            if rmg_env else ["python", "-c", script, *arguments])


def _marked_json(stdout: str, marker: str) -> dict:
    matches = [line[len(marker):] for line in stdout.splitlines() if line.startswith(marker)]
    if not matches:
        raise RuntimeError(f"RMG helper did not emit {marker.rstrip('=')}")
    return json.loads(matches[-1])


def validate_arkane_input_syntax(input_file: Path, *, rmg_env: str | None = "rmg_env") -> dict:
    """Ask the selected Arkane installation to parse the job without executing chemistry."""
    marker = "STORCA_ARKANE_INPUT_JSON="
    script = r'''
import json, sys
import rmgpy
from arkane.input import load_input_file
jobs, reactions, species, transition_states, networks, level = load_input_file(sys.argv[1])
print("STORCA_ARKANE_INPUT_JSON=" + json.dumps({
    "rmg_version": getattr(rmgpy, "__version__", None),
    "job_types": [job.__class__.__name__ for job in jobs],
    "reaction_labels": sorted(reactions),
    "species_labels": sorted(species),
    "transition_state_labels": sorted(transition_states),
    "network_labels": sorted(networks),
    "model_chemistry": level.to_model_chem() if level is not None else None,
}, sort_keys=True))
'''
    command = _rmg_python_command(rmg_env, script, str(Path(input_file).resolve()))
    result = run_bounded_kinetics_command(command, cwd=Path(input_file).parent)
    if result.returncode != 0:
        raise RuntimeError(
            f"Arkane input failed syntax/API validation in {rmg_env or 'the active environment'}: "
            f"{result.stderr.strip() or result.stdout.strip()}"
        )
    return {**_marked_json(result.stdout, marker), "command": command}


def _literal_keyword(call: ast.Call, name: str) -> object:
    for keyword in call.keywords:
        if keyword.arg == name:
            return ast.literal_eval(keyword.value)
    raise ValueError(f"Arkane output call lacks {name!r}")


def _kinetics_expression(call: ast.Call, source: str) -> tuple[str, str]:
    for keyword in call.keywords:
        if keyword.arg == "kinetics":
            expression = ast.get_source_segment(source, keyword.value)
            if not expression or not isinstance(keyword.value, ast.Call) or not isinstance(keyword.value.func, ast.Name):
                raise ValueError("Arkane output kinetics expression is not a direct RMG kinetics constructor")
            return expression, keyword.value.func.id
    raise ValueError("Arkane output call lacks a kinetics model")


def _extract_arkane_rate_model(spec: VerifiedArkaneRouteSpec, output_file: Path) -> dict:
    """Extract exactly the intended Arkane-generated model from output.py."""
    source = Path(output_file).read_text(errors="replace")
    tree = ast.parse(source, filename=str(output_file))
    matches: list[dict] = []
    if spec.use_pressure_dependence:
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name) \
                    or node.func.id != "pdepreaction":
                continue
            reactants = tuple(_literal_keyword(node, "reactants"))
            products = tuple(_literal_keyword(node, "products"))
            if Counter(reactants) == Counter(spec.reactant_labels) and Counter(products) == Counter(spec.product_labels):
                expression, model_class = _kinetics_expression(node, source)
                matches.append({"expression": expression, "model_class": model_class,
                                "reactants": list(reactants), "products": list(products)})
        allowed = {"Chebyshev", "PDepArrhenius"}
    else:
        # The TST table is only written when Arkane calculated the rate from
        # the validated TS rather than accepting pre-existing kinetics.
        if "k (TST)" not in source or "k (TST+T)" not in source:
            raise ValueError("Arkane output lacks evidence that high-pressure TST was executed")
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name) or node.func.id != "kinetics":
                continue
            if _literal_keyword(node, "label") == spec.label:
                expression, model_class = _kinetics_expression(node, source)
                matches.append({"expression": expression, "model_class": model_class,
                                "reactants": list(spec.reactant_labels), "products": list(spec.product_labels)})
        allowed = {"Arrhenius", "MultiArrhenius"}
    if len(matches) != 1:
        raise ValueError(f"Arkane output contains {len(matches)} matching condition-rate models; expected exactly one")
    if matches[0]["model_class"] not in allowed:
        raise ValueError(
            f"Arkane returned {matches[0]['model_class']} for the requested "
            f"{'pressure-dependent' if spec.use_pressure_dependence else 'high-pressure TST'} branch"
        )
    return matches[0]


def _evaluate_kinetics_expression(
    expression: str,
    *,
    temperature_K: float,
    pressure_bar: float,
    pressure_dependent: bool,
    rmg_env: str | None,
) -> dict:
    """Evaluate an Arkane/RMG model at the immutable condition in its owning environment."""
    marker = "STORCA_ARKANE_RATE_JSON="
    script = r'''
import json, math, sys
import rmgpy
import rmgpy.kinetics as kinetics_module
context = {name: getattr(kinetics_module, name) for name in dir(kinetics_module) if not name.startswith("_")}
model = eval(sys.argv[1], {"__builtins__": {}}, context)
temperature = float(sys.argv[2])
pressure_pa = float(sys.argv[3]) * 1.0e5
pressure_dependent = sys.argv[4] == "1"
value = model.get_rate_coefficient(temperature, pressure_pa) if pressure_dependent else model.get_rate_coefficient(temperature)
if not math.isfinite(value) or value <= 0.0:
    raise ValueError("Arkane model returned a non-finite or non-positive condition rate")
print("STORCA_ARKANE_RATE_JSON=" + json.dumps({
    "value_si": value, "model_class": model.__class__.__name__,
    "canonical_expression": repr(model), "rmg_version": getattr(rmgpy, "__version__", None),
}, sort_keys=True))
'''
    command = _rmg_python_command(
        rmg_env, script, expression, repr(temperature_K), repr(pressure_bar), "1" if pressure_dependent else "0",
    )
    result = run_bounded_kinetics_command(command)
    if result.returncode != 0:
        raise RuntimeError(
            "Could not evaluate the Arkane condition-rate model: "
            f"{result.stderr.strip() or result.stdout.strip()}"
        )
    return {**_marked_json(result.stdout, marker), "command": command[:-4] + ["<kinetics>", *command[-3:]]}


def _tst_output_block(spec: VerifiedArkaneRouteSpec, output_file: Path) -> str:
    """Select comments belonging to this kinetics call, not another reaction."""
    source = Path(output_file).read_text(errors="replace")
    calls = sorted((node for node in ast.walk(ast.parse(source))
                    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "kinetics"), key=lambda node: node.lineno)
    matching = [index for index, call in enumerate(calls)
                if _literal_keyword(call, "label") == spec.label]
    if len(matching) != 1:
        raise ValueError("Retained Arkane output does not identify exactly one intended TST reaction")
    index = matching[0]
    start = calls[index - 1].end_lineno if index else 0
    return "\n".join(source.splitlines()[start:calls[index].lineno - 1])


def _validated_fit_expression(expression: str, *, molecularity: int | None = None) -> str:
    """Permit only serialized Arrhenius constructors and literal parameters."""
    tree = ast.parse(expression, mode="eval")
    if not isinstance(tree.body, ast.Call) or not isinstance(tree.body.func, ast.Name) \
            or tree.body.func.id not in {"Arrhenius", "MultiArrhenius"}:
        raise ValueError("Retained TST fit is not an allowed Arrhenius constructor")
    allowed = (ast.Expression, ast.Call, ast.Name, ast.Load, ast.keyword, ast.Constant,
               ast.Tuple, ast.List, ast.Dict, ast.UnaryOp, ast.USub, ast.UAdd)
    for node in ast.walk(tree):
        if not isinstance(node, allowed) or isinstance(node, ast.Name) \
                and node.id not in {"Arrhenius", "MultiArrhenius"}:
            raise ValueError("Retained TST fit contains a non-literal or unsupported expression")
        if isinstance(node, ast.Call) and (not isinstance(node.func, ast.Name)
                                          or node.func.id not in {"Arrhenius", "MultiArrhenius"}):
            raise ValueError("Retained TST fit contains an unsupported nested constructor")
        if isinstance(node, ast.Call) and node.func.id == "Arrhenius" and molecularity is not None:
            preexponential = _literal_keyword(node, "A")
            allowed_units = {1: {"s^-1"}, 2: {"cm^3/(mol*s)", "m^3/(mol*s)"}}[molecularity]
            if not isinstance(preexponential, (tuple, list)) or len(preexponential) < 2 \
                    or preexponential[1] not in allowed_units:
                raise ValueError("Retained TST fit rate units differ from the reaction molecularity")
    return expression


def _retained_arkane_artifact_basis(spec: VerifiedArkaneRouteSpec, workdir: Path) -> dict:
    """Bind a read-only diagnostic to the actual inputs, minima and TS artifacts."""
    input_file, provenance_file = workdir / "input.py", workdir / "provenance.json"
    provenance = json.loads(provenance_file.read_text())
    if provenance.get("kind") != "verified_arkane_route":
        raise ValueError("Retained Arkane provenance is not a verified elementary-route contract")
    retained = provenance.get("route") or {}
    current = json.loads(json.dumps(spec.__dict__, default=str))
    required = {"label", "reactant_labels", "reactant_smiles", "reactant_orca_outputs",
                "reactant_multiplicities", "product_labels", "product_smiles", "product_orca_outputs",
                "product_multiplicities", "transition_state_orca_output", "transition_state_multiplicity",
                "model_chemistry", "frequency_scale_factor", "tunneling", "condition_temperature_K",
                "condition_pressure_bar", "temperatures_K"}
    if not required <= retained.keys() or any(retained[key] != current[key] for key in retained):
        raise ValueError("Retained Arkane route differs from the requested condition/artifact contract")
    if provenance.get("pressure_dependent") is not spec.use_pressure_dependence:
        raise ValueError("Retained Arkane provenance has a different pressure model")
    tree = ast.parse(input_file.read_text())
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name)]
    reactions = [node for node in calls if node.func.id == "reaction" and node.args
                 and ast.literal_eval(node.args[0]) == spec.label]
    if len(reactions) != 1 or tuple(_literal_keyword(reactions[0], "reactants")) != spec.reactant_labels \
            or tuple(_literal_keyword(reactions[0], "products")) != spec.product_labels \
            or _literal_keyword(reactions[0], "tunneling") != spec.tunneling:
        raise ValueError("Retained Arkane input changes reaction identity or tunneling model")
    ts_label = _literal_keyword(reactions[0], "transitionState")
    assignments = {node.targets[0].id: ast.literal_eval(node.value) for node in tree.body
                   if isinstance(node, ast.Assign) and len(node.targets) == 1
                   and isinstance(node.targets[0], ast.Name)}
    if assignments.get("modelChemistry") != spec.model_chemistry \
            or assignments.get("frequencyScaleFactor") != spec.frequency_scale_factor \
            or assignments.get("useHinderedRotors") is not False \
            or assignments.get("useAtomCorrections") is not False \
            or assignments.get("useBondCorrections") is not False:
        raise ValueError("Retained Arkane input differs from the electronic/statistical-mechanics model")
    kinetics_calls = [node for node in calls if node.func.id == "kinetics" and node.args
                      and ast.literal_eval(node.args[0]) == spec.label]
    if len(kinetics_calls) != 1 or tuple(_literal_keyword(kinetics_calls[0], "Tlist")) \
            != (list(spec.temperatures_K), "K"):
        raise ValueError("Retained Arkane input changes the sampled temperature grid")
    expected = [(label, Path(output), multiplicity, "species") for label, output, multiplicity in zip(
        (*spec.reactant_labels, *spec.product_labels),
        (*spec.reactant_orca_outputs, *spec.product_orca_outputs),
        (*spec.reactant_multiplicities, *spec.product_multiplicities))]
    expected.append((ts_label, Path(spec.transition_state_orca_output),
                     spec.transition_state_multiplicity, "transition_state"))
    basis: dict[str, dict] = {}
    for label, output, multiplicity, kind in expected:
        records = [record for record in provenance.get("orca_artifacts", ())
                   if record.get("label") == label and record.get("kind") == kind]
        if len(records) != 1:
            raise ValueError(f"Retained Arkane artifact basis is ambiguous or missing for {label}")
        record = records[0]
        if Path(record["output"]).resolve() != output.resolve() or record.get("multiplicity") != multiplicity:
            raise ValueError(f"Retained Arkane endpoint/state basis changed for {label}")
        frequencies = record.get("frequency_check") or {}
        if kind == "species" and (frequencies.get("IsMinimum") is not True
                                  or frequencies.get("NumImaginary") != 0) \
                or kind == "transition_state" and frequencies.get("NumImaginary") != 1:
            raise ValueError(f"Retained Arkane stationary-point frequency evidence is incomplete for {label}")
        hashes = record.get("checksums") or {}
        artifacts = {}
        for field, checksum in (("output", "output"), ("hessian", "hessian"),
                                ("staged_output", "output"), ("staged_hessian", "hessian")):
            path = Path(record[field])
            digest = _sha256(path)
            if not hashes.get(checksum) or digest != hashes[checksum]:
                raise ValueError(f"Retained Arkane artifact hash changed: {label}/{field}")
            artifacts[field] = {"path": str(path.resolve()), "sha256": digest}
        statmech_file = Path(record["statmech_file"])
        declarations = [call for call in calls if call.func.id == ("transitionState" if kind == "transition_state"
                                                                  else "species") and call.args
                        and ast.literal_eval(call.args[0]) == label]
        if len(declarations) != 1 or len(declarations[0].args) < 2 \
                or (workdir / ast.literal_eval(declarations[0].args[1])).resolve() != statmech_file.resolve():
            raise ValueError(f"Retained Arkane input does not use the claimed statmech file for {label}")
        if kind == "species":
            smiles_call = next((keyword.value for keyword in declarations[0].keywords
                                if keyword.arg == "structure"), None)
            declared_smiles = dict(zip((*spec.reactant_labels, *spec.product_labels),
                                      (*spec.reactant_smiles, *spec.product_smiles)))[label]
            if not isinstance(smiles_call, ast.Call) or not isinstance(smiles_call.func, ast.Name) \
                    or smiles_call.func.id != "SMILES" or len(smiles_call.args) != 1 \
                    or ast.literal_eval(smiles_call.args[0]) != declared_smiles:
                raise ValueError(f"Retained Arkane input changes the declared graph for {label}")
        statmech = ast.parse(statmech_file.read_text())
        terms = {node.targets[0].id: node.value for node in statmech.body
                 if isinstance(node, ast.Assign) and len(node.targets) == 1
                 and isinstance(node.targets[0], ast.Name)}
        for field in ("geometry", "frequencies"):
            value = terms.get(field)
            if not isinstance(value, ast.Call) or not isinstance(value.func, ast.Name) \
                    or value.func.id != "Log" or len(value.args) != 1 \
                    or Path(ast.literal_eval(value.args[0])).resolve() != Path(record["staged_output"]).resolve():
                raise ValueError(f"Retained Arkane {field} does not use the claimed reached endpoint for {label}")
        if ast.literal_eval(terms.get("spin_multiplicity")) != multiplicity:
            raise ValueError(f"Retained Arkane statmech multiplicity changed for {label}")
        reference = record.get("energy_reference") or {}
        source_by_label = {item["label"]: item for item in provenance["orca_artifacts"]}
        atom_count = sum(sum(source_by_label[item]["formula"].values()) for item in spec.reactant_labels)
        species_atom_count = sum(record["formula"].values())
        expected_reference = {"policy_id": "balanced-per-atom-electronic-gauge-v1",
                              "reference_reactant_labels": list(spec.reactant_labels),
                              "reference_atom_count": atom_count, "species_atom_count": species_atom_count,
                              "zero_point_energy_included": False,
                              "absolute_formation_thermochemistry_claim_allowed": False}
        if reference != expected_reference:
            raise ValueError(f"Retained Arkane balanced electronic-energy reference changed for {label}")
        reference_terms = " + ".join(f"Log({source_by_label[item]['staged_output']!r}).load_energy()"
                                     for item in spec.reactant_labels)
        expected_energy = ast.parse(
            f"(Log({record['staged_output']!r}).load_energy() - "
            f"(({reference_terms}) / {atom_count}) * {species_atom_count}, 'J/mol', 'e_electronic')",
            mode="eval").body
        if ast.dump(terms.get("energy")) != ast.dump(expected_energy):
            raise ValueError(f"Retained Arkane energy does not use the claimed balanced endpoint basis for {label}")
        artifacts["statmech_file"] = {"path": str(statmech_file.resolve()), "sha256": _sha256(statmech_file)}
        basis[label] = {"label": label, "kind": kind, "multiplicity": multiplicity, "artifacts": artifacts}
    return {"reactants": [basis[label] for label in spec.reactant_labels],
            "products": [basis[label] for label in spec.product_labels],
            "transition_state": basis[ts_label],
            "product_basis_used_by_tunneling": str(spec.tunneling).casefold() == "eckart",
            "all_conformers_covered": False}


def extract_retained_arkane_rate_evidence(
    spec: VerifiedArkaneRouteSpec,
    workdir: Path,
    output_file: Path | None = None,
    *,
    rmg_env: str | None = "rmg_env",
    forward_evaluation: Mapping[str, Any] | None = None,
) -> dict:
    """Reopen central forward/reverse fits without Arkane or electronic-structure jobs.

    A fitted reverse rate belongs to the reached endpoint, not the complete
    conformer pool. Printed tables are rounded diagnostics, never uncertainty
    bounds. This function is read-only and does not admit a network replacement.
    """
    workdir = Path(workdir)
    output_file = Path(output_file) if output_file is not None else workdir / "output.py"
    common: dict[str, Any] = {
        "schema_version": 1, "route_id": spec.label,
        "temperature_K": spec.condition_temperature_K, "pressure_bar": spec.condition_pressure_bar,
        "model_scope": "master_equation_condition_rate" if spec.use_pressure_dependence else "high_pressure_limit",
        "upper_bound_claim_allowed": False, "network_reinjection_claim_allowed": False,
        "all_conformer_return_bound_allowed": False,
        "fit_uncertainty_is_physical_uncertainty": False,
        "uncertainty": {"status": "calibration_required", "lower_rate_multiplier": None,
                        "upper_rate_multiplier": None},
    }
    try:
        basis = _retained_arkane_artifact_basis(spec, workdir)
        common["endpoint_artifact_basis"] = basis
        common["source_artifacts"] = {key: {"path": str(path.resolve()), "sha256": _sha256(path)}
                                      for key, path in (("input", workdir / "input.py"),
                                                        ("output", output_file),
                                                        ("provenance", workdir / "provenance.json"))}
    except (OSError, ValueError, KeyError, TypeError, SyntaxError) as error:
        unresolved = {**common, "status": "unresolved_artifact_basis", "failure_reason": str(error)}
        return {"conditional_forward_rate_evidence": {**unresolved, "kind": "conditional_forward_rate_evidence"},
                "reverse_rate_evidence": {**unresolved, "kind": "reached_endpoint_reverse_rate_evidence"}}
    extracted = _extract_arkane_rate_model(spec, output_file)
    expression = extracted["expression"]
    if not spec.use_pressure_dependence:
        _validated_fit_expression(expression, molecularity=len(spec.reactant_labels))
    evaluated = dict(forward_evaluation) if forward_evaluation is not None else _evaluate_kinetics_expression(
        expression, temperature_K=spec.condition_temperature_K, pressure_bar=spec.condition_pressure_bar,
        pressure_dependent=spec.use_pressure_dependence, rmg_env=rmg_env)
    forward = {**common, "kind": "conditional_forward_rate_evidence", "status": "conditional_central_model_rate",
               "method": "native_RMG_evaluation_of_retained_Arkane_fit", "fit_expression": expression,
               "condition_rate": {"value": _finite_positive(evaluated.get("value_si"), "forward fitted rate"),
                                  "units": "s^-1" if len(spec.reactant_labels) == 1 else "m^3/(mol*s)",
                                  "temperature_K": spec.condition_temperature_K,
                                  "pressure_bar": spec.condition_pressure_bar},
               "native_evaluation": evaluated}
    reverse = {**common, "kind": "reached_endpoint_reverse_rate_evidence", "status": "unresolved_reverse_model",
               "model_scope": "high_pressure_limit", "reverse_reactants": list(spec.product_labels),
               "reverse_products": list(spec.reactant_labels),
               "basis_origin": "Arkane reverse fit to forward TST+tunneling divided by reached-endpoint RRHO Kc",
               "limitation": "Reached-endpoint RRHO detailed balance only; no conformer ensemble, rate calibration, "
                             "or low-pressure applicability is established."}
    if spec.use_pressure_dependence:
        reverse["failure_reason"] = "High-pressure reverse TST diagnostics cannot bound the pressure-dependent net rate"
        return {"conditional_forward_rate_evidence": forward, "reverse_rate_evidence": reverse}
    block = _tst_output_block(spec, output_file)
    fits = re.findall(r"^\s*#\s*k_rev \(TST\+T\)\s*=\s*(.+)$", block, re.MULTILINE)
    if len(fits) != 1:
        reverse["failure_reason"] = "Selected reaction does not contain exactly one retained reverse TST+tunneling fit"
        return {"conditional_forward_rate_evidence": forward, "reverse_rate_evidence": reverse}
    try:
        reverse_expression = _validated_fit_expression(fits[0].strip(), molecularity=len(spec.product_labels))
        reverse_evaluated = _evaluate_kinetics_expression(
            reverse_expression, temperature_K=spec.condition_temperature_K, pressure_bar=spec.condition_pressure_bar,
            pressure_dependent=False, rmg_env=rmg_env)
        reverse_value = _finite_positive(reverse_evaluated.get("value_si"), "reverse fitted rate")
    except (OSError, RuntimeError, ValueError, SyntaxError) as error:
        reverse["failure_reason"] = str(error)
        return {"conditional_forward_rate_evidence": forward, "reverse_rate_evidence": reverse}
    reverse.update({"status": "conditional_central_model_rate",
                    "method": "native_RMG_evaluation_of_retained_Arkane_reverse_fit",
                    "fit_expression": reverse_expression, "native_evaluation": reverse_evaluated,
                    "condition_rate": {"value": reverse_value,
                                       "units": "s^-1" if len(spec.product_labels) == 1 else "m^3/(mol*s)",
                                       "temperature_K": spec.condition_temperature_K,
                                       "pressure_bar": spec.condition_pressure_bar}})
    reverse_table = block.split("Kc (eq)", 1)[-1] if "Kc (eq)" in block else ""
    rows = re.findall(r"^\s*#\s*([\d.eE+-]+) K\s+([\d.eE+-]+)\s+"
                      r"(?:mol\^2/cm\^6|mol/cm\^3|cm\^3/mol|cm\^6/mol\^2)?\s*"
                      r"([\d.eE+-]+)\s+([\d.eE+-]+)\s+(\S+)\s*$", reverse_table, re.MULTILINE)
    condition_rows = [row for row in rows if float(row[0]) == spec.condition_temperature_K]
    if len(condition_rows) == 1:
        row = condition_rows[0]
        reverse["rounded_table_diagnostic"] = {
            "temperature_K": float(row[0]), "equilibrium_constant_printed": float(row[1]),
            "reverse_TST_printed": float(row[2]), "reverse_TST_tunneling_printed": float(row[3]),
            "printed_units": row[4], "precision": "Arkane %.3e rate/equilibrium table formatting",
            "upper_bound_claim_allowed": False,
            "limitation": "Rounded sampled-model values, not confidence bounds or the evaluated fitted rate."}
    return {"conditional_forward_rate_evidence": forward, "reverse_rate_evidence": reverse}


def load_retained_arkane_rate_evidence(workdir: Path, *, rmg_env: str | None = "rmg_env") -> dict:
    """Read-only compatibility adapter for jobs predating structured reverse evidence."""
    workdir = Path(workdir)
    provenance = json.loads((workdir / "provenance.json").read_text())
    route = dict(provenance["route"])
    for name in ("reactant_orca_outputs", "product_orca_outputs"):
        route[name] = tuple(Path(value) for value in route[name])
    route["transition_state_orca_output"] = Path(route["transition_state_orca_output"])
    for name in ("reactant_labels", "reactant_smiles", "reactant_multiplicities", "product_labels",
                 "product_smiles", "product_multiplicities", "temperatures_K", "pressures_bar", "interpolation_model"):
        route[name] = tuple(route[name])
    route["bath_gas"] = tuple(tuple(value) for value in route["bath_gas"])
    return extract_retained_arkane_rate_evidence(VerifiedArkaneRouteSpec(**route), workdir, rmg_env=rmg_env)


def _write_verified_kinetics_library(
    spec: VerifiedArkaneRouteSpec,
    workdir: Path,
    *,
    kinetics_expression: str,
    species_dictionary: Path | None = None,
) -> Path:
    """Build the one-reaction RMG library Arkane itself does not export for pdep net rates."""
    native = Path(workdir) / "RMG_libraries" / "kinetics"
    native_dictionary = species_dictionary or (native / "dictionary.txt")
    if not native_dictionary.is_file():
        raise ValueError("Arkane did not write the species dictionary needed by the verified kinetics library")
    # RMG keys external kinetics libraries by the final directory name.  A
    # generic ``kinetics`` basename for every iteration causes later repairs to
    # overwrite earlier ones in ``database.kinetics.libraries``.  Bind the
    # basename to the stable route ID so all verified replacements survive a
    # multi-route full-mechanism rerun.
    library = Path(workdir) / "verified_RMG_library" / f"storca-verified-{_safe_label(spec.label)}"
    library.mkdir(parents=True, exist_ok=True)
    equation = f"{' + '.join(spec.reactant_labels)} <=> {' + '.join(spec.product_labels)}"
    (library / "reactions.py").write_text(
        "#!/usr/bin/env python\n# encoding: utf-8\n\n"
        "name = 'STORCA verified Arkane route'\n"
        "shortDesc = 'One ORCA/Arkane-verified elementary route'\n"
        "longDesc = '''Generated from validated ORCA stationary points and retained Arkane output.'''\n"
        "autoGenerated = True\n\n"
        "entry(\n"
        "    index = 1,\n"
        f"    label = {equation!r},\n"
        f"    kinetics = {kinetics_expression},\n"
        "    shortDesc = 'STORCA condition-rate repair',\n"
        "    longDesc = '''Do not detach this rate from generated-kinetics.json provenance.''',\n"
        ")\n"
    )
    shutil.copy2(native_dictionary, library / "dictionary.txt")
    return library


def _prepare_arkane_species_dictionary(
    spec: VerifiedArkaneRouteSpec, workdir: Path, *, rmg_env: str | None,
) -> Path:
    """Export exact declared species when Arkane omits its native dictionary.

    This is graph serialization in the owning RMG environment, not discovery,
    thermo estimation, or rate generation. Existing exports are not overwritten.
    """
    native = Path(workdir) / "RMG_libraries" / "kinetics" / "dictionary.txt"
    if native.is_file():
        return native
    rows = list(zip((*spec.reactant_labels, *spec.product_labels),
                    (*spec.reactant_smiles, *spec.product_smiles),
                    (*spec.reactant_multiplicities, *spec.product_multiplicities)))
    script = r'''
import json, sys
import rmgpy
from rmgpy.molecule import Molecule
records = {}
for label, smiles, multiplicity in json.loads(sys.argv[1]):
    molecule = Molecule(smiles=smiles)
    if molecule.multiplicity != multiplicity:
        raise ValueError("Declared spin differs from RMG graph for " + label)
    adjacency = molecule.to_adjacency_list(label=label, remove_h=False)
    restored = Molecule().from_adjacency_list(adjacency)
    if not restored.is_isomorphic(molecule) or restored.multiplicity != multiplicity:
        raise ValueError("Species dictionary roundtrip failed for " + label)
    if label in records and records[label] != adjacency:
        raise ValueError("Conflicting species dictionary label " + label)
    records[label] = adjacency
print("STORCA_ARKANE_DICTIONARY_JSON=" + json.dumps({
    "dictionary": "\n\n".join(records.values()) + "\n",
    "labels": list(records), "rmg_version": getattr(rmgpy, "__version__", None),
}, sort_keys=True))
'''
    command = _rmg_python_command(rmg_env, script, json.dumps(rows))
    result = run_bounded_kinetics_command(command)
    if result.returncode:
        raise ValueError("RMG species dictionary export failed: " + result.stderr.strip())
    report = _marked_json(result.stdout, "STORCA_ARKANE_DICTIONARY_JSON=")
    if set(report["labels"]) != {row[0] for row in rows}:
        raise ValueError("Species dictionary export omitted or added labels")
    directory = Path(workdir) / "storca_species_dictionary"
    directory.mkdir(exist_ok=True)
    destination = directory / "dictionary.txt"
    contents = report["dictionary"]
    if destination.exists() and destination.read_text() != contents:
        raise ValueError("Retained species dictionary differs from current declared species")
    if not destination.exists():
        destination.write_text(contents)
    evidence = {"source": "declared_species_RMG_graph_roundtrip", "species": rows,
                "rmg_version": report["rmg_version"], "dictionary_sha256": _sha256(destination),
                "rate_claim_allowed": False}
    evidence_path = directory / "provenance.json"
    serialized = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    if evidence_path.exists() and evidence_path.read_text() != serialized:
        raise ValueError("Retained species dictionary provenance changed")
    if not evidence_path.exists():
        evidence_path.write_text(serialized)
    return destination


def _validate_library_with_rmg(
    library: Path,
    *,
    temperature_K: float,
    pressure_bar: float,
    pressure_dependent: bool,
    rmg_env: str | None,
) -> dict:
    marker = "STORCA_RMG_LIBRARY_JSON="
    script = r'''
import json, math, sys
import rmgpy
from rmgpy.data.kinetics.database import KineticsDatabase
database = KineticsDatabase()
database.load_libraries(sys.argv[1], libraries=[sys.argv[2]])
library = database.libraries[sys.argv[2]]
reactions = library.get_library_reactions()
if len(reactions) != 1:
    raise ValueError("Verified library must contain exactly one loadable reaction")
reaction = reactions[0]
temperature = float(sys.argv[3])
pressure_pa = float(sys.argv[4]) * 1.0e5
is_pdep = bool(reaction.kinetics.is_pressure_dependent())
value = reaction.get_rate_coefficient(temperature, pressure_pa) if is_pdep else reaction.get_rate_coefficient(temperature)
if not math.isfinite(value) or value <= 0.0:
    raise ValueError("Verified library returned an invalid condition rate")
print("STORCA_RMG_LIBRARY_JSON=" + json.dumps({
    "reaction": str(reaction), "reactants": [species.label for species in reaction.reactants],
    "products": [species.label for species in reaction.products], "value_si": value,
    "model_class": reaction.kinetics.__class__.__name__, "pressure_dependent": is_pdep,
    "rmg_version": getattr(rmgpy, "__version__", None),
}, sort_keys=True))
'''
    command = _rmg_python_command(
        rmg_env, script, str(Path(library).parent.resolve()), Path(library).name,
        repr(temperature_K), repr(pressure_bar),
    )
    result = run_bounded_kinetics_command(command)
    if result.returncode != 0:
        raise RuntimeError(
            f"Generated Arkane library failed RMG loading: {result.stderr.strip() or result.stdout.strip()}"
        )
    report = _marked_json(result.stdout, marker)
    if bool(report["pressure_dependent"]) != pressure_dependent:
        raise ValueError("Generated RMG library rate class does not match the requested Arkane branch")
    return {**report, "command": command}


def validate_verified_arkane_library(
    spec: VerifiedArkaneRouteSpec,
    workdir: Path,
    output_file: Path,
    *,
    rmg_env: str | None = "rmg_env",
) -> dict:
    """Validate Arkane's generated RMG library for the requested rate class."""
    from .generated_kinetics import _checksum, validate_rate_replacement_evidence

    # A syntactically valid Arkane fit is not sufficient: check the native
    # loaded external modes and their hash-bound execution receipt before
    # producing any replacement library. Historical unaudited rates remain
    # diagnostic artifacts, not silently re-certified replacements.
    if Path(output_file).resolve() != (Path(workdir) / "output.py").resolve():
        raise ValueError("Arkane rate output differs from the audited workspace output")
    statmech_validation = validate_rotor_receipt(Path(workdir))
    rate_source_artifacts = {
        name: {"path": str(path.resolve()), "sha256": _sha256(path)}
        for name, path in (("input", Path(workdir) / "input.py"),
                           ("output", Path(output_file)),
                           ("provenance", Path(workdir) / "provenance.json"))
    }

    if spec.conformer_rate_basis is not None:
        retained_manifest = (Path(workdir) / "verified_RMG_library"
                             / f"storca-verified-{_safe_label(spec.label)}" / "generated-kinetics.json")
        if retained_manifest.exists() or retained_manifest.is_symlink():
            raise ValueError(
                "Reached-conformer rate generation requires a new library workspace; "
                "retained generated-kinetics evidence will not be rewritten"
            )
    pdep = parse_arkane_pdep_output(output_file)
    explicit_pressure_model = (
        _validated_explicit_pressure_model(spec) if spec.use_pressure_dependence else None
    )
    if spec.use_pressure_dependence and pdep["status"] != "completed":
        return {"verification_status": "incomplete", "failure_reason": "Arkane output lacks pdep net kinetics",
                "pdep": pdep}
    extracted = _extract_arkane_rate_model(spec, output_file)
    evaluated = _evaluate_kinetics_expression(
        extracted["expression"], temperature_K=spec.condition_temperature_K,
        pressure_bar=spec.condition_pressure_bar, pressure_dependent=spec.use_pressure_dependence,
        rmg_env=rmg_env,
    )
    library = _write_verified_kinetics_library(
        spec, workdir, kinetics_expression=evaluated["canonical_expression"],
        species_dictionary=_prepare_arkane_species_dictionary(spec, workdir, rmg_env=rmg_env),
    )
    library_check = _validate_library_with_rmg(
        library, temperature_K=spec.condition_temperature_K, pressure_bar=spec.condition_pressure_bar,
        pressure_dependent=spec.use_pressure_dependence, rmg_env=rmg_env,
    )
    if Counter(library_check["reactants"]) != Counter(spec.reactant_labels) \
            or Counter(library_check["products"]) != Counter(spec.product_labels):
        raise ValueError("Generated RMG library did not preserve the original reaction direction and stoichiometry")
    if library_check["model_class"] != evaluated["model_class"]:
        raise ValueError("Generated RMG library changed Arkane's kinetics model class")
    if not math.isclose(evaluated["value_si"], library_check["value_si"], rel_tol=1e-10, abs_tol=0.0):
        raise ValueError("Generated RMG library does not reproduce Arkane's condition-specific rate")
    reactions, dictionary = library / "reactions.py", library / "dictionary.txt"
    verification = "verified_arkane_pdep" if spec.use_pressure_dependence else "verified_arkane_tst"
    model = "Arkane pressure dependent" if spec.use_pressure_dependence else "Arkane high-pressure-limit TST"
    molecularity = len(spec.reactant_labels)
    units = "s^-1" if molecularity == 1 else "m^3/(mol*s)"
    library_equation = f"{' + '.join(spec.reactant_labels)} <=> {' + '.join(spec.product_labels)}"
    equation = str(spec.reaction_equation or library_equation)
    manifest = {
        "schema_version": 2,
        "kind": "storca_generated_kinetics",
        "source": "arkane_pdep" if spec.use_pressure_dependence else "arkane_tst",
        "arkane_statmech_validation": statmech_validation,
        "arkane_rate_source_artifacts": rate_source_artifacts,
        "route_id": spec.label,
        "reactants": list(spec.reactant_labels),
        "products": list(spec.product_labels),
        "reaction_equation": equation,
        "replaces_reaction": equation,
        "library_reaction_equation": library_equation,
        "purpose": "kinetic_repair",
        "molecularity": molecularity,
        "kinetics_model": model,
        "rate_model_class": library_check["model_class"],
        "temperature_range_K": [min(spec.temperatures_K), max(spec.temperatures_K)],
        "pressure_range_bar": (
            [min(spec.pressures_bar), max(spec.pressures_bar)] if spec.use_pressure_dependence else None
        ),
        "condition_rate": {
            "value": library_check["value_si"], "units": units,
            "temperature_K": spec.condition_temperature_K, "pressure_bar": spec.condition_pressure_bar,
        },
        "collision_limit": (
            {"value": spec.collision_limit_m3_mol_s, "units": units, "source": spec.collision_limit_source}
            if molecularity == 2 else {"applicable": False}
        ),
        "pressure_interpretation": (
            "master_equation_condition_rate" if spec.use_pressure_dependence
            else "high_pressure_limit_TST; recorded pressure is the condition contract but is not an input to k(T)"
        ),
        "pressure_dependence_model": (
            {
                "status": (
                    "explicit_applicable_collisional_model"
                    if explicit_pressure_model is not None else
                    "central_model_with_unvalidated_collisional_estimates"
                ),
                "master_equation_method": spec.method,
                "bath_gas": dict(spec.bath_gas),
                "energy_transfer_model": (
                    explicit_pressure_model["energy_transfer_model"]
                    if explicit_pressure_model is not None else
                    "RMG generic SingleExponentialDown estimate"
                ),
                "transport_model": (
                    explicit_pressure_model["transport_models"]
                    if explicit_pressure_model is not None else
                    "RMG estimated transport from declared structures"
                ),
                "maximum_grain_size_kcal_mol": spec.maximum_grain_size_kcal_mol,
                "minimum_grain_count": spec.minimum_grain_count,
                "applicability_claimed": explicit_pressure_model is not None,
                "explicit_model_sha256": (
                    explicit_pressure_model["model_sha256"]
                    if explicit_pressure_model is not None else None
                ),
                "required_for_prediction": (
                    "explicit applicable uncertainty or a future calibrated collisional model"
                ),
            }
            if spec.use_pressure_dependence else {"applicable": False}
        ),
        "rate_uncertainty": {
            "status": "calibration_required",
            "lower_rate_multiplier": None,
            "upper_rate_multiplier": None,
            "reason": (
                "Arkane produced a central model rate, but STORCA has no calibrated "
                "default uncertainty for this reaction, electronic-structure method, "
                "statistical-mechanics treatment, and pressure-dependence model."
            ),
            "claim_effect": (
                "A nonzero-flux calculated rate requires explicit condition-matched "
                "uncertainty propagation before an observable can converge."
            ),
        },
        "library": str(library.resolve()),
        "checksums": {
            "reactions.py": _checksum(reactions),
            "dictionary.txt": _checksum(dictionary),
        },
        "verification_status": verification,
        "pdep": pdep,
        "arkane_rate_evaluation": evaluated,
        "rmg_library_validation": library_check,
    }
    manifest.update(extract_retained_arkane_rate_evidence(
        spec, workdir, output_file, rmg_env=rmg_env, forward_evaluation=evaluated,
    ))
    physical = validate_rate_replacement_evidence(
        manifest, temperature_K=spec.condition_temperature_K, pressure_bar=spec.condition_pressure_bar,
    )
    if physical["status"] != "passed":
        rejected = {**manifest, "verification_status": "incomplete_physical_validation",
                    "physical_validation": physical}
        evidence_path = Path(workdir) / "arkane-rate-validation.json"
        evidence_path.write_text(json.dumps(rejected, indent=2, sort_keys=True) + "\n")
        return {**rejected, "failure_reason": ", ".join(physical["blocking_reasons"]),
                "validation_evidence": str(evidence_path)}
    manifest["physical_validation"] = physical
    if spec.conformer_rate_basis is not None:
        # Do not prevent Arkane from retaining the conditional elementary
        # rate.  Network admission is a distinct scientific decision, checked
        # both by the verification engine and direct generated-library use.
        manifest["conformer_rate_basis"] = dict(spec.conformer_rate_basis)
        network_validation = validate_rate_replacement_evidence(
            manifest, temperature_K=spec.condition_temperature_K,
            pressure_bar=spec.condition_pressure_bar,
        )
        manifest["network_reinjection_validation"] = network_validation
        manifest["network_reinjection_allowed"] = network_validation["status"] == "passed"
        manifest["conformer_thermochemistry_limitation"] = (
            "Coordinate agreement with a reference minimum is not a calibration of "
            "RMG species thermochemistry. A distinct reached conformer is not assumed "
            "to equilibrate rapidly with the graph-only network's species pool."
        )
    manifest_path = library / "generated-kinetics.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {**manifest, "manifest": str(manifest_path)}


def run_verified_arkane_route(
    spec: VerifiedArkaneRouteSpec, workdir: Path, *, rmg_env: str | None = "rmg_env"
) -> dict:
    """Execute family-neutral TST/pdep kinetics and retain fail-closed evidence."""
    workdir = Path(workdir)
    artifacts = {
        "input": str(workdir / "input.py"),
        "output": str(workdir / "output.py"),
        "stdout": str(workdir / "arkane.stdout.log"),
        "stderr": str(workdir / "arkane.stderr.log"),
    }
    try:
        if any((workdir / name).exists() or (workdir / name).is_symlink() for name in
               ("input.py", "provenance.json", "output.py", "arkane-rotor-contract.json",
                "arkane-statmech-validation.json", "arkane-statmech-validation.failed.json")):
            raise ValueError(
                "Arkane has retained result or statmech evidence; use a new rate workspace. "
                "Existing calculations and evidence will not be rewritten."
            )
        input_file = create_verified_arkane_input(spec, workdir)
        syntax_validation = validate_arkane_input_syntax(input_file, rmg_env=rmg_env)
        contract = create_rotor_contract(workdir)
        artifacts["statmech_contract"] = str(contract)
        artifacts["statmech_validation"] = str(workdir / "arkane-statmech-validation.json")
        execution = run_arkane(input_file, rmg_env=rmg_env, statmech_contract=contract)
        Path(artifacts["stdout"]).write_text(execution["stdout"])
        Path(artifacts["stderr"]).write_text(execution["stderr"])
        output = Path(artifacts["output"])
        if not output.is_file():
            return {"status": "incomplete", "artifacts": artifacts,
                    "failure_reason": "Arkane completed without output.py"}
        library = validate_verified_arkane_library(spec, workdir, output, rmg_env=rmg_env)
        status = "completed" if library["verification_status"].startswith("verified") else "incomplete"
        return {"schema_version": 2, "status": status, "artifacts": artifacts,
                "input_syntax_validation": syntax_validation, "execution": execution,
                "generated_kinetics_library": library}
    except Exception as error:
        if getattr(error, "stdout", None) is not None:
            Path(artifacts["stdout"]).write_text(error.stdout)
        if getattr(error, "stderr", None) is not None:
            Path(artifacts["stderr"]).write_text(error.stderr)
        return {"schema_version": 2, "status": "failed", "artifacts": artifacts,
                "failure_reason": str(error)}


def validate_arkane_library(spec: ArkaneRouteSpec, workdir: Path, output_file: Path) -> dict:
    """Validate Arkane's native RMG library before RMG may load it."""
    from .generated_kinetics import _checksum
    library = Path(workdir) / "RMG_libraries" / "kinetics"
    reactions, dictionary = library / "reactions.py", library / "dictionary.txt"
    pdep = parse_arkane_pdep_output(output_file)
    if not reactions.is_file() or not dictionary.is_file():
        return {"verification_status": "incomplete", "failure_reason": "Arkane did not write an RMG kinetics library",
                "library": str(library), "pdep": pdep}
    text = reactions.read_text()
    # Do not silently admit a high-pressure-only export for a pdep stage.
    if pdep["status"] != "completed" or not ("Chebyshev(" in text or "PDepArrhenius(" in text):
        return {"verification_status": "incomplete", "failure_reason": "Arkane library lacks an exported pressure-dependent kinetics model",
                "library": str(library), "pdep": pdep}
    manifest = {
        "kind": "storca_generated_kinetics", "source": "arkane_pdep", "route_id": spec.label,
        "reactants": [spec.reactant_label], "products": list(spec.product_labels),
        "kinetics_model": "Arkane pressure dependent", "temperature_range_K": [min(spec.temperatures_K), max(spec.temperatures_K)],
        "pressure_range_bar": [min(spec.pressures_bar), max(spec.pressures_bar)], "library": str(library.resolve()),
        "checksums": {"reactions.py": _checksum(reactions), "dictionary.txt": _checksum(dictionary)},
        "verification_status": "verified_arkane_pdep", "pdep": pdep,
    }
    manifest_path = library / "generated-kinetics.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {**manifest, "manifest": str(manifest_path)}


def run_arkane_route(spec: ArkaneRouteSpec, workdir: Path, *, rmg_env: str | None = "rmg_env") -> dict:
    """Execute an Arkane route calculation and retain all artifacts and failures."""
    artifacts = {"input": str(Path(workdir) / "input.py"), "output": str(Path(workdir) / "output.py"),
                 "stdout": str(Path(workdir) / "arkane.stdout.log"), "stderr": str(Path(workdir) / "arkane.stderr.log")}
    try:
        input_file = create_arkane_input(spec, workdir)
        artifacts["input"] = str(input_file)
        execution = run_arkane(input_file, rmg_env=rmg_env)
        Path(artifacts["stdout"]).write_text(execution["stdout"])
        Path(artifacts["stderr"]).write_text(execution["stderr"])
        output = Path(artifacts["output"])
        if not output.is_file():
            return {"status": "incomplete", "artifacts": artifacts,
                    "failure_reason": "Arkane completed without output.py"}
        library = validate_arkane_library(spec, workdir, output)
        status = "completed" if library["verification_status"].startswith("verified") else "incomplete"
        return {"status": status, "artifacts": artifacts, "execution": execution,
                "pdep": library.get("pdep"), "generated_kinetics_library": library}
    except Exception as error:
        if getattr(error, "stdout", None) is not None:
            Path(artifacts["stdout"]).write_text(error.stdout)
        if getattr(error, "stderr", None) is not None:
            Path(artifacts["stderr"]).write_text(error.stderr)
        return {"status": "failed", "artifacts": artifacts, "failure_reason": str(error)}
