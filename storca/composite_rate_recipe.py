"""Explicit, source-bound composite energies for isolated gas-phase Arkane work.

Electronic energy and geometry/vibration evidence are independent bindings.
This module never launches programs, modifies native outputs, infers frequency
scales/symmetries, or promotes a diagnostic into a repaired network library.
Local consistency checks do not certify an IRC, network coverage, a physical
rate, harmonic validity, or electronic-state accuracy.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from src.stability.freq_check import parse_orca_frequencies
from .arkane_rotor import POLICY_ID, create_rotor_contract
from .arkane_runner import _smiles_identity
from .ir_modes import parse_orca_normal_modes
from .route_geometry import _coupled_multiplicities
from .route_thermochemistry import HARTREE_TO_J_MOL


CONTRACT_ID = "isolated_gas_composite_rate_recipe_v1"
ATOMIC_CONTRACT_ID = "isolated_gas_composite_rate_recipe_v2_atomic_hydrogen"
ATOMIC_ENERGY_POLICY = "one_electron_hydrogen_native_hf_same_basis_v1"
BOHR_TO_ANGSTROM = 0.529177210903
CM_INV_TO_J_MOL = 11.9626565638697
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


@dataclass(frozen=True)
class ArtifactBinding:
    path: str
    sha256: str

    def record(self) -> dict:
        return asdict(self)


def bind_artifact(path: Path) -> ArtifactBinding:
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError("Composite recipe artifact is not a file")
    return ArtifactBinding(str(path), hashlib.sha256(path.read_bytes()).hexdigest())


def reopen(binding: ArtifactBinding) -> Path:
    if (not isinstance(binding, ArtifactBinding) or not isinstance(binding.path, str)
            or not Path(binding.path).is_absolute() or not isinstance(binding.sha256, str)
            or not re.fullmatch(r"[0-9a-f]{64}", binding.sha256)):
        raise ValueError("Invalid immutable composite artifact binding")
    try:
        actual = bind_artifact(Path(binding.path))
    except OSError as error:
        raise ValueError("Bound composite artifact is absent") from error
    if actual != binding:
        raise ValueError("Composite source artifact hash changed: " + binding.path)
    return Path(binding.path)


@dataclass(frozen=True)
class CompositeStationaryPoint:
    label: str
    kind: str  # species or transition_state
    frequency_input: ArtifactBinding
    frequency_output: ArtifactBinding
    hessian: ArtifactBinding
    energy_input: ArtifactBinding
    energy_output: ArtifactBinding
    energy_hartree: float  # Electronic energy only; independently checked below.
    charge: int
    multiplicity: int
    external_symmetry: int
    optical_isomers: int
    smiles: str | None


@dataclass(frozen=True)
class CompositeAtomicSpecies:
    """One explicitly declared H atom, with no optimizer or Hessian fiction.

    The one-electron HF calculation uses the molecular composite's *same*
    orbital basis. Its finite-basis native energy, not -0.5 Eh, is retained.
    No automatic fallback from a failed correlated calculation is permitted.
    """
    label: str
    energy_input: ArtifactBinding
    energy_output: ArtifactBinding
    energy_hartree: float
    charge: int
    multiplicity: int
    external_symmetry: int
    optical_isomers: int
    smiles: str
    atomic_mass_amu: float
    isotope: int
    atomic_energy_method: str
    kind: str = "species"
    atomic_energy_policy: str = ATOMIC_ENERGY_POLICY


def bind_atomic_species(label: str, energy_output: Path, *, energy_hartree: float,
                        atomic_mass_amu: float, isotope: int = 1,
                        atomic_energy_method: str = "HF/def2-QZVPP",
                        smiles: str = "[H]") -> CompositeAtomicSpecies:
    """Bind a separately authorized native atomic single point; never execute."""
    energy_output = Path(energy_output)
    return CompositeAtomicSpecies(label, bind_artifact(energy_output.with_suffix(".inp")),
                                  bind_artifact(energy_output), energy_hartree, 0, 2,
                                  1, 1, smiles, atomic_mass_amu, isotope,
                                  atomic_energy_method)


def bind_stationary_point(
    label: str, kind: str, frequency_output: Path, energy_output: Path, *,
    energy_hartree: float, charge: int, multiplicity: int,
    external_symmetry: int, optical_isomers: int, smiles: str | None,
) -> CompositeStationaryPoint:
    """Bind native files once; subsequent validation refuses changed bytes."""
    frequency_output, energy_output = Path(frequency_output), Path(energy_output)
    return CompositeStationaryPoint(
        label, kind, bind_artifact(frequency_output.with_suffix(".inp")),
        bind_artifact(frequency_output), bind_artifact(frequency_output.with_suffix(".hess")),
        bind_artifact(energy_output.with_suffix(".inp")), bind_artifact(energy_output),
        energy_hartree, charge, multiplicity, external_symmetry, optical_isomers, smiles,
    )


@dataclass(frozen=True)
class CompositeRateRecipe:
    label: str
    points: tuple[CompositeStationaryPoint | CompositeAtomicSpecies, ...]
    reactant_labels: tuple[str, ...]
    product_labels: tuple[str, ...]
    transition_state_label: str
    electronic_method: str  # Exact input method/basis, e.g. CCSD(T)/def2-QZVPP.
    vibrational_method: str  # Exact input method/basis, e.g. WB97X-D3/def2-TZVPP.
    frequency_scale_factor: float  # Arkane scale applied once to raw frequencies.
    temperatures_K: tuple[float, ...]
    tunneling: str  # none, Wigner or Eckart; deliberately no inferred default.
    environment: str
    surface_model: str
    electronic_state_assumption: str
    channel_counting_assumption: str
    soft_mode_assumption: str
    route_evidence: tuple[ArtifactBinding, ...]


def _finite(value, label, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) \
            or (positive and value <= 0):
        raise ValueError("Invalid explicit composite " + label)
    return float(value)


def _integer(value, label, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, int) or (positive and value < 1):
        raise ValueError("Invalid explicit composite " + label)


def _input(path: Path, method: str) -> dict:
    text = re.sub(r"#[^\n]*", "", path.read_text(errors="replace"))
    if not isinstance(method, str):
        raise ValueError("Composite method must explicitly name method/basis")
    parts = method.split("/")
    if len(parts) != 2 or any(not part or re.search(r"\s", part) for part in parts):
        raise ValueError("Composite method must explicitly name method/basis")
    keywords = " ".join(re.findall(r"(?m)^\s*!([^\n]*)", text)).split()
    lowered = {word.casefold() for word in keywords}
    if not all(part.casefold() in lowered for part in parts):
        raise ValueError("Native input differs from declared composite method: " + str(path))
    if (any(word.casefold().startswith(("opt", "neb", "irc")) for word in keywords)
            or re.search(r"(?i)\b(?:CPCM|SMD|ALPB|COSMO|QMMM|PBC)\b|%\s*(?:cpcm|tddft|casscf|pointcharges|qmmm)\b", text)):
        raise ValueError("Composite diagnostic requires separate ground-state gas single-point/frequency jobs")
    blocks = re.findall(r"(?ims)^\s*\*\s+xyz\s+([-+]?\d+)\s+(\d+)\s*\n(.*?)^\s*\*\s*$", text)
    if len(blocks) != 1:
        raise ValueError("Exact one-block inline Cartesian native input required")
    charge, spin, coordinates = blocks[0]
    rows = [line.split() for line in coordinates.splitlines() if line.strip()]
    if not rows or any(len(row) != 4 or not re.fullmatch(r"[A-Z][a-z]?", row[0]) for row in rows):
        raise ValueError("Malformed native input geometry")
    try:
        xyz = np.asarray([[float(value) for value in row[1:]] for row in rows])
    except ValueError as error:
        raise ValueError("Non-numeric native input geometry") from error
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite native input geometry")
    return {"charge": int(charge), "multiplicity": int(spin), "labels": tuple(row[0] for row in rows),
            "xyz": xyz, "keywords": lowered, "text": text}


def _native(path: Path, reference: dict) -> tuple[str, float]:
    text = path.read_text(errors="replace")
    if (text.count("ORCA TERMINATED NORMALLY") != 1 or "SCF CONVERGED AFTER" not in text
            or re.search(r"ORCA.*(?:error termination|TERMINATED ABNORMALLY)", text, re.I)):
        raise ValueError("Composite source requires unique normal termination and SCF convergence")
    versions = re.findall(r"(?m)^\s*Program Version\s+(\d+\.\d+\.\d+)\b", text)
    if len(versions) != 1 or not versions[0].startswith("6.1."):
        raise ValueError("Composite native evidence parser is restricted to ORCA 6.1")
    charges = re.findall(r"(?m)^\s*Total Charge\s+Charge\s+\.{2,}\s+([-+]?\d+)\s*$", text)
    spins = re.findall(r"(?m)^\s*Multiplicity\s+Mult\s+\.{2,}\s+(\d+)\s*$", text)
    if (not charges or not spins or {int(value) for value in charges} != {reference["charge"]}
            or {int(value) for value in spins} != {reference["multiplicity"]}):
        raise ValueError("Native output/input charge or multiplicity mismatch")
    geometries = []
    sections = text.split("CARTESIAN COORDINATES (ANGSTROEM)")[1:]
    for section in sections:
        rows = []
        for line in section.splitlines():
            fields = line.split()
            if len(fields) == 4 and re.fullmatch(r"[A-Z][a-z]?", fields[0]):
                try:
                    rows.append((fields[0], [float(value) for value in fields[1:]]))
                except ValueError as error:
                    raise ValueError("Malformed native output geometry") from error
            elif rows:
                break
        if rows:
            geometries.append(rows)
    if not geometries:
        raise ValueError("Native output geometry is absent")
    for rows in geometries:
        if (tuple(row[0] for row in rows) != reference["labels"]
                or not np.allclose([row[1] for row in rows], reference["xyz"], atol=1e-5, rtol=0)):
            raise ValueError("Native output/input geometry or nucleus order mismatch")
    energies = re.findall(rf"(?m)^\s*FINAL SINGLE POINT ENERGY\s+({NUMBER})\s*$", text)
    if len(energies) != 1 or not math.isfinite(float(energies[0].replace("D", "E").replace("d", "e"))):
        raise ValueError("Unique finite native electronic energy required")
    return text, float(energies[0].replace("D", "E").replace("d", "e"))


def _full_hessian(path: Path, count: int):
    """Reject missing/truncated matrices before native statmech sees them."""
    lines = path.read_text(errors="replace").splitlines()
    if lines.count("$hessian") != 1:
        raise ValueError("Unique full native Hessian matrix required")
    index = lines.index("$hessian") + 1
    if index >= len(lines):
        raise ValueError("Missing full Hessian dimension")
    dimension = int(lines[index]); index += 1
    if dimension != 3 * count:
        raise ValueError("Full native Hessian dimension mismatch")
    matrix = np.full((dimension, dimension), np.nan)
    while index < len(lines) and not lines[index].lstrip().startswith("$"):
        if not lines[index].strip():
            index += 1
            continue
        columns = [int(value) for value in lines[index].split()]; index += 1
        if not columns or len(set(columns)) != len(columns) or any(column < 0 or column >= dimension for column in columns):
            raise ValueError("Invalid full Hessian column block")
        for _ in range(dimension):
            if index >= len(lines) or lines[index].lstrip().startswith("$"):
                raise ValueError("Incomplete full Hessian row")
            fields = lines[index].split(); index += 1
            if len(fields) != len(columns) + 1:
                raise ValueError("Incomplete full Hessian row")
            row = int(fields[0])
            if row < 0 or row >= dimension or np.isfinite(matrix[row, columns]).any():
                raise ValueError("Duplicate/invalid full Hessian row")
            matrix[row, columns] = [float(value.replace("D", "E").replace("d", "e")) for value in fields[1:]]
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T, atol=1e-8, rtol=0):
        raise ValueError("Incomplete/nonfinite/asymmetric full native Hessian")


def _atomic_point(point: CompositeAtomicSpecies, recipe: CompositeRateRecipe) -> dict:
    _integer(point.charge, "atomic charge")
    for value, label in ((point.multiplicity, "atomic multiplicity"),
                         (point.external_symmetry, "atomic external symmetry"),
                         (point.optical_isomers, "atomic optical-isomer count")):
        _integer(value, label, positive=True)
    if (point.kind != "species" or not isinstance(point.label, str) or not point.label.strip()
            or point.atomic_energy_policy != ATOMIC_ENERGY_POLICY
            or (point.charge, point.multiplicity, point.external_symmetry, point.optical_isomers) != (0, 2, 1, 1)
            or any(isinstance(value, bool) for value in (point.charge, point.multiplicity,
                    point.external_symmetry, point.optical_isomers, point.isotope))
            or point.isotope not in (0, 1)
            or _smiles_identity(point.smiles) != ({"H": 1}, 0)):
        raise ValueError("Atomic composite evidence currently admits only explicitly declared neutral ground-state protium")
    mass = _finite(point.atomic_mass_amu, "atomic mass", positive=True)
    # Native Arkane's ORCA geometry reader chooses the most abundant isotope.
    # Refuse arbitrary/natural-average masses instead of silently changing it.
    if not math.isclose(mass, 1.00782503224, abs_tol=1e-8, rel_tol=0):
        raise ValueError("Atomic mass differs from the supported native protium mass")
    energy = _finite(point.energy_hartree, "atomic electronic energy")
    ei, eo = reopen(point.energy_input), reopen(point.energy_output)
    if eo.with_suffix(".inp") != ei:
        raise ValueError("Atomic composite log/input basename binding mismatch")
    declared = point.atomic_energy_method.split("/")
    composite = recipe.electronic_method.split("/")
    if (len(declared) != 2 or len(composite) != 2 or declared[0].casefold() != "hf"
            or composite[0].casefold() != "ccsd(t)"
            or declared[1].casefold() != composite[1].casefold()):
        raise ValueError("Atomic one-electron policy requires explicitly declared HF in the identical CCSD(T) composite basis")
    reference = _input(ei, point.atomic_energy_method)
    if (reference["labels"] != ("H",) or (reference["charge"], reference["multiplicity"]) != (0, 2)
            or {"freq", "anfreq", "numfreq", "ccsd(t)", "mp2"} & reference["keywords"]):
        raise ValueError("Atomic evidence must be a one-electron HF single point, not a fictitious frequency/correlated job")
    native, observed = _native(eo, reference)
    electrons = re.findall(r"Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)\b", native)
    if not electrons or set(electrons) != {"1"}:
        raise ValueError("Atomic native one-electron inventory is absent or inconsistent")
    if not math.isclose(observed, energy, abs_tol=1e-10, rel_tol=0):
        raise ValueError("Atomic declared energy differs from the actual finite-basis native energy")
    reopen(point.energy_input); reopen(point.energy_output)
    return {**asdict(point), "evidence_kind": "monatomic_single_point",
            "formula": {"H": 1}, "atom_count": 1,
            "observed_electronic_energy_hartree": observed,
            "external_projected_mode_count": 3,
            "negative_physical_mode_indices": [], "real_physical_frequencies_cm_inv": [],
            "imaginary_frequency_cm_inv": None, "raw_harmonic_ZPE_J_mol": 0.0,
            "frequency_calculation": "not_applicable_monatomic",
            "coordinates_angstrom": reference["xyz"].tolist(),
            "statmech_policy": "native_SP_geometry_translation_and_spin_only_no_Hessian",
            "stationary_character_scope": "A single nucleus has no internal geometry or vibrational stationary-point test."}


def _point(point: CompositeStationaryPoint, recipe: CompositeRateRecipe) -> dict:
    if isinstance(point, CompositeAtomicSpecies):
        return _atomic_point(point, recipe)
    if not isinstance(point, CompositeStationaryPoint) or point.kind not in {"species", "transition_state"} \
            or not isinstance(point.label, str) or not point.label.strip():
        raise ValueError("Invalid composite stationary point declaration")
    _integer(point.charge, "charge")
    for value, name in ((point.multiplicity, "multiplicity"), (point.external_symmetry, "external symmetry"),
                        (point.optical_isomers, "optical-isomer count")):
        _integer(value, name, positive=True)
    energy = _finite(point.energy_hartree, "electronic energy")
    for name in ("frequency_input", "frequency_output", "hessian", "energy_input", "energy_output"):
        reopen(getattr(point, name))
    fi, fo, hess = map(Path, (point.frequency_input.path, point.frequency_output.path, point.hessian.path))
    ei, eo = map(Path, (point.energy_input.path, point.energy_output.path))
    if fo.with_suffix(".inp") != fi or fo.with_suffix(".hess") != hess or eo.with_suffix(".inp") != ei:
        raise ValueError("Composite log/input/Hessian basename binding mismatch")
    freq = _input(fi, recipe.vibrational_method)
    elec = _input(ei, recipe.electronic_method)
    if (not {"freq", "anfreq", "numfreq"} & freq["keywords"]
            or not re.search(r"(?is)%freq\b.*?\bProjectTR\s+true\b.*?\bend\b", freq["text"])):
        raise ValueError("Explicit projected native frequency input required")
    if (freq["labels"] != elec["labels"] or freq["xyz"].shape != elec["xyz"].shape
            or not np.allclose(freq["xyz"], elec["xyz"], atol=1e-10, rtol=0)):
        raise ValueError("Composite electronic and vibration input geometry/order mismatch")
    for source in (freq, elec):
        if (source["charge"], source["multiplicity"]) != (point.charge, point.multiplicity):
            raise ValueError("Composite source/declaration charge or multiplicity mismatch")
    native_freq, frequency_energy = _native(fo, freq)
    native_energy, observed_energy = _native(eo, elec)
    if not math.isclose(energy, observed_energy, abs_tol=1e-10, rel_tol=0):
        raise ValueError("Declared composite electronic energy differs from bound native energy")
    if recipe.electronic_method.split("/")[0].casefold() == "ccsd(t)":
        totals = re.findall(rf"E\(CCSD\(T\)\)\s+\.{{2,}}\s+({NUMBER})", native_energy)
        if ("--- The Coupled-Cluster iterations have converged ---" not in native_energy
                or len(totals) != 1 or not math.isclose(float(totals[0].replace("D", "E").replace("d", "e")), energy, abs_tol=1e-8, rel_tol=0)):
            raise ValueError("Composite CCSD(T) completion/total-energy evidence missing")
    modes = parse_orca_normal_modes(hess)
    count = len(freq["labels"])
    hessian_xyz = modes.coordinates_bohr * BOHR_TO_ANGSTROM
    if (modes.atom_labels != freq["labels"] or len(modes.frequencies_cm_1) != 3 * count
            or not np.allclose(hessian_xyz - hessian_xyz.mean(0), freq["xyz"] - freq["xyz"].mean(0), atol=1e-5, rtol=0)
            or not np.allclose(parse_orca_frequencies(fo), modes.frequencies_cm_1, atol=.05, rtol=0)):
        raise ValueError("Composite frequency log/Hessian/geometry/order mismatch")
    _full_hessian(hess, count)
    hessian_scales = re.findall(r"\$frequency_scale_factor\s+([\d.]+)", hess.read_text())
    printed_scales = re.findall(rf"Scaling factor for frequencies\s*=\s*({NUMBER})\s+\(already applied!\)", native_freq)
    if len(hessian_scales) != 1 or float(hessian_scales[0]) != 1 or len(printed_scales) != 1 or float(printed_scales[0]) != 1:
        raise ValueError("Raw scale1 native frequencies required; Arkane scaling must be applied once")
    rank = np.linalg.svd(freq["xyz"] - freq["xyz"].mean(0), compute_uv=False)
    external = 3 if count == 1 else (5 if len(rank) < 2 or rank[1] < 1e-8 else 6)
    physical = modes.frequencies_cm_1[external:]
    negative = (np.flatnonzero(physical < 0) + external).tolist()
    expected = int(point.kind == "transition_state")
    if (not np.all(np.abs(modes.frequencies_cm_1[:external]) < .01) or np.any(physical == 0)
            or len(negative) != expected):
        raise ValueError("Composite physical stationary-mode inventory unresolved; no negative-mode waiver")
    formula = dict(sorted(Counter(freq["labels"]).items()))
    if point.kind == "species":
        if not isinstance(point.smiles, str) or _smiles_identity(point.smiles) != (formula, point.charge):
            raise ValueError("Composite declared species composition/charge differs from native nuclei")
    elif point.smiles is not None:
        raise ValueError("Transition-state SMILES is not a stationary-point identity certificate")
    # Recheck bindings after every parser read: a changed source cannot supply a recipe.
    for name in ("frequency_input", "frequency_output", "hessian", "energy_input", "energy_output"):
        reopen(getattr(point, name))
    return {**asdict(point), "formula": formula, "atom_count": count,
            "observed_frequency_energy_hartree": frequency_energy,
            "observed_electronic_energy_hartree": observed_energy,
            "external_projected_mode_count": external, "negative_physical_mode_indices": negative,
            "real_physical_frequencies_cm_inv": physical[physical > 0].tolist(),
            "imaginary_frequency_cm_inv": float(physical[physical < 0][0]) if expected else None,
            "raw_harmonic_ZPE_J_mol": float(np.maximum(physical, 0).sum() * CM_INV_TO_J_MOL / 2),
            "frequency_calculation": "numerical" if "numfreq" in freq["keywords"] else "analytical",
            "geometry_energy_matching_tolerance_angstrom": 1e-10,
            "native_log_geometry_rounding_tolerance_angstrom": 1e-5,
            "hessian_geometry_comparison": "Exact nucleus order and centered Cartesian frame; uniform translation only is ignored.",
            "state_evidence": "Bound native input/output declarations; Hessian multiplicity metadata is not authoritative.",
            "stationary_character_scope": "Local physical mode counts only; meaningful reaction coordinate and IRC remain independent route evidence."}


def validate_composite_rate_recipe(recipe: CompositeRateRecipe) -> dict:
    """Pure admission of native-data consistency, not scientific rate accuracy."""
    if (not isinstance(recipe, CompositeRateRecipe) or not isinstance(recipe.label, str)
            or not recipe.label.strip()):
        raise ValueError("Explicit composite recipe declaration required")
    if recipe.environment != "gas" or recipe.surface_model != "single_adiabatic_surface":
        raise ValueError("Composite diagnostic currently supports declared gas single-adiabatic-surface models only")
    for name in ("electronic_state_assumption", "channel_counting_assumption", "soft_mode_assumption"):
        if not isinstance(getattr(recipe, name), str) or not getattr(recipe, name).strip():
            raise ValueError("Explicit composite assumption required: " + name)
    _finite(recipe.frequency_scale_factor, "frequency scale", positive=True)
    if (not isinstance(recipe.points, tuple) or not isinstance(recipe.reactant_labels, tuple)
            or not isinstance(recipe.product_labels, tuple) or not isinstance(recipe.temperatures_K, tuple)
            or not isinstance(recipe.route_evidence, tuple) or not recipe.route_evidence):
        raise ValueError("Immutable tuples and bound route evidence required")
    if str(recipe.tunneling).casefold() not in {"none", "wigner", "eckart"}:
        raise ValueError("Explicit tunneling must be none, Wigner or Eckart")
    if len(recipe.temperatures_K) < 3 or tuple(sorted(set(recipe.temperatures_K))) != recipe.temperatures_K:
        raise ValueError("Composite rate grid needs at least three strictly increasing temperatures")
    for value in recipe.temperatures_K:
        _finite(value, "temperature", positive=True)
    points = [_point(point, recipe) for point in recipe.points]
    if any(isinstance(point, CompositeAtomicSpecies) for point in recipe.points):
        # The new policy is intentionally stronger than the historical v1
        # contract, whose source bytes/admission semantics remain unchanged.
        from .native_rate_audit import read_electronic_single_point
        for point in recipe.points:
            read_electronic_single_point(point.energy_input.path, point.energy_output.path,
                method=point.atomic_energy_method if isinstance(point, CompositeAtomicSpecies) else recipe.electronic_method,
                charge=point.charge, multiplicity=point.multiplicity)
    by_label = {point["label"]: point for point in points}
    if len(by_label) != len(points) or not recipe.reactant_labels or not recipe.product_labels:
        raise ValueError("Unique stationary labels and nonempty reaction sides required")
    if set(by_label) != set((*recipe.reactant_labels, *recipe.product_labels, recipe.transition_state_label)):
        raise ValueError("Composite reaction/stationary-point inventory mismatch")
    if by_label[recipe.transition_state_label]["kind"] != "transition_state" \
            or sum(point["kind"] == "transition_state" for point in points) != 1:
        raise ValueError("Exactly one declared composite transition state required")
    for side in (recipe.reactant_labels, recipe.product_labels):
        if any(by_label[label]["kind"] != "species" for label in side):
            raise ValueError("Composite reaction sides must contain species")
    totals = {}
    for name, labels in (("reactants", recipe.reactant_labels), ("products", recipe.product_labels)):
        formula = Counter()
        for label in labels:
            formula.update(by_label[label]["formula"])
        totals[name] = (dict(sorted(formula.items())), sum(by_label[label]["charge"] for label in labels))
    ts = by_label[recipe.transition_state_label]
    if totals["reactants"] != totals["products"] or totals["reactants"] != (ts["formula"], ts["charge"]):
        raise ValueError("Composite reaction and TS are not atom/charge balanced")
    common = set(_coupled_multiplicities(tuple(by_label[label]["multiplicity"] for label in recipe.reactant_labels))) \
        & set(_coupled_multiplicities(tuple(by_label[label]["multiplicity"] for label in recipe.product_labels)))
    if ts["multiplicity"] not in common:
        raise ValueError("Composite TS is not on a common endpoint spin surface")
    for evidence in recipe.route_evidence:
        reopen(evidence)
    reactant_energy = sum(by_label[label]["energy_hartree"] for label in recipe.reactant_labels)
    gauge = reactant_energy / ts["atom_count"]
    shifted = {label: (point["energy_hartree"] - gauge * point["atom_count"]) * HARTREE_TO_J_MOL
               for label, point in by_label.items()}
    barriers = {}
    for direction, labels in (("forward", recipe.reactant_labels), ("reverse", recipe.product_labels)):
        electronic = (ts["energy_hartree"] - sum(by_label[label]["energy_hartree"] for label in labels)) * HARTREE_TO_J_MOL
        zpe = (ts["raw_harmonic_ZPE_J_mol"] - sum(by_label[label]["raw_harmonic_ZPE_J_mol"] for label in labels)) * recipe.frequency_scale_factor
        barriers[direction] = {"electronic_J_mol": electronic, "scaled_raw_harmonic_ZPE_delta_J_mol": zpe,
                               "raw_harmonic_E0_barrier_J_mol": electronic + zpe,
                               "native_arkane_ZPE_scale_factor": recipe.frequency_scale_factor / 1.014,
                               "native_arkane_expected_E0_barrier_J_mol": electronic + zpe / 1.014}
        gauged = shifted[ts["label"]] - sum(shifted[label] for label in labels)
        if not math.isclose(gauged, electronic, rel_tol=0, abs_tol=1e-6):
            raise ValueError("Balanced per-atom energy gauge changed a barrier")
    plain = {"contract": ATOMIC_CONTRACT_ID if any(isinstance(point, CompositeAtomicSpecies) for point in recipe.points) else CONTRACT_ID,
             "recipe": asdict(recipe), "points": points,
             "balanced_electronic_gauge_hartree_per_atom": gauge, "gauged_electronic_energies_J_mol": shifted,
             "diagnostic_barriers": barriers}
    return {**plain, "recipe_id": "sha256:" + hashlib.sha256(json.dumps(plain, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest(),
            "status": "validated_isolated_composite_recipe", "zero_point_energy_in_electronic_sources": False,
            "conditional_model_execution_allowed": True, "network_reinjection_allowed": False,
            "route_verification_claim_allowed": False, "physical_rate_accuracy_claim_allowed": False,
            "limitations": ["Source/recipe consistency is not a route or physical-accuracy certificate; supplied route artifacts remain independent evidence.",
                            "Gas single-surface RRHO/TST only; no solvent, pressure falloff, spin crossing or electronic-state partition correction is inferred.",
                            "A positive soft vibration is not proof of harmonic validity or a hindered rotor; no frequency floor or mode deletion is applied.",
                            "Explicit symmetry, optical counts, spin and channel assumptions are declarations, not automatic calibration.",
                            "Raw harmonic ZPE diagnostics do not equal native Arkane E0: installed Arkane statmech applies frequencyScaleFactor/1.014 to ZPE; native execution must confirm its convention.",
                            "No experimental fitting, confidence interval, lifetime, absolute formation thermochemistry or production method promotion is supplied."]}


def reopen_composite_recipe_assessment(path: Path) -> tuple[CompositeRateRecipe, dict]:
    """Deserialize a staged recipe and repeat native admission, never trust flags.

    Historical v1 records and the explicit atomic v2 variant remain distinct.
    No scientific binding is refreshed or silently rewritten by this reader.
    """
    def unique(pairs):
        value = {}
        for key, item in pairs:
            if key in value:
                raise ValueError("Duplicate composite recipe JSON key")
            value[key] = item
        return value
    try:
        observed = json.loads(Path(path).read_text(), object_pairs_hook=unique)
        record = dict(observed["recipe"])
        points = []
        for raw in record["points"]:
            raw = dict(raw)
            atomic = "atomic_energy_policy" in raw
            fields = ("energy_input", "energy_output") if atomic else (
                "frequency_input", "frequency_output", "hessian", "energy_input", "energy_output")
            for field in fields:
                raw[field] = ArtifactBinding(**raw[field])
            points.append((CompositeAtomicSpecies if atomic else CompositeStationaryPoint)(**raw))
        record["points"] = tuple(points)
        record["route_evidence"] = tuple(ArtifactBinding(**item) for item in record["route_evidence"])
        for field in ("reactant_labels", "product_labels", "temperatures_K"):
            record[field] = tuple(record[field])
        recipe = CompositeRateRecipe(**record)
        current = validate_composite_rate_recipe(recipe)
        # JSON arrays and dataclass tuples are semantically the same wire type.
        wire = json.loads(json.dumps(current, allow_nan=False))
        if observed != wire:
            raise ValueError("Stored composite assessment differs from revalidated native evidence")
        return recipe, current
    except (OSError, TypeError, KeyError, json.JSONDecodeError) as error:
        raise ValueError("Stored composite recipe assessment is invalid") from error


def _write_new(path: Path, text: str):
    with path.open("x") as handle:
        handle.write(text)


def write_isolated_arkane_recipe(recipe: CompositeRateRecipe, workdir: Path) -> dict:
    """Stage exact native bytes plus separate scalar energies; never execute.

    Validation precedes directory creation. The destination must not exist.
    Failed partial output is preserved, never silently overwritten or deleted.
    The existing diatomic rotor adapter is reused without changing its policy.
    """
    assessment = validate_composite_rate_recipe(recipe)
    workdir = Path(workdir).absolute()
    if workdir.exists() or workdir.is_symlink():
        raise ValueError("Isolated composite recipe requires a fresh nonexistent workspace")
    workdir.mkdir(parents=True, exist_ok=False)
    species_dir = workdir / "species"
    species_dir.mkdir()
    by_label, records, files = {point.label: point for point in recipe.points}, [], {}
    for index, observed in enumerate(assessment["points"], 1):
        point = by_label[observed["label"]]
        artifact_dir = species_dir / f"point-{index:03d}-artifacts"
        artifact_dir.mkdir()
        staged_output, staged_hessian = artifact_dir / "calculation.out", artifact_dir / "calculation.hess"
        atomic = isinstance(point, CompositeAtomicSpecies)
        sources = ((point.energy_output, staged_output),) if atomic else ((point.frequency_output, staged_output), (point.hessian, staged_hessian))
        for binding, destination in sources:
            data = reopen(binding).read_bytes()
            if hashlib.sha256(data).hexdigest() != binding.sha256:
                raise ValueError("Composite source changed during staging")
            with destination.open("xb") as handle:
                handle.write(data)
        file = species_dir / f"point-{index:03d}.py"
        files[point.label] = file
        energy = assessment["gauged_electronic_energies_J_mol"][point.label]
        _write_new(file, f"energy = ({energy!r}, 'J/mol', 'e_electronic')\n"
                   + f"geometry = Log({str(staged_output.resolve())!r})\nfrequencies = Log({str(staged_output.resolve())!r})\n"
                   + f"spin_multiplicity = {point.multiplicity}\nexternalSymmetry = {point.external_symmetry}\nopticalIsomers = {point.optical_isomers}\n")
        record = {"label": point.label, "kind": point.kind, "formula": observed["formula"],
                        "charge": point.charge, "multiplicity": point.multiplicity,
                        "output": point.energy_output.path if atomic else point.frequency_output.path,
                        "staged_output": str(staged_output.resolve()),
                        "statmech_file": str(file.resolve()),
                        "checksums": {"output": point.energy_output.sha256} if atomic else {"output": point.frequency_output.sha256, "hessian": point.hessian.sha256},
                        "electronic_energy_source": {"input": point.energy_input.record(), "output": point.energy_output.record(),
                                                     "energy_hartree": point.energy_hartree, "zero_point_energy_included": False},
                        "energy_reference": {"policy_id": "balanced-per-atom-electronic-gauge-v1",
                                             "species_atom_count": observed["atom_count"], "zero_point_energy_included": False}}
        if atomic:
            record.update({"evidence_kind": "monatomic_single_point", "atomic_mass_amu": point.atomic_mass_amu,
                           "atomic_energy_policy": point.atomic_energy_policy, "isotope": point.isotope})
        else:
            record.update({"hessian": point.hessian.path, "staged_hessian": str(staged_hessian.resolve())})
        records.append(record)
    def level(method):
        electronic, basis = method.split("/")
        return f"LevelOfTheory(method={electronic!r}, basis={basis!r}, software='orca')"

    model = (level(recipe.vibrational_method) if recipe.electronic_method == recipe.vibrational_method
             else f"CompositeLevelOfTheory(freq={level(recipe.vibrational_method)}, energy={level(recipe.electronic_method)})")
    lines = ["# Isolated diagnostic; native outputs are staged unchanged.", f"modelChemistry = {model}",
             f"frequencyScaleFactor = {recipe.frequency_scale_factor!r}", "useHinderedRotors = False",
             "useAtomCorrections = False", "useBondCorrections = False"]
    for point in recipe.points:
        call = "transitionState" if point.kind == "transition_state" else "species"
        extra = "" if point.kind == "transition_state" else f", structure=SMILES({point.smiles!r})"
        lines.append(f"{call}({point.label!r}, {str(files[point.label].relative_to(workdir))!r}{extra})")
    tunneling = None if recipe.tunneling.casefold() == "none" else recipe.tunneling.capitalize()
    lines.extend([f"reaction({recipe.label!r}, reactants={list(recipe.reactant_labels)!r}, products={list(recipe.product_labels)!r}, "
                  f"transitionState={recipe.transition_state_label!r}, tunneling={tunneling!r})",
                  f"kinetics({recipe.label!r}, Tlist=({list(recipe.temperatures_K)!r}, 'K'))"])
    _write_new(workdir / "input.py", "\n".join(lines) + "\n")
    _write_new(workdir / "composite-rate-recipe.json", json.dumps(assessment, indent=2, sort_keys=True, allow_nan=False) + "\n")
    # The kind/policy fields are the existing statmech adapter's wire format,
    # not permission to promote this isolated diagnostic as a verified library.
    provenance = {"schema_version": 2, "kind": "verified_arkane_route", "arkane_statmech_policy": POLICY_ID,
                  "route": {"label": recipe.label}, "orca_artifacts": records,
                  "composite_recipe": bind_artifact(workdir / "composite-rate-recipe.json").record(),
                  "execution_scope": "isolated_composite_rate_diagnostic", "network_reinjection_allowed": False,
                  "physical_rate_accuracy_claim_allowed": False, "recipe_id": assessment["recipe_id"]}
    _write_new(workdir / "provenance.json", json.dumps(provenance, indent=2, sort_keys=True, allow_nan=False) + "\n")
    validate_composite_rate_recipe(recipe)  # Catch mutation during staging.
    rotor_contract = create_rotor_contract(workdir)
    return {"input": bind_artifact(workdir / "input.py").record(),
            "composite_recipe": bind_artifact(workdir / "composite-rate-recipe.json").record(),
            "provenance": bind_artifact(workdir / "provenance.json").record(),
            "statmech_contract": bind_artifact(rotor_contract).record(), "recipe_id": assessment["recipe_id"],
            "native_jobs_dispatched": 0, "network_reinjection_allowed": False}


__all__ = ["ArtifactBinding", "CompositeStationaryPoint", "CompositeAtomicSpecies", "CompositeRateRecipe", "bind_artifact",
           "bind_stationary_point", "bind_atomic_species", "reopen", "validate_composite_rate_recipe", "write_isolated_arkane_recipe",
           "reopen_composite_recipe_assessment"]
