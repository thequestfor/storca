"""Fail-closed preparation and parsing boundary for the external MESS solver.

This module does not estimate statistical-mechanical or collisional inputs.  It
renders only an explicit, source-backed model and parses one named native MESS
temperature/pressure table.  Execution and scientific calibration remain
separate provider responsibilities.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Mapping, Sequence


MESS_INPUT_SCHEMA_VERSION = 1
MESS_RATE_ADMISSION_SCHEMA_VERSION = 1
AVOGADRO_MOL_1 = 6.02214076e23
MESS_BIMOLECULAR_TO_RMG_SI = 1.0e-6 * AVOGADRO_MOL_1
_MESS_UNCERTAINTY_COMPONENTS = {
    "electronic_structure",
    "conformer_sampling",
    "statistical_mechanics",
    "capture_or_transition_state",
    "energy_transfer",
    "collision_transport",
    "master_equation_numerics",
}


def _positive(value: object, label: str, *, allow_zero: bool = False) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    threshold = 0.0 if allow_zero else math.nextafter(0.0, 1.0)
    if not math.isfinite(number) or number < threshold:
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{label} must be finite and {qualifier}")
    return number


def _canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _scientific_identity(value: Any) -> Any:
    """Remove bundle locations while retaining checksums in cache identity."""
    if isinstance(value, Mapping):
        return {
            key: _scientific_identity(item)
            for key, item in value.items()
            if key != "path"
        }
    if isinstance(value, (list, tuple)):
        return [_scientific_identity(item) for item in value]
    return value


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_immutable(path: Path, content: str) -> None:
    encoded = content.encode()
    if path.is_file():
        if path.read_bytes() != encoded:
            raise ValueError(f"Refusing to overwrite mutated retained MESS artifact: {path}")
        return
    path.write_bytes(encoded)


def _float_list(values: Sequence[object], label: str) -> list[float]:
    result = [_positive(item, label) for item in values]
    if not result:
        raise ValueError(f"{label} must not be empty")
    if result != sorted(set(result)):
        raise ValueError(f"{label} must be unique and sorted")
    return result


def _label(value: object, label: str) -> str:
    text = str(value or "").strip()
    if not text or re.search(r"\s", text):
        raise ValueError(f"{label} must be a nonempty token without whitespace")
    return text


def _validate_artifact(item: Mapping[str, Any], label: str) -> Path:
    path = Path(str(item.get("path") or "")).expanduser().resolve()
    expected = str(item.get("sha256") or "")
    if not path.is_file() or not expected or _file_sha256(path) != expected:
        raise ValueError(f"{label} artifact is missing or its SHA-256 does not match")
    return path


def _validate_rrho(value: Mapping[str, Any], label: str) -> dict[str, Any]:
    geometry = []
    for index, atom in enumerate(value.get("geometry_angstrom") or ()):
        if len(atom) != 4 or not str(atom[0]).strip():
            raise ValueError(f"{label} geometry atom {index} is invalid")
        geometry.append([
            str(atom[0]).strip(),
            float(atom[1]), float(atom[2]), float(atom[3]),
        ])
        if not all(math.isfinite(number) for number in geometry[-1][1:]):
            raise ValueError(f"{label} geometry atom {index} is nonfinite")
    if not geometry:
        raise ValueError(f"{label} geometry is missing")
    frequencies = [
        _positive(item, f"{label} frequency")
        for item in value.get("frequencies_cm1") or ()
    ]
    rotor_kind = str(value.get("rotor_kind") or "").strip()
    expected_modes = {
        "atomic": 0,
        "linear": max(0, 3 * len(geometry) - 5),
        "nonlinear": max(0, 3 * len(geometry) - 6),
    }.get(rotor_kind)
    if expected_modes is None:
        raise ValueError(f"{label} rotor_kind must be atomic, linear, or nonlinear")
    if rotor_kind == "atomic" and len(geometry) != 1:
        raise ValueError(f"{label} atomic rotor_kind requires one atom")
    if rotor_kind == "linear" and len(geometry) < 2:
        raise ValueError(f"{label} linear rotor_kind requires at least two atoms")
    if rotor_kind == "nonlinear" and len(geometry) < 3:
        raise ValueError(f"{label} nonlinear rotor_kind requires at least three atoms")
    if len(frequencies) != expected_modes:
        raise ValueError(
            f"{label} has {len(frequencies)} frequencies; {rotor_kind} "
            f"{len(geometry)}-atom RRHO requires {expected_modes}"
        )
    levels = []
    for item in value.get("electronic_levels") or ():
        levels.append({
            "energy_cm1": _positive(
                item.get("energy_cm1"), f"{label} electronic energy", allow_zero=True,
            ),
            "degeneracy": int(_positive(item.get("degeneracy"), f"{label} degeneracy")),
        })
    if not levels:
        raise ValueError(f"{label} electronic levels are missing")
    artifacts = value.get("evidence_artifacts") or ()
    if not artifacts:
        raise ValueError(f"{label} needs source evidence artifacts")
    for index, artifact in enumerate(artifacts):
        _validate_artifact(artifact, f"{label} evidence {index}")
    return {
        "geometry_angstrom": geometry,
        "rotor_kind": rotor_kind,
        "frequencies_cm1": frequencies,
        "symmetry_factor": _positive(value.get("symmetry_factor"), f"{label} symmetry"),
        "zero_energy_kcal_mol": float(value.get("zero_energy_kcal_mol")),
        "electronic_levels": levels,
        "evidence_artifacts": [dict(item) for item in artifacts],
    }


def validate_mess_model_spec(spec: Mapping[str, Any]) -> dict[str, Any]:
    """Validate a narrow one-well/one-channel MESS input without inventing data."""
    value = dict(spec)
    if value.get("schema_version") != MESS_INPUT_SCHEMA_VERSION:
        raise ValueError("Unsupported STORCA MESS input schema")
    if value.get("kind") != "storca_mess_one_well_model":
        raise ValueError("MESS model kind must be storca_mess_one_well_model")
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", str(value.get("request_sha256") or "")):
        raise ValueError("MESS model must be bound to a specialist request SHA-256")
    temperatures = _float_list(value.get("temperatures_K") or (), "MESS temperatures")
    pressures = _float_list(value.get("pressures_bar") or (), "MESS pressures")

    collision = value.get("collision_model") or {}
    if not isinstance(collision, Mapping):
        raise ValueError("MESS collision model is missing")
    relaxation = collision.get("energy_relaxation") or {}
    if relaxation.get("kind") != "Exponential":
        raise ValueError("Only explicit MESS Exponential energy relaxation is supported")
    normalized_relaxation = {
        "kind": "Exponential",
        "factor_cm1": _positive(relaxation.get("factor_cm1"), "energy relaxation factor"),
        "power": float(relaxation.get("power")),
        "exponent_cutoff": _positive(relaxation.get("exponent_cutoff"), "energy relaxation cutoff"),
        "source": str(relaxation.get("source") or "").strip(),
        "applicability": str(relaxation.get("applicability") or "").strip(),
    }
    if not math.isfinite(normalized_relaxation["power"]):
        raise ValueError("energy relaxation power must be finite")
    if not normalized_relaxation["source"] or not normalized_relaxation["applicability"]:
        raise ValueError("energy relaxation source and applicability are required")

    lj = collision.get("lennard_jones") or {}
    if lj.get("kind") != "LennardJones":
        raise ValueError("Only explicit MESS LennardJones collision frequency is supported")
    normalized_lj = {
        "kind": "LennardJones",
        "epsilons_cm1": [_positive(item, "LJ epsilon") for item in lj.get("epsilons_cm1") or ()],
        "sigmas_angstrom": [_positive(item, "LJ sigma") for item in lj.get("sigmas_angstrom") or ()],
        "masses_amu": [_positive(item, "collision mass") for item in lj.get("masses_amu") or ()],
        "bath_gas": str(lj.get("bath_gas") or "").strip(),
        "source": str(lj.get("source") or "").strip(),
        "applicability": str(lj.get("applicability") or "").strip(),
    }
    if not all(len(normalized_lj[key]) == 2 for key in (
        "epsilons_cm1", "sigmas_angstrom", "masses_amu",
    )):
        raise ValueError("MESS Lennard-Jones model needs well and bath-gas pairs")
    if not all(normalized_lj[key] for key in ("bath_gas", "source", "applicability")):
        raise ValueError("MESS bath gas, collision source, and applicability are required")
    collision_artifacts = collision.get("evidence_artifacts") or ()
    if not collision_artifacts:
        raise ValueError("MESS collision model needs evidence artifacts")
    for index, artifact in enumerate(collision_artifacts):
        _validate_artifact(artifact, f"collision evidence {index}")

    well = value.get("well") or {}
    well_label = _label(well.get("label"), "MESS well label")
    well_rrho = _validate_rrho(well.get("rrho") or {}, "MESS well")
    bimolecular = value.get("bimolecular") or {}
    bimolecular_label = _label(bimolecular.get("label"), "MESS bimolecular label")
    fragments = []
    for index, fragment in enumerate(bimolecular.get("fragments") or ()):
        fragments.append({
            "label": _label(fragment.get("label"), f"MESS fragment {index} label"),
            "rrho": _validate_rrho(fragment.get("rrho") or {}, f"MESS fragment {index}"),
        })
    if len(fragments) != 2:
        raise ValueError("Narrow MESS adapter requires exactly two bimolecular fragments")
    barrier = value.get("barrier") or {}
    barrier_label = _label(barrier.get("label"), "MESS barrier label")
    if barrier.get("model_kind") != "external_rotd":
        raise ValueError("Narrow MESS adapter requires external_rotD microcanonical flux evidence")
    rotd_path = _validate_artifact(barrier.get("rotd_file") or {}, "MESS Rotd")
    barrier_frequencies = [
        _positive(item, "MESS barrier frequency")
        for item in barrier.get("frequencies_cm1") or ()
    ]
    if not barrier_frequencies:
        raise ValueError("MESS Rotd barrier needs retained fragment frequencies")
    barrier_levels = barrier.get("electronic_levels") or ()
    if not barrier_levels:
        raise ValueError("MESS barrier electronic levels are missing")
    normalized_levels = [{
        "energy_cm1": _positive(item.get("energy_cm1"), "barrier electronic energy", allow_zero=True),
        "degeneracy": int(_positive(item.get("degeneracy"), "barrier electronic degeneracy")),
    } for item in barrier_levels]

    controls = value.get("calculation_controls") or {}
    normalized = {
        "schema_version": MESS_INPUT_SCHEMA_VERSION,
        "kind": value["kind"],
        "request_sha256": value["request_sha256"],
        "temperatures_K": temperatures,
        "pressures_bar": pressures,
        "requested_channel": str(value.get("requested_channel") or "").strip(),
        "calculation_controls": {
            "energy_step_over_temperature": _positive(
                controls.get("energy_step_over_temperature", 0.5), "MESS energy step ratio",
            ),
            "excess_energy_over_temperature": _positive(
                controls.get("excess_energy_over_temperature", 30.0), "MESS excess energy ratio",
            ),
            "model_energy_limit_kcal_mol": _positive(
                controls.get("model_energy_limit_kcal_mol"), "MESS model energy limit",
            ),
            "well_cutoff": _positive(controls.get("well_cutoff"), "MESS well cutoff"),
            "calculation_method": str(controls.get("calculation_method") or "").strip(),
        },
        "collision_model": {
            "energy_relaxation": normalized_relaxation,
            "lennard_jones": normalized_lj,
            "evidence_artifacts": [dict(item) for item in collision_artifacts],
        },
        "well": {"label": well_label, "rrho": well_rrho},
        "bimolecular": {
            "label": bimolecular_label,
            "fragments": fragments,
            "ground_energy_kcal_mol": float(bimolecular.get("ground_energy_kcal_mol")),
        },
        "barrier": {
            "label": barrier_label,
            "model_kind": "external_rotd",
            "stoichiometry": _label(barrier.get("stoichiometry"), "MESS barrier stoichiometry"),
            "rotd_file": {"path": str(rotd_path), "sha256": _file_sha256(rotd_path)},
            "symmetry_factor": _positive(barrier.get("symmetry_factor"), "MESS barrier symmetry"),
            "frequencies_cm1": barrier_frequencies,
            "zero_energy_kcal_mol": float(barrier.get("zero_energy_kcal_mol")),
            "electronic_levels": normalized_levels,
        },
    }
    if not normalized["calculation_controls"]["calculation_method"]:
        raise ValueError("MESS calculation method is required")
    if normalized["requested_channel"] not in {
        "well_to_bimolecular", "bimolecular_to_well",
    }:
        raise ValueError(
            "MESS requested_channel must be well_to_bimolecular or bimolecular_to_well"
        )
    for label, number in (
        ("well zero energy", well_rrho["zero_energy_kcal_mol"]),
        ("bimolecular ground energy", normalized["bimolecular"]["ground_energy_kcal_mol"]),
        ("barrier zero energy", normalized["barrier"]["zero_energy_kcal_mol"]),
    ):
        if not math.isfinite(number):
            raise ValueError(f"MESS {label} must be finite")
    return {
        **normalized,
        "model_sha256": _canonical_hash(_scientific_identity(normalized)),
    }


def _format_numbers(values: Sequence[float]) -> str:
    return " ".join(f"{value:.12g}" for value in values)


def _render_rrho(rrho: Mapping[str, Any], indent: str) -> list[str]:
    lines = [f"{indent}RRHO", f"{indent}  Geometry[angstrom] {len(rrho['geometry_angstrom'])}"]
    lines.extend(
        f"{indent}    {atom[0]} {atom[1]:.12g} {atom[2]:.12g} {atom[3]:.12g}"
        for atom in rrho["geometry_angstrom"]
    )
    lines.extend([
        f"{indent}  Core RigidRotor",
        f"{indent}    SymmetryFactor {rrho['symmetry_factor']:.12g}",
        f"{indent}  End",
        f"{indent}  Frequencies[1/cm] {len(rrho['frequencies_cm1'])}",
        f"{indent}    {_format_numbers(rrho['frequencies_cm1'])}",
        f"{indent}  ZeroEnergy[kcal/mol] {rrho['zero_energy_kcal_mol']:.12g}",
        f"{indent}  ElectronicLevels[1/cm] {len(rrho['electronic_levels'])}",
    ])
    lines.extend(
        f"{indent}    {item['energy_cm1']:.12g} {item['degeneracy']}"
        for item in rrho["electronic_levels"]
    )
    lines.append(f"{indent}End")
    return lines


def prepare_mess_input(spec: Mapping[str, Any], workdir: Path) -> dict[str, Any]:
    """Write deterministic MESS input and provenance; never execute the solver."""
    model = validate_mess_model_spec(spec)
    root = Path(workdir)
    root.mkdir(parents=True, exist_ok=True)
    rotd_source = Path(model["barrier"]["rotd_file"]["path"])
    rotd_name = "capture-microcanonical.dat"
    rotd_destination = root / rotd_name
    if rotd_destination.is_file():
        if _file_sha256(rotd_destination) != _file_sha256(rotd_source):
            raise ValueError("Refusing to overwrite mutated retained MESS Rotd evidence")
    elif rotd_source != rotd_destination:
        shutil.copy2(rotd_source, rotd_destination)
    controls = model["calculation_controls"]
    relaxation = model["collision_model"]["energy_relaxation"]
    lj = model["collision_model"]["lennard_jones"]
    well = model["well"]
    bimolecular = model["bimolecular"]
    barrier = model["barrier"]
    lines = [
        f"TemperatureList[K] {_format_numbers(model['temperatures_K'])}",
        f"PressureList[bar] {_format_numbers(model['pressures_bar'])}",
        f"EnergyStepOverTemperature {controls['energy_step_over_temperature']:.12g}",
        f"ExcessEnergyOverTemperature {controls['excess_energy_over_temperature']:.12g}",
        f"ModelEnergyLimit[kcal/mol] {controls['model_energy_limit_kcal_mol']:.12g}",
        f"WellCutoff {controls['well_cutoff']:.12g}",
        f"CalculationMethod {controls['calculation_method']}",
        "Model",
        "  EnergyRelaxation",
        "    Exponential",
        f"      Factor[1/cm] {relaxation['factor_cm1']:.12g}",
        f"      Power {relaxation['power']:.12g}",
        f"      ExponentCutoff {relaxation['exponent_cutoff']:.12g}",
        "    End",
        "  CollisionFrequency",
        "    LennardJones",
        f"      Epsilons[1/cm] {_format_numbers(lj['epsilons_cm1'])}",
        f"      Sigmas[angstrom] {_format_numbers(lj['sigmas_angstrom'])}",
        f"      Masses[amu] {_format_numbers(lj['masses_amu'])}",
        "    End",
        f"  Well {well['label']}",
        "    Species",
        *_render_rrho(well["rrho"], "      "),
        "    End",
        "  End",
        f"  Bimolecular {bimolecular['label']}",
    ]
    for fragment in bimolecular["fragments"]:
        lines.extend([
            f"    Fragment {fragment['label']}",
            *_render_rrho(fragment["rrho"], "      "),
            "    End",
        ])
    lines.extend([
        f"    GroundEnergy[kcal/mol] {bimolecular['ground_energy_kcal_mol']:.12g}",
        "  End",
        f"  Barrier {barrier['label']} {well['label']} {bimolecular['label']}",
        "    RRHO",
        f"      Stoichiometry {barrier['stoichiometry']}",
        "      Core Rotd",
        f"        File {rotd_name}",
        f"        SymmetryFactor {barrier['symmetry_factor']:.12g}",
        "      End",
        f"      Frequencies[1/cm] {len(barrier['frequencies_cm1'])}",
        f"        {_format_numbers(barrier['frequencies_cm1'])}",
        f"      ZeroEnergy[kcal/mol] {barrier['zero_energy_kcal_mol']:.12g}",
        f"      ElectronicLevels[1/cm] {len(barrier['electronic_levels'])}",
    ])
    lines.extend(
        f"        {item['energy_cm1']:.12g} {item['degeneracy']}"
        for item in barrier["electronic_levels"]
    )
    lines.extend(["    End", "  End", "End", ""])
    input_path = root / "mess.inp"
    _write_immutable(input_path, "\n".join(lines))
    manifest_body = {
        "schema_version": 1,
        "kind": "storca_prepared_mess_input",
        "status": "prepared_not_executed",
        "request_sha256": model["request_sha256"],
        "model_sha256": model["model_sha256"],
        "solver": {"name": "MESS", "version": "unresolved_until_execution"},
        "input": input_path.name,
        "input_sha256": _file_sha256(input_path),
        "copied_rotd_file": rotd_destination.name,
        "copied_rotd_sha256": _file_sha256(rotd_destination),
        "model_summary": {
            "well_label": well["label"],
            "bimolecular_label": bimolecular["label"],
            "fragment_labels": [item["label"] for item in bimolecular["fragments"]],
            "requested_channel": model["requested_channel"],
            "temperatures_K": model["temperatures_K"],
            "pressures_bar": model["pressures_bar"],
            "pressure_dependence_model": {
                "calculation_method": controls["calculation_method"],
                "bath_gas": lj["bath_gas"],
                "energy_relaxation": relaxation,
                "lennard_jones": lj,
            },
        },
        "claims_permitted": ["mess_input_prepared"],
        "claims_prohibited": ["rate", "falloff", "lifetime"],
    }
    manifest = {
        **manifest_body,
        "manifest_sha256": _canonical_hash(manifest_body),
    }
    manifest_path = root / "mess-input-manifest.json"
    _write_immutable(
        manifest_path, json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return {**manifest, "manifest": str(manifest_path)}


def validate_prepared_mess_input(workdir: Path) -> dict[str, Any]:
    """Validate an immutable prepared bundle before any external execution."""
    root = Path(workdir).resolve()
    manifest_path = root / "mess-input-manifest.json"
    if not manifest_path.is_file():
        raise ValueError("Prepared MESS input manifest is missing")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("kind") != "storca_prepared_mess_input" \
            or manifest.get("status") != "prepared_not_executed":
        raise ValueError("Prepared MESS input manifest has an invalid kind or status")
    unsigned = {key: value for key, value in manifest.items() if key != "manifest_sha256"}
    if manifest.get("manifest_sha256") != _canonical_hash(unsigned):
        raise ValueError("Prepared MESS input manifest content hash does not match")
    resolved: dict[str, Path] = {}
    for key, hash_key in (
        ("input", "input_sha256"),
        ("copied_rotd_file", "copied_rotd_sha256"),
    ):
        relative = Path(str(manifest.get(key) or ""))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Prepared MESS {key} path is unsafe")
        path = root / relative
        if not path.is_file() or _file_sha256(path) != manifest.get(hash_key):
            raise ValueError(f"Prepared MESS {key} is missing or mutated")
        resolved[key] = path
    summary = manifest.get("model_summary") or {}
    _label(summary.get("well_label"), "prepared MESS well label")
    _label(summary.get("bimolecular_label"), "prepared MESS bimolecular label")
    # Bundles prepared before admission support was introduced did not retain
    # fragment labels.  They remain executable/reusable, but downstream rate
    # admission will fail closed until the input is re-prepared from its model.
    if "fragment_labels" in summary:
        fragment_labels = [
            _label(item, "prepared MESS fragment label")
            for item in summary.get("fragment_labels") or ()
        ]
        if len(fragment_labels) != 2 or len(set(fragment_labels)) != 2:
            raise ValueError("Prepared MESS input must retain two unique fragment labels")
    if summary.get("requested_channel") not in {
        "well_to_bimolecular", "bimolecular_to_well",
    }:
        raise ValueError("Prepared MESS requested channel is invalid")
    _float_list(summary.get("temperatures_K") or (), "prepared MESS temperatures")
    _float_list(summary.get("pressures_bar") or (), "prepared MESS pressures")
    return {
        **manifest,
        "manifest": str(manifest_path),
        "resolved_input": str(resolved["input"]),
        "resolved_rotd_file": str(resolved["copied_rotd_file"]),
    }


def _resolve_executable(executable: str | Path) -> Path:
    raw = str(executable)
    candidate = Path(raw).expanduser()
    if candidate.is_absolute() or candidate.parent != Path("."):
        resolved = candidate.resolve()
    else:
        found = shutil.which(raw)
        if not found:
            raise FileNotFoundError(
                "MESS executable was not found. Install MESS separately or put `mess` on PATH."
            )
        resolved = Path(found).resolve()
    if not resolved.is_file() or not resolved.stat().st_mode & 0o111:
        raise ValueError(f"MESS executable is not an executable file: {resolved}")
    return resolved


def _execution_artifacts(root: Path) -> list[dict[str, Any]]:
    result = []
    for name in (
        "mess.out", "mess.log", "mess.aux", "mess.stdout", "mess.stderr",
        "mess-rate-table.json",
    ):
        path = root / name
        if path.is_file():
            result.append({"path": name, "sha256": _file_sha256(path), "size_bytes": path.stat().st_size})
    return result


def _load_retained_execution(
    root: Path, executable_sha256: str, prepared: Mapping[str, Any],
) -> dict[str, Any] | None:
    path = root / "mess-execution-manifest.json"
    if not path.is_file():
        return None
    retained = json.loads(path.read_text())
    unsigned = {key: value for key, value in retained.items() if key != "execution_sha256"}
    if retained.get("execution_sha256") != _canonical_hash(unsigned):
        raise ValueError("Retained MESS execution manifest content hash does not match")
    if retained.get("executable", {}).get("sha256") != executable_sha256:
        raise ValueError("Retained MESS execution used a different executable binary")
    for key, prepared_key in (
        ("request_sha256", "request_sha256"),
        ("model_sha256", "model_sha256"),
        ("prepared_input_sha256", "input_sha256"),
    ):
        if retained.get(key) != prepared.get(prepared_key):
            raise ValueError("Retained MESS execution does not match the prepared model")
    for artifact in retained.get("artifacts") or ():
        relative = Path(str(artifact.get("path") or ""))
        candidate = root / relative
        if relative.is_absolute() or ".." in relative.parts or not candidate.is_file() \
                or _file_sha256(candidate) != artifact.get("sha256"):
            raise ValueError("Retained MESS execution artifact is missing or mutated")
    return {**retained, "checkpoint_reused": True, "manifest": str(path)}


def run_prepared_mess_input(
    workdir: Path,
    *,
    executable: str | Path = "mess",
    timeout_seconds: float = 3600.0,
) -> dict[str, Any]:
    """Run one prepared MESS bundle with a timeout and retain unadmitted output.

    A successful parse is still not a verified rate. Detailed balance, reference
    calibration, uncertainty, and generated-library admission remain downstream.
    """
    prepared = validate_prepared_mess_input(workdir)
    root = Path(workdir).resolve()
    resolved_executable = _resolve_executable(executable)
    executable_sha256 = _file_sha256(resolved_executable)
    retained = _load_retained_execution(root, executable_sha256, prepared)
    if retained is not None:
        return retained
    orphaned = [item for item in _execution_artifacts(root) if item["path"] in {
        "mess.out", "mess.log", "mess.aux",
    }]
    if orphaned:
        raise ValueError("Unmanifested MESS output exists; use a new attempt directory")
    timeout = _positive(timeout_seconds, "MESS timeout")
    command = [str(resolved_executable), Path(prepared["resolved_input"]).name]
    started = time.monotonic()
    timed_out = False
    try:
        process = subprocess.run(
            command, cwd=root, capture_output=True, text=True, timeout=timeout,
        )
        returncode = process.returncode
        stdout = process.stdout
        stderr = process.stderr
    except subprocess.TimeoutExpired as error:
        timed_out = True
        returncode = None
        stdout = error.stdout.decode(errors="replace") if isinstance(error.stdout, bytes) else (error.stdout or "")
        stderr = error.stderr.decode(errors="replace") if isinstance(error.stderr, bytes) else (error.stderr or "")
    duration = time.monotonic() - started
    summary = prepared["model_summary"]
    if summary["requested_channel"] == "well_to_bimolecular":
        reactant = summary["well_label"]
        product = summary["bimolecular_label"]
        native_rate_units = "1/sec"
    else:
        reactant = summary["bimolecular_label"]
        product = summary["well_label"]
        native_rate_units = "cm^3/sec"
    blockers: list[str] = []
    parsed = None
    output_path = root / "mess.out"
    if timed_out:
        blockers.append("mess_execution_timeout")
    elif returncode != 0:
        blockers.append("mess_execution_nonzero_exit")
    elif not output_path.is_file():
        blockers.append("mess_rate_output_missing")
    else:
        output_text = output_path.read_text(errors="replace")
        if "Unimolecular Rate Units: 1/sec;  Bimolecular Rate Units: cm^3/sec" not in output_text:
            blockers.append("mess_native_rate_unit_declaration_missing")
        try:
            parsed = parse_mess_temperature_pressure_table(
                output_text, reactant_label=reactant, product_label=product,
                pressure_units="bar", rate_units=native_rate_units,
            )
        except ValueError as error:
            blockers.append("mess_requested_rate_table_invalid")
            stderr = stderr + ("\n" if stderr else "") + str(error)
            _write_immutable(root / "mess.stderr", stderr)
        if parsed is not None:
            expected_temperatures = [float(item) for item in summary["temperatures_K"]]
            expected_pressures = [float(item) for item in summary["pressures_bar"]]
            if parsed["temperatures_K"] != expected_temperatures:
                blockers.append("mess_temperature_grid_mismatch")
            if parsed["pressures"] != expected_pressures:
                blockers.append("mess_pressure_grid_mismatch")
            parsed_path = root / "mess-rate-table.json"
            _write_immutable(parsed_path, json.dumps(parsed, indent=2, sort_keys=True) + "\n")
    _write_immutable(root / "mess.stdout", stdout)
    _write_immutable(root / "mess.stderr", stderr)
    status = "completed_unadmitted_rate_table" if not blockers else (
        "execution_timeout" if timed_out else "execution_or_output_validation_failed"
    )
    artifacts = _execution_artifacts(root)
    manifest_body = {
        "schema_version": 1,
        "kind": "storca_mess_execution",
        "status": status,
        "blocking_reasons": blockers,
        "request_sha256": prepared["request_sha256"],
        "model_sha256": prepared["model_sha256"],
        "prepared_input_sha256": prepared["input_sha256"],
        "command": command,
        "returncode": returncode,
        "timed_out": timed_out,
        "timeout_seconds": timeout,
        "wall_time_seconds": duration,
        "executable": {
            "name": "MESS", "path": str(resolved_executable),
            "sha256": executable_sha256,
            "version": "not_reported_by_executable",
        },
        "requested_channel": {"reactant": reactant, "product": product},
        "native_rate_units": native_rate_units,
        "parsed_rate_table": "mess-rate-table.json" if parsed is not None else None,
        "artifacts": artifacts,
        "checkpoint_reused": False,
        "claims_permitted": ["mess_execution_completed", "unadmitted_native_rate_table"] if not blockers else ["mess_execution_failure"],
        "claims_prohibited": ["verified_rate", "falloff_accuracy", "lifetime"],
    }
    manifest = {
        **manifest_body,
        "execution_sha256": _canonical_hash(manifest_body),
    }
    manifest_path = root / "mess-execution-manifest.json"
    _write_immutable(manifest_path, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {**manifest, "manifest": str(manifest_path)}


def parse_mess_temperature_pressure_table(
    output: str,
    *,
    reactant_label: str,
    product_label: str,
    pressure_units: str,
    rate_units: str,
) -> dict[str, Any]:
    """Parse one exact native MESS T/P table; ambiguous or incomplete tables fail."""
    reactant = _label(reactant_label, "MESS output reactant")
    product = _label(product_label, "MESS output product")
    target = f"{reactant}->{product}"
    lines = output.splitlines()
    indices = [index for index, line in enumerate(lines) if line.strip() == target]
    if len(indices) != 1:
        raise ValueError(f"Expected exactly one MESS table named {target}; found {len(indices)}")
    index = indices[0] + 1
    while index < len(lines) and not lines[index].strip():
        index += 1
    if index >= len(lines):
        raise ValueError(f"MESS table {target} has no header")
    header = lines[index].split()
    if not header or header[0] != "P\\T":
        raise ValueError(f"MESS table {target} has an invalid P\\T header")
    temperatures = [_positive(item, "MESS output temperature") for item in header[1:]]
    if not temperatures:
        raise ValueError(f"MESS table {target} has no temperatures")
    index += 1
    points = []
    seen_pressures: set[float] = set()
    while index < len(lines):
        tokens = lines[index].split()
        index += 1
        if not tokens:
            if points:
                break
            continue
        if tokens[0] == "O-O":
            break
        if len(tokens) != len(temperatures) + 1:
            raise ValueError(f"MESS table {target} pressure row has the wrong width")
        pressure = _positive(tokens[0], "MESS output pressure")
        if pressure in seen_pressures:
            raise ValueError(f"MESS table {target} repeats pressure {pressure}")
        seen_pressures.add(pressure)
        for temperature, token in zip(temperatures, tokens[1:]):
            if token == "***":
                raise ValueError(
                    f"MESS table {target} has no rate at {temperature} K and {pressure} {pressure_units}"
                )
            points.append({
                "temperature_K": temperature,
                "pressure": pressure,
                "pressure_units": pressure_units,
                "rate": _positive(token, "MESS output rate"),
                "rate_units": rate_units,
            })
    if not points:
        raise ValueError(f"MESS table {target} contains no rates")
    return {
        "schema_version": 1,
        "kind": "storca_parsed_mess_rate_table",
        "reactant_label": reactant,
        "product_label": product,
        "temperatures_K": temperatures,
        "pressures": sorted(seen_pressures),
        "pressure_units": pressure_units,
        "rate_units": rate_units,
        "points": points,
    }


def validate_mess_execution_bundle(workdir: Path) -> dict[str, Any]:
    """Revalidate a retained successful execution without rerunning MESS."""
    root = Path(workdir).resolve()
    prepared = validate_prepared_mess_input(root)
    path = root / "mess-execution-manifest.json"
    if not path.is_file():
        raise ValueError("MESS execution manifest is missing")
    manifest = json.loads(path.read_text())
    unsigned = {key: value for key, value in manifest.items() if key != "execution_sha256"}
    if manifest.get("execution_sha256") != _canonical_hash(unsigned):
        raise ValueError("MESS execution manifest content hash does not match")
    if manifest.get("kind") != "storca_mess_execution" \
            or manifest.get("status") != "completed_unadmitted_rate_table":
        raise ValueError("MESS execution has no complete unadmitted rate table")
    for key, prepared_key in (
        ("request_sha256", "request_sha256"),
        ("model_sha256", "model_sha256"),
        ("prepared_input_sha256", "input_sha256"),
    ):
        if manifest.get(key) != prepared.get(prepared_key):
            raise ValueError("MESS execution does not match the prepared model")
    declared: dict[str, dict[str, Any]] = {}
    for artifact in manifest.get("artifacts") or ():
        relative = Path(str((artifact or {}).get("path") or ""))
        if relative.is_absolute() or ".." in relative.parts or str(relative) in declared:
            raise ValueError("MESS execution artifact path is unsafe or duplicated")
        candidate = root / relative
        if not candidate.is_file() or _file_sha256(candidate) != artifact.get("sha256"):
            raise ValueError("MESS execution artifact is missing or mutated")
        declared[str(relative)] = dict(artifact)
    for required in ("mess.out", "mess-rate-table.json"):
        if required not in declared:
            raise ValueError(f"MESS execution manifest does not retain {required}")
    table = json.loads((root / "mess-rate-table.json").read_text())
    if table.get("kind") != "storca_parsed_mess_rate_table":
        raise ValueError("Retained MESS rate table has an invalid kind")
    return {
        **manifest,
        "manifest": str(path),
        "prepared": prepared,
        "rate_table": table,
    }


def _evidence_artifact_checks(
    artifacts: object, root: Path, blockers: list[str], prefix: str,
) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    if not isinstance(artifacts, (list, tuple)) or not artifacts:
        blockers.append(f"{prefix}_evidence_artifacts_missing")
        return checks
    for raw in artifacts:
        item = dict(raw or {})
        relative = Path(str(item.get("path") or ""))
        path = relative if relative.is_absolute() else (root / relative).resolve()
        expected = str(item.get("sha256") or "")
        actual = _file_sha256(path) if path.is_file() else None
        valid = bool(actual and actual == expected)
        checks.append({
            "path": str(path), "expected_sha256": expected,
            "actual_sha256": actual, "valid": valid,
        })
        if not valid:
            blockers.append(f"{prefix}_evidence_artifact_missing_or_mutated")
    return checks


def _exact_grid_record(
    records: object, temperatures: Sequence[float], *, value_key: str, label: str,
) -> dict[float, dict[str, Any]]:
    if not isinstance(records, (list, tuple)):
        raise ValueError(f"{label} must be a temperature-indexed list")
    indexed: dict[float, dict[str, Any]] = {}
    for raw in records:
        item = dict(raw or {})
        temperature = _positive(item.get("temperature_K"), f"{label} temperature")
        _positive(item.get(value_key), f"{label} {value_key}")
        if temperature in indexed:
            raise ValueError(f"{label} repeats temperature {temperature}")
        indexed[temperature] = item
    if set(indexed) != set(float(item) for item in temperatures):
        raise ValueError(f"{label} does not exactly cover the MESS temperature grid")
    return indexed


def _point_index(table: Mapping[str, Any]) -> dict[tuple[float, float], dict[str, Any]]:
    result: dict[tuple[float, float], dict[str, Any]] = {}
    for raw in table.get("points") or ():
        item = dict(raw)
        key = (float(item["temperature_K"]), float(item["pressure"]))
        if key in result:
            raise ValueError("MESS rate table contains a duplicate T/P point")
        result[key] = item
    return result


def admit_mess_rate_table(
    workdir: Path,
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    """Admit a MESS T/P table only with explicit units and physical evidence.

    The amount basis hidden by native ``cm^3/sec`` output is never inferred.
    For the reversible one-well model, an independently sourced equilibrium
    rate-ratio table and the opposite MESS direction are checked at every grid
    point.  Interpolation is deliberately outside this admission boundary.
    """
    root = Path(workdir).resolve()
    execution = validate_mess_execution_bundle(root)
    value = dict(evidence)
    blockers: list[str] = []
    if value.get("schema_version") != MESS_RATE_ADMISSION_SCHEMA_VERSION:
        blockers.append("mess_admission_schema_unsupported")
    if value.get("kind") != "storca_mess_rate_admission_evidence":
        blockers.append("mess_admission_kind_invalid")
    if value.get("request_sha256") != execution.get("request_sha256"):
        blockers.append("mess_admission_request_hash_mismatch")
    route_id = str(value.get("route_id") or "").strip()
    equation = str(value.get("reaction_equation") or "").strip()
    if not route_id:
        blockers.append("mess_admission_route_id_missing")
    if not equation:
        blockers.append("mess_admission_reaction_equation_missing")

    summary = execution["prepared"]["model_summary"]
    direction = summary["requested_channel"]
    fragment_labels = list(summary["fragment_labels"])
    if direction == "bimolecular_to_well":
        expected_reactants = fragment_labels
        expected_products = [summary["well_label"]]
        molecularity = 2
        si_units = "m^3/(mol*s)"
        ratio_units = "m^3/mol"
    else:
        expected_reactants = [summary["well_label"]]
        expected_products = fragment_labels
        molecularity = 1
        si_units = "s^-1"
        ratio_units = "mol/m^3"
    reactants = list(value.get("reactant_labels") or ())
    products = list(value.get("product_labels") or ())
    reactant_smiles = list(value.get("reactant_smiles") or ())
    product_smiles = list(value.get("product_smiles") or ())
    if reactants != expected_reactants or products != expected_products:
        blockers.append("mess_admission_reaction_direction_or_labels_mismatch")
    if len(reactants) != len(reactant_smiles) or len(products) != len(product_smiles) \
            or not reactant_smiles or not product_smiles:
        blockers.append("mess_admission_species_identity_incomplete")

    table = execution["rate_table"]
    temperatures = [float(item) for item in table.get("temperatures_K") or ()]
    pressures = [float(item) for item in table.get("pressures") or ()]
    condition = dict(value.get("condition") or {})
    try:
        condition_temperature = _positive(condition.get("temperature_K"), "admission temperature")
        condition_pressure = _positive(condition.get("pressure_bar"), "admission pressure")
    except ValueError:
        condition_temperature = condition_pressure = math.nan
        blockers.append("mess_admission_condition_invalid")
    if condition_temperature not in temperatures or condition_pressure not in pressures:
        blockers.append("mess_admission_condition_not_on_computed_grid")

    unit = dict(value.get("bimolecular_unit_interpretation") or {})
    unit_checks = _evidence_artifact_checks(
        unit.get("evidence_artifacts"), root, blockers, "mess_unit_basis",
    )
    if unit.get("status") != "verified":
        blockers.append("mess_bimolecular_unit_basis_not_verified")
    if unit.get("native_units") != "cm^3/sec" \
            or unit.get("interpreted_units") != "cm^3 molecule^-1 s^-1" \
            or unit.get("si_units") != "m^3/(mol*s)":
        blockers.append("mess_bimolecular_unit_basis_unsupported_or_ambiguous")
    try:
        conversion = float(unit.get("conversion_factor_to_si"))
    except (TypeError, ValueError):
        conversion = math.nan
    if not math.isclose(conversion, MESS_BIMOLECULAR_TO_RMG_SI, rel_tol=1e-12):
        blockers.append("mess_bimolecular_unit_conversion_invalid")
    if unit.get("executable_sha256") != execution.get("executable", {}).get("sha256"):
        blockers.append("mess_unit_basis_executable_mismatch")
    if not all(str(unit.get(key) or "").strip() for key in (
        "source", "source_kind", "applicability",
    )):
        blockers.append("mess_unit_basis_provenance_incomplete")

    uncertainty = dict(value.get("uncertainty") or {})
    try:
        lower = float(uncertainty.get("lower_rate_multiplier"))
        upper = float(uncertainty.get("upper_rate_multiplier"))
    except (TypeError, ValueError):
        lower = upper = math.nan
    if not (math.isfinite(lower) and math.isfinite(upper) and 0 < lower <= 1 <= upper):
        blockers.append("mess_rate_uncertainty_bounds_invalid")
    if not all(str(uncertainty.get(key) or "").strip() for key in (
        "source", "source_kind", "applicability",
    )):
        blockers.append("mess_rate_uncertainty_provenance_incomplete")
    components = {str(item) for item in uncertainty.get("uncertainty_components") or ()}
    missing_components = sorted(_MESS_UNCERTAINTY_COMPONENTS - components)
    if missing_components:
        blockers.append("mess_rate_uncertainty_components_incomplete")

    requested_points = _point_index(table)
    converted_points = [{
        "temperature_K": temperature,
        "pressure_bar": pressure,
        "rate": point["rate"] * (conversion if molecularity == 2 else 1.0),
        "rate_units": si_units,
    } for (temperature, pressure), point in sorted(requested_points.items())]

    collision_checks: list[dict[str, Any]] = []
    collision_artifacts: list[dict[str, Any]] = []
    collision_at_condition = {"applicable": False}
    if molecularity == 2:
        collision = dict(value.get("collision_limit") or {})
        collision_artifacts = _evidence_artifact_checks(
            collision.get("evidence_artifacts"), root, blockers, "mess_collision_limit",
        )
        if not all(str(collision.get(key) or "").strip() for key in (
            "source", "source_kind", "applicability",
        )):
            blockers.append("mess_collision_limit_provenance_incomplete")
        try:
            collision_grid = _exact_grid_record(
                collision.get("values"), temperatures,
                value_key="value", label="MESS collision-limit grid",
            )
        except ValueError:
            collision_grid = {}
            blockers.append("mess_collision_limit_grid_invalid")
        for point in converted_points:
            record = collision_grid.get(point["temperature_K"])
            limit = float(record["value"]) if record else math.nan
            units_match = bool(record and record.get("units") == si_units)
            passed = units_match and math.isfinite(limit) and point["rate"] <= limit * (1 + 1e-9)
            collision_checks.append({
                "temperature_K": point["temperature_K"],
                "pressure_bar": point["pressure_bar"],
                "rate": point["rate"], "limit": limit,
                "units": si_units, "passed": passed,
            })
            if not units_match:
                blockers.append("mess_collision_limit_units_mismatch")
            elif not passed:
                blockers.append("mess_rate_exceeds_collision_limit")
        condition_record = collision_grid.get(condition_temperature)
        if condition_record:
            collision_at_condition = {
                "value": float(condition_record["value"]), "units": si_units,
                "source": collision.get("source"),
            }

    thermodynamics = dict(value.get("thermodynamic_consistency") or {})
    thermodynamic_artifacts = _evidence_artifact_checks(
        thermodynamics.get("evidence_artifacts"), root, blockers, "mess_thermodynamic",
    )
    if thermodynamics.get("status") != "verified" or not all(
        str(thermodynamics.get(key) or "").strip()
        for key in ("source", "source_kind", "applicability")
    ):
        blockers.append("mess_thermodynamic_provenance_incomplete")
    try:
        tolerance = float(thermodynamics.get("maximum_relative_error"))
    except (TypeError, ValueError):
        tolerance = math.nan
    if not math.isfinite(tolerance) or not 0 < tolerance <= 0.05:
        blockers.append("mess_detailed_balance_tolerance_invalid")
    try:
        expected_ratios = _exact_grid_record(
            thermodynamics.get("expected_requested_over_reverse_rate_ratios"),
            temperatures, value_key="value", label="MESS equilibrium rate-ratio grid",
        )
    except ValueError:
        expected_ratios = {}
        blockers.append("mess_equilibrium_rate_ratio_grid_invalid")
    if any(item.get("units") != ratio_units for item in expected_ratios.values()):
        blockers.append("mess_equilibrium_rate_ratio_units_mismatch")

    output_text = (root / "mess.out").read_text(errors="replace")
    reverse_native_units = "1/sec" if molecularity == 2 else "cm^3/sec"
    try:
        reverse_table = parse_mess_temperature_pressure_table(
            output_text,
            reactant_label=execution["requested_channel"]["product"],
            product_label=execution["requested_channel"]["reactant"],
            pressure_units="bar", rate_units=reverse_native_units,
        )
        reverse_points = _point_index(reverse_table)
        if set(reverse_points) != set(requested_points):
            blockers.append("mess_reverse_rate_grid_mismatch")
    except ValueError:
        reverse_table = None
        reverse_points = {}
        blockers.append("mess_reverse_rate_table_missing_or_invalid")

    balance_checks: list[dict[str, Any]] = []
    for key, forward in sorted(requested_points.items()):
        reverse = reverse_points.get(key)
        expected = expected_ratios.get(key[0])
        if reverse is None or expected is None:
            continue
        forward_si = float(forward["rate"]) * (conversion if molecularity == 2 else 1.0)
        reverse_si = float(reverse["rate"]) * (1.0 if molecularity == 2 else conversion)
        actual_ratio = forward_si / reverse_si
        expected_ratio = float(expected["value"])
        relative_error = abs(actual_ratio / expected_ratio - 1.0)
        passed = math.isfinite(tolerance) and relative_error <= tolerance
        balance_checks.append({
            "temperature_K": key[0], "pressure_bar": key[1],
            "actual_requested_over_reverse_rate_ratio": actual_ratio,
            "expected_requested_over_reverse_rate_ratio": expected_ratio,
            "ratio_units": ratio_units, "relative_error": relative_error,
            "passed": passed,
        })
        if not passed:
            blockers.append("mess_detailed_balance_mismatch")

    condition_point = next((
        item for item in converted_points
        if item["temperature_K"] == condition_temperature
        and item["pressure_bar"] == condition_pressure
    ), None)
    blockers = list(dict.fromkeys(blockers))
    status = "admitted" if not blockers else "rejected"
    body = {
        "schema_version": MESS_RATE_ADMISSION_SCHEMA_VERSION,
        "kind": "storca_mess_rate_admission",
        "status": status,
        "blocking_reasons": blockers,
        "request_sha256": execution.get("request_sha256"),
        "route_id": route_id or None,
        "reaction_equation": equation or None,
        "reactant_labels": reactants,
        "product_labels": products,
        "reactant_smiles": reactant_smiles,
        "product_smiles": product_smiles,
        "molecularity": molecularity,
        "requested_channel": direction,
        "temperature_range_K": [min(temperatures), max(temperatures)],
        "pressure_range_bar": [min(pressures), max(pressures)],
        "condition_rate": ({
            "value": condition_point["rate"], "units": si_units,
            "temperature_K": condition_temperature, "pressure_bar": condition_pressure,
        } if condition_point else None),
        "converted_rate_table": {
            "units": si_units, "points": converted_points,
        },
        "bimolecular_unit_interpretation": {
            **unit, "evidence_artifact_checks": unit_checks,
        },
        "collision_limit": collision_at_condition,
        "collision_limit_validation": {
            "applicable": molecularity == 2, "checks": collision_checks,
            "evidence_artifact_checks": collision_artifacts,
        },
        "thermodynamic_consistency": {
            "status": "passed" if balance_checks and all(item["passed"] for item in balance_checks) else "failed",
            "maximum_relative_error": tolerance,
            "checks": balance_checks,
            "evidence_artifact_checks": thermodynamic_artifacts,
        },
        "uncertainty": uncertainty,
        "execution_manifest": execution["manifest"],
        "execution_sha256": execution["execution_sha256"],
        "evidence_sha256": _canonical_hash(value),
        "verification_status": (
            "verified_mess_pressure_dependent_rate_table" if not blockers
            else "incomplete_physical_validation"
        ),
        "claims_permitted": (
            ["condition_rate", "pressure_dependent_rate_table"] if not blockers
            else ["rate_table_rejection"]
        ),
        "claims_prohibited": (
            ["lifetime_before_network_repropagation"] if not blockers
            else ["verified_rate", "network_repair", "lifetime"]
        ),
    }
    result = {**body, "admission_sha256": _canonical_hash(body)}
    path = root / "mess-rate-admission.json"
    _write_immutable(path, json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {**result, "manifest": str(path)}


def validate_mess_rate_admission(workdir: Path) -> dict[str, Any]:
    """Revalidate an admitted rate table and all retained evidence hashes."""
    root = Path(workdir).resolve()
    path = root / "mess-rate-admission.json"
    if not path.is_file():
        raise ValueError("MESS rate-admission manifest is missing")
    admission = json.loads(path.read_text())
    unsigned = {key: value for key, value in admission.items() if key != "admission_sha256"}
    if admission.get("admission_sha256") != _canonical_hash(unsigned):
        raise ValueError("MESS rate-admission manifest content hash does not match")
    if admission.get("kind") != "storca_mess_rate_admission" \
            or admission.get("status") != "admitted" \
            or not str(admission.get("verification_status") or "").startswith("verified"):
        raise ValueError("MESS rate table has not passed scientific admission")
    execution = validate_mess_execution_bundle(root)
    if admission.get("execution_sha256") != execution.get("execution_sha256") \
            or admission.get("request_sha256") != execution.get("request_sha256"):
        raise ValueError("MESS admission does not match the retained execution")
    for section in (
        admission.get("bimolecular_unit_interpretation") or {},
        admission.get("collision_limit_validation") or {},
        admission.get("thermodynamic_consistency") or {},
    ):
        for check in section.get("evidence_artifact_checks") or ():
            artifact = Path(str(check.get("path") or ""))
            if not check.get("valid") or not artifact.is_file() \
                    or _file_sha256(artifact) != check.get("actual_sha256"):
                raise ValueError("MESS admission evidence artifact is missing or mutated")
    return {**admission, "manifest": str(path), "execution": execution}


def _safe_library_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip("-") or "route"


def _rmg_command(rmg_env: str | None, script: str, *arguments: str) -> list[str]:
    return (["conda", "run", "-n", rmg_env, "python", "-c", script, *arguments]
            if rmg_env else ["python", "-c", script, *arguments])


def _marked_result(stdout: str, marker: str) -> dict[str, Any]:
    matches = [line[len(marker):] for line in stdout.splitlines() if line.startswith(marker)]
    if not matches:
        raise RuntimeError(f"RMG helper did not emit {marker.rstrip('=')}")
    return json.loads(matches[-1])


def _fit_mess_rate_table_with_rmg(
    admission: Mapping[str, Any], *, rmg_env: str | None,
    fit_model: str, degree_temperature: int | None,
    degree_pressure: int | None,
) -> dict[str, Any]:
    temperatures = sorted({
        float(item["temperature_K"])
        for item in admission["converted_rate_table"]["points"]
    })
    pressures = sorted({
        float(item["pressure_bar"])
        for item in admission["converted_rate_table"]["points"]
    })
    indexed = {
        (float(item["temperature_K"]), float(item["pressure_bar"])): float(item["rate"])
        for item in admission["converted_rate_table"]["points"]
    }
    matrix = [[indexed[(temperature, pressure)] for pressure in pressures]
              for temperature in temperatures]
    if fit_model not in {"Chebyshev", "PDepArrhenius"}:
        raise ValueError("MESS RMG fit_model must be Chebyshev or PDepArrhenius")
    degree_t = degree_temperature if degree_temperature is not None else min(6, len(temperatures) - 1)
    degree_p = degree_pressure if degree_pressure is not None else min(4, len(pressures) - 1)
    if fit_model == "Chebyshev" and not (
        1 <= degree_t < len(temperatures) and 1 <= degree_p < len(pressures)
    ):
        raise ValueError("Chebyshev degrees must be positive and smaller than the MESS grid dimensions")
    marker = "STORCA_MESS_RMG_FIT_JSON="
    script = r'''
import json, math, sys
import numpy as np
import rmgpy
from rmgpy.kinetics import Chebyshev, PDepArrhenius
data = json.loads(sys.argv[1])
temperatures = np.array(data["temperatures_K"], dtype=float)
pressures_pa = np.array(data["pressures_bar"], dtype=float) * 1.0e5
rates = np.array(data["rates"], dtype=float)
if data["fit_model"] == "Chebyshev":
    model = Chebyshev().fit_to_data(
        temperatures, pressures_pa, rates, data["units"],
        int(data["degree_temperature"]), int(data["degree_pressure"]),
        float(temperatures.min()), float(temperatures.max()),
        float(pressures_pa.min()), float(pressures_pa.max()),
    )
else:
    model = PDepArrhenius().fit_to_data(
        temperatures, pressures_pa, rates, data["units"],
    )
evaluated = []
for i, temperature in enumerate(temperatures):
    for j, pressure_pa in enumerate(pressures_pa):
        predicted = model.get_rate_coefficient(float(temperature), float(pressure_pa))
        if not math.isfinite(predicted) or predicted <= 0:
            raise ValueError("Fitted RMG model returned an invalid rate")
        evaluated.append({
            "temperature_K": float(temperature),
            "pressure_bar": float(pressure_pa / 1.0e5),
            "source_rate": float(rates[i, j]),
            "fitted_rate": float(predicted),
            "log10_error": abs(math.log10(predicted / rates[i, j])),
        })
print("STORCA_MESS_RMG_FIT_JSON=" + json.dumps({
    "canonical_expression": repr(model),
    "model_class": model.__class__.__name__,
    "pressure_dependent": bool(model.is_pressure_dependent()),
    "rmg_version": getattr(rmgpy, "__version__", None),
    "grid_evaluation": evaluated,
}, sort_keys=True))
'''
    payload = {
        "temperatures_K": temperatures, "pressures_bar": pressures,
        "rates": matrix, "units": admission["converted_rate_table"]["units"],
        "fit_model": fit_model, "degree_temperature": degree_t,
        "degree_pressure": degree_p,
    }
    command = _rmg_command(rmg_env, script, json.dumps(payload, separators=(",", ":")))
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        raise RuntimeError(
            "Could not fit the admitted MESS grid with RMG: "
            + (process.stderr.strip() or process.stdout.strip())
        )
    result = _marked_result(process.stdout, marker)
    if not result.get("pressure_dependent"):
        raise ValueError("RMG fit did not produce a pressure-dependent model")
    return {
        **result,
        "fit_request": payload,
        "command": command[:-1] + ["<rate-grid-json>"],
    }


def _validate_mess_library_grid_with_rmg(
    library: Path, admission: Mapping[str, Any], *, rmg_env: str | None,
) -> dict[str, Any]:
    marker = "STORCA_MESS_LIBRARY_GRID_JSON="
    script = r'''
import json, math, sys
import rmgpy
from rmgpy.data.kinetics.database import KineticsDatabase
payload = json.loads(sys.argv[3])
database = KineticsDatabase()
database.load_libraries(sys.argv[1], libraries=[sys.argv[2]])
library = database.libraries[sys.argv[2]]
reactions = library.get_library_reactions()
if len(reactions) != 1:
    raise ValueError("MESS repair library must contain exactly one reaction")
reaction = reactions[0]
if not reaction.kinetics.is_pressure_dependent():
    raise ValueError("MESS repair library did not retain pressure-dependent kinetics")
checks = []
for point in payload:
    fitted = reaction.kinetics.get_rate_coefficient(
        float(point["temperature_K"]), float(point["pressure_bar"]) * 1.0e5,
    )
    source = float(point["rate"])
    if not math.isfinite(fitted) or fitted <= 0:
        raise ValueError("Saved MESS library returned an invalid grid rate")
    checks.append({
        "temperature_K": float(point["temperature_K"]),
        "pressure_bar": float(point["pressure_bar"]),
        "source_rate": source, "saved_library_rate": float(fitted),
        "log10_error": abs(math.log10(fitted / source)),
    })
print("STORCA_MESS_LIBRARY_GRID_JSON=" + json.dumps({
    "reaction": str(reaction),
    "reactants": [species.label for species in reaction.reactants],
    "products": [species.label for species in reaction.products],
    "model_class": reaction.kinetics.__class__.__name__,
    "pressure_dependent": True,
    "rmg_version": getattr(rmgpy, "__version__", None),
    "checks": checks,
}, sort_keys=True))
'''
    payload = json.dumps(admission["converted_rate_table"]["points"], separators=(",", ":"))
    command = _rmg_command(
        rmg_env, script, str(library.parent.resolve()), library.name, payload,
    )
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        raise RuntimeError(
            "Saved MESS library failed full-grid RMG validation: "
            + (process.stderr.strip() or process.stdout.strip())
        )
    return {**_marked_result(process.stdout, marker),
            "command": command[:-1] + ["<rate-grid-json>"]}


def write_admitted_mess_rmg_library(
    workdir: Path,
    *,
    output_dir: Path | None = None,
    rmg_env: str | None = "rmg_env",
    fit_model: str = "Chebyshev",
    degree_temperature: int | None = None,
    degree_pressure: int | None = None,
    maximum_log10_error: float = 0.05,
) -> dict[str, Any]:
    """Fit, load, and grid-validate an admitted MESS table as an RMG library."""
    root = Path(workdir).resolve()
    admission = validate_mess_rate_admission(root)
    tolerance = _positive(maximum_log10_error, "MESS RMG fit tolerance")
    if tolerance > 0.3:
        raise ValueError("MESS RMG fit tolerance cannot exceed 0.3 log10 units")
    fit = _fit_mess_rate_table_with_rmg(
        admission, rmg_env=rmg_env, fit_model=fit_model,
        degree_temperature=degree_temperature, degree_pressure=degree_pressure,
    )
    maximum_error = max(float(item["log10_error"]) for item in fit["grid_evaluation"])
    if maximum_error > tolerance:
        raise ValueError(
            f"RMG {fit_model} fit exceeds the declared MESS grid tolerance: "
            f"{maximum_error:.6g} > {tolerance:.6g} log10 units"
        )

    from .generated_kinetics import _adjacency_list, _checksum
    from .arkane_runner import _validate_library_with_rmg

    library = (
        Path(output_dir).resolve() if output_dir is not None else
        root / "verified_RMG_library" / f"storca-verified-{_safe_library_label(admission['route_id'])}"
    )
    library.mkdir(parents=True, exist_ok=True)
    equation = f"{' + '.join(admission['reactant_labels'])} <=> {' + '.join(admission['product_labels'])}"
    reaction_text = (
        "#!/usr/bin/env python\n# encoding: utf-8\n\n"
        "name = 'STORCA admitted MESS pressure-dependent route'\n"
        "shortDesc = 'One evidence-admitted MESS rate table'\n"
        "longDesc = '''Generated only after unit, collision, uncertainty, and detailed-balance checks.'''\n"
        "autoGenerated = True\n\n"
        "entry(\n"
        "    index = 1,\n"
        f"    label = {equation!r},\n"
        f"    kinetics = {fit['canonical_expression']},\n"
        "    shortDesc = 'STORCA MESS pressure-dependent repair',\n"
        "    longDesc = '''Do not detach this model from generated-kinetics.json provenance.''',\n"
        ")\n"
    )
    _write_immutable(library / "reactions.py", reaction_text)
    species = list(zip(
        [*admission["reactant_labels"], *admission["product_labels"]],
        [*admission["reactant_smiles"], *admission["product_smiles"]],
    ))
    seen: dict[str, str] = {}
    for label, smiles in species:
        if label in seen and seen[label] != smiles:
            raise ValueError("MESS library maps one species label to multiple structures")
        seen[label] = smiles
    dictionary_text = "\n\n".join(
        f"{label}\n{_adjacency_list(smiles, rmg_env=rmg_env)}"
        for label, smiles in seen.items()
    ) + "\n"
    _write_immutable(library / "dictionary.txt", dictionary_text)
    condition = admission["condition_rate"]
    loaded = _validate_library_with_rmg(
        library,
        temperature_K=float(condition["temperature_K"]),
        pressure_bar=float(condition["pressure_bar"]),
        pressure_dependent=True,
        rmg_env=rmg_env,
    )
    if loaded["model_class"] != fit["model_class"]:
        raise ValueError("Saved MESS library changed the fitted RMG model class")
    saved_grid = _validate_mess_library_grid_with_rmg(
        library, admission, rmg_env=rmg_env,
    )
    saved_maximum_error = max(float(item["log10_error"]) for item in saved_grid["checks"])
    if saved_maximum_error > tolerance:
        raise ValueError(
            "Saved MESS library exceeds the declared full-grid tolerance: "
            f"{saved_maximum_error:.6g} > {tolerance:.6g} log10 units"
        )
    reactions = library / "reactions.py"
    dictionary = library / "dictionary.txt"
    manifest = {
        "schema_version": 2,
        "kind": "storca_generated_kinetics",
        "source": "mess_pdep",
        "route_id": admission["route_id"],
        "reactants": admission["reactant_labels"],
        "products": admission["product_labels"],
        "reactant_smiles": admission["reactant_smiles"],
        "product_smiles": admission["product_smiles"],
        "reaction_equation": admission["reaction_equation"],
        "replaces_reaction": admission["reaction_equation"],
        "library_reaction_equation": equation,
        "purpose": "kinetic_repair",
        "molecularity": admission["molecularity"],
        "kinetics_model": "MESS master equation with RMG interpolation",
        "rate_model_class": fit["model_class"],
        "temperature_range_K": admission["temperature_range_K"],
        "pressure_range_bar": admission["pressure_range_bar"],
        "condition_rate": condition,
        "collision_limit": admission["collision_limit"],
        "pressure_interpretation": "MESS master-equation condition rate",
        "pressure_dependence_model": admission["execution"]["prepared"]["model_summary"].get(
            "pressure_dependence_model"
        ),
        "rate_uncertainty": admission["uncertainty"],
        "mess_admission": {
            "path": admission["manifest"],
            "sha256": admission["admission_sha256"],
        },
        "rmg_grid_fit_validation": {
            "fit_model": fit_model,
            "maximum_allowed_log10_error": tolerance,
            "maximum_observed_log10_error_before_serialization": maximum_error,
            "maximum_observed_log10_error_saved_library": saved_maximum_error,
            "all_grid_points_passed": True,
            "fit_points": fit["grid_evaluation"],
            "saved_library_points": saved_grid["checks"],
            "rmg_version": fit.get("rmg_version"),
        },
        "rmg_library_validation": loaded,
        "library": str(library),
        "checksums": {
            "reactions.py": _checksum(reactions),
            "dictionary.txt": _checksum(dictionary),
        },
        "verification_status": "verified_mess_pdep",
    }
    from .generated_kinetics import validate_rate_replacement_evidence
    physical = validate_rate_replacement_evidence(
        manifest,
        temperature_K=float(condition["temperature_K"]),
        pressure_bar=float(condition["pressure_bar"]),
    )
    if physical["status"] != "passed":
        raise ValueError(
            "Admitted MESS library failed replacement validation: "
            + ", ".join(physical["blocking_reasons"])
        )
    manifest["physical_validation"] = physical
    manifest_path = library / "generated-kinetics.json"
    _write_immutable(manifest_path, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {**manifest, "manifest": str(manifest_path)}


def finalize_mess_specialist_rate_result(
    workdir: Path,
    specialist_request: Mapping[str, Any],
    *,
    library_dir: Path | None = None,
) -> dict[str, Any]:
    """Package an admitted MESS library for the existing repair/repropagation loop."""
    from .generated_kinetics import validate_generated_library
    from .specialist_rates import (
        validate_specialist_rate_request,
        validate_specialist_rate_result,
    )

    root = Path(workdir).resolve()
    request = dict(specialist_request)
    request_validation = validate_specialist_rate_request(request)
    if request_validation["status"] != "passed":
        raise ValueError(
            "Invalid specialist request: "
            + ", ".join(request_validation["blocking_reasons"])
        )
    admission = validate_mess_rate_admission(root)
    if request.get("request_sha256") != admission.get("request_sha256"):
        raise ValueError("MESS admission does not match the specialist request")
    if "validated_external_falloff_model" not in set(request.get("allowed_rate_theories") or ()):
        raise ValueError("Diagnosed path does not permit an external falloff model")
    for field in ("route_id", "reactant_labels", "product_labels", "reactant_smiles", "product_smiles"):
        if request.get(field) != admission.get(field):
            raise ValueError(f"MESS admission {field} does not match the specialist request")
    condition = request["condition_contract"]
    if not math.isclose(
        float(condition["temperature_K"]), float(admission["condition_rate"]["temperature_K"]),
        rel_tol=1e-12,
    ) or not math.isclose(
        float(condition["pressure_bar"]), float(admission["condition_rate"]["pressure_bar"]),
        rel_tol=1e-12,
    ):
        raise ValueError("MESS admission condition does not match the specialist request")

    library = (
        Path(library_dir).resolve() if library_dir is not None else
        root / "verified_RMG_library" / f"storca-verified-{_safe_library_label(admission['route_id'])}"
    )
    manifest = validate_generated_library(
        library,
        temperature_K=float(condition["temperature_K"]),
        pressure_bar=float(condition["pressure_bar"]),
    )
    if manifest.get("source") != "mess_pdep" or manifest.get("route_id") != request.get("route_id"):
        raise ValueError("Generated library is not the admitted MESS model for this route")
    execution = admission["execution"]
    pressure_model = manifest.get("pressure_dependence_model") or {}
    lj = pressure_model.get("lennard_jones") or {}
    relaxation = pressure_model.get("energy_relaxation") or {}
    admission_path = Path(admission["manifest"])
    execution_path = Path(execution["manifest"])
    library_manifest_path = library / "generated-kinetics.json"
    result = {
        "schema_version": 1,
        "kind": "storca_specialist_rate_result",
        "status": "verified",
        "request_sha256": request["request_sha256"],
        "condition_contract_sha256": request["condition_contract_sha256"],
        "route_id": request["route_id"],
        "rate_theory": "validated_external_falloff_model",
        "provider": {
            "name": "MESS",
            "version": "binary-sha256:" + str(execution["executable"]["sha256"]),
            "executable_sha256": execution["executable"]["sha256"],
            "rmg_interpolation_version": manifest["rmg_grid_fit_validation"].get("rmg_version"),
        },
        "applicability": {
            "status": "applicable",
            "diagnostics": [
                "exact specialist-request and reaction identity match",
                "explicit MESS binary-bound bimolecular amount basis",
                "full-grid forward/reverse detailed-balance check",
                "full-grid collision-ceiling check where bimolecular",
                "full-grid serialized RMG interpolation check",
                "explicit source-backed collisional and uncertainty models",
            ],
        },
        "uncertainty": admission["uncertainty"],
        "pressure_dependence": {
            "pressure_range_bar": manifest["pressure_range_bar"],
            "temperature_range_K": manifest["temperature_range_K"],
            "bath_gas": lj.get("bath_gas"),
            "energy_transfer_model": relaxation,
            "collision_transport_model": lj,
            "master_equation_method": pressure_model.get("calculation_method"),
        },
        "evidence_artifacts": [{
            "path": str(path), "sha256": _file_sha256(path),
        } for path in (admission_path, execution_path, library_manifest_path)],
        "generated_kinetics_library": manifest,
        "claims_permitted": [
            "condition_rate", "pressure_dependent_rate_model", "network_repair_input",
        ],
        "claims_prohibited": [
            "lifetime_before_network_repropagation",
            "applicability_outside_declared_temperature_pressure_and_bath_gas",
        ],
    }
    validation = validate_specialist_rate_result(result, request)
    if validation["status"] != "passed":
        raise ValueError(
            "Final MESS specialist result failed repair-loop validation: "
            + ", ".join(validation["blocking_reasons"])
        )
    result["validation"] = validation
    path = root / "specialist-rate-result.json"
    _write_immutable(path, json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {**result, "result": str(path)}


__all__ = [
    "MESS_INPUT_SCHEMA_VERSION",
    "MESS_RATE_ADMISSION_SCHEMA_VERSION",
    "admit_mess_rate_table",
    "finalize_mess_specialist_rate_result",
    "MESS_BIMOLECULAR_TO_RMG_SI",
    "parse_mess_temperature_pressure_table",
    "prepare_mess_input",
    "run_prepared_mess_input",
    "validate_mess_rate_admission",
    "validate_mess_execution_bundle",
    "validate_prepared_mess_input",
    "validate_mess_model_spec",
    "write_admitted_mess_rmg_library",
]
