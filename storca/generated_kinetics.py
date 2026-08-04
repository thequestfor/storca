"""Validated local RMG libraries for calculated thermal and photolysis routes."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess


def _checksum(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _adjacency_list(smiles: str, *, rmg_env: str | None = "rmg_env") -> str:
    """Use RMG itself to generate the dictionary representation it will load."""
    script = "from rmgpy.molecule import Molecule; import sys; print(Molecule().from_smiles(sys.argv[1]).to_adjacency_list())"
    command = (["conda", "run", "-n", rmg_env, "python", "-c", script, smiles]
               if rmg_env else ["python", "-c", script, smiles])
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Could not generate RMG adjacency list: {result.stderr.strip()}")
    return result.stdout.strip()


def write_photolysis_library(
    workdir: Path, *, route_id: str, reactant_label: str, reactant_smiles: str,
    products: list[tuple[str, str]], photolysis_rate_constant_s_1: float,
    photolysis_evidence: Path, rmg_env: str | None = "rmg_env",
) -> dict:
    """Write an RMG first-order photon-source library from a validated J rate."""
    if photolysis_rate_constant_s_1 <= 0:
        raise ValueError("Photolysis rate constant must be positive")
    library = Path(workdir) / "kinetics"
    library.mkdir(parents=True, exist_ok=True)
    labels = [(reactant_label, reactant_smiles), *products]
    if len({label for label, _ in labels}) != len(labels):
        raise ValueError("Generated library species labels must be unique")
    # A photon source has no thermal reverse process.  Writing ``=>`` rather
    # than a reversible Chemkin arrow prevents RMG from regenerating the
    # chromophore through an invented detailed-balance reverse reaction.
    reaction_label = f"{reactant_label} => {' + '.join(label for label, _ in products)}"
    reaction_file = library / "reactions.py"
    dictionary_file = library / "dictionary.txt"
    reaction_file.write_text(f'''#!/usr/bin/env python
# encoding: utf-8
name = "storca-photolysis-{route_id}"
shortDesc = "STORCA generated photon source reaction"
longDesc = u"""The Arrhenius A is a declared photolysis rate J, not thermal kinetics.
Evidence: {Path(photolysis_evidence).resolve()}
"""

entry(
    index = 1,
    label = "{reaction_label}",
    kinetics = Arrhenius(A = ({photolysis_rate_constant_s_1!r}, 's^-1'), n = 0, Ea = (0, 'kcal/mol'), T0 = (1, 'K')),
    shortDesc = u"Calculated spectral photolysis source",
    longDesc = u"""Calculated from retained spectrum, cross section, quantum yield, and transmission evidence.""",
)
''')
    dictionary_file.write_text("\n\n".join(f"{label}\n{_adjacency_list(smiles, rmg_env=rmg_env)}" for label, smiles in labels) + "\n")
    manifest = {
        "kind": "storca_generated_kinetics", "source": "sunlight_photolysis", "route_id": route_id,
        "reactants": [reactant_label], "products": [label for label, _ in products], "kinetics_model": "Arrhenius",
        "photolysis_rate_constant_s-1": photolysis_rate_constant_s_1, "thermal_interpretation": "not_a_thermal_arrhenius_rate",
        "photolysis_evidence": str(Path(photolysis_evidence).resolve()), "library": str(library.resolve()),
        "checksums": {"reactions.py": _checksum(reaction_file), "dictionary.txt": _checksum(dictionary_file)},
        "verification_status": "verified_spectral_rate",
    }
    manifest_file = library / "generated-kinetics.json"
    manifest_file.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {**manifest, "manifest": str(manifest_file)}


def validate_generated_library(path: Path, *, temperature_K: float | None = None,
                               pressure_bar: float | None = None) -> dict:
    """Validate a generated-library manifest before it is admitted to RMG."""
    path = Path(path)
    manifest_path = path / "generated-kinetics.json"
    if not manifest_path.is_file():
        raise ValueError(f"Generated kinetics library lacks manifest: {path}")
    manifest = json.loads(manifest_path.read_text())
    if not manifest.get("verification_status", "").startswith("verified"):
        raise ValueError("Generated kinetics library is not verified")
    for name, digest in manifest.get("checksums", {}).items():
        if _checksum(path / name) != digest:
            raise ValueError(f"Generated kinetics checksum mismatch: {name}")
    if manifest.get("source") == "arkane_pdep":
        minimum, maximum = manifest["temperature_range_K"]
        if temperature_K is not None and not minimum <= temperature_K <= maximum:
            raise ValueError("Requested temperature lies outside Arkane kinetics grid")
        minimum, maximum = manifest["pressure_range_bar"]
        if pressure_bar is not None and not minimum <= pressure_bar <= maximum:
            raise ValueError("Requested pressure lies outside Arkane kinetics grid")
    if manifest.get("replaces_reaction") or manifest.get("purpose") == "kinetic_repair":
        validation = validate_rate_replacement_evidence(
            manifest, temperature_K=temperature_K, pressure_bar=pressure_bar,
        )
        if validation["status"] != "passed":
            raise ValueError(
                "Generated replacement kinetics failed physical validation: "
                + ", ".join(validation["blocking_reasons"])
            )
    return manifest


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _rate_record(evidence: dict) -> dict:
    """Resolve the condition-specific rate record without evaluating a model."""
    record = evidence.get("condition_rate") or evidence.get("rate_constant") or {}
    if isinstance(record, (int, float)):
        return {"value": record, "units": evidence.get("rate_units")}
    return record if isinstance(record, dict) else {}


def validate_rate_replacement_evidence(
    evidence: dict,
    *,
    temperature_K: float | None = None,
    pressure_bar: float | None = None,
    collision_tolerance: float = 1e-9,
) -> dict:
    """Validate a condition-specific verified rate before model replacement.

    This function does not evaluate Arrhenius, Chebyshev, or capture models.
    Their producing workflow must retain the evaluated rate and its physical
    ceiling at the immutable condition.  Comparing unlike or absent units is
    prohibited rather than silently converted.
    """
    blockers: list[str] = []
    verification_status = str(evidence.get("verification_status") or "")
    if not verification_status.startswith("verified"):
        blockers.append("rate_source_not_verified")
    route_id = evidence.get("route_id")
    if not route_id:
        blockers.append("replacement_route_id_missing")
    if not (evidence.get("reaction_equation") or evidence.get("replaces_reaction")):
        blockers.append("replacement_reaction_equation_missing")
    try:
        molecularity = int(evidence.get("molecularity"))
    except (TypeError, ValueError):
        molecularity = 0
        blockers.append("reaction_molecularity_missing")
    rate = _rate_record(evidence)
    value = _finite(rate.get("value"))
    units = str(rate.get("units") or "").strip()
    if value is None or value < 0.0:
        blockers.append("condition_rate_missing_or_invalid")
    if not units:
        blockers.append("condition_rate_units_missing")
    rate_temperature = _finite(rate.get("temperature_K"))
    rate_pressure = _finite(rate.get("pressure_bar"))
    if temperature_K is not None and (rate_temperature is None or
                                      not math.isclose(rate_temperature, temperature_K, rel_tol=1e-9, abs_tol=1e-9)):
        blockers.append("condition_rate_temperature_mismatch")
    if pressure_bar is not None and (rate_pressure is None or
                                    not math.isclose(rate_pressure, pressure_bar, rel_tol=1e-9, abs_tol=1e-12)):
        blockers.append("condition_rate_pressure_mismatch")

    collision = evidence.get("collision_limit") or {}
    collision_value = _finite(collision.get("value"))
    collision_units = str(collision.get("units") or "").strip()
    ratio = None
    if molecularity == 2:
        if collision_value is None or collision_value <= 0.0 or not collision.get("source"):
            blockers.append("bimolecular_collision_limit_missing")
        elif collision_units != units:
            blockers.append("collision_limit_units_do_not_match_rate")
        elif value is not None:
            ratio = value / collision_value
            if ratio > 1.0 + collision_tolerance:
                blockers.append("replacement_rate_exceeds_collision_limit")
    elif molecularity == 1:
        # A gas collision ceiling is not the relevant limit for a first-order
        # elementary step; endpoint/TS/capture validation is still mandatory.
        pass
    elif molecularity > 2:
        blockers.append("higher_order_replacement_rate_unsupported")

    return {
        "status": "passed" if not blockers else "failed",
        "route_id": route_id,
        "verification_status": verification_status or None,
        "molecularity": molecularity or None,
        "condition_rate": {"value": value, "units": units or None,
                           "temperature_K": rate_temperature, "pressure_bar": rate_pressure},
        "collision_limit": ({"value": collision_value, "units": collision_units or None,
                             "source": collision.get("source"), "rate_to_limit_ratio": ratio}
                            if molecularity == 2 else {"applicable": False}),
        "blocking_reasons": blockers,
    }
