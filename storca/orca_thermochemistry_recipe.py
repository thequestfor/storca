"""Observe ORCA's free-energy recipe; never infer it from program defaults.

This is a provenance/admission check, not an accuracy model. ORCA harmonic
frequencies do not imply purely harmonic vibrational *entropy*: ORCA 6.1 can
apply its Grimme quasi-RRHO treatment. Arkane's own partition functions are a
different recipe and are intentionally not altered by this helper.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re


THERMOCHEMISTRY_RECIPE_CONTRACT = "observed_orca_6_1_thermochemistry_v1"
MANUAL_REFERENCE = "https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html"
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def _one(text: str, pattern: str, label: str) -> str:
    values = re.findall(pattern, text, flags=re.MULTILINE | re.IGNORECASE)
    if len(values) != 1:
        raise ValueError(f"ORCA thermochemistry recipe requires one observed {label}")
    return values[0]


def _positive(value: str, label: str, *, allow_zero: bool = False) -> float:
    result = float(value)
    if not math.isfinite(result) or (result < 0.0 if allow_zero else result <= 0.0):
        raise ValueError(f"Invalid observed ORCA thermochemistry {label}")
    return result


def observed_orca_thermochemistry_recipe(output: Path) -> dict:
    """Bind one normally terminated frequency input/output to observed settings.

    No documented/default value is substituted for missing output evidence.
    Version restriction concerns parser/recipe applicability, not DFT validity.
    Printed scalar matching tolerances below are rounding tolerances, not errors
    or confidence bounds. Species-specific rotational symmetry is retained but
    excluded from the shared recipe identity: reactant and TS may differ.
    """
    output = Path(output)
    input_path = output.with_suffix(".inp")
    if not output.is_file() or not input_path.is_file():
        raise ValueError("ORCA thermochemistry recipe requires exact input/output evidence")
    output_data, input_data = output.read_bytes(), input_path.read_bytes()
    text = output_data.decode("utf-8", errors="replace")
    if text.count("ORCA TERMINATED NORMALLY") != 1 or re.search(
        r"ORCA.*(?:error termination|TERMINATED ABNORMALLY)", text, flags=re.IGNORECASE,
    ):
        raise ValueError("ORCA thermochemistry recipe requires one normal termination")
    version = _one(text, r"^\s*Program Version\s+(\d+\.\d+\.\d+)\b", "program version")
    if not version.startswith("6.1."):
        raise ValueError("Observed thermochemistry recipe parser is restricted to ORCA 6.1")
    quasi = _one(text, r"^\s*Quasi RRHO\s+\.\.\.\s+(True|False)\s*$", "Quasi RRHO selection").casefold() == "true"
    cutoff = _positive(_one(
        text, rf"^\s*Cut-Off Frequency\s+\.\.\.\s+({_NUMBER})\s+cm(?:\^-1|\*\*-1|-1)\s*$",
        "frequency cutoff",
    ), "frequency cutoff", allow_zero=True)
    scale = _positive(_one(
        text, rf"^\s*Scaling factor for frequencies\s*=\s*({_NUMBER})\s+\(already applied!\)\s*$",
        "applied frequency scaling factor",
    ), "frequency scaling factor")
    pressure = _positive(_one(
        text, rf"^\s*Pressure\s+\.\.\.\s+({_NUMBER})\s+atm\s*$", "thermochemistry pressure",
    ), "pressure")
    references = re.findall(
        rf"using a reference frequency of\s+({_NUMBER})\s+cm(?:\^-1|\*\*-1|-1)\b", text,
        flags=re.IGNORECASE,
    )
    if quasi:
        if len(references) != 1 or not re.search(
            r"Vibrational entropy computed according to the QRRHO of S\. Grimme", text,
            flags=re.IGNORECASE,
        ):
            raise ValueError("Quasi-RRHO recipe lacks its observed entropy/reference-frequency evidence")
        reference = _positive(references[0], "QRRHO reference frequency")
    else:
        if references or re.search(r"Vibrational entropy computed according to the QRRHO", text, flags=re.IGNORECASE):
            raise ValueError("RRHO selection contradicts observed QRRHO entropy evidence")
        reference = None
    if not re.search(r"equations for an ideal gas", text, flags=re.IGNORECASE):
        raise ValueError("ORCA thermochemistry recipe lacks its ideal-gas statistical-mechanics declaration")

    input_text = re.sub(r"#[^\n]*", "", input_data.decode("utf-8", errors="replace"))
    keywords = " ".join(re.findall(r"(?m)^\s*!([^\n]*)", input_text))
    if not re.search(r"\b(?:Freq|AnFreq|NumFreq)\b", keywords, flags=re.IGNORECASE):
        raise ValueError("Observed recipe requires the exact frequency calculation input")
    blocks = re.findall(r"%freq\b(.*?)\bend\b", input_text, flags=re.IGNORECASE | re.DOTALL)
    if len(blocks) > 1 or input_text.casefold().count("%freq") != len(blocks):
        raise ValueError("ORCA thermochemistry recipe has conflicting or incomplete frequency blocks")
    block = blocks[0] if blocks else ""
    overrides = {}
    declarations = re.findall(
        r"\b(QuasiRRHO|QRRHORefFreq|CutOffFreq|ScalFreq|Pressure|ProjectTR)\s+(\S+)",
        block, flags=re.IGNORECASE,
    )
    for name, value in declarations:
        key = name.casefold()
        if key in overrides:
            raise ValueError(f"Conflicting ORCA thermochemistry override: {name}")
        if key in {"quasirrho", "projecttr"}:
            if value.casefold() not in {"true", "false"}:
                raise ValueError(f"Invalid ORCA thermochemistry override: {name}")
            overrides[key] = value.casefold() == "true"
        else:
            if not re.fullmatch(_NUMBER, value):
                raise ValueError(f"Scalar ORCA thermochemistry recipe cannot admit a {name} grid/override")
            overrides[key] = _positive(value, name, allow_zero=key == "cutofffreq")
    if "quasirrho" in overrides and overrides["quasirrho"] != quasi:
        raise ValueError("Input QuasiRRHO override contradicts observed thermochemistry")
    if overrides.get("projecttr") is False:
        raise ValueError("Unprojected translation/rotation frequencies cannot supply this solution thermochemistry recipe")
    observed_numbers = {
        "scalfreq": (scale, 0.000000000500001),
        "cutofffreq": (cutoff, 0.00500001),
        "pressure": (pressure, 0.00500001),
        "qrrhoreffreq": (reference, 0.05000001),
    }
    for key, (observed, tolerance) in observed_numbers.items():
        # An explicit reference is inactive when QuasiRRHO is false; retain the
        # override but do not pretend it contributes to RRHO free energies.
        if key == "qrrhoreffreq" and not quasi:
            continue
        if key in overrides and not math.isclose(overrides[key], observed, rel_tol=0.0, abs_tol=tolerance):
            raise ValueError(f"Input {key} override contradicts observed thermochemistry")
    recipe = {
        "contract": THERMOCHEMISTRY_RECIPE_CONTRACT,
        "program_version": version,
        "vibrational_entropy_model": "orca_grimme2012_qrrho" if quasi else "rrho",
        "qrrho_reference_frequency_cm_inv": reference,
        "frequency_cutoff_cm_inv": cutoff,
        "frequency_scale_factor": scale,
        "thermochemistry_pressure_atm": pressure,
        "translation_rotation_model": "orca_ideal_gas_statistical_mechanics",
    }
    recipe_id = "sha256:" + hashlib.sha256(json.dumps(
        recipe, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode()).hexdigest()
    symmetry = re.findall(r"Symmetry Number:\s+(\d+)\b", text, flags=re.IGNORECASE)
    return {
        "schema_version": 1, "kind": "observed_orca_thermochemistry_recipe",
        "status": "observed_and_input_consistent", "program_version": version,
        "recipe": recipe, "recipe_id": recipe_id,
        "source_input": {"path": str(input_path), "sha256": hashlib.sha256(input_data).hexdigest()},
        "source_output": {"path": str(output), "sha256": hashlib.sha256(output_data).hexdigest()},
        "explicit_input_overrides": overrides,
        "observed_rotational_symmetry_number": int(symmetry[0]) if len(symmetry) == 1 else None,
        "limitations": [
            "Recipe observation is not a thermochemistry accuracy or applicability calibration.",
            "Quasi-RRHO entropy is not an explicit hindered-rotor or anharmonic calculation.",
            "Ideal-gas translational/rotational statistical mechanics remains an approximation when combined with implicit-solvent energies.",
        ],
        "manual_reference": MANUAL_REFERENCE,
    }


def matched_orca_thermochemistry_recipes(reactant: dict, transition_state: dict) -> dict:
    """Bind like-with-like free energies while leaving species symmetry distinct."""
    if reactant["recipe_id"] != transition_state["recipe_id"]:
        raise ValueError("Reactant and transition-state ORCA thermochemistry recipes differ")
    return {
        "status": "matched_observed_recipe", "recipe_id": reactant["recipe_id"],
        "recipe": reactant["recipe"],
        "reactant": reactant, "transition_state": transition_state,
    }
