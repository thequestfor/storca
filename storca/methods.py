"""Named, explicitly provisional harmonic-spectrum method profiles."""

from __future__ import annotations


HARMONIC_METHOD_PROFILES = {
    "b3lyp-def2-svp": {
        "method": "B3LYP",
        "basis": "def2-SVP",
        "orca_keywords": ["B3LYP", "def2-SVP", "NoAutoStart"],
        "harmonic_scale_factor": 0.970,
        "scale_factor_status": "project baseline; validate against the target chemical class",
    },
    "b3lyp-3c": {
        "method": "B3LYP-3c",
        "basis": "composite def2-mSVP",
        "orca_keywords": ["B3LYP-3c"],
        "harmonic_scale_factor": 1.0,
        "scale_factor_status": "unscaled during the initial bakeoff; no project calibration accepted yet",
        "manual_reference": "ORCA Manual 6.1, section 3.6.5: composite method introduced for gas-phase infrared spectra.",
    },
    "r2scan-3c": {
        "method": "r2SCAN-3c",
        "basis": "composite def2-mTZVPP",
        "orca_keywords": ["r2SCAN-3c"],
        "harmonic_scale_factor": 1.0,
        "scale_factor_status": "unscaled during the initial bakeoff; no project calibration accepted yet",
        "manual_reference": "ORCA Manual 6.1, section 3.6.3: robust composite method for thermochemistry, geometries, and noncovalent interactions.",
    },
}

ORCA_IR_BAKEOFF_CANDIDATES = {
    "b3lyp-def2-svp": "Current project baseline.",
    "b3lyp-3c": "ORCA Manual 6.1, section 3.6.5: composite method introduced for gas-phase IR spectra.",
    "r2scan-3c": "ORCA Manual 6.1, section 3.6.3: robust composite method for thermochemistry, geometries, and noncovalent interactions.",
}


def harmonic_method_profile(name: str = "b3lyp-def2-svp") -> dict:
    """Return a copy of a supported profile so metadata cannot mutate defaults."""
    try:
        return {**HARMONIC_METHOD_PROFILES[name], "name": name}
    except KeyError as error:
        choices = ", ".join(sorted(HARMONIC_METHOD_PROFILES))
        raise ValueError(f"Unknown harmonic method profile '{name}'. Available: {choices}") from error


def arkane_model_chemistry(profile: dict) -> str:
    """Return Arkane's model-chemistry spelling for a retained ORCA profile."""
    method = str(profile.get("method") or "").strip()
    basis = str(profile.get("basis") or "").strip()
    if not method:
        raise ValueError("The ORCA method profile does not identify its method")
    # ORCA composite methods carry their basis/dispersion recipe in the method
    # name and Arkane/RMG recognizes that name as one level of theory.
    if method.casefold().endswith("-3c"):
        return method
    if not basis:
        raise ValueError("A non-composite ORCA method needs a basis for Arkane")
    return f"{method}/{basis}"


def resolve_scale_factor(profile_name: str, explicit_scale_factor: float | None) -> tuple[float, str]:
    """Resolve an override while recording whether it came from a named profile."""
    profile = harmonic_method_profile(profile_name)
    if explicit_scale_factor is not None:
        return explicit_scale_factor, "user_override"
    return float(profile["harmonic_scale_factor"]), "method_profile"
