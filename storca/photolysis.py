"""Explicit photolysis evidence contracts.

Ground-state electronic-structure and RMG calculations do not supply a
photolysis rate.  This module therefore records a named illumination source
and refuses to manufacture a lifetime without absorption and quantum-yield
evidence.
"""

from __future__ import annotations

from pathlib import Path
import json
import csv
import math


ASTM_G173_23_AM15G = {
    "name": "ASTM G173-03 AM1.5 global hemispherical",
    "source_url": "https://www.nrel.gov/grid/solar-resource/spectra-am1.5",
    "geometry": "37 degree sun-facing tilted surface",
    "spectrum_kind": "terrestrial reference spectral irradiance",
    "exposure_schedule": "continuous_reference_exposure",
    "container_transmission": "directly_exposed_gas; no container attenuation model",
    "status": "reference_spectrum_declared",
}


def sunlight_contract(*, exposure_schedule: str = "continuous_reference_exposure",
                      container_transmission: str = "directly_exposed_gas; no container attenuation model") -> dict:
    """Return the immutable named sunlight contract used by ladder stages."""
    return {**ASTM_G173_23_AM15G, "exposure_schedule": exposure_schedule,
            "container_transmission": container_transmission}


def integrate_photolysis_rate(spectrum_csv: Path, absorption_csv: Path, *, quantum_yield: float | None = None,
                              transmission: float = 1.0) -> dict:
    """Integrate a wavelength-resolved photolysis rate in s^-1.

    ``spectrum_csv`` must have ``wavelength_nm`` and
    ``irradiance_W_m2_nm`` columns. ``absorption_csv`` has
    ``wavelength_nm``, ``absorption_cross_section_cm2_molecule`` and an
    optional ``quantum_yield``.  The wavelength grids must match exactly so
    that a provenance error cannot silently interpolate an arbitrary spectrum.
    """
    if not math.isfinite(transmission) or not 0.0 <= transmission <= 1.0:
        raise ValueError("Transmission must be from zero through one")
    with Path(spectrum_csv).open(newline="") as handle:
        spectrum = list(csv.DictReader(handle))
    with Path(absorption_csv).open(newline="") as handle:
        absorption = list(csv.DictReader(handle))
    required_spectrum = {"wavelength_nm", "irradiance_W_m2_nm"}
    required_absorption = {"wavelength_nm", "absorption_cross_section_cm2_molecule"}
    if not spectrum or required_spectrum - set(spectrum[0]):
        raise ValueError("Spectrum CSV requires wavelength_nm and irradiance_W_m2_nm")
    if not absorption or required_absorption - set(absorption[0]):
        raise ValueError("Absorption CSV requires wavelength_nm and absorption_cross_section_cm2_molecule")
    if len(spectrum) < 2 or len(absorption) < 2:
        raise ValueError("Photolysis integration requires at least two spectral grid points")
    try:
        source_points = [(float(row["wavelength_nm"]), float(row["irradiance_W_m2_nm"])) for row in spectrum]
        absorption_points = [(float(row["wavelength_nm"]),
                              float(row["absorption_cross_section_cm2_molecule"]), row) for row in absorption]
    except (TypeError, ValueError) as error:
        raise ValueError("Spectrum and absorption CSV values must be numeric") from error
    if not all(math.isfinite(wavelength) and wavelength > 0 and math.isfinite(irradiance) and irradiance >= 0
               for wavelength, irradiance in source_points):
        raise ValueError("Spectrum wavelengths and irradiances must be finite and nonnegative")
    if not all(math.isfinite(wavelength) and wavelength > 0 and math.isfinite(cross_section) and cross_section >= 0
               for wavelength, cross_section, _ in absorption_points):
        raise ValueError("Absorption wavelengths and cross sections must be finite and nonnegative")
    wavelengths = [item[0] for item in source_points]
    absorption_wavelengths = [item[0] for item in absorption_points]
    if len(set(wavelengths)) != len(wavelengths) or len(set(absorption_wavelengths)) != len(absorption_wavelengths):
        raise ValueError("Spectrum and absorption wavelength grids must not contain duplicates")
    cross_sections = {wavelength: row for wavelength, _, row in absorption_points}
    if set(wavelengths) != set(cross_sections):
        raise ValueError("Spectrum and absorption wavelength grids must match exactly")
    h = 6.62607015e-34
    c = 299792458.0
    contributions = []
    irradiance_by_wavelength = dict(source_points)
    for wavelength_nm in wavelengths:
        absorption_row = cross_sections[wavelength_nm]
        try:
            phi = float(absorption_row.get("quantum_yield", quantum_yield if quantum_yield is not None else "nan"))
        except (TypeError, ValueError) as error:
            raise ValueError("Each quantum yield must be numeric") from error
        if not math.isfinite(phi) or not 0.0 <= phi <= 1.0:
            raise ValueError("Each quantum yield must be from zero through one")
        wavelength_m = wavelength_nm * 1e-9
        photon_flux = irradiance_by_wavelength[wavelength_nm] * wavelength_m / (h * c)
        # cm2 to m2 conversion is 1e-4.
        rate_density = photon_flux * float(absorption_row["absorption_cross_section_cm2_molecule"]) * 1e-4 * phi * transmission
        contributions.append({"wavelength_nm": wavelength_nm, "rate_density_s-1_nm-1": rate_density})
    contributions.sort(key=lambda item: item["wavelength_nm"])
    integral = sum((left["rate_density_s-1_nm-1"] + right["rate_density_s-1_nm-1"]) *
                   (right["wavelength_nm"] - left["wavelength_nm"]) / 2.0
                   for left, right in zip(contributions, contributions[1:]))
    return {"photolysis_rate_constant_s-1": integral, "integration": "trapezoidal",
            "spectrum": str(spectrum_csv), "absorption": str(absorption_csv),
            "container_transmission": transmission, "wavelength_contributions": contributions}


def evaluate_sunlight_photolysis(*, condition: dict, photolysis_evidence: Path | None = None) -> dict:
    """Validate external absorption/quantum-yield evidence for a light stage.

    A future provider can calculate a spectral integral from the retained
    spectrum and molecular photophysics.  Until then, this function makes the
    missing evidence a first-class result instead of treating sunlight as dark.
    """
    result = {
        "kind": "sunlight_photolysis_screen",
        "condition": condition,
        "status": "photochemical_evidence_incomplete",
        "estimated_time_to_retention_seconds": None,
        "required_evidence": [
            "wavelength-resolved absorption cross section for the target",
            "quantum yield or verified excited-state dissociation pathway",
            "spectral integral using the declared sunlight and transmission contracts",
            "declared radical or molecular photoproduct structures",
        ],
    }
    if photolysis_evidence is None:
        return result
    payload = json.loads(Path(photolysis_evidence).read_text())
    if {"spectrum_csv", "absorption_csv"} <= set(payload):
        evidence_dir = Path(photolysis_evidence).resolve().parent
        spectrum_path = Path(payload["spectrum_csv"])
        absorption_path = Path(payload["absorption_csv"])
        if not spectrum_path.is_absolute():
            spectrum_path = evidence_dir / spectrum_path
        if not absorption_path.is_absolute():
            absorption_path = evidence_dir / absorption_path
        integral = integrate_photolysis_rate(
            spectrum_path, absorption_path,
            quantum_yield=payload.get("quantum_yield"),
            transmission=float(payload.get("container_transmission", 1.0)),
        )
        payload = {**payload, **integral,
                   "absorption_cross_section": {"source": str(absorption_path)},
                   "quantum_yield": {"source": str(absorption_path)}}
    required = {"absorption_cross_section", "quantum_yield", "photolysis_rate_constant_s-1", "photoproducts"}
    missing = sorted(required - set(payload))
    if missing:
        result["provided_evidence"] = str(photolysis_evidence)
        result["missing_fields"] = missing
        return result
    rate = float(payload["photolysis_rate_constant_s-1"])
    if not math.isfinite(rate) or rate <= 0:
        result["provided_evidence"] = str(photolysis_evidence)
        result["failure_reason"] = "photolysis_rate_constant_s-1 must be positive"
        return result
    retention = float(condition["retention_fraction"])
    products = payload["photoproducts"]
    if not isinstance(products, list) or not products or any(not isinstance(item, dict) or not {"label", "smiles"} <= set(item) for item in products):
        result["provided_evidence"] = str(photolysis_evidence)
        result["failure_reason"] = "photoproducts must be a non-empty list of {label, smiles} records"
        return result
    result.update(
        status="completed",
        provided_evidence=str(photolysis_evidence),
        estimated_time_to_retention_seconds=-math.log(retention) / rate,
        model="first_order_photolysis_from_user_supplied_spectral_evidence",
        photoproducts=products,
        photolysis_rate_constant_s_1=rate,
    )
    return result
