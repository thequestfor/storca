"""ORCA TD-DFT absorption screening with explicit photolysis boundaries."""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import numpy as np

from src.molecule_tools import sanitize_smiles, smiles_to_xyz
from src.orca_runner import run_orca

from .runs import write_metadata
from .workflow import run_optimization_and_frequency


EV_NM = 1239.841984


def _validated_transitions(transitions: list[dict]) -> list[tuple[float, float]]:
    """Return finite ``(wavelength_nm, oscillator_strength)`` pairs."""
    if not transitions:
        raise ValueError("At least one electronic transition is required")
    validated: list[tuple[float, float]] = []
    for index, state in enumerate(transitions, start=1):
        try:
            wavelength = float(state["wavelength_nm"])
            oscillator_strength = float(state["oscillator_strength"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"Transition {index} lacks numeric wavelength and oscillator strength") from error
        if not math.isfinite(wavelength) or wavelength <= 0:
            raise ValueError(f"Transition {index} has an invalid wavelength")
        if not math.isfinite(oscillator_strength) or oscillator_strength < 0:
            raise ValueError(f"Transition {index} has an invalid oscillator strength")
        validated.append((wavelength, oscillator_strength))
    return validated


def oscillator_strength_cross_sections(transitions: list[dict], *, minimum_nm: float = 200.0,
                                       maximum_nm: float = 800.0, step_nm: float = 1.0,
                                       fwhm_nm: float = 20.0, wavelengths_nm: list[float] | None = None) -> list[dict]:
    """Convert oscillator strengths to an approximate absolute cross section.

    The integrated oscillator-strength sum rule fixes the area of each band;
    the Gaussian width remains a declared modelling assumption.  Values are
    appropriate for a screening photon-absorption calculation, not a
    substitute for a measured condensed-phase UV–Vis spectrum.
    """
    if (not all(math.isfinite(value) for value in (minimum_nm, maximum_nm, step_nm, fwhm_nm))
            or minimum_nm <= 0 or maximum_nm <= minimum_nm or step_nm <= 0 or fwhm_nm <= 0):
        raise ValueError("Invalid cross-section grid or broadening")
    # Integral sigma(nu) dnu = e^2/(4 epsilon_0 m_e c) f in SI.  This is
    # equivalent to pi r_e c f; the former expression used an extra 4*pi.
    elementary_charge, electron_mass, epsilon_0, speed_of_light = 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12, 299792458.0
    area_m2_hz = elementary_charge ** 2 / (4 * epsilon_0 * electron_mass * speed_of_light)
    states = _validated_transitions(transitions)
    points = []
    grid = wavelengths_nm
    supplied_grid = wavelengths_nm is not None
    if grid is None:
        grid = []
        wavelength = minimum_nm
        while wavelength <= maximum_nm + step_nm * 1e-8:
            grid.append(wavelength)
            wavelength += step_nm
    if not grid:
        raise ValueError("Cross-section wavelength grid cannot be empty")
    try:
        grid = [float(wavelength) for wavelength in grid]
    except (TypeError, ValueError) as error:
        raise ValueError("Cross-section wavelength grid must be numeric") from error
    if not all(math.isfinite(wavelength) and wavelength > 0 for wavelength in grid):
        raise ValueError("Cross-section wavelengths must be finite and greater than zero")
    for wavelength in grid:
        frequency = speed_of_light / (wavelength * 1e-9)
        sigma_m2 = 0.0
        for state_wavelength, oscillator_strength in states:
            center_frequency = speed_of_light / (state_wavelength * 1e-9)
            # Convert the declared wavelength FWHM locally to frequency FWHM.
            frequency_fwhm = speed_of_light * (fwhm_nm * 1e-9) / (state_wavelength * 1e-9) ** 2
            stddev = frequency_fwhm / (2 * math.sqrt(2 * math.log(2)))
            gaussian_hz_inverse = math.exp(-0.5 * ((frequency - center_frequency) / stddev) ** 2) / (stddev * math.sqrt(2 * math.pi))
            sigma_m2 += area_m2_hz * oscillator_strength * gaussian_hz_inverse
        points.append({"wavelength_nm": wavelength if supplied_grid else round(wavelength, 6),
                       "absorption_cross_section_cm2_molecule": sigma_m2 * 1e4})
    return points


def write_tddft_input(xyz_file: Path, *, charge: int = 0, multiplicity: int = 1,
                      ncores: int = 1, nroots: int = 20,
                      method_keywords: list[str] | None = None) -> Path:
    """Write a self-contained ORCA TD-DFT vertical-excitation input."""
    if nroots < 1 or ncores < 1 or multiplicity < 1:
        raise ValueError("nroots, ncores, and multiplicity must be at least one")
    xyz_file = Path(xyz_file)
    if not xyz_file.is_file():
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")
    # ORCA enables linear-response TD-DFT through the %tddft block.  `TDDFT`
    # is not a simple-input keyword in ORCA 6 and causes an input error.
    # Near-UV photodissociation commonly contains diffuse/Rydberg character.
    # The former B3LYP/def2-SVP default placed H2O2's first bright state near
    # 205 nm and supplied no representation of its weak long-wavelength tail.
    # A range-separated functional plus an augmented basis is a still-bounded
    # ORCA TD-DFT screen, but is materially more appropriate for this task.
    keywords = [item for item in (method_keywords or ["CAM-B3LYP", "aug-cc-pVDZ"])
                if item.upper() != "TDDFT"]
    lines = ["! " + " ".join([*keywords, "TightSCF"]), f"%pal\n  nprocs {ncores}\nend",
             f"%tddft\n  nroots {nroots}\n  tda false\nend", f"* xyz {charge} {multiplicity}"]
    raw = xyz_file.read_text().splitlines()
    try:
        int(raw[0].strip())
        raw = raw[2:]
    except (ValueError, IndexError):
        pass
    lines.extend(line.strip() for line in raw if line.strip())
    lines.append("*")
    output = xyz_file.parent / "tddft.inp"
    output.write_text("\n".join(lines) + "\n")
    return output


def parse_orca_tddft(output_file: Path) -> list[dict]:
    """Parse ORCA absorption-table transitions without relying on column width."""
    lines = Path(output_file).read_text(errors="replace").splitlines()
    start = next((i for i, line in enumerate(lines)
                  if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line.upper()), None)
    if start is None:
        return []
    transitions = []
    for line in lines[start + 1:]:
        if not line.strip():
            if transitions:
                break
            continue
        values = line.split()
        if not values:
            continue
        arrow = values.index("->") if "->" in values else None
        if arrow is not None and len(values) >= arrow + 6:
            state_match = re.match(r"(\d+)-", values[arrow + 1])
            if not state_match:
                continue
            try:
                energy_eV = float(values[arrow + 2].replace("D", "E"))
                energy_cm_1 = float(values[arrow + 3].replace("D", "E"))
                wavelength_nm = float(values[arrow + 4].replace("D", "E"))
                oscillator_strength = float(values[arrow + 5].replace("D", "E"))
            except ValueError:
                continue
            if (not all(math.isfinite(value) for value in (energy_eV, energy_cm_1, wavelength_nm, oscillator_strength))
                    or energy_eV <= 0 or energy_cm_1 <= 0 or wavelength_nm <= 0 or oscillator_strength < 0):
                continue
            transitions.append({"state": int(state_match.group(1)),
                                "initial_state": values[0], "final_state": values[arrow + 1],
                                "oscillator_strength_gauge": "electric_length",
                                "energy_cm-1": energy_cm_1,
                                "energy_eV": energy_eV, "wavelength_nm": wavelength_nm,
                                "oscillator_strength": oscillator_strength})
            continue
        if not re.fullmatch(r"\d+", values[0]):
            continue
        numbers = []
        for value in values[1:]:
            try:
                numbers.append(float(value.replace("D", "E")))
            except ValueError:
                pass
        # ORCA's electric-dipole table convention begins State, Energy (cm-1),
        # wavelength (nm), fosc.  Keep this conservative and reject malformed
        # rows rather than guessing a spectrum.
        if len(numbers) < 3:
            continue
        energy_cm_1, wavelength_nm, oscillator_strength = numbers[:3]
        if (not all(math.isfinite(value) for value in (energy_cm_1, wavelength_nm, oscillator_strength))
                or energy_cm_1 <= 0 or wavelength_nm <= 0 or oscillator_strength < 0):
            continue
        transitions.append({"state": int(values[0]), "oscillator_strength_gauge": "electric_length",
                            "energy_cm-1": energy_cm_1,
                            "energy_eV": EV_NM / wavelength_nm,
                            "wavelength_nm": wavelength_nm,
                            "oscillator_strength": oscillator_strength})
    return transitions


def predicted_absorption_spectrum(transitions: list[dict], *, minimum_nm: float = 200.0,
                                  maximum_nm: float = 800.0, step_nm: float = 1.0,
                                  fwhm_nm: float = 20.0) -> list[dict]:
    """Broaden oscillator strengths into a relative—not absolute—spectrum."""
    if (not all(math.isfinite(value) for value in (minimum_nm, maximum_nm, step_nm, fwhm_nm))
            or minimum_nm <= 0 or maximum_nm <= minimum_nm or step_nm <= 0 or fwhm_nm <= 0):
        raise ValueError("Invalid absorption-spectrum grid or broadening")
    # Use the same frequency-domain, oscillator-sum-rule-preserving line shape
    # as the absolute cross sections, then normalize only the display trace.
    absolute = oscillator_strength_cross_sections(
        transitions, minimum_nm=minimum_nm, maximum_nm=maximum_nm,
        step_nm=step_nm, fwhm_nm=fwhm_nm,
    )
    peak = max(item["absorption_cross_section_cm2_molecule"] for item in absolute)
    return [
        {
            "wavelength_nm": item["wavelength_nm"],
            "relative_absorption": (item["absorption_cross_section_cm2_molecule"] / peak if peak > 0 else 0.0),
        }
        for item in absolute
    ]


def parse_orca_excited_state_model(output_file: Path) -> str:
    """Identify whether ORCA actually ran full LR-TDDFT or TDA-TDDFT."""
    text = Path(output_file).read_text(errors="replace")
    match = re.search(
        r"Tamm-Dancoff approximation\s+\.\.\.\s+(not operative|deactivated|operative|activated)",
        text, re.IGNORECASE,
    )
    if not match:
        return "unresolved"
    return "full_lr_tddft" if match.group(1).lower() in {"not operative", "deactivated"} else "tda_tddft"


def solar_overlap(spectrum: list[dict], sunlight_csv: Path) -> dict:
    """Calculate relative irradiance-weighted overlap on the source grid.

    Predicted absorption is linearly interpolated onto the sunlight grid.  An
    exact floating-point wavelength match is neither physically meaningful nor
    required, and previously caused valid spectra on different grids to appear
    to have no overlap.
    """
    with Path(sunlight_csv).open(newline="") as handle:
        sunlight = list(csv.DictReader(handle))
    required = {"wavelength_nm", "irradiance_W_m2_nm"}
    if not sunlight or required - set(sunlight[0]):
        raise ValueError("Sunlight CSV requires wavelength_nm and irradiance_W_m2_nm")
    if len(spectrum) < 2:
        raise ValueError("Predicted absorption spectrum requires at least two points")
    try:
        absorption_points = sorted(
            (float(item["wavelength_nm"]), float(item["relative_absorption"]))
            for item in spectrum
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Predicted spectrum requires numeric wavelength_nm and relative_absorption") from error
    if not all(math.isfinite(wavelength) and wavelength > 0 and math.isfinite(value) and value >= 0
               for wavelength, value in absorption_points):
        raise ValueError("Predicted spectrum values must be finite and nonnegative")
    absorption_wavelengths = np.asarray([item[0] for item in absorption_points], dtype=float)
    absorption_values = np.asarray([item[1] for item in absorption_points], dtype=float)
    if np.any(np.diff(absorption_wavelengths) <= 0):
        raise ValueError("Predicted spectrum wavelengths must be unique")

    source_points: list[tuple[float, float]] = []
    try:
        for item in sunlight:
            source_points.append((float(item["wavelength_nm"]), float(item["irradiance_W_m2_nm"])))
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Sunlight spectrum contains a nonnumeric wavelength or irradiance") from error
    source_points.sort()
    if not all(math.isfinite(wavelength) and wavelength > 0 and math.isfinite(irradiance) and irradiance >= 0
               for wavelength, irradiance in source_points):
        raise ValueError("Sunlight wavelengths and irradiances must be finite and nonnegative")
    if any(right[0] <= left[0] for left, right in zip(source_points, source_points[1:])):
        raise ValueError("Sunlight wavelengths must be unique")

    contributions = []
    lower, upper = float(absorption_wavelengths[0]), float(absorption_wavelengths[-1])
    for wavelength, irradiance in source_points:
        if lower <= wavelength <= upper:
            interpolated_absorption = float(np.interp(wavelength, absorption_wavelengths, absorption_values))
            contributions.append({"wavelength_nm": wavelength,
                                  "relative_absorption": interpolated_absorption,
                                  "irradiance_W_m2_nm": irradiance,
                                  "relative_overlap_density": interpolated_absorption * irradiance})
    if len(contributions) < 2:
        raise ValueError("Sunlight and predicted absorption ranges need at least two overlapping source points")
    overlap = sum((left["relative_overlap_density"] + right["relative_overlap_density"]) *
                  (right["wavelength_nm"] - left["wavelength_nm"]) / 2
                  for left, right in zip(contributions, contributions[1:]))
    return {"relative_solar_overlap": overlap, "sunlight_spectrum": str(Path(sunlight_csv)),
            "grid_alignment": "linear_absorption_interpolation_to_sunlight_grid",
            "overlap_range_nm": [contributions[0]["wavelength_nm"], contributions[-1]["wavelength_nm"]],
            "wavelength_contributions": contributions}


def run_photostability_screen(smiles: str, run_dir: Path, *, charge: int = 0, multiplicity: int = 1,
                              ncores: int = 1, nroots: int = 20, sunlight_csv: Path | None = None,
                              method_keywords: list[str] | None = None) -> dict:
    """Run a vertical-excitation sunlight-absorption screen and retain all evidence."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    xyz = run_dir / "input.xyz"
    canonical_smiles = sanitize_smiles(smiles)
    smiles_to_xyz(canonical_smiles, xyz)
    result = {"schema_version": 2, "kind": "orca_tddft_photostability_screen", "smiles": canonical_smiles,
              "status": "failed", "assessment": "incomplete", "photolysis_rate_available": False,
              "limitations": ["TD-DFT oscillator strengths predict relative absorption, not an absolute absorption cross section.",
                              "No photolysis lifetime is reported without quantum yield and resolved product evidence.",
                              "Vertical excitation does not verify an excited-state decomposition pathway."],
              "artifacts": {"input_xyz": str(xyz), "tddft_input": str(run_dir / "tddft.inp"),
                            "tddft_output": str(run_dir / "tddft.out"), "predicted_absorption": str(run_dir / "predicted-absorption.csv")}}
    try:
        ground = run_optimization_and_frequency(xyz, run_dir, charge=charge, multiplicity=multiplicity,
                                                 ncores=ncores, run_frequency=True, method_keywords=method_keywords)
        if not ground["frequency_check"].get("IsMinimum"):
            raise RuntimeError("TD-DFT ground-state geometry is not a verified local minimum")
        tddft_input = write_tddft_input(ground["optimized_xyz"], charge=charge, multiplicity=multiplicity,
                                        ncores=ncores, nroots=nroots, method_keywords=method_keywords)
        output = run_orca(tddft_input)["out"]
        response_model = parse_orca_excited_state_model(output)
        if response_model != "full_lr_tddft":
            raise RuntimeError(
                "ORCA output did not confirm full LR-TDDFT"
                + ("; TDA remained operative" if response_model == "tda_tddft" else "")
            )
        transitions = parse_orca_tddft(output)
        if not transitions:
            raise RuntimeError("ORCA completed without a parseable TD-DFT absorption table")
        spectrum = predicted_absorption_spectrum(transitions)
        absorption_csv = Path(result["artifacts"]["predicted_absorption"])
        with absorption_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["wavelength_nm", "relative_absorption"])
            writer.writeheader(); writer.writerows(spectrum)
        solar_states = [state for state in transitions if 280.0 <= state["wavelength_nm"] <= 800.0]
        bright_solar_states = [state for state in solar_states if state["oscillator_strength"] >= 0.01]
        solar_oscillator_strength_sum = sum(state["oscillator_strength"] for state in solar_states)
        material_solar_absorption = solar_oscillator_strength_sum >= 0.01
        shortest_wavelength = min(state["wavelength_nm"] for state in transitions)
        window_covered = shortest_wavelength <= 280.0
        status = "completed" if material_solar_absorption or window_covered else "incomplete"
        assessment = (
            "solar_absorption_risk_identified" if material_solar_absorption
            else "no_material_oscillator_strength_in_default_solar_window" if window_covered
            else "spectral_window_incomplete"
        )
        result.update(status=status, transitions=transitions, bright_solar_states=bright_solar_states,
                      solar_window_oscillator_strength_sum=solar_oscillator_strength_sum,
                      assessment=assessment, orca_response_model=response_model,
                      ground_state_frequency_check=ground["frequency_check"],
                      conformer_scope="single_rdkit_seed_conformer_ground_state_minimum",
                      spectral_coverage={"window_nm": [280.0, 800.0],
                                         "shortest_computed_transition_nm": shortest_wavelength,
                                         "high_energy_edge_covered": window_covered,
                                         "nroots": nroots,
                                         "material_oscillator_strength_sum_threshold": 0.01},
                      absorption_model={"kind": "relative_display_from_oscillator_strength_cross_sections",
                                        "oscillator_strength_gauge": "electric_length",
                                        "broadening_domain": "frequency_hz",
                                        "fwhm_nm_local_at_band_center": 20.0,
                                        "normalization": "unit_peak_after_sum_rule_area_normalized_bands"})
        if status != "completed":
            result["failure_reason"] = "Requested roots did not reach the 280 nm high-energy edge; a negative bright-state conclusion is unavailable"
        if sunlight_csv:
            result["solar_overlap"] = solar_overlap(spectrum, sunlight_csv)
        result["required_for_photolysis_t95"] = ["absolute absorption cross section or calibrated spectrum", "quantum yield", "resolved photoproduct SMILES"]
    except Exception as error:
        result["failure_reason"] = str(error)
    output_json = run_dir / "photostability-screen.json"
    output_json.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
    write_metadata(run_dir, workflow="orca_tddft_photostability_screen", result_json=str(output_json),
                   charge=charge, multiplicity=multiplicity, nroots=nroots,
                   response_model_requested="full_lr_tddft", oscillator_strength_gauge="electric_length",
                   broadening_model="frequency_domain_sum_rule_gaussian_v1",
                   sunlight_csv=str(sunlight_csv) if sunlight_csv else None)
    return {**result, "result_json": output_json}
