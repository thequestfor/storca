"""Peak-level, non-fitting analysis of calculated and reference IR spectra."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from .benchmark import read_spectrum_csv


def _column_name(path: Path) -> str:
    with Path(path).open(newline="") as handle:
        header = next(csv.reader(handle), None)
    if not header or len(header) < 2:
        raise ValueError("Spectrum CSV requires a two-column header")
    return header[1].strip().lower()


def _absorption_signal(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return x values and an absorption-like signal, regardless of display style."""
    x_values, values = read_spectrum_csv(path)
    column_name = _column_name(path)
    if "transmittance" in column_name:
        transmission = values / 100.0 if "percent" in column_name or float(values.max()) > 1.5 else values
        if np.any(transmission <= 0):
            raise ValueError("Transmittance must be greater than zero for Beer-Lambert conversion")
        absorbance = -np.log10(transmission)
        signal = absorbance - np.percentile(absorbance, 5)
    else:
        signal = values - np.percentile(values, 5)
    return x_values, np.clip(signal, 0.0, None)


def detect_ir_peaks(
    path: Path, *, min_prominence: float = 0.05, minimum_separation_cm_1: float = 12.0
) -> list[dict]:
    """Detect prominent absorption bands from a spectrum CSV.

    This intentionally uses a transparent local-extrema rule rather than peak
    fitting.  The returned intensity is normalized only for comparison within
    the same trace.
    """
    if (not all(math.isfinite(value) for value in (min_prominence, minimum_separation_cm_1))
            or not 0 < min_prominence < 1 or minimum_separation_cm_1 <= 0):
        raise ValueError("min_prominence must be between zero and one and separation must be positive")
    x_values, signal = _absorption_signal(path)
    maximum = float(signal.max())
    if maximum <= 0:
        return []
    candidates = np.where((signal[1:-1] >= signal[:-2]) & (signal[1:-1] > signal[2:]))[0] + 1
    candidates = [index for index in candidates if signal[index] >= maximum * min_prominence]
    # Keep the strongest candidate when nearby digitization steps create a
    # cluster of local extrema.
    selected: list[int] = []
    for index in sorted(candidates, key=lambda item: signal[item], reverse=True):
        if all(abs(x_values[index] - x_values[existing]) >= minimum_separation_cm_1 for existing in selected):
            selected.append(index)
    return [
        {
            "wavenumber_cm-1": round(float(x_values[index]), 2),
            "relative_intensity": round(float(signal[index] / maximum), 5),
        }
        for index in sorted(selected, key=lambda item: x_values[item], reverse=True)
    ]


def match_ir_peaks(reference_peaks: list[dict], calculated_peaks: list[dict], *, window_cm_1: float = 40.0) -> dict:
    """Maximum-cardinality, minimum-position-error one-to-one peak match."""
    if not math.isfinite(window_cm_1) or window_cm_1 <= 0:
        raise ValueError("window_cm_1 must be positive")
    references = sorted(enumerate(reference_peaks), key=lambda item: item[1]["wavenumber_cm-1"])
    calculations = sorted(enumerate(calculated_peaks), key=lambda item: item[1]["wavenumber_cm-1"])

    # On a one-dimensional ordered axis an optimal absolute-error assignment
    # has no crossing pairs. Dynamic programming therefore gives maximum
    # cardinality first and minimum total position error second without a
    # heavyweight assignment dependency.
    empty: tuple[int, float, tuple[tuple[int, int], ...]] = (0, 0.0, ())
    table = [[empty for _ in range(len(calculations) + 1)] for _ in range(len(references) + 1)]

    def better(left, right):
        if left[0] != right[0]:
            return left if left[0] > right[0] else right
        if not math.isclose(left[1], right[1], abs_tol=1e-12):
            return left if left[1] < right[1] else right
        return left if left[2] <= right[2] else right

    for ref_position in range(1, len(references) + 1):
        for calc_position in range(1, len(calculations) + 1):
            best = better(table[ref_position - 1][calc_position], table[ref_position][calc_position - 1])
            ref_index, reference = references[ref_position - 1]
            calc_index, calculated = calculations[calc_position - 1]
            error = abs(calculated["wavenumber_cm-1"] - reference["wavenumber_cm-1"])
            if error <= window_cm_1:
                previous = table[ref_position - 1][calc_position - 1]
                matched = (previous[0] + 1, previous[1] + error, previous[2] + ((ref_index, calc_index),))
                best = better(best, matched)
            table[ref_position][calc_position] = best

    used_reference, used_calculated, matches = set(), set(), []
    for ref_index, calc_index in table[-1][-1][2]:
        difference = calculated_peaks[calc_index]["wavenumber_cm-1"] - reference_peaks[ref_index]["wavenumber_cm-1"]
        used_reference.add(ref_index)
        used_calculated.add(calc_index)
        reference, calculated = reference_peaks[ref_index], calculated_peaks[calc_index]
        matches.append({
            "reference_cm-1": reference["wavenumber_cm-1"],
            "calculated_cm-1": calculated["wavenumber_cm-1"],
            "position_error_cm-1": round(float(difference), 2),
            "reference_relative_intensity": reference["relative_intensity"],
            "calculated_relative_intensity": calculated["relative_intensity"],
            "relative_intensity_error": round(calculated["relative_intensity"] - reference["relative_intensity"], 5),
        })
    errors = [abs(match["position_error_cm-1"]) for match in matches]
    return {
        "matches": sorted(matches, key=lambda item: item["reference_cm-1"], reverse=True),
        "unmatched_reference": [peak for index, peak in enumerate(reference_peaks) if index not in used_reference],
        "unmatched_calculated": [peak for index, peak in enumerate(calculated_peaks) if index not in used_calculated],
        "mean_absolute_position_error_cm-1": round(float(np.mean(errors)), 2) if errors else None,
    }


def analyze_ir_spectra(
    calculated: Path,
    reference: Path,
    *,
    min_prominence: float = 0.05,
    minimum_separation_cm_1: float = 12.0,
    match_window_cm_1: float = 40.0,
) -> dict:
    """Create a peak-level report without shifting, rescaling, or fitting spectra."""
    calculated_peaks = detect_ir_peaks(calculated, min_prominence=min_prominence, minimum_separation_cm_1=minimum_separation_cm_1)
    reference_peaks = detect_ir_peaks(reference, min_prominence=min_prominence, minimum_separation_cm_1=minimum_separation_cm_1)
    report = match_ir_peaks(reference_peaks, calculated_peaks, window_cm_1=match_window_cm_1)
    return {
        "calculated_spectrum": str(calculated),
        "reference_spectrum": str(reference),
        "settings": {
            "min_prominence_fraction": min_prominence,
            "minimum_separation_cm-1": minimum_separation_cm_1,
            "match_window_cm-1": match_window_cm_1,
            "fitting_or_frequency_shift_applied": False,
        },
        "calculated_peaks": calculated_peaks,
        "reference_peaks": reference_peaks,
        **report,
    }


def write_peak_analysis(path: Path, report: dict) -> Path:
    Path(path).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return Path(path)


def write_peak_analysis_plot(path: Path, *, calculated: Path, reference: Path, report: dict) -> Path:
    """Write an annotated trace overlay. Labels identify matched reference bands."""
    import matplotlib.pyplot as plt

    calc_x, calc_y = read_spectrum_csv(calculated)
    ref_x, ref_y = read_spectrum_csv(reference)
    figure, axis = plt.subplots(figsize=(11, 5.5))
    axis.plot(ref_x, ref_y, color="#202020", linewidth=1.1, label="reference")
    axis.plot(calc_x, calc_y, color="#be123c", linewidth=1.0, label="calculated")
    for match in report["matches"]:
        position = match["reference_cm-1"]
        axis.axvline(position, color="#2563eb", alpha=0.2, linewidth=0.8)
        axis.annotate(
            f"{position:.0f}", xy=(position, 0.98), xycoords=("data", "axes fraction"),
            rotation=90, ha="right", va="top", fontsize=7, color="#1d4ed8",
        )
    axis.invert_xaxis()
    axis.set_xlabel("Wavenumber (cm$^{-1}$)")
    axis.set_ylabel("Spectrum value")
    axis.set_title("Peak-level IR comparison (no fitted shift)")
    axis.legend(frameon=False)
    axis.grid(alpha=0.16)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return Path(path)
