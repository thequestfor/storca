"""Numerical comparison of predicted and reference IR spectra."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


def _spectrum_convention(path: Path) -> dict:
    with Path(path).open(newline="") as handle:
        header = next(csv.reader(handle), None)
    if not header or len(header) < 2:
        raise ValueError("Spectrum CSV requires a two-column header")
    axis_header, intensity_header = (value.strip().lower() for value in header[:2])
    axis = "wavelength" if "wavelength" in axis_header else "wavenumber" if "wavenumber" in axis_header else "unspecified"
    if "transmittance" in intensity_header:
        intensity = "transmittance"
    elif "absorbance" in intensity_header:
        intensity = "absorbance"
    elif "relative" in intensity_header:
        intensity = "relative"
    else:
        intensity = "unspecified"
    return {"axis": axis, "intensity": intensity, "headers": header[:2]}


def read_spectrum_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 3 or len(rows[0]) < 2:
        raise ValueError("Spectrum CSV needs a header and at least two numeric rows")
    try:
        data = np.asarray([[float(row[0]), float(row[1])] for row in rows[1:] if row and any(cell.strip() for cell in row)], dtype=float)
    except (IndexError, TypeError, ValueError) as error:
        raise ValueError("Spectrum CSV rows must contain two numeric columns") from error
    if data.ndim != 2 or data.shape[0] < 2 or data.shape[1] != 2 or not np.all(np.isfinite(data)):
        raise ValueError("Spectrum CSV requires at least two finite two-column data rows")
    order = np.argsort(data[:, 0])
    x_values, y_values = data[order, 0], data[order, 1]
    if np.any(np.diff(x_values) <= 0):
        raise ValueError("Spectrum coordinates must be unique")
    return x_values, y_values


def _absorption_like(values: np.ndarray, convention: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if convention == "transmittance":
        signal = float(np.max(values)) - values
    else:
        signal = values - float(np.min(values))
    return np.maximum(signal, 0.0)


def _shape_metrics(
    x_values: np.ndarray, predicted_values: np.ndarray, reference_values: np.ndarray,
    predicted_convention: str, reference_convention: str,
) -> dict:
    predicted = _absorption_like(predicted_values, predicted_convention)
    reference = _absorption_like(reference_values, reference_convention)
    predicted_norm = predicted / max(float(np.max(predicted)), 1e-300)
    reference_norm = reference / max(float(np.max(reference)), 1e-300)
    denominator = float(np.linalg.norm(predicted_norm) * np.linalg.norm(reference_norm))
    cosine = None if denominator <= 0 else float(np.dot(predicted_norm, reference_norm) / denominator)

    def moments(signal: np.ndarray) -> tuple[float | None, float | None]:
        area = float(np.trapezoid(signal, x_values))
        if area <= 0:
            return None, None
        center = float(np.trapezoid(signal * x_values, x_values) / area)
        variance = float(np.trapezoid(signal * np.square(x_values - center), x_values) / area)
        if variance <= 0:
            return 0.0, 0.0
        skew = float(
            np.trapezoid(signal * np.power(x_values - center, 3), x_values)
            / (area * variance ** 1.5)
        )
        return 2.0 * np.sqrt(2.0 * np.log(2.0)) * np.sqrt(variance), skew

    predicted_width, predicted_asymmetry = moments(predicted_norm)
    reference_width, reference_asymmetry = moments(reference_norm)
    regions = []
    for lower, upper in ((400.0, 1800.0), (1800.0, 2800.0), (2800.0, 4000.0)):
        mask = (x_values >= lower) & (x_values <= upper)
        if np.count_nonzero(mask) < 2:
            continue
        predicted_area = float(np.trapezoid(predicted_norm[mask], x_values[mask]))
        reference_area = float(np.trapezoid(reference_norm[mask], x_values[mask]))
        regions.append({
            "range_cm-1": [lower, upper],
            "predicted_normalized_area": predicted_area,
            "reference_normalized_area": reference_area,
            "area_error": predicted_area - reference_area,
        })
    return {
        "whole_spectrum_cosine_overlap": cosine,
        "equivalent_envelope_fwhm_cm-1": {
            "predicted": predicted_width, "reference": reference_width,
            "error": None if None in {predicted_width, reference_width} else predicted_width - reference_width,
        },
        "envelope_asymmetry": {
            "predicted": predicted_asymmetry, "reference": reference_asymmetry,
            "error": None if None in {predicted_asymmetry, reference_asymmetry}
            else predicted_asymmetry - reference_asymmetry,
        },
        "integrated_regions": regions,
    }


def _metrics(predicted_values: np.ndarray, reference_values: np.ndarray) -> dict:
    """Return comparable shape and absolute-error metrics for two equal arrays."""
    if len(predicted_values) < 3:
        raise ValueError("At least three overlapping spectrum points are required")
    rmse = float(np.sqrt(np.mean((predicted_values - reference_values) ** 2)))
    # A flat trace has no meaningful Pearson correlation.
    if np.std(predicted_values) == 0 or np.std(reference_values) == 0:
        correlation = None
    else:
        correlation = float(np.corrcoef(predicted_values, reference_values)[0, 1])
    return {"rmse": rmse, "pearson_correlation": correlation}


def _uniform_axis(px: np.ndarray, rx: np.ndarray, *, lower: float, upper: float) -> np.ndarray:
    if not lower < upper:
        return np.asarray([], dtype=float)
    step = min(float(np.median(np.diff(px))), float(np.median(np.diff(rx))))
    if not np.isfinite(step) or step <= 0:
        raise ValueError("Spectrum coordinate spacing must be finite and positive")
    intervals = max(1, int(np.floor((upper - lower) / step + 1e-12)))
    return np.linspace(lower, upper, intervals + 1)


def _shift_grid(window: float, step: float) -> list[float]:
    if window == 0:
        return [0.0]
    values = {-float(window), 0.0, float(window)}
    value = -window
    while value <= window + step * 1e-12:
        values.add(max(-window, min(window, float(value))))
        value += step
    return sorted(values)


def compare_spectra(
    predicted: Path,
    reference: Path,
    *,
    shift_window: float = 0.0,
    shift_step: float = 1.0,
) -> dict:
    """Compare spectra and optionally find the best global band-position shift.

    The shift search is a diagnostic, not a correction to the calculated spectrum:
    it shows how much mismatch is explained by a uniform frequency offset before
    attributing the remaining difference to line shape, intensity, or environment.
    """
    if not np.isfinite(shift_window) or shift_window < 0:
        raise ValueError("shift_window must be zero or positive")
    if not np.isfinite(shift_step) or shift_step <= 0:
        raise ValueError("shift_step must be positive")
    predicted_convention = _spectrum_convention(predicted)
    reference_convention = _spectrum_convention(reference)
    if "wavelength" in {predicted_convention["axis"], reference_convention["axis"]}:
        raise ValueError("IR benchmark inputs must use a wavenumber axis, not wavelength")
    for key in ("axis", "intensity"):
        left, right = predicted_convention[key], reference_convention[key]
        if left != "unspecified" and right != "unspecified" and left != right:
            raise ValueError(f"Predicted and reference spectrum {key} conventions differ ({left} versus {right})")
    px, py = read_spectrum_csv(predicted)
    rx, ry = read_spectrum_csv(reference)
    raw_lower, raw_upper = max(px.min(), rx.min()), min(px.max(), rx.max())
    x = _uniform_axis(px, rx, lower=float(raw_lower), upper=float(raw_upper))
    if len(x) < 3:
        raise ValueError("Predicted and reference spectra have insufficient overlapping range")
    interpolated = np.interp(x, px, py)
    reference_values = np.interp(x, rx, ry)
    raw = _metrics(interpolated, reference_values)
    shape = _shape_metrics(
        x, interpolated, reference_values,
        predicted_convention["intensity"], reference_convention["intensity"],
    )
    result = {
        "overlap_cm-1": [float(x.min()), float(x.max())],
        "points": int(len(x)),
        **raw,
        **shape,
        "raw": raw,
        "shift_search_cm-1": {"window": float(shift_window), "step": float(shift_step)},
        "sampling": "uniform_common_wavenumber_grid",
        "conventions": {"predicted": predicted_convention, "reference": reference_convention},
    }
    shifts = _shift_grid(float(shift_window), float(shift_step))
    search_lower = max(float(rx.min()), float(px.min() + shift_window))
    search_upper = min(float(rx.max()), float(px.max() - shift_window))
    search_x = _uniform_axis(px, rx, lower=search_lower, upper=search_upper)
    if len(search_x) < 3:
        raise ValueError("Requested shift window leaves fewer than three common comparison points")
    search_reference = np.interp(search_x, rx, ry)
    candidates = []
    for shift in shifts:
        shifted_predicted = np.interp(search_x, px + float(shift), py)
        candidate = _metrics(shifted_predicted, search_reference)
        candidates.append((float(shift), candidate))
    if not candidates:
        raise ValueError("No valid points available for the requested shift search")
    best_shift, aligned = max(
        candidates,
        key=lambda item: (
            item[1]["pearson_correlation"] is not None,
            item[1]["pearson_correlation"] if item[1]["pearson_correlation"] is not None else -float("inf"),
            -item[1]["rmse"],
            -abs(item[0]),
        ),
    )
    result["best_global_shift_cm-1"] = best_shift
    result["shift_aligned"] = {"overlap_cm-1": [float(search_x.min()), float(search_x.max())],
                               "points": int(len(search_x)), **aligned}
    return result


def write_comparison_plot(
    path: Path, *, predicted: Path, reference: Path, best_shift_cm_1: float = 0.0
) -> Path:
    """Write a QA overlay; the calculated trace is shifted only for diagnosis."""
    import matplotlib.pyplot as plt

    px, py = read_spectrum_csv(predicted)
    rx, ry = read_spectrum_csv(reference)
    figure, axis = plt.subplots(figsize=(10, 5))
    axis.plot(rx, ry, color="#202020", linewidth=1.2, label="reference")
    label = "calculated"
    if best_shift_cm_1:
        label += f" (diagnostic shift {best_shift_cm_1:+.0f} cm$^{{-1}}$)"
    axis.plot(px + best_shift_cm_1, py, color="#be123c", linewidth=1.1, label=label)
    axis.set_xlim(max(px.max() + best_shift_cm_1, rx.max()), min(px.min() + best_shift_cm_1, rx.min()))
    axis.set_xlabel("Wavenumber (cm$^{-1}$)")
    axis.set_ylabel("Spectrum value")
    axis.set_title("Calculated vs. reference IR spectrum")
    axis.legend(frameon=False)
    axis.grid(alpha=0.18)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return Path(path)


def write_benchmark_result(path: Path, *, predicted: Path, reference: Path, conditions: dict, metrics: dict) -> Path:
    result = {"predicted_spectrum": str(predicted), "reference_spectrum": str(reference), "conditions": conditions, "metrics": metrics}
    Path(path).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return Path(path)
