"""Calibrated extraction of line spectra from simple transmittance plot images."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


@dataclass(frozen=True)
class PlotCalibration:
    left: int
    right: int
    top: int
    bottom: int
    x_left_value: float = 4000.0
    x_right_value: float = 400.0
    y_top_value: float = 100.0
    y_bottom_value: float = 0.0
    # Optional (pixel_x, wavenumber) pairs for a non-linear printed axis.
    x_anchors: tuple[tuple[float, float], ...] = ()


# Calibration for the supplied 715×553 SDBS-style images. SDBS deliberately
# expands the fingerprint region: its printed 4000, 3000, 2000, 1500, 1000,
# 500 cm^-1 ticks are not uniformly spaced by wavenumber. These anchors were
# read from the actual tick marks, rather than inferred from the plot edges.
SDBS_IR_CALIBRATION = PlotCalibration(
    left=29,
    right=713,
    top=96,
    bottom=417,
    x_anchors=((29, 4000), (161, 3000), (293, 2000), (425, 1500), (556, 1000), (688, 500), (713, 400)),
)


def digitize_transmittance(image_path: Path, calibration: PlotCalibration = SDBS_IR_CALIBRATION) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Trace a single dark transmittance curve inside a calibrated plot area.

    Returns wavenumber, transmittance, and pixel y-coordinate arrays. This is
    appropriate for clean monochrome line plots, not photographs of plots.
    """
    rgb = np.asarray(Image.open(image_path).convert("RGB"))
    grayscale = np.asarray(Image.open(image_path).convert("L"))
    left, right, top, bottom = calibration.left, calibration.right, calibration.top, calibration.bottom
    if not (0 <= left < right < grayscale.shape[1] and 0 <= top < bottom < grayscale.shape[0]):
        raise ValueError("Plot calibration is outside image bounds")
    rgb_values = rgb.astype(np.int16)
    red_dominant = (rgb_values[:, :, 0] > rgb_values[:, :, 1] + 25) & (rgb_values[:, :, 0] > rgb_values[:, :, 2] + 25)
    crop_area = (right - left + 1) * (bottom - top + 1)
    # NIST-style exports commonly use a red trace. Prefer that colour mask
    # when it is substantial, avoiding black plot borders and tick marks.
    trace_mask = red_dominant if red_dominant[top:bottom + 1, left:right + 1].sum() > crop_area * 0.001 else grayscale < 90
    trace_y: list[float] = []
    observed: list[bool] = []
    previous = top + 18.0
    for x in range(left + 2, right - 1):
        candidates = np.where(trace_mask[top + 2:bottom - 1, x])[0] + top + 2
        if not len(candidates):
            trace_y.append(float("nan"))
            observed.append(False)
            continue
        # The trace is continuous; choosing the nearest dark ink suppresses
        # labels, borders, and sparse scan artefacts within the plot area.
        selected = float(candidates[np.argmin(np.abs(candidates - previous))])
        trace_y.append(selected)
        observed.append(True)
        previous = selected
    pixels = np.asarray(trace_y, dtype=float)
    observed_mask = np.asarray(observed, dtype=bool)
    if observed_mask.mean() < 0.10:
        raise ValueError("Too little trace ink was detected to digitize this spectrum reliably")
    sample_positions = np.arange(len(pixels), dtype=float)
    pixels = np.interp(sample_positions, sample_positions[observed_mask], pixels[observed_mask])
    xs = np.arange(left + 2, right - 1)
    if calibration.x_anchors:
        anchor_pixels, anchor_wavenumbers = zip(*calibration.x_anchors)
        if (not all(np.isfinite(value) for value in (*anchor_pixels, *anchor_wavenumbers))
                or any(right <= left for left, right in zip(anchor_pixels, anchor_pixels[1:]))):
            raise ValueError("x_anchors must be finite and strictly ordered from left to right")
        wavenumbers = np.interp(xs, anchor_pixels, anchor_wavenumbers)
    else:
        wavenumbers = calibration.x_left_value + (xs - left) * (calibration.x_right_value - calibration.x_left_value) / (right - left)
    transmittance = calibration.y_top_value + (pixels - top) * (calibration.y_bottom_value - calibration.y_top_value) / (bottom - top)
    return wavenumbers, np.clip(transmittance, calibration.y_bottom_value, calibration.y_top_value), pixels


def write_digitized_csv(path: Path, wavenumbers: np.ndarray, transmittance: np.ndarray) -> Path:
    wavenumbers = np.asarray(wavenumbers, dtype=float)
    transmittance = np.asarray(transmittance, dtype=float)
    if (wavenumbers.ndim != 1 or transmittance.ndim != 1 or len(wavenumbers) == 0
            or len(wavenumbers) != len(transmittance)
            or not np.all(np.isfinite(wavenumbers)) or not np.all(np.isfinite(transmittance))):
        raise ValueError("Digitized coordinates must be equal-length finite one-dimensional arrays")
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber_cm-1", "transmittance_percent"])
        writer.writerows(zip(wavenumbers, transmittance))
    return Path(path)


def write_trace_overlay(image_path: Path, output_path: Path, calibration: PlotCalibration, trace_y: np.ndarray) -> Path:
    image = Image.open(image_path).convert("RGB")
    draw = ImageDraw.Draw(image)
    points = [(calibration.left + 2 + index, float(y)) for index, y in enumerate(trace_y)]
    draw.line(points, fill=(220, 20, 60), width=1)
    draw.rectangle((calibration.left, calibration.top, calibration.right, calibration.bottom), outline=(30, 90, 220), width=1)
    image.save(output_path)
    return Path(output_path)
