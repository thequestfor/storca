# src/stability/freq_check.py
from pathlib import Path
import re

FREQ_HEADER = re.compile(r"VIBRATIONAL FREQUENCIES", re.IGNORECASE)
FREQ_VALUE = re.compile(r"\b\d+[:\s]+(-?\d+\.\d+)\b")


def parse_orca_frequencies(out_file: Path) -> list[float]:
    """
    Parse vibrational frequencies from an ORCA output file.
    Returns a list of frequencies in cm^-1.
    """
    if not out_file.exists():
        raise FileNotFoundError(f"ORCA output file not found: {out_file}")

    freqs: list[float] = []
    in_freq_block = False

    with out_file.open("r") as f:
        for line in f:
            if FREQ_HEADER.search(line):
                in_freq_block = True
                continue

            if in_freq_block:
                match = FREQ_VALUE.search(line)
                if match:
                    freqs.append(float(match.group(1)))

                if "NORMAL MODES" in line or "THERMOCHEMISTRY" in line:
                    break

    if not freqs:
        raise RuntimeError(
            f"No vibrational frequencies found in ORCA output: {out_file}"
        )

    return freqs


def frequency_stability_check(out_file: Path, imag_tol: float = 20.0) -> dict:
    """
    Determines whether a structure is a true minimum.
    imag_tol: frequencies more negative than -imag_tol (cm^-1)
              are treated as imaginary
    """
    freqs = parse_orca_frequencies(out_file)

    imaginary = [f for f in freqs if f < -imag_tol]
    near_zero = [f for f in freqs if -imag_tol <= f <= imag_tol]

    return {
        "IsMinimum": len(imaginary) == 0,
        "NumImaginary": len(imaginary),
        "LowestFrequency": min(freqs),
        "ImaginaryFrequencies": imaginary,
        "NearZeroModes": near_zero,
        "AllFrequencies": freqs,
    }
