from pathlib import Path
from typing import List, Tuple
from collections import defaultdict
import os
import re

def parse_orca_orbitals(out_file: Path) -> dict:
    orbitals = []
    reading = False
    table_started = False

    with open(out_file) as f:
        for line in f:
            line_strip = line.strip()
            if "ORBITAL ENERGIES" in line_strip:
                reading = True
                continue
            if reading:
                if set(line_strip) == {"-"}:
                    continue  # skip separator
                if not table_started:
                    if "NO" in line_strip and "OCC" in line_strip:
                        table_started = True
                    continue
                # stop at empty line or notes
                if line_strip == "" or "Only the first" in line_strip:
                    break
                # parse orbital line
                parts = line_strip.split()
                if len(parts) >= 4:
                    try:
                        no = int(parts[0])
                        occ = float(parts[1])
                        eh = float(parts[2])
                        ev = float(parts[3])
                        orbitals.append({"no": no, "occ": occ, "eh": eh, "ev": ev})
                    except ValueError:
                        continue

    if not orbitals:
        return {"homo_number": None, "homo_energy": None,
                "lumo_number": None, "lumo_energy": None}

    homo_orb = max([o for o in orbitals if o["occ"] > 0], key=lambda x: x["no"])
    lumo_orb = min([o for o in orbitals if o["occ"] == 0], key=lambda x: x["no"])

    return {
        "homo_number": homo_orb["no"],
        "homo_energy": homo_orb["ev"],
        "lumo_number": lumo_orb["no"],
        "lumo_energy": lumo_orb["ev"]
    }

def parse_orca_energy(out_file: Path) -> float:
    """
    Extract the final SCF energy from an ORCA output file.
    Looks for the last 'Energy :' line and returns the value in Hartree.
    """
    energy = None
    with open(out_file) as f:
        for line in f:
            if "Energy" in line and "Eh" in line:
                try:
                    # Split on ':' and take the numeric part before 'Eh'
                    energy_str = line.split(':')[1].split()[0]
                    energy = float(energy_str)
                except Exception:
                    continue
    if energy is None:
        raise RuntimeError(f"Could not find SCF energy in {out_file}")
    return energy

# --------------------------
# GOAT OUTPUT PARSER
# --------------------------
def parse_goat_out(out_file):
    """
    Parse GOAT .out file for conformers and Boltzmann weighings.
    Returns:
        conformers: list of all conformers with their weights
        sig_conformers: list of conformers with weight > 0.05
    """
    conformers = []

    with open(out_file, "r") as f:
        lines = f.readlines()

    # Find the start of the Boltzmann chart
    chart_start = None
    for i, line in enumerate(lines):
        if re.match(r'\s*Conformer\s+Energy\s+Degen\.', line):
            chart_start = i + 2  # skip header/dashed line
            break

    if chart_start is None:
        raise ValueError("Boltzmann chart not found in .out file")

    # Parse the chart lines
    for line in lines[chart_start:]:
        if line.strip() == "":
            break
        parts = line.split()
        if len(parts) >= 5:
            idx = int(parts[0])
            energy = float(parts[1])
            degeneracy = int(parts[2])
            weight = float(parts[3]) / 100.0  # convert percent -> fraction
            cum_weight = float(parts[4]) / 100.0
            conformers.append({
                "idx": idx,
                "energy": energy,
                "degeneracy": degeneracy,
                "weight": weight,
                "cum_weight": cum_weight
            })

    # Significant conformers: weight > 0.05
    sig_conformers = [c for c in conformers if c["weight"] > 0.05]

    return conformers, sig_conformers


# --------------------------
# XYZ PARSER
# --------------------------
def parse_xyz_ensemble(xyz_file: Path):
    """
    Parse an XYZ file containing multiple frames.
    Returns a list of frames:
        [{'natoms': int, 'comment': str, 'xyz_lines': list[str]}]
    """
    frames = []
    with open(xyz_file, "r") as f:
        lines = [line.rstrip("\n") for line in f]

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            natoms = int(line)
        except ValueError:
            raise ValueError(f"Expected number of atoms at line {i+1}, got: {line}")

        if i + 1 + natoms >= len(lines):
            raise ValueError(f"Incomplete frame starting at line {i+1}")

        comment = lines[i + 1]
        xyz_lines = lines[i + 2 : i + 2 + natoms]
        frames.append({
            "natoms": natoms,
            "comment": comment,
            "xyz_lines": xyz_lines
        })
        i += 2 + natoms

    return frames



# --------------------------
# HELPER: COMBINE OUT & XYZ
# --------------------------
def attach_xyz_to_conformers(conformers, xyz_frames):
    """
    Attach XYZ frames to parsed conformers.
    Assumes frames are in the same order as conformers.
    """
    if len(conformers) != len(xyz_frames):
        raise ValueError(
            f"Mismatch: {len(conformers)} conformers vs {len(xyz_frames)} xyz frames"
        )

    for c, frame in zip(conformers, xyz_frames):
        c["xyz"] = frame

    return conformers



def merge_weighted_ir(
    conformers,
    freq_tol: float = 10.0
):
    """
    Merge Boltzmann-weighted IR modes into a single stick spectrum.

    Returns:
        List of (freq_cm1, intensity_km_mol)
    """
    all_modes = []

    for conf in conformers:
        if conf["freqs"] is None or conf["ir_intensities"] is None:
            continue
        w = conf["weight"]
        for f, I in zip(conf["freqs"], conf["ir_intensities"]):
            all_modes.append((f, I * w))

    # Sort by frequency
    all_modes.sort(key=lambda x: x[0])

    merged = []
    current = []

    for f, I in all_modes:
        if not current:
            current = [(f, I)]
            continue
        if abs(f - current[-1][0]) <= freq_tol:
            current.append((f, I))
        else:
            merged.append(current)
            current = [(f, I)]

    if current:
        merged.append(current)

    # Weighted average frequency, summed intensity
    final_spectrum = []
    for group in merged:
        total_I = sum(I for _, I in group)
        if total_I <= 0.0:
            continue
        avg_f = sum(f * I for f, I in group) / total_I
        final_spectrum.append((avg_f, total_I))

    return final_spectrum

from pathlib import Path

def parse_orca_ir(file_path: str):
    """
    Parse ORCA IR output file for frequencies, epsilon, and intensity.

    Returns:
        List of dicts:
            [{'mode': int, 'freq': float, 'eps': float, 'intensity': float}, ...]
    """
    vibrations = []

    with open(file_path) as f:
        lines = f.readlines()

    # Find the line that says "IR SPECTRUM"
    for i, line in enumerate(lines):
        if "IR SPECTRUM" in line:
            start_data = i + 6  # skip 6 lines to get to actual data
            break
    else:
        raise RuntimeError("IR SPECTRUM section not found in ORCA output")

    # Parse table
    for line in lines[start_data:]:
        line = line.strip()
        if not line or line.startswith("*"):
            # Stop at blank line or footer
            break

        # Split tokens
        tokens = line.split()
        if len(tokens) < 4:
            # Not enough columns; skip
            continue

        # Extract data
        mode = int(tokens[0].rstrip(":"))
        freq = float(tokens[1])
        eps = float(tokens[2])
        intensity = float(tokens[3])

        vibrations.append({
            "mode": mode,
            "freq": freq,
            "eps": eps,
            "intensity": intensity
        })

    return vibrations


