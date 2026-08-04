# src/utils/orca_xyz_clean.py
from pathlib import Path

def clean_xyz_for_orca(xyz_file: Path) -> list[str]:
    """
    Reads an XYZ file and returns only the atomic coordinates,
    skipping the first two lines (atom count + comment).
    """
    with open(xyz_file, "r") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise RuntimeError(f"XYZ file {xyz_file} is too short.")

    # Skip first two lines, return only coordinates
    coords = [line.strip() for line in lines[2:] if line.strip()]
    return coords
