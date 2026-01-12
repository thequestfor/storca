# src/utils/orca_xyz_clean.py
from pathlib import Path

def clean_xyz_lines(xyz_file: Path) -> Path:
    """
    Ensure ORCA-friendly XYZ: remove first two lines if present.
    Returns path to cleaned temporary XYZ.
    """
    cleaned_xyz = xyz_file.with_name(f"{xyz_file.stem}_clean.xyz")
    with open(xyz_file) as f_in, open(cleaned_xyz, "w") as f_out:
        lines = f_in.readlines()
        # Detect and skip first 2 lines if they are atom count + comment
        start_idx = 0
        if len(lines) >= 2:
            try:
                int(lines[0].strip())  # atom count line
                start_idx = 2  # skip atom count + comment
            except ValueError:
                start_idx = 0
        f_out.writelines(lines[start_idx:])
    return cleaned_xyz
