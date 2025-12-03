import subprocess
from pathlib import Path

def gbw_to_wfn(gbw_file: Path) -> Path:
    """
    Convert an ORCA .gbw file to a .wfn file using orca_2aim.
    Returns the path to the generated .wfn file.
    """
    if gbw_file.suffix != ".gbw":
        raise ValueError("Input file must have a .gbw extension")

    base_name = gbw_file.with_suffix("")  # remove .gbw
    wfn_file = base_name.with_suffix(".wfn")

    print(f"Converting {gbw_file} → {wfn_file} using orca_2aim...")

    # Pass the base name (without .gbw) to orca_2aim
    result = subprocess.run(["orca_2aim", str(base_name)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"orca_2aim failed:\n{result.stderr}")

    if not wfn_file.exists():
        raise FileNotFoundError(f"Expected WFN file not found: {wfn_file}")

    return wfn_file

