import subprocess
from pathlib import Path

def run_orca(inp_file: Path, capture_out: bool = True) -> dict:
    """
    Run ORCA on the given input file.
    Returns a dictionary with:
        - gbw: Path to .gbw file
        - xyz: Path to optimized .xyz
        - out: Path to a captured ORCA output file (optional)
    """
    print(f"Running ORCA on {inp_file}...")

    out_file = inp_file.with_suffix(".out") if capture_out else None

    if capture_out:
        # Open a file to capture ORCA terminal output
        with open(out_file, "w") as f:
            result = subprocess.run(["orca", str(inp_file)], stdout=f, stderr=subprocess.STDOUT)
    else:
        result = subprocess.run(["orca", str(inp_file)])

    if result.returncode != 0:
        raise RuntimeError(f"ORCA failed: check {out_file if capture_out else 'terminal output'}")

    gbw_file = inp_file.with_suffix(".gbw")
    xyz_file = inp_file.with_suffix(".xyz")  # optimized geometry

    return {"gbw": gbw_file, "xyz": xyz_file, "out": out_file}

