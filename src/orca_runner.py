from pathlib import Path
import shutil
import subprocess

from pathlib import Path
import subprocess

def find_orca():
    """Return absolute path to ORCA executable or exit with error."""
    result = subprocess.run(
        ["which", "orca"],
        capture_output=True,
        text=True
    )
    if result.returncode != 0 or not result.stdout.strip():
        print("❌ ORCA executable not found in PATH")
        sys.exit(1)
    return result.stdout.strip()

def run_orca(inp_file: Path, capture_out: bool = True) -> dict:
    """
    Run ORCA on the given input file in serial mode.
    Produces .gbw, .xyz, and optional .out file.
    """
    parresult = subprocess.run(
        ["which", "orca"],
        capture_output=True,
        text=True
    )
    if parresult.returncode != 0 or not parresult.stdout.strip():
        print("❌ ORCA executable not found in PATH")
        sys.exit(1)

    ORCA_CMD = parresult.stdout.strip()  # this is the path to ORCA
    out_file = inp_file.with_suffix(".out") if capture_out else None

    if capture_out:
        with open(out_file, "w") as f:
            result = subprocess.run([ORCA_CMD, str(inp_file)], stdout=f, stderr=subprocess.STDOUT)
    else:
        result = subprocess.run([ORCA_CMD, str(inp_file)])

    if result.returncode != 0:
        raise RuntimeError(f"ORCA failed: check {out_file if capture_out else 'terminal output'}")

    gbw_file = inp_file.with_suffix(".gbw")
    xyz_file = inp_file.with_suffix(".xyz")

    return {"gbw": gbw_file, "xyz": xyz_file, "out": out_file}
