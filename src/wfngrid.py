import subprocess
from pathlib import Path

def run_multiwfn_fukui(neutral_wfn: Path, cation_wfn: Path, anion_wfn: Path):
    """
    Run Multiwfn to compute Fukui grid and export cube files.
    
    Args:
        neutral_wfn: Path to neutral .wfn file
        cation_wfn: Path to cation (+1) .wfn file
        anion_wfn: Path to anion (-1) .wfn file
    """
    neutral_wfn = Path(neutral_wfn).resolve()
    cation_wfn = Path(cation_wfn).resolve()
    anion_wfn = Path(anion_wfn).resolve()

    print("Running Fukui grid calculation in Multiwfn (this may take a while)...")

    # -------------------------------
    # Step 1: Compute Fukui grid
    # -------------------------------
    cmds_grid = f"22\n3\n{neutral_wfn}\n{cation_wfn}\n{anion_wfn}\n3\n"
    subprocess.run(
        f'echo -e "{cmds_grid}" | Multiwfn {neutral_wfn}',
        shell=True,
        check=True
    )

    print("Fukui grid calculation finished.")

    # -------------------------------
    # Step 2: Export grids as .cub
    # -------------------------------
    print("Exporting grid data to .cub files...")

    cmds_export = "22\n5\n6\n7\n8\n0\n"
    subprocess.run(
        f'echo -e "{cmds_export}" | Multiwfn {neutral_wfn}',
        shell=True,
        check=True
    )

    print("Export finished. Renaming .cub → .cube for Avogadro 2...")

    # -------------------------------
    # Step 3: Rename .cub → .cube
    # -------------------------------
    for cub_file in Path(".").glob("*.cub"):
        cube_file = cub_file.with_suffix(".cube")
        cub_file.rename(cube_file)
        print(f"Renamed: {cub_file} → {cube_file}")

    print("All Fukui grids are ready for visualization.")

