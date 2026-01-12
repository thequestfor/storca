# src/inputgen.py
from pathlib import Path

def create_orca_input(
    xyz_file: Path,
    charge: int,
    multiplicity: int,
    opt: bool = False,
    freq: bool = False,
    label: str = "job",
    ncores: int = 1  # NEW: number of CPU cores
) -> Path:
    inp_file = Path(f"{label}.inp")
    lines = []

    # Method and basis
    keywords = ["B3LYP", "def2-SVP"]
    if opt:
        keywords.append("Opt")
    if freq:
        keywords.append("Freq")
    lines.append("! " + " ".join(keywords))

    # PAL block for parallel execution
    lines.append(f"%pal\n  nprocs {ncores}\nend")

    # Coordinates
    lines.append(f"* xyz {charge} {multiplicity}")

    with open(xyz_file) as f:
        xyz_lines = f.readlines()
        # Skip first two lines ONLY if they are atom count + comment
        try:
            int(xyz_lines[0].strip())
            start_idx = 2
        except (ValueError, IndexError):
            start_idx = 0
        lines.extend([line.strip() for line in xyz_lines[start_idx:] if line.strip()])

    lines.append("*")

    with open(inp_file, "w") as f:
        f.write("\n".join(lines))

    return inp_file

