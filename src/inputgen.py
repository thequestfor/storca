# src/inputgen.py
from pathlib import Path
import subprocess
import textwrap
import re
import shutil


class XTBNotFoundError(RuntimeError):
    """Raised when xTB executable is required but not found."""
    pass


def create_orca_input(
    xyz_file: Path,
    charge: int,
    multiplicity: int,
    opt: bool = False,
    freq: bool = False,
    keepdens: bool = False,
    label: str = "job",
    ncores: int = 1,      # number of CPU cores
    use_goat: bool = False,  # if True, generate GOAT/xTB input
    b97: bool = False # if true, use b973c
) -> Path:
    """
    Generate an ORCA input file.
    If use_goat=True, generates a GOAT input using xTB as the method.

    Args:
        xyz_file: Path to initial geometry file.
        charge: Molecular charge.
        multiplicity: Spin multiplicity.
        opt: Add optimization keyword.
        freq: Add frequency calculation keyword.
        keepdens: Add KeepDens keyword.
        label: Output file base name.
        ncores: Number of CPU cores.
        use_goat: If True, generate GOAT/xTB input.
        b97: If true, use b97-3c

    Returns:
        Path to the generated ORCA input file.

    Raises:
        XTBNotFoundError: if use_goat=True but xTB executable is not found.
        FileNotFoundError: if xyz_file does not exist.
    """

    if not xyz_file.exists():
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")

    inp_file = Path(f"{label}.inp")
    lines = []

    # Method and basis / keywords
    if use_goat:
        xtb_exe = shutil.which("xtb")
        if xtb_exe is None:
            raise XTBNotFoundError(
                "xTB executable not found in PATH. GOAT input requires xTB. "
                "Please install xTB and ensure it is in your PATH."
            )
        keywords = ["GOAT", "xTB"]
    elif b97:
        keywords = ["B97-3c"]
        if opt:
            keywords.append("Opt")
        if freq:
            keywords.append("Freq")
        if keepdens:
            keywords.append("KeepDens")
    else:
        keywords = ["B3LYP", "def2-SVP", "NoAutoStart"]
        if opt:
            keywords.append("Opt")
        if freq:
            keywords.append("Freq")
        if keepdens:
            keywords.append("KeepDens")

    lines.append("! " + " ".join(keywords))

    # PAL block for parallel execution
    lines.append(f"%pal\n  nprocs {ncores}\nend")

    # Coordinates
    lines.append(f"* xyz {charge} {multiplicity}")

    with xyz_file.open() as f:
        xyz_lines = f.readlines()
        # Skip first two lines ONLY if they are atom count + comment
        try:
            int(xyz_lines[0].strip())
            start_idx = 2
        except (ValueError, IndexError):
            start_idx = 0
        lines.extend([line.strip() for line in xyz_lines[start_idx:] if line.strip()])

    lines.append("*")

    inp_file.write_text("\n".join(lines))

    return inp_file



import textwrap
from pathlib import Path

def create_rmg_input(
    label: str,
    smiles: str,
    workdir: Path,
    temperature: float = 298.0,
    pressure: float = 1.0,
    termination_time: float = 1e7,
) -> Path:
    """
    Generate an RMG input.py file for stability testing.
    This version does NOT set an output directory, so logs go directly
    into the current working directory (workdir).

    Args:
        label: Species label
        smiles: SMILES string
        workdir: Directory where input/output live
        temperature: K
        pressure: bar
        termination_time: seconds

    Returns:
        Path to input.py
    """
    workdir.mkdir(parents=True, exist_ok=True)
    input_file = workdir / "input.py"

    content = textwrap.dedent(f"""
    database(
        thermoLibraries=['primaryThermoLibrary'],
        reactionLibraries=[],
        seedMechanisms=[],
        kineticsDepositories=['training'],
        kineticsFamilies='default',
        kineticsEstimator='rate rules',
    )

    species(
        label='{label}',
        reactive=True,
        structure=SMILES('{smiles}'),
    )

    generatedSpeciesConstraints(
        allowed=('input species',)
    )

    simpleReactor(
        temperature=({temperature}, 'K'),
        pressure=({pressure}, 'bar'),
        initialMoleFractions={{
            '{label}': 1.0,
        }},
        terminationTime=({termination_time}, 's'),
    )

    simulator(
        atol=1e-16,
        rtol=1e-8,
    )

    model(
        toleranceMoveToCore=1e-1,
        toleranceInterruptSimulation=1e8,
        maximumEdgeSpecies=1000,
    )

    options(
        units='si',
        saveSimulationProfiles=True,
        generateOutputHTML=True,
        generatePlots=False
    )
    """)

    input_file.write_text(content.strip() + "\n")
    return input_file
