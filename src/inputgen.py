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
    b97: bool = False, # if true, use b973c
    method_keywords: list[str] | None = None,
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

    if ncores < 1:
        raise ValueError("ncores must be at least 1")

    # Keep all ORCA input and output files with the geometry that spawned them.
    inp_file = xyz_file.parent / f"{label}.inp"
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
    elif method_keywords:
        keywords = list(method_keywords)
        if opt:
            keywords.append("Opt")
        if freq:
            keywords.append("Freq")
        if keepdens:
            keywords.append("KeepDens")
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


def _orca_method_line(method_keywords: list[str] | None, *job_keywords: str) -> str:
    """Return one deterministic ORCA simple-keyword line."""
    method = list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"])
    return "! " + " ".join([*method, *job_keywords])


def _write_orca_xyzfile_job(
    xyz_file: Path,
    *,
    label: str,
    charge: int,
    multiplicity: int,
    ncores: int,
    keyword_line: str,
    blocks: list[str] | None = None,
) -> Path:
    xyz_file = Path(xyz_file)
    if not xyz_file.is_file():
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")
    if ncores < 1:
        raise ValueError("ncores must be at least 1")
    path = xyz_file.parent / f"{label}.inp"
    lines = [keyword_line, f"%pal\n  nprocs {ncores}\nend", *(blocks or []),
             f"* xyzfile {charge} {multiplicity} {xyz_file.name}"]
    path.write_text("\n".join(lines) + "\n")
    return path


def create_orca_neb_ts_input(
    reactant_xyz: Path,
    product_xyz: Path,
    *,
    charge: int,
    multiplicity: int,
    label: str = "neb-ts",
    ncores: int = 1,
    nimages: int = 8,
    method_keywords: list[str] | None = None,
    preopt_ends: bool = False,
) -> Path:
    """Create an ORCA 6 double-ended NEB-TS input.

    ``preopt_ends`` is appropriate only when both supplied endpoints are
    physical bound minima.  It must remain false for assembled encounter
    geometries of separated collision partners, because relaxing those ends
    would silently change the declared channel.
    """
    product_xyz = Path(product_xyz)
    if not product_xyz.is_file():
        raise FileNotFoundError(f"Product XYZ file not found: {product_xyz}")
    if Path(reactant_xyz).parent.resolve() != product_xyz.parent.resolve():
        raise ValueError("NEB reactant and product XYZ files must share one working directory")
    if nimages < 3:
        raise ValueError("NEB-TS needs at least three movable images")
    return _write_orca_xyzfile_job(
        reactant_xyz,
        label=label,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        keyword_line=_orca_method_line(method_keywords, "TightSCF", "NEB-TS"),
        blocks=[
            f'%neb\n  Product "{product_xyz.name}"\n  NImages {nimages}'
            + ("\n  PreOptEnds true" if preopt_ends else "")
            + "\nend"
        ],
    )


def create_orca_irc_input(
    ts_xyz: Path,
    *,
    charge: int,
    multiplicity: int,
    hessian_file: Path | None = None,
    label: str = "irc",
    ncores: int = 1,
    max_iterations: int = 40,
    method_keywords: list[str] | None = None,
) -> Path:
    """Create a bidirectional ORCA IRC input from a validated TS geometry."""
    blocks = [f"%irc\n  MaxIter {max_iterations}\n  Direction both"]
    if hessian_file is not None:
        hessian_file = Path(hessian_file)
        if not hessian_file.is_file():
            raise FileNotFoundError(f"TS Hessian file not found: {hessian_file}")
        if hessian_file.parent.resolve() != Path(ts_xyz).parent.resolve():
            raise ValueError("IRC Hessian and TS XYZ must share one working directory")
        blocks[0] += f'\n  InitHess read\n  Hess_Filename "{hessian_file.name}"'
    blocks[0] += "\nend"
    return _write_orca_xyzfile_job(
        ts_xyz,
        label=label,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        keyword_line=_orca_method_line(method_keywords, "TightSCF", "IRC"),
        blocks=blocks,
    )


def create_orca_relaxed_scan_input(
    reactant_xyz: Path,
    *,
    bond_atom_indices: tuple[int, int],
    start_distance_angstrom: float,
    end_distance_angstrom: float,
    charge: int,
    multiplicity: int,
    label: str = "dissociation-scan",
    ncores: int = 1,
    steps: int = 20,
    method_keywords: list[str] | None = None,
) -> Path:
    """Create a full relaxed bond-distance scan for a dissociation route."""
    if steps < 3:
        raise ValueError("A relaxed path scan needs at least three points")
    left, right = bond_atom_indices
    scan = (
        "%geom\n"
        "  FullScan true\n"
        "  Scan\n"
        f"    B {left} {right} = {start_distance_angstrom:.8f}, {end_distance_angstrom:.8f}, {steps}\n"
        "  end\n"
        "end"
    )
    return _write_orca_xyzfile_job(
        reactant_xyz,
        label=label,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        keyword_line=_orca_method_line(method_keywords, "TightSCF", "Opt"),
        blocks=[scan],
    )


def create_orca_constrained_opt_input(
    xyz_file: Path,
    *,
    bond_atom_indices: tuple[int, int],
    distance_angstrom: float,
    charge: int,
    multiplicity: int,
    label: str = "constrained-opt",
    ncores: int = 1,
    max_iterations: int = 300,
    fresh_hessian: bool = False,
    method_keywords: list[str] | None = None,
) -> Path:
    """Create one independently resumable constrained ORCA optimization."""
    if max_iterations < 1:
        raise ValueError("max_iterations must be positive")
    left, right = bond_atom_indices
    hessian = "\n  Calc_Hess true" if fresh_hessian else ""
    geom = (
        "%geom\n"
        f"  MaxIter {max_iterations}{hessian}\n"
        "  Constraints\n"
        f"    {{ B {left} {right} {distance_angstrom:.8f} C }}\n"
        "  end\n"
        "end"
    )
    return _write_orca_xyzfile_job(
        xyz_file, label=label, charge=charge, multiplicity=multiplicity, ncores=ncores,
        keyword_line=_orca_method_line(method_keywords, "TightSCF", "Opt"), blocks=[geom],
    )



import textwrap
from pathlib import Path

def create_rmg_input(
    label: str,
    smiles: str,
    workdir: Path,
    temperature: float = 298.0,
    pressure: float = 1.0,
    termination_time: float = 1e7,
    max_edge_species: int = 250,
    additional_species: list[dict] | None = None,
    initial_mole_fractions: dict[str, float] | None = None,
    reaction_libraries: list[Path] | None = None,
    database_reaction_libraries: list[str] | None = None,
    cap_generated_carbon_at_target: bool = True,
    maximum_heavy_atoms: int | None = None,
    maximum_radical_electrons: int | None = None,
    filter_reactions: bool = False,
    max_objects_per_iteration: int = 1,
    restart_from_seed: Path | None = None,
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
        additional_species: Optional RMG species dictionaries with ``label``,
            ``smiles``, and optional ``reactive`` keys.
        initial_mole_fractions: Explicit reactor composition. Defaults to the
            target species at mole fraction one.
        reaction_libraries: Validated local STORCA-generated library paths.
        database_reaction_libraries: Named libraries shipped with the selected
            RMG database. Local calculated libraries are listed first so an
            ORCA/Arkane repair can replace a reference-library rate.

    Returns:
        Path to input.py
    """
    workdir.mkdir(parents=True, exist_ok=True)
    input_file = workdir / "input.py"
    additional_species = [dict(species) for species in (additional_species or [])]
    reaction_libraries = reaction_libraries or []
    database_reaction_libraries = list(database_reaction_libraries or [])
    if any(not str(name).strip() for name in database_reaction_libraries):
        raise ValueError("RMG database reaction-library names must be nonempty")
    database_reaction_libraries = list(dict.fromkeys(str(name).strip() for name in database_reaction_libraries))
    for library in reaction_libraries:
        library = Path(library)
        if not (library / "reactions.py").is_file() or not (library / "dictionary.txt").is_file():
            raise ValueError(f"RMG reaction library must contain reactions.py and dictionary.txt: {library}")
    mole_fractions = initial_mole_fractions or {label: 1.0}
    # RMG requires every declared species graph to be unique.  A stability
    # target can itself be water/oxygen/a bath-gas component, in which case
    # the environmental fraction belongs to the same molecular pool.
    try:
        from rdkit import Chem

        def canonical(value: str) -> str:
            molecule = Chem.MolFromSmiles(value)
            if molecule is None:
                raise ValueError(f"Invalid RMG species SMILES: {value}")
            return Chem.MolToSmiles(molecule, canonical=True)

        target_canonical = canonical(smiles)
        retained_species = []
        seen = {target_canonical: label}
        for species in additional_species:
            if not {"label", "smiles"} <= species.keys():
                raise ValueError("Each additional RMG species needs label and smiles keys")
            species_label = str(species["label"])
            species_canonical = canonical(str(species["smiles"]))
            if species_canonical == target_canonical:
                amount = float(mole_fractions.get(species_label, 0.0))
                mole_fractions = dict(mole_fractions)
                mole_fractions.pop(species_label, None)
                mole_fractions[label] = float(mole_fractions.get(label, 0.0)) + amount
                continue
            if species_canonical in seen:
                raise ValueError(
                    f"RMG additional species {species_label} duplicates {seen[species_canonical]} structurally"
                )
            seen[species_canonical] = species_label
            retained_species.append(species)
        additional_species = retained_species
    except ImportError:
        # The normal STORCA environment includes RDKit. Keep the prior input
        # validation behavior if a minimal environment does not.
        pass
    declared_labels = {label, *(species["label"] for species in additional_species)}
    if set(mole_fractions) - declared_labels:
        raise ValueError("RMG initial mole fractions reference an undeclared species")
    if any(value < 0 for value in mole_fractions.values()) or abs(sum(mole_fractions.values()) - 1.0) > 1e-9:
        raise ValueError("RMG initial mole fractions must be non-negative and sum to one")
    if max_objects_per_iteration < 1:
        raise ValueError("max_objects_per_iteration must be at least one")
    if maximum_heavy_atoms is not None and maximum_heavy_atoms < 1:
        raise ValueError("maximum_heavy_atoms must be positive when supplied")
    if maximum_radical_electrons is not None and maximum_radical_electrons < 0:
        raise ValueError("maximum_radical_electrons must be non-negative when supplied")
    if restart_from_seed is not None and not Path(restart_from_seed).is_dir():
        raise ValueError(f"RMG restart seed directory not found: {restart_from_seed}")
    species_blocks = [
        f"species(\n    label={label!r},\n    reactive=True,\n    structure=SMILES({smiles!r}),\n)"
    ]
    for species in additional_species:
        species_blocks.append(
            f"species(\n    label={species['label']!r},\n    reactive={bool(species.get('reactive', True))!r},\n"
            f"    structure=SMILES({species['smiles']!r}),\n)"
        )
    species_text = "\n\n".join(species_blocks)
    fractions_text = textwrap.indent(
        "\n".join(f"{name!r}: {fraction}," for name, fraction in mole_fractions.items()),
        "            ",
    )
    # A stability screen follows loss of one declared target molecule.  For an
    # organic target, species containing more carbon than the target must have
    # arisen through target-target association or an undeclared organic
    # co-reactant.  Excluding those generated species prevents dimer/polymer
    # explosions from consuming the bounded screen, without blocking oxidation,
    # hydration, fragmentation, or isomerization (which do not add carbon).
    maximum_carbon = None
    if cap_generated_carbon_at_target:
        try:
            from rdkit import Chem
            target = Chem.MolFromSmiles(smiles)
            if target is not None:
                maximum_carbon = sum(atom.GetAtomicNum() == 6 for atom in target.GetAtoms()) or None
        except ImportError:
            # The input remains usable without RDKit; it simply cannot apply
            # the organic scaffold cap in that environment.
            pass
    carbon_constraint = f"\n    maximumCarbonAtoms={maximum_carbon}," if maximum_carbon is not None else ""
    heavy_atom_constraint = (
        f"\n    maximumHeavyAtoms={int(maximum_heavy_atoms)}," if maximum_heavy_atoms is not None else ""
    )
    radical_constraint = (
        f"\n    maximumRadicalElectrons={int(maximum_radical_electrons)},"
        if maximum_radical_electrons is not None else ""
    )

    library_entries = [
        *(str(Path(library).resolve()) for library in reaction_libraries),
        *database_reaction_libraries,
    ]

    content = f"""database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries={library_entries!r},
    seedMechanisms=[],
    kineticsDepositories=['training'],
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

{species_text}

generatedSpeciesConstraints(
    allowed=('input species',),{carbon_constraint}{heavy_atom_constraint}{radical_constraint}
)

simpleReactor(
    temperature=({temperature}, 'K'),
    pressure=({pressure}, 'bar'),
    initialMoleFractions={{
{fractions_text}
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
    maximumEdgeSpecies={max_edge_species},
    maxNumObjsPerIter={max_objects_per_iteration},
    filterReactions={bool(filter_reactions)!r},
)

options(
    units='si',
    saveSimulationProfiles=True,
    generateOutputHTML=True,
    generatePlots=False
)
"""
    if restart_from_seed is not None:
        content += f"\nrestartFromSeed(path={str(Path(restart_from_seed).resolve())!r})\n"

    input_file.write_text(content.strip() + "\n")
    return input_file
