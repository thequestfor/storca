# src/inputgen.py
from pathlib import Path
import math
import subprocess
import textwrap
import re
import shutil
import hashlib

from .orca_job_settings import NativeJobSettings, apply_native_job_settings


class XTBNotFoundError(RuntimeError):
    """Raised when xTB executable is required but not found."""
    pass


def normalize_thermochemistry_temperature(value: float) -> float:
    """One finite positive temperature; grids need a separate result contract."""
    if isinstance(value, bool):
        raise ValueError("Thermochemistry temperature must be one finite positive number")
    try:
        temperature = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("Thermochemistry temperature must be one finite positive number") from error
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("Thermochemistry temperature must be one finite positive number")
    return temperature


def orca_frequency_temperature_input_evidence(input_path: Path, temperature_K: float) -> dict:
    """Prove one explicit %freq Temp request without executing ORCA.

    This narrow producer emits a single thermochemistry temperature. Multiple
    blocks/temperatures and implicit defaults cannot supply explicit evidence.
    """
    temperature = normalize_thermochemistry_temperature(temperature_K)
    input_path = Path(input_path)
    if not input_path.is_file():
        raise ValueError("Requested-temperature thermochemistry input evidence is missing")
    data = input_path.read_bytes()
    text = data.decode("utf-8", errors="replace")
    keyword_lines = " ".join(re.findall(r"(?m)^\s*!([^\n]*)", text))
    if not re.search(r"\b(?:Freq|AnFreq|NumFreq)\b", keyword_lines, flags=re.IGNORECASE):
        raise ValueError("Thermochemistry temperature evidence requires a frequency input")
    blocks = re.findall(r"%freq\b(.*?)\bend\b", text, flags=re.IGNORECASE | re.DOTALL)
    if len(blocks) != 1:
        raise ValueError("Thermochemistry requires exactly one explicit %freq temperature block")
    declarations = re.findall(r"\b(?:Temp|T)\s+([^\n#]+)", blocks[0], flags=re.IGNORECASE)
    if len(declarations) != 1 or not re.fullmatch(r"\s*[-+0-9.Ee]+\s*", declarations[0]):
        raise ValueError("Thermochemistry requires one temperature, not a grid or conflicting overrides")
    retained = normalize_thermochemistry_temperature(declarations[0].strip())
    if not math.isclose(retained, temperature, rel_tol=1e-14, abs_tol=1e-10):
        raise ValueError("Retained frequency input temperature differs from the requested temperature")
    return {
        "status": "matched", "source_input": str(input_path),
        "source_input_sha256": hashlib.sha256(data).hexdigest(),
        "requested_temperature_K": temperature, "input_temperature_K": retained,
        "orca_control": "%freq Temp", "temperature_count": 1,
    }


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
    preserve_existing: bool = False,
    thermochemistry_temperature_K: float | None = None,
    native_settings: NativeJobSettings | None = None,
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
        thermochemistry_temperature_K: Explicit single thermochemistry temperature
            for Freq only. None preserves the legacy ORCA default input bytes.

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
    temperature = None
    if thermochemistry_temperature_K is not None:
        if not freq or use_goat:
            raise ValueError("Explicit thermochemistry temperature requires a frequency job, not GOAT")
        temperature = normalize_thermochemistry_temperature(thermochemistry_temperature_K)

    # Keep all ORCA input and output files with the geometry that spawned them.
    inp_file = xyz_file.parent / f"{label}.inp"
    lines = []

    # Method and basis / keywords
    if use_goat:
        try:
            from src.orca_runner import find_xtb
            find_xtb()
        except RuntimeError as error:
            raise XTBNotFoundError(
                "GOAT requires xTB, but no executable was found in the active Python environment or PATH. "
                "Install xTB or set STORCA_XTB_BIN."
            ) from error
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
    if temperature is not None:
        # ORCA 6.1 manual: Frequencies keyword list, %freq Temp.
        # One T avoids ambiguity in downstream scalar thermochemistry parsing.
        lines.append(f"%freq\n  Temp {temperature:.15g}\nend")

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

    input_text = apply_native_job_settings("\n".join(lines), native_settings)
    if preserve_existing and inp_file.exists():
        if inp_file.read_text() != input_text:
            raise RuntimeError(f"Retained ORCA input differs; refusing to overwrite it: {inp_file}")
    else:
        inp_file.write_text(input_text)

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
    native_settings: NativeJobSettings | None = None,
) -> Path:
    xyz_file = Path(xyz_file)
    if not xyz_file.is_file():
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")
    if ncores < 1:
        raise ValueError("ncores must be at least 1")
    path = xyz_file.parent / f"{label}.inp"
    lines = [keyword_line, f"%pal\n  nprocs {ncores}\nend", *(blocks or []),
             f"* xyzfile {charge} {multiplicity} {xyz_file.name}"]
    text = apply_native_job_settings("\n".join(lines) + "\n", native_settings)
    if path.exists():
        if path.read_text() == text:
            return path
        if "NEB-TS" in keyword_line.upper():
            from .orca_input_repair import prepare_neb_preopt_input_repair
            return prepare_neb_preopt_input_repair(path, text)
        raise RuntimeError(f"Conflicting retained ORCA input: {path}")
    with path.open("x") as stream:
        stream.write(text)
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
    prepare_fragments: bool | None = None,
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
            + ("\n  PreOpt true" if preopt_ends else "")
            + (f"\n  Prepare_Frags {'true' if prepare_fragments else 'false'}"
               if prepare_fragments is not None else "")
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


def create_orca_guided_scan_input(
    reactant_xyz: Path,
    *,
    coordinates: list[dict],
    charge: int,
    multiplicity: int,
    label: str = "guided-scan",
    ncores: int = 1,
    steps: int = 12,
    method_keywords: list[str] | None = None,
) -> Path:
    """Create a bounded simultaneous one- or two-bond relaxed scan."""
    if steps < 3:
        raise ValueError("A guided relaxed scan needs at least three points")
    if not 1 <= len(coordinates) <= 2:
        raise ValueError("The initial guided scan supports one or two coordinates")
    scan_lines = ["%geom", "  FullScan true"]
    for coordinate in coordinates:
        atom_indices = coordinate.get("atom_indices") or ()
        if len(atom_indices) != 2:
            raise ValueError("A guided bond coordinate needs two atom indices")
        left, right = (int(atom_indices[0]), int(atom_indices[1]))
        start = float(coordinate["start_distance_angstrom"])
        end = float(coordinate["end_distance_angstrom"])
        if min(left, right) < 0 or left == right or min(start, end) <= 0.0:
            raise ValueError("Invalid guided bond-scan coordinate")
        scan_lines.append(
            f"  Scan B {left} {right} = {start:.8f}, {end:.8f}, {int(steps)} end"
        )
    if len(coordinates) > 1:
        scan_lines.append("  Simul_Scan true")
    scan_lines.append("end")
    return _write_orca_xyzfile_job(
        reactant_xyz,
        label=label,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        keyword_line=_orca_method_line(method_keywords, "TightSCF", "Opt"),
        blocks=["\n".join(scan_lines)],
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


def create_orca_optimization_recovery_input(
    xyz_file: Path,
    *,
    charge: int,
    multiplicity: int,
    label: str = "opt",
    ncores: int = 1,
    max_iterations: int = 300,
    coordinate_system: str = "redundant",
    calculate_hessian: bool = True,
    recalculate_hessian_every: int | None = None,
    trust_radius: float | None = None,
    method_keywords: list[str] | None = None,
) -> Path:
    """Create one explicit, append-only minimum-optimization recovery job.

    The options are deliberately limited to documented ORCA optimizer
    controls.  They are scientific inputs and therefore belong in the caller's
    retained recovery contract and cache identity.
    """
    if max_iterations < 1:
        raise ValueError("max_iterations must be positive")
    normalized_coordinate_system = str(coordinate_system).strip().casefold()
    if normalized_coordinate_system not in {"redundant", "cartesian"}:
        raise ValueError("coordinate_system must be redundant or cartesian")
    if recalculate_hessian_every is not None and recalculate_hessian_every < 1:
        raise ValueError("recalculate_hessian_every must be positive when supplied")
    if trust_radius is not None and not 0.0 < float(trust_radius) <= 1.0:
        raise ValueError("trust_radius must be in (0, 1] when supplied")
    lines = [
        "%geom",
        f"  MaxIter {int(max_iterations)}",
        f"  CoordSys {normalized_coordinate_system}",
        "  CartFallback true",
    ]
    if calculate_hessian:
        lines.append("  Calc_Hess true")
    if recalculate_hessian_every is not None:
        lines.append(f"  Recalc_Hess {int(recalculate_hessian_every)}")
    if trust_radius is not None:
        lines.append(f"  Trust {float(trust_radius):.8f}")
    lines.append("end")
    return _write_orca_xyzfile_job(
        xyz_file,
        label=label,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        keyword_line=_orca_method_line(method_keywords, "TightSCF", "Opt"),
        blocks=["\n".join(lines)],
    )


def create_orca_environment_refinement_input(
    xyz_file: Path,
    *,
    interactions: list[dict],
    charge: int,
    multiplicity: int,
    label: str = "environment-refine",
    ncores: int = 1,
    max_iterations: int = 120,
    constrain_angles: bool = True,
    method_keywords: list[str] | None = None,
) -> Path:
    """Create a DFT optimization that preserves explicit H-bond coordinates.

    Constraints define the environment during optimization only.  Frequencies
    must be evaluated in a separate input without this ``%geom`` block.
    """
    if max_iterations < 1 or not interactions:
        raise ValueError("Environment refinement needs interactions and positive MaxIter")
    constraints = []
    for interaction in interactions:
        donor = int(interaction["donor_atom"])
        hydrogen = int(interaction["donor_hydrogen"])
        acceptor = int(interaction["acceptor_atom"])
        distance = float(interaction["h_bond_distance_angstrom"])
        angle = float(interaction["donor_h_acceptor_angle_degrees"])
        if min(donor, hydrogen, acceptor) < 0 or not 0.5 < distance < 5.0 or not 30.0 < angle <= 180.0:
            raise ValueError("Invalid environment-refinement interaction geometry")
        constraints.append(f"    {{ B {hydrogen} {acceptor} {distance:.8f} C }}")
        if constrain_angles:
            # ORCA's redundant-coordinate optimizer treats near-linear bends
            # specially and only permits freezing their current value.
            constraints.append(
                f"    {{ A {donor} {hydrogen} {acceptor} C }}"
                if angle >= 170.0 else
                f"    {{ A {donor} {hydrogen} {acceptor} {angle:.8f} C }}"
            )
    geom = "\n".join([
        "%geom", f"  MaxIter {max_iterations}", "  Constraints",
        *constraints, "  end", "end",
    ])
    return _write_orca_xyzfile_job(
        xyz_file, label=label, charge=charge, multiplicity=multiplicity, ncores=ncores,
        keyword_line=_orca_method_line(method_keywords, "TightSCF", "Opt"), blocks=[geom],
    )


def create_orca_gradient_input(
    xyz_file: Path, *, charge: int, multiplicity: int, label: str = "gradient",
    ncores: int = 1, method_keywords: list[str] | None = None,
    point_charge_file: Path | None = None,
) -> Path:
    """Create an unrestrained analytic-gradient job for a retained geometry."""
    blocks = None
    if point_charge_file is not None:
        point_charge_file = Path(point_charge_file)
        if not point_charge_file.is_file():
            raise FileNotFoundError(f"Point-charge file not found: {point_charge_file}")
        if point_charge_file.parent.resolve() != Path(xyz_file).parent.resolve():
            raise ValueError("Gradient geometry and point-charge file must share one directory")
        blocks = [f'%pointcharges "{point_charge_file.name}"']
    return _write_orca_xyzfile_job(
        xyz_file, label=label, charge=charge, multiplicity=multiplicity, ncores=ncores,
        keyword_line=_orca_method_line(method_keywords, "TightSCF", "EnGrad"),
        blocks=blocks,
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
    seed_mechanisms: list[Path] | None = None,
    database_reaction_libraries: list[str] | None = None,
    cap_generated_carbon_at_target: bool = True,
    maximum_heavy_atoms: int | None = None,
    maximum_radical_electrons: int | None = None,
    filter_reactions: bool = False,
    max_objects_per_iteration: int = 1,
    restart_from_seed: Path | None = None,
    simulator_atol: float = 1e-16,
    simulator_rtol: float = 1e-8,
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
        seed_mechanisms: Local mechanisms to place directly in the RMG core.
            Scientific admission is the caller's responsibility.
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
    seed_mechanisms = list(dict.fromkeys(Path(path).resolve() for path in (seed_mechanisms or [])))
    database_reaction_libraries = list(database_reaction_libraries or [])
    if any(not str(name).strip() for name in database_reaction_libraries):
        raise ValueError("RMG database reaction-library names must be nonempty")
    database_reaction_libraries = list(dict.fromkeys(str(name).strip() for name in database_reaction_libraries))
    for library in [*reaction_libraries, *seed_mechanisms]:
        library = Path(library)
        if not (library / "reactions.py").is_file() or not (library / "dictionary.txt").is_file():
            raise ValueError(f"RMG reaction library must contain reactions.py and dictionary.txt: {library}")
    if seed_mechanisms:
        # RMG registers an external library by basename. Distinct paths with
        # the same registered key can otherwise silently replace a core seed.
        registered_sources = {}
        for entry in [*seed_mechanisms, *reaction_libraries, *database_reaction_libraries]:
            path = Path(entry)
            source = str(path.resolve()) if path.is_absolute() or path.is_dir() else str(entry)
            key = Path(source).name if Path(source).is_absolute() else source
            if key in registered_sources and registered_sources[key] != source:
                raise ValueError(f"RMG library key collision for seeded mechanism: {key}")
            registered_sources[key] = source
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
    if not math.isfinite(float(simulator_atol)) or float(simulator_atol) <= 0.0:
        raise ValueError("simulator_atol must be positive and finite")
    if not math.isfinite(float(simulator_rtol)) or not 0.0 < float(simulator_rtol) < 1.0:
        raise ValueError("simulator_rtol must be finite and between zero and one")
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
        *(str(Path(library).resolve()) for library in reaction_libraries
          if Path(library).resolve() not in seed_mechanisms),
        *database_reaction_libraries,
    ]

    content = f"""database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries={library_entries!r},
    seedMechanisms={[str(path) for path in seed_mechanisms]!r},
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
    atol={float(simulator_atol)!r},
    rtol={float(simulator_rtol)!r},
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
