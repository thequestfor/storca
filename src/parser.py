from pathlib import Path
from typing import List, Tuple
from collections import defaultdict
import math
import os
import re

def parse_orca_orbitals(out_file: Path) -> dict:
    orbitals = []
    reading = False
    table_started = False

    with open(out_file) as f:
        for line in f:
            line_strip = line.strip()
            if "ORBITAL ENERGIES" in line_strip:
                reading = True
                continue
            if reading:
                if set(line_strip) == {"-"}:
                    continue  # skip separator
                if not table_started:
                    if "NO" in line_strip and "OCC" in line_strip:
                        table_started = True
                    continue
                # stop at empty line or notes
                if line_strip == "" or "Only the first" in line_strip:
                    break
                # parse orbital line
                parts = line_strip.split()
                if len(parts) >= 4:
                    try:
                        no = int(parts[0])
                        occ = float(parts[1])
                        eh = float(parts[2])
                        ev = float(parts[3])
                        orbitals.append({"no": no, "occ": occ, "eh": eh, "ev": ev})
                    except ValueError:
                        continue

    if not orbitals:
        return {"homo_number": None, "homo_energy": None,
                "lumo_number": None, "lumo_energy": None}

    occupied = [o for o in orbitals if o["occ"] > 0]
    virtual = [o for o in orbitals if o["occ"] == 0]
    homo_orb = max(occupied, key=lambda x: x["no"]) if occupied else None
    lumo_orb = min(virtual, key=lambda x: x["no"]) if virtual else None

    return {
        "homo_number": homo_orb["no"] if homo_orb else None,
        "homo_energy": homo_orb["ev"] if homo_orb else None,
        "lumo_number": lumo_orb["no"] if lumo_orb else None,
        "lumo_energy": lumo_orb["ev"] if lumo_orb else None,
    }

def parse_orca_energy(out_file: Path) -> float:
    """
    Extract the final ORCA energy in Hartree.
    """
    energy = None
    with open(out_file) as f:
        for line in f:
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    energy = float(line.split()[-1])
                    continue
                except ValueError:
                    continue
            if "Energy" in line and "Eh" in line:
                try:
                    # Split on ':' and take the numeric part before 'Eh'
                    energy_str = line.split(':')[1].split()[0]
                    energy = float(energy_str)
                except Exception:
                    continue
    if energy is None:
        raise RuntimeError(f"Could not find SCF energy in {out_file}")
    return energy

# --------------------------
# GOAT OUTPUT PARSER
# --------------------------
def parse_goat_out(out_file):
    """
    Parse GOAT .out file for conformers and Boltzmann weighings.
    Returns:
        conformers: list of all conformers with their weights
        sig_conformers: list of conformers with weight > 0.05
    """
    conformers = []

    with open(out_file, "r") as f:
        lines = f.readlines()

    # Find the start of the Boltzmann chart
    chart_start = None
    for i, line in enumerate(lines):
        if re.match(r'\s*Conformer\s+Energy\s+Degen\.', line):
            chart_start = i + 2  # skip header/dashed line
            break

    if chart_start is None:
        raise ValueError("Boltzmann chart not found in .out file")

    # Parse the chart lines
    for line in lines[chart_start:]:
        if line.strip() == "":
            break
        parts = line.split()
        if len(parts) >= 5:
            idx = int(parts[0])
            energy = float(parts[1])
            degeneracy = int(parts[2])
            weight = float(parts[3]) / 100.0  # convert percent -> fraction
            cum_weight = float(parts[4]) / 100.0
            conformers.append({
                "idx": idx,
                "energy": energy,
                "degeneracy": degeneracy,
                "weight": weight,
                "cum_weight": cum_weight
            })

    # Significant conformers: weight > 0.05
    sig_conformers = [c for c in conformers if c["weight"] > 0.05]

    return conformers, sig_conformers


# --------------------------
# XYZ PARSER
# --------------------------
def parse_xyz_ensemble(xyz_file: Path):
    """
    Parse an XYZ file containing multiple frames.
    Returns a list of frames:
        [{'natoms': int, 'comment': str, 'xyz_lines': list[str]}]
    """
    frames = []
    with open(xyz_file, "r") as f:
        lines = [line.rstrip("\n") for line in f]

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            natoms = int(line)
        except ValueError:
            raise ValueError(f"Expected number of atoms at line {i+1}, got: {line}")

        if i + 1 + natoms >= len(lines):
            raise ValueError(f"Incomplete frame starting at line {i+1}")

        comment = lines[i + 1]
        xyz_lines = lines[i + 2 : i + 2 + natoms]
        frames.append({
            "natoms": natoms,
            "comment": comment,
            "xyz_lines": xyz_lines
        })
        i += 2 + natoms

    return frames



# --------------------------
# HELPER: COMBINE OUT & XYZ
# --------------------------
def attach_xyz_to_conformers(conformers, xyz_frames):
    """
    Attach XYZ frames to parsed conformers.
    Assumes frames are in the same order as conformers.
    """
    if len(conformers) != len(xyz_frames):
        raise ValueError(
            f"Mismatch: {len(conformers)} conformers vs {len(xyz_frames)} xyz frames"
        )

    for c, frame in zip(conformers, xyz_frames):
        c["xyz"] = frame

    return conformers



def merge_weighted_ir(
    conformers,
    freq_tol: float = 10.0
):
    """
    Merge Boltzmann-weighted IR modes into a single stick spectrum.

    Returns:
        List of (freq_cm1, intensity_km_mol)
    """
    all_modes = []

    for conf in conformers:
        if conf["freqs"] is None or conf["ir_intensities"] is None:
            continue
        w = conf["weight"]
        for f, I in zip(conf["freqs"], conf["ir_intensities"]):
            all_modes.append((f, I * w))

    # Sort by frequency
    all_modes.sort(key=lambda x: x[0])

    merged = []
    current = []

    for f, I in all_modes:
        if not current:
            current = [(f, I)]
            continue
        if abs(f - current[-1][0]) <= freq_tol:
            current.append((f, I))
        else:
            merged.append(current)
            current = [(f, I)]

    if current:
        merged.append(current)

    # Weighted average frequency, summed intensity
    final_spectrum = []
    for group in merged:
        total_I = sum(I for _, I in group)
        if total_I <= 0.0:
            continue
        avg_f = sum(f * I for f, I in group) / total_I
        final_spectrum.append((avg_f, total_I))

    return final_spectrum

from pathlib import Path

def parse_orca_ir(file_path: str | Path):
    """
    Parse ORCA IR output file for frequencies, epsilon, and intensity.

    Returns:
        List of dicts:
            [{'mode': int, 'freq': float, 'eps': float, 'intensity': float}, ...]
    """
    vibrations = []

    with open(file_path) as f:
        lines = f.readlines()

    # Find the line that says "IR SPECTRUM"
    for i, line in enumerate(lines):
        if "IR SPECTRUM" in line:
            start_data = i + 1
            break
    else:
        raise RuntimeError("IR SPECTRUM section not found in ORCA output")

    # Parse table. Header layout differs between ORCA versions, so only rows
    # whose first four fields are numeric are accepted. Blank lines before the
    # first mode are header spacing; blank lines after modes end the table.
    table_started = False
    for line in lines[start_data:]:
        line = line.strip()
        if not line or line.startswith("*"):
            if not table_started:
                continue
            # Stop at a blank line or footer after at least one data row.
            break

        # Split tokens
        tokens = line.split()
        if len(tokens) < 4:
            # Not enough columns; skip
            continue

        # Extract data
        try:
            mode = int(tokens[0].rstrip(":"))
            freq = float(tokens[1].replace("D", "E").replace("d", "e"))
            eps = float(tokens[2].replace("D", "E").replace("d", "e"))
            intensity = float(tokens[3].replace("D", "E").replace("d", "e"))
        except ValueError:
            continue

        if (mode < 0 or not all(math.isfinite(value) for value in (freq, eps, intensity))
                or intensity < 0):
            continue

        table_started = True

        vibrations.append({
            "mode": mode,
            "freq": freq,
            "eps": eps,
            "intensity": intensity
        })

    if not vibrations:
        raise RuntimeError("IR SPECTRUM section contained no parseable finite modes")
    return vibrations


def classify_chemkin_route(reaction_equation: str, label: str = "stability") -> dict:
    """Classify the reactant context of an RMG/Chemkin route.

    The classification concerns the forward direction as written by RMG.  It
    describes whether the target can react alone or requires another species;
    it is not proof that a co-reactant is present in a real storage system.
    """
    try:
        lhs, _ = re.split(r"<=>|=>", reaction_equation, maxsplit=1)
    except ValueError as error:
        raise ValueError(f"Chemkin reaction has no supported arrow: {reaction_equation}") from error

    reactants = []
    for item in lhs.split("+"):
        token = item.strip()
        if not token:
            continue
        coefficient_match = re.match(r"^(?P<coefficient>\d+(?:\.\d+)?)\s+(?P<label>.+)$", token)
        coefficient = float(coefficient_match.group("coefficient")) if coefficient_match else 1.0
        species_label = coefficient_match.group("label") if coefficient_match else token
        # RMG Chemkin labels species with a numeric index, e.g. stability(1).
        species_label = re.sub(r"\(\d+\)$", "", species_label).strip()
        reactants.append({"label": species_label, "stoichiometry": coefficient})

    target_stoichiometry = sum(item["stoichiometry"] for item in reactants if item["label"] == label)
    co_reactants = [item for item in reactants if item["label"] != label]
    radical_co_reactants = [item["label"] for item in co_reactants if re.search(r"\[[^\]]+\]", item["label"])]
    if target_stoichiometry <= 0:
        raise ValueError(f"Target label '{label}' is not a reactant in: {reaction_equation}")
    if not co_reactants and target_stoichiometry == 1:
        context = "unimolecular_decomposition"
    elif not co_reactants and target_stoichiometry > 1:
        context = "self_reaction"
    elif radical_co_reactants:
        context = "radical_or_impurity_initiated"
    else:
        context = "co_reactant_dependent"
    return {
        "route_context": context,
        "target_stoichiometry": target_stoichiometry,
        "co_reactants": co_reactants,
        "radical_co_reactants": radical_co_reactants,
    }


def parse_chemkin_annotated(
    chemkin_file: Path,
    barrier_threshold: float = 50.0,
    label: str = "stability",
) -> dict:
    """
    Determine kinetic stability from chem_annotated.inp.
    
    A species is considered unstable if it appears as a reactant
    in any decomposition reaction with barrier below the threshold.
    The investigated species can be selected with *label*.
    """
    text = chemkin_file.read_text()

    # Chemkin's third Arrhenius coefficient uses the energy unit declared on
    # the REACTIONS line.  RMG normally writes this declaration, and silently
    # treating an undeclared value as kcal/mol would make the screen unsafe.
    energy_units = {
        "KCAL/MOLE": 1.0,
        "CAL/MOLE": 0.001,
        "KJOULES/MOLE": 0.239005736,
        "JOULES/MOLE": 0.000239005736,
    }
    header = re.search(r"^\s*REACTIONS\b(?P<units>[^\n]*)", text, re.IGNORECASE | re.MULTILINE)
    declared_units = header.group("units").upper().replace(" ", "") if header else ""
    source_unit = next((unit for unit in energy_units if unit in declared_units), None)

    decomposition_reactions = []
    reaction_barriers = []
    candidate_routes = []

    reactions_section = re.search(r"REACTIONS.*?END", text, re.DOTALL | re.IGNORECASE)
    if not reactions_section:
        return {
            "stable": True,
            "decomposition_reactions": [],
            "reaction_barriers": [],
            "candidate_routes": [],
            "activation_energy_source_unit": source_unit,
            "activation_energy_unit": "kcal/mol",
        }
    if source_unit is None:
        raise ValueError(
            "Chemkin REACTIONS header must declare activation-energy units "
            "(for example KCAL/MOLE or CAL/MOLE)"
        )

    for line in reactions_section.group(0).splitlines():
        line = line.strip()
        if not line or line.startswith("!") or "<=>" not in line:
            continue

        match = re.match(
            r"^(?P<reaction>.+?(?:<=>|=>).+?)\s+"
            r"(?P<a>[-+0-9.Ee]+)\s+(?P<n>[-+0-9.Ee]+)\s+(?P<ea>[-+0-9.Ee]+)"
            r"(?:\s*!.*)?$",
            line,
        )
        if not match:
            continue
        reaction_part = match.group("reaction")
        Ea = float(match.group("ea")) * energy_units[source_unit]
        lhs, rhs = re.split(r"<=>|=>", reaction_part, maxsplit=1)
        clean_species = lambda item: re.sub(r"^\d+(?:\.\d+)?\s*|\(\d+\)", "", item).strip()
        lhs_species = [clean_species(species) for species in lhs.split("+")]
        rhs_species = [clean_species(species) for species in rhs.split("+")]

        # Chemkin may write a reversible reaction in either direction.  The
        # Arrhenius parameters belong only to the printed forward direction;
        # never reuse them as a reverse decomposition barrier.
        if label in lhs_species:
            if rhs_species == [label]:
                continue
            if Ea >= barrier_threshold:
                continue
            decomposition_reactions.append(line)
            reaction_barriers.append(Ea)
            candidate_routes.append({
                "reaction": line,
                "reaction_equation": reaction_part,
                "activation_energy_kcal_mol": Ea,
                "direction": "chemkin_forward",
                "kinetics_representation": "printed_arrhenius",
                **classify_chemkin_route(reaction_part, label=label),
            })
        elif "<=>" in reaction_part and label in rhs_species:
            reverse_equation = f"{rhs.strip()}<=>{lhs.strip()}"
            candidate_routes.append({
                "reaction": line,
                "reaction_equation": reverse_equation,
                "printed_reaction_equation": reaction_part,
                "activation_energy_kcal_mol": None,
                "direction": "reverse_of_chemkin_direction",
                "kinetics_representation": "reversible_rate_requires_thermodynamic_reverse",
                **classify_chemkin_route(reverse_equation, label=label),
            })

    stable = len(decomposition_reactions) == 0

    return {
        "stable": stable,
        "decomposition_reactions": decomposition_reactions,
        "reaction_barriers": reaction_barriers,
        "candidate_routes": candidate_routes,
        "activation_energy_source_unit": source_unit,
        "activation_energy_unit": "kcal/mol",
    }
