# src/inputgen.py
from pathlib import Path
from rdkit import Chem


def create_orca_input(xyz_file: Path, charge=0, multiplicity=1, opt=True, label: str = "", intel: dict = None):
    """
    Smarter ORCA input generator using structure intelligence.
    intel: optional dict containing info from analyze_molecule (e.g., HasMetal, NumAtoms).
    """

    # Default method
    method_line = "! B3LYP def2-SVP"

    # Adjust for molecule size or metals
    if intel:
        num_atoms = intel.get("NumAtoms", 0)
        has_metal = intel.get("HasMetal", False)

        if has_metal:
            method_line = "! B3LYP def2-TZVP ECP"
        elif num_atoms > 50:
            method_line = "! RI-B3LYP def2-SVP TightSCF"

    # Optimization flag
    method_line += " Opt TightSCF" if opt else " TightSCF"

    # Build input content
    inp_content = f"""{method_line}

* xyzfile {charge} {multiplicity} {xyz_file}
"""

    # Build filename with label
    if label:
        new_name = f"{xyz_file.stem}_{label}.inp"
    else:
        new_name = f"{xyz_file.stem}.inp"

    inp_file = xyz_file.with_name(new_name)
    inp_file.write_text(inp_content)

    return inp_file
