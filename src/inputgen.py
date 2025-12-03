from pathlib import Path

def create_orca_input(xyz_file: Path, charge=0, multiplicity=1, opt=True, label: str = ""):
    """
    Create an ORCA input file from an XYZ file.
    label: optional string to append to the base name (e.g., 'opt', 'sp0', 'sp+1')
    """
    method_line = "! B3LYP def2-SVP"
    method_line += " Opt TightSCF" if opt else " TightSCF"

    inp_content = f"""{method_line}

* xyzfile {charge} {multiplicity} {xyz_file}
"""

    # Build filename with label
    if label:
        new_name = f"{xyz_file.stem}_{label}.inp"
    else:
        new_name = f"{xyz_file.stem}.inp"

    inp_file = xyz_file.with_name(new_name)

    with open(inp_file, "w") as f:
        f.write(inp_content)

    return inp_file

