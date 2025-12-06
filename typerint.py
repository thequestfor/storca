#!/usr/bin/env python3
from pathlib import Path
import typer
from rich.prompt import Prompt

from src.inputgen import create_orca_input
from src.orca_runner import run_orca
from src.parser import parse_orca_orbitals
from src.converter import gbw_to_wfn
from src.wfngrid import run_multiwfn_fukui
from src.quickprops import analyze_xyz

app = typer.Typer(help="Automated ORCA + Multiwfn workflow with quick property estimates")


def display_descriptors(info):
    typer.echo("\n=== Quick Descriptors ===")
    descriptors = {k: info[k] for k in ["SMILES", "IUPACName", "MolecularWeight", "TPSA", "XLogP"]}
    for k, v in descriptors.items():
        typer.echo(f"  {k}: {v}")


def display_pubchem(info):
    typer.echo("\n=== PubChem / Additional Info ===")
    properties = {k: info[k] for k in ["CID", "Synonyms", "Description"]}
    for k, v in properties.items():
        if isinstance(v, list):
            v = ", ".join(v) if v else "N/A"
        typer.echo(f"  {k}: {v if v else 'N/A'}")


def run_orca_workflow(xyz_file: Path, base_name: str):
    typer.echo("\n=== Running ORCA + Multiwfn Workflow ===")

    # Step 1: Geometry Optimization
    opt_label = "opt"
    opt_inp = create_orca_input(xyz_file, charge=0, multiplicity=1, opt=True, label=opt_label)
    typer.echo(f"Created ORCA input for optimization: {opt_inp}")

    opt_outputs = run_orca(opt_inp)
    optimized_xyz = opt_outputs["xyz"]
    new_xyz_name = optimized_xyz.with_name(f"{base_name}_{opt_label}.xyz")
    optimized_xyz.rename(new_xyz_name)
    optimized_xyz = new_xyz_name

    typer.echo("\nGeometry optimization finished. Files generated:")
    for key, path in opt_outputs.items():
        if key == "xyz":
            path = optimized_xyz
        typer.echo(f"  {key}: {path}")

    # Step 2: Single-Point Calculations
    charges = [0, +1, -1]
    sp_gbws = {}
    for charge in charges:
        mult = 1 if charge == 0 else 2
        label = f"sp{charge:+}"
        sp_inp = create_orca_input(optimized_xyz, charge=charge, multiplicity=mult, opt=False, label=label)
        typer.echo(f"\nRunning single-point calculation for charge {charge} (multiplicity {mult}): {sp_inp}")
        sp_outputs = run_orca(sp_inp)
        sp_gbws[charge] = sp_outputs["gbw"]
        typer.echo(f"  GBW generated: {sp_outputs['gbw']}")

    # Step 3: Parse HOMO/LUMO
    neutral_label = f"sp{0:+}"
    neutral_out = optimized_xyz.with_name(f"{base_name}_opt_{neutral_label}.out")
    orbitals = parse_orca_orbitals(neutral_out)
    typer.echo("\nParsed orbitals for neutral species:")
    typer.echo(f"  HOMO: {orbitals['homo_number']}  energy: {orbitals['homo_energy']} eV")
    typer.echo(f"  LUMO: {orbitals['lumo_number']}  energy: {orbitals['lumo_energy']} eV")

    # Step 4: Convert GBW → WFN
    wfn_files = {}
    for charge, gbw in sp_gbws.items():
        wfn = gbw_to_wfn(gbw)
        wfn_files[charge] = wfn
        typer.echo(f"GBW → WFN conversion done for charge {charge}: {wfn}")

    # Step 5: Rename WFNs
    new_names = {0: Path("N.wfn"), 1: Path("N+1.wfn"), -1: Path("N-1.wfn")}
    for charge, old_path in wfn_files.items():
        new_path = new_names[charge]
        old_path.rename(new_path)
        wfn_files[charge] = new_path

    typer.echo("\nRenamed WFN files:")
    for charge, path in wfn_files.items():
        typer.echo(f"  Charge {charge}: {path}")


@app.command()
def run(xyz_file: Path = typer.Argument(..., help="Input XYZ file")):
    if not xyz_file.exists():
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")
    base_name = xyz_file.stem

    # Step 0: Quick properties
    try:
        info = analyze_xyz(xyz_file)
    except Exception as e:
        typer.echo(f"Failed to estimate properties: {e}")
        info = None

    # Main interactive menu
    while True:
        choices = ["descriptors", "pubchem", "ORCA", "exit"]
        choice = Prompt.ask("\nSelect a section to view", choices=choices)

        if choice == "descriptors":
            if info:
                display_descriptors(info)
            else:
                typer.echo("No descriptor info available.")
        elif choice == "pubchem":
            if info:
                display_pubchem(info)
                typer.echo("\nDescription:\n" + info["Description"])
            else:
                typer.echo("No PubChem info available.")
        elif choice == "ORCA":
            run_orca_workflow(xyz_file, base_name)
        elif choice == "exit":
            typer.echo("Exiting menu.")
            break


if __name__ == "__main__":
    app()

