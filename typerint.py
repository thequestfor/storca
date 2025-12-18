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
from src.molecule_tools import xyz_to_smiles, smiles_to_xyz, smiles_to_png, xyz_to_png

app = typer.Typer(
    help="Automated ORCA + Multiwfn workflow with quick property and hazard estimates"
)


# ===========================
# DISPLAY: DESCRIPTORS + HAZARDS
# ===========================
def display_descriptors(info):
    typer.echo("\n=== Quick Descriptors ===")

    descriptors = [
        ("SMILES", info.get("SMILES")),
        ("IUPAC Name", info.get("IUPACName")),
        ("Molecular Weight", info.get("MolecularWeight")),
        ("TPSA", info.get("TPSA")),
        ("XLogP", info.get("XLogP")),
    ]

    for key, value in descriptors:
        typer.echo(f"  {key}: {value if value is not None else 'N/A'}")

    hazards_info = info.get("Hazards")
    if hazards_info:
        typer.echo("\n=== Hazard Summary ===")
        typer.echo(f"Source: {hazards_info.get('Source')}")
        typer.echo(f"Signal Word: {hazards_info.get('SignalWord')}")

        ghs = hazards_info.get("GHS", {})
        codes = ghs.get("Codes") or []
        descs = ghs.get("Descriptions") or {}

        typer.echo("GHS Hazard Statements:")
        if codes:
            for code in codes:
                meaning = descs.get(code) if isinstance(descs, dict) else "See SDS"
                typer.echo(f"  {code}: {meaning}")
        else:
            typer.echo("  None")

        pictos = ghs.get("Pictograms") or []
        if pictos:
            typer.echo(f"Pictograms: {', '.join(pictos)}")

        practical = hazards_info.get("Practical", {})
        typer.echo("\n=== Practical Lab Guidance ===")
        typer.echo(f"Flammable: {practical.get('Flammable')}")
        typer.echo(f"Toxicity: {practical.get('Toxicity')}")
        typer.echo(f"Environmental Hazard: {practical.get('EnvironmentalHazard')}")
        typer.echo(f"Advice: {practical.get('Advice')}")


# ===========================
# DISPLAY: PUBCHEM / DESCRIPTION
# ===========================
def display_pubchem(info):
    typer.echo("\n=== Compound Information ===")
    typer.echo(f"  CID: {info.get('CID') or 'N/A'}")

    synonyms = info.get("Synonyms", [])
    if synonyms:
        typer.echo(f"  Synonyms: {', '.join(synonyms[:10])}")
    else:
        typer.echo("  Synonyms: N/A")

    typer.echo("\n=== Description ===")
    typer.echo(info.get("Description") or "No description available.")


# ===========================
# ORCA WORKFLOW (UNCHANGED)
# ===========================
def run_orca_workflow(xyz_file: Path, base_name: str):
    typer.echo("\n=== Running ORCA + Multiwfn Workflow ===")

    opt_label = "opt"
    opt_inp = create_orca_input(
        xyz_file, charge=0, multiplicity=1, opt=True, label=opt_label
    )
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

    charges = [0, +1, -1]
    sp_gbws = {}

    for charge in charges:
        mult = 1 if charge == 0 else 2
        label = f"sp{charge:+}"
        sp_inp = create_orca_input(
            optimized_xyz, charge=charge, multiplicity=mult, opt=False, label=label
        )
        typer.echo(
            f"\nRunning single-point calculation for charge {charge} (multiplicity {mult}): {sp_inp}"
        )
        sp_outputs = run_orca(sp_inp)
        sp_gbws[charge] = sp_outputs["gbw"]
        typer.echo(f"  GBW generated: {sp_outputs['gbw']}")

    neutral_out = optimized_xyz.with_name(f"{base_name}_opt_sp+0.out")
    orbitals = parse_orca_orbitals(neutral_out)

    typer.echo("\nParsed orbitals for neutral species:")
    typer.echo(
        f"  HOMO: {orbitals['homo_number']}  energy: {orbitals['homo_energy']} eV"
    )
    typer.echo(
        f"  LUMO: {orbitals['lumo_number']}  energy: {orbitals['lumo_energy']} eV"
    )

    wfn_files = {}
    for charge, gbw in sp_gbws.items():
        wfn = gbw_to_wfn(gbw)
        wfn_files[charge] = wfn
        typer.echo(f"GBW → WFN conversion done for charge {charge}: {wfn}")

    new_names = {0: Path("N.wfn"), 1: Path("N+1.wfn"), -1: Path("N-1.wfn")}
    for charge, old_path in wfn_files.items():
        new_path = new_names[charge]
        old_path.rename(new_path)
        wfn_files[charge] = new_path

    typer.echo("\nRenamed WFN files:")
    for charge, path in wfn_files.items():
        typer.echo(f"  Charge {charge}: {path}")


# ===========================
# MAIN COMMAND (XYZ or SMILES)
# ===========================
@app.command()
def run(
    xyz_file: Path = typer.Option(None, help="Input XYZ file"),
    smiles: str = typer.Option(None, help="Input SMILES string"),
    generate_png: bool = typer.Option(True, help="Generate 2D PNG of molecule")
):
    """Analyze molecule from either XYZ or SMILES and optionally generate PNG."""

    if not xyz_file and not smiles:
        raise typer.BadParameter("You must provide either --xyz or --smiles input.")

    # Convert XYZ → SMILES if needed
    if xyz_file:
        if not xyz_file.exists():
            raise FileNotFoundError(f"XYZ file not found: {xyz_file}")
        smiles_value = xyz_to_smiles(xyz_file)
        info = analyze_xyz(xyz_file)
        base_name = xyz_file.stem
    else:
        smiles_value = smiles
        base_name = "from_smiles"
        # generate temporary XYZ for analysis
        xyz_file = smiles_to_xyz(smiles_value, Path(f"{base_name}.xyz"))
        info = analyze_xyz(xyz_file)

    # Optional PNG
    if generate_png:
        png_file = Path(f"{base_name}.png")
        smiles_to_png(smiles_value, png_file)
        typer.echo(f"2D PNG saved at: {png_file}")

    # Interactive display
    while True:
        choice = Prompt.ask(
            "\nSelect a section to view",
            choices=["descriptors", "pubchem", "ORCA", "exit"],
        )

        if choice == "descriptors":
            if info:
                display_descriptors(info)
            else:
                typer.echo("No descriptor info available.")
        elif choice == "pubchem":
            if info:
                display_pubchem(info)
            else:
                typer.echo("No compound information available.")
        elif choice == "ORCA":
            run_orca_workflow(xyz_file, base_name)
        elif choice == "exit":
            typer.echo("Exiting.")
            break


if __name__ == "__main__":
    app()

