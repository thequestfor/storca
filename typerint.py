#!/usr/bin/env python3
from pathlib import Path
import typer
from rich.prompt import Prompt
from typing import Optional
import subprocess

from src.inputgen import create_orca_input
from src.orca_runner import run_orca
from src.parser import parse_orca_orbitals
from src.converter import gbw_to_wfn
from src.wfngrid import run_multiwfn_fukui
from src.quickprops import analyze_molecule
from src.molecule_tools import (
    smiles_to_xyz,
    smiles_to_png,
    xyz_to_png,
    sanitize_smiles,
)
from src.stability.freq_check import frequency_stability_check
from src.utils.orca_xyz_clean import clean_xyz_for_orca

app = typer.Typer(
    help="Automated ORCA + Multiwfn workflow with quick property and hazard estimates"
)

# --- display functions ---
def display_descriptors(info):
    typer.echo("\n=== Quick Descriptors ===")
    for key in ["SMILES", "IUPACName", "MolecularWeight", "TPSA", "XLogP"]:
        typer.echo(f"  {key}: {info.get(key) or 'N/A'}")

    hazards = info.get("Hazards")
    if hazards:
        typer.echo("\n=== Hazards ===")
        typer.echo(f"Source: {hazards.get('Source')}")
        typer.echo(f"Signal Word: {hazards.get('SignalWord')}")
        ghs = hazards.get("GHS", {})
        codes = ghs.get("Codes") or []
        descs = ghs.get("Descriptions") or {}
        if codes:
            for code in codes:
                typer.echo(f"  {code}: {descs.get(code) or 'See SDS'}")
        else:
            typer.echo("  None")

        practical = hazards.get("Practical", {})
        if practical:
            typer.echo("\n=== Practical Lab Advice ===")
            for k in ["Flammable", "Toxicity", "EnvironmentalHazard", "Advice"]:
                typer.echo(f"{k}: {practical.get(k)}")

    intel = info.get("StructureIntelligence", {})
    if intel:
        typer.echo("\n=== Structure Intelligence ===")
        for k in ["CompoundClass", "HasMetal", "NumAtoms", "FormalCharge", "Multiplicity"]:
            if k in intel:
                typer.echo(f"  {k}: {intel.get(k)}")


def display_pubchem(info):
    typer.echo("\n=== Compound Info ===")
    typer.echo(f"CID: {info.get('CID') or 'N/A'}")
    synonyms = info.get("Synonyms") or []
    typer.echo(f"Synonyms: {', '.join(synonyms[:10]) if synonyms else 'N/A'}")
    source = info.get("DescriptionSource", "Unknown")
    typer.echo(f"\nDescription (source: {source}):")
    typer.echo(info.get("Description") or "No description available.")


def clean_xyz_lines(xyz_file: Path) -> Path:
    cleaned_xyz = xyz_file.with_name(f"{xyz_file.stem}_clean.xyz")
    with open(xyz_file) as f_in, open(cleaned_xyz, "w") as f_out:
        lines = f_in.readlines()
        start_idx = 0
        if len(lines) >= 2:
            try:
                int(lines[0].strip())
                start_idx = 2
            except ValueError:
                start_idx = 0
        f_out.writelines(lines[start_idx:])
    return cleaned_xyz


def check_imaginary_frequencies(freq_out: Path):
    imag_freqs = []
    with open(freq_out) as f:
        for line in f:
            if "Frequencies" in line:
                parts = line.split()[1:]
                for val in parts:
                    try:
                        freq = float(val)
                        if freq < 0.0:
                            imag_freqs.append(freq)
                    except ValueError:
                        continue
    if imag_freqs:
        raise RuntimeError(
            f"Imaginary frequencies detected: {imag_freqs}\n"
            "This compound does not exist as a stable minimum."
        )


def run_orca_workflow(xyz_file: Path, base_name: str, intel: dict = None, ncores: int = 1):
    """
    Full ORCA + Multiwfn workflow with checks for imaginary frequencies.
    """
    typer.echo("\n=== Running ORCA + Multiwfn Workflow ===")
    intel = intel or {}
    charge0 = intel.get("Charge", 0)
    mult0 = intel.get("Multiplicity", 1)

    # --- Geometry optimization ---
    opt_label = "opt"
    opt_inp = create_orca_input(
        xyz_file, charge=charge0, multiplicity=mult0, opt=True, label=opt_label, ncores=ncores
    )
    typer.echo(f"Created ORCA input for optimization: {opt_inp}")
    opt_outputs = run_orca(opt_inp)
    optimized_xyz = opt_outputs["xyz"]
    typer.echo("\nGeometry optimization finished. Files generated:")
    for key, path in opt_outputs.items():
        typer.echo(f"  {key}: {path}")

    # --- Frequency calculation ---
    freq_label = "freq"
    freq_inp = create_orca_input(
        optimized_xyz, charge=charge0, multiplicity=mult0, opt=False, freq=True, label=freq_label, ncores=ncores
    )
    typer.echo(f"\nRunning frequency calculation: {freq_inp}")
    freq_outputs = run_orca(freq_inp)
    typer.echo("\nFrequency calculation finished. Files generated:")
    for key, path in freq_outputs.items():
        typer.echo(f"  {key}: {path}")

    check_imaginary_frequencies(freq_outputs["out"])
    typer.echo("No imaginary frequencies detected. Compound is a true minimum.")

    # --- Single-point calculations ---
    charges = [charge0, charge0 + 1, charge0 - 1]
    sp_gbws = {}
    sp_outs = {}
    for charge in charges:
        multiplicity = 1 if charge == charge0 else 2
        label = f"sp{charge - charge0:+}"
        sp_inp = create_orca_input(
            optimized_xyz, charge=charge, multiplicity=multiplicity, opt=False, keepdens=True, label=label, ncores=ncores
        )
        typer.echo(f"\nRunning single-point calculation for charge {charge} (multiplicity {multiplicity}): {sp_inp}")
        sp_outputs = run_orca(sp_inp)
        sp_gbws[charge - charge0] = sp_outputs["gbw"]
        sp_outs[charge - charge0] = sp_outputs["out"]
        typer.echo(f"  GBW generated: {sp_outputs['gbw']}")

    # --- Parse orbitals ---
    orbitals = parse_orca_orbitals(sp_outs[0])
    typer.echo("\nParsed orbitals for neutral species:")
    typer.echo(f"  HOMO: {orbitals['homo_number']}  energy: {orbitals['homo_energy']} eV")
    typer.echo(f"  LUMO: {orbitals['lumo_number']}  energy: {orbitals['lumo_energy']} eV")

    # --- Convert GBW → WFN ---
    wfn_files = {}
    for charge_key, gbw_path in sp_gbws.items():
        wfn = gbw_to_wfn(gbw_path)
        wfn_files[charge_key] = wfn
        typer.echo(f"GBW → WFN conversion done for charge {charge_key}: {wfn}")

    # Rename to N.wfn, N+1.wfn, N-1.wfn
    new_names = {0: Path("N.wfn"), 1: Path("N+1.wfn"), -1: Path("N-1.wfn")}
    for charge_key, old_path in wfn_files.items():
        new_path = new_names[charge_key]
        old_path.rename(new_path)
        wfn_files[charge_key] = new_path

    typer.echo("\nRenamed WFN files:")
    for charge_key, path in wfn_files.items():
        typer.echo(f"  Charge {charge_key}: {path}")

    typer.echo("\nAll ORCA + Multiwfn steps completed successfully.")


@app.command()
def run(
    xyz_file: Optional[Path] = typer.Argument(None, help="Input XYZ file (if using SMILES, omit this)"),
    smiles_value: Optional[str] = typer.Option(None, "--smiles", help="SMILES string"),
    generate_png: bool = typer.Option(True, help="Generate PNG representation"),
    cores: int = typer.Option(1, "--cores", help="Number of CPU cores for ORCA")
):
    if not xyz_file and not smiles_value:
        raise typer.BadParameter("Must provide either an XYZ file or a SMILES string.")

    base_name = "from_smiles" if smiles_value else xyz_file.stem
    info = None

    if smiles_value:
        try:
            canonical_smiles = sanitize_smiles(smiles_value)
            xyz_file = smiles_to_xyz(canonical_smiles, Path(f"{base_name}.xyz"))
            info = analyze_molecule(smiles_value=canonical_smiles)
            if generate_png:
                png_path = smiles_to_png(canonical_smiles, Path(f"{base_name}.png"))
                typer.echo(f"Generated PNG: {png_path}")
        except Exception as e:
            typer.echo(f"Analysis partially failed: {e}")
            info = {
                "SMILES": smiles_value,
                "StructureIntelligence": {
                    "Charge": 0,
                    "Multiplicity": 2 if "Cr" in smiles_value else 1,
                    "HasMetal": any(m in smiles_value for m in ["Cr","Fe","Ni","Cu","Zn"]),
                    "NumAtoms": 5,
                }
            }
    elif xyz_file and xyz_file.exists():
        try:
            info = analyze_molecule(xyz_path=xyz_file)
            if generate_png:
                png_path = xyz_to_png(xyz_file, Path(f"{base_name}.png"))
                typer.echo(f"Generated PNG: {png_path}")
        except Exception as e:
            typer.echo(f"Analysis partially failed: {e}")
            info = {
                "SMILES": None,
                "StructureIntelligence": {"Charge": 0, "Multiplicity": 1, "HasMetal": False, "NumAtoms": 5}
            }
    else:
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")

    while True:
        choice = Prompt.ask("\nSelect a section to view", choices=["descriptors", "pubchem", "ORCA", "exit"])
        if choice == "descriptors":
            if info:
                display_descriptors(info)
            else:
                typer.echo("No descriptor info available.")
        elif choice == "pubchem":
            if info:
                display_pubchem(info)
            else:
                typer.echo("No compound info available.")
        elif choice == "ORCA":
            if xyz_file.exists() and info:
                intel = info.get("StructureIntelligence", {}) if info else {}
                run_orca_workflow(xyz_file, base_name, intel=intel, ncores=cores)
            else:
                typer.echo("No valid XYZ structure available for ORCA.")
        elif choice == "exit":
            typer.echo("Exiting.")
            break


if __name__ == "__main__":
    app()
