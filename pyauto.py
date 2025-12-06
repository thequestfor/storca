#!/usr/bin/env python3
from pathlib import Path
import typer

from src.inputgen import create_orca_input
from src.orca_runner import run_orca
from src.parser import parse_orca_orbitals
from src.converter import gbw_to_wfn
from src.wfngrid import run_multiwfn_fukui
from src.quickprops import analyze_xyz

app = typer.Typer(help="Automated ORCA + Multiwfn workflow with quick property estimates")

@app.command()
def run(
    xyz_file: Path = typer.Argument(..., help="Input XYZ file"),
    run_orca_flag: bool = typer.Option(True, help="Whether to run ORCA + Multiwfn workflow"),
    run_quickprops: bool = typer.Option(True, help="Whether to run quick property estimations"),
    run_fukui: bool = typer.Option(True, help="Whether to run Multiwfn Fukui calculation")
):
    # ----------------------
    # Check XYZ file
    # ----------------------
    if not xyz_file.exists():
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")

    base_name = xyz_file.stem

    # ----------------------
    # Quick property estimations
    # ----------------------
    if run_quickprops:
        print("\n=== Quick Molecular Property Estimation ===")
        from src.quickprops import analyze_xyz
    
        try:
            info = analyze_xyz(xyz_file)  # single dictionary returned
    
            # Separate descriptors vs PubChem properties (optional)
            descriptors = {k: info[k] for k in ["SMILES", "IUPACName", "MolecularWeight", "TPSA", "XLogP"]}
            properties = {k: info[k] for k in ["CID", "Synonyms", "Description"]}
    
            # Print descriptors
            print("\nDescriptors:")
            for k, v in descriptors.items():
                print(f"  {k}: {v}")
    
            # Print PubChem properties
            print("\nProperties from PubChem / Additional Info:")
            for k, v in properties.items():
                print(f"  {k}: {v if v else 'N/A'}")
    
        except Exception as e:
            print(f"Failed to estimate properties for {xyz_file}: {e}")
    


    # ----------------------
    # ORCA + Multiwfn workflow
    # ----------------------
    if run_orca_flag:
        print("\n=== Running ORCA + Multiwfn Workflow ===")

        # Step 1: Geometry Optimization
        opt_label = "opt"
        opt_inp = create_orca_input(
            xyz_file, charge=0, multiplicity=1, opt=True, label=opt_label
        )
        print(f"Created ORCA input for optimization: {opt_inp}")

        opt_outputs = run_orca(opt_inp)
        optimized_xyz = opt_outputs["xyz"]
        new_xyz_name = optimized_xyz.with_name(f"{base_name}_{opt_label}.xyz")
        optimized_xyz.rename(new_xyz_name)
        optimized_xyz = new_xyz_name

        print("\nGeometry optimization finished. Files generated:")
        for key, path in opt_outputs.items():
            if key == "xyz":
                path = optimized_xyz
            print(f"  {key}: {path}")

        # Step 2: Single-Point Calculations (neutral, +1, -1)
        charges = [0, +1, -1]
        sp_gbws = {}

        for charge in charges:
            mult = 1 if charge == 0 else 2
            label = f"sp{charge:+}"
            sp_inp = create_orca_input(
                optimized_xyz, charge=charge, multiplicity=mult, opt=False, label=label
            )
            print(f"\nRunning single-point calculation for charge {charge} (multiplicity {mult}): {sp_inp}")
            sp_outputs = run_orca(sp_inp)
            sp_gbws[charge] = sp_outputs["gbw"]
            print(f"  GBW generated: {sp_outputs['gbw']}")

        # Step 3: Parse HOMO/LUMO
        neutral_label = f"sp{0:+}"
        neutral_out = optimized_xyz.with_name(f"{base_name}_opt_{neutral_label}.out")
        orbitals = parse_orca_orbitals(neutral_out)
        print("\nParsed orbitals for neutral species:")
        print(f"  HOMO: {orbitals['homo_number']}  energy: {orbitals['homo_energy']} eV")
        print(f"  LUMO: {orbitals['lumo_number']}  energy: {orbitals['lumo_energy']} eV")

        # Step 4: Convert GBW → WFN
        wfn_files = {}
        for charge, gbw in sp_gbws.items():
            wfn = gbw_to_wfn(gbw)
            wfn_files[charge] = wfn
            print(f"GBW → WFN conversion done for charge {charge}: {wfn}")

        # Step 5: Rename WFNs for Multiwfn
        new_names = {
            0: Path("N.wfn"),
            1: Path("N+1.wfn"),
           -1: Path("N-1.wfn")
        }
        
        for charge, old_path in wfn_files.items():
            new_path = new_names[charge]
            old_path.rename(new_path)
            wfn_files[charge] = new_path
        
        print("\nRenamed WFN files:")
        for charge, path in wfn_files.items():
            print(f"  Charge {charge}: {path}")


if __name__ == "__main__":
    app()

