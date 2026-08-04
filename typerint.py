#!/usr/bin/env python3
from pathlib import Path
import typer
from rich.prompt import Prompt
from typing import Optional
import subprocess
import time
import threading
import numpy as np
from src.inputgen import (
    create_orca_input,
    create_rmg_input
)
from src.orca_runner import (
    run_orca,
    run_rmg
)
from src.parser import (
    parse_orca_orbitals,
    parse_orca_energy,
    parse_orca_ir,
    parse_goat_out,
    parse_xyz_ensemble,
    attach_xyz_to_conformers,
    parse_chemkin_annotated
)
from src.converter import gbw_to_wfn
from src.wfngrid import run_multiwfn_fukui
from src.quickprops import analyze_molecule
from src.molecule_tools import (
    write_weighted_orca_freq,
    smiles_to_xyz,
    smiles_to_png,
    generate_unique_conformers,
    smiles_to_conformers,
    xyz_to_png,
    sanitize_smiles,
)
from src.plotter import build_ir_spectrum
from src.stability.freq_check import frequency_stability_check
from src.utils.orca_xyz_clean import clean_xyz_for_orca
import csv 
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


# -------------------------------
# Frequency compression function
# -------------------------------
def compress_wavenumber(freq_raw):
    """
    Smoothly compress harmonic frequencies to experimental scale.
    
    Parameters
    ----------
    freq_raw : float or ndarray
        Harmonic frequency in cm^-1

    Returns
    -------
    freq_scaled : float or ndarray
        Compressed frequency
    """
    f0 = 2000.0     # midpoint of compression (cm^-1)
    delta = 1200.0  # controls how sharply compression ramps
    min_scale = 0.9  # scale for very high freq
    max_scale = 0.995  # scale for low freq
    
    scale = min_scale + (max_scale - min_scale) / (1 + np.exp((freq_raw - f0)/delta))
    
    return freq_raw * scale

def run_stability_analysis(
    smiles: str,
    rmg_env: str = "rmg_env",
    rmg_command: str = "rmg.py",
):
    typer.echo("\n=== Running RMG Stability Analysis ===")

    # Main working directory
    workdir = Path("rmg_stability")
    workdir.mkdir(parents=True, exist_ok=True)

    # Generate RMG input
    species_label = "stability"
    input_py = create_rmg_input(species_label, smiles, workdir)

    typer.echo(f"RMG input generated: {input_py}")
    typer.echo("Launching RMG... (this may take time)")

    # Run RMG
    try:
        run_rmg(input_file=input_py, rmg_env=rmg_env)
    except subprocess.CalledProcessError as e:
        typer.echo(f"⚠ RMG failed:\n{e.stdout}\n{e.stderr}")
        return

    typer.echo("RMG finished. Parsing results...")

    # Parse annotated Chemkin output
    chemkin_file = workdir / "chemkin" / "chem_annotated.inp"
    result = parse_chemkin_annotated(chemkin_file)

    typer.echo("\n=== Stability Assessment ===")
    typer.echo(f"Stable: {result['stable']}")

    if not result["stable"]:
        typer.echo("\nDecomposition reactions detected:")
        for rxn in result["decomposition_reactions"]:
            typer.echo(f"  {rxn}")

        if result["reaction_barriers"]:  # updated key
            typer.echo("\nReaction barriers (kcal/mol):")
            for barrier in result["reaction_barriers"]:
                typer.echo(f"  {barrier:.3f}")
def run_conformers(
    xyz_initial: str = typer.Option(..., help="filename of XYZ geometry of molecule"),
    out_dir: Path = typer.Option(Path("conformers"), help="Output directory"),
    ncores: int = typer.Option(1, help="Number of CPU cores for ORCA"),
    temperature: float = typer.Option(298.15, help="Temperature for Boltzmann weighting in K")
):
    """
    Run GOAT conformer search, optimize significant conformers, and compute Boltzmann weights.
    """

    out_dir.mkdir(parents=True, exist_ok=True)
    if hasattr(temperature, "default"):
        temperature = float(temperature.default)

    typer.echo("\nGenerating input file for conformer search (GOAT)...")
    inp_file = create_orca_input(
        xyz_file=xyz_initial,
        charge=0,
        multiplicity=1,
        use_goat=True,
        label=out_dir / "goat",
        ncores=ncores
    )
    typer.echo(f"   Input file generated at {inp_file}")
    typer.echo("\nLaunching conformer search...\n")

    # -------------------------------
    # Helper: run ORCA/GOAT with monitoring
    # -------------------------------
    def run_goat_with_monitor(inp_file: Path, out_dir: Path, stall_threshold: int = 60):
        goat_out = inp_file.with_suffix(".out")
        out_dir.mkdir(parents=True, exist_ok=True)

        def run_orca_thread():
            try:
                run_orca(inp_file, capture_out=True, stdout_file=goat_out)
            except Exception as e:
                run_orca_thread.error = e

        run_orca_thread.error = None
        t = threading.Thread(target=run_orca_thread)
        t.start()

        stalled_count = 0
        last_sizes = {f: f.stat().st_size for f in out_dir.glob("*") if f.is_file()}

        while t.is_alive():
            total_size = 0
            changed = False

            for f in out_dir.glob("*"):
                if not f.is_file():
                    continue
                try:
                    size = f.stat().st_size
                except FileNotFoundError:
                    continue

                total_size += size
                if f not in last_sizes or size != last_sizes[f]:
                    last_sizes[f] = size
                    changed = True

            if changed:
                stalled_count = 0
            else:
                stalled_count += 1

            status = f"GOAT running... thread alive, total output size: {total_size} bytes"
            if stalled_count > 0:
                status += f" | {stalled_count}s since last change"
            typer.echo(f"\r{status}", nl=False)

            if stalled_count >= stall_threshold:
                typer.echo(f"\n⚠ No file changes for {stall_threshold}s, computation may be slow or frozen.")

            time.sleep(1)

        t.join()
        typer.echo("\nGOAT run finished.")

        if run_orca_thread.error is not None:
            typer.echo(f"⚠ ORCA crashed with error: {run_orca_thread.error}")
            typer.echo("Check goat.out for details. Consider lowering the number of cores.")
            return False

        if not goat_out.exists() or goat_out.stat().st_size == 0:
            typer.echo("⚠ ORCA finished but goat.out is missing or empty. Check input and environment.")
            return False

        typer.echo(f"Output directory: {out_dir}")
        typer.echo(f"Main output file: {goat_out}")

        xyz_file = inp_file.with_name("goat.finalensemble.xyz")
        if xyz_file.exists():
            typer.echo(f"XYZ file: {xyz_file}")
        else:
            typer.echo("Warning: XYZ file not found!")

        return True

    # -------------------------------
    # Run GOAT
    # -------------------------------
    success = run_goat_with_monitor(inp_file, out_dir, stall_threshold=60)
    if not success:
        return

    # -------------------------------
    # Parse GOAT output
    # -------------------------------
    goat_out = inp_file.with_suffix(".out")
    conformers, sig_conformers = parse_goat_out(goat_out)
    typer.echo(f"\nTotal conformers found: {len(conformers)}")
    typer.echo(f"Significant conformers (>0.05 weight): {len(sig_conformers)}")

    # -------------------------------
    # Parse XYZ ensemble
    # -------------------------------
    xyz_file = inp_file.with_name("goat.finalensemble.xyz")
    typer.echo(f"Parsing XYZ ensemble: {xyz_file}")
    xyz_frames = parse_xyz_ensemble(xyz_file)
    typer.echo(f"XYZ frames found: {len(xyz_frames)}")

    # -------------------------------
    # Attach XYZ frames to conformers
    # -------------------------------
    conformers = attach_xyz_to_conformers(conformers, xyz_frames)
    typer.echo("Conformers successfully matched with XYZ frames.")

    # -------------------------------
    # Optimize significant conformers using run_orca
    # -------------------------------
    kB_Hartree_per_K = 3.1668114e-6
    opt_dir = out_dir / "optimization"
    opt_dir.mkdir(exist_ok=True)

    for c in sig_conformers:
        idx = c["idx"]
        typer.echo(f"\nOptimizing conformer {idx}...")

        xyz_path = opt_dir / f"conformer_{idx}.xyz"
        xyz_lines = c["xyz"]["xyz_lines"]  # Get the list of atom lines.
        xyz_content = "\n".join(xyz_lines) + "\n"  # Join lines into a single string
        with open(xyz_path, "w") as f:
            f.write(xyz_content)


        orca_inp = create_orca_input(
            xyz_file=xyz_path,
            charge=0,
            multiplicity=1,
            label=opt_dir / f"conformer_{idx}",
            ncores=ncores
        )

        orca_out = run_orca(orca_inp, capture_out=True)["out"]

        # Parse energy
        energy = None
        with open(orca_out) as f:
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:
                    energy = float(line.split()[-1])
        if energy is None:
            typer.echo(f"⚠ Could not find energy for conformer {idx}, skipping...")
            continue

        c["energy"] = energy
        typer.echo(f"Conformer {idx} optimized: Energy = {energy:.6f} Hartree")
    # -------------------------------
    # Boltzmann weights (after ALL optimizations)
    # -------------------------------
    kB_Hartree_per_K = 3.1668114e-6
    
    energies = np.array([c["energy"] for c in sig_conformers])
    E_min = energies.min()
    delta_E = energies - E_min
    
    boltz_factors = np.exp(-delta_E / (kB_Hartree_per_K * temperature))
    weights = boltz_factors / boltz_factors.sum()
    
    for c, w in zip(sig_conformers, weights):
        c["weight"] = w
        typer.echo(
            f"Conformer {c['idx']}: Energy = {c['energy']:.6f} Hartree, Weight = {w:.4f}"
        )

    # -------------------------------
    # Save CSV summary
    # -------------------------------
    summary_file = opt_dir / "conformer_summary.csv"
    with open(summary_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["idx", "energy", "weight"])
        writer.writeheader()
        for c in sig_conformers:
            writer.writerow({
                "idx": c["idx"],
                "energy": c["energy"],
                "weight": c["weight"]
            })
    
    typer.echo(f"Summary saved to {summary_file}")
    
    # -------------------------------
    # Frequency calculations
    # -------------------------------
    freq_dir = opt_dir / "frequencies"
    freq_dir.mkdir(exist_ok=True)
    
    for c in sig_conformers:
        idx = c["idx"]
        typer.echo(f"\nRunning frequency calculation for conformer {idx}...")
    
        # ---- Write XYZ ----
        xyz_path = freq_dir / f"conformer_{idx}.xyz"
        xyz_lines = c["xyz"]["xyz_lines"]
        xyz_content = "\n".join(xyz_lines) + "\n"
    
        with open(xyz_path, "w") as f:
            f.write(xyz_content)
    
        # ---- ORCA frequency input ----
        freq_inp = create_orca_input(
            xyz_file=xyz_path,
            charge=0,
            multiplicity=1,
            label=freq_dir / f"conformer_{idx}_freq",
            ncores=ncores,
            freq=True
        )
    
        freq_out = freq_inp.with_suffix(".out")
    
        # ---- Run ORCA ----
        try:
            run_orca(freq_inp, capture_out=True, stdout_file=freq_out)
        except Exception as e:
            typer.echo(f"⚠ Frequency calculation failed for conformer {idx}: {e}")
            c["freq_failed"] = True
            continue
    
        if not freq_out.exists() or freq_out.stat().st_size == 0:
            typer.echo(f"⚠ Frequency output missing for conformer {idx}")
            c["freq_failed"] = True
            continue
    
        typer.echo(f"   Frequency calculation finished for conformer {idx}")
    
        # Mark success (parsing comes later)
        c["freq_out"] = freq_out
        c["freq_failed"] = False
    # -------------------------------
    # Parse IR frequencies
    # -------------------------------
    for c in sig_conformers:
        if c.get("freq_failed", False):
            continue
    
        freq_out = c.get("freq_out")
        if freq_out is None:
            continue
    
        typer.echo(f"Parsing IR data for conformer {c['idx']}...")
    
        vibrations = parse_orca_ir(freq_out)  # returns list of dicts
    
        # Extract separate lists
        freqs = [v["freq"] for v in vibrations]
        intensities = [v["intensity"] for v in vibrations]
    
        c["ir_freqs"] = freqs
        c["ir_intensities"] = intensities
    
    # -------------------------------
    # Build Boltzmann-weighted IR spectrum in transmittance
    # -------------------------------
    
    typer.echo("\nBuilding Boltzmann-weighted IR spectrum in transmittance...")
    
    # -------------------------------
    # Frequency grid
    # -------------------------------
    freq_min, freq_max = 400, 3800  # cm^-1
    resolution = 1.0
    fwhm = 25.0  # cm^-1
    freq_grid = np.arange(freq_min, freq_max + resolution, resolution)
    
    # Accumulate molar absorptivity (L mol^-1 cm^-1)
    spectrum_eps = np.zeros_like(freq_grid, dtype=float)
    
    # -------------------------------
    # Renormalize Boltzmann weights
    # -------------------------------
    W = sum(
        c["weight"] for c in sig_conformers
        if not c.get("freq_failed", False)
    )
    for c in sig_conformers:
        if not c.get("freq_failed", False):
            c["weight"] /= W
    
    # -------------------------------
    # Loop over conformers
    # -------------------------------
    sigma = fwhm / 2.35482  # Gaussian sigma (cm^-1)
    norm = sigma * np.sqrt(2.0 * np.pi)
    
    for c in sig_conformers:
        if c.get("freq_failed", False):
            continue
    
        vibrations = parse_orca_ir(c["freq_out"])  # dicts: freq, eps, intensity
        weight = c["weight"]
    
        for v in vibrations:
            # Compressed fundamental frequency
            freq = compress_wavenumber(v["freq"])
    
            # Boltzmann-weighted molar absorptivity
            eps = v["eps"] * weight
    
            # Normalized Gaussian broadening
            gauss = np.exp(-(freq_grid - freq) ** 2 / (2.0 * sigma ** 2)) / norm
            spectrum_eps += eps * gauss
    
            # ---- Optional first overtone for strong C–H stretches ----
            if v["freq"] > 2800.0 and v["eps"] > 5.0:
                overtone_freq = compress_wavenumber(2.0 * v["freq"])
                overtone_eps = 0.03 * eps
                gauss_ot = np.exp(-(freq_grid - overtone_freq) ** 2 / (2.0 * sigma ** 2)) / norm
                spectrum_eps += overtone_eps * gauss_ot
    
    # -------------------------------
    # Experimental absorbance mapping
    # -------------------------------
    
    # Global absorbance scaling (effective path length × concentration)
    # Neat liquids typically require much larger values than solutions
    ABS_SCALE = 2000.0   # reasonable starting point for neat liquid IR
    
    absorbance_raw = ABS_SCALE * spectrum_eps
    
    # Detector saturation / finite dynamic range
    A_max = 3.0  # ~0.1% transmittance floor
    absorbance = A_max * (1.0 - np.exp(-absorbance_raw / A_max))
    
    # Convert to transmittance
    transmittance = 10.0 ** (-absorbance)
    
    # -------------------------------
    # Save CSV
    # -------------------------------
    csv_path = out_dir / "ir_spectrum_transmittance.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Wavenumber_cm-1", "Transmittance"])
        for x, y in zip(freq_grid, transmittance):
            writer.writerow([x, y])
    
    typer.echo(f"IR transmittance spectrum saved to {csv_path}")
    
    # -------------------------------
    # Plot
    # -------------------------------
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(8, 4))
    plt.plot(freq_grid, transmittance, color="black", linewidth=1.2)
    plt.gca().invert_xaxis()
    plt.ylim(0.0, 1.02)
    plt.xlabel("Wavenumber (cm$^{-1}$)")
    plt.ylabel("Transmittance")
    plt.title("Boltzmann-weighted IR Spectrum (Neat Liquid, Compressed)")
    plt.tight_layout()
    
    png_path = out_dir / "ir_spectrum_transmittance.png"
    plt.savefig(png_path, dpi=300)
    plt.close()
    
    typer.echo(f"IR spectrum plot saved to {png_path}")


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
        choice = Prompt.ask("\nSelect a section to view", choices=["descriptors", "pubchem", "ORCA", "spec", "stability", "exit"])
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
        elif choice == "spec":
            if smiles_value:
                run_conformers(
                    xyz_initial=xyz_file,
                    out_dir=Path("conformers"),
                    ncores=cores
                )
            else:
                typer.echo("Conformer generation currently only supported from SMILES input.")
        elif choice == "stability":
            if smiles_value:
                typer.echo("Running stability analysis...")

                # Provide all required arguments
                run_stability_analysis(
                    smiles=smiles_value,
                    rmg_env="rmg_env",  # default, can also make this configurable
                )
        elif choice == "exit":
            typer.echo("Exiting.")
            break


if __name__ == "__main__":
    typer.echo(
        "DEPRECATED: typerint.py is archived prototype code. Use `storca --help`; "
        "the supported STORCA 2.0 path does not use Multiwfn.",
        err=True,
    )
    app()
