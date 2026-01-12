from pathlib import Path
import subprocess
import math
import sys
import re

# ==============================
# USER SETTINGS
# ==============================

def find_orca():
    """Return absolute path to ORCA executable or exit with error."""
    result = subprocess.run(
        ["which", "orca"],
        capture_output=True,
        text=True
    )
    if result.returncode != 0 or not result.stdout.strip():
        print("❌ ORCA executable not found in PATH")
        sys.exit(1)
    return result.stdout.strip()

ORCA_CMD = find_orca()

METHOD_LINE = "! B3LYP def2-SVP TightSCF Opt PAL4"
MAX_STRETCH = 1.8            # Å beyond equilibrium
STEPS = 8                    # number of scan points (NOT step size)
ENERGY_DROP_THRESHOLD = 0.5  # kcal/mol (downhill = unstable)

# ==============================
# BASIC PARSING
# ==============================

def read_xyz(xyz_path: Path):
    lines = xyz_path.read_text().splitlines()
    atoms = []
    for line in lines[2:]:
        if not line.strip():
            continue
        el, x, y, z = line.split()
        atoms.append([el, float(x), float(y), float(z)])
    return atoms

def distance(a, b):
    return math.sqrt(
        (a[1] - b[1])**2 +
        (a[2] - b[2])**2 +
        (a[3] - b[3])**2
    )

# ==============================
# BOND DETECTION
# ==============================

COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66
}

def find_bonds(atoms):
    bonds = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            ai, aj = atoms[i], atoms[j]

            # skip H–H
            if ai[0] == "H" and aj[0] == "H":
                continue

            ri = COVALENT_RADII.get(ai[0], 0.8)
            rj = COVALENT_RADII.get(aj[0], 0.8)

            if distance(ai, aj) < 1.25 * (ri + rj):
                bonds.append((i, j))

    return bonds

# ==============================
# ORCA INPUT GENERATION
# ==============================

def write_scan_input(atoms, bond, start, end, nsteps, label):
    i, j = bond
    inp = Path(f"scan_{label}.inp")

    lines = []
    lines.append(METHOD_LINE)
    lines.append("")
    lines.append("%geom")
    lines.append("  Scan")
    lines.append(f"    B {i+1} {j+1} = {start:.4f}, {end:.4f}, {nsteps}")
    lines.append("  end")
    lines.append("end")
    lines.append("")
    lines.append("* xyz 0 1")

    for sym, x, y, z in atoms:
        lines.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")

    lines.append("*")

    inp.write_text("\n".join(lines))
    return inp

# ==============================
# ORCA RUN + SCAN ENERGY PARSING
# ==============================

def run_orca_scan(inp: Path):
    """
    Run an ORCA scan and parse the energies from the output.

    Returns:
        energies (list of float): Scan energies in kcal/mol.
    """
    out = inp.with_suffix(".out")

    with open(out, "w") as f:
        subprocess.run(
            [ORCA_CMD, str(inp)],
            stdout=f,
            stderr=subprocess.STDOUT,
            check=False
        )

    text = out.read_text()
    energies = []
    in_table = False

    for line in text.splitlines():
        if "RELAXED SURFACE SCAN RESULTS" in line:
            in_table = True
            continue

        if in_table:
            if not line.strip():
                if energies:
                    break
                else:
                    continue

            m = re.match(r"\s*([0-9eE\.\+\-]+)\s+(-?[0-9eE\.\+\-]+)", line)
            if m:
                try:
                    e = float(m.group(2)) * 627.509  # Eh → kcal/mol
                    energies.append(e)
                except ValueError:
                    continue

    if len(energies) < 2:
        raise RuntimeError(f"Scan failed or produced insufficient points: {inp}")

    return energies

# ==============================
# SCAN LOGIC
# ==============================

def scan_bond(atoms, bond):
    i, j = bond
    r0 = distance(atoms[i], atoms[j])
    start = r0
    end = r0 + MAX_STRETCH

    label = f"bond_{i+1}_{j+1}"
    inp = write_scan_input(atoms, bond, start, end, STEPS, label)

    energies = run_orca_scan(inp)
    
    # Look for a local minimum
    min_index = energies.index(min(energies))
    
    # If the energy keeps dropping after the minimum, it's barrierless
    for e in energies[min_index+1:]:
        if e < energies[min_index] - ENERGY_DROP_THRESHOLD:
            return False  # barrierless dissociation
    
    return True  # metastable


# ==============================
# ENTRY POINT
# ==============================

def main(xyz_file):
    atoms = read_xyz(Path(xyz_file))
    bonds = find_bonds(atoms)

    print(f"Using ORCA at: {ORCA_CMD}")
    print(f"Found {len(bonds)} bonds to scan")

    for bond in bonds:
        print(f"Scanning bond {bond[0]+1}–{bond[1]+1}...")
        stable = scan_bond(atoms, bond)
        if not stable:
            print("❌ Molecule is NOT metastable (barrierless bond breaking detected)")
            sys.exit(1)

    print("✅ Molecule passed Phase 1 (no barrierless bond breaking detected)")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python stabtest.py molecule.xyz")
        sys.exit(1)

    main(sys.argv[1])
