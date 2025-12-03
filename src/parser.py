from pathlib import Path

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

    homo_orb = max([o for o in orbitals if o["occ"] > 0], key=lambda x: x["no"])
    lumo_orb = min([o for o in orbitals if o["occ"] == 0], key=lambda x: x["no"])

    return {
        "homo_number": homo_orb["no"],
        "homo_energy": homo_orb["ev"],
        "lumo_number": lumo_orb["no"],
        "lumo_energy": lumo_orb["ev"]
    }

