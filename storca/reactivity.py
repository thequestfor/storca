"""ORCA-native, explicitly qualitative frontier-orbital reactivity summaries."""

from __future__ import annotations

import json
from pathlib import Path

from src.parser import parse_orca_orbitals


def frontier_reactivity_summary(orca_output: Path) -> dict:
    """Summarize frontier orbital energies without claiming atom-resolved reactivity.

    The derived quantities use a Koopmans-like frontier-orbital approximation.
    They are useful for comparing consistently calculated molecules, not as
    measured ionization energies, electron affinities, redox potentials, or
    site-selectivity predictions.
    """
    orca_output = Path(orca_output)
    orbitals = parse_orca_orbitals(orca_output)
    homo, lumo = orbitals["homo_energy"], orbitals["lumo_energy"]
    if homo is None or lumo is None:
        raise ValueError("ORCA output must contain both occupied and virtual frontier orbitals")
    gap = float(lumo - homo)
    if gap <= 0:
        raise ValueError("Frontier orbital gap is not positive; cannot form a qualitative closed-shell summary")
    chemical_potential = -0.5 * (homo + lumo)
    hardness = 0.5 * gap
    softness = 1.0 / (2.0 * hardness)
    electrophilicity = (chemical_potential**2) / (2.0 * hardness)
    return {
        "schema_version": 1,
        "kind": "qualitative_frontier_orbital_summary",
        "orca_output": str(orca_output),
        "frontier_orbitals": {
            "homo_number": orbitals["homo_number"],
            "homo_energy_eV": homo,
            "lumo_number": orbitals["lumo_number"],
            "lumo_energy_eV": lumo,
            "gap_eV": gap,
        },
        "derived_frontier_proxies": {
            "chemical_potential_proxy_eV": chemical_potential,
            "hardness_proxy_eV": hardness,
            "softness_proxy_eV_inverse": softness,
            "electrophilicity_proxy_eV": electrophilicity,
        },
        "provenance": {
            "source": "Parsed directly from an ORCA orbital-energy table",
            "model": "Koopmans-like frontier-orbital approximation",
            "network_accessed": False,
        },
        "limitations": [
            "Kohn-Sham orbital energies are not measured ionization energies, electron affinities, or redox potentials.",
            "This summary is molecule-level only; it does not identify reactive atoms or provide Fukui functions.",
            "Compare values only across calculations with compatible charge, spin, geometry, method, and basis settings.",
        ],
    }


def write_frontier_reactivity(path: Path, report: dict) -> Path:
    path = Path(path)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return path
