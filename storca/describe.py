"""Offline, provenance-labelled molecule descriptions from RDKit."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FunctionalGroup:
    name: str
    smarts: str


FUNCTIONAL_GROUPS = (
    FunctionalGroup("alcohol_or_phenol", "[OX2H]"),
    FunctionalGroup("primary_or_secondary_amine", "[NX3;H1,H2]"),
    FunctionalGroup("carbonyl", "[CX3]=[OX1]"),
    FunctionalGroup("aldehyde", "[CH1X3](=O)"),
    FunctionalGroup("carboxylic_acid", "C(=O)[OX2H1]"),
    FunctionalGroup("ester", "C(=O)O[#6]"),
    FunctionalGroup("ether", "[OD2]([#6])[#6]"),
    FunctionalGroup("halogen", "[F,Cl,Br,I]"),
)


def _inchi_values(molecule) -> dict:
    from rdkit import Chem

    try:
        inchi = Chem.MolToInchi(molecule)
        return {"inchi": inchi, "inchi_key": Chem.InchiToInchiKey(inchi)} if inchi else {"inchi": None, "inchi_key": None}
    except Exception:
        return {"inchi": None, "inchi_key": None}


def _functional_groups(molecule) -> list[str]:
    from rdkit import Chem

    detected = []
    for group in FUNCTIONAL_GROUPS:
        pattern = Chem.MolFromSmarts(group.smarts)
        if pattern and molecule.HasSubstructMatch(pattern):
            detected.append(group.name)
    if any(atom.GetIsAromatic() for atom in molecule.GetAtoms()):
        detected.append("aromatic_ring")
    return detected


def _compound_class(molecule) -> str:
    from src.chem_intel import METALS

    symbols = {atom.GetSymbol() for atom in molecule.GetAtoms()}
    if symbols & METALS:
        return "organometallic" if "C" in symbols else "inorganic"
    return "organic" if "C" in symbols else "inorganic"


def _generated_description(result: dict) -> str:
    groups = result["functional_groups"]
    group_text = ", ".join(groups).replace("_", " ") if groups else "no named functional-group pattern"
    article = "an" if result["compound_class"][:1].lower() in "aeiou" else "a"
    return (
        f"This is {article} {result['compound_class']} molecule with formula {result['formula']} "
        f"and molar mass {result['molar_mass_g_mol']:.3f} g/mol. "
        f"Detected structural patterns: {group_text}."
    )


def describe_smiles(smiles: str) -> dict:
    """Return local, deterministic descriptors and their provenance."""
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    result = {
        "input_smiles": smiles,
        "canonical_smiles": canonical_smiles,
        "formula": formula,
        "compound_class": _compound_class(molecule),
        "formal_charge": Chem.GetFormalCharge(molecule),
        "radical_electrons": sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()),
        "molar_mass_g_mol": round(float(Descriptors.MolWt(molecule)), 6),
        "exact_mass_da": round(float(Descriptors.ExactMolWt(molecule)), 6),
        "heavy_atom_count": molecule.GetNumHeavyAtoms(),
        "atom_count_with_implicit_hydrogens": Chem.AddHs(molecule).GetNumAtoms(),
        "bond_count": molecule.GetNumBonds(),
        "ring_count": int(rdMolDescriptors.CalcNumRings(molecule)),
        "aromatic_ring_count": int(rdMolDescriptors.CalcNumAromaticRings(molecule)),
        "rotatable_bond_count": int(Lipinski.NumRotatableBonds(molecule)),
        "hydrogen_bond_donors": int(Lipinski.NumHDonors(molecule)),
        "hydrogen_bond_acceptors": int(Lipinski.NumHAcceptors(molecule)),
        "tpsa_a2": round(float(rdMolDescriptors.CalcTPSA(molecule)), 6),
        "xlogp_rdkit": round(float(Crippen.MolLogP(molecule)), 6),
        "fraction_csp3": round(float(rdMolDescriptors.CalcFractionCSP3(molecule)), 6),
        "functional_groups": _functional_groups(molecule),
        **_inchi_values(molecule),
        "provenance": {
            "source": "RDKit local calculation",
            "network_accessed": False,
            "notes": "Descriptors are structure-derived estimates or counts, not measured properties.",
        },
    }
    result["description"] = _generated_description(result)
    result["description_source"] = "generated from local RDKit descriptors"
    return result


def description_lines(result: dict) -> list[str]:
    """Format the most useful local fields for the command line."""
    fields = (
        ("Canonical SMILES", result["canonical_smiles"]),
        ("Formula", result["formula"]),
        ("Molar mass", f"{result['molar_mass_g_mol']:.3f} g/mol"),
        ("Exact mass", f"{result['exact_mass_da']:.6f} Da"),
        ("Class", result["compound_class"]),
        ("Formal charge", result["formal_charge"]),
        ("H-bond donors / acceptors", f"{result['hydrogen_bond_donors']} / {result['hydrogen_bond_acceptors']}"),
        ("TPSA", f"{result['tpsa_a2']:.2f} Å²"),
        ("XlogP (RDKit)", result["xlogp_rdkit"]),
        ("Rotatable bonds", result["rotatable_bond_count"]),
        ("Functional groups", ", ".join(result["functional_groups"]) or "none detected"),
    )
    return [f"{label}: {value}" for label, value in fields]
