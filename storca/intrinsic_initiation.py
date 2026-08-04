"""Generic radical-fragment discovery without molecule-class alerts."""

from __future__ import annotations


def enumerate_homolytic_cleavages(smiles: str) -> dict:
    """Enumerate every covalent single-bond homolysis as radical fragments.

    These are route-discovery candidates only.  Bond dissociation energy and
    pressure-dependent kinetics are deliberately not inferred from topology.
    """
    try:
        from rdkit import Chem
    except ImportError as error:
        raise ImportError("Intrinsic initiation discovery requires RDKit") from error
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES")
    candidates = []
    for bond in molecule.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        editable = Chem.RWMol(molecule)
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        editable.RemoveBond(left, right)
        editable.GetAtomWithIdx(left).SetNumRadicalElectrons(
            editable.GetAtomWithIdx(left).GetNumRadicalElectrons() + 1)
        editable.GetAtomWithIdx(right).SetNumRadicalElectrons(
            editable.GetAtomWithIdx(right).GetNumRadicalElectrons() + 1)
        fragments = Chem.GetMolFrags(editable.GetMol(), asMols=True, sanitizeFrags=True)
        candidates.append({
            "bond_atom_indices": [left, right],
            "bond_type": "single",
            "fragment_radical_smiles": sorted(Chem.MolToSmiles(fragment, canonical=True) for fragment in fragments),
            "kinetic_status": "requires_orca_and_pressure_dependent_rate",
        })
    return {
        "kind": "generic_intrinsic_radical_initiation_discovery",
        "canonical_smiles": Chem.MolToSmiles(molecule, canonical=True),
        "candidates": candidates,
        "interpretation": "Every covalent single-bond homolysis is enumerated uniformly; no structural-alert classes are used.",
    }
