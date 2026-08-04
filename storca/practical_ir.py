"""Explicit, calculation-first rules for experiment-like IR display spectra.

The rules here are deliberately small and versioned.  They transform ORCA
normal modes; they do not predict a spectrum directly from SMILES.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class MolecularIRFeatures:
    smiles: str
    hydroxyl_donors: int
    amine_donors: int
    hydrogen_bond_acceptors: int
    rotatable_bonds: int


def describe_molecule(smiles: str) -> MolecularIRFeatures:
    """Extract only transparent descriptors used by the practical rules."""
    from rdkit import Chem
    from rdkit.Chem import Lipinski

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    hydroxyl = Chem.MolFromSmarts("[OX2H]")
    amine = Chem.MolFromSmarts("[NX3;H1,H2]")
    return MolecularIRFeatures(
        smiles=Chem.MolToSmiles(molecule, canonical=True),
        hydroxyl_donors=len(molecule.GetSubstructMatches(hydroxyl)),
        amine_donors=len(molecule.GetSubstructMatches(amine)),
        hydrogen_bond_acceptors=Lipinski.NumHAcceptors(molecule),
        rotatable_bonds=Lipinski.NumRotatableBonds(molecule),
    )


def practical_band_transform(
    frequency_cm_1: float, *, baseline_fwhm: float, features: MolecularIRFeatures
) -> tuple[float, float, list[str]]:
    """Apply the provisional ``ambient-organic-v0.3`` display rules.

    Rules apply only after the selected harmonic scale factor. They aim for a
    useful mixed-method organic reference display, not a particular solvent or
    sample preparation. Their values are intentionally recorded in metadata
    and will be revised only through held-out benchmark validation.
    """
    frequency = frequency_cm_1
    fwhm = baseline_fwhm
    rules = ["instrumental_baseline"]
    hbond_capable = features.hydrogen_bond_acceptors > 0

    if 800 <= frequency_cm_1 <= 1800 and features.rotatable_bonds >= 3:
        fwhm = max(fwhm, 25.0 + min(features.rotatable_bonds, 5) * 2.0)
        rules.append("flexible_fingerprint_width")
    return frequency, fwhm, rules


def practical_profile_metadata(features: MolecularIRFeatures) -> dict:
    return {
        "name": "ambient-organic-v0.3",
        "kind": "explicit_rule_calibration",
        "description": "Calculation-first broadening proxy for common organic experimental IR displays",
        "line_shape_model": "unit_area_gaussian_v2",
        "features": asdict(features),
        "limitations": "Does not model a specified solvent, crystal, concentration, or first-principles vibrational lifetime.",
    }
