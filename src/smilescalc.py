# quickprops/modules/smiles_calculations.py
from rdkit import Chem
from rdkit.Chem import Descriptors

def compute_rdkit_properties(smiles):
    """
    Compute molar mass and logP from SMILES.
    """
    mol = Chem.MolFromSmiles(smiles)
    molar_mass = Descriptors.MolWt(mol)
    logP = Descriptors.MolLogP(mol)
    return molar_mass, logP

