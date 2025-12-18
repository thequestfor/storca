# src/molecule_tools.py
from pathlib import Path
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import Draw

def xyz_to_smiles(xyz_file: Path) -> str:
    """Convert XYZ file to SMILES using Open Babel."""
    conv = ob.OBConversion()
    conv.SetInFormat("xyz")
    conv.SetOutFormat("smi")
    mol = ob.OBMol()
    if not conv.ReadFile(mol, str(xyz_file)):
        raise ValueError(f"Failed to read XYZ file: {xyz_file}")
    smiles = conv.WriteString(mol).strip().split()[0]
    return smiles

def smiles_to_xyz(smiles: str, output_file: Path) -> Path:
    """Generate approximate XYZ from SMILES using Open Babel."""
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    conv.SetOutFormat("xyz")
    mol = ob.OBMol()
    conv.ReadString(mol, smiles)
    builder = ob.OBBuilder()
    builder.Build(mol)
    conv.WriteFile(mol, str(output_file))
    return output_file

def smiles_to_png(smiles: str, output_file: Path, size=(300, 300)) -> Path:
    """Generate 2D PNG from SMILES using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    Chem.rdDepictor.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=size)
    img.save(output_file)
    return output_file

def xyz_to_png(xyz_file: Path, output_file: Path, size=(300, 300)) -> Path:
    """Generate 2D PNG from XYZ by converting XYZ → SMILES first."""
    smiles = xyz_to_smiles(xyz_file)
    return smiles_to_png(smiles, output_file, size=size)

