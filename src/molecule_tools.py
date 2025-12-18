from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def sanitize_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid or ambiguous SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)

def smiles_to_xyz(smiles: str, output_path: Path) -> Path:
    smiles = sanitize_smiles(smiles)
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    conf = mol.GetConformer()
    with output_path.open("w") as f:
        f.write(f"{mol.GetNumAtoms()}\n\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
    return output_path

def smiles_to_png(smiles: str, output_path: Path) -> Path:
    smiles = sanitize_smiles(smiles)
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(output_path)
    return output_path

def xyz_to_png(xyz_path: Path, output_path: Path) -> Path:
    with xyz_path.open() as f:
        lines = f.readlines()[2:]
    mol = Chem.RWMol()
    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        mol.AddAtom(Chem.Atom(parts[0]))
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(output_path)
    return output_path
