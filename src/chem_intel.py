# src/chem_intel.py
from rdkit import Chem

METALS = set([
    "Li","Be","Na","Mg","Al","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co",
    "Ni","Cu","Zn","Ga","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
    "Ag","Cd","In","Sn","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd",
    "Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
    "Au","Hg","Tl","Pb","Bi","Po","Fr","Ra","Ac","Th","Pa","U"
])

def analyze_structure(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "Valid": False,
            "CompoundClass": "unknown",
            "HasMetal": False,
            "Atoms": [],
        }

    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    atom_set = set(atoms)
    has_metal = bool(atom_set & METALS)

    # Radical check for multiplicity
    num_radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    multiplicity = 2 if num_radicals else 1

    return {
        "Valid": True,
        "CompoundClass": "inorganic/metal" if has_metal else "organic",
        "HasMetal": has_metal,
        "Atoms": sorted(atom_set),
        "Charge": Chem.GetFormalCharge(mol),
        "Multiplicity": multiplicity,
        "NumAtoms": mol.GetNumAtoms(),
        "NumHeavyAtoms": mol.GetNumHeavyAtoms(),
        "Aromatic": any(a.GetIsAromatic() for a in mol.GetAtoms()),
    }


