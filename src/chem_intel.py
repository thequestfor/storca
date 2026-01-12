# src/chem_intel.py
from rdkit import Chem

METALS = {
    "Li","Be","Na","Mg","Al","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co",
    "Ni","Cu","Zn","Ga","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
    "Ag","Cd","In","Sn","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd",
    "Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
    "Au","Hg","Tl","Pb","Bi","Po","Fr","Ra","Ac","Th","Pa","U"
}

def analyze_structure(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "Valid": False,
            "CompoundClass": "unknown",
            "HasMetal": False,
            "Atoms": [],
        }

    # --- Add explicit hydrogens for correct atom counting ---
    mol_h = Chem.AddHs(mol)

    atoms = [a.GetSymbol() for a in mol_h.GetAtoms()]
    atom_set = set(atoms)

    has_metal = bool(atom_set & METALS)
    has_carbon = "C" in atom_set

    # --- Correct chemical classification ---
    if has_metal and has_carbon:
        compound_class = "organometallic"
    elif has_carbon:
        compound_class = "organic"
    else:
        compound_class = "inorganic"

    # --- Radical check (very coarse; marked as estimated) ---
    num_radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    multiplicity = 2 if num_radicals else 1

    return {
        "Valid": True,
        "CompoundClass": compound_class,
        "HasMetal": has_metal,
        "Atoms": sorted(atom_set),
        "Charge": Chem.GetFormalCharge(mol),
        "Multiplicity": multiplicity,
        "MultiplicityEstimated": True,
        "NumAtoms": mol_h.GetNumAtoms(),          # total atoms (incl. H)
        "NumHeavyAtoms": mol.GetNumHeavyAtoms(),  # heavy atoms only
        "Aromatic": any(a.GetIsAromatic() for a in mol.GetAtoms()),
    }
