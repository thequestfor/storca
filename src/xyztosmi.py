import openbabel

def xyz_to_smiles(xyz_file):
    """
    Convert XYZ file to SMILES using Open Babel.
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "smi")
    mol = openbabel.OBMol()
    if not obConversion.ReadFile(mol, xyz_file):
        raise ValueError("Failed to read XYZ file.")
    smiles = obConversion.WriteString(mol).strip()
    return smiles

