# quickprops/modules/formatter.py

def format_description(pubchem_props=None, smiles=None, molar_mass=None, logP=None, cid=None):
    """
    Create a Wikipedia-style macroscopic description.
    """
    description = {}
    if pubchem_props:
        description["Chemical formula"] = pubchem_props.get('MolecularFormula', 'N/A')
        description["IUPAC name"] = pubchem_props.get('IUPACName', 'N/A')
        description["SMILES"] = pubchem_props.get('CanonicalSMILES', smiles)
    else:
        description["SMILES"] = smiles
    
    description["Molar mass (g/mol)"] = molar_mass
    description["logP"] = logP
    if cid:
        description["PubChem CID"] = cid

    return description

