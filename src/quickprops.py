# src/quickprops.py

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import subprocess
import tempfile
import os

def xyz_to_smiles(xyz_file):
    """
    Converts an XYZ file to SMILES using Open Babel.
    RDKit cannot read XYZ directly, so we rely on obabel CLI.
    """
    with tempfile.NamedTemporaryFile(suffix=".mol", delete=False) as tmpmol:
        tmpmol_name = tmpmol.name

    # Use Open Babel to convert XYZ -> MOL
    cmd = ["obabel", xyz_file, "-O", tmpmol_name]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Read MOL with RDKit
    mol = rdmolfiles.MolFromMolFile(tmpmol_name)
    os.remove(tmpmol_name)
    if mol is None:
        raise ValueError(f"Cannot convert {xyz_file} to a valid molecule.")
    
    smiles = Chem.MolToSmiles(mol)
    return smiles

def fetch_pubchem_info(identifier):
    """
    Query PubChem using SMILES, InChI, or CID to get full substance info.
    Returns a dictionary of descriptors and macroscopic properties.
    """
    compounds = pcp.get_compounds(identifier, 'smiles')
    if not compounds:
        return None, None
    
    compound = compounds[0]

    # Basic molecular descriptors
    descriptors = {
        "SMILES": compound.isomeric_smiles,
        "CID": compound.cid,
        "MolecularWeight": compound.molecular_weight,
        "IUPACName": compound.iupac_name,
        "XLogP": compound.xlogp,
        "TPSA": compound.tpsa
    }

    # Fetch full substance records for macroscopic info
    substances = pcp.get_substances(compound.cid)
    macroscopic = {}
    if substances:
        sub = substances[0]
        macroscopic = {
            "CID": compound.cid,
            "MolecularWeight": compound.molecular_weight,
            "IUPACName": compound.iupac_name,
            "XLogP": compound.xlogp,
            "TPSA": compound.tpsa,
            "PhysicalDescription": getattr(sub, "physical_description", None),
            "RecordTitle": getattr(sub, "title", None),
            "Synonyms": getattr(sub, "synonyms", [])
        }

    return descriptors, macroscopic

def analyze_xyz(xyz_file):
    """
    Main entry point: convert XYZ -> SMILES -> query PubChem.
    Returns descriptors and macroscopic properties.
    """
    try:
        smiles = xyz_to_smiles(xyz_file)
    except Exception as e:
        print(f"Error converting XYZ to SMILES: {e}")
        return None, None
    
    try:
        descriptors, macroscopic = fetch_pubchem_info(smiles)
        return descriptors, macroscopic
    except Exception as e:
        print(f"Error fetching PubChem data: {e}")
        return None, None

