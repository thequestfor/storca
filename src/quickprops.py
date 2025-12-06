import json
import requests
import pubchempy as pcp
from openbabel import openbabel as ob
import os

# ===========================
# XYZ → SMILES CONVERSION
# ===========================
def xyz_to_smiles(xyz_file):
    """Convert an XYZ file to SMILES using Open Babel."""
    xyz_file = str(xyz_file)  # ensure string path
    if not os.path.exists(xyz_file):
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")

    conv = ob.OBConversion()
    conv.SetInFormat("xyz")
    conv.SetOutFormat("smi")

    mol = ob.OBMol()
    success = conv.ReadFile(mol, xyz_file)
    if not success:
        raise ValueError(f"Failed to read XYZ file: {xyz_file}")

    smiles = conv.WriteString(mol).strip()
    if not smiles:
        raise ValueError(f"OpenBabel failed to convert {xyz_file} to SMILES")

    return smiles
# ===========================
# NIH CACTUS API: SMILES → IUPAC
# ===========================
def smiles_to_iupac(smiles):
    """Convert a SMILES string to IUPAC using NIH Cactus API."""
    url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/iupac_name"
    try:
        r = requests.get(url, timeout=5)
        if r.ok and r.text.strip():
            return r.text.strip()
    except requests.RequestException:
        pass
    return None

# ===========================
# PUG-VIEW DESCRIPTION PARSER
# ===========================
def extract_relevant_descriptions(data, keywords=None, max_paragraphs=3):
    """
    Extract descriptive text from PubChem PUG-View JSON, filter by length and keyword relevance.
    """
    if keywords is None:
        keywords = [
            "boiling", "melting", "solubility", "density", "logP", "XLogP",
            "odor", "taste", "toxicity", "hazard", "flammable", "reactivity"
        ]

    descriptions = []

    def walk(node):
        if isinstance(node, dict):
            if "StringWithMarkup" in node:
                for item in node["StringWithMarkup"]:
                    if "String" in item:
                        descriptions.append(item["String"])
            if "String" in node and isinstance(node["String"], str):
                descriptions.append(node["String"])
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(data)

    # Filter by length
    long_descriptions = [d for d in descriptions if len(d.split()) >= 8]

    # Rank by relevance
    def relevance_score(text):
        text_lower = text.lower()
        return sum(1 for kw in keywords if kw.lower() in text_lower)

    ranked = sorted(
        long_descriptions,
        key=lambda d: (relevance_score(d), len(d)),
        reverse=True
    )

    # Return only top N
    return "\n\n".join(ranked[:max_paragraphs]) if ranked else "No relevant description available."

# ===========================
# FETCH PUG-VIEW DESCRIPTION
# ===========================
def fetch_full_description(cid, max_paragraphs=3):
    """Download and extract the long-form narrative description from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    try:
        r = requests.get(url, timeout=5)
        if not r.ok:
            return None
        data = r.json()
        return extract_relevant_descriptions(data, max_paragraphs=max_paragraphs)
    except requests.RequestException:
        return None

# ===========================
# MASTER FUNCTION
# ===========================
def analyze_xyz(xyz_path):
    """
    Returns:
      {
        "SMILES": ...,
        "CID": ...,
        "IUPACName": ...,
        "MolecularWeight": ...,
        "TPSA": ...,
        "XLogP": ...,
        "Synonyms": [...],
        "Description": "long descriptive text from PubChem"
      }
    """

    # --- 1. XYZ → SMILES ---
    smiles = xyz_to_smiles(xyz_path)

    # --- 2. Attempt PubChem lookup ---
    compounds = pcp.get_compounds(smiles, namespace="smiles")
    if compounds:
        compound = compounds[0]
        cid = compound.cid
        mw = compound.molecular_weight
        xlogp = compound.xlogp
        tpsa = compound.tpsa
        synonyms = compound.synonyms or []
        # Prefer PubChem IUPAC if available, else fallback to Cactus
        name = compound.iupac_name or smiles_to_iupac(smiles)
        description = fetch_full_description(cid) or "No description available."
    else:
        # Molecule not found in PubChem
        cid = None
        mw = None
        xlogp = None
        tpsa = None
        synonyms = []
        name = smiles_to_iupac(smiles) or "Unknown"
        description = "No description available."

    # --- 3. Return structured block ---
    return {
        "SMILES": smiles,
        "CID": cid,
        "IUPACName": name,
        "MolecularWeight": mw,
        "TPSA": tpsa,
        "XLogP": xlogp,
        "Synonyms": synonyms,
        "Description": description
    }

# ===========================
# DEBUG/STANDALONE TEST
# ===========================
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python quick_props.py <file.xyz>")
        sys.exit(1)
    xyz = sys.argv[1]
    info = analyze_xyz(xyz)
    print(json.dumps(info, indent=2))

