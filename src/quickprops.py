# src/quickprops.py
import os
import json
import requests
import pubchempy as pcp
from openbabel import openbabel as ob
from rdkit import Chem
from . import hazards  # access via hazards.classify_hazards / hazards.estimate_practical_hazards
from .chem_intel import analyze_structure

GOOD_SECTIONS = {"description", "physical description", "chemical and physical properties", "safety and hazards", "uses"}

# ------------------------------
# XYZ → SMILES using OpenBabel
# ------------------------------
def xyz_to_smiles(xyz_path):
    if not os.path.exists(str(xyz_path)):
        raise FileNotFoundError(f"XYZ file not found: {xyz_path}")

    conv = ob.OBConversion()
    conv.SetInFormat("xyz")
    conv.SetOutFormat("smi")
    mol = ob.OBMol()
    if not conv.ReadFile(mol, str(xyz_path)):
        raise ValueError(f"Failed to read XYZ file: {xyz_path}")

    smiles = conv.WriteString(mol).strip().split()[0]
    if not smiles:
        raise ValueError("OpenBabel failed to generate SMILES")

    return smiles

# ------------------------------
# Canonical SMILES using RDKit
# ------------------------------
def canonicalize_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)

# ------------------------------
# IUPAC Name via NIH CACTUS
# ------------------------------
def smiles_to_iupac(smiles: str):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/iupac_name"
    try:
        r = requests.get(url, timeout=5)
        if r.ok and r.text.strip():
            return r.text.strip()
    except requests.RequestException:
        pass
    return None

# ------------------------------
# PubChem description helpers
# ------------------------------
def extract_pubchem_description(data, max_paragraphs=3):
    collected = []
    def walk(node, current_heading=None):
        if isinstance(node, dict):
            heading = node.get("TOCHeading", current_heading)
            heading_norm = heading.lower() if isinstance(heading, str) else None
            if heading_norm in GOOD_SECTIONS and "StringWithMarkup" in node:
                for item in node["StringWithMarkup"]:
                    text = item.get("String")
                    if text and len(text.split()) >= 8:
                        collected.append(text)
            for v in node.values():
                walk(v, heading)
        elif isinstance(node, list):
            for item in node:
                walk(item, current_heading)
    walk(data)
    return "\n\n".join(collected[:max_paragraphs]) if collected else None

def fetch_pubchem_description(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    try:
        r = requests.get(url, timeout=5)
        if r.ok:
            return extract_pubchem_description(r.json())
    except requests.RequestException:
        pass
    return None

# ------------------------------
# ChEBI description
# ------------------------------
def fetch_chebi_description(inchikey):
    search_url = f"https://www.ebi.ac.uk/webservices/chebi/2.0/test/getLiteEntity?search={inchikey}&searchCategory=INCHIKEY"
    try:
        r = requests.get(search_url, timeout=5)
        if r.ok:
            entities = r.json().get("ListElement", [])
            if entities:
                chebi_id = entities[0]["chebiId"]
                entity_url = f"https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId={chebi_id}"
                r2 = requests.get(entity_url, timeout=5)
                if r2.ok:
                    return r2.json().get("chebiAsciiName")
    except requests.RequestException:
        pass
    return None

# ------------------------------
# Wikipedia summary
# ------------------------------
def fetch_wikipedia_summary(name):
    url = f"https://en.wikipedia.org/api/rest_v1/page/summary/{name}"
    try:
        r = requests.get(url, timeout=5)
        if r.ok:
            data = r.json()
            extract = data.get("extract")
            if extract and len(extract.split()) >= 20:
                return extract
    except requests.RequestException:
        pass
    return None

# ------------------------------
# Fallback description
# ------------------------------
def generate_description(name, mw, xlogp, tpsa):
    parts = [f"{name} is an organic compound"]
    if mw:
        parts.append(f"with a molecular weight of {mw:.2f} g/mol")
    if xlogp is not None:
        parts.append(f"and moderate hydrophobicity (XLogP ≈ {xlogp})")
    if tpsa is not None:
        parts.append(f"with a topological polar surface area of {tpsa:.1f} Å²")
    return " ".join(parts) + "."

# ------------------------------
# MASTER FUNCTION
# ------------------------------
def analyze_molecule(xyz_path=None, smiles_value=None):
    """Analyze a molecule (XYZ or SMILES). Returns descriptors, hazards, description."""
    if smiles_value:
        smiles = canonicalize_smiles(smiles_value)
    elif xyz_path:
        smiles = xyz_to_smiles(xyz_path)
        smiles = canonicalize_smiles(smiles)
    else:
        raise ValueError("Must provide either xyz_path or smiles_value")

    # --- Structure intelligence ---
    intel = analyze_structure(smiles)

    # Defaults (always defined)
    cid = None
    mw = None
    xlogp = None
    tpsa = None
    synonyms = []

    # --- PubChem lookup (only for valid organic compounds) ---
    compound = None
    if intel["Valid"] and not intel["HasMetal"]:
        try:
            compounds = pcp.get_compounds(smiles, namespace="smiles")
            compound = compounds[0] if compounds else None
        except Exception:
            compound = None

    if compound:
        cid = compound.cid
        mw = compound.molecular_weight
        xlogp = compound.xlogp
        tpsa = compound.tpsa
        synonyms = compound.synonyms or []

    # Hazards
    try:
        ghs = hazards.classify_hazards(cid=cid, smiles=smiles, molecular_weight=mw, xlogp=xlogp)
    except Exception:
        ghs = {"Source": "Estimated", "SignalWord": None, "GHSCodes": [], "HCodeDescriptions": {}, "Pictograms": []}

    try:
        practical = hazards.estimate_practical_hazards(smiles, intel=intel)
    except Exception:
        practical = {"Flammable": None, "Toxicity": None, "EnvironmentalHazard": None, "Advice": None}

    hazards_summary = {
        "Source": ghs.get("Source", "Estimated"),
        "SignalWord": ghs.get("SignalWord"),
        "GHS": {
            "Codes": ghs.get("GHSCodes", []),
            "Descriptions": ghs.get("HCodeDescriptions", {}),
            "Pictograms": ghs.get("Pictograms", []),
        },
        "Practical": practical,
        "StructureIntelligence": intel,
    }

    # IUPAC
    iupac = compound.iupac_name if compound and compound.iupac_name else smiles_to_iupac(smiles) or "Unknown compound"

    # Description
    description = (
        fetch_wikipedia_summary(iupac)
        or (compound and fetch_chebi_description(compound.inchikey))
        or (cid and fetch_pubchem_description(cid))
        or generate_description(iupac, mw, xlogp, tpsa)
    )

    return {
        "SMILES": smiles,
        "CID": cid,
        "IUPACName": iupac,
        "MolecularWeight": mw,
        "TPSA": tpsa,
        "XLogP": xlogp,
        "Synonyms": synonyms,
        "Description": description,
        "Hazards": hazards_summary,
        "StructureIntelligence": intel,
    }
