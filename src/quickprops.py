# src/quickprops.py
import os
import requests
import pubchempy as pcp
from openbabel import openbabel as ob
from rdkit import Chem

from . import hazards
from .chem_intel import analyze_structure

GOOD_SECTIONS = {
    "description",
    "physical description",
    "chemical and physical properties",
    "safety and hazards",
    "uses",
}

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
# Canonical SMILES
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
# Wikipedia summary
# ------------------------------
def fetch_wikipedia_summary(name):
    url = f"https://en.wikipedia.org/api/rest_v1/page/summary/{name}"
    try:
        r = requests.get(url, timeout=5)
        if r.ok:
            extract = r.json().get("extract")
            if extract and len(extract.split()) >= 20:
                return extract
    except requests.RequestException:
        pass
    return None


# ------------------------------
# Generated description
# ------------------------------
def generate_description(name, mw, xlogp, tpsa, intel=None, hazards=None):
    compound_class = intel.get("CompoundClass", "compound") if intel else "compound"

    if compound_class == "organometallic":
        parts = [f"{name} is an organometallic compound"]
    elif compound_class == "inorganic":
        parts = [f"{name} is an inorganic compound"]
    else:
        parts = [f"{name} is an organic compound"]

    if mw:
        parts.append(f"with a molecular weight of {mw:.2f} g/mol")

    if xlogp is not None:
        parts.append(f"and an estimated XLogP of {xlogp}")

    if tpsa is not None:
        parts.append(f"(TPSA ≈ {tpsa:.1f} Å²)")

    toxicity = hazards.get("Practical", {}).get("Toxicity") if hazards else None
    if toxicity in {"Moderate (systemic)", "High (systemic)"}:
        parts.append("and requires appropriate safety precautions during handling")

    return " ".join(parts) + "."


# ------------------------------
# MASTER FUNCTION
# ------------------------------
def analyze_molecule(xyz_path=None, smiles_value=None):
    """Analyze a molecule (XYZ or SMILES)."""

    # --- Resolve SMILES ---
    if smiles_value:
        smiles = canonicalize_smiles(smiles_value)
    elif xyz_path:
        smiles = canonicalize_smiles(xyz_to_smiles(xyz_path))
    else:
        raise ValueError("Must provide either xyz_path or smiles_value")

    # --- Structural intelligence ---
    intel = analyze_structure(smiles)

    # --- Defaults ---
    cid = mw = xlogp = tpsa = None
    synonyms = []
    compound = None

    # --- PubChem lookup ---
    try:
        results = pcp.get_compounds(smiles, namespace="smiles")
        if results:
            compound = results[0]
            cid = compound.cid
            mw = compound.molecular_weight
            xlogp = compound.xlogp
            tpsa = compound.tpsa
            synonyms = compound.synonyms or []
    except Exception:
        pass

    # --- IUPAC name ---
    iupac = (
        compound.iupac_name
        if compound and compound.iupac_name
        else smiles_to_iupac(smiles)
        or "Unknown compound"
    )

    # --- Hazards (single source of truth) ---
    hazards_result = hazards.classify_hazards(
        cid=cid,
        smiles=smiles,
        molecular_weight=mw,
        xlogp=xlogp,
    )

    hazards_summary = {
        "Source": hazards_result.get("Source"),
        "SignalWord": hazards_result.get("SignalWord"),
        "GHS": {
            "Codes": hazards_result.get("GHSCodes", []),
            "Descriptions": hazards_result.get("HCodeDescriptions", {}),
            "Pictograms": hazards_result.get("Pictograms", []),
        },
        "Practical": {
            "Flammable": hazards_result.get("Flammable"),
            "Toxicity": hazards_result.get("Toxicity"),
            "PhysicalHazard": hazards_result.get("PhysicalHazard"),
            "EnvironmentalHazard": hazards_result.get("EnvironmentalHazard"),
            "Advice": hazards_result.get("Advice"),
        },
        "StructureIntelligence": intel,
    }

    # --- Description cascade ---
    description = description_source = None

    if cid:
        description = fetch_pubchem_description(cid)
        if description:
            description_source = "PubChem"

    if not description:
        description = fetch_wikipedia_summary(iupac)
        if description:
            description_source = "Wikipedia"

    if not description:
        description = generate_description(
            name=iupac,
            mw=mw,
            xlogp=xlogp,
            tpsa=tpsa,
            intel=intel,
            hazards=hazards_summary,
        )
        description_source = "Generated"

    return {
        "SMILES": smiles,
        "CID": cid,
        "IUPACName": iupac,
        "MolecularWeight": mw,
        "TPSA": tpsa,
        "XLogP": xlogp,
        "Synonyms": synonyms,
        "Description": description,
        "DescriptionSource": description_source,
        "Hazards": hazards_summary,
        "StructureIntelligence": intel,
    }
