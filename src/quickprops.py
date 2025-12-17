import json
import os
import requests
import pubchempy as pcp
from openbabel import openbabel as ob

# ============================================================
# XYZ → SMILES
# ============================================================
def xyz_to_smiles(xyz_file):
    xyz_file = str(xyz_file)
    if not os.path.exists(xyz_file):
        raise FileNotFoundError(f"XYZ file not found: {xyz_file}")

    conv = ob.OBConversion()
    conv.SetInFormat("xyz")
    conv.SetOutFormat("smi")

    mol = ob.OBMol()
    if not conv.ReadFile(mol, xyz_file):
        raise ValueError(f"Failed to read XYZ file: {xyz_file}")

    smiles = conv.WriteString(mol).strip()
    if not smiles:
        raise ValueError("Open Babel failed to generate SMILES")

    return smiles


# ============================================================
# SMILES → IUPAC (NIH CACTUS fallback)
# ============================================================
def smiles_to_iupac(smiles):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/iupac_name"
    try:
        r = requests.get(url, timeout=5)
        if r.ok and r.text.strip():
            return r.text.strip()
    except requests.RequestException:
        pass
    return None


# ============================================================
# PUBCHEM PUG-VIEW (SECTION-AWARE EXTRACTION)
# ============================================================
GOOD_SECTIONS = {
    "description",
    "physical description",
    "chemical and physical properties",
    "safety and hazards",
    "uses"
}


def extract_pubchem_description(data, max_paragraphs=3):
    collected = []

    def walk(node, current_heading=None):
        if isinstance(node, dict):
            heading = node.get("TOCHeading", current_heading)
            heading_norm = heading.lower() if isinstance(heading, str) else None

            if heading_norm in GOOD_SECTIONS:
                if "StringWithMarkup" in node:
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
        if not r.ok:
            return None
        return extract_pubchem_description(r.json())
    except requests.RequestException:
        return None


# ============================================================
# ChEBI (CURATED DESCRIPTIONS)
# ============================================================
def fetch_chebi_description(inchikey):
    search_url = (
        "https://www.ebi.ac.uk/webservices/chebi/2.0/test/"
        f"getLiteEntity?search={inchikey}&searchCategory=INCHIKEY"
    )
    try:
        r = requests.get(search_url, timeout=5)
        if not r.ok:
            return None
        entities = r.json().get("ListElement", [])
        if not entities:
            return None

        chebi_id = entities[0]["chebiId"]
        entity_url = (
            "https://www.ebi.ac.uk/webservices/chebi/2.0/test/"
            f"getCompleteEntity?chebiId={chebi_id}"
        )
        r2 = requests.get(entity_url, timeout=5)
        if r2.ok:
            return r2.json().get("chebiAsciiName")
    except requests.RequestException:
        pass

    return None


# ============================================================
# Wikipedia Summary
# ============================================================
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


# ============================================================
# Generated Fallback Description
# ============================================================
def generate_description(name, mw, xlogp, tpsa):
    parts = [f"{name} is an organic compound"]
    if mw:
        parts.append(f"with a molecular weight of {mw:.2f} g/mol")
    if xlogp is not None:
        parts.append(f"and moderate hydrophobicity (XLogP ≈ {xlogp})")
    if tpsa is not None:
        parts.append(f"with a topological polar surface area of {tpsa:.1f} Å²")
    return " ".join(parts) + "."


# ============================================================
# MASTER FUNCTION (API UNCHANGED)
# ============================================================
def analyze_xyz(xyz_path):
    smiles = xyz_to_smiles(xyz_path)

    compounds = pcp.get_compounds(smiles, namespace="smiles")
    compound = compounds[0] if compounds else None

    cid = compound.cid if compound else None
    mw = compound.molecular_weight if compound else None
    xlogp = compound.xlogp if compound else None
    tpsa = compound.tpsa if compound else None
    synonyms = compound.synonyms if compound and compound.synonyms else []

    iupac = (
        compound.iupac_name
        if compound and compound.iupac_name
        else smiles_to_iupac(smiles)
        or "Unknown compound"
    )

    description = None

    # 1. Wikipedia (best prose)
    description = fetch_wikipedia_summary(iupac)

    # 2. ChEBI (curated)
    if not description and compound and compound.inchikey:
        description = fetch_chebi_description(compound.inchikey)

    # 3. PubChem narrative
    if not description and cid:
        description = fetch_pubchem_description(cid)

    # 4. Generated fallback
    if not description:
        description = generate_description(iupac, mw, xlogp, tpsa)

    return {
        "SMILES": smiles,
        "CID": cid,
        "IUPACName": iupac,
        "MolecularWeight": mw,
        "TPSA": tpsa,
        "XLogP": xlogp,
        "Synonyms": synonyms,
        "Description": description,
    }


# ============================================================
# CLI TEST
# ============================================================
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python quickprops.py <file.xyz>")
        sys.exit(1)

    info = analyze_xyz(sys.argv[1])
    print(json.dumps(info, indent=2))

