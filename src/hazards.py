"""
hazards.py

Hazard classification module.

Phase 1: fully operational GHS extraction with H-code explanations.
- Uses PubChem PUG-View GHS data when available
- Falls back to estimated hazards if PubChem data unavailable
"""

from __future__ import annotations
import re
import requests
from typing import Dict, Optional, Set


PUBCHEM_PUGVIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"


# ===========================
# Public API
# ===========================

def classify_hazards(
    *,
    cid: Optional[int],
    smiles: str,
    molecular_weight: Optional[float],
    xlogp: Optional[float],
) -> Dict:
    """
    Determine hazard classification for a compound.

    Returns a structured hazard dictionary with:
      - Source
      - SignalWord
      - GHSCodes
      - HazardStatements
      - HCodeDescriptions
      - Pictograms
      - Flammable/Toxicity/EnvironmentalHazard (if estimated)
    """

    pubchem_checked = False
    pubchem_available = False

    # --- Attempt PubChem GHS extraction ---
    if cid:
        pubchem_checked = True
        data = _fetch_pubchem_pugview(cid)
        if data:
            pubchem_available = _has_ghs_section(data)
            hazards = _extract_pubchem_ghs(data)
            if hazards:
                hazards.update({
                    "PubChemChecked": True,
                    "PubChemGHSAvailable": True,
                    "Flammable": None,
                    "Toxicity": None,
                    "EnvironmentalHazard": None
                })
                return hazards

    # --- Fallback estimation ---
    hazards = _estimate_hazards(
        smiles=smiles,
        molecular_weight=molecular_weight,
        xlogp=xlogp,
    )
    hazards.update({
        "PubChemChecked": pubchem_checked,
        "PubChemGHSAvailable": pubchem_available,
    })
    return hazards


# ===========================
# PubChem PUG-View Parsing
# ===========================

def _fetch_pubchem_pugview(cid: int) -> Optional[dict]:
    try:
        r = requests.get(PUBCHEM_PUGVIEW.format(cid=cid), timeout=5)
        if r.ok:
            return r.json()
    except requests.RequestException:
        pass
    return None


def _has_ghs_section(data: dict) -> bool:
    """
    Check if a PubChem PUG-View JSON contains a GHS section.
    """
    found = False

    def walk(node):
        nonlocal found
        if isinstance(node, dict):
            heading = node.get("TOCHeading", "").lower()
            if "hazard" in heading or "ghs" in heading:
                found = True
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(data)
    return found


def _extract_pubchem_ghs(data: dict) -> Optional[Dict]:
    """
    Extract GHS hazard data from PubChem PUG-View JSON.
    Includes H-code descriptions.
    """
    hazard_texts: Set[str] = set()
    signal_words: Set[str] = set()
    pictograms: Set[str] = set()

    def walk(node):
        if isinstance(node, dict):
            name = node.get("Name")
            value = node.get("Value")

            # Extract full hazard statements
            if name == "GHS Hazard Statements" and value:
                for item in value.get("StringWithMarkup", []):
                    text = item.get("String", "").strip()
                    if text:
                        hazard_texts.add(text)

            # Signal words (Danger/Warning)
            elif name == "Signal" and value:
                for item in value.get("StringWithMarkup", []):
                    sw = item.get("String", "").strip()
                    if sw:
                        signal_words.add(sw)

            # Pictograms (Irritant, Health Hazard, Environmental)
            elif name == "Pictogram(s)" and value:
                for item in value.get("StringWithMarkup", []):
                    for markup in item.get("Markup", []):
                        if markup.get("Type") == "Icon":
                            pictograms.add(markup.get("Extra"))

            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(data)

    if not hazard_texts:
        return None

    # Extract H-codes
    h_codes = sorted(set(re.findall(r"H\d{3}", " ".join(hazard_texts))))

    # Map H-codes → full description
    hcode_map = {}
    for text in hazard_texts:
        match = re.match(r"(H\d{3}):?\s*(.*)", text)
        if match:
            h, desc = match.groups()
            hcode_map[h] = desc

    return {
        "Source": "PubChem GHS",
        "SignalWord": _select_signal_word(signal_words),
        "GHSCodes": h_codes,
        "HazardStatements": sorted(hazard_texts),
        "HCodeDescriptions": hcode_map,
        "Pictograms": sorted(pictograms),
    }


def _select_signal_word(words: Set[str]) -> Optional[str]:
    """
    Prefer 'Danger' over 'Warning' if both appear.
    """
    if "Danger" in words:
        return "Danger"
    if "Warning" in words:
        return "Warning"
    return None

def estimate_practical_hazards(smiles):
    """
    Returns a practical lab-oriented hazard summary.
    Flammable, Toxicity, Environmental Hazard, and short advice.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    hazards = {
        "Flammable": False,
        "Toxicity": "Low",
        "EnvironmentalHazard": False,
        "Advice": "Use standard PPE (gloves, goggles), avoid inhalation."
    }

    # Rough heuristics:
    smarts_flammable = ["[OH]", "[CH3]", "[CH2]"]  # simple example
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in smarts_flammable):
        hazards["Flammable"] = True

    # Functional groups that could increase toxicity
    smarts_toxic = ["[N+]", "[C#N]", "N=O"]
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in smarts_toxic):
        hazards["Toxicity"] = "Moderate"

    # Environmental hazards: halogens
    smarts_env = ["F", "Cl", "Br", "I"]
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in smarts_env):
        hazards["EnvironmentalHazard"] = True

    # Tailor advice
    advice = []
    if hazards["Flammable"]:
        advice.append("Keep away from open flame or heat.")
    if hazards["Toxicity"] in ["Moderate", "High"]:
        advice.append("Avoid ingestion and skin contact.")
    if hazards["EnvironmentalHazard"]:
        advice.append("Prevent environmental release.")
    if advice:
        hazards["Advice"] = " ".join(advice)

    return hazards

# ===========================
# Fallback Estimation
# ===========================

def estimate_practical_hazards(smiles):
    """
    Returns a practical lab-oriented hazard summary.
    Flammable, Toxicity, Environmental Hazard, and short advice.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # If RDKit cannot parse the SMILES, return conservative defaults
        return {
            "Flammable": False,
            "Toxicity": "Unknown",
            "EnvironmentalHazard": False,
            "Advice": "SMILES could not be parsed. Exercise general laboratory safety precautions."
        }

    hazards = {
        "Flammable": False,
        "Toxicity": "Low",
        "EnvironmentalHazard": False,
        "Advice": "Use standard PPE (gloves, goggles), avoid inhalation."
    }

    # --- Flammable heuristics ---
    smarts_flammable = ["[OH]", "[CH3]", "[CH2]"]
    for s in smarts_flammable:
        query = Chem.MolFromSmarts(s)
        if query and mol.HasSubstructMatch(query):
            hazards["Flammable"] = True
            break

    # --- Toxicity heuristics ---
    smarts_toxic = ["[N+]", "[C#N]", "N=O"]
    for s in smarts_toxic:
        query = Chem.MolFromSmarts(s)
        if query and mol.HasSubstructMatch(query):
            hazards["Toxicity"] = "Moderate"
            break

    # --- Environmental hazards (halogens) ---
    smarts_env = ["F", "Cl", "Br", "I"]
    for s in smarts_env:
        query = Chem.MolFromSmarts(s)
        if query and mol.HasSubstructMatch(query):
            hazards["EnvironmentalHazard"] = True
            break

    # --- Tailor advice ---
    advice = []
    if hazards["Flammable"]:
        advice.append("Keep away from open flame or heat.")
    if hazards["Toxicity"] in ["Moderate", "High"]:
        advice.append("Avoid ingestion and skin contact.")
    if hazards["EnvironmentalHazard"]:
        advice.append("Prevent environmental release.")
    if advice:
        hazards["Advice"] = " ".join(advice)

    return hazards

