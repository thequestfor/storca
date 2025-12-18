"""
hazards.py

Hazard classification module.

- Uses PubChem PUG-View GHS data when available.
- Falls back to estimated practical hazards if PubChem data unavailable.
- Practical hazards now cover aromatics, ethers, alcohols, aldehydes, halogens, etc.
"""

from __future__ import annotations
import re
import requests
from typing import Dict, Optional, Set
from rdkit import Chem

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
      - Flammable/Toxicity/EnvironmentalHazard
      - Advice
    """

    # --- Attempt PubChem GHS extraction ---
    if cid:
        data = _fetch_pubchem_pugview(cid)
        if data and _has_ghs_section(data):
            hazards = _extract_pubchem_ghs(data)
            if hazards:
                # Merge practical advice instead of overwriting
                practical = estimate_practical_hazards(smiles)
                hazards.update({
                    "Flammable": practical["Flammable"],
                    "Toxicity": practical["Toxicity"],
                    "EnvironmentalHazard": practical["EnvironmentalHazard"],
                    "Advice": practical["Advice"]
                })
                return hazards

    # --- Fallback estimation ---
    hazards = estimate_practical_hazards(smiles)
    hazards.update({
        "Source": "Estimated",
        "SignalWord": None,
        "GHSCodes": [],
        "HCodeDescriptions": {},
        "Pictograms": []
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
    hazard_texts: Set[str] = set()
    signal_words: Set[str] = set()
    pictograms: Set[str] = set()

    def walk(node):
        if isinstance(node, dict):
            name = node.get("Name")
            value = node.get("Value")

            if name == "GHS Hazard Statements" and value:
                for item in value.get("StringWithMarkup", []):
                    text = item.get("String", "").strip()
                    if text:
                        hazard_texts.add(text)

            elif name == "Signal" and value:
                for item in value.get("StringWithMarkup", []):
                    sw = item.get("String", "").strip()
                    if sw:
                        signal_words.add(sw)

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
    if "Danger" in words:
        return "Danger"
    if "Warning" in words:
        return "Warning"
    return None


# ===========================
# Practical hazard estimation
# ===========================

def estimate_practical_hazards(smiles: str, intel: dict = None) -> dict:
    """Returns practical lab-oriented hazard summary."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "Flammable": None,
            "Toxicity": "Unknown",
            "EnvironmentalHazard": None,
            "Advice": "Could not parse SMILES. Exercise general lab safety precautions."
        }

    hazards = {
        "Flammable": False,
        "Toxicity": "Low",
        "EnvironmentalHazard": False,
        "Advice": "Use standard PPE (gloves, goggles), avoid inhalation."
    }

    # --- Flammable heuristics ---
    smarts_flammable = ["[CH3]", "[CH2]", "[OH]", "c", "[O][C]"]
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

    # --- Environmental hazards ---
    smarts_env = ["F", "Cl", "Br", "I"]
    for s in smarts_env:
        query = Chem.MolFromSmarts(s)
        if query and mol.HasSubstructMatch(query):
            hazards["EnvironmentalHazard"] = True
            break

    # --- Metal compounds or special functional groups ---
    if intel:
        if intel.get("HasMetal"):
            hazards["Toxicity"] = "High"
            hazards["Advice"] += " Treat as potentially toxic metal compound."
        if any(atom in intel.get("Atoms", []) for atom in ["O"]):
            # crude oxidizer detection
            if ">=O" in smiles or "=O)=O" in smiles:
                hazards["Advice"] += " Strong oxidizer – avoid contact with organics."

    # --- Tailor advice ---
    advice_list = []
    if hazards["Flammable"]:
        advice_list.append("Keep away from open flame or heat.")
    if hazards["Toxicity"] in ["Moderate", "High"]:
        advice_list.append("Avoid ingestion and skin contact.")
    if hazards["EnvironmentalHazard"]:
        advice_list.append("Prevent environmental release.")

    if advice_list:
        hazards["Advice"] = " ".join(advice_list)

    return hazards
