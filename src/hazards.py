"""
hazards.py

Hazard classification module.

- Uses PubChem PUG-View GHS data when available.
- Falls back to structural heuristics if PubChem data unavailable.
- Maps GHS H-codes to practical hazard flags: Flammable, Toxicity, Corrosive, EnvironmentalHazard, PhysicalHazard.
"""

from __future__ import annotations
import re
import requests
from typing import Dict, Optional, Set
from rdkit import Chem

PUBCHEM_PUGVIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"

# --- Mapping of H-codes to practical hazard flags ---
H_TO_PRACTICAL = {
    # Flammable
    "H220": "Flammable",
    "H221": "Flammable",
    "H222": "Flammable",
    "H223": "Flammable",
    "H224": "Flammable",
    "H225": "Flammable",
    "H226": "Flammable",
    "H227": "Flammable",

    # Acute toxicity
    "H300": "ToxicityHigh",
    "H301": "ToxicityHigh",
    "H302": "ToxicityModerate",
    "H310": "ToxicityHigh",
    "H311": "ToxicityHigh",
    "H312": "ToxicityModerate",
    "H330": "ToxicityHigh",
    "H331": "ToxicityHigh",
    "H332": "ToxicityModerate",
    "H333": "ToxicityModerate",

    # Corrosive / Irritation
    "H314": "Corrosive",
    "H315": "Corrosive",
    "H318": "Corrosive",
    "H319": "Corrosive",

    # Environmental hazards
    "H400": "EnvironmentalHazard",
    "H401": "EnvironmentalHazard",
    "H410": "EnvironmentalHazard",
    "H411": "EnvironmentalHazard",
}

# ===========================
# Public API
# ===========================

def classify_hazards(
    *,
    cid: Optional[int],
    smiles: str,
    molecular_weight: Optional[float] = None,
    xlogp: Optional[float] = None,
) -> Dict:
    """
    Return hazard dictionary with:
        - GHS H-codes, signal words, pictograms
        - Practical lab flags: Flammable, Toxicity, Corrosive, EnvironmentalHazard
        - Tailored Advice
    """
    # Start with structural heuristics
    practical = estimate_practical_hazards(smiles)

    # Attempt PubChem GHS extraction
    ghs = None
    if cid:
        data = _fetch_pubchem_pugview(cid)
        if data and _has_ghs_section(data):
            ghs = _extract_pubchem_ghs(data)

    # Apply H-code mapping to practical hazards
    if ghs:
        for hcode in ghs.get("GHSCodes", []):
            flag = H_TO_PRACTICAL.get(hcode)
            if flag == "Flammable":
                practical["Flammable"] = True
                practical["PhysicalHazard"] = max_hazard(practical.get("PhysicalHazard"), "High")
            elif flag == "ToxicityHigh":
                practical["Toxicity"] = max_toxicity(practical.get("Toxicity"), "High (systemic)")
            elif flag == "ToxicityModerate":
                practical["Toxicity"] = max_toxicity(practical.get("Toxicity"), "Moderate (systemic)")
            elif flag == "Corrosive":
                practical["PhysicalHazard"] = max_hazard(practical.get("PhysicalHazard"), "High")
                practical["Advice"] += " Handle with gloves and eye protection."
            elif flag == "EnvironmentalHazard":
                practical["EnvironmentalHazard"] = True

        # Merge GHS info with practical flags
        result = {**ghs, **practical}
    else:
        # No PubChem GHS, fallback
        result = {
            "Source": "Estimated",
            "SignalWord": None,
            "GHSCodes": [],
            "HCodeDescriptions": {},
            "Pictograms": [],
            **practical
        }

    # Build advice based on updated flags
    advice = []
    if result.get("Flammable"):
        advice.append("Keep away from heat, sparks, and open flame.")
    if result.get("PhysicalHazard") in ["High", "Extreme"]:
        advice.append("Use in well-ventilated area with appropriate precautions.")
    if result.get("EnvironmentalHazard"):
        advice.append("Prevent environmental release.")
    if "Corrosive" in H_TO_PRACTICAL.values() and "Handle with gloves" in result.get("Advice", ""):
        advice.append("Handle with gloves and eye protection.")
    if advice:
        result["Advice"] = " ".join(advice)

    return result

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
    def walk(node):
        if isinstance(node, dict):
            heading = node.get("TOCHeading", "").lower()
            if "hazard" in heading or "ghs" in heading:
                return True
            return any(walk(v) for v in node.values())
        if isinstance(node, list):
            return any(walk(i) for i in node)
        return False
    return walk(data)

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
                    txt = item.get("String", "").strip()
                    if txt:
                        hazard_texts.add(txt)

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

    h_codes = sorted(set(re.findall(r"H\d{3}", " ".join(hazard_texts))))
    hcode_map = {}
    for text in hazard_texts:
        m = re.match(r"(H\d{3}):?\s*(.*)", text)
        if m:
            h, desc = m.groups()
            hcode_map[h] = desc

    return {
        "Source": "PubChem GHS",
        "SignalWord": _select_signal_word(signal_words),
        "GHSCodes": h_codes,
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
# Practical hazard estimation (structural fallback)
# ===========================

def estimate_practical_hazards(smiles: str, intel: dict = None) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "Flammable": None,
            "Toxicity": "Unknown",
            "PhysicalHazard": "Unknown",
            "EnvironmentalHazard": None,
            "Advice": "Could not parse SMILES. Exercise general lab safety precautions."
        }

    hazards = {
        "Flammable": False,
        "Toxicity": "Low (systemic)",
        "PhysicalHazard": "Low",
        "EnvironmentalHazard": False,
        "Advice": "Use standard PPE (gloves, goggles), avoid inhalation."
    }

    atoms = {a.GetSymbol() for a in mol.GetAtoms()}

    # --- Inorganic flammable gases ---
    inorganic_flammable = {"H", "Si", "B", "P", "As"}
    if atoms <= inorganic_flammable:
        hazards["Flammable"] = True
        hazards["PhysicalHazard"] = "Extreme"

    # --- Organic flammability heuristics ---
    smarts_flammable = ["[CH3]", "[CH2]", "[OH]", "c", "[O][C]"]
    for s in smarts_flammable:
        q = Chem.MolFromSmarts(s)
        if q and mol.HasSubstructMatch(q):
            hazards["Flammable"] = True
            hazards["PhysicalHazard"] = max_hazard(hazards.get("PhysicalHazard"), "High")
            break

    # --- Toxicity heuristics ---
    smarts_toxic = ["[N+]", "[C#N]", "N=O"]
    for s in smarts_toxic:
        q = Chem.MolFromSmarts(s)
        if q and mol.HasSubstructMatch(q):
            hazards["Toxicity"] = max_toxicity(hazards.get("Toxicity"), "Moderate (systemic)")
            break

    # --- Environmental hazards ---
    for hal in ["F", "Cl", "Br", "I"]:
        if any(a.GetSymbol() == hal for a in mol.GetAtoms()):
            hazards["EnvironmentalHazard"] = True
            break

    # --- Tailored advice ---
    advice = []
    if hazards["Flammable"]:
        advice.append("Keep away from heat, sparks, and open flame.")
    if hazards["PhysicalHazard"] in ["High", "Extreme"]:
        advice.append("Use in well-ventilated area with appropriate precautions.")
    if hazards["EnvironmentalHazard"]:
        advice.append("Prevent environmental release.")

    if advice:
        hazards["Advice"] = " ".join(advice)

    return hazards

# ===========================
# Utilities
# ===========================

def max_hazard(current: str, new: str) -> str:
    """Return the more severe hazard level."""
    levels = ["Low", "Moderate", "High", "Extreme"]
    try:
        return max([current, new], key=lambda x: levels.index(x))
    except ValueError:
        return new

def max_toxicity(current: str, new: str) -> str:
    """Return the more severe toxicity level."""
    levels = ["Low (systemic)", "Moderate (systemic)", "High (systemic)"]
    try:
        return max([current, new], key=lambda x: levels.index(x))
    except ValueError:
        return new
