"""Explicit, source-labelled external molecule enrichment."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


PUBCHEM_PROPERTIES = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/Title,IUPACName,MolecularFormula,MolecularWeight,XLogP,TPSA/JSON"
PUBCHEM_SYNONYMS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"


def _get_json(url: str, *, timeout_seconds: float) -> dict:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "STORCA/0.1 (+local chemistry workflow)"})
    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        raise RuntimeError(f"PubChem request failed with HTTP {error.code}") from error
    except (URLError, TimeoutError) as error:
        raise RuntimeError(f"PubChem request failed: {error.reason if isinstance(error, URLError) else error}") from error


def enrich_pubchem(smiles: str, *, timeout_seconds: float = 10.0, synonym_limit: int = 10) -> dict:
    """Look up PubChem identifiers/properties only when the user invokes this API."""
    if timeout_seconds <= 0 or synonym_limit < 0:
        raise ValueError("Timeout must be positive and synonym limit cannot be negative")
    from src.molecule_tools import sanitize_smiles

    canonical_smiles = sanitize_smiles(smiles)
    property_url = PUBCHEM_PROPERTIES.format(smiles=quote(canonical_smiles, safe=""))
    payload = _get_json(property_url, timeout_seconds=timeout_seconds)
    properties = payload.get("PropertyTable", {}).get("Properties", [])
    if not properties:
        raise RuntimeError("PubChem returned no compound properties for this structure")
    record = properties[0]
    cid = record.get("CID")
    synonyms: list[str] = []
    synonym_url = None
    if cid is not None and synonym_limit:
        synonym_url = PUBCHEM_SYNONYMS.format(cid=cid)
        synonym_payload = _get_json(synonym_url, timeout_seconds=timeout_seconds)
        synonyms = synonym_payload.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])[:synonym_limit]
    return {
        "schema_version": 1,
        "kind": "external_compound_enrichment",
        "input_smiles": smiles,
        "canonical_smiles": canonical_smiles,
        "pubchem": {
            "cid": cid,
            "title": record.get("Title"),
            "iupac_name": record.get("IUPACName"),
            "molecular_formula": record.get("MolecularFormula"),
            "molecular_weight_g_mol": record.get("MolecularWeight"),
            "xlogp": record.get("XLogP"),
            "tpsa_a2": record.get("TPSA"),
            "synonyms": synonyms,
        },
        "provenance": {
            "source": "PubChem PUG REST",
            "network_accessed": True,
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "property_url": property_url,
            "synonym_url": synonym_url,
        },
        "limitations": [
            "Values are returned by an external database and may change after retrieval.",
            "This command does not provide a safety assessment, SDS, experimental property guarantee, or identity confirmation.",
        ],
    }


def write_enrichment(path: Path, result: dict) -> Path:
    path = Path(path)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return path
