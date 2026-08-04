"""Evidence-first molecular persistence dossiers; never an optimizer-only verdict."""

from __future__ import annotations

import json
from pathlib import Path

from .describe import describe_smiles


AMBIENT_DEFAULTS = {
    "temperature_K": 298.15,
    "pressure_bar": 1.0,
    "atmosphere": "air",
    "phase": "ordinary condensed-phase handling",
    "persistence_target": "at least 24 hours",
}


def create_plausibility_dossier(
    smiles: str,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    conditions: dict | None = None,
) -> dict:
    """Create a validated request dossier without overstating persistence.

    This first stage intentionally makes no stability prediction. It establishes
    the chemical state and records exactly which RMG/ORCA evidence is still
    required before any non-abstaining category can be assigned.
    """
    if multiplicity < 1:
        raise ValueError("Multiplicity must be at least one")
    from rdkit import Chem

    description = describe_smiles(smiles)
    molecule = Chem.MolFromSmiles(description["canonical_smiles"])
    electron_count = sum(atom.GetAtomicNum() for atom in molecule.GetAtoms()) - charge
    parity_consistent = (electron_count % 2 == 0) == (multiplicity % 2 == 1)
    resolved_conditions = {**AMBIENT_DEFAULTS, **(conditions or {})}
    return {
        "schema_version": 1,
        "kind": "molecular_persistence_dossier",
        "species": {
            "canonical_smiles": description["canonical_smiles"],
            "charge": charge,
            "multiplicity": multiplicity,
            "electron_count": electron_count,
            "spin_parity_consistent": parity_consistent,
        },
        "conditions": resolved_conditions,
        "evidence": [{
            "kind": "local_structure_validation",
            "status": "passed",
            "source": "RDKit local structure parsing",
            "details": {"description": description, "spin_parity_consistent": parity_consistent},
        }],
        "assessment": {
            "persistence_category": "insufficient_evidence",
            "reason": "A valid structure and spin-parity check do not establish a local minimum, decomposition kinetics, or ambient persistence.",
        },
        "required_evidence": [
            "ORCA conformer/optimization/frequency evidence with identity-retention check",
            "RMG condition-specific pathway search with retained coverage and logs",
            "ORCA verification of leading decomposition energetics or barriers when RMG identifies material risk",
        ],
        "provenance": {"network_accessed": False, "source": "STORCA local plausibility request"},
    }


def attach_rmg_evidence(dossier: dict, evidence: dict) -> dict:
    """Attach bounded RMG evidence while preserving an abstaining assessment."""
    dossier["evidence"].append(evidence)
    risk_flags = dossier.setdefault("risk_flags", [])
    if evidence.get("search_outcome") == "candidate_decomposition_pathway_identified":
        risk_flags.append("candidate_decomposition_pathway_identified")
    elif evidence.get("status") == "incomplete":
        risk_flags.append("incomplete_rmg_search")
    elif evidence.get("status") == "failed":
        risk_flags.append("rmg_unavailable_or_failed")
    else:
        risk_flags.append("no_route_found_within_bounded_rmg_model")
    dossier["assessment"]["persistence_category"] = "insufficient_evidence"
    dossier["assessment"]["reason"] = "RMG evidence is a bounded pathway screen and cannot establish persistence without ORCA structure and pathway verification."
    return dossier


def write_plausibility_dossier(path: Path, dossier: dict) -> Path:
    path = Path(path)
    path.write_text(json.dumps(dossier, indent=2, sort_keys=True) + "\n")
    return path
