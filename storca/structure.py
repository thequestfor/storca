"""Deterministic local structure artifacts for STORCA molecule records."""

from __future__ import annotations

import json
from pathlib import Path

from .describe import describe_smiles


def build_structure_artifacts(smiles: str, output_dir: Path) -> dict:
    """Write a 2D depiction, a deterministic 3D seed geometry, and JSON record.

    The XYZ geometry is an RDKit/UFF seed for subsequent calculations, not a
    calculated equilibrium geometry or a structure determination.
    """
    from src.molecule_tools import sanitize_smiles, smiles_to_png, smiles_to_xyz

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    canonical_smiles = sanitize_smiles(smiles)
    description = describe_smiles(canonical_smiles)
    xyz_path = smiles_to_xyz(canonical_smiles, output_dir / "initial.xyz")
    image_path = smiles_to_png(canonical_smiles, output_dir / "structure.png")
    result = {
        "schema_version": 1,
        "kind": "local_structure_record",
        "description": description,
        "artifacts": {
            "initial_xyz": str(xyz_path),
            "structure_png": str(image_path),
        },
        "coordinate_provenance": {
            "source": "RDKit conformer embedding followed by UFF minimization",
            "random_seed": 42,
            "intended_use": "initial geometry only; optimize with ORCA before interpreting geometry-dependent properties",
        },
        "network_accessed": False,
    }
    record_path = output_dir / "molecule.json"
    record_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    result["record_json"] = str(record_path)
    return result
