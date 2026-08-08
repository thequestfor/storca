"""Materialize small, condition-relevant subsets of installed RMG libraries."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess


def build_rmg_reference_subset(
    target_smiles: str,
    declared_reactive_smiles: list[str],
    library_names: list[str],
    output_dir: Path,
    *,
    rmg_env: str | None,
    timeout_seconds: float = 300.0,
) -> dict:
    """Extract directly reachable reference reactions in the RMG runtime.

    A complete combustion library must not be injected into a bounded storage
    screen: doing so admits unrelated species and can violate RMG species
    constraints before generation starts.  The companion script imports RMG
    only inside its selected Conda environment and saves a normal, tiny local
    kinetics library.
    """
    names = list(dict.fromkeys(name.strip() for name in library_names if name.strip()))
    if not names:
        return {
            "schema_version": 1,
            "kind": "storca_rmg_reference_subset",
            "status": "not_requested",
            "source_libraries": [],
            "selected_reaction_count": 0,
        }
    output_dir = Path(output_dir).resolve()
    script = Path(__file__).with_name("rmg_reference_subset.py").resolve()
    command = (["conda", "run", "-n", rmg_env, "python", str(script)]
               if rmg_env else ["python", str(script)])
    command += ["--target-smiles", target_smiles, "--output", str(output_dir)]
    for smiles in list(dict.fromkeys([target_smiles, *declared_reactive_smiles])):
        command += ["--declared-smiles", smiles]
    for name in names:
        command += ["--library", name]
    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise RuntimeError(
            "RMG reference-library subset extraction failed in "
            f"{rmg_env or 'the active environment'}: {detail[-2000:]}"
        )
    manifest_file = output_dir / "reference-subset.json"
    if not manifest_file.is_file():
        raise RuntimeError("RMG reference-library extractor did not write its manifest")
    manifest = json.loads(manifest_file.read_text())
    if manifest.get("source_libraries") != names:
        raise RuntimeError("RMG reference-library manifest does not match the requested libraries")
    required = (output_dir / "reactions.py", output_dir / "dictionary.txt")
    if not all(path.is_file() for path in required):
        raise RuntimeError("RMG reference-library extractor produced an incomplete local library")
    manifest["status"] = "completed"
    manifest["runtime_environment"] = rmg_env or "active"
    manifest["artifacts"] = {
        "manifest": str(manifest_file),
        "reactions": str(required[0]),
        "dictionary": str(required[1]),
    }
    return manifest
