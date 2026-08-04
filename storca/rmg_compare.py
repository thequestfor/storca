"""Controlled, artifact-retaining comparisons between RMG installations."""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

from .conditions import build_condition_spec
from .runs import write_metadata
from .stability import collect_rmg_evidence, resolve_stability_configuration


def rmg_version(rmg_env: str | None) -> dict:
    """Ask the selected RMG environment for its installed RMG-Py version."""
    command = (["conda", "run", "-n", rmg_env, "python", "-c", "import rmgpy; print(rmgpy.__version__)"]
               if rmg_env else ["python", "-c", "import rmgpy; print(rmgpy.__version__)"])
    result = subprocess.run(command, capture_output=True, text=True)
    return {"environment": rmg_env, "command": command, "returncode": result.returncode,
            "version": result.stdout.strip() if result.returncode == 0 else None,
            "stderr": result.stderr.strip() or None}


def compare_rmg_versions(
    smiles: str, run_dir: Path, *, rmg3_env: str = "rmg_env", rmg4_env: str = "rmg4_env",
    scenario: str = "ambient-air-gas-screen", screen_tier: str = "quick-screen",
    temperature: float = 298.15, pressure: float = 1.0, target_duration_hours: float = 24.0,
    retention_fraction: float = 0.95, rmg_walltime: str | None = None,
    rmg_max_processes: int | None = None, rmg_max_iterations: int | None = None,
    rmg_max_edge_species: int | None = None,
) -> dict:
    """Run two installations with one identical condition and resource contract.

    This is a model-comparison tool, not an accuracy ranking: database and
    solver changes can legitimately yield different mechanisms.
    """
    run_dir = Path(run_dir)
    scenario_config, resources = resolve_stability_configuration(
        scenario, screen_tier, rmg_walltime=rmg_walltime, rmg_max_processes=rmg_max_processes,
        rmg_max_iterations=rmg_max_iterations, rmg_max_edge_species=rmg_max_edge_species,
    )
    conditions = build_condition_spec(
        scenario_config, temperature_K=temperature, pressure_bar=pressure,
        target_duration_hours=target_duration_hours, retention_fraction=retention_fraction,
    )
    common = {
        "barrier_threshold": 50.0, "temperature": temperature, "pressure": pressure,
        "rmg_walltime": resources["walltime"], "rmg_max_processes": resources["max_processes"],
        "rmg_max_iterations": resources["max_iterations"], "rmg_max_edge_species": resources["max_edge_species"],
        "requested_phase": conditions.phase_model, "scenario": scenario_config, "conditions": conditions,
    }
    runs = []
    for label, environment in (("rmg3", rmg3_env), ("rmg4", rmg4_env)):
        version = rmg_version(environment)
        started = time.monotonic()
        evidence = collect_rmg_evidence(smiles, run_dir / label, rmg_env=environment, **common)
        wall_seconds = time.monotonic() - started
        profile = evidence.get("solver_profile") or {}
        runs.append({
            "label": label, "version": version, "wall_seconds": wall_seconds,
            "status": evidence["status"], "search_outcome": evidence.get("search_outcome"),
            "candidate_route_count": len(evidence.get("candidate_routes", [])),
            "species_count": len(evidence.get("species_dictionary", {})),
            "target_loss_fraction": profile.get("target_loss_fraction"),
            "target_fraction_remaining": profile.get("target_fraction_remaining"),
            "end_time_seconds": profile.get("end_time_seconds"), "artifacts": evidence["artifacts"],
            "failure_reason": evidence.get("failure_reason"),
        })
    report = {
        "schema_version": 1, "kind": "rmg_version_comparison", "smiles": smiles,
        "comparison_contract": {"conditions": conditions.as_dict(), "resource_limits": resources,
                                "same_input_policy": "Each run is independently generated from this identical contract; no generated kinetics libraries are included."},
        "runs": runs,
        "interpretation": "Compare completion, wall time, model size, target loss, and retained Chemkin artifacts. Different predicted chemistry is a result to investigate, not proof that either installation is correct.",
    }
    output = run_dir / "rmg-comparison.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    write_metadata(run_dir, workflow="rmg_version_comparison", result_json=str(output),
                   condition_contract=conditions.as_dict(), resource_limits=resources)
    return {**report, "result_json": output}
