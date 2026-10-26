"""Orchestration and retained artifacts for the verified-input solution lane."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .lanes import (
    ImplicitSolutionLaneExecutor,
    LaneExecutionRequest,
    execute_lane,
)
from .solution_kinetics import ImplicitSolutionEnvironment, VerifiedSolutionNetwork


def _canonical_hash(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def run_implicit_solution_stability_bundle(
    bundle: Mapping[str, Any],
    workdir: Path,
    *,
    profile_points: int = 160,
) -> dict[str, Any]:
    """Validate and execute one restricted implicit-solution evidence bundle."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    retained_input = workdir / "solution-input-bundle.json"
    result_path = workdir / "solution-stability.json"
    try:
        if not isinstance(bundle, Mapping):
            raise TypeError("Solution stability bundle must be a mapping")
        if bundle.get("kind") != "verified_implicit_solution_stability_bundle":
            raise ValueError("Bundle kind must be 'verified_implicit_solution_stability_bundle'")
        serializable = dict(bundle)
        bundle_id = _canonical_hash(serializable)
        retained_input.write_text(json.dumps(serializable, indent=2, sort_keys=True) + "\n")
        environment = ImplicitSolutionEnvironment.from_dict(bundle.get("environment") or {})
        network = VerifiedSolutionNetwork.from_dict(bundle.get("network") or {})
        from .reaction_candidates import canonical_smiles

        declared_target = str(bundle.get("target_smiles") or "").strip()
        if not declared_target:
            raise ValueError("Solution stability bundle requires target_smiles")
        network_target = next(
            item for item in network.species if item.species_id == network.target_species_id
        )
        if canonical_smiles(declared_target) != canonical_smiles(network_target.smiles):
            raise ValueError("Bundle target_smiles does not match the network target species")
        environment_spec = environment.as_environment_spec()
        request = LaneExecutionRequest(
            environment=environment_spec,
            requested_observable="t95",
            inputs={
                "solution_environment": environment,
                "verified_solution_network": network,
            },
            compute_budget={"profile_points": int(profile_points)},
        )
        lane_run = execute_lane(ImplicitSolutionLaneExecutor(), request, workdir / "lane")
        result = {
            "schema_version": 1,
            "kind": "implicit_solution_stability_result",
            "status": lane_run["status"],
            "input_bundle_sha256": bundle_id,
            "lane_run": lane_run,
            "estimated_time_to_retention_seconds": (
                lane_run.get("lane_result", {}).get("observable_value")
                if lane_run.get("lane_result", {}).get("observable_status") == "verified"
                else None
            ),
            "claim_scope": (
                "Dark neutral homogeneous implicit solution with externally verified first-order rates "
                "and a supplied restricted network-coverage claim."
            ),
            "result_json": str(result_path),
        }
    except Exception as error:
        result = {
            "schema_version": 1,
            "kind": "implicit_solution_stability_result",
            "status": "invalid_or_incomplete_evidence",
            "reason": f"{type(error).__name__}: {error}",
            "estimated_time_to_retention_seconds": None,
            "claim_scope": "No solution observable claim; the evidence bundle failed validation.",
            "result_json": str(result_path),
        }
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def load_and_run_implicit_solution_stability(
    bundle_path: Path,
    workdir: Path,
    *,
    profile_points: int = 160,
) -> dict[str, Any]:
    payload = json.loads(Path(bundle_path).read_text())
    return run_implicit_solution_stability_bundle(payload, workdir, profile_points=profile_points)


__all__ = [
    "load_and_run_implicit_solution_stability", "run_implicit_solution_stability_bundle",
]
