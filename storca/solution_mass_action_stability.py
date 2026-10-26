"""Retained-bundle orchestration for the solution mass-action lane."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .lanes.execution import LaneExecutionRequest, execute_lane
from .lanes.solution_mass_action import ImplicitSolutionMassActionLaneExecutor
from .solution_mass_action import MassActionSolutionEnvironment, VerifiedMassActionNetwork


def _canonical_hash(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def run_implicit_solution_mass_action_bundle(
    bundle: Mapping[str, Any],
    workdir: Path,
    *,
    profile_points: int = 180,
    relative_tolerance: float = 1.0e-7,
    absolute_tolerance: float = 1.0e-12,
    maximum_steps: int = 200_000,
    requested_observable: str = "t95",
    observable_parameters: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    input_path = workdir / "solution-mass-action-input-bundle.json"
    result_path = workdir / "solution-mass-action-stability.json"
    try:
        if not isinstance(bundle, Mapping):
            raise TypeError("Solution mass-action stability bundle must be a mapping")
        if bundle.get("kind") != "verified_implicit_solution_mass_action_stability_bundle":
            raise ValueError(
                "Bundle kind must be 'verified_implicit_solution_mass_action_stability_bundle'"
            )
        serializable = dict(bundle)
        bundle_id = _canonical_hash(serializable)
        input_path.write_text(json.dumps(serializable, indent=2, sort_keys=True) + "\n")
        environment = MassActionSolutionEnvironment.from_dict(bundle.get("environment") or {})
        network = VerifiedMassActionNetwork.from_dict(bundle.get("network") or {})
        from .reaction_candidates import canonical_smiles

        declared_target = str(bundle.get("target_smiles") or "").strip()
        if not declared_target:
            raise ValueError("Solution mass-action bundle requires target_smiles")
        target = next(item for item in network.species if item.species_id == network.target_species_id)
        if canonical_smiles(declared_target) != canonical_smiles(target.smiles):
            raise ValueError("Bundle target_smiles does not match the network target species")
        lane_run = execute_lane(
            ImplicitSolutionMassActionLaneExecutor(),
            LaneExecutionRequest(
                environment=environment.as_environment_spec(),
                requested_observable=str(requested_observable),
                inputs={
                    "solution_environment": environment,
                    "verified_solution_mass_action_network": network,
                    "observable_parameters": dict(observable_parameters or {}),
                },
                compute_budget={
                    "profile_points": int(profile_points),
                    "relative_tolerance": float(relative_tolerance),
                    "absolute_tolerance": float(absolute_tolerance),
                    "maximum_steps": int(maximum_steps),
                },
            ),
            workdir / "lane",
        )
        result = {
            "schema_version": 1,
            "kind": "implicit_solution_mass_action_stability_result",
            "status": lane_run["status"],
            "input_bundle_sha256": bundle_id,
            "lane_run": lane_run,
            "estimated_time_to_retention_seconds": (
                lane_run.get("lane_result", {}).get("observable_value")
                if (
                    requested_observable == "t95"
                    and lane_run.get("lane_result", {}).get("observable_status") == "verified"
                ) else None
            ),
            "claim_scope": (
                "Dark neutral ideal-activity homogeneous implicit solution with verified dimensional "
                "first/second-order rates, diffusion ceilings, and supplied restricted network coverage."
            ),
            "result_json": str(result_path),
        }
    except Exception as error:
        result = {
            "schema_version": 1,
            "kind": "implicit_solution_mass_action_stability_result",
            "status": "invalid_or_incomplete_evidence",
            "reason": f"{type(error).__name__}: {error}",
            "estimated_time_to_retention_seconds": None,
            "claim_scope": "No solution observable claim; the mass-action evidence bundle failed validation.",
            "result_json": str(result_path),
        }
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def load_and_run_implicit_solution_mass_action(
    bundle_path: Path, workdir: Path, **kwargs: Any,
) -> dict[str, Any]:
    return run_implicit_solution_mass_action_bundle(
        json.loads(Path(bundle_path).read_text()), workdir, **kwargs,
    )


__all__ = [
    "load_and_run_implicit_solution_mass_action", "run_implicit_solution_mass_action_bundle",
]
