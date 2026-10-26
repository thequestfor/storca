"""Run one preregistered, bounded production RMG discovery preflight.

This driver does not run ORCA or Arkane, inject a pathway, score experimental
agreement, or establish a stability prediction. Use a new output directory for
every invocation. The registered expected product is never passed to discovery.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback


REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

BUDGET = {
    "overall_timeout_seconds": 600.0,
    "rmg_retry_attempts": 1,
    "rmg_walltime": "00:00:02:00",
    "rmg_hard_timeout_seconds": 180.0,
    "rmg_heartbeat_timeout_seconds": 120.0,
    "rmg_max_processes": 1,
    "rmg_max_iterations": 20,
    "rmg_max_edge_species": 100,
    "rmg_max_objects_per_iteration": 1,
}
THREAD_LIMITS = {
    "OMP_NUM_THREADS": "1",
    "OMP_THREAD_LIMIT": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def save(directory: Path, name: str, value: object) -> None:
    with (directory / name).open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def manifest_conditions(manifest: dict):
    """Translate only the explicitly registered pure neutral gas conditions."""
    from storca.conditions import build_condition_spec

    if (manifest.get("schema_version") != 1
            or manifest.get("kind") != "storca_p1_development_preregistration"
            or manifest.get("predictions_inspected_at_registration") is not False):
        raise ValueError("A pre-prediction P1 preregistration is required")
    system = manifest["chemical_system"]
    parent = system["parent"]
    declared = manifest["conditions"]
    if (system.get("pathway_injection_allowed") is not False
            or parent.get("charge") != 0 or parent.get("multiplicity") != 1
            or declared.get("phase") != "gas"
            or declared.get("added_bath_gas") is not None
            or declared.get("environment") != "homogeneous thermal gas; no irradiation"
            or declared.get("reactor") != "isothermal ideal-gas batch model"
            or declared.get("initial_composition") != [
                {"smiles": parent["smiles"], "mole_fraction": 1.0}
            ]):
        raise ValueError("This preflight requires the registered pure, dark, neutral-singlet gas model")
    for name in ("temperature_K", "pressure_bar", "duration_s", "retention_fraction"):
        value = declared[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"Registered {name} must be a finite number")
    scenario = {
        "name": "p1-pure-gas-development",
        "phase": "homogeneous gas-phase pure-target model",
        "atmosphere": "pure target; no added bath gas",
        "additional_species": [],
        "initial_mole_fractions": {"stability": 1.0},
        "model_applicability": (
            "Isothermal ideal-gas batch, pure parent initially, dark; no surface, "
            "condensed-phase, sampling-withdrawal or excited-state model. "
            + str(declared.get("model_vs_experiment_caveat") or "")
        ),
    }
    conditions = build_condition_spec(
        scenario, temperature_K=declared["temperature_K"], pressure_bar=declared["pressure_bar"],
        target_duration_hours=declared["duration_s"] / 3600.0,
        retention_fraction=declared["retention_fraction"], light_condition="dark",
    )
    return parent["smiles"], scenario, conditions


def worker(output: Path) -> int:
    """Collect production evidence in the externally bounded child process."""
    started = time.monotonic()
    try:
        contract = json.loads((output / "execution-contract.json").read_text())
        raw = (output / "preregistration.json").read_bytes()
        if digest(raw) != contract["manifest_sha256"]:
            raise ValueError("Retained preregistration hash mismatch")
        if digest(Path(__file__).read_bytes()) != contract["driver_sha256"]:
            raise ValueError("Preflight driver changed after dispatch")
        if contract["budget"] != BUDGET or contract["thread_limits"] != THREAD_LIMITS:
            raise ValueError("Frozen preflight resource contract mismatch")
        os.environ.update(THREAD_LIMITS)
        smiles, scenario, conditions = manifest_conditions(json.loads(raw))
        if conditions.as_dict() != contract["condition_contract"]:
            raise ValueError("Retained conditions do not match dispatch contract")
        # Import only after capping numerical-library threads, before chemistry.
        from storca.stability import collect_rmg_evidence

        arguments = {key: value for key, value in BUDGET.items() if key != "overall_timeout_seconds"}
        evidence = collect_rmg_evidence(
            smiles, output / "discovery", rmg_env=contract["rmg_environment"],
            temperature=conditions.temperature_K, pressure=conditions.pressure_bar,
            scenario=scenario, conditions=conditions, requested_phase="gas", **arguments,
        )
        save(output, "rmg-evidence.json", evidence)
        save(output, "worker-result.json", {
            "status": "completed",
            "discovery_status": evidence.get("status"),
            "enumeration_status": evidence.get("enumeration_status"),
            "candidate_route_count": len(evidence.get("candidate_routes") or []),
            "retained_partial_route_count": len(evidence.get("retained_partial_routes") or []),
            "failure_reason": evidence.get("failure_reason"),
            "elapsed_seconds": time.monotonic() - started,
            "scientific_prediction_claim_allowed": False,
        })
        return 0
    except Exception as error:
        save(output, "worker-failure.json", {
            "status": "failed", "failure_type": type(error).__name__, "failure_reason": str(error),
            "traceback": traceback.format_exc(), "elapsed_seconds": time.monotonic() - started,
            "scientific_prediction_claim_allowed": False,
        })
        traceback.print_exc()
        return 1


def run_preflight(manifest_path: Path, output: Path, rmg_env: str) -> dict:
    """Reserve an exclusive directory and supervise exactly one child run."""
    from src.orca_runner import run_bounded_kinetics_command

    manifest_path, output = manifest_path.resolve(), output.resolve()
    raw = manifest_path.read_bytes()
    manifest = json.loads(raw)
    _, scenario, conditions = manifest_conditions(manifest)
    if not rmg_env or not rmg_env.strip():
        raise ValueError("A named --rmg-env is required")
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    with (output / "preregistration.json").open("xb") as stream:
        stream.write(raw)
    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--output-dir", str(output)]
    save(output, "execution-contract.json", {
        "schema_version": 1, "kind": "storca_p1_discovery_preflight",
        "case_id": manifest["case_id"], "manifest_source": str(manifest_path),
        "manifest_sha256": digest(raw), "driver_sha256": digest(Path(__file__).read_bytes()),
        "command": command, "rmg_environment": rmg_env, "budget": BUDGET,
        "thread_limits": THREAD_LIMITS, "maximum_cores": 1,
        "condition_contract": conditions.as_dict(), "scenario": scenario,
        "expected_product_supplied_to_discovery": False, "external_pathways_supplied": False,
        "quantum_jobs_allowed": False, "scientific_prediction_claim_allowed": False,
        "budget_note": (
            "The self-child ceiling includes imports, RMG, parsing and bridge calls; no retries. "
            "RMG also has a 180 s external stop and at most 30 s shutdown grace. "
            "Overall process-group cleanup may take up to 5 additional seconds."
        ),
    })
    result = {
        "status": "failed", "case_id": manifest["case_id"],
        "scientific_prediction_claim_allowed": False,
        "manifest_sha256": digest(raw), "budget": BUDGET,
    }
    stdout, stderr = "", ""
    try:
        execution = run_bounded_kinetics_command(
            command, cwd=REPOSITORY, timeout_seconds=BUDGET["overall_timeout_seconds"],
        )
        stdout, stderr = execution.stdout, execution.stderr
        result["child_returncode"] = execution.returncode
        if execution.returncode:
            result["failure_reason"] = "preflight_child_failed"
        elif not (output / "worker-result.json").is_file() or not (output / "rmg-evidence.json").is_file():
            result["failure_reason"] = "preflight_child_missing_evidence"
        else:
            result["status"] = "completed"
            result["worker"] = json.loads((output / "worker-result.json").read_text())
            result["rmg_evidence_sha256"] = digest((output / "rmg-evidence.json").read_bytes())
    except subprocess.TimeoutExpired as error:
        stdout, stderr = error.output or "", error.stderr or ""
        result["failure_reason"] = "overall_preflight_timeout"
    except Exception as error:
        result.update(failure_reason=str(error), failure_type=type(error).__name__)
        stderr = traceback.format_exc()
    finally:
        for name, value in (("stdout.log", stdout), ("stderr.log", stderr)):
            with (output / name).open("x", encoding="utf-8") as stream:
                stream.write(value.decode(errors="replace") if isinstance(value, bytes) else value)
        result["elapsed_seconds"] = time.monotonic() - started
        try:
            result["manifest_source_unchanged"] = digest(manifest_path.read_bytes()) == digest(raw)
        except OSError:
            result["manifest_source_unchanged"] = False
        if not result["manifest_source_unchanged"]:
            result.update(status="failed", failure_reason="manifest_source_changed_during_preflight")
        save(output, "result.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rmg-env", default="rmg4_env")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    arguments = parser.parse_args()
    if arguments.worker:
        return worker(arguments.output_dir.resolve())
    if arguments.manifest is None:
        parser.error("--manifest is required")
    result = run_preflight(arguments.manifest, arguments.output_dir, arguments.rmg_env)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
