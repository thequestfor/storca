"""Bounded production stability run for a preregistered development case.

No routes, products, rates or uncertainty calibration are injected. Parent
minimum evidence is passed through the production reuse validator. An ordinary
workflow return is not a scientific acceptance decision. Use a new directory.
"""

import argparse
import json
import math
import os
from pathlib import Path
import signal
import time
import traceback

from run_p1_preflight import THREAD_LIMITS, digest, manifest_conditions, save


BUDGET = {
    # A per-dispatch cap, never a renewal of the shared 2400-second campaign.
    "overall_timeout_seconds": 2340.0,
    "cleanup_grace_seconds": 20.0,
    "maximum_cores": 4,
    "route_timeout_seconds_per_job": 600.0,
    "verification_max_iterations": 3,
    "verification_orientations": 3,
    "verification_neb_images": 8,
    "rmg_retry_attempts": 1,
    "rmg_walltime": "00:00:02:00",
    "rmg_hard_timeout_seconds": 180.0,
    "rmg_heartbeat_timeout_seconds": 120.0,
    "rmg_max_processes": 1,
    "rmg_max_iterations": 20,
    "rmg_max_edge_species": 100,
    "rmg_max_objects_per_iteration": 1,
}

CAMPAIGN_CEILING_SECONDS = 2400.0
SIGINT_HANDLER_POLICY_ID = "p1-driver-explicit-sigint-v1"


def signal_handler_label(handler):
    """Serialize the disposition actually returned by ``signal.getsignal``."""
    if handler == signal.SIG_IGN:
        return "SIG_IGN"
    if handler == signal.SIG_DFL:
        return "SIG_DFL"
    return f"{getattr(handler, '__module__', '')}.{getattr(handler, '__qualname__', repr(handler))}"


def production_budget(history_root, manifest_sha256, *, exclude=None):
    """Charge every terminal production attempt, including failed/resumed ones.

    A missing result is not evidence of termination; do not start a duplicate.
    This diagnostic driver's ledger is distinct from candidate preflight costs.
    """
    records = []
    for directory in sorted(Path(history_root).glob("p1-cyclobutene-workflow-*")):
        if exclude is not None and directory.resolve() == Path(exclude).resolve():
            continue
        contract_path = directory / "execution-contract.json"
        if not contract_path.is_file():
            raise ValueError(f"Production history lacks a dispatch contract: {directory}")
        contract = json.loads(contract_path.read_text())
        if contract.get("manifest_sha256") != manifest_sha256:
            raise ValueError(f"Production history has a different preregistration: {directory}")
        result_path = directory / "result.json"
        if not result_path.is_file():
            raise ValueError(f"Production history is not authoritatively terminal: {directory}")
        result = json.loads(result_path.read_text())
        elapsed = result.get("elapsed_seconds")
        if isinstance(elapsed, bool) or not isinstance(elapsed, (int, float)) \
                or not math.isfinite(elapsed) or elapsed < 0:
            raise ValueError(f"Production history has invalid measured cost: {directory}")
        records.append({"contract": artifact(contract_path), "result": artifact(result_path),
                        "elapsed_seconds": elapsed})
    consumed = math.fsum(record["elapsed_seconds"] for record in records)
    remaining = max(0.0, CAMPAIGN_CEILING_SECONDS - consumed)
    if remaining <= 1.0:
        raise ValueError("Shared P1 production budget exhausted; additional authority is required")
    budget = {**BUDGET, "overall_timeout_seconds": min(BUDGET["overall_timeout_seconds"], remaining)}
    return budget, {"kind": "p1_shared_production_budget", "ceiling_seconds": CAMPAIGN_CEILING_SECONDS,
                    "consumed_seconds": consumed, "remaining_seconds": remaining, "attempts": records}


def artifact(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": digest(path.read_bytes())}


def run(manifest_path, parent_path, output, *, rmg_env, resume_run=None):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    raw = Path(manifest_path).read_bytes()
    smiles, scenario, conditions = manifest_conditions(json.loads(raw))
    parent_record = artifact(parent_path)
    parent = json.loads(Path(parent_path).read_text())
    parent_sources = {role: artifact(path) for role, path in (parent.get("artifacts") or {}).items()}
    budget, accounting = production_budget(output.parent, digest(raw), exclude=output)
    resume_source = None
    if resume_run is not None:
        resume_run = Path(resume_run).resolve(strict=True)
        source = resume_run / "verification-execution-contract.json"
        retained = json.loads(source.read_text())
        if (retained.get("condition_contract") != conditions.as_dict()
                or retained.get("method_profile") != parent.get("method_profile")
                or retained.get("parent") != {"smiles": smiles, "charge": 0, "multiplicity": 1}
                or (retained.get("rmg") or {}).get("environment") != rmg_env
                or (retained.get("resources") or {}).get("ncores") != budget["maximum_cores"]):
            raise ValueError("Retained workflow differs from the frozen P1 dispatch")
        resume_source = artifact(source)
    os.environ.update(THREAD_LIMITS)
    os.environ["STORCA_ORCA_LIVE_OUTPUT"] = "1"
    from storca.stability import run_stability_screen

    contract = {
        "kind": "p1_bounded_production_workflow", "manifest_sha256": digest(raw),
        "driver_sha256": digest(Path(__file__).read_bytes()),
        "condition_contract": conditions.as_dict(), "method_profile": parent.get("method_profile"),
        "parent_evidence": parent_record, "parent_source_artifacts": parent_sources,
        "rmg_environment": rmg_env, "budget": budget, "shared_budget_accounting": accounting,
        "resume_execution_contract": resume_source, "resume_run": str(resume_run) if resume_run else None,
        "thread_limits": THREAD_LIMITS,
        "signal_handling_policy_id": SIGINT_HANDLER_POLICY_ID,
        "graph_xtb_discovery": False, "auto_reference_libraries": False,
        "rate_uncertainty_evidence_supplied": False,
        "interpretation": (
            "Production path/rate/repair exercise. Missing independent coverage and rate "
            "uncertainty must remain unresolved; a normal return is not a P1 scientific pass."
        ),
    }
    save(output, "execution-contract.json", contract)
    with (output / "preregistration.json").open("xb") as stream:
        stream.write(raw)
    workflow = resume_run if resume_run else output / "workflow"
    if resume_run is None:
        workflow.mkdir(exist_ok=False)
    started = time.monotonic()
    outcome = {"status": "failed", "budget": budget, "p1_scientific_gate_achieved": False}
    interruption = {"received_signal": None}

    def controlled_sigint(signum, _frame):
        # Explicitly replace inherited SIG_IGN. Raising into the existing
        # ORCA/RMG runners invokes their owned-process-group cleanup. A second
        # request must not interrupt that cleanup while it is already running.
        if interruption["received_signal"] is not None:
            return
        interruption.update(received_signal=signal.Signals(signum).name,
                            received_at_epoch=time.time())
        raise KeyboardInterrupt("P1 controlled SIGINT interruption requested")

    def deadline(_signum, _frame):
        raise KeyboardInterrupt("P1 production workflow exhausted its recorded overall budget")

    previous_sigint = signal.signal(signal.SIGINT, controlled_sigint)
    previous = signal.signal(signal.SIGALRM, deadline)
    previous_timer = signal.setitimer(signal.ITIMER_REAL, budget["overall_timeout_seconds"])
    try:
        save(output, "driver-signal-handlers.json", {
            "kind": "p1_driver_signal_handlers", "handler_policy_id": SIGINT_HANDLER_POLICY_ID,
            "driver_pid": os.getpid(), "attempt_dir": str(output),
            "driver_sha256": contract["driver_sha256"],
            "sigint": {"previous": signal_handler_label(previous_sigint),
                       "installed": signal_handler_label(signal.getsignal(signal.SIGINT)),
                       "installed_callable": callable(signal.getsignal(signal.SIGINT))},
            "sigalrm": {"installed": signal_handler_label(signal.getsignal(signal.SIGALRM))},
            "cleanup": "KeyboardInterrupt propagates through existing ORCA/RMG owned-process-group cleanup",
        })
        if resume_run is not None:
            from storca.verification_resume import resume_stability_verification
            result = resume_stability_verification(resume_run)
        else:
            rmg_resources = {key: value for key, value in budget.items() if key.startswith("rmg_")}
            result = run_stability_screen(
                smiles, workflow, rmg_env=rmg_env,
                scenario_config=scenario, temperature=conditions.temperature_K,
                pressure=conditions.pressure_bar, target_duration_hours=conditions.target_duration_seconds / 3600,
                retention_fraction=conditions.retention_fraction, light_condition="dark",
                charge=0, multiplicity=1, ncores=budget["maximum_cores"],
                method_profile=parent["method_profile"], precomputed_orca_evidence=parent,
                auto_reference_libraries=False, auto_verify_routes=True,
                verification_max_iterations=budget["verification_max_iterations"],
                verification_timeout_seconds=budget["route_timeout_seconds_per_job"],
                verification_orientations=budget["verification_orientations"],
                verification_neb_images=budget["verification_neb_images"],
                graph_xtb_discovery=False, **rmg_resources,
            )
        outcome.update(status="workflow_returned", result_json=str(result["result_json"]),
                       assessment=result.get("assessment"),
                       verification_status=(result.get("verification_summary") or {}).get("status"),
                       observable=(result.get("kinetic_lifetime") or {}).get("estimated_time_to_retention_seconds"))
    except BaseException as error:
        outcome.update(failure_reason=f"{type(error).__name__}: {error}",
                       traceback=traceback.format_exc(),
                       stdout=getattr(error, "stdout", ""), stderr=getattr(error, "stderr", ""))
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
        signal.signal(signal.SIGINT, previous_sigint)
        elapsed = time.monotonic() - started
        if previous_timer[0]:
            signal.setitimer(signal.ITIMER_REAL, max(0.001, previous_timer[0] - elapsed), previous_timer[1])
        outcome["elapsed_seconds"] = elapsed
        outcome["interruption_control"] = interruption
        outcome["signal_handler_evidence"] = artifact(output / "driver-signal-handlers.json")
        outcome["manifest_unchanged"] = digest(Path(manifest_path).read_bytes()) == digest(raw)
        outcome["parent_evidence_unchanged"] = artifact(parent_path) == parent_record
        outcome["parent_artifacts_unchanged"] = all(artifact(record["path"]) == record
                                                    for record in parent_sources.values())
        save(output, "result.json", outcome)
    print(json.dumps(outcome, sort_keys=True))
    return outcome


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--parent-evidence", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rmg-env", required=True)
    parser.add_argument("--resume-run", type=Path, help="Continue a retained production workflow, without new discovery")
    args = parser.parse_args()
    result = run(args.manifest, args.parent_evidence, args.output_dir, rmg_env=args.rmg_env,
                 resume_run=args.resume_run)
    raise SystemExit(0 if result["status"] == "workflow_returned" else 1)
