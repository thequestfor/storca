"""Run one reference-free OH/methane development workflow through STORCA.

This is a selected-rate development attempt, not an experimental flow-reactor
reconstruction, methane lifetime prediction, or benchmark acceptance decision.
Use a fresh directory. No product, atom mapping, route, or reference rate is
provided to production discovery. The independently retained scoring evidence
is deliberately not read by this driver.
"""

from __future__ import annotations

import argparse
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time
import traceback


THREAD_LIMITS = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
}
POLICY_ID = "oh-methane-reference-free-production-driver-v1"
BUDGET = {
    "overall_timeout_seconds": 86400.0,
    "maximum_cores": 4,
    "verification_timeout_seconds_per_job": 14400.0,
    "verification_max_iterations": 3,
    "verification_orientations": 3,
    "verification_neb_images": 8,
    "rmg_walltime": "00:00:30:00",
    "rmg_hard_timeout_seconds": 2100.0,
    "rmg_heartbeat_timeout_seconds": 300.0,
    "rmg_max_processes": 1,
    "rmg_max_iterations": 60,
    "rmg_max_edge_species": 300,
    "rmg_max_objects_per_iteration": 2,
    "rmg_maximum_heavy_atoms": 3,
    "rmg_retry_attempts": 2,
}
TEMPERATURE_K = 340.3
PRESSURE_BAR = 1.0
INITIAL_OH_NUMBER_DENSITY_CM3 = 1.0e9
BOLTZMANN_J_K = 1.380649e-23


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def save(output, name, value):
    path = Path(output) / name
    with path.open("x") as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return path


def artifact(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": digest_bytes(path.read_bytes())}


def scenario_configuration():
    # P = n k T; the arbitrary initial OH concentration is used solely to
    # obtain a finite, dilute discovery trajectory, not as a measured value.
    hydroxyl = INITIAL_OH_NUMBER_DENSITY_CM3 * 1.0e6 * BOLTZMANN_J_K * TEMPERATURE_K / (PRESSURE_BAR * 1.0e5)
    dry_remainder = 1.0 - hydroxyl - 63.0e-6 - 0.013
    return {
        "name": "oh-methane-rate-only-development-proxy",
        "phase": "homogeneous gas-phase surrogate",
        "atmosphere": "explicit humid synthetic-air proxy",
        "model_applicability": (
            "A dark, closed homogeneous batch initialized with OH, methane and a declared humid-air proxy; "
            "not the experimental illuminated flow tube. Initial OH density, dry-air ratio, absence of ozone "
            "and the 0.5-second discovery window are modeling assumptions, not measured conditions. "
            "No ongoing photochemical or radical source, experimental trajectory or methane lifetime is modeled. "
            "A generated gas TST coefficient requires a separate finite-pressure, tunneling, electronic-state "
            "and statistical-counting applicability assessment before experimental rate comparison."
        ),
        "additional_species": [
            {"label": "methane", "smiles": "C", "charge": 0, "multiplicity": 1, "reactive": True},
            {"label": "water", "smiles": "O", "charge": 0, "multiplicity": 1, "reactive": True},
            {"label": "oxygen", "smiles": "[O][O]", "charge": 0, "multiplicity": 3, "reactive": True},
            {"label": "nitrogen", "smiles": "N#N", "charge": 0, "multiplicity": 1, "reactive": False},
        ],
        "initial_mole_fractions": {
            "stability": hydroxyl, "methane": 63.0e-6, "water": 0.013,
            "oxygen": dry_remainder * 0.2095, "nitrogen": dry_remainder * 0.7905,
        },
    }


class _Tee:
    def __init__(self, terminal, logfile):
        self.terminal, self.logfile = terminal, logfile

    def write(self, value):
        self.terminal.write(value)
        self.logfile.write(value)
        self.flush()
        return len(value)

    def flush(self):
        self.terminal.flush()
        self.logfile.flush()

    def isatty(self):
        return False


def run(output, *, rmg_env):
    """Dispatch exactly one production workflow and preserve every outcome."""
    if not isinstance(rmg_env, str) or not rmg_env.strip():
        raise ValueError("Specify the existing RMG/Arkane conda environment with --rmg-env")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    os.environ.update(THREAD_LIMITS)
    os.environ["STORCA_ORCA_LIVE_OUTPUT"] = "1"
    from storca.conditions import build_condition_spec
    from storca.methods import harmonic_method_profile
    from storca.stability import run_stability_screen

    scenario = scenario_configuration()
    profile = harmonic_method_profile()
    conditions = build_condition_spec(
        scenario, temperature_K=TEMPERATURE_K, pressure_bar=PRESSURE_BAR,
        target_duration_hours=0.5 / 3600.0, retention_fraction=0.95, light_condition="dark",
    ).as_dict()
    request = {
        "schema_version": 1, "kind": "storca_oh_methane_development_execution_request",
        "policy_id": POLICY_ID,
        "chemical_system": {"tracked_parent": {"smiles": "[OH]", "charge": 0, "multiplicity": 2}},
        "scenario_config": scenario, "condition_contract": conditions,
        "requested_observables": ["retained_discovered_bimolecular_rate_evidence"],
        "modeling_assumptions": {
            "initial_OH_number_density_cm^-3": INITIAL_OH_NUMBER_DENSITY_CM3,
            "initial_OH_density_measured": False,
            "dry_air_oxygen_fraction": 0.2095, "dry_air_ratio_measured": False,
            "reactor": "closed_batch_proxy_not_experimental_flow_tube",
            "ongoing_radical_generation": False, "ozone_included": False,
            "discovery_window_seconds": 0.5, "window_measured": False,
        },
        "reference_answers_supplied": False, "reference_sources_read": False,
        "benchmark_acceptance_thresholds": None, "support_claim_allowed": False,
        "lifetime_claim_allowed": False, "calibration_claim_allowed": False,
    }
    request_path = save(output, "request.json", request)
    workflow = output / "workflow"
    workflow.mkdir(exist_ok=False)
    contract = {
        "schema_version": 1, "kind": "storca_oh_methane_development_execution_contract",
        "policy_id": POLICY_ID, "request": artifact(request_path),
        "driver": artifact(Path(__file__)), "condition_contract": conditions,
        "method_profile": profile, "rmg_environment": rmg_env,
        "budget": dict(BUDGET), "thread_limits": dict(THREAD_LIMITS),
        "driver_pid": os.getpid(), "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "production_api": "storca.stability.run_stability_screen",
        "graph_xtb_discovery": False, "auto_reference_libraries": False,
        "rate_uncertainty_evidence_supplied": False,
        "retry_scope": "Only declared existing bounded RMG retry; no wrapper retry or additional cases.",
        "interruption_policy": "First SIGINT or overall SIGALRM raises through existing owned-job cleanup; repeated requests do not interrupt cleanup.",
        "scientific_acceptance_pass": None, "physical_prediction_claim_allowed": False,
    }
    contract_path = save(output, "execution-contract.json", contract)
    started = time.monotonic()
    interruption = {"received_signal": None}
    outcome = {
        "schema_version": 1, "kind": "storca_oh_methane_development_outcome",
        "status": "failed", "execution_contract": artifact(contract_path),
        "scientific_acceptance_pass": None, "benchmark_scoring_performed": False,
        "support_claim_allowed": False, "lifetime_claim_allowed": False,
        "calibration_claim_allowed": False,
    }

    def stop(signum, _frame):
        if interruption["received_signal"] is not None:
            return
        interruption.update(received_signal=signal.Signals(signum).name,
                            received_at_epoch=time.time())
        raise KeyboardInterrupt("OH/methane development run interrupted by " + signal.Signals(signum).name)

    previous_sigint = signal.signal(signal.SIGINT, stop)
    previous_sigalrm = signal.signal(signal.SIGALRM, stop)
    previous_timer = signal.setitimer(signal.ITIMER_REAL, BUDGET["overall_timeout_seconds"])
    try:
        with (output / "console.log").open("x", buffering=1) as logfile:
            with redirect_stdout(_Tee(sys.stdout, logfile)), redirect_stderr(_Tee(sys.stderr, logfile)):
                print(f"[STORCA benchmark] Driver PID {os.getpid()}; output {output}", flush=True)
                print("[STORCA benchmark] Reference-free OH/methane discovery; 4 ORCA cores, one RMG worker.", flush=True)
                print("[STORCA benchmark] 24-hour overall ceiling; 4-hour per-route job ceiling. No conclusive lifetime is requested.", flush=True)
                try:
                    result = run_stability_screen(
                        "[OH]", workflow, rmg_env=rmg_env,
                        scenario_config=scenario, temperature=TEMPERATURE_K, pressure=PRESSURE_BAR,
                        target_duration_hours=0.5 / 3600.0, retention_fraction=0.95,
                        light_condition="dark", charge=0, multiplicity=2,
                        ncores=BUDGET["maximum_cores"], method_profile=profile,
                        auto_reference_libraries=False, auto_verify_routes=True,
                        verification_max_iterations=BUDGET["verification_max_iterations"],
                        verification_timeout_seconds=BUDGET["verification_timeout_seconds_per_job"],
                        verification_orientations=BUDGET["verification_orientations"],
                        verification_neb_images=BUDGET["verification_neb_images"],
                        graph_xtb_discovery=False,
                        **{key: value for key, value in BUDGET.items() if key.startswith("rmg_")},
                    )
                    result_path = result.get("result_json")
                    outcome.update(
                        status="workflow_returned", result_json=str(result_path) if result_path else None,
                        assessment=result.get("assessment"),
                        verification_status=(result.get("verification_summary") or {}).get("status"),
                    )
                    if result_path and Path(result_path).is_file():
                        outcome["production_result"] = artifact(result_path)
                    print("[STORCA benchmark] Production workflow returned. This is not a chemical-accuracy pass.", flush=True)
                except BaseException as error:
                    outcome.update(
                        status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                        failure_reason=f"{type(error).__name__}: {error}",
                        traceback=traceback.format_exc(),
                        stdout=str(getattr(error, "stdout", "") or ""),
                        stderr=str(getattr(error, "stderr", "") or ""),
                    )
                    print(outcome["traceback"], flush=True)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_sigalrm)
        signal.signal(signal.SIGINT, previous_sigint)
        elapsed = time.monotonic() - started
        if previous_timer[0]:
            signal.setitimer(signal.ITIMER_REAL, max(0.001, previous_timer[0] - elapsed), previous_timer[1])
        outcome.update(elapsed_seconds=elapsed, interruption_control=interruption,
                       request_unchanged=artifact(request_path) == contract["request"],
                       driver_unchanged=artifact(Path(__file__)) == contract["driver"],
                       driver_pid=os.getpid(),
                       completed_at_utc=datetime.now(timezone.utc).isoformat())
        if (output / "console.log").is_file():
            outcome["console_log"] = artifact(output / "console.log")
        result_path = save(output, "result.json", outcome)
    print(f"[STORCA benchmark] Retained outcome: {result_path}", flush=True)
    return outcome


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New development-run directory")
    parser.add_argument("--rmg-env", required=True, help="Existing RMG/Arkane conda environment")
    args = parser.parse_args()
    outcome = run(args.output, rmg_env=args.rmg_env)
    return 0 if outcome["status"] == "workflow_returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
