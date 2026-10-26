"""Interrupt only an exact live P1 driver after a validated NEB checkpoint.

This produces controlled execution evidence, not a path/rate/resume-success claim.
"""

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import signal
import shlex
import subprocess
import time

from run_p1_preflight import THREAD_LIMITS, save
import run_p1_workflow as driver
from run_p1_workflow import artifact


def _absent(identifier, *, group=False):
    try:
        os.kill(-identifier if group else identifier, 0)
    except ProcessLookupError:
        return True
    except (PermissionError, OSError):
        return False
    return False


def _terminal_confirmation(attempt, neb_input, driver_pid, orca_pid, group_id, requested_at):
    """Require acknowledgment AND terminal process ownership, not a sent signal."""
    attempt, neb_input = Path(attempt).resolve(), Path(neb_input).resolve()
    result_path, lease_path = attempt / "result.json", neb_input.with_suffix(".process.json")
    if not result_path.is_file() or not lease_path.is_file():
        return None
    try:
        result, lease = json.loads(result_path.read_text()), json.loads(lease_path.read_text())
        ack = result.get("interruption_control") or {}
        received, ended = ack.get("received_at_epoch"), lease.get("ended_at_epoch")
        if (result.get("status") != "failed" or ack.get("received_signal") != "SIGINT"
                or not str(result.get("failure_reason", "")).startswith("KeyboardInterrupt:")
                or isinstance(received, bool) or not isinstance(received, (int, float))
                or not math.isfinite(received) or received < requested_at
                or isinstance(ended, bool) or not isinstance(ended, (int, float))
                or not math.isfinite(ended) or ended < received
                or lease.get("status") != "terminated" or type(lease.get("returncode")) is not int
                or lease.get("pid") != orca_pid or lease.get("process_group_id") != group_id
                or Path(lease.get("input", "")).resolve() != neb_input
                or lease.get("input_sha256") != artifact(neb_input)["sha256"]
                or result.get("signal_handler_evidence") != artifact(attempt / "driver-signal-handlers.json")
                or not all((_absent(driver_pid), _absent(orca_pid), _absent(group_id, group=True)))):
            return None
        return {"status": "controlled_interruption_confirmed", "result": artifact(result_path),
                "terminal_lease": artifact(lease_path), "received_at_epoch": ack["received_at_epoch"],
                "recorded_driver_pid_absent": True, "recorded_orca_pid_absent": True,
                "recorded_orca_group_absent": True, "path_or_rate_claim_allowed": False,
                "resume_equivalence_demonstrated": False}
    except (OSError, ValueError, TypeError):
        return None


def run(attempt, neb_input, output):
    os.environ.update(THREAD_LIMITS)
    from storca.reaction_path import _NEB_LBFGS_ITERATION, _validate_orca_neb_restart_band
    attempt, neb_input, output = Path(attempt).resolve(strict=True), Path(neb_input).resolve(strict=True), Path(output).resolve()
    if (attempt / "result.json").exists():
        raise ValueError("Attempt is already terminal; do not interrupt it")
    contract = json.loads((attempt / "execution-contract.json").read_text())
    if contract.get("budget", {}).get("maximum_cores") != 4:
        raise ValueError("Expected the frozen four-core P1 attempt")
    handler_path = attempt / "driver-signal-handlers.json"
    if not handler_path.is_file():
        raise ValueError("Installed SIGINT handler evidence is missing")
    handlers = json.loads(handler_path.read_text())
    if (handlers.get("handler_policy_id") != driver.SIGINT_HANDLER_POLICY_ID
            or contract.get("signal_handling_policy_id") != driver.SIGINT_HANDLER_POLICY_ID
            or handlers.get("driver_sha256") != contract.get("driver_sha256")
            or contract.get("driver_sha256") != artifact(Path(driver.__file__))["sha256"]
            or Path(handlers.get("attempt_dir", "")).resolve() != attempt
            or (handlers.get("sigint") or {}).get("installed_callable") is not True
            or not str((handlers.get("sigint") or {}).get("installed", "")).endswith("run.<locals>.controlled_sigint")):
        raise ValueError("Installed SIGINT handler evidence differs from the exact driver contract")
    grace = contract.get("budget", {}).get("cleanup_grace_seconds")
    if isinstance(grace, bool) or not isinstance(grace, (int, float)) or not math.isfinite(grace) or not 0 < grace <= 20:
        raise ValueError("A bounded cleanup confirmation interval is required")
    lease = json.loads(neb_input.with_suffix(".process.json").read_text())
    if (lease.get("status") != "running" or Path(lease["input"]).resolve() != neb_input
            or lease.get("input_sha256") != artifact(neb_input)["sha256"]):
        raise ValueError("The exact NEB process is not authoritatively running")
    orca_pid = int(lease["pid"])
    process = subprocess.check_output(["ps", "-p", str(orca_pid), "-o", "ppid=,command="], text=True).strip()
    driver_pid = int(process.split(maxsplit=1)[0])
    orca_command = shlex.split(process.split(maxsplit=1)[1])
    if len(orca_command) != 2 or Path(orca_command[0]).name != "orca" or orca_command[1] != neb_input.name:
        raise ValueError("NEB lease process identity differs from the live command")
    command = subprocess.check_output(["ps", "-p", str(driver_pid), "-o", "command="], text=True).strip()
    args = shlex.split(command)
    output_positions = [i for i, arg in enumerate(args) if arg == "--output-dir"]
    if (not any(Path(arg).resolve() == Path(driver.__file__).resolve() for arg in args)
            or len(output_positions) != 1 or output_positions[0] + 1 >= len(args)
            or Path(args[output_positions[0] + 1]).resolve() != attempt
            or handlers.get("driver_pid") != driver_pid):
        raise ValueError("NEB parent is not the exact selected P1 driver")
    group_id = lease.get("process_group_id")
    try:
        live_group = os.getpgid(orca_pid)
    except (ProcessLookupError, PermissionError) as error:
        raise ValueError("The selected NEB process is not provably live") from error
    if type(group_id) is not int or group_id <= 0 or live_group != group_id:
        raise ValueError("The recorded NEB process group differs from its live owner")
    band = neb_input.with_name(f"{neb_input.stem}_MEP.allxyz")
    check = _validate_orca_neb_restart_band(band, expected_images=10)
    if not check.get("valid"):
        raise ValueError(f"Checkpoint is not restartable: {check}")
    neb_output = neb_input.with_suffix(".out")
    try:
        output_bytes = neb_output.read_bytes()
    except OSError as error:
        raise ValueError("No completed NEB iteration output is available for checkpoint continuation") from error
    completed_rows = [line for line in output_bytes.decode(errors="replace").splitlines()
                      if _NEB_LBFGS_ITERATION.search(line)]
    if not completed_rows:
        raise ValueError("No completed NEB iteration exists; a band alone cannot justify interruption/resume")
    output.mkdir(parents=True, exist_ok=False)
    retained_output = output / "pre-interruption-output.out"
    with retained_output.open("xb") as stream:
        stream.write(output_bytes)
    retained_band = output / "pre-interruption-band.allxyz"
    shutil.copy2(band, retained_band)
    retained_check = _validate_orca_neb_restart_band(retained_band, expected_images=10)
    if not retained_check.get("valid") or artifact(band)["sha256"] != artifact(retained_band)["sha256"]:
        raise ValueError("The checkpoint changed while retaining interruption evidence")
    requested_at = time.time()
    save(output, "control.json", {
        "kind": "p1_controlled_inflight_interruption", "status": "signal_requested_not_terminal_proof",
        "attempt_contract": artifact(attempt / "execution-contract.json"),
        "source_input": artifact(neb_input), "source_lease": artifact(neb_input.with_suffix(".process.json")),
        "source_output": {"path": str(neb_output), "sha256": driver.digest(output_bytes)},
        "retained_output": artifact(retained_output),
        "completed_neb_iteration_count": len(completed_rows),
        "last_completed_neb_iteration_row": completed_rows[-1].strip(),
        "retained_band": artifact(retained_band), "band_integrity": check,
        "driver_pid": driver_pid, "driver_command": command, "orca_pid": orca_pid,
        "process_group_id": group_id, "handler_evidence": artifact(handler_path),
        "requested_at_epoch": requested_at, "confirmation_timeout_seconds": grace,
        "signal": "SIGINT_to_exact_driver", "path_or_rate_claim_allowed": False,
        "resume_equivalence_demonstrated": False,
    })
    if (attempt / "result.json").exists() or json.loads(neb_input.with_suffix(".process.json").read_text()).get("status") != "running":
        raise ValueError("The selected execution became terminal before the signal request")
    os.kill(driver_pid, signal.SIGINT)
    deadline, confirmation = time.monotonic() + grace, None
    while time.monotonic() < deadline:
        confirmation = _terminal_confirmation(attempt, neb_input, driver_pid, orca_pid, group_id, requested_at)
        if confirmation is not None:
            break
        time.sleep(min(.05, max(0, deadline - time.monotonic())))
    if confirmation is None:
        confirmation = {"status": "signal_request_not_confirmed_within_cleanup_bound",
                        "path_or_rate_claim_allowed": False, "resume_equivalence_demonstrated": False}
    save(output, "confirmation.json", confirmation)
    return {"driver_pid": driver_pid, "control_artifact": str(output / "control.json"),
            "confirmation": str(output / "confirmation.json"), "status": confirmation["status"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempt", type=Path, required=True)
    parser.add_argument("--neb-input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.attempt, args.neb_input, args.output_dir)))
