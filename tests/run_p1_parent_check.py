"""Bounded real parent-minimum check for the preregistered P1 case.

Uses the production collector; no pathway, product, rate, or experimental
answer is supplied. Completed evidence is eligible for the production reuse
validator. This stage alone cannot establish the P1 prediction gate.
"""

import argparse
import json
import os
from pathlib import Path
import signal
import time
import traceback

from run_p1_preflight import THREAD_LIMITS, digest, manifest_conditions, save


BUDGET = {"overall_timeout_seconds": 500.0, "orca_timeout_seconds_per_job": 180.0,
          "maximum_cores": 4, "cleanup_grace_seconds": 20.0}


def run(manifest_path, output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    raw = Path(manifest_path).read_bytes()
    manifest = json.loads(raw)
    smiles, _, conditions = manifest_conditions(manifest)
    os.environ.update(THREAD_LIMITS)
    # Streaming mode owns and cleans ORCA's separate process group on timeout
    # and interruption; do not use buffered subprocess execution here.
    os.environ["STORCA_ORCA_LIVE_OUTPUT"] = "1"
    from storca.methods import harmonic_method_profile
    from storca.stability import collect_orca_evidence

    profile = harmonic_method_profile()
    save(output, "execution-contract.json", {
        "kind": "p1_parent_minimum_check", "manifest_sha256": digest(raw),
        "driver_sha256": digest(Path(__file__).read_bytes()),
        "condition_contract": conditions.as_dict(), "method_profile": profile,
        "budget": BUDGET, "thread_limits": THREAD_LIMITS,
        "interpretation": "Production parent structure check only; no rate or lifetime claim.",
    })
    with (output / "preregistration.json").open("xb") as stream:
        stream.write(raw)
    started = time.monotonic()
    outcome = {"status": "failed", "budget": BUDGET}

    def deadline(_signum, _frame):
        raise KeyboardInterrupt("P1 parent check exhausted its overall wall budget")

    previous = signal.signal(signal.SIGALRM, deadline)
    signal.setitimer(signal.ITIMER_REAL, BUDGET["overall_timeout_seconds"])
    try:
        evidence = collect_orca_evidence(
            smiles, output, charge=0, multiplicity=1, ncores=BUDGET["maximum_cores"],
            method_profile=profile, timeout_seconds=BUDGET["orca_timeout_seconds_per_job"],
        )
        save(output, "parent-evidence.json", evidence)
        outcome.update(status=evidence["status"], local_minimum=evidence.get("local_minimum"),
                       failure_reason=evidence.get("failure_reason"))
    except BaseException as error:
        outcome.update(failure_reason=f"{type(error).__name__}: {error}",
                       traceback=traceback.format_exc())
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
        outcome["elapsed_seconds"] = time.monotonic() - started
        outcome["manifest_unchanged"] = digest(Path(manifest_path).read_bytes()) == digest(raw)
        save(output, "result.json", outcome)
    print(json.dumps(outcome, sort_keys=True))
    return outcome


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.manifest, args.output_dir)
    raise SystemExit(0 if result["status"] == "completed" and result.get("local_minimum") else 1)
