"""Prepare/freeze/run a bounded known-reaction diagnostic; never a fate job.

Preparation and freezing launch no calculations. Native execution requires a
separately authorized immutable plan and the user's existing --rmg-env.
"""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from storca.rate_transfer_execution import (prepare_execution, freeze_execution,
                                            reopen_execution, run_execution, prepare_h_ch4_campaign)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare", help="zero-calculation execution draft")
    prepare.add_argument("preparation", type=Path)
    prepare.add_argument("case", choices=("h_ch4", "oh_h2"))
    prepare.add_argument("destination", type=Path)
    candidate = commands.add_parser("prepare-h-ch4", help="freeze conditional counting/runtime, no authorization or calculations")
    candidate.add_argument("preparation", type=Path)
    candidate.add_argument("runtime_source", type=Path)
    candidate.add_argument("destination", type=Path)
    freeze = commands.add_parser("freeze", help="seal an explicitly completed plan; no execution")
    freeze.add_argument("plan", type=Path)
    freeze.add_argument("destination", type=Path)
    inspect = commands.add_parser("inspect", help="reopen source bindings; no execution")
    inspect.add_argument("workdir", type=Path)
    archive = commands.add_parser("archive-sources", help="archive exact historical sources BEFORE editing bound code; no execution")
    archive.add_argument("original_campaign", type=Path)
    archive.add_argument("source_capsule", type=Path)
    continuation = commands.add_parser("prepare-continuation", help="prepare same-ledger adaptive continuation, UNAUTHORIZED and zero native work")
    continuation.add_argument("original_campaign", type=Path)
    continuation.add_argument("source_capsule", type=Path)
    continuation.add_argument("destination", type=Path)
    continuation_freeze = commands.add_parser("freeze-continuation", help="seal an explicit continuation plan without running")
    continuation_freeze.add_argument("plan", type=Path)
    continuation_freeze.add_argument("destination", type=Path)
    continuation_inspect = commands.add_parser("inspect-continuation", help="reopen archived/live lineage and original native journal")
    continuation_inspect.add_argument("workdir", type=Path)
    for name in ("run", "resume"):
        run = commands.add_parser(name, help="only a separately authorized sealed campaign")
        run.add_argument("workdir", type=Path)
        run.add_argument("--rmg-env", required=True)
    for name in ("run-continuation", "resume-continuation"):
        run = commands.add_parser(name, help="only a separately authorized same-ledger continuation")
        run.add_argument("workdir", type=Path)
        run.add_argument("--rmg-env", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        result = prepare_execution(args.preparation, args.case, args.destination)
    elif args.command == "prepare-h-ch4":
        result = prepare_h_ch4_campaign(args.preparation, args.runtime_source, args.destination)
    elif args.command == "freeze":
        result = freeze_execution(args.plan, args.destination)
    elif args.command == "inspect":
        contract = reopen_execution(args.workdir)
        result = {"status": "bindings_reopened", "authorized": contract["plan"]["control"]["authorized"],
                  "native_jobs_launched": 0}
    elif args.command == "archive-sources":
        from archive_rate_transfer_sources import archive_sources
        value = archive_sources(args.original_campaign, args.source_capsule)
        result = {"status": "historical_sources_archived", "files": len(value["files"]),
                  "native_jobs_launched": 0}
    elif args.command in {"prepare-continuation", "freeze-continuation", "inspect-continuation",
                           "run-continuation", "resume-continuation"}:
        from storca.rate_transfer_continuation import (prepare_continuation, freeze_continuation,
                                                       reopen_continuation, run_continuation)
        if args.command == "prepare-continuation":
            result = prepare_continuation(args.original_campaign, args.source_capsule, args.destination)
        elif args.command == "freeze-continuation":
            result = freeze_continuation(args.plan, args.destination)
        elif args.command == "inspect-continuation":
            contract = reopen_continuation(args.workdir)
            result = {"status": "continuation_bindings_reopened",
                      "authorized": contract["plan"]["authorized"], "native_jobs_launched": 0,
                      "original_control_contract": contract["plan"]["native_control_contract"]}
        else:
            result = run_continuation(args.workdir, rmg_env=args.rmg_env)
    else:
        result = run_execution(args.workdir, rmg_env=args.rmg_env)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
