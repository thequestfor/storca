"""Freeze a preparation-only known-reaction benchmark; never launch chemistry."""

from pathlib import Path
import argparse
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from storca.rate_transfer_preparation import prepare_transfer_benchmark


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    report = prepare_transfer_benchmark(args.request, args.output)
    print(json.dumps({key: report[key] for key in (
        "status", "execution_ready", "native_jobs_launched", "readiness_gaps",
    )}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
