"""Parent-side client for the isolated RMG/Cantera artifact bridge."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


def run_rmg_bridge(payload: dict, *, rmg_env: str | None, timeout_seconds: float = 300.0) -> dict:
    """Run the bridge under the selected RMG environment and decode JSON only."""
    script = Path(__file__).with_name("rmg_bridge.py").resolve()
    # conda run does not reliably relay stdin on every supported platform.
    # A retained temporary JSON file keeps the protocol byte-for-byte stable.
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", encoding="utf-8") as handle:
        json.dump(payload, handle)
        handle.flush()
        command = (["conda", "run", "-n", rmg_env, "python", str(script), "--payload", handle.name] if rmg_env
                   else [sys.executable, str(script), "--payload", handle.name])
        completed = subprocess.run(command, text=True, capture_output=True, timeout=timeout_seconds)
    if completed.returncode:
        raise RuntimeError(f"RMG bridge failed with exit code {completed.returncode}: {completed.stderr.strip()}")
    try:
        result = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(f"RMG bridge emitted invalid JSON: {completed.stdout[-500:]}") from error
    if result.get("status") == "failed":
        raise RuntimeError(f"RMG bridge failed: {result.get('failure_reason', 'unknown failure')}")
    return result
