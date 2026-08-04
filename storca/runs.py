"""Creation and bookkeeping of isolated calculation runs."""

from __future__ import annotations

import json
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def _safe_name(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip(".-")
    return name or "calculation"


def create_run_directory(output_root: Path, name: str) -> Path:
    """Create a unique run directory below *output_root*."""
    output_root = Path(output_root)
    timestamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    candidate = output_root / f"{_safe_name(name)}-{timestamp}"
    suffix = 1
    while candidate.exists():
        candidate = output_root / f"{_safe_name(name)}-{timestamp}-{suffix}"
        suffix += 1
    candidate.mkdir(parents=True)
    return candidate


def write_metadata(run_dir: Path, **values: Any) -> Path:
    """Merge JSON-serializable values into persistent run metadata."""
    path = Path(run_dir) / "metadata.json"
    metadata = {}
    if path.is_file():
        metadata = json.loads(path.read_text())
    metadata.update(values)
    metadata.setdefault("created_at", datetime.now(UTC).isoformat())
    metadata["updated_at"] = datetime.now(UTC).isoformat()
    path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return path
