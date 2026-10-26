"""Content-addressed evidence resolution shared by campaign adapters."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from .stability_campaign import CampaignLaneRequest


def evidence_content_hash(path: Path) -> str:
    """Hash one evidence file or a deterministic nonempty directory tree."""
    path = Path(path)
    digest = hashlib.sha256()
    if path.is_file():
        digest.update(path.read_bytes())
    elif path.is_dir():
        files = sorted(item for item in path.rglob("*") if item.is_file())
        if not files:
            raise ValueError(f"Evidence directory is empty: {path}")
        for item in files:
            digest.update(item.relative_to(path).as_posix().encode("utf-8"))
            digest.update(b"\0")
            digest.update(item.read_bytes())
            digest.update(b"\0")
    else:
        raise ValueError(f"Evidence input does not exist: {path}")
    return "sha256:" + digest.hexdigest()


def resolve_lane_evidence(
    lane: "CampaignLaneRequest",
    *,
    base_dir: Path,
    require_paths: bool = True,
) -> dict[str, list[dict[str, Any]]]:
    """Resolve and verify every retained evidence reference for one lane."""
    resolved: dict[str, list[dict[str, Any]]] = {}
    for reference in lane.evidence_inputs:
        if reference.path is None:
            if require_paths:
                raise ValueError(
                    "Campaign execution needs a resolvable path for evidence input "
                    f"{reference.input_id}"
                )
            continue
        path = Path(reference.path)
        if not path.is_absolute():
            path = Path(base_dir) / path
        actual = evidence_content_hash(path)
        if actual != reference.content_sha256:
            raise ValueError(
                f"Evidence input {reference.input_id} content hash does not match: "
                f"declared={reference.content_sha256}, actual={actual}"
            )
        resolved.setdefault(reference.input_kind, []).append({
            "input_id": reference.input_id,
            "input_kind": reference.input_kind,
            "path": path,
            "content_sha256": actual,
            "claim_id": reference.claim_id,
        })
    return resolved


def evidence_contract_records(
    evidence: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    """Return path-independent evidence identities for execution contracts."""
    return [
        {
            "input_id": item["input_id"],
            "input_kind": item["input_kind"],
            "content_sha256": item["content_sha256"],
            "claim_id": item["claim_id"],
        }
        for kind in sorted(evidence)
        for item in evidence[kind]
    ]


__all__ = [
    "evidence_content_hash", "evidence_contract_records", "resolve_lane_evidence",
]
