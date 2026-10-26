"""Immutable, source-bound records of individual bounded RMG attempts.

A retained attempt is historical evidence, not a convergence certificate.
Loading it preserves its original failure, rate and coverage states; callers
must independently decide whether a candidate is eligible for verification.
Seed directories are recorded only as locations and are never authorized for
restart by this module.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


SCHEMA_VERSION = 1
SNAPSHOT_KIND = "storca_rmg_attempt_evidence"
REFERENCE_KIND = "storca_rmg_attempt_reference"
FILE_ROLES = (
    "input", "stdout", "stderr", "log", "chemkin", "species_dictionary",
    "transport", "collision_rate_violators",
)


class RMGAttemptEvidenceError(ValueError):
    """The retained record or the source evidence cannot be admitted."""


def _encoded(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _path(value: object) -> Path:
    if not isinstance(value, str) or not value:
        raise RMGAttemptEvidenceError("rmg_attempt_artifact_path_invalid")
    return Path(value).absolute()


def _file_record(value: object) -> dict:
    if value is None:
        return {"status": "not_declared"}
    path = _path(value)
    record = {"path": str(path), "resolved_path": str(path.resolve())}
    if not path.exists():
        return {**record, "status": "missing"}
    if not path.is_file():
        return {**record, "status": "not_a_file"}
    return {**record, "status": "present", "sha256": _file_digest(path)}


def _directory_record(value: object, *, solver: bool = False) -> dict:
    if value is None:
        return {"status": "not_declared"}
    path = _path(value)
    record = {"path": str(path), "resolved_path": str(path.resolve())}
    if not path.exists():
        return {**record, "status": "missing"}
    if not path.is_dir():
        return {**record, "status": "not_a_directory"}
    if solver:
        return {**record, "status": "present", "scope": "simulation_csv_files",
                "files": [_file_record(str(item)) for item in sorted(path.glob("simulation_*.csv"))]}
    return {**record, "status": "present", "scope": "location_only_not_restart_evidence"}


def _normalized_evidence(evidence: dict) -> dict:
    if not isinstance(evidence, dict):
        raise RMGAttemptEvidenceError("rmg_attempt_evidence_invalid")
    # A strict JSON roundtrip makes an independent copy and refuses NaN/Inf.
    copied = json.loads(_encoded({key: value for key, value in evidence.items() if key != "controller"}))
    artifacts = copied.get("artifacts") or {}
    if not isinstance(artifacts, dict):
        raise RMGAttemptEvidenceError("rmg_attempt_artifacts_invalid")
    for role, value in artifacts.items():
        if value is not None:
            artifacts[role] = str(_path(value))
    if "artifacts" in copied:
        copied["artifacts"] = artifacts
    for field in ("solver_profile", "kinetics_validation"):
        record = copied.get(field)
        if isinstance(record, dict) and record.get("path") is not None:
            record["path"] = str(_path(record["path"]))
    for record in _identity_records(copied).values():
        if not isinstance(record, dict):
            raise RMGAttemptEvidenceError("rmg_attempt_identity_artifact_invalid")
        _path(record.get("path"))
        # Identity sidecar bytes already bind the path spelling. Rewriting a
        # relative proof path here would invalidate that earlier evidence.
        if not Path(record["path"]).is_absolute():
            raise RMGAttemptEvidenceError("rmg_attempt_identity_path_not_absolute")
    for _, manifest, _ in _library_manifests(copied):
        if manifest.get("library") is not None:
            manifest["library"] = str(_path(manifest["library"]))
        for role, value in (manifest.get("artifacts") or {}).items():
            manifest["artifacts"][role] = str(_path(value))
    return copied


def _identity_records(evidence: dict) -> dict:
    inspection = evidence.get("mechanism_inspection") or {}
    if not isinstance(inspection, dict):
        raise RMGAttemptEvidenceError("rmg_attempt_inspection_invalid")
    records = {}
    sidecar = inspection.get("replacement_identity_artifact")
    if sidecar is not None:
        records["sidecar"] = sidecar
    report = inspection.get("replacement_identity")
    if isinstance(report, dict):
        for name, record in (report.get("mechanism_artifacts") or {}).items():
            records["mechanism/" + name] = record
        for index, correspondence in enumerate(report.get("correspondences") or []):
            for name, record in (correspondence.get("source_artifacts") or {}).items():
                records[f"source/{index}/{name}"] = record
    return records


def _library_manifests(evidence: dict) -> list[tuple[str, dict, str]]:
    manifests = [(f"generated/{index}", manifest, "generated-kinetics.json")
                 for index, manifest in enumerate(evidence.get("generated_kinetics_libraries") or [])]
    subset = (evidence.get("reference_kinetics_libraries") or {}).get("direct_subset")
    if subset is not None:
        manifests.append(("reference_subset", subset, "reference-subset.json"))
    if any(not isinstance(manifest, dict) for _, manifest, _ in manifests):
        raise RMGAttemptEvidenceError("rmg_attempt_library_manifest_invalid")
    return manifests


def _library_records(evidence: dict) -> dict:
    records = {}
    for prefix, manifest, manifest_name in _library_manifests(evidence):
        for role, value in (manifest.get("artifacts") or {}).items():
            records[prefix + "/artifact/" + role] = {"path": value}
        if manifest.get("library") is None:
            if manifest.get("checksums"):
                raise RMGAttemptEvidenceError("rmg_attempt_library_path_missing")
            continue
        library = _path(manifest["library"])
        records[prefix + "/manifest"] = {"path": str(library / manifest_name)}
        for name in {"reactions.py", "dictionary.txt", *(manifest.get("checksums") or {})}:
            if not isinstance(name, str) or Path(name).name != name or name in {".", ".."}:
                raise RMGAttemptEvidenceError("rmg_attempt_library_checksum_path_invalid")
            records[prefix + "/" + name] = {"path": str(library / name)}
            if name in (manifest.get("checksums") or {}):
                records[prefix + "/" + name]["sha256"] = manifest["checksums"][name]
    return records


def _inventory(evidence: dict) -> dict:
    artifacts = evidence.get("artifacts") or {}
    inventory = {"artifacts/" + role: _file_record(artifacts.get(role)) for role in FILE_ROLES}
    for role, value in artifacts.items():
        if role not in {*FILE_ROLES, "solver", "seed"}:
            inventory["artifacts/" + role] = _file_record(value)
    inventory["artifacts/solver"] = _directory_record(artifacts.get("solver"), solver=True)
    inventory["artifacts/seed"] = _directory_record(artifacts.get("seed"))
    for field in ("solver_profile", "kinetics_validation"):
        record = evidence.get(field)
        if isinstance(record, dict) and record.get("path") is not None:
            inventory[field + "/path"] = _file_record(record["path"])
    for role, record in _identity_records(evidence).items():
        if not isinstance(record, dict) or not record.get("path"):
            raise RMGAttemptEvidenceError("rmg_attempt_identity_artifact_invalid")
        inventory["identity/" + role] = _file_record(record["path"])
    for role, record in _library_records(evidence).items():
        inventory["library/" + role] = _file_record(record["path"])
    return inventory


def _validate_identity(evidence: dict) -> None:
    for record in _identity_records(evidence).values():
        current = _file_record(record["path"])
        if current["status"] != "present" or current["sha256"] != record.get("sha256"):
            raise RMGAttemptEvidenceError("rmg_attempt_identity_source_changed")
    inspection = evidence.get("mechanism_inspection") or {}
    if "replacement_identity" in inspection:
        sidecar = inspection.get("replacement_identity_artifact")
        if not isinstance(sidecar, dict):
            raise RMGAttemptEvidenceError("rmg_attempt_identity_sidecar_missing")
        saved = json.loads(Path(sidecar["path"]).read_bytes())
        if _encoded(saved) != _encoded(inspection["replacement_identity"]):
            raise RMGAttemptEvidenceError("rmg_attempt_identity_sidecar_mismatch")
    for record in _library_records(evidence).values():
        if "sha256" in record:
            current = _file_record(record["path"])
            if current["status"] != "present" or current["sha256"] != record["sha256"]:
                raise RMGAttemptEvidenceError("rmg_attempt_library_source_changed")


def retain_rmg_attempt(
    evidence: dict,
    directory: Path,
    *,
    attempt_number: int,
    condition_contract: dict | None,
    predecessor: dict | None = None,
) -> dict:
    """Retain a complete parsed attempt before another attempt can start.

    Missing outputs are explicit states and can be retained for failed jobs.
    Repeating this operation with identical evidence reuses the same bytes.
    The caller owns the attempt sequence; predecessor is a compact provenance
    reference rather than recursively embedded historical evidence.
    """
    try:
        if type(attempt_number) is not int or attempt_number < 1:
            raise RMGAttemptEvidenceError("rmg_attempt_number_invalid")
        if condition_contract is not None and not isinstance(condition_contract, dict):
            raise RMGAttemptEvidenceError("rmg_attempt_condition_invalid")
        copied = _normalized_evidence(evidence)
        condition_hash = _digest(_encoded(condition_contract))
        if predecessor is not None:
            _validate_reference(predecessor)
            if (predecessor["attempt_number"] >= attempt_number
                    or predecessor["condition_sha256"] != condition_hash):
                raise RMGAttemptEvidenceError("rmg_attempt_predecessor_invalid")
        snapshot = {
            "kind": SNAPSHOT_KIND, "schema_version": SCHEMA_VERSION,
            "attempt_number": attempt_number, "condition_contract": condition_contract,
            "condition_sha256": condition_hash, "predecessor": predecessor,
            "source_artifacts": _inventory(copied), "evidence": copied,
        }
        content = _encoded(snapshot) + b"\n"
        digest = _digest(content)
        path = Path(directory).absolute() / f"rmg-attempt-{attempt_number:02d}-{digest}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with path.open("xb") as stream:
                stream.write(content)
                stream.flush()
                # The reference must not precede durable evidence on resume.
                os.fsync(stream.fileno())
        except FileExistsError:
            if path.read_bytes() != content:
                raise RMGAttemptEvidenceError("rmg_attempt_snapshot_mutated")
        return {"kind": REFERENCE_KIND, "schema_version": SCHEMA_VERSION,
                "attempt_number": attempt_number, "path": str(path), "sha256": digest,
                "condition_sha256": condition_hash}
    except RMGAttemptEvidenceError:
        raise
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as error:
        raise RMGAttemptEvidenceError(f"rmg_attempt_retention_failed: {error}") from error


def _validate_reference(reference: dict) -> None:
    if (not isinstance(reference, dict)
            or reference.get("kind") != REFERENCE_KIND
            or reference.get("schema_version") != SCHEMA_VERSION
            or type(reference.get("attempt_number")) is not int
            or reference["attempt_number"] < 1
            or not isinstance(reference.get("sha256"), str)
            or not isinstance(reference.get("condition_sha256"), str)
            or not Path(reference.get("path", "")).is_absolute()):
        raise RMGAttemptEvidenceError("rmg_attempt_reference_invalid")


def load_rmg_attempt(reference: dict, *, condition_contract: dict | None) -> dict:
    """Verify a retained snapshot, exact conditions and all inventoried files.

    A file appearing after a failed attempt is also drift: it cannot silently
    upgrade that attempt's evidence. Invalid or stale records raise a named
    error; callers must not replace them with unverified parsed data.
    """
    try:
        _validate_reference(reference)
        content = Path(reference["path"]).read_bytes()
        if _digest(content) != reference["sha256"]:
            raise RMGAttemptEvidenceError("rmg_attempt_snapshot_hash_mismatch")
        snapshot = json.loads(content)
        if (not isinstance(snapshot, dict) or snapshot.get("kind") != SNAPSHOT_KIND
                or snapshot.get("schema_version") != SCHEMA_VERSION
                or snapshot.get("attempt_number") != reference["attempt_number"]
                or _encoded(snapshot) + b"\n" != content):
            raise RMGAttemptEvidenceError("rmg_attempt_snapshot_invalid")
        condition_hash = _digest(_encoded(condition_contract))
        if (snapshot.get("condition_sha256") != condition_hash
                or reference["condition_sha256"] != condition_hash
                or _encoded(snapshot.get("condition_contract")) != _encoded(condition_contract)):
            raise RMGAttemptEvidenceError("rmg_attempt_condition_mismatch")
        evidence = snapshot["evidence"]
        if not isinstance(evidence, dict) or "controller" in evidence:
            raise RMGAttemptEvidenceError("rmg_attempt_evidence_invalid")
        if _inventory(evidence) != snapshot.get("source_artifacts"):
            raise RMGAttemptEvidenceError("rmg_attempt_source_artifacts_changed")
        _validate_identity(evidence)
        return evidence
    except RMGAttemptEvidenceError:
        raise
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as error:
        raise RMGAttemptEvidenceError(f"rmg_attempt_load_failed: {error}") from error
