"""Exact-contract continuation for retained gas stability verification runs."""

from __future__ import annotations

from datetime import UTC, datetime
import ast
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time
from typing import Any, Iterable

from .campaign_evidence import evidence_content_hash
from .conditions import build_condition_spec
from .kinetic_uncertainty import coerce_rate_uncertainty_evidence
from .runs import write_metadata
from .verification_engine import VerificationEngineConfig, continue_verification_engine
from .surface_resolution import SurfaceResolutionPolicy
from .verification_policy import (
    DISSOCIATION_ESCALATION_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    ENDPOINT_RETURN_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    LEGACY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    PREVIOUS_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    PRE_IRC_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    PRE_RING_CONNECTIVITY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    PRE_GEOMETRY_AUDIT_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    PRE_REACHED_MINIMA_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
)


EXECUTION_CONTRACT_NAME = "verification-execution-contract.json"
EXECUTION_HEARTBEAT_QUIET_SECONDS = 15 * 60
EXECUTION_RESOURCE_AMENDMENT_PREFIX = "verification-execution-resource-amendment"
PATH_POLICY_AMENDMENT_PREFIX = "verification-path-policy-amendment"
SESSION_RECONCILIATION_DIRECTORY = "verification-session-reconciliations"


def _canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _evidence_input_hash(path: Path) -> str:
    """Hash a library input, including an explicitly empty directory."""
    path = Path(path)
    if path.is_dir() and not any(item.is_file() for item in path.rglob("*")):
        return "sha256:" + hashlib.sha256(b"").hexdigest()
    return evidence_content_hash(path)


def _write_json_atomic(path: Path, value: Any) -> Path:
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


def build_verification_execution_contract(
    *,
    smiles: str,
    charge: int,
    multiplicity: int,
    ncores: int,
    method_profile: dict,
    condition_contract: dict,
    scenario_config: dict,
    resource_config: dict,
    barrier_threshold: float,
    rmg_env: str | None,
    rmg_command: str,
    rmg_maximum_heavy_atoms: int | None,
    rmg_retry_attempts: int,
    rmg_hard_timeout_seconds: float | None,
    rmg_heartbeat_timeout_seconds: float | None,
    rmg_start_method: str,
    database_reaction_libraries: Iterable[str],
    reaction_libraries: Iterable[Path],
    verification_max_iterations: int,
    verification_timeout_seconds: float | None,
    verification_orientations: int,
    verification_neb_images: int,
) -> dict:
    """Capture every setting needed to rebuild the verification callbacks."""
    base_libraries = []
    for raw_path in reaction_libraries:
        path = Path(raw_path).resolve()
        base_libraries.append({
            "path": str(path),
            "content_sha256": _evidence_input_hash(path),
        })
    path_policy = SurfaceResolutionPolicy().as_dict()
    body = {
        "schema_version": 1,
        "kind": "stability_verification_execution_contract",
        "contract_status": "fully_retained",
        "parent": {
            "smiles": smiles,
            "charge": int(charge),
            "multiplicity": int(multiplicity),
        },
        "condition_contract": condition_contract,
        "condition_contract_sha256": _canonical_hash(condition_contract),
        "method_profile": method_profile,
        "resources": {
            "ncores": int(ncores),
            "route_job_timeout_seconds": verification_timeout_seconds,
            "orientations": int(verification_orientations),
            "neb_images": int(verification_neb_images),
        },
        "path_verifier_policy_id": GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        "path_verifier_policy": path_policy,
        "path_verifier_policy_sha256": path_policy["scientific_hash"],
        "verification_engine_config": {
            "max_iterations": int(verification_max_iterations),
        },
        "rmg": {
            "environment": rmg_env,
            "command": rmg_command,
            "barrier_threshold_kcal_mol": float(barrier_threshold),
            "resource_config": resource_config,
            "maximum_heavy_atoms": rmg_maximum_heavy_atoms,
            "retry_attempts": int(rmg_retry_attempts),
            "hard_timeout_seconds": rmg_hard_timeout_seconds,
            "heartbeat_timeout_seconds": rmg_heartbeat_timeout_seconds,
            "start_method": rmg_start_method,
            "scenario_config": scenario_config,
            "database_reaction_libraries": list(database_reaction_libraries),
            "base_reaction_libraries": base_libraries,
        },
    }
    return {**body, "execution_contract_sha256": _canonical_hash(body)}


def write_verification_execution_contract(run_dir: Path, contract: dict) -> Path:
    return _write_json_atomic(Path(run_dir) / EXECUTION_CONTRACT_NAME, contract)


def _effective_execution_contract(run_dir: Path, contract: dict) -> tuple[dict, list[dict]]:
    effective = copy.deepcopy(contract)
    amendments = []
    previous_hash = None
    for path in sorted(Path(run_dir).glob(f"{EXECUTION_RESOURCE_AMENDMENT_PREFIX}-*.json")):
        amendment = json.loads(path.read_text())
        retained_hash = amendment.get("amendment_sha256")
        unsigned = {key: value for key, value in amendment.items() if key != "amendment_sha256"}
        if retained_hash != _canonical_hash(unsigned):
            raise ValueError(f"Execution resource amendment has been mutated: {path}")
        if amendment.get("execution_contract_sha256") != contract.get("execution_contract_sha256"):
            raise ValueError(f"Execution resource amendment targets another contract: {path}")
        if amendment.get("predecessor_amendment_sha256") != previous_hash:
            raise ValueError(f"Execution resource amendment chain is broken: {path}")
        current = int(effective["resources"]["ncores"])
        if int(amendment.get("prior_ncores")) != current:
            raise ValueError(f"Execution resource amendment prior core count is inconsistent: {path}")
        new_ncores = int(amendment.get("new_ncores"))
        if new_ncores < 1:
            raise ValueError(f"Execution resource amendment has an invalid core count: {path}")
        effective["resources"]["ncores"] = new_ncores
        amendments.append(amendment)
        previous_hash = retained_hash
    return effective, amendments


def write_execution_resource_amendment(run_dir: Path, contract: dict, *, ncores: int) -> dict:
    """Append an auditable execution-only core-count change."""
    if int(ncores) < 1:
        raise ValueError("ORCA core count must be at least one")
    effective, amendments = _effective_execution_contract(run_dir, contract)
    prior = int(effective["resources"]["ncores"])
    if prior == int(ncores):
        return amendments[-1] if amendments else {
            "kind": "no_execution_resource_change", "new_ncores": prior,
        }
    number = len(amendments) + 1
    body = {
        "schema_version": 1,
        "kind": "verification_execution_resource_amendment",
        "amendment_number": number,
        "created_at": datetime.now(UTC).isoformat(),
        "execution_contract_sha256": contract["execution_contract_sha256"],
        "predecessor_amendment_sha256": (
            amendments[-1]["amendment_sha256"] if amendments else None
        ),
        "condition_contract_sha256": contract["condition_contract_sha256"],
        "method_profile_sha256": _canonical_hash(contract["method_profile"]),
        "prior_ncores": prior,
        "new_ncores": int(ncores),
        "scope": "execution_parallelism_only",
        "scientific_method_changed": False,
        "reason": "explicit_user_requested_resume_resource_change",
    }
    amendment = {**body, "amendment_sha256": _canonical_hash(body)}
    _write_json_atomic(
        Path(run_dir) / f"{EXECUTION_RESOURCE_AMENDMENT_PREFIX}-{number:03d}.json",
        amendment,
    )
    return amendment


def _effective_path_policy(run_dir: Path, contract: dict) -> tuple[str, list[dict]]:
    """Validate the append-only scientific-policy chain and return its tip."""
    effective = str(contract.get("path_verifier_policy_id") or "").strip()
    if not effective:
        raise ValueError("Execution contract has no path-verifier policy ID")
    amendments: list[dict] = []
    previous_hash = None
    for path in sorted(Path(run_dir).glob(f"{PATH_POLICY_AMENDMENT_PREFIX}-*.json")):
        amendment = json.loads(path.read_text())
        retained_hash = amendment.get("amendment_sha256")
        unsigned = {key: value for key, value in amendment.items() if key != "amendment_sha256"}
        if retained_hash != _canonical_hash(unsigned):
            raise ValueError(f"Path-policy amendment has been mutated: {path}")
        if amendment.get("execution_contract_sha256") != contract.get("execution_contract_sha256"):
            raise ValueError(f"Path-policy amendment targets another contract: {path}")
        if amendment.get("predecessor_amendment_sha256") != previous_hash:
            raise ValueError(f"Path-policy amendment chain is broken: {path}")
        if str(amendment.get("prior_path_verifier_policy_id")) != effective:
            raise ValueError(f"Path-policy amendment predecessor is inconsistent: {path}")
        successor = str(amendment.get("new_path_verifier_policy_id") or "").strip()
        if not successor or successor == effective:
            raise ValueError(f"Path-policy amendment has an invalid successor: {path}")
        effective = successor
        amendments.append(amendment)
        previous_hash = retained_hash
    return effective, amendments


def write_path_policy_amendment(
    run_dir: Path,
    contract: dict,
    *,
    path_verifier_policy_id: str,
) -> dict:
    """Append an explicit, hashed scientific path-policy migration."""
    successor = str(path_verifier_policy_id).strip()
    if not successor:
        raise ValueError("A path-verifier policy amendment needs a successor policy ID")
    effective, amendments = _effective_path_policy(run_dir, contract)
    if successor == effective:
        return amendments[-1] if amendments else {
            "kind": "no_path_policy_change",
            "new_path_verifier_policy_id": effective,
        }
    if effective in {
        LEGACY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        PREVIOUS_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        ENDPOINT_RETURN_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        DISSOCIATION_ESCALATION_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        PRE_IRC_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        PRE_RING_CONNECTIVITY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        PRE_GEOMETRY_AUDIT_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        PRE_REACHED_MINIMA_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    } and successor == GENERIC_ORCA_PATH_VERIFIER_POLICY_ID:
        successor_policy = SurfaceResolutionPolicy().as_dict()
    else:
        raise ValueError(
            "Unsupported automatic scientific-policy migration: "
            f"{effective} -> {successor}"
        )
    number = len(amendments) + 1
    body = {
        "schema_version": 1,
        "kind": "verification_path_policy_amendment",
        "amendment_number": number,
        "created_at": datetime.now(UTC).isoformat(),
        "execution_contract_sha256": contract["execution_contract_sha256"],
        "predecessor_amendment_sha256": (
            amendments[-1]["amendment_sha256"] if amendments else None
        ),
        "condition_contract_sha256": contract["condition_contract_sha256"],
        "method_profile_sha256": _canonical_hash(contract["method_profile"]),
        "prior_path_verifier_policy_id": effective,
        "new_path_verifier_policy_id": successor,
        "new_path_verifier_policy_sha256": successor_policy["scientific_hash"],
        "new_path_verifier_policy": successor_policy,
        "scope": "reaction_surface_resolution_policy",
        "scientific_method_changed": True,
        "reason": "bounded_surface_resolution_ladder_policy_upgrade",
        "claims_permitted": ["retry_previously_deferred_route_under_new_policy"],
        "claims_prohibited": ["mechanism", "rate", "lifetime"],
    }
    amendment = {**body, "amendment_sha256": _canonical_hash(body)}
    _write_json_atomic(
        Path(run_dir) / f"{PATH_POLICY_AMENDMENT_PREFIX}-{number:03d}.json",
        amendment,
    )
    return amendment


def _legacy_execution_contract(run_dir: Path, result: dict, metadata: dict) -> dict:
    """Make old retained defaults explicit; never infer local input libraries."""
    reference = (result.get("rmg_evidence") or {}).get("reference_kinetics_libraries") or {}
    if reference.get("requested_source_libraries"):
        raise ValueError(
            "Legacy run used local reaction libraries whose paths/hashes were not retained; "
            "exact-contract continuation is unavailable"
        )
    required = (
        "smiles", "charge", "multiplicity", "ncores", "harmonic_method_profile",
        "condition_contract", "stability_scenario", "rmg_resource_tier",
        "rmg_command", "route_energy_ceiling_kcal_mol",
    )
    missing = [key for key in required if metadata.get(key) is None]
    if missing:
        raise ValueError("Legacy run lacks required continuation metadata: " + ", ".join(missing))
    if metadata["condition_contract"] != result.get("condition_contract"):
        raise ValueError("Legacy run metadata and result condition contracts disagree")

    endpoint_metadata = sorted(run_dir.glob(
        "verification-engine*/iteration-*/path/stationary-points/*/*/metadata.json"
    ))
    if not endpoint_metadata:
        raise ValueError("Legacy verification has no retained endpoint metadata")
    endpoint = json.loads(endpoint_metadata[0].read_text())
    if endpoint.get("ncores") != metadata["ncores"] \
            or endpoint.get("orca_method_keywords") != metadata["harmonic_method_profile"].get("orca_keywords"):
        raise ValueError("Legacy endpoint execution does not match retained parent method/resources")

    evidence = result.get("rmg_evidence") or {}
    limits = evidence.get("resource_limits") or metadata["rmg_resource_tier"]
    body = build_verification_execution_contract(
        smiles=metadata["smiles"], charge=metadata["charge"],
        multiplicity=metadata["multiplicity"], ncores=metadata["ncores"],
        method_profile=metadata["harmonic_method_profile"],
        condition_contract=metadata["condition_contract"],
        scenario_config=metadata["stability_scenario"],
        resource_config={
            "name": metadata["rmg_resource_tier"].get("name", "retained-legacy-tier"),
            "walltime": limits.get("walltime"),
            "max_processes": limits.get("max_processes"),
            "max_iterations": limits.get("max_iterations"),
            "max_edge_species": limits.get("max_edge_species"),
            "max_objects_per_iteration": metadata["rmg_resource_tier"].get("max_objects_per_iteration", 1),
        },
        barrier_threshold=metadata["route_energy_ceiling_kcal_mol"],
        rmg_env=metadata.get("rmg_env"), rmg_command=metadata["rmg_command"],
        rmg_maximum_heavy_atoms=metadata.get("rmg_maximum_heavy_atoms"),
        rmg_retry_attempts=2, rmg_hard_timeout_seconds=None,
        rmg_heartbeat_timeout_seconds=None, rmg_start_method="auto",
        database_reaction_libraries=metadata.get("rmg_database_reaction_libraries") or (),
        reaction_libraries=(), verification_max_iterations=12,
        verification_timeout_seconds=endpoint.get("timeout_seconds"),
        verification_orientations=3, verification_neb_images=8,
    )
    migrated = {
        **body,
        "contract_status": "retained_legacy_contract_with_explicit_default_inference",
        "legacy_inferred_fields": {
            "verification_engine_config.max_iterations": 12,
            "resources.orientations": 3,
            "resources.neb_images": 8,
            "rmg.retry_attempts": 2,
            "rmg.hard_timeout_seconds": None,
            "rmg.heartbeat_timeout_seconds": None,
            "rmg.start_method": "auto",
        },
        "legacy_inference_basis": (
            "Values are the CLI defaults at creation; method, cores, timeout, condition, "
            "parent state, RMG resource limits, and absence of local libraries are retained independently."
        ),
    }
    unsigned = {key: value for key, value in migrated.items() if key != "execution_contract_sha256"}
    migrated["execution_contract_sha256"] = _canonical_hash(unsigned)
    return migrated


def _load_execution_contract(run_dir: Path, result: dict, metadata: dict) -> dict:
    path = run_dir / EXECUTION_CONTRACT_NAME
    if path.is_file():
        contract = json.loads(path.read_text())
    else:
        contract = _legacy_execution_contract(run_dir, result, metadata)
        write_verification_execution_contract(run_dir, contract)
    retained_hash = contract.get("execution_contract_sha256")
    unsigned = {key: value for key, value in contract.items() if key != "execution_contract_sha256"}
    if retained_hash != _canonical_hash(unsigned):
        raise ValueError("Verification execution contract has been mutated")
    if contract.get("condition_contract") != result.get("condition_contract"):
        raise ValueError("Execution contract does not match stability.json conditions")
    parent = contract.get("parent") or {}
    for key in ("smiles", "charge", "multiplicity"):
        if parent.get(key) != metadata.get(key):
            raise ValueError(f"Execution contract parent {key} does not match metadata")
    for record in (contract.get("rmg") or {}).get("base_reaction_libraries") or ():
        path = Path(record["path"])
        if _evidence_input_hash(path) != record.get("content_sha256"):
            raise ValueError(f"Base reaction library content changed: {path}")
    return contract


def _resolve_artifact_path(raw_path: object, run_dir: Path) -> Path:
    path = Path(str(raw_path))
    if path.is_absolute() or path.exists():
        return path
    run_dir = Path(run_dir).resolve()
    candidates = [run_dir / path]
    if run_dir.name in path.parts:
        marker = len(path.parts) - 1 - tuple(reversed(path.parts)).index(run_dir.name)
        candidates.append(run_dir / Path(*path.parts[marker + 1:]))
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return path


def _started_execution_heartbeat(
    state_path: Path,
    *,
    quiet_seconds: float = EXECUTION_HEARTBEAT_QUIET_SECONDS,
    observed_at_epoch: float | None = None,
) -> dict:
    """Classify a retained started calculation using artifact activity only.

    Process IDs are not portable across detached shells or machines.  A fresh
    calculation artifact is therefore positive evidence that an external
    engine may still be writing; STORCA fails closed and refuses a duplicate.
    """
    state_path = Path(state_path)
    state = json.loads(state_path.read_text())
    started = [
        item for item in (state.get("iterations") or [])
        if item.get("status") == "started" and item.get("route_id")
    ]
    if not started:
        raise ValueError("Retained running session has no started route to recover")
    observed = time.time() if observed_at_epoch is None else float(observed_at_epoch)
    newest: tuple[float, Path] | None = None
    workspaces = []
    live_leases = []
    dead_leases = []
    unproved_terminal_leases = []
    for iteration in started:
        raw = (iteration.get("artifacts") or {}).get("calculation_workspace") \
            or (iteration.get("artifacts") or {}).get("directory")
        workspace = _resolve_artifact_path(raw, state_path.parent.parent)
        if not workspace.is_dir():
            raise FileNotFoundError(f"Started calculation workspace is missing: {workspace}")
        workspaces.append(str(workspace))
        for artifact in workspace.rglob("*"):
            if not artifact.is_file():
                continue
            modified = artifact.stat().st_mtime
            if newest is None or modified > newest[0]:
                newest = (modified, artifact)
            if artifact.name.endswith(".process.json"):
                try:
                    lease = json.loads(artifact.read_text())
                    pid = int(lease.get("pid"))
                    process_group_id = int(lease.get("process_group_id") or pid)
                except (OSError, TypeError, ValueError, json.JSONDecodeError):
                    continue
                terminal_interruption = (lease.get("status") == "terminated"
                                         and lease.get("returncode") in {-2, -9})
                if lease.get("status") not in {
                    "running", "runner_detached_process_may_be_active",
                } and not terminal_interruption or pid < 1:
                    continue
                if terminal_interruption:
                    source = Path(str(lease.get("input") or ""))
                    ended = lease.get("ended_at_epoch")
                    if (lease.get("kind") != "storca_external_process_lease" or lease.get("engine") != "orca"
                            or not source.is_file() or lease.get("input_sha256") != hashlib.sha256(source.read_bytes()).hexdigest()
                            or isinstance(ended, bool) or not isinstance(ended, (int, float)) or not math.isfinite(ended)):
                        unproved_terminal_leases.append(str(artifact))
                        continue
                try:
                    os.kill(pid, 0)
                    alive = True
                except ProcessLookupError:
                    try:
                        os.killpg(process_group_id, 0)
                        alive = True
                    except ProcessLookupError:
                        alive = False
                    except PermissionError:
                        alive = True
                except PermissionError:
                    alive = True
                record = {
                    "lease": str(artifact), "pid": pid,
                    "process_group_id": process_group_id,
                    "status": lease.get("status"),
                    "ended_at_epoch": lease.get("ended_at_epoch"),
                }
                (live_leases if alive else dead_leases).append(record)
    if newest is None:
        raise ValueError("Started calculation workspace contains no heartbeat artifacts")
    age = max(0.0, observed - newest[0])
    # A missing recorded PID *and* process group is stronger evidence than a
    # fresh file timestamp: no retained process remains capable of writing.
    # PID reuse errs safely in the other direction because it appears live.
    stopped_lease_covers_latest = any(
        record.get("status") != "terminated" or record.get("ended_at_epoch", 0) + 1 >= newest[0]
        for record in dead_leases
    )
    confirmed_stopped_process = stopped_lease_covers_latest and not live_leases and not unproved_terminal_leases
    return {
        "kind": "calculation_artifact_heartbeat",
        "status": (
            "active" if live_leases or unproved_terminal_leases else
            "stale" if confirmed_stopped_process else
            "active" if age < quiet_seconds else "stale"
        ),
        "route_ids": [str(item["route_id"]) for item in started],
        "calculation_workspaces": workspaces,
        "newest_artifact": str(newest[1]),
        "newest_artifact_mtime_epoch": newest[0],
        "observed_at_epoch": observed,
        "artifact_age_seconds": age,
        "quiet_threshold_seconds": float(quiet_seconds),
        "classification_basis": (
            "live_or_unproved_external_process_lease" if unproved_terminal_leases
            else "live_external_process_lease" if live_leases
            else "confirmed_stopped_external_process_lease"
            if confirmed_stopped_process
            else "newest_retained_calculation_artifact_mtime"
        ),
        "live_process_leases": live_leases,
        "dead_process_leases": dead_leases,
        "unproved_terminal_process_leases": unproved_terminal_leases,
    }


def _recoverable_orca_startup_failure(state_path: Path, state: dict) -> dict | None:
    """Recognize a retained execution-environment failure, not chemistry."""
    latest_by_route = {}
    for iteration in state.get("iterations") or []:
        route_id = str(iteration.get("route_id") or "")
        if route_id:
            latest_by_route[route_id] = iteration
    eligible = []
    outputs = []
    markers = (
        "ORCA finished by error termination in Startup",
        "bind() failed for port 0: Operation not permitted",
        "No sockets were able to be opened on the available protocols",
    )
    for route_id, iteration in latest_by_route.items():
        if iteration.get("status") != "deferred_unverified":
            continue
        raw_workspace = (iteration.get("artifacts") or {}).get("calculation_workspace")
        workspace = _resolve_artifact_path(raw_workspace, state_path.parent.parent)
        failed_outputs = []
        retained_outputs: set[Path] = set()
        current_result_is_startup_failure = (
            iteration.get("path_verification_status") == "execution_infrastructure_failure"
        )
        raw_path_result = (iteration.get("artifacts") or {}).get("path_result")
        if raw_path_result:
            path_result = _resolve_artifact_path(raw_path_result, state_path.parent.parent)
            if path_result.is_file():
                route_result = json.loads(path_result.read_text())
                guided = (
                    (route_result.get("surface_resolution_guided_recovery") or {}).get("result")
                    or {}
                )
                # Compatibility for the first retained v3 run, before startup
                # failures received their own categorical route state.
                legacy_startup_misclassification = (
                    guided.get("status") == "coordinate_scan_no_ts_seed"
                    and "without normal-termination marker" in str(guided.get("reason") or "")
                )
                current_result_is_startup_failure = (
                    current_result_is_startup_failure or legacy_startup_misclassification
                )
                budget = route_result.get("surface_resolution_budget") or {}
                for item in (
                    budget.get("retained_process_leases") or ()
                    if current_result_is_startup_failure else ()
                ):
                    lease_path = Path(str(item.get("lease") or ""))
                    if not lease_path.is_absolute():
                        lease_path = workspace / lease_path
                    if not lease_path.is_file():
                        continue
                    lease = json.loads(lease_path.read_text())
                    raw_output = lease.get("output")
                    if raw_output:
                        retained_outputs.add(Path(raw_output))
        if not retained_outputs and iteration.get("path_verification_status") \
                == "stationary_point_validation_failed":
            retained_outputs.update(workspace.glob("path/stationary-points/*/*/opt.out"))
        for output in sorted(retained_outputs):
            if not output.is_file():
                continue
            text = output.read_text(errors="replace")
            if any(marker in text for marker in markers):
                failed_outputs.append(str(output))
        if failed_outputs:
            eligible.append(route_id)
            outputs.extend(failed_outputs)
    if not eligible:
        return None
    return {
        "kind": "orca_startup_infrastructure_execution_failure",
        "status": "recoverable_execution_failure",
        "route_ids": eligible,
        "failed_outputs": outputs,
        "scientific_path_failure": False,
        "reason": "ORCA/OpenMPI startup was denied before electronic-structure execution",
    }


def _recoverable_budget_serialization_failure(state_path: Path, state: dict) -> dict | None:
    """Recognize the pre-calculation tuple/list policy comparison defect."""
    expected = "Retained surface-resolution budget policy hash was mutated"
    eligible = []
    failures = []
    for iteration in state.get("iterations") or ():
        artifacts = iteration.get("artifacts") or {}
        if iteration.get("status") != "failed" \
                or iteration.get("scientific_failure_state") != "verification_iteration_failed" \
                or expected not in str(iteration.get("failure_reason") or "") \
                or artifacts.get("path_result"):
            continue
        raw_failure = artifacts.get("failure")
        failure = _resolve_artifact_path(raw_failure, state_path.parent.parent)
        if not failure.is_file():
            continue
        retained = json.loads(failure.read_text())
        if retained.get("error_type") != "RuntimeError" \
                or retained.get("message") != expected:
            continue
        eligible.append(str(iteration["route_id"]))
        failures.append(str(failure))
    if not eligible:
        return None
    return {
        "kind": "surface_budget_policy_serialization_execution_failure",
        "status": "recoverable_precalculation_software_failure",
        "route_ids": eligible,
        "failure_artifacts": failures,
        "scientific_path_failure": False,
        "calculations_started": False,
        "reason": "Canonical JSON policy equals the retained policy; tuple/list comparison was invalid",
    }


def _recoverable_surface_context_amendment_failure(
    state_path: Path, state: dict,
) -> dict | None:
    """Recognize the retired v4 context-comparison defect exactly.

    This is deliberately not a generic exception retry.  Eligibility requires
    the old message, no returned path result, and a retained v2/v3 -> current
    scientific-policy amendment governing the failed iteration.
    """
    expected = "Retained surface-resolution scientific context differs from this run"
    continuation = state.get("continuation_evidence") or {}
    if (
        (continuation.get("execution_recovery_evidence") or {}).get("kind")
        == "surface_context_policy_amendment_execution_failure"
    ):
        # One repaired execution segment is the bound.  Seeing the retired
        # defect again after that segment indicates a new defect, not license
        # for an unbounded same-policy retry loop.
        return None
    amendment = continuation.get("scientific_policy_amendment") or {}
    current_policy = str(state.get("path_verifier_policy_id") or "")
    if (
        current_policy != ENDPOINT_RETURN_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID
        or amendment.get("new_path_verifier_policy_id") != current_policy
        or amendment.get("prior_path_verifier_policy_id") not in {
            LEGACY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
            PREVIOUS_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        }
    ):
        return None
    eligible: list[str] = []
    failures: list[str] = []
    for iteration in state.get("iterations") or ():
        artifacts = iteration.get("artifacts") or {}
        if (
            iteration.get("status") != "failed"
            or iteration.get("scientific_failure_state") != "verification_iteration_failed"
            or iteration.get("path_verifier_policy_id") != current_policy
            or str(iteration.get("failure_reason") or "") != f"RuntimeError: {expected}"
            or artifacts.get("path_result")
        ):
            continue
        failure = _resolve_artifact_path(
            artifacts.get("failure"), state_path.parent.parent,
        )
        if not failure.is_file():
            continue
        retained = json.loads(failure.read_text())
        if (
            retained.get("error_type") != "RuntimeError"
            or retained.get("message") != expected
        ):
            continue
        eligible.append(str(iteration["route_id"]))
        failures.append(str(failure))
    if not eligible:
        return None
    return {
        "kind": "surface_context_policy_amendment_execution_failure",
        "status": "recoverable_precalculation_software_failure",
        "route_ids": eligible,
        "failure_artifacts": failures,
        "scientific_path_failure": False,
        "calculations_started_in_failed_segment": False,
        "reason": (
            "The retired context comparator rejected a valid, retained path-policy "
            "amendment before starting a new electronic-structure calculation"
        ),
    }


def _recoverable_surface_attempt_context_identity_failure(
    state_path: Path, state: dict,
) -> dict | None:
    """Recognize the retired v6 attempt-ID/context omission exactly once."""
    expected = "Retained surface-resolution attempt evidence was mutated"
    recovery_kind = "surface_attempt_policy_context_identity_execution_failure"
    continuation = state.get("continuation_evidence") or {}
    if (
        (continuation.get("execution_recovery_evidence") or {}).get("kind")
        == recovery_kind
    ):
        return None
    amendment = continuation.get("scientific_policy_amendment") or {}
    current_policy = str(state.get("path_verifier_policy_id") or "")
    if (
        current_policy != GENERIC_ORCA_PATH_VERIFIER_POLICY_ID
        or amendment.get("new_path_verifier_policy_id") != current_policy
        or amendment.get("prior_path_verifier_policy_id")
        != DISSOCIATION_ESCALATION_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID
    ):
        return None
    eligible: list[str] = []
    failures: list[str] = []
    for iteration in state.get("iterations") or ():
        artifacts = iteration.get("artifacts") or {}
        if (
            iteration.get("status") != "failed"
            or iteration.get("scientific_failure_state") != "verification_iteration_failed"
            or iteration.get("path_verifier_policy_id") != current_policy
            or str(iteration.get("failure_reason") or "") != f"RuntimeError: {expected}"
            or artifacts.get("path_result")
        ):
            continue
        failure = _resolve_artifact_path(
            artifacts.get("failure"), state_path.parent.parent,
        )
        if not failure.is_file():
            continue
        retained = json.loads(failure.read_text())
        if (
            retained.get("error_type") != "RuntimeError"
            or retained.get("message") != expected
        ):
            continue
        eligible.append(str(iteration["route_id"]))
        failures.append(str(failure))
    if not eligible:
        return None
    return {
        "kind": recovery_kind,
        "status": "recoverable_precalculation_software_failure",
        "route_ids": eligible,
        "failure_artifacts": failures,
        "scientific_path_failure": False,
        "calculations_started_in_failed_segment": False,
        "reason": (
            "The retired attempt identifier omitted the append-only scientific "
            "context hash and falsely classified a v5-to-v6 evidence successor as mutation"
        ),
    }


def _recoverable_execution_failure(state_path: Path, state: dict) -> dict | None:
    return (
        _recoverable_neb_preopt_parser_failure(state_path, state)
        or _recoverable_neb_restart_preopt_parser_failure(state_path, state)
        or _recoverable_neb_time_limit_checkpoint(state_path, state)
        or _recoverable_orca_startup_failure(state_path, state)
        or _recoverable_budget_serialization_failure(state_path, state)
        or _recoverable_surface_context_amendment_failure(state_path, state)
        or _recoverable_surface_attempt_context_identity_failure(state_path, state)
    )


def _recoverable_neb_time_limit_checkpoint(state_path: Path, state: dict) -> dict | None:
    """Resume only a distinct native band from an owned, terminal time limit."""
    from .reaction_path import inspect_neb_time_limit_checkpoint, NEB_TIME_LIMIT_CHECKPOINT_POLICY

    prior = ((state.get("continuation_evidence") or {}).get("execution_recovery_evidence") or {})
    def segment_key(item, route_id=None):
        source = item.get("source_input") or {}
        return (str(item.get("route_id") or route_id or ""),
                str(Path(source.get("path", "")).resolve()),
                (item.get("source_job_contract") or {}).get("sha256"))
    already_used = {segment_key(item)
                    for item in prior.get("failures") or ()} if prior.get("kind") == NEB_TIME_LIMIT_CHECKPOINT_POLICY else set()
    latest = {str(item.get("route_id")): item for item in state.get("iterations") or ()}
    route_ids, failures = [], []
    for route_id, iteration in latest.items():
        if iteration.get("status") != "deferred_unverified":
            continue
        raw_path = (iteration.get("artifacts") or {}).get("path_result")
        if not raw_path:
            continue
        path = _resolve_artifact_path(raw_path, state_path.parent.parent)
        if not path.is_file():
            continue
        result = json.loads(path.read_text())
        for orientation in result.get("orientations") or ():
            neb = orientation.get("neb") or {}
            match = re.fullmatch(r"Command '(.*)' timed out after ([0-9.]+) seconds", str(neb.get("failure_reason", "")))
            if neb.get("status") != "failed" or not neb.get("input") or not match:
                continue
            try:
                command = ast.literal_eval(match.group(1))
                inp = Path(neb["input"])
                if (not isinstance(command, list) or len(command) != 2
                        or not all(isinstance(arg, str) for arg in command)
                        or command[1] != inp.name or Path(command[0]).name != "orca"):
                    continue
                evidence = inspect_neb_time_limit_checkpoint(inp, timeout_seconds=float(match.group(2)))
            except (ValueError, TypeError, SyntaxError):
                continue
            if evidence is not None and segment_key(evidence, route_id) not in already_used:
                route_ids.append(route_id)
                failures.append({"route_id": route_id, "path_result_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                 "path_result": str(path), **evidence})
                break
    return ({"kind": NEB_TIME_LIMIT_CHECKPOINT_POLICY, "status": "restartable_stopped_time_limit",
             "route_ids": route_ids, "failures": failures, "rate_claim_allowed": False,
             "scientific_path_failure": False, "renews_compute_budget": False,
             "reason": "A terminal owned NEB time limit retained a finite native band; the original budget still applies"}
            if route_ids else None)


def _recoverable_neb_preopt_parser_failure(state_path: Path, state: dict) -> dict | None:
    """Recover the exact generated-keyword defect once, not a chemical failure."""
    from src.orca_input_repair import inspect_neb_preopt_parser_failure, NEB_PREOPT_REPAIR_POLICY

    if ((state.get("continuation_evidence") or {}).get("execution_recovery_evidence") or {}).get(
        "kind"
    ) == NEB_PREOPT_REPAIR_POLICY:
        return None
    latest = {str(item.get("route_id")): item for item in state.get("iterations") or ()}
    route_ids, failures = [], []
    for route_id, iteration in latest.items():
        if iteration.get("status") != "deferred_unverified":
            continue
        raw_path = (iteration.get("artifacts") or {}).get("path_result")
        if not raw_path:
            continue
        path = _resolve_artifact_path(raw_path, state_path.parent.parent)
        if not path.is_file():
            continue
        result = json.loads(path.read_text())
        for orientation in result.get("orientations") or ():
            neb = orientation.get("neb") or {}
            if neb.get("status") != "failed" or not neb.get("input"):
                continue
            evidence = inspect_neb_preopt_parser_failure(Path(neb["input"]))
            if evidence is not None:
                route_ids.append(route_id)
                failures.append({"path_result_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                 "path_result": str(path), **evidence})
                break
    return ({"kind": NEB_PREOPT_REPAIR_POLICY, "status": "recoverable_precalculation_software_failure",
             "route_ids": route_ids, "failures": failures, "scientific_path_failure": False,
             "reason": "ORCA 6.1 rejected the generated PreOptEnds keyword before any NEB calculation"}
            if route_ids else None)


def _recoverable_neb_restart_preopt_parser_failure(state_path: Path, state: dict) -> dict | None:
    """Recover only the observed stopped ORCA restart/preoptimization conflict."""
    from src.orca_input_repair import inspect_neb_restart_preopt_parser_failure, NEB_RESTART_PREOPT_REPAIR_POLICY

    prior = ((state.get("continuation_evidence") or {}).get("execution_recovery_evidence") or {})
    latest = {str(item.get("route_id")): item for item in state.get("iterations") or ()}
    route_ids, failures = [], []
    for route_id, iteration in latest.items():
        if iteration.get("status") != "deferred_unverified":
            continue
        raw_path = (iteration.get("artifacts") or {}).get("path_result")
        if not raw_path:
            continue
        path = _resolve_artifact_path(raw_path, state_path.parent.parent)
        if not path.is_file():
            continue
        result = json.loads(path.read_text())
        for orientation in result.get("orientations") or ():
            neb = orientation.get("neb") or {}
            if neb.get("status") != "failed" or not neb.get("input"):
                continue
            evidence = inspect_neb_restart_preopt_parser_failure(Path(neb["input"]))
            if evidence is None:
                continue
            if prior.get("kind") == NEB_RESTART_PREOPT_REPAIR_POLICY and any(
                    item.get("route_id") == route_id and item.get("source_input") == evidence["source_input"]
                    and item.get("source_job_contract") == evidence["source_job_contract"]
                    for item in prior.get("failures") or ()):
                continue
            route_ids.append(route_id)
            failures.append({"route_id": route_id, "path_result_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                             "path_result": str(path), **evidence})
            break
    return ({"kind": NEB_RESTART_PREOPT_REPAIR_POLICY, "status": "recoverable_precalculation_software_failure",
             "route_ids": route_ids, "failures": failures, "scientific_path_failure": False,
             "rate_claim_allowed": False, "renews_compute_budget": False,
             "reason": "ORCA 6.1 rejected endpoint preoptimization combined with its native restart band before new NEB iterations"}
            if route_ids else None)


def write_verification_session_reconciliation(run_dir: Path) -> dict:
    """Append an operational interpretation of immutable engine-state records.

    Historical state files are evidence and are never rewritten.  A retained
    ``running`` record with a hashed continuation descendant is therefore
    reported as superseded, while only the leaf session is assessed for a live
    or interrupted execution.
    """
    run_dir = Path(run_dir)
    state_paths = sorted(run_dir.glob("verification-engine*/engine-state.json"))
    records: list[tuple[Path, dict, str]] = []
    for state_path in state_paths:
        raw = state_path.read_bytes()
        state = json.loads(raw)
        records.append((state_path, state, hashlib.sha256(raw).hexdigest()))
    hashes = {digest: path for path, _, digest in records}
    descendants: dict[Path, list[Path]] = {path: [] for path, _, _ in records}
    for path, state, _ in records:
        predecessor_hash = str(
            (state.get("continuation_evidence") or {}).get("predecessor_state_sha256")
            or ""
        )
        predecessor = hashes.get(predecessor_hash)
        if predecessor is not None:
            descendants[predecessor].append(path)

    sessions: list[dict[str, Any]] = []
    for path, state, digest in records:
        successors = sorted(descendants[path])
        raw_status = str(state.get("status") or "unknown")
        heartbeat = None
        if successors:
            effective_status = "superseded_by_continuation"
            reason = "A hash-linked continuation retains the effective later execution state."
        elif raw_status == "running":
            heartbeat = _started_execution_heartbeat(path)
            if heartbeat.get("status") == "active":
                effective_status = "active"
                reason = "The leaf session retains live or indeterminate execution evidence."
            else:
                effective_status = "interrupted_execution_recoverable"
                reason = "The leaf session is marked running but retained process evidence is stopped."
        else:
            effective_status = raw_status
            reason = "The immutable leaf session already retains a terminal state."
        sessions.append({
            "session": path.parent.name,
            "state_path": str(path),
            "state_sha256": digest,
            "raw_status": raw_status,
            "effective_status": effective_status,
            "reason": reason,
            "successor_state_paths": [str(item) for item in successors],
            "heartbeat": heartbeat,
        })
    leaves = [item for item in sessions if not item["successor_state_paths"]]
    snapshot = {
        "schema_version": 1,
        "kind": "verification_session_reconciliation",
        "status": (
            "complete" if len(leaves) == 1 else
            "empty" if not sessions else "ambiguous_multiple_leaves"
        ),
        "immutable_state_files_modified": False,
        "session_count": len(sessions),
        "leaf_sessions": [item["session"] for item in leaves],
        "sessions": sessions,
        "interpretation": (
            "This artifact reconciles operational session state only. It does not alter "
            "reaction, rate, convergence, or lifetime evidence."
        ),
    }
    folder = run_dir / SESSION_RECONCILIATION_DIRECTORY
    folder.mkdir(parents=True, exist_ok=True)
    existing = sorted(folder.glob("*.json"))
    for path in existing:
        retained = json.loads(path.read_text())
        comparable = {
            key: value for key, value in retained.items()
            if key not in {"revision", "created_at", "predecessor_revision_sha256", "revision_sha256"}
        }
        if comparable == snapshot:
            return {**retained, "path": str(path)}
    body = {
        **snapshot,
        "revision": len(existing) + 1,
        "created_at": datetime.now(UTC).isoformat(),
        "predecessor_revision_sha256": (
            json.loads(existing[-1].read_text()).get("revision_sha256") if existing else None
        ),
    }
    body["revision_sha256"] = _canonical_hash(body)
    path = folder / f"{len(existing) + 1:06d}.json"
    _write_json_atomic(path, body)
    return {**body, "path": str(path)}


def _apply_verification_summary(result: dict, summary: dict) -> None:
    previous = result.get("verification_summary")
    if isinstance(previous, dict):
        result.setdefault("verification_history", []).append({
            "recorded_at": datetime.now(UTC).isoformat(),
            "status": previous.get("status"),
            "reason": previous.get("reason"),
            "artifacts": previous.get("artifacts"),
        })
    result["verification_summary"] = summary
    status = summary.get("status")
    if status == "verified_t95_converged":
        result["assessment"] = {
            "status": "orca_verified_condition_dependent_instability",
            "reason": (
                "Flux-controlling ORCA/Arkane rates were exactly reinjected and the repaired "
                "condition model converged to the retention threshold."
            ),
        }
        result.setdefault("kinetic_lifetime", {})["estimated_time_to_retention_seconds"] = \
            summary.get("estimated_time_to_retention_seconds")
    elif status == "verified_below_loss_threshold":
        result["assessment"] = {
            "status": "verified_route_below_loss_threshold",
            "reason": (
                "Verified controlling-route rates and conservative unresolved-route bounds "
                "remain below the retention threshold in this condition model."
            ),
        }
    elif status == "no_target_loss_in_completed_rmg_model":
        result["assessment"] = {"status": "no_loss_detected_in_rmg_model", "reason": summary.get("reason")}
    elif status == "conditional_model_prediction":
        result["assessment"] = {
            "status": "conditional_kinetic_model_prediction",
            "reason": summary.get("reason") or (
                "A condition-bound fixed kinetic model was propagated, but pressure "
                "applicability, rate uncertainty, and mechanism coverage remain unresolved."
            ),
        }
        result["conditional_model_prediction"] = summary.get("conditional_model_prediction")
        result.setdefault("kinetic_lifetime", {})["estimated_time_to_retention_seconds"] = None
    else:
        result["assessment"] = {
            "status": "orca_verification_incomplete",
            "reason": summary.get("reason") or "Continued verification did not converge.",
        }


def export_specialist_rate_requests(run_dir: Path) -> dict:
    """Export exact requests for retained paths that need specialist kinetics."""
    from .specialist_rates import build_specialist_rate_request

    run_dir = Path(run_dir).resolve()
    result_path = run_dir / "stability.json"
    if not result_path.is_file():
        raise FileNotFoundError(f"Retained stability result not found: {result_path}")
    result = json.loads(result_path.read_text())
    summary = result.get("verification_summary") or {}
    if summary.get("status") != "verification_incomplete":
        raise ValueError("Specialist requests require an incomplete retained verification")
    state_raw = (summary.get("artifacts") or {}).get("engine_state")
    if not state_raw:
        raise ValueError("Incomplete verification does not retain an engine state")
    state_path = _resolve_artifact_path(state_raw, run_dir)
    if not state_path.is_file():
        raise FileNotFoundError(f"Retained verification state is missing: {state_path}")
    state = json.loads(state_path.read_text())
    condition = state.get("condition_contract") or {}
    if state.get("condition_contract_sha256") != hashlib.sha256(
        json.dumps(condition, sort_keys=True, separators=(",", ":"), default=str).encode()
    ).hexdigest():
        raise ValueError("Retained verification condition contract has been mutated")

    eligible = {
        route_id: item
        for route_id, item in (state.get("deferred_unverified_routes") or {}).items()
        if str((item or {}).get("status")) in {
            "barrierless_rate_model_required", "dissociation_rate_model_required",
        }
    }
    if not eligible:
        raise ValueError("No retained route currently requires a specialist rate model")
    exported = []
    root = run_dir / "specialist-rate-requests"
    for route_id, deferred in sorted(eligible.items()):
        matching = [
            item for item in state.get("iterations") or ()
            if str(item.get("route_id")) == str(route_id)
        ]
        if not matching:
            raise ValueError(f"Deferred route {route_id} has no retained iteration")
        artifacts = matching[-1].get("artifacts") or {}
        route_path = _resolve_artifact_path(artifacts.get("route"), run_dir)
        path_result_path = _resolve_artifact_path(
            deferred.get("path_result") or artifacts.get("path_result"), run_dir,
        )
        if not route_path.is_file() or not path_result_path.is_file():
            raise FileNotFoundError(
                f"Specialist request source artifacts are missing for route {route_id}"
            )
        route = json.loads(route_path.read_text())
        path_result = json.loads(path_result_path.read_text())
        request = build_specialist_rate_request(route, path_result, condition)
        safe_route = re.sub(r"[^A-Za-z0-9_.-]+", "-", str(route_id)).strip("-") or "route"
        request_dir = root / safe_route / str(request["request_sha256"]).removeprefix("sha256:")
        request_dir.mkdir(parents=True, exist_ok=True)
        request_path = request_dir / "request.json"
        if request_path.exists() and json.loads(request_path.read_text()) != request:
            raise ValueError(f"Existing exported specialist request was mutated: {request_path}")
        _write_json_atomic(request_path, request)
        source = {
            "schema_version": 1,
            "kind": "storca_specialist_rate_request_export",
            "route_id": route_id,
            "request": str(request_path),
            "request_sha256": request["request_sha256"],
            "source_engine_state": str(state_path),
            "source_engine_state_sha256": _evidence_input_hash(state_path),
            "source_route": str(route_path),
            "source_route_sha256": _evidence_input_hash(route_path),
            "source_path_result": str(path_result_path),
            "source_path_result_sha256": _evidence_input_hash(path_result_path),
            "allowed_rate_theories": request["allowed_rate_theories"],
            "ordinary_harmonic_tst_authorized": False,
        }
        source["export_sha256"] = _canonical_hash(source)
        manifest_path = _write_json_atomic(request_dir / "export.json", source)
        exported.append({
            "route_id": route_id,
            "request": str(request_path),
            "manifest": str(manifest_path),
            "request_sha256": request["request_sha256"],
            "allowed_rate_theories": request["allowed_rate_theories"],
        })
    return {
        "schema_version": 1,
        "kind": "storca_specialist_rate_request_set",
        "status": "requests_exported",
        "run_directory": str(run_dir),
        "requests": exported,
        "next_action": (
            "Run an applicable external specialist model, create a schema-v1 "
            "storca_specialist_rate_result with immutable evidence and a checksummed RMG "
            "library, then resume with that result file or result directory."
        ),
    }


def resume_stability_verification(
    run_dir: Path,
    *,
    ncores: int | None = None,
    specialist_rate_result: Path | None = None,
) -> dict:
    """Resume deferred route verification without accepting scientific overrides."""
    from .stability import collect_rmg_evidence
    from .verification_adapters import make_arkane_rate_builder, make_orca_route_verifier
    from .specialist_rates import RetainedSpecialistRateProvider

    run_dir = Path(run_dir)
    result_path = run_dir / "stability.json"
    metadata_path = run_dir / "metadata.json"
    if not result_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError("A retained stability.json and metadata.json are required")
    result = json.loads(result_path.read_text())
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("workflow") != "orca_rmg_stability_screen":
        raise ValueError("Run is not a gas stability-screen workflow")
    if (result.get("verification_summary") or {}).get("status") != "verification_incomplete":
        raise ValueError("Only an incomplete retained verification can be continued")
    retained_contract = _load_execution_contract(run_dir, result, metadata)
    resource_amendment = None
    if ncores is not None:
        resource_amendment = write_execution_resource_amendment(
            run_dir, retained_contract, ncores=int(ncores),
        )
    contract, resource_amendments = _effective_execution_contract(run_dir, retained_contract)
    if resource_amendments:
        resource_amendment = resource_amendments[-1]
    effective_path_policy, path_policy_amendments = _effective_path_policy(
        run_dir, retained_contract,
    )
    path_policy_amendment = path_policy_amendments[-1] if path_policy_amendments else None
    retained_policy_artifact = (
        path_policy_amendment.get("new_path_verifier_policy")
        if path_policy_amendment else retained_contract.get("path_verifier_policy")
    )
    if effective_path_policy == GENERIC_ORCA_PATH_VERIFIER_POLICY_ID \
            and retained_policy_artifact:
        current_policy = SurfaceResolutionPolicy().as_dict()
        if (
            retained_policy_artifact.get("scientific_hash")
            != current_policy["scientific_hash"]
            or retained_policy_artifact != current_policy
        ):
            raise ValueError(
                "Retained path-policy definition differs from the current implementation; "
                "a new scientific policy ID is required before continuation"
            )
    # An existing scientific-policy chain is already part of the retained
    # history and must govern descendant discovery.  A first migration is
    # still deferred until the original-policy chain is audited below.
    if path_policy_amendments:
        contract["path_verifier_policy_id"] = effective_path_policy
    condition_contract = contract["condition_contract"]
    if _canonical_hash(condition_contract) != contract.get("condition_contract_sha256"):
        raise ValueError("Execution contract condition hash does not match its payload")

    predecessor_raw = ((result.get("verification_summary") or {}).get("artifacts") or {}).get("engine_state")
    if not predecessor_raw:
        raise ValueError("Incomplete verification summary does not retain an engine state")
    predecessor = _resolve_artifact_path(predecessor_raw, run_dir)
    if not predecessor.is_file():
        raise FileNotFoundError(f"Retained verification engine state is missing: {predecessor}")
    state = json.loads(predecessor.read_text())
    retained_summary = state.get("verification_summary") or {}
    current_summary = result.get("verification_summary") or {}
    if retained_summary.get("status") != current_summary.get("status") \
            or retained_summary.get("reason") != current_summary.get("reason"):
        raise ValueError("stability.json and predecessor engine terminal summaries disagree")
    state_deferred = state.get("deferred_unverified_routes") or {}
    summary_deferred = {
        str(item.get("route_id")): str(item.get("status"))
        for item in current_summary.get("deferred_unverified_routes") or ()
    }
    if summary_deferred != {
        str(route_id): str(item.get("status"))
        for route_id, item in state_deferred.items()
    }:
        raise ValueError("stability.json and predecessor deferred-route evidence disagree")

    # Follow the append-only continuation chain.  A detached shell is not
    # evidence that its ORCA child died: fresh artifacts make the session
    # active and categorically prohibit a duplicate process.  Only a quiet,
    # retained, started execution may reuse the same scientific policy attempt.
    recovery_route_ids: tuple[str, ...] = ()
    recovery_heartbeat: dict | None = None
    continuation_predecessor = predecessor
    predecessor_infrastructure_failure = _recoverable_execution_failure(predecessor, state)
    if predecessor_infrastructure_failure:
        recovery_heartbeat = predecessor_infrastructure_failure
        recovery_route_ids = tuple(predecessor_infrastructure_failure["route_ids"])
    seen_states = {predecessor.resolve()}
    while True:
        continuation_digest = hashlib.sha256(continuation_predecessor.read_bytes()).hexdigest()
        descendants = []
        for candidate in sorted(run_dir.glob("verification-engine*/engine-state.json")):
            if candidate.resolve() in seen_states:
                continue
            candidate_state = json.loads(candidate.read_text())
            continuation = candidate_state.get("continuation_evidence") or {}
            if continuation.get("predecessor_state_sha256") == continuation_digest \
                    and continuation.get("path_verifier_policy_id") == contract["path_verifier_policy_id"]:
                descendants.append((candidate, candidate_state))
        if not descendants:
            break
        if len(descendants) > 1:
            raise ValueError(
                "Ambiguous verification history: multiple continuations share predecessor "
                f"{continuation_predecessor}"
            )
        candidate, candidate_state = descendants[0]
        seen_states.add(candidate.resolve())
        infrastructure_failure = _recoverable_execution_failure(candidate, candidate_state)
        if infrastructure_failure:
            continuation_predecessor = candidate
            recovery_heartbeat = infrastructure_failure
            recovery_route_ids = tuple(infrastructure_failure["route_ids"])
            continue
        interrupted_terminal = (
            candidate_state.get("status") == "failed_closed"
            and (candidate_state.get("terminal_error") or {}).get("error_type")
            in {"KeyboardInterrupt", "InterruptedError", "SystemExit"}
        )
        if candidate_state.get("status") != "running" and not interrupted_terminal:
            raise ValueError(
                "This scientific continuation is already retained and is not an interrupted "
                f"running session: {candidate.parent} (status={candidate_state.get('status')})"
            )
        heartbeat = _started_execution_heartbeat(candidate)
        if heartbeat["status"] == "active":
            raise ValueError(
                "A verification continuation is still active; refusing to launch a duplicate. "
                f"Newest artifact {heartbeat['newest_artifact']} is "
                f"{heartbeat['artifact_age_seconds']:.1f} seconds old"
            )
        continuation_predecessor = candidate
        recovery_heartbeat = heartbeat
        recovery_route_ids = tuple(heartbeat["route_ids"])

    # Scientific policy changes occur only after the retained-policy chain has
    # been audited.  An interrupted execution must be recovered under the same
    # policy; a terminal v2/v3 deferral may enter the current policy through a
    # hashed amendment.
    if not recovery_route_ids and specialist_rate_result is None \
            and effective_path_policy in {
                LEGACY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                PREVIOUS_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                ENDPOINT_RETURN_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                DISSOCIATION_ESCALATION_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                PRE_IRC_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                PRE_RING_CONNECTIVITY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                PRE_GEOMETRY_AUDIT_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                PRE_REACHED_MINIMA_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
            }:
        path_policy_amendment = write_path_policy_amendment(
            run_dir,
            retained_contract,
            path_verifier_policy_id=GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
        )
        effective_path_policy, path_policy_amendments = _effective_path_policy(
            run_dir, retained_contract,
        )
    contract["path_verifier_policy_id"] = effective_path_policy

    rmg = contract["rmg"]
    resources = rmg["resource_config"]
    scenario_config = rmg["scenario_config"]
    conditions = build_condition_spec(
        scenario_config,
        temperature_K=float(condition_contract["temperature_K"]),
        pressure_bar=float(condition_contract["pressure_bar"]),
        target_duration_hours=float(condition_contract["target_duration_seconds"]) / 3600.0,
        retention_fraction=float(condition_contract["retention_fraction"]),
        light_condition=str(condition_contract["light_condition"]),
        light_model=condition_contract.get("light_model"),
        relative_humidity=condition_contract.get("relative_humidity"),
    )
    if conditions.as_dict() != condition_contract:
        raise ValueError("Retained scenario cannot reproduce the immutable condition contract")

    base_libraries = [Path(item["path"]) for item in rmg.get("base_reaction_libraries") or ()]
    verify_route = make_orca_route_verifier(
        contract["parent"]["smiles"], ncores=int(contract["resources"]["ncores"]),
        method_profile=contract["method_profile"], condition_contract=condition_contract,
        timeout_seconds=contract["resources"]["route_job_timeout_seconds"],
        orientations=int(contract["resources"]["orientations"]),
        neb_images=int(contract["resources"]["neb_images"]),
        rmg_env=rmg.get("environment"),
    )
    build_rate = make_arkane_rate_builder(
        condition_contract, method_profile=contract["method_profile"],
        rmg_env=rmg.get("environment"),
    )
    specialist_provider = (
        RetainedSpecialistRateProvider(Path(specialist_rate_result))
        if specialist_rate_result is not None else None
    )

    def rerun_model(library_paths: list[Path], verification_dir: Path) -> dict:
        libraries: list[Path] = []
        for library in [*library_paths, *base_libraries]:
            resolved = Path(library)
            if resolved not in libraries:
                libraries.append(resolved)
        repaired = collect_rmg_evidence(
            contract["parent"]["smiles"], verification_dir,
            rmg_env=rmg.get("environment"), rmg_command=rmg["command"],
            barrier_threshold=float(rmg["barrier_threshold_kcal_mol"]),
            temperature=float(condition_contract["temperature_K"]),
            pressure=float(condition_contract["pressure_bar"]),
            rmg_walltime=resources["walltime"],
            rmg_max_processes=resources["max_processes"],
            rmg_max_iterations=resources["max_iterations"],
            rmg_max_edge_species=resources["max_edge_species"],
            rmg_max_objects_per_iteration=resources["max_objects_per_iteration"],
            rmg_maximum_heavy_atoms=rmg.get("maximum_heavy_atoms"),
            scenario=scenario_config, conditions=conditions,
            reaction_libraries=libraries,
            database_reaction_libraries=rmg.get("database_reaction_libraries") or (),
            rmg_retry_attempts=int(rmg["retry_attempts"]),
            rmg_hard_timeout_seconds=rmg.get("hard_timeout_seconds"),
            rmg_heartbeat_timeout_seconds=rmg.get("heartbeat_timeout_seconds"),
            rmg_start_method=rmg["start_method"],
        )
        return {"rmg_evidence": repaired, "condition_contract": condition_contract}

    retained_deferred = state.get("deferred_unverified_routes") or {}
    if specialist_provider is not None:
        deferred_ids = tuple(
            route_id for route_id, item in retained_deferred.items()
            if str((item or {}).get("status")) in {
                "barrierless_rate_model_required", "dissociation_rate_model_required",
            }
        )
        if not deferred_ids:
            raise ValueError(
                "A specialist result was supplied, but no retained route requires a specialist rate"
            )
    else:
        deferred_ids = tuple(retained_deferred.keys())
    continuation_kwargs = (
        {
            "recover_in_progress_route_ids": recovery_route_ids,
            "execution_recovery_evidence": {
                **(recovery_heartbeat or {}),
                "execution_resource_amendment": resource_amendment,
                "effective_ncores": int(contract["resources"]["ncores"]),
            },
        }
        if recovery_route_ids else
        {"retry_deferred_route_ids": deferred_ids}
    )
    summary = continue_verification_engine(
        result["rmg_evidence"], condition_contract, run_dir, continuation_predecessor,
        verify_route=verify_route, build_rate=build_rate, rerun_model=rerun_model,
        path_verifier_policy_id=contract["path_verifier_policy_id"],
        scientific_policy_amendment=path_policy_amendment,
        config=VerificationEngineConfig(
            max_iterations=int(contract["verification_engine_config"]["max_iterations"]),
        ),
        build_specialist_rate=(specialist_provider.evaluate if specialist_provider else None),
        **continuation_kwargs,
    )
    _apply_verification_summary(result, summary)
    reconciliation = write_verification_session_reconciliation(run_dir)
    result["verification_session_reconciliation"] = {
        "status": reconciliation["status"],
        "path": reconciliation["path"],
        "leaf_sessions": reconciliation["leaf_sessions"],
    }
    _write_json_atomic(result_path, result)
    assessment = result["assessment"]["status"]
    terminal = assessment in {
        "no_loss_detected_in_rmg_model", "verified_route_below_loss_threshold",
        "orca_verified_condition_dependent_instability",
    }
    write_metadata(
        run_dir, result_json=str(result_path), assessment_status=assessment,
        status="completed" if terminal else "completed_with_risk_or_incomplete_evidence",
        verification_resume_count=len(result.get("verification_history") or ()),
        latest_verification_session=(summary.get("artifacts") or {}).get("session_directory"),
        verification_session_reconciliation=reconciliation["path"],
        verification_effective_ncores=int(contract["resources"]["ncores"]),
        verification_execution_resource_amendment=(
            resource_amendment.get("amendment_sha256") if isinstance(resource_amendment, dict) else None
        ),
        verification_path_policy_amendment=(
            path_policy_amendment.get("amendment_sha256")
            if isinstance(path_policy_amendment, dict) else None
        ),
    )
    return {**result, "result_json": result_path}


__all__ = [
    "EXECUTION_CONTRACT_NAME", "PATH_POLICY_AMENDMENT_PREFIX",
    "build_verification_execution_contract",
    "export_specialist_rate_requests",
    "resume_stability_verification", "write_execution_resource_amendment",
    "write_path_policy_amendment",
    "write_verification_session_reconciliation",
    "write_verification_execution_contract",
]
