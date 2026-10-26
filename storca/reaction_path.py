"""Resumable ORCA decomposition-path calculations and explanation artifacts."""

from __future__ import annotations

import json
import hashlib
import math
import os
from pathlib import Path
import re
import shutil
import time
from typing import Any, Callable, Iterable

import numpy as np

from src.inputgen import (
    create_orca_guided_scan_input,
    create_orca_input,
    create_orca_irc_input,
    create_orca_neb_ts_input,
    create_orca_optimization_recovery_input,
    create_orca_relaxed_scan_input,
    normalize_thermochemistry_temperature,
    orca_frequency_temperature_input_evidence,
)
from src.molecule_tools import smiles_to_xyz
from src.orca_runner import run_orca
from src.orca_job_settings import active_native_settings, apply_native_job_settings
from src.stability.freq_check import frequency_stability_check

from .decomposition_visuals import (
    render_candidate_storyboard,
    render_energy_profile,
    render_xyz_trajectory_gif,
    write_visual_manifest,
    xyz_frames,
)
from .ir_modes import parse_orca_normal_modes
from .intermediate_identity import resolve_validated_intermediate_identity, propose_unexpected_irc_candidate
from .reaction_candidates import ReactionCandidate
from .route_geometry import (
    build_endpoint_complex_seeds,
    aligned_rmsd,
    build_dissociation_endpoint,
    build_mapped_endpoint_seeds,
    interpolate_xyz,
    read_xyz,
    resolve_route_atom_mapping,
    validate_declared_connectivity,
    validate_route_mapping,
    validate_route_endpoint_connectivity,
    validate_trajectory_endpoints,
    validate_trajectory_frames,
    write_xyz,
    write_declared_fragment_seed,
)
from .route_verify import RouteSpec, prepare_explanation, route_spec_from_candidate_route
from .route_thermochemistry import (
    assemble_route_thermochemistry,
    condition_specific_forward_loss_bound,
    parse_orca_thermochemistry,
    atomic_hydrogen_identity,
    read_atomic_hydrogen_energy,
    atomic_hydrogen_thermochemistry,
)
from .surface_resolution import (
    EndpointReferenceAssessment,
    SurfaceDiagnosis,
    SurfaceResolutionAttempt,
    SurfaceResolutionPolicy,
    aggregate_surface_diagnoses,
    bond_change_mode_overlap,
    budget_status,
    build_bond_change_scan_plan,
    diagnose_surface,
    select_interior_scan_maximum,
)
from .workflow import run_optimization_and_frequency
from .rate_transfer_control import TransferControlError
from .verification_policy import (
    GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    PRE_RING_CONNECTIVITY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    PRE_GEOMETRY_AUDIT_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
    PRE_REACHED_MINIMA_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
)


def _write_json(path: Path, value: dict) -> Path:
    """Replace a JSON status file only after the new document is complete."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


class SurfaceResolutionBudgetExhausted(RuntimeError):
    """Raised before a new external calculation would exceed its route budget."""

    def __init__(self, evidence: dict[str, Any]):
        self.evidence = evidence
        super().__init__(
            "surface_resolution_budget_exhausted: "
            + ", ".join(evidence.get("exceeded") or ["requested calculation exceeds limit"])
        )


def _orca_input_ncores(path: Path) -> int:
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return 1
    match = re.search(r"(?im)^\s*nprocs\s+(\d+)\s*$", text)
    return max(1, int(match.group(1))) if match else 1


def _initialize_surface_budget(
    workdir: Path, policy: SurfaceResolutionPolicy, *, now_epoch: float | None = None,
) -> dict[str, Any]:
    """Create one immutable policy-attempt budget contract and an active pointer."""
    workdir = Path(workdir)
    policy_artifact = policy.as_dict()
    policy_hash = str(policy_artifact["scientific_hash"]).split(":", 1)[-1]
    attempt_dir = workdir / "surface-resolution-budgets" / policy_hash[:16]
    attempt_dir.mkdir(parents=True, exist_ok=True)
    contract_path = attempt_dir / "budget-contract.json"
    if contract_path.is_file():
        contract = json.loads(contract_path.read_text())
        if (
            contract.get("policy_scientific_hash") != policy_artifact["scientific_hash"]
            or contract.get("policy") != policy_artifact
        ):
            raise RuntimeError("Retained surface-resolution budget policy hash was mutated")
    else:
        contract = {
            "schema_version": 1,
            "kind": "storca_surface_resolution_budget_contract",
            "policy": policy_artifact,
            "policy_scientific_hash": policy_artifact["scientific_hash"],
            "started_at_epoch": float(now_epoch if now_epoch is not None else time.time()),
            "claims_permitted": [],
            "claims_prohibited": ["rate", "lifetime"],
        }
        _write_json(contract_path, contract)
    active = {
        "schema_version": 1,
        "kind": "storca_active_surface_resolution_budget",
        "contract": str(contract_path),
        "policy_scientific_hash": policy_artifact["scientific_hash"],
    }
    _write_json(workdir / ".surface-resolution-budget.json", active)
    return {**contract, "contract_path": str(contract_path)}


def _active_surface_budget(anchor: Path) -> tuple[Path, dict[str, Any]] | None:
    path = Path(anchor).resolve()
    directory = path if path.is_dir() else path.parent
    for candidate_root in (directory, *directory.parents):
        pointer = candidate_root / ".surface-resolution-budget.json"
        if not pointer.is_file():
            continue
        active = json.loads(pointer.read_text())
        contract_path = Path(str(active.get("contract") or ""))
        if not contract_path.is_file():
            raise RuntimeError("Active surface-resolution budget contract is missing")
        contract = json.loads(contract_path.read_text())
        if active.get("policy_scientific_hash") != contract.get("policy_scientific_hash"):
            raise RuntimeError("Active surface-resolution budget pointer does not match its contract")
        return candidate_root, contract
    return None


def _surface_budget_snapshot(
    anchor: Path,
    *,
    requested_calculations: int = 0,
    requested_ncores: int = 1,
    requested_timeout_seconds: float | None = None,
    now_epoch: float | None = None,
) -> dict[str, Any] | None:
    active = _active_surface_budget(anchor)
    if active is None:
        return None
    root, contract = active
    now = float(now_epoch if now_epoch is not None else time.time())
    started = float(contract["started_at_epoch"])
    policy_values = dict(contract["policy"])
    policy_values.pop("kind", None)
    policy_values.pop("scientific_hash", None)
    policy = SurfaceResolutionPolicy(**policy_values)
    from .rate_transfer_control import active_execution_control
    transfer_control = active_execution_control()
    leases: list[dict[str, Any]] = []
    total_core_seconds = 0.0
    lease_paths = root.rglob("*.process.json") if transfer_control is None else ()
    for lease_path in lease_paths:
        try:
            lease = json.loads(lease_path.read_text())
            lease_started = float(lease["started_at_epoch"])
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
        if lease_started + 1.0e-9 < started:
            continue
        ended = lease.get("ended_at_epoch")
        lease_end = min(now, float(ended)) if isinstance(ended, (int, float)) else now
        duration = max(0.0, lease_end - lease_started)
        input_path = Path(str(lease.get("input") or ""))
        ncores = int(lease.get("ncores") or _orca_input_ncores(input_path))
        total_core_seconds += duration * max(1, ncores)
        leases.append({
            "lease": str(lease_path),
            "input": str(input_path),
            "status": lease.get("status"),
            "ncores": max(1, ncores),
            "elapsed_seconds": duration,
        })
    calculation_count = len(leases)
    core_hours = total_core_seconds / 3600.0
    wall_seconds = max(0.0, now - started)
    if transfer_control is not None:
        # The opt-in transfer campaign has immutable attempt receipts. Idle
        # calendar time and unrelated process files cannot spend its compute.
        owned = transfer_control.snapshot()
        reservations = {
            key: value for key, value in owned["reservations"].items()
            if Path(value["input"]).resolve().is_relative_to(root.resolve())
        }
        leases = []
        total_core_seconds = 0.0
        wall_seconds = 0.0
        for key, reservation in reservations.items():
            terminal = owned["terminals"].get(key)
            duration = terminal["elapsed_seconds"] if terminal is not None else 0.0
            wall_seconds += duration
            total_core_seconds += duration * reservation["ncores"]
            leases.append({
                "input": reservation["input"], "ncores": reservation["ncores"],
                "elapsed_seconds": duration,
                "status": terminal["outcome"] if terminal is not None else "reservation_unresolved",
                "attempt_id": key, "accounting": "owned_transfer_execution_receipt",
            })
        calculation_count = len(reservations)
        core_hours = total_core_seconds / 3600.0
    gate = budget_status(
        policy,
        calculation_count=calculation_count,
        core_hours=core_hours,
        wall_seconds=wall_seconds,
    )
    exceeded = list(gate["exceeded"])
    if calculation_count + requested_calculations > policy.maximum_calculations:
        exceeded.append("maximum_calculations")
    remaining_core_seconds = max(0.0, policy.maximum_core_hours * 3600.0 - total_core_seconds)
    remaining_wall_seconds = max(0.0, policy.maximum_wall_seconds - wall_seconds)
    if requested_calculations and (remaining_core_seconds <= 0.0 or remaining_wall_seconds <= 0.0):
        exceeded.append(
            "maximum_core_hours" if remaining_core_seconds <= 0.0 else "maximum_wall_seconds"
        )
    requested_slots = max(1, int(requested_calculations))
    maximum_safe_timeout = min(
        remaining_wall_seconds / requested_slots,
        remaining_core_seconds / (requested_slots * max(1, int(requested_ncores))),
    )
    effective_timeout = (
        min(float(requested_timeout_seconds), maximum_safe_timeout)
        if requested_timeout_seconds is not None else maximum_safe_timeout
    )
    snapshot = {
        **gate,
        "status": "surface_resolution_budget_exhausted" if exceeded else "within_budget",
        "exceeded": sorted(set(exceeded)),
        "requested": {
            "calculation_count": int(requested_calculations),
            "ncores": int(requested_ncores),
            "timeout_seconds_per_calculation": requested_timeout_seconds,
        },
        "effective_timeout_seconds_per_calculation": max(0.0, effective_timeout),
        "retained_process_leases": leases,
        **({"time_accounting": "owned_native_execution_elapsed"} if transfer_control is not None else {}),
        "budget_contract": str(
            root / "surface-resolution-budgets"
            / str(contract["policy_scientific_hash"]).split(":", 1)[-1][:16]
            / "budget-contract.json"
        ),
    }
    snapshot_dir = Path(snapshot["budget_contract"]).parent / "snapshots"
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    existing = [
        int(path.stem) for path in snapshot_dir.glob("*.json") if path.stem.isdigit()
    ]
    _write_json(snapshot_dir / f"{max(existing, default=0) + 1:06d}.json", snapshot)
    return snapshot


def _surface_budget_gate(
    anchor: Path,
    *,
    requested_calculations: int,
    ncores: int,
    timeout_seconds: float | None,
) -> float | None:
    snapshot = _surface_budget_snapshot(
        anchor,
        requested_calculations=requested_calculations,
        requested_ncores=ncores,
        requested_timeout_seconds=timeout_seconds,
    )
    if snapshot is None:
        return timeout_seconds
    if snapshot["status"] == "surface_resolution_budget_exhausted":
        raise SurfaceResolutionBudgetExhausted(snapshot)
    return float(snapshot["effective_timeout_seconds_per_calculation"])


def _surface_context(anchor: Path) -> dict[str, Any] | None:
    """Read the immutable route context governing a nested strategy attempt."""
    path = Path(anchor).resolve()
    directory = path if path.is_dir() else path.parent
    for candidate in (directory, *directory.parents):
        context_paths = sorted(
            (candidate / "surface-resolution-contexts").glob("*.json")
        ) if (candidate / "surface-resolution-contexts").is_dir() else []
        legacy = candidate / ".surface-resolution-context.json"
        if legacy.is_file():
            context_paths.insert(0, legacy)
        if context_paths:
            validated: list[tuple[Path, dict[str, Any]]] = []
            for context_path in context_paths:
                value = json.loads(context_path.read_text())
                retained_hash = value.get("scientific_context_sha256")
                unhashed = {
                    key: item for key, item in value.items()
                    if key != "scientific_context_sha256"
                }
                calculated = hashlib.sha256(
                    json.dumps(unhashed, sort_keys=True, separators=(",", ":")).encode("utf-8")
                ).hexdigest()
                if retained_hash != calculated:
                    raise RuntimeError("Retained surface-resolution scientific context was mutated")
                validated.append((context_path, value))
            context_path, value = validated[-1]
            return {**value, "path": str(context_path)}
    return None


def _initialize_surface_context(
    workdir: Path,
    *,
    route: RouteSpec,
    mapping: dict[str, Any],
    environment_applicability: dict[str, Any],
    policy: SurfaceResolutionPolicy,
    ncores: int,
    method_keywords: list[str] | None,
    thermochemistry_temperature_K: float | None = None,
) -> dict[str, Any]:
    """Retain one immutable scientific context shared by every nested attempt."""
    def equivalent(left: Any, right: Any) -> bool:
        return json.dumps(
            left, sort_keys=True, separators=(",", ":"), default=str,
        ) == json.dumps(
            right, sort_keys=True, separators=(",", ":"), default=str,
        )

    def scientific_spec(record: dict[str, Any]) -> dict[str, Any]:
        return {
            key: record.get(key) for key in (
                "kind", "route", "mapping", "environment_applicability",
                "electronic_structure", "surface_resolution_policy",
                "claims_permitted", "claims_prohibited",
            )
        }

    context_path = Path(workdir) / ".surface-resolution-context.json"
    scientific_body: dict[str, Any] = {
        "schema_version": 1,
        "kind": "storca_surface_resolution_scientific_context",
        "route": route.as_dict(),
        "mapping": mapping,
        "environment_applicability": environment_applicability,
        "electronic_structure": {
            "backend": "orca",
            "method_keywords": list(
                method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]
            ),
            "ncores": int(ncores),
            "charge": int(route.charge),
            "multiplicity": int(route.multiplicity),
        },
        "surface_resolution_policy": policy.as_dict(),
        "claims_permitted": [],
        "claims_prohibited": ["rate", "lifetime"],
    }
    value = dict(scientific_body)
    if thermochemistry_temperature_K is not None:
        scientific_body["electronic_structure"]["thermochemistry_temperature_K"] = (
            normalize_thermochemistry_temperature(thermochemistry_temperature_K)
        )
    value["scientific_context_sha256"] = hashlib.sha256(
        json.dumps(scientific_body, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    retained = _surface_context(workdir)
    if retained is not None:
        if equivalent(scientific_spec(retained), scientific_spec(value)):
            return retained
        prior_electronic = dict(retained.get("electronic_structure") or {})
        next_electronic = dict(value.get("electronic_structure") or {})
        prior_ncores = prior_electronic.pop("ncores", None)
        next_ncores = next_electronic.pop("ncores", None)
        prior_policy = retained.get("surface_resolution_policy") or {}
        next_policy = value.get("surface_resolution_policy") or {}
        invariant_mismatch = (
            not equivalent(retained.get("route"), value.get("route"))
            or not equivalent(retained.get("mapping"), value.get("mapping"))
            or not equivalent(
                retained.get("environment_applicability"),
                value.get("environment_applicability"),
            )
            or not equivalent(prior_electronic, next_electronic)
        )
        if invariant_mismatch:
            raise RuntimeError(
                "Retained surface-resolution scientific context differs outside an auditable policy/resource amendment"
            )
        if prior_policy.get("policy_id") == next_policy.get("policy_id"):
            raise RuntimeError(
                "Retained surface-resolution policy definition changed without a new policy ID"
            )
        context_dir = Path(workdir) / "surface-resolution-contexts"
        context_dir.mkdir(parents=True, exist_ok=True)
        existing = sorted(context_dir.glob("*.json"))
        scientific_body.update({
            "schema_version": 2,
            "predecessor_context": retained.get("path"),
            "predecessor_scientific_context_sha256": retained.get(
                "scientific_context_sha256"
            ),
            "scientific_change": {
                "kind": "surface_resolution_policy_amendment",
                "prior_policy_id": prior_policy.get("policy_id"),
                "new_policy_id": next_policy.get("policy_id"),
                "execution_resource_change": (
                    {"prior_ncores": prior_ncores, "new_ncores": next_ncores}
                    if prior_ncores != next_ncores else None
                ),
            },
        })
        value = dict(scientific_body)
        value["scientific_context_sha256"] = hashlib.sha256(
            json.dumps(scientific_body, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        for existing_path in existing:
            existing_value = json.loads(existing_path.read_text())
            if equivalent(existing_value, value):
                return {**value, "path": str(existing_path)}
        retained_indexes = [
            int(path.stem) for path in existing if path.stem.isdigit()
        ]
        context_path = context_dir / f"{max(retained_indexes, default=1) + 1:06d}.json"
    elif context_path.is_file():
        # _surface_context would have returned this path; retain this guard for
        # malformed filesystem races.
        raise RuntimeError("Retained surface-resolution context could not be validated")
    _write_json(context_path, value)
    return {**value, "path": str(context_path)}


def _write_surface_attempt_evidence(
    attempt_root: Path,
    *,
    strategy: str,
    status: str,
    result: dict[str, Any],
    parameters: dict[str, Any] | None = None,
    input_paths: list[Path] | None = None,
    failure_reason: str | None = None,
    claims_permitted: tuple[str, ...] = (),
) -> dict[str, Any]:
    """Write an append-only, uniform evidence envelope for one strategy attempt."""
    attempt_root = Path(attempt_root)
    attempt_root.mkdir(parents=True, exist_ok=True)
    context = _surface_context(attempt_root)
    input_hashes: list[str] = []
    input_records: list[dict[str, str]] = []
    for candidate in input_paths or []:
        path = Path(candidate)
        if not path.is_file():
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        input_hashes.append("sha256:" + digest)
        input_records.append({"path": str(path), "sha256": digest})
    budget = _surface_budget_snapshot(attempt_root)
    leases: list[dict[str, Any]] = []
    if budget:
        root_resolved = attempt_root.resolve()
        for lease in budget.get("retained_process_leases") or []:
            try:
                Path(str(lease["lease"])).resolve().relative_to(root_resolved)
            except (KeyError, ValueError):
                continue
            leases.append(dict(lease))
    core_seconds = sum(
        float(item.get("elapsed_seconds") or 0.0) * int(item.get("ncores") or 1)
        for item in leases
    )
    result_hash = hashlib.sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    ).hexdigest()
    scientific_context = (
        {
            "path": str(context.get("path") or ""),
            "scientific_context_sha256": context["scientific_context_sha256"],
            "route_id": (context.get("route") or {}).get("route_id"),
            "reaction_equation": (context.get("route") or {}).get("reaction_equation"),
            "environment_applicability": context.get("environment_applicability"),
            "mapping_status": (context.get("mapping") or {}).get("status"),
            "mapping": context.get("mapping"),
            "electronic_structure": context.get("electronic_structure"),
            "surface_resolution_policy": context.get("surface_resolution_policy"),
        }
        if context else
        {"status": "direct_strategy_invocation_without_parent_route_context"}
    )
    attempt_id_seed = {
        "strategy": strategy,
        "attempt_root": str(attempt_root.resolve()),
        # The same retained calculation may be reinterpreted only under an
        # explicit append-only scientific-policy context.  Bind the envelope
        # identity to that context so a vN -> vN+1 policy amendment appends a
        # successor record instead of masquerading as mutation of vN.
        "scientific_context_sha256": scientific_context.get(
            "scientific_context_sha256"
        ),
        "input_evidence_hashes": input_hashes,
        "result_sha256": result_hash,
    }
    attempt_id = "surface-attempt-" + hashlib.sha256(
        json.dumps(attempt_id_seed, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:20]
    detailed_artifact_paths: list[str] = []
    for value in result.values():
        if not isinstance(value, (str, Path)):
            continue
        try:
            if Path(str(value)).exists():
                detailed_artifact_paths.append(str(value))
        except OSError:
            continue
    evidence = SurfaceResolutionAttempt(
        attempt_id=attempt_id,
        strategy=strategy,
        status=str(status),
        input_evidence_hashes=tuple(input_hashes),
        scientific_context=scientific_context,
        calculation_parameters={**(parameters or {}), "input_files": input_records},
        budget_consumed={
            "calculation_count": len(leases),
            "core_hours": core_seconds / 3600.0,
            "route_budget_status": budget.get("status") if budget else "not_available",
            "route_budget_contract": budget.get("budget_contract") if budget else None,
        },
        process_leases=tuple(leases),
        artifacts={
            "result_sha256_before_evidence_envelope": result_hash,
            "detailed_result_kind": result.get("kind"),
            "detailed_artifact_paths": sorted(set(detailed_artifact_paths)),
        },
        failure_reason=failure_reason,
        claims_permitted=claims_permitted,
        claims_prohibited=("rate", "lifetime"),
    ).as_dict()
    evidence_root = attempt_root / "surface-resolution-attempt-evidence"
    evidence_root.mkdir(parents=True, exist_ok=True)
    matching = sorted(evidence_root.glob("*.json"))
    for retained in matching:
        value = json.loads(retained.read_text())
        if value.get("attempt_id") == attempt_id:
            if value != evidence:
                raise RuntimeError("Retained surface-resolution attempt evidence was mutated")
            return {**evidence, "path": str(retained)}
    numbers = [int(path.stem) for path in matching if path.stem.isdigit()]
    evidence_path = evidence_root / f"{max(numbers, default=0) + 1:06d}.json"
    _write_json(evidence_path, evidence)
    return {**evidence, "path": str(evidence_path)}


def _normal_orca_output(path: Path) -> bool:
    return Path(path).is_file() and "ORCA TERMINATED NORMALLY" in Path(path).read_text(errors="replace")


_ORCA_RESTART_QUIET_SECONDS = 15 * 60


def _validate_orca_neb_restart_band(path: Path, *, expected_images: int | None = None) -> dict:
    """Validate ORCA's ``.allxyz`` band without interpreting its chemistry."""
    path = Path(path)
    if not path.is_file():
        return {"status": "missing", "valid": False, "path": str(path)}
    lines = [line for line in path.read_text(errors="replace").splitlines()
             if line.strip() and line.strip() != ">"]
    frames: list[list[str]] = []
    index = 0
    try:
        while index < len(lines):
            count = int(lines[index].strip())
            if count < 1 or index + count + 1 >= len(lines):
                raise ValueError("incomplete frame")
            atoms = []
            for row in lines[index + 2:index + 2 + count]:
                fields = row.split()
                if len(fields) < 4:
                    raise ValueError("incomplete atom row")
                if not all(math.isfinite(float(value)) for value in fields[1:4]):
                    raise ValueError("nonfinite checkpoint coordinates")
                atoms.append(fields[0])
            frames.append(atoms)
            index += count + 2
    except (TypeError, ValueError) as error:
        return {
            "status": "invalid_xyz_band", "valid": False, "path": str(path),
            "reason": str(error), "complete_frame_count": len(frames),
        }
    if not frames or any(elements != frames[0] for elements in frames[1:]):
        return {
            "status": "atom_order_changed_or_empty", "valid": False,
            "path": str(path), "complete_frame_count": len(frames),
        }
    if expected_images is not None and len(frames) != expected_images:
        return {
            "status": "image_count_mismatch", "valid": False, "path": str(path),
            "complete_frame_count": len(frames), "expected_frame_count": expected_images,
        }
    return {
        "status": "validated", "valid": True, "path": str(path),
        "complete_frame_count": len(frames), "atom_count": len(frames[0]),
        "expected_frame_count": expected_images,
        "elements": frames[0],
    }


NEB_TIME_LIMIT_CHECKPOINT_POLICY = "orca-neb-time-limit-checkpoint-continuation-v1"


def inspect_neb_time_limit_checkpoint(input_path: Path, *, timeout_seconds: float) -> dict | None:
    """Admit a hash-bound stopped NEB segment, never a path/rate claim.

    The caller must independently enforce its retained aggregate compute
    ceiling. This inspection does not grant a fresh budget or chemical retry.
    """
    inp = Path(input_path)
    output, contract_path, lease_path = (inp.with_suffix(suffix) for suffix in
                                         (".out", ".contract.json", ".process.json"))
    band = inp.parent / f"{inp.stem}_MEP.allxyz"
    if not all(p.is_file() for p in (inp, output, contract_path, lease_path, band)):
        return None
    try:
        timeout = float(timeout_seconds)
        if not math.isfinite(timeout) or timeout <= 0:
            return None
        text, out = inp.read_text(), output.read_text(errors="replace")
        if "NEB-TS" not in text.upper() or not any(
                _NEB_LBFGS_ITERATION.search(line) for line in out.splitlines()):
            return None
        if (re.search(r"ORCA TERMINATED NORMALLY|THE NEB OPTIMIZATION HAS CONVERGED", out, re.I)
                or _neb_intermediate_image_indices(output) or _neb_angle_warning_image_indices(output)
                or _neb_monotonic_endpoint_diagnostic(output)["valid"]
                or _neb_persistent_endpoint_highest_image_diagnostic(output, inp)["valid"]
                or _neb_submerged_climbing_image_diagnostic(output)["valid"]):
            return None
        lease = json.loads(lease_path.read_text())
        contract = json.loads(contract_path.read_text())
        if (lease.get("kind") != "storca_external_process_lease" or lease.get("engine") != "orca"
                or lease.get("status") != "terminated" or type(lease.get("returncode")) is not int
                or lease["returncode"] not in {-2, -9}
                or lease.get("input_sha256") != hashlib.sha256(inp.read_bytes()).hexdigest()
                or Path(lease.get("input", "")).resolve() != inp.resolve()
                or Path(lease.get("output", "")).resolve() != output.resolve()):
            return None
        started, ended = float(lease["started_at_epoch"]), float(lease["ended_at_epoch"])
        if not all(math.isfinite(x) for x in (started, ended)) or started <= 0 \
                or not timeout - .5 <= ended - started <= timeout + 20:
            return None
        for key, group in (("pid", False), ("process_group_id", True)):
            identifier = lease.get(key)
            if type(identifier) is not int or identifier <= 0:
                return None
            try:
                os.kill(-identifier if group else identifier, 0)
            except ProcessLookupError:
                pass
            else:
                return None
        if output.stat().st_mtime > ended + 1 or band.stat().st_mtime > ended + 1:
            return None
        names = set(re.findall(r'(?:Product|Hess_Filename|Restart_ALLXYZFile)\s+"([^"\n]+)"', text, re.I))
        xyz = re.search(r"\*\s+xyzfile\s+[-+]?\d+\s+\d+\s+(\S+)", text, re.I)
        images = re.search(r"(?im)^\s*NImages\s+(\d+)\s*$", text)
        if not xyz or not images:
            return None
        names.add(xyz.group(1).strip('"'))
        digest = hashlib.sha256(inp.read_bytes())
        dependencies = []
        def artifact(path):
            return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for name in sorted(names):
            path = inp.parent / name
            if not path.is_file():
                return None
            digest.update(name.encode())
            digest.update(path.read_bytes())
            dependencies.append(artifact(path))
        if contract.get("sha256") != digest.hexdigest():
            return None
        validation = _validate_orca_neb_restart_band(band, expected_images=int(images.group(1)) + 2)
        elements, coordinates, _ = read_xyz(inp.parent / xyz.group(1).strip('"'))
        if not validation["valid"] or validation["elements"] != elements or not np.isfinite(coordinates).all():
            return None
        return {"policy_id": NEB_TIME_LIMIT_CHECKPOINT_POLICY,
                "source_input": artifact(inp), "source_output": artifact(output),
                "source_job_contract": artifact(contract_path), "source_process_lease": artifact(lease_path),
                "source_band": artifact(band), "dependencies": dependencies, "band_validation": validation,
                "timeout_seconds": timeout, "elapsed_seconds": ended - started,
                "scientific_model_changed": False, "rate_claim_allowed": False,
                "renews_compute_budget": False}
    except (OSError, TypeError, ValueError, KeyError):
        return None


HISTORICAL_NEB_TS_SEED_POLICY = "historical-partial-neb-optts-trial-v1"


def _read_orca_neb_band_frames(path: Path) -> list[tuple[list[str], np.ndarray, str]]:
    """Read native allxyz separators without changing ordinary XYZ parsing."""
    lines = Path(path).read_text().splitlines()
    frames = []
    index = 0
    while index < len(lines):
        if not lines[index].strip() or lines[index].strip() == ">":
            index += 1
            continue
        count = int(lines[index].strip())
        if count <= 0 or index + count + 1 >= len(lines):
            raise ValueError("Incomplete native NEB band frame")
        rows = [line.split() for line in lines[index + 2:index + count + 2]]
        if any(len(row) != 4 for row in rows):
            raise ValueError("Invalid native NEB band coordinate row")
        coordinates = np.asarray([[float(value) for value in row[1:]] for row in rows])
        if not np.isfinite(coordinates).all():
            raise ValueError("Nonfinite native NEB band coordinates")
        frames.append(([row[0] for row in rows], coordinates, lines[index + 1]))
        index += count + 2
    integrity = validate_trajectory_frames(frames)
    if not integrity.get("valid"):
        raise ValueError("Native NEB band frame integrity failed: " + str(integrity.get("status")))
    return frames


def inspect_historical_neb_ts_seed(
    route: RouteSpec, request: dict, *, ncores: int,
    method_keywords: list[str] | None, condition_contract: dict | None,
) -> dict:
    """Read-only admission of an explicitly pinned *new* TS-search hypothesis.

    An incomplete historical band is not a converged NEB, barrier, or resume
    source for a newer execution.  This narrow strategy admits one bound
    species per side and never creates a budget or upgrades a prepared table.
    """
    result: dict[str, Any] = {
        "schema_version": 1, "kind": "historical_partial_neb_ts_seed_admission",
        "strategy_id": HISTORICAL_NEB_TS_SEED_POLICY,
        "status": "historical_ts_seed_refused", "valid": False,
        "rate_claim_allowed": False, "lifetime_claim_allowed": False,
        "neb_converged": False, "renews_compute_budget": False,
        "claims_prohibited": ["converged_neb", "transition_state", "barrier", "rate", "lifetime"],
    }
    def equivalent(left: Any, right: Any) -> bool:
        return json.dumps(left, sort_keys=True, default=str) == json.dumps(right, sort_keys=True, default=str)

    try:
        if request.get("schema_version") != 1 or request.get("kind") != "historical_partial_neb_ts_seed_request":
            raise ValueError("Explicit historical partial-band request required")
        references: dict[str, dict] = {}
        checkpoint = request["source_checkpoint"]
        for reference in [request["source_path_result"], request["source_execution_contract"],
                          *request["source_artifacts"],
                          *(checkpoint[key] for key in ("source_input", "source_output", "source_job_contract",
                                                       "source_process_lease", "source_band")),
                          *checkpoint["dependencies"]]:
            path = Path(reference["path"]).resolve(strict=True)
            previous = references.get(str(path))
            if previous and previous["sha256"] != reference["sha256"]:
                raise ValueError("Conflicting source hashes in historical request")
            if hashlib.sha256(path.read_bytes()).hexdigest() != reference["sha256"]:
                raise ValueError(f"Historical source changed: {path}")
            references[str(path)] = {"path": str(path), "sha256": reference["sha256"]}

        def pinned(path: Path | str) -> Path:
            path = Path(path).resolve(strict=True)
            if str(path) not in references:
                raise ValueError(f"Historical source lacks an immutable hash: {path}")
            return path

        source = json.loads(pinned(request["source_path_result"]["path"]).read_text())
        execution_path = pinned(request["source_execution_contract"]["path"])
        execution = json.loads(execution_path.read_text())
        method = list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"])
        if (not equivalent(source.get("route"), route.as_dict())
                or not equivalent(execution.get("condition_contract"), condition_contract)
                or (execution.get("method_profile") or {}).get("orca_keywords") != method
                or (execution.get("resources") or {}).get("ncores") != ncores
                or execution.get("parent") != {"smiles": route.parent_smiles,
                                               "charge": route.charge, "multiplicity": route.multiplicity}):
            raise ValueError("Historical route, state, method, conditions or resources differ")
        if len(route.reactant_smiles) != 1 or len(route.product_smiles) != 1:
            raise ValueError("Historical band strategy currently requires one bound species per side")
        inp = pinned(checkpoint["source_input"]["path"])
        if not inp.is_relative_to(execution_path.parent):
            raise ValueError("Historical NEB source is outside its execution workspace")
        timeout = (execution.get("resources") or {}).get("route_job_timeout_seconds")
        if checkpoint.get("timeout_seconds") != timeout:
            raise ValueError("Historical checkpoint timeout differs from the execution contract")
        current = inspect_neb_time_limit_checkpoint(inp, timeout_seconds=timeout)
        if current is None or not equivalent(current, checkpoint):
            raise ValueError("Historical NEB lacks its exact terminal time-limit checkpoint")
        if _orca_input_ncores(inp) != ncores:
            raise ValueError("Historical NEB core count differs from the execution contract")
        lease = json.loads(pinned(checkpoint["source_process_lease"]["path"]).read_text())
        neb_text = inp.read_text()
        neb_tokens = " ".join(re.findall(r"^\s*!\s*(.*)$", neb_text, re.M)).lower().split()
        neb_states = re.findall(r"^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)", neb_text, re.M | re.I)
        if (lease.get("ncores") != ncores or "neb-ts" not in neb_tokens
                or sorted(token for token in neb_tokens if token not in {"neb-ts", "tightscf"})
                != sorted(token.lower() for token in method if token.lower() != "tightscf")
                or neb_states != [(str(route.charge), str(route.multiplicity))]):
            raise ValueError("Historical NEB method, state or lease resources differ")
        orientations = [item for item in source.get("orientations") or []
                        if Path((item.get("neb") or {}).get("input") or "").resolve() == inp]
        if len(orientations) != 1:
            raise ValueError("Historical path does not uniquely own the requested NEB execution")
        orientation = orientations[0]
        if ((orientation.get("neb") or {}).get("status") != "failed"
                or Path(orientation["neb"]["output"]).resolve() != inp.with_suffix(".out")
                or not (orientation.get("endpoint_validation") or {}).get("valid")):
            raise ValueError("Historical endpoint/NEB execution ownership is unresolved")
        mapping = resolve_route_atom_mapping(route)
        context_path = pinned(source["surface_resolution_scientific_context"]["path"])
        context = json.loads(context_path.read_text())
        electronic = context.get("electronic_structure") or {}
        if (not mapping.get("valid") or not equivalent(context.get("route"), route.as_dict())
                or not equivalent(context.get("mapping"), mapping)
                or electronic.get("method_keywords") != method or electronic.get("ncores") != ncores
                or electronic.get("charge") != route.charge or electronic.get("multiplicity") != route.multiplicity):
            raise ValueError("Historical mapping or electronic scientific context differs")
        for reference in (mapping.get("source_evidence") or {}).values():
            if isinstance(reference, dict) and "path" in reference:
                path = pinned(reference["path"])
                if references[str(path)]["sha256"] != reference.get("sha256"):
                    raise ValueError("Historical atom-map evidence hash differs")
        applicability = _surface_environment_applicability(route, condition_contract, SurfaceResolutionPolicy(), method)
        if (applicability.get("status") != "applicable"
                or not equivalent(applicability, context.get("environment_applicability"))):
            raise ValueError("Historical environment model differs or is inapplicable")

        stationary = {"status": "validated", "valid": True, "reactants": [], "products": []}
        endpoints = {}
        validations = {"mode": "revalidated_historical_bound_endpoint_minima", "valid": True}
        for side, singular, smiles_values, labels in (
            ("reactants", "reactant", route.reactant_smiles, route.reactant_labels),
            ("products", "product", route.product_smiles, route.product_labels),
        ):
            old = orientation["endpoint_validation"][singular]
            xyz = pinned(old["endpoint_xyz"])
            folder = xyz.parent
            paths = {name: pinned(folder / name) for name in (
                "metadata.json", "endpoint-contract.json", "input.xyz", "opt.inp", "opt.out", "opt.xyz",
                "freq.inp", "freq.out", "opt.contract.json", "freq.contract.json",
            )}
            if (pinned(old["optimization_output"]) != paths["opt.out"]
                    or pinned(old["frequency_output"]) != paths["freq.out"]):
                raise ValueError("Historical endpoint output ownership differs")
            endpoint_contract = json.loads(paths["endpoint-contract.json"].read_text())
            if (endpoint_contract.get("charge") != route.charge
                    or endpoint_contract.get("multiplicity") != route.multiplicity
                    or endpoint_contract.get("ncores") != ncores
                    or endpoint_contract.get("method_keywords") != method
                    or endpoint_contract.get("input_xyz_sha256") != references[str(paths["input.xyz"])]["sha256"]):
                raise ValueError("Historical endpoint contract differs")
            retained = _retained_endpoint_workflow(folder, charge=route.charge, multiplicity=route.multiplicity,
                                                   ncores=ncores, method_keywords=method)
            if (not retained or not (retained.get("frequency_check") or {}).get("IsMinimum")
                    or not _orca_optimization_converged(paths["opt.out"])):
                raise ValueError("Historical endpoint is not a converged, frequency-validated minimum")
            for name, job, geometry in (("opt.inp", "opt", paths["input.xyz"]),
                                        ("freq.inp", "freq", xyz)):
                text = paths[name].read_text()
                actual = " ".join(re.findall(r"^\s*!\s*(.*)$", text, re.M)).lower().split()
                states = re.findall(r"^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)", text, re.M | re.I)
                if (job not in actual or sorted(token for token in actual if token not in {job, "tightscf"})
                        != sorted(token.lower() for token in method if token.lower() != "tightscf")
                        or states != [(str(route.charge), str(route.multiplicity))]
                        or _orca_input_ncores(paths[name]) != ncores):
                    raise ValueError("Historical endpoint job method, state or resources differ")
                _validate_retained_orca_job_input(paths[name])
                inline = re.search(r"^\s*\*\s+xyz\s+-?\d+\s+\d+\s*\n(.*?)^\s*\*",
                                   text, re.M | re.S | re.I)
                rows = [line.split() for line in inline[1].strip().splitlines()] if inline else []
                elements, coordinates, _ = read_xyz(geometry)
                if (not rows or any(len(row) != 4 for row in rows) or [row[0] for row in rows] != elements
                        or not np.allclose([[float(v) for v in row[1:]] for row in rows],
                                           coordinates, atol=1e-5, rtol=0)):
                    raise ValueError("Historical endpoint input geometry differs")
            connectivity = validate_route_endpoint_connectivity(route, xyz, side=side, mapping_result=mapping)
            if not connectivity.get("valid"):
                raise ValueError("Historical endpoint connectivity or atom order differs")
            label = str(labels[0] if labels else smiles_values[0])
            record = {"label": label, "original_label": label, "stoichiometric_copy": 1,
                      "stoichiometric_copy_count": 1, "smiles": smiles_values[0], "charge": route.charge,
                      "multiplicity": route.multiplicity, "input_xyz": str(paths["input.xyz"]),
                      "status": "validated_minimum", "optimized_xyz": str(xyz),
                      "optimization_output": str(paths["opt.out"]), "frequency_output": str(paths["freq.out"]),
                      "frequency_check": retained["frequency_check"], "connectivity_check": connectivity,
                      "revalidated_historical_source": True}
            stationary[side].append(record)
            endpoints[singular] = xyz
            validations[singular] = {**record, "endpoint_xyz": str(xyz)}
        band = pinned(checkpoint["source_band"]["path"])
        frames = _read_orca_neb_band_frames(band)
        reactant_elements, reactant_coordinates, _ = read_xyz(endpoints["reactant"])
        product_elements, product_coordinates, _ = read_xyz(endpoints["product"])
        integrity = validate_trajectory_frames(frames, expected_elements=reactant_elements)
        tolerance = SurfaceResolutionPolicy().endpoint_rmsd_tolerance_angstrom
        distances = (aligned_rmsd(reactant_coordinates, frames[0][1]),
                     aligned_rmsd(product_coordinates, frames[-1][1]))
        band_endpoints = {"valid": integrity.get("valid") and product_elements == reactant_elements
                         and all(value <= tolerance for value in distances),
                         "orientation": "forward", "reactant_rmsd_angstrom": distances[0],
                         "product_rmsd_angstrom": distances[1], "per_endpoint_rmsd_tolerance_angstrom": tolerance,
                         "trajectory_integrity": integrity, "claim_scope": "coordinate_correspondence_only"}
        if not band_endpoints["valid"]:
            raise ValueError("Historical band endpoints differ from validated, ordered minima")
        energies = []
        for _, _, comment in frames:
            matches = re.findall(r"(?:energy|\bE)\s*(?:=|:)?\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
                                 comment, re.I)
            if not matches:
                raise ValueError("Native NEB band energy profile is incomplete")
            energies.append(float(matches[-1]))
        selection = select_interior_scan_maximum(energies)
        if selection.get("status") != "interior_ts_seed_candidate" or len(energies or []) != len(frames):
            raise ValueError("Historical band lacks a finite, unique bracketed interior maximum")
        plan = build_bond_change_scan_plan(reactant_coordinates, product_coordinates,
                                           source["endpoint_complexes"]["bond_changes"])
        if not plan.get("coordinates"):
            raise ValueError("Historical band lacks mapped intended reaction coordinates")
        from .route_geometry import _build_endpoint_side, _mapped_bond_changes
        expected_changes = _mapped_bond_changes(
            _build_endpoint_side(route.reactant_smiles), _build_endpoint_side(route.product_smiles),
            {int(left): int(right) for left, right in mapping["all_atom_mapping"]},
        )
        if not equivalent(source["endpoint_complexes"]["bond_changes"], expected_changes):
            raise ValueError("Historical intended coordinates differ from the exact atom-mapped reaction")
        elements, coordinates, comment = frames[selection["index"]]
        result.update(valid=True, status="historical_ts_seed_admitted", mapping=mapping,
                      source_checkpoint=current, source_artifacts=list(references.values()),
                      source_execution_contract=request["source_execution_contract"], source_band=current["source_band"],
                      seed_index=selection["index"], seed_elements=elements, seed_coordinates=coordinates.tolist(),
                      seed_comment=comment, seed_selection=selection,
                      reactant_xyz=str(endpoints["reactant"]), product_xyz=str(endpoints["product"]),
                      stationary_points=stationary, endpoint_validation=validations,
                      intended_coordinates=plan["coordinates"], band_endpoint_correspondence=band_endpoints)
        if request.get("retained_ts_frequency") is not None:
            reuse = inspect_historical_ts_frequency_reuse(
                route, request["retained_ts_frequency"], historical_admission=result,
                ncores=ncores, method_keywords=method_keywords,
            )
            result["retained_ts_frequency"] = reuse
            if not reuse["valid"]:
                result.update(valid=False, status="historical_ts_seed_refused",
                              failure_reason=reuse.get("failure_reason"))
    except (KeyError, IndexError, TypeError, ValueError, OSError, RuntimeError) as error:
        result["failure_reason"] = f"{type(error).__name__}: {error}"
    return result


def inspect_historical_ts_frequency_reuse(
    route: RouteSpec, request: dict, *, historical_admission: dict,
    ncores: int, method_keywords: list[str] | None,
    surface_resolution_policy: SurfaceResolutionPolicy | None = None,
) -> dict:
    """Reopen completed TS stages without borrowing an incomplete IRC claim.

    This is a read-only stage-reuse gate. A new IRC still needs its own bounded
    execution and ordinary bidirectional/endpoint checks. No source is copied,
    rewritten or reexecuted here; resource authorization remains the caller's.
    """
    result = {"schema_version": 1, "kind": "completed_ts_frequency_reuse_admission",
              "valid": False, "status": "ts_frequency_reuse_refused",
              "claims_prohibited": ["verified_path", "rate", "lifetime"],
              "renews_compute_budget": False}
    policy = surface_resolution_policy or SurfaceResolutionPolicy()
    equivalent = lambda left, right: json.dumps(left, sort_keys=True) == json.dumps(right, sort_keys=True)
    try:
        if (not historical_admission.get("valid") or request.get("schema_version") != 1
                or request.get("kind") != "completed_ts_frequency_reuse_request"):
            raise ValueError("Explicit source-bound completed-stage reuse request required")
        references = {}
        for reference in [request["source_path_result"], request["source_execution_contract"],
                          request["source_run_result"], *request["source_artifacts"]]:
            path = Path(reference["path"]).resolve(strict=True)
            key = str(path)
            if key in references and references[key]["sha256"] != reference["sha256"]:
                raise ValueError("Conflicting retained TS source hashes")
            if hashlib.sha256(path.read_bytes()).hexdigest() != reference["sha256"]:
                raise ValueError(f"Retained TS source changed: {path}")
            references[key] = {"path": key, "sha256": reference["sha256"]}

        def pinned(value):
            path = Path(value).resolve(strict=True)
            if str(path) not in references:
                raise ValueError(f"Retained TS source lacks a hash: {path}")
            return path

        source_path = pinned(request["source_path_result"]["path"])
        contract_path = pinned(request["source_execution_contract"]["path"])
        run_path = pinned(request["source_run_result"]["path"])
        source, contract, run = (json.loads(path.read_text()) for path in (source_path, contract_path, run_path))
        original = json.loads(Path(historical_admission["source_execution_contract"]["path"]).read_text())
        method = list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"])
        if (source.get("kind") != "generic_orca_route_verification"
                or not equivalent(source.get("route"), route.as_dict())
                or not source_path.is_relative_to(contract_path.parent)
                or run_path.parent != contract_path.parent
                or contract.get("condition_contract") != original.get("condition_contract")
                or contract.get("method_profile") != original.get("method_profile")
                or (contract.get("budget") or {}).get("maximum_cores") != ncores
                or run.get("status") not in {"workflow_returned", "failed"}
                or run.get("original_sources_unchanged") is not True
                or run.get("manifest_unchanged") is not True):
            raise ValueError("Retained TS execution does not own the exact terminal scientific contract")
        prior = source["historical_neb_seed_admission"]
        if (not prior.get("valid") or prior.get("seed_index") != historical_admission["seed_index"]
                or not equivalent(prior.get("source_checkpoint"), historical_admission["source_checkpoint"])):
            raise ValueError("Retained TS was not generated from this admitted historical image")
        context = _surface_context(source_path.parent)
        if context is None or not equivalent(context.get("route"), route.as_dict()) \
                or not equivalent(context.get("mapping"), historical_admission["mapping"]):
            raise ValueError("Retained TS scientific mapping/context differs")
        declared_context = source["surface_resolution_scientific_context"]
        if (pinned(declared_context["path"]) != pinned(context["path"])
                or declared_context.get("scientific_context_sha256") != context.get("scientific_context_sha256")):
            raise ValueError("Retained TS declared context pointer or scientific hash differs")
        electronic = context["electronic_structure"]
        if (electronic.get("method_keywords") != method or electronic.get("ncores") != ncores
                or electronic.get("charge") != route.charge or electronic.get("multiplicity") != route.multiplicity):
            raise ValueError("Retained TS method, state or execution resources differ")
        ts = source["historical_neb_seed_trial"]["transition_state"]
        optimization, frequency = ts["optimization"], ts["frequency_execution"]
        paths = {"seed": pinned(ts["seed_xyz"]), "xyz": pinned(ts["optimized_xyz"]),
                 "opt": pinned(optimization["input"]), "freq": pinned(frequency["input"]),
                 "opt_output": pinned(optimization["output"]), "freq_output": pinned(ts["frequency_output"])}
        folder = paths["opt"].parent
        if any(path.parent != folder for path in paths.values()) or not folder.is_relative_to(source_path.parent):
            raise ValueError("Retained TS artifacts do not share their owned execution workspace")
        paths["hessian"] = pinned(paths["freq"].with_suffix(".hess"))
        for inp, out, job in ((paths["opt"], paths["opt_output"], "optts"),
                              (paths["freq"], paths["freq_output"], "freq"),
                              (pinned(folder / "irc.inp"), pinned(folder / "irc.out"), "irc")):
            if out != inp.with_suffix(".out"):
                raise ValueError("Retained TS job output ownership differs")
            text = inp.read_text()
            tokens = " ".join(re.findall(r"^\s*!\s*(.*)$", text, re.M)).lower().split()
            if (job not in tokens or sorted(token for token in tokens if token not in {job, "tightscf"})
                    != sorted(token.lower() for token in method if token.lower() != "tightscf")
                    or re.findall(r"^\s*\*\s+xyz(?:file)?\s+(-?\d+)\s+(\d+)", text, re.M | re.I)
                    != [(str(route.charge), str(route.multiplicity))]
                    or _orca_input_ncores(inp) != ncores):
                raise ValueError("Retained TS job method, state or resource declaration differs")
            job_contract = json.loads(pinned(inp.with_suffix(".contract.json")).read_text())
            _, dependencies = _orca_job_input_evidence(inp)
            if (job_contract.get("schema_version") != 1 or job_contract.get("kind") != "orca_job_input_contract"
                    or Path(job_contract.get("input", "")).resolve() != inp
                    or job_contract.get("referenced_files") != dependencies):
                raise ValueError("Retained TS job contract ownership or dependencies differ")
            for dependency in dependencies:
                pinned(dependency)
            _validate_retained_orca_job_input(inp)
            lease = json.loads(pinned(inp.with_suffix(".process.json")).read_text())
            if (lease.get("kind") != "storca_external_process_lease" or lease.get("engine") != "orca"
                    or lease.get("status") != "terminated" or type(lease.get("returncode")) is not int
                    or lease.get("input_sha256") != references[str(inp)]["sha256"]
                    or Path(lease.get("input", "")).resolve() != inp
                    or Path(lease.get("output", "")).resolve() != out or lease.get("ncores") != ncores
                    or not all(math.isfinite(float(lease[key])) for key in ("started_at_epoch", "ended_at_epoch"))
                    or float(lease["ended_at_epoch"]) < float(lease["started_at_epoch"])):
                raise ValueError("Retained TS job lease is not terminal with exact ownership")
            for key, group in (("pid", False), ("process_group_id", True)):
                identifier = lease.get(key)
                if type(identifier) is not int or identifier <= 0:
                    raise ValueError("Retained TS lease lacks positive process identifiers")
                try:
                    os.kill(-identifier if group else identifier, 0)
                except ProcessLookupError:
                    pass
                else:
                    raise ValueError("Retained TS source process may still be alive")
            if job != "irc" and (lease["returncode"] != 0 or not _normal_orca_output(out)):
                raise ValueError("Retained TS stage did not terminate normally")
        if not _orca_optimization_converged(paths["opt_output"]):
            raise ValueError("Retained TS optimization is unconverged")
        _, ts_dependencies = _orca_job_input_evidence(paths["opt"])
        if (ts_dependencies != [str(paths["seed"])] or paths["xyz"] != paths["opt"].with_suffix(".xyz")
                or Path(frequency.get("output", "")).resolve() != paths["freq_output"]):
            raise ValueError("Retained TS optimization seed/geometry or frequency output ownership differs")
        labels, coordinates, _ = read_xyz(paths["xyz"])
        seed_labels, seed_coordinates, _ = read_xyz(paths["seed"])
        if (seed_labels != historical_admission["seed_elements"]
                or not np.allclose(seed_coordinates, historical_admission["seed_coordinates"], atol=1e-7, rtol=0)
                or labels != seed_labels or not np.isfinite(coordinates).all()):
            raise ValueError("Retained TS seed does not match the admitted historical image")
        modes = parse_orca_normal_modes(paths["hessian"])
        imaginary = np.flatnonzero(modes.frequencies_cm_1 < policy.meaningful_imaginary_frequency_cm_1)
        check = frequency_stability_check(paths["freq_output"], imag_tol=abs(policy.meaningful_imaginary_frequency_cm_1))
        printed_frequencies = np.asarray(check.get("AllFrequencies") or [], dtype=float)
        hessian_xyz = modes.coordinates_bohr * 0.529177210903
        if (check.get("NumImaginary") != 1 or len(imaginary) != 1 or tuple(labels) != modes.atom_labels
                or printed_frequencies.shape != modes.frequencies_cm_1.shape
                or not np.allclose(printed_frequencies, modes.frequencies_cm_1, atol=0.05, rtol=0)
                or not np.isfinite(modes.frequencies_cm_1).all() or not np.isfinite(modes.vectors).all()
                or not np.isfinite(hessian_xyz).all()
                or not np.allclose(hessian_xyz + coordinates.mean(axis=0) - hessian_xyz.mean(axis=0),
                                   coordinates, atol=1e-5, rtol=0)):
            raise ValueError("Retained TS Hessian, geometry or imaginary mode differs")
        inline = re.search(r"^\s*\*\s+xyz\s+-?\d+\s+\d+\s*\n(.*?)^\s*\*",
                           paths["freq"].read_text(), re.M | re.S | re.I)
        rows = [line.split() for line in inline[1].strip().splitlines()] if inline else []
        if (len(rows) != len(labels) or any(len(row) != 4 for row in rows) or [row[0] for row in rows] != labels
                or not np.allclose([[float(v) for v in row[1:]] for row in rows], coordinates, atol=1e-5, rtol=0)):
            raise ValueError("Retained frequency input geometry differs from the optimized TS")
        overlap = bond_change_mode_overlap(coordinates, modes.vectors[int(imaginary[0])],
                                           historical_admission["intended_coordinates"])
        if overlap["status"] != "matched":
            raise ValueError("Retained TS imaginary mode does not match the intended mapped reaction")
        reopened = {**ts, "optimization": {**optimization, "resumed": True},
                    "frequency_execution": {**frequency, "resumed": True},
                    "frequency_check": check, "reaction_coordinate_overlap": overlap}
        result.update(valid=True, status="ts_frequency_stages_reusable", transition_state=reopened,
                      source_artifacts=list(references.values()), source_path_result=request["source_path_result"],
                      reused_stages=["ts_optimization", "ts_frequency"],
                      irc_reused=False, irc_completion=_orca_irc_completion(folder / "irc.out"))
    except (KeyError, IndexError, TypeError, ValueError, OSError, RuntimeError) as error:
        result["failure_reason"] = f"{type(error).__name__}: {error}"
    return result


def _prepare_neb_restart_input(input_path: Path, band_path: Path) -> tuple[Path, dict]:
    """Create an append-only ORCA NEB restart input from a validated band."""
    input_path, band_path = Path(input_path), Path(band_path)
    text = input_path.read_text(errors="replace")
    match = re.search(r"(?im)^\s*NImages\s+(\d+)\s*$", text)
    expected = int(match.group(1)) + 2 if match else None
    validation = _validate_orca_neb_restart_band(band_path, expected_images=expected)
    if not validation["valid"]:
        raise ValueError(f"ORCA NEB restart band is invalid: {validation['status']}")
    # ORCA 6.1 forbids endpoint preoptimization when reading a saved band.
    # The band already retains the optimized, fixed endpoint structures.
    if len(re.findall(r"(?im)^\s*PreOpt\s+\S+\s*$", text)) > 1:
        raise ValueError("ORCA NEB restart has ambiguous endpoint preoptimization")
    text = re.sub(r"(?im)^(\s*PreOpt\s+)true(\s*)$", r"\1false\2", text)
    root_stem = re.sub(r"-restart-\d{3}$", "", input_path.stem)
    existing = []
    for path in input_path.parent.glob(f"{root_stem}-restart-*.inp"):
        suffix = re.search(r"-restart-(\d{3})$", path.stem)
        if suffix:
            existing.append(int(suffix.group(1)))
    restart_input = input_path.parent / f"{root_stem}-restart-{max(existing, default=0) + 1:03d}.inp"
    relative_band = band_path.relative_to(input_path.parent)
    if re.search(r"(?im)^\s*Restart_ALLXYZFile\b", text):
        text = re.sub(
            r'(?im)^(\s*Restart_ALLXYZFile\s+)"[^"]+"\s*$',
            rf'\1"{relative_band}"', text, count=1,
        )
    else:
        text, count = re.subn(
            r"(?im)^(\s*%neb\s*)$",
            rf'\1\n  Restart_ALLXYZFile "{relative_band}"', text, count=1,
        )
        if count != 1:
            raise ValueError("ORCA NEB input has no unique %neb block")
    with restart_input.open("x") as stream:
        stream.write(text)
    evidence = {
        "schema_version": 1,
        "kind": "orca_neb_checkpoint_restart",
        "source_input": str(input_path),
        "source_output": str(input_path.with_suffix(".out")),
        "source_band": str(band_path),
        "source_band_sha256": hashlib.sha256(band_path.read_bytes()).hexdigest(),
        "source_input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "source_output_sha256": hashlib.sha256(input_path.with_suffix(".out").read_bytes()).hexdigest(),
        "restart_input": str(restart_input),
        "band_validation": validation,
        "endpoint_preoptimization_disabled_for_restart": True,
        "scientific_model_changed": False,
    }
    for suffix, key in ((".contract.json", "source_job_contract"), (".process.json", "source_process_lease")):
        source = input_path.with_suffix(suffix)
        if source.is_file():
            evidence[key] = {"path": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
    _write_json(restart_input.with_suffix(".restart.json"), evidence)
    return restart_input, evidence


def _orca_optimization_converged(path: Path) -> bool:
    if not Path(path).is_file():
        return False
    text = Path(path).read_text(errors="replace").upper()
    return "THE OPTIMIZATION HAS CONVERGED" in text


def _orca_neb_no_barrier_outcome(path: Path) -> dict:
    """Recognize ORCA's semantically successful no-saddle NEB outcome.

    ORCA 6 can return a nonzero process status after a converged elastic band
    when the highest-energy image is an endpoint.  In that case it explicitly
    says that no barrier was found and skips the NEB-TS refinement.  The process
    status is not normal termination, but the retained MEP is still usable
    evidence once its images, energies, and endpoints are validated below.
    """
    output = Path(path)
    if not output.is_file():
        return {"status": "not_detected", "valid": False}
    text = output.read_text(errors="replace")
    converged = "THE NEB OPTIMIZATION HAS CONVERGED" in text
    no_barrier = bool(re.search(
        r"No\s+barrier\s+was\s+found\.\s*Skipping\s+NEB-TS\s+run\s+here\.",
        text,
        re.IGNORECASE,
    ))
    if not (converged and no_barrier):
        return {
            "status": "not_detected",
            "valid": False,
            "neb_converged": converged,
            "explicit_no_barrier_message": no_barrier,
        }
    return {
        "status": "converged_no_interior_barrier",
        "valid": True,
        "neb_converged": True,
        "explicit_no_barrier_message": True,
        "normal_termination": _normal_orca_output(output),
        "climbing_image_activated": "climbing image was never activated" not in text.lower(),
        "interpretation": (
            "ORCA converged the minimum-energy path and found no interior maximum, so it skipped "
            "NEB-TS. Reaction direction and any forward energetic threshold must be determined from "
            "the validated endpoint energy profile; this does not by itself establish a forward rate."
        ),
    }


_NEB_INTERMEDIATE_WARNING = re.compile(
    r"Possible\s+intermediate\s+minimum\s+found\s+at\s+image\(s\):",
    re.IGNORECASE,
)
_NEB_ANGLE_WARNING = re.compile(
    r"Warning:\s*maximum\s+angle\s+between\s+images\s*\(([^)]*)\)",
    re.IGNORECASE,
)
_NEB_MONOTONIC_WARNING = re.compile(
    r"Warning:\s*path\s+is\s+strictly\s+increasing\s+in\s+energy\s+"
    r"LBFGS\s+(\d+)\s+(\d+)",
    re.IGNORECASE,
)
_NEB_LBFGS_ITERATION = re.compile(r"^\s*(?:Warning:.*?\s+)?LBFGS\s+\d+\s+\d+", re.IGNORECASE)
_NEB_MONOTONIC_STOP_STREAK = 12
_NEB_ENDPOINT_MAX_STOP_STREAK = 24
_NEB_ENDPOINT_MAX_LOW_FORCE_COUNT = 8
_NEB_SUBMERGED_CI_STOP_STREAK = 12
_ORCA_STARTUP_FAILURE_MARKERS = (
    "ORCA finished by error termination in Startup",
    "bind() failed for port 0: Operation not permitted",
    "No sockets were able to be opened on the available protocols",
    "ptl_tool: problems getting address",
    "no available interfaces were found",
    "You requested support for tool connections",
)
_ORCA_EXECUTION_RETRY_LIMIT = 2
_ORCA_SAFE_MPI_ENVIRONMENT = {
    "OMPI_MCA_pml": "ob1",
    "OMPI_MCA_btl": "self,tcp",
}


def _orca_startup_infrastructure_failure(output: Path) -> bool:
    output = Path(output)
    if not output.is_file():
        return False
    text = output.read_text(errors="replace")
    return any(marker in text for marker in _ORCA_STARTUP_FAILURE_MARKERS)


def _orca_serial_execution_input(input_text: str) -> str:
    """Return the same ORCA model with its execution parallelism set to one."""
    updated, count = re.subn(
        r"(?im)^(\s*nprocs\s+)\d+(\s*)$", r"\g<1>1\g<2>", input_text,
    )
    if count:
        return updated
    # Inputs without a PAL block already execute serially.  Do not invent a
    # scientific or syntactic change merely to label the fallback.
    return input_text


def _orca_parallel_execution_input(input_text: str, ncores: int) -> str:
    """Restore a retained PAL contract without changing scientific keywords."""
    if int(ncores) < 1:
        raise ValueError("ORCA execution recovery needs at least one core")
    updated, count = re.subn(
        r"(?im)^(\s*nprocs\s+)\d+(\s*)$",
        rf"\g<1>{int(ncores)}\g<2>",
        input_text,
    )
    if not count and int(ncores) != 1:
        raise RuntimeError("Cannot restore parallelism: retained ORCA input has no PAL nprocs line")
    return updated


def _prepare_orca_execution_retry(
    input_path: Path,
    *,
    output: Path,
    existing_retry_inputs: list[Path],
) -> tuple[Path, dict]:
    """Create one append-only infrastructure retry with an auditable policy."""
    number = len(existing_retry_inputs) + 1
    if number > _ORCA_EXECUTION_RETRY_LIMIT:
        raise RuntimeError(
            "ORCA startup infrastructure retry budget exhausted after "
            f"{_ORCA_EXECUTION_RETRY_LIMIT} attempts: {output}"
        )
    root_stem = re.sub(r"-execution-retry-\d{3}$", "", input_path.stem)
    retry_input = input_path.parent / f"{root_stem}-execution-retry-{number:03d}.inp"
    if retry_input.exists():
        raise RuntimeError(f"Execution retry input collision: {retry_input}")
    source_text = input_path.read_text(errors="replace")
    serial_fallback = number == _ORCA_EXECUTION_RETRY_LIMIT
    retry_text = _orca_serial_execution_input(source_text) if serial_fallback else source_text
    retry_input.write_text(retry_text)
    execution_environment = dict(_ORCA_SAFE_MPI_ENVIRONMENT)
    policy = {
        "schema_version": 1,
        "kind": "storca_orca_execution_retry_policy",
        "retry_number": number,
        "strategy": (
            "serial_execution_fallback" if serial_fallback
            else "safe_openmpi_transport_retry"
        ),
        "extra_environment": execution_environment,
        "requested_ncores": _orca_input_ncores(input_path),
        "effective_ncores": _orca_input_ncores(retry_input),
        "scientific_model_changed": False,
        "execution_resources_changed": retry_text != source_text,
    }
    policy_path = retry_input.with_suffix(".execution-policy.json")
    _write_json(policy_path, policy)
    retry_evidence = {
        "schema_version": 2,
        "kind": "orca_startup_infrastructure_execution_retry",
        "source_input": str(input_path),
        "source_output": str(output),
        "source_output_sha256": (
            hashlib.sha256(output.read_bytes()).hexdigest() if output.is_file() else None
        ),
        "retry_input": str(retry_input),
        "execution_policy": str(policy_path),
        "strategy": policy["strategy"],
        "scientific_input_changed": False,
        "execution_resources_changed": policy["execution_resources_changed"],
        "reason": "retained_orca_startup_infrastructure_failure",
    }
    _write_json(retry_input.with_suffix(".restart.json"), retry_evidence)
    return retry_input, retry_evidence


def _detached_orca_execution_is_stopped(input_path: Path) -> bool:
    lease_path = Path(input_path).with_suffix(".process.json")
    if not lease_path.is_file():
        return False
    try:
        lease = json.loads(lease_path.read_text())
        if lease.get("status") not in {
            "running", "runner_detached_process_may_be_active",
        }:
            return False
        pid = int(lease.get("pid"))
        process_group_id = int(lease.get("process_group_id") or pid)
        try:
            os.kill(pid, 0)
            return False
        except ProcessLookupError:
            try:
                os.killpg(process_group_id, 0)
                return False
            except ProcessLookupError:
                return True
    except (OSError, PermissionError, TypeError, ValueError, json.JSONDecodeError):
        return False


def _prepare_orca_resource_contract_recovery(
    source_input: Path,
    *,
    desired_ncores: int,
    execution_root_stem: str,
    existing_retry_inputs: list[Path],
    surface_context: dict[str, Any],
) -> tuple[Path, dict[str, Any]]:
    """Create an append-only execution segment restoring retained resources."""
    number = max(
        [
            int(match.group(1))
            for path in existing_retry_inputs
            if (match := re.search(r"-execution-retry-(\d{3})$", path.stem))
        ],
        default=0,
    ) + 1
    recovered_input = Path(source_input).parent / (
        f"{execution_root_stem}-execution-retry-{number:03d}.inp"
    )
    if recovered_input.exists():
        raise RuntimeError(f"Execution resource recovery input collision: {recovered_input}")
    source_text = Path(source_input).read_text(errors="replace")
    recovered_text = _orca_parallel_execution_input(source_text, desired_ncores)
    recovered_input.write_text(recovered_text)
    source_policy_path = Path(source_input).with_suffix(".execution-policy.json")
    source_policy = (
        json.loads(source_policy_path.read_text())
        if source_policy_path.is_file() else {}
    )
    policy = {
        "schema_version": 1,
        "kind": "storca_orca_execution_resource_contract_recovery",
        "strategy": "restore_retained_parallelism_after_serial_fallback",
        "source_input": str(source_input),
        "source_input_sha256": hashlib.sha256(Path(source_input).read_bytes()).hexdigest(),
        "prior_ncores": _orca_input_ncores(source_input),
        "effective_ncores": int(desired_ncores),
        "surface_context": surface_context.get("path"),
        "surface_context_sha256": surface_context.get("scientific_context_sha256"),
        "extra_environment": dict(source_policy.get("extra_environment") or {}),
        "scientific_model_changed": False,
        "execution_resources_changed": True,
        "claims_permitted": ["resume_same_orca_job_under_retained_resource_contract"],
        "claims_prohibited": ["mechanism", "rate", "lifetime"],
    }
    policy_path = recovered_input.with_suffix(".execution-policy.json")
    _write_json(policy_path, policy)
    evidence = {
        "schema_version": 1,
        "kind": "orca_execution_resource_contract_recovery",
        "source_input": str(source_input),
        "recovery_input": str(recovered_input),
        "execution_policy": str(policy_path),
        "scientific_input_changed": False,
        "execution_resources_changed": True,
        "reason": "serial_fallback_process_stopped_and_retained_parallel_contract_restored",
    }
    _write_json(recovered_input.with_suffix(".restart.json"), evidence)
    return recovered_input, evidence


def _neb_intermediate_image_indices(output: Path) -> list[int]:
    """Read ORCA's retained intermediate-image diagnostic without guessing chemistry."""
    output = Path(output)
    if not output.is_file():
        return []
    return _neb_intermediate_indices_from_text(output.read_text(errors="replace"))


def _neb_intermediate_indices_from_text(text: str) -> list[int]:
    """Extract ORCA intermediate image labels from retained or streamed text."""
    indices: list[int] = []
    for match in re.finditer(r"Possible\s+intermediate\s+minimum\s+found\s+at\s+image\(s\):\s*([^\n]+)", text, re.I):
        indices.extend(int(value) for value in re.findall(r"\d+", match.group(1)))
    return sorted(set(indices))


def _neb_angle_warning_image_indices(output: Path) -> list[int]:
    """Extract the adjacent image pair from ORCA's path-cusp warning.

    ORCA 6.1 can print an invalid first integer in this diagnostic.  The last
    two integers are the meaningful adjacent image labels; the malformed
    value is deliberately ignored.
    """
    output = Path(output)
    if not output.is_file():
        return []
    indices: set[int] = set()
    text = output.read_text(errors="replace")
    for match in _NEB_ANGLE_WARNING.finditer(text):
        values = re.findall(r"[-+]?\d+", match.group(1))
        if len(values) < 2:
            continue
        for value in values[-2:]:
            index = abs(int(value))
            if index > 0:
                indices.add(index)
    return sorted(indices)


def _neb_monotonic_endpoint_diagnostic(output: Path) -> dict:
    """Recognize a persistent endpoint maximum, not a transient initial band."""
    output = Path(output)
    if not output.is_file():
        return {"status": "not_detected", "valid": False, "trailing_iteration_count": 0}
    streak = 0
    maximum_streak = 0
    highest_energy_image = None
    last_iteration = None
    for line in output.read_text(errors="replace").splitlines():
        match = _NEB_MONOTONIC_WARNING.search(line)
        if match:
            streak += 1
            maximum_streak = max(maximum_streak, streak)
            last_iteration = int(match.group(1))
            highest_energy_image = int(match.group(2))
        elif _NEB_LBFGS_ITERATION.search(line):
            streak = 0
    valid = streak >= _NEB_MONOTONIC_STOP_STREAK
    return {
        "status": "persistent_monotonic_endpoint_maximum" if valid else "not_detected",
        "valid": valid,
        "trailing_iteration_count": streak,
        "maximum_iteration_count": maximum_streak,
        "last_iteration": last_iteration,
        "highest_energy_image": highest_energy_image,
        "required_consecutive_iterations": _NEB_MONOTONIC_STOP_STREAK,
        "interpretation": (
            "The highest-energy image remained an endpoint for the bounded warning streak; "
            "no interior transition-state candidate was available."
            if valid else None
        ),
    }


def _neb_persistent_endpoint_highest_image_diagnostic(
    output: Path, input_path: Path,
) -> dict[str, Any]:
    """Bound a band that never produces an interior climbing-image seed.

    This is an allocation diagnostic, not a no-barrier claim.  It requires a
    long trailing endpoint-HEI streak and repeated force excursions below the
    retained ORCA climbing-image switch threshold before redirecting the route
    to surface diagnosis.
    """
    output = Path(output)
    input_path = Path(input_path)
    if not output.is_file() or not input_path.is_file():
        return {"status": "not_detected", "valid": False}
    input_text = input_path.read_text(errors="replace")
    image_match = re.search(r"(?im)^\s*NImages\s+(\d+)\s*$", input_text)
    if image_match is None:
        return {"status": "not_detected", "valid": False, "reason": "nimages_unavailable"}
    endpoint_indices = {0, int(image_match.group(1)) + 1}
    text = output.read_text(errors="replace")
    threshold_matches = re.findall(
        r"Switch-on\s+CI\s+threshold\s+([-+0-9.Ee]+)", text, re.IGNORECASE,
    )
    if not threshold_matches:
        return {
            "status": "not_detected", "valid": False,
            "reason": "climbing_image_switch_threshold_unavailable",
        }
    switch_threshold = float(threshold_matches[-1])
    row_pattern = re.compile(
        r"LBFGS\s+(\d+)\s+(\d+)\s+([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
        re.IGNORECASE,
    )
    streak = 0
    low_force_count = 0
    highest_energy_image: int | None = None
    last_iteration: int | None = None
    maximum_force: float | None = None
    for match in row_pattern.finditer(text):
        iteration = int(match.group(1))
        hei = int(match.group(2))
        force = float(match.group(4))
        if hei in endpoint_indices:
            streak += 1
            low_force_count += int(force <= switch_threshold)
            highest_energy_image = hei
            last_iteration = iteration
            maximum_force = force
        else:
            streak = 0
            low_force_count = 0
            highest_energy_image = None
            last_iteration = iteration
            maximum_force = force
    valid = bool(
        streak >= _NEB_ENDPOINT_MAX_STOP_STREAK
        and low_force_count >= _NEB_ENDPOINT_MAX_LOW_FORCE_COUNT
    )
    return {
        "kind": "persistent_endpoint_highest_image",
        "status": "persistent_endpoint_highest_image" if valid else "not_detected",
        "valid": valid,
        "trailing_iteration_count": streak,
        "low_force_iteration_count": low_force_count,
        "required_trailing_iterations": _NEB_ENDPOINT_MAX_STOP_STREAK,
        "required_low_force_iterations": _NEB_ENDPOINT_MAX_LOW_FORCE_COUNT,
        "climbing_image_switch_threshold_Eh_per_bohr": switch_threshold,
        "highest_energy_image": highest_energy_image,
        "endpoint_image_indices": sorted(endpoint_indices),
        "last_iteration": last_iteration,
        "last_maximum_perpendicular_force_Eh_per_bohr": maximum_force,
        "interpretation": (
            "No interior highest-energy image emerged during the bounded streak; "
            "the retained band is redirected to surface diagnosis. This does not "
            "establish a barrierless path or a rate."
            if valid else None
        ),
    }


def _neb_submerged_climbing_image_diagnostic(output: Path) -> dict[str, Any]:
    """Detect a CI persistently below the higher fixed endpoint asymptote.

    Such an image may still be a real stationary point, but it is not the
    energetic bottleneck for motion toward the higher endpoint.  The result is
    therefore a method-selection diagnostic, never a barrier or rate claim.
    """
    output = Path(output)
    if not output.is_file():
        return {"status": "not_detected", "valid": False}
    text = output.read_text(errors="replace")
    endpoint_match = re.search(
        r"Properties\s+of\s+fixed\s+NEB\s+end\s+points:.*?"
        r"Reactant:.*?\bE\s+\.{2,}\s+([-+0-9.Ee]+).*?"
        r"Product:.*?\bE\s+\.{2,}\s+([-+0-9.Ee]+)",
        text, re.IGNORECASE | re.DOTALL,
    )
    ci_header = re.search(r"Optim\.\s+Iteration\s+CI\s+E\(CI\)-E\(0\)", text)
    if endpoint_match is None or ci_header is None:
        return {"status": "not_detected", "valid": False}
    reactant_energy = float(endpoint_match.group(1))
    product_energy = float(endpoint_match.group(2))
    product_relative = product_energy - reactant_energy
    higher_endpoint_relative = max(0.0, product_relative)
    row_pattern = re.compile(
        r"LBFGS\s+(\d+)\s+(\d+)\s+([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
        re.IGNORECASE,
    )
    streak = 0
    last: tuple[int, int, float, float] | None = None
    tolerance = 1.0e-4
    for match in row_pattern.finditer(text[ci_header.start():]):
        iteration = int(match.group(1))
        image = int(match.group(2))
        ci_relative = float(match.group(3))
        ci_force = float(match.group(7))
        if ci_relative < higher_endpoint_relative - tolerance:
            streak += 1
        else:
            streak = 0
        last = (iteration, image, ci_relative, ci_force)
    valid = streak >= _NEB_SUBMERGED_CI_STOP_STREAK
    return {
        "kind": "submerged_climbing_image_below_endpoint_asymptote",
        "status": (
            "persistent_submerged_climbing_image" if valid else "not_detected"
        ),
        "valid": valid,
        "trailing_iteration_count": streak,
        "required_trailing_iterations": _NEB_SUBMERGED_CI_STOP_STREAK,
        "reactant_energy_hartree": reactant_energy,
        "product_energy_hartree": product_energy,
        "product_relative_energy_hartree": product_relative,
        "higher_endpoint_relative_energy_hartree": higher_endpoint_relative,
        "last_iteration": last[0] if last else None,
        "climbing_image_index": last[1] if last else None,
        "climbing_image_relative_energy_hartree": last[2] if last else None,
        "climbing_image_max_force_Eh_per_bohr": last[3] if last else None,
        "interpretation": (
            "The climbing image remained below the higher fixed endpoint energy. "
            "This is a path/reference diagnostic: finite encounter energies are not "
            "automatically asymptotes. No mechanism or rate is established."
            if valid else None
        ),
    }


def _validate_neb_intermediates(
    route: RouteSpec,
    *,
    neb_output: Path,
    path_dir: Path,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    candidate_indices: list[int] | None = None,
    detection_source: str = "explicit_intermediate_warning",
    surface_resolution_policy: SurfaceResolutionPolicy | None = None,
) -> dict:
    """Optimize/frequency-check ORCA-flagged NEB intermediate images.

    This intentionally establishes only a retained intermediate/cusp minimum
    and two elementary segment tasks.  It does not assign a rate or infer a
    reaction family from a geometry; each segment still needs its own TS/IRC.
    """
    surface_resolution_policy = surface_resolution_policy or SurfaceResolutionPolicy()
    neb_output, path_dir = Path(neb_output), Path(path_dir)
    all_indices = sorted(set(candidate_indices or _neb_intermediate_image_indices(neb_output)))
    indices = all_indices[:surface_resolution_policy.maximum_intermediates]
    search_truncated = len(all_indices) > len(indices)
    if not indices:
        return {"status": "not_detected", "image_indices": [], "candidates": []}
    source_execution = {
        "kind": "orca_neb_intermediate_source_execution",
        "execution_stem": neb_output.stem,
        "output": str(neb_output.resolve()),
        "sha256": (
            hashlib.sha256(neb_output.read_bytes()).hexdigest()
            if neb_output.is_file() else None
        ),
        "detection_source": detection_source,
        "requested_image_indices": indices,
        "valid": bool(
            neb_output.is_file() and neb_output.suffix == ".out"
            and neb_output.parent.resolve() == path_dir.resolve()
        ),
    }
    expected_elements, _, _ = read_xyz(path_dir / "reactant.xyz")
    candidates = []
    for image_index in indices:
        # Retries/restarts retain their predecessors in the same directory.
        # An image index is meaningful only within the output that reported it.
        image = path_dir / f"{neb_output.stem}_im{image_index}.xyz"
        record: dict[str, Any] = {"image_index": image_index, "image_xyz": str(image)}
        if not source_execution["valid"]:
            record.update(status="source_output_unresolved")
            candidates.append(record)
            continue
        if not image.is_file():
            record.update(status="image_missing")
            candidates.append(record)
            continue
        try:
            source_binding = {
                "schema_version": 1,
                "kind": "orca_neb_intermediate_source_binding",
                "execution_stem": neb_output.stem,
                "image_index": image_index,
                "source_output": {
                    "path": source_execution["output"],
                    "sha256": source_execution["sha256"],
                },
                "source_image": {
                    "path": str(image.resolve()),
                    "sha256": hashlib.sha256(image.read_bytes()).hexdigest(),
                },
            }
            binding_id = hashlib.sha256(json.dumps(
                source_binding, sort_keys=True, separators=(",", ":"),
            ).encode()).hexdigest()
            validation_workspace = (
                path_dir / "intermediates" / neb_output.stem
                / f"image-{image_index:02d}-{binding_id}"
            )
            record.update(
                source_binding=source_binding,
                source_binding_id=binding_id,
                validation_workspace=str(validation_workspace),
            )
            elements, _, _ = read_xyz(image)
            if elements != expected_elements:
                raise RuntimeError("NEB intermediate image changed atom count or order")
            calculation = _optimize_endpoint(
                image, validation_workspace,
                charge=route.charge, multiplicity=route.multiplicity,
                ncores=ncores, method_keywords=method_keywords,
                timeout_seconds=timeout_seconds,
                surface_resolution_policy=surface_resolution_policy,
            )
            if (
                hashlib.sha256(image.read_bytes()).hexdigest()
                != source_binding["source_image"]["sha256"]
                or hashlib.sha256(neb_output.read_bytes()).hexdigest()
                != source_binding["source_output"]["sha256"]
            ):
                raise RuntimeError("NEB intermediate source image or output changed during validation")
            frequency = calculation["frequency_check"]
            optimized = Path(calculation["optimized_xyz"])
            optimized_elements, _, _ = read_xyz(optimized)
            if optimized_elements != expected_elements:
                raise RuntimeError("optimized intermediate changed atom count or order")
            recovery = calculation.get("stationary_point_recovery") or {}
            recovery_status = recovery.get("status")
            status = (
                "intermediate_optimization_unconverged"
                if recovery_status == "intermediate_optimization_unconverged" else
                "validated_intermediate_minimum" if frequency.get("IsMinimum") else
                "intermediate_nonminimum"
            )
            record.update(
                status=status,
                optimized_xyz=str(optimized),
                optimization_output=str(calculation["optimization"]["out"]),
                frequency_output=(
                    str(calculation["frequency"]["out"])
                    if calculation.get("frequency") else None
                ),
                frequency_check=frequency,
                stationary_point_recovery=recovery or None,
                minimum_recovery=calculation.get("minimum_recovery"),
            )
            if status == "validated_intermediate_minimum":
                identity = resolve_validated_intermediate_identity(
                    route,
                    optimized,
                    image_index=image_index,
                    parent_evidence={
                        "detection_source": detection_source,
                        "source_binding": source_binding,
                        "source_binding_id": binding_id,
                        "optimization_output": str(calculation["optimization"]["out"]),
                        "frequency_output": (
                            str(calculation["frequency"]["out"])
                            if calculation.get("frequency") else None
                        ),
                    },
                )
                record["intermediate_identity"] = identity
                _write_json(
                    optimized.parent / "intermediate-identity.json", identity,
                )
        except Exception as error:
            if isinstance(error, (SurfaceResolutionBudgetExhausted, TransferControlError)):
                raise
            record.update(status="intermediate_validation_failed", failure_reason=str(error))
        candidates.append(record)
    validated = [item for item in candidates if item.get("status") == "validated_intermediate_minimum"]
    identity_resolved = [
        item for item in validated
        if (item.get("intermediate_identity") or {}).get("status") == "resolved"
    ]
    endpoint_returns = [
        item for item in validated
        if (item.get("intermediate_identity") or {}).get("status") == "endpoint_return"
    ]
    identity_unresolved = [
        item for item in validated
        if item not in identity_resolved and item not in endpoint_returns
    ]
    return {
        "schema_version": 1,
        "kind": "orca_neb_intermediate_evidence",
        "status": (
            "surface_resolution_budget_exhausted" if search_truncated else
            "validated_intermediate_detected" if identity_resolved else
            "intermediate_identity_unresolved" if identity_unresolved else
            "intermediate_warning_dismissed_endpoint_return" if endpoint_returns else
            "intermediate_detected_unvalidated"
        ),
        "detection_source": detection_source,
        "source_execution": source_execution,
        "image_indices": indices,
        "all_reported_image_indices": all_indices,
        "candidate_search_truncated": search_truncated,
        "maximum_intermediates": surface_resolution_policy.maximum_intermediates,
        "candidates": candidates,
        "validated_candidate_count": len(validated),
        "identity_resolved_candidate_count": len(identity_resolved),
        "identity_unresolved_candidate_count": len(identity_unresolved),
        "endpoint_return_candidate_count": len(endpoint_returns),
        "segmented_verification": {
            "status": "prepared" if identity_resolved else "not_prepared",
            "segments": [
                segment
                for intermediate in identity_resolved
                for segment in (
                    {
                        "from": "declared_reactants",
                        "to": intermediate["optimized_xyz"],
                        "native_candidate": intermediate["intermediate_identity"]["native_candidates"][0],
                    },
                    {
                        "from": intermediate["optimized_xyz"],
                        "to": "declared_products",
                        "native_candidate": intermediate["intermediate_identity"]["native_candidates"][1],
                    },
                )
            ],
            "requirements": [
                "optimize/validate each segment transition state",
                "require one relevant imaginary frequency per transition state",
                "run IRC or endpoint-connection validation for each segment",
                "derive the overall rate from the validated elementary sequence",
            ],
        },
        "interpretation": (
            "An ORCA-flagged image that relaxes to a declared endpoint is retained as a dismissed warning, "
            "not a distinct intermediate. Identity-resolved distinct minima are split into two explicitly "
            "retained segment-verification tasks; unresolved identities remain fail-closed."
        ),
    }


def _run_segmented_intermediate_nebs(
    route: RouteSpec,
    *,
    reactant_xyz: Path,
    product_xyz: Path,
    intermediate_evidence: dict,
    path_dir: Path,
    separated_fragments: bool,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    neb_images: int,
    surface_resolution_policy: SurfaceResolutionPolicy | None = None,
    segmentation_depth: int = 1,
    visited_endpoint_pair_hashes: set[str] | None = None,
    parent_endpoint_stationary: dict[str, list[dict]] | None = None,
    require_reached_minima: bool = False,
) -> dict:
    """Verify identity-resolved elementary segments under explicit recursion bounds."""
    surface_resolution_policy = surface_resolution_policy or SurfaceResolutionPolicy()
    visited = set(visited_endpoint_pair_hashes or ())
    if segmentation_depth > surface_resolution_policy.maximum_segmentation_depth:
        return {
            "status": "surface_resolution_budget_exhausted",
            "reason": "maximum_segmentation_depth_exceeded",
            "segmentation_depth": segmentation_depth,
            "maximum_segmentation_depth": surface_resolution_policy.maximum_segmentation_depth,
            "segments": [],
            "rate_claim_allowed": False,
        }
    candidates = [
        item for item in intermediate_evidence.get("candidates", [])
        if item.get("status") == "validated_intermediate_minimum"
        and item.get("optimized_xyz")
        and (item.get("intermediate_identity") or {}).get("status") == "resolved"
    ]
    if not candidates:
        identity_unresolved = any(
            item.get("status") == "validated_intermediate_minimum"
            and (item.get("intermediate_identity") or {}).get("status")
                == "intermediate_identity_unresolved"
            for item in intermediate_evidence.get("candidates", [])
        )
        return {
            "status": (
                "intermediate_identity_unresolved" if identity_unresolved else "not_prepared"
            ),
            "segments": [],
            "reason": (
                "validated_minimum_has_no_defensible_molecular_identity"
                if identity_unresolved else "no_validated_intermediate_geometry"
            ),
            "rate_claim_allowed": False,
        }
    segment_records: list[dict] = []
    # A shorter band is sufficient for each local segment and avoids repeating
    # the full expensive parent path resolution.
    segment_image_count = max(3, min(int(neb_images), 6))
    for intermediate in candidates:
        intermediate_xyz = Path(intermediate["optimized_xyz"])
        native_candidates = list(
            (intermediate.get("intermediate_identity") or {}).get("native_candidates") or []
        )
        if len(native_candidates) != 2:
            segment_records.append({
                "label": "intermediate-segmentation",
                "status": "intermediate_identity_unresolved",
                "reason": "resolved identity did not retain exactly two native segment candidates",
                "rate_claim_allowed": False,
            })
            continue
        endpoints = [
            (
                "reactants-to-intermediate", Path(reactant_xyz), intermediate_xyz,
                len(route.reactant_smiles), 1,
                list(route.reactant_multiplicities), [route.multiplicity], native_candidates[0],
            ),
            (
                "intermediate-to-products", intermediate_xyz, Path(product_xyz),
                1, len(route.product_smiles),
                [route.multiplicity], list(route.product_multiplicities), native_candidates[1],
            ),
        ]
        for (label, source, target, source_fragment_count, target_fragment_count,
             source_fragment_multiplicities, target_fragment_multiplicities,
             native_candidate) in endpoints:
            segment_dir = path_dir / "segments" / label
            segment_dir.mkdir(parents=True, exist_ok=True)
            segment_reactant = segment_dir / "reactant.xyz"
            segment_product = segment_dir / "product.xyz"
            shutil.copy2(source, segment_reactant)
            shutil.copy2(target, segment_product)
            channel_kind = (
                "association" if source_fragment_count > target_fragment_count else
                "dissociation" if source_fragment_count < target_fragment_count else
                "bound_rearrangement"
            )
            endpoint_pair_hash = str(native_candidate.get("candidate_id") or "")
            if not endpoint_pair_hash:
                segment_records.append({
                    "label": label,
                    "status": "intermediate_identity_unresolved",
                    "reason": "native segment candidate has no stable endpoint-pair identity",
                    "native_candidate": native_candidate,
                    "rate_claim_allowed": False,
                })
                continue
            if endpoint_pair_hash in visited:
                segment_records.append({
                    "label": label,
                    "status": "segment_path_unresolved",
                    "reason": "visited_endpoint_pair_recursion_prevented",
                    "endpoint_pair_hash": endpoint_pair_hash,
                    "segmentation_depth": segmentation_depth,
                    "native_candidate": native_candidate,
                    "rate_claim_allowed": False,
                })
                continue
            visited.add(endpoint_pair_hash)
            if source_fragment_count > 1 or target_fragment_count > 1:
                segment_records.append({
                    "label": label,
                    "status": "channel_specific_solver_required",
                    "channel_kind": channel_kind,
                    "source_fragment_count": source_fragment_count,
                    "target_fragment_count": target_fragment_count,
                    "source_fragment_multiplicities": source_fragment_multiplicities,
                    "target_fragment_multiplicities": target_fragment_multiplicities,
                    "selected_total_multiplicity": route.multiplicity,
                    "native_candidate": native_candidate,
                    "reactant_xyz": str(segment_reactant),
                    "product_xyz": str(segment_product),
                    "required_solver": (
                        "fragment_spin_validated_capture_recrossing_model"
                        if channel_kind == "association" else
                        "fragment_spin_validated_constrained_dissociation_scan"
                    ),
                    "reason": (
                        "A separated-channel endpoint preserves only the selected total multiplicity. "
                        "An unconstrained endpoint optimization or fixed-end NEB cannot prove that "
                        "the declared fragment-local spin states and connectivity were retained."
                    ),
                    "rate_claim_allowed": False,
                })
                continue
            try:
                candidate = ReactionCandidate.from_dict(native_candidate)
                segment_route_dict = candidate.as_route_dict(
                    parent_smiles=route.parent_smiles,
                    target_label=(candidate.reactant_labels[0] if candidate.reactant_labels else "segment"),
                )
                segment_route = route_spec_from_candidate_route(
                    route.parent_smiles, segment_route_dict,
                    stage=f"segment-depth-{segmentation_depth}",
                )
                segment_mapping = resolve_route_atom_mapping(segment_route)
                if not segment_mapping.get("valid"):
                    raise ValueError(
                        "segment atom mapping unresolved: " + str(segment_mapping.get("status"))
                    )
                endpoint_records = dict(parent_endpoint_stationary or {})
                def bind_parent_records(
                    records: list[dict], smiles_values: tuple[str, ...],
                    labels: tuple[str, ...], charges: tuple[int, ...],
                    multiplicities: tuple[int, ...],
                ) -> list[dict]:
                    if len(records) != len(smiles_values):
                        return records
                    bound: list[dict] = []
                    for index, (record, smiles) in enumerate(zip(records, smiles_values)):
                        item = dict(record)
                        label_value = labels[index] if len(labels) == len(smiles_values) else item.get("label")
                        item.update(
                            label=label_value,
                            original_label=label_value,
                            smiles=smiles,
                            charge=(
                                charges[index] if len(charges) == len(smiles_values)
                                else item.get("charge")
                            ),
                            multiplicity=(
                                multiplicities[index]
                                if len(multiplicities) == len(smiles_values)
                                else item.get("multiplicity")
                            ),
                        )
                        bound.append(item)
                    return bound

                parent_reactants = bind_parent_records(
                    list(endpoint_records.get("reactants") or []),
                    candidate.reactants,
                    candidate.reactant_labels,
                    candidate.reactant_charges,
                    candidate.reactant_multiplicities,
                )
                parent_products = bind_parent_records(
                    list(endpoint_records.get("products") or []),
                    candidate.products,
                    candidate.product_labels,
                    candidate.product_charges,
                    candidate.product_multiplicities,
                )
                intermediate_record = {
                    "label": (
                        candidate.product_labels[0]
                        if label == "reactants-to-intermediate" else candidate.reactant_labels[0]
                    ),
                    "original_label": (
                        candidate.product_labels[0]
                        if label == "reactants-to-intermediate" else candidate.reactant_labels[0]
                    ),
                    "smiles": (
                        candidate.products[0]
                        if label == "reactants-to-intermediate" else candidate.reactants[0]
                    ),
                    "charge": route.charge,
                    "multiplicity": route.multiplicity,
                    "status": "validated_minimum",
                    "optimized_xyz": str(intermediate_xyz),
                    "frequency_output": intermediate.get("frequency_output"),
                    "intermediate_evidence_id": (
                        (intermediate.get("intermediate_identity") or {}).get("evidence_id")
                    ),
                }
                if label == "reactants-to-intermediate":
                    segment_reactant_records = parent_reactants
                    segment_product_records = [intermediate_record]
                else:
                    segment_reactant_records = [intermediate_record]
                    segment_product_records = parent_products
                rate_inputs_ready = bool(
                    len(segment_reactant_records) == len(candidate.reactants)
                    and len(segment_product_records) == len(candidate.products)
                    and all(
                        item.get("status") == "validated_minimum"
                        and item.get("frequency_output")
                        and Path(str(item["frequency_output"])).is_file()
                        for item in (*segment_reactant_records, *segment_product_records)
                    )
                )
            except Exception as error:
                if isinstance(error, (SurfaceResolutionBudgetExhausted, TransferControlError)):
                    raise
                segment_records.append({
                    "label": label,
                    "status": "segment_path_unresolved",
                    "reason": f"segment_route_construction_failed: {type(error).__name__}: {error}",
                    "endpoint_pair_hash": endpoint_pair_hash,
                    "native_candidate": native_candidate,
                    "rate_claim_allowed": False,
                })
                continue
            neb_input = create_orca_neb_ts_input(
                segment_reactant, segment_product,
                charge=route.charge, multiplicity=route.multiplicity,
                label="neb-ts", ncores=ncores, nimages=segment_image_count,
                method_keywords=method_keywords,
                preopt_ends=not separated_fragments,
            )
            try:
                execution = _run_or_resume(neb_input, timeout_seconds=timeout_seconds)
                execution_stem = Path(execution.get("input") or neb_input).stem
                trajectory = _locate_neb_trajectory(segment_dir, execution_stem)
                ts = _locate_neb_ts(segment_dir, execution_stem)
                execution_status = execution.get("status")
                if execution_status == "completed" and ts is not None:
                    intended = None
                    if require_reached_minima:
                        from .route_geometry import _build_endpoint_side, _mapped_bond_changes
                        edits = _mapped_bond_changes(
                            _build_endpoint_side(segment_route.reactant_smiles),
                            _build_endpoint_side(segment_route.product_smiles),
                            dict(segment_mapping["all_atom_mapping"]))
                        intended = [{"atom_indices": pair, "kind": "mapped_segment_bond_change"}
                                    for pair in edits["formed_bonds"] + edits["broken_bonds"]]
                        intended += [{"atom_indices": item["atom_indices"], "kind": "bond_order_change"}
                                     for item in edits["bond_order_changes"]]
                    verification = _verify_ts_seed(
                        segment_route,
                        ts,
                        segment_dir / "verification",
                        reactant_xyz=segment_reactant,
                        product_xyz=segment_product,
                        mapping_result=segment_mapping,
                        separated_fragments=False,
                        ncores=ncores,
                        method_keywords=method_keywords,
                        timeout_seconds=timeout_seconds,
                        surface_resolution_policy=surface_resolution_policy,
                        **({"intended_coordinates": intended} if require_reached_minima else {}),
                    )
                    if require_reached_minima and verification.get("status") == "verified_transition_state_path":
                        verification.update(_refine_reached_irc_stationary_points(
                            segment_route, verification, segment_dir / "reached-stationary-points",
                            segment_mapping, ncores=ncores, method_keywords=method_keywords,
                            timeout_seconds=timeout_seconds, surface_resolution_policy=surface_resolution_policy,
                            reference_endpoints={"reactants": segment_reactant, "products": segment_product}))
                        reached = verification.get("stationary_points") or {}
                        if not reached.get("valid"):
                            verification["status"] = "reached_minimum_validation_failed"
                        else:
                            segment_reactant_records = reached["reactants"]
                            segment_product_records = reached["products"]
                    verified_status = str(verification.get("status") or "segment_path_unresolved")
                    segment_status = (
                        "verified_transition_state_path"
                        if verified_status == "verified_transition_state_path" else
                        verified_status
                        if verified_status in {
                            "ts_mode_mismatch", "irc_endpoint_validation_failed", "irc_path_incomplete",
                            "ts_frequency_validation_failed", "ts_optimization_unconverged",
                        } else
                        "segment_path_unresolved"
                    )
                    transition_state = verification.get("transition_state") or {}
                    rate_inputs_ready = bool(
                        rate_inputs_ready
                        and transition_state.get("frequency_output")
                        and Path(str(transition_state["frequency_output"])).is_file()
                    )
                elif execution_status == "completed_no_interior_barrier":
                    verification = {
                        "status": "barrierless_rate_model_required",
                        "reason": "segment_neb_has_no_conventional_interior_saddle",
                        "method_escalation": "barrierless_or_variational_rate_model",
                    }
                    segment_status = "barrierless_rate_model_required"
                else:
                    verification = {
                        "status": "segment_path_unresolved",
                        "reason": (
                            "segment_neb_produced_no_ts_candidate"
                            if ts is None else f"segment_neb_status_{execution_status}"
                        ),
                    }
                    segment_status = "segment_path_unresolved"
                segment_records.append({
                    "label": label,
                    "status": segment_status,
                    "channel_kind": channel_kind,
                    "native_candidate": native_candidate,
                    "segment_route": segment_route.as_dict(),
                    "segment_mapping": segment_mapping,
                    "endpoint_pair_hash": endpoint_pair_hash,
                    "segmentation_depth": segmentation_depth,
                    "input": str(neb_input),
                    "output": execution.get("output"),
                    "execution": execution,
                    "trajectory": str(trajectory) if trajectory else None,
                    "ts_candidate": str(ts) if ts else None,
                    "verification": verification,
                    "rate_input_status": (
                        "ready" if segment_status == "verified_transition_state_path"
                        and rate_inputs_ready else "unavailable"
                    ),
                    "arkane_inputs": (
                        {
                            "reactants": segment_reactant_records,
                            "products": segment_product_records,
                            "transition_state": verification.get("transition_state"),
                        }
                        if segment_status == "verified_transition_state_path"
                        and rate_inputs_ready else None
                    ),
                    "rate_claim_allowed": False,
                })
            except Exception as error:
                if isinstance(error, (SurfaceResolutionBudgetExhausted, TransferControlError)):
                    raise
                segment_records.append({
                    "label": label,
                    "status": "failed",
                    "input": str(neb_input),
                    "native_candidate": native_candidate,
                    "failure_reason": f"{type(error).__name__}: {error}",
                    "rate_claim_allowed": False,
                })
    for record in segment_records:
        segment_root = (
            Path(str(record["input"])).parent
            if record.get("input") else path_dir / "segments" / str(record.get("label") or "unresolved")
        )
        input_paths = [
            Path(str(value)) for value in (
                record.get("input"),
                record.get("reactant_xyz"),
                record.get("product_xyz"),
            ) if value
        ]
        record["surface_resolution_attempt_evidence"] = _write_surface_attempt_evidence(
            segment_root,
            strategy="segmented_path_verification",
            status=str(record.get("status") or "segment_path_unresolved"),
            result=record,
            parameters={
                "segment_label": record.get("label"),
                "channel_kind": record.get("channel_kind"),
                "segmentation_depth": record.get("segmentation_depth", segmentation_depth),
                "endpoint_pair_hash": record.get("endpoint_pair_hash"),
            },
            input_paths=input_paths,
            failure_reason=record.get("failure_reason") or record.get("reason"),
            claims_permitted=(
                ("verified_transition_state_path",)
                if record.get("status") == "verified_transition_state_path" else ()
            ),
        )
    completed = sum(
        item.get("status") == "verified_transition_state_path" for item in segment_records
    )
    rate_ready = sum(
        item.get("status") == "verified_transition_state_path"
        and item.get("rate_input_status") == "ready"
        for item in segment_records
    )
    deferred = sum(item.get("status") == "channel_specific_solver_required" for item in segment_records)
    barrierless = sum(item.get("status") == "barrierless_rate_model_required" for item in segment_records)
    budget_exhausted = any(
        item.get("status") == "surface_resolution_budget_exhausted" for item in segment_records
    )
    return {
        "status": (
            "verified_segments_rate_generation_required"
            if segment_records and completed == len(segment_records)
            and rate_ready == len(segment_records) else
            "segment_rate_model_unavailable"
            if segment_records and completed == len(segment_records) else
            "surface_resolution_budget_exhausted" if budget_exhausted else
            "barrierless_rate_model_required" if barrierless else
            "segment_rate_model_unavailable" if deferred else
            "segment_path_unresolved"
        ),
        "segment_count": len(segment_records),
        "completed_segment_count": completed,
        "rate_ready_segment_count": rate_ready,
        "channel_specific_segment_count": deferred,
        "barrierless_segment_count": barrierless,
        "segmentation_depth": segmentation_depth,
        "maximum_segmentation_depth": surface_resolution_policy.maximum_segmentation_depth,
        "visited_endpoint_pair_hashes": sorted(visited),
        "image_count_per_segment": segment_image_count,
        "segments": segment_records,
        "rate_claim_allowed": False,
        "requirements": [
            "validate each segment endpoint connectivity",
            "optimize each segment TS and require one relevant imaginary frequency",
            "run IRC or endpoint-connection validation for each segment",
            "derive the overall rate from the validated elementary sequence",
        ],
    }


def _orca_job_input_evidence(input_path: Path) -> tuple[Any, list[str]]:
    """Hash the exact input and referenced evidence without dispatch or writes."""
    input_path = Path(input_path)
    input_text = input_path.read_text(errors="replace")
    referenced_names = set(re.findall(
        r'(?:Product|Hess_Filename|Restart_ALLXYZFile)\s+"([^"]+)"', input_text,
        re.IGNORECASE,
    ))
    xyzfile = re.search(r"\*\s+xyzfile\s+[-+]?\d+\s+\d+\s+(\S+)", input_text, re.IGNORECASE)
    if xyzfile:
        referenced_names.add(xyzfile.group(1).strip('"'))
    digest = hashlib.sha256(input_path.read_bytes())
    referenced = []
    for name in sorted(referenced_names):
        path = input_path.parent / name
        if not path.is_file():
            raise FileNotFoundError(f"ORCA input references a missing artifact: {path}")
        digest.update(name.encode())
        digest.update(path.read_bytes())
        referenced.append(str(path))
    return digest, referenced


def _validate_retained_orca_job_input(input_path: Path) -> None:
    """Read-only cache check; absent legacy contracts receive no new proof."""
    contract_path = Path(input_path).with_suffix(".contract.json")
    if not contract_path.is_file():
        if active_native_settings() is not None:
            raise RuntimeError("Scoped native reuse requires its retained complete job contract")
        return
    contract = json.loads(contract_path.read_text())
    settings = active_native_settings()
    if settings is not None and contract.get("native_settings_sha256") != settings.sha256:
        raise RuntimeError("Retained ORCA input lacks the active complete native settings identity")
    digest, _ = _orca_job_input_evidence(input_path)
    if contract.get("sha256") != digest.hexdigest():
        raise RuntimeError(f"Retained ORCA input or referenced evidence changed: {input_path}")
    from .rate_transfer_control import active_execution_control
    control = active_execution_control()
    if control is not None:
        control.validate_reuse(input_path, digest.hexdigest())


def _run_or_resume(
    input_path: Path,
    *,
    timeout_seconds: float | None,
    ignore_intermediate_warning_indices: Iterable[int] = (),
) -> dict:
    input_path = Path(input_path)
    settings = active_native_settings()
    if settings is not None:
        retained_text = input_path.read_text()
        complete_text = apply_native_job_settings(retained_text, settings)
        if complete_text != retained_text:
            if input_path.with_suffix(".contract.json").exists() or input_path.with_suffix(".out").exists():
                raise RuntimeError("Retained ORCA input does not contain the complete frozen native settings")
            input_path.write_text(complete_text)
    native_settings_identity = ({"native_settings_sha256": settings.sha256} if settings is not None else {})
    from .rate_transfer_control import active_execution_control
    execution_control = active_execution_control()
    ignored_intermediate_indices = {
        int(value) for value in ignore_intermediate_warning_indices
    }
    output = input_path.with_suffix(".out")
    input_text = input_path.read_text(errors="replace")
    is_neb = "NEB-TS" in input_text.upper()
    digest, referenced = _orca_job_input_evidence(input_path)
    contract_path = input_path.with_suffix(".contract.json")
    contract = json.loads(contract_path.read_text()) if contract_path.is_file() else {}
    if settings is not None and contract and contract.get("native_settings_sha256") != settings.sha256:
        raise RuntimeError("Retained ORCA contract has a different native settings identity")
    if contract and contract.get("sha256") != digest.hexdigest():
        raise RuntimeError(
            f"Retained ORCA input or referenced evidence changed; refusing to replace its job contract: {input_path}"
        )
    if execution_control is not None:
        execution_control.validate_reuse(input_path, digest.hexdigest())
    no_barrier_outcome = _orca_neb_no_barrier_outcome(output)
    existing_intermediates = _neb_intermediate_image_indices(output)
    existing_angle_warning = _neb_angle_warning_image_indices(output)
    existing_monotonic_warning = _neb_monotonic_endpoint_diagnostic(output)
    existing_endpoint_highest = _neb_persistent_endpoint_highest_image_diagnostic(
        output, input_path,
    )
    existing_submerged_ci = _neb_submerged_climbing_image_diagnostic(output)
    legacy_retainable_completion = not contract and (
        _normal_orca_output(output) or no_barrier_outcome["valid"]
        or bool(existing_intermediates) or bool(existing_angle_warning)
        or existing_monotonic_warning["valid"] or existing_endpoint_highest["valid"]
        or existing_submerged_ci["valid"]
    )
    if legacy_retainable_completion:
        _write_json(contract_path, {
            "schema_version": 1,
            "kind": "orca_job_input_contract",
            "sha256": digest.hexdigest(),
            "input": str(input_path),
            "referenced_files": referenced,
            "status": "retained_legacy_completion",
            **native_settings_identity,
        })
        contract = {"sha256": digest.hexdigest()}
    if contract.get("sha256") == digest.hexdigest():
        if _normal_orca_output(output):
            return {"status": "completed", "resumed": True, "input": str(input_path), "output": str(output)}
        if no_barrier_outcome["valid"]:
            return {
                "status": "completed_no_interior_barrier", "resumed": True,
                "input": str(input_path), "output": str(output),
                "orca_neb_outcome": no_barrier_outcome,
            }
        if existing_intermediates and not set(existing_intermediates).issubset(
            ignored_intermediate_indices
        ):
            return {
                "status": "completed_intermediate_detected", "resumed": True,
                "input": str(input_path), "output": str(output),
                "intermediate_image_indices": existing_intermediates,
            }
        if existing_angle_warning:
            return {
                "status": "completed_path_quality_warning", "resumed": True,
                "input": str(input_path), "output": str(output),
                "path_quality_warning": {
                    "kind": "maximum_angle_between_images",
                    "image_indices": existing_angle_warning,
                },
            }
        if existing_monotonic_warning["valid"]:
            return {
                "status": "completed_path_quality_warning", "resumed": True,
                "input": str(input_path), "output": str(output),
                "path_quality_warning": existing_monotonic_warning,
            }
        if existing_endpoint_highest["valid"]:
            return {
                "status": "completed_path_quality_warning", "resumed": True,
                "input": str(input_path), "output": str(output),
                "path_quality_warning": existing_endpoint_highest,
            }
        if existing_submerged_ci["valid"]:
            return {
                "status": "completed_path_quality_warning", "resumed": True,
                "input": str(input_path), "output": str(output),
                "path_quality_warning": existing_submerged_ci,
            }

    # A retained native restart may have failed before any iteration because
    # the old generator inherited incompatible endpoint preoptimization.
    # Gate this before any generic execution retry: inability to establish
    # stopped, immutable ownership is a refusal, not permission to rerun.
    if is_neb and re.search(r"(?im)^\s*Restart_ALLXYZFile\b", input_text):
        from src.orca_input_repair import (
            inspect_neb_restart_preopt_parser_failure,
            prepare_neb_restart_preopt_input_repair,
        )
        restart_preopt_evidence = inspect_neb_restart_preopt_parser_failure(input_path)
        if restart_preopt_evidence is not None:
            repaired_input = prepare_neb_restart_preopt_input_repair(input_path)
            result = _run_or_resume(
                repaired_input, timeout_seconds=timeout_seconds,
                ignore_intermediate_warning_indices=ignored_intermediate_indices,
            )
            result["restart_preopt_input_repair"] = str(repaired_input.with_suffix(".syntax-repair.json"))
            return result
        if output.is_file() and (
            "Preoptimization of end points requested. This option is not possible with restart."
            in output.read_text(errors="replace")
        ):
            raise RuntimeError(
                "Retained ORCA restart/preoptimization conflict lacks a stopped, immutable "
                f"source-bound failure proof; refusing execution: {input_path}"
            )

    # A launcher denial is neither chemistry nor a reusable completion.
    # Preserve it under the original basename and create an exact-input,
    # append-only execution retry rather than overwriting its evidence.
    execution_root_stem = re.sub(r"-execution-retry-\d{3}$", "", input_path.stem)
    execution_retry_inputs = sorted(
        input_path.parent.glob(f"{execution_root_stem}-execution-retry-*.inp"),
        key=lambda path: int(re.search(r"-execution-retry-(\d{3})$", path.stem).group(1))
        if re.search(r"-execution-retry-(\d{3})$", path.stem) else -1,
    )
    if execution_retry_inputs and input_path.resolve() != execution_retry_inputs[-1].resolve():
        latest_retry = execution_retry_inputs[-1]
        retained_surface_context = _surface_context(input_path)
        desired_ncores = int(
            ((retained_surface_context or {}).get("electronic_structure") or {}).get(
                "ncores", _orca_input_ncores(latest_retry),
            )
        )
        resource_recovery: dict[str, Any] | None = None
        if (
            desired_ncores != _orca_input_ncores(latest_retry)
            and not _normal_orca_output(latest_retry.with_suffix(".out"))
            and _detached_orca_execution_is_stopped(latest_retry)
        ):
            latest_retry, resource_recovery = _prepare_orca_resource_contract_recovery(
                latest_retry,
                desired_ncores=desired_ncores,
                execution_root_stem=execution_root_stem,
                existing_retry_inputs=execution_retry_inputs,
                surface_context=retained_surface_context or {},
            )
        result = _run_or_resume(
            latest_retry, timeout_seconds=timeout_seconds,
            ignore_intermediate_warning_indices=ignored_intermediate_indices,
        )
        if resource_recovery is not None:
            result["execution_resource_contract_recovery"] = resource_recovery
        result.setdefault("execution_retry_chain_root_input", str(input_path))
        return result
    startup_failure = _orca_startup_infrastructure_failure(output)
    stopped_interruption = _detached_orca_execution_is_stopped(input_path)
    if contract.get("sha256") == digest.hexdigest() \
            and (startup_failure or stopped_interruption):
        if startup_failure:
            retry_input, retry_evidence = _prepare_orca_execution_retry(
                input_path,
                output=output,
                existing_retry_inputs=execution_retry_inputs,
            )
        else:
            number = len(execution_retry_inputs) + 1
            retry_input = input_path.parent / f"{execution_root_stem}-execution-retry-{number:03d}.inp"
            if retry_input.exists():
                raise RuntimeError(f"Execution retry input collision: {retry_input}")
            retry_input.write_text(input_text)
            retry_evidence = {
                "schema_version": 1,
                "kind": "orca_interrupted_execution_retry",
                "source_input": str(input_path),
                "source_output": str(output),
                "source_output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
                "retry_input": str(retry_input),
                "scientific_input_changed": False,
                "reason": "retained_interrupted_execution_confirmed_stopped",
            }
            _write_json(retry_input.with_suffix(".restart.json"), retry_evidence)
        result = _run_or_resume(
            retry_input, timeout_seconds=timeout_seconds,
            ignore_intermediate_warning_indices=ignored_intermediate_indices,
        )
        result["execution_startup_retry"] = retry_evidence
        return result

    # Prefer the newest append-only restart segment if one already exists.
    # This makes a high-level workflow resume converge on the actual ORCA job
    # basename and prevents silently overwriting an interrupted segment.
    root_stem = re.sub(r"-restart-\d{3}$", "", input_path.stem)
    restart_inputs = sorted(
        input_path.parent.glob(f"{root_stem}-restart-*.inp"),
        key=lambda path: int(re.search(r"-restart-(\d{3})$", path.stem).group(1))
        if re.search(r"-restart-(\d{3})$", path.stem) else -1,
    ) if is_neb else []
    if restart_inputs and input_path.resolve() != restart_inputs[-1].resolve():
        result = _run_or_resume(
            restart_inputs[-1], timeout_seconds=timeout_seconds,
            ignore_intermediate_warning_indices=ignored_intermediate_indices,
        )
        result.setdefault("restart_chain_root_input", str(input_path))
        return result

    # ORCA documents the MEP .allxyz file as the NEB restart checkpoint.  A
    # partial band is reusable only after its entire artifact family has been
    # quiet long enough to rule out a detached writer.
    has_completed_iteration = output.is_file() and any(
        _NEB_LBFGS_ITERATION.search(line) for line in output.read_text(errors="replace").splitlines()
    )
    if is_neb and has_completed_iteration:
        process_lease = input_path.with_suffix(".process.json")
        confirmed_stopped_at: float | None = None
        if process_lease.is_file():
            try:
                lease = json.loads(process_lease.read_text())
                pid = int(lease.get("pid"))
                process_group_id = int(lease.get("process_group_id") or pid)
                live_status = lease.get("status") in {
                    "running", "runner_detached_process_may_be_active", "terminated",
                }
                if live_status and pid > 0:
                    alive = True
                    try:
                        os.kill(pid, 0)
                    except ProcessLookupError:
                        try:
                            os.killpg(process_group_id, 0)
                        except ProcessLookupError:
                            alive = False
                    if alive:
                        raise RuntimeError(
                            "ORCA NEB process lease is still live; refusing a duplicate restart "
                            f"(pid={pid}, lease={process_lease})"
                        )
                    ended = lease.get("ended_at_epoch")
                    if isinstance(ended, (int, float)):
                        confirmed_stopped_at = float(ended)
            except ProcessLookupError:
                pass
            except PermissionError as error:
                raise RuntimeError(
                    f"ORCA NEB process lease cannot be disproved: {process_lease}"
                ) from error
            except (OSError, TypeError, ValueError, json.JSONDecodeError):
                # A corrupt lease cannot prove liveness; the conservative
                # artifact quiet gate below remains mandatory.
                pass
        # Input generation is deterministic but may rewrite the identical
        # ``.inp`` immediately before this check.  Inputs and provenance JSON
        # are not engine heartbeats; only ORCA runtime artifacts count here.
        related = [
            path for path in input_path.parent.glob(f"{input_path.stem}*")
            if path.is_file()
            and path != input_path
            and not path.name.endswith((".contract.json", ".restart.json", ".process.json"))
        ]
        if not related:
            raise RuntimeError("Interrupted ORCA NEB has no retained runtime artifacts")
        newest = max(related, key=lambda path: path.stat().st_mtime)
        age = max(0.0, time.time() - newest.stat().st_mtime)
        stopped_lease_covers_latest = (
            confirmed_stopped_at is not None
            and confirmed_stopped_at + 30.0 >= newest.stat().st_mtime
        )
        if age < _ORCA_RESTART_QUIET_SECONDS and not stopped_lease_covers_latest:
            raise RuntimeError(
                "ORCA NEB artifacts still have a fresh heartbeat; refusing a duplicate "
                f"restart ({newest} is {age:.1f} seconds old)"
            )
        band = input_path.parent / f"{input_path.stem}_MEP.allxyz"
        if not contract_path.exists():
            _write_json(contract_path, {
                "schema_version": 1, "kind": "orca_job_input_contract",
                "sha256": digest.hexdigest(), "input": str(input_path),
                "referenced_files": referenced, "status": "retained_interrupted_neb_checkpoint",
            })
        restart_input, restart_evidence = _prepare_neb_restart_input(input_path, band)
        result = _run_or_resume(
            restart_input, timeout_seconds=timeout_seconds,
            ignore_intermediate_warning_indices=ignored_intermediate_indices,
        )
        result["neb_checkpoint_restart"] = restart_evidence
        if ignored_intermediate_indices:
            result["validated_endpoint_return_warning_policy"] = {
                "status": "active",
                "ignored_image_indices": sorted(ignored_intermediate_indices),
                "scope": "stream_stop_policy_only",
                "scientific_input_changed": False,
                "rule": (
                    "Only warnings already shown by optimized identity evidence to be "
                    "declared endpoint returns may be ignored. Any new image warning stops."
                ),
            }
        return result

    monotonic_streak = 0
    stopped_for_monotonic = False
    endpoint_highest_streak = 0
    endpoint_low_force_count = 0
    stopped_for_endpoint_highest = False
    stopped_for_submerged_ci = False
    nimages_match = re.search(r"(?im)^\s*NImages\s+(\d+)\s*$", input_text)
    endpoint_image_indices = (
        {0, int(nimages_match.group(1)) + 1} if nimages_match else set()
    )
    climbing_switch_threshold: float | None = None

    def _stream(line: str) -> bool | None:
        # Preserve ORCA's native output verbatim while allowing a NEB warning
        # to interrupt the expensive single-band calculation immediately.
        nonlocal monotonic_streak, stopped_for_monotonic
        nonlocal endpoint_highest_streak, endpoint_low_force_count
        nonlocal stopped_for_endpoint_highest, climbing_switch_threshold
        nonlocal stopped_for_submerged_ci
        print(line, end="", flush=True)
        threshold_match = re.search(
            r"Switch-on\s+CI\s+threshold\s+([-+0-9.Ee]+)", line, re.IGNORECASE,
        )
        if threshold_match:
            climbing_switch_threshold = float(threshold_match.group(1))
        monotonic_match = _NEB_MONOTONIC_WARNING.search(line) if is_neb else None
        if monotonic_match:
            monotonic_streak += 1
        elif is_neb and _NEB_LBFGS_ITERATION.search(line):
            monotonic_streak = 0
        stopped_for_monotonic = monotonic_streak >= _NEB_MONOTONIC_STOP_STREAK
        row_match = re.search(
            r"LBFGS\s+(\d+)\s+(\d+)\s+([-+0-9.Ee]+)\s+"
            r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
            line, re.IGNORECASE,
        ) if is_neb else None
        if row_match and int(row_match.group(2)) in endpoint_image_indices:
            endpoint_highest_streak += 1
            if (
                climbing_switch_threshold is not None
                and float(row_match.group(4)) <= climbing_switch_threshold
            ):
                endpoint_low_force_count += 1
        elif row_match:
            endpoint_highest_streak = 0
            endpoint_low_force_count = 0
        stopped_for_endpoint_highest = bool(
            endpoint_highest_streak >= _NEB_ENDPOINT_MAX_STOP_STREAK
            and endpoint_low_force_count >= _NEB_ENDPOINT_MAX_LOW_FORCE_COUNT
        )
        if row_match:
            stopped_for_submerged_ci = _neb_submerged_climbing_image_diagnostic(
                output,
            )["valid"]
        warning_indices = (
            _neb_intermediate_indices_from_text(line) if is_neb else []
        )
        unvalidated_intermediate_warning = bool(
            warning_indices
            and not set(warning_indices).issubset(ignored_intermediate_indices)
        )
        return bool(is_neb and (
            unvalidated_intermediate_warning
            or _NEB_ANGLE_WARNING.search(line)
            or stopped_for_monotonic
            or stopped_for_endpoint_highest
            or stopped_for_submerged_ci
        ))

    timeout_seconds = _surface_budget_gate(
        input_path,
        requested_calculations=1,
        ncores=_orca_input_ncores(input_path),
        timeout_seconds=timeout_seconds,
    )
    _write_json(contract_path, {
        "schema_version": 1,
        "kind": "orca_job_input_contract",
        "sha256": digest.hexdigest(),
        "input": str(input_path),
        "referenced_files": referenced,
        "status": "started",
        **native_settings_identity,
    })
    execution_policy_path = input_path.with_suffix(".execution-policy.json")
    execution_policy = (
        json.loads(execution_policy_path.read_text())
        if execution_policy_path.is_file() else {}
    )
    try:
        def native_dispatch(effective_timeout, controlled_environment):
            environment = dict(execution_policy.get("extra_environment") or {})
            if execution_control is not None and any(
                key not in controlled_environment or controlled_environment[key] != value
                for key, value in environment.items()
            ):
                raise RuntimeError("Native per-job environment is outside the frozen runtime contract")
            for key, value in controlled_environment.items():
                if key in environment and environment[key] != value:
                    raise RuntimeError("Native execution policy conflicts with the frozen runtime environment")
                environment[key] = value
            return run_orca(input_path, timeout_seconds=effective_timeout, stream_output=_stream,
                            live_output=False, extra_env=environment)
        if execution_control is None:
            artifacts = native_dispatch(timeout_seconds, {})
        else:
            artifacts = execution_control.dispatch(
                input_path, digest.hexdigest(), timeout_seconds=timeout_seconds,
                ncores=_orca_input_ncores(input_path), operation=native_dispatch,
            )
    except RuntimeError as error:
        error.storca_failed_execution = {
            "input": str(input_path), "output": str(output),
            "job_contract_sha256": hashlib.sha256(contract_path.read_bytes()).hexdigest(),
        }
        if _orca_startup_infrastructure_failure(output):
            retry_input, retry_evidence = _prepare_orca_execution_retry(
                input_path,
                output=output,
                existing_retry_inputs=execution_retry_inputs,
            )
            result = _run_or_resume(
                retry_input, timeout_seconds=timeout_seconds,
                ignore_intermediate_warning_indices=ignored_intermediate_indices,
            )
            result["execution_startup_retry"] = retry_evidence
            result.setdefault("execution_retry_chain_root_input", str(input_path))
            return result
        no_barrier_outcome = _orca_neb_no_barrier_outcome(output)
        existing_intermediates = _neb_intermediate_image_indices(output)
        existing_angle_warning = _neb_angle_warning_image_indices(output)
        existing_monotonic_warning = _neb_monotonic_endpoint_diagnostic(output)
        existing_endpoint_highest = _neb_persistent_endpoint_highest_image_diagnostic(
            output, input_path,
        )
        existing_submerged_ci = _neb_submerged_climbing_image_diagnostic(output)
        if (not no_barrier_outcome["valid"] and not existing_intermediates
                and not existing_angle_warning and not existing_monotonic_warning["valid"]
                and not existing_endpoint_highest["valid"]
                and not existing_submerged_ci["valid"]):
            raise
        artifacts = {
            "out": output,
            "stopped_early": bool(
                existing_intermediates or existing_angle_warning
                or existing_monotonic_warning["valid"]
                or existing_endpoint_highest["valid"]
                or existing_submerged_ci["valid"]
            ),
        }
    except BaseException as error:
        # Recursion executes an append-only basename. Preserve that actual
        # failing segment when the outer route records its incomplete path.
        error.storca_failed_execution = {
            "input": str(input_path), "output": str(output),
            "job_contract_sha256": hashlib.sha256(contract_path.read_bytes()).hexdigest(),
        }
        raise
    normal_termination = _normal_orca_output(Path(artifacts["out"]))
    no_barrier_outcome = _orca_neb_no_barrier_outcome(Path(artifacts["out"]))
    intermediate_indices = _neb_intermediate_image_indices(Path(artifacts["out"]))
    angle_warning_indices = _neb_angle_warning_image_indices(Path(artifacts["out"]))
    monotonic_warning = _neb_monotonic_endpoint_diagnostic(Path(artifacts["out"]))
    endpoint_highest_warning = _neb_persistent_endpoint_highest_image_diagnostic(
        Path(artifacts["out"]), input_path,
    )
    submerged_ci_warning = _neb_submerged_climbing_image_diagnostic(
        Path(artifacts["out"]),
    )
    interrupted_for_intermediate = bool(artifacts.get("stopped_early")) and bool(intermediate_indices)
    interrupted_for_angle = bool(artifacts.get("stopped_early")) and bool(angle_warning_indices)
    interrupted_for_monotonic = bool(artifacts.get("stopped_early")) and monotonic_warning["valid"]
    interrupted_for_endpoint_highest = bool(
        artifacts.get("stopped_early")
    ) and endpoint_highest_warning["valid"]
    interrupted_for_submerged_ci = bool(
        artifacts.get("stopped_early")
    ) and submerged_ci_warning["valid"]
    if (not normal_termination and not no_barrier_outcome["valid"] and not interrupted_for_intermediate
            and not interrupted_for_angle and not interrupted_for_monotonic
            and not interrupted_for_endpoint_highest
            and not interrupted_for_submerged_ci):
        raise RuntimeError(f"ORCA returned without its normal-termination marker: {artifacts['out']}")
    _write_json(contract_path, {
        "schema_version": 1,
        "kind": "orca_job_input_contract",
        "sha256": digest.hexdigest(),
        "input": str(input_path),
        "referenced_files": referenced,
        "status": "completed",
        **native_settings_identity,
    })
    status = (
        "completed_intermediate_detected" if interrupted_for_intermediate else
        "completed_path_quality_warning" if (
            interrupted_for_angle or interrupted_for_monotonic
            or interrupted_for_endpoint_highest
            or interrupted_for_submerged_ci
        ) else
        "completed" if normal_termination else "completed_no_interior_barrier"
    )
    result = {
        "status": status, "resumed": False,
        "input": str(input_path), "output": str(artifacts["out"]),
    }
    if no_barrier_outcome["valid"]:
        result["orca_neb_outcome"] = no_barrier_outcome
    if intermediate_indices:
        result["intermediate_image_indices"] = intermediate_indices
    if angle_warning_indices:
        result["path_quality_warning"] = {
            "kind": "maximum_angle_between_images",
            "image_indices": angle_warning_indices,
        }
    elif monotonic_warning["valid"]:
        result["path_quality_warning"] = monotonic_warning
    elif endpoint_highest_warning["valid"]:
        result["path_quality_warning"] = endpoint_highest_warning
    elif submerged_ci_warning["valid"]:
        result["path_quality_warning"] = submerged_ci_warning
    if interrupted_for_intermediate:
        result["early_stop"] = {
            "status": "stopped_after_orca_intermediate_warning",
            "reason": "ORCA reported an intermediate minimum during the live NEB stream",
            "image_indices": intermediate_indices,
        }
    elif interrupted_for_angle:
        result["early_stop"] = {
            "status": "stopped_after_neb_path_quality_warning",
            "reason": "ORCA reported a maximum-angle path cusp during the live NEB stream",
            "image_indices": angle_warning_indices,
        }
    elif interrupted_for_monotonic:
        result["early_stop"] = {
            "status": "stopped_after_persistent_monotonic_endpoint_maximum",
            "reason": (
                "ORCA retained an endpoint as the highest-energy image for "
                f"{_NEB_MONOTONIC_STOP_STREAK} consecutive NEB iterations"
            ),
            "highest_energy_image": monotonic_warning.get("highest_energy_image"),
        }
    elif interrupted_for_endpoint_highest:
        result["early_stop"] = {
            "status": "stopped_after_persistent_endpoint_highest_image",
            "reason": (
                "No interior highest-energy image emerged during the bounded "
                "endpoint streak despite repeated low-force iterations"
            ),
            "highest_energy_image": endpoint_highest_warning.get("highest_energy_image"),
            "trailing_iteration_count": endpoint_highest_warning.get(
                "trailing_iteration_count"
            ),
            "low_force_iteration_count": endpoint_highest_warning.get(
                "low_force_iteration_count"
            ),
        }
    elif interrupted_for_submerged_ci:
        result["early_stop"] = {
            "status": "stopped_after_persistent_submerged_climbing_image",
            "reason": (
                "The climbing image remained below the higher endpoint asymptote "
                "for the bounded streak"
            ),
            "climbing_image_index": submerged_ci_warning.get("climbing_image_index"),
            "trailing_iteration_count": submerged_ci_warning.get(
                "trailing_iteration_count"
            ),
        }
    return result


def _combine_xyz_files(files: list[Path], output: Path) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(path.read_text().rstrip() + "\n" for path in files))
    return output


def _scan_frames(workdir: Path, label: str) -> list[Path]:
    pattern = re.compile(rf"^{re.escape(label)}\.(\d+)\.xyz$")
    matches = []
    for path in Path(workdir).glob(f"{label}.*.xyz"):
        match = pattern.match(path.name)
        if match:
            matches.append((int(match.group(1)), path))
    return [path for _, path in sorted(matches)]


def _scan_energies(workdir: Path, label: str, expected: int, scan_files: list[Path] | None = None) -> list[float] | None:
    candidates = [Path(workdir) / f"{label}.relaxscanact.dat", Path(workdir) / f"{label}_relaxscanact.dat"]
    path = next((candidate for candidate in candidates if candidate.is_file()), None)
    energies = []
    if path is not None:
        for line in path.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            values = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", line)
            if len(values) >= 2:
                energies.append(float(values[-1]))
    elif scan_files:
        for scan_file in scan_files:
            lines = scan_file.read_text().splitlines()
            match = re.search(r"\bE\s+(-?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)", lines[1] if len(lines) > 1 else "")
            if match:
                energies.append(float(match.group(1)))
    return energies if len(energies) == expected else None


def _validate_dissociation_path(trajectory: Path, bond: tuple[int, int],
                                final_distance: float) -> dict:
    frames = xyz_frames(trajectory)
    elements = frames[0][0]
    if any(frame[0] != elements for frame in frames):
        return {"status": "atom_order_changed", "valid": False}
    left, right = bond
    distances = [
        float(((coordinates[left] - coordinates[right]) ** 2).sum() ** 0.5)
        for _, coordinates, _ in frames
    ]
    decreases = sum(current + 0.05 < previous for previous, current in zip(distances, distances[1:]))
    valid = len(frames) >= 3 and distances[-1] >= final_distance - 0.2 and decreases <= max(1, len(distances) // 10)
    return {
        "status": "validated" if valid else "dissociation_coordinate_mismatch",
        "valid": valid,
        "frame_count": len(frames),
        "bond_atom_indices": list(bond),
        "bond_distances_angstrom": distances,
    }


def _locate_neb_ts(workdir: Path, label: str) -> Path | None:
    # A restart shares a directory and prefix with its predecessors. Prefix
    # globs can borrow a different execution's saddle, so use exact native
    # artifact names only; missing evidence is not a directory-wide search.
    candidates = [
        Path(workdir) / f"{label}{suffix}"
        for suffix in ("_NEB-TS_converged.xyz", "_NEB-TS.xyz",
                       "_TS_converged.xyz", "_TS.xyz", ".xyz")
    ]
    return next((path for path in candidates if path.is_file()), None)


def _locate_irc_trajectory(workdir: Path, label: str) -> Path | None:
    trajectory = Path(workdir) / f"{label}_IRC_Full_trj.xyz"
    return trajectory if trajectory.is_file() else None


def _frequency_temperature_kwargs(anchor: Path) -> dict[str, float]:
    """Inherit only an explicit, hash-validated parent frequency contract."""
    context = _surface_context(anchor)
    value = (context or {}).get("electronic_structure", {}).get("thermochemistry_temperature_K")
    return ({"thermochemistry_temperature_K": normalize_thermochemistry_temperature(value)}
            if value is not None else {})


def validated_orca_temperature_thermochemistry(output: Path, temperature_K: float) -> dict:
    """Bind scalar totals to the exact frequency request and printed T evidence."""
    temperature = normalize_thermochemistry_temperature(temperature_K)
    output = Path(output)
    input_evidence = orca_frequency_temperature_input_evidence(output.with_suffix(".inp"), temperature)
    record = parse_orca_thermochemistry(output)
    if not record.get("valid"):
        raise ValueError("Requested-temperature ORCA thermochemistry is incomplete")
    temperatures = re.findall(
        r"(?mi)^\s*Temperature\s+\.\.\.\s+([-+0-9.Ee]+)\s+K", output.read_text(errors="replace"),
    )
    if len(temperatures) != 1:
        raise ValueError("Scalar thermochemistry requires one output temperature section")
    # ORCA prints this line to 0.01 K. Preserve printed versus exact requested
    # values separately; this tolerance is output rounding, not uncertainty.
    if not math.isclose(float(record["temperature_K"]), temperature, rel_tol=0.0, abs_tol=0.00500001):
        raise ValueError("ORCA thermochemistry temperature does not match the requested temperature")
    return {
        **record, "requested_temperature_K": temperature,
        "temperature_input_evidence": input_evidence,
        "temperature_output_rounding_tolerance_K": 0.00500001,
        "source_output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }


def _retained_endpoint_workflow(
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    ncores: int,
    method_keywords: list[str] | None,
) -> dict | None:
    """Adopt an exact retained pre-contract workflow without overwriting it."""
    workdir = Path(workdir)
    metadata_path = workdir / "metadata.json"
    optimized = workdir / "opt.xyz"
    optimization_output = workdir / "opt.out"
    if not (metadata_path.is_file() and optimized.is_file() and _normal_orca_output(optimization_output)):
        return None
    try:
        metadata = json.loads(metadata_path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    expected_method = list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"])
    settings = active_native_settings()
    temperature_contract = _frequency_temperature_kwargs(workdir)
    if (
        int(metadata.get("charge", 10**9)) != int(charge)
        or int(metadata.get("multiplicity", 0)) != int(multiplicity)
        or int(metadata.get("ncores", 0)) != int(ncores)
        or metadata.get("orca_method_keywords") != expected_method
        or metadata.get("run_frequency") is not True
        or (settings is not None and metadata.get("native_job_settings") != settings.as_dict())
        or (temperature_contract and metadata.get("thermochemistry_temperature_K")
            != temperature_contract["thermochemistry_temperature_K"])
    ):
        return None
    result: dict[str, Any] = {
        "optimized_xyz": optimized,
        "optimization": {"out": optimization_output},
        "resumed": True,
        "legacy_evidence_adopted": True,
    }
    frequency = workdir / "freq.out"
    if _normal_orca_output(frequency):
        if temperature_contract:
            validated_orca_temperature_thermochemistry(frequency, temperature_contract["thermochemistry_temperature_K"])
        result.update(
            frequency={"out": frequency},
            frequency_check=frequency_stability_check(frequency),
        )
    return result


def _run_stationary_point_recovery(
    initial_result: dict,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    policy: SurfaceResolutionPolicy,
) -> dict:
    """Continue a normally terminated, unconverged optimization append-only."""
    workdir = Path(workdir)
    recovery_root = workdir / "stationary-point-recovery"
    recovery_root.mkdir(parents=True, exist_ok=True)
    manifest = recovery_root / "stationary-point-recovery.json"
    source_xyz = Path(initial_result["optimized_xyz"])
    source_elements, _, _ = read_xyz(source_xyz)
    strategies = (
        {
            "coordinate_system": "redundant", "max_iterations": 300,
            "calculate_hessian": True, "recalculate_hessian_every": None,
            "trust_radius": 0.10,
        },
        {
            "coordinate_system": "cartesian", "max_iterations": 300,
            "calculate_hessian": True, "recalculate_hessian_every": None,
            "trust_radius": 0.10,
        },
        {
            "coordinate_system": "redundant", "max_iterations": 400,
            "calculate_hessian": True, "recalculate_hessian_every": 5,
            "trust_radius": 0.05,
        },
    )[:policy.maximum_stationary_point_recoveries]
    attempts: list[dict[str, Any]] = []
    for position, strategy in enumerate(strategies, start=1):
        attempt_dir = recovery_root / f"attempt-{position:02d}"
        attempt_dir.mkdir(parents=True, exist_ok=True)
        seed = attempt_dir / "input.xyz"
        if not seed.is_file():
            shutil.copy2(source_xyz, seed)
        else:
            retained_elements, _, _ = read_xyz(seed)
            if retained_elements != source_elements:
                raise RuntimeError("Retained stationary-point recovery seed changed atom order")
        contract = {
            "schema_version": 1,
            "kind": "orca_stationary_point_recovery_contract",
            "attempt": position,
            "source_xyz": str(source_xyz),
            "source_xyz_sha256": hashlib.sha256(source_xyz.read_bytes()).hexdigest(),
            "seed_xyz_sha256": hashlib.sha256(seed.read_bytes()).hexdigest(),
            "charge": int(charge),
            "multiplicity": int(multiplicity),
            "ncores": int(ncores),
            "method_keywords": list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]),
            "strategy": strategy,
            "surface_resolution_policy_id": policy.policy_id,
        }
        contract_path = attempt_dir / "attempt-contract.json"
        if contract_path.is_file() and json.loads(contract_path.read_text()) != contract:
            raise RuntimeError("Retained stationary-point recovery contract was mutated")
        _write_json(contract_path, contract)
        attempt: dict[str, Any] = {
            "attempt": position,
            "status": "started",
            "contract": str(contract_path),
            "contract_sha256": hashlib.sha256(contract_path.read_bytes()).hexdigest(),
            "strategy": strategy,
            "claims_permitted": [],
            "claims_prohibited": ["intermediate_minimum", "rate", "lifetime"],
        }

        def retain_attempt_evidence() -> None:
            inputs = [seed, contract_path]
            for key in ("optimization", "frequency_execution"):
                input_value = (attempt.get(key) or {}).get("input")
                if input_value:
                    inputs.append(Path(str(input_value)))
            attempt["surface_resolution_attempt_evidence"] = _write_surface_attempt_evidence(
                attempt_dir,
                strategy="stationary_point_recovery",
                status=str(attempt.get("status") or "recovery_execution_failed"),
                result=attempt,
                parameters={
                    "attempt": position,
                    "optimizer_strategy": strategy,
                    "charge": int(charge),
                    "multiplicity": int(multiplicity),
                    "ncores": int(ncores),
                    "method_keywords": list(
                        method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]
                    ),
                    "surface_resolution_policy": policy.as_dict(),
                },
                input_paths=inputs,
                failure_reason=attempt.get("failure_reason"),
                claims_permitted=(
                    ("validated_stationary_minimum",)
                    if attempt.get("status") == "validated_minimum" else ()
                ),
            )
        try:
            opt_input = create_orca_optimization_recovery_input(
                seed,
                charge=charge,
                multiplicity=multiplicity,
                label="opt",
                ncores=ncores,
                method_keywords=method_keywords,
                **strategy,
            )
            execution = _run_or_resume(opt_input, timeout_seconds=timeout_seconds)
            executed_input = Path(execution.get("input") or opt_input)
            opt_output = Path(execution.get("output") or executed_input.with_suffix(".out"))
            optimized = executed_input.with_suffix(".xyz")
            attempt["optimization"] = execution
            if not optimized.is_file():
                raise RuntimeError("ORCA recovery produced no retained geometry")
            optimized_elements, _, _ = read_xyz(optimized)
            if optimized_elements != source_elements:
                raise RuntimeError("ORCA recovery changed composition or atom order")
            # Continue later attempts from the most recent valid geometry even
            # when this attempt reaches its iteration bound.
            source_xyz = optimized
            if not _orca_optimization_converged(opt_output):
                attempt.update(
                    status="intermediate_optimization_unconverged",
                    last_geometry=str(optimized),
                )
                retain_attempt_evidence()
                attempts.append(attempt)
                _write_json(manifest, {
                    "schema_version": 1,
                    "kind": "orca_stationary_point_recovery",
                    "status": "running",
                    "attempts": attempts,
                })
                continue
            freq_input = create_orca_input(
                optimized,
                charge=charge,
                multiplicity=multiplicity,
                freq=True,
                label="freq",
                ncores=ncores,
                method_keywords=method_keywords,
                preserve_existing=bool(_frequency_temperature_kwargs(workdir)),
                **_frequency_temperature_kwargs(workdir),
            )
            frequency_execution = _run_or_resume(freq_input, timeout_seconds=timeout_seconds)
            frequency_output = Path(
                frequency_execution.get("output") or freq_input.with_suffix(".out")
            )
            temperature_contract = _frequency_temperature_kwargs(workdir)
            if temperature_contract:
                validated_orca_temperature_thermochemistry(
                    frequency_output, temperature_contract["thermochemistry_temperature_K"],
                )
            check = frequency_stability_check(
                frequency_output,
                imag_tol=abs(policy.meaningful_imaginary_frequency_cm_1),
            )
            attempt.update(
                status="validated_minimum" if check.get("IsMinimum") else "intermediate_nonminimum",
                optimized_xyz=str(optimized),
                frequency_execution=frequency_execution,
                frequency_output=str(frequency_output),
                frequency_check=check,
            )
            retain_attempt_evidence()
            attempts.append(attempt)
            result = {
                "optimized_xyz": optimized,
                "optimization": {"out": opt_output, **execution},
                "frequency": {"out": frequency_output, **frequency_execution},
                "frequency_check": check,
                "stationary_point_recovery": {
                    "schema_version": 1,
                    "kind": "orca_stationary_point_recovery",
                    "status": attempt["status"],
                    "attempts": attempts,
                },
                "resumed": bool(execution.get("resumed")),
            }
            _write_json(manifest, result["stationary_point_recovery"])
            if check.get("IsMinimum"):
                return result
            recovered = _recover_endpoint_minimum(
                result,
                workdir,
                charge=charge,
                multiplicity=multiplicity,
                ncores=ncores,
                method_keywords=method_keywords,
                timeout_seconds=timeout_seconds,
                imaginary_frequency_threshold_cm_1=policy.meaningful_imaginary_frequency_cm_1,
            )
            recovered["stationary_point_recovery"] = result["stationary_point_recovery"]
            return recovered
        except Exception as error:
            if isinstance(error, (SurfaceResolutionBudgetExhausted, TransferControlError)):
                raise
            attempt.update(
                status="recovery_execution_failed",
                failure_reason=f"{type(error).__name__}: {error}",
            )
            retain_attempt_evidence()
            attempts.append(attempt)
            _write_json(manifest, {
                "schema_version": 1,
                "kind": "orca_stationary_point_recovery",
                "status": "running",
                "attempts": attempts,
            })
    result = {
        **initial_result,
        "initial_nonstationary_frequency_check": initial_result.get("frequency_check"),
        "frequency_check": {
            "status": "invalid_nonstationary_frequency",
            "IsMinimum": False,
            "reason": "No stationary-point recovery optimization reached its convergence marker.",
        },
        "stationary_point_recovery": {
            "schema_version": 1,
            "kind": "orca_stationary_point_recovery",
            "status": "intermediate_optimization_unconverged",
            "attempts": attempts,
            "claims_permitted": [],
            "claims_prohibited": ["intermediate_minimum", "rate", "lifetime"],
        },
    }
    _write_json(manifest, result["stationary_point_recovery"])
    return result


def _orca_irc_completion(
    output: Path, *, required_directions: tuple[str, ...] = ("forward", "backward"),
) -> dict:
    """Require both IRC direction tables to meet their printed tolerances."""
    output = Path(output)
    text = output.read_text(errors="replace") if output.is_file() else ""
    if required_directions not in {("forward", "backward"), ("downhill",)}:
        raise ValueError("Unsupported IRC completion direction contract")
    sections = list(re.finditer(r"\b(FORWARD|BACKWARD|DOWNHILL)\s+IRC\b", text))
    directions = {}
    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
    for index, heading in enumerate(sections):
        section = text[heading.end():sections[index + 1].start() if index + 1 < len(sections) else len(text)]
        section = section.split("IRC PATH SUMMARY", 1)[0]
        thresholds = re.search(rf"Convergence thresholds\s+({number})\s+({number})", section)
        rows = list(re.finditer(
            rf"(?m)^\s*(\d+)\s+({number})\s+({number})\s+({number})\s+({number})\s*$", section,
        ))
        limited = "MAXIMUM NUMBER OF ITERATIONS REACHED" in section
        record = {"iteration_limit_reached": limited, "valid": False}
        if thresholds and rows:
            values = [float(value.replace("D", "E").replace("d", "e")) for value in (
                thresholds.group(1), thresholds.group(2), rows[-1].group(4), rows[-1].group(5),
            )]
            max_tol, rms_tol, max_gradient, rms_gradient = values
            record.update(
                last_iteration=int(rows[-1].group(1)),
                maximum_gradient_Eh_per_bohr=max_gradient,
                rms_gradient_Eh_per_bohr=rms_gradient,
                maximum_gradient_tolerance=max_tol,
                rms_gradient_tolerance=rms_tol,
                valid=bool(
                    all(np.isfinite(values)) and min(max_tol, rms_tol) > 0
                    and 0 <= max_gradient <= max_tol and 0 <= rms_gradient <= rms_tol
                    and not limited
                ),
            )
        directions[heading.group(1).lower()] = record
    normal = _normal_orca_output(output)
    valid = bool(normal and set(directions) == set(required_directions)
                 and all(item["valid"] for item in directions.values()))
    limited = "MAXIMUM NUMBER OF ITERATIONS REACHED" in text
    return {
        "schema_version": 1,
        "kind": ("orca_bidirectional_irc_completion" if required_directions == ("forward", "backward")
                 else "orca_directional_irc_completion"),
        "required_directions": list(required_directions),
        "status": ("irc_iteration_limit_reached" if limited else
                   ("irc_downhill_converged" if required_directions == ("downhill",)
                    else "irc_bidirectional_converged") if valid else "irc_completion_unresolved"),
        "valid": valid and not limited,
        "normal_termination": normal,
        "directions": directions,
        "source_output": str(output),
        "source_output_sha256": hashlib.sha256(output.read_bytes()).hexdigest() if output.is_file() else None,
        "interpretation": "Program termination and endpoint connectivity alone do not establish completed bidirectional IRC tracing.",
    }


def _downhill_irc_input_text(source_text: str, seed_name: str) -> str:
    blocks = list(re.finditer(r"(?ims)^\s*%irc\b.*?^\s*end\s*$", source_text))
    if len(blocks) != 1:
        raise ValueError("IRC continuation requires one IRC block")
    block = re.sub(r"(?im)^(\s*Direction\s+)both\s*$", r"\1down", blocks[0].group())
    block = re.sub(r"(?im)^\s*(?:InitHess|Hess_Filename)\b[^\n]*\n?", "", block)
    text = source_text[:blocks[0].start()] + block + source_text[blocks[0].end():]
    return re.sub(
        r"(?im)^(\*\s+xyzfile\s+[-+]?\d+\s+\d+\s+)\S+\s*$",
        rf"\g<1>{seed_name}\n", text,
    )


def _validate_irc_continuation_blocks(source_input: Path, source_text: str) -> None:
    """Relocate only legacy blocks or the exact source-bound scoped profile."""
    names = [value.casefold() for value in re.findall(r"(?im)^\s*%(\w+)", source_text)]
    settings = active_native_settings()
    allowed = {"pal", "irc"}
    if settings is not None:
        contract = json.loads(Path(source_input).with_suffix(".contract.json").read_text())
        if contract.get("native_settings_sha256") != settings.sha256:
            raise ValueError("IRC continuation native settings identity differs from its source contract")
        if names.count("scf") != 1 or names.count("maxcore") != 1:
            raise ValueError("IRC continuation requires the complete scoped SCF and MaxCore controls")
        if apply_native_job_settings(source_text, settings) != source_text:
            raise ValueError("IRC continuation source settings do not match the frozen profile")
        allowed.update({"scf", "maxcore"})
    if set(names) - allowed:
        raise ValueError("IRC continuation cannot relocate unsupported source input blocks")


def _irc_downhill_execution_sources(root_input: Path, execution: dict) -> dict[str, str]:
    """Bind a downhill continuation to its actual append-only execution chain."""
    root_input = Path(root_input).resolve()
    selected = Path(execution.get("input") or "").resolve()
    output = Path(execution.get("output") or "").resolve()
    pattern = re.compile(re.escape(root_input.stem) + r"-execution-retry-\d{3}$")
    if (execution.get("status") != "completed" or selected.parent != root_input.parent
            or selected.suffix != ".inp" or output != selected.with_suffix(".out")
            or (selected != root_input and not pattern.fullmatch(selected.stem))):
        raise ValueError("IRC continuation execution ownership is inconsistent")
    sources = {}

    def retain(path):
        path = Path(path).resolve()
        if not path.is_file():
            raise ValueError(f"IRC continuation execution evidence is missing: {path}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if str(path) in sources and sources[str(path)] != digest:
            raise ValueError("IRC continuation execution has conflicting retained source hashes")
        sources[str(path)] = digest

    visited = set()
    current = selected
    while True:
        if current in visited:
            raise ValueError("IRC continuation execution ownership contains a retry cycle")
        visited.add(current)
        retain(current)
        job_contract_path = current.with_suffix(".contract.json")
        retain(job_contract_path)
        job_contract = json.loads(job_contract_path.read_text())
        digest, referenced = _orca_job_input_evidence(current)
        if (job_contract.get("sha256") != digest.hexdigest()
                or Path(job_contract.get("input") or "").resolve() != current):
            raise ValueError("IRC continuation execution input contract hash changed")
        for reference in referenced:
            retain(reference)
        if _orca_serial_execution_input(current.read_text()) != _orca_serial_execution_input(root_input.read_text()):
            raise ValueError("IRC continuation execution method or geometry changed")
        if current == root_input:
            break
        restart_path = current.with_suffix(".restart.json")
        retain(restart_path)
        restart = json.loads(restart_path.read_text())
        prior = Path(restart.get("source_input") or "").resolve()
        if (restart.get("kind") not in {
                "orca_startup_infrastructure_execution_retry", "orca_interrupted_execution_retry",
                "orca_execution_resource_contract_recovery",
            } or restart.get("scientific_input_changed") is not False
                or prior.parent != root_input.parent or prior.suffix != ".inp"
                or (prior != root_input and not pattern.fullmatch(prior.stem))
                or Path(restart.get("retry_input") or restart.get("recovery_input") or "").resolve() != current):
            raise ValueError("IRC continuation execution retry ownership is inconsistent")
        if restart.get("kind") != "orca_execution_resource_contract_recovery":
            prior_output = prior.with_suffix(".out")
            if Path(restart.get("source_output") or "").resolve() != prior_output:
                raise ValueError("IRC continuation execution retry output ownership changed")
            retain(prior_output)
            if restart.get("source_output_sha256") != sources[str(prior_output)]:
                raise ValueError("IRC continuation execution retry source output hash changed")
        policy_value = restart.get("execution_policy")
        if policy_value:
            policy_path = Path(policy_value).resolve()
            if policy_path != current.with_suffix(".execution-policy.json"):
                raise ValueError("IRC continuation execution policy ownership changed")
            retain(policy_path)
            policy = json.loads(policy_path.read_text())
            if policy.get("scientific_model_changed") is not False:
                raise ValueError("IRC continuation execution policy changed the scientific model")
            prior_field = ("prior_ncores" if restart.get("kind") == "orca_execution_resource_contract_recovery"
                           else "requested_ncores")
            for field, expected_ncores in (("effective_ncores", _orca_input_ncores(current)),
                                           (prior_field, _orca_input_ncores(prior))):
                declared = policy.get(field)
                if type(declared) is not int or declared < 1 or declared != expected_ncores:
                    raise ValueError("IRC continuation execution policy core counts differ from retained inputs")
            if restart.get("kind") == "orca_execution_resource_contract_recovery":
                if (Path(policy.get("source_input") or "").resolve() != prior
                        or policy.get("source_input_sha256") != hashlib.sha256(prior.read_bytes()).hexdigest()):
                    raise ValueError("IRC continuation execution resource-recovery source hash changed")
        elif restart.get("kind") == "orca_execution_resource_contract_recovery":
            raise ValueError("IRC continuation execution resource recovery lacks its policy")
        current = prior
    retain(output)
    return sources


def run_irc_branch_continuation(
    source_input: Path, *, branch: str, timeout_seconds: float | None,
) -> dict:
    """Extend one capped IRC branch downhill once, without admitting a path.

    ORCA 6.1 documents Direction down for a previous unconverged IRC endpoint.
    This bounded diagnostic retains the source method, charge, spin, resources,
    and MaxIter. It is not an endpoint-identity or kinetic admission mechanism.
    """
    if branch not in {"forward", "backward"}:
        raise ValueError("IRC continuation branch must be forward or backward")
    source_input = Path(source_input).resolve()
    source_output = source_input.with_suffix(".out")
    completion = _orca_irc_completion(source_output)
    direction = completion["directions"].get(branch) or {}
    if not completion["normal_termination"] or not direction.get("iteration_limit_reached"):
        raise ValueError("Only a normally terminated, iteration-limited IRC branch can be continued")
    if not source_input.with_suffix(".contract.json").is_file():
        raise ValueError("IRC continuation requires the original immutable job contract")
    source_execution = _run_or_resume(source_input, timeout_seconds=timeout_seconds)
    if not source_execution.get("resumed"):
        raise RuntimeError("IRC continuation unexpectedly reexecuted its source")
    source_text = source_input.read_text()
    _validate_irc_continuation_blocks(source_input, source_text)
    blocks = list(re.finditer(r"(?ims)^\s*%irc\b.*?^\s*end\s*$", source_text))
    if len(blocks) != 1 or not re.search(r"(?im)^\s*Direction\s+both\s*$", blocks[0].group()):
        raise ValueError("Continuation requires a single original bidirectional IRC block")
    block = blocks[0].group()
    bound = re.search(r"(?im)^\s*MaxIter\s+(\d+)\s*$", block)
    if bound is None or not 1 <= int(bound.group(1)) <= 200:
        raise ValueError("IRC continuation requires an original MaxIter bound between 1 and 200")
    suffix = "F" if branch == "forward" else "B"
    endpoint = source_input.parent / f"{source_input.stem}_IRC_{suffix}.xyz"
    trajectory = source_input.parent / f"{source_input.stem}_IRC_{suffix}_trj.xyz"
    endpoint_elements, endpoint_coordinates, _ = read_xyz(endpoint)
    frames = xyz_frames(trajectory)
    if endpoint_elements != frames[-1][0] or not np.allclose(
        endpoint_coordinates, frames[-1][1], atol=1e-6, rtol=0,
    ):
        raise ValueError("IRC branch endpoint does not match its retained trajectory tail")
    geometry = re.search(r"(?im)^\*\s+xyzfile\s+([-+]?\d+)\s+(\d+)\s+(\S+)\s*$", source_text)
    if geometry is None:
        raise ValueError("IRC source lacks a uniquely retained XYZ geometry reference")
    original_geometry = source_input.parent / geometry.group(3).strip('"')
    if read_xyz(original_geometry)[0] != endpoint_elements:
        raise ValueError("IRC branch changed atom composition or order")
    workdir = source_input.parent / "irc-branch-continuation" / branch
    contract = {
        "schema_version": 1, "kind": "bounded_irc_downhill_continuation",
        "policy_id": "irc-downhill-single-extension-v1", "branch": branch,
        "maximum_iterations": int(bound.group(1)),
        "source_input": str(source_input),
        "source_input_sha256": hashlib.sha256(source_input.read_bytes()).hexdigest(),
        "source_output_sha256": completion["source_output_sha256"],
        "source_endpoint_sha256": hashlib.sha256(endpoint.read_bytes()).hexdigest(),
        "source_trajectory_sha256": hashlib.sha256(trajectory.read_bytes()).hexdigest(),
        "source_completion": completion,
        **({"native_settings_sha256": active_native_settings().sha256}
           if active_native_settings() is not None else {}),
        "claims_permitted": [], "claims_prohibited": ["verified_path", "product_identity", "rate", "lifetime"],
    }
    contract_path = workdir / "continuation-contract.json"
    if contract_path.is_file() and json.loads(contract_path.read_text()) != contract:
        raise RuntimeError("Retained IRC continuation source evidence changed")
    workdir.mkdir(parents=True, exist_ok=True)
    seed = workdir / "branch-endpoint.xyz"
    if seed.is_file() and seed.read_bytes() != endpoint.read_bytes():
        raise RuntimeError("Retained IRC continuation seed changed")
    if not seed.is_file():
        shutil.copy2(endpoint, seed)
    input_text = _downhill_irc_input_text(source_text, seed.name)
    inp = workdir / "irc-downhill.inp"
    if inp.is_file() and inp.read_text() != input_text:
        raise RuntimeError("Retained downhill IRC input changed")
    if not inp.is_file():
        inp.write_text(input_text)
    if not contract_path.exists():
        _write_json(contract_path, contract)
    retained_output = inp.with_suffix(".out")
    if retained_output.is_file() and not _normal_orca_output(retained_output) \
            and not _orca_startup_infrastructure_failure(retained_output) \
            and not _detached_orca_execution_is_stopped(inp):
        raise RuntimeError("IRC continuation may still be active; refusing duplicate execution")
    execution = _run_or_resume(inp, timeout_seconds=timeout_seconds)
    result = {**contract, "execution": execution,
              "completion": _orca_irc_completion(Path(execution["output"]), required_directions=("downhill",)),
              "status": "continuation_executed_unadmitted"}
    result["execution_source_hashes"] = _irc_downhill_execution_sources(inp, execution)
    executed_input = Path(execution["input"])
    for suffix in ("_IRC_D_trj.xyz", "_IRC_D.xyz"):
        artifact = executed_input.parent / f"{executed_input.stem}{suffix}"
        if artifact.is_file():
            result["execution_source_hashes"][str(artifact.resolve())] = hashlib.sha256(artifact.read_bytes()).hexdigest()
    _write_json(workdir / "continuation-result.json", result)
    return result


def assemble_completed_irc_branches(source_input: Path) -> dict:
    """Join converged, hash-linked downhill segments; do not admit chemistry."""
    source_input = Path(source_input).resolve()
    root = source_input.parent
    source_output = source_input.with_suffix(".out")
    completion = _orca_irc_completion(source_output)
    incomplete = {
        "kind": "storca_irc_branch_assembly", "schema_version": 1,
        "status": "irc_branch_assembly_incomplete", "valid": False,
        "claims_permitted": [], "claims_prohibited": ["verified_path", "product_identity", "rate", "lifetime"],
    }
    if not completion["normal_termination"] or set(completion["directions"]) != {"forward", "backward"}:
        return {**incomplete, "reason": "original_direction_evidence_incomplete"}
    full = root / f"{source_input.stem}_IRC_Full_trj.xyz"
    if not full.is_file():
        return {**incomplete, "reason": "original_full_trajectory_missing"}
    sources = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
               for path in (source_input, source_output, full)}
    base_frames = xyz_frames(full)
    elements = base_frames[0][0]
    tails = {}
    extras = {}
    diagnostics = {}

    def same_frame(left, right):
        return left[0] == right[0] and np.allclose(left[1], right[1], atol=1e-6, rtol=0)

    for branch, suffix in (("forward", "F"), ("backward", "B")):
        endpoint = root / f"{source_input.stem}_IRC_{suffix}.xyz"
        branch_trajectory = root / f"{source_input.stem}_IRC_{suffix}_trj.xyz"
        if not endpoint.is_file() or not branch_trajectory.is_file():
            return {**incomplete, "reason": "original_branch_artifacts_missing", "branch": branch}
        for path in (endpoint, branch_trajectory):
            sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        tail = read_xyz(endpoint)
        if not same_frame(tail, xyz_frames(branch_trajectory)[-1]):
            raise ValueError("Original IRC endpoint and trajectory disagree")
        tails[branch] = tail
        direction = completion["directions"][branch]
        extras[branch] = []
        diagnostics[branch] = direction
        if direction["valid"]:
            continue
        if not direction["iteration_limit_reached"]:
            return {**incomplete, "reason": "branch_not_eligible_for_downhill_extension", "branch": branch}
        directory = root / "irc-branch-continuation" / branch
        contract_path = directory / "continuation-contract.json"
        if not contract_path.is_file():
            return {**incomplete, "reason": "branch_continuation_missing", "branch": branch}
        sources[str(contract_path)] = hashlib.sha256(contract_path.read_bytes()).hexdigest()
        contract = json.loads(contract_path.read_text())
        expected = {
            "policy_id": "irc-downhill-single-extension-v1", "branch": branch,
            "source_input_sha256": sources[str(source_input)],
            "source_output_sha256": sources[str(source_output)],
            "source_endpoint_sha256": sources[str(endpoint)],
            "source_trajectory_sha256": sources[str(branch_trajectory)],
        }
        if any(contract.get(key) != value for key, value in expected.items()):
            raise ValueError("IRC continuation source hashes or branch contract changed")
        inp = directory / "irc-downhill.inp"
        seed = directory / "branch-endpoint.xyz"
        if not inp.is_file() or not seed.is_file():
            return {**incomplete, "reason": "continuation_input_missing", "branch": branch}
        for path in (inp, seed):
            sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        if inp.read_text() != _downhill_irc_input_text(source_input.read_text(), seed.name) \
                or seed.read_bytes() != endpoint.read_bytes():
            raise ValueError("IRC continuation method or seed changed")
        result_path = directory / "continuation-result.json"
        if not result_path.is_file():
            return {**incomplete, "reason": "continuation_execution_evidence_missing", "branch": branch}
        sources[str(result_path)] = hashlib.sha256(result_path.read_bytes()).hexdigest()
        continuation_result = json.loads(result_path.read_text())
        if any(continuation_result.get(key) != value for key, value in contract.items()):
            raise ValueError("IRC continuation execution source contract changed")
        execution = continuation_result.get("execution") or {}
        execution_sources = _irc_downhill_execution_sources(inp, execution)
        recorded = continuation_result.get("execution_source_hashes") or {}
        if (Path(execution["input"]).resolve() != inp
                and (not recorded or any(recorded.get(path) != digest
                                         for path, digest in execution_sources.items()))):
            raise ValueError("IRC continuation execution retry lacks complete source hashes")
        for path, digest in recorded.items():
            retained = Path(path).resolve()
            if (retained.parent != directory or not retained.is_file()
                    or hashlib.sha256(retained.read_bytes()).hexdigest() != digest):
                raise ValueError("IRC continuation execution source hash changed")
        executed_input = Path(execution["input"]).resolve()
        output = Path(execution["output"]).resolve()
        diagnostic = _orca_irc_completion(output, required_directions=("downhill",))
        if diagnostic["source_output_sha256"] != execution_sources[str(output)]:
            raise ValueError("IRC continuation execution output hash changed during admission")
        diagnostics[branch] = diagnostic
        if not diagnostic["valid"]:
            return {**incomplete, "reason": "downhill_branch_not_converged", "branch": branch,
                    "completion": diagnostic}
        continued_trajectory = directory / f"{executed_input.stem}_IRC_D_trj.xyz"
        final_endpoint = directory / f"{executed_input.stem}_IRC_D.xyz"
        if not continued_trajectory.is_file() or not final_endpoint.is_file():
            return {**incomplete, "reason": "downhill_trajectory_missing", "branch": branch}
        for path in (continued_trajectory, final_endpoint):
            sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            if executed_input != inp and recorded.get(str(path)) != sources[str(path)]:
                raise ValueError("IRC continuation execution trajectory source hash changed")
        frames = xyz_frames(continued_trajectory)
        if not same_frame(frames[0], tail) or not same_frame(frames[-1], read_xyz(final_endpoint)):
            raise ValueError("Downhill IRC trajectory has a discontinuous junction or endpoint")
        extras[branch] = frames[1:]
        if any(path in sources and sources[path] != digest
               for path, digest in execution_sources.items()):
            raise ValueError("IRC assembly has conflicting retained source hashes")
        sources.update(execution_sources)

    if not (same_frame(base_frames[0], tails["backward"]) and same_frame(base_frames[-1], tails["forward"])):
        raise ValueError("Original full IRC does not connect its declared directional tails")
    frames = list(reversed(extras["backward"])) + base_frames + extras["forward"]
    if any(frame[0] != elements or not np.isfinite(frame[1]).all() for frame in frames):
        raise ValueError("IRC assembly changed atom order or contains nonfinite coordinates")
    if any(not Path(path).is_file() or hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest
           for path, digest in sources.items()):
        raise ValueError("IRC assembly source evidence changed during admission")
    identity = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    directory = root / "irc-branch-continuation" / f"assembled-{identity[:16]}"
    directory.mkdir(parents=True, exist_ok=True)
    trajectory = directory / "trajectory.xyz"
    text = "".join(
        f"{len(labels)}\n{comment}\n" + "".join(
            f"{label} {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f}\n"
            for label, xyz in zip(labels, coordinates)
        ) for labels, coordinates, comment in frames
    )
    if trajectory.is_file() and trajectory.read_text() != text:
        raise ValueError("Retained assembled IRC trajectory changed")
    if not trajectory.is_file():
        trajectory.write_text(text)
    result = {**incomplete, "status": "irc_branches_joined", "valid": True,
              "claims_permitted": ["completed_continuous_irc_tracing"],
              "source_hashes": sources, "direction_completion": diagnostics,
              "trajectory": str(trajectory), "frame_count": len(frames),
              "trajectory_sha256": hashlib.sha256(trajectory.read_bytes()).hexdigest()}
    manifest = directory / "assembly.json"
    if manifest.is_file() and json.loads(manifest.read_text()) != result:
        raise ValueError("Retained IRC assembly evidence changed")
    if not manifest.is_file():
        _write_json(manifest, result)
    return result


def _run_endpoint_job(input_path: Path, *, timeout_seconds: float | None) -> dict:
    """Use the append-only path executor for endpoint Opt/Freq jobs too."""
    input_path = Path(input_path)
    output = input_path.with_suffix(".out")
    if output.is_file() and not _normal_orca_output(output):
        if not _orca_startup_infrastructure_failure(output) and not \
                _detached_orca_execution_is_stopped(input_path):
            raise RuntimeError(
                f"Retained endpoint output is incomplete; refusing to overwrite it: {output}"
            )
        # Old endpoint workflows did not create job contracts. Preserve their
        # startup failure and record its hash before an execution-only retry.
        if not input_path.with_suffix(".contract.json").exists():
            if not _orca_startup_infrastructure_failure(output):
                raise RuntimeError("Interrupted legacy endpoint lacks an immutable job contract")
            retries = sorted(input_path.parent.glob(
                f"{input_path.stem}-execution-retry-*.inp"
            ))
            if retries:
                original_evidence = json.loads(retries[0].with_suffix(".restart.json").read_text())
                if retries[0].read_bytes() != input_path.read_bytes() or \
                        original_evidence.get("source_output_sha256") != hashlib.sha256(output.read_bytes()).hexdigest():
                    raise RuntimeError("Retained endpoint execution retry evidence no longer matches its source")
                retry_input = retries[-1]
                evidence = json.loads(retry_input.with_suffix(".restart.json").read_text())
            else:
                retry_input, evidence = _prepare_orca_execution_retry(
                    input_path, output=output, existing_retry_inputs=[],
                )
            execution = _run_or_resume(retry_input, timeout_seconds=timeout_seconds)
            execution["execution_startup_retry"] = evidence
        else:
            execution = _run_or_resume(input_path, timeout_seconds=timeout_seconds)
    else:
        execution = _run_or_resume(input_path, timeout_seconds=timeout_seconds)
    executed = Path(execution.get("input") or input_path)
    return {
        **execution,
        "out": Path(execution.get("output") or executed.with_suffix(".out")),
        "xyz": executed.with_suffix(".xyz"),
    }


def _optimize_endpoint(seed: Path, workdir: Path, *, charge: int, multiplicity: int,
                       ncores: int, method_keywords: list[str] | None,
                       timeout_seconds: float | None,
                       surface_resolution_policy: SurfaceResolutionPolicy | None = None) -> dict:
    surface_resolution_policy = surface_resolution_policy or SurfaceResolutionPolicy()
    temperature_contract = _frequency_temperature_kwargs(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    local_seed = workdir / "input.xyz"
    if seed.resolve() != local_seed.resolve():
        if local_seed.is_file() and local_seed.read_bytes() != seed.read_bytes():
            raise RuntimeError("Retained endpoint seed differs; a new calculation workspace is required")
        if not local_seed.is_file():
            shutil.copy2(seed, local_seed)
    contract_path = workdir / "endpoint-contract.json"
    contract_value = {
        "schema_version": 2,
        "kind": "orca_endpoint_contract",
        "input_xyz_sha256": hashlib.sha256(local_seed.read_bytes()).hexdigest(),
        "charge": int(charge),
        "multiplicity": int(multiplicity),
        "ncores": int(ncores),
        "method_keywords": list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]),
        **temperature_contract,
        **({"native_settings_sha256": active_native_settings().sha256}
           if active_native_settings() is not None else {}),
        "minimum_recovery": {
            "maximum_attempts": 2,
            "rms_displacement_angstrom": 0.15,
            "imaginary_frequency_threshold_cm-1": surface_resolution_policy.meaningful_imaginary_frequency_cm_1,
        },
        "stationary_point_recovery_policy": {
            "maximum_attempts": surface_resolution_policy.maximum_stationary_point_recoveries,
            "surface_resolution_policy_id": surface_resolution_policy.policy_id,
        },
    }
    retained_contract = json.loads(contract_path.read_text()) if contract_path.is_file() else None
    if temperature_contract and isinstance(retained_contract, dict) and (
        retained_contract.get("thermochemistry_temperature_K") != temperature_contract["thermochemistry_temperature_K"]
    ):
        raise RuntimeError("Retained endpoint thermochemistry temperature differs; a new workspace is required")
    optimized = workdir / "opt.xyz"
    frequency = workdir / "freq.out"
    optimization_output = workdir / "opt.out"
    # This cache branch does not call the executor, so it must independently
    # recheck contracts that bind retained job inputs and referenced geometry.
    for retained_output in (optimization_output, frequency):
        if _normal_orca_output(retained_output):
            _validate_retained_orca_job_input(retained_output.with_suffix(".inp"))
    contract_core_fields = (
        "kind", "input_xyz_sha256", "charge", "multiplicity", "ncores", "method_keywords",
    )
    if active_native_settings() is not None:
        contract_core_fields += ("native_settings_sha256",)
    if temperature_contract:
        contract_core_fields += ("thermochemistry_temperature_K",)
        if _normal_orca_output(frequency):
            validated_orca_temperature_thermochemistry(frequency, temperature_contract["thermochemistry_temperature_K"])
    retained_contract_matches = (
        isinstance(retained_contract, dict)
        and all(retained_contract.get(key) == contract_value.get(key) for key in contract_core_fields)
    )
    cached_base = (
        retained_contract_matches
        and optimized.is_file()
        and _normal_orca_output(optimization_output)
        and _orca_optimization_converged(optimization_output)
        and _normal_orca_output(frequency)
    )
    retained_base = _retained_endpoint_workflow(
        workdir,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        method_keywords=method_keywords,
    )
    if cached_base:
        input_elements, _, _ = read_xyz(local_seed)
        optimized_elements, _, _ = read_xyz(optimized)
        if input_elements != optimized_elements:
            raise RuntimeError("Cached ORCA endpoint changed composition or atom order")
        result = {
            "optimized_xyz": optimized,
            "optimization": {"out": optimization_output},
            "frequency": {"out": frequency},
            "frequency_check": frequency_stability_check(frequency),
            "resumed": True,
        }
        if retained_contract != contract_value:
            # Schema-v1 endpoint contracts describe the same electronic job.
            # Upgrade the small contract in place after validating every
            # retained output instead of invalidating costly ORCA evidence.
            _write_json(contract_path, contract_value)
    elif retained_base is not None:
        result = retained_base
        input_elements, _, _ = read_xyz(local_seed)
        optimized_elements, _, _ = read_xyz(Path(result["optimized_xyz"]))
        if input_elements != optimized_elements:
            raise RuntimeError("Retained ORCA endpoint changed composition or atom order")
        _write_json(contract_path, {
            **contract_value,
            "retained_legacy_evidence_adopted": True,
            "retained_metadata": str(workdir / "metadata.json"),
        })
    else:
        timeout_seconds = _surface_budget_gate(
            workdir,
            requested_calculations=2,
            ncores=ncores,
            timeout_seconds=timeout_seconds,
        )
        result = run_optimization_and_frequency(
            local_seed,
            workdir,
            charge=charge,
            multiplicity=multiplicity,
            ncores=ncores,
            run_frequency=True,
            method_keywords=method_keywords,
            timeout_seconds=timeout_seconds,
            job_runner=_run_endpoint_job,
            **temperature_contract,
        )
        input_elements, _, _ = read_xyz(local_seed)
        optimized_elements, _, _ = read_xyz(Path(result["optimized_xyz"]))
        if input_elements != optimized_elements:
            raise RuntimeError("ORCA optimized geometry changed endpoint composition or atom order")
        _write_json(contract_path, contract_value)

    if not _orca_optimization_converged(Path(result["optimization"]["out"])):
        return _run_stationary_point_recovery(
            result,
            workdir,
            charge=charge,
            multiplicity=multiplicity,
            ncores=ncores,
            method_keywords=method_keywords,
            timeout_seconds=timeout_seconds,
            policy=surface_resolution_policy,
        )

    if "frequency_check" not in result:
        raise RuntimeError("A converged endpoint lacks retained frequency evidence")
    if temperature_contract:
        validated_orca_temperature_thermochemistry(
            Path(result["frequency"]["out"]), temperature_contract["thermochemistry_temperature_K"],
        )

    if result["frequency_check"].get("IsMinimum"):
        return result
    return _recover_endpoint_minimum(
        result, workdir, charge=charge, multiplicity=multiplicity, ncores=ncores,
        method_keywords=method_keywords, timeout_seconds=timeout_seconds,
        imaginary_frequency_threshold_cm_1=surface_resolution_policy.meaningful_imaginary_frequency_cm_1,
    )


def _imaginary_mode_displacement_seed(
    optimized_xyz: Path,
    hessian: Path,
    output: Path,
    *,
    direction: int,
    rms_displacement_angstrom: float = 0.15,
    imaginary_frequency_threshold_cm_1: float = -20.0,
) -> dict:
    """Displace a stationary point along its most negative normal mode."""
    if direction not in {-1, 1}:
        raise ValueError("Imaginary-mode displacement direction must be -1 or +1")
    modes = parse_orca_normal_modes(Path(hessian))
    if imaginary_frequency_threshold_cm_1 >= 0.0:
        raise ValueError("Imaginary-frequency threshold must be negative")
    imaginary = np.flatnonzero(
        modes.frequencies_cm_1 < float(imaginary_frequency_threshold_cm_1)
    )
    if not len(imaginary):
        raise ValueError("ORCA Hessian contains no significant imaginary mode to follow")
    mode_index = int(imaginary[np.argmin(modes.frequencies_cm_1[imaginary])])
    elements, coordinates, _ = read_xyz(Path(optimized_xyz))
    if tuple(elements) != modes.atom_labels:
        raise ValueError("ORCA Hessian atom order does not match the optimized endpoint")
    vector = np.asarray(modes.vectors[mode_index], dtype=float)
    vector_rms = float(np.sqrt(np.mean(np.sum(vector * vector, axis=1))))
    if not np.isfinite(vector_rms) or vector_rms <= 0.0:
        raise ValueError("ORCA imaginary mode has no finite Cartesian displacement")
    displacement = direction * float(rms_displacement_angstrom) * vector / vector_rms
    write_xyz(
        output, elements, coordinates + displacement,
        f"Minimum recovery from mode {mode_index} ({modes.frequencies_cm_1[mode_index]:.3f} cm-1), direction {direction:+d}",
    )
    return {
        "mode_index": mode_index,
        "frequency_cm-1": float(modes.frequencies_cm_1[mode_index]),
        "direction": int(direction),
        "rms_displacement_angstrom": float(rms_displacement_angstrom),
        "seed_xyz": str(output),
    }


def _recover_endpoint_minimum(
    initial_result: dict,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    imaginary_frequency_threshold_cm_1: float = -20.0,
) -> dict:
    """Try both directions of the lowest imaginary mode before failing closed."""
    workdir = Path(workdir)
    recovery_dir = workdir / "minimum-recovery"
    manifest = recovery_dir / "minimum-recovery.json"
    source_xyz = Path(initial_result["optimized_xyz"])
    hessian = Path(initial_result["frequency"]["out"]).with_suffix(".hess")
    attempts: list[dict] = []
    results: list[dict] = [initial_result]
    selected_result: dict | None = None
    temperature_contract = _frequency_temperature_kwargs(workdir)
    if temperature_contract:
        for retained_frequency in recovery_dir.glob("attempt-*/freq.out"):
            if _normal_orca_output(retained_frequency):
                validated_orca_temperature_thermochemistry(retained_frequency, temperature_contract["thermochemistry_temperature_K"])
    for position, direction in enumerate((1, -1), start=1):
        attempt_dir = recovery_dir / f"attempt-{position:02d}"
        attempt_dir.mkdir(parents=True, exist_ok=True)
        seed = attempt_dir / "input.xyz"
        attempt: dict[str, Any] = {"attempt": position, "direction": direction, "status": "started"}
        try:
            displacement = _imaginary_mode_displacement_seed(
                source_xyz, hessian, seed, direction=direction,
                imaginary_frequency_threshold_cm_1=imaginary_frequency_threshold_cm_1,
            )
            attempt["displacement"] = displacement
            cached_optimized = attempt_dir / "opt.xyz"
            cached_optimization = attempt_dir / "opt.out"
            cached_frequency = attempt_dir / "freq.out"
            if (
                cached_optimized.is_file()
                and _normal_orca_output(cached_optimization)
                and _orca_optimization_converged(cached_optimization)
                and _normal_orca_output(cached_frequency)
            ):
                calculation = {
                    "optimized_xyz": cached_optimized,
                    "optimization": {"out": cached_optimization},
                    "frequency": {"out": cached_frequency},
                    "frequency_check": frequency_stability_check(cached_frequency),
                    "resumed": True,
                }
            else:
                bounded_timeout = _surface_budget_gate(
                    attempt_dir,
                    requested_calculations=2,
                    ncores=ncores,
                    timeout_seconds=timeout_seconds,
                )
                calculation = run_optimization_and_frequency(
                    seed, attempt_dir, charge=charge, multiplicity=multiplicity,
                    ncores=ncores, run_frequency=True, method_keywords=method_keywords,
                    timeout_seconds=bounded_timeout,
                    job_runner=_run_endpoint_job,
                    **temperature_contract,
                )
            if temperature_contract:
                validated_orca_temperature_thermochemistry(
                    Path(calculation["frequency"]["out"]), temperature_contract["thermochemistry_temperature_K"],
                )
            if not _orca_optimization_converged(Path(calculation["optimization"]["out"])):
                raise RuntimeError("ORCA recovery optimization did not reach its convergence marker")
            source_elements, _, _ = read_xyz(source_xyz)
            recovered_elements, _, _ = read_xyz(Path(calculation["optimized_xyz"]))
            if source_elements != recovered_elements:
                raise RuntimeError("ORCA recovery optimization changed endpoint composition or atom order")
            results.append(calculation)
            check = calculation["frequency_check"]
            attempt.update(
                status="validated_minimum" if check.get("IsMinimum") else "nonminimum_endpoint",
                optimized_xyz=str(calculation["optimized_xyz"]),
                optimization_output=str(calculation["optimization"]["out"]),
                frequency_output=str(calculation["frequency"]["out"]),
                frequency_check=check,
            )
            if check.get("IsMinimum"):
                selected_result = calculation
        except Exception as error:
            if isinstance(error, (SurfaceResolutionBudgetExhausted, TransferControlError)):
                raise
            attempt.update(status="failed", failure_reason=str(error))
        attempts.append(attempt)
        _write_json(manifest, {
            "schema_version": 1,
            "kind": "orca_endpoint_minimum_recovery",
            "initial_frequency_check": initial_result.get("frequency_check"),
            "attempts": attempts,
            "status": "validated_minimum" if selected_result else "running",
        })
        if selected_result is not None:
            break

    if selected_result is None:
        selected_result = min(
            results,
            key=lambda item: (
                int((item.get("frequency_check") or {}).get("NumImaginary", 10**6)),
                abs(float((item.get("frequency_check") or {}).get("LowestFrequency", -10**6))),
            ),
        )
    recovery = {
        "schema_version": 1,
        "kind": "orca_endpoint_minimum_recovery",
        "initial_frequency_check": initial_result.get("frequency_check"),
        "attempts": attempts,
        "status": "validated_minimum" if selected_result["frequency_check"].get("IsMinimum") else "recovery_exhausted",
        "selected_optimized_xyz": str(selected_result["optimized_xyz"]),
    }
    _write_json(manifest, recovery)
    return {**selected_result, "minimum_recovery": recovery}


def _render_preview(route: RouteSpec, endpoints: dict, visual_dir: Path, storyboard: Path) -> dict:
    trajectory = interpolate_xyz(
        Path(endpoints["reactant_xyz"]),
        Path(endpoints["product_xyz"]),
        visual_dir / "candidate-path-preview.xyz",
        comment_prefix="Candidate interpolation; not an ORCA reaction path",
    )
    animation = render_xyz_trajectory_gif(
        trajectory,
        visual_dir / "candidate-path-preview.gif",
        title=f"{route.route_id}: candidate geometry interpolation",
        evidence_label="CANDIDATE PREVIEW — NOT A COMPUTED REACTION PATH",
    )
    manifest = write_visual_manifest(
        visual_dir / "visuals.json",
        route_id=route.route_id,
        evidence_level="candidate_preview",
        storyboard=storyboard,
        animation=animation,
        trajectory=trajectory,
    )
    return {"storyboard": str(storyboard), "animation": str(animation),
            "trajectory": str(trajectory), "manifest": str(manifest)}


def _run_dissociation(route: RouteSpec, explanation_dir: Path, *, prepare_only: bool,
                      ncores: int, method_keywords: list[str] | None,
                      timeout_seconds: float | None, scan_steps: int) -> dict:
    endpoint_dir = explanation_dir / "endpoints"
    endpoint_dir.mkdir(parents=True, exist_ok=True)
    seed = endpoint_dir / "parent-seed.xyz"
    if not seed.is_file():
        smiles_to_xyz(route.parent_smiles, seed)
    if prepare_only:
        reactant = seed
    else:
        endpoint_result = _optimize_endpoint(
            seed, endpoint_dir / "reactant", charge=route.charge, multiplicity=route.multiplicity,
            ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
        )
        if not endpoint_result["frequency_check"]["IsMinimum"]:
            raise RuntimeError("The parent endpoint is not an ORCA local minimum")
        reactant = Path(endpoint_result["optimized_xyz"])
    endpoints = build_dissociation_endpoint(
        route, reactant, endpoint_dir / "product-mapped.xyz"
    )
    if prepare_only:
        return {"status": "prepared", "endpoints": endpoints}

    path_dir = explanation_dir / "path"
    path_dir.mkdir(parents=True, exist_ok=True)
    path_reactant = path_dir / "reactant.xyz"
    shutil.copy2(reactant, path_reactant)
    start = endpoints["initial_distance_angstrom"]
    end = endpoints["final_distance_angstrom"]
    scan_input = create_orca_relaxed_scan_input(
        path_reactant,
        bond_atom_indices=route.broken_bonds[0],
        start_distance_angstrom=start,
        end_distance_angstrom=end,
        charge=route.charge,
        multiplicity=route.multiplicity,
        ncores=ncores,
        steps=scan_steps,
        method_keywords=method_keywords,
    )
    scan_files = _scan_frames(path_dir, scan_input.stem)
    abnormal_reason = None
    if len(scan_files) >= scan_steps and scan_input.with_suffix(".out").is_file():
        normal = _normal_orca_output(scan_input.with_suffix(".out"))
        execution = {
            "status": "completed" if normal else "incomplete",
            "resumed": True,
            "input": str(scan_input),
            "output": str(scan_input.with_suffix(".out")),
        }
        if not normal:
            abnormal_reason = "ORCA retained every requested scan geometry but did not terminate normally"
    else:
        try:
            execution = _run_or_resume(scan_input, timeout_seconds=timeout_seconds)
        except Exception as error:
            if isinstance(error, (SurfaceResolutionBudgetExhausted, TransferControlError)):
                raise
            abnormal_reason = str(error)
            execution = {
                "status": "incomplete",
                "resumed": False,
                "input": str(scan_input),
                "output": str(scan_input.with_suffix(".out")),
                "failure_reason": abnormal_reason,
            }
        executed_scan_input = Path(execution.get("input") or scan_input)
        scan_files = _scan_frames(path_dir, executed_scan_input.stem)
    if len(scan_files) < 3:
        raise RuntimeError(abnormal_reason or "ORCA completed the relaxed scan without at least three retained scan geometries")
    trajectory = _combine_xyz_files(scan_files, path_dir / "trajectory.xyz")
    validation = _validate_dissociation_path(trajectory, route.broken_bonds[0], end)
    if not validation["valid"]:
        raise RuntimeError("The retained ORCA scan does not follow the declared dissociation coordinate")
    energies = _scan_energies(
        path_dir, Path(execution.get("input") or scan_input).stem, len(scan_files), scan_files,
    )
    complete = abnormal_reason is None and len(scan_files) >= scan_steps
    return {
        "status": "computed_dissociation_path" if complete else "incomplete_dissociation_path",
        "endpoints": endpoints,
        "execution": execution,
        "trajectory": str(trajectory),
        "energies_hartree": energies,
        "validation": validation,
        "failure_reason": abnormal_reason,
        "interpretation": (
            "A relaxed ground-state bond-distance scan was retained. No transition state or excited-state "
            "photochemical dynamics is claimed. An incomplete label means at least one requested point or "
            "ORCA normal termination/convergence check was missing."
        ),
    }


def classify_path_energy_profile(
    energies_hartree: list[float] | None,
    *,
    reactant_fragment_count: int,
    product_fragment_count: int,
    barrier_tolerance_hartree: float = 5.0e-4,
) -> dict:
    """Classify a complete path-energy profile without assuming a family.

    This is a hook-level physical classification, not a rate model.  In
    particular, a downhill association curve still requires capture and
    recrossing kinetics before it can contribute to a lifetime.
    """
    import math

    if not energies_hartree or len(energies_hartree) < 5:
        return {"status": "surface_unresolved", "reason": "fewer_than_five_energy_points"}
    energies = [float(value) for value in energies_hartree]
    if any(not math.isfinite(value) for value in energies):
        return {"status": "surface_unresolved", "reason": "nonfinite_energy"}
    if barrier_tolerance_hartree <= 0 or not math.isfinite(barrier_tolerance_hartree):
        raise ValueError("Barrier tolerance must be finite and positive")
    largest_step = max(abs(right - left) for left, right in zip(energies, energies[1:]))
    if largest_step > 0.2:
        return {
            "status": "surface_unresolved",
            "reason": "discontinuous_energy_profile",
            "largest_adjacent_step_hartree": largest_step,
        }
    interior_maximum = max(energies[1:-1])
    endpoint_maximum = max(energies[0], energies[-1])
    if interior_maximum > endpoint_maximum + barrier_tolerance_hartree:
        maximum_index = max(range(1, len(energies) - 1), key=energies.__getitem__)
        return {
            "status": "barriered_ts_candidate",
            "interior_maximum_index": maximum_index,
            "interior_barrier_hartree": interior_maximum - energies[0],
            "largest_adjacent_step_hartree": largest_step,
        }
    direction_change = energies[-1] - energies[0]
    if reactant_fragment_count > product_fragment_count and direction_change <= barrier_tolerance_hartree:
        status = "barrierless_capture_candidate"
    elif reactant_fragment_count < product_fragment_count and direction_change >= -barrier_tolerance_hartree:
        status = "barrierless_dissociation_candidate"
    elif reactant_fragment_count == product_fragment_count and direction_change <= barrier_tolerance_hartree:
        # A two-to-two profile can be barrierless only in the asserted forward
        # direction.  Calling an uphill encounter path "barrierless" would
        # incorrectly send it to a capture-rate model.
        status = "barrierless_path_candidate"
    else:
        status = "surface_unresolved"
    downhill_direction = (
        "forward" if direction_change < -barrier_tolerance_hartree
        else "reverse" if direction_change > barrier_tolerance_hartree
        else "approximately_thermoneutral"
    )
    return {
        "status": status,
        "reaction_energy_hartree": direction_change,
        "downhill_direction": downhill_direction,
        "forward_energetic_threshold_hartree": max(0.0, direction_change),
        "forward_target_loss_barrierless": bool(direction_change <= barrier_tolerance_hartree),
        "largest_adjacent_step_hartree": largest_step,
        "reason": (
            None if status != "surface_unresolved" else
            "forward_path_uphill_without_interior_ts_candidate"
            if reactant_fragment_count == product_fragment_count and direction_change > barrier_tolerance_hartree
            else "profile_direction_inconsistent_with_stoichiometry"
        ),
    }


def _trajectory_energies(path: Path) -> list[float] | None:
    values: list[float] = []
    for _, _, comment in xyz_frames(path):
        matches = re.findall(
            r"(?:energy|\bE)\s*(?:=|:)?\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
            comment,
            re.IGNORECASE,
        )
        if not matches:
            return None
        values.append(float(matches[-1]))
    return values


def _locate_neb_trajectory(workdir: Path, label: str) -> Path | None:
    candidates = [
        Path(workdir) / f"{label}_MEP_trj.xyz",
        Path(workdir) / f"{label}_NEB_trj.xyz",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _create_ts_optimization_input(
    seed_xyz: Path,
    *,
    charge: int,
    multiplicity: int,
    ncores: int,
    method_keywords: list[str] | None,
) -> Path:
    """Create an explicit ORCA OptTS refinement with a fresh Hessian."""
    if ncores < 1:
        raise ValueError("ncores must be at least one")
    seed_xyz = Path(seed_xyz)
    method = list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"])
    output = seed_xyz.parent / "ts-opt.inp"
    input_text = apply_native_job_settings(
        "! " + " ".join([*method, "TightSCF", "OptTS"]) + "\n"
        + f"%pal\n  nprocs {ncores}\nend\n"
        + "%geom\n  Calc_Hess true\n  Recalc_Hess 5\nend\n"
        + f"* xyzfile {charge} {multiplicity} {seed_xyz.name}\n"
    )
    if active_native_settings() is not None and output.exists() and output.read_text() != input_text:
        raise RuntimeError("Retained TS optimization input differs from the frozen native settings")
    output.write_text(input_text)
    return output


def retained_ts_source_evidence(
    route: RouteSpec, transition_state: dict, irc: dict,
    *, method_keywords: list[str] | None,
) -> dict:
    """Bind reusable electronic calculations to their original route and trace.

    This records candidate source evidence only. It deliberately cannot admit
    a new reaction, its product states, or a rate.
    """
    base = {"kind": "retained_ts_source_evidence", "schema_version": 1,
            "status": "unresolved", "valid": False, "parent_route_id": route.route_id,
            "environment": route.environment, "charge": route.charge,
            "multiplicity": route.multiplicity,
            "method_keywords": list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]),
            "claims_prohibited": ["verified_new_path", "validated_product_states", "rate", "lifetime"]}
    try:
        optimization = transition_state["optimization"]
        frequency = transition_state["frequency_execution"]
        execution = irc["execution"]
        paths = {
            "ts_input": Path(optimization["input"]), "ts_output": Path(optimization["output"]),
            "ts_xyz": Path(transition_state["optimized_xyz"]),
            "frequency_input": Path(frequency["input"]),
            "frequency_output": Path(transition_state["frequency_output"]),
            "hessian": Path(frequency["input"]).with_suffix(".hess"),
            "irc_input": Path(execution["input"]), "irc_output": Path(execution["output"]),
            "trajectory": Path(irc["trajectory"]),
        }
        references = {key: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                      for key, path in paths.items()}
        expected_method = " ".join(base["method_keywords"]).lower().split()
        for key, job in (("ts_input", "optts"), ("frequency_input", "freq"), ("irc_input", "irc")):
            text = paths[key].read_text()
            keyword_lines = re.findall(r"^\s*!\s*(.*)$", text, re.M)
            actual = " ".join(keyword_lines).lower().split()
            if job not in actual or sorted(token for token in actual if token not in {job, "tightscf"}) != sorted(
                    token for token in expected_method if token != "tightscf"):
                raise ValueError("Retained electronic-structure method differs from requested method")
            states = re.findall(r"^\s*\*\s+xyz(?:file)?\s+(-?\d+)\s+(\d+)", text, re.M | re.I)
            if states != [(str(route.charge), str(route.multiplicity))]:
                raise ValueError("Retained charge or multiplicity differs from route")
        if not _normal_orca_output(paths["ts_output"]) or not _orca_optimization_converged(paths["ts_output"]):
            raise ValueError("Retained TS optimization did not converge normally")
        if not _normal_orca_output(paths["frequency_output"]):
            raise ValueError("Retained TS frequency calculation did not terminate normally")
        check = frequency_stability_check(paths["frequency_output"], imag_tol=20.0)
        modes = parse_orca_normal_modes(paths["hessian"])
        if (not np.isfinite(modes.frequencies_cm_1).all() or check.get("NumImaginary") != 1
                or int(np.count_nonzero(modes.frequencies_cm_1 < -20.0)) != 1):
            raise ValueError("Retained TS lacks exactly one meaningful imaginary mode")
        labels, coordinates, _ = read_xyz(paths["ts_xyz"])
        if tuple(labels) != modes.atom_labels:
            raise ValueError("TS geometry and Hessian atom order differ")
        hessian_xyz = modes.coordinates_bohr * 0.529177210903
        translation = coordinates.mean(axis=0) - hessian_xyz.mean(axis=0)
        if not np.allclose(hessian_xyz + translation, coordinates, atol=1e-5, rtol=0):
            raise ValueError("Hessian geometry differs from retained TS")
        base["hessian_geometry_alignment"] = {
            "policy": "translation_only_no_atom_permutation_or_rotation",
            "translation_angstrom": translation.tolist(),
            "maximum_residual_angstrom": float(np.max(np.abs(hessian_xyz + translation - coordinates))),
        }
        frequency_text = paths["frequency_input"].read_text()
        inline = re.search(r"^\s*\*\s+xyz\s+-?\d+\s+\d+\s*\n(.*?)^\s*\*", frequency_text, re.M | re.S | re.I)
        if inline is None:
            raise ValueError("Retained frequency geometry requires an explicit inline geometry contract")
        rows = [line.split() for line in inline[1].strip().splitlines()]
        if ([row[0] for row in rows] != labels or any(len(row) != 4 for row in rows)
                or not np.allclose([[float(v) for v in row[1:]] for row in rows], coordinates, atol=1e-5, rtol=0)):
            raise ValueError("Frequency input geometry differs from retained TS")
        frames = xyz_frames(paths["trajectory"])
        if not any(elements == labels and np.allclose(xyz, coordinates, atol=1e-5, rtol=0)
                   for elements, xyz, _ in frames):
            raise ValueError("Retained IRC does not contain the optimized TS geometry")
        completion = _orca_irc_completion(paths["irc_output"])
        if not completion["valid"]:
            assembly = assemble_completed_irc_branches(paths["irc_input"])
            if (not assembly["valid"] or
                    assembly["trajectory_sha256"] != references["trajectory"]["sha256"]):
                raise ValueError("Retained IRC continuation is incomplete or differs from the cited trace")
        elif paths["trajectory"].read_bytes() != Path(
            _locate_irc_trajectory(paths["irc_input"].parent, paths["irc_input"].stem)
        ).read_bytes():
            raise ValueError("Retained trace differs from the completed IRC output trajectory")
        base.update(valid=True, status="retained_ts_source_bound", references=references,
                    imaginary_frequency_threshold_cm_1=-20.0)
    except (KeyError, IndexError, ValueError, TypeError, OSError, RuntimeError) as error:
        base["reason"] = f"{type(error).__name__}: {error}"
    return base


def assess_unexpected_retained_path(
    route: RouteSpec, proposal: dict, workdir: Path, *, method_keywords: list[str] | None,
) -> dict:
    """Recheck a source TS against a new channel without executing chemistry."""
    from .route_geometry import align_mapping_to_retained_endpoints

    result = {"kind": "unexpected_retained_path_assessment", "valid": False,
              "status": "retained_path_unresolved", "route_id": route.route_id,
              "claims_prohibited": ["verified_path", "validated_product_states", "rate", "lifetime"]}
    try:
        candidate = ReactionCandidate.from_dict(proposal["native_candidates"][0])
        expected = route_spec_from_candidate_route(route.parent_smiles, candidate.as_route_dict(parent_smiles=route.parent_smiles))
        if expected.as_dict() != route.as_dict():
            raise ValueError("Requested route differs from retained native candidate")
        source = proposal["retained_ts_source"]
        if not source.get("valid") or source.get("parent_route_id") != candidate.provenance.get("parent_route_id"):
            raise ValueError("Retained TS source parent identity unresolved")
        if (source.get("environment") != route.environment or source.get("charge") != route.charge
                or source.get("multiplicity") != route.multiplicity):
            raise ValueError("Retained TS source environment or electronic state differs")
        refs = source["references"]
        for reference in refs.values():
            if hashlib.sha256(Path(reference["path"]).read_bytes()).hexdigest() != reference["sha256"]:
                raise ValueError("Retained TS source hash mismatch")
        paths = {key: Path(value["path"]) for key, value in refs.items()}
        transition_state = {
            "optimization": {"input": str(paths["ts_input"]), "output": str(paths["ts_output"])},
            "optimized_xyz": str(paths["ts_xyz"]),
            "frequency_execution": {"input": str(paths["frequency_input"])},
            "frequency_output": str(paths["frequency_output"]),
        }
        irc = {"execution": {"input": str(paths["irc_input"]), "output": str(paths["irc_output"])},
               "trajectory": str(paths["trajectory"])}
        fresh = retained_ts_source_evidence(route, transition_state, irc, method_keywords=method_keywords)
        if not fresh["valid"] or fresh["references"] != refs:
            raise ValueError("Retained TS source no longer passes source validation")
        if refs["trajectory"]["sha256"] != candidate.provenance.get("source_trajectory_sha256"):
            raise ValueError("Candidate and TS cite different IRC trajectories")
        source_endpoint = Path(candidate.provenance["source_xyz"])
        if hashlib.sha256(source_endpoint.read_bytes()).hexdigest() != candidate.provenance["source_xyz_sha256"]:
            raise ValueError("Candidate product source hash mismatch")
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        retained = workdir / "trajectory.xyz"
        if retained.exists() and retained.read_bytes() != paths["trajectory"].read_bytes():
            raise ValueError("Retained assessment workspace contains a different trajectory")
        if not retained.exists():
            shutil.copy2(paths["trajectory"], retained)
        frames = xyz_frames(retained)
        endpoints = []
        for name, frame in (("first", frames[0]), ("last", frames[-1])):
            endpoint = workdir / f"source-{name}.xyz"
            if endpoint.exists():
                labels, xyz, _ = read_xyz(endpoint)
                if labels != frame[0] or not np.allclose(xyz, frame[1], atol=1e-9, rtol=0):
                    raise ValueError("Retained assessment endpoint was changed")
            else:
                write_xyz(endpoint, *frame)
            endpoints.append(endpoint)
        labels, xyz, _ = read_xyz(source_endpoint)
        product_indices = [index for index, frame in enumerate((frames[0], frames[-1]))
                           if labels == frame[0] and np.allclose(xyz, frame[1], atol=1e-6, rtol=0)]
        if len(product_indices) != 1:
            raise ValueError("Unique candidate product endpoint is not established")
        product_index = product_indices[0]
        mapping = align_mapping_to_retained_endpoints(route, endpoints[1-product_index], endpoints[product_index])
        if not mapping.get("valid"):
            raise ValueError("New reaction cannot be mapped onto the retained path")
        changes = candidate.provenance["observed_bond_changes"]
        if changes.get("index_convention") != "zero_based_original_irc_atom_order":
            raise ValueError("Unknown observed bond-change atom convention")
        checks = mapping["trace_alignment"]
        reactant_bonds = {tuple(sorted(pair)) for pair in checks["reactant_check"]["declared_bonds"]}
        product_bonds = {tuple(sorted(pair)) for pair in checks["product_check"]["inferred_covalent_bonds"]}
        expected_changes = {"broken_bonds": reactant_bonds - product_bonds,
                            "formed_bonds": product_bonds - reactant_bonds}
        for kind, expected_pairs in expected_changes.items():
            if {tuple(pair) for pair in changes.get(kind, [])} != expected_pairs:
                raise ValueError("Retained bond changes differ from rechecked endpoint connectivity")
        coordinates = [{"atom_indices": pair, "kind": kind}
                       for kind in ("broken_bonds", "formed_bonds") for pair in changes.get(kind, [])]
        modes = parse_orca_normal_modes(paths["hessian"])
        imaginary = np.flatnonzero(modes.frequencies_cm_1 < -20.0)
        _, ts_xyz, _ = read_xyz(paths["ts_xyz"])
        mode_check = bond_change_mode_overlap(ts_xyz, modes.vectors[int(imaginary[0])], coordinates)
        result.update(mapping=mapping, mode_check=mode_check,
                      reactant_xyz=str(endpoints[1-product_index]), product_xyz=str(endpoints[product_index]),
                      source_references=refs)
        if mode_check["status"] != "matched":
            raise ValueError("New reaction bond changes do not match the imaginary mode")
        result.update(valid=True, status="retained_path_geometry_supported",
                      limitations=["Mode overlap tests at least one changed bond, not simultaneity of all changes.",
                                   "Product minima, electronic-state applicability and rate theory remain unvalidated."])
    except (KeyError, IndexError, ValueError, TypeError, OSError, RuntimeError) as error:
        result["reason"] = f"{type(error).__name__}: {error}"
    return result


def _prepare_reused_ts_irc(
    route: RouteSpec, admission: dict, workdir: Path, *, ncores: int,
    method_keywords: list[str] | None,
) -> tuple[Path, dict]:
    """Copy only hash-bound geometry/Hessian into a new owned IRC workspace."""
    if (admission.get("kind") != "completed_ts_frequency_reuse_admission"
            or admission.get("valid") is not True):
        raise ValueError("Completed TS stage reuse requires source admission")
    for reference in admission["source_artifacts"]:
        if hashlib.sha256(Path(reference["path"]).read_bytes()).hexdigest() != reference["sha256"]:
            raise ValueError("Retained TS source changed before IRC preparation")
    workdir = Path(workdir)
    ts = admission["transition_state"]
    origins = {"geometry": Path(ts["optimized_xyz"]),
               "hessian": Path(ts["frequency_execution"]["input"]).with_suffix(".hess")}
    if workdir.resolve().is_relative_to(origins["geometry"].parent.resolve()):
        raise ValueError("Reused TS IRC requires a workspace outside the source calculation")
    workdir.mkdir(parents=True, exist_ok=True)
    targets = {"geometry": workdir / "retained-ts.xyz", "hessian": workdir / "retained-ts.hess"}
    copies = {}
    for key, origin in origins.items():
        target = targets[key]
        data = origin.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        expected = next(ref["sha256"] for ref in admission["source_artifacts"]
                        if Path(ref["path"]).resolve() == origin.resolve())
        if digest != expected:
            raise ValueError("Retained TS copy source changed")
        if target.exists() and target.read_bytes() != data:
            raise ValueError("Derived IRC workspace contains conflicting retained TS evidence")
        if not target.exists():
            shutil.copy2(origin, target)
        copies[key] = {"source": str(origin), "copy": str(target), "sha256": digest}
    inp = create_orca_irc_input(targets["geometry"], charge=route.charge, multiplicity=route.multiplicity,
                                hessian_file=targets["hessian"], label="irc", ncores=ncores,
                                method_keywords=method_keywords)
    evidence = {"kind": "source_bound_ts_frequency_reuse_for_new_irc", "copies": copies,
                "source_path_result": admission["source_path_result"],
                "reused_stages": admission["reused_stages"], "irc_reused": False,
                "claims_prohibited": ["verified_path", "rate", "lifetime"]}
    _write_json(workdir / "ts-frequency-reuse.json", evidence)
    return inp, evidence


def _verify_ts_seed(
    route: RouteSpec,
    ts_seed_source: Path,
    workdir: Path,
    *,
    reactant_xyz: Path,
    product_xyz: Path,
    mapping_result: dict,
    separated_fragments: bool,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    surface_resolution_policy: SurfaceResolutionPolicy,
    intended_coordinates: list[dict] | None = None,
    retained_ts_frequency: dict | None = None,
    retained_completed_irc: dict | None = None,
) -> dict:
    """Apply the ordinary TS/frequency/mode/IRC gates to one retained seed."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    ts_seed = workdir / "ts-seed.xyz"
    if Path(ts_seed_source).resolve() != ts_seed.resolve():
        if ts_seed.is_file() and ts_seed.read_bytes() != Path(ts_seed_source).read_bytes():
            raise RuntimeError("Retained transition-state seed differs from the requested seed")
        if not ts_seed.is_file():
            shutil.copy2(ts_seed_source, ts_seed)
    result: dict[str, Any] = {
        "status": "running",
        "path_classification": "surface_unresolved",
        "seed_xyz": str(ts_seed),
        "intended_coordinates": intended_coordinates or [],
    }
    try:
        if retained_completed_irc is not None and retained_ts_frequency is None:
            raise ValueError("Completed IRC reuse also requires admitted TS/frequency evidence")
        if retained_completed_irc is not None:
            if retained_completed_irc.get("valid") is not True:
                raise ValueError("Completed IRC evidence is not admitted")
            ts_owner = retained_ts_frequency["transition_state"]
            irc_ts_owner = retained_completed_irc["transition_state"]
            paths = [("optimized_xyz",), ("frequency_output",),
                     ("optimization", "input"), ("optimization", "output"),
                     ("frequency_execution", "input")]
            for keys in paths:
                left, right = ts_owner, irc_ts_owner
                for key in keys:
                    left, right = left[key], right[key]
                if Path(left).resolve(strict=True) != Path(right).resolve(strict=True):
                    raise ValueError("Completed IRC and retained TS evidence have different owners")
        if retained_ts_frequency is None:
            ts_input = _create_ts_optimization_input(
                ts_seed, charge=route.charge, multiplicity=route.multiplicity,
                ncores=ncores, method_keywords=method_keywords,
            )
            ts_execution = _run_or_resume(ts_input, timeout_seconds=timeout_seconds)
        else:
            if retained_ts_frequency.get("valid") is not True:
                raise ValueError("Retained TS stages are not admitted")
            ts_execution = retained_ts_frequency["transition_state"]["optimization"]
            ts_input = Path(ts_execution["input"])
        executed_ts_input = Path(ts_execution.get("input") or ts_input)
        if not _orca_optimization_converged(Path(ts_execution["output"])):
            result.update(status="ts_optimization_unconverged", optimization=ts_execution)
            return result
        optimized_ts = (Path(retained_ts_frequency["transition_state"]["optimized_xyz"])
                        if retained_ts_frequency is not None else executed_ts_input.with_suffix(".xyz"))
        if not optimized_ts.is_file():
            raise RuntimeError("ORCA OptTS completed without an optimized TS geometry")
        if retained_ts_frequency is None:
            frequency_input = create_orca_input(
                optimized_ts, route.charge, route.multiplicity, freq=True, label="ts-frequency",
                ncores=ncores, method_keywords=method_keywords,
                preserve_existing=bool(_frequency_temperature_kwargs(workdir)),
                **_frequency_temperature_kwargs(workdir),
            )
            frequency_execution = _run_or_resume(frequency_input, timeout_seconds=timeout_seconds)
        else:
            frequency_execution = retained_ts_frequency["transition_state"]["frequency_execution"]
            frequency_input = Path(frequency_execution["input"])
        executed_frequency_input = Path(frequency_execution.get("input") or frequency_input)
        frequency_output = Path(
            frequency_execution.get("output") or executed_frequency_input.with_suffix(".out")
        )
        temperature_contract = _frequency_temperature_kwargs(workdir)
        if temperature_contract:
            validated_orca_temperature_thermochemistry(
                frequency_output, temperature_contract["thermochemistry_temperature_K"],
            )
        frequency_check = frequency_stability_check(
            frequency_output,
            imag_tol=abs(surface_resolution_policy.meaningful_imaginary_frequency_cm_1),
        )
        transition_state = {
            "seed_xyz": str(ts_seed),
            "optimization": ts_execution,
            "optimized_xyz": str(optimized_ts),
            "frequency_execution": frequency_execution,
            "frequency_output": str(frequency_output),
            "frequency_check": frequency_check,
        }
        result["transition_state"] = transition_state
        if frequency_check.get("NumImaginary") != 1:
            result.update(status="ts_frequency_validation_failed")
            return result
        if intended_coordinates:
            hessian = executed_frequency_input.with_suffix(".hess")
            modes = parse_orca_normal_modes(hessian)
            imaginary = np.flatnonzero(
                modes.frequencies_cm_1
                < surface_resolution_policy.meaningful_imaginary_frequency_cm_1
            )
            if len(imaginary) != 1:
                result.update(status="ts_frequency_validation_failed")
                return result
            elements, coordinates, _ = read_xyz(optimized_ts)
            if tuple(elements) != modes.atom_labels:
                raise RuntimeError("TS Hessian atom order does not match the optimized geometry")
            overlap = bond_change_mode_overlap(
                coordinates,
                modes.vectors[int(imaginary[0])],
                intended_coordinates,
            )
            transition_state["reaction_coordinate_overlap"] = overlap
            if overlap["status"] != "matched":
                result.update(status="ts_mode_mismatch")
                return result
        hessian = executed_frequency_input.with_suffix(".hess")
        if retained_completed_irc is not None:
            from .completed_irc_reuse import prepare_completed_irc_reuse
            (irc_input, irc_execution, irc_completion, trajectory, reuse) = prepare_completed_irc_reuse(
                retained_completed_irc, workdir,
            )
            result["completed_irc_reuse"] = reuse
        elif retained_ts_frequency is None:
            irc_input = create_orca_irc_input(
                optimized_ts, charge=route.charge, multiplicity=route.multiplicity,
                hessian_file=hessian if hessian.is_file() else None, label="irc",
                ncores=ncores, method_keywords=method_keywords,
            )
        else:
            irc_input, reuse = _prepare_reused_ts_irc(
                route, retained_ts_frequency, workdir, ncores=ncores, method_keywords=method_keywords,
            )
            result["ts_frequency_reuse"] = reuse
        if retained_completed_irc is None:
            irc_execution = _run_or_resume(irc_input, timeout_seconds=timeout_seconds)
            irc_completion = _orca_irc_completion(Path(irc_execution["output"]))
            trajectory = _locate_irc_trajectory(
                workdir, Path(irc_execution.get("input") or irc_input).stem,
            )
        result["irc"] = {
            "execution": irc_execution,
            "completion": irc_completion,
            "trajectory": str(trajectory) if trajectory else None,
        }
        if not irc_completion["valid"] and (
            surface_resolution_policy.policy_id in {
                PRE_RING_CONNECTIVITY_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                PRE_GEOMETRY_AUDIT_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                PRE_REACHED_MINIMA_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
                GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
            }
            and surface_resolution_policy.irc_continuation_policy_id == "irc-downhill-single-extension-v1"
            and irc_completion["normal_termination"]
            and set(irc_completion["directions"]) == {"forward", "backward"}
            and all(item["valid"] or item["iteration_limit_reached"]
                    for item in irc_completion["directions"].values())
        ):
            extensions = []
            result["irc"]["branch_continuations"] = extensions
            executed_irc_input = Path(irc_execution.get("input") or irc_input)
            for branch in ("forward", "backward"):
                if not irc_completion["directions"][branch]["valid"]:
                    extensions.append(run_irc_branch_continuation(
                        executed_irc_input, branch=branch, timeout_seconds=timeout_seconds,
                    ))
            assembly = assemble_completed_irc_branches(executed_irc_input)
            result["irc"]["branch_assembly"] = assembly
            if assembly["valid"]:
                trajectory = Path(assembly["trajectory"])
                irc_completion = {
                    "kind": "continued_bidirectional_irc_completion", "valid": True,
                    "original_completion": irc_completion,
                    "assembly": assembly,
                }
        if not irc_completion["valid"]:
            result.update(status="irc_path_incomplete")
            return result
        if trajectory is None:
            raise RuntimeError("ORCA IRC completed without a retained bidirectional trajectory")
        retained = (trajectory if (result["irc"].get("branch_assembly") or {}).get("valid")
                    else workdir / "validated-trajectory.xyz")
        if trajectory.resolve() != retained.resolve():
            shutil.copy2(trajectory, retained)
        trajectory_validation = _validate_irc_route_endpoints(
            route,
            retained,
            reactant_xyz,
            product_xyz,
            mapping_result=mapping_result,
            separated_fragments=separated_fragments,
            rmsd_tolerance_angstrom=surface_resolution_policy.endpoint_rmsd_tolerance_angstrom,
        )
        result["irc"] = {
            **result["irc"],
            "execution": irc_execution,
            "completion": irc_completion,
            "trajectory": str(retained),
            "endpoint_validation": trajectory_validation,
        }
        if not trajectory_validation["valid"]:
            if trajectory_validation.get("status") == "irc_endpoint_identity_unresolved":
                result.update(status="irc_endpoint_identity_unresolved")
                return result
            result["unexpected_endpoint_discovery"] = propose_unexpected_irc_candidate(
                route, trajectory_validation, {
                    "valid": irc_completion["valid"], "trajectory": str(retained),
                    "trajectory_sha256": hashlib.sha256(retained.read_bytes()).hexdigest(),
                },
            )
            if result["unexpected_endpoint_discovery"].get("status") == "candidate_proposed":
                result["unexpected_endpoint_discovery"]["retained_ts_source"] = retained_ts_source_evidence(
                    route, transition_state, result["irc"], method_keywords=method_keywords,
                )
            result.update(status="irc_endpoint_validation_failed")
            return result
        result.update(
            status="verified_transition_state_path",
            path_classification="verified_barriered_path",
            trajectory=str(retained),
            endpoint_reference_geometries={"reactants": str(reactant_xyz), "products": str(product_xyz)},
            path_evidence={
                "kind": "orca_ts_frequency_irc",
                "one_meaningful_imaginary_frequency": True,
                "intended_mode_match": (
                    transition_state.get("reaction_coordinate_overlap", {}).get("status") == "matched"
                    if intended_coordinates else None
                ),
                "irc_matches_declared_endpoints": True,
                "trajectory": str(retained),
            },
        )
        return result
    except Exception as error:
        if isinstance(error, (SurfaceResolutionBudgetExhausted, TransferControlError)):
            raise
        result.update(
            status="path_execution_failed",
            failure_reason=f"{type(error).__name__}: {error}",
        )
        return result


def _run_bond_change_guided_ts_recovery(
    route: RouteSpec,
    *,
    reactant_xyz: Path,
    product_xyz: Path,
    bond_changes: dict,
    mapping_result: dict,
    path_dir: Path,
    separated_fragments: bool,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    surface_resolution_policy: SurfaceResolutionPolicy,
) -> dict:
    """Run one bounded simultaneous scan, then the ordinary TS/IRC gates."""
    recovery_root = Path(path_dir) / "bond-change-guided-recovery"
    recovery_root.mkdir(parents=True, exist_ok=True)
    reactant_elements, reactant_coordinates, _ = read_xyz(reactant_xyz)
    product_elements, product_coordinates, _ = read_xyz(product_xyz)

    def finish_preflight(value: dict[str, Any]) -> dict[str, Any]:
        value["surface_resolution_attempt_evidence"] = _write_surface_attempt_evidence(
            recovery_root / "preflight",
            strategy="bond_change_guided_scan",
            status=str(value.get("status") or "coordinate_scan_no_ts_seed"),
            result=value,
            parameters={
                "phase": "scan_plan_preflight",
                "charge": int(route.charge),
                "multiplicity": int(route.multiplicity),
                "ncores": int(ncores),
                "method_keywords": list(
                    method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]
                ),
                "surface_resolution_policy": surface_resolution_policy.as_dict(),
            },
            input_paths=[Path(reactant_xyz), Path(product_xyz)],
            failure_reason=value.get("reason"),
        )
        return value

    if reactant_elements != product_elements:
        return finish_preflight({
            "status": "coordinate_scan_no_ts_seed",
            "reason": "endpoint_atom_order_mismatch",
            "rate_claim_allowed": False,
        })
    try:
        plan = build_bond_change_scan_plan(
            reactant_coordinates,
            product_coordinates,
            bond_changes,
            policy=surface_resolution_policy,
            requested_points=min(12, surface_resolution_policy.maximum_scan_points),
        )
    except Exception as error:
        return finish_preflight({
            "status": "coordinate_scan_no_ts_seed",
            "reason": f"{type(error).__name__}: {error}",
            "rate_claim_allowed": False,
        })
    if plan["status"] != "prepared":
        return finish_preflight({**plan, "rate_claim_allowed": False})
    identity = {
        "schema_version": 1,
        "kind": "orca_bond_change_guided_ts_recovery_contract",
        "route_id": route.route_id,
        "reactant_xyz_sha256": hashlib.sha256(Path(reactant_xyz).read_bytes()).hexdigest(),
        "product_xyz_sha256": hashlib.sha256(Path(product_xyz).read_bytes()).hexdigest(),
        "charge": route.charge,
        "multiplicity": route.multiplicity,
        "ncores": ncores,
        "method_keywords": list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]),
        "surface_resolution_policy_id": surface_resolution_policy.policy_id,
        "scan_plan_sha256": plan["scientific_hash"],
    }
    identity_hash = hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    attempt_dir = None
    for retained in sorted(recovery_root.glob("attempt-*/attempt-contract.json")):
        value = json.loads(retained.read_text())
        if value.get("scientific_identity_sha256") == identity_hash:
            attempt_dir = retained.parent
            break
    if attempt_dir is None:
        number = 1 + max(
            [int(path.name.split("-")[-1]) for path in recovery_root.glob("attempt-*") if path.is_dir()]
            or [0]
        )
        attempt_dir = recovery_root / f"attempt-{number:03d}"
        attempt_dir.mkdir(parents=True, exist_ok=False)
    contract = {**identity, "scientific_identity_sha256": identity_hash, "scan_plan": plan}
    contract_path = attempt_dir / "attempt-contract.json"
    if contract_path.is_file() and json.loads(contract_path.read_text()) != contract:
        raise RuntimeError("Retained bond-change-guided recovery contract was mutated")
    _write_json(contract_path, contract)
    local_reactant = attempt_dir / "reactant.xyz"
    if not local_reactant.is_file():
        shutil.copy2(reactant_xyz, local_reactant)
    result: dict[str, Any] = {
        "schema_version": 1,
        "kind": "orca_bond_change_guided_ts_recovery",
        "status": "running",
        "contract": str(contract_path),
        "scan_plan": plan,
        "rate_claim_allowed": False,
        "lifetime_claim_allowed": False,
    }
    result_path = attempt_dir / "result.json"

    def finish() -> dict[str, Any]:
        scan_input_value = ((result.get("scan") or {}).get("execution") or {}).get("input")
        result["surface_resolution_attempt_evidence"] = _write_surface_attempt_evidence(
            attempt_dir,
            strategy="bond_change_guided_scan",
            status=str(result.get("status") or "coordinate_scan_no_ts_seed"),
            result=result,
            parameters={
                "scan_plan": plan,
                "charge": int(route.charge),
                "multiplicity": int(route.multiplicity),
                "ncores": int(ncores),
                "method_keywords": list(
                    method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]
                ),
                "surface_resolution_policy": surface_resolution_policy.as_dict(),
            },
            input_paths=[
                Path(reactant_xyz), Path(product_xyz), contract_path,
                *([Path(str(scan_input_value))] if scan_input_value else []),
            ],
            failure_reason=result.get("reason"),
            claims_permitted=(
                ("verified_transition_state_path",)
                if result.get("status") == "verified_transition_state_path" else ()
            ),
        )
        _write_json(result_path, result)
        return result

    try:
        scan_input = create_orca_guided_scan_input(
            local_reactant,
            coordinates=plan["coordinates"],
            charge=route.charge,
            multiplicity=route.multiplicity,
            label="guided-scan",
            ncores=ncores,
            steps=plan["scan_points"],
            method_keywords=method_keywords,
        )
        execution = _run_or_resume(scan_input, timeout_seconds=timeout_seconds)
        executed_scan_input = Path(execution.get("input") or scan_input)
        scan_files = _scan_frames(attempt_dir, executed_scan_input.stem)
        energies = _scan_energies(
            attempt_dir, executed_scan_input.stem, len(scan_files), scan_files,
        )
        result["scan"] = {
            "execution": execution,
            "geometry_count": len(scan_files),
            "geometry_files": [str(path) for path in scan_files],
            "energies_hartree": energies,
        }
        if execution.get("status") != "completed" or len(scan_files) != plan["scan_points"]:
            result.update(
                status="coordinate_scan_no_ts_seed",
                reason="guided_scan_incomplete_or_abnormally_terminated",
            )
            return finish()
        selection = select_interior_scan_maximum(
            energies,
            tolerance_hartree=surface_resolution_policy.energy_extremum_tolerance_hartree,
        )
        result["ts_seed_selection"] = selection
        if selection["status"] != "interior_ts_seed_candidate":
            result.update(status="coordinate_scan_no_ts_seed", reason=selection.get("reason"))
            return finish()
        seed_source = scan_files[int(selection["index"])]
        ts_result = _verify_ts_seed(
            route,
            seed_source,
            attempt_dir / "transition-state",
            reactant_xyz=reactant_xyz,
            product_xyz=product_xyz,
            mapping_result=mapping_result,
            separated_fragments=separated_fragments,
            ncores=ncores,
            method_keywords=method_keywords,
            timeout_seconds=timeout_seconds,
            surface_resolution_policy=surface_resolution_policy,
            intended_coordinates=plan["coordinates"],
        )
        result["transition_state_verification"] = ts_result
        result.update(
            status=ts_result["status"],
            path_classification=ts_result.get("path_classification", "surface_unresolved"),
            transition_state=ts_result.get("transition_state"),
            irc=ts_result.get("irc"),
            trajectory=ts_result.get("trajectory"),
            path_evidence=ts_result.get("path_evidence"),
        )
        return finish()
    except Exception as error:
        if isinstance(error, TransferControlError):
            raise
        scan_output = attempt_dir / "guided-scan.out"
        infrastructure_failure = _orca_startup_infrastructure_failure(scan_output)
        result.update(
            status=(
                "execution_infrastructure_failure"
                if infrastructure_failure else "coordinate_scan_no_ts_seed"
            ),
            reason=f"{type(error).__name__}: {error}",
            infrastructure_failure={
                "kind": "orca_startup_infrastructure_execution_failure",
                "scientific_path_failure": False,
                "failed_output": str(scan_output),
            } if infrastructure_failure else None,
        )
        return finish()


def _atomic_stationary_fragment(
    seed: Path, folder: Path, *, smiles: str, charge: int, multiplicity: int,
    ncores: int, method_keywords: list[str] | None, timeout_seconds: float | None,
    execute: bool,
) -> dict:
    """Validate an atom with a controlled energy job, never an invented Hessian."""
    evidence_path = Path(folder) / "atomic-stationary.json"
    if execute:
        native_input = create_orca_input(
            Path(seed), charge, multiplicity, opt=False, freq=False, label="atom-energy",
            ncores=ncores, method_keywords=method_keywords, preserve_existing=True,
        )
        execution = _run_or_resume(native_input, timeout_seconds=timeout_seconds)
        native_input = Path(execution.get("input") or native_input)
        native_output = Path(execution.get("output") or native_input.with_suffix(".out"))
    else:
        if not evidence_path.is_file():
            raise ValueError("Retained atomic stationary evidence is missing")
        if not (Path(folder) / "atom-energy.inp").is_file():
            raise ValueError("Retained atomic method contract is missing")
        # With preserve_existing and a checked existing path this is read-only;
        # it also rejects a different requested method/core/settings contract.
        create_orca_input(
            Path(seed), charge, multiplicity, opt=False, freq=False, label="atom-energy",
            ncores=ncores, method_keywords=method_keywords, preserve_existing=True,
        )
        retained = json.loads(evidence_path.read_text())
        native_input, native_output = Path(retained["energy_input"]), Path(retained["energy_output"])
        _validate_retained_orca_job_input(native_input)
    evidence = read_atomic_hydrogen_energy(
        seed, native_input, native_output, smiles=smiles, charge=charge, multiplicity=multiplicity,
    )
    if evidence_path.is_file():
        if json.loads(evidence_path.read_text()) != evidence:
            raise ValueError("Retained atomic stationary evidence changed; use a new workspace")
    elif execute:
        _write_json(evidence_path, evidence)
    temperature = _frequency_temperature_kwargs(folder).get("thermochemistry_temperature_K")
    if temperature is not None:
        evidence["atomic_thermochemistry"] = atomic_hydrogen_thermochemistry(evidence, temperature)
    return evidence


def _stationary_fragment_records(
    route: RouteSpec,
    workdir: Path,
    mapping: dict,
    *,
    execute: bool,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    retained_endpoints: dict[str, Path] | None = None,
    surface_resolution_policy: SurfaceResolutionPolicy | None = None,
) -> dict:
    """Prepare/validate separated species needed by Arkane."""
    retained_seeds = None
    if retained_endpoints is not None:
        from .route_geometry import retained_endpoint_fragment_geometries
        retained_seeds = retained_endpoint_fragment_geometries(
            route, retained_endpoints["reactants"], retained_endpoints["products"],
        )
        if not retained_seeds.get("valid"):
            return {"valid": False, "status": "retained_endpoint_mapping_unresolved",
                    "retained_seed_evidence": retained_seeds, "reactants": [], "products": []}
        if retained_seeds["mapping"]["surface"] != mapping["surface"]:
            raise ValueError("Retained endpoint states differ from the stationary-point contract")
    surface = mapping["surface"]
    records: dict[str, list[dict]] = {"reactants": [], "products": []}
    for side, smiles_values, declared_labels, charges, multiplicities in (
        ("reactants", route.reactant_smiles, route.reactant_labels, surface["reactant_fragment_charges"],
         surface["reactant_fragment_multiplicities"]),
        ("products", route.product_smiles, route.product_labels, surface["product_fragment_charges"],
         surface["product_fragment_multiplicities"]),
    ):
        labels = list(declared_labels) if len(declared_labels) == len(smiles_values) else list(smiles_values)
        totals = {label: labels.count(label) for label in set(labels)}
        seen: dict[str, int] = {}
        label_signatures: dict[str, list[tuple[str, int, int]]] = {}
        cached: dict[tuple[str, str, int, int], dict] = {}
        for index, (smiles, charge, multiplicity) in enumerate(
            zip(smiles_values, charges, multiplicities), start=1,
        ):
            original_label = str(labels[index - 1])
            seen[original_label] = seen.get(original_label, 0) + 1
            signature = (smiles, int(charge), int(multiplicity))
            distinct = label_signatures.setdefault(original_label, [])
            if signature not in distinct:
                distinct.append(signature)
            # Repeated stoichiometric copies intentionally reuse one species
            # label/job.  A suffix is reserved for a genuine label collision
            # between different structures or electronic states.
            collision_index = distinct.index(signature) + 1
            label = original_label if collision_index == 1 else f"{original_label}__collision_{collision_index}"
            cache_key = (label, *signature)
            if cache_key in cached and retained_seeds is None:
                record = {
                    **cached[cache_key],
                    "stoichiometric_copy": seen[original_label],
                    "reused_stationary_point": True,
                }
                records[side].append(record)
                continue
            folder = Path(workdir) / side / f"{side[:-1]}-{index:02d}"
            folder.mkdir(parents=True, exist_ok=True)
            seed = folder / "input.xyz"
            seed_evidence = None
            if retained_seeds is not None:
                seed_evidence = retained_seeds["fragments"][side][index - 1]
                if seed.is_file():
                    seed_elements, coordinates, _ = read_xyz(seed)
                    if (seed_elements != seed_evidence["elements"] or
                            not np.allclose(coordinates, seed_evidence["coordinates"], atol=1e-9, rtol=0)):
                        raise ValueError("Retained fragment seed differs; use a new calculation workspace")
                else:
                    write_xyz(seed, seed_evidence["elements"], np.asarray(seed_evidence["coordinates"]),
                              "Retained IRC fragment seed; minimum and electronic state unvalidated")
            elif not seed.is_file():
                if mapping.get("source") == "complete_SMILES_nuclear_atom_map_numbers":
                    write_declared_fragment_seed(smiles, seed)
                else:
                    smiles_to_xyz(smiles, seed)
            record = {
                "label": label,
                "original_label": original_label,
                "stoichiometric_copy": seen[original_label],
                "stoichiometric_copy_count": totals[original_label],
                "smiles": smiles,
                "charge": int(charge),
                "multiplicity": int(multiplicity),
                "input_xyz": str(seed),
                "status": "prepared",
                **({"retained_seed_evidence": seed_evidence} if seed_evidence is not None else {}),
            }
            try:
                atom = atomic_hydrogen_identity(smiles, int(charge), int(multiplicity))
            except ValueError as error:
                record.update(status="unsupported_atomic_electronic_state", failure_reason=str(error))
                records[side].append(record)
                continue
            if atom is not None:
                record.update(kind="monatomic_stationary_evidence", **atom)
            if execute:
                print(
                    f"[STORCA] Route {route.route_id}: validating {side[:-1]} species "
                    f"{index}/{len(smiles_values)} ({label})", flush=True,
                )
                try:
                    if atom is not None:
                        record.update(_atomic_stationary_fragment(
                            seed, folder, smiles=smiles, charge=int(charge), multiplicity=int(multiplicity),
                            ncores=ncores, method_keywords=method_keywords,
                            timeout_seconds=timeout_seconds, execute=True,
                        ))
                        record["connectivity_check"] = validate_declared_connectivity(seed, (smiles,))
                        if not record["connectivity_check"].get("valid"):
                            record["status"] = "endpoint_connectivity_validation_failed"
                        print(
                            f"[STORCA] Route {route.route_id}: {side[:-1]} species "
                            f"{index}/{len(smiles_values)} finished ({record['status']})", flush=True,
                        )
                        cached[cache_key] = dict(record)
                        records[side].append(record)
                        continue
                    calculation = _optimize_endpoint(
                        seed, folder, charge=int(charge), multiplicity=int(multiplicity),
                        ncores=ncores, method_keywords=method_keywords,
                        timeout_seconds=timeout_seconds,
                        **({"surface_resolution_policy": surface_resolution_policy}
                           if surface_resolution_policy is not None else {}),
                    )
                    check = calculation["frequency_check"]
                    connectivity = validate_declared_connectivity(
                        Path(calculation["optimized_xyz"]), (smiles,),
                    )
                    if not check.get("IsMinimum"):
                        endpoint_status = "nonminimum_endpoint"
                    elif not connectivity.get("valid"):
                        endpoint_status = "endpoint_connectivity_validation_failed"
                    else:
                        endpoint_status = "validated_minimum"
                    record.update(
                        status=endpoint_status,
                        optimized_xyz=str(calculation["optimized_xyz"]),
                        optimization_output=str(calculation["optimization"]["out"]),
                        frequency_output=str(calculation["frequency"]["out"]),
                        frequency_check=check,
                        connectivity_check=connectivity,
                    )
                except Exception as error:
                    if isinstance(error, (SurfaceResolutionBudgetExhausted, TransferControlError)):
                        raise
                    record.update(status="failed", failure_reason=str(error))
                print(
                    f"[STORCA] Route {route.route_id}: {side[:-1]} species "
                    f"{index}/{len(smiles_values)} finished ({record['status']})", flush=True,
                )
            cached[cache_key] = dict(record)
            records[side].append(record)
    records["valid"] = bool(execute) and all(
        record["status"] in {"validated_minimum", "validated_atom"}
        for side in ("reactants", "products")
        for record in records[side]
    )
    records["status"] = (
        "validated" if records["valid"] else "prepared" if not execute else "stationary_point_validation_failed"
    )
    return records


def reuse_stationary_fragment_records(
    route: RouteSpec, workdir: Path, mapping: dict, *, ncores: int,
    method_keywords: list[str] | None, retained_endpoints: dict[str, Path] | None = None,
) -> dict:
    """Revalidate completed fragment calculations without launching or copying jobs."""
    for side, values in (("reactants", route.reactant_smiles), ("products", route.product_smiles)):
        for index, smiles in enumerate(values, start=1):
            folder = Path(workdir) / side / f"{side[:-1]}-{index:02d}"
            atom = atomic_hydrogen_identity(
                smiles, mapping["surface"][f"{side[:-1]}_fragment_charges"][index - 1],
                mapping["surface"][f"{side[:-1]}_fragment_multiplicities"][index - 1],
            )
            required = (("input.xyz", "atomic-stationary.json") if atom is not None else
                        ("input.xyz", "endpoint-contract.json", "metadata.json", "opt.out", "opt.xyz", "freq.out"))
            for name in required:
                if not (folder / name).is_file():
                    raise ValueError(f"Completed stationary evidence missing: {folder / name}")
    records = _stationary_fragment_records(
        route, workdir, mapping, execute=False, ncores=ncores,
        method_keywords=method_keywords, timeout_seconds=None, retained_endpoints=retained_endpoints,
    )
    if records.get("status") != "prepared":
        return records
    for side in ("reactants", "products"):
        for record in records[side]:
            folder = Path(record["input_xyz"]).parent
            if record.get("kind") == "monatomic_stationary_evidence":
                record.update(_atomic_stationary_fragment(
                    Path(record["input_xyz"]), folder, smiles=record["smiles"], charge=record["charge"],
                    multiplicity=record["multiplicity"], ncores=ncores, method_keywords=method_keywords,
                    timeout_seconds=None, execute=False,
                ), reused_stationary_point=True)
                record["connectivity_check"] = validate_declared_connectivity(
                    Path(record["optimized_xyz"]), (record["smiles"],),
                )
                if not record["connectivity_check"].get("valid"):
                    raise ValueError("Retained atomic stationary nucleus identity changed")
                continue
            contract = json.loads((folder / "endpoint-contract.json").read_text())
            if contract.get("input_xyz_sha256") != hashlib.sha256(Path(record["input_xyz"]).read_bytes()).hexdigest():
                raise ValueError("Retained stationary seed contract was changed")
            calculation = _retained_endpoint_workflow(
                folder, charge=record["charge"], multiplicity=record["multiplicity"],
                ncores=ncores, method_keywords=method_keywords,
            )
            if not calculation or not calculation.get("frequency_check", {}).get("IsMinimum") \
                    or not _orca_optimization_converged(folder / "opt.out"):
                raise ValueError("Retained stationary method/state/completion contract is unresolved")
            check = validate_declared_connectivity(Path(calculation["optimized_xyz"]), (record["smiles"],))
            if not check["valid"]:
                raise ValueError("Retained stationary geometry changed chemical identity")
            record.update(status="validated_minimum", optimized_xyz=str(calculation["optimized_xyz"]),
                          optimization_output=str(folder / "opt.out"), frequency_output=str(folder / "freq.out"),
                          frequency_check=calculation["frequency_check"], connectivity_check=check,
                          reused_stationary_point=True,
                          retained_evidence_hashes={name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
                              for name in ("input.xyz", "endpoint-contract.json", "metadata.json", "opt.inp", "opt.out", "opt.xyz", "freq.inp", "freq.out", "freq.hess", "opt.process.json", "freq.process.json")
                              if (folder / name).is_file()})
    records.update(valid=True, status="validated", reuse_mode="read_only_no_electronic_jobs")
    return records


def _route_uses_separated_fragment_endpoints(route: RouteSpec) -> bool:
    """Whether either reaction side is an explicitly separated channel."""
    return len(route.reactant_smiles) > 1 or len(route.product_smiles) > 1


def _stationary_fragment_xyz(stationary_points: dict, *, optimized: bool) -> dict[str, list[str]]:
    """Extract one atom-ordered geometry for every stoichiometric fragment."""
    key = "optimized_xyz" if optimized else "input_xyz"
    paths: dict[str, list[str]] = {"reactants": [], "products": []}
    for side in paths:
        for record in stationary_points.get(side) or []:
            value = record.get(key)
            if not value or not Path(value).is_file():
                raise FileNotFoundError(
                    f"Stationary-point record for {record.get('label') or side} has no retained {key}"
                )
            paths[side].append(str(value))
    return paths


def _validate_irc_route_endpoints(
    route: RouteSpec,
    trajectory: Path,
    reactant_xyz: Path,
    product_xyz: Path,
    *,
    mapping_result: dict,
    separated_fragments: bool,
    rmsd_tolerance_angstrom: float = 0.75,
) -> dict:
    """Keep mapped chemistry mandatory; reference RMSD diagnoses conformers.

    A failed reference RMSD is rescued only by independent, exact supported
    bond/state assignment, never by the distance connectivity screen alone.
    Reached ends still need local optimization and frequencies before rates.
    """
    coordinate_check = None if separated_fragments else validate_trajectory_endpoints(
        trajectory, reactant_xyz, product_xyz, rmsd_tolerance_angstrom=rmsd_tolerance_angstrom,
    )
    try:
        frames = xyz_frames(trajectory)
        expected_elements, _, _ = read_xyz(reactant_xyz)
    except (OSError, ValueError) as error:
        return {"status": "irc_endpoint_frames_unreadable", "valid": False, "reason": str(error)}
    integrity = validate_trajectory_frames(frames, expected_elements=expected_elements)
    if not integrity["valid"]:
        return {**integrity, "status": "irc_endpoint_frames_incomplete"}
    first = Path(trajectory).with_name("irc-endpoint-first.xyz")
    last = Path(trajectory).with_name("irc-endpoint-last.xyz")
    write_xyz(first, frames[0][0], frames[0][1], "IRC endpoint; first retained frame")
    write_xyz(last, frames[-1][0], frames[-1][1], "IRC endpoint; last retained frame")
    forward = {
        "reactant": validate_route_endpoint_connectivity(
            route, first, side="reactants", mapping_result=mapping_result,
        ),
        "product": validate_route_endpoint_connectivity(
            route, last, side="products", mapping_result=mapping_result,
        ),
    }
    reverse = {
        "reactant": validate_route_endpoint_connectivity(
            route, last, side="reactants", mapping_result=mapping_result,
        ),
        "product": validate_route_endpoint_connectivity(
            route, first, side="products", mapping_result=mapping_result,
        ),
    }
    forward_valid = forward["reactant"].get("valid") and forward["product"].get("valid")
    reverse_valid = reverse["reactant"].get("valid") and reverse["product"].get("valid")
    conformer_checks = {}
    rescued = []
    if coordinate_check is not None:
        directions = coordinate_check.get("direction_checks") or {}
        original_forward = forward_valid and bool((directions.get("forward") or {}).get("valid"))
        original_reverse = reverse_valid and bool((directions.get("reverse") or {}).get("valid"))
        from .endpoint_identity import validate_endpoint_identity
        product_map = {int(product): int(reactant)
                       for reactant, product in mapping_result["all_atom_mapping"]}
        for direction, chemistry_valid, endpoints in (
            ("forward", forward_valid, (first, last)),
            ("reverse", reverse_valid, (last, first)),
        ):
            if chemistry_valid and not (directions.get(direction) or {}).get("valid"):
                checks = {
                    "reactant": validate_endpoint_identity(endpoints[0], route.reactant_smiles),
                    "product": validate_endpoint_identity(endpoints[1], route.product_smiles,
                                                          expected_to_xyz=product_map),
                }
                conformer_checks[direction] = checks
                if all(check.get("valid") for check in checks.values()):
                    rescued.append(direction)
        # Never orient a path arbitrarily when chemistry cannot distinguish
        # both directions and neither matches the original mapped reference.
        unique_rescue = rescued[0] if len(rescued) == 1 else None
        forward_valid = original_forward or unique_rescue == "forward"
        reverse_valid = original_reverse or unique_rescue == "reverse"
    return {
        "schema_version": 1,
        "kind": "separated_fragment_irc_endpoint_check" if separated_fragments else "bound_fragment_irc_endpoint_check",
        "status": ("validated" if forward_valid or reverse_valid else
                   "irc_endpoint_identity_unresolved" if conformer_checks else "declared_fragment_channels_not_reached"),
        "valid": bool(forward_valid or reverse_valid),
        "matched_direction": "forward" if forward_valid else "reverse" if reverse_valid else None,
        "first_endpoint_xyz": str(first),
        "last_endpoint_xyz": str(last),
        "forward_checks": forward,
        "reverse_checks": reverse,
        "coordinate_correspondence_required": False,
        "coordinate_validation": coordinate_check,
        "conformer_identity_checks": conformer_checks,
        "conformer_rescue_orientation_ambiguous": len(rescued) > 1,
        "reference_conformer_mismatch": bool(rescued and (forward_valid or reverse_valid)),
        "minimum_refinement_required": bool(forward_valid or reverse_valid),
        "stationary_point_claim_allowed": False,
        "trajectory_integrity": integrity,
        "interpretation": (
            "Both IRC ends must match mapped declared chemistry and stereochemistry. A bound reference "
            "RMSD mismatch is a conformer diagnostic, not proof of a different molecule; rescuing it "
            "requires supported independent bond/state assignment. Reached ends are not certified "
            "minima until separately optimized and frequency-checked. Electronic method applicability remains separate."
        ),
    }


def _refine_reached_irc_stationary_points(
    route: RouteSpec, trial: dict, workdir: Path, mapping: dict, *,
    ncores: int, method_keywords: list[str] | None, timeout_seconds: float | None,
    reference_endpoints: dict[str, str] | None = None,
    surface_resolution_policy: SurfaceResolutionPolicy | None = None,
) -> dict:
    """Bind rate artifacts to locally validated minima reached by THIS IRC.

    Molecular identity is not conformer thermodynamics. A distinct refined
    conformer may support a conditional elementary rate, but it cannot silently
    replace a graph-species rate in a thermochemical model with another basis.
    """
    endpoint_check = (trial.get("irc") or {}).get("endpoint_validation") or {}
    policy = surface_resolution_policy or SurfaceResolutionPolicy()
    if not endpoint_check.get("valid") or endpoint_check.get("matched_direction") not in {"forward", "reverse"}:
        return {"stationary_points": {"valid": False, "status": "reached_endpoint_evidence_missing",
                                      "reactants": [], "products": []}}
    first, last = endpoint_check["first_endpoint_xyz"], endpoint_check["last_endpoint_xyz"]
    reached = {"reactants": Path(first), "products": Path(last)}
    if endpoint_check["matched_direction"] == "reverse":
        reached = {"reactants": Path(last), "products": Path(first)}
    print(f"[STORCA] Route {route.route_id}: refining actual IRC species (molecules: Opt + Freq; atoms: energy/state)", flush=True)
    stationary = _stationary_fragment_records(
        route, workdir, mapping, execute=True, ncores=ncores, method_keywords=method_keywords,
        timeout_seconds=timeout_seconds, retained_endpoints=reached,
        surface_resolution_policy=policy,
    )
    stationary["geometry_origin"] = "locally_refined_reached_irc_endpoints"
    stationary["irc_endpoint_validation"] = endpoint_check
    value = {"stationary_points": stationary}
    if not stationary.get("valid"):
        return value
    from .endpoint_identity import validate_endpoint_identity
    rescued = endpoint_check.get("reference_conformer_mismatch") is True
    for side in ("reactants", "products"):
        for record in stationary[side]:
            identity = validate_endpoint_identity(Path(record["optimized_xyz"]), (record["smiles"],))
            record["independent_identity_check"] = identity
            if ((rescued and not identity.get("valid"))
                    or identity.get("status") == "mapped_chemical_identity_mismatch"):
                record["status"] = "refined_endpoint_identity_unresolved"
    if any(record["status"] not in {"validated_minimum", "validated_atom"}
           for side in ("reactants", "products") for record in stationary[side]):
        stationary.update(valid=False, status="stationary_point_validation_failed")
        return value

    references = reference_endpoints or trial.get("endpoint_reference_geometries")
    if not references:
        endpoint_validation = trial.get("endpoint_validation") or {}
        references = {side: (endpoint_validation.get(name) or {}).get("endpoint_xyz")
                      for side, name in (("reactants", "reactant"), ("products", "product"))}
    reference_seeds = None
    if references and all(references.get(side) for side in ("reactants", "products")):
        from .route_geometry import retained_endpoint_fragment_geometries
        reference_seeds = retained_endpoint_fragment_geometries(
            route, Path(references["reactants"]), Path(references["products"]),
        )
    comparisons = {"reactants": [], "products": []}
    comparable = bool(reference_seeds and reference_seeds.get("valid"))
    for side in ("reactants", "products"):
        for index, record in enumerate(stationary[side]):
            comparison = {"label": record["label"], "status": "reference_conformer_unresolved"}
            if comparable:
                reference = reference_seeds["fragments"][side][index]
                elements, coordinates, _ = read_xyz(Path(record["optimized_xyz"]))
                if elements == reference["elements"]:
                    rmsd = aligned_rmsd(coordinates, np.asarray(reference["coordinates"]))
                    comparison.update(
                        status="reference_geometry_correspondence" if rmsd <= policy.endpoint_rmsd_tolerance_angstrom else "distinct_reference_conformer",
                        mapped_rmsd_angstrom=rmsd, diagnostic_tolerance_angstrom=policy.endpoint_rmsd_tolerance_angstrom,
                        reference_xyz=reference["source_xyz"],
                        reference_xyz_sha256=reference["source_xyz_sha256"],
                    )
            comparisons[side].append(comparison)
    same_basis = all(item["status"] == "reference_geometry_correspondence"
                     for side in comparisons for item in comparisons[side])
    basis = {
        "schema_version": 1, "kind": "reached_endpoint_conformer_rate_basis",
        "status": "reference_geometry_correspondence" if same_basis else "conformer_thermo_basis_unresolved",
        "network_reinjection_allowed": same_basis,
        "comparisons": comparisons,
        "fast_conformer_equilibration_assumed": False,
        "interpretation": (
            "The rate uses locally refined IRC minima. Reference geometry correspondence is a conservative "
            "handoff diagnostic, not proof of conformer completeness or absolute RMG thermochemistry accuracy. "
            "A distinct or unresolved conformer basis requires explicit thermochemical/equilibration evidence "
            "before graph-species reinjection; a conditional elementary rate may still be calculated."
        ),
    }
    def pinned(path):
        path = Path(path).resolve(strict=True)
        return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    artifact_contract = {"schema_version": 1, "kind": "reached_endpoint_artifact_contract"}
    for side in ("reactants", "products"):
        artifact_contract[side] = [({"kind": "monatomic_stationary_evidence",
                                    "optimized_xyz": pinned(record["optimized_xyz"]),
                                    "energy_input": pinned(record["energy_input"]),
                                    "energy_output": pinned(record["energy_output"]),
                                    "mass_amu": record["mass_amu"], "isotope": record["isotope"]}
                                   if record.get("kind") == "monatomic_stationary_evidence" else
                                  {"optimized_xyz": pinned(record["optimized_xyz"]),
                                   "optimization_output": pinned(record["optimization_output"]),
                                   "frequency_output": pinned(record["frequency_output"]),
                                   "frequency_hessian": pinned(Path(record["frequency_output"]).with_suffix(".hess"))})
                                  for record in stationary[side]]
    artifact_contract["transition_state"] = {
        "frequency_output": pinned(trial["transition_state"]["frequency_output"]),
        "frequency_hessian": pinned(Path(trial["transition_state"]["frequency_output"]).with_suffix(".hess")),
    }
    value.update(conformer_rate_basis=basis, reached_endpoint_artifact_contract=artifact_contract)
    _write_json(Path(workdir) / "reached-endpoint-rate-handoff.json", value)
    return value


def _run_generic_orientation(
    route: RouteSpec,
    seed_record: dict,
    *,
    mapping_result: dict,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    neb_images: int,
    barrierless_classifier: Callable[[dict], dict] | None,
    surface_resolution_policy: SurfaceResolutionPolicy | None = None,
    surface_bond_changes: dict | None = None,
    separated_fragments: bool = False,
    stationary_points: dict | None = None,
) -> dict:
    """Run one endpoint orientation through NEB-TS, TS/frequency, and IRC."""
    surface_resolution_policy = surface_resolution_policy or SurfaceResolutionPolicy()
    orientation = int(seed_record["orientation"])
    folder = Path(seed_record["reactant_xyz"]).parent
    record: dict[str, Any] = {"orientation": orientation, "status": "running"}
    try:
        if separated_fragments:
            if not stationary_points or not stationary_points.get("valid"):
                raise RuntimeError("Separated-fragment NEB requires validated stationary species")
            reactant_source = Path(seed_record["reactant_xyz"])
            product_source = Path(seed_record["product_xyz"])
            reactant_check = {
                "status": "not_applicable_to_separated_channel",
                "IsMinimum": None,
                "reason": "The encounter geometry is assembled from independently validated fragment minima.",
            }
            product_check = dict(reactant_check)
        else:
            reactant_result = _optimize_endpoint(
                Path(seed_record["reactant_xyz"]), folder / "reactant-complex",
                charge=route.charge, multiplicity=route.multiplicity, ncores=ncores,
                method_keywords=method_keywords, timeout_seconds=timeout_seconds,
            )
            product_result = _optimize_endpoint(
                Path(seed_record["product_xyz"]), folder / "product-complex",
                charge=route.charge, multiplicity=route.multiplicity, ncores=ncores,
                method_keywords=method_keywords, timeout_seconds=timeout_seconds,
            )
            reactant_source = Path(reactant_result["optimized_xyz"])
            product_source = Path(product_result["optimized_xyz"])
            reactant_check = reactant_result["frequency_check"]
            product_check = product_result["frequency_check"]
        reactant_connectivity = validate_route_endpoint_connectivity(
            route, reactant_source,
            side="reactants", mapping_result=mapping_result,
        )
        product_connectivity = validate_route_endpoint_connectivity(
            route, product_source,
            side="products", mapping_result=mapping_result,
        )
        endpoint_validation = {
            "mode": (
                "assembled_from_validated_separated_fragment_minima"
                if separated_fragments else "optimized_bound_endpoint_minima"
            ),
            "full_complex_minimum_required": not separated_fragments,
            "stationary_fragment_validation_status": (
                stationary_points.get("status") if stationary_points else "not_applicable"
            ),
            "reactant": {
                "endpoint_xyz": str(reactant_source),
                "frequency_check": reactant_check,
                "connectivity_check": reactant_connectivity,
            },
            "product": {
                "endpoint_xyz": str(product_source),
                "frequency_check": product_check,
                "connectivity_check": product_connectivity,
            },
            "valid": bool(
                (separated_fragments or (
                    reactant_check.get("IsMinimum") and product_check.get("IsMinimum")
                ))
                and reactant_connectivity.get("valid") and product_connectivity.get("valid")
            ),
        }
        if not separated_fragments:
            endpoint_validation["reactant"].update(
                optimization_output=str(reactant_result["optimization"]["out"]),
                frequency_output=str(reactant_result["frequency"]["out"]),
            )
            endpoint_validation["product"].update(
                optimization_output=str(product_result["optimization"]["out"]),
                frequency_output=str(product_result["frequency"]["out"]),
            )
        record["endpoint_validation"] = endpoint_validation
        endpoint_sources = tuple(str(value) for item in (
            endpoint_validation["reactant"], endpoint_validation["product"])
            for key, value in item.items() if key in {"optimization_output", "frequency_output"})
        reference = EndpointReferenceAssessment(
            "finite_encounter" if separated_fragments else
            "bound_minimum" if endpoint_validation["valid"] else "unresolved",
            sources=endpoint_sources,
            reason=("Independently validated fragments at finite separation are not an asymptote."
                    if separated_fragments else "Native optimized/frequency-validated bound endpoints."),
        ).as_dict()
        record["endpoint_reference"] = reference
        if not endpoint_validation["valid"]:
            connectivity_failed = not (
                reactant_connectivity.get("valid") and product_connectivity.get("valid")
            )
            record.update(
                status=(
                    "endpoint_connectivity_validation_failed"
                    if connectivity_failed else "endpoint_validation_failed"
                ),
                path_classification="surface_unresolved",
            )
            return record
        path_dir = folder / "path"
        path_dir.mkdir(parents=True, exist_ok=True)
        reactant_xyz = path_dir / "reactant.xyz"
        product_xyz = path_dir / "product.xyz"
        shutil.copy2(reactant_source, reactant_xyz)
        shutil.copy2(product_source, product_xyz)
        # These endpoints have already been independently validated as
        # minima.  Asking ORCA to preoptimize them once more makes the NEB
        # robust to small residual forces.  Encounter geometries deliberately
        # combine separated fragments at a selected orientation, so preopting
        # those ends would change the reaction channel and is prohibited.
        preopt_ends = not separated_fragments
        neb_input = create_orca_neb_ts_input(
            reactant_xyz, product_xyz, charge=route.charge, multiplicity=route.multiplicity,
            label="neb-ts", ncores=ncores, nimages=neb_images, method_keywords=method_keywords,
            preopt_ends=preopt_ends,
        )
        record["neb_endpoint_policy"] = {
            "preopt_ends": preopt_ends,
            "free_end_neb": False,
            "reason": (
                "both declared endpoints are bound, validated minima"
                if preopt_ends else
                "assembled separated-fragment encounter endpoints must retain their declared channel"
            ),
        }
        try:
            neb_execution = _run_or_resume(neb_input, timeout_seconds=timeout_seconds)
        except Exception as error:
            if isinstance(error, TransferControlError):
                raise
            failed_execution = getattr(error, "storca_failed_execution", {})
            neb_execution = {
                "status": "failed", "input": failed_execution.get("input", str(neb_input)),
                "output": failed_execution.get("output", str(neb_input.with_suffix(".out"))),
                "failure_reason": str(error), "failed_execution": failed_execution or None,
            }
        record["neb"] = neb_execution
        semantic_no_barrier = neb_execution["status"] == "completed_no_interior_barrier"
        execution_stem = Path(neb_execution.get("input") or neb_input).stem
        ts_xyz = _locate_neb_ts(path_dir, execution_stem) if neb_execution["status"] == "completed" else None
        if ts_xyz is None:
            path_quality_warning = neb_execution.get("path_quality_warning") or {}
            angle_candidates = list(path_quality_warning.get("image_indices") or [])
            explicit_candidates = _neb_intermediate_image_indices(Path(neb_execution["output"]))
            # For a cusp warning the latter image is the turning-point side
            # of the reported pair.  It is only a candidate until it passes
            # the same optimization/frequency gate as an explicit ORCA
            # intermediate.
            candidate_indices = (
                None if explicit_candidates else ([max(angle_candidates)] if angle_candidates else None)
            )
            intermediate_evidence = _validate_neb_intermediates(
                route, neb_output=Path(neb_execution["output"]), path_dir=path_dir,
                ncores=ncores, method_keywords=method_keywords,
                timeout_seconds=timeout_seconds,
                candidate_indices=candidate_indices,
                detection_source=(
                    "maximum_angle_path_cusp_warning" if candidate_indices
                    else "explicit_intermediate_warning"
                ),
                surface_resolution_policy=surface_resolution_policy,
            )
            if intermediate_evidence.get("status") == "intermediate_warning_dismissed_endpoint_return":
                ignored_indices = sorted({
                    int(item["image_index"])
                    for item in intermediate_evidence.get("candidates") or []
                    if (item.get("intermediate_identity") or {}).get("status") == "endpoint_return"
                    and item.get("image_index") is not None
                })
                if not ignored_indices:
                    raise RuntimeError(
                        "Endpoint-return recovery has no validated warning image identity"
                    )
                recovery = _run_or_resume(
                    neb_input,
                    timeout_seconds=timeout_seconds,
                    ignore_intermediate_warning_indices=ignored_indices,
                )
                recovery["validated_endpoint_return_identity"] = {
                    "image_indices": ignored_indices,
                    "identity_evidence": [
                        item.get("intermediate_identity")
                        for item in intermediate_evidence.get("candidates") or []
                        if (item.get("intermediate_identity") or {}).get("status") == "endpoint_return"
                    ],
                }
                record["endpoint_return_neb_recovery"] = recovery
                record["neb"] = recovery
                neb_execution = recovery
                semantic_no_barrier = neb_execution["status"] == "completed_no_interior_barrier"
                execution_stem = Path(neb_execution.get("input") or neb_input).stem
                ts_xyz = (
                    _locate_neb_ts(path_dir, execution_stem)
                    if neb_execution["status"] == "completed" else None
                )
                if ts_xyz is not None:
                    ts_result = _verify_ts_seed(
                        route,
                        ts_xyz,
                        folder / "transition-state",
                        reactant_xyz=reactant_xyz,
                        product_xyz=product_xyz,
                        mapping_result=mapping_result,
                        separated_fragments=separated_fragments,
                        ncores=ncores,
                        method_keywords=method_keywords,
                        timeout_seconds=timeout_seconds,
                        surface_resolution_policy=surface_resolution_policy,
                    )
                    ts_result["intermediate_warning_evidence"] = intermediate_evidence
                    record.update(ts_result)
                    return record
                path_quality_warning = neb_execution.get("path_quality_warning") or {}
                explicit_candidates = _neb_intermediate_image_indices(
                    Path(neb_execution["output"])
                )
                new_candidate_indices = [
                    index for index in explicit_candidates if index not in ignored_indices
                ]
                if new_candidate_indices:
                    intermediate_evidence = _validate_neb_intermediates(
                        route,
                        neb_output=Path(neb_execution["output"]),
                        path_dir=path_dir,
                        ncores=ncores,
                        method_keywords=method_keywords,
                        timeout_seconds=timeout_seconds,
                        candidate_indices=new_candidate_indices,
                        detection_source="post_endpoint_return_restart_intermediate_warning",
                        surface_resolution_policy=surface_resolution_policy,
                    )
            if intermediate_evidence.get("status") == "validated_intermediate_detected":
                parent_endpoint_stationary = None
                if not separated_fragments:
                    def endpoint_record(
                        side: str, smiles: str, label: str, charge: int,
                        multiplicity: int, calculation: dict,
                    ) -> dict:
                        return {
                            "label": label,
                            "original_label": label,
                            "smiles": smiles,
                            "charge": charge,
                            "multiplicity": multiplicity,
                            "status": "validated_minimum",
                            "optimized_xyz": str(calculation["optimized_xyz"]),
                            "optimization_output": str(calculation["optimization"]["out"]),
                            "frequency_output": str(calculation["frequency"]["out"]),
                            "source": f"parent_route_{side}_endpoint",
                        }

                    parent_endpoint_stationary = {
                        "reactants": [endpoint_record(
                            "reactant",
                            route.reactant_smiles[0],
                            (route.reactant_labels or ("reactant",))[0],
                            (route.reactant_charges or (route.charge,))[0],
                            (route.reactant_multiplicities or (route.multiplicity,))[0],
                            reactant_result,
                        )],
                        "products": [endpoint_record(
                            "product",
                            route.product_smiles[0],
                            (route.product_labels or ("product",))[0],
                            (route.product_charges or (route.charge,))[0],
                            (route.product_multiplicities or (route.multiplicity,))[0],
                            product_result,
                        )],
                    }
                intermediate_evidence["segmented_verification"] = _run_segmented_intermediate_nebs(
                    route,
                    reactant_xyz=reactant_xyz,
                    product_xyz=product_xyz,
                    intermediate_evidence=intermediate_evidence,
                    path_dir=path_dir,
                    separated_fragments=separated_fragments,
                    ncores=ncores,
                    method_keywords=method_keywords,
                    timeout_seconds=timeout_seconds,
                    neb_images=neb_images,
                    surface_resolution_policy=surface_resolution_policy,
                    parent_endpoint_stationary=parent_endpoint_stationary,
                )
                segmented_status = str(
                    (intermediate_evidence.get("segmented_verification") or {}).get("status")
                    or "segment_path_unresolved"
                )
                route_status = {
                    "verified_segments_rate_generation_required": "verified_multistep_path",
                    "segment_rate_model_unavailable": "segment_rate_model_unavailable",
                    "barrierless_rate_model_required": "barrierless_rate_model_required",
                    "surface_resolution_budget_exhausted": "surface_resolution_budget_exhausted",
                    "intermediate_identity_unresolved": "intermediate_identity_unresolved",
                }.get(segmented_status, "segment_path_unresolved")
                record.update(
                    status=route_status,
                    path_classification="multistep_intermediate_detected",
                    classification_evidence={
                        "status": route_status,
                        "reason": (
                            "validated_intermediate_segment_status_" + segmented_status
                        ),
                    },
                    intermediate_evidence=intermediate_evidence,
                    path_evidence={
                        "kind": "orca_neb_intermediate_path",
                        "trajectory": str(_locate_neb_trajectory(path_dir, execution_stem) or "") or None,
                        "rate_claim_allowed": False,
                    },
                    trajectory=None,
                )
                if route_status == "verified_multistep_path":
                    record["path_evidence"] = {
                        "kind": "verified_segmented_ts_frequency_irc_paths",
                        "segment_count": len(
                            intermediate_evidence["segmented_verification"].get("segments") or []
                        ),
                        "rate_claim_allowed": False,
                    }
                return record
            if intermediate_evidence.get("status") == "intermediate_identity_unresolved":
                record.update(
                    status="intermediate_identity_unresolved",
                    path_classification="surface_unresolved",
                    classification_evidence={
                        "status": "intermediate_identity_unresolved",
                        "reason": "validated_stationary_minimum_has_unresolved_chemical_identity",
                    },
                    intermediate_evidence=intermediate_evidence,
                    path_evidence={
                        "kind": "orca_neb_intermediate_path",
                        "trajectory": str(_locate_neb_trajectory(path_dir, execution_stem) or "") or None,
                        "rate_claim_allowed": False,
                    },
                    trajectory=None,
                )
                return record
            if intermediate_evidence.get("status") != "not_detected":
                record["intermediate_evidence"] = intermediate_evidence
            if path_quality_warning:
                record["path_quality_warning"] = path_quality_warning
            trajectory = _locate_neb_trajectory(path_dir, execution_stem)
            energies = _trajectory_energies(trajectory) if trajectory else None
            frames = xyz_frames(trajectory) if trajectory else []
            frame_count = len(frames)
            trajectory_validation = (
                validate_trajectory_endpoints(trajectory, reactant_xyz, product_xyz)
                if trajectory else None
            )
            diagnostic_evidence_complete = bool(
                trajectory is not None
                and frame_count >= neb_images
                and energies is not None
                and len(energies) == frame_count
                and trajectory_validation
                and trajectory_validation["valid"]
            )
            surface_diagnosis = diagnose_surface(
                energies_hartree=energies,
                coordinates=[coordinates for _, coordinates, _ in frames] if frames else None,
                atom_labels=frames[0][0] if frames else None,
                bond_changes=surface_bond_changes or {
                    "broken_bonds": [list(item) for item in route.broken_bonds],
                    "formed_bonds": [list(item) for item in route.formed_bonds],
                },
                fragment_atom_indices=_route_fragment_atom_indices(route, mapping_result),
                explicit_intermediate_indices=explicit_candidates,
                output_text=(
                    Path(neb_execution["output"]).read_text(errors="replace")
                    if Path(neb_execution["output"]).is_file() else None
                ),
                evidence_complete=diagnostic_evidence_complete,
                endpoints_valid=bool(trajectory_validation and trajectory_validation.get("valid")),
                electronic_state_status="applicable_single_surface",
                reactant_fragment_count=len(route.reactant_smiles),
                product_fragment_count=len(route.product_smiles),
                policy=surface_resolution_policy,
            )
            record["surface_diagnosis"] = surface_diagnosis.as_dict()
            if path_quality_warning.get("kind") == (
                "submerged_climbing_image_below_endpoint_asymptote"
            ) and path_quality_warning.get("valid"):
                if not reference["energetic_classification_permitted"]:
                    record.update(
                        status="endpoint_reference_unqualified",
                        path_classification="surface_unresolved",
                        classification_evidence={
                            **path_quality_warning, "endpoint_reference": reference,
                            "interpretation": "A finite encounter is not a separated-fragment energy reference.",
                            "mechanism_claim_allowed": False, "rate_claim_allowed": False,
                        },
                        path_evidence={"output": neb_execution.get("output"),
                                       "trajectory": str(trajectory) if trajectory else None},
                    )
                    return record
                record.update(
                    status="dissociation_rate_model_required",
                    path_classification="submerged_saddle_dissociation_candidate",
                    classification_evidence={
                        **path_quality_warning,
                        "status": "dissociation_rate_model_required",
                        "rate_claim_allowed": False,
                        "lifetime_claim_allowed": False,
                    },
                    method_selection={
                        "status": "specialized_rate_model_required",
                        "detection": "automatic_retained_energy_diagnostic",
                        "ordinary_harmonic_tst_applicable": False,
                        "candidate_methods": [
                            "variational_transition_state_theory",
                            "capture_theory_for_reverse_association",
                            "rrkm_master_equation_for_pressure_dependent_dissociation",
                        ],
                        "required_evidence": [
                            "validated_fragment_thermochemistry",
                            "long_range_or_variational_dividing_surface",
                            "pressure_and_bath_gas_contract_when_falloff_is_relevant",
                        ],
                        "selection_reason": (
                            "the retained climbing image is below the higher endpoint "
                            "asymptote and cannot control forward passage to that endpoint"
                        ),
                        "rate_claim_allowed": False,
                    },
                    path_evidence={
                        "kind": "orca_submerged_climbing_image_dissociation_diagnostic",
                        "output": neb_execution.get("output"),
                        "trajectory": str(trajectory) if trajectory else None,
                        "energies_hartree": energies,
                        "endpoint_validation": trajectory_validation,
                        "rate_claim_allowed": False,
                    },
                    trajectory=(
                        str(trajectory)
                        if trajectory and trajectory_validation
                        and trajectory_validation.get("valid") else None
                    ),
                )
                return record
            evidence_complete = bool(
                neb_execution["status"] in {"completed", "completed_no_interior_barrier"}
                and diagnostic_evidence_complete
            )
            if not evidence_complete:
                record.update(
                    status="surface_unresolved",
                    path_classification="surface_unresolved",
                    classification_evidence={
                        "status": "surface_unresolved",
                        "reason": "incomplete_neb_path_evidence",
                        "normal_termination": neb_execution["status"] == "completed",
                        "semantic_no_barrier_completion": semantic_no_barrier,
                        "orca_neb_outcome": neb_execution.get("orca_neb_outcome"),
                        "expected_minimum_image_count": neb_images,
                        "retained_image_count": frame_count,
                        "energies_complete": energies is not None and len(energies) == frame_count,
                        "endpoint_validation": trajectory_validation,
                    },
                    path_evidence={
                        "kind": "incomplete_neb_path",
                        "trajectory": str(trajectory) if trajectory else None,
                        "energies_hartree": energies,
                        "rate_claim_allowed": False,
                    },
                    trajectory=None,
                )
                return record
            context = {
                "route": route.as_dict(), "orientation": orientation,
                "trajectory": str(trajectory) if trajectory else None,
                "energies_hartree": energies,
                "orca_neb_outcome": neb_execution.get("orca_neb_outcome"),
                "default_classification": classify_path_energy_profile(
                    energies,
                    reactant_fragment_count=len(route.reactant_smiles),
                    product_fragment_count=len(route.product_smiles),
                ),
            }
            classification = barrierless_classifier(context) if barrierless_classifier else context["default_classification"]
            allowed = {
                "barriered_ts_candidate", "barrierless_capture_candidate",
                "barrierless_dissociation_candidate", "barrierless_path_candidate", "surface_unresolved",
            }
            if not isinstance(classification, dict) or classification.get("status") not in allowed:
                raise ValueError("Barrierless-classification hook returned an unsupported status")
            record.update(
                status=("classified_without_transition_state" if classification["status"] != "surface_unresolved"
                        else "surface_unresolved"),
                path_classification=classification["status"],
                classification_evidence=classification,
                path_evidence={
                    "kind": (
                        "orca_converged_neb_without_interior_barrier"
                        if semantic_no_barrier else "neb_energy_path_without_validated_ts"
                    ),
                    "trajectory": str(trajectory) if trajectory else None,
                    "energies_hartree": energies,
                    "endpoint_validation": trajectory_validation,
                    "orca_neb_outcome": neb_execution.get("orca_neb_outcome"),
                    "rate_claim_allowed": False,
                },
                trajectory=str(trajectory) if trajectory and trajectory_validation and trajectory_validation["valid"] else None,
            )
            return record

        ts_result = _verify_ts_seed(
            route,
            ts_xyz,
            folder / "transition-state",
            reactant_xyz=reactant_xyz,
            product_xyz=product_xyz,
            mapping_result=mapping_result,
            separated_fragments=separated_fragments,
            ncores=ncores,
            method_keywords=method_keywords,
            timeout_seconds=timeout_seconds,
            surface_resolution_policy=surface_resolution_policy,
        )
        record.update(ts_result)
        return record
    except Exception as error:
        if isinstance(error, TransferControlError):
            raise
        if isinstance(error, SurfaceResolutionBudgetExhausted):
            record.update(
                status="surface_resolution_budget_exhausted",
                path_classification="surface_unresolved",
                failure_reason=str(error),
                surface_resolution_budget=error.evidence,
            )
        else:
            record.update(
                status="path_execution_failed", path_classification="surface_unresolved",
                failure_reason=str(error),
            )
        return record


def _surface_diagnoses_from_records(records: list[dict]) -> list[SurfaceDiagnosis]:
    diagnoses: list[SurfaceDiagnosis] = []
    for item in records:
        artifact = item.get("surface_diagnosis")
        if not isinstance(artifact, dict) or artifact.get("kind") != "storca_surface_diagnosis":
            continue
        try:
            diagnoses.append(SurfaceDiagnosis(
                status=str(artifact["status"]),
                reasons=tuple(str(value) for value in artifact.get("reasons") or ()),
                evidence=artifact.get("evidence") or {},
                corroborated=bool(artifact.get("corroborated")),
                recommended_strategy=artifact.get("recommended_strategy"),
                claims_permitted=tuple(artifact.get("claims_permitted") or ()),
                claims_prohibited=tuple(artifact.get("claims_prohibited") or ()),
            ))
        except (KeyError, TypeError, ValueError):
            continue
    return diagnoses


def _route_fragment_atom_indices(route: RouteSpec, mapping: dict) -> list[list[int]]:
    """Return a separated-side fragment partition in reactant XYZ atom order."""
    from rdkit import Chem

    def counts(values: tuple[str, ...]) -> list[int]:
        result = []
        for value in values:
            molecule = Chem.AddHs(Chem.MolFromSmiles(value))
            if molecule is None:
                raise ValueError(f"invalid fragment SMILES: {value}")
            result.append(molecule.GetNumAtoms())
        return result

    if len(route.reactant_smiles) >= 2:
        groups = []
        offset = 0
        for count in counts(route.reactant_smiles):
            groups.append(list(range(offset, offset + count)))
            offset += count
        return groups
    if len(route.product_smiles) >= 2:
        product_groups = []
        offset = 0
        for count in counts(route.product_smiles):
            product_groups.append(list(range(offset, offset + count)))
            offset += count
        product_to_reactant = {
            int(product): int(reactant)
            for reactant, product in mapping.get("all_atom_mapping") or ()
        }
        if len(product_to_reactant) != offset:
            return []
        return [sorted(product_to_reactant[index] for index in group) for group in product_groups]
    return []


def _surface_environment_applicability(
    route: RouteSpec,
    condition_contract: dict | None,
    policy: SurfaceResolutionPolicy,
    method_keywords: list[str] | None,
) -> dict[str, Any]:
    environment = dict(route.environment or {})
    raw = (
        environment.get("model")
        or environment.get("phase_model")
        or environment.get("phase")
        or (condition_contract or {}).get("phase_model")
    )
    source = "route_environment" if environment else (
        "condition_contract" if condition_contract else "legacy_generic_gas_default"
    )
    if raw is None:
        canonical = "homogeneous_gas"
    else:
        normalized = str(raw).strip().lower().replace("-", "_")
        if "gas" in normalized:
            canonical = "homogeneous_gas"
        elif (
            normalized in {"solution", "implicit_solution", "implicit_solvent", "pcm", "cpcm", "smd"}
            or ("implicit" in normalized and ("solution" in normalized or "solvent" in normalized))
        ):
            canonical = "implicit_solution"
        else:
            canonical = normalized
    keywords = list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"])
    solvent_tokens = re.findall(
        r"\b(CPCM|SMD)\s*\(([^\n)]*)\)", " ".join(keywords), flags=re.IGNORECASE,
    )
    blockers: list[str] = []
    declared_solvent = environment.get("solvent") or (condition_contract or {}).get("solvent")
    declared_solvent_model = (
        environment.get("solvent_model") or (condition_contract or {}).get("solvent_model")
    )
    if canonical == "implicit_solution":
        if not declared_solvent:
            blockers.append("implicit_solution_solvent_identity_missing")
        if str(declared_solvent_model or "").upper() not in {"CPCM", "SMD"}:
            blockers.append("implicit_solution_model_must_be_cpcm_or_smd")
        expected = (
            f"{str(declared_solvent_model).upper()}({declared_solvent})".casefold()
            if declared_solvent and declared_solvent_model else None
        )
        retained = [f"{model.upper()}({solvent})".casefold() for model, solvent in solvent_tokens]
        if len(retained) != 1 or expected not in retained:
            blockers.append("orca_implicit_solvent_keyword_missing_or_mismatched")
    elif solvent_tokens:
        blockers.append("gas_or_nonsolution_route_contains_implicit_solvent_keyword")
    applicable = canonical in set(policy.allowed_environment_models) and not blockers
    return {
        "schema_version": 1,
        "kind": "surface_resolution_environment_applicability",
        "status": "applicable" if applicable else "environment_model_escalation_required",
        "declared_environment": environment,
        "condition_phase_model": (condition_contract or {}).get("phase_model"),
        "raw_model": raw,
        "canonical_model": canonical,
        "declaration_source": source,
        "allowed_environment_models": list(policy.allowed_environment_models),
        "declared_solvent": declared_solvent,
        "declared_solvent_model": declared_solvent_model,
        "orca_method_keywords": keywords,
        "orca_implicit_solvent_tokens": [
            f"{model.upper()}({solvent})" for model, solvent in solvent_tokens
        ],
        "blockers": blockers,
        "rate_claim_allowed": False,
        "lifetime_claim_allowed": False,
    }


def run_generic_reaction_path(
    route: RouteSpec,
    workdir: Path,
    *,
    prepare_only: bool = False,
    ncores: int = 1,
    method_keywords: list[str] | None = None,
    timeout_seconds: float | None = 14400.0,
    orientations: int = 3,
    neb_images: int = 8,
    barrierless_classifier: Callable[[dict], dict] | None = None,
    condition_contract: dict | None = None,
    surface_resolution_policy: SurfaceResolutionPolicy | None = None,
    retained_unexpected_proposal: dict | None = None,
    retained_stationary_workdir: Path | None = None,
    thermochemistry_temperature_K: float | None = None,
    historical_neb_seed_request: dict | None = None,
    historical_trial_resource_budget: dict | None = None,
    route_search_policy: Any | None = None,
    retained_fragment_request: dict | None = None,
    retained_partial_neb_request: dict | None = None,
    retained_planner_context: dict | None = None,
) -> dict:
    """Verify a resolved route without dispatching on a reaction family."""
    if route_search_policy is not None:
        if any(value is not None for value in (
                historical_neb_seed_request, historical_trial_resource_budget,
                retained_unexpected_proposal, retained_stationary_workdir)):
            raise ValueError("Adaptive search requires its own explicit retained-fragment context")
        from .route_search_execution import run_adaptive_reaction_path
        return run_adaptive_reaction_path(
            route, workdir, policy=route_search_policy, prepare_only=prepare_only,
            ncores=ncores, method_keywords=method_keywords,
            timeout_seconds=timeout_seconds, neb_images=neb_images,
            condition_contract=condition_contract,
            surface_resolution_policy=surface_resolution_policy,
            retained_fragment_request=retained_fragment_request,
            retained_partial_neb_request=retained_partial_neb_request,
            retained_planner_context=retained_planner_context,
            thermochemistry_temperature_K=thermochemistry_temperature_K,
        )
    if any(value is not None for value in (retained_fragment_request, retained_partial_neb_request,
                                           retained_planner_context)):
        raise ValueError("Independent retained evidence requires an explicit adaptive policy")
    workdir = Path(workdir)
    historical_assessment = None
    completed_irc_assessment = None
    if historical_neb_seed_request is not None:
        # Historical seed trials are new hypotheses, not in-place recovery.
        # Refuse before creating even a budget/context pointer in old sources.
        execution_source = Path(historical_neb_seed_request["source_execution_contract"]["path"]).resolve()
        reuse_request = historical_neb_seed_request.get("retained_ts_frequency")
        reuse_source = (Path(reuse_request["source_execution_contract"]["path"]).resolve()
                        if reuse_request is not None else None)
        completed_request = historical_neb_seed_request.get("retained_completed_irc")
        completed_source = (Path(completed_request["ts_frequency_request"]["source_execution_contract"]["path"]).resolve()
                            if completed_request is not None else None)
        if (workdir.resolve().is_relative_to(execution_source.parent)
                or (reuse_source is not None and workdir.resolve().is_relative_to(reuse_source.parent))
                or (completed_source is not None and workdir.resolve().is_relative_to(completed_source.parent))
                or (workdir.exists() and (not workdir.is_dir() or any(workdir.iterdir())))):
            raise ValueError("Historical NEB seed trial requires a new workspace outside the source workflow")
        if retained_unexpected_proposal is not None or retained_stationary_workdir is not None:
            raise ValueError("Historical NEB seed trial cannot combine distinct retained-path strategies")
        historical_assessment = inspect_historical_neb_ts_seed(
            route, historical_neb_seed_request, ncores=ncores,
            method_keywords=method_keywords, condition_contract=condition_contract,
        )
        if historical_neb_seed_request.get("retained_completed_irc") is not None:
            from .completed_irc_reuse import inspect_completed_irc_reuse
            completed_irc_assessment = inspect_completed_irc_reuse(
                route, historical_neb_seed_request["retained_completed_irc"],
                historical_admission=historical_assessment, ncores=ncores,
                method_keywords=method_keywords, surface_resolution_policy=surface_resolution_policy,
            )
            if not completed_irc_assessment.get("valid"):
                historical_assessment.update(valid=False, status="completed_irc_reuse_refused",
                    failure_reason=completed_irc_assessment.get("failure_reason"))
        if historical_assessment.get("valid"):
            source_cap = historical_assessment["source_checkpoint"]["timeout_seconds"]
            if (isinstance(timeout_seconds, bool) or timeout_seconds is None
                    or not math.isfinite(float(timeout_seconds)) or float(timeout_seconds) <= 0):
                raise ValueError("Historical TS trial requires a positive finite per-job timeout")
            if historical_trial_resource_budget is None:
                if float(timeout_seconds) > source_cap:
                    raise ValueError("Historical TS trial must retain a finite per-job timeout no larger than the source cap")
            else:
                # A new, explicitly budgeted hypothesis may have more execution
                # time. Never rewrite the source cap or relax any physical gate.
                budget = historical_trial_resource_budget
                duration = budget.get("overall_timeout_seconds")
                authorization = budget.get("authorization")
                if (budget.get("schema_version") != 1
                        or budget.get("kind") != "historical_trial_resource_budget"
                        or budget.get("source_execution_contract") != historical_assessment["source_execution_contract"]
                        or budget.get("previous_timeout_seconds_per_job") != source_cap
                        or type(budget.get("maximum_cores")) is not int
                        or budget.get("maximum_cores") != ncores
                        or isinstance(budget.get("timeout_seconds_per_job"), bool)
                        or budget.get("timeout_seconds_per_job") != timeout_seconds
                        or isinstance(duration, bool) or not isinstance(duration, (int, float))
                        or not math.isfinite(duration) or duration < float(timeout_seconds)
                        or not isinstance(authorization, str) or not authorization.strip()):
                    raise ValueError("Historical trial resource budget does not match its source and requested resources")
    elif historical_trial_resource_budget is not None:
        raise ValueError("Historical trial resource budget requires a historical source request")
    if thermochemistry_temperature_K is not None:
        thermochemistry_temperature_K = normalize_thermochemistry_temperature(thermochemistry_temperature_K)
        if not condition_contract or not math.isclose(
            float(condition_contract.get("temperature_K", float("nan"))), thermochemistry_temperature_K,
            rel_tol=0.0, abs_tol=1e-10,
        ):
            raise ValueError("Explicit thermochemistry temperature must match the reaction condition contract")
    workdir.mkdir(parents=True, exist_ok=True)
    surface_resolution_policy = surface_resolution_policy or SurfaceResolutionPolicy()
    budget_contract = _initialize_surface_budget(workdir, surface_resolution_policy)
    if orientations > surface_resolution_policy.maximum_orientations:
        raise ValueError(
            "Requested endpoint orientations exceed the immutable surface-resolution policy"
        )
    mapping = resolve_route_atom_mapping(route)
    retained_assessment = None
    if retained_unexpected_proposal is not None:
        retained_assessment = assess_unexpected_retained_path(
            route, retained_unexpected_proposal, workdir / "retained-path",
            method_keywords=method_keywords,
        )
        if retained_assessment.get("valid"):
            mapping = retained_assessment["mapping"]
    canonical: dict[str, Any] = {
        "status": "mapping_unresolved",
        "route_id": route.route_id,
        "reaction_equation": route.reaction_equation,
        "path_classification": "surface_unresolved",
        "endpoint_validation": {"status": "not_run", "valid": False},
        "mapping": mapping,
        "path_evidence": None,
        "rate_model": {"status": "rate_model_unavailable"},
        "propagation": {"status": "not_run"},
        "convergence": {"status": "not_converged"},
        "trajectory": None,
        "stationary_points": {"status": "not_run", "reactants": [], "products": [], "valid": False},
    }
    result: dict[str, Any] = {
        "schema_version": 1,
        "kind": "generic_orca_route_verification",
        "route": route.as_dict(),
        "route_verification": canonical,
        "orientations": [],
        "surface_resolution_policy": surface_resolution_policy.as_dict(),
        "surface_resolution_budget_contract": budget_contract,
    }
    if historical_trial_resource_budget is not None:
        result["historical_trial_resource_budget"] = dict(historical_trial_resource_budget)
    environment_applicability = _surface_environment_applicability(
        route, condition_contract, surface_resolution_policy, method_keywords,
    )
    result["environment_applicability"] = environment_applicability
    result["surface_resolution_scientific_context"] = _initialize_surface_context(
        workdir,
        route=route,
        mapping=mapping,
        environment_applicability=environment_applicability,
        policy=surface_resolution_policy,
        ncores=ncores,
        method_keywords=method_keywords,
        **({"thermochemistry_temperature_K": thermochemistry_temperature_K}
           if thermochemistry_temperature_K is not None else {}),
    )
    manifest = workdir / "route-verification.json"
    if retained_assessment is not None:
        result["retained_path_assessment"] = retained_assessment
        if not retained_assessment.get("valid"):
            canonical.update(status="irc_endpoint_validation_failed",
                             convergence={"status": "retained_path_reuse_unresolved"})
            _write_json(manifest, result)
            return {**result, "path": manifest}
    if environment_applicability["status"] != "applicable":
        canonical.update(
            status="environment_model_escalation_required",
            convergence={
                "status": "blocked_before_orca",
                "reason": "declared_environment_outside_surface_resolution_policy",
            },
        )
        _write_json(manifest, result)
        return {**result, "path": manifest}
    if not mapping.get("valid"):
        mapping_status = str(mapping.get("status") or "mapping_unresolved")
        electronic_state_failures = {
            "fragment_spin_contract_incomplete",
            "no_common_adiabatic_spin_surface",
            "multiple_common_spin_surfaces_require_selection",
            "declared_multiplicity_not_common",
        }
        canonical["status"] = (
            "electronic_state_escalation_required"
            if mapping_status in electronic_state_failures else mapping_status
        )
        if mapping_status in electronic_state_failures:
            canonical["electronic_state_evidence"] = {
                "status": mapping_status,
                "mapping_surface_assessment": mapping.get("surface") or mapping,
                "required_action": (
                    "declare_or_verify_an_applicable_adiabatic_surface_or_use_a_"
                    "specialized_spin_crossing_nonadiabatic_method"
                ),
                "rate_claim_allowed": False,
                "lifetime_claim_allowed": False,
            }
        canonical["convergence"] = {"status": "blocked_before_orca"}
        _write_json(manifest, result)
        return {**result, "path": manifest}
    selected_multiplicity = int(mapping["surface"]["selected_multiplicity"])
    execution_route = (
        route if route.multiplicity == selected_multiplicity
        else RouteSpec(**{**route.as_dict(), "multiplicity": selected_multiplicity})
    )
    canonical["selected_charge"] = int(mapping["surface"]["charge"])
    canonical["selected_multiplicity"] = selected_multiplicity
    separated_fragments = _route_uses_separated_fragment_endpoints(execution_route)
    canonical["endpoint_mode"] = (
        "separated_fragment_channels" if separated_fragments else "bound_endpoint_minima"
    )
    if historical_assessment is not None:
        result["historical_neb_seed_admission"] = historical_assessment
        if completed_irc_assessment is not None:
            result["completed_irc_reuse_admission"] = completed_irc_assessment
        if not historical_assessment.get("valid"):
            canonical.update(status="historical_ts_seed_refused",
                             convergence={"status": "blocked_before_orca"})
        else:
            canonical["stationary_points"] = historical_assessment["stationary_points"]
            canonical["endpoint_validation"] = historical_assessment["endpoint_validation"]
            seed = workdir / "historical-neb-seed.xyz"
            write_xyz(seed, historical_assessment["seed_elements"],
                      np.asarray(historical_assessment["seed_coordinates"]),
                      "Historical partial NEB image; new OptTS hypothesis, not a verified saddle")
            _write_json(workdir / "historical-neb-seed-request.json", historical_neb_seed_request)
            canonical.update(status="historical_ts_seed_prepared",
                             convergence={"status": "path_unresolved"})
            if prepare_only and completed_irc_assessment is not None:
                from .completed_irc_reuse import prepare_completed_irc_reuse
                prepared = prepare_completed_irc_reuse(completed_irc_assessment, workdir / "historical-seed-ts-trial")
                result["prepared_completed_irc_reuse"] = {"trajectory": str(prepared[3]), "reuse": prepared[4]}
            elif prepare_only and historical_assessment.get("retained_ts_frequency") is not None:
                irc_input, reuse = _prepare_reused_ts_irc(
                    execution_route, historical_assessment["retained_ts_frequency"],
                    workdir / "historical-seed-ts-trial", ncores=ncores, method_keywords=method_keywords,
                )
                result["prepared_reused_ts_irc"] = {"input": str(irc_input), "reuse": reuse}
            if not prepare_only:
                trial = _verify_ts_seed(
                    execution_route, seed, workdir / "historical-seed-ts-trial",
                    reactant_xyz=Path(historical_assessment["reactant_xyz"]),
                    product_xyz=Path(historical_assessment["product_xyz"]),
                    mapping_result=mapping, separated_fragments=False, ncores=ncores,
                    method_keywords=method_keywords, timeout_seconds=timeout_seconds,
                    surface_resolution_policy=surface_resolution_policy,
                    intended_coordinates=historical_assessment["intended_coordinates"],
                    retained_ts_frequency=historical_assessment.get("retained_ts_frequency"),
                    retained_completed_irc=completed_irc_assessment,
                )
                result["historical_neb_seed_trial"] = trial
                # A source mutation during new calculations cannot be repaired
                # by a final driver flag after rates have already been made.
                reassessment = inspect_historical_neb_ts_seed(
                    route, historical_neb_seed_request, ncores=ncores,
                    method_keywords=method_keywords, condition_contract=condition_contract,
                )
                result["historical_neb_seed_post_trial_admission"] = reassessment
                source_unchanged = reassessment.get("valid") is True
                if source_unchanged and completed_irc_assessment is not None:
                    from .completed_irc_reuse import inspect_completed_irc_reuse
                    reopened = inspect_completed_irc_reuse(
                        route, historical_neb_seed_request["retained_completed_irc"],
                        historical_admission=reassessment, ncores=ncores, method_keywords=method_keywords,
                        surface_resolution_policy=surface_resolution_policy,
                    )
                    result["completed_irc_post_trial_admission"] = reopened
                    source_unchanged = reopened.get("valid") is True
                canonical.update(status=trial["status"], path_classification=trial["path_classification"],
                                 transition_state=trial.get("transition_state"), irc=trial.get("irc"),
                                 trajectory=trial.get("trajectory"), path_evidence=trial.get("path_evidence"),
                                 convergence={"status": "path_unresolved"})
                if not source_unchanged:
                    canonical.update(status="historical_ts_seed_source_changed", path_classification="surface_unresolved",
                                     convergence={"status": "blocked_before_kinetics"})
                elif trial["status"] == "verified_transition_state_path":
                    refinement = _refine_reached_irc_stationary_points(
                        execution_route, trial, workdir / "reached-stationary-points", mapping,
                        ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
                        reference_endpoints={"reactants": historical_assessment["reactant_xyz"],
                                             "products": historical_assessment["product_xyz"]},
                        surface_resolution_policy=surface_resolution_policy,
                    )
                    canonical.update(refinement)
                    # Endpoint work can be lengthy. Reopen source contracts
                    # again after it, before exposing any rate input envelope.
                    final_source = inspect_historical_neb_ts_seed(
                        route, historical_neb_seed_request, ncores=ncores,
                        method_keywords=method_keywords, condition_contract=condition_contract,
                    )
                    source_unchanged = final_source.get("valid") is True
                    result["historical_neb_seed_post_refinement_admission"] = final_source
                    if source_unchanged and completed_irc_assessment is not None:
                        final_irc = inspect_completed_irc_reuse(
                            route, historical_neb_seed_request["retained_completed_irc"],
                            historical_admission=final_source, ncores=ncores, method_keywords=method_keywords,
                            surface_resolution_policy=surface_resolution_policy,
                        )
                        source_unchanged = final_irc.get("valid") is True
                        result["completed_irc_post_refinement_admission"] = final_irc
                    if not source_unchanged:
                        canonical.update(status="historical_ts_seed_source_changed", path_classification="surface_unresolved",
                                         convergence={"status": "blocked_before_kinetics"})
                    elif not canonical["stationary_points"].get("valid"):
                        canonical.update(status="stationary_point_validation_failed",
                                         convergence={"status": "blocked_before_kinetics"})
                    else:
                        canonical.update(
                            convergence={"status": "path_verified_rate_pending"},
                            arkane_inputs={"reactants": canonical["stationary_points"]["reactants"],
                                           "products": canonical["stationary_points"]["products"],
                                           "transition_state": trial["transition_state"]},
                        )
                result["historical_neb_seed_attempt_evidence"] = _write_surface_attempt_evidence(
                    workdir / "historical-seed-ts-trial", strategy=HISTORICAL_NEB_TS_SEED_POLICY,
                    status=canonical["status"], result=trial,
                    parameters={"seed_index": historical_assessment["seed_index"],
                                "historical_source_is_new_hypothesis_not_resume": True,
                                "source_unchanged_after_trial": source_unchanged,
                                "ncores": ncores, "timeout_seconds": timeout_seconds},
                    input_paths=[workdir / "historical-neb-seed-request.json", seed],
                    claims_permitted=(("verified_transition_state_path",)
                                      if canonical["status"] == "verified_transition_state_path" else ()),
                )
        result["surface_resolution_budget"] = _surface_budget_snapshot(workdir)
        _write_json(manifest, result)
        return {**result, "path": manifest}
    stationary: dict | None = None
    fragment_xyz = None
    if separated_fragments:
        # Each physical species is a real stationary point.  The assembled
        # collision/encounter geometry is only an atom-ordered NEB endpoint and
        # must not be required to optimize as a fictitious bound molecule.
        if retained_stationary_workdir is not None and retained_assessment is None:
            raise ValueError("External stationary reuse requires a revalidated retained path")
        stationary = (reuse_stationary_fragment_records(
            execution_route, retained_stationary_workdir, mapping, ncores=ncores,
            method_keywords=method_keywords, retained_endpoints={
                "reactants": Path(retained_assessment["reactant_xyz"]),
                "products": Path(retained_assessment["product_xyz"]),
            },
        ) if retained_stationary_workdir is not None and not prepare_only else _stationary_fragment_records(
            execution_route, workdir / "stationary-points", mapping,
            execute=not prepare_only, ncores=ncores, method_keywords=method_keywords,
            timeout_seconds=timeout_seconds,
            **({"retained_endpoints": {
                "reactants": Path(retained_assessment["reactant_xyz"]),
                "products": Path(retained_assessment["product_xyz"]),
            }} if retained_assessment else {}),
        ))
        canonical["stationary_points"] = stationary
        if not prepare_only and not stationary["valid"]:
            canonical.update(
                status="stationary_point_validation_failed",
                endpoint_validation={"status": "blocked_before_endpoint_assembly", "valid": False},
                convergence={"status": "blocked_before_neb"},
            )
            _write_json(manifest, result)
            return {**result, "path": manifest}
        if not prepare_only:
            thermochemistry = assemble_route_thermochemistry(
                execution_route, stationary, condition_contract,
            )
            loss_bound = condition_specific_forward_loss_bound(
                execution_route, thermochemistry, condition_contract,
            )
            canonical["endpoint_thermochemistry"] = thermochemistry
            canonical["forward_loss_upper_bound"] = loss_bound
            if loss_bound.get("status") == "forward_loss_below_retention_threshold_upper_bound":
                canonical.update(
                    status="thermochemistry_bounded_below_loss_threshold",
                    path_classification="no_path_required_for_bounded_route",
                    endpoint_validation={
                        "status": "validated_stationary_species",
                        "valid": True,
                        "mode": "separated_fragment_channels",
                    },
                    convergence={"status": "not_required_after_conservative_loss_bound"},
                    rate_model={
                        "status": "forward_rate_upper_bound_only",
                        "rate_claim_allowed": False,
                    },
                )
                _write_json(manifest, result)
                return {**result, "path": manifest}
        fragment_xyz = _stationary_fragment_xyz(stationary, optimized=not prepare_only)
    if retained_assessment is not None:
        # No new endpoint orientation/NEB search is warranted for a completed,
        # rechecked IRC. Stationary-species and downstream rate gates still apply.
        if stationary is None or prepare_only or not stationary.get("valid"):
            canonical.update(status="prepared" if prepare_only else "stationary_point_validation_failed")
        else:
            refs = retained_assessment["source_references"]
            transition_state = {
                "optimized_xyz": refs["ts_xyz"]["path"],
                "frequency_output": refs["frequency_output"]["path"],
                "reaction_coordinate_overlap": retained_assessment["mode_check"],
                "source_references": refs,
            }
            canonical.update(
                status="verified_transition_state_path", path_classification="verified_barriered_path",
                endpoint_validation={"status": "validated_stationary_species_and_retained_irc", "valid": True},
                transition_state=transition_state,
                trajectory=str(workdir / "retained-path" / "trajectory.xyz"),
                path_evidence={"kind": "revalidated_retained_ts_frequency_irc", "assessment": retained_assessment,
                               "scope": "declared_single_surface_method_not_electronic_structure_ground_truth"},
                arkane_inputs={"reactants": stationary["reactants"], "products": stationary["products"],
                               "transition_state": transition_state},
                convergence={"status": "path_verified_rate_pending"},
            )
        result["surface_resolution_budget"] = _surface_budget_snapshot(workdir)
        _write_json(manifest, result)
        return {**result, "path": manifest}
    seeds = build_endpoint_complex_seeds(
        execution_route, workdir / "endpoint-complexes", orientations=orientations,
        mapping_result=mapping, fragment_xyz=fragment_xyz,
    )
    result["endpoint_complexes"] = seeds
    if prepare_only:
        if stationary is None:
            stationary = _stationary_fragment_records(
                execution_route, workdir / "stationary-points", mapping, execute=False,
                ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
            )
        canonical.update(
            status="prepared", path_classification="not_computed",
            endpoint_validation={"status": "prepared", "valid": False},
            stationary_points=stationary,
            convergence={"status": "not_run"},
        )
        result["orientations"] = [
            {**record, "status": "prepared", "path_classification": "not_computed"}
            for record in seeds["orientations"]
        ]
        _write_json(manifest, result)
        return {**result, "path": manifest}

    def _resume_priority(seed_record: dict) -> tuple[int, float, int]:
        folder = Path(seed_record["reactant_xyz"]).parent / "path"
        output = folder / "neb-ts.out"
        band = folder / "neb-ts_MEP.allxyz"
        interrupted = (
            output.is_file() and not _normal_orca_output(output)
            and any(_NEB_LBFGS_ITERATION.search(line)
                    for line in output.read_text(errors="replace").splitlines())
        )
        validation = _validate_orca_neb_restart_band(
            band, expected_images=neb_images + 2,
        ) if interrupted else {"valid": False}
        if validation.get("valid"):
            return (0, -band.stat().st_mtime, int(seed_record["orientation"]))
        return (1, 0.0, int(seed_record["orientation"]))

    ordered_seeds = sorted(seeds["orientations"], key=_resume_priority)
    result["orientation_execution_order"] = [
        int(item["orientation"]) for item in ordered_seeds
    ]
    for seed_record in ordered_seeds:
        print(
            f"[STORCA] Route {route.route_id}: ORCA path orientation "
            f"{seed_record['orientation']}/{seeds['orientation_count']}", flush=True,
        )
        orientation_result = _run_generic_orientation(
            execution_route, seed_record, mapping_result=mapping,
            ncores=ncores, method_keywords=method_keywords,
            timeout_seconds=timeout_seconds, neb_images=neb_images,
            barrierless_classifier=barrierless_classifier,
            surface_resolution_policy=surface_resolution_policy,
            surface_bond_changes=seeds.get("bond_changes"),
            separated_fragments=separated_fragments,
            stationary_points=stationary,
        )
        orientation_root = Path(seed_record["reactant_xyz"]).parent / "path"
        neb_input_value = (orientation_result.get("neb") or {}).get("input")
        orientation_result["surface_resolution_attempt_evidence"] = (
            _write_surface_attempt_evidence(
                orientation_root,
                strategy="standard_neb_ts",
                status=str(orientation_result.get("status") or "surface_unresolved"),
                result=orientation_result,
                parameters={
                    "orientation": int(seed_record["orientation"]),
                    "neb_images": int(neb_images),
                    "separated_fragments": bool(separated_fragments),
                },
                input_paths=[
                    Path(seed_record["reactant_xyz"]),
                    Path(seed_record["product_xyz"]),
                    *([Path(str(neb_input_value))] if neb_input_value else []),
                ],
                failure_reason=orientation_result.get("failure_reason"),
                claims_permitted=(
                    ("verified_transition_state_path",)
                    if orientation_result.get("status") == "verified_transition_state_path"
                    else ()
                ),
            )
        )
        result["orientations"].append(orientation_result)
        _write_json(manifest, result)
        if orientation_result.get("status") == "verified_transition_state_path":
            break
        if (
            orientation_result.get("status")
            in {
                "intermediate_detected_requires_segmented_verification",
                "segment_path_unresolved",
                "segment_rate_model_unavailable",
                "verified_multistep_path",
                "intermediate_identity_unresolved",
                "surface_resolution_budget_exhausted",
                "dissociation_rate_model_required",
            }
            or (
                orientation_result.get("status") == "barrierless_rate_model_required"
                and orientation_result.get("path_classification")
                    == "multistep_intermediate_detected"
            )
        ):
            # A validated stationary intermediate is stronger evidence than
            # another equivalent full-endpoint NEB.  Continue with its bounded
            # segments instead of spending another orientation automatically.
            break
        current_diagnoses = _surface_diagnoses_from_records(result["orientations"])
        if len(current_diagnoses) >= surface_resolution_policy.independent_orientation_agreement:
            current_aggregate = aggregate_surface_diagnoses(
                current_diagnoses, policy=surface_resolution_policy,
            )
            intermediate_recovery_exhausted = (
                current_aggregate.status == "intermediate_candidate"
                and all(
                    (item.get("intermediate_evidence") or {}).get("status")
                    != "validated_intermediate_detected"
                    for item in result["orientations"]
                    if (item.get("surface_diagnosis") or {}).get("status") == "intermediate_candidate"
                )
            )
            if (
                intermediate_recovery_exhausted
                or (
                    current_aggregate.status in {
                        "endpoint_orientation_artifact", "single_saddle_candidate",
                    }
                    and current_aggregate.recommended_strategy == "bond_change_guided_scan"
                )
            ):
                matching = next(
                    item for item in result["orientations"]
                    if (item.get("surface_diagnosis") or {}).get("status") == current_aggregate.status
                )
                seed = next(
                    item for item in ordered_seeds
                    if int(item["orientation"]) == int(matching["orientation"])
                )
                selected_path_dir = Path(seed["reactant_xyz"]).parent / "path"
                guided = _run_bond_change_guided_ts_recovery(
                    execution_route,
                    reactant_xyz=selected_path_dir / "reactant.xyz",
                    product_xyz=selected_path_dir / "product.xyz",
                    bond_changes=seeds.get("bond_changes") or {},
                    mapping_result=mapping,
                    path_dir=selected_path_dir,
                    separated_fragments=separated_fragments,
                    ncores=ncores,
                    method_keywords=method_keywords,
                    timeout_seconds=timeout_seconds,
                    surface_resolution_policy=surface_resolution_policy,
                )
                guided_contract = guided.get("contract")
                guided_root = (
                    Path(str(guided_contract)).parent
                    if guided_contract else selected_path_dir / "bond-change-guided-recovery"
                )
                guided_scan_input = ((guided.get("scan") or {}).get("execution") or {}).get("input")
                if "surface_resolution_attempt_evidence" not in guided:
                    guided["surface_resolution_attempt_evidence"] = _write_surface_attempt_evidence(
                        guided_root,
                        strategy="bond_change_guided_scan",
                        status=str(guided.get("status") or "coordinate_scan_no_ts_seed"),
                        result=guided,
                        parameters={
                            "scan_plan": guided.get("scan_plan"),
                            "selected_orientation": int(matching["orientation"]),
                        },
                        input_paths=[
                            selected_path_dir / "reactant.xyz",
                            selected_path_dir / "product.xyz",
                            *([Path(str(guided_scan_input))] if guided_scan_input else []),
                        ],
                        failure_reason=guided.get("reason"),
                        claims_permitted=(
                            ("verified_transition_state_path",)
                            if guided.get("status") == "verified_transition_state_path" else ()
                        ),
                    )
                matching["bond_change_guided_recovery"] = guided
                result["surface_resolution_guided_recovery"] = {
                    "trigger": current_aggregate.as_dict(),
                    "selected_orientation": matching["orientation"],
                    "result": guided,
                    "equivalent_full_neb_stopped": True,
                }
                if guided.get("status") == "verified_transition_state_path":
                    matching.update(
                        status="verified_transition_state_path",
                        path_classification="verified_barriered_path",
                        transition_state=guided.get("transition_state"),
                        irc=guided.get("irc"),
                        trajectory=guided.get("trajectory"),
                        path_evidence=guided.get("path_evidence"),
                    )
                else:
                    matching.update(
                        status=str(guided.get("status") or "coordinate_scan_no_ts_seed"),
                        path_classification="surface_unresolved",
                    )
                _write_json(manifest, result)
                break
        classified = [
            item for item in result["orientations"]
            if item.get("path_classification") in {
                "barrierless_capture_candidate", "barrierless_dissociation_candidate", "barrierless_path_candidate",
            }
        ]
        reproducible_pair_found = False
        for left_index, left in enumerate(classified):
            left_energies = (left.get("path_evidence") or {}).get("energies_hartree") or []
            for right in classified[left_index + 1:]:
                right_energies = (right.get("path_evidence") or {}).get("energies_hartree") or []
                if (left.get("path_classification") == right.get("path_classification")
                        and len(left_energies) >= 2 and len(right_energies) >= 2
                        and abs((left_energies[-1] - left_energies[0])
                                - (right_energies[-1] - right_energies[0])) <= 0.02):
                    reproducible_pair_found = True
                    break
            if reproducible_pair_found:
                break
        if reproducible_pair_found:
            break
    verified = [item for item in result["orientations"] if item.get("status") == "verified_transition_state_path"]
    segmented = [
        item for item in result["orientations"]
        if (item.get("intermediate_evidence") or {}).get("status")
            == "validated_intermediate_detected"
    ]
    barrierless = [
        item for item in result["orientations"]
        if item.get("path_classification") in {
            "barrierless_capture_candidate", "barrierless_dissociation_candidate", "barrierless_path_candidate",
        }
    ]
    reproducible_barrierless: list[dict] = []
    barrierless_reproducibility: dict[str, Any] = {
        "status": "not_applicable", "required_independent_orientations": 2,
    }
    if barrierless and not verified:
        grouped: dict[str, list[dict]] = {}
        for item in barrierless:
            grouped.setdefault(item["path_classification"], []).append(item)
        for classification, group in sorted(grouped.items()):
            for left_index, left in enumerate(group):
                left_profile = (left.get("path_evidence") or {}).get("energies_hartree") or []
                if len(left_profile) < 2:
                    continue
                left_energy = float(left_profile[-1]) - float(left_profile[0])
                for right in group[left_index + 1:]:
                    right_profile = (right.get("path_evidence") or {}).get("energies_hartree") or []
                    if len(right_profile) < 2:
                        continue
                    right_energy = float(right_profile[-1]) - float(right_profile[0])
                    if abs(left_energy - right_energy) <= 0.02:
                        reproducible_barrierless = [left, right]
                        barrierless_reproducibility = {
                            "status": "reproducible",
                            "classification": classification,
                            "orientation_count": 2,
                            "orientation_ids": [left["orientation"], right["orientation"]],
                            "reaction_energy_spread_hartree": abs(left_energy - right_energy),
                            "required_independent_orientations": 2,
                        }
                        break
                if reproducible_barrierless:
                    break
            if reproducible_barrierless:
                break
        if not reproducible_barrierless:
            barrierless_reproducibility = {
                "status": "surface_unresolved",
                "reason": "barrierless_profile_not_reproduced_in_two_independent_orientations",
                "classified_orientation_count": len(barrierless),
                "required_independent_orientations": 2,
            }
    canonical["barrierless_reproducibility"] = barrierless_reproducibility
    diagnoses = _surface_diagnoses_from_records(result["orientations"])
    aggregate_diagnosis = aggregate_surface_diagnoses(
        diagnoses, policy=surface_resolution_policy,
    )
    canonical["surface_resolution"] = {
        "status": "diagnosed" if diagnoses else "not_diagnosed",
        "policy_id": surface_resolution_policy.policy_id,
        "orientation_diagnosis_count": len(diagnoses),
        "aggregate_diagnosis": aggregate_diagnosis.as_dict(),
        "rate_claim_allowed": False,
        "lifetime_claim_allowed": False,
    }
    precise_terminal = next((
        item for item in result["orientations"]
        if item.get("status") in surface_resolution_policy.terminal_abstention_states
    ), None)
    chosen = verified[0] if verified else reproducible_barrierless[0] if reproducible_barrierless else None
    if verified:
        refinement = _refine_reached_irc_stationary_points(
            execution_route, chosen, workdir / "reached-stationary-points", mapping,
            ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
            surface_resolution_policy=surface_resolution_policy,
        )
        canonical.update(refinement)
        stationary = refinement["stationary_points"]
    elif stationary is None:
        stationary = _stationary_fragment_records(
            execution_route, workdir / "stationary-points", mapping,
            execute=chosen is not None, ncores=ncores, method_keywords=method_keywords,
            timeout_seconds=timeout_seconds,
        )
    canonical["stationary_points"] = stationary
    if chosen is None:
        guided_summary = result.get("surface_resolution_guided_recovery") or {}
        guided_result = guided_summary.get("result") or {}
        if guided_result and guided_result.get("status") != "verified_transition_state_path":
            guided_status = str(guided_result.get("status") or "surface_unresolved")
            canonical.update(
                status=guided_status,
                path_classification="surface_unresolved",
                endpoint_validation={"status": "validated_for_guided_recovery", "valid": True},
                path_evidence=guided_result.get("path_evidence"),
                surface_resolution_attempt=guided_summary,
                selected_orientation=guided_summary.get("selected_orientation"),
                convergence={"status": "surface_resolution_incomplete", "reason": guided_status},
            )
        elif segmented:
            selected = segmented[0]
            segmented_status = str(
                selected.get("status") or "segment_path_unresolved"
            )
            canonical.update(
                status=segmented_status,
                path_classification="multistep_intermediate_detected",
                endpoint_validation=selected.get("endpoint_validation", {"valid": False}),
                path_evidence=selected.get("path_evidence"),
                intermediate_evidence=selected.get("intermediate_evidence"),
                selected_orientation=selected.get("orientation"),
                convergence={
                    "status": (
                        "segmented_paths_verified_rate_pending"
                        if (
                            ((selected.get("intermediate_evidence") or {})
                             .get("segmented_verification") or {})
                            .get("status") == "verified_segments_rate_generation_required"
                        ) else
                        "segmented_path_verification_incomplete"
                    ),
                    "reason": segmented_status,
                },
            )
            if segmented_status == "verified_multistep_path":
                canonical["multistep_path"] = (
                    (selected.get("intermediate_evidence") or {})
                    .get("segmented_verification")
                )
        elif precise_terminal:
            precise_status = str(precise_terminal["status"])
            canonical.update(
                status=precise_status,
                path_classification=str(
                    precise_terminal.get("path_classification") or "surface_unresolved"
                ),
                endpoint_validation=precise_terminal.get(
                    "endpoint_validation", {"status": "not_completed", "valid": False},
                ),
                path_evidence=precise_terminal.get("path_evidence"),
                classification_evidence=precise_terminal.get("classification_evidence"),
                method_selection=precise_terminal.get("method_selection"),
                intermediate_evidence=precise_terminal.get("intermediate_evidence"),
                surface_resolution_budget=precise_terminal.get("surface_resolution_budget"),
                selected_orientation=precise_terminal.get("orientation"),
                convergence={"status": "surface_resolution_incomplete", "reason": precise_status},
            )
        else:
            canonical.update(
                status="surface_unresolved",
                endpoint_validation={"status": "no_validated_orientation", "valid": False},
                convergence={"status": "path_unresolved"},
            )
    else:
        canonical.update(
            path_classification=chosen["path_classification"],
            endpoint_validation=chosen.get("endpoint_validation", {"valid": False}),
            path_evidence=chosen.get("path_evidence"),
            trajectory=chosen.get("trajectory"),
            selected_orientation=chosen["orientation"],
            transition_state=chosen.get("transition_state"),
            irc=chosen.get("irc"),
        )
        if not stationary["valid"]:
            canonical.update(
                status="stationary_point_validation_failed",
                convergence={"status": "blocked_before_kinetics"},
            )
        elif verified:
            canonical.update(
                status="verified_transition_state_path",
                convergence={"status": "path_verified_rate_pending"},
                arkane_inputs={
                    "reactants": stationary["reactants"],
                    "products": stationary["products"],
                    "transition_state": canonical["transition_state"],
                },
            )
        else:
            canonical.update(
                status="barrierless_rate_model_required",
                convergence={"status": "route_physics_classified_rate_pending"},
            )
    result["surface_resolution_budget"] = _surface_budget_snapshot(workdir)
    _write_json(manifest, result)
    return {**result, "path": manifest}


def _run_neb(route: RouteSpec, explanation_dir: Path, *, prepare_only: bool,
             ncores: int, method_keywords: list[str] | None,
             timeout_seconds: float | None, neb_images: int) -> dict:
    endpoint_dir = explanation_dir / "endpoints"
    seeds = build_mapped_endpoint_seeds(route, endpoint_dir)
    if prepare_only:
        return {"status": "prepared", "endpoints": seeds}
    reactant_result = _optimize_endpoint(
        Path(seeds["reactant_xyz"]), endpoint_dir / "reactant", charge=route.charge,
        multiplicity=route.multiplicity, ncores=ncores, method_keywords=method_keywords,
        timeout_seconds=timeout_seconds,
    )
    product_result = _optimize_endpoint(
        Path(seeds["product_xyz"]), endpoint_dir / "product", charge=route.charge,
        multiplicity=route.multiplicity, ncores=ncores, method_keywords=method_keywords,
        timeout_seconds=timeout_seconds,
    )
    if not reactant_result["frequency_check"]["IsMinimum"] or not product_result["frequency_check"]["IsMinimum"]:
        raise RuntimeError("Both NEB endpoints must be ORCA local minima")
    path_dir = explanation_dir / "path"
    path_dir.mkdir(parents=True, exist_ok=True)
    reactant_xyz = path_dir / "reactant.xyz"
    product_xyz = path_dir / "product.xyz"
    shutil.copy2(reactant_result["optimized_xyz"], reactant_xyz)
    shutil.copy2(product_result["optimized_xyz"], product_xyz)
    neb_input = create_orca_neb_ts_input(
        reactant_xyz, product_xyz, charge=route.charge, multiplicity=route.multiplicity,
        ncores=ncores, nimages=neb_images, method_keywords=method_keywords, preopt_ends=True,
    )
    neb_execution = _run_or_resume(neb_input, timeout_seconds=timeout_seconds)
    ts_xyz = _locate_neb_ts(path_dir, Path(neb_execution.get("input") or neb_input).stem)
    if ts_xyz is None:
        raise RuntimeError("ORCA NEB-TS completed without a retained TS geometry")
    retained_ts = path_dir / "transition-state.xyz"
    if ts_xyz.resolve() != retained_ts.resolve():
        shutil.copy2(ts_xyz, retained_ts)
    frequency_input = create_orca_input(
        retained_ts, route.charge, route.multiplicity, freq=True, label="ts-frequency",
        ncores=ncores, method_keywords=method_keywords,
    )
    frequency_execution = _run_or_resume(frequency_input, timeout_seconds=timeout_seconds)
    executed_frequency_input = Path(frequency_execution.get("input") or frequency_input)
    frequency_check = frequency_stability_check(
        Path(frequency_execution.get("output") or executed_frequency_input.with_suffix(".out")),
    )
    if frequency_check["NumImaginary"] != 1:
        raise RuntimeError(
            f"Transition-state frequency check found {frequency_check['NumImaginary']} significant imaginary modes; exactly one is required"
        )
    hessian = executed_frequency_input.with_suffix(".hess")
    irc_input = create_orca_irc_input(
        retained_ts, charge=route.charge, multiplicity=route.multiplicity,
        hessian_file=hessian if hessian.is_file() else None, ncores=ncores,
        method_keywords=method_keywords,
    )
    irc_execution = _run_or_resume(irc_input, timeout_seconds=timeout_seconds)
    trajectory = _locate_irc_trajectory(
        path_dir, Path(irc_execution.get("input") or irc_input).stem,
    )
    if trajectory is None:
        raise RuntimeError("ORCA IRC completed without the combined forward/backward trajectory")
    retained_trajectory = path_dir / "trajectory.xyz"
    if trajectory.resolve() != retained_trajectory.resolve():
        shutil.copy2(trajectory, retained_trajectory)
    endpoint_validation = validate_trajectory_endpoints(retained_trajectory, reactant_xyz, product_xyz)
    if not endpoint_validation["valid"]:
        raise RuntimeError("The IRC endpoints do not match the declared reactant and product")
    return {
        "status": "verified_transition_state_path",
        "endpoints": {
            **seeds,
            "optimized_reactant_xyz": str(reactant_xyz),
            "optimized_product_xyz": str(product_xyz),
        },
        "neb_execution": neb_execution,
        "frequency_execution": frequency_execution,
        "frequency_check": frequency_check,
        "irc_execution": irc_execution,
        "trajectory": str(retained_trajectory),
        "endpoint_validation": endpoint_validation,
    }


def run_decomposition_explanation(
    report_path: Path,
    *,
    output_dir: Path | None = None,
    route_id: str | None = None,
    prepare_only: bool = False,
    charge: int | None = None,
    multiplicity: int | None = None,
    ncores: int = 1,
    method_keywords: list[str] | None = None,
    timeout_seconds: float | None = 14400.0,
    scan_steps: int = 20,
    neb_images: int = 8,
) -> dict:
    """Select, calculate, validate, and visualize one decomposition route."""
    report_path = Path(report_path)
    source_dir = report_path if report_path.is_dir() else report_path.parent
    explanation_dir = Path(output_dir) if output_dir else source_dir / "explanation"
    explanation_dir.mkdir(parents=True, exist_ok=True)
    prepared = prepare_explanation(report_path, route_id)
    route_data = prepared["selected_route"]
    if charge is not None:
        route_data["charge"] = charge
    if multiplicity is not None:
        route_data["multiplicity"] = multiplicity
    route = RouteSpec(**route_data)
    result: dict[str, Any] = {
        **prepared,
        "selected_route": route.as_dict(),
        "status": "running",
        "method_keywords": method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"],
        "resources": {
            "ncores": ncores,
            "timeout_seconds_per_orca_job": timeout_seconds,
            "scan_steps": scan_steps,
            "neb_images": neb_images,
        },
    }
    manifest = explanation_dir / "decomposition-explanation.json"
    _write_json(manifest, result)
    visual_dir = explanation_dir / "visuals"
    storyboard = render_candidate_storyboard(
        route.parent_smiles, list(route.product_smiles), visual_dir / "candidate-storyboard.png",
        title=f"{route.route_id} ({route.source})",
    ) if route.product_smiles else None
    mapping = validate_route_mapping(route)
    result["mapping_validation"] = mapping
    if route.protocol == "unsupported_atom_mapping" or not mapping["valid"]:
        result["status"] = "unsupported_route_physics"
        result["failure_reason"] = route.limitation or mapping["status"]
        if storyboard:
            result["visuals"] = {
                "storyboard": str(storyboard),
                "manifest": str(write_visual_manifest(
                    visual_dir / "visuals.json", route_id=route.route_id,
                    evidence_level="candidate_storyboard", storyboard=storyboard,
                )),
            }
        _write_json(manifest, result)
        return {**result, "result_json": manifest}
    try:
        if route.protocol == "relaxed_dissociation_scan":
            path_result = _run_dissociation(
                route, explanation_dir, prepare_only=prepare_only, ncores=ncores,
                method_keywords=method_keywords, timeout_seconds=timeout_seconds, scan_steps=scan_steps,
            )
        elif route.protocol in {"neb_ts_irc", "generic_endpoint_path"}:
            generic = run_generic_reaction_path(
                route, explanation_dir / "generic-path", prepare_only=prepare_only,
                ncores=ncores, method_keywords=method_keywords,
                timeout_seconds=timeout_seconds, orientations=3, neb_images=neb_images,
            )
            canonical = generic["route_verification"]
            first_seed = generic["endpoint_complexes"]["orientations"][0]
            path_result = {
                "status": canonical["status"],
                "endpoints": {
                    "reactant_xyz": first_seed["reactant_xyz"],
                    "product_xyz": first_seed["product_xyz"],
                },
                "trajectory": canonical.get("trajectory"),
                "generic_verification": generic,
            }
            if not prepare_only and not path_result["trajectory"]:
                raise RuntimeError(
                    f"Generic route verification stopped at {canonical['status']} without a validated trajectory"
                )
        else:
            raise ValueError(f"Unsupported route protocol: {route.protocol}")
        result["path_calculation"] = path_result
        if prepare_only:
            preview = _render_preview(route, path_result["endpoints"], visual_dir, storyboard)
            result.update(status="prepared", visuals=preview)
        else:
            trajectory = Path(path_result["trajectory"])
            evidence = path_result["status"]
            animation = render_xyz_trajectory_gif(
                trajectory, visual_dir / "decomposition.gif",
                title=f"{route.route_id}: decomposition path", evidence_label=evidence.replace("_", " ").upper(),
                energies=path_result.get("energies_hartree"),
            )
            energy_profile = None
            if path_result.get("energies_hartree"):
                energy_profile = render_energy_profile(
                    path_result["energies_hartree"], visual_dir / "energy-profile.png"
                )
            visual_manifest = write_visual_manifest(
                visual_dir / "visuals.json", route_id=route.route_id, evidence_level=evidence,
                storyboard=storyboard, animation=animation, energy_profile=energy_profile,
                trajectory=trajectory,
            )
            result.update(
                status=evidence,
                visuals={
                    "storyboard": str(storyboard) if storyboard else None,
                    "animation": str(animation),
                    "energy_profile": str(energy_profile) if energy_profile else None,
                    "trajectory": str(trajectory),
                    "manifest": str(visual_manifest),
                },
            )
    except Exception as error:
        result["status"] = "failed"
        result["failure_reason"] = str(error)
        # A failed computation may still provide a clearly labelled preview if
        # endpoint generation completed before the failure.
        endpoints = (result.get("path_calculation") or {}).get("endpoints")
        if storyboard and endpoints:
            result["visuals"] = _render_preview(route, endpoints, visual_dir, storyboard)
    _write_json(manifest, result)
    return {**result, "result_json": manifest}
