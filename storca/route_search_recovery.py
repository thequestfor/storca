"""Narrow recovery evidence adapters; no runner, ledger, or relaxed validators.

Paid-stage replay is allowed to read/reanalyze existing owned work. It cannot
reserve a new calculation. When a composite existing validator reaches an
unpaid descendant, the caller gets an explicit request to plan that remainder.
"""
from __future__ import annotations

from pathlib import Path

from .rate_transfer_control import (TransferControlError, active_execution_control,
                                     classify_native_input)
from .rate_transfer_preparation import binding, _reopen


class UnpaidStageRequired(TransferControlError):
    """A retained composite stage reached an unperformed native descendant."""
    def __init__(self, request):
        self.request = request
        super().__init__("Paid-stage replay reached an unpaid native descendant")


def owned_execution_facts(control, reservation_ids=None):
    """Reopen owned native sources, not a wrapper's mutable success label."""
    from .reaction_path import _orca_job_input_evidence
    snapshot = control.snapshot()
    if snapshot["pending"]:
        raise TransferControlError("Paid-stage ownership remains unresolved")
    ids = sorted(snapshot["reservations"] if reservation_ids is None else reservation_ids)
    if len(set(ids)) != len(ids) or any(type(item) is not int or item not in snapshot["reservations"] for item in ids):
        raise TransferControlError("Paid-stage reservation identity is missing or duplicated")
    facts, sources = [], {}
    for identifier in ids:
        reservation = snapshot["reservations"][identifier]
        native_input = Path(reservation["input"])
        digest = _orca_job_input_evidence(native_input)[0].hexdigest()
        control.validate_reuse(native_input, digest)
        terminal = snapshot["terminals"][identifier]
        for source in terminal["sources"]:
            _reopen(source)
            sources[source["path"]] = source
        facts.append({"attempt_id": identifier, "category": reservation["category"],
                      "input": reservation["input"], "input_sha256": reservation["input_sha256"],
                      "native_outcome": terminal["outcome"],
                      "elapsed_seconds": terminal["elapsed_seconds"],
                      "physically_validated_by_execution_receipt": False})
    return {"native_reservation_ids": ids, "owned_execution_facts": facts,
            "sources": [sources[path] for path in sorted(sources)]}


def safe_replay_stage(execute, control, *, reservation_ids=None):
    """Run an existing composite reader/reuser, refusing EVERY new reservation.

    Requires the caller's already active, exclusively owned native control.
    Completed native caches remain usable. A stopped failed job is not relabelled
    successful; any existing append-only retry becomes an unpaid remainder.
    No attempt is started, refunded, or removed by this helper.
    """
    if active_execution_control() is not control:
        raise TransferControlError("Paid replay requires the existing active native controller")
    before = control.snapshot()
    facts = owned_execution_facts(control, reservation_ids)
    requested = []
    sentinel = object()
    previous = control.__dict__.get("reserve", sentinel)

    def deny_reservation(input_path, input_sha256, *, timeout_seconds, ncores):
        native_input = Path(input_path).resolve(strict=True)
        matches = [item for item in before["reservations"].values()
                   if item["input"] == str(native_input)]
        if matches:
            raise TransferControlError("Paid replay cannot repeat an attempted native basename")
        request = {"input": binding(native_input), "input_sha256": input_sha256,
                   "category": classify_native_input(native_input.read_text()),
                   "ncores": ncores, "requested_timeout_seconds": timeout_seconds,
                   "reason": "existing_stage_requires_unpaid_descendant",
                   "native_jobs_launched": 0}
        requested.append(request)
        raise UnpaidStageRequired(request)

    control.reserve = deny_reservation
    try:
        try:
            trial = execute()
        except UnpaidStageRequired as error:
            result = {"status": "paid_stage_partial", "trial": None,
                      "unpaid_remainder_request": error.request}
        else:
            # Some legacy physical adapters catch broad Exception. Even if they
            # swallow our refusal, never call that a completed stage replay.
            result = {"status": "paid_stage_partial" if requested else "paid_stage_replayed",
                      "trial": trial,
                      "unpaid_remainder_request": requested[0] if requested else None}
    finally:
        if previous is sentinel:
            del control.reserve
        else:
            control.reserve = previous
    after = control.snapshot()
    if after != before:
        raise TransferControlError("Native journal changed during dispatch-forbidden paid replay")
    owned_execution_facts(control, facts["native_reservation_ids"])
    return {**result, **facts, "new_native_jobs": 0,
            "physical_validation_passed": False}


def _local_barriered_gates(trial, *, require_reached):
    ts, irc = trial.get("transition_state") or {}, trial.get("irc") or {}
    evidence = trial.get("path_evidence") or {}
    gates = bool(ts.get("frequency_check", {}).get("NumImaginary") == 1
                 and ts.get("reaction_coordinate_overlap", {}).get("status") == "matched"
                 and evidence.get("intended_mode_match") is True
                 and irc.get("completion", {}).get("valid") is True
                 and irc.get("endpoint_validation", {}).get("valid") is True)
    if not gates or not require_reached:
        return gates
    stationary = trial.get("stationary_points") or {}
    rows = [row for side in ("reactants", "products") for row in stationary.get(side, [])]
    return bool(trial.get("status") == "verified_transition_state_path"
                and trial.get("path_classification") == "verified_barriered_path"
                and stationary.get("valid") is True
                and stationary.get("geometry_origin") == "locally_refined_reached_irc_endpoints"
                and stationary.get("reactants") and stationary.get("products") and rows
                and all(row.get("status") in {"validated_minimum", "validated_atom"} for row in rows))


def classify_recovery(trial, action_kind):
    """Suggest evidence-eligible next actions, never grant route/rate admission.

    Results are proposal records. Caller still checks exact sources, method,
    environment, physical applicability, meaningful difference and native cost.
    """
    if not isinstance(trial, dict) or not isinstance(action_kind, str):
        raise ValueError("A retained trial and explicit action kind are required")
    actions = []
    def propose(kind, reason):
        actions.append({"kind": kind, "reason": reason})
    status = trial.get("status")
    ts, irc = trial.get("transition_state") or {}, trial.get("irc") or {}
    if _local_barriered_gates(trial, require_reached=True):
        return {"status": "complete_barriered_evidence", "actions": [],
                "rate_claim_allowed": False}
    if _local_barriered_gates(trial, require_reached=False):
        propose("validate_reached_minima", "TS, mapped mode and both IRC branches passed; actual reached species remain unvalidated")
    elif status == "irc_path_incomplete":
        completion = irc.get("completion") or {}
        directions = completion.get("directions") or {}
        if (completion.get("normal_termination") is True
                and set(directions) == {"forward", "backward"}
                and all(item.get("valid") is True or item.get("iteration_limit_reached") is True
                        for item in directions.values())
                and any(item.get("iteration_limit_reached") is True and item.get("valid") is not True
                        for item in directions.values())):
            propose("continue_irc", "Normally terminated capped branch qualifies for the existing bounded downhill continuation")
    if not actions and trial.get("ts_seed") and status not in {
            "ts_mode_mismatch", "ts_frequency_validation_failed", "irc_endpoint_validation_failed"}:
        propose("verify_ts", "Retained candidate can enter the existing full TS/mode/IRC validator")
    if not actions and ts.get("optimized_xyz") and status in {
            "running", "path_execution_failed", "native_execution_failed", "awaiting_verification"}:
        propose("continue_verification", "Reanalyze owned completed stages and explicitly plan missing descendants")
    intermediate = trial.get("intermediate_evidence") or trial.get("intermediate_validation") or {}
    if any(row.get("status") == "validated_intermediate_minimum"
           and (row.get("intermediate_identity") or {}).get("status") == "resolved"
           for row in intermediate.get("candidates", [])):
        propose("verify_segments", "A physically validated identity-resolved intermediate requires independently verified segments")
    qualification = trial.get("band_restart_qualification") or {}
    warning = (trial.get("neb") or {}).get("path_quality_warning") or {}
    if (qualification.get("valid") is True and qualification.get("useful_restart") is True
            and qualification.get("source_bindings") and not warning):
        propose("resume_band", "Retained band has explicit useful-restart qualification and exact source bindings")
    return {"status": "recoverable" if actions else "distinct_hypothesis_or_escalation_required",
            "actions": actions, "rate_claim_allowed": False}


def validate_multistep_native_route(trial):
    """Refuse aggregate success labels without every segment's full evidence.

    Existing segmentation often establishes only TS/IRC and initial endpoint
    rates. That is not actual reached-species validation and does not pass here.
    Native scientific readers must have produced and source-bound these facts.
    """
    segmented = (trial.get("segmented_verification") or
                 (trial.get("intermediate_evidence") or {}).get("segmented_verification") or {})
    segments = segmented.get("segments") or []
    failures = []
    if len(segments) != 2:
        failures.append("exactly_two_supported_elementary_segments_required")
    elif [item.get("label") for item in segments] != [
            "reactants-to-intermediate", "intermediate-to-products"]:
        failures.append("ordered_intermediate_chain_required")
    else:
        first, second = [item.get("native_candidate") or {} for item in segments]
        if not first.get("products") or first.get("products") != second.get("reactants"):
            failures.append("disconnected_intermediate_identity")
    pairs = set()
    for index, segment in enumerate(segments):
        prefix = f"segment_{index + 1}:"
        pair = segment.get("endpoint_pair_hash")
        if not pair or pair in pairs:
            failures.append(prefix + "missing_or_repeated_endpoint_identity")
        pairs.add(pair)
        if segment.get("status") != "verified_transition_state_path":
            failures.append(prefix + "elementary_route_not_verified")
        if (segment.get("segment_mapping") or {}).get("valid") is not True:
            failures.append(prefix + "mapping_unresolved")
        verification = segment.get("verification") or {}
        if not _local_barriered_gates(verification, require_reached=True):
            failures.append(prefix + "full_mode_IRC_actual_reached_evidence_required")
    return {"valid": not failures, "status": "verified_multistep_route" if not failures else
            "multistep_route_unresolved", "segment_count": len(segments),
            "failures": failures, "rate_claim_allowed": False,
            "all_required_segments_verified": not failures}


def continue_qualifying_irc(source_input, *, control, timeout_seconds):
    """Use the established single-extension implementation, then its assembler.

    Completion is not endpoint identity, a verified route, or kinetic evidence.
    This helper is an explicitly paid action unless enclosed in safe replay.
    """
    if active_execution_control() is not control:
        raise TransferControlError("IRC recovery requires the active native controller")
    from . import reaction_path as rp
    source_input = Path(source_input).resolve(strict=True)
    snapshot = control.snapshot()
    matches = [identifier for identifier, row in snapshot["reservations"].items()
               if row["input"] == str(source_input)]
    if not matches:
        raise TransferControlError("IRC recovery source has no owned native reservation")
    if (classify_native_input(source_input.read_text()) != "irc"
            or any(snapshot["reservations"][identifier]["category"] != "irc" for identifier in matches)):
        raise ValueError("IRC recovery requires an actual owned IRC source")
    facts = owned_execution_facts(control, matches)
    completion = rp._orca_irc_completion(source_input.with_suffix(".out"))
    directions = completion.get("directions") or {}
    if (completion.get("normal_termination") is not True
            or set(directions) != {"forward", "backward"}
            or not all(row.get("valid") is True or row.get("iteration_limit_reached") is True
                       for row in directions.values())):
        raise ValueError("IRC branches do not qualify for bounded downhill continuation")
    branches = [branch for branch in ("forward", "backward")
                if directions[branch].get("valid") is not True]
    if not branches:
        raise ValueError("IRC is already complete; use offline endpoint validation")
    continuations = [rp.run_irc_branch_continuation(
        source_input, branch=branch, timeout_seconds=timeout_seconds) for branch in branches]
    assembly = rp.assemble_completed_irc_branches(source_input)
    return {"status": "continued_irc_awaiting_endpoint_validation", "source_completion": completion,
            "branch_continuations": continuations, "branch_assembly": assembly,
            "source_execution_facts": facts, "route_verified": False, "rate_claim_allowed": False}


__all__ = ["UnpaidStageRequired", "owned_execution_facts", "safe_replay_stage",
           "classify_recovery", "validate_multistep_native_route", "continue_qualifying_irc"]
