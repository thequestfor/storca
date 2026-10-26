"""Resume the frozen angular experiment with loader tracing diagnostic-only.

This explicit reader adapter changes no frozen source, native input, runtime
environment, physical test, attempt budget or schedule. Its immutable extension
is bound into every new runtime analysis and final result. Historical receipts
are still replayed with their original reader. Completed +120 work is reused.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from storca.configured_mpi_runtime import configured_mpi_receipt

spec = importlib.util.spec_from_file_location(
    "runtime_policy_frozen_candidate", Path(__file__).with_name("continue_oh_methane_angular_profile.py"))
candidate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)
old, pilot, state = candidate.old, candidate.pilot, candidate.state
REVISION = "configured_runtime_loader_diagnostic_v1"
EXTENSION = "runtime-policy-extension-v1.json"
POLICY = {
    "revision": REVISION,
    "loader_trace_required_for_acceptance": False,
    "selected_runtime_inventory_required": True,
    "owned_four_core_normal_completion_required": True,
    "physical_policy_changed": False,
    "native_inputs_or_environment_changed": False,
    "budget_or_schedule_changed": False,
    "historical_receipts_reinterpreted": False,
    "rate_claim_allowed": False,
}


def validate_extension(folder, contract, ledger):
    """An extension anchors a journal prefix, not an unchanging future head."""
    binding = pilot.read_json(folder / "runtime-policy-extension-binding.json")
    path = pilot.reopen(binding)
    if path != folder / EXTENSION:
        raise ValueError("Runtime extension location differs")
    extension = pilot.read_json(path)
    if (extension.get("schema_version") != 1 or extension.get("policy") != POLICY
            or extension.get("contract") != pilot.artifact(folder / "contract.json")):
        raise ValueError("Runtime policy extension differs")
    old.reopen_bindings(extension)
    state.reopen_inventory(extension["source_inventory"])
    state.reopen_inventory(contract["source_inventory"])
    state.reopen_inventory(contract["runtime"]["inventory"])
    expected_sources = {str(Path(__file__).resolve()), str(ROOT / "storca/configured_mpi_runtime.py")}
    if set(extension["source_inventory"]) != expected_sources:
        raise ValueError("Exact runtime reader sources required")
    anchor = extension["attempt_evidence_at_revision"]
    if not anchor or ledger.evidence()[:len(anchor)] != anchor:
        raise ValueError("Runtime extension's retained journal prefix differs")
    if extension["whole_increment_attempts_at_revision"] != {"saddles": 8, "NumFreq": 3}:
        raise ValueError("Runtime extension cannot reset charged attempts")
    if Path(contract["orca"]["path"]).resolve() != Path(pilot.find_orca()).resolve():
        raise ValueError("Configured ORCA executable differs")
    pilot.reopen(contract["orca"])
    return binding


def prepare_extension(folder, contract, prepared, ledger):
    """Publish the narrow reader correction once, without changing a journal."""
    with prepared["parent_ledger"].controller(), ledger.controller():
        path = folder / EXTENSION
        if not path.exists():
            if (folder / "result.json").exists() or ledger.counts() != {"saddles": 8, "NumFreq": 3}:
                raise ValueError("Expected retained +120 attention state required")
            attempts = ledger.reservations()
            execution = ledger.execution(attempts[0]) if len(attempts) == 1 else None
            if (len(attempts) != 1 or attempts[0]["point_id"] != "point-p120"
                    or attempts[0]["kind"] != "saddles" or ledger.analysis(attempts[0]) is not None
                    or execution is None or execution.get("status") != "normal"):
                raise ValueError("One paid normal, unanalyzed +120 attempt required")
            attention = folder / "checkpoints/0000.json"
            checkpoint = pilot.read_json(attention)
            if (checkpoint.get("failure_reason") != "ValueError: Private MPI loading at actual gradient workers not proved"
                    or checkpoint.get("whole_increment_attempts") != ledger.counts()):
                raise ValueError("Expected loader-reader attention checkpoint required")
            old.reopen_bindings(ledger.execution(attempts[0]))
            old.atomic_save(path, {
                "schema_version": 1, "kind": "OH_CH4_explicit_runtime_reader_extension",
                "policy": POLICY, "contract": pilot.artifact(folder / "contract.json"),
                "attention_checkpoint": pilot.artifact(attention),
                "attempt_evidence_at_revision": ledger.evidence(),
                "whole_increment_attempts_at_revision": ledger.counts(),
                "source_inventory": state.inventory({Path(__file__), ROOT / "storca/configured_mpi_runtime.py"}),
                "reason": "Native +120 completed; interleaved diagnostic stderr is not an execution or chemistry failure",
            })
        seal = folder / "runtime-policy-extension-binding.json"
        if not seal.exists():
            # Exclusive publication supports an interruption between these two
            # immutable files; all semantic/source checks still run below.
            old.atomic_save(seal, pilot.artifact(path))
        return validate_extension(folder, contract, ledger)


@contextmanager
def runtime_reader_adapter(folder, contract, ledger):
    """Explicit scoped callbacks: diagnostics, dispatch, checkpoint, closeout.

    Only current journal-owned normal outputs get the revised runtime reader.
    The old physical analyzers and historical source replay run unmodified.
    This is deliberately not a process-wide change to production STORCA.
    """
    previous_witness, previous_dispatch = old.private_mpi_runtime_witness, candidate.dispatch
    previous_checkpoint, previous_finish = old.checkpoint, old.finish
    binding = validate_extension(folder, contract, ledger)

    def witness(native, prefix):
        digest = hashlib.sha256(native.encode()).hexdigest()
        for attempt in ledger.reservations():
            execution = ledger.execution(attempt)
            inp = Path(attempt["metadata"]["input_path"])
            outputs = [] if execution is None else [record for record in execution.get("artifacts", [])
                        if record.get("path") == str(inp.with_suffix(".out"))]
            if len(outputs) != 1 or outputs[0].get("sha256") != digest:
                continue
            validate_extension(folder, contract, ledger)
            old.reopen_bindings(execution)
            if execution.get("status") != "normal" or str(prefix) != contract["runtime"]["prefix"]:
                raise ValueError("Current native execution/runtime binding differs")
            if pilot.artifact(inp)["sha256"] != attempt["metadata"]["input_sha256"]:
                raise ValueError("Current native reserved input differs")
            outcome = old.native_outcome(inp, inp.read_text())
            if outcome["status"] != "normal":
                raise ValueError("Owned four-core native completion required")
            receipt = configured_mpi_receipt(native, prefix, contract["selected_environment"],
                                             trace_reader=previous_witness)
            return {**receipt, "runtime_policy_extension": binding,
                    "native_outcome": outcome,
                    "selected_runtime_inventory": contract["runtime"]["inventory"]}
        # In particular, keep the reused +60 receipt exactly identical.
        return previous_witness(native, prefix)

    def dispatch(*args, **kwargs):
        validate_extension(folder, contract, ledger)
        return previous_dispatch(*args, **kwargs)

    def checkpoint(destination, payload):
        if Path(destination) == folder:
            payload = {**payload, "runtime_policy_extension": binding}
        return previous_checkpoint(destination, payload)

    def finish(destination, run_contract, run_ledger, points, *args, **kwargs):
        validate_extension(folder, contract, ledger)
        return previous_finish(destination, run_contract, run_ledger, points, *args,
                               runtime_policy_extension=binding, **kwargs)

    old.private_mpi_runtime_witness, candidate.dispatch = witness, dispatch
    old.checkpoint, old.finish = checkpoint, finish
    try:
        yield binding
    finally:
        old.private_mpi_runtime_witness, candidate.dispatch = previous_witness, previous_dispatch
        old.checkpoint, old.finish = previous_checkpoint, previous_finish


def replay_paid_point(folder, contract, prepared, ledger):
    """No native launch, no reservation, no analysis event: replay paid +120."""
    attempt = next(row for row in ledger.reservations()
                   if row["point_id"] == "point-p120" and row["kind"] == "saddles")
    metadata = attempt["metadata"]
    inp = Path(metadata["input_path"])
    if (inp != folder / "point-p120/optimization/opt.inp"
            or metadata.get("from_offset_deg") != 60 or metadata.get("offset_deg") != 120):
        raise ValueError("Retained +120 mapping/seed differs")
    before = ledger.evidence()
    with runtime_reader_adapter(folder, contract, ledger):
        node = candidate.analyze_candidate(inp, metadata["seed_A"], 120, prepared["origin"], contract)
    if ledger.evidence() != before:
        raise ValueError("Zero-job replay unexpectedly changed the journal")
    receipt = {"schema_version": 1, "kind": "OH_CH4_paid_plus120_runtime_replay",
               "new_native_jobs": 0, "runtime_policy_extension": pilot.artifact(folder / EXTENSION),
               "point": node, "rate_claim_allowed": False}
    path = folder / "runtime-policy-plus120-replay.json"
    if path.exists():
        if not old.continuation.replay_equal(pilot.read_json(path), receipt):
            raise ValueError("Retained runtime replay receipt differs")
    else:
        old.atomic_save(path, receipt)
    return node


def run(source, *, preflight_only=False):
    folder = old.owned_run(source)
    contract, prepared = candidate.reopen_campaign(folder)
    ledger = old.ledger_for(folder, contract)
    ledger.initialize()
    prepare_extension(folder, contract, prepared, ledger)
    with prepared["parent_ledger"].controller(), ledger.controller():
        node = replay_paid_point(folder, contract, prepared, ledger)
    print("REUSED +120: native and physical candidate checks passed; zero new jobs; "
          "loader trace is diagnostic-only", flush=True)
    print("ATTEMPTS " + str(ledger.counts()) + "; unchanged caps " + str(contract["caps"]), flush=True)
    if preflight_only:
        return node
    print("NEXT NEW TARGET +170; four cores, six hours per job; no rate claim", flush=True)
    with runtime_reader_adapter(folder, contract, ledger):
        return candidate.run(folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--preflight-only", action="store_true", help="Replay the paid point without native dispatch")
    arguments = parser.parse_args()
    run(arguments.source, preflight_only=arguments.preflight_only)
