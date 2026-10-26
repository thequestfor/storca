"""Reuse the completed +15 OptTS after an append-only loader-reader correction.

Exactly one pending NumFreq is permitted. Old code, logs, receipts and physical
gates stay frozen; no new saddle, automatic retry or rate admission is allowed.
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys
import tempfile
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import storca.mpi_runtime_evidence as mpi_evidence

spec = importlib.util.spec_from_file_location("frozen_matched_angular_restart", Path(__file__).with_name("continue_oh_methane_matched_mpi.py"))
frozen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frozen)
pilot, state, continuation = frozen.pilot, frozen.state, frozen.continuation
PRIOR_ATTEMPTS = {"saddles": 6, "NumFreq": 2}
LOCAL_CAPS = {"saddles": 0, "NumFreq": 1}
EXPECTED_FAILURE = "ValueError: Private MPI loading at actual gradient workers not proved"


def claim_path(source):
    return ROOT / "runs" / ("oh-methane-angular-frequency-claim-"
                            + pilot.artifact(source / "result.json")["sha256"] + ".json")


def preflight(source):
    source = source.resolve(strict=True)
    if source.parent != (ROOT / "runs").resolve():
        raise ValueError("Owned repository continuation required")
    if claim_path(source).exists():
        raise ValueError("This completed optimization already has a frequency continuation")
    result = pilot.read_json(source / "result.json")
    if (result.get("status") != "unresolved" or result.get("failure_reason") != EXPECTED_FAILURE
            or result.get("attempts") != {"saddles": 1, "NumFreq": 0}
            or result.get("whole_increment_attempts") != PRIOR_ATTEMPTS
            or result.get("old_sources_unchanged") is not True
            or any(result.get(key) is not False for key in ("rate_claim_allowed", "lifetime_claim_allowed", "network_reinjection_allowed"))):
        raise ValueError("Exact stopped runtime-reader predecessor required")
    old = pilot.read_json(pilot.reopen(result["contract"]))
    state.reopen_inventory(old["source_inventory"])
    state.reopen_inventory(old["runtime"]["inventory"])
    if old["local_caps"] != {"saddles": 1, "NumFreq": 1} or old["physical_policy_unchanged"] is not True:
        raise ValueError("Retained finite runtime-only restart policy required")
    diagnostic_file = pilot.reopen(old["source"])
    diagnostic_result = pilot.read_json(diagnostic_file)
    diagnostic_contract_file = pilot.reopen(diagnostic_result["contract"])
    diagnostic_contract = pilot.read_json(diagnostic_contract_file)
    stopped_file = pilot.reopen(diagnostic_contract["source"])
    stopped = pilot.read_json(stopped_file)
    stopped_contract_file = pilot.reopen(stopped["contract"])
    stopped_contract = pilot.read_json(stopped_contract_file)
    original_file = pilot.reopen(stopped_contract["predecessor_contract"])
    original = pilot.read_json(original_file)
    origin = pilot.load_source(Path(original["source_replay"]))
    if old["policy"] != stopped_contract["policy"] or old["candidate_policy"] != stopped_contract["candidate_policy"]:
        raise ValueError("Retained candidate/physical gates differ")
    if Path(pilot.find_orca()).resolve() != pilot.reopen(diagnostic_contract["orca"]):
        raise ValueError("Bound native ORCA differs")
    prefix = Path(old["runtime"]["prefix"]).resolve(strict=True)
    if prefix != (diagnostic_file.parent / "mpi-build/install").resolve(strict=True):
        raise ValueError("Previously witnessed private runtime identity differs")
    env, runtime = frozen.runtime_environment(prefix)
    if runtime["inventory"] != old["runtime"]["inventory"]:
        raise ValueError("Private runtime files changed")
    opt = source / "point-02/optimization/opt.inp"
    continuation.completed_job(opt)
    if opt.read_text() != pilot.input_text(old["restart_xyz_A"], "opt", pilot.reopen(old["initial_hessian"])):
        raise ValueError("Retained optimization input differs")
    witness = mpi_evidence.private_mpi_runtime_witness(opt.with_suffix(".out").read_text(), prefix)
    xyz, candidate = continuation.candidate_receipt(opt, old["original_plus_seed_A"], old["policy"])
    restart_gap = abs(continuation.wrap_angle(continuation.signed_dihedral(xyz, pilot.DIHEDRAL)
                                            - continuation.signed_dihedral(old["restart_xyz_A"], pilot.DIHEDRAL)))
    if restart_gap > old["policy"]["angular_tolerance_rad"]:
        raise ValueError("Completed geometry differs from saved restart target")
    if [row["offset_deg"] for row in result["points"]] != [0.0, -15.0, 15.0]:
        raise ValueError("Exact three-point predecessor required")
    points = []
    paths = {Path(path) for path in old["source_inventory"]}
    paths.update(state.bound_paths(result))
    paths.update({Path(__file__), Path(mpi_evidence.__file__), source / "result.json", source / "contract.json",
                  source / "point-02/point-result.json", opt, opt.with_suffix(".out"), opt.with_suffix(".xyz"),
                  opt.with_suffix(".process.json"), opt.parent / "execution.json", diagnostic_file,
                  diagnostic_contract_file, stopped_file, stopped_contract_file, original_file})
    for row in result["points"][:2]:
        inp = pilot.reopen(row["frequency"]["input"])
        continuation.completed_job(inp)
        point_xyz = pilot.inline_geometry(inp)
        if inp.read_text() != pilot.input_text(point_xyz, "freq"):
            raise ValueError("Retained frequency input differs")
        replay = pilot.frequency_receipt(inp, point_xyz, origin, old["policy"])
        if not continuation.replay_equal(replay, row["frequency"]) or replay["pilot_point_gate_passed"] is not True:
            raise ValueError("Retained center/minus physical evidence replay differs")
        paths.add(inp.parent / "execution.json")
        points.append({**row, "derivatives_reused": True})
    return {"source": source, "old": old, "origin": origin, "runtime": runtime, "environment": env,
            "paths": paths, "points": points, "xyz": xyz, "candidate": candidate,
            "restart_gap": restart_gap, "optimization": opt, "witness": witness, "prefix": prefix}


def dispatch_frequency(folder, xyz, attempts, contract):
    if (attempts != {"saddles": 0, "NumFreq": 0}
            or contract["local_caps"] != LOCAL_CAPS
            or contract["prior_attempts"]["NumFreq"] >= contract["caps"]["whole_increment_NumFreq_jobs"]):
        raise ValueError("One pending frequency attempt only; no saddle/retry dispatch")
    state.reopen_inventory(contract["source_inventory"])
    inp = pilot.execute(folder, "freq", pilot.input_text(xyz, "freq"), attempts)
    continuation.completed_job(inp)
    witness = mpi_evidence.private_mpi_runtime_witness(inp.with_suffix(".out").read_text(), Path(contract["runtime"]["prefix"]))
    pilot.save(folder / "runtime-witness.json", witness)
    state.reopen_inventory(contract["source_inventory"])
    return inp


def run(source):
    prepared = preflight(source)
    source = prepared["source"]
    folder = Path(tempfile.mkdtemp(prefix="oh-methane-angular-frequency-20261013-", dir=ROOT / "runs"))
    pilot.save(claim_path(source), {"source": pilot.artifact(source / "result.json"),
                                  "continuation_directory": str(folder), "automatic_retry_allowed": False})
    print("WORKDIR " + str(folder), flush=True)
    old = prepared["old"]
    contract = {"schema_version": 1, "kind": "bounded_OH_CH4_angular_frequency_reader_continuation",
                "source": pilot.artifact(source / "result.json"), "source_inventory": state.inventory(prepared["paths"]),
                "runtime": prepared["runtime"], "selected_environment": prepared["environment"],
                "prior_attempts": PRIOR_ATTEMPTS, "local_caps": LOCAL_CAPS, "caps": old["caps"],
                "prior_EnGrad_execution_diagnostic_attempts": 2,
                "policy": old["policy"], "candidate_policy": old["candidate_policy"], "physical_policy_unchanged": True,
                "reused_optimization": pilot.artifact(prepared["optimization"]),
                "fixed_frequency_xyz_A": prepared["xyz"].tolist(),
                "reader_revision": {"source": pilot.artifact(Path(mpi_evidence.__file__)),
                                    "retained_failure": EXPECTED_FAILURE,
                                    "reason": "Stop path capture before the next intact dyld token; do not consume then discard an adjacent worker header",
                                    "physical_gate_changed": False, "old_artifacts_rewritten": False},
                "cores": 4, "per_job_seconds": pilot.JOB_SECONDS, "minimum_disk_reserve_bytes": pilot.MIN_FREE,
                "automatic_retries": 0, "rate_claim_allowed": False, "lifetime_claim_allowed": False,
                "network_reinjection_allowed": False, "production_defaults_changed": False}
    pilot.save(folder / "contract.json", contract)
    pilot.save(folder / "replayed-optimization-runtime.json", prepared["witness"])
    pilot.save(folder / "replayed-candidate-eligibility.json", prepared["candidate"])
    attempts = {"saddles": 0, "NumFreq": 0}
    points, failure, unchanged = prepared["points"], None, True
    row = {"offset_deg": 15.0, "status": "started", "optimization_reused": True,
           "candidate_eligibility": prepared["candidate"], "restart_constraint_gap_rad": prepared["restart_gap"]}
    points.append(row)
    point_folder = folder / "point-02"
    point_folder.mkdir()
    started = time.monotonic()
    try:
        with frozen.scoped_environment(prepared["environment"]):
            inp = dispatch_frequency(point_folder / "numerical-Hessian", prepared["xyz"], attempts, contract)
            row["frequency"] = pilot.frequency_receipt(inp, prepared["xyz"], prepared["origin"], contract["policy"])
            if row["frequency"]["pilot_point_gate_passed"] is not True:
                raise ValueError("Unchanged constrained saddle/transfer/stationarity gate unresolved")
            row["status"] = "passed_pilot_point"
        comparisons = pilot.orthogonal_comparison(points, prepared["origin"]["temperature_grid_K"])
        pilot.save(folder / "orthogonal-comparison.json", comparisons)
        if any(abs(value["ZPE_shift_J_mol"]) > contract["policy"]["maximum_orthogonal_ZPE_shift_J_mol"]
               or value["maximum_abs_log_partition_shift"] > contract["policy"]["maximum_orthogonal_partition_log_shift"] for value in comparisons):
            raise ValueError("Orthogonal coupling diagnostic unresolved; separable angular model not admitted")
    except BaseException as error:
        failure = f"{type(error).__name__}: {error}"
        if row["status"] != "passed_pilot_point":
            row["status"] = "unresolved"
            row["failure_reason"] = failure
        print("FREQUENCY CONTINUATION UNRESOLVED: " + failure, flush=True)
    finally:
        try:
            state.reopen_inventory(contract["source_inventory"])
        except (OSError, ValueError) as error:
            unchanged = False
            failure = failure or str(error)
        pilot.save(point_folder / "point-result.json", row)
        result = {"schema_version": 1, "kind": contract["kind"],
                  "status": "unresolved" if failure else "pilot_passed_profile_pending", "failure_reason": failure,
                  "contract": pilot.artifact(folder / "contract.json"), "points": points, "attempts": attempts,
                  "whole_increment_attempts": {key: PRIOR_ATTEMPTS[key]+value for key, value in attempts.items()},
                  "prior_EnGrad_execution_diagnostic_attempts": 2, "elapsed_seconds": time.monotonic()-started,
                  "old_sources_unchanged": unchanged, "full_periodic_domain_covered": False,
                  "angular_partition_admitted": False, "mirror_basin_count_resolved": False,
                  "verification_status": "diagnostic_not_admitted", "rate_claim_allowed": False,
                  "lifetime_claim_allowed": False, "production_defaults_changed": False, "network_reinjection_allowed": False}
        pilot.save(folder / "result.json", result)
        print(result["status"] + ": " + str(folder / "result.json"), flush=True)
    if failure:
        raise RuntimeError(failure)
    return folder


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Completed OptTS stopped only at runtime-evidence parsing")
    run(parser.parse_args().source)
