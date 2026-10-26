"""One append-only geometry restart and fresh NumFreq, using witnessed MPI4.1.1.

Reuse the frozen pilot's inputs, owned runner, candidate and physical readers.
Do not modify the failed job, production defaults, rate model or angular gates.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import importlib.util
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


continuation = load("matched_frozen_continuation", "continue_oh_methane_angular_audit.py")
diagnostic = load("matched_frozen_gradient_diagnostic", "run_oh_methane_gradient_diagnostic.py")
pilot, state = continuation.pilot, continuation.state
PRIOR_ATTEMPTS = {"saddles": 5, "NumFreq": 2}
LOCAL_CAPS = {"saddles": 1, "NumFreq": 1}


def claim_path(source):
    return ROOT / "runs" / ("oh-methane-angular-matched-mpi-claim-"
                            + pilot.artifact(source / "result.json")["sha256"] + ".json")


def runtime_witness(native, prefix):
    """Parse complete dyld tokens even when parallel stderr interleaves lines."""
    entries = re.findall(r"dyld\[(\d+)\]:\s*<[0-9A-Fa-f-]{36}>\s+(/[^\s<>]+)", native)
    entries = [(pid, path.split("dyld[", 1)[0]) for pid, path in entries]
    workers = {pid for pid, path in entries if Path(path).name == "orca_scfgrad_mpi"}
    mpi = [(pid, path) for pid, path in entries if Path(path).name == "libmpi.40.dylib"]
    expected = (prefix / "lib/libmpi.40.dylib").resolve(strict=True)
    if (len(workers) < 4 or not mpi or not workers.issubset({pid for pid, _ in mpi})
            or any(Path(path).resolve(strict=True) != expected for _, path in mpi)):
        raise ValueError("Private MPI loading at actual gradient workers not proved")
    return {"gradient_worker_PIDs": sorted(workers), "MPI_loading_witnesses": len(mpi),
            "loaded_MPI_library": str(expected), "all_identified_gradient_workers_witnessed": True,
            "scope": "Complete native loader tokens; not a proof of original crash causation"}


def runtime_environment(prefix):
    env, receipt = diagnostic.matched_runtime(prefix)
    # Explicit values survive a protected shell stripping inherited DYLD_*.
    env.update({"DYLD_LIBRARY_PATH": str(prefix / "lib"),
                "DYLD_FALLBACK_LIBRARY_PATH": str(prefix / "lib"),
                "DYLD_PRINT_LIBRARIES": "1"})
    env["OMPI_MCA_mca_base_env_list"] = ";".join(f"{name}={env[name]}" for name in (
        "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH", "DYLD_PRINT_LIBRARIES"))
    receipt["selected_environment"] = env
    paths = {Path(path) for path in receipt["inventory"]}
    paths.add(prefix / "bin/orted")
    paths.update((prefix / "lib").rglob("*.so"))
    paths.update(path for path in (prefix / "etc").rglob("*") if path.is_file())
    receipt["inventory"] = state.inventory(paths)
    return env, receipt


@contextmanager
def scoped_environment(values):
    """Process-local environment only; restored even after native failure."""
    old = {key: os.environ.get(key) for key in values}
    os.environ.update(values)
    try:
        yield
    finally:
        for key, value in old.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def preflight(source):
    source = source.resolve(strict=True)
    if source.parent != (ROOT / "runs").resolve():
        raise ValueError("Retained repository diagnostic required")
    if claim_path(source).exists():
        raise ValueError("Diagnostic already has a matched-MPI continuation; no redispatch")
    result = pilot.read_json(source / "result.json")
    if (result.get("status") != "completed_comparison" or result.get("old_sources_unchanged") is not True
            or result.get("native_attempts") != {"EnGrad": 2}
            or any(result.get(key) is not False for key in ("rate_claim_allowed", "lifetime_claim_allowed", "saddle_claim_allowed"))):
        raise ValueError("Closed unadmitted gradient comparison required")
    contract = pilot.read_json(pilot.reopen(result["contract"]))
    state.reopen_inventory(contract["source_inventory"])
    prefix = (source / "mpi-build/install").resolve(strict=True)
    env, runtime = runtime_environment(prefix)
    paths = {Path(path) for path in contract["source_inventory"]}
    paths.update({source / "result.json", source / "contract.json", Path(__file__),
                  Path(continuation.__file__), Path(diagnostic.__file__)})
    paths.update(Path(path) for path in runtime["inventory"])
    xyz = np.asarray(contract["xyz_A"])
    pair = []
    for case, record in zip(diagnostic.CASES, result["cases"], strict=True):
        case_file = pilot.reopen(record)
        if case_file != source / case / "case-result.json":
            raise ValueError("Diagnostic case identity differs")
        saved = pilot.read_json(case_file)
        inp = pilot.reopen(saved["input"])
        if (inp != source / case / "gradient.inp" or saved["status"] != "completed_execution_diagnostic"
                or inp.read_text() != diagnostic.input_text(xyz, 1 if case == "serial" else 4)):
            raise ValueError("Diagnostic input/completion differs")
        for artifact in [saved["runtime"], *saved["artifacts"]]:
            paths.add(pilot.reopen(artifact))
        paths.update({case_file, inp})
        lease = pilot.read_json(inp.with_suffix(".process.json"))
        native = inp.with_suffix(".out").read_text()
        if (lease["status"] != "terminated" or lease["returncode"] != 0
                or lease["input_sha256"] != pilot.artifact(inp)["sha256"]
                or native.count("ORCA TERMINATED NORMALLY") != 1):
            raise ValueError("Owned diagnostic completion missing")
        energy, gradient = pilot.read_engrad(inp.with_suffix(".engrad"), xyz / pilot.BOHR_A)
        printed_E, printed_g = diagnostic.printed_derivatives(native)
        if (abs(energy-printed_E) > 1e-12 or np.max(abs(gradient-printed_g)) > 5.1e-10
                or not np.allclose(diagnostic.native_frame(native), xyz, atol=6e-7, rtol=0)
                or energy != saved["receipt"]["energy_Eh"]
                or not np.array_equal(gradient, saved["receipt"]["gradient_Eh_bohr"])):
            raise ValueError("Raw diagnostic derivatives differ")
        pilot.native_state(inp)
        if case == "matched-mpi":
            retained_runtime = pilot.read_json(pilot.reopen(saved["runtime"]))
            if Path(retained_runtime["prefix"]).resolve() != prefix:
                raise ValueError("Previously witnessed runtime prefix differs")
            state.reopen_inventory(retained_runtime["inventory"])
            runtime_witness(native, prefix)
        pair.append((energy, gradient))
    if (abs(pair[0][0]-pair[1][0]) > contract["agreement_tolerances"]["energy_Eh"]
            or np.max(abs(pair[0][1]-pair[1][1])) > contract["agreement_tolerances"]["max_gradient_component_Eh_bohr"]):
        raise ValueError("Reparsed diagnostic comparison does not pass")
    stopped_file = pilot.reopen(contract["source"])
    stopped = pilot.read_json(stopped_file)
    previous = pilot.read_json(pilot.reopen(stopped["contract"]))
    state.reopen_inventory(previous["source_inventory"])
    if (stopped.get("status") != "unresolved" or stopped.get("old_sources_unchanged") is not True
            or stopped.get("whole_increment_attempts") != PRIOR_ATTEMPTS
            or [row["offset_deg"] for row in stopped["points"]] != [0.0, -15.0, 15.0]):
        raise ValueError("Exact stopped three-point predecessor required")
    original = pilot.read_json(pilot.reopen(previous["predecessor_contract"]))
    origin = pilot.load_source(Path(original["source_replay"]))
    if Path(pilot.find_orca()).resolve() != pilot.reopen(contract["orca"]):
        raise ValueError("Bound native ORCA differs")
    orca_folder = pilot.reopen(contract["orca"]).parent
    paths.update(orca_folder / name for name in ("orca_numfreq", "orca_startup_mpi", "orca_prop_mpi"))
    paths.update((orca_folder / "lib").glob("liborca*.dylib"))
    paths.update(state.bound_paths(stopped))
    paths.update({stopped_file, pilot.reopen(stopped["contract"]), pilot.reopen(previous["predecessor_contract"])})
    points = []
    for row in stopped["points"][:2]:
        inp = pilot.reopen(row["frequency"]["input"])
        continuation.completed_job(inp)
        point_xyz = pilot.inline_geometry(inp)
        if inp.read_text() != pilot.input_text(point_xyz, "freq"):
            raise ValueError("Retained derivative method differs")
        replay = pilot.frequency_receipt(inp, point_xyz, origin, previous["policy"])
        if not continuation.replay_equal(replay, row["frequency"]) or replay["pilot_point_gate_passed"] is not True:
            raise ValueError("Retained center/minus physical receipt replay differs")
        paths.add(inp.parent / "execution.json")
        points.append({**row, "derivatives_reused": True})
    failed = stopped_file.parent / "point-02/optimization/opt.inp"
    failed_native = failed.with_suffix(".out").read_text()
    if ("error termination in SCF gradient" not in failed_native
            or "ORCA TERMINATED NORMALLY" in failed_native
            or failed_native.count("GEOMETRY OPTIMIZATION CYCLE") != 6):
        raise ValueError("Exact interrupted cycle-six optimization required")
    saved_xyz_file = pilot.reopen(contract["coordinates_source"])
    labels, saved_xyz, _ = pilot.read_xyz(saved_xyz_file)
    if tuple(labels) != pilot.LABELS or not np.array_equal(saved_xyz, xyz):
        raise ValueError("Restart is not the exact tested full-precision geometry")
    seed = continuation.rotate_dihedral_seed(origin["xyz_A"] / pilot.BOHR_A, pilot.DIHEDRAL, math.radians(15))*pilot.BOHR_A
    if failed.read_text() != pilot.input_text(seed, "opt", origin["hessian"]):
        raise ValueError("Original plus-angle method/seed differs")
    angle_gap = abs(continuation.wrap_angle(continuation.signed_dihedral(xyz, pilot.DIHEDRAL)
                                          - continuation.signed_dihedral(seed, pilot.DIHEDRAL)))
    if angle_gap > previous["policy"]["angular_tolerance_rad"]:
        raise ValueError("Saved restart target differs from original plus-angle target")
    paths.update({failed, failed.with_suffix(".out"), saved_xyz_file, failed.with_suffix(".process.json"),
                  failed.parent / "execution.json", origin["hessian"]})
    return {"source": source, "origin": origin, "points": points, "xyz": xyz, "original_seed": seed,
            "angle_gap": angle_gap, "previous": previous, "paths": paths, "environment": env,
            "runtime": runtime, "prefix": prefix, "diagnostic_contract": contract}


def dispatch(folder, stem, text, attempts, contract):
    if stem not in ("opt", "freq"):
        raise ValueError("Only opt/freq continuation tasks")
    key = "saddles" if stem == "opt" else "NumFreq"
    global_key = "whole_increment_saddle_attempts" if stem == "opt" else "whole_increment_NumFreq_jobs"
    if (attempts[key] >= contract["local_caps"][key]
            or attempts[key] + contract["prior_attempts"][key] >= contract["caps"][global_key]):
        raise ValueError("Local/cumulative continuation attempt cap reached")
    state.reopen_inventory(contract["source_inventory"])
    inp = pilot.execute(folder, stem, text, attempts)
    continuation.completed_job(inp)
    witness = runtime_witness(inp.with_suffix(".out").read_text(), Path(contract["runtime"]["prefix"]))
    pilot.save(folder / "runtime-witness.json", witness)
    state.reopen_inventory(contract["source_inventory"])
    return inp


def run(source):
    prepared = preflight(source)
    source = prepared["source"]
    folder = Path(tempfile.mkdtemp(prefix="oh-methane-angular-matched-mpi-20261013-", dir=ROOT / "runs"))
    pilot.save(claim_path(source), {"source": pilot.artifact(source / "result.json"),
                                  "continuation_directory": str(folder), "automatic_retry_allowed": False})
    print("WORKDIR " + str(folder), flush=True)
    contract = {"schema_version": 1, "kind": "bounded_OH_CH4_angular_matched_MPI_geometry_restart",
                "source": pilot.artifact(source / "result.json"), "source_inventory": state.inventory(prepared["paths"]),
                "runtime": prepared["runtime"], "selected_environment": prepared["environment"],
                "prior_attempts": PRIOR_ATTEMPTS, "local_caps": LOCAL_CAPS, "caps": prepared["previous"]["caps"],
                "prior_EnGrad_execution_diagnostic_attempts": 2, "policy": prepared["previous"]["policy"],
                "candidate_policy": prepared["previous"]["candidate_policy"], "physical_policy_unchanged": True,
                "restart_xyz_A": prepared["xyz"].tolist(), "original_plus_seed_A": prepared["original_seed"].tolist(),
                "restart_target_gap_rad": prepared["angle_gap"], "initial_hessian": pilot.artifact(prepared["origin"]["hessian"]),
                "restart_scope": "Geometry-only fresh SCF; verified original initial Hessian reused as an optimizer estimate, not a checkpoint/point verification",
                "cores": 4, "per_job_seconds": pilot.JOB_SECONDS, "minimum_disk_reserve_bytes": pilot.MIN_FREE,
                "automatic_retries": 0, "rate_claim_allowed": False, "lifetime_claim_allowed": False,
                "network_reinjection_allowed": False, "production_defaults_changed": False}
    pilot.save(folder / "contract.json", contract)
    attempts = {"saddles": 0, "NumFreq": 0}
    points, failure, unchanged = prepared["points"], None, True
    row = {"offset_deg": 15.0, "status": "started", "geometry_only_restart": True}
    points.append(row)
    started = time.monotonic()
    point_folder = folder / "point-02"
    point_folder.mkdir()
    try:
        with scoped_environment(prepared["environment"]):
            opt = dispatch(point_folder / "optimization", "opt", pilot.input_text(prepared["xyz"], "opt", prepared["origin"]["hessian"]), attempts, contract)
            xyz, candidate = continuation.candidate_receipt(opt, prepared["original_seed"], contract["policy"])
            restart_gap = abs(continuation.wrap_angle(continuation.signed_dihedral(xyz, pilot.DIHEDRAL)
                                                    - continuation.signed_dihedral(prepared["xyz"], pilot.DIHEDRAL)))
            if restart_gap > contract["policy"]["angular_tolerance_rad"]:
                raise ValueError("Optimized angle differs from saved restart target")
            row["candidate_eligibility"] = candidate
            row["restart_constraint_gap_rad"] = restart_gap
            pilot.save(point_folder / "candidate-eligibility.json", candidate)
            freq = dispatch(point_folder / "numerical-Hessian", "freq", pilot.input_text(xyz, "freq"), attempts, contract)
            row["frequency"] = pilot.frequency_receipt(freq, xyz, prepared["origin"], contract["policy"])
            if row["frequency"]["pilot_point_gate_passed"] is not True:
                raise ValueError("Unchanged constrained saddle/transfer/stationarity gate unresolved")
            row["status"] = "passed_pilot_point"
        comparisons = pilot.orthogonal_comparison(points, prepared["origin"]["temperature_grid_K"])
        pilot.save(folder / "orthogonal-comparison.json", comparisons)
        if any(abs(value["ZPE_shift_J_mol"]) > contract["policy"]["maximum_orthogonal_ZPE_shift_J_mol"]
               or value["maximum_abs_log_partition_shift"] > contract["policy"]["maximum_orthogonal_partition_log_shift"] for value in comparisons):
            raise ValueError("Orthogonal coupling diagnostic unresolved; no separable angular admission")
    except BaseException as error:
        failure = f"{type(error).__name__}: {error}"
        if row["status"] != "passed_pilot_point":
            row["status"] = "unresolved"
            row["failure_reason"] = failure
        print("MATCHED-MPI CONTINUATION UNRESOLVED: " + failure, flush=True)
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
    parser.add_argument("source", type=Path, help="Completed serial/matched-MPI gradient diagnostic")
    run(parser.parse_args().source)
