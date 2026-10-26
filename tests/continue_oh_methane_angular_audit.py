"""Append-only continuation of the stopped OH/CH4 three-point angular pilot.

Reuse its verified center derivatives and native-converged -15 degree candidate.
Only candidate acquisition changes; the frozen physical derivative gates do not.
No rotor partition, rate correction, or production admission is implemented here.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path
import re
import sys
import tempfile
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np

from storca.ts_candidate_eligibility import orca_optimization_candidate
from storca.ts_angular_audit import rotate_dihedral_seed, signed_dihedral, wrap_angle

spec = importlib.util.spec_from_file_location("frozen_angular_pilot", Path(__file__).with_name("run_oh_methane_angular_audit.py"))
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)
state = pilot.state
PHYSICAL_KEYS = ("angular_tolerance_rad", "tangent_gradient_RMS_Eh_bohr", "tangent_gradient_max_Eh_bohr",
                 "negative_tangent_eigenvalues", "positive_tangent_eigenvalues", "minimum_transfer_mode_overlap",
                 "maximum_S_squared_change", "maximum_orthogonal_ZPE_shift_J_mol", "maximum_orthogonal_partition_log_shift")
LOCAL_CAPS = {"saddles": 1, "NumFreq": 2}
PRIOR_ATTEMPTS = {"saddles": 4, "NumFreq": 1}


def completed_job(inp):
    """Reopen owned input/output/lease, never trust a detached success flag."""
    execution = pilot.read_json(inp.parent / "execution.json")
    if execution["failure"] is not None or pilot.reopen(execution["input"]) != inp:
        raise ValueError("Retained job execution differs")
    paths = {pilot.reopen(record) for record in execution["artifacts"]}
    if not {inp.with_suffix(".out"), inp.with_suffix(".process.json")}.issubset(paths):
        raise ValueError("Retained output/lease binding missing")
    lease = pilot.read_json(inp.with_suffix(".process.json"))
    if (lease.get("status") != "terminated" or lease.get("returncode") != 0
            or lease.get("engine") != "orca" or lease.get("ncores") != 4
            or lease.get("input") != str(inp) or lease.get("output") != str(inp.with_suffix(".out"))
            or lease.get("input_sha256") != pilot.artifact(inp)["sha256"]
            or inp.with_suffix(".out").read_text().count("ORCA TERMINATED NORMALLY") != 1):
        raise ValueError("Owned normal native completion required")
    return execution


def candidate_receipt(inp, seed, policy):
    completed_job(inp)
    text = inp.with_suffix(".out").read_text()
    receipt = orca_optimization_candidate(text, max_step_multiple=2.0)
    primitive = r"D\((?:H\s+1,O\s+0,C\s+2,H\s+4|H\s+4,C\s+2,O\s+0,H\s+1)\)[^\n]*\bC\s*$"
    if not re.search(primitive, text, re.M):
        raise ValueError("Native mapped constrained primitive missing")
    labels, xyz, _ = pilot.read_xyz(inp.with_suffix(".xyz"))
    if tuple(labels) != pilot.LABELS or xyz.shape != (7, 3) or not np.isfinite(xyz).all():
        raise ValueError("Candidate mapped geometry differs")
    # Native final Cartesian coordinates have six decimals, unlike the XYZ.
    blocks = re.findall(r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\s*\n(.*?)(?:\n\s*\n)", text, re.S)
    fields = [row.split() for row in blocks[-1].splitlines()] if blocks else []
    if (len(fields) != 7 or tuple(row[0] for row in fields) != pilot.LABELS
            or any(len(row) != 4 for row in fields)
            or not np.allclose([[float(v) for v in row[1:]] for row in fields], xyz, atol=6e-7, rtol=0)):
        raise ValueError("Candidate XYZ does not match final native frame")
    gap = abs(wrap_angle(signed_dihedral(xyz, pilot.DIHEDRAL)-signed_dihedral(seed, pilot.DIHEDRAL)))
    if gap > policy["angular_tolerance_rad"]:
        raise ValueError("Candidate constrained mapped angle differs")
    receipt.update({"input": pilot.artifact(inp), "output": pilot.artifact(inp.with_suffix(".out")),
                    "xyz": pilot.artifact(inp.with_suffix(".xyz")), "native_state": pilot.native_state(inp),
                    "constraint_gap_rad": gap})
    return xyz, receipt


def replay_equal(actual, expected):
    """Tight numerical replay equality; booleans/identities remain exact."""
    if isinstance(actual, dict) and isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(replay_equal(actual[k], expected[k]) for k in actual)
    if isinstance(actual, (list, tuple)) and isinstance(expected, (list, tuple)):
        return len(actual) == len(expected) and all(replay_equal(a, b) for a, b in zip(actual, expected))
    if isinstance(actual, float) and isinstance(expected, float):
        return math.isfinite(actual) and math.isfinite(expected) and math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12)
    return type(actual) is type(expected) and actual == expected


def claim_path(source):
    record = pilot.artifact(source / "contract.json")
    return ROOT / "runs" / ("oh-methane-angular-continuation-claim-" + record["sha256"] + ".json")


def load_predecessor(source):
    source = source.resolve(strict=True)
    if source.parent != (ROOT / "runs").resolve():
        raise ValueError("Retained repository runs child required")
    if claim_path(source).exists():
        raise ValueError("This source already has a continuation; inspect its claim instead of redispatching")
    result = pilot.read_json(source / "result.json")
    if (result.get("kind") != "bounded_OH_CH4_TS_angular_saddle_pilot" or result.get("status") != "unresolved"
            or result.get("failure_reason") != "ValueError: Native constrained convergence candidate gate unresolved"
            or result.get("old_sources_unchanged") is not True or result.get("rate_claim_allowed") is not False
            or result.get("whole_increment_attempts") != PRIOR_ATTEMPTS
            or result.get("attempts") != {"saddles": 1, "NumFreq": 1}):
        raise ValueError("Exact stopped candidate-eligibility predecessor required")
    contract = pilot.read_json(pilot.reopen(result["contract"]))
    state.reopen_inventory(contract["source_inventory"])
    for record in contract["source_snapshots"]:
        if pilot.reopen(record["source"]).read_bytes() != pilot.reopen(record["snapshot"]).read_bytes():
            raise ValueError("Frozen predecessor implementation snapshot differs")
    origin = pilot.load_source(Path(contract["source_replay"]))
    center = result["points"][0]
    if center["offset_deg"] != 0.0 or center["status"] != "passed_pilot_point":
        raise ValueError("Completed center point required")
    center_opt = pilot.reopen(center["optimization"]["input"])
    if center_opt.read_text() != pilot.input_text(origin["xyz_A"], "opt", origin["hessian"]):
        raise ValueError("Center acquisition input differs")
    center_xyz, center_candidate = candidate_receipt(center_opt, origin["xyz_A"], contract["policy"])
    center_freq = pilot.reopen(center["frequency"]["input"])
    completed_job(center_freq)
    if center_freq.read_text() != pilot.input_text(center_xyz, "freq"):
        raise ValueError("Center derivative input differs")
    replay = pilot.frequency_receipt(center_freq, center_xyz, origin, contract["policy"])
    if not replay_equal(replay, center["frequency"]) or replay["pilot_point_gate_passed"] is not True:
        raise ValueError("Retained center physical derivative replay differs")
    minus = source / "point-01/optimization/opt.inp"
    seed = rotate_dihedral_seed(origin["xyz_A"] / pilot.BOHR_A, pilot.DIHEDRAL, math.radians(-15.0))*pilot.BOHR_A
    if minus.read_text() != pilot.input_text(seed, "opt", origin["hessian"]):
        raise ValueError("Minus acquisition input differs from exact retained requested seed")
    xyz, candidate = candidate_receipt(minus, seed, contract["policy"])
    paths = {Path(path) for path in contract["source_inventory"]}
    paths.update(state.bound_paths(contract))
    paths.update(state.bound_paths(result))
    paths.update({source / "contract.json", source / "result.json", source / "point-00/point-result.json",
                  minus, minus.with_suffix(".out"), minus.with_suffix(".xyz"), minus.with_suffix(".process.json"),
                  minus.parent / "execution.json", center_freq.parent / "execution.json",
                  center_freq.with_suffix(".process.json"), center_opt.parent / "execution.json"})
    return {"origin": origin, "contract": contract, "result": result, "paths": paths,
            "center": {**center, "candidate_eligibility": center_candidate, "derivatives_reused": True},
            "minus_xyz": xyz, "minus_candidate": candidate, "minus_seed": seed}


def build_contract(source, prepared, folder):
    import storca.ts_candidate_eligibility as eligibility
    paths = prepared["paths"] | {Path(__file__), Path(eligibility.__file__)}
    inventory = state.inventory(paths)
    snapshots = folder / "source-snapshots"
    snapshots.mkdir()
    saved = []
    for index, record in enumerate(inventory.values()):
        path = pilot.reopen(record)
        if path.suffix == ".py":
            copy = snapshots / f"{index}-{path.name}"
            with copy.open("xb") as stream:
                stream.write(path.read_bytes())
            saved.append({"source": record, "snapshot": pilot.artifact(copy)})
    old = prepared["contract"]
    return {"schema_version": 1, "kind": "bounded_OH_CH4_TS_angular_candidate_continuation",
            "predecessor": pilot.artifact(source / "result.json"), "predecessor_contract": pilot.artifact(source / "contract.json"),
            "source_inventory": inventory, "source_snapshots": saved,
            "prior_attempts": PRIOR_ATTEMPTS, "local_caps": LOCAL_CAPS, "caps": old["caps"],
            "prior_candidate_policy": {key: value for key, value in old["policy"].items() if key not in PHYSICAL_KEYS},
            "candidate_policy": {"gradient_fraction_of_requested_tolerance": 1.0,
                                 "max_step_tolerance_multiplier": 2.0, "native_final_cycle_markers_required": True,
                                 "independent_derivative_validation_required": True},
            "policy": {key: old["policy"][key] for key in PHYSICAL_KEYS},
            "physical_policy_unchanged": True,
            "policy_revision_reason": "The extra 5%-of-tolerance screen denied native-converged candidates before independent physical validation; remove only that acquisition requirement.",
            "reused_center_frequency": prepared["center"]["frequency"]["input"],
            "reused_minus_optimization": prepared["minus_candidate"]["input"],
            "cores": old["cores"], "per_job_seconds": old["per_job_seconds"], "minimum_disk_reserve_bytes": old["minimum_disk_reserve_bytes"],
            "automatic_retries": 0, "rate_claim_allowed": False, "lifetime_claim_allowed": False,
            "production_defaults_changed": False, "network_reinjection_allowed": False}


def dispatch(folder, stem, text, attempts, contract):
    if stem not in ("opt", "freq"):
        raise ValueError("Unknown continuation task")
    key = "saddles" if stem == "opt" else "NumFreq"
    cap_key = "whole_increment_saddle_attempts" if stem == "opt" else "whole_increment_NumFreq_jobs"
    if (attempts[key] >= contract["local_caps"][key]
            or attempts[key] + contract["prior_attempts"][key] >= contract["caps"][cap_key]):
        raise ValueError("Continuation/cumulative attempt cap reached")
    state.reopen_inventory(contract["source_inventory"])
    return pilot.execute(folder, stem, text, attempts)


def run(source):
    source = source.resolve(strict=True)
    prepared = load_predecessor(source)  # All reuse and policy checks before native dispatch.
    folder = Path(tempfile.mkdtemp(prefix="oh-methane-angular-continuation-20261012-", dir=ROOT / "runs"))
    # Exclusive claim prevents silent repeats, concurrent launches and budget resets.
    pilot.save(claim_path(source), {"source_contract": pilot.artifact(source / "contract.json"),
                                  "continuation_directory": str(folder), "automatic_retry_allowed": False})
    print("WORKDIR " + str(folder), flush=True)
    started = time.monotonic()
    attempts = {"saddles": 0, "NumFreq": 0}
    points, contract, failure, unchanged = [prepared["center"]], None, None, True
    try:
        contract = build_contract(source, prepared, folder)
        pilot.save(folder / "contract.json", contract)
        origin = prepared["origin"]
        for index, offset in ((1, -15.0), (2, 15.0)):
            point_folder = folder / f"point-{index:02d}"
            point_folder.mkdir()
            seed = (prepared["minus_seed"] if index == 1 else
                    rotate_dihedral_seed(origin["xyz_A"] / pilot.BOHR_A, pilot.DIHEDRAL, math.radians(offset))*pilot.BOHR_A)
            row = {"offset_deg": offset, "seed_A": seed.tolist(), "status": "started"}
            points.append(row)
            try:
                if index == 1:
                    xyz, candidate = prepared["minus_xyz"], prepared["minus_candidate"]
                    row["optimization_reused"] = True
                else:
                    opt = dispatch(point_folder / "optimization", "opt", pilot.input_text(seed, "opt", origin["hessian"]), attempts, contract)
                    xyz, candidate = candidate_receipt(opt, seed, contract["policy"])
                row["candidate_eligibility"] = candidate
                pilot.save(point_folder / "candidate-eligibility.json", candidate)
                freq = dispatch(point_folder / "numerical-Hessian", "freq", pilot.input_text(xyz, "freq"), attempts, contract)
                row["frequency"] = pilot.frequency_receipt(freq, xyz, origin, contract["policy"])
                row["status"] = "passed_pilot_point" if row["frequency"]["pilot_point_gate_passed"] else "unresolved"
                if row["status"] != "passed_pilot_point":
                    raise ValueError("Constraint-projected saddle/transfer/tangent-gradient physical gate unresolved")
            except BaseException as error:
                row["status"] = "unresolved"
                row["failure_reason"] = f"{type(error).__name__}: {error}"
                raise
            finally:
                pilot.save(point_folder / "point-result.json", row)
        comparisons = pilot.orthogonal_comparison(points, origin["temperature_grid_K"])
        pilot.save(folder / "orthogonal-comparison.json", comparisons)
        if any(abs(row["ZPE_shift_J_mol"]) > contract["policy"]["maximum_orthogonal_ZPE_shift_J_mol"]
               or row["maximum_abs_log_partition_shift"] > contract["policy"]["maximum_orthogonal_partition_log_shift"]
               for row in comparisons):
            raise ValueError("Orthogonal contributions vary materially; separable angular model unresolved")
    except BaseException as error:
        failure = f"{type(error).__name__}: {error}"
        print("CONTINUATION UNRESOLVED: " + failure, flush=True)
    if contract is not None:
        try:
            state.reopen_inventory(contract["source_inventory"])
        except (ValueError, OSError) as error:
            unchanged = False
            failure = failure or str(error)
    result = {"schema_version": 1, "kind": "bounded_OH_CH4_TS_angular_candidate_continuation",
              "status": "unresolved" if failure else "pilot_passed_profile_pending", "failure_reason": failure,
              "points": points, "attempts": attempts,
              "whole_increment_attempts": {key: PRIOR_ATTEMPTS[key]+value for key, value in attempts.items()},
              "elapsed_seconds": time.monotonic()-started, "old_sources_unchanged": unchanged,
              "contract": pilot.artifact(folder / "contract.json") if contract is not None else None,
              "full_periodic_domain_covered": False, "angular_partition_admitted": False, "mirror_basin_count_resolved": False,
              "verification_status": "diagnostic_not_admitted", "rate_claim_allowed": False, "lifetime_claim_allowed": False,
              "production_defaults_changed": False, "network_reinjection_allowed": False}
    pilot.save(folder / "result.json", result)
    print(result["status"] + ": " + str(folder / "result.json"), flush=True)
    if failure:
        raise RuntimeError(failure)
    return folder


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Stopped three-point angular pilot directory")
    run(parser.parse_args().source)
