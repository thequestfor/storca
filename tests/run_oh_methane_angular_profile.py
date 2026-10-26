"""Bounded OH/CH4 periodic constrained-profile continuation, not a rate model.

One source/resume argument. Completed native jobs are replayed, not repeated.
The unreviewed reduced Hamiltonian is deliberately diagnostic-only. Every
scientific decision and every attempted native calculation remains recorded.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np

from storca import ts_angular_profile as profile
from storca.angular_attempts import AttemptLedger
from storca.mpi_runtime_evidence import private_mpi_runtime_witness
from storca.ts_angular_audit import rotate_dihedral_seed, signed_dihedral, wrap_angle
from storca.ts_candidate_eligibility import orca_optimization_candidate


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


completion = load("periodic_frozen_frequency", Path(__file__).with_name("complete_oh_methane_angular_frequency.py"))
frozen, pilot, state, continuation = completion.frozen, completion.pilot, completion.state, completion.continuation
KIND = "bounded_OH_CH4_periodic_profile"
PRIOR = {"saddles": 6, "NumFreq": 3}
CAPS = {"saddles": 17, "NumFreq": 4}
CLAIM_FLAGS = ("rate_claim_allowed", "lifetime_claim_allowed", "network_reinjection_allowed", "production_defaults_changed")


class ScientificGateFailure(ValueError):
    """A physical diagnostic failed; do not silently retry or relax it."""


class AnalysisPending(ValueError):
    """Execution is retained but analysis needs attention; never redispatch."""


@contextmanager
def scientific_boundary():
    """Known applicability failures close; unknown reader failures stay replayable."""
    try:
        yield
    except ValueError as error:
        known_physical = (
            "Angular chart has a zero-length", "Angular chart is near-collinear",
            "Angular gradient/Hessian finite-step convergence unresolved",
            "Angular derivative violates rigid-motion invariance",
            "Independent angular constraint normal missing",
            "Native optimization candidate eligibility unresolved",
            "Native doublet spin continuity warning", "Electronic spin continuity warning",
        )
        if isinstance(error, ScientificGateFailure):
            raise
        if str(error).startswith(known_physical):
            raise ScientificGateFailure(str(error)) from error
        raise


def owned_run(path):
    path = Path(path).resolve(strict=True)
    if path.parent != (ROOT / "runs").resolve() or not path.is_dir():
        raise ValueError("Owned repository runs child required")
    return path


def assert_unadmitted(value):
    if any(value.get(key) is not False for key in CLAIM_FLAGS):
        raise ValueError("Unadmitted diagnostic evidence required")


def origin_from_contract(contract):
    """Follow only explicit source/contract bindings, not a runs-tree search."""
    paths, seen, current = set(), set(), contract
    for _ in range(12):
        if "source_replay" in current:
            return pilot.load_source(Path(current["source_replay"])), paths
        if "predecessor_contract" in current:
            path = pilot.reopen(current["predecessor_contract"])
        else:
            source = pilot.reopen(current["source"])
            paths.add(source)
            path = pilot.reopen(pilot.read_json(source)["contract"])
        if path in seen:
            raise ValueError("Cyclic predecessor chain")
        seen.add(path)
        paths.add(path)
        current = pilot.read_json(path)
        if "source_inventory" in current:
            state.reopen_inventory(current["source_inventory"])
    raise ValueError("Predecessor chain does not resolve to retained source")


def preflight(source):
    """Recompute all three local receipts with zero ORCA launches."""
    source = owned_run(source)
    result = pilot.read_json(source / "result.json")
    if (result.get("kind") != "bounded_OH_CH4_angular_frequency_reader_continuation"
            or result.get("status") != "pilot_passed_profile_pending"
            or result.get("attempts") != {"saddles": 0, "NumFreq": 1}
            or result.get("whole_increment_attempts") != PRIOR
            or result.get("old_sources_unchanged") is not True
            or result.get("angular_partition_admitted") is not False
            or result.get("mirror_basin_count_resolved") is not False):
        raise ValueError("Exact closed three-point frequency predecessor required")
    assert_unadmitted(result)
    contract_path = pilot.reopen(result["contract"])
    if contract_path != source / "contract.json":
        raise ValueError("Predecessor contract location differs")
    old = pilot.read_json(contract_path)
    assert_unadmitted(old)
    state.reopen_inventory(old["source_inventory"])
    state.reopen_inventory(old["runtime"]["inventory"])
    if (old["caps"]["whole_increment_saddle_attempts"] != CAPS["saddles"]
            or old["caps"]["whole_increment_NumFreq_jobs"] != CAPS["NumFreq"]
            or old["physical_policy_unchanged"] is not True):
        raise ValueError("Retained whole-campaign caps/physical policy differ")
    origin, paths = origin_from_contract(old)
    if Path(pilot.find_orca()).resolve() not in {Path(path) for path in old["source_inventory"]}:
        raise ValueError("Selected ORCA not bound by predecessor")
    prefix = Path(old["runtime"]["prefix"]).resolve(strict=True)
    environment, runtime = frozen.runtime_environment(prefix)
    if runtime["inventory"] != old["runtime"]["inventory"]:
        raise ValueError("Private runtime identity differs")
    if [row["offset_deg"] for row in result["points"]] != [0.0, -15.0, 15.0]:
        raise ValueError("Three retained mapped angular points required")
    points = []
    for row in result["points"]:
        freq = pilot.reopen(row["frequency"]["input"])
        continuation.completed_job(freq)
        xyz = pilot.inline_geometry(freq)
        if freq.read_text() != pilot.input_text(xyz, "freq"):
            raise ValueError("Retained frequency input differs")
        replay = pilot.frequency_receipt(freq, xyz, origin, old["policy"])
        if not continuation.replay_equal(replay, row["frequency"]) or replay["pilot_point_gate_passed"] is not True:
            raise ValueError("Retained constrained physical evidence differs")
        opt = pilot.reopen(row["candidate_eligibility"]["input"])
        seed = pilot.inline_geometry(opt)
        if opt.read_text() != pilot.input_text(seed, "opt", origin["hessian"]):
            raise ValueError("Retained acquisition input differs")
        candidate_xyz, candidate = continuation.candidate_receipt(opt, seed, old["policy"])
        if not np.allclose(candidate_xyz, xyz, atol=1e-9, rtol=0):
            raise ValueError("Retained optimization/frequency geometry differs")
        energy, gradient = pilot.read_engrad(freq.with_suffix(".engrad"), xyz / pilot.BOHR_A)
        diagnostics = profile.candidate_gradient_diagnostics(
            xyz / pilot.BOHR_A, origin["masses_amu"], gradient,
            target_angle_rad=signed_dihedral(seed, pilot.DIHEDRAL))
        if diagnostics["candidate_gate_passed"] is not True:
            raise ValueError("Retained gradient/coordinate replay gate differs")
        points.append({"offset_deg": row["offset_deg"], "xyz_A": xyz.tolist(),
                       "energy_Eh": energy, "gradient_Eh_bohr": gradient.tolist(),
                       "status": "derivative_verified_constrained", "frequency": replay,
                       "candidate_eligibility": candidate, "gradient_diagnostics": diagnostics,
                       "derivatives_reused": True,
                       "runtime_scope": "historical bound runtime; not relabelled as new private-MPI evidence"})
        paths.update({freq.parent / "execution.json", freq.with_suffix(".process.json"),
                      opt.parent / "execution.json", opt.with_suffix(".process.json")})
        if row["offset_deg"] == 15.0:
            for inp in (freq, opt):
                private_mpi_runtime_witness(inp.with_suffix(".out").read_text(), prefix)
    coupling = pilot.orthogonal_comparison(points, origin["temperature_grid_K"])
    require_coupling(coupling, old["policy"])
    paths.update(Path(path) for path in old["source_inventory"])
    paths.update(state.bound_paths(result))
    paths.update({source / "result.json", source / "contract.json"})
    return {"source": source, "old": old, "origin": origin, "points": points,
            "paths": paths, "runtime": runtime, "environment": environment}


def require_coupling(comparison, policy):
    if any(abs(row["ZPE_shift_J_mol"]) > policy["maximum_orthogonal_ZPE_shift_J_mol"]
           or row["maximum_abs_log_partition_shift"] > policy["maximum_orthogonal_partition_log_shift"]
           for row in comparison):
        raise ScientificGateFailure("Orthogonal coupling exceeds unchanged local-model gates")


def claim_path(source):
    return ROOT / "runs" / ("oh-methane-periodic-claim-" + pilot.artifact(source / "result.json")["sha256"] + ".json")


def make_campaign(prepared):
    """Exclusive source claim prevents parallel fresh directories/reset budgets."""
    source = prepared["source"]
    claim = claim_path(source)
    if claim.exists():
        recorded = pilot.read_json(claim)
        if recorded["source"] != pilot.artifact(source / "result.json"):
            raise ValueError("Source claim binding differs")
        return owned_run(recorded["continuation_directory"])
    folder = Path(tempfile.mkdtemp(prefix="oh-methane-periodic-20261015-", dir=ROOT / "runs"))
    # An interrupted creation stays claimed and requires inspection, not redispatch.
    pilot.save(claim, {"source": pilot.artifact(source / "result.json"),
                       "continuation_directory": str(folder), "automatic_retries": 0})
    paths = prepared["paths"] | {Path(__file__), Path(profile.__file__),
                                 ROOT / "storca/angular_attempts.py",
                                 Path(__file__).with_name("replay_oh_methane_angular_rates.py")}
    inventory = state.inventory(paths)
    snapshot_dir = folder / "source-snapshots"
    snapshot_dir.mkdir()
    snapshots = []
    for index, record in enumerate(inventory.values()):
        source_path = pilot.reopen(record)
        if source_path.suffix == ".py":
            dest = snapshot_dir / f"{index}-{source_path.name}"
            with dest.open("xb") as stream:
                stream.write(source_path.read_bytes())
            snapshots.append({"source": record, "snapshot": pilot.artifact(dest)})
    origin, old = prepared["origin"], prepared["old"]
    contract = {"schema_version": 1, "kind": KIND,
                "source": pilot.artifact(source / "result.json"), "source_inventory": inventory,
                "source_snapshots": snapshots, "runtime": prepared["runtime"],
                "selected_environment": prepared["environment"],
                "orca": pilot.artifact(pilot.find_orca()), "initial_hessian": pilot.artifact(origin["hessian"]),
                "baseline_recipe": pilot.artifact(Path(origin["state_contract"]["source_workdir"]) / "composite/composite-rate-recipe.json"),
                "previous_conditional_branches": pilot.artifact(origin["source"] / "result.json"),
                "prior_attempts": PRIOR, "caps": CAPS, "prior_EnGrad_execution_diagnostic_attempts": 2,
                "policy": old["policy"], "candidate_policy": old["candidate_policy"],
                "grid_offsets_deg": list(profile.GRID),
                "acquisition_schedule": [list(pair) for pair in profile.ACQUISITION_SCHEDULE],
                "temperature_grid_K": origin["temperature_grid_K"], "masses_amu": origin["masses_amu"].tolist(),
                "center_angle_rad": signed_dihedral(np.asarray(prepared["points"][0]["xyz_A"]), pilot.DIHEDRAL),
                "mapped_labels": list(pilot.LABELS), "dihedral_atom_indices": list(pilot.DIHEDRAL),
                "method": pilot.METHOD, "cores": 4, "maxcore_MB_per_rank": 1000,
                "per_job_seconds": pilot.JOB_SECONDS, "minimum_disk_reserve_bytes": pilot.MIN_FREE,
                "automatic_retries": 0, "physical_policy_unchanged": True,
                "closure_policy": {"proper_mapped_mass_weighted_RMSD_A": .02,
                                   "maximum_displacement_A": .05, "energy_gap_J_mol": 100.0},
                "counting_policy": {"domain": "full labelled signed 2pi chart",
                                    "fixed_reacting_atoms": [0, 1, 2, 6], "spectator_atoms": [3, 4, 5],
                                    "quotient_proved": False, "additional_optical_factor_applied": False,
                                    "additional_abstraction_factor_applied": False,
                                    "scope": "external symmetry, channel and optical counts unresolved jointly"},
                "partition_policy": {"mode": "diagnostic_only", "supported_quantitative_hamiltonian": False,
                                     "reason": "No reviewed angular/reaction/rotation measure and transverse-bath model; sparse Hessian coverage",
                                     "refinement_policy": "No quantitative-partition refinement without an admitted measure; preserve unused budget"},
                **{key: False for key in CLAIM_FLAGS}}
    pilot.save(folder / "contract.json", contract)
    return folder


def ledger_for(folder, contract):
    return AttemptLedger(folder / "attempts", prior=contract["prior_attempts"], caps=contract["caps"],
                         contract_sha256=pilot.artifact(folder / "contract.json")["sha256"])


def reopen_campaign(folder):
    folder = owned_run(folder)
    contract = pilot.read_json(folder / "contract.json")
    if contract.get("kind") != KIND or contract.get("schema_version") != 1:
        raise ValueError("Periodic-profile contract required")
    assert_unadmitted(contract)
    if (contract["prior_attempts"] != PRIOR or contract["caps"] != CAPS
            or contract["grid_offsets_deg"] != list(profile.GRID)
            or contract["acquisition_schedule"] != [list(pair) for pair in profile.ACQUISITION_SCHEDULE]
            or contract["partition_policy"]["mode"] != "diagnostic_only"
            or contract["partition_policy"]["supported_quantitative_hamiltonian"] is not False):
        raise ValueError("Frozen bounded diagnostic policy differs")
    if (contract["method"] != pilot.METHOD or contract["cores"] != 4
            or contract["maxcore_MB_per_rank"] != 1000
            or contract["per_job_seconds"] != pilot.JOB_SECONDS
            or contract["minimum_disk_reserve_bytes"] != pilot.MIN_FREE
            or contract["automatic_retries"] != 0
            or contract["mapped_labels"] != list(pilot.LABELS)
            or contract["dihedral_atom_indices"] != list(pilot.DIHEDRAL)
            or contract["selected_environment"] != contract["runtime"]["selected_environment"]):
        raise ValueError("Frozen method/mapping/resource/runtime contract differs")
    state.reopen_inventory(contract["source_inventory"])
    state.reopen_inventory(contract["runtime"]["inventory"])
    for record in contract["source_snapshots"]:
        if pilot.reopen(record["source"]).read_bytes() != pilot.reopen(record["snapshot"]).read_bytes():
            raise ValueError("Participating implementation snapshot differs")
    prepared = preflight(pilot.reopen(contract["source"]).parent)
    if contract["policy"] != prepared["old"]["policy"] or contract["runtime"]["inventory"] != prepared["runtime"]["inventory"]:
        raise ValueError("Reopened physical/runtime policy differs")
    return contract, prepared


def checkpoint(folder, payload):
    directory = folder / "checkpoints"
    directory.mkdir(exist_ok=True)
    index = len(list(directory.glob("*.json")))
    pilot.save(directory / f"{index:04d}.json", payload)


def reopen_bindings(payload):
    def visit(value):
        if isinstance(value, dict):
            if "path" in value and "sha256" in value:
                pilot.reopen(value)
            for child in value.values():
                visit(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                visit(child)
    visit(payload)


def freeze_rate_assessment(folder):
    consumer = load("periodic_angular_rate_consumer", Path(__file__).with_name("replay_oh_methane_angular_rates.py"))
    destination = folder / "angular-rate-assessment.json"
    if destination.exists():
        current = consumer.assess(folder)
        if not continuation.replay_equal(current, pilot.read_json(destination)):
            raise ValueError("Frozen angular rate assessment differs on replay")
    else:
        consumer.freeze_assessment(folder, destination)


def native_outcome(inp, expected_text):
    """Distinguish terminal native failure from a wrapper/analysis interruption."""
    if inp.read_text() != expected_text:
        raise AnalysisPending("Reserved native input differs; no redispatch")
    output, lease_path = inp.with_suffix(".out"), inp.with_suffix(".process.json")
    lease = pilot.read_json(lease_path)
    if (lease.get("status") != "terminated" or type(lease.get("returncode")) is not int
            or lease.get("engine") != "orca" or type(lease.get("ncores")) is not int or lease["ncores"] != 4
            or lease.get("input") != str(inp) or lease.get("output") != str(output)
            or lease.get("input_sha256") != pilot.artifact(inp)["sha256"]):
        raise AnalysisPending("Reservation has running or ambiguous native evidence; no redispatch")
    normal_markers = output.read_text().count("ORCA TERMINATED NORMALLY")
    return {"status": "normal" if lease["returncode"] == 0 and normal_markers == 1 else "failed",
            "returncode": lease["returncode"], "normal_termination_markers": normal_markers,
            "input": pilot.artifact(inp), "output": pilot.artifact(output), "lease": pilot.artifact(lease_path)}


def recover_execution(inp, expected_text):
    """Append an adoption receipt, never overwrite a wrapper's failure receipt."""
    outcome = native_outcome(inp, expected_text)
    if outcome["status"] != "normal":
        raise ScientificGateFailure("Owned native job did not complete normally")
    output, lease_path = inp.with_suffix(".out"), inp.with_suffix(".process.json")
    receipt = inp.parent / "execution.json"
    if not receipt.exists():
        atomic_save(receipt, {"input": pilot.artifact(inp), "failure": None,
                              "recovered_from_owned_terminal_lease": True,
                              "artifacts": [pilot.artifact(output), pilot.artifact(lease_path)],
                              "rate_claim_allowed": False, "lifetime_claim_allowed": False})
    try:
        execution = pilot.read_json(receipt)
    except (ValueError, UnicodeError):
        execution = None  # Preserve the malformed bytes; raw lease/output were validated above.
    if execution is not None:
        if pilot.reopen(execution["input"]) != inp:
            raise AnalysisPending("Retained execution input binding differs")
        reopen_bindings(execution)
        paths = {pilot.reopen(record) for record in execution["artifacts"]}
        if not {output, lease_path}.issubset(paths):
            raise AnalysisPending("Retained output/lease bindings missing")
    if execution is None or execution["failure"] is not None:
        adoption = inp.parent / "execution-adoption.json"
        payload = {"kind": "owned_normal_native_adoption", "original_execution": pilot.artifact(receipt),
                   "json_readable": execution is not None,
                   "native_outcome": outcome, "reason": "Native normal completion retained despite wrapper exception",
                   "old_execution_rewritten": False, "rate_claim_allowed": False}
        if adoption.exists():
            if not continuation.replay_equal(pilot.read_json(adoption), payload):
                raise AnalysisPending("Native adoption receipt differs")
        else:
            atomic_save(adoption, payload)
        return pilot.artifact(adoption)
    return pilot.artifact(receipt)


def owned_candidate_receipt(inp, seed, contract):
    """Same frozen candidate gates, using this controller's owned-job adoption.

    Historical receipts continue through the unmodified historical reader.
    This adapter is needed only because a new wrapper can stop after normal
    native completion; its failed execution receipt must not be overwritten.
    """
    recover_execution(inp, pilot.input_text(seed, "opt", pilot.reopen(contract["initial_hessian"])))
    text = inp.with_suffix(".out").read_text()
    if text.count("THE OPTIMIZATION HAS CONVERGED") == 0:
        raise ScientificGateFailure("Native optimization ended without a converged constrained candidate")
    candidate = orca_optimization_candidate(text, max_step_multiple=2.0)
    primitive = r"D\((?:H\s+1,O\s+0,C\s+2,H\s+4|H\s+4,C\s+2,O\s+0,H\s+1)\)[^\n]*\bC\s*$"
    if not re.search(primitive, text, re.M):
        raise ValueError("Native mapped constrained primitive missing")
    labels, xyz, _ = pilot.read_xyz(inp.with_suffix(".xyz"))
    if tuple(labels) != pilot.LABELS or xyz.shape != (7, 3) or not np.isfinite(xyz).all():
        raise ValueError("Candidate mapped geometry differs")
    blocks = re.findall(r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\s*\n(.*?)(?:\n\s*\n)", text, re.S)
    fields = [row.split() for row in blocks[-1].splitlines()] if blocks else []
    if (len(fields) != 7 or tuple(row[0] for row in fields) != pilot.LABELS
            or any(len(row) != 4 for row in fields)
            or not np.allclose([[float(v) for v in row[1:]] for row in fields], xyz, atol=6e-7, rtol=0)):
        raise ValueError("Candidate XYZ does not match final native frame")
    gap = abs(wrap_angle(signed_dihedral(xyz, pilot.DIHEDRAL) - signed_dihedral(seed, pilot.DIHEDRAL)))
    if gap > contract["policy"]["angular_tolerance_rad"]:
        raise ScientificGateFailure("Candidate constrained mapped angle differs")
    candidate.update({"input": pilot.artifact(inp), "output": pilot.artifact(inp.with_suffix(".out")),
                      "xyz": pilot.artifact(inp.with_suffix(".xyz")), "native_state": pilot.native_state(inp),
                      "constraint_gap_rad": gap})
    return xyz, candidate


def native_job(folder, point_id, stem, text, metadata, ledger, contract):
    """Reserve once; replay normally completed jobs after any analysis interruption."""
    job_folder = folder / point_id / ("optimization" if stem == "opt" else "numerical-Hessian")
    inp = job_folder / f"{stem}.inp"
    kind = "saddles" if stem == "opt" else "NumFreq"
    metadata = {**metadata, "input_path": str(inp), "input_sha256": hashlib.sha256(text.encode()).hexdigest()}
    existing = [r for r in ledger.reservations() if r["point_id"] == point_id and r["kind"] == kind]
    fresh = not existing
    if fresh:
        reservation = ledger.reserve(point_id, kind, metadata)
    else:
        if len(existing) != 1 or existing[0]["metadata"] != metadata:
            raise AnalysisPending("Reservation metadata differs; no redispatch")
        reservation = existing[0]
    attempt_id = reservation["id"]
    previous = ledger.execution(attempt_id)
    if previous is not None:
        reopen_bindings(previous)
        if previous["status"] != "normal":
            raise ScientificGateFailure(previous["failure_reason"])
    elif fresh:
        (folder / point_id).mkdir(exist_ok=True)
        state.reopen_inventory(contract["source_inventory"])
        state.reopen_inventory(contract["runtime"]["inventory"])
        try:
            pilot.execute(job_folder, stem, text, {"saddles": 0, "NumFreq": 0})
        except BaseException as error:
            try:
                outcome = native_outcome(inp, text)
            except (OSError, ValueError) as evidence_error:
                raise AnalysisPending(f"Executor interrupted; native state ambiguous ({evidence_error}); no redispatch") from error
            if outcome["status"] == "failed":
                reason = f"Native {stem} attempt failed: {type(error).__name__}: {error}; no automatic retry"
                artifacts = [pilot.artifact(path) for path in job_folder.glob("*")
                             if path.is_file() and path.suffix in (".json", ".inp", ".out")]
                ledger.seal_execution(attempt_id, {"status": "failed", "failure_reason": reason, "artifacts": artifacts})
                raise ScientificGateFailure(reason) from error
            # A wrapper failure is not native failure. Preserve its receipt and
            # adopt only the exact normally completed owned input/output/lease.
    # Includes recovery of a completed native job whose parent never sealed it.
    admitted_execution = recover_execution(inp, text)
    if previous is None:
        suffixes = (".inp", ".out", ".process.json", ".xyz", ".engrad") + ((".hess",) if stem == "freq" else ())
        ledger.seal_execution(attempt_id, {"status": "normal",
                                          "artifacts": [pilot.artifact(inp.with_suffix(ext)) for ext in suffixes if inp.with_suffix(ext).is_file()],
                                          "execution": admitted_execution})
    return inp, attempt_id


def topology_guard(xyz_A):
    rows = []
    for i, j, low, high in ((0, 1, .65, 1.25), (2, 3, .7, 1.35), (2, 4, .7, 1.35), (2, 5, .7, 1.35)):
        distance = float(np.linalg.norm(xyz_A[i] - xyz_A[j]))
        rows.append({"atoms": [i, j], "distance_A": distance, "guard_A": [low, high]})
        if not low <= distance <= high:
            raise ScientificGateFailure("Spectator topology acquisition guard failed; branch unresolved")
    return {"bonds": rows, "scope": "Mapped spectator acquisition guard, not chemical identity or saddle proof"}


def path_metric_diagnostics(points, masses_amu):
    """Illustrate path/constraint metric differences without admitting a measure."""
    ordered = sorted(points, key=lambda row: row["offset_deg"])
    angles = [math.radians(row["offset_deg"]) for row in ordered]
    frames = [np.asarray(row["xyz_A"]) / pilot.BOHR_A for row in ordered]
    extended_angles = [angles[-1] - 2 * math.pi, *angles, angles[0] + 2 * math.pi]
    extended_frames = [frames[-1], *frames, frames[0]]
    diagnostics = []
    for index, node in enumerate(ordered, 1):
        try:
            row = profile.relaxed_path_inertia(extended_angles, extended_frames, masses_amu, index)
            diagnostics.append({"offset_deg": node["offset_deg"], **row,
                                "step_refinement_verified": False,
                                "scope": "Coarse finite differences, reaction-normal/rotation coupling not resolved; not a partition measure"})
        except ValueError as error:
            diagnostics.append({"offset_deg": node["offset_deg"], "status": "unresolved",
                                "reason": str(error), "hamiltonian_admitted": False})
    return diagnostics


def analyze_candidate(inp, seed, offset, prepared, contract):
    xyz, candidate = owned_candidate_receipt(inp, seed, contract)
    energy, gradient = pilot.read_engrad(inp.with_suffix(".engrad"), xyz / pilot.BOHR_A)
    energies = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?[\d.]+)", inp.with_suffix(".out").read_text())
    if not energies or abs(float(energies[-1]) - energy) > 1e-9:
        raise ValueError("Final full EnGrad energy differs from native final energy")
    diagnostics = profile.candidate_gradient_diagnostics(
        xyz / pilot.BOHR_A, prepared["origin"]["masses_amu"], gradient,
        target_angle_rad=contract["center_angle_rad"] + math.radians(offset))
    native = candidate["native_state"]
    if max(abs(value - .758947) for value in native["S_squared"]) > contract["policy"]["maximum_S_squared_change"]:
        raise ScientificGateFailure("Native spin continuity gate unresolved")
    passed = (diagnostics["candidate_gate_passed"] is True
              and diagnostics["tangent_gradient_RMS_Eh_bohr"] <= contract["policy"]["tangent_gradient_RMS_Eh_bohr"]
              and diagnostics["tangent_gradient_max_Eh_bohr"] <= contract["policy"]["tangent_gradient_max_Eh_bohr"])
    if not passed:
        raise ScientificGateFailure("Final constrained tangent-gradient gate unresolved")
    topology = topology_guard(xyz)
    witness = private_mpi_runtime_witness(inp.with_suffix(".out").read_text(), Path(contract["runtime"]["prefix"]))
    return {"offset_deg": float(offset), "xyz_A": xyz.tolist(), "energy_Eh": energy,
            "gradient_Eh_bohr": gradient.tolist(), "gradient_diagnostics": diagnostics,
            "candidate_eligibility": candidate, "engrad": pilot.artifact(inp.with_suffix(".engrad")),
            "runtime_witness": witness, "spectator_topology": topology,
            "status": "constrained_optimizer_candidate", "saddle_character_verified_by_optimization_alone": False}


def seal_or_replay_analysis(ledger, attempt_id, value):
    saved = ledger.analysis(attempt_id)
    if saved is None:
        ledger.seal_analysis(attempt_id, value)
    elif not continuation.replay_equal(saved, value):
        raise AnalysisPending("Derived analysis differs from immutable receipt")


def atomic_save(path, value):
    """Publish a complete derived JSON exclusively; leave no partial final name."""
    content = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    descriptor, temporary = tempfile.mkstemp(prefix=".derived-json-", dir=path.parent)
    temporary = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path)  # Atomic no-overwrite publication on this filesystem.
        except FileExistsError:
            if path.read_bytes() != content:
                raise AnalysisPending("Previously published derived record differs: " + str(path))
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        temporary.unlink(missing_ok=True)


def publish_result(folder, result):
    """Bind exact proposed bytes before publication; resume can replay a seal-only stop."""
    path, binding_path = folder / "result.json", folder / "result-binding.json"
    content = (json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    binding = {"path": str(path), "sha256": hashlib.sha256(content).hexdigest()}
    if binding_path.exists():
        if pilot.read_json(binding_path) != binding:
            raise AnalysisPending("Previously sealed closeout differs from replay")
    else:
        atomic_save(binding_path, binding)
    if path.exists():
        if path.read_bytes() != content:
            raise AnalysisPending("Published result differs from its proposed closeout")
    else:
        atomic_save(path, result)
    pilot.reopen(binding)


def finish(folder, contract, ledger, points, failure=None, **details):
    state.reopen_inventory(contract["source_inventory"])
    state.reopen_inventory(contract["runtime"]["inventory"])
    result = {"schema_version": 1, "kind": KIND,
              "status": "profile_unresolved" if failure else "profile_complete_partition_unsupported",
              "failure_reason": failure, "contract": pilot.artifact(folder / "contract.json"),
              "source_replay": contract["source"], "points": points,
              "whole_increment_attempts": ledger.counts(), "attempt_reservations": ledger.reservations(),
              "attempt_evidence": ledger.evidence(),
              "prior_EnGrad_execution_diagnostic_attempts": 2,
              "coverage": profile.profile_coverage([p["offset_deg"] for p in points]),
              "full_periodic_domain_covered": False,
              "coverage_scope": "Discrete coarse candidates only; no continuous saddle/bath completeness proof",
              "angular_partition_admitted": False, "mirror_basin_count_resolved": False,
              "partition_assessment": profile.partition_admission({
                  "periodic_coverage": profile.profile_coverage([p["offset_deg"] for p in points])["coarse_grid_complete"],
                  "branch_closure": details.get("closure", {}).get("strict_labelled_closure_passed", False),
                  "candidate_gates_accepted": failure is None,
                  "selected_hessian_passed": any(abs(p["offset_deg"]) >= 60
                      and p.get("frequency", {}).get("pilot_point_gate_passed") is True for p in points),
                  "orthogonal_coupling_accepted": failure is None and "selected_orthogonal_comparison" in details}),
              "baseline_recipe": contract["baseline_recipe"],
              "previous_conditional_branches": contract["previous_conditional_branches"],
              "old_sources_unchanged": True, "verification_status": "diagnostic_not_admitted",
              "failed_diagnostic_does_not_establish_absent_path": True,
              "refinement_dispatched": False,
              "refinement_reason": "No admitted quantitative measure; no defensible partition convergence target",
              **{key: False for key in CLAIM_FLAGS}, **details}
    publish_result(folder, result)
    print(result["status"] + ": " + str(folder / "result.json"), flush=True)
    return result


def run(source):
    source = owned_run(source)
    maybe = pilot.read_json(source / "contract.json")
    folder = source if maybe.get("kind") == KIND else make_campaign(preflight(source))
    contract, prepared = reopen_campaign(folder)
    ledger = ledger_for(folder, contract)
    ledger.initialize()
    print("WORKDIR " + str(folder), flush=True)
    with ledger.controller():
        if (folder / "result.json").exists():
            pilot.reopen(pilot.read_json(folder / "result-binding.json"))
            result = pilot.read_json(folder / "result.json")
            if (result["contract"] != pilot.artifact(folder / "contract.json")
                    or result["whole_increment_attempts"] != ledger.counts()
                    or result["attempt_evidence"] != ledger.evidence()):
                raise ValueError("Closed periodic result/ledger binding differs")
            reopen_bindings(result)
            freeze_rate_assessment(folder)
            print("REUSED CLOSED RESULT " + result["status"], flush=True)
            return folder
        points = prepared["points"].copy()
        by_offset = {row["offset_deg"]: row for row in points}
        details = {}
        try:
            with scientific_boundary(), frozen.scoped_environment(contract["selected_environment"]):
                for previous_offset, offset in profile.ACQUISITION_SCHEDULE:
                    prior_xyz = np.asarray(by_offset[previous_offset]["xyz_A"])
                    desired = contract["center_angle_rad"] + math.radians(offset)
                    delta = wrap_angle(desired - signed_dihedral(prior_xyz, pilot.DIHEDRAL))
                    seed = rotate_dihedral_seed(prior_xyz / pilot.BOHR_A, pilot.DIHEDRAL, delta) * pilot.BOHR_A
                    point_id = f"point-{'p' if offset >= 0 else 'm'}{abs(offset):03d}"
                    inp, attempt_id = native_job(folder, point_id, "opt", pilot.input_text(seed, "opt", pilot.reopen(contract["initial_hessian"])),
                                                 {"offset_deg": offset, "from_offset_deg": previous_offset,
                                                  "seed_A": seed.tolist()}, ledger, contract)
                    row = analyze_candidate(inp, seed, offset, prepared, contract)
                    seal_or_replay_analysis(ledger, attempt_id, row)
                    if offset == 190:
                        details["closure_control"] = row
                        reference = by_offset[-170]
                        closure = profile.branch_closure(
                            reference["xyz_A"], row["xyz_A"], contract["masses_amu"],
                            reference["energy_Eh"], row["energy_Eh"], True, True)
                        details["closure"] = closure
                        details["spectator_permutation_diagnostics"] = profile.spectator_permutation_diagnostics(
                            reference["xyz_A"], row["xyz_A"], contract["masses_amu"])
                        if not closure["strict_labelled_closure_passed"]:
                            raise ScientificGateFailure("Independent periodic branch closure unresolved")
                    else:
                        points.append(row)
                        by_offset[offset] = row
                    checkpoint(folder, {"stage": "coarse_acquisition", "completed_offset_deg": offset,
                                        "whole_increment_attempts": ledger.counts(), "rate_claim_allowed": False})
                details["relaxed_path_metric_diagnostics"] = path_metric_diagnostics(points, contract["masses_amu"])
                selection = profile.select_final_hessian(points, contract["temperature_grid_K"], contract["masses_amu"])
                details["hessian_selection"] = selection
                if selection["selected_offset_deg"] is None:
                    raise ScientificGateFailure("Final-Hessian diagnostic selection unresolved; no safe allocation made")
                selected = by_offset[selection["selected_offset_deg"]]
                offset = selected["offset_deg"]
                point_id = f"point-{'p' if offset >= 0 else 'm'}{abs(int(offset)):03d}"
                inp, attempt_id = native_job(folder, point_id, "freq", pilot.input_text(selected["xyz_A"], "freq"),
                                             {"offset_deg": offset, "selection": selection}, ledger, contract)
                frequency = pilot.frequency_receipt(inp, np.asarray(selected["xyz_A"]), prepared["origin"], contract["policy"])
                witness = private_mpi_runtime_witness(inp.with_suffix(".out").read_text(), Path(contract["runtime"]["prefix"]))
                frequency["runtime_witness"] = witness
                seal_or_replay_analysis(ledger, attempt_id, frequency)
                selected["frequency"] = frequency
                if frequency["pilot_point_gate_passed"] is not True:
                    raise ScientificGateFailure("Selected wider-angle constrained spectrum/transfer gate unresolved")
                coupling = pilot.orthogonal_comparison([points[0], selected], contract["temperature_grid_K"])
                details["selected_orthogonal_comparison"] = coupling
                require_coupling(coupling, contract["policy"])
                selected["status"] = "derivative_verified_constrained"
            finish(folder, contract, ledger, points, **details)
        except ScientificGateFailure as error:
            finish(folder, contract, ledger, points, failure=str(error), **details)
        except BaseException as error:
            # Do not make a parser interruption irreversible or spend another job.
            try:
                counts = ledger.counts()
            except (OSError, ValueError):
                counts = {"status": "unavailable_retained_journal_requires_inspection"}
            checkpoint(folder, {"stage": "analysis_or_execution_attention", "status": "attention_required",
                                "failure_reason": f"{type(error).__name__}: {error}",
                                "whole_increment_attempts": counts, "rate_claim_allowed": False})
            print("ATTENTION REQUIRED; retained native evidence will not be redispatched: " + str(error), flush=True)
            raise
        freeze_rate_assessment(folder)
    return folder


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Completed three-point pilot, or this profile's resume directory")
    args = parser.parse_args(argv)
    run(args.source)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
