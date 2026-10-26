"""Correct candidate acquisition and continue the retained periodic profile.

Never rewrite the closed 2026-10-15 run or repeat its +60-degree optimization.
Select the actual indexed native gradient geometry, not a relabelled EnGrad
header. Native convergence admits a candidate; unchanged physical checks and
selected full Hessians remain separate evidence. No corrected rate is enabled.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import math
from pathlib import Path
import re
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np

from storca.constrained_candidate_eligibility import orca_constrained_candidate
from storca.orca_evaluated_gradient import evaluated_gradient_receipt

spec = importlib.util.spec_from_file_location(
    "retained_periodic_profile", Path(__file__).with_name("run_oh_methane_angular_profile.py"))
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
pilot, state, profile = old.pilot, old.state, old.profile
KIND = old.KIND
REVISION = "native_convergence_evaluated_gradient_v1"
EXPECTED_FAILURE = "Native optimization candidate eligibility unresolved"
PRIOR = {"saddles": 7, "NumFreq": 3}
SCHEDULE = tuple(profile.ACQUISITION_SCHEDULE[1:])
INHERITED_POLICY_KEYS = (
    "schema_version", "kind", "caps", "policy", "mapped_labels", "dihedral_atom_indices", "method",
    "cores", "maxcore_MB_per_rank", "per_job_seconds", "minimum_disk_reserve_bytes", "automatic_retries",
    "temperature_grid_K", "masses_amu", "center_angle_rad", "grid_offsets_deg", "runtime",
    "selected_environment", "initial_hessian", "closure_policy", "counting_policy", "partition_policy",
    "baseline_recipe", "previous_conditional_branches", "physical_policy_unchanged",
    "prior_EnGrad_execution_diagnostic_attempts",
)
CANDIDATE_POLICY = {
    "native_convergence_required": True, "extra_displacement_cutoff_applied": False,
    "independent_evaluated_frame_gradient_required": True,
    "gradient_geometry_policy": "same native property GeometryIndex for E/g/geometry",
    "independent_hessian_required_for_saddle_claim": True,
}


def validate_inherited_policy(contract, parent_contract):
    if any(not old.continuation.replay_equal(contract[key], parent_contract[key]) for key in INHERITED_POLICY_KEYS):
        raise ValueError("Inherited scientific/execution policy differs from parent")


def dispatch(folder, point_id, stem, text, metadata, ledger, contract, parent_ledger):
    """Refuse changed evidence before reserving; old executor checks again after."""
    state.reopen_inventory(contract["source_inventory"])
    state.reopen_inventory(contract["runtime"]["inventory"])
    if parent_ledger.evidence() != contract["parent_attempt_evidence"]:
        raise ValueError("Parent ledger changed before dispatch")
    return old.native_job(folder, point_id, stem, text, metadata, ledger, contract)


def analyze_candidate(inp, seed, offset, origin, contract):
    """Admit the same-index evaluated E/g/geometry, never a mixed native file."""
    old.recover_execution(inp, pilot.input_text(seed, "opt", pilot.reopen(contract["initial_hessian"])))
    text = inp.with_suffix(".out").read_text()
    candidate = orca_constrained_candidate(text)
    primitive = r"D\((?:H\s+1,O\s+0,C\s+2,H\s+4|H\s+4,C\s+2,O\s+0,H\s+1)\)[^\n]*\bC\s*$"
    if not re.search(primitive, text, re.M):
        raise ValueError("Native mapped constrained primitive missing")
    labels, final_xyz, _ = pilot.read_xyz(inp.with_suffix(".xyz"))
    if tuple(labels) != pilot.LABELS or final_xyz.shape != (7, 3) or not np.isfinite(final_xyz).all():
        raise ValueError("Native final mapped geometry differs")
    property_path = inp.with_suffix(".property.txt")
    evaluation = evaluated_gradient_receipt(property_path.read_text(), text, final_xyz, pilot.BOHR_A)
    xyz = np.asarray(evaluation["xyz_A"])
    xyz_bohr = np.asarray(evaluation["xyz_bohr"])
    gradient = np.asarray(evaluation["gradient_Eh_bohr"])
    diagnostics = profile.candidate_gradient_diagnostics(
        xyz_bohr, origin["masses_amu"], gradient,
        target_angle_rad=contract["center_angle_rad"] + math.radians(offset))
    native = pilot.native_state(inp)
    if max(abs(value - .758947) for value in native["S_squared"]) > contract["policy"]["maximum_S_squared_change"]:
        raise old.ScientificGateFailure("Native spin continuity gate unresolved")
    if (diagnostics["candidate_gate_passed"] is not True
            or diagnostics["tangent_gradient_RMS_Eh_bohr"] > contract["policy"]["tangent_gradient_RMS_Eh_bohr"]
            or diagnostics["tangent_gradient_max_Eh_bohr"] > contract["policy"]["tangent_gradient_max_Eh_bohr"]):
        raise old.ScientificGateFailure("Evaluated-frame constrained tangent-gradient gate unresolved")
    topology = old.topology_guard(xyz)
    witness = old.private_mpi_runtime_witness(text, Path(contract["runtime"]["prefix"]))
    candidate.update({"input": pilot.artifact(inp), "output": pilot.artifact(inp.with_suffix(".out")),
                      "native_final_xyz": pilot.artifact(inp.with_suffix(".xyz")), "native_state": native,
                      "evaluated_geometry_index": evaluation["selected_geometry_index"],
                      "constraint_gap_rad": diagnostics["constraint_gap_rad"],
                      "force_geometry_source": "same-index native property geometry/gradient/dispersion-inclusive energy"})
    return {"offset_deg": float(offset), "xyz_A": xyz.tolist(), "xyz_bohr": xyz_bohr.tolist(),
            "energy_Eh": evaluation["energy_Eh"], "gradient_Eh_bohr": gradient.tolist(),
            "gradient_diagnostics": diagnostics, "candidate_eligibility": candidate,
            "evaluated_gradient": {"property_file": pilot.artifact(property_path), **evaluation},
            "original_engrad": pilot.artifact(inp.with_suffix(".engrad")),
            "original_engrad_header_used_for_gradient_geometry": False,
            "spectator_topology": topology, "runtime_witness": witness,
            "status": "constrained_optimizer_candidate", "saddle_character_verified_by_optimization_alone": False}


def parent_claim(parent, contract):
    source = pilot.reopen(contract["source"]).parent
    path = old.claim_path(source)
    claim = pilot.read_json(path)
    if claim.get("source") != contract["source"] or Path(claim.get("continuation_directory", "")).resolve() != parent:
        raise ValueError("Original periodic source claim differs")
    return path


def preflight_parent(parent):
    parent = old.owned_run(parent)
    seal = pilot.read_json(parent / "result-binding.json")
    if pilot.reopen(seal) != parent / "result.json":
        raise ValueError("Sealed parent result location differs")
    result = pilot.read_json(parent / "result.json")
    if (result.get("kind") != KIND or result.get("status") != "profile_unresolved"
            or result.get("failure_reason") != EXPECTED_FAILURE
            or result.get("whole_increment_attempts") != PRIOR
            or result.get("old_sources_unchanged") is not True
            or [row["offset_deg"] for row in result.get("points", [])] != [0.0, -15.0, 15.0]):
        raise ValueError("Exact closed displacement-screen predecessor required")
    old.assert_unadmitted(result)
    contract, prepared = old.reopen_campaign(parent)
    if pilot.reopen(result["contract"]) != parent / "contract.json":
        raise ValueError("Parent result/contract binding differs")
    ledger = old.ledger_for(parent, contract)
    if ledger.counts() != PRIOR or result["attempt_evidence"] != ledger.evidence():
        raise ValueError("Parent charged attempts/history differ")
    old.reopen_bindings(result)
    reservations = ledger.reservations()
    if (len(reservations) != 1 or reservations[0]["point_id"] != "point-p060"
            or reservations[0]["kind"] != "saddles" or ledger.analysis(reservations[0]) is not None):
        raise ValueError("One paid, unanalyzed +60 optimization required")
    attempt = reservations[0]
    execution = ledger.execution(attempt)
    if execution is None or execution.get("status") != "normal":
        raise ValueError("Normally completed retained optimization required")
    old.reopen_bindings(execution)
    metadata = attempt["metadata"]
    inp = Path(metadata["input_path"])
    if inp != parent / "point-p060/optimization/opt.inp" or metadata.get("offset_deg") != 60 or metadata.get("from_offset_deg") != 15:
        raise ValueError("Retained +60 input/seed identity differs")
    if pilot.artifact(inp)["sha256"] != metadata["input_sha256"]:
        raise ValueError("Retained +60 input hash differs")
    seed = np.asarray(metadata["seed_A"])
    old.recover_execution(inp, pilot.input_text(seed, "opt", pilot.reopen(contract["initial_hessian"])))
    adopted = analyze_candidate(inp, seed, 60, prepared["origin"], contract)
    adopted.update({"optimization_reused": True, "new_optimization_attempts": 0,
                    "source_attempt": attempt, "source_execution": execution})
    paths = {Path(path) for path in contract["source_inventory"]}
    paths.update(state.bound_paths(result))
    paths.update(state.bound_paths(adopted))
    paths.update(Path(record["path"]) for record in ledger.evidence())
    paths.update({parent / "result.json", parent / "result-binding.json", parent / "contract.json",
                  parent / "angular-rate-assessment.json", parent_claim(parent, contract)})
    return {"parent": parent, "parent_contract": contract, "parent_result": result,
            "origin": prepared["origin"], "points": prepared["points"], "adopted": adopted,
            "parent_ledger": ledger, "paths": paths}


def correction_claim(parent):
    return ROOT / "runs" / ("oh-methane-periodic-candidate-claim-" + pilot.artifact(parent / "result.json")["sha256"] + ".json")


def make_campaign(prepared):
    parent = prepared["parent"]
    claim = correction_claim(parent)
    if claim.exists():
        saved = pilot.read_json(claim)
        if saved.get("source") != pilot.artifact(parent / "result.json"):
            raise ValueError("Candidate correction source claim differs")
        return old.owned_run(saved["continuation_directory"])
    folder = Path(tempfile.mkdtemp(prefix="oh-methane-periodic-candidate-20261016-", dir=ROOT / "runs"))
    old.atomic_save(claim, {"source": pilot.artifact(parent / "result.json"),
                            "continuation_directory": str(folder), "automatic_retries": 0})
    helpers = {Path(__file__), ROOT / "storca/constrained_candidate_eligibility.py",
               ROOT / "storca/orca_evaluated_gradient.py"}
    inventory = state.inventory(prepared["paths"] | helpers)
    snapshots = list(prepared["parent_contract"]["source_snapshots"])
    directory = folder / "source-snapshots"
    directory.mkdir()
    for index, source in enumerate(sorted(helpers)):
        snapshot = directory / f"{index}-{source.name}"
        with snapshot.open("xb") as stream:
            stream.write(source.read_bytes())
        snapshots.append({"source": pilot.artifact(source), "snapshot": pilot.artifact(snapshot)})
    old.atomic_save(folder / "reused-plus60.json", prepared["adopted"])
    contract = copy.deepcopy(prepared["parent_contract"])
    contract.update({"candidate_revision": REVISION, "source": pilot.artifact(parent / "result.json"),
                     "predecessor_contract": pilot.artifact(parent / "contract.json"),
                     "source_inventory": inventory, "source_snapshots": snapshots,
                     "prior_attempts": dict(prepared["parent_ledger"].counts()),
                     "parent_attempt_evidence": prepared["parent_ledger"].evidence(),
                     "acquisition_schedule": [list(pair) for pair in SCHEDULE],
                     "reused_candidate": pilot.artifact(folder / "reused-plus60.json"),
                     "candidate_policy": dict(CANDIDATE_POLICY),
                     "reader_revision": {"new_candidate_reader": pilot.artifact(ROOT / "storca/constrained_candidate_eligibility.py"),
                                         "new_gradient_frame_reader": pilot.artifact(ROOT / "storca/orca_evaluated_gradient.py"),
                                         "reason": "Native step size is diagnostic; require a genuinely evaluated-frame gradient instead",
                                         "physical_stationarity_thresholds_changed": False,
                                         "old_files_rewritten": False}})
    old.atomic_save(folder / "contract.json", contract)
    return folder


def reopen_campaign(folder):
    folder = old.owned_run(folder)
    contract = pilot.read_json(folder / "contract.json")
    old.assert_unadmitted(contract)
    if (contract.get("kind") != KIND or contract.get("candidate_revision") != REVISION
            or contract["prior_attempts"] != PRIOR or contract["caps"] != old.CAPS
            or contract["acquisition_schedule"] != [list(pair) for pair in SCHEDULE]
            or not old.continuation.replay_equal(contract["candidate_policy"], CANDIDATE_POLICY)
            or contract["partition_policy"]["mode"] != "diagnostic_only"
            or contract["partition_policy"]["supported_quantitative_hamiltonian"] is not False):
        raise ValueError("Bounded evaluated-gradient continuation policy differs")
    state.reopen_inventory(contract["source_inventory"])
    state.reopen_inventory(contract["runtime"]["inventory"])
    for pair in contract["source_snapshots"]:
        if pilot.reopen(pair["source"]).read_bytes() != pilot.reopen(pair["snapshot"]).read_bytes():
            raise ValueError("Retained implementation snapshot differs")
    prepared = preflight_parent(pilot.reopen(contract["source"]).parent)
    validate_inherited_policy(contract, prepared["parent_contract"])
    if prepared["parent_ledger"].evidence() != contract["parent_attempt_evidence"]:
        raise ValueError("Parent history/physical/runtime policy differs")
    if not old.continuation.replay_equal(pilot.read_json(pilot.reopen(contract["reused_candidate"])), prepared["adopted"]):
        raise ValueError("Reused +60 evaluated candidate differs")
    return contract, prepared


def run(source):
    source = old.owned_run(source)
    maybe = pilot.read_json(source / "contract.json")
    folder = source if maybe.get("candidate_revision") == REVISION else make_campaign(preflight_parent(source))
    contract, prepared = reopen_campaign(folder)
    ledger = old.ledger_for(folder, contract)
    ledger.initialize()
    print("WORKDIR " + str(folder), flush=True)
    with prepared["parent_ledger"].controller(), ledger.controller():
        if (folder / "result.json").exists():
            pilot.reopen(pilot.read_json(folder / "result-binding.json"))
            result = pilot.read_json(folder / "result.json")
            if (result["contract"] != pilot.artifact(folder / "contract.json")
                    or result["whole_increment_attempts"] != ledger.counts()
                    or result["attempt_evidence"] != ledger.evidence()):
                raise ValueError("Closed correction result/journal differs")
            old.reopen_bindings(result)
            old.freeze_rate_assessment(folder)
            print("REUSED CLOSED RESULT " + result["status"], flush=True)
            return folder
        points = [*copy.deepcopy(prepared["points"]), copy.deepcopy(prepared["adopted"])]
        by_offset = {point["offset_deg"]: point for point in points}
        details = {"candidate_revision": REVISION, "reused_plus60": contract["reused_candidate"],
                   "parent_attempt_evidence": contract["parent_attempt_evidence"],
                   "old_optimization_rerun": False}
        try:
            with old.scientific_boundary(), old.frozen.scoped_environment(contract["selected_environment"]):
                for previous_offset, offset in SCHEDULE:
                    xyz = np.asarray(by_offset[previous_offset]["xyz_A"])
                    desired = contract["center_angle_rad"] + math.radians(offset)
                    delta = old.wrap_angle(desired - old.signed_dihedral(xyz, pilot.DIHEDRAL))
                    seed = old.rotate_dihedral_seed(xyz / pilot.BOHR_A, pilot.DIHEDRAL, delta) * pilot.BOHR_A
                    point_id = f"point-{'p' if offset >= 0 else 'm'}{abs(offset):03d}"
                    # Parent controller remains locked; inventory and exact
                    # ledger head are also rechecked before any native dispatch.
                    inp, attempt = dispatch(folder, point_id, "opt",
                        pilot.input_text(seed, "opt", pilot.reopen(contract["initial_hessian"])),
                        {"offset_deg": offset, "from_offset_deg": previous_offset, "seed_A": seed.tolist()}, ledger, contract,
                        prepared["parent_ledger"])
                    node = analyze_candidate(inp, seed, offset, prepared["origin"], contract)
                    old.seal_or_replay_analysis(ledger, attempt, node)
                    if offset == 190:
                        details["closure_control"] = node
                        reference = by_offset[-170]
                        details["closure"] = profile.branch_closure(reference["xyz_A"], node["xyz_A"], contract["masses_amu"],
                            reference["energy_Eh"], node["energy_Eh"], True, True)
                        details["spectator_permutation_diagnostics"] = profile.spectator_permutation_diagnostics(
                            reference["xyz_A"], node["xyz_A"], contract["masses_amu"])
                        if not details["closure"]["strict_labelled_closure_passed"]:
                            raise old.ScientificGateFailure("Independent periodic branch closure unresolved")
                    else:
                        points.append(node)
                        by_offset[offset] = node
                    old.checkpoint(folder, {"stage": "evaluated_frame_acquisition", "completed_offset_deg": offset,
                                            "whole_increment_attempts": ledger.counts(), "rate_claim_allowed": False})
                details["relaxed_path_metric_diagnostics"] = old.path_metric_diagnostics(points, contract["masses_amu"])
                selection = profile.select_final_hessian(points, contract["temperature_grid_K"], contract["masses_amu"])
                details["hessian_selection"] = selection
                if selection["selected_offset_deg"] is None:
                    raise old.ScientificGateFailure("Final-Hessian diagnostic selection unresolved")
                selected = by_offset[selection["selected_offset_deg"]]
                offset = selected["offset_deg"]
                point_id = f"point-{'p' if offset >= 0 else 'm'}{abs(int(offset)):03d}"
                inp, attempt = dispatch(folder, point_id, "freq", pilot.input_text(selected["xyz_A"], "freq"),
                                       {"offset_deg": offset, "selection": selection}, ledger, contract,
                                       prepared["parent_ledger"])
                receipt = pilot.frequency_receipt(inp, np.asarray(selected["xyz_A"]), prepared["origin"], contract["policy"])
                receipt["runtime_witness"] = old.private_mpi_runtime_witness(inp.with_suffix(".out").read_text(), Path(contract["runtime"]["prefix"]))
                old.seal_or_replay_analysis(ledger, attempt, receipt)
                selected["frequency"] = receipt
                if receipt["pilot_point_gate_passed"] is not True:
                    raise old.ScientificGateFailure("Selected wider-angle constrained spectrum/transfer gate unresolved")
                comparison = pilot.orthogonal_comparison([points[0], selected], contract["temperature_grid_K"])
                details["selected_orthogonal_comparison"] = comparison
                old.require_coupling(comparison, contract["policy"])
                selected["status"] = "derivative_verified_constrained"
            old.finish(folder, contract, ledger, points, **details)
        except old.ScientificGateFailure as error:
            old.finish(folder, contract, ledger, points, failure=str(error), **details)
        except BaseException as error:
            try:
                counts = ledger.counts()
            except (OSError, ValueError):
                counts = {"status": "unavailable_retained_journal_requires_inspection"}
            old.checkpoint(folder, {"stage": "analysis_or_execution_attention", "status": "attention_required",
                                    "failure_reason": f"{type(error).__name__}: {error}",
                                    "whole_increment_attempts": counts, "rate_claim_allowed": False})
            print("ATTENTION REQUIRED; native results will not be redispatched: " + str(error), flush=True)
            raise
        old.freeze_rate_assessment(folder)
    return folder


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Closed +60 displacement-screen stop, or this correction's resume directory")
    run(parser.parse_args().source)
