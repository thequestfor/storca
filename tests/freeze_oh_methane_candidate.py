"""Freeze this bounded OH/CH4 diagnostic, never launch chemistry or admit rates.

The small bundle preserves participating Python bytes independently of future
checkout edits. Reopening is byte-integrity verification, not re-execution or a
new scientific certificate. Native artifacts remain external, hash-bound files.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
KIND = "OH_CH4_frozen_conditional_rate_candidate"
FLAGS = ("rate_claim_allowed", "lifetime_claim_allowed", "network_reinjection_allowed",
         "production_defaults_changed", "calibrated_uncertainty_claim_allowed",
         "physical_rate_accuracy_claim_allowed", "angular_partition_admitted")


def read(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate JSON key: " + key)
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError("Nonfinite JSON")))


def artifact(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError("Expected source file")
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def reopen(binding):
    if (not isinstance(binding, dict) or not isinstance(binding.get("path"), str)
            or not Path(binding["path"]).is_absolute()
            or not isinstance(binding.get("sha256"), str)
            or re.fullmatch(r"[0-9a-f]{64}", binding["sha256"]) is None
            or artifact(binding["path"]) != {k: binding[k] for k in ("path", "sha256")}):
        raise ValueError("Frozen candidate source hash changed or malformed")
    return Path(binding["path"])


def bindings(value):
    """Only explicit records; never scan the runs tree."""
    if isinstance(value, dict):
        if {"path", "sha256"} <= value.keys() and Path(value["path"]).is_absolute():
            yield {k: value[k] for k in ("path", "sha256")}
        for item in value.values():
            yield from bindings(item)
    elif isinstance(value, list):
        for item in value:
            yield from bindings(item)


def journal_payload_bindings(result):
    """Expand only the closed result's explicitly bound attempt-event files."""
    for record in result["attempt_evidence"]:
        yield from bindings(read(reopen(record)))


def load(name):
    spec = importlib.util.spec_from_file_location("candidate_freeze_" + name,
                                                Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_campaign(folder):
    """Replay source/native admission without calling run, execute or dispatch."""
    runtime = load("resume_oh_methane_angular_runtime")
    contract, _ = runtime.candidate.reopen_campaign(folder)
    ledger = runtime.old.ledger_for(folder, contract)
    extension = runtime.validate_extension(folder, contract, ledger)
    if reopen(read(folder / "result-binding.json")) != folder / "result.json":
        raise ValueError("Closed result binding differs")
    result = read(folder / "result.json")
    if (result.get("contract") != artifact(folder / "contract.json")
            or result.get("whole_increment_attempts") != ledger.counts()
            or result.get("attempt_evidence") != ledger.evidence()
            or result.get("runtime_policy_extension") != extension):
        raise ValueError("Closed campaign/journal/runtime extension differs")
    runtime.old.reopen_bindings(result)
    consumer = load("replay_oh_methane_angular_rates")
    assessment = consumer.assess(folder)
    if not runtime.old.continuation.replay_equal(assessment, read(folder / "angular-rate-assessment.json")):
        raise ValueError("Frozen angular-rate refusal differs")
    return contract, result, assessment


def cost_records(records):
    """All identifiable angular jobs, including failed attempts; no inferred zero costs."""
    by_path = {row["path"]: row for row in records}
    rows = []
    for name in sorted(by_path):
        inp = Path(name)
        if inp.suffix != ".inp":
            continue
        angular = any(part.startswith(("oh-methane-angular-", "oh-methane-periodic-")) for part in inp.parts)
        diagnostic = any(part.startswith("oh-methane-gradient-mpi-diagnostic-") for part in inp.parts)
        if not (angular or diagnostic):
            continue
        text = reopen(by_path[name]).read_text()
        kind = "EnGrad" if diagnostic else ("NumFreq" if re.search(r"\bNumFreq\b", text, re.I) else "saddles")
        output, lease = inp.with_suffix(".out"), inp.with_suffix(".process.json")
        execution = inp.parent / ("case-result.json" if diagnostic else "execution.json")
        row = {"input": by_path[name], "kind": kind, "elapsed_seconds": None,
               "native_seconds": None, "ncores": None, "requested_core_seconds": None,
               "native_normal_marker_present": None, "native_normal_termination": None,
               "execution_failure": None}
        for path, key in ((execution, "execution"), (lease, "lease"), (output, "output")):
            if str(path) not in by_path:
                continue
            row[key] = by_path[str(path)]
            source = reopen(row[key])
            if key == "output":
                native = source.read_text(errors="replace")
                row["native_normal_marker_present"] = native.count("ORCA TERMINATED NORMALLY") == 1
                match = re.findall(r"TOTAL RUN TIME:\s*(\d+) days\s*(\d+) hours\s*(\d+) minutes\s*(\d+) seconds\s*(\d+) msec", native)
                if len(match) == 1:
                    d, h, m, s, ms = map(int, match[0])
                    row["native_seconds"] = d * 86400 + h * 3600 + m * 60 + s + ms / 1000
            else:
                try:
                    data = read(source)
                except (ValueError, UnicodeError) as error:
                    row[key + "_reader_warning"] = type(error).__name__
                    continue
                if key == "lease":
                    cores = data.get("ncores")
                    row["ncores"] = cores if type(cores) is int and cores > 0 else None
                    row["returncode"] = data.get("returncode")
                    row["lease_status"] = data.get("status")
                else:
                    elapsed = data.get("elapsed_seconds")
                    if type(elapsed) in (float, int) and math.isfinite(elapsed) and elapsed >= 0:
                        row["elapsed_seconds"] = elapsed
                    row["execution_failure"] = data.get("failure", data.get("failure_reason"))
        if row["elapsed_seconds"] is not None and row["ncores"] is not None:
            row["requested_core_seconds"] = row["elapsed_seconds"] * row["ncores"]
        if row["native_normal_marker_present"] is False:
            row["native_normal_termination"] = False
        elif row["native_normal_marker_present"] is True and type(row.get("returncode")) is int:
            row["native_normal_termination"] = (row["returncode"] == 0
                and row.get("lease_status") in (None, "terminated"))
        rows.append(row)
    return rows


def scientific_summary(result):
    points = result.get("points", [])
    energies = [point["energy_Eh"] for point in points if type(point.get("energy_Eh")) in (int, float)]
    return {"profile_status": result["status"], "failure_reason": result.get("failure_reason"),
            "points": [{"offset_deg": point["offset_deg"], "status": point["status"],
                        "energy_Eh": point.get("energy_Eh"), "frequency": point.get("frequency")}
                       for point in points],
            "sampled_electronic_spread_Eh": max(energies) - min(energies) if energies else None,
            "closure": result.get("closure"), "coverage": result.get("coverage"),
            "selected_orthogonal_comparison": result.get("selected_orthogonal_comparison"),
            "partition_assessment": result.get("partition_assessment"),
            "scope": "Discrete constrained candidates and selected tangent Hessians, not an ordinary TS inventory or partition."}


def freeze(campaignfolder, outputfolder):
    folder, destination = Path(campaignfolder).resolve(strict=True), Path(outputfolder).absolute()
    if destination.exists() or destination.is_symlink():
        raise ValueError("Candidate freeze requires a fresh nonexistent output folder")
    contract, result, assessment = validate_campaign(folder)
    recipe_path = reopen(contract["baseline_recipe"])
    recipe = read(recipe_path)
    state_path = reopen(contract["previous_conditional_branches"])
    state_result = read(state_path)
    state_contract = read(reopen(state_result["contract"]))
    state_native = read(reopen(state_result["native_replay"]))
    state_reader = load("replay_oh_methane_electronic_states")
    state_reader.validate_native_receipt(state_contract, state_native)
    _, closure, _ = state_reader.load_reporter(recipe_path.parent.parent)
    direct = [folder / name for name in ("contract.json", "result.json", "result-binding.json",
              "angular-rate-assessment.json", "runtime-policy-extension-v1.json",
              "runtime-policy-extension-binding.json")]
    records = (list(bindings([result, assessment, state_result]))
               + list(journal_payload_bindings(result)) + [artifact(path) for path in direct])
    for inventory in (contract["source_inventory"], contract["runtime"]["inventory"],
                      state_contract["retained_source_inventory"],
                      read(folder / "runtime-policy-extension-v1.json")["source_inventory"]):
        records.extend({k: row[k] for k in ("path", "sha256")} for row in inventory.values())
    records.extend(bindings([recipe, state_contract["historic_generator_snapshot"], closure]))
    records.append(artifact(__file__))
    # Historic generator bindings can name an intentionally changed live file;
    # its exact original bytes already have an admitted reconstruction snapshot.
    historical = state_contract["historic_generator_snapshot"]
    historic_original = historical["original_binding"]
    records = [row for row in records if row != historic_original]
    unique = {(row["path"], row["sha256"]): row for row in records}
    records = list(unique.values())
    for row in records:
        reopen(row)
    costs = cost_records(records)
    destination.mkdir(parents=True, exist_ok=False)
    snapshots = destination / "source-snapshots"
    snapshots.mkdir()
    source_map, external = [], []
    for index, row in enumerate(sorted(records, key=lambda row: (row["path"], row["sha256"]))):
        if Path(row["path"]).suffix != ".py":
            external.append(row)
            continue
        target = snapshots / (str(index) + "-" + Path(row["path"]).name)
        with target.open("xb") as stream:
            stream.write(reopen(row).read_bytes())
        source_map.append({"original": row, "snapshot": artifact(target)})
    source_map.append({"original": historic_original,
                       "snapshot": next(pair["snapshot"] for pair in source_map
                                        if pair["original"] == historical["snapshot"])})
    branch_specs = state_contract["branches"]
    if {row["TS_optical_isomers"] for row in branch_specs} != {1, 2}:
        raise ValueError("Both unresolved counting alternatives must remain visible")
    document = {"schema_version": 1, "kind": KIND, "status": "diagnostic_not_admitted",
                **{flag: False for flag in FLAGS}, "recipe_candidate_has_unresolved_assumptions": True,
                "new_native_jobs": 0, "candidate_revision": "oh_ch4_bound_conditional_recipe_v1",
                "campaign": artifact(folder / "result.json"), "recipe": artifact(recipe_path),
                "electronic_state_replay": artifact(state_path), "source_bindings": external,
                "source_snapshots": source_map, "scientific_summary": scientific_summary(result),
                "conditional_branches": branch_specs, "selected_branch": None,
                "method_assumptions": {key: recipe["recipe"][key] for key in (
                    "electronic_method", "vibrational_method", "frequency_scale_factor", "environment",
                    "surface_model", "electronic_state_assumption", "channel_counting_assumption",
                    "soft_mode_assumption", "tunneling", "temperatures_K")},
                "whole_increment_attempts": result["whole_increment_attempts"],
                "caps": contract["caps"], "separate_EnGrad_diagnostic_attempts": result["prior_EnGrad_execution_diagnostic_attempts"],
                "cost_records": costs,
                "identified_attempt_counts": {kind: sum(row["kind"] == kind for row in costs)
                                              for kind in ("saddles", "NumFreq", "EnGrad")},
                "cost_scope": "Bound angular increment only. All charged attempts remain; missing timings are null. Requested core-seconds are not measured CPU time. Earlier geometry/composite development is excluded.",
                "known_elapsed_seconds_subtotal": sum(row["elapsed_seconds"] for row in costs if row["elapsed_seconds"] is not None),
                "elapsed_coverage_jobs": sum(row["elapsed_seconds"] is not None for row in costs),
                "identified_jobs": len(costs), "limitations": [
                    "No reviewed global angular Hamiltonian/measure, continuous transverse-bath proof or reconciled basin/symmetry counting.",
                    "A constrained tangent saddle, including +170 degrees, is not an ordinary first-order transition state.",
                    "FrozenCore canonical CCSD(T), finite QZ basis and single-reference DFT geometries/frequencies are not ground truth.",
                    "Isolated ideal-gas elementary bimolecular k(T), not a pressure/solvent model, network completeness proof or lifetime."]
                    + [item["reason"] for item in assessment["missing_evidence"]]
                    + state_result["limitations"],
                "reopen_scope": "Byte integrity of native sources and copied implementation versions; exact historical semantic execution is separate."}
    for name, value in (("candidate.json", document),):
        with (destination / name).open("x") as stream:
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    with (destination / "candidate-binding.json").open("x") as stream:
        json.dump(artifact(destination / "candidate.json"), stream, indent=2, sort_keys=True)
        stream.write("\n")
    return verify(destination)


def verify(folder_or_candidate_json):
    source = Path(folder_or_candidate_json).resolve(strict=True)
    folder = source if source.is_dir() else source.parent
    if reopen(read(folder / "candidate-binding.json")) != folder / "candidate.json":
        raise ValueError("Candidate seal differs")
    document = read(folder / "candidate.json")
    if (document.get("schema_version") != 1 or document.get("kind") != KIND
            or document.get("status") != "diagnostic_not_admitted"
            or any(document.get(flag) is not False for flag in FLAGS)
            or document.get("recipe_candidate_has_unresolved_assumptions") is not True
            or type(document.get("new_native_jobs")) is not int or document["new_native_jobs"] != 0
            or document.get("selected_branch") is not None
            or {row["TS_optical_isomers"] for row in document.get("conditional_branches", [])} != {1, 2}):
        raise ValueError("Frozen conditional claim boundary differs")
    for row in document["source_bindings"]:
        reopen(row)
    for pair in document["source_snapshots"]:
        if pair["original"]["sha256"] != pair["snapshot"]["sha256"]:
            raise ValueError("Implementation snapshot differs from original revision")
        path = reopen(pair["snapshot"])
        if folder / "source-snapshots" not in path.parents:
            raise ValueError("Source snapshot is outside this bundle")
    return document


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "reopen"))
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", nargs="?", type=Path)
    args = parser.parse_args()
    if args.action == "freeze" and args.destination is None:
        parser.error("freeze requires a fresh destination")
    result = freeze(args.source, args.destination) if args.action == "freeze" else verify(args.source)
    print(json.dumps({"status": result["status"], "whole_increment_attempts": result["whole_increment_attempts"],
                      "new_native_jobs": result["new_native_jobs"]}, indent=2))
