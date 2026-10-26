"""Authorized energy-only PMIx hash continuation; original native ledger retained."""
import argparse
from contextlib import nullcontext
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import resume_h_ch4_remaining_work as previous
from storca import rate_transfer_preparation as prep
from storca.rate_transfer_control import NativeExecutionControl
from storca.scientific_panel import canonical_hash

KIND = "h_ch4_pmix_hash_energy_continuation"
CHANGED = {"storca/rate_transfer_control.py", "storca/rate_transfer_execution.py"}


def _sealed(binding, kind):
    value = prep._json(prep._reopen(binding))
    if (set(value) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or value["schema_version"] != 1 or value["kind"] != kind + "_sealed"
            or value["plan"]["authorized"] is not True
            or canonical_hash(prep._json(prep._reopen(value["plan_source"]))) != canonical_hash(value["plan"])):
        raise ValueError("Exact authorized historical contract required")
    return value


def _files(inventory, rows, folder):
    required = {r["path"]: r["sha256"] for r in inventory}
    if len(required) != len(inventory):
        raise ValueError("Duplicate historical implementation")
    found, paths = {}, set()
    for row in rows:
        relative = Path(row["relative_path"])
        if (relative.is_absolute() or ".." in relative.parts or str(relative) in paths
                or row["original_path"] in found or required.get(row["original_path"]) != row["sha256"]):
            raise ValueError("Historical source inventory differs")
        path = folder / relative
        if path.is_symlink() or not path.resolve().is_relative_to(folder.resolve()):
            raise ValueError("Historical source escaped archive")
        prep._reopen({"path": str(path), "sha256": row["sha256"]})
        found[row["original_path"]] = row["sha256"]
        paths.add(str(relative))
    if found != required:
        raise ValueError("Incomplete historical source archive")


def archive(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    sealed = previous.reopen(source)
    if destination.exists():
        raise ValueError("Archive into a fresh directory")
    destination.mkdir(parents=True)
    rows = []
    for item in sealed["implementation"]:
        path = prep._reopen(item)
        relative = path.relative_to(ROOT)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        rows.append({"original_path": str(path), "relative_path": str(relative), "sha256": item["sha256"]})
    shutil.copy2(source / "resume-contract.json", destination / "resume-contract.json")
    capsule = {"schema_version": 1, "kind": KIND + "_source_capsule",
               "previous_contract": prep.binding(source / "resume-contract.json"), "files": rows,
               "native_jobs_launched": 0}
    prep._publish(destination / "source-capsule.json", capsule)
    return {"status": "sources_archived", "files": len(rows), "native_jobs_launched": 0}


def _history(binding):
    capsule_path = prep._reopen(binding)
    capsule = prep._json(capsule_path)
    if (set(capsule) != {"schema_version", "kind", "previous_contract", "files", "native_jobs_launched"}
            or capsule["schema_version"] != 1 or capsule["kind"] != KIND + "_source_capsule"
            or type(capsule["native_jobs_launched"]) is not int or capsule["native_jobs_launched"] != 0):
        raise ValueError("Exact previous energy-continuation source archive required")
    old = _sealed(capsule["previous_contract"], previous.KIND)
    if prep.binding(capsule_path.parent / "resume-contract.json")["sha256"] != capsule["previous_contract"]["sha256"]:
        raise ValueError("Archived previous execution contract changed")
    _files(old["implementation"], capsule["files"], capsule_path.parent)
    for row in old["implementation"]:
        if Path(row["path"]).relative_to(ROOT).as_posix() not in CHANGED:
            prep._reopen(row)
    recovery = previous.recovery
    prior = _sealed(old["plan"]["previous_contract"], recovery.KIND)
    archive_path = prep._reopen(old["plan"]["source_capsule"])
    archived_recovery = prep._json(archive_path)
    if (archived_recovery.get("kind") != previous.KIND + "_source_capsule"
            or archived_recovery.get("previous_contract") != old["plan"]["previous_contract"]
            or archived_recovery.get("native_jobs_launched") != 0):
        raise ValueError("Prior recovery source lineage changed")
    _files(prior["implementation"], archived_recovery["files"], archive_path.parent)
    _, original, settings, control, _ = recovery._validate_plan(
        prior["plan"], Path(old["plan"]["previous_contract"]["path"]).parent,
        diagnostic_implementation=prior["implementation"])
    recovery._prefix(old["plan"], control)
    old_workdir = Path(old["plan"]["previous_contract"]["path"]).parent
    return old_workdir, original, settings, control


def _inventory():
    return sorted([*previous._inventory(), prep.binding(Path(__file__))], key=lambda r: r["path"])


def prepare(capsule, destination, *, authorized=False):
    capsule, destination = Path(capsule).resolve(), Path(destination).resolve()
    if capsule.is_dir():
        capsule /= "source-capsule.json"
    old_workdir, original, settings, control = _history(prep.binding(capsule))
    snapshot = control.snapshot()
    if snapshot["pending"] or snapshot["counts"] != {"optimization_path": 9, "frequency": 7, "irc": 2, "single_point": 6}:
        raise ValueError("Expected closed 24-attempt checkpoint required")
    if destination.exists():
        raise ValueError("Prepare a fresh runtime continuation")
    destination.mkdir(parents=True)
    execution = previous.recovery.execution
    route = execution.prediction_to_route(execution._prediction(prep._json(prep._reopen(original["plan"]["prediction"]))))
    source_path = old_workdir / "path/route-verification.json"
    path_copy = destination / "verified-path.json"
    prep._publish(path_copy, prep._json(source_path))
    points, _ = execution._verified_points({"path": path_copy}, route)
    retry = old_workdir / "energies/r2/pmix-hash-retry-001"
    folders = {label: str(retry if label == "r2" else old_workdir / "energies" / label)
               for label in execution.LABELS}
    inputs = [prep.binding(execution._single_point_input(points[label], Path(folders[label]), original["plan"])[0])
              for label in execution.LABELS]
    failed = old_workdir / "energies/r2/energy.inp"
    if Path(inputs[1]["path"]).read_bytes() != failed.read_bytes():
        raise ValueError("Retry changed the frozen molecular input")
    owned = [r for r in snapshot["reservations"].values() if r["input"] == str(failed)]
    if len(owned) != 1 or snapshot["terminals"][owned[0]["attempt_id"]]["outcome"] != "failed":
        raise ValueError("Exact owned failed methane attempt required")
    plugin = Path(control.contract["extra_environment"]["DYLD_LIBRARY_PATH"]) / "pmix/mca_gds_hash.so"
    amendment = {"schema_version": 1, "kind": "storca_pmix_hash_single_point_amendment",
                 "authorized": authorized, "original_control_contract": prep.binding(control.root / "contract.json"),
                 "original_single_point_cap": 9, "revised_single_point_cap": 10,
                 "journal_prefix": [prep.binding(p) for p in sorted((control.root / "events").iterdir())],
                 "prefix_spending": {k: snapshot[k] for k in ("counts", "wall_seconds", "core_seconds")},
                 "extra_environment": {"PMIX_MCA_gds": "hash"}, "runtime_source": prep.binding(plugin),
                 "allowed_inputs": inputs, "failed_input": prep.binding(failed),
                 "failed_output": prep.binding(failed.with_suffix(".out"))}
    prep._publish(destination / "execution-amendment.json", amendment)
    plan = {"schema_version": 1, "kind": KIND, "authorized": authorized,
            "source_capsule": prep.binding(capsule), "execution_amendment": prep.binding(destination / "execution-amendment.json"),
            "retained_verified_path": {"path_summary": prep.binding(path_copy),
                                       "handoff": prep.binding(old_workdir / "rate-source-contract.json")},
            "energy_directories": folders, "inputs": inputs}
    prep._publish(destination / "runtime-plan.json", plan)
    prep._publish(destination / "runtime-contract.json", {"schema_version": 1, "kind": KIND + "_sealed",
                  "plan": plan, "plan_source": prep.binding(destination / "runtime-plan.json"),
                  "implementation": _inventory()})
    return inspect(destination)


def reopen(destination):
    destination = Path(destination).resolve()
    sealed = _sealed(prep.binding(destination / "runtime-contract.json"), KIND)
    for row in sealed["implementation"]:
        prep._reopen(row)
    if sealed["implementation"] != _inventory():
        raise ValueError("Runtime-continuation implementation changed")
    plan = sealed["plan"]
    if (set(plan) != {"schema_version", "kind", "authorized", "source_capsule", "execution_amendment",
                      "retained_verified_path", "energy_directories", "inputs"}
            or plan["schema_version"] != 1 or plan["kind"] != KIND or plan["authorized"] is not True):
        raise ValueError("Invalid runtime continuation plan")
    old_workdir, original, settings, control = _history(plan["source_capsule"])
    amended = NativeExecutionControl(control.root, control.contract, execution_amendment=plan["execution_amendment"])
    amended.validate_execution_amendment()
    execution = previous.recovery.execution
    route = execution.prediction_to_route(execution._prediction(prep._json(prep._reopen(original["plan"]["prediction"]))))
    path_result = execution._retained_verified_path(plan["retained_verified_path"], route, amended)
    points, _ = execution._verified_points(path_result, route)
    expected = {label: str(old_workdir / "energies" / label / "pmix-hash-retry-001" if label == "r2"
                           else old_workdir / "energies" / label) for label in execution.LABELS}
    if plan["energy_directories"] != expected or plan["inputs"] != amended.validate_execution_amendment()["allowed_inputs"]:
        raise ValueError("Matched-energy destinations changed")
    for label, source in zip(execution.LABELS, plan["inputs"]):
        prep._reopen(source)
        inp, _ = execution._single_point_input(points[label], Path(expected[label]), original["plan"])
        if prep.binding(inp) != source:
            raise ValueError("Frozen matched-energy input changed")
    return sealed, original, settings, amended


def inspect(destination):
    sealed, _, _, control = reopen(destination)
    s = control.snapshot()
    return {"status": "pmix_hash_energy_continuation_ready", "authorized": sealed["plan"]["authorized"],
            "counts": s["counts"], "pending": len(s["pending"]), "native_jobs_launched": 0,
            "contract": prep.binding(Path(destination).resolve() / "runtime-contract.json"),
            **previous.recovery.execution.FLAGS}


def run(destination, *, rmg_env):
    destination = Path(destination).resolve()
    sealed, original, settings, control = reopen(destination)
    plan = sealed["plan"]
    from storca.rate_transfer_runtime import runtime_environment
    runtime = runtime_environment(original["plan"]["runtime_inventory"]) if original["plan"]["schema_version"] == 2 else nullcontext()
    with runtime:
        (destination / "campaign").mkdir(exist_ok=True)
        return previous.recovery.execution._execute_validated_context(
            destination / "campaign", plan=original["plan"], settings=settings, rmg_env=rmg_env,
            control_root=control.root, context_contract_path=destination / "runtime-contract.json",
            validate_context=lambda: reopen(destination), execution_amendment=plan["execution_amendment"],
            retained_verified_path=plan["retained_verified_path"], energy_directories=plan["energy_directories"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    a = sub.add_parser("archive")
    a.add_argument("source", type=Path)
    a.add_argument("destination", type=Path)
    p = sub.add_parser("prepare")
    p.add_argument("capsule", type=Path)
    p.add_argument("destination", type=Path)
    p.add_argument("--authorized", action="store_true")
    sub.add_parser("inspect").add_argument("destination", type=Path)
    r = sub.add_parser("run")
    r.add_argument("destination", type=Path)
    r.add_argument("--rmg-env", required=True)
    args = parser.parse_args()
    if args.command == "archive":
        result = archive(args.source, args.destination)
    elif args.command == "prepare":
        result = prepare(args.capsule, args.destination, authorized=args.authorized)
    else:
        result = inspect(args.destination) if args.command == "inspect" else run(args.destination, rmg_env=args.rmg_env)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
