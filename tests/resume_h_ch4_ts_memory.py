"""One TS-only memory retry, with all four completed energies reused."""
import argparse
from contextlib import nullcontext
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import resume_h_ch4_h2_serial as previous
from storca import rate_transfer_preparation as prep
from storca.rate_transfer_control import NativeExecutionControl

KIND = "h_ch4_ts_memory_energy_continuation"
EXECUTION = previous.EXECUTION
CHANGED = {"storca/rate_transfer_control.py", "storca/rate_transfer_execution.py"}


def archive(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    sealed, _, _, _ = previous.reopen(source)
    if destination.exists():
        raise ValueError("Archive into a fresh directory")
    destination.mkdir(parents=True)
    rows = []
    for row in sealed["implementation"]:
        path = prep._reopen(row)
        relative = path.relative_to(ROOT)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        rows.append({"original_path": str(path), "relative_path": str(relative), "sha256": row["sha256"]})
    shutil.copy2(source / "runtime-contract.json", destination / "runtime-contract.json")
    prep._publish(destination / "source-capsule.json", {
        "schema_version": 1, "kind": KIND + "_source_capsule", "files": rows,
        "previous_contract": prep.binding(source / "runtime-contract.json"), "native_jobs_launched": 0})
    return {"status": "sources_archived", "files": len(rows), "native_jobs_launched": 0}


def _history(capsule_binding):
    path = prep._reopen(capsule_binding)
    capsule = prep._json(path)
    if (set(capsule) != {"schema_version", "kind", "previous_contract", "files", "native_jobs_launched"}
            or capsule["schema_version"] != 1 or capsule["kind"] != KIND + "_source_capsule"
            or type(capsule["native_jobs_launched"]) is not int or capsule["native_jobs_launched"] != 0):
        raise ValueError("Exact serial-H2 continuation source archive required")
    old = previous.previous._sealed(capsule["previous_contract"], previous.KIND)
    if prep.binding(path.parent / "runtime-contract.json")["sha256"] != capsule["previous_contract"]["sha256"]:
        raise ValueError("Archived execution contract differs")
    previous.previous._files(old["implementation"], capsule["files"], path.parent)
    for row in old["implementation"]:
        if Path(row["path"]).relative_to(ROOT).as_posix() not in CHANGED:
            prep._reopen(row)
    _, work, original, settings, control = previous._history(old["plan"]["source_capsule"])
    control = NativeExecutionControl(control.root, control.contract,
                                    execution_amendment=old["plan"]["execution_amendment"])
    if old["plan"]["inputs"] != control.validate_execution_amendment()["allowed_inputs"]:
        raise ValueError("Prior matched-energy bindings differ")
    for row in old["plan"]["inputs"]:
        prep._reopen(row)
    return old, work, original, settings, control


def _inventory():
    return sorted([*previous._inventory(), prep.binding(Path(__file__))], key=lambda r: r["path"])


def _verify_inputs(plan, old, original, control):
    parent = old["plan"]
    folders = {**parent["energy_directories"], "TS": str(Path(parent["energy_directories"]["TS"]) / "memory-retry-001")}
    memory = {label: 1500 if label == "TS" else 1000 for label in EXECUTION.LABELS}
    amendment = control.validate_execution_amendment()
    if (amendment["previous_amendment"] != parent["execution_amendment"]
            or plan["energy_directories"] != folders or plan["energy_ncores"] != parent["energy_ncores"]
            or plan["energy_maxcore_mb_per_rank"] != memory or plan["inputs"] != amendment["allowed_inputs"]
            or amendment["memory_input"] != plan["inputs"][4]):
        raise ValueError("TS memory exception escaped its exact input")
    route = EXECUTION.prediction_to_route(EXECUTION._prediction(prep._json(prep._reopen(original["plan"]["prediction"]))))
    points, _ = EXECUTION._verified_points(EXECUTION._retained_verified_path(plan["retained_verified_path"], route, control), route)
    for label, row in zip(EXECUTION.LABELS, plan["inputs"]):
        prep._reopen(row)
        kwargs = {"ncores": 1} if label == "p1" else {"maxcore_mb_per_rank": 1500} if label == "TS" else {}
        inp, _ = EXECUTION._single_point_input(points[label], Path(folders[label]), original["plan"], **kwargs)
        if prep.binding(inp) != row:
            raise ValueError("Frozen matched-energy input changed")
        EXECUTION._preflight_energy_parallelism(inp)


def prepare(capsule, destination, *, authorized=False):
    if not authorized:
        raise ValueError("This TS-only memory/budget continuation requires authorization")
    capsule, destination = Path(capsule).resolve(), Path(destination).resolve()
    if capsule.is_dir():
        capsule /= "source-capsule.json"
    old, work, original, _, control = _history(prep.binding(capsule))
    snapshot = control.snapshot()
    if snapshot["pending"] or snapshot["counts"] != {"optimization_path": 9, "frequency": 7, "irc": 2, "single_point": 11}:
        raise ValueError("Expected closed 29-attempt checkpoint required")
    if destination.exists():
        raise ValueError("Prepare a fresh TS-memory continuation")
    destination.mkdir(parents=True)
    (destination / "campaign").mkdir()
    parent = old["plan"]
    path_copy = destination / "verified-path.json"
    prep._publish(path_copy, prep._json(prep._reopen(parent["retained_verified_path"]["path_summary"])))
    route = EXECUTION.prediction_to_route(EXECUTION._prediction(prep._json(prep._reopen(original["plan"]["prediction"]))))
    points, _ = EXECUTION._verified_points({"path": path_copy}, route)
    folders = {**parent["energy_directories"], "TS": str(work / "energies/TS/memory-retry-001")}
    inputs = [prep.binding(EXECUTION._single_point_input(points[label], Path(folders[label]), original["plan"],
        **({"ncores": 1} if label == "p1" else {"maxcore_mb_per_rank": 1500} if label == "TS" else {}))[0])
        for label in EXECUTION.LABELS]
    failed = work / "energies/TS/energy.inp"
    amendment = {**prep._json(prep._reopen(parent["execution_amendment"])),
        "kind": "storca_ts_memory_single_point_amendment", "authorized": True,
        "original_single_point_cap": 11, "revised_single_point_cap": 12,
        "previous_amendment": parent["execution_amendment"], "memory_input": inputs[4], "maxcore_mb_per_rank": 1500,
        "journal_prefix": [prep.binding(p) for p in sorted((control.root / "events").iterdir())],
        "prefix_spending": {k: snapshot[k] for k in ("counts", "wall_seconds", "core_seconds")},
        "allowed_inputs": inputs, "failed_input": prep.binding(failed), "failed_output": prep.binding(failed.with_suffix(".out"))}
    prep._publish(destination / "execution-amendment.json", amendment)
    plan = {"schema_version": 1, "kind": KIND, "authorized": True, "source_capsule": prep.binding(capsule),
        "execution_amendment": prep.binding(destination / "execution-amendment.json"),
        "retained_verified_path": {"path_summary": prep.binding(path_copy), "handoff": parent["retained_verified_path"]["handoff"]},
        "energy_directories": folders, "energy_ncores": parent["energy_ncores"], "inputs": inputs,
        "energy_maxcore_mb_per_rank": {label: 1500 if label == "TS" else 1000 for label in EXECUTION.LABELS}}
    prep._publish(destination / "runtime-plan.json", plan)
    prep._publish(destination / "runtime-contract.json", {"schema_version": 1, "kind": KIND + "_sealed",
        "plan": plan, "plan_source": prep.binding(destination / "runtime-plan.json"), "implementation": _inventory()})
    return inspect(destination)


def reopen(destination):
    destination = Path(destination).resolve()
    sealed = previous.previous._sealed(prep.binding(destination / "runtime-contract.json"), KIND)
    for row in sealed["implementation"]:
        prep._reopen(row)
    if sealed["implementation"] != _inventory():
        raise ValueError("TS-memory continuation implementation changed")
    plan = sealed["plan"]
    if (set(plan) != {"schema_version", "kind", "authorized", "source_capsule", "execution_amendment",
                     "retained_verified_path", "energy_directories", "energy_ncores", "energy_maxcore_mb_per_rank", "inputs"}
            or plan["kind"] != KIND or plan["schema_version"] != 1 or plan["authorized"] is not True):
        raise ValueError("Exact authorized TS-memory plan required")
    old, _, original, settings, control = _history(plan["source_capsule"])
    control = NativeExecutionControl(control.root, control.contract, execution_amendment=plan["execution_amendment"])
    _verify_inputs(plan, old, original, control)
    return sealed, original, settings, control


def inspect(destination):
    sealed, _, _, control = reopen(destination)
    snapshot = control.snapshot()
    return {"status": "ts_memory_energy_continuation_ready", "authorized": sealed["plan"]["authorized"],
        "counts": snapshot["counts"], "pending": len(snapshot["pending"]), "native_jobs_launched": 0,
        "contract": prep.binding(Path(destination).resolve() / "runtime-contract.json"), **EXECUTION.FLAGS}


def run(destination, *, rmg_env):
    destination = Path(destination).resolve()
    sealed, original, settings, control = reopen(destination)
    plan = sealed["plan"]
    from storca.rate_transfer_runtime import runtime_environment
    runtime = runtime_environment(original["plan"]["runtime_inventory"]) if original["plan"]["schema_version"] == 2 else nullcontext()
    with runtime:
        (destination / "campaign").mkdir(exist_ok=True)
        return EXECUTION._execute_validated_context(destination / "campaign", plan=original["plan"], settings=settings,
            rmg_env=rmg_env, control_root=control.root, context_contract_path=destination / "runtime-contract.json",
            validate_context=lambda: reopen(destination), execution_amendment=plan["execution_amendment"],
            retained_verified_path=plan["retained_verified_path"], energy_directories=plan["energy_directories"],
            energy_ncores=plan["energy_ncores"], energy_maxcore_mb_per_rank=plan["energy_maxcore_mb_per_rank"])


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
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
