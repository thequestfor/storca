"""One authorized H2 serial-energy exception; reuse the original H/CH4 ledger."""
import argparse
from contextlib import nullcontext
from pathlib import Path
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import resume_h_ch4_pmix_runtime as previous
from storca import rate_transfer_preparation as prep
from storca.rate_transfer_control import NativeExecutionControl

KIND = "h_ch4_h2_serial_energy_continuation"
CHANGED = {"storca/rate_transfer_control.py", "storca/rate_transfer_execution.py"}
EXECUTION = previous.previous.recovery.execution


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
        "schema_version": 1, "kind": KIND + "_source_capsule",
        "previous_contract": prep.binding(source / "runtime-contract.json"),
        "files": rows, "native_jobs_launched": 0})
    return {"status": "sources_archived", "files": len(rows), "native_jobs_launched": 0}


def _history(capsule_binding):
    path = prep._reopen(capsule_binding)
    capsule = prep._json(path)
    if (set(capsule) != {"schema_version", "kind", "previous_contract", "files", "native_jobs_launched"}
            or capsule["schema_version"] != 1 or capsule["kind"] != KIND + "_source_capsule"
            or type(capsule["native_jobs_launched"]) is not int or capsule["native_jobs_launched"] != 0):
        raise ValueError("Exact prior PMIx continuation source archive required")
    old = previous._sealed(capsule["previous_contract"], previous.KIND)
    if prep.binding(path.parent / "runtime-contract.json")["sha256"] != capsule["previous_contract"]["sha256"]:
        raise ValueError("Archived execution contract differs")
    previous._files(old["implementation"], capsule["files"], path.parent)
    for row in old["implementation"]:
        if Path(row["path"]).relative_to(ROOT).as_posix() not in CHANGED:
            prep._reopen(row)
    work, original, settings, control = previous._history(old["plan"]["source_capsule"])
    amended = NativeExecutionControl(control.root, control.contract,
                                    execution_amendment=old["plan"]["execution_amendment"])
    amendment = amended.validate_execution_amendment()
    if old["plan"]["inputs"] != amendment["allowed_inputs"]:
        raise ValueError("Prior matched-energy bindings differ")
    for source in old["plan"]["inputs"]:
        prep._reopen(source)
    route = EXECUTION.prediction_to_route(EXECUTION._prediction(prep._json(prep._reopen(original["plan"]["prediction"]))))
    EXECUTION._retained_verified_path(old["plan"]["retained_verified_path"], route, amended)
    return old, work, original, settings, control


def _inventory():
    return sorted([*previous._inventory(), prep.binding(Path(__file__))], key=lambda row: row["path"])


def prepare(capsule, destination, *, authorized=False):
    if not authorized:
        raise ValueError("This specific serial H2/budget continuation requires authorization")
    capsule, destination = Path(capsule).resolve(), Path(destination).resolve()
    if capsule.is_dir():
        capsule /= "source-capsule.json"
    old, work, original, _, control = _history(prep.binding(capsule))
    snapshot = control.snapshot()
    if snapshot["pending"] or snapshot["counts"] != {"optimization_path": 9, "frequency": 7, "irc": 2, "single_point": 8}:
        raise ValueError("Expected closed 26-attempt checkpoint required")
    if destination.exists():
        raise ValueError("Prepare a fresh serial-energy continuation")
    destination.mkdir(parents=True)
    (destination / "campaign").mkdir()
    path_copy = destination / "verified-path.json"
    prep._publish(path_copy, prep._json(prep._reopen(old["plan"]["retained_verified_path"]["path_summary"])))
    route = EXECUTION.prediction_to_route(EXECUTION._prediction(prep._json(prep._reopen(original["plan"]["prediction"]))))
    points, _ = EXECUTION._verified_points({"path": path_copy}, route)
    folders = {**old["plan"]["energy_directories"], "p1": str(work / "energies/p1/serial-retry-001")}
    cores = {label: (1 if label == "p1" else 4) for label in EXECUTION.LABELS}
    inputs = [prep.binding(EXECUTION._single_point_input(points[label], Path(folders[label]), original["plan"],
                  **({"ncores": 1} if label == "p1" else {}))[0]) for label in EXECUTION.LABELS]
    for row in inputs:
        EXECUTION._preflight_energy_parallelism(Path(row["path"]))
    failed = work / "energies/p1/energy.inp"
    parent = prep._json(prep._reopen(old["plan"]["execution_amendment"]))
    amendment = {**parent, "kind": "storca_h2_serial_single_point_amendment",
        "authorized": True, "original_single_point_cap": 10, "revised_single_point_cap": 11,
        "previous_amendment": old["plan"]["execution_amendment"], "serial_input": inputs[2],
        "journal_prefix": [prep.binding(p) for p in sorted((control.root / "events").iterdir())],
        "prefix_spending": {k: snapshot[k] for k in ("counts", "wall_seconds", "core_seconds")},
        "allowed_inputs": inputs, "failed_input": prep.binding(failed),
        "failed_output": prep.binding(failed.with_suffix(".out"))}
    prep._publish(destination / "execution-amendment.json", amendment)
    plan = {"schema_version": 1, "kind": KIND, "authorized": True,
        "source_capsule": prep.binding(capsule), "execution_amendment": prep.binding(destination / "execution-amendment.json"),
        "retained_verified_path": {"path_summary": prep.binding(path_copy),
                                   "handoff": old["plan"]["retained_verified_path"]["handoff"]},
        "energy_directories": folders, "energy_ncores": cores, "inputs": inputs}
    prep._publish(destination / "runtime-plan.json", plan)
    prep._publish(destination / "runtime-contract.json", {"schema_version": 1, "kind": KIND + "_sealed",
        "plan": plan, "plan_source": prep.binding(destination / "runtime-plan.json"), "implementation": _inventory()})
    return inspect(destination)


def reopen(destination):
    destination = Path(destination).resolve()
    sealed = previous._sealed(prep.binding(destination / "runtime-contract.json"), KIND)
    for row in sealed["implementation"]:
        prep._reopen(row)
    if sealed["implementation"] != _inventory():
        raise ValueError("Serial-energy continuation implementation changed")
    plan = sealed["plan"]
    if (set(plan) != {"schema_version", "kind", "authorized", "source_capsule", "execution_amendment",
                     "retained_verified_path", "energy_directories", "energy_ncores", "inputs"}
            or plan["schema_version"] != 1 or plan["kind"] != KIND or plan["authorized"] is not True):
        raise ValueError("Exact serial-energy continuation plan required")
    old, work, original, settings, control = _history(plan["source_capsule"])
    control = NativeExecutionControl(control.root, control.contract, execution_amendment=plan["execution_amendment"])
    amendment = control.validate_execution_amendment()
    if amendment["previous_amendment"] != old["plan"]["execution_amendment"]:
        raise ValueError("Serial-energy amendment predecessor differs")
    folders = {**old["plan"]["energy_directories"], "p1": str(work / "energies/p1/serial-retry-001")}
    cores = {label: (1 if label == "p1" else 4) for label in EXECUTION.LABELS}
    if (plan["energy_directories"] != folders or plan["energy_ncores"] != cores
            or plan["inputs"] != amendment["allowed_inputs"]
            or amendment["serial_input"] != plan["inputs"][2]):
        raise ValueError("Serial-energy exception escaped its frozen H2 input")
    route = EXECUTION.prediction_to_route(EXECUTION._prediction(prep._json(prep._reopen(original["plan"]["prediction"]))))
    points, _ = EXECUTION._verified_points(EXECUTION._retained_verified_path(plan["retained_verified_path"], route, control), route)
    for label, row in zip(EXECUTION.LABELS, plan["inputs"]):
        prep._reopen(row)
        inp, _ = EXECUTION._single_point_input(points[label], Path(folders[label]), original["plan"],
                                              **({"ncores": 1} if label == "p1" else {}))
        if prep.binding(inp) != row:
            raise ValueError("Frozen energy input changed")
        EXECUTION._preflight_energy_parallelism(inp)
    return sealed, original, settings, control


def inspect(destination):
    sealed, _, _, control = reopen(destination)
    snapshot = control.snapshot()
    return {"status": "serial_h2_energy_continuation_ready", "authorized": sealed["plan"]["authorized"],
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
        return EXECUTION._execute_validated_context(destination / "campaign", plan=original["plan"],
            settings=settings, rmg_env=rmg_env, control_root=control.root,
            context_contract_path=destination / "runtime-contract.json", validate_context=lambda: reopen(destination),
            execution_amendment=plan["execution_amendment"], retained_verified_path=plan["retained_verified_path"],
            energy_directories=plan["energy_directories"], energy_ncores=plan["energy_ncores"])


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
