"""Explicit source-bound correction of an interrupted H + methane resume.

Keeps the original ledger, workspace, decisions and native artifacts. The old
implementation is archived, not relabelled as the new implementation. This
campaign-specific wrapper grants no extra jobs or physical-model permissions.
"""
import argparse
from contextlib import nullcontext
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import run_h_ch4_partial_neb_recovery as recovery
from storca import rate_transfer_preparation as prep
from storca.scientific_panel import canonical_hash

KIND = "h_ch4_remaining_work_resume"
CHANGED_FILES = {"storca/route_search_execution.py", "storca/reaction_path.py",
                 "tests/run_h_ch4_partial_neb_recovery.py"}


def archive(workdir, destination):
    """Run BEFORE editing bound sources; no native work or journal changes."""
    workdir, destination = Path(workdir).resolve(), Path(destination).resolve()
    sealed = recovery.reopen_recovery(workdir)
    if destination.exists():
        raise ValueError("Archive requires a fresh destination")
    destination.mkdir(parents=True)
    files = []
    for item in sealed["implementation"]:
        source = prep._reopen(item)
        relative = source.relative_to(ROOT)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if prep.binding(target)["sha256"] != item["sha256"]:
            raise ValueError("Source archive copy differs")
        files.append({"original_path": item["path"], "relative_path": str(relative),
                      "sha256": item["sha256"]})
    previous = prep.binding(workdir / "recovery-contract.json")
    shutil.copy2(prep._reopen(previous), destination / "recovery-contract.json")
    shutil.copy2(prep._reopen(sealed["plan_source"]), destination / "recovery-plan.json")
    capsule = {"schema_version": 1, "kind": KIND + "_source_capsule",
               "previous_contract": previous, "files": files, "native_jobs_launched": 0}
    prep._publish(destination / "source-capsule.json", capsule)
    return {"status": "original_sources_archived", "capsule": prep.binding(destination / "source-capsule.json"),
            "files": len(files), "native_jobs_launched": 0}


def _archive(previous_binding, capsule_binding):
    previous = prep._json(prep._reopen(previous_binding))
    if (set(previous) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or previous["schema_version"] != 1 or previous["kind"] != recovery.KIND + "_sealed"
            or previous["plan"]["authorized"] is not True
            or canonical_hash(prep._json(prep._reopen(previous["plan_source"]))) != canonical_hash(previous["plan"])):
        raise ValueError("Exact authorized prior recovery required")
    capsule_path = prep._reopen(capsule_binding)
    capsule = prep._json(capsule_path)
    if (set(capsule) != {"schema_version", "kind", "previous_contract", "files", "native_jobs_launched"}
            or capsule["schema_version"] != 1 or capsule["kind"] != KIND + "_source_capsule"
            or capsule["previous_contract"] != previous_binding
            or type(capsule["native_jobs_launched"]) is not int or capsule["native_jobs_launched"] != 0):
        raise ValueError("Exact recovery source capsule required")
    if (prep.binding(capsule_path.parent / "recovery-contract.json")["sha256"] != previous_binding["sha256"]
            or prep.binding(capsule_path.parent / "recovery-plan.json")["sha256"] != previous["plan_source"]["sha256"]):
        raise ValueError("Archived prior contract or plan changed")
    required = {row["path"]: row["sha256"] for row in previous["implementation"]}
    if len(required) != len(previous["implementation"]):
        raise ValueError("Duplicate prior implementation source")
    found, archived, relatives = {}, {}, set()
    for row in capsule["files"]:
        if set(row) != {"original_path", "relative_path", "sha256"}:
            raise ValueError("Malformed recovery source capsule")
        relative = Path(row["relative_path"])
        if (relative.is_absolute() or ".." in relative.parts or str(relative) in relatives
                or row["original_path"] in found or required.get(row["original_path"]) != row["sha256"]):
            raise ValueError("Recovery archive inventory differs")
        source = capsule_path.parent / relative
        if source.is_symlink() or not source.resolve().is_relative_to(capsule_path.parent):
            raise ValueError("Recovery source escaped archive")
        archived[row["original_path"]] = prep.binding(prep._reopen({"path": str(source), "sha256": row["sha256"]}))
        found[row["original_path"]] = row["sha256"]
        relatives.add(str(relative))
    if found != required:
        raise ValueError("Incomplete recovery source archive")
    # Only these concrete implementation fixes may differ from the old run.
    for original, digest in required.items():
        if Path(original).relative_to(ROOT).as_posix() not in CHANGED_FILES:
            prep._reopen({"path": original, "sha256": digest})
    return previous, archived


def _inventory():
    rows = recovery._implementation_bindings()
    return sorted([*rows, prep.binding(Path(__file__))], key=lambda row: row["path"])


def _context(plan):
    previous, archived = _archive(plan["previous_contract"], plan["source_capsule"])
    workdir = Path(plan["previous_contract"]["path"]).parent
    prior, original, settings, control, request = recovery._validate_plan(
        previous["plan"], workdir, diagnostic_implementation=previous["implementation"])
    # The old scientific decision identities remain intact, but the outer
    # revision binds the exact NEW executing sources and the exact OLD archive.
    planner_path = workdir / "path/route-search/planner/contract.json"
    if plan["planner_contract"] != prep.binding(planner_path):
        raise ValueError("Prior planner contract changed")
    planner = prep._json(planner_path)
    sources = [{"original": row, "archived": archived[row["path"]]}
               for row in planner["context"]["implementation"]]
    recovery._prefix(plan, control)
    return workdir, prior, original, settings, control, request, {
        "contract": plan["planner_contract"], "implementation_sources": sources}


def prepare(workdir, capsule, destination, *, authorized=False):
    workdir, destination = Path(workdir).resolve(), Path(destination).resolve()
    capsule = Path(capsule).resolve()
    if capsule.is_dir():
        capsule /= "source-capsule.json"
    previous = prep.binding(workdir / "recovery-contract.json")
    old, _ = _archive(previous, prep.binding(capsule))
    _, _, _, control, _ = recovery._validate_plan(old["plan"], workdir,
                                                diagnostic_implementation=old["implementation"])
    snapshot = control.snapshot()
    if snapshot["pending"]:
        raise ValueError("Owned native work is still pending; do not duplicate")
    if destination.exists():
        raise ValueError("Resume revision requires a fresh destination")
    if destination == workdir or destination.is_relative_to(workdir):
        raise ValueError("Keep the new revision outside historical campaign evidence")
    plan = {"schema_version": 1, "kind": KIND, "authorized": authorized,
            "previous_contract": previous, "source_capsule": prep.binding(capsule),
            "planner_contract": prep.binding(workdir / "path/route-search/planner/contract.json"),
            "journal_prefix": [prep.binding(p) for p in sorted((control.root / "events").iterdir())],
            "prefix_spending": {k: snapshot[k] for k in ("counts", "wall_seconds", "core_seconds")}}
    _context(plan)
    destination.mkdir(parents=True)
    prep._publish(destination / "resume-plan.json", plan)
    prep._publish(destination / "resume-contract.json", {
        "schema_version": 1, "kind": KIND + "_sealed", "plan": plan,
        "plan_source": prep.binding(destination / "resume-plan.json"), "implementation": _inventory()})
    return inspect(destination)


def reopen(destination):
    destination = Path(destination).resolve()
    sealed = prep._json(destination / "resume-contract.json")
    if (set(sealed) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or sealed["schema_version"] != 1 or sealed["kind"] != KIND + "_sealed"):
        raise ValueError("Invalid remaining-work revision contract")
    plan = sealed["plan"]
    if (set(plan) != {"schema_version", "kind", "authorized", "previous_contract", "source_capsule",
                      "planner_contract", "journal_prefix", "prefix_spending"}
            or plan["schema_version"] != 1 or plan["kind"] != KIND or type(plan["authorized"]) is not bool
            or canonical_hash(prep._json(prep._reopen(sealed["plan_source"]))) != canonical_hash(plan)):
        raise ValueError("Frozen remaining-work revision changed")
    for row in sealed["implementation"]:
        prep._reopen(row)
    if sealed["implementation"] != _inventory():
        raise ValueError("Remaining-work implementation inventory changed")
    _context(plan)
    return sealed


def inspect(destination):
    sealed = reopen(destination)
    _, _, _, _, control, _, _ = _context(sealed["plan"])
    snapshot = control.snapshot()
    return {"status": "remaining_work_revision_reopened", "authorized": sealed["plan"]["authorized"],
            "contract": prep.binding(Path(destination).resolve() / "resume-contract.json"),
            "counts": snapshot["counts"], "pending": len(snapshot["pending"]),
            "native_jobs_launched": 0, **recovery.execution.FLAGS}


def run(destination, *, rmg_env):
    destination = Path(destination).resolve()
    sealed = reopen(destination)
    if sealed["plan"]["authorized"] is not True:
        raise ValueError("Remaining-work revision has not been authorized")
    workdir, prior, original, settings, control, request, planner = _context(sealed["plan"])
    plan = original["plan"]
    from storca.rate_transfer_runtime import runtime_environment
    runtime = runtime_environment(plan["runtime_inventory"]) if plan["schema_version"] == 2 else nullcontext()
    fragments = prep._json(prep._reopen(prior["plan"]["retained_fragments"]))
    with runtime:
        return recovery.execution._execute_validated_context(
            workdir, plan=plan, settings=settings, rmg_env=rmg_env, control_root=control.root,
            context_contract_path=destination / "resume-contract.json",
            validate_context=lambda: reopen(destination), path_options={
                "route_search_policy": recovery.continuation._policy(prior["plan"]["search_policy"]),
                "retained_fragment_request": {"schema_version": 1,
                    "workdir": fragments["stationary_points_root"], "sources": fragments["sources"]},
                "retained_partial_neb_request": request, "retained_planner_context": planner})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    a = sub.add_parser("archive")
    a.add_argument("workdir", type=Path)
    a.add_argument("destination", type=Path)
    p = sub.add_parser("prepare")
    p.add_argument("workdir", type=Path)
    p.add_argument("capsule", type=Path)
    p.add_argument("destination", type=Path)
    p.add_argument("--authorized", action="store_true")
    sub.add_parser("inspect").add_argument("destination", type=Path)
    r = sub.add_parser("run")
    r.add_argument("destination", type=Path)
    r.add_argument("--rmg-env", required=True)
    args = parser.parse_args()
    if args.command == "archive":
        result = archive(args.workdir, args.destination)
    elif args.command == "prepare":
        result = prepare(args.workdir, args.capsule, args.destination, authorized=args.authorized)
    elif args.command == "inspect":
        result = inspect(args.destination)
    else:
        result = run(args.destination, rmg_env=args.rmg_env)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
