"""Explicit historical-to-live transfer continuation, with no free child ledger.

The historical execution contract stays strict and unchanged. Exact archived
source bytes establish its implementation lineage; a separately sealed live
inventory establishes the new search. Preparation is zero-dispatch and never
inherits native continuation authorization from the old campaign.
"""
from __future__ import annotations

import ast
from contextlib import nullcontext
import hashlib
from pathlib import Path
import re
from typing import Any

from src.orca_job_settings import active_native_settings, native_settings_scope
from . import rate_transfer_execution as execution
from . import rate_transfer_preparation as prep
from .rate_transfer_control import NativeExecutionControl, TransferControlError
from .scientific_panel import canonical_hash


KIND = "storca_rate_transfer_adaptive_continuation"


def _source_capsule(binding, original_contract_binding, original_contract):
    path = prep._reopen(binding)
    capsule = prep._json(path)
    if (set(capsule) != {"schema_version", "kind", "execution_contract", "files", "native_jobs_launched"}
            or type(capsule["schema_version"]) is not int or capsule["schema_version"] != 1
            or capsule["kind"] != "rate_transfer_historical_source_capsule"
            or capsule["execution_contract"] != original_contract_binding
            or type(capsule["native_jobs_launched"]) is not int or capsule["native_jobs_launched"] != 0):
        raise ValueError("Exact historical source capsule required")
    required = {item["path"]: item["sha256"] for item in original_contract["implementation"]}
    found, originals, relatives = {}, set(), set()
    for item in capsule["files"]:
        if (not isinstance(item, dict)
                or set(item) != {"original_path", "relative_path", "sha256", "originally_bound"}
                or type(item["originally_bound"]) is not bool):
            raise ValueError("Malformed historical source inventory")
        relative = Path(item["relative_path"])
        original = item["original_path"]
        if (relative.is_absolute() or ".." in relative.parts or not relative.parts
                or not Path(original).is_absolute() or original in originals
                or str(relative) in relatives):
            raise ValueError("Historical source path is unsafe or duplicated")
        archived = path.parent / relative
        if archived.is_symlink() or not archived.resolve().is_relative_to(path.parent.resolve()):
            raise ValueError("Historical source escaped its capsule")
        prep._reopen({"path": str(archived.resolve()), "sha256": item["sha256"]})
        if item["originally_bound"] != (original in required):
            raise ValueError("Original bound source inventory differs")
        if original in required:
            if item["sha256"] != required[original]:
                raise ValueError("Archived source does not match its historical binding")
            found[original] = archived
        originals.add(original)
        relatives.add(str(relative))
    if set(found) != set(required):
        raise ValueError("Historical source capsule omitted bound implementation")
    return capsule, found


def _historical_input_reader(sources):
    """Load only the exact archived pure reader; never import a mixed old package.

    The historical function reads native input plus referenced XYZ/Hessian/band
    bytes. Compiling its original AST avoids importing historical orchestration
    into live modules or replacing its digest semantics with a file-only hash.
    """
    paths = [path for path in sources.values() if path.name == "reaction_path.py"]
    if len(paths) != 1:
        raise ValueError("One archived reaction-path input reader is required")
    tree = ast.parse(paths[0].read_text())
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                 and node.name == "_orca_job_input_evidence"]
    if len(functions) != 1:
        raise ValueError("Historical native input reader is missing or ambiguous")
    isolated = ast.Module(body=[ast.ImportFrom(module="__future__", names=[
        ast.alias(name="annotations")], level=0), functions[0]], type_ignores=[])
    namespace = {"Path": Path, "re": re, "hashlib": hashlib, "Any": Any}
    exec(compile(ast.fix_missing_locations(isolated), str(paths[0]), "exec"), namespace)
    return namespace["_orca_job_input_evidence"]


def _original_context(campaign, capsule_binding):
    campaign = Path(campaign).resolve()
    contract_binding = prep.binding(campaign / "execution-contract.json")
    contract = prep._json(campaign / "execution-contract.json")
    if (set(contract) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or type(contract["schema_version"]) is not int or contract["schema_version"] != 1
            or contract["kind"] != execution.KIND + "_sealed"):
        raise ValueError("Invalid historical execution contract")
    source = prep._reopen(contract["plan_source"])
    if canonical_hash(prep._json(source)) != canonical_hash(contract["plan"]):
        raise ValueError("Historical raw plan and sealed plan differ")
    _, archived_sources = _source_capsule(capsule_binding, contract_binding, contract)
    prediction, settings = execution._validate_plan(contract["plan"])
    if contract["plan"]["control"]["authorized"] is not True:
        raise ValueError("The historical campaign was not authorized")
    control = NativeExecutionControl(campaign / "native-control", contract["plan"]["control"])
    # No initialize: missing old evidence must fail, never create a zero ledger.
    snapshot = control.snapshot()
    if snapshot["pending"]:
        raise TransferControlError("Historical native ownership is unresolved")
    settings_binding = prep.binding(control.root / "native-settings.json")
    if prep._json(Path(settings_binding["path"])) != {
            "sha256": settings.sha256, "settings": settings.as_dict()}:
        raise ValueError("Historical native settings differ")
    old_reader = _historical_input_reader(archived_sources)
    from .reaction_path import _orca_job_input_evidence
    for item in snapshot["reservations"].values():
        native_input = Path(item["input"])
        historical_digest = old_reader(native_input)[0].hexdigest()
        if historical_digest != item["input_sha256"]:
            raise ValueError("Historical native input or referenced evidence changed")
        if _orca_job_input_evidence(native_input)[0].hexdigest() != historical_digest:
            raise ValueError("Live native-input digest semantics differ from the historical reader")
        control.validate_reuse(native_input, historical_digest)
    return contract, prediction, settings, control, snapshot


def _fragment_evidence(campaign, prediction, settings):
    from .reaction_path import reuse_stationary_fragment_records
    from .route_geometry import resolve_route_atom_mapping
    route = execution.prediction_to_route(prediction)
    mapping = resolve_route_atom_mapping(route)
    if not mapping.get("valid"):
        raise ValueError("Historical reaction mapping is unresolved")
    root = Path(campaign) / "path" / "stationary-points"
    active = active_native_settings()
    if active is not None and active.sha256 != settings.sha256:
        raise ValueError("Active native settings differ from retained fragment settings")
    # Final source revalidation may run inside the same owned controller.
    # Reuse only an identical active profile; never nest or replace a different one.
    with nullcontext() if active is not None else native_settings_scope(settings):
        records = reuse_stationary_fragment_records(
            route, root, mapping, ncores=4, method_keywords=list(execution.LOW_KEYWORDS),
            retained_endpoints=None)
    if records.get("valid") is not True:
        raise ValueError("Historical independent fragment evidence did not revalidate")
    names = ("input.xyz", "endpoint-contract.json", "metadata.json", "opt.inp", "opt.out",
             "opt.xyz", "opt.contract.json", "opt.process.json", "freq.inp", "freq.out",
             "freq.hess", "freq.contract.json", "freq.process.json", "atomic-stationary.json",
             "atom-energy.inp", "atom-energy.out", "atom-energy.contract.json", "atom-energy.process.json")
    sources = {}
    for side in ("reactants", "products"):
        for row in records[side]:
            folder = Path(row["input_xyz"]).parent
            for name in names:
                path = folder / name
                if path.is_file():
                    sources[str(path)] = prep.binding(path)
    return {"schema_version": 1, "kind": KIND + "_retained_fragments",
            "stationary_points_root": str(root.resolve()), "records": records,
            "sources": [sources[key] for key in sorted(sources)],
            "native_jobs_launched": 0, "full_historical_path_admitted": False}


def _policy(value=None):
    from .route_search import RouteSearchPolicy
    policy = (RouteSearchPolicy(extra_completion_cost={"single_point": 5})
              if value is None else RouteSearchPolicy.from_dict(value))
    if policy.extra_completion_cost != {"single_point": 5}:
        raise ValueError("The original rate-transfer campaign must retain five downstream matched energies")
    return policy


def _implementation_bindings():
    root = Path(__file__).resolve().parents[1]
    sources = {item["path"]: item for item in execution._implementation_bindings()}
    for name in ("storca/rate_transfer_continuation.py", "storca/route_search.py",
                 "storca/route_search_execution.py",
                 "storca/route_search_recovery.py",
                 "storca/reaction_seed_geometry.py", "tests/run_rate_transfer_benchmark.py"):
        item = prep.binding(root / name)
        sources[item["path"]] = item
    return [sources[key] for key in sorted(sources)]


def prepare_continuation(original_campaign, source_capsule, destination):
    """Freeze a zero-native-work draft preserving the original spent budget."""
    original_campaign, destination = Path(original_campaign).resolve(), Path(destination).resolve()
    capsule_path = Path(source_capsule).resolve()
    if capsule_path.is_dir():
        capsule_path /= "source-capsule.json"
    capsule_binding = prep.binding(capsule_path)
    contract, prediction, settings, control, snapshot = _original_context(original_campaign, capsule_binding)
    if destination == original_campaign or destination.is_relative_to(original_campaign):
        raise ValueError("Continuation workspace must be outside the historical campaign")
    evidence = _fragment_evidence(original_campaign, prediction, settings)
    destination.mkdir(parents=True, exist_ok=True)
    prep._publish(destination / "retained-fragments.json", evidence)
    plan = {"schema_version": 1, "kind": KIND, "authorized": False,
            "original_campaign": str(original_campaign),
            "original_execution_contract": prep.binding(original_campaign / "execution-contract.json"),
            "source_capsule": capsule_binding,
            "native_control_contract": prep.binding(control.root / "contract.json"),
            "native_settings": prep.binding(control.root / "native-settings.json"),
            "journal_prefix": [prep.binding(path) for path in sorted((control.root / "events").iterdir())],
            "prefix_spending": {key: snapshot[key] for key in ("counts", "wall_seconds", "core_seconds")},
            "retained_fragments": prep.binding(destination / "retained-fragments.json"),
            "search_policy": _policy().as_dict(),
            "search_workdir": str(destination / "campaign" / "path"),
            "scope": "Original conditional known-reaction model; new route search only, unchanged physical gates and cumulative native caps."}
    prep._publish(destination / "continuation-plan.json", plan)
    sealed = freeze_continuation(destination / "continuation-plan.json", destination / "campaign")
    result = {"status": "adaptive_continuation_prepared_not_authorized", "authorized": False,
              "native_jobs_launched": 0, "historical_charged_attempts": sum(snapshot["counts"].values()),
              "spent": plan["prefix_spending"], "control_root": str(control.root),
              "contract": sealed["contract"], "plan": prep.binding(destination / "continuation-plan.json"),
              **execution.FLAGS}
    prep._publish(destination / "readiness.json", result)
    return result


def _validate_plan(plan, workdir):
    required = {"schema_version", "kind", "authorized", "original_campaign", "original_execution_contract",
                "source_capsule", "native_control_contract", "native_settings", "journal_prefix",
                "prefix_spending", "retained_fragments", "search_policy", "search_workdir", "scope"}
    if (set(plan) != required or type(plan["schema_version"]) is not int
            or plan["schema_version"] != 1 or plan["kind"] != KIND
            or type(plan["authorized"]) is not bool):
        raise ValueError("Invalid adaptive continuation plan")
    campaign = Path(plan["original_campaign"])
    if not campaign.is_absolute() or campaign != campaign.resolve():
        raise ValueError("Exact absolute historical campaign required")
    workdir = Path(workdir).resolve()
    if (Path(plan["search_workdir"]) != workdir / "path"
            or workdir == campaign or workdir.is_relative_to(campaign)):
        raise ValueError("Fresh separately bound search workspace required")
    prep._reopen(plan["original_execution_contract"])
    if plan["original_execution_contract"] != prep.binding(campaign / "execution-contract.json"):
        raise ValueError("Historical contract belongs to a different campaign")
    contract, prediction, settings, control, snapshot = _original_context(campaign, plan["source_capsule"])
    for key, path in (("native_control_contract", control.root / "contract.json"),
                      ("native_settings", control.root / "native-settings.json")):
        prep._reopen(plan[key])
        if plan[key] != prep.binding(path):
            raise ValueError("Continuation cannot replace the original native ledger/settings")
    current = sorted((control.root / "events").iterdir())
    prefix = plan["journal_prefix"]
    if (not isinstance(prefix, list) or len(current) < len(prefix)
            or [prep.binding(path) for path in current[:len(prefix)]] != prefix):
        raise ValueError("Original immutable native journal prefix changed")
    # Validate the prefix's own chain/counts as well as the current full chain.
    counts = {key: 0 for key in control.contract["limits"] if key in snapshot["counts"]}
    wall, cores, ranks = 0, 0, {}
    for item in prefix:
        row = prep._json(prep._reopen(item))
        if row["event"] == "reserved":
            counts[row["category"]] += 1
            ranks[row["attempt_id"]] = row["ncores"]
        else:
            wall += row["elapsed_seconds"]
            cores += row["elapsed_seconds"] * ranks[row["attempt_id"]]
    if canonical_hash(plan["prefix_spending"]) != canonical_hash({
            "counts": counts, "wall_seconds": wall, "core_seconds": cores}):
        raise ValueError("Historical charges cannot be reset or refunded")
    fragments = prep._json(prep._reopen(plan["retained_fragments"]))
    for binding in fragments["sources"]:
        prep._reopen(binding)
    fresh = _fragment_evidence(campaign, prediction, settings)
    if canonical_hash(fragments) != canonical_hash(fresh):
        raise ValueError("Retained fragment evidence or its scientific assessment changed")
    _policy(plan["search_policy"])
    return contract, settings, control


def freeze_continuation(plan_path, workdir):
    """Seal only explicit consent; never modify the historical plan or ledger."""
    plan_path, workdir = Path(plan_path).resolve(), Path(workdir).resolve()
    plan = prep._json(plan_path)
    _validate_plan(plan, workdir)
    target = workdir / "continuation-contract.json"
    if workdir.exists() and not target.exists() and (not workdir.is_dir() or any(workdir.iterdir())):
        raise ValueError("Freeze into a fresh continuation workspace")
    workdir.mkdir(parents=True, exist_ok=True)
    contract = {"schema_version": 1, "kind": KIND + "_sealed", "plan": plan,
                "plan_source": prep.binding(plan_path), "implementation": _implementation_bindings()}
    prep._publish(target, contract)
    return {"status": "continuation_frozen_not_run", "authorized": plan["authorized"],
            "contract": prep.binding(target), "native_jobs_launched": 0, **execution.FLAGS}


def reopen_continuation(workdir):
    workdir = Path(workdir).resolve()
    contract = prep._json(workdir / "continuation-contract.json")
    if (set(contract) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or type(contract["schema_version"]) is not int
            or contract["schema_version"] != 1 or contract["kind"] != KIND + "_sealed"):
        raise ValueError("Invalid sealed continuation contract")
    source = prep._reopen(contract["plan_source"])
    if canonical_hash(prep._json(source)) != canonical_hash(contract["plan"]):
        raise ValueError("Frozen continuation plan changed")
    for binding in contract["implementation"]:
        prep._reopen(binding)
    if contract["implementation"] != _implementation_bindings():
        raise ValueError("Continuation implementation inventory changed")
    _validate_plan(contract["plan"], workdir)
    return contract


def run_continuation(workdir, *, rmg_env):
    """Run/resume separately authorized search under the SAME native controller."""
    workdir = Path(workdir).resolve()
    sealed = reopen_continuation(workdir)
    plan = sealed["plan"]
    if plan["authorized"] is not True:
        raise ValueError("This continuation has not been separately authorized")
    historical, settings, control = _validate_plan(plan, workdir)
    original_plan = historical["plan"]
    from .rate_transfer_runtime import runtime_environment
    from contextlib import nullcontext
    runtime = (runtime_environment(original_plan["runtime_inventory"])
               if original_plan["schema_version"] == 2 else nullcontext())
    with runtime:
        retained = prep._json(prep._reopen(plan["retained_fragments"]))
        return execution._execute_validated_context(
            workdir, plan=original_plan, settings=settings, rmg_env=rmg_env,
            control_root=control.root, context_contract_path=workdir / "continuation-contract.json",
            validate_context=lambda: reopen_continuation(workdir),
            path_options={"route_search_policy": _policy(plan["search_policy"]),
                          "retained_fragment_request": {
                              "schema_version": 1,
                              "workdir": str(Path(plan["original_campaign"]) / "path" / "stationary-points"),
                              "sources": retained["sources"]}})


__all__ = ["prepare_continuation", "freeze_continuation", "reopen_continuation", "run_continuation"]
