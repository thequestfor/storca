"""Same-ledger, source-bound TS-only recovery of the paid H + methane band.

This wrapper does not weaken the historical continuation's live-source checks.
It establishes a new, explicit execution lineage from the exact archived prior
implementation. Preparation and inspection are always zero-calculation reads.
"""
from __future__ import annotations

import argparse
from contextlib import nullcontext
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from storca import rate_transfer_continuation as continuation
from storca import rate_transfer_execution as execution
from storca import rate_transfer_preparation as prep
from storca.partial_neb_recovery import inspect_partial_neb_recovery
from storca.scientific_panel import canonical_hash


KIND = "storca_h_ch4_owned_partial_neb_recovery"
NUMERICAL_ENERGY_EQUIVALENCE_EH = 1e-8
NUMERICAL_SPIN_EQUIVALENCE = 1e-4
NUMERICAL_ELECTRON_INTEGRATION_TOLERANCE = 1e-5
DOUBLET_SPIN_DISTANCE_LIMIT = 0.05
SCF_RECOVERY_REASON = ("Explicit TRAH bypasses native DIIS path that aborted; matched electronic "
                       "evidence conditional, original cause unresolved")


def _archive(prior_binding, capsule_binding):
    prior_path = prep._reopen(prior_binding)
    prior = prep._json(prior_path)
    if (set(prior) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or prior["schema_version"] != 1 or prior["kind"] != continuation.KIND + "_sealed"):
        raise ValueError("A sealed historical adaptive continuation is required")
    source = prep._json(prep._reopen(prior["plan_source"]))
    if canonical_hash(source) != canonical_hash(prior["plan"]) or source["authorized"] is not True:
        raise ValueError("Historical raw/sealed authorization differs")
    capsule_path = prep._reopen(capsule_binding)
    capsule = prep._json(capsule_path)
    if (set(capsule) != {"schema_version", "kind", "continuation_contract", "files", "native_jobs_launched"}
            or capsule["schema_version"] != 1
            or capsule["kind"] != "rate_transfer_continuation_source_capsule"
            or capsule["continuation_contract"] != prior_binding
            or type(capsule["native_jobs_launched"]) is not int or capsule["native_jobs_launched"] != 0):
        raise ValueError("Exact prior-continuation source capsule required")
    archived_contract = capsule_path.parent / "continuation-contract.json"
    if prep.binding(archived_contract)["sha256"] != prior_binding["sha256"]:
        raise ValueError("Archived continuation contract changed")
    required = {item["path"]: item["sha256"] for item in prior["implementation"]}
    if len(required) != len(prior["implementation"]):
        raise ValueError("Duplicated historical implementation binding")
    found, relatives = {}, set()
    for item in capsule["files"]:
        if set(item) != {"original_path", "relative_path", "sha256", "originally_bound"}:
            raise ValueError("Malformed prior source inventory")
        relative = Path(item["relative_path"])
        original = item["original_path"]
        if (relative.is_absolute() or ".." in relative.parts or not relative.parts
                or str(relative) in relatives or original in found or not Path(original).is_absolute()
                or item["originally_bound"] is not True or required.get(original) != item["sha256"]):
            raise ValueError("Prior source inventory differs or contains unsafe paths")
        archived = capsule_path.parent / relative
        if archived.is_symlink() or not archived.resolve().is_relative_to(capsule_path.parent.resolve()):
            raise ValueError("Prior archived source escaped its capsule")
        prep._reopen({"path": str(archived.resolve()), "sha256": item["sha256"]})
        found[original] = item["sha256"]
        relatives.add(str(relative))
    if found != required:
        raise ValueError("Prior implementation source capsule is incomplete")
    # Validate the original model, settings, immutable journal and raw sources;
    # never call the prior strict live-source reopen against edited live code.
    original, settings, control = continuation._validate_plan(prior["plan"], prior_path.parent)
    prediction = execution._prediction(prep._json(prep._reopen(original["plan"]["prediction"])))
    if prediction["case_id"] != "h_ch4":
        raise ValueError("This recovery wrapper supports the locked H + methane benchmark only")
    return prior, original, settings, control


def _implementation_bindings():
    sources = {item["path"]: item for item in continuation._implementation_bindings()}
    for relative in ("storca/partial_neb_recovery.py", "tests/run_h_ch4_partial_neb_recovery.py",
                     "tests/run_h_ch4_scf_diagnostic.py"):
        item = prep.binding(ROOT / relative)
        sources[item["path"]] = item
    return [sources[key] for key in sorted(sources)]


def _prefix(plan, control):
    snapshot = control.snapshot()
    current = sorted((control.root / "events").iterdir())
    prefix = plan["journal_prefix"]
    if (not isinstance(prefix, list) or len(current) < len(prefix)
            or [prep.binding(path) for path in current[:len(prefix)]] != prefix):
        raise ValueError("The original native journal prefix changed")
    counts = {key: 0 for key in snapshot["counts"]}
    wall, cores, ranks = 0, 0, {}
    for item in prefix:
        event = prep._json(prep._reopen(item))
        if event["event"] == "reserved":
            counts[event["category"]] += 1
            ranks[event["attempt_id"]] = event["ncores"]
        else:
            wall += event["elapsed_seconds"]
            cores += event["elapsed_seconds"] * ranks[event["attempt_id"]]
    if canonical_hash(plan["prefix_spending"]) != canonical_hash({
            "counts": counts, "wall_seconds": wall, "core_seconds": cores}):
        raise ValueError("Historical paid work cannot be reset or refunded")


def _diagnostic_selection(report_binding, control, prior, prior_capsule_binding, *, implementation_inventory=None):
    # The comparison adapter reopens both owned raw electronic solves and
    # reruns its fixed, non-experimental selection policy; no success-label trust.
    from run_h_ch4_scf_diagnostic import read_comparison
    path = prep._reopen(report_binding)
    if path.name != "comparison.json":
        raise ValueError("The immutable two-arm comparison report is required")
    result = read_comparison(path.parent)
    if result["plan"]["native_control_root"] != str(control.root):
        raise ValueError("Electronic comparison belongs to another native ledger")
    capsule = result["plan"].get("prior_source_capsule")
    if capsule is not None:
        if capsule != prior_capsule_binding:
            raise ValueError("Electronic comparison belongs to a different archived continuation")
        # This explicit diagnostic revision was run after the recovery changes.
        # Its full live implementation must still be exactly the sealed one.
        # A separately bound code-repair revision may supply the exact prior
        # inventory only AFTER validating every archived byte against its old
        # sealed recovery contract. Legacy reopening still requires live bytes.
        required = {item["path"]: item for item in (
            _implementation_bindings() if implementation_inventory is None else implementation_inventory)}
    else:
        # Legacy diagnostics may be reused only with the exact old source
        # capsule and their still-exact diagnostic reader, never a hash waiver.
        required = {item["path"]: item for item in prior["implementation"]}
        diagnostic = prep.binding(ROOT / "tests" / "run_h_ch4_scf_diagnostic.py")
        required[diagnostic["path"]] = diagnostic
    historical = result["historical_implementation"]
    if (len(historical) != len(required)
            or {item["path"]: item for item in historical} != required):
        raise ValueError("Electronic comparison's exact archived/current source lineage differs")
    return result


def _numerical_equivalence(assessments):
    """Numerical equivalence screens only; not a rate-accuracy calibration."""
    if len(assessments) != 2:
        raise ValueError("Both original and TRAH raw electronic solves are required")
    for item in assessments:
        if (item.get("status") != "electronic_solve_completed_not_TS_validation"
                or item.get("normal_termination") is not True or item.get("scf_converged") is not True
                or item.get("charge") != 0 or item.get("multiplicity") != 2
                or item.get("total_electrons") != 11):
            raise ValueError("Converged neutral-doublet electronic state evidence is required from both arms")
        for name, expected in (("alpha_electrons", 6), ("beta_electrons", 5)):
            value = item.get(name)
            if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
                    or abs(value - expected) > NUMERICAL_ELECTRON_INTEGRATION_TOLERANCE):
                raise ValueError("The frozen electronic population screen failed")
        for name in ("energy_hartree", "S_squared"):
            value = item.get(name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError("Nonfinite or missing electronic comparison quantity")
        if abs(item["S_squared"] - 0.75) > DOUBLET_SPIN_DISTANCE_LIMIT:
            raise ValueError("Electronic comparison fails the frozen doublet-spin screen")
    if abs(assessments[1]["energy_hartree"] - assessments[0]["energy_hartree"]) > NUMERICAL_ENERGY_EQUIVALENCE_EH:
        raise ValueError("Electronic energies differ beyond the frozen numerical-equivalence tolerance")
    if abs(assessments[1]["S_squared"] - assessments[0]["S_squared"]) > NUMERICAL_SPIN_EQUIVALENCE:
        raise ValueError("Electronic spins differ beyond the frozen numerical-equivalence tolerance")


def _validate_request(request_binding, report_binding, prior_path, prior, original, control, prior_capsule_binding,
                      *, diagnostic_implementation=None):
    request = prep._json(prep._reopen(request_binding))
    required = {"schema_version", "kind", "input_path", "context_binding", "expected_context",
                "source_bindings", "selection", "scf_solver", "solver_evidence", "selection_reason"}
    if (set(request) != required or type(request["schema_version"]) is not int
            or request["schema_version"] != 1 or request["kind"] != "owned_partial_neb_ts_recovery_request"
            or request["selection"] not in {"converged_ci", "last_evaluated_opt"}
            or request["scf_solver"] not in {"original", "trah"}
            or request["solver_evidence"] != report_binding
            or request["selection_reason"] != SCF_RECOVERY_REASON):
        raise ValueError("Explicit candidate-only partial-NEB recovery request required")
    parent = Path(prior_path).parent / "path"
    if (not Path(request["input_path"]).resolve().is_relative_to(parent.resolve())
            or not Path(request["context_binding"]["path"]).resolve().is_relative_to(parent.resolve())):
        raise ValueError("Partial native evidence belongs to a different continuation")
    prediction = execution._prediction(prep._json(prep._reopen(original["plan"]["prediction"])))
    context = request["expected_context"]
    electronic = context.get("electronic_structure") or {}
    if (canonical_hash(context.get("route")) != canonical_hash(execution.prediction_to_route(prediction).as_dict())
            or electronic.get("method_keywords") != list(execution.LOW_KEYWORDS)
            or electronic.get("ncores") != 4
            or electronic.get("thermochemistry_temperature_K") != prediction["temperature_grid_K"][0]):
        raise ValueError("Recovery cannot alter the locked reaction, model, temperature or resources")
    comparison = _diagnostic_selection(report_binding, control, prior, prior_capsule_binding,
        **({"implementation_inventory": diagnostic_implementation} if diagnostic_implementation is not None else {}))
    if (comparison["plan"]["source_continuation"] != prep.binding(prior_path)
            or comparison["plan"]["failed_source"]["input"]["path"] != request["input_path"]
            or comparison["plan"]["failed_source"]["geometry"] != prep.binding(
                Path(request["input_path"]).with_suffix(".xyz"))):
        raise ValueError("Electronic comparison belongs to a different failed native source")
    _numerical_equivalence(comparison["assessments"])
    admission = inspect_partial_neb_recovery(request["input_path"], control=control,
        context_binding=request["context_binding"], expected_context=context,
        source_bindings=request["source_bindings"], selection=request["selection"])
    if admission.get("valid") is not True or admission.get("eligible_for_ts_recovery") is not True:
        raise ValueError("Paid partial-NEB candidate was not admitted: " + str(admission.get("failure_reason")))
    return request, admission


def _validate_plan(plan, workdir, *, diagnostic_implementation=None):
    if (set(plan) != {"schema_version", "kind", "authorized", "prior_continuation", "prior_source_capsule",
                      "partial_request", "diagnostic_report", "search_workdir", "journal_prefix", "prefix_spending",
                      "native_control_contract", "native_settings", "candidate_admission_sha256"}
            or type(plan["schema_version"]) is not int or plan["schema_version"] != 1
            or plan["kind"] != KIND or type(plan["authorized"]) is not bool):
        raise ValueError("Invalid TS-only recovery plan")
    prior, original, settings, control = _archive(plan["prior_continuation"], plan["prior_source_capsule"])
    workdir = Path(workdir).resolve()
    old = Path(plan["prior_continuation"]["path"]).parent.resolve()
    initial = Path(prior["plan"]["original_campaign"]).resolve()
    if (Path(plan["search_workdir"]) != workdir / "path" or workdir.is_relative_to(old)
            or workdir.is_relative_to(initial)):
        raise ValueError("Recovery requires a fresh separately bound workspace")
    for key, path in (("native_control_contract", control.root / "contract.json"),
                      ("native_settings", control.root / "native-settings.json")):
        if plan[key] != prep.binding(path):
            raise ValueError("Recovery must retain the original ledger and native settings")
    _prefix(plan, control)
    request, admission = _validate_request(plan["partial_request"], plan["diagnostic_report"],
                                         plan["prior_continuation"]["path"], prior, original, control,
                                         plan["prior_source_capsule"], diagnostic_implementation=diagnostic_implementation)
    if canonical_hash(admission) != plan["candidate_admission_sha256"]:
        raise ValueError("Retained candidate admission changed")
    return prior, original, settings, control, request


def prepare_recovery(prior_campaign, prior_capsule, partial_request, diagnostic_report, destination, *, authorized=False):
    """Freeze a fresh wrapper; explicit authorization never edits old contracts."""
    prior_path = Path(prior_campaign).resolve() / "continuation-contract.json"
    capsule = Path(prior_capsule).resolve()
    if capsule.is_dir():
        capsule /= "source-capsule.json"
    prior_binding, capsule_binding = prep.binding(prior_path), prep.binding(capsule)
    prior, original, settings, control = _archive(prior_binding, capsule_binding)
    request_binding, report_binding = prep.binding(partial_request), prep.binding(diagnostic_report)
    _, admission = _validate_request(request_binding, report_binding, prior_path, prior, original, control,
                                     capsule_binding)
    destination = Path(destination).resolve()
    workdir = destination / "campaign"
    snapshot = control.snapshot()
    plan = {"schema_version": 1, "kind": KIND, "authorized": authorized,
            "prior_continuation": prior_binding, "prior_source_capsule": capsule_binding,
            "partial_request": request_binding, "diagnostic_report": report_binding,
            "search_workdir": str(workdir / "path"),
            "journal_prefix": [prep.binding(path) for path in sorted((control.root / "events").iterdir())],
            "prefix_spending": {key: snapshot[key] for key in ("counts", "wall_seconds", "core_seconds")},
            "native_control_contract": prep.binding(control.root / "contract.json"),
            "native_settings": prep.binding(control.root / "native-settings.json"),
            "candidate_admission_sha256": canonical_hash(admission)}
    _validate_plan(plan, workdir)
    target = workdir / "recovery-contract.json"
    if workdir.exists() and not target.exists() and any(workdir.iterdir()):
        raise ValueError("Seal recovery only into an empty workspace")
    workdir.mkdir(parents=True, exist_ok=True)
    prep._publish(destination / "recovery-plan.json", plan)
    contract = {"schema_version": 1, "kind": KIND + "_sealed", "plan": plan,
                "plan_source": prep.binding(destination / "recovery-plan.json"),
                "implementation": _implementation_bindings()}
    prep._publish(target, contract)
    return {"status": "ts_only_recovery_frozen_not_run", "authorized": authorized,
            "contract": prep.binding(target), "native_jobs_launched": 0,
            "historical_charged_attempts": sum(snapshot["counts"].values()),
            "control_root": str(control.root), **execution.FLAGS}


def reopen_recovery(workdir):
    workdir = Path(workdir).resolve()
    contract = prep._json(workdir / "recovery-contract.json")
    if (set(contract) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or contract["schema_version"] != 1 or contract["kind"] != KIND + "_sealed"):
        raise ValueError("Invalid source-bound recovery contract")
    if canonical_hash(prep._json(prep._reopen(contract["plan_source"]))) != canonical_hash(contract["plan"]):
        raise ValueError("Frozen recovery plan changed")
    for item in contract["implementation"]:
        prep._reopen(item)
    if contract["implementation"] != _implementation_bindings():
        raise ValueError("Live recovery implementation inventory changed")
    _validate_plan(contract["plan"], workdir)
    return contract


def run_recovery(workdir, *, rmg_env):
    workdir = Path(workdir).resolve()
    sealed = reopen_recovery(workdir)
    if sealed["plan"]["authorized"] is not True:
        raise ValueError("TS-only recovery has not been separately authorized")
    prior, original, settings, control, request = _validate_plan(sealed["plan"], workdir)
    plan = original["plan"]
    from storca.rate_transfer_runtime import runtime_environment
    runtime = runtime_environment(plan["runtime_inventory"]) if plan["schema_version"] == 2 else nullcontext()
    fragments = prep._json(prep._reopen(prior["plan"]["retained_fragments"]))
    with runtime:
        return execution._execute_validated_context(workdir, plan=plan, settings=settings, rmg_env=rmg_env,
            control_root=control.root, context_contract_path=workdir / "recovery-contract.json",
            validate_context=lambda: reopen_recovery(workdir), path_options={
                "route_search_policy": continuation._policy(prior["plan"]["search_policy"]),
                "retained_fragment_request": {"schema_version": 1,
                    "workdir": fragments["stationary_points_root"], "sources": fragments["sources"]},
                "retained_partial_neb_request": request})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare", help="freeze a fresh TS-only recovery, zero native jobs")
    for name in ("prior_campaign", "prior_capsule", "partial_request", "diagnostic_report", "destination"):
        prepare.add_argument(name, type=Path)
    prepare.add_argument("--authorized", action="store_true", help="record explicit user consent in this NEW plan")
    for name in ("inspect", "run", "resume"):
        child = commands.add_parser(name)
        child.add_argument("workdir", type=Path)
        if name != "inspect":
            child.add_argument("--rmg-env", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        result = prepare_recovery(args.prior_campaign, args.prior_capsule, args.partial_request,
                                  args.diagnostic_report, args.destination, authorized=args.authorized)
    elif args.command == "inspect":
        contract = reopen_recovery(args.workdir)
        result = {"status": "ts_only_recovery_bindings_reopened", "authorized": contract["plan"]["authorized"],
                  "native_jobs_launched": 0, "control_contract": contract["plan"]["native_control_contract"]}
    else:
        result = run_recovery(args.workdir, rmg_env=args.rmg_env)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
