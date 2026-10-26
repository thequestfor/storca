"""Bounded, explicit historical-band hypothesis in a NEW P1 workspace.

The default is a read-only admission/preparation replay, not chemistry. With
--execute, the actual generic path, Arkane and repaired-model callbacks run
under the ORIGINAL shared P1 ledger unless --generous-budget explicitly records
a separately authorized extended run. No reference answer, barrier or uncertainty
is supplied. A central model calculation is not scientific P1 acceptance.
"""

import argparse
from contextlib import ExitStack
import json
import os
from pathlib import Path
import signal
import sys
import time
import traceback
from unittest.mock import patch

from run_p1_preflight import REPOSITORY, THREAD_LIMITS, digest, manifest_conditions, save
from run_p1_workflow import artifact, production_budget


SOURCE_WORKFLOW = REPOSITORY / "runs/p1-cyclobutene-workflow-002/workflow"
SOURCE_PATH_RESULT = SOURCE_WORKFLOW / "verification-engine-003/iteration-001-rmg-02be421ca3cb/path-result.json"
POLICY_ID = "generic-orca-path-v10-reached-endpoint-minima+historical-partial-neb-ts-seed-v1"
REUSE_POLICY_ID = POLICY_ID + "+completed-ts-frequency-reuse-v1"
COMPLETED_IRC_POLICY_ID = POLICY_ID + "+completed-irc-reuse-v1"
EXTENDED_AUTHORIZATION = (
    "User explicitly requested another visible Kitty run with enough time, "
    "after clarifying that compute is not constrained (2026-10-07)."
)
EXTENDED_BUDGET = {
    "overall_timeout_seconds": 48.0 * 3600.0,
    "cleanup_grace_seconds": 20.0,
    "maximum_cores": 4,
    "route_timeout_seconds_per_job": 4.0 * 3600.0,
    "verification_neb_images": 8,
    "verification_max_iterations": 12,
    "rmg_walltime": "00:00:30:00",
    "rmg_max_processes": 1,
    "rmg_max_iterations": 300,
    "rmg_max_edge_species": 1000,
    "rmg_max_objects_per_iteration": 1,
    "rmg_retry_attempts": 2,
    "rmg_hard_timeout_seconds": None,
    "rmg_heartbeat_timeout_seconds": None,
}


class ForbiddenExecution(BaseException):
    pass


def json_paths(value):
    """Serialize native Path fields without hiding other invalid evidence types."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: json_paths(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_paths(item) for item in value]
    return value


def inventory(root):
    """Retain bytes AND timestamps; inspection must not modify source evidence."""
    root = Path(root).resolve(strict=True)
    files = {}
    directories = {str(root): root.stat().st_mtime_ns}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError("Historical evidence must not contain symbolic links")
        if path.is_file():
            files[str(path)] = {**artifact(path), "mtime_ns": path.stat().st_mtime_ns}
        elif path.is_dir():
            directories[str(path)] = path.stat().st_mtime_ns
    return files, directories


def historical_request(source_workflow, source_path_result, timeout_seconds):
    from storca.reaction_path import inspect_neb_time_limit_checkpoint

    source_workflow = Path(source_workflow).resolve(strict=True)
    source_path_result = Path(source_path_result).resolve(strict=True)
    if not source_path_result.is_relative_to(source_workflow):
        raise ValueError("Historical result must belong to the retained workflow")
    result = json.loads(source_path_result.read_text())
    eligible = [item for item in result.get("orientations", [])
                if (item.get("neb") or {}).get("status") == "failed"
                and (item.get("path_evidence") or {}).get("kind") == "incomplete_neb_path"]
    if len(eligible) != 1:
        raise ValueError("This bounded trial requires one unambiguous historical incomplete band")
    source_input = Path(eligible[0]["neb"]["input"])
    checkpoint = inspect_neb_time_limit_checkpoint(source_input, timeout_seconds=timeout_seconds)
    if checkpoint is None:
        raise ValueError("Historical native time-limit checkpoint is inadmissible")
    files, directories = inventory(source_workflow)
    request = {
        "schema_version": 1, "kind": "historical_partial_neb_ts_seed_request",
        "source_path_result": artifact(source_path_result),
        "source_execution_contract": artifact(source_workflow / "verification-execution-contract.json"),
        "source_checkpoint": checkpoint,
        "source_artifacts": [{"path": key, "sha256": value["sha256"]}
                             for key, value in files.items()],
    }
    return request, result, (files, directories)


def run(output, *, rmg_env, execute=False, source_workflow=SOURCE_WORKFLOW,
        source_path_result=SOURCE_PATH_RESULT,
        manifest=REPOSITORY / "benchmarks/stability/p1_cyclobutene_development.json",
        resume_run=None, generous_budget=False):
    started = time.monotonic()
    policy_id = REUSE_POLICY_ID if resume_run is not None else POLICY_ID
    from storca.reaction_path import inspect_historical_neb_ts_seed, run_generic_reaction_path
    from storca.route_verify import RouteSpec, route_spec_from_rmg_route

    output = Path(output).resolve()
    source_workflow = Path(source_workflow).resolve(strict=True)
    if output.exists() or output.is_relative_to(source_workflow):
        raise ValueError("A new workspace outside the original workflow is required")
    if execute:
        prefix = "p1-cyclobutene-extended-" if generous_budget else "p1-cyclobutene-workflow-"
        if output.parent != REPOSITORY / "runs" or not output.name.startswith(prefix):
            raise ValueError("Executed hypotheses must join the explicitly selected budget ledger; "
                             f"required directory prefix is {prefix}; "
                             "ordinary trials require the existing P1 production ledger")
    raw = Path(manifest).read_bytes()
    smiles, scenario, conditions = manifest_conditions(json.loads(raw))
    contract = json.loads((source_workflow / "verification-execution-contract.json").read_text())
    if (contract["condition_contract"] != conditions.as_dict()
            or contract["parent"] != {"smiles": smiles, "charge": 0, "multiplicity": 1}
            or contract["rmg"]["environment"] != rmg_env
            or contract["resources"]["ncores"] != 4):
        raise ValueError("Historical trial differs from the frozen P1 scientific conditions")
    if generous_budget:
        budget = dict(EXTENDED_BUDGET)
        ledger = {"kind": "p1_user_authorized_extended_budget", "authorization": EXTENDED_AUTHORIZATION,
                  "ceiling_seconds": budget["overall_timeout_seconds"], "consumed_seconds": 0.0,
                  "remaining_seconds": budget["overall_timeout_seconds"],
                  "legacy_ledger_renewed": False, "original_workflow": artifact(
                      source_workflow / "verification-execution-contract.json")}
    else:
        budget, ledger = production_budget(REPOSITORY / "runs", digest(raw))
    source_timeout = contract["resources"]["route_job_timeout_seconds"]
    request, historical, before = historical_request(
        source_workflow, source_path_result, source_timeout)
    resource_budget = None
    if generous_budget:
        resource_budget = {
            "schema_version": 1, "kind": "historical_trial_resource_budget",
            "source_execution_contract": request["source_execution_contract"],
            "previous_timeout_seconds_per_job": source_timeout,
            "timeout_seconds_per_job": budget["route_timeout_seconds_per_job"],
            "overall_timeout_seconds": budget["overall_timeout_seconds"],
            "maximum_cores": budget["maximum_cores"], "authorization": EXTENDED_AUTHORIZATION,
        }
    resumed_before = None
    if resume_run is not None:
        resume_run = Path(resume_run).resolve(strict=True)
        if output.is_relative_to(resume_run):
            raise ValueError("Derived continuation must be outside the retained run")
        prior_contract = resume_run / "execution-contract.json"
        prior_result = resume_run / "result.json"
        candidates = [path for path in resume_run.glob(
            "verification-engine/**/historical-path/route-verification.json")
            if json.loads(path.read_text()).get("route") == historical["route"]]
        if len(candidates) != 1:
            raise ValueError("Continuation requires one unambiguous retained historical TS trial")
        resumed_before = inventory(resume_run)
        container_request = {
            "schema_version": 1, "kind": "completed_ts_frequency_reuse_request",
            "source_path_result": artifact(candidates[0]),
            "source_execution_contract": artifact(prior_contract),
            "source_run_result": artifact(prior_result),
            "source_artifacts": [{"path": key, "sha256": value["sha256"]}
                                 for key, value in resumed_before[0].items()],
        }
        prior_path = json.loads(candidates[0].read_text())
        if ((prior_path.get("historical_neb_seed_trial") or {}).get("irc") or {}).get("completion", {}).get("valid") is True:
            # Completed IRC owner may itself have reused TS stages owned by an
            # older run. Keep their original ownership/hashes; the new gate
            # independently reopens the completed trace's separate container.
            prior_request_path = Path(json.loads(prior_contract.read_text())["source_request"]["path"])
            prior_request = json.loads(prior_request_path.read_text())
            request["retained_ts_frequency"] = prior_request.get("retained_ts_frequency", container_request)
            request["retained_completed_irc"] = {
                "schema_version": 1, "kind": "completed_irc_reuse_request",
                "ts_frequency_request": container_request,
            }
            policy_id = COMPLETED_IRC_POLICY_ID
        else:
            request["retained_ts_frequency"] = container_request
    route_value = dict(historical["route"])
    for field in ("reactant_smiles", "product_smiles", "reactant_labels", "product_labels",
                  "reactant_charges", "product_charges", "reactant_multiplicities",
                  "product_multiplicities", "surface_multiplicities"):
        if field in route_value:
            route_value[field] = tuple(route_value[field])
    for field in ("broken_bonds", "formed_bonds", "atom_mapping"):
        if field in route_value:
            route_value[field] = tuple(tuple(pair) for pair in route_value[field])
    spec = RouteSpec(**route_value)
    keywords = contract["method_profile"]["orca_keywords"]
    admission = inspect_historical_neb_ts_seed(
        spec, request, ncores=4, method_keywords=keywords, condition_contract=conditions.as_dict())
    output.mkdir(parents=True, exist_ok=False)
    save(output, "historical-seed-request.json", request)
    save(output, "execution-contract.json", {
        "kind": "p1_bounded_historical_neb_seed_hypothesis", "manifest_sha256": digest(raw),
        "driver_sha256": digest(Path(__file__).read_bytes()), "budget": budget,
        "shared_budget_accounting": ledger, "condition_contract": conditions.as_dict(),
        "method_profile": contract["method_profile"], "rmg_environment": rmg_env,
        "source_request": artifact(output / "historical-seed-request.json"),
        "path_verifier_policy_id": policy_id, "execute": execute,
        "interpretation": ("New reached-endpoint validation from a revalidated completed IRC; no reused endpoint/rate claim."
                           if "retained_completed_irc" in request else
                           "New TS-search hypothesis; historical band is neither a saddle nor a barrier."),
        "rate_uncertainty_evidence_supplied": False, "reference_answers_supplied": False,
        "resume_source": str(resume_run) if resume_run is not None else None,
        "historical_trial_resource_budget": resource_budget,
    })
    save(output, "admission.json", admission)
    outcome = {"status": "failed", "execute": execute, "p1_scientific_gate_achieved": False,
               "estimated_time_to_retention_seconds": None, "native_dispatches": [],
               "budget": budget, "path_verifier_policy_id": policy_id}
    previous_handlers, previous_timer = {}, (0.0, 0.0)
    timer_started = None
    interruption = {"received_signal": None}

    def deadline(_signal, _frame):
        raise KeyboardInterrupt("P1 historical trial reached its shared overall deadline")

    def controlled_sigint(signum, _frame):
        if interruption["received_signal"] is not None:
            return
        interruption["received_signal"] = signal.Signals(signum).name
        # Permit existing owned-process cleanup to finish without a second
        # signal or the campaign alarm interrupting that same cleanup.
        signal.setitimer(signal.ITIMER_REAL, 0)
        raise KeyboardInterrupt("P1 historical trial received controlled SIGINT")

    def verify_route(route, workdir):
        current = route_spec_from_rmg_route(smiles, route)
        for field in ("route_id", "reaction_equation", "parent_smiles", "charge", "multiplicity",
                      "reactant_smiles", "product_smiles", "reactant_charges", "product_charges",
                      "reactant_multiplicities", "product_multiplicities"):
            left, right = getattr(current, field), getattr(spec, field)
            if isinstance(left, (list, tuple)):
                left, right = tuple(left), tuple(right)
            if left != right:
                raise ValueError("Native selected route differs from historical hypothesis: " + field)
        result = run_generic_reaction_path(spec, workdir / "historical-path", ncores=4, method_keywords=keywords,
            timeout_seconds=budget["route_timeout_seconds_per_job"], orientations=1,
            neb_images=budget["verification_neb_images"], condition_contract=conditions.as_dict(),
            historical_neb_seed_request=request, historical_trial_resource_budget=resource_budget)
        return result

    try:
        if not admission.get("valid"):
            outcome.update(status="source_admission_refused", reason=admission.get("failure_reason") or admission.get("reason"))
        elif not execute:
            attempts = []
            def forbid(*args, **kwargs):
                attempts.append("external_execution")
                raise ForbiddenExecution("Read-only preparation attempted an external job")
            with ExitStack() as stack:
                for name in ("subprocess.run", "subprocess.Popen", "storca.reaction_path._run_or_resume",
                             "storca.reaction_path.run_optimization_and_frequency", "src.orca_runner.run_orca"):
                    stack.enter_context(patch(name, side_effect=forbid))
                prepared = run_generic_reaction_path(spec, output / "prepared-path", prepare_only=True,
                    ncores=4, method_keywords=keywords, orientations=1, neb_images=8,
                    timeout_seconds=budget["route_timeout_seconds_per_job"],
                    condition_contract=conditions.as_dict(), historical_neb_seed_request=request,
                    historical_trial_resource_budget=resource_budget)
            save(output, "prepared-path-result.json", json_paths(prepared))
            path_status = prepared["route_verification"]["status"]
            outcome.update(status=("prepared" if path_status == "historical_ts_seed_prepared"
                                   else "preparation_refused"), path_status=path_status,
                           forbidden_execution_attempts=attempts, chemistry_jobs_executed=0)
            if path_status != "historical_ts_seed_prepared":
                outcome["failure_reason"] = "Native preparation refused: " + path_status
        else:
            os.environ.update(THREAD_LIMITS)
            os.environ["STORCA_ORCA_LIVE_OUTPUT"] = "1"
            print("STORCA cyclobutene: reached IRC endpoints -> validated minima -> Arkane -> conditional network handoff", flush=True)
            print(f"Conditions: {conditions.temperature_K} K, {conditions.pressure_bar} bar, dark gas", flush=True)
            print(f"Budget: {budget['route_timeout_seconds_per_job'] / 3600:g} hours per ORCA job; "
                  f"{budget['overall_timeout_seconds'] / 3600:g} hours overall; four aggregate cores", flush=True)
            print(f"Output: {output}", flush=True)
            if resume_run is not None:
                if "retained_completed_irc" in request:
                    print(f"Revalidating completed TS/frequency/IRC from {resume_run}; no NEB, TS search or IRC rerun", flush=True)
                    print("New endpoint Opt/Freq calculations; rates/network claims retain all scientific gates", flush=True)
                else:
                    print(f"Reusing completed TS/frequency from {resume_run}; the incomplete IRC is NOT reused", flush=True)
            print("Scientific checks unchanged. A completed job is not automatically a defensible lifetime.\n", flush=True)
            timer_started = time.monotonic()
            remaining = budget["overall_timeout_seconds"] - (timer_started - started)
            if remaining <= 0:
                raise TimeoutError("P1 budget exhausted during source admission; no chemistry dispatched")
            previous_handlers[signal.SIGINT] = signal.signal(signal.SIGINT, controlled_sigint)
            previous_handlers[signal.SIGALRM] = signal.signal(signal.SIGALRM, deadline)
            previous_timer = signal.setitimer(signal.ITIMER_REAL, remaining)
            from storca.stability import collect_rmg_evidence
            from storca.verification_adapters import make_arkane_rate_builder
            from storca.verification_engine import VerificationEngineConfig, run_verification_engine
            initial = json.loads((source_workflow / "stability.json").read_text())["rmg_evidence"]
            rmg = contract["rmg"]
            resources = dict(rmg["resource_config"])
            if generous_budget:
                resources.update({name: budget["rmg_" + name] for name in (
                    "walltime", "max_processes", "max_iterations", "max_edge_species", "max_objects_per_iteration")})
            def rerun_model(library_paths, workdir):
                print("\nSTORCA: propagating the repaired RMG model", flush=True)
                libraries = list(library_paths)
                for item in rmg["base_reaction_libraries"]:
                    path = Path(item["path"])
                    if path not in libraries:
                        libraries.append(path)
                repaired = collect_rmg_evidence(smiles, workdir, rmg_env=rmg_env,
                    rmg_command=rmg["command"], barrier_threshold=rmg["barrier_threshold_kcal_mol"],
                    temperature=conditions.temperature_K, pressure=conditions.pressure_bar,
                    rmg_walltime=resources["walltime"], rmg_max_processes=resources["max_processes"],
                    rmg_max_iterations=resources["max_iterations"], rmg_max_edge_species=resources["max_edge_species"],
                    rmg_max_objects_per_iteration=resources["max_objects_per_iteration"],
                    rmg_maximum_heavy_atoms=rmg["maximum_heavy_atoms"], scenario=scenario, conditions=conditions,
                    reaction_libraries=libraries, database_reaction_libraries=rmg["database_reaction_libraries"],
                    rmg_retry_attempts=budget["rmg_retry_attempts"] if generous_budget else rmg["retry_attempts"],
                    rmg_hard_timeout_seconds=budget["rmg_hard_timeout_seconds"] if generous_budget else rmg["hard_timeout_seconds"],
                    rmg_heartbeat_timeout_seconds=budget["rmg_heartbeat_timeout_seconds"] if generous_budget else rmg["heartbeat_timeout_seconds"],
                    rmg_start_method=rmg["start_method"])
                return {"rmg_evidence": repaired, "condition_contract": conditions.as_dict()}
            native_rate_builder = make_arkane_rate_builder(conditions.as_dict(),
                method_profile=contract["method_profile"], rmg_env=rmg_env)
            def build_rate(*args, **kwargs):
                print("\nSTORCA: building rates from the verified path using Arkane", flush=True)
                return native_rate_builder(*args, **kwargs)
            summary = run_verification_engine(initial, conditions.as_dict(), output,
                verify_route=verify_route, build_rate=build_rate if generous_budget else native_rate_builder,
                rerun_model=rerun_model,
                config=VerificationEngineConfig(max_iterations=budget.get("verification_max_iterations", 3)), rmg_env=rmg_env,
                path_verifier_policy_id=policy_id)
            save(output, "verification-result.json", summary)
            outcome.update(status="workflow_returned", verification_status=summary.get("status"),
                           estimated_time_to_retention_seconds=summary.get("estimated_time_to_retention_seconds"))
    except BaseException as error:
        outcome.update(failure_reason=f"{type(error).__name__}: {error}", traceback=traceback.format_exc())
    finally:
        if previous_handlers:
            signal.setitimer(signal.ITIMER_REAL, 0)
            for sig, handler in previous_handlers.items():
                signal.signal(sig, handler)
            if previous_timer[0]:
                signal.setitimer(signal.ITIMER_REAL, max(.001, previous_timer[0] - (time.monotonic() - timer_started)), previous_timer[1])
        after = inventory(source_workflow)
        outcome.update(elapsed_seconds=time.monotonic() - started,
            original_sources_unchanged=before == after, original_file_count=len(before[0]),
            original_directory_count=len(before[1]), manifest_unchanged=Path(manifest).read_bytes() == raw)
        if resumed_before is not None:
            outcome.update(retained_stage_sources_unchanged=resumed_before == inventory(resume_run),
                           retained_stage_file_count=len(resumed_before[0]),
                           retained_stage_directory_count=len(resumed_before[1]))
        if not outcome["original_sources_unchanged"] or not outcome["manifest_unchanged"]:
            outcome.update(status="source_changed", estimated_time_to_retention_seconds=None,
                           failure_reason="Original evidence or preregistration changed during the trial")
        if outcome.get("retained_stage_sources_unchanged") is False:
            outcome.update(status="source_changed", estimated_time_to_retention_seconds=None,
                           failure_reason="Retained TS stage sources changed during continuation")
        outcome["interruption_control"] = interruption
        outcome["native_dispatches"] = [artifact(path) for path in sorted(output.rglob("*.process.json"))]
        save(output, "result.json", outcome)
    print(json.dumps(outcome, sort_keys=True))
    return outcome


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rmg-env", required=True)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--resume-run", type=Path,
                        help="Reuse admitted completed TS/frequency and, when present, completed IRC in a new workspace")
    parser.add_argument("--generous-budget", action="store_true",
                        help="Record a separate authorized run: four hours per ORCA job, 48 hours overall, four cores")
    args = parser.parse_args()
    result = run(args.output_dir, rmg_env=args.rmg_env, execute=args.execute, resume_run=args.resume_run,
                 generous_budget=args.generous_budget)
    raise SystemExit(0 if result["status"] in {"prepared", "workflow_returned"} else 1)
