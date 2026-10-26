"""Adaptive route-search adapters over the established physical executors.

This module creates no process runner or compute ledger.  Non-preparation work
requires the existing owned NativeExecutionControl. Geometry is proposal
evidence; only the existing final TS/IRC/reached-species gates grant success.
"""
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import re
import shutil

import numpy as np

from .route_search import (RouteSearchPolicy, RouteSearchController, SearchAction,
                           ActionOutcome, RouteValidation)


def _json(value):
    return json.loads(json.dumps(value, default=str, allow_nan=False))


def _binding(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _reopen(sources):
    for source in sources:
        if _binding(source["path"]) != source:
            raise ValueError("Adaptive route source changed")


def _publish(path, value):
    """Immutable action receipts; mutable route summaries use the old writer."""
    path = Path(path)
    # Match the canonical bytes written by _atomic, including on first resume.
    text = json.dumps(_json(value), sort_keys=True, separators=(",", ":")) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text() != text:
            raise ValueError("Retained adaptive action receipt changed")
    else:
        from .route_search import _atomic
        _atomic(path, json.loads(text))


def _geometry_portfolio(route, root, mapping, fragment_xyz, environment):
    from .route_geometry import _build_endpoint_side, _fragment_coordinates_from_xyz
    from .reaction_seed_geometry import (derive_geometry_context, propose_seed_pairs,
                                         rank_diverse_candidates, seal_seed_portfolio,
                                         SeedGeometryPolicy)
    reactant, product = (_build_endpoint_side(values) for values in
                         (route.reactant_smiles, route.product_smiles))
    r = np.vstack(_fragment_coordinates_from_xyz(reactant, fragment_xyz["reactants"]))
    p = np.vstack(_fragment_coordinates_from_xyz(product, fragment_xyz["products"]))
    mapping_pairs = dict(mapping["all_atom_mapping"])
    inverse = {right: left for left, right in mapping_pairs.items()}
    declared_maps = [atom.GetAtomMapNum() for molecule in reactant.explicit_molecules
                     for atom in molecule.GetAtoms()]
    exact_maps = all(declared_maps) and len(set(declared_maps)) == len(declared_maps)
    atoms = [{"map_number": declared_maps[index] if exact_maps else index + 1, "element": element,
              "isotope": reactant.isotopes[index]} for index, element in enumerate(reactant.elements)]
    context = derive_geometry_context(
        atoms=atoms,
        reactants={"coordinates": r.tolist(), "fragment_ids": reactant.explicit_fragment,
                   "bonds": [[a, b, order] for (a, b), order in reactant.bonds.items()],
                   "fragment_states": list(zip(mapping["surface"]["reactant_fragment_charges"],
                                               mapping["surface"]["reactant_fragment_multiplicities"])),
                   "source_bindings": [_binding(x) for x in fragment_xyz["reactants"]]},
        products={"coordinates": [p[mapping_pairs[i]].tolist() for i in range(len(r))],
                  "fragment_ids": [product.explicit_fragment[mapping_pairs[i]] for i in range(len(r))],
                  "bonds": [[inverse[a], inverse[b], order] for (a, b), order in product.bonds.items()],
                  "fragment_states": list(zip(mapping["surface"]["product_fragment_charges"],
                                              mapping["surface"]["product_fragment_multiplicities"])),
                  "source_bindings": [_binding(x) for x in fragment_xyz["products"]]},
        charge=route.charge, multiplicity=route.multiplicity, environment=environment,
        parent_evidence=({"route": route.as_dict(), "mapping": mapping,
                          "map_namespace": "declared_nuclear_maps" if exact_maps else "resolved_route_index_namespace"},),
    )
    policy = SeedGeometryPolicy()
    candidates = propose_seed_pairs(context, policy)
    # Keep rejected alternatives in the immutable portfolio, not just winners.
    portfolio = seal_seed_portfolio(context, candidates, root, policy)
    ranked = rank_diverse_candidates(candidates, context, policy)
    ranks = {item.candidate_id: index for index, item in enumerate(ranked)}
    records = []
    from .route_geometry import validate_route_endpoint_connectivity
    for item in portfolio["candidates"]:
        if item["candidate_id"] in ranks:
            checks = {side: validate_route_endpoint_connectivity(
                route, Path(item[name + "_xyz"]), side=side, mapping_result=mapping)
                for side, name in (("reactants", "reactant"), ("products", "product"))}
            if all(check.get("valid") for check in checks.values()):
                records.append({**item, "rank": ranks[item["candidate_id"]],
                                "map_namespace": "declared_nuclear_maps" if exact_maps else "resolved_route_index_namespace",
                                "endpoint_connectivity": checks})
    records.sort(key=lambda item: item["rank"])
    return portfolio, records


def _complete_route(trial):
    """No preliminary TS/IRC status can substitute for reached-species gates."""
    if trial.get("status") != "verified_transition_state_path":
        return False
    ts = trial.get("transition_state") or {}
    irc = trial.get("irc") or {}
    stationary = trial.get("stationary_points") or {}
    records = [item for side in ("reactants", "products") for item in stationary.get(side, [])]
    return bool(ts.get("frequency_check", {}).get("NumImaginary") == 1
                and ts.get("reaction_coordinate_overlap", {}).get("status") == "matched"
                and (trial.get("path_evidence") or {}).get("intended_mode_match") is True
                and irc.get("completion", {}).get("valid")
                and irc.get("endpoint_validation", {}).get("valid")
                and stationary.get("valid")
                and stationary.get("geometry_origin") == "locally_refined_reached_irc_endpoints"
                and stationary.get("reactants")
                and stationary.get("products") and records
                and all(item.get("status") in {"validated_minimum", "validated_atom"}
                        for item in records))


def _complete_multistep(trial):
    from .route_search_recovery import validate_multistep_native_route
    if trial.get("status") != "verified_multistep_path":
        return False
    chain = trial.get("chain_sources") or {}
    if set(chain) != {"reactants", "intermediate", "products"}:
        return False
    _reopen(chain.values())
    segmented = trial.get("segmented_verification") or {}
    segments = segmented.get("segments") or []
    if len(segments) != 2 or [item.get("label") for item in segments] != [
            "reactants-to-intermediate", "intermediate-to-products"]:
        return False
    first, second = [item.get("native_candidate") or {} for item in segments]
    if not first.get("products") or first.get("products") != second.get("reactants"):
        return False
    return validate_multistep_native_route(trial)["valid"]


def _remaining_interrupted_cost(original, outcome, control, workspace):
    """Forecast only unfinished stages; never refund actual ledger charges.

    Credits require exact owned, normally returned native evidence inside the
    interrupted workspace. They affect planning only: the existing stage
    executors still reopen and physically validate every reused calculation.
    Failed attempts, pending work and arbitrary side files earn no credit.
    """
    from . import reaction_path as rp
    from .rate_transfer_control import classify_native_input
    evidence = outcome["evidence"]
    _reopen(evidence["owned_sources"])
    snapshot = control.snapshot()
    if snapshot["pending"]:
        raise ValueError("Interrupted native ownership is still pending")
    ids = evidence.get("reusable_reservation_ids", outcome["native_reservation_ids"])
    if (not isinstance(ids, (list, tuple)) or len(set(ids)) != len(ids)
            or any(type(value) is not int or value <= 0 for value in ids)):
        raise ValueError("Invalid interrupted native reservation inventory")
    covered = {row["path"]: row for row in evidence["owned_sources"]}
    credit = {}
    root = Path(workspace).resolve(strict=True)
    for identifier in ids:
        reservation = snapshot["reservations"].get(identifier)
        terminal = snapshot["terminals"].get(identifier)
        if reservation is None or terminal is None:
            raise ValueError("Interrupted native ownership is unresolved")
        if terminal.get("outcome") != "returned_native_work_not_physically_validated":
            continue
        inp = Path(reservation["input"]).resolve(strict=True)
        if not inp.is_relative_to(root):
            # The existing branch executor owns its original IRC directory,
            # while reached-species refinement uses this action's directory.
            # Admit only that exact hash-bound branch input, not its siblings.
            source = original.get("inputs", {}).get("source_irc")
            branch = original.get("parameters", {}).get("branch")
            if original.get("kind") != "continue_verification" or source is None or branch not in {"forward", "backward"}:
                raise ValueError("Interrupted native input belongs to a different workspace")
            _reopen([source])
            expected = Path(source["path"]).parent / "irc-branch-continuation" / branch / "irc-downhill.inp"
            if reservation["category"] != "irc" or inp != expected.resolve(strict=True):
                raise ValueError("Interrupted native input belongs to a different workspace")
        if any(covered.get(row["path"]) != row for row in terminal["sources"]):
            raise ValueError("Interrupted receipt omits owned native sources")
        control.validate_reuse(inp, reservation["input_sha256"])
        category = classify_native_input(inp.read_text())
        if category != reservation["category"]:
            raise ValueError("Interrupted native category differs from its reservation")
        output = inp.with_suffix(".out")
        if not rp._normal_orca_output(output):
            continue
        if category == "optimization_path" and not rp._orca_optimization_converged(output):
            continue
        credit[category] = credit.get(category, 0) + 1
    dispatch = dict(original["cost_envelope"])
    completion = dict(original["completion_cost"] or {})
    for category, count in credit.items():
        if count > dispatch.get(category, 0) + completion.get(category, 0):
            raise ValueError("Owned stages exceed the interrupted action forecast")
        take = min(count, dispatch.get(category, 0))
        dispatch[category] = dispatch.get(category, 0) - take
        completion[category] = completion.get(category, 0) - (count - take)
    return ({k: v for k, v in dispatch.items() if v},
            {k: v for k, v in completion.items() if v}, list(ids))


def _retained_planner_context(request, current, policy):
    """Keep original decision IDs across an explicitly archived code repair.

    The caller's NEW execution contract binds current sources. This evidence
    binds the OLD producer sources without pretending their paths still hold
    those bytes. Every non-implementation scientific/policy field must match.
    """
    if set(request) != {"contract", "implementation_sources"}:
        raise ValueError("Explicit archived planner lineage required")
    _reopen([request["contract"]])
    previous = json.loads(Path(request["contract"]["path"]).read_text())
    if (set(previous) != {"schema_version", "kind", "policy", "context"}
            or previous["schema_version"] != 1 or previous["kind"] != "storca_adaptive_route_search"
            or _json(previous["policy"]) != _json(policy.as_dict())):
        raise ValueError("Retained planner policy differs")
    old = previous["context"]
    if _json({k: v for k, v in old.items() if k != "implementation"}) != _json(
            {k: v for k, v in current.items() if k != "implementation"}):
        raise ValueError("Retained planner scientific context differs")
    expected = {row["path"]: row for row in old["implementation"]}
    found = {}
    for row in request["implementation_sources"]:
        if set(row) != {"original", "archived"}:
            raise ValueError("Malformed archived planner source")
        original, archived = row["original"], row["archived"]
        if (original["path"] in found or expected.get(original["path"]) != original
                or original["sha256"] != archived["sha256"]):
            raise ValueError("Archived planner source inventory differs")
        _reopen([archived])
        found[original["path"]] = original
    if found != expected or len(expected) != len(old["implementation"]):
        raise ValueError("Incomplete archived planner implementation")
    _reopen(current["implementation"])
    return old


def _admit_retained_partial_neb(request, control, scientific_context, root):
    """Readmit paid candidate evidence, never the failed job as a verified TS."""
    from .partial_neb_recovery import inspect_partial_neb_recovery
    required = {"schema_version", "kind", "input_path", "context_binding", "expected_context",
                "source_bindings", "selection", "scf_solver", "solver_evidence", "selection_reason"}
    if (set(request) != required or type(request["schema_version"]) is not int
            or request["schema_version"] != 1 or request["kind"] != "owned_partial_neb_ts_recovery_request"
            or request["selection"] not in {"converged_ci", "last_evaluated_opt"}
            or request["scf_solver"] not in {"original", "trah"}):
        raise ValueError("Explicit owned partial-NEB candidate request required")
    current = {key: value for key, value in scientific_context.items() if key != "path"}
    if _json(current) != _json(request["expected_context"]):
        raise ValueError("Retained candidate differs from current scientific context")
    _reopen([request["solver_evidence"]])
    comparison = json.loads(Path(request["solver_evidence"]["path"]).read_text())
    if comparison.get("both_electronic_solves_completed") is not True:
        raise ValueError("Partial recovery needs the completed controlled electronic comparison")
    assessment = inspect_partial_neb_recovery(
        request["input_path"], control=control, context_binding=request["context_binding"],
        expected_context=request["expected_context"], source_bindings=request["source_bindings"],
        selection=request["selection"])
    if assessment.get("valid") is not True or assessment.get("eligible_for_ts_recovery") is not True:
        raise ValueError("Paid partial NEB not admitted: " + str(assessment.get("failure_reason")))
    native = Path(request["input_path"]).resolve(strict=True)
    root = Path(root).resolve()
    if root == native.parent or root.is_relative_to(native.parent):
        raise ValueError("Partial recovery must not write into historical native evidence")
    text = native.read_text()
    reactant = re.findall(r"(?im)^\s*\*\s+xyzfile\s+-?\d+\s+\d+\s+(\S+)", text)
    product = re.findall(r'(?im)^\s*Product\s+"([^"\n]+)"\s*$', text)
    if len(reactant) != 1 or len(product) != 1:
        raise ValueError("Retained mapped endpoints are ambiguous")
    _publish(root / "admission.json", {"request": request, "assessment": assessment})
    geometry = assessment["candidate_geometry"]
    content = "\n".join([str(len(geometry["elements"])), "Paid partial NEB candidate; not a verified TS",
        *(f"{element} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}"
          for element, point in zip(geometry["elements"], geometry["coordinates"]))]) + "\n"
    candidate = root / "candidate.xyz"
    if candidate.exists():
        if candidate.read_text() != content:
            raise ValueError("Retained partial candidate geometry changed")
    else:
        with candidate.open("x") as stream:
            stream.write(content)
    identifier = hashlib.sha256(json.dumps({
        "native_input_sha256": _binding(native)["sha256"],
        "candidate_sha256": _binding(candidate)["sha256"],
        "selection": request["selection"], "scf_solver": request["scf_solver"]}, sort_keys=True).encode()).hexdigest()
    seed = {"candidate_id": "partial-" + identifier, "rank": -1,
            "reactant_xyz": str((native.parent / reactant[0].strip('"')).resolve(strict=True)),
            "product_xyz": str((native.parent / product[0]).resolve(strict=True)),
            "ts_seed": _binding(candidate), "admission": _binding(root / "admission.json")}
    _reopen(assessment["sources"])
    return assessment, seed


def run_adaptive_reaction_path(route, workdir, *, policy, prepare_only=False,
                              ncores=1, method_keywords=None, timeout_seconds=14400,
                              neb_images=8, condition_contract=None,
                              surface_resolution_policy=None, retained_fragment_request=None,
                              retained_partial_neb_request=None,
                              retained_planner_context=None,
                              thermochemistry_temperature_K=None):
    from . import reaction_path as rp
    from .rate_transfer_control import active_execution_control, TransferControlError
    from .route_geometry import resolve_route_atom_mapping, validate_route_endpoint_connectivity
    from .surface_resolution import SurfaceResolutionPolicy

    workdir = Path(workdir).resolve()
    policy = policy if isinstance(policy, RouteSearchPolicy) else RouteSearchPolicy.from_dict(policy)
    surface = surface_resolution_policy or SurfaceResolutionPolicy()
    control = active_execution_control()
    if not prepare_only and control is None:
        raise TransferControlError("Adaptive native work requires an owned campaign controller")
    mapping = resolve_route_atom_mapping(route)
    environment = rp._surface_environment_applicability(route, condition_contract, surface, method_keywords)
    canonical = {"status": "surface_unresolved", "route_id": route.route_id,
                 "mapping": mapping, "reaction_equation": route.reaction_equation,
                 "path_classification": "surface_unresolved",
                 "convergence": {"status": "path_unresolved"},
                 "rate_model": {"status": "rate_model_unavailable"},
                 "propagation": {"status": "not_run"}}
    result = {"schema_version": 1, "kind": "generic_orca_route_verification",
              "route": route.as_dict(), "route_verification": canonical,
              "environment_applicability": environment, "orientations": [],
              "surface_resolution_policy": surface.as_dict()}
    manifest = workdir / "route-verification.json"
    if not mapping.get("valid") or environment["status"] != "applicable":
        canonical["status"] = ("mapping_unresolved" if not mapping.get("valid") else
                               "environment_model_escalation_required")
        rp._write_json(manifest, result)
        return {**result, "path": manifest}
    selected = mapping["surface"]["selected_multiplicity"]
    if route.multiplicity != selected:
        from .route_verify import RouteSpec
        route = RouteSpec(**{**route.as_dict(), "multiplicity": selected})
    rp._initialize_surface_budget(workdir, surface)
    result["surface_resolution_scientific_context"] = rp._initialize_surface_context(
        workdir, route=route, mapping=mapping, environment_applicability=environment,
        policy=surface, ncores=ncores, method_keywords=method_keywords,
        **({"thermochemistry_temperature_K": thermochemistry_temperature_K}
           if thermochemistry_temperature_K is not None else {}))
    if retained_fragment_request is not None:
        request = retained_fragment_request
        if request.get("schema_version") != 1 or not request.get("sources"):
            raise ValueError("Independent fragment reuse needs exact source bindings")
        _reopen(request["sources"])
        source_root = Path(request["workdir"]).resolve(strict=True)
        if workdir == source_root or workdir.is_relative_to(source_root):
            raise ValueError("Adaptive search cannot write into retained fragment evidence")
        stationary = rp.reuse_stationary_fragment_records(
            route, source_root, mapping, ncores=ncores,
            method_keywords=method_keywords, retained_endpoints=None)
        _reopen(request["sources"])
    else:
        stationary = rp._stationary_fragment_records(
            route, workdir / "stationary-points", mapping, execute=not prepare_only,
            ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
            surface_resolution_policy=surface)
    canonical["stationary_points"] = stationary
    if not prepare_only and not stationary.get("valid"):
        canonical["status"] = "stationary_point_validation_failed"
        rp._write_json(manifest, result)
        return {**result, "path": manifest}
    fragment_xyz = rp._stationary_fragment_xyz(stationary, optimized=bool(stationary.get("valid")))
    portfolio, seeds = _geometry_portfolio(
        route, workdir / "route-search" / "candidates", mapping, fragment_xyz,
        {**dict(route.environment or {}), "model": environment["canonical_model"]})
    canonical["endpoint_reference"] = {
        "kind": "finite_encounter" if rp._route_uses_separated_fragment_endpoints(route) else "bound_minimum",
        "rate_claim_allowed": False, "mechanism_claim_allowed": False}
    result["endpoint_complexes"] = portfolio
    retained_seed = None
    if retained_partial_neb_request is not None:
        if control is None:
            raise TransferControlError("Partial native evidence requires the existing owned controller")
        assessment, retained_seed = _admit_retained_partial_neb(
            retained_partial_neb_request, control, result["surface_resolution_scientific_context"],
            workdir / "route-search" / "retained-partial-neb")
        result["retained_partial_neb"] = assessment
        seeds = [retained_seed]  # TS-only recovery never authorizes a new NEB sweep.
    if prepare_only:
        canonical["status"] = "prepared"
        result["search_preparation"] = {"qualified_candidates": len(seeds), "native_jobs_launched": 0,
                                        "policy": policy.as_dict()}
        rp._write_json(manifest, result)
        return {**result, "path": manifest}

    context = {"route": route.as_dict(), "mapping": mapping, "environment": environment,
               "method_keywords": method_keywords, "ncores": ncores,
               "surface_policy": surface.as_dict(), "portfolio": _binding(
                   workdir / "route-search" / "candidates" / "portfolio.json")}
    context["implementation"] = [_binding(Path(__file__).with_name(name)) for name in (
        "route_search_execution.py", "route_search.py", "route_search_recovery.py",
        "reaction_seed_geometry.py", "reaction_path.py", "surface_resolution.py", "partial_neb_recovery.py")]
    if retained_seed is not None:
        context["retained_partial_neb"] = {
            "admission": retained_seed["admission"],
            "selection": retained_partial_neb_request["selection"],
            "scf_solver": retained_partial_neb_request["scf_solver"],
            "solver_evidence": retained_partial_neb_request["solver_evidence"]}
    if retained_planner_context is not None:
        original_context = _retained_planner_context(retained_planner_context, context, policy)
        result["planner_implementation_continuation"] = {
            "previous_contract": retained_planner_context["contract"],
            "archived_sources": retained_planner_context["implementation_sources"],
            "executing_implementation": context["implementation"]}
        context = original_context
    separated = rp._route_uses_separated_fragment_endpoints(route)
    molecules = sum(item.get("kind") != "monatomic_stationary_evidence"
                    for side in ("reactants", "products") for item in stationary[side])
    atoms = sum(item.get("kind") == "monatomic_stationary_evidence"
                for side in ("reactants", "products") for item in stationary[side])
    reached_cost = {"optimization_path": molecules, "frequency": molecules, "single_point": atoms}
    completion = {"optimization_path": molecules + 1, "frequency": molecules + 1,
                  "single_point": atoms, "irc": 2}
    seed_by_id = {item["candidate_id"]: item for item in seeds}

    def budget():
        snapshot = control.snapshot()
        return {"remaining": {key: control.contract["limits"][key] - value
                              for key, value in snapshot["counts"].items()},
                "pending": snapshot["pending"], "reservation_ids": list(snapshot["reservations"]),
                "wall_seconds_remaining": control.contract["limits"]["total_wall_seconds"] - snapshot["wall_seconds"],
                "core_seconds_remaining": control.contract["limits"]["total_core_seconds"] - snapshot["core_seconds"]}

    def action_for(kind, seed, *, parent=None):
        inputs = {"reactant": _binding(seed["reactant_xyz"]), "product": _binding(seed["product_xyz"])}
        if parent is not None:
            trial = parent["outcome"]["evidence"].get("trial", {})
            if kind == "verify_ts":
                inputs["ts_seed"] = _binding(trial["ts_seed"])
            inputs["parent_receipt"] = parent["outcome"]["evidence"]["receipt"]
        return SearchAction(
            kind, inputs, context=context,
            parameters={"candidate_id": seed["candidate_id"]},
            artifacts={"seed": seed}, rank=(seed["rank"],),
            cost_envelope=({"optimization_path": 1} if kind == "neb" else
                           {"optimization_path": 1, "frequency": 1, "irc": 1} if kind == "verify_ts" else
                           {"optimization_path": 2, "frequency": 1, "irc": 1}),
            completion_cost=completion if kind == "neb" else reached_cost,
            reason=("Best distinct qualified paired geometry" if parent is None else
                    "Use retained candidate evidence rather than repeat the same search"))

    def actions(state):
        if state.completed:
            last = state.completed[-1]["outcome"]
            if last["execution"] == "failed" and last["physical"] in {
                    "native_execution_failed", "neb_execution_failed", "path_execution_failed", "operational_block"}:
                # A broken reader/runtime is not a reason to spend every new
                # geometry on the same operational failure.
                return []
        if retained_seed is None:
            values = [action_for("neb", item) for item in seeds]
        else:
            values = [SearchAction("verify_ts", {
                "reactant": _binding(retained_seed["reactant_xyz"]),
                "product": _binding(retained_seed["product_xyz"]),
                "ts_seed": retained_seed["ts_seed"], "partial_admission": retained_seed["admission"],
                "solver_evidence": retained_partial_neb_request["solver_evidence"]}, context=context,
                parameters={"candidate_id": retained_seed["candidate_id"],
                            "scf_solver": retained_partial_neb_request["scf_solver"]},
                artifacts={"seed": retained_seed}, rank=(-1,),
                cost_envelope={"optimization_path": 1, "frequency": 1, "irc": 1},
                completion_cost=reached_cost,
                reason="Refine the owned partial-stage candidate without repeating its converged NEB")]
        for previous in state.completed:
            evidence = previous["outcome"]["evidence"]
            seed = seed_by_id.get(previous["action"]["parameters"].get("candidate_id"))
            if seed is None:
                continue
            trial = evidence.get("trial") or {}
            if trial.get("status") == "irc_path_incomplete":
                irc = trial.get("irc") or {}
                complete = irc.get("completion") or {}
                done = {item.get("branch") for item in irc.get("branch_continuations", [])}
                for branch in ("forward", "backward"):
                    direction = complete.get("directions", {}).get(branch, {})
                    if (complete.get("normal_termination") and direction.get("iteration_limit_reached")
                            and branch not in done):
                        inputs = {"source_receipt": evidence["receipt"],
                                  "source_irc": _binding(irc["execution"]["input"])}
                        values.append(SearchAction(
                            "continue_verification", inputs, context=context,
                            parameters={"candidate_id": seed["candidate_id"], "branch": branch},
                            artifacts={"seed": seed, "parent_trial": trial},
                            cost_envelope={"irc": 1}, completion_cost=reached_cost,
                            rank=(seed["rank"],), reason="Extend one normally completed capped IRC branch"))
                        break
            intermediate = trial.get("intermediate_evidence") or {}
            qualified = [item for item in intermediate.get("candidates", [])
                         if item.get("status") == "validated_intermediate_minimum"
                         and (item.get("intermediate_identity") or {}).get("status") == "resolved"]
            if qualified:
                values.append(SearchAction(
                    "verify_segments", {"source_receipt": evidence["receipt"],
                                        "intermediate": _binding(qualified[0]["optimized_xyz"])},
                    context=context, parameters={"candidate_id": seed["candidate_id"]},
                    artifacts={"seed": seed, "intermediate_evidence": {
                        **intermediate, "candidates": qualified[:1]}},
                    cost_envelope={"optimization_path": 8, "frequency": 6, "irc": 2},
                    completion_cost={}, reason="Verify both segments through one qualified intermediate"))
            if evidence.get("interrupted_action"):
                original = evidence["interrupted_action"]
                dispatch, completion_remaining, reusable_ids = _remaining_interrupted_cost(
                    original, previous["outcome"], control, evidence["interrupted_workspace"])
                values.append(SearchAction(
                    "continue_verification", {"stopped_native_evidence": evidence["owned_sources"],
                                              "original_inputs": original["inputs"]},
                    context=context, parameters={"candidate_id": seed["candidate_id"],
                                                  "interrupted_fingerprint": evidence["interrupted_fingerprint"]},
                    artifacts={"seed": seed, "resume_action": original,
                               "resume_root": evidence["interrupted_workspace"],
                               "resume_reservation_ids": reusable_ids},
                    cost_envelope=dispatch, completion_cost=completion_remaining,
                    reason="Resume stopped paid stages without repeating completed native work"))
            effective_kind = (previous["action"].get("artifacts", {}).get("resume_action") or
                              previous["action"])["kind"]
            if effective_kind != "neb":
                continue
            if trial.get("intermediate_indices"):
                values.append(SearchAction(
                    "validate_intermediate", {"source_receipt": evidence["receipt"],
                                              "native_output": _binding(trial["neb"]["output"])},
                    context=context, parameters={"candidate_id": seed["candidate_id"],
                                                  "indices": trial["intermediate_indices"]},
                    artifacts={"seed": seed, "source_path": str(Path(trial["neb"]["output"]).parent)},
                    cost_envelope={"optimization_path": 2, "frequency": 2},
                    completion_cost={"optimization_path": 8, "frequency": 6, "irc": 2},
                    reason="Qualify an observed intermediate before choosing segment verification"))
            if trial.get("ts_seed"):
                values.append(action_for("verify_ts", seed, parent=previous))
            elif trial.get("status") in {"endpoint_reference_unqualified", "neb_path_unresolved"}:
                values.append(action_for("guided_scan", seed, parent=previous))
        return values

    from .route_geometry import _build_endpoint_side, _mapped_bond_changes
    changes = _mapped_bond_changes(_build_endpoint_side(route.reactant_smiles),
                                  _build_endpoint_side(route.product_smiles),
                                  dict(mapping["all_atom_mapping"]))
    intended_coordinates = [
        {"atom_indices": pair, "start_angstrom": 1.0, "end_angstrom": 2.0}
        for pair in changes["broken_bonds"]] + [
        {"atom_indices": pair, "start_angstrom": 2.0, "end_angstrom": 1.0}
        for pair in changes["formed_bonds"]]
    intended_coordinates += [
        {"atom_indices": item["atom_indices"], "start_angstrom": 1.0,
         "end_angstrom": 0.9 if item["product_order"] > item["reactant_order"] else 1.1}
        for item in changes["bond_order_changes"]]

    def execute(action):
        _reopen([value for value in action.inputs.values() if isinstance(value, dict) and "path" in value])
        identifier = action.fingerprint(policy.validation_policy, context)
        root = workdir / "route-search" / "actions" / identifier
        receipt = root / "outcome.json"
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            _reopen(saved["sources"])
            old = saved["outcome"]
            if old["execution"] not in {"blocked_resources", "not_started"} or old["native_reservation_ids"]:
                old["evidence"]["receipt"] = _binding(receipt)
                return old
            # A resource block before reservation is not a scientific attempt.
            # Keep its receipt and append a new execution receipt after retry.
            sequence = 1
            while (root / f"outcome-retry-{sequence:04d}.json").exists():
                retry = root / f"outcome-retry-{sequence:04d}.json"
                previous = json.loads(retry.read_text())
                _reopen(previous["sources"])
                if previous["outcome"]["execution"] not in {"blocked_resources", "not_started"}:
                    value = previous["outcome"]
                    value["evidence"]["receipt"] = _binding(retry)
                    return value
                sequence += 1
            receipt = root / f"outcome-retry-{sequence:04d}.json"
        root.mkdir(parents=True, exist_ok=True)
        operation = action
        execution_root = root
        if action.artifacts.get("resume_action"):
            operation = SearchAction.from_dict(action.artifacts["resume_action"])
            execution_root = Path(action.artifacts["resume_root"])
        seed = action.artifacts["seed"]
        before = set(control.snapshot()["reservations"])
        trial = {}
        try:
            if operation.kind == "neb":
                endpoints = []
                for name in ("reactant", "product"):
                    source = Path(seed[name + "_xyz"])
                    target = execution_root / (name + ".xyz")
                    if target.exists() and target.read_bytes() != source.read_bytes():
                        raise ValueError("Attempted endpoint geometry cannot be overwritten")
                    if not target.exists():
                        shutil.copy2(source, target)
                    endpoints.append(target)
                inp = rp.create_orca_neb_ts_input(
                    *endpoints, charge=route.charge, multiplicity=route.multiplicity,
                    ncores=ncores, nimages=neb_images, method_keywords=method_keywords,
                    preopt_ends=False, prepare_fragments=False)
                native = rp._run_or_resume(inp, timeout_seconds=timeout_seconds)
                stem = Path(native.get("input") or inp).stem
                ts_seed = rp._locate_neb_ts(execution_root, stem) if native["status"] == "completed" else None
                warning = native.get("path_quality_warning") or {}
                trial = {"status": "neb_execution_failed" if native["status"] == "failed" else
                         "ts_seed_available" if ts_seed else
                         "endpoint_reference_unqualified" if separated and warning.get("kind") ==
                         "submerged_climbing_image_below_endpoint_asymptote" else "neb_path_unresolved",
                         "neb": native, "ts_seed": str(ts_seed) if ts_seed else None,
                         "endpoint_reference": canonical["endpoint_reference"]}
                trial["intermediate_indices"] = rp._neb_intermediate_image_indices(Path(native["output"]))
                initial = execution_root / f"{stem}_initial_path_trj.xyz"
                if initial.is_file():
                    frames = rp.xyz_frames(initial)
                    labels = rp.read_xyz(endpoints[0])[0]
                    intact = bool(len(frames) >= 2 and frames[0][0] == labels and frames[-1][0] == labels)
                    changes_native = []
                    if intact:
                        for frame, endpoint in ((frames[0][1], endpoints[0]), (frames[-1][1], endpoints[1])):
                            original = rp.read_xyz(endpoint)[1]
                            distances = lambda x: np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
                            changes_native.append(float(np.max(np.abs(distances(frame) - distances(original)))))
                    trial["native_endpoint_preparation"] = {
                        "source": _binding(initial), "frame_identity_intact": intact,
                        "maximum_distance_changes_angstrom": changes_native,
                        "nominal_reference_preserved": intact and max(changes_native) <= 1e-4,
                        "energetic_asymptote_claim_allowed": False}
            elif operation.kind == "validate_intermediate":
                intermediate = rp._validate_neb_intermediates(
                    route, neb_output=Path(operation.inputs["native_output"]["path"]),
                    path_dir=Path(operation.artifacts["source_path"]),
                    ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
                    candidate_indices=operation.parameters["indices"],
                    surface_resolution_policy=surface)
                trial = {"status": intermediate["status"], "intermediate_evidence": intermediate}
            elif operation.kind == "verify_segments":
                intermediate = operation.artifacts["intermediate_evidence"]
                segmented = rp._run_segmented_intermediate_nebs(
                    route, reactant_xyz=Path(seed["reactant_xyz"]), product_xyz=Path(seed["product_xyz"]),
                    intermediate_evidence=intermediate, path_dir=execution_root,
                    separated_fragments=separated, ncores=ncores, method_keywords=method_keywords,
                    timeout_seconds=timeout_seconds, neb_images=neb_images,
                    surface_resolution_policy=replace(surface, irc_continuation_policy_id="adaptive-planner-owned-v1"),
                    parent_endpoint_stationary=stationary, require_reached_minima=True)
                trial = {"status": "verified_multistep_path", "segmented_verification": segmented,
                         "chain_sources": {"reactants": _binding(seed["reactant_xyz"]),
                                           "intermediate": _binding(intermediate["candidates"][0]["optimized_xyz"]),
                                           "products": _binding(seed["product_xyz"])}}
                if not _complete_multistep(trial):
                    trial["status"] = "segment_path_unresolved"
            elif operation.kind == "continue_verification" and operation.parameters.get("branch"):
                trial = _json(operation.artifacts["parent_trial"])
                irc = trial["irc"]
                source = Path(operation.inputs["source_irc"]["path"])
                assembly = (rp.assemble_completed_irc_branches(source)
                            if action.artifacts.get("resume_action") else {})
                if assembly.get("valid"):
                    # The raw assembly reader has checked this exact old
                    # execution. Do not rewrite its historical result merely
                    # to change a transient "resumed" flag.
                    extension = json.loads((source.parent / "irc-branch-continuation" /
                        operation.parameters["branch"] / "continuation-result.json").read_text())
                else:
                    extension = rp.run_irc_branch_continuation(
                        source, branch=operation.parameters["branch"], timeout_seconds=timeout_seconds)
                    assembly = rp.assemble_completed_irc_branches(source)
                irc.setdefault("branch_continuations", []).append(extension)
                irc["branch_assembly"] = assembly
                if assembly.get("valid"):
                    trajectory = Path(assembly["trajectory"])
                    identity = rp._validate_irc_route_endpoints(
                        route, trajectory, Path(seed["reactant_xyz"]), Path(seed["product_xyz"]),
                        mapping_result=mapping, separated_fragments=separated,
                        rmsd_tolerance_angstrom=surface.endpoint_rmsd_tolerance_angstrom)
                    irc.update(completion={"valid": True, "assembly": assembly},
                               endpoint_validation=identity, trajectory=str(trajectory))
                    trial.update(status="verified_transition_state_path" if identity.get("valid") else
                                 "irc_endpoint_validation_failed", trajectory=str(trajectory))
                    if identity.get("valid"):
                        ts = trial.get("transition_state") or {}
                        trial.update(path_classification="verified_barriered_path",
                            endpoint_reference_geometries={"reactants": seed["reactant_xyz"],
                                                           "products": seed["product_xyz"]},
                            path_evidence={"kind": "orca_ts_frequency_assembled_irc",
                                "one_meaningful_imaginary_frequency": ts.get("frequency_check", {}).get("NumImaginary") == 1,
                                "intended_mode_match": ts.get("reaction_coordinate_overlap", {}).get("status") == "matched",
                                "irc_matches_declared_endpoints": True, "trajectory": str(trajectory)})
            else:
                kwargs = dict(reactant_xyz=Path(seed["reactant_xyz"]), product_xyz=Path(seed["product_xyz"]),
                              mapping_result=mapping, separated_fragments=separated,
                              ncores=ncores, method_keywords=method_keywords,
                              timeout_seconds=timeout_seconds, surface_resolution_policy=replace(
                                  surface, irc_continuation_policy_id="adaptive-planner-owned-v1"))
                if operation.kind == "verify_ts":
                    if operation.parameters.get("scf_solver") == "trah":
                        kwargs["method_keywords"] = [*(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]), "TRAH"]
                    trial = rp._verify_ts_seed(route, Path(operation.inputs["ts_seed"]["path"]), execution_root,
                                               intended_coordinates=intended_coordinates, **kwargs)
                    if retained_seed is not None:
                        trial["partial_neb_recovery"] = {
                            "admission": retained_seed["admission"],
                            "solver_evidence": retained_partial_neb_request["solver_evidence"],
                            "solver": operation.parameters["scf_solver"],
                            "scope": "TS refinement, independent frequency and IRC only; unchanged physical validation"}
                else:
                    trial = rp._run_bond_change_guided_ts_recovery(route, path_dir=execution_root,
                                                                  bond_changes=changes, **kwargs)
            if trial.get("status") == "verified_transition_state_path":
                    refined = rp._refine_reached_irc_stationary_points(
                        route, trial, execution_root / "reached-stationary-points", mapping,
                        ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
                        reference_endpoints={"reactants": Path(seed["reactant_xyz"]),
                                             "products": Path(seed["product_xyz"])},
                        surface_resolution_policy=surface)
                    trial.update(refined)
                    if not _complete_route(trial):
                        trial["status"] = "reached_minimum_validation_failed"
            status = trial.get("status", "unresolved")
            operational_failure = status in {"neb_execution_failed", "path_execution_failed"}
            outcome = ActionOutcome(execution="failed" if operational_failure else "completed",
                                    reader="not_run" if operational_failure else "passed",
                                    physical=status, evidence={"trial": _json(trial)},
                                    reason=trial.get("failure_reason", ""))
        except (rp.SurfaceResolutionBudgetExhausted, TransferControlError) as error:
            message = str(error)
            if isinstance(error, TransferControlError) and any(word in message.lower() for word in (
                    "changed", "mutat", "integrity", "differs", "unowned", "no owned")):
                raise
            status = "blocked_resources" if any(word in message.lower() for word in ("memory", "disk")) else "ownership_unresolved" if any(word in message.lower() for word in ("pending", "owner", "unclosed")) else "failed"
            physical = ("budget_exhausted" if isinstance(error, rp.SurfaceResolutionBudgetExhausted)
                        or any(word in message.lower() for word in ("cap reached", "remaining native timeout", "exceed"))
                        else "operational_block")
            outcome = ActionOutcome(execution=status, reader="not_run", physical=physical,
                                    evidence={"trial": _json(trial)}, reason=message)
        except (RuntimeError, OSError) as error:
            outcome = ActionOutcome(execution="failed", reader="not_run", physical="native_execution_failed",
                                    evidence={"trial": _json(trial)}, reason=str(error))
        after = control.snapshot()
        spent = sorted(set(after["reservations"]) - before)
        sources = {}
        for attempt in spent:
            for source in after["terminals"].get(attempt, {}).get("sources", []):
                sources[source["path"]] = source
        preparation = trial.get("native_endpoint_preparation") or {}
        if preparation.get("source"):
            sources[preparation["source"]["path"]] = preparation["source"]
        value = outcome.as_dict()
        value["native_reservation_ids"] = spent
        if (spent or action.artifacts.get("resume_action")) and outcome.execution in {"blocked_resources", "ownership_unresolved"}:
            reusable_ids = sorted(set(action.artifacts.get("resume_reservation_ids", [])) | set(spent))
            for attempt in reusable_ids:
                for source in after["terminals"].get(attempt, {}).get("sources", []):
                    sources[source["path"]] = source
            value["evidence"].update(
                interrupted_action=operation.as_dict(), interrupted_fingerprint=identifier,
                interrupted_workspace=str(execution_root), owned_sources=list(sources.values()),
                reusable_reservation_ids=reusable_ids)
        # Receipt is not self-referential; the journal binds it after publication.
        _publish(receipt, {"outcome": value, "sources": list(sources.values())})
        value["evidence"]["receipt"] = _binding(receipt)
        return value

    def validate(state):
        for item in state.completed:
            evidence = item["outcome"]["evidence"]
            if evidence.get("receipt"):
                _reopen([evidence["receipt"]])
                saved = json.loads(Path(evidence["receipt"]["path"]).read_text())
                _reopen(saved["sources"])
            trial = evidence.get("trial") or {}
            if _complete_route(trial):
                return RouteValidation(status="verified_route", reached_minima_validated=True,
                                       all_required_segments_verified=True,
                                       evidence={"trial": trial, "receipt": evidence["receipt"]})
            if _complete_multistep(trial):
                return RouteValidation(status="verified_multistep_route", reached_minima_validated=True,
                                       all_required_segments_verified=True,
                                       evidence={"trial": trial, "receipt": evidence["receipt"]})
        return RouteValidation()

    def recover(decision):
        action = SearchAction.from_dict(decision["action"])
        root = workdir / "route-search" / "actions" / decision["fingerprint"]
        receipt = root / "outcome.json"
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            _reopen(saved["sources"])
            outcome = saved["outcome"]
            outcome["evidence"]["receipt"] = _binding(receipt)
            return outcome
        previous = set(decision["budget_snapshot"].get("reservation_ids", []))
        snapshot = control.snapshot()
        if not (set(snapshot["reservations"]) - previous):
            return ActionOutcome(execution="not_started", reader="not_run", physical="not_attempted")
        spent = sorted(set(snapshot["reservations"]) - previous)
        if snapshot["pending"]:
            return None
        sources = {source["path"]: source for attempt in spent
                   for source in snapshot["terminals"][attempt]["sources"]}
        _reopen(sources.values())
        return ActionOutcome(physical="stopped_paid_stage_pending_completion",
                             native_reservation_ids=tuple(spent), evidence={
                                 "interrupted_action": action.as_dict(),
                                 "interrupted_fingerprint": decision["fingerprint"],
                                 "interrupted_workspace": str(root),
                                 "owned_sources": list(sources.values())})

    controller = RouteSearchController(workdir / "route-search" / "planner", policy, context)
    search = controller.run(actions, execute, budget, validate, recover_pending=recover)
    result["route_search"] = search
    for item in search["completed"]:
        result["orientations"].append(item["outcome"]["evidence"].get("trial", {}))
    final = search["validation"].get("evidence", {}).get("trial")
    if final and _complete_multistep(final):
        canonical.update(status="verified_multistep_path", path_classification="multistep_intermediate_detected",
                         multistep_path=final["segmented_verification"],
                         convergence={"status": "segmented_paths_verified_rate_pending"},
                         chain_sources=final["chain_sources"])
    elif final and _complete_route(final):
        canonical.update(final, convergence={"status": "path_verified_rate_pending"},
                         path_classification="verified_barriered_path",
                         arkane_inputs={"reactants": final["stationary_points"]["reactants"],
                                        "products": final["stationary_points"]["products"],
                                        "transition_state": final["transition_state"]})
    else:
        status = search["status"]
        if (search["completed"] and search["completed"][-1]["outcome"]["execution"] == "failed"
                and search["completed"][-1]["outcome"]["physical"] != "budget_exhausted"):
            status = "native_or_reader_failure_requires_recovery"
        canonical.update(status=status, convergence={"status": "path_unresolved",
                         "reason": status})
    rp._write_json(manifest, result)
    return {**result, "path": manifest}
