"""Opt-in, bounded known-reaction transfer; never an autonomous fate claim.

The preparation/scoring package is deliberately not a prediction backend. This
adapter consumes only its answer-free reaction, freezes a separate execution
contract, and reuses the existing path, owned native runner and composite code.
Preparing a template does not authorize a calculation. No experimental rates
or OH/methane numerical counting factors enter this module.
"""
from __future__ import annotations

from contextlib import nullcontext
from pathlib import Path

from src.orca_job_settings import NativeJobSettings
from . import rate_transfer_preparation as prep
from .scientific_panel import canonical_hash
from .route_verify import RouteSpec
from .composite_rate_recipe import ATOMIC_ENERGY_POLICY


KIND = "storca_known_reaction_execution"
LABELS = ("r1", "r2", "p1", "p2", "TS")
LOW_KEYWORDS = ["UKS", "WB97X-D3", "def2-TZVPP", "NoAutoStart", "ExtremeSCF",
                "DefGrid3", "RIJONX", "def2/J"]
# Separate from, and explicitly superseding only a future execution proposal:
# the original 6/4 preparation is preserved verbatim. Minimum generic flow is
# 8 Opt/path + 7 molecular frequencies + 1 bidirectional IRC + 7 SP (two low
# level atomic H points and five matched high-level points). Recovery is capped.
PROPOSED_CONTROL = {
    "schema_version": 1, "kind": "storca_rate_transfer_execution_control",
    "campaign_id": "REPLACE_WITH_UNIQUE_CAMPAIGN_ID", "authorized": False,
    "resources": {"ncores": 4, "maxcore_mb_per_rank": 1000,
                  "reserve_disk_bytes": 4 * 1024**3, "reserve_memory_bytes": 1024**3},
    "limits": {"optimization_path": 12, "frequency": 10, "irc": 2, "single_point": 9,
               "per_job_seconds": 21600, "total_wall_seconds": 172800,
               "total_core_seconds": 691200},
    "extra_environment": {},
}
FLAGS = {"physical_rate_accuracy_claim_allowed": False, "rate_claim_allowed": False,
         "lifetime_claim_allowed": False, "calibrated_uncertainty_claim_allowed": False,
         "network_reinjection_allowed": False, "production_defaults_changed": False}


def _publish(path, value):
    prep._publish(Path(path), value)


def _prediction(value):
    """Exact allowlist, not arbitrary nested metadata carrying reference answers."""
    if not isinstance(value, dict) or value.get("case_id") not in prep.CASE_IDS:
        raise ValueError("A supported answer-free known-reaction request is required")
    expected = prep.prediction_request(value["case_id"], value["method_recipe"]["candidate"])
    if canonical_hash(expected) != canonical_hash(value):
        raise ValueError("Prediction inputs differ from the frozen known-reaction protocol")
    return value


def prediction_to_route(prediction):
    value = _prediction(prediction)
    r, p = value["reactants"], value["products"]
    return RouteSpec(
        route_id="known-reaction:" + value["case_id"], source="known_reaction",
        parent_smiles=".".join(s["mapped_smiles"] for s in r),
        reactant_smiles=tuple(s["mapped_smiles"] for s in r),
        product_smiles=tuple(s["mapped_smiles"] for s in p),
        reactant_labels=LABELS[:2], product_labels=LABELS[2:4],
        charge=value["total_charge"], multiplicity=value["surface_multiplicity"],
        reactant_charges=tuple(s["charge"] for s in r),
        product_charges=tuple(s["charge"] for s in p),
        reactant_multiplicities=tuple(s["multiplicity"] for s in r),
        product_multiplicities=tuple(s["multiplicity"] for s in p),
        surface_multiplicities=(value["surface_multiplicity"],),
        protocol="generic_endpoint_path", kinetic_relevance="known_reaction_not_network_discovery",
        evidence_level="user_declared_known_reaction", environment=dict(value["environment"]),
    )


def prepare_execution(preparation_dir, case_id, output_dir):
    """Write an unapproved draft, zero native work; never alter the preparation."""
    preparation_dir, output_dir = Path(preparation_dir), Path(output_dir)
    manifest = prep.validate_preparation(preparation_dir)
    predictions = prep._json(preparation_dir / "prediction-requests.json")
    try:
        prediction = next(p for p in predictions["cases"] if p["case_id"] == case_id)
    except StopIteration as error:
        raise ValueError("Unsupported transfer case") from error
    _prediction(prediction)
    output_dir.mkdir(parents=True, exist_ok=True)
    _publish(output_dir / "prediction.json", prediction)
    draft = {
        "schema_version": 1, "kind": KIND,
        "prediction": prep.binding(output_dir / "prediction.json"),
        "native_settings": NativeJobSettings().as_dict(),
        "control": PROPOSED_CONTROL,
        "orca_binary": None,
        "counting": {label: {"external_symmetry": None, "optical_isomers": None}
                     for label in LABELS},
        "assumptions": {"electronic_state": None, "channel_counting": None, "soft_modes": None},
        "model_scope": "isolated_canonical_gas_conditional_diagnostic_not_experiment_equivalence",
        "tunneling": "Eckart", "frequency_scale_factor": 1.0,
        "electronic_method": "CCSD(T)/def2-QZVPP",
        "vibrational_method": "WB97X-D3/def2-TZVPP",
        "atomic_energy_policy": ATOMIC_ENERGY_POLICY,
        "original_budget_proposal_sha256": canonical_hash(manifest["proposed_budget"]),
    }
    _publish(output_dir / "execution-plan.json", draft)
    return {"status": "draft_not_authorized", "native_jobs_launched": 0, **FLAGS,
            "plan": prep.binding(output_dir / "execution-plan.json"),
            "required_decisions": ["separate_compute_authorization", "explicit_ORCA_runtime",
                                   "case_specific_symmetry_optical_and_channel_counts",
                                   "explicit_electronic_state_and_harmonic_assumptions"],
            "minimum_unretried_flow": {"optimization_path": 8, "frequency": 7,
                                       "irc": 1, "single_point": 7},
            "comparison_ready": False}


def _validate_plan(plan):
    required = {"schema_version", "kind", "prediction", "native_settings", "control", "orca_binary",
                "counting", "assumptions", "model_scope", "tunneling", "frequency_scale_factor",
                "electronic_method", "vibrational_method", "atomic_energy_policy",
                "original_budget_proposal_sha256"}
    version = plan.get("schema_version")
    if version == 2:
        required.add("runtime_inventory")
    if (set(plan) != required or type(version) is not int
            or version not in {1, 2} or plan["kind"] != KIND):
        raise ValueError("Unsupported execution contract or non-allowlisted fields")
    prediction = _prediction(prep._json(prep._reopen(plan["prediction"])))
    settings = NativeJobSettings.from_dict(plan["native_settings"])
    if (settings != NativeJobSettings()
            or canonical_hash(plan["native_settings"]) != canonical_hash(NativeJobSettings().as_dict())):
        raise ValueError("Transfer requires the frozen native procedure, not a method sweep")
    fixed = {"model_scope": "isolated_canonical_gas_conditional_diagnostic_not_experiment_equivalence",
             "tunneling": "Eckart", "frequency_scale_factor": 1.0,
             "electronic_method": "CCSD(T)/def2-QZVPP", "vibrational_method": "WB97X-D3/def2-TZVPP",
             "atomic_energy_policy": ATOMIC_ENERGY_POLICY}
    if canonical_hash({k: plan[k] for k in fixed}) != canonical_hash(fixed):
        raise ValueError("Frozen transfer method/model differs")
    if plan["original_budget_proposal_sha256"] != canonical_hash(prep.PROPOSED_BUDGET):
        raise ValueError("Original preparation budget identity differs")
    if set(plan["counting"]) != set(LABELS):
        raise ValueError("Explicit counts required for both asymptotes and TS")
    for label, counts in plan["counting"].items():
        if set(counts) != {"external_symmetry", "optical_isomers"} or any(
                type(v) is not int or v < 1 for v in counts.values()):
            raise ValueError("Positive integer case-specific symmetry/optical counts required: " + label)
    if set(plan["assumptions"]) != {"electronic_state", "channel_counting", "soft_modes"} or any(
            not isinstance(s, str) or not s.strip() for s in plan["assumptions"].values()):
        raise ValueError("Explicit case-specific physical assumptions required")
    binary = prep._reopen(plan["orca_binary"])
    if not binary.is_file():
        raise ValueError("Bound ORCA executable missing")
    control = plan["control"]
    # Validate the whole control schema without creating its workspace.
    from .rate_transfer_control import validate_control_contract
    validate_control_contract(control)
    if control["resources"]["ncores"] != 4 or control["resources"]["maxcore_mb_per_rank"] != 1000:
        raise ValueError("Transfer retains four cores and 1000 MB/rank")
    if (control["resources"]["reserve_disk_bytes"] < 4 * 1024**3
            or control["resources"]["reserve_memory_bytes"] < 1024**3
            or control["limits"]["per_job_seconds"] != 21600
            or control["campaign_id"] == PROPOSED_CONTROL["campaign_id"]):
        raise ValueError("Explicit campaign identity and unchanged per-job/disk safeguards required")
    minimum = {"optimization_path": 8, "frequency": 7, "irc": 1, "single_point": 7}
    if any(control["limits"][key] < n for key, n in minimum.items()):
        raise ValueError("Execution proposal cannot cover the existing mandatory path checks")
    if any(control["limits"][key] > n for key, n in PROPOSED_CONTROL["limits"].items()):
        raise ValueError("Expanded compute campaigns require a separate revised protocol")
    if prediction["case_id"] == "h_ch4":
        from .h_ch4_counting import expected_counts
        if plan["counting"] != expected_counts():
            raise ValueError("The frozen H + CH4 model requires its source-supported symmetry counts")
    if version == 2:
        from .rate_transfer_runtime import reopen_runtime
        reopen_runtime(plan["runtime_inventory"], orca_binary=plan["orca_binary"],
                       extra_environment=control["extra_environment"])
    return prediction, settings


def _implementation_bindings():
    root = Path(__file__).resolve().parents[1]
    return [prep.binding(root / name) for name in (
        "storca/rate_transfer_execution.py", "storca/rate_transfer_control.py",
        "src/orca_job_settings.py", "src/inputgen.py", "src/orca_runner.py",
        "storca/workflow.py", "storca/route_geometry.py", "storca/reaction_path.py",
        "storca/route_thermochemistry.py", "storca/composite_rate_recipe.py",
        "storca/arkane_launcher.py", "storca/arkane_rotor.py", "storca/native_rate_audit.py",
        "storca/route_verify.py", "storca/surface_resolution.py", "storca/endpoint_identity.py",
        "storca/ir_modes.py", "src/stability/freq_check.py", "storca/h_ch4_counting.py",
        "storca/rate_transfer_runtime.py", "storca/configured_mpi_runtime.py")]


def freeze_execution(plan_path, workdir):
    """Freeze a separately authored plan. Authorization is never inferred here."""
    plan_path, workdir = Path(plan_path), Path(workdir)
    plan = prep._json(plan_path)
    _validate_plan(plan)
    if (workdir.exists() and not (workdir / "execution-contract.json").exists()
            and (not workdir.is_dir() or any(workdir.iterdir()))):
        raise ValueError("Freeze into a fresh campaign, not an existing calculation directory")
    workdir.mkdir(parents=True, exist_ok=True)
    contract = {"schema_version": 1, "kind": KIND + "_sealed", "plan": plan,
                "plan_source": prep.binding(plan_path), "implementation": _implementation_bindings()}
    _publish(workdir / "execution-contract.json", contract)
    return {"status": "frozen_not_run", "authorized": plan["control"]["authorized"],
            "contract": prep.binding(workdir / "execution-contract.json"),
            "native_jobs_launched": 0, **FLAGS}


def reopen_execution(workdir):
    contract = prep._json(Path(workdir) / "execution-contract.json")
    if (set(contract) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or contract["schema_version"] != 1 or contract["kind"] != KIND + "_sealed"):
        raise ValueError("Invalid sealed execution contract")
    prep._reopen(contract["plan_source"])
    if prep._json(Path(contract["plan_source"]["path"])) != contract["plan"]:
        raise ValueError("Frozen plan changed")
    for source in contract["implementation"]:
        prep._reopen(source)
    if contract["implementation"] != _implementation_bindings():
        raise ValueError("Execution implementation inventory changed")
    _validate_plan(contract["plan"])
    return contract


def _verified_points(path_result, route):
    """Path code, not this adapter, establishes chemistry and reached minima."""
    path = Path(path_result["path"])
    recorded = prep._json(path)
    canonical = recorded["route_verification"]
    if canonical_hash(recorded["route"]) != canonical_hash(route.as_dict()):
        raise ValueError("Path belongs to another mapped reaction")
    stationary = canonical.get("stationary_points") or {}
    ts = canonical.get("transition_state") or {}
    irc = canonical.get("irc") or {}
    evidence = canonical.get("path_evidence") or {}
    if (canonical.get("status") != "verified_transition_state_path"
            or canonical.get("path_classification") != "verified_barriered_path"
            or stationary.get("valid") is not True
            or stationary.get("geometry_origin") != "locally_refined_reached_irc_endpoints"
            or ts.get("frequency_check", {}).get("NumImaginary") != 1
            or ts.get("reaction_coordinate_overlap", {}).get("status") != "matched"
            or irc.get("completion", {}).get("valid") is not True
            or irc.get("endpoint_validation", {}).get("valid") is not True
            or evidence.get("intended_mode_match") is not True):
        raise ValueError("Mapped mode, both IRC directions and actual reached minima must pass")
    points = {}
    for side, labels, smiles, charges, spins in (
        ("reactants", LABELS[:2], route.reactant_smiles, route.reactant_charges, route.reactant_multiplicities),
        ("products", LABELS[2:4], route.product_smiles, route.product_charges, route.product_multiplicities),
    ):
        rows = stationary.get(side, [])
        if len(rows) != 2:
            raise ValueError("Exactly two independently validated species required per side")
        for record, label, identity, charge, spin in zip(rows, labels, smiles, charges, spins):
            if (record.get("smiles") != identity or record.get("charge") != charge
                    or record.get("multiplicity") != spin
                    or record.get("status") not in {"validated_minimum", "validated_atom"}):
                raise ValueError("Reached stationary species identity/state differs")
            points[label] = {**record, "kind": "species", "smiles": identity}
    points["TS"] = {**ts, "kind": "transition_state", "smiles": None,
                    "charge": route.charge, "multiplicity": route.multiplicity}
    # Seal all raw scientific sources used in the rate handoff, not just a
    # mutable 'passed' label in the route summary.
    # The path summary is intentionally NOT a source identity: rerunning the
    # reader changes 'resumed' runtime metadata even when all paid native work
    # is identical. Freeze a stable handoff of raw bindings below instead.
    sources = []
    for point in points.values():
        sources.append(prep.binding(Path(point["optimized_xyz"])))
        if point.get("status") == "validated_atom":
            sources += [prep.binding(Path(point[k])) for k in ("energy_input", "energy_output")]
        else:
            output = Path(point["frequency_output"])
            sources += [prep.binding(output), prep.binding(output.with_suffix(".inp")),
                        prep.binding(output.with_suffix(".hess"))]
            optimization = point.get("optimization") or {}
            optimization_output = point.get("optimization_output") or optimization.get("output")
            if not optimization_output:
                raise ValueError("Native molecular/TS optimization evidence is missing")
            optimized_output = Path(optimization_output)
            sources += [prep.binding(optimized_output),
                        prep.binding(Path(optimization.get("input") or optimized_output.with_suffix(".inp")))]
    def execution_sources(execution):
        if not execution or not all(execution.get(key) for key in ("input", "output")):
            raise ValueError("Native IRC execution evidence is missing")
        sources.extend(prep.binding(Path(execution[key])) for key in ("input", "output"))
    execution_sources(irc.get("execution"))
    for continuation in irc.get("branch_continuations", []):
        execution_sources(continuation.get("execution"))
    sources.append(prep.binding(Path(irc["trajectory"])))
    return points, sorted({s["path"]: s for s in sources}.values(), key=lambda s: s["path"])


def _retained_verified_path(request, route, control):
    """Reopen an immutable verified handoff, never rerun or adopt native work."""
    if not isinstance(request, dict) or set(request) != {"path_summary", "handoff"}:
        raise ValueError("Exact retained path and raw handoff bindings required")
    path = prep._reopen(request["path_summary"])
    handoff_path = prep._reopen(request["handoff"])
    handoff = prep._json(handoff_path)
    _, sources = _verified_points({"path": path}, route)
    counting_path = handoff_path.parent / "counting-admission.json"
    if canonical_hash(route.as_dict()) != canonical_hash(handoff.get("route")) or (
            handoff.get("kind") != "known_reaction_verified_path_native_handoff"
            or handoff.get("schema_version") != 1
            or handoff.get("mapped_reaction_mode_and_bidirectional_IRC_passed") is not True
            or handoff.get("actual_reached_stationary_species_validated") is not True
            or handoff.get("physical_accuracy_claim_allowed") is not False
            or handoff.get("sources") != [*sources, prep.binding(counting_path)]
            or prep._json(counting_path).get("status") != "passed"):
        raise ValueError("Retained verified handoff or counting evidence differs")
    from .reaction_path import _orca_job_input_evidence
    owned = {row["input"] for row in control.snapshot()["reservations"].values()}
    for source in handoff["sources"]:
        raw = prep._reopen(source)
        if raw.suffix == ".inp":
            if str(raw) not in owned:
                raise ValueError("Retained native input has no owned receipt")
            control.validate_reuse(raw, _orca_job_input_evidence(raw)[0].hexdigest())
    return {"path": path}


def _single_point_input(point, folder, plan, *, ncores=4, maxcore_mb_per_rank=1000):
    from .route_geometry import read_xyz
    folder.mkdir(parents=True, exist_ok=True)
    elements, coordinates, _ = read_xyz(Path(point["optimized_xyz"]))
    atom = len(elements) == 1
    if type(ncores) is not int or ncores not in {1, 4}:
        raise ValueError("Explicit matched-energy rank count must be one or four")
    if ncores == 1 and (list(elements) != ["H", "H"] or point["charge"] != 0 or point["multiplicity"] != 1):
        raise ValueError("One-rank exception is restricted to neutral singlet H2")
    if (type(maxcore_mb_per_rank) is not int or maxcore_mb_per_rank not in {1000, 1500}
            or (maxcore_mb_per_rank != 1000 and (point.get("kind") != "transition_state" or ncores != 4))):
        raise ValueError("Raised MaxCore is restricted to the four-core TS energy")
    method = "HF" if atom else "CCSD(T)"
    if atom and (elements != ["H"] and tuple(elements) != ("H",)):
        raise ValueError("Only the explicit one-electron hydrogen atomic energy policy is supported")
    reference = "UHF" if point["multiplicity"] > 1 else "RHF"
    text = f"! {reference} {method} def2-QZVPP NoAutoStart VeryTightSCF FrozenCore\n"
    text += f"%pal\n  nprocs {ncores}\nend\n%maxcore {maxcore_mb_per_rank}\n%scf\n  MaxIter 500\nend\n"
    if not atom:
        text += "%mdci\n  DoSDiagnostic true\n  UseQROs false\nend\n"
    text += f"* xyz {point['charge']} {point['multiplicity']}\n"
    text += "".join(f"{e} {x:.12f} {y:.12f} {z:.12f}\n" for e, (x, y, z) in zip(elements, coordinates))
    text += "*\n"
    # High-level SPs have their own explicit HF/CC SCF controls, not the DFT
    # numerical-frequency profile. They use the same controlled dispatcher.
    target = folder / "energy.inp"
    if target.exists():
        if target.read_text() != text:
            raise ValueError("Retained matched single-point input changed")
    else:
        with target.open("x") as stream:
            stream.write(text)
    return target, ("HF/def2-QZVPP" if atom else plan["electronic_method"])


def _preflight_energy_parallelism(input_path):
    """Known ORCA MDCI restriction, not a general parallel-scaling predictor."""
    import re
    from .rate_transfer_control import _is_h2_ccsd_input
    text = Path(input_path).read_text()
    if _is_h2_ccsd_input(text) and re.findall(r"(?im)^\s*nprocs\s+(\d+)\s*$", text) != ["1"]:
        raise ValueError("H2 canonical CCSD(T) has one electron pair; this ORCA backend requires one process")


def _composite(points, energy, plan, sources):
    from .composite_rate_recipe import (CompositeRateRecipe, bind_stationary_point,
                                       bind_atomic_species, bind_artifact)
    prediction = prep._json(prep._reopen(plan["prediction"]))
    if prediction["case_id"] == "h_ch4":
        from .h_ch4_counting import assess_h_ch4_counting
        if assess_h_ch4_counting(points, plan["counting"])["status"] != "passed":
            raise ValueError("Actual geometry does not support the frozen H + CH4 counting model")
    declarations = []
    for label in LABELS:
        point, result = points[label], energy[label]
        prep._reopen(result["input"])
        prep._reopen(result["output"])
        counts = plan["counting"][label]
        if point.get("status") == "validated_atom":
            if counts != {"external_symmetry": 1, "optical_isomers": 1}:
                raise ValueError("Atomic H has no external rotor or optical multiplicity")
            declarations.append(bind_atomic_species(
                label, Path(result["output"]["path"]), energy_hartree=result["energy_hartree"],
                atomic_mass_amu=point["atomic_mass_amu"], smiles=point["smiles"],
                isotope=point.get("isotope", 0)))
        else:
            declarations.append(bind_stationary_point(
                label, point["kind"], Path(point["frequency_output"]), Path(result["output"]["path"]),
                energy_hartree=result["energy_hartree"], charge=point["charge"],
                multiplicity=point["multiplicity"], smiles=point["smiles"], **counts))
    return CompositeRateRecipe(
        label="transfer-" + prediction["case_id"], points=tuple(declarations),
        reactant_labels=LABELS[:2], product_labels=LABELS[2:4], transition_state_label="TS",
        electronic_method=plan["electronic_method"], vibrational_method=plan["vibrational_method"],
        frequency_scale_factor=plan["frequency_scale_factor"],
        temperatures_K=tuple(prediction["temperature_grid_K"]), tunneling=plan["tunneling"],
        environment="gas", surface_model="single_adiabatic_surface",
        electronic_state_assumption=plan["assumptions"]["electronic_state"],
        channel_counting_assumption=plan["assumptions"]["channel_counting"],
        soft_mode_assumption=plan["assumptions"]["soft_modes"],
        route_evidence=tuple(bind_artifact(prep._reopen(s)) for s in sources))


def run_execution(workdir, *, rmg_env):
    """Run/resume only a separately frozen authorized campaign; no rate reinjection."""
    contract = reopen_execution(workdir)
    plan = contract["plan"]
    if plan["control"]["authorized"] is not True:
        raise ValueError("This execution proposal has not been separately authorized")
    if plan["schema_version"] == 2:
        from .rate_transfer_runtime import runtime_environment
        selection = runtime_environment(plan["runtime_inventory"])
    else:
        selection = nullcontext()
    with selection:
        return _run_execution(workdir, rmg_env=rmg_env)


def _run_execution(workdir, *, rmg_env):
    workdir = Path(workdir).resolve()
    contract = reopen_execution(workdir)
    plan = contract["plan"]
    _, settings = _validate_plan(plan)
    if plan["control"]["authorized"] is not True:
        raise ValueError("This execution proposal has not been separately authorized")
    return _execute_validated_context(
        workdir, plan=plan, settings=settings, rmg_env=rmg_env,
        control_root=workdir / "native-control",
        context_contract_path=workdir / "execution-contract.json",
        validate_context=lambda: reopen_execution(workdir),
    )


def _execute_validated_context(workdir, *, plan, settings, rmg_env, control_root,
                               context_contract_path, validate_context, path_options=None,
                               execution_amendment=None, retained_verified_path=None,
                               energy_directories=None, energy_ncores=None,
                               energy_maxcore_mb_per_rank=None):
    """Shared stages behind a strict original or continuation lineage validator.

    This private body grants no authorization or source-hash exemption. Its
    caller validates the complete context, including separate continuation
    consent when applicable, and supplies the original native control ledger.
    """
    from src.orca_runner import find_orca
    from .reaction_path import run_generic_reaction_path
    from .rate_transfer_control import NativeExecutionControl, controlled_run_or_resume
    from .native_rate_audit import read_electronic_single_point, run_native_rate_audit
    from .composite_rate_recipe import write_isolated_arkane_recipe, validate_composite_rate_recipe

    workdir = Path(workdir).resolve()
    validate_context()
    prediction = _prediction(prep._json(prep._reopen(plan["prediction"])))
    if plan["control"]["authorized"] is not True:
        raise ValueError("This execution proposal has not been separately authorized")
    if prep.binding(Path(find_orca())) != plan["orca_binary"]:
        raise ValueError("Selected ORCA differs from the frozen runtime; set STORCA_ORCA_BIN explicitly")
    if not isinstance(rmg_env, str) or not rmg_env.strip():
        raise ValueError("Select the working Arkane environment with --rmg-env")
    _publish(workdir / "rmg-environment.json", {"rmg_env": rmg_env})
    route = prediction_to_route(prediction)
    if energy_directories is not None and set(energy_directories) != set(LABELS):
        raise ValueError("Exact matched-energy label destinations required")
    if energy_ncores is not None and (set(energy_ncores) != set(LABELS) or execution_amendment is None
            or any(type(n) is not int or n not in {1, 4} for n in energy_ncores.values())):
        raise ValueError("Exact authorized matched-energy core counts required")
    if energy_maxcore_mb_per_rank is not None and (execution_amendment is None
            or energy_maxcore_mb_per_rank != {label: 1500 if label == "TS" else 1000 for label in LABELS}
            or any(type(n) is not int for n in energy_maxcore_mb_per_rank.values())):
        raise ValueError("Exact authorized TS-only matched-energy memory settings required")
    if retained_verified_path is not None and (path_options or execution_amendment is None):
        raise ValueError("Energy-only reuse requires an explicit amendment and no path search options")
    control = NativeExecutionControl(control_root, plan["control"],
                                    **({"execution_amendment": execution_amendment} if execution_amendment else {}))
    with control.activate(settings):
        path_kwargs = dict(ncores=4, timeout_seconds=21600,
            method_keywords=list(LOW_KEYWORDS), orientations=1,
            condition_contract={"phase_model": "homogeneous_gas",
                                "temperature_K": prediction["temperature_grid_K"][0]},
            thermochemistry_temperature_K=prediction["temperature_grid_K"][0])
        path_kwargs.update(path_options or {})
        path_result = (_retained_verified_path(retained_verified_path, route, control)
                       if retained_verified_path is not None else
                       run_generic_reaction_path(route, workdir / "path", **path_kwargs))
        try:
            points, sources = _verified_points(path_result, route)
        except ValueError as error:
            result = {"status": "path_unresolved", "failure_reason": str(error),
                      "path": prep.binding(Path(path_result["path"])), **FLAGS}
            if (path_options or {}).get("route_search_policy") is not None:
                # Adaptive progress can legitimately produce another unresolved
                # assessment on resume. Preserve the assessed summary itself,
                # not a binding to the mutable canonical path manifest, and
                # retain every closeout rather than overwriting the first one.
                summary = prep._json(Path(path_result["path"]))
                receipt = workdir / "path-closeouts" / canonical_hash({
                    "summary": summary, "failure_reason": str(error)})
                receipt.mkdir(parents=True, exist_ok=True)
                _publish(receipt / "path-summary.json", summary)
                result["path"] = prep.binding(receipt / "path-summary.json")
                _publish(receipt / "closeout.json", result)
            else:
                _publish(workdir / "path-closeout.json", result)
            return result
        if prediction["case_id"] == "h_ch4":
            from .h_ch4_counting import assess_h_ch4_counting
            counting = assess_h_ch4_counting(points, plan["counting"])
            _publish(workdir / "counting-admission.json", counting)
            sources.append(prep.binding(workdir / "counting-admission.json"))
            if counting["status"] != "passed":
                result = {"status": "counting_model_unresolved", "counting": sources[-1], **FLAGS}
                _publish(workdir / "counting-closeout.json", result)
                return result
        _publish(workdir / "rate-source-contract.json", {
            "schema_version": 1, "kind": "known_reaction_verified_path_native_handoff",
            "route": route.as_dict(), "sources": sources,
            "mapped_reaction_mode_and_bidirectional_IRC_passed": True,
            "actual_reached_stationary_species_validated": True,
            "physical_accuracy_claim_allowed": False})
        sources.append(prep.binding(workdir / "rate-source-contract.json"))
        # Validate every planned input before spending even the first energy.
        # Historic four-rank H2 inputs remain reopenable, but cannot be launched.
        energy_inputs = {}
        for label, point in points.items():
            folder = Path(energy_directories[label]) if energy_directories is not None else workdir / "energies" / label
            ncores = (energy_ncores or {}).get(label, 4)
            memory = (energy_maxcore_mb_per_rank or {}).get(label, 1000)
            options = {**({"ncores": ncores} if ncores != 4 else {}),
                       **({"maxcore_mb_per_rank": memory} if memory != 1000 else {})}
            inp, method = _single_point_input(point, folder, plan, **options)
            _preflight_energy_parallelism(inp)
            energy_inputs[label] = (inp, method)
        energies = {}
        for label, point in points.items():
            for source in sources:
                prep._reopen(source)
            inp, method = energy_inputs[label]
            executed = controlled_run_or_resume(inp, timeout_seconds=21600,
                                               explicit_single_point_settings=True)
            executed_input = Path(executed.get("input") or inp)
            output = Path(executed.get("output") or executed_input.with_suffix(".out"))
            energies[label] = read_electronic_single_point(
                executed_input, output, method=method,
                charge=point["charge"], multiplicity=point["multiplicity"])
        _publish(workdir / "energies.json", energies)
        for source in sources:
            prep._reopen(source)
        recipe = _composite(points, energies, plan, sources)
        validate_composite_rate_recipe(recipe)
        staged = workdir / "composite"
        if not staged.exists():
            write_isolated_arkane_recipe(recipe, staged)
        else:
            # A staged recipe is reusable only as exactly the same recipe.
            if canonical_hash(prep._json(staged / "composite-rate-recipe.json")) != canonical_hash(validate_composite_rate_recipe(recipe)):
                raise ValueError("Retained composite recipe differs")
        audit = run_native_rate_audit(staged, rmg_env=rmg_env)
        for source in sources:
            prep._reopen(source)
        validate_context()
        result = {"status": "conditional_rate_available", "case_id": prediction["case_id"],
                  "rate_points": audit["rate_points"], "model_scope": plan["model_scope"],
                  "audit": prep.binding(staged / "native-rate-audit.json"),
                  "execution_contract": prep.binding(context_contract_path),
                  "comparison_ready": False, "scientific_pass": False, **FLAGS}
        _publish(workdir / "result.json", result)
        return result


def prepare_h_ch4_campaign(preparation_dir, runtime_source, output_dir):
    """Resolve only the H/CH4 conditional model, freeze it UNAUTHORIZED.

    No chemical calculation or experimental-answer ingestion takes place.
    Numerical geometry symmetry remains a required future evidence gate, not
    an assumed outcome of the native search.
    """
    from .h_ch4_counting import expected_counts, POLICY_ID
    from .rate_transfer_runtime import prepare_runtime
    output_dir = Path(output_dir)
    prepare_execution(preparation_dir, "h_ch4", output_dir / "original-draft")
    runtime = prepare_runtime(Path(runtime_source), output_dir / "runtime-selection.json")
    selected = prep._json(prep._reopen(runtime))
    plan = prep._json(output_dir / "original-draft" / "execution-plan.json")
    plan.update(schema_version=2, runtime_inventory=runtime, orca_binary=selected["orca_binary"],
                counting=expected_counts())
    identity = "h-ch4-conditional-transfer-20261017-" + canonical_hash({"workdir": str(output_dir.resolve())})[:10]
    plan["control"].update(campaign_id=identity, authorized=False,
                           extra_environment=selected["selected_environment"])
    plan["assumptions"] = {
        "electronic_state": "Ground-state unpolarized single adiabatic doublet surface. H, CH3 and TS doublets; CH4 and H2 singlets. Spin-only electronic partition; no extra electronic or nuclear-spin factor. SCF spin/root continuity and single-reference applicability remain independent diagnostics, not certified by multiplicities.",
        "channel_counting": "One mapped path, ordinary protium nuclei. H1, CH4(Td)12, H2(Dinfh)2, CH3(D3h)6, TS(C3v)3; optical1 throughout. Rotational partitions already account for four equivalent abstraction sites: no additional times-four multiplier, native reaction degeneracy1. Geometry/achirality must pass " + POLICY_ID + "; no symmetrization or posthoc count change. Source: Fernandez-Ramos et al., DOI10.1007/s00214-007-0328-0, Eq6/Table2.",
        "soft_modes": "Conventional RRHO/TST with classical rigid rotation, high-temperature spin-averaged nuclear statistics at 1000/1200/1500K, and all positive physical vibrations harmonic. No nuclear-spin-resolved ortho/para kinetics, angular correction, mode floor or extra TS basin. Exactly one physical H-transfer negative mode and both mapped IRC branches required. Eckart is conditional on positive forward and reverse E0 barriers; no claimed RRHO, tunneling or physical accuracy calibration.",
    }
    _publish(output_dir / "execution-plan.json", plan)
    frozen = freeze_execution(output_dir / "execution-plan.json", output_dir / "campaign")
    result = {"schema_version": 1, "status": "conditional_model_frozen_not_authorized",
              "contract": frozen["contract"], "plan": prep.binding(output_dir / "execution-plan.json"),
              "counting_policy": POLICY_ID, "counting_geometry_validated": False,
              "runtime_selection_bound": True, "actual_worker_libraries_proved": False,
              "native_jobs_launched": 0, "comparison_ready": False,
              "next_decision": "Separate authorization of H + CH4 only under the frozen proposed budget; primary-reference audit is independently required before scoring.",
              **FLAGS}
    _publish(output_dir / "readiness.json", result)
    return result


__all__ = ["prepare_execution", "freeze_execution", "reopen_execution", "run_execution",
           "prediction_to_route", "prepare_h_ch4_campaign", "PROPOSED_CONTROL"]
