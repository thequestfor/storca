"""Two owned fixed-geometry SCF diagnostics; no optimizer, gradient or rate claim.

Preparation is zero-native work. ``prepare --authorize`` records the user's
explicit consent for exactly two SP attempts in the EXISTING transfer ledger.
Failed attempts remain charged and are never automatically retried.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import nullcontext
import math
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from storca import rate_transfer_continuation as continuation
from storca import rate_transfer_execution as execution
from storca import rate_transfer_preparation as prep
from storca.rate_transfer_control import NativeExecutionControl, TransferControlError, controlled_run_or_resume
from storca.scientific_panel import canonical_hash
from src.orca_job_settings import active_native_settings, explicit_native_input_scope


KIND = "h_ch4_owned_two_arm_scf_diagnostic"
ARMS = ("original", "explicit-trah")
DOWNSTREAM_SP_RESERVE = 6  # Five matched energies plus one reached atomic point.
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def _inventory(prior_source_capsule=None):
    records = {row["path"]: row for row in continuation._implementation_bindings()}
    row = prep.binding(Path(__file__))
    records[row["path"]] = row
    if prior_source_capsule is not None:
        for relative in ("storca/partial_neb_recovery.py", "tests/run_h_ch4_partial_neb_recovery.py"):
            row = prep.binding(Path(__file__).resolve().parents[1] / relative)
            records[row["path"]] = row
    return [records[path] for path in sorted(records)]


def _context(workdir, prior_source_capsule=None):
    # Continuation reopening rechecks retained fragments under its OWN frozen
    # profile. During an owned run the same profile is already active. Suspend
    # only implicit rendering while this read-only revalidation opens that
    # scope, then restore it and require exact profile equality. The native
    # controller and exclusive ownership remain active throughout.
    active = active_native_settings()
    with explicit_native_input_scope():
        result = _read_context(workdir, prior_source_capsule)
    if active is not None and active.sha256 != result[2].sha256:
        raise TransferControlError("Active native settings differ from diagnostic frozen settings")
    return result


def _read_context(workdir, prior_source_capsule=None):
    if prior_source_capsule is None:
        sealed = continuation.reopen_continuation(workdir)
        if sealed["plan"]["authorized"] is not True:
            raise ValueError("Source continuation must have been explicitly authorized")
        historical, prediction, settings, control, snapshot = continuation._original_context(
            sealed["plan"]["original_campaign"], sealed["plan"]["source_capsule"])
    else:
        # Exact old source bytes, old raw plan, original model/settings/charges
        # are revalidated by the independently bound recovery archive reader.
        # No live-source mismatch is ignored and no ledger is replaced.
        from run_h_ch4_partial_neb_recovery import _archive
        _, historical, settings, control = _archive(
            prep.binding(Path(workdir) / "continuation-contract.json"), prior_source_capsule)
        prediction = execution._prediction(prep._json(prep._reopen(historical["plan"]["prediction"])))
        snapshot = control.snapshot()
    if (prediction["case_id"] != "h_ch4" or prediction["total_charge"] != 0
            or prediction["surface_multiplicity"] != 2):
        raise ValueError("This diagnostic is only for the existing neutral doublet H + methane case")
    resource = control.contract["resources"]
    if (resource["ncores"] != 4 or resource["maxcore_mb_per_rank"] != 1000
            or control.contract["limits"]["per_job_seconds"] != 21600):
        raise ValueError("Four cores, 1000 MB/rank and six-hour cap must remain unchanged")
    return historical, prediction, settings, control, snapshot


def _geometry(path):
    lines = path.read_text().splitlines()
    if not lines or lines[0].strip() != "6" or len(lines) != 8:
        raise ValueError("Exactly one retained six-nucleus XYZ frame is required")
    labels = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError("Malformed retained XYZ coordinate")
        labels.append(fields[0])
        if not all(math.isfinite(float(number)) for number in fields[1:]):
            raise ValueError("Nonfinite retained XYZ coordinate")
    if Counter(labels) != {"C": 1, "H": 5}:
        raise ValueError("Retained H + methane nucleus inventory differs")
    return "\n".join(lines[2:]) + "\n"


def _failed_source(snapshot):
    matches = []
    for attempt, row in snapshot["reservations"].items():
        if row["category"] != "optimization_path":
            continue
        inp = Path(row["input"])
        out = inp.with_suffix(".out")
        if not out.is_file():
            continue
        text = out.read_text(errors="replace")
        if "Diagonalization failed with error code: 203" in text:
            if not re.search(r"NEB.*CONVERGED|NEB.*converged|NEB convergence", text):
                raise ValueError("Expected retained converged NEB stage is absent")
            if "ORCA TERMINATED NORMALLY" in text:
                raise ValueError("Expected failed compound native job instead terminated normally")
            terminal = snapshot["terminals"].get(attempt)
            if terminal is None:
                raise TransferControlError("Failed native ownership is unresolved")
            xyz = inp.with_suffix(".xyz")
            _geometry(xyz)
            bound = {source["path"]: source for source in terminal["sources"]}
            required = [prep.binding(path) for path in (inp, out, xyz)]
            if any(bound.get(source["path"]) != source for source in required):
                raise ValueError("Failing geometry and native input/output must be owned terminal evidence")
            matches.append({"attempt_id": attempt, "input": required[0],
                            "output": required[1], "geometry": required[2],
                            "terminal_sources": terminal["sources"]})
    if len(matches) != 1:
        raise ValueError("Exactly one owned converged-NEB/203-failure source is required")
    return matches[0]


def _input(coordinates, arm):
    keywords = list(execution.LOW_KEYWORDS)
    if arm == "explicit-trah":
        keywords.append("TRAH")
    elif arm != "original":
        raise ValueError("Unsupported diagnostic arm")
    return ("! " + " ".join(keywords) + "\n%pal\n  nprocs 4\nend\n%maxcore 1000\n"
            "%scf\n  CNVSOSCF false\n  ConvCheckMode 0\n  ConvForced 1\n"
            "  MaxIter 500\nend\n* xyz 0 2\n" + coordinates + "*\n")


def _available(control, snapshot, jobs):
    remaining = control.contract["limits"]["single_point"] - snapshot["counts"]["single_point"]
    if remaining < jobs + DOWNSTREAM_SP_RESERVE:
        raise TransferControlError("Two-arm diagnostics must preserve six downstream single-point slots")


def _capsule_binding(value):
    if value is None:
        return None
    if isinstance(value, dict):
        prep._reopen(value)
        return value
    path = Path(value).resolve()
    return prep.binding(path / "source-capsule.json" if path.is_dir() else path)


def prepare(source_continuation, destination, *, authorized=False, prior_source_capsule=None):
    source_continuation, destination = Path(source_continuation).resolve(), Path(destination).resolve()
    capsule = _capsule_binding(prior_source_capsule)
    contract_path = destination / "diagnostic-contract.json"
    if contract_path.exists():
        value = inspect(destination)
        if (value["authorized"] != authorized or value["source_continuation"]["path"] != str(
                source_continuation / "continuation-contract.json")
                or value["prior_source_capsule"] != capsule):
            raise ValueError("Prepared authorization/source cannot be changed in place")
        return value
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise ValueError("Prepare into a fresh diagnostic directory")
    _, _, _, control, snapshot = _context(source_continuation, capsule)
    if destination.is_relative_to(source_continuation):
        raise ValueError("Diagnostic evidence must not be written inside historical continuation")
    _available(control, snapshot, 2)
    failed = _failed_source(snapshot)
    coordinates = _geometry(prep._reopen(failed["geometry"]))
    destination.mkdir(parents=True, exist_ok=True)
    arms = []
    for arm in ARMS:
        folder = destination / arm
        folder.mkdir()
        inp = folder / "scf.inp"
        # Immutable creation: inputs cannot adopt an existing GBW or output.
        with inp.open("x") as stream:
            stream.write(_input(coordinates, arm))
        arms.append({"arm": arm, "input": prep.binding(inp)})
    plan = {"schema_version": 1, "kind": KIND, "authorized": authorized,
            "source_continuation": prep.binding(source_continuation / "continuation-contract.json"),
            "prior_source_capsule": capsule,
            "source_workdir": str(source_continuation), "failed_source": failed,
            "native_control_contract": prep.binding(control.root / "contract.json"),
            "native_control_root": str(control.root), "native_settings": prep.binding(control.root / "native-settings.json"),
            "journal_prefix": [prep.binding(path) for path in sorted((control.root / "events").iterdir())],
            "initial_counts": snapshot["counts"], "arms": arms,
            "downstream_single_point_reserve": DOWNSTREAM_SP_RESERVE,
            "guess_policy": "fresh_NoAutoStart_no_old_GBW_same_default_guess_both_arms",
            "scope": "Fixed failing geometry; original SCF versus explicit TRAH only; no causal proof or TS/rate admission."}
    prep._publish(destination / "diagnostic-plan.json", plan)
    prep._publish(contract_path, {"schema_version": 1, "kind": KIND + "_sealed",
                  "plan": plan, "plan_source": prep.binding(destination / "diagnostic-plan.json"),
                  "implementation": _inventory(capsule)})
    return inspect(destination)


def _reopen(destination):
    destination = Path(destination).resolve()
    sealed = prep._json(destination / "diagnostic-contract.json")
    if (set(sealed) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or sealed["schema_version"] != 1 or sealed["kind"] != KIND + "_sealed"):
        raise ValueError("Invalid diagnostic contract")
    plan = prep._json(prep._reopen(sealed["plan_source"]))
    if canonical_hash(plan) != canonical_hash(sealed["plan"]):
        raise ValueError("Diagnostic plan changed")
    for binding in sealed["implementation"]:
        prep._reopen(binding)
    if sealed["implementation"] != _inventory(plan.get("prior_source_capsule")):
        raise ValueError("Diagnostic implementation inventory changed")
    prep._reopen(plan["source_continuation"])
    historical, prediction, settings, control, snapshot = _context(
        plan["source_workdir"], plan.get("prior_source_capsule"))
    if (plan["schema_version"] != 1 or plan["kind"] != KIND or type(plan["authorized"]) is not bool
            or plan["downstream_single_point_reserve"] != DOWNSTREAM_SP_RESERVE
            or plan["native_control_root"] != str(control.root)):
        raise ValueError("Diagnostic scope or same-ledger binding changed")
    for key, name in (("native_control_contract", "contract.json"), ("native_settings", "native-settings.json")):
        if plan[key] != prep.binding(control.root / name):
            raise ValueError("Diagnostic cannot replace the original native ledger/settings")
    current = sorted((control.root / "events").iterdir())
    prefix = plan["journal_prefix"]
    if [prep.binding(path) for path in current[:len(prefix)]] != prefix:
        raise ValueError("Historical native journal prefix changed")
    if _failed_source(snapshot) != plan["failed_source"]:
        raise ValueError("Failed native source changed")
    coordinates = _geometry(prep._reopen(plan["failed_source"]["geometry"]))
    if len(plan["arms"]) != 2 or [row["arm"] for row in plan["arms"]] != list(ARMS):
        raise ValueError("Exactly the two declared diagnostic arms are required")
    for row in plan["arms"]:
        inp = prep._reopen(row["input"])
        if inp != destination / row["arm"] / "scf.inp" or inp.read_text() != _input(coordinates, row["arm"]):
            raise ValueError("Diagnostic input differs from the controlled comparison")
    return plan, historical, settings, control, snapshot


def read_output(path):
    """Read electronic completion, not exit codes or thermochemical validity."""
    path = Path(path)
    text = path.read_text(errors="replace") if path.is_file() else ""
    energies = re.findall(r"FINAL SINGLE POINT ENERGY\s+(" + NUMBER + r")", text)
    spins = re.findall(r"Expectation value of <S\*\*2>\s*:\s*(" + NUMBER + r")", text)
    number = lambda value: float(value.replace("D", "E").replace("d", "e"))
    energy = number(energies[-1]) if energies else None
    spin = number(spins[-1]) if spins else None
    failures = [line.strip() for line in text.splitlines() if re.search(
        r"Diagonalization failed|SCF NOT CONVERGED|SCF did not converge|error termination|OUT OF MEMORY|MPI_ABORT", line, re.I)]
    normal = "ORCA TERMINATED NORMALLY" in text
    converged = bool(re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", text))
    def last(pattern, convert=float):
        found = re.findall(pattern, text)
        return convert(found[-1]) if found else None
    passed = normal and converged and energy is not None and spin is not None and not failures \
        and math.isfinite(energy) and math.isfinite(spin)
    return {"status": "electronic_solve_completed_not_TS_validation" if passed else "electronic_solve_unresolved",
            "normal_termination": normal, "scf_converged": converged,
            "energy_hartree": energy, "S_squared": spin, "reported_failures": failures,
            "charge": last(r"Total Charge\s+Charge\s+\.{2,}\s+([-+]?\d+)", int),
            "multiplicity": last(r"Multiplicity\s+Mult\s+\.{2,}\s+(\d+)", int),
            "total_electrons": last(r"Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)", int),
            "alpha_electrons": last(r"N\(Alpha\)\s*:\s*(" + NUMBER + r")", number),
            "beta_electrons": last(r"N\(Beta\)\s*:\s*(" + NUMBER + r")", number),
            "output": prep.binding(path) if path.is_file() else None,
            "physical_accuracy_claim_allowed": False, "rate_claim_allowed": False,
            "crash_cause_established": False}


def _owned(snapshot, input_path):
    rows = [row for row in snapshot["reservations"].values() if row["input"] == str(input_path)]
    if len(rows) > 1:
        raise TransferControlError("Diagnostic arm was dispatched more than once")
    if not rows:
        return None
    reservation = rows[0]
    terminal = snapshot["terminals"].get(reservation["attempt_id"])
    if terminal is None:
        raise TransferControlError("Owned diagnostic work remains pending; do not duplicate")
    return {"reservation": reservation, "terminal": terminal}


def read_comparison(destination):
    """Reopen exact report/native bytes without claiming live source equivalence.

    This is a pure evidence reader for a later separately source-bound recovery.
    Its caller MUST independently establish the archived diagnostic implementation
    lineage. It does not authorize work, select a solver, or waive source guards
    for ``run``. The old implementation inventory is returned, not replaced.
    """
    destination = Path(destination).resolve()
    sealed = prep._json(destination / "diagnostic-contract.json")
    if (set(sealed) != {"schema_version", "kind", "plan", "plan_source", "implementation"}
            or sealed["schema_version"] != 1 or sealed["kind"] != KIND + "_sealed"):
        raise ValueError("Invalid diagnostic contract")
    plan = prep._json(prep._reopen(sealed["plan_source"]))
    if canonical_hash(plan) != canonical_hash(sealed["plan"]):
        raise ValueError("Diagnostic plan changed")
    control_contract = prep._json(prep._reopen(plan["native_control_contract"]))
    control = NativeExecutionControl(plan["native_control_root"], control_contract)
    if plan["native_control_contract"] != prep.binding(control.root / "contract.json"):
        raise ValueError("Diagnostic control belongs to another ledger")
    snapshot = control.snapshot()
    prefix = plan["journal_prefix"]
    current = sorted((control.root / "events").iterdir())
    if [prep.binding(path) for path in current[:len(prefix)]] != prefix:
        raise ValueError("Historical diagnostic journal prefix changed")
    failed = plan["failed_source"]
    failure_reservation = snapshot["reservations"].get(failed["attempt_id"])
    failure_terminal = snapshot["terminals"].get(failed["attempt_id"])
    if (failure_reservation is None or failure_terminal is None
            or failure_reservation["input"] != failed["input"]["path"]
            or failure_terminal["sources"] != failed["terminal_sources"]):
        raise ValueError("Diagnostic source does not match original owned failed attempt")
    sources = [prep.binding(destination / "diagnostic-contract.json"), sealed["plan_source"],
               plan["source_continuation"], plan["native_control_contract"], plan["native_settings"]]
    if plan.get("prior_source_capsule") is not None:
        sources.append(plan["prior_source_capsule"])
    for source in sources + plan["failed_source"]["terminal_sources"]:
        prep._reopen(source)
    coordinates = _geometry(prep._reopen(plan["failed_source"]["geometry"]))
    if len(plan["arms"]) != 2 or [row["arm"] for row in plan["arms"]] != list(ARMS):
        raise ValueError("Exactly the two declared diagnostic arms are required")
    results = []
    for row in plan["arms"]:
        inp = prep._reopen(row["input"])
        if inp != destination / row["arm"] / "scf.inp" or inp.read_text() != _input(coordinates, row["arm"]):
            raise ValueError("Controlled diagnostic input changed")
        owned = _owned(snapshot, inp)
        if owned is None:
            raise ValueError("Diagnostic report has no owned attempt")
        control.validate_reuse(inp, owned["reservation"]["input_sha256"])
        result_path = inp.parent / "result.json"
        result = prep._json(result_path)
        if (result["schema_version"] != 1 or result["arm"] != row["arm"]
                or result["input"] != row["input"] or result["owned_attempt"] != owned
                or result["assessment"] != read_output(inp.with_suffix(".out"))):
            raise ValueError("Diagnostic receipt or native assessment changed")
        sources.extend([prep.binding(result_path), row["input"], *owned["terminal"]["sources"]])
        results.append(result)
    report_path = destination / "comparison.json"
    report = prep._json(report_path)
    last_attempt = max(result["owned_attempt"]["reservation"]["attempt_id"] for result in results)
    counts = {category: sum(row["category"] == category and key <= last_attempt
                           for key, row in snapshot["reservations"].items())
              for category in snapshot["counts"]}
    if report != _comparison(destination, results, counts):
        raise ValueError("Comparison differs from exact owned raw evidence")
    sources.extend([prep.binding(report_path), *plan["failed_source"]["terminal_sources"]])
    unique = {row["path"]: row for row in sources}
    return {"comparison": report, "assessments": [row["assessment"] for row in results],
            "plan": plan, "sources": [unique[path] for path in sorted(unique)],
            "historical_implementation": sealed["implementation"],
            "archived_implementation_lineage_must_be_validated_by_caller": True}


def _comparison(destination, results, counts):
    a, b = [result["assessment"] for result in results]
    comparable = all(item["status"] == "electronic_solve_completed_not_TS_validation" for item in (a, b))
    return {"schema_version": 1, "status": "two_arm_diagnostic_complete",
            "arms": [prep.binding(destination / arm / "result.json") for arm in ARMS],
            "both_electronic_solves_completed": comparable,
            "TRAH_minus_original_energy_hartree": b["energy_hartree"] - a["energy_hartree"] if comparable else None,
            "TRAH_minus_original_S_squared": b["S_squared"] - a["S_squared"] if comparable else None,
            "interpretation": "Fresh-guess fixed-geometry evidence only; does not reproduce lost optimizer DIIS history or prove original crash cause.",
            "crash_cause_established": False, "verified_TS": False,
            "rate_claim_allowed": False, "counts": counts}


def inspect(destination):
    plan, _, _, control, snapshot = _reopen(destination)
    return {"status": "diagnostic_bindings_reopened", "authorized": plan["authorized"],
            "source_continuation": plan["source_continuation"], "native_control_root": str(control.root),
            "prior_source_capsule": plan.get("prior_source_capsule"),
            "counts": snapshot["counts"], "pending": len(snapshot["pending"]),
            "native_jobs_launched": 0, "rate_claim_allowed": False}


def run(destination):
    destination = Path(destination).resolve()
    plan, historical, settings, control, _ = _reopen(destination)
    if plan["authorized"] is not True:
        raise ValueError("Diagnostic comparison has not been separately authorized")
    if (destination / "comparison.json").exists():
        return read_comparison(destination)["comparison"]
    from storca.rate_transfer_runtime import runtime_environment
    runtime = runtime_environment(historical["plan"]["runtime_inventory"]) \
        if historical["plan"]["schema_version"] == 2 else nullcontext()
    with runtime, control.activate(settings):
        results = []
        for index, row in enumerate(plan["arms"]):
            _reopen(destination)
            inp = prep._reopen(row["input"])
            result_path = inp.parent / "result.json"
            owned = _owned(control.snapshot(), inp)
            if result_path.exists():
                result = prep._json(result_path)
                if owned is None or result["owned_attempt"] != owned:
                    raise ValueError("Diagnostic result lacks exact owned native evidence")
                if result["assessment"] != read_output(inp.with_suffix(".out")):
                    raise ValueError("Retained diagnostic output assessment changed")
                results.append(result)
                continue
            exception = None
            if owned is None:
                if any(inp.parent.glob("*.gbw")) or inp.with_suffix(".out").exists():
                    raise TransferControlError("Unowned orbital/output evidence cannot seed a fresh diagnostic")
                _available(control, control.snapshot(), len(ARMS) - index)
                try:
                    controlled_run_or_resume(inp, timeout_seconds=21600,
                                             explicit_single_point_settings=True)
                except Exception as error:
                    exception = {"type": type(error).__name__, "message": str(error)}
                    if _owned(control.snapshot(), inp) is None:
                        return {"status": "blocked_before_charge", "error": exception,
                                "completed_arms": len(results), "counts": control.snapshot()["counts"],
                                "rate_claim_allowed": False}
                owned = _owned(control.snapshot(), inp)
            if owned is None:
                raise TransferControlError("Diagnostic dispatcher returned without owned native work")
            result = {"schema_version": 1, "arm": row["arm"], "input": row["input"],
                      "owned_attempt": owned, "assessment": read_output(inp.with_suffix(".out")),
                      "dispatch_exception": exception}
            prep._publish(result_path, result)
            results.append(result)
        _reopen(destination)
        comparison = _comparison(destination, results, control.snapshot()["counts"])
        prep._publish(destination / "comparison.json", comparison)
        return comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("prepare")
    create.add_argument("source_continuation", type=Path)
    create.add_argument("destination", type=Path)
    create.add_argument("--authorize", action="store_true", help="record explicit consent for exactly two SP attempts")
    create.add_argument("--prior-source-capsule", type=Path,
                        help="exact archived prior continuation implementation, never a source-check waiver")
    for name in ("inspect", "run"):
        commands.add_parser(name).add_argument("destination", type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        value = prepare(args.source_continuation, args.destination, authorized=args.authorize,
                        prior_source_capsule=args.prior_source_capsule)
    else:
        value = (inspect if args.command == "inspect" else run)(args.destination)
    import json
    print(json.dumps(value, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
