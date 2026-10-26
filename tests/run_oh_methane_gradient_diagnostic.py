"""Two fixed-gradient execution diagnostics; not saddle/rate verification.

Prepare from the stopped angular continuation, run serial and package-matched
MPI cases in fresh directories, then close with a numerical comparison. No
old job, shell configuration, MPI installation or production default is edited.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from src.orca_runner import find_orca, run_orca

spec = importlib.util.spec_from_file_location("gradient_frozen_pilot", Path(__file__).with_name("run_oh_methane_angular_audit.py"))
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)
state = pilot.state
CASES = ("serial", "matched-mpi")
E_TOL = 1e-8
G_TOL = 1e-7


def native_frame(text):
    blocks = re.findall(r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\s*\n(.*?)(?:\n\s*\n)", text, re.S)
    rows = [row.split() for row in blocks[-1].splitlines()] if blocks else []
    if len(rows) != 7 or tuple(row[0] for row in rows) != pilot.LABELS or any(len(row) != 4 for row in rows):
        raise ValueError("Seven-atom mapped native Cartesian frame required")
    xyz = np.array([[float(value) for value in row[1:]] for row in rows])
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite native Cartesian frame")
    return xyz


def printed_derivatives(text):
    energies = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?[\d.]+)", text)
    blocks = re.findall(r"CARTESIAN GRADIENT\s*\n-+\s*\n(.*?)Difference to translation invariance", text, re.S)
    rows = re.findall(r"(?m)^\s*(\d+)\s+([A-Z][a-z]?)\s*:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$", blocks[-1]) if blocks else []
    if (not energies or len(rows) != 7 or tuple(row[1] for row in rows) != pilot.LABELS
            or [int(row[0]) for row in rows] != list(range(1, 8))):
        raise ValueError("Complete mapped native energy/gradient table required")
    energy = float(energies[-1])
    gradient = np.array([[float(v) for v in row[2:]] for row in rows]).ravel()
    if not math.isfinite(energy) or not np.isfinite(gradient).all():
        raise ValueError("Nonfinite native derivatives")
    return energy, gradient


def input_text(xyz, cores):
    xyz = np.asarray(xyz, dtype=float)
    if xyz.shape != (7, 3) or not np.isfinite(xyz).all() or cores not in (1, 4):
        raise ValueError("Mapped fixed seven-atom geometry and one/four cores required")
    coordinates = "".join(f"{label} {x:.14f} {y:.14f} {z:.14f}\n"
                          for label, (x, y, z) in zip(pilot.LABELS, xyz, strict=True))
    return (f"! {pilot.METHOD} EnGrad\n%pal\n nprocs {cores}\nend\n%maxcore 1000\n"
            + pilot.SCF + "* xyz 0 2\n" + coordinates + "*\n")


def prepare(source):
    source = source.resolve(strict=True)
    if source.parent != (ROOT / "runs").resolve():
        raise ValueError("Owned repository run required")
    result = pilot.read_json(source / "result.json")
    if (result.get("status") != "unresolved" or result.get("old_sources_unchanged") is not True
            or result.get("whole_increment_attempts") != {"saddles": 5, "NumFreq": 2}
            or result.get("points", [])[-1].get("offset_deg") != 15.0):
        raise ValueError("Exact stopped plus-angle continuation required")
    previous = pilot.read_json(pilot.reopen(result["contract"]))
    state.reopen_inventory(previous["source_inventory"])
    opt = source / "point-02/optimization/opt.inp"
    execution = pilot.read_json(opt.parent / "execution.json")
    pilot.reopen(execution["input"])
    for record in execution["artifacts"]:
        pilot.reopen(record)
    native = opt.with_suffix(".out").read_text()
    if ("error termination in SCF gradient" not in native or "exiting\nimproperly" not in native
            or "ORCA TERMINATED NORMALLY" in native or native.count("GEOMETRY OPTIMIZATION CYCLE") != 6):
        raise ValueError("Exact cycle-six native MPI gradient failure required")
    if not opt.read_text().startswith(f"! {pilot.METHOD} OptTS VeryTightOpt\n") or pilot.SCF not in opt.read_text():
        raise ValueError("Retained chemical method/SCF settings differ")
    labels, xyz, _ = pilot.read_xyz(opt.with_suffix(".xyz"))
    if (tuple(labels) != pilot.LABELS or xyz.shape != (7, 3) or not np.isfinite(xyz).all()
            or not np.allclose(native_frame(native), xyz, atol=6e-7, rtol=0)):
        raise ValueError("Saved full-precision geometry differs from failed native frame")
    failed_energy, failed_gradient = printed_derivatives(native)
    orca = Path(find_orca()).resolve(strict=True)
    paths = {Path(path) for path in previous["source_inventory"]}
    paths.update(state.bound_paths(result))
    paths.update({source / "result.json", source / "contract.json", opt, opt.with_suffix(".out"),
                  opt.with_suffix(".xyz"), opt.with_suffix(".process.json"), opt.parent / "execution.json",
                  Path(__file__), Path(pilot.__file__), Path(state.__file__), ROOT / "src/orca_runner.py",
                  orca, orca.parent / "orca_scfgrad", orca.parent / "orca_scfgrad_mpi", orca.parent / "orca_leanscf", orca.parent / "orca_leanscf_mpi"})
    folder = Path(tempfile.mkdtemp(prefix="oh-methane-gradient-mpi-diagnostic-20261012-", dir=ROOT / "runs"))
    contract = {"schema_version": 1, "kind": "fixed_geometry_serial_matched_MPI_gradient_diagnostic",
                "source_inventory": state.inventory(paths), "source": pilot.artifact(source / "result.json"),
                "xyz_A": xyz.tolist(), "orca": pilot.artifact(orca), "cases": CASES,
                "per_case_attempt_cap": 1, "per_job_seconds": pilot.JOB_SECONDS, "minimum_disk_reserve_bytes": pilot.MIN_FREE,
                "method": pilot.METHOD, "SCF": pilot.SCF, "coordinates_source": pilot.artifact(opt.with_suffix(".xyz")),
                "failed_native_reference": {"energy_Eh": failed_energy, "gradient_Eh_bohr": failed_gradient.tolist(),
                                             "scope": "printed incomplete-job diagnostic only; not admitted derivatives"},
                "agreement_tolerances": {"energy_Eh": E_TOL, "max_gradient_component_Eh_bohr": G_TOL},
                "rate_claim_allowed": False, "lifetime_claim_allowed": False, "saddle_claim_allowed": False,
                "angular_budget_unchanged": {"saddles": 5, "NumFreq": 2}, "new_budget": {"EnGrad": 2},
                "scope": "fresh SCF comparisons at one fixed geometry, not literal replay of failed orbitals or proof of crash causation"}
    pilot.save(folder / "contract.json", contract)
    print("WORKDIR " + str(folder), flush=True)
    return folder


def matched_runtime(prefix):
    prefix = prefix.resolve(strict=True)
    env = {"PATH": str(prefix / "bin") + os.pathsep + os.environ.get("PATH", ""),
           "DYLD_LIBRARY_PATH": str(prefix / "lib") + os.pathsep + os.environ.get("DYLD_LIBRARY_PATH", ""),
           "DYLD_FALLBACK_LIBRARY_PATH": str(prefix / "lib") + os.pathsep + os.environ.get("DYLD_FALLBACK_LIBRARY_PATH", "")}
    launcher = prefix / "bin/mpirun"
    version = subprocess.check_output([str(launcher), "--version"], env={**os.environ, **env}, text=True)
    if not re.search(r"mpirun \(Open MPI\) 4\.1\.1\b", version):
        raise ValueError("Package-matched Open MPI 4.1.1 required")
    libs = sorted((prefix / "lib").glob("*.dylib"))
    if not (prefix / "lib/libmpi.40.dylib").is_file() or not libs:
        raise ValueError("Matched native MPI libraries missing")
    records = state.inventory({launcher, prefix / "bin/ompi_info", *libs})
    return env, {"prefix": str(prefix), "mpirun_version": version, "inventory": records, "selected_environment": env}


def run_case(folder, case, prefix=None):
    if case not in CASES or (case == "matched-mpi") != (prefix is not None):
        raise ValueError("Exactly serial or matched-MPI runtime selection required")
    contract = pilot.read_json(folder / "contract.json")
    state.reopen_inventory(contract["source_inventory"])
    if Path(find_orca()).resolve() != pilot.reopen(contract["orca"]):
        raise ValueError("Bound ORCA executable differs")
    target = folder / case
    if target.exists():
        raise ValueError("Case already dispatched; no automatic retry or overwrite")
    if shutil.disk_usage(folder).free < contract["minimum_disk_reserve_bytes"]:
        raise ValueError("Minimum disk reserve reached before dispatch")
    extra, runtime = matched_runtime(prefix) if prefix is not None else ({}, {"scope": "serial; no MPI launcher"})
    extra.update({name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")})
    extra["DYLD_PRINT_LIBRARIES"] = "1"  # Prove actual MPI library selection, not PATH alone.
    target.mkdir()
    inp = target / "gradient.inp"
    with inp.open("x") as stream:
        stream.write(input_text(contract["xyz_A"], 1 if case == "serial" else 4))
    pilot.save(target / "runtime.json", runtime)
    failure, receipt, resource_stop = None, None, []

    def progress(line):
        if re.search(r"FINAL SINGLE POINT|SCF CONVERGED|ORCA TERMINATED|error termination|aborting", line):
            print(line.rstrip(), flush=True)
        if shutil.disk_usage(folder).free < contract["minimum_disk_reserve_bytes"]:
            resource_stop.append("Minimum disk reserve reached")
            return True
        return False

    print("START " + str(inp), flush=True)
    try:
        run_orca(inp, timeout_seconds=contract["per_job_seconds"], extra_env=extra, stream_output=progress)
        if resource_stop:
            raise ValueError(resource_stop[0])
        native = inp.with_suffix(".out").read_text()
        lease = pilot.read_json(inp.with_suffix(".process.json"))
        if (lease["status"] != "terminated" or lease["returncode"] != 0
                or lease["input_sha256"] != pilot.artifact(inp)["sha256"]
                or native.count("ORCA TERMINATED NORMALLY") != 1):
            raise ValueError("Owned unique normal native completion required")
        energy, gradient = pilot.read_engrad(inp.with_suffix(".engrad"), np.asarray(contract["xyz_A"])/pilot.BOHR_A)
        printed_energy, printed_gradient = printed_derivatives(native)
        if (abs(energy-printed_energy) > 1e-12 or np.max(abs(gradient-printed_gradient)) > 5.1e-10
                or not np.allclose(native_frame(native), contract["xyz_A"], atol=6e-7, rtol=0)):
            raise ValueError("Native fixed-frame EnGrad/output disagreement")
        if prefix is not None:
            state.reopen_inventory(runtime["inventory"])
            expected = str((prefix / "lib/libmpi.40.dylib").resolve())
            loaded = re.findall(r"(?m)^dyld\[\d+\]:.*?(/\S*libmpi\.40\.dylib)\s*$", native)
            if not loaded or any(str(Path(path).resolve()) != expected for path in loaded):
                raise ValueError("Matched MPI library loading not proved or mixed runtime observed")
        else:
            loaded = []
        receipt = {"energy_Eh": energy, "gradient_Eh_bohr": gradient.tolist(), "native_state": pilot.native_state(inp),
                   "loaded_MPI_libraries": loaded, "fixed_geometry_derivative_job_completed": True}
    except BaseException as error:
        failure = f"{type(error).__name__}: {error}"
        print("DIAGNOSTIC UNRESOLVED: " + failure, flush=True)
    finally:
        try:
            state.reopen_inventory(contract["source_inventory"])
            unchanged = True
        except (OSError, ValueError) as error:
            unchanged = False
            failure = failure or str(error)
        output = {"case": case, "status": "completed_execution_diagnostic" if failure is None else "unresolved",
                  "failure_reason": failure, "receipt": receipt, "input": pilot.artifact(inp),
                  "artifacts": [pilot.artifact(path) for path in (inp.with_suffix(".out"), inp.with_suffix(".engrad"), inp.with_suffix(".process.json")) if path.is_file()],
                  "runtime": pilot.artifact(target / "runtime.json"), "old_sources_unchanged": unchanged,
                  "rate_claim_allowed": False, "lifetime_claim_allowed": False, "saddle_claim_allowed": False}
        pilot.save(target / "case-result.json", output)
    if failure:
        raise RuntimeError(failure)
    return output


def close(folder):
    contract = pilot.read_json(folder / "contract.json")
    state.reopen_inventory(contract["source_inventory"])
    cases = [pilot.read_json(folder / case / "case-result.json") for case in CASES]
    for case in cases:
        for record in [case["input"], case["runtime"], *case["artifacts"]]:
            pilot.reopen(record)
    comparison = None
    if all(case["status"] == "completed_execution_diagnostic" for case in cases):
        left, right = (case["receipt"] for case in cases)
        delta_E = abs(left["energy_Eh"]-right["energy_Eh"])
        delta_g = float(np.max(abs(np.asarray(left["gradient_Eh_bohr"])-np.asarray(right["gradient_Eh_bohr"]))))
        comparison = {"absolute_energy_difference_Eh": delta_E, "max_gradient_component_difference_Eh_bohr": delta_g,
                      "within_declared_engineering_tolerances": delta_E <= E_TOL and delta_g <= G_TOL,
                      "scope": "numerical reproducibility at one fixed frame; not physical uncertainty or root-cause proof"}
    result = {"status": "completed_comparison" if comparison and comparison["within_declared_engineering_tolerances"] else "unresolved",
              "contract": pilot.artifact(folder / "contract.json"), "cases": [pilot.artifact(folder / case / "case-result.json") for case in CASES],
              "comparison": comparison, "old_sources_unchanged": True, "native_attempts": {"EnGrad": 2},
              "crash_root_cause_proved": False, "rate_claim_allowed": False, "lifetime_claim_allowed": False, "saddle_claim_allowed": False}
    pilot.save(folder / "result.json", result)
    print(str(result), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "serial", "matched-mpi", "close"))
    parser.add_argument("directory", type=Path)
    parser.add_argument("mpi_prefix", type=Path, nargs="?")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.directory)
    elif args.action == "close":
        close(args.directory)
    else:
        run_case(args.directory.resolve(strict=True), args.action, args.mpi_prefix)
