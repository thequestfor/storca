"""Fresh, bounded constrained-saddle pilot; never a production rate correction.

The source argument is the completed electronic-state replay directory. Native
ORCA runs use the retained electronic surface, not a fitted angular potential.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import os
from pathlib import Path
import re
import shutil
import signal
import sys
import tempfile
import threading
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np

from src.orca_runner import find_orca, run_orca
from storca.ir_modes import parse_orca_normal_modes
from storca.route_geometry import read_xyz
from storca.surface_resolution import bond_change_mode_overlap


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


state = load("angular_state_source_reader", Path(__file__).with_name("replay_oh_methane_electronic_states.py"))
artifact, reopen, read_json, save = state.artifact, state.reopen, state.read_json, state.write_new
LABELS = ("O", "H", "C", "H", "H", "H", "H")
ATOMIC_NUMBERS = (8, 1, 6, 1, 1, 1, 1)
DIHEDRAL = (1, 0, 2, 4)
# ORCA6.1's native Cartesian conversion, checked against the bound input/Hessian.
# Keep this distinct from exact-SI constants used for energy/frequency units.
BOHR_A = 0.5291772083
EH_J_MOL = 4.3597447222071e-18 * 6.02214076e23
CM_J_MOL = 6.62607015e-34 * 299792458.0 * 100 * 6.02214076e23
JOB_SECONDS = 21600
MIN_FREE = 4 * 1024**3
PILOT_OFFSETS_DEG = (0.0, -15.0, 15.0)
# Append-only predecessors: no SCF, neither a chemistry conclusion nor a waiver.
PRIOR_FAILURES = (("oh-methane-angular-pilot-20261011-r91v4lhr", "Scan_TSMode",
                   "TS_Mode subsection terminator omitted"),
                  ("oh-methane-angular-pilot-20261011-479rd0_z", "oob_tcp_component.c",
                   "Sandbox blocked native local MPI startup"))
REUSE_CENTER = ROOT / "runs/oh-methane-angular-pilot-20261011-xqavjjg0/point-00/optimization/opt.inp"
METHOD = "UKS WB97X-D3 def2-TZVPP NoAutoStart ExtremeSCF DefGrid3 RIJONX def2/J"
SCF = "%scf\n CNVSOSCF false\n ConvCheckMode 0\n ConvForced 1\n MaxIter 500\nend\n"
TRANSFER = ({"atom_indices": [2, 6], "kind": "broken"},
            {"atom_indices": [0, 6], "kind": "formed"})


def inline_geometry(path):
    text = Path(path).read_text()
    blocks = re.findall(r"(?im)^\s*\*\s*xyz\s+0\s+2\s*\n(.*?)^\s*\*\s*$", text, re.S)
    if len(blocks) != 1:
        raise ValueError("Exactly one inline neutral-doublet geometry required")
    fields = [row.split() for row in blocks[0].splitlines() if row.strip()]
    if tuple(row[0] for row in fields) != LABELS or any(len(row) != 4 for row in fields):
        raise ValueError("Retained mapped atom order differs")
    xyz = np.array([[float(v) for v in row[1:]] for row in fields])
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite seed coordinates")
    return xyz


def input_text(xyz_A, task, hessian=None):
    xyz_A = np.asarray(xyz_A, dtype=float)
    if xyz_A.shape != (7, 3) or not np.isfinite(xyz_A).all():
        raise ValueError("Finite seven-atom mapped seed required")
    if task == "opt":
        if not hessian:
            raise ValueError("Bound numerical initial Hessian required")
        extra = ("%geom\n CoordSys redundant\n CartFallback false\n InHess Read\n"
                 f' InHessName "{hessian}"\n TS_Mode {{ B 2 6 }} end\n'
                 " Constraints\n  { D 1 0 2 4 C }\n end\n MaxIter 100\nend\n")
        keyword = "OptTS VeryTightOpt"
    elif task == "freq":
        if hessian is not None:
            raise ValueError("Full NumFreq must not inherit constraints/Hessian")
        extra = ("%freq\n CentralDiff true\n Increment 0.001\n Restart false\n"
                 " ScalFreq 1.0\n ProjectTR true\n QuasiRRHO false\n Temp 340.3\nend\n")
        keyword = "NumFreq"
    else:
        raise ValueError("Only bounded OptTS and full fixed-geometry NumFreq tasks")
    coordinates = "".join(f"{label} {x:.14f} {y:.14f} {z:.14f}\n"
                          for label, (x, y, z) in zip(LABELS, xyz_A, strict=True))
    return (f"! {METHOD} {keyword}\n%pal\n nprocs 4\nend\n%maxcore 1000\n"
            + SCF + extra + "* xyz 0 2\n" + coordinates + "*\n")


def load_source(source):
    source = source.resolve(strict=True)
    if source.parent != (ROOT / "runs").resolve():
        raise ValueError("Source must be a retained repository runs child")
    result = read_json(source / "result.json")
    if (result.get("status") != "completed_conditional_state_diagnostic"
            or result.get("retained_artifacts_unchanged") is not True
            or result.get("network_reinjection_allowed") is not False):
        raise ValueError("Closed conditional electronic-state replay required")
    contract = read_json(reopen(result["contract"]))
    state.reopen_inventory(contract["retained_source_inventory"])
    state.validate_native_receipt(contract, read_json(reopen(result["native_replay"])))
    state.load_reporter(Path(contract["source_workdir"]))  # Read-only source closure.
    audit = read_json(reopen(contract["retained_audits"]["composite"]))
    recipe = read_json(Path(contract["source_workdir"]) / "composite/composite-rate-recipe.json")
    point = next(point for point in recipe["points"] if point["label"] == "TS")
    freq_input = reopen(point["frequency_input"])
    hessian = reopen(point["hessian"])
    if hessian.parent != freq_input.parent or hessian.name != "freq.hess":
        raise ValueError("Exact retained numerical-Hessian source required")
    geometry_folder = hessian.parents[3]
    admission = read_json(geometry_folder / "numerical-continuation-result.json")
    if (admission.get("status") != "completed_numerical_route_diagnostic"
            or admission.get("bound_source_bytes_and_mtime_unchanged") is not True
            or admission.get("prior_scientific_artifacts_unchanged") is not True
            or admission.get("rate_claim_allowed") is not False):
        raise ValueError("Closed conditional numerical route source required")
    matrix_reader = load("angular_retained_matrix_reader", geometry_folder / "matrix_replay.py")
    H, masses, xyz_bohr = matrix_reader.read_matrix(artifact(hessian))
    xyz = inline_geometry(freq_input)
    modes = parse_orca_normal_modes(hessian)
    if (modes.atom_labels != LABELS or H.shape != (21, 21)
            or not np.allclose(xyz / BOHR_A, xyz_bohr, atol=1e-9, rtol=0)
            or not np.allclose(modes.coordinates_bohr, xyz_bohr, atol=1e-9, rtol=0)
            or len(modes.frequencies_cm_1) != 21
            or np.count_nonzero(modes.frequencies_cm_1[6:] < 0) != 1):
        raise ValueError("Retained geometry/state/Hessian/mode source mismatch")
    return {"source": source, "state_contract": contract, "frequency_input": freq_input,
            "hessian": hessian, "geometry_folder": geometry_folder, "xyz_A": xyz,
            "masses_amu": masses, "H": H, "modes": modes,
            "matrix_reader": matrix_reader, "temperature_grid_K": audit["evaluation_temperature_grid_K"]}


def build_contract(origin, folder):
    from storca import ts_angular_audit as angular
    # Old inputs and helper versions stay bound; no old experiment is modified.
    paths = set(Path(path) for path in origin["state_contract"]["retained_source_inventory"])
    paths.update(origin["source"].glob("*.json"))
    paths.update({Path(__file__), Path(state.__file__), Path(angular.__file__),
                  ROOT / "src/orca_runner.py", ROOT / "storca/ir_modes.py",
                  ROOT / "storca/route_geometry.py", ROOT / "storca/surface_resolution.py",
                  origin["geometry_folder"] / "matrix_replay.py",
                  origin["geometry_folder"] / "numerical-continuation-contract.json",
                  origin["geometry_folder"] / "numerical-continuation-result.json",
                  origin["hessian"], origin["frequency_input"],
                  origin["frequency_input"].with_suffix(".engrad"), Path(find_orca())})
    predecessors = []
    for name, marker, reason in PRIOR_FAILURES:
        previous = ROOT / "runs" / name
        prior_result = read_json(previous / "result.json")
        prior_output = previous / "point-00/optimization/opt.out"
        prior_native = prior_output.read_text()
        if (prior_result["attempts"] != {"saddles": 1, "NumFreq": 0}
                or marker not in prior_native or "SCF ITERATIONS" in prior_native):
            raise ValueError("Predecessor no longer matches its one pre-SCF attempt")
        paths.update({previous / "result.json", previous / "contract.json", prior_output,
                      prior_output.with_suffix(".inp"), prior_output.with_suffix(".process.json"),
                      prior_output.parent / "execution.json"})
        predecessors.append({"result": artifact(previous / "result.json"), "reason": reason,
                             "saddles": 1, "NumFreq": 0, "SCF_started": False})
    center_execution = read_json(REUSE_CENTER.parent / "execution.json")
    center_result = read_json(REUSE_CENTER.parents[2] / "result.json")
    reopen(center_execution["input"])
    for record in center_execution["artifacts"]:
        reopen(record)
    if (center_execution["failure"] is not None or center_result["attempts"] != {"saddles": 1, "NumFreq": 0}
            or center_result["failure_reason"] != "ValueError: Constrained geometry all-five-criteria gate unresolved"
            or REUSE_CENTER.read_text() != input_text(origin["xyz_A"], "opt", origin["hessian"])):
        raise ValueError("Completed center optimization does not match exact requested input")
    paths.update({REUSE_CENTER, REUSE_CENTER.with_suffix(".out"), REUSE_CENTER.with_suffix(".xyz"),
                  REUSE_CENTER.with_suffix(".process.json"), REUSE_CENTER.parent / "execution.json",
                  REUSE_CENTER.parents[2] / "result.json", REUSE_CENTER.parents[2] / "contract.json"})
    predecessors.append({"result": artifact(REUSE_CENTER.parents[2] / "result.json"),
                         "reason": "Native overachieved-gradient convergence; all-five audit refused candidate before derivative validation",
                         "saddles": 1, "NumFreq": 0, "SCF_started": True})
    records = state.inventory(paths)
    snapshots = []
    destination = folder / "source-snapshots"
    destination.mkdir()
    for i, record in enumerate(records.values()):
        path = reopen(record)
        if path.suffix == ".py":
            saved = destination / f"{i}-{path.name}"
            with saved.open("xb") as handle:
                handle.write(path.read_bytes())
            snapshots.append({"source": record, "snapshot": artifact(saved)})
    return {"schema_version": 1, "kind": "bounded_OH_CH4_TS_angular_saddle_pilot",
            "source_inventory": records, "source_snapshots": snapshots,
            "source_replay": str(origin["source"]), "initial_hessian": artifact(origin["hessian"]),
            "orca": artifact(find_orca()), "mapped_labels": LABELS,
            "atom_mapping": {"OH_O": 0, "OH_H": 1, "CH4_C": 2,
                             "spectator_H": [3, 4, 5], "transferring_H": 6},
            "spectator_permutations": "none within pilot; domain quotient/counting remains unresolved",
            "dihedral_atom_indices": DIHEDRAL, "coordinate_units": "radians internally; seed offsets in degrees",
            "coordinate": "signed H_OH-O-C-Hspect4 dihedral; nonbonded O-C geometric axis",
            "periodic_domain_rad": [-math.pi, math.pi], "domain_coverage_in_pilot": False,
            "pilot_offsets_deg": PILOT_OFFSETS_DEG, "method": METHOD,
            "native_ORCA_bohr_to_angstrom": BOHR_A,
            "cores": 4, "maxcore_MB_per_rank": 1000, "per_job_seconds": JOB_SECONDS,
            "minimum_disk_reserve_bytes": MIN_FREE,
            "caps": {"whole_increment_saddle_attempts": 17, "whole_increment_NumFreq_jobs": 4,
                     "pilot_saddle_attempts": 3, "pilot_NumFreq_jobs": 3, "automatic_retries": 0},
            "prior_attempts": {"saddles": len(predecessors), "NumFreq": 0, "predecessors": predecessors},
            "policy": {"native_geometry_all_five_criteria_reported_separately": True,
                       "overachieved_gradient_candidate_requires_fixed_derivative_validation": True,
                       "overachieved_gradient_fraction_of_tolerance": 0.05,
                       "overachieved_max_step_tolerance_multiplier": 2.0,
                       "angular_tolerance_rad": 2e-4,
                       "tangent_gradient_RMS_Eh_bohr": 8e-6, "tangent_gradient_max_Eh_bohr": 3e-5,
                       "negative_tangent_eigenvalues": 1, "positive_tangent_eigenvalues": 13,
                       "minimum_transfer_mode_overlap": 0.90, "maximum_S_squared_change": 0.05,
                       "maximum_orthogonal_ZPE_shift_J_mol": 500.0,
                       "maximum_orthogonal_partition_log_shift": math.log(1.20)},
            "metric_scope": "local inverse mass metric; not a proved reduced relaxed-path kinetic metric",
            "reused_center_optimization": artifact(REUSE_CENTER),
            "policy_revision": "All-five printed criteria are not falsely claimed for native alternate convergence. Both candidate routes require independently fixed geometry/gradient and Lagrangian-Hessian saddle validation; no acceptance by marker alone.",
            "SCF_scope": "native ConvCheckMode0/ConvForced1; not proof of every printed SCF residual inequality",
            "rate_claim_allowed": False, "lifetime_claim_allowed": False,
            "production_defaults_changed": False, "network_reinjection_allowed": False,
            "manual_sources": ["https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations_TS.html",
                               "https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html#constrained-optimization",
                               "https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/geomopt.html"]}


def execute(folder, stem, text, attempts):
    if stem not in ("opt", "freq"):
        raise ValueError("Unknown bounded native task")
    key, cap = ("saddles", 3) if stem == "opt" else ("NumFreq", 3)
    if attempts[key] >= cap:
        raise ValueError("Pilot attempt cap reached")
    if shutil.disk_usage(folder.parent).free < MIN_FREE:
        raise ValueError("Disk reserve reached before dispatch")
    folder.mkdir(exist_ok=False)
    inp = folder / f"{stem}.inp"
    with inp.open("x") as handle:
        handle.write(text)
    input_hash = artifact(inp)["sha256"]
    done = threading.Event()
    resource_stop = []

    def watch():
        while not done.wait(2):
            if shutil.disk_usage(folder).free >= MIN_FREE:
                continue
            resource_stop.append("Owned job reached minimum disk reserve")
            try:
                lease = read_json(inp.with_suffix(".process.json"))
                if (lease["status"] == "running" and lease["input_sha256"] == input_hash
                        and Path(lease["input"]) == inp.resolve()
                        and os.getpgid(lease["pid"]) == lease["process_group_id"]):
                    os.killpg(lease["process_group_id"], signal.SIGINT)
            except (OSError, ValueError, KeyError):
                pass
            return

    def progress(line):
        if re.search(r"OPTIMIZATION HAS CONVERGED|FINAL SINGLE POINT|ORCA TERMINATED|"
                     r"VIBRATIONAL FREQUENCIES|GEOMETRY OPTIMIZATION CYCLE|ERROR|error termination", line):
            print(line.rstrip(), flush=True)
        return bool(resource_stop)

    attempts[key] += 1  # Dispatch failures also consume the finite attempt budget.
    started = time.monotonic()
    print("START " + str(inp), flush=True)
    watcher = threading.Thread(target=watch, daemon=True)
    watcher.start()
    failure = None
    try:
        run_orca(inp, stream_output=progress, timeout_seconds=JOB_SECONDS,
                 extra_env={key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")})
        if resource_stop:
            raise ValueError(resource_stop[0])
        lease = read_json(inp.with_suffix(".process.json"))
        native = inp.with_suffix(".out").read_text()
        if (lease["returncode"] != 0 or lease["input_sha256"] != input_hash
                or native.count("ORCA TERMINATED NORMALLY") != 1):
            raise ValueError("Owned unique native completion required")
    except BaseException as error:
        failure = f"{type(error).__name__}: {error}"
        raise
    finally:
        done.set()
        watcher.join(timeout=3)
        save(folder / "execution.json", {"input": artifact(inp), "failure": failure,
             "elapsed_seconds": time.monotonic() - started,
             "artifacts": [artifact(path) for path in (inp.with_suffix(".out"), inp.with_suffix(".process.json")) if path.is_file()],
             "rate_claim_allowed": False, "lifetime_claim_allowed": False})
    return inp


def native_state(inp):
    text = inp.with_suffix(".out").read_text()
    for regex, required in ((r"Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)", 0),
                            (r"Multiplicity\s+Mult\s+\.{2,}\s+(\d+)", 2),
                            (r"Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)", 19)):
        values = re.findall(regex, text)
        if not values or any(int(value) != required for value in values):
            raise ValueError("Native electronic state/atom count differs")
    for pattern in (r"Exchange Functional\s+Exchange\s+\.{2,}\s+WB97X-D3\b",
                    r"Correlation Functional\s+Correlation\s+\.{2,}\s+WB97X-D3\b",
                    r"Hartree-Fock type\s+HFTyp\s+\.{2,}\s+UHF\b",
                    r"Angular Grid \(max\. ang\.\)\s+AngularGrid\s+\.{2,}\s+6\b",
                    r"General Integration Accuracy\s+IntAcc\s+\.{2,}\s+4\.959\b"):
        if not re.search(pattern, text):
            raise ValueError("Native method/ground-state grid differs")
    for phrase in ("SCF CONVERGED AFTER", "USING zero damping", "COSX micro iterations for RIJONX calc"):
        if phrase not in text:
            raise ValueError("Native method/convergence evidence missing: " + phrase)
    spins = [float(value) for value in re.findall(r"Expectation value of <S\*\*2>\s*:\s*([\d.]+)", text)]
    if not spins or any(not 0.70 <= value <= 0.85 for value in spins):
        raise ValueError("Native doublet spin continuity warning")
    return {"charge": 0, "multiplicity": 2, "electrons": 19, "S_squared": spins,
            "continuity_scope": "same declared state and small spin diagnostic change; not an orbital-state proof"}


def optimization_receipt(inp, seed_A, policy):
    from storca.ts_angular_audit import signed_dihedral, wrap_angle
    text = inp.with_suffix(".out").read_text()
    if text.count("THE OPTIMIZATION HAS CONVERGED") != 1:
        raise ValueError("Constrained native optimization did not uniquely converge")
    tables = re.findall(r"\|Geometry convergence\|(.*?)(?:\.\.\.\.\.\.|Max\(Bonds\))", text, re.S)
    if not tables:
        raise ValueError("Native five-criterion convergence table missing")
    rows = re.findall(r"(Energy change|RMS gradient|MAX gradient|RMS step|MAX step)\s+"
                      r"(-?[\d.]+)\s+([\d.]+)\s+(YES|NO)", tables[-1])
    if len(rows) not in (4, 5) or len({row[0] for row in rows}) != len(rows):
        raise ValueError("Malformed native geometry convergence inventory")
    all_five = len(rows) == 5 and all(row[3] == "YES" and abs(float(row[1])) <= float(row[2]) for row in rows)
    criteria = {name: (abs(float(value)), float(tolerance), passed) for name, value, tolerance, passed in rows}
    alternative = ("The gradient convergence is overachieved" in text
                   and "Convergence will therefore be signaled now" in text
                   and all(name in criteria and criteria[name][2] == "YES"
                           and criteria[name][0] <= policy["overachieved_gradient_fraction_of_tolerance"]*criteria[name][1]
                           for name in ("RMS gradient", "MAX gradient"))
                   and all(name in criteria and criteria[name][0] <= policy["overachieved_max_step_tolerance_multiplier"]*criteria[name][1]
                           for name in ("RMS step", "MAX step")))
    if not all_five and not alternative:
        raise ValueError("Native constrained convergence candidate gate unresolved")
    # Full reversal preserves the signed dihedral; ORCA canonicalizes primitives.
    primitive = r"D\((?:H\s+1,O\s+0,C\s+2,H\s+4|H\s+4,C\s+2,O\s+0,H\s+1)\)[^\n]*\bC\s*$"
    if not re.search(primitive, text, re.M):
        raise ValueError("Native constrained mapped dihedral primitive not found")
    labels, xyz, _ = read_xyz(inp.with_suffix(".xyz"))
    if tuple(labels) != LABELS:
        raise ValueError("Optimized mapped atom order changed")
    gap = abs(wrap_angle(signed_dihedral(xyz, DIHEDRAL) - signed_dihedral(seed_A, DIHEDRAL)))
    if gap > policy["angular_tolerance_rad"]:
        raise ValueError("Constrained dihedral not preserved")
    return xyz, {"input": artifact(inp), "output": artifact(inp.with_suffix(".out")),
                 "xyz": artifact(inp.with_suffix(".xyz")), "native_state": native_state(inp),
                 "convergence_rows": rows, "constraint_gap_rad": gap,
                 "all_five_printed_criteria_passed": all_five,
                 "native_convergence_route": "all_five" if all_five else "overachieved_gradient_candidate",
                 "saddle_character_verified_by_optimization_alone": False,
                 "fixed_geometry_derivative_validation_required": True,
                 "rate_claim_allowed": False}


def read_engrad(path, xyz_bohr):
    fields = [line.strip().split() for line in Path(path).read_text().splitlines()
              if line.strip() and not line.lstrip().startswith("#")]
    if len(fields) != 30 or int(fields[0][0]) != 7:
        raise ValueError("Complete seven-atom fixed-geometry EnGrad required")
    energy = float(fields[1][0])
    gradient = np.array([float(row[0]) for row in fields[2:23]])
    rows = fields[23:]
    if tuple(int(row[0]) for row in rows) != ATOMIC_NUMBERS:
        raise ValueError("EnGrad mapped elements changed")
    coordinates = np.array([[float(v) for v in row[1:]] for row in rows])
    if (not math.isfinite(energy) or not np.isfinite(gradient).all()
            or not np.isfinite(coordinates).all()
            or not np.allclose(coordinates, xyz_bohr, atol=1e-7, rtol=0)):
        raise ValueError("EnGrad is not the exact fixed Hessian geometry")
    return energy, gradient


def aligned_mode_overlap(source_xyz, candidate_xyz, source_vector, candidate_vector, masses):
    left, _, right = np.linalg.svd((candidate_xyz-candidate_xyz.mean(0)).T @ (source_xyz-source_xyz.mean(0)))
    if np.linalg.det(left @ right) < 0:
        left[:, -1] *= -1
    rotation = left @ right  # Proper rotation only, never mirror a mapped saddle.
    ref = (source_vector * np.sqrt(masses[:, None])).ravel()
    new = ((candidate_vector @ rotation) * np.sqrt(masses[:, None])).ravel()
    return float(abs(np.dot(ref, new)) / np.linalg.norm(ref) / np.linalg.norm(new))


def frequency_receipt(inp, xyz_A, origin, policy):
    from storca.ts_angular_audit import constrained_saddle_diagnostics
    H, masses, xyz = origin["matrix_reader"].read_matrix(artifact(inp.with_suffix(".hess")))
    modes = parse_orca_normal_modes(inp.with_suffix(".hess"))
    if (modes.atom_labels != LABELS or not np.array_equal(masses, origin["masses_amu"])
            or not np.allclose(xyz, xyz_A / BOHR_A, atol=1e-9, rtol=0)
            or not np.allclose(modes.coordinates_bohr, xyz, atol=1e-9, rtol=0)
            or not np.all(np.abs(modes.frequencies_cm_1[:6]) < 0.01)):
        raise ValueError("Fresh fixed-geometry NumFreq binding differs")
    text = inp.with_suffix(".hess").read_text()
    if re.findall(r"\$frequency_scale_factor\s+([\d.]+)", text) != ["1.000000000000"]:
        # Native formatting can differ without changing the exact numerical value.
        scales = re.findall(r"\$frequency_scale_factor\s+([\d.]+)", text)
        if len(scales) != 1 or float(scales[0]) != 1.0:
            raise ValueError("Unscaled numerical Hessian required")
    energy, gradient = read_engrad(inp.with_suffix(".engrad"), xyz)
    native = native_state(inp)
    if max(abs(value - 0.758947) for value in native["S_squared"]) > policy["maximum_S_squared_change"]:
        raise ValueError("Electronic spin continuity warning")
    spectrum = constrained_saddle_diagnostics(xyz, masses, H, gradient, DIHEDRAL)
    vector = np.asarray(spectrum["unstable_cartesian_mode"])
    transfer = bond_change_mode_overlap(xyz, vector, TRANSFER)
    signed = [row["signed_stretch_component"] for row in transfer["coordinates"]]
    overlap = aligned_mode_overlap(origin["xyz_A"] / BOHR_A, xyz,
                                   origin["modes"].vectors[6], vector, masses)
    passed = (spectrum["negative_count"] == 1 and spectrum["positive_count"] == 13
              and spectrum["tangent_gradient_RMS_Eh_bohr"] <= policy["tangent_gradient_RMS_Eh_bohr"]
              and spectrum["tangent_gradient_max_Eh_bohr"] <= policy["tangent_gradient_max_Eh_bohr"]
              and transfer["status"] == "matched" and signed[0]*signed[1] < 0
              and overlap >= policy["minimum_transfer_mode_overlap"])
    return {"input": artifact(inp), "output": artifact(inp.with_suffix(".out")),
            "hessian": artifact(inp.with_suffix(".hess")), "engrad": artifact(inp.with_suffix(".engrad")),
            "energy_Eh": energy, "native_state": native, "constrained_spectrum": spectrum,
            "transfer_mode": transfer, "retained_transfer_mode_overlap": overlap,
            "pilot_point_gate_passed": bool(passed),
            "raw_full_frequencies_cm_inv": modes.frequencies_cm_1.tolist(),
            "rate_claim_allowed": False, "lifetime_claim_allowed": False}


def orthogonal_comparison(points, temperatures):
    center = points[0]["frequency"]["constrained_spectrum"]["positive_frequencies_cm_inv"]
    if len(center) != 13:
        raise ValueError("Thirteen positive tangent modes required")

    def logQ(frequencies, T):
        return float(-np.log(-np.expm1(-np.asarray(frequencies) * 1.4387768775039338 / T)).sum())

    results = []
    for point in points[1:]:
        frequencies = point["frequency"]["constrained_spectrum"]["positive_frequencies_cm_inv"]
        if len(frequencies) != 13 or any(v <= 0 for v in frequencies):
            raise ValueError("Orthogonal mode inventory changed")
        results.append({"offset_deg": point["offset_deg"],
                        "ZPE_shift_J_mol": float((sum(frequencies)-sum(center)) * CM_J_MOL/2),
                        "maximum_abs_log_partition_shift": max(abs(logQ(frequencies, T)-logQ(center, T)) for T in temperatures),
                        "energy_shift_J_mol": (point["frequency"]["energy_Eh"]-points[0]["frequency"]["energy_Eh"])*EH_J_MOL,
                        "scope": "local constrained spectra coupling diagnostic; not a rotor contribution or rate"})
    return results


def run(source):
    from storca import ts_angular_audit as angular
    origin = load_source(source)  # Refuse before creating a campaign or dispatching.
    folder = Path(tempfile.mkdtemp(prefix="oh-methane-angular-pilot-20261011-", dir=ROOT / "runs"))
    print("WORKDIR " + str(folder), flush=True)
    attempts = {"saddles": 0, "NumFreq": 0}
    points, failure, contract = [], None, None
    started = time.monotonic()
    try:
        contract = build_contract(origin, folder)
        save(folder / "contract.json", contract)
        geometry = angular.angular_chart_diagnostics(origin["xyz_A"] / BOHR_A, origin["masses_amu"],
                      DIHEDRAL, soft_mode=origin["modes"].vectors[7], unstable_mode=origin["modes"].vectors[6])
        save(folder / "coordinate-diagnostic.json", geometry)
        if geometry["soft_mode_overlap"] < 0.75 or geometry["unstable_mode_overlap"] > 0.05:
            raise ValueError("Angular chart does not isolate the intended soft direction sufficiently for pilot")
        for index, offset in enumerate(PILOT_OFFSETS_DEG):
            state.reopen_inventory(contract["source_inventory"])
            seed = angular.rotate_dihedral_seed(origin["xyz_A"] / BOHR_A, DIHEDRAL, math.radians(offset)) * BOHR_A
            point_folder = folder / f"point-{index:02d}"
            point_folder.mkdir()
            row = {"offset_deg": offset, "seed_A": seed.tolist(), "status": "started"}
            points.append(row)
            if index == 0:
                # Exact input/output/lease/XYZ are bound above. Reuse the job,
                # not its old failed audit decision; new derivatives remain fresh.
                opt = REUSE_CENTER
                lease = read_json(opt.with_suffix(".process.json"))
                if (lease["returncode"] != 0 or lease["input_sha256"] != artifact(opt)["sha256"]
                        or opt.with_suffix(".out").read_text().count("ORCA TERMINATED NORMALLY") != 1):
                    raise ValueError("Reused center lacks owned unique native completion")
                row["optimization_reused"] = True
            else:
                opt = execute(point_folder / "optimization", "opt", input_text(seed, "opt", origin["hessian"]), attempts)
            xyz, receipt = optimization_receipt(opt, seed, contract["policy"])
            row["optimization"] = receipt
            freq = execute(point_folder / "numerical-Hessian", "freq", input_text(xyz, "freq"), attempts)
            row["frequency"] = frequency_receipt(freq, xyz, origin, contract["policy"])
            row["status"] = "passed_pilot_point" if row["frequency"]["pilot_point_gate_passed"] else "unresolved"
            save(point_folder / "point-result.json", row)
            if row["status"] == "unresolved":
                raise ValueError("Constraint-projected saddle/transfer/tangent-gradient pilot gate unresolved")
        comparisons = orthogonal_comparison(points, origin["temperature_grid_K"])
        save(folder / "orthogonal-comparison.json", comparisons)
        if any(abs(row["ZPE_shift_J_mol"]) > contract["policy"]["maximum_orthogonal_ZPE_shift_J_mol"]
               or row["maximum_abs_log_partition_shift"] > contract["policy"]["maximum_orthogonal_partition_log_shift"]
               for row in comparisons):
            raise ValueError("Pilot orthogonal contribution varies materially; separable one-angle correction unresolved")
    except BaseException as error:
        failure = f"{type(error).__name__}: {error}"
        print("PILOT UNRESOLVED: " + failure, flush=True)
    unchanged = True
    if contract is not None:
        try:
            state.reopen_inventory(contract["source_inventory"])
        except (ValueError, OSError) as error:
            unchanged = False
            failure = failure or str(error)
    result = {"schema_version": 1, "kind": "bounded_OH_CH4_TS_angular_saddle_pilot",
              "status": "unresolved" if failure else "pilot_passed_profile_pending",
              "failure_reason": failure, "points": points, "attempts": attempts,
              "whole_increment_attempts": {key: value + (contract["prior_attempts"][key] if contract else 0)
                                           for key, value in attempts.items()},
              "elapsed_seconds": time.monotonic()-started, "old_sources_unchanged": unchanged,
              "contract": artifact(folder / "contract.json") if contract is not None else None,
              "full_periodic_domain_covered": False, "angular_partition_admitted": False,
              "mirror_basin_count_resolved": False, "verification_status": "diagnostic_not_admitted",
              "rate_claim_allowed": False, "lifetime_claim_allowed": False,
              "production_defaults_changed": False, "network_reinjection_allowed": False,
              "next_step": ("Freeze unresolved angular assumptions; no automatic method or multidimensional search."
                            if failure else "Pilot supports a separate bounded periodic-profile continuation, not a rate correction.")}
    save(folder / "result.json", result)
    print(result["status"] + ": " + str(folder / "result.json"), flush=True)
    if failure:
        raise RuntimeError(failure)
    return folder


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Completed electronic-state replay directory")
    run(parser.parse_args().source)
