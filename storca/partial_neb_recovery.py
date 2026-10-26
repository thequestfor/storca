"""Read-only admission of a converged NEB stage followed by failed OptTS.

This is deliberately not a restart runner or a transition-state validator.
Retained side XYZ files are only witnesses: coordinates must also agree, in
the original mapped order, with the terminal-owned native output and band.
Every charge, failed parent and scientific validation gate remains intact.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np


POLICY_ID = "owned-converged-neb-failed-optts-candidate-v1"
_NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"
_COORDINATE_TOLERANCE = 1.1e-6  # Native output/MEP print six decimal places.
_ENERGY_TOLERANCE = 1e-8


def _hash(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Missing regular source: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path):
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("Duplicate source JSON key")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def _equivalent(left, right):
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(right, sort_keys=True, allow_nan=False)


def _float(value):
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError("Nonfinite native coordinate or energy")
    return result


def _frames(path, *, require_energy=True):
    lines = Path(path).read_text().splitlines()
    frames, index = [], 0
    while index < len(lines):
        if not lines[index].strip() or lines[index].strip() == ">":
            index += 1
            continue
        count = int(lines[index].strip())
        if count <= 0 or index + count + 1 >= len(lines):
            raise ValueError("Incomplete native XYZ frame")
        rows = [row.split() for row in lines[index + 2:index + count + 2]]
        if any(len(row) != 4 or not re.fullmatch(r"[A-Z][a-z]?", row[0]) for row in rows):
            raise ValueError("Malformed native XYZ atom row")
        energy = re.search(r"\bE\s*(?:=\s*)?(" + _NUMBER + r")\s*$", lines[index + 1])
        if not energy and require_energy:
            raise ValueError("Native XYZ frame lacks an explicit finite energy")
        frames.append({"elements": [row[0] for row in rows],
                       "coordinates": [[_float(value) for value in row[1:]] for row in rows],
                       "energy_hartree": _float(energy[1]) if energy else None, "comment": lines[index + 1]})
        index += count + 2
    if not frames:
        raise ValueError("Empty native XYZ trajectory")
    return frames


def _printed_coordinates(text, title, count):
    matches = list(re.finditer(re.escape(title), text))
    if len(matches) != 1:
        raise ValueError(f"Missing or ambiguous native coordinate block: {title}")
    rows = []
    for line in text[matches[0].end():].splitlines():
        fields = line.split()
        if len(fields) == 4 and re.fullmatch(r"[A-Z][a-z]?", fields[0]):
            rows.append(fields)
            if len(rows) == count:
                break
        elif rows:
            raise ValueError("Incomplete printed native coordinate block")
        elif line.strip() and not set(line.strip()) <= {"-"}:
            raise ValueError("Unexpected text before native coordinates")
    if len(rows) != count:
        raise ValueError("Truncated printed native coordinates")
    return {"elements": [row[0] for row in rows],
            "coordinates": [[_float(value) for value in row[1:]] for row in rows]}


def _same_coordinates(left, right, elements):
    if left["elements"] != elements or right["elements"] != elements:
        raise ValueError("Native atom identity/order differs from the mapped route")
    if not np.allclose(left["coordinates"], right["coordinates"], rtol=0,
                       atol=_COORDINATE_TOLERANCE):
        raise ValueError("Retained native coordinates lack exact mapped-order correspondence")


def inspect_partial_neb_recovery(input_path, *, control, context_binding,
                                 expected_context, source_bindings,
                                 selection="last_evaluated_opt"):
    """Admit only a candidate from an owned failed compound NEB-TS attempt.

    ``control`` is the existing NativeExecutionControl, not a replacement
    budget. ``context_binding`` pins the original surface-resolution scientific
    context, which must equal the caller's independently frozen
    ``expected_context``. ``source_bindings`` must pin the job contract and
    named CI/OptTS trajectory, which older terminal receipts did not hash.
    Their contents are independently matched against terminal-owned output.

    Selections are ``converged_ci``, ``last_evaluated_opt`` (a successfully
    evaluated but unconverged geometry), or ``failing_geometry_diagnostic``
    (never eligible for TS recovery). No files, ledger events or jobs are made.
    """
    result = {"schema_version": 1, "kind": "partial_neb_recovery_admission",
              "policy_id": POLICY_ID, "status": "partial_neb_recovery_refused", "valid": False,
              "neb_converged": False, "transition_state_verified": False,
              "route_claim_allowed": False, "rate_claim_allowed": False,
              "lifetime_claim_allowed": False, "renews_compute_budget": False,
              "native_jobs_launched": 0,
              "claims_prohibited": ["verified_transition_state", "verified_route", "barrier", "rate", "lifetime"]}
    before = None
    try:
        if selection not in {"converged_ci", "last_evaluated_opt", "failing_geometry_diagnostic"}:
            raise ValueError("Unsupported partial-stage candidate selection")
        inp = Path(input_path).absolute()
        _hash(inp)
        inp = inp.resolve(strict=True)
        references = {}
        def pin(binding):
            raw = Path(binding["path"])
            digest = _hash(raw)
            path = raw.resolve(strict=True)
            if digest != binding["sha256"] or (str(path) in references and references[str(path)]["sha256"] != digest):
                raise ValueError(f"Source hash differs: {path}")
            references[str(path)] = {"path": str(path), "sha256": digest}
            return path
        context_path = pin(context_binding)
        for reference in source_bindings:
            pin(reference)
        context = _json(context_path)
        if not _equivalent(context, expected_context):
            raise ValueError("Frozen route/mapping/environment/method context differs")
        from .route_verify import RouteSpec
        from .route_geometry import _build_endpoint_side, resolve_route_atom_mapping
        route = RouteSpec(**context["route"])
        mapping = resolve_route_atom_mapping(route)
        if (not mapping.get("valid") or not _equivalent(mapping, context["mapping"])
                or context.get("kind") != "storca_surface_resolution_scientific_context"
                or context["environment_applicability"].get("status") != "applicable"):
            raise ValueError("Mapped route or environmental applicability is unresolved")
        electronic = context["electronic_structure"]
        if electronic["charge"] != route.charge or electronic["multiplicity"] != route.multiplicity:
            raise ValueError("Declared route and electronic state differ")
        elements = _build_endpoint_side(tuple(route.reactant_smiles)).elements
        before = control.snapshot()
        if before["pending"]:
            raise ValueError("Native ownership has a pending reservation")
        matches = [row for row in before["reservations"].values() if Path(row["input"]).resolve() == inp]
        if len(matches) != 1:
            raise ValueError("Exactly one owned native attempt required")
        reservation = matches[0]
        terminal = before["terminals"].get(reservation["attempt_id"])
        if (terminal is None or terminal.get("outcome") != "failed"
                or reservation["category"] != "optimization_path"
                or reservation["ncores"] != electronic["ncores"]):
            raise ValueError("Failed terminal optimization/path ownership required")
        owned = set()
        for reference in terminal["sources"]:
            owned.add(str(pin(reference)))
        def required(path, *, terminal_owned=False):
            path = Path(path).resolve(strict=True)
            if str(path) not in references or (terminal_owned and str(path) not in owned):
                raise ValueError(f"Required source lacks {'terminal ownership' if terminal_owned else 'a pinned hash'}: {path}")
            return path
        out = required(inp.with_suffix(".out"), terminal_owned=True)
        lease_path = required(inp.with_suffix(".process.json"), terminal_owned=True)
        contract_path = required(inp.with_suffix(".contract.json"))
        band_path = required(inp.parent / f"{inp.stem}_MEP.allxyz", terminal_owned=True)
        ci_path = required(inp.parent / f"{inp.stem}_NEB-CI_converged.xyz")
        trajectory_path = required(inp.parent / f"{inp.stem}_trj.xyz")
        required(inp, terminal_owned=True)
        text = inp.read_text()
        keywords = [token.casefold() for line in re.findall(r"(?m)^\s*!([^\n]*)", text) for token in line.split()]
        if (keywords.count("neb-ts") != 1
                or sorted(token for token in keywords if token not in {"neb-ts", "verytightopt"})
                != sorted(token.casefold() for token in electronic["method_keywords"])):
            raise ValueError("Compound NEB method differs from the frozen electronic method")
        xyz = re.findall(r"(?im)^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)\s*$", text)
        images = re.findall(r"(?im)^\s*NImages\s+(\d+)\s*$", text)
        cores = re.findall(r"(?im)^\s*nprocs\s+(\d+)\s*$", text)
        if (len(xyz) != 1 or (int(xyz[0][0]), int(xyz[0][1])) != (route.charge, route.multiplicity)
                or len(images) != 1 or cores != [str(electronic["ncores"]) ]):
            raise ValueError("Native state/image count/cores differ")
        dependencies = set(re.findall(r'(?:Product|Hess_Filename|Restart_ALLXYZFile)\s+"([^"]+)"', text, re.I))
        dependencies.add(xyz[0][2].strip('"'))
        digest = hashlib.sha256(inp.read_bytes())
        for name in sorted(dependencies):
            path = inp.parent / name
            _hash(path)
            digest.update(name.encode())
            digest.update(path.read_bytes())
            pin({"path": str(path), "sha256": _hash(path)})
        contract = _json(contract_path)
        if (contract.get("sha256") != digest.hexdigest()
                or reservation["input_sha256"] != digest.hexdigest()
                or reservation["native_input_file_sha256"] != _hash(inp)):
            raise ValueError("Owned native input/dependency digest differs")
        control.validate_reuse(inp, digest.hexdigest())
        lease = _json(lease_path)
        if (lease.get("kind") != "storca_external_process_lease" or lease.get("engine") != "orca"
                or lease.get("status") != "terminated" or type(lease.get("returncode")) is not int
                or lease.get("input_sha256") != _hash(inp) or Path(lease.get("input", "")).resolve() != inp
                or Path(lease.get("output", "")).resolve() != out or lease.get("ncores") != electronic["ncores"]
                or not all(math.isfinite(lease[key]) for key in ("started_at_epoch", "ended_at_epoch"))
                or lease["ended_at_epoch"] < lease["started_at_epoch"]):
            raise ValueError("Native process lacks stopped, matching terminal evidence")
        output = out.read_text(errors="strict")
        converged = "THE NEB OPTIMIZATION HAS CONVERGED"
        opt_start = "Next a TS optimization is performed on the CI from the NEB run."
        if (output.count(converged) != 1 or output.count(opt_start) != 1
                or output.index(converged) >= output.index(opt_start)
                or "ORCA TERMINATED NORMALLY" in output or "THE OPTIMIZATION HAS CONVERGED" in output
                or "ORCA finished by error termination" not in output):
            raise ValueError("Only native-converged NEB followed by failed OptTS is eligible")
        summary, optimization = output.split(opt_start)
        ci_index_match = re.findall(r"Climbing image\s*\.{2,}\s*(\d+)", summary)
        ci_energy_match = re.findall(r"Energy\s*\.{2,}\s*(" + _NUMBER + r")\s+Eh", summary)
        if len(ci_index_match) != 1 or len(ci_energy_match) != 1:
            raise ValueError("Native CI identity/energy ambiguous")
        ci_index = int(ci_index_match[0])
        band = _frames(band_path)
        if len(band) != int(images[0]) + 2 or not 0 < ci_index < len(band) - 1:
            raise ValueError("Native band size or interior CI index differs")
        for index, frame in enumerate(band):
            if frame["elements"] != elements or not re.search(r"\bNEB Path Image\s+" + str(index + 1) + r"\b", frame["comment"]):
                raise ValueError("Native band image ordering differs")
        ci = _frames(ci_path)
        if len(ci) != 1:
            raise ValueError("Native CI is not a single retained geometry")
        printed_ci = _printed_coordinates(summary, "SADDLE POINT (ANGSTROEM)", len(elements))
        _same_coordinates(ci[0], printed_ci, elements)
        _same_coordinates(ci[0], band[ci_index], elements)
        if abs(ci[0]["energy_hartree"] - _float(ci_energy_match[0])) > _ENERGY_TOLERANCE:
            raise ValueError("Native CI energy differs from terminal-owned output")
        # MEP comments may retain prior-cycle energies. Coordinate agreement
        # establishes image correspondence only, never an energetic threshold.
        cycles = list(re.finditer(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)", optimization))
        if len(cycles) < 2 or [int(match[1]) for match in cycles] != list(range(1, len(cycles) + 1)):
            raise ValueError("OptTS cycle sequence missing or reordered")
        trajectory = _frames(trajectory_path)
        successful = []
        failing = None
        for index, match in enumerate(cycles):
            block = optimization[match.end():cycles[index + 1].start() if index + 1 < len(cycles) else len(optimization)]
            geometry = printed_ci if index == 0 else _printed_coordinates(block, "CARTESIAN COORDINATES (ANGSTROEM)", len(elements))
            energies = re.findall(r"FINAL SINGLE POINT ENERGY\s+(" + _NUMBER + r")", block)
            gradient = re.findall(r"(?m)^\s*MAX gradient\s+(" + _NUMBER + r")\s+(" + _NUMBER + r")\s+(YES|NO)\s*$", block)
            completed = bool("SCF CONVERGED AFTER" in block and len(energies) == 1
                             and "ORCA SCF GRADIENT CALCULATION" in block and len(gradient) == 1)
            if completed:
                if failing is not None:
                    raise ValueError("Evaluated OptTS cycle follows an incomplete cycle")
                successful.append({**geometry, "optimization_cycle": index + 1,
                                   "energy_hartree": _float(energies[0]),
                                   "maximum_gradient_hartree_per_bohr": _float(gradient[0][0])})
            elif index == len(cycles) - 1 and not energies and not gradient:
                failing = {**geometry, "optimization_cycle": index + 1, "energy_hartree": None}
            else:
                raise ValueError("Incomplete or inconsistent OptTS evaluation evidence")
        if not successful or failing is None or len(trajectory) != len(successful):
            raise ValueError("OptTS trajectory does not match completed native evaluations")
        for frame, evaluated in zip(trajectory, successful):
            _same_coordinates(frame, evaluated, elements)
            if abs(frame["energy_hartree"] - evaluated["energy_hartree"]) > _ENERGY_TOLERANCE:
                raise ValueError("OptTS trajectory energy differs from completed native evaluation")
        failing_path = required(inp.with_suffix(".xyz"), terminal_owned=True)
        failed_frames = _frames(failing_path, require_energy=False)
        if len(failed_frames) != 1:
            raise ValueError("Failed native geometry is ambiguous")
        _same_coordinates(failed_frames[0], failing, elements)
        geometries = {"converged_ci": {**ci[0], "origin": "native_converged_climbing_image"},
                      "last_evaluated_opt": {**successful[-1], "origin": "last_completed_opt_evaluation"},
                      "failing_geometry_diagnostic": {**failing, "origin": "failed_scf_geometry_not_evaluated"}}
        for reference in references.values():
            if _hash(reference["path"]) != reference["sha256"]:
                raise ValueError("Native source changed during admission")
        if control.snapshot() != before:
            raise ValueError("Native journal changed during read-only admission")
        result.update(valid=True, status="partial_neb_candidate_admitted", neb_converged=True,
                      native_reservation_id=reservation["attempt_id"], failed_parent_remains_charged=True,
                      selection=selection, candidate_geometry=geometries[selection],
                      eligible_for_ts_recovery=selection != "failing_geometry_diagnostic",
                      converged_ci=geometries["converged_ci"],
                      last_evaluated_opt=geometries["last_evaluated_opt"],
                      failing_geometry_diagnostic=geometries["failing_geometry_diagnostic"],
                      completed_opt_evaluations=len(successful), native_ci_index=ci_index,
                      context=context_binding, mapping=mapping,
                      sources=[references[key] for key in sorted(references)],
                      band_energy_interpretation="coordinate_correspondence_only; native comments may contain stale energies",
                      requirements=["converged TS refinement and electronic-state checks", "independent Hessian and mapped reaction mode",
                                    "both IRC directions", "validated actual reached minima"])
    except (OSError, ValueError, TypeError, KeyError, IndexError, RuntimeError) as error:
        result["failure_reason"] = f"{type(error).__name__}: {error}"
    return result
