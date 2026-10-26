"""Read-only admission of completed IRC evidence into a new owned workspace.

Reusing a trace is not reusing an endpoint, stationary-point, or rate claim.
Endpoint chemistry and local minima must be checked under the caller's current
policy.  In particular, a historical RMSD rejection is never relabeled here.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil

import numpy as np

from .decomposition_visuals import xyz_frames
from .route_geometry import read_xyz


def _digest(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _equivalent(left, right) -> bool:
    return json.dumps(left, sort_keys=True) == json.dumps(right, sort_keys=True)


def _references(request: dict) -> dict[str, dict]:
    references = {}
    for reference in [request["source_path_result"], request["source_execution_contract"],
                      request["source_run_result"], *request["source_artifacts"]]:
        path = Path(reference["path"]).resolve(strict=True)
        if str(path) in references and references[str(path)]["sha256"] != reference["sha256"]:
            raise ValueError("Conflicting completed IRC source hashes")
        if _digest(path) != reference["sha256"]:
            raise ValueError(f"Completed IRC source changed: {path}")
        references[str(path)] = {"path": str(path), "sha256": reference["sha256"]}
    return references


def _pinned(value, references):
    path = Path(value).resolve(strict=True)
    if str(path) not in references or _digest(path) != references[str(path)]["sha256"]:
        raise ValueError(f"Completed IRC source missing or changed hash: {path}")
    return path


def _terminal_job(inp, references, *, require_success=True):
    """Validate owned, hash-bound jobs without dispatching or renewing leases."""
    from . import reaction_path as rp

    inp = _pinned(inp, references)
    out = _pinned(inp.with_suffix(".out"), references)
    contract = json.loads(_pinned(inp.with_suffix(".contract.json"), references).read_text())
    digest, dependencies = rp._orca_job_input_evidence(inp)
    if (contract.get("kind") != "orca_job_input_contract"
            or Path(contract.get("input", "")).resolve() != inp
            or contract.get("sha256") != digest.hexdigest()
            or contract.get("referenced_files") != dependencies):
        raise ValueError("Completed IRC job input contract or dependency ownership differs")
    for dependency in dependencies:
        _pinned(dependency, references)
    lease = json.loads(_pinned(inp.with_suffix(".process.json"), references).read_text())
    if (lease.get("kind") != "storca_external_process_lease" or lease.get("engine") != "orca"
            or lease.get("status") != "terminated" or type(lease.get("returncode")) is not int
            or lease.get("input_sha256") != _digest(inp)
            or Path(lease.get("input", "")).resolve() != inp
            or Path(lease.get("output", "")).resolve() != out
            or lease.get("ncores") != rp._orca_input_ncores(inp)
            or not all(math.isfinite(float(lease[key])) for key in ("started_at_epoch", "ended_at_epoch"))
            or float(lease["ended_at_epoch"]) < float(lease["started_at_epoch"])):
        raise ValueError("Completed IRC job lease is not terminal with exact ownership")
    for field, sign in (("pid", 1), ("process_group_id", -1)):
        identifier = lease.get(field)
        if type(identifier) is not int or identifier <= 0:
            raise ValueError("Completed IRC lease lacks positive process identifiers")
        try:
            os.kill(sign * identifier, 0)
        except ProcessLookupError:
            pass
        else:
            raise ValueError("Completed IRC source process may still be alive")
    if require_success and (lease["returncode"] != 0 or not rp._normal_orca_output(out)):
        raise ValueError("Completed IRC job did not terminate normally")
    return dependencies


def _same_frame(left, right):
    return left[0] == right[0] and np.allclose(left[1], right[1], atol=1e-6, rtol=0)


def _reopen_assembly(inp: Path, retained_trajectory: Path, references: dict) -> dict:
    """Reconstruct the existing assembly in memory; never call its writer."""
    from . import reaction_path as rp

    root = inp.parent
    original = rp._orca_irc_completion(inp.with_suffix(".out"))
    if not original["normal_termination"] or set(original["directions"]) != {"forward", "backward"}:
        raise ValueError("Original IRC directional evidence is incomplete")
    manifest = _pinned(retained_trajectory.parent / "assembly.json", references)
    assembly = json.loads(manifest.read_text())
    if (assembly.get("kind") != "storca_irc_branch_assembly" or assembly.get("valid") is not True
            or Path(assembly.get("trajectory", "")).resolve() != retained_trajectory
            or assembly.get("trajectory_sha256") != _digest(retained_trajectory)):
        raise ValueError("Completed IRC assembly manifest differs from its trajectory")
    sources = {}

    def retain(value):
        path = _pinned(value, references)
        sources[str(path)] = _digest(path)
        return path

    for path in (inp, inp.with_suffix(".out"), root / f"{inp.stem}_IRC_Full_trj.xyz"):
        retain(path)
    base = xyz_frames(root / f"{inp.stem}_IRC_Full_trj.xyz")
    if not base:
        raise ValueError("Completed IRC original trajectory is empty")
    tails, extras, diagnostics = {}, {}, {}
    for branch, suffix in (("forward", "F"), ("backward", "B")):
        endpoint = retain(root / f"{inp.stem}_IRC_{suffix}.xyz")
        branch_trajectory = retain(root / f"{inp.stem}_IRC_{suffix}_trj.xyz")
        tail = read_xyz(endpoint)
        branch_frames = xyz_frames(branch_trajectory)
        if not branch_frames or not _same_frame(tail, branch_frames[-1]):
            raise ValueError("Original IRC endpoint and trajectory disagree")
        tails[branch], extras[branch] = tail, []
        direction = original["directions"][branch]
        diagnostics[branch] = direction
        if direction["valid"]:
            continue
        if not direction["iteration_limit_reached"]:
            raise ValueError("IRC branch is not eligible for downhill continuation")
        directory = root / "irc-branch-continuation" / branch
        contract = json.loads(retain(directory / "continuation-contract.json").read_text())
        expected = {"policy_id": "irc-downhill-single-extension-v1", "branch": branch,
                    "source_input_sha256": sources[str(inp)],
                    "source_output_sha256": sources[str(inp.with_suffix(".out"))],
                    "source_endpoint_sha256": sources[str(endpoint)],
                    "source_trajectory_sha256": sources[str(branch_trajectory)]}
        if any(contract.get(key) != value for key, value in expected.items()):
            raise ValueError("IRC continuation source contract changed")
        downhill = retain(directory / "irc-downhill.inp")
        seed = retain(directory / "branch-endpoint.xyz")
        if (downhill.read_text() != rp._downhill_irc_input_text(inp.read_text(), seed.name)
                or seed.read_bytes() != endpoint.read_bytes()):
            raise ValueError("IRC continuation method or seed changed")
        result = json.loads(retain(directory / "continuation-result.json").read_text())
        if any(result.get(key) != value for key, value in contract.items()):
            raise ValueError("IRC continuation execution source contract changed")
        execution = result["execution"]
        execution_sources = rp._irc_downhill_execution_sources(downhill, execution)
        recorded = result.get("execution_source_hashes") or {}
        selected = Path(execution["input"]).resolve()
        if selected != downhill and (not recorded or any(recorded.get(path) != digest
                                                         for path, digest in execution_sources.items())):
            raise ValueError("IRC continuation retry lacks complete source hashes")
        for path, digest in recorded.items():
            if Path(path).resolve().parent != directory or _digest(retain(path)) != digest:
                raise ValueError("IRC continuation execution source hash changed")
        for path, digest in execution_sources.items():
            if _digest(retain(path)) != digest:
                raise ValueError("IRC continuation execution source changed during admission")
            if Path(path).suffix == ".inp":
                _terminal_job(path, references, require_success=Path(path) == selected)
        diagnostic = rp._orca_irc_completion(Path(execution["output"]), required_directions=("downhill",))
        if not diagnostic["valid"]:
            raise ValueError("Downhill IRC branch did not converge")
        diagnostics[branch] = diagnostic
        continued = retain(directory / f"{selected.stem}_IRC_D_trj.xyz")
        final = retain(directory / f"{selected.stem}_IRC_D.xyz")
        frames = xyz_frames(continued)
        if not frames or not _same_frame(frames[0], tail) or not _same_frame(frames[-1], read_xyz(final)):
            raise ValueError("Downhill IRC has a discontinuous junction or endpoint")
        extras[branch] = frames[1:]
    if not _same_frame(base[0], tails["backward"]) or not _same_frame(base[-1], tails["forward"]):
        raise ValueError("Original full IRC does not connect its directional tails")
    frames = list(reversed(extras["backward"])) + base + extras["forward"]
    if any(frame[0] != base[0][0] or not np.isfinite(frame[1]).all() for frame in frames):
        raise ValueError("Completed IRC changes atom order or contains nonfinite coordinates")
    expected_text = "".join(
        f"{len(labels)}\n{comment}\n" + "".join(
            f"{label} {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f}\n"
            for label, xyz in zip(labels, coordinates)) for labels, coordinates, comment in frames)
    identity = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    if (retained_trajectory.read_text() != expected_text or assembly.get("source_hashes") != sources
            or assembly.get("frame_count") != len(frames)
            or not _equivalent(assembly.get("direction_completion"), diagnostics)
            or retained_trajectory.parent != root / "irc-branch-continuation" / f"assembled-{identity[:16]}"):
        raise ValueError("Completed IRC assembly cannot be reconstructed from its owned sources")
    return {"kind": "continued_bidirectional_irc_completion", "valid": True,
            "original_completion": original, "assembly": assembly,
            "revalidated_read_only": True}


def inspect_completed_irc_reuse(route, request: dict, *, historical_admission: dict,
                                ncores: int, method_keywords: list[str] | None,
                                surface_resolution_policy=None) -> dict:
    """Reopen source evidence, not previous endpoint/path/rate conclusions."""
    from . import reaction_path as rp

    result = {"schema_version": 1, "kind": "completed_irc_reuse_admission", "valid": False,
              "status": "completed_irc_reuse_refused", "renews_compute_budget": False,
              "claims_permitted": [], "claims_prohibited": ["verified_path", "validated_minima", "rate", "lifetime"]}
    try:
        if (request.get("schema_version") != 1 or request.get("kind") != "completed_irc_reuse_request"
                or historical_admission.get("valid") is not True):
            raise ValueError("Explicit source-bound completed IRC reuse request required")
        container = request["ts_frequency_request"]
        if container.get("kind") != "completed_ts_frequency_reuse_request" or container.get("schema_version") != 1:
            raise ValueError("Completed IRC source container request differs")
        references = _references(container)
        source_path = _pinned(container["source_path_result"]["path"], references)
        contract_path = _pinned(container["source_execution_contract"]["path"], references)
        run_path = _pinned(container["source_run_result"]["path"], references)
        source, contract, run = (json.loads(path.read_text()) for path in (source_path, contract_path, run_path))
        original = json.loads(Path(historical_admission["source_execution_contract"]["path"]).read_text())
        if (source.get("kind") != "generic_orca_route_verification"
                or not _equivalent(source.get("route"), route.as_dict())
                or not source_path.is_relative_to(contract_path.parent) or run_path.parent != contract_path.parent
                or contract.get("execute") is not True
                or contract.get("condition_contract") != original.get("condition_contract")
                or contract.get("method_profile") != original.get("method_profile")
                or (contract.get("budget") or {}).get("maximum_cores") != ncores
                or run.get("status") not in {"workflow_returned", "failed"}
                or run.get("original_sources_unchanged") is not True or run.get("manifest_unchanged") is not True):
            raise ValueError("Completed IRC owner differs from the terminal scientific execution contract")
        prior = source["historical_neb_seed_admission"]
        if (not prior.get("valid") or prior.get("seed_index") != historical_admission["seed_index"]
                or not _equivalent(prior.get("source_checkpoint"), historical_admission["source_checkpoint"])):
            raise ValueError("Completed IRC was not generated from this admitted historical image")
        context = rp._surface_context(source_path.parent)
        declared = source["surface_resolution_scientific_context"]
        method = list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"])
        if (context is None or not _equivalent(context.get("route"), route.as_dict())
                or not _equivalent(context.get("mapping"), historical_admission["mapping"])
                or _pinned(context["path"], references) != _pinned(declared["path"], references)
                or context["scientific_context_sha256"] != declared["scientific_context_sha256"]
                or context["electronic_structure"].get("method_keywords") != method
                or context["electronic_structure"].get("ncores") != ncores
                or context["electronic_structure"].get("charge") != route.charge
                or context["electronic_structure"].get("multiplicity") != route.multiplicity):
            raise ValueError("Completed IRC scientific context, mapping, method or state differs")
        ts_request = container
        if prior.get("retained_ts_frequency") is not None:
            source_request = contract["source_request"]
            source_request_path = _pinned(source_request["path"], references)
            if _digest(source_request_path) != source_request["sha256"]:
                raise ValueError("Completed IRC original TS reuse request changed")
            ts_request = json.loads(source_request_path.read_text())["retained_ts_frequency"]
            if ts_request["source_path_result"] != prior["retained_ts_frequency"]["source_path_result"]:
                raise ValueError("Completed IRC original TS owner differs from recorded admission")
        ts_admission = rp.inspect_historical_ts_frequency_reuse(
            route, ts_request, historical_admission=historical_admission, ncores=ncores,
            method_keywords=method_keywords, surface_resolution_policy=surface_resolution_policy)
        if not ts_admission["valid"]:
            raise ValueError(f"Completed IRC TS source refused: {ts_admission.get('failure_reason')}")
        for reference in ts_admission["source_artifacts"]:
            key = reference["path"]
            if key in references and references[key] != reference:
                raise ValueError("Completed IRC and TS source references conflict")
            references[key] = reference
        trial = source["historical_neb_seed_trial"]
        ts = ts_admission["transition_state"]
        stored_ts = trial["transition_state"]
        for field in ("optimized_xyz", "frequency_output"):
            if Path(stored_ts[field]).resolve() != Path(ts[field]).resolve():
                raise ValueError("Completed IRC references a different retained TS")
        for field in ("optimization", "frequency_execution"):
            if any(Path(stored_ts[field][key]).resolve() != Path(ts[field][key]).resolve()
                   for key in ("input", "output")):
                raise ValueError("Completed IRC TS execution ownership differs")
        seed = _pinned(stored_ts["seed_xyz"], references)
        labels, coordinates, _ = read_xyz(seed)
        if (labels != historical_admission["seed_elements"]
                or not np.allclose(coordinates, historical_admission["seed_coordinates"], atol=1e-7, rtol=0)):
            raise ValueError("Completed IRC historical seed differs")
        irc = trial["irc"]
        inp = _pinned(irc["execution"]["input"], references)
        output = _pinned(irc["execution"]["output"], references)
        trajectory = _pinned(irc["trajectory"], references)
        if output != inp.with_suffix(".out") or not inp.parent.is_relative_to(source_path.parent):
            raise ValueError("Completed IRC execution is outside its owned workspace")
        text = inp.read_text()
        tokens = " ".join(re.findall(r"^\s*!\s*(.*)$", text, re.M)).lower().split()
        if ("irc" not in tokens
                or sorted(token for token in tokens if token not in {"irc", "tightscf"})
                != sorted(token.lower() for token in method if token.lower() != "tightscf")
                or re.findall(r"^\s*\*\s+xyz(?:file)?\s+(-?\d+)\s+(\d+)", text, re.M | re.I)
                != [(str(route.charge), str(route.multiplicity))] or rp._orca_input_ncores(inp) != ncores):
            raise ValueError("Completed IRC input method, state or resource declaration differs")
        dependencies = _terminal_job(inp, references)
        geometry = re.search(r"\*\s+xyzfile\s+[-+]?\d+\s+\d+\s+(\S+)", text, re.I)
        hessian = re.search(r'Hess_Filename\s+"([^"]+)"', text, re.I)
        if geometry is None or hessian is None:
            raise ValueError("Completed IRC requires explicit retained geometry and Hessian")
        if (_digest(inp.parent / geometry[1].strip('"')) != _digest(Path(ts["optimized_xyz"]))
                or _digest(inp.parent / hessian[1]) != _digest(Path(ts["frequency_execution"]["input"]).with_suffix(".hess"))):
            raise ValueError("Completed IRC geometry or Hessian differs from retained TS")
        completion = rp._orca_irc_completion(output)
        if completion["valid"]:
            native = _pinned(rp._locate_irc_trajectory(inp.parent, inp.stem), references)
            if native.read_bytes() != trajectory.read_bytes():
                raise ValueError("Completed IRC trace differs from its native trajectory")
        else:
            completion = _reopen_assembly(inp, trajectory, references)
        frames = xyz_frames(trajectory)
        ts_labels, ts_xyz, _ = read_xyz(Path(ts["optimized_xyz"]))
        if (not frames or any(frame[0] != ts_labels or not np.isfinite(frame[1]).all() for frame in frames)
                or not any(frame[0] == ts_labels and np.allclose(frame[1], ts_xyz, atol=1e-5, rtol=0)
                           for frame in frames)):
            raise ValueError("Completed IRC lacks its admitted TS or valid ordered trajectory")
        for reference in references.values():
            _pinned(reference["path"], references)
        result.update(valid=True, status="completed_irc_stages_reusable", transition_state=ts,
                      ts_frequency_admission=ts_admission, completion=completion,
                      execution={**irc["execution"], "resumed": True}, trajectory=str(trajectory),
                      dependencies=dependencies, source_artifacts=list(references.values()),
                      source_path_result=container["source_path_result"], source_workspace=str(contract_path.parent),
                      source_policy_id=source.get("surface_resolution_policy", {}).get("policy_id"),
                      endpoint_claims_reused=False, reused_stages=["ts_optimization", "ts_frequency", "irc"],
                      claims_permitted=["completed_continuous_irc_tracing"], irc_reused=True)
    except (KeyError, IndexError, TypeError, ValueError, OSError, RuntimeError) as error:
        result["failure_reason"] = f"{type(error).__name__}: {error}"
    return result


def prepare_completed_irc_reuse(admission: dict, workdir: Path):
    """Copy admitted evidence only; this function never invokes an executor."""
    if admission.get("kind") != "completed_irc_reuse_admission" or admission.get("valid") is not True:
        raise ValueError("Completed IRC copy requires valid source admission")
    references = {reference["path"]: reference for reference in admission["source_artifacts"]}
    for reference in references.values():
        _pinned(reference["path"], references)
    workdir = Path(workdir).resolve()
    source_root = Path(admission["source_workspace"]).resolve()
    if workdir.is_relative_to(source_root) or source_root.is_relative_to(workdir):
        raise ValueError("Completed IRC reuse requires a separate owned workspace")
    evidence_dir = workdir / "retained-irc-evidence"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    original_input = Path(admission["execution"]["input"])
    origins = {"input": original_input, "output": Path(admission["execution"]["output"]),
               "trajectory": Path(admission["trajectory"])}
    targets = {"input": evidence_dir / original_input.name,
               "output": evidence_dir / original_input.with_suffix(".out").name,
               "trajectory": workdir / "validated-trajectory.xyz"}
    for index, value in enumerate(admission["dependencies"]):
        origin = Path(value)
        relative = origin.relative_to(original_input.parent)
        if ".." in relative.parts:
            raise ValueError("Completed IRC evidence dependency leaves its owned directory")
        origins[f"dependency-{index}"] = origin
        targets[f"dependency-{index}"] = evidence_dir / relative
    copies = {}
    for key, origin in origins.items():
        _pinned(origin, references)
        target = targets[key]
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.read_bytes() != origin.read_bytes():
            raise ValueError("Completed IRC derived workspace contains conflicting evidence")
        if not target.exists():
            shutil.copy2(origin, target)
        if _digest(target) != references[str(origin.resolve())]["sha256"]:
            raise ValueError("Completed IRC evidence copy differs from admitted source")
        copies[key] = {"source": str(origin), "copy": str(target), "sha256": _digest(target)}
    evidence = {"schema_version": 1, "kind": "source_bound_completed_irc_reuse",
                "source_path_result": admission["source_path_result"], "copies": copies,
                "source_policy_id": admission.get("source_policy_id"),
                "reused_stages": admission["reused_stages"], "irc_reused": True,
                "endpoint_claims_reused": False, "renews_compute_budget": False,
                "claims_permitted": ["completed_continuous_irc_tracing"],
                "claims_prohibited": ["verified_path", "validated_minima", "rate", "lifetime"]}
    manifest = workdir / "completed-irc-reuse.json"
    if manifest.exists() and not _equivalent(json.loads(manifest.read_text()), evidence):
        raise ValueError("Completed IRC reuse evidence changed")
    if not manifest.exists():
        manifest.write_text(json.dumps(evidence, sort_keys=True, indent=2) + "\n")
    for reference in references.values():
        _pinned(reference["path"], references)
    execution = {**admission["execution"], "input": str(targets["input"]),
                 "output": str(targets["output"]), "resumed": True,
                 "read_only_reuse": True, "source_execution": admission["execution"]}
    return targets["input"], execution, admission["completion"], targets["trajectory"], evidence
