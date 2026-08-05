"""Resumable ORCA decomposition-path calculations and explanation artifacts."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path
import re
import shutil
from typing import Any, Callable

from src.inputgen import (
    create_orca_input,
    create_orca_irc_input,
    create_orca_neb_ts_input,
    create_orca_relaxed_scan_input,
)
from src.molecule_tools import smiles_to_xyz
from src.orca_runner import run_orca
from src.stability.freq_check import frequency_stability_check

from .decomposition_visuals import (
    render_candidate_storyboard,
    render_energy_profile,
    render_xyz_trajectory_gif,
    write_visual_manifest,
    xyz_frames,
)
from .route_geometry import (
    build_endpoint_complex_seeds,
    build_dissociation_endpoint,
    build_mapped_endpoint_seeds,
    interpolate_xyz,
    read_xyz,
    resolve_route_atom_mapping,
    validate_declared_connectivity,
    validate_route_mapping,
    validate_route_endpoint_connectivity,
    validate_trajectory_endpoints,
    write_xyz,
)
from .route_verify import RouteSpec, prepare_explanation
from .route_thermochemistry import (
    assemble_route_thermochemistry,
    condition_specific_forward_loss_bound,
)
from .workflow import run_optimization_and_frequency


def _write_json(path: Path, value: dict) -> Path:
    """Replace a JSON status file only after the new document is complete."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)
    return path


def _normal_orca_output(path: Path) -> bool:
    return Path(path).is_file() and "ORCA TERMINATED NORMALLY" in Path(path).read_text(errors="replace")


def _orca_optimization_converged(path: Path) -> bool:
    if not Path(path).is_file():
        return False
    text = Path(path).read_text(errors="replace").upper()
    return "THE OPTIMIZATION HAS CONVERGED" in text


def _orca_neb_no_barrier_outcome(path: Path) -> dict:
    """Recognize ORCA's semantically successful no-saddle NEB outcome.

    ORCA 6 can return a nonzero process status after a converged elastic band
    when the highest-energy image is an endpoint.  In that case it explicitly
    says that no barrier was found and skips the NEB-TS refinement.  The process
    status is not normal termination, but the retained MEP is still usable
    evidence once its images, energies, and endpoints are validated below.
    """
    output = Path(path)
    if not output.is_file():
        return {"status": "not_detected", "valid": False}
    text = output.read_text(errors="replace")
    converged = "THE NEB OPTIMIZATION HAS CONVERGED" in text
    no_barrier = bool(re.search(
        r"No\s+barrier\s+was\s+found\.\s*Skipping\s+NEB-TS\s+run\s+here\.",
        text,
        re.IGNORECASE,
    ))
    if not (converged and no_barrier):
        return {
            "status": "not_detected",
            "valid": False,
            "neb_converged": converged,
            "explicit_no_barrier_message": no_barrier,
        }
    return {
        "status": "converged_no_interior_barrier",
        "valid": True,
        "neb_converged": True,
        "explicit_no_barrier_message": True,
        "normal_termination": _normal_orca_output(output),
        "climbing_image_activated": "climbing image was never activated" not in text.lower(),
        "interpretation": (
            "ORCA converged the minimum-energy path and found no interior maximum, so it skipped "
            "NEB-TS. Reaction direction and any forward energetic threshold must be determined from "
            "the validated endpoint energy profile; this does not by itself establish a forward rate."
        ),
    }


_NEB_INTERMEDIATE_WARNING = re.compile(
    r"Possible\s+intermediate\s+minimum\s+found\s+at\s+image\(s\):",
    re.IGNORECASE,
)


def _neb_intermediate_image_indices(output: Path) -> list[int]:
    """Read ORCA's retained intermediate-image diagnostic without guessing chemistry."""
    output = Path(output)
    if not output.is_file():
        return []
    text = output.read_text(errors="replace")
    indices: list[int] = []
    for match in re.finditer(r"Possible\s+intermediate\s+minimum\s+found\s+at\s+image\(s\):\s*([^\n]+)", text, re.I):
        indices.extend(int(value) for value in re.findall(r"\d+", match.group(1)))
    return sorted(set(indices))


def _validate_neb_intermediates(
    route: RouteSpec,
    *,
    neb_output: Path,
    path_dir: Path,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
) -> dict:
    """Optimize/frequency-check ORCA-flagged NEB intermediate images.

    This intentionally establishes only a retained intermediate minimum and
    two *prepared* elementary segments.  It does not assign a rate or infer a
    reaction family from a geometry; each segment still needs its own TS/IRC.
    """
    indices = _neb_intermediate_image_indices(neb_output)
    if not indices:
        return {"status": "not_detected", "image_indices": [], "candidates": []}
    expected_elements, _, _ = read_xyz(path_dir / "reactant.xyz")
    candidates = []
    for image_index in indices:
        image = path_dir / f"neb-ts_im{image_index}.xyz"
        record: dict[str, Any] = {"image_index": image_index, "image_xyz": str(image)}
        if not image.is_file():
            record.update(status="image_missing")
            candidates.append(record)
            continue
        try:
            elements, _, _ = read_xyz(image)
            if elements != expected_elements:
                raise RuntimeError("NEB intermediate image changed atom count or order")
            calculation = _optimize_endpoint(
                image, path_dir / "intermediates" / f"image-{image_index:02d}",
                charge=route.charge, multiplicity=route.multiplicity,
                ncores=ncores, method_keywords=method_keywords,
                timeout_seconds=timeout_seconds,
            )
            frequency = calculation["frequency_check"]
            optimized = Path(calculation["optimized_xyz"])
            optimized_elements, _, _ = read_xyz(optimized)
            if optimized_elements != expected_elements:
                raise RuntimeError("optimized intermediate changed atom count or order")
            record.update(
                status="validated_intermediate_minimum" if frequency.get("IsMinimum") else "nonminimum_intermediate",
                optimized_xyz=str(optimized),
                optimization_output=str(calculation["optimization"]["out"]),
                frequency_output=str(calculation["frequency"]["out"]),
                frequency_check=frequency,
            )
        except Exception as error:
            record.update(status="intermediate_validation_failed", failure_reason=str(error))
        candidates.append(record)
    validated = [item for item in candidates if item.get("status") == "validated_intermediate_minimum"]
    return {
        "schema_version": 1,
        "kind": "orca_neb_intermediate_evidence",
        "status": "validated_intermediate_detected" if validated else "intermediate_detected_unvalidated",
        "image_indices": indices,
        "candidates": candidates,
        "validated_candidate_count": len(validated),
        "segmented_verification": {
            "status": "prepared" if validated else "not_prepared",
            "segments": [
                {"from": "declared_reactants", "to": item["optimized_xyz"]},
                {"from": item["optimized_xyz"], "to": "declared_products"},
            ] if validated else [],
            "requirements": [
                "optimize/validate each segment transition state",
                "require one relevant imaginary frequency per transition state",
                "run IRC or endpoint-connection validation for each segment",
                "derive the overall rate from the validated elementary sequence",
            ],
        },
        "interpretation": (
            "The intermediate is a validated ORCA stationary structure, not an automatically accepted mechanism or rate. "
            "The parent route is split into two explicitly retained segment-verification tasks."
        ),
    }


def _run_segmented_intermediate_nebs(
    route: RouteSpec,
    *,
    reactant_xyz: Path,
    product_xyz: Path,
    intermediate_evidence: dict,
    path_dir: Path,
    separated_fragments: bool,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    neb_images: int,
) -> dict:
    """Immediately launch smaller NEBs after ORCA flags an intermediate.

    These calculations are path evidence only.  Their endpoint/TS/IRC/rate
    validation remains a separate gate, so a segmented NEB can never by
    itself create an instability or lifetime claim.
    """
    candidates = [
        item for item in intermediate_evidence.get("candidates", [])
        if item.get("status") == "validated_intermediate_minimum"
        and item.get("optimized_xyz")
    ]
    if not candidates:
        return {
            "status": "not_prepared",
            "segments": [],
            "reason": "no_validated_intermediate_geometry",
        }
    segment_records: list[dict] = []
    # A shorter band is sufficient for each local segment and avoids repeating
    # the full expensive parent path resolution.
    segment_image_count = max(3, min(int(neb_images), 6))
    for intermediate in candidates:
        intermediate_xyz = Path(intermediate["optimized_xyz"])
        endpoints = [
            ("reactants-to-intermediate", Path(reactant_xyz), intermediate_xyz),
            ("intermediate-to-products", intermediate_xyz, Path(product_xyz)),
        ]
        for label, source, target in endpoints:
            segment_dir = path_dir / "segments" / label
            segment_dir.mkdir(parents=True, exist_ok=True)
            segment_reactant = segment_dir / "reactant.xyz"
            segment_product = segment_dir / "product.xyz"
            shutil.copy2(source, segment_reactant)
            shutil.copy2(target, segment_product)
            neb_input = create_orca_neb_ts_input(
                segment_reactant, segment_product,
                charge=route.charge, multiplicity=route.multiplicity,
                label="neb-ts", ncores=ncores, nimages=segment_image_count,
                method_keywords=method_keywords,
                preopt_ends=not separated_fragments,
            )
            try:
                execution = _run_or_resume(neb_input, timeout_seconds=timeout_seconds)
                trajectory = _locate_neb_trajectory(segment_dir, neb_input.stem)
                ts = _locate_neb_ts(segment_dir, neb_input.stem)
                segment_records.append({
                    "label": label,
                    "status": "computed",
                    "input": str(neb_input),
                    "output": execution.get("output"),
                    "execution": execution,
                    "trajectory": str(trajectory) if trajectory else None,
                    "ts_candidate": str(ts) if ts else None,
                    "rate_claim_allowed": False,
                })
            except Exception as error:
                segment_records.append({
                    "label": label,
                    "status": "failed",
                    "input": str(neb_input),
                    "failure_reason": f"{type(error).__name__}: {error}",
                    "rate_claim_allowed": False,
                })
    completed = sum(item.get("status") == "computed" for item in segment_records)
    return {
        "status": "computed" if completed == len(segment_records) else "incomplete",
        "segment_count": len(segment_records),
        "completed_segment_count": completed,
        "image_count_per_segment": segment_image_count,
        "segments": segment_records,
        "rate_claim_allowed": False,
        "requirements": [
            "validate each segment endpoint connectivity",
            "optimize each segment TS and require one relevant imaginary frequency",
            "run IRC or endpoint-connection validation for each segment",
            "derive the overall rate from the validated elementary sequence",
        ],
    }


def _run_or_resume(input_path: Path, *, timeout_seconds: float | None) -> dict:
    input_path = Path(input_path)
    output = Path(input_path).with_suffix(".out")
    input_text = input_path.read_text(errors="replace")
    referenced_names = set(re.findall(r'(?:Product|Hess_Filename)\s+"([^"]+)"', input_text))
    xyzfile = re.search(r"\*\s+xyzfile\s+[-+]?\d+\s+\d+\s+(\S+)", input_text, re.IGNORECASE)
    if xyzfile:
        referenced_names.add(xyzfile.group(1).strip('"'))
    digest = hashlib.sha256(input_path.read_bytes())
    referenced = []
    for name in sorted(referenced_names):
        path = input_path.parent / name
        if not path.is_file():
            raise FileNotFoundError(f"ORCA input references a missing artifact: {path}")
        digest.update(name.encode())
        digest.update(path.read_bytes())
        referenced.append(str(path))
    contract_path = input_path.with_suffix(".contract.json")
    contract = json.loads(contract_path.read_text()) if contract_path.is_file() else {}
    no_barrier_outcome = _orca_neb_no_barrier_outcome(output)
    existing_intermediates = _neb_intermediate_image_indices(output)
    if contract.get("sha256") == digest.hexdigest():
        if _normal_orca_output(output):
            return {"status": "completed", "resumed": True, "input": str(input_path), "output": str(output)}
        if no_barrier_outcome["valid"]:
            return {
                "status": "completed_no_interior_barrier", "resumed": True,
                "input": str(input_path), "output": str(output),
                "orca_neb_outcome": no_barrier_outcome,
            }
        if existing_intermediates:
            return {
                "status": "completed_intermediate_detected", "resumed": True,
                "input": str(input_path), "output": str(output),
                "intermediate_image_indices": existing_intermediates,
            }

    is_neb = "NEB-TS" in input_text.upper()

    def _stream(line: str) -> bool | None:
        # Preserve ORCA's native output verbatim while allowing a NEB warning
        # to interrupt the expensive single-band calculation immediately.
        print(line, end="", flush=True)
        return bool(is_neb and _NEB_INTERMEDIATE_WARNING.search(line))

    try:
        artifacts = run_orca(
            input_path,
            timeout_seconds=timeout_seconds,
            stream_output=_stream,
            live_output=False,
        )
    except RuntimeError:
        no_barrier_outcome = _orca_neb_no_barrier_outcome(output)
        existing_intermediates = _neb_intermediate_image_indices(output)
        if not no_barrier_outcome["valid"] and not existing_intermediates:
            raise
        artifacts = {"out": output}
    normal_termination = _normal_orca_output(Path(artifacts["out"]))
    no_barrier_outcome = _orca_neb_no_barrier_outcome(Path(artifacts["out"]))
    intermediate_indices = _neb_intermediate_image_indices(Path(artifacts["out"]))
    interrupted_for_intermediate = bool(artifacts.get("stopped_early")) and bool(intermediate_indices)
    if not normal_termination and not no_barrier_outcome["valid"] and not interrupted_for_intermediate:
        raise RuntimeError(f"ORCA returned without its normal-termination marker: {artifacts['out']}")
    _write_json(contract_path, {
        "schema_version": 1,
        "kind": "orca_job_input_contract",
        "sha256": digest.hexdigest(),
        "input": str(input_path),
        "referenced_files": referenced,
    })
    status = (
        "completed_intermediate_detected" if interrupted_for_intermediate else
        "completed" if normal_termination else "completed_no_interior_barrier"
    )
    result = {
        "status": status, "resumed": False,
        "input": str(input_path), "output": str(artifacts["out"]),
    }
    if no_barrier_outcome["valid"]:
        result["orca_neb_outcome"] = no_barrier_outcome
    if intermediate_indices:
        result["intermediate_image_indices"] = intermediate_indices
    if interrupted_for_intermediate:
        result["early_stop"] = {
            "status": "stopped_after_orca_intermediate_warning",
            "reason": "ORCA reported an intermediate minimum during the live NEB stream",
            "image_indices": intermediate_indices,
        }
    return result


def _combine_xyz_files(files: list[Path], output: Path) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(path.read_text().rstrip() + "\n" for path in files))
    return output


def _scan_frames(workdir: Path, label: str) -> list[Path]:
    pattern = re.compile(rf"^{re.escape(label)}\.(\d+)\.xyz$")
    matches = []
    for path in Path(workdir).glob(f"{label}.*.xyz"):
        match = pattern.match(path.name)
        if match:
            matches.append((int(match.group(1)), path))
    return [path for _, path in sorted(matches)]


def _scan_energies(workdir: Path, label: str, expected: int, scan_files: list[Path] | None = None) -> list[float] | None:
    candidates = [Path(workdir) / f"{label}.relaxscanact.dat", Path(workdir) / f"{label}_relaxscanact.dat"]
    path = next((candidate for candidate in candidates if candidate.is_file()), None)
    energies = []
    if path is not None:
        for line in path.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            values = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", line)
            if len(values) >= 2:
                energies.append(float(values[-1]))
    elif scan_files:
        for scan_file in scan_files:
            lines = scan_file.read_text().splitlines()
            match = re.search(r"\bE\s+(-?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)", lines[1] if len(lines) > 1 else "")
            if match:
                energies.append(float(match.group(1)))
    return energies if len(energies) == expected else None


def _validate_dissociation_path(trajectory: Path, bond: tuple[int, int],
                                final_distance: float) -> dict:
    frames = xyz_frames(trajectory)
    elements = frames[0][0]
    if any(frame[0] != elements for frame in frames):
        return {"status": "atom_order_changed", "valid": False}
    left, right = bond
    distances = [
        float(((coordinates[left] - coordinates[right]) ** 2).sum() ** 0.5)
        for _, coordinates, _ in frames
    ]
    decreases = sum(current + 0.05 < previous for previous, current in zip(distances, distances[1:]))
    valid = len(frames) >= 3 and distances[-1] >= final_distance - 0.2 and decreases <= max(1, len(distances) // 10)
    return {
        "status": "validated" if valid else "dissociation_coordinate_mismatch",
        "valid": valid,
        "frame_count": len(frames),
        "bond_atom_indices": list(bond),
        "bond_distances_angstrom": distances,
    }


def _locate_neb_ts(workdir: Path, label: str) -> Path | None:
    candidates = [
        *Path(workdir).glob(f"{label}*NEB-TS*converged*.xyz"),
        *Path(workdir).glob(f"{label}*TS*.xyz"),
        Path(workdir) / f"{label}.xyz",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _locate_irc_trajectory(workdir: Path, label: str) -> Path | None:
    candidates = [
        *Path(workdir).glob(f"{label}*_IRC_Full_trj.xyz"),
        *Path(workdir).glob("*_IRC_Full_trj.xyz"),
    ]
    return next((path for path in candidates if path.is_file()), None)


def _optimize_endpoint(seed: Path, workdir: Path, *, charge: int, multiplicity: int,
                       ncores: int, method_keywords: list[str] | None,
                       timeout_seconds: float | None) -> dict:
    workdir.mkdir(parents=True, exist_ok=True)
    local_seed = workdir / "input.xyz"
    if seed.resolve() != local_seed.resolve():
        shutil.copy2(seed, local_seed)
    contract_path = workdir / "endpoint-contract.json"
    contract_value = {
        "schema_version": 1,
        "kind": "orca_endpoint_contract",
        "input_xyz_sha256": hashlib.sha256(local_seed.read_bytes()).hexdigest(),
        "charge": int(charge),
        "multiplicity": int(multiplicity),
        "ncores": int(ncores),
        "method_keywords": list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"]),
    }
    retained_contract = json.loads(contract_path.read_text()) if contract_path.is_file() else None
    optimized = workdir / "opt.xyz"
    frequency = workdir / "freq.out"
    optimization_output = workdir / "opt.out"
    if (retained_contract == contract_value and optimized.is_file() and _normal_orca_output(optimization_output)
            and _orca_optimization_converged(optimization_output) and _normal_orca_output(frequency)):
        input_elements, _, _ = read_xyz(local_seed)
        optimized_elements, _, _ = read_xyz(optimized)
        if input_elements != optimized_elements:
            raise RuntimeError("Cached ORCA endpoint changed composition or atom order")
        return {
            "optimized_xyz": optimized,
            "optimization": {"out": optimization_output},
            "frequency": {"out": frequency},
            "frequency_check": frequency_stability_check(frequency),
            "resumed": True,
        }
    result = run_optimization_and_frequency(
        local_seed,
        workdir,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        run_frequency=True,
        method_keywords=method_keywords,
        timeout_seconds=timeout_seconds,
    )
    if not _orca_optimization_converged(Path(result["optimization"]["out"])):
        raise RuntimeError(f"ORCA optimization did not reach its convergence marker: {result['optimization']['out']}")
    input_elements, _, _ = read_xyz(local_seed)
    optimized_elements, _, _ = read_xyz(Path(result["optimized_xyz"]))
    if input_elements != optimized_elements:
        raise RuntimeError("ORCA optimized geometry changed endpoint composition or atom order")
    _write_json(contract_path, contract_value)
    return result


def _render_preview(route: RouteSpec, endpoints: dict, visual_dir: Path, storyboard: Path) -> dict:
    trajectory = interpolate_xyz(
        Path(endpoints["reactant_xyz"]),
        Path(endpoints["product_xyz"]),
        visual_dir / "candidate-path-preview.xyz",
        comment_prefix="Candidate interpolation; not an ORCA reaction path",
    )
    animation = render_xyz_trajectory_gif(
        trajectory,
        visual_dir / "candidate-path-preview.gif",
        title=f"{route.route_id}: candidate geometry interpolation",
        evidence_label="CANDIDATE PREVIEW — NOT A COMPUTED REACTION PATH",
    )
    manifest = write_visual_manifest(
        visual_dir / "visuals.json",
        route_id=route.route_id,
        evidence_level="candidate_preview",
        storyboard=storyboard,
        animation=animation,
        trajectory=trajectory,
    )
    return {"storyboard": str(storyboard), "animation": str(animation),
            "trajectory": str(trajectory), "manifest": str(manifest)}


def _run_dissociation(route: RouteSpec, explanation_dir: Path, *, prepare_only: bool,
                      ncores: int, method_keywords: list[str] | None,
                      timeout_seconds: float | None, scan_steps: int) -> dict:
    endpoint_dir = explanation_dir / "endpoints"
    endpoint_dir.mkdir(parents=True, exist_ok=True)
    seed = endpoint_dir / "parent-seed.xyz"
    if not seed.is_file():
        smiles_to_xyz(route.parent_smiles, seed)
    if prepare_only:
        reactant = seed
    else:
        endpoint_result = _optimize_endpoint(
            seed, endpoint_dir / "reactant", charge=route.charge, multiplicity=route.multiplicity,
            ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
        )
        if not endpoint_result["frequency_check"]["IsMinimum"]:
            raise RuntimeError("The parent endpoint is not an ORCA local minimum")
        reactant = Path(endpoint_result["optimized_xyz"])
    endpoints = build_dissociation_endpoint(
        route, reactant, endpoint_dir / "product-mapped.xyz"
    )
    if prepare_only:
        return {"status": "prepared", "endpoints": endpoints}

    path_dir = explanation_dir / "path"
    path_dir.mkdir(parents=True, exist_ok=True)
    path_reactant = path_dir / "reactant.xyz"
    shutil.copy2(reactant, path_reactant)
    start = endpoints["initial_distance_angstrom"]
    end = endpoints["final_distance_angstrom"]
    scan_input = create_orca_relaxed_scan_input(
        path_reactant,
        bond_atom_indices=route.broken_bonds[0],
        start_distance_angstrom=start,
        end_distance_angstrom=end,
        charge=route.charge,
        multiplicity=route.multiplicity,
        ncores=ncores,
        steps=scan_steps,
        method_keywords=method_keywords,
    )
    scan_files = _scan_frames(path_dir, scan_input.stem)
    abnormal_reason = None
    if len(scan_files) >= scan_steps and scan_input.with_suffix(".out").is_file():
        normal = _normal_orca_output(scan_input.with_suffix(".out"))
        execution = {
            "status": "completed" if normal else "incomplete",
            "resumed": True,
            "input": str(scan_input),
            "output": str(scan_input.with_suffix(".out")),
        }
        if not normal:
            abnormal_reason = "ORCA retained every requested scan geometry but did not terminate normally"
    else:
        try:
            execution = _run_or_resume(scan_input, timeout_seconds=timeout_seconds)
        except Exception as error:
            abnormal_reason = str(error)
            execution = {
                "status": "incomplete",
                "resumed": False,
                "input": str(scan_input),
                "output": str(scan_input.with_suffix(".out")),
                "failure_reason": abnormal_reason,
            }
        scan_files = _scan_frames(path_dir, scan_input.stem)
    if len(scan_files) < 3:
        raise RuntimeError(abnormal_reason or "ORCA completed the relaxed scan without at least three retained scan geometries")
    trajectory = _combine_xyz_files(scan_files, path_dir / "trajectory.xyz")
    validation = _validate_dissociation_path(trajectory, route.broken_bonds[0], end)
    if not validation["valid"]:
        raise RuntimeError("The retained ORCA scan does not follow the declared dissociation coordinate")
    energies = _scan_energies(path_dir, scan_input.stem, len(scan_files), scan_files)
    complete = abnormal_reason is None and len(scan_files) >= scan_steps
    return {
        "status": "computed_dissociation_path" if complete else "incomplete_dissociation_path",
        "endpoints": endpoints,
        "execution": execution,
        "trajectory": str(trajectory),
        "energies_hartree": energies,
        "validation": validation,
        "failure_reason": abnormal_reason,
        "interpretation": (
            "A relaxed ground-state bond-distance scan was retained. No transition state or excited-state "
            "photochemical dynamics is claimed. An incomplete label means at least one requested point or "
            "ORCA normal termination/convergence check was missing."
        ),
    }


def classify_path_energy_profile(
    energies_hartree: list[float] | None,
    *,
    reactant_fragment_count: int,
    product_fragment_count: int,
    barrier_tolerance_hartree: float = 5.0e-4,
) -> dict:
    """Classify a complete path-energy profile without assuming a family.

    This is a hook-level physical classification, not a rate model.  In
    particular, a downhill association curve still requires capture and
    recrossing kinetics before it can contribute to a lifetime.
    """
    import math

    if not energies_hartree or len(energies_hartree) < 5:
        return {"status": "surface_unresolved", "reason": "fewer_than_five_energy_points"}
    energies = [float(value) for value in energies_hartree]
    if any(not math.isfinite(value) for value in energies):
        return {"status": "surface_unresolved", "reason": "nonfinite_energy"}
    if barrier_tolerance_hartree <= 0 or not math.isfinite(barrier_tolerance_hartree):
        raise ValueError("Barrier tolerance must be finite and positive")
    largest_step = max(abs(right - left) for left, right in zip(energies, energies[1:]))
    if largest_step > 0.2:
        return {
            "status": "surface_unresolved",
            "reason": "discontinuous_energy_profile",
            "largest_adjacent_step_hartree": largest_step,
        }
    interior_maximum = max(energies[1:-1])
    endpoint_maximum = max(energies[0], energies[-1])
    if interior_maximum > endpoint_maximum + barrier_tolerance_hartree:
        maximum_index = max(range(1, len(energies) - 1), key=energies.__getitem__)
        return {
            "status": "barriered_ts_candidate",
            "interior_maximum_index": maximum_index,
            "interior_barrier_hartree": interior_maximum - energies[0],
            "largest_adjacent_step_hartree": largest_step,
        }
    direction_change = energies[-1] - energies[0]
    if reactant_fragment_count > product_fragment_count and direction_change <= barrier_tolerance_hartree:
        status = "barrierless_capture_candidate"
    elif reactant_fragment_count < product_fragment_count and direction_change >= -barrier_tolerance_hartree:
        status = "barrierless_dissociation_candidate"
    elif reactant_fragment_count == product_fragment_count and direction_change <= barrier_tolerance_hartree:
        # A two-to-two profile can be barrierless only in the asserted forward
        # direction.  Calling an uphill encounter path "barrierless" would
        # incorrectly send it to a capture-rate model.
        status = "barrierless_path_candidate"
    else:
        status = "surface_unresolved"
    downhill_direction = (
        "forward" if direction_change < -barrier_tolerance_hartree
        else "reverse" if direction_change > barrier_tolerance_hartree
        else "approximately_thermoneutral"
    )
    return {
        "status": status,
        "reaction_energy_hartree": direction_change,
        "downhill_direction": downhill_direction,
        "forward_energetic_threshold_hartree": max(0.0, direction_change),
        "forward_target_loss_barrierless": bool(direction_change <= barrier_tolerance_hartree),
        "largest_adjacent_step_hartree": largest_step,
        "reason": (
            None if status != "surface_unresolved" else
            "forward_path_uphill_without_interior_ts_candidate"
            if reactant_fragment_count == product_fragment_count and direction_change > barrier_tolerance_hartree
            else "profile_direction_inconsistent_with_stoichiometry"
        ),
    }


def _trajectory_energies(path: Path) -> list[float] | None:
    values: list[float] = []
    for _, _, comment in xyz_frames(path):
        matches = re.findall(
            r"(?:energy|\bE)\s*(?:=|:)?\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
            comment,
            re.IGNORECASE,
        )
        if not matches:
            return None
        values.append(float(matches[-1]))
    return values


def _locate_neb_trajectory(workdir: Path, label: str) -> Path | None:
    candidates = [
        *Path(workdir).glob(f"{label}*_MEP_trj.xyz"),
        *Path(workdir).glob("*_MEP_trj.xyz"),
        *Path(workdir).glob(f"{label}*NEB*trj.xyz"),
    ]
    return next((path for path in candidates if path.is_file()), None)


def _create_ts_optimization_input(
    seed_xyz: Path,
    *,
    charge: int,
    multiplicity: int,
    ncores: int,
    method_keywords: list[str] | None,
) -> Path:
    """Create an explicit ORCA OptTS refinement with a fresh Hessian."""
    if ncores < 1:
        raise ValueError("ncores must be at least one")
    seed_xyz = Path(seed_xyz)
    method = list(method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"])
    output = seed_xyz.parent / "ts-opt.inp"
    output.write_text(
        "! " + " ".join([*method, "TightSCF", "OptTS"]) + "\n"
        + f"%pal\n  nprocs {ncores}\nend\n"
        + "%geom\n  Calc_Hess true\n  Recalc_Hess 5\nend\n"
        + f"* xyzfile {charge} {multiplicity} {seed_xyz.name}\n"
    )
    return output


def _stationary_fragment_records(
    route: RouteSpec,
    workdir: Path,
    mapping: dict,
    *,
    execute: bool,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
) -> dict:
    """Prepare/validate separated species needed by Arkane."""
    surface = mapping["surface"]
    records: dict[str, list[dict]] = {"reactants": [], "products": []}
    for side, smiles_values, declared_labels, charges, multiplicities in (
        ("reactants", route.reactant_smiles, route.reactant_labels, surface["reactant_fragment_charges"],
         surface["reactant_fragment_multiplicities"]),
        ("products", route.product_smiles, route.product_labels, surface["product_fragment_charges"],
         surface["product_fragment_multiplicities"]),
    ):
        labels = list(declared_labels) if len(declared_labels) == len(smiles_values) else list(smiles_values)
        totals = {label: labels.count(label) for label in set(labels)}
        seen: dict[str, int] = {}
        label_signatures: dict[str, list[tuple[str, int, int]]] = {}
        cached: dict[tuple[str, str, int, int], dict] = {}
        for index, (smiles, charge, multiplicity) in enumerate(
            zip(smiles_values, charges, multiplicities), start=1,
        ):
            original_label = str(labels[index - 1])
            seen[original_label] = seen.get(original_label, 0) + 1
            signature = (smiles, int(charge), int(multiplicity))
            distinct = label_signatures.setdefault(original_label, [])
            if signature not in distinct:
                distinct.append(signature)
            # Repeated stoichiometric copies intentionally reuse one species
            # label/job.  A suffix is reserved for a genuine label collision
            # between different structures or electronic states.
            collision_index = distinct.index(signature) + 1
            label = original_label if collision_index == 1 else f"{original_label}__collision_{collision_index}"
            cache_key = (label, *signature)
            if cache_key in cached:
                record = {
                    **cached[cache_key],
                    "stoichiometric_copy": seen[original_label],
                    "reused_stationary_point": True,
                }
                records[side].append(record)
                continue
            folder = Path(workdir) / side / f"{side[:-1]}-{index:02d}"
            folder.mkdir(parents=True, exist_ok=True)
            seed = folder / "input.xyz"
            if not seed.is_file():
                smiles_to_xyz(smiles, seed)
            record = {
                "label": label,
                "original_label": original_label,
                "stoichiometric_copy": seen[original_label],
                "stoichiometric_copy_count": totals[original_label],
                "smiles": smiles,
                "charge": int(charge),
                "multiplicity": int(multiplicity),
                "input_xyz": str(seed),
                "status": "prepared",
            }
            if execute:
                print(
                    f"[STORCA] Route {route.route_id}: validating {side[:-1]} species "
                    f"{index}/{len(smiles_values)} ({label})", flush=True,
                )
                try:
                    calculation = _optimize_endpoint(
                        seed, folder, charge=int(charge), multiplicity=int(multiplicity),
                        ncores=ncores, method_keywords=method_keywords,
                        timeout_seconds=timeout_seconds,
                    )
                    check = calculation["frequency_check"]
                    connectivity = validate_declared_connectivity(
                        Path(calculation["optimized_xyz"]), (smiles,),
                    )
                    if not check.get("IsMinimum"):
                        endpoint_status = "nonminimum_endpoint"
                    elif not connectivity.get("valid"):
                        endpoint_status = "endpoint_connectivity_validation_failed"
                    else:
                        endpoint_status = "validated_minimum"
                    record.update(
                        status=endpoint_status,
                        optimized_xyz=str(calculation["optimized_xyz"]),
                        optimization_output=str(calculation["optimization"]["out"]),
                        frequency_output=str(calculation["frequency"]["out"]),
                        frequency_check=check,
                        connectivity_check=connectivity,
                    )
                except Exception as error:
                    record.update(status="failed", failure_reason=str(error))
                print(
                    f"[STORCA] Route {route.route_id}: {side[:-1]} species "
                    f"{index}/{len(smiles_values)} finished ({record['status']})", flush=True,
                )
            cached[cache_key] = dict(record)
            records[side].append(record)
    records["valid"] = bool(execute) and all(
        record["status"] == "validated_minimum"
        for side in ("reactants", "products")
        for record in records[side]
    )
    records["status"] = (
        "validated" if records["valid"] else "prepared" if not execute else "stationary_point_validation_failed"
    )
    return records


def _route_uses_separated_fragment_endpoints(route: RouteSpec) -> bool:
    """Whether either reaction side is an explicitly separated channel."""
    return len(route.reactant_smiles) > 1 or len(route.product_smiles) > 1


def _stationary_fragment_xyz(stationary_points: dict, *, optimized: bool) -> dict[str, list[str]]:
    """Extract one atom-ordered geometry for every stoichiometric fragment."""
    key = "optimized_xyz" if optimized else "input_xyz"
    paths: dict[str, list[str]] = {"reactants": [], "products": []}
    for side in paths:
        for record in stationary_points.get(side) or []:
            value = record.get(key)
            if not value or not Path(value).is_file():
                raise FileNotFoundError(
                    f"Stationary-point record for {record.get('label') or side} has no retained {key}"
                )
            paths[side].append(str(value))
    return paths


def _validate_irc_route_endpoints(
    route: RouteSpec,
    trajectory: Path,
    reactant_xyz: Path,
    product_xyz: Path,
    *,
    mapping_result: dict,
    separated_fragments: bool,
) -> dict:
    """Validate IRC ends by geometry for bound routes and topology for channels."""
    if not separated_fragments:
        return validate_trajectory_endpoints(trajectory, reactant_xyz, product_xyz)
    frames = xyz_frames(trajectory)
    if len(frames) < 2 or frames[0][0] != frames[-1][0]:
        return {"status": "irc_endpoint_frames_incomplete", "valid": False}
    first = Path(trajectory).with_name("irc-endpoint-first.xyz")
    last = Path(trajectory).with_name("irc-endpoint-last.xyz")
    write_xyz(first, frames[0][0], frames[0][1], "IRC endpoint; first retained frame")
    write_xyz(last, frames[-1][0], frames[-1][1], "IRC endpoint; last retained frame")
    forward = {
        "reactant": validate_route_endpoint_connectivity(
            route, first, side="reactants", mapping_result=mapping_result,
        ),
        "product": validate_route_endpoint_connectivity(
            route, last, side="products", mapping_result=mapping_result,
        ),
    }
    reverse = {
        "reactant": validate_route_endpoint_connectivity(
            route, last, side="reactants", mapping_result=mapping_result,
        ),
        "product": validate_route_endpoint_connectivity(
            route, first, side="products", mapping_result=mapping_result,
        ),
    }
    forward_valid = forward["reactant"].get("valid") and forward["product"].get("valid")
    reverse_valid = reverse["reactant"].get("valid") and reverse["product"].get("valid")
    return {
        "schema_version": 1,
        "kind": "separated_fragment_irc_endpoint_check",
        "status": "validated" if forward_valid or reverse_valid else "declared_fragment_channels_not_reached",
        "valid": bool(forward_valid or reverse_valid),
        "matched_direction": "forward" if forward_valid else "reverse" if reverse_valid else None,
        "first_endpoint_xyz": str(first),
        "last_endpoint_xyz": str(last),
        "forward_checks": forward,
        "reverse_checks": reverse,
        "interpretation": (
            "Separated channels are matched by declared fragment connectivity; "
            "a freely rotating/translating encounter geometry is not required to reproduce an arbitrary seed RMSD."
        ),
    }


def _run_generic_orientation(
    route: RouteSpec,
    seed_record: dict,
    *,
    mapping_result: dict,
    ncores: int,
    method_keywords: list[str] | None,
    timeout_seconds: float | None,
    neb_images: int,
    barrierless_classifier: Callable[[dict], dict] | None,
    separated_fragments: bool = False,
    stationary_points: dict | None = None,
) -> dict:
    """Run one endpoint orientation through NEB-TS, TS/frequency, and IRC."""
    orientation = int(seed_record["orientation"])
    folder = Path(seed_record["reactant_xyz"]).parent
    record: dict[str, Any] = {"orientation": orientation, "status": "running"}
    try:
        if separated_fragments:
            if not stationary_points or not stationary_points.get("valid"):
                raise RuntimeError("Separated-fragment NEB requires validated stationary species")
            reactant_source = Path(seed_record["reactant_xyz"])
            product_source = Path(seed_record["product_xyz"])
            reactant_check = {
                "status": "not_applicable_to_separated_channel",
                "IsMinimum": None,
                "reason": "The encounter geometry is assembled from independently validated fragment minima.",
            }
            product_check = dict(reactant_check)
        else:
            reactant_result = _optimize_endpoint(
                Path(seed_record["reactant_xyz"]), folder / "reactant-complex",
                charge=route.charge, multiplicity=route.multiplicity, ncores=ncores,
                method_keywords=method_keywords, timeout_seconds=timeout_seconds,
            )
            product_result = _optimize_endpoint(
                Path(seed_record["product_xyz"]), folder / "product-complex",
                charge=route.charge, multiplicity=route.multiplicity, ncores=ncores,
                method_keywords=method_keywords, timeout_seconds=timeout_seconds,
            )
            reactant_source = Path(reactant_result["optimized_xyz"])
            product_source = Path(product_result["optimized_xyz"])
            reactant_check = reactant_result["frequency_check"]
            product_check = product_result["frequency_check"]
        reactant_connectivity = validate_route_endpoint_connectivity(
            route, reactant_source,
            side="reactants", mapping_result=mapping_result,
        )
        product_connectivity = validate_route_endpoint_connectivity(
            route, product_source,
            side="products", mapping_result=mapping_result,
        )
        endpoint_validation = {
            "mode": (
                "assembled_from_validated_separated_fragment_minima"
                if separated_fragments else "optimized_bound_endpoint_minima"
            ),
            "full_complex_minimum_required": not separated_fragments,
            "stationary_fragment_validation_status": (
                stationary_points.get("status") if stationary_points else "not_applicable"
            ),
            "reactant": {
                "endpoint_xyz": str(reactant_source),
                "frequency_check": reactant_check,
                "connectivity_check": reactant_connectivity,
            },
            "product": {
                "endpoint_xyz": str(product_source),
                "frequency_check": product_check,
                "connectivity_check": product_connectivity,
            },
            "valid": bool(
                (separated_fragments or (
                    reactant_check.get("IsMinimum") and product_check.get("IsMinimum")
                ))
                and reactant_connectivity.get("valid") and product_connectivity.get("valid")
            ),
        }
        if not separated_fragments:
            endpoint_validation["reactant"].update(
                optimization_output=str(reactant_result["optimization"]["out"]),
                frequency_output=str(reactant_result["frequency"]["out"]),
            )
            endpoint_validation["product"].update(
                optimization_output=str(product_result["optimization"]["out"]),
                frequency_output=str(product_result["frequency"]["out"]),
            )
        record["endpoint_validation"] = endpoint_validation
        if not endpoint_validation["valid"]:
            connectivity_failed = not (
                reactant_connectivity.get("valid") and product_connectivity.get("valid")
            )
            record.update(
                status=(
                    "endpoint_connectivity_validation_failed"
                    if connectivity_failed else "endpoint_validation_failed"
                ),
                path_classification="surface_unresolved",
            )
            return record
        path_dir = folder / "path"
        path_dir.mkdir(parents=True, exist_ok=True)
        reactant_xyz = path_dir / "reactant.xyz"
        product_xyz = path_dir / "product.xyz"
        shutil.copy2(reactant_source, reactant_xyz)
        shutil.copy2(product_source, product_xyz)
        # These endpoints have already been independently validated as
        # minima.  Asking ORCA to preoptimize them once more makes the NEB
        # robust to small residual forces.  Encounter geometries deliberately
        # combine separated fragments at a selected orientation, so preopting
        # those ends would change the reaction channel and is prohibited.
        preopt_ends = not separated_fragments
        neb_input = create_orca_neb_ts_input(
            reactant_xyz, product_xyz, charge=route.charge, multiplicity=route.multiplicity,
            label="neb-ts", ncores=ncores, nimages=neb_images, method_keywords=method_keywords,
            preopt_ends=preopt_ends,
        )
        record["neb_endpoint_policy"] = {
            "preopt_ends": preopt_ends,
            "free_end_neb": False,
            "reason": (
                "both declared endpoints are bound, validated minima"
                if preopt_ends else
                "assembled separated-fragment encounter endpoints must retain their declared channel"
            ),
        }
        try:
            neb_execution = _run_or_resume(neb_input, timeout_seconds=timeout_seconds)
        except Exception as error:
            neb_execution = {
                "status": "failed", "input": str(neb_input),
                "output": str(neb_input.with_suffix(".out")), "failure_reason": str(error),
            }
        record["neb"] = neb_execution
        semantic_no_barrier = neb_execution["status"] == "completed_no_interior_barrier"
        ts_xyz = _locate_neb_ts(path_dir, neb_input.stem) if neb_execution["status"] == "completed" else None
        if ts_xyz is None:
            intermediate_evidence = _validate_neb_intermediates(
                route, neb_output=Path(neb_execution["output"]), path_dir=path_dir,
                ncores=ncores, method_keywords=method_keywords,
                timeout_seconds=timeout_seconds,
            )
            if intermediate_evidence.get("status") == "validated_intermediate_detected":
                intermediate_evidence["segmented_verification"] = _run_segmented_intermediate_nebs(
                    route,
                    reactant_xyz=reactant_xyz,
                    product_xyz=product_xyz,
                    intermediate_evidence=intermediate_evidence,
                    path_dir=path_dir,
                    separated_fragments=separated_fragments,
                    ncores=ncores,
                    method_keywords=method_keywords,
                    timeout_seconds=timeout_seconds,
                    neb_images=neb_images,
                )
                record.update(
                    status="intermediate_detected_requires_segmented_verification",
                    path_classification="multistep_intermediate_detected",
                    classification_evidence={
                        "status": "multistep_intermediate_detected",
                        "reason": "orca_neb_detected_and_validated_intermediate_minimum",
                    },
                    intermediate_evidence=intermediate_evidence,
                    path_evidence={
                        "kind": "orca_neb_intermediate_path",
                        "trajectory": str(_locate_neb_trajectory(path_dir, neb_input.stem) or "") or None,
                        "rate_claim_allowed": False,
                    },
                    trajectory=None,
                )
                return record
            trajectory = _locate_neb_trajectory(path_dir, neb_input.stem)
            energies = _trajectory_energies(trajectory) if trajectory else None
            frame_count = len(xyz_frames(trajectory)) if trajectory else 0
            trajectory_validation = (
                validate_trajectory_endpoints(trajectory, reactant_xyz, product_xyz)
                if trajectory else None
            )
            evidence_complete = bool(
                neb_execution["status"] in {"completed", "completed_no_interior_barrier"}
                and trajectory is not None
                and frame_count >= neb_images
                and energies is not None
                and len(energies) == frame_count
                and trajectory_validation
                and trajectory_validation["valid"]
            )
            if not evidence_complete:
                record.update(
                    status="surface_unresolved",
                    path_classification="surface_unresolved",
                    classification_evidence={
                        "status": "surface_unresolved",
                        "reason": "incomplete_neb_path_evidence",
                        "normal_termination": neb_execution["status"] == "completed",
                        "semantic_no_barrier_completion": semantic_no_barrier,
                        "orca_neb_outcome": neb_execution.get("orca_neb_outcome"),
                        "expected_minimum_image_count": neb_images,
                        "retained_image_count": frame_count,
                        "energies_complete": energies is not None and len(energies) == frame_count,
                        "endpoint_validation": trajectory_validation,
                    },
                    path_evidence={
                        "kind": "incomplete_neb_path",
                        "trajectory": str(trajectory) if trajectory else None,
                        "energies_hartree": energies,
                        "rate_claim_allowed": False,
                    },
                    trajectory=None,
                )
                return record
            context = {
                "route": route.as_dict(), "orientation": orientation,
                "trajectory": str(trajectory) if trajectory else None,
                "energies_hartree": energies,
                "orca_neb_outcome": neb_execution.get("orca_neb_outcome"),
                "default_classification": classify_path_energy_profile(
                    energies,
                    reactant_fragment_count=len(route.reactant_smiles),
                    product_fragment_count=len(route.product_smiles),
                ),
            }
            classification = barrierless_classifier(context) if barrierless_classifier else context["default_classification"]
            allowed = {
                "barriered_ts_candidate", "barrierless_capture_candidate",
                "barrierless_dissociation_candidate", "barrierless_path_candidate", "surface_unresolved",
            }
            if not isinstance(classification, dict) or classification.get("status") not in allowed:
                raise ValueError("Barrierless-classification hook returned an unsupported status")
            record.update(
                status=("classified_without_transition_state" if classification["status"] != "surface_unresolved"
                        else "surface_unresolved"),
                path_classification=classification["status"],
                classification_evidence=classification,
                path_evidence={
                    "kind": (
                        "orca_converged_neb_without_interior_barrier"
                        if semantic_no_barrier else "neb_energy_path_without_validated_ts"
                    ),
                    "trajectory": str(trajectory) if trajectory else None,
                    "energies_hartree": energies,
                    "endpoint_validation": trajectory_validation,
                    "orca_neb_outcome": neb_execution.get("orca_neb_outcome"),
                    "rate_claim_allowed": False,
                },
                trajectory=str(trajectory) if trajectory and trajectory_validation and trajectory_validation["valid"] else None,
            )
            return record

        ts_dir = folder / "transition-state"
        ts_dir.mkdir(parents=True, exist_ok=True)
        ts_seed = ts_dir / "ts-seed.xyz"
        shutil.copy2(ts_xyz, ts_seed)
        ts_input = _create_ts_optimization_input(
            ts_seed, charge=route.charge, multiplicity=route.multiplicity,
            ncores=ncores, method_keywords=method_keywords,
        )
        ts_execution = _run_or_resume(ts_input, timeout_seconds=timeout_seconds)
        if not _orca_optimization_converged(Path(ts_execution["output"])):
            raise RuntimeError("ORCA OptTS terminated normally without converging the transition state")
        optimized_ts = ts_input.with_suffix(".xyz")
        if not optimized_ts.is_file():
            raise RuntimeError("ORCA OptTS completed without an optimized TS geometry")
        frequency_input = create_orca_input(
            optimized_ts, route.charge, route.multiplicity, freq=True,
            label="ts-frequency", ncores=ncores, method_keywords=method_keywords,
        )
        frequency_execution = _run_or_resume(frequency_input, timeout_seconds=timeout_seconds)
        frequency_output = frequency_input.with_suffix(".out")
        frequency_check = frequency_stability_check(frequency_output)
        transition_state = {
            "seed_xyz": str(ts_seed),
            "optimization": ts_execution,
            "optimized_xyz": str(optimized_ts),
            "frequency_execution": frequency_execution,
            "frequency_output": str(frequency_output),
            "frequency_check": frequency_check,
        }
        record["transition_state"] = transition_state
        if frequency_check.get("NumImaginary") != 1:
            record.update(status="ts_frequency_validation_failed", path_classification="surface_unresolved")
            return record
        hessian = frequency_input.with_suffix(".hess")
        irc_input = create_orca_irc_input(
            optimized_ts, charge=route.charge, multiplicity=route.multiplicity,
            hessian_file=hessian if hessian.is_file() else None, label="irc",
            ncores=ncores, method_keywords=method_keywords,
        )
        irc_execution = _run_or_resume(irc_input, timeout_seconds=timeout_seconds)
        trajectory = _locate_irc_trajectory(ts_dir, irc_input.stem)
        if trajectory is None:
            raise RuntimeError("ORCA IRC completed without a retained bidirectional trajectory")
        retained = folder / "validated-trajectory.xyz"
        if trajectory.resolve() != retained.resolve():
            shutil.copy2(trajectory, retained)
        trajectory_validation = _validate_irc_route_endpoints(
            route, retained, reactant_xyz, product_xyz,
            mapping_result=mapping_result, separated_fragments=separated_fragments,
        )
        record["irc"] = {
            "execution": irc_execution,
            "trajectory": str(retained),
            "endpoint_validation": trajectory_validation,
        }
        if not trajectory_validation["valid"]:
            record.update(status="irc_endpoint_validation_failed", path_classification="surface_unresolved")
            return record
        record.update(
            status="verified_transition_state_path",
            path_classification="verified_barriered_path",
            trajectory=str(retained),
            path_evidence={
                "kind": "orca_ts_frequency_irc",
                "one_imaginary_frequency": True,
                "irc_matches_declared_endpoints": True,
                "trajectory": str(retained),
            },
        )
        return record
    except Exception as error:
        record.update(
            status="path_execution_failed", path_classification="surface_unresolved",
            failure_reason=str(error),
        )
        return record


def run_generic_reaction_path(
    route: RouteSpec,
    workdir: Path,
    *,
    prepare_only: bool = False,
    ncores: int = 1,
    method_keywords: list[str] | None = None,
    timeout_seconds: float | None = 14400.0,
    orientations: int = 3,
    neb_images: int = 8,
    barrierless_classifier: Callable[[dict], dict] | None = None,
    condition_contract: dict | None = None,
) -> dict:
    """Verify a resolved route without dispatching on a reaction family."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    mapping = resolve_route_atom_mapping(route)
    canonical: dict[str, Any] = {
        "status": "mapping_unresolved",
        "route_id": route.route_id,
        "reaction_equation": route.reaction_equation,
        "path_classification": "surface_unresolved",
        "endpoint_validation": {"status": "not_run", "valid": False},
        "mapping": mapping,
        "path_evidence": None,
        "rate_model": {"status": "rate_model_unavailable"},
        "propagation": {"status": "not_run"},
        "convergence": {"status": "not_converged"},
        "trajectory": None,
        "stationary_points": {"status": "not_run", "reactants": [], "products": [], "valid": False},
    }
    result: dict[str, Any] = {
        "schema_version": 1,
        "kind": "generic_orca_route_verification",
        "route": route.as_dict(),
        "route_verification": canonical,
        "orientations": [],
    }
    manifest = workdir / "route-verification.json"
    if not mapping.get("valid"):
        canonical["status"] = str(mapping.get("status") or "mapping_unresolved")
        canonical["convergence"] = {"status": "blocked_before_orca"}
        _write_json(manifest, result)
        return {**result, "path": manifest}
    selected_multiplicity = int(mapping["surface"]["selected_multiplicity"])
    execution_route = (
        route if route.multiplicity == selected_multiplicity
        else RouteSpec(**{**route.as_dict(), "multiplicity": selected_multiplicity})
    )
    canonical["selected_charge"] = int(mapping["surface"]["charge"])
    canonical["selected_multiplicity"] = selected_multiplicity
    separated_fragments = _route_uses_separated_fragment_endpoints(execution_route)
    canonical["endpoint_mode"] = (
        "separated_fragment_channels" if separated_fragments else "bound_endpoint_minima"
    )
    stationary: dict | None = None
    fragment_xyz = None
    if separated_fragments:
        # Each physical species is a real stationary point.  The assembled
        # collision/encounter geometry is only an atom-ordered NEB endpoint and
        # must not be required to optimize as a fictitious bound molecule.
        stationary = _stationary_fragment_records(
            execution_route, workdir / "stationary-points", mapping,
            execute=not prepare_only, ncores=ncores, method_keywords=method_keywords,
            timeout_seconds=timeout_seconds,
        )
        canonical["stationary_points"] = stationary
        if not prepare_only and not stationary["valid"]:
            canonical.update(
                status="stationary_point_validation_failed",
                endpoint_validation={"status": "blocked_before_endpoint_assembly", "valid": False},
                convergence={"status": "blocked_before_neb"},
            )
            _write_json(manifest, result)
            return {**result, "path": manifest}
        if not prepare_only:
            thermochemistry = assemble_route_thermochemistry(
                execution_route, stationary, condition_contract,
            )
            loss_bound = condition_specific_forward_loss_bound(
                execution_route, thermochemistry, condition_contract,
            )
            canonical["endpoint_thermochemistry"] = thermochemistry
            canonical["forward_loss_upper_bound"] = loss_bound
            if loss_bound.get("status") == "forward_loss_below_retention_threshold_upper_bound":
                canonical.update(
                    status="thermochemistry_bounded_below_loss_threshold",
                    path_classification="no_path_required_for_bounded_route",
                    endpoint_validation={
                        "status": "validated_stationary_species",
                        "valid": True,
                        "mode": "separated_fragment_channels",
                    },
                    convergence={"status": "not_required_after_conservative_loss_bound"},
                    rate_model={
                        "status": "forward_rate_upper_bound_only",
                        "rate_claim_allowed": False,
                    },
                )
                _write_json(manifest, result)
                return {**result, "path": manifest}
        fragment_xyz = _stationary_fragment_xyz(stationary, optimized=not prepare_only)
    seeds = build_endpoint_complex_seeds(
        execution_route, workdir / "endpoint-complexes", orientations=orientations,
        mapping_result=mapping, fragment_xyz=fragment_xyz,
    )
    result["endpoint_complexes"] = seeds
    if prepare_only:
        if stationary is None:
            stationary = _stationary_fragment_records(
                execution_route, workdir / "stationary-points", mapping, execute=False,
                ncores=ncores, method_keywords=method_keywords, timeout_seconds=timeout_seconds,
            )
        canonical.update(
            status="prepared", path_classification="not_computed",
            endpoint_validation={"status": "prepared", "valid": False},
            stationary_points=stationary,
            convergence={"status": "not_run"},
        )
        result["orientations"] = [
            {**record, "status": "prepared", "path_classification": "not_computed"}
            for record in seeds["orientations"]
        ]
        _write_json(manifest, result)
        return {**result, "path": manifest}

    for seed_record in seeds["orientations"]:
        print(
            f"[STORCA] Route {route.route_id}: ORCA path orientation "
            f"{seed_record['orientation']}/{seeds['orientation_count']}", flush=True,
        )
        orientation_result = _run_generic_orientation(
            execution_route, seed_record, mapping_result=mapping,
            ncores=ncores, method_keywords=method_keywords,
            timeout_seconds=timeout_seconds, neb_images=neb_images,
            barrierless_classifier=barrierless_classifier,
            separated_fragments=separated_fragments,
            stationary_points=stationary,
        )
        result["orientations"].append(orientation_result)
        _write_json(manifest, result)
        if orientation_result.get("status") == "verified_transition_state_path":
            break
        classified = [
            item for item in result["orientations"]
            if item.get("path_classification") in {
                "barrierless_capture_candidate", "barrierless_dissociation_candidate", "barrierless_path_candidate",
            }
        ]
        reproducible_pair_found = False
        for left_index, left in enumerate(classified):
            left_energies = (left.get("path_evidence") or {}).get("energies_hartree") or []
            for right in classified[left_index + 1:]:
                right_energies = (right.get("path_evidence") or {}).get("energies_hartree") or []
                if (left.get("path_classification") == right.get("path_classification")
                        and len(left_energies) >= 2 and len(right_energies) >= 2
                        and abs((left_energies[-1] - left_energies[0])
                                - (right_energies[-1] - right_energies[0])) <= 0.02):
                    reproducible_pair_found = True
                    break
            if reproducible_pair_found:
                break
        if reproducible_pair_found:
            break
    verified = [item for item in result["orientations"] if item.get("status") == "verified_transition_state_path"]
    barrierless = [
        item for item in result["orientations"]
        if item.get("path_classification") in {
            "barrierless_capture_candidate", "barrierless_dissociation_candidate", "barrierless_path_candidate",
        }
    ]
    reproducible_barrierless: list[dict] = []
    barrierless_reproducibility: dict[str, Any] = {
        "status": "not_applicable", "required_independent_orientations": 2,
    }
    if barrierless and not verified:
        grouped: dict[str, list[dict]] = {}
        for item in barrierless:
            grouped.setdefault(item["path_classification"], []).append(item)
        for classification, group in sorted(grouped.items()):
            for left_index, left in enumerate(group):
                left_profile = (left.get("path_evidence") or {}).get("energies_hartree") or []
                if len(left_profile) < 2:
                    continue
                left_energy = float(left_profile[-1]) - float(left_profile[0])
                for right in group[left_index + 1:]:
                    right_profile = (right.get("path_evidence") or {}).get("energies_hartree") or []
                    if len(right_profile) < 2:
                        continue
                    right_energy = float(right_profile[-1]) - float(right_profile[0])
                    if abs(left_energy - right_energy) <= 0.02:
                        reproducible_barrierless = [left, right]
                        barrierless_reproducibility = {
                            "status": "reproducible",
                            "classification": classification,
                            "orientation_count": 2,
                            "orientation_ids": [left["orientation"], right["orientation"]],
                            "reaction_energy_spread_hartree": abs(left_energy - right_energy),
                            "required_independent_orientations": 2,
                        }
                        break
                if reproducible_barrierless:
                    break
            if reproducible_barrierless:
                break
        if not reproducible_barrierless:
            barrierless_reproducibility = {
                "status": "surface_unresolved",
                "reason": "barrierless_profile_not_reproduced_in_two_independent_orientations",
                "classified_orientation_count": len(barrierless),
                "required_independent_orientations": 2,
            }
    canonical["barrierless_reproducibility"] = barrierless_reproducibility
    chosen = verified[0] if verified else reproducible_barrierless[0] if reproducible_barrierless else None
    if stationary is None:
        stationary = _stationary_fragment_records(
            execution_route, workdir / "stationary-points", mapping,
            execute=chosen is not None, ncores=ncores, method_keywords=method_keywords,
            timeout_seconds=timeout_seconds,
        )
    canonical["stationary_points"] = stationary
    if chosen is None:
        canonical.update(
            status="surface_unresolved",
            endpoint_validation={"status": "no_validated_orientation", "valid": False},
            convergence={"status": "path_unresolved"},
        )
    else:
        canonical.update(
            path_classification=chosen["path_classification"],
            endpoint_validation=chosen.get("endpoint_validation", {"valid": False}),
            path_evidence=chosen.get("path_evidence"),
            trajectory=chosen.get("trajectory"),
            selected_orientation=chosen["orientation"],
            transition_state=chosen.get("transition_state"),
            irc=chosen.get("irc"),
        )
        if not stationary["valid"]:
            canonical.update(
                status="stationary_point_validation_failed",
                convergence={"status": "blocked_before_kinetics"},
            )
        elif verified:
            canonical.update(
                status="verified_transition_state_path",
                convergence={"status": "path_verified_rate_pending"},
                arkane_inputs={
                    "reactants": stationary["reactants"],
                    "products": stationary["products"],
                    "transition_state": canonical["transition_state"],
                },
            )
        else:
            canonical.update(
                status="barrierless_rate_model_required",
                convergence={"status": "route_physics_classified_rate_pending"},
            )
    _write_json(manifest, result)
    return {**result, "path": manifest}


def _run_neb(route: RouteSpec, explanation_dir: Path, *, prepare_only: bool,
             ncores: int, method_keywords: list[str] | None,
             timeout_seconds: float | None, neb_images: int) -> dict:
    endpoint_dir = explanation_dir / "endpoints"
    seeds = build_mapped_endpoint_seeds(route, endpoint_dir)
    if prepare_only:
        return {"status": "prepared", "endpoints": seeds}
    reactant_result = _optimize_endpoint(
        Path(seeds["reactant_xyz"]), endpoint_dir / "reactant", charge=route.charge,
        multiplicity=route.multiplicity, ncores=ncores, method_keywords=method_keywords,
        timeout_seconds=timeout_seconds,
    )
    product_result = _optimize_endpoint(
        Path(seeds["product_xyz"]), endpoint_dir / "product", charge=route.charge,
        multiplicity=route.multiplicity, ncores=ncores, method_keywords=method_keywords,
        timeout_seconds=timeout_seconds,
    )
    if not reactant_result["frequency_check"]["IsMinimum"] or not product_result["frequency_check"]["IsMinimum"]:
        raise RuntimeError("Both NEB endpoints must be ORCA local minima")
    path_dir = explanation_dir / "path"
    path_dir.mkdir(parents=True, exist_ok=True)
    reactant_xyz = path_dir / "reactant.xyz"
    product_xyz = path_dir / "product.xyz"
    shutil.copy2(reactant_result["optimized_xyz"], reactant_xyz)
    shutil.copy2(product_result["optimized_xyz"], product_xyz)
    neb_input = create_orca_neb_ts_input(
        reactant_xyz, product_xyz, charge=route.charge, multiplicity=route.multiplicity,
        ncores=ncores, nimages=neb_images, method_keywords=method_keywords, preopt_ends=True,
    )
    neb_execution = _run_or_resume(neb_input, timeout_seconds=timeout_seconds)
    ts_xyz = _locate_neb_ts(path_dir, neb_input.stem)
    if ts_xyz is None:
        raise RuntimeError("ORCA NEB-TS completed without a retained TS geometry")
    retained_ts = path_dir / "transition-state.xyz"
    if ts_xyz.resolve() != retained_ts.resolve():
        shutil.copy2(ts_xyz, retained_ts)
    frequency_input = create_orca_input(
        retained_ts, route.charge, route.multiplicity, freq=True, label="ts-frequency",
        ncores=ncores, method_keywords=method_keywords,
    )
    frequency_execution = _run_or_resume(frequency_input, timeout_seconds=timeout_seconds)
    frequency_check = frequency_stability_check(frequency_input.with_suffix(".out"))
    if frequency_check["NumImaginary"] != 1:
        raise RuntimeError(
            f"Transition-state frequency check found {frequency_check['NumImaginary']} significant imaginary modes; exactly one is required"
        )
    hessian = frequency_input.with_suffix(".hess")
    irc_input = create_orca_irc_input(
        retained_ts, charge=route.charge, multiplicity=route.multiplicity,
        hessian_file=hessian if hessian.is_file() else None, ncores=ncores,
        method_keywords=method_keywords,
    )
    irc_execution = _run_or_resume(irc_input, timeout_seconds=timeout_seconds)
    trajectory = _locate_irc_trajectory(path_dir, irc_input.stem)
    if trajectory is None:
        raise RuntimeError("ORCA IRC completed without the combined forward/backward trajectory")
    retained_trajectory = path_dir / "trajectory.xyz"
    if trajectory.resolve() != retained_trajectory.resolve():
        shutil.copy2(trajectory, retained_trajectory)
    endpoint_validation = validate_trajectory_endpoints(retained_trajectory, reactant_xyz, product_xyz)
    if not endpoint_validation["valid"]:
        raise RuntimeError("The IRC endpoints do not match the declared reactant and product")
    return {
        "status": "verified_transition_state_path",
        "endpoints": {
            **seeds,
            "optimized_reactant_xyz": str(reactant_xyz),
            "optimized_product_xyz": str(product_xyz),
        },
        "neb_execution": neb_execution,
        "frequency_execution": frequency_execution,
        "frequency_check": frequency_check,
        "irc_execution": irc_execution,
        "trajectory": str(retained_trajectory),
        "endpoint_validation": endpoint_validation,
    }


def run_decomposition_explanation(
    report_path: Path,
    *,
    output_dir: Path | None = None,
    route_id: str | None = None,
    prepare_only: bool = False,
    charge: int | None = None,
    multiplicity: int | None = None,
    ncores: int = 1,
    method_keywords: list[str] | None = None,
    timeout_seconds: float | None = 14400.0,
    scan_steps: int = 20,
    neb_images: int = 8,
) -> dict:
    """Select, calculate, validate, and visualize one decomposition route."""
    report_path = Path(report_path)
    source_dir = report_path if report_path.is_dir() else report_path.parent
    explanation_dir = Path(output_dir) if output_dir else source_dir / "explanation"
    explanation_dir.mkdir(parents=True, exist_ok=True)
    prepared = prepare_explanation(report_path, route_id)
    route_data = prepared["selected_route"]
    if charge is not None:
        route_data["charge"] = charge
    if multiplicity is not None:
        route_data["multiplicity"] = multiplicity
    route = RouteSpec(**route_data)
    result: dict[str, Any] = {
        **prepared,
        "selected_route": route.as_dict(),
        "status": "running",
        "method_keywords": method_keywords or ["B3LYP", "def2-SVP", "NoAutoStart"],
        "resources": {
            "ncores": ncores,
            "timeout_seconds_per_orca_job": timeout_seconds,
            "scan_steps": scan_steps,
            "neb_images": neb_images,
        },
    }
    manifest = explanation_dir / "decomposition-explanation.json"
    _write_json(manifest, result)
    visual_dir = explanation_dir / "visuals"
    storyboard = render_candidate_storyboard(
        route.parent_smiles, list(route.product_smiles), visual_dir / "candidate-storyboard.png",
        title=f"{route.route_id} ({route.source})",
    ) if route.product_smiles else None
    mapping = validate_route_mapping(route)
    result["mapping_validation"] = mapping
    if route.protocol == "unsupported_atom_mapping" or not mapping["valid"]:
        result["status"] = "unsupported_route_physics"
        result["failure_reason"] = route.limitation or mapping["status"]
        if storyboard:
            result["visuals"] = {
                "storyboard": str(storyboard),
                "manifest": str(write_visual_manifest(
                    visual_dir / "visuals.json", route_id=route.route_id,
                    evidence_level="candidate_storyboard", storyboard=storyboard,
                )),
            }
        _write_json(manifest, result)
        return {**result, "result_json": manifest}
    try:
        if route.protocol == "relaxed_dissociation_scan":
            path_result = _run_dissociation(
                route, explanation_dir, prepare_only=prepare_only, ncores=ncores,
                method_keywords=method_keywords, timeout_seconds=timeout_seconds, scan_steps=scan_steps,
            )
        elif route.protocol in {"neb_ts_irc", "generic_endpoint_path"}:
            generic = run_generic_reaction_path(
                route, explanation_dir / "generic-path", prepare_only=prepare_only,
                ncores=ncores, method_keywords=method_keywords,
                timeout_seconds=timeout_seconds, orientations=3, neb_images=neb_images,
            )
            canonical = generic["route_verification"]
            first_seed = generic["endpoint_complexes"]["orientations"][0]
            path_result = {
                "status": canonical["status"],
                "endpoints": {
                    "reactant_xyz": first_seed["reactant_xyz"],
                    "product_xyz": first_seed["product_xyz"],
                },
                "trajectory": canonical.get("trajectory"),
                "generic_verification": generic,
            }
            if not prepare_only and not path_result["trajectory"]:
                raise RuntimeError(
                    f"Generic route verification stopped at {canonical['status']} without a validated trajectory"
                )
        else:
            raise ValueError(f"Unsupported route protocol: {route.protocol}")
        result["path_calculation"] = path_result
        if prepare_only:
            preview = _render_preview(route, path_result["endpoints"], visual_dir, storyboard)
            result.update(status="prepared", visuals=preview)
        else:
            trajectory = Path(path_result["trajectory"])
            evidence = path_result["status"]
            animation = render_xyz_trajectory_gif(
                trajectory, visual_dir / "decomposition.gif",
                title=f"{route.route_id}: decomposition path", evidence_label=evidence.replace("_", " ").upper(),
                energies=path_result.get("energies_hartree"),
            )
            energy_profile = None
            if path_result.get("energies_hartree"):
                energy_profile = render_energy_profile(
                    path_result["energies_hartree"], visual_dir / "energy-profile.png"
                )
            visual_manifest = write_visual_manifest(
                visual_dir / "visuals.json", route_id=route.route_id, evidence_level=evidence,
                storyboard=storyboard, animation=animation, energy_profile=energy_profile,
                trajectory=trajectory,
            )
            result.update(
                status=evidence,
                visuals={
                    "storyboard": str(storyboard) if storyboard else None,
                    "animation": str(animation),
                    "energy_profile": str(energy_profile) if energy_profile else None,
                    "trajectory": str(trajectory),
                    "manifest": str(visual_manifest),
                },
            )
    except Exception as error:
        result["status"] = "failed"
        result["failure_reason"] = str(error)
        # A failed computation may still provide a clearly labelled preview if
        # endpoint generation completed before the failure.
        endpoints = (result.get("path_calculation") or {}).get("endpoints")
        if storyboard and endpoints:
            result["visuals"] = _render_preview(route, endpoints, visual_dir, storyboard)
    _write_json(manifest, result)
    return {**result, "result_json": manifest}
