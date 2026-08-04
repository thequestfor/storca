"""Repair preparation for RMG rates that exceed the gas-collision limit.

These routes cannot be accepted as kinetics.  This module turns the retained
RMG artifact into auditable ORCA endpoint jobs, while refusing to pretend that
a barrierless radical association has a transition state.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from src.inputgen import create_orca_input, create_orca_relaxed_scan_input
from src.molecule_tools import smiles_to_xyz
from src.orca_runner import run_orca
from src.parser import parse_orca_energy

from .rmg_execution import parse_collision_rate_violators
from .route_geometry import read_xyz, write_xyz
from .workflow import run_optimization_and_frequency


def _thermochemistry(output: Path) -> dict:
    """Retain ORCA's reported thermochemistry verbatim enough for audit."""
    text = Path(output).read_text(errors="replace")
    lines = [line.strip() for line in text.splitlines()
             if any(term in line.upper() for term in ("THERMAL ENTHALPY", "GIBBS FREE ENERGY", "ZERO POINT ENERGY"))]
    return {"electronic_energy_hartree": parse_orca_energy(output), "reported_thermal_lines": lines,
            "source_output": str(output)}


def _violators(path: Path) -> list[dict]:
    """Read the compact RMG collision report without depending on formatting."""
    text = Path(path).read_text(errors="replace") if Path(path).is_file() else ""
    blocks = re.split(r"(?=^! Template reaction:)", text, flags=re.MULTILINE)
    result = []
    for block in blocks:
        match = re.search(r"^([^!\n][^\n]*?(?:<=>|=>)[^\n]+?)\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s*$", block, re.MULTILINE)
        factor = re.search(r"Violation factor:\s*([0-9.Ee+-]+)", block)
        template = re.search(r"Template reaction:\s*([^\n]+)", block)
        if match and factor:
            result.append({"reaction_equation": match.group(1).strip(), "family": template.group(1).strip() if template else None,
                           "violation_factor": float(factor.group(1))})
    return result


def _endpoint_records(equation: str, species: list[dict]) -> tuple[list[dict], list[dict]]:
    by_label = {str(item.get("label")): item for item in species}
    by_label.update({str(item.get("chemkin_identifier")): item for item in species if item.get("chemkin_identifier")})
    arrow = "<=>" if "<=>" in equation else "=>"
    def parse(side: str) -> list[dict]:
        records = []
        for token in side.split("+"):
            label = re.sub(r"^\s*\d+(?:\.\d+)?\s*", "", token).strip()
            item = by_label.get(label, {})
            records.append({"label": label, "smiles": item.get("smiles"), "multiplicity": item.get("multiplicity"), "charge": 0})
        return records
    left, right = equation.split(arrow, 1)
    return parse(left), parse(right)


def _allowed_total_multiplicities(endpoints: list[dict]) -> list[int]:
    """Couple retained endpoint spins; no arbitrary spin surface is assumed."""
    spins = {0.0}
    for endpoint in endpoints:
        multiplicity = endpoint.get("multiplicity")
        if not multiplicity:
            return []
        spin = (int(multiplicity) - 1) / 2.0
        spins = {value for prior in spins for value in
                 [abs(prior - spin) + step for step in range(int(2 * min(prior, spin)) + 1)]}
    return sorted({int(round(2 * spin + 1)) for spin in spins})


def prepare_collision_rate_repairs(run_dir: Path, *, execute_endpoints: bool = False,
                                   ncores: int = 1, method_keywords: list[str] | None = None) -> dict:
    """Create, and optionally execute, ORCA endpoint checks for every violator.

    ``execute_endpoints`` intentionally does not claim a TS or rate.  It only
    establishes spin-labelled stationary-point evidence needed before a scan,
    TS/IRC, and Arkane step can be selected.
    """
    run_dir = Path(run_dir)
    screen_path = run_dir / "stability.json"
    screen = json.loads(screen_path.read_text())
    rmg = screen.get("rmg_evidence", {})
    collision = Path((rmg.get("artifacts") or {}).get("collision_rate_violators", run_dir / "rmg" / "collision_rate_violators.log"))
    validation = parse_collision_rate_violators(collision)
    if validation["status"] != "kinetics_unreliable":
        raise ValueError("This screen has no collision-rate-limit violator to repair.")
    species = (rmg.get("mechanism_inspection") or {}).get("species", [])
    output = run_dir / "collision-repair"
    output.mkdir(exist_ok=True)
    repairs = []
    for index, item in enumerate(_violators(collision), 1):
        print(f"[STORCA] Collision-route {index}: resolving ORCA endpoints", flush=True)
        reactants, products = _endpoint_records(item["reaction_equation"], species)
        resolved = all(x["smiles"] and x["multiplicity"] for x in reactants + products)
        barrierless = item.get("family") in {"R_Recombination", "R_Addition_MultipleBond"}
        reactant_surfaces = _allowed_total_multiplicities(reactants) if resolved else []
        product_surfaces = _allowed_total_multiplicities(products) if resolved else []
        common_surfaces = sorted(set(reactant_surfaces) & set(product_surfaces))
        folder = output / f"route-{index:02d}"
        folder.mkdir(exist_ok=True)
        endpoint_jobs = []
        if resolved:
            for side, endpoints in (("reactant", reactants), ("product", products)):
                for endpoint_index, endpoint in enumerate(endpoints, 1):
                    name = f"{side}-{endpoint_index:02d}"
                    job_dir = folder / name
                    job_dir.mkdir(exist_ok=True)
                    xyz = job_dir / "input.xyz"
                    smiles_to_xyz(endpoint["smiles"], xyz)
                    if execute_endpoints:
                        print(f"[STORCA] ORCA endpoint {name}: optimization/frequency started", flush=True)
                        calculation = run_optimization_and_frequency(xyz, job_dir, charge=endpoint["charge"],
                                                                     multiplicity=endpoint["multiplicity"], ncores=ncores,
                                                                     method_keywords=method_keywords)
                        print(f"[STORCA] ORCA endpoint {name}: local-minimum check completed", flush=True)
                        endpoint_jobs.append({"name": name, "endpoint": endpoint, "status": "completed",
                                              "optimized_xyz": str(calculation["optimized_xyz"]),
                                              "frequency_check": calculation.get("frequency_check"),
                                              "thermochemistry": _thermochemistry(calculation["frequency"]["out"])})
                    else:
                        endpoint_jobs.append({"name": name, "endpoint": endpoint, "status": "prepared",
                                              "input_xyz": str(xyz), "orca_input": str(create_orca_input(
                                                  xyz, charge=endpoint["charge"], multiplicity=endpoint["multiplicity"],
                                                  opt=True, freq=True, label="endpoint", ncores=ncores, method_keywords=method_keywords))})
        repair = {
            "route_id": f"collision-route-{index}", **item, "reactants": reactants, "products": products,
            "endpoint_structures_resolved": resolved, "reaction_class": "barrierless_association_candidate" if barrierless else "barriered_or_unknown",
            "spin_surfaces": {"reactant_total_multiplicities": reactant_surfaces,
                              "product_total_multiplicities": product_surfaces,
                              "common_total_multiplicities": common_surfaces,
                              "status": "surface_available" if common_surfaces else "spin_crossing_or_unresolved"},
            "next_orca_protocol": ("optimize endpoints; perform spin-surface approach scan; only seek a TS if scan shows a barrier"
                                    if barrierless else "optimize endpoints; locate TS; require one imaginary mode and IRC"),
            "endpoint_jobs": endpoint_jobs,
            "rate_status": "unverified", "can_repair_lifetime": False,
        }
        repairs.append(repair)
        (folder / "repair-route.json").write_text(json.dumps(repair, indent=2, sort_keys=True) + "\n")
    dossier = {"schema_version": 1, "kind": "collision_rate_repair", "status": "endpoint_jobs_completed" if execute_endpoints else "prepared_for_orca_endpoint_checks",
               "source_stability_json": str(screen_path), "collision_validation": validation,
               "limitations": ["Endpoint calculations alone do not verify a reaction rate.", "No t95 may be reported until a route-specific rate passes physical validation and is propagated again."],
               "repairs": repairs}
    path = output / "collision-repair.json"
    path.write_text(json.dumps(dossier, indent=2, sort_keys=True) + "\n")
    return {**dossier, "path": path}


def _radical_atom_index(smiles: str) -> int:
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid radical endpoint SMILES: {smiles}")
    radicals = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetNumRadicalElectrons()]
    if len(radicals) != 1:
        raise ValueError("Approach scan requires exactly one retained radical centre")
    return radicals[0]


def _association_seed(oxygen_smiles: str, radical_smiles: str, output: Path, *, orientation: int) -> tuple[Path, tuple[int, int]]:
    """Place O2 near the radical centre; this is an approach seed, never a TS."""
    import numpy as np
    output = Path(output)
    oxygen_xyz, radical_xyz = output.with_name("oxygen.xyz"), output.with_name("radical.xyz")
    smiles_to_xyz(oxygen_smiles, oxygen_xyz)
    smiles_to_xyz(radical_smiles, radical_xyz)
    oxygen_elements, oxygen_coordinates, _ = read_xyz(oxygen_xyz)
    radical_elements, radical_coordinates, _ = read_xyz(radical_xyz)
    radical_index = _radical_atom_index(radical_smiles)
    # Orient the O--O vector perpendicular to the approach.  Two opposing
    # placements protect against a single arbitrary initial orientation.
    sign = 1.0 if orientation == 1 else -1.0
    oxygen_coordinates -= oxygen_coordinates[0]
    oxygen_coordinates[:, 0] += radical_coordinates[radical_index, 0] + 4.5
    oxygen_coordinates[:, 1] += radical_coordinates[radical_index, 1] + sign * 0.8
    oxygen_coordinates[:, 2] += radical_coordinates[radical_index, 2]
    write_xyz(output, oxygen_elements + radical_elements,
              np.vstack((oxygen_coordinates, radical_coordinates)),
              "O2-radical approach seed; not a transition state")
    return output, (0, len(oxygen_elements) + radical_index)


def _scan_energies(path: Path) -> list[float] | None:
    candidates = list(Path(path).glob("*.relaxscanact.dat")) + list(Path(path).glob("*_relaxscanact.dat"))
    if not candidates:
        return None
    energies = []
    for line in candidates[0].read_text(errors="replace").splitlines():
        values = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", line)
        if len(values) >= 2:
            energies.append(float(values[-1]))
    return energies or None


def prepare_association_scans(run_dir: Path, *, execute: bool = False, ncores: int = 1,
                              steps: int = 16, method_keywords: list[str] | None = None) -> dict:
    """Prepare/execute doublet association scans for repaired barrierless routes."""
    from .association_scan import run_resumable_association_scans
    return run_resumable_association_scans(
        run_dir, execute=execute, ncores=ncores, points=steps,
        method_keywords=method_keywords,
    )

    # Legacy monolithic FullScan implementation retained below temporarily for
    # source-history readability; the return above intentionally supersedes it.
    run_dir = Path(run_dir)
    repair_path = run_dir / "collision-repair" / "collision-repair.json"
    dossier = json.loads(repair_path.read_text())
    scans = []
    for repair in dossier.get("repairs", []):
        if repair.get("reaction_class") != "barrierless_association_candidate":
            continue
        if execute:
            endpoint_jobs = repair.get("endpoint_jobs", [])
            endpoint_minima = all(
                job.get("status") == "completed"
                and (job.get("frequency_check") or {}).get("IsMinimum") is True
                for job in endpoint_jobs
            )
            if not endpoint_jobs or not endpoint_minima:
                raise ValueError(
                    "Run stability-repair-collision --execute-endpoints first and require every endpoint to be an ORCA local minimum before executing an association scan."
                )
        if 2 not in (repair.get("spin_surfaces") or {}).get("common_total_multiplicities", []):
            scans.append({"route_id": repair["route_id"], "status": "spin_surface_unresolved"})
            continue
        oxygen = next((x for x in repair["reactants"] if x.get("smiles") == "[O][O]"), None)
        radical = next((x for x in repair["reactants"] if x.get("multiplicity") == 2), None)
        if not oxygen or not radical:
            scans.append({"route_id": repair["route_id"], "status": "unsupported_approach_endpoints"})
            continue
        for orientation in (1, 2):
            folder = run_dir / "collision-repair" / repair["route_id"] / f"approach-{orientation}"
            folder.mkdir(parents=True, exist_ok=True)
            seed, bond = _association_seed(oxygen["smiles"], radical["smiles"], folder / "approach.xyz", orientation=orientation)
            scan_input = create_orca_relaxed_scan_input(seed, bond_atom_indices=bond,
                start_distance_angstrom=4.5, end_distance_angstrom=1.35, charge=0, multiplicity=2,
                label="doublet-association-scan", ncores=ncores, steps=steps, method_keywords=method_keywords)
            record = {"route_id": repair["route_id"], "orientation": orientation, "multiplicity": 2,
                      "coordinate_atom_indices": list(bond), "input": str(scan_input), "status": "prepared",
                      "interpretation": "A constrained approach scan tests for a ground-state barrier. It is not a TS, IRC, rate, or lifetime."}
            if execute:
                output = run_orca(scan_input)["out"]
                normal = "ORCA TERMINATED NORMALLY" in Path(output).read_text(errors="replace")
                energies = _scan_energies(folder)
                record.update(status="completed" if normal and energies else "incomplete", output=str(output), energies_hartree=energies)
                if normal and energies and max(energies[1:-1] or energies) > max(energies[0], energies[-1]):
                    record["classification"] = "barriered_ts_candidate"
                elif normal and energies:
                    record["classification"] = "barrierless_capture_candidate"
            scans.append(record)
    result = {"schema_version": 1, "kind": "collision_association_scan", "status": "computed" if execute else "prepared",
              "source_repair": str(repair_path), "scans": scans,
              "limitations": ["A scan must be complete and physically coherent before selecting a TS workflow.", "A barrierless result requires a capture-compatible rate model; it is not an Arkane TST result."]}
    path = run_dir / "collision-repair" / "association-scans.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {**result, "path": path}


def finalize_association_verification(run_dir: Path) -> dict:
    """Turn completed ORCA scans into a precise route-physics conclusion.

    This deliberately does *not* convert barrierless capture into a t95.  The
    conclusion is limited to the scanned, declared charge and spin surface.
    """
    run_dir = Path(run_dir)
    scan_path = run_dir / "collision-repair" / "association-scans.json"
    scan_report = json.loads(scan_path.read_text())
    scans = scan_report.get("scans", [])
    if scan_report.get("orientations"):
        scans = []
        for orientation in scan_report["orientations"]:
            points = [point for point in orientation.get("points", []) if point.get("status") == "completed"]
            energies = [point.get("energy_hartree") for point in points]
            complete = len(points) >= 5 and orientation.get("coverage_angstrom", 0.0) >= 2.0 and all(value is not None for value in energies)
            classification = None
            if complete:
                interior = energies[1:-1]
                classification = ("barriered_ts_candidate" if interior and max(interior) > max(energies[0], energies[-1])
                                  else "barrierless_capture_candidate")
            scans.append({"orientation": orientation.get("orientation"), "status": "completed" if complete else "incomplete",
                          "classification": classification, "energies_hartree": energies})
    completed = [scan for scan in scans if scan.get("status") == "completed"]
    classifications = {scan.get("classification") for scan in completed}
    if not completed:
        status = "association_scan_not_completed"
        reason = "No ORCA association scan has completed normally with a retained energy profile."
    elif "barriered_ts_candidate" in classifications:
        status = "barriered_route_requires_ts_irc"
        reason = "At least one orientation has an interior scan maximum; a TS optimization, one imaginary mode, and IRC are required."
    elif classifications == {"barrierless_capture_candidate"}:
        status = "barrierless_on_scanned_doublet_surface"
        reason = "Every completed doublet approach scan is downhill or lacks an interior maximum on the declared coordinate. No saddle point is claimed or required by this result."
    else:
        status = "association_surface_unresolved"
        reason = "The completed scan orientations do not yield one consistent route classification."
    result = {
        "schema_version": 1,
        "kind": "orca_association_route_verification",
        "status": status,
        "reason": reason,
        "source_scans": str(scan_path),
        "completed_orientation_count": len(completed),
        "route_classifications": sorted(value for value in classifications if value),
        "rate_status": "not_calculated",
        "t95_status": "not_calculated",
        "next_step": (
            "Use a collision-capped capture model only if a quantitative rate is required."
            if status == "barrierless_on_scanned_doublet_surface" else
            "Optimize the scan maximum as a TS, then require frequency and IRC validation."
            if status == "barriered_route_requires_ts_irc" else
            "Complete or repair the ORCA scan evidence before selecting kinetics."
        ),
        "limitations": [
            "This conclusion applies only to the scanned doublet ground-state coordinate and sampled orientations.",
            "Barrierless association is route-physics evidence, not a kinetic lifetime or a proof of irreversible decomposition.",
        ],
    }
    path = run_dir / "collision-repair" / "association-verification.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {**result, "path": path}


def write_capture_rate_bounds(run_dir: Path, *, temperature_K: float = 298.15,
                              collision_diameter_angstrom: float = 4.0) -> dict:
    """Write an explicitly conservative capture-rate bound for a verified scan.

    This is not TST or an Arkane result.  Its only quantitative claim is a
    hard-sphere gas-collision ceiling; the lower bound is deliberately zero
    until a capture/recrossing model is supplied.
    """
    import math
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    run_dir = Path(run_dir)
    verified = json.loads((run_dir / "collision-repair" / "association-verification.json").read_text())
    if verified.get("status") != "barrierless_on_scanned_doublet_surface":
        raise ValueError("Capture bounds require a completed barrierless doublet-surface classification.")
    repairs = json.loads((run_dir / "collision-repair" / "collision-repair.json").read_text()).get("repairs", [])
    repair = next((x for x in repairs if x.get("reaction_class") == "barrierless_association_candidate"), None)
    if not repair:
        raise ValueError("No barrierless collision-repair route is available.")
    reactants = repair["reactants"]
    masses = [Descriptors.MolWt(Chem.MolFromSmiles(item["smiles"])) for item in reactants]
    if len(masses) != 2:
        raise ValueError("The conservative capture model currently supports bimolecular association only.")
    amu_kg = 1.66053906660e-27
    reduced_mass = masses[0] * masses[1] / sum(masses) * amu_kg
    boltzmann = 1.380649e-23
    diameter_m = collision_diameter_angstrom * 1e-10
    # m3/s -> cm3 molecule-1 s-1
    ceiling = math.pi * diameter_m ** 2 * math.sqrt(8 * boltzmann * temperature_K / (math.pi * reduced_mass)) * 1e6
    result = {
        "schema_version": 1, "kind": "conservative_barrierless_capture_bounds",
        "status": "capture_rate_range_not_sufficient_for_t95",
        "temperature_K": temperature_K, "reactant_molecular_weights_g_mol": masses,
        "assumed_collision_diameter_angstrom": collision_diameter_angstrom,
        "rate_constant_cm3_molecule_s": {"lower_bound": 0.0, "upper_collision_bound": ceiling},
        "assumptions": ["Hard-sphere collision upper bound.", "Unknown capture probability and recrossing are represented by a zero lower bound."],
        "prohibited_uses": ["Do not replace an RMG rate with this interval.", "Do not calculate t95 from this interval."],
    }
    path = run_dir / "collision-repair" / "capture-rate-bounds.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {**result, "path": path}


def run_collision_repair_pipeline(run_dir: Path, *, ncores: int = 1, scan_steps: int = 16,
                                  method_keywords: list[str] | None = None) -> dict:
    """Execute the available ORCA evidence stages, never fabricating a rate."""
    endpoints = prepare_collision_rate_repairs(run_dir, execute_endpoints=True, ncores=ncores, method_keywords=method_keywords)
    scans = prepare_association_scans(run_dir, execute=True, ncores=ncores, steps=scan_steps, method_keywords=method_keywords)
    verification = finalize_association_verification(run_dir)
    capture = None
    if verification["status"] == "barrierless_on_scanned_doublet_surface":
        capture = write_capture_rate_bounds(run_dir)
    status = (
        "orca_supported_barrierless_air_reactivity_lifetime_unverified"
        if verification["status"] == "barrierless_on_scanned_doublet_surface" else
        "orca_barriered_route_requires_ts_irc"
        if verification["status"] == "barriered_route_requires_ts_irc" else
        verification["status"]
    )
    result = {"schema_version": 1, "kind": "collision_repair_pipeline", "status": status,
              "endpoint_dossier": str(endpoints["path"]), "scan_dossier": str(scans["path"]),
              "verification_dossier": str(verification["path"]), "capture_bounds": str(capture["path"]) if capture else None,
              "assessment": "No final rate, repaired mechanism, flux attribution, or t95 is issued until a validated rate model is available."}
    path = Path(run_dir) / "collision-repair" / "pipeline.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {**result, "path": path}
