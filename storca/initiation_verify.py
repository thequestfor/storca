"""Prepare and execute endpoint evidence for a dependency-selected initiation route."""

from __future__ import annotations

import json
from pathlib import Path

from src.molecule_tools import smiles_to_xyz

from .collision_repair import _allowed_total_multiplicities, _thermochemistry
from .workflow import run_optimization_and_frequency


def prepare_initiation_verification(run_dir: Path, route_index: int, *, execute_endpoints: bool = True,
                                    ncores: int = 1, method_keywords: list[str] | None = None) -> dict:
    run_dir = Path(run_dir)
    screen = json.loads((run_dir / "stability.json").read_text())
    route = screen["rmg_evidence"]["candidate_routes"][route_index]
    endpoints = route.get("resolved_endpoints") or {}
    reactants, products = endpoints.get("reactants", []), endpoints.get("products", [])
    if not reactants or not products or not all(item.get("smiles") and item.get("multiplicity") for item in reactants + products):
        raise ValueError("Selected initiation route does not have resolved endpoint structures and multiplicities")
    reactant_surfaces = _allowed_total_multiplicities(reactants)
    product_surfaces = _allowed_total_multiplicities(products)
    common = sorted(set(reactant_surfaces) & set(product_surfaces))
    folder = run_dir / "initiation-verification"
    folder.mkdir(exist_ok=True)
    jobs = []
    for side, records in (("reactant", reactants), ("product", products)):
        for number, endpoint in enumerate(records, 1):
            job_dir = folder / f"{side}-{number:02d}"
            job_dir.mkdir(exist_ok=True)
            xyz = job_dir / "input.xyz"
            smiles_to_xyz(endpoint["smiles"], xyz)
            record = {"side": side, "index": number, "endpoint": endpoint, "input_xyz": str(xyz)}
            if execute_endpoints:
                print(f"[STORCA] Upstream initiation endpoint {side}-{number:02d}: ORCA optimization/frequency started", flush=True)
                result = run_optimization_and_frequency(xyz, job_dir, charge=0, multiplicity=int(endpoint["multiplicity"]),
                                                        ncores=ncores, method_keywords=method_keywords)
                check = result["frequency_check"]
                record.update(status="completed" if check.get("IsMinimum") else "not_a_local_minimum",
                              optimized_xyz=str(result["optimized_xyz"]), frequency_check=check,
                              thermochemistry=_thermochemistry(result["frequency"]["out"]))
            else:
                record["status"] = "prepared"
            jobs.append(record)
    minima = execute_endpoints and all(job["status"] == "completed" for job in jobs)
    dossier = {
        "schema_version": 1, "kind": "orca_upstream_initiation_verification",
        "status": "requires_h_transfer_ts_verification" if minima else ("endpoint_jobs_prepared" if not execute_endpoints else "endpoint_validation_failed"),
        "selected_route_index": route_index, "reaction_equation": route.get("reaction_equation"),
        "route": route, "endpoint_jobs": jobs,
        "spin_surfaces": {"reactant_total_multiplicities": reactant_surfaces,
                          "product_total_multiplicities": product_surfaces,
                          "common_total_multiplicities": common,
                          "status": "surface_available" if common else "spin_crossing_required"},
        "downstream_status": "blocked_until_initiation_rate_verified",
        "rate_status": "not_calculated", "t95_status": "not_calculated",
        "next_step": "Generate H-transfer TS guesses on each common spin surface; require one relevant imaginary mode and bidirectional IRC.",
    }
    path = folder / "initiation-verification.json"
    path.write_text(json.dumps(dossier, indent=2, sort_keys=True, default=str) + "\n")
    return {**dossier, "path": path}
