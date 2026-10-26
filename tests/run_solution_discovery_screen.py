"""One preregistered native screening diagnostic; never computes rates/lifetimes."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from unittest.mock import patch

from storca.graph_xtb_discovery import GraphXtbBudget, GraphXtbDiscoveryBackend, run_xtb_endpoint_screen
from storca.reaction_candidates import canonical_smiles


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inventory(root):
    return {str(path.relative_to(root)): {"sha256": digest(path), "mtime_ns": path.stat().st_mtime_ns}
            for path in root.rglob("*") if path.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    if root.exists():
        raise ValueError("Choose a new diagnostic output directory; prior evidence is immutable")
    root.mkdir(parents=True)
    budget = GraphXtbBudget(max_generated_graph_candidates=8, max_structurally_valid_candidates=1,
                           max_products_per_reactant_state=1, max_xtb_optimizations=2, max_xtb_failures=1,
                           max_candidates_for_verification=1, timeout_seconds=15.0, overall_timeout_seconds=30.0,
                           allow_dissociation=False, allow_bond_formation=False, allow_bond_rearrangement=False,
                           allow_hydrogen_transfer=False, allow_bimolecular_radical_association=False,
                           allow_bimolecular_hydrogen_abstraction=False)
    cases = (("ring-water", "C1=CCC1", "water", "C=CC=C"),
             ("diene-water", "C=CC=C", "water", "C1=CCC1"),
             ("ring-acetonitrile", "C1=CCC1", "acetonitrile", "C=CC=C"))
    plan = {"kind": "native_solution_discovery_screening_plan", "cases": cases,
            "maximum_endpoint_optimizations": 6, "aggregate_cores": 1,
            "batch_hard_timeout_seconds": 90, "budget_per_case": budget.as_dict(),
            "prediction_proved": False, "orca_jobs_allowed": 0, "propagation_allowed": False,
            "reference_isolation": "Expected products are diagnostic comparisons only; never supplied to discovery"}
    (root / "plan.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    started = time.monotonic()
    records, native_jobs = [], 0
    for case_id, smiles, solvent, expected_product in cases:
        if time.monotonic() - started >= 90:
            break
        environment = {"phase": "solution", "solvent": solvent, "solvent_model": "SMD", "temperature_K": 300.0}
        result = GraphXtbDiscoveryBackend(ncores=1).discover(
            {"smiles": smiles, "charge": 0, "multiplicity": 1, "environment": environment},
            {"contract_id": "development-screen:" + case_id}, budget, workdir=root / case_id)
        record = {"case_id": case_id, "status": result.status, "screening_environment": result.provenance["screening_environment"],
                  "candidate_count": len(result.candidates), "expected_product_recovered": False,
                  "native_endpoint_jobs": result.provenance["xtb_endpoint_optimizations"], "cache_replay": None}
        native_jobs += record["native_endpoint_jobs"]
        if len(result.candidates) == 1:
            candidate = result.candidates[0]
            record["candidate_id"] = candidate.candidate_id
            record["screening_status"] = candidate.screening_status
            record["expected_product_recovered"] = tuple(map(canonical_smiles, candidate.products)) == (canonical_smiles(expected_product),)
            screen = candidate.provenance.get("xtb_screen") or {}
            record["estimated_reaction_screen_energy_kcal_mol"] = candidate.estimated_reaction_energy_kcal_mol
            if candidate.screening_status == "candidate_unverified":
                endpoint_root = root / case_id / candidate.candidate_id
                before = inventory(endpoint_root)
                forbidden, calls = [], []
                native_run = subprocess.run
                def version_only(command, **kwargs):
                    if command[1:] != ["--version"]:
                        forbidden.append(command)
                        raise AssertionError("Cache replay attempted new chemistry")
                    calls.append(command)
                    return native_run(command, **kwargs)
                reused = []
                with patch("storca.graph_xtb_discovery.subprocess.run", version_only):
                    for endpoint in screen["endpoint_records"]:
                        folder = Path(endpoint["artifacts"]["input_xyz"]).parent
                        replay = run_xtb_endpoint_screen(endpoint["smiles"], folder,
                                                         charge=0, multiplicity=1, timeout_seconds=15.0,
                                                         environment=environment, ncores=1)
                        reused.append(replay["status"] == "completed" and replay["energy_hartree"] == endpoint["energy_hartree"])
                record["cache_replay"] = {"all_endpoints_reused": all(reused), "native_job_attempts": len(forbidden),
                                           "version_reads": len(calls), "files_unchanged": before == inventory(endpoint_root),
                                           "file_count": len(before)}
        records.append(record)
    result = {"kind": "native_solution_discovery_screening_diagnostic", "plan_sha256": digest(root / "plan.json"),
              "elapsed_seconds": time.monotonic() - started, "aggregate_cores": 1, "native_endpoint_jobs": native_jobs,
              "cases": records, "completed_cases": len(records), "orca_jobs": 0, "rates_generated": 0,
              "prediction_proved": False, "uncertainty_calibration_proved": False, "mechanism_completeness_proved": False,
              "implementation_hashes": {str(path): digest(path) for path in
                                        (Path("storca/graph_xtb_discovery.py"), Path(__file__))}}
    (root / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
