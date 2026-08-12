"""Resumable, coverage-balanced extension of static xTB environment ensembles."""

from __future__ import annotations

import json
from pathlib import Path

from .frequency_transfer import build_frequency_transfer_artifacts
from .xtb_frequencies import sample_xtb_snapshot_frequencies
from .xtb_sampling import sample_xtb_dimer_environments


EXTENSION_SCHEMA_VERSION = 1


def extension_targets(initial: int, maximum: int, batch_size: int = 20) -> list[int]:
    if initial < 0 or maximum < initial or batch_size < 1:
        raise ValueError("Invalid adaptive xTB extension limits")
    targets = list(range(initial + batch_size, maximum + 1, batch_size))
    if maximum > initial and (not targets or targets[-1] != maximum):
        targets.append(maximum)
    return targets


def _retained_sampling_records(run_dir: Path) -> list[dict]:
    path = Path(run_dir) / "environment-sampling.json"
    if not path.is_file():
        return []
    payload = json.loads(path.read_text())
    return list((payload.get("sampling") or {}).get("candidates") or [])


def extend_xtb_environment_ensemble(
    run_dir: Path, *, smiles: str, monomer_xyz: Path, maximum_candidates: int = 100,
    batch_size: int = 20, charge: int = 0, multiplicity: int = 1,
    ncores: int = 1, include_trimers: bool = True, progress=None,
) -> dict:
    """Extend in real acquisition rounds and stop after two stable comparisons."""
    run_dir = Path(run_dir)
    existing = _retained_sampling_records(run_dir)
    if not existing:
        raise RuntimeError("Adaptive extension requires a retained initial xTB ensemble")
    retained_plan_path = run_dir / "xtb-extension-plan.json"
    retained_plan = (
        json.loads(retained_plan_path.read_text()) if retained_plan_path.is_file() else {}
    )
    initial = int(retained_plan.get("initial_candidates", len(existing)))
    all_targets = extension_targets(initial, maximum_candidates, batch_size)
    targets = [target for target in all_targets if target > len(existing)]
    plan = {
        "schema_version": EXTENSION_SCHEMA_VERSION,
        "kind": "adaptive_xtb_environment_extension",
        "initial_candidates": initial,
        "maximum_candidates": maximum_candidates,
        "batch_size": batch_size,
        "cumulative_targets": all_targets,
        "round_model": "independently_seeded_balanced_dimer_trimer_strata",
        "population_model": "stratified_equal_weights_not_liquid_occupancies",
    }
    plan_path = retained_plan_path
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    retained_report_path = run_dir / "xtb-extension.json"
    retained_report = (
        json.loads(retained_report_path.read_text())
        if retained_report_path.is_file() else {}
    )
    rounds = list(retained_report.get("rounds") or [])
    for target in targets:
        round_index = all_targets.index(target) + 1
        records = _retained_sampling_records(run_dir)
        if len(records) < target:
            requested = target - len(records)
            if progress:
                progress(
                    f"Extending xTB environment ensemble to {target} candidates "
                    f"(round {round_index}, {requested} new)"
                )
            records, _ = sample_xtb_dimer_environments(
                smiles, run_dir, fidelity="accurate", monomer_xyz=monomer_xyz,
                charge=charge, multiplicity=multiplicity, ncores=ncores,
                seed=42 + 1009 * round_index, candidate_count=requested,
                include_trimers=include_trimers,
                candidate_id_prefix=f"extension-r{round_index:02d}",
                acquisition_round=round_index, merge_existing=True,
                progress=progress,
            )
        sample_xtb_snapshot_frequencies(
            records, run_dir, charge=charge, multiplicity=multiplicity,
            ncores=ncores, progress=progress,
        )
        transfer = build_frequency_transfer_artifacts(run_dir)
        convergence_path = run_dir / "frequency-transfer-convergence.json"
        convergence = json.loads(convergence_path.read_text())
        acquisition = convergence.get("acquisition_round_diagnostics") or {}
        rounds.append({
            "round": round_index,
            "cumulative_candidates": len(records),
            "transfer_status": transfer.get("status"),
            "acquisition_round_converged": acquisition.get("converged", False),
            "ensemble_converged": convergence.get("converged", False),
            "bootstrap_precision_pass": convergence.get(
                "bootstrap_precision_pass", False
            ),
            "consecutive_passing_comparisons": acquisition.get(
                "consecutive_passing_comparisons", 0
            ),
            "out_of_domain_candidates": convergence.get("out_of_domain_candidates"),
            "coverage": {
                "by_topology": {
                    topology: sum(
                        item.get("topology", "dimer") == topology
                        and int(item.get("acquisition_round", 0)) == round_index
                        for item in records
                    )
                    for topology in sorted({
                        item.get("topology", "dimer") for item in records
                        if int(item.get("acquisition_round", 0)) == round_index
                    })
                },
                "association_classes": {
                    label: sum(
                        int(item.get("acquisition_round", 0)) == round_index
                        and (
                            (float((item.get("environment_features") or {}).get("h_bond_distance_angstrom", 9.0)) < 1.85 and label == "strong")
                            or (1.85 <= float((item.get("environment_features") or {}).get("h_bond_distance_angstrom", 9.0)) <= 2.20 and label == "intermediate")
                            or (float((item.get("environment_features") or {}).get("h_bond_distance_angstrom", 9.0)) > 2.20 and label == "weak")
                        )
                        for item in records
                    )
                    for label in ("strong", "intermediate", "weak")
                },
            },
        })
        if convergence.get("converged") is True:
            break
    final_convergence_path = run_dir / "frequency-transfer-convergence.json"
    final_convergence = (
        json.loads(final_convergence_path.read_text())
        if final_convergence_path.is_file() else {}
    )
    final_converged = final_convergence.get("converged") is True
    numerical_stability = (
        (final_convergence.get("acquisition_round_diagnostics") or {}).get("converged")
        is True
    )
    report = {
        **plan,
        "status": (
            "converged" if final_converged
            else "maximum_candidates_reached_statistically_imprecise"
            if numerical_stability
            else "maximum_candidates_reached_unconverged"
        ),
        "numerically_stable_across_acquisition_rounds": numerical_stability,
        "bootstrap_precision_pass": final_convergence.get(
            "bootstrap_precision_pass", False
        ),
        "rounds": rounds,
        "final_candidates": len(_retained_sampling_records(run_dir)),
        "plan_artifact": str(plan_path),
    }
    report_path = retained_report_path
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(report_path)
    return report
