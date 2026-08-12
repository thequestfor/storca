"""Budget-aware local-mode fallback orchestration for environment representatives."""

from __future__ import annotations

import json
from pathlib import Path

from .local_mode_fd import calculate_orca_local_modes


ENVIRONMENT_LOCAL_MODE_SCHEMA_VERSION = 1
MINIMUM_REPRESENTATIVES_PER_CLASS = 3
INVOCATIONS_PER_BOND = 4


def _mode_class(bond: dict) -> str:
    return ":".join((
        str(bond.get("bond_class", "X-H")),
        str(bond.get("spectral_band_class", "xh_stretch")),
        str(bond.get("coordination_class", "unclassified")),
    ))


def plan_environment_local_mode_fallbacks(
    conformers: list[dict], *, maximum_orca_invocations: int,
) -> dict:
    """Choose three distinct poor-stationarity representatives per mode class."""
    if maximum_orca_invocations < 0:
        raise ValueError("Local-mode ORCA invocation allowance cannot be negative")
    available: dict[str, list[tuple[int, dict]]] = {}
    for position, conformer in enumerate(conformers):
        reliability = conformer.get("snapshot_hessian_reliability") or {}
        if reliability.get("stationarity_status") != "poor":
            continue
        seen_classes = set()
        for bond in conformer.get("local_stretch_bonds") or []:
            mode_class = _mode_class(bond)
            # One bond per class and representative is enough for independent
            # transfer coverage; cyclic trimers must not count three times.
            if mode_class in seen_classes:
                continue
            seen_classes.add(mode_class)
            available.setdefault(mode_class, []).append((position, dict(bond)))
    selected: dict[int, list[dict]] = {}
    class_reports = {}
    remaining = maximum_orca_invocations
    for mode_class in sorted(available):
        candidates = available[mode_class]
        required = MINIMUM_REPRESENTATIVES_PER_CLASS
        cost = required * INVOCATIONS_PER_BOND
        if len(candidates) < required:
            status = "insufficient_distinct_representatives"
            chosen = []
        elif remaining < cost:
            status = "insufficient_orca_invocation_allowance"
            chosen = []
        else:
            chosen = candidates[:required]
            remaining -= cost
            status = "selected"
            for position, bond in chosen:
                selected.setdefault(position, []).append(bond)
        class_reports[mode_class] = {
            "status": status,
            "available_distinct_representatives": len(candidates),
            "required_distinct_representatives": required,
            "selected_representative_positions": [position + 1 for position, _ in chosen],
            "orca_invocations": len(chosen) * INVOCATIONS_PER_BOND,
        }
    jobs = []
    for position, bonds in sorted(selected.items()):
        conformer = conformers[position]
        jobs.append({
            "conformer_position": position + 1,
            "independent_environment_id": conformer.get("independent_environment_id"),
            "frequency_output": conformer.get("frequency_output"),
            "xyz": str(conformer.get("optimized_xyz") or ""),
            "bonds": bonds,
            "orca_invocations": len(bonds) * INVOCATIONS_PER_BOND,
        })
    used = maximum_orca_invocations - remaining
    return {
        "schema_version": ENVIRONMENT_LOCAL_MODE_SCHEMA_VERSION,
        "kind": "environment_local_mode_fallback_plan",
        "maximum_orca_invocations": maximum_orca_invocations,
        "planned_orca_invocations": used,
        "unallocated_orca_invocations": remaining,
        "invocations_per_bond": INVOCATIONS_PER_BOND,
        "minimum_representatives_per_class": MINIMUM_REPRESENTATIVES_PER_CLASS,
        "classes": class_reports,
        "jobs": jobs,
        "status": (
            "planned" if jobs and all(item["status"] == "selected" for item in class_reports.values())
            else "partial_plan" if jobs else "insufficient_budget_or_coverage"
        ),
    }


def run_environment_local_mode_fallbacks(
    run_dir: Path, *, maximum_orca_invocations: int, ncores: int = 1,
    method_keywords: list[str] | None = None,
    progress=None,
) -> dict:
    """Execute selected local modes and retain a run-level process ledger."""
    run_dir = Path(run_dir)
    conformers_path = run_dir / "clusters" / "conformers.json"
    if not conformers_path.is_file():
        return {"status": "missing_environment_conformers"}
    conformers = json.loads(conformers_path.read_text())
    plan = plan_environment_local_mode_fallbacks(
        conformers, maximum_orca_invocations=maximum_orca_invocations,
    )
    plan_path = run_dir / "environment-local-mode-plan.json"
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    ledger_path = run_dir / "environment-local-mode-orca-invocations.json"
    results = []
    for job in plan["jobs"]:
        position = int(job["conformer_position"])
        xyz = Path(job["xyz"])
        frequency_output = Path(job["frequency_output"])
        if not xyz.is_file() or not frequency_output.is_file():
            results.append({**job, "status": "missing_representative_artifact"})
            continue
        job_dir = frequency_output.parent / "local-mode-finite-differences"
        if progress:
            progress(
                f"Calculating {len(job['bonds'])} local modes for environment "
                f"{position} ({job['orca_invocations']} ORCA gradient invocations)"
            )
        try:
            result = calculate_orca_local_modes(
                xyz, job["bonds"], job_dir,
                hard_orca_invocation_cap=maximum_orca_invocations,
                ncores=ncores, method_keywords=method_keywords,
                orca_ledger_path=ledger_path,
                ledger_namespace=f"environment-{position:03d}",
            )
            results.append({
                **job, "status": result["status"],
                "local_mode_output": str(job_dir / "local-modes.json"),
            })
        except Exception as error:
            results.append({**job, "status": "failed", "error": str(error)})
            break
    completed_by_class: dict[str, set[int]] = {}
    for result in results:
        if result.get("status") != "completed":
            continue
        for bond in result["bonds"]:
            completed_by_class.setdefault(_mode_class(bond), set()).add(
                int(result["conformer_position"])
            )
    completed_classes = sorted(
        mode_class for mode_class, positions in completed_by_class.items()
        if len(positions) >= MINIMUM_REPRESENTATIVES_PER_CLASS
    )
    report = {
        "schema_version": ENVIRONMENT_LOCAL_MODE_SCHEMA_VERSION,
        "kind": "environment_local_mode_fallback_execution",
        "status": (
            "completed" if results and all(item.get("status") == "completed" for item in results)
            else "partial_or_failed"
        ),
        "plan_artifact": str(plan_path),
        "orca_invocation_ledger": str(ledger_path),
        "completed_mode_classes": completed_classes,
        "results": results,
    }
    report_path = run_dir / "environment-local-modes.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(report_path)
    return report
