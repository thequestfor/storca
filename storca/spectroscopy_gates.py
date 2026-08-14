"""Promotion controller for direct, adaptive, and bulk-liquid spectroscopy."""

from __future__ import annotations

import json
from pathlib import Path


def evaluate_spectroscopy_gate_sequence(run_dir: Path) -> dict:
    """Read retained artifacts and authorize at most the next expensive stage."""
    run_dir = Path(run_dir)
    direct_path = run_dir / "direct-local-dft-ensemble.json"
    adaptive_path = run_dir / "local-mode-acquisition.json"
    bulk_path = run_dir / "bulk-embedding-acceptance.json"
    direct = json.loads(direct_path.read_text()) if direct_path.is_file() else None
    adaptive = json.loads(adaptive_path.read_text()) if adaptive_path.is_file() else None
    bulk = json.loads(bulk_path.read_text()) if bulk_path.is_file() else None

    if direct is None:
        state, authorized, reason = "direct_local_dft_required", "direct_local_dft", "no_direct_local_dft_artifact"
    elif direct.get("status") != "accepted":
        # A failed accuracy gate justifies a small bulk prototype, while a
        # missing validation/coverage gate authorizes no extra expense.
        gates = direct.get("acceptance_gates") or {}
        accuracy_evaluated = gates.get("nist_oh_center_improves", {}).get("status") != "not_evaluated_missing_reference"
        coverage_passed = gates.get("minimum_three_valid_environments_per_coordination_class", {}).get("passed") is True
        modes_passed = gates.get("all_expected_local_modes_step_validated", {}).get("passed") is True
        if accuracy_evaluated and coverage_passed and modes_passed:
            state, authorized, reason = "direct_accuracy_not_improved", "bulk_prototype", "direct_cluster_model_did_not_improve_oh_accuracy"
        else:
            state, authorized, reason = "direct_gate_failed_closed", None, "direct_validation_or_coverage_incomplete"
    elif adaptive is None:
        state, authorized, reason = "adaptive_expansion_justified", "adaptive_local_mode_acquisition", "direct_oh_accuracy_improved"
    else:
        convergence = adaptive.get("convergence") or {}
        summary = convergence.get("summary") or {}
        if summary.get("converged") is True:
            state, authorized, reason = "cluster_model_retained", None, "direct_distribution_converged"
        elif summary.get("promotion_to_bulk_embedding") is True:
            state, authorized, reason = "bulk_embedding_justified", "bulk_embedding", "cluster_distribution_failed_to_converge_at_hard_budget"
        else:
            state, authorized, reason = "adaptive_batch_required", "adaptive_local_mode_acquisition", "adaptive_hard_budget_not_reached"
    if bulk is not None:
        if bulk.get("status") == "accepted":
            state, authorized, reason = "bulk_embedding_accepted", None, "independent_seed_and_heldout_gates_passed"
        else:
            state, authorized, reason = "bulk_embedding_failed_closed", None, "bulk_validation_or_heldout_gain_incomplete"
    report = {
        "schema_version": 1, "kind": "gated_spectroscopy_implementation_sequence",
        "state": state, "authorized_next_expensive_stage": authorized,
        "reason": reason,
        "artifacts": {"direct": str(direct_path), "adaptive": str(adaptive_path), "bulk": str(bulk_path)},
    }
    path = run_dir / "spectroscopy-gate-sequence.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(path)
    return report
