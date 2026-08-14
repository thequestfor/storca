import json

from storca.spectroscopy_gates import evaluate_spectroscopy_gate_sequence


def test_sequence_starts_at_direct_stage(tmp_path):
    report = evaluate_spectroscopy_gate_sequence(tmp_path)
    assert report["authorized_next_expensive_stage"] == "direct_local_dft"


def test_accepted_direct_authorizes_adaptive_only(tmp_path):
    (tmp_path / "direct-local-dft-ensemble.json").write_text(json.dumps({
        "status": "accepted", "acceptance_gates": {},
    }))
    report = evaluate_spectroscopy_gate_sequence(tmp_path)
    assert report["authorized_next_expensive_stage"] == "adaptive_local_mode_acquisition"


def test_unconverged_hard_budget_authorizes_bulk(tmp_path):
    (tmp_path / "direct-local-dft-ensemble.json").write_text(json.dumps({"status": "accepted"}))
    (tmp_path / "local-mode-acquisition.json").write_text(json.dumps({
        "convergence": {"summary": {"converged": False, "promotion_to_bulk_embedding": True}},
    }))
    report = evaluate_spectroscopy_gate_sequence(tmp_path)
    assert report["authorized_next_expensive_stage"] == "bulk_embedding"


def test_failed_bulk_attempt_does_not_reauthorize_itself(tmp_path):
    (tmp_path / "direct-local-dft-ensemble.json").write_text(json.dumps({"status": "accepted"}))
    (tmp_path / "bulk-embedding-acceptance.json").write_text(json.dumps({"status": "failed_closed"}))
    report = evaluate_spectroscopy_gate_sequence(tmp_path)
    assert report["state"] == "bulk_embedding_failed_closed"
    assert report["authorized_next_expensive_stage"] is None
