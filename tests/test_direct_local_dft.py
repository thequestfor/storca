import json
from pathlib import Path

import numpy as np

from storca.direct_local_dft import (
    DirectLocalDFTConfig, assemble_direct_local_dft,
    coordination_class_coverage, plan_six_representative_direct_local_dft,
    write_direct_local_dft_artifacts,
)


def _conformer(folder: Path, position: int, *, valid=True, cluster_size=2):
    job = folder / f"conf-{position}"
    job.mkdir()
    xyz = job / "input.xyz"
    xyz.write_text("2\ntest\nO 0 0 0\nH 0 0 1\n")
    output = job / "freq.out"
    output.write_text("retained")
    bond = {
        "bond_class": "O-H", "spectral_band_class": "hydrogen_bonded_oh",
        "coordination_class": "donor", "heavy_atom": 0, "hydrogen_atom": 1,
    }
    local = job / "local-mode-finite-differences" / "local-modes.json"
    local.parent.mkdir()
    local.write_text(json.dumps({"modes": [{
        "bond": bond, "status": "validated" if valid else "failed_validation",
        "frequency_cm-1": 3300.0 + position, "intensity_km_mol": 100.0,
        "frequency_step_disagreement_cm-1": 2.0,
    }]}))
    return {
        "optimized_xyz": str(xyz), "frequency_output": str(output),
        "local_mode_output": str(local), "local_stretch_bonds": [bond],
        "population_weight": 1.0 / 6, "cluster_size": cluster_size,
        "independent_environment_id": f"env-{position}",
        "ir_modes": [
            {"mode": 1, "freq": 1000.0, "intensity": 10.0},
            {"mode": 2, "freq": 3400.0, "intensity": 9999.0},
        ],
    }


def test_assembler_replaces_only_oh_and_never_uses_xtb(tmp_path):
    conformers = [_conformer(tmp_path, index) for index in range(1, 4)]
    grid, signal, modes, provenance = assemble_direct_local_dft(conformers)
    assert signal[np.argmin(abs(grid - 1000.0))] > 0
    assert provenance["xtb_frequency_correction_used_in_oh_region"] is False
    assert all("xtb" not in row["frequency_source"] for row in modes)
    assert sum(row.get("frequency_cm-1") == 3400.0 for row in modes) == 0


def test_coverage_counts_independent_environments_and_ess(tmp_path):
    conformers = [_conformer(tmp_path, index) for index in range(1, 4)]
    _, _, modes, _ = assemble_direct_local_dft(conformers)
    coverage = coordination_class_coverage(
        row for row in modes if row.get("validation_status")
    )
    detail = next(iter(coverage.values()))
    assert detail["valid_independent_environments"] == 3
    assert np.isclose(detail["effective_sample_size"], 3.0)


def test_direct_artifacts_fail_closed_without_nist(tmp_path):
    conformers = [_conformer(tmp_path, index) for index in range(1, 4)]
    report = write_direct_local_dft_artifacts(
        tmp_path, conformers,
        config=DirectLocalDFTConfig(expected_valid_local_modes=3),
    )
    assert report["status"] == "failed_closed"
    assert Path(report["artifacts"]["intrinsic_spectrum"]).is_file()
    assert Path(report["artifacts"]["transmission_spectrum"]).is_file()
    assert report["acceptance_gates"]["nist_oh_center_improves"]["passed"] is False


def test_24_call_plan_exposes_nine_mode_inconsistency(tmp_path):
    conformers = [
        _conformer(tmp_path, index, valid=index <= 3) for index in range(1, 7)
    ]
    for conformer in conformers[3:]:
        conformer["local_stretch_bonds"] *= 3
    plan = plan_six_representative_direct_local_dft(conformers)
    assert plan["required_additional_orca_invocations"] == 36
    assert plan["status"] == "failed_closed_budget_or_coverage"
