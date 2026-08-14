import numpy as np

from storca.adaptive_local_modes import (
    AdaptiveLocalModeConfig, cluster_bootstrap_center_ci,
    evaluate_adaptive_batches, run_adaptive_local_mode_acquisition,
    select_adaptive_local_mode_representatives,
)


def _record(index):
    return {
        "candidate_id": f"c{index}", "population_weight": 0.5 if index == 2 else 0.1,
        "topology": "trimer" if index >= 2 else "dimer",
        "environment_features": {
            "h_bond_distance_angstrom": 1.7 + index * 0.1,
            "donor_h_acceptor_angle_degrees": 140 + index * 5,
        },
        "local_stretch_bonds": [{
            "bond_class": "O-H", "spectral_band_class": "bonded",
            "coordination_class": "network" if index >= 2 else "donor",
        }],
    }


def test_selector_uses_geometry_and_occupancy_not_frequencies():
    records = [_record(index) for index in range(4)]
    matrix = np.abs(np.arange(4)[:, None] - np.arange(4)[None, :])
    selected, report = select_adaptive_local_mode_representatives(records, matrix, 2)
    assert 2 in selected
    assert report["xtb_frequency_used"] is False
    assert report["prohibited_selection_input"] == "xtb_frequency"


def test_bootstrap_resamples_trajectory_clusters():
    observations = [
        {"trajectory_cluster_id": f"cluster-{index}", "frequency_cm-1": 3300 + index, "weight": 1}
        for index in range(6)
    ]
    report = cluster_bootstrap_center_ci(observations, samples=20)
    assert report["resampling_unit"] == "trajectory_cluster"
    assert report["clusters"] == 6


def test_convergence_requires_two_consecutive_batches():
    grid = np.linspace(2800, 3900, 1101)
    spectrum = np.exp(-0.5 * ((grid - 3300) / 30) ** 2)
    classes = {"network": {"occupancy_fraction": 1, "effective_sample_size": 6, "occupancy_coverage_fraction": 1}}
    observations = [{"trajectory_cluster_id": f"c{i}", "frequency_cm-1": 3300, "weight": 1} for i in range(6)]
    batches = [{"representatives": endpoint, "grid": grid, "spectrum": spectrum, "coordination_classes": classes, "observations": observations} for endpoint in (6, 9, 12)]
    report = evaluate_adaptive_batches(batches, config=AdaptiveLocalModeConfig(bootstrap_samples=20))
    assert report["summary"]["converged"] is True
    assert len(report["batches"]) == 3


def test_runner_writes_converged_manifest_uncertainty_and_spectrum(tmp_path):
    selected = [_record(index) for index in range(12)]
    grid = np.linspace(2800, 3900, 101)
    spectrum = np.exp(-0.5 * ((grid - 3300) / 30) ** 2)

    def acquire(_, endpoint):
        return {
            "representatives": endpoint, "grid": grid.tolist(), "spectrum": spectrum.tolist(),
            "coordination_classes": {"network": {"occupancy_fraction": 1, "effective_sample_size": 6, "occupancy_coverage_fraction": 1}},
            "observations": [{"trajectory_cluster_id": f"c{i}", "frequency_cm-1": 3300, "weight": 1} for i in range(6)],
        }

    # Seed the retained initial-six result; acquisition itself adds batches of three.
    (tmp_path / "local-mode-acquisition.json").write_text(__import__("json").dumps({
        "completed_batches": [acquire(selected[:6], 6)],
    }))
    report = run_adaptive_local_mode_acquisition(
        tmp_path, selected, acquire_batch=acquire,
        config=AdaptiveLocalModeConfig(maximum_representatives=12, bootstrap_samples=20),
    )
    assert report["convergence"]["summary"]["converged"] is True
    assert (tmp_path / "local-mode-uncertainty-bands.json").is_file()
    assert (tmp_path / "final-local-mode-representative-manifest.json").is_file()
    assert (tmp_path / "spectrum_direct_local_dft_converged.csv").is_file()
