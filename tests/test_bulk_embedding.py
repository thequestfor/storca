import numpy as np
from pathlib import Path

from storca.bulk_embedding import (
    BulkEmbeddingConfig, build_periodic_methanol_box,
    evaluate_bulk_embedding_acceptance, extract_central_solvation_shell,
    methanol_box_length_angstrom, validate_nvt_trajectory,
    write_electrostatic_embedding_input,
)


def test_periodic_box_matches_density_and_is_reproducible(tmp_path):
    config = BulkEmbeddingConfig(molecule_count=8)
    first = build_periodic_methanol_box(tmp_path / "one", config=config, seed=7)
    second = build_periodic_methanol_box(tmp_path / "two", config=config, seed=7)
    assert first["box_vectors_angstrom"] == second["box_vectors_angstrom"]
    assert methanol_box_length_angstrom(8, config.density_g_cm3) == first["box_vectors_angstrom"][0][0]


def test_shell_and_embedding_input_are_explicit(tmp_path):
    symbols = ["C", "O", "H", "H", "H", "H"] * 2
    coordinates = np.zeros((12, 3))
    coordinates[2, 0] = 1.0
    coordinates[6:, 0] = 2.5
    coordinates[8, 0] = 3.5
    shell = extract_central_solvation_shell(symbols, coordinates, central_molecule=0, box_length_angstrom=20)
    assert shell["qm_molecules"] == [0, 1]
    manifest = write_electrostatic_embedding_input(tmp_path, symbols, coordinates, shell)
    assert manifest["ncores"] == 8
    assert "pointcharges" in (tmp_path / "embedded-gradient.inp").read_text().lower()


def test_radial_qm_shell_includes_non_hbonded_oxygen_neighbor():
    symbols = ["C", "O", "H", "H", "H", "H"] * 2
    coordinates = np.zeros((12, 3))
    coordinates[2] = [1.0, 0.0, 0.0]
    coordinates[6:] += [0.0, 3.2, 0.0]
    coordinates[8] += [1.0, 0.0, 0.0]
    shell = extract_central_solvation_shell(
        symbols, coordinates, central_molecule=0, box_length_angstrom=20,
        config=BulkEmbeddingConfig(qm_shell_policy="radial_first_shell"),
    )
    assert shell["hydrogen_bonds"] == []
    assert shell["qm_molecules"] == [0, 1]


def test_embedding_rebuilds_molecule_split_across_periodic_boundary(tmp_path):
    symbols = ["C", "O", "H", "H", "H", "H"] * 2
    coordinates = np.zeros((12, 3))
    coordinates[1] = [1.0, 1.0, 1.0]
    coordinates[2] = [2.0, 1.0, 1.0]
    coordinates[7] = [9.8, 1.0, 1.0]
    coordinates[8] = [0.8, 1.0, 1.0]
    shell = extract_central_solvation_shell(
        symbols, coordinates, central_molecule=0, box_length_angstrom=10,
        config=BulkEmbeddingConfig(qm_shell_policy="radial_first_shell"),
    )
    manifest = write_electrostatic_embedding_input(tmp_path, symbols, coordinates, shell)
    rows = Path(manifest["qm_geometry"]).read_text().splitlines()[2:]
    points = np.asarray([[float(value) for value in row.split()[1:]] for row in rows])
    assert np.isclose(np.linalg.norm(points[8] - points[7]), 1.0)


def test_trajectory_and_two_seed_acceptance():
    records = [{
        "temperature_K": 298.15 + (-1) ** i, "density_g_cm3": 0.7866,
        "potential_energy": -100 + (-1) ** i, "h_bonds_per_molecule": 1.8,
    } for i in range(20)]
    validation = validate_nvt_trajectory(records)
    assert validation["status"] == "validated"
    seeds = [{
        "seed": seed, "trajectory_validation": validation,
        "finite_difference_validation_passed": True, "held_out_gain_reproduced": True,
        "maximum_block_center_change_cm-1": 2, "maximum_block_width_change_cm-1": 5,
    } for seed in (1, 2)]
    assert evaluate_bulk_embedding_acceptance(seeds)["status"] == "accepted"
