import json
from pathlib import Path

import numpy as np

from storca.bulk_embedding import BulkEmbeddingConfig, build_periodic_methanol_box
from storca.gromacs_bulk import (
    GromacsMethanolConfig, METHANOL_ATOMS, prepare_bulk_seed,
    write_mdp_files, write_methanol_gro, write_oplsaa_methanol_topology,
)


def test_force_field_contract_is_neutral_and_uses_pinned_opls_types():
    assert np.isclose(sum(item[2] for item in METHANOL_ATOMS), 0.0)
    assert [item[1] for item in METHANOL_ATOMS] == [
        "opls_157", "opls_154", "opls_155", "opls_156", "opls_156", "opls_156",
    ]


def test_box_conversion_preserves_molecules_and_density(tmp_path):
    box = build_periodic_methanol_box(
        tmp_path, config=BulkEmbeddingConfig(molecule_count=8), seed=4,
    )
    gro = write_methanol_gro(Path(box["initial_geometry"]), tmp_path / "initial.gro")
    lines = gro.read_text().splitlines()
    assert int(lines[1]) == 48
    assert lines[2][5:10].strip() == "MEOH"
    assert np.isclose(float(lines[-1].split()[0]) * 10, box["box_vectors_angstrom"][0][0], atol=0.006)


def test_topology_declares_all_bonds_angles_and_torsions(tmp_path):
    result = write_oplsaa_methanol_topology(tmp_path, molecule_count=216)
    text = Path(result["molecule_include"]).read_text()
    assert text.count("opls_15") >= 6
    assert "[ bonds ]" in text and "[ angles ]" in text and "[ dihedrals ]" in text
    assert "MEOH 216" in Path(result["topology"]).read_text()


def test_mdp_schedule_uses_flexible_one_fs_contract(tmp_path):
    paths = write_mdp_files(tmp_path, seed=9, config=GromacsMethanolConfig(molecule_count=8))
    production = Path(paths["production"]).read_text()
    assert "dt = 0.001" in production
    assert "constraints = none" in production
    assert "coulombtype = PME" in production


def test_real_grompp_preflight_for_small_box(tmp_path):
    gmx = Path("/home/entity/storca/.miniforge/envs/storca_gromacs/bin/gmx")
    if not gmx.is_file():
        return
    build_periodic_methanol_box(
        tmp_path, config=BulkEmbeddingConfig(molecule_count=8), seed=3,
    )
    report = prepare_bulk_seed(
        tmp_path, gmx, config=GromacsMethanolConfig(molecule_count=8, cutoff_nm=0.35),
    )
    assert report["status"] == "prepared"
    assert (tmp_path / "minim.tpr").is_file()
    retained = json.loads((tmp_path / "gromacs-preflight.json").read_text())
    assert retained["force_field_contract"]["combination_rule"] == 3
