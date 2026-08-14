import json
from pathlib import Path
import tempfile

import numpy as np

from storca.cluster_ir import generate_stratified_dimer_poses
from storca.xtb_sampling import run_xtb_candidate
from storca.environment_selection import select_xtb_environment_representatives
from storca.xtb_trajectory import (
    TRAJECTORY_POPULATION_MODEL,
    decorrelated_frame_indices,
    read_xyz_trajectory,
    run_restrained_xtb_trajectory,
    statistical_inefficiency,
    xtb_trajectory_defaults,
)
from storca.xtb_frequencies import run_xtb_snapshot_frequency


def _fake_xtb(path: Path) -> Path:
    script = path / "xtb"
    script.write_text("""#!/usr/bin/env python3
import json, pathlib, shutil, sys
if '--version' in sys.argv:
    print('xtb version 6.7.1')
    raise SystemExit(0)
source = pathlib.Path(sys.argv[1])
if '--opt' in sys.argv:
    shutil.copyfile(source, 'xtbopt.xyz')
    pathlib.Path('xtbout.json').write_text(json.dumps({'total energy': -16.5}))
else:
    text = source.read_text()
    pathlib.Path('xtb.trj').write_text(text * 8)
    pathlib.Path('xtbmdok').touch()
print('normal termination of xtb')
""")
    script.chmod(0o755)
    return script


def test_trajectory_defaults_are_298_k_and_occupancy_weighted():
    settings = xtb_trajectory_defaults("auto")
    assert settings["temperature_K"] == 298.15
    assert settings["trajectory_seed_count"] == 4
    assert settings["time_step_fs"] == 0.5
    assert settings["hydrogen_mass_amu"] == 4
    assert settings["population_model"] == TRAJECTORY_POPULATION_MODEL


def test_decorrelation_detects_correlated_series_and_discards_equilibration():
    rng = np.random.default_rng(5)
    values = np.cumsum(rng.normal(size=200))
    assert statistical_inefficiency(values) > 1.0
    frames = [
        {"coordinates": np.asarray([[0.0, 0.0, 0.0], [0.0, 0.0, value]])}
        for value in values
    ]
    indices, report = decorrelated_frame_indices(
        frames, equilibration_fraction=0.25, maximum_snapshots=10,
    )
    assert indices[0] >= 50
    assert len(indices) <= 10
    assert report["decorrelation_stride_frames"] >= 2


def test_restrained_trajectory_is_resumable_and_extracts_snapshots():
    pose = generate_stratified_dimer_poses("CO", candidate_count=1)[0]
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        executable = _fake_xtb(root)
        seed = run_xtb_candidate(
            pose, root / "seed", executable=str(executable), executable_version="6.7.1",
        )
        settings = {**xtb_trajectory_defaults("auto"), "production_time_ps": 0.08}
        snapshots, report = run_restrained_xtb_trajectory(
            seed, root / "trajectory", settings=settings, executable=str(executable),
            executable_version="6.7.1", maximum_snapshots=3,
        )
        reused, reused_report = run_restrained_xtb_trajectory(
            seed, root / "trajectory", settings=settings, executable=str(executable),
            executable_version="6.7.1", maximum_snapshots=3,
        )
        parsed = read_xyz_trajectory(root / "trajectory" / "xtb.trj")
        contract = json.loads((root / "trajectory" / "calculation-contract.json").read_text())
    assert report["trajectory_status"] == "completed"
    assert len(parsed) == 8
    assert 1 <= len(snapshots) <= 3
    assert reused == snapshots
    assert reused_report == report
    assert contract["settings"]["temperature_K"] == 298.15
    assert all(item["population_model"] == TRAJECTORY_POPULATION_MODEL for item in snapshots)


def test_representative_clusters_preserve_trajectory_occupancy_mass():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        paths = []
        for index in range(3):
            path = root / f"candidate-{index}.xyz"
            path.write_text(f"2\nframe\nO 0 0 0\nH 0 0 {1.0 + 0.01 * index}\n")
            paths.append(path)
        records = []
        for index, (distance, weight) in enumerate(((1.60, 0.1), (1.65, 0.2), (2.50, 0.7))):
            records.append({
                "candidate_id": f"candidate-{index}",
                "sampling_status": "retained",
                "optimized_xyz": str(paths[index]),
                "cluster_size": 2,
                "topology": "dimer",
                "site_identity": {"site": 1},
                "relative_xtb_energy_kcal_mol": 0.0,
                "population_weight": weight,
                "population_model": TRAJECTORY_POPULATION_MODEL,
                "population_warning": "conditional occupancy",
                "environment_features": {
                    "h_bond_distance_angstrom": distance,
                    "donor_h_acceptor_angle_degrees": 170.0,
                    "donor_acceptor_distance_angstrom": distance + 1.0,
                    "intermolecular_orientation_degrees": 0.0,
                    "heavy_atom_rmsd_angstrom": 0.0,
                    "estimated_local_frequency_cm-1": None,
                },
            })
        _, manifest = select_xtb_environment_representatives(
            records, root, representative_count=2,
        )
        entries = json.loads(manifest.read_text())["conformers"]
    assert sorted(round(item["population_weight"], 8) for item in entries) == [0.3, 0.7]
    assert all(item["population_model"] == TRAJECTORY_POPULATION_MODEL for item in entries)


def test_cached_frequency_record_refreshes_population_provenance():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        executable = _fake_xtb(root)
        source = root / "snapshot" / "optimized.xyz"
        source.parent.mkdir()
        source.write_text("2\nsnapshot\nO 0 0 0\nH 0 0 1\n")
        record = {
            "candidate_id": "snapshot-1",
            "sampling_status": "retained",
            "optimized_xyz": str(source),
            "cluster_size": 1,
            "topology": "monomer",
            "population_weight": 1.0,
            "population_model": TRAJECTORY_POPULATION_MODEL,
            "population_warning": "conditional occupancy",
            "population_weight_basis": "decorrelated frames",
        }
        frequency_dir = source.parent / "xtb-frequency"
        frequency_dir.mkdir()
        # Exercise the resume path without needing to synthesize an xTB g98 Hessian.
        from storca.xtb_frequencies import XTB_FREQUENCY_SCHEMA_VERSION, _sha256_file
        input_path = frequency_dir / "input.xyz"
        input_path.write_text(source.read_text())
        contract = {
            "schema_version": XTB_FREQUENCY_SCHEMA_VERSION,
            "candidate_id": "snapshot-1",
            "input_xyz_sha256": _sha256_file(input_path),
            "xtb_executable": str(executable.resolve()),
            "xtb_version": "6.7.1",
            "method": "GFN2-xTB",
            "calculation": "unrestrained_snapshot_hessian",
            "charge": 0,
            "unpaired_electrons": 0,
            "ncores": 1,
        }
        (frequency_dir / "calculation-contract.json").write_text(json.dumps(contract))
        (frequency_dir / "frequency-record.json").write_text(json.dumps({
            "candidate_id": "snapshot-1", "frequency_status": "completed",
        }))
        refreshed = run_xtb_snapshot_frequency(
            record, executable=str(executable), executable_version="6.7.1",
        )
    assert refreshed["population_weight"] == 1.0
    assert refreshed["population_model"] == TRAJECTORY_POPULATION_MODEL
    assert refreshed["population_warning"] == "conditional occupancy"
    assert refreshed["population_weight_basis"] == "decorrelated frames"
