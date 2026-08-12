import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from storca.frequency_transfer import (apply_finite_difference_pair_overrides,
                                       evaluate_transferred_ensemble_convergence,
                                       transfer_snapshot_modes, validate_frequency_transfer,
                                       _write_or_activate_hybrid_spectrum)
from storca.xtb_frequencies import parse_xtb_g98, run_xtb_snapshot_frequency


def minimal_g98() -> str:
    return """ Entering Gaussian System
                         Standard orientation:
 --------------------------------------------------------------------
  Center     Atomic     Atomic              Coordinates (Angstroms)
  Number     Number      Type              X           Y           Z
 --------------------------------------------------------------------
    1          8             0        0.000000    0.000000    0.000000
    2          1             0        0.960000    0.000000    0.000000
    3          1             0       -0.240000    0.930000    0.000000
 --------------------------------------------------------------------
 Frequencies --  1600.0000              3500.0000              3650.0000
 Red. masses --     1.0000                 1.0000                 1.0000
 Frc consts  --     1.0000                 1.0000                 1.0000
 IR Inten    --    10.0000                20.0000                30.0000
 Atom AN      X      Y      Z        X      Y      Z        X      Y      Z
   1   8     0.00   0.10   0.00    -0.50   0.00   0.00     0.20  -0.40   0.00
   2   1     0.00  -0.10   0.00     1.00   0.00   0.00    -0.10   0.20   0.00
   3   1     0.00  -0.10   0.00     0.00   0.00   0.00    -0.10   0.20   0.00
"""


class XTBFrequencyTests(unittest.TestCase):
    def test_g98_parser_returns_geometry_modes_and_intensities(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "g98.out"
            path.write_text(minimal_g98())
            parsed = parse_xtb_g98(path)
        self.assertEqual(parsed.mode_set.atom_labels, ("O", "H", "H"))
        self.assertEqual(parsed.mode_set.vectors.shape, (3, 3, 3))
        self.assertEqual(parsed.mode_set.frequencies_cm_1.tolist(), [1600.0, 3500.0, 3650.0])
        self.assertEqual(parsed.intensities_km_mol.tolist(), [10.0, 20.0, 30.0])

    def test_snapshot_hessian_is_unrestrained_and_resumable(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            candidate = folder / "candidate-0001"
            candidate.mkdir()
            xyz = candidate / "optimized.xyz"
            xyz.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
            executable = folder / "fake-xtb"
            executable.write_text(
                "#!/bin/sh\n"
                "if [ \"$1\" = \"--version\" ]; then echo 'xtb version 6.7.1'; exit 0; fi\n"
                "cp ../fixture.g98 g98.out\n"
                "echo 'normal termination of xtb'\n"
            )
            executable.chmod(0o755)
            (candidate / "fixture.g98").write_text(minimal_g98())
            sampling = {
                "candidate_id": "candidate-0001", "sampling_status": "retained",
                "optimized_xyz": str(xyz), "environment_features": {},
                "local_stretch_bonds": [{
                    "bond_class": "O-H:0-1", "heavy_atom": 0, "hydrogen_atom": 1,
                    "molecule_index": 0, "spectral_band_class": "non_donating_oh",
                    "coordination_class": "free",
                }],
            }
            first = run_xtb_snapshot_frequency(sampling, executable=str(executable))
            second = run_xtb_snapshot_frequency(sampling, executable=str(executable))
            contract = json.loads((candidate / "xtb-frequency" / "calculation-contract.json").read_text())
        self.assertEqual(first, second)
        self.assertEqual(first["frequency_status"], "completed")
        self.assertEqual(contract["calculation"], "unrestrained_snapshot_hessian")
        self.assertNotIn("xcontrol", json.dumps(contract).lower())


class FrequencyTransferTests(unittest.TestCase):
    @staticmethod
    def pairs():
        output = []
        for index, correction in enumerate((98.0, 100.0, 102.0), start=1):
            output.append({
                "representative_id": f"rep-{index}", "mode_class": "O-H:hb:donor",
                "xtb_frequency_cm-1": 3400.0 + 10.0 * index,
                "dft_frequency_cm-1": 3400.0 + 10.0 * index + correction,
                "additive_correction_cm-1": correction,
                "mode_character_similarity": 0.9,
                "environment_features": {"h_bond_distance_angstrom": 1.7 + 0.1 * index},
            })
        return output

    def test_leave_one_out_validation_requires_improvement(self):
        validation = validate_frequency_transfer(self.pairs())
        self.assertEqual(validation["status"], "validated")
        result = validation["classes"]["O-H:hb:donor"]
        self.assertLess(result["corrected_loo_mae_cm-1"], result["raw_xtb_loo_mae_cm-1"])

    def test_transfer_is_fail_closed_with_too_few_representatives(self):
        pairs = self.pairs()[:2]
        validation = validate_frequency_transfer(pairs)
        self.assertEqual(validation["status"], "insufficient_dft_transfer_validation")
        records = [{
            "candidate_id": "snapshot", "frequency_status": "completed",
            "environment_features": {"h_bond_distance_angstrom": 1.8},
            "modes": [{
                "mode": 1, "freq": 3400.0, "intensity": 10.0,
                "bond_class": "O-H", "spectral_band_class": "hb",
                "coordination_class": "donor",
            }],
        }]
        transferred = transfer_snapshot_modes(records, pairs, validation)
        self.assertEqual(transferred[0]["modes"], [])

    def test_validated_local_finite_difference_overrides_snapshot_hessian_pair(self):
        pair = self.pairs()[0]
        pair.update({
            "bond": {"molecule_index": 0, "heavy_atom": 1, "hydrogen_atom": 5},
            "dft_intensity_km_mol": 12.0,
        })
        payload = {"modes": [{
            "status": "validated", "frequency_cm-1": 3333.0,
            "intensity_km_mol": 45.0, "heavy_atom": 1, "hydrogen_atom": 5,
            "displacement_stability_status": "passed",
            "frequency_step_disagreement_cm-1": 1.2,
            "bond": {"molecule_index": 0, "heavy_atom": 1, "hydrogen_atom": 5},
        }]}
        updated = apply_finite_difference_pair_overrides([pair], payload)[0]
        self.assertEqual(updated["dft_frequency_cm-1"], 3333.0)
        self.assertEqual(updated["dft_intensity_km_mol"], 45.0)
        self.assertEqual(
            updated["dft_frequency_source"],
            "orca_projected_local_mode_finite_difference",
        )

    def test_poor_snapshot_pairs_without_finite_difference_are_excluded(self):
        selected = self.pairs()[0]
        selected["bond"] = {"molecule_index": 0, "heavy_atom": 1, "hydrogen_atom": 5}
        unselected = {**self.pairs()[1], "bond": {
            "molecule_index": 1, "heavy_atom": 7, "hydrogen_atom": 11,
        }}
        payload = {"modes": [{
            "status": "validated", "frequency_cm-1": 3333.0,
            "intensity_km_mol": 45.0, "heavy_atom": 1, "hydrogen_atom": 5,
            "displacement_stability_status": "passed",
            "frequency_step_disagreement_cm-1": 1.0,
            "bond": {"molecule_index": 0, "heavy_atom": 1, "hydrogen_atom": 5},
        }]}
        updated = apply_finite_difference_pair_overrides(
            [selected, unselected], payload,
            exclude_unvalidated_snapshot_pairs=True,
        )
        self.assertEqual(len(updated), 1)
        self.assertEqual(updated[0]["bond"]["molecule_index"], 0)

    def test_transferred_ensemble_requires_two_stable_batch_comparisons(self):
        records = [{
            "candidate_id": f"candidate-{index}", "cluster_size": 2,
            "acquisition_round": index // 10,
            "modes": [{
                "mode_class": "O-H:hb:donor", "freq": 3400.0,
                "intensity": 10.0,
            }],
        } for index in range(30)]
        report = evaluate_transferred_ensemble_convergence(records, batch_size=10)
        self.assertTrue(report["converged"])
        self.assertEqual(report["batches"][-1]["comparison"]["consecutive_passes"], 2)
        final_band = report["batches"][-1]["bands"]["O-H:hb:donor"]
        self.assertEqual(final_band["independent_environments"], 30)
        self.assertEqual(final_band["samples"], 300)
        self.assertIsNotNone(final_band["center_95_ci_cm-1"])
        self.assertTrue(report["acquisition_round_diagnostics"]["converged"])
        self.assertTrue(report["bootstrap_precision_pass"])

    def test_multiple_local_modes_do_not_inflate_independent_environment_count(self):
        records = [{
            "candidate_id": f"candidate-{index}", "cluster_size": 3,
            "modes": [
                {"mode_class": "O-H:hb:donor", "freq": 3390.0, "intensity": 10.0},
                {"mode_class": "O-H:hb:donor", "freq": 3410.0, "intensity": 12.0},
            ],
        } for index in range(30)]
        report = evaluate_transferred_ensemble_convergence(records, batch_size=10)
        band = report["batches"][-1]["bands"]["O-H:hb:donor"]
        self.assertEqual(band["independent_environments"], 30)
        self.assertEqual(band["center_cm-1"], 3400.0)

    def test_unconverged_hybrid_is_diagnostic_and_does_not_replace_display(self):
        with tempfile.TemporaryDirectory() as temp:
            run_dir = Path(temp)
            displayed = run_dir / "spectrum.csv"
            displayed.write_text("retained display\n")
            grid = np.asarray([1000.0, 1001.0])
            intrinsic = np.asarray([0.2, 0.3])
            result = _write_or_activate_hybrid_spectrum(
                run_dir, grid, intrinsic, activate=False,
            )
            self.assertEqual(displayed.read_text(), "retained display\n")
            self.assertEqual(result["hybrid_spectrum_role"], "diagnostic_not_displayed")
            self.assertTrue((run_dir / "spectrum_hybrid_multifidelity_intrinsic.csv").is_file())


if __name__ == "__main__":
    unittest.main()
