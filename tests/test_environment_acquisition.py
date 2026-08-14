import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from storca.environment_acquisition import (
    EnvironmentAcquisitionConfig, build_candidate_mode_profiles,
    record_acquisition_batch, reprioritize_pending_representatives,
    select_mode_coverage_representatives,
    significant_mode_classes, write_environment_acquisition_report,
)


def frequency_record(candidate_id, classes):
    return {
        "candidate_id": candidate_id, "frequency_status": "completed",
        "modes": [
            {
                "mode": index, "freq": frequency, "intensity": 10.0,
                "bond_class": "O-H", "spectral_band_class": spectral,
                "coordination_class": coordination,
            }
            for index, (spectral, coordination, frequency) in enumerate(classes)
        ],
    }


class EnvironmentAcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.records = [
            {
                "candidate_id": f"candidate-{index}",
                "topology": "linear_trimer" if index >= 3 else "dimer",
                "environment_features": {"h_bond_distance_angstrom": 1.6 + 0.1 * index},
            }
            for index in range(6)
        ]
        common = ("non_donating_oh", "acceptor_only", 3550.0)
        network = ("hydrogen_bonded_oh", "donor_acceptor", 3300.0)
        self.frequencies = [
            frequency_record(f"candidate-{index}", [common] + ([network] if index >= 3 else []))
            for index in range(6)
        ]
        self.matrix = np.abs(np.arange(6)[:, None] - np.arange(6)[None, :]).astype(float)

    def test_profiles_count_candidates_not_oscillators(self):
        duplicated = frequency_record("candidate", [
            ("hydrogen_bonded_oh", "donor_acceptor", 3300.0),
            ("hydrogen_bonded_oh", "donor_acceptor", 3320.0),
        ])
        profiles = build_candidate_mode_profiles([duplicated])
        classes = significant_mode_classes(
            profiles,
            config=EnvironmentAcquisitionConfig(
                minimum_sampled_environments=1, minimum_sampled_fraction=0.0,
            ),
        )
        key = "O-H:hydrogen_bonded_oh:donor_acceptor"
        self.assertEqual(classes[key]["sampled_environments"], 1)
        self.assertEqual(profiles["candidate"]["class_frequencies_cm-1"][key], 3310.0)

    def test_selector_closes_underrepresented_network_class_first(self):
        selected, report = select_mode_coverage_representatives(
            self.records, self.frequencies, self.matrix, 4,
        )
        selected_ids = {self.records[index]["candidate_id"] for index in selected}
        self.assertTrue({"candidate-3", "candidate-4", "candidate-5"}.issubset(selected_ids))
        key = "O-H:hydrogen_bonded_oh:donor_acceptor"
        self.assertEqual(report["final_independent_dft_coverage"][key], 3)
        self.assertEqual(report["status"], "coverage_targets_satisfied")

    def test_hard_budget_reports_remaining_class_deficit(self):
        _, report = select_mode_coverage_representatives(
            self.records, self.frequencies, self.matrix, 2,
        )
        key = "O-H:hydrogen_bonded_oh:donor_acceptor"
        self.assertGreater(report["remaining_class_deficits"][key], 0)
        self.assertEqual(report["status"], "orca_budget_insufficient_for_mode_class_coverage")

    def test_selector_prefers_usable_snapshot_when_coverage_is_equal(self):
        records = self.records[:2]
        frequencies = [
            frequency_record(record["candidate_id"], [
                ("non_donating_oh", "acceptor_only", 3550.0),
            ])
            for record in records
        ]
        frequencies[0]["snapshot_reliability"] = {"status": "usable_snapshot"}
        frequencies[1]["snapshot_reliability"] = {"status": "large_imaginary_curvature"}
        selected, report = select_mode_coverage_representatives(
            records, frequencies, np.zeros((2, 2)), 1,
            config=EnvironmentAcquisitionConfig(
                minimum_sampled_environments=1, minimum_sampled_fraction=0.0,
                target_dft_environments_per_class=1,
            ),
        )
        self.assertEqual(selected, [0])
        self.assertEqual(
            report["decisions"][0]["snapshot_reliability"]["status"],
            "usable_snapshot",
        )

    def test_report_records_each_orca_feedback_batch_resumably(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            (folder / "environment-sampling.json").write_text("{}")
            report = {
                "status": "coverage_targets_satisfied", "stop_reason": "selected",
                "representative_budget": 4, "remaining_class_deficits": {},
            }
            write_environment_acquisition_report(folder, report)
            record_acquisition_batch(folder, 2, {"status": "insufficient_dft_transfer_validation"})
            record_acquisition_batch(folder, 2, {"status": "validated"})
            retained = json.loads((folder / "environment-acquisition.json").read_text())
        self.assertEqual(len(retained["orca_batches"]), 1)
        self.assertEqual(retained["status"], "dft_transfer_validated")

    def test_failed_class_reorders_only_pending_representatives(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            clusters = folder / "clusters"
            clusters.mkdir()
            entries = [
                {"selected_position": index, "execution_position": index,
                 "source_xtb_candidate_id": candidate, "xyz": str(folder / f"{candidate}.xyz")}
                for index, candidate in enumerate(("done", "other", "network"), start=1)
            ]
            (clusters / "selected-conformers.json").write_text(json.dumps({"conformers": entries}))
            (folder / "environment-acquisition.json").write_text(json.dumps({
                "status": "additional_orca_representative_required", "stop_reason": "test",
                "representative_budget": 3, "remaining_class_deficits": {},
                "candidate_profiles": {
                    "done": {"mode_classes": ["common"]},
                    "other": {"mode_classes": ["common"]},
                    "network": {"mode_classes": ["network-class"]},
                },
            }))
            (folder / "frequency-transfer-validation.json").write_text(json.dumps({
                "classes": {
                    "common": {"status": "validated"},
                    "network-class": {"status": "insufficient_dft_transfer_validation"},
                },
            }))
            paths = reprioritize_pending_representatives(folder, 1)
            retained = json.loads((clusters / "selected-conformers.json").read_text())
        self.assertEqual([path.name for path in paths], ["done.xyz", "network.xyz", "other.xyz"])
        self.assertEqual(retained["conformers"][0]["source_xtb_candidate_id"], "done")


if __name__ == "__main__":
    unittest.main()
