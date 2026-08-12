import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from storca.benchmark import compare_spectra
from storca.experimental_ir import resolve_experimental_profile
from storca.ir_benchmark import load_ir_manifest
from storca.ir_contracts import (assert_partition_separation, build_experimental_condition,
                                 condition_compatibility, sha256_file)
from storca.ir_modes import NormalModeSet
from storca.mode_character import mode_character_fingerprints, match_mode_fingerprints
from storca.spectrum import _write_ir_artifacts


class IRFoundationTests(unittest.TestCase):
    def test_condition_contract_is_versioned_and_fail_closed(self):
        condition = build_experimental_condition(
            phase="liquid", measurement="atr", temperature_K=298.15,
            pressure_bar=1.0, solvent="neat", composition={"CO": 1.0},
            resolution_cm_1=4.0, apodization="gaussian", atr_crystal="diamond",
            atr_incidence_angle_degrees=45.0, sample_refractive_index=1.33,
        )
        self.assertEqual(condition.as_dict()["schema_version"], 1)
        self.assertEqual(condition.measurement.geometry, "atr")
        compatible = condition_compatibility(condition.as_dict(), condition.as_dict())
        self.assertEqual(compatible["status"], "compatible")
        unknown = condition_compatibility(
            condition.as_dict(),
            {"phase": "unspecified", "temperature_K": None,
             "measurement": {"geometry": "unspecified"}},
        )
        self.assertEqual(unknown["status"], "insufficient_metadata")
        self.assertFalse(unknown["quantitative_comparison_allowed"])
        with self.assertRaisesRegex(ValueError, "incompatible"):
            build_experimental_condition(phase="gas", measurement="atr")

    def test_reference_partitions_cannot_split_one_identity(self):
        with self.assertRaisesRegex(ValueError, "multiple benchmark partitions"):
            assert_partition_separation([
                {"id": "one", "smiles": "CO", "partition": "training"},
                {"id": "two", "smiles": "CO", "partition": "holdout"},
            ])

    def test_schema_v2_manifest_validates_hash_provenance_and_partition(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            reference = folder / "reference.csv"
            reference.write_text("wavenumber,value\n400,0\n500,1\n600,0\n")
            manifest = folder / "manifest.json"
            manifest.write_text(json.dumps({
                "schema_version": 2,
                "entries": [{
                    "id": "methanol", "smiles": "CO",
                    "reference_spectrum": reference.name,
                    "reference_sha256": sha256_file(reference),
                    "partition": "training",
                    "condition": {"phase": "liquid", "measurement": {"geometry": "atr"}},
                    "provenance": {"source": "fixture"},
                }],
            }))
            _, entries = load_ir_manifest(manifest)
            self.assertEqual(entries[0]["partition"], "training")
            entries[0]["reference_sha256"] = "bad"
            manifest.write_text(json.dumps({"schema_version": 2, "entries": entries}))
            with self.assertRaisesRegex(ValueError, "content hash changed"):
                load_ir_manifest(manifest)

    @staticmethod
    def _water_modes(reordered=False):
        coordinates = np.asarray([
            [0.0, 0.0, 0.0], [1.80, 0.0, 0.0], [-0.45, 1.74, 0.0],
        ])
        stretch = np.zeros((3, 3))
        stretch[0, 0], stretch[1, 0] = -0.4, 1.0
        bend = np.zeros((3, 3))
        bend[1, 1], bend[2, 0] = 1.0, 1.0
        vectors = np.asarray([stretch, bend])
        frequencies = np.asarray([3650.0, 1600.0])
        if reordered:
            vectors, frequencies = vectors[[1, 0]], frequencies[[1, 0]] + np.asarray([2.0, -3.0])
        return NormalModeSet(frequencies, vectors, ("O", "H", "H"), coordinates)

    def test_internal_coordinate_fingerprints_track_reordered_stretch_and_bend(self):
        reference = mode_character_fingerprints(self._water_modes())
        candidate = mode_character_fingerprints(self._water_modes(reordered=True))
        self.assertTrue(any("stretch" in item["dominant_coordinate"] for item in reference))
        self.assertTrue(any("bend" in item["dominant_coordinate"] for item in reference))
        matches = match_mode_fingerprints(reference, candidate)
        self.assertEqual([item["candidate_mode"] for item in matches], [1, 0])
        self.assertTrue(all(item["confidence"] in {"high", "medium"} for item in matches))

    def test_explicit_atom_mapping_makes_permuted_fingerprints_comparable(self):
        original = self._water_modes()
        permutation = [1, 0, 2]
        permuted = NormalModeSet(
            original.frequencies_cm_1.copy(), original.vectors[:, permutation, :],
            tuple(original.atom_labels[index] for index in permutation),
            original.coordinates_bohr[permutation],
        )
        reference = mode_character_fingerprints(original)
        candidate = mode_character_fingerprints(
            permuted, atom_mapping={0: 1, 1: 0, 2: 2},
        )
        matches = match_mode_fingerprints(reference, candidate)
        self.assertTrue(all(item["mode_character_similarity"] > 0.99 for item in matches))

    def test_short_carbon_oxygen_bond_is_labelled_as_carbonyl(self):
        coordinates = np.asarray([[0.0, 0.0, 0.0], [1.22 * 1.889726125, 0.0, 0.0]])
        vector = np.zeros((1, 2, 3))
        vector[0, 0, 0], vector[0, 1, 0] = -1.0, 1.0
        modes = NormalModeSet(np.asarray([1720.0]), vector, ("C", "O"), coordinates)
        fingerprint = mode_character_fingerprints(modes)[0]
        self.assertTrue(fingerprint["dominant_coordinate"].startswith("C=O stretch"))

    def test_intrinsic_spectrum_is_independent_of_measurement_response(self):
        records = [
            {"index": 1, "energy": -100.0, "temperature": 298.15,
             "ir_modes": [{"mode": 1, "freq": 1000.0, "intensity": 10.0},
                          {"mode": 2, "freq": 3000.0, "intensity": 10.0}]},
        ]
        with tempfile.TemporaryDirectory() as temp, patch(
            "storca.spectrum.write_spectrum_plot", side_effect=lambda path, *args, **kwargs: Path(path),
        ):
            root = Path(temp)
            (root / "first").mkdir()
            (root / "second").mkdir()
            first = _write_ir_artifacts(
                copy.deepcopy(records), root / "first", scale_factor=1.0, fwhm=15.0,
                spectrum_style="relative", max_absorbance=1.0,
                spectrum_model="experimental", practical_smiles=None,
                phase="liquid", measurement="atr", instrument_resolution=2.0,
                apodization="gaussian", residual_fwhm=3.0,
            )
            second = _write_ir_artifacts(
                copy.deepcopy(records), root / "second", scale_factor=1.0, fwhm=15.0,
                spectrum_style="relative", max_absorbance=1.0,
                spectrum_model="experimental", practical_smiles=None,
                phase="liquid", measurement="transmission", instrument_resolution=12.0,
                apodization="gaussian", residual_fwhm=3.0,
            )
            intrinsic_first = Path(first["intrinsic_spectrum_csv"]).read_bytes()
            intrinsic_second = Path(second["intrinsic_spectrum_csv"]).read_bytes()
            measured_first = Path(first["spectrum_csv"]).read_bytes()
            measured_second = Path(second["spectrum_csv"]).read_bytes()
            condition = json.loads(Path(first["experimental_condition"]).read_text())
        self.assertEqual(intrinsic_first, intrinsic_second)
        self.assertNotEqual(measured_first, measured_second)
        self.assertEqual(condition["measurement"]["geometry"], "atr")

    def test_existing_condition_contract_cannot_be_mutated_in_place(self):
        records = [{
            "index": 1, "energy": -100.0, "temperature": 298.15,
            "ir_modes": [{"mode": 1, "freq": 1000.0, "intensity": 10.0}],
        }]
        with tempfile.TemporaryDirectory() as temp, patch(
            "storca.spectrum.write_spectrum_plot", side_effect=lambda path, *args, **kwargs: Path(path),
        ):
            folder = Path(temp)
            _write_ir_artifacts(
                copy.deepcopy(records), folder, scale_factor=1.0, fwhm=15.0,
                spectrum_style="relative", max_absorbance=1.0,
                spectrum_model="experimental", practical_smiles=None,
                phase="liquid", measurement="atr", instrument_resolution=4.0,
                apodization="gaussian", residual_fwhm=3.0,
            )
            with self.assertRaisesRegex(RuntimeError, "condition contract differs"):
                _write_ir_artifacts(
                    copy.deepcopy(records), folder, scale_factor=1.0, fwhm=15.0,
                    spectrum_style="relative", max_absorbance=1.0,
                    spectrum_model="experimental", practical_smiles=None,
                    phase="liquid", measurement="atr", instrument_resolution=8.0,
                    apodization="gaussian", residual_fwhm=3.0,
                )

    def test_benchmark_reports_overlap_width_asymmetry_and_regions(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            predicted = folder / "predicted.csv"
            reference = folder / "reference.csv"
            x = np.arange(400.0, 4001.0, 10.0)
            y = np.exp(-np.square((x - 1700.0) / 120.0))
            content = "wavenumber_cm-1,relative_intensity\n" + "".join(
                f"{left},{right}\n" for left, right in zip(x, y)
            )
            predicted.write_text(content)
            reference.write_text(content)
            result = compare_spectra(predicted, reference)
        self.assertAlmostEqual(result["whole_spectrum_cosine_overlap"], 1.0)
        self.assertAlmostEqual(result["equivalent_envelope_fwhm_cm-1"]["error"], 0.0)
        self.assertAlmostEqual(result["envelope_asymmetry"]["error"], 0.0)
        self.assertEqual(len(result["integrated_regions"]), 3)


if __name__ == "__main__":
    unittest.main()
