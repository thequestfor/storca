import json
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw

from src.inputgen import create_orca_input, create_orca_irc_input, create_orca_neb_ts_input, create_orca_relaxed_scan_input
from src.inputgen import create_rmg_input
from src.orca_runner import find_orca, run_arkane, run_rmg
from src.parser import classify_chemkin_route, parse_chemkin_annotated, parse_orca_ir, parse_orca_orbitals
from storca.runs import create_run_directory, write_metadata
from storca.benchmark import compare_spectra
from storca.analysis import analyze_ir_spectra, detect_ir_peaks, match_ir_peaks
from storca.digitize import PlotCalibration, SDBS_IR_CALIBRATION, digitize_transmittance
from storca.methods import harmonic_method_profile, resolve_scale_factor
from storca.practical_ir import describe_molecule, practical_band_transform
from storca.describe import describe_smiles
from storca.ir_benchmark import evaluate_ir_manifest
from storca.calibration import scale_grid
from storca.structure import build_structure_artifacts
from storca.reactivity import frontier_reactivity_summary
from storca.enrich import enrich_pubchem
from storca.plausibility_benchmark import evaluate_plausibility_manifest
from storca.plausibility import attach_rmg_evidence, create_plausibility_dossier
from storca.spectrum import (boltzmann_weights, build_ir_spectrum, format_spectrum,
                             resume_ir_spectrum, select_goat_conformers, _completed_optimization,
                             _write_ir_artifacts)
from storca.stability import resolve_stability_configuration, run_stability_screen
from storca.conditions import build_condition_spec, normalize_target_environment_species
from storca.reachability import enrich_candidate_route
from storca.rmg_evidence import parse_final_solver_profile, parse_species_dictionary, time_to_retention_seconds
from storca.rmg_execution import assess_rmg_execution, requested_time_coverage
from storca.rmg_controller import convergence_comparison, restart_is_justified, staged_budgets
from storca.condition_ladder import build_default_ladder, humid_air_scenario, run_condition_ladder
from storca.intrinsic_initiation import enumerate_homolytic_cleavages
from storca.photolysis import evaluate_sunlight_photolysis, integrate_photolysis_rate, sunlight_contract
from storca.arkane_runner import ArkaneRouteSpec, create_arkane_input, parse_arkane_pdep_output
from storca.generated_kinetics import validate_generated_library
from storca.photophysics import (oscillator_strength_cross_sections, parse_orca_excited_state_model,
                                 parse_orca_tddft, predicted_absorption_spectrum, solar_overlap,
                                 write_tddft_input)
from storca.light_model import ComputationalLightModel, distribute_reactive_branch, energetic_accessibility
from storca.computational_light import simulate_computational_light_profiles
from storca.photo_routes import generic_photo_candidates, rank_photo_routes
from storca.decomposition_visuals import render_xyz_trajectory_gif
from storca.route_geometry import build_dissociation_endpoint, read_xyz
from storca.route_verify import RouteSpec, discover_routes, select_explanatory_route
from storca.reaction_path import run_decomposition_explanation

FIXTURES = Path(__file__).parent / "fixtures"


class FoundationTests(unittest.TestCase):
    def test_orca_input_stays_with_geometry(self):
        with tempfile.TemporaryDirectory() as temp:
            run_dir = Path(temp)
            xyz = run_dir / "input.xyz"
            xyz.write_text("2\nwater\nH 0 0 0\nH 0 0 0.7\n")
            inp = create_orca_input(xyz, charge=0, multiplicity=1, opt=True, label="opt")
            self.assertEqual(inp, run_dir / "opt.inp")
            self.assertIn("* xyz 0 1", inp.read_text())

    def test_rmg_stability_input_caps_generated_organic_dimers(self):
        with tempfile.TemporaryDirectory() as temp:
            input_file = create_rmg_input("stability", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", Path(temp))
            content = input_file.read_text()
        self.assertIn("maximumCarbonAtoms=13", content)

    def test_orca_input_accepts_profile_keywords(self):
        with tempfile.TemporaryDirectory() as temp:
            xyz = Path(temp) / "input.xyz"
            xyz.write_text("1\nhydrogen\nH 0 0 0\n")
            inp = create_orca_input(xyz, charge=0, multiplicity=1, freq=True, method_keywords=["B3LYP-3c"])
            self.assertIn("! B3LYP-3c Freq", inp.read_text())

    def test_run_directory_is_unique_and_has_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = create_run_directory(root, "ethanol run")
            second = create_run_directory(root, "ethanol run")
            self.assertNotEqual(first, second)
            metadata = write_metadata(first, charge=0)
            self.assertIn('"charge": 0', metadata.read_text())

    def test_metadata_merges_updates(self):
        with tempfile.TemporaryDirectory() as temp:
            run_dir = Path(temp)
            write_metadata(run_dir, charge=0)
            metadata = write_metadata(run_dir, temperature=298.15)
            content = metadata.read_text()
            self.assertIn('"charge": 0', content)
            self.assertIn('"temperature": 298.15', content)

    def test_benchmark_compares_overlapping_csvs(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            predicted = folder / "predicted.csv"
            reference = folder / "reference.csv"
            predicted.write_text("wavenumber,value\n400,1\n500,2\n600,3\n")
            reference.write_text("wavenumber,value\n400,2\n500,4\n600,6\n")
            metrics = compare_spectra(predicted, reference)
            self.assertEqual(metrics["points"], 3)
            self.assertAlmostEqual(metrics["pearson_correlation"], 1.0)

    def test_benchmark_recovers_known_global_band_shift(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            predicted = folder / "predicted.csv"
            reference = folder / "reference.csv"
            x_values = np.arange(400.0, 801.0, 10.0)
            predicted_values = np.exp(-((x_values - 600.0) / 35.0) ** 2)
            reference_values = np.exp(-((x_values - 620.0) / 35.0) ** 2)
            predicted.write_text("wavenumber,value\n" + "".join(f"{x},{y}\n" for x, y in zip(x_values, predicted_values)))
            reference.write_text("wavenumber,value\n" + "".join(f"{x},{y}\n" for x, y in zip(x_values, reference_values)))
            metrics = compare_spectra(predicted, reference, shift_window=30, shift_step=1)
            self.assertEqual(metrics["best_global_shift_cm-1"], 20.0)
            self.assertGreater(metrics["shift_aligned"]["pearson_correlation"], 0.999)

    def test_benchmark_handles_flat_traces_and_keeps_zero_shift(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            predicted = folder / "predicted.csv"
            reference = folder / "reference.csv"
            content = "wavenumber,value\n0,1\n1,1\n2,1\n3,1\n4,1\n"
            predicted.write_text(content)
            reference.write_text(content)
            metrics = compare_spectra(predicted, reference, shift_window=0.2, shift_step=1.0)
        self.assertIsNone(metrics["pearson_correlation"])
        self.assertEqual(metrics["rmse"], 0.0)
        self.assertEqual(metrics["best_global_shift_cm-1"], 0.0)

    def test_benchmark_rejects_duplicate_and_nonfinite_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            good = folder / "good.csv"
            bad = folder / "bad.csv"
            good.write_text("wavenumber,value\n0,0\n1,1\n2,0\n")
            bad.write_text("wavenumber,value\n0,0\n0,nan\n2,0\n")
            with self.assertRaises(ValueError):
                compare_spectra(good, bad)

    def test_peak_matching_maximizes_cardinality_before_distance(self):
        reference = [{"wavenumber_cm-1": 0.0, "relative_intensity": 1.0},
                     {"wavenumber_cm-1": 3.0, "relative_intensity": 1.0}]
        calculated = [{"wavenumber_cm-1": 1.0, "relative_intensity": 1.0},
                      {"wavenumber_cm-1": -2.0, "relative_intensity": 1.0}]
        result = match_ir_peaks(reference, calculated, window_cm_1=2.0)
        self.assertEqual(len(result["matches"]), 2)

    def test_transmittance_peak_intensity_uses_beer_lambert_absorbance(self):
        absorbance = np.array([0.0, 0.5, 0.0, 1.0, 0.0])
        transmittance = 100.0 * np.power(10.0, -absorbance)
        with tempfile.TemporaryDirectory() as temp:
            spectrum = Path(temp) / "spectrum.csv"
            spectrum.write_text("wavenumber_cm-1,transmittance_percent\n" + "".join(
                f"{index},{value}\n" for index, value in enumerate(transmittance)
            ))
            peaks = detect_ir_peaks(spectrum, min_prominence=0.1, minimum_separation_cm_1=1.0)
        intensities = sorted(item["relative_intensity"] for item in peaks)
        self.assertEqual(intensities, [0.5, 1.0])

    def test_peak_analysis_matches_without_applying_a_shift(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            calculated = folder / "calculated.csv"
            reference = folder / "reference.csv"
            x_values = np.arange(400.0, 801.0, 2.0)
            calculated_values = np.exp(-((x_values - 600.0) / 10.0) ** 2)
            reference_values = np.exp(-((x_values - 610.0) / 10.0) ** 2)
            calculated.write_text("wavenumber_cm-1,relative_intensity\n" + "".join(f"{x},{y}\n" for x, y in zip(x_values, calculated_values)))
            reference.write_text("wavenumber_cm-1,relative_intensity\n" + "".join(f"{x},{y}\n" for x, y in zip(x_values, reference_values)))
            report = analyze_ir_spectra(calculated, reference, match_window_cm_1=20)
            self.assertEqual(len(report["matches"]), 1)
            self.assertEqual(report["matches"][0]["position_error_cm-1"], -10.0)
            self.assertFalse(report["settings"]["fitting_or_frequency_shift_applied"])

    def test_harmonic_profile_scale_factor_and_override_are_explicit(self):
        self.assertEqual(resolve_scale_factor("b3lyp-def2-svp", None), (0.97, "method_profile"))
        self.assertEqual(resolve_scale_factor("b3lyp-def2-svp", 0.98), (0.98, "user_override"))

    def test_scale_grid_is_inclusive_and_stable(self):
        self.assertEqual(scale_grid(0.97, 0.974, 0.002), [0.97, 0.972, 0.974])
        self.assertEqual(scale_grid(0.94, 1.0, 0.04), [0.94, 0.98, 1.0])
        with self.assertRaises(ValueError):
            scale_grid(0.94, 1.0, 1e-10)

    def test_local_molecule_description_is_deterministic_and_offline(self):
        result = describe_smiles("CCO")
        self.assertEqual(result["canonical_smiles"], "CCO")
        self.assertEqual(result["formula"], "C2H6O")
        self.assertEqual(result["hydrogen_bond_donors"], 1)
        self.assertIn("alcohol_or_phenol", result["functional_groups"])
        self.assertFalse(result["provenance"]["network_accessed"])

    def test_structure_record_writes_local_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            result = build_structure_artifacts("CCO", Path(temp))
            self.assertEqual(result["description"]["canonical_smiles"], "CCO")
            self.assertFalse(result["network_accessed"])
            self.assertTrue(Path(result["artifacts"]["initial_xyz"]).is_file())
            self.assertTrue(Path(result["artifacts"]["structure_png"]).is_file())
            self.assertTrue(Path(result["record_json"]).is_file())

    @patch("storca.enrich._get_json")
    def test_pubchem_enrichment_is_explicitly_source_labelled(self, get_json):
        get_json.side_effect = [
            {"PropertyTable": {"Properties": [{"CID": 702, "Title": "Ethanol", "IUPACName": "ethanol", "MolecularFormula": "C2H6O", "MolecularWeight": "46.07"}]}},
            {"InformationList": {"Information": [{"Synonym": ["ethanol", "ethyl alcohol"]}]}},
        ]
        result = enrich_pubchem("CCO")
        self.assertEqual(result["pubchem"]["cid"], 702)
        self.assertEqual(result["pubchem"]["synonyms"], ["ethanol", "ethyl alcohol"])
        self.assertTrue(result["provenance"]["network_accessed"])

    def test_ir_benchmark_aggregates_registered_profile(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            reference = folder / "reference.csv"
            prediction = folder / "prediction.csv"
            reference.write_text("wavenumber,value\n400,1\n500,2\n600,3\n")
            prediction.write_text("wavenumber,value\n400,2\n500,4\n600,6\n")
            manifest = folder / "manifest.json"
            manifest.write_text(json.dumps({"schema_version": 1, "entries": [{"id": "test", "smiles": "CCO", "reference_spectrum": "reference.csv", "baseline_spectra": {"test-profile": "prediction.csv"}}]}))
            report = evaluate_ir_manifest(manifest, "test-profile")
            self.assertEqual(report["evaluated_entries"], 1)
            self.assertAlmostEqual(report["aggregate"]["mean_pearson_correlation"], 1.0)

    def test_plausibility_benchmark_prioritizes_false_reassurance(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            manifest = folder / "manifest.json"
            manifest.write_text(json.dumps({"schema_version": 1, "entries": [{
                "id": "reactive", "smiles": "OO", "charge": 0, "multiplicity": 1,
                "evidence_status": "accepted", "expected_category": "reactive_but_isolable",
                "evidence_note": "test", "evidence_source": "test source",
            }]}))
            (folder / "reactive.json").write_text(json.dumps({"assessment": {"persistence_category": "ordinary_condition_persistent"}}))
            report = evaluate_plausibility_manifest(manifest, folder)
            self.assertEqual(report["evaluated_entries"], 1)
            self.assertEqual(report["aggregate"]["false_reassurance_count"], 1)

    def test_plausibility_dossier_abstains_without_rmg_or_orca_evidence(self):
        dossier = create_plausibility_dossier("OO")
        self.assertEqual(dossier["assessment"]["persistence_category"], "insufficient_evidence")
        self.assertTrue(dossier["species"]["spin_parity_consistent"])

    def test_plausibility_dossier_keeps_abstaining_after_rmg_route(self):
        dossier = create_plausibility_dossier("OO")
        attach_rmg_evidence(dossier, {"status": "completed", "search_outcome": "candidate_decomposition_pathway_identified"})
        self.assertEqual(dossier["assessment"]["persistence_category"], "insufficient_evidence")
        self.assertIn("candidate_decomposition_pathway_identified", dossier["risk_flags"])

    def test_practical_profile_does_not_apply_unvalidated_oh_rule(self):
        features = describe_molecule("OCCO")
        frequency, width, rules = practical_band_transform(3700.0, baseline_fwhm=15.0, features=features)
        self.assertEqual(frequency, 3700.0)
        self.assertEqual(width, 15.0)
        self.assertEqual(rules, ["instrumental_baseline"])

    def test_practical_profile_does_not_apply_unvalidated_amine_rule(self):
        features = describe_molecule("CCN")
        frequency, width, rules = practical_band_transform(3400.0, baseline_fwhm=15.0, features=features)
        self.assertEqual(frequency, 3400.0)
        self.assertEqual(width, 15.0)
        self.assertEqual(rules, ["instrumental_baseline"])

    def test_digitizer_recovers_simple_trace(self):
        with tempfile.TemporaryDirectory() as temp:
            image_path = Path(temp) / "trace.png"
            image = Image.new("L", (100, 100), 255)
            draw = ImageDraw.Draw(image)
            draw.line([(12, 25), (87, 60)], fill=0, width=1)
            image.save(image_path)
            config = PlotCalibration(10, 90, 10, 90, 4000, 400, 100, 0)
            x, y, pixels = digitize_transmittance(image_path, config)
            self.assertEqual(len(x), 77)
            self.assertTrue(np.isfinite(pixels).all())
            self.assertTrue(np.isfinite(y).all())

    def test_digitizer_traces_red_curve_without_following_black_axes(self):
        with tempfile.TemporaryDirectory() as temp:
            image_path = Path(temp) / "red-trace.png"
            image = Image.new("RGB", (100, 100), "white")
            draw = ImageDraw.Draw(image)
            draw.rectangle((10, 10, 90, 90), outline="black", width=1)
            draw.line([(12, 25), (87, 60)], fill="red", width=1)
            image.save(image_path)
            config = PlotCalibration(10, 90, 10, 90, 4000, 400, 100, 0)
            _, y, pixels = digitize_transmittance(image_path, config)
            self.assertTrue(np.isfinite(y).all())
            self.assertLess(float(pixels.min()), 40.0)
            self.assertLess(float(pixels.max()), 70.0)

    def test_sdbs_calibration_uses_expanded_fingerprint_axis(self):
        anchors = dict(SDBS_IR_CALIBRATION.x_anchors)
        self.assertEqual(anchors[293], 2000)
        self.assertEqual(anchors[425], 1500)

    def test_orbital_parser_handles_missing_virtual_orbitals(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "orca.out"
            output.write_text("ORBITAL ENERGIES\nNO   OCC   E(Eh) E(eV)\n0 2.0 -0.2 -5.4\n\n")
            parsed = parse_orca_orbitals(output)
            self.assertEqual(parsed["homo_number"], 0)
            self.assertIsNone(parsed["lumo_number"])

    def test_orbital_parser_reads_orca_fixture(self):
        parsed = parse_orca_orbitals(FIXTURES / "orca_orbitals.out")
        self.assertEqual(parsed["homo_energy"], -13.6057)
        self.assertEqual(parsed["lumo_energy"], 2.7211)

    def test_frontier_reactivity_is_explicitly_qualitative(self):
        result = frontier_reactivity_summary(FIXTURES / "orca_orbitals.out")
        self.assertAlmostEqual(result["frontier_orbitals"]["gap_eV"], 16.3268)
        self.assertAlmostEqual(result["derived_frontier_proxies"]["hardness_proxy_eV"], 8.1634)
        self.assertFalse(result["provenance"]["network_accessed"])
        self.assertIn("does not identify reactive atoms", " ".join(result["limitations"]))

    def test_chemkin_parser_detects_low_barrier_decomposition(self):
        with tempfile.TemporaryDirectory() as temp:
            chemkin = Path(temp) / "chem_annotated.inp"
            chemkin.write_text((FIXTURES / "chem_annotated.inp").read_text())
            result = parse_chemkin_annotated(chemkin, label="stability")
            self.assertFalse(result["stable"])
            self.assertEqual(result["reaction_barriers"], [12.5])
            self.assertEqual(result["activation_energy_source_unit"], "KCAL/MOLE")
            self.assertEqual(result["activation_energy_unit"], "kcal/mol")

    def test_chemkin_parser_normalizes_cal_per_mole_to_kcal_per_mole(self):
        with tempfile.TemporaryDirectory() as temp:
            chemkin = Path(temp) / "chem_annotated.inp"
            chemkin.write_text(
                "REACTIONS CAL/MOLE\n"
                "stability(1)<=>product(2) 1.0E+12 0.0 12500\n"
                "END\n"
            )
            result = parse_chemkin_annotated(chemkin, label="stability")
            self.assertEqual(result["reaction_barriers"], [12.5])
            self.assertEqual(result["activation_energy_source_unit"], "CAL/MOLE")

    def test_chemkin_parser_reports_target_consuming_reverse_without_reusing_forward_barrier(self):
        with tempfile.TemporaryDirectory() as temp:
            chemkin = Path(temp) / "chem_annotated.inp"
            chemkin.write_text(
                "REACTIONS KCAL/MOLE\n"
                "[O]O(4)+[O]O(4)<=>oxygen(2)+stability(1) 1.0E+10 0.0 -3.275\n"
                "END\n"
            )
            route = parse_chemkin_annotated(chemkin, label="stability")["candidate_routes"][0]
            self.assertEqual(route["direction"], "reverse_of_chemkin_direction")
            self.assertIsNone(route["activation_energy_kcal_mol"])
            self.assertIn("stability(1)", route["reaction_equation"].split("<=>")[0])
            self.assertIn("oxygen(2)", route["reaction_equation"].split("<=>")[0])

    def test_chemkin_route_context_distinguishes_intrinsic_and_conditional_routes(self):
        self.assertEqual(
            classify_chemkin_route("stability(1)<=>product(2)")["route_context"],
            "unimolecular_decomposition",
        )
        self.assertEqual(
            classify_chemkin_route("2 stability(1)<=>dimer(2)")["route_context"],
            "self_reaction",
        )
        radical = classify_chemkin_route("[OH](2)+stability(1)<=>product(3)")
        self.assertEqual(radical["route_context"], "radical_or_impurity_initiated")
        self.assertEqual(radical["radical_co_reactants"], ["[OH]"])
        self.assertEqual(
            classify_chemkin_route("O2(2)+stability(1)<=>product(3)")["route_context"],
            "co_reactant_dependent",
        )

    @patch("storca.stability.collect_rmg_evidence")
    @patch("storca.stability.collect_orca_evidence")
    def test_combined_stability_requires_orca_and_rmg_evidence(self, collect_orca, collect_rmg):
        collect_orca.return_value = {
            "status": "completed", "local_minimum": True,
            "artifacts": {"frequency_output": "freq.out"},
        }
        collect_rmg.return_value = {
            "status": "completed", "candidate_routes": [],
            "kinetics_validation": {"status": "passed", "violation_count": 0},
            "artifacts": {"chemkin": "chemkin/chem_annotated.inp"},
        }
        with tempfile.TemporaryDirectory() as temp:
            result = run_stability_screen("CCO", Path(temp), auto_verify_routes=False)
            self.assertEqual(result["assessment"]["status"], "no_loss_detected_in_rmg_model")
            self.assertEqual(result["risk_flags"], [])
            saved = json.loads(Path(result["result_json"]).read_text())
            self.assertIn("orca_evidence", saved)
            self.assertIn("rmg_evidence", saved)

    @patch("storca.stability.collect_rmg_evidence")
    @patch("storca.stability.collect_orca_evidence")
    def test_ladder_can_reuse_exact_validated_parent_orca_evidence(self, collect_orca, collect_rmg):
        collect_rmg.return_value = {
            "status": "completed", "candidate_routes": [],
            "kinetics_validation": {"status": "passed", "violation_count": 0},
            "artifacts": {"chemkin": "chemkin/chem_annotated.inp"},
        }
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            optimized = folder / "source-opt.xyz"
            frequency = folder / "source-freq.out"
            optimized.write_text("9\nethanol\n" + "H 0 0 0\n" * 9)
            frequency.write_text("ORCA TERMINATED NORMALLY\n")
            profile = harmonic_method_profile()
            reusable = {
                "kind": "orca_geometry_and_frequency_check",
                "status": "completed",
                "canonical_smiles": "CCO",
                "charge": 0,
                "multiplicity": 1,
                "method_profile": profile,
                "local_minimum": True,
                "element_inventory_retained": True,
                "frequency_mode_coverage": {"adequate": True},
                "artifacts": {"optimized_xyz": str(optimized), "frequency_output": str(frequency)},
            }
            (folder / "next-stage").mkdir()
            result = run_stability_screen(
                "CCO", folder / "next-stage", method_profile=profile,
                precomputed_orca_evidence=reusable, auto_verify_routes=False,
            )
        collect_orca.assert_not_called()
        self.assertEqual(result["orca_evidence"]["reuse"]["status"], "validated_reuse")
        self.assertEqual(result["assessment"]["status"], "no_loss_detected_in_rmg_model")

    @patch("storca.stability.collect_rmg_evidence")
    @patch("storca.stability.collect_orca_evidence")
    def test_combined_stability_flags_intrinsic_rmg_route_even_for_orca_minimum(self, collect_orca, collect_rmg):
        collect_orca.return_value = {"status": "completed", "local_minimum": True, "artifacts": {}}
        collect_rmg.return_value = {
            "status": "completed",
            "kinetics_validation": {"status": "passed", "violation_count": 0},
            "candidate_routes": [{"reaction": "stability(1)<=>product(2)", "activation_energy_kcal_mol": 12.5,
                                  "route_context": "unimolecular_decomposition"}],
            "artifacts": {"chemkin": "chemkin/chem_annotated.inp"},
        }
        with tempfile.TemporaryDirectory() as temp:
            result = run_stability_screen("CCO", Path(temp), auto_verify_routes=False)
            self.assertEqual(result["assessment"]["status"], "candidate_intrinsic_instability")
            self.assertIn("candidate_unimolecular_decomposition", result["risk_flags"])

    @patch("storca.stability.collect_rmg_evidence")
    @patch("storca.stability.collect_orca_evidence")
    def test_incomplete_rmg_cannot_be_overridden_by_below_threshold_route(self, collect_orca, collect_rmg):
        collect_orca.return_value = {"status": "completed", "local_minimum": True, "artifacts": {}}
        collect_rmg.return_value = {
            "status": "incomplete", "candidate_routes": [], "artifacts": {"chemkin": "chemkin/chem_annotated.inp"},
            "kinetic_relevance": {"status": "reachable_but_below_loss_threshold"},
        }
        with tempfile.TemporaryDirectory() as temp:
            result = run_stability_screen("CCO", Path(temp))
        self.assertEqual(result["assessment"]["status"], "incomplete_evidence")

    @patch("storca.stability.collect_rmg_evidence")
    @patch("storca.stability.collect_orca_evidence")
    def test_combined_stability_marks_radical_route_as_condition_dependent(self, collect_orca, collect_rmg):
        collect_orca.return_value = {"status": "completed", "local_minimum": True, "artifacts": {}}
        collect_rmg.return_value = {
            "status": "completed",
            "kinetics_validation": {"status": "passed", "violation_count": 0},
            "candidate_routes": [{"reaction": "[OH](2)+stability(1)<=>product(3)", "activation_energy_kcal_mol": 6.325,
                                  "route_context": "radical_or_impurity_initiated", "co_reactants": [{"label": "[OH]", "stoichiometry": 1.0}]}],
            "artifacts": {"chemkin": "chemkin/chem_annotated.inp"},
        }
        with tempfile.TemporaryDirectory() as temp:
            result = run_stability_screen("OO", Path(temp), auto_verify_routes=False)
            self.assertEqual(result["assessment"]["status"], "condition_dependent_risk")
            self.assertIn("candidate_radical_or_impurity_initiated_pathway", result["risk_flags"])

    def test_ir_parser_reads_fixture(self):
        modes = parse_orca_ir(FIXTURES / "orca_ir.out")
        self.assertEqual(len(modes), 3)
        self.assertEqual(modes[1]["freq"], 1595.40)

    def test_ir_parser_supports_fortran_exponents_and_rejects_empty_tables(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "freq.out"
            output.write_text("IR SPECTRUM\n\n0: 1.000D+03 2.0D+00 3.0D+00\n\n")
            self.assertEqual(parse_orca_ir(output)[0]["intensity"], 3.0)
            output.write_text("IR SPECTRUM\n\nnot a data row\n")
            with self.assertRaises(RuntimeError):
                parse_orca_ir(output)

    def test_spectrum_uses_only_successful_conformers(self):
        conformers = [
            {"energy": -100.0, "temperature": 298.15, "ir_modes": [{"freq": 1000.0, "intensity": 10.0}]},
            {"energy": -99.999, "temperature": 298.15, "ir_modes": [{"freq": 1200.0, "intensity": 10.0}]},
            {"energy": None, "temperature": 298.15, "ir_modes": []},
        ]
        grid, spectrum, successful = build_ir_spectrum(conformers, frequency_min=900, frequency_max=1300, resolution=10)
        self.assertEqual(len(successful), 2)
        self.assertAlmostEqual(sum(item["weight"] for item in successful), 1.0)
        self.assertEqual(len(grid), len(spectrum))
        self.assertGreater(float(spectrum.max()), 0.0)

    def test_ir_spectrum_excludes_nonminimum_conformers(self):
        conformers = [
            {"index": 1, "energy": -100.0, "temperature": 298.15,
             "frequency_check": {"IsMinimum": True},
             "ir_modes": [{"freq": 1000.0, "intensity": 10.0}]},
            {"index": 2, "energy": -101.0, "temperature": 298.15,
             "frequency_check": {"IsMinimum": False},
             "ir_modes": [{"freq": 1200.0, "intensity": 100.0}]},
        ]
        _, _, successful = build_ir_spectrum(
            conformers, frequency_min=900, frequency_max=1300, resolution=1, scale_factor=1.0,
        )
        self.assertEqual([item["index"] for item in successful], [1])
        self.assertEqual(successful[0]["weight"], 1.0)

    def test_ir_gaussian_preserves_integrated_band_strength_across_widths(self):
        conformers = [{"energy": -100.0, "temperature": 298.15,
                       "ir_modes": [{"freq": 1000.0, "intensity": 10.0}]}]
        grid_narrow, narrow, _ = build_ir_spectrum(
            conformers, frequency_min=700, frequency_max=1300, resolution=0.1,
            fwhm=10.0, scale_factor=1.0,
        )
        grid_wide, wide, _ = build_ir_spectrum(
            conformers, frequency_min=700, frequency_max=1300, resolution=0.1,
            fwhm=60.0, scale_factor=1.0,
        )
        self.assertAlmostEqual(float(np.trapezoid(narrow, grid_narrow)), 10.0, places=5)
        self.assertAlmostEqual(float(np.trapezoid(wide, grid_wide)), 10.0, places=5)

    def test_relative_spectrum_is_unit_peak_and_rejects_nonfinite_values(self):
        values, label = format_spectrum(np.array([0.0, 2.0, 4.0]), style="relative")
        self.assertEqual(label, "relative_intensity")
        self.assertEqual(float(values.max()), 1.0)
        with self.assertRaises(ValueError):
            format_spectrum(np.array([0.0, np.nan]), style="relative")

    def test_ir_artifact_builder_propagates_weights_to_summary_records(self):
        records = [
            {"index": 1, "energy": -100.0, "temperature": 298.15,
             "ir_modes": [{"freq": 1000.0, "intensity": 10.0}]},
            {"index": 2, "energy": -99.999, "temperature": 298.15,
             "ir_modes": [{"freq": 1100.0, "intensity": 10.0}]},
        ]
        with tempfile.TemporaryDirectory() as temp, patch("storca.spectrum.write_spectrum_plot", side_effect=lambda path, *args, **kwargs: Path(path)):
            _write_ir_artifacts(records, Path(temp), scale_factor=1.0, fwhm=15.0,
                                spectrum_style="relative", max_absorbance=1.0,
                                spectrum_model="raw", practical_smiles=None)
        self.assertAlmostEqual(sum(item["weight"] for item in records), 1.0)

    def test_completed_optimization_requires_convergence_marker(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "opt.out"
            output.write_text("ORCA TERMINATED NORMALLY\n")
            self.assertFalse(_completed_optimization(output))
            output.write_text("THE OPTIMIZATION HAS CONVERGED\nORCA TERMINATED NORMALLY\n")
            self.assertTrue(_completed_optimization(output))

    @patch("storca.spectrum.run_ir_spectrum")
    def test_spectrum_resume_restores_retained_method_and_model(self, run_ir):
        run_ir.return_value = {"spectrum_csv": Path("spectrum.csv")}
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            (folder / "selected-conformers").mkdir()
            (folder / "selected-conformers" / "conf.xyz").write_text("1\ntest\nH 0 0 0\n")
            (folder / "metadata.json").write_text(json.dumps({
                "charge": 0, "multiplicity": 1, "ncores": 4, "temperature": 310.0,
                "scale_factor": 1.0, "fwhm_cm_1": 20.0, "spectrum_style": "relative",
                "spectrum_model": "practical", "smiles": "CCO",
                "harmonic_method_profile": {"orca_keywords": ["r2SCAN-3c"]},
            }))
            resume_ir_spectrum(folder)
        kwargs = run_ir.call_args.kwargs
        self.assertEqual(kwargs["method_keywords"], ["r2SCAN-3c"])
        self.assertEqual(kwargs["spectrum_model"], "practical")
        self.assertEqual(kwargs["temperature"], 310.0)

    def test_spectrum_resume_rejects_incompatible_completed_cache(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            selected = folder / "selected-conformers"
            selected.mkdir()
            xyz_text = "1\ntest\nH 0 0 0\n"
            (selected / "conf.xyz").write_text(xyz_text)
            job = folder / "conformers" / "conf-001"
            job.mkdir(parents=True)
            (job / "input.xyz").write_text(xyz_text)
            (job / "opt.out").write_text("THE OPTIMIZATION HAS CONVERGED\nORCA TERMINATED NORMALLY\n")
            (folder / "metadata.json").write_text(json.dumps({
                "charge": 0, "multiplicity": 1,
                "harmonic_method_profile": {"orca_keywords": ["B3LYP", "def2-SVP"]},
            }))
            with self.assertRaisesRegex(RuntimeError, "different or unverified"):
                resume_ir_spectrum(folder, charge=1)

    def test_boltzmann_weights_reject_invalid_energy(self):
        with self.assertRaises(ValueError):
            boltzmann_weights([float("inf")], 298.15)
        self.assertEqual(boltzmann_weights([-1.0, -1.0], 298.15, degeneracies=[1, 3]), [0.25, 0.75])

    def test_goat_selection_uses_cumulative_population(self):
        selected = select_goat_conformers([
            {"idx": 0, "weight": 0.60, "cum_weight": 0.60},
            {"idx": 1, "weight": 0.35, "cum_weight": 0.95},
            {"idx": 2, "weight": 0.05, "cum_weight": 1.00},
        ], population_cutoff=0.95, max_conformers=10)
        self.assertEqual([item["idx"] for item in selected], [0, 1])

    def test_transmittance_display_turns_absorptions_into_downward_bands(self):
        values, label = format_spectrum(np.array([0.0, 5.0]), style="transmittance", max_absorbance=1.0)
        self.assertEqual(label, "transmittance_percent")
        self.assertAlmostEqual(values[0], 100.0)
        self.assertAlmostEqual(values[1], 10.0)

    def test_orca_path_can_be_configured_per_machine(self):
        with tempfile.TemporaryDirectory() as temp:
            executable = Path(temp) / "orca"
            executable.write_text("#!/bin/sh\n")
            executable.chmod(0o755)
            with patch.dict("os.environ", {"STORCA_ORCA_BIN": str(executable)}):
                self.assertEqual(find_orca(), str(executable.resolve()))

    def test_rmg_input_records_requested_conditions(self):
        with tempfile.TemporaryDirectory() as temp:
            input_file = create_rmg_input(
                "stability", "CCO", Path(temp), temperature=350.0, pressure=2.0
            )
            content = input_file.read_text()
            self.assertIn("SMILES('CCO')", content)
            self.assertIn("temperature=(350.0, 'K')", content)
            self.assertIn("pressure=(2.0, 'bar')", content)

    def test_rmg_input_can_declare_an_inert_atmosphere(self):
        with tempfile.TemporaryDirectory() as temp:
            input_file = create_rmg_input(
                "stability", "CCO", Path(temp),
                additional_species=[{"label": "nitrogen", "smiles": "N#N", "reactive": False}],
                initial_mole_fractions={"stability": 0.01, "nitrogen": 0.99},
            )
            content = input_file.read_text()
            self.assertIn("label='nitrogen'", content)
            self.assertIn("reactive=False", content)
            self.assertIn("'stability': 0.01", content)
            self.assertIn("'nitrogen': 0.99", content)

    def test_rmg_input_merges_target_identical_to_environment_species(self):
        with tempfile.TemporaryDirectory() as temp:
            input_file = create_rmg_input(
                "stability", "O", Path(temp),
                additional_species=[
                    {"label": "water", "smiles": "O", "reactive": True},
                    {"label": "nitrogen", "smiles": "N#N", "reactive": False},
                ],
                initial_mole_fractions={"stability": 0.01, "water": 0.02, "nitrogen": 0.97},
            )
            content = input_file.read_text()
        self.assertEqual(content.count("structure=SMILES('O')"), 1)
        self.assertNotIn("label='water'", content)
        self.assertIn("'stability': 0.03", content)

    def test_rmg_input_loads_an_explicit_local_reaction_library(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            library = folder / "kinetics"
            library.mkdir()
            (library / "reactions.py").write_text("# library\n")
            (library / "dictionary.txt").write_text("# dictionary\n")
            input_file = create_rmg_input("stability", "CCO", folder / "run", reaction_libraries=[library])
            self.assertIn(str(library.resolve()), input_file.read_text())

    def test_stability_configuration_resolves_named_profiles_and_rejects_unknown_ones(self):
        scenario, tier = resolve_stability_configuration("ambient-air-gas-screen", "review-screen")
        self.assertEqual(scenario["atmosphere"], "dry air (oxygen/nitrogen)")
        self.assertEqual(tier["max_iterations"], 300)
        with self.assertRaises(ValueError):
            resolve_stability_configuration("humid-liquid-storage", "quick-screen")

    def test_generated_radical_route_requires_initiation_when_profile_concentration_is_negligible(self):
        scenario, _ = resolve_stability_configuration("ambient-air-gas-screen", "quick-screen")
        conditions = build_condition_spec(scenario, temperature_K=298.15, pressure_bar=1.0)
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            dictionary = folder / "species_dictionary.txt"
            dictionary.write_text(
                "stability(1)\n1 C u0 p0 c0\n\n"
                "[O]O(15)\nmultiplicity 2\n1 O u1 p2 c0\n"
            )
            solver = folder / "solver"
            solver.mkdir()
            (solver / "simulation_1_1.csv").write_text(
                "Time (s),stability(1),[O]O(15)\n0,0.01,0\n86400,0.01,1e-25\n"
            )
            profile = parse_final_solver_profile(solver, "stability")
            route = enrich_candidate_route(
                {"co_reactants": [{"label": "[O]O", "stoichiometry": 1.0}]},
                species_dictionary=parse_species_dictionary(dictionary), profile=profile, conditions=conditions,
            )
            self.assertEqual(route["initiation_status"], "requires_generated_intermediate")
            self.assertEqual(route["kinetic_relevance"], "no_credible_initiation_demonstrated")

    def test_solver_profile_interpolates_time_to_retention(self):
        profile = {"target_time_series": [
            {"time_seconds": 0.0, "fraction_remaining": 1.0},
            {"time_seconds": 10.0, "fraction_remaining": 0.90},
        ]}
        self.assertAlmostEqual(time_to_retention_seconds(profile, 0.95), 5.0)
        self.assertIsNone(time_to_retention_seconds(profile, 0.85))

    def test_generic_intrinsic_discovery_finds_h2o2_oo_homolysis_without_an_alert_rule(self):
        result = enumerate_homolytic_cleavages("OO")
        self.assertEqual(len(result["candidates"]), 1)
        self.assertEqual(result["candidates"][0]["fragment_radical_smiles"], ["[OH]", "[OH]"])

    def test_humid_air_composition_records_water_vapour_and_sums_to_one(self):
        scenario = humid_air_scenario(temperature_K=298.15, pressure_bar=1.0, relative_humidity=0.5)
        self.assertGreater(scenario["initial_mole_fractions"]["water"], 0.0)
        self.assertAlmostEqual(sum(scenario["initial_mole_fractions"].values()), 1.0)

    def test_water_target_merges_humidity_into_one_target_pool(self):
        scenario = humid_air_scenario(temperature_K=298.15, pressure_bar=1.0, relative_humidity=0.5)
        normalized = normalize_target_environment_species("O", scenario)
        condition = build_condition_spec(normalized, temperature_K=298.15, pressure_bar=1.0)
        self.assertNotIn("water", normalized["initial_mole_fractions"])
        self.assertGreater(condition.target_mole_fraction, 0.01)
        self.assertEqual(condition.as_dict()["target_identity_merges"][0]["environment_label"], "water")

    def test_sunlight_requires_photochemical_evidence_before_claiming_a_lifetime(self):
        result = evaluate_sunlight_photolysis(
            condition={"retention_fraction": 0.95, "light_model": sunlight_contract()},
        )
        self.assertEqual(result["status"], "photochemical_evidence_incomplete")
        self.assertIsNone(result["estimated_time_to_retention_seconds"])

    def test_no_reactive_photon_profiles_are_completed_negative_evidence(self):
        light = {
            "status": "completed",
            "smiles": "CCC",
            "model_contract": {"name": "test-light"},
            "profiles": {
                "low": {"reaction_libraries": []},
                "central": {"reaction_libraries": []},
                "high": {"reaction_libraries": []},
            },
        }
        with tempfile.TemporaryDirectory() as temp:
            result = simulate_computational_light_profiles(light, Path(temp), rmg_env="rmg_env")
        self.assertEqual(result["status"], "completed_no_reactive_photon_channel")
        self.assertEqual(result["outcome"], "no_modeled_sunlight_loss")

    @patch("storca.computational_light.simulate_computational_light_profiles")
    @patch("storca.computational_light.run_computational_light_model")
    def test_completed_no_channel_sunlight_stages_remain_photochemically_incomplete(self, run_light, profiles):
        run_light.return_value = {
            "status": "completed",
            "smiles": "CCC",
            "model_contract": {"name": "test-light"},
            "profiles": {"low": {}, "central": {}, "high": {}},
        }
        profiles.return_value = {
            "status": "completed_no_reactive_photon_channel",
            "profiles": {
                "low": {"status": "no_reactive_photon_channel"},
                "central": {"status": "no_reactive_photon_channel"},
                "high": {"status": "no_reactive_photon_channel"},
            },
        }

        def stage_runner(*args, **kwargs):
            return {
                "assessment": {"status": "no_loss_detected_in_rmg_model", "reason": "test"},
                "orca_evidence": {"status": "completed"},
                "rmg_evidence": {"status": "completed"},
                "kinetic_lifetime": {"estimated_time_to_retention_seconds": None},
            }

        with tempfile.TemporaryDirectory() as temp:
            result = run_condition_ladder(
                "CCC", Path(temp), stage_runner=stage_runner,
                computational_light_spectrum=Path(temp) / "sun.csv",
            )
        self.assertEqual(result["verdict"], "no_verified_t95_with_incomplete_evidence")
        self.assertEqual(result["stages"][4]["status"], "completed_with_incomplete_evidence")

    def test_spectral_photolysis_integral_uses_absorption_and_quantum_yield(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            spectrum = folder / "spectrum.csv"
            absorption = folder / "absorption.csv"
            spectrum.write_text("wavelength_nm,irradiance_W_m2_nm\n300,1\n400,1\n")
            absorption.write_text("wavelength_nm,absorption_cross_section_cm2_molecule,quantum_yield\n300,1e-20,1\n400,1e-20,1\n")
            result = integrate_photolysis_rate(spectrum, absorption)
        self.assertGreater(result["photolysis_rate_constant_s-1"], 0.0)
        self.assertEqual(len(result["wavelength_contributions"]), 2)

    def test_photolysis_integral_rejects_singleton_duplicate_and_nonfinite_grids(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            spectrum = folder / "spectrum.csv"
            absorption = folder / "absorption.csv"
            spectrum.write_text("wavelength_nm,irradiance_W_m2_nm\n300,1\n")
            absorption.write_text("wavelength_nm,absorption_cross_section_cm2_molecule,quantum_yield\n300,1e-18,0.1\n")
            with self.assertRaises(ValueError):
                integrate_photolysis_rate(spectrum, absorption)
            spectrum.write_text("wavelength_nm,irradiance_W_m2_nm\n300,1\n300,2\n")
            absorption.write_text("wavelength_nm,absorption_cross_section_cm2_molecule,quantum_yield\n300,1e-18,0.1\n400,1e-18,0.1\n")
            with self.assertRaises(ValueError):
                integrate_photolysis_rate(spectrum, absorption)
            spectrum.write_text("wavelength_nm,irradiance_W_m2_nm\n300,nan\n400,1\n")
            with self.assertRaises(ValueError):
                integrate_photolysis_rate(spectrum, absorption)

    def test_spectral_photolysis_evidence_can_complete_a_sunlight_stage(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            spectrum = folder / "spectrum.csv"
            absorption = folder / "absorption.csv"
            evidence = folder / "evidence.json"
            spectrum.write_text("wavelength_nm,irradiance_W_m2_nm\n300,1\n400,1\n")
            absorption.write_text("wavelength_nm,absorption_cross_section_cm2_molecule,quantum_yield\n300,1e-20,1\n400,1e-20,1\n")
            evidence.write_text(json.dumps({"spectrum_csv": "spectrum.csv", "absorption_csv": "absorption.csv",
                                             "photoproducts": [{"label": "radical", "smiles": "[OH]"}]}))
            result = evaluate_sunlight_photolysis(
                condition={"retention_fraction": 0.95, "light_model": sunlight_contract()}, photolysis_evidence=evidence,
            )
        self.assertEqual(result["status"], "completed")
        self.assertGreater(result["estimated_time_to_retention_seconds"], 0.0)

    def test_arkane_input_retains_orca_artifact_paths_and_pressure_grid(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            outputs = []
            for name in ("reactant.out", "product_a.out", "product_b.out", "ts.out"):
                path = folder / name
                path.write_text("ORCA TERMINATED NORMALLY\n")
                outputs.append(path)
            spec = ArkaneRouteSpec(
                label="homolysis", reactant_label="parent", reactant_smiles="OO", reactant_orca_output=outputs[0],
                product_labels=("oh_a", "oh_b"), product_smiles=("[OH]", "[OH]"),
                product_orca_outputs=(outputs[1], outputs[2]), transition_state_orca_output=outputs[3],
            )
            with self.assertRaisesRegex(ValueError, "validated profile"):
                create_arkane_input(spec, folder / "arkane")

    def test_arkane_pdep_parser_requires_pdep_blocks(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "output.py"
            output.write_text("pdepreaction('route')\n")
            self.assertEqual(parse_arkane_pdep_output(output)["status"], "completed")

    def test_generated_library_rejects_checksum_tampering(self):
        with tempfile.TemporaryDirectory() as temp:
            library = Path(temp) / "kinetics"
            library.mkdir()
            reactions, dictionary = library / "reactions.py", library / "dictionary.txt"
            reactions.write_text("reaction\n")
            dictionary.write_text("dictionary\n")
            import hashlib
            checksum = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
            (library / "generated-kinetics.json").write_text(json.dumps({
                "verification_status": "verified_spectral_rate",
                "checksums": {"reactions.py": checksum(reactions), "dictionary.txt": checksum(dictionary)},
            }))
            validate_generated_library(library)
            reactions.write_text("changed\n")
            with self.assertRaises(ValueError):
                validate_generated_library(library)

    def test_condition_ladder_stops_only_after_modelled_t95(self):
        def stage_runner(smiles, stage_dir, **kwargs):
            if kwargs["scenario"] == "ambient-air-gas-screen":
                return {"assessment": {"status": "orca_verified_condition_dependent_instability"},
                        "orca_evidence": {"status": "completed"}, "rmg_evidence": {"status": "completed"},
                        "condition_contract": {"scenario": "ambient-air-gas-screen", "retention_fraction": 0.95},
                        "verification_summary": {"status": "verified_t95_converged",
                                                 "estimated_time_to_retention_seconds": 10.0,
                                                 "condition_contract_match": True,
                                                 "target_loss_fraction": 0.05,
                                                 "verified_flux_coverage": 1.0,
                                                 "unresolved_loss_upper_bound": 0.0},
                        "kinetic_lifetime": {"estimated_time_to_retention_seconds": 10.0}}
            return {"assessment": {"status": "orca_verified_condition_dependent_instability"}, "orca_evidence": {"status": "completed"}, "rmg_evidence": {"status": "completed"},
                    "kinetic_lifetime": {"estimated_time_to_retention_seconds": None}}
        with tempfile.TemporaryDirectory() as temp:
            result = run_condition_ladder("CCO", Path(temp), stage_runner=stage_runner)
        self.assertEqual(result["verdict"], "condition_scoped_t95_identified")
        self.assertEqual(result["first_t95_stage"], "dark-dry-air")
        self.assertEqual(result["stages"][3]["status"], "not_run_after_verified_t95")

    def test_condition_ladder_serializes_stage_artifact_paths(self):
        def stage_runner(smiles, stage_dir, **kwargs):
            return {"kinetic_lifetime": {"estimated_time_to_retention_seconds": None},
                    "artifact": Path(stage_dir) / "artifact.json"}
        with tempfile.TemporaryDirectory() as temp:
            result = run_condition_ladder("CCO", Path(temp), stage_runner=stage_runner)
            saved = json.loads(Path(result["result_json"]).read_text())
        self.assertEqual(saved["stages"][1]["artifact"], str(Path(temp) / "ladder" / "02-dark-dry-inert" / "artifact.json"))
        self.assertEqual(saved["stages"][1]["status"], "completed_with_incomplete_evidence")

    @patch("storca.computational_light.run_computational_light_model")
    def test_condition_ladder_forwards_electronic_contract_to_light_screen(self, run_light):
        run_light.return_value = {"status": "failed", "failure_reason": "test"}
        def stage_runner(*args, **kwargs):
            return {"orca_evidence": {"status": "completed"},
                    "rmg_evidence": {"status": "completed"},
                    "kinetic_lifetime": {"estimated_time_to_retention_seconds": None}}
        with tempfile.TemporaryDirectory() as temp:
            run_condition_ladder(
                "[NH3+]O", Path(temp) / "run", stage_runner=stage_runner,
                computational_light_spectrum=Path(temp) / "sun.csv",
                charge=1, multiplicity=2, method_keywords=["r2SCAN-3c"],
                light_nroots=24, ncores=4,
            )
        kwargs = run_light.call_args.kwargs
        self.assertEqual(kwargs["charge"], 1)
        self.assertEqual(kwargs["multiplicity"], 2)
        self.assertEqual(kwargs["method_keywords"], ["r2SCAN-3c"])
        self.assertEqual(kwargs["nroots"], 24)

    @patch("storca.generated_kinetics.write_photolysis_library")
    def test_sunlight_j_does_not_stop_ladder_without_modelled_loss(self, write_library):
        def stage_runner(smiles, stage_dir, **kwargs):
            return {"rmg_evidence": {"status": "completed"}, "orca_evidence": {"status": "completed"},
                    "kinetic_lifetime": {"estimated_time_to_retention_seconds": None}}
        write_library.return_value = {"library": "/tmp/kinetics"}
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            spectrum, absorption, evidence = folder / "spectrum.csv", folder / "absorption.csv", folder / "evidence.json"
            spectrum.write_text("wavelength_nm,irradiance_W_m2_nm\n300,1\n400,1\n")
            absorption.write_text("wavelength_nm,absorption_cross_section_cm2_molecule,quantum_yield\n300,1e-20,1\n400,1e-20,1\n")
            evidence.write_text(json.dumps({"spectrum_csv": str(spectrum), "absorption_csv": str(absorption),
                                             "photoproducts": [{"label": "radical", "smiles": "[OH]"}]}))
            result = run_condition_ladder("CCO", folder / "run", stage_runner=stage_runner,
                                          photolysis_evidence=evidence)
        self.assertEqual(result["verdict"], "no_verified_t95_within_completed_models")
        self.assertIsNone(result["first_t95_stage"])

    @patch("storca.generated_kinetics.write_photolysis_library")
    def test_condition_ladder_stops_at_completed_sunlight_t95(self, write_library):
        def stage_runner(smiles, stage_dir, **kwargs):
            verified = kwargs.get("light_condition") == "sunlight"
            return {"orca_evidence": {"status": "completed"}, "rmg_evidence": {"status": "completed"},
                    "assessment": {"status": "orca_verified_condition_dependent_instability"},
                    "condition_contract": {"scenario": "ambient-air-gas-screen", "retention_fraction": 0.95},
                    "verification_summary": ({"status": "verified_t95_converged",
                                              "estimated_time_to_retention_seconds": 10.0,
                                              "condition_contract_match": True,
                                              "target_loss_fraction": 0.05,
                                              "verified_flux_coverage": 1.0,
                                              "unresolved_loss_upper_bound": 0.0}
                                             if verified else {"status": "not_required"}),
                    "kinetic_lifetime": {"estimated_time_to_retention_seconds": 10.0 if verified else None}}
        write_library.return_value = {"library": "/tmp/kinetics"}
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            spectrum = folder / "spectrum.csv"
            absorption = folder / "absorption.csv"
            evidence = folder / "evidence.json"
            spectrum.write_text("wavelength_nm,irradiance_W_m2_nm\n300,1\n400,1\n")
            absorption.write_text("wavelength_nm,absorption_cross_section_cm2_molecule,quantum_yield\n300,1e-20,1\n400,1e-20,1\n")
            evidence.write_text(json.dumps({"spectrum_csv": str(spectrum), "absorption_csv": str(absorption),
                                             "photoproducts": [{"label": "radical", "smiles": "[OH]"}]}))
            result = run_condition_ladder("CCO", folder / "run", stage_runner=stage_runner,
                                          photolysis_evidence=evidence)
        self.assertEqual(result["first_t95_stage"], "sunlight-dry-air")
        self.assertEqual(result["stages"][5]["status"], "not_run_after_verified_t95")

    def test_nno_ladder_does_not_turn_model_no_loss_into_stability(self):
        fixture = json.loads((FIXTURES / "nno_ladder_regression.json").read_text())
        def stage_runner(smiles, stage_dir, **kwargs):
            if kwargs.get("scenario") == "dry-inert-gas-screen":
                return fixture["dark_inert"]
            return fixture["dark_air"]
        with tempfile.TemporaryDirectory() as temp:
            result = run_condition_ladder(fixture["smiles"], Path(temp), stage_runner=stage_runner)
        inert = result["stages"][1]
        self.assertEqual(inert["assessment"]["status"], "no_loss_detected_in_rmg_model")
        self.assertEqual(inert["status"], "completed_model_screen_no_loss")
        self.assertEqual(result["verdict"], "no_verified_t95_with_incomplete_evidence")
        gaps = {(item["stage_id"], item["code"]) for item in result["evidence_summary"]["missing_evidence"]}
        self.assertIn(("dark-low-pressure", "initiation_ts_kinetics_required"), gaps)
        self.assertIn(("dark-low-pressure", "screening_t95_unverified"), gaps)
        self.assertIn(("dark-dry-air", "initiation_ts_kinetics_required"), gaps)
        self.assertIn(("dark-dry-air", "screening_t95_unverified"), gaps)
        self.assertNotIn("stable", json.dumps(inert["assessment"]))

    def test_orca_assessment_and_number_do_not_stop_without_condition_contract(self):
        calls = []
        def stage_runner(smiles, stage_dir, **kwargs):
            calls.append(stage_dir)
            return {"status": "completed", "orca_evidence": {"status": "completed"},
                    "rmg_evidence": {"status": "completed"},
                    "assessment": {"status": "orca_verified_condition_dependent_instability"},
                    "kinetic_lifetime": {"estimated_time_to_retention_seconds": 12.0}}
        with tempfile.TemporaryDirectory() as temp:
            result = run_condition_ladder("NNO", Path(temp), stage_runner=stage_runner)
        self.assertIsNone(result["first_t95_stage"])
        self.assertEqual(len(calls), 4)
        self.assertTrue(any(item["code"] == "reported_t95_not_verified_converged"
                            for item in result["evidence_summary"]["missing_evidence"]))

    def test_verified_below_threshold_is_distinct_from_rmg_no_loss(self):
        contract = {"scenario": "dry-inert-gas-screen", "retention_fraction": 0.95}
        stage = {"status": "completed", "condition_contract": contract,
                 "verification_summary": {"status": "verified_below_loss_threshold",
                                          "condition_contract_match": True,
                                          "estimated_time_to_retention_seconds": None,
                                          "target_loss_fraction": 0.01,
                                          "verified_flux_coverage": 0.99,
                                          "unresolved_loss_upper_bound": 0.005}}
        from storca.condition_ladder import _normalize_stage_result
        result = _normalize_stage_result(stage)
        self.assertEqual(result["status"], "completed_verified_below_loss_threshold")
        self.assertEqual(result["assessment"]["status"], "verified_route_below_loss_threshold")
        self.assertEqual(result["missing_evidence"], [])


    @patch("src.orca_runner.subprocess.run")
    def test_rmg_runner_applies_bounded_execution_options(self, run):
        with tempfile.TemporaryDirectory() as temp:
            input_file = Path(temp) / "input.py"
            input_file.write_text("# test\n")
            run.return_value = __import__("subprocess").CompletedProcess([], 0, stdout="ok", stderr="")
            result = run_rmg(
                input_file, rmg_env="rmg_env", walltime="00:00:10:00",
                max_processes=2, max_iterations=100,
            )
            command = run.call_args.args[0]
            self.assertIn("-t", command)
            self.assertIn("00:00:10:00", command)
            self.assertIn("-n", command)
            self.assertEqual(result["stdout"], "ok")

    def test_rmg_execution_rejects_saved_chemkin_after_resource_termination(self):
        with tempfile.TemporaryDirectory() as temp:
            log = Path(temp) / "RMG.log"
            log.write_text("Saving current model core to Chemkin file...\nMODEL GENERATION TERMINATED\n"
                           "The output model may be incomplete.\nRMG execution terminated\n")
            result = assess_rmg_execution(log)
        self.assertEqual(result["status"], "incomplete")
        self.assertEqual(result["termination_reason"], "resource_or_interrupt_termination")

    def test_time_coverage_does_not_extrapolate_a_short_solver_profile(self):
        result = requested_time_coverage({"end_time_seconds": 0.03}, 365 * 24 * 3600)
        self.assertFalse(result["complete"])
        self.assertLess(result["fraction"], 1e-8)

    def test_rmg_controller_restarts_only_a_real_budget_exhaustion_with_seed(self):
        evidence = {"status": "incomplete", "execution_assessment": {"termination_reason": "resource_or_interrupt_termination"},
                    "artifacts": {"seed": "/tmp/seed"}}
        self.assertTrue(restart_is_justified(evidence))
        evidence["execution_assessment"]["termination_reason"] = "missing_rmg_log"
        self.assertFalse(restart_is_justified(evidence))

    def test_rmg_controller_expands_only_declared_retry_budgets(self):
        budgets = staged_budgets({"walltime": "00:00:10:00", "max_processes": 4,
                                  "max_iterations": 100, "max_edge_species": 250}, 2)
        self.assertEqual([item.walltime for item in budgets], ["00:00:10:00", "00:00:20:00"])
        self.assertEqual(budgets[1].max_edge_species, 500)
        comparison = convergence_comparison([{"candidate_routes": []}, {"candidate_routes": []}])
        self.assertEqual(comparison["new_candidate_route_count"], 0)

    @patch("src.orca_runner.subprocess.run")
    def test_arkane_runner_uses_the_rmg_environment(self, run):
        with tempfile.TemporaryDirectory() as temp:
            input_file = Path(temp) / "input.py"
            input_file.write_text("# test\n")
            run.return_value = __import__("subprocess").CompletedProcess([], 0, stdout="ok", stderr="")
            result = run_arkane(input_file, rmg_env="rmg_env")
        self.assertEqual(run.call_args.args[0][:5], ["conda", "run", "-n", "rmg_env", "Arkane.py"])
        self.assertEqual(result["stdout"], "ok")

    def test_tddft_input_and_absorption_parser(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            xyz = folder / "input.xyz"
            xyz.write_text("2\ntest\nH 0 0 0\nH 0 0 0.7\n")
            inp = write_tddft_input(xyz, nroots=12)
            self.assertIn("%tddft", inp.read_text())
            self.assertIn("nroots 12", inp.read_text())
            self.assertIn("tda false", inp.read_text())
            out = folder / "tddft.out"
            out.write_text("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n\n"
                           " 1  25000.0  400.0  0.120\n\n")
            transitions = parse_orca_tddft(out)
        self.assertEqual(transitions[0]["state"], 1)
        self.assertAlmostEqual(transitions[0]["energy_eV"], 3.0996, places=3)
        self.assertEqual(transitions[0]["oscillator_strength_gauge"], "electric_length")

    def test_tddft_parser_uses_electric_not_velocity_gauge(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "tddft.out"
            output.write_text(
                "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n\n"
                "0-1A -> 1-1A 4.0 32262.0 310.0 0.100\n\n"
                "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\n\n"
                "0-1A -> 1-1A 4.0 32262.0 310.0 0.900\n\n"
            )
            transitions = parse_orca_tddft(output)
        self.assertEqual(len(transitions), 1)
        self.assertEqual(transitions[0]["oscillator_strength"], 0.1)
        self.assertEqual(transitions[0]["initial_state"], "0-1A")

    def test_orca_excited_state_model_distinguishes_tda(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "tddft.out"
            output.write_text("Tamm-Dancoff approximation     ... operative\n")
            self.assertEqual(parse_orca_excited_state_model(output), "tda_tddft")
            output.write_text("Tamm-Dancoff approximation     ... not operative\n")
            self.assertEqual(parse_orca_excited_state_model(output), "full_lr_tddft")
            output.write_text("Tamm-Dancoff approximation     ... deactivated\n")
            self.assertEqual(parse_orca_excited_state_model(output), "full_lr_tddft")

    def test_relative_solar_overlap_requires_declared_spectrum(self):
        spectrum = predicted_absorption_spectrum([{"state": 1, "wavelength_nm": 400.0, "oscillator_strength": 1.0}],
                                                 minimum_nm=300, maximum_nm=500, step_nm=100)
        with tempfile.TemporaryDirectory() as temp:
            sunlight = Path(temp) / "sun.csv"
            sunlight.write_text("wavelength_nm,irradiance_W_m2_nm\n300,1\n400,2\n500,1\n")
            result = solar_overlap(spectrum, sunlight)
        self.assertGreater(result["relative_solar_overlap"], 0.0)

    def test_relative_solar_overlap_interpolates_different_grids(self):
        spectrum = predicted_absorption_spectrum(
            [{"wavelength_nm": 400.0, "oscillator_strength": 1.0}],
            minimum_nm=300, maximum_nm=500, step_nm=100,
        )
        with tempfile.TemporaryDirectory() as temp:
            sunlight = Path(temp) / "sun.csv"
            sunlight.write_text(
                "wavelength_nm,irradiance_W_m2_nm\n325,1\n375,2\n425,2\n475,1\n"
            )
            result = solar_overlap(spectrum, sunlight)
        self.assertEqual(len(result["wavelength_contributions"]), 4)
        self.assertEqual(result["grid_alignment"], "linear_absorption_interpolation_to_sunlight_grid")
        self.assertGreater(result["relative_solar_overlap"], 0.0)

    def test_computational_light_model_uses_ordered_branch_profiles(self):
        model = ComputationalLightModel()
        self.assertEqual(model.profiles(), {"low": 0.001, "central": 0.01, "high": 0.1})
        self.assertGreater(energetic_accessibility(4.0, 3.0), energetic_accessibility(2.0, 3.0))
        distributed = distribute_reactive_branch([{"accessibility": 1.0}, {"accessibility": 3.0}], 0.1)
        self.assertAlmostEqual(sum(route["quantum_yield"] for route in distributed), 0.1)

    def test_photo_route_ranking_is_generic_and_energy_scored(self):
        candidates = generic_photo_candidates("OO", maximum_routes=3)
        ranked = rank_photo_routes(candidates,
                                   [{"state": 1, "wavelength_nm": 300.0, "energy_eV": 4.13, "oscillator_strength": 0.2}],
                                   {candidates[0]["route_id"]: 3.0})
        self.assertEqual(ranked[0]["kind"], "homolysis")
        self.assertEqual(ranked[0]["status"], "computed")
        self.assertGreater(ranked[0]["accessibility"], 0.5)

    def test_oscillator_strength_cross_sections_are_positive(self):
        cross_sections = oscillator_strength_cross_sections(
            [{"wavelength_nm": 300.0, "oscillator_strength": 0.1}], wavelengths_nm=[280.0, 300.0, 320.0])
        self.assertEqual(len(cross_sections), 3)
        self.assertGreater(cross_sections[1]["absorption_cross_section_cm2_molecule"], cross_sections[0]["absorption_cross_section_cm2_molecule"])

    def test_oscillator_strength_cross_section_conserves_sum_rule_area(self):
        speed_of_light = 299792458.0
        center_wavelength_nm = 300.0
        fwhm_nm = 20.0
        center_frequency = speed_of_light / (center_wavelength_nm * 1e-9)
        frequency_fwhm = speed_of_light * (fwhm_nm * 1e-9) / (center_wavelength_nm * 1e-9) ** 2
        sigma = frequency_fwhm / (2 * np.sqrt(2 * np.log(2)))
        frequencies = np.linspace(center_frequency - 7 * sigma, center_frequency + 7 * sigma, 30001)
        wavelengths = speed_of_light / frequencies * 1e9
        points = oscillator_strength_cross_sections(
            [{"wavelength_nm": center_wavelength_nm, "oscillator_strength": 0.1}],
            fwhm_nm=fwhm_nm, wavelengths_nm=list(wavelengths),
        )
        area = np.trapezoid(
            np.asarray([item["absorption_cross_section_cm2_molecule"] for item in points]) * 1e-4,
            frequencies,
        )
        expected = (1.602176634e-19 ** 2 / (4 * 8.8541878128e-12 * 9.1093837015e-31 * speed_of_light)) * 0.1
        self.assertAlmostEqual(float(area / expected), 1.0, places=5)

    def test_inaccessible_photo_routes_do_not_consume_full_branch_prior(self):
        distributed = distribute_reactive_branch(
            [{"accessibility": 1e-6}, {"accessibility": 2e-6}], 0.1,
        )
        self.assertAlmostEqual(sum(route["quantum_yield"] for route in distributed), 3e-7)

    def test_photo_route_score_changes_smoothly_with_oscillator_strength(self):
        candidate = [{"route_id": "r1", "kind": "homolysis"}]
        base_states = [
            {"state": 1, "wavelength_nm": 300.0, "energy_eV": 4.0, "oscillator_strength": 0.051},
            {"state": 2, "wavelength_nm": 500.0, "energy_eV": 2.0, "oscillator_strength": 0.050},
        ]
        changed_states = [{**base_states[0], "oscillator_strength": 0.053}, base_states[1]]
        first = rank_photo_routes(candidate, base_states, {"r1": 3.0})[0]["accessibility"]
        second = rank_photo_routes(candidate, changed_states, {"r1": 3.0})[0]["accessibility"]
        self.assertLess(abs(second - first), 0.02)

    def test_xyz_trajectory_renderer_writes_a_gif(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            xyz = folder / "trajectory.xyz"
            xyz.write_text("2\nframe 1\nH 0 0 0\nH 0 0 0.7\n2\nframe 2\nH 0 0 0\nH 0 0 1.2\n")
            output = render_xyz_trajectory_gif(xyz, folder / "route.gif")
            self.assertTrue(output.is_file())

    def test_orca_reaction_path_inputs_retain_endpoints_and_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            reactant = folder / "reactant.xyz"
            product = folder / "product.xyz"
            reactant.write_text("2\nreactant\nH 0 0 0\nH 0 0 0.7\n")
            product.write_text("2\nproduct\nH 0 0 0\nH 0 0 3.5\n")
            neb = create_orca_neb_ts_input(reactant, product, charge=0, multiplicity=1)
            preoptimized_neb = create_orca_neb_ts_input(
                reactant, product, charge=0, multiplicity=1, label="preoptimized-neb",
                preopt_ends=True,
            )
            scan = create_orca_relaxed_scan_input(
                reactant, bond_atom_indices=(0, 1), start_distance_angstrom=0.7,
                end_distance_angstrom=3.5, charge=0, multiplicity=1,
            )
            hessian = folder / "ts.hess"
            hessian.write_text("synthetic")
            irc = create_orca_irc_input(reactant, charge=0, multiplicity=1, hessian_file=hessian)
            self.assertIn('Product "product.xyz"', neb.read_text())
            self.assertNotIn("PreOptEnds true", neb.read_text())
            self.assertIn("PreOptEnds true", preoptimized_neb.read_text())
            self.assertIn("B 0 1 = 0.70000000, 3.50000000, 20", scan.read_text())
            self.assertIn('Hess_Filename "ts.hess"', irc.read_text())

    def test_dissociation_endpoint_preserves_atom_order_and_separates_bond(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            reactant = folder / "reactant.xyz"
            from src.molecule_tools import smiles_to_xyz
            smiles_to_xyz("OO", reactant)
            route = RouteSpec(
                route_id="oo-break", source="test", parent_smiles="OO",
                reactant_smiles=("OO",), product_smiles=("[OH]", "[OH]"),
                broken_bonds=((0, 1),), atom_mapping=((0, 0), (1, 1)),
                product_multiplicities=(2, 2), protocol="relaxed_dissociation_scan",
            )
            result = build_dissociation_endpoint(route, reactant, folder / "product.xyz")
            reactant_elements, _, _ = read_xyz(reactant)
            product_elements, coordinates, _ = read_xyz(folder / "product.xyz")
            distance = np.linalg.norm(coordinates[0] - coordinates[1])
            self.assertEqual(reactant_elements, product_elements)
            self.assertAlmostEqual(distance, result["final_distance_angstrom"], places=5)

    def test_computational_light_route_can_prepare_labelled_animation_preview(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            report = folder / "computational-light.json"
            report.write_text(json.dumps({
                "schema_version": 1,
                "kind": "computational_light_model",
                "status": "completed",
                "smiles": "OO",
                "candidate_routes": [{
                    "route_id": "homolysis-1",
                    "kind": "homolysis",
                    "bond_atom_indices": [0, 1],
                    "fragment_radical_smiles": ["[OH]", "[OH]"],
                    "accessibility": 0.9,
                    "status": "computed",
                }],
            }))
            _, routes = discover_routes(report)
            selected = select_explanatory_route(routes)
            self.assertEqual(selected.protocol, "relaxed_dissociation_scan")
            result = run_decomposition_explanation(report, prepare_only=True)
            self.assertEqual(result["status"], "prepared")
            self.assertTrue(Path(result["visuals"]["animation"]).is_file())
            manifest = json.loads(Path(result["visuals"]["manifest"]).read_text())
            self.assertEqual(manifest["evidence_level"], "candidate_preview")

    def test_verification_plan_selects_upstream_radical_formation_before_collision_route(self):
        from storca.verification_plan import build_verification_dependency_plan
        rmg = {
            "kinetics_validation": {"violators": [{
                "reaction_equation": "oxygen(2)+radical(12)<=>adduct(42)",
            }]},
            "candidate_routes": [{
                "route_id": "rmg-1",
                "reaction_equation": "oxygen(2)+stability(1)<=>peroxy(13)+radical(12)",
                "resolved_endpoints": {
                    "reactants": [{"label": "oxygen(2)"}, {"label": "stability(1)"}],
                    "products": [{"label": "peroxy(13)"}, {"label": "radical(12)"}],
                },
            }],
        }
        plan = build_verification_dependency_plan(
            rmg, {"composition": {"stability": 0.01, "oxygen": 0.21, "nitrogen": 0.78}},
        )
        self.assertEqual(plan["status"], "upstream_initiation_verification_required")
        self.assertEqual(plan["selected_route_index"], 0)
        self.assertEqual(plan["blocked_downstream_routes"][0]["generated_reactants"], ["radical(12)"])


if __name__ == "__main__":
    unittest.main()
