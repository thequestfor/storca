from __future__ import annotations

import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from storca.reaction_path import (
    _neb_angle_warning_image_indices,
    _neb_intermediate_image_indices,
    _orca_neb_no_barrier_outcome,
    _run_generic_orientation,
    _run_or_resume,
    _stationary_fragment_records,
    classify_path_energy_profile,
    run_generic_reaction_path,
)
from storca.route_geometry import (
    build_endpoint_complex_seeds,
    read_xyz,
    resolve_route_atom_mapping,
    validate_declared_connectivity,
    write_xyz,
)
from storca.route_verify import RouteSpec, route_spec_from_rmg_route


class GenericReactionPathTests(unittest.TestCase):
    def test_orca_intermediate_image_diagnostic_is_parsed(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "neb-ts.out"
            output.write_text(
                "Possible intermediate minimum found at image(s): 1 4\n"
                "Possible intermediate minimum found at image(s): 4 7\n"
            )
            self.assertEqual(_neb_intermediate_image_indices(output), [1, 4, 7])

    def test_orca_angle_warning_ignores_corrupt_first_integer(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "neb-ts.out"
            output.write_text(
                "Warning: maximum angle between images (-1659107148 - 5 - 6): 0.00\n"
                "Warning: maximum angle between images (745032783 - 5 - 6): 0.00\n"
            )
            self.assertEqual(_neb_angle_warning_image_indices(output), [5, 6])

    @patch("storca.reaction_path.run_orca")
    def test_neb_intermediate_warning_stops_and_is_retained(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "neb-ts.inp"
            inp.write_text(
                "! B3LYP def2-SVP NEB-TS\n"
                '%neb\n  Product "product.xyz"\nend\n'
                "* xyzfile 0 1 reactant.xyz\n"
            )
            (folder / "reactant.xyz").write_text("1\nreactant\nH 0 0 0\n")
            (folder / "product.xyz").write_text("1\nproduct\nH 0 0 1\n")
            output = inp.with_suffix(".out")

            def stopped(input_path, **kwargs):
                output.write_text(
                    "Possible intermediate minimum found at image(s): 2\n"
                )
                self.assertFalse(kwargs["live_output"])
                return {"out": output, "stopped_early": True}

            run_orca.side_effect = stopped
            result = _run_or_resume(inp, timeout_seconds=10.0)
            self.assertEqual(result["status"], "completed_intermediate_detected")
            self.assertEqual(result["intermediate_image_indices"], [2])
            self.assertEqual(result["early_stop"]["status"], "stopped_after_orca_intermediate_warning")

    @staticmethod
    def _association_route() -> RouteSpec:
        return RouteSpec(
            route_id="association",
            source="test",
            parent_smiles="[CH3]",
            reactant_smiles=("[CH3]", "[OH]"),
            product_smiles=("CO",),
            multiplicity=1,
            reactant_multiplicities=(2, 2),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )

    def test_unique_multifragment_association_maps_and_builds_orientations(self):
        route = self._association_route()
        mapping = resolve_route_atom_mapping(route)
        self.assertTrue(mapping["valid"])
        self.assertEqual(mapping["surface"]["common_total_multiplicities"], [1])
        with tempfile.TemporaryDirectory() as temp:
            seeds = build_endpoint_complex_seeds(route, Path(temp), mapping_result=mapping)
            self.assertEqual(seeds["orientation_count"], 3)
            self.assertEqual(seeds["bond_changes"]["formed_bonds"], [[0, 4]])
            first = seeds["orientations"][0]
            left, _, _ = read_xyz(Path(first["reactant_xyz"]))
            right, _, _ = read_xyz(Path(first["product_xyz"]))
            self.assertEqual(left, right)

    def test_automatic_mapping_canonicalizes_true_endpoint_symmetry(self):
        route = RouteSpec(
            route_id="symmetric-identity",
            source="test",
            parent_smiles="CC",
            reactant_smiles=("CC",),
            product_smiles=("CC",),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["source"], "automatic_minimum_graph_edit_symmetry_canonicalized")
        self.assertEqual(result["symmetry_resolution"]["equivalent_mapping_count"], 2)

    def test_automatic_mapping_canonicalizes_o2_symmetry_for_nno_initiation(self):
        route = RouteSpec(
            route_id="nno-o2-initiation",
            source="test",
            parent_smiles="NNO",
            reactant_smiles=("[O][O]", "NNO"),
            product_smiles=("[O]O", "NN[O]"),
            multiplicity=3,
            reactant_multiplicities=(3, 1),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["surface"]["selected_multiplicity"], 3)
        self.assertEqual(result["source"], "automatic_minimum_graph_edit_symmetry_canonicalized")
        self.assertEqual(result["symmetry_resolution"]["equivalent_mapping_count"], 2)
        self.assertEqual(len(result["all_atom_mapping"]), 9)

    def test_automatic_mapping_canonicalizes_identical_reactant_exchange(self):
        """Two identical reactant molecules may exchange without ambiguity."""
        route = RouteSpec(
            route_id="hno-disproportionation",
            source="test",
            parent_smiles="N=O",
            reactant_smiles=("N=O", "N=O"),
            product_smiles=("[N]=O", "N[O]"),
            multiplicity=1,
            reactant_multiplicities=(1, 1),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["source"], "automatic_minimum_graph_edit_symmetry_canonicalized")
        self.assertEqual(result["symmetry_resolution"]["equivalent_mapping_count"], 2)

    def test_automatic_mapping_keeps_distinct_isomerization_maps_ambiguous(self):
        route = RouteSpec(
            route_id="non-equivalent-maps",
            source="test",
            parent_smiles="CCCO",
            reactant_smiles=("CCCO",),
            product_smiles=("CC(C)O",),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertFalse(result["valid"])
        self.assertEqual(result["status"], "ambiguous_atom_mapping")

    def test_supplied_mapping_extends_across_hydrogen_transfer(self):
        route = RouteSpec(
            route_id="abstraction",
            source="test",
            parent_smiles="CCO",
            reactant_smiles=("CCO", "[O][O]"),
            product_smiles=("CC[O]", "[O]O"),
            multiplicity=3,
            reactant_multiplicities=(1, 3),
            product_multiplicities=(2, 2),
            atom_mapping=((0, 0), (1, 1), (2, 2), (3, 3), (4, 4)),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"])
        self.assertEqual(len(result["hydrogen_transfers"]), 1)
        self.assertEqual(result["surface"]["common_total_multiplicities"], [3])

    def test_standalone_hydrogen_fragment_uses_all_hydrogen_mapping_layer(self):
        route = RouteSpec(
            route_id="nno-h-loss",
            source="test",
            parent_smiles="NNO",
            reactant_smiles=("NNO",),
            product_smiles=("[H]", "[NH]NO"),
            reactant_labels=("stability(1)",),
            product_labels=("H(2)", "radical(3)"),
            multiplicity=1,
            reactant_multiplicities=(1,),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["reactant_skeleton_atom_count"], 3)
        self.assertEqual(result["product_skeleton_atom_count"], 3)
        self.assertEqual(result["explicit_atom_count"], 7)
        self.assertEqual(len(result["all_atom_mapping"]), 7)
        self.assertEqual(len(result["hydrogen_transfers"]), 1)
        transfer = result["hydrogen_transfers"][0]
        self.assertEqual(transfer["reactant_donor_index"], 0)
        self.assertIsNone(transfer["product_acceptor_index"])
        self.assertEqual(transfer["product_fragment_index"], 0)
        self.assertEqual(result["surface"]["charge"], 0)
        self.assertEqual(result["surface"]["selected_multiplicity"], 1)

        with tempfile.TemporaryDirectory() as temp:
            seeds = build_endpoint_complex_seeds(route, Path(temp), mapping_result=result)
            self.assertEqual(seeds["orientation_count"], 3)
            first = seeds["orientations"][0]
            reactant_elements, _, _ = read_xyz(Path(first["reactant_xyz"]))
            product_elements, _, _ = read_xyz(Path(first["product_xyz"]))
            self.assertEqual(reactant_elements, product_elements)
            self.assertEqual(len(product_elements), 7)

            prepared = run_generic_reaction_path(
                route, Path(temp) / "prepared", prepare_only=True,
            )
            h_record = prepared["route_verification"]["stationary_points"]["products"][0]
            self.assertEqual(read_xyz(Path(h_record["input_xyz"]))[0], ["H"])

    def test_standalone_charged_hydrogen_preserves_fragment_state(self):
        for smiles, charge in (("[H+]", 1), ("[H-]", -1)):
            with self.subTest(smiles=smiles):
                route = RouteSpec(
                    route_id=f"identity-{charge}",
                    source="test",
                    parent_smiles=smiles,
                    reactant_smiles=(smiles,),
                    product_smiles=(smiles,),
                    charge=charge,
                    multiplicity=1,
                    reactant_charges=(charge,),
                    product_charges=(charge,),
                    reactant_multiplicities=(1,),
                    product_multiplicities=(1,),
                    protocol="generic_endpoint_path",
                )
                result = resolve_route_atom_mapping(route)
                self.assertTrue(result["valid"], result)
                self.assertEqual(result["reactant_skeleton_atom_count"], 0)
                self.assertEqual(result["product_skeleton_atom_count"], 0)
                self.assertEqual(result["explicit_atom_count"], 1)
                self.assertEqual(result["all_atom_mapping"], [[0, 0]])
                self.assertEqual(result["surface"]["charge"], charge)
                self.assertEqual(result["surface"]["selected_multiplicity"], 1)

    def test_standalone_hydrogen_isotope_imbalance_fails_closed(self):
        route = RouteSpec(
            route_id="isotope-imbalance",
            source="test",
            parent_smiles="[2H]",
            reactant_smiles=("[2H]",),
            product_smiles=("[H]",),
            multiplicity=2,
            reactant_multiplicities=(2,),
            product_multiplicities=(2,),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertFalse(result["valid"])
        self.assertEqual(result["status"], "unbalanced_elements")

    def test_connectivity_gate_distinguishes_bound_and_separated_hydrogen(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            bound = write_xyz(
                folder / "bound.xyz", ["H", "H"],
                np.asarray([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]),
                "mock optimized H2",
            )
            separated = write_xyz(
                folder / "separated.xyz", ["H", "H"],
                np.asarray([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
                "mock separated H atoms",
            )
            self.assertTrue(validate_declared_connectivity(bound, ("[H][H]",))["valid"])
            self.assertTrue(validate_declared_connectivity(separated, ("[H]", "[H]"))["valid"])
            collapsed = validate_declared_connectivity(bound, ("[H]", "[H]"))
            dissociated = validate_declared_connectivity(separated, ("[H][H]",))
            self.assertFalse(collapsed["valid"])
            self.assertEqual(collapsed["status"], "connectivity_changed")
            self.assertEqual(collapsed["unexpected_covalent_bonds"][0]["atom_indices"], [0, 1])
            self.assertFalse(dissociated["valid"])
            self.assertEqual(dissociated["missing_declared_bonds"][0]["atom_indices"], [0, 1])

    @patch("storca.reaction_path._optimize_endpoint")
    def test_collapsed_product_complex_is_rejected_before_neb(self, optimize):
        route = self._association_route()
        mapping = resolve_route_atom_mapping(route)
        with tempfile.TemporaryDirectory() as temp:
            seeds = build_endpoint_complex_seeds(
                route, Path(temp) / "seeds", orientations=1, mapping_result=mapping,
            )
            reactant_seed = Path(seeds["orientations"][0]["reactant_xyz"])
            optimize.return_value = {
                # Both optimizations are mocked as the separated reactant
                # complex: the declared bonded product has collapsed to the
                # wrong endpoint channel.
                "optimized_xyz": reactant_seed,
                "optimization": {"out": Path(temp) / "opt.out"},
                "frequency": {"out": Path(temp) / "freq.out"},
                "frequency_check": {"IsMinimum": True},
            }
            result = _run_generic_orientation(
                route,
                seeds["orientations"][0],
                mapping_result=mapping,
                ncores=1,
                method_keywords=None,
                timeout_seconds=1.0,
                neb_images=4,
                barrierless_classifier=None,
            )
        self.assertEqual(result["status"], "endpoint_connectivity_validation_failed")
        self.assertTrue(result["endpoint_validation"]["reactant"]["connectivity_check"]["valid"])
        self.assertFalse(result["endpoint_validation"]["product"]["connectivity_check"]["valid"])

    @patch("storca.reaction_path._optimize_endpoint")
    def test_separated_stationary_species_must_retain_connectivity(self, optimize):
        route = RouteSpec(
            route_id="mock-h2",
            source="test",
            parent_smiles="[H][H]",
            reactant_smiles=("[H][H]",),
            product_smiles=("[H][H]",),
            multiplicity=1,
            reactant_multiplicities=(1,),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        mapping = resolve_route_atom_mapping(route)
        with tempfile.TemporaryDirectory() as temp:
            dissociated = write_xyz(
                Path(temp) / "dissociated.xyz", ["H", "H"],
                np.asarray([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
                "mock wrong minimum",
            )
            optimize.return_value = {
                "optimized_xyz": dissociated,
                "optimization": {"out": Path(temp) / "opt.out"},
                "frequency": {"out": Path(temp) / "freq.out"},
                "frequency_check": {"IsMinimum": True},
            }
            records = _stationary_fragment_records(
                route, Path(temp) / "stationary", mapping,
                execute=True, ncores=1, method_keywords=None, timeout_seconds=1.0,
            )
        self.assertFalse(records["valid"])
        self.assertEqual(records["status"], "stationary_point_validation_failed")
        self.assertTrue(all(
            record["status"] == "endpoint_connectivity_validation_failed"
            for side in ("reactants", "products") for record in records[side]
        ))

    def test_noncommon_declared_spin_surface_is_rejected(self):
        route = RouteSpec(
            route_id="wrong-spin",
            source="test",
            parent_smiles="[CH3]",
            reactant_smiles=("[CH3]", "[OH]"),
            product_smiles=("CO",),
            multiplicity=3,
            reactant_multiplicities=(2, 2),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertFalse(result["valid"])
        self.assertEqual(result["status"], "declared_multiplicity_not_common")

    def test_profile_classifier_separates_barriered_and_barrierless(self):
        barriered = classify_path_energy_profile(
            [-100.0, -99.99, -99.96, -99.99, -100.01],
            reactant_fragment_count=2,
            product_fragment_count=1,
        )
        capture = classify_path_energy_profile(
            [-100.0, -100.01, -100.02, -100.03, -100.04],
            reactant_fragment_count=2,
            product_fragment_count=1,
        )
        self.assertEqual(barriered["status"], "barriered_ts_candidate")
        self.assertEqual(capture["status"], "barrierless_capture_candidate")

    def test_equal_fragment_no_saddle_profile_retains_reaction_direction(self):
        uphill = classify_path_energy_profile(
            [-336.942205, -336.950666, -336.944031, -336.929571, -336.903440],
            reactant_fragment_count=2,
            product_fragment_count=2,
        )
        self.assertEqual(uphill["status"], "surface_unresolved")
        self.assertEqual(uphill["downhill_direction"], "reverse")
        self.assertFalse(uphill["forward_target_loss_barrierless"])
        self.assertAlmostEqual(uphill["forward_energetic_threshold_hartree"], 0.038765)
        self.assertEqual(uphill["reason"], "forward_path_uphill_without_interior_ts_candidate")

    def test_orca_converged_no_barrier_output_is_semantic_completion(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "neb-ts.inp"
            out = folder / "neb-ts.out"
            inp.write_text("! B3LYP def2-SVP NEB-TS\n")

            def failed_orca(input_path, **kwargs):
                out.write_text(
                    "*** THE NEB OPTIMIZATION HAS CONVERGED ***\n"
                    "The elastic band has converged successfully to a MEP. "
                    "However, climbing image was never activated!\n"
                    "No barrier was found. Skipping NEB-TS run here.\n"
                    "ORCA finished with error return - aborting the run\n"
                )
                raise RuntimeError("ORCA failed")

            with patch("storca.reaction_path.run_orca", side_effect=failed_orca):
                result = _run_or_resume(inp, timeout_seconds=10.0)

            self.assertEqual(result["status"], "completed_no_interior_barrier")
            self.assertTrue(result["orca_neb_outcome"]["valid"])
            self.assertFalse(result["orca_neb_outcome"]["climbing_image_activated"])
            self.assertTrue(inp.with_suffix(".contract.json").is_file())

            resumed = _run_or_resume(inp, timeout_seconds=10.0)
            self.assertTrue(resumed["resumed"])
            self.assertEqual(resumed["status"], "completed_no_interior_barrier")

    def test_orca_no_barrier_parser_requires_convergence_and_explicit_message(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "neb-ts.out"
            output.write_text("No barrier was found. Skipping NEB-TS run here.\n")
            parsed = _orca_neb_no_barrier_outcome(output)
        self.assertFalse(parsed["valid"])
        self.assertFalse(parsed["neb_converged"])

    def test_prepare_only_returns_canonical_fail_closed_contract(self):
        route = RouteSpec(
            route_id="prepare",
            source="test",
            parent_smiles="[CH3]",
            reactant_smiles=("[CH3]", "[OH]"),
            product_smiles=("CO",),
            multiplicity=1,
            reactant_multiplicities=(2, 2),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(route, Path(temp), prepare_only=True)
            verification = result["route_verification"]
            self.assertEqual(verification["status"], "prepared")
            self.assertEqual(verification["rate_model"]["status"], "rate_model_unavailable")
            self.assertIsNone(verification["trajectory"])
            self.assertEqual(len(result["orientations"]), 3)
            self.assertEqual(
                result["endpoint_complexes"]["fragment_geometry_source"],
                "supplied_atom_ordered_fragment_xyz",
            )
            self.assertEqual(verification["endpoint_mode"], "separated_fragment_channels")

    @patch("storca.reaction_path.build_endpoint_complex_seeds")
    @patch("storca.reaction_path.condition_specific_forward_loss_bound")
    @patch("storca.reaction_path.assemble_route_thermochemistry")
    @patch("storca.reaction_path._stationary_fragment_records")
    def test_thermochemical_loss_bound_stops_before_neb(
        self, stationary, assemble_thermo, loss_bound, build_seeds,
    ):
        stationary.return_value = {
            "status": "validated", "valid": True, "reactants": [], "products": [],
        }
        assemble_thermo.return_value = {"status": "validated", "valid": True}
        loss_bound.return_value = {
            "status": "forward_loss_below_retention_threshold_upper_bound",
            "applicable": True,
            "target_loss_fraction_upper_bound": 0.001,
        }
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(
                self._association_route(), Path(temp), condition_contract={
                    "temperature_K": 298.15, "pressure_bar": 1.0,
                    "target_duration_seconds": 1000.0, "retention_fraction": 0.95,
                },
            )
        verification = result["route_verification"]
        self.assertEqual(verification["status"], "thermochemistry_bounded_below_loss_threshold")
        self.assertEqual(verification["convergence"]["status"], "not_required_after_conservative_loss_bound")
        build_seeds.assert_not_called()

    def test_public_rmg_route_normalizer_retains_fragment_states(self):
        route = route_spec_from_rmg_route("[CH3]", {
            "route_id": "rmg-1",
            "reaction_equation": "methyl(1) + hydroxyl(2) <=> methanol(3)",
            "multiplicity": 1,
            "resolved_endpoints": {
                "reactants": [
                    {"label": "methyl(1)", "smiles": "[CH3]", "charge": 0, "multiplicity": 2},
                    {"label": "hydroxyl(2)", "smiles": "[OH]", "charge": 0, "multiplicity": 2},
                ],
                "products": [{"label": "methanol(3)", "smiles": "CO", "charge": 0, "multiplicity": 1}],
            },
        })
        self.assertEqual(route.protocol, "generic_endpoint_path")
        self.assertEqual(route.reactant_charges, (0, 0))
        self.assertEqual(route.reactant_multiplicities, (2, 2))
        self.assertEqual(route.product_multiplicities, (1,))
        self.assertEqual(route.reactant_labels, ("methyl(1)", "hydroxyl(2)"))
        self.assertEqual(route.reaction_equation, "methyl(1) + hydroxyl(2) <=> methanol(3)")

    def test_stoichiometric_duplicates_reuse_original_stationary_label(self):
        route = RouteSpec(
            route_id="homolysis",
            source="test",
            parent_smiles="OO",
            reactant_smiles=("OO",),
            product_smiles=("[OH]", "[OH]"),
            reactant_labels=("peroxide",),
            product_labels=("hydroxyl", "hydroxyl"),
            multiplicity=1,
            reactant_multiplicities=(1,),
            product_multiplicities=(2, 2),
            atom_mapping=((0, 0), (1, 1)),
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(route, Path(temp), prepare_only=True)
            products = result["route_verification"]["stationary_points"]["products"]
            self.assertEqual([record["label"] for record in products], ["hydroxyl", "hydroxyl"])
            self.assertEqual(products[0]["input_xyz"], products[1]["input_xyz"])
            self.assertTrue(products[1]["reused_stationary_point"])

    @patch("storca.reaction_path._stationary_fragment_xyz", return_value=None)
    @patch("storca.reaction_path._stationary_fragment_records")
    @patch("storca.reaction_path._run_generic_orientation")
    def test_reproducible_barrierless_path_stops_after_two_orientations(
        self, run_orientation, stationary, stationary_xyz,
    ):
        def result(seed):
            orientation = seed["orientation"]
            return {
                "orientation": orientation,
                "status": "classified_without_transition_state",
                "path_classification": "barrierless_capture_candidate",
                "endpoint_validation": {"valid": True},
                "trajectory": f"orientation-{orientation}.xyz",
                "path_evidence": {
                    "energies_hartree": [-100.0, -100.01, -100.02, -100.03, -100.04],
                },
            }

        run_orientation.side_effect = lambda route, seed, **kwargs: result(seed)
        stationary.return_value = {"status": "validated", "valid": True, "reactants": [], "products": []}
        with tempfile.TemporaryDirectory() as temp:
            verification = run_generic_reaction_path(self._association_route(), Path(temp))["route_verification"]
        self.assertEqual(run_orientation.call_count, 2)
        self.assertEqual(verification["status"], "barrierless_rate_model_required")
        self.assertEqual(verification["barrierless_reproducibility"]["status"], "reproducible")

    @patch("storca.reaction_path._stationary_fragment_xyz", return_value=None)
    @patch("storca.reaction_path._stationary_fragment_records")
    @patch("storca.reaction_path._run_generic_orientation")
    def test_verified_ts_path_stops_after_first_orientation(
        self, run_orientation, stationary, stationary_xyz,
    ):
        run_orientation.return_value = {
            "orientation": 1,
            "status": "verified_transition_state_path",
            "path_classification": "verified_barriered_path",
            "endpoint_validation": {"valid": True},
            "transition_state": {"frequency_output": "ts-frequency.out"},
            "irc": {"endpoint_validation": {"valid": True}},
            "trajectory": "validated-trajectory.xyz",
            "path_evidence": {"kind": "orca_ts_frequency_irc"},
        }
        stationary.return_value = {"status": "validated", "valid": True, "reactants": [], "products": []}
        with tempfile.TemporaryDirectory() as temp:
            verification = run_generic_reaction_path(self._association_route(), Path(temp))["route_verification"]
        self.assertEqual(run_orientation.call_count, 1)
        self.assertEqual(verification["status"], "verified_transition_state_path")
        self.assertEqual(verification["trajectory"], "validated-trajectory.xyz")

    @patch("storca.reaction_path._run_or_resume", side_effect=RuntimeError("mock NEB stop"))
    @patch("storca.reaction_path._optimize_endpoint")
    def test_separated_nno_o2_channel_bypasses_encounter_complex_optimization(
        self, optimize_endpoint, run_or_resume,
    ):
        route = RouteSpec(
            route_id="nno-o2-initiation",
            source="test",
            parent_smiles="NNO",
            reactant_smiles=("[O][O]", "NNO"),
            product_smiles=("[O]O", "NN[O]"),
            multiplicity=3,
            reactant_multiplicities=(3, 1),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        mapping = resolve_route_atom_mapping(route)
        with tempfile.TemporaryDirectory() as temp:
            seeds = build_endpoint_complex_seeds(
                route, Path(temp) / "seeds", orientations=1, mapping_result=mapping,
            )
            result = _run_generic_orientation(
                route, seeds["orientations"][0], mapping_result=mapping,
                ncores=1, method_keywords=None, timeout_seconds=1.0,
                neb_images=4, barrierless_classifier=None,
                separated_fragments=True,
                stationary_points={"status": "validated", "valid": True},
            )
        optimize_endpoint.assert_not_called()
        run_or_resume.assert_called_once()
        self.assertEqual(
            result["endpoint_validation"]["mode"],
            "assembled_from_validated_separated_fragment_minima",
        )
        self.assertFalse(result["endpoint_validation"]["full_complex_minimum_required"])
        self.assertEqual(result["status"], "surface_unresolved")

    @patch("storca.reaction_path.build_endpoint_complex_seeds")
    @patch("storca.reaction_path._stationary_fragment_records")
    def test_invalid_separated_species_blocks_before_encounter_assembly(self, stationary, build_seeds):
        stationary.return_value = {
            "status": "stationary_point_validation_failed",
            "valid": False,
            "reactants": [{"status": "failed", "failure_reason": "mock failure"}],
            "products": [],
        }
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(self._association_route(), Path(temp))
        build_seeds.assert_not_called()
        verification = result["route_verification"]
        self.assertEqual(verification["status"], "stationary_point_validation_failed")
        self.assertEqual(verification["convergence"]["status"], "blocked_before_neb")


if __name__ == "__main__":
    unittest.main()
