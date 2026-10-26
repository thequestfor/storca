import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from propagation_fixtures import isothermal_model
from arkane_rotor_fixtures import rotor_rate_evidence_fixture

from storca.rmg_execution import (
    merge_collision_rate_validation,
    validate_inspected_collision_rates,
)
from storca.stability import run_stability_screen
from storca.verification_engine import _planning_replacements, run_verification_engine


def _route():
    return {
        "route_id": "route-stable",
        "reaction": "stability(1)+nitrogen(2)<=>P(3)",
        "reaction_equation": "stability(1)+nitrogen(2)<=>P(3)",
        "activation_energy_kcal_mol": 2.0,
        "route_context": "co_reactant_dependent",
        "initiation_status": "directly_present",
        "co_reactants": [{"label": "nitrogen", "origin": "initial_scenario_species"}],
        "resolved_endpoints": {
            "reactants": [{"label": "stability(1)"}, {"label": "nitrogen(2)"}],
            "products": [{"label": "P(3)"}],
        },
    }


def _propagation(*, loss=0.20, t95=50.0, equation=None):
    equation = equation or _route()["reaction_equation"]
    return {
        "status": "completed",
        "coverage": {"complete": True, "requested_seconds": 86400.0, "simulated_seconds": 86400.0},
        "estimated_time_to_retention_seconds": t95,
        "reactor_model": isothermal_model(298.0),
        "target_loss_fraction": loss,
        "target_profile": [
            {"time_seconds": 0.0, "target_fraction_remaining": 1.0},
            {"time_seconds": 86400.0, "target_fraction_remaining": 1.0 - loss},
        ],
        "reaction_flux_attribution": {
            "status": "completed",
            "total_integrated_target_destruction_kmol": 1.0,
            "numerical_closure_relative_error": 0.0,
            "reactions": [{
                "reaction_index": 0,
                "reaction_equation": equation,
                "integrated_target_destruction_kmol": 1.0,
                "integrated_forward_extent_kmol": 1.0,
                "integrated_reverse_extent_kmol": 0.0,
            }],
        },
    }


def _rmg_evidence(*, loss=0.20, t95=50.0, routes=None, kinetics=None):
    routes = [_route()] if routes is None else routes
    propagation = _propagation(loss=loss, t95=t95)
    if not routes:
        propagation["reaction_flux_attribution"].update(
            total_integrated_target_destruction_kmol=0.0,
            reactions=[],
        )
    return {
        "status": "completed",
        "execution_assessment": {"status": "completed"},
        "time_coverage": {"complete": True},
        "candidate_routes": routes,
        "network_routes": [],
        "independent_cantera_propagation": propagation,
        "kinetic_relevance": {
            "status": "kinetically_relevant_candidate" if loss >= 0.05 else "reachable_but_below_loss_threshold",
            "estimated_time_to_retention_seconds": t95,
        },
        "kinetics_validation": kinetics or {"status": "passed", "violation_count": 0},
        "generated_kinetics_libraries": [],
        "mechanism_inspection": {"status": "completed", "reactions": []},
        "artifacts": {"chemkin": "chemkin/chem_annotated.inp"},
    }


def _orca(*, minimum=True):
    return {"status": "completed", "local_minimum": minimum, "artifacts": {}}


def _inspection(*, missing_reverse_limit=False):
    return {
        "status": "completed",
        "reactions": [{
            "reaction_index": 0,
            "equation": "A(1)+B(2)<=>C(3)+D(4)",
            "molecularity": 2,
            "reverse_molecularity": 2,
            "evaluated_forward_rate_coefficient_si": 2.0,
            "evaluated_reverse_rate_coefficient_si": 3.0,
            "forward_collision_limit_si": 4.0,
            "reverse_collision_limit_si": None if missing_reverse_limit else 5.0,
            "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
            "collision_limit_source": "RMG retained transport calculation",
        }],
    }


class CollisionValidationIntegrationTests(unittest.TestCase):
    def test_absent_optional_log_can_pass_only_from_complete_inspection(self):
        result = merge_collision_rate_validation(
            {"status": "not_reported", "violation_count": 0, "path": "missing.log"},
            _inspection(),
            temperature_K=298.15,
            pressure_bar=1.0,
        )
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["condition_mechanism_validation"]["applicable_direction_count"], 2)
        self.assertEqual(result["condition_mechanism_validation"]["validated_direction_count"], 2)

    def test_one_missing_directional_ceiling_remains_incomplete(self):
        result = validate_inspected_collision_rates(
            _inspection(missing_reverse_limit=True), temperature_K=298.15, pressure_bar=1.0,
        )
        self.assertEqual(result["status"], "incomplete")
        self.assertEqual(result["validated_direction_count"], 1)
        self.assertEqual(result["missing_or_invalid_directions"][0]["direction"], "reverse")

    def test_explicit_log_violation_is_never_erased_by_clean_inspection(self):
        reported = {
            "status": "kinetics_unreliable",
            "violation_count": 1,
            "violators": [{
                "reaction_equation": "X+Y<=>Z", "direction": "forward", "violation_factor": 7.0,
            }],
        }
        result = merge_collision_rate_validation(
            reported, _inspection(), temperature_K=298.15, pressure_bar=1.0,
        )
        self.assertEqual(result["status"], "kinetics_unreliable")
        self.assertEqual(result["violation_count"], 1)
        self.assertEqual(result["maximum_violation_factor"], 7.0)

    def test_irreversible_reaction_does_not_require_nonexistent_reverse_bound(self):
        inspection = _inspection(missing_reverse_limit=True)
        inspection["reactions"][0]["reversible"] = False
        result = validate_inspected_collision_rates(
            inspection, temperature_K=298.15, pressure_bar=1.0,
        )
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["applicable_direction_count"], 1)
        self.assertEqual(result["validated_direction_count"], 1)

    def test_missing_molecularity_cannot_vacuously_pass(self):
        inspection = _inspection()
        inspection["reactions"][0].pop("molecularity")
        result = validate_inspected_collision_rates(
            inspection, temperature_K=298.15, pressure_bar=1.0,
        )
        self.assertEqual(result["status"], "incomplete")
        self.assertIn(
            "reaction_molecularity_missing_or_invalid",
            result["missing_or_invalid_directions"][0]["reasons"],
        )


class StabilityDefaultVerificationWiringTests(unittest.TestCase):
    def test_conditional_model_prediction_is_not_an_orca_failure_or_lifetime(self):
        prediction = {"target_profile": [{"time_seconds": 86400.0, "target_fraction_remaining": 0.88}],
                      "conclusive_observable_allowed": False}

        def conditional_engine(*args, **kwargs):
            return {"status": "conditional_model_prediction", "condition_contract": args[1],
                    "reason": "Fixed-model profile only; physical applicability and uncertainty unresolved.",
                    "estimated_time_to_retention_seconds": None,
                    "conditional_model_prediction": prediction}

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=_rmg_evidence()), \
                patch("storca.verification_engine.run_verification_engine", side_effect=conditional_engine), \
                patch("storca.verification_adapters.make_orca_route_verifier", return_value=lambda *_: {}), \
                patch("storca.verification_adapters.make_arkane_rate_builder", return_value=lambda *_: {}):
            result = run_stability_screen("NNO", Path(temp),
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0})
        self.assertEqual(result["assessment"]["status"], "conditional_kinetic_model_prediction")
        self.assertEqual(result["conditional_model_prediction"], prediction)
        self.assertIsNone(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"])
        self.assertIn("conditional_model_only", result["risk_flags"])

    def test_rate_uncertainty_is_validated_retained_and_forwarded_to_engine(self):
        initial = _rmg_evidence()
        supplied = [{
            "route_id": "route-stable",
            "lower_rate_multiplier": 0.2,
            "upper_rate_multiplier": 5.0,
            # This test exercises forwarding explicit user input, not
            # condition-specific calibration. The screen creates its contract
            # later, so a global fixture must not claim calibration authority.
            "source": "explicit test user allocation bound",
            "source_kind": "user_supplied",
            "applicability": "test fixture chemistry",
        }]
        captured = {}

        def fake_engine(*args, **kwargs):
            captured["records"] = kwargs["rate_uncertainties"]
            return {
                "status": "verification_incomplete",
                "reason": "fixture stopped after wiring",
                "condition_contract": args[1],
                "estimated_time_to_retention_seconds": None,
            }

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial), \
                patch("storca.verification_engine.run_verification_engine", side_effect=fake_engine), \
                patch("storca.verification_adapters.make_orca_route_verifier", return_value=lambda *_: {}), \
                patch("storca.verification_adapters.make_arkane_rate_builder", return_value=lambda *_: {}):
            result = run_stability_screen(
                "NNO", Path(temp), rate_uncertainty_evidence=supplied,
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )
        self.assertEqual(len(captured["records"]), 1)
        self.assertEqual(captured["records"][0].route_id, "route-stable")
        self.assertEqual(result["rate_uncertainty_evidence"]["record_count"], 1)

    def test_malformed_rate_uncertainty_stops_before_orca_or_rmg(self):
        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence") as collect_orca, \
                patch("storca.stability.collect_rmg_evidence") as collect_rmg:
            with self.assertRaisesRegex(ValueError, "applicability"):
                run_stability_screen(
                    "NNO", Path(temp),
                    rate_uncertainty_evidence=[{
                        "route_id": "route-stable",
                        "lower_rate_multiplier": 0.2,
                        "upper_rate_multiplier": 5.0,
                        "source": "test",
                        "source_kind": "user_supplied",
                    }],
                    method_profile={"orca_keywords": "TEST"},
                )
        collect_orca.assert_not_called()
        collect_rmg.assert_not_called()

    def test_unbound_method_calibration_stops_before_orca_or_rmg(self):
        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence") as collect_orca, \
                patch("storca.stability.collect_rmg_evidence") as collect_rmg:
            with self.assertRaisesRegex(ValueError, "method_calibration.*condition_contract_id"):
                run_stability_screen(
                    "NNO", Path(temp), rate_uncertainty_evidence=[{
                        "route_id": "route-stable", "lower_rate_multiplier": 0.2,
                        "upper_rate_multiplier": 5.0, "source": "synthetic test calibration",
                        "source_kind": "method_calibration",
                        "applicability": "unexecutable gas-only fixture prose",
                    }], method_profile={"orca_keywords": "TEST"},
                )
        collect_orca.assert_not_called()
        collect_rmg.assert_not_called()

    def test_wrong_electron_parity_stops_before_orca_or_rmg(self):
        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence") as collect_orca, \
                patch("storca.stability.collect_rmg_evidence") as collect_rmg:
            with self.assertRaisesRegex(ValueError, "wrong electron parity"):
                run_stability_screen(
                    "[CH3]", Path(temp), multiplicity=1,
                    method_profile={"orca_keywords": "TEST"},
                )
        collect_orca.assert_not_called()
        collect_rmg.assert_not_called()

    def test_opt_in_graph_discovery_passes_only_graph_only_survivors_to_engine(self):
        from storca.reaction_candidates import ReactionCandidate
        from storca.reaction_discovery import DiscoveryResult

        initial = _rmg_evidence(loss=0.001, t95=None, routes=[])
        candidate = ReactionCandidate.create(
            ("[CH3:1][CH3:2]",), ("[CH3:1]", "[CH3:2]"),
            total_charge=0, multiplicity=1,
            discovery_assertions=({"source": "graph_xtb"},),
            screening_status="candidate_unverified", verification_status="candidate_unverified",
            provenance={"atom_mapping_representation": "smiles_atom_map_numbers"},
            reactant_charges=(0,), product_charges=(0, 0),
            reactant_multiplicities=(1,), product_multiplicities=(2, 2),
        )
        graph = DiscoveryResult("graph_xtb", "completed", (candidate,))
        captured = {}

        def fake_engine(*args, **kwargs):
            captured["supplemental"] = kwargs["supplemental_candidate_routes"]
            return {
                "status": "verification_incomplete",
                "reason": "fixture graph candidate retained",
                "condition_contract": args[1],
                "estimated_time_to_retention_seconds": None,
            }

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial), \
                patch("storca.graph_xtb_discovery.GraphXtbDiscoveryBackend.discover", return_value=graph), \
                patch("storca.verification_engine.run_verification_engine", side_effect=fake_engine), \
                patch("storca.verification_adapters.make_orca_route_verifier", return_value=lambda *_: {}), \
                patch("storca.verification_adapters.make_arkane_rate_builder", return_value=lambda *_: {}):
            result = run_stability_screen(
                "CC", Path(temp), graph_xtb_discovery=True,
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )
        self.assertEqual(len(captured["supplemental"]), 1)
        self.assertEqual(captured["supplemental"][0]["kinetic_relevance"], "unknown")
        self.assertEqual(result["reaction_discovery"]["backend"], "merged")
        self.assertEqual(result["reaction_discovery"]["candidates"][0]["discovery_sources"], ["graph_xtb"])

    def test_air_graph_discovery_receives_declared_partner_states(self):
        from storca.reaction_discovery import DiscoveryResult

        initial = _rmg_evidence(loss=0.001, t95=None, routes=[])
        captured = {}

        def discover(system, conditions, budget, **kwargs):
            captured["system"] = system
            return DiscoveryResult(
                "graph_xtb", "completed", (),
                provenance={
                    "target_smiles": "CC",
                    "target_charge": 0,
                    "target_multiplicity": 1,
                    "search_scope": (
                        "single_target_unimolecular_and_condition_specific_"
                        "bimolecular_graph_edit"
                    ),
                    "reaction_classes": [
                        "bimolecular_radical_association",
                        "bimolecular_hydrogen_abstraction",
                    ],
                    "reactive_partner_assessments": [{
                        "label": "oxygen",
                        "status": "completed_declared_search",
                    }],
                },
            )

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial), \
                patch(
                    "storca.graph_xtb_discovery.GraphXtbDiscoveryBackend.discover",
                    side_effect=discover,
                ):
            result = run_stability_screen(
                "CC", Path(temp), scenario="ambient-air-gas-screen",
                graph_xtb_discovery=True, auto_verify_routes=False,
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )
        self.assertEqual(captured["system"]["reactive_partners"], [{
            "label": "oxygen", "smiles": "[O][O]", "charge": 0,
            "multiplicity": 3, "reactive": True,
        }])
        codes = {
            item["code"] for item in result["mechanism_coverage"]["blocking_findings"]
        }
        self.assertNotIn("independent_bimolecular_discovery_scope_missing", codes)

    def test_precomputed_graph_discovery_is_reused_without_rerunning_xtb(self):
        from storca.graph_xtb_discovery import GraphXtbBudget
        from storca.reaction_candidates import ReactionCandidate
        from storca.reaction_discovery import DiscoveryResult

        initial = _rmg_evidence(loss=0.001, t95=None, routes=[])
        candidate = ReactionCandidate.create(
            ("[CH3:1][CH3:2]",), ("[CH3:1]", "[CH3:2]"),
            total_charge=0, multiplicity=1,
            discovery_assertions=({"source": "graph_xtb"},),
            screening_status="candidate_unverified", verification_status="candidate_unverified",
            provenance={"atom_mapping_representation": "smiles_atom_map_numbers"},
            reactant_charges=(0,), product_charges=(0, 0),
            reactant_multiplicities=(1,), product_multiplicities=(2, 2),
        )
        shared = DiscoveryResult(
            "graph_xtb", "completed", (candidate,),
            provenance={
                "target_smiles": "CC", "target_charge": 0, "target_multiplicity": 1,
                "search_scope": "single_target_unimolecular_graph_edit",
                "budget": GraphXtbBudget().as_dict(),
            },
        ).as_dict()
        captured = {}

        def fake_engine(*args, **kwargs):
            captured["supplemental"] = kwargs["supplemental_candidate_routes"]
            return {
                "status": "verification_incomplete",
                "reason": "fixture graph candidate retained",
                "condition_contract": args[1],
                "estimated_time_to_retention_seconds": None,
            }

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial), \
                patch("storca.graph_xtb_discovery.GraphXtbDiscoveryBackend.discover", side_effect=AssertionError("must reuse")), \
                patch("storca.verification_engine.run_verification_engine", side_effect=fake_engine), \
                patch("storca.verification_adapters.make_orca_route_verifier", return_value=lambda *_: {}), \
                patch("storca.verification_adapters.make_arkane_rate_builder", return_value=lambda *_: {}):
            result = run_stability_screen(
                "CC", Path(temp), precomputed_graph_xtb_discovery=shared,
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )
        self.assertEqual(len(captured["supplemental"]), 1)
        self.assertEqual(
            captured["supplemental"][0]["discovery_condition_scope"],
            "target_only_connectivity_and_endpoint_screening",
        )
        self.assertTrue(result["reaction_discovery"]["reuse"]["status"].startswith("reused_precomputed"))

    def test_incomplete_requested_discovery_vetoes_converged_lifetime(self):
        from storca.reaction_discovery import DiscoveryResult

        initial = _rmg_evidence(loss=0.1, t95=None, routes=[])
        truncated = DiscoveryResult(
            "graph_xtb", "discovery_truncated",
            provenance={
                "target_smiles": "CC", "target_charge": 0,
                "target_multiplicity": 1,
            },
        )

        def falsely_converged(*args, **kwargs):
            return {
                "status": "verified_t95_converged",
                "reason": "fixture repaired model convergence",
                "condition_contract": args[1],
                "estimated_time_to_retention_seconds": 12.0,
            }

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial), \
                patch("storca.graph_xtb_discovery.GraphXtbDiscoveryBackend.discover", return_value=truncated), \
                patch("storca.verification_engine.run_verification_engine", side_effect=falsely_converged), \
                patch("storca.verification_adapters.make_orca_route_verifier", return_value=lambda *_: {}), \
                patch("storca.verification_adapters.make_arkane_rate_builder", return_value=lambda *_: {}):
            result = run_stability_screen(
                "CC", Path(temp), graph_xtb_discovery=True,
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )
        verification = result["verification_summary"]
        self.assertEqual(result["reaction_discovery"]["assembly_status"], "discovery_truncated")
        self.assertEqual(verification["status"], "verification_incomplete")
        self.assertEqual(
            verification["rate_network_verification_status"],
            "verified_t95_converged",
        )
        self.assertIsNone(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"])
        self.assertEqual(result["assessment"]["status"], "orca_verification_incomplete")

    def test_incomplete_discovery_blocks_preliminary_no_loss_without_verifier(self):
        from storca.reaction_discovery import DiscoveryResult

        initial = _rmg_evidence(loss=0.001, t95=None, routes=[])
        truncated = DiscoveryResult("graph_xtb", "discovery_truncated")
        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial), \
                patch("storca.graph_xtb_discovery.GraphXtbDiscoveryBackend.discover", return_value=truncated):
            result = run_stability_screen(
                "CC", Path(temp), graph_xtb_discovery=True,
                auto_verify_routes=False,
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )
        self.assertEqual(result["assessment"]["status"], "incomplete_evidence")
        self.assertIn(
            "requested_reaction_discovery_incomplete", result["risk_flags"],
        )
        self.assertIsNone(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"])

    def test_unqueued_graph_survivor_vetoes_model_convergence(self):
        from storca.reaction_candidates import ReactionCandidate
        from storca.reaction_discovery import DiscoveryResult

        initial = _rmg_evidence(loss=0.1, t95=None, routes=[])
        candidate = ReactionCandidate.create(
            ("[CH3:1][CH3:2]",), ("[CH3:1]", "[CH3:2]"),
            total_charge=0, multiplicity=1,
            discovery_assertions=({
                "source": "graph_xtb",
                "edit": {"kind": "single_bond_break"},
            },),
            screening_status="candidate_unverified",
            verification_status="candidate_unverified",
            provenance={"atom_mapping_representation": "smiles_atom_map_numbers"},
            reactant_charges=(0,), product_charges=(0, 0),
            reactant_multiplicities=(1,), product_multiplicities=(2, 2),
        )
        graph = DiscoveryResult(
            "graph_xtb", "completed", (candidate,),
            provenance={
                "reaction_classes": ["single_bond_break"],
                "search_scope": "single_target_unimolecular_graph_edit",
            },
        )

        def falsely_converged(*args, **kwargs):
            self.assertEqual(kwargs["supplemental_candidate_routes"], [])
            return {
                "status": "verified_t95_converged",
                "reason": "fixture retained-model convergence",
                "condition_contract": args[1],
                "estimated_time_to_retention_seconds": 12.0,
                "verified_route_ids": [],
                "deferred_unverified_routes": [],
            }

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial), \
                patch("storca.graph_xtb_discovery.GraphXtbDiscoveryBackend.discover", return_value=graph), \
                patch("storca.verification_engine.run_verification_engine", side_effect=falsely_converged), \
                patch("storca.verification_adapters.make_orca_route_verifier", return_value=lambda *_: {}), \
                patch("storca.verification_adapters.make_arkane_rate_builder", return_value=lambda *_: {}):
            result = run_stability_screen(
                "CC", Path(temp), graph_xtb_discovery=True,
                graph_xtb_budget={"max_candidates_for_verification": 0},
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )
        self.assertEqual(result["verification_summary"]["status"], "verification_incomplete")
        self.assertEqual(
            result["verification_summary"]["rate_network_verification_status"],
            "verified_t95_converged",
        )
        self.assertEqual(
            result["mechanism_coverage"]["status"],
            "unresolved_candidate_kinetics",
        )
        self.assertFalse(result["mechanism_coverage"]["supports_observable_claim"])
        self.assertIsNone(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"])

    def test_uncalibrated_xtb_endpoint_exclusion_vetoes_model_convergence(self):
        from storca.reaction_candidates import ReactionCandidate
        from storca.reaction_discovery import DiscoveryResult

        initial = _rmg_evidence(loss=0.1, t95=None, routes=[])
        candidate = ReactionCandidate.create(
            ("[CH3:1][CH3:2]",), ("[CH3:1]", "[CH3:2]"),
            total_charge=0, multiplicity=1,
            discovery_assertions=({
                "source": "graph_xtb",
                "edit": {"kind": "single_bond_break"},
            },),
            screening_status="xtb_screened_out",
            verification_status="candidate_unverified",
            provenance={"atom_mapping_representation": "smiles_atom_map_numbers"},
            reactant_charges=(0,), product_charges=(0, 0),
            reactant_multiplicities=(1,), product_multiplicities=(2, 2),
        )
        graph = DiscoveryResult(
            "graph_xtb", "completed", (candidate,),
            provenance={
                "reaction_classes": ["single_bond_break"],
                "search_scope": "single_target_unimolecular_graph_edit",
            },
        )

        def falsely_converged(*args, **kwargs):
            self.assertEqual(kwargs["supplemental_candidate_routes"], [])
            return {
                "status": "verified_t95_converged",
                "reason": "fixture retained-model convergence",
                "condition_contract": args[1],
                "estimated_time_to_retention_seconds": 12.0,
                "verified_route_ids": [],
                "deferred_unverified_routes": [],
            }

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial), \
                patch("storca.graph_xtb_discovery.GraphXtbDiscoveryBackend.discover", return_value=graph), \
                patch("storca.verification_engine.run_verification_engine", side_effect=falsely_converged), \
                patch("storca.verification_adapters.make_orca_route_verifier", return_value=lambda *_: {}), \
                patch("storca.verification_adapters.make_arkane_rate_builder", return_value=lambda *_: {}):
            result = run_stability_screen(
                "CC", Path(temp), graph_xtb_discovery=True,
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )

        coverage = result["mechanism_coverage"]
        self.assertEqual(coverage["status"], "unresolved_mechanism_coverage_risk")
        self.assertIn(
            "uncalibrated_endpoint_screen_exclusion",
            {item["code"] for item in coverage["blocking_findings"]},
        )
        self.assertEqual(result["verification_summary"]["status"], "verification_incomplete")
        self.assertEqual(
            result["verification_summary"]["rate_network_verification_status"],
            "verified_t95_converged",
        )
        self.assertIsNone(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"])

    def test_default_tiny_loss_model_stops_before_orca_callbacks(self):
        initial = _rmg_evidence(loss=0.001, t95=None)
        callback_calls = []

        def should_not_run(*args, **kwargs):
            callback_calls.append((args, kwargs))
            raise AssertionError("below-threshold completed model must not start route chemistry")

        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial) as collect_rmg, \
                patch("storca.verification_adapters.make_orca_route_verifier", return_value=should_not_run), \
                patch("storca.verification_adapters.make_arkane_rate_builder", return_value=should_not_run):
            result = run_stability_screen(
                "NNO", Path(temp), target_duration_hours=24.0,
                method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
            )
        self.assertEqual(callback_calls, [])
        self.assertEqual(collect_rmg.call_count, 1)
        self.assertEqual(result["verification_summary"]["status"], "verification_incomplete")
        self.assertEqual(
            result["verification_summary"]["rate_network_verification_status"],
            "no_target_loss_in_completed_rmg_model",
        )
        self.assertEqual(result["assessment"]["status"], "orca_verification_incomplete")
        self.assertIsNone(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"])
        self.assertIsNotNone(result["reaction_discovery"])
        self.assertEqual(
            result["mechanism_coverage"]["status"],
            "unresolved_mechanism_coverage_risk",
        )
        self.assertIn(
            "required_backend_missing",
            {
                item["code"]
                for item in result["mechanism_coverage"]["blocking_findings"]
            },
        )

    def test_default_full_rerun_extracts_manifest_and_matches_canonical_reverse(self):
        initial = _rmg_evidence()
        repaired = _rmg_evidence()
        repaired["generated_kinetics_libraries"] = [{"route_id": "route-stable"}]
        repaired["mechanism_inspection"] = {
            "status": "completed",
            "reactions": [{
                "reaction_index": 0,
                "equation": "P(9)<=>nitrogen(8)+stability(7)",
                "molecularity": 1,
                "reverse_molecularity": 2,
                "evaluated_forward_rate_coefficient_si": 9.0,
                "evaluated_reverse_rate_coefficient_si": 2.0,
                "evaluated_at": {"temperature_K": 298.0, "pressure_bar": 1.0},
                "source_library": None,
                "kinetics_comment": None,
            }],
        }
        verified_calls = []

        def verify(route, workdir):
            verified_calls.append(route["route_id"])
            return {"route_verification": {"status": "verified_transition_state_path"}}

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "run").mkdir()
            verified_library = root / "verified-library"
            base_library = root / "base-library"
            verified_library.mkdir()
            base_library.mkdir()

            def build(route, path_result, workdir):
                return {
                    "status": "completed",
                    "generated_kinetics_library": {
                        "route_id": route["route_id"],
                        "reaction_equation": route["reaction_equation"],
                        "replaces_reaction": route["reaction_equation"],
                        "purpose": "kinetic_repair",
                        "verification_status": "verified_arkane_tst",
                        **rotor_rate_evidence_fixture(Path(workdir) / "statmech-fixture", route_id=route["route_id"]),
                        "molecularity": 2,
                        "condition_rate": {
                            "value": 2.0, "units": "m^3/(mol*s)",
                            "temperature_K": 298.0, "pressure_bar": 1.0,
                        },
                        "collision_limit": {
                            "value": 3.0, "units": "m^3/(mol*s)", "source": "retained ceiling",
                        },
                        "library": str(verified_library),
                    },
                }

            with patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                    patch("storca.stability.collect_rmg_evidence", side_effect=[initial, repaired]) as collect_rmg, \
                    patch("storca.verification_adapters.make_orca_route_verifier", return_value=verify), \
                    patch("storca.verification_adapters.make_arkane_rate_builder", return_value=build):
                result = run_stability_screen(
                    "NNO", root / "run", temperature=298.0, target_duration_hours=24.0,
                    reaction_libraries=[base_library],
                    method_profile={"orca_keywords": "TEST", "harmonic_scale_factor": 1.0},
                )

            rerun_call = collect_rmg.call_args_list[1]
            self.assertEqual(
                rerun_call.kwargs["reaction_libraries"], [verified_library, base_library],
            )
            self.assertEqual(rerun_call.kwargs["conditions"].as_dict(), result["condition_contract"])
        self.assertEqual(verified_calls, ["route-stable"])
        self.assertEqual(result["verification_summary"]["status"], "verification_incomplete")
        self.assertEqual(
            result["verification_summary"]["rate_network_verification_status"],
            "verified_t95_converged",
        )
        self.assertEqual(result["verification_summary"]["condition_contract"], result["condition_contract"])
        self.assertEqual(result["assessment"]["status"], "orca_verification_incomplete")
        self.assertIsNone(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"])
        self.assertEqual(
            result["mechanism_coverage"]["status"],
            "unresolved_mechanism_coverage_risk",
        )

    def test_invalid_parent_minimum_cannot_be_overwritten_by_model_no_loss(self):
        initial = _rmg_evidence(loss=0.0, t95=None, routes=[])
        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca(minimum=False)), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial):
            result = run_stability_screen(
                "NNO", Path(temp), method_profile={"orca_keywords": "TEST"},
            )
        self.assertEqual(result["assessment"]["status"], "incomplete_evidence")
        self.assertEqual(result["verification_summary"]["status"], "verification_incomplete")
        self.assertIn("not_a_validated_local_minimum", result["verification_summary"]["reason"])

    def test_incomplete_collision_bounds_do_not_report_model_no_loss_when_verification_disabled(self):
        initial = _rmg_evidence(
            loss=0.001,
            t95=None,
            kinetics={"status": "incomplete", "violation_count": 0},
        )
        with tempfile.TemporaryDirectory() as temp, \
                patch("storca.stability.collect_orca_evidence", return_value=_orca()), \
                patch("storca.stability.collect_rmg_evidence", return_value=initial):
            result = run_stability_screen(
                "NNO", Path(temp), auto_verify_routes=False,
                method_profile={"orca_keywords": "TEST"},
            )
        self.assertEqual(result["assessment"]["status"], "incomplete_evidence")


class StableIdentityAndDeduplicationTests(unittest.TestCase):
    def test_reverse_network_copy_does_not_make_replacement_rebind_ambiguous(self):
        conditions = {"temperature_K": 298.15, "pressure_bar": 1.0}
        evidence = {
            "candidate_routes": [{
                "route_id": "current-id",
                "reaction_equation": "stability(8)+A(2)<=>P(4)",
            }],
            "network_routes": [{
                "route_id": "reverse-network-copy",
                "reaction_equation": "P(9)<=>A(7)+stability(3)",
            }],
            "mechanism_inspection": {"reactions": [{
                "reaction_index": 0,
                "equation": "stability(8)+A(2)<=>P(4)",
                "molecularity": 2,
                "evaluated_forward_rate_coefficient_si": 4.0,
                "evaluated_at": conditions,
            }]},
        }
        replacement = {
            "route_id": "old-id",
            "reaction_equation": "A(1)+stability(1)<=>P(2)",
            "molecularity": 2,
            "condition_rate": {"value": 4.0, "units": "m^3/(mol*s)", **conditions},
        }
        rebound = _planning_replacements([replacement], evidence, conditions)
        self.assertEqual(rebound[0]["route_id"], "current-id")
        self.assertEqual(rebound[0]["planning_applied_rate_validation"]["status"], "passed")
        self.assertEqual(replacement["route_id"], "old-id")
        # Legacy label deduplication is not permission to skip the new rate
        # evidence gate, even when the route labels otherwise match exactly.
        del evidence["mechanism_inspection"]
        unresolved = _planning_replacements([replacement], evidence, conditions)
        self.assertEqual(unresolved[0]["route_id"], "old-id")
        self.assertEqual(unresolved[0]["planning_applied_rate_validation"]["status"], "failed")

    def test_incomplete_solver_evidence_cannot_use_below_threshold_shortcut(self):
        evidence = _rmg_evidence(loss=0.001, t95=None)
        evidence["time_coverage"] = {"complete": False}

        def fail_closed(*args, **kwargs):
            raise AssertionError("no chemistry should run after the solver gate fails")

        contract = {
            "temperature_K": 298.15,
            "pressure_bar": 1.0,
            "composition": {"stability": 0.01, "nitrogen": 0.99},
            "target_duration_seconds": 86400.0,
            "retention_fraction": 0.95,
        }
        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                evidence, contract, Path(temp),
                verify_route=fail_closed, build_rate=fail_closed, rerun_model=fail_closed,
            )
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("insufficient_solver_evidence", summary["reason"])
        self.assertEqual(summary["claim_scope"], "verification_incomplete_no_condition_lifetime")


if __name__ == "__main__":
    unittest.main()
