import json
import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from propagation_fixtures import isothermal_model
from arkane_rotor_fixtures import rotor_rate_evidence_fixture

from storca.verification_engine import (
    VerificationEngine,
    VerificationEngineConfig,
    continue_verification_engine,
    run_verification_engine,
    _planning_replacements,
)
from storca.flux_verification import build_flux_ranked_verification_plan


CONTRACT = {
    "contract_id": "sha256:verification-engine-fixture-condition",
    "scenario": "ambient-air-gas-screen",
    "phase_model": "homogeneous gas-phase surrogate",
    "temperature_K": 298.15,
    "pressure_bar": 1.0,
    "composition": {"stability": 0.01, "A": 0.99},
    "target_label": "stability",
    "target_duration_seconds": 1000.0,
    "retention_fraction": 0.95,
    "light_condition": "dark",
}


ROUTES = [
    {
        "route_id": "init", "reaction_equation": "stability(1) + A(2) => R(3) + P1(4)",
        "initiation_status": "directly_present",
        "co_reactants": [{"label": "A(2)", "origin": "initial_scenario_species"}],
        "resolved_endpoints": {
            "reactants": [{"label": "stability(1)"}, {"label": "A(2)"}],
            "products": [{"label": "R(3)"}, {"label": "P1(4)"}],
        },
    },
    {
        "route_id": "main", "reaction_equation": "stability(1) + R(3) => P2(5)",
        "initiation_status": "thermally_reachable",
        "co_reactants": [{"label": "R(3)", "origin": "generated_intermediate"}],
        "resolved_endpoints": {
            "reactants": [{"label": "stability(1)"}, {"label": "R(3)"}],
            "products": [{"label": "P2(5)"}],
        },
    },
]


def _propagation(*, t95=50.0, loss=0.20, weights=(1.0, 9.0)):
    final = 1.0 - loss
    reactions = [
        {"reaction_index": index, "reaction_equation": route["reaction_equation"],
         "integrated_target_destruction_kmol": float(weight)}
        for index, (route, weight) in enumerate(zip(ROUTES, weights)) if weight > 0.0
    ]
    return {
        "status": "completed",
        "coverage": {"complete": True, "requested_seconds": 1000.0, "simulated_seconds": 1000.0},
        "estimated_time_to_retention_seconds": t95,
        "reactor_model": isothermal_model(),
        "target_loss_fraction": loss,
        "target_profile": [
            {"time_seconds": 0.0, "target_fraction_remaining": 1.0},
            {"time_seconds": 1000.0, "target_fraction_remaining": final},
        ],
        "reaction_flux_attribution": {
            "status": "completed",
            "total_integrated_target_destruction_kmol": float(sum(weights)),
            "numerical_closure_relative_error": 0.01,
            "reactions": reactions,
        },
    }


def _replacement(route, library, *, molecularity=2):
    replacement = {
        "route_id": route["route_id"],
        "reaction_equation": route["reaction_equation"],
        "verification_status": "verified_test_rate",
        "molecularity": molecularity,
        "condition_rate": {
            "value": 2.0,
            "units": "m^3/(mol*s)" if molecularity == 2 else "s^-1",
            "temperature_K": 298.15, "pressure_bar": 1.0,
        },
        "library": str(library),
    }
    if molecularity == 2:
        replacement["collision_limit"] = {
            "value": 3.0, "units": "m^3/(mol*s)",
            "source": "retained_test_ceiling",
        }
    return replacement


def _model(applied=(), *, weights=(1.0, 9.0), t95=50.0, loss=0.20,
           applied_rate=2.0, envelope=None, routes=None):
    routes = ROUTES if routes is None else routes
    propagation = _propagation(t95=t95, loss=loss, weights=weights)
    if not routes:
        propagation["reaction_flux_attribution"].update(
            total_integrated_target_destruction_kmol=0.0, reactions=[]
        )
    evidence = {
        "status": "completed",
        "execution_assessment": {"status": "completed"},
        "time_coverage": {"complete": True},
        "candidate_routes": routes,
        "independent_cantera_propagation": propagation,
        "kinetics_validation": {"status": "passed", "violation_count": 0},
        "generated_kinetics_libraries": [{"route_id": route_id} for route_id in applied],
        "mechanism_inspection": {
            "status": "completed",
            "reactions": [
                {
                    "reaction_index": index,
                    "equation": next(route["reaction_equation"] for route in ROUTES if route["route_id"] == route_id),
                    "molecularity": 2,
                    "evaluated_forward_rate_coefficient_si": applied_rate,
                    "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
                    "source_library": None,
                    "kinetics_comment": None,
                }
                for index, route_id in enumerate(applied)
            ],
        },
    }
    wrapper = {"condition_contract": CONTRACT, "rmg_evidence": evidence}
    if envelope is not None:
        wrapper["unresolved_sensitivity_envelope"] = envelope
    return wrapper


def _bound_propagation(t95, loss=0.20):
    return {
        "status": "completed",
        "coverage": {"complete": True},
        "estimated_time_to_retention_seconds": t95,
        "target_loss_fraction": loss,
        "target_profile": [{"time_seconds": 1000.0, "target_fraction_remaining": 1.0 - loss}],
    }


class VerificationEngineTests(unittest.TestCase):
    @staticmethod
    def _uncertainty(route_id):
        return {
            "route_id": route_id,
            "lower_rate_multiplier": 0.5,
            "upper_rate_multiplier": 2.0,
            "source": "test calibration",
            "source_kind": "method_calibration",
            "applicability": "test fixture chemistry",
            "condition_contract_id": CONTRACT["contract_id"],
        }

    @patch("storca.reaction_correspondence.resolve_replacement_correspondence")
    def test_planning_reuses_correspondence_and_preserves_source_identity_and_rate_orientation(self, resolve):
        source = _replacement(
            {"route_id": "native-original", "reaction_equation": "parent => product1 + product2"},
            "verified-library", molecularity=1,
        )
        source["native_candidate_id"] = "native-original"
        original = copy.deepcopy(source)
        source_oriented = "stability(9) <=> fragment(4) + resonance_form(8)"
        candidate = {
            "reaction_index": 17, "equation": "fragment(4) + resonance_form(8) <=> stability(9)",
            "molecularity": 2, "reverse_molecularity": 1,
            "evaluated_forward_rate_coefficient_si": 9.0,
            "evaluated_reverse_rate_coefficient_si": 2.0,
            "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
        }
        route = {
            "route_id": "model-route-3", "reaction_equation": source_oriented,
            "retained_rate_evidence": {"reaction_index": 17, "direction": "reverse"},
            "initiation_status": "directly_present",
        }
        evidence = {
            "status": "completed", "execution_assessment": {"status": "completed"},
            "time_coverage": {"complete": True},
            "kinetics_validation": {"status": "passed", "violation_count": 0},
            "candidate_routes": [route],
            "mechanism_inspection": {"status": "completed", "reactions": [candidate]},
            "independent_cantera_propagation": {
                **_propagation(weights=(10.0, 0.0)),
                "reaction_flux_attribution": {
                    "status": "completed", "total_integrated_target_destruction_kmol": 10.0,
                    "numerical_closure_relative_error": 0.0,
                    "reactions": [{
                        "reaction_index": 17, "reaction_equation": source_oriented,
                        "integrated_target_destruction_kmol": 10.0,
                    }],
                },
            },
        }
        correspondence = {"source_route_id": "native-original", "reaction_index": 17}
        resolve.return_value = {
            "status": "matched", "candidate": candidate, "applied_direction": "reverse",
            "source_oriented_equation": source_oriented, "source_route_id": "native-original",
            "correspondence": correspondence, "blocking_reasons": [],
        }
        planned = _planning_replacements([source], evidence, CONTRACT)
        resolve.assert_called_once_with(evidence, source)
        self.assertEqual(source, original)
        self.assertEqual(planned[0]["route_id"], "model-route-3")
        self.assertEqual(planned[0]["identity_source_route_id"], "native-original")
        self.assertEqual(planned[0]["native_candidate_id"], "native-original")
        self.assertEqual(planned[0]["source_reaction_equation"], original["reaction_equation"])
        self.assertEqual(planned[0]["reaction_equation"], source_oriented)
        self.assertEqual(planned[0]["condition_rate"], original["condition_rate"])
        self.assertEqual(planned[0]["molecularity"], 1)
        self.assertEqual(planned[0]["planning_applied_rate_validation"]["applied_direction"], "reverse")
        plan = build_flux_ranked_verification_plan(evidence, CONTRACT, verified_replacements=planned)
        self.assertTrue(next(row for row in plan["ranked_routes"] if row["route_id"] == "model-route-3")["verified_replacement"])
        candidate["evaluated_reverse_rate_coefficient_si"] = 8.0
        invalid = _planning_replacements([source], evidence, CONTRACT)
        self.assertEqual(invalid[0]["planning_applied_rate_validation"]["status"], "failed")
        invalid_plan = build_flux_ranked_verification_plan(evidence, CONTRACT, verified_replacements=invalid)
        self.assertFalse(next(row for row in invalid_plan["ranked_routes"] if row["route_id"] == "model-route-3")["verified_replacement"])

    @patch("storca.reaction_correspondence.resolve_replacement_correspondence")
    def test_planning_failed_correspondence_cannot_credit_even_same_route_id_and_labels(self, resolve):
        source = _replacement(ROUTES[0], "verified-library")
        original = copy.deepcopy(source)
        evidence = _model(applied=("init",), weights=(9.0, 1.0))["rmg_evidence"]
        resolve.return_value = {
            "status": "failed", "blocking_reasons": ["correspondence_artifact_hash_mismatch"],
        }
        planned = _planning_replacements([source], evidence, CONTRACT)
        plan = build_flux_ranked_verification_plan(evidence, CONTRACT, verified_replacements=planned)
        self.assertFalse(next(row for row in plan["ranked_routes"] if row["route_id"] == "init")["verified_replacement"])
        self.assertEqual(source, original)
        self.assertEqual(planned[0]["verification_status"], original["verification_status"])
        self.assertEqual(planned[0]["planning_applied_rate_validation"]["status"], "failed")

    @patch("storca.reaction_correspondence.resolve_replacement_correspondence")
    def test_planning_structural_match_does_not_guess_between_parallel_route_views(self, resolve):
        source = _replacement(ROUTES[0], "verified-library")
        evidence = _model(applied=("init",))["rmg_evidence"]
        evidence["candidate_routes"] = [
            {**ROUTES[0], "route_id": "parallel-a"},
            {**ROUTES[0], "route_id": "parallel-b"},
        ]
        candidate = evidence["mechanism_inspection"]["reactions"][0]
        resolve.return_value = {
            "status": "matched", "candidate": candidate, "applied_direction": "forward",
            "source_oriented_equation": candidate["equation"], "source_route_id": source["route_id"],
            "correspondence": {}, "blocking_reasons": [],
        }
        planned = _planning_replacements([source], evidence, CONTRACT)
        self.assertEqual(planned[0]["route_id"], source["route_id"])
        self.assertIn("planning_model_route_missing_or_ambiguous", planned[0]["planning_applied_rate_validation"]["blocking_reasons"])

        # An explicit retained model index can distinguish the route views;
        # a matching equation alone cannot choose a parallel channel.
        evidence["candidate_routes"][0]["retained_rate_evidence"] = {"reaction_index": 0}
        evidence["candidate_routes"][1]["retained_rate_evidence"] = {"reaction_index": 9}
        planned = _planning_replacements([source], evidence, CONTRACT)
        self.assertEqual(planned[0]["route_id"], "parallel-a")
        self.assertEqual(planned[0]["planning_applied_rate_validation"]["status"], "passed")

    def test_planning_legacy_wrong_rate_gets_no_credit_for_matching_route_id(self):
        source = _replacement(ROUTES[0], "verified-library")
        evidence = _model(applied=("init",), weights=(9.0, 1.0), applied_rate=8.0)["rmg_evidence"]
        planned = _planning_replacements([source], evidence, CONTRACT)
        self.assertEqual(planned[0]["planning_applied_rate_validation"]["reaction_identity_basis"], "legacy_normalized_equation_labels")
        self.assertIn("applied_rate_does_not_match_replacement", planned[0]["planning_applied_rate_validation"]["blocking_reasons"])
        plan = build_flux_ranked_verification_plan(evidence, CONTRACT, verified_replacements=planned)
        self.assertFalse(next(row for row in plan["ranked_routes"] if row["route_id"] == "init")["verified_replacement"])

    def test_uncertainty_evidence_is_persisted_and_dependency_order_is_preserved(self):
        initial = _model(weights=(9.0, 1.0))["rmg_evidence"]
        initial["independent_cantera_propagation"]["reaction_sensitivity_attribution"] = {
            "status": "completed",
            "reactions": [
                {"reaction_equation": ROUTES[0]["reaction_equation"],
                 "absolute_log_t95_sensitivity": 0.01},
                {"reaction_equation": ROUTES[1]["reaction_equation"],
                 "absolute_log_t95_sensitivity": 10.0},
            ],
        }
        calls = []

        def unresolved(route, workdir):
            calls.append(route["route_id"])
            return {"route_verification": {"status": "surface_unresolved"}}

        def should_not_run(*args):
            raise AssertionError("An unresolved path cannot create a rate")

        with tempfile.TemporaryDirectory() as temp:
            engine = VerificationEngine(
                initial, CONTRACT, Path(temp), verify_route=unresolved,
                build_rate=should_not_run, rerun_model=should_not_run,
                rate_uncertainties=[self._uncertainty("init"), self._uncertainty("main")],
            )
            summary = engine.run()
            first_plan = json.loads((engine.session / "plan-001.json").read_text())
            state = json.loads((engine.session / "engine-state.json").read_text())
        self.assertEqual(first_plan["allocation_policy"]["status"], "observable_uncertainty_priority_active")
        self.assertEqual(first_plan["controlling_route_id"], "main")
        self.assertEqual(first_plan["selected_route_id"], "init")
        self.assertEqual(calls[0], "init")
        self.assertEqual(state["rate_uncertainty_evidence"]["record_count"], 2)
        self.assertEqual(summary["status"], "verification_incomplete")

    def test_malformed_uncertainty_fails_before_session_directory_is_created(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with self.assertRaisesRegex(ValueError, "applicability"):
                VerificationEngine(
                    _model()["rmg_evidence"], CONTRACT, root,
                    verify_route=lambda *_: {}, build_rate=lambda *_: {}, rerun_model=lambda *_: {},
                    rate_uncertainties=[{
                        "route_id": "init",
                        "lower_rate_multiplier": 0.5,
                        "upper_rate_multiplier": 2.0,
                        "source": "test",
                        "source_kind": "user_supplied",
                    }],
                )
            self.assertEqual(list(root.iterdir()), [])

    def test_unbound_calibration_fails_before_session_or_callbacks(self):
        supplied = self._uncertainty("init")
        supplied.pop("condition_contract_id")
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with self.assertRaisesRegex(ValueError, "method_calibration.*condition_contract_id"):
                VerificationEngine(
                    _model()["rmg_evidence"], CONTRACT, root,
                    verify_route=lambda *_: self.fail("unbound calibration reached path"),
                    build_rate=lambda *_: self.fail("unbound calibration reached rate"),
                    rerun_model=lambda *_: self.fail("unbound calibration reached propagation"),
                    rate_uncertainties=[supplied],
                )
            self.assertEqual(list(root.iterdir()), [])

    @patch("storca.rmg_bridge_client.run_rmg_bridge")
    def test_engine_generates_missing_sensitivity_then_replans(self, bridge):
        bridge.side_effect = [
            _propagation(t95=105.0), _propagation(t95=95.0),
            _propagation(t95=150.0), _propagation(t95=50.0),
        ]
        initial = _model(weights=(9.0, 1.0))["rmg_evidence"]
        initial["artifacts"] = {
            "chemkin": "chem.inp",
            "species_dictionary": "dictionary.txt",
        }
        calls = []

        def unresolved(route, workdir):
            calls.append(route["route_id"])
            return {"route_verification": {"status": "surface_unresolved"}}

        def should_not_run(*args):
            raise AssertionError("An unresolved path cannot create a rate")

        with tempfile.TemporaryDirectory() as temp:
            engine = VerificationEngine(
                initial, CONTRACT, Path(temp), verify_route=unresolved,
                build_rate=should_not_run, rerun_model=should_not_run,
                rate_uncertainties=[self._uncertainty("init"), self._uncertainty("main")],
                rmg_env="rmg_env",
            )
            engine.run()
            attribution = json.loads(
                (engine.session / "t95-sensitivity-001.json").read_text()
            )
            plan = json.loads((engine.session / "plan-001.json").read_text())
        self.assertEqual(bridge.call_count, 4)
        self.assertEqual(attribution["status"], "completed")
        self.assertEqual(plan["allocation_policy"]["status"], "observable_uncertainty_priority_active")
        self.assertEqual(plan["controlling_route_id"], "main")
        self.assertEqual(calls[0], "init")

    def test_failed_optional_sensitivity_diagnostic_falls_back_without_aborting(self):
        initial = _model()["rmg_evidence"]

        def unresolved(route, workdir):
            return {"route_verification": {"status": "surface_unresolved"}}

        def should_not_run(*args):
            raise AssertionError("An unresolved path cannot create a rate")

        with tempfile.TemporaryDirectory() as temp:
            engine = VerificationEngine(
                initial, CONTRACT, Path(temp), verify_route=unresolved,
                build_rate=should_not_run, rerun_model=should_not_run,
                rate_uncertainties=[self._uncertainty("init"), self._uncertainty("main")],
            )
            summary = engine.run()
            attribution = json.loads(
                (engine.session / "t95-sensitivity-001.json").read_text()
            )
            plan = json.loads((engine.session / "plan-001.json").read_text())
        self.assertEqual(attribution["status"], "incomplete")
        self.assertIn("mechanism artifacts are missing", attribution["reason"])
        self.assertEqual(plan["allocation_policy"]["status"], "fallback_incomplete_evidence")
        self.assertEqual(summary["status"], "verification_incomplete")

    def test_graph_only_candidate_uses_existing_path_gate_and_cannot_restore_no_loss_claim(self):
        evidence = _model(routes=[], loss=0.001, t95=None)
        calls = []
        candidate_route = {
            "native_candidate_id": "graph-route-1",
            "route_id": "graph-route-1",
            "source": "graph_xtb",
            "reaction_equation": "stability <=> graph-product",
            "kinetic_relevance": "unknown",
            "screening_status": "candidate_unverified",
            "verification_queue_class": "graph_xtb_screened_kinetic_relevance_unknown",
        }

        def verify(route, workdir):
            calls.append(route["route_id"])
            return {"route_verification": {"status": "mapping_unresolved"}}

        def should_not_build_or_rerun(*args, **kwargs):
            raise AssertionError("a deferred graph path must not produce a rate or rerun")

        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                evidence, CONTRACT, Path(temp),
                verify_route=verify,
                build_rate=should_not_build_or_rerun,
                rerun_model=should_not_build_or_rerun,
                supplemental_candidate_routes=[candidate_route],
            )
        self.assertEqual(calls, ["graph-route-1"])
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("supplemental_graph_candidates_deferred", summary["reason"])
        self.assertEqual(summary["claim_scope"], "verification_incomplete_no_condition_lifetime")

    def test_deferred_dissociation_retains_physical_rate_uncertainty_bound(self):
        evidence = _model(routes=[], loss=0.001, t95=None)
        route = {
            "route_id": "bounded-dissociation",
            "reaction_equation": "stability <=> R1 + R2",
            "reactant_smiles": ["CC"],
            "product_smiles": ["[CH3]", "[CH3]"],
            "kinetic_relevance": "unknown",
            "screening_status": "candidate_unverified",
            "retained_rate_evidence": {
                "rate_coefficient_si": 2.0, "rate_units": "s^-1",
            },
        }

        def verify(candidate, workdir):
            return {"route_verification": {
                "status": "dissociation_rate_model_required",
                "path_classification": "submerged_saddle_dissociation_candidate",
                "forward_loss_upper_bound": {
                    "applicable": True,
                    "bound_basis": "fixture_physical_bound",
                    "forward_rate_upper_bound_si": 20.0,
                    "forward_rate_upper_bound_units": "s^-1",
                    "assumptions": ["fixture"],
                },
            }}

        with tempfile.TemporaryDirectory() as temp:
            engine = VerificationEngine(
                evidence, CONTRACT, Path(temp),
                verify_route=verify,
                build_rate=lambda *_: self.fail("deferred route cannot build a rate"),
                rerun_model=lambda *_: self.fail("deferred route cannot rerun"),
                supplemental_candidate_routes=[route],
            )
            engine.run()
            state = json.loads((engine.session / "engine-state.json").read_text())

        retained = state["rate_uncertainty_evidence"]["records"][0]
        self.assertEqual(retained["route_id"], "bounded-dissociation")
        self.assertEqual(retained["uncertainty_status"], "unbounded_downward")
        self.assertEqual(retained["upper_rate_multiplier"], 10.0)
        self.assertEqual(
            state["deferred_unverified_routes"]["bounded-dissociation"]
            ["rate_uncertainty"]["source_kind"],
            "physical_bound",
        )

    def test_partial_rmg_exhaustion_attempts_direct_graph_path_without_claiming_lifetime(self):
        for partner, condition_id, expected in (
            ("A", None, ["graph-direct"]),
            ("generated-radical", None, []),
            ("A", "another-condition", []),
        ):
            with self.subTest(partner=partner, condition=condition_id):
                initial = _model(routes=[], loss=0.001, t95=None)["rmg_evidence"]
                initial.update(status="incomplete", enumeration_status="partial_enumeration")
                route = {
                    "route_id": "graph-direct",
                    "reaction_equation": f"stability + {partner} <=> product",
                    "source": "graph_xtb",
                    "screening_status": "candidate_unverified",
                    "atom_balance": {"status": "passed"},
                    "resolved_endpoints": {
                        "reactants": [{"label": "stability"}, {"label": partner}],
                        "products": [{"label": "product"}],
                    },
                    "verification_condition_contract_id": condition_id,
                }
                calls = []

                def verify(candidate, workdir):
                    calls.append(candidate["route_id"])
                    return {"route_verification": {"status": "mapping_unresolved"}}

                with tempfile.TemporaryDirectory() as temp:
                    engine = VerificationEngine(
                        initial, CONTRACT, Path(temp), verify_route=verify,
                        build_rate=lambda *_: self.fail("Unverified path cannot produce a rate"),
                        rerun_model=lambda *_: self.fail("Unverified path cannot repair the network"),
                        supplemental_candidate_routes=[route],
                    )
                    summary = engine.run()
                    if expected:
                        plan = json.loads((engine.session / "supplemental-plan-001.json").read_text())
                        self.assertTrue(plan["partial_enumeration"])
                        self.assertFalse(plan["observable_claim_allowed"])
                self.assertEqual(calls, expected)
                self.assertEqual(summary["status"], "verification_incomplete")
                self.assertIsNone(summary.get("estimated_time_to_retention_seconds"))

    def test_graph_only_candidate_verified_rate_is_inserted_and_repropagated(self):
        initial = _model(routes=[], loss=0.001, t95=None)["rmg_evidence"]
        candidate_route = {
            "native_candidate_id": "reaction-graph-verified",
            "route_id": "reaction-graph-verified",
            "source": "graph_xtb",
            "discovery_sources": ["graph_xtb"],
            "reaction_equation": "stability <=> methyl_1 + methyl_2",
            "reactant_smiles": ["CC"],
            "product_smiles": ["[CH3]", "[CH3]"],
            "reactant_labels": ["stability"],
            "product_labels": ["methyl_1", "methyl_2"],
            "kinetic_relevance": "unknown",
            "screening_status": "candidate_unverified",
            "verification_queue_class": "graph_xtb_screened_kinetic_relevance_unknown",
        }
        rerun_calls = []

        def verify(route, workdir):
            return {"route_verification": {"status": "verified_transition_state_path"}}

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            library = root / "verified-graph-library"

            def build_rate(route, path_result, workdir):
                library.mkdir()
                return {
                    "rate_verification_status": "validated",
                    "replacement": _replacement(route, library, molecularity=1),
                }

            def rerun(library_paths, workdir):
                rerun_calls.append(list(library_paths))
                repaired = _model(routes=[], loss=0.001, t95=None)
                repaired_evidence = repaired["rmg_evidence"]
                repaired_evidence["generated_kinetics_libraries"] = [{
                    "route_id": candidate_route["route_id"],
                }]
                repaired_evidence["mechanism_inspection"]["reactions"] = [{
                    "reaction_index": 0,
                    "equation": candidate_route["reaction_equation"],
                    "molecularity": 1,
                    "evaluated_forward_rate_coefficient_si": 2.0,
                    "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
                    "source_library": None,
                    "kinetics_comment": None,
                }]
                return repaired

            engine = VerificationEngine(
                initial, CONTRACT, root,
                verify_route=verify, build_rate=build_rate, rerun_model=rerun,
                supplemental_candidate_routes=[candidate_route],
                config=VerificationEngineConfig(max_iterations=1),
            )
            summary = engine.run()
            state = json.loads((engine.session / "engine-state.json").read_text())
            handoff_path = Path(
                state["iterations"][0]["artifacts"]["native_candidate_kinetic_handoff"]
            )
            handoff = json.loads(handoff_path.read_text())

        self.assertEqual(rerun_calls, [[library]])
        self.assertEqual(handoff["status"], "verified_rate_inserted_and_repropagated")
        self.assertEqual(handoff["native_candidate_id"], "reaction-graph-verified")
        self.assertEqual(handoff["discovery_sources"], ["graph_xtb"])
        self.assertEqual(
            state["verified_replacements"][0]["native_candidate_id"],
            "reaction-graph-verified",
        )
        self.assertEqual(
            summary["native_candidate_handoffs"][0]["status"],
            "verified_rate_inserted_and_repropagated",
        )

    def test_specialist_dissociation_rate_enters_same_repair_loop(self):
        import hashlib

        initial = _model(routes=[], loss=0.001, t95=None)["rmg_evidence"]
        route = {
            "route_id": "specialist-dissociation",
            "reaction_equation": "stability <=> radical-1 + radical-2",
            "reactant_smiles": ["CC"],
            "product_smiles": ["[CH3]", "[CH3]"],
            "reactant_labels": ["stability"],
            "product_labels": ["radical-1", "radical-2"],
            "charge": 0,
            "multiplicity": 1,
            "kinetic_relevance": "unknown",
            "screening_status": "candidate_unverified",
            "verification_queue_class": "specialist-fixture",
        }
        normal_rate_calls = []
        rerun_calls = []

        def verify(candidate, workdir):
            return {"route_verification": {
                "status": "dissociation_rate_model_required",
                "path_classification": "submerged_saddle_dissociation_candidate",
                "method_selection": {"ordinary_harmonic_tst_applicable": False},
            }}

        def normal_rate(*args):
            normal_rate_calls.append(args)
            raise AssertionError("Harmonic-TST builder must not receive a dissociation escalation")

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            library = root / "specialist-library"
            evidence_file = root / "master-equation.out"

            def specialist(request, workdir):
                library.mkdir()
                evidence_file.write_text("verified fixture output\n")
                reactions = library / "reactions.py"
                dictionary = library / "dictionary.txt"
                reactions.write_text("# fixture specialist kinetics\n")
                dictionary.write_text("fixture\n")
                return {
                    "schema_version": 1,
                    "kind": "storca_specialist_rate_result",
                    "status": "verified",
                    "request_sha256": request["request_sha256"],
                    "condition_contract_sha256": request["condition_contract_sha256"],
                    "rate_theory": "rrkm_master_equation",
                    "provider": {"name": "fixture", "version": "1.0"},
                    "applicability": {"status": "applicable", "diagnostics": ["fixture"]},
                    "uncertainty": {
                        "lower_rate_multiplier": 0.1, "upper_rate_multiplier": 10.0,
                        "source": "fixture", "source_kind": "method_calibration",
                        "applicability": "fixture dissociation",
                        "uncertainty_components": ["energy_transfer", "rate_theory"],
                    },
                    "pressure_dependence": {
                        "pressure_range_bar": [0.01, 10.0],
                        "bath_gas": {"A": 1.0},
                        "energy_transfer_model": "fixture_delta_down",
                    },
                    "evidence_artifacts": [{
                        "path": str(evidence_file),
                        "sha256": hashlib.sha256(evidence_file.read_bytes()).hexdigest(),
                    }],
                    "generated_kinetics_library": {
                        "route_id": request["route_id"],
                        "reaction_equation": request["reaction_equation"],
                        "replaces_reaction": request["reaction_equation"],
                        "verification_status": "verified_rrkm_master_equation",
                        "molecularity": 1,
                        "condition_rate": {
                            "value": 2.0, "units": "s^-1",
                            "temperature_K": 298.15, "pressure_bar": 1.0,
                        },
                        "collision_limit": {"applicable": False},
                        "library": str(library),
                        "checksums": {
                            "reactions.py": hashlib.sha256(reactions.read_bytes()).hexdigest(),
                            "dictionary.txt": hashlib.sha256(dictionary.read_bytes()).hexdigest(),
                        },
                    },
                }

            def rerun(library_paths, workdir):
                rerun_calls.append(list(library_paths))
                repaired = _model(routes=[], loss=0.001, t95=None)
                evidence = repaired["rmg_evidence"]
                evidence["generated_kinetics_libraries"] = [{"route_id": route["route_id"]}]
                evidence["mechanism_inspection"]["reactions"] = [{
                    "reaction_index": 0,
                    "equation": route["reaction_equation"],
                    "molecularity": 1,
                    "evaluated_forward_rate_coefficient_si": 2.0,
                    "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
                    "source_library": None,
                    "kinetics_comment": None,
                }]
                return repaired

            engine = VerificationEngine(
                initial, CONTRACT, root,
                verify_route=verify, build_rate=normal_rate, rerun_model=rerun,
                build_specialist_rate=specialist,
                supplemental_candidate_routes=[route],
                config=VerificationEngineConfig(max_iterations=1),
            )
            engine.run()
            state = json.loads((engine.session / "engine-state.json").read_text())
            iteration = state["iterations"][0]
            validation = json.loads(
                Path(iteration["artifacts"]["specialist_rate_validation"]).read_text()
            )

        self.assertEqual(normal_rate_calls, [])
        self.assertEqual(rerun_calls, [[library]])
        self.assertEqual(validation["status"], "passed")
        self.assertEqual(iteration["path_verification_status"], "verified_specialist_rate_model")
        self.assertEqual(
            state["rate_uncertainty_evidence"]["records"][0]["route_id"],
            "specialist-dissociation",
        )

    def test_specialist_continuation_reuses_retained_path_without_orca_retry(self):
        import hashlib

        initial = _model(routes=[], loss=0.001, t95=None)["rmg_evidence"]
        route = {
            "route_id": "retained-dissociation",
            "reaction_equation": "stability <=> radical-1 + radical-2",
            "reactant_smiles": ["CC"], "product_smiles": ["[CH3]", "[CH3]"],
            "reactant_labels": ["stability"],
            "product_labels": ["radical-1", "radical-2"],
            "charge": 0, "multiplicity": 1,
            "kinetic_relevance": "unknown", "screening_status": "candidate_unverified",
            "verification_queue_class": "specialist-fixture",
        }
        path_calls = []

        def first_path(candidate, workdir):
            path_calls.append(candidate["route_id"])
            return {"route_verification": {
                "status": "dissociation_rate_model_required",
                "path_classification": "submerged_saddle_dissociation_candidate",
                "method_selection": {"ordinary_harmonic_tst_applicable": False},
            }}

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = VerificationEngine(
                initial, CONTRACT, root,
                verify_route=first_path, build_rate=lambda *_: {}, rerun_model=lambda *_: {},
                supplemental_candidate_routes=[route],
                config=VerificationEngineConfig(max_iterations=1),
                path_verifier_policy_id="retained-policy",
            )
            first.run()
            predecessor = first.session / "engine-state.json"

            library = root / "specialist-library"
            evidence_file = root / "specialist.out"

            def specialist(request, workdir):
                library.mkdir()
                evidence_file.write_text("specialist evidence\n")
                reactions = library / "reactions.py"
                dictionary = library / "dictionary.txt"
                reactions.write_text("# retained specialist kinetics\n")
                dictionary.write_text("fixture\n")
                return {
                    "schema_version": 1, "kind": "storca_specialist_rate_result",
                    "status": "verified", "request_sha256": request["request_sha256"],
                    "condition_contract_sha256": request["condition_contract_sha256"],
                    "rate_theory": "rrkm_master_equation",
                    "provider": {"name": "fixture", "version": "1"},
                    "applicability": {"status": "applicable", "diagnostics": ["fixture"]},
                    "uncertainty": {
                        "lower_rate_multiplier": 0.1, "upper_rate_multiplier": 10.0,
                        "source": "fixture", "source_kind": "method_calibration",
                        "applicability": "fixture dissociation",
                        "uncertainty_components": ["energy_transfer", "rate_theory"],
                    },
                    "pressure_dependence": {
                        "pressure_range_bar": [0.1, 10.0], "bath_gas": {"A": 1.0},
                        "energy_transfer_model": "fixture",
                    },
                    "evidence_artifacts": [{
                        "path": str(evidence_file),
                        "sha256": hashlib.sha256(evidence_file.read_bytes()).hexdigest(),
                    }],
                    "generated_kinetics_library": {
                        "route_id": route["route_id"],
                        "reaction_equation": route["reaction_equation"],
                        "replaces_reaction": route["reaction_equation"],
                        "verification_status": "verified_rrkm_master_equation",
                        "molecularity": 1,
                        "condition_rate": {"value": 2.0, "units": "s^-1",
                                           "temperature_K": 298.15, "pressure_bar": 1.0},
                        "collision_limit": {"applicable": False}, "library": str(library),
                        "checksums": {
                            "reactions.py": hashlib.sha256(reactions.read_bytes()).hexdigest(),
                            "dictionary.txt": hashlib.sha256(dictionary.read_bytes()).hexdigest(),
                        },
                    },
                }

            def forbidden_path(*_):
                raise AssertionError("Retained specialist continuation must not rerun ORCA")

            def rerun(libraries, workdir):
                repaired = _model(routes=[], loss=0.001, t95=None)
                evidence = repaired["rmg_evidence"]
                evidence["generated_kinetics_libraries"] = [{"route_id": route["route_id"]}]
                evidence["mechanism_inspection"]["reactions"] = [{
                    "reaction_index": 0, "equation": route["reaction_equation"],
                    "molecularity": 1, "evaluated_forward_rate_coefficient_si": 2.0,
                    "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
                    "source_library": None, "kinetics_comment": None,
                }]
                return repaired

            continue_verification_engine(
                initial, CONTRACT, root, predecessor,
                verify_route=forbidden_path, build_rate=lambda *_: {}, rerun_model=rerun,
                path_verifier_policy_id="retained-policy",
                build_specialist_rate=specialist,
                config=VerificationEngineConfig(max_iterations=1),
            )
            sessions = sorted(root.glob("verification-engine*/engine-state.json"))
            continued = json.loads(sessions[-1].read_text())

        self.assertEqual(path_calls, [route["route_id"]])
        self.assertTrue(continued["iterations"][0]["retained_path_reused"])
        self.assertFalse(continued["continuation_evidence"]["path_verifier_reexecuted"])

    def test_arkane_central_rate_without_uncertainty_cannot_converge_observable(self):
        initial = _model(weights=(1.0, 0.0))["rmg_evidence"]
        replacement_records = []

        def verified(route, workdir):
            return {"route_verification": {"status": "irc_endpoint_verified"}}

        def arkane_rate(route, path_result, workdir):
            library = Path(workdir) / "library"
            library.mkdir(parents=True)
            item = _replacement(route, library)
            item["source"] = "arkane_tst"
            item.update(rotor_rate_evidence_fixture(Path(workdir) / "statmech-fixture", route_id=route["route_id"]))
            replacement_records.append(item)
            return item

        def rerun(libraries, workdir):
            return _model(["init"], weights=(1.0, 0.0))

        with tempfile.TemporaryDirectory() as temp:
            engine = VerificationEngine(
                initial, CONTRACT, Path(temp), verify_route=verified,
                build_rate=arkane_rate, rerun_model=rerun,
                config=VerificationEngineConfig(max_iterations=1),
            )
            summary = engine.run()
            state = json.loads((engine.session / "engine-state.json").read_text())

        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("calculated_rate_uncertainty_missing_for_init", summary["reason"])
        self.assertTrue(
            state["verified_replacements"][0]["rate_uncertainty_propagation_required"]
        )

    def test_newly_discovered_rate_is_propagated_but_not_a_converged_lifetime(self):
        initial = _model(routes=[], loss=0.0, t95=None)["rmg_evidence"]
        candidate = {
            **ROUTES[0], "native_candidate_id": "init",
            "source": "unexpected_completed_irc_endpoint",
            "discovery_sources": ["unexpected_completed_irc_endpoint"],
            "reactant_smiles": ["C", "[OH]"], "product_smiles": ["[CH3]", "O"],
            "reactant_labels": ["stability", "A"], "product_labels": ["R", "P1"],
            "kinetic_relevance": "unknown", "screening_status": "candidate_unverified",
            "verification_queue_class": "unexpected_channel_unverified",
        }
        calls = []
        def verify(route, workdir):
            calls.append("path")
            return {"route_verification": {"status": "verified_transition_state_path"}}
        def rate(route, path_result, workdir):
            calls.append("rate")
            library = Path(workdir) / "library"
            library.mkdir(parents=True)
            result = _replacement(route, library)
            result["source"] = "arkane_tst"
            result.update(rotor_rate_evidence_fixture(Path(workdir) / "statmech-fixture", route_id=route["route_id"]))
            result["rate_uncertainty"] = {"status": "calibration_required"}
            return result
        def rerun(libraries, workdir):
            calls.append("propagation")
            return _model(["init"], weights=(1.0, 0.0))
        with tempfile.TemporaryDirectory() as temp:
            engine = VerificationEngine(
                initial, CONTRACT, Path(temp), verify_route=verify,
                build_rate=rate, rerun_model=rerun,
                supplemental_candidate_routes=[candidate],
                config=VerificationEngineConfig(max_iterations=1),
            )
            summary = engine.run()
            state = json.loads((engine.session / "engine-state.json").read_text())
        self.assertEqual(calls, ["path", "rate", "propagation"])
        self.assertEqual(state["iterations"][0]["native_candidate_handoff_status"],
                         "verified_rate_inserted_and_repropagated")
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("calculated_rate_uncertainty_missing_for_init", summary["reason"])
        self.assertIsNone(summary.get("estimated_time_to_retention_seconds"))

    @patch("storca.rmg_bridge_client.run_rmg_bridge")
    def test_arkane_rate_uncertainty_is_propagated_before_convergence(self, bridge):
        def with_sensitivity(wrapper):
            wrapper["rmg_evidence"]["artifacts"] = {
                "chemkin": "chem.inp", "species_dictionary": "dictionary.txt",
            }
            wrapper["rmg_evidence"]["independent_cantera_propagation"][
                "reaction_sensitivity_attribution"
            ] = {
                "status": "completed",
                "reactions": [{
                    "reaction_equation": ROUTES[0]["reaction_equation"],
                    "absolute_log_t95_sensitivity": 1.0,
                }],
            }
            return wrapper

        initial_wrapper = with_sensitivity(_model(weights=(1.0, 0.0)))

        def verified(route, workdir):
            return {"route_verification": {"status": "irc_endpoint_verified"}}

        def arkane_rate(route, path_result, workdir):
            library = Path(workdir) / "library"
            library.mkdir(parents=True)
            return {**_replacement(route, library), "source": "arkane_tst",
                    **rotor_rate_evidence_fixture(Path(workdir) / "statmech-fixture", route_id=route["route_id"])}

        def rerun(libraries, workdir):
            return with_sensitivity(_model(["init"], weights=(1.0, 0.0)))

        bridge.side_effect = [
            _bound_propagation(52.0), _bound_propagation(48.0),
        ]
        calibrated_uncertainty = self._uncertainty("init")
        calibrated_uncertainty["uncertainty_components"] = ["total_rate_model"]
        with tempfile.TemporaryDirectory() as temp:
            engine = VerificationEngine(
                initial_wrapper["rmg_evidence"], CONTRACT, Path(temp),
                verify_route=verified, build_rate=arkane_rate, rerun_model=rerun,
                rate_uncertainties=[calibrated_uncertainty],
                config=VerificationEngineConfig(max_iterations=1),
            )
            summary = engine.run()
            envelope = json.loads(
                (engine.session / "combined-prediction-rate-sensitivity-002.json").read_text()
            )

        self.assertEqual(summary["status"], "verified_t95_converged", summary)
        self.assertEqual(envelope["uncertain_route_ids"], ["init"])
        self.assertEqual(bridge.call_count, 2)

    def test_unhandled_interrupt_closes_engine_state(self):
        with tempfile.TemporaryDirectory() as temp:
            engine = VerificationEngine(
                {}, CONTRACT, Path(temp),
                verify_route=lambda *_: {}, build_rate=lambda *_: {}, rerun_model=lambda *_: {},
            )

            def interrupt():
                raise KeyboardInterrupt("test interrupt")

            engine._run_started = interrupt
            with self.assertRaises(KeyboardInterrupt):
                engine.run()
            state = json.loads((engine.session / "engine-state.json").read_text())
        self.assertEqual(state["status"], "failed_closed")
        self.assertEqual(state["terminal_error"]["error_type"], "KeyboardInterrupt")

    def test_continuation_retries_deferred_route_only_under_new_policy_and_reuses_workspace(self):
        initial = _model(weights=(1.0, 0.0))["rmg_evidence"]
        first_workspaces = []

        def deferred(route, workdir):
            first_workspaces.append(Path(workdir))
            Path(workdir).mkdir(parents=True)
            return {"route_verification": {"status": "stationary_point_validation_failed"}}

        def should_not_run(*args):
            raise AssertionError("A deferred route cannot create a rate")

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = VerificationEngine(
                initial, CONTRACT, root, verify_route=deferred,
                build_rate=should_not_run, rerun_model=should_not_run,
                path_verifier_policy_id="path-policy-v1",
            )
            first.run()
            predecessor = first.session / "engine-state.json"
            second_workspaces = []

            def verified(route, workdir):
                second_workspaces.append(Path(workdir))
                return {"route_verification": {"status": "irc_endpoint_verified"}}

            library = root / "verified-library"

            def rate(route, path_result, workdir):
                library.mkdir()
                item = _replacement(route, library)
                return {"rate_verification_status": "validated", "replacement": item}

            def rerun(library_paths, workdir):
                return _model(["init"], weights=(1.0, 0.0))

            summary = continue_verification_engine(
                initial, CONTRACT, root, predecessor,
                verify_route=verified, build_rate=rate, rerun_model=rerun,
                path_verifier_policy_id="path-policy-v2",
                retry_deferred_route_ids=["init"],
            )
            second_state = json.loads(
                Path(summary["artifacts"]["engine_state"]).read_text()
            )

        self.assertEqual(second_workspaces, first_workspaces)
        self.assertNotEqual(
            summary["artifacts"]["session_directory"], str(first.session),
        )
        self.assertEqual(
            second_state["continuation_evidence"]["path_verifier_policy_id"],
            "path-policy-v2",
        )
        self.assertEqual(
            second_state["route_policy_attempts"]["init"],
            ["path-policy-v1", "path-policy-v2"],
        )

    def test_continuation_rejects_condition_mutation_and_unchanged_policy_retry(self):
        initial = _model(weights=(1.0, 0.0))["rmg_evidence"]

        def deferred(route, workdir):
            Path(workdir).mkdir(parents=True)
            return {"route_verification": {"status": "surface_unresolved"}}

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = VerificationEngine(
                initial, CONTRACT, root, verify_route=deferred,
                build_rate=lambda *_: {}, rerun_model=lambda *_: {},
                path_verifier_policy_id="same-policy",
            )
            first.run()
            predecessor = first.session / "engine-state.json"
            callbacks = {
                "verify_route": lambda *_: {}, "build_rate": lambda *_: {},
                "rerun_model": lambda *_: {},
            }
            with self.assertRaisesRegex(ValueError, "does not exactly match"):
                continue_verification_engine(
                    initial, {**CONTRACT, "pressure_bar": 2.0}, root, predecessor,
                    path_verifier_policy_id="new-policy", **callbacks,
                )
            with self.assertRaisesRegex(ValueError, "already attempted"):
                continue_verification_engine(
                    initial, CONTRACT, root, predecessor,
                    path_verifier_policy_id="same-policy", **callbacks,
                )

    def test_stale_execution_recovery_reuses_policy_attempt_and_workspace(self):
        initial = _model(weights=(1.0, 0.0))["rmg_evidence"]
        first_workspaces = []

        def interrupted(route, workdir):
            first_workspaces.append(Path(workdir).parent)
            Path(workdir).mkdir(parents=True)
            return {"route_verification": {"status": "path_execution_failed"}}

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = VerificationEngine(
                initial, CONTRACT, root, verify_route=interrupted,
                build_rate=lambda *_: {}, rerun_model=lambda *_: {},
                path_verifier_policy_id="same-policy",
            )
            first.run()
            predecessor = first.session / "engine-state.json"
            state = json.loads(predecessor.read_text())
            state["status"] = "running"
            state["deferred_unverified_routes"] = {}
            state.pop("verification_summary", None)
            state["iterations"][-1]["status"] = "started"
            predecessor.write_text(json.dumps(state))

            recovered_workspaces = []

            def verified(route, workdir):
                recovered_workspaces.append(Path(workdir).parent)
                return {"route_verification": {"status": "irc_endpoint_verified"}}

            def rate(route, path_result, workdir):
                library = Path(workdir) / "library"
                library.mkdir(parents=True)
                return {
                    "rate_verification_status": "validated",
                    "replacement": _replacement(route, library),
                }

            summary = continue_verification_engine(
                initial, CONTRACT, root, predecessor,
                verify_route=verified, build_rate=rate,
                rerun_model=lambda *_: _model(["init"], weights=(1.0, 0.0)),
                path_verifier_policy_id="same-policy",
                recover_in_progress_route_ids=["init"],
                execution_recovery_evidence={"status": "stale"},
            )
            recovered_state = json.loads(Path(summary["artifacts"]["engine_state"]).read_text())

        self.assertEqual(recovered_workspaces, first_workspaces)
        self.assertEqual(recovered_state["route_policy_attempts"]["init"], ["same-policy"])
        self.assertTrue(recovered_state["iterations"][0]["same_policy_execution_recovery"])
        self.assertEqual(
            recovered_state["continuation_evidence"]["kind"],
            "verification_engine_stale_execution_recovery",
        )

    def _callbacks(self, root, *, weights=(1.0, 9.0), envelope_after=None,
                   bad_rate=False, bad_contract=False):
        calls = []
        replacements = []

        def verify_route(route, workdir):
            calls.append(route["route_id"])
            Path(workdir).mkdir(parents=True)
            return {"route_verification": {"status": "irc_endpoint_verified"}}

        def build_rate(route, path_result, workdir):
            library = Path(workdir) / f"library-{route['route_id']}"
            library.mkdir(parents=True)
            item = _replacement(route, library)
            replacements.append(item)
            return {"rate_verification_status": "validated", "replacement": item, "library": str(library)}

        def rerun_model(library_paths, workdir):
            result = _model(
                [item["route_id"] for item in replacements], weights=weights,
                applied_rate=4.0 if bad_rate else 2.0,
                envelope=envelope_after if len(replacements) == 1 else None,
            )
            if bad_contract:
                result["condition_contract"] = {**CONTRACT, "pressure_bar": 2.0}
            return result

        return calls, verify_route, build_rate, rerun_model

    def test_engine_verifies_generated_dependency_then_controlling_route(self):
        with tempfile.TemporaryDirectory() as temp:
            calls, verify, build, rerun = self._callbacks(Path(temp))
            summary = run_verification_engine(
                _model()["rmg_evidence"], CONTRACT, Path(temp),
                verify_route=verify, build_rate=build, rerun_model=rerun,
            )
        self.assertEqual(calls, ["init", "main"])
        self.assertEqual(summary["status"], "verified_t95_converged")
        self.assertEqual(summary["estimated_time_to_retention_seconds"], 50.0)
        self.assertEqual(summary["iterations_completed"], 2)
        self.assertEqual(summary["condition_contract"], CONTRACT)

    @patch("storca.verification_engine.run_combined_unresolved_sensitivity")
    def test_engine_automatically_propagates_complete_unresolved_rate_bounds(
        self, combined,
    ):
        calls = []
        combined.return_value = {
            "schema_version": 1,
            "kind": "combined_unresolved_route_sensitivity",
            "status": "completed",
            "combined_perturbation": True,
            "unresolved_route_ids": ["main"],
            "rate_bounds": [{
                "route_id": "main", "baseline_rate": 4.0,
                "lower_rate": 2.0, "upper_rate": 8.0,
                "units": "m^3/(mol*s)", "bound_source": "test calibration",
            }],
            "lower_bound_propagation": _bound_propagation(55.0),
            "upper_bound_propagation": _bound_propagation(48.0),
        }

        def add_allocation_evidence(wrapper):
            evidence = wrapper["rmg_evidence"] if "rmg_evidence" in wrapper else wrapper
            evidence["independent_cantera_propagation"]["reaction_sensitivity_attribution"] = {
                "status": "completed",
                "reactions": [
                    {"reaction_equation": route["reaction_equation"],
                     "absolute_log_t95_sensitivity": 1.0}
                    for route in ROUTES
                ],
            }
            for route in evidence["candidate_routes"]:
                route["retained_rate_evidence"] = {
                    "rate_coefficient_si": 4.0,
                    "rate_units": "m^3/(mol*s)",
                }
            return wrapper

        initial = add_allocation_evidence(_model()["rmg_evidence"])
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            library = root / "verified-init"

            def verify(route, workdir):
                calls.append(route["route_id"])
                return {"route_verification": {"status": "irc_endpoint_verified"}}

            def build(route, path_result, workdir):
                library.mkdir()
                return {
                    "rate_verification_status": "validated",
                    "replacement": _replacement(route, library),
                }

            def rerun(library_paths, workdir):
                return add_allocation_evidence(_model(["init"]))

            summary = run_verification_engine(
                initial, CONTRACT, root,
                verify_route=verify, build_rate=build, rerun_model=rerun,
                rate_uncertainties=[self._uncertainty("init"), self._uncertainty("main")],
            )

        self.assertEqual(calls, ["init"])
        self.assertEqual(summary["status"], "verified_t95_converged")
        self.assertEqual(summary["t95_seconds"], 50.0)
        combined.assert_called_once()

    def test_verified_multistep_bundle_inserts_all_rates_before_one_repropagation(self):
        rerun_calls = []

        def verify_multistep(route, workdir):
            return {"route_verification": {
                "status": "verified_multistep_path",
                "path_classification": "verified_multistep_path",
            }}

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            def build_multistep(route, path_result, workdir):
                segments = []
                for segment_route in ROUTES:
                    library = root / f"multistep-{segment_route['route_id']}"
                    library.mkdir()
                    rate = {
                        "rate_verification_status": "validated",
                        "replacement": _replacement(segment_route, library),
                        "library": str(library),
                    }
                    segments.append({"route": segment_route, "rate_result": rate})
                return {
                    "kind": "storca_verified_multistep_rate_bundle",
                    "status": "completed_multistep",
                    "segments": segments,
                }

            def rerun(library_paths, workdir):
                rerun_calls.append([Path(item) for item in library_paths])
                return _model(["init", "main"])

            summary = run_verification_engine(
                _model()["rmg_evidence"], CONTRACT, root,
                verify_route=verify_multistep,
                build_rate=build_multistep,
                rerun_model=rerun,
            )
            state = json.loads(Path(summary["artifacts"]["engine_state"]).read_text())

        self.assertEqual(summary["status"], "verified_t95_converged")
        self.assertEqual(len(rerun_calls), 1)
        self.assertEqual(len(rerun_calls[0]), 2)
        self.assertEqual(
            state["iterations"][0]["replacement_route_ids"], ["init", "main"],
        )
        self.assertEqual(
            state["iterations"][0]["parent_multistep_route_resolved"], "init",
        )

    def test_failed_multistep_repropagation_retains_precise_network_pending_state(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            def build_multistep(route, path_result, workdir):
                segments = []
                for segment_route in ROUTES:
                    library = root / f"pending-{segment_route['route_id']}"
                    library.mkdir()
                    segments.append({
                        "route": segment_route,
                        "rate_result": {
                            "rate_verification_status": "validated",
                            "replacement": _replacement(segment_route, library),
                            "library": str(library),
                        },
                    })
                return {
                    "kind": "storca_verified_multistep_rate_bundle",
                    "status": "completed_multistep",
                    "segments": segments,
                }

            summary = run_verification_engine(
                _model()["rmg_evidence"], CONTRACT, root,
                verify_route=lambda *_: {"route_verification": {
                    "status": "verified_multistep_path",
                }},
                build_rate=build_multistep,
                rerun_model=lambda *_: (_ for _ in ()).throw(
                    RuntimeError("test repaired network rejected intermediate")
                ),
            )
            state = json.loads(Path(summary["artifacts"]["engine_state"]).read_text())

        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("verified_multistep_path_network_insertion_pending", summary["reason"])
        self.assertEqual(
            state["iterations"][0]["scientific_failure_state"],
            "verified_multistep_path_network_insertion_pending",
        )
        self.assertIsNone(summary["t95_seconds"])

    def test_engine_defers_unverified_route_and_continues_ranked_routes(self):
        independent = [
            dict(ROUTES[0], co_reactants=[], initiation_status="directly_present"),
            dict(ROUTES[1], co_reactants=[], initiation_status="directly_present"),
        ]
        initial = _model(weights=(9.0, 1.0), routes=independent)["rmg_evidence"]
        initial["candidate_routes"] = independent
        calls = []

        def unresolved(route, workdir):
            calls.append(route["route_id"])
            return {"route_verification": {
                "status": "surface_unresolved", "path_classification": "surface_unresolved",
            }}

        def should_not_run(*args):
            raise AssertionError("No rate or rerun is allowed for an unresolved path")

        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=unresolved,
                build_rate=should_not_run, rerun_model=should_not_run,
                config=VerificationEngineConfig(max_iterations=3),
            )
            state = json.loads(Path(summary["artifacts"]["engine_state"]).read_text())
        self.assertEqual(calls, ["init", "main"])
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("all_unverified_routes_deferred", summary["reason"])
        self.assertEqual(len(summary["deferred_unverified_routes"]), 2)
        self.assertEqual([item["status"] for item in state["iterations"]], [
            "deferred_unverified", "deferred_unverified",
        ])

    def test_engine_defers_route_without_common_adiabatic_spin_surface(self):
        initial = _model(weights=(10.0, 0.0), routes=[ROUTES[0]])["rmg_evidence"]
        initial["candidate_routes"] = [ROUTES[0]]

        def incompatible_spin_surface(route, workdir):
            return {"route_verification": {
                "status": "no_common_adiabatic_spin_surface",
                "path_classification": "surface_unresolved",
            }}

        def should_not_run(*args):
            raise AssertionError("No rate or rerun is allowed across incompatible spin surfaces")

        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=incompatible_spin_surface,
                build_rate=should_not_run, rerun_model=should_not_run,
            )

        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("all_unverified_routes_deferred", summary["reason"])
        self.assertEqual(
            summary["deferred_unverified_routes"][0]["status"],
            "no_common_adiabatic_spin_surface",
        )

    def test_engine_defers_bounded_mapping_search_exhaustion(self):
        initial = _model(weights=(10.0, 0.0), routes=[ROUTES[0]])["rmg_evidence"]
        initial["candidate_routes"] = [ROUTES[0]]

        def exhausted_mapping_search(route, workdir):
            return {"route_verification": {
                "status": "mapping_search_exhausted",
                "path_classification": "surface_unresolved",
                "mapping": {"status": "mapping_search_exhausted", "search_states": 200000},
            }}

        def should_not_run(*args):
            raise AssertionError("No rate or rerun is allowed without a resolved atom mapping")

        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=exhausted_mapping_search,
                build_rate=should_not_run, rerun_model=should_not_run,
            )

        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("all_unverified_routes_deferred", summary["reason"])
        self.assertEqual(
            summary["deferred_unverified_routes"][0]["status"],
            "mapping_search_exhausted",
        )

    def test_engine_backtracks_to_alternative_generated_intermediate_source(self):
        routes = [
            {
                "route_id": "strong-source",
                "reaction_equation": "stability(1) + A(2) => R(3) + P0(4)",
                "initiation_status": "directly_present",
                "co_reactants": [{"label": "A(2)", "origin": "initial_scenario_species"}],
                "resolved_endpoints": {
                    "reactants": [{"label": "stability(1)"}, {"label": "A(2)"}],
                    "products": [{"label": "R(3)"}, {"label": "P0(4)"}],
                },
            },
            {
                "route_id": "alternate-source",
                "reaction_equation": "stability(1) + B(6) => R(3) + P1(7)",
                "initiation_status": "directly_present",
                "co_reactants": [{"label": "B(6)", "origin": "initial_scenario_species"}],
                "resolved_endpoints": {
                    "reactants": [{"label": "stability(1)"}, {"label": "B(6)"}],
                    "products": [{"label": "R(3)"}, {"label": "P1(7)"}],
                },
            },
            {
                "route_id": "main",
                "reaction_equation": "stability(1) + R(3) => P2(5)",
                "initiation_status": "thermally_reachable",
                "co_reactants": [{"label": "R(3)", "origin": "generated_intermediate"}],
                "resolved_endpoints": {
                    "reactants": [{"label": "stability(1)"}, {"label": "R(3)"}],
                    "products": [{"label": "P2(5)"}],
                },
            },
        ]
        weights = (1.0, 0.5, 8.5)

        def model(applied=()):
            flux_reactions = [
                {
                    "reaction_index": index,
                    "reaction_equation": route["reaction_equation"],
                    "integrated_target_destruction_kmol": weight,
                    "integrated_forward_extent_kmol": weight,
                    "integrated_reverse_extent_kmol": 0.0,
                }
                for index, (route, weight) in enumerate(zip(routes, weights))
            ]
            evidence = {
                "status": "completed",
                "execution_assessment": {"status": "completed"},
                "time_coverage": {"complete": True},
                "candidate_routes": routes,
                "independent_cantera_propagation": {
                    "status": "completed",
                    "reactor_model": isothermal_model(),
                    "coverage": {"complete": True, "requested_seconds": 1000.0, "simulated_seconds": 1000.0},
                    "estimated_time_to_retention_seconds": 50.0,
                    "target_loss_fraction": 0.20,
                    "target_profile": [
                        {"time_seconds": 0.0, "target_fraction_remaining": 1.0},
                        {"time_seconds": 1000.0, "target_fraction_remaining": 0.8},
                    ],
                    "reaction_flux_attribution": {
                        "status": "completed",
                        "total_integrated_target_destruction_kmol": sum(weights),
                        "numerical_closure_relative_error": 0.01,
                        "reactions": flux_reactions,
                    },
                },
                "kinetics_validation": {"status": "passed", "violation_count": 0},
                "generated_kinetics_libraries": [{"route_id": route_id} for route_id in applied],
                "mechanism_inspection": {
                    "status": "completed",
                    "reactions": [
                        {
                            "reaction_index": index,
                            "equation": next(
                                route["reaction_equation"] for route in routes
                                if route["route_id"] == route_id
                            ),
                            "molecularity": 2,
                            "evaluated_forward_rate_coefficient_si": 2.0,
                            "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
                            "source_library": None,
                            "kinetics_comment": None,
                        }
                        for index, route_id in enumerate(applied)
                    ],
                },
            }
            return {"condition_contract": CONTRACT, "rmg_evidence": evidence}

        calls = []
        replacements = []

        def verify(route, workdir):
            calls.append(route["route_id"])
            if route["route_id"] == "strong-source":
                return {"route_verification": {"status": "surface_unresolved"}}
            return {"route_verification": {"status": "irc_endpoint_verified"}}

        def build(route, path_result, workdir):
            library = Path(workdir) / f"library-{route['route_id']}"
            library.mkdir(parents=True)
            item = _replacement(route, library)
            replacements.append(item)
            return {"rate_verification_status": "validated", "replacement": item}

        def rerun(library_paths, workdir):
            return model([item["route_id"] for item in replacements])

        contract = {
            **CONTRACT,
            "composition": {"stability": 0.01, "A": 0.49, "B": 0.50},
        }
        # The rerun must retain the exact contract supplied to this session.
        def rerun_with_contract(library_paths, workdir):
            result = rerun(library_paths, workdir)
            result["condition_contract"] = contract
            return result

        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                model()["rmg_evidence"], contract, Path(temp),
                verify_route=verify, build_rate=build, rerun_model=rerun_with_contract,
                config=VerificationEngineConfig(max_iterations=4),
            )
        self.assertEqual(calls, ["strong-source", "alternate-source", "main"])
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("all_unverified_routes_deferred", summary["reason"])

    def test_flux_coverage_alone_does_not_stop_without_robustness(self):
        independent = [dict(ROUTES[0], co_reactants=[], initiation_status="directly_present"),
                       dict(ROUTES[1], co_reactants=[], initiation_status="directly_present")]
        initial = _model(weights=(9.6, 0.4), routes=independent)["rmg_evidence"]
        initial["candidate_routes"] = independent
        with tempfile.TemporaryDirectory() as temp:
            calls, verify, build, rerun_base = self._callbacks(Path(temp), weights=(9.6, 0.4))
            def rerun(paths, workdir):
                result = rerun_base(paths, workdir)
                result["rmg_evidence"]["candidate_routes"] = independent
                return result
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=verify, build_rate=build, rerun_model=rerun,
            )
        self.assertEqual(calls, ["init", "main"])
        self.assertEqual(summary["status"], "verified_t95_converged")

    def test_downstream_collision_warning_does_not_block_upstream_reranking(self):
        with tempfile.TemporaryDirectory() as temp:
            calls, verify, build, rerun_base = self._callbacks(Path(temp))
            def rerun(paths, workdir):
                result = rerun_base(paths, workdir)
                if len(paths) == 1:
                    result["rmg_evidence"]["kinetics_validation"] = {
                        "status": "kinetics_unreliable", "violation_count": 1,
                    }
                return result
            summary = run_verification_engine(
                _model()["rmg_evidence"], CONTRACT, Path(temp),
                verify_route=verify, build_rate=build, rerun_model=rerun,
            )
        self.assertEqual(calls, ["init", "main"])
        self.assertEqual(summary["status"], "verified_t95_converged")

    def test_combined_physical_bounds_can_stop_after_controlling_coverage(self):
        independent = [dict(ROUTES[0], co_reactants=[], initiation_status="directly_present"),
                       dict(ROUTES[1], co_reactants=[], initiation_status="directly_present")]
        envelope = {
            "status": "completed", "combined_perturbation": True,
            "unresolved_route_ids": ["main"],
            "rate_bounds": [{"route_id": "main", "lower_rate": 0.0, "upper_rate": 3.0,
                             "bound_source": "retained_collision_ceiling"}],
            "lower_bound_propagation": _bound_propagation(54.0),
            "upper_bound_propagation": _bound_propagation(47.0),
        }
        initial = _model(weights=(9.6, 0.4), routes=independent)["rmg_evidence"]
        initial["candidate_routes"] = independent
        with tempfile.TemporaryDirectory() as temp:
            calls, verify, build, rerun_base = self._callbacks(
                Path(temp), weights=(9.6, 0.4), envelope_after=envelope,
            )
            def rerun(paths, workdir):
                result = rerun_base(paths, workdir)
                result["rmg_evidence"]["candidate_routes"] = independent
                return result
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=verify, build_rate=build, rerun_model=rerun,
            )
        self.assertEqual(calls, ["init"])
        self.assertEqual(summary["status"], "verified_t95_converged")
        self.assertAlmostEqual(summary["verified_flux_coverage"], 0.96)

    def test_completed_no_loss_flux_stays_model_scoped_without_verified_rates(self):
        initial = _model(routes=[], weights=(0.0, 0.0), t95=None, loss=0.0)["rmg_evidence"]
        def should_not_run(*args):
            raise AssertionError("No chemistry callback should run")
        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=should_not_run,
                build_rate=should_not_run, rerun_model=should_not_run,
            )
        self.assertEqual(summary["status"], "no_target_loss_in_completed_rmg_model")
        self.assertEqual(summary["target_loss_fraction"], 0.0)
        self.assertEqual(summary["unresolved_loss_upper_bound"], 0.0)
        self.assertIsNone(summary["verified_flux_coverage"])
        self.assertEqual(summary["route_verification_status"], "not_required_for_retained_model")
        self.assertEqual(summary["rate_verification_status"], "not_required_for_retained_model")
        self.assertEqual(summary["propagation_status"], "completed_model")

    def test_below_threshold_completed_model_skips_tiny_flux_orca_work(self):
        initial = _model(t95=None, loss=0.01, weights=(0.001, 0.002))["rmg_evidence"]
        def should_not_run(*args):
            raise AssertionError("Below-threshold passed model should not start ORCA")
        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=should_not_run,
                build_rate=should_not_run, rerun_model=should_not_run,
            )
        self.assertEqual(summary["status"], "no_target_loss_in_completed_rmg_model")
        self.assertIsNone(summary["estimated_time_to_retention_seconds"])
        self.assertTrue(summary["screening_estimate_prohibited"])

    def test_below_threshold_model_with_invalid_collision_rate_does_not_skip_verification(self):
        initial = _model(t95=None, loss=0.01)["rmg_evidence"]
        initial["kinetics_validation"] = {"status": "kinetics_unreliable", "violation_count": 1}
        with tempfile.TemporaryDirectory() as temp:
            calls, verify, build, rerun = self._callbacks(Path(temp))
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=verify, build_rate=build, rerun_model=rerun,
            )
        self.assertTrue(calls)
        self.assertEqual(summary["status"], "verified_t95_converged")

    @patch("storca.verification_engine.run_collision_limit_rate_envelope")
    def test_collision_envelope_can_certify_only_the_bounded_observable(self, envelope):
        initial = _model(routes=[], weights=(0.0, 0.0), t95=None, loss=0.01)["rmg_evidence"]
        initial["kinetics_validation"] = {
            "status": "kinetics_unreliable",
            "violation_count": 1,
            "condition_mechanism_validation": {"violators": [{"reaction_index": 4}]},
        }
        envelope.return_value = {
            "schema_version": 1,
            "kind": "collision_limit_rate_envelope",
            "status": "completed",
            "observable_classification_robust": True,
            "maximum_loss_upper_bound_fraction": 0.02,
            "loss_threshold_fraction": 0.05,
            "rate_bounds": [{"reaction_index": 4, "lower_multiplier": 0.0,
                             "upper_multiplier": 0.1}],
        }

        def should_not_run(*args):
            raise AssertionError("Collision-bounded no-loss model should not start ORCA")

        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=should_not_run,
                build_rate=should_not_run, rerun_model=should_not_run,
                rmg_env="rmg4_env",
            )
            session = Path(summary["artifacts"]["session_directory"])
            self.assertTrue((session / "collision-limit-rate-envelope-001.json").is_file())

        self.assertEqual(summary["status"], "verified_below_loss_threshold")
        self.assertEqual(
            summary["rate_verification_status"],
            "collision_bounded_for_observable_not_rate_verified",
        )
        self.assertEqual(
            summary["claim_scope"],
            "condition_specific_retained_network_collision_bounded_observable",
        )
        self.assertAlmostEqual(summary["unresolved_loss_upper_bound"], 0.01)
        envelope.assert_called_once()

    @patch("storca.verification_engine.run_collision_limit_rate_envelope")
    def test_collision_envelope_crossing_threshold_fails_closed(self, envelope):
        initial = _model(routes=[], weights=(0.0, 0.0), t95=None, loss=0.01)["rmg_evidence"]
        initial["kinetics_validation"] = {
            "status": "kinetics_unreliable",
            "violation_count": 1,
            "condition_mechanism_validation": {"violators": [{"reaction_index": 4}]},
        }
        envelope.return_value = {
            "schema_version": 1,
            "kind": "collision_limit_rate_envelope",
            "status": "completed",
            "observable_classification_robust": False,
            "maximum_loss_upper_bound_fraction": 0.08,
            "loss_threshold_fraction": 0.05,
        }

        def should_not_run(*args):
            raise AssertionError("No route exists to verify in this fixture")

        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=should_not_run,
                build_rate=should_not_run, rerun_model=should_not_run,
                config=VerificationEngineConfig(max_iterations=1),
            )

        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("robustness_not_established", summary["reason"])

    def test_applied_rate_mismatch_fails_closed_and_retains_iteration(self):
        with tempfile.TemporaryDirectory() as temp:
            calls, verify, build, rerun = self._callbacks(Path(temp), bad_rate=True)
            summary = run_verification_engine(
                _model()["rmg_evidence"], CONTRACT, Path(temp),
                verify_route=verify, build_rate=build, rerun_model=rerun,
            )
            session = Path(summary["artifacts"]["session_directory"])
            state = json.loads((session / "engine-state.json").read_text())
            iteration = next(session.glob("iteration-001-*"))
            self.assertTrue((iteration / "rerun-result.json").is_file())
            self.assertTrue((iteration / "repaired-model-validation.json").is_file())
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIsNone(summary["estimated_time_to_retention_seconds"])
        self.assertEqual(state["status"], "failed_closed")
        self.assertIn("applied_replacement_not_rate_matched", summary["reason"])

    def test_rerun_must_return_exact_condition_contract(self):
        with tempfile.TemporaryDirectory() as temp:
            calls, verify, build, rerun = self._callbacks(Path(temp), bad_contract=True)
            summary = run_verification_engine(
                _model()["rmg_evidence"], CONTRACT, Path(temp),
                verify_route=verify, build_rate=build, rerun_model=rerun,
            )
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("exact immutable condition contract", summary["reason"])

    def test_max_iterations_is_checked_after_final_rerun(self):
        with tempfile.TemporaryDirectory() as temp:
            calls, verify, build, rerun = self._callbacks(Path(temp))
            summary = run_verification_engine(
                _model()["rmg_evidence"], CONTRACT, Path(temp),
                verify_route=verify, build_rate=build, rerun_model=rerun,
                config=VerificationEngineConfig(max_iterations=1),
            )
        self.assertEqual(calls, ["init"])
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("maximum_verification_iterations_reached", summary["reason"])

    def test_engine_can_verify_non_target_network_prerequisite(self):
        main = {
            "route_id": "target-main", "reaction_equation": "stability(1) + R(3) => P2(5)",
            "initiation_status": "thermally_reachable",
            "co_reactants": [{"label": "R(3)", "origin": "generated_intermediate"}],
            "resolved_endpoints": {
                "reactants": [{"label": "stability(1)"}, {"label": "R(3)"}],
                "products": [{"label": "P2(5)"}],
            },
        }
        source = {
            "route_id": "network-source", "reaction_equation": "A(2) => R(3)",
            "resolved_endpoints": {
                "reactants": [{"label": "A(2)"}], "products": [{"label": "R(3)"}],
            },
        }
        initial = _model()["rmg_evidence"]
        initial["candidate_routes"] = [main]
        initial["network_routes"] = [source]
        flux = initial["independent_cantera_propagation"]["reaction_flux_attribution"]
        flux.update(total_integrated_target_destruction_kmol=10.0)
        flux["reactions"] = [
            {"reaction_index": 0, "reaction_equation": main["reaction_equation"],
             "integrated_target_destruction_kmol": 10.0,
             "integrated_forward_extent_kmol": 10.0, "integrated_reverse_extent_kmol": 0.0},
            {"reaction_index": 1, "reaction_equation": source["reaction_equation"],
             "integrated_target_destruction_kmol": 0.0,
             "integrated_forward_extent_kmol": 1.0, "integrated_reverse_extent_kmol": 0.0},
        ]
        calls = []
        replacement = None

        def verify(route, workdir):
            calls.append(route["route_id"])
            return {"status": "irc_endpoint_verified"}

        def build(route, path_result, workdir):
            nonlocal replacement
            library = Path(workdir) / "network-library"
            library.mkdir(parents=True)
            replacement = {
                "route_id": route["route_id"], "reaction_equation": route["reaction_equation"],
                "verification_status": "verified_network_rate", "molecularity": 1,
                "condition_rate": {"value": 2.0, "units": "s^-1", "temperature_K": 298.15,
                                   "pressure_bar": 1.0},
            }
            return {"replacement": replacement, "library": str(library), "rate_verification_status": "validated"}

        def rerun(paths, workdir):
            evidence = dict(initial)
            evidence["generated_kinetics_libraries"] = [{"route_id": "network-source"}]
            evidence["mechanism_inspection"] = {"reactions": [{
                "reaction_index": 1, "equation": source["reaction_equation"], "molecularity": 1,
                "evaluated_forward_rate_coefficient_si": 2.0,
                "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
                "source_library": None, "kinetics_comment": None,
            }]}
            return {"condition_contract": CONTRACT, "rmg_evidence": evidence}

        with tempfile.TemporaryDirectory() as temp:
            summary = run_verification_engine(
                initial, CONTRACT, Path(temp), verify_route=verify, build_rate=build, rerun_model=rerun,
                config=VerificationEngineConfig(max_iterations=1),
            )
        self.assertEqual(calls, ["network-source"])
        self.assertEqual(summary["status"], "verification_incomplete")
        self.assertIn("maximum_verification_iterations_reached", summary["reason"])


if __name__ == "__main__":
    unittest.main()
