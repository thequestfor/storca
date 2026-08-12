import json
import tempfile
import unittest
from pathlib import Path

from storca.verification_engine import VerificationEngine, VerificationEngineConfig, run_verification_engine


CONTRACT = {
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


def _replacement(route, library):
    return {
        "route_id": route["route_id"],
        "reaction_equation": route["reaction_equation"],
        "verification_status": "verified_test_rate",
        "molecularity": 2,
        "condition_rate": {
            "value": 2.0, "units": "m^3/(mol*s)",
            "temperature_K": 298.15, "pressure_bar": 1.0,
        },
        "collision_limit": {
            "value": 3.0, "units": "m^3/(mol*s)", "source": "retained_test_ceiling",
        },
        "library": str(library),
    }


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
