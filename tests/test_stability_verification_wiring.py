import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

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
        self.assertEqual(result["verification_summary"]["status"], "no_target_loss_in_completed_rmg_model")
        self.assertEqual(result["assessment"]["status"], "no_loss_detected_in_rmg_model")
        self.assertIsNone(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"])

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
        self.assertEqual(result["verification_summary"]["status"], "verified_t95_converged")
        self.assertEqual(result["verification_summary"]["condition_contract"], result["condition_contract"])
        self.assertEqual(result["assessment"]["status"], "orca_verified_condition_dependent_instability")
        self.assertEqual(result["kinetic_lifetime"]["estimated_time_to_retention_seconds"], 50.0)

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
        evidence = {
            "candidate_routes": [{
                "route_id": "current-id",
                "reaction_equation": "stability(8)+A(2)<=>P(4)",
            }],
            "network_routes": [{
                "route_id": "reverse-network-copy",
                "reaction_equation": "P(9)<=>A(7)+stability(3)",
            }],
        }
        replacement = {
            "route_id": "old-id",
            "reaction_equation": "A(1)+stability(1)<=>P(2)",
        }
        rebound = _planning_replacements([replacement], evidence)
        self.assertEqual(rebound[0]["route_id"], "current-id")

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
