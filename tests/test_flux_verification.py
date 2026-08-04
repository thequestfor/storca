import unittest
from unittest.mock import patch

from storca.flux_verification import (
    assess_repaired_mechanism,
    assess_target_loss_commitment,
    assess_t95_robustness,
    build_flux_ranked_verification_plan,
    reaction_signature,
    run_combined_unresolved_sensitivity,
    verify_applied_replacement,
)
from storca.generated_kinetics import validate_rate_replacement_evidence


CONDITIONS = {
    "temperature_K": 298.15,
    "pressure_bar": 1.0,
    "retention_fraction": 0.95,
}


def replacement(route_id, *, rate=1.0e9, collision=2.0e9):
    return {
        "route_id": route_id,
        "reaction_equation": f"stability + X <=> P-{route_id}",
        "verification_status": "verified_orca_irc_arkane",
        "molecularity": 2,
        "condition_rate": {
            "value": rate, "units": "cm3 molecule-1 s-1",
            "temperature_K": 298.15, "pressure_bar": 1.0,
        },
        "collision_limit": {
            "value": collision, "units": "cm3 molecule-1 s-1", "source": "retained_collision_model",
        },
    }


def propagation(*, t95=100.0, final=0.80, flux=True):
    result = {
        "status": "completed",
        "coverage": {"complete": True, "requested_seconds": 1000.0, "simulated_seconds": 1000.0},
        "estimated_time_to_retention_seconds": t95,
        "target_profile": [
            {"time_seconds": 0.0, "target_fraction_remaining": 1.0},
            {"time_seconds": 1000.0, "target_fraction_remaining": final},
        ],
    }
    if flux:
        result["reaction_flux_attribution"] = {
            "status": "completed",
            "total_integrated_target_destruction_kmol": 10.0,
            "numerical_closure_relative_error": 0.01,
            "reactions": [
                {"reaction_index": 0, "reaction_equation": "stability(1) + A(2) <=> P1(3)",
                 "integrated_target_destruction_kmol": 8.0},
                {"reaction_index": 1, "reaction_equation": "B(4) + stability(1) <=> P2(5)",
                 "integrated_target_destruction_kmol": 2.0},
            ],
        }
    return result


def rmg_evidence(*, flux=True):
    return {
        "status": "completed",
        "execution_assessment": {"status": "completed"},
        "time_coverage": {"complete": True},
        "candidate_routes": [
            {"route_id": "r1", "reaction_equation": "A(2) + stability(1)<=>P1(3)",
             "initiation_status": "directly_present"},
            {"route_id": "r2", "reaction_equation": "stability(1)+B(4)<=>P2(5)",
             "initiation_status": "thermally_reachable"},
        ],
        "independent_cantera_propagation": propagation(flux=flux),
        "kinetics_validation": {"status": "passed", "violation_count": 0},
    }


class FluxVerificationTests(unittest.TestCase):
    def test_partial_rmg_enumeration_queues_only_direct_atom_balanced_route(self):
        evidence = {
            "status": "incomplete",
            "enumeration_status": "partial_enumeration",
            "candidate_routes": [
                {
                    "route_id": "direct",
                    "reaction_equation": "stability(1) + oxygen(2) <=> P(3)",
                    "target_stoichiometry": 1.0,
                    "initiation_status": "directly_present",
                    "co_reactants": [{"origin": "initial_scenario_species"}],
                    "resolved_endpoints": {"reactants": [{"label": "stability"}], "products": [{"label": "P"}]},
                    "atom_balance": {"status": "passed"},
                },
                {
                    "route_id": "generated-radical",
                    "reaction_equation": "stability(1) + R(2) <=> P(3)",
                    "target_stoichiometry": 1.0,
                    "initiation_status": "requires_generated_intermediate",
                    "co_reactants": [{"origin": "generated_intermediate"}],
                    "resolved_endpoints": {"reactants": [{"label": "stability"}], "products": [{"label": "P"}]},
                    "atom_balance": {"status": "passed"},
                },
            ],
        }
        plan = build_flux_ranked_verification_plan(evidence, CONDITIONS)
        self.assertEqual(plan["status"], "verification_required")
        self.assertTrue(plan["partial_enumeration"])
        self.assertEqual(plan["selected_route_id"], "direct")
        self.assertEqual([item["route_id"] for item in plan["ranked_routes"]], ["direct"])

    def test_reaction_signature_is_orientation_and_order_independent(self):
        self.assertEqual(
            reaction_signature("2 A(1) + stability(2) <=> P(3)"),
            reaction_signature("P(3) <=> stability(2) + 2 A(1)"),
        )

    def test_alternate_parent_return_is_credited_to_radical_loss_route(self):
        flux = {
            "integrated_net_target_loss_kmol": 1.0,
            "reactions": [
                {
                    "reaction_index": 0,
                    "reaction_equation": "stability(1) + A(2) <=> R(3) + B(4)",
                    "integrated_net_extent_kmol": 100.0,
                    "integrated_net_target_loss_kmol": 100.0,
                },
                {
                    "reaction_index": 1,
                    "reaction_equation": "R(3) + B(4) <=> stability(1) + A(2)",
                    "integrated_net_extent_kmol": 99.0,
                    "integrated_net_target_loss_kmol": -99.0,
                },
            ],
        }
        result = assess_target_loss_commitment(flux, "stability")
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["committed_target_loss_kmol"], 1.0)
        self.assertEqual(result["reaction_commitment"][0]["allocated_parent_return_credit_kmol"], 99.0)
        self.assertEqual(result["reaction_commitment"][0]["commitment_fraction"], 0.01)

    def test_complete_radical_cycle_has_no_committed_target_loss(self):
        flux = {
            "integrated_net_target_loss_kmol": 0.0,
            "reactions": [
                {
                    "reaction_index": 0,
                    "reaction_equation": "stability(1) + A(2) <=> R(3) + B(4)",
                    "integrated_net_extent_kmol": 100.0,
                    "integrated_net_target_loss_kmol": 100.0,
                },
                {
                    "reaction_index": 1,
                    "reaction_equation": "R(3) + B(4) <=> stability(1) + A(2)",
                    "integrated_net_extent_kmol": 100.0,
                    "integrated_net_target_loss_kmol": -100.0,
                },
            ],
        }
        result = assess_target_loss_commitment(flux, "stability")
        self.assertEqual(result["status"], "no_committed_target_loss")
        self.assertEqual(result["committed_target_loss_kmol"], 0.0)

    def test_complete_radical_cycle_is_not_scheduled_for_orca_verification(self):
        evidence = rmg_evidence()
        flux = evidence["independent_cantera_propagation"]["reaction_flux_attribution"]
        flux["integrated_net_target_loss_kmol"] = 0.0
        flux["reactions"] = [
            {
                "reaction_index": 0,
                "reaction_equation": "stability(1) + A(2) <=> P1(3)",
                "integrated_target_destruction_kmol": 100.0,
                "integrated_net_extent_kmol": 100.0,
                "integrated_net_target_loss_kmol": 100.0,
            },
            {
                "reaction_index": 1,
                "reaction_equation": "P1(3) <=> stability(1) + A(2)",
                "integrated_target_destruction_kmol": 0.0,
                "integrated_net_extent_kmol": 100.0,
                "integrated_net_target_loss_kmol": -100.0,
            },
        ]
        plan = build_flux_ranked_verification_plan(evidence, CONDITIONS)
        self.assertEqual(plan["status"], "no_target_destruction_flux")
        self.assertEqual(
            plan["target_loss_commitment"]["status"], "no_committed_target_loss",
        )
        self.assertIsNone(plan["selected_route_id"])

    def test_unrelated_parent_reformation_cannot_create_false_reassurance(self):
        flux = {
            "integrated_net_target_loss_kmol": 1.0,
            "reactions": [
                {
                    "reaction_index": 0,
                    "reaction_equation": "stability(1) + A(2) <=> R(3) + B(4)",
                    "integrated_net_extent_kmol": 100.0,
                    "integrated_net_target_loss_kmol": 100.0,
                },
                {
                    "reaction_index": 1,
                    "reaction_equation": "C(5) + D(6) <=> stability(1) + E(7)",
                    "integrated_net_extent_kmol": 99.0,
                    "integrated_net_target_loss_kmol": -99.0,
                },
            ],
        }
        result = assess_target_loss_commitment(flux, "stability")
        self.assertEqual(result["status"], "unresolved_parent_return_attribution")
        self.assertEqual(result["unallocated_parent_reformation_kmol"], 99.0)

    def test_plan_fails_closed_without_quantitative_flux(self):
        plan = build_flux_ranked_verification_plan(rmg_evidence(flux=False), CONDITIONS)
        self.assertEqual(plan["status"], "insufficient_solver_evidence")
        self.assertIn("quantitative_reaction_flux_missing_or_invalid", plan["blocking_reasons"])
        self.assertIsNone(plan["selected_route_id"])

    def test_plan_selects_only_the_largest_controlling_route_then_reranks(self):
        plan = build_flux_ranked_verification_plan(rmg_evidence(), CONDITIONS)
        self.assertEqual(plan["status"], "verification_required")
        self.assertEqual(plan["selected_route_id"], "r1")
        self.assertAlmostEqual(plan["ranked_routes"][0]["gross_target_destruction_fraction"], 0.8)
        self.assertEqual(plan["deferred_route_ids"], ["r2"])
        self.assertTrue(plan["requires_repropagation_after_each_verification"])

    def test_generated_reactant_prerequisite_precedes_high_flux_downstream_route(self):
        evidence = rmg_evidence()
        evidence["candidate_routes"][0].update({
            "resolved_endpoints": {
                "reactants": [{"label": "stability(1)"}, {"label": "oxygen(2)"}],
                "products": [{"label": "P1(3)"}, {"label": "R(6)"}],
            },
            "co_reactants": [{"label": "oxygen(2)", "origin": "initial_scenario_species"}],
        })
        evidence["candidate_routes"][1].update({
            "resolved_endpoints": {
                "reactants": [{"label": "stability(1)"}, {"label": "R(6)"}],
                "products": [{"label": "P2(5)"}],
            },
            "co_reactants": [{"label": "R(6)", "origin": "generated_intermediate"}],
        })
        # Make the downstream reaction controlling and its direct initiator small.
        flux = evidence["independent_cantera_propagation"]["reaction_flux_attribution"]
        flux["reactions"][0]["integrated_target_destruction_kmol"] = 2.0
        flux["reactions"][1]["integrated_target_destruction_kmol"] = 8.0
        conditions = {**CONDITIONS, "composition": {"stability": 0.01, "oxygen": 0.99}}
        plan = build_flux_ranked_verification_plan(evidence, conditions)
        self.assertEqual(plan["controlling_route_id"], "r2")
        self.assertEqual(plan["selected_route_id"], "r1")
        self.assertEqual(plan["ordered_dependency_route_ids"], ["r1", "r2"])

    def test_non_target_network_route_can_supply_generated_prerequisite(self):
        evidence = rmg_evidence()
        main = {
            "route_id": "target-main", "reaction_equation": "stability(1) + R(6) => P2(5)",
            "initiation_status": "thermally_reachable",
            "co_reactants": [{"label": "R(6)", "origin": "generated_intermediate"}],
            "resolved_endpoints": {
                "reactants": [{"label": "stability(1)"}, {"label": "R(6)"}],
                "products": [{"label": "P2(5)"}],
            },
        }
        evidence["candidate_routes"] = [main]
        evidence["network_routes"] = [
            {**main, "route_id": "duplicate-target-core"},
            {
                "route_id": "network-radical-source", "reaction_equation": "A(2) => R(6)",
                "resolved_endpoints": {
                    "reactants": [{"label": "A(2)"}], "products": [{"label": "R(6)"}],
                },
            },
        ]
        flux = evidence["independent_cantera_propagation"]["reaction_flux_attribution"]
        flux.update(total_integrated_target_destruction_kmol=10.0)
        flux["reactions"] = [
            {"reaction_index": 0, "reaction_equation": main["reaction_equation"],
             "integrated_target_destruction_kmol": 10.0,
             "integrated_forward_extent_kmol": 10.0, "integrated_reverse_extent_kmol": 0.0},
            {"reaction_index": 1, "reaction_equation": "A(2) => R(6)",
             "integrated_target_destruction_kmol": 0.0,
             "integrated_forward_extent_kmol": 1.0, "integrated_reverse_extent_kmol": 0.0},
        ]
        conditions = {**CONDITIONS, "composition": {"stability": 0.01, "A": 0.99}}
        plan = build_flux_ranked_verification_plan(evidence, conditions)
        self.assertEqual(plan["controlling_route_id"], "target-main")
        self.assertEqual(plan["selected_route_id"], "network-radical-source")
        self.assertEqual(plan["selected_route_source"], "network_routes")
        self.assertEqual(plan["selected_route_index"], 1)
        self.assertEqual(plan["ordered_dependency_route_ids"], ["network-radical-source", "target-main"])

    def test_deferred_producer_backtracks_to_alternative_radical_source(self):
        evidence = rmg_evidence()
        main = {
            "route_id": "target-main",
            "reaction_equation": "stability(1) + R(6) => P2(5)",
            "initiation_status": "thermally_reachable",
            "co_reactants": [{"label": "R(6)", "origin": "generated_intermediate"}],
            "resolved_endpoints": {
                "reactants": [{"label": "stability(1)"}, {"label": "R(6)"}],
                "products": [{"label": "P2(5)"}],
            },
        }
        evidence["candidate_routes"] = [main]
        evidence["network_routes"] = [
            {
                "route_id": "strong-source",
                "reaction_equation": "A(2) => R(6)",
                "resolved_endpoints": {
                    "reactants": [{"label": "A(2)"}],
                    "products": [{"label": "R(6)"}],
                },
            },
            {
                "route_id": "alternate-source",
                "reaction_equation": "B(7) => R(6)",
                "resolved_endpoints": {
                    "reactants": [{"label": "B(7)"}],
                    "products": [{"label": "R(6)"}],
                },
            },
        ]
        flux = evidence["independent_cantera_propagation"]["reaction_flux_attribution"]
        flux.update(total_integrated_target_destruction_kmol=10.0)
        flux["reactions"] = [
            {
                "reaction_index": 0,
                "reaction_equation": main["reaction_equation"],
                "integrated_target_destruction_kmol": 10.0,
                "integrated_forward_extent_kmol": 10.0,
                "integrated_reverse_extent_kmol": 0.0,
            },
            {
                "reaction_index": 1,
                "reaction_equation": "A(2) => R(6)",
                "integrated_target_destruction_kmol": 0.0,
                "integrated_forward_extent_kmol": 2.0,
                "integrated_reverse_extent_kmol": 0.0,
            },
            {
                "reaction_index": 2,
                "reaction_equation": "B(7) => R(6)",
                "integrated_target_destruction_kmol": 0.0,
                "integrated_forward_extent_kmol": 1.0,
                "integrated_reverse_extent_kmol": 0.0,
            },
        ]
        conditions = {
            **CONDITIONS,
            "composition": {"stability": 0.01, "A": 0.49, "B": 0.50},
        }
        plan = build_flux_ranked_verification_plan(
            evidence, conditions, deferred_route_ids=["strong-source"],
        )
        self.assertEqual(plan["status"], "verification_required")
        self.assertEqual(plan["controlling_route_id"], "target-main")
        self.assertEqual(plan["selected_route_id"], "alternate-source")
        self.assertEqual(
            plan["ordered_dependency_route_ids"], ["alternate-source", "target-main"],
        )
        producer_search = plan["dependency_branch_search"][0]["producer_search"][0]
        self.assertEqual(producer_search["deferred_producer_route_ids"], ["strong-source"])

    def test_blocked_high_flux_branch_falls_back_to_independent_loss_route(self):
        evidence = rmg_evidence()
        blocked = {
            "route_id": "blocked-main",
            "reaction_equation": "stability(1) + R(6) => P1(3)",
            "initiation_status": "thermally_reachable",
            "co_reactants": [{"label": "R(6)", "origin": "generated_intermediate"}],
            "resolved_endpoints": {
                "reactants": [{"label": "stability(1)"}, {"label": "R(6)"}],
                "products": [{"label": "P1(3)"}],
            },
        }
        independent = {
            "route_id": "independent-loss",
            "reaction_equation": "stability(1) + A(2) => P2(5)",
            "initiation_status": "directly_present",
            "co_reactants": [{"label": "A(2)", "origin": "initial_scenario_species"}],
            "resolved_endpoints": {
                "reactants": [{"label": "stability(1)"}, {"label": "A(2)"}],
                "products": [{"label": "P2(5)"}],
            },
        }
        evidence["candidate_routes"] = [blocked, independent]
        evidence["network_routes"] = [{
            "route_id": "failed-source",
            "reaction_equation": "A(2) => R(6)",
            "resolved_endpoints": {
                "reactants": [{"label": "A(2)"}],
                "products": [{"label": "R(6)"}],
            },
        }]
        flux = evidence["independent_cantera_propagation"]["reaction_flux_attribution"]
        flux.update(total_integrated_target_destruction_kmol=10.0)
        flux["reactions"] = [
            {
                "reaction_index": 0,
                "reaction_equation": blocked["reaction_equation"],
                "integrated_target_destruction_kmol": 9.0,
            },
            {
                "reaction_index": 1,
                "reaction_equation": independent["reaction_equation"],
                "integrated_target_destruction_kmol": 1.0,
            },
            {
                "reaction_index": 2,
                "reaction_equation": "A(2) => R(6)",
                "integrated_target_destruction_kmol": 0.0,
                "integrated_forward_extent_kmol": 1.0,
                "integrated_reverse_extent_kmol": 0.0,
            },
        ]
        conditions = {**CONDITIONS, "composition": {"stability": 0.01, "A": 0.99}}
        plan = build_flux_ranked_verification_plan(
            evidence, conditions, deferred_route_ids=["failed-source"],
        )
        self.assertEqual(plan["status"], "verification_required")
        self.assertEqual(plan["highest_flux_unverified_route_id"], "blocked-main")
        self.assertEqual(plan["controlling_route_id"], "independent-loss")
        self.assertEqual(plan["selected_route_id"], "independent-loss")
        self.assertEqual(
            plan["dependency_branch_search"][0]["status"], "branch_unresolved",
        )
        self.assertEqual(plan["dependency_branch_search"][1]["status"], "viable")

    def test_verified_flux_coverage_is_recorded_without_pretending_unresolved_is_robust(self):
        plan = build_flux_ranked_verification_plan(
            rmg_evidence(), CONDITIONS, verified_replacements=[replacement("r1")],
            flux_coverage_target=0.75,
        )
        self.assertEqual(plan["status"], "flux_coverage_target_met")
        self.assertAlmostEqual(plan["verified_flux_fraction"], 0.8)
        self.assertIsNone(plan["selected_route_id"])
        robust = assess_t95_robustness(rmg_evidence(), plan, retention_fraction=0.95)
        self.assertEqual(robust["status"], "not_robust")
        self.assertEqual(robust["reason"], "combined_unresolved_route_sensitivity_missing")

    def test_unattributed_flux_prevents_false_coverage(self):
        evidence = rmg_evidence()
        evidence["independent_cantera_propagation"]["reaction_flux_attribution"]["reactions"][1][
            "reaction_equation"
        ] = "stability(1) + UNKNOWN(9) <=> P3(10)"
        plan = build_flux_ranked_verification_plan(
            evidence, CONDITIONS, verified_replacements=[replacement("r1")], flux_coverage_target=0.95,
        )
        self.assertEqual(plan["status"], "insufficient_route_attribution")
        self.assertAlmostEqual(plan["unattributed_flux_fraction"], 0.2)
        self.assertIsNone(plan["selected_route_id"])

    def test_bad_flux_numerical_closure_fails_closed(self):
        evidence = rmg_evidence()
        evidence["independent_cantera_propagation"]["reaction_flux_attribution"][
            "numerical_closure_relative_error"
        ] = 0.9
        plan = build_flux_ranked_verification_plan(evidence, CONDITIONS)
        self.assertEqual(plan["status"], "insufficient_solver_evidence")
        self.assertIn("reaction_flux_numerical_closure_failed", plan["blocking_reasons"])

    def test_negligible_absolute_closure_error_does_not_fail_reversible_no_loss_model(self):
        evidence = rmg_evidence()
        flux = evidence["independent_cantera_propagation"]["reaction_flux_attribution"]
        flux.update({
            "numerical_closure_relative_error": 0.9,
            "initial_target_amount_kmol": 1.0e-2,
            "amount_profile_target_loss_kmol": 2.0e-15,
            "integrated_net_target_loss_kmol": 1.0e-15,
            "numerical_closure_absolute_error_kmol": 1.0e-15,
        })
        plan = build_flux_ranked_verification_plan(evidence, CONDITIONS)
        self.assertNotEqual(plan["status"], "insufficient_solver_evidence")

    def test_invalid_rmg_screening_promotes_only_rmg_direct_unimolecular_route(self):
        evidence = rmg_evidence()
        evidence["kinetics_validation"] = {"status": "kinetics_unreliable", "violation_count": 1}
        evidence["candidate_routes"].append({
            "route_id": "direct-rmg", "reaction_equation": "stability(1) <=> A(2) + B(3)",
            "initiation_status": "directly_present",
            "resolved_endpoints": {
                "reactants": [{"label": "stability(1)"}],
                "products": [{"label": "A(2)"}, {"label": "B(3)"}],
            },
        })
        plan = build_flux_ranked_verification_plan(evidence, CONDITIONS)
        self.assertEqual(plan["status"], "verification_required")
        self.assertEqual(plan["selected_route_id"], "direct-rmg")
        self.assertTrue(plan["intrinsic_failover_required"])

    def test_bimolecular_replacement_must_be_collision_bounded_at_condition(self):
        passed = validate_rate_replacement_evidence(
            replacement("r1"), temperature_K=298.15, pressure_bar=1.0,
        )
        self.assertEqual(passed["status"], "passed")
        failed = validate_rate_replacement_evidence(
            replacement("r1", rate=3.0e9, collision=2.0e9), temperature_K=298.15, pressure_bar=1.0,
        )
        self.assertEqual(failed["status"], "failed")
        self.assertIn("replacement_rate_exceeds_collision_limit", failed["blocking_reasons"])

    def test_combined_sensitivity_bounds_can_establish_robust_t95(self):
        evidence = rmg_evidence()
        plan = build_flux_ranked_verification_plan(
            evidence, CONDITIONS, verified_replacements=[replacement("r1")], flux_coverage_target=0.75,
        )
        envelope = {
            "status": "completed", "combined_perturbation": True,
            "unresolved_route_ids": ["r2"],
            "rate_bounds": [{"route_id": "r2", "lower_rate": 0.0, "upper_rate": 2.0e9,
                             "bound_source": "retained_collision_model"}],
            "lower_bound_propagation": propagation(t95=108.0),
            "upper_bound_propagation": propagation(t95=94.0),
        }
        result = assess_t95_robustness(
            evidence, plan, retention_fraction=0.95,
            unresolved_sensitivity_envelope=envelope, relative_tolerance=0.20,
        )
        self.assertEqual(result["status"], "robust_verified_t95")
        self.assertEqual(result["t95_seconds"], 100.0)

    @patch("storca.rmg_bridge_client.run_rmg_bridge")
    def test_combined_sensitivity_uses_only_explicit_retained_rate_bounds(self, bridge):
        evidence = rmg_evidence()
        evidence["artifacts"] = {"chemkin": "chem.inp", "species_dictionary": "dictionary.txt"}
        plan = build_flux_ranked_verification_plan(
            evidence, CONDITIONS, verified_replacements=[replacement("r1")], flux_coverage_target=0.75,
        )
        bridge.side_effect = [propagation(t95=108.0), propagation(t95=94.0)]
        condition = {**CONDITIONS, "composition": {"stability": 0.01, "A": 0.99},
                     "target_label": "stability", "target_duration_seconds": 1000.0}
        result = run_combined_unresolved_sensitivity(
            evidence, plan, condition, rmg_env="rmg_env",
            rate_bounds=[{"route_id": "r2", "baseline_rate": 4.0, "lower_rate": 0.0,
                          "upper_rate": 8.0, "units": "s-1", "bound_source": "retained_bound"}],
        )
        self.assertEqual(result["status"], "completed")
        self.assertEqual(bridge.call_args_list[0].args[0]["reaction_multipliers"], {"1": 0.0})
        self.assertEqual(bridge.call_args_list[1].args[0]["reaction_multipliers"], {"1": 2.0})

    def test_repaired_mechanism_rejects_unloaded_replacement_and_collision_warning(self):
        evidence = rmg_evidence()
        evidence["kinetics_validation"] = {"status": "kinetics_unreliable", "violation_count": 1}
        result = assess_repaired_mechanism(evidence, [replacement("r1")], CONDITIONS)
        self.assertEqual(result["status"], "rejected")
        self.assertIn("validated_replacement_library_not_present_in_repaired_model", result["blocking_reasons"])
        self.assertIn("collision_limit_violation_remains", result["blocking_reasons"])
        self.assertIsNone(result["estimated_time_to_retention_seconds"])

    def test_exact_upstream_repair_can_rerank_with_downstream_collision_warning(self):
        evidence = rmg_evidence()
        item = replacement("r1")
        evidence["generated_kinetics_libraries"] = [{"route_id": "r1"}]
        evidence["kinetics_validation"] = {"status": "kinetics_unreliable", "violation_count": 1}
        evidence["mechanism_inspection"] = {"reactions": [{
            "reaction_index": 0,
            "equation": item["reaction_equation"],
            "molecularity": 2,
            "evaluated_forward_rate_coefficient_si": 1.0e9 * 1e-6 * 6.02214076e23,
            "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
            "source_library": None,
            "kinetics_comment": None,
        }]}
        result = assess_repaired_mechanism(evidence, [item], CONDITIONS)
        self.assertEqual(
            result["status"],
            "accepted_for_flux_reranking_with_unrepaired_collision_violations",
        )
        self.assertFalse(result["terminal_kinetics_eligible"])
        self.assertIn("collision_limit_violation_remains", result["blocking_reasons"])

    def test_applied_replacement_requires_matching_signature_direction_and_rate(self):
        item = {
            "route_id": "uni", "reaction_equation": "stability(1) => P(2)",
            "verification_status": "verified_orca_irc_arkane", "molecularity": 1,
            "condition_rate": {"value": 4.0, "units": "s^-1", "temperature_K": 298.15,
                               "pressure_bar": 1.0},
        }
        evidence = {"mechanism_inspection": {"reactions": [{
            "reaction_index": 3, "equation": "stability(1) => P(2)", "molecularity": 1,
            "evaluated_forward_rate_coefficient_si": 4.0,
            "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
            "source_library": "repair", "kinetics_comment": None,
        }]}}
        result = verify_applied_replacement(evidence, item, CONDITIONS)
        self.assertEqual(result["status"], "passed")
        evidence["mechanism_inspection"]["reactions"][0]["evaluated_forward_rate_coefficient_si"] = 8.0
        result = verify_applied_replacement(evidence, item, CONDITIONS)
        self.assertEqual(result["status"], "failed")
        self.assertIn("applied_rate_does_not_match_replacement", result["blocking_reasons"])

    def test_applied_replacement_can_match_rmg_canonical_reverse_rate(self):
        item = {
            "route_id": "reverse", "reaction_equation": "stability(1) => P(2)",
            "verification_status": "verified_orca_irc_arkane", "molecularity": 1,
            "condition_rate": {"value": 4.0, "units": "s^-1", "temperature_K": 298.15,
                               "pressure_bar": 1.0},
        }
        evidence = {"mechanism_inspection": {"reactions": [{
            "reaction_index": 3, "equation": "P(2) <=> stability(1)",
            "molecularity": 1, "reverse_molecularity": 1,
            "evaluated_forward_rate_coefficient_si": 9.0,
            "evaluated_reverse_rate_coefficient_si": 4.0,
            "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
            "source_library": None, "kinetics_comment": None,
        }]}}
        result = verify_applied_replacement(evidence, item, CONDITIONS)
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["applied_direction"], "reverse")


if __name__ == "__main__":
    unittest.main()
