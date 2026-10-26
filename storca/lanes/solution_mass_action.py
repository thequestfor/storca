"""Executable verified-input homogeneous-solution mass-action lane."""

from __future__ import annotations

import json
from pathlib import Path

from .base import LaneDefinition, LaneResult
from .execution import LaneExecutionOutcome, LaneExecutionRequest
from ..solution_mass_action import (
    MassActionSolutionEnvironment,
    VerifiedMassActionNetwork,
    assess_mass_action_applicability,
    propagate_verified_mass_action_network,
)


SOLUTION_MASS_ACTION_CAPABILITIES = frozenset({
    "implicit_solution_environment_contract",
    "ideal_molar_activity_contract",
    "rdkit_solution_species_identity_validation",
    "verified_solution_first_and_second_order_rate_evidence",
    "dimensional_solution_rate_validation",
    "bimolecular_diffusion_ceiling_validation",
    "reversible_nonlinear_mass_action_propagation",
    "stiffness_diagnostic_and_linearly_implicit_propagation",
    "target_reformation_and_product_metrics",
    "element_charge_nonnegativity_extent_and_flux_closure_gates",
    "deterministic_rate_bound_sensitivity",
    "solution_t95_observable_gate",
})


def implicit_solution_mass_action_lane_definition() -> LaneDefinition:
    return LaneDefinition(
        lane_id="solution-neutral-implicit-mass-action-v1",
        environment_kind="homogeneous_solution",
        support_level="observable_supported",
        stage_kind="verified_input_solution_mass_action_kinetics",
        applicability=(
            "Dark, neutral, homogeneous constant-volume implicit solution with CPCM or SMD, "
            "ideal molar activities, fixed supplied species/concentrations, verified first- or "
            "second-order elementary rates, evidenced bimolecular diffusion ceilings, and a "
            "complete-network claim within the declared mass-action scope."
        ),
        discovery_backends=("external_verified_mass_action_network_bundle",),
        path_verifier="retained_external_condition_matched_solution_path_or_rate_claims",
        rate_model="verified_first_or_second_order_solution_rate_with_dimensions_and_diffusion_diagnostic",
        propagator="native_reversible_nonlinear_mass_action_with_audited_stiff_method_escalation",
        convergence_policy="network_coverage_applicability_rate_bound_scenario_element_charge_extent_flux_and_numerical_closure",
        supported_observables=(
            "t95", "retention_at_duration", "concentration_at_time",
            "product_yield", "dominant_product", "branching_fraction", "pathway_flux",
        ),
        required_capabilities=tuple(sorted(SOLUTION_MASS_ACTION_CAPABILITIES)),
        required_inputs=("solution_environment", "verified_solution_mass_action_network"),
        benchmark_profile="implicit_solution_mass_action_analytical_regression_v1",
        metadata={
            "scientific_support_scope": (
                "No autonomous solution discovery, gas-rate conversion, activity correction, dynamic "
                "speciation, explicit solvent participation, photochemistry, surface chemistry, or "
                "unresolved electronic-state chemistry. Schema v1 standardizes t95, retention, "
                "concentration, product yield, dominant product, gross target-consumption branching, "
                "and time-integrated pathway extent; half-life, selectivity, and induction time remain open."
            ),
        },
    )


class ImplicitSolutionMassActionLaneExecutor:
    @property
    def definition(self) -> LaneDefinition:
        return implicit_solution_mass_action_lane_definition()

    @property
    def capabilities(self) -> frozenset[str]:
        return SOLUTION_MASS_ACTION_CAPABILITIES

    def execute(self, request: LaneExecutionRequest, workdir: Path) -> LaneExecutionOutcome:
        environment = request.inputs.get("solution_environment")
        network = request.inputs.get("verified_solution_mass_action_network")
        if not isinstance(environment, MassActionSolutionEnvironment):
            raise TypeError("solution_environment must be a MassActionSolutionEnvironment")
        if not isinstance(network, VerifiedMassActionNetwork):
            raise TypeError("verified_solution_mass_action_network must be a VerifiedMassActionNetwork")
        if environment.as_environment_spec().as_dict() != request.environment.as_dict():
            raise ValueError("Lane request environment does not match mass-action environment evidence")
        contract = self.definition.as_dict()
        applicability = assess_mass_action_applicability(environment, network)
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        evidence_path = workdir / "solution-mass-action-lane-evidence.json"

        if applicability["status"] != "supported":
            result = LaneResult(
                lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
                admission_status="supported", support_level=self.definition.support_level,
                execution_status="outside_method_applicability", observable_status="unresolved",
                requested_observable=request.requested_observable,
                reason="The supplied network contains unresolved or unsupported solution chemistry.",
            )
            evidence = {
                "applicability": applicability,
                "scientific_claim": {
                    "status": "unresolved", "claim": "No mass-action solution t95 is defensible.",
                    "blocking_findings": applicability["blocking_findings"],
                },
            }
            evidence_path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
            return LaneExecutionOutcome(result, evidence, {"lane_evidence": str(evidence_path)})

        try:
            propagation = propagate_verified_mass_action_network(
                environment, network,
                profile_points=int(request.compute_budget.get("profile_points", 180)),
                relative_tolerance=float(request.compute_budget.get("relative_tolerance", 1.0e-7)),
                absolute_tolerance=float(request.compute_budget.get("absolute_tolerance", 1.0e-12)),
                maximum_steps=int(request.compute_budget.get("maximum_steps", 200_000)),
            )
        except Exception as error:
            result = LaneResult(
                lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
                admission_status="supported", support_level=self.definition.support_level,
                execution_status="numerical_failure", observable_status="unresolved",
                requested_observable=request.requested_observable,
                reason=f"Nonlinear concentration propagation failed closed: {type(error).__name__}: {error}",
            )
            evidence = {
                "applicability": applicability,
                "propagation": {"status": "failed", "reason": result.reason},
                "scientific_claim": {"status": "unresolved", "claim": "No mass-action solution t95 is available."},
            }
            evidence_path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
            return LaneExecutionOutcome(result, evidence, {"lane_evidence": str(evidence_path)})

        uncertainty = propagation["rate_uncertainty_sensitivity"]
        t95 = propagation["estimated_time_to_retention_seconds"]
        if request.requested_observable != "t95":
            try:
                result, claim = self._non_t95_result(
                    request, network, propagation, uncertainty, contract,
                )
            except Exception as error:
                result = LaneResult(
                    lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
                    admission_status="supported", support_level=self.definition.support_level,
                    execution_status="invalid_observable_request", observable_status="unresolved",
                    requested_observable=request.requested_observable,
                    reason=f"Observable request failed closed: {type(error).__name__}: {error}",
                )
                claim = {
                    "status": "unresolved", "claim": "No requested solution observable is defensible.",
                    "reason": result.reason,
                }
            evidence = {
                "applicability": applicability, "propagation": propagation,
                "scientific_claim": claim,
            }
            evidence_path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
            return LaneExecutionOutcome(result, evidence, {"lane_evidence": str(evidence_path)})
        if uncertainty["status"] == "crossing_depends_on_rate_uncertainty":
            result = LaneResult(
                lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
                admission_status="supported", support_level=self.definition.support_level,
                execution_status="observable_unresolved_rate_uncertainty", observable_status="unresolved",
                requested_observable=request.requested_observable,
                reason="Whether t95 occurs changes across retained rate-bound scenarios.",
            )
            claim = {
                "status": "unresolved", "claim": "No within-horizon t95 is robust across retained rate bounds.",
                "rate_uncertainty_sensitivity": uncertainty,
            }
        elif t95 is None:
            result = LaneResult(
                lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
                admission_status="supported", support_level=self.definition.support_level,
                execution_status="completed_no_t95_within_horizon",
                observable_status="not_reached_within_horizon",
                requested_observable=request.requested_observable,
                reason="No first t95 crossing occurred in any retained rate-bound scenario.",
            )
            claim = {
                "status": "model_supported_no_crossing_within_horizon",
                "claim": "No t95 crossing occurred within the declared duration in the admitted model.",
                "target_fraction_remaining": propagation["target_fraction_remaining"],
                "conditional_on_network_coverage_claim_id": network.coverage_claim_id,
            }
        elif uncertainty["status"] != "crossing_within_horizon_for_all_retained_scenarios":
            result = LaneResult(
                lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
                admission_status="supported", support_level=self.definition.support_level,
                execution_status="observable_unresolved_uncertainty_inconsistency",
                observable_status="unresolved", requested_observable=request.requested_observable,
                reason="Central t95 is inconsistent with the retained uncertainty classification.",
            )
            claim = {"status": "unresolved", "claim": "No robust mass-action solution t95 is available."}
        else:
            result = LaneResult(
                lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
                admission_status="supported", support_level=self.definition.support_level,
                execution_status="completed_verified_t95", observable_status="verified",
                requested_observable=request.requested_observable,
                observable_value=float(t95), observable_units="s",
                reason="Verified mass-action network crossed t95 in every retained rate-bound scenario.",
            )
            claim = {
                "status": "verified_condition_specific_t95",
                "claim": "The admitted verified-input mass-action model predicts the retained t95.",
                "t95_seconds": float(t95),
                "t95_bounds_seconds": [
                    uncertainty["lower_crossing_seconds"], uncertainty["upper_crossing_seconds"],
                ],
                "environment_contract_id": request.environment.as_dict()["contract_id"],
                "network_coverage_claim_id": network.coverage_claim_id,
                "rate_claim_ids": [item.rate_claim_id for item in network.reactions],
                "path_claim_ids": [item.path_claim_id for item in network.reactions if item.path_claim_id],
                "limitations": self.definition.metadata["scientific_support_scope"],
            }
        evidence = {"applicability": applicability, "propagation": propagation, "scientific_claim": claim}
        evidence_path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
        return LaneExecutionOutcome(result, evidence, {"lane_evidence": str(evidence_path)})

    def _non_t95_result(
        self,
        request: LaneExecutionRequest,
        network: VerifiedMassActionNetwork,
        propagation: dict,
        uncertainty: dict,
        contract: dict,
    ) -> tuple[LaneResult, dict]:
        kind = request.requested_observable
        parameters = request.inputs.get("observable_parameters") or {}
        if not isinstance(parameters, dict):
            raise TypeError("observable_parameters must be a mapping")
        scenarios = uncertainty["scenarios"]
        central_value = None
        units = None
        observable_data = {}
        scenario_values = []

        if kind == "retention_at_duration":
            species_id = str(parameters.get("species_id") or network.target_species_id)
            if species_id != network.target_species_id:
                raise ValueError("retention_at_duration is defined for the stability target only")
            central_value = float(propagation["target_fraction_remaining"])
            scenario_values = [float(item["target_fraction_remaining"]) for item in scenarios]
            units = "fraction"
        elif kind == "concentration_at_time":
            species_id = str(parameters.get("species_id") or "").strip()
            if species_id not in propagation["final_concentrations_mol_L"]:
                raise ValueError("concentration_at_time requires a retained species_id")
            central_value = float(propagation["final_concentrations_mol_L"][species_id])
            scenario_values = [float(item["final_concentrations_mol_L"][species_id]) for item in scenarios]
            units = "mol L^-1"
            observable_data = {
                "species_id": species_id,
                **({"component_id": str(parameters["component_id"])} if parameters.get("component_id") else {}),
            }
        elif kind == "product_yield":
            species_id = str(parameters.get("species_id") or "").strip()
            yields = propagation["product_yields_relative_to_initial_target"]
            if species_id == network.target_species_id or species_id not in yields:
                raise ValueError("product_yield requires a retained non-target species_id")
            central_value = float(yields[species_id])
            scenario_values = [
                float(item["product_yields_relative_to_initial_target"][species_id])
                for item in scenarios
            ]
            units = "mol product per mol initial target"
            observable_data = {
                "species_id": species_id,
                **({"component_id": str(parameters["component_id"])} if parameters.get("component_id") else {}),
            }
        elif kind == "dominant_product":
            candidates = tuple(str(item) for item in parameters.get("candidate_species_ids") or ())
            available = set(propagation["product_yields_relative_to_initial_target"])
            if candidates and (set(candidates) - available):
                raise ValueError("dominant_product candidate_species_ids include an unknown species")
            considered = candidates or tuple(sorted(available))

            def dominant(values: dict) -> str | None:
                if not considered or not any(float(values[item]) > 0.0 for item in considered):
                    return None
                return max(considered, key=lambda item: (float(values[item]), item))

            central_product = dominant(propagation["product_yields_relative_to_initial_target"])
            scenario_products = [
                dominant(item["product_yields_relative_to_initial_target"])
                for item in scenarios
            ]
            if any(item != central_product for item in scenario_products):
                result = LaneResult(
                    lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
                    admission_status="supported", support_level=self.definition.support_level,
                    execution_status="observable_unresolved_rate_uncertainty",
                    observable_status="unresolved", requested_observable=kind,
                    reason="Dominant product changes across retained rate-bound scenarios.",
                )
                return result, {
                    "status": "unresolved", "claim": result.reason,
                    "central_species_id": central_product,
                    "scenario_species_ids": scenario_products,
                }
            observable_data = {
                "species_id": central_product,
                "component_id": (
                    dict(parameters.get("species_to_component_id") or {}).get(central_product)
                    if central_product else None
                ),
                "status": "dominant_positive_gain" if central_product else "no_positive_product_gain",
                "candidate_species_ids": list(considered),
            }
        elif kind in {"branching_fraction", "pathway_flux"}:
            reaction_id = str(parameters.get("reaction_id") or "").strip()
            central_reactions = {
                item["reaction_id"]: item
                for item in propagation["reaction_flux_attribution"]["reactions"]
            }
            if reaction_id not in central_reactions:
                raise ValueError(f"{kind} requires a retained reaction_id")
            if kind == "branching_fraction":
                central_value = float(
                    central_reactions[reaction_id]["gross_target_consumption_fraction"],
                )
                scenario_values = [
                    float(item["gross_target_consumption_fractions"][reaction_id])
                    for item in scenarios
                ]
                units = "fraction of gross target consumption"
            else:
                if parameters.get("metric", "integrated_extent") != "integrated_extent":
                    raise ValueError("pathway_flux schema v1 supports metric='integrated_extent' only")
                central_value = float(central_reactions[reaction_id]["integrated_extent_mol_L"])
                scenario_values = [
                    float(item["integrated_reaction_extents_mol_L"][reaction_id])
                    for item in scenarios
                ]
                units = "mol L^-1 integrated extent"
            observable_data = {"reaction_id": reaction_id}
        else:
            raise ValueError(f"Unsupported mass-action observable: {kind}")

        if central_value is not None:
            observable_data = {
                **observable_data,
                "rate_bound_scenario_min": min(scenario_values),
                "rate_bound_scenario_max": max(scenario_values),
                "scenario_count": len(scenario_values),
            }
        result = LaneResult(
            lane_id=self.definition.lane_id, lane_contract_id=contract["contract_id"],
            admission_status="supported", support_level=self.definition.support_level,
            execution_status=f"completed_verified_{kind}", observable_status="verified",
            requested_observable=kind, observable_value=central_value,
            observable_units=units, observable_data=observable_data,
            reason="Requested observable is retained with every explicit rate-bound scenario.",
        )
        claim = {
            "status": f"verified_condition_specific_{kind}",
            "claim": "The admitted verified-input mass-action model predicts the requested observable.",
            "observable_kind": kind, "value": central_value, "units": units,
            "observable_data": observable_data,
            "environment_contract_id": request.environment.as_dict()["contract_id"],
            "network_coverage_claim_id": network.coverage_claim_id,
            "rate_claim_ids": [item.rate_claim_id for item in network.reactions],
            "path_claim_ids": [item.path_claim_id for item in network.reactions if item.path_claim_id],
            "uncertainty_interpretation": uncertainty["interpretation"],
            "limitations": self.definition.metadata["scientific_support_scope"],
        }
        return result, claim


__all__ = [
    "ImplicitSolutionMassActionLaneExecutor", "SOLUTION_MASS_ACTION_CAPABILITIES",
    "implicit_solution_mass_action_lane_definition",
]
