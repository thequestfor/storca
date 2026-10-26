"""First executable homogeneous-solution stability lane."""

from __future__ import annotations

import json
from pathlib import Path

from .base import LaneDefinition, LaneResult
from .execution import LaneExecutionOutcome, LaneExecutionRequest
from ..solution_kinetics import (
    ImplicitSolutionEnvironment,
    VerifiedSolutionNetwork,
    assess_solution_network_applicability,
    propagate_verified_solution_network,
)


SOLUTION_UNIMOLECULAR_CAPABILITIES = frozenset({
    "implicit_solution_environment_contract",
    "rdkit_solution_species_identity_validation",
    "verified_solution_unimolecular_rate_evidence",
    "native_concentration_network",
    "linear_solution_concentration_propagation",
    "element_and_flux_closure_gates",
    "solution_t95_observable_gate",
})


def implicit_solution_lane_definition() -> LaneDefinition:
    return LaneDefinition(
        lane_id="solution-neutral-implicit-v1",
        environment_kind="homogeneous_solution",
        support_level="observable_supported",
        stage_kind="verified_input_solution_kinetics",
        applicability=(
            "Dark, neutral, homogeneous constant-volume implicit solution with a declared CPCM or SMD "
            "solvent, one-molar standard state, monotonic target loss, externally verified first-order "
            "elementary rates, and verified network coverage within the declared unimolecular scope."
        ),
        discovery_backends=(
            "external_verified_network_bundle", "storca_native_candidate_adapter",
            "rmg_resolved_candidate_adapter", "graph_xtb_candidate_adapter",
        ),
        path_verifier="storca_solvent_consistent_orca_path_or_retained_external_solution_path_claims",
        rate_model="verified_solution_eyring_tst_or_experimental_first_order",
        propagator="native_linear_concentration_matrix_exponential",
        convergence_policy="observable_conditioned_rate_verification_then_network_coverage_applicability_element_flux_and_numerical_closure",
        supported_observables=("t95",),
        required_capabilities=tuple(sorted(SOLUTION_UNIMOLECULAR_CAPABILITIES)),
        required_inputs=("solution_environment", "verified_solution_network"),
        benchmark_profile="implicit_solution_unimolecular_analytical_regression_v1",
        metadata={
            "required_input_alternatives": [
                ["solution_verification_session_request"],
            ],
            "scientific_support_scope": (
                "This lane accepts retained external rates or STORCA's solvent-consistent ORCA path/rate evidence; "
                "it is not autonomous solution reaction discovery. "
                "Any influential bimolecular, solvent-assisted, ionic, photochemical, spin-crossing, "
                "surface, or unresolved multireference chemistry blocks the observable."
            ),
        },
    )


class ImplicitSolutionLaneExecutor:
    @property
    def definition(self) -> LaneDefinition:
        return implicit_solution_lane_definition()

    @property
    def capabilities(self) -> frozenset[str]:
        return SOLUTION_UNIMOLECULAR_CAPABILITIES

    def execute(
        self,
        request: LaneExecutionRequest,
        workdir: Path,
    ) -> LaneExecutionOutcome:
        environment = request.inputs.get("solution_environment")
        network = request.inputs.get("verified_solution_network")
        if not isinstance(environment, ImplicitSolutionEnvironment):
            raise TypeError("solution_environment must be an ImplicitSolutionEnvironment")
        if not isinstance(network, VerifiedSolutionNetwork):
            raise TypeError("verified_solution_network must be a VerifiedSolutionNetwork")
        if environment.as_environment_spec().as_dict() != request.environment.as_dict():
            raise ValueError("Lane request environment does not match the solution environment evidence")

        contract = self.definition.as_dict()
        applicability = assess_solution_network_applicability(environment, network)
        Path(workdir).mkdir(parents=True, exist_ok=True)
        if applicability["status"] != "supported":
            result = LaneResult(
                lane_id=self.definition.lane_id,
                lane_contract_id=contract["contract_id"],
                admission_status="supported",
                support_level=self.definition.support_level,
                execution_status="outside_method_applicability",
                observable_status="unresolved",
                requested_observable=request.requested_observable,
                reason="The supplied network contains unresolved or unsupported solution chemistry.",
            )
            evidence = {
                "applicability": applicability,
                "scientific_claim": {
                    "status": "unresolved",
                    "claim": "No solution t95 is defensible for this evidence bundle.",
                    "blocking_findings": applicability["blocking_findings"],
                },
            }
            path = Path(workdir) / "solution-lane-evidence.json"
            path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
            return LaneExecutionOutcome(result, evidence, {"lane_evidence": str(path)})

        try:
            propagation = propagate_verified_solution_network(
                environment,
                network,
                profile_points=int(request.compute_budget.get("profile_points", 160)),
            )
        except Exception as error:
            result = LaneResult(
                lane_id=self.definition.lane_id,
                lane_contract_id=contract["contract_id"],
                admission_status="supported",
                support_level=self.definition.support_level,
                execution_status="numerical_failure",
                observable_status="unresolved",
                requested_observable=request.requested_observable,
                reason=f"Concentration propagation failed closed: {type(error).__name__}: {error}",
            )
            evidence = {
                "applicability": applicability,
                "propagation": {"status": "failed", "reason": result.reason},
                "scientific_claim": {"status": "unresolved", "claim": "No solution t95 is available."},
            }
            path = Path(workdir) / "solution-lane-evidence.json"
            path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
            return LaneExecutionOutcome(result, evidence, {"lane_evidence": str(path)})

        t95 = propagation.get("estimated_time_to_retention_seconds")
        uncertainty = propagation.get("t95_uncertainty") or {}
        if uncertainty.get("status") == "crossing_depends_on_rate_uncertainty":
            result = LaneResult(
                lane_id=self.definition.lane_id,
                lane_contract_id=contract["contract_id"],
                admission_status="supported",
                support_level=self.definition.support_level,
                execution_status="observable_unresolved_rate_uncertainty",
                observable_status="unresolved",
                requested_observable=request.requested_observable,
                reason="Whether t95 occurs within the horizon changes across the explicit rate bounds.",
            )
            evidence = {
                "applicability": applicability,
                "propagation": propagation,
                "scientific_claim": {
                    "status": "unresolved",
                    "claim": "No within-horizon solution t95 is defensible across the retained rate bounds.",
                    "t95_uncertainty": uncertainty,
                },
            }
            path = Path(workdir) / "solution-lane-evidence.json"
            path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
            return LaneExecutionOutcome(result, evidence, {"lane_evidence": str(path)})
        if t95 is None:
            result = LaneResult(
                lane_id=self.definition.lane_id,
                lane_contract_id=contract["contract_id"],
                admission_status="supported",
                support_level=self.definition.support_level,
                execution_status="completed_no_t95_within_horizon",
                observable_status="not_reached_within_horizon",
                requested_observable=request.requested_observable,
                reason="The verified-input model did not cross the requested retention threshold within the horizon.",
            )
            claim = {
                "status": "model_supported_no_crossing_within_horizon",
                "claim": "No t95 crossing occurred within the declared duration in the admitted model.",
                "target_fraction_remaining": propagation["target_fraction_remaining"],
                "conditional_on_network_coverage_claim_id": network.coverage_claim_id,
            }
        else:
            result = LaneResult(
                lane_id=self.definition.lane_id,
                lane_contract_id=contract["contract_id"],
                admission_status="supported",
                support_level=self.definition.support_level,
                execution_status="completed_verified_t95",
                observable_status="verified",
                requested_observable=request.requested_observable,
                observable_value=float(t95),
                observable_units="s",
                reason="Verified first-order solution network crossed the declared retention threshold.",
            )
            claim = {
                "status": "verified_condition_specific_t95",
                "claim": "The admitted verified-input implicit-solution model predicts the retained t95.",
                "t95_seconds": float(t95),
                "environment_contract_id": request.environment.as_dict()["contract_id"],
                "network_coverage_claim_id": network.coverage_claim_id,
                "t95_bounds_seconds": [uncertainty.get("lower_seconds"), uncertainty.get("upper_seconds")],
                "rate_claim_ids": [reaction.rate_claim_id for reaction in network.reactions],
                "path_claim_ids": [
                    reaction.path_claim_id for reaction in network.reactions if reaction.path_claim_id
                ],
                "limitations": self.definition.metadata["scientific_support_scope"],
            }
        evidence = {
            "applicability": applicability,
            "propagation": propagation,
            "scientific_claim": claim,
        }
        path = Path(workdir) / "solution-lane-evidence.json"
        path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
        return LaneExecutionOutcome(result, evidence, {"lane_evidence": str(path)})


__all__ = [
    "ImplicitSolutionLaneExecutor", "SOLUTION_UNIMOLECULAR_CAPABILITIES",
    "implicit_solution_lane_definition",
]
