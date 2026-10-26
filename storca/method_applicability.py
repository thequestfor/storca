"""Deterministic applicability decisions around the existing path/rate gates.

This policy neither chooses a quantum method nor estimates its error. It keeps
an executable conditional TST model separate from evidence that its physical
approximations apply to an observed reaction. Unavailable diagnostics are
retained as assumptions; known contradictions cannot be waived by a warning.
"""

from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Mapping


METHOD_APPLICABILITY_POLICY_ID = "ordinary-single-surface-tst-applicability-v1"
ENDPOINT_STATE_AUDIT_POLICY_ID = "prepath-coupled-endpoint-state-audit-v1"
_DIAGNOSTICS = (
    "electronic_structure", "conformer_sampling", "statistical_mechanics",
    "tunneling", "pressure_dependence",
)
_SOURCES = {
    "orca_electronic_structure": (
        "https://www.faccts.de/docs/orca/6.1/manual/contents/quickstartguide/troubleshooting.html"
        "#strange-results-for-energies-and-properties"
    ),
    "arkane_rate_models": (
        "https://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/introduction.html"
    ),
    "arkane_inputs": (
        "https://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/input.html"
    ),
}
_ENVIRONMENT_ALIASES = {
    "gas": "homogeneous_gas",
    "gas phase": "homogeneous_gas",
    "homogeneous gas": "homogeneous_gas",
    "homogeneous gas phase": "homogeneous_gas",
    "homogeneous gas phase pure target model": "homogeneous_gas",
    "homogeneous gas phase surrogate": "homogeneous_gas",
    "homogeneous gas phase low pressure surrogate": "homogeneous_gas",
    "implicit solution": "implicit_solution",
    "homogeneous implicit solution": "implicit_solution",
    "implicit solvent": "implicit_solution",
    "homogeneous implicit solvent": "implicit_solution",
    "pcm": "implicit_solution", "cpcm": "implicit_solution", "smd": "implicit_solution",
}
_SPECIALIST_FLAGS = {
    "explicit_solvent_participation": "solvent_participating_rate_model",
    "acid_base_speciation_required": "speciation_and_condition_matched_rates",
    "ion_pairing_required": "ion_pair_speciation_and_rates",
    "photochemical": "excited_state_or_photochemical_rate_model",
    "spin_crossing": "spin_crossing_or_nonadiabatic_rate_model",
    "surface_or_interface": "surface_or_interface_rate_model",
    "multireference_unresolved": "electronic_structure_assessment_or_escalation",
    "diffusion_limited": "encounter_diffusion_rate_model",
}


def _hash(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return "sha256:" + hashlib.sha256(encoded.encode()).hexdigest()


def method_applicability_policy() -> dict[str, Any]:
    """Return the versioned decision rules, distinct from accuracy calibration."""
    result = {
        "schema_version": 1,
        "kind": "storca_method_applicability_policy",
        "policy_id": METHOD_APPLICABILITY_POLICY_ID,
        "scope": "ordinary_ground_state_single_surface_harmonic_tst",
        "required_diagnostic_categories": list(_DIAGNOSTICS),
        "missing_diagnostic_action": "explicit_conditional_assumption",
        "known_contradiction_action": "withhold_ordinary_tst_rate_and_escalate",
        "warnings_are_proof_of_inapplicability": False,
        "absence_of_warnings_establishes_applicability": False,
        "rate_accuracy_calibration_provided": False,
        "sources": dict(_SOURCES),
    }
    return {**result, "scientific_hash": _hash(result)}


def _environment_model(value: Any) -> str | None:
    if value is None or not str(value).strip():
        return None
    normalized = " ".join(str(value).strip().casefold().replace("_", " ").replace("-", " ").split())
    return _ENVIRONMENT_ALIASES.get(normalized, normalized)


def _side_states(route: Mapping[str, Any], verification: Mapping[str, Any], side: str) -> tuple[int, ...]:
    records = (verification.get("arkane_inputs") or {}).get(side) or ()
    if records:
        try:
            return tuple(int(item["multiplicity"]) for item in records)
        except (KeyError, TypeError, ValueError):
            return ()
    key = "reactant_multiplicities" if side == "reactants" else "product_multiplicities"
    if route.get(key):
        try:
            return tuple(int(value) for value in route[key])
        except (TypeError, ValueError):
            return ()
    records = (route.get("resolved_endpoints") or {}).get(side) or ()
    try:
        return tuple(int(item["multiplicity"]) for item in records)
    except (KeyError, TypeError, ValueError):
        return ()


def endpoint_state_audit_policy() -> dict[str, Any]:
    """Expose the rules for a state audit that cannot establish a valid rate."""
    result = {
        "schema_version": 1,
        "kind": "storca_endpoint_state_audit_policy",
        "policy_id": ENDPOINT_STATE_AUDIT_POLICY_ID,
        "scope": "declared_fragment_spin_compatibility_before_path_search",
        "multiplicity_validation": "positive_integer_excluding_boolean_without_coercion",
        "coupling_rule": "route_geometry._coupled_multiplicities",
        "no_common_surface_action": "spin_changing_rate_model_unresolved",
        "multiple_common_surfaces_action": "require_explicit_surface_selection",
        "single_common_surface_action": "eligible_for_further_ordinary_path_review",
        "rate_accuracy_calibration_provided": False,
        "limitations": [
            "Endpoint multiplicities do not locate a crossing or establish its rate.",
            "Common total spin does not establish single-reference validity or path connectivity.",
            "No path, selected electronic state energy, or specialist calculation is fabricated.",
        ],
    }
    return {**result, "scientific_hash": _hash(result)}


def _positive_integer(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def _fragment_states(values: Any, side: str) -> tuple[tuple[int, ...], list[dict[str, Any]]]:
    """Retain validation evidence without coercing floats, strings, or booleans."""
    if not isinstance(values, (list, tuple)) or not values:
        return (), [{"side": side, "reason": "nonempty_fragment_state_sequence_required",
                     "input_type": type(values).__name__}]
    invalid = []
    for index, value in enumerate(values):
        if not _positive_integer(value):
            invalid.append({"side": side, "fragment_index": index,
                            "reason": "positive_integer_multiplicity_required",
                            "input_type": type(value).__name__,
                            "input_value": str(value) if isinstance(value, (str, int, float, bool)) else None})
    return (() if invalid else tuple(values)), invalid


def audit_endpoint_electronic_states(
    reactant_multiplicities: Any,
    product_multiplicities: Any,
    selected_multiplicity: Any = None,
) -> dict[str, Any]:
    """Audit fragment spin compatibility before choosing an ordinary path.

    A compatible result only makes further path/applicability review possible.
    It never admits an ordinary TST rate, validates a path, or authorizes a job.
    Fragment spin coupling uses the same rules as the route-geometry gate.
    """
    reactants, reactant_errors = _fragment_states(reactant_multiplicities, "reactants")
    products, product_errors = _fragment_states(product_multiplicities, "products")
    invalid = reactant_errors + product_errors
    selected_valid = selected_multiplicity is None or _positive_integer(selected_multiplicity)
    if not selected_valid:
        invalid.append({"side": "selected_surface", "reason": "positive_integer_multiplicity_required",
                        "input_type": type(selected_multiplicity).__name__,
                        "input_value": str(selected_multiplicity)
                        if isinstance(selected_multiplicity, (str, int, float, bool)) else None})
    reactant_surfaces: list[int] = []
    product_surfaces: list[int] = []
    if reactants and products:
        from .route_geometry import _coupled_multiplicities

        reactant_surfaces = _coupled_multiplicities(reactants)
        product_surfaces = _coupled_multiplicities(products)
    common = sorted(set(reactant_surfaces) & set(product_surfaces))
    specialist_flags: dict[str, bool] = {}
    required_actions: list[str] = []
    if invalid:
        status = "electronic_state_identity_unresolved"
        required_actions = ["retain_valid_endpoint_and_selected_state_declarations"]
    elif not common:
        status = "spin_changing_rate_model_unresolved"
        specialist_flags["spin_crossing"] = True
        required_actions = ["spin_crossing_or_nonadiabatic_rate_model"]
    elif selected_multiplicity is not None and selected_multiplicity not in common:
        status = "selected_surface_not_common_to_endpoints"
        required_actions = ["electronic_surface_resolution"]
    elif len(common) > 1 and selected_multiplicity is None:
        status = "surface_selection_required"
        required_actions = ["select_and_justify_common_electronic_surface"]
    else:
        status = "eligible_for_ordinary_path_review"
    eligible = status == "eligible_for_ordinary_path_review"
    policy = endpoint_state_audit_policy()
    result = {
        "schema_version": 1,
        "kind": "storca_prepath_endpoint_state_audit",
        "policy_id": policy["policy_id"], "policy_hash": policy["scientific_hash"],
        "status": status,
        "reactant_fragment_multiplicities": list(reactants) if reactants else None,
        "product_fragment_multiplicities": list(products) if products else None,
        "reactant_coupled_multiplicities": reactant_surfaces,
        "product_coupled_multiplicities": product_surfaces,
        "common_multiplicities": common,
        "selected_multiplicity": selected_multiplicity if selected_valid else None,
        "selection_basis": "explicit_declaration" if selected_multiplicity is not None and selected_valid
        else "unique_common_surface_for_review_only" if eligible else "not_selected",
        "invalid_state_evidence": invalid,
        "specialist_flags": specialist_flags,
        "required_actions": required_actions,
        "ordinary_path_eligible": eligible,
        "central_rate_allowed": False,
        "rate_model_validated": False,
        "validated_prediction_allowed": False,
        "lifetime_claim_allowed": False,
        "automatic_specialist_job_allowed": False,
        "limitations": list(policy["limitations"]),
    }
    return {**result, "scientific_hash": _hash(result)}


def audit_retained_network_states(network: Mapping[str, Any]) -> dict[str, Any]:
    """Audit exact snapshot-local species states and group native reaction refs.

    The caller must separately validate the retained network and its sources.
    This helper preserves full native identifiers and does not establish stereo,
    atom balance, thermochemistry, source attribution, or rate applicability.
    """
    policy = endpoint_state_audit_policy()
    problems: list[str] = []
    if not isinstance(network, Mapping):
        network = {}
        problems.append("network_mapping_required")
    try:
        network_input_hash = _hash(network)
    except (TypeError, ValueError):
        network_input_hash = None
        problems.append("network_input_hash_unavailable")
    species = network.get("species")
    if not isinstance(species, (list, tuple)) or not species:
        species = []
        problems.append("native_species_inventory_required")
    states: dict[str, Any] = {}
    duplicate_ids: set[str] = set()
    for item in species:
        if not isinstance(item, Mapping):
            problems.append("native_species_record_malformed")
            continue
        identifier = item.get("chemkin_identifier")
        if not isinstance(identifier, str) or not identifier:
            problems.append("full_native_species_identifier_required")
            continue
        if identifier in states:
            duplicate_ids.add(identifier)
        states[identifier] = item.get("multiplicity")
    for identifier in duplicate_ids:
        states[identifier] = None
    if duplicate_ids:
        problems.append("duplicate_native_species_identifiers")
    channels = network.get("channels")
    if not isinstance(channels, (list, tuple)) or not channels:
        channels = []
        problems.append("directed_channel_inventory_required")

    def side_states(side: Any) -> list[Any] | None:
        if not isinstance(side, Mapping) or not side:
            return None
        values: list[Any] = []
        for identifier, count in side.items():
            if not isinstance(identifier, str) or not _positive_integer(count):
                return None
            values.extend([states.get(identifier)] * count)
        return values

    channel_audits = []
    groups: dict[int, list[dict[str, Any]]] = {}
    for index, channel in enumerate(channels):
        if not isinstance(channel, Mapping):
            channel = {}
            problems.append("directed_channel_record_malformed")
        reaction_index = channel.get("reaction_index", channel.get("native_reaction_index"))
        valid_reference = (isinstance(reaction_index, int) and not isinstance(reaction_index, bool)
                           and reaction_index >= 0)
        if ("reaction_index" in channel and "native_reaction_index" in channel
                and channel["reaction_index"] != channel["native_reaction_index"]):
            valid_reference = False
        if not valid_reference:
            reaction_index = None
            problems.append("native_reaction_reference_unresolved")
        channel_id = channel.get("channel_id")
        if not isinstance(channel_id, str) or not channel_id:
            channel_id = None
            problems.append("directed_channel_identifier_required")
        direction = channel.get("direction")
        if not isinstance(direction, str) or direction not in {"forward", "reverse"}:
            direction = None
            problems.append("directed_channel_direction_required")
        audit = audit_endpoint_electronic_states(
            side_states(channel.get("reactants")), side_states(channel.get("products")),
            channel.get("selected_multiplicity"),
        )
        record = {
            "channel_id": channel_id,
            "channel_position": index,
            "reaction_index": reaction_index,
            "direction": direction,
            "reactants": dict(channel["reactants"]) if side_states(channel.get("reactants")) is not None else None,
            "products": dict(channel["products"]) if side_states(channel.get("products")) is not None else None,
            "status": audit["status"],
            "ordinary_path_eligible": audit["ordinary_path_eligible"] and valid_reference and not problems,
            "audit": audit,
        }
        channel_audits.append(record)
        if valid_reference:
            groups.setdefault(reaction_index, []).append(record)
    if problems:
        for record in channel_audits:
            record["ordinary_path_eligible"] = False
    reaction_audits = []
    for reaction_index, records in sorted(groups.items()):
        flags = {flag: True for record in records for flag, value in record["audit"]["specialist_flags"].items() if value}
        reaction_audits.append({
            "reaction_index": reaction_index,
            "channel_ids": [record["channel_id"] for record in records],
            "statuses": sorted({record["status"] for record in records}),
            "ordinary_path_eligible": all(record["ordinary_path_eligible"] for record in records),
            "specialist_flags": flags,
            "required_actions": sorted({action for record in records for action in record["audit"]["required_actions"]}),
        })
    if problems:
        status = "electronic_state_identity_unresolved"
    elif any(record["status"] == "electronic_state_identity_unresolved" for record in channel_audits):
        status = "electronic_state_identity_unresolved"
    elif any(record["status"] == "spin_changing_rate_model_unresolved" for record in channel_audits):
        status = "spin_changing_rate_model_unresolved"
    elif any(not record["ordinary_path_eligible"] for record in channel_audits):
        status = "surface_resolution_required"
    else:
        status = "eligible_for_ordinary_path_review"
    result = {
        "schema_version": 1,
        "kind": "storca_retained_network_endpoint_state_audit",
        "policy_id": policy["policy_id"], "policy_hash": policy["scientific_hash"],
        "network_input_hash": network_input_hash,
        "input_network_scientific_hash": network.get("scientific_hash")
        if isinstance(network.get("scientific_hash"), str) else None,
        "status": status,
        "blocking_reasons": sorted(set(problems)),
        "channels": channel_audits,
        "reactions": reaction_audits,
        "ordinary_path_eligible": status == "eligible_for_ordinary_path_review",
        "source_validation_required": True,
        "central_rate_allowed": False,
        "validated_prediction_allowed": False,
        "lifetime_claim_allowed": False,
    }
    return {**result, "scientific_hash": _hash(result)}


def evaluate_method_applicability(
    route: Mapping[str, Any],
    *,
    condition_contract: Mapping[str, Any],
    method_profile: Mapping[str, Any],
    path_result: Mapping[str, Any] | None = None,
    rate_theory: str = "harmonic_gas_tst",
    tunneling_model: str | None = None,
    pressure_dependent: bool = False,
    diagnostics: Mapping[str, Any] | None = None,
    environment_assessment: Mapping[str, Any] | None = None,
    surface_assessment: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit conventional TST without replacing existing scientific gates.

    Supply the actual selected path and its environment/surface assessments.
    A ``supported`` diagnostic needs a retained source or evidence reference.
    ``warning`` remains conditional; ``inapplicable`` or ``escalation_required``
    withholds this model. Numerical execution/calibration is outside this API.
    """
    route = dict(route)
    conditions = dict(condition_contract)
    path = dict(path_result or {})
    verification = dict(path.get("route_verification") or {})
    retained_route = dict(path.get("route") or {})
    retained_environment = dict(retained_route.get("environment") or {})
    environment = dict(route.get("environment") or retained_environment)
    retained_context = dict(path.get("surface_resolution_scientific_context") or {})
    mapping = dict(verification.get("mapping") or path.get("mapping") or {})
    environment_gate = dict(environment_assessment or path.get("environment_applicability") or {})
    surface_gate = dict(surface_assessment or mapping.get("surface") or {})
    diagnostic_inputs = dict(diagnostics or route.get("method_diagnostics") or {})
    blockers: list[str] = []
    assumptions: list[str] = []
    warnings: list[str] = []
    escalation: list[str] = []

    route_model = _environment_model(
        environment.get("model") or environment.get("phase_model") or environment.get("phase")
    )
    condition_model = _environment_model(conditions.get("phase_model") or conditions.get("phase"))
    if route_model and condition_model and route_model != condition_model:
        blockers.append("route_and_condition_environment_mismatch")
    retained_model = _environment_model(
        retained_environment.get("model") or retained_environment.get("phase_model")
        or retained_environment.get("phase") or environment_gate.get("canonical_model")
    )
    if retained_model and (route_model or condition_model) and retained_model != (route_model or condition_model):
        blockers.append("retained_path_environment_mismatch")
    model = route_model or condition_model
    if model is None:
        model = "homogeneous_gas"
        assumptions.append("undeclared_environment_assumed_homogeneous_gas")
    if model not in {"homogeneous_gas", "implicit_solution"}:
        blockers.append("environment_outside_ordinary_homogeneous_tst_scope")
        escalation.append("environment_specific_path_and_rate_model")
    theory = str(rate_theory).strip().casefold()
    if theory not in {"harmonic_gas_tst", "harmonic_solution_eyring_tst", "solution_eyring_tst"}:
        blockers.append("rate_theory_outside_ordinary_tst_policy")
    elif (theory == "harmonic_gas_tst") != (model == "homogeneous_gas"):
        blockers.append("rate_theory_environment_mismatch")
    if environment_gate and environment_gate.get("status") != "applicable":
        blockers.append("existing_environment_gate_failed")
    for key in ("solvent", "solvent_model"):
        declared, requested = environment.get(key), conditions.get(key)
        if declared and requested and str(declared).casefold() != str(requested).casefold():
            blockers.append(f"route_and_condition_{key}_mismatch")
    if model == "implicit_solution":
        if not (environment.get("solvent") or conditions.get("solvent")):
            blockers.append("implicit_solution_solvent_missing")
        if str(environment.get("solvent_model") or conditions.get("solvent_model") or "").upper() not in {"CPCM", "SMD"}:
            blockers.append("implicit_solution_model_missing_or_unsupported")
        if not environment_gate:
            assumptions.append("orca_solvent_keyword_gate_evidence_unavailable")

    retained_keywords = (retained_context.get("electronic_structure") or {}).get("method_keywords")
    expected_keywords = method_profile.get("orca_keywords")
    if retained_keywords and expected_keywords:
        def electronic_keywords(values):
            return sorted(str(value).strip().casefold() for value in values
                          if not re.fullmatch(r"(?:cpcm|smd)\s*\([^)]*\)", str(value), re.IGNORECASE))
        if electronic_keywords(retained_keywords) != electronic_keywords(expected_keywords):
            blockers.append("retained_path_electronic_structure_method_mismatch")
    else:
        assumptions.append("retained_path_electronic_structure_method_evidence_unavailable")

    flags = {**dict(environment.get("applicability_flags") or {}),
             **dict(route.get("applicability_flags") or {})}
    for flag, required_method in _SPECIALIST_FLAGS.items():
        if flags.get(flag) is True or route.get(flag) is True or environment.get(flag) is True:
            blockers.append("declared_" + flag)
            escalation.append(required_method)
    light = str(environment.get("light_condition") or conditions.get("light_condition") or "dark").casefold()
    if light != "dark":
        assumptions.append("thermal_rate_is_only_one_component_of_illuminated_fate")
        escalation.append("complementary_photochemical_evidence_for_fate_claim")

    if surface_gate and not surface_gate.get("valid"):
        blockers.append("existing_surface_gate_failed")
        escalation.append("electronic_surface_resolution")
    reactant_states = _side_states(route, verification, "reactants")
    product_states = _side_states(route, verification, "products")
    selected = verification.get("selected_multiplicity", route.get("multiplicity"))
    declared_surfaces = (
        route.get("multiplicity"), verification.get("selected_multiplicity"),
        surface_gate.get("selected_multiplicity"),
        (retained_context.get("electronic_structure") or {}).get("multiplicity"),
    )
    try:
        positive_selections = {int(value) for value in declared_surfaces
                               if value is not None and int(value) > 0}
    except (TypeError, ValueError):
        positive_selections = set()
        blockers.append("selected_surface_declaration_malformed")
    if len(positive_selections) > 1:
        blockers.append("selected_surface_declarations_disagree")
    selected = surface_gate.get("selected_multiplicity", selected)
    surface_record: dict[str, Any] = {"existing_assessment": surface_gate or None}
    if reactant_states and product_states and selected is not None:
        from .route_geometry import _coupled_multiplicities

        common = sorted(set(_coupled_multiplicities(reactant_states)) & set(_coupled_multiplicities(product_states)))
        try:
            selected = int(selected)
        except (TypeError, ValueError):
            selected = 0
        surface_record.update(
            selected_multiplicity=selected, common_multiplicities=common,
            reactant_fragment_multiplicities=list(reactant_states),
            product_fragment_multiplicities=list(product_states),
        )
        if selected not in common:
            blockers.append("selected_surface_not_common_to_endpoints")
            escalation.append("spin_crossing_or_electronic_surface_resolution")
        elif len(common) > 1:
            assumptions.append("only_selected_common_spin_surface_is_represented")
    elif surface_gate.get("valid") and surface_gate.get("selected_multiplicity"):
        surface_record["status"] = "reused_existing_selected_common_surface"
    else:
        blockers.append("selected_common_surface_evidence_missing")

    path_status = verification.get("status")
    classification = verification.get("path_classification")
    if path_status != "verified_transition_state_path":
        blockers.append("verified_transition_state_path_required")
    if classification != "verified_barriered_path":
        blockers.append("verified_barriered_path_required")
    if str(classification or "").startswith("barrierless") or path_status in {
        "barrierless_rate_model_required", "dissociation_rate_model_required",
    }:
        escalation.append("capture_variational_or_dissociation_rate_model")

    reactant_count = len((verification.get("arkane_inputs") or {}).get("reactants") or ())
    product_count = len((verification.get("arkane_inputs") or {}).get("products") or ())
    reactant_count = reactant_count or len(route.get("reactant_smiles") or ()) or len(reactant_states)
    product_count = product_count or len(route.get("product_smiles") or ()) or len(product_states)
    if model == "homogeneous_gas" and {reactant_count, product_count} == {1, 2} and not pressure_dependent:
        blockers.append("well_channel_pressure_model_required")
        escalation.append("condition_matched_master_equation_model")
    if flags.get("collisional_stabilization_required") is True and not pressure_dependent:
        blockers.append("declared_collisional_stabilization_missing_rate_model")
        escalation.append("condition_matched_master_equation_model")
    if model == "homogeneous_gas" and not pressure_dependent:
        assumptions.append("high_pressure_limit_rate_does_not_establish_requested_pressure_adequacy")
    if model == "implicit_solution" and pressure_dependent:
        blockers.append("gas_master_equation_not_solution_rate_model")

    tunneling = str(tunneling_model or "none").strip().casefold() or "none"
    if tunneling not in {"none", "wigner", "eckart"}:
        blockers.append("unsupported_conventional_tunneling_model")
    if tunneling == "none":
        assumptions.append("transmission_coefficient_one_without_tunneling_correction")
    else:
        assumptions.append("declared_one_dimensional_tunneling_approximation")
    hydrogen_transfers = mapping.get("hydrogen_transfers")
    hydrogen_record = {
        "status": ("detected" if hydrogen_transfers else "not_detected_by_atom_mapping")
        if "hydrogen_transfers" in mapping else "not_available",
        "evidence": hydrogen_transfers,
        "interpretation": "Atom correspondence detects H transfer; it does not measure tunneling importance.",
    }
    if hydrogen_transfers:
        warnings.append("hydrogen_transfer_requires_tunneling_treatment_assessment")

    diagnostic_records: dict[str, Any] = {}
    for category in _DIAGNOSTICS:
        raw = diagnostic_inputs.get(category)
        if raw is not None and not isinstance(raw, Mapping):
            blockers.append(category + "_diagnostic_malformed")
            diagnostic_records[category] = {"status": "invalid"}
            continue
        record = dict(raw or {})
        status = str(record.get("status") or "not_available")
        if status not in {
            "not_available", "not_evaluated", "conditional_assumption", "supported",
            "warning", "inapplicable", "escalation_required",
        }:
            blockers.append(category + "_diagnostic_status_unrecognized")
        elif status in {"inapplicable", "escalation_required"}:
            blockers.append(category + "_diagnostic_requires_escalation")
            escalation.append(str(record.get("required_method") or category + "_model_assessment"))
        elif status == "warning":
            warnings.append(category + "_diagnostic_warning")
            assumptions.append(category + "_adequacy_unestablished")
        elif status != "supported":
            assumptions.append(category + "_adequacy_unestablished")
        elif not (record.get("source") or record.get("evidence_reference") or record.get("evidence_artifacts")):
            status = "not_available"
            assumptions.append(category + "_supported_diagnostic_source_missing")
        diagnostic_records[category] = {**record, "status": status}

    policy = method_applicability_policy()
    missing_path_or_state = {
        "verified_transition_state_path_required", "verified_barriered_path_required",
        "selected_common_surface_evidence_missing",
    }
    status = (
        "path_or_state_evidence_required" if blockers and set(blockers) <= missing_path_or_state and not escalation
        else "method_escalation_required" if blockers else "conditional_model"
    )
    result = {
        "schema_version": 1,
        "kind": "storca_method_applicability_decision",
        "policy_id": policy["policy_id"], "policy_hash": policy["scientific_hash"],
        "status": status,
        "route_id": route.get("route_id"),
        "scope": policy["scope"],
        "conditions": conditions,
        "method_profile": dict(method_profile),
        "rate_theory": rate_theory, "tunneling_model": tunneling,
        "pressure_dependent": bool(pressure_dependent),
        "environment": {"route_model": route_model, "condition_model": condition_model,
                        "retained_path_model": retained_model,
                        "selected_model": model, "existing_assessment": environment_gate or None},
        "surface": surface_record,
        "path": {"status": path_status, "classification": classification},
        "hydrogen_transfer": hydrogen_record,
        "diagnostics": diagnostic_records,
        "blocking_reasons": sorted(set(blockers)),
        "conditional_assumptions": sorted(set(assumptions)),
        "warnings": sorted(set(warnings)), "required_actions": sorted(set(escalation)),
        "central_rate_allowed": not blockers,
        "method_accuracy_validated": False,
        "validated_prediction_allowed": False, "lifetime_claim_allowed": False,
        "rate_uncertainty": {"status": "calibration_required",
                             "lower_rate_multiplier": None, "upper_rate_multiplier": None},
        "sources": dict(_SOURCES),
    }
    return {**result, "scientific_hash": _hash(result)}


__all__ = [
    "METHOD_APPLICABILITY_POLICY_ID", "method_applicability_policy", "evaluate_method_applicability",
    "ENDPOINT_STATE_AUDIT_POLICY_ID", "endpoint_state_audit_policy",
    "audit_endpoint_electronic_states", "audit_retained_network_states",
]
