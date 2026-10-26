"""Trace inherited rates without promoting provenance to physical validity.

This bounded audit consumes a retained network and read-only source snapshots.
It does not construct, adjust, remove or calibrate rates. Source documentation
of spin treatment is distinct from its applicability to the requested system.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

from .retained_linear_network import reopen, scientific_hash


POLICY_ID = "retained-rate-source-applicability-review-v1"
AUDIT_KIND = "retained_network_rate_applicability_audit"
CLAIMS = {
    "physical_rates_validated": False,
    "physical_rate_bounds_established": False,
    "uncertainty_calibrated": False,
    "rate_correction_allowed": False,
    "branch_removal_allowed": False,
    "ordinary_path_override_allowed": False,
    "lifetime_claim_allowed": False,
    "conclusive_observable_allowed": False,
}


def _source_spin_documentation(source):
    """Explicit textual source triage, not a nonadiabatic rate certificate.

Require both a spin-crossing description and a lifetime-to-rate convention
in the extracted family description. Citation presence, rank, and absence of
an activation barrier cannot by themselves meet this descriptive criterion.
    """
    matches = []
    for document in source.get("source_documents", []):
        parsed = document.get("parsed") or {}
        metadata = parsed.get("metadata") or {}
        text = metadata.get("longDesc")
        if not isinstance(text, str):
            continue
        normalized = " ".join(text.casefold().split())
        compact = normalized.replace(" ", "")
        if "intersystem crossing" in normalized and "k=1/tau" in compact:
            matches.append({"role": document.get("role"), "artifact": document["artifact"],
                            "documentation": text,
                            "criterion": "explicit_intersystem_crossing_and_k_inverse_lifetime"})
    return {
        "status": ("empirical_spin_crossing_basis_documented_in_source" if matches
                   else "spin_treatment_not_established_by_extracted_sources"),
        "documented_basis": matches,
        "exact_system_applicability_established": False,
        "original_primary_reference_audited": False,
        "state_specific_rate_validated": False,
        "matched_rule_correlation_association_validated": False,
        "interpretation": "Source-text classification only; matching a documented family does not validate its rate.",
    }


def _native_sources(network):
    snapshot = json.loads(reopen(network["input_bindings"]["rate_snapshot_file"]).read_text())
    native = snapshot.get("rate_source_snapshot")
    if native is None:
        return None
    from .native_rate_sources import validate_rate_source_snapshot
    native = validate_rate_source_snapshot(native)
    if native.get("condition") != network["condition"] \
            or native.get("condition_contract_id") != network["condition_contract"]["contract_id"]:
        raise ValueError("retained_rate_audit_native_source_condition_mismatch")
    # Every model export used for source attribution must be one of the exact
    # model exports already admitted by the operative-rate snapshot.
    exports = native.get("input_bindings", {}).get("model_exports")
    if not isinstance(exports, dict):
        raise ValueError("retained_rate_audit_native_source_model_bindings_missing")
    records = list(exports.values())
    if any(not isinstance(item, dict) for item in records) \
            or sorted(records, key=lambda item: item["path"]) != sorted(
                snapshot["source_artifacts"], key=lambda item: item["path"]):
        raise ValueError("retained_rate_audit_native_source_model_mismatch")
    identities = [(item["chemkin_identifier"], item["multiplicity"]) for item in native.get("species", [])]
    expected = [(item["chemkin_identifier"], item["multiplicity"]) for item in network["species"]]
    if sorted(identities) != sorted(expected):
        raise ValueError("retained_rate_audit_native_source_electronic_state_mismatch")
    return native


def audit_retained_network_rates(network, state_audit):
    """Reopen evidence and distinguish source facts from unresolved physics."""
    from .retained_rate_recipe import audit_retained_arkane_recipe
    from .retained_linear_network import validate_retained_linear_network
    from .method_applicability import audit_retained_network_states

    network = validate_retained_linear_network(network)
    if state_audit != audit_retained_network_states(network):
        raise ValueError("retained_rate_audit_state_network_mismatch")
    document = json.loads(reopen(network["input_bindings"]["model_file"]).read_text())
    inspection = document["rmg_evidence"]["mechanism_inspection"]
    reactions = inspection["reactions"]
    native = _native_sources(network)
    sources = {} if native is None else {item["reaction_index"]: item
                                         for item in native.get("reactions", [])}
    if native is not None and len(sources) != len(native.get("reactions", [])):
        raise ValueError("retained_rate_audit_duplicate_native_sources")
    if native is not None and set(sources) != {item["reaction_index"] for item in reactions}:
        raise ValueError("retained_rate_audit_complete_native_source_inventory_required")
    audits = {item["channel_id"]: item for item in state_audit["channels"]}
    species = {item["chemkin_identifier"]: item for item in network["species"]}
    result = {"schema_version": 1, "kind": AUDIT_KIND, "policy_id": POLICY_ID,
              "status": "physical_rate_applicability_unresolved",
              "condition_contract": copy.deepcopy(network["condition_contract"]),
              "network_scientific_hash": network["scientific_hash"],
              "native_source_snapshot": native, "channels": [],
              "source_snapshot_status": "not_retained" if native is None else native.get("status"),
              "limitations": [
                  "Database annotations are provenance, not experimental or theoretical accuracy certificates.",
                  "A documented spin-crossing correlation does not establish exact-system state or pressure applicability.",
                  "A recycled training rate is not an independent experimental observation.",
                  "Nominal flux or a removal scenario cannot justify deleting an unresolved reaction.",
                  "No rate, electronic state, condition or original artifact is altered by this audit."],
              **CLAIMS}
    replacement_id = network["replacement_forward_channel_id"]
    for channel in network["channels"]:
        index = channel["reaction_index"]
        if type(index) is not int or not 0 <= index < len(reactions) \
                or reactions[index].get("reaction_index") != index:
            raise ValueError("retained_rate_audit_native_reaction_index_mismatch")
        inspected = reactions[index]
        source = sources.get(index)
        left, right = (list(channel["reactants"]), list(channel["products"]))
        forward_left, forward_right = (left, right) if channel["direction"] == "forward" else (right, left)
        if inspected["reactant_labels"] != forward_left or inspected["product_labels"] != forward_right:
            raise ValueError("retained_rate_audit_inspected_channel_mismatch")
        if source is not None and (source.get("reactants") != forward_left
                                  or source.get("products") != forward_right
                                  or source.get("reversible") != inspected.get("reversible")):
            raise ValueError("retained_rate_audit_source_channel_mismatch")
        spin = _source_spin_documentation(source or {})
        ordinary = audits[channel["channel_id"]]["ordinary_path_eligible"]
        if channel["channel_id"] == replacement_id:
            category = "retained_arkane_central_model_rate"
        elif not ordinary and spin["documented_basis"]:
            category = "documented_spin_aware_empirical_source_applicability_unresolved"
        elif not ordinary:
            category = "state_changing_estimate_spin_treatment_unestablished"
        else:
            category = "inherited_rate_applicability_unresolved"
        result["channels"].append({
            "channel_id": channel["channel_id"], "reaction_index": index,
            "direction": channel["direction"], "reactants": copy.deepcopy(channel["reactants"]),
            "products": copy.deepcopy(channel["products"]),
            "endpoint_multiplicities": {item: species[item]["multiplicity"] for item in left + right},
            "endpoint_state_basis": "exact_retained_native_graph_state_inventory",
            "operative_rate_s^-1": channel["rate_s^-1"], "rate_category": category,
            "ordinary_path_eligible": ordinary, "physical_rate_applicability": "unresolved",
            "rate_origin_annotation": {key: copy.deepcopy(inspected.get(key)) for key in
                                       ("source_library", "kinetics_comment")},
            "source_ancestry": copy.deepcopy(source), "source_spin_documentation": spin,
            "owning_storca_library_ancestry_reopened": channel["channel_id"] == replacement_id,
            "reverse_rate_basis": ("retained_native_thermochemistry_and_operative_forward_rate"
                                   if channel["direction"] == "reverse" else None),
            "independent_measurement_established": False,
            "physical_rate_validated": False, "automatic_rate_change_allowed": False,
        })
    replacement = network["input_bindings"]["replacement"]
    result["generated_library_ancestry"] = {
        "status": "existing_source_bound_insertion_validation_reopened",
        "channel_id": replacement_id,
        "native_source_snapshot_library_lookup_not_required": True,
        "evidence": copy.deepcopy(network["applied_replacement"]),
        "physical_rate_validated": False,
    }
    result["dominant_rate_recipe_audit"] = audit_retained_arkane_recipe(replacement, network["condition_contract"])
    result["thermochemistry_sources"] = copy.deepcopy([] if native is None else native.get("species", []))
    result["requested_pressure_review"] = {
        "status": "requested_pressure_applicability_unresolved",
        "temperature_K": network["condition"]["temperature_K"],
        "pressure_bar": network["condition"]["pressure_bar"],
        "quantitative_high_pressure_limit_error_bound": None,
        "reason": "A recorded pressure, pressure grid, or local experimental plateau is not a quantified high-pressure-limit error bound.",
        "experimental_reference_used_to_adjust_rate": False,
        "pressure_model_correction_allowed": False,
    }
    records = [network["input_bindings"]["model_file"], network["input_bindings"]["rate_snapshot_file"]]
    if native is not None:
        records.extend(native.get("source_artifacts", []))
    records.extend(result["dominant_rate_recipe_audit"].get("source_artifacts", []))
    unique = {}
    for original in records:
        record = {key: original[key] for key in ("path", "sha256")}
        reopen(record)
        if record["path"] in unique and unique[record["path"]] != record:
            raise ValueError("retained_rate_audit_conflicting_source_hashes")
        unique[record["path"]] = record
    result["source_artifacts"] = sorted(unique.values(), key=lambda item: item["path"])
    result["scientific_hash"] = scientific_hash(result)
    return result


def validate_retained_network_rate_audit(audit, network, state_audit):
    if not isinstance(audit, dict) or audit.get("kind") != AUDIT_KIND \
            or any(audit.get(key) is not value for key, value in CLAIMS.items()):
        raise ValueError("retained_rate_audit_claim_or_kind_mismatch")
    replay = audit_retained_network_rates(network, state_audit)
    if replay != audit:
        raise ValueError("retained_rate_audit_not_supported_by_retained_sources")
    return replay
