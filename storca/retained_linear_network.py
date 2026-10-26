"""Source-bound, snapshot-local first-order networks, not chemical universals.

Full native identifiers distinguish structures/electronic states in one exact
export. They are deliberately not persistent IDs across regenerated mechanisms.
This narrow diagnostic domain does not validate estimated rates or establish
physical network completeness, pressure applicability or a lifetime.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path

from .reverse_relevance import _conditions


RETAINED_NETWORK_ALGORITHM = "source_bound_retained_linear_network_v2"
RETAINED_NETWORK_KIND = "source_bound_retained_linear_network"


def artifact(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def reopen(record):
    if not isinstance(record, dict) or artifact(record["path"]) != record:
        raise ValueError("retained_network_source_missing_or_mutated")
    return Path(record["path"])


def scientific_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def _number(value, field, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise ValueError("retained_network_" + field + "_invalid")
    value = float(value)
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError("retained_network_" + field + "_invalid")
    return value


def _sources(records):
    unique = {}
    for record in records:
        reopen(record)
        if record["path"] in unique and unique[record["path"]] != record:
            raise ValueError("retained_network_conflicting_source_artifacts")
        unique[record["path"]] = record
    return sorted(unique.values(), key=lambda item: item["path"])


def _scope(evidence, contract):
    _conditions(contract)
    def normalized(value):
        return " ".join(str(value).strip().casefold().replace("_", " ").replace("-", " ").split())
    supported = {"gas", "gas phase", "homogeneous gas", "homogeneous gas phase",
                 "homogeneous gas phase pure target model", "homogeneous gas phase ideal gas",
                 "homogeneous gas phase surrogate", "homogeneous gas phase low pressure surrogate"}
    # Inspect all declarations, not the first truthy phase alias. A gas-phase
    # label cannot license a simultaneously declared surface/solution model.
    aliases = [contract[key] for key in ("phase_model", "phase", "model", "environment_kind", "environment_model")
               if key in contract and contract[key] not in (None, "")]
    specialized = ("solvent", "solvent_model", "explicit_solvent", "microsolvation", "qmmm", "qm_mm",
                   "surface", "surface_model", "electrode", "electrochemical_model", "electrode_potential",
                   "interface", "interface_model", "bulk_model", "solid_model", "liquid_model",
                   "environment", "irradiance", "photon_flux", "light_exposure", "light_spectrum",
                   "temperature_profile", "pressure_profile", "transport_model")
    if not aliases or any(normalized(value) not in supported for value in aliases) \
            or contract.get("light_condition") != "dark" \
            or contract.get("light_model") is not None \
            or any(contract.get(key) is not None and contract.get(key) != "" and contract.get(key) is not False
                   for key in specialized) \
            or any(key in contract and contract[key] is not True for key in ("closed", "isothermal")) \
            or ("open_system" in contract and contract["open_system"] is not False) \
            or contract.get("relative_humidity") is not None \
            or contract.get("radical_sources") not in (None, [], ["none"]) \
            or contract.get("target_mole_fraction", 1.0) != 1.0 \
            or contract.get("observable") not in (None, "parent_retention"):
        raise ValueError("retained_network_environment_outside_scope")
    if contract.get("composition") != {contract["target_label"]: 1.0}:
        raise ValueError("retained_network_requires_pure_parent_inventory")
    reactor = (evidence.get("independent_cantera_propagation") or {}).get("reactor_model") or {}
    if reactor.get("policy_id") != "isothermal-ideal-gas-constant-pressure-v1" \
            or reactor.get("energy_enabled") is not False or reactor.get("condition_validation") != "passed" \
            or any(reactor.get(key) != contract[key] for key in ("temperature_K", "pressure_bar")):
        raise ValueError("retained_network_isothermal_reactor_evidence_unresolved")


def _composition(value, field, *, allow_empty=False):
    if not isinstance(value, dict) or (not value and not allow_empty) \
            or any(not isinstance(key, str) or not key or type(count) is not int or count <= 0
                   for key, count in value.items()):
        raise ValueError("retained_network_" + field + "_invalid")
    return value


def _record_list(value, field):
    if not isinstance(value, list) or not value or any(not isinstance(item, dict) for item in value):
        raise ValueError("retained_network_" + field + "_invalid")
    return value


def build_retained_linear_network(model_file, replacement, condition_contract, rate_snapshot_file):
    """Reopen native evidence and admit only its complete first-order inventory.

The snapshot must be produced read-only by the native ``linear_model_snapshot``
operation. No label normalization, legacy identity fallback, source mutation or
simulation is performed here. Reversible channels use Cantera operative rates,
not an inspected fitted reverse Arrhenius approximation.
    """
    from .flux_verification import verify_applied_replacement
    from .reaction_correspondence import resolve_replacement_correspondence

    model_ref, snapshot_ref = artifact(model_file), artifact(rate_snapshot_file)
    document = json.loads(Path(model_file).read_text())
    if document.get("condition_contract") != condition_contract:
        raise ValueError("retained_network_exact_condition_contract_mismatch")
    evidence = document["rmg_evidence"]
    _scope(evidence, condition_contract)
    inspection = evidence.get("mechanism_inspection") or {}
    if inspection.get("status") != "completed":
        raise ValueError("retained_network_native_inspection_incomplete")
    identity = resolve_replacement_correspondence(evidence, replacement)
    if identity.get("status") != "matched":
        raise ValueError("retained_network_modern_correspondence_required: "
                         + ",".join(identity.get("blocking_reasons") or [identity.get("status", "missing")]))
    applied = verify_applied_replacement(evidence, replacement, condition_contract)
    if applied.get("status") != "passed":
        raise ValueError("retained_network_applied_replacement_unresolved: "
                         + ",".join(applied.get("blocking_reasons") or []))
    maps = identity["correspondence"]["species_correspondence"]
    if any(len(maps.get(side) or []) != 1 for side in ("reactants", "products")):
        raise ValueError("retained_network_replacement_must_be_one_to_one")
    parent, product = (maps[side][0]["model_label"] for side in ("reactants", "products"))
    if maps["reactants"][0]["source_label"] != condition_contract["target_label"]:
        raise ValueError("retained_network_replacement_is_not_declared_parent_loss")

    snapshot = json.loads(Path(rate_snapshot_file).read_text())
    condition = {key: condition_contract[key] for key in ("temperature_K", "pressure_bar")}
    if snapshot.get("status") != "completed" \
            or snapshot.get("kind") != "native_fixed_first_order_rate_snapshot" \
            or snapshot.get("rate_basis") != "cantera_operative_rate_constants" \
            or snapshot.get("condition") != condition:
        raise ValueError("retained_network_native_rate_snapshot_unresolved")
    snapshot_sources = _sources(snapshot.get("source_artifacts") or [])
    exports = evidence.get("artifacts") or {}
    required = [artifact(exports[field]) for field in ("chemkin", "species_dictionary")]
    if exports.get("transport"):
        required.append(artifact(exports["transport"]))
    if any(record not in snapshot_sources for record in required):
        raise ValueError("retained_network_snapshot_exports_mismatch")
    species = _record_list(snapshot.get("species"), "species")
    ids = [item.get("chemkin_identifier") for item in species]
    inspected = _record_list(inspection.get("species"), "inspected_species")
    if not ids or any(not isinstance(item, str) or not item for item in ids) \
            or len(set(ids)) != len(ids) \
            or set(ids) != {item.get("chemkin_identifier") for item in inspected} \
            or len(ids) != len(inspected) or parent == product or not {parent, product}.issubset(ids):
        raise ValueError("retained_network_species_inventory_invalid")
    by_id = {item["chemkin_identifier"]: item for item in species}
    inspected_by_id = {item["chemkin_identifier"]: item for item in inspected}
    for native_id, item in by_id.items():
        if type(item.get("multiplicity")) is not int or item["multiplicity"] < 1 \
                or type(item.get("net_charge")) is not int \
                or not isinstance(item.get("adjacency_list"), str) or not item["adjacency_list"].strip():
            raise ValueError("retained_network_species_state_or_graph_missing")
        _composition(item.get("atom_composition"), "atom_composition")
        _composition(item.get("isotopic_composition"), "isotopic_composition", allow_empty=True)
        for key in ("multiplicity", "net_charge"):
            if inspected_by_id[native_id].get(key) != item[key]:
                raise ValueError("retained_network_species_state_changed")
    initial = {item: float(item == parent) for item in ids}
    if not isinstance(snapshot.get("initial_fractions"), dict):
        raise ValueError("retained_network_native_initial_inventory_mismatch")
    for value in snapshot["initial_fractions"].values():
        _number(value, "initial_fraction")
    if snapshot.get("parent_id") != parent or snapshot.get("initial_fractions") != initial:
        raise ValueError("retained_network_native_initial_inventory_mismatch")

    channels = _record_list(snapshot.get("channels"), "channels")
    expected = {}
    for position, reaction in enumerate(inspection.get("reactions") or []):
        left, right = reaction.get("reactant_labels") or [], reaction.get("product_labels") or []
        if len(left) != 1 or len(right) != 1 or reaction.get("molecularity") != 1 \
                or (reaction.get("reversible") is True and reaction.get("reverse_molecularity") != 1):
            raise ValueError("retained_network_nonlinear_or_nonconservative_model")
        if reaction.get("evaluated_at") != condition:
            raise ValueError("retained_network_native_inspection_condition_mismatch")
        for direction in ("forward", "reverse") if reaction.get("reversible") is True else ("forward",):
            a, b = (left[0], right[0]) if direction == "forward" else (right[0], left[0])
            expected[f"native-{position}-{direction}"] = (reaction.get("reaction_index"), direction, a, b)
    if not channels or len(channels) != len(expected) \
            or len({item.get("channel_id") for item in channels}) != len(channels):
        raise ValueError("retained_network_complete_directed_inventory_required")
    for item in channels:
        signature = expected.get(item.get("channel_id"))
        if signature is None or type(item.get("reaction_index")) is not int \
                or (item.get("reaction_index"), item.get("direction"),
                next(iter(item.get("reactants") or {}), None), next(iter(item.get("products") or {}), None)) != signature \
                or item.get("reactants") != {signature[2]: 1} or item.get("products") != {signature[3]: 1} \
                or any(type(value) is not int or value != 1 for side in ("reactants", "products")
                       for value in (item.get(side) or {}).values()) \
                or item.get("molecularity") != 1 or isinstance(item.get("molecularity"), bool) \
                or item.get("rate_units", "s^-1") != "s^-1":
            raise ValueError("retained_network_channel_mapping_or_units_invalid")
        a, b = signature[2:]
        if a == b or a not in by_id or b not in by_id:
            raise ValueError("retained_network_channel_species_invalid")
        _number(item.get("rate_s^-1"), "operative_rate")
        if any(by_id[a][key] != by_id[b][key]
               for key in ("atom_composition", "isotopic_composition", "net_charge")):
            raise ValueError("retained_network_atom_isotope_or_charge_unbalanced")
    candidate = identity["candidate"]
    if candidate.get("reversible") is not False or identity["applied_direction"] != "forward":
        raise ValueError("retained_network_replacement_must_remain_irreversible")
    position = (inspection.get("reactions") or []).index(candidate)
    forward_id = f"native-{position}-forward"
    selected = next(item for item in channels if item["channel_id"] == forward_id)
    expected_rate = _number((replacement.get("condition_rate") or {}).get("value"), "replacement_rate", positive=True)
    if selected["reactants"] != {parent: 1} or selected["products"] != {product: 1} \
            or not math.isclose(selected["rate_s^-1"], expected_rate, rel_tol=1e-4, abs_tol=0):
        raise ValueError("retained_network_operative_replacement_rate_mismatch")
    if any(item["reactants"] == {product: 1} and item["products"] == {parent: 1} for item in channels):
        raise ValueError("retained_network_omitted_direct_return_reintroduced")
    library = Path(replacement["library"])
    sources = _sources([model_ref, snapshot_ref, *snapshot_sources,
                        inspection["replacement_identity_artifact"],
                        *(artifact(library / name) for name in ("reactions.py", "dictionary.txt", "generated-kinetics.json"))])
    result = {"schema_version": 2, "kind": RETAINED_NETWORK_KIND,
              "algorithm": RETAINED_NETWORK_ALGORITHM, "closed": True, "isothermal": True,
              "condition_contract": copy.deepcopy(condition_contract), "condition": condition,
              "species_ids": ids, "species": copy.deepcopy(species), "channels": copy.deepcopy(channels),
              "initial_fractions": initial, "parent_id": parent, "product_id": product,
              "replacement_forward_channel_id": forward_id,
              "replacement_route_id": replacement["route_id"], "applied_replacement": applied,
              "operative_forward_relative_difference_from_retained_fit":
                  abs(selected["rate_s^-1"] - expected_rate) / max(selected["rate_s^-1"], expected_rate),
              "rate_serialization_difference_not_bounded_by_reverse_omission": True,
              "rate_basis": snapshot["rate_basis"], "source_artifacts": sources,
              "input_bindings": {"model_file": model_ref, "rate_snapshot_file": snapshot_ref,
                                 "replacement": copy.deepcopy(replacement)},
              "limitations": ["Snapshot-local native IDs, not persistent chemical identities.",
                               "Declared fixed rates only; estimated electronic-state-changing channels remain unverified.",
                               "No physical uncertainty, mechanism completeness, pressure or lifetime claim."],
              "lifetime_claim_allowed": False, "conclusive_observable_allowed": False}
    result["scientific_hash"] = scientific_hash(result)
    return result


def validate_retained_linear_network(network, condition_contract=None):
    """Reopen and recompute; accepting an opaque hash is not sufficient."""
    if not isinstance(network, dict) or network.get("kind") != RETAINED_NETWORK_KIND \
            or network.get("algorithm") != RETAINED_NETWORK_ALGORITHM:
        raise ValueError("retained_network_policy_unsupported")
    contract = network.get("condition_contract") if condition_contract is None else condition_contract
    if contract != network.get("condition_contract"):
        raise ValueError("retained_network_exact_condition_contract_mismatch")
    bindings = network["input_bindings"]
    model = reopen(bindings["model_file"])
    snapshot = reopen(bindings["rate_snapshot_file"])
    replay = build_retained_linear_network(model, bindings["replacement"], contract, snapshot)
    if replay != network:
        raise ValueError("retained_network_not_supported_by_retained_sources")
    return replay
