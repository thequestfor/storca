"""Auditable finite-horizon omission bounds for a narrowly declared model.

This module does not estimate missing chemistry, rate uncertainty, or pressure
applicability.  It bounds the effect of omitting declared first-order product
return channels relative to a supplied, closed isothermal kinetic model.  In
particular, a small fitted reverse rate alone is not an upper-bound certificate.

For an inventory bound B_j relative to the initial parent inventory, the
integrated reverse source is at most t * sum_j(n_j * k_j_upper * B_j).  Variation
of constants and positivity/nonexpansiveness of a conservative linear network
make this a conservative absolute parent-retention error bound.  Each omitted
return injects opposite parent/product source terms; the difference between two
transition probabilities for the parent observable is at most one.  Nonlinear
chemistry cannot use this argument.  Alternative descriptions
of the *same* return channel use their maximum; distinct channels are summed.
The result is conditional on the declared kinetic model, never a confidence
interval for a physical lifetime or proof of conformer/network completeness.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping


REVERSE_RELEVANCE_POLICY_ID = "finite-horizon-reverse-return-v1"
MODEL_INVENTORY_ALGORITHM = "closed_first_order_linear_model_inventory_v1"
_TWO_SPECIES_INVENTORY_ALGORITHM = "closed_first_order_two_species_model_inventory_v1"
_SNAPSHOT_FIELDS = (
    "temperature_K", "pressure_bar", "target_duration_seconds",
    "retention_fraction", "target_label",
)
_UPPER_BOUND_METHODS = {
    "declared_fixed_model_envelope", "physical_bound", "interval_rate_bound",
}


def _canonical(value: Any) -> str:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as error:
        raise ValueError("Reverse relevance evidence must be finite JSON data") from error


def _copy(value: Any) -> Any:
    return json.loads(_canonical(value))


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be a mapping")
    return value


def _text(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a nonempty string")
    return value


def _number(value: Any, field: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise ValueError(f"{field} must be a finite number")
    number = float(value)
    if not math.isfinite(number) or number < 0 or (positive and number == 0):
        raise ValueError(f"{field} must be {'positive' if positive else 'nonnegative'} and finite")
    return number


def _flag(value: Any, field: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{field} must be a boolean")
    return value


def _list(value: Any, field: str, *, nonempty: bool = True) -> list:
    if not isinstance(value, list) or (nonempty and not value):
        raise ValueError(f"{field} must be {'a nonempty' if nonempty else 'a'} list")
    return value


def _assumptions(value: Any, field: str) -> None:
    for item in _list(value, field):
        _text(item, field)


def _sources(record: Mapping[str, Any], field: str) -> None:
    """Reopen immutable source records without launching or modifying anything."""
    for raw in _list(record.get("source_artifacts"), f"{field}.source_artifacts"):
        artifact = _mapping(raw, f"{field}.source_artifact")
        path = Path(_text(artifact.get("path"), f"{field}.source_artifact.path"))
        expected = _text(artifact.get("sha256"), f"{field}.source_artifact.sha256")
        if len(expected) != 64 or any(character not in "0123456789abcdef" for character in expected):
            raise ValueError(f"{field} source artifact requires a lowercase SHA-256")
        if not path.is_absolute() or not path.is_file():
            raise ValueError(f"{field} source artifact must be an existing absolute file: {path}")
        try:
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
        except OSError as error:
            raise ValueError(f"{field} source artifact is unreadable: {path}") from error
        if digest.hexdigest() != expected:
            raise ValueError(f"{field} source artifact hash changed: {path}")


def reverse_relevance_condition_snapshot(condition_contract: Mapping[str, Any]) -> dict:
    """Exact values required on every rate/inventory evidence record."""
    contract = _mapping(condition_contract, "condition_contract")
    try:
        return {field: contract[field] for field in _SNAPSHOT_FIELDS}
    except KeyError as error:
        raise ValueError(f"Condition contract is missing {error.args[0]}") from error


def _conditions(condition_contract: Mapping[str, Any]) -> tuple[dict, str]:
    contract = _mapping(condition_contract, "condition_contract")
    snapshot = reverse_relevance_condition_snapshot(contract)
    for field in ("temperature_K", "pressure_bar", "target_duration_seconds"):
        _number(snapshot[field], f"condition_contract.{field}", positive=True)
    retention = _number(snapshot["retention_fraction"], "condition_contract.retention_fraction", positive=True)
    if retention >= 1:
        raise ValueError("condition_contract.retention_fraction must be in (0, 1)")
    _text(snapshot["target_label"], "condition_contract.target_label")
    condition_id = _text(contract.get("contract_id"), "condition_contract.contract_id")
    payload = {key: value for key, value in contract.items() if key not in {"schema_version", "contract_id"}}
    if condition_id != "sha256:" + _digest(payload):
        raise ValueError("condition_contract hash does not match its exact scientific payload")
    return snapshot, condition_id


def _bind(record: Mapping[str, Any], snapshot: dict, condition_id: str, field: str) -> None:
    if record.get("condition_contract_id") != condition_id:
        raise ValueError(f"{field} condition_contract_id mismatch")
    if record.get("condition_snapshot") != snapshot:
        raise ValueError(f"{field} condition_snapshot mismatch")
    _sources(record, field)


def _inventory_model(envelope: Mapping[str, Any], *, parent: str, product: str,
                     route_id: str, forward_value: float, channel_uppers: dict[str, float]) -> bool:
    """Validate an algorithmic inventory, not an opaque completeness boolean."""
    model = _mapping(envelope.get("model_inventory"), "reverse_envelope.model_inventory")
    two_species = model.get("kind") == "closed_first_order_two_species_model_inventory"
    if model.get("kind") not in {
        "closed_first_order_two_species_model_inventory", "closed_first_order_linear_model_inventory",
    }:
        raise ValueError("Unsupported reverse-envelope model inventory kind")
    expected_algorithm = _TWO_SPECIES_INVENTORY_ALGORITHM if two_species else MODEL_INVENTORY_ALGORITHM
    if model.get("algorithm") != expected_algorithm:
        raise ValueError("Unsupported reverse-envelope model inventory algorithm")
    if model.get("isothermal") is not True or model.get("closed") is not True:
        raise ValueError("Model inventory must be closed and isothermal")
    species = _list(model.get("species_ids"), "model_inventory.species_ids")
    for species_id in species:
        _text(species_id, "model_inventory.species_id")
    if len(set(species)) != len(species) or not {parent, product}.issubset(species):
        raise ValueError("Model inventory requires distinct species including the parent and product")
    if two_species and set(species) != {parent, product}:
        raise ValueError("Two-species model inventory requires exactly parent and product")
    fractions = _mapping(model.get("initial_fractions"), "model_inventory.initial_fractions")
    if set(fractions) != set(species):
        raise ValueError("Model inventory initial fractions must name every species")
    if (_number(fractions[parent], "model_inventory.parent_fraction") != 1.0
            or any(_number(value, "model_inventory.initial_fraction") != 0.0
                   for species_id, value in fractions.items() if species_id != parent)):
        raise ValueError("Model inventory must start from a pure parent inventory")
    forward_ids, return_ids, seen = set(), set(), set()
    isolated_two_species = set(species) == {parent, product}
    for raw in _list(model.get("channels"), "model_inventory.channels"):
        channel = _mapping(raw, "model_inventory.channel")
        channel_id = _text(channel.get("channel_id"), "model_inventory.channel_id")
        if channel_id in seen:
            raise ValueError("Duplicate model inventory channel_id")
        seen.add(channel_id)
        if channel.get("molecularity") != 1 or isinstance(channel.get("molecularity"), bool):
            raise ValueError("Model inventory contains a non-first-order channel")
        model_rate = _number(channel.get("rate_s^-1"), "model_inventory.rate_s^-1")
        reactants = _mapping(channel.get("reactants"), "model_inventory.reactants")
        products = _mapping(channel.get("products"), "model_inventory.products")
        for stoichiometry in (reactants, products):
            if len(stoichiometry) != 1 or any(
                isinstance(value, bool) or _number(value, "model_inventory.stoichiometry") != 1.0
                for value in stoichiometry.values()
            ):
                raise ValueError("Model inventory channel must be one-to-one")
            if not set(stoichiometry).issubset(species):
                raise ValueError("Model inventory channel references an undeclared species")
        if reactants == products:
            raise ValueError("Model inventory contains a self-reaction")
        if reactants == {parent: 1} and products == {product: 1}:
            forward_ids.add(channel_id)
            if channel_id == route_id and model_rate != forward_value:
                raise ValueError("Forward rate and declared model inventory rate differ")
        elif reactants == {product: 1} and products == {parent: 1}:
            return_ids.add(channel_id)
            if channel_id in channel_uppers and model_rate > channel_uppers[channel_id]:
                raise ValueError("Reverse upper envelope is below the declared model inventory rate")
        else:
            isolated_two_species = False
    if forward_ids != {route_id} or return_ids != set(channel_uppers):
        raise ValueError("Model inventory and forward/reverse envelope channels differ")
    _sources(model, "model_inventory")
    return isolated_two_species


def _policy(value: Mapping[str, Any] | None, allowed_loss: float) -> dict:
    policy = {} if value is None else dict(_mapping(value, "policy"))
    if set(policy) - {"policy_id", "max_abs_retention_error"}:
        raise ValueError("Unknown reverse relevance policy field")
    policy_id = policy.get("policy_id", REVERSE_RELEVANCE_POLICY_ID)
    if policy_id != REVERSE_RELEVANCE_POLICY_ID:
        raise ValueError("Unsupported reverse relevance policy_id")
    tolerance = _number(policy.get("max_abs_retention_error", 1e-4),
                        "policy.max_abs_retention_error", positive=True)
    if tolerance > allowed_loss:
        raise ValueError("Absolute retention-error tolerance cannot exceed the allowed parent loss")
    return {"policy_id": policy_id, "max_abs_retention_error": tolerance}


def assess_reverse_return_relevance(
    condition_contract: Mapping[str, Any],
    forward_rate: Mapping[str, Any],
    reverse_envelope: Mapping[str, Any],
    inventory_bound: Mapping[str, Any],
    policy: Mapping[str, Any] | None = None,
) -> dict:
    """Assess a model-scoped omission error; malformed evidence raises ValueError.

    Unsupported but well-formed upper/inventory bounds yield diagnostic results
    with ``bounded_approximation_allowed=False``.  Neither a small reverse rate
    nor this assessment can promote a conditional rate to a conclusive lifetime.
    Inputs are retained verbatim as JSON-compatible evidence for replay.
    """
    inputs = _copy({"condition_contract": condition_contract, "forward_rate": forward_rate,
                    "reverse_envelope": reverse_envelope, "inventory_bound": inventory_bound})
    snapshot, condition_id = _conditions(inputs["condition_contract"])
    normalized_policy = _policy(policy, 1.0 - float(snapshot["retention_fraction"]))
    forward = _mapping(inputs["forward_rate"], "forward_rate")
    envelope = _mapping(inputs["reverse_envelope"], "reverse_envelope")
    inventory = _mapping(inputs["inventory_bound"], "inventory_bound")
    for name, record, kind in (
        ("forward_rate", forward, "verified_conditional_forward_rate"),
        ("reverse_envelope", envelope, "first_order_reverse_rate_upper_envelope"),
        ("inventory_bound", inventory, "finite_horizon_product_inventory_bound"),
    ):
        if record.get("kind") != kind:
            raise ValueError(f"{name} kind must be {kind}")
        _bind(record, snapshot, condition_id, name)
        _assumptions(record.get("model_assumptions"), f"{name}.model_assumptions")
    parent = _text(forward.get("parent_id"), "forward_rate.parent_id")
    product = _text(forward.get("product_id"), "forward_rate.product_id")
    route_id = _text(forward.get("route_id"), "forward_rate.route_id")
    if parent == product:
        raise ValueError("Parent and product must be distinct declared species")
    if forward.get("verification_status") != "verified_arkane_tst":
        raise ValueError("Forward rate must be a verified conditional Arkane rate")
    if forward.get("molecularity") != 1 or isinstance(forward.get("molecularity"), bool):
        raise ValueError("Only first-order forward rates are in scope")
    stoichiometry = _mapping(forward.get("stoichiometry"), "forward_rate.stoichiometry")
    if (set(stoichiometry) != {"parent_consumed", "product_formed"}
            or any(isinstance(value, bool) or _number(value, "forward_rate.stoichiometry") != 1.0
                   for value in stoichiometry.values())):
        raise ValueError("Forward channel must be one parent to one product")
    forward_value = _number(forward.get("value_s^-1"), "forward_rate.value_s^-1", positive=True)
    scope = _mapping(envelope.get("scope"), "reverse_envelope.scope")
    expected_scope = {"closed": True, "isothermal": True, "observable": "parent_retention",
                      "parent_id": parent, "forward_route_id": route_id,
                      "first_order_single_parent_product": True}
    if any(scope.get(key) != value or (isinstance(value, bool) and scope.get(key) is not value)
           for key, value in expected_scope.items()):
        raise ValueError("Reverse return relevance is outside its declared model/observable scope")
    upper_supported = _flag(envelope.get("upper_bound_supported"), "reverse_envelope.upper_bound_supported")
    inventory_supported = _flag(inventory.get("supported"), "inventory_bound.supported")
    if inventory.get("basis") != "initial_parent_inventory":
        raise ValueError("Product inventory must be normalized to initial parent inventory")
    inventory_by_product = {}
    for raw in _list(inventory.get("products"), "inventory_bound.products"):
        record = _mapping(raw, "inventory_bound.product")
        product_id = _text(record.get("product_id"), "inventory_bound.product_id")
        if product_id in inventory_by_product or product_id != product:
            raise ValueError("Inventory must name the single declared product exactly once")
        fraction = _number(record.get("maximum_fraction"), "inventory_bound.maximum_fraction")
        if fraction > 1.0:
            raise ValueError("Pure-parent one-to-one inventory fraction cannot exceed one")
        inventory_by_product[product_id] = fraction
    if inventory_supported and (
        inventory.get("bound_method") != "closed_one_to_one_conserved_inventory"
        or inventory_by_product[product] != 1.0
    ):
        raise ValueError("Supported inventory bound requires the conserved pure-parent ceiling of one")
    terms, channel_ids = [], set()
    all_alternatives_supported = True
    for raw in _list(envelope.get("channels"), "reverse_envelope.channels"):
        channel = _mapping(raw, "reverse_envelope.channel")
        channel_id = _text(channel.get("channel_id"), "reverse_envelope.channel_id")
        if channel_id in channel_ids:
            raise ValueError("Duplicate reverse envelope channel_id")
        channel_ids.add(channel_id)
        if channel.get("product_id") != product:
            raise ValueError("Reverse channel must return the declared product to parent")
        stoich = _number(channel.get("parent_stoichiometry"), "reverse_envelope.parent_stoichiometry", positive=True)
        if stoich != 1:
            raise ValueError("Only one-parent return channels are in scope")
        upper = _number(channel.get("upper_rate_s^-1"), "reverse_envelope.upper_rate_s^-1")
        alternatives, basis_ids = [], set()
        for raw_alternative in _list(channel.get("rate_basis_alternatives"), "reverse_envelope.rate_basis_alternatives"):
            alternative = _mapping(raw_alternative, "reverse_envelope.alternative")
            basis = _text(alternative.get("basis"), "reverse_envelope.alternative.basis")
            if basis in basis_ids:
                raise ValueError("Duplicate reverse rate basis alternative")
            basis_ids.add(basis)
            alternative_upper = _number(alternative.get("upper_rate_s^-1"), "reverse_envelope.alternative.upper_rate_s^-1")
            supported = _flag(alternative.get("upper_bound_supported"), "reverse_envelope.alternative.upper_bound_supported")
            method = alternative.get("upper_bound_method")
            if supported and method not in _UPPER_BOUND_METHODS:
                raise ValueError("Supported rate alternative requires an explicit upper-bound method")
            _text(alternative.get("bound_justification"), "reverse_envelope.alternative.bound_justification")
            _sources(alternative, "reverse_envelope.alternative")
            all_alternatives_supported = all_alternatives_supported and supported
            alternatives.append(alternative_upper)
        if not {"reached_conformer", "network_species"}.issubset(basis_ids):
            raise ValueError("Reverse envelope must include reached-conformer and network-species alternatives")
        if upper < max(alternatives):
            raise ValueError("Reverse upper envelope is below a declared basis alternative")
        coverage = _mapping(channel.get("coverage"), "reverse_envelope.coverage")
        if (coverage.get("reached_conformer") is not True or coverage.get("network_species") is not True
                or coverage.get("other_return_channels") != "excluded_by_model_inventory"):
            raise ValueError("Reverse envelope coverage must be tied to the explicit model inventory")
        source_upper = upper * stoich * inventory_by_product[product]
        if not math.isfinite(source_upper):
            raise ValueError("Reverse source upper bound overflows")
        terms.append({"channel_id": channel_id, "product_id": product,
                      "upper_rate_s^-1": upper, "parent_stoichiometry": stoich,
                      "product_inventory_fraction_bound": inventory_by_product[product],
                      "parent_return_fraction_per_second_upper": source_upper})
    isolated_two_species = _inventory_model(
        envelope, parent=parent, product=product, route_id=route_id,
        forward_value=forward_value,
        channel_uppers={term["channel_id"]: term["upper_rate_s^-1"] for term in terms},
    )
    horizon = float(snapshot["target_duration_seconds"])
    try:
        integrated_upper = math.fsum(term["parent_return_fraction_per_second_upper"] for term in terms) * horizon
    except OverflowError as error:
        raise ValueError("Integrated reverse source upper bound overflows") from error
    if not math.isfinite(integrated_upper):
        raise ValueError("Integrated reverse source upper bound overflows")
    error_bound = min(1.0, integrated_upper)
    reasons = []
    if not upper_supported:
        reasons.append("reverse_upper_bound_not_supported")
    if not all_alternatives_supported:
        reasons.append("rate_basis_alternative_upper_bound_not_supported")
    if not inventory_supported:
        reasons.append("product_inventory_bound_not_supported")
    bound_supported = not reasons
    negligible = error_bound <= normalized_policy["max_abs_retention_error"]
    if bound_supported and not negligible:
        reasons.append("reverse_return_bound_exceeds_observable_tolerance")
    allowed = bound_supported and negligible
    parent_without_return = math.exp(-forward_value * horizon) if isolated_two_species else None
    result = {
        "schema_version": 1, "kind": "finite_horizon_reverse_return_assessment",
        "policy": normalized_policy, "condition_contract_id": condition_id,
        "condition_snapshot": snapshot, "route_id": route_id,
        "status": ("reverse_return_negligible_within_declared_model" if allowed else
                   "reverse_return_requires_resolution" if bound_supported else "unresolved"),
        "bound_scope": "declared_kinetic_model_only",
        "observable": "parent_retention", "observable_horizon_seconds": horizon,
        "bound_supported": bound_supported, "bounded_approximation_allowed": allowed,
        "maximum_abs_parent_retention_error": error_bound if bound_supported else None,
        "diagnostic_maximum_abs_parent_retention_error": error_bound,
        "integrated_reverse_source_upper_bound": integrated_upper,
        "parent_retention_without_return_under_declared_model": parent_without_return,
        "parent_retention_interval_under_declared_model": (
            [parent_without_return, min(1.0, parent_without_return + error_bound)]
            if bound_supported and parent_without_return is not None else None
        ),
        "terms": terms, "blocking_reasons": reasons,
        "pressure_applicability_verified": False, "rate_uncertainty_resolved": False,
        "conformer_completeness_verified": False, "network_completeness_verified": False,
        "conclusive_observable_allowed": False, "lifetime_claim_allowed": False,
        "limitations": [
            "The omission bound is relative to the declared fixed kinetic model, not physical truth.",
            "Source hashes prove artifact integrity, not accuracy or authenticity of the rate model.",
            "No pressure applicability, rate uncertainty, conformer completeness, or mechanism completeness is established.",
            "The bound applies only to parent retention at the exact supplied finite horizon, not lifetime or other observables.",
        ],
        "evidence_inputs": inputs,
    }
    result["scientific_hash"] = _digest(result)
    return result


def validate_reverse_return_assessment(assessment: Mapping[str, Any]) -> dict:
    """Reopen source hashes and exactly recompute a persisted assessment.

    A recomputed checksum alone cannot bless edited conclusions: all derived
    fields must agree with a fresh evaluation under the versioned policy.
    """
    record = _copy(_mapping(assessment, "assessment"))
    expected_hash = record.pop("scientific_hash", None)
    if expected_hash != _digest(record):
        raise ValueError("Reverse return assessment scientific hash changed")
    inputs = _mapping(record.get("evidence_inputs"), "assessment.evidence_inputs")
    if set(inputs) != {"condition_contract", "forward_rate", "reverse_envelope", "inventory_bound"}:
        raise ValueError("Reverse return assessment evidence inputs are incomplete")
    replay = assess_reverse_return_relevance(**inputs, policy=record.get("policy"))
    if {**record, "scientific_hash": expected_hash} != replay:
        raise ValueError("Reverse return assessment differs from scientific recomputation")
    return replay


__all__ = [
    "REVERSE_RELEVANCE_POLICY_ID", "MODEL_INVENTORY_ALGORITHM",
    "reverse_relevance_condition_snapshot", "assess_reverse_return_relevance",
    "validate_reverse_return_assessment",
]
