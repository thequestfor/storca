"""Condition-bound forward-only approximations, never reversible rate repairs.

The immutable source library remains inadmissible when its conformer basis is
unresolved. A *different* library may omit a reverse channel only relative to
an explicitly declared closed, conservative, first-order model. This module
does not establish high-pressure applicability, physical uncertainty, conformer
equilibration, or a lifetime. Its certificate is for parent retention only.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
import math
from pathlib import Path
import re

from .reverse_relevance import (
    MODEL_INVENTORY_ALGORITHM, assess_reverse_return_relevance,
    reverse_relevance_condition_snapshot, validate_reverse_return_assessment,
)


CERTIFICATE_KIND = "bounded_parent_retention_forward_only_approximation"


class BoundedHandoffUnresolved(ValueError):
    """A valid assessment can refuse admission without discarding its evidence."""

    def __init__(self, message, assessment):
        super().__init__(message)
        self.assessment = assessment


def _artifact(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _reopen(record):
    if not isinstance(record, dict) or _artifact(record["path"]) != record:
        raise ValueError("bounded_handoff_source_missing_or_mutated")
    return Path(record["path"])


def _save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Append-only scientific evidence; callers must select a new workspace.
    with path.open("x") as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return path


def _label(value):
    return re.sub(r"\(\d+\)$", "", str(value).strip())


def _positive(value, field, *, zero=False):
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise ValueError(f"{field} is not a finite rate/constant")
    value = float(value)
    if not math.isfinite(value) or value < 0 or (not zero and value == 0):
        raise ValueError(f"{field} is not a finite positive rate/constant")
    return value


def _entry(path):
    """Read a single declarative entry; never execute an RMG Python library."""
    tree = ast.parse(Path(path).read_text())
    calls = [node.value for node in tree.body if isinstance(node, ast.Expr)
             and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name)
             and node.value.func.id == "entry"]
    if len(calls) != 1 or calls[0].args:
        raise ValueError("bounded_handoff_requires_one_declarative_entry")
    fields = {item.arg: item.value for item in calls[0].keywords}
    if None in fields or len(fields) != len(calls[0].keywords):
        raise ValueError("bounded_handoff_entry_keywords_invalid")
    kinetics = fields.get("kinetics")
    if not isinstance(kinetics, ast.Call) or not isinstance(kinetics.func, ast.Name) \
            or kinetics.func.id != "Arrhenius":
        raise ValueError("bounded_handoff_requires_retained_arrhenius")
    # Reject executable expressions even though we do not execute the source.
    for node in ast.walk(kinetics):
        if isinstance(node, (ast.Call, ast.Name)) and node is not kinetics and node is not kinetics.func:
            raise ValueError("bounded_handoff_kinetics_expression_not_declarative")
    return fields


def _source_library(path):
    path = Path(path).resolve(strict=True)
    manifest = json.loads((path / "generated-kinetics.json").read_text())
    if manifest.get("source") != "arkane_tst" or manifest.get("verification_status") != "verified_arkane_tst" \
            or manifest.get("molecularity") != 1:
        raise ValueError("bounded_handoff_requires_verified_first_order_arkane_source")
    for name in ("reactions.py", "dictionary.txt"):
        if _artifact(path / name)["sha256"] != (manifest.get("checksums") or {}).get(name):
            raise ValueError("bounded_handoff_original_library_checksum_mismatch")
    sides = [manifest.get(side) for side in ("reactants", "products")]
    if any(not isinstance(side, list) or len(side) != 1 for side in sides):
        raise ValueError("bounded_handoff_requires_one_to_one_route")
    entry = _entry(path / "reactions.py")
    expected = f"{sides[0][0]} <=> {sides[1][0]}"
    if ast.literal_eval(entry.get("label")) != expected:
        raise ValueError("bounded_handoff_source_direction_mismatch")
    basis = manifest.get("conformer_rate_basis") or {}
    if basis.get("kind") != "reached_endpoint_conformer_rate_basis" \
            or basis.get("network_reinjection_allowed") is not False:
        raise ValueError("bounded_handoff_requires_explicit_unresolved_conformer_basis")
    # Product conformer uncertainty does not authorize ignoring reactant
    # population uncertainty. The current extension is product-side only.
    reactants = (basis.get("comparisons") or {}).get("reactants") or []
    if len(reactants) != 1 or reactants[0].get("status") != "reference_geometry_correspondence":
        raise ValueError("bounded_handoff_reactant_conformer_basis_unresolved")
    for records in (basis.get("comparisons") or {}).values():
        for comparison in records:
            if comparison.get("reference_xyz") and comparison.get("reference_xyz_sha256"):
                _reopen({"path": comparison["reference_xyz"], "sha256": comparison["reference_xyz_sha256"]})
    return manifest, entry


def _model_sources(evidence, source_file):
    artifacts = evidence.get("artifacts") or {}
    sources = [_artifact(source_file)]
    for field in ("chemkin", "species_dictionary"):
        if not artifacts.get(field):
            raise ValueError("bounded_handoff_native_model_artifacts_missing")
        sources.append(_artifact(artifacts[field]))
    for field in ("transport", "input"):
        if artifacts.get(field):
            sources.append(_artifact(artifacts[field]))
    native_sources = (evidence.get("mechanism_inspection") or {}).get("equilibrium_model_evidence")
    if native_sources is not None:
        if native_sources.get("status") != "completed" or native_sources.get("source_artifacts_unchanged") is not True:
            raise ValueError("bounded_handoff_native_model_source_integrity_unresolved")
        for record in native_sources.get("source_artifacts", {}).values():
            _reopen(record)
    return sources


def _inventory(evidence, condition_contract, sources):
    phase = " ".join(str(condition_contract.get("phase_model", "")).lower().replace("-", " ").split())
    supported_phases = {"gas", "gas phase", "homogeneous gas", "homogeneous gas phase",
                        "homogeneous gas phase pure target model", "homogeneous gas phase ideal gas",
                        "homogeneous gas phase surrogate", "homogeneous gas phase low pressure surrogate"}
    if phase not in supported_phases \
            or condition_contract.get("light_condition") != "dark" \
            or condition_contract.get("light_model") is not None \
            or any(condition_contract.get(key) not in (None, "") for key in ("solvent", "solvent_model")) \
            or condition_contract.get("closed") is False or condition_contract.get("open_system") is True \
            or condition_contract.get("observable") not in (None, "parent_retention"):
        raise ValueError("bounded_handoff_environment_outside_scope")
    target = _label(condition_contract.get("target_label"))
    composition = condition_contract.get("composition")
    if composition != {condition_contract["target_label"]: 1.0}:
        raise ValueError("bounded_handoff_requires_closed_pure_parent_initial_inventory")
    reactor = (evidence.get("independent_cantera_propagation") or {}).get("reactor_model") or {}
    if reactor.get("policy_id") != "isothermal-ideal-gas-constant-pressure-v1" \
            or reactor.get("energy_enabled") is not False or reactor.get("condition_validation") != "passed" \
            or any(reactor.get(key) != condition_contract[key] for key in ("temperature_K", "pressure_bar")):
        raise ValueError("bounded_handoff_closed_isothermal_reactor_evidence_unresolved")
    inspection = evidence.get("mechanism_inspection") or {}
    if inspection.get("status") != "completed":
        raise ValueError("bounded_handoff_native_inspection_incomplete")
    species = [_label(item["chemkin_identifier"]) for item in inspection.get("species") or []]
    if target not in species or len(set(species)) != len(species):
        raise ValueError("bounded_handoff_species_inventory_invalid")
    channels = []
    for index, reaction in enumerate(inspection.get("reactions") or []):
        left, right = reaction.get("reactant_labels") or [], reaction.get("product_labels") or []
        if len(left) != 1 or len(right) != 1 or reaction.get("molecularity") != 1 \
                or (reaction.get("reversible") is True and reaction.get("reverse_molecularity") != 1):
            raise ValueError("bounded_handoff_nonlinear_or_nonconservative_model")
        left, right = _label(left[0]), _label(right[0])
        if left == right or left not in species or right not in species:
            raise ValueError("bounded_handoff_channel_species_invalid")
        expected = {key: condition_contract[key] for key in ("temperature_K", "pressure_bar")}
        if reaction.get("evaluated_at") != expected:
            raise ValueError("bounded_handoff_native_rates_condition_mismatch")
        for direction in ("forward", "reverse") if reaction.get("reversible") is True else ("forward",):
            if reaction.get(f"evaluated_{direction}_rate_units", reaction.get("evaluated_rate_units")) != "s^-1":
                raise ValueError("bounded_handoff_native_rate_units_invalid")
            value = _positive(reaction.get(f"evaluated_{direction}_rate_coefficient_si"), "native rate", zero=True)
            a, b = (left, right) if direction == "forward" else (right, left)
            channels.append({"channel_id": f"native-{index}-{direction}", "molecularity": 1,
                             "reactants": {a: 1}, "products": {b: 1}, "rate_s^-1": value})
    if not channels:
        raise ValueError("bounded_handoff_empty_model")
    return {"kind": "closed_first_order_linear_model_inventory", "algorithm": MODEL_INVENTORY_ALGORITHM,
            "isothermal": True, "closed": True, "species_ids": species,
            "initial_fractions": {item: float(item == target) for item in species},
            "channels": channels, "source_artifacts": sources}


def _rate_sources(record):
    sources = list((record.get("source_artifacts") or {}).values())
    basis = record.get("endpoint_artifact_basis") or {}
    for endpoint in [*(basis.get("reactants") or []), *(basis.get("products") or []),
                     *([basis["transition_state"]] if basis.get("transition_state") else [])]:
        sources.extend((endpoint.get("artifacts") or {}).values())
    if not sources:
        raise ValueError("bounded_handoff_retained_rate_sources_missing")
    unique = {item["path"]: item for item in sources}
    for item in unique.values():
        _reopen(item)
    return list(unique.values())


def _assessment(source_library, model_file, rate_file, condition_contract, policy=None):
    manifest, _ = _source_library(source_library)
    model_doc = json.loads(Path(model_file).read_text())
    if model_doc.get("condition_contract") != condition_contract:
        raise ValueError("bounded_handoff_model_condition_contract_mismatch")
    evidence = model_doc["rmg_evidence"]
    sources = _model_sources(evidence, model_file)
    model = _inventory(evidence, condition_contract, sources)
    rates = json.loads(Path(rate_file).read_text())
    forward, reverse = (rates.get(key) or {} for key in ("conditional_forward_rate_evidence", "reverse_rate_evidence"))
    for rate in (forward, reverse):
        if rate.get("status") != "conditional_central_model_rate" or rate.get("route_id") != manifest["route_id"] \
                or rate.get("model_scope") != "high_pressure_limit":
            raise ValueError("bounded_handoff_retained_fixed_rate_unresolved")
        condition = rate.get("condition_rate") or {}
        if condition.get("units") != "s^-1" or any(condition.get(key) != condition_contract[key]
                                                  for key in ("temperature_K", "pressure_bar")):
            raise ValueError("bounded_handoff_retained_fixed_rate_condition_mismatch")
    if forward["condition_rate"] != manifest["condition_rate"]:
        raise ValueError("bounded_handoff_original_forward_rate_changed")
    kf = _positive(forward["condition_rate"]["value"], "forward rate")
    kr = _positive(reverse["condition_rate"]["value"], "reached reverse rate", zero=True)
    parent, product = _label(manifest["reactants"][0]), _label(manifest["products"][0])
    if parent != _label(condition_contract["target_label"]):
        raise ValueError("bounded_handoff_forward_is_not_parent_loss")
    matches = [item for item in (evidence["mechanism_inspection"].get("reactions") or [])
               if {_label(x) for x in item.get("reactant_labels", []) + item.get("product_labels", [])}
               == {parent, product}]
    if len(matches) != 1 or matches[0].get("reversible") is not True:
        raise ValueError("bounded_handoff_native_reversible_channel_missing_or_ambiguous")
    native = matches[0]
    Kc = _positive(native.get("evaluated_equilibrium_constant_Kc"), "native equilibrium constant")
    if native.get("equilibrium_constant_evaluated_at") != {"temperature_K": condition_contract["temperature_K"]}:
        raise ValueError("bounded_handoff_native_equilibrium_condition_missing_or_mismatched")
    # Kc is read from species Gibbs thermochemistry, NOT a fitted reverse/forward ratio.
    network_reverse = kf / Kc if _label(native["reactant_labels"][0]) == parent else kf * Kc
    _positive(network_reverse, "network-species reverse rate", zero=True)
    forward_matches = [item for item in model["channels"] if item["reactants"] == {parent: 1}
                       and item["products"] == {product: 1}]
    return_matches = [item for item in model["channels"] if item["reactants"] == {product: 1}
                      and item["products"] == {parent: 1}]
    if len(forward_matches) != 1 or len(return_matches) != 1:
        raise ValueError("bounded_handoff_model_channel_count_ambiguous")
    forward_matches[0].update(channel_id=manifest["route_id"], **{"rate_s^-1": kf})
    return_matches[0].update(channel_id=manifest["route_id"] + ":return", **{"rate_s^-1": network_reverse})
    assumptions = ["Closed isothermal conservative first-order gas model, pure parent initially.",
                   "Fixed forward Arkane fit and two explicitly named reverse rate models only.",
                   "No claim about unknown conformers, missing chemistry, rate error, or pressure applicability."]
    binding = {"condition_contract_id": condition_contract["contract_id"],
               "condition_snapshot": reverse_relevance_condition_snapshot(condition_contract),
               "source_artifacts": [*sources, _artifact(rate_file), _artifact(Path(source_library) / "generated-kinetics.json")],
               "model_assumptions": assumptions}
    alternatives = [{"basis": name, "upper_rate_s^-1": value, "upper_bound_supported": True,
                     "upper_bound_method": "declared_fixed_model_envelope",
                     "bound_justification": justification, "source_artifacts": refs}
                    for name, value, justification, refs in (
                        ("reached_conformer", kr, "Exact evaluated retained reverse fit, bounded only as a fixed model.",
                         [*_rate_sources(reverse), _artifact(rate_file)]),
                        ("network_species", network_reverse, "New fixed forward fit divided by native species Kc at exact T.", sources))]
    inputs = {
        "condition_contract": condition_contract,
        "forward_rate": {**binding, "kind": "verified_conditional_forward_rate", "route_id": manifest["route_id"],
                         "parent_id": parent, "product_id": product, "molecularity": 1,
                         "stoichiometry": {"parent_consumed": 1, "product_formed": 1}, "value_s^-1": kf,
                         "verification_status": "verified_arkane_tst",
                         "source_artifacts": [*_rate_sources(forward), *binding["source_artifacts"]]},
        "reverse_envelope": {**binding, "kind": "first_order_reverse_rate_upper_envelope",
                             "upper_bound_supported": True, "bound_scope": "declared_kinetic_model_only",
                             "scope": {"closed": True, "isothermal": True, "observable": "parent_retention",
                                       "parent_id": parent, "forward_route_id": manifest["route_id"],
                                       "first_order_single_parent_product": True},
                             "channels": [{"channel_id": return_matches[0]["channel_id"], "product_id": product,
                                           "parent_stoichiometry": 1, "upper_rate_s^-1": max(kr, network_reverse),
                                           "rate_basis_alternatives": alternatives,
                                           "coverage": {"reached_conformer": True, "network_species": True,
                                                        "other_return_channels": "excluded_by_model_inventory"}}],
                             "model_inventory": model},
        "inventory_bound": {**binding, "kind": "finite_horizon_product_inventory_bound", "supported": True,
                            "basis": "initial_parent_inventory", "bound_method": "closed_one_to_one_conserved_inventory",
                            "products": [{"product_id": product, "maximum_fraction": 1.0}]},
    }
    return assess_reverse_return_relevance(**inputs, policy=policy), manifest


def write_forward_only_approximation_library(source_library, model_file, rate_file,
                                             condition_contract, workdir, *, policy=None):
    """Write a NEW irreversible library when a recomputable model bound permits it.

    model_file contains {condition_contract,rmg_evidence}, including a native
    Kc inspection. rate_file is saved read-only Arkane extraction evidence.
    Sources are not changed. Unsupported scopes raise ValueError, not a bypass.
    """
    assessment, original = _assessment(source_library, model_file, rate_file, condition_contract, policy)
    if not assessment["bounded_approximation_allowed"]:
        raise BoundedHandoffUnresolved("reverse_return_requires_resolution: "
                                      + ",".join(assessment["blocking_reasons"]), assessment)
    source_library = Path(source_library).resolve(strict=True)
    root = Path(workdir).resolve()
    if root == source_library or root.is_relative_to(source_library):
        raise ValueError("bounded_handoff_workspace_must_not_modify_source")
    library = root / ("storca-forward-only-" + assessment["scientific_hash"][:20])
    library.mkdir(parents=True, exist_ok=False)
    _, entry = _source_library(source_library)
    kinetics = ast.unparse(entry["kinetics"])
    equation = f"{original['reactants'][0]} => {original['products'][0]}"
    reaction_text = ("name = " + repr(library.name) + "\nshortDesc = 'Condition-bound model approximation; not a lifetime'\n"
                     "entry(index=1, label=" + repr(equation) + ", reversible=False, kinetics=" + kinetics + ",\n"
                     " shortDesc='Bounded omission of model reverse return',\n"
                     " longDesc='Conditional fixed model only; see generated-kinetics.json.')\n")
    with (library / "reactions.py").open("x") as handle:
        handle.write(reaction_text)
    with (library / "dictionary.txt").open("xb") as handle:
        handle.write((source_library / "dictionary.txt").read_bytes())
    certificate = {"schema_version": 1, "kind": CERTIFICATE_KIND, "assessment": assessment,
                   "maximum_abs_parent_retention_error": assessment["maximum_abs_parent_retention_error"],
                   "condition_contract": copy.deepcopy(condition_contract),
                   "original_rate_manifest": _artifact(source_library / "generated-kinetics.json"),
                   "source_model": _artifact(model_file), "retained_rate_evidence": _artifact(rate_file),
                   "conclusive_observable_allowed": False, "lifetime_claim_allowed": False}
    rates = json.loads(Path(rate_file).read_text())
    manifest = {**original, "source": "arkane_tst_forward_only_approximation", "purpose": "kinetic_repair",
                "library": str(library), "library_name": library.name, "source_library": library.name,
                "reaction_equation": equation, "kinetic_reversibility": "irreversible_model_approximation",
                "library_reaction_equation": equation,
                "claim_scope": "conditional_fixed_kinetic_model_no_physical_lifetime",
                "forward_only_approximation_evidence": certificate,
                "reverse_relevance_assessment": assessment, **rates,
                "checksums": {name: _artifact(library / name)["sha256"] for name in ("reactions.py", "dictionary.txt")}}
    if "rmg_library_validation" in manifest:
        manifest["original_rate_library_validation"] = manifest.pop("rmg_library_validation")
    _save(library / "generated-kinetics.json", manifest)
    validate_forward_only_approximation_evidence(manifest, condition_contract)
    return {**manifest, "manifest": str(library / "generated-kinetics.json")}


def validate_forward_only_approximation_evidence(evidence, condition_contract):
    """Reopen original evidence, recompute the bound and inspect derived direction."""
    certificate = evidence.get("forward_only_approximation_evidence") or {}
    if evidence.get("source") != "arkane_tst_forward_only_approximation" \
            or evidence.get("kinetic_reversibility") != "irreversible_model_approximation" \
            or evidence.get("claim_scope") != "conditional_fixed_kinetic_model_no_physical_lifetime" \
            or evidence.get("purpose") != "kinetic_repair":
        raise ValueError("bounded_handoff_semantic_claim_scope_invalid")
    if certificate.get("kind") != CERTIFICATE_KIND or condition_contract is None \
            or certificate.get("condition_contract") != condition_contract \
            or certificate.get("conclusive_observable_allowed") is not False \
            or certificate.get("lifetime_claim_allowed") is not False:
        raise ValueError("bounded_handoff_exact_condition_or_claim_scope_invalid")
    original_file = _reopen(certificate["original_rate_manifest"])
    model_file, rate_file = _reopen(certificate["source_model"]), _reopen(certificate["retained_rate_evidence"])
    assessment = validate_reverse_return_assessment(certificate["assessment"])
    replay, original = _assessment(original_file.parent, model_file, rate_file, condition_contract, assessment["policy"])
    if replay != assessment or not assessment["bounded_approximation_allowed"] \
            or certificate.get("maximum_abs_parent_retention_error") != assessment["maximum_abs_parent_retention_error"]:
        raise ValueError("bounded_handoff_assessment_not_supported_by_original_model")
    for field in ("condition_rate", "route_id", "reactants", "products", "molecularity",
                  "conformer_rate_basis", "verification_status", "replaces_reaction"):
        if evidence.get(field) != original.get(field):
            raise ValueError("bounded_handoff_original_rate_or_basis_changed: " + field)
    library = Path(evidence["library"]).resolve(strict=True)
    if library == original_file.parent:
        raise ValueError("bounded_handoff_original_library_cannot_be_relabelled")
    for name in ("reactions.py", "dictionary.txt"):
        if _artifact(library / name)["sha256"] != evidence["checksums"].get(name):
            raise ValueError("bounded_handoff_derived_library_checksum_mismatch")
    if _artifact(library / "dictionary.txt")["sha256"] != original["checksums"]["dictionary.txt"]:
        raise ValueError("bounded_handoff_species_dictionary_changed")
    derived, source = _entry(library / "reactions.py"), _entry(original_file.parent / "reactions.py")
    expected = f"{original['reactants'][0]} => {original['products'][0]}"
    if ast.literal_eval(derived.get("label")) != expected or ast.literal_eval(derived.get("reversible")) is not False \
            or ast.dump(derived["kinetics"]) != ast.dump(source["kinetics"]) \
            or evidence.get("reaction_equation") != expected \
            or evidence.get("library_reaction_equation") != expected:
        raise ValueError("bounded_handoff_requires_unchanged_forward_fit_and_irreversible_entry")
    return assessment


def validate_forward_only_repaired_model(repaired_evidence, replacement, condition_contract):
    """Do not carry a narrow certificate into a newly expanded/reversible model."""
    result = {"kind": "conditional_forward_only_repaired_model_validation", "status": "failed",
              "blocking_reasons": [], "conclusive_observable_allowed": False}
    try:
        assessment = validate_forward_only_approximation_evidence(replacement, condition_contract)
        sources = _model_sources(repaired_evidence, replacement["library"] + "/generated-kinetics.json")
        original = assessment["evidence_inputs"]["reverse_envelope"]["model_inventory"]
        forward = assessment["evidence_inputs"]["forward_rate"]
        parent, product = forward["parent_id"], forward["product_id"]
        inspection = repaired_evidence.get("mechanism_inspection") or {}
        expected_channel_count = len(original["channels"]) - len(
            assessment["evidence_inputs"]["reverse_envelope"]["channels"])
        actual_channel_count = sum(2 if item.get("reversible") is True else 1
                                   for item in inspection.get("reactions") or [])
        if len(inspection.get("species") or []) != len(original["species_ids"]) \
                or actual_channel_count != expected_channel_count:
            raise ValueError("bounded_handoff_repaired_inventory_changed_new_species_or_channels")
        inventory = _inventory(repaired_evidence, condition_contract, sources)
        from .flux_verification import verify_applied_replacement
        applied = verify_applied_replacement(repaired_evidence, replacement, condition_contract)
        if applied["status"] != "passed":
            raise ValueError("bounded_handoff_applied_replacement_unresolved: " + ",".join(applied["blocking_reasons"]))
        # A changed network needs a new certificate, even if still linear.
        # This prevents hidden new return channels or changed spectator rates.
        def signatures(model, omit_return=False):
            items = []
            for channel in model["channels"]:
                if omit_return and channel["reactants"] == {product: 1} and channel["products"] == {parent: 1}:
                    continue
                items.append((next(iter(channel["reactants"])), next(iter(channel["products"])), channel["rate_s^-1"]))
            return sorted(items)
        expected, actual = signatures(original, True), signatures(inventory)
        if set(inventory["species_ids"]) != set(original["species_ids"]) or len(expected) != len(actual) \
                or any(a[:2] != b[:2] or not math.isclose(a[2], b[2], rel_tol=1e-12, abs_tol=0)
                       for a, b in zip(expected, actual)):
            raise ValueError("bounded_handoff_repaired_inventory_changed_or_reverse_reintroduced")
        result.update(status="passed", derived_inventory=inventory, applied_replacement=applied,
                      maximum_abs_parent_retention_error=assessment["maximum_abs_parent_retention_error"],
                      rate_serialization_maximum_relative_difference=max(
                          abs(a[2] - b[2]) / max(abs(a[2]), abs(b[2]), 1e-300)
                          for a, b in zip(expected, actual)),
                      reverse_omission_bound_not_total_prediction_error=True,
                      source_evidence=sources)
    except (ValueError, KeyError, TypeError, OSError, SyntaxError) as error:
        result["blocking_reasons"].append(str(error))
    return result


EXPANDED_ASSESSMENT_KIND = "expanded_forward_only_fixed_model_assessment"
EXPANDED_ASSESSMENT_POLICY = "expanded-retained-first-order-forward-only-v2"


def assess_expanded_forward_only_model(model_file, retained_rate_file, replacement,
                                      condition_contract, rate_snapshot_file, policy=None):
    """Assess a NEW complete model without broadening the old library certificate.

The retained network already omits the selected direct product return. The
comparison inventory adds exactly that one virtual channel; every existing
entry, recycling/return and escape channel retains its operative fixed rate.
The resulting error bound is relative to two named reverse-rate models only.
No estimated branch acquires verification or uncertainty credit.
    """
    from .retained_linear_network import build_retained_linear_network, scientific_hash

    original = validate_forward_only_approximation_evidence(replacement, condition_contract)
    certificate = replacement["forward_only_approximation_evidence"]
    if _artifact(retained_rate_file) != certificate["retained_rate_evidence"]:
        raise ValueError("expanded_handoff_retained_rate_ancestry_mismatch")
    network = build_retained_linear_network(model_file, replacement, condition_contract, rate_snapshot_file)
    rates = json.loads(Path(retained_rate_file).read_text())
    forward, reverse = (rates.get(key) or {} for key in
                        ("conditional_forward_rate_evidence", "reverse_rate_evidence"))
    for rate in (forward, reverse):
        if rate.get("status") != "conditional_central_model_rate" \
                or rate.get("route_id") != replacement["route_id"] \
                or rate.get("model_scope") != "high_pressure_limit":
            raise ValueError("expanded_handoff_retained_fixed_rate_unresolved")
        condition = rate.get("condition_rate") or {}
        if condition.get("units") != "s^-1" or any(condition.get(key) != condition_contract[key]
                                                   for key in ("temperature_K", "pressure_bar")):
            raise ValueError("expanded_handoff_retained_fixed_rate_condition_mismatch")
    if forward["condition_rate"] != replacement["condition_rate"]:
        raise ValueError("expanded_handoff_original_forward_rate_changed")
    parent, product = network["parent_id"], network["product_id"]
    route_id = replacement["route_id"]
    selected = next(item for item in network["channels"]
                    if item["channel_id"] == network["replacement_forward_channel_id"])
    kf = _positive(selected["rate_s^-1"], "operative forward rate")
    kr = _positive(reverse["condition_rate"]["value"], "reached reverse rate", zero=True)
    native = network["applied_replacement"]["identity_resolution"]["candidate"]
    Kc = _positive(native.get("evaluated_equilibrium_constant_Kc"), "native equilibrium constant")
    if native.get("equilibrium_constant_evaluated_at") != {"temperature_K": condition_contract["temperature_K"]}:
        raise ValueError("expanded_handoff_native_equilibrium_condition_unresolved")
    network_reverse = _positive(kf / Kc, "native-species reverse rate", zero=True)
    virtual_id = route_id + ":virtual-direct-return"
    if virtual_id in {item["channel_id"] for item in network["channels"]}:
        raise ValueError("expanded_handoff_virtual_return_id_collision")
    # This projection changes only the selected channel's local identifier.
    # Rates, other channels, inert species and all initial fractions are kept.
    channels = [{key: copy.deepcopy(item[key]) for key in
                 ("channel_id", "molecularity", "reactants", "products", "rate_s^-1")}
                for item in network["channels"]]
    next(item for item in channels if item["channel_id"] == network["replacement_forward_channel_id"])["channel_id"] = route_id
    channels.append({"channel_id": virtual_id, "molecularity": 1,
                     "reactants": {product: 1}, "products": {parent: 1}, "rate_s^-1": network_reverse})
    sources = [*network["source_artifacts"], _artifact(retained_rate_file),
               certificate["original_rate_manifest"], certificate["source_model"]]
    unique = {item["path"]: item for item in sources}
    sources = sorted(unique.values(), key=lambda item: item["path"])
    model = {"kind": "closed_first_order_linear_model_inventory", "algorithm": MODEL_INVENTORY_ALGORITHM,
             "closed": True, "isothermal": True, "species_ids": network["species_ids"],
             "initial_fractions": network["initial_fractions"], "channels": channels,
             "source_artifacts": sources}
    assumptions = ["Complete retained closed, isothermal, conservative first-order fixed gas model; pure parent initially.",
                   "Only the selected direct product return is omitted; existing recycling channels remain explicit.",
                   "Two named fixed reverse-rate alternatives, not physical uncertainty bounds.",
                   "No pressure applicability, conformer equilibration, network completeness or lifetime claim."]
    binding = {"condition_contract_id": condition_contract["contract_id"],
               "condition_snapshot": reverse_relevance_condition_snapshot(condition_contract),
               "source_artifacts": sources, "model_assumptions": assumptions}
    alternatives = [
        {"basis": "reached_conformer", "upper_rate_s^-1": kr, "upper_bound_supported": True,
         "upper_bound_method": "declared_fixed_model_envelope",
         "bound_justification": "Evaluated retained reached-conformer reverse fit, fixed model only.",
         "source_artifacts": [*_rate_sources(reverse), _artifact(retained_rate_file)]},
        {"basis": "network_species", "upper_rate_s^-1": network_reverse, "upper_bound_supported": True,
         "upper_bound_method": "declared_fixed_model_envelope",
         "bound_justification": "Operative retained forward rate divided by direct native species Kc at exact T; named model only.",
         "source_artifacts": sources},
    ]
    inputs = {
        "condition_contract": condition_contract,
        "forward_rate": {**binding, "kind": "verified_conditional_forward_rate", "route_id": route_id,
                         "parent_id": parent, "product_id": product, "molecularity": 1,
                         "stoichiometry": {"parent_consumed": 1, "product_formed": 1}, "value_s^-1": kf,
                         "verification_status": "verified_arkane_tst",
                         "source_artifacts": [*_rate_sources(forward), *sources]},
        "reverse_envelope": {**binding, "kind": "first_order_reverse_rate_upper_envelope",
                             "upper_bound_supported": True, "bound_scope": "declared_kinetic_model_only",
                             "scope": {"closed": True, "isothermal": True, "observable": "parent_retention",
                                       "parent_id": parent, "forward_route_id": route_id,
                                       "first_order_single_parent_product": True},
                             "channels": [{"channel_id": virtual_id, "product_id": product,
                                           "parent_stoichiometry": 1, "upper_rate_s^-1": max(kr, network_reverse),
                                           "rate_basis_alternatives": alternatives,
                                           "coverage": {"reached_conformer": True, "network_species": True,
                                                        "other_return_channels": "excluded_by_model_inventory"}}],
                             "model_inventory": model},
        "inventory_bound": {**binding, "kind": "finite_horizon_product_inventory_bound", "supported": True,
                            "basis": "initial_parent_inventory", "bound_method": "closed_one_to_one_conserved_inventory",
                            "products": [{"product_id": product, "maximum_fraction": 1.0}]},
    }
    omission = assess_reverse_return_relevance(**inputs, policy=policy)
    result = {"schema_version": 2, "kind": EXPANDED_ASSESSMENT_KIND,
              "policy_id": EXPANDED_ASSESSMENT_POLICY,
              "status": "completed" if omission["bounded_approximation_allowed"] else "unresolved",
              "bounded_approximation_allowed": omission["bounded_approximation_allowed"],
              "blocking_reasons": omission["blocking_reasons"], "condition_contract": copy.deepcopy(condition_contract),
              "retained_network": network, "reverse_relevance_assessment": omission,
              "maximum_abs_parent_retention_error": omission["maximum_abs_parent_retention_error"],
              "omitted_channel_id": virtual_id,
              "original_certificate_scientific_hash": original["scientific_hash"],
              "source_artifacts": sources,
              "input_bindings": {"model_file": _artifact(model_file), "retained_rate_file": _artifact(retained_rate_file),
                                 "rate_snapshot_file": _artifact(rate_snapshot_file),
                                 "replacement": copy.deepcopy(replacement)},
              "reverse_omission_bound_not_total_prediction_error": True,
              "physical_rate_bounds_established": False, "uncertainty_calibrated": False,
              "rate_uncertainty_resolved": False, "pressure_applicability_verified": False,
              "conclusive_observable_allowed": False, "lifetime_claim_allowed": False,
              "terminal_robustness_eligible": False,
              "old_certificate_coverage_extended": False}
    result["scientific_hash"] = scientific_hash(result)
    return result


def validate_expanded_forward_only_model_assessment(assessment, condition_contract=None):
    """Recompute the expanded sidecar; never write it back into its ancestry."""
    if not isinstance(assessment, dict) or assessment.get("kind") != EXPANDED_ASSESSMENT_KIND \
            or assessment.get("policy_id") != EXPANDED_ASSESSMENT_POLICY:
        raise ValueError("expanded_handoff_assessment_policy_unsupported")
    contract = assessment.get("condition_contract") if condition_contract is None else condition_contract
    if contract != assessment.get("condition_contract"):
        raise ValueError("expanded_handoff_exact_condition_contract_mismatch")
    bindings = assessment["input_bindings"]
    replay = assess_expanded_forward_only_model(
        _reopen(bindings["model_file"]), _reopen(bindings["retained_rate_file"]),
        bindings["replacement"], contract, _reopen(bindings["rate_snapshot_file"]),
        policy=assessment["reverse_relevance_assessment"]["policy"])
    if replay != assessment:
        raise ValueError("expanded_handoff_assessment_not_supported_by_retained_sources")
    return replay
