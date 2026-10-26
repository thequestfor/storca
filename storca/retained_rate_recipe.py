"""Read-only, hash-bound recipe facts; never execute retained Arkane Python.

An identified recipe is not a physically validated rate.  In particular, a
pressure grid is not an executed master equation, and optical, rotational and
reaction-path factors cannot be corrected independently from channel counting.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
import math
from pathlib import Path


POLICY_ID = "retained-arkane-recipe-audit-v1"
MAXIMUM_ARTIFACT_BYTES = 32 * 1024 * 1024
MAXIMUM_RECIPE_BYTES = 2 * 1024 * 1024
MAXIMUM_AST_NODES = 20000
MAXIMUM_AST_DEPTH = 96


class _Missing(ValueError):
    pass


class _Conflict(ValueError):
    pass


def _primitive(node, description):
    try:
        value = ast.literal_eval(node)
    except (ValueError, TypeError, SyntaxError) as error:
        raise _Conflict(f"nonliteral_recipe_metadata:{description}") from error
    return value


def _name(call):
    return call.func.id if isinstance(call.func, ast.Name) else None


def _keywords(call):
    result = {}
    for item in call.keywords:
        if item.arg is None or item.arg in result:
            raise _Conflict("duplicate_or_expanded_recipe_keywords")
        result[item.arg] = item.value
    return result


def _assignments(tree):
    result = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                raise _Conflict("ambiguous_recipe_assignment")
            name = node.targets[0].id
            if name in result:
                raise _Conflict(f"duplicate_recipe_declaration:{name}")
            result[name] = node.value
    return result


def _calls(tree, name):
    return [node.value for node in tree.body
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
            and _name(node.value) == name]


def _positive_integer(value, description):
    if type(value) is not int or value < 1:
        raise _Conflict(f"invalid_positive_integer:{description}")
    return value


def _positive_number(value, description):
    if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
        raise _Conflict(f"invalid_positive_number:{description}")
    return value


def _label(value, description):
    if not isinstance(value, str) or not value.strip():
        raise _Conflict(f"invalid_recipe_label:{description}")
    return value


def _unique_json(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise _Conflict(f"duplicate_provenance_key:{key}")
        result[key] = value
    return result


def _reopen(record, role, reopened):
    if not isinstance(record, dict) or not record.get("path") or not record.get("sha256"):
        raise _Missing(f"missing_source_binding:{role}")
    expected = record["sha256"]
    if not isinstance(expected, str) or len(expected) != 64 \
            or any(character not in "0123456789abcdef" for character in expected):
        raise _Conflict(f"invalid_source_hash:{role}")
    path = Path(record["path"]).resolve()
    try:
        if path.stat().st_size > MAXIMUM_ARTIFACT_BYTES:
            raise _Conflict(f"source_artifact_size_budget_exhausted:{role}")
        data = path.read_bytes()
    except OSError as error:
        raise _Missing(f"source_unavailable:{role}") from error
    actual = hashlib.sha256(data).hexdigest()
    if len(data) > MAXIMUM_ARTIFACT_BYTES:
        raise _Conflict(f"source_artifact_size_budget_exhausted:{role}")
    if actual != expected:
        raise _Conflict(f"source_hash_mismatch:{role}")
    reopened.append({"role": role, "path": str(path), "sha256": actual})
    return data, path


def _tree(data, role):
    if len(data) > MAXIMUM_RECIPE_BYTES:
        raise _Conflict(f"recipe_size_budget_exhausted:{role}")
    try:
        tree = ast.parse(data.decode("utf-8"))
    except (UnicodeError, SyntaxError, RecursionError, MemoryError) as error:
        raise _Conflict(f"malformed_retained_python:{role}") from error
    stack = [(tree, 0)]
    count = 0
    while stack:
        node, depth = stack.pop()
        count += 1
        if count > MAXIMUM_AST_NODES or depth > MAXIMUM_AST_DEPTH:
            raise _Conflict(f"recipe_ast_budget_exhausted:{role}")
        stack.extend((child, depth + 1) for child in ast.iter_child_nodes(node))
    allowed_calls = ({"species", "transitionState", "reaction", "kinetics", "network", "pressureDependence"}
                     if role == "input" else {"conformer", "kinetics"} if role == "output" else set())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            continue
        if isinstance(node, ast.Expr):
            if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                continue
            if isinstance(node.value, ast.Call) and _name(node.value) in allowed_calls:
                continue
            raise _Conflict(f"unknown_executable_recipe_statement:{role}")
        # Declarations hidden behind executable control flow are not facts
        # about the effective recipe. Imports and mutation statements are also
        # outside this bounded data grammar, even though nothing is executed.
        raise _Conflict(f"dynamic_recipe_control_flow:{role}")
    return tree


def _input_recipe(tree, route_id):
    assignments = _assignments(tree)
    facts = {}
    for name in ("modelChemistry", "frequencyScaleFactor", "useHinderedRotors",
                 "useAtomCorrections", "useBondCorrections"):
        if name not in assignments:
            raise _Missing(f"missing_recipe_declaration:{name}")
        facts[name] = _primitive(assignments[name], name)
    _label(facts["modelChemistry"], "modelChemistry")
    _positive_number(facts["frequencyScaleFactor"], "frequencyScaleFactor")
    for name in ("useHinderedRotors", "useAtomCorrections", "useBondCorrections"):
        if type(facts[name]) is not bool:
            raise _Conflict(f"invalid_boolean_recipe_declaration:{name}")
    definitions = {}
    for kind in ("species", "transitionState"):
        for call in _calls(tree, kind):
            if len(call.args) < 2:
                # Bath-gas species without statmech files are not endpoints.
                continue
            label = _label(_primitive(call.args[0], kind), kind)
            filename = _label(_primitive(call.args[1], label), label)
            if label in definitions:
                raise _Conflict(f"duplicate_stationary_label:{label}")
            definitions[label] = {"kind": kind, "statmech_file": filename}
    reactions = {}
    for call in _calls(tree, "reaction"):
        if len(call.args) != 1:
            raise _Conflict("ambiguous_reaction_label")
        label = _label(_primitive(call.args[0], "reaction"), "reaction")
        if label in reactions:
            raise _Conflict(f"duplicate_reaction_label:{label}")
        reactions[label] = _keywords(call)
    if route_id not in reactions:
        raise _Missing("selected_reaction_recipe_missing")
    declaration = reactions[route_id]
    reaction = {}
    for name in ("reactants", "products", "transitionState", "tunneling"):
        if name not in declaration:
            raise _Missing(f"selected_reaction_metadata_missing:{name}")
        reaction[name] = _primitive(declaration[name], name)
    for name in ("reactants", "products"):
        if not isinstance(reaction[name], (list, tuple)) or not reaction[name]:
            raise _Conflict(f"invalid_endpoint_labels:{name}")
        reaction[name] = [_label(label, name) for label in reaction[name]]
    _label(reaction["transitionState"], "transitionState")
    if not isinstance(reaction["tunneling"], str):
        raise _Conflict("invalid_tunneling_declaration")
    degeneracy = (None if "degeneracy" not in declaration else
                  _primitive(declaration["degeneracy"], "degeneracy"))
    if degeneracy is not None:
        _positive_number(degeneracy, "degeneracy")
    networks = _calls(tree, "network")
    pressure = _calls(tree, "pressureDependence")
    if bool(networks) != bool(pressure):
        raise _Conflict("network_and_pressure_job_metadata_inconsistent")
    # Do not infer a pressure job from a grid retained only in provenance.
    scope = "master_equation_recipe_declared" if pressure else "high_pressure_limit_recipe"
    return facts, definitions, reaction, degeneracy, scope, len(networks), len(pressure)


def _output_conformers(tree):
    result = {}
    for call in _calls(tree, "conformer"):
        kwargs = _keywords(call)
        if "label" not in kwargs:
            raise _Missing("output_conformer_label_missing")
        label = _label(_primitive(kwargs["label"], "output_label"), "output_label")
        if label in result:
            raise _Conflict(f"duplicate_output_conformer_label:{label}")
        values = {}
        for name in ("spin_multiplicity", "optical_isomers"):
            if name not in kwargs:
                raise _Missing(f"output_conformer_metadata_missing:{label}:{name}")
            values[name] = _positive_integer(_primitive(kwargs[name], name), name)
        modes = kwargs.get("modes")
        if not isinstance(modes, (ast.List, ast.Tuple)):
            raise _Missing(f"output_conformer_modes_missing:{label}")
        values["mode_types"] = []
        values["external_rotational_symmetry"] = None
        for mode in modes.elts:
            if not isinstance(mode, ast.Call) or _name(mode) is None:
                raise _Conflict(f"nonliteral_output_mode:{label}")
            mode_name = _name(mode)
            values["mode_types"].append(mode_name)
            if mode_name in {"NonlinearRotor", "LinearRotor"}:
                if values["external_rotational_symmetry"] is not None:
                    raise _Conflict(f"duplicate_external_rotor:{label}")
                symmetry = _keywords(mode).get("symmetry")
                if symmetry is None:
                    raise _Missing(f"rotational_symmetry_missing:{label}")
                values["external_rotational_symmetry"] = _positive_integer(
                    _primitive(symmetry, "rotational_symmetry"), "rotational_symmetry")
        result[label] = values
    return result


def audit_retained_arkane_recipe(replacement, condition_contract):
    """Return source-bound facts and explicit unresolved scientific questions.

    Missing legacy evidence is diagnostic, not an execution exception.  Supplied
    but contradictory or mutated evidence is refused as metadata conflict.
    No default path degeneracy is guessed from the presence of an Arkane file.
    """
    report = {
        "schema_version": 1, "kind": "retained_arkane_recipe_audit",
        "policy_id": POLICY_ID, "condition_contract": copy.deepcopy(condition_contract),
        "status": "evidence_missing", "blocking_reasons": [], "source_artifacts": [],
        "recipe_facts": None, "stationary_structures": [],
        "rate_validated": False, "lifetime_claim_allowed": False,
        "confidence_claim_allowed": False, "pressure_applicability_verified": False,
        "conformer_completeness_verified": False, "counting_basis_verified": False,
        "rate_correction_allowed": False, "network_reinjection_claim_allowed": False,
    }
    try:
        if not isinstance(replacement, dict) or not isinstance(condition_contract, dict):
            raise _Missing("replacement_or_condition_contract_missing")
        forward = replacement.get("conditional_forward_rate_evidence")
        if not isinstance(forward, dict):
            raise _Missing("conditional_forward_rate_evidence_missing")
        condition = forward.get("condition_rate")
        if not isinstance(condition, dict):
            raise _Missing("retained_condition_rate_missing")
        for name in ("temperature_K", "pressure_bar"):
            expected = _positive_number(condition_contract.get(name), name)
            actual = _positive_number(condition.get(name), name)
            if actual != expected:
                raise _Conflict(f"condition_contract_mismatch:{name}")
        sources = forward.get("source_artifacts")
        if not isinstance(sources, dict):
            raise _Missing("retained_recipe_sources_missing")
        data = {}
        paths = {}
        for role in ("input", "output", "provenance"):
            data[role], paths[role] = _reopen(sources.get(role), role, report["source_artifacts"])
        try:
            provenance = json.loads(data["provenance"].decode("utf-8"), object_pairs_hook=_unique_json)
        except (UnicodeError, json.JSONDecodeError) as error:
            raise _Conflict("malformed_recipe_provenance") from error
        if not isinstance(provenance, dict) or not isinstance(provenance.get("route"), dict):
            raise _Missing("retained_route_provenance_missing")
        route = provenance["route"]
        route_id = _label(forward.get("route_id"), "route_id")
        if route.get("label") != route_id:
            raise _Conflict("route_label_metadata_conflict")
        facts, definitions, reaction, degeneracy, scope, networks, pressure_jobs = _input_recipe(
            _tree(data["input"], "input"), route_id)
        for declared, key in (("modelChemistry", "model_chemistry"),
                              ("frequencyScaleFactor", "frequency_scale_factor")):
            if key not in route:
                raise _Missing(f"route_metadata_missing:{key}")
            if route[key] != facts[declared]:
                raise _Conflict(f"recipe_provenance_conflict:{key}")
        _positive_number(route["frequency_scale_factor"], "provenance_frequency_scale_factor")
        for name in ("temperature_K", "pressure_bar"):
            key = "condition_" + name
            if key not in route:
                raise _Missing(f"route_metadata_missing:{key}")
            if route[key] != condition_contract[name]:
                raise _Conflict(f"recipe_provenance_condition_conflict:{name}")
        for key, declared in (("reactant_labels", "reactants"), ("product_labels", "products")):
            if key not in route or list(route[key]) != reaction[declared]:
                raise _Conflict(f"recipe_provenance_endpoint_conflict:{key}")
        if route.get("tunneling") != reaction["tunneling"]:
            raise _Conflict("recipe_provenance_tunneling_conflict")
        if "pressure_dependent" in provenance and provenance["pressure_dependent"] is not bool(pressure_jobs):
            raise _Conflict("recipe_provenance_pressure_job_conflict")
        expected_scope = "high_pressure_limit" if not pressure_jobs else "master_equation_condition_rate"
        if forward.get("model_scope") not in (None, expected_scope):
            raise _Conflict("retained_rate_model_scope_conflict")
        conformers = _output_conformers(_tree(data["output"], "output"))
        basis = forward.get("endpoint_artifact_basis")
        if not isinstance(basis, dict):
            raise _Missing("endpoint_statmech_bindings_missing")
        endpoints = []
        for side in ("reactants", "products"):
            records = basis.get(side)
            if not isinstance(records, list) or not records:
                raise _Missing(f"endpoint_statmech_bindings_missing:{side}")
            if [record.get("label") for record in records if isinstance(record, dict)] != reaction[side]:
                raise _Conflict(f"endpoint_basis_label_conflict:{side}")
            endpoints.extend((side, record) for record in records)
        transition = basis.get("transition_state")
        if not isinstance(transition, dict):
            raise _Missing("transition_state_statmech_binding_missing")
        if transition.get("label") != reaction["transitionState"]:
            raise _Conflict("transition_state_basis_label_conflict")
        endpoints.append(("transition_state", transition))
        stationary = {}
        for side, endpoint in endpoints:
            label = _label(endpoint.get("label"), side)
            if label not in definitions or label not in conformers:
                raise _Missing(f"stationary_structure_recipe_missing:{label}")
            expected_kind = "transitionState" if side == "transition_state" else "species"
            if definitions[label]["kind"] != expected_kind:
                raise _Conflict(f"stationary_definition_kind_conflict:{label}")
            artifacts = endpoint.get("artifacts")
            if not isinstance(artifacts, dict) or "statmech_file" not in artifacts:
                raise _Missing(f"endpoint_statmech_source_missing:{label}")
            statmech = None
            statmech_path = None
            for name, record in artifacts.items():
                body, path = _reopen(record, f"{side}:{label}:{name}", report["source_artifacts"])
                if name == "statmech_file":
                    statmech, statmech_path = body, path
            expected_path = (paths["input"].parent / definitions[label]["statmech_file"]).resolve()
            if expected_path != statmech_path:
                raise _Conflict(f"endpoint_statmech_path_conflict:{label}")
            declarations = _assignments(_tree(statmech, f"statmech:{label}"))
            row = {"label": label, "sides": [side], **conformers[label]}
            for declared, key in (("opticalIsomers", "optical_isomers"),
                                  ("spin_multiplicity", "spin_multiplicity")):
                if declared not in declarations:
                    raise _Missing(f"stationary_declaration_missing:{label}:{declared}")
                value = _positive_integer(_primitive(declarations[declared], declared), declared)
                if value != row[key]:
                    raise _Conflict(f"stationary_output_declaration_conflict:{label}:{declared}")
            if _positive_integer(endpoint.get("multiplicity"), "endpoint_multiplicity") != row["spin_multiplicity"]:
                raise _Conflict(f"endpoint_multiplicity_metadata_conflict:{label}")
            external = declarations.get("externalSymmetry")
            row["external_symmetry_explicitly_declared"] = external is not None
            if external is not None and _positive_integer(_primitive(external, "externalSymmetry"),
                                                         "externalSymmetry") != row["external_rotational_symmetry"]:
                raise _Conflict(f"stationary_external_symmetry_conflict:{label}")
            row["optical_isomers_explicitly_declared"] = True
            row["independently_validated"] = False
            if label in stationary:
                previous = stationary[label]
                if {key: value for key, value in previous.items() if key != "sides"} != \
                        {key: value for key, value in row.items() if key != "sides"}:
                    raise _Conflict(f"reused_stationary_label_conflict:{label}")
                previous["sides"].append(side)
            else:
                stationary[label] = row
        for side, key in (("reactants", "reactant_multiplicities"), ("products", "product_multiplicities")):
            if key in route:
                values = route[key]
                if not isinstance(values, (list, tuple)) or len(values) != len(reaction[side]):
                    raise _Conflict(f"provenance_spin_inventory_conflict:{side}")
                for label, multiplicity in zip(reaction[side], values):
                    if _positive_integer(multiplicity, "provenance_multiplicity") != stationary[label]["spin_multiplicity"]:
                        raise _Conflict(f"provenance_spin_inventory_conflict:{label}")
        if "transition_state_multiplicity" in route:
            if _positive_integer(route["transition_state_multiplicity"], "provenance_ts_multiplicity") != \
                    stationary[reaction["transitionState"]]["spin_multiplicity"]:
                raise _Conflict("provenance_transition_state_spin_conflict")
        rrho = all(set(row["mode_types"]) in (
            {"IdealGasTranslation", "NonlinearRotor", "HarmonicOscillator"},
            {"IdealGasTranslation", "LinearRotor", "HarmonicOscillator"})
            and len(row["mode_types"]) == 3 for row in stationary.values())
        report.update({
            "status": "unresolved", "recipe_facts": {
                "electronic_structure_model": facts["modelChemistry"],
                "frequency_scale_factor": facts["frequencyScaleFactor"],
                "statistical_mechanics": ("retained_RRHO_mode_inventory" if rrho else
                                           "retained_mode_inventory_not_classified"),
                "hindered_rotors_used": facts["useHinderedRotors"],
                "atom_corrections_used": facts["useAtomCorrections"],
                "bond_corrections_used": facts["useBondCorrections"],
                "tunneling_model": reaction["tunneling"], "pressure_model": scope,
                "network_declaration_count": networks,
                "pressure_dependence_job_declaration_count": pressure_jobs,
                "pressure_grid_in_provenance_establishes_applicability": False,
                "bath_gas_applicable_to_recipe": bool(pressure_jobs),
                "reaction_path_degeneracy": {
                    "explicitly_declared": degeneracy is not None,
                    "declared_value": degeneracy, "effective_native_default": None,
                    "native_default_source_bound": False, "independently_validated": False,
                },
            },
            "stationary_structures": list(stationary.values()),
            "blocking_reasons": ["statistical_counting_basis_unresolved",
                                 "conformer_and_transition_state_channel_coverage_unresolved",
                                 "physical_rate_uncertainty_unresolved",
                                 "requested_pressure_applicability_unresolved"],
            "limitations": [
                "Retained declarations and native output facts are not independent scientific validation.",
                "Missing path degeneracy is unknown; no native default is guessed without version/source-bound evidence.",
                "Optical, external rotational symmetry and channel degeneracy need a joint no-double-counting audit.",
                "Fit-parameter errors do not establish physical uncertainty.",
                "A declared network/pressure job is not proof of execution or collisional applicability.",
                "No historical partition diagnostic, rate multiplier or experimental comparator is automatically imported.",
            ],
        })
    except _Missing as error:
        report["status"] = "evidence_missing"
        report["blocking_reasons"] = [str(error)]
    except (_Conflict, TypeError, KeyError, AttributeError) as error:
        report["status"] = "recipe_metadata_conflict"
        report["blocking_reasons"] = [str(error)]
    # Binding includes all supplied metadata as well as reopened file hashes.
    try:
        binding = json.dumps({"replacement": replacement, "condition_contract": condition_contract},
                             sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        report["input_binding_hash"] = hashlib.sha256(binding).hexdigest()
    except (ValueError, TypeError):
        # Invalid caller metadata cannot become a reproducible scientific
        # binding; diagnostic refusal must still not crash legacy reporting.
        report["status"] = "recipe_metadata_conflict"
        report["blocking_reasons"] = ["unhashable_or_nonfinite_recipe_metadata"]
        report["recipe_facts"] = None
        report["stationary_structures"] = []
        report["input_binding_hash"] = None
        report["condition_contract"] = None
    report["scientific_hash"] = hashlib.sha256(json.dumps(
        report, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    return report
