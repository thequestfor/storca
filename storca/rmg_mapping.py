"""Lazy RMG recipe evidence handoff to the existing all-H path mapper.

Original endpoints and their atom order are preserved. Exact endpoint graph
isomorphisms translate recipe IDs, never a graph-edit score. Multiple distinct
recipe maps and hydrogen assignments unsupported by the existing H layer fail
closed. A family trace is discovery evidence, not TS or rate evidence.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
import hashlib
import json
from pathlib import Path

from .route_verify import RouteSpec


POLICY_ID = "rmg-family-atom-trace-v1"


def _encoded(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _artifact(path: Path) -> dict:
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _check_artifacts(records: dict) -> None:
    if not isinstance(records, dict) or not records:
        raise ValueError("recipe_source_artifacts_missing")
    for record in records.values():
        if not isinstance(record, dict) or _artifact(Path(record["path"])) != record:
            raise ValueError("recipe_source_artifact_missing_or_mutated")


def _endpoint_contract(route: RouteSpec) -> dict:
    return {side: [{"label": label, "smiles": smiles, "charge": charge, "multiplicity": spin}
                   for label, smiles, charge, spin in zip(labels, smiles_values, charges, spins)]
            for side, labels, smiles_values, charges, spins in (
                ("reactants", route.reactant_labels, route.reactant_smiles,
                 route.reactant_charges, route.reactant_multiplicities),
                ("products", route.product_labels, route.product_smiles,
                 route.product_charges, route.product_multiplicities))}


def bind_rmg_mapping_source(route: dict, evidence: dict, condition_contract: dict | None) -> dict:
    """Attach only exact structural source inputs to an engine-selected route."""
    if route.get("native_candidate_id") or route.get("atom_mapping"):
        return route
    source, original = evidence, route
    if route.get("source_attempt"):
        from .rmg_attempts import load_rmg_attempt
        source = load_rmg_attempt(route["source_attempt"], condition_contract=condition_contract or {})
        index = route.get("source_route_index")
        if type(index) is not int or index < 0 or index >= len(source.get("candidate_routes") or []):
            raise ValueError("recipe_earlier_source_route_index_invalid")
        original = source["candidate_routes"][index]
        if original.get("resolved_endpoints") != route.get("resolved_endpoints"):
            # Earlier sanitized endpoints omit optional producer metadata.
            def structural(items):
                return [{key: item.get(key) for key in ("label", "smiles", "charge", "multiplicity")}
                        for item in items]
            if any(structural((original.get("resolved_endpoints") or {}).get(side) or []) !=
                   structural((route.get("resolved_endpoints") or {}).get(side) or [])
                   for side in ("reactants", "products")):
                raise ValueError("recipe_earlier_source_endpoint_mismatch")
    reaction_index = (original.get("retained_rate_evidence") or {}).get("reaction_index")
    artifacts = source.get("artifacts") or {}
    if (type(reaction_index) is not int or reaction_index < 0
            or not artifacts.get("chemkin") or not artifacts.get("species_dictionary")):
        return route  # Legacy source unavailable; no fabricated provenance.
    direction = "reverse" if original.get("direction") == "reverse_of_chemkin_direction" else "forward"
    return {**route, "rmg_mapping_source": {
        "reaction_index": reaction_index, "direction": direction,
        "reaction_family": original.get("reaction_family"),
        "source_artifacts": {"chemkin": _artifact(Path(artifacts["chemkin"])),
                             "dictionary": _artifact(Path(artifacts["species_dictionary"]))},
    }}


def _translated_trace(route: RouteSpec, trace: dict) -> tuple[dict[int, int], Counter]:
    from rdkit import Chem
    from .route_geometry import _build_endpoint_side

    ids_by_side = []
    endpoint_sides = []
    for side, smiles_values, charges, spins in (
            ("reactants", route.reactant_smiles, route.reactant_charges, route.reactant_multiplicities),
            ("products", route.product_smiles, route.product_charges, route.product_multiplicities)):
        originals = _build_endpoint_side(smiles_values)
        endpoint_sides.append(originals)
        remaining = list(trace[side])
        atom_ids = []
        for fragment, charge, spin in zip(originals.explicit_molecules, charges, spins):
            found = None
            for index, record in enumerate(remaining):
                if record["charge"] != charge or record["multiplicity"] != spin:
                    continue
                params = Chem.SmilesParserParams()
                params.removeHs = False
                mapped = Chem.MolFromSmiles(record["smiles"], params)
                if mapped is None:
                    raise ValueError("recipe_mapped_endpoint_invalid")
                mapped = Chem.AddHs(mapped)
                if fragment.GetNumAtoms() != mapped.GetNumAtoms():
                    continue
                match = mapped.GetSubstructMatch(fragment, useChirality=True)
                if not match:
                    continue
                # RDKit substructure matches are not sufficient by themselves
                # for isotope/radical/H-count identity. Check complete labeled
                # graphs in both directions before using this alignment.
                if not fragment.GetSubstructMatch(mapped, useChirality=True):
                    continue
                if any((fragment.GetAtomWithIdx(i).GetAtomicNum(), fragment.GetAtomWithIdx(i).GetIsotope(),
                        fragment.GetAtomWithIdx(i).GetFormalCharge(), fragment.GetAtomWithIdx(i).GetNumRadicalElectrons()) !=
                       (mapped.GetAtomWithIdx(j).GetAtomicNum(), mapped.GetAtomWithIdx(j).GetIsotope(),
                        mapped.GetAtomWithIdx(j).GetFormalCharge(), mapped.GetAtomWithIdx(j).GetNumRadicalElectrons())
                       for i, j in enumerate(match)):
                    continue
                found = (index, [mapped.GetAtomWithIdx(i).GetAtomMapNum() for i in match])
                break
            if found is None:
                raise ValueError("recipe_mapped_endpoint_identity_mismatch")
            remaining.pop(found[0])
            atom_ids.extend(found[1])
        if remaining or len(atom_ids) != len(originals.elements) or not all(atom_ids) \
                or len(set(atom_ids)) != len(atom_ids):
            raise ValueError("recipe_mapped_endpoint_nuclei_incomplete")
        ids_by_side.append(atom_ids)
    left, right = endpoint_sides
    if set(ids_by_side[0]) != set(ids_by_side[1]):
        raise ValueError("recipe_trace_atom_ids_not_bijective")
    product_by_id = {atom_id: index for index, atom_id in enumerate(ids_by_side[1])}
    all_mapping = {index: product_by_id[atom_id] for index, atom_id in enumerate(ids_by_side[0])}
    if any((left.elements[i], left.isotopes[i]) != (right.elements[j], right.isotopes[j])
           for i, j in all_mapping.items()):
        raise ValueError("recipe_trace_atom_or_isotope_mismatch")
    skeleton = {left.explicit_to_skeleton[i]: right.explicit_to_skeleton[j]
                for i, j in all_mapping.items() if left.elements[i] != "H"}
    transfers = Counter()
    for i, j in all_mapping.items():
        if left.elements[i] != "H":
            continue
        donor, acceptor = left.hydrogen_parent.get(i), right.hydrogen_parent.get(j)
        if (donor is None) != (acceptor is None) or (donor is not None and all_mapping[donor] != acceptor):
            transfers[(donor, acceptor, left.isotopes[i])] += 1
    return skeleton, transfers


def _resolve_report(route: RouteSpec, report: dict) -> tuple[tuple[int, int], ...]:
    from .route_geometry import (_build_endpoint_side, _extend_hydrogen_mapping,
                                 _skeleton_permutation_is_automorphism, _validate_skeleton_mapping)

    if (report.get("kind") != "rmg_family_atom_trace" or report.get("policy_id") != POLICY_ID
            or report.get("rate_claim_allowed") is not False or report.get("path_claim_allowed") is not False
            or report.get("status") != "completed"):
        raise ValueError("recipe_trace_not_completed_discovery_evidence")
    if route.source != "rmg" or not isinstance(report.get("reaction_family"), str) \
            or not report["reaction_family"]:
        raise ValueError("recipe_trace_source_or_family_missing")
    request = report.get("request") or {}
    if request.get("reaction_family") and request["reaction_family"] != report["reaction_family"]:
        raise ValueError("recipe_trace_family_contract_mismatch")
    if request.get("route_id") != route.route_id or request.get("endpoints") != _endpoint_contract(route):
        raise ValueError("recipe_trace_route_or_endpoint_contract_mismatch")
    if not all(request["endpoints"][side] for side in ("reactants", "products")):
        raise ValueError("recipe_endpoint_state_contract_incomplete")
    _check_artifacts(report.get("source_artifacts"))
    for role in ("family_sources", "implementation_sources"):
        sources = report.get(role)
        if not isinstance(sources, dict) or not sources:
            raise ValueError("recipe_family_or_implementation_sources_missing")
        for path, digest in sources.items():
            if _artifact(Path(path))["sha256"] != digest:
                raise ValueError("recipe_family_or_implementation_source_mutated")
    traces = report.get("traces")
    if not isinstance(traces, list) or not 1 <= len(traces) <= 64:
        raise ValueError("recipe_trace_count_invalid")
    left, right = _build_endpoint_side(route.reactant_smiles), _build_endpoint_side(route.product_smiles)

    def exact_symmetry(side, permutation, spins):
        if not _skeleton_permutation_is_automorphism(side, permutation):
            return False
        for i, j in permutation.items():
            first, second = side.skeleton_to_explicit[i], side.skeleton_to_explicit[j]
            if (spins[side.explicit_fragment[first]] != spins[side.explicit_fragment[second]]
                    or Counter(side.isotopes[h] for h in side.attached_hydrogens[first]) !=
                    Counter(side.isotopes[h] for h in side.attached_hydrogens[second])):
                return False
        return True

    def equivalent(previous, mapping):
        inverse = {j: i for i, j in mapping.items()}
        if exact_symmetry(left, {i: inverse[j] for i, j in previous.items()}, route.reactant_multiplicities):
            return True
        return exact_symmetry(right, {j: mapping[i] for i, j in previous.items()}, route.product_multiplicities)

    classes = []
    for trace in traces:
        mapping, transfers = _translated_trace(route, trace)
        reason = _validate_skeleton_mapping(left, right, mapping)
        if reason:
            raise ValueError(reason)
        extended = _extend_hydrogen_mapping(left, right, mapping)
        if not extended["valid"]:
            raise ValueError(extended["status"])
        inferred = Counter((item["reactant_donor_index"], item["product_acceptor_index"], item["isotope"])
                           for item in extended["hydrogen_transfers"])
        if inferred != transfers:
            raise ValueError("recipe_hydrogen_trace_not_supported_by_existing_hydrogen_layer")
        if not any(equivalent(previous, mapping) for previous in classes):
            classes.append(mapping)
    if len(classes) != 1:
        raise ValueError("recipe_multiple_inequivalent_atom_mappings")
    return tuple(sorted(classes[0].items()))


def validate_recipe_mapping(route: RouteSpec) -> None:
    """Reopen immutable evidence at geometry/reuse boundaries, not just dispatch."""
    evidence = route.atom_mapping_evidence
    record = evidence.get("artifact") or {}
    if _artifact(Path(record["path"])) != record:
        raise ValueError("recipe_mapping_evidence_missing_or_mutated")
    report = json.loads(Path(record["path"]).read_text())
    request_record = evidence.get("request_artifact") or {}
    if _artifact(Path(request_record["path"])) != request_record:
        raise ValueError("recipe_mapping_request_missing_or_mutated")
    payload = json.loads(Path(request_record["path"]).read_text())
    if (payload.get("request") != report.get("request")
            or payload.get("source_artifacts") != report.get("source_artifacts")
            or payload.get("operation") != "trace_family_recipe"):
        raise ValueError("recipe_mapping_request_result_mismatch")
    mapping = route.atom_mapping
    if (not isinstance(mapping, (tuple, list))
            or any(not isinstance(pair, (tuple, list)) or len(pair) != 2
                   or any(type(index) is not int for index in pair) for pair in mapping)):
        raise ValueError("recipe_mapping_indices_must_be_integer_pairs")
    # RouteSpec.as_dict() deliberately serializes tuples as JSON lists. Reopen
    # that representation without relaxing index, ordering or evidence checks.
    if _resolve_report(route, report) != tuple(tuple(pair) for pair in mapping):
        raise ValueError("recipe_mapping_does_not_match_evidence")


def apply_rmg_recipe_mapping(route: RouteSpec, source: dict, directory: Path, *,
                             rmg_env: str | None, timeout_seconds: float) -> tuple[RouteSpec, dict]:
    """Retain and consume one bounded native trace, reusing exact prior evidence."""
    from .rmg_bridge_client import run_rmg_bridge

    request = {"route_id": route.route_id, "endpoints": _endpoint_contract(route),
               **{key: source.get(key) for key in ("reaction_index", "direction", "reaction_family")}}
    for side, smiles_values in (("reactants", route.reactant_smiles), ("products", route.product_smiles)):
        if len(request["endpoints"][side]) != len(smiles_values):
            raise ValueError("recipe_endpoint_state_contract_incomplete")
    artifacts = source["source_artifacts"]
    _check_artifacts(artifacts)
    payload = {"operation": "trace_family_recipe", "request": request, "source_artifacts": artifacts,
               **{key: artifacts[key]["path"] for key in ("chemkin", "dictionary")}}
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(_encoded(payload)).hexdigest()[:24]
    request_path = directory / ("rmg-recipe-request-" + key + ".json")
    result_path = directory / ("rmg-recipe-result-" + key + ".json")
    encoded_request = _encoded(payload) + b"\n"
    try:
        with request_path.open("xb") as stream:
            stream.write(encoded_request)
    except FileExistsError:
        if request_path.read_bytes() != encoded_request:
            raise ValueError("recipe_retained_request_mutated")
    reused = result_path.is_file()
    if reused:
        report = json.loads(result_path.read_text())
    else:
        try:
            report = run_rmg_bridge(payload, rmg_env=rmg_env, timeout_seconds=timeout_seconds)
        except Exception as error:
            failure_path = directory / ("rmg-recipe-failure-" + key + ".json")
            if not failure_path.exists():
                with failure_path.open("xb") as stream:
                    stream.write(_encoded({"status": "failed", "stage": "native_family_recipe_trace",
                                           "request_artifact": _artifact(request_path),
                                           "reason": str(error), "exception_type": type(error).__name__,
                                           "stdout": getattr(error, "stdout", ""),
                                           "stderr": getattr(error, "stderr", ""),
                                           "path_claim_allowed": False, "rate_claim_allowed": False}) + b"\n")
            raise
        with result_path.open("xb") as stream:
            stream.write(_encoded(report) + b"\n")
    # Bind every reopened result to the exact dispatch, including reaction
    # index/direction; route endpoints alone cannot authorize source reuse.
    if report.get("request") != request or report.get("source_artifacts") != artifacts:
        raise ValueError("recipe_result_dispatch_source_mismatch")
    admission = {"policy_id": POLICY_ID, "artifact": _artifact(result_path),
                 "request_artifact": _artifact(request_path), "reused": reused,
                 "path_claim_allowed": False, "rate_claim_allowed": False}
    if report.get("status") == "unresolved":
        return route, {**admission, "status": "unresolved", "reason": report.get("reason")}
    mapping = _resolve_report(route, report)
    evidence = {"policy_id": POLICY_ID, "artifact": admission["artifact"],
                "request_artifact": admission["request_artifact"]}
    return replace(route, atom_mapping=mapping, atom_mapping_evidence=evidence), {
        **admission, "status": "resolved", "atom_mapping": [list(pair) for pair in mapping],
    }
