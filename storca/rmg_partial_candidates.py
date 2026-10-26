"""Recover structural discoveries from immutable earlier RMG attempts.

The final attempt owns all propagation and convergence evidence. Earlier
attempts may supply only direct, atom-resolved hypotheses for path verification.
"""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
import math
from pathlib import Path
import re

from .rmg_evidence import normalize_rmg_label


def _endpoint(item: dict) -> tuple[dict, Counter]:
    from rdkit import Chem

    label, smiles = item.get("label"), item.get("smiles")
    charge, multiplicity = item.get("charge"), item.get("multiplicity")
    if not isinstance(label, str) or not label or not isinstance(smiles, str) or not smiles:
        raise ValueError("endpoint_identity_missing")
    if type(charge) is not int or type(multiplicity) is not int or multiplicity < 1:
        raise ValueError("endpoint_electronic_state_missing")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or Chem.GetFormalCharge(molecule) != charge:
        raise ValueError("endpoint_smiles_or_charge_invalid")
    inventory = Counter((atom.GetAtomicNum(), atom.GetIsotope())
                        for atom in Chem.AddHs(molecule).GetAtoms())
    if any(number <= 0 for number, _ in inventory):
        raise ValueError("endpoint_contains_unspecified_atom")
    clean = {"label": label, "smiles": smiles, "charge": charge, "multiplicity": multiplicity}
    if item.get("stoichiometric_copy") is not None:
        clean["stoichiometric_copy"] = item["stoichiometric_copy"]
    return clean, inventory


def _sanitize(route: dict, reference: dict, index: int, conditions: dict) -> dict:
    if (route.get("atom_balance") or {}).get("status") != "passed":
        raise ValueError("earlier_route_atom_balance_not_established")
    if route.get("initiation_status") != "directly_present":
        raise ValueError("earlier_route_not_directly_present")
    if not isinstance(route.get("route_id"), str) or not route["route_id"]:
        raise ValueError("earlier_route_id_missing")
    if not route.get("reaction_equation"):
        raise ValueError("earlier_route_equation_missing")
    endpoints, inventories = {}, {}
    for side in ("reactants", "products"):
        items = (route.get("resolved_endpoints") or {}).get(side)
        if not isinstance(items, list) or not items:
            raise ValueError("earlier_route_endpoints_missing")
        prepared = [_endpoint(item) for item in items]
        endpoints[side] = [item for item, _ in prepared]
        inventories[side] = sum((inventory for _, inventory in prepared), Counter())
    if inventories["reactants"] != inventories["products"]:
        raise ValueError("earlier_route_atom_or_isotope_imbalance")
    charges = {side: sum(item["charge"] for item in items) for side, items in endpoints.items()}
    if charges["reactants"] != charges["products"]:
        raise ValueError("earlier_route_charge_imbalance")
    target = normalize_rmg_label(str(conditions.get("target_label") or ""))
    composition = conditions.get("composition") or {}
    initial = {normalize_rmg_label(str(label)) for label, fraction in composition.items()
               if math.isfinite(float(fraction)) and float(fraction) > 0.0}
    if not target or target not in initial:
        raise ValueError("declared_initial_target_missing")
    labels = {side: Counter(normalize_rmg_label(item["label"]) for item in items)
              for side, items in endpoints.items()}
    equation_sides = re.split(r"<=>|=>", route["reaction_equation"], maxsplit=1)
    if len(equation_sides) != 2:
        raise ValueError("earlier_route_equation_invalid")
    for side, text in zip(("reactants", "products"), equation_sides):
        equation_counts = Counter()
        for token in text.split("+"):
            match = re.fullmatch(r"\s*(?:(\d+)\s+)?(.+?)\s*", token)
            if match is None:
                raise ValueError("earlier_route_equation_invalid")
            equation_counts[normalize_rmg_label(match.group(2))] += int(match.group(1) or 1)
        if equation_counts != labels[side]:
            raise ValueError("earlier_route_equation_endpoint_mismatch")
    if any(label not in initial for label in labels["reactants"]):
        raise ValueError("earlier_route_requires_undeclared_reactant")
    target_count = labels["reactants"][target]
    if target_count <= labels["products"][target]:
        raise ValueError("earlier_route_not_net_target_loss")
    if float(route.get("target_stoichiometry") or 0) != target_count:
        raise ValueError("earlier_route_target_stoichiometry_mismatch")
    if (route.get("electronic_state_validation") or {}).get("status") == "inconsistent":
        raise ValueError("earlier_route_electronic_surface_inconsistent")
    # Deliberate structural whitelist: no historical rate, flux, concentration,
    # ORCA acceptance, barrier, or inferred kinetic relevance crosses this boundary.
    clean = {key: deepcopy(route[key]) for key in (
        "route_id", "reaction_equation", "stable_route_identity", "reaction_family",
        "direction", "charge", "multiplicity", "surface_multiplicities",
        "electronic_state_validation",
    ) if key in route}
    clean.update({
        "source_attempt": deepcopy(reference), "source_route_index": index,
        "discovery_only": True, "kinetic_relevance": "unknown",
        "initiation_status": "directly_present", "target_stoichiometry": target_count,
        "resolved_endpoints": endpoints,
        "atom_balance": {"status": "passed", "basis": "recomputed_element_and_isotope_inventory"},
        "co_reactants": [{"label": label, "normalized_label": label,
                          "stoichiometry": count, "origin": "initial_scenario_species"}
                         for label, count in sorted(labels["reactants"].items()) if label != target],
        "orca_verification": {"status": "candidate_unverified", "endpoint_structures_resolved": True,
                              "atom_mapping_resolved": False, "can_execute": True},
        "discovery_provenance": {"source": "rmg", "scope": "earlier_bounded_attempt",
                                 "source_attempt": deepcopy(reference), "source_route_index": index},
    })
    for side, prefix in (("reactants", "reactant"), ("products", "product")):
        clean[prefix + "_smiles"] = [item["smiles"] for item in endpoints[side]]
        clean[prefix + "_charges"] = [item["charge"] for item in endpoints[side]]
        clean[prefix + "_multiplicities"] = [item["multiplicity"] for item in endpoints[side]]
    return clean


def build_retained_partial_routes(evidence: dict, condition_contract: dict) -> dict:
    """Rebuild discoveries from source records; an invalid source rejects the view."""
    from .rmg_attempts import load_rmg_attempt

    controller = evidence.get("controller") or {}
    references = controller.get("attempt_evidence") or []
    result = {"routes": [], "blocking_reasons": [], "findings": []}
    if not references:
        return result
    final_number = controller.get("final_attempt_number")
    if type(final_number) is not int or final_number < 1:
        result["blocking_reasons"].append("final_rmg_attempt_number_missing")
        return result
    routes, seen_attempts = [], set()
    for reference in references:
        try:
            number = reference.get("attempt_number")
            if type(number) is not int or number < 1 or number in seen_attempts or number > final_number:
                raise ValueError("invalid_or_duplicate_rmg_attempt_number")
            seen_attempts.add(number)
            if number == final_number:
                continue
            source = load_rmg_attempt(reference, condition_contract=condition_contract)
        except Exception as error:
            result["blocking_reasons"].append("earlier_attempt_evidence_invalid")
            result["findings"].append({"status": "rejected_source", "source_attempt": reference,
                                       "reason": f"{type(error).__name__}: {error}"})
            continue
        if not source.get("candidate_routes"):
            continue
        artifacts = source.get("artifacts") or {}
        if ((source.get("mechanism_inspection") or {}).get("status") != "completed"
                or any(not artifacts.get(role) or not Path(artifacts[role]).is_file()
                       for role in ("input", "chemkin", "species_dictionary"))):
            result["findings"].append({"status": "ineligible_source", "source_attempt": reference,
                                       "reason": "earlier_route_source_artifacts_or_inspection_missing"})
            continue
        for index, route in enumerate(source.get("candidate_routes") or []):
            try:
                clean = _sanitize(route, reference, index, condition_contract)
            except Exception as error:
                result["findings"].append({"status": "ineligible_route", "source_attempt": reference,
                                           "source_route_index": index,
                                           "reason": f"{type(error).__name__}: {error}"})
            else:
                routes.append(clean)
    if not result["blocking_reasons"]:
        result["routes"] = routes
    result["blocking_reasons"] = sorted(set(result["blocking_reasons"]))
    return result


def validated_retained_partial_routes(evidence: dict, condition_contract: dict) -> dict:
    """Admit only the exact current source-bound view stored by the controller."""
    rebuilt = build_retained_partial_routes(evidence, condition_contract)
    stored = evidence.get("retained_partial_routes") or []
    if rebuilt["routes"] != stored:
        rebuilt["routes"] = []
        rebuilt["blocking_reasons"] = sorted(set([
            *rebuilt["blocking_reasons"], "retained_partial_route_view_mismatch",
        ]))
    return rebuilt


def _structural_channel_identity(route: dict) -> tuple:
    """Compare physical channels independent of labels and reversible direction."""
    from .reaction_candidates import canonical_smiles

    sides = []
    for side in ("reactants", "products"):
        items = (route.get("resolved_endpoints") or {}).get(side)
        if not isinstance(items, list) or not items:
            raise ValueError("final_route_resolved_endpoints_missing")
        clean = [_endpoint(item)[0] for item in items]
        sides.append(tuple(sorted((canonical_smiles(item["smiles"]), item["charge"], item["multiplicity"])
                                  for item in clean)))
    multiplicity = route.get("multiplicity")
    if type(multiplicity) is not int or multiplicity < 1:
        raise ValueError("route_surface_multiplicity_unresolved")
    surfaces = route.get("surface_multiplicities") or [multiplicity]
    if any(type(value) is not int or value < 1 for value in surfaces):
        raise ValueError("route_surface_multiplicities_invalid")
    return (tuple(sorted(sides)), multiplicity, tuple(sorted(set(surfaces))))


def retained_earlier_rmg_verification_queue(
    evidence: dict, condition_contract: dict, *, maximum_candidates: int,
) -> dict:
    """Expose earlier-only hypotheses to the existing bounded supplemental queue.

    The completed final model keeps ownership of all flux and propagation. No
    earlier numerical or verification result is copied. A catalog with missing
    endpoint/state identity cannot establish that a channel is earlier-only.
    """
    if type(maximum_candidates) is not int or maximum_candidates < 0:
        raise ValueError("Earlier-RMG verification ceiling must be a nonnegative integer")
    report = {"routes": [], "blocking_reasons": [], "findings": [],
              "eligible_route_count": 0, "queue_truncated": False,
              "unqueued_route_ids": [], "historical_kinetic_credit_allowed": False}
    if evidence.get("status") != "completed" or maximum_candidates == 0:
        return report
    retained = validated_retained_partial_routes(evidence, condition_contract)
    report["blocking_reasons"] = retained["blocking_reasons"]
    report["findings"] = list(retained["findings"])
    if retained["blocking_reasons"] or not retained["routes"]:
        return report
    if not isinstance(evidence.get("candidate_routes"), list) or not isinstance(evidence.get("network_routes"), list):
        report["blocking_reasons"].append("final_route_catalog_unavailable")
        return report
    final_identities = set()
    for route in [*evidence["candidate_routes"], *evidence["network_routes"]]:
        try:
            final_identities.add(_structural_channel_identity(route))
        except Exception as error:
            report["blocking_reasons"].append("final_route_identity_unresolved")
            report["findings"].append({"status": "unresolved_final_route_identity",
                                       "route_id": route.get("route_id") if isinstance(route, dict) else None,
                                       "reason": f"{type(error).__name__}: {error}"})
    if report["blocking_reasons"]:
        report["blocking_reasons"] = sorted(set(report["blocking_reasons"]))
        return report
    unique, seen = [], set()
    for index, route in enumerate(retained["routes"]):
        try:
            identity = _structural_channel_identity(route)
        except Exception as error:
            report["blocking_reasons"].append("earlier_route_surface_identity_unresolved")
            report["findings"].append({"status": "ineligible_earlier_route_identity",
                                       "route_id": route.get("route_id"), "reason": f"{type(error).__name__}: {error}"})
            continue
        if identity in final_identities or identity in seen:
            report["findings"].append({"status": "already_present_in_final_model" if identity in final_identities else "duplicate_earlier_discovery",
                                       "route_id": route["route_id"], "retained_partial_route_index": index})
            continue
        seen.add(identity)
        queued = deepcopy(route)
        queued.update({
            "retained_partial_route_index": index,
            "screening_status": "candidate_unverified",
            "verification_queue_class": "earlier_rmg_attempt_discovery_kinetic_relevance_unknown",
            "verification_condition_contract_id": condition_contract.get("contract_id"),
            "selection_reason": "Source-bound earlier RMG channel absent from the completed final model; kinetic relevance remains unknown.",
        })
        unique.append(queued)
    if report["blocking_reasons"]:
        report["blocking_reasons"] = sorted(set(report["blocking_reasons"]))
        return report
    # Stable identity ordering; no historical kinetics or a made-up cost ranks work.
    unique.sort(key=lambda route: (route["route_id"], route["source_attempt"]["attempt_number"],
                                   route["source_route_index"]))
    report["eligible_route_count"] = len(unique)
    report["routes"] = unique[:maximum_candidates]
    report["queue_truncated"] = len(unique) > maximum_candidates
    report["unqueued_route_ids"] = [route["route_id"] for route in unique[maximum_candidates:]]
    return report
