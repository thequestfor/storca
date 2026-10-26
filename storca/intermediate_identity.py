"""Fail-closed molecular identity for validated path intermediates.

An optimized minimum is a stationary geometry, not yet a chemical species in
the native reaction network.  This module assigns a graph only when RDKit can
construct one molecule with the retained total charge and a spin-compatible
radical inventory.  The resulting graph is explicitly mapped back to the
parent route atom order so that the two elementary segment candidates can use
the ordinary route-verification machinery.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from .reaction_candidates import ReactionCandidate, canonical_smiles
from .route_geometry import resolve_route_atom_mapping, validate_route_endpoint_connectivity


INTERMEDIATE_IDENTITY_SCHEMA_VERSION = 1
INTERMEDIATE_IDENTITY_POLICY_ID = "rdkit-connectivity-first-charge-spin-map-v3"


def propose_unexpected_irc_candidate(route: Any, endpoint_check: dict, trace_evidence: dict) -> dict:
    """Retain one neutral closed-shell graph hypothesis, never a verified rate.

    This is discovery from a completed IRC, not a minimum-identity assertion.
    Geometry-derived bond orders and formal electronic states remain provisional.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from .decomposition_visuals import xyz_frames
    from .route_geometry import read_xyz
    import numpy as np

    base = {
        "schema_version": 1, "kind": "unexpected_irc_endpoint_discovery",
        "policy_id": "neutral_closed_shell_irc_graph_hypothesis-v1",
        "status": "unexpected_endpoint_unresolved", "native_candidates": [],
        "claims_permitted": [],
        "claims_prohibited": ["validated_product_identity", "verified_path", "rate", "lifetime"],
    }
    if not trace_evidence.get("valid") or endpoint_check.get("valid"):
        return {**base, "reason": "requires_completed_irc_and_failed_declared_endpoint_match"}
    if route.charge != 0 or route.multiplicity != 1:
        return {**base, "reason": "outside_neutral_singlet_graph_proposal_scope"}
    trajectory = Path(str(trace_evidence.get("trajectory") or ""))
    if not trajectory.is_file() or _sha256(trajectory) != trace_evidence.get("trajectory_sha256"):
        return {**base, "reason": "completed_trajectory_hash_missing_or_mismatched"}
    mapping = resolve_route_atom_mapping(route)
    if not mapping.get("valid"):
        return {**base, "reason": "parent_mapping_unresolved"}
    paths = [Path(str(endpoint_check.get(key) or ""))
             for key in ("first_endpoint_xyz", "last_endpoint_xyz")]
    if any(not path.is_file() for path in paths):
        return {**base, "reason": "irc_endpoint_geometry_missing"}
    try:
        frames = xyz_frames(trajectory)
        if len(frames) < 2 or any(
            labels != frames[0][0] or coordinates.shape != (len(labels), 3)
            or not np.isfinite(coordinates).all()
            for labels, coordinates, _ in frames
        ):
            return {**base, "reason": "completed_trajectory_geometry_invalid"}
        for path, frame in zip(paths, (frames[0], frames[-1])):
            labels, coordinates, _ = read_xyz(path)
            if (labels != frame[0] or coordinates.shape != frame[1].shape
                    or not np.allclose(coordinates, frame[1], atol=1e-6, rtol=0)):
                return {**base, "reason": "endpoint_does_not_belong_to_completed_trace"}
    except (ValueError, IndexError, OSError) as error:
        return {**base, "reason": "completed_trajectory_geometry_invalid",
                "diagnostic": f"{type(error).__name__}: {error}"}
    anchors = [validate_route_endpoint_connectivity(route, path, side="reactants", mapping_result=mapping)
               for path in paths]
    matching = [index for index, check in enumerate(anchors) if check.get("valid")]
    if len(matching) != 1:
        return {**base, "reason": "unique_declared_reactant_anchor_not_established"}
    unknown = paths[1 - matching[0]]
    try:
        molecule = Chem.MolFromXYZBlock(unknown.read_text())
        if molecule is None:
            raise ValueError("invalid endpoint XYZ")
        rdDetermineBonds.DetermineConnectivity(molecule, charge=0)
        rdDetermineBonds.DetermineBondOrders(molecule, charge=0, allowChargedFragments=True, embedChiral=False)
        Chem.SanitizeMol(molecule)
        fragments = Chem.GetMolFrags(molecule, asMols=True)
        if any(Chem.GetFormalCharge(fragment) != 0 or
               any(atom.GetNumRadicalElectrons() for atom in fragment.GetAtoms())
               or not any(atom.GetAtomicNum() > 1 for atom in fragment.GetAtoms())
               for fragment in fragments):
            return {**base, "reason": "fragment_charge_or_spin_assignment_requires_escalation"}
        mapped, inventory = _mapped_intermediate_smiles(molecule)
        if inventory != tuple(range(1, int(mapping["reactant_skeleton_atom_count"]) + 1)):
            raise ValueError("inferred endpoint changed parent heavy-atom inventory")
        products = tuple(mapped.split("."))
        reactants = _mapped_endpoint_smiles(tuple(route.reactant_smiles))
        if _side_identity(products) in (_side_identity(reactants), _side_identity(route.product_smiles)):
            return {**base, "reason": "no_distinct_product_hypothesis"}
        evidence = {
            "source": "unexpected_completed_irc_endpoint", "parent_route_id": route.route_id,
            "source_xyz_sha256": _sha256(unknown),
            "source_trajectory_sha256": _sha256(trajectory),
            "assertion": "geometry_derived_candidate_only",
        }
        parent_bonds = {tuple(sorted(pair)) for pair in anchors[matching[0]]["declared_bonds"]}
        observed_bonds = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
                          for bond in molecule.GetBonds()}
        observed_changes = {
            "index_convention": "zero_based_original_irc_atom_order",
            "broken_bonds": sorted(parent_bonds - observed_bonds),
            "formed_bonds": sorted(observed_bonds - parent_bonds),
        }
        candidate = ReactionCandidate.create(
            reactants, products, total_charge=0, multiplicity=1,
            discovery_assertions=(evidence,), screening_status="candidate_unverified",
            verification_status="unexpected_channel_unverified", kinetic_relevance="unknown",
            reactant_labels=tuple(route.reactant_labels), reactant_charges=tuple(route.reactant_charges),
            reactant_multiplicities=tuple(route.reactant_multiplicities),
            product_charges=(0,) * len(products), product_multiplicities=(1,) * len(products),
            environment=dict(route.environment or {}),
            warnings=("geometry_derived_bond_orders_unvalidated", "fragment_electronic_states_unvalidated",
                      "stationary_minima_not_yet_verified", "original_reaction_remains_unresolved",
                      "all_atom_trace_mapping_not_yet_admitted"),
            provenance={**evidence, "atom_mapping_representation": "smiles_atom_map_numbers",
                        "source_xyz": str(unknown), "source_trajectory": str(trajectory),
                        "observed_bond_changes": observed_changes,
                        "source_atom_order": [atom.GetSymbol() for atom in molecule.GetAtoms()],
                        "state_assignment_scope": "one_neutral_closed_shell_graph_hypothesis_not_exhaustive",
                        "parent_rate_reuse_allowed": False},
        )
        from .route_verify import route_spec_from_candidate_route
        from .route_geometry import align_mapping_to_retained_endpoints
        candidate_route = route_spec_from_candidate_route(
            route.parent_smiles, candidate.as_route_dict(parent_smiles=route.parent_smiles),
        )
        trace_mapping = align_mapping_to_retained_endpoints(candidate_route, paths[matching[0]], unknown)
        return {**base, "status": "candidate_proposed", "native_candidates": [candidate.as_dict()],
                "trace_atom_mapping": trace_mapping,
                "claims_permitted": ["unverified_alternative_reaction_candidate"],
                "reactant_anchor": anchors[matching[0]],
                "source_endpoint": str(unknown), "source_endpoint_sha256": _sha256(unknown)}
    except Exception as error:
        return {**base, "reason": "endpoint_graph_proposal_failed",
                "diagnostic": f"{type(error).__name__}: {error}"}


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _side_identity(smiles: Iterable[str]) -> Counter[str]:
    return Counter(canonical_smiles(item) for item in smiles)


def _allowed_multiplicities(unpaired_electrons: int) -> tuple[int, ...]:
    if unpaired_electrons < 0:
        return ()
    return tuple(range(unpaired_electrons + 1, 0, -2))


def _mapped_endpoint_smiles(
    smiles: tuple[str, ...], *, product_to_reactant: dict[int, int] | None = None,
) -> tuple[str, ...]:
    """Add common skeleton map numbers without changing endpoint chemistry."""
    from rdkit import Chem

    mapped: list[str] = []
    compact_index = 0
    for value in smiles:
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            raise ValueError(f"invalid endpoint SMILES: {value}")
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue
            common_index = (
                product_to_reactant.get(compact_index)
                if product_to_reactant is not None else compact_index
            )
            if common_index is None:
                raise ValueError("incomplete parent-route skeleton mapping")
            atom.SetAtomMapNum(int(common_index) + 1)
            compact_index += 1
        mapped.append(Chem.MolToSmiles(molecule, canonical=False, isomericSmiles=True))
    return tuple(mapped)


def _mapped_intermediate_smiles(molecule: Any) -> tuple[str, tuple[int, ...]]:
    """Return a mapped heavy-atom graph and its common atom-map inventory."""
    from rdkit import Chem

    common_index = 0
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            atom.SetAtomMapNum(0)
            continue
        atom.SetAtomMapNum(common_index + 1)
        common_index += 1
    heavy = Chem.RemoveHs(molecule, sanitize=True)
    mapped_smiles = Chem.MolToSmiles(heavy, canonical=False, isomericSmiles=True)
    reparsed = Chem.MolFromSmiles(mapped_smiles)
    if reparsed is None:
        raise ValueError("inferred mapped intermediate cannot be reparsed")
    inventory = tuple(sorted(
        atom.GetAtomMapNum() for atom in reparsed.GetAtoms() if atom.GetAtomicNum() != 1
    ))
    return mapped_smiles, inventory


def _failure(
    xyz_path: Path, status_detail: str, reason: str, *, diagnostics: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema_version": INTERMEDIATE_IDENTITY_SCHEMA_VERSION,
        "kind": "storca_intermediate_identity",
        "identity_policy_id": INTERMEDIATE_IDENTITY_POLICY_ID,
        "status": "intermediate_identity_unresolved",
        "status_detail": status_detail,
        "source_xyz": str(xyz_path),
        "source_xyz_sha256": _sha256(xyz_path),
        "reason": reason,
        "diagnostics": diagnostics or {},
        "native_candidates": [],
        "claims_permitted": ["validated_stationary_minimum"],
        "claims_prohibited": [
            "intermediate_species_identity", "elementary_mechanism", "rate", "lifetime",
        ],
    }


def resolve_validated_intermediate_identity(
    route: Any,
    xyz_path: Path,
    *,
    image_index: int | None = None,
    parent_evidence: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Resolve a validated minimum and create its two native segment candidates.

    Bond perception is deliberately conservative.  A disconnected geometry,
    an incompatible charge/spin assignment, a non-bijective parent mapping, or
    a non-bijective parent mapping remains unresolved. A graph equal to a
    declared endpoint is positively identified as an endpoint return, not
    misreported as an unresolved intermediate.
    Geometry-derived bond orders are evidence, not a substitute for subsequent
    TS/frequency/IRC verification of either segment.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    xyz_path = Path(xyz_path)
    if not xyz_path.is_file():
        return _failure(xyz_path, "source_geometry_missing", "validated intermediate XYZ is missing")
    mapping = resolve_route_atom_mapping(route)
    if not mapping.get("valid"):
        return _failure(
            xyz_path, "parent_mapping_unresolved",
            "the parent atom mapping is required to preserve intermediate identity",
            diagnostics={"parent_mapping": mapping},
        )
    try:
        molecule = Chem.MolFromXYZBlock(xyz_path.read_text())
        if molecule is None:
            raise ValueError("RDKit could not parse the XYZ geometry")
        rdDetermineBonds.DetermineConnectivity(
            molecule,
            charge=int(route.charge),
        )
    except Exception as error:
        return _failure(
            xyz_path, "connectivity_assignment_failed",
            f"geometry connectivity assignment failed: {type(error).__name__}: {error}",
        )

    fragments = Chem.GetMolFrags(molecule, asMols=False, sanitizeFrags=False)
    if len(fragments) != 1:
        return _failure(
            xyz_path, "fragmented_stationary_geometry",
            "a disconnected stationary geometry needs fragment-local charge and spin resolution",
            diagnostics={"fragment_count": len(fragments), "fragment_atom_indices": fragments},
        )

    try:
        rdDetermineBonds.DetermineBondOrders(
            molecule,
            charge=int(route.charge),
            allowChargedFragments=True,
            embedChiral=False,
        )
        Chem.SanitizeMol(molecule)
    except Exception as error:
        return _failure(
            xyz_path, "bond_assignment_failed",
            f"geometry-to-graph bond-order assignment failed: {type(error).__name__}: {error}",
        )

    formal_charge = int(Chem.GetFormalCharge(molecule))
    radical_electrons = int(sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()))
    allowed_multiplicities = _allowed_multiplicities(radical_electrons)
    declared_multiplicity = int(route.multiplicity)
    if formal_charge != int(route.charge):
        return _failure(
            xyz_path, "charge_mismatch",
            "the inferred graph does not retain the declared total charge",
            diagnostics={"declared_charge": route.charge, "inferred_charge": formal_charge},
        )
    if declared_multiplicity not in allowed_multiplicities:
        return _failure(
            xyz_path, "multiplicity_incompatible",
            "the inferred radical inventory is incompatible with the selected total multiplicity",
            diagnostics={
                "declared_multiplicity": declared_multiplicity,
                "radical_electrons": radical_electrons,
                "allowed_multiplicities": allowed_multiplicities,
            },
        )

    try:
        mapped_intermediate, map_inventory = _mapped_intermediate_smiles(molecule)
        skeleton_count = int(mapping["reactant_skeleton_atom_count"])
        if map_inventory != tuple(range(1, skeleton_count + 1)):
            raise ValueError("intermediate heavy-atom inventory differs from the parent route")
        parent_map = {
            int(reactant): int(product)
            for reactant, product in mapping["skeleton_atom_mapping"]
        }
        product_to_reactant = {product: reactant for reactant, product in parent_map.items()}
        if len(product_to_reactant) != skeleton_count:
            raise ValueError("parent route has a non-bijective skeleton mapping")
        mapped_reactants = _mapped_endpoint_smiles(tuple(route.reactant_smiles))
        mapped_products = _mapped_endpoint_smiles(
            tuple(route.product_smiles), product_to_reactant=product_to_reactant,
        )
    except Exception as error:
        return _failure(
            xyz_path, "mapped_identity_failed",
            f"atom-mapped identity construction failed: {type(error).__name__}: {error}",
            diagnostics={"parent_mapping": mapping},
        )

    intermediate_identity = _side_identity((mapped_intermediate,))
    endpoint_side = (
        "reactants" if intermediate_identity == _side_identity(mapped_reactants)
        else "products" if intermediate_identity == _side_identity(mapped_products)
        else None
    )
    if endpoint_side is not None:
        return {
            "schema_version": INTERMEDIATE_IDENTITY_SCHEMA_VERSION,
            "kind": "storca_intermediate_identity",
            "identity_policy_id": INTERMEDIATE_IDENTITY_POLICY_ID,
            "status": "endpoint_return",
            "status_detail": f"matches_declared_{endpoint_side}",
            "matched_endpoint_side": endpoint_side,
            "source_xyz": str(xyz_path),
            "source_xyz_sha256": _sha256(xyz_path),
            "reason": (
                "The ORCA-flagged candidate relaxed to a declared endpoint graph; "
                "it is not a distinct intermediate species."
            ),
            "diagnostics": {"mapped_intermediate_smiles": mapped_intermediate},
            "native_candidates": [],
            "claims_permitted": [
                "validated_stationary_minimum", "declared_endpoint_return",
            ],
            "claims_prohibited": [
                "intermediate_species_identity", "elementary_segment_candidates", "rate", "lifetime",
            ],
        }

    source_hash = _sha256(xyz_path)
    evidence_seed = {
        "parent_route_id": route.route_id,
        "parent_source": route.source,
        "source_image_index": image_index,
        "source_xyz_sha256": source_hash,
        "parent_mapping_source": mapping.get("source"),
        "bond_assignment_backend": "rdkit.DetermineBonds",
        "charge": formal_charge,
        "multiplicity": declared_multiplicity,
    }
    evidence_id = "intermediate-" + hashlib.sha256(
        json.dumps(evidence_seed, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()[:16]
    environment = dict(getattr(route, "environment", {}) or {})
    assertion = {
        "source": "validated_neb_intermediate",
        "evidence_id": evidence_id,
        "parent_route_id": route.route_id,
        "source_image_index": image_index,
        "assertion": "elementary_segment_candidate_only",
    }
    provenance = {
        **evidence_seed,
        "intermediate_evidence_id": evidence_id,
        "atom_mapping_representation": "smiles_atom_map_numbers",
        "parent_evidence": dict(parent_evidence or {}),
    }
    first = ReactionCandidate.create(
        mapped_reactants,
        (mapped_intermediate,),
        total_charge=formal_charge,
        multiplicity=declared_multiplicity,
        discovery_assertions=(assertion,),
        screening_status="candidate_unverified",
        verification_status="segment_path_unresolved",
        kinetic_relevance="unknown",
        reactant_labels=tuple(route.reactant_labels),
        product_labels=(f"{route.route_id}-intermediate-{image_index or 0}",),
        reactant_charges=tuple(route.reactant_charges),
        product_charges=(formal_charge,),
        reactant_multiplicities=tuple(route.reactant_multiplicities),
        product_multiplicities=(declared_multiplicity,),
        environment=environment,
        provenance={**provenance, "segment": "reactants_to_intermediate"},
    )
    second = ReactionCandidate.create(
        (mapped_intermediate,),
        mapped_products,
        total_charge=formal_charge,
        multiplicity=declared_multiplicity,
        discovery_assertions=(assertion,),
        screening_status="candidate_unverified",
        verification_status="segment_path_unresolved",
        kinetic_relevance="unknown",
        reactant_labels=(f"{route.route_id}-intermediate-{image_index or 0}",),
        product_labels=tuple(route.product_labels),
        reactant_charges=(formal_charge,),
        product_charges=tuple(route.product_charges),
        reactant_multiplicities=(declared_multiplicity,),
        product_multiplicities=tuple(route.product_multiplicities),
        environment=environment,
        provenance={**provenance, "segment": "intermediate_to_products"},
    )
    return {
        "schema_version": INTERMEDIATE_IDENTITY_SCHEMA_VERSION,
        "kind": "storca_intermediate_identity",
        "identity_policy_id": INTERMEDIATE_IDENTITY_POLICY_ID,
        "status": "resolved",
        "evidence_id": evidence_id,
        "source_xyz": str(xyz_path),
        "source_xyz_sha256": source_hash,
        "mapped_smiles": mapped_intermediate,
        "canonical_smiles": canonical_smiles(mapped_intermediate),
        "charge": formal_charge,
        "multiplicity": declared_multiplicity,
        "radical_electrons": radical_electrons,
        "allowed_multiplicities": list(allowed_multiplicities),
        "fragment_count": 1,
        "atom_map_numbers": list(map_inventory),
        "parent_mapping": mapping,
        "native_candidates": [first.as_dict(), second.as_dict()],
        "claims_permitted": ["intermediate_species_identity", "elementary_segment_candidates"],
        "claims_prohibited": ["verified_elementary_mechanism", "rate", "lifetime"],
        "interpretation": (
            "A charge-, spin-, and atom-map-consistent intermediate graph was resolved. "
            "Both elementary reactions remain unverified candidates until their own path and rate gates pass."
        ),
    }
