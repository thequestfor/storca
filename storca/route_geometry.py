"""Atom-consistent geometries for executable decomposition routes."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .route_verify import RouteSpec


@dataclass
class _EndpointSide:
    """One side of a reaction in stable, concatenated fragment atom order."""

    smiles: tuple[str, ...]
    skeleton_molecules: list[Any]
    explicit_molecules: list[Any]
    skeleton_to_explicit: dict[int, int]
    explicit_to_skeleton: dict[int, int]
    explicit_fragment: list[int]
    elements: list[str]
    isotopes: list[int]
    bonds: dict[tuple[int, int], float]
    attached_hydrogens: dict[int, list[int]]
    hydrogen_parent: dict[int, int | None]
    skeleton_count: int
    formal_charges: tuple[int, ...]
    inferred_multiplicities: tuple[int, ...]


def _bond_value(bond: Any) -> float:
    """Return a stable numerical bond order for graph-edit comparison."""
    return 1.5 if bond.GetIsAromatic() else float(bond.GetBondTypeAsDouble())


def _build_endpoint_side(smiles: tuple[str, ...]) -> _EndpointSide:
    """Construct explicit-H fragments without losing SMILES atom indices."""
    from rdkit import Chem

    if not smiles:
        raise ValueError("A reaction endpoint must contain at least one species")
    skeleton_molecules: list[Any] = []
    explicit_molecules: list[Any] = []
    skeleton_to_explicit: dict[int, int] = {}
    explicit_to_skeleton: dict[int, int] = {}
    explicit_fragment: list[int] = []
    elements: list[str] = []
    isotopes: list[int] = []
    bonds: dict[tuple[int, int], float] = {}
    attached_hydrogens: dict[int, list[int]] = {}
    hydrogen_parent: dict[int, int | None] = {}
    skeleton_offset = 0
    explicit_offset = 0
    charges: list[int] = []
    multiplicities: list[int] = []
    for fragment_index, value in enumerate(smiles):
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            raise ValueError(f"Invalid endpoint SMILES: {value}")
        explicit = Chem.AddHs(molecule)
        if explicit.GetNumAtoms() == 0:
            raise ValueError(f"Empty endpoint species: {value}")
        skeleton_molecules.append(molecule)
        explicit_molecules.append(explicit)
        # Skeleton indices are heavy-atom-only on both sides.  Explicit H in a
        # SMILES string (including a standalone [H], [H+], or [H-] fragment)
        # belongs to the same all-hydrogen layer as atoms introduced by AddHs.
        local_heavy_indices = [
            atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1
        ]
        for compact_index, local_index in enumerate(local_heavy_indices):
            skeleton_to_explicit[skeleton_offset + compact_index] = explicit_offset + local_index
            explicit_to_skeleton[explicit_offset + local_index] = skeleton_offset + compact_index
        for atom in explicit.GetAtoms():
            global_index = explicit_offset + atom.GetIdx()
            elements.append(atom.GetSymbol())
            isotopes.append(int(atom.GetIsotope()))
            explicit_fragment.append(fragment_index)
            if atom.GetAtomicNum() != 1:
                attached_hydrogens[global_index] = sorted(
                    explicit_offset + neighbour.GetIdx()
                    for neighbour in atom.GetNeighbors()
                    if neighbour.GetAtomicNum() == 1
                )
            else:
                heavy_neighbours = [
                    explicit_offset + neighbour.GetIdx()
                    for neighbour in atom.GetNeighbors()
                    if neighbour.GetAtomicNum() != 1
                ]
                hydrogen_parent[global_index] = heavy_neighbours[0] if len(heavy_neighbours) == 1 else None
        for bond in explicit.GetBonds():
            left = explicit_offset + bond.GetBeginAtomIdx()
            right = explicit_offset + bond.GetEndAtomIdx()
            bonds[tuple(sorted((left, right)))] = _bond_value(bond)
        charges.append(sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()))
        radicals = sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms())
        # Closed shells and monoradicals have one unambiguous minimal spin.
        # Two or more radical electrons can couple on multiple surfaces, so a
        # SMILES string alone is not an electronic-state declaration.
        multiplicities.append(radicals + 1 if radicals <= 1 else 0)
        skeleton_offset += len(local_heavy_indices)
        explicit_offset += explicit.GetNumAtoms()
    return _EndpointSide(
        smiles=smiles,
        skeleton_molecules=skeleton_molecules,
        explicit_molecules=explicit_molecules,
        skeleton_to_explicit=skeleton_to_explicit,
        explicit_to_skeleton=explicit_to_skeleton,
        explicit_fragment=explicit_fragment,
        elements=elements,
        isotopes=isotopes,
        bonds=bonds,
        attached_hydrogens=attached_hydrogens,
        hydrogen_parent=hydrogen_parent,
        skeleton_count=skeleton_offset,
        formal_charges=tuple(charges),
        inferred_multiplicities=tuple(multiplicities),
    )


def _skeleton_atoms(side: _EndpointSide) -> list[Any]:
    return [
        atom
        for molecule in side.skeleton_molecules
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() != 1
    ]


def _skeleton_bonds(side: _EndpointSide) -> dict[tuple[int, int], float]:
    bonds: dict[tuple[int, int], float] = {}
    offset = 0
    for molecule in side.skeleton_molecules:
        local_to_compact = {
            atom.GetIdx(): compact
            for compact, atom in enumerate(
                atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1
            )
        }
        for bond in molecule.GetBonds():
            if bond.GetBeginAtomIdx() not in local_to_compact or bond.GetEndAtomIdx() not in local_to_compact:
                continue
            pair = tuple(sorted((
                offset + local_to_compact[bond.GetBeginAtomIdx()],
                offset + local_to_compact[bond.GetEndAtomIdx()],
            )))
            bonds[pair] = _bond_value(bond)
        offset += len(local_to_compact)
    return bonds


def _neighbour_shells(side: _EndpointSide, index: int, depth: int = 2) -> tuple[Counter, ...]:
    atoms = _skeleton_atoms(side)
    bonds = _skeleton_bonds(side)
    adjacency: dict[int, set[int]] = {position: set() for position in range(len(atoms))}
    for left, right in bonds:
        adjacency[left].add(right)
        adjacency[right].add(left)
    shells: list[Counter] = []
    visited = {index}
    frontier = {index}
    for _ in range(depth):
        frontier = {neighbour for atom in frontier for neighbour in adjacency[atom]} - visited
        shells.append(Counter(atoms[item].GetSymbol() for item in frontier))
        visited.update(frontier)
    return tuple(shells)


def _counter_distance(left: Counter, right: Counter) -> int:
    return sum((left - right).values()) + sum((right - left).values())


def _atom_correspondence_cost(
    reactants: _EndpointSide,
    products: _EndpointSide,
    reactant_index: int,
    product_index: int,
) -> float:
    """Local graph-edit cost; finite only for element/isotope-compatible atoms."""
    left = _skeleton_atoms(reactants)[reactant_index]
    right = _skeleton_atoms(products)[product_index]
    if left.GetAtomicNum() != right.GetAtomicNum() or left.GetIsotope() != right.GetIsotope():
        return math.inf
    cost = 5.0 * abs(left.GetDegree() - right.GetDegree())
    cost += 6.0 * int(left.GetIsAromatic() != right.GetIsAromatic())
    cost += 4.0 * int(left.IsInRing() != right.IsInRing())
    cost += 2.0 * abs(left.GetFormalCharge() - right.GetFormalCharge())
    cost += 1.0 * abs(left.GetNumRadicalElectrons() - right.GetNumRadicalElectrons())
    cost += 0.5 * abs(left.GetTotalNumHs(includeNeighbors=True) - right.GetTotalNumHs(includeNeighbors=True))
    for radius, (first, second) in enumerate(
        zip(_neighbour_shells(reactants, reactant_index), _neighbour_shells(products, product_index)),
        start=1,
    ):
        cost += (4.0 / radius) * _counter_distance(first, second)
    return cost


def _validate_skeleton_mapping(
    reactants: _EndpointSide,
    products: _EndpointSide,
    mapping: dict[int, int],
) -> str | None:
    if set(mapping) != set(range(reactants.skeleton_count)):
        return "incomplete_atom_mapping"
    if len(set(mapping.values())) != products.skeleton_count:
        return "non_bijective_atom_mapping"
    if any(index < 0 or index >= products.skeleton_count for index in mapping.values()):
        return "atom_mapping_index_out_of_range"
    left_atoms = _skeleton_atoms(reactants)
    right_atoms = _skeleton_atoms(products)
    if any(
        left_atoms[left].GetAtomicNum() != right_atoms[right].GetAtomicNum()
        or left_atoms[left].GetIsotope() != right_atoms[right].GetIsotope()
        for left, right in mapping.items()
    ):
        return "mapped_atom_element_mismatch"
    return None


def _skeleton_fragment_index(side: _EndpointSide, skeleton_index: int) -> int:
    """Return the declared fragment containing a compact heavy-atom index."""
    explicit_index = side.skeleton_to_explicit[skeleton_index]
    return side.explicit_fragment[explicit_index]


def _fragments_are_exchangeable(side: _EndpointSide, left_fragment: int, right_fragment: int) -> bool:
    """Whether two separately declared fragments are truly interchangeable.

    Reactant order is bookkeeping, not atom identity.  Two identical HNO (or
    O2) molecules may exchange wholesale in an elementary reaction without
    creating a second chemical pathway.  This deliberately permits only an
    exact canonical molecular-graph match; similar fragments remain distinct.
    """
    if left_fragment == right_fragment:
        return True
    try:
        from rdkit import Chem

        left = Chem.MolToSmiles(side.skeleton_molecules[left_fragment], canonical=True, isomericSmiles=True)
        right = Chem.MolToSmiles(side.skeleton_molecules[right_fragment], canonical=True, isomericSmiles=True)
        return left == right
    except Exception:
        return False


def _skeleton_permutation_is_automorphism(
    side: _EndpointSide,
    permutation: dict[int, int],
) -> bool:
    """Check a conservative atom-identity symmetry of one endpoint.

    This is deliberately stricter than an element-only permutation.  It must
    preserve the complete heavy-atom graph, atomic electronic annotations, and
    the declared fragment boundary.  A boundary may be crossed only when the
    complete source and destination fragments are exact graph-identical
    molecules, which is ordinary reactant-exchange symmetry rather than a new
    reaction channel.
    """
    atom_count = side.skeleton_count
    if set(permutation) != set(range(atom_count)) or set(permutation.values()) != set(range(atom_count)):
        return False
    atoms = _skeleton_atoms(side)
    for source, target in permutation.items():
        left, right = atoms[source], atoms[target]
        source_fragment = _skeleton_fragment_index(side, source)
        target_fragment = _skeleton_fragment_index(side, target)
        if not _fragments_are_exchangeable(side, source_fragment, target_fragment):
            return False
        if (
            left.GetAtomicNum() != right.GetAtomicNum()
            or left.GetIsotope() != right.GetIsotope()
            or left.GetFormalCharge() != right.GetFormalCharge()
            or left.GetNumRadicalElectrons() != right.GetNumRadicalElectrons()
            or left.GetIsAromatic() != right.GetIsAromatic()
        ):
            return False
    bonds = _skeleton_bonds(side)
    permuted_bonds = {
        tuple(sorted((permutation[left], permutation[right]))): order
        for (left, right), order in bonds.items()
    }
    return bonds == permuted_bonds


def _maps_are_symmetry_equivalent(
    reactants: _EndpointSide,
    products: _EndpointSide,
    reference: dict[int, int],
    candidate: dict[int, int],
) -> bool:
    """Whether two maps differ only by a true endpoint graph symmetry.

    Either the reactant or product endpoint may supply the symmetry.  For
    example, swapping the two oxygen atoms of a single O2 reactant before a
    hydrogen abstraction does not define a second chemical path.
    """
    if set(reference) != set(candidate) or set(reference.values()) != set(candidate.values()):
        return False
    inverse_candidate = {product: reactant for reactant, product in candidate.items()}
    reactant_permutation = {
        reactant: inverse_candidate[product]
        for reactant, product in reference.items()
    }
    if _skeleton_permutation_is_automorphism(reactants, reactant_permutation):
        return True
    inverse_reference = {product: reactant for reactant, product in reference.items()}
    product_permutation = {
        product: candidate[inverse_reference[product]]
        for product in inverse_reference
    }
    return _skeleton_permutation_is_automorphism(products, product_permutation)


def _mapped_smiles_numbers(reactants: _EndpointSide, products: _EndpointSide) -> dict[int, int] | None:
    left_atoms = _skeleton_atoms(reactants)
    right_atoms = _skeleton_atoms(products)
    left = [atom.GetAtomMapNum() for atom in left_atoms]
    right = [atom.GetAtomMapNum() for atom in right_atoms]
    if not left or not all(left) or not all(right):
        return None
    if len(set(left)) != len(left) or set(left) != set(right):
        return None
    product_by_number = {number: index for index, number in enumerate(right)}
    return {index: product_by_number[number] for index, number in enumerate(left)}


def _automatic_skeleton_mapping(
    reactants: _EndpointSide,
    products: _EndpointSide,
    *,
    max_search_states: int,
) -> dict:
    """Find a unique minimum graph-edit mapping with a bounded exact search.

    The search is deliberately conservative.  It returns ambiguity or search
    exhaustion rather than resolving atom identity by RDKit's arbitrary
    canonical ordering.
    """
    left_atoms = _skeleton_atoms(reactants)
    right_atoms = _skeleton_atoms(products)
    if len(left_atoms) != len(right_atoms):
        return {"status": "unbalanced_smiles_atom_count", "valid": False}
    left_shells = {index: _neighbour_shells(reactants, index) for index in range(len(left_atoms))}
    right_shells = {index: _neighbour_shells(products, index) for index in range(len(right_atoms))}

    def local_cost(left_index: int, right_index: int) -> float:
        left, right = left_atoms[left_index], right_atoms[right_index]
        if left.GetAtomicNum() != right.GetAtomicNum() or left.GetIsotope() != right.GetIsotope():
            return math.inf
        cost = 5.0 * abs(left.GetDegree() - right.GetDegree())
        cost += 6.0 * int(left.GetIsAromatic() != right.GetIsAromatic())
        cost += 4.0 * int(left.IsInRing() != right.IsInRing())
        cost += 2.0 * abs(left.GetFormalCharge() - right.GetFormalCharge())
        cost += abs(left.GetNumRadicalElectrons() - right.GetNumRadicalElectrons())
        cost += 0.5 * abs(left.GetTotalNumHs(includeNeighbors=True) - right.GetTotalNumHs(includeNeighbors=True))
        for radius, (first, second) in enumerate(zip(left_shells[left_index], right_shells[right_index]), start=1):
            cost += (4.0 / radius) * _counter_distance(first, second)
        return cost

    candidates: dict[int, list[tuple[float, int]]] = {}
    for left in range(len(left_atoms)):
        scored = [
            (local_cost(left, right), right)
            for right in range(len(right_atoms))
        ]
        scored = [(cost, right) for cost, right in scored if math.isfinite(cost)]
        if not scored:
            return {"status": "unbalanced_elements", "valid": False}
        # Do not prune a chemically possible correspondence by a heuristic
        # local score: that could make a wrong mapping appear unique.  Large
        # same-element search spaces fail closed and require supplied atom maps.
        retained = sorted(scored)
        if len(retained) > 14:
            return {
                "status": "mapping_search_exhausted",
                "valid": False,
                "reason": f"atom {left} has {len(retained)} indistinguishable candidates",
            }
        candidates[left] = retained
    order = sorted(
        range(len(left_atoms)),
        key=lambda index: (len(candidates[index]), -left_atoms[index].GetDegree(), index),
    )
    left_bonds = _skeleton_bonds(reactants)
    right_bonds = _skeleton_bonds(products)
    best_cost = math.inf
    # We need every equal-cost candidate to distinguish a harmless endpoint
    # symmetry from genuinely different reaction channels.  The cap is a
    # safety bound: if the equivalence class itself is too large, require an
    # explicit map rather than silently choosing from a combinatorial space.
    max_equal_best_mappings = 64
    best: list[dict[int, int]] = []
    best_overflow = False
    states = 0
    exhausted = False

    def visit(position: int, mapping: dict[int, int], used: set[int], cost: float) -> None:
        nonlocal best_cost, best, best_overflow, states, exhausted
        if exhausted:
            return
        states += 1
        if states > max_search_states:
            exhausted = True
            return
        if cost > best_cost + 1e-9:
            return
        if position == len(order):
            if cost < best_cost - 1e-9:
                best_cost, best = cost, [dict(mapping)]
            elif abs(cost - best_cost) <= 1e-9:
                if len(best) < max_equal_best_mappings:
                    best.append(dict(mapping))
                else:
                    best_overflow = True
            return
        left = order[position]
        for atom_cost, right in candidates[left]:
            if right in used:
                continue
            incremental = atom_cost
            for prior_left, prior_right in mapping.items():
                left_bond = left_bonds.get(tuple(sorted((left, prior_left))))
                right_bond = right_bonds.get(tuple(sorted((right, prior_right))))
                if (left_bond is None) != (right_bond is None):
                    incremental += 18.0
                elif left_bond is not None and abs(left_bond - right_bond) > 0.1:
                    incremental += 7.0
            if cost + incremental <= best_cost + 1e-9:
                mapping[left] = right
                used.add(right)
                visit(position + 1, mapping, used, cost + incremental)
                used.remove(right)
                del mapping[left]

    visit(0, {}, set(), 0.0)
    if exhausted:
        return {
            "status": "mapping_search_exhausted",
            "valid": False,
            "search_states": states,
            "max_search_states": max_search_states,
        }
    if not best:
        return {"status": "atom_mapping_unresolved", "valid": False, "search_states": states}
    if best_overflow:
        return {
            "status": "mapping_symmetry_class_exhausted",
            "valid": False,
            "mapping_cost": best_cost,
            "search_states": states,
            "retained_equal_cost_mappings": len(best),
            "maximum_retained_equal_cost_mappings": max_equal_best_mappings,
        }
    if len(best) > 1:
        canonical = min(best, key=lambda item: tuple(item[index] for index in range(len(left_atoms))))
        if all(_maps_are_symmetry_equivalent(reactants, products, canonical, item) for item in best):
            return {
                "status": "resolved",
                "valid": True,
                "source": "automatic_minimum_graph_edit_symmetry_canonicalized",
                "mapping_cost": best_cost,
                "search_states": states,
                "mapping": canonical,
                "symmetry_resolution": {
                    "status": "equivalent_endpoint_symmetry_canonicalized",
                    "equivalent_mapping_count": len(best),
                    "canonical_mapping": sorted([left, right] for left, right in canonical.items()),
                    "equivalent_mappings": [
                        sorted([left, right] for left, right in item.items()) for item in best
                    ],
                },
            }
        return {
            "status": "ambiguous_atom_mapping",
            "valid": False,
            "mapping_cost": best_cost,
            "search_states": states,
            "candidate_mappings": [sorted([left, right] for left, right in item.items()) for item in best],
        }
    return {
        "status": "resolved",
        "valid": True,
        "source": "automatic_minimum_graph_edit",
        "mapping_cost": best_cost,
        "search_states": states,
        "mapping": best[0],
    }


def _extend_hydrogen_mapping(
    reactants: _EndpointSide,
    products: _EndpointSide,
    skeleton_mapping: dict[int, int],
) -> dict:
    """Extend a skeleton map across retained and transferred explicit hydrogens."""
    all_mapping = {
        reactants.skeleton_to_explicit[left]: products.skeleton_to_explicit[right]
        for left, right in skeleton_mapping.items()
    }
    consumed_left: set[int] = set(all_mapping)
    consumed_right: set[int] = set(all_mapping.values())
    transfers: list[dict] = []
    for left_skeleton, right_skeleton in skeleton_mapping.items():
        left_parent = reactants.skeleton_to_explicit[left_skeleton]
        right_parent = products.skeleton_to_explicit[right_skeleton]
        left_h = [value for value in reactants.attached_hydrogens.get(left_parent, []) if value not in consumed_left]
        right_h = [value for value in products.attached_hydrogens.get(right_parent, []) if value not in consumed_right]
        for isotope in sorted({reactants.isotopes[value] for value in left_h}
                              | {products.isotopes[value] for value in right_h}):
            left_isotope = sorted(value for value in left_h if reactants.isotopes[value] == isotope)
            right_isotope = sorted(value for value in right_h if products.isotopes[value] == isotope)
            retained = min(len(left_isotope), len(right_isotope))
            for left_hydrogen, right_hydrogen in zip(left_isotope[:retained], right_isotope[:retained]):
                all_mapping[left_hydrogen] = right_hydrogen
                consumed_left.add(left_hydrogen)
                consumed_right.add(right_hydrogen)

    # Collect every remaining hydrogen, including explicit standalone hydrogen
    # fragments that have no heavy-atom parent on that endpoint.
    unmatched_left = [
        index for index, element in enumerate(reactants.elements)
        if element == "H" and index not in consumed_left
    ]
    unmatched_right = [
        index for index, element in enumerate(products.elements)
        if element == "H" and index not in consumed_right
    ]
    unmatched_by_left_parent: dict[int | None, list[int]] = {}
    unmatched_by_right_parent: dict[int | None, list[int]] = {}
    for value in unmatched_left:
        unmatched_by_left_parent.setdefault(reactants.hydrogen_parent.get(value), []).append(value)
    for value in unmatched_right:
        unmatched_by_right_parent.setdefault(products.hydrogen_parent.get(value), []).append(value)
    if len(unmatched_left) != len(unmatched_right):
        return {"status": "unbalanced_explicit_hydrogens", "valid": False}
    # A transfer is chemically ambiguous only when the same isotope can be
    # assigned between multiple distinct donor *and* acceptor centres.  When
    # one side has a single centre, permutations on that centre are endpoint
    # symmetry (for example 2 HNO -> N2O + H2O) and may be canonicalized.
    ambiguous_isotopes = []
    for isotope in sorted({reactants.isotopes[value] for value in unmatched_left}
                          | {products.isotopes[value] for value in unmatched_right}):
        donor_parents = {
            reactants.hydrogen_parent.get(value)
            for value in unmatched_left if reactants.isotopes[value] == isotope
        }
        acceptor_parents = {
            products.hydrogen_parent.get(value)
            for value in unmatched_right if products.isotopes[value] == isotope
        }
        if len(donor_parents) > 1 and len(acceptor_parents) > 1:
            ambiguous_isotopes.append(isotope)
    if ambiguous_isotopes:
        return {
            "status": "ambiguous_hydrogen_transfer_mapping",
            "valid": False,
            "ambiguous_isotopes": ambiguous_isotopes,
            "reactant_donor_atom_indices": [
                value for value in sorted(
                    unmatched_by_left_parent, key=lambda item: (-1 if item is None else item)
                )
            ],
            "product_acceptor_atom_indices": [
                value for value in sorted(
                    unmatched_by_right_parent, key=lambda item: (-1 if item is None else item)
                )
            ],
        }
    symmetry_canonicalized = bool(unmatched_left) and (
        len(unmatched_by_left_parent) > 1 or len(unmatched_by_right_parent) > 1
    )
    for isotope in sorted({reactants.isotopes[value] for value in unmatched_left}
                          | {products.isotopes[value] for value in unmatched_right}):
        left_isotope = sorted(value for value in unmatched_left if reactants.isotopes[value] == isotope)
        right_isotope = sorted(value for value in unmatched_right if products.isotopes[value] == isotope)
        if len(left_isotope) != len(right_isotope):
            return {"status": "unbalanced_hydrogen_isotopes", "valid": False, "isotope": isotope}
        for left_hydrogen, right_hydrogen in zip(left_isotope, right_isotope):
            all_mapping[left_hydrogen] = right_hydrogen
            donor = reactants.hydrogen_parent.get(left_hydrogen)
            acceptor = products.hydrogen_parent.get(right_hydrogen)
            transfers.append({
                "reactant_hydrogen_index": left_hydrogen,
                "product_hydrogen_index": right_hydrogen,
                "reactant_donor_index": donor,
                "product_acceptor_index": acceptor,
                "reactant_fragment_index": reactants.explicit_fragment[left_hydrogen],
                "product_fragment_index": products.explicit_fragment[right_hydrogen],
                "isotope": isotope,
            })
    if set(all_mapping) != set(range(len(reactants.elements))) or set(all_mapping.values()) != set(range(len(products.elements))):
        return {"status": "incomplete_all_atom_mapping", "valid": False}
    if any(
        reactants.elements[left] != products.elements[right]
        or reactants.isotopes[left] != products.isotopes[right]
        for left, right in all_mapping.items()
    ):
        return {"status": "all_atom_element_mismatch", "valid": False}
    result = {"status": "resolved", "valid": True, "mapping": all_mapping, "hydrogen_transfers": transfers}
    if symmetry_canonicalized:
        result["hydrogen_symmetry_resolution"] = {
            "status": "equivalent_hydrogen_permutations_canonicalized",
            "reactant_donor_atom_indices": sorted(
                unmatched_by_left_parent, key=lambda item: (-1 if item is None else item)
            ),
            "product_acceptor_atom_indices": sorted(
                unmatched_by_right_parent, key=lambda item: (-1 if item is None else item)
            ),
        }
    return result


def _coupled_multiplicities(multiplicities: tuple[int, ...]) -> list[int]:
    if not multiplicities or any(value < 1 for value in multiplicities):
        return []
    twice_spins = {0}
    for multiplicity in multiplicities:
        fragment_twice_spin = multiplicity - 1
        twice_spins = {
            coupled
            for prior in twice_spins
            for coupled in range(abs(prior - fragment_twice_spin), prior + fragment_twice_spin + 1, 2)
        }
    return sorted(value + 1 for value in twice_spins)


def resolve_route_surface(route: RouteSpec, reactants: _EndpointSide | None = None,
                          products: _EndpointSide | None = None) -> dict:
    """Resolve charge and one explicitly common total-spin surface."""
    reactants = reactants or _build_endpoint_side(route.reactant_smiles)
    products = products or _build_endpoint_side(route.product_smiles)
    reactant_charges = route.reactant_charges or reactants.formal_charges
    product_charges = route.product_charges or products.formal_charges
    if len(reactant_charges) != len(route.reactant_smiles) or len(product_charges) != len(route.product_smiles):
        return {"status": "fragment_charge_contract_incomplete", "valid": False}
    if sum(reactant_charges) != sum(product_charges):
        return {"status": "unbalanced_total_charge", "valid": False,
                "reactant_total_charge": sum(reactant_charges), "product_total_charge": sum(product_charges)}
    if sum(reactant_charges) != route.charge:
        return {"status": "declared_total_charge_mismatch", "valid": False,
                "declared_total_charge": route.charge, "endpoint_total_charge": sum(reactant_charges)}
    reactant_multiplicities = route.reactant_multiplicities or reactants.inferred_multiplicities
    product_multiplicities = route.product_multiplicities or products.inferred_multiplicities
    if len(reactant_multiplicities) != len(route.reactant_smiles) or len(product_multiplicities) != len(route.product_smiles):
        return {"status": "fragment_spin_contract_incomplete", "valid": False}
    reactant_surfaces = _coupled_multiplicities(tuple(reactant_multiplicities))
    product_surfaces = _coupled_multiplicities(tuple(product_multiplicities))
    common = sorted(set(reactant_surfaces) & set(product_surfaces))
    if route.surface_multiplicities:
        common = sorted(set(common) & set(route.surface_multiplicities))
    if not common:
        return {
            "status": "no_common_adiabatic_spin_surface",
            "valid": False,
            "reactant_total_multiplicities": reactant_surfaces,
            "product_total_multiplicities": product_surfaces,
        }
    if route.multiplicity < 1:
        if len(common) != 1:
            return {
                "status": "multiple_common_spin_surfaces_require_selection",
                "valid": False,
                "common_total_multiplicities": common,
            }
        selected_multiplicity = common[0]
        selected_provenance = "unique_common_surface_from_fragment_states"
    elif route.multiplicity not in common:
        return {
            "status": "declared_multiplicity_not_common",
            "valid": False,
            "declared_multiplicity": route.multiplicity,
            "common_total_multiplicities": common,
        }
    else:
        selected_multiplicity = route.multiplicity
        selected_provenance = "RouteSpec.multiplicity"
    return {
        "status": "selected_common_surface",
        "valid": True,
        "charge": route.charge,
        "selected_multiplicity": selected_multiplicity,
        "reactant_fragment_charges": list(reactant_charges),
        "product_fragment_charges": list(product_charges),
        "reactant_fragment_multiplicities": list(reactant_multiplicities),
        "product_fragment_multiplicities": list(product_multiplicities),
        "reactant_total_multiplicities": reactant_surfaces,
        "product_total_multiplicities": product_surfaces,
        "common_total_multiplicities": common,
        "unattempted_common_surfaces": [value for value in common if value != selected_multiplicity],
        "state_provenance": {
            "charges": "route" if route.reactant_charges or route.product_charges else "SMILES_formal_charge",
            "multiplicities": "route" if route.reactant_multiplicities or route.product_multiplicities else "SMILES_radical_count",
            "selected_total_surface": selected_provenance,
        },
    }


def resolve_route_atom_mapping(route: RouteSpec, *, max_search_states: int = 200_000) -> dict:
    """Resolve a complete all-atom correspondence or return a fail-closed status."""
    if max_search_states < 1:
        raise ValueError("max_search_states must be positive")
    try:
        reactants = _build_endpoint_side(route.reactant_smiles)
        products = _build_endpoint_side(route.product_smiles)
    except ValueError as error:
        return {"status": "invalid_endpoint", "valid": False, "reason": str(error)}
    reactant_nuclei = Counter(zip(reactants.elements, reactants.isotopes))
    product_nuclei = Counter(zip(products.elements, products.isotopes))
    if reactant_nuclei != product_nuclei:
        return {
            "status": "unbalanced_elements",
            "valid": False,
            "reactant_elements": {
                f"{isotope or ''}{element}": count
                for (element, isotope), count in sorted(reactant_nuclei.items())
            },
            "product_elements": {
                f"{isotope or ''}{element}": count
                for (element, isotope), count in sorted(product_nuclei.items())
            },
        }
    surface = resolve_route_surface(route, reactants, products)
    if not surface["valid"]:
        return {**surface, "surface": surface}
    source = None
    if route.atom_mapping:
        skeleton_mapping = dict(route.atom_mapping)
        source = "RouteSpec.atom_mapping"
    else:
        skeleton_mapping = _mapped_smiles_numbers(reactants, products)
        source = "SMILES_atom_map_numbers" if skeleton_mapping is not None else None
    if skeleton_mapping is None:
        automatic = _automatic_skeleton_mapping(
            reactants, products, max_search_states=max_search_states,
        )
        if not automatic["valid"]:
            return {**automatic, "surface": surface}
        skeleton_mapping = automatic.pop("mapping")
        source = automatic.get("source")
        automatic_metadata = automatic
    else:
        automatic_metadata = {}
    mapping_error = _validate_skeleton_mapping(reactants, products, skeleton_mapping)
    if mapping_error:
        return {"status": mapping_error, "valid": False, "surface": surface, "source": source}
    extended = _extend_hydrogen_mapping(reactants, products, skeleton_mapping)
    if not extended["valid"]:
        return {**extended, "surface": surface, "source": source}
    all_mapping = extended["mapping"]
    return {
        "schema_version": 1,
        "kind": "resolved_reaction_atom_mapping",
        "status": "resolved",
        "valid": True,
        "source": source,
        "skeleton_atom_mapping": sorted([left, right] for left, right in skeleton_mapping.items()),
        "all_atom_mapping": sorted([left, right] for left, right in all_mapping.items()),
        "reactant_skeleton_atom_count": reactants.skeleton_count,
        "product_skeleton_atom_count": products.skeleton_count,
        "explicit_atom_count": len(reactants.elements),
        "hydrogen_transfers": extended["hydrogen_transfers"],
        "surface": surface,
        **({"hydrogen_symmetry_resolution": extended["hydrogen_symmetry_resolution"]}
           if "hydrogen_symmetry_resolution" in extended else {}),
        **{key: value for key, value in automatic_metadata.items() if key not in {"valid", "status", "source"}},
    }


def read_xyz(path: Path) -> tuple[list[str], np.ndarray, str]:
    lines = Path(path).read_text().splitlines()
    if not lines:
        raise ValueError(f"Empty XYZ file: {path}")
    count = int(lines[0])
    rows = [line.split() for line in lines[2:2 + count]]
    if len(rows) != count or any(len(row) < 4 for row in rows):
        raise ValueError(f"Incomplete XYZ geometry: {path}")
    return [row[0] for row in rows], np.asarray([[float(value) for value in row[1:4]] for row in rows]), (
        lines[1] if len(lines) > 1 else ""
    )


def write_xyz(path: Path, elements: list[str], coordinates: np.ndarray, comment: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(elements)), comment]
    lines.extend(
        f"{element:<2} {x: .10f} {y: .10f} {z: .10f}"
        for element, (x, y, z) in zip(elements, coordinates)
    )
    path.write_text("\n".join(lines) + "\n")
    return path


def _graph_component_count(atom_count: int, bonds: set[tuple[int, int]]) -> int:
    adjacency = {index: set() for index in range(atom_count)}
    for left, right in bonds:
        adjacency[left].add(right)
        adjacency[right].add(left)
    count = 0
    remaining = set(adjacency)
    while remaining:
        count += 1
        stack = [remaining.pop()]
        while stack:
            neighbours = adjacency[stack.pop()] & remaining
            stack.extend(neighbours)
            remaining -= neighbours
    return count


def validate_declared_connectivity(
    xyz_path: Path,
    smiles: tuple[str, ...],
    *,
    expected_to_xyz: dict[int, int] | None = None,
    expected_bond_max_covalent_ratio: float = 1.45,
    unexpected_bond_ratio: float = 1.20,
    ambiguous_contact_ratio: float = 1.45,
) -> dict:
    """Check that an optimized XYZ still represents a declared fragment graph.

    XYZ files contain no bond table, so this is deliberately a conservative
    geometry gate rather than a bond-order assignment.  Declared bonds must
    remain within a liberal covalent-radius envelope, while undeclared pairs
    inside the covalent or near-covalent envelope fail closed.  Atom order is
    explicit: ``expected_to_xyz`` maps the declared explicit-H order to the XYZ
    order used by an endpoint complex.
    """
    from rdkit import Chem

    side = _build_endpoint_side(tuple(smiles))
    elements, coordinates, _ = read_xyz(Path(xyz_path))
    atom_count = len(side.elements)
    if len(elements) != atom_count:
        return {
            "status": "atom_count_mismatch", "valid": False,
            "declared_atom_count": atom_count, "xyz_atom_count": len(elements),
        }
    mapping = dict(expected_to_xyz or {index: index for index in range(atom_count)})
    if set(mapping) != set(range(atom_count)) or set(mapping.values()) != set(range(atom_count)):
        return {
            "status": "invalid_expected_to_xyz_mapping", "valid": False,
            "declared_atom_count": atom_count,
        }
    expected_elements = [""] * atom_count
    for declared_index, xyz_index in mapping.items():
        expected_elements[xyz_index] = side.elements[declared_index]
    if elements != expected_elements:
        return {
            "status": "element_order_mismatch", "valid": False,
            "declared_elements_in_xyz_order": expected_elements,
            "xyz_elements": elements,
        }
    expected_bonds = {
        tuple(sorted((mapping[left], mapping[right])))
        for left, right in side.bonds
    }
    periodic_table = Chem.GetPeriodicTable()
    missing_bonds: list[dict] = []
    unexpected_bonds: list[dict] = []
    ambiguous_contacts: list[dict] = []
    implausible_overlaps: list[dict] = []
    inferred_bonds: set[tuple[int, int]] = set()
    for left in range(atom_count):
        for right in range(left + 1, atom_count):
            pair = (left, right)
            distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
            covalent_sum = float(
                periodic_table.GetRcovalent(periodic_table.GetAtomicNumber(elements[left]))
                + periodic_table.GetRcovalent(periodic_table.GetAtomicNumber(elements[right]))
            )
            if covalent_sum <= 0.0:
                return {
                    "status": "covalent_radius_unavailable", "valid": False,
                    "atom_indices": [left, right], "elements": [elements[left], elements[right]],
                }
            ratio = distance / covalent_sum
            detail = {
                "atom_indices": [left, right],
                "elements": [elements[left], elements[right]],
                "distance_angstrom": distance,
                "covalent_radius_ratio": ratio,
            }
            if ratio < 0.45:
                implausible_overlaps.append(detail)
            if pair in expected_bonds:
                if ratio <= expected_bond_max_covalent_ratio:
                    inferred_bonds.add(pair)
                else:
                    missing_bonds.append(detail)
            elif ratio <= unexpected_bond_ratio:
                unexpected_bonds.append(detail)
                inferred_bonds.add(pair)
            elif ratio <= ambiguous_contact_ratio:
                # A near-covalent undeclared contact is not enough to assign a
                # bond, but it is also not safe evidence that the channel was
                # retained.  Fail closed and request better endpoint evidence.
                ambiguous_contacts.append(detail)
    if implausible_overlaps or missing_bonds or unexpected_bonds:
        status = "connectivity_changed"
    elif ambiguous_contacts:
        status = "connectivity_ambiguous"
    else:
        status = "validated"
    valid = status == "validated"
    return {
        "schema_version": 1,
        "kind": "declared_endpoint_connectivity_check",
        "status": status,
        "valid": valid,
        "xyz": str(Path(xyz_path)),
        "declared_smiles": list(smiles),
        "atom_order_aware": True,
        "bond_orders_checked": False,
        "declared_bonds": [list(pair) for pair in sorted(expected_bonds)],
        "inferred_covalent_bonds": [list(pair) for pair in sorted(inferred_bonds)],
        "declared_fragment_count": _graph_component_count(atom_count, expected_bonds),
        "inferred_fragment_count": _graph_component_count(atom_count, inferred_bonds),
        "missing_declared_bonds": missing_bonds,
        "unexpected_covalent_bonds": unexpected_bonds,
        "ambiguous_near_covalent_contacts": ambiguous_contacts,
        "implausible_overlaps": implausible_overlaps,
        "thresholds": {
            "expected_bond_max_covalent_ratio": expected_bond_max_covalent_ratio,
            "unexpected_bond_ratio": unexpected_bond_ratio,
            "ambiguous_contact_ratio": ambiguous_contact_ratio,
        },
    }


def validate_route_endpoint_connectivity(
    route: RouteSpec,
    xyz_path: Path,
    *,
    side: str,
    mapping_result: dict,
) -> dict:
    """Validate a mapped endpoint-complex XYZ against one side of a route."""
    if side == "reactants":
        smiles = route.reactant_smiles
        expected_to_xyz = None
    elif side == "products":
        smiles = route.product_smiles
        # Endpoint-complex product XYZ files are written in reactant all-atom
        # order, while the declared product graph is in product order.
        expected_to_xyz = {
            int(product): int(reactant)
            for reactant, product in mapping_result["all_atom_mapping"]
        }
    else:
        raise ValueError("side must be 'reactants' or 'products'")
    return validate_declared_connectivity(
        xyz_path, tuple(smiles), expected_to_xyz=expected_to_xyz,
    )


def validate_route_mapping(route: RouteSpec) -> dict:
    """Validate or uniquely resolve a route's complete atom/state contract."""
    result = resolve_route_atom_mapping(route)
    if not result.get("valid"):
        return result
    return {
        **result,
        "status": "validated",
        "skeleton_atom_count": result["reactant_skeleton_atom_count"],
        "atom_count_with_explicit_hydrogens": result["explicit_atom_count"],
    }


def _fragment_indices(parent_smiles: str, bond: tuple[int, int]) -> tuple[list[int], list[int], list[str]]:
    """Return explicit-H fragment indices after deleting one heavy-atom bond."""
    from rdkit import Chem

    molecule = Chem.AddHs(Chem.MolFromSmiles(parent_smiles))
    if molecule is None:
        raise ValueError(f"Invalid parent SMILES: {parent_smiles}")
    left, right = bond
    existing = molecule.GetBondBetweenAtoms(left, right)
    if existing is None:
        raise ValueError(f"Retained bond {bond} does not exist in the parent")
    editable = Chem.RWMol(molecule)
    editable.RemoveBond(left, right)
    components = [list(component) for component in Chem.GetMolFrags(editable.GetMol())]
    if len(components) != 2:
        raise ValueError("Dissociation v2 requires one bond that separates the parent into exactly two fragments")
    left_component = next(component for component in components if left in component)
    right_component = next(component for component in components if right in component)
    elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    return left_component, right_component, elements


def build_dissociation_endpoint(route: RouteSpec, reactant_xyz: Path, product_xyz: Path,
                                *, final_distance_angstrom: float = 3.5) -> dict:
    """Translate the two fragments while preserving atom order for ORCA scans."""
    if route.protocol != "relaxed_dissociation_scan" or len(route.broken_bonds) != 1:
        raise ValueError("Dissociation endpoint construction requires one broken-bond route")
    elements, coordinates, _ = read_xyz(reactant_xyz)
    left_component, right_component, expected_elements = _fragment_indices(
        route.parent_smiles, route.broken_bonds[0]
    )
    if elements != expected_elements:
        raise ValueError(
            "Reactant XYZ atom order does not match the explicit-H parent order; "
            "the endpoint cannot be mapped safely"
        )
    left, right = route.broken_bonds[0]
    vector = coordinates[right] - coordinates[left]
    distance = float(np.linalg.norm(vector))
    if distance < 1e-8:
        raise ValueError("Broken-bond atoms occupy the same coordinates")
    requested = max(final_distance_angstrom, distance + 0.5)
    displacement = (requested - distance) * vector / distance
    product_coordinates = coordinates.copy()
    # Split the displacement evenly to retain an approximately fixed centre.
    product_coordinates[left_component] -= displacement * 0.5
    product_coordinates[right_component] += displacement * 0.5
    write_xyz(
        product_xyz,
        elements,
        product_coordinates,
        f"STORCA mapped dissociation endpoint; bond {left}-{right} = {requested:.4f} Angstrom",
    )
    return {
        "status": "completed",
        "reactant_xyz": str(Path(reactant_xyz)),
        "product_xyz": str(Path(product_xyz)),
        "broken_bond": [left, right],
        "initial_distance_angstrom": distance,
        "final_distance_angstrom": requested,
        "atom_order_preserved": True,
        "fragment_atom_indices": [left_component, right_component],
    }


def interpolate_xyz(start_xyz: Path, end_xyz: Path, output: Path, *, frames: int = 24,
                    comment_prefix: str = "STORCA path preview") -> Path:
    """Create a geometry interpolation for candidate previews only."""
    if frames < 2:
        raise ValueError("A trajectory needs at least two frames")
    start_elements, start, _ = read_xyz(start_xyz)
    end_elements, end, _ = read_xyz(end_xyz)
    if start_elements != end_elements or start.shape != end.shape:
        raise ValueError("Trajectory endpoints must have identical atom order")
    blocks = []
    for index, fraction in enumerate(np.linspace(0.0, 1.0, frames)):
        coordinates = start + fraction * (end - start)
        blocks.append(
            "\n".join(
                [str(len(start_elements)), f"{comment_prefix}; frame={index + 1}; fraction={fraction:.6f}"]
                + [
                    f"{element:<2} {x: .10f} {y: .10f} {z: .10f}"
                    for element, (x, y, z) in zip(start_elements, coordinates)
                ]
            )
        )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(blocks) + "\n")
    return output


def _rotation_matrix(orientation: int, fragment_index: int, side_index: int) -> np.ndarray:
    """Return one of a deterministic set of proper 3-D rotations."""
    angles = (
        (0.0, 0.0, 0.0),
        (math.pi, 0.0, math.pi / 2.0),
        (math.pi / 2.0, math.pi, math.pi / 3.0),
        (3.0 * math.pi / 2.0, math.pi / 2.0, math.pi),
        (math.pi / 3.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0),
        (5.0 * math.pi / 3.0, math.pi / 3.0, 2.0 * math.pi / 3.0),
    )
    alpha, beta, gamma = angles[(orientation + 2 * fragment_index + side_index) % len(angles)]
    x = np.asarray([[1, 0, 0], [0, math.cos(alpha), -math.sin(alpha)],
                    [0, math.sin(alpha), math.cos(alpha)]], dtype=float)
    y = np.asarray([[math.cos(beta), 0, math.sin(beta)], [0, 1, 0],
                    [-math.sin(beta), 0, math.cos(beta)]], dtype=float)
    z = np.asarray([[math.cos(gamma), -math.sin(gamma), 0],
                    [math.sin(gamma), math.cos(gamma), 0], [0, 0, 1]], dtype=float)
    return z @ y @ x


def _orientation_direction(orientation: int, fragment_index: int) -> np.ndarray:
    directions = np.asarray([
        [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0],
        [1.0, 1.0, 1.0], [-1.0, 1.0, -1.0],
    ])
    value = directions[(orientation + fragment_index) % len(directions)]
    return value / np.linalg.norm(value)


def _embed_fragment(molecule: Any, *, random_seed: int) -> np.ndarray:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    embedded = Chem.Mol(molecule)
    # Distance geometry has no internal coordinates to solve for a monatomic
    # fragment.  Handle standalone H/H+/H- (and other atoms) directly so a
    # valid reaction endpoint does not depend on an RDKit embedding edge case.
    if embedded.GetNumAtoms() == 1:
        return np.zeros((1, 3), dtype=float)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = int(random_seed)
    parameters.useRandomCoords = False
    if AllChem.EmbedMolecule(embedded, parameters) != 0:
        # Random coordinates are slower but rescue linear/highly ionic species.
        parameters.useRandomCoords = True
        if AllChem.EmbedMolecule(embedded, parameters) != 0:
            raise RuntimeError("RDKit could not generate an endpoint fragment conformer")
    try:
        if AllChem.UFFHasAllMoleculeParams(embedded):
            AllChem.UFFOptimizeMolecule(embedded, maxIters=500)
    except Exception:
        # UFF is a seed refinement, not evidence; ORCA performs validation.
        pass
    coordinates = np.asarray(embedded.GetConformer().GetPositions(), dtype=float)
    return coordinates - coordinates.mean(axis=0)


def _explicit_fragment_offsets(side: _EndpointSide) -> list[int]:
    offsets, value = [], 0
    for molecule in side.explicit_molecules:
        offsets.append(value)
        value += molecule.GetNumAtoms()
    return offsets


def _pack_endpoint_fragments(
    side: _EndpointSide,
    *,
    orientation: int,
    side_index: int,
    contact_pairs: list[tuple[int, int]],
    encounter_distance_angstrom: float,
    fragment_coordinates: list[np.ndarray] | None = None,
) -> np.ndarray:
    """Pack fragments, prioritising topology-derived interfragment contacts."""
    offsets = _explicit_fragment_offsets(side)
    if fragment_coordinates is not None and len(fragment_coordinates) != len(side.explicit_molecules):
        raise ValueError("Supplied fragment coordinates do not match the declared endpoint stoichiometry")
    blocks = []
    for index, molecule in enumerate(side.explicit_molecules):
        if fragment_coordinates is None:
            coordinates = _embed_fragment(
                molecule, random_seed=1729 + 101 * side_index + 17 * orientation + index,
            )
        else:
            coordinates = np.asarray(fragment_coordinates[index], dtype=float).copy()
            if coordinates.shape != (molecule.GetNumAtoms(), 3) or not np.isfinite(coordinates).all():
                raise ValueError(f"Invalid coordinates for endpoint fragment {index + 1}")
            coordinates -= coordinates.mean(axis=0)
        blocks.append(coordinates @ _rotation_matrix(orientation, index, side_index).T)
    anchor = max(range(len(blocks)), key=lambda index: side.explicit_molecules[index].GetNumHeavyAtoms())
    placed: dict[int, np.ndarray] = {anchor: blocks[anchor]}
    remaining = [index for index in range(len(blocks)) if index != anchor]

    def local(global_index: int) -> tuple[int, int]:
        fragment = side.explicit_fragment[global_index]
        return fragment, global_index - offsets[fragment]

    while remaining:
        fragment = remaining.pop(0)
        relevant = []
        for left, right in contact_pairs:
            left_fragment, left_local = local(left)
            right_fragment, right_local = local(right)
            if left_fragment in placed and right_fragment == fragment:
                relevant.append((left_fragment, left_local, right_local))
            elif right_fragment in placed and left_fragment == fragment:
                relevant.append((right_fragment, right_local, left_local))
        direction = _orientation_direction(orientation, fragment)
        if relevant:
            placed_fragment, placed_atom, moving_atom = relevant[0]
            target = placed[placed_fragment][placed_atom] + direction * encounter_distance_angstrom
            moved = blocks[fragment] + (target - blocks[fragment][moving_atom])
        else:
            anchor_radius = max(np.linalg.norm(value, axis=1).max(initial=0.0) for value in placed.values())
            radius = np.linalg.norm(blocks[fragment], axis=1).max(initial=0.0)
            moved = blocks[fragment] + direction * (anchor_radius + radius + 1.5)
        # A contact atom can point the rest of a randomly embedded fragment
        # through the placed complex.  Translate outward until no atom pair is
        # pathologically close, retaining the same deterministic orientation.
        for _ in range(12):
            shortest_to_placed = min(
                float(np.linalg.norm(left - right))
                for previous in placed.values()
                for left in moved
                for right in previous
            )
            if shortest_to_placed >= 0.75:
                break
            moved += direction * (0.85 - shortest_to_placed)
        placed[fragment] = moved
    coordinates = np.vstack([placed[index] for index in range(len(blocks))])
    # Fail rather than hand ORCA a malformed complex with overlapping fragments.
    shortest = math.inf
    for left in range(len(coordinates)):
        for right in range(left + 1, len(coordinates)):
            if side.explicit_fragment[left] != side.explicit_fragment[right]:
                shortest = min(shortest, float(np.linalg.norm(coordinates[left] - coordinates[right])))
    if shortest < 0.65:
        raise RuntimeError(f"Endpoint packing produced an interfragment contact of {shortest:.3f} Å")
    return coordinates


def _fragment_coordinates_from_xyz(
    side: _EndpointSide,
    paths: list[str | Path],
) -> list[np.ndarray]:
    """Load atom-order-preserving geometries for every declared fragment."""
    if len(paths) != len(side.explicit_molecules):
        raise ValueError("Fragment XYZ count does not match the declared endpoint stoichiometry")
    result = []
    for index, (molecule, path) in enumerate(zip(side.explicit_molecules, paths), start=1):
        elements, coordinates, _ = read_xyz(Path(path))
        expected = [atom.GetSymbol() for atom in molecule.GetAtoms()]
        if elements != expected:
            raise ValueError(
                f"Fragment {index} XYZ changed atom composition/order: expected {expected}, found {elements}"
            )
        result.append(coordinates)
    return result


def _mapped_bond_changes(
    reactants: _EndpointSide,
    products: _EndpointSide,
    all_mapping: dict[int, int],
) -> dict:
    inverse = {right: left for left, right in all_mapping.items()}
    product_in_reactant_order = {
        tuple(sorted((inverse[left], inverse[right]))): order
        for (left, right), order in products.bonds.items()
    }
    broken, formed, changed_order = [], [], []
    for pair in sorted(set(reactants.bonds) | set(product_in_reactant_order)):
        reactant_order = reactants.bonds.get(pair)
        product_order = product_in_reactant_order.get(pair)
        if reactant_order is not None and product_order is None:
            broken.append(list(pair))
        elif reactant_order is None and product_order is not None:
            formed.append(list(pair))
        elif abs(float(reactant_order) - float(product_order)) > 0.1:
            changed_order.append({"atom_indices": list(pair), "reactant_order": reactant_order,
                                  "product_order": product_order})
    return {"broken_bonds": broken, "formed_bonds": formed, "bond_order_changes": changed_order}


def _largest_unchanged_component(atom_count: int, changes: dict,
                                 bonds: dict[tuple[int, int], float]) -> list[int]:
    changed = {
        atom for pair in [*changes["broken_bonds"], *changes["formed_bonds"]] for atom in pair
    }
    changed.update(atom for item in changes["bond_order_changes"] for atom in item["atom_indices"])
    adjacency = {index: set() for index in range(atom_count) if index not in changed}
    for left, right in bonds:
        if left in adjacency and right in adjacency:
            adjacency[left].add(right)
            adjacency[right].add(left)
    components: list[list[int]] = []
    remaining = set(adjacency)
    while remaining:
        seed = min(remaining)
        stack, component = [seed], []
        remaining.remove(seed)
        while stack:
            current = stack.pop()
            component.append(current)
            neighbours = adjacency[current] & remaining
            stack.extend(sorted(neighbours, reverse=True))
            remaining -= neighbours
        components.append(sorted(component))
    return max(components, key=lambda value: (len(value), [-item for item in value]), default=[])


def _align_product_core(reference: np.ndarray, candidate: np.ndarray,
                        core_indices: list[int]) -> np.ndarray:
    if not core_indices:
        return candidate
    reference_core = reference[core_indices]
    candidate_core = candidate[core_indices]
    reference_center = reference_core.mean(axis=0)
    candidate_center = candidate_core.mean(axis=0)
    if len(core_indices) < 2:
        return candidate - candidate_center + reference_center
    covariance = (candidate_core - candidate_center).T @ (reference_core - reference_center)
    left, _, right = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right))
    rotation = left @ correction @ right
    return (candidate - candidate_center) @ rotation + reference_center


def build_endpoint_complex_seeds(
    route: RouteSpec,
    output_dir: Path,
    *,
    orientations: int = 3,
    encounter_distance_angstrom: float = 2.8,
    mapping_result: dict | None = None,
    fragment_xyz: dict[str, list[str | Path]] | None = None,
) -> dict:
    """Build atom-consistent endpoint complexes for any resolved stoichiometry.

    Multiple deterministic orientations are generated whenever either endpoint
    has more than one fragment.  They are geometry seeds, never reaction-path
    evidence; ORCA optimization/frequency and TS/IRC validation remain required.
    """
    if orientations < 1 or orientations > 12:
        raise ValueError("Endpoint construction supports between 1 and 12 orientations")
    if not math.isfinite(encounter_distance_angstrom) or encounter_distance_angstrom <= 1.0:
        raise ValueError("Encounter distance must be finite and greater than 1 Å")
    mapping_result = mapping_result or resolve_route_atom_mapping(route)
    if not mapping_result.get("valid"):
        raise ValueError(f"Cannot build endpoint complexes: {mapping_result.get('status')}")
    reactants = _build_endpoint_side(route.reactant_smiles)
    products = _build_endpoint_side(route.product_smiles)
    reactant_fragment_coordinates = None
    product_fragment_coordinates = None
    if fragment_xyz is not None:
        reactant_fragment_coordinates = _fragment_coordinates_from_xyz(
            reactants, list(fragment_xyz.get("reactants") or []),
        )
        product_fragment_coordinates = _fragment_coordinates_from_xyz(
            products, list(fragment_xyz.get("products") or []),
        )
    all_mapping = {int(left): int(right) for left, right in mapping_result["all_atom_mapping"]}
    inverse = {right: left for left, right in all_mapping.items()}
    changes = _mapped_bond_changes(reactants, products, all_mapping)
    product_formed = [tuple(all_mapping[index] for index in pair) for pair in changes["broken_bonds"]]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    rejected_orientations = []
    count = orientations if max(len(route.reactant_smiles), len(route.product_smiles)) > 1 else 1
    # Encounter orientations are seeds, not evidence.  Try additional
    # deterministic rotations and gradually widen the encounter if a seed
    # creates an undeclared covalent contact.  Invalid seed geometry must not
    # consume one of the independent ORCA orientations required for a
    # barrierless classification.
    maximum_attempts = max(count * 8, 8)
    for attempt in range(maximum_attempts):
        orientation = attempt
        encounter_distance = encounter_distance_angstrom + 0.35 * (attempt // count)
        folder = output_dir / f"candidate-{attempt + 1:02d}"
        folder.mkdir(parents=True, exist_ok=True)
        reactant_coordinates = _pack_endpoint_fragments(
            reactants,
            orientation=orientation,
            side_index=0,
            contact_pairs=[tuple(pair) for pair in changes["formed_bonds"]],
            encounter_distance_angstrom=encounter_distance,
            fragment_coordinates=reactant_fragment_coordinates,
        )
        raw_product = _pack_endpoint_fragments(
            products,
            orientation=orientation,
            side_index=1,
            contact_pairs=product_formed,
            encounter_distance_angstrom=encounter_distance,
            fragment_coordinates=product_fragment_coordinates,
        )
        product_coordinates = np.zeros_like(reactant_coordinates)
        for reactant_index, product_index in all_mapping.items():
            product_coordinates[reactant_index] = raw_product[product_index]
        core = _largest_unchanged_component(len(reactant_coordinates), changes, reactants.bonds)
        product_coordinates = _align_product_core(reactant_coordinates, product_coordinates, core)
        reactant_xyz = write_xyz(
            folder / "reactant-seed.xyz", reactants.elements, reactant_coordinates,
            f"STORCA endpoint-complex seed; route={route.route_id}; orientation={orientation + 1}",
        )
        product_xyz = write_xyz(
            folder / "product-seed.xyz", reactants.elements, product_coordinates,
            f"STORCA mapped endpoint-complex seed; route={route.route_id}; candidate={attempt + 1}",
        )
        reactant_connectivity = validate_route_endpoint_connectivity(
            route, reactant_xyz, side="reactants", mapping_result=mapping_result,
        )
        product_connectivity = validate_route_endpoint_connectivity(
            route, product_xyz, side="products", mapping_result=mapping_result,
        )
        if not (reactant_connectivity.get("valid") and product_connectivity.get("valid")):
            rejected_orientations.append({
                "candidate": attempt + 1,
                "encounter_distance_angstrom": encounter_distance,
                "reactant_connectivity": reactant_connectivity,
                "product_connectivity": product_connectivity,
            })
            continue
        records.append({
            "orientation": len(records) + 1,
            "generation_candidate": attempt + 1,
            "encounter_distance_angstrom": encounter_distance,
            "reactant_xyz": str(reactant_xyz),
            "product_xyz": str(product_xyz),
            "alignment_core_atom_indices": core,
        })
        if len(records) >= count:
            break
    result = {
        "schema_version": 1,
        "kind": "reaction_endpoint_complex_seeds",
        "status": "completed",
        "route_id": route.route_id,
        "orientation_count": len(records),
        "requested_orientation_count": count,
        "generation_attempt_count": min(maximum_attempts, attempt + 1),
        "rejected_orientations": rejected_orientations,
        "orientations": records,
        "mapping": mapping_result,
        "bond_changes": changes,
        "encounter_distance_angstrom": encounter_distance_angstrom,
        "fragment_geometry_source": (
            "supplied_atom_ordered_fragment_xyz" if fragment_xyz is not None else "rdkit_seed_conformers"
        ),
        "interpretation": "Deterministic endpoint seeds; not optimized minima, transition states, or dynamics.",
    }
    write_geometry_manifest(route, result, output_dir / "endpoint-complexes.json")
    return result


def build_mapped_endpoint_seeds(route: RouteSpec, output_dir: Path) -> dict:
    """Backward-compatible one-orientation wrapper around the generic builder."""
    result = build_endpoint_complex_seeds(route, output_dir, orientations=1)
    orientation = result["orientations"][0]
    return {
        "status": "completed",
        "reactant_xyz": orientation["reactant_xyz"],
        "product_xyz": orientation["product_xyz"],
        "all_atom_mapping": result["mapping"]["all_atom_mapping"],
        "mapping": result["mapping"],
        "bond_changes": result["bond_changes"],
    }


def aligned_rmsd(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Kabsch-aligned RMSD without changing either input."""
    if reference.shape != candidate.shape:
        raise ValueError("RMSD arrays must have the same shape")
    ref = reference - reference.mean(axis=0)
    trial = candidate - candidate.mean(axis=0)
    covariance = trial.T @ ref
    left, _, right = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right))
    rotation = left @ correction @ right
    aligned = trial @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - ref) ** 2, axis=1))))


def validate_trajectory_endpoints(trajectory_xyz: Path, reactant_xyz: Path, product_xyz: Path,
                                  *, rmsd_tolerance_angstrom: float = 0.75) -> dict:
    """Require the two trajectory ends to pair with the declared endpoints."""
    from .decomposition_visuals import xyz_frames

    frames = xyz_frames(trajectory_xyz)
    reactant_elements, reactant, _ = read_xyz(reactant_xyz)
    product_elements, product, _ = read_xyz(product_xyz)
    if frames[0][0] != reactant_elements or frames[-1][0] != reactant_elements:
        return {"status": "element_order_mismatch", "valid": False}
    if product_elements != reactant_elements:
        return {"status": "reference_element_order_mismatch", "valid": False}
    forward = aligned_rmsd(reactant, frames[0][1]) + aligned_rmsd(product, frames[-1][1])
    reverse = aligned_rmsd(reactant, frames[-1][1]) + aligned_rmsd(product, frames[0][1])
    selected = "forward" if forward <= reverse else "reverse"
    pair = forward if forward <= reverse else reverse
    return {
        "status": "validated" if pair <= 2 * rmsd_tolerance_angstrom else "endpoint_mismatch",
        "valid": pair <= 2 * rmsd_tolerance_angstrom,
        "orientation": selected,
        "paired_rmsd_sum_angstrom": pair,
        "per_endpoint_rmsd_tolerance_angstrom": rmsd_tolerance_angstrom,
    }


def write_geometry_manifest(route: RouteSpec, result: dict, output: Path) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({
        "schema_version": 1,
        "kind": "decomposition_route_geometry",
        "route": route.as_dict(),
        "result": result,
    }, indent=2, sort_keys=True) + "\n")
    return output
