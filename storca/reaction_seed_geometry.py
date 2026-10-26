"""Cheap, map-preserving *proposals* for paired reaction endpoints.

Nothing in this module runs an electronic calculation or establishes a minimum,
transition state, mechanism, population, symmetry number, or rate.  In particular,
an exhausted/rejected portfolio is not evidence that a reaction cannot occur.
All atom indices on both sides use one common, explicit-nucleus order.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


POLICY_ID = "paired_rigid_reaction_geometry_v1"
_RADII = {
    "H": .31, "He": .28, "Li": 1.28, "Be": .96, "B": .84, "C": .76,
    "N": .71, "O": .66, "F": .57, "Ne": .58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Br": 1.20, "I": 1.39,
}
_ELEMENTS = frozenset("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split())


class GeometryPreparationError(ValueError):
    """Malformed or inconsistent input, not a chemical nonexistence claim."""


class UnsupportedGeometryContext(GeometryPreparationError):
    """The declared context needs a different allowed-motion model."""


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(value: Any) -> str:
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _coordinates(value: Any, count: int) -> tuple[tuple[float, float, float], ...]:
    array = np.asarray(value, dtype=float)
    if array.shape != (count, 3) or not np.isfinite(array).all():
        raise GeometryPreparationError("Coordinates must be finite N by 3 arrays")
    return tuple(tuple(float(v) for v in row) for row in array)


def _index(value: Any, count: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise GeometryPreparationError("Atom indices must be integers")
    if not 0 <= int(value) < count:
        raise GeometryPreparationError("Atom index is outside the common nuclear order")
    return int(value)


def _bonds(value: Sequence[Any], count: int) -> tuple[tuple[int, int, float], ...]:
    bonds: dict[tuple[int, int], float] = {}
    for item in value:
        if isinstance(item, Mapping):
            left, right = item["atom_indices"]
            order = item.get("order", item.get("bond_order", 1.0))
        else:
            if len(item) not in (2, 3):
                raise GeometryPreparationError("Bonds must contain two indices and an optional order")
            left, right = item[:2]
            order = item[2] if len(item) == 3 else 1.0
        pair = tuple(sorted((_index(left, count), _index(right, count))))
        if pair[0] == pair[1] or pair in bonds:
            raise GeometryPreparationError("Self bonds and duplicate bonds are not allowed")
        order = float(order)
        if not math.isfinite(order) or order <= 0:
            raise GeometryPreparationError("Bond orders must be positive and finite")
        bonds[pair] = order
    return tuple((left, right, order) for (left, right), order in sorted(bonds.items()))


@dataclass(frozen=True)
class AtomIdentity:
    map_number: int
    element: str
    isotope: int = 0
    stereochemistry: str | None = None

    def to_dict(self) -> dict:
        return {"map_number": self.map_number, "element": self.element,
                "isotope": self.isotope, "stereochemistry": self.stereochemistry}


@dataclass(frozen=True)
class EndpointGeometry:
    coordinates: tuple[tuple[float, float, float], ...]
    fragment_ids: tuple[int, ...]
    bonds: tuple[tuple[int, int, float], ...]
    source_bindings_json: str = "[]"
    conformer_ids_json: str = "[]"
    fragment_states_json: str = "[]"

    def to_dict(self) -> dict:
        return {"coordinates": [list(row) for row in self.coordinates],
                "fragment_ids": list(self.fragment_ids),
                "bonds": [list(row) for row in self.bonds],
                "source_bindings": json.loads(self.source_bindings_json),
                "conformer_ids": json.loads(self.conformer_ids_json),
                "fragment_states": json.loads(self.fragment_states_json)}


@dataclass(frozen=True)
class ReactionGeometryContext:
    atoms: tuple[AtomIdentity, ...]
    reactants: EndpointGeometry
    products: EndpointGeometry
    charge: int
    multiplicity: int
    environment_json: str
    allowed_motions: str = "free_molecular_fragments"
    parent_evidence_json: str = "[]"

    @property
    def context_id(self) -> str:
        return _sha(self.to_dict())

    @property
    def environment(self) -> dict:
        return json.loads(self.environment_json)

    @property
    def bond_changes(self) -> dict:
        first = {(a, b): order for a, b, order in self.reactants.bonds}
        second = {(a, b): order for a, b, order in self.products.bonds}
        return {
            "broken_bonds": [list(pair) for pair in sorted(first.keys() - second.keys())],
            "formed_bonds": [list(pair) for pair in sorted(second.keys() - first.keys())],
            "bond_order_changes": [
                {"atom_indices": list(pair), "reactant_order": first[pair], "product_order": second[pair]}
                for pair in sorted(first.keys() & second.keys())
                if abs(first[pair] - second[pair]) > 1e-8
            ],
        }

    def to_dict(self) -> dict:
        return {"atoms": [atom.to_dict() for atom in self.atoms],
                "reactants": self.reactants.to_dict(), "products": self.products.to_dict(),
                "charge": self.charge, "multiplicity": self.multiplicity,
                "environment": self.environment, "allowed_motions": self.allowed_motions,
                "parent_evidence": json.loads(self.parent_evidence_json)}


def derive_geometry_context(*, atoms: Sequence[Mapping[str, Any]],
                            reactants: Mapping[str, Any], products: Mapping[str, Any],
                            charge: int = 0, multiplicity: int = 1,
                            environment: Mapping[str, Any] | None = None,
                            allowed_motions: str = "free_molecular_fragments",
                            parent_evidence: Sequence[Any] = ()) -> ReactionGeometryContext:
    """Normalize already mapped graph/coordinate inputs without guessing atom maps.

    Fragment IDs may be any integer labels. Coordinates for different fragments
    may initially overlap: this is their input pose, not a qualified encounter.
    Within each fragment the graph must be connected. Rigid preparation cannot
    repair an invalid or conformationally obstructed individual molecule.
    """
    identities = []
    for value in atoms:
        mapping = value.get("map_number")
        isotope = value.get("isotope", 0)
        if isinstance(mapping, bool) or not isinstance(mapping, int) or mapping <= 0:
            raise GeometryPreparationError("Every nucleus needs a positive explicit map number")
        if isinstance(isotope, bool) or not isinstance(isotope, int) or isotope < 0:
            raise GeometryPreparationError("Isotope labels must be nonnegative integers")
        if value.get("element") not in _ELEMENTS:
            raise GeometryPreparationError("Unknown atomic element")
        stereo = value.get("stereochemistry")
        if stereo is not None and not isinstance(stereo, str):
            raise GeometryPreparationError("Stereochemistry must be an explicit label or null")
        identities.append(AtomIdentity(mapping, value["element"], isotope, stereo))
    if not identities or len({atom.map_number for atom in identities}) != len(identities):
        raise GeometryPreparationError("A nonempty, unique complete nuclear map is required")
    if isinstance(charge, bool) or not isinstance(charge, int):
        raise GeometryPreparationError("Charge must be an integer")
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) or multiplicity < 1:
        raise GeometryPreparationError("Multiplicity must be a positive integer")
    environment = dict(environment or {})
    if allowed_motions != "free_molecular_fragments":
        raise UnsupportedGeometryContext("This proposal policy supports only free molecular fragment motions")
    # Carry solution/gas declarations, but never rotate a fixed boundary away.
    if (environment.get("fixed_atom_indices") or environment.get("periodic_cell")
            or environment.get("participating_solvent")
            or str(environment.get("kind", environment.get("phase", ""))).lower()
            in {"surface", "solid", "periodic", "interface", "electrochemical", "qm/mm", "qmmm"}):
        raise UnsupportedGeometryContext("Declared environment requires a specialized allowed-motion model")
    count = len(identities)

    def endpoint(value: Mapping[str, Any]) -> EndpointGeometry:
        labels = value.get("fragment_ids")
        if labels is None or len(labels) != count:
            raise GeometryPreparationError("Each endpoint needs one fragment ID per nucleus")
        if any(isinstance(v, bool) or not isinstance(v, (int, np.integer)) for v in labels):
            raise GeometryPreparationError("Fragment IDs must be integers")
        labels = tuple(int(v) for v in labels)
        bonds = _bonds(value.get("bonds", ()), count)
        adjacency = {index: set() for index in range(count)}
        for a, b, _ in bonds:
            if labels[a] != labels[b]:
                raise GeometryPreparationError("A declared covalent bond crosses fragment IDs")
            adjacency[a].add(b)
            adjacency[b].add(a)
        for label in set(labels):
            members = {i for i, v in enumerate(labels) if v == label}
            reached, pending = set(), [min(members)]
            while pending:
                index = pending.pop()
                if index not in reached:
                    reached.add(index)
                    pending.extend(adjacency[index] - reached)
            if reached != members:
                raise GeometryPreparationError("Each declared fragment must be graph-connected")
        bindings = value.get("source_bindings", ())
        for binding in bindings:
            if not isinstance(binding, Mapping) or not binding.get("path") or len(str(binding.get("sha256", ""))) != 64:
                raise GeometryPreparationError("Source bindings require a path and SHA256")
        return EndpointGeometry(_coordinates(value["coordinates"], count), labels, bonds,
                                _json(list(bindings)), _json(value.get("conformer_ids", [])),
                                _json(value.get("fragment_states", [])))

    return ReactionGeometryContext(tuple(identities), endpoint(reactants), endpoint(products),
                                   charge, multiplicity, _json(environment), allowed_motions,
                                   _json(list(parent_evidence)))


@dataclass(frozen=True)
class SeedGeometryPolicy:
    max_candidates: int = 24
    encounter_distance_angstrom: float = 2.8
    separation_increments_angstrom: tuple[float, ...] = (0.0, 0.5, 1.0)
    bend_angles_degrees: tuple[float, ...] = (0.0, 15.0, -15.0)
    azimuth_angles_degrees: tuple[float, ...] = (0.0, 120.0, 240.0)
    rigid_iterations: int = 48
    duplicate_rmsd_angstrom: float = 0.035

    def __post_init__(self) -> None:
        if isinstance(self.max_candidates, bool) or not isinstance(self.max_candidates, int) or not 1 <= self.max_candidates <= 96:
            raise GeometryPreparationError("Proposal count must be between 1 and 96")
        if not math.isfinite(self.encounter_distance_angstrom) or self.encounter_distance_angstrom <= 1:
            raise GeometryPreparationError("Encounter distance must exceed 1 Angstrom")
        for name in ("separation_increments_angstrom", "bend_angles_degrees", "azimuth_angles_degrees"):
            values = getattr(self, name)
            if not values or not all(math.isfinite(v) for v in values):
                raise GeometryPreparationError("Proposal parameter grids must be finite and nonempty")
            object.__setattr__(self, name, tuple(float(value) for value in values))
        if any(v < 0 for v in self.separation_increments_angstrom):
            raise GeometryPreparationError("Separation increments cannot be negative")
        if (isinstance(self.rigid_iterations, bool) or not isinstance(self.rigid_iterations, int)
                or not 0 <= self.rigid_iterations <= 200
                or not math.isfinite(self.duplicate_rmsd_angstrom) or self.duplicate_rmsd_angstrom <= 0):
            raise GeometryPreparationError("Invalid iteration or deduplication bound")

    def to_dict(self) -> dict:
        return {"policy_id": POLICY_ID, "max_candidates": self.max_candidates,
                "encounter_distance_angstrom": self.encounter_distance_angstrom,
                "separation_increments_angstrom": list(self.separation_increments_angstrom),
                "bend_angles_degrees": list(self.bend_angles_degrees),
                "azimuth_angles_degrees": list(self.azimuth_angles_degrees),
                "rigid_iterations": self.rigid_iterations,
                "duplicate_rmsd_angstrom": self.duplicate_rmsd_angstrom}


@dataclass(frozen=True)
class SeedCandidate:
    candidate_id: str
    context_id: str
    generator: str
    reactant_coordinates: tuple[tuple[float, float, float], ...]
    product_coordinates: tuple[tuple[float, float, float], ...]
    parameters_json: str
    metrics_json: str
    rejection_reasons: tuple[str, ...]
    diversity_aliases: tuple[str, ...] = ()

    @property
    def metrics(self) -> dict:
        return json.loads(self.metrics_json)

    @property
    def parameters(self) -> dict:
        return json.loads(self.parameters_json)

    @property
    def valid(self) -> bool:
        return not self.rejection_reasons

    def to_dict(self) -> dict:
        return {"candidate_id": self.candidate_id, "context_id": self.context_id,
                "generator": self.generator,
                "reactant_coordinates": [list(row) for row in self.reactant_coordinates],
                "product_coordinates": [list(row) for row in self.product_coordinates],
                "parameters": self.parameters, "metrics": self.metrics,
                "valid": self.valid, "rejection_reasons": list(self.rejection_reasons),
                "diversity_aliases": list(self.diversity_aliases)}


def _groups(side: EndpointGeometry) -> list[list[int]]:
    return [[i for i, value in enumerate(side.fragment_ids) if value == label]
            for label in sorted(set(side.fragment_ids), key=lambda label: min(
                i for i, value in enumerate(side.fragment_ids) if value == label))]


def _unit(vector: np.ndarray) -> np.ndarray | None:
    length = float(np.linalg.norm(vector))
    return vector / length if length > 1e-10 else None


def _rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = _unit(axis)
    if axis is None or abs(angle) < 1e-14:
        return np.eye(3)
    x, y, z = axis
    skew = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    return np.eye(3) + math.sin(angle) * skew + (1 - math.cos(angle)) * (skew @ skew)


def _perpendicular(axis: np.ndarray, vectors: Sequence[np.ndarray]) -> tuple[np.ndarray, bool]:
    for vector in vectors:
        projected = vector - np.dot(vector, axis) * axis
        value = _unit(projected)
        if value is not None:
            return value, False
    # A perfectly atomic/linear context has no physically distinguished azimuth.
    # This laboratory-frame choice fixes a gauge only; it is explicitly recorded,
    # and no unique equivariant frame or additional physical diversity is claimed.
    gauge = np.eye(3)[int(np.argmin(np.abs(axis)))]
    return _unit(gauge - np.dot(gauge, axis) * axis), True


def _fit(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Proper rigid fit only; never reflect a stereocenter."""
    left_center, right_center = source.mean(axis=0), target.mean(axis=0)
    if len(source) == 1:
        return source - left_center + right_center
    u, _, vt = np.linalg.svd((source - left_center).T @ (target - right_center))
    correction = np.eye(3)
    correction[-1, -1] = 1. if np.linalg.det(u @ vt) >= 0 else -1.
    return (source - left_center) @ (u @ correction @ vt) + right_center


def _rotate_to(vector: np.ndarray, target: np.ndarray,
               transverse: np.ndarray) -> np.ndarray:
    source, target = _unit(vector), _unit(target)
    if source is None or target is None:
        return np.eye(3)
    cross = np.cross(source, target)
    length = float(np.linalg.norm(cross))
    if length > 1e-10:
        return _rotation(cross / length, math.atan2(length, np.dot(source, target)))
    return _rotation(transverse, math.pi) if np.dot(source, target) < 0 else np.eye(3)


def _transfer_motifs(context: ReactionGeometryContext, side_name: str) -> list[tuple[int, int, int]]:
    changes = context.bond_changes
    broken = changes["broken_bonds"] if side_name == "reactants" else changes["formed_bonds"]
    formed = changes["formed_bonds"] if side_name == "reactants" else changes["broken_bonds"]
    result = []
    for hydrogen, atom in enumerate(context.atoms):
        if atom.element != "H":
            continue
        donors = [b if a == hydrogen else a for a, b in broken if hydrogen in (a, b)]
        acceptors = [b if a == hydrogen else a for a, b in formed if hydrogen in (a, b)]
        if len(donors) == len(acceptors) == 1:
            result.append((donors[0], hydrogen, acceptors[0]))
    return result


def _contact_pairs(context: ReactionGeometryContext, side_name: str) -> list[tuple[int, int]]:
    side = getattr(context, side_name)
    changes = context.bond_changes
    values = changes["formed_bonds"] if side_name == "reactants" else changes["broken_bonds"]
    # Include every changed contact, but only interfragment ones are adjustable
    # with rigid motions. Internal conformational exploration is not implemented.
    return [tuple(pair) for pair in values if side.fragment_ids[pair[0]] != side.fragment_ids[pair[1]]]


def _outward(side: EndpointGeometry, coordinates: np.ndarray, atom: int,
             counterpart: int, reference: np.ndarray | None) -> np.ndarray | None:
    neighbours = [b if a == atom else a for a, b, _ in side.bonds if atom in (a, b)]
    vectors = [_unit(coordinates[atom] - coordinates[i]) for i in neighbours]
    value = _unit(sum((v for v in vectors if v is not None), np.zeros(3)))
    if value is not None:
        return value
    if reference is not None:
        value = _unit(reference[counterpart] - reference[atom])
        if value is not None:
            return value
    return None


def _place_side(context: ReactionGeometryContext, side_name: str, *, distance: float,
                bend: float, azimuth: float, iterations: int,
                reference: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    side = getattr(context, side_name)
    coordinates = np.array(side.coordinates, dtype=float)
    atom_order = sorted(range(len(context.atoms)), key=lambda atom: context.atoms[atom].map_number)
    groups = sorted(_groups(side), key=lambda group: min(context.atoms[atom].map_number for atom in group))
    groups = [sorted(group, key=lambda atom: context.atoms[atom].map_number) for group in groups]
    if reference is not None:
        for group in groups:
            coordinates[group] = _fit(coordinates[group], reference[group])
    if len(groups) == 1:
        return coordinates, {"adjusted_interfragment_contacts": 0, "azimuth_gauge_only": False,
                             "internal_conformer_search": "not_performed"}
    anchor = max(range(len(groups)), key=lambda i: (sum(context.atoms[a].element != "H" for a in groups[i]),
                                                   len(groups[i]), -min(context.atoms[a].map_number for a in groups[i])))
    contacts = sorted(_contact_pairs(context, side_name), key=lambda pair: sorted(context.atoms[a].map_number for a in pair))
    motifs = _transfer_motifs(context, side_name)
    placed = set(groups[anchor])
    pending = [i for i in range(len(groups)) if i != anchor]
    gauge_only = False
    while pending:
        next_group = max(pending, key=lambda i: (sum((a in groups[i] and b in placed) or
                                                     (b in groups[i] and a in placed) for a, b in contacts), -i))
        pending.remove(next_group)
        moving = groups[next_group]
        relevant = [(a, b) if a in placed else (b, a) for a, b in contacts
                    if (a in placed and b in moving) or (b in placed and a in moving)]
        if not relevant:
            # No connecting edit: this fragment is an uninvolved spectator. Keep
            # a finite, nonoverlapping proposal without calling it an asymptote.
            center = coordinates[[atom for atom in atom_order if atom in placed]].mean(axis=0)
            vectors = [coordinates[a] - center for a in atom_order if a in placed]
            axis = next((v for vector in vectors if (v := _unit(vector)) is not None), np.array([1., 0., 0.]))
            radius = max(np.linalg.norm(coordinates[a] - center) for a in placed)
            moving_center = coordinates[moving].mean(axis=0)
            moving_radius = max(np.linalg.norm(coordinates[a] - moving_center) for a in moving)
            coordinates[moving] += center + axis * (radius + moving_radius + distance) - moving_center
            placed.update(moving)
            continue
        targets, sites, directions = [], [], []
        common_direction = None
        fixed_sites = sorted({fixed for fixed, _ in relevant}, key=lambda atom: context.atoms[atom].map_number)
        if len(fixed_sites) > 1:
            span = _unit(coordinates[fixed_sites[1]] - coordinates[fixed_sites[0]])
            if span is not None:
                internal_vectors = [coordinates[a] - coordinates[fixed_sites[0]] for a in atom_order if a in placed]
                internal_vectors.extend(coordinates[a] - coordinates[moving[0]] for a in moving)
                common_direction, is_gauge = _perpendicular(span, internal_vectors)
                gauge_only |= is_gauge
                # Three noncollinear contacting sites define an approach normal.
                # Do not make independent outward rays that pull the moving
                # sites in mutually incompatible directions.
                normal = next((normal for atom in fixed_sites[2:]
                               if (normal := _unit(np.cross(span, coordinates[atom] - coordinates[fixed_sites[0]]))) is not None), None)
                if normal is not None:
                    common_direction = normal
        for fixed, mobile in relevant:
            direction = common_direction
            for donor, hydrogen, acceptor in motifs:
                if direction is None and fixed == hydrogen and mobile == acceptor:
                    direction = _unit(coordinates[hydrogen] - coordinates[donor])
                    break
            if direction is None:
                direction = _outward(side, coordinates, fixed, mobile, reference)
            if direction is None:
                # Derive a frame from internal vectors, not the input fragment
                # centers (which are frequently all at zero).
                direction = next((value for group in groups for a in group for b in group
                                  if a < b and (value := _unit(coordinates[b] - coordinates[a])) is not None),
                                 np.array([1., 0., 0.]))
            vectors = [coordinates[a] - coordinates[fixed] for a in atom_order if a in placed]
            vectors.extend(coordinates[a] - coordinates[mobile] for a in moving)
            transverse, is_gauge = _perpendicular(direction, vectors)
            gauge_only |= is_gauge
            transverse = _rotation(direction, math.radians(azimuth)) @ transverse
            direction = math.cos(math.radians(bend)) * direction + math.sin(math.radians(bend)) * transverse
            targets.append(coordinates[fixed] + distance * direction)
            sites.append(mobile)
            directions.append(direction)
        if len(sites) > 1:
            # Repeated moving sites can legitimately have multiple contacts.
            # Their targets are averaged, then the full contact objective below
            # resolves distances jointly rather than choosing relevant[0].
            unique_sites = sorted(set(sites), key=lambda atom: context.atoms[atom].map_number)
            averaged = np.array([np.mean([target for site, target in zip(sites, targets) if site == atom], axis=0)
                                 for atom in unique_sites])
            source = coordinates[unique_sites]
            center = source.mean(axis=0)
            fitted = _fit(source, averaged)
            # Obtain the same proper transform for every atom in this fragment.
            if len(unique_sites) >= 2:
                u, _, vt = np.linalg.svd((source - center).T @ (averaged - averaged.mean(axis=0)))
                correction = np.eye(3)
                correction[-1, -1] = 1. if np.linalg.det(u @ vt) >= 0 else -1.
                coordinates[moving] = (coordinates[moving] - center) @ (u @ correction @ vt) + averaged.mean(axis=0)
            else:
                coordinates[moving] += fitted[0] - source[0]
        else:
            site = sites[0]
            outward = _outward(side, coordinates, site, relevant[0][0], reference)
            if outward is not None:
                transverse, is_gauge = _perpendicular(outward, [coordinates[a] - coordinates[site] for a in moving])
                gauge_only |= is_gauge
                rotation = _rotate_to(outward, -directions[0], transverse)
                coordinates[moving] = (coordinates[moving] - coordinates[site]) @ rotation.T + coordinates[site]
            coordinates[moving] += targets[0] - coordinates[site]
        placed.update(moving)
    coordinates = _relax_rigid(context, side, coordinates, groups, anchor, contacts, distance, iterations)
    return coordinates, {"adjusted_interfragment_contacts": len(contacts),
                         "azimuth_gauge_only": gauge_only,
                         "internal_conformer_search": "not_performed"}


def _rigid_objective(context: ReactionGeometryContext, side: EndpointGeometry,
                     coordinates: np.ndarray, contacts: list[tuple[int, int]],
                     distance: float) -> tuple[float, np.ndarray]:
    gradient = np.zeros_like(coordinates)
    value = 0.0
    contact_set = {tuple(sorted(pair)) for pair in contacts}
    for left in range(len(coordinates)):
        for right in range(left + 1, len(coordinates)):
            if side.fragment_ids[left] == side.fragment_ids[right]:
                continue
            delta = coordinates[left] - coordinates[right]
            separation = float(np.linalg.norm(delta))
            if separation < 1e-10:
                # A collapsed pair cannot have a useful radial derivative, but
                # must never disappear from the objective or look inexpensive.
                value += distance ** 2 if (left, right) in contact_set else 100.
                continue
            if (left, right) in contact_set:
                error, weight = separation - distance, 1.0
            else:
                threshold = .95 * (_RADII.get(context.atoms[left].element, .77) + _RADII.get(context.atoms[right].element, .77))
                error, weight = min(0., separation - threshold), 4.0
            value += weight * error ** 2
            force = 2 * weight * error * delta / separation
            gradient[left] += force
            gradient[right] -= force
    return value, gradient


def _relax_rigid(context: ReactionGeometryContext, side: EndpointGeometry,
                 coordinates: np.ndarray, groups: list[list[int]], anchor: int,
                 contacts: list[tuple[int, int]], distance: float, iterations: int) -> np.ndarray:
    """Bounded rigid-body descent; all contacts contribute and no bonds distort."""
    for _ in range(iterations):
        old_value, gradient = _rigid_objective(context, side, coordinates, contacts, distance)
        if old_value < 1e-12:
            break
        updates = []
        for index, group in enumerate(groups):
            if index == anchor:
                continue
            centered = coordinates[group] - coordinates[group].mean(axis=0)
            translation = -gradient[group].mean(axis=0)
            torque = -np.cross(centered, gradient[group]).sum(axis=0)
            torque /= max(1., float(np.sum(centered * centered)))
            updates.append((group, centered, translation, torque))
        accepted = False
        for step in (.25, .125, .0625, .03125, .015625):
            trial = coordinates.copy()
            for group, centered, translation, torque in updates:
                turn = min(.15, float(np.linalg.norm(torque)) * step)
                trial[group] = centered @ _rotation(torque, turn).T + coordinates[group].mean(axis=0) + step * translation
            new_value, _ = _rigid_objective(context, side, trial, contacts, distance)
            if new_value < old_value - 1e-12:
                coordinates, accepted = trial, True
                break
        if not accepted:
            break
    return coordinates


def assess_seed_pair(candidate: SeedCandidate, context: ReactionGeometryContext,
                     policy: SeedGeometryPolicy | None = None) -> dict:
    """Geometry-only integrity and ranking; never qualifies a physical endpoint."""
    policy = policy or SeedGeometryPolicy()
    if candidate.context_id != context.context_id:
        raise GeometryPreparationError("Candidate belongs to a different mapped/state/environment context")
    arrays = [np.array(_coordinates(candidate.reactant_coordinates, len(context.atoms))),
              np.array(_coordinates(candidate.product_coordinates, len(context.atoms)))]
    severe, crowding, contact_errors = 0, 0.0, []
    reasons = []
    for side_name, coordinates in zip(("reactants", "products"), arrays):
        side = getattr(context, side_name)
        original = np.array(side.coordinates)
        contacts = set(_contact_pairs(context, side_name))
        bonded = {(a, b) for a, b, _ in side.bonds}
        for group in _groups(side):
            # Proper rigid superposition detects reflection as well as unwanted
            # internal deformation; no ideal angle or strained-ring cutoff.
            if np.max(np.linalg.norm(_fit(original[group], coordinates[group]) - coordinates[group], axis=1)) > 1e-6:
                reasons.append(f"{side_name}:fragment_internal_geometry_or_stereochemistry_changed")
                break
        for a in range(len(coordinates)):
            for b in range(a + 1, len(coordinates)):
                separation = float(np.linalg.norm(coordinates[a] - coordinates[b]))
                radii = _RADII.get(context.atoms[a].element, .77) + _RADII.get(context.atoms[b].element, .77)
                if side.fragment_ids[a] == side.fragment_ids[b]:
                    # Only actual collapse, not unusual legitimate geometry.
                    if separation < (.18 if (a, b) in bonded else .25) * radii:
                        severe += 1
                elif (a, b) in contacts:
                    if separation < .30 * radii:
                        severe += 1
                    contact_errors.append(abs(separation - candidate.parameters.get(
                        "encounter_distance_angstrom", policy.encounter_distance_angstrom)))
                else:
                    severe += separation < .45 * radii
                    crowding += max(0., .95 * radii - separation) ** 2
    if severe:
        reasons.append("collapsed_or_impossible_contact")
    # This measures mapped endpoint motion, not a calculated reaction path.
    aligned_product = _fit(arrays[1], arrays[0])
    displacement = np.linalg.norm(arrays[0] - aligned_product, axis=1)
    motifs = _transfer_motifs(context, "reactants")
    misalignment = []
    for donor, hydrogen, attacker in motifs:
        donor_axis = _unit(arrays[0][hydrogen] - arrays[0][donor])
        approach = _unit(arrays[0][attacker] - arrays[0][hydrogen])
        if donor_axis is not None and approach is not None:
            misalignment.append(1. - float(np.clip(np.dot(donor_axis, approach), -1, 1)))
    return {"valid": not reasons, "rejection_reasons": sorted(set(reasons)),
            "severe_contact_count": int(severe), "spectator_crowding_penalty": crowding,
            "maximum_reactive_contact_error_angstrom": max(contact_errors, default=0.),
            "hydrogen_transfer_misalignment": max(misalignment, default=0.),
            "paired_mapped_displacement_rmsd_angstrom": float(np.sqrt(np.mean(displacement ** 2))),
            "paired_maximum_displacement_angstrom": float(max(displacement, default=0.)),
            "unknown_covalent_radius_elements": sorted({a.element for a in context.atoms if a.element not in _RADII}),
            "interpretation": "Proposal screening only; no optimized minimum, TS, mechanism, or rate claim."}


def propose_seed_pairs(context: ReactionGeometryContext,
                       policy: SeedGeometryPolicy | None = None) -> tuple[SeedCandidate, ...]:
    policy = policy or SeedGeometryPolicy()
    generated = []
    motifs = _transfer_motifs(context, "reactants")
    is_multifragment = max(len(_groups(context.reactants)), len(_groups(context.products))) > 1
    grid = [(increment, bend, azimuth)
            for increment in policy.separation_increments_angstrom
            for bend in policy.bend_angles_degrees
            for azimuth in policy.azimuth_angles_degrees]
    if not is_multifragment:
        grid = [(0., 0., 0.)]
    for increment, bend, azimuth in grid[:policy.max_candidates]:
        distance = policy.encounter_distance_angstrom + increment
        reactant, reactant_record = _place_side(context, "reactants", distance=distance, bend=bend,
                                                azimuth=azimuth, iterations=policy.rigid_iterations)
        product, product_record = _place_side(context, "products", distance=distance, bend=bend,
                                              azimuth=azimuth, iterations=policy.rigid_iterations,
                                              reference=reactant)
        # Preserve relative paired pose while removing avoidable whole-system
        # product drift with a proper mapped alignment.
        product = _fit(product, reactant)
        parameters = {"policy": policy.to_dict(), "encounter_distance_angstrom": distance,
                      "bend_degrees": bend, "azimuth_degrees": azimuth,
                      "reactant_preconditioning": reactant_record,
                      "product_preconditioning": product_record,
                      "internal_conformer_coverage": "not_established"}
        generator = "hydrogen_transfer_and_paired_rigid_contacts" if motifs else "paired_rigid_all_contacts"
        raw = SeedCandidate("", context.context_id, generator,
                            _coordinates(reactant, len(context.atoms)), _coordinates(product, len(context.atoms)),
                            _json(parameters), "{}", ())
        metrics = assess_seed_pair(raw, context, policy)
        generated.append(replace(raw, candidate_id=_candidate_hash(raw), metrics_json=_json(metrics),
                                 rejection_reasons=tuple(metrics["rejection_reasons"])))
    return tuple(generated)


def _pair_rmsd(first: SeedCandidate, second: SeedCandidate) -> float:
    left = np.vstack((first.reactant_coordinates, first.product_coordinates))
    right = np.vstack((second.reactant_coordinates, second.product_coordinates))
    return float(np.sqrt(np.mean(np.sum((_fit(right, left) - left) ** 2, axis=1))))


def _candidate_hash(candidate: SeedCandidate) -> str:
    return _sha({"context_id": candidate.context_id, "generator": candidate.generator,
                 "reactants": candidate.reactant_coordinates, "products": candidate.product_coordinates,
                 "parameters": candidate.parameters})


def rank_diverse_candidates(candidates: Sequence[SeedCandidate], evidence: Any = None,
                            policy: SeedGeometryPolicy | None = None) -> tuple[SeedCandidate, ...]:
    """Lexicographic proposal ranking with proper, map-aware geometric dedup.

    ``evidence`` is reserved for explicit future low-level assessments. No energy,
    learned probability, rate, degeneracy, or experimental target is inferred.
    Rejected candidates remain in the sealed full portfolio, not this frontier.
    """
    policy = policy or SeedGeometryPolicy()
    if len({candidate.context_id for candidate in candidates}) > 1:
        raise GeometryPreparationError("Cannot deduplicate different maps, isotopes, conformers, states, or environments")
    # Numerical SVD/rotation noise must not outrank physically meaningful access.
    ranked = sorted((candidate for candidate in candidates if candidate.valid), key=lambda candidate: (
        candidate.metrics["severe_contact_count"], round(candidate.metrics["spectator_crowding_penalty"], 8),
        round(candidate.metrics["maximum_reactive_contact_error_angstrom"], 8),
        round(candidate.metrics["hydrogen_transfer_misalignment"], 8),
        round(candidate.metrics["paired_mapped_displacement_rmsd_angstrom"], 8), candidate.candidate_id))
    selected: list[SeedCandidate] = []
    for candidate in ranked:
        duplicate = next((i for i, previous in enumerate(selected)
                          if _pair_rmsd(previous, candidate) < policy.duplicate_rmsd_angstrom), None)
        if duplicate is None:
            selected.append(candidate)
        else:
            previous = selected[duplicate]
            selected[duplicate] = replace(previous, diversity_aliases=tuple(sorted(set(
                previous.diversity_aliases + (candidate.candidate_id,) + candidate.diversity_aliases))))
    return tuple(selected)


def _xyz(context: ReactionGeometryContext, coordinates: Sequence[Sequence[float]], candidate_id: str, side: str) -> str:
    return "\n".join([str(len(context.atoms)), f"STORCA geometry proposal only; candidate={candidate_id}; side={side}"] + [
        f"{atom.element:<2} {row[0]: .15f} {row[1]: .15f} {row[2]: .15f}"
        for atom, row in zip(context.atoms, coordinates)]) + "\n"


def _immutable_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = content.encode()
    try:
        with path.open("xb") as stream:
            stream.write(data)
    except FileExistsError:
        if path.read_bytes() != data:
            raise GeometryPreparationError(f"Immutable geometry artifact differs: {path}")


def _verify_sources(context: ReactionGeometryContext) -> None:
    for side in (context.reactants, context.products):
        for binding in json.loads(side.source_bindings_json):
            path = Path(binding["path"])
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != binding["sha256"]:
                raise GeometryPreparationError(f"Fragment source binding changed: {path}")


def seal_seed_portfolio(context: ReactionGeometryContext, candidates: Sequence[SeedCandidate],
                        output_dir: str | Path, policy: SeedGeometryPolicy | None = None) -> dict:
    """Write or exactly reopen content-addressed proposals; never overwrite seeds."""
    policy = policy or SeedGeometryPolicy()
    output_dir = Path(output_dir).resolve()
    _verify_sources(context)
    frontier = rank_diverse_candidates(candidates, policy=policy)
    aliases = {candidate.candidate_id: list(candidate.diversity_aliases) for candidate in frontier}
    records = []
    for candidate in candidates:
        assessed = assess_seed_pair(candidate, context, policy)
        if (candidate.candidate_id != _candidate_hash(candidate)
                or candidate.parameters.get("policy") != policy.to_dict()
                or assessed != candidate.metrics or tuple(assessed["rejection_reasons"]) != candidate.rejection_reasons):
            raise GeometryPreparationError("Candidate assessment was altered or uses a different policy")
        record = candidate.to_dict()
        record["diversity_aliases"] = aliases.get(candidate.candidate_id, [])
        folder = output_dir / "candidates" / candidate.candidate_id
        for side, coordinates in (("reactant", candidate.reactant_coordinates), ("product", candidate.product_coordinates)):
            path = folder / f"{side}-seed.xyz"
            content = _xyz(context, coordinates, candidate.candidate_id, side)
            _immutable_write(path, content)
            record[f"{side}_xyz"] = str(path)
            record[f"{side}_xyz_sha256"] = hashlib.sha256(content.encode()).hexdigest()
        _immutable_write(folder / "candidate.json", json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
        records.append(record)
    manifest = {"schema_version": 1, "kind": "immutable_paired_reaction_seed_portfolio",
                "policy": policy.to_dict(), "context_id": context.context_id, "context": context.to_dict(),
                "candidates": records, "frontier_candidate_ids": [candidate.candidate_id for candidate in frontier],
                "generated_count": len(candidates), "qualified_proposal_count": len(frontier),
                "interpretation": "Finite geometry proposal space, not reaction absence, route verification, or coverage proof."}
    _immutable_write(output_dir / "portfolio.json", json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return manifest


def reopen_seed_portfolio(output_dir: str | Path, *, context: ReactionGeometryContext | None = None,
                          policy: SeedGeometryPolicy | None = None) -> dict:
    output_dir = Path(output_dir).resolve()
    original = (output_dir / "portfolio.json").read_bytes()
    manifest = json.loads(original)
    if manifest.get("schema_version") != 1 or manifest.get("kind") != "immutable_paired_reaction_seed_portfolio":
        raise GeometryPreparationError("Unsupported geometry portfolio schema")
    stored_context = derive_geometry_context(**manifest["context"])
    if context is not None and context.context_id != stored_context.context_id:
        raise GeometryPreparationError("Retained portfolio has a different scientific context")
    policy_data = dict(manifest["policy"])
    if policy_data.pop("policy_id", None) != POLICY_ID:
        raise GeometryPreparationError("Unsupported geometry proposal policy")
    for name in ("separation_increments_angstrom", "bend_angles_degrees", "azimuth_angles_degrees"):
        policy_data[name] = tuple(policy_data[name])
    stored_policy = SeedGeometryPolicy(**policy_data)
    if policy is not None and policy.to_dict() != stored_policy.to_dict():
        raise GeometryPreparationError("Retained geometry policy changed")
    candidates = []
    for record in manifest["candidates"]:
        candidates.append(SeedCandidate(
            record["candidate_id"], record["context_id"], record["generator"],
            _coordinates(record["reactant_coordinates"], len(stored_context.atoms)),
            _coordinates(record["product_coordinates"], len(stored_context.atoms)),
            _json(record["parameters"]), _json(record["metrics"]), tuple(record["rejection_reasons"])))
    # Byte-exact comparison includes canonical candidate manifests, XYZ, source
    # bindings and assessment recomputation. Partial writes may be completed by
    # an explicit sealing operation, but reopen itself must not create files.
    _verify_sources(stored_context)
    for candidate, record in zip(candidates, manifest["candidates"]):
        if (candidate.candidate_id != _candidate_hash(candidate)
                or candidate.parameters.get("policy") != stored_policy.to_dict()
                or assess_seed_pair(candidate, stored_context, stored_policy) != candidate.metrics
                or candidate.valid != record.get("valid")):
            raise GeometryPreparationError("Retained candidate assessment changed")
        for side, coordinates in (("reactant", candidate.reactant_coordinates), ("product", candidate.product_coordinates)):
            path = output_dir / "candidates" / candidate.candidate_id / f"{side}-seed.xyz"
            expected = _xyz(stored_context, coordinates, candidate.candidate_id, side).encode()
            if path.read_bytes() != expected or str(path) != record[f"{side}_xyz"] or hashlib.sha256(expected).hexdigest() != record[f"{side}_xyz_sha256"]:
                raise GeometryPreparationError("Retained candidate XYZ or binding changed")
        path = output_dir / "candidates" / candidate.candidate_id / "candidate.json"
        if path.read_bytes() != (json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n").encode():
            raise GeometryPreparationError("Retained candidate manifest changed")
    # A forged frontier is not accepted merely because the XYZ files match.
    frontier = rank_diverse_candidates(candidates, policy=stored_policy)
    if manifest["frontier_candidate_ids"] != [candidate.candidate_id for candidate in frontier]:
        raise GeometryPreparationError("Retained proposal frontier changed")
    aliases = {candidate.candidate_id: list(candidate.diversity_aliases) for candidate in frontier}
    if (manifest["generated_count"] != len(candidates) or manifest["qualified_proposal_count"] != len(frontier)
            or any(record["diversity_aliases"] != aliases.get(record["candidate_id"], [])
                   for record in manifest["candidates"])):
        raise GeometryPreparationError("Retained proposal counts or aliases changed")
    if manifest["context_id"] != stored_context.context_id or original != (json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n").encode():
        raise GeometryPreparationError("Retained portfolio identity or serialization changed")
    return manifest
