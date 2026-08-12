"""Bounded, deterministic self-dimer sampling for condensed-phase IR spectra."""

from __future__ import annotations

import json
import math
from pathlib import Path
import random

import numpy as np


def classify_local_stretch_bonds(
    bonds: list[dict], interactions: list[dict], monomer_atom_count: int,
) -> list[dict]:
    """Annotate each local X-H oscillator by its molecule's H-bond coordination."""
    annotated = []
    for bond in bonds:
        molecule_index = int(bond.get("molecule_index", int(bond["heavy_atom"]) // monomer_atom_count))
        start, stop = molecule_index * monomer_atom_count, (molecule_index + 1) * monomer_atom_count
        donates = any(int(item["donor_hydrogen"]) == int(bond["hydrogen_atom"]) for item in interactions)
        accepts = any(start <= int(item["acceptor_atom"]) < stop for item in interactions)
        coordination = (
            "donor_acceptor" if donates and accepts else
            "donor_only" if donates else
            "acceptor_only" if accepts else "free"
        )
        annotated.append({
            **bond,
            "donates_hydrogen_bond": donates,
            "accepts_hydrogen_bond": accepts,
            "coordination_class": coordination,
            "spectral_band_class": "hydrogen_bonded_oh" if donates else "non_donating_oh",
        })
    return annotated


def parse_xcontrol_interactions(path: Path) -> list[dict]:
    """Recover one-based xTB angle restraints as zero-based H-bond identities."""
    interactions = []
    for line in Path(path).read_text().splitlines():
        if "angle:" not in line:
            continue
        fields = line.split("angle:", 1)[1].strip().split(",")
        if len(fields) >= 3:
            interactions.append({
                "donor_atom": int(fields[0]) - 1,
                "donor_hydrogen": int(fields[1]) - 1,
                "acceptor_atom": int(fields[2]) - 1,
            })
    return interactions


def _rotation_between(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    source = source / np.linalg.norm(source)
    target = target / np.linalg.norm(target)
    cross = np.cross(source, target)
    cosine = float(np.dot(source, target))
    sine = float(np.linalg.norm(cross))
    if sine < 1e-12:
        if cosine > 0:
            return np.eye(3)
        axis = np.array([1.0, 0.0, 0.0])
        if abs(source[0]) > 0.8:
            axis = np.array([0.0, 1.0, 0.0])
        axis -= np.dot(axis, source) * source
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    skew = np.array([[0.0, -cross[2], cross[1]], [cross[2], 0.0, -cross[0]], [-cross[1], cross[0], 0.0]])
    return np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / (sine * sine))


def _axis_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    cross = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
    return np.eye(3) + math.sin(angle) * cross + (1.0 - math.cos(angle)) * (cross @ cross)


def _interaction_sites(molecule) -> tuple[list[tuple[int, int]], list[int]]:
    """Return donor-heavy/H pairs and conservative acceptor atom indices."""
    from rdkit.Chem import Lipinski

    donor_atoms = {match[0] for match in molecule.GetSubstructMatches(Lipinski.HDonorSmarts)}
    donor_atoms.update(
        atom.GetIdx() for atom in molecule.GetAtoms()
        if atom.GetSymbol() == "O" and atom.GetFormalCharge() == 0
        and any(neighbor.GetSymbol() == "H" for neighbor in atom.GetNeighbors())
        and not any(neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors())
    )
    acceptors = [match[0] for match in molecule.GetSubstructMatches(Lipinski.HAcceptorSmarts)]
    # Lipinski's conservative SMARTS intentionally excludes isolated water,
    # although its oxygen is the required acceptor in a water self-cluster.
    acceptors.extend(
        atom.GetIdx() for atom in molecule.GetAtoms()
        if atom.GetSymbol() == "O" and atom.GetFormalCharge() == 0
        and not any(neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors())
    )
    donors: list[tuple[int, int]] = []
    for atom_index in donor_atoms:
        atom = molecule.GetAtomWithIdx(atom_index)
        donors.extend(
            (atom_index, neighbor.GetIdx())
            for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == "H"
        )
    return donors, sorted(set(acceptors))


def _read_xyz_coordinates(path: Path, expected_symbols: list[str]) -> np.ndarray:
    lines = Path(path).read_text().splitlines()
    try:
        atom_count = int(lines[0].strip())
    except (IndexError, ValueError) as error:
        raise ValueError(f"Invalid XYZ atom-count line: {path}") from error
    rows = [line.split() for line in lines[2:2 + atom_count]]
    if len(rows) != atom_count or any(len(row) < 4 for row in rows):
        raise ValueError(f"Incomplete XYZ coordinates: {path}")
    symbols = [row[0] for row in rows]
    if symbols != expected_symbols:
        raise ValueError("Monomer XYZ atom ordering does not match the hydrogen-complete SMILES ordering")
    try:
        coordinates = np.asarray([[float(value) for value in row[1:4]] for row in rows], dtype=float)
    except ValueError as error:
        raise ValueError(f"Non-numeric XYZ coordinates: {path}") from error
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"Non-finite XYZ coordinates: {path}")
    return coordinates


def _diverse_parameter_order(
    parameters: list[dict], *, seed: int, maximum_items: int,
) -> list[dict]:
    """Return a deterministic maximin ordering over interaction coordinates."""
    shuffled = list(parameters)
    random.Random(seed).shuffle(shuffled)
    selected_indices: list[int] = []
    available = np.ones(len(shuffled), dtype=bool)
    site_pairs = sorted({item["site_pair"] for item in shuffled})
    for site_pair in site_pairs:
        site_items = [item for item in shuffled if item["site_pair"] == site_pair]
        coordinate_values = {
            key: sorted({item[key] for item in site_items})
            for key in (
                "h_bond_distance_angstrom", "donor_h_acceptor_angle_degrees",
                "axis_rotation_degrees", "framework_rotation_degrees",
            )
        }
        anchor_count = max(len(values) for values in coordinate_values.values())
        for offset in range(anchor_count):
            if len(selected_indices) >= maximum_items:
                break
            targets = {
                key: values[offset % len(values)]
                for key, values in coordinate_values.items()
            }
            anchor_index = next(
                index for index, item in enumerate(shuffled)
                if available[index] and item["site_pair"] == site_pair
                and all(item[key] == value for key, value in targets.items())
            )
            selected_indices.append(anchor_index)
            available[anchor_index] = False

    site_positions = {site_pair: position for position, site_pair in enumerate(site_pairs)}
    features = np.zeros((len(shuffled), 6 + len(site_pairs)), dtype=float)
    for index, item in enumerate(shuffled):
        azimuth = math.radians(item["axis_rotation_degrees"])
        framework = math.radians(item["framework_rotation_degrees"])
        features[index, :6] = (
            (item["h_bond_distance_angstrom"] - 1.6) / 0.9,
            (item["donor_h_acceptor_angle_degrees"] - 120.0) / 60.0,
            math.cos(azimuth), math.sin(azimuth),
            math.cos(framework), math.sin(framework),
        )
        features[index, 6 + site_positions[item["site_pair"]]] = 1.0
    if selected_indices:
        differences = features[:, None, :] - features[np.asarray(selected_indices), :][None, :, :]
        minimum_squared_distance = np.square(differences).sum(axis=2).min(axis=1)
    else:
        minimum_squared_distance = np.full(len(shuffled), np.inf)
    minimum_squared_distance[~available] = -1.0
    while available.any() and len(selected_indices) < maximum_items:
        candidate_index = int(np.argmax(minimum_squared_distance))
        selected_indices.append(candidate_index)
        available[candidate_index] = False
        squared_distance = np.square(features - features[candidate_index]).sum(axis=1)
        minimum_squared_distance = np.minimum(minimum_squared_distance, squared_distance)
        minimum_squared_distance[~available] = -1.0
    return [shuffled[index] for index in selected_indices]


def generate_stratified_dimer_poses(
    smiles: str,
    *,
    candidate_count: int,
    seed: int = 42,
    monomer_xyz: Path | None = None,
) -> list[dict]:
    """Build a deterministic, coverage-oriented set of directed H-bond poses."""
    if candidate_count < 1:
        raise ValueError("candidate_count must be at least one")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    molecule = Chem.AddHs(base)
    if Chem.GetFormalCharge(molecule) != 0:
        raise ValueError("Automatic environment poses are limited to neutral molecules")
    symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    if monomer_xyz is None:
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = seed
        if AllChem.EmbedMolecule(molecule, parameters) < 0:
            raise RuntimeError("Could not embed the molecule for environment sampling")
        try:
            AllChem.MMFFOptimizeMolecule(molecule, mmffVariant="MMFF94s")
        except Exception:
            AllChem.UFFOptimizeMolecule(molecule)
        coordinates = np.asarray([
            list(molecule.GetConformer().GetAtomPosition(index))
            for index in range(molecule.GetNumAtoms())
        ])
    else:
        coordinates = _read_xyz_coordinates(Path(monomer_xyz), symbols)

    donors, acceptors = _interaction_sites(molecule)
    if not donors or not acceptors:
        raise RuntimeError("No donor/acceptor pair supports automatic environment sampling")
    site_pairs = [
        (donor_atom, donor_hydrogen, acceptor)
        for donor_atom, donor_hydrogen in donors for acceptor in acceptors
    ]
    distances = (1.60, 1.75, 1.90, 2.10, 2.30, 2.50)
    angles = (120.0, 135.0, 150.0, 165.0, 180.0)
    axis_rotations = (0.0, 60.0, 120.0, 180.0, 240.0, 300.0)
    framework_rotations = (0.0, 90.0, 180.0, 270.0)
    parameter_pool = [
        {
            "site_pair": site_pair,
            "h_bond_distance_angstrom": distance,
            "donor_h_acceptor_angle_degrees": angle,
            "axis_rotation_degrees": axis_rotation,
            "framework_rotation_degrees": framework_rotation,
        }
        for site_pair in site_pairs
        for distance in distances
        for angle in angles
        for axis_rotation in axis_rotations
        for framework_rotation in framework_rotations
    ]
    ordered_parameters = _diverse_parameter_order(
        parameter_pool, seed=seed,
        maximum_items=min(len(parameter_pool), max(candidate_count * 3, candidate_count + 20)),
    )
    atom_count = molecule.GetNumAtoms()
    target_axis = np.array([1.0, 0.0, 0.0])
    transverse_y = np.array([0.0, 1.0, 0.0])
    transverse_z = np.array([0.0, 0.0, 1.0])
    poses: list[dict] = []
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in molecule.GetBonds()]

    for parameters in ordered_parameters:
        donor_atom, donor_hydrogen, acceptor = parameters["site_pair"]
        source_direction = coordinates[donor_atom] - coordinates[donor_hydrogen]
        if np.linalg.norm(source_direction) < 1e-8:
            continue
        theta = math.radians(parameters["donor_h_acceptor_angle_degrees"])
        phi = math.radians(parameters["axis_rotation_degrees"])
        transverse = math.cos(phi) * transverse_y + math.sin(phi) * transverse_z
        target_donor_direction = -math.cos(theta) * target_axis + math.sin(theta) * transverse
        alignment = _rotation_between(source_direction, target_donor_direction)
        donor_coordinates = (coordinates - coordinates[donor_hydrogen]) @ alignment.T
        framework_rotation = _axis_rotation(
            target_donor_direction, math.radians(parameters["framework_rotation_degrees"]),
        )
        donor_coordinates = donor_coordinates @ framework_rotation.T
        target_hydrogen = (
            coordinates[acceptor]
            + parameters["h_bond_distance_angstrom"] * target_axis
        )
        donor_coordinates += target_hydrogen
        pose = np.vstack([coordinates, donor_coordinates])
        inter_distances = np.linalg.norm(
            pose[:atom_count, None, :] - pose[None, atom_count:, :], axis=2,
        )
        if float(inter_distances.min()) < 0.65:
            continue
        acceptor_neighbors = [
            neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(acceptor).GetNeighbors()
            if neighbor.GetIdx() != donor_hydrogen
        ]
        orientation_atoms = None
        if parameters["donor_h_acceptor_angle_degrees"] < 175.0 and acceptor_neighbors:
            orientation_atoms = (
                acceptor_neighbors[0], acceptor,
                atom_count + donor_hydrogen, atom_count + donor_atom,
            )
        pose_index = len(poses) + 1
        poses.append({
            "candidate_id": f"candidate-{pose_index:04d}",
            "symbols": symbols + symbols,
            "coordinates": pose,
            "monomer_atom_count": atom_count,
            "bonds": bonds + [(left + atom_count, right + atom_count) for left, right in bonds],
            "donor_atom": atom_count + donor_atom,
            "donor_hydrogen": atom_count + donor_hydrogen,
            "acceptor_atom": acceptor,
            "orientation_atoms": orientation_atoms,
            "target_geometry": {
                key: value for key, value in parameters.items() if key != "site_pair"
            },
            "site_identity": {
                "donor_atom_in_monomer": donor_atom,
                "donor_hydrogen_in_monomer": donor_hydrogen,
                "acceptor_atom_in_monomer": acceptor,
            },
            "cluster_size": 2,
            "topology": "dimer",
            "molecule_atom_ranges": [[0, atom_count], [atom_count, 2 * atom_count]],
            "interactions": [{
                "donor_atom": atom_count + donor_atom,
                "donor_hydrogen": atom_count + donor_hydrogen,
                "acceptor_atom": acceptor,
                "orientation_atoms": list(orientation_atoms) if orientation_atoms else None,
                "target_geometry": {
                    key: value for key, value in parameters.items() if key != "site_pair"
                },
            }],
            "local_stretch_bonds": [
                {
                    "heavy_atom": molecule_index * atom_count + local_donor,
                    "hydrogen_atom": molecule_index * atom_count + local_hydrogen,
                    "bond_class": f"{symbols[local_donor]}-H:{local_donor}-{local_hydrogen}",
                    "molecule_index": molecule_index,
                }
                for molecule_index in range(2)
                for local_donor, local_hydrogen in donors
            ],
        })
        if len(poses) >= candidate_count:
            break
    if len(poses) < candidate_count:
        raise RuntimeError(
            f"Only {len(poses)} non-colliding poses were generated for {candidate_count} requested candidates"
        )
    return poses


def _place_donor_molecule(
    coordinates: np.ndarray, *, donor_atom: int, donor_hydrogen: int,
    acceptor_position: np.ndarray, outward_axis: np.ndarray,
    distance: float, angle_degrees: float, rotation_degrees: float,
) -> np.ndarray:
    axis = np.asarray(outward_axis, dtype=float)
    axis /= np.linalg.norm(axis)
    helper = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.8 else np.array([0.0, 1.0, 0.0])
    transverse = np.cross(axis, helper)
    transverse /= np.linalg.norm(transverse)
    theta = math.radians(angle_degrees)
    donor_direction = -math.cos(theta) * axis + math.sin(theta) * transverse
    source_direction = coordinates[donor_atom] - coordinates[donor_hydrogen]
    alignment = _rotation_between(source_direction, donor_direction)
    placed = (coordinates - coordinates[donor_hydrogen]) @ alignment.T
    placed = placed @ _axis_rotation(donor_direction, math.radians(rotation_degrees)).T
    placed += acceptor_position + distance * axis
    return placed


def generate_stratified_trimer_poses(
    smiles: str, *, candidate_count: int, seed: int = 43,
    monomer_xyz: Path | None = None, cyclic_fraction: float = 0.25,
) -> list[dict]:
    """Extend diverse dimers into restrained donor-acceptor chains and plausible rings."""
    if candidate_count < 1:
        return []
    from rdkit import Chem

    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    atom_count = len(symbols)
    donors, acceptors = _interaction_sites(molecule)
    if not donors or not acceptors:
        return []
    seeds = generate_stratified_dimer_poses(
        smiles, candidate_count=max(candidate_count * 2, candidate_count + 8),
        seed=seed, monomer_xyz=monomer_xyz,
    )
    base_coordinates = np.asarray(seeds[0]["coordinates"][:atom_count], dtype=float)
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in molecule.GetBonds()]
    poses: list[dict] = []
    cyclic_target = int(round(candidate_count * cyclic_fraction))
    for cycle_index in range(cyclic_target):
        donor, hydrogen = donors[cycle_index % len(donors)]
        acceptor = acceptors[cycle_index % len(acceptors)]
        local_span = base_coordinates[hydrogen] - base_coordinates[acceptor]
        if np.linalg.norm(local_span) < 1e-8:
            continue
        distance = (1.70, 1.90, 2.15, 2.35)[cycle_index % 4]
        side = float(np.linalg.norm(local_span)) + distance
        vertices = np.asarray([
            [0.0, side / math.sqrt(3.0), 0.0],
            [-side / 2.0, -side / (2.0 * math.sqrt(3.0)), 0.0],
            [side / 2.0, -side / (2.0 * math.sqrt(3.0)), 0.0],
        ])
        molecules = []
        for molecule_index in range(3):
            edge = vertices[(molecule_index + 1) % 3] - vertices[molecule_index]
            rotation = _rotation_between(local_span, edge)
            placed = (base_coordinates - base_coordinates[acceptor]) @ rotation.T
            framework = _axis_rotation(edge, math.radians((cycle_index % 3) * 120.0))
            placed = placed @ framework.T + vertices[molecule_index]
            molecules.append(placed)
        combined = np.vstack(molecules)
        interactions = []
        valid_ring = True
        for molecule_index in range(3):
            next_index = (molecule_index + 1) % 3
            donor_atom = molecule_index * atom_count + donor
            donor_hydrogen = molecule_index * atom_count + hydrogen
            acceptor_atom = next_index * atom_count + acceptor
            left = combined[donor_atom] - combined[donor_hydrogen]
            right = combined[acceptor_atom] - combined[donor_hydrogen]
            denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
            if denominator <= 0:
                valid_ring = False
                break
            angle = math.degrees(math.acos(np.clip(np.dot(left, right) / denominator, -1.0, 1.0)))
            actual_distance = float(np.linalg.norm(right))
            if angle < 110.0 or not 1.4 <= actual_distance <= 2.7:
                valid_ring = False
                break
            interactions.append({
                "donor_atom": donor_atom, "donor_hydrogen": donor_hydrogen,
                "acceptor_atom": acceptor_atom, "orientation_atoms": None,
                "target_geometry": {
                    "h_bond_distance_angstrom": actual_distance,
                    "donor_h_acceptor_angle_degrees": angle,
                },
            })
        pair_minima = [
            float(np.linalg.norm(molecules[left][:, None, :] - molecules[right][None, :, :], axis=2).min())
            for left, right in ((0, 1), (0, 2), (1, 2))
        ]
        if not valid_ring or min(pair_minima) < 0.65:
            continue
        pose_index = len(poses) + 1
        poses.append({
            "candidate_id": f"trimer-cyclic-{pose_index:04d}",
            "symbols": symbols * 3, "coordinates": combined,
            "monomer_atom_count": atom_count,
            "bonds": [(left + i * atom_count, right + i * atom_count) for i in range(3) for left, right in bonds],
            "donor_atom": interactions[0]["donor_atom"],
            "donor_hydrogen": interactions[0]["donor_hydrogen"],
            "acceptor_atom": interactions[0]["acceptor_atom"],
            "orientation_atoms": None, "target_geometry": interactions[0]["target_geometry"],
            "site_identity": {"topology": "cyclic_trimer", "interaction_sites": [[donor, hydrogen, acceptor]] * 3},
            "cluster_size": 3, "topology": "cyclic_trimer",
            "molecule_atom_ranges": [[i * atom_count, (i + 1) * atom_count] for i in range(3)],
            "interactions": interactions,
            "local_stretch_bonds": [
                {"heavy_atom": i * atom_count + local_donor,
                 "hydrogen_atom": i * atom_count + local_hydrogen,
                 "bond_class": f"{symbols[local_donor]}-H:{local_donor}-{local_hydrogen}",
                 "molecule_index": i}
                for i in range(3) for local_donor, local_hydrogen in donors
            ],
        })
    for seed_index, dimer in enumerate(seeds):
        first = np.asarray(dimer["coordinates"][:atom_count])
        second = np.asarray(dimer["coordinates"][atom_count:2 * atom_count])
        first_centroid, second_centroid = first.mean(axis=0), second.mean(axis=0)
        outward = second_centroid - first_centroid
        if np.linalg.norm(outward) < 1e-8:
            continue
        local_donor, local_hydrogen = donors[seed_index % len(donors)]
        local_acceptor = acceptors[(seed_index // len(donors)) % len(acceptors)]
        distance = (1.65, 1.85, 2.10, 2.35, 2.50)[seed_index % 5]
        angle = (125.0, 145.0, 165.0, 180.0)[seed_index % 4]
        rotation = (0.0, 90.0, 180.0, 270.0)[seed_index % 4]
        third = _place_donor_molecule(
            base_coordinates, donor_atom=local_donor, donor_hydrogen=local_hydrogen,
            acceptor_position=second[local_acceptor], outward_axis=outward,
            distance=distance, angle_degrees=angle, rotation_degrees=rotation,
        )
        combined = np.vstack([first, second, third])
        pair_distances = [
            np.linalg.norm(combined[left * atom_count:(left + 1) * atom_count, None, :]
                           - combined[None, right * atom_count:(right + 1) * atom_count, :], axis=2)
            for left, right in ((0, 1), (0, 2), (1, 2))
        ]
        if min(float(values.min()) for values in pair_distances) < 0.65:
            continue
        interactions = [dict(dimer["interactions"][0])]
        interactions.append({
            "donor_atom": 2 * atom_count + local_donor,
            "donor_hydrogen": 2 * atom_count + local_hydrogen,
            "acceptor_atom": atom_count + local_acceptor,
            "orientation_atoms": None,
            "target_geometry": {
                "h_bond_distance_angstrom": distance,
                "donor_h_acceptor_angle_degrees": angle,
                "axis_rotation_degrees": rotation,
                "framework_rotation_degrees": rotation,
            },
        })
        topology = "linear_trimer"
        # A ring is only requested when the pre-relaxed closing contact is
        # already chemically plausible; xTB is not asked to pull fragments
        # across a large vacuum gap.
        close_donor, close_hydrogen = donors[(seed_index + 1) % len(donors)]
        close_acceptor = acceptors[(seed_index + 1) % len(acceptors)]
        close_distance = float(np.linalg.norm(first[close_hydrogen] - third[close_acceptor]))
        if sum(pose.get("topology") == "cyclic_trimer" for pose in poses) < cyclic_target and close_distance <= 3.0:
            close_angle_vector = None
            try:
                left = first[close_donor] - first[close_hydrogen]
                right = third[close_acceptor] - first[close_hydrogen]
                close_angle_vector = math.degrees(math.acos(np.clip(np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right)), -1.0, 1.0)))
            except Exception:
                pass
            if close_angle_vector is not None and close_angle_vector >= 110.0:
                topology = "cyclic_trimer"
                interactions.append({
                    "donor_atom": close_donor,
                    "donor_hydrogen": close_hydrogen,
                    "acceptor_atom": 2 * atom_count + close_acceptor,
                    "orientation_atoms": None,
                    "target_geometry": {
                        "h_bond_distance_angstrom": close_distance,
                        "donor_h_acceptor_angle_degrees": close_angle_vector,
                    },
                })
        pose_index = len(poses) + 1
        primary = interactions[0]
        poses.append({
            "candidate_id": f"trimer-{topology.replace('_trimer', '')}-{pose_index:04d}",
            "symbols": symbols * 3,
            "coordinates": combined,
            "monomer_atom_count": atom_count,
            "bonds": [
                (left + molecule_index * atom_count, right + molecule_index * atom_count)
                for molecule_index in range(3) for left, right in bonds
            ],
            "donor_atom": primary["donor_atom"],
            "donor_hydrogen": primary["donor_hydrogen"],
            "acceptor_atom": primary["acceptor_atom"],
            "orientation_atoms": primary.get("orientation_atoms"),
            "target_geometry": primary["target_geometry"],
            "site_identity": {
                "topology": topology,
                "interaction_sites": [
                    [item["donor_atom"] % atom_count, item["donor_hydrogen"] % atom_count,
                     item["acceptor_atom"] % atom_count] for item in interactions
                ],
            },
            "cluster_size": 3,
            "topology": topology,
            "molecule_atom_ranges": [[i * atom_count, (i + 1) * atom_count] for i in range(3)],
            "interactions": interactions,
            "local_stretch_bonds": [
                {
                    "heavy_atom": molecule_index * atom_count + donor,
                    "hydrogen_atom": molecule_index * atom_count + hydrogen,
                    "bond_class": f"{symbols[donor]}-H:{donor}-{hydrogen}",
                    "molecule_index": molecule_index,
                }
                for molecule_index in range(3) for donor, hydrogen in donors
            ],
        })
        if len(poses) >= candidate_count:
            break
    if len(poses) < candidate_count:
        raise RuntimeError(f"Only {len(poses)} non-colliding trimer poses were generated for {candidate_count} requested")
    return poses


def dimer_sampling_plan(smiles: str, *, phase: str, multiplicity: int, charge: int) -> dict:
    from rdkit import Chem

    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    donors, acceptors = _interaction_sites(molecule)
    eligible = phase in {"liquid", "solid"} and charge == 0 and multiplicity == 1 and bool(donors and acceptors)
    reasons = []
    if phase not in {"liquid", "solid"}:
        reasons.append("self-dimer sampling is restricted to neat condensed phases")
    if charge != 0 or multiplicity != 1:
        reasons.append("charged or open-shell self-dimers require an explicit spin/counterion model")
    if not donors or not acceptors:
        reasons.append("no matched donor/acceptor interaction sites were found")
    heavy_atoms = sum(atom.GetAtomicNum() > 1 for atom in molecule.GetAtoms())
    trimer_atom_count = 3 * molecule.GetNumAtoms()
    trimer_eligible = bool(eligible and heavy_atoms <= 4 and trimer_atom_count <= 36)
    return {
        "eligible": eligible,
        "donor_sites": len(donors),
        "acceptor_sites": len(acceptors),
        "reason_not_eligible": "; ".join(reasons) if reasons else None,
        "trimer_eligible": trimer_eligible,
        "trimer_topologies": ["linear_trimer", "cyclic_trimer"] if trimer_eligible else [],
        "monomer_heavy_atoms": heavy_atoms,
        "trimer_atom_count": trimer_atom_count,
        "maximum_trimer_heavy_atoms": 4,
        "maximum_trimer_atoms": 36,
    }


def _force_field_relax(molecule, coordinates: np.ndarray) -> tuple[np.ndarray, float]:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.Mol(molecule)
    molecule.RemoveAllConformers()
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for atom, position in enumerate(coordinates):
        conformer.SetAtomPosition(atom, tuple(float(value) for value in position))
    molecule.AddConformer(conformer, assignId=True)
    force_field = None
    properties = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant="MMFF94s")
    if properties is not None:
        force_field = AllChem.MMFFGetMoleculeForceField(
            molecule, properties, confId=0, ignoreInterfragInteractions=False,
        )
    if force_field is None:
        force_field = AllChem.UFFGetMoleculeForceField(
            molecule, confId=0, ignoreInterfragInteractions=False,
        )
    if force_field is None:
        raise RuntimeError("No RDKit force field is available for the dimer")
    force_field.Minimize(maxIts=500)
    optimized = molecule.GetConformer(0)
    positions = np.asarray([list(optimized.GetAtomPosition(atom)) for atom in range(molecule.GetNumAtoms())])
    return positions, float(force_field.CalcEnergy())


def generate_dimer_candidates(
    smiles: str,
    output_dir: Path,
    *,
    max_clusters: int = 4,
    energy_window_kcal_mol: float = 5.0,
) -> tuple[list[Path], Path]:
    """Generate a small, force-field-filtered set of directed H-bonded dimers."""
    if max_clusters < 1 or not math.isfinite(energy_window_kcal_mol) or energy_window_kcal_mol <= 0:
        raise ValueError("Dimer limits must be positive")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    molecule = Chem.AddHs(base)
    if Chem.GetFormalCharge(molecule) != 0:
        raise ValueError("Automatic self-dimers are limited to neutral molecules")
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 42
    if AllChem.EmbedMolecule(molecule, parameters) < 0:
        raise RuntimeError("Could not embed the molecule for dimer sampling")
    try:
        AllChem.MMFFOptimizeMolecule(molecule, mmffVariant="MMFF94s")
    except Exception:
        AllChem.UFFOptimizeMolecule(molecule)
    coordinates = np.asarray([list(molecule.GetConformer().GetAtomPosition(index)) for index in range(molecule.GetNumAtoms())])
    donors, acceptors = _interaction_sites(molecule)
    if not donors or not acceptors:
        raise RuntimeError("No donor/acceptor pair supports automatic self-dimer sampling")

    combined = Chem.CombineMols(molecule, molecule)
    Chem.SanitizeMol(combined)
    atom_count = molecule.GetNumAtoms()
    candidates: list[dict] = []
    target_direction = np.array([1.0, 0.0, 0.0])
    for donor_atom, donor_hydrogen in donors:
        source_direction = coordinates[donor_atom] - coordinates[donor_hydrogen]
        if np.linalg.norm(source_direction) < 1e-8:
            continue
        alignment = _rotation_between(source_direction, target_direction)
        rotated = coordinates @ alignment.T
        for acceptor in acceptors:
            for angle in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0):
                around_axis = _axis_rotation(target_direction, angle)
                oriented = (rotated - rotated[donor_hydrogen]) @ around_axis.T + rotated[donor_hydrogen]
                target_hydrogen = coordinates[acceptor] + target_direction * 1.85
                oriented += target_hydrogen - oriented[donor_hydrogen]
                pose = np.vstack([coordinates, oriented])
                inter_distances = np.linalg.norm(
                    pose[:atom_count, None, :] - pose[None, atom_count:, :], axis=2,
                )
                if float(inter_distances.min()) < 0.65:
                    continue
                try:
                    optimized, energy = _force_field_relax(combined, pose)
                except Exception:
                    continue
                centroid_distance = float(np.linalg.norm(
                    optimized[:atom_count].mean(axis=0) - optimized[atom_count:].mean(axis=0)
                ))
                inter_distances = np.linalg.norm(
                    optimized[:atom_count, None, :] - optimized[None, atom_count:, :], axis=2,
                )
                if not 1.0 < centroid_distance < 10.0 or float(inter_distances.min()) < 0.65:
                    continue
                signature = np.sort(inter_distances.ravel())
                if any(float(np.sqrt(np.mean(np.square(signature - item["signature"])))) < 0.12 for item in candidates):
                    continue
                candidates.append({
                    "coordinates": optimized,
                    "energy": energy,
                    "signature": signature,
                    "donor_atom_second_molecule": donor_atom,
                    "donor_hydrogen_second_molecule": donor_hydrogen,
                    "acceptor_atom_first_molecule": acceptor,
                    "initial_azimuth_degrees": math.degrees(angle),
                })
    if not candidates:
        raise RuntimeError("Dimer pose generation produced no non-colliding bound candidates")
    candidates.sort(key=lambda item: item["energy"])
    minimum = candidates[0]["energy"]
    selected = [item for item in candidates if item["energy"] - minimum <= energy_window_kcal_mol][:max_clusters]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    manifest_entries = []
    symbols = [atom.GetSymbol() for atom in combined.GetAtoms()]
    for index, item in enumerate(selected, start=1):
        path = output_dir / f"cluster_{index:03d}.xyz"
        lines = [str(len(symbols)), f"Self-dimer candidate; MMFF/UFF energy {item['energy']:.8f}"]
        lines.extend(
            f"{symbol} {position[0]:.8f} {position[1]:.8f} {position[2]:.8f}"
            for symbol, position in zip(symbols, item["coordinates"])
        )
        path.write_text("\n".join(lines) + "\n")
        paths.append(path)
        manifest_entries.append({
            "index": index,
            "xyz": str(path),
            "force_field_energy": item["energy"],
            "relative_force_field_energy_kcal_mol": item["energy"] - minimum,
            "donor_atom_second_molecule": item["donor_atom_second_molecule"],
            "donor_hydrogen_second_molecule": item["donor_hydrogen_second_molecule"],
            "acceptor_atom_first_molecule": item["acceptor_atom_first_molecule"],
            "initial_azimuth_degrees": item["initial_azimuth_degrees"],
        })
    manifest = output_dir.parent / "clusters.json"
    manifest.write_text(json.dumps({
        "schema_version": 1,
        "kind": "neutral_self_dimer_candidates",
        "screening_method": "RDKit MMFF94s with UFF fallback and intermolecular terms enabled",
        "population_warning": "Force-field energies filter poses; they are not bulk-liquid populations.",
        "monomer_atom_count": atom_count,
        "generated_unique_candidates": len(candidates),
        "retained_candidates": manifest_entries,
    }, indent=2, sort_keys=True) + "\n")
    return paths, manifest


def normalize_self_dimer_records(records: list[dict]) -> list[dict]:
    """Convert dimer band strengths to a per-target-molecule basis."""
    normalized = []
    fallback_weight = 1.0 / len(records) if records else None
    for record in records:
        copied = {**record}
        copied["ir_modes"] = [
            {**mode, "intensity": float(mode["intensity"]) / 2.0}
            for mode in record.get("ir_modes", [])
        ]
        copied["environment_kind"] = "neutral_self_dimer"
        if copied.get("population_weight") is None and fallback_weight is not None:
            copied["population_weight"] = fallback_weight
            copied["population_model"] = "stratified_equal_representative_weights"
            copied["population_warning"] = "Weights represent geometry strata, not liquid occupancies."
        normalized.append(copied)
    return normalized


def normalize_environment_cluster_records(records: list[dict]) -> list[dict]:
    """Normalize every cluster Hessian to one target molecule."""
    normalized = []
    fallback_weight = 1.0 / len(records) if records else None
    for record in records:
        copied = {**record}
        size = int(copied.get("cluster_size") or 2)
        copied["ir_modes"] = [
            {**mode, "intensity": float(mode["intensity"]) / size}
            for mode in record.get("ir_modes", [])
        ]
        copied["environment_kind"] = str(copied.get("topology") or "dimer")
        copied["intensity_normalization_divisor"] = size
        if copied.get("population_weight") is None and fallback_weight is not None:
            copied["population_weight"] = fallback_weight
            copied["population_model"] = "stratified_equal_representative_weights"
            copied["population_warning"] = "Weights represent geometry strata, not liquid occupancies."
        normalized.append(copied)
    return normalized
