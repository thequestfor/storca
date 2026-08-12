"""Diversity clustering and representative selection for sampled environments."""

from __future__ import annotations

import json
import math
from pathlib import Path
import shutil

import numpy as np


SELECTION_SCHEMA_VERSION = 1
FIDELITY_ENVIRONMENT_TARGETS = {"fast": 0, "auto": 4, "balanced": 6, "accurate": 10}


def allocate_orca_budget(
    *,
    fidelity: str,
    max_conformers: int,
    max_environments: int,
    max_orca_jobs: int,
    environment_eligible: bool,
) -> dict:
    """Reserve environment jobs first while retaining at least one monomer."""
    fidelity = fidelity.strip().lower()
    if fidelity not in FIDELITY_ENVIRONMENT_TARGETS:
        raise ValueError("Unsupported environment fidelity")
    if min(max_conformers, max_environments, max_orca_jobs) < 1:
        raise ValueError("ORCA budget limits must be positive")
    requested = FIDELITY_ENVIRONMENT_TARGETS[fidelity] if environment_eligible else 0
    environment_jobs = min(requested, max_environments, max(0, max_orca_jobs - 1))
    monomer_jobs = min(max_conformers, max_orca_jobs - environment_jobs)
    return {
        "schema_version": 1,
        "fidelity": fidelity,
        "environment_eligible": bool(environment_eligible),
        "requested_environment_jobs": requested,
        "reserved_environment_jobs": environment_jobs,
        "maximum_monomer_jobs": monomer_jobs,
        "maximum_total_orca_jobs": max_orca_jobs,
        "minimum_environment_gate_reachable": environment_jobs >= 3,
        "environment_gate_warning": (
            None if environment_jobs >= 3 or requested == 0
            else "The hard ORCA budget cannot provide three independent environments."
        ),
    }


def _read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = Path(path).read_text().splitlines()
    try:
        atom_count = int(lines[0].strip())
    except (IndexError, ValueError) as error:
        raise ValueError(f"Invalid XYZ file: {path}") from error
    rows = [line.split() for line in lines[2:2 + atom_count]]
    if len(rows) != atom_count or any(len(row) < 4 for row in rows):
        raise ValueError(f"Incomplete XYZ file: {path}")
    symbols = [row[0] for row in rows]
    try:
        coordinates = np.asarray([[float(value) for value in row[1:4]] for row in rows], dtype=float)
    except ValueError as error:
        raise ValueError(f"Non-numeric XYZ coordinates: {path}") from error
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"Non-finite XYZ coordinates: {path}")
    return symbols, coordinates


def aligned_heavy_atom_rmsd(left_path: Path, right_path: Path) -> float:
    """Return atom-mapped heavy-atom RMSD after proper Kabsch alignment."""
    left_symbols, left = _read_xyz(left_path)
    right_symbols, right = _read_xyz(right_path)
    if left_symbols != right_symbols or left.shape != right.shape:
        raise ValueError("Environment RMSD requires identical atom identities and ordering")
    indices = [index for index, symbol in enumerate(left_symbols) if symbol != "H"]
    if not indices:
        indices = list(range(len(left_symbols)))
    reference = left[indices] - left[indices].mean(axis=0)
    candidate = right[indices] - right[indices].mean(axis=0)
    covariance = candidate.T @ reference
    u_matrix, _, v_matrix = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u_matrix @ v_matrix))
    aligned = candidate @ (u_matrix @ correction @ v_matrix)
    return float(np.sqrt(np.mean(np.square(aligned - reference))))


def _wrapped_angle_gap(left: float | None, right: float | None) -> float:
    if left is None and right is None:
        return 0.0
    if left is None or right is None:
        return 180.0
    difference = abs(float(left) - float(right)) % 360.0
    return min(difference, 360.0 - difference)


def association_class(record: dict) -> str:
    distance = float(record["environment_features"]["h_bond_distance_angstrom"])
    if distance < 1.85:
        return "strong"
    if distance <= 2.20:
        return "intermediate"
    return "weak"


def _pair_components(left: dict, right: dict) -> dict:
    left_features = left["environment_features"]
    right_features = right["environment_features"]
    same_topology = (
        int(left.get("cluster_size", 2)) == int(right.get("cluster_size", 2))
        and left.get("topology", "dimer") == right.get("topology", "dimer")
    )
    return {
        "h_bond_distance_angstrom": abs(
            float(left_features["h_bond_distance_angstrom"])
            - float(right_features["h_bond_distance_angstrom"])
        ),
        "donor_h_acceptor_angle_degrees": abs(
            float(left_features["donor_h_acceptor_angle_degrees"])
            - float(right_features["donor_h_acceptor_angle_degrees"])
        ),
        "donor_acceptor_distance_angstrom": abs(
            float(left_features["donor_acceptor_distance_angstrom"])
            - float(right_features["donor_acceptor_distance_angstrom"])
        ),
        "intermolecular_orientation_degrees": _wrapped_angle_gap(
            left_features.get("intermolecular_orientation_degrees"),
            right_features.get("intermolecular_orientation_degrees"),
        ),
        "heavy_atom_rmsd_angstrom": (
            aligned_heavy_atom_rmsd(Path(left["optimized_xyz"]), Path(right["optimized_xyz"]))
            if same_topology else abs(
                float(left_features.get("heavy_atom_rmsd_angstrom", 0.0))
                - float(right_features.get("heavy_atom_rmsd_angstrom", 0.0))
            )
        ),
        "relative_xtb_energy_kcal_mol": min(
            10.0,
            abs(
                float(left.get("relative_xtb_energy_kcal_mol") or 0.0)
                - float(right.get("relative_xtb_energy_kcal_mol") or 0.0)
            ),
        ),
        "estimated_local_frequency_cm-1": (
            abs(
                float(left_features["estimated_local_frequency_cm-1"])
                - float(right_features["estimated_local_frequency_cm-1"])
            )
            if all(
                isinstance(features.get("estimated_local_frequency_cm-1"), (int, float))
                and math.isfinite(float(features["estimated_local_frequency_cm-1"]))
                for features in (left_features, right_features)
            ) else 0.0
        ),
        "same_interaction_sites": left.get("site_identity") == right.get("site_identity"),
        "same_topology": same_topology,
    }


def environment_distance(left: dict, right: dict) -> tuple[float, dict]:
    components = _pair_components(left, right)
    squared = (
        (components["h_bond_distance_angstrom"] / 0.15) ** 2
        + (components["donor_h_acceptor_angle_degrees"] / 15.0) ** 2
        + (components["donor_acceptor_distance_angstrom"] / 0.20) ** 2
        + (components["intermolecular_orientation_degrees"] / 30.0) ** 2
        + (components["heavy_atom_rmsd_angstrom"] / 0.35) ** 2
        + 0.25 * (components["relative_xtb_energy_kcal_mol"] / 3.0) ** 2
        + 0.50 * (components["estimated_local_frequency_cm-1"] / 50.0) ** 2
        + (0.0 if components["same_interaction_sites"] else 4.0)
        + (0.0 if components["same_topology"] else 9.0)
    )
    return math.sqrt(squared), components


def pairwise_environment_distances(records: list[dict]) -> tuple[np.ndarray, dict[tuple[int, int], dict]]:
    size = len(records)
    matrix = np.zeros((size, size), dtype=float)
    components: dict[tuple[int, int], dict] = {}
    for left in range(size):
        for right in range(left + 1, size):
            distance, detail = environment_distance(records[left], records[right])
            matrix[left, right] = matrix[right, left] = distance
            components[(left, right)] = detail
    return matrix, components


def _is_duplicate(components: dict) -> bool:
    return bool(
        components["same_interaction_sites"]
        and components["same_topology"]
        and components["h_bond_distance_angstrom"] <= 0.08
        and components["donor_h_acceptor_angle_degrees"] <= 8.0
        and components["donor_acceptor_distance_angstrom"] <= 0.10
        and components["intermolecular_orientation_degrees"] <= 15.0
        and components["heavy_atom_rmsd_angstrom"] <= 0.20
        and components["estimated_local_frequency_cm-1"] <= 20.0
    )


def deduplicate_environments(
    records: list[dict], matrix: np.ndarray, components: dict[tuple[int, int], dict],
) -> tuple[list[int], list[dict]]:
    """Collapse connected duplicate groups and retain their geometry medoids."""
    parents = list(range(len(records)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parents[max(left_root, right_root)] = min(left_root, right_root)

    for pair, detail in components.items():
        if _is_duplicate(detail):
            union(*pair)
    groups: dict[int, list[int]] = {}
    for index in range(len(records)):
        groups.setdefault(find(index), []).append(index)
    retained: list[int] = []
    report: list[dict] = []
    for members in sorted(groups.values(), key=lambda group: min(records[index]["candidate_id"] for index in group)):
        medoid = min(
            members,
            key=lambda index: (
                float(matrix[index, members].sum()), records[index]["candidate_id"],
            ),
        )
        retained.append(medoid)
        report.append({
            "retained_candidate_id": records[medoid]["candidate_id"],
            "member_candidate_ids": sorted(records[index]["candidate_id"] for index in members),
            "duplicates_removed": len(members) - 1,
        })
    retained.sort(key=lambda index: records[index]["candidate_id"])
    return retained, report


def deterministic_k_medoids(
    records: list[dict], matrix: np.ndarray, cluster_count: int,
) -> tuple[list[int], list[int]]:
    if not 1 <= cluster_count <= len(records):
        raise ValueError("Invalid environment cluster count")
    medoids: list[int] = []
    seeded_topologies: list[str] = []
    for label in ("dimer", "linear_trimer", "cyclic_trimer"):
        members = [index for index, record in enumerate(records) if record.get("topology", "dimer") == label]
        if members and len(medoids) < cluster_count:
            medoids.append(min(
                members,
                key=lambda index: (float(matrix[index, members].sum()), records[index]["candidate_id"]),
            ))
            seeded_topologies.append(label)
    for label in ("strong", "intermediate", "weak"):
        members = [index for index, record in enumerate(records) if association_class(record) == label]
        members = [index for index in members if index not in medoids]
        if members and len(medoids) < cluster_count:
            medoids.append(min(
                members,
                key=lambda index: (float(matrix[index, members].sum()), records[index]["candidate_id"]),
            ))
    while len(medoids) < cluster_count:
        remaining = [index for index in range(len(records)) if index not in medoids]
        medoids.append(max(
            remaining,
            key=lambda index: (
                min(float(matrix[index, medoid]) for medoid in medoids) if medoids else math.inf,
                records[index]["candidate_id"],
            ),
        ))
    assignments = [-1] * len(records)
    for _ in range(100):
        new_assignments = [
            min(
                range(len(medoids)),
                key=lambda cluster: (float(matrix[index, medoids[cluster]]), records[medoids[cluster]]["candidate_id"]),
            )
            for index in range(len(records))
        ]
        new_medoids = []
        for cluster in range(len(medoids)):
            members = [index for index, assignment in enumerate(new_assignments) if assignment == cluster]
            if cluster < len(seeded_topologies):
                topology_members = [
                    index for index in members
                    if records[index].get("topology", "dimer") == seeded_topologies[cluster]
                ]
                if topology_members:
                    members = topology_members
            new_medoids.append(min(
                members,
                key=lambda index: (float(matrix[index, members].sum()), records[index]["candidate_id"]),
            ))
        if new_assignments == assignments and new_medoids == medoids:
            break
        assignments, medoids = new_assignments, new_medoids
    order = sorted(range(len(medoids)), key=lambda cluster: records[medoids[cluster]]["candidate_id"])
    remap = {old: new for new, old in enumerate(order)}
    return [medoids[old] for old in order], [remap[assignment] for assignment in assignments]


def incremental_diversity_order(
    distance_matrix: np.ndarray, topology_labels: list[str] | None = None,
) -> list[int]:
    """Order representatives central-first, then by farthest-point diversity."""
    matrix = np.asarray(distance_matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] < 1:
        raise ValueError("Representative distance matrix must be nonempty and square")
    first = min(range(len(matrix)), key=lambda index: (float(matrix[index].sum()), index))
    order = [first]
    if topology_labels is not None:
        if len(topology_labels) != len(matrix):
            raise ValueError("Topology labels must match the representative matrix")
        for topology in ("dimer", "linear_trimer", "cyclic_trimer"):
            if topology == topology_labels[first] or topology not in topology_labels:
                continue
            choices = [index for index, label in enumerate(topology_labels) if label == topology and index not in order]
            if choices:
                order.append(max(choices, key=lambda index: min(float(matrix[index, selected]) for selected in order)))
    while len(order) < len(matrix):
        remaining = [index for index in range(len(matrix)) if index not in order]
        order.append(max(
            remaining,
            key=lambda index: (min(float(matrix[index, selected]) for selected in order), -index),
        ))
    return order


def _merge_environment_report(path: Path, additions: dict) -> Path:
    retained = json.loads(path.read_text()) if path.is_file() else {}
    retained["schema_version"] = max(2, int(retained.get("schema_version", 0)))
    retained.update(additions)
    path.write_text(json.dumps(retained, indent=2, sort_keys=True) + "\n")
    return path


def select_xtb_environment_representatives(
    records: list[dict], run_dir: Path, *, representative_count: int,
    frequency_records: list[dict] | None = None,
) -> tuple[list[Path], Path]:
    """Deduplicate and export coverage-constrained, diverse xTB representatives."""
    eligible = sorted(
        (
            record for record in records
            if record.get("sampling_status") == "retained"
            and record.get("environment_features")
            and record.get("optimized_xyz")
            and Path(record["optimized_xyz"]).is_file()
        ),
        key=lambda record: record["candidate_id"],
    )
    if not eligible:
        raise RuntimeError("No retained xTB environments are available for representative selection")
    full_matrix, full_components = pairwise_environment_distances(eligible)
    retained_indices, deduplication = deduplicate_environments(
        eligible, full_matrix, full_components,
    )
    unique = [eligible[index] for index in retained_indices]
    unique_matrix = full_matrix[np.ix_(retained_indices, retained_indices)]
    selected_count = min(int(representative_count), len(unique))
    if selected_count < 1:
        raise RuntimeError("The ORCA budget does not permit an environment representative")
    acquisition_report = None
    if frequency_records is None:
        frequency_manifest = Path(run_dir) / "xtb-snapshot-frequencies.json"
        if frequency_manifest.is_file():
            frequency_records = json.loads(frequency_manifest.read_text()).get("candidates") or []
    if frequency_records:
        from .environment_acquisition import select_mode_coverage_representatives
        medoids, acquisition_report = select_mode_coverage_representatives(
            unique, frequency_records, unique_matrix, selected_count,
        )
        assignments = [
            min(
                range(len(medoids)),
                key=lambda cluster: (
                    float(unique_matrix[index, medoids[cluster]]),
                    unique[medoids[cluster]]["candidate_id"],
                ),
            )
            for index in range(len(unique))
        ]
    else:
        medoids, assignments = deterministic_k_medoids(unique, unique_matrix, selected_count)
    cluster_dir = Path(run_dir) / "clusters"
    selected_dir = cluster_dir / "selected-conformers"
    selected_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    manifest_entries = []
    cluster_reports = []
    equal_weight = 1.0 / selected_count
    for cluster_index, medoid in enumerate(medoids, start=1):
        members = [index for index, assignment in enumerate(assignments) if assignment == cluster_index - 1]
        representative = unique[medoid]
        cluster_id = f"environment-cluster-{cluster_index:03d}"
        destination = selected_dir / f"cluster_{cluster_index:03d}.xyz"
        shutil.copy2(representative["optimized_xyz"], destination)
        paths.append(destination)
        features = {**representative["environment_features"], "geometry_cluster_id": cluster_id}
        entry = {
            "selected_position": cluster_index,
            "xyz": str(destination),
            "source_xtb_candidate_id": representative["candidate_id"],
            "source_xtb_xyz": representative["optimized_xyz"],
            "geometry_cluster_id": cluster_id,
            "independent_environment_id": cluster_id,
            "environment_features": features,
            "association_class": association_class(representative),
            "sampling_support_count": len(members),
            "sampling_support_candidate_ids": sorted(unique[index]["candidate_id"] for index in members),
            "population_weight": equal_weight,
            "population_model": "stratified_equal_representative_weights",
            "population_warning": "Weights represent geometry strata, not liquid occupancies.",
            "geometry_role": "environment_snapshot",
            "cluster_size": int(representative.get("cluster_size", 2)),
            "topology": representative.get("topology", "dimer"),
            "molecule_atom_ranges": representative.get("molecule_atom_ranges"),
            "local_stretch_bonds": representative.get("local_stretch_bonds", []),
            "hydrogen_bond_interactions": representative.get("hydrogen_bond_interactions", []),
        }
        manifest_entries.append(entry)
        cluster_reports.append({
            "geometry_cluster_id": cluster_id,
            "medoid_candidate_id": representative["candidate_id"],
            "association_class": entry["association_class"],
            "cluster_size": entry["cluster_size"],
            "topology": entry["topology"],
            "sampling_support_count": len(members),
            "member_candidate_ids": entry["sampling_support_candidate_ids"],
        })
    representative_matrix = unique_matrix[np.ix_(medoids, medoids)]
    execution_order = (
        list(range(len(medoids))) if acquisition_report is not None
        else incremental_diversity_order(
            representative_matrix,
            [unique[medoid].get("topology", "dimer") for medoid in medoids],
        )
    )
    paths = [paths[index] for index in execution_order]
    manifest_entries = [manifest_entries[index] for index in execution_order]
    for execution_position, entry in enumerate(manifest_entries, start=1):
        entry["selected_position"] = execution_position
        entry["execution_position"] = execution_position
    execution_by_cluster = {
        entry["geometry_cluster_id"]: entry["execution_position"] for entry in manifest_entries
    }
    for report in cluster_reports:
        report["execution_position"] = execution_by_cluster[report["geometry_cluster_id"]]
    manifest = cluster_dir / "selected-conformers.json"
    manifest.write_text(json.dumps({
        "schema_version": SELECTION_SCHEMA_VERSION,
        "kind": "xtb_diversity_selected_environment_representatives",
        "population_model": "stratified_equal_representative_weights",
        "population_warning": "Weights represent geometry strata, not liquid occupancies.",
        "execution_order_model": (
            "mode_class_deficits_then_frequency_geometry_diversity"
            if acquisition_report is not None else "central_then_farthest_point_diversity"
        ),
        "cluster_model": "mixed_dimer_trimer_network" if any(entry["cluster_size"] == 3 for entry in manifest_entries) else "dimer_only",
        "conformers": manifest_entries,
    }, indent=2, sort_keys=True) + "\n")
    compatibility_manifest = cluster_dir / "clusters.json"
    compatibility_manifest.write_text(json.dumps({
        "schema_version": 2,
        "kind": "xtb_diversity_selected_environment_representatives",
        "generated_candidates": len(records),
        "eligible_candidates": len(eligible),
        "deduplicated_candidates": len(unique),
        "retained_candidates": manifest_entries,
    }, indent=2, sort_keys=True) + "\n")
    report_path = Path(run_dir) / "environment-sampling.json"
    _merge_environment_report(report_path, {
        "deduplication": {
            "input_candidates": len(eligible),
            "unique_candidates": len(unique),
            "groups": deduplication,
        },
        "clustering": {
            "method": (
                "mode_class_coverage_constrained_diversity_acquisition"
                if acquisition_report is not None else "deterministic_k_medoids"
            ),
            "distance_model": (
                "interaction_geometry_orientation_rmsd_local_frequency_"
                "with_low_weight_xtb_energy"
            ),
            "requested_clusters": representative_count,
            "selected_clusters": selected_count,
            "clusters": cluster_reports,
        },
        "selection": {
            "manifest": str(manifest),
            "representatives": manifest_entries,
        },
        "population_model": {
            "kind": "stratified_equal_representative_weights",
            "warning": "Weights represent geometry strata, not liquid occupancies.",
            "vacuum_energies_used_as_liquid_populations": False,
        },
    })
    if acquisition_report is not None:
        from .environment_acquisition import write_environment_acquisition_report
        acquisition_report["representative_manifest"] = str(manifest)
        write_environment_acquisition_report(run_dir, acquisition_report)
    return paths, manifest
