"""Reproducible periodic-methanol and electrostatic-embedding prototypes."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from .direct_local_dft import DirectLocalDFTConfig, assemble_direct_local_dft
from .local_mode_fd import LocalModeFiniteDifferenceConfig, calculate_orca_local_modes
from .spectrum import write_spectrum_csv
from src.orca_runner import run_orca


METHANOL_MOLAR_MASS_G_MOL = 32.04186
AVOGADRO_MOL_1 = 6.02214076e23
DEFAULT_METHANOL_CHARGES = {"C": 0.145, "O": -0.683, "HO": 0.418, "H": 0.040}


@dataclass(frozen=True)
class BulkEmbeddingConfig:
    schema_version: int = 1
    temperature_K: float = 298.15
    density_g_cm3: float = 0.7866
    molecule_count: int = 216
    shell_cutoff_angstrom: float = 3.5
    embedding_cutoff_angstrom: float = 10.0
    qm_shell_policy: str = "hydrogen_bond"
    minimum_independent_seeds: int = 2
    maximum_temperature_relative_error: float = 0.03
    maximum_density_relative_error: float = 0.03
    maximum_block_center_change_cm_1: float = 5.0
    maximum_block_width_change_cm_1: float = 10.0


def methanol_box_length_angstrom(molecule_count: int, density_g_cm3: float) -> float:
    if molecule_count < 1 or not math.isfinite(density_g_cm3) or density_g_cm3 <= 0:
        raise ValueError("Periodic box needs a positive molecule count and density")
    volume_cm3 = molecule_count * METHANOL_MOLAR_MASS_G_MOL / (AVOGADRO_MOL_1 * density_g_cm3)
    return float(np.cbrt(volume_cm3 * 1.0e24))


def _rotation(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(3, 3))
    q, _ = np.linalg.qr(matrix)
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def build_periodic_methanol_box(
    output_dir: Path, *, config: BulkEmbeddingConfig | None = None, seed: int = 1,
) -> dict:
    """Write a deterministic initial periodic box at the experimental density."""
    resolved = config or BulkEmbeddingConfig()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    length = methanol_box_length_angstrom(resolved.molecule_count, resolved.density_g_cm3)
    # C, O, hydroxyl H, then three methyl H.  Coordinates are a compact gas
    # geometry used only as the initial condition for force-field equilibration.
    symbols = ["C", "O", "H", "H", "H", "H"]
    template = np.asarray([
        [0.000, 0.000, 0.000], [1.430, 0.000, 0.000], [1.800, 0.900, 0.000],
        [-0.360, 1.030, 0.000], [-0.360, -0.515, 0.892], [-0.360, -0.515, -0.892],
    ])
    template -= template.mean(axis=0)
    side = int(math.ceil(resolved.molecule_count ** (1.0 / 3.0)))
    spacing = length / side
    rng = np.random.default_rng(seed)
    atoms = []
    for molecule in range(resolved.molecule_count):
        ix = molecule % side
        iy = (molecule // side) % side
        iz = molecule // (side * side)
        center = (np.asarray([ix, iy, iz], dtype=float) + 0.5) * spacing
        coordinates = template @ _rotation(rng).T + center
        atoms.extend((symbol, point, molecule) for symbol, point in zip(symbols, coordinates))
    xyz_path = output_dir / "periodic-methanol-box.xyz"
    rows = [str(len(atoms)), f"box={length:.10f},{length:.10f},{length:.10f} Angstrom seed={seed}"]
    rows.extend(f"{symbol} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}" for symbol, point, _ in atoms)
    xyz_path.write_text("\n".join(rows) + "\n")
    manifest = {
        "schema_version": 1, "kind": "periodic_methanol_box",
        "seed": seed, "temperature_K": resolved.temperature_K,
        "target_density_g_cm3": resolved.density_g_cm3,
        "molecule_count": resolved.molecule_count, "atoms_per_molecule": 6,
        "box_vectors_angstrom": [[length, 0, 0], [0, length, 0], [0, 0, length]],
        "initial_geometry": str(xyz_path),
        "required_next_stage": "nvt_equilibration_and_sampling",
    }
    path = output_dir / "periodic-box.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    manifest["artifact"] = str(path)
    return manifest


def statistical_inefficiency(values: Iterable[float]) -> dict:
    data = np.asarray(list(values), dtype=float)
    if data.ndim != 1 or len(data) < 2 or not np.all(np.isfinite(data)):
        return {"integrated_autocorrelation_time_frames": None, "effective_sample_size": 0.0}
    centered = data - data.mean()
    variance = float(np.dot(centered, centered))
    if variance <= 0:
        return {"integrated_autocorrelation_time_frames": 1.0, "effective_sample_size": float(len(data))}
    # FFT autocorrelation avoids the quadratic cost of long MD observables.
    fft_size = 1 << (2 * len(data) - 1).bit_length()
    transformed = np.fft.rfft(centered, n=fft_size)
    correlations = np.fft.irfft(transformed * transformed.conjugate(), n=fft_size)[:len(data)] / variance
    positive = []
    for value in correlations[1:]:
        if value <= 0:
            break
        positive.append(float(value))
    tau = 1.0 + 2.0 * sum(positive)
    return {
        "integrated_autocorrelation_time_frames": tau,
        "effective_sample_size": float(len(data) / tau),
    }


def validate_nvt_trajectory(
    records: list[dict], *, config: BulkEmbeddingConfig | None = None,
) -> dict:
    """Validate temperature, density, energy, H bonds, and their sampling ESS."""
    resolved = config or BulkEmbeddingConfig()
    keys = ("temperature_K", "density_g_cm3", "potential_energy", "h_bonds_per_molecule")
    series = {key: np.asarray([float(row[key]) for row in records], dtype=float) for key in keys} if records else {}
    statistics = {
        key: {
            "mean": float(values.mean()), "stddev": float(values.std()),
            **statistical_inefficiency(values),
        } for key, values in series.items()
    }
    temperature_ok = bool(records) and abs(statistics["temperature_K"]["mean"] - resolved.temperature_K) / resolved.temperature_K <= resolved.maximum_temperature_relative_error
    density_ok = bool(records) and abs(statistics["density_g_cm3"]["mean"] - resolved.density_g_cm3) / resolved.density_g_cm3 <= resolved.maximum_density_relative_error
    energy_ok = bool(records) and np.isfinite(series["potential_energy"]).all() and abs(np.polyfit(np.arange(len(records)), series["potential_energy"], 1)[0]) <= max(float(series["potential_energy"].std()), 1.0) / max(len(records), 1)
    hbonds_ok = bool(records) and statistics["h_bonds_per_molecule"]["mean"] > 0
    gates = {"temperature": temperature_ok, "density": density_ok, "energy_stationarity": energy_ok, "hydrogen_bond_statistics": hbonds_ok}
    return {"kind": "nvt_trajectory_validation", "status": "validated" if all(gates.values()) else "failed_closed", "frames": len(records), "statistics": statistics, "gates": gates}


def minimum_image(displacement: np.ndarray, box_length: float) -> np.ndarray:
    return displacement - box_length * np.round(displacement / box_length)


def extract_central_solvation_shell(
    symbols: list[str], coordinates: np.ndarray, *, central_molecule: int,
    box_length_angstrom: float, config: BulkEmbeddingConfig | None = None,
) -> dict:
    """Split a six-atom/molecule frame into central+first-shell QM and MM atoms."""
    resolved = config or BulkEmbeddingConfig()
    coordinates = np.asarray(coordinates, dtype=float)
    if len(symbols) != len(coordinates) or len(symbols) % 6:
        raise ValueError("Methanol shell extraction requires contiguous six-atom molecules")
    molecule_count = len(symbols) // 6
    if not 0 <= central_molecule < molecule_count:
        raise ValueError("Central molecule is outside the periodic frame")
    central_oxygen = coordinates[central_molecule * 6 + 1]
    central_hydroxyl_h = coordinates[central_molecule * 6 + 2]
    distances = []
    image_shifts = {}
    hydrogen_bonds = []

    def angle_at_hydrogen(donor_o: np.ndarray, hydrogen: np.ndarray, acceptor_o: np.ndarray) -> float:
        left = minimum_image(donor_o - hydrogen, box_length_angstrom)
        right = minimum_image(acceptor_o - hydrogen, box_length_angstrom)
        denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
        if denominator <= 1.0e-12:
            return 180.0
        cosine = float(np.clip(np.dot(left, right) / denominator, -1.0, 1.0))
        return math.degrees(math.acos(cosine))

    for molecule in range(molecule_count):
        oxygen = coordinates[molecule * 6 + 1]
        displacement = minimum_image(oxygen - central_oxygen, box_length_angstrom)
        distance = float(np.linalg.norm(displacement))
        distances.append(distance)
        image_shifts[str(molecule)] = (central_oxygen + displacement - oxygen).tolist()
        if molecule == central_molecule or distance > resolved.shell_cutoff_angstrom:
            continue
        acceptor_o = oxygen
        neighbor_h = coordinates[molecule * 6 + 2]
        central_donor_distance = float(np.linalg.norm(minimum_image(acceptor_o - central_hydroxyl_h, box_length_angstrom)))
        neighbor_donor_distance = float(np.linalg.norm(minimum_image(central_oxygen - neighbor_h, box_length_angstrom)))
        central_angle = angle_at_hydrogen(central_oxygen, central_hydroxyl_h, acceptor_o)
        neighbor_angle = angle_at_hydrogen(acceptor_o, neighbor_h, central_oxygen)
        if central_donor_distance <= 2.6 and central_angle >= 140.0:
            hydrogen_bonds.append({"donor_molecule": central_molecule, "acceptor_molecule": molecule, "h_acceptor_distance_angstrom": central_donor_distance, "angle_degrees": central_angle})
        if neighbor_donor_distance <= 2.6 and neighbor_angle >= 140.0:
            hydrogen_bonds.append({"donor_molecule": molecule, "acceptor_molecule": central_molecule, "h_acceptor_distance_angstrom": neighbor_donor_distance, "angle_degrees": neighbor_angle})
    if resolved.qm_shell_policy == "hydrogen_bond":
        selected = {central_molecule}
        for bond in hydrogen_bonds:
            selected.update((bond["donor_molecule"], bond["acceptor_molecule"]))
        shell_definition = "central_methanol_plus_geometric_first_hydrogen_bond_shell"
    elif resolved.qm_shell_policy == "radial_first_shell":
        selected = {index for index, distance in enumerate(distances) if distance <= resolved.shell_cutoff_angstrom}
        shell_definition = "central_methanol_plus_all_oxygen_neighbors_within_radial_first_shell"
    else:
        raise ValueError(f"Unknown QM shell policy: {resolved.qm_shell_policy}")
    qm_molecules = [central_molecule] + sorted(selected - {central_molecule})
    embedding_molecules = [
        index for index, distance in enumerate(distances)
        if index not in qm_molecules and distance <= resolved.embedding_cutoff_angstrom
    ]
    return {
        "central_molecule": central_molecule, "qm_molecules": qm_molecules,
        "box_length_angstrom": float(box_length_angstrom),
        "embedding_molecules": embedding_molecules,
        "qm_atom_indices": [index for molecule in qm_molecules for index in range(molecule * 6, molecule * 6 + 6)],
        "embedding_atom_indices": [index for molecule in embedding_molecules for index in range(molecule * 6, molecule * 6 + 6)],
        "oxygen_distances_angstrom": distances,
        "hydrogen_bonds": hydrogen_bonds,
        "shell_definition": shell_definition,
        "molecule_image_shifts_angstrom": image_shifts,
    }


def write_electrostatic_embedding_input(
    output_dir: Path, symbols: list[str], coordinates: np.ndarray, shell: dict, *,
    method_keywords: Iterable[str] = ("B3LYP", "def2-SVP", "EnGrad", "NoAutoStart"),
    ncores: int = 8,
) -> dict:
    """Write reproducible QM geometry, point charges, and an ORCA EnGrad input."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    coordinates = np.asarray(coordinates, dtype=float)
    qm = list(shell["qm_atom_indices"])
    mm = list(shell["embedding_atom_indices"])
    box_length = float(shell["box_length_angstrom"])

    def imaged_point(atom: int) -> np.ndarray:
        molecule = atom // 6
        oxygen = coordinates[molecule * 6 + 1]
        # XTC frames can wrap individual atoms across a box face. Rebuild each
        # molecule around its oxygen before moving that oxygen into the image
        # nearest the central QM molecule.
        rebuilt = oxygen + minimum_image(coordinates[atom] - oxygen, box_length)
        return rebuilt + np.asarray(
            (shell.get("molecule_image_shifts_angstrom") or {}).get(str(molecule), [0.0, 0.0, 0.0])
        )
    point_charge_path = output_dir / "embedding.pc"
    charge_rows = [str(len(mm))]
    for atom in mm:
        local = atom % 6
        charge = DEFAULT_METHANOL_CHARGES["C" if local == 0 else "O" if local == 1 else "HO" if local == 2 else "H"]
        point = imaged_point(atom)
        charge_rows.append(f"{charge:.8f} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}")
    point_charge_path.write_text("\n".join(charge_rows) + "\n")
    xyz_path = output_dir / "qm-region.xyz"
    xyz_rows = [str(len(qm)), "central methanol plus first hydrogen-bond shell"]
    for index in qm:
        point = imaged_point(index)
        xyz_rows.append(f"{symbols[index]} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}")
    xyz_path.write_text("\n".join(xyz_rows) + "\n")
    input_path = output_dir / "embedded-gradient.inp"
    input_path.write_text(
        f"! {' '.join(method_keywords)}\n%pal nprocs {int(ncores)} end\n%pointcharges \"{point_charge_path.name}\"\n* xyzfile 0 1 {xyz_path.name}\n"
    )
    manifest = {
        "schema_version": 1, "kind": "electrostatically_embedded_orca_input",
        "qm_atoms": len(qm), "point_charges": len(mm), "ncores": int(ncores),
        "qm_geometry": str(xyz_path), "point_charge_file": str(point_charge_path),
        "orca_input": str(input_path), "charge_model": DEFAULT_METHANOL_CHARGES,
        "shell": shell,
    }
    path = output_dir / "embedding-input.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    manifest["artifact"] = str(path)
    return manifest


def write_shell_extraction_artifact(
    output_dir: Path, frames: list[dict], *, config: BulkEmbeddingConfig | None = None,
) -> dict:
    """Extract central shells from decorrelated frames and retain their identities."""
    resolved = config or BulkEmbeddingConfig()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    snapshots = []
    for position, frame in enumerate(frames, start=1):
        shell = extract_central_solvation_shell(
            list(frame["symbols"]), np.asarray(frame["coordinates"], dtype=float),
            central_molecule=int(frame.get("central_molecule", 0)),
            box_length_angstrom=float(frame["box_length_angstrom"]), config=resolved,
        )
        snapshot_dir = output_dir / f"snapshot-{position:04d}"
        embedding = write_electrostatic_embedding_input(
            snapshot_dir, list(frame["symbols"]), np.asarray(frame["coordinates"], dtype=float), shell,
            ncores=8,
        )
        snapshots.append({
            "snapshot_id": str(frame.get("snapshot_id") or f"snapshot-{position:04d}"),
            "trajectory_seed": frame.get("trajectory_seed"),
            "trajectory_block": frame.get("trajectory_block"),
            "trajectory_cluster_id": frame.get("trajectory_cluster_id"),
            "trajectory_occupancy_weight": float(frame.get("trajectory_occupancy_weight", 0.0)),
            "shell": shell, "embedding_input": embedding["artifact"],
        })
    report = {
        "schema_version": 1, "kind": "decorrelated_bulk_methanol_shell_extraction",
        "configuration": asdict(resolved), "snapshots": snapshots,
        "status": "completed" if snapshots and all(item["shell"]["hydrogen_bonds"] for item in snapshots) else "failed_closed_missing_first_hydrogen_bond_shell",
    }
    path = output_dir / "shell-extraction.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(path)
    return report


def run_embedded_local_mode_snapshot(
    snapshot: dict, *, maximum_orca_invocations: int = 4,
    method_keywords: list[str] | None = None,
    finite_difference_config: LocalModeFiniteDifferenceConfig | None = None,
    job_label: str = "embedded-local-mode-finite-differences",
) -> dict:
    """Run the central O--H projected mode with fixed retained point charges."""
    embedding = json.loads(Path(snapshot["embedding_input"]).read_text())
    xyz = Path(embedding["qm_geometry"])
    point_charge_value = embedding.get("point_charge_file")
    point_charges = Path(point_charge_value) if point_charge_value else None
    job_dir = xyz.parent / job_label
    result = calculate_orca_local_modes(
        xyz, [{
            "bond_class": "O-H", "spectral_band_class": "bulk_liquid_oh",
            "coordination_class": "embedded_first_shell", "heavy_atom": 1,
            "hydrogen_atom": 2, "molecule_index": 0,
        }], job_dir, hard_orca_invocation_cap=maximum_orca_invocations,
        ncores=8, method_keywords=method_keywords,
        config=finite_difference_config,
        point_charge_file=point_charges,
    )
    return {
        **snapshot, "status": result["status"],
        "local_mode_output": str(job_dir / "local-modes.json"),
    }


def optimize_embedded_qm_shell(
    snapshot: dict, *, ncores: int = 8,
    method_keywords: Iterable[str] = ("B3LYP", "def2-SVP", "Opt", "TightSCF", "NoAutoStart"),
    constrain_heavy_atoms: bool = True,
    job_label: str = "embedded-qm-relaxation",
) -> dict:
    """Relax the extracted QM first shell in its frozen electrostatic environment."""
    embedding = json.loads(Path(snapshot["embedding_input"]).read_text())
    source_xyz = Path(embedding["qm_geometry"])
    point_charge_value = embedding.get("point_charge_file")
    source_pc = Path(point_charge_value) if point_charge_value else None
    output_dir = source_xyz.parent / job_label
    output_dir.mkdir(parents=True, exist_ok=True)
    xyz = output_dir / "input.xyz"
    point_charges = output_dir / "embedding.pc"
    if not xyz.exists():
        xyz.write_bytes(source_xyz.read_bytes())
    if source_pc is not None and not point_charges.exists():
        point_charges.write_bytes(source_pc.read_bytes())
    inp = output_dir / "embedded-opt.inp"
    xyz_rows = source_xyz.read_text().splitlines()[2:]
    heavy_atoms = [index for index, row in enumerate(xyz_rows) if row.split()[0] not in {"H", "D"}]
    geometry_block = ""
    if constrain_heavy_atoms:
        constraints = "\n".join(f"  {{ C {index} C }}" for index in heavy_atoms)
        geometry_block = f"%geom\n Constraints\n{constraints}\n end\nend\n"
    point_charge_block = f"%pointcharges \"{point_charges.name}\"\n" if source_pc is not None else ""
    inp.write_text(
        f"! {' '.join(method_keywords)}\n%pal nprocs {int(ncores)} end\n"
        f"{point_charge_block}{geometry_block}* xyzfile 0 1 {xyz.name}\n"
    )
    run_orca(inp)
    optimized = inp.with_suffix(".xyz")
    if not optimized.is_file():
        raise RuntimeError("Embedded ORCA optimization produced no final XYZ geometry")
    relaxed_embedding = {**embedding, "qm_geometry": str(optimized), "optimization_input": str(inp), "heavy_atoms_constrained": constrain_heavy_atoms}
    manifest = output_dir / "relaxed-embedding-input.json"
    manifest.write_text(json.dumps(relaxed_embedding, indent=2, sort_keys=True) + "\n")
    return {**snapshot, "embedding_input": str(manifest), "relaxation_status": "completed"}


def write_embedded_local_dft_ensemble(
    run_dir: Path, snapshot_results: list[dict], *,
    isolated_cluster_baseline: Path | None = None,
) -> dict:
    """Assemble an occupancy-weighted embedded O--H spectrum and attribution baseline."""
    conformers = []
    for result in snapshot_results:
        embedding = json.loads(Path(result["embedding_input"]).read_text())
        conformers.append({
            "optimized_xyz": embedding["qm_geometry"],
            "frequency_output": str(Path(embedding["orca_input"]).with_suffix(".out")),
            "local_mode_output": result.get("local_mode_output"),
            "local_stretch_bonds": [{
                "bond_class": "O-H", "spectral_band_class": "bulk_liquid_oh",
                "coordination_class": "embedded_first_shell", "heavy_atom": 1,
                "hydrogen_atom": 2,
            }],
            "population_weight": result.get("trajectory_occupancy_weight", 0.0),
            "cluster_size": 1, "independent_environment_id": result["snapshot_id"],
        })
    grid, spectrum, modes, provenance = assemble_direct_local_dft(
        conformers, config=DirectLocalDFTConfig(expected_valid_local_modes=len(conformers)),
    )
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    spectrum_path = write_spectrum_csv(
        run_dir / "spectrum_embedded_direct_local_dft_intrinsic.csv", grid, spectrum,
        "absorption_strength",
    )
    valid = [item for item in modes if item.get("validation_status") == "validated"]
    report = {
        "schema_version": 1, "kind": "embedded_bulk_liquid_direct_local_dft_ensemble",
        "status": "completed" if len(valid) == len(conformers) and conformers else "failed_closed",
        "population_model": "decorrelated_bulk_trajectory_cluster_occupancy",
        "provenance": provenance, "modes": modes,
        "isolated_cluster_baseline": str(isolated_cluster_baseline or ""),
        "artifacts": {"intrinsic_spectrum": str(spectrum_path)},
    }
    path = run_dir / "embedded-local-dft-ensemble.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(path)
    return report


def evaluate_bulk_embedding_acceptance(
    seed_reports: list[dict], *, config: BulkEmbeddingConfig | None = None,
) -> dict:
    """Require independent-seed, held-out-block, convergence, and FD evidence."""
    resolved = config or BulkEmbeddingConfig()
    valid_trajectories = [item for item in seed_reports if (item.get("trajectory_validation") or {}).get("status") == "validated"]
    finite_difference = all(bool(item.get("finite_difference_validation_passed")) for item in seed_reports)
    heldout = all(bool(item.get("held_out_gain_reproduced")) for item in seed_reports)
    block_convergence = all(
        float(item.get("maximum_block_center_change_cm-1", math.inf)) <= resolved.maximum_block_center_change_cm_1
        and float(item.get("maximum_block_width_change_cm-1", math.inf)) <= resolved.maximum_block_width_change_cm_1
        for item in seed_reports
    )
    gates = {
        "stable_trajectories": len(valid_trajectories) == len(seed_reports) and bool(seed_reports),
        "at_least_two_independent_seeds": len({item.get("seed") for item in seed_reports}) >= resolved.minimum_independent_seeds,
        "block_center_and_width_convergence": block_convergence,
        "embedded_finite_difference_validation": finite_difference,
        "gain_on_held_out_trajectory_blocks": heldout,
    }
    return {
        "schema_version": 1, "kind": "explicit_bulk_liquid_embedding_acceptance",
        "configuration": asdict(resolved), "seed_reports": seed_reports,
        "acceptance_gates": gates,
        "status": "accepted" if all(gates.values()) else "failed_closed",
    }


def write_bulk_embedding_acceptance(
    run_dir: Path, seed_reports: list[dict], *, config: BulkEmbeddingConfig | None = None,
) -> dict:
    report = evaluate_bulk_embedding_acceptance(seed_reports, config=config)
    path = Path(run_dir) / "bulk-embedding-acceptance.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(path)
    return report
