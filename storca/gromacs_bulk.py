"""Pinned GROMACS orchestration for periodic OPLS-AA liquid methanol."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess

import numpy as np

from .bulk_embedding import BulkEmbeddingConfig, extract_central_solvation_shell, statistical_inefficiency


@dataclass(frozen=True)
class GromacsMethanolConfig:
    schema_version: int = 1
    temperature_K: float = 298.15
    pressure_bar: float = 1.0
    molecule_count: int = 216
    timestep_ps: float = 0.001
    heating_ps: float = 100.0
    npt_ps: float = 1000.0
    nvt_equilibration_ps: float = 1000.0
    production_ps: float = 5000.0
    coordinate_interval_ps: float = 0.5
    observable_interval_ps: float = 0.1
    cutoff_nm: float = 1.0


METHANOL_ATOMS = (
    ("C", "opls_157", 0.145, 12.01100),
    ("O", "opls_154", -0.683, 15.99940),
    ("HO", "opls_155", 0.418, 1.00800),
    ("H1", "opls_156", 0.040, 1.00800),
    ("H2", "opls_156", 0.040, 1.00800),
    ("H3", "opls_156", 0.040, 1.00800),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _run(command: list[str], cwd: Path, *, stdin: str | None = None) -> dict:
    result = subprocess.run(command, cwd=cwd, input=stdin, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n"
            f"{result.stdout[-4000:]}\n{result.stderr[-4000:]}"
        )
    return {"command": command, "stdout": result.stdout, "stderr": result.stderr}


def gromacs_version(gmx: Path) -> dict:
    gmx = Path(gmx).resolve()
    result = _run([str(gmx), "--version"], Path.cwd())
    text = result["stdout"] + result["stderr"]
    def field(label: str) -> str | None:
        match = re.search(rf"^{re.escape(label)}:\s*(.+)$", text, re.MULTILINE)
        return match.group(1).strip() if match else None
    return {
        "executable": str(gmx), "sha256": _sha256(gmx),
        "version": field("GROMACS version"), "precision": field("Precision"),
        "mpi_library": field("MPI library"), "openmp_support": field("OpenMP support"),
        "simd": field("SIMD instructions"), "gpu_support": field("GPU support"),
    }


def audit_gromacs_oplsaa(gmx: Path) -> dict:
    """Verify the installed OPLS-AA defaults and methanol nonbonded parameters."""
    gmx = Path(gmx).resolve()
    prefix = gmx.parents[1]
    force_field = prefix / "share" / "gromacs" / "top" / "oplsaa.ff"
    defaults_path = force_field / "forcefield.itp"
    nonbonded_path = force_field / "ffnonbonded.itp"
    bonded_path = force_field / "ffbonded.itp"
    if not all(path.is_file() for path in (defaults_path, nonbonded_path, bonded_path)):
        raise FileNotFoundError("Pinned GROMACS OPLS-AA force-field files are unavailable")
    defaults = defaults_path.read_text()
    match = re.search(r"^\s*1\s+3\s+yes\s+0\.5\s+0\.5\s*$", defaults, re.MULTILINE)
    expected = {item[1]: (item[2], item[3]) for item in METHANOL_ATOMS}
    found = {}
    for line in nonbonded_path.read_text().splitlines():
        fields = line.split()
        if fields and fields[0] in expected and len(fields) >= 8:
            found[fields[0]] = {"charge_e": float(fields[4]), "mass_u": float(fields[3]), "sigma_nm": float(fields[6]), "epsilon_kj_mol": float(fields[7])}
    unique_expected = set(expected)
    passed = bool(match) and unique_expected.issubset(found) and all(
        math.isclose(found[key]["charge_e"], expected[key][0], abs_tol=1e-12)
        and math.isclose(found[key]["mass_u"], expected[key][1], abs_tol=1e-5)
        for key in unique_expected
    )
    report = {
        "kind": "installed_gromacs_oplsaa_methanol_audit",
        "status": "validated" if passed else "failed_validation",
        "force_field_directory": str(force_field),
        "defaults": {"nbfunc": 1, "combination_rule": 3, "generate_pairs": True, "fudge_lj": 0.5, "fudge_qq": 0.5},
        "methanol_atom_types": found,
        "files": {
            "forcefield.itp": _sha256(defaults_path),
            "ffnonbonded.itp": _sha256(nonbonded_path),
            "ffbonded.itp": _sha256(bonded_path),
        },
    }
    if not passed:
        raise RuntimeError("Installed GROMACS OPLS-AA parameters differ from the frozen methanol contract")
    return report


def _read_xyz_box(path: Path) -> tuple[list[str], np.ndarray, float]:
    lines = Path(path).read_text().splitlines()
    count = int(lines[0])
    match = re.search(r"box=([0-9.eE+-]+),", lines[1])
    if not match:
        raise ValueError("Periodic methanol XYZ lacks a cubic box declaration")
    rows = [line.split() for line in lines[2:2 + count]]
    return [row[0] for row in rows], np.asarray([[float(x) for x in row[1:4]] for row in rows]), float(match.group(1))


def write_methanol_gro(xyz: Path, output: Path) -> Path:
    symbols, coordinates, box_angstrom = _read_xyz_box(xyz)
    if len(symbols) % 6 or any(
        symbols[index:index + 6] != ["C", "O", "H", "H", "H", "H"]
        for index in range(0, len(symbols), 6)
    ):
        raise ValueError("Expected contiguous C,O,HO,H,H,H methanol atom order")
    rows = ["Periodic OPLS-AA methanol", str(len(symbols))]
    names = [item[0] for item in METHANOL_ATOMS]
    for index, point in enumerate(coordinates):
        residue = index // 6 + 1
        atom = index % 6
        rows.append(
            f"{residue:5d}{'MEOH':<5}{names[atom]:>5}{index + 1:5d}"
            f"{point[0] / 10:8.3f}{point[1] / 10:8.3f}{point[2] / 10:8.3f}"
        )
    length = box_angstrom / 10.0
    rows.append(f"{length:10.5f}{length:10.5f}{length:10.5f}")
    Path(output).write_text("\n".join(rows) + "\n")
    return Path(output)


def write_oplsaa_methanol_topology(output_dir: Path, *, molecule_count: int) -> dict:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    itp = output_dir / "methanol-oplsaa.itp"
    lines = [
        "; Methanol using canonical GROMACS OPLS-AA atom and bonded types",
        "[ moleculetype ]", "; name nrexcl", "MEOH 3", "",
        "[ atoms ]", "; nr type resnr residue atom cgnr charge mass",
    ]
    for index, (name, atom_type, charge, mass) in enumerate(METHANOL_ATOMS, 1):
        lines.append(f"{index:3d} {atom_type:<10} 1 MEOH {name:<4} 1 {charge: .6f} {mass:.5f}")
    lines += [
        "", "[ bonds ]", "; ai aj funct; parameters from oplsaa.ff/ffbonded.itp",
        "1 2 1", "2 3 1", "1 4 1", "1 5 1", "1 6 1", "",
        "[ angles ]", "3 2 1 1", "4 1 2 1", "5 1 2 1", "6 1 2 1",
        "4 1 5 1", "4 1 6 1", "5 1 6 1", "",
        "[ dihedrals ]", "; alcohol torsions from OPLS-AA",
        "4 1 2 3 3", "5 1 2 3 3", "6 1 2 3 3", "",
    ]
    itp.write_text("\n".join(lines))
    top = output_dir / "topol.top"
    top.write_text(
        '#include "oplsaa.ff/forcefield.itp"\n'
        '#include "methanol-oplsaa.itp"\n\n'
        '[ system ]\nPeriodic liquid methanol\n\n'
        f'[ molecules ]\nMEOH {int(molecule_count)}\n'
    )
    return {"topology": str(top), "molecule_include": str(itp), "topology_sha256": _sha256(top), "molecule_include_sha256": _sha256(itp)}


def write_mdp_files(output_dir: Path, *, seed: int, config: GromacsMethanolConfig | None = None) -> dict:
    resolved = config or GromacsMethanolConfig()
    output_dir = Path(output_dir)
    dt = resolved.timestep_ps
    common = f"""cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = {resolved.cutoff_nm}
vdwtype = Cut-off
rvdw = {resolved.cutoff_nm}
DispCorr = EnerPres
pbc = xyz
constraints = none
dt = {dt}
nstlist = 20
nstcalcenergy = {max(1, round(resolved.observable_interval_ps / dt))}
nstenergy = {max(1, round(resolved.observable_interval_ps / dt))}
nstxout-compressed = {max(1, round(resolved.coordinate_interval_ps / dt))}
compressed-x-precision = 1000
"""
    stages = {
        "minim": "integrator = steep\nemtol = 100.0\nemstep = 0.001\nnsteps = 50000\n" + common,
        "heat": f"integrator = md\nnsteps = {round(resolved.heating_ps / dt)}\ngen-vel = yes\ngen-temp = 50\ngen-seed = {seed}\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 1.0\nref-t = {resolved.temperature_K}\nannealing = single\nannealing-npoints = 2\nannealing-time = 0 {resolved.heating_ps}\nannealing-temp = 50 {resolved.temperature_K}\npcoupl = no\n" + common,
        "npt": f"integrator = md\nnsteps = {round(resolved.npt_ps / dt)}\ncontinuation = yes\ngen-vel = no\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 1.0\nref-t = {resolved.temperature_K}\npcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 5.0\nref-p = {resolved.pressure_bar}\ncompressibility = 1.2e-4\n" + common,
        "nvt": f"integrator = md\nnsteps = {round(resolved.nvt_equilibration_ps / dt)}\ncontinuation = yes\ngen-vel = no\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 1.0\nref-t = {resolved.temperature_K}\npcoupl = no\n" + common,
        "production": f"integrator = md\nnsteps = {round(resolved.production_ps / dt)}\ncontinuation = yes\ngen-vel = no\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 1.0\nref-t = {resolved.temperature_K}\npcoupl = no\n" + common,
    }
    paths = {}
    for name, content in stages.items():
        path = output_dir / f"{name}.mdp"
        path.write_text(content)
        paths[name] = str(path)
    return paths


def prepare_bulk_seed(
    seed_dir: Path, gmx: Path, *, config: GromacsMethanolConfig | None = None,
) -> dict:
    resolved = config or GromacsMethanolConfig()
    seed_dir = Path(seed_dir)
    box = json.loads((seed_dir / "periodic-box.json").read_text())
    if int(box["molecule_count"]) != resolved.molecule_count:
        raise ValueError("Periodic box molecule count differs from MD contract")
    gro = write_methanol_gro(Path(box["initial_geometry"]), seed_dir / "initial.gro")
    topology = write_oplsaa_methanol_topology(seed_dir, molecule_count=resolved.molecule_count)
    mdp = write_mdp_files(seed_dir, seed=int(box["seed"]), config=resolved)
    _run([str(Path(gmx).resolve()), "grompp", "-f", "minim.mdp", "-c", gro.name, "-p", "topol.top", "-o", "minim.tpr", "-po", "minim-expanded.mdp", "-maxwarn", "0"], seed_dir)
    report = {
        "schema_version": 1, "kind": "gromacs_oplsaa_methanol_seed_preflight",
        "status": "prepared", "configuration": asdict(resolved), "box": box,
        "gromacs": gromacs_version(gmx), "topology": topology, "mdp": mdp,
        "installed_force_field_audit": audit_gromacs_oplsaa(gmx),
        "initial_gro": str(gro), "minimization_tpr": str(seed_dir / "minim.tpr"),
        "force_field_contract": {
            "name": "GROMACS bundled OPLS-AA",
            "atom_types": [item[1] for item in METHANOL_ATOMS],
            "charges_e": [item[2] for item in METHANOL_ATOMS],
            "net_molecular_charge_e": sum(item[2] for item in METHANOL_ATOMS),
            "combination_rule": 3, "fudge_lj_1_4": 0.5, "fudge_qq_1_4": 0.5,
        },
    }
    path = seed_dir / "gromacs-preflight.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(path)
    return report


def run_stage(
    seed_dir: Path,
    gmx: Path,
    stage: str,
    *,
    threads: int = 8,
    accelerator: str = "auto",
    resume: bool = True,
) -> dict:
    seed_dir = Path(seed_dir)
    predecessors = {"minim": "initial.gro", "heat": "minim.gro", "npt": "heat.gro", "nvt": "npt.gro", "production": "nvt.gro"}
    checkpoints = {"npt": "heat.cpt", "nvt": "npt.cpt", "production": "nvt.cpt"}
    command = [str(Path(gmx).resolve()), "grompp", "-f", f"{stage}.mdp", "-c", predecessors[stage], "-p", "topol.top", "-o", f"{stage}.tpr", "-po", f"{stage}-expanded.mdp", "-maxwarn", "0"]
    if stage in checkpoints:
        command += ["-t", checkpoints[stage]]
    _run(command, seed_dir)
    env = os.environ.copy()
    mdrun = [
        str(Path(gmx).resolve()), "mdrun", "-deffnm", stage,
        "-ntmpi", "1", "-ntomp", str(threads), "-pin", "on",
    ]
    gpu_support = (gromacs_version(gmx).get("gpu_support") or "").lower()
    use_gpu = accelerator == "gpu" or (accelerator == "auto" and gpu_support not in {"", "disabled", "none"})
    if use_gpu:
        mdrun += ["-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "cpu"]
    checkpoint = seed_dir / f"{stage}.cpt"
    resumed = bool(resume and checkpoint.is_file())
    if resumed:
        mdrun += ["-cpi", checkpoint.name, "-append"]
    result = subprocess.run(mdrun, cwd=seed_dir, capture_output=True, text=True, env=env)
    if result.returncode:
        raise RuntimeError(result.stdout[-4000:] + result.stderr[-4000:])
    payload = {"stage": stage, "status": "completed", "threads": threads, "accelerator": "gpu" if use_gpu else "cpu", "resumed": resumed, "command": mdrun, "log": str(seed_dir / f"{stage}.log"), "coordinates": str(seed_dir / f"{stage}.gro"), "checkpoint": str(seed_dir / f"{stage}.cpt"), "energy": str(seed_dir / f"{stage}.edr"), "trajectory": str(seed_dir / f"{stage}.xtc")}
    path = seed_dir / f"{stage}-execution.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def _read_xvg(path: Path) -> np.ndarray:
    rows = [line.split() for line in Path(path).read_text().splitlines() if line and line[0] not in "#@"]
    return np.asarray(rows, dtype=float)


def analyze_production_trajectory(
    seed_dir: Path,
    gmx: Path,
    *,
    sample_interval_ps: float = 5.0,
    config: GromacsMethanolConfig | None = None,
) -> dict:
    """Write thermodynamic, H-bond, block, autocorrelation, and candidate records."""
    import MDAnalysis as mda

    resolved = config or GromacsMethanolConfig()
    seed_dir = Path(seed_dir)
    xvg = seed_dir / "production-observables.xvg"
    _run(
        [str(Path(gmx).resolve()), "energy", "-f", "production.edr", "-o", xvg.name],
        seed_dir, stdin="Potential\nTemperature\n0\n",
    )
    energy = _read_xvg(xvg)
    universe = mda.Universe(str(seed_dir / "nvt.gro"), str(seed_dir / "production.xtc"))
    stride = max(1, round(sample_interval_ps / float(universe.trajectory.dt)))
    symbols = ["H" if atom.name.startswith("H") else atom.name for atom in universe.atoms]
    seed_match = re.search(r"seed-(\d+)$", seed_dir.name)
    seed = int(seed_match.group(1)) if seed_match else 0
    records, candidates = [], []
    total_mass_kg = resolved.molecule_count * 32.04186 / 1000.0 / 6.02214076e23
    shell_config = BulkEmbeddingConfig()
    for sampled_index, ts in enumerate(universe.trajectory[::stride]):
        coordinates = np.asarray(universe.atoms.positions, dtype=float).copy()
        box = float(ts.dimensions[0])
        oxygens = coordinates[1::6]
        hydroxyl_h = coordinates[2::6]
        donor_vectors = minimum_image_array(oxygens - hydroxyl_h, box)
        h_to_o = minimum_image_array(oxygens[None, :, :] - hydroxyl_h[:, None, :], box)
        distances = np.linalg.norm(h_to_o, axis=2)
        left = donor_vectors[:, None, :]
        denom = np.linalg.norm(left, axis=2) * np.maximum(distances, 1.0e-12)
        cosines = np.clip(np.sum(left * h_to_o, axis=2) / denom, -1.0, 1.0)
        angles = np.degrees(np.arccos(cosines))
        hbonds = (distances <= 2.6) & (angles >= 140.0)
        np.fill_diagonal(hbonds, False)
        count = int(hbonds.sum())
        volume_m3 = float(np.prod(ts.dimensions[:3])) * 1.0e-30
        time_ps = float(ts.time)
        energy_row = energy[int(np.argmin(np.abs(energy[:, 0] - time_ps)))]
        records.append({
            "time_ps": time_ps,
            "temperature_K": float(energy_row[2]),
            "density_g_cm3": total_mass_kg / volume_m3 / 1000.0,
            "potential_energy": float(energy_row[1]),
            "h_bonds_per_molecule": count / resolved.molecule_count,
        })
        block = min(4, int(time_ps // max(resolved.production_ps / 5.0, 1.0)))
        if sampled_index % 5 == 0:
            central = (sampled_index * 37 + seed * 17) % resolved.molecule_count
            shell = extract_central_solvation_shell(
                symbols, coordinates, central_molecule=central,
                box_length_angstrom=box, config=shell_config,
            )
            if shell["hydrogen_bonds"]:
                candidates.append({
                    "snapshot_id": f"seed-{seed}-t-{time_ps:07.1f}-mol-{central:03d}",
                    "trajectory_seed": seed, "trajectory_block": block,
                    "time_ps": time_ps, "central_molecule": central,
                    "box_length_angstrom": box, "symbols": symbols,
                    "coordinates": coordinates.tolist(),
                    "donor_count": sum(x["donor_molecule"] == central for x in shell["hydrogen_bonds"]),
                    "acceptor_count": sum(x["acceptor_molecule"] == central for x in shell["hydrogen_bonds"]),
                    "hydrogen_bonds": shell["hydrogen_bonds"],
                })
    blocks = []
    for block in range(5):
        rows = records[block * len(records) // 5:(block + 1) * len(records) // 5]
        blocks.append({"block": block, **{
            key: float(np.mean([row[key] for row in rows]))
            for key in ("temperature_K", "density_g_cm3", "potential_energy", "h_bonds_per_molecule")
        }})
    statistics = {
        key: {"mean": float(np.mean(values)), "stddev": float(np.std(values)), **statistical_inefficiency(values)}
        for key in ("temperature_K", "density_g_cm3", "potential_energy", "h_bonds_per_molecule")
        for values in [[row[key] for row in records]]
    }
    potential = np.asarray([row["potential_energy"] for row in records])
    energy_slope = float(np.polyfit(np.arange(len(potential)), potential, 1)[0])
    energy_stationary = abs(energy_slope) <= max(float(potential.std()), 1.0) / len(potential)
    statistics["potential_energy"]["linear_drift_per_sample"] = energy_slope
    report = {
        "schema_version": 1, "kind": "gromacs_bulk_methanol_production_analysis",
        "seed": seed, "frames_analyzed": len(records), "sample_interval_ps": sample_interval_ps,
        "statistics": statistics, "blocks": blocks, "environment_candidates": candidates,
        "gates": {
            "temperature": abs(statistics["temperature_K"]["mean"] - resolved.temperature_K) / resolved.temperature_K <= 0.03,
            "density": abs(statistics["density_g_cm3"]["mean"] - 0.7866) / 0.7866 <= 0.03,
            "energy_stationarity": energy_stationary,
            "hydrogen_bonds": statistics["h_bonds_per_molecule"]["mean"] > 0.5,
        },
    }
    report["status"] = "validated" if all(report["gates"].values()) else "failed_closed"
    path = seed_dir / "production-analysis.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["artifact"] = str(path)
    return report


def minimum_image_array(displacement: np.ndarray, box_length: float) -> np.ndarray:
    return displacement - box_length * np.round(displacement / box_length)


def select_embedded_environments(
    reports: list[dict], *, acquisition_count: int = 12, heldout_count: int = 4,
) -> dict:
    """Select seed/block/coordination-balanced environments without spectral data."""
    candidates = [item for report in reports for item in report["environment_candidates"]]
    class_counts = {}
    for item in candidates:
        key = f"D{item['donor_count']}A{item['acceptor_count']}"
        class_counts[key] = class_counts.get(key, 0) + 1
    acquisition_pool = [item for item in candidates if int(item["trajectory_block"]) < 4]
    heldout_pool = [item for item in candidates if int(item["trajectory_block"]) == 4]

    def choose(pool: list[dict], count: int) -> list[dict]:
        selected = []
        groups = {}
        for item in pool:
            key = (int(item["trajectory_seed"]), int(item["donor_count"]), int(item["acceptor_count"]))
            groups.setdefault(key, []).append(item)
        ordered = sorted(groups, key=lambda key: (len(groups[key]), key))
        while len(selected) < count and any(groups.values()):
            for key in ordered:
                rows = groups[key]
                if not rows or len(selected) >= count:
                    continue
                target = rows[len(rows) // 2]
                selected.append(target)
                groups[key] = [row for row in rows if abs(float(row["time_ps"]) - float(target["time_ps"])) >= 100.0]
        return selected

    def annotate(rows: list[dict], role: str) -> list[dict]:
        result = []
        for item in rows:
            key = f"D{item['donor_count']}A{item['acceptor_count']}"
            result.append({
                **item, "selection_role": role, "coordination_class": key,
                "trajectory_cluster_id": f"seed-{item['trajectory_seed']}-block-{item['trajectory_block']}",
                "trajectory_occupancy_weight": class_counts[key] / len(candidates),
            })
        return result

    acquisition = annotate(choose(acquisition_pool, acquisition_count), "acquisition")
    heldout = annotate(choose(heldout_pool, heldout_count), "heldout")
    return {
        "schema_version": 1, "kind": "bulk_embedded_environment_selection",
        "selection_features": ["trajectory_seed", "trajectory_block", "donor_count", "acceptor_count", "time_separation"],
        "spectral_information_used": False, "candidate_count": len(candidates),
        "coordination_class_occupancy": {key: value / len(candidates) for key, value in class_counts.items()},
        "acquisition": acquisition, "heldout": heldout,
        "status": "completed" if len(acquisition) == acquisition_count and len(heldout) == heldout_count else "failed_closed",
    }
