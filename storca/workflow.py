"""Non-interactive ORCA optimization and frequency workflow."""

from __future__ import annotations

from pathlib import Path

from src.inputgen import create_orca_input
from src.orca_runner import run_orca
from src.stability.freq_check import frequency_stability_check

from .runs import write_metadata


def run_optimization_and_frequency(
    xyz_file: Path,
    run_dir: Path,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    ncores: int = 1,
    run_frequency: bool = True,
    method_keywords: list[str] | None = None,
    timeout_seconds: float | None = None,
) -> dict:
    """Run ORCA in *run_dir* and return paths and parsed frequency status.

    The frequency job always consumes ORCA's optimized XYZ geometry.  A failed
    optimization or a missing optimized geometry stops the workflow early.
    """
    xyz_file = Path(xyz_file)
    run_dir = Path(run_dir)
    if xyz_file.parent != run_dir:
        raise ValueError("Input XYZ must be located inside the run directory")

    opt_input = create_orca_input(
        xyz_file,
        charge=charge,
        multiplicity=multiplicity,
        opt=True,
        label="opt",
        ncores=ncores,
        method_keywords=method_keywords,
    )
    opt_outputs = run_orca(opt_input, timeout_seconds=timeout_seconds)
    optimized_xyz = opt_outputs["xyz"]
    if not optimized_xyz.exists():
        raise RuntimeError(f"ORCA completed without an optimized geometry: {optimized_xyz}")

    result: dict = {
        "run_dir": run_dir,
        "input_xyz": xyz_file,
        "optimized_xyz": optimized_xyz,
        "optimization": opt_outputs,
    }
    if run_frequency:
        freq_input = create_orca_input(
            optimized_xyz,
            charge=charge,
            multiplicity=multiplicity,
            freq=True,
            label="freq",
            ncores=ncores,
            method_keywords=method_keywords,
        )
        freq_outputs = run_orca(freq_input, timeout_seconds=timeout_seconds)
        result["frequency"] = freq_outputs
        result["frequency_check"] = frequency_stability_check(freq_outputs["out"])

    write_metadata(
        run_dir,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        run_frequency=run_frequency,
        orca_method_keywords=method_keywords,
        timeout_seconds=timeout_seconds,
        optimized_xyz=str(optimized_xyz),
    )
    return result
