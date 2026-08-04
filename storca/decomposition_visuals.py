"""Portable visual assets for modelled and verified decomposition routes."""

from __future__ import annotations

import json
from pathlib import Path
import re

import numpy as np


def render_candidate_storyboard(parent_smiles: str, product_smiles: list[str], output: Path, *, title: str) -> Path:
    """Render a labelled 2-D candidate-route storyboard, not a trajectory claim."""
    from rdkit import Chem
    from rdkit.Chem import Draw

    molecules = [Chem.MolFromSmiles(parent_smiles), *(Chem.MolFromSmiles(item) for item in product_smiles)]
    if any(molecule is None for molecule in molecules):
        raise ValueError("Could not render one or more route SMILES")
    labels = ["parent", *[f"product {index}" for index in range(1, len(product_smiles) + 1)]]
    image = Draw.MolsToGridImage(molecules, molsPerRow=len(molecules), subImgSize=(320, 260), legends=labels)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(str(output))
    output.with_suffix(".txt").write_text(
        f"{title}\nComputed candidate-route storyboard; this is not an IRC or molecular-dynamics trajectory.\n"
    )
    return output


def xyz_frames(path: Path) -> list[tuple[list[str], np.ndarray, str]]:
    """Read a multi-frame XYZ file while retaining per-frame comments."""
    lines = Path(path).read_text().splitlines()
    frames, index = [], 0
    while index < len(lines):
        if not lines[index].strip():
            index += 1
            continue
        count = int(lines[index])
        comment = lines[index + 1] if index + 1 < len(lines) else ""
        block = [row.split() for row in lines[index + 2:index + 2 + count]]
        if len(block) != count or any(len(row) < 4 for row in block):
            raise ValueError("Incomplete XYZ trajectory frame")
        frames.append((
            [row[0] for row in block],
            np.asarray([[float(value) for value in row[1:4]] for row in block]),
            comment,
        ))
        index += count + 2
    if not frames:
        raise ValueError("XYZ trajectory has no frames")
    return frames


def _aligned_frames(frames: list[tuple[list[str], np.ndarray, str]]) -> list[tuple[list[str], np.ndarray, str]]:
    reference = frames[0][1] - frames[0][1].mean(axis=0)
    aligned = []
    for elements, coordinates, comment in frames:
        centered = coordinates - coordinates.mean(axis=0)
        covariance = centered.T @ reference
        left, _, right = np.linalg.svd(covariance)
        correction = np.eye(3)
        correction[-1, -1] = np.sign(np.linalg.det(left @ right))
        aligned.append((elements, centered @ left @ correction @ right, comment))
    return aligned


def _comment_energy(comment: str) -> float | None:
    matches = re.findall(r"(?:energy|E)\s*[=:]\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)", comment, re.IGNORECASE)
    return float(matches[-1]) if matches else None


def _bonds(elements: list[str], coordinates: np.ndarray) -> list[tuple[int, int, float]]:
    radii = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
             "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39}
    bonds = []
    for left in range(len(elements)):
        for right in range(left + 1, len(elements)):
            distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
            cutoff = 1.25 * (radii.get(elements[left], 0.8) + radii.get(elements[right], 0.8))
            if 0.2 < distance <= cutoff:
                bonds.append((left, right, distance))
    return bonds


def render_energy_profile(energies: list[float], output: Path, *,
                          title: str = "Computed reaction path energy") -> Path:
    """Write a compact relative-energy plot for a computed path."""
    import matplotlib.pyplot as plt

    if len(energies) < 2:
        raise ValueError("An energy profile needs at least two points")
    relative = np.asarray(energies, dtype=float) - float(energies[0])
    figure, axis = plt.subplots(figsize=(6, 3.5))
    axis.plot(range(1, len(relative) + 1), relative, marker="o", linewidth=1.5, markersize=3)
    axis.set(xlabel="Path frame", ylabel="Relative energy", title=title)
    axis.grid(alpha=0.25)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(figure)
    return output


def render_xyz_trajectory_gif(
    trajectory_xyz: Path,
    output: Path,
    *,
    title: str = "Computed IRC trajectory",
    evidence_label: str = "ORCA-computed path",
    energies: list[float] | None = None,
) -> Path:
    """Render an aligned XYZ trajectory with dynamic connectivity as a GIF."""
    import io
    import matplotlib.pyplot as plt
    from PIL import Image

    frames = _aligned_frames(xyz_frames(trajectory_xyz))
    expected_elements = frames[0][0]
    if any(elements != expected_elements for elements, _, _ in frames):
        raise ValueError("All trajectory frames must preserve atom order")
    coordinates = np.concatenate([frame[1] for frame in frames], axis=0)
    limits = [
        (float(coordinates[:, axis].min() - 0.7), float(coordinates[:, axis].max() + 0.7))
        for axis in range(3)
    ]
    colours = {"H": "#eeeeee", "C": "#333333", "N": "#3366cc", "O": "#cc3333",
               "F": "#55aa55", "P": "#ef8f22", "S": "#d6bb22", "Cl": "#44aa44"}
    parsed_energies = [_comment_energy(comment) for _, _, comment in frames]
    plot_energies = energies if energies is not None else (
        [float(value) for value in parsed_energies] if all(value is not None for value in parsed_energies) else None
    )
    if plot_energies is not None and len(plot_energies) != len(frames):
        raise ValueError("Energy and trajectory frame counts differ")
    relative = (
        np.asarray(plot_energies) - float(plot_energies[0])
        if plot_energies is not None else None
    )
    images = []
    for frame_index, (elements, points, _) in enumerate(frames):
        figure = plt.figure(figsize=(7.2 if relative is not None else 5.4, 5.2))
        axis = figure.add_subplot(121 if relative is not None else 111, projection="3d")
        for left, right, _ in _bonds(elements, points):
            axis.plot(*zip(points[left], points[right]), color="#777777", linewidth=3, alpha=0.75)
        for element, (x, y, z) in zip(elements, points):
            axis.scatter(x, y, z, s=180, color=colours.get(element, "#aa55aa"), edgecolors="black")
        axis.set(xlim=limits[0], ylim=limits[1], zlim=limits[2], title=title)
        axis.text2D(0.02, 0.02, f"{evidence_label}\nframe {frame_index + 1}/{len(frames)}",
                    transform=axis.transAxes, fontsize=8)
        axis.set_axis_off()
        if relative is not None:
            energy_axis = figure.add_subplot(122)
            energy_axis.plot(range(len(relative)), relative, color="#777777", linewidth=1.4)
            energy_axis.scatter([frame_index], [relative[frame_index]], color="#cc3333", s=35, zorder=3)
            energy_axis.set(xlabel="Path frame", ylabel="Relative energy", title="Energy profile")
            energy_axis.grid(alpha=0.2)
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
        plt.close(figure)
        buffer.seek(0)
        images.append(Image.open(buffer).convert("P"))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    images[0].save(output, save_all=True, append_images=images[1:], duration=160, loop=0)
    return output


def write_visual_manifest(output: Path, *, route_id: str, evidence_level: str,
                          storyboard: Path | None = None, animation: Path | None = None,
                          energy_profile: Path | None = None, trajectory: Path | None = None) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({
        "schema_version": 1,
        "kind": "decomposition_visuals",
        "route_id": route_id,
        "evidence_level": evidence_level,
        "artifacts": {
            "candidate_storyboard": str(storyboard) if storyboard else None,
            "animation": str(animation) if animation else None,
            "energy_profile": str(energy_profile) if energy_profile else None,
            "trajectory": str(trajectory) if trajectory else None,
        },
        "interpretation": (
            "candidate_storyboard is topology-only; computed_dissociation_scan is a ground-state coordinate path; "
            "verified_transition_state_path requires one imaginary mode and an IRC matching both endpoints."
        ),
    }, indent=2, sort_keys=True) + "\n")
    return output
