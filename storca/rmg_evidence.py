"""Read retained RMG artifacts into explicit evidence records."""

from __future__ import annotations

import csv
import re
from pathlib import Path


def normalize_rmg_label(label: str) -> str:
    """Remove RMG's trailing numeric species index from a label."""
    return re.sub(r"\(\d+\)$", "", label).strip()


def parse_species_dictionary(path: Path) -> dict[str, dict]:
    """Parse RMG adjacency-list records by exact and compatibility labels.

    Exact Chemkin labels retain their RMG index.  A normalized alias is kept
    only when unambiguous enough for legacy callers; it never overwrites a
    prior species record.
    """
    records: dict[str, dict] = {}
    blocks = re.split(r"\n\s*\n", Path(path).read_text().strip())
    for block in blocks:
        lines = [line.rstrip() for line in block.splitlines() if line.strip()]
        if not lines:
            continue
        label = lines[0].strip()
        body = "\n".join(lines[1:])
        multiplicity_match = re.search(r"^multiplicity\s+(\d+)", body, re.MULTILINE)
        radical_electrons = sum(int(value) for value in re.findall(r"\bu(\d+)", body))
        normalized = normalize_rmg_label(label)
        record = {
            "label": label,
            "normalized_label": normalized,
            "adjacency_list": body,
            "multiplicity": int(multiplicity_match.group(1)) if multiplicity_match else 1,
            "radical_electrons": radical_electrons,
            "is_radical": radical_electrons > 0,
        }
        records[label] = record
        records.setdefault(normalized, record)
    return records


def parse_annotated_reaction_provenance(path: Path) -> dict[str, dict]:
    """Associate annotated Chemkin comments with each reaction equation."""
    records: dict[str, dict] = {}
    comments: list[str] = []
    reaction_pattern = re.compile(r"^(?P<equation>.+?(?:<=>|=>).+?)\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s*$")
    for raw_line in Path(path).read_text().splitlines():
        line = raw_line.strip()
        if line.startswith("!"):
            comments.append(line[1:].strip())
            continue
        match = reaction_pattern.match(line)
        if not match:
            if line and line != "REACTIONS":
                comments = []
            continue
        equation = match.group("equation")
        source_text = "\n".join(comments)
        if "Matched reaction" in source_text or "From Training reaction" in source_text:
            rate_source = "training"
        elif "Estimated using" in source_text or "Estimated from node" in source_text:
            rate_source = "rate_rule"
        else:
            rate_source = "unknown"
        template = re.search(r"Template reaction:\s*(.+)", source_text)
        family = re.search(r"family:\s*(.+)", source_text)
        records[equation] = {
            "reaction_equation": equation,
            "reaction_family": (family.group(1).strip() if family else (template.group(1).strip() if template else None)),
            "rate_source": rate_source,
            "provenance_comments": comments,
        }
        comments = []
    return records


def _final_profile_path(solver_dir: Path) -> Path | None:
    paths = list(Path(solver_dir).glob("simulation_*.csv"))
    if not paths:
        return None
    def order(path: Path) -> int:
        values = re.findall(r"(\d+)", path.stem)
        return int(values[-1]) if values else -1
    return max(paths, key=order)


def parse_final_solver_profile(solver_dir: Path, target_label: str) -> dict | None:
    """Return final-profile mole-fraction extrema and target loss.

    RMG simple-reactor CSVs report mole fractions (alongside volume), not
    amounts.  The field names state that basis explicitly so downstream code
    cannot accidentally treat a presence floor as an absolute mole amount.
    """
    path = _final_profile_path(solver_dir)
    if path is None:
        return None
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    resolved_target_label = next((name for name in rows[0] if normalize_rmg_label(name) == target_label), None)
    if resolved_target_label is None:
        return None
    concentrations = {
        name: {"initial_mole_fraction": float(rows[0][name]), "final_mole_fraction": float(rows[-1][name]),
               "maximum_mole_fraction": max(float(row[name]) for row in rows)}
        for name in rows[0] if name not in {"Time (s)", "Volume (m^3)"}
    }
    target = concentrations[resolved_target_label]
    initial = target["initial_mole_fraction"]
    target_loss = 0.0 if initial == 0 else max(0.0, 1.0 - target["final_mole_fraction"] / initial)
    target_time_series = [
        {"time_seconds": float(row["Time (s)"]), "fraction_remaining": float(row[resolved_target_label]) / initial}
        for row in rows if initial > 0
    ]
    return {
        "path": str(path),
        "end_time_seconds": float(rows[-1]["Time (s)"]),
        "species": concentrations,
        "target_loss_fraction": target_loss,
        "target_fraction_remaining": 1.0 - target_loss,
        "target_profile_label": resolved_target_label,
        "profile_basis": "mole_fraction",
        "target_time_series": target_time_series,
    }


def time_to_retention_seconds(profile: dict | None, retention_fraction: float) -> float | None:
    """Linearly interpolate the first modeled crossing of a retention target."""
    if not profile:
        return None
    series = profile.get("target_time_series", [])
    for previous, current in zip(series, series[1:]):
        if previous["fraction_remaining"] >= retention_fraction >= current["fraction_remaining"]:
            change = previous["fraction_remaining"] - current["fraction_remaining"]
            if change == 0:
                return current["time_seconds"]
            fraction = (previous["fraction_remaining"] - retention_fraction) / change
            return previous["time_seconds"] + fraction * (current["time_seconds"] - previous["time_seconds"])
    return None
