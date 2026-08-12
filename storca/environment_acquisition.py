"""Mode-class-aware acquisition of environment representatives under a hard budget."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np


ACQUISITION_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class EnvironmentAcquisitionConfig:
    schema_version: int = ACQUISITION_SCHEMA_VERSION
    minimum_sampled_environments: int = 3
    minimum_sampled_fraction: float = 0.05
    target_dft_environments_per_class: int = 3
    geometry_diversity_weight: float = 0.15
    frequency_diversity_weight: float = 0.35
    topology_diversity_bonus: float = 3.00


def local_mode_class(mode: dict) -> str | None:
    """Return a stable chemical/coordination identity for a localized mode."""
    if not mode.get("bond_class"):
        return None
    return ":".join((
        str(mode.get("bond_class")),
        str(mode.get("spectral_band_class", "xh_stretch")),
        str(mode.get("coordination_class", "unclassified")),
    ))


def build_candidate_mode_profiles(frequency_records: Iterable[dict]) -> dict[str, dict]:
    """Collapse multiple local oscillators into one contribution per candidate/class."""
    profiles: dict[str, dict] = {}
    for record in frequency_records:
        if record.get("frequency_status") != "completed":
            continue
        grouped: dict[str, list[float]] = {}
        for mode in record.get("modes") or []:
            mode_class = local_mode_class(mode)
            frequency = mode.get("freq")
            if mode_class and isinstance(frequency, (int, float)) and math.isfinite(frequency):
                grouped.setdefault(mode_class, []).append(float(frequency))
        profiles[str(record["candidate_id"])] = {
            "candidate_id": str(record["candidate_id"]),
            "mode_classes": sorted(grouped),
            "class_frequencies_cm-1": {
                key: float(np.mean(values)) for key, values in grouped.items()
            },
            "snapshot_reliability": record.get("snapshot_reliability"),
        }
    return profiles


def significant_mode_classes(
    profiles: dict[str, dict], *, config: EnvironmentAcquisitionConfig | None = None,
) -> dict[str, dict]:
    """Identify prevalent sampled classes without promoting isolated outliers."""
    resolved = config or EnvironmentAcquisitionConfig()
    total = len(profiles)
    support: dict[str, list[str]] = {}
    frequencies: dict[str, list[float]] = {}
    for candidate_id, profile in profiles.items():
        for mode_class in profile["mode_classes"]:
            support.setdefault(mode_class, []).append(candidate_id)
            frequencies.setdefault(mode_class, []).append(
                profile["class_frequencies_cm-1"][mode_class]
            )
    output = {}
    for mode_class in sorted(support):
        count = len(set(support[mode_class]))
        fraction = count / total if total else 0.0
        significant = bool(
            count >= resolved.minimum_sampled_environments
            and fraction >= resolved.minimum_sampled_fraction
        )
        output[mode_class] = {
            "sampled_environments": count,
            "sampled_fraction": fraction,
            "frequency_min_cm-1": min(frequencies[mode_class]),
            "frequency_max_cm-1": max(frequencies[mode_class]),
            "frequency_span_cm-1": max(frequencies[mode_class]) - min(frequencies[mode_class]),
            "significant": significant,
            "significance_reason": (
                "prevalent_sampled_mode_class" if significant
                else "below_sampled_support_threshold"
            ),
        }
    return output


def select_mode_coverage_representatives(
    records: list[dict], frequency_records: Iterable[dict], distance_matrix: np.ndarray,
    representative_count: int, *, already_selected_ids: Iterable[str] = (),
    config: EnvironmentAcquisitionConfig | None = None,
) -> tuple[list[int], dict]:
    """Greedily close class deficits, then maximize frequency/geometry diversity."""
    resolved = config or EnvironmentAcquisitionConfig()
    if not 0 <= representative_count <= len(records):
        raise ValueError("Invalid mode-aware representative count")
    matrix = np.asarray(distance_matrix, dtype=float)
    if matrix.shape != (len(records), len(records)) or not np.all(np.isfinite(matrix)):
        raise ValueError("Mode-aware acquisition requires a finite square distance matrix")
    profiles = build_candidate_mode_profiles(frequency_records)
    classes = significant_mode_classes(profiles, config=resolved)
    significant = {key for key, value in classes.items() if value["significant"]}
    index_by_id = {str(record["candidate_id"]): index for index, record in enumerate(records)}
    selected = [
        index_by_id[candidate_id] for candidate_id in already_selected_ids
        if candidate_id in index_by_id
    ]
    selected = list(dict.fromkeys(selected))

    def coverage(indices: list[int]) -> dict[str, int]:
        return {
            mode_class: sum(
                mode_class in profiles.get(str(records[index]["candidate_id"]), {}).get("mode_classes", [])
                for index in indices
            )
            for mode_class in significant
        }

    decisions = []
    while len(selected) < representative_count:
        counts = coverage(selected)
        remaining = [index for index in range(len(records)) if index not in selected]
        if not remaining:
            break

        def score(index: int) -> tuple[float, float, str]:
            candidate_id = str(records[index]["candidate_id"])
            profile = profiles.get(candidate_id, {"mode_classes": [], "class_frequencies_cm-1": {}})
            deficit_score = sum(
                max(0, resolved.target_dft_environments_per_class - counts.get(mode_class, 0))
                for mode_class in profile["mode_classes"] if mode_class in significant
            )
            frequency_score = 0.0
            for mode_class in profile["mode_classes"]:
                if mode_class not in significant:
                    continue
                selected_frequencies = [
                    profiles.get(str(records[chosen]["candidate_id"]), {}).get(
                        "class_frequencies_cm-1", {}
                    ).get(mode_class)
                    for chosen in selected
                ]
                selected_frequencies = [value for value in selected_frequencies if value is not None]
                if selected_frequencies:
                    frequency_score += min(
                        abs(profile["class_frequencies_cm-1"][mode_class] - value)
                        for value in selected_frequencies
                    ) / 50.0
            geometry_score = min((float(matrix[index, chosen]) for chosen in selected), default=0.0)
            selected_topologies = {records[chosen].get("topology", "dimer") for chosen in selected}
            topology_bonus = (
                resolved.topology_diversity_bonus
                if records[index].get("topology", "dimer") not in selected_topologies else 0.0
            )
            total = (
                deficit_score
                + resolved.frequency_diversity_weight * frequency_score
                + resolved.geometry_diversity_weight * geometry_score
                + topology_bonus
            )
            return total, geometry_score, candidate_id

        chosen = max(remaining, key=score)
        before = coverage(selected)
        selected.append(chosen)
        after = coverage(selected)
        decisions.append({
            "acquisition_position": len(selected),
            "candidate_id": records[chosen]["candidate_id"],
            "topology": records[chosen].get("topology", "dimer"),
            "association_class": _association_class(records[chosen]),
            "mode_classes": profiles.get(str(records[chosen]["candidate_id"]), {}).get("mode_classes", []),
            "coverage_before": before,
            "coverage_after": after,
            "selection_score": score(chosen)[0],
            "reason": (
                "closes_significant_mode_class_deficit"
                if any(after.get(key, 0) > before.get(key, 0) and before.get(key, 0) < resolved.target_dft_environments_per_class for key in significant)
                else "remaining_geometry_frequency_diversity"
            ),
        })
    final_coverage = coverage(selected)
    deficits = {
        mode_class: max(
            0, resolved.target_dft_environments_per_class - final_coverage.get(mode_class, 0)
        )
        for mode_class in sorted(significant)
    }
    return selected, {
        "schema_version": ACQUISITION_SCHEMA_VERSION,
        "kind": "mode_class_coverage_constrained_orca_acquisition",
        "configuration": asdict(resolved),
        "representative_budget": representative_count,
        "candidate_profiles": profiles,
        "mode_classes": classes,
        "significant_mode_classes": sorted(significant),
        "decisions": decisions,
        "selected_candidate_ids": [records[index]["candidate_id"] for index in selected],
        "final_independent_dft_coverage": final_coverage,
        "remaining_class_deficits": deficits,
        "status": (
            "coverage_targets_satisfied" if not any(deficits.values())
            else "orca_budget_insufficient_for_mode_class_coverage"
        ),
        "stop_reason": (
            "all_significant_classes_have_target_coverage" if not any(deficits.values())
            else "hard_environment_orca_budget_reached"
        ),
    }


def _association_class(record: dict) -> str:
    distance = float(record["environment_features"]["h_bond_distance_angstrom"])
    return "strong" if distance < 1.85 else "intermediate" if distance <= 2.20 else "weak"


def write_environment_acquisition_report(run_dir: Path, report: dict) -> Path:
    """Persist acquisition state and merge its summary into related diagnostics."""
    run_dir = Path(run_dir)
    path = run_dir / "environment-acquisition.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    for name in ("environment-sampling.json", "environment-convergence.json"):
        target = run_dir / name
        if not target.is_file():
            continue
        payload = json.loads(target.read_text())
        payload["environment_acquisition"] = {
            "artifact": str(path),
            "status": report["status"],
            "stop_reason": report["stop_reason"],
            "representative_budget": report["representative_budget"],
            "remaining_class_deficits": report["remaining_class_deficits"],
        }
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def record_acquisition_batch(run_dir: Path, endpoint: int, transfer_artifacts: dict) -> dict:
    """Record transfer feedback after one cumulative ORCA acquisition batch."""
    path = Path(run_dir) / "environment-acquisition.json"
    if not path.is_file():
        return {"status": "environment_acquisition_not_configured"}
    report = json.loads(path.read_text())
    batches = report.setdefault("orca_batches", [])
    entry = {
        "completed_representatives": int(endpoint),
        "transfer_status": transfer_artifacts.get("status"),
        "frequency_transfer_validation": str(
            transfer_artifacts.get("frequency_transfer_validation") or ""
        ),
    }
    validation_path = transfer_artifacts.get("frequency_transfer_validation")
    if validation_path and Path(validation_path).is_file():
        validation = json.loads(Path(validation_path).read_text())
        entry["limiting_mode_classes"] = sorted(
            key for key, value in validation.get("classes", {}).items()
            if value.get("status") != "validated"
        )
    batches = [item for item in batches if int(item["completed_representatives"]) != int(endpoint)]
    batches.append(entry)
    report["orca_batches"] = sorted(batches, key=lambda item: item["completed_representatives"])
    if transfer_artifacts.get("status") == "validated":
        report["status"] = "dft_transfer_validated"
        report["stop_reason"] = "all_significant_classes_passed_withheld_validation"
    elif endpoint >= int(report.get("representative_budget", 0)):
        report["status"] = "orca_budget_reached_before_transfer_validation"
        report["stop_reason"] = "hard_environment_orca_budget_reached"
    else:
        report["status"] = "additional_orca_representative_required"
        report["stop_reason"] = "dft_transfer_validation_incomplete"
    write_environment_acquisition_report(run_dir, report)
    return report


def reprioritize_pending_representatives(run_dir: Path, completed_count: int) -> list[Path]:
    """Move pending representatives covering failed classes to the next batch."""
    run_dir = Path(run_dir)
    acquisition_path = run_dir / "environment-acquisition.json"
    manifest_path = run_dir / "clusters" / "selected-conformers.json"
    validation_path = run_dir / "frequency-transfer-validation.json"
    if not all(path.is_file() for path in (acquisition_path, manifest_path, validation_path)):
        return []
    acquisition = json.loads(acquisition_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    validation = json.loads(validation_path.read_text())
    entries = manifest.get("conformers") or []
    if not 0 <= completed_count <= len(entries):
        raise ValueError("Completed acquisition count is outside the representative manifest")
    limiting = {
        key for key, value in validation.get("classes", {}).items()
        if value.get("status") != "validated"
    }
    profiles = acquisition.get("candidate_profiles") or {}
    completed, pending = entries[:completed_count], entries[completed_count:]
    original = [entry["source_xtb_candidate_id"] for entry in pending]
    pending.sort(key=lambda entry: (
        -sum(
            mode_class in limiting
            for mode_class in profiles.get(
                str(entry["source_xtb_candidate_id"]), {}
            ).get("mode_classes", [])
        ),
        -float((entry.get("environment_features") or {}).get(
            "estimated_local_frequency_cm-1"
        ) or 0.0),
        str(entry["source_xtb_candidate_id"]),
    ))
    reordered = completed + pending
    for position, entry in enumerate(reordered, start=1):
        entry["selected_position"] = position
        entry["execution_position"] = position
    manifest["conformers"] = reordered
    manifest["execution_order_model"] = "adaptive_failed_mode_classes_then_frequency_geometry_diversity"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    updated = [entry["source_xtb_candidate_id"] for entry in pending]
    if updated != original:
        acquisition.setdefault("reprioritization_events", []).append({
            "after_completed_representatives": completed_count,
            "limiting_mode_classes": sorted(limiting),
            "previous_pending_candidate_ids": original,
            "updated_pending_candidate_ids": updated,
        })
        write_environment_acquisition_report(run_dir, acquisition)
    return [Path(entry["xyz"]) for entry in reordered]
