"""Admission and immutable storage of backend-proven reaction correspondence.

No graph matching happens here. The owning RMG environment establishes the
correspondence; these checks bind its evidence to the source library, exported
network and inspection used by the requesting consumer. Native IDs stay intact.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import re


IDENTITY_POLICY = "rmg-resonance-strict-v1"
SOURCE_FILES = ("reactions.py", "dictionary.txt", "generated-kinetics.json")


def _encoded(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _digest(value: object) -> str:
    return hashlib.sha256(_encoded(value)).hexdigest()


def _artifact(path: Path) -> dict:
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _check_artifact(record: dict, expected_path: Path | None = None) -> None:
    path = Path(record["path"]).resolve(strict=True)
    if expected_path is not None and path != Path(expected_path).resolve(strict=True):
        raise ValueError("identity_artifact_path_mismatch")
    if _artifact(path) != record:
        raise ValueError("identity_artifact_missing_or_mutated")


def build_replacement_identity_requests(
    libraries: list[Path], *, temperature_K: float, pressure_bar: float,
    condition_contract: dict | None = None,
) -> list[dict]:
    """Prepare exact, validated repair libraries for the owning graph matcher."""
    from .generated_kinetics import validate_generated_library

    requests = []
    seen = set()
    for library in dict.fromkeys(Path(path).resolve() for path in libraries):
        manifest = validate_generated_library(
            library, temperature_K=temperature_K, pressure_bar=pressure_bar,
            condition_contract=condition_contract,
        )
        if manifest.get("kind") != "storca_generated_kinetics" or not manifest.get("route_id"):
            raise ValueError("Replacement identity requires a generated, route-bound library")
        if manifest["route_id"] in seen:
            raise ValueError("Replacement identity has duplicate source route IDs")
        seen.add(manifest["route_id"])
        artifacts = {name: _artifact(library / name) for name in SOURCE_FILES}
        for name in SOURCE_FILES[:2]:
            if artifacts[name]["sha256"] != (manifest.get("checksums") or {}).get(name):
                raise ValueError("Replacement identity requires exact library checksums")
        requests.append({"route_id": manifest["route_id"], "library": str(library),
                         "source_artifacts": artifacts})
    return requests


def retain_replacement_identity(inspection: dict, directory: Path) -> dict:
    """Append a content-addressed sidecar; never rewrite a historical report."""
    report = inspection.get("replacement_identity")
    if report is None:
        return inspection
    content = _encoded(report) + b"\n"
    path = Path(directory) / ("reaction-correspondence-" + _digest(report)[:20] + ".json")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("xb") as handle:
            handle.write(content)
    except FileExistsError:
        if path.read_bytes() != content:
            raise ValueError("Retained reaction correspondence was mutated")
    return {**inspection, "replacement_identity_artifact": _artifact(path)}


def _failed(reason: str, *, correspondence: dict | None = None) -> dict:
    return {"status": "failed", "blocking_reasons": [reason],
            "correspondence": correspondence}


def resolve_replacement_correspondence(evidence: dict, replacement: dict) -> dict:
    """Admit a current snapshot's structural proof, or explicitly withhold it.

    Historical evidence with no structural policy remains distinguishable as
    ``not_available``. An expected, failed or stale structural report may never
    fall back to label matching.
    """
    inspection = evidence.get("mechanism_inspection") or {}
    report = inspection.get("replacement_identity")
    required = bool(evidence.get("replacement_identity_required"))
    if "replacement_identity" not in inspection:
        return (_failed("replacement_identity_report_missing") if required else
                {"status": "not_available", "blocking_reasons": []})
    if not isinstance(report, dict):
        return _failed("replacement_identity_report_invalid")
    entry = None
    try:
        if (report.get("schema_version") != 1
                or report.get("kind") != "storca_replacement_correspondences"
                or report.get("policy_id") != IDENTITY_POLICY):
            return _failed("replacement_identity_policy_unsupported")
        saved = inspection.get("replacement_identity_artifact") or {}
        _check_artifact(saved)
        if json.loads(Path(saved["path"]).read_text()) != report:
            return _failed("replacement_identity_sidecar_mismatch")
        if report.get("inspection_sha256") != _digest(inspection.get("reactions") or []):
            return _failed("replacement_identity_inspection_changed")
        model_artifacts = report.get("mechanism_artifacts") or {}
        requested_artifacts = evidence.get("artifacts") or {}
        for name, field in (("chemkin", "chemkin"), ("dictionary", "species_dictionary")):
            if not requested_artifacts.get(field):
                return _failed("replacement_identity_model_artifacts_missing")
            _check_artifact(model_artifacts[name], Path(requested_artifacts[field]))
        if model_artifacts.get("transport"):
            _check_artifact(model_artifacts["transport"], Path(requested_artifacts["transport"]))
        route_id = str(replacement.get("identity_source_route_id") or replacement.get("route_id") or "")
        matches = [item for item in report.get("correspondences") or [] if item.get("route_id") == route_id]
        if not route_id or len(matches) != 1:
            return _failed("replacement_identity_missing_or_ambiguous")
        entry = matches[0]
        library = Path(replacement["library"]).resolve(strict=True)
        if Path(entry["library"]).resolve(strict=True) != library:
            return _failed("replacement_identity_source_library_mismatch", correspondence=entry)
        sources = entry["source_artifacts"]
        for name in SOURCE_FILES:
            _check_artifact(sources[name], library / name)
        manifest = json.loads((library / "generated-kinetics.json").read_text())
        if (manifest.get("kind") != "storca_generated_kinetics"
                or manifest.get("route_id") != route_id
                or not str(manifest.get("verification_status", "")).startswith("verified")):
            return _failed("replacement_identity_source_manifest_mismatch", correspondence=entry)
        for name in SOURCE_FILES[:2]:
            if (sources[name]["sha256"] != (manifest.get("checksums") or {}).get(name)
                    or sources[name]["sha256"] != (replacement.get("checksums") or {}).get(name)):
                return _failed("replacement_identity_source_checksums_mismatch", correspondence=entry)
        if (replacement.get("condition_rate") != manifest.get("condition_rate")
                or replacement.get("molecularity") != manifest.get("molecularity")):
            return _failed("replacement_identity_source_rate_contract_changed", correspondence=entry)
        if entry.get("status") != "matched" or entry.get("blocking_reasons"):
            return {"status": "failed", "blocking_reasons": entry.get("blocking_reasons") or
                    ["replacement_identity_unresolved"], "correspondence": entry}
        index = entry.get("matched_reaction_index")
        candidates = [item for item in inspection.get("reactions") or [] if item.get("reaction_index") == index]
        if type(index) is not int or len(candidates) != 1:
            return _failed("replacement_identity_model_reaction_ambiguous", correspondence=entry)
        candidate = candidates[0]
        direction = entry.get("applied_direction")
        if direction not in {"forward", "reverse"} or (direction == "reverse" and not candidate.get("reversible")):
            return _failed("replacement_identity_direction_invalid", correspondence=entry)
        if not entry.get("source_library") or candidate.get("source_library") != entry["source_library"]:
            return _failed("replacement_identity_library_provenance_mismatch", correspondence=entry)
        left, right = candidate["reactant_labels"], candidate["product_labels"]
        if direction == "reverse":
            left, right = right, left
        arrow = " <=> " if candidate.get("reversible") else " => "
        equation = " + ".join(left) + arrow + " + ".join(right)
        if re.sub(r"\s+", "", entry.get("source_oriented_model_equation", "")) != re.sub(r"\s+", "", equation):
            return _failed("replacement_identity_model_equation_mismatch", correspondence=entry)
        for side, labels in (("reactants", left), ("products", right)):
            mappings = (entry.get("species_correspondence") or {}).get(side) or []
            if (Counter(item["model_label"] for item in mappings) != Counter(labels)
                    or sorted(item["model_position"] for item in mappings) != list(range(len(labels)))
                    or sorted(item["source_position"] for item in mappings) != list(range(len(labels)))):
                return _failed("replacement_identity_stoichiometry_mismatch", correspondence=entry)
            if manifest.get(side) is not None and Counter(item["source_label"] for item in mappings) != Counter(manifest[side]):
                return _failed("replacement_identity_source_stoichiometry_mismatch", correspondence=entry)
        return {"status": "matched", "candidate": candidate, "applied_direction": direction,
                "source_oriented_equation": equation, "source_route_id": route_id,
                "correspondence": entry, "blocking_reasons": []}
    except (AttributeError, KeyError, TypeError, ValueError, OSError) as error:
        return {"status": "failed", "blocking_reasons": ["replacement_identity_evidence_invalid"],
                "failure_reason": str(error), "correspondence": entry}
