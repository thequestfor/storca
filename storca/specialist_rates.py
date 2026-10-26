"""Auditable boundary for non-harmonic-TST specialist rate providers.

The provider may wrap Arkane, MESMER, MultiWell, a VTST program, or retained
user/external evidence.  STORCA does not infer that an external result is valid
from its filename or solver name: the condition, applicability diagnostics,
uncertainty, artifacts, and generated network library must all be explicit.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
import shutil
from typing import Any, Mapping, Protocol

from .generated_kinetics import validate_rate_replacement_evidence
from .kinetic_uncertainty import RateUncertaintyRecord


SPECIALIST_RATE_REQUEST_SCHEMA_VERSION = 1
SPECIALIST_RATE_RESULT_SCHEMA_VERSION = 1


class SpecialistRateProvider(Protocol):
    """One condition-bound provider; implementations remain independently testable."""

    name: str

    def evaluate(self, request: Mapping[str, Any], workdir: Path) -> dict[str, Any]:
        ...


def _canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _condition_hash(condition_contract: Mapping[str, Any]) -> str:
    return _canonical_hash(dict(condition_contract))


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _resolve_bundle_path(raw_path: object, bundle_root: Path) -> Path:
    path = Path(str(raw_path or ""))
    return path if path.is_absolute() else (bundle_root / path).resolve()


def _normalized_result_paths(result: Mapping[str, Any], bundle_root: Path) -> dict[str, Any]:
    """Resolve only declared artifact paths relative to a retained result bundle."""
    value = dict(result)
    artifacts = []
    for raw in value.get("evidence_artifacts") or ():
        item = dict(raw or {})
        item["path"] = str(_resolve_bundle_path(item.get("path"), bundle_root))
        artifacts.append(item)
    value["evidence_artifacts"] = artifacts
    replacement = dict(value.get("generated_kinetics_library") or {})
    if replacement.get("library"):
        replacement["library"] = str(
            _resolve_bundle_path(replacement["library"], bundle_root)
        )
    value["generated_kinetics_library"] = replacement
    return value


@dataclass(frozen=True)
class RetainedSpecialistRateProvider:
    """Import one or more explicit specialist results without rerunning chemistry.

    A file supplies one result.  A directory may contain multiple JSON results;
    exactly one must match the immutable request hash.  Relative evidence and
    library paths are resolved against the submitted result file, never against
    the caller's current working directory.
    """

    source: Path
    name: str = "retained-specialist-rate-result-v1"

    def _candidates(self) -> tuple[Path, ...]:
        source = Path(self.source).expanduser().resolve()
        if source.is_file():
            return (source,)
        if source.is_dir():
            return tuple(sorted(path for path in source.rglob("*.json") if path.is_file()))
        raise FileNotFoundError(f"Specialist-rate result source does not exist: {source}")

    def evaluate(self, request: Mapping[str, Any], workdir: Path) -> dict[str, Any]:
        matches: list[tuple[Path, dict[str, Any]]] = []
        malformed: list[str] = []
        for path in self._candidates():
            try:
                parsed = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError) as error:
                malformed.append(f"{path}: {type(error).__name__}")
                continue
            if not isinstance(parsed, dict):
                continue
            if parsed.get("kind") == "storca_specialist_rate_result" \
                    and parsed.get("request_sha256") == request.get("request_sha256"):
                matches.append((path, parsed))
        if len(matches) != 1:
            detail = f"; malformed candidates: {', '.join(malformed)}" if malformed else ""
            raise ValueError(
                "Expected exactly one specialist result matching request "
                f"{request.get('request_sha256')}; found {len(matches)}{detail}"
            )
        source_path, raw_result = matches[0]
        result = _normalized_result_paths(raw_result, source_path.parent)
        result["import_provenance"] = {
            "kind": "retained_specialist_rate_import",
            "provider_adapter": self.name,
            "source_result": str(source_path),
            "source_result_sha256": _file_sha256(source_path),
            "path_resolution_root": str(source_path.parent),
        }
        validation = validate_specialist_rate_result(result, request)
        if validation["status"] != "passed":
            raise ValueError(
                "Submitted specialist result failed validation: "
                + ", ".join(validation["blocking_reasons"])
            )
        destination = Path(workdir)
        destination.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination / "submitted-result.json")
        (destination / "import-provenance.json").write_text(
            json.dumps(result["import_provenance"], indent=2, sort_keys=True) + "\n"
        )
        return result


def build_specialist_rate_request(
    route: Mapping[str, Any],
    path_result: Mapping[str, Any],
    condition_contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Build an immutable request from a fail-closed path-method escalation."""
    verification = dict(path_result.get("route_verification") or path_result)
    status = str(verification.get("status") or "")
    classification = str(verification.get("path_classification") or "")
    if status == "dissociation_rate_model_required":
        problem_class = "unimolecular_dissociation"
        allowed = (
            "variational_transition_state_theory",
            "rrkm_master_equation",
            "validated_external_falloff_model",
        )
    elif status == "barrierless_rate_model_required" and "capture" in classification:
        problem_class = "barrierless_association"
        allowed = (
            "capture_theory_with_recrossing",
            "variational_transition_state_theory",
            "validated_external_falloff_model",
        )
    elif status == "barrierless_rate_model_required":
        problem_class = "barrierless_or_loose_path"
        allowed = (
            "variational_transition_state_theory",
            "rrkm_master_equation",
            "validated_external_rate_model",
        )
    else:
        raise ValueError(
            "A specialist rate request requires an explicit barrierless or dissociation path escalation"
        )
    route_id = str(route.get("route_id") or "").strip()
    if not route_id:
        raise ValueError("Specialist rate request needs a route ID")
    reactants = list(route.get("reactant_smiles") or ())
    products = list(route.get("product_smiles") or ())
    if not reactants or not products:
        raise ValueError("Specialist rate request needs explicit reactant and product identities")
    body = {
        "schema_version": SPECIALIST_RATE_REQUEST_SCHEMA_VERSION,
        "kind": "storca_specialist_rate_request",
        "route_id": route_id,
        "reaction_equation": route.get("reaction_equation"),
        "reactant_smiles": reactants,
        "product_smiles": products,
        "reactant_labels": list(route.get("reactant_labels") or ()),
        "product_labels": list(route.get("product_labels") or ()),
        "charge": route.get("charge"),
        "multiplicity": route.get("multiplicity"),
        "problem_class": problem_class,
        "path_status": status,
        "path_classification": classification,
        "allowed_rate_theories": list(allowed),
        "condition_contract": dict(condition_contract),
        "condition_contract_sha256": _condition_hash(condition_contract),
        "path_evidence_sha256": _canonical_hash(verification),
        "method_selection": verification.get("method_selection"),
        "claims_permitted": ["request_specialist_rate_evidence"],
        "claims_prohibited": ["rate", "lifetime"],
    }
    return {**body, "request_sha256": _canonical_hash(body)}


def validate_specialist_rate_request(request: Mapping[str, Any]) -> dict[str, Any]:
    """Verify that neither the request body nor its condition was altered."""
    value = dict(request)
    blockers: list[str] = []
    if value.get("kind") != "storca_specialist_rate_request":
        blockers.append("specialist_request_kind_invalid")
    if value.get("schema_version") != SPECIALIST_RATE_REQUEST_SCHEMA_VERSION:
        blockers.append("specialist_request_schema_unsupported")
    unsigned = {key: item for key, item in value.items() if key != "request_sha256"}
    if value.get("request_sha256") != _canonical_hash(unsigned):
        blockers.append("specialist_request_content_hash_mismatch")
    condition = value.get("condition_contract")
    if not isinstance(condition, dict) \
            or value.get("condition_contract_sha256") != _condition_hash(condition):
        blockers.append("specialist_request_condition_hash_mismatch")
    if not str(value.get("route_id") or "").strip():
        blockers.append("specialist_request_route_missing")
    if not value.get("allowed_rate_theories"):
        blockers.append("specialist_request_allowed_theories_missing")
    return {
        "schema_version": 1,
        "kind": "storca_specialist_rate_request_validation",
        "status": "passed" if not blockers else "failed",
        "request_sha256": value.get("request_sha256"),
        "blocking_reasons": blockers,
    }


def validate_specialist_rate_result(
    result: Mapping[str, Any],
    request: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate specialist evidence before it can enter the ordinary repair loop."""
    value = dict(result)
    request_validation = validate_specialist_rate_request(request)
    blockers: list[str] = list(request_validation["blocking_reasons"])
    if value.get("kind") != "storca_specialist_rate_result":
        blockers.append("specialist_result_kind_invalid")
    if value.get("schema_version") != SPECIALIST_RATE_RESULT_SCHEMA_VERSION:
        blockers.append("specialist_result_schema_unsupported")
    if value.get("status") != "verified":
        blockers.append("specialist_result_not_verified")
    if value.get("request_sha256") != request.get("request_sha256"):
        blockers.append("specialist_request_hash_mismatch")
    if value.get("condition_contract_sha256") != request.get("condition_contract_sha256"):
        blockers.append("specialist_condition_contract_mismatch")
    theory = str(value.get("rate_theory") or "")
    if theory not in set(request.get("allowed_rate_theories") or ()):
        blockers.append("rate_theory_not_allowed_for_diagnosed_path")
    provider = value.get("provider") or {}
    if not all(str(provider.get(key) or "").strip() for key in ("name", "version")):
        blockers.append("specialist_provider_identity_incomplete")
    applicability = value.get("applicability") or {}
    if applicability.get("status") != "applicable":
        blockers.append("specialist_method_applicability_not_established")
    if not applicability.get("diagnostics"):
        blockers.append("specialist_applicability_diagnostics_missing")
    uncertainty = value.get("uncertainty") or {}
    try:
        lower = float(uncertainty.get("lower_rate_multiplier"))
        upper = float(uncertainty.get("upper_rate_multiplier"))
    except (TypeError, ValueError):
        lower = upper = math.nan
    if not (math.isfinite(lower) and math.isfinite(upper) and 0.0 < lower <= 1.0 <= upper):
        blockers.append("specialist_rate_uncertainty_bounds_invalid")
    if not str(uncertainty.get("source") or "").strip():
        blockers.append("specialist_rate_uncertainty_source_missing")
    if not str(uncertainty.get("applicability") or "").strip():
        blockers.append("specialist_rate_uncertainty_applicability_missing")
    if not str(uncertainty.get("source_kind") or "").strip():
        blockers.append("specialist_rate_uncertainty_source_kind_missing")
    components = uncertainty.get("uncertainty_components") or ()
    if not isinstance(components, (list, tuple)) or not components \
            or any(not str(item).strip() for item in components):
        blockers.append("specialist_rate_uncertainty_components_missing")

    artifact_checks: list[dict[str, Any]] = []
    artifacts = value.get("evidence_artifacts") or ()
    if not artifacts:
        blockers.append("specialist_evidence_artifacts_missing")
    for item in artifacts:
        path = Path(str((item or {}).get("path") or ""))
        expected = str((item or {}).get("sha256") or "")
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        valid = bool(actual and expected == actual)
        artifact_checks.append({"path": str(path), "expected_sha256": expected, "actual_sha256": actual, "valid": valid})
        if not valid:
            blockers.append("specialist_evidence_artifact_missing_or_mutated")

    replacement = dict(value.get("generated_kinetics_library") or {})
    replacement_validation = validate_rate_replacement_evidence(
        replacement,
        temperature_K=float(request["condition_contract"]["temperature_K"]),
        pressure_bar=float(request["condition_contract"]["pressure_bar"]),
    )
    if replacement_validation.get("status") != "passed":
        blockers.extend(
            "replacement_" + reason
            for reason in replacement_validation.get("blocking_reasons") or ()
        )
    if replacement.get("route_id") != request.get("route_id"):
        blockers.append("specialist_replacement_route_mismatch")
    library = Path(str(replacement.get("library") or ""))
    if not library.is_dir():
        blockers.append("specialist_replacement_library_missing")
    checksums = replacement.get("checksums") or {}
    if not isinstance(checksums, dict) or not checksums:
        blockers.append("specialist_replacement_library_checksums_missing")
    else:
        for name, expected in checksums.items():
            relative = Path(str(name))
            if relative.is_absolute() or ".." in relative.parts:
                blockers.append("specialist_replacement_library_checksum_path_invalid")
                continue
            artifact = library / relative
            actual = _file_sha256(artifact) if artifact.is_file() else None
            if not actual or actual != str(expected):
                blockers.append("specialist_replacement_library_missing_or_mutated")

    if theory in {"rrkm_master_equation", "validated_external_falloff_model"}:
        pdep = value.get("pressure_dependence") or {}
        pressure_range = pdep.get("pressure_range_bar") or ()
        condition_pressure = float(request["condition_contract"]["pressure_bar"])
        if (
            len(pressure_range) != 2
            or not float(pressure_range[0]) <= condition_pressure <= float(pressure_range[1])
        ):
            blockers.append("specialist_pressure_grid_does_not_cover_condition")
        if not pdep.get("bath_gas"):
            blockers.append("specialist_bath_gas_missing")
        if not pdep.get("energy_transfer_model"):
            blockers.append("specialist_energy_transfer_model_missing")
    if theory == "variational_transition_state_theory":
        surface = value.get("variational_dividing_surface") or {}
        if not isinstance(surface, Mapping) or not surface:
            blockers.append("variational_dividing_surface_evidence_missing")
        else:
            if not str(surface.get("reaction_coordinate") or "").strip():
                blockers.append("variational_reaction_coordinate_missing")
            if not surface.get("sampled_dividing_surfaces"):
                blockers.append("variational_surface_scan_missing")
            if not str(surface.get("variational_criterion") or "").strip():
                blockers.append("variational_criterion_missing")

    if theory == "capture_theory_with_recrossing":
        capture = value.get("capture_model") or {}
        if not isinstance(capture, Mapping) or not capture:
            blockers.append("capture_model_evidence_missing")
        else:
            long_range = capture.get("long_range_interaction") or {}
            if not isinstance(long_range, Mapping) or not str(
                long_range.get("model") or ""
            ).strip():
                blockers.append("capture_long_range_interaction_model_missing")
            elif not str(long_range.get("parameter_source") or "").strip():
                blockers.append("capture_long_range_parameter_source_missing")

            surface = capture.get("dividing_surface") or {}
            if not isinstance(surface, Mapping) or not surface.get("sampled_surfaces"):
                blockers.append("capture_dividing_surface_scan_missing")
            elif not str(surface.get("variational_criterion") or "").strip():
                blockers.append("capture_variational_criterion_missing")

            recrossing = capture.get("recrossing") or {}
            if not isinstance(recrossing, Mapping) or not str(
                recrossing.get("model") or ""
            ).strip():
                blockers.append("capture_recrossing_model_missing")
            elif not str(recrossing.get("evidence") or "").strip():
                blockers.append("capture_recrossing_evidence_missing")

            for field, reason in (
                ("electronic_degeneracy_factor", "capture_electronic_degeneracy_missing"),
                ("symmetry_factor", "capture_symmetry_factor_missing"),
            ):
                try:
                    factor = float(capture.get(field))
                except (TypeError, ValueError):
                    factor = math.nan
                if not math.isfinite(factor) or factor <= 0.0:
                    blockers.append(reason)

            spin = capture.get("electronic_state_assessment") or {}
            if not isinstance(spin, Mapping) or spin.get("status") != "applicable":
                blockers.append("capture_electronic_state_applicability_missing")
            elif not str(spin.get("basis") or "").strip():
                blockers.append("capture_electronic_state_basis_missing")

    blockers = list(dict.fromkeys(blockers))
    return {
        "schema_version": 1,
        "kind": "storca_specialist_rate_validation",
        "status": "passed" if not blockers else "failed",
        "request_sha256": request.get("request_sha256"),
        "route_id": request.get("route_id"),
        "rate_theory": theory or None,
        "artifact_checks": artifact_checks,
        "request_validation": request_validation,
        "replacement_validation": replacement_validation,
        "blocking_reasons": blockers,
        "rate_claim_allowed": not blockers,
        "lifetime_claim_allowed": False,
    }


def specialist_rate_uncertainty_record(
    result: Mapping[str, Any], request: Mapping[str, Any],
) -> RateUncertaintyRecord:
    """Convert already-validated provider uncertainty into allocation evidence."""
    uncertainty = dict(result.get("uncertainty") or {})
    return RateUncertaintyRecord(
        route_id=str(request["route_id"]),
        reaction_equation=str(request.get("reaction_equation") or "") or None,
        lower_rate_multiplier=uncertainty["lower_rate_multiplier"],
        upper_rate_multiplier=uncertainty["upper_rate_multiplier"],
        source=uncertainty["source"],
        source_kind=uncertainty["source_kind"],
        applicability=uncertainty["applicability"],
        # Preserve the declaration's type so the record's method-calibration
        # admission cannot be bypassed by stringifying an invalid request ID.
        condition_contract_id=(request.get("condition_contract") or {}).get("contract_id"),
        confidence_level=uncertainty.get("confidence_level"),
        uncertainty_components=tuple(uncertainty.get("uncertainty_components") or ()),
    )


__all__ = [
    "SPECIALIST_RATE_REQUEST_SCHEMA_VERSION",
    "SPECIALIST_RATE_RESULT_SCHEMA_VERSION",
    "RetainedSpecialistRateProvider",
    "SpecialistRateProvider",
    "build_specialist_rate_request",
    "specialist_rate_uncertainty_record",
    "validate_specialist_rate_request",
    "validate_specialist_rate_result",
]
