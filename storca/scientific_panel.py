"""Bounded, offline development-panel scoring; never a chemistry dispatcher.

Calculation requests are allowlisted projections, separate from scoring-only
rates, products and mechanisms. Pending references, failed calculations and
abstentions stay in the ledger. These descriptive scores grant neither support
nor uncertainty calibration, even when every available comparator agrees.
"""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import median
from typing import Any, Mapping, Sequence


SCHEMA_VERSION = 1
MAXIMUM_CASE_COUNT = 12
MAXIMUM_RATE_OBSERVATIONS = 12
_DISPOSITIONS = {"prediction", "unresolved", "outside_applicability", "crashed"}
_CONDITION_FIELDS = {
    "phase", "environment", "temperature_K", "pressure_bar", "pressure_regime", "reactor",
    "initial_composition", "added_bath_gas", "duration_s", "retention_fraction",
    "solvent", "solvent_model", "standard_state", "light_condition",
    "initial_concentrations_mol_L",
}
_PARENT_FIELDS = {"name", "smiles", "charge", "multiplicity"}
_BUDGET_FIELDS = {
    "overall_timeout_seconds", "cleanup_grace_seconds", "maximum_cores",
    "maximum_high_level_jobs", "route_timeout_seconds_per_job",
}


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode()).hexdigest()


def _number(value: Any, *, positive: bool = False) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    try:
        number = float(value)
    except (OverflowError, ValueError):
        return None
    return number if math.isfinite(number) and (number > 0 if positive else number >= 0) else None


def _project_fields(value: Mapping[str, Any], fields: set[str]) -> dict[str, Any]:
    return {field: deepcopy(value[field]) for field in sorted(fields & value.keys())}


def _project_scalars(value: Mapping[str, Any], fields: set[str]) -> dict[str, Any]:
    projected = _project_fields(value, fields)
    if any(isinstance(item, (dict, list, tuple)) for item in projected.values()):
        raise ValueError("Request scalar fields cannot contain nested reference records")
    return projected


def prediction_request(case: Mapping[str, Any]) -> dict[str, Any]:
    """Export inputs, never reference products/rates, acceptance answers or sources.

    Nested composition/species records are projected too: a reference concealed
    inside a conditions record must not travel to a calculation backend.
    """
    if case.get("rate_only_scoring") is not None:
        return _rate_prediction_request(case)
    supplied = case.get("prediction_request") or {}
    declared = supplied.get("conditions") or {}
    conditions = _project_scalars(declared, _CONDITION_FIELDS - {
        "initial_composition", "added_bath_gas", "initial_concentrations_mol_L",
    })
    if "initial_composition" in declared:
        conditions["initial_composition"] = [
            _project_scalars(item, {"smiles", "mole_fraction", "charge", "multiplicity"})
            for item in declared["initial_composition"]
        ]
    if isinstance(declared.get("added_bath_gas"), list):
        conditions["added_bath_gas"] = [
            _project_scalars(item, {"smiles", "mole_fraction", "charge", "multiplicity"})
            for item in declared["added_bath_gas"]
        ]
    elif "added_bath_gas" in declared:
        conditions.update(_project_scalars(declared, {"added_bath_gas"}))
    if "initial_concentrations_mol_L" in declared:
        concentrations = declared["initial_concentrations_mol_L"]
        if not isinstance(concentrations, dict) or any(_number(value) is None for value in concentrations.values()):
            raise ValueError("Concentration requests require a species-to-number mapping")
        conditions["initial_concentrations_mol_L"] = deepcopy(concentrations)
    return {
        "schema_version": SCHEMA_VERSION,
        "kind": "storca_scientific_panel_prediction_request",
        "case_id": case["case_id"],
        "chemical_system": {"parent": _project_scalars(
            (supplied.get("chemical_system") or {}).get("parent") or {}, _PARENT_FIELDS,
        )},
        "conditions": conditions,
        "requested_observables": list(supplied.get("requested_observables") or ()),
        "budget": _project_scalars(supplied.get("budget") or {}, _BUDGET_FIELDS),
        # Projection is preparation; existing lane/driver gates still admit execution.
        "execution_admitted": False,
    }


def case_from_p1_preregistration(path: Path) -> dict[str, Any]:
    """Reuse the frozen P1 comparator without claiming its pending audit passed."""
    payload = json.loads(Path(path).read_text())
    if (payload.get("schema_version") != 1
            or payload.get("kind") != "storca_p1_development_preregistration"
            or payload.get("predictions_inspected_at_registration") is not False
            or payload.get("partition") != "development"
            or payload["chemical_system"].get("pathway_injection_allowed") is not False):
        raise ValueError("A pre-prediction, development-only P1 preregistration is required")
    reference = payload["reference"]
    comparator = reference["primary_numeric_comparator"]
    accuracy = payload["acceptance"]["accuracy"]
    rate = comparator["k_s^-1"]
    if rate != accuracy["reference"] or _number(rate, positive=True) is None:
        raise ValueError("P1 comparator and registered acceptance reference disagree")
    derived = reference["derived_comparators"]
    audited = (reference.get("transcription") or {}).get("visual_page_verified") is True
    return {
        "case_id": payload["case_id"], "partition": "development",
        "classes": ["gas", "ordinary_thermal", "unimolecular_isomerization"],
        "prediction_request": {
            "chemical_system": {"parent": payload["chemical_system"]["parent"]},
            "conditions": payload["conditions"],
            "requested_observables": payload["requested_observables"],
        },
        "scoring": {
            "reference_status": "audited" if audited else "provisional",
            "reference_audit_pending": not audited,
            "expected_disposition": "prediction",
            # An expected P1 chain does not establish calibrated method applicability.
            "applicability_assessment": "unresolved",
            "rate": {"value": rate, "units": "s^-1", "convention": reference["rate_convention"]},
            "lifetime": {"value": derived["t95_s"], "units": "s", "kind": "derived_comparator",
                         "assumption": derived["assumption"]},
            "parent_fraction": {"value": derived["parent_fraction_at_3600_s"],
                                "kind": "derived_comparator", "assumption": derived["assumption"]},
            "dominant_product": payload["chemical_system"]["expected_major_product_for_scoring_only"]["smiles"],
            "maximum_rate_multiplicative_error": accuracy["maximum_multiplicative_error"],
            "experimental_uncertainty": comparator.get("experimental_total_uncertainty"),
            "provenance": {key: reference[key] for key in (
                "citation", "publisher_url", "source_kind", "source_location",
            )},
            # Product formation is not independent elementary-mechanism ground truth.
            "mechanism_reactions": None, "consequential_reactions": None,
        },
    }


def _validate_reference_audit(case: Mapping[str, Any], manifest_path: Path) -> dict[str, Any] | None:
    """Bind optional frozen transcriptions to retained primary sources.

    A source audit is not scoring readiness or execution admission. This only
    checks the declared immutable dependencies and refuses promotion of an
    explicitly unready observation set into a scoring record.
    """
    record = case.get("reference_audit")
    if record is None:
        return
    try:
        if not isinstance(record, dict):
            raise ValueError("Invalid reference audit record")
        audit_path = Path(record["path"])
        if not audit_path.is_absolute():
            audit_path = Path(manifest_path).parent / audit_path
        content = audit_path.read_bytes()
        if hashlib.sha256(content).hexdigest() != record.get("sha256"):
            raise ValueError("Changed reference audit artifact")
        audit = json.loads(content)
        if (audit.get("schema_version") != 1
                or audit.get("kind") != "storca_scientific_panel_reference_audit"
                or any(audit.get(key) is not False for key in (
                    "prediction_inputs_allowed", "support_claim_allowed", "calibration_claim_allowed",
                ))):
            raise ValueError("Unsupported reference audit claims or schema")
        matches = [item for item in audit["audited_observation_sets"]
                   if item["case_id"] == case["case_id"]]
        if len(matches) != 1:
            raise ValueError("Unique reference audit observation identity required")
        observation = matches[0]
        if case.get("scoring") is not None and observation.get("quantitative_scoring_ready") is not True:
            raise ValueError("Audited observations are not scoring ready")
        binding = observation.get("scoring_binding")
        if case.get("scoring") is not None and binding is not None:
            # Freeze a source transcription, not merely a readiness boolean.
            # Optional bindings leave earlier audits/gates unchanged.
            if not isinstance(binding, dict):
                raise ValueError("Invalid reference scoring binding")
            scoring = case["scoring"]
            request = prediction_request(case)
            if scoring.get("rate") != binding.get("rate"):
                raise ValueError("Reference scoring rate differs from frozen audit")
            if ("reference_uncertainty" in observation
                    and scoring.get("reference_uncertainty") != observation["reference_uncertainty"]):
                raise ValueError("Reference uncertainty differs from frozen audit")
            if request["conditions"] != binding.get("conditions"):
                raise ValueError("Reference prediction conditions differ from frozen audit")
            if ("chemical_system" in binding
                    and request["chemical_system"] != binding["chemical_system"]):
                raise ValueError("Reference prediction identity differs from frozen audit")
            if request["requested_observables"] != binding.get("requested_observables"):
                raise ValueError("Reference requested observables differ from frozen audit")
            permitted = binding.get("permitted_scoring_observables")
            if not isinstance(permitted, list) or any(
                scoring.get(name) is not None and name not in permitted
                for name in ("rate", "lifetime", "parent_fraction", "dominant_product", "branching",
                             "mechanism_reactions", "consequential_reactions")
            ):
                raise ValueError("Reference scoring observables exceed frozen audit scope")
            constraints = binding.get("scoring_constraints")
            if not isinstance(constraints, dict) or any(
                scoring.get(name) != value for name, value in constraints.items()
            ):
                raise ValueError("Reference scoring constraints differ from frozen audit")
        source_ids = observation["source_ids"]
        sources = audit["sources"]
        if (not isinstance(source_ids, list) or not source_ids
                or len(set(source_ids)) != len(source_ids)
                or len({item["source_id"] for item in sources}) != len(sources)):
            raise ValueError("Unique nonempty reference source identities required")
        for identifier in source_ids:
            source_matches = [item for item in sources if item["source_id"] == identifier]
            if len(source_matches) != 1:
                raise ValueError("Missing reference source identity")
            source = source_matches[0]
            source_path = Path(source["path"])
            if not source_path.is_absolute():
                source_path = audit_path.parent / source_path
            if hashlib.sha256(source_path.read_bytes()).hexdigest() != source.get("sha256"):
                raise ValueError("Changed reference source artifact")
        return observation
    except (OSError, KeyError, TypeError, AttributeError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid or unavailable reference audit for {case['case_id']}: {error}") from error


def _rate_prediction_request(case: Mapping[str, Any]) -> dict[str, Any]:
    """Only genuine starting identities and physical rate conditions may leave scoring."""
    supplied = case["rate_only_request"]
    declared = supplied.get("conditions") or {}
    fields = {"phase", "solvent", "solvent_isotopic_formula", "pressure_bar",
              "atmosphere", "light_condition", "standard_state_mol_L"}
    conditions = _project_scalars(declared, fields)
    grid = declared.get("temperature_grid_K")
    if (not isinstance(grid, list) or not 0 < len(grid) <= MAXIMUM_RATE_OBSERVATIONS
            or any(_number(value, positive=True) is None for value in grid)
            or len(set(grid)) != len(grid)):
        raise ValueError("Rate requests need a unique finite positive temperature grid")
    conditions["temperature_grid_K"] = list(grid)
    species = (supplied.get("chemical_system") or {}).get("initial_species")
    if not isinstance(species, list) or not 0 < len(species) <= MAXIMUM_CASE_COUNT:
        raise ValueError("Rate requests need declared starting molecular identities")
    initial = [_project_scalars(item, _PARENT_FIELDS) for item in species]
    if any(not item.get("smiles") or type(item.get("charge")) is not int or item["charge"] != 0
           or type(item.get("multiplicity")) is not int or item["multiplicity"] != 1 for item in initial):
        raise ValueError("Rate request starting identities need neutral singlet SMILES")
    return {
        "schema_version": SCHEMA_VERSION, "kind": "storca_scientific_panel_rate_request",
        "case_id": case["case_id"], "chemical_system": {"initial_species": initial},
        "conditions": conditions, "requested_observables": ["condition_matched_rate_coefficients"],
        # No budget is granted by this offline reference protocol.
        "budget": {}, "execution_admitted": False,
    }


def _attach_rate_preregistrations(cases: list[dict[str, Any]], paths: Sequence[Path]) -> None:
    """Attach immutable descriptive overlays, preserving ordinary/full-case gates."""
    if any("rate_only_scoring" in case or "rate_only_request" in case for case in cases):
        raise ValueError("Rate-only scoring requires a separate preregistration, not a manifest injection")
    attached: set[str] = set()
    for path in paths:
        path = Path(path)
        raw = path.read_bytes()
        protocol = json.loads(raw)
        identifier = protocol.get("case_id")
        matches = [case for case in cases if case.get("case_id") == identifier]
        if (protocol.get("schema_version") != SCHEMA_VERSION
                or protocol.get("kind") != "storca_rate_only_development_preregistration"
                or protocol.get("partition") != "development"
                or protocol.get("predictions_inspected_at_registration") is not False
                or any(protocol.get(key) is not False for key in (
                    "holdout_locked", "support_claim_allowed", "calibration_claim_allowed"))
                or len(matches) != 1 or identifier in attached):
            raise ValueError("Unique pre-prediction development-only rate preregistration required")
        case = matches[0]
        if case.get("scoring") is not None or case.get("partition") != "development":
            raise ValueError("Rate-only overlays cannot replace full-case scoring or holdouts")
        audit_record = protocol.get("reference_audit") or {}
        if audit_record.get("sha256") != (case.get("reference_audit") or {}).get("sha256"):
            raise ValueError("Rate preregistration source audit mismatch")
        observation = _validate_reference_audit({**case, "reference_audit": audit_record}, path)
        if (not observation or observation.get("transcription_frozen") is not True
                or observation.get("source_fitted_rate_only_comparison_available") is not True):
            raise ValueError("Source-fitted rate-only reference readiness required")
        references = observation.get("observations") or []
        identifiers = protocol.get("observation_ids")
        if (not isinstance(identifiers, list) or not 0 < len(identifiers) <= MAXIMUM_RATE_OBSERVATIONS
                or len(set(identifiers)) != len(identifiers)
                or set(identifiers) != {item["observation_id"] for item in references}):
            raise ValueError("All audited rate observations must remain uniquely registered")
        if any(_number(ref.get("rate"), positive=True) is None
               or _number(ref.get("joint_fit_interval_half_width")) is None
               or ref.get("direction") not in {"forward", "reverse"}
               or ref.get("rate_units") != ("L mol^-1 s^-1" if ref.get("direction") == "forward" else "s^-1")
               for ref in references):
            raise ValueError("Positive source coefficients, molecularity units and fit widths required")
        policy = protocol.get("scoring_policy") or {}
        if (policy.get("kind") != "descriptive_rate_only_development"
                or policy.get("rate_convention") != "concentration_mass_action"
                or policy.get("exact_isotopic_solvent_required") is not True
                or policy.get("all_registered_rates_remain_in_denominator") is not True
                or any(policy.get(key) is not None for key in (
                    "accuracy_threshold", "coverage_threshold", "cost_threshold", "scientific_acceptance_pass"))
                or any(policy.get(key) is not False for key in (
                    "experimental_fit_intervals_may_define_prediction_error", "macroscopic_observables_admitted"))):
            raise ValueError("Explicit descriptive-only rate scoring policy required")
        case["rate_only_request"] = deepcopy(protocol["prediction_request"])
        request = _rate_prediction_request(case)
        conditions = observation["source_conditions"]
        if (request["conditions"].get("phase") != "solution"
                or request["conditions"].get("solvent_isotopic_formula") != conditions["solvent_isotopic_formula"]
                or request["conditions"].get("solvent") != conditions["solvent"]
                or request["conditions"]["temperature_grid_K"] != conditions["temperature_grid_K"]
                or (protocol["prediction_request"].get("chemical_system") or {}).get("pathway_injection_allowed") is not False):
            raise ValueError("Rate preregistration must preserve source conditions and prohibit pathway injection")
        case["rate_only_scoring"] = {
            "references": deepcopy(references), "source_fit_uncertainty": deepcopy(observation["fit_uncertainty"]),
            "experimental_total_uncertainty": observation.get("experimental_total_uncertainty"),
            "policy": deepcopy(policy), "preregistration_sha256": hashlib.sha256(raw).hexdigest(),
            "source_audit_sha256": audit_record["sha256"],
        }
        attached.add(identifier)


def load_scientific_panel_manifest(
    path: Path, *, preregistration_paths: Sequence[Path] = (),
    rate_preregistration_paths: Sequence[Path] = (),
) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    if (payload.get("schema_version") != SCHEMA_VERSION
            or payload.get("kind") != "storca_scientific_panel_curation"):
        raise ValueError("Unsupported scientific-panel curation manifest")
    cases = list(payload.get("cases") or ())
    cases.extend(case_from_p1_preregistration(item) for item in preregistration_paths)
    if len(cases) > MAXIMUM_CASE_COUNT:
        raise ValueError(f"Offline panel slice is limited to {MAXIMUM_CASE_COUNT} cases")
    _attach_rate_preregistrations(cases, rate_preregistration_paths)
    identifiers: set[str] = set()
    for case in cases:
        identifier = str(case.get("case_id") or "")
        if not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_.-]*", identifier) or identifier in identifiers:
            raise ValueError("Panel case IDs must be unique filename-safe identifiers")
        identifiers.add(identifier)
        if case.get("partition") not in {"development", "holdout", "challenge"}:
            raise ValueError(f"Invalid partition for {identifier}")
        _validate_reference_audit(case, Path(path))
        scoring = case.get("scoring")
        if scoring is None:
            continue
        if case["partition"] == "holdout":
            raise ValueError("Prospective holdout scoring requires a separate locked blinded protocol")
        if scoring.get("reference_status") not in {"audited", "provisional"}:
            raise ValueError(f"Explicit reference readiness required for {identifier}")
        if scoring.get("expected_disposition") not in {"prediction", "outside_applicability"}:
            raise ValueError(f"Expected applicability required for {identifier}")
        if scoring.get("applicability_assessment", "unresolved") not in {
            "applicable", "outside_applicability", "unresolved",
        }:
            raise ValueError(f"Explicit applicability assessment invalid for {identifier}")
        if not scoring.get("provenance"):
            raise ValueError(f"Reference provenance required for {identifier}")
        request = prediction_request(case)
        conditions = request["conditions"]
        if not request["chemical_system"]["parent"].get("smiles") or not conditions:
            raise ValueError(f"Prediction identity and conditions required for {identifier}")
        if _number(conditions.get("temperature_K"), positive=True) is None:
            raise ValueError(f"Finite positive temperature required for {identifier}")
        for name, value in request["budget"].items():
            if _number(value, positive=name != "cleanup_grace_seconds") is None:
                raise ValueError(f"Finite bounded resource required for {identifier}/{name}")
        for name in ("rate", "lifetime"):
            if scoring.get(name):
                item = scoring[name]
                if _number(item.get("value"), positive=True) is None or not item.get("units"):
                    raise ValueError(f"Positive {name} and explicit units required for {identifier}")
        if scoring.get("rate") and not scoring["rate"].get("convention"):
            raise ValueError(f"Rate convention required for {identifier}")
        if scoring.get("lifetime") and scoring["lifetime"].get("kind") not in {
            "derived_comparator", "experimental_measurement",
        }:
            raise ValueError(f"Lifetime reference origin required for {identifier}")
        if scoring.get("parent_fraction"):
            fraction = _number(scoring["parent_fraction"].get("value"))
            retention = _number(conditions.get("retention_fraction"), positive=True)
            if fraction is None or fraction > 1 or retention is None or retention > 1:
                raise ValueError(f"Reference fraction and retention threshold required for {identifier}")
        if scoring.get("branching"):
            fractions = scoring["branching"].values()
            if (any(_number(value) is None or value > 1 for value in fractions)
                    or not math.isclose(sum(scoring["branching"].values()), 1.0, abs_tol=1e-6)):
                raise ValueError(f"Reference branching must sum to one for {identifier}")
        threshold = scoring.get("maximum_rate_multiplicative_error")
        if threshold is not None and (_number(threshold, positive=True) is None or threshold < 1):
            raise ValueError(f"Multiplicative error threshold must be at least one for {identifier}")
    return {**payload, "cases": cases}


def _artifacts_valid(records: Sequence[Mapping[str, Any]], root: Path) -> bool:
    if not isinstance(records, (list, tuple)) or not records:
        return False
    for record in records:
        if not isinstance(record, dict):
            return False
        path = Path(str(record.get("path") or ""))
        path = path if path.is_absolute() else root / path
        try:
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != record.get("sha256"):
                return False
        except OSError:
            return False
    return True


def _compute_score(result: Mapping[str, Any], root: Path) -> dict[str, Any]:
    compute = result.get("compute") or {}
    if not isinstance(compute, dict):
        compute = {}
    basis = compute.get("measurement_basis") or {}
    if not isinstance(basis, dict):
        basis = {}
    valid_evidence = _artifacts_valid(compute.get("evidence_artifacts") or (), root)
    measurements, missing = {}, []
    for name in ("wall_time_seconds", "core_time_seconds", "executed_high_level_jobs"):
        value = _number(compute.get(name))
        if (value is None or not valid_evidence
                or not isinstance(basis.get(name), str) or not basis[name].strip()
                or (name == "executed_high_level_jobs" and value != int(value))):
            missing.append(name)
        else:
            measurements[name] = int(value) if name == "executed_high_level_jobs" else value
    return {"measurements": measurements, "missing_measurements": missing,
            "evidence_valid": valid_evidence,
            "evidence_artifacts": [{"path": str(ref["path"]), "sha256": ref["sha256"]}
                                   for ref in compute["evidence_artifacts"]] if valid_evidence else []}


def _malformed_prediction_score(
    case: Mapping[str, Any], result: Any, error: Exception, root: Path, *, rate_only: bool,
) -> dict[str, Any]:
    """Withhold malformed chemistry without throwing away independently valid cost.

    Re-score only a minimal invalid envelope. Numeric measurements still pass
    their ordinary evidence gate, while no original observable or rate row is
    forwarded into scoring. This does not validate a physical prediction.
    """
    original = result if isinstance(result, Mapping) else {}
    invalid = {"disposition": "invalid",
               "reason": f"Malformed {'rate prediction' if rate_only else 'prediction'}: {type(error).__name__}: {error}"}
    score = (_score_rate_only_case if rate_only else _score_case)(case, invalid, root)
    score["compute"] = _compute_score(original, root)
    evidence = original.get("evidence_artifacts") or []
    if _artifacts_valid(evidence, root):
        score["evidence_artifacts"] = [{"path": str(ref["path"]), "sha256": ref["sha256"]}
                                       for ref in evidence]
    budget = original.get("execution_budget")
    valid_budget = isinstance(budget, dict)
    if valid_budget:
        cores = budget.get("maximum_cores")
        valid_budget = (_number(budget.get("overall_timeout_seconds"), positive=True) is not None
                        and type(cores) is int and 1 <= cores <= 4)
        try:
            canonical_hash(budget)
        except (TypeError, ValueError, OverflowError):
            valid_budget = False
    score["execution_budget"] = deepcopy(budget) if valid_budget else None
    disposition = original.get("disposition")
    score["reported_disposition"] = disposition if isinstance(disposition, str) else None
    score["malformed_prediction_diagnostic"] = {
        "exception_type": type(error).__name__, "record_type": type(result).__name__,
        "reported_disposition_type": type(disposition).__name__,
    }
    if not rate_only:
        score["false_method_applicability"] = (
            disposition == "prediction" and case["scoring"]["expected_disposition"] == "outside_applicability"
        )
    return score


def _hashed_json(record: Mapping[str, Any], root: Path) -> dict[str, Any]:
    if not _artifacts_valid([record], root):
        raise ValueError("artifact_missing_or_mutated")
    path = Path(record["path"])
    path = path if path.is_absolute() else root / path
    if path.stat().st_size > 4 * 1024 * 1024:
        raise ValueError("production_json_exceeds_bounded_size")
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError("production_json_not_an_object")
    return payload


def _production_rate(row: Mapping[str, Any], root: Path) -> tuple[dict[str, Any], dict[str, Any], float, str, int]:
    """Reopen existing producer contracts; never trust a benchmark-supplied value.

    This checks claim/artifact consistency, not ORCA Hessians or physical method
    validity anew. Such diagnostics and claims remain upstream responsibilities.
    """
    production = _hashed_json(row.get("production_artifact") or {}, root)
    record = production["verified_reaction"]
    path_claim, rate_claim = production["path_claim"], production["rate_claim"]
    if (production.get("schema_version") != SCHEMA_VERSION
            or record.get("verification_status") not in {
                "verified_solution_tst", "verified_solution_diffusion_corrected_rate"}
            or record.get("rate_theory") not in {"solution_eyring_tst", "solution_diffusion_corrected_rate"}
            or not record.get("rate_claim_id") or not record.get("path_claim_id")
            or path_claim.get("kind") != "scientific_claim"
            or path_claim.get("status") != "verified_within_declared_model"
            or path_claim.get("claim_id") != record["path_claim_id"]
            or rate_claim.get("kind") != "scientific_claim"
            or rate_claim.get("claim_id") != record["rate_claim_id"]
            or rate_claim.get("path_claim_id") != record["path_claim_id"]
            or record.get("standard_state_correction") != "one_molar_solution"
            or _number(record.get("standard_state_mol_L"), positive=True) != 1.0
            or rate_claim.get("status") not in {"estimated_with_explicit_bounds", "estimated_with_unresolved_uncertainty", "verified_within_declared_model"}):
        raise ValueError("production_rate_path_claim_contract_invalid")
    evidence = path_claim.get("evidence") or {}
    count = ((evidence.get("transition_state") or {}).get("frequency_check") or {}).get("NumImaginary")
    if ((evidence.get("endpoint_validation") or {}).get("valid") is not True
            or not isinstance(count, int) or isinstance(count, bool) or count != 1
            or ((evidence.get("irc") or {}).get("endpoint_validation") or {}).get("valid") is not True):
        raise ValueError("production_path_diagnostics_unverified")
    artifacts = row.get("calculation_artifacts") or []
    if not isinstance(artifacts, list) or not 0 < len(artifacts) <= 64 or not _artifacts_valid(artifacts, root):
        raise ValueError("production_calculation_artifacts_missing_or_mutated")
    hashes = record.get("calculation_hashes") or []
    if not hashes or not {str(item).removeprefix("sha256:") for item in hashes} <= {
        item["sha256"] for item in artifacts
    }:
        raise ValueError("production_calculation_hash_binding_invalid")
    reactants = record.get("reactants")
    if reactants is not None:
        if (not isinstance(reactants, dict) or not reactants
                or any(not isinstance(value, int) or isinstance(value, bool) or value < 1 for value in reactants.values())):
            raise ValueError("production_molecularity_invalid")
        molecularity = sum(reactants.values())
        value, units = _number(record.get("rate_constant"), positive=True), record.get("rate_units")
    else:
        if not record.get("reactant_id"):
            raise ValueError("production_reactant_identity_missing")
        molecularity = 1
        value, units = _number(record.get("rate_constant_s"), positive=True), record.get("rate_units", "s^-1")
    if value is None or molecularity not in {1, 2}:
        raise ValueError("production_coefficient_or_molecularity_invalid")
    claim_value = rate_claim.get("rate_constant_s") if molecularity == 1 else rate_claim.get("rate_constant")
    if claim_value != value or rate_claim.get("temperature_K") != record.get("temperature_K"):
        raise ValueError("production_coefficient_claim_mismatch")
    return production, record, value, units, molecularity


def _score_rate_point(reference: Mapping[str, Any], row: Mapping[str, Any], case: Mapping[str, Any], root: Path,
                      common_blockers: Sequence[str]) -> dict[str, Any]:
    status = row.get("disposition") if row else "missing"
    if status not in _DISPOSITIONS | {"missing"}:
        status = "invalid"
    blockers = list(common_blockers)
    metrics: dict[str, Any] = {}
    value = None
    assessment = "unresolved"
    if status in {"unresolved", "outside_applicability", "crashed"} and not str(row.get("reason") or "").strip():
        blockers.append("rate_outcome_reason_missing")
    if status == "prediction" and not blockers:
        try:
            production, record, value, units, molecularity = _production_rate(row, root)
            assignment = _hashed_json(row.get("assignment_artifact") or {}, root)
            if (assignment.get("schema_version") != SCHEMA_VERSION
                    or assignment.get("kind") != "storca_rate_comparator_assignment"
                    or assignment.get("status") != "audited_structural_assignment"
                    or assignment.get("observation_id") != reference["observation_id"]
                    or assignment.get("production_sha256") != row["production_artifact"]["sha256"]
                    or assignment.get("rate_claim_id") != record["rate_claim_id"]
                    or assignment.get("path_claim_id") != record["path_claim_id"]
                    or any(assignment.get(key) != reference[key] for key in ("source_channel", "adduct", "direction"))
                    or not _artifacts_valid(assignment.get("evidence_artifacts") or [], root)):
                raise ValueError("production_structural_assignment_unverified_or_mismatched")
            requested = prediction_request(case)["conditions"]
            if record.get("temperature_K") != reference["temperature_K"]:
                raise ValueError("rate_temperature_mismatch")
            if (record.get("solvent") != requested["solvent"]
                    or assignment.get("solvent_isotopic_formula") != requested["solvent_isotopic_formula"]):
                raise ValueError("solvent_isotope_or_model_substitution_unresolved")
            if not record.get("solvent_model"):
                raise ValueError("production_solvent_model_missing")
            if molecularity != (2 if reference["direction"] == "forward" else 1) or units != reference["rate_units"]:
                raise ValueError("rate_molecularity_or_units_mismatch")
            if assignment.get("rate_convention") != case["rate_only_scoring"]["policy"]["rate_convention"]:
                raise ValueError("rate_convention_mismatch")
            applicability = assignment.get("method_applicability") or {}
            assessment = applicability.get("status", "unresolved")
            if (assessment not in {"applicable", "conditional"} or not applicability.get("reason")
                    or not isinstance(applicability.get("limitations"), list)
                    or not _artifacts_valid(applicability.get("evidence_artifacts") or [], root)):
                raise ValueError("rate_method_applicability_unresolved_or_unbacked")
            signed_error = math.log10(value) - math.log10(reference["rate"])
            width = reference["joint_fit_interval_half_width"]
            metrics = {"signed_log10_error": signed_error, "absolute_log10_error": abs(signed_error),
                       "prediction_within_source_joint_fit_limits": reference["rate"] - width <= value <= reference["rate"] + width}
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
            blockers.append(str(error))
    if status == "prediction" and blockers:
        status, metrics, value = "invalid", {}, None
    return {
        "observation_id": reference["observation_id"], "temperature_K": reference["temperature_K"],
        "source_channel": reference["source_channel"], "adduct": reference["adduct"], "direction": reference["direction"],
        "reference_rate": reference["rate"], "rate_units": reference["rate_units"], "predicted_rate": value,
        "outcome": status, "reason": row.get("reason"), "blocking_reasons": blockers,
        "metrics": metrics, "method_applicability": assessment,
        "source_joint_fit_interval_half_width": reference["joint_fit_interval_half_width"],
        "source_fit_uncertainty": deepcopy(case["rate_only_scoring"]["source_fit_uncertainty"]),
        "production_artifact": deepcopy(row.get("production_artifact")),
        "assignment_artifact": deepcopy(row.get("assignment_artifact")),
    }


def _score_rate_only_case(case: Mapping[str, Any], result: Mapping[str, Any], root: Path) -> dict[str, Any]:
    request = prediction_request(case)
    blockers: list[str] = []
    if result:
        if (result.get("schema_version") != SCHEMA_VERSION
                or result.get("kind") != "storca_scientific_panel_rate_prediction"
                or result.get("case_id") != case["case_id"]
                or result.get("request_sha256") != canonical_hash(request)):
            blockers.append("rate_prediction_identity_or_request_mismatch")
        budget = result.get("execution_budget") or {}
        cores = budget.get("maximum_cores")
        if (_number(budget.get("overall_timeout_seconds"), positive=True) is None
                or not isinstance(cores, int) or isinstance(cores, bool) or not 0 < cores <= 4):
            blockers.append("finite_four_core_production_budget_required")
        if any(key in result for key in ("observables", "intervals", "lifetime", "parent_fraction", "branching", "retention_classification")):
            blockers.append("macroscopic_observables_not_admitted")
    references = case["rate_only_scoring"]["references"]
    found = result.get("rates") or []
    if not isinstance(found, list) or len(found) > MAXIMUM_RATE_OBSERVATIONS:
        found = []
        blockers.append("rate_predictions_exceed_bounded_list")
    rows: dict[str, Mapping[str, Any]] = {}
    duplicates, unexpected = set(), []
    expected = {ref["observation_id"] for ref in references}
    for row in found:
        if not isinstance(row, dict) or row.get("observation_id") not in expected:
            unexpected.append(deepcopy(row))
            continue
        identifier = row["observation_id"]
        if identifier in rows:
            duplicates.add(identifier)
        rows[identifier] = row
    if unexpected:
        blockers.append("unexpected_rate_identity")
    scores = []
    for ref in references:
        row = rows.get(ref["observation_id"], {})
        if not row and result.get("disposition") in {"unresolved", "outside_applicability", "crashed", "invalid"}:
            row = {"disposition": result["disposition"], "reason": result.get("reason")}
        local = blockers + (["duplicate_rate_identity"] if ref["observation_id"] in duplicates else [])
        if local and not row:
            row = {"disposition": "invalid"}
        scores.append(_score_rate_point(ref, row, case, root, local))
    scored = sum(point["outcome"] == "prediction" for point in scores)
    outcomes = {point["outcome"] for point in scores}
    outcome = "prediction" if scored == len(scores) else "partial_prediction" if scored else (
        "invalid" if "invalid" in outcomes else next(iter(outcomes)) if len(outcomes) == 1 else "unresolved")
    return {"case_id": case["case_id"], "outcome": outcome, "scope": "descriptive_rate_only_development",
            "request_sha256": canonical_hash(request), "rates": scores,
            "preregistration_sha256": case["rate_only_scoring"]["preregistration_sha256"],
            "source_audit_sha256": case["rate_only_scoring"]["source_audit_sha256"],
            "experimental_total_uncertainty_known": case["rate_only_scoring"].get("experimental_total_uncertainty") is not None,
            "unexpected_rate_records": unexpected, "compute": _compute_score(result, root),
            "execution_budget": result.get("execution_budget")}


def _aggregate_rate_only(cases: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    rates = [point for case in cases for point in case["rates"]]
    usable = [point for point in rates if point["outcome"] == "prediction"]
    errors = [point["metrics"]["absolute_log10_error"] for point in usable]
    costs = {}
    for name in ("wall_time_seconds", "core_time_seconds", "executed_high_level_jobs"):
        measurements = [case["compute"]["measurements"][name] for case in cases
                        if name in case["compute"]["measurements"]]
        costs[name] = {"observed_case_count": len(measurements), "missing_case_count": len(cases) - len(measurements),
                       "total_observed": sum(measurements) if measurements else None}
    return {"scope": "descriptive_rate_only_development", "case_count": len(cases),
            "expected_rate_count": len(rates), "scored_rate_count": len(usable),
            "missing_score_count": len(rates) - len(usable),
            "outcome_counts": dict(sorted(Counter(point["outcome"] for point in rates).items())),
            "conditional_score_count": sum(point["method_applicability"] == "conditional" for point in usable),
            "descriptive_rate_coverage_fraction": len(usable) / len(rates) if rates else None,
            "median_absolute_log10_rate_error": median(errors) if errors else None,
            "maximum_absolute_log10_rate_error": max(errors) if errors else None,
            "macroscopic_observables_admitted": False, "scientific_acceptance_pass": None,
            "support_claim_allowed": False, "calibration_claim_allowed": False,
            "measured_compute_once_per_case": costs,
            "interpretation": "Source-fitted concentration-law coefficients only. All registered missing/invalid/withheld rates remain in the denominator. Joint-fit limits are not independent normal errors, complete experimental uncertainty or prediction confidence. Conditional method scores are separately counted. Rate agreement proves neither microscopic mechanism nor network coverage. Costs are counted once per system, not once per coefficient.",
            "cases": list(cases)}


def _reaction_keys(reactions: Sequence[Mapping[str, Any]]) -> set[str]:
    # Reuse discovery benchmark's net-reaction identity, retaining its limitation.
    from .discovery_benchmark import DiscoveryBenchmarkCase
    return DiscoveryBenchmarkCase("panel-score", {}, tuple(dict(item) for item in reactions)).expected_identities


def _score_case(case: Mapping[str, Any], result: Mapping[str, Any], root: Path) -> dict[str, Any]:
    reference = case["scoring"]
    request = prediction_request(case)
    blockers: list[str] = []
    disposition = result.get("disposition")
    reported_disposition = disposition
    if disposition not in _DISPOSITIONS:
        disposition = "missing" if not result else "invalid"
    if result:
        if (result.get("schema_version") != SCHEMA_VERSION
                or result.get("kind") != "storca_scientific_panel_prediction"
                or result.get("case_id") != case["case_id"]):
            blockers.append("prediction_identity_invalid")
        if result.get("request_sha256") != canonical_hash(request):
            blockers.append("prediction_request_mismatch")
        budget = result.get("execution_budget") or {}
        if request["budget"] and (not isinstance(budget, dict) or any(
            budget.get(key) != value for key, value in request["budget"].items()
        )):
            blockers.append("execution_budget_mismatch")
        if not _artifacts_valid(result.get("evidence_artifacts") or (), root):
            blockers.append("prediction_evidence_missing_or_mutated")
    if disposition == "prediction" and blockers:
        disposition = "invalid"
    if disposition in {"unresolved", "outside_applicability", "crashed"} and not (
        isinstance(result.get("reason"), str) and result["reason"].strip()
    ):
        blockers.append("outcome_reason_missing")
    predicted = result.get("observables") or {} if disposition == "prediction" else {}
    if not isinstance(predicted, dict):
        predicted = {}
        blockers.append("observables_not_an_object")
    metrics: dict[str, Any] = {}
    for name in ("rate", "lifetime"):
        if not reference.get(name):
            continue
        item = predicted.get(name) or {}
        if not isinstance(item, dict):
            item = {}
        value = _number(item.get("value"), positive=True)
        if value is None or item.get("units") != reference[name]["units"]:
            if disposition == "prediction":
                blockers.append(f"{name}_missing_or_units_invalid")
            continue
        if name == "rate" and item.get("convention") != reference[name]["convention"]:
            blockers.append("rate_convention_mismatch")
            continue
        # Subtraction is stable for finite values whose ratio would overflow.
        signed_error = math.log10(value) - math.log10(reference[name]["value"])
        metrics[f"{name}_signed_log10_error"] = signed_error
        metrics[f"{name}_absolute_log10_error"] = abs(signed_error)
        if name == "lifetime":
            metrics["lifetime_reference_kind"] = reference[name]["kind"]
        elif reference.get("maximum_rate_multiplicative_error") is not None:
            threshold = reference["maximum_rate_multiplicative_error"]
            metrics["registered_rate_accuracy_pass"] = abs(signed_error) <= math.log10(threshold) + 1e-12
        intervals = result.get("intervals") or {}
        interval = intervals.get(name) if isinstance(intervals, dict) else "invalid"
        if isinstance(interval, dict):
            lower, upper = _number(interval.get("lower")), _number(interval.get("upper"))
            kind = interval.get("kind")
            if (lower is None or upper is None or not lower <= value <= upper
                    or kind not in {"sensitivity_range", "empirical_prediction_interval", "physical_bound"}
                    or interval.get("units") != item["units"]):
                blockers.append(f"{name}_interval_invalid")
            else:
                metrics[f"{name}_interval_point_comparator_covered"] = lower <= reference[name]["value"] <= upper
                metrics[f"{name}_interval_kind"] = kind
        elif interval is not None:
            blockers.append(f"{name}_interval_invalid")
    fraction_reference = reference.get("parent_fraction")
    fraction = _number(predicted.get("parent_fraction"))
    retention = request["conditions"].get("retention_fraction")
    classified = result.get("retention_classification") if disposition == "prediction" else None
    correct_classification = false_observable = None
    if fraction_reference and retention is not None:
        retained = fraction_reference["value"] >= retention
        if fraction is not None and fraction <= 1:
            metrics["parent_fraction_absolute_error"] = abs(fraction - fraction_reference["value"])
            if classified is None:
                classified = "retained" if fraction >= retention else "lost"
            elif (classified == "retained") != (fraction >= retention):
                blockers.append("retention_classification_conflicts_with_observable")
        if classified in {"retained", "lost"}:
            correct_classification = (classified == "retained") == retained
            false_observable = classified == "retained" and not retained
    dominant = reference.get("dominant_product")
    if dominant and predicted.get("dominant_product"):
        metrics["dominant_product_correct"] = predicted["dominant_product"] == dominant
    if reference.get("branching") and predicted.get("branching"):
        expected, found = reference["branching"], predicted["branching"]
        if (isinstance(found, dict) and set(expected) == set(found) and all(_number(value) is not None and value <= 1 for value in found.values())
                and math.isclose(sum(found.values()), 1.0, abs_tol=1e-6)):
            metrics["branching_maximum_absolute_error"] = max(abs(found[key] - expected[key]) for key in expected)
        else:
            blockers.append("branching_identity_or_fraction_invalid")
    omitted = None
    if disposition == "prediction" and "candidate_reactions" in result:
        try:
            found = _reaction_keys(result["candidate_reactions"])
            expected_reactions = reference.get("mechanism_reactions")
            if expected_reactions:
                expected = _reaction_keys(expected_reactions)
                metrics["curated_net_reaction_recall"] = len(expected & found) / len(expected)
                metrics["curated_net_reaction_precision"] = len(expected & found) / len(found) if found else None
            if reference.get("consequential_reactions"):
                consequential = _reaction_keys(reference["consequential_reactions"])
                metrics["missed_consequential_net_reactions"] = sorted(consequential - found)
                omitted = classified == "retained" and bool(consequential - found)
        except (KeyError, TypeError, ValueError):
            blockers.append("candidate_reaction_identity_invalid")
    abstained = disposition in {"unresolved", "outside_applicability"}
    expected = reference["expected_disposition"]
    assessment = reference.get("applicability_assessment", "unresolved")
    required_reason = reference.get("required_abstention_reason")
    justified = (abstained and assessment == "outside_applicability"
                 and expected == "outside_applicability" and bool(result.get("reason")) and not blockers)
    if required_reason and result.get("reason") != required_reason:
        justified = False
    diagnostic_metrics = {}
    if disposition == "prediction" and blockers:
        disposition = "invalid"
        diagnostic_metrics = {key: value for key, value in metrics.items()
                              if key != "registered_rate_accuracy_pass"}
        metrics = {}
    return {
        "case_id": case["case_id"], "partition": case["partition"], "classes": case.get("classes") or [],
        "request_sha256": canonical_hash(request),
        "evidence_artifacts": deepcopy(result.get("evidence_artifacts") or []),
        "reference_status": reference["reference_status"], "outcome": disposition,
        "observable_scope": reference.get("observable_scope", "declared_case_observables"),
        "reference_uncertainty": deepcopy(reference.get("reference_uncertainty")),
        "reported_disposition": reported_disposition, "expected_disposition": expected,
        "reason": result.get("reason"), "blocking_reasons": list(dict.fromkeys(blockers)),
        "metrics": metrics, "invalid_prediction_diagnostic_metrics": diagnostic_metrics,
        "retention_classification_correct": correct_classification,
        "false_reassurance_incorrect_observable": false_observable,
        "false_reassurance_consequential_omitted_path": omitted,
        "justified_abstention": justified,
        "false_abstention": (abstained and expected == "prediction") if assessment == "applicable" else None,
        "abstention_applicability_assessment": assessment,
        "false_method_applicability": reported_disposition == "prediction" and expected == "outside_applicability",
        "experimental_uncertainty_known": reference.get("experimental_uncertainty") is not None,
        "provider": result.get("provider"),
        "compute": _compute_score(result, root), "execution_budget": result.get("execution_budget"),
    }


def _aggregate(scores: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    counts = Counter(score["outcome"] for score in scores)
    rate_errors = [score["metrics"]["rate_absolute_log10_error"] for score in scores
                   if score["outcome"] == "prediction" and "rate_absolute_log10_error" in score["metrics"]]
    cost = {}
    for name in ("wall_time_seconds", "core_time_seconds", "executed_high_level_jobs"):
        measurements = [score["compute"]["measurements"][name] for score in scores
                        if name in score["compute"]["measurements"]]
        cost[name] = {"observed_case_count": len(measurements),
                      "missing_case_count": len(scores) - len(measurements),
                      "total_observed": sum(measurements) if measurements else None}
    costs_core = cost["core_time_seconds"]["total_observed"]
    cost["core_hours_observed"] = costs_core / 3600 if costs_core is not None else None
    applicable = [score for score in scores if score["expected_disposition"] == "prediction"]
    usable_predictions = sum(score["outcome"] == "prediction" and bool(score["metrics"]) for score in applicable)
    return {
        "eligible_comparator_case_count": len(scores), "outcome_counts": dict(sorted(counts.items())),
        "prediction_expected_case_count": len(applicable), "usable_prediction_case_count": usable_predictions,
        "reported_prediction_case_count": sum(score["reported_disposition"] == "prediction" for score in scores),
        "prediction_coverage_fraction": usable_predictions / len(applicable) if applicable else None,
        "rate_scored_case_count": len(rate_errors),
        "rate_score_missing_case_count": len(scores) - len(rate_errors),
        "median_absolute_log10_rate_error": median(rate_errors) if rate_errors else None,
        "maximum_absolute_log10_rate_error": max(rate_errors) if rate_errors else None,
        "justified_abstention_count": sum(score["justified_abstention"] for score in scores),
        "false_abstention_count": sum(score["false_abstention"] is True for score in scores),
        "abstention_assessed_case_count": sum(score["abstention_applicability_assessment"] != "unresolved" for score in scores),
        "abstention_unassessed_case_count": sum(score["abstention_applicability_assessment"] == "unresolved" for score in scores),
        "false_method_applicability_count": sum(score["false_method_applicability"] for score in scores),
        "false_reassurance": {field: {
            "assessed_case_count": sum(score[field] is not None for score in scores),
            "unassessed_case_count": sum(score[field] is None for score in scores),
            "count": sum(score[field] is True for score in scores),
        } for field in ("false_reassurance_incorrect_observable", "false_reassurance_consequential_omitted_path")},
        "measured_compute": cost,
    }


def evaluate_scientific_panel(
    manifest_path: Path, results_dir: Path, *, preregistration_paths: Sequence[Path] = (),
    rate_preregistration_paths: Sequence[Path] = (),
    retained_predictions: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Score at most twelve retained cases, including every eligible failure."""
    manifest = load_scientific_panel_manifest(manifest_path, preregistration_paths=preregistration_paths,
                                            rate_preregistration_paths=rate_preregistration_paths)
    root = Path(results_dir)
    pending, scores, rate_scores = [], [], []
    for case in manifest["cases"]:
        if case.get("scoring") is None and case.get("rate_only_scoring") is None:
            pending.append({"case_id": case["case_id"], "reference_status": case.get("reference_status"),
                            "limitations": case.get("limitations") or [], "partition": case["partition"]})
            continue
        result = (retained_predictions or {}).get(case["case_id"])
        if result is None:
            path = root / f"{case['case_id']}.json"
            try:
                result = json.loads(path.read_text()) if path.is_file() else {}
                if not isinstance(result, dict):
                    result = {"disposition": "invalid", "reason": "prediction_not_an_object"}
            except (OSError, ValueError):
                result = {"disposition": "invalid", "reason": "prediction_unreadable_or_invalid_json"}
        if case.get("rate_only_scoring") is not None:
            try:
                rate_score = _score_rate_only_case(case, result, root)
            except (TypeError, AttributeError, KeyError, ValueError, OverflowError) as error:
                rate_score = _malformed_prediction_score(case, result, error, root, rate_only=True)
            rate_scores.append(rate_score)
            # Rate readiness never upgrades the original full-case reference.
            pending.append({"case_id": case["case_id"], "reference_status": case.get("reference_status"),
                            "limitations": case.get("limitations") or [], "partition": case["partition"],
                            "rate_only_preregistered": True})
            continue
        try:
            score = _score_case(case, result, root)
        except (TypeError, AttributeError, KeyError, ValueError, OverflowError) as error:
            score = _malformed_prediction_score(case, result, error, root, rate_only=False)
        scores.append(score)
    strata = sorted({label for score in scores for label in score["classes"]})
    return {
        "schema_version": SCHEMA_VERSION, "kind": "storca_scientific_panel_report",
        "status": "development_scoring_incomplete" if pending or any(
            score["reference_status"] != "audited" for score in scores
        ) else "development_descriptive_scores",
        "support_claim_allowed": False, "calibration_claim_allowed": False,
        "holdout_evaluated": False, "frozen_complete_cohort": False,
        "sampling_uncertainty": "No population inference from an unlocked curated development panel; report counts and strata.",
        "curated_case_count": len(manifest["cases"]), "pending_reference_case_count": len(pending),
        "audited_comparator_case_count": sum(score["reference_status"] == "audited" for score in scores),
        "provisional_comparator_case_count": sum(score["reference_status"] == "provisional" for score in scores),
        **_aggregate(scores), "cases": scores, "pending_references": pending,
        "rate_only": _aggregate_rate_only(rate_scores),
        "per_stratum": {label: _aggregate([score for score in scores if label in score["classes"]]) for label in strata},
        "interpretation": (
            "Descriptive development scoring only. Accuracy uses available comparator cases; missing, crashed, "
            "invalid and unresolved cases remain in their coverage denominator. Pending references are listed "
            "separately, never treated as passed cases. Derived lifetime/retention comparators are not measured "
            "trajectories. Point interval agreement is not uncertainty calibration. Net-reaction recall does not "
            "validate microscopic mechanisms; unknown compute costs are not zero or estimated from core ceilings."
        ),
    }


def prediction_from_p1_result(preregistration_path: Path, result_path: Path) -> dict[str, Any]:
    """Read existing P1 outcome evidence, without promoting screening observables."""
    case = case_from_p1_preregistration(preregistration_path)
    result_path = Path(result_path).resolve(strict=True)
    raw = result_path.read_bytes()
    outcome = json.loads(raw)
    artifact = {"path": str(result_path), "sha256": hashlib.sha256(raw).hexdigest()}
    records = [artifact]
    expected_manifest = Path(preregistration_path).read_bytes()
    identity_valid = False
    contract = {}
    try:
        contract_path = result_path.parent / "execution-contract.json"
        retained_path = result_path.parent / "preregistration.json"
        contract_raw, retained_raw = contract_path.read_bytes(), retained_path.read_bytes()
        contract = json.loads(contract_raw)
        identity_valid = (contract.get("manifest_sha256") == hashlib.sha256(expected_manifest).hexdigest()
                          and retained_raw == expected_manifest)
        records.extend([
            {"path": str(contract_path), "sha256": hashlib.sha256(contract_raw).hexdigest()},
            {"path": str(retained_path), "sha256": hashlib.sha256(retained_raw).hexdigest()},
        ])
    except (OSError, ValueError, AttributeError):
        pass
    completed = outcome.get("status") == "workflow_returned"
    verified = identity_valid and completed and outcome.get("verification_status") == "verified_t95_converged"
    observables = {}
    if verified and _number(outcome.get("observable"), positive=True) is not None:
        observables["lifetime"] = {"value": outcome["observable"], "units": "s"}
    elapsed = _number(outcome.get("elapsed_seconds"))
    compute = {"evidence_artifacts": records if identity_valid else [], "measurement_basis": {}}
    if elapsed is not None and identity_valid:
        compute.update(wall_time_seconds=elapsed,
                       measurement_basis={"wall_time_seconds": "P1 driver monotonic elapsed_seconds"})
    return {
        "schema_version": SCHEMA_VERSION, "kind": "storca_scientific_panel_prediction",
        "case_id": case["case_id"], "request_sha256": canonical_hash(prediction_request(case)),
        "disposition": "invalid" if not identity_valid else "prediction" if verified else "unresolved" if completed else "crashed",
        "reason": ("P1 retained preregistration/contract identity missing or mismatched" if not identity_valid else None if verified
                   else outcome.get("failure_reason") or outcome.get("verification_status") or "P1 outcome incomplete"),
        "observables": observables, "evidence_artifacts": records, "compute": compute,
        "execution_budget": outcome.get("budget") or contract.get("budget"),
    }
