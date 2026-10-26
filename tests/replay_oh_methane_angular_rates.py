"""Offline refusal/freeze of unsupported OH/CH4 angular-rate corrections.

This consumer reopens explicit source bindings; it does not dispatch programs,
evaluate a supplied partition function, or recertify retained conditional rates.
The periodic-profile experiment currently has no reviewed quantitative angular
Hamiltonian. Success flags and arbitrary Q(T) tables cannot supply one.
"""
from __future__ import annotations

import argparse
from enum import Enum
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
_spec = importlib.util.spec_from_file_location(
    "angular_rate_retained_binding_reader",
    Path(__file__).with_name("replay_oh_methane_electronic_states.py"),
)
_retained = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_retained)

PROFILE_KIND = "bounded_OH_CH4_periodic_profile"
CONSUMER_KIND = "bounded_OH_CH4_angular_rate_admission"


class RateAdmissionStatus(str, Enum):
    PROFILE_UNRESOLVED = "profile_unresolved"
    PARTITION_UNSUPPORTED = "profile_complete_partition_unsupported"


_REQUIRED_EVIDENCE = (
    ("reviewed_angular_overlay", "No reviewed quantitative angular-overlay adapter is implemented."),
    ("global_domain_and_branch", "A complete, single-branch physical integration domain is not certified by this consumer."),
    ("kinetic_measure", "A derived relaxed-surface Hamiltonian/kinetic measure is required; the inverse-constraint metric is insufficient."),
    ("rotation_and_reaction_normal", "Overall-rotation coupling and exclusion of the reactive direction need an explicit physical derivation."),
    ("harmonic_and_rotor_limits", "Harmonic/free-rotor limits and coordinate-reparameterization invariance need replayable model tests."),
    ("transverse_bath", "Sparse selected Hessians do not establish transverse stability, ZPE and separability throughout the domain."),
    ("counting", "External symmetry, channel, spectator-permutation, basin and optical counting must be reconciled together."),
    ("energy_zero_and_ZPE", "Exactly-once angular/ZPE replacement and an explicit common energy zero have not been implemented."),
    ("forward_reverse_and_tunneling", "Both barriers, transmission assumptions, energy-zero invariance and detailed balance need a validated overlay replay."),
)
_CLAIM_FLAGS = (
    "angular_partition_admitted", "rate_claim_allowed", "lifetime_claim_allowed",
    "network_reinjection_allowed", "production_defaults_changed",
    "physical_rate_accuracy_claim_allowed", "calibrated_uncertainty_claim_allowed",
)


def _binding(record):
    if (not isinstance(record, dict) or not isinstance(record.get("path"), str)
            or not Path(record["path"]).is_absolute()
            or not isinstance(record.get("sha256"), str)
            or re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is None):
        raise ValueError("An absolute, SHA256-bound source artifact is required")
    if "mtime_ns" in record and (type(record["mtime_ns"]) is not int or record["mtime_ns"] < 0):
        raise ValueError("Invalid source modification-time binding")
    return _retained.reopen(record)


def _source(value):
    if isinstance(value, dict):
        if not {"path", "sha256"} <= value.keys():
            raise ValueError("Use an immutable artifact binding, not an unbound result dictionary")
        path = _binding(value)
    elif isinstance(value, (str, Path)):
        path = Path(value).resolve(strict=True)
        if path.is_dir():
            path = path / "result.json"
    else:
        raise ValueError("Use a profile result path or an immutable artifact binding, not an unbound result dictionary")
    return _retained.artifact(path)


class _Sources:
    """Reopen current receipts; hash historical JSON without reinterpreting it.

    The campaign's predecessor replay owns historical semantic/revision policy.
    No directory is scanned: only explicitly bound files can be opened.
    """

    def __init__(self, current_directory):
        self.current_directory = current_directory
        self.bindings = []
        self.documents = {}
        self.hashes = {}
        self.opaque_json = set()

    def add(self, record, *, expand=False, opaque=False):
        path = _binding(record)
        previous = self.hashes.setdefault(str(path), record["sha256"])
        if previous != record["sha256"]:
            raise ValueError("Conflicting source hashes for the same artifact")
        self.bindings.append(dict(record))
        if path.suffix.lower() == ".json":
            current_receipt = self.current_directory in path.parents
            if not opaque and (expand or current_receipt):
                self.opaque_json.discard(str(path))
                if str(path) not in self.documents:
                    document = _retained.read_json(path)
                    self.documents[str(path)] = document
                    self.walk(document, self._reference_pdfs(document, path))
            elif str(path) not in self.documents:
                self.opaque_json.add(str(path))
        return path

    @staticmethod
    def _reference_pdfs(document, audit_path):
        """Interpret only the retained reference-audit schema's local PDFs.

        This does not relax the absolute-artifact contract elsewhere. Both the
        containing audit and the resolved PDF remain independently hash bound.
        """
        if (not isinstance(document, dict) or type(document.get("schema_version")) is not int
                or document["schema_version"] != 1
                or document.get("kind") != "storca_scientific_panel_reference_audit"):
            return {}
        entries = document.get("sources")
        if not isinstance(entries, list):
            raise ValueError("Reference-audit sources must be a list")
        interpreted = {}
        for entry in entries:
            if (not isinstance(entry, dict) or not isinstance(entry.get("path"), str)
                    or not {"path", "sha256"} <= entry.keys()):
                raise ValueError("Reference-audit source lacks its artifact binding")
            relative = Path(entry["path"])
            if relative.is_absolute():
                continue
            if (relative.name != entry["path"] or relative.suffix.lower() != ".pdf"
                    or any(not isinstance(entry.get(key), str) or not entry[key].strip()
                           for key in ("source_id", "source_kind"))):
                raise ValueError("Only simple local PDF filenames are allowed in reference-audit sources")
            resolved = (audit_path.parent / relative).resolve(strict=True)
            if resolved.parent != audit_path.parent:
                raise ValueError("Reference-audit PDF escapes its bound audit directory")
            interpreted[id(entry)] = {**entry, "path": str(resolved)}
        return interpreted

    def walk(self, value, reference_pdfs):
        if isinstance(value, dict):
            if "path" in value and "sha256" in value:
                self.add(reference_pdfs.get(id(value), value))
            opaque_execution = None
            if value.get("kind") == "owned_normal_native_adoption" and value.get("json_readable") is False:
                outcome = value.get("native_outcome")
                if (not isinstance(outcome, dict) or outcome.get("status") != "normal"
                        or type(outcome.get("returncode")) is not int or outcome["returncode"] != 0
                        or type(outcome.get("normal_termination_markers")) is not int
                        or outcome["normal_termination_markers"] != 1
                        or value.get("old_execution_rewritten") is not False
                        or value.get("rate_claim_allowed") is not False):
                    raise ValueError("Unreadable execution adoption requires diagnostic-only owned-normal native evidence")
                for name in ("input", "output", "lease"):
                    _binding(outcome.get(name))
                opaque_execution = value.get("original_execution")
                # Preserve and verify damaged wrapper bytes. The controller,
                # not this rate-refusal consumer, validates owned native state.
                self.add(opaque_execution, opaque=True)
            for key, item in value.items():
                if key == "original_execution" and opaque_execution is not None:
                    continue
                explicit_inventory = (key in {"source_inventory", "retained_source_inventory"}
                                      or (key == "inventory" and isinstance(item, dict)
                                          and ("prefix" in value or any(isinstance(record, dict)
                                               and "path" in record and "sha256" in record
                                               for record in item.values()))))
                if explicit_inventory and isinstance(item, dict):
                    for source_path, source_record in item.items():
                        if not isinstance(source_record, dict) or source_record.get("path") != source_path:
                            raise ValueError("Inventory key differs from its bound source path")
                self.walk(item, reference_pdfs)
        elif isinstance(value, list):
            for item in value:
                self.walk(item, reference_pdfs)

    def finish(self):
        # Reopen after parsing too: a changed source must not supply a receipt.
        for record in self.bindings:
            _binding(record)


def _point_summary(points):
    if not isinstance(points, list):
        raise ValueError("Profile points must be a list")
    counts = {"constrained_optimizer_candidate": 0, "derivative_verified_constrained": 0,
              "other_or_unresolved": 0, "extra_full_negative_modes": 0}
    refused = []
    for index, point in enumerate(points):
        if not isinstance(point, dict):
            raise ValueError("Each profile point must be an object")
        offset = point.get("offset_deg")
        if type(offset) not in (int, float) or not math.isfinite(offset):
            raise ValueError("A finite mapped angular offset is required")
        status = point.get("status")
        if not isinstance(status, str):
            raise ValueError("A profile point must declare its evidence status")
        counts[status if status in {"constrained_optimizer_candidate", "derivative_verified_constrained"}
               else "other_or_unresolved"] += 1
        frequency = point.get("frequency", {})
        if not isinstance(frequency, dict):
            raise ValueError("Profile frequency evidence must be an object")
        values = frequency.get("raw_full_frequencies_cm_inv")
        if values is not None:
            if (not isinstance(values, list) or len(values) != 21
                    or any(type(v) not in (int, float) or not math.isfinite(v) for v in values)):
                raise ValueError("Full-frequency inventory must contain 21 finite numbers")
            physical = values[6:]
            negative = sum(value < 0 for value in physical)
            if negative > 1:
                counts["extra_full_negative_modes"] += 1
                ordinary = (point.get("kind") in {"transition_state", "ordinary_transition_state"}
                            or point.get("ordinary_first_order_TS") is True
                            or point.get("ordinary_transition_state") is True)
                if ordinary:
                    refused.append({"code": "non_first_order_point_claimed_as_ordinary_TS",
                                    "point_index": index, "offset_deg": offset,
                                    "full_negative_count": negative})
    return counts, refused


def assess(profile_result, profile_contract=None):
    """Return an immutable-source-bound refusal, never corrected rate constants.

    `profile_result` is a result path/directory or {path, sha256} binding. An
    optional explicit contract must match the result's contract binding. No
    success boolean or external Q(T) payload can enable an unimplemented model.
    """
    result_binding = _source(profile_result)
    sources = _Sources(Path(result_binding["path"]).parent)
    result_path = sources.add(result_binding, expand=True)
    result = sources.documents.get(str(result_path))
    if (not isinstance(result, dict) or type(result.get("schema_version")) is not int
            or result["schema_version"] != 1 or result.get("kind") != PROFILE_KIND):
        raise ValueError("Expected the schema-1 bounded OH/CH4 periodic-profile result")
    contract_record = result.get("contract")
    contract_path = sources.add(contract_record, expand=True)
    if profile_contract is not None and _source(profile_contract) != _retained.artifact(contract_path):
        raise ValueError("Explicit profile contract differs from the result's bound contract")
    contract = sources.documents.get(str(contract_path))
    if (not isinstance(contract, dict) or type(contract.get("schema_version")) is not int
            or contract["schema_version"] != 1 or contract.get("kind") != PROFILE_KIND):
        raise ValueError("Expected the matching schema-1 periodic-profile contract")
    if not isinstance(contract.get("source_inventory"), dict):
        raise ValueError("Periodic-profile contract requires an explicit source inventory")
    sources.add(contract.get("source"))
    policy = contract.get("partition_policy")
    if (not isinstance(policy, dict) or policy.get("mode") != "diagnostic_only"
            or policy.get("supported_quantitative_hamiltonian") is not False):
        raise ValueError("Unsupported profile partition policy: only the explicit diagnostic-only protocol is implemented")
    source_status = result.get("status")
    if not isinstance(source_status, str) or source_status not in {item.value for item in RateAdmissionStatus}:
        raise ValueError("Unsupported periodic-profile disposition")
    counts, refused = _point_summary(result.get("points"))
    for location, document in (("result", result), ("contract", contract)):
        for flag in _CLAIM_FLAGS:
            if flag in document and document[flag] is not False:
                refused.append({"code": "unsupported_admission_claim", "location": location, "field": flag})

    retained = {}
    for name in ("baseline_recipe", "previous_conditional_branches"):
        records = [document[name] for document in (contract, result) if name in document]
        if records:
            paths = [sources.add(record) for record in records]
            if len({str(path) for path in paths}) != 1:
                raise ValueError("Conflicting retained " + name + " bindings")
            retained[name] = {"source": _retained.artifact(paths[0]),
                              "disposition": "preserved_source_only_not_recomputed_or_recertified"}

    missing = [{"code": code, "reason": reason} for code, reason in _REQUIRED_EVIDENCE]
    if counts["constrained_optimizer_candidate"]:
        missing.append({"code": "optimizer_only_nodes_not_independently_verified",
                        "reason": "A read-in or quasi-Newton optimization Hessian is not a final independent saddle-index certificate.",
                        "count": counts["constrained_optimizer_candidate"]})
    if source_status == RateAdmissionStatus.PROFILE_UNRESOLVED.value:
        missing.insert(0, {"code": "source_profile_unresolved",
                           "reason": result.get("failure_reason") or "The bounded source profile remains unresolved."})
    receipt = {
        "schema_version": 1, "kind": CONSUMER_KIND, "status": source_status,
        "source_profile_result": result_binding, "source_profile_contract": _retained.artifact(contract_path),
        "consumer": _retained.artifact(Path(__file__)),
        "binding_reader": _retained.artifact(Path(_retained.__file__)),
        "point_evidence_summary": counts, "missing_evidence": missing, "refused_claims": refused,
        "retained_conditional_evidence": retained,
        "sources_reopened": len(sources.hashes),
        "current_receipt_directory": str(sources.current_directory),
        "opaque_provenance_JSON": sorted(sources.opaque_json),
        "source_integrity_scope": "Current profile result/contract and explicitly bound JSON receipts beneath the current profile directory are traversed. Referenced historical JSON bytes and all explicit current artifact hashes/optional mtimes are checked; historical JSON contents, legacy list inventories and scientific semantics are not reinterpreted. Only owned-normal adoption's explicitly unreadable original execution is hash-checked without JSON decoding. Campaign preflight owns native semantic replay and exact historical revision policy. Known current reference-audit PDFs are parent-relative, never arbitrary JSON path strings.",
        "historical_semantics_replayed": False,
        "prior_conditional_rates_preserved": True, "prior_conditional_rates_recomputed": False,
        "corrected_rates_produced": False, "angular_partition_admitted": False,
        "rate_claim_allowed": False, "lifetime_claim_allowed": False,
        "physical_rate_accuracy_claim_allowed": False, "calibrated_uncertainty_claim_allowed": False,
        "network_reinjection_allowed": False, "production_defaults_changed": False,
        "native_jobs_dispatched": 0, "native_launch_allowed": False,
        "future_overlay_requirement": "A separately typed, reviewed physical overlay with source-bound semantic replay must be implemented; ordinary stationary-point validation remains unchanged.",
    }
    sources.finish()
    return receipt


def freeze_assessment(profile_result, destination, profile_contract=None):
    """Write a new refusal sidecar exclusively; never replace historical output."""
    receipt = assess(profile_result, profile_contract)
    _retained.write_new(Path(destination), receipt)
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Periodic-profile directory or result JSON; no native environment is needed")
    args = parser.parse_args(argv)
    print(json.dumps(assess(args.source), indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
