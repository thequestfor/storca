"""Ordered, retained condition-ladder orchestration for stability screens."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

from .conditions import ConditionSpec, build_condition_spec
from .photolysis import evaluate_sunlight_photolysis, sunlight_contract
from .runs import write_metadata


SECONDS_PER_DAY = 86400.0


_MODEL_ONLY_STABILITY_STATUSES = {
    "stable",
    "no_loss_detected_in_rmg_model",
    "model_supported_persistence",
    "verified_route_below_loss_threshold",
}
_NO_LOSS_MODEL_OUTCOMES = {
    "no_target_loss_above_threshold_in_retained_rmg_model",
    "verified_below_loss_threshold_in_condition_model",
}
_VERIFIED_T95_COMPOSITE_STATUSES = {
    "completed_verified_converged",
    "orca_verified_condition_specific_t95",
    "verified_t95_converged",
    "verified_condition_specific_t95",
}
_VERIFIED_ROUTE_STATUSES = {
    "barrierless_capture_verified",
    "irc_endpoint_verified",
    "orca_verified",
    "verified",
}
_VERIFIED_RATE_STATUSES = {
    "arkane_verified",
    "collision_bounded_verified",
    "rate_verified",
    "validated",
    "verified",
}
_CONVERGED_PROPAGATION_STATUSES = {
    "completed_converged",
    "converged",
    "verified_converged",
}
_VERIFICATION_BLOCK_KEYS = (
    "verification_summary",
    "generic_route_verification",
    "generic_verification",
    "route_verification_summary",
    "verified_kinetics",
    "kinetic_verification",
)


def _positive_finite(value: Any) -> float | None:
    """Return a physically usable lifetime, never a bool/NaN/infinity."""
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and number > 0.0 else None


def _verification_blocks(stage_result: dict) -> list[dict]:
    blocks = [stage_result[key] for key in _VERIFICATION_BLOCK_KEYS
              if isinstance(stage_result.get(key), dict)]
    lifetime = stage_result.get("kinetic_lifetime")
    if isinstance(lifetime, dict):
        blocks.append(lifetime)
    return blocks


def _nested_value(block: dict, *paths: tuple[str, ...]) -> Any:
    for path in paths:
        value: Any = block
        for key in path:
            if not isinstance(value, dict) or key not in value:
                break
            value = value[key]
        else:
            if value is not None:
                return value
    return None


def _condition_contract_matches(stage_result: dict, block: dict) -> bool:
    explicit = _nested_value(
        block,
        ("condition_contract_match",),
        ("condition_specific",),
        ("condition_validation", "matches_stage_contract"),
    )
    if explicit is False:
        return False
    expected = stage_result.get("condition_contract")
    retained = _nested_value(block, ("condition_contract",), ("conditions",))
    if expected is not None and retained is not None:
        return expected == retained
    return explicit is True


def verified_t95_evidence(stage_result: dict) -> dict | None:
    """Return a terminal lifetime only with explicit verification and convergence.

    An assessment label or an RMG screening estimate is deliberately
    insufficient.  A generic verifier may supply one composite status, or the
    equivalent route/rate/propagation fields.  In both forms it must also bind
    the result to this stage's immutable condition contract.
    """
    stage_status = str(stage_result.get("status", ""))
    if stage_status.startswith("not_run") or stage_status in {
        "failed", "incomplete", "intrinsic_kinetics_incomplete",
        "computational_light_incomplete", "photochemical_evidence_incomplete",
    }:
        return None
    for evidence_name in ("orca_evidence", "rmg_evidence"):
        evidence = stage_result.get(evidence_name)
        if isinstance(evidence, dict) and evidence.get("status") != "completed":
            return None
    for block in _verification_blocks(stage_result):
        lifetime = _positive_finite(_nested_value(
            block,
            ("estimated_time_to_retention_seconds",),
            ("t95_seconds",),
            ("kinetic_lifetime", "estimated_time_to_retention_seconds"),
        ))
        if lifetime is None:
            continue
        status = str(block.get("status") or block.get("verification_status") or "")
        composite = status in _VERIFIED_T95_COMPOSITE_STATUSES
        route_status = str(_nested_value(
            block,
            ("route_verification_status",),
            ("route_status",),
            ("route_verification", "status"),
        ) or "")
        rate_status = str(_nested_value(
            block,
            ("rate_verification_status",),
            ("rate_status",),
            ("kinetics", "status"),
        ) or "")
        propagation_status = str(_nested_value(
            block,
            ("propagation_status",),
            ("kinetic_model_status",),
            ("propagation", "status"),
        ) or "")
        split_verified = (
            route_status in _VERIFIED_ROUTE_STATUSES
            and rate_status in _VERIFIED_RATE_STATUSES
            and propagation_status in _CONVERGED_PROPAGATION_STATUSES
        )
        # The canonical status names convergence itself.  Older composite
        # records still need their separate propagation field; this prevents a
        # route-only ORCA result from being mistaken for a propagated t95.
        converged = (
            status == "verified_t95_converged"
            or (composite and propagation_status in _CONVERGED_PROPAGATION_STATUSES)
        )
        if (split_verified or converged) and _condition_contract_matches(stage_result, block):
            return {
                "estimated_time_to_retention_seconds": lifetime,
                "verification_status": status or "split_verification_contract",
                "propagation_status": propagation_status,
                "source": next((key for key in _VERIFICATION_BLOCK_KEYS
                                if stage_result.get(key) is block), "kinetic_lifetime"),
            }
    # Compatibility for retained v1 reports: their final lifetime field was
    # contractually null until ORCA route verification and repaired propagation
    # completed.  Require the explicit ORCA-verified alias, both completed
    # engines, and a retained condition contract; an assessment string alone is
    # never sufficient.
    assessment = str((stage_result.get("assessment") or {}).get("status") or "")
    lifetime = _positive_finite((stage_result.get("kinetic_lifetime") or {}).get(
        "estimated_time_to_retention_seconds"
    ))
    if (
        assessment in {"orca_verified_intrinsic_instability", "orca_verified_condition_dependent_instability"}
        and lifetime is not None
        and isinstance(stage_result.get("condition_contract"), dict)
        and all(
            not isinstance(stage_result.get(name), dict)
            or stage_result[name].get("status") == "completed"
            for name in ("orca_evidence", "rmg_evidence")
        )
    ):
        return {
            "estimated_time_to_retention_seconds": lifetime,
            "verification_status": assessment,
            "propagation_status": "legacy_contract_implied_converged",
            "source": "legacy_v1_assessment_and_lifetime_contract",
        }
    return None


def _verified_below_threshold_evidence(stage_result: dict) -> dict | None:
    """Validate the canonical completed/below-threshold verifier summary."""
    summary = stage_result.get("verification_summary")
    if not isinstance(summary, dict) or summary.get("status") != "verified_below_loss_threshold":
        return None
    if not _condition_contract_matches(stage_result, summary):
        return None
    target_loss = _positive_finite(summary.get("target_loss_fraction"))
    if summary.get("target_loss_fraction") == 0:
        target_loss = 0.0
    unresolved = _positive_finite(summary.get("unresolved_loss_upper_bound"))
    if summary.get("unresolved_loss_upper_bound") == 0:
        unresolved = 0.0
    coverage = summary.get("verified_flux_coverage")
    try:
        coverage_value = float(coverage)
    except (TypeError, ValueError):
        return None
    if (target_loss is None or unresolved is None or not math.isfinite(coverage_value)
            or not 0.0 <= coverage_value <= 1.0):
        return None
    retention = (stage_result.get("condition_contract") or {}).get("retention_fraction")
    try:
        loss_threshold = 1.0 - float(retention)
    except (TypeError, ValueError):
        return None
    if not 0.0 <= loss_threshold <= 1.0 or target_loss + unresolved >= loss_threshold:
        return None
    return {
        "status": "verified_below_loss_threshold",
        "target_loss_fraction": target_loss,
        "unresolved_loss_upper_bound": unresolved,
        "verified_flux_coverage": coverage_value,
    }


def _gap(code: str, detail: str) -> dict:
    return {"code": code, "detail": detail}


def _stage_evidence_gaps(result: dict) -> list[dict]:
    """Collect actionable omissions without treating a model no-loss as one."""
    gaps: list[dict] = []
    status = str(result.get("status") or "")
    assessment = str((result.get("assessment") or {}).get("status") or "")
    verification = result.get("verification_summary")
    if isinstance(verification, dict):
        verification_status = verification.get("status")
        if verification_status == "verification_incomplete":
            gaps.append(_gap(
                "generic_route_verification_incomplete",
                verification.get("reason") or "The generic route verifier did not complete.",
            ))
        elif verification_status == "verification_disabled":
            gaps.append(_gap(
                "generic_route_verification_disabled",
                verification.get("reason") or "Automatic generic route verification was disabled.",
            ))
        elif verification_status == "verified_t95_converged" and verified_t95_evidence(result) is None:
            gaps.append(_gap(
                "verified_t95_contract_invalid",
                "The verifier claims a converged t95, but the lifetime or condition-contract binding is invalid.",
            ))
        elif (verification_status == "verified_below_loss_threshold"
              and _verified_below_threshold_evidence(result) is None):
            gaps.append(_gap(
                "below_threshold_coverage_invalid",
                "The below-threshold result lacks valid flux coverage, an unresolved-loss bound, or a matching condition contract.",
            ))
    if status == "intrinsic_kinetics_incomplete":
        gaps.append(_gap(
            "intrinsic_kinetics_incomplete",
            result.get("reason") or "Intrinsic routes lack verified pressure-dependent kinetics.",
        ))
    elif status in {"computational_light_incomplete", "photochemical_evidence_incomplete"}:
        gaps.append(_gap("photochemical_evidence_incomplete", "The photochemical model did not complete."))
    elif (result.get("model_outcome") == "no_reactive_photon_channel_in_computational_light_model"):
        gaps.append(_gap(
            "photochemical_channel_not_validated",
            "The computational screen admitted no reactive photon channel; this is not evidence of photostability.",
        ))
    elif status in {"failed", "incomplete", "completed_with_incomplete_evidence"}:
        gaps.append(_gap("stage_evidence_incomplete", "One or more required stage calculations did not complete."))

    for name, label in (("orca_evidence", "ORCA"), ("rmg_evidence", "RMG")):
        evidence = result.get(name)
        if isinstance(evidence, dict) and evidence.get("status") != "completed":
            gaps.append(_gap(f"{name}_incomplete", f"{label} evidence did not complete."))

    assessment_gaps = {
        "incomplete_evidence": ("required_evidence_incomplete", "The combined stage evidence is incomplete."),
        "kinetics_unreliable": ("rmg_kinetics_unreliable", "Retained RMG kinetics failed physical validation."),
        "likely_air_reactive_lifetime_unverified": ("route_rate_unverified", "The controlling route rate has not been physically verified."),
        "orca_verification_required": ("route_verification_required", "The controlling route still requires ORCA verification."),
        "orca_verification_incomplete": ("route_verification_incomplete", "ORCA route verification did not complete."),
        "orca_initiation_ts_verification_required": ("initiation_ts_kinetics_required", "The initiating route still needs TS/IRC and rate verification."),
        "orca_supported_barrierless_air_reactivity_lifetime_unverified": ("barrierless_rate_unverified", "Barrierless route physics is available, but its condition-specific rate is not."),
    }
    if assessment in assessment_gaps:
        code, detail = assessment_gaps[assessment]
        reason = (result.get("assessment") or {}).get("reason")
        gaps.append(_gap(code, reason or detail))

    lifetime = result.get("kinetic_lifetime") or {}
    screening_t95 = _positive_finite(lifetime.get("screening_estimated_time_to_retention_seconds"))
    reported_t95 = _positive_finite(lifetime.get("estimated_time_to_retention_seconds"))
    if screening_t95 is not None:
        gaps.append(_gap(
            "screening_t95_unverified",
            "RMG produced a screening t95, but verified route kinetics have not been propagated.",
        ))
    if reported_t95 is not None and verified_t95_evidence(result) is None:
        gaps.append(_gap(
            "reported_t95_not_verified_converged",
            "A numerical t95 is present without an explicit condition-matched, verified, converged propagation record.",
        ))

    unique: dict[str, dict] = {}
    for item in gaps:
        unique.setdefault(item["code"], item)
    return list(unique.values())


def _normalize_stage_result(result: dict) -> dict:
    """Make model scope and unresolved evidence visible at the stage boundary."""
    normalized = dict(result)
    assessment = normalized.get("assessment")
    if isinstance(assessment, dict) and assessment.get("status") in _MODEL_ONLY_STABILITY_STATUSES:
        normalized["source_assessment"] = dict(assessment)
        normalized["assessment"] = {
            "status": "no_loss_detected_in_rmg_model",
            "reason": (
                f"{assessment.get('reason', 'The retained model did not cross the loss threshold.')} "
                "This is a model-scoped no-loss result, not an overall stability conclusion."
            ),
        }
        normalized["model_outcome"] = "no_target_loss_above_threshold_in_retained_rmg_model"

    verified = verified_t95_evidence(normalized)
    verified_below = _verified_below_threshold_evidence(normalized)
    gaps = _stage_evidence_gaps(normalized)
    if verified is not None:
        superseded = {
            "barrierless_rate_unverified", "initiation_ts_kinetics_required",
            "reported_t95_not_verified_converged", "route_rate_unverified",
            "route_verification_incomplete", "route_verification_required",
            "screening_t95_unverified",
        }
        gaps = [item for item in gaps if item["code"] not in superseded]
        prior = normalized.get("assessment")
        if isinstance(prior, dict):
            normalized["source_assessment"] = dict(prior)
        normalized["assessment"] = {
            "status": "orca_verified_condition_dependent_instability",
            "reason": "Verified route kinetics were propagated to a condition-matched, converged t95.",
        }
        normalized.pop("model_outcome", None)
    elif verified_below is not None:
        prior = normalized.get("assessment")
        if isinstance(prior, dict):
            normalized["source_assessment"] = dict(prior)
        normalized["assessment"] = {
            "status": "verified_route_below_loss_threshold",
            "reason": (
                "Verified controlling-route kinetics and the unresolved-loss bound remain below "
                "the declared retention threshold in this condition model."
            ),
        }
        normalized["verified_below_threshold_evidence"] = verified_below
        normalized["model_outcome"] = "verified_below_loss_threshold_in_condition_model"
    normalized["missing_evidence"] = gaps
    if verified is not None:
        normalized["verified_t95_evidence"] = verified
        normalized["status"] = "completed_with_verified_t95"
    elif verified_below is not None and not gaps:
        normalized["status"] = "completed_verified_below_loss_threshold"
    elif normalized.get("model_outcome") and not gaps:
        normalized["status"] = "completed_model_screen_no_loss"
    elif gaps and normalized.get("status") not in {
        "intrinsic_kinetics_incomplete", "computational_light_incomplete", "photochemical_evidence_incomplete",
    }:
        normalized["status"] = "completed_with_incomplete_evidence"
    return normalized


def water_vapor_pressure_bar(temperature_K: float) -> float:
    """Buck saturation-vapour-pressure correlation for liquid water (0–50 C)."""
    celsius = temperature_K - 273.15
    if not 0.0 <= celsius <= 50.0:
        raise ValueError("Humid-air stage currently supports temperatures from 273.15 through 323.15 K")
    import math
    return 0.0061121 * math.exp((18.678 - celsius / 234.5) * (celsius / (257.14 + celsius)))


def humid_air_scenario(*, temperature_K: float, pressure_bar: float, relative_humidity: float) -> dict:
    """Return an RMG homogeneous-gas composition for humid air."""
    if not 0.0 <= relative_humidity <= 1.0:
        raise ValueError("relative_humidity must be a fraction from zero through one")
    water_fraction = relative_humidity * water_vapor_pressure_bar(temperature_K) / pressure_bar
    if water_fraction >= 0.99:
        raise ValueError("Humidity/pressure combination leaves no meaningful dry-air fraction")
    target_fraction = 0.01
    remaining = 1.0 - target_fraction - water_fraction
    if remaining <= 0:
        raise ValueError("Humidity leaves no room for the target composition")
    return {
        "name": "humid-air-gas-screen",
        "phase": "homogeneous gas-phase surrogate",
        "atmosphere": "humid air (water vapour/oxygen/nitrogen)",
        "model_applicability": "A dilute target in homogeneous humid gas. It does not model droplets, condensed water, dissolved ions, surfaces, or container compatibility.",
        "additional_species": [
            {"label": "water", "smiles": "O", "reactive": True},
            {"label": "oxygen", "smiles": "[O][O]", "reactive": True},
            {"label": "nitrogen", "smiles": "N#N", "reactive": False},
        ],
        "initial_mole_fractions": {"stability": target_fraction, "water": water_fraction,
                                   "oxygen": remaining * 0.2095, "nitrogen": remaining * 0.7905},
        "relative_humidity": relative_humidity,
        "water_vapor_partial_pressure_bar": water_fraction * pressure_bar,
    }


def build_default_ladder(*, temperature_K: float = 298.15, pressure_bar: float = 1.0,
                         relative_humidity: float = 0.5, retention_fraction: float = 0.95,
                         maximum_duration_days: float = 365.0) -> list[dict]:
    """Build fixed-order stages; later stages are retained as untested after t95."""
    humid = humid_air_scenario(temperature_K=temperature_K, pressure_bar=pressure_bar,
                               relative_humidity=relative_humidity)
    common = {"temperature_K": temperature_K, "pressure_bar": pressure_bar,
              "retention_fraction": retention_fraction,
              "maximum_duration_seconds": maximum_duration_days * SECONDS_PER_DAY}
    return [
        {"id": "dark-low-pressure", "kind": "rmg", "scenario": "low-pressure-intrinsic-gas-screen",
         "pressure_bar": 1e-6, **{k: v for k, v in common.items() if k != "pressure_bar"}},
        {"id": "dark-dry-inert", "kind": "rmg", "scenario": "dry-inert-gas-screen", **common},
        {"id": "dark-dry-air", "kind": "rmg", "scenario": "ambient-air-gas-screen", **common},
        {"id": "dark-humid-air", "kind": "rmg", "scenario_config": humid, "relative_humidity": relative_humidity, **common},
        {"id": "sunlight-dry-air", "kind": "photolysis", "scenario": "ambient-air-gas-screen", "light_model": sunlight_contract(), **common},
        {"id": "sunlight-humid-air", "kind": "photolysis", "scenario_config": humid, "relative_humidity": relative_humidity, "light_model": sunlight_contract(), **common},
    ]


def _t95(stage_result: dict) -> float | None:
    # A standalone RMG t95 or photolysis J is screening evidence only.  The
    # terminal gate lives in one place so every stage must meet the identical
    # ORCA/rate/condition/converged-propagation contract.
    evidence = verified_t95_evidence(stage_result)
    return evidence["estimated_time_to_retention_seconds"] if evidence else None


def run_condition_ladder(smiles: str, run_dir: Path, *, stage_runner: Callable[..., dict],
                         temperature_K: float = 298.15, pressure_bar: float = 1.0,
                         relative_humidity: float = 0.5, retention_fraction: float = 0.95,
                         maximum_duration_days: float = 365.0,
                         photolysis_evidence: Path | None = None,
                         computational_light_spectrum: Path | None = None, **runner_kwargs) -> dict:
    """Run ordered stages and stop only after a verified, converged t95."""
    run_dir = Path(run_dir)
    stages = build_default_ladder(temperature_K=temperature_K, pressure_bar=pressure_bar,
                                  relative_humidity=relative_humidity, retention_fraction=retention_fraction,
                                  maximum_duration_days=maximum_duration_days)
    results, terminal_index = [], None
    for index, stage in enumerate(stages):
        if terminal_index is not None:
            results.append({"id": stage["id"], "status": "not_run_after_verified_t95", "condition": stage,
                            "missing_evidence": []})
            continue
        stage_dir = run_dir / "ladder" / f"{index + 1:02d}-{stage['id']}"
        stage_dir.mkdir(parents=True, exist_ok=True)
        print(f"[STORCA] Ladder stage {index + 1}/{len(stages)} started: {stage['id']}", flush=True)
        if stage["kind"] == "photolysis":
            if photolysis_evidence is None and computational_light_spectrum is not None:
                from .computational_light import run_computational_light_model, simulate_computational_light_profiles
                light = run_computational_light_model(
                    smiles, stage_dir / "computational-light", sunlight_spectrum=computational_light_spectrum,
                    rmg_env=runner_kwargs.get("rmg_env") or "rmg_env",
                    charge=runner_kwargs.get("charge", 0),
                    multiplicity=runner_kwargs.get("multiplicity", 1),
                    ncores=runner_kwargs.get("ncores", 1),
                    nroots=runner_kwargs.get("light_nroots", 20),
                    method_keywords=runner_kwargs.get("method_keywords"),
                )
                if light["status"] != "completed":
                    result = {"id": stage["id"], "status": "computational_light_incomplete", "condition": stage,
                              "computational_light": light}
                else:
                    kinetics = simulate_computational_light_profiles(
                        light, stage_dir / "computational-light-rmg", rmg_env=runner_kwargs.get("rmg_env") or "rmg_env",
                        scenario=stage.get("scenario") or "ambient-air-gas-screen",
                        scenario_config=stage.get("scenario_config"),
                        temperature=temperature_K, pressure=pressure_bar,
                        target_duration_hours=maximum_duration_days * 24, retention_fraction=retention_fraction,
                    )
                    central = kinetics["profiles"].get("central", {})
                    if kinetics.get("status") == "completed_no_reactive_photon_channel":
                        result = {
                            "id": stage["id"],
                            "status": "completed_with_incomplete_evidence",
                            "condition": stage,
                            "computational_light": light,
                            "computational_light_kinetics": kinetics,
                            "assessment": {
                                "status": "photochemical_channel_not_validated",
                                "reason": (
                                    "TD-DFT and all configured photon profiles completed, but no route met "
                                    "the declared source-window and accessibility admission criteria. This does not prove photostability."
                                ),
                            },
                            "model_outcome": "no_reactive_photon_channel_in_computational_light_model",
                            "kinetic_lifetime": {
                                "estimated_time_to_retention_seconds": None,
                                "interpretation": (
                                    "No admitted photon-driven channel was available for propagation under "
                                    "the immutable computational-light model."
                                ),
                            },
                        }
                    else:
                        result = {
                            "id": stage["id"],
                            "status": (
                                "completed"
                                if kinetics.get("status") == "completed" and central.get("status") == "completed"
                                else "completed_with_incomplete_evidence"
                            ),
                            "condition": stage,
                            "computational_light": light,
                            "computational_light_kinetics": kinetics,
                            "kinetic_lifetime": {
                                "estimated_time_to_retention_seconds": central.get("estimated_time_to_retention_seconds"),
                                "interpretation": "Central branch-profile projection; consult low/high profiles for sensitivity.",
                            },
                        }
                result = _normalize_stage_result(result)
                results.append(result)
                if _t95(result) is not None:
                    terminal_index = index
                continue
            scenario = stage.get("scenario_config")
            if scenario is None:
                from .stability import resolve_stability_configuration
                scenario, _ = resolve_stability_configuration(stage["scenario"], "quick-screen")
            condition = build_condition_spec(scenario, temperature_K=temperature_K, pressure_bar=pressure_bar,
                                             target_duration_hours=maximum_duration_days * 24,
                                             retention_fraction=retention_fraction, light_condition="sunlight",
                                             relative_humidity=stage.get("relative_humidity"), light_model=stage["light_model"])
            photolysis = evaluate_sunlight_photolysis(condition=condition.as_dict(), photolysis_evidence=photolysis_evidence)
            if photolysis["status"] != "completed":
                result = {"id": stage["id"], **photolysis}
            else:
                from .generated_kinetics import write_photolysis_library
                library = write_photolysis_library(
                    stage_dir / "generated-photolysis", route_id=stage["id"], reactant_label="stability", reactant_smiles=smiles,
                    products=[(item["label"], item["smiles"]) for item in photolysis["photoproducts"]],
                    photolysis_rate_constant_s_1=photolysis["photolysis_rate_constant_s_1"],
                    photolysis_evidence=photolysis_evidence,
                    rmg_env=runner_kwargs.get("rmg_env") or "rmg_env",
                )
                stage_result = stage_runner(
                    smiles, stage_dir, scenario=stage.get("scenario"), scenario_config=stage.get("scenario_config"),
                    target_duration_hours=maximum_duration_days * 24, retention_fraction=retention_fraction,
                    light_condition="sunlight", light_model=stage["light_model"], reaction_libraries=[Path(library["library"])],
                    temperature=stage["temperature_K"], pressure=stage["pressure_bar"], **runner_kwargs,
                )
                complete = (stage_result.get("orca_evidence", {}).get("status") == "completed"
                            and stage_result.get("rmg_evidence", {}).get("status") == "completed")
                result = {"id": stage["id"], "status": "completed" if complete else "completed_with_incomplete_evidence",
                          "condition": stage, "photolysis": photolysis, "generated_kinetics_library": library, **stage_result}
        else:
            stage_result = stage_runner(smiles, stage_dir, scenario=stage.get("scenario"),
                                        scenario_config=stage.get("scenario_config"),
                                        target_duration_hours=maximum_duration_days * 24,
                                        retention_fraction=retention_fraction,
                                        temperature=stage["temperature_K"], pressure=stage["pressure_bar"],
                                        **runner_kwargs)
            complete = (stage_result.get("orca_evidence", {}).get("status") == "completed"
                        and stage_result.get("rmg_evidence", {}).get("status") == "completed")
            result = {"id": stage["id"], "status": "completed" if complete else "completed_with_incomplete_evidence",
                      "condition": stage, **stage_result}
        result = _normalize_stage_result(result)
        results.append(result)
        print(f"[STORCA] Ladder stage {index + 1}/{len(stages)} finished: {stage['id']} ({result.get('status')})", flush=True)
        if _t95(result) is not None:
            terminal_index = index
    missing_evidence = [
        {"stage_id": stage["id"], **gap}
        for stage in results for gap in stage.get("missing_evidence", [])
    ]
    model_no_loss_stages = [
        stage["id"] for stage in results
        if stage.get("model_outcome") in _NO_LOSS_MODEL_OUTCOMES
    ]
    verified_below_threshold_stages = [
        stage["id"] for stage in results if stage.get("verified_below_threshold_evidence")
    ]
    if terminal_index is not None:
        verdict = "condition_scoped_t95_identified"
        evidence_status = ("condition_specific_t95_verified_with_other_evidence_gaps"
                           if missing_evidence else "condition_specific_t95_verified")
        overall_reason = (
            "A condition-matched ORCA-verified rate was propagated with a converged kinetic model to the retention threshold."
        )
    elif missing_evidence:
        verdict = "no_verified_t95_with_incomplete_evidence"
        evidence_status = "incomplete"
        overall_reason = (
            "No verified t95 is available, and at least one required intrinsic, route, kinetic, or model-coverage step is incomplete."
        )
    elif len(model_no_loss_stages) == len(stages) and len(results) == len(stages):
        verdict = "stable_under_tested_conditions"
        evidence_status = "completed_model_no_loss"
        overall_reason = (
            "Every applicable ladder condition completed without a modeled loss above the declared "
            "retention threshold. This is stability under the retained condition contracts, not a universal claim."
        )
    elif verified_below_threshold_stages:
        verdict = "verified_below_loss_threshold_within_completed_models"
        evidence_status = "completed_verified_below_loss_threshold"
        overall_reason = (
            "Verified route coverage and conservative unresolved-loss bounds remained below the retention threshold in the completed condition models."
        )
    else:
        verdict = "no_verified_t95_within_completed_models"
        evidence_status = "completed_without_verified_t95"
        overall_reason = (
            "The completed supported models did not produce a verified condition-specific t95; this is not a universal stability claim."
        )
    report = {"schema_version": 2, "kind": "condition_ladder_stability_screen", "smiles": smiles,
              "verdict": verdict, "first_t95_stage": results[terminal_index]["id"] if terminal_index is not None else None,
              "overall_assessment": {"status": evidence_status, "reason": overall_reason},
              "evidence_summary": {"status": evidence_status,
                                   "missing_evidence": missing_evidence,
                                   "model_no_loss_stages": model_no_loss_stages,
                                   "verified_below_threshold_stages": verified_below_threshold_stages,
                                   "verified_t95": (results[terminal_index].get("verified_t95_evidence")
                                                    if terminal_index is not None else None)},
              "stages": results,
              "limitations": ["A stage with incomplete intrinsic or photochemical kinetics is not a stability finding.",
                              "A stable-under-tested-conditions verdict is scoped to the retained condition contracts and is not a universal stability finding.",
                              "Humid-air stages model water vapour only, not liquid or surface chemistry."],}
    path = run_dir / "stability-ladder.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    write_metadata(run_dir, workflow="condition_ladder_stability_screen", result_json=str(path),
                   ladder_verdict=verdict)
    return {**report, "result_json": path}
