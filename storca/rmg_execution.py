"""Fail-closed interpretation of retained RMG execution artifacts."""

from __future__ import annotations

import math
import re
from pathlib import Path


def assess_rmg_execution(log_file: Path, execution: dict | None = None) -> dict:
    """Classify whether RMG actually completed its requested model/simulation.

    A successful process return merely proves that the wrapper exited normally.
    RMG logs explicit resource/interruption messages which must not become
    persistence evidence.
    """
    execution = execution or {}
    result = {
        "wrapper_returned_success": (execution or {}).get("returncode", 0) == 0,
        "log_present": Path(log_file).is_file(),
        "status": "incomplete",
        "termination_reason": "missing_rmg_log",
        "evidence": [],
    }
    if not result["wrapper_returned_success"]:
        result.update(status="incomplete", termination_reason="rmg_process_failed")
        return result
    if not result["log_present"]:
        return result
    text = Path(log_file).read_text(errors="replace")
    lowered = text.lower()
    resurrections = lowered.count("resurrecting model")
    initial_step_failures = lowered.count("daspk returned with an idid")
    # RMG may recover once from a stiff start. Repeated recovery means the
    # reactor is not providing a numerically trustworthy bounded-network
    # result, even if the outer process later exits normally.
    if resurrections >= 2 or initial_step_failures >= 2:
        result.update(
            status="incomplete", termination_reason="repeated_solver_resurrection",
            evidence=[f"resurrections={resurrections}", f"daspk_initial_step_failures={initial_step_failures}"],
        )
        return result
    # ``Reached max number of objects`` and the following simulator interrupt
    # are normal RMG enlargement control flow when maxNumObjsPerIter is hit.
    # They become a failure only when RMG subsequently says that the model is
    # incomplete or that a real execution limit was reached.
    interrupted = [
        "output model may be incomplete",
        "there is not enough time to complete the next iteration",
        "maximum execution time reached",
        "time limit reached",
        "keyboardinterrupt",
    ]
    matches = [marker for marker in interrupted if marker in lowered]
    if matches:
        result.update(status="incomplete", termination_reason="resource_or_interrupt_termination", evidence=matches)
        return result
    # RMG's normal closing line has varied by version; accept either explicit
    # completion wording, but deliberately do not infer completion from a log.
    complete_markers = [
        "rmg execution terminated",
        "rmg execution completed",
        "finished rmg job",
        "execution complete",
    ]
    matches = [marker for marker in complete_markers if marker in lowered]
    if matches:
        result.update(status="completed", termination_reason="normal_completion", evidence=matches)
    else:
        result.update(status="incomplete", termination_reason="no_recognized_completion_marker")
    return result


def requested_time_coverage(profile: dict | None, requested_seconds: float) -> dict:
    """Return explicit simulated-time coverage; never extrapolate it."""
    end = (profile or {}).get("end_time_seconds")
    if end is None:
        return {"status": "not_available", "requested_seconds": requested_seconds, "simulated_seconds": None,
                "fraction": None, "complete": False}
    fraction = min(1.0, max(0.0, float(end) / requested_seconds)) if requested_seconds else 0.0
    return {"status": "completed" if end >= requested_seconds * (1 - 1e-9) else "incomplete",
            "requested_seconds": requested_seconds, "simulated_seconds": float(end), "fraction": fraction,
            "complete": end >= requested_seconds * (1 - 1e-9)}


def parse_collision_rate_violators(path: Path) -> dict:
    """Treat RMG collision-limit warnings as kinetic-validity failures."""
    path = Path(path)
    if not path.is_file():
        return {"status": "not_reported", "violation_count": 0, "path": str(path)}
    text = path.read_text(errors="replace")
    block_pattern = re.compile(
        r"^(?P<equation>[^!\n][^\n]*?(?:<=>|=>)[^\n]+?)\s+"
        r"[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s*$"
        r"(?P<trailer>.*?)(?=^[^!\n][^\n]*?(?:<=>|=>)[^\n]+?\s+"
        r"[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s*$|\Z)",
        re.MULTILINE | re.DOTALL,
    )
    violators = []
    for match in block_pattern.finditer(text):
        trailer = match.group("trailer")
        factor = re.search(r"Violation factor:\s*([0-9.Ee+-]+)", trailer)
        direction = re.search(r"Direction:\s*(forward|reverse)", trailer, re.IGNORECASE)
        condition = re.search(
            r"Violation condition:\s*([0-9.Ee+-]+)\s*K,\s*([0-9.Ee+-]+)\s*bar",
            trailer, re.IGNORECASE,
        )
        if factor:
            violators.append({
                "reaction_equation": match.group("equation").strip(),
                "direction": direction.group(1).lower() if direction else None,
                "violation_factor": float(factor.group(1)),
                "temperature_K": float(condition.group(1)) if condition else None,
                "pressure_bar": float(condition.group(2)) if condition else None,
            })
    factors = [item["violation_factor"] for item in violators]
    return {"status": "kinetics_unreliable" if factors else "passed", "violation_count": len(factors),
            "maximum_violation_factor": max(factors) if factors else None, "violators": violators, "path": str(path)}


def _finite_number(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def validate_inspected_collision_rates(
    mechanism_inspection: dict | None,
    *,
    temperature_K: float,
    pressure_bar: float,
    collision_tolerance: float = 1e-9,
) -> dict:
    """Validate every retained bimolecular direction at one condition.

    RMG does not always create ``collision_rate_violators.log``.  Absence of
    that optional report is not evidence that the mechanism is collision
    safe, while the retained bridge inspection already contains the evaluated
    forward/reverse rates and RMG collision ceilings.  This check is
    deliberately exhaustive and fail-closed: one missing rate, ceiling, or
    condition stamp keeps the result incomplete.
    """
    if collision_tolerance < 0.0:
        raise ValueError("Collision tolerance must be nonnegative")
    inspection = mechanism_inspection or {}
    if inspection.get("status") != "completed":
        return {
            "status": "incomplete",
            "applicable_direction_count": None,
            "validated_direction_count": 0,
            "violators": [],
            "missing_or_invalid_directions": [{"reason": "mechanism_inspection_missing_or_incomplete"}],
            "evaluated_at": {"temperature_K": temperature_K, "pressure_bar": pressure_bar},
        }

    checked: list[dict] = []
    missing: list[dict] = []
    violators: list[dict] = []
    applicable = 0
    for reaction_index, reaction in enumerate(inspection.get("reactions") or []):
        equation = str(reaction.get("equation") or reaction.get("reaction_equation") or "")
        evaluated_at = reaction.get("evaluated_at") or {}
        actual_temperature = _finite_number(evaluated_at.get("temperature_K"))
        actual_pressure = _finite_number(evaluated_at.get("pressure_bar"))
        condition_matches = (
            actual_temperature is not None
            and actual_pressure is not None
            and math.isclose(actual_temperature, float(temperature_K), rel_tol=1e-9, abs_tol=1e-9)
            and math.isclose(actual_pressure, float(pressure_bar), rel_tol=1e-9, abs_tol=1e-12)
        )
        directions = [(
            "forward", reaction.get("molecularity"),
            reaction.get("evaluated_forward_rate_coefficient_si"),
            reaction.get("forward_collision_limit_si"),
        )]
        # Generated photolysis and other irreversible source reactions have no
        # physical reverse kinetics to validate, even if a product count happens
        # to be two.  Older inspection artifacts lacked this flag and came only
        # from reversible Chemkin entries, so preserve that compatibility.
        if reaction.get("reversible", True):
            directions.append((
                "reverse", reaction.get("reverse_molecularity"),
                reaction.get("evaluated_reverse_rate_coefficient_si"),
                reaction.get("reverse_collision_limit_si"),
            ))
        for direction, raw_molecularity, raw_rate, raw_limit in directions:
            try:
                molecularity = int(raw_molecularity)
            except (TypeError, ValueError):
                missing.append({
                    "reaction_index": reaction.get("reaction_index", reaction_index),
                    "reaction_equation": equation,
                    "direction": direction,
                    "molecularity": None,
                    "rate_coefficient_si": _finite_number(raw_rate),
                    "collision_limit_si": _finite_number(raw_limit),
                    "temperature_K": actual_temperature,
                    "pressure_bar": actual_pressure,
                    "collision_limit_source": reaction.get("collision_limit_source"),
                    "reasons": ["reaction_molecularity_missing_or_invalid"],
                })
                continue
            if molecularity != 2:
                continue
            applicable += 1
            rate = _finite_number(raw_rate)
            limit = _finite_number(raw_limit)
            reasons = []
            if not condition_matches:
                reasons.append("evaluated_condition_missing_or_mismatched")
            if rate is None or rate < 0.0:
                reasons.append("condition_rate_missing_or_invalid")
            if limit is None or limit <= 0.0:
                reasons.append("collision_limit_missing_or_invalid")
            record = {
                "reaction_index": reaction.get("reaction_index", reaction_index),
                "reaction_equation": equation,
                "direction": direction,
                "molecularity": molecularity,
                "rate_coefficient_si": rate,
                "collision_limit_si": limit,
                "temperature_K": actual_temperature,
                "pressure_bar": actual_pressure,
                "collision_limit_source": reaction.get("collision_limit_source"),
            }
            if reasons:
                missing.append({**record, "reasons": reasons})
                continue
            ratio = rate / limit
            checked_record = {**record, "rate_to_collision_limit_ratio": ratio}
            checked.append(checked_record)
            if ratio > 1.0 + collision_tolerance:
                violators.append({
                    **checked_record,
                    "violation_factor": ratio,
                })

    if violators:
        status = "kinetics_unreliable"
    elif missing:
        status = "incomplete"
    else:
        status = "passed"
    return {
        "status": status,
        "applicable_direction_count": applicable,
        "validated_direction_count": len(checked),
        "maximum_rate_to_collision_limit_ratio": (
            max((item["rate_to_collision_limit_ratio"] for item in checked), default=None)
        ),
        "violators": violators,
        "missing_or_invalid_directions": missing,
        "evaluated_at": {"temperature_K": temperature_K, "pressure_bar": pressure_bar},
        "rule": (
            "Pass only when every retained forward/reverse direction with bimolecular reactants has a "
            "finite condition rate and finite positive RMG collision ceiling, and none exceeds that ceiling."
        ),
    }


def merge_collision_rate_validation(
    reported_validation: dict | None,
    mechanism_inspection: dict | None,
    *,
    temperature_K: float,
    pressure_bar: float,
) -> dict:
    """Merge RMG's optional warning log with exhaustive condition checks.

    Explicit RMG warnings always survive the merge.  A clean or absent log can
    become ``passed`` only when the retained mechanism inspection independently
    validates every applicable bimolecular direction.
    """
    reported = dict(reported_validation or {})
    inspected = validate_inspected_collision_rates(
        mechanism_inspection,
        temperature_K=temperature_K,
        pressure_bar=pressure_bar,
    )
    reported_violators = list(reported.get("violators") or [])
    inspected_violators = list(inspected.get("violators") or [])
    explicit_violation = (
        reported.get("status") == "kinetics_unreliable"
        or int(reported.get("violation_count") or 0) > 0
        or bool(reported_violators)
    )
    if explicit_violation or inspected.get("status") == "kinetics_unreliable":
        status = "kinetics_unreliable"
    elif inspected.get("status") == "passed":
        status = "passed"
    else:
        status = "incomplete"
    violators = [*reported_violators, *inspected_violators]
    factors = [
        value for value in (_finite_number(item.get("violation_factor")) for item in violators)
        if value is not None
    ]
    return {
        **reported,
        "status": status,
        "violation_count": len(violators) if violators else int(reported.get("violation_count") or 0),
        "maximum_violation_factor": max(factors) if factors else reported.get("maximum_violation_factor"),
        "violators": violators,
        "reported_log_validation": reported,
        "condition_mechanism_validation": inspected,
        "validation_basis": "reported_collision_warnings_plus_exhaustive_condition_mechanism_inspection",
    }
