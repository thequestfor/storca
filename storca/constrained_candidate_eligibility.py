"""Native constrained-optimizer eligibility without a second displacement rule.

This source revision leaves historical acquisition policies untouched. Native
convergence supplies a candidate for independent final-frame gradient validation,
never evidence of the final saddle index, a verified reaction, or a rate.
"""
from __future__ import annotations

import json
import math
import re

from storca.ts_candidate_eligibility import _NAMES, _NUMBER

_LABEL = r"(?:Energy change|RMS gradient|MAX gradient|RMS step|MAX step)"
_ROW = re.compile(rf"\s*({_LABEL})\s+({_NUMBER})\s+({_NUMBER})\s+(YES|NO)\s*")
_GRADIENTS = ("RMS gradient", "MAX gradient")
_STEPS = ("RMS step", "MAX step")


def orca_constrained_candidate(text: str) -> dict:
    """Admit only native all-five or final-cycle alternate convergence.

    Both requested native gradient tests must pass in the alternate route.
    Printed step ratios are retained without an additional STORCA threshold.
    The caller must still bind ownership, state, mapping and final geometry,
    then independently evaluate its full constrained tangent gradient. Only
    separate full-Hessian evidence can establish a constrained saddle index.
    """
    if not isinstance(text, str):
        raise ValueError("Native optimizer text required")
    marker = "THE OPTIMIZATION HAS CONVERGED"
    normal = "ORCA TERMINATED NORMALLY"
    if (text.count(marker) != 1 or text.count(normal) != 1
            or text.index(normal) < text.index(marker)
            or "ORCA finished by error termination" in text):
        raise ValueError("Unique normal native optimization completion required")
    tables = list(re.finditer(r"\|Geometry convergence\|(.*?)(?:\.{6}|Max\(Bonds\))", text, re.S))
    if not tables or tables[-1].end() >= text.index(marker):
        raise ValueError("Final native convergence table must precede completion")
    table = tables[-1]
    rows = []
    for line in table.group(1).splitlines():
        if not re.match(rf"\s*{_LABEL}\b", line):
            continue
        match = _ROW.fullmatch(line)
        if match is None:
            raise ValueError("Malformed native geometry convergence row")
        rows.append(match.groups())
    if (len(rows) not in (4, 5) or len({row[0] for row in rows}) != len(rows)
            or not (_NAMES-{"Energy change"}).issubset(row[0] for row in rows)):
        raise ValueError("Malformed native geometry convergence inventory")
    criteria = {}
    for name, value, tolerance, passed in rows:
        actual = float(value.replace("D", "E").replace("d", "e"))
        tol = float(tolerance.replace("D", "E").replace("d", "e"))
        if (not math.isfinite(actual) or not math.isfinite(tol) or tol <= 0
                or (passed == "YES" and abs(actual) > tol)):
            raise ValueError("Nonfinite or contradictory native convergence criterion")
        criteria[name] = {"value": actual, "tolerance": tol, "native_passed": passed == "YES"}
    all_five = len(rows) == 5 and all(row["native_passed"] for row in criteria.values())
    decision = text[table.end():text.index(marker)]
    alternate = ("The gradient convergence is overachieved" in decision
                 and "Convergence will therefore be signaled now" in decision
                 and all(criteria[name]["native_passed"] for name in _GRADIENTS))
    if not all_five and not alternate:
        raise ValueError("Native optimization candidate eligibility unresolved")
    ratios = {name: abs(criteria[name]["value"])/criteria[name]["tolerance"] for name in _STEPS}
    if not all(math.isfinite(value) for value in ratios.values()):
        raise ValueError("Nonfinite reported displacement ratio")
    unmet = [name for name in _STEPS if not criteria[name]["native_passed"]]
    warnings = (["Native alternate convergence did not satisfy every printed displacement criterion; "
                 "independent final-frame tangent-gradient validation remains mandatory"] if unmet else [])
    result = {"candidate_eligible_for_derivative_validation": True,
              "scientifically_verified": False,
              "all_five_printed_criteria_passed": all_five,
              "native_convergence_route": "all_five" if all_five else "overachieved_gradient_candidate",
              "convergence_rows": rows, "criteria": criteria, "step_ratios": ratios,
              "native_displacement_criteria_unmet": unmet, "warnings": warnings,
              "extra_displacement_cutoff_applied": False,
              "displacement_bound_scope": "Reported diagnostics only; no additional STORCA displacement cutoff",
              "independent_final_geometry_validation_required": True,
              "independent_tangent_gradient_validation_required": True,
              "fixed_geometry_derivative_validation_required": True,
              "saddle_character_verified_by_optimization_alone": False,
              "independent_hessian_verified": False,
              "rate_claim_allowed": False, "lifetime_claim_allowed": False}
    json.dumps(result, allow_nan=False)
    return result
