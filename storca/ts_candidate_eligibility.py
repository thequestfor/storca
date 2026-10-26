"""Native optimizer receipts admit candidates, never verify physical saddles.

This parser is deliberately separate from fixed-geometry derivative validation.
Its current adapter recognizes ORCA's all-five and overachieved-gradient routes.
It does not certify a rate, a reaction path, or the quality of an approximate
optimizer Hessian. Electronic state, geometry and process ownership belong to
the caller's bound acquisition contract.
"""
from __future__ import annotations

import math
import re

_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
_ROW = re.compile(r"(Energy change|RMS gradient|MAX gradient|RMS step|MAX step)\s+"
                  rf"({_NUMBER})\s+({_NUMBER})\s+(YES|NO)\b")
_NAMES = {"Energy change", "RMS gradient", "MAX gradient", "RMS step", "MAX step"}


def orca_optimization_candidate(text: str, *, max_step_multiple: float = 2.0) -> dict:
    """Return eligibility for independent derivatives, with truthful criteria.

    The optional displacement bound is an acquisition guard, not stationarity
    evidence. Both final-cycle native gradient tests must pass for the alternate
    route. No extra percentage of an already requested gradient tolerance is
    imposed. Alternate markers from earlier cycles cannot authorize this one.
    """
    if not math.isfinite(max_step_multiple) or max_step_multiple < 1:
        raise ValueError("Finite displacement acquisition bound >=1 required")
    marker = "THE OPTIMIZATION HAS CONVERGED"
    if (text.count(marker) != 1 or text.count("ORCA TERMINATED NORMALLY") != 1
            or text.index("ORCA TERMINATED NORMALLY") < text.index(marker)):
        raise ValueError("Unique normal native optimization completion required")
    tables = list(re.finditer(r"\|Geometry convergence\|(.*?)(?:\.{6}|Max\(Bonds\))", text, re.S))
    if not tables or tables[-1].end() >= text.index(marker):
        raise ValueError("Final native convergence table must precede completion")
    table = tables[-1]
    rows = _ROW.findall(table.group(1))
    if (len(rows) not in (4, 5) or len({row[0] for row in rows}) != len(rows)
            or len(re.findall(r"(?m)^\s*(?:Energy change|RMS gradient|MAX gradient|RMS step|MAX step)\b", table.group(1))) != len(rows)
            or not _NAMES.difference({"Energy change"}).issubset(row[0] for row in rows)):
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
    final_decision = text[table.end():text.index(marker)]
    alternate = ("The gradient convergence is overachieved" in final_decision
                 and "Convergence will therefore be signaled now" in final_decision
                 and all(criteria[name]["native_passed"] for name in ("RMS gradient", "MAX gradient"))
                 and all(abs(criteria[name]["value"]) <= max_step_multiple*criteria[name]["tolerance"]
                         for name in ("RMS step", "MAX step")))
    if not all_five and not alternate:
        raise ValueError("Native optimization candidate eligibility unresolved")
    return {"candidate_eligible_for_derivative_validation": True,
            "scientifically_verified": False,
            "all_five_printed_criteria_passed": all_five,
            "native_convergence_route": "all_five" if all_five else "overachieved_gradient_candidate",
            "convergence_rows": rows, "criteria": criteria,
            "displacement_bound_scope": "candidate acquisition guard, not physical stationarity",
            "saddle_character_verified_by_optimization_alone": False,
            "fixed_geometry_derivative_validation_required": True,
            "rate_claim_allowed": False, "lifetime_claim_allowed": False}
