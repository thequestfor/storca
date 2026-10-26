"""Explicit evidence contracts for rate uncertainty and verification cost.

These records are allocation evidence, not rate coefficients and not a claim
that the resulting lifetime uncertainty has been propagated.  STORCA never
invents a default uncertainty or cost when a record is absent.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping


_SOURCE_KINDS = {
    "database",
    "method_calibration",
    "model_ensemble",
    "physical_bound",
    "user_supplied",
}


def _finite(value: Any, *, field: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{field} must be a finite number") from error
    if not math.isfinite(number):
        raise ValueError(f"{field} must be a finite number")
    return number


def _optional_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


@dataclass(frozen=True)
class RateUncertaintyRecord:
    """A provenance-bearing multiplicative uncertainty bound for one route.

    Method calibration is admitted only with an explicit exact-condition
    binding; descriptive applicability prose is never an executable domain.
    The binding is not authentication of the source or method recipe. Explicit
    user/database allocation records retain their optional condition scope.
    """

    lower_rate_multiplier: float
    upper_rate_multiplier: float
    source: str
    source_kind: str
    applicability: str
    route_id: str | None = None
    reaction_equation: str | None = None
    condition_contract_id: str | None = None
    confidence_level: float | None = None
    uncertainty_components: tuple[str, ...] = ()
    estimated_cost: float | None = None
    cost_units: str | None = None
    cost_model: str | None = None

    def __post_init__(self) -> None:
        lower = _finite(self.lower_rate_multiplier, field="lower_rate_multiplier")
        upper = _finite(self.upper_rate_multiplier, field="upper_rate_multiplier")
        if lower < 0.0 or lower > 1.0 or upper < 1.0:
            raise ValueError(
                "Rate bounds must satisfy 0 <= lower_rate_multiplier <= 1 <= upper_rate_multiplier"
            )
        route_id = _optional_text(self.route_id)
        equation = _optional_text(self.reaction_equation)
        if route_id is None and equation is None:
            raise ValueError("Rate uncertainty requires route_id or reaction_equation")
        source = _optional_text(self.source)
        applicability = _optional_text(self.applicability)
        source_kind = _optional_text(self.source_kind)
        if source is None:
            raise ValueError("Rate uncertainty source is required")
        if applicability is None:
            raise ValueError("Rate uncertainty applicability is required")
        if source_kind not in _SOURCE_KINDS:
            raise ValueError(
                "source_kind must be one of " + ", ".join(sorted(_SOURCE_KINDS))
            )
        condition_id = self.condition_contract_id
        if source_kind == "method_calibration" and (
            not isinstance(condition_id, str) or not condition_id.strip()
        ):
            raise ValueError(
                "method_calibration requires an explicit nonempty string "
                "condition_contract_id; applicability prose does not bind a condition"
            )
        condition_id = _optional_text(condition_id)
        confidence = self.confidence_level
        if confidence is not None:
            confidence = _finite(confidence, field="confidence_level")
            if not 0.0 < confidence < 1.0:
                raise ValueError("confidence_level must be in (0, 1)")
        components = tuple(
            text for text in (_optional_text(item) for item in self.uncertainty_components)
            if text is not None
        )
        cost = self.estimated_cost
        units = _optional_text(self.cost_units)
        model = _optional_text(self.cost_model)
        if cost is not None:
            cost = _finite(cost, field="estimated_cost")
            if cost <= 0.0:
                raise ValueError("estimated_cost must be positive")
            if units is None or model is None:
                raise ValueError("estimated_cost requires cost_units and cost_model")
        elif units is not None or model is not None:
            raise ValueError("cost_units and cost_model require estimated_cost")
        object.__setattr__(self, "lower_rate_multiplier", lower)
        object.__setattr__(self, "upper_rate_multiplier", upper)
        object.__setattr__(self, "route_id", route_id)
        object.__setattr__(self, "reaction_equation", equation)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "source_kind", source_kind)
        object.__setattr__(self, "applicability", applicability)
        object.__setattr__(self, "condition_contract_id", condition_id)
        object.__setattr__(self, "confidence_level", confidence)
        object.__setattr__(self, "uncertainty_components", components)
        object.__setattr__(self, "estimated_cost", cost)
        object.__setattr__(self, "cost_units", units)
        object.__setattr__(self, "cost_model", model)

    @property
    def unbounded_downward(self) -> bool:
        return self.lower_rate_multiplier == 0.0

    @property
    def max_log_rate_deviation(self) -> float | None:
        """Worst finite |ln(k/k0)|; ``None`` means the lower bound is zero."""
        if self.unbounded_downward:
            return None
        return max(
            abs(math.log(self.lower_rate_multiplier)),
            abs(math.log(self.upper_rate_multiplier)),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RateUncertaintyRecord":
        if not isinstance(value, Mapping):
            raise TypeError("Each rate uncertainty record must be a mapping")
        required = (
            "lower_rate_multiplier", "upper_rate_multiplier", "source",
            "source_kind", "applicability",
        )
        missing = [field for field in required if field not in value]
        if missing:
            raise ValueError("Rate uncertainty record is missing: " + ", ".join(missing))
        return cls(
            lower_rate_multiplier=value["lower_rate_multiplier"],
            upper_rate_multiplier=value["upper_rate_multiplier"],
            source=value["source"],
            source_kind=value["source_kind"],
            applicability=value["applicability"],
            route_id=value.get("route_id"),
            reaction_equation=value.get("reaction_equation"),
            condition_contract_id=value.get("condition_contract_id"),
            confidence_level=value.get("confidence_level"),
            uncertainty_components=tuple(value.get("uncertainty_components") or ()),
            estimated_cost=value.get("estimated_cost"),
            cost_units=value.get("cost_units"),
            cost_model=value.get("cost_model"),
        )

    def as_dict(self) -> dict:
        return {
            "schema_version": 1,
            "kind": "rate_uncertainty_record",
            "route_id": self.route_id,
            "reaction_equation": self.reaction_equation,
            "lower_rate_multiplier": self.lower_rate_multiplier,
            "upper_rate_multiplier": self.upper_rate_multiplier,
            "max_log_rate_deviation": self.max_log_rate_deviation,
            "uncertainty_status": (
                "unbounded_downward" if self.unbounded_downward else "finite_bound"
            ),
            "source": self.source,
            "source_kind": self.source_kind,
            "applicability": self.applicability,
            "condition_contract_id": self.condition_contract_id,
            "confidence_level": self.confidence_level,
            "uncertainty_components": list(self.uncertainty_components),
            "estimated_cost": self.estimated_cost,
            "cost_units": self.cost_units,
            "cost_model": self.cost_model,
        }


def coerce_rate_uncertainty_evidence(
    value: Mapping[str, Any] | Iterable[Mapping[str, Any] | RateUncertaintyRecord] | None,
) -> tuple[RateUncertaintyRecord, ...]:
    """Validate an evidence artifact or iterable without filling missing data."""
    if value is None:
        return ()
    if isinstance(value, Mapping):
        kind = value.get("kind")
        if kind not in {None, "rate_uncertainty_evidence"}:
            raise ValueError("Rate uncertainty artifact kind must be 'rate_uncertainty_evidence'")
        raw_records = value.get("records")
        if raw_records is None:
            raise ValueError("Rate uncertainty artifact requires a records array")
    elif isinstance(value, (str, bytes)):
        raise TypeError("Rate uncertainty evidence must be a mapping or iterable of records")
    else:
        raw_records = value
    records = tuple(
        item if isinstance(item, RateUncertaintyRecord) else RateUncertaintyRecord.from_dict(item)
        for item in raw_records
    )
    seen_ids: set[tuple[str | None, str]] = set()
    seen_equations: set[tuple[str | None, str]] = set()
    for record in records:
        if record.route_id is not None:
            identity = (record.condition_contract_id, record.route_id)
            if identity in seen_ids:
                raise ValueError(f"Duplicate rate uncertainty route_id: {record.route_id}")
            seen_ids.add(identity)
        if record.route_id is None and record.reaction_equation is not None:
            equation_key = " ".join(record.reaction_equation.split())
            identity = (record.condition_contract_id, equation_key)
            if identity in seen_equations:
                raise ValueError(
                    f"Duplicate reaction-equation-only rate uncertainty: {record.reaction_equation}"
                )
            seen_equations.add(identity)
    return records


def rate_uncertainty_artifact(records: Iterable[RateUncertaintyRecord]) -> dict:
    serialized = [record.as_dict() for record in records]
    return {
        "schema_version": 1,
        "kind": "rate_uncertainty_evidence",
        "record_count": len(serialized),
        "records": serialized,
        "interpretation": (
            "Explicit multiplicative rate bounds used only for verification allocation; "
            "they are not propagated lifetime confidence intervals."
        ),
    }


def load_rate_uncertainty_evidence(path: Path) -> tuple[RateUncertaintyRecord, ...]:
    payload = json.loads(Path(path).read_text())
    return coerce_rate_uncertainty_evidence(payload)


def uncertainty_from_forward_rate_upper_bound(
    route: Mapping[str, Any],
    path_result: Mapping[str, Any],
    condition_contract: Mapping[str, Any],
) -> RateUncertaintyRecord | None:
    """Convert a retained physical rate ceiling into an explicitly one-sided bound."""
    verification = path_result.get("route_verification") or path_result
    bound = verification.get("forward_loss_upper_bound") or {}
    if not bound.get("applicable"):
        return None
    baseline = route.get("retained_rate_evidence") or {}
    try:
        upper_rate = float(bound.get("forward_rate_upper_bound_si"))
        baseline_rate = float(baseline.get("rate_coefficient_si"))
    except (TypeError, ValueError):
        return None
    upper_units = str(bound.get("forward_rate_upper_bound_units") or "")
    baseline_units = str(baseline.get("rate_units") or "")
    if (
        not math.isfinite(upper_rate) or upper_rate < 0.0
        or not math.isfinite(baseline_rate) or baseline_rate <= 0.0
        or not upper_units or upper_units != baseline_units
    ):
        return None
    assumptions = bound.get("assumptions") or ()
    return RateUncertaintyRecord(
        route_id=str(route.get("route_id") or "") or None,
        reaction_equation=str(
            route.get("reaction_equation") or route.get("printed_reaction_equation") or ""
        ) or None,
        lower_rate_multiplier=0.0,
        upper_rate_multiplier=max(1.0, upper_rate / baseline_rate),
        source=str(bound.get("bound_basis") or "condition_specific_forward_rate_upper_bound"),
        source_kind="physical_bound",
        applicability=(
            "Exact retained condition; " + "; ".join(str(item) for item in assumptions)
        ),
        condition_contract_id=(
            str(condition_contract.get("contract_id"))
            if condition_contract.get("contract_id") else None
        ),
        uncertainty_components=(
            "rate_theory", "capture_probability", "recrossing",
            "pressure_dependence", "mechanism_applicability",
        ),
    )


__all__ = [
    "RateUncertaintyRecord",
    "coerce_rate_uncertainty_evidence",
    "load_rate_uncertainty_evidence",
    "rate_uncertainty_artifact",
    "uncertainty_from_forward_rate_upper_bound",
]
