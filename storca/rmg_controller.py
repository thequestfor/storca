"""Bounded, restartable RMG mechanism-search policy."""

from __future__ import annotations

from dataclasses import dataclass

from src.orca_runner import rmg_walltime_seconds


@dataclass(frozen=True)
class RmgBudget:
    attempt: int
    walltime: str
    max_processes: int
    max_iterations: int
    max_edge_species: int

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def _walltime(seconds: float) -> str:
    seconds = int(seconds)
    days, seconds = divmod(seconds, 86400)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return f"{days:02d}:{hours:02d}:{minutes:02d}:{seconds:02d}"


def staged_budgets(resources: dict, attempts: int) -> list[RmgBudget]:
    """Create bounded restarts; each later attempt grows only declared limits."""
    if attempts < 1:
        raise ValueError("RMG attempts must be at least one")
    base_seconds = rmg_walltime_seconds(resources["walltime"])
    if base_seconds is None:
        raise ValueError("Staged RMG control requires an explicit walltime")
    budgets = []
    for index in range(attempts):
        factor = 2 ** index
        budgets.append(RmgBudget(
            attempt=index + 1, walltime=_walltime(base_seconds * factor),
            max_processes=int(resources["max_processes"]),
            max_iterations=int(resources["max_iterations"]) * factor,
            max_edge_species=int(resources["max_edge_species"]) * factor,
        ))
    return budgets


def restart_is_justified(evidence: dict) -> bool:
    """Restart only genuine resource exhaustion with a retained seed."""
    if evidence.get("status") != "incomplete":
        return False
    execution = evidence.get("execution_assessment") or {}
    reason = execution.get("termination_reason")
    seed = evidence.get("artifacts", {}).get("seed")
    return reason == "resource_or_interrupt_termination" and bool(seed)


def convergence_comparison(attempts: list[dict]) -> dict:
    """Compare retained attempts without asserting convergence by omission."""
    if len(attempts) < 2:
        return {"status": "not_comparable", "reason": "fewer_than_two_attempts"}
    previous, current = attempts[-2], attempts[-1]
    previous_loss = ((previous.get("independent_cantera_propagation") or previous.get("solver_profile") or {})
                     .get("target_loss_fraction"))
    current_loss = ((current.get("independent_cantera_propagation") or current.get("solver_profile") or {})
                    .get("target_loss_fraction"))
    previous_routes = {route.get("reaction_equation") for route in previous.get("candidate_routes", [])}
    current_routes = {route.get("reaction_equation") for route in current.get("candidate_routes", [])}
    return {
        "status": "compared",
        "target_loss_change": None if previous_loss is None or current_loss is None else abs(current_loss - previous_loss),
        "new_candidate_route_count": len(current_routes - previous_routes),
        "previous_candidate_route_count": len(previous_routes),
        "current_candidate_route_count": len(current_routes),
        "interpretation": "Agreement is diagnostic only; RMG normal completion and requested-time coverage remain required for convergence.",
    }
