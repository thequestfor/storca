from __future__ import annotations

import math

from storca.rmg_bridge import (
    _flux_closure,
    _integrate_flux_snapshots,
    _logarithmic_times,
)


def _snapshot(rate: float) -> list[dict]:
    return [{
        "destruction_kmol_s": rate,
        "production_kmol_s": 0.0,
        "forward_extent_kmol_s": rate,
        "reverse_extent_kmol_s": 0.0,
    }]


def test_logarithmic_grid_honors_fast_chemistry_floor() -> None:
    times = _logarithmic_times(31_536_000.0, 320, 1.0e-10)
    assert times[0] == 0.0
    assert math.isclose(times[1], 1.0e-10, rel_tol=1e-12)
    assert times[-1] == 31_536_000.0
    assert all(right > left for left, right in zip(times, times[1:]))


def test_refinement_grid_can_resolve_a_crossing_before_the_coarse_first_point() -> None:
    horizon = 31_536_000.0
    coarse_first = horizon * 1e-12
    coarse_crossing = 1.6e-6
    refined_first = min(coarse_first * 1e-3, coarse_crossing * 1e-4)
    assert refined_first < coarse_crossing * 1e-2


def test_refined_fast_decay_flux_closes_against_inventory() -> None:
    initial = 1.0
    rate_constant = 1.0e6
    horizon = 1.0
    times = _logarithmic_times(horizon, 1200, 1.0e-12)
    profile = [
        {"time_seconds": time, "target_amount_kmol": math.exp(-rate_constant * time)}
        for time in times
    ]
    snapshots = [_snapshot(rate_constant * item["target_amount_kmol"]) for item in profile]
    totals = _integrate_flux_snapshots(profile, snapshots, 1)
    _, _, _, relative_error = _flux_closure(
        initial,
        profile[-1]["target_amount_kmol"],
        totals["destruction_kmol_s"],
        totals["production_kmol_s"],
    )
    assert relative_error < 0.01
