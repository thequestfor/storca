from __future__ import annotations

import io
import json
import math
import sys
import types
import unittest
from unittest.mock import patch

from storca import rmg_bridge
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
        "net_target_loss_kmol_s": rate,
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
        totals["net_target_loss_kmol_s"],
    )
    assert relative_error < 0.01


def test_direct_net_integration_avoids_gross_reversible_cancellation() -> None:
    profile = [
        {"time_seconds": 0.0, "target_amount_kmol": 1.0},
        {"time_seconds": 1.0, "target_amount_kmol": 1.0 - 1.0e-9},
    ]
    snapshots = [[{
        "destruction_kmol_s": 1.0e8,
        "production_kmol_s": 1.0e8 - 1.0e-9,
        "forward_extent_kmol_s": 1.0e8,
        "reverse_extent_kmol_s": 1.0e8 - 1.0e-9,
        "net_target_loss_kmol_s": 1.0e-9,
    }]] * 2
    totals = _integrate_flux_snapshots(profile, snapshots, 1)
    _, direct_net, _, relative_error = _flux_closure(
        1.0,
        1.0 - 1.0e-9,
        totals["destruction_kmol_s"],
        totals["production_kmol_s"],
        totals["net_target_loss_kmol_s"],
    )
    assert math.isclose(direct_net, 1.0e-9, rel_tol=1e-12)
    assert relative_error < 1.0e-6


class TestInspectionNumericalAdmission(unittest.TestCase):
    def _inspect(self, forward, *, reverse=0.5, collision=1.0e8):
        # Deliberately independent of the installed owning RMG environment.
        # Pmin=0 mirrors the invalid low-pressure Chemkin export found live.
        class Chebyshev:
            Pmin = types.SimpleNamespace(value_si=0.0)
            comment = None

        def species(label):
            return types.SimpleNamespace(label=label, molecule=[], multiplicity=1, reactive=True)

        reactants = [species("A"), species("B")]
        products = [species("C"), species("D")]
        reaction = types.SimpleNamespace(
            reactants=reactants, products=products, reversible=True,
            kinetics=Chebyshev(), library="retained-library",
            get_rate_coefficient=lambda *args: forward,
            generate_reverse_rate_coefficient=lambda: types.SimpleNamespace(
                get_rate_coefficient=lambda *args: reverse,
            ),
            calculate_coll_limit=lambda **kwargs: collision,
        )
        module = types.ModuleType("rmgpy.chemkin")
        module.load_chemkin_file = lambda *args, **kwargs: (reactants + products, [reaction])
        module.get_species_identifier = lambda item: item.label
        with patch.dict(sys.modules, {"rmgpy": types.ModuleType("rmgpy"), "rmgpy.chemkin": module}):
            return rmg_bridge.inspect({
                "chemkin": "retained-chemkin.inp", "temperature_K": 298.15,
                "pressure_bar": 1.0e-6,
            })

    def test_nonfinite_or_negative_forward_rate_fails_with_reaction_context(self):
        for value in (float("nan"), float("inf"), -float("inf"), -1.0, None):
            with self.subTest(value=value):
                with self.assertRaisesRegex(
                    ValueError, r"Reaction index 0 \(A\+B<=>C\+D\), kinetics Chebyshev: forward rate",
                ):
                    self._inspect(value)

    def test_finite_and_zero_forward_rates_remain_valid_numerical_evidence(self):
        for value in (0.0, 3.235611525078417e-7):
            with self.subTest(value=value):
                result = self._inspect(value, reverse=0.0)
                self.assertEqual(result["status"], "completed")
                self.assertEqual(result["reactions"][0]["evaluated_forward_rate_coefficient_si"], value)
                self.assertEqual(result["reactions"][0]["evaluated_reverse_rate_coefficient_si"], 0.0)
                json.dumps(result, allow_nan=False)

    def test_invalid_optional_reverse_and_collision_values_are_unavailable(self):
        for value in (float("nan"), float("inf"), -1.0, None):
            with self.subTest(value=value):
                result = self._inspect(1.0, reverse=value, collision=value)
                record = result["reactions"][0]
                self.assertIsNone(record["evaluated_reverse_rate_coefficient_si"])
                self.assertEqual(record["reverse_rate_unavailable_reason"],
                                 "reverse_rate_nonfinite_negative_or_unavailable")
                self.assertIsNone(record["forward_collision_limit_si"])
                self.assertIsNone(record["reverse_collision_limit_si"])
                self.assertIsNone(record["collision_limit_source"])
                json.dumps(result, allow_nan=False)

    def test_zero_collision_estimate_is_not_a_physical_bound(self):
        record = self._inspect(1.0, collision=0.0)["reactions"][0]
        self.assertIsNone(record["forward_collision_limit_si"])
        self.assertIsNone(record["reverse_collision_limit_si"])

    def test_main_withholds_nonfinite_json_payload_for_either_operation(self):
        for operation in ("inspect", "propagate"):
            for invalid in (float("nan"), float("inf"), -float("inf")):
                with self.subTest(operation=operation, invalid=invalid):
                    result = {"status": "completed", "nested": {"unvalidated": invalid}}
                    output = io.StringIO()
                    with patch.object(rmg_bridge, operation, return_value=result), \
                            patch.object(sys, "argv", ["rmg_bridge.py"]), \
                            patch.object(sys, "stdin", io.StringIO(json.dumps({"operation": operation}))), \
                            patch.object(sys, "stdout", output):
                        rmg_bridge.main()
                    payload = json.loads(output.getvalue())
                    self.assertEqual(payload["status"], "failed")
                    self.assertIn("Out of range float values", payload["failure_reason"])
                    self.assertNotIn("nested", payload)
                    self.assertNotIn(": NaN", output.getvalue())
                    self.assertNotIn(": Infinity", output.getvalue())

    def test_main_reports_contextual_inspection_failure_as_valid_json(self):
        output = io.StringIO()
        with patch.object(rmg_bridge, "inspect", side_effect=ValueError(
            "Reaction index 0 (A+B<=>C+D), kinetics Chebyshev: non-finite forward rate",
        )), patch.object(sys, "argv", ["rmg_bridge.py"]), \
                patch.object(sys, "stdin", io.StringIO('{"operation": "inspect"}')), \
                patch.object(sys, "stdout", output):
            rmg_bridge.main()
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["status"], "failed")
        self.assertIn("Reaction index 0", payload["failure_reason"])
        self.assertIn("Chebyshev", payload["failure_reason"])
