from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

from storca.route_thermochemistry import (
    assemble_route_thermochemistry,
    condition_specific_forward_loss_bound,
    parse_orca_thermochemistry,
)
from storca.route_verify import RouteSpec


def _orca_frequency(path: Path, *, electronic: float, enthalpy: float, gibbs: float,
                    temperature: float = 298.15, normal: bool = True) -> Path:
    path.write_text(
        f"FINAL SINGLE POINT ENERGY      {electronic:.12f}\n"
        f"THERMOCHEMISTRY AT {temperature:.2f}K\n"
        f"Temperature         ...   {temperature:.2f} K\n"
        "Zero point energy                ...      0.01000000 Eh       6.28 kcal/mol\n"
        f"Total Enthalpy                    ...   {enthalpy:.12f} Eh\n"
        f"Final Gibbs free energy         ...   {gibbs:.12f} Eh\n"
        + ("ORCA TERMINATED NORMALLY\n" if normal else "ORCA finished with error\n")
    )
    return path


class RouteThermochemistryTests(unittest.TestCase):
    @staticmethod
    def _route() -> RouteSpec:
        return RouteSpec(
            route_id="thermo-route", source="test", parent_smiles="NNO",
            reactant_smiles=("NNO", "[O][O]"), product_smiles=("NN[O]", "[O]O"),
            reactant_labels=("stability(1)", "oxygen(2)"),
            product_labels=("NN[O](3)", "[O]O(4)"),
            reactant_charges=(0, 0), product_charges=(0, 0),
            multiplicity=3, reactant_multiplicities=(1, 3), product_multiplicities=(2, 2),
        )

    @staticmethod
    def _condition() -> dict:
        return {
            "phase_model": "homogeneous gas-phase surrogate",
            "temperature_K": 298.15, "pressure_bar": 1.0,
            "composition": {"stability": 0.01, "oxygen": 0.20, "nitrogen": 0.79},
            "target_label": "stability", "target_duration_seconds": 31536000.0,
            "retention_fraction": 0.95,
        }

    def test_parser_requires_complete_normal_frequency_thermochemistry(self):
        with tempfile.TemporaryDirectory() as temp:
            valid_path = _orca_frequency(
                Path(temp) / "valid.out", electronic=-100.0, enthalpy=-99.99, gibbs=-100.02,
            )
            invalid_path = _orca_frequency(
                Path(temp) / "invalid.out", electronic=-100.0, enthalpy=-99.99,
                gibbs=-100.02, normal=False,
            )
            valid = parse_orca_thermochemistry(valid_path)
            invalid = parse_orca_thermochemistry(invalid_path)
        self.assertTrue(valid["valid"])
        self.assertEqual(valid["gibbs_free_energy_hartree"], -100.02)
        self.assertFalse(invalid["valid"])

    def test_strongly_endergonic_two_to_two_route_gets_below_loss_upper_bound(self):
        route = self._route()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            files = [
                _orca_frequency(root / "r1.out", electronic=-100.0, enthalpy=-99.99, gibbs=-100.02),
                _orca_frequency(root / "r2.out", electronic=-200.0, enthalpy=-199.99, gibbs=-200.02),
                _orca_frequency(root / "p1.out", electronic=-99.95, enthalpy=-99.94, gibbs=-99.97),
                _orca_frequency(root / "p2.out", electronic=-199.95, enthalpy=-199.94, gibbs=-199.97),
            ]
            stationary = {
                "valid": True,
                "reactants": [
                    {"label": "stability(1)", "smiles": "NNO", "frequency_output": str(files[0])},
                    {"label": "oxygen(2)", "smiles": "[O][O]", "frequency_output": str(files[1])},
                ],
                "products": [
                    {"label": "NN[O](3)", "smiles": "NN[O]", "frequency_output": str(files[2])},
                    {"label": "[O]O(4)", "smiles": "[O]O", "frequency_output": str(files[3])},
                ],
            }
            thermo = assemble_route_thermochemistry(route, stationary, self._condition())
            bound = condition_specific_forward_loss_bound(route, thermo, self._condition())
        self.assertTrue(thermo["valid"])
        self.assertAlmostEqual(thermo["reaction"]["delta_g_hartree"], 0.1)
        self.assertTrue(bound["applicable"])
        self.assertEqual(bound["status"], "forward_loss_below_retention_threshold_upper_bound")
        self.assertLess(bound["target_loss_fraction_upper_bound"], 0.05)
        self.assertEqual(bound["co_reactant_concentration_source"], "declared_condition_composition")

    def test_condition_temperature_mismatch_fails_closed(self):
        route = self._route()
        with tempfile.TemporaryDirectory() as temp:
            output = _orca_frequency(
                Path(temp) / "freq.out", electronic=-10.0, enthalpy=-9.9, gibbs=-10.1,
                temperature=350.0,
            )
            stationary = {
                "reactants": [{"frequency_output": str(output)}, {"frequency_output": str(output)}],
                "products": [{"frequency_output": str(output)}, {"frequency_output": str(output)}],
            }
            thermo = assemble_route_thermochemistry(route, stationary, self._condition())
        self.assertFalse(thermo["valid"])
        self.assertEqual(thermo["status"], "incomplete_or_condition_mismatch")


if __name__ == "__main__":
    unittest.main()
