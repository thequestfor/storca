import tempfile
import unittest
from pathlib import Path

from src.inputgen import (create_orca_environment_refinement_input,
                          create_orca_gradient_input)
from storca.environment_refinement import (parse_orca_engrad,
                                           representative_vibrational_route)


class EnvironmentRefinementTests(unittest.TestCase):
    def test_refinement_constraints_are_absent_from_gradient_input(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            xyz = folder / "input.xyz"
            xyz.write_text("3\nwater contact\nO 0 0 0\nH 1 0 0\nO 2.8 0 0\n")
            refinement = create_orca_environment_refinement_input(
                xyz,
                interactions=[{
                    "donor_atom": 0, "donor_hydrogen": 1, "acceptor_atom": 2,
                    "h_bond_distance_angstrom": 1.8,
                    "donor_h_acceptor_angle_degrees": 180.0,
                }],
                charge=0, multiplicity=1,
            )
            gradient = create_orca_gradient_input(
                xyz, charge=0, multiplicity=1, label="unrestrained-gradient",
            )
            refinement_text, gradient_text = refinement.read_text(), gradient.read_text()
        self.assertIn("{ B 1 2 1.80000000 C }", refinement_text)
        self.assertIn("{ A 0 1 2 180.00000000 C }", refinement_text)
        self.assertIn("EnGrad", gradient_text)
        self.assertNotIn("Constraints", gradient_text)

    def test_engrad_parser_reports_rms_and_maximum_component(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "gradient.engrad"
            path.write_text("""# Number of atoms
2
# energy
-10.0
# gradient
1.0D-4
-1.0D-4
2.0D-4
-2.0D-4
0.0
0.0
# coordinates follow
1 0 0 0
1 1 0 0
""")
            parsed = parse_orca_engrad(path)
        self.assertEqual(parsed["atom_count"], 2)
        self.assertAlmostEqual(parsed["gradient_maximum_component_hartree_per_bohr"], 2e-4)
        self.assertGreater(parsed["gradient_rms_hartree_per_bohr"], 0.0)

    def test_failed_stationarity_gate_routes_to_local_mode_finite_differences(self):
        route = representative_vibrational_route({
            "status": "completed", "stationarity_status": "poor",
            "full_hessian_use": "diagnostic_only",
        })
        self.assertEqual(route["route"], "orca_projected_local_mode_finite_differences")
        self.assertEqual(route["full_hessian_use"], "not_permitted")

    def test_passed_stationarity_gate_routes_to_full_hessian(self):
        route = representative_vibrational_route({
            "status": "completed", "stationarity_status": "usable",
            "full_hessian_use": "permitted",
        })
        self.assertEqual(route["route"], "full_orca_hessian")


if __name__ == "__main__":
    unittest.main()
