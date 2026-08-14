import tempfile
import unittest
from unittest.mock import patch
import json
from pathlib import Path

from src.inputgen import (create_orca_environment_refinement_input,
                          create_orca_gradient_input)
from storca.environment_refinement import (parse_orca_engrad,
                                           refine_selected_orca_environments,
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
        self.assertIn("{ A 0 1 2 C }", refinement_text)
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

    def test_selected_refinement_updates_execution_geometry_and_gate(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            clusters = folder / "clusters"
            clusters.mkdir()
            source = folder / "source.xyz"
            source.write_text("3\ncontact\nO 0 0 0\nH 1 0 0\nO 2.8 0 0\n")
            manifest = clusters / "selected-conformers.json"
            manifest.write_text(json.dumps({"conformers": [{
                "selected_position": 1,
                "source_xtb_candidate_id": "candidate-1",
                "source_xtb_xyz": str(source),
                "xyz": str(source),
                "hydrogen_bond_interactions": [{
                    "donor_atom": 0, "donor_hydrogen": 1, "acceptor_atom": 2,
                }],
            }]}))
            refined = folder / "refined.xyz"
            refined.write_text(source.read_text())
            record = {
                "status": "completed", "refined_xyz": str(refined),
                "stationarity_status": "usable", "full_hessian_use": "permitted",
                "gradient": {
                    "gradient_rms_hartree_per_bohr": 1e-4,
                    "gradient_maximum_component_hartree_per_bohr": 2e-4,
                },
            }
            with patch(
                "storca.environment_refinement.refine_orca_environment",
                return_value=record,
            ):
                paths, report = refine_selected_orca_environments(folder)
            retained = json.loads(manifest.read_text())
        self.assertEqual(paths, [refined])
        self.assertTrue(report.name == "environment-refinement.json")
        self.assertEqual(retained["conformers"][0]["xyz"], str(refined))
        self.assertEqual(
            retained["conformers"][0]["vibrational_route"]["route"],
            "full_orca_hessian",
        )


if __name__ == "__main__":
    unittest.main()
