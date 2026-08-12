import unittest

from storca.environment_local_modes import plan_environment_local_mode_fallbacks


class EnvironmentLocalModePlanningTests(unittest.TestCase):
    @staticmethod
    def conformers():
        output = []
        for index in range(4):
            bonds = [{
                "bond_class": "O-H", "spectral_band_class": "hydrogen_bonded_oh",
                "coordination_class": "donor", "molecule_index": 0,
                "heavy_atom": 0, "hydrogen_atom": 1,
            }]
            if index < 3:
                bonds.append({
                    "bond_class": "O-H", "spectral_band_class": "non_donating_oh",
                    "coordination_class": "acceptor", "molecule_index": 1,
                    "heavy_atom": 3, "hydrogen_atom": 4,
                })
            output.append({
                "independent_environment_id": f"env-{index}",
                "optimized_xyz": f"env-{index}.xyz", "frequency_output": f"env-{index}.out",
                "snapshot_hessian_reliability": {"stationarity_status": "poor"},
                "local_stretch_bonds": bonds,
            })
        return output

    def test_planner_selects_three_distinct_representatives_per_class(self):
        plan = plan_environment_local_mode_fallbacks(
            self.conformers(), maximum_orca_invocations=24,
        )
        self.assertEqual(plan["status"], "planned")
        self.assertEqual(plan["planned_orca_invocations"], 24)
        for detail in plan["classes"].values():
            self.assertEqual(len(detail["selected_representative_positions"]), 3)

    def test_planner_never_starts_a_class_without_complete_coverage_budget(self):
        plan = plan_environment_local_mode_fallbacks(
            self.conformers(), maximum_orca_invocations=12,
        )
        selected = [item for item in plan["classes"].values() if item["status"] == "selected"]
        rejected = [item for item in plan["classes"].values() if item["status"] != "selected"]
        self.assertEqual(len(selected), 1)
        self.assertEqual(len(rejected), 1)
        self.assertEqual(plan["planned_orca_invocations"], 12)


if __name__ == "__main__":
    unittest.main()
