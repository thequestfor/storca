import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from storca.route_verify import RouteSpec
from storca.verification_adapters import make_arkane_rate_builder, make_orca_route_verifier
from storca.verification_engine import _normalize_replacement


def _condition() -> dict:
    return {
        "temperature_K": 310.0,
        "pressure_bar": 0.8,
        "target_label": "stability",
        "composition": {
            "stability": 0.01,
            "nitrogen": 0.74,
            "oxygen": 0.20,
            "water": 0.05,
        },
    }


def _profile() -> dict:
    return {
        "method": "r2SCAN-3c", "basis": "composite def2-mTZVPP",
        "orca_keywords": ["r2SCAN-3c"], "harmonic_scale_factor": 0.99,
    }


class VerificationAdapterTests(unittest.TestCase):
    def test_default_orca_adapter_passes_route_and_execution_contract(self):
        normalized = RouteSpec(
            route_id="route-1", source="rmg", parent_smiles="[H]",
            reactant_smiles=("[H]",), product_smiles=("[H]",), multiplicity=2,
        )
        path_result = {
            "route_verification": {
                "route_id": "route-1", "status": "verified_transition_state_path",
                "path_classification": "verified_barriered_path",
            }
        }
        with tempfile.TemporaryDirectory() as temporary, \
                patch("storca.verification_adapters.route_spec_from_rmg_route", return_value=normalized) as normalize, \
                patch("storca.verification_adapters.run_generic_reaction_path", return_value=path_result) as run, \
                patch("storca.verification_adapters._validated_trajectory_visual", return_value=None):
            callback = make_orca_route_verifier(
                "[H]", ncores=4, method_profile=_profile(), timeout_seconds=12.0,
                orientations=2, neb_images=6,
            )
            result = callback({"route_id": "route-1", "reaction_equation": "A<=>B"}, Path(temporary))
        self.assertIs(result, path_result)
        normalize.assert_called_once()
        self.assertEqual(normalize.call_args.args[0], "[H]")
        self.assertIs(run.call_args.args[0], normalized)
        self.assertEqual(run.call_args.kwargs, {
            "ncores": 4, "method_keywords": ["r2SCAN-3c"], "timeout_seconds": 12.0,
            "orientations": 2, "neb_images": 6, "condition_contract": None,
        })

    def test_arkane_adapter_preserves_repeated_species_spin_collision_bath_and_nested_result(self):
        route = {
            "route_id": "recombination",
            "reaction_equation": "2 hydrogen(7) <=> hydrogen2(8)",
            "resolved_endpoints": {
                "reactants": [
                    {"label": "hydrogen(7)", "smiles": "[H]"},
                    {"label": "hydrogen(7)", "smiles": "[H]"},
                ],
                "products": [{"label": "hydrogen2(8)", "smiles": "[H][H]"}],
            },
            "retained_rate_evidence": {
                "collision_limit_si": 4.2,
                "collision_limit_units": "m^3/(mol*s)",
                "collision_limit_source": "RMG retained transport collision bound",
                "evaluated_at": {"temperature_K": 310.0, "pressure_bar": 0.8},
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            hydrogen = folder / "hydrogen.out"
            product = folder / "hydrogen2.out"
            ts = folder / "ts.out"
            path_result = {
                "route_verification": {
                    "route_id": "recombination",
                    "status": "verified_transition_state_path",
                    "selected_multiplicity": 3,
                    "arkane_inputs": {
                        "reactants": [
                            {"label": "hydrogen(7)", "smiles": "[H]", "frequency_output": str(hydrogen),
                             "multiplicity": 2, "status": "validated_minimum"},
                            {"label": "hydrogen(7)", "smiles": "[H]", "frequency_output": str(hydrogen),
                             "multiplicity": 2, "status": "validated_minimum"},
                        ],
                        "products": [
                            {"label": "hydrogen2(8)", "smiles": "[H][H]", "frequency_output": str(product),
                             "multiplicity": 1, "status": "validated_minimum"},
                        ],
                        "transition_state": {"frequency_output": str(ts)},
                    },
                }
            }
            library = folder / "verified-library" / "kinetics"
            library.mkdir(parents=True)
            nested = {
                "verification_status": "verified_arkane_pdep",
                "route_id": "recombination",
                "reaction_equation": route["reaction_equation"],
                "molecularity": 2,
                "condition_rate": {"value": 1.1, "units": "m^3/(mol*s)",
                                   "temperature_K": 310.0, "pressure_bar": 0.8},
                "collision_limit": {"value": 4.2, "units": "m^3/(mol*s)",
                                    "source": "RMG retained transport collision bound"},
                "library": str(library),
            }
            captured = {}

            def fake_run(spec, workdir, *, rmg_env):
                captured.update(spec=spec, workdir=workdir, rmg_env=rmg_env)
                return {"status": "completed", "generated_kinetics_library": nested}

            with patch("storca.verification_adapters.run_verified_arkane_route", side_effect=fake_run):
                callback = make_arkane_rate_builder(
                    _condition(), method_profile=_profile(), rmg_env="rmg4_env",
                )
                result = callback(route, path_result, folder / "rate")

            spec = captured["spec"]
            self.assertEqual(spec.reactant_labels, ("hydrogen", "hydrogen"))
            self.assertEqual(spec.reactant_orca_outputs, (hydrogen, hydrogen))
            self.assertEqual(spec.transition_state_multiplicity, 3)
            self.assertEqual(spec.reaction_equation, route["reaction_equation"])
            self.assertEqual(spec.collision_limit_m3_mol_s, 4.2)
            self.assertEqual(spec.collision_limit_source, "RMG retained transport collision bound")
            self.assertEqual(dict(spec.bath_gas), {"nitrogen": 0.74 / 0.99, "oxygen": 0.20 / 0.99,
                                                   "water": 0.05 / 0.99})
            self.assertIn(310.0, spec.temperatures_K)
            self.assertIn(0.8, spec.pressures_bar)
            self.assertEqual(captured["rmg_env"], "rmg4_env")

            replacement, paths = _normalize_replacement(result, route)
            self.assertEqual(replacement, nested)
            self.assertEqual(paths, [library])

    def test_arkane_adapter_rejects_collision_bound_from_wrong_condition(self):
        route = {
            "route_id": "association", "reaction_equation": "A+B<=>C",
            "resolved_endpoints": {
                "reactants": [{"label": "A", "smiles": "[H]"}, {"label": "B", "smiles": "[OH]"}],
                "products": [{"label": "C", "smiles": "O"}],
            },
            "retained_rate_evidence": {
                "collision_limit_si": 1.0, "collision_limit_units": "m^3/(mol*s)",
                "collision_limit_source": "retained bound",
                "evaluated_at": {"temperature_K": 298.15, "pressure_bar": 1.0},
            },
        }
        path = {
            "route_verification": {
                "status": "verified_transition_state_path", "selected_multiplicity": 1,
                "arkane_inputs": {
                    "reactants": [{"label": "A", "smiles": "[H]", "frequency_output": "a", "multiplicity": 2,
                                   "status": "validated_minimum"},
                                  {"label": "B", "smiles": "[OH]", "frequency_output": "b", "multiplicity": 2,
                                   "status": "validated_minimum"}],
                    "products": [{"label": "C", "smiles": "O", "frequency_output": "c", "multiplicity": 1,
                                  "status": "validated_minimum"}],
                    "transition_state": {"frequency_output": "ts"},
                },
            }
        }
        callback = make_arkane_rate_builder(_condition(), method_profile=_profile(), rmg_env="rmg_env")
        with self.assertRaisesRegex(ValueError, "does not match the immutable condition"):
            callback(route, path, Path("unused"))


if __name__ == "__main__":
    unittest.main()
