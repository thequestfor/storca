import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from storca.reaction_candidates import ReactionCandidate
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


def _ordinary_route_and_path(folder: Path) -> tuple[dict, dict]:
    route = {
        "route_id": "ordinary-callback-fixture", "multiplicity": 1,
        "reaction_equation": "cyclobutene <=> butadiene",
        "environment": {"model": "homogeneous_gas"},
        "resolved_endpoints": {
            "reactants": [{"label": "cyclobutene", "smiles": "C1=CCC1", "multiplicity": 1}],
            "products": [{"label": "butadiene", "smiles": "C=CC=C", "multiplicity": 1}],
        },
    }
    inputs = {"transition_state": {"frequency_output": str(folder / "ts.out")}}
    for side, records in route["resolved_endpoints"].items():
        inputs[side] = [{**record, "status": "validated_minimum",
                         "frequency_output": str(folder / (record["label"] + ".out"))}
                        for record in records]
    path = {"route_verification": {
        "status": "verified_transition_state_path", "path_classification": "verified_barriered_path",
        "selected_multiplicity": 1, "arkane_inputs": inputs,
    }}
    return route, path


class VerificationAdapterTests(unittest.TestCase):
    def test_arkane_callback_rejects_solution_as_gas_before_execution(self):
        with tempfile.TemporaryDirectory() as temporary, patch(
            "storca.verification_adapters.run_verified_arkane_route", return_value={"status": "completed"},
        ) as run:
            folder = Path(temporary)
            route, path = _ordinary_route_and_path(folder)
            route["environment"] = {"model": "implicit_solution", "solvent": "hexane", "solvent_model": "SMD"}
            condition = {**_condition(), "phase_model": "homogeneous implicit solution",
                         "solvent": "hexane", "solvent_model": "SMD"}
            callback = make_arkane_rate_builder(condition, method_profile=_profile(), rmg_env="rmg_env")
            with self.assertRaisesRegex(ValueError, "rate_theory_environment_mismatch"):
                callback(route, path, folder / "rate")
            run.assert_not_called()

    def test_arkane_callback_rejects_contradictory_barrierless_classification(self):
        with tempfile.TemporaryDirectory() as temporary, patch(
            "storca.verification_adapters.run_verified_arkane_route", return_value={"status": "completed"},
        ) as run:
            folder = Path(temporary)
            route, path = _ordinary_route_and_path(folder)
            path["route_verification"]["path_classification"] = "barrierless_capture_candidate"
            callback = make_arkane_rate_builder(_condition(), method_profile=_profile(), rmg_env="rmg_env")
            with self.assertRaisesRegex(ValueError, "verified_barriered_path_required"):
                callback(route, path, folder / "rate")
            run.assert_not_called()

    def test_arkane_callback_retains_matching_evidence_without_rewriting_library(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            route, path = _ordinary_route_and_path(folder)
            library = folder / "old-library"
            library.mkdir()
            manifest = library / "generated-kinetics.json"
            retained_bytes = b'{"fixture":"retained old library; no computed rate"}\n'
            manifest.write_bytes(retained_bytes)
            nested = {"library": str(library), "verification_status": "verified_arkane_tst",
                      "rate_uncertainty": {"status": "calibration_required"}}
            native_result = {"status": "completed", "generated_kinetics_library": nested}
            callback = make_arkane_rate_builder(_condition(), method_profile=_profile(), rmg_env="rmg_env")
            with patch("storca.verification_adapters.run_verified_arkane_route", return_value=native_result) as run:
                first = callback(route, path, folder / "rate")
                decision_path = folder / "rate" / "method-applicability.json"
                first_bytes, first_mtime = decision_path.read_bytes(), decision_path.stat().st_mtime_ns
                second = callback(route, path, folder / "rate")
                self.assertEqual(run.call_count, 2)
            self.assertEqual(decision_path.read_bytes(), first_bytes)
            self.assertEqual(decision_path.stat().st_mtime_ns, first_mtime)
            self.assertEqual(manifest.read_bytes(), retained_bytes)
            self.assertEqual(first["method_applicability_evidence"], second["method_applicability_evidence"])
            self.assertEqual(first["method_applicability_evidence"]["sha256"], hashlib.sha256(first_bytes).hexdigest())
            decision = json.loads(first_bytes)
            self.assertEqual(decision, first["method_applicability"])
            self.assertEqual(decision["tunneling_model"], run.call_args.args[0].tunneling.casefold())
            self.assertEqual(decision["pressure_dependent"], run.call_args.args[0].use_pressure_dependence)
            self.assertTrue(decision["central_rate_allowed"])
            self.assertFalse(decision["validated_prediction_allowed"])
            self.assertEqual(first["generated_kinetics_library"]["rate_uncertainty"], nested["rate_uncertainty"])
            self.assertNotIn("method_applicability", native_result)
            self.assertNotIn("method_applicability", nested)

    def test_arkane_callback_retains_rejection_and_refuses_corrupt_or_changed_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            route, path = _ordinary_route_and_path(folder)
            path["route_verification"]["path_classification"] = "barrierless_capture_candidate"
            callback = make_arkane_rate_builder(_condition(), method_profile=_profile(), rmg_env="rmg_env")
            with patch("storca.verification_adapters.run_verified_arkane_route") as run:
                with self.assertRaisesRegex(ValueError, "verified_barriered_path_required"):
                    callback(route, path, folder / "rate")
                decision_path = folder / "rate" / "method-applicability.json"
                rejected_bytes = decision_path.read_bytes()
                self.assertFalse(json.loads(rejected_bytes)["central_rate_allowed"])
                path["route_verification"]["path_classification"] = "verified_barriered_path"
                with self.assertRaisesRegex(ValueError, "decision differs"):
                    callback(route, path, folder / "rate")
                self.assertEqual(decision_path.read_bytes(), rejected_bytes)
                decision_path.write_text("corrupt retained decision\n")
                with self.assertRaisesRegex(ValueError, "unreadable or corrupt"):
                    callback(route, path, folder / "rate")
                self.assertEqual(decision_path.read_text(), "corrupt retained decision\n")
                run.assert_not_called()

    def test_multistep_adapter_cannot_overwrite_a_contradictory_segment_classification(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            route, single_path = _ordinary_route_and_path(folder)
            candidate = ReactionCandidate.create(
                ("C1=CCC1",), ("C=CC=C",), total_charge=0, multiplicity=1,
                reactant_labels=("cyclobutene",), product_labels=("butadiene",),
                reactant_charges=(0,), product_charges=(0,),
                reactant_multiplicities=(1,), product_multiplicities=(1,),
            )
            segment = {"status": "verified_transition_state_path", "rate_input_status": "ready",
                       "native_candidate": candidate.as_dict(),
                       "verification": {"status": "verified_transition_state_path",
                                        "path_classification": "barrierless_capture_candidate"},
                       "arkane_inputs": single_path["route_verification"]["arkane_inputs"]}
            path = {"route_verification": {"status": "verified_multistep_path",
                                           "multistep_path": {"segments": [segment]}}}
            callback = make_arkane_rate_builder(_condition(), method_profile=_profile(), rmg_env="rmg_env")
            with patch("storca.verification_adapters.run_verified_arkane_route",
                       return_value={"status": "completed"}) as run:
                with self.assertRaisesRegex(ValueError, "every elementary segment"):
                    callback(route, path, folder / "rate")
                run.assert_not_called()

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

    def test_default_orca_adapter_uses_native_candidate_normalizer(self):
        normalized = RouteSpec(
            route_id="graph-route", source="graph_xtb", parent_smiles="CC",
            reactant_smiles=("CC",), product_smiles=("[CH3]", "[CH3]"), multiplicity=1,
        )
        path_result = {"route_verification": {"status": "mapping_unresolved"}}
        route = {
            "native_candidate_id": "graph-route", "route_id": "graph-route",
            "source": "graph_xtb", "reaction_equation": "stability <=> p1 + p2",
        }
        with tempfile.TemporaryDirectory() as temporary, \
                patch("storca.verification_adapters.route_spec_from_candidate_route", return_value=normalized) as normalize, \
                patch("storca.verification_adapters.route_spec_from_rmg_route") as rmg_normalize, \
                patch("storca.verification_adapters.run_generic_reaction_path", return_value=path_result), \
                patch("storca.verification_adapters._validated_trajectory_visual", return_value=None):
            callback = make_orca_route_verifier("CC", ncores=1, method_profile=_profile())
            self.assertIs(callback(route, Path(temporary)), path_result)
        normalize.assert_called_once_with("CC", route)
        rmg_normalize.assert_not_called()

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
        route["pressure_dependence_model"] = {"sentinel": "explicit-route-model"}
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            hydrogen = folder / "hydrogen.out"
            product = folder / "hydrogen2.out"
            ts = folder / "ts.out"
            path_result = {
                "route_verification": {
                    "route_id": "recombination",
                    "status": "verified_transition_state_path",
                    "path_classification": "verified_barriered_path",
                    "selected_multiplicity": 1,
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
            self.assertEqual(spec.transition_state_multiplicity, 1)
            self.assertEqual(spec.reaction_equation, route["reaction_equation"])
            self.assertEqual(
                spec.explicit_pressure_dependence_model,
                {"sentinel": "explicit-route-model"},
            )
            self.assertEqual(spec.collision_limit_m3_mol_s, 4.2)
            self.assertEqual(spec.collision_limit_source, "RMG retained transport collision bound")
            self.assertEqual(dict(spec.bath_gas), {"nitrogen": 0.74 / 0.99, "oxygen": 0.20 / 0.99,
                                                   "water": 0.05 / 0.99})
            self.assertIn(310.0, spec.temperatures_K)
            self.assertIn(0.8, spec.pressures_bar)
            self.assertEqual(captured["rmg_env"], "rmg4_env")

            replacement, paths = _normalize_replacement(result, route)
            for field, value in nested.items():
                self.assertEqual(replacement[field], value)
            self.assertNotIn("method_applicability", nested)
            self.assertEqual(replacement["method_applicability"], result["method_applicability"])
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
                "path_classification": "verified_barriered_path",
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

    def test_arkane_adapter_builds_each_verified_multistep_segment_without_effective_barrier(self):
        shared = dict(
            total_charge=0, multiplicity=1,
            discovery_assertions=({"source": "validated_neb_intermediate"},),
            verification_status="segment_path_unresolved",
            provenance={"atom_mapping_representation": "smiles_atom_map_numbers"},
            reactant_charges=(0,), product_charges=(0,),
            reactant_multiplicities=(1,), product_multiplicities=(1,),
        )
        candidates = [
            ReactionCandidate.create(
                ("[CH3:1][CH2:2][CH2:3][OH:4]",),
                ("[CH3:1][CH:2]([CH3:3])[OH:4]",),
                reactant_labels=("propanol",), product_labels=("isopropanol",), **shared,
            ),
            ReactionCandidate.create(
                ("[CH3:1][CH:2]([CH3:3])[OH:4]",),
                ("[CH3:1][CH2:2][O:4][CH3:3]",),
                reactant_labels=("isopropanol",), product_labels=("methoxyethane",), **shared,
            ),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            endpoint_outputs = []
            ts_outputs = []
            for index in range(3):
                path = folder / f"endpoint-{index}.out"
                path.write_text("retained endpoint frequency\n")
                endpoint_outputs.append(path)
            for index in range(2):
                path = folder / f"ts-{index}.out"
                path.write_text("retained TS frequency\n")
                ts_outputs.append(path)
            segments = []
            for index, candidate in enumerate(candidates):
                reactant_label = candidate.reactant_labels[0]
                product_label = candidate.product_labels[0]
                segments.append({
                    "label": f"segment-{index + 1}",
                    "status": "verified_transition_state_path",
                    "rate_input_status": "ready",
                    "native_candidate": candidate.as_dict(),
                    "verification": {
                        "status": "verified_transition_state_path",
                        "path_classification": "verified_barriered_path",
                        "transition_state": {"frequency_output": str(ts_outputs[index])},
                    },
                    "arkane_inputs": {
                        "reactants": [{
                            "label": reactant_label, "original_label": reactant_label,
                            "smiles": candidate.reactants[0], "multiplicity": 1,
                            "status": "validated_minimum",
                            "frequency_output": str(endpoint_outputs[index]),
                        }],
                        "products": [{
                            "label": product_label, "original_label": product_label,
                            "smiles": candidate.products[0], "multiplicity": 1,
                            "status": "validated_minimum",
                            "frequency_output": str(endpoint_outputs[index + 1]),
                        }],
                        "transition_state": {"frequency_output": str(ts_outputs[index])},
                    },
                })
            path_result = {"route_verification": {
                "status": "verified_multistep_path",
                "multistep_path": {"segments": segments},
            }}
            calls = []

            def fake_run(spec, workdir, *, rmg_env):
                calls.append(spec)
                library = folder / f"library-{len(calls)}"
                library.mkdir()
                return {
                    "status": "completed",
                    "generated_kinetics_library": {
                        "verification_status": "verified_arkane_pdep",
                        "route_id": spec.label,
                        "library": str(library),
                    },
                }

            with patch("storca.verification_adapters.run_verified_arkane_route", side_effect=fake_run):
                callback = make_arkane_rate_builder(
                    _condition(), method_profile=_profile(), rmg_env="rmg4_env",
                )
                result = callback(
                    {"route_id": "parent-multistep"}, path_result, folder / "rates",
                )

        self.assertEqual(result["status"], "completed_multistep")
        self.assertEqual(len(result["segments"]), 2)
        self.assertEqual(len(calls), 2)
        self.assertFalse(result["effective_single_barrier_used"])
        self.assertTrue(result["network_insertion_required"])
        self.assertFalse(result["lifetime_claim_allowed"])


if __name__ == "__main__":
    unittest.main()
