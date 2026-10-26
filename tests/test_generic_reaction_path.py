from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from storca.reaction_path import (
    _imaginary_mode_displacement_seed,
    _initialize_surface_context,
    _initialize_surface_budget,
    _neb_angle_warning_image_indices,
    _neb_intermediate_image_indices,
    _neb_monotonic_endpoint_diagnostic,
    _neb_persistent_endpoint_highest_image_diagnostic,
    _neb_submerged_climbing_image_diagnostic,
    _orca_neb_no_barrier_outcome,
    _optimize_endpoint,
    _run_generic_orientation,
    _run_bond_change_guided_ts_recovery,
    _run_or_resume,
    _run_endpoint_job,
    _run_segmented_intermediate_nebs,
    _run_stationary_point_recovery,
    _validate_orca_neb_restart_band,
    _stationary_fragment_records,
    _surface_budget_snapshot,
    _surface_environment_applicability,
    _validate_neb_intermediates,
    _write_surface_attempt_evidence,
    classify_path_energy_profile,
    run_generic_reaction_path,
)
from storca.ir_modes import NormalModeSet
from storca.reaction_candidates import ReactionCandidate
from storca.route_geometry import (
    build_endpoint_complex_seeds,
    read_xyz,
    resolve_route_atom_mapping,
    validate_declared_connectivity,
    write_xyz,
)
from storca.route_verify import RouteSpec, route_spec_from_rmg_route
from storca.surface_resolution import SurfaceResolutionPolicy
from storca.verification_policy import PREVIOUS_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID


class GenericReactionPathTests(unittest.TestCase):
    @staticmethod
    def _bound_segment_identity() -> dict:
        shared = dict(
            total_charge=0,
            multiplicity=1,
            discovery_assertions=({"source": "validated_neb_intermediate"},),
            verification_status="segment_path_unresolved",
            provenance={"atom_mapping_representation": "smiles_atom_map_numbers"},
            reactant_charges=(0,),
            product_charges=(0,),
            reactant_multiplicities=(1,),
            product_multiplicities=(1,),
        )
        first = ReactionCandidate.create(
            ("[CH3:1][CH2:2][CH2:3][OH:4]",),
            ("[CH3:1][CH:2]([CH3:3])[OH:4]",),
            reactant_labels=("propan-1-ol",),
            product_labels=("propan-2-ol",),
            **shared,
        )
        second = ReactionCandidate.create(
            ("[CH3:1][CH:2]([CH3:3])[OH:4]",),
            ("[CH3:1][CH2:2][O:4][CH3:3]",),
            reactant_labels=("propan-2-ol",),
            product_labels=("methoxyethane",),
            **shared,
        )
        return {"status": "resolved", "native_candidates": [first.as_dict(), second.as_dict()]}

    @staticmethod
    def _write_neb_restart_fixture(folder: Path) -> Path:
        inp = folder / "neb-ts.inp"
        inp.write_text(
            "! B3LYP def2-SVP NoAutoStart NEB-TS\n"
            '%neb\n  Product "product.xyz"\n  NImages 2\nend\n'
            "* xyzfile 0 1 reactant.xyz\n"
        )
        (folder / "reactant.xyz").write_text("1\nreactant\nH 0 0 0\n")
        (folder / "product.xyz").write_text("1\nproduct\nH 0 0 1\n")
        (folder / "neb-ts.out").write_text(
            "Starting iterations:\n   LBFGS     0      1    0.100000    0.050000\n"
        )
        frames = []
        for index in range(4):
            frames.append(f"1\nimage {index}\nH 0 0 {index / 3:.6f}\n>\n")
        (folder / "neb-ts_MEP.allxyz").write_text("".join(frames))
        return inp

    def test_neb_restart_band_requires_complete_consistent_image_set(self):
        with tempfile.TemporaryDirectory() as temp:
            band = Path(temp) / "band.allxyz"
            band.write_text("1\none\nH 0 0 0\n>\n1\ntwo\nH 0 0 1\n")
            valid = _validate_orca_neb_restart_band(band, expected_images=2)
            mismatch = _validate_orca_neb_restart_band(band, expected_images=3)
        self.assertTrue(valid["valid"])
        self.assertEqual(mismatch["status"], "image_count_mismatch")

    def test_surface_context_policy_upgrade_is_append_only(self):
        route = RouteSpec(
            route_id="policy-amendment", source="test", parent_smiles="CC",
            reactant_smiles=("CC",), product_smiles=("CC",),
            protocol="generic_endpoint_path",
        )
        mapping = {"status": "mapped", "valid": True, "pairs": [[0, 0], [1, 1]]}
        environment = {"status": "applicable", "environment_model": "homogeneous_gas"}
        previous = SurfaceResolutionPolicy(
            policy_id=PREVIOUS_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
            intermediate_identity_policy_id="rdkit-connectivity-first-charge-spin-map-v2",
        )
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            original = _initialize_surface_context(
                root, route=route, mapping=mapping,
                environment_applicability=environment, policy=previous,
                ncores=2, method_keywords=None,
            )
            original_bytes = Path(original["path"]).read_bytes()
            amended = _initialize_surface_context(
                root, route=route, mapping=mapping,
                environment_applicability=environment,
                policy=SurfaceResolutionPolicy(), ncores=4,
                method_keywords=None,
            )
            reused = _initialize_surface_context(
                root, route=route, mapping=mapping,
                environment_applicability=environment,
                policy=SurfaceResolutionPolicy(), ncores=4,
                method_keywords=None,
            )

            self.assertEqual(Path(original["path"]).read_bytes(), original_bytes)
            self.assertEqual(Path(amended["path"]).name, "000002.json")
            self.assertEqual(Path(reused["path"]).resolve(), Path(amended["path"]).resolve())
            self.assertEqual(
                amended["predecessor_scientific_context_sha256"],
                original["scientific_context_sha256"],
            )
            self.assertEqual(
                amended["scientific_change"]["kind"],
                "surface_resolution_policy_amendment",
            )
            self.assertEqual(
                amended["scientific_change"]["execution_resource_change"],
                {"prior_ncores": 2, "new_ncores": 4},
            )

    def test_surface_context_rejects_same_policy_definition_mutation(self):
        route = RouteSpec(
            route_id="policy-mutation", source="test", parent_smiles="CC",
            reactant_smiles=("CC",), product_smiles=("CC",),
            protocol="generic_endpoint_path",
        )
        mapping = {"status": "mapped", "valid": True}
        environment = {"status": "applicable", "environment_model": "homogeneous_gas"}
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _initialize_surface_context(
                root, route=route, mapping=mapping,
                environment_applicability=environment,
                policy=SurfaceResolutionPolicy(), ncores=2,
                method_keywords=None,
            )
            with self.assertRaisesRegex(RuntimeError, "without a new policy ID"):
                _initialize_surface_context(
                    root, route=route, mapping=mapping,
                    environment_applicability=environment,
                    policy=SurfaceResolutionPolicy(maximum_scan_points=17),
                    ncores=2, method_keywords=None,
                )

    def test_surface_attempt_policy_upgrade_appends_instead_of_mutating(self):
        route = RouteSpec(
            route_id="attempt-policy-amendment", source="test", parent_smiles="CC",
            reactant_smiles=("CC",), product_smiles=("CC",),
            protocol="generic_endpoint_path",
        )
        mapping = {"status": "mapped", "valid": True}
        environment = {"status": "applicable", "environment_model": "homogeneous_gas"}
        previous = SurfaceResolutionPolicy(
            policy_id=PREVIOUS_GENERIC_ORCA_PATH_VERIFIER_POLICY_ID,
            intermediate_identity_policy_id="rdkit-connectivity-first-charge-spin-map-v2",
        )
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _initialize_surface_context(
                root, route=route, mapping=mapping,
                environment_applicability=environment, policy=previous,
                ncores=2, method_keywords=None,
            )
            first = _write_surface_attempt_evidence(
                root, strategy="standard_neb_ts", status="surface_unresolved",
                result={"kind": "fixture", "status": "surface_unresolved"},
            )
            _initialize_surface_context(
                root, route=route, mapping=mapping,
                environment_applicability=environment,
                policy=SurfaceResolutionPolicy(), ncores=2, method_keywords=None,
            )
            second = _write_surface_attempt_evidence(
                root, strategy="standard_neb_ts", status="surface_unresolved",
                result={"kind": "fixture", "status": "surface_unresolved"},
            )

        self.assertNotEqual(first["attempt_id"], second["attempt_id"])
        self.assertEqual(Path(first["path"]).name, "000001.json")
        self.assertEqual(Path(second["path"]).name, "000002.json")
        self.assertNotEqual(
            first["scientific_context"]["scientific_context_sha256"],
            second["scientific_context"]["scientific_context_sha256"],
        )

    def test_unsupported_environment_escalates_before_orca(self):
        route = RouteSpec(
            route_id="surface-route",
            source="test",
            parent_smiles="CC",
            reactant_smiles=("CC",),
            product_smiles=("CC",),
            charge=0,
            multiplicity=1,
            protocol="generic_endpoint_path",
            environment={"phase": "periodic_surface"},
        )
        with tempfile.TemporaryDirectory() as temp, patch(
            "storca.reaction_path._run_or_resume",
        ) as execute:
            result = run_generic_reaction_path(route, Path(temp))

        self.assertEqual(
            result["route_verification"]["status"],
            "environment_model_escalation_required",
        )
        self.assertFalse(result["environment_applicability"]["rate_claim_allowed"])
        execute.assert_not_called()

    def test_implicit_solution_requires_exact_solvent_method_contract(self):
        route = RouteSpec(
            route_id="solution-route", source="test", parent_smiles="CC",
            reactant_smiles=("CC",), product_smiles=("CC",), charge=0,
            multiplicity=1, protocol="generic_endpoint_path",
        )
        condition = {
            "phase_model": "homogeneous implicit solution",
            "solvent": "water", "solvent_model": "SMD",
        }
        accepted = _surface_environment_applicability(
            route, condition, SurfaceResolutionPolicy(), ["r2SCAN-3c", "SMD(water)"],
        )
        rejected = _surface_environment_applicability(
            route, condition, SurfaceResolutionPolicy(), ["r2SCAN-3c"],
        )
        mismatched = _surface_environment_applicability(
            route, condition, SurfaceResolutionPolicy(), ["r2SCAN-3c", "CPCM(water)"],
        )
        self.assertEqual(accepted["status"], "applicable")
        self.assertEqual(rejected["status"], "environment_model_escalation_required")
        self.assertIn("orca_implicit_solvent_keyword_missing_or_mismatched", rejected["blockers"])
        self.assertEqual(mismatched["status"], "environment_model_escalation_required")

    def test_gas_route_rejects_accidental_solvent_keyword(self):
        route = RouteSpec(
            route_id="gas-route", source="test", parent_smiles="CC",
            reactant_smiles=("CC",), product_smiles=("CC",), charge=0,
            multiplicity=1, protocol="generic_endpoint_path",
        )
        result = _surface_environment_applicability(
            route, {"phase_model": "homogeneous gas"}, SurfaceResolutionPolicy(),
            ["r2SCAN-3c", "SMD(water)"],
        )
        self.assertEqual(result["status"], "environment_model_escalation_required")
        self.assertIn(
            "gas_or_nonsolution_route_contains_implicit_solvent_keyword", result["blockers"],
        )

    @patch("storca.reaction_path.build_endpoint_complex_seeds")
    @patch("storca.reaction_path.resolve_route_atom_mapping")
    def test_spin_surface_failure_uses_top_level_escalation_state(self, mapping, build_seeds):
        mapping.return_value = {
            "status": "no_common_adiabatic_spin_surface", "valid": False,
            "reactant_total_multiplicities": [1],
            "product_total_multiplicities": [2],
            "surface": {
                "status": "no_common_adiabatic_spin_surface", "valid": False,
                "reactant_total_multiplicities": [1],
                "product_total_multiplicities": [2],
            },
        }
        route = RouteSpec(
            route_id="spin-route", source="test", parent_smiles="CC",
            reactant_smiles=("CC",), product_smiles=("CC",), charge=0,
            multiplicity=1, protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(route, Path(temp))
        verification = result["route_verification"]
        self.assertEqual(verification["status"], "electronic_state_escalation_required")
        self.assertEqual(
            verification["electronic_state_evidence"]["status"],
            "no_common_adiabatic_spin_surface",
        )
        self.assertIsNone(verification["trajectory"])
        build_seeds.assert_not_called()

    @patch("storca.reaction_path.run_orca")
    def test_fresh_partial_neb_refuses_duplicate_restart(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            inp = self._write_neb_restart_fixture(Path(temp))
            with self.assertRaisesRegex(RuntimeError, "fresh heartbeat"):
                _run_or_resume(inp, timeout_seconds=10.0)
        run_orca.assert_not_called()

    @patch("storca.reaction_path.run_orca")
    def test_quiet_partial_neb_restarts_from_retained_allxyz_band(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = self._write_neb_restart_fixture(folder)
            for artifact in folder.iterdir():
                os.utime(artifact, (1, 1))
            # High-level route preparation may regenerate the byte-identical
            # input immediately before resume; that is not ORCA activity.
            inp.write_text(inp.read_text())

            def completed(restart_input, **kwargs):
                output = Path(restart_input).with_suffix(".out")
                output.write_text("ORCA TERMINATED NORMALLY\n")
                return {"out": output}

            run_orca.side_effect = completed
            result = _run_or_resume(inp, timeout_seconds=10.0)
            restart_input = folder / "neb-ts-restart-001.inp"
            restart_record = folder / "neb-ts-restart-001.restart.json"

            self.assertEqual(Path(result["input"]), restart_input)
            self.assertIn('Restart_ALLXYZFile "neb-ts_MEP.allxyz"', restart_input.read_text())
            self.assertTrue(restart_record.is_file())
            self.assertTrue(result["neb_checkpoint_restart"]["band_validation"]["valid"])
            run_orca.assert_called_once()

    @patch("storca.reaction_path.run_orca")
    def test_neb_restart_does_not_request_unsupported_endpoint_preoptimization(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = self._write_neb_restart_fixture(folder)
            inp.write_text(inp.read_text().replace("  NImages 2", "  PreOpt true\n  NImages 2"))
            original = inp.read_bytes()
            for artifact in folder.iterdir():
                os.utime(artifact, (1, 1))
            def completed(restart_input, **kwargs):
                self.assertIn("PreOpt false", restart_input.read_text())
                self.assertNotIn("PreOpt true", restart_input.read_text())
                out = restart_input.with_suffix(".out")
                out.write_text("ORCA TERMINATED NORMALLY\n")
                return {"out": out}
            run_orca.side_effect = completed
            _run_or_resume(inp, timeout_seconds=600)
            self.assertEqual(inp.read_bytes(), original)

    @patch("storca.reaction_path.run_orca")
    def test_validated_endpoint_return_warning_gets_one_checkpoint_continuation(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = self._write_neb_restart_fixture(folder)
            inp.with_suffix(".out").write_text(
                inp.with_suffix(".out").read_text()
                + "Possible intermediate minimum found at image(s): 2\n"
            )
            for artifact in folder.iterdir():
                os.utime(artifact, (1, 1))

            def completed(restart_input, **kwargs):
                stop = kwargs["stream_output"](
                    "Possible intermediate minimum found at image(s): 2\n"
                )
                self.assertFalse(stop)
                output = Path(restart_input).with_suffix(".out")
                output.write_text(
                    "Possible intermediate minimum found at image(s): 2\n"
                    "ORCA TERMINATED NORMALLY\n"
                )
                return {"out": output}

            run_orca.side_effect = completed
            result = _run_or_resume(
                inp,
                timeout_seconds=10.0,
                ignore_intermediate_warning_indices=[2],
            )

        self.assertEqual(result["status"], "completed")
        self.assertEqual(
            result["validated_endpoint_return_warning_policy"]["ignored_image_indices"],
            [2],
        )
        self.assertFalse(
            result["validated_endpoint_return_warning_policy"]["scientific_input_changed"]
        )

    @patch("storca.reaction_path.run_orca")
    def test_job_contract_rejects_changed_input_or_geometry_without_reexecution(self, run_orca):
        for mutation in ("input", "referenced_geometry"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                folder = Path(temp)
                inp = folder / "opt.inp"
                geometry = folder / "seed.xyz"
                geometry.write_text("1\nseed\nH 0 0 0\n")
                inp.write_text("! B3LYP Opt\n* xyzfile 0 2 seed.xyz\n")

                def completed(input_path, **kwargs):
                    output = Path(input_path).with_suffix(".out")
                    output.write_text("ORCA TERMINATED NORMALLY\n")
                    return {"out": output}

                run_orca.reset_mock()
                run_orca.side_effect = completed
                _run_or_resume(inp, timeout_seconds=10.)
                retained_contract = inp.with_suffix(".contract.json").read_bytes()
                retained_output = inp.with_suffix(".out").read_bytes()
                if mutation == "input":
                    inp.write_text(inp.read_text().replace("B3LYP", "PBE0"))
                else:
                    geometry.write_text("1\nchanged\nH 0.1 0 0\n")
                with self.assertRaisesRegex(RuntimeError, "refusing to replace its job contract"):
                    _run_or_resume(inp, timeout_seconds=10.)
                self.assertEqual(run_orca.call_count, 1)
                self.assertEqual(inp.with_suffix(".contract.json").read_bytes(), retained_contract)
                self.assertEqual(inp.with_suffix(".out").read_bytes(), retained_output)

    @patch("storca.reaction_path.run_orca")
    def test_orca_startup_failure_retries_append_only_without_scientific_change(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "opt.inp"
            inp.write_text("! B3LYP def2-SVP Opt\n* xyz 0 1\nH 0 0 0\n*\n")
            original_output = inp.with_suffix(".out")
            startup_text = (
                "bind() failed for port 0: Operation not permitted\n"
                "ORCA finished by error termination in Startup\n"
            )

            def startup_failure(input_path, **kwargs):
                Path(input_path).with_suffix(".out").write_text(startup_text)
                raise RuntimeError("ORCA exited without normal termination")

            def completed(input_path, **kwargs):
                output = Path(input_path).with_suffix(".out")
                output.write_text("ORCA TERMINATED NORMALLY\n")
                return {"out": output}

            calls = iter((startup_failure, completed))
            run_orca.side_effect = lambda input_path, **kwargs: next(calls)(input_path, **kwargs)
            result = _run_or_resume(inp, timeout_seconds=10.0)

            retry = folder / "opt-execution-retry-001.inp"
            retry_record = retry.with_suffix(".restart.json")
            retry_policy = retry.with_suffix(".execution-policy.json")
            self.assertEqual(original_output.read_text(), startup_text)
            self.assertEqual(Path(result["input"]), retry)
            self.assertEqual(retry.read_text(), inp.read_text())
            self.assertTrue(retry_record.is_file())
            self.assertTrue(retry_policy.is_file())
            self.assertFalse(result["execution_startup_retry"]["scientific_input_changed"])
            self.assertEqual(
                json.loads(retry_policy.read_text())["strategy"],
                "safe_openmpi_transport_retry",
            )
            self.assertEqual(
                run_orca.call_args_list[1].kwargs["extra_env"]["OMPI_MCA_btl"],
                "self,tcp",
            )

    @patch("storca.reaction_path.run_orca")
    def test_repeated_tool_connection_failure_falls_back_to_serial_execution(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "opt.inp"
            inp.write_text(
                "! B3LYP def2-SVP Opt\n%pal\n  nprocs 4\nend\n"
                "* xyz 0 1\nH 0 0 0\n*\n"
            )
            startup_text = (
                "ptl_tool: problems getting address for index 0\n"
                "You requested support for tool connections, but no available "
                "interfaces were found\n"
                "ORCA finished by error termination in Startup\n"
            )

            def fail_startup(input_path, **kwargs):
                Path(input_path).with_suffix(".out").write_text(startup_text)
                raise RuntimeError("ORCA exited without normal termination")

            def completed(input_path, **kwargs):
                output = Path(input_path).with_suffix(".out")
                output.write_text("ORCA TERMINATED NORMALLY\n")
                return {"out": output}

            functions = iter((fail_startup, fail_startup, completed))
            run_orca.side_effect = (
                lambda input_path, **kwargs: next(functions)(input_path, **kwargs)
            )
            result = _run_or_resume(inp, timeout_seconds=10.0)

            safe = folder / "opt-execution-retry-001.inp"
            serial = folder / "opt-execution-retry-002.inp"
            self.assertEqual(Path(result["input"]), serial)
            self.assertEqual(safe.read_text(), inp.read_text())
            self.assertIn("nprocs 1", serial.read_text())
            policy = json.loads(serial.with_suffix(".execution-policy.json").read_text())
            self.assertEqual(policy["strategy"], "serial_execution_fallback")
            self.assertEqual(policy["requested_ncores"], 4)
            self.assertEqual(policy["effective_ncores"], 1)
            self.assertFalse(policy["scientific_model_changed"])
            self.assertTrue(policy["execution_resources_changed"])

    @patch("storca.reaction_path.run_orca")
    def test_stopped_serial_fallback_restores_retained_parallel_contract(self, run_orca):
        route = RouteSpec(
            route_id="resource-recovery", source="test", parent_smiles="CC",
            reactant_smiles=("CC",), product_smiles=("CC",),
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            _initialize_surface_context(
                folder, route=route,
                mapping={"status": "mapped", "valid": True},
                environment_applicability={
                    "status": "applicable", "environment_model": "homogeneous_gas",
                },
                policy=SurfaceResolutionPolicy(), ncores=4,
                method_keywords=None,
            )
            original = folder / "opt.inp"
            original.write_text(
                "! B3LYP def2-SVP Opt\n%pal\n  nprocs 4\nend\n"
                "* xyz 0 1\nH 0 0 0\n*\n"
            )
            safe = folder / "opt-execution-retry-001.inp"
            safe.write_text(original.read_text())
            serial = folder / "opt-execution-retry-002.inp"
            serial.write_text(original.read_text().replace("nprocs 4", "nprocs 1"))
            serial.with_suffix(".out").write_text("interrupted serial fallback\n")
            serial.with_suffix(".contract.json").write_text(json.dumps({
                "schema_version": 1, "kind": "orca_job_input_contract",
                "sha256": hashlib.sha256(serial.read_bytes()).hexdigest(),
                "input": str(serial), "referenced_files": [], "status": "started",
            }))
            serial.with_suffix(".process.json").write_text(json.dumps({
                "kind": "storca_external_process_lease",
                "status": "runner_detached_process_may_be_active",
                "pid": 99999991, "process_group_id": 99999991,
            }))
            serial.with_suffix(".execution-policy.json").write_text(json.dumps({
                "extra_environment": {"OMPI_MCA_btl": "self,tcp"},
            }))

            def completed(input_path, **kwargs):
                output = Path(input_path).with_suffix(".out")
                output.write_text("ORCA TERMINATED NORMALLY\n")
                return {"out": output}

            run_orca.side_effect = completed
            result = _run_or_resume(original, timeout_seconds=10.0)
            recovered = folder / "opt-execution-retry-003.inp"
            policy = json.loads(
                recovered.with_suffix(".execution-policy.json").read_text()
            )
            recovered_text = recovered.read_text()

        self.assertEqual(Path(result["input"]), recovered)
        self.assertIn("nprocs 4", recovered_text)
        self.assertEqual(policy["prior_ncores"], 1)
        self.assertEqual(policy["effective_ncores"], 4)
        self.assertFalse(policy["scientific_model_changed"])
        self.assertEqual(
            result["execution_resource_contract_recovery"]["kind"],
            "orca_execution_resource_contract_recovery",
        )

    @patch("storca.reaction_path.run_orca")
    def test_confirmed_stopped_interruption_gets_new_execution_basename(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "opt.inp"
            inp.write_text("! B3LYP def2-SVP Opt\n* xyz 0 1\nH 0 0 0\n*\n")
            partial = "interrupted electronic-structure output\n"
            inp.with_suffix(".out").write_text(partial)
            inp.with_suffix(".contract.json").write_text(json.dumps({
                "schema_version": 1,
                "kind": "orca_job_input_contract",
                "sha256": hashlib.sha256(inp.read_bytes()).hexdigest(),
                "input": str(inp),
                "referenced_files": [],
                "status": "started",
            }))
            inp.with_suffix(".process.json").write_text(json.dumps({
                "kind": "storca_external_process_lease",
                "status": "runner_detached_process_may_be_active",
                "pid": 99999991,
                "process_group_id": 99999991,
            }))

            def completed(input_path, **kwargs):
                output = Path(input_path).with_suffix(".out")
                output.write_text("ORCA TERMINATED NORMALLY\n")
                return {"out": output}

            run_orca.side_effect = completed
            result = _run_or_resume(inp, timeout_seconds=10.0)

            self.assertEqual(inp.with_suffix(".out").read_text(), partial)
            self.assertTrue(str(result["input"]).endswith("opt-execution-retry-001.inp"))
            self.assertEqual(
                result["execution_startup_retry"]["kind"],
                "orca_interrupted_execution_retry",
            )

    @patch("storca.reaction_path.run_orca")
    def test_endpoint_legacy_startup_retry_preserves_output_and_reuses_completion(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "opt.inp"
            inp.write_text("! B3LYP def2-SVP Opt\n%pal\n nprocs 4\nend\n* xyz 0 1\nH 0 0 0\n*\n")
            failed = "ORCA finished by error termination in Startup\n"
            inp.with_suffix(".out").write_text(failed)

            def completed(input_path, **kwargs):
                output = Path(input_path).with_suffix(".out")
                output.write_text("ORCA TERMINATED NORMALLY\n")
                Path(input_path).with_suffix(".xyz").write_text("1\nrecovered\nH 0 0 0\n")
                return {"out": output}

            run_orca.side_effect = completed
            result = _run_endpoint_job(inp, timeout_seconds=10.0)
            resumed = _run_endpoint_job(inp, timeout_seconds=10.0)
            self.assertEqual(inp.with_suffix(".out").read_text(), failed)
            self.assertEqual(result["out"].name, "opt-execution-retry-001.out")
            self.assertEqual(result["xyz"].name, "opt-execution-retry-001.xyz")
            self.assertEqual(resumed["xyz"], result["xyz"])
            self.assertTrue(resumed["resumed"])
            self.assertEqual(run_orca.call_count, 1)
            evidence = json.loads((folder / "opt-execution-retry-001.restart.json").read_text())
            self.assertEqual(evidence["source_output_sha256"], hashlib.sha256(failed.encode()).hexdigest())
            inp.write_text(inp.read_text().replace("B3LYP", "PBE0"))
            with self.assertRaisesRegex(RuntimeError, "no longer matches"):
                _run_endpoint_job(inp, timeout_seconds=10.0)
            self.assertEqual(run_orca.call_count, 1)

    @patch("storca.reaction_path.run_orca")
    @patch("storca.workflow.frequency_stability_check", return_value={"IsMinimum": True})
    def test_endpoint_workflow_frequencies_use_recovered_geometry(self, frequency_check, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            seed = write_xyz(folder / "input.xyz", ["H"], np.asarray([[0., 0., 0.]]), "seed")
            failed = "ORCA finished by error termination in Startup\n"
            (folder / "opt.out").write_text(failed)

            def completed(input_path, **kwargs):
                input_path = Path(input_path)
                output = input_path.with_suffix(".out")
                output.write_text("THE OPTIMIZATION HAS CONVERGED\nORCA TERMINATED NORMALLY\n")
                write_xyz(input_path.with_suffix(".xyz"), ["H"], np.asarray([[0.123, 0., 0.]]), "optimized")
                return {"out": output}

            run_orca.side_effect = completed
            for _ in range(2):
                result = _optimize_endpoint(
                    seed, folder, charge=0, multiplicity=2, ncores=4,
                    method_keywords=None, timeout_seconds=10.,
                )
                self.assertTrue(result["frequency_check"]["IsMinimum"])
            self.assertEqual(run_orca.call_count, 2)
            self.assertEqual((folder / "opt.out").read_text(), failed)
            self.assertEqual(Path(result["optimized_xyz"]).name, "opt-execution-retry-001.xyz")
            self.assertIn("0.1230000000", (folder / "freq.inp").read_text())
            retained_input = (folder / "opt.inp").read_bytes()
            with self.assertRaisesRegex(RuntimeError, "refusing to overwrite"):
                _optimize_endpoint(
                    seed, folder, charge=0, multiplicity=2, ncores=4,
                    method_keywords=["PBE0", "def2-SVP"], timeout_seconds=10.,
                )
            self.assertEqual((folder / "opt.inp").read_bytes(), retained_input)
            self.assertEqual(run_orca.call_count, 2)

    @patch("storca.reaction_path.run_optimization_and_frequency")
    def test_changed_endpoint_seed_is_rejected_before_copying(self, workflow):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            workdir = folder / "endpoint"
            workdir.mkdir()
            retained = workdir / "input.xyz"
            retained.write_text("1\noriginal\nH 0 0 0\n")
            changed = folder / "new.xyz"
            changed.write_text("1\nchanged\nH 0.1 0 0\n")
            before = retained.read_bytes()
            with self.assertRaisesRegex(RuntimeError, "Retained endpoint seed differs"):
                _optimize_endpoint(
                    changed, workdir, charge=0, multiplicity=2, ncores=4,
                    method_keywords=None, timeout_seconds=10.,
                )
            self.assertEqual(retained.read_bytes(), before)
            workflow.assert_not_called()

    @patch("storca.reaction_path.run_orca")
    def test_endpoint_unknown_incomplete_output_is_not_overwritten(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            inp = Path(temp) / "opt.inp"
            inp.write_text("! B3LYP Opt\n")
            inp.with_suffix(".out").write_text("unresolved execution\n")
            with self.assertRaisesRegex(RuntimeError, "refusing to overwrite"):
                _run_endpoint_job(inp, timeout_seconds=10.0)
            self.assertEqual(inp.with_suffix(".out").read_text(), "unresolved execution\n")
            run_orca.assert_not_called()

    @patch("storca.reaction_path.run_optimization_and_frequency")
    @patch("storca.reaction_path.frequency_stability_check", return_value={"IsMinimum": True})
    @patch("storca.reaction_path._orca_optimization_converged", return_value=True)
    @patch("storca.reaction_path._normal_orca_output", return_value=True)
    def test_legacy_endpoint_contract_is_upgraded_without_rerunning_orca(
        self, normal, converged, frequency_check, run_calculation,
    ):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            seed = write_xyz(
                folder / "input.xyz", ["H"], np.asarray([[0.0, 0.0, 0.0]]), "seed",
            )
            write_xyz(folder / "opt.xyz", ["H"], np.asarray([[0.0, 0.0, 0.0]]), "optimized")
            (folder / "opt.out").write_text("retained optimization")
            (folder / "freq.out").write_text("retained frequencies")
            (folder / "endpoint-contract.json").write_text(json.dumps({
                "schema_version": 1,
                "kind": "orca_endpoint_contract",
                "input_xyz_sha256": hashlib.sha256(seed.read_bytes()).hexdigest(),
                "charge": 0,
                "multiplicity": 2,
                "ncores": 1,
                "method_keywords": ["B3LYP", "def2-SVP", "NoAutoStart"],
            }))
            result = _optimize_endpoint(
                seed, folder, charge=0, multiplicity=2, ncores=1,
                method_keywords=None, timeout_seconds=10.0,
            )
            upgraded = json.loads((folder / "endpoint-contract.json").read_text())

        self.assertTrue(result["resumed"])
        self.assertEqual(upgraded["schema_version"], 2)
        run_calculation.assert_not_called()

    @patch("storca.reaction_path._run_or_resume")
    def test_unconverged_stationary_point_continues_append_only(self, run_or_resume):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            source = write_xyz(
                folder / "base-opt.xyz", ["H"], np.asarray([[0.0, 0.0, 0.0]]), "base",
            )
            base_out = folder / "base-opt.out"
            base_out.write_text("ORCA TERMINATED NORMALLY\n")

            def execute(inp, **kwargs):
                inp = Path(inp)
                out = inp.with_suffix(".out")
                if inp.name == "freq.inp":
                    out.write_text(
                        "VIBRATIONAL FREQUENCIES\n  0: 0.00\n  1: 0.00\n"
                        "  2: 0.00\n  3: 0.00\n  4: 0.00\n  5: 0.00\n"
                        "  6: 100.00\nNORMAL MODES\nORCA TERMINATED NORMALLY\n"
                    )
                    return {"status": "completed", "out": out, "output": str(out), "resumed": False}
                write_xyz(
                    inp.with_suffix(".xyz"), ["H"], np.asarray([[0.01, 0.0, 0.0]]), "recovered",
                )
                attempt = int(inp.parent.name.split("-")[-1])
                executed_input = inp
                if attempt == 2:
                    executed_input = inp.with_name("opt-execution-retry-001.inp")
                    executed_input.write_text(inp.read_text())
                    write_xyz(
                        executed_input.with_suffix(".xyz"), ["H"],
                        np.asarray([[0.02, 0.0, 0.0]]), "append-only retry",
                    )
                    out = executed_input.with_suffix(".out")
                out.write_text(
                    ("THE OPTIMIZATION HAS CONVERGED\n" if attempt == 2 else "")
                    + "ORCA TERMINATED NORMALLY\n"
                )
                return {
                    "status": "completed", "input": str(executed_input),
                    "out": out, "output": str(out), "resumed": False,
                }

            run_or_resume.side_effect = execute
            result = _run_stationary_point_recovery(
                {
                    "optimized_xyz": source,
                    "optimization": {"out": base_out},
                    "frequency_check": {"IsMinimum": False},
                },
                folder,
                charge=0,
                multiplicity=2,
                ncores=1,
                method_keywords=None,
                timeout_seconds=1.0,
                policy=SurfaceResolutionPolicy(maximum_stationary_point_recoveries=2),
            )
            manifest = json.loads(
                (folder / "stationary-point-recovery" / "stationary-point-recovery.json").read_text()
            )
            recovery_evidence = json.loads(Path(
                manifest["attempts"][1]["surface_resolution_attempt_evidence"]["path"]
            ).read_text())

        self.assertTrue(result["frequency_check"]["IsMinimum"])
        self.assertEqual(manifest["status"], "validated_minimum")
        self.assertEqual(
            [item["status"] for item in manifest["attempts"]],
            ["intermediate_optimization_unconverged", "validated_minimum"],
        )
        self.assertTrue(manifest["attempts"][0]["claims_prohibited"])
        self.assertTrue(str(result["optimized_xyz"]).endswith("opt-execution-retry-001.xyz"))
        self.assertEqual(recovery_evidence["strategy"], "stationary_point_recovery")
        self.assertEqual(recovery_evidence["status"], "validated_minimum")
        self.assertIn("process_leases", recovery_evidence)

    @patch("storca.reaction_path._verify_ts_seed")
    @patch("storca.reaction_path._scan_energies")
    @patch("storca.reaction_path._scan_frames")
    @patch("storca.reaction_path._run_or_resume")
    def test_guided_scan_seed_still_requires_ts_frequency_and_irc(
        self, run_or_resume, scan_frames, scan_energies, verify_ts,
    ):
        route = RouteSpec(
            route_id="guided",
            source="test",
            parent_smiles="CCO",
            reactant_smiles=("CCO",),
            product_smiles=("C", "C=O"),
            charge=0,
            multiplicity=1,
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            reactant = write_xyz(
                folder / "reactant.xyz", ["C", "C", "O"],
                np.asarray([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.8, 0.0, 0.0]]),
                "reactant",
            )
            product = write_xyz(
                folder / "product.xyz", ["C", "C", "O"],
                np.asarray([[0.0, 0.0, 0.0], [2.8, 0.0, 0.0], [3.9, 0.0, 0.0]]),
                "product",
            )
            frames = []
            for index in range(3):
                frames.append(write_xyz(
                    folder / f"scan-{index}.xyz", ["C", "C", "O"],
                    np.asarray([[0.0, 0.0, 0.0], [1.5 + 0.5 * index, 0.0, 0.0],
                                [2.8 + 0.3 * index, 0.0, 0.0]]), f"frame {index}",
                ))
            run_or_resume.return_value = {"status": "completed", "output": str(folder / "scan.out")}
            scan_frames.return_value = frames
            scan_energies.return_value = [-100.0, -99.8, -100.1]
            verify_ts.return_value = {
                "status": "verified_transition_state_path",
                "path_classification": "verified_barriered_path",
                "transition_state": {"frequency_check": {"NumImaginary": 1}},
                "irc": {"endpoint_validation": {"valid": True}},
                "trajectory": "irc.xyz",
                "path_evidence": {"kind": "orca_ts_frequency_irc"},
            }
            result = _run_bond_change_guided_ts_recovery(
                route,
                reactant_xyz=reactant,
                product_xyz=product,
                bond_changes={
                    "broken_bonds": [[0, 1]],
                    "bond_order_changes": [{"atom_indices": [1, 2]}],
                },
                mapping_result={"valid": True},
                path_dir=folder / "path",
                separated_fragments=True,
                ncores=1,
                method_keywords=None,
                timeout_seconds=1.0,
                surface_resolution_policy=SurfaceResolutionPolicy(maximum_scan_points=3),
            )
            guided_evidence = json.loads(Path(
                result["surface_resolution_attempt_evidence"]["path"]
            ).read_text())

        self.assertEqual(result["status"], "verified_transition_state_path")
        self.assertEqual(result["ts_seed_selection"]["index"], 1)
        verify_ts.assert_called_once()
        self.assertEqual(guided_evidence["strategy"], "bond_change_guided_scan")
        self.assertEqual(guided_evidence["claims_permitted"], ["verified_transition_state_path"])

    @patch("storca.reaction_path.parse_orca_normal_modes")
    def test_minimum_recovery_seed_follows_most_negative_mode(self, parse_modes):
        parse_modes.return_value = NormalModeSet(
            frequencies_cm_1=np.asarray([-120.0, -35.0, 100.0]),
            vectors=np.asarray([
                [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
                [[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
                [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
            ]),
            atom_labels=("H", "H"),
            coordinates_bohr=np.zeros((2, 3)),
        )
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            source = write_xyz(
                folder / "source.xyz", ["H", "H"],
                np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), "source",
            )
            detail = _imaginary_mode_displacement_seed(
                source, folder / "freq.hess", folder / "displaced.xyz", direction=1,
            )
            _, displaced, _ = read_xyz(folder / "displaced.xyz")

        self.assertEqual(detail["mode_index"], 0)
        self.assertEqual(detail["frequency_cm-1"], -120.0)
        self.assertTrue(np.allclose(displaced[:, 0], [0.15, 0.85]))
        self.assertTrue(np.allclose(displaced[:, 1:], 0.0))

    @patch("storca.reaction_path.resolve_validated_intermediate_identity")
    @patch("storca.reaction_path._optimize_endpoint")
    def test_validated_neb_intermediate_prepares_both_segments(self, optimize, resolve_identity):
        route = self._association_route()
        resolve_identity.return_value = {
            "status": "resolved",
            "native_candidates": [{"candidate_id": "segment-1"}, {"candidate_id": "segment-2"}],
        }
        with tempfile.TemporaryDirectory() as temp:
            path_dir = Path(temp)
            write_xyz(
                path_dir / "reactant.xyz", ["H"], np.asarray([[0.0, 0.0, 0.0]]), "reactant",
            )
            write_xyz(
                path_dir / "neb-ts_im3.xyz", ["H"], np.asarray([[0.0, 0.0, 0.5]]), "intermediate",
            )
            optimized = write_xyz(
                path_dir / "intermediate.xyz", ["H"], np.asarray([[0.0, 0.0, 0.4]]), "optimized",
            )
            output = path_dir / "neb-ts.out"
            output.write_text("Possible intermediate minimum found at image(s): 3\n")
            optimize.return_value = {
                "optimized_xyz": optimized,
                "optimization": {"out": path_dir / "opt.out"},
                "frequency": {"out": path_dir / "freq.out"},
                "frequency_check": {"IsMinimum": True},
            }

            result = _validate_neb_intermediates(
                route, neb_output=output, path_dir=path_dir, ncores=1,
                method_keywords=None, timeout_seconds=1.0,
            )

        self.assertEqual(result["status"], "validated_intermediate_detected")
        self.assertEqual(result["validated_candidate_count"], 1)
        self.assertEqual(
            [
                {"from": item["from"], "to": item["to"]}
                for item in result["segmented_verification"]["segments"]
            ],
            [
                {"from": "declared_reactants", "to": str(optimized)},
                {"from": str(optimized), "to": "declared_products"},
            ],
        )
        self.assertEqual(
            [item["native_candidate"]["candidate_id"] for item in result["segmented_verification"]["segments"]],
            ["segment-1", "segment-2"],
        )

    @patch("storca.reaction_path.resolve_validated_intermediate_identity")
    @patch("storca.reaction_path._optimize_endpoint")
    def test_neb_warning_that_relaxes_to_endpoint_does_not_block_path(self, optimize, resolve_identity):
        route = self._association_route()
        resolve_identity.return_value = {
            "status": "endpoint_return",
            "status_detail": "matches_declared_reactants",
            "matched_endpoint_side": "reactants",
            "native_candidates": [],
        }
        with tempfile.TemporaryDirectory() as temp:
            path_dir = Path(temp)
            write_xyz(
                path_dir / "reactant.xyz", ["H"], np.asarray([[0.0, 0.0, 0.0]]), "reactant",
            )
            write_xyz(
                path_dir / "neb-ts_im3.xyz", ["H"], np.asarray([[0.0, 0.0, 0.5]]), "warning",
            )
            optimized = write_xyz(
                path_dir / "endpoint.xyz", ["H"], np.asarray([[0.0, 0.0, 0.0]]), "endpoint",
            )
            output = path_dir / "neb-ts.out"
            output.write_text("Possible intermediate minimum found at image(s): 3\n")
            optimize.return_value = {
                "optimized_xyz": optimized,
                "optimization": {"out": path_dir / "opt.out"},
                "frequency": {"out": path_dir / "freq.out"},
                "frequency_check": {"IsMinimum": True},
            }

            result = _validate_neb_intermediates(
                route, neb_output=output, path_dir=path_dir, ncores=1,
                method_keywords=None, timeout_seconds=1.0,
            )

        self.assertEqual(result["status"], "intermediate_warning_dismissed_endpoint_return")
        self.assertEqual(result["endpoint_return_candidate_count"], 1)
        self.assertEqual(result["identity_unresolved_candidate_count"], 0)
        self.assertEqual(result["segmented_verification"]["status"], "not_prepared")

    def test_orca_intermediate_image_diagnostic_is_parsed(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "neb-ts.out"
            output.write_text(
                "Possible intermediate minimum found at image(s): 1 4\n"
                "Possible intermediate minimum found at image(s): 4 7\n"
            )
            self.assertEqual(_neb_intermediate_image_indices(output), [1, 4, 7])

    def test_orca_angle_warning_ignores_corrupt_first_integer(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "neb-ts.out"
            output.write_text(
                "Warning: maximum angle between images (-1659107148 - 5 - 6): 0.00\n"
                "Warning: maximum angle between images (745032783 - 5 - 6): 0.00\n"
            )
            self.assertEqual(_neb_angle_warning_image_indices(output), [5, 6])

    @patch("storca.reaction_path.run_orca")
    def test_neb_intermediate_warning_stops_and_is_retained(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "neb-ts.inp"
            inp.write_text(
                "! B3LYP def2-SVP NEB-TS\n"
                '%neb\n  Product "product.xyz"\nend\n'
                "* xyzfile 0 1 reactant.xyz\n"
            )
            (folder / "reactant.xyz").write_text("1\nreactant\nH 0 0 0\n")
            (folder / "product.xyz").write_text("1\nproduct\nH 0 0 1\n")
            output = inp.with_suffix(".out")

            def stopped(input_path, **kwargs):
                output.write_text(
                    "Possible intermediate minimum found at image(s): 2\n"
                )
                self.assertFalse(kwargs["live_output"])
                return {"out": output, "stopped_early": True}

            run_orca.side_effect = stopped
            result = _run_or_resume(inp, timeout_seconds=10.0)
            self.assertEqual(result["status"], "completed_intermediate_detected")
            self.assertEqual(result["intermediate_image_indices"], [2])
            self.assertEqual(result["early_stop"]["status"], "stopped_after_orca_intermediate_warning")

    @patch("storca.reaction_path.run_orca")
    def test_persistent_monotonic_endpoint_path_stops_early(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "neb-ts.inp"
            inp.write_text(
                "! B3LYP def2-SVP NEB-TS\n"
                '%neb\n  Product "product.xyz"\nend\n'
                "* xyzfile 0 1 reactant.xyz\n"
            )
            (folder / "reactant.xyz").write_text("1\nreactant\nH 0 0 0\n")
            (folder / "product.xyz").write_text("1\nproduct\nH 0 0 1\n")
            output = inp.with_suffix(".out")

            def stopped(input_path, **kwargs):
                lines = [
                    f"Warning: path is strictly increasing in energy   LBFGS {index:5d}      7    0.08\n"
                    for index in range(12)
                ]
                decisions = [kwargs["stream_output"](line) for line in lines]
                output.write_text("".join(lines))
                self.assertFalse(any(decisions[:-1]))
                self.assertTrue(decisions[-1])
                return {"out": output, "stopped_early": True}

            run_orca.side_effect = stopped
            result = _run_or_resume(inp, timeout_seconds=10.0)

        self.assertEqual(result["status"], "completed_path_quality_warning")
        self.assertEqual(
            result["early_stop"]["status"],
            "stopped_after_persistent_monotonic_endpoint_maximum",
        )
        self.assertEqual(result["path_quality_warning"]["highest_energy_image"], 7)

    def test_monotonic_warning_parser_requires_a_trailing_streak(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "neb-ts.out"
            output.write_text("".join(
                f"Warning: path is strictly increasing in energy LBFGS {index} 7 0.08\n"
                for index in range(12)
            ) + "   LBFGS    12      4    0.06\n")
            result = _neb_monotonic_endpoint_diagnostic(output)
        self.assertFalse(result["valid"])
        self.assertEqual(result["maximum_iteration_count"], 12)

    @patch("storca.reaction_path.run_orca")
    def test_persistent_endpoint_highest_image_stops_as_unresolved_allocation(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "neb-ts.inp"
            inp.write_text(
                "! B3LYP def2-SVP NEB-TS\n"
                '%neb\n  Product "product.xyz"\n  NImages 6\nend\n'
                "* xyzfile 0 1 reactant.xyz\n"
            )
            (folder / "reactant.xyz").write_text("1\nreactant\nH 0 0 0\n")
            (folder / "product.xyz").write_text("1\nproduct\nH 0 0 1\n")
            output = inp.with_suffix(".out")

            def stopped(input_path, **kwargs):
                lines = ["Switch-on CI threshold               0.020000\n"]
                decisions = [kwargs["stream_output"](lines[0])]
                for index in range(24):
                    force = 0.010 if index < 8 else 0.030
                    line = (
                        f"   LBFGS {index:5d}      7    0.080000    "
                        f"{force:.6f}   0.001000  10.0000\n"
                    )
                    lines.append(line)
                    decisions.append(kwargs["stream_output"](line))
                output.write_text("".join(lines))
                self.assertFalse(any(decisions[:-1]))
                self.assertTrue(decisions[-1])
                return {"out": output, "stopped_early": True}

            run_orca.side_effect = stopped
            result = _run_or_resume(inp, timeout_seconds=10.0)
            diagnostic = _neb_persistent_endpoint_highest_image_diagnostic(
                output, inp,
            )

        self.assertTrue(diagnostic["valid"])
        self.assertEqual(result["status"], "completed_path_quality_warning")
        self.assertEqual(
            result["early_stop"]["status"],
            "stopped_after_persistent_endpoint_highest_image",
        )
        self.assertIn("does not establish a barrierless path", diagnostic["interpretation"])

    @patch("storca.reaction_path.run_orca")
    def test_submerged_climbing_image_is_a_neutral_path_reference_warning(self, run_orca):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "neb-ts.inp"
            inp.write_text(
                "! B3LYP def2-SVP NEB-TS\n"
                '%neb\n  Product "product.xyz"\n  NImages 6\nend\n'
                "* xyzfile 0 1 reactant.xyz\n"
            )
            (folder / "reactant.xyz").write_text("1\nreactant\nH 0 0 0\n")
            (folder / "product.xyz").write_text("1\nproduct\nH 0 0 1\n")
            output = inp.with_suffix(".out")
            prefix = (
                "Properties of fixed NEB end points:\n"
                "  Reactant:\n    E .... -100.000000 Eh\n"
                "  Product:\n    E ....  -99.920000 Eh\n"
                "Optim.  Iteration  CI   E(CI)-E(0)   max(|Fp|)   RMS(Fp)    "
                "dS     max(|FCI|)   RMS(FCI)\n"
            )

            def stopped(input_path, **kwargs):
                output.write_text(prefix)
                decisions = []
                for index in range(12):
                    line = (
                        f"   LBFGS {index:5d}      3    0.060000    0.010000   "
                        "0.002000  10.0000    0.008000    0.003000\n"
                    )
                    with output.open("a") as stream:
                        stream.write(line)
                    decisions.append(kwargs["stream_output"](line))
                self.assertFalse(any(decisions[:-1]))
                self.assertTrue(decisions[-1])
                return {"out": output, "stopped_early": True}

            run_orca.side_effect = stopped
            result = _run_or_resume(inp, timeout_seconds=10.0)
            diagnostic = _neb_submerged_climbing_image_diagnostic(output)

        self.assertTrue(diagnostic["valid"])
        self.assertAlmostEqual(diagnostic["product_relative_energy_hartree"], 0.08)
        self.assertEqual(result["status"], "completed_path_quality_warning")
        self.assertEqual(
            result["early_stop"]["status"],
            "stopped_after_persistent_submerged_climbing_image",
        )
        self.assertIn("finite encounter energies are not", diagnostic["interpretation"])
        self.assertIn("No mechanism or rate is established", diagnostic["interpretation"])

    def test_surface_budget_counts_process_leases_and_blocks_the_next_job(self):
        policy = SurfaceResolutionPolicy(
            maximum_calculations=2,
            maximum_core_hours=1.0,
            maximum_wall_seconds=100.0,
        )
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contract = _initialize_surface_budget(root, policy, now_epoch=1000.0)
            input_path = root / "job.inp"
            input_path.write_text("%pal\n nprocs 4\nend\n")
            first = _surface_budget_snapshot(
                input_path,
                requested_calculations=2,
                requested_ncores=4,
                requested_timeout_seconds=80.0,
                now_epoch=1001.0,
            )
            for index in range(2):
                (root / f"job-{index}.process.json").write_text(json.dumps({
                    "kind": "storca_external_process_lease",
                    "status": "terminated",
                    "input": str(input_path),
                    "ncores": 4,
                    "started_at_epoch": 1001.0,
                    "ended_at_epoch": 1002.0,
                }))
            exhausted = _surface_budget_snapshot(
                input_path,
                requested_calculations=1,
                requested_ncores=4,
                requested_timeout_seconds=80.0,
                now_epoch=1003.0,
            )
            snapshots = list(
                (Path(contract["contract_path"]).parent / "snapshots").glob("*.json")
            )

        self.assertEqual(first["status"], "within_budget")
        self.assertLessEqual(first["effective_timeout_seconds_per_calculation"], 49.5)
        self.assertEqual(exhausted["status"], "surface_resolution_budget_exhausted")
        self.assertIn("maximum_calculations", exhausted["exceeded"])
        self.assertEqual(exhausted["consumed"]["calculation_count"], 2)
        self.assertEqual(len(snapshots), 2)

    def test_surface_budget_policy_change_is_append_only_and_mutation_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = _initialize_surface_budget(
                root, SurfaceResolutionPolicy(maximum_calculations=4), now_epoch=10.0,
            )
            second = _initialize_surface_budget(
                root, SurfaceResolutionPolicy(maximum_calculations=5), now_epoch=20.0,
            )
            self.assertNotEqual(first["contract_path"], second["contract_path"])
            self.assertTrue(Path(first["contract_path"]).is_file())
            retained = json.loads(Path(second["contract_path"]).read_text())
            retained["policy"]["maximum_calculations"] = 999
            Path(second["contract_path"]).write_text(json.dumps(retained))
            with self.assertRaisesRegex(RuntimeError, "mutated"):
                _initialize_surface_budget(
                    root, SurfaceResolutionPolicy(maximum_calculations=5), now_epoch=30.0,
                )

    @patch("storca.reaction_path._run_or_resume")
    def test_separated_intermediate_segments_require_channel_specific_solvers(self, run_or_resume):
        route = RouteSpec(
            route_id="two-to-two",
            source="test",
            parent_smiles="N=O",
            reactant_smiles=("N=O", "N=O"),
            product_smiles=("[N]=O", "N[O]"),
            multiplicity=1,
            reactant_multiplicities=(1, 1),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            reactant = folder / "reactant.xyz"
            product = folder / "product.xyz"
            intermediate = folder / "intermediate.xyz"
            geometry = "6\ntest\nN 0 0 0\nO 0 0 1\nH 0 1 0\nN 3 0 0\nO 3 0 1\nH 3 1 0\n"
            for path in (reactant, product, intermediate):
                path.write_text(geometry)
            result = _run_segmented_intermediate_nebs(
                route,
                reactant_xyz=reactant,
                product_xyz=product,
                intermediate_evidence={"candidates": [{
                    "status": "validated_intermediate_minimum",
                    "optimized_xyz": str(intermediate),
                    "intermediate_identity": {
                        "status": "resolved",
                        "native_candidates": [
                            {"candidate_id": "reactants-to-intermediate"},
                            {"candidate_id": "intermediate-to-products"},
                        ],
                    },
                }]},
                path_dir=folder / "path",
                separated_fragments=True,
                ncores=1,
                method_keywords=None,
                timeout_seconds=1.0,
                neb_images=8,
            )
        run_or_resume.assert_not_called()
        self.assertEqual(result["status"], "segment_rate_model_unavailable")
        self.assertEqual(result["channel_specific_segment_count"], 2)
        self.assertEqual(
            [item["channel_kind"] for item in result["segments"]],
            ["association", "dissociation"],
        )
        self.assertTrue(all(
            item["surface_resolution_attempt_evidence"]["strategy"]
            == "segmented_path_verification"
            for item in result["segments"]
        ))

    @patch("storca.reaction_path._verify_ts_seed")
    @patch("storca.reaction_path._locate_neb_ts")
    @patch("storca.reaction_path._run_or_resume")
    def test_bound_intermediate_segments_each_require_ts_frequency_and_irc(
        self, run_or_resume, locate_ts, verify_ts,
    ):
        route = RouteSpec(
            route_id="c3h8o",
            source="test",
            parent_smiles="CCCO",
            reactant_smiles=("[CH3:1][CH2:2][CH2:3][OH:4]",),
            product_smiles=("[CH3:1][CH2:2][O:4][CH3:3]",),
            charge=0,
            multiplicity=1,
            reactant_charges=(0,),
            product_charges=(0,),
            reactant_multiplicities=(1,),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            geometry = (
                "4\ntest\nC 0 0 0\nC 1.5 0 0\nC 2.8 0 0\nO 4.0 0 0\n"
            )
            reactant = folder / "reactant.xyz"
            product = folder / "product.xyz"
            intermediate = folder / "intermediate.xyz"
            ts = folder / "ts.xyz"
            for path in (reactant, product, intermediate, ts):
                path.write_text(geometry)
            reactant_frequency = folder / "reactant-frequency.out"
            intermediate_frequency = folder / "intermediate-frequency.out"
            product_frequency = folder / "product-frequency.out"
            ts_frequency = folder / "ts-frequency.out"
            for path in (
                reactant_frequency, intermediate_frequency, product_frequency, ts_frequency,
            ):
                path.write_text("retained test frequency evidence\n")
            run_or_resume.return_value = {
                "status": "completed", "input": str(folder / "neb-ts.inp"),
                "output": str(folder / "neb-ts.out"),
            }
            locate_ts.return_value = ts
            verify_ts.return_value = {
                "status": "verified_transition_state_path",
                "path_classification": "verified_barriered_path",
                "transition_state": {
                    "frequency_check": {"NumImaginary": 1},
                    "frequency_output": str(ts_frequency),
                },
                "irc": {"endpoint_validation": {"valid": True}},
            }
            result = _run_segmented_intermediate_nebs(
                route,
                reactant_xyz=reactant,
                product_xyz=product,
                intermediate_evidence={"candidates": [{
                    "status": "validated_intermediate_minimum",
                    "optimized_xyz": str(intermediate),
                    "frequency_output": str(intermediate_frequency),
                    "intermediate_identity": self._bound_segment_identity(),
                }]},
                path_dir=folder / "path",
                separated_fragments=False,
                ncores=1,
                method_keywords=None,
                timeout_seconds=1.0,
                neb_images=6,
                parent_endpoint_stationary={
                    "reactants": [{
                        "label": "propan-1-ol", "original_label": "propan-1-ol",
                        "smiles": "[CH3:1][CH2:2][CH2:3][OH:4]",
                        "charge": 0, "multiplicity": 1,
                        "status": "validated_minimum",
                        "frequency_output": str(reactant_frequency),
                    }],
                    "products": [{
                        "label": "methoxyethane", "original_label": "methoxyethane",
                        "smiles": "[CH3:1][CH2:2][O:4][CH3:3]",
                        "charge": 0, "multiplicity": 1,
                        "status": "validated_minimum",
                        "frequency_output": str(product_frequency),
                    }],
                },
            )

        self.assertEqual(result["status"], "verified_segments_rate_generation_required")
        self.assertEqual(result["completed_segment_count"], 2)
        self.assertEqual(verify_ts.call_count, 2)
        self.assertFalse(result["rate_claim_allowed"])
        self.assertTrue(all(item["segment_mapping"]["valid"] for item in result["segments"]))
        self.assertTrue(all(
            item["surface_resolution_attempt_evidence"]["claims_permitted"]
            == ["verified_transition_state_path"]
            for item in result["segments"]
        ))

    @patch("storca.reaction_path._run_or_resume")
    def test_segment_endpoint_pair_recursion_is_prevented(self, run_or_resume):
        identity = self._bound_segment_identity()
        visited = {
            item["candidate_id"] for item in identity["native_candidates"]
        }
        route = RouteSpec(
            route_id="c3h8o", source="test", parent_smiles="CCCO",
            reactant_smiles=("CCCO",), product_smiles=("CCOC",),
            charge=0, multiplicity=1,
            reactant_multiplicities=(1,), product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            geometry = "4\ntest\nC 0 0 0\nC 1 0 0\nC 2 0 0\nO 3 0 0\n"
            paths = [folder / name for name in ("reactant.xyz", "product.xyz", "intermediate.xyz")]
            for path in paths:
                path.write_text(geometry)
            result = _run_segmented_intermediate_nebs(
                route,
                reactant_xyz=paths[0], product_xyz=paths[1],
                intermediate_evidence={"candidates": [{
                    "status": "validated_intermediate_minimum",
                    "optimized_xyz": str(paths[2]),
                    "intermediate_identity": identity,
                }]},
                path_dir=folder / "path", separated_fragments=False,
                ncores=1, method_keywords=None, timeout_seconds=1.0, neb_images=6,
                visited_endpoint_pair_hashes=visited,
            )

        run_or_resume.assert_not_called()
        self.assertEqual(result["status"], "segment_path_unresolved")
        self.assertTrue(all(
            item["reason"] == "visited_endpoint_pair_recursion_prevented"
            for item in result["segments"]
        ))

    @staticmethod
    def _association_route() -> RouteSpec:
        return RouteSpec(
            route_id="association",
            source="test",
            parent_smiles="[CH3]",
            reactant_smiles=("[CH3]", "[OH]"),
            product_smiles=("CO",),
            multiplicity=1,
            reactant_multiplicities=(2, 2),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )

    def test_unique_multifragment_association_maps_and_builds_orientations(self):
        route = self._association_route()
        mapping = resolve_route_atom_mapping(route)
        self.assertTrue(mapping["valid"])
        self.assertEqual(mapping["surface"]["common_total_multiplicities"], [1])
        with tempfile.TemporaryDirectory() as temp:
            seeds = build_endpoint_complex_seeds(route, Path(temp), mapping_result=mapping)
            self.assertEqual(seeds["orientation_count"], 3)
            self.assertEqual(seeds["bond_changes"]["formed_bonds"], [[0, 4]])
            first = seeds["orientations"][0]
            left, _, _ = read_xyz(Path(first["reactant_xyz"]))
            right, _, _ = read_xyz(Path(first["product_xyz"]))
            self.assertEqual(left, right)

    def test_automatic_mapping_canonicalizes_true_endpoint_symmetry(self):
        route = RouteSpec(
            route_id="symmetric-identity",
            source="test",
            parent_smiles="CC",
            reactant_smiles=("CC",),
            product_smiles=("CC",),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["source"], "automatic_minimum_graph_edit_symmetry_canonicalized")
        self.assertEqual(result["symmetry_resolution"]["equivalent_mapping_count"], 2)

    def test_automatic_mapping_canonicalizes_o2_symmetry_for_nno_initiation(self):
        route = RouteSpec(
            route_id="nno-o2-initiation",
            source="test",
            parent_smiles="NNO",
            reactant_smiles=("[O][O]", "NNO"),
            product_smiles=("[O]O", "NN[O]"),
            multiplicity=3,
            reactant_multiplicities=(3, 1),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["surface"]["selected_multiplicity"], 3)
        self.assertEqual(result["source"], "automatic_minimum_graph_edit_symmetry_canonicalized")
        self.assertEqual(result["symmetry_resolution"]["equivalent_mapping_count"], 2)
        self.assertEqual(len(result["all_atom_mapping"]), 9)

    def test_automatic_mapping_canonicalizes_identical_reactant_exchange(self):
        """Two identical reactant molecules may exchange without ambiguity."""
        route = RouteSpec(
            route_id="hno-disproportionation",
            source="test",
            parent_smiles="N=O",
            reactant_smiles=("N=O", "N=O"),
            product_smiles=("[N]=O", "N[O]"),
            multiplicity=1,
            reactant_multiplicities=(1, 1),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["source"], "automatic_minimum_graph_edit_symmetry_canonicalized")
        self.assertEqual(result["symmetry_resolution"]["equivalent_mapping_count"], 2)

    def test_large_symmetric_radical_recombination_does_not_hit_mapping_storage_cap(self):
        """Equivalent methyl/copy permutations are one chemical mapping class."""
        route = RouteSpec(
            route_id="tert-butoxy-recombination",
            source="rmg",
            parent_smiles="CC(C)(C)[O]",
            reactant_smiles=("CC(C)(C)[O]", "CC(C)(C)[O]"),
            product_smiles=("CC(C)(C)OOC(C)(C)C",),
            multiplicity=1,
            reactant_multiplicities=(2, 2),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["source"], "automatic_minimum_graph_edit_symmetry_canonicalized")
        self.assertGreater(result["symmetry_resolution"]["equivalent_mapping_count"], 64)
        self.assertEqual(result["symmetry_resolution"]["symmetry_class_count"], 1)

    def test_hno_water_route_canonicalizes_equivalent_hydrogen_transfers(self):
        route = RouteSpec(
            route_id="hno-water-disproportionation",
            source="test",
            parent_smiles="N=O",
            reactant_smiles=("N=O", "N=O"),
            product_smiles=("[N-]=[N+]=O", "O"),
            multiplicity=1,
            reactant_multiplicities=(1, 1),
            product_multiplicities=(1, 1),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(len(result["hydrogen_transfers"]), 2)
        self.assertEqual(
            result["hydrogen_symmetry_resolution"]["status"],
            "equivalent_hydrogen_permutations_canonicalized",
        )
        self.assertEqual(len({item["product_acceptor_index"] for item in result["hydrogen_transfers"]}), 1)

    def test_multiple_distinct_hydrogen_donors_and_acceptors_remain_ambiguous(self):
        route = RouteSpec(
            route_id="distinct-two-centre-transfer",
            source="test",
            parent_smiles="CO",
            reactant_smiles=("CO", "CO", "[NH2]", "[NH2]"),
            product_smiles=("C[O]", "C[O]", "N", "N"),
            multiplicity=1,
            reactant_multiplicities=(1, 1, 2, 2),
            product_multiplicities=(2, 2, 1, 1),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertFalse(result["valid"], result)
        self.assertEqual(result["status"], "ambiguous_hydrogen_transfer_mapping")

    def test_automatic_mapping_keeps_distinct_isomerization_maps_ambiguous(self):
        route = RouteSpec(
            route_id="non-equivalent-maps",
            source="test",
            parent_smiles="CCCO",
            reactant_smiles=("CCCO",),
            product_smiles=("CC(C)O",),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertFalse(result["valid"])
        self.assertEqual(result["status"], "ambiguous_atom_mapping")

    def test_supplied_mapping_extends_across_hydrogen_transfer(self):
        route = RouteSpec(
            route_id="abstraction",
            source="test",
            parent_smiles="CCO",
            reactant_smiles=("CCO", "[O][O]"),
            product_smiles=("CC[O]", "[O]O"),
            multiplicity=3,
            reactant_multiplicities=(1, 3),
            product_multiplicities=(2, 2),
            atom_mapping=((0, 0), (1, 1), (2, 2), (3, 3), (4, 4)),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"])
        self.assertEqual(len(result["hydrogen_transfers"]), 1)
        self.assertEqual(result["surface"]["common_total_multiplicities"], [3])

    def test_standalone_hydrogen_fragment_uses_all_hydrogen_mapping_layer(self):
        route = RouteSpec(
            route_id="nno-h-loss",
            source="test",
            parent_smiles="NNO",
            reactant_smiles=("NNO",),
            product_smiles=("[H]", "[NH]NO"),
            reactant_labels=("stability(1)",),
            product_labels=("H(2)", "radical(3)"),
            multiplicity=1,
            reactant_multiplicities=(1,),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["reactant_skeleton_atom_count"], 3)
        self.assertEqual(result["product_skeleton_atom_count"], 3)
        self.assertEqual(result["explicit_atom_count"], 7)
        self.assertEqual(len(result["all_atom_mapping"]), 7)
        self.assertEqual(len(result["hydrogen_transfers"]), 1)
        transfer = result["hydrogen_transfers"][0]
        self.assertEqual(transfer["reactant_donor_index"], 0)
        self.assertIsNone(transfer["product_acceptor_index"])
        self.assertEqual(transfer["product_fragment_index"], 0)
        self.assertEqual(result["surface"]["charge"], 0)
        self.assertEqual(result["surface"]["selected_multiplicity"], 1)

        with tempfile.TemporaryDirectory() as temp:
            seeds = build_endpoint_complex_seeds(route, Path(temp), mapping_result=result)
            self.assertEqual(seeds["orientation_count"], 3)
            first = seeds["orientations"][0]
            reactant_elements, _, _ = read_xyz(Path(first["reactant_xyz"]))
            product_elements, _, _ = read_xyz(Path(first["product_xyz"]))
            self.assertEqual(reactant_elements, product_elements)
            self.assertEqual(len(product_elements), 7)

            prepared = run_generic_reaction_path(
                route, Path(temp) / "prepared", prepare_only=True,
            )
            h_record = prepared["route_verification"]["stationary_points"]["products"][0]
            self.assertEqual(read_xyz(Path(h_record["input_xyz"]))[0], ["H"])

    def test_standalone_charged_hydrogen_preserves_fragment_state(self):
        for smiles, charge in (("[H+]", 1), ("[H-]", -1)):
            with self.subTest(smiles=smiles):
                route = RouteSpec(
                    route_id=f"identity-{charge}",
                    source="test",
                    parent_smiles=smiles,
                    reactant_smiles=(smiles,),
                    product_smiles=(smiles,),
                    charge=charge,
                    multiplicity=1,
                    reactant_charges=(charge,),
                    product_charges=(charge,),
                    reactant_multiplicities=(1,),
                    product_multiplicities=(1,),
                    protocol="generic_endpoint_path",
                )
                result = resolve_route_atom_mapping(route)
                self.assertTrue(result["valid"], result)
                self.assertEqual(result["reactant_skeleton_atom_count"], 0)
                self.assertEqual(result["product_skeleton_atom_count"], 0)
                self.assertEqual(result["explicit_atom_count"], 1)
                self.assertEqual(result["all_atom_mapping"], [[0, 0]])
                self.assertEqual(result["surface"]["charge"], charge)
                self.assertEqual(result["surface"]["selected_multiplicity"], 1)

    def test_standalone_hydrogen_isotope_imbalance_fails_closed(self):
        route = RouteSpec(
            route_id="isotope-imbalance",
            source="test",
            parent_smiles="[2H]",
            reactant_smiles=("[2H]",),
            product_smiles=("[H]",),
            multiplicity=2,
            reactant_multiplicities=(2,),
            product_multiplicities=(2,),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertFalse(result["valid"])
        self.assertEqual(result["status"], "unbalanced_elements")

    def test_connectivity_gate_distinguishes_bound_and_separated_hydrogen(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            bound = write_xyz(
                folder / "bound.xyz", ["H", "H"],
                np.asarray([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]),
                "mock optimized H2",
            )
            separated = write_xyz(
                folder / "separated.xyz", ["H", "H"],
                np.asarray([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
                "mock separated H atoms",
            )
            self.assertTrue(validate_declared_connectivity(bound, ("[H][H]",))["valid"])
            self.assertTrue(validate_declared_connectivity(separated, ("[H]", "[H]"))["valid"])
            collapsed = validate_declared_connectivity(bound, ("[H]", "[H]"))
            dissociated = validate_declared_connectivity(separated, ("[H][H]",))
            self.assertFalse(collapsed["valid"])
            self.assertEqual(collapsed["status"], "connectivity_changed")
            self.assertEqual(collapsed["unexpected_covalent_bonds"][0]["atom_indices"], [0, 1])
            self.assertFalse(dissociated["valid"])
            self.assertEqual(dissociated["missing_declared_bonds"][0]["atom_indices"], [0, 1])

    @patch("storca.reaction_path._optimize_endpoint")
    def test_collapsed_product_complex_is_rejected_before_neb(self, optimize):
        route = self._association_route()
        mapping = resolve_route_atom_mapping(route)
        with tempfile.TemporaryDirectory() as temp:
            seeds = build_endpoint_complex_seeds(
                route, Path(temp) / "seeds", orientations=1, mapping_result=mapping,
            )
            reactant_seed = Path(seeds["orientations"][0]["reactant_xyz"])
            optimize.return_value = {
                # Both optimizations are mocked as the separated reactant
                # complex: the declared bonded product has collapsed to the
                # wrong endpoint channel.
                "optimized_xyz": reactant_seed,
                "optimization": {"out": Path(temp) / "opt.out"},
                "frequency": {"out": Path(temp) / "freq.out"},
                "frequency_check": {"IsMinimum": True},
            }
            result = _run_generic_orientation(
                route,
                seeds["orientations"][0],
                mapping_result=mapping,
                ncores=1,
                method_keywords=None,
                timeout_seconds=1.0,
                neb_images=4,
                barrierless_classifier=None,
            )
        self.assertEqual(result["status"], "endpoint_connectivity_validation_failed")
        self.assertTrue(result["endpoint_validation"]["reactant"]["connectivity_check"]["valid"])
        self.assertFalse(result["endpoint_validation"]["product"]["connectivity_check"]["valid"])

    @patch("storca.reaction_path._optimize_endpoint")
    def test_separated_stationary_species_must_retain_connectivity(self, optimize):
        route = RouteSpec(
            route_id="mock-h2",
            source="test",
            parent_smiles="[H][H]",
            reactant_smiles=("[H][H]",),
            product_smiles=("[H][H]",),
            multiplicity=1,
            reactant_multiplicities=(1,),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        mapping = resolve_route_atom_mapping(route)
        with tempfile.TemporaryDirectory() as temp:
            dissociated = write_xyz(
                Path(temp) / "dissociated.xyz", ["H", "H"],
                np.asarray([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
                "mock wrong minimum",
            )
            optimize.return_value = {
                "optimized_xyz": dissociated,
                "optimization": {"out": Path(temp) / "opt.out"},
                "frequency": {"out": Path(temp) / "freq.out"},
                "frequency_check": {"IsMinimum": True},
            }
            records = _stationary_fragment_records(
                route, Path(temp) / "stationary", mapping,
                execute=True, ncores=1, method_keywords=None, timeout_seconds=1.0,
            )
        self.assertFalse(records["valid"])
        self.assertEqual(records["status"], "stationary_point_validation_failed")
        self.assertTrue(all(
            record["status"] == "endpoint_connectivity_validation_failed"
            for side in ("reactants", "products") for record in records[side]
        ))

    def test_noncommon_declared_spin_surface_is_rejected(self):
        route = RouteSpec(
            route_id="wrong-spin",
            source="test",
            parent_smiles="[CH3]",
            reactant_smiles=("[CH3]", "[OH]"),
            product_smiles=("CO",),
            multiplicity=3,
            reactant_multiplicities=(2, 2),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        result = resolve_route_atom_mapping(route)
        self.assertFalse(result["valid"])
        self.assertEqual(result["status"], "declared_multiplicity_not_common")

    def test_profile_classifier_separates_barriered_and_barrierless(self):
        barriered = classify_path_energy_profile(
            [-100.0, -99.99, -99.96, -99.99, -100.01],
            reactant_fragment_count=2,
            product_fragment_count=1,
        )
        capture = classify_path_energy_profile(
            [-100.0, -100.01, -100.02, -100.03, -100.04],
            reactant_fragment_count=2,
            product_fragment_count=1,
        )
        self.assertEqual(barriered["status"], "barriered_ts_candidate")
        self.assertEqual(capture["status"], "barrierless_capture_candidate")

    def test_equal_fragment_no_saddle_profile_retains_reaction_direction(self):
        uphill = classify_path_energy_profile(
            [-336.942205, -336.950666, -336.944031, -336.929571, -336.903440],
            reactant_fragment_count=2,
            product_fragment_count=2,
        )
        self.assertEqual(uphill["status"], "surface_unresolved")
        self.assertEqual(uphill["downhill_direction"], "reverse")
        self.assertFalse(uphill["forward_target_loss_barrierless"])
        self.assertAlmostEqual(uphill["forward_energetic_threshold_hartree"], 0.038765)
        self.assertEqual(uphill["reason"], "forward_path_uphill_without_interior_ts_candidate")

    def test_orca_converged_no_barrier_output_is_semantic_completion(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder / "neb-ts.inp"
            out = folder / "neb-ts.out"
            inp.write_text("! B3LYP def2-SVP NEB-TS\n")

            def failed_orca(input_path, **kwargs):
                out.write_text(
                    "*** THE NEB OPTIMIZATION HAS CONVERGED ***\n"
                    "The elastic band has converged successfully to a MEP. "
                    "However, climbing image was never activated!\n"
                    "No barrier was found. Skipping NEB-TS run here.\n"
                    "ORCA finished with error return - aborting the run\n"
                )
                raise RuntimeError("ORCA failed")

            with patch("storca.reaction_path.run_orca", side_effect=failed_orca):
                result = _run_or_resume(inp, timeout_seconds=10.0)

            self.assertEqual(result["status"], "completed_no_interior_barrier")
            self.assertTrue(result["orca_neb_outcome"]["valid"])
            self.assertFalse(result["orca_neb_outcome"]["climbing_image_activated"])
            self.assertTrue(inp.with_suffix(".contract.json").is_file())

            resumed = _run_or_resume(inp, timeout_seconds=10.0)
            self.assertTrue(resumed["resumed"])
            self.assertEqual(resumed["status"], "completed_no_interior_barrier")

    def test_orca_no_barrier_parser_requires_convergence_and_explicit_message(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "neb-ts.out"
            output.write_text("No barrier was found. Skipping NEB-TS run here.\n")
            parsed = _orca_neb_no_barrier_outcome(output)
        self.assertFalse(parsed["valid"])
        self.assertFalse(parsed["neb_converged"])

    def test_prepare_only_returns_canonical_fail_closed_contract(self):
        route = RouteSpec(
            route_id="prepare",
            source="test",
            parent_smiles="[CH3]",
            reactant_smiles=("[CH3]", "[OH]"),
            product_smiles=("CO",),
            multiplicity=1,
            reactant_multiplicities=(2, 2),
            product_multiplicities=(1,),
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(route, Path(temp), prepare_only=True)
            verification = result["route_verification"]
            self.assertEqual(verification["status"], "prepared")
            self.assertEqual(verification["rate_model"]["status"], "rate_model_unavailable")
            self.assertIsNone(verification["trajectory"])
            self.assertEqual(len(result["orientations"]), 3)
            self.assertEqual(
                result["endpoint_complexes"]["fragment_geometry_source"],
                "supplied_atom_ordered_fragment_xyz",
            )
            self.assertEqual(verification["endpoint_mode"], "separated_fragment_channels")

    @patch("storca.reaction_path.build_endpoint_complex_seeds")
    @patch("storca.reaction_path.condition_specific_forward_loss_bound")
    @patch("storca.reaction_path.assemble_route_thermochemistry")
    @patch("storca.reaction_path._stationary_fragment_records")
    def test_thermochemical_loss_bound_stops_before_neb(
        self, stationary, assemble_thermo, loss_bound, build_seeds,
    ):
        stationary.return_value = {
            "status": "validated", "valid": True, "reactants": [], "products": [],
        }
        assemble_thermo.return_value = {"status": "validated", "valid": True}
        loss_bound.return_value = {
            "status": "forward_loss_below_retention_threshold_upper_bound",
            "applicable": True,
            "target_loss_fraction_upper_bound": 0.001,
        }
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(
                self._association_route(), Path(temp), condition_contract={
                    "temperature_K": 298.15, "pressure_bar": 1.0,
                    "target_duration_seconds": 1000.0, "retention_fraction": 0.95,
                },
            )
        verification = result["route_verification"]
        self.assertEqual(verification["status"], "thermochemistry_bounded_below_loss_threshold")
        self.assertEqual(verification["convergence"]["status"], "not_required_after_conservative_loss_bound")
        build_seeds.assert_not_called()

    def test_public_rmg_route_normalizer_retains_fragment_states(self):
        route = route_spec_from_rmg_route("[CH3]", {
            "route_id": "rmg-1",
            "reaction_equation": "methyl(1) + hydroxyl(2) <=> methanol(3)",
            "multiplicity": 1,
            "resolved_endpoints": {
                "reactants": [
                    {"label": "methyl(1)", "smiles": "[CH3]", "charge": 0, "multiplicity": 2},
                    {"label": "hydroxyl(2)", "smiles": "[OH]", "charge": 0, "multiplicity": 2},
                ],
                "products": [{"label": "methanol(3)", "smiles": "CO", "charge": 0, "multiplicity": 1}],
            },
        })
        self.assertEqual(route.protocol, "generic_endpoint_path")
        self.assertEqual(route.reactant_charges, (0, 0))
        self.assertEqual(route.reactant_multiplicities, (2, 2))
        self.assertEqual(route.product_multiplicities, (1,))
        self.assertEqual(route.reactant_labels, ("methyl(1)", "hydroxyl(2)"))
        self.assertEqual(route.reaction_equation, "methyl(1) + hydroxyl(2) <=> methanol(3)")

    def test_stoichiometric_duplicates_reuse_original_stationary_label(self):
        route = RouteSpec(
            route_id="homolysis",
            source="test",
            parent_smiles="OO",
            reactant_smiles=("OO",),
            product_smiles=("[OH]", "[OH]"),
            reactant_labels=("peroxide",),
            product_labels=("hydroxyl", "hydroxyl"),
            multiplicity=1,
            reactant_multiplicities=(1,),
            product_multiplicities=(2, 2),
            atom_mapping=((0, 0), (1, 1)),
            protocol="generic_endpoint_path",
        )
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(route, Path(temp), prepare_only=True)
            products = result["route_verification"]["stationary_points"]["products"]
            self.assertEqual([record["label"] for record in products], ["hydroxyl", "hydroxyl"])
            self.assertEqual(products[0]["input_xyz"], products[1]["input_xyz"])
            self.assertTrue(products[1]["reused_stationary_point"])

    @patch("storca.reaction_path._stationary_fragment_xyz", return_value=None)
    @patch("storca.reaction_path._stationary_fragment_records")
    @patch("storca.reaction_path._run_generic_orientation")
    def test_reproducible_barrierless_path_stops_after_two_orientations(
        self, run_orientation, stationary, stationary_xyz,
    ):
        def result(seed):
            orientation = seed["orientation"]
            return {
                "orientation": orientation,
                "status": "classified_without_transition_state",
                "path_classification": "barrierless_capture_candidate",
                "endpoint_validation": {"valid": True},
                "trajectory": f"orientation-{orientation}.xyz",
                "path_evidence": {
                    "energies_hartree": [-100.0, -100.01, -100.02, -100.03, -100.04],
                },
            }

        run_orientation.side_effect = lambda route, seed, **kwargs: result(seed)
        stationary.return_value = {"status": "validated", "valid": True, "reactants": [], "products": []}
        with tempfile.TemporaryDirectory() as temp:
            verification = run_generic_reaction_path(self._association_route(), Path(temp))["route_verification"]
        self.assertEqual(run_orientation.call_count, 2)
        self.assertEqual(verification["status"], "barrierless_rate_model_required")
        self.assertEqual(verification["barrierless_reproducibility"]["status"], "reproducible")

    @patch("storca.reaction_path._stationary_fragment_xyz", return_value=None)
    @patch("storca.reaction_path._stationary_fragment_records")
    @patch("storca.reaction_path._run_generic_orientation")
    def test_explicitly_qualified_submerged_candidate_escalation_is_preserved(
        self, run_orientation, stationary, stationary_xyz,
    ):
        run_orientation.return_value = {
            "orientation": 1,
            "status": "dissociation_rate_model_required",
            "path_classification": "submerged_saddle_dissociation_candidate",
            "endpoint_validation": {"valid": True},
            "endpoint_reference": {
                "kind": "bound_minimum",
                "sources": ["validated-bound-opt.out", "validated-bound-freq.out"],
                "energetic_classification_permitted": True,
                "mechanism_claim_allowed": False,
                "rate_claim_allowed": False,
            },
            "path_evidence": {
                "kind": "orca_submerged_climbing_image_dissociation_diagnostic",
                "rate_claim_allowed": False,
            },
            "method_selection": {
                "ordinary_harmonic_tst_applicable": False,
                "rate_claim_allowed": False,
            },
        }
        stationary.return_value = {
            "status": "validated", "valid": True, "reactants": [], "products": [],
        }

        with tempfile.TemporaryDirectory() as temp:
            verification = run_generic_reaction_path(
                self._association_route(), Path(temp),
            )["route_verification"]

        self.assertEqual(run_orientation.call_count, 1)
        self.assertEqual(verification["status"], "dissociation_rate_model_required")
        self.assertEqual(
            verification["path_classification"],
            "submerged_saddle_dissociation_candidate",
        )
        self.assertEqual(verification["convergence"]["status"], "surface_resolution_incomplete")

    @patch("storca.reaction_path._stationary_fragment_xyz", return_value=None)
    @patch("storca.reaction_path._stationary_fragment_records")
    @patch("storca.reaction_path._run_generic_orientation")
    def test_validated_intermediate_status_is_preserved_at_route_level(
        self, run_orientation, stationary, stationary_xyz,
    ):
        run_orientation.return_value = {
            "orientation": 1,
            "status": "intermediate_detected_requires_segmented_verification",
            "path_classification": "multistep_intermediate_detected",
            "endpoint_validation": {"valid": True},
            "path_evidence": {"kind": "orca_neb_intermediate_path", "rate_claim_allowed": False},
            "intermediate_evidence": {
                "status": "validated_intermediate_detected",
                "segmented_verification": {"status": "computed"},
            },
        }
        stationary.return_value = {"status": "validated", "valid": True, "reactants": [], "products": []}

        with tempfile.TemporaryDirectory() as temp:
            verification = run_generic_reaction_path(self._association_route(), Path(temp))["route_verification"]

        # A validated stationary intermediate changes the scientific task to
        # segment verification; equivalent full-endpoint orientations would
        # spend compute without resolving that new task.
        self.assertEqual(run_orientation.call_count, 1)
        self.assertEqual(verification["status"], "intermediate_detected_requires_segmented_verification")
        self.assertEqual(verification["path_classification"], "multistep_intermediate_detected")
        self.assertEqual(verification["selected_orientation"], 1)
        self.assertEqual(verification["intermediate_evidence"]["segmented_verification"]["status"], "computed")
        self.assertEqual(verification["convergence"]["status"], "segmented_path_verification_incomplete")

    @patch("storca.reaction_path._run_bond_change_guided_ts_recovery")
    @patch("storca.reaction_path._refine_reached_irc_stationary_points", return_value={
        "stationary_points": {"status": "validated", "valid": True, "reactants": [], "products": []}})
    @patch("storca.reaction_path._stationary_fragment_xyz", return_value=None)
    @patch("storca.reaction_path._stationary_fragment_records")
    @patch("storca.reaction_path._run_generic_orientation")
    def test_two_unvalidated_intermediate_paths_trigger_one_guided_recovery(
        self, run_orientation, stationary, stationary_xyz, refine, guided,
    ):
        def unresolved(route, seed, **kwargs):
            orientation = int(seed["orientation"])
            return {
                "orientation": orientation,
                "status": "surface_unresolved",
                "path_classification": "surface_unresolved",
                "endpoint_validation": {"valid": True},
                "intermediate_evidence": {
                    "status": "intermediate_detected_unvalidated",
                    "candidates": [{"status": "intermediate_optimization_unconverged"}],
                },
                "surface_diagnosis": {
                    "schema_version": 1,
                    "kind": "storca_surface_diagnosis",
                    "status": "intermediate_candidate",
                    "reasons": ["orca_reported_possible_intermediate"],
                    "evidence": {"orientation": orientation},
                    "corroborated": True,
                    "recommended_strategy": "stationary_point_recovery",
                    "claims_permitted": [],
                    "claims_prohibited": ["mechanism", "rate", "lifetime"],
                },
            }

        run_orientation.side_effect = unresolved
        guided.return_value = {
            "status": "verified_transition_state_path",
            "path_classification": "verified_barriered_path",
            "transition_state": {"frequency_check": {"NumImaginary": 1}},
            "irc": {"endpoint_validation": {"valid": True}},
            "trajectory": "guided-irc.xyz",
            "path_evidence": {"kind": "orca_ts_frequency_irc"},
        }
        stationary.return_value = {
            "status": "validated", "valid": True, "reactants": [], "products": [],
        }

        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(self._association_route(), Path(temp))

        self.assertEqual(run_orientation.call_count, 2)
        guided.assert_called_once()
        refine.assert_called_once()
        self.assertTrue(result["surface_resolution_guided_recovery"]["equivalent_full_neb_stopped"])
        self.assertEqual(result["route_verification"]["status"], "verified_transition_state_path")

    @patch("storca.reaction_path._refine_reached_irc_stationary_points", return_value={
        "stationary_points": {"status": "validated", "valid": True, "reactants": [], "products": []}})
    @patch("storca.reaction_path._stationary_fragment_xyz", return_value=None)
    @patch("storca.reaction_path._stationary_fragment_records")
    @patch("storca.reaction_path._run_generic_orientation")
    def test_verified_ts_path_stops_after_first_orientation(
        self, run_orientation, stationary, stationary_xyz, refine,
    ):
        run_orientation.return_value = {
            "orientation": 1,
            "status": "verified_transition_state_path",
            "path_classification": "verified_barriered_path",
            "endpoint_validation": {"valid": True},
            "transition_state": {"frequency_output": "ts-frequency.out"},
            "irc": {"endpoint_validation": {"valid": True}},
            "trajectory": "validated-trajectory.xyz",
            "path_evidence": {"kind": "orca_ts_frequency_irc"},
        }
        stationary.return_value = {"status": "validated", "valid": True, "reactants": [], "products": []}
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(self._association_route(), Path(temp))
            verification = result["route_verification"]
            attempt_evidence = json.loads(Path(
                result["orientations"][0]["surface_resolution_attempt_evidence"]["path"]
            ).read_text())
        self.assertEqual(run_orientation.call_count, 1)
        refine.assert_called_once()
        self.assertEqual(verification["status"], "verified_transition_state_path")
        self.assertEqual(verification["trajectory"], "validated-trajectory.xyz")
        self.assertEqual(attempt_evidence["strategy"], "standard_neb_ts")
        self.assertEqual(
            attempt_evidence["scientific_context"]["environment_applicability"]["canonical_model"],
            "homogeneous_gas",
        )
        self.assertIn("mapping", attempt_evidence["scientific_context"])

    @patch("storca.reaction_path._run_or_resume", side_effect=RuntimeError("mock NEB stop"))
    @patch("storca.reaction_path._optimize_endpoint")
    def test_separated_nno_o2_channel_bypasses_encounter_complex_optimization(
        self, optimize_endpoint, run_or_resume,
    ):
        route = RouteSpec(
            route_id="nno-o2-initiation",
            source="test",
            parent_smiles="NNO",
            reactant_smiles=("[O][O]", "NNO"),
            product_smiles=("[O]O", "NN[O]"),
            multiplicity=3,
            reactant_multiplicities=(3, 1),
            product_multiplicities=(2, 2),
            protocol="generic_endpoint_path",
        )
        mapping = resolve_route_atom_mapping(route)
        with tempfile.TemporaryDirectory() as temp:
            seeds = build_endpoint_complex_seeds(
                route, Path(temp) / "seeds", orientations=1, mapping_result=mapping,
            )
            result = _run_generic_orientation(
                route, seeds["orientations"][0], mapping_result=mapping,
                ncores=1, method_keywords=None, timeout_seconds=1.0,
                neb_images=4, barrierless_classifier=None,
                separated_fragments=True,
                stationary_points={"status": "validated", "valid": True},
            )
        optimize_endpoint.assert_not_called()
        run_or_resume.assert_called_once()
        self.assertEqual(
            result["endpoint_validation"]["mode"],
            "assembled_from_validated_separated_fragment_minima",
        )
        self.assertFalse(result["endpoint_validation"]["full_complex_minimum_required"])
        self.assertEqual(result["status"], "surface_unresolved")

    @patch("storca.reaction_path.build_endpoint_complex_seeds")
    @patch("storca.reaction_path._stationary_fragment_records")
    def test_invalid_separated_species_blocks_before_encounter_assembly(self, stationary, build_seeds):
        stationary.return_value = {
            "status": "stationary_point_validation_failed",
            "valid": False,
            "reactants": [{"status": "failed", "failure_reason": "mock failure"}],
            "products": [],
        }
        with tempfile.TemporaryDirectory() as temp:
            result = run_generic_reaction_path(self._association_route(), Path(temp))
        build_seeds.assert_not_called()
        verification = result["route_verification"]
        self.assertEqual(verification["status"], "stationary_point_validation_failed")
        self.assertEqual(verification["convergence"]["status"], "blocked_before_neb")


if __name__ == "__main__":
    unittest.main()
