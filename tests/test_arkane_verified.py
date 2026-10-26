import tempfile
import unittest
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

from storca.arkane_runner import (
    VerifiedArkaneRouteSpec,
    _write_verified_kinetics_library,
    _prepare_arkane_species_dictionary,
    _evaluate_kinetics_expression,
    _extract_arkane_rate_model,
    create_verified_arkane_input,
    validate_verified_arkane_library,
    validate_arkane_input_syntax,
)
from storca.generated_kinetics import validate_generated_library, validate_rate_replacement_evidence
from arkane_rotor_fixtures import conditional_rate_evidence_fixture, rotor_source_artifacts, rotor_validation_fixture


def _artifact(folder: Path, name: str) -> Path:
    output = folder / f"{name}.out"
    output.write_text(name)
    output.with_suffix(".hess").write_text("hessian")
    return output


def _fake_validation(path: Path, **kwargs) -> dict:
    path = Path(path)
    return {
        "output": str(path.resolve()), "hessian": str(path.with_suffix(".hess").resolve()),
        "frequency_check": {"NumImaginary": int(kwargs["transition_state"])},
        "formula": kwargs["expected_formula"], "charge": kwargs["expected_charge"],
        "multiplicity": kwargs["expected_multiplicity"], "checksums": {},
    }


class VerifiedArkaneTests(unittest.TestCase):
    def test_transition_state_must_share_both_endpoint_spin_surfaces(self):
        # H(doublet) + H(doublet) can couple to singlet or triplet, but the
        # declared H2 singlet endpoint cannot connect through a triplet TS on
        # the same adiabatic surface. No ORCA artifacts or jobs are needed to
        # detect this inconsistency at rate admission.
        for selected in (2, 3):
            with self.subTest(transition_state_multiplicity=selected), \
                    self.assertRaisesRegex(ValueError, "common endpoint spin surface"):
                VerifiedArkaneRouteSpec(
                    label="invalid-spin-recombination",
                    reactant_labels=("hydrogen", "hydrogen"),
                    reactant_smiles=("[H]", "[H]"),
                    reactant_orca_outputs=(Path("hydrogen.out"), Path("hydrogen.out")),
                    reactant_multiplicities=(2, 2),
                    product_labels=("hydrogen2",), product_smiles=("[H][H]",),
                    product_orca_outputs=(Path("hydrogen2.out"),),
                    product_multiplicities=(1,),
                    transition_state_orca_output=Path("ts.out"),
                    transition_state_multiplicity=selected, model_chemistry="r2SCAN-3c",
                    collision_limit_m3_mol_s=1.0, collision_limit_source="fixture ceiling",
                )

    def test_missing_native_dictionary_export_preserves_and_checks_evidence(self):
        from types import SimpleNamespace
        import subprocess
        spec = SimpleNamespace(reactant_labels=("a",), product_labels=("b",),
                               reactant_smiles=("[H]",), product_smiles=("[H]",),
                               reactant_multiplicities=(2,), product_multiplicities=(2,))
        report = {"dictionary": "fixture dictionary\n", "labels": ["a", "b"], "rmg_version": "fixture"}
        response = subprocess.CompletedProcess([], 0,
            stdout="STORCA_ARKANE_DICTIONARY_JSON=" + json.dumps(report), stderr="")
        with tempfile.TemporaryDirectory() as temp, \
             patch("storca.arkane_runner.run_bounded_kinetics_command", return_value=response) as run:
            root = Path(temp)
            path = _prepare_arkane_species_dictionary(spec, root, rmg_env="rmg4_env")
            self.assertEqual(path.read_text(), report["dictionary"])
            self.assertIn("rmg4_env", run.call_args.args[0])
            self.assertEqual(_prepare_arkane_species_dictionary(spec, root, rmg_env="rmg4_env"), path)
            path.write_text("mutated\n")
            with self.assertRaisesRegex(ValueError, "differs"):
                _prepare_arkane_species_dictionary(spec, root, rmg_env="rmg4_env")
            self.assertEqual(path.read_text(), "mutated\n")

    @staticmethod
    def _explicit_pressure_model(root: Path, *, well_label: str = "methane") -> dict:
        evidence = root / "collision-model-reference.txt"
        evidence.write_text("fixture calibrated collision model\n")
        common = {"source": "fixture calibration", "applicability": "fixture methane in nitrogen"}
        return {
            "schema_version": 1,
            "kind": "storca_explicit_gas_pressure_dependence_model",
            "status": "applicable",
            "temperature_range_K": [250.0, 1000.0],
            "pressure_range_bar": [1.0e-6, 10.0],
            "energy_transfer_model": {
                "model": "SingleExponentialDown", "well_label": well_label,
                "alpha0_kJ_mol": 3.5, "T0_K": 300.0, "n": 0.85, **common,
            },
            "transport_models": {
                well_label: {
                    "model": "LennardJones", "sigma_angstrom": 3.8,
                    "epsilon_K": 150.0, "shape_index": 2, **common,
                },
                "nitrogen": {
                    "model": "LennardJones", "sigma_angstrom": 3.7,
                    "epsilon_K": 100.0, "shape_index": 1, **common,
                },
            },
            "evidence_artifacts": [{
                "path": str(evidence),
                "sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
            }],
        }

    def test_verified_repair_libraries_have_route_unique_rmg_labels(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workdirs = [root / "first", root / "second"]
            specs = []
            for label, workdir in zip(("route-a", "route-b"), workdirs):
                native = workdir / "RMG_libraries" / "kinetics"
                native.mkdir(parents=True)
                (native / "dictionary.txt").write_text("placeholder\n1 H u1 p0 c0\n")
                specs.append(VerifiedArkaneRouteSpec(
                    label=label,
                    reactant_labels=("a",), reactant_smiles=("[H]",),
                    reactant_orca_outputs=(Path("a"),), reactant_multiplicities=(2,),
                    product_labels=("b",), product_smiles=("[H]",),
                    product_orca_outputs=(Path("b"),), product_multiplicities=(2,),
                    transition_state_orca_output=Path("ts"), transition_state_multiplicity=2,
                    model_chemistry="r2SCAN-3c",
                ))
            libraries = [
                _write_verified_kinetics_library(
                    spec, workdir,
                    kinetics_expression="Arrhenius(A=(1,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'))",
                )
                for spec, workdir in zip(specs, workdirs)
            ]
        self.assertEqual(len({path.name for path in libraries}), 2)
        self.assertEqual([path.name for path in libraries], [
            "storca-verified-route-a", "storca-verified-route-b",
        ])

    def test_pdep_dissociation_and_association_are_oriented_around_the_well(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            methane = _artifact(folder, "methane")
            methyl = _artifact(folder, "methyl")
            hydrogen = _artifact(folder, "hydrogen")
            ts = _artifact(folder, "ts")
            common = dict(
                label="route", model_chemistry="r2SCAN-3c", transition_state_orca_output=ts,
                transition_state_multiplicity=1,
            )
            dissociation = VerifiedArkaneRouteSpec(
                reactant_labels=("methane",), reactant_smiles=("C",), reactant_orca_outputs=(methane,),
                reactant_multiplicities=(1,), product_labels=("methyl", "hydrogen"),
                product_smiles=("[CH3]", "[H]"), product_orca_outputs=(methyl, hydrogen),
                product_multiplicities=(2, 2), **common,
            )
            association = VerifiedArkaneRouteSpec(
                reactant_labels=("methyl", "hydrogen"), reactant_smiles=("[CH3]", "[H]"),
                reactant_orca_outputs=(methyl, hydrogen), reactant_multiplicities=(2, 2),
                product_labels=("methane",), product_smiles=("C",), product_orca_outputs=(methane,),
                product_multiplicities=(1,), collision_limit_m3_mol_s=1.0e6,
                collision_limit_source="test hard-sphere bound", **common,
            )
            with patch("storca.arkane_runner._validate_orca_frequency_artifact", side_effect=_fake_validation):
                dissociation_input = create_verified_arkane_input(dissociation, folder / "dissociation").read_text()
                association_input = create_verified_arkane_input(association, folder / "association").read_text()
            # Execute the generated statmech declarations with a fake energy
            # reader, testing the actual adapter expressions (not a duplicate
            # reference-shift implementation). No electronic jobs run here.
            raw = {"methane": -100000000.0, "methyl": -80000000.0,
                   "hydrogen": -19900000.0, "ts": -99800000.0}
            class EnergyLog:
                def __init__(self, path):
                    self.name = Path(path).read_text()
                def load_energy(self):
                    return raw[self.name]
            for branch in ("dissociation", "association"):
                provenance = json.loads((folder / branch / "provenance.json").read_text())
                energies = {}
                for record in provenance["orca_artifacts"]:
                    context = {"Log": EnergyLog}
                    exec(Path(record["statmech_file"]).read_text(), context)
                    energy, units, kind = context["energy"]
                    self.assertEqual((units, kind), ("J/mol", "e_electronic"))
                    self.assertFalse(record["energy_reference"]["zero_point_energy_included"])
                    energies[record["label"]] = energy
                self.assertAlmostEqual(energies["methyl"] + energies["hydrogen"]
                                       - energies["methane"], 100000.0, places=6)
                self.assertAlmostEqual(energies["__storca_ts__"] - energies["methane"],
                                       200000.0, places=6)
            self.assertIn("isomers=['methane']", dissociation_input)
            self.assertIn("products=[['methyl', 'hydrogen']]", dissociation_input)
            self.assertIn("isomers=['methane']", association_input)
            self.assertIn("reactants=[['methyl', 'hydrogen']]", association_input)
            self.assertIn("products=[]", association_input)

    def test_stoichiometric_duplicate_is_defined_once_but_retained_in_reaction(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            peroxide = _artifact(folder, "peroxide")
            hydroxyl = _artifact(folder, "hydroxyl")
            ts = _artifact(folder, "ts")
            spec = VerifiedArkaneRouteSpec(
                label="homolysis", reactant_labels=("peroxide",), reactant_smiles=("OO",),
                reactant_orca_outputs=(peroxide,), reactant_multiplicities=(1,),
                product_labels=("hydroxyl", "hydroxyl"), product_smiles=("[OH]", "[OH]"),
                product_orca_outputs=(hydroxyl, hydroxyl), product_multiplicities=(2, 2),
                transition_state_orca_output=ts, transition_state_multiplicity=1,
                model_chemistry="r2SCAN-3c",
            )
            with patch("storca.arkane_runner._validate_orca_frequency_artifact", side_effect=_fake_validation):
                text = create_verified_arkane_input(spec, folder / "arkane").read_text()
            self.assertEqual(text.count("species('hydroxyl'"), 1)
            self.assertIn("products=['hydroxyl', 'hydroxyl']", text)

    def test_bimolecular_route_fails_without_collision_bound(self):
        with self.assertRaisesRegex(ValueError, "collision limit"):
            VerifiedArkaneRouteSpec(
                label="exchange", reactant_labels=("a", "b"), reactant_smiles=("[H]", "[OH]"),
                reactant_orca_outputs=(Path("a"), Path("b")), reactant_multiplicities=(2, 2),
                product_labels=("c", "d"), product_smiles=("[H]", "[OH]"),
                product_orca_outputs=(Path("c"), Path("d")), product_multiplicities=(2, 2),
                transition_state_orca_output=Path("ts"), transition_state_multiplicity=1,
                model_chemistry="r2SCAN-3c",
            )

    def test_pdep_extractor_requires_original_direction(self):
        spec = VerifiedArkaneRouteSpec(
            label="association", reactant_labels=("methyl", "hydrogen"),
            reactant_smiles=("[CH3]", "[H]"), reactant_orca_outputs=(Path("a"), Path("b")),
            reactant_multiplicities=(2, 2), product_labels=("methane",), product_smiles=("C",),
            product_orca_outputs=(Path("c"),), product_multiplicities=(1,),
            transition_state_orca_output=Path("ts"), transition_state_multiplicity=1,
            model_chemistry="r2SCAN-3c", collision_limit_m3_mol_s=1.0e6,
            collision_limit_source="test hard-sphere bound",
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output.py"
            output.write_text(
                "pdepreaction(reactants=['methane'], products=['methyl', 'hydrogen'], "
                "kinetics=Chebyshev(coeffs=[[1]], kunits='s^-1', Tmin=(250,'K'), Tmax=(1000,'K'), "
                "Pmin=(1e-6,'bar'), Pmax=(10,'bar')))\n"
            )
            with self.assertRaisesRegex(ValueError, "0 matching"):
                _extract_arkane_rate_model(spec, output)

    def test_rate_expression_evaluates_in_both_installed_rmg_versions(self):
        expression = "Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'))"
        for environment in ("rmg_env", "rmg4_env"):
            with self.subTest(environment=environment):
                report = _evaluate_kinetics_expression(
                    expression, temperature_K=298.15, pressure_bar=1.0,
                    pressure_dependent=False, rmg_env=environment,
                )
                self.assertAlmostEqual(report["value_si"], 1.0e10)
                self.assertIn(report["rmg_version"], {"3.3.0", "4.0.0"})

    def test_generated_pdep_input_parses_in_both_installed_arkane_versions(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            methane = _artifact(folder, "methane")
            methyl = _artifact(folder, "methyl")
            hydrogen = _artifact(folder, "hydrogen")
            ts = _artifact(folder, "ts")
            spec = VerifiedArkaneRouteSpec(
                label="dissociation", reactant_labels=("methane",), reactant_smiles=("C",),
                reactant_orca_outputs=(methane,), reactant_multiplicities=(1,),
                product_labels=("methyl", "hydrogen"), product_smiles=("[CH3]", "[H]"),
                product_orca_outputs=(methyl, hydrogen), product_multiplicities=(2, 2),
                transition_state_orca_output=ts, transition_state_multiplicity=1,
                model_chemistry="r2SCAN-3c",
            )
            with patch("storca.arkane_runner._validate_orca_frequency_artifact", side_effect=_fake_validation):
                input_file = create_verified_arkane_input(spec, folder / "arkane")
            for environment in ("rmg_env", "rmg4_env"):
                with self.subTest(environment=environment):
                    report = validate_arkane_input_syntax(input_file, rmg_env=environment)
                    self.assertEqual(report["reaction_labels"], ["dissociation"])
                    self.assertIn("PressureDependenceJob", report["job_types"])

    def test_explicit_collisional_model_is_written_and_parses(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            methane = _artifact(folder, "methane")
            methyl = _artifact(folder, "methyl")
            hydrogen = _artifact(folder, "hydrogen")
            ts = _artifact(folder, "ts")
            spec = VerifiedArkaneRouteSpec(
                label="explicit-dissociation",
                reactant_labels=("methane",), reactant_smiles=("C",),
                reactant_orca_outputs=(methane,), reactant_multiplicities=(1,),
                product_labels=("methyl", "hydrogen"),
                product_smiles=("[CH3]", "[H]"),
                product_orca_outputs=(methyl, hydrogen), product_multiplicities=(2, 2),
                transition_state_orca_output=ts, transition_state_multiplicity=1,
                model_chemistry="r2SCAN-3c",
                explicit_pressure_dependence_model=self._explicit_pressure_model(folder),
            )
            with patch(
                "storca.arkane_runner._validate_orca_frequency_artifact",
                side_effect=_fake_validation,
            ):
                input_file = create_verified_arkane_input(spec, folder / "arkane")
            text = input_file.read_text()
            self.assertIn("energyTransferModel=SingleExponentialDown", text)
            self.assertIn("collisionModel=TransportData", text)
            provenance = json.loads((input_file.parent / "provenance.json").read_text())
            self.assertEqual(
                provenance["explicit_pressure_dependence_model"]["status"], "applicable",
            )
            self.assertTrue(
                provenance["explicit_pressure_dependence_model"]["model_sha256"].startswith(
                    "sha256:"
                )
            )
            report = validate_arkane_input_syntax(input_file, rmg_env="rmg4_env")
            self.assertIn("PressureDependenceJob", report["job_types"])

    def test_explicit_collisional_model_rejects_mutated_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            model = self._explicit_pressure_model(folder)
            Path(model["evidence_artifacts"][0]["path"]).write_text("mutated\n")
            spec = VerifiedArkaneRouteSpec(
                label="mutated", reactant_labels=("methane",), reactant_smiles=("C",),
                reactant_orca_outputs=(Path("a"),), reactant_multiplicities=(1,),
                product_labels=("methyl", "hydrogen"), product_smiles=("[CH3]", "[H]"),
                product_orca_outputs=(Path("b"), Path("c")), product_multiplicities=(2, 2),
                transition_state_orca_output=Path("ts"), transition_state_multiplicity=1,
                model_chemistry="r2SCAN-3c", explicit_pressure_dependence_model=model,
            )
            with self.assertRaisesRegex(ValueError, "missing or mutated"):
                create_verified_arkane_input(spec, folder / "arkane")

    def test_generated_manifest_satisfies_replacement_evidence_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            native = folder / "RMG_libraries" / "kinetics"
            native.mkdir(parents=True)
            (native / "dictionary.txt").write_text("placeholder\n1 H u1 p0 c0\n")
            (native / "reactions.py").write_text("name = 'native'\n")
            output = folder / "output.py"
            rate_output = (
                "# k (TST)\n# k (TST+T)\n"
                "kinetics(label='exchange', kinetics=Arrhenius(A=(1,'m^3/(mol*s)'), n=0, "
                "Ea=(0,'kJ/mol'), T0=(1,'K')))\n"
            )
            spec = VerifiedArkaneRouteSpec(
                label="exchange", reactant_labels=("h_left", "oh_left"),
                reactant_smiles=("[H]", "[OH]"), reactant_orca_outputs=(Path("a"), Path("b")),
                reactant_multiplicities=(2, 2), product_labels=("oh_right", "h_right"),
                product_smiles=("[OH]", "[H]"), product_orca_outputs=(Path("c"), Path("d")),
                product_multiplicities=(2, 2), transition_state_orca_output=Path("ts"),
                transition_state_multiplicity=1, model_chemistry="r2SCAN-3c",
                reaction_equation="h_left(7) + oh_left(8) <=> oh_right(9) + h_right(10)",
                collision_limit_m3_mol_s=2.0, collision_limit_source="test collision bound",
            )
            evaluated = {"value_si": 1.0, "model_class": "Arrhenius",
                         "canonical_expression": "Arrhenius(A=(1,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'))"}
            library_check = {
                "value_si": 1.0, "model_class": "Arrhenius", "pressure_dependent": False,
                "reactants": list(spec.reactant_labels), "products": list(spec.product_labels),
            }
            statmech = rotor_validation_fixture(folder, route_id=spec.label, rate_output=rate_output)
            with patch("storca.arkane_runner.extract_retained_arkane_rate_evidence",
                          return_value=conditional_rate_evidence_fixture(statmech)), \
                    patch("storca.arkane_runner._evaluate_kinetics_expression", return_value=evaluated), \
                    patch("storca.arkane_runner._validate_library_with_rmg", return_value=library_check):
                manifest = validate_verified_arkane_library(spec, folder, output)
            self.assertEqual(manifest["arkane_statmech_validation"], statmech)
            self.assertEqual(manifest["arkane_rate_source_artifacts"], rotor_source_artifacts(statmech))
            self.assertEqual(validate_rate_replacement_evidence(
                manifest, temperature_K=298.15, pressure_bar=1.0,
            )["status"], "passed")
            self.assertEqual(
                manifest["reaction_equation"],
                "h_left(7) + oh_left(8) <=> oh_right(9) + h_right(10)",
            )
            self.assertEqual(
                manifest["library_reaction_equation"],
                "h_left + oh_left <=> oh_right + h_right",
            )
            self.assertEqual(
                manifest["rate_uncertainty"]["status"], "calibration_required",
            )
            self.assertEqual(validate_generated_library(
                Path(manifest["library"]), temperature_K=298.15, pressure_bar=1.0,
            )["verification_status"], "verified_arkane_tst")


if __name__ == "__main__":
    unittest.main()
