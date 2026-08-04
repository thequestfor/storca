import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from storca.arkane_runner import (
    VerifiedArkaneRouteSpec,
    _write_verified_kinetics_library,
    _evaluate_kinetics_expression,
    _extract_arkane_rate_model,
    create_verified_arkane_input,
    validate_verified_arkane_library,
    validate_arkane_input_syntax,
)
from storca.generated_kinetics import validate_generated_library, validate_rate_replacement_evidence


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
                transition_state_multiplicity=2,
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
            transition_state_orca_output=Path("ts"), transition_state_multiplicity=2,
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
                transition_state_orca_output=ts, transition_state_multiplicity=2,
                model_chemistry="r2SCAN-3c",
            )
            with patch("storca.arkane_runner._validate_orca_frequency_artifact", side_effect=_fake_validation):
                input_file = create_verified_arkane_input(spec, folder / "arkane")
            for environment in ("rmg_env", "rmg4_env"):
                with self.subTest(environment=environment):
                    report = validate_arkane_input_syntax(input_file, rmg_env=environment)
                    self.assertEqual(report["reaction_labels"], ["dissociation"])
                    self.assertIn("PressureDependenceJob", report["job_types"])

    def test_generated_manifest_satisfies_replacement_evidence_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            native = folder / "RMG_libraries" / "kinetics"
            native.mkdir(parents=True)
            (native / "dictionary.txt").write_text("placeholder\n1 H u1 p0 c0\n")
            (native / "reactions.py").write_text("name = 'native'\n")
            output = folder / "output.py"
            output.write_text(
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
            with patch("storca.arkane_runner._evaluate_kinetics_expression", return_value=evaluated), \
                    patch("storca.arkane_runner._validate_library_with_rmg", return_value=library_check):
                manifest = validate_verified_arkane_library(spec, folder, output)
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
            self.assertEqual(validate_generated_library(
                Path(manifest["library"]), temperature_K=298.15, pressure_bar=1.0,
            )["verification_status"], "verified_arkane_tst")


if __name__ == "__main__":
    unittest.main()
