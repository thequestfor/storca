import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from storca.local_mode_fd import (
    ANGSTROM_PER_BOHR,
    ATOMIC_MASS_UNIT_IN_ELECTRON_MASSES,
    DEBYE_PER_ATOMIC_DIPOLE,
    HARTREE_VIBRATIONAL_WAVENUMBER_CM_1,
    IR_INTENSITY_CONVERSION_KM_MOL,
    LocalModeHessianValidationConfig,
    OrcaInvocationLedger,
    local_stretch_displacement,
    parse_orca_dipole,
    projected_local_mode_from_differences,
    validate_local_modes_against_harmonic_subspace,
)


class LocalModeFiniteDifferenceTests(unittest.TestCase):
    def test_displacement_changes_only_bond_length_and_preserves_pair_center_of_mass(self):
        symbols = ["O", "H", "C"]
        coordinates = np.asarray([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [2.0, 1.0, 0.0]])
        displaced, detail = local_stretch_displacement(
            symbols, coordinates, 0, 1, 0.01,
        )
        original_center = (
            detail["heavy_mass_u"] * coordinates[0]
            + detail["hydrogen_mass_u"] * coordinates[1]
        ) / (detail["heavy_mass_u"] + detail["hydrogen_mass_u"])
        displaced_center = (
            detail["heavy_mass_u"] * displaced[0]
            + detail["hydrogen_mass_u"] * displaced[1]
        ) / (detail["heavy_mass_u"] + detail["hydrogen_mass_u"])
        self.assertTrue(np.allclose(original_center, displaced_center, atol=1e-14))
        self.assertAlmostEqual(np.linalg.norm(displaced[1] - displaced[0]), 0.97, places=12)
        self.assertTrue(np.array_equal(displaced[2], coordinates[2]))

    def test_dipole_parser_uses_last_orca_total_vector(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "gradient.out"
            output.write_text(
                "Total Dipole Moment : 1.0 2.0 3.0\n"
                "Total Dipole Moment    : -4.0D-1 5.0E-1 6.0e-1\n"
            )
            dipole = parse_orca_dipole(output)
        self.assertTrue(np.allclose(dipole, [-0.4, 0.5, 0.6]))

    def test_projected_frequency_matches_stationary_harmonic_fixture(self):
        heavy_mass, hydrogen_mass = 15.999, 1.00794
        reduced_mass = heavy_mass * hydrogen_mass / (heavy_mass + hydrogen_mass)
        curvature = 0.55
        step_angstrom = 0.005
        step_bohr = step_angstrom / ANGSTROM_PER_BOHR
        generalized_force = curvature * step_bohr
        plus = np.zeros((2, 3))
        minus = np.zeros((2, 3))
        # Equal and opposite Cartesian forces project exactly onto the
        # center-of-mass-preserving local coordinate.
        plus[1, 0], plus[0, 0] = generalized_force, -generalized_force
        minus[1, 0], minus[0, 0] = -generalized_force, generalized_force
        dipole_derivative = np.asarray([0.25, -0.10, 0.05])
        plus_dipole = dipole_derivative * step_bohr
        minus_dipole = -dipole_derivative * step_bohr
        result = projected_local_mode_from_differences(
            plus, minus, plus_dipole, minus_dipole,
            unit_vector=np.asarray([1.0, 0.0, 0.0]),
            heavy_atom=0, hydrogen_atom=1,
            heavy_mass_u=heavy_mass, hydrogen_mass_u=hydrogen_mass,
            step_angstrom=step_angstrom,
        )
        expected_frequency = HARTREE_VIBRATIONAL_WAVENUMBER_CM_1 * math.sqrt(
            curvature / (reduced_mass * ATOMIC_MASS_UNIT_IN_ELECTRON_MASSES)
        )
        derivative_debye_per_angstrom = (
            dipole_derivative * DEBYE_PER_ATOMIC_DIPOLE / ANGSTROM_PER_BOHR
        )
        expected_intensity = (
            IR_INTENSITY_CONVERSION_KM_MOL
            * float(np.dot(derivative_debye_per_angstrom, derivative_debye_per_angstrom))
            / reduced_mass
        )
        self.assertAlmostEqual(result["curvature_hartree_per_bohr2"], curvature, places=12)
        self.assertAlmostEqual(result["frequency_cm-1"], expected_frequency, places=9)
        self.assertAlmostEqual(result["intensity_km_mol"], expected_intensity, places=9)

    def test_invocation_ledger_counts_attempts_and_enforces_hard_cap(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "ledger.json"
            ledger = OrcaInvocationLedger(path, hard_cap=2)
            first = ledger.begin("job-a", "hash-a")
            ledger.finish(first)
            self.assertIsNone(ledger.begin("job-a", "hash-a"))
            second = ledger.begin("job-b", "hash-b")
            ledger.finish(second, error="failure")
            with self.assertRaisesRegex(RuntimeError, "budget exhausted"):
                ledger.begin("job-b", "hash-b")
            payload = json.loads(path.read_text())
        self.assertEqual(payload["invocations_used"], 2)
        self.assertEqual(payload["invocations_remaining"], 0)
        self.assertEqual([item["status"] for item in payload["invocations"]], ["completed", "failed"])

    def test_shared_ledger_namespaces_keep_representatives_independent(self):
        with tempfile.TemporaryDirectory() as temp:
            ledger = OrcaInvocationLedger(Path(temp) / "ledger.json", hard_cap=2)
            first = ledger.begin("environment-001:bond-000:minus", "same-input")
            ledger.finish(first)
            second = ledger.begin("environment-002:bond-000:minus", "same-input")
            self.assertIsNotNone(second)
            ledger.finish(second)

    def test_localized_oscillators_validate_against_coupled_hessian_subspace(self):
        local = [
            {"status": "validated", "frequency_cm-1": 3832.0, "intensity_km_mol": 17.0},
            {"status": "validated", "frequency_cm-1": 3832.0, "intensity_km_mol": 17.0},
        ]
        harmonic = [
            {"frequency_cm-1": 3787.3, "intensity_km_mol": 4.58},
            {"frequency_cm-1": 3881.77, "intensity_km_mol": 26.66},
        ]
        report = validate_local_modes_against_harmonic_subspace(local, harmonic)
        self.assertEqual(report["status"], "validated")
        self.assertLess(report["center_error_cm-1"], 3.0)
        self.assertLess(report["total_intensity_relative_error"], 0.10)

    def test_hessian_subspace_validation_fails_dimension_mismatch(self):
        report = validate_local_modes_against_harmonic_subspace(
            [{"status": "validated", "frequency_cm-1": 3500.0, "intensity_km_mol": 10.0}],
            [
                {"frequency_cm-1": 3500.0, "intensity_km_mol": 5.0},
                {"frequency_cm-1": 3600.0, "intensity_km_mol": 5.0},
            ],
            config=LocalModeHessianValidationConfig(
                maximum_subspace_center_error_cm_1=100.0,
                maximum_total_intensity_relative_error=1.0,
            ),
        )
        self.assertEqual(report["status"], "failed_validation")
        self.assertIn("subspace_dimension_mismatch", report["failures"])


if __name__ == "__main__":
    unittest.main()
