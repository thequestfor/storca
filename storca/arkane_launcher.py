"""Run native Arkane with a process-local, admitted ORCA rotor adapter."""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import json
import math
import os
from pathlib import Path
import sys
import tempfile

# Direct absolute-path launch in an RMG env does not require STORCA installed.
if __package__:
    from .arkane_rotor import (COLLINEAR_ROUNDOFF_RATIO, FAILURE_NAME, POLICY_ID, RECEIPT_NAME, SYMBOLS, artifact,
                              assess_geometry, mint_rotor_receipt, reopen, validate_contract)
else:
    from arkane_rotor import (COLLINEAR_ROUNDOFF_RATIO, FAILURE_NAME, POLICY_ID, RECEIPT_NAME, SYMBOLS, artifact,
                             assess_geometry, mint_rotor_receipt, reopen, validate_contract)


def register_adapter(contract):
    """Register a native subclass only in the invoking Arkane process."""
    contract = validate_contract(contract)
    from arkane.ess.orca import OrcaLog
    from arkane.ess.factory import register_ess_adapter
    from rmgpy.statmech import LinearRotor, NonlinearRotor

    points = {point["sources"]["staged_output"]["path"]: point for point in contract["stationary_points"]}

    def rotor_snapshot(conformer):
        rotors = [mode for mode in conformer.modes if isinstance(mode, (LinearRotor, NonlinearRotor))]
        if len(rotors) > 1:
            raise ValueError("Multiple native external rotors are not supported")
        mode = rotors[0] if rotors else None
        values = None if mode is None else mode.inertia.value
        if values is not None:
            values = [float(values)] if isinstance(mode, LinearRotor) else [float(value) for value in values]
            literal, units = ast.literal_eval(repr(mode.inertia))
            if units != "amu*angstrom^2":
                raise ValueError("Unsupported native inertia units/serialization")
            serialized = [float(literal)] if isinstance(mode, LinearRotor) else [float(value) for value in literal]
        return {"mode_types": [type(item).__name__ for item in conformer.modes],
                "rotor": None if mode is None else {
                    "class": type(mode).__name__, "inertia_amu_angstrom2": values,
                    "serialized_inertia_amu_angstrom2": serialized,
                    "symmetry": int(mode.symmetry), "quantum": bool(mode.quantum)}}

    class StorcaOrcaLog(OrcaLog):
        storca_audits = {}
        storca_contract = contract

        def __init__(self, path, check_for_errors=True, scratch_directory=None):
            path = str(Path(path).resolve())
            if path not in points:
                raise ValueError("ORCA log is outside the admitted STORCA rotor contract")
            self.storca_point = points[path]
            for source in self.storca_point["sources"].values():
                reopen(source)
            # Symmetry calls allocate their own scratch below. Caller-owned
            # scratch is never handed to the native deletion routine.
            super().__init__(path, check_for_errors=check_for_errors, scratch_directory=None)

        def get_software(self):
            return "Orca"

        def get_symmetry_properties(self):
            with tempfile.TemporaryDirectory(prefix="storca-orca-symmetry-") as temporary:
                previous = self.scratch_directory
                self.scratch_directory = str(Path(temporary) / "owned-symmetry")
                try:
                    return super().get_symmetry_properties()
                finally:
                    self.scratch_directory = previous

        def load_geometry(self):
            reopen(self.storca_point["sources"]["staged_output"])
            coords, numbers, masses = super().load_geometry()
            geometry = assess_geometry(coords, numbers, masses, allow_linear=True)
            composition = dict(Counter(SYMBOLS[int(number)] for number in numbers))
            if geometry["atom_count"] != self.storca_point["atom_count"] \
                    or composition != self.storca_point["formula"]:
                raise ValueError("Native ORCA geometry differs from the admitted atom composition")
            if self.storca_point.get("evidence_kind") == "monatomic_single_point" and (
                    len(masses) != 1 or not math.isclose(float(masses[0]), self.storca_point["atomic_mass_amu"],
                                                       rel_tol=0, abs_tol=1e-8)):
                raise ValueError("Native atomic isotope mass differs from its explicit declaration")
            return coords, numbers, masses

        def load_force_constant_matrix(self):
            if self.storca_point.get("evidence_kind") == "monatomic_single_point":
                # A single nucleus has no internal force matrix. This is a
                # typed not-applicable branch, never a missing-Hessian waiver.
                self.load_geometry()
                return None
            # Arkane 3.3's ORCA reader opens the discovered Hessian basename
            # relative to cwd. Keep that native read bound to this admitted
            # directory rather than any unrelated job-root Hessian. The
            # launcher is an isolated process and native loading is serial.
            hessian = reopen(self.storca_point["sources"]["staged_hessian"])
            previous = Path.cwd()
            try:
                os.chdir(hessian.parent)
                matrix = super().load_force_constant_matrix()
            finally:
                os.chdir(previous)
            reopen(self.storca_point["sources"]["staged_hessian"])
            return matrix

        def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=""):
            point = self.storca_point
            if label and label != point["label"]:
                raise ValueError("Native ORCA conformer label differs from its admitted source")
            conformer, frequencies = super().load_conformer(
                symmetry=symmetry, spin_multiplicity=spin_multiplicity,
                optical_isomers=optical_isomers, label=label,
            )
            coords, numbers, masses = self.load_geometry()
            geometry = assess_geometry(coords, numbers, masses, allow_linear=True)
            before = rotor_snapshot(conformer)
            nonrotor_modes = [mode for mode in conformer.modes if not isinstance(mode, (LinearRotor, NonlinearRotor))]
            energy, spin, optical = conformer.E0, conformer.spin_multiplicity, conformer.optical_isomers
            original_frequencies = list(frequencies)
            native_moments = None
            rotor = before["rotor"]
            if spin != point["multiplicity"] or optical != point["optical_isomers"]:
                raise ValueError("Native ORCA conformer spin/optical count differs from admitted declarations")
            if geometry["classification"] == "atom":
                if rotor is not None or before["mode_types"] != ["IdealGasTranslation"] or frequencies:
                    raise ValueError("Native atomic mode inventory is invalid")
                action = "unchanged_atom"
            else:
                if rotor is None or before["mode_types"] != ["IdealGasTranslation", rotor["class"], "HarmonicOscillator"]:
                    raise ValueError("Native ORCA RRHO mode inventory is incomplete or unsupported")
                if rotor["symmetry"] < 1 or not all(math.isfinite(value) for value in rotor["inertia_amu_angstrom2"]):
                    raise ValueError("Native ORCA external rotor metadata is invalid")
                if any(not math.isfinite(float(value)) or value <= 0 for value in frequencies):
                    raise ValueError("Native ORCA vibrational modes are invalid")
                if geometry["classification"] == "diatomic":
                    expected = 0 if point["kind"] == "transition_state" else 1
                    if len(frequencies) != expected:
                        raise ValueError("Native ORCA diatomic physical vibration inventory is invalid")
                    replacement = LinearRotor(inertia=(geometry["analytic_inertia_amu_angstrom2"], "amu*angstrom^2"),
                                              symmetry=rotor["symmetry"], quantum=rotor["quantum"])
                    conformer.modes = [replacement if isinstance(mode, (LinearRotor, NonlinearRotor)) else mode
                                       for mode in conformer.modes]
                    action = "diatomic_linear_rotor"
                else:
                    values = rotor["inertia_amu_angstrom2"]
                    linear = geometry["classification"] == "linear_polyatomic"
                    expected = 3 * geometry["atom_count"] - (5 if linear else 6) - int(point["kind"] == "transition_state")
                    if linear:
                        from arkane.common import get_principal_moments_of_inertia
                        native_moments = [float(value) for value in get_principal_moments_of_inertia(coords, numbers=numbers)[0]]
                        smallest, first, second = sorted(native_moments)
                        if rotor["class"] != "LinearRotor" or len(values) != 1 \
                                or not all(math.isfinite(value) for value in native_moments) \
                                or first <= 0 or second <= 0 \
                                or abs(smallest) > COLLINEAR_ROUNDOFF_RATIO * first \
                                or not math.isclose(first, second, rel_tol=1e-12, abs_tol=0) \
                                or not math.isclose(values[0], first, rel_tol=1e-12, abs_tol=0) \
                                or len(frequencies) != expected:
                            raise ValueError("Native linear polyatomic rotor classification is unresolved")
                        action = "unchanged_linear_polyatomic"
                    else:
                        if rotor["class"] != "NonlinearRotor" or len(values) != 3 \
                                or min(values) <= 0 or min(values) / max(values) <= 1e-8 \
                                or len(frequencies) != expected:
                            raise ValueError("Native polyatomic rotor or vibration classification is unresolved")
                        action = "unchanged_nonlinear_polyatomic"
            after = rotor_snapshot(conformer)
            preserved = {"nonrotor_modes": len(nonrotor_modes) == len([mode for mode in conformer.modes
                         if not isinstance(mode, (LinearRotor, NonlinearRotor))]) and all(left is right for left, right in
                         zip(nonrotor_modes, [mode for mode in conformer.modes
                              if not isinstance(mode, (LinearRotor, NonlinearRotor))])),
                         "energy": conformer.E0 is energy, "unscaled_frequencies": list(frequencies) == original_frequencies,
                         "spin_multiplicity": conformer.spin_multiplicity == spin,
                         "optical_isomers": conformer.optical_isomers == optical}
            if not all(preserved.values()):
                raise ValueError("ORCA rotor adapter changed an unrelated statmech quantity")
            self.storca_rotor_audit = {"label": point["label"], "source_path": self.path,
                                       "geometry": geometry, "before": before, "after": after,
                                       "spin_multiplicity": int(spin), "optical_isomers": int(optical),
                                       "preserved": preserved, "status": "passed", "action": action}
            if native_moments is not None:
                self.storca_rotor_audit["native_principal_inertia_amu_angstrom2"] = native_moments
            self.storca_audits[point["label"]] = self.storca_rotor_audit
            return conformer, frequencies

    register_ess_adapter("OrcaLog", StorcaOrcaLog)
    return StorcaOrcaLog


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("contract")
    arguments = parser.parse_args()
    workdir = Path(arguments.input).resolve().parent
    audits = {}
    try:
        contract = validate_contract(Path(arguments.contract).resolve())
        if Path(arguments.input).resolve() != reopen(contract["input"]) \
                or Path(arguments.contract).resolve() != workdir / "arkane-rotor-contract.json":
            raise ValueError("Arkane launcher input/contract path mismatch")
        for name in ("output.py", RECEIPT_NAME, FAILURE_NAME):
            if (workdir / name).exists() or (workdir / name).is_symlink():
                raise ValueError("Arkane launcher requires a fresh output/receipt workspace")
        adapter = register_adapter(contract)
        audits = adapter.storca_audits
        import rmgpy
        import arkane.common
        import arkane.ess.orca
        import arkane.ess.factory
        import arkane.input
        import arkane.statmech
        import rmgpy.statmech.rotation
        import rmgpy.quantity
        modules = (arkane.common, arkane.ess.orca, arkane.ess.factory, arkane.input,
                   arkane.statmech, rmgpy.statmech.rotation, rmgpy.quantity)
        native = {"rmg_version": str(rmgpy.__version__),
                  "source_artifacts": {module.__name__: artifact(module.__file__) for module in modules}}
        from arkane.__main__ import main as arkane_main
        sys.argv = ["Arkane.py", str(Path(arguments.input).resolve())]
        arkane_main()
        receipt = mint_rotor_receipt(workdir, contract, list(audits.values()), native)
        print("STORCA_ARKANE_STATMECH_JSON=" + json.dumps(receipt, sort_keys=True), flush=True)
    except BaseException as error:
        # Never turn partial output into a success receipt. Failure evidence is
        # append-only and cannot overwrite the result of an earlier attempt.
        failure = {"schema_version": 1, "policy_id": POLICY_ID, "status": "failed",
                   "failure_reason": str(error), "stationary_points": list(audits.values())}
        try:
            with (workdir / FAILURE_NAME).open("x") as handle:
                handle.write(json.dumps(failure, indent=2, sort_keys=True, allow_nan=False) + "\n")
        except (OSError, ValueError):
            pass
        raise


if __name__ == "__main__":
    main()
