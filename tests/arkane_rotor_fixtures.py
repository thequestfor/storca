"""Synthetic statmech receipts for isolated unit tests; no native chemistry runs."""

import json
from pathlib import Path

from storca.arkane_rotor import (
    POLICY_ID, SYMBOLS, artifact, assess_geometry, create_rotor_contract, mint_rotor_receipt,
)


def rotor_validation_fixture(root, *, route_id="oh-methane", rate_output=""):
    """Retain real hash bindings while supplying a synthetic hydrogen audit.

    Tests may mock native execution, extraction, and model evaluation.
    Receipt acquisition and replacement admission use the portable validators.
    Native rotor loading and scientific mode preservation have separate tests.
    """
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    source_dir, staged_dir = root / "source", root / "species" / "hydrogen-artifacts"
    source_dir.mkdir()
    staged_dir.mkdir(parents=True)
    source = source_dir / "hydrogen.out"
    hessian = source.with_suffix(".hess")
    source.write_text(
        "O   R   C   A\nCARTESIAN COORDINATES (ANGSTROEM)\n"
        "---------------------------------\nH 0.0 0.0 0.0\n\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    hessian.write_text("synthetic atom Hessian fixture\n")
    staged = staged_dir / "calculation.out"
    staged_hessian = staged.with_suffix(".hess")
    staged.write_bytes(source.read_bytes())
    staged_hessian.write_bytes(hessian.read_bytes())
    statmech = root / "species" / "hydrogen.py"
    statmech.write_text(
        "energy = (0.0, 'kJ/mol')\n"
        f"geometry = Log({str(staged.resolve())!r})\n"
        f"frequencies = Log({str(staged.resolve())!r})\n"
        "spin_multiplicity = 2\nexternalSymmetry = 1\nopticalIsomers = 1\n"
    )
    (root / "input.py").write_text(
        "modelChemistry = 'mock native unit fixture'\n"
        f"species('hydrogen', {str(statmech.resolve())!r}, structure=SMILES('[H]'))\n"
        "statmech('hydrogen')\n"
    )
    provenance = {
        "schema_version": 1, "kind": "verified_arkane_route",
        "arkane_statmech_policy": POLICY_ID, "route": {"label": route_id},
        "orca_artifacts": [{
            "label": "hydrogen", "kind": "species", "formula": {"H": 1},
            "charge": 0, "multiplicity": 2,
            "output": str(source.resolve()), "hessian": str(hessian.resolve()),
            "staged_output": str(staged.resolve()), "staged_hessian": str(staged_hessian.resolve()),
            "statmech_file": str(statmech.resolve()),
            "checksums": {"output": artifact(source)["sha256"],
                          "hessian": artifact(hessian)["sha256"]},
        }],
    }
    (root / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    contract = create_rotor_contract(root)
    (root / "output.py").write_text(
        "conformer(label='hydrogen', E0=(0.0, 'kJ/mol'),\n"
        "          modes=[IdealGasTranslation(mass=(1.00782503223, 'amu'))],\n"
        "          spin_multiplicity=2, optical_isomers=1)\n"
        + rate_output
    )
    modes = {"mode_types": ["IdealGasTranslation"], "rotor": None}
    audit = {
        "label": "hydrogen", "status": "passed", "source_path": str(staged.resolve()),
        "geometry": assess_geometry([[0.0, 0.0, 0.0]], [1], [1.00782503223]),
        "action": "unchanged_atom", "optical_isomers": 1,
        "before": modes, "after": modes,
        "preserved": {name: True for name in (
            "nonrotor_modes", "energy", "unscaled_frequencies",
            "spin_multiplicity", "optical_isomers",
        )},
    }
    return mint_rotor_receipt(
        root, contract, [audit],
        {"rmg_version": "mock-native-unit-fixture",
         "source_artifacts": {"fixture": artifact(Path(__file__))}},
    )


def rotor_source_artifacts(statmech):
    """Return the actual rate artifacts bound into an existing fixture receipt."""
    receipt = json.loads(Path(statmech["receipt"]["path"]).read_text())
    return {name: receipt[name] for name in ("input", "output", "provenance")}


def rotor_rate_evidence_fixture(root, *, route_id="oh-methane"):
    statmech = rotor_validation_fixture(root, route_id=route_id)
    return {"arkane_statmech_validation": statmech,
            "arkane_rate_source_artifacts": rotor_source_artifacts(statmech)}


def conditional_rate_evidence_fixture(statmech):
    receipt = json.loads(Path(statmech["receipt"]["path"]).read_text())
    common = {"route_id": receipt["route_id"], "source_artifacts": rotor_source_artifacts(statmech),
              "upper_bound_claim_allowed": False}
    return {"conditional_forward_rate_evidence": {**common, "status": "conditional_central_model_rate"},
            "reverse_rate_evidence": {**common, "status": "unresolved_reverse_model"}}


def rotor_receipt_for_existing_input(root, rate_output):
    """Audit synthetic modes for already generated, fully retained unit inputs."""
    root = Path(root)
    contract = create_rotor_contract(root)
    provenance = json.loads((root / "provenance.json").read_text())
    declarations, audits = [], []
    for point in provenance["orca_artifacts"]:
        numbers = [SYMBOLS.index(symbol) for symbol, count in point["formula"].items() for _ in range(count)]
        masses = [float(number + 1) for number in numbers]
        coords = [[float(index), float(index * index), float(index % 2)] for index in range(len(numbers))]
        geometry = assess_geometry(coords, numbers, masses)
        rotor = None
        mode_text = "IdealGasTranslation(mass=(1.0, 'amu'))"
        mode_types = ["IdealGasTranslation"]
        if len(numbers) > 1:
            name = "LinearRotor" if len(numbers) == 2 else "NonlinearRotor"
            inertia = ([geometry["analytic_inertia_amu_angstrom2"]] if len(numbers) == 2
                       else [1.0, 2.0, 3.0])
            rotor = {"class": name, "inertia_amu_angstrom2": inertia,
                     "serialized_inertia_amu_angstrom2": inertia, "symmetry": 1, "quantum": False}
            mode_text += f", {name}(inertia=({inertia[0] if len(numbers) == 2 else inertia!r}, 'amu*angstrom^2'), symmetry=1)"
            mode_types.append(name)
        modes = {"mode_types": mode_types, "rotor": rotor}
        declarations.append(
            f"conformer(label={point['label']!r}, E0=(0.0, 'kJ/mol'), modes=[{mode_text}], "
            f"spin_multiplicity={point['multiplicity']}, optical_isomers=1)\n"
        )
        audits.append({
            "label": point["label"], "status": "passed", "source_path": point["staged_output"],
            "geometry": geometry, "optical_isomers": 1, "before": modes, "after": modes,
            "action": ("unchanged_atom" if len(numbers) == 1 else "diatomic_linear_rotor"
                       if len(numbers) == 2 else "unchanged_nonlinear_polyatomic"),
            "preserved": {name: True for name in (
                "nonrotor_modes", "energy", "unscaled_frequencies", "spin_multiplicity", "optical_isomers",
            )},
        })
    (root / "output.py").write_text("".join(declarations) + rate_output)
    return mint_rotor_receipt(
        root, contract, audits,
        {"rmg_version": "mock-native-unit-fixture", "source_artifacts": {"fixture": artifact(Path(__file__))}},
    )
