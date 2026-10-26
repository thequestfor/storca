"""Mapped endpoint identity evidence, distinct from reference conformation.

An XYZ has no bond table and does not certify an electronic state. This module
combines the existing conservative connectivity/stereo gate with an independent
coordinate-only bond assignment for ordinary neutral closed-shell fragments.
It can corroborate a changed reference conformation; it cannot certify a
stationary point, a reaction path, or a rate. Unsupported assignments remain
unresolved, rather than being declared a different molecule.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

from .route_geometry import (
    _bond_value,
    _build_endpoint_side,
    aligned_rmsd,
    read_xyz,
    validate_declared_connectivity,
)


ENDPOINT_IDENTITY_POLICY_ID = "mapped-coordinate-graph-conformer-distinction-v1"


def _reference_conformer(xyz_path: Path, reference_xyz: Path | None, maximum_rmsd: float) -> dict:
    record = {
        "kind": "mapped_reference_conformer_diagnostic",
        "identity_claim_allowed": False,
        "maximum_reference_rmsd_angstrom": maximum_rmsd,
        "same_reference_conformer": None,
    }
    if reference_xyz is None:
        return {**record, "status": "reference_not_supplied"}
    record["reference_xyz"] = str(reference_xyz)
    try:
        elements, coordinates, _ = read_xyz(xyz_path)
        reference_elements, reference_coordinates, _ = read_xyz(reference_xyz)
        if elements != reference_elements:
            return {**record, "status": "reference_atom_order_mismatch"}
        rmsd = aligned_rmsd(reference_coordinates, coordinates)
        if not math.isfinite(rmsd):
            raise ValueError("Aligned RMSD is not finite")
        matches = rmsd <= maximum_rmsd
        return {
            **record,
            "status": "reference_conformer_matched" if matches else "different_reference_conformer",
            "same_reference_conformer": matches,
            "aligned_rmsd_angstrom": rmsd,
            "reference_sha256": hashlib.sha256(reference_xyz.read_bytes()).hexdigest(),
            "interpretation": "Atom-ordered proper-rotation RMSD is conformation evidence, not chemical identity.",
        }
    except (OSError, ValueError) as error:
        return {**record, "status": "reference_conformer_unresolved", "reason": str(error)}


def validate_endpoint_identity(
    xyz_path: Path,
    smiles: tuple[str, ...],
    *,
    expected_to_xyz: dict[int, int] | None = None,
    reference_xyz: Path | None = None,
    maximum_reference_rmsd_angstrom: float = 0.75,
) -> dict:
    """Corroborate exact mapped chemistry without demanding one conformer.

The caller must supply the fixed declared-to-XYZ atom map. No graph isomorphism
or hydrogen reassignment is performed here. All declared stereochemistry still
passes the existing geometric check. Coordinate-only perception must agree in
every atom-ordered bond/order, formal charge and radical count. Charge/spin are
inherited declarations, not experimentally or electronically inferred here.
"""
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdDetermineBonds

    xyz_path = Path(xyz_path)
    reference_xyz = None if reference_xyz is None else Path(reference_xyz)
    result = {
        "schema_version": 1,
        "kind": "mapped_endpoint_identity_evidence",
        "policy_id": ENDPOINT_IDENTITY_POLICY_ID,
        "xyz": str(xyz_path),
        "declared_smiles": list(smiles),
        "valid": False,
        "chemistry_valid": False,
        "minimum_refinement_required": True,
        "stationary_point_claim_allowed": False,
        "electronic_state_claim_allowed": False,
        "path_claim_allowed": False,
        "rate_claim_allowed": False,
        "claim_scope": "coordinate_graph_corroboration_of_declared_mapped_chemistry",
        "fragment_checks": [],
    }
    if (isinstance(maximum_reference_rmsd_angstrom, bool)
            or not isinstance(maximum_reference_rmsd_angstrom, (int, float))
            or not math.isfinite(maximum_reference_rmsd_angstrom)
            or maximum_reference_rmsd_angstrom <= 0):
        return {**result, "status": "invalid_reference_conformer_threshold"}
    try:
        connectivity = validate_declared_connectivity(
            xyz_path, tuple(smiles), expected_to_xyz=expected_to_xyz,
        )
        result["connectivity_check"] = connectivity
        result["reference_conformer"] = _reference_conformer(
            xyz_path, reference_xyz, maximum_reference_rmsd_angstrom,
        )
        if not connectivity.get("valid"):
            return {**result, "status": "declared_chemistry_not_corroborated",
                    "reason": connectivity.get("status")}
        elements, coordinates, _ = read_xyz(xyz_path)
        result["xyz_sha256"] = hashlib.sha256(xyz_path.read_bytes()).hexdigest()
        side = _build_endpoint_side(tuple(smiles))
        mapping = ({index: index for index in range(len(elements))}
                   if expected_to_xyz is None else dict(expected_to_xyz))
        for fragment, declared in enumerate(side.explicit_molecules):
            declared_indices = [index for index, owner in enumerate(side.explicit_fragment)
                                if owner == fragment]
            check = {
                "fragment_index": fragment,
                "declared_smiles": side.smiles[fragment],
                "backend": "rdkit.DetermineBonds.connect_the_dots",
                "rdkit_version": rdBase.rdkitVersion,
                "maximum_bond_order_iterations": 1000,
                "declared_charge": side.formal_charges[fragment],
                "declared_multiplicity": side.inferred_multiplicities[fragment],
                "mapped_xyz_atom_indices": [mapping[index] for index in declared_indices],
                "isotopes_inherited_from_declared_atom_mapping": True,
                "electronic_state_claim_allowed": False,
            }
            result["fragment_checks"].append(check)
            if (side.formal_charges[fragment] != 0 or side.inferred_multiplicities[fragment] != 1
                    or any(atom.GetFormalCharge() != 0 or atom.GetNumRadicalElectrons() != 0
                           for atom in declared.GetAtoms())):
                check.update(status="assignment_outside_supported_domain",
                             reason="Coordinate-only corroboration is restricted to neutral non-zwitterionic closed-shell fragments.")
                continue
            block = [str(len(declared_indices)), "Independent mapped endpoint identity diagnostic"]
            for index in declared_indices:
                x, y, z = coordinates[mapping[index]]
                block.append(f"{side.elements[index]} {x:.12f} {y:.12f} {z:.12f}")
            try:
                perceived = Chem.MolFromXYZBlock("\n".join(block) + "\n")
                if perceived is None:
                    raise ValueError("Coordinate-only fragment could not be parsed")
                rdDetermineBonds.DetermineBonds(
                    perceived, charge=0, useHueckel=False, useVdw=False,
                    allowChargedFragments=True, embedChiral=False, maxIterations=1000,
                )
                Chem.SanitizeMol(perceived)
                expected_bonds = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))): _bond_value(bond)
                                  for bond in declared.GetBonds()}
                observed_bonds = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))): _bond_value(bond)
                                  for bond in perceived.GetBonds()}
                check["declared_atom_ordered_bonds"] = [
                    {"atom_indices": list(pair), "order": order}
                    for pair, order in sorted(expected_bonds.items())]
                check["perceived_atom_ordered_bonds"] = [
                    {"atom_indices": list(pair), "order": order}
                    for pair, order in sorted(observed_bonds.items())]
                states = []
                for expected, observed in zip(declared.GetAtoms(), perceived.GetAtoms()):
                    matches = ((expected.GetAtomicNum(), expected.GetFormalCharge(), expected.GetNumRadicalElectrons())
                               == (observed.GetAtomicNum(), observed.GetFormalCharge(), observed.GetNumRadicalElectrons()))
                    states.append({"atom_index": expected.GetIdx(),
                                   "declared_formal_charge": expected.GetFormalCharge(),
                                   "observed_formal_charge": observed.GetFormalCharge(),
                                   "declared_radical_electrons": expected.GetNumRadicalElectrons(),
                                   "observed_radical_electrons": observed.GetNumRadicalElectrons(),
                                   "matches": matches})
                check["atom_state_checks"] = states
                matches = (observed_bonds == expected_bonds
                           and perceived.GetNumAtoms() == declared.GetNumAtoms()
                           and all(item["matches"] for item in states))
                check["status"] = "exact_atom_ordered_graph_match" if matches else "perceived_graph_mismatch"
            except Exception as error:
                check.update(status="graph_assignment_unresolved", reason=f"{type(error).__name__}: {error}")
        statuses = {check["status"] for check in result["fragment_checks"]}
        if "perceived_graph_mismatch" in statuses:
            return {**result, "status": "mapped_chemical_identity_mismatch"}
        if statuses != {"exact_atom_ordered_graph_match"}:
            return {**result, "status": "mapped_chemical_identity_unresolved"}
        return {**result, "valid": True, "chemistry_valid": True,
                "status": "mapped_chemical_identity_corroborated",
                "interpretation": "Declared mapped chemistry is independently corroborated. Reference RMSD is only a conformer diagnostic; local optimization and frequency validation are still required."}
    except (OSError, ValueError) as error:
        return {**result, "status": "mapped_chemical_identity_unresolved", "reason": str(error)}


def validate_route_endpoint_identity(route, xyz_path: Path, *, side: str, mapping_result: dict,
                                     reference_xyz: Path | None = None) -> dict:
    """Route adapter using the same mapped atom order as existing endpoint gates."""
    if side == "reactants":
        smiles, expected_to_xyz = route.reactant_smiles, None
    elif side == "products":
        smiles = route.product_smiles
        expected_to_xyz = {int(product): int(reactant)
                           for reactant, product in mapping_result["all_atom_mapping"]}
    else:
        raise ValueError("side must be 'reactants' or 'products'")
    return validate_endpoint_identity(xyz_path, tuple(smiles), expected_to_xyz=expected_to_xyz,
                                      reference_xyz=reference_xyz)
