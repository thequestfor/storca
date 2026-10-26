"""Bounded, labelled angular-profile diagnostics; no rotor or rate admission.

Lengths are explicit in each API. Gradient projection reuses the frozen pilot's
mass-weighted rigid/constraint conventions. An optimizer candidate is never
promoted to a verified saddle by these diagnostics or by an approximate Hessian.
"""
from __future__ import annotations

from functools import wraps
from itertools import permutations
import json
import math

import numpy as np

from storca.ts_angular_audit import _tangent, signed_dihedral, wrap_angle

GRID = (-170, -120, -60, -15, 0, 15, 60, 120, 170)
ACQUISITION_SCHEDULE = ((15, 60), (60, 120), (120, 170),
                        (-15, -60), (-60, -120), (-120, -170), (170, 190))
DIHEDRAL = (1, 0, 2, 4)
SPECTATORS = (3, 4, 5)
ACCEPTED_STATES = ("constrained_optimizer_candidate", "derivative_verified_constrained")
EH_J_MOL = 4.3597447222071e-18 * 6.02214076e23
GAS_CONSTANT_J_MOL_K = 8.31446261815324
ANGLE_TOLERANCE_RAD = 2e-4
TANGENT_RMS_TOLERANCE = 8e-6
TANGENT_MAX_TOLERANCE = 3e-5
CLOSURE_RMSD_A = .02
CLOSURE_MAX_DISPLACEMENT_A = .05
CLOSURE_ENERGY_KJ_MOL = .10
REQUIRED_PARTITION_EVIDENCE = (
    "periodic_coverage", "branch_closure", "candidate_gates_accepted",
    "selected_hessian_passed", "counting_reconciled",
    "kinetic_measure_derived", "reaction_normal_treatment", "transverse_bath_model",
    "harmonic_limit", "free_rotor_limit", "reparameterization_invariance",
    "orthogonal_coupling_accepted", "refinement_converged",
)


def _json_result(function):
    @wraps(function)
    def checked(*args, **kwargs):
        result = function(*args, **kwargs)
        # Reject overflow/nonfinite derived values rather than publish a NaN
        # receipt. This is serialization hygiene, not scientific admission.
        json.dumps(result, allow_nan=False)
        return result
    return checked


def _number(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f"Finite numeric {name} required")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"Finite {name} required")
    return value


def _geometry(xyz):
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 4 or not np.isfinite(xyz).all():
        raise ValueError("Finite mapped N x 3 geometry with at least four atoms required")
    return xyz


def _masses(masses_amu, count):
    masses = np.asarray(masses_amu, dtype=float)
    if masses.shape != (count,) or not np.isfinite(masses).all() or np.any(masses <= 0):
        raise ValueError("Matching finite positive atomic masses required")
    return masses


def _align(reference, candidate, masses):
    ref_center = np.average(reference, axis=0, weights=masses)
    new_center = np.average(candidate, axis=0, weights=masses)
    left, _, right = np.linalg.svd(((candidate-new_center)*masses[:, None]).T @ (reference-ref_center))
    correction = np.eye(3)
    correction[-1, -1] = 1 if np.linalg.det(left @ right) >= 0 else -1
    rotation = left @ correction @ right
    aligned = (candidate-new_center) @ rotation + ref_center
    return aligned, rotation


@_json_result
def proper_mapped_alignment(reference_xyz, candidate_xyz, masses_amu):
    """Proper mapped Kabsch alignment; input/output lengths are in angstrom.

    Atom correspondence is positional. No permutation or improper rotation is
    searched. A mismatch concerns branch geometry, not chemical identity.
    """
    reference, candidate = _geometry(reference_xyz), _geometry(candidate_xyz)
    if candidate.shape != reference.shape:
        raise ValueError("Matching mapped geometries required")
    masses = _masses(masses_amu, len(reference))
    aligned, rotation = _align(reference, candidate, masses)
    residual = aligned-reference
    displacement = np.linalg.norm(residual, axis=1)
    return {"aligned_candidate_xyz_A": aligned.tolist(), "rotation": rotation.tolist(),
            "rotation_determinant": float(np.linalg.det(rotation)),
            "residual_vectors_A": residual.tolist(), "atomic_displacements_A": displacement.tolist(),
            "mass_weighted_RMSD_A": float(np.sqrt(np.average(displacement**2, weights=masses))),
            "maximum_atomic_displacement_A": float(displacement.max()),
            "reflection_allowed": False, "permutation_applied": False,
            "chemical_identity_claim": False}


@_json_result
def branch_closure(reference_xyz, candidate_xyz, masses_amu, reference_energy_Eh,
                   candidate_energy_Eh, reference_gradient_passed, candidate_gradient_passed):
    """Strict labelled same-branch control, not a molecule-identity validator."""
    if type(reference_gradient_passed) is not bool or type(candidate_gradient_passed) is not bool:
        raise ValueError("Explicit boolean tangent-gradient decisions required")
    result = proper_mapped_alignment(reference_xyz, candidate_xyz, masses_amu)
    difference = abs(_number(candidate_energy_Eh, "candidate energy")
                     - _number(reference_energy_Eh, "reference energy")) * EH_J_MOL / 1000
    passed = (reference_gradient_passed and candidate_gradient_passed
              and result["mass_weighted_RMSD_A"] <= CLOSURE_RMSD_A
              and result["maximum_atomic_displacement_A"] <= CLOSURE_MAX_DISPLACEMENT_A
              and difference <= CLOSURE_ENERGY_KJ_MOL)
    return {**result, "absolute_energy_difference_kJ_mol": difference,
            "reference_gradient_passed": reference_gradient_passed,
            "candidate_gradient_passed": candidate_gradient_passed,
            "thresholds": {"mass_weighted_RMSD_A": CLOSURE_RMSD_A,
                           "maximum_atomic_displacement_A": CLOSURE_MAX_DISPLACEMENT_A,
                           "absolute_energy_difference_kJ_mol": CLOSURE_ENERGY_KJ_MOL},
            "strict_labelled_closure_passed": bool(passed), "domain_quotient_admitted": False,
            "scope": "Same labelled branch; failure does not mean a wrong molecule",
            "rate_claim_allowed": False}


@_json_result
def spectator_permutation_diagnostics(reference_xyz, candidate_xyz, masses_amu):
    """Enumerate only H3/H4/H5 relabellings, never admit a symmetry quotient."""
    reference, candidate = _geometry(reference_xyz), _geometry(candidate_xyz)
    if reference.shape != (7, 3) or candidate.shape != (7, 3):
        raise ValueError("Seven-atom OH/methane mapping required")
    masses = _masses(masses_amu, 7)
    if not np.allclose(masses[list(SPECTATORS)], masses[SPECTATORS[0]], atol=1e-12, rtol=0):
        raise ValueError("Spectator permutations must preserve isotope masses")
    records = []
    reference_angle = signed_dihedral(reference, DIHEDRAL)
    original_angle = signed_dihedral(candidate, DIHEDRAL)
    for ordering in permutations(SPECTATORS):
        mapping = list(range(7))
        for target, source in zip(SPECTATORS, ordering):
            mapping[target] = source
        permuted = candidate[mapping]
        alignment = proper_mapped_alignment(reference, permuted, masses)
        angle = signed_dihedral(permuted, DIHEDRAL)
        records.append({"reference_to_candidate_atom_indices": mapping,
                        "spectator_ordering": list(ordering),
                        "mass_weighted_RMSD_A": alignment["mass_weighted_RMSD_A"],
                        "maximum_atomic_displacement_A": alignment["maximum_atomic_displacement_A"],
                        "reference_angle_rad": reference_angle,
                        "original_candidate_angle_rad": original_angle,
                        "permuted_candidate_angle_rad": angle,
                        "angle_action_rad": wrap_angle(angle-original_angle),
                        "mapped_angle_gap_rad": abs(wrap_angle(angle-reference_angle)),
                        "rotation_determinant": alignment["rotation_determinant"]})
    best = min(records, key=lambda row: (row["mass_weighted_RMSD_A"], row["maximum_atomic_displacement_A"],
                                       row["spectator_ordering"]))
    return {"permutations": records, "best_spectator_ordering": best["spectator_ordering"],
            "reacting_atoms_fixed": [0, 1, 2, 6], "reflection_allowed": False,
            "domain_quotient_admitted": False, "counting_obligation_remaining": True,
            "rate_claim_allowed": False}


@_json_result
def candidate_gradient_diagnostics(xyz_bohr, masses_amu, gradient_Eh_bohr,
                                   indices=DIHEDRAL, *, target_angle_rad=None):
    """Frozen pilot projection without claiming a Hessian or native completion."""
    xyz = _geometry(xyz_bohr)
    masses = _masses(masses_amu, len(xyz))
    gradient = np.asarray(gradient_Eh_bohr, dtype=float)
    if gradient.shape not in ((xyz.size,), xyz.shape) or not np.isfinite(gradient).all():
        raise ValueError("Matching finite full Cartesian gradient required")
    gradient = gradient.ravel()
    _, _, normal, tangent, external, root, evidence = _tangent(xyz, masses, indices)
    weighted = gradient/root
    internal = weighted-external @ (external.T @ weighted)
    multiplier = float(normal @ internal/(normal @ normal))
    residual = (tangent @ (tangent.T @ internal))*root
    rms, maximum = float(np.sqrt(np.mean(residual**2))), float(np.max(abs(residual)))
    angle = signed_dihedral(xyz, indices)
    gap = (None if target_angle_rad is None else
           abs(wrap_angle(angle-_number(target_angle_rad, "target angle"))))
    passed = rms <= TANGENT_RMS_TOLERANCE and maximum <= TANGENT_MAX_TOLERANCE
    passed = passed and (gap is None or gap <= ANGLE_TOLERANCE_RAD)
    return {**evidence, "angle_rad": angle, "constraint_gap_rad": gap,
            "constraint_multiplier_Eh_rad": multiplier,
            "constraint_tangent_dimension": int(tangent.shape[1]),
            "tangent_gradient_RMS_Eh_bohr": rms, "tangent_gradient_max_Eh_bohr": maximum,
            "tangent_gradient_vectors_Eh_bohr": residual.reshape(xyz.shape).tolist(),
            "rigid_weighted_gradient_norm": float(np.linalg.norm(external.T @ weighted)),
            "candidate_gate_passed": bool(passed), "independent_hessian_verified": False,
            "scope": "Gradient/chart gate only; caller must bind native method/state/completion",
            "rate_claim_allowed": False}


def node_evidence_state(candidate_gate_passed, *, hessian_gate_passed=None):
    """The Hessian flag must come from the bound independent physical reader."""
    if type(candidate_gate_passed) is not bool or (hessian_gate_passed is not None and type(hessian_gate_passed) is not bool):
        raise ValueError("Explicit candidate and independent-Hessian decisions required")
    if not candidate_gate_passed or hessian_gate_passed is False:
        return "unresolved"
    return "derivative_verified_constrained" if hessian_gate_passed is True else "constrained_optimizer_candidate"


@_json_result
def evidence_counts(nodes):
    counts = {key: 0 for key in (*ACCEPTED_STATES, "unresolved")}
    accepted = []
    controls = 0
    for node in nodes:
        status = node.get("status")
        if status not in counts:
            raise ValueError("Explicit candidate/verified/unresolved node state required")
        counts[status] += 1
        control = node.get("closure_control", False)
        if type(control) is not bool:
            raise ValueError("Explicit boolean closure-control flag required")
        controls += int(control)
        if status in ACCEPTED_STATES and not control:
            accepted.append(node["offset_deg"])
    return {**counts, "closure_controls": controls,
            "accepted_unique_periodic_nodes": len(profile_coverage(accepted)["unique_offsets_deg"]),
            "verified_count_scope": "Counts labelled evidence; no inference from optimizer Hessians",
            "rate_claim_allowed": False}


def _canonical_offset(offset):
    value = (_number(offset, "angular offset")+180) % 360-180
    return 0.0 if abs(value) < 1e-10 else value


@_json_result
def profile_coverage(offsets_deg):
    groups = []
    for offset in offsets_deg:
        canonical = _canonical_offset(offset)
        group = next((group for group in groups if abs(wrap_angle(math.radians(canonical-group[0]))) < 1e-10), None)
        if group is None:
            groups.append([canonical, [float(offset)]])
        else:
            group[1].append(float(offset))
    groups.sort(key=lambda group: group[0])
    unique = [group[0] for group in groups]
    gaps = [(unique[(i+1) % len(unique)]-value) % 360 or 360.0 for i, value in enumerate(unique)]
    missing = [offset for offset in GRID if not any(abs(wrap_angle(math.radians(offset-value))) < 1e-10 for value in unique)]
    duplicates = [{"offset_deg": group[0], "supplied_offsets_deg": group[1]} for group in groups if len(group[1]) > 1]
    return {"unique_offsets_deg": unique, "duplicate_periodic_nodes": duplicates,
            "missing_grid_offsets_deg": missing,
            "coarse_grid_complete": not missing and not duplicates,
            "periodic_gaps_deg": gaps, "maximum_gap_deg": max(gaps, default=360.0),
            "seam_gap_deg": gaps[-1] if gaps else 360.0,
            "scope": "Labelled sampling coverage, not proof of narrow-feature absence or branch closure",
            "periodic_domain_quotiented": False, "rate_claim_allowed": False}


@_json_result
def select_final_hessian(nodes, temperature_grid_K, masses_amu):
    """Deterministic diagnostic ranking, never a physical thermal probability.

    Unit angular measure and electronic energies provide *provisional* Voronoi
    interval weights. Rank consequential wide candidates by branch ambiguity,
    nearest-pilot geometry/constraint-metric departure, then maximum interval
    weight over the whole temperature grid and ascending labelled angle.
    """
    accepted = [node for node in nodes if node.get("status") in ACCEPTED_STATES and not node.get("closure_control", False)]
    coverage = profile_coverage([node["offset_deg"] for node in accepted])
    temperatures = sorted(_number(T, "temperature") for T in temperature_grid_K)
    if not temperatures or min(temperatures) <= 0 or len(set(temperatures)) != len(temperatures):
        raise ValueError("Unique positive temperature grid required")
    result = {"selected_offset_deg": None, "provisional_electronic_interval_weights": [],
              "diagnostic_rank": [], "temperature_grid_K": temperatures,
              "rule": "Consequential wide unverified candidates: branch ambiguity, max(normalized geometry departure, inverse-metric log departure), maximum interval relevance, ascending angle",
              "weight_scope": "Electronic Boltzmann/unit-angle Voronoi heuristic, not a validated measure or physical probability",
              "potential_relevance_cutoff": .001, "rate_claim_allowed": False}
    if not coverage["coarse_grid_complete"]:
        return {**result, "status": "unresolved", "reason": "Complete unique coarse grid required before final-Hessian selection"}
    accepted.sort(key=lambda node: _canonical_offset(node["offset_deg"]))
    masses = _masses(masses_amu, len(_geometry(accepted[0]["xyz_A"])))
    energies = np.array([_number(node["energy_Eh"], "node energy") for node in accepted])
    gaps = np.asarray(coverage["periodic_gaps_deg"])
    widths = (gaps+np.roll(gaps, 1))/2
    weights = []
    for T in temperatures:
        raw = widths*np.exp(-(energies-energies.min())*EH_J_MOL/(GAS_CONSTANT_J_MOL_K*T))
        weights.append(raw/raw.sum())
    weights = np.asarray(weights)
    pilot = [node for node in accepted if abs(_canonical_offset(node["offset_deg"])) <= 15+1e-8]

    def metric(node):
        receipt = node.get("gradient_diagnostics") or node.get("frequency", {}).get("constrained_spectrum", {})
        value = _number(receipt["local_inverse_mass_metric_amu_bohr2"], "inverse constraint metric")
        if value <= 0:
            raise ValueError("Positive inverse-constraint metric required")
        return value

    ranks = []
    for index, node in enumerate(accepted):
        offset = _canonical_offset(node["offset_deg"])
        maximum = float(weights[:, index].max())
        result["provisional_electronic_interval_weights"].append({
            "offset_deg": offset, "weights_over_temperature_grid": weights[:, index].tolist(),
            "maximum_over_temperature_grid": maximum, "interval_width_deg": float(widths[index])})
        if abs(offset) < 60-1e-8 or node["status"] != "constrained_optimizer_candidate":
            continue
        geometry_departure = min(proper_mapped_alignment(ref["xyz_A"], node["xyz_A"], masses)["mass_weighted_RMSD_A"] for ref in pilot)/CLOSURE_RMSD_A
        metric_departure = min(abs(math.log(metric(node))-math.log(metric(ref))) for ref in pilot)
        ambiguity = node.get("branch_ambiguity", False)
        if type(ambiguity) is not bool:
            raise ValueError("Explicit branch-ambiguity flag required")
        ranks.append({"offset_deg": offset, "branch_ambiguity": ambiguity,
                      "normalized_geometry_departure": geometry_departure,
                      "inverse_constraint_metric_log_departure": metric_departure,
                      "departure_rank": max(geometry_departure, metric_departure),
                      "maximum_provisional_interval_weight": maximum,
                      "potentially_consequential": maximum >= result["potential_relevance_cutoff"]})
    ranks.sort(key=lambda row: (-int(row["potentially_consequential"]), -int(row["branch_ambiguity"]),
                               -row["departure_rank"], -row["maximum_provisional_interval_weight"], row["offset_deg"]))
    result["diagnostic_rank"] = ranks
    consequential = [row for row in ranks if row["potentially_consequential"]]
    if not consequential:
        return {**result, "status": "unresolved", "reason": "No wide unverified candidate passes provisional relevance screen"}
    return {**result, "status": "selected_for_diagnostic_hessian", "selected_offset_deg": consequential[0]["offset_deg"]}


@_json_result
def relaxed_path_inertia(angles_rad, geometries_bohr, masses_amu, index, *, reaction_mode=None, indices=DIHEDRAL):
    """Local three-point aligned path-tangent estimator in amu bohr².

    Angles must already be unwrapped. Neighbours are independently *properly*
    aligned to the central mapped frame; rigid velocity is removed. Optional
    Cartesian reaction-mode removal is explicit. This is neither the inverse
    constraint metric nor an admitted Hamiltonian/measure.
    """
    angles = np.array([_number(value, "unwrapped path angle") for value in angles_rad])
    frames = np.asarray(geometries_bohr, dtype=float)
    if (len(angles) < 3 or frames.ndim != 3 or frames.shape[0] != len(angles)
            or not np.isfinite(frames).all() or np.any(np.diff(angles) <= 0)
            or type(index) is not int or not 0 < index < len(angles)-1):
        raise ValueError("Matching finite path and interior unwrapped increasing angle required")
    center = _geometry(frames[index])
    masses = _masses(masses_amu, len(center))
    before, _ = _align(center, _geometry(frames[index-1]), masses)
    after, _ = _align(center, _geometry(frames[index+1]), masses)
    minus, plus = angles[index]-angles[index-1], angles[index+1]-angles[index]
    derivative = (-plus/(minus*(minus+plus))*before
                  + (plus-minus)/(minus*plus)*center
                  + minus/(plus*(minus+plus))*after)
    q, _, _, _, external, root, evidence = _tangent(center, masses, indices)
    weighted = derivative.ravel()*root
    rigid_removed = weighted-external @ (external.T @ weighted)
    selected = rigid_removed.copy()
    reaction_fraction = None
    if reaction_mode is not None:
        mode = np.asarray(reaction_mode, dtype=float)
        if mode.shape != center.shape or not np.isfinite(mode).all():
            raise ValueError("Finite mapped Cartesian reaction mode required")
        vector = mode.ravel()*root
        vector -= external @ (external.T @ vector)
        length = np.linalg.norm(vector)
        if length < 1e-12:
            raise ValueError("Nonrigid reaction mode required")
        vector /= length
        component = float(vector @ selected)
        denominator = float(selected @ selected)
        reaction_fraction = component**2/denominator if denominator > 1e-24 else 0.0
        selected -= component*vector
    inertia = float(selected @ selected)
    coordinate_derivative = float(q @ (selected/root))
    return {"inertia_amu_bohr2": inertia,
            "rigid_removed_inertia_amu_bohr2": float(rigid_removed @ rigid_removed),
            "local_inverse_mass_metric_amu_bohr2": evidence["local_inverse_mass_metric_amu_bohr2"],
            "aligned_path_derivative_bohr_per_rad": (selected/root).reshape(center.shape).tolist(),
            "angular_coordinate_derivative": coordinate_derivative,
            "coordinate_derivative_error": abs(coordinate_derivative-1),
            "unwrapped_neighbor_steps_rad": [float(minus), float(plus)],
            "reaction_mode_removed": reaction_mode is not None,
            "removed_reaction_projection_fraction": reaction_fraction,
            "metric_scope": "Finite-difference relaxed-path tangent estimate, not a validated kinetic measure",
            "hamiltonian_admitted": False, "rate_claim_allowed": False}


@_json_result
def partition_admission(evidence):
    """Report missing necessary evidence; never manufacture a rotor partition.

    Passing a checklist cannot replace a reviewed, implemented Hamiltonian.
    This increment can always close a useful profile-only assessment.
    """
    if not isinstance(evidence, dict):
        raise ValueError("Explicit partition evidence mapping required")
    missing = [name for name in REQUIRED_PARTITION_EVIDENCE if evidence.get(name) is not True]
    profile_ready = all(evidence.get(name) is True for name in (
        "periodic_coverage", "branch_closure", "candidate_gates_accepted",
        "selected_hessian_passed", "orthogonal_coupling_accepted"))
    return {"status": "profile_complete_partition_unsupported" if profile_ready else "profile_unresolved",
            "partition_admitted": False, "missing_necessary_evidence": missing,
            "quantitative_partition_model_implemented": False,
            "reason": "No reviewed reduced Hamiltonian/measure is implemented; evidence flags alone cannot admit one",
            "rate_claim_allowed": False, "lifetime_claim_allowed": False,
            "network_reinjection_allowed": False}
