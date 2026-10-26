"""Local constrained-saddle diagnostics, not a rotor or rate-model adapter.

Coordinates/derivatives use bohr and radians, energy derivatives Eh. Atomic
masses are in amu. A nonlinear angular constraint needs the Lagrangian Hessian
H(E) - lambda H(phi), not a deleted harmonic mode or projected H(E) alone.
"""
from __future__ import annotations

import math
import numpy as np

FREQUENCY_FACTOR_CM = (math.sqrt(4.3597447222071e-18 /
    (1.66053906660e-27 * 5.29177210903e-11**2)) / (2 * math.pi * 2.99792458e10))
DERIVATIVE_STEPS_BOHR = (2e-4, 1e-4)
EIGENVALUE_TOLERANCE = 1e-8  # Eh/(amu bohr²); small modes remain unresolved.


def wrap_angle(angle):
    if not math.isfinite(float(angle)):
        raise ValueError("Angular coordinate must be finite")
    return math.atan2(math.sin(float(angle)), math.cos(float(angle)))


def _geometry(xyz, indices):
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 4 or not np.isfinite(xyz).all():
        raise ValueError("Finite N x 3 geometry required")
    if (len(indices) != 4 or len(set(indices)) != 4
            or any(type(i) is not int or not 0 <= i < len(xyz) for i in indices)):
        raise ValueError("Four distinct mapped atom indices required")
    return xyz


def _chart(xyz, indices):
    xyz = _geometry(xyz, indices)
    i, j, k, l = indices
    left, axis, right = xyz[i]-xyz[j], xyz[k]-xyz[j], xyz[l]-xyz[k]
    lengths = np.array([np.linalg.norm(v) for v in (left, axis, right)])
    if min(lengths) <= 1e-8:
        raise ValueError("Angular chart has a zero-length axis or arm")
    axis /= lengths[1]
    u = left - np.dot(left, axis)*axis
    v = right - np.dot(right, axis)*axis
    sines = (np.linalg.norm(u)/lengths[0], np.linalg.norm(v)/lengths[2])
    if min(sines) < 0.15:
        raise ValueError("Angular chart is near-collinear or singular")
    angle = math.atan2(float(np.dot(np.cross(axis, u), v)), float(np.dot(u, v)))
    return angle, axis, sines


def signed_dihedral(xyz_bohr, indices):
    return _chart(xyz_bohr, indices)[0]


def rotate_dihedral_seed(xyz_bohr, indices, offset_rad):
    xyz = _geometry(xyz_bohr, indices).copy()
    angle, axis, _ = _chart(xyz, indices)
    offset = wrap_angle(offset_rad)
    i, j, _, _ = indices
    vector = xyz[i]-xyz[j]
    # phi measures the rotation from the first arm to the last arm, so rotate
    # the first arm by -offset to increase phi by offset.
    xyz[i] = xyz[j] + (vector*math.cos(offset) - np.cross(axis, vector)*math.sin(offset)
                        + axis*np.dot(axis, vector)*(1-math.cos(offset)))
    if abs(wrap_angle(signed_dihedral(xyz, indices)-angle-offset)) > 1e-10:
        raise ValueError("Angular seed rotation failed its coordinate witness")
    return xyz


def _derivatives(xyz, indices, increment):
    shape, flat = xyz.shape, xyz.ravel().copy()
    base = signed_dihedral(xyz, indices)
    gradient = np.zeros(len(flat))
    hessian = np.zeros((len(flat), len(flat)))

    def difference(displacement):
        return wrap_angle(signed_dihedral((flat+displacement).reshape(shape), indices)-base)

    # Only these twelve Cartesian components can affect this mapped dihedral.
    active = [3*atom+axis for atom in indices for axis in range(3)]
    shifts = {}
    for i in active:
        shift = np.zeros(len(flat))
        shift[i] = increment
        shifts[i] = shift
        plus, minus = difference(shift), difference(-shift)
        gradient[i] = (plus-minus)/(2*increment)
        hessian[i, i] = (plus+minus)/increment**2
    for index, i in enumerate(active):
        for j in active[:index]:
            u, v = shifts[i], shifts[j]
            value = (difference(u+v)-difference(u-v)-difference(-u+v)+difference(-u-v))/(4*increment**2)
            hessian[i, j] = hessian[j, i] = value
    return gradient, hessian


def angle_derivatives(xyz_bohr, indices):
    xyz = _geometry(xyz_bohr, indices)
    coarse_g, coarse_H = _derivatives(xyz, indices, DERIVATIVE_STEPS_BOHR[0])
    g, H = _derivatives(xyz, indices, DERIVATIVE_STEPS_BOHR[1])
    g_gap, H_gap = float(np.max(abs(g-coarse_g))), float(np.max(abs(H-coarse_H)))
    if (g_gap > 2e-7 + 2e-5*np.max(abs(g))
            or H_gap > 2e-6 + 2e-4*np.max(abs(H))):
        raise ValueError("Angular gradient/Hessian finite-step convergence unresolved")
    return g, H, {"steps_bohr": list(DERIVATIVE_STEPS_BOHR),
                  "gradient_max_step_difference_rad_bohr": g_gap,
                  "Hessian_max_step_difference_rad_bohr2": H_gap,
                  "finite_step_gate_passed": True}


def rigid_basis(xyz_bohr, masses_amu):
    xyz = np.asarray(xyz_bohr, dtype=float)
    masses = np.asarray(masses_amu, dtype=float)
    if (xyz.ndim != 2 or xyz.shape[1] != 3 or masses.shape != (len(xyz),)
            or not np.isfinite(xyz).all() or not np.isfinite(masses).all() or min(masses) <= 0):
        raise ValueError("Finite matching geometry and positive atomic masses required")
    centered = xyz - (xyz*masses[:, None]).sum(0)/masses.sum()
    root = np.sqrt(masses[:, None])
    external = [(np.tile(axis, (len(xyz), 1))*root).ravel() for axis in np.eye(3)]
    external += [(np.cross(np.tile(axis, (len(xyz), 1)), centered)*root).ravel() for axis in np.eye(3)]
    U, singular, _ = np.linalg.svd(np.asarray(external).T, full_matrices=False)
    if np.count_nonzero(singular > 1e-8*singular[0]) != 6:
        raise ValueError("Nonlinear geometry with six rigid directions required")
    return U[:, :6]


def _tangent(xyz, masses, indices):
    g, H, derivative = angle_derivatives(xyz, indices)
    root = np.repeat(np.sqrt(masses), 3)
    external = rigid_basis(xyz, masses)
    normal = g/root
    contamination = float(np.linalg.norm(external.T@normal)/np.linalg.norm(normal))
    if contamination > 1e-5:
        raise ValueError("Angular derivative violates rigid-motion invariance")
    normal -= external@(external.T@normal)
    magnitude = np.linalg.norm(normal)
    if magnitude < 1e-10:
        raise ValueError("Independent angular constraint normal missing")
    U, _, _ = np.linalg.svd(np.column_stack((external, normal/magnitude)), full_matrices=True)
    tangent = U[:, 7:]
    return g, H, normal, tangent, external, root, {**derivative,
        "external_contamination_fraction": contamination,
        "local_inverse_mass_metric_amu_bohr2": float(1/magnitude**2),
        "metric_scope": "local constraint inverse metric; not a proved relaxed-path kinetic metric"}


def angular_chart_diagnostics(xyz_bohr, masses_amu, indices, *, soft_mode=None, unstable_mode=None):
    xyz = _geometry(xyz_bohr, indices)
    masses = np.asarray(masses_amu, dtype=float)
    _, _, normal, _, _, root, evidence = _tangent(xyz, masses, indices)
    angle, _, sines = _chart(xyz, indices)
    result = {**evidence, "angle_rad": angle, "angle_deg": math.degrees(angle),
              "adjacent_angle_sines": list(sines), "constraint_tangent_dimension": 3*len(xyz)-7,
              "rate_claim_allowed": False}
    for label, mode in (("soft", soft_mode), ("unstable", unstable_mode)):
        if mode is None:
            continue
        vector = np.asarray(mode, dtype=float)
        if vector.shape != xyz.shape or not np.isfinite(vector).all() or np.linalg.norm(vector) == 0:
            raise ValueError("Finite mapped nonzero normal-mode vector required")
        weighted = vector.ravel()*root
        result[label+"_mode_overlap"] = float(abs(normal@weighted)/np.linalg.norm(normal)/np.linalg.norm(weighted))
    return result


def constrained_saddle_diagnostics(xyz_bohr, masses_amu, H_E_Eh_bohr2, gradient_Eh_bohr, indices):
    xyz = _geometry(xyz_bohr, indices)
    masses = np.asarray(masses_amu, dtype=float)
    size = xyz.size
    H = np.asarray(H_E_Eh_bohr2, dtype=float)
    gradient = np.asarray(gradient_Eh_bohr, dtype=float).reshape(-1)
    if H.shape != (size, size) or gradient.shape != (size,) or not np.isfinite(H).all() or not np.isfinite(gradient).all():
        raise ValueError("Complete finite Cartesian energy gradient/Hessian required")
    if np.max(abs(H-H.T)) > 1e-7:
        raise ValueError("Raw energy Hessian asymmetry exceeds the diagnostic tolerance")
    q, q_H, normal, tangent, external, root, evidence = _tangent(xyz, masses, indices)
    weighted_gradient = gradient/root
    internal_gradient = weighted_gradient-external@(external.T@weighted_gradient)
    multiplier = float(normal@internal_gradient/(normal@normal))
    residual = tangent@(tangent.T@internal_gradient)
    cartesian_residual = residual*root
    H_L = (H+H.T)/2-multiplier*q_H
    reduced = tangent.T@(H_L/np.outer(root, root))@tangent
    values, vectors = np.linalg.eigh((reduced+reduced.T)/2)
    frequencies = np.sign(values)*np.sqrt(abs(values))*FREQUENCY_FACTOR_CM
    negative = values < -EIGENVALUE_TOLERANCE
    positive = values > EIGENVALUE_TOLERANCE
    # Return the lowest direction even when it is not unstable. The negative
    # count gate must pass before any caller interprets this as a TS mode.
    unstable = ((tangent@vectors[:, 0])/root).reshape(xyz.shape)
    unstable /= np.linalg.norm(unstable)
    raw_reduced = tangent.T@((H+H.T)/2/np.outer(root, root))@tangent
    return {**evidence, "constraint_tangent_dimension": len(values),
            "constraint_multiplier_Eh_rad": multiplier,
            "gradient_convention": "g(E)=lambda g(phi); H(L)=H(E)-lambda H(phi)",
            "tangent_gradient_RMS_Eh_bohr": float(np.sqrt(np.mean(cartesian_residual**2))),
            "tangent_gradient_max_Eh_bohr": float(np.max(abs(cartesian_residual))),
            "rigid_weighted_gradient_norm": float(np.linalg.norm(external.T@weighted_gradient)),
            "eigenvalues_Eh_amu_bohr2": values.tolist(),
            "projected_frequencies_cm_inv": frequencies.tolist(),
            "positive_frequencies_cm_inv": frequencies[positive].tolist(),
            "negative_count": int(np.count_nonzero(negative)), "positive_count": int(np.count_nonzero(positive)),
            "unresolved_small_count": int(np.count_nonzero(~negative & ~positive)),
            "eigenvalue_tolerance_Eh_amu_bohr2": EIGENVALUE_TOLERANCE,
            "raw_vs_Lagrangian_reduced_matrix_max_difference": float(np.max(abs(reduced-raw_reduced))),
            "unstable_cartesian_mode": unstable.tolist(),
            "rate_claim_allowed": False, "lifetime_claim_allowed": False}
