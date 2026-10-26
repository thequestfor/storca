# OH + methane: direct curvature audit — 2026-10-10

Status: five fixed-frame energy/gradient calculations completed normally. Their
local curvature supports the two retained numerical Hessians, not the three
analytical Hessians. This local derivative diagnosis does **not** admit the
entire TS, an IRC, a rate or a lifetime. No production source or method default
was changed; all prior scientific artifacts remain unchanged.

## Question and bounded experiment

The [geometry/Hessian audit](OH_METHANE_GEOMETRY_HESSIAN_AUDIT_2026-10-10.md)
found the same weak OH angular displacement at approximately −31 cm⁻¹ in
analytical Hessians and +76 cm⁻¹ in two numerical Hessians. Tighter analytical
response convergence and a larger response grid did not settle the sign.
Independent raw-matrix reconstruction excluded a STORCA mode-label or matching
artifact, but did not decide which derivative treatment agrees with the energy.

This experiment uses the same refined neutral-doublet seven-atom geometry,
UKS WB97X-D3/def2-TZVPP with native D3(0), RIJONX/def2/J, ExtremeSCF and DefGrid3.
It computes a fresh center and displacements ±0.01 and ±0.02 bohr along the
**same Cartesian-normalized** Num001 weak vector. Cartesian displacement is not
mass-weighted displacement or a nonlinear angular-coordinate scan. Coordinates
remain in the retained frame; no alignment, constraint or optimization changes
the displaced structures. [ORCA's EnGrad run type](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basics.html)
supplies single-point energies and gradients, not stationary-point verification.

The original sandbox launch stopped at local MPI socket binding before SCF.
Its failed input/output/lease and closed result remain preserved. A fresh,
source-bound resource-only retry allowed local MPI communication with identical
scientific settings and center input bytes. Five jobs then completed normally
in 87.8 seconds total. They ran sequentially on four workers, with MaxCore
1000 MB/worker, six hours per job and a 4 GiB disk reserve. No guard fired;
about 60 GiB remained free. No additional user data were deleted.

## Results

For unit Cartesian vector u and amplitude h, the independent estimates are:

```text
κE(h) = [E(+h) + E(−h) − 2E(0)] / h²
κg(h) = u · [g(+h) − g(−h)] / (2h)
```

All curvatures below have units Eh/bohr².

| Evidence | Energy curvature | Gradient curvature |
| --- | ---: | ---: |
| Direct ±0.01 bohr | +2.324800×10⁻⁴ | +2.349439×10⁻⁴ |
| Direct ±0.02 bohr | +2.399000×10⁻⁴ | +2.498141×10⁻⁴ |
| Two-size local extrapolation | +2.300067×10⁻⁴ | +2.299872×10⁻⁴ |
| Retained Num005 full Hessian, uᵀHu | +2.302471×10⁻⁴ | — |
| Retained Num001 full Hessian, uᵀHu | +2.300081×10⁻⁴ | — |
| Retained analytical full Hessian, uᵀHu | −3.945097×10⁻⁵ | — |

The two-size extrapolation is `(4κ(0.01) − κ(0.02))/3`, assuming a smooth
leading h² truncation term. It is a post-run descriptive local estimate, not a
calibrated uncertainty interval or an independently demonstrated asymptotic
limit. Its energy and gradient values differ from Num001 by approximately
−1.4×10⁻⁹ and −2.1×10⁻⁸ Eh/bohr² respectively. Raw finite-amplitude estimates
also pass the preregistered sign and max(5×10⁻⁶ Eh/bohr², 10%) review threshold
against both numerical Hessians, while all analytical controls fail. These
thresholds allocate review; they are not physical accuracy certificates.

The fresh center reproduces the retained total energy, −116.259386991363 Eh,
and native S², 0.758947. Its Cartesian gradient norm is 3.026×10⁻⁶ Eh/bohr,
with directional projection −5.8785×10⁻⁷ Eh/bohr. The corresponding symmetric
odd energy slopes are −5.872×10⁻⁷ and −5.848×10⁻⁷ Eh/bohr. This linear drift
is why one-sided “energy goes up/down” reasoning is insufficient. Native
optimizer convergence uses projected/internal gradients, so the raw Cartesian
residual is not a contradiction of the earlier printed optimization criteria.
Finite-amplitude even-gradient drift is recorded, not silently discarded.

Input charge/multiplicity, nineteen electrons, integrated alpha/beta counts
near 10/9, atom order, geometry, native method/grid/dispersion, total energy,
printed Cartesian gradients and full-precision `.engrad` data agree within
their declared checks. S² remains 0.758947–0.758948 and atomic population screens
pass. These are electronic-state continuity diagnostics, not proof of unique
wavefunction identity, absence of multireference character or model accuracy.
All printed SCF residuals and tolerances are retained; native convergence is
not an assertion that every individual printed ExtremeSCF target passes.

## What this resolves—and what it does not

There is now independent energy **and** gradient evidence for positive local
curvature along the disputed displacement at the tested model/settings.
The analytical weak-mode curvature is inconsistent with that local evidence.
The result does not implicate RMG, Arkane or STORCA's normal-mode labels.
It also does not prove that this DFT energy surface accurately describes the
experimental chemistry.

Native D3 energy and printed dispersion gradients supply a separate component
check: the D3 curvature is approximately −9.63×10⁻⁶ Eh/bohr² and its two energy
estimates agree closely. That is much smaller than the approximately
2.69×10⁻⁴ analytical–numerical gap and has the opposite sign to the positive
total curvature. Simply missing that D3 contribution cannot explain the gap.
The retained analytical output includes a van der Waals Hessian completion
line, but no isolated component matrix; no specific native bug is identified
from this observation.

Full TS stationary-character/derivative admission remains false. One tested
direction does not validate all Hessian directions or remove the need for
fresh mapped bidirectional IRC evidence. No negative mode was deleted, no
rotor correction inferred, no original acceptance gate bypassed, and no
IRC, new correlated energies, Arkane rate, propagation or lifetime ran here.
The original approximately 5200-fold rate discrepancy remains unrecomputed.

## Next smallest repair step

Use this source-bound case as a minimal native derivative reproducer. Isolate
the analytical integration/derivative path with a bounded implementation or
native-version comparison before attributing the cause. The manual
[distinguishes SCF, response and Hessian grids](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/numericalintegration.html#changing-td-dft-cp-scf-and-hessian-grids);
the preceding response-only experiment was not an exhaustive Hessian-grid
convergence study. No package upgrade or new issue ticket was performed.

A future numerical-Hessian continuation must be an explicit, evidence-backed
alternative method policy, not “ignore the analytical negative”: repeated full
Hessians, meaningful reaction-coordinate character, justified stationarity,
fresh native/mapped IRC, and limits on applicable configurations remain gates.
Do not weaken all frequency checks or replace defaults project-wide from one
case. Consistent five-species energy/statistical-mechanical/rate construction
and experimental-condition matching are subsequent separate tasks.

## Reproducible evidence

- [Immutable five-point driver and read-only self-test/preflight](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/curvature.py)
- [Preserved sandbox startup failure](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/curvature-result.json)
- [Resource-only retry driver](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/curvature_retry.py)
- [Bound retry contract, settings, vector and inventories](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/curvature-mpi-retry-contract.json)
- [Closed results, raw data, SCF diagnostics and preservation receipts](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/curvature-mpi-retry-result.json)

The closed retry result has SHA256
`1d5be8b606b16e14785963a8690447af16063bb72e664bd8161fc4745a0caafc`.
Prior-artifact and bound-source byte/mtime preservation both pass. The isolated
driver self-tests and read-only preflight pass; this is not a full production
test-suite or benchmark calibration result.
