# OH + methane: bounded periodic-profile implementation plan

Closeout (2026-10-17): all coarse nodes, independent labelled closure and the
selected +170° numerical Hessian completed and passed the unchanged checks.
The result is `profile_complete_partition_unsupported`, consuming 13/4 under
the original 17/4 caps. No quantitative partition refinement or new rate ran.
See [the closeout](OH_METHANE_CLOSEOUT_2026-10-17.md). Historical launch status
and pre-dispatch allocations below are retained as dated planning evidence.

Status: **bounded profile implementation complete; native validation running**
(2026-10-15). See the
[implementation checkpoint](OH_METHANE_PERIODIC_PROFILE_IMPLEMENTATION_2026-10-15.md).
The reviewed quantitative Hamiltonian is not available: the explicit
profile-diagnostics-only branch in section 3 is selected. The new rate consumer
records this refusal and preserves prior conditional evidence; it does not
invent a partition correction.

This implements the remaining angular-profile part of
[the bounded rate-validation plan](BOUNDED_RATE_VALIDATION_PLAN_2026-10-11.md),
then freezes the OH + methane assessment. It is not a new conformer platform,
a production rate admission, or a replacement for the multi-environment
[stability roadmap](../../STABILITY_ROADMAP.md).

## Scientific question and finish line

For OH + CH4 → H2O + CH3, determine whether the soft TS motion can be treated
as one connected angular degree of freedom without losing the H-transfer
saddle, miscounting equivalent configurations, or assuming an invalid thermal
partition. A successful implementation must also be able to conclude that
the evidence is insufficient for a corrected rate.

Mandatory deliverable: a bounded, reproducible angular-profile assessment,
its coverage and unresolved assumptions, followed by a frozen candidate recipe.
Conditional deliverable: a rate comparison **only if** the statistical-mechanical
model below passes its explicit admission tests. Agreement with experiment
is not an admission test. No lifetime or network-completeness claim follows.

## 1. Reuse the actual completed evidence

Start from `runs/oh-methane-angular-frequency-20261013-ivb3p6ba/result.json`,
status `pilot_passed_profile_pending`. Reopen its bound predecessor chain,
not just the top-level status. The 0°, −15° and +15° constrained points passed
the retained local physical/coupling checks. All three numerical Hessians
are reused; none is repeated merely to simplify orchestration.

The +15° **full** Hessian has an additional angular-like negative mode. Its
constraint-tangent Lagrangian Hessian has one negative/thirteen positive
curvatures. This is a constrained profile point, not an ordinary independently
verified first-order TS. The original IRC establishes one mapped reaction
route; it does not verify the whole angular surface.

Retain the atom map: O0, OH-H1, C2, spectator H3/H4/H5, transferred H6.
The coordinate is the signed dihedral `D(1,0,2,4)`; all offsets below are
relative to the retained center, not absolute ORCA dihedral values. The O–C
axis is geometric/nonbonded, not a conventional covalent torsional bond.

### Budget, including historical failures

| Resource | Whole increment cap | Already consumed | Remaining allocation |
| --- | ---: | ---: | --- |
| Constrained optimization attempts | 17 | 6 | 6 coarse + 1 closure + at most 4 refinement/retry |
| Full numerical-Hessian attempts | 4 | 3 | 1 selected wider-angle point |

Failures consume attempts. An unlaunched, durably reserved attempt remains
charged conservatively until its execution state is established. The two
earlier fixed-geometry serial/MPI EnGrad diagnostics remain separate recorded
costs. Do not reset either ledger when creating a new continuation directory.

Keep four cores, private witnessed Open MPI 4.1.1, MaxCore 1000 MB/rank,
six hours per native job and the existing 4 GiB free-space reserve. A timeout
is an execution limit, not chemical evidence. No new MPI build, benchmark or
resource-flag family is needed. Closing the laptop can still suspend progress.

## 2. Implement the smallest safe continuation

Proposed additions, with names finalized during implementation:

- `storca/ts_angular_profile.py`: pure profile/domain, candidate-gradient,
  branch-closure, metric and conditional partition diagnostics.
- `tests/run_oh_methane_angular_profile.py`: this experiment's bounded native
  controller, durable attempt reservations and evidence receipts.
- `tests/replay_oh_methane_angular_rates.py`: offline conditional rate overlay,
  only after the partition gate; use the existing declared `--rmg-env` when
  invoking the native Arkane environment.
- Corresponding focused unit/replay tests and one closeout checkpoint.

Reuse `ts_angular_audit.py` for angle derivatives, rigid-motion removal and
Lagrangian-Hessian diagnostics; existing pilot input/execution/readers;
`mpi_runtime_evidence.py`; and the composite/electronic-state replay helpers.
Do not edit frozen source versions to make old manifests pass. Bind the new
consumer revision and participating sources explicitly, preserving historical
hashes, mtimes, logs and failures. A source snapshot is not permission to
substitute arbitrary code for an old calculation.

The existing pilot executor has a local three-optimization cap. Reuse it as
a **one-job executor only after** the new durable outer reservation; its
one-job counter must never replace/reset the whole-campaign accounting.
Keep the owned streaming process/lease path, not a new unowned subprocess.

Separate immutable reservation, native execution and derived analysis receipts:

```text
reopen predecessors → reserve attempt exclusively → execute one owned job
                    → seal native receipt → derive/replay scientific receipt
```

On resume: completed native job + missing/broken reader receipt means offline
reanalysis, not redispatch. An ambiguous/running lease blocks redispatch;
never kill a process on PID alone. A native failure stays charged, with no
automatic retry. A second controller must not reserve the same point/Hessian.
Keep the user interface to the source/resume location and necessary environment
selection; methods, thresholds, atom maps and budgets belong in the contract.

## 3. Freeze the domain/counting policy before new calculations

Use a full signed 2π labelled chart initially; do not assume a 120° period.
Hold the reacting atoms fixed. Enumerate candidate spectator-H permutations
and record how each acts on the coordinate, but do not quotient the domain
without demonstrating the equivalence of the mapped branch. Proper rotations,
reflections and permutations are distinct operations.

One counting ledger must reconcile the existing reactant/product/TS external
symmetry numbers, four equivalent methane abstraction sites, any angular-domain
quotient, and mirror configurations. Do not multiply by a new factor of four,
divide by three, or add optical ×2 simply from visual symmetry. If a full-domain
integral already includes both mirrors, they cannot be counted again as
independent TS basins. Preserve the earlier one-/two-basin alternatives until
there is evidence to supersede them; never choose by experimental agreement.

Before dispatch, document whether a reduced partition is mathematically
supportable. If not, explicitly authorize **profile diagnostics only**, with
no corrected-rate output. This is not a reason to expand the campaign.

## 4. Acquire six wider points and one genuine closure control

Target the following nine unique initial nodes:

```text
−170, −120, −60, −15, 0, +15, +60, +120, +170 degrees
```

The center and ±15° are reused. The six new nodes are ±60°, ±120° and ±170°.
Grow two independent arms:
`+15 → +60 → +120 → +170` and `−15 → −60 → −120 → −170`.
Use the preceding accepted mapped geometry as the seed and the retained
numerical Hessian only as an optimizer acquisition estimate. Keep the native
method/state/constraint/transfer-mode prescription; do not use minimum
optimization or GOAT as a substitute for a constrained saddle search.

Reserve a seventh attempt: continue **+170 → +190°**, then compare with the
independently grown −170° point. The targets are periodically identical,
but the histories are different. This tests branch closure; rotating an old
geometry by 360° and rerunning it would only test repeatability.

This nine-node grid has a maximum 60° gap and a 20° seam gap. It is a bounded
screen, not proof that narrow features cannot exist between samples. An
alternative ±180° seam saves an optimization but supplies one fewer unique
wide node; prefer the explicit nine-node grid plus independent closure here.

Predeclare practical branch-closure guards: proper, mapped mass-weighted
alignment RMSD ≤0.02 Å, maximum atomic displacement ≤0.05 Å, and electronic
energy mismatch ≤0.10 kJ/mol, with both points separately passing the same
tangent-gradient gate. These are same-branch diagnostics, not energy accuracy
estimates. Test their numerical implementation before dispatch; do not loosen
them after inspecting a failed scan. Record raw mismatches as well as decisions.

Test strict labelled closure first. If only a spectator permutation matches,
retain both results and the corresponding counting obligation; do not silently
declare strict closure. Same connectivity but different relaxed shapes can
mean distinct branches, **not** that ORCA made the wrong molecule. An unresolved
seam blocks the single-branch partition, not a chemical identity claim.

## 5. Keep optimizer candidates separate from verified constrained points

At every new node require owned native completion, the same mapped atoms,
charge 0/multiplicity 2/19 electrons, method and resource contract; the intended
angle within 2e−4 rad; nonsingular angle derivatives; and preserved spectator
connectivity. Bind the final input/output/XYZ and final full native EnGrad at
the **same** geometry. Do not use an earlier cycle or dispersion-only gradient.

Use that gradient to test constrained stationarity: tangent RMS ≤8e−6 and
maximum ≤3e−5 Eh/bohr, with the existing derivative/rigid-motion conventions.
Retain native convergence-route and transfer-mode warnings. The existing
all-five/overachieved-gradient candidate rule remains unchanged; do not add
another stricter acquisition screen or claim every SCF residual met its
printed target. S² continuity is not proof of an unchanged electronic state.

An admitted optimization is `constrained_optimizer_candidate`. Its initial
or quasi-Newton Hessian cannot establish the final saddle index. Only an
independent full-Hessian receipt can label a node
`derivative_verified_constrained`.

Allocate the final numerical Hessian **after** coarse acquisition, using a
rule frozen before dispatch: among wide, potentially thermally consequential
points, prioritize unresolved branch/geometry/metric departures from the
three-point pilot, then interval-weight relevance over the whole retained
temperature grid, with a deterministic angular-order tie-break. If relevance
depends on an unvalidated measure, show that ranking as a provisional sampling
heuristic, not a physical probability. Never rank by experimental rate fit.

At this selected point enforce the unchanged physical checks: tangent index
one negative/thirteen positive, meaningful opposing C–H/O–H transfer mode and
overlap ≥0.90, ΔS² ≤0.05, orthogonal ZPE shift ≤500 J/mol and maximum absolute
Δlog Qorth ≤log(1.20) over the retained grid. A failure ends the admitted
one-coordinate model; do not weaken the gates or automatically buy more jobs.
A pass checks **one** wider location, not the unsampled spectra everywhere.

## 6. Gate any angular partition on an explicit physical model

Do not compute a corrected rate from a scalar `∫exp(−V/RT)dφ` alone. That lacks
the momentum normalization and kinetic measure. The existing inverse-constraint
metric is explicitly **not** a demonstrated relaxed-path inertia.

For a reduced model, specify its Hamiltonian/measure, boundary conditions,
rigid-rotation and reaction-normal treatment, aligned relaxed-path derivatives,
and interpolation. Distinguish it from a surface/coarea measure; these are
not interchangeable assumptions. Reject unsupported cross-couplings rather
than calling the result an exact rotor. Validate with synthetic coupled
coordinates before using chemistry output.

For seven nonlinear atoms there are fifteen internal degrees: one reactive
direction and fourteen thermal directions. Fixing this angle leaves fourteen
tangent directions; a valid constrained branch must have one unstable and
thirteen stable directions. Restore the angular thermal
degree once, exclude the reactive degree once, and recover the **full local
stable harmonic partition**, not merely the frequency nearest 76 cm⁻¹.
The observed mixed soft motion is not permission to delete a mode by index.
Specify electronic potential, orthogonal ZPE, angular ground energy, thermal
excitation and rotational changes separately; prevent duplicate ZPE/counting.

Arkane's documented `HinderedRotor` uses a torsional bond and a rotating top;
that is not automatically this nonbonded angular chart. Its external-rotor
handling is also distinct from the scoped `storca/arkane_rotor.py` admission
logic. See the [official Arkane input documentation](https://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/input.html).
Coupled torsional treatments are a relevant methodological precedent, not
validation of this proposed reduced coordinate; the derivation must be
reviewed before admitting a quantitative partition.

Require harmonic/free-rotor limits, units and coordinate-reparameterization
invariance, periodic seam uniqueness and a declared transverse-bath treatment.
With only four total Hessians, interpolation/frozen-bath assumptions at other
nodes remain explicitly conditional. Unresolved stability/coupling or an
unjustified metric means **diagnostic profile only**, not an invented interval.

Use at most four remaining optimization attempts, reduced by every failure
or retry, in one adaptive refinement batch. Select gaps using energy/gradient
shape, branch ambiguity and disagreement between admissible interpolants and
coarse/refined integrals across **all** retained temperatures. Never bridge
a missing branch by extrapolation. Require <10% partition change under the
declared refinement/interpolant checks. This is numerical convergence only;
it neither bounds mechanism error nor establishes 10% physical rate accuracy.

## 7. Reevaluate rates only after those gates, then freeze

Reuse the composite/electronic-state replay rather than write a second rate
pipeline. Keep the matched CCSD(T)/def2-QZVPP central energy; call the new
angular shape a **DFT-relative, QZ-anchored dual-level approximation**, not a
QZ angular surface. Leave ordinary `CompositeStationaryPoint` validation
unchanged: a constrained point with an extra full-Hessian negative mode cannot
be smuggled in as an ordinary TS.

Use a separately typed angular-partition overlay. Remove/replace the baseline
angular harmonic/ZPE contribution exactly once. Preserve native constants,
standard states, the retained frequency/ZPE convention and declared OH-level
energy offsets. Recompute both forward/reverse ground-energy barriers and
classical TST/Wigner/Eckart comparisons consistently; retaining an old Eckart
curve after changing its barriers is not acceptable. A central transmission
coefficient factorized over the angular model is an explicit **conditional**
approximation, not angle-resolved Eckart, VTST or multidimensional tunneling.

Check forward/reverse detailed balance and energy-zero invariance. Report
every retained temperature and counting/electronic-state branch. Compare to
the fourteen measurements only after freezing physical choices; those rows
remain one development study. Preserve the humid-air/finite-pressure versus
isolated-gas applicability warning. No basis/functional/counting choice is
selected by which curve is closest to experiment.

Freeze one of these outcomes, without automatic scope expansion:

1. `profile_unresolved`: domain, branch, execution or physical checks failed.
2. `profile_complete_partition_unsupported`: useful profile, no defensible
   reduced thermal model; retain prior conditional rates.
3. `conditional_angular_rate`: numerically converged declared approximation,
   with sparse Hessian coverage and physical assumptions visible.

None promotes production defaults, calibrated uncertainty, network reinjection,
environment support or a methane lifetime. The next separately authorized
campaign is condition-locked H + CH4, then OH + H2 transfer tests, followed by
the existing P1 fate regression—not indefinite OH + methane adjustment.

## Ordered implementation and acceptance checklist

- [x] Add the pure domain/closure/candidate types and the shared durable ledger.
- [x] Offline replay the actual three-point predecessor with native launch
  forbidden; preserve source/runtime identities and historical runtime labels.
- [x] Test interruptions before/after reservation, native completion and receipt
  generation; two controllers; ambiguous lease; failures; exhausted budget;
  final-Hessian double reservation; changed geometry/method/source/runtime.
- [x] Test mapped proper alignment versus reflection/permutation, seam closure,
  branch switching, angle singularity, final-EnGrad identity and candidate
  versus independently verified evidence states.
- [ ] Validate the proposed partition mathematics offline: harmonic/free-rotor
  and coupled-coordinate limits, relaxed versus inverse-constraint metric,
  rigid-motion removal, units, reparameterization and no unstable thermal mode.
  **Not admitted:** diagnostic path-metric invariance/unit tests are implemented,
  but no quantitative partition Hamiltonian or harmonic/free-rotor partition
  replay has been established. This intentionally blocks corrected rates.
- [x] Freeze contract/counting, closure thresholds, Hessian selection and
  refinement policy before any new job. Explicitly choose profile-only output
  if the partition derivation is not supportable.
- [x] Acquire six wide candidates and the independent closure control; select
  one remaining Hessian; refine only within the cumulative 17/4 caps.
  Completed with no refinement: the explicitly chosen profile-only branch
  does not admit quantitative partition refinement.
- [ ] If admitted, test exactly-once angular/ZPE/electronic counting, both
  barriers, tunneling consistency, energy-zero invariance and detailed balance
  against the existing replay. Never edit the ordinary-TS validator to pass.
- [x] Freeze result, all-outcome compute ledger, coverage table, source-bound
  recipe, raw evidence and unresolved assumptions; update the roadmap.

No new production CLI, discovery engine, multidimensional scan, solvent lane,
electronic-method sweep, extra QZ profile or generic rotor/conformer framework
belongs in this increment. Broader environments remain on the main roadmap;
this gas-phase benchmark does not validate them.
