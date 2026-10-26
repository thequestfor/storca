# OH/methane constrained angular pilot

Status: **−15° independent derivative check passed; append-only continuation
stopped at +15° native MPI execution failure (2026-10-12). Fixed-geometry
serial/matched-MPI diagnostics passed. A geometry-only matched-MPI restart of
+15° converged at cycle 35 (2026-10-13). Its frequency-only continuation passed
the unchanged constrained physical/coupling gates after an append-only reader
correction. The local pilot passed; the full periodic profile remains pending**.
This is the pilot portion of step 2 in the
[bounded validation plan](BOUNDED_RATE_VALIDATION_PLAN_2026-10-11.md).
It does not complete the periodic angular/basin audit or promote a production rate.

## Question and implementation

The retained TS's 76.15 cm⁻¹ mode has substantial OH-H angular motion. The
experiment-specific [driver](../../tests/run_oh_methane_angular_audit.py) tests
whether nearby constrained hydrogen-transfer saddles remain well behaved on
the **same UKS WB97X-D3/def2-TZVPP, RIJONX/def2-J surface**. It does not fit a
potential or select a counting factor using the experiment.

Mapped atoms are O0, OH-H1, methane-C2, spectator-H3/H4/H5 and transferring-H6.
Use D(1,0,2,4), a signed OH-H–O–C–spectator-H dihedral. Its adjacent angle
sines are 0.99626 and 0.97119, mass-weighted overlap with the retained soft
mode is 0.81640, and overlap with the unstable transfer mode is 0.000760.
These are local coordinate diagnostics, not proof of a separable physical
rotor. The O–C axis is geometrical, not a covalent-bond declaration. Using
transfer-H6 as the outer reference is near-collinear and is rejected.

The [pure diagnostic module](../../storca/ts_angular_audit.py) checks angular
derivatives at two finite steps, removes rigid directions and constructs the
14-dimensional constrained tangent space. At a constrained stationary point
g(E)=λg(φ), it uses **H(L)=H(E)−λH(φ)**. Simply projecting the raw Hessian or
deleting the 76 cm⁻¹ mode would omit nonlinear-constraint curvature.

Each point requires an owned native completion, preserved mapped angle/state,
a fresh no-constraint central NumFreq Hessian at its exact optimized geometry,
a matching central gradient, one negative plus thirteen positive tangent
curvatures, and the intended opposing C–H/O–H transfer motion. Small eigenvalues
remain unresolved. Spin checks are continuity diagnostics, not proof that no
electronic-state switch occurred.

The pilot offsets are 0/−15/+15 degrees, not a complete periodic domain.
Orthogonal ZPE and modal partition changes are coupling diagnostics. No angular
partition, free-energy profile, optical/mirror count or rate is admitted.

## Attempts and auditable policy correction

All prior attempts are preserved:

- `r91v4lhr`: new input omitted the TS_Mode subsection `end`; ORCA rejected it
  before SCF. Corrected syntax has a regression test. The compact 6.1 keyword
  table omitted the terminator; the [full official syntax](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/geomopt.html)
  includes it.
- `479rd0_z`: corrected input was accepted but the sandbox blocked MPI startup;
  no SCF. The retry permits ORCA's local MPI workers without changing the method.
- `xqavjjg0`: center optimization completed normally, preserved the angle and
  reached ORCA's overachieved-gradient convergence. The first audit demanded
  all five printed criteria and refused before launching NumFreq. Four rows were
  printed: both gradient criteria passed strongly; both displacement criteria
  did not. **Do not report this as all-five convergence.**

The original revised, pre-off-center contract distinguished two candidate routes:
all-five convergence, or native overachieved-gradient convergence with both
gradient values ≤5% of their requested tolerance and displacements ≤2 times
their requested tolerance. Neither route verifies a saddle on its own. Both
still require the independently evaluated fixed-geometry gradient and numerical
Lagrangian-Hessian checks. This changes candidate acquisition, not the physical
stationarity/saddle gate or a production admission rule.

The center job is exactly input/output/XYZ/process bound and reused; it is not
recalculated to obtain a preferred convergence message. ORCA's reversed
D(4,2,0,1) primitive is accepted only because full atom-order reversal preserves
the signed dihedral and the native primitive is explicitly marked constrained.
No arbitrary atom permutation or reflection is admitted.

Preserved original campaign:
[zswvx1r1 contract](../../runs/oh-methane-angular-pilot-20261011-zswvx1r1/contract.json).
Participating source snapshots, original scientific files, predecessor failures,
native binary, initial Hessian and exact coordinate/method policy are bound.
The native ORCA bohr-to-Å conversion is recorded separately from SI constants.

The three earlier saddle attempts count toward the whole-increment limit.
The fresh pilot may dispatch two additional optimizations and three NumFreq
jobs: cumulative **five saddle attempts and three numerical Hessians**, within
the plan's 17/four limits. There are no automatic retries. Each native job has
four cores, 1000 MB MaxCore per rank, six hours and a 4 GiB disk reserve. A
resource stop returns unresolved, not chemical stability.

## Tests and remaining gate

Nineteen original focused tests pass: coordinate and derivative analytic checks,
rigid/mirror/reversal invariance, nonlinear Lagrangian curvature, small-mode
refusal, mapped transfer input, native convergence alternatives, no-constraint
frequency input, attempt accounting and retained failures. This is software
validation; native pilot results are still required.

## Candidate eligibility correction and continuation (2026-10-12)

The `zswvx1r1` center's independent numerical Hessian passed the physical
constrained-saddle gate. The −15° optimization then completed normally under
ORCA's native overachieved-gradient route, but the extra ≤5%-of-tolerance
candidate screen refused it **before** a numerical Hessian was evaluated.
Its RMS/MAX gradients were 1.5174e−6/6.4641e−6 Eh/bohr, both within the
requested 8e−6/3e−5 tolerances. Its RMS/MAX displacements were 1.42/1.62 times
their tolerances. Neither its optimizer Hessian nor this stopping decision
establishes whether the physical constrained saddle is valid.

The new [candidate receipt parser](../../storca/ts_candidate_eligibility.py)
removes only the extra 5% gradient screen. It requires final-cycle native
convergence evidence, gradients within their requested tolerances, and the
existing ≤2× displacement **acquisition guard**. Earlier-cycle alternate
markers cannot authorize the final geometry. It reports native failed step
criteria truthfully and returns eligible-for-derivatives, **not scientifically
verified**. Twenty focused parser tests include malformed/nonfinite rows,
completion ordering, alternate-marker locality and the actual −15° gradients.
Eleven continuation regressions additionally cover unchanged physical policy,
source-mutation refusal, exclusive redispatch protection, cumulative budgets and
failure accounting. Together with the original nineteen angular tests and
twenty-two electronic-state regressions, **72 focused tests pass**; this is not
a claim of a fresh full-suite run or scientific rate validation.

The [append-only continuation driver](../../tests/continue_oh_methane_angular_audit.py)
reopens the original contract, implementation snapshots, owned native leases
and exact inputs. It independently recomputes and compares the retained center
frequency receipt without new chemistry, reuses the −15° optimization, and
launches only its fresh fixed-geometry NumFreq. If that physical check passes,
it may complete the original +15° optimization/NumFreq point. It calls the
**unchanged frozen physical validator**: tangent RMS/MAX gradients ≤8e−6/3e−5,
one negative/thirteen positive constrained curvatures, matching H-transfer mode,
state/geometry binding and the same orthogonal-coupling limits. Failure stops
unresolved; no tolerance relaxation, method sweep or automatic retry follows.

The old driver and all old scientific artifacts remain unchanged. New
[continuation contract](../../runs/oh-methane-angular-continuation-20261012-ws7fj_vt/contract.json)
records both candidate policies, identical physical gates and the cumulative
starting ledger of **four saddle attempts/one NumFreq**. Successful completion
would total **five saddles/three NumFreq**, not a reset budget. An exclusive
source-contract claim prevents repeating this continuation or launching it
concurrently. All jobs retain four cores, six hours/job and the 4 GiB reserve.
Periodic coverage, angular partition, basin counting, rate and lifetime claims
remain unadmitted regardless of this three-point outcome.

If the pilot fails, retain unresolved angular assumptions and proceed to the
candidate-recipe freeze—no automatic multidimensional search or functional sweep.
If it passes, report `pilot_passed_profile_pending`, not `verified_rate`. A
separately bound continuation must establish periodic coverage, basin/counting,
the relaxed-path metric, orthogonal contributions, harmonic-limit recovery and
the plan's refinement criterion before considering any reduced partition model.

## Actual continuation outcome and execution diagnostic

The −15° fresh NumFreq passed the unchanged physical gate: one negative and
thirteen positive constrained curvatures, tangent RMS/MAX gradients
2.34138e−6/5.42299e−6 Eh/bohr, and retained H-transfer overlap 0.999593.
The center also replayed successfully. The +15° optimization then stopped at
cycle six: its SCF and printed gradient completed, but an MPI worker exited
improperly and ORCA reported error termination in the SCF gradient. The outer
driver's zero return code did not override missing native normal completion.
No +15° NumFreq ran. The cumulative angular ledger is **five saddle attempts
and two NumFreq**, not a completed three-point pilot.

The separate [fixed-gradient diagnostic](OH_METHANE_GRADIENT_MPI_DIAGNOSTIC_2026-10-12.md)
preserved this failed job and ran exactly two EnGrad calculations at its saved
full-precision geometry: serial and private package-matched Open MPI 4.1.1.
Both completed normally and agreed numerically. This is execution evidence,
not completion of the +15° optimization, an admitted saddle, or proof of crash
causation. No angular partition, counting factor, rate or lifetime was changed.

The separately authorized [matched-MPI continuation](OH_METHANE_MATCHED_MPI_CONTINUATION_2026-10-13.md)
reuses the exact tested cycle-six geometry and original verified initial Hessian
as an optimizer estimate. It does not reuse interrupted optimizer/orbital state.
It may dispatch one optimization and one fresh NumFreq, retaining all physical
gates and both original/restart angular targets. Center and −15° are replayed,
not recalculated. Maximum cumulative counts for this continuation are six
saddle attempts/three NumFreq; the two diagnostic EnGrad jobs stay separate.
