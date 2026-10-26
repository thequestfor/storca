# STORCA Stability and Chemical-Fate Roadmap

## Current implementation sequence

H + methane TS memory continuation (2026-10-23): the user authorized one
TS-only retry at 1500 MB/rank, four cores and cumulative SP cap 11 → 12.
The exact memory-only input is bound to the retained failed attempt, all four
other completed energies are reused, and only this TS may launch under the
new amendment. The memory guard rises with the allowance to 6.859 GiB;
all other time/resource/physical gates and historical charges remain.
460 combined regression tests passed. The prior 33 sources were archived;
the new 34-source revision is `runs/h-ch4-ts-memory-continuation-20261023`.
Post-test RAM was only 5.740 GiB, below the guard. The actual kitty invocation
refused at the memory guard before reservation: **zero new native jobs/charges**,
still 29 total attempts / 11 SP, none pending/running. Free roughly 1.2 GiB more
available memory plus headroom, then resume the same authorized revision; no
further cap increase is needed for this unspent retry. No rate is admitted.
See [the TS-memory handoff](docs/stability/H_CH4_TS_MEMORY_CONTINUATION_2026-10-23.md).

H + methane serial-H2 continuation (2026-10-23): the user authorized one core
for the exact H2 matched-energy retry and cumulative SP cap 10 → 11. All failed
charges and original resource/time/physical gates are retained; the conservative
four-rank memory guard is unchanged even for H2. All planned energy inputs now
undergo the known H2/CCSD process-count preflight before native dispatch.
451 combined regressions and the final three-test parallelism subset passed.
The prior 32 sources are archived; the new 33-source authorized revision is
`runs/h-ch4-h2-serial-continuation-20261023`. Native attempts 27 (H2, one core)
and 28 (methyl, four cores) completed and passed the strict energy reader.
Attempt 29 (TS) refused integral transformation: ORCA explicitly requires at
least 1359.7 MB/rank against the frozen 1000 MB. This is a native allocation
limit, not established macOS OOM or another geometry/IRC failure. Final counts
are 9 Opt/path, 7 Freq, 2 IRC, 11 SP, with none pending/running. Only the TS
matched energy remains; the cap is spent. Proposed next decision: TS-only
1500 MB/rank, four cores, one extra SP slot (11 → 12), retaining the memory
guard, all failed charges, cached energies and physical checks. Not yet
authorized; no rate is admitted.
See [the serial-H2 handoff](docs/stability/H_CH4_H2_SERIAL_CONTINUATION_2026-10-23.md).

H + methane PMIx continuation (2026-10-23): the user authorized the scoped
`PMIX_MCA_gds=hash` workaround and one extra single-point slot (9 → 10), retaining
all failed charges and original resource/time safeguards. The original verified
path, counting result and completed atomic energy are reused; only four matched
energies and the existing isolated rate audit remain. The control amendment is
explicit, input-bound and recorded per reservation, not a production default.
441 combined tests passed; a new wrapper output-directory startup omission was
then corrected and 56 focused tests passed. Both initial startup invocations
spent zero chemistry attempts. The executable revision is now
`runs/h-ch4-pmix-hash-continuation-20261023-r2`; the earlier preparation and its
driver bytes are retained as superseded zero-native evidence. Methane retry 25
terminated normally and passed the strict energy reader. H2 attempt 26 then
failed because four MDCI ranks exceed its one electron pair; normal-termination
admission correctly refused the top-level returncode 0. No methyl/TS energy ran.
Final counts: 9 Opt/path, 7 Freq, 2 IRC, 8 SP; none pending/running. This reveals
a rank-selection defect for a two-electron CCSD input. Continuing requires
separate authority for a one-core H2 exception and one extra SP slot (10 → 11),
not another four-core retry or a refund. No rate is admitted.
See [the PMIx handoff](docs/stability/H_CH4_PMIX_ENERGY_CONTINUATION_2026-10-23.md).

H + methane remaining-work repair (2026-10-23): the RAM-restored resume exposed
a planning defect: it reserved another IRC and atomic point even though both
were already paid. Remaining-cost forecasts now reopen owned native evidence,
without refunding any charges, and preserve the five matched-energy slots.
Resuming a branch action also retains its original operation and reuses its
existing assembly without rewriting historical continuation results. The
assembled path gets explicit, evidence-derived final validation fields.
All **424 combined regression tests pass**, including the actual retained
18-attempt interruption and source-lineage refusal tests. The authorized
revision is `runs/h-ch4-remaining-work-20261023-r3`; it retains the original
ledger and workspace through an exact archived-old/current-new source binding.
Native resume succeeded: both products passed reached-minimum validation, the
full route is verified and geometric counting passed. The atomic matched energy
completed (23); methane CCSD(T) then failed (24) during MPI/PMIx finalization.
The matching macOS report identifies recursive unfair-lock failure in the ds12
cleanup stack, not an established OOM event. Printed energy is not admitted
without normal termination. Nothing remains running/pending. Counts are now
9 Opt/path, 7 Freq, 2 IRC, 6 SP: three SP slots remain but four are needed.
No NEB/TS/IRC/reactant calculation was repeated and no rate is admitted. Further
native work requires a separately authorized budget/runtime continuation;
do not refund attempt 24 or retry automatically.
See [the repair handoff](docs/stability/H_CH4_REMAINING_WORK_REPAIR_2026-10-23.md).

H + methane recovery resume (2026-10-20): both owned fixed-geometry diagnostic
single points completed normally and passed all frozen electronic-equivalence
screens. TRAH minus original energy is -9.9476e-13 Eh, with identical printed
S² = 0.759080. This does not reproduce the lost DIIS history or establish an
ORCA bug. A separately sealed, authorized TS-only recovery is running under
`runs/h-ch4-partial-neb-recovery-20261020/campaign`, reusing the last evaluated
OptTS geometry without another NEB. OptTS converged in four cycles, and the
independent NumFreq completed with one imaginary mode at -1317.42 cm-1; the
mapped reaction-mode gate passed and bidirectional IRC is running. The same
ledger now includes attempt 14 (6 optimization/path, 4 frequency, 1 IRC,
3 single points). Four native workers and the unchanged six-hour per-job
allowance are verified. This is a running checkpoint, not a verified route/rate;
inspect live evidence before resuming.
See [the recovery resume handoff](docs/stability/H_CH4_RECOVERY_RESUME_2026-10-20.md).

Partial-NEB recovery increment (2026-10-19): attempt nine's NEB converged, but
ORCA's internal TS refinement aborted in electronic diagonalization. Pure paid-
stage admission, opt-in TS-only adaptive recovery and a two-point same-geometry
SCF comparison are implemented; 409 combined regression tests pass. The
original nine charges and all physical gates remain. Native comparison launches
were refused by the unchanged RAM guard, with zero new charges. Current
authorized diagnostic is `runs/h-ch4-scf-diagnostic-20261019-r3`; free memory,
then resume that package on the same ledger. Solver equivalence and TS recovery
remain unvalidated; an ORCA bug is not yet established. See
[the implementation handoff](docs/stability/H_CH4_PARTIAL_NEB_RECOVERY_2026-10-19.md).

Adaptive H + methane resume (2026-10-18): the same authorized campaign passed
the unchanged resource guard and started native NEB in kitty. The owned ORCA
lease and live MPI process inventory confirm four workers, with 21,600 seconds
allowed for the job. Original-ledger attempt nine is now charged (5 optimization/
path, 3 frequency, 0 IRC, 1 single point); the eight earlier charges remain.
The new native initial path and endpoint energies/gradients completed, and path
iterations are starting. No TS, route or rate is yet verified. The dated launch
snapshot is `runs/h-ch4-adaptive-authorized-20261018/resume-status-20261018.json`;
inspect live artifacts for subsequent status. The earlier memory refusal below
is preserved history, not the current running state.

Authorized adaptive H + methane launch (2026-10-18): the user authorized the
same-ledger continuation. The separately sealed authorized workspace is
`runs/h-ch4-adaptive-authorized-20261018`. Kitty opened and the controller
reused fragments/prepared the corrected seeds, but the unchanged memory reserve
refused the first native dispatch. **Zero new ORCA jobs or charged attempts**;
all eight historical charges remain, with no pending reservation. The last
available-memory check was 5,211,422,720 bytes versus the required 5,268,045,824
bytes, before another controller's overhead. Free memory, then resume this
same authorized campaign; do not lower safeguards or initialize another ledger.
Its `launch-status.json` contains the exact resume command. This is a resource
block, not evidence for or against the chemical route.

Adaptive route-search implementation (2026-10-18): the opt-in sequence in
[the saved plan](docs/stability/ADAPTIVE_VERIFIED_ROUTE_SEARCH_PLAN_2026-10-18.md)
is **implemented; native route success and efficiency remain unvalidated**.
Mapped paired geometry, persistent evidence-driven strategy selection,
paid-stage/IRC recovery and mandatory reached-species finalization use the
existing native executor and physical validators. The real zero-native H + CH4
replay corrects the inward/crowded approach: nearest spectator-H distance
changes from 1.0304 to 4.3728 Å, with both mapped endpoint identity screens
passing. This is geometry evidence, not a verified route or rate. All eight
historical charges remain. The current source-bound continuation is prepared,
**unauthorized**, under `runs/h-ch4-adaptive-preparation-20261018-r2` on the
original ledger. No new chemistry or production-default change occurred. See
[the implementation handoff](docs/stability/ADAPTIVE_VERIFIED_ROUTE_SEARCH_IMPLEMENTATION_2026-10-18.md)
for tests, remaining budget, unresolved claims and the next authorization.

H + methane stop (2026-10-18): the retry completed fragment work, then STORCA
interrupted NEB after 3033 seconds at its submerged-climbing-image diagnostic.
The finite encounter is crowded and is not a qualified asymptote, so the saved
dissociation-model classification is not established by that comparison. No
verified route/rate exists. Retain all eight charges (4 optimization/path,
3 frequency, 0 IRC, 1 single point), all native artifacts and the strict frozen
contract. The plan begins by preserving bound source bytes, qualifying endpoint
references, preparing better mapped arrangements and enabling an explicitly
source-bound continuation on the same native ledger. Earlier launch snapshots
below describe history, not current running status.

Authorized H + methane retry (2026-10-18): the user-authorized retry passed
the unchanged memory guard in the same campaign/ledger. The initial atomic-H
single point completed with native return code 0, and methane optimization is
running in kitty with four verified MPI workers and the six-hour job ceiling.
This is execution progress, not a validated path or rate. See
`runs/h-ch4-authorized-20261017/retry-status-20261018.json` for the launch
snapshot; all original attempt and cumulative caps remain unchanged.

Authorized H + methane launch (2026-10-18): a fresh authorized copy is retained
under `runs/h-ch4-authorized-20261017`, preserving the original preflight. Kitty
opened, but both bounded starts refused the first native dispatch at the
unchanged memory reserve; **zero ORCA jobs and zero charged attempts**. The
missing non-chemistry `psutil` dependency was declared/installed normally and
48 focused regression tests pass. The remaining blocker is available RAM, not
a chemical failure. Free memory, then resume the same campaign/ledger; do not
relax safeguards or open duplicate controllers. See its `launch-status.json`.

H + methane preflight followup (2026-10-17): the conditional H + CH4 model and
selected existing ORCA/MPI runtime are now frozen, **unauthorized**, at
`runs/h-ch4-preflight-20261017/campaign/execution-contract.json`. Case-specific
geometric symmetry/achirality is checked before expensive single points; no
coordinates/counts are repaired to pass. All 662 combined tests pass; actual
zero-job preparation/reopening is reproducible. No new chemistry ran. The next decision
is separate authorization of this bounded H + CH4 campaign, not another method
sweep or environment integration. The bounded primary-reference followup still
lacks full measured tables/correction, so experimental scoring remains blocked
independently of conditional prediction. See
[the preflight handoff](docs/stability/H_CH4_PREFLIGHT_2026-10-17.md).

Known-reaction executor increment (2026-10-17): the H + CH4 → H2 + CH3,
then OH + H2 adapter is implemented with preserved explicit hydrogen maps,
typed atomic-H evidence, complete opt-in native settings, durable categorized
compute control and reference-free composite/native-rate auditing. Synthetic
stage/restart tests and actual composite-reader/staging tests pass; the combined
631-test campaign and overlapping 458-test compatibility suites pass. No transfer
chemistry ran. The original 6 Opt/path / 4 frequency budget was insufficient
for mandatory reached-endpoint refinement; a separately documented 12/10/2/9
attempt proposal replaces it only for future authorized execution, without
altering the old preparation. Case-specific counting/state declarations,
runtime freeze, primary-reference/condition audits and separate native compute
authorization remain open. No physical accuracy/support gate is promoted. See
[the executor checkpoint](docs/stability/RATE_TRANSFER_EXECUTION_2026-10-17.md).

Latest bounded closeout (2026-10-17): the OH + methane periodic scan, labelled
closure control and selected +170° NumFreq completed normally. All unchanged
constrained-point/closure/local coupling gates pass, with cumulative spending
13 optimization attempts/4 NumFreq. The sealed outcome is
`profile_complete_partition_unsupported`: a reviewed angular kinetic measure
and global counting model are still absent, so no corrected rate, calibrated
uncertainty or lifetime is admitted. OH/methane tuning stops here. The
source-bound conditional candidate is frozen at
`runs/oh-methane-candidate-freeze-20261017-complete/candidate.json`; zero-dispatch
H + CH4, then OH + H2 preparation is retained at
`runs/rate-transfer-preparation-20261017`. All 424 focused tests pass. The next
benchmark is prepared, not executable/admitted: primary-table/conditions/
correction audits, case-specific counting/state assumptions, a known-reaction
execution adapter and separate compute authorization remain explicit gaps. See
[the closeout](docs/stability/OH_METHANE_CLOSEOUT_2026-10-17.md).

The accepted 2026-09-09 execution plan is
[Comprehensive, Working STORCA Stability](docs/stability/COMPREHENSIVE_IMPLEMENTATION_PLAN.md).
Use it for future implementation goals: close the prediction loop and start
scientific benchmarks first, then expand independent discovery and autonomous
solution chemistry with observable-driven allocation and method gates. It
preserves the milestones below while refining their execution order and making
neat-liquid and bulk-solid scope explicit. Existing milestone support states
are not promoted by the plan itself.

The plan's **Prediction-first execution order — 2026-10-01** section is the
current ordering: repair the real shared prediction loop while freezing and
running scientific benchmarks, produce applicable uncertainty rather than
requiring invented user bounds, then make ordinary solution chemistry autonomous.
Independent physics discovery and broader matter/exposure lanes advance through
measured complete applications. The section includes the first ten work items,
their dependencies and closure evidence. It changes planning priorities only;
the focused P1 acceptance gate and runtime support levels are unchanged.

The latest bounded shared-loop increment is
[Expanded-network assessment and next-evidence plan — 2026-10-08](docs/stability/EXPANDED_NETWORK_IMPLEMENTATION_PLAN_2026-10-08.md):
preserve full native identities, assess the saved expanded network separately
from its older certificate, diagnose branch influence and electronic-state
concerns, and retain explicit next-evidence requests. Its status is
**implemented_unvalidated**: software/native saved-model replay passes, but
physical rate accuracy, applicability and coverage remain unresolved. See
[the implementation checkpoint](docs/stability/EXPANDED_NETWORK_ASSESSMENT_2026-10-08.md).
It requires no new chemistry jobs and does not replace the wider environment
roadmap or close physical prediction gates.

Native-rate ingestion repair (2026-10-08): the STORCA-owned Arkane diatomic rotor
adapter passes a real saved OH + methane rate replay without new ORCA jobs or
installed-package edits. Source-bound native mode receipts are now required for
Arkane rate replacement. The corrected central rate passes the collision gate
but still strongly disagrees with the retained descriptive experimental
comparator; calibration and applicability remain unresolved. See
[the rotor checkpoint](docs/stability/ARKANE_DIATOMIC_ROTOR_FIX_2026-10-08.md).
The next distinct scientific task is a bounded barrier/rate-construction audit,
not support promotion or a new environment feature.

Basis/reference audit (2026-10-09): four new native correlated points completed.
The alternative open-shell reference changes the fixed-geometry TZ barrier by
−0.443893 kJ/mol; this is sensitivity evidence, not a calibrated error bound.
The canonical QZ transition state remains unresolved: the first attempt lacked
native-required memory, and a resource-only retry passed that gate but reached
the disk reserve. Completed points and original evidence are preserved. See
[the partial scientific checkpoint](docs/stability/OH_METHANE_BASIS_REFERENCE_AUDIT_2026-10-09.md).
Next: adequate scratch space for the missing TS/QZ point, then a bounded
geometry/Hessian/ZPE, soft-mode and electronic-partition assessment before
consistent five-species rate reconstruction. No rate, lifetime, calibration or
environment support gate was promoted.

Basis blocker resolved (2026-10-10): the missing canonical CCSD(T)/def2-QZVPP
OH + methane transition-state single point completed normally in 669 seconds
with two MPI workers and 2500 MB per worker. Only that point was dispatched;
completed inputs/outputs and all prior attempts remain unchanged. The electronic
barrier is 27.508540 kJ/mol, a −3.167869 kJ/mol TZ→QZ shift; adding the unchanged
low-level ZPE gives the diagnostic composite 20.969660 kJ/mol. These are
sensitivity results, not a calibrated accuracy certificate or a corrected rate.
QRO/QZ interaction remains untested. See
[the completed scientific checkpoint](docs/stability/OH_METHANE_BASIS_REFERENCE_AUDIT_2026-10-10.md).
Next remains bounded geometry/Hessian/ZPE, soft-mode and electronic-partition
assessment, then consistent five-species rate reconstruction. No support gate
or production method default was promoted.

Geometry/Hessian audit (2026-10-10): nineteen bounded native calculations
completed normally, including five B3LYP controls, independently refined
WB97X-D3/def2-TZVPP minima/TS and analytical/numerical derivative repeats.
The TS geometry converges, but analytical Hessians retain a weak negative
mode near −31 cm⁻¹ while two central numerical Hessians reproduce the same
mode at +76.19/+76.15 cm⁻¹. Tighter response residuals and a larger response
XC grid do not resolve the disagreement. First-order derivative/route
admission remains **unresolved**; conditional IRC preflight refuses before
dispatch. All prior scientific artifacts remain unchanged. See
[the geometry/Hessian checkpoint](docs/stability/OH_METHANE_GEOMETRY_HESSIAN_AUDIT_2026-10-10.md).
Next: a small independent energy/gradient curvature test along the disputed
mode, not mode deletion, an automatic rotor correction or a new environment
feature. No revised rate, lifetime, accuracy calibration, supported-lane gate
or production default is promoted.

Direct-curvature checkpoint (2026-10-10): five unchanged-model, fixed-frame
EnGrad calculations completed normally after a preserved pre-SCF sandbox MPI
startup failure and resource-only retry. Symmetric energy and gradient changes
at ±0.01/±0.02 bohr support positive curvature along the disputed angular mode;
two-size local extrapolations are approximately +2.30007/+2.29987×10⁻⁴
Eh/bohr², matching retained numerical Hessians rather than negative analytical
curvature. D3's approximately −9.63×10⁻⁶ contribution cannot explain the full
gap by simple omission. All prior scientific artifacts remain unchanged. See
[the direct-curvature checkpoint](docs/stability/OH_METHANE_CURVATURE_AUDIT_2026-10-10.md).
Next: isolate the native analytical derivative discrepancy, then define any
evidence-backed numerical-Hessian continuation explicitly. One-direction
diagnosis does not promote whole-TS/IRC/rate admission, calibration, supported
lanes or production defaults; no mode deletion or rotor/rate correction ran.

Numerical-route follow-up (2026-10-10): the conventional-integral NoRI control
completed six native jobs, with an append-only reader repair reusing its
normally completed frequency calculation rather than rerunning it. Stored
analytical curvature remains negative while same-control direct energy/gradient
curvature is positive. The NoRI geometry is nonstationary on its changed
surface, so this is diagnostic evidence, not a replacement TS/rate model or a
unique native-component diagnosis. Original-model full numerical repeats,
direct E/G corroboration, all-five native optimization criteria and fresh
Cartesian gradients pass a scoped alternative derivative admission. One fresh
bidirectional IRC completed with exact native numerical-Hessian read proof:
forward converged after 28 steps to OH + CH4, backward after 27 steps to
H2O + CH3, with the declared atom-mapped transfer verified. See
[the bounded continuation checkpoint](docs/stability/OH_METHANE_NUMERICAL_ROUTE_AUDIT_2026-10-10.md).
No production fallback/default, revised rate, lifetime, calibrated uncertainty
or environment support gate is promoted. All prior attempts remain preserved.

Consistent-rate reconstruction (2026-10-11): completed using
the five refined structures and admitted numerical TS frequencies, native
statistical-mechanics replay, and five matched
canonical CCSD(T)/def2-QZVPP single-point energies. Separate hash-bound energy
and geometry/vibration sources prevent old-geometry energy substitution or
native-log rewriting. The initial string-only omegaB97X-D3 Arkane declaration
failed before chemistry because the native method object requires an explicit
software identity. A fresh isolated recipe supplies ORCA explicitly; the
failed attempt is retained. All five QZ points terminate normally; native
composite E0 forward barrier is22.051678kJ/mol. At340.3K, direct Eckart k
is3.06232× the frozen measured comparator versus99.9241× for the matched
DFT-only baseline; across all14 measured temperatures it remains2.92–3.27×.
Source-bound composite energy/frequency handling and final hash-reopening
closeout pass. A precisely bounded parser revision is append-only admitted;
the original contract/native evidence are unchanged. See
[the completed rate audit](docs/stability/OH_METHANE_CONSISTENT_RATE_AUDIT_2026-10-11.md).
Next is physically consistent electronic-state/TS-basin counting and justified
soft-mode/tunneling assessment, not selection by experimental agreement.
Rate accuracy, electronic-state counting, soft-mode
harmonic validity, tunneling and experimental-environment applicability remain
independent questions. This work does not promote environment support or
production defaults.

Current focused increment (2026-10-11):
[Bounded rate validation, then return to chemical fate](docs/stability/BOUNDED_RATE_VALIDATION_PLAN_2026-10-11.md).
Complete a zero-ORCA electronic-state replay and a capped TS angular/basin
audit, freeze a candidate recipe rather than tune OH + methane indefinitely,
then test H + methane and OH + H2 before any scoped production handoff and
the existing P1 fate regression. Reference curation, all-outcome/cost accounting
and explicit stop gates are prerequisites. Step 1's conditional implementation
checkpoint follows; remaining steps are planned. Unresolved assumptions,
support levels and the broader multi-environment objective remain unchanged.

Electronic-state replay checkpoint (2026-10-11): the first bounded-plan step
is implemented and passes a fresh zero-ORCA native replay. Explicit two-level
OH populations, ground offsets, both barriers, recalculated Eckart and
corrected Gibbs equilibrium are consistent. At 340.3 K the conditional
composite Eckart comparison is 1.4827× for one harmonic TS basin or 2.9653×
for two assumed independent mirror basins; neither branch is admitted as the
true rate. Source snapshots and the unchanged historical evidence graph are
retained. See [the replay checkpoint](docs/stability/OH_METHANE_ELECTRONIC_STATE_REPLAY_2026-10-11.md).
Next is the capped angular/basin audit, then independent rate-transfer tests.
No production rate/default, uncertainty calibration, lifetime or environment
support gate is promoted.

Angular pilot checkpoint (2026-10-11): the experiment-specific constrained
TS pilot and numerical Lagrangian-Hessian diagnostics are implemented; nineteen
focused tests pass. The original pilot stopped before −15° NumFreq because
an extra candidate-gradient screen rejected native alternate convergence.
The 2026-10-12 append-only continuation removes that screen, replays/reuses
the center derivatives and reuses −15° optimization; fresh −15° NumFreq passed
the unchanged physical gates. The +15° optimization stopped at cycle six when
an MPI gradient worker exited improperly; no +15° NumFreq ran. The cumulative
angular ledger is five saddle attempts/two NumFreq. Twenty candidate-parser and eleven
continuation regressions pass; seventy-two focused tests pass including the
original angular and electronic-state regressions.
Prior syntax/MPI
failures and the overachieved-gradient center candidate remain preserved, with
attempt accounting and exact completed-center reuse. Full periodic coverage,
basin counting and a reduced partition/rate correction are not admitted. See
[the angular pilot checkpoint](docs/stability/OH_METHANE_ANGULAR_PILOT_2026-10-11.md).

Execution diagnostic checkpoint (2026-10-12): two separate fixed-geometry
EnGrad jobs at the saved +15° cycle-six frame passed in serial and with private
package-matched Open MPI 4.1.1. Actual worker library loading was witnessed;
energy and gradient components agree to approximately 1e−12 in native units.
Eighteen focused diagnostic tests pass. All prior artifacts/source inventories
remain unchanged. This does not establish the original crash cause, complete
the angular pilot, validate runtime performance, or admit a revised rate.
See [the gradient/MPI checkpoint](docs/stability/OH_METHANE_GRADIENT_MPI_DIAGNOSTIC_2026-10-12.md).

Matched-runtime continuation (2026-10-13): an append-only geometry-only +15°
restart is running from the exact gradient-tested frame, under private MPI4.1.1.
Center/−15° derivative evidence replays without new chemistry. One optimization
and one fresh NumFreq are allowed, for maximum cumulative six saddles/three
NumFreq; unchanged candidate/physical/coupling gates still apply. Seventeen
new focused tests pass, with prior 18 diagnostic/11 continuation tests rerun.
No periodic profile, rate, lifetime or support claim is promoted. See
[the matched-MPI checkpoint](docs/stability/OH_METHANE_MATCHED_MPI_CONTINUATION_2026-10-13.md).

Frequency-only closeout (2026-10-13): the +15° optimization completed normally
at cycle 35, passing all five candidate criteria. Its next-stage refusal was
a loader-parser defect, not a worker crash. A new pure reader recovers the
intact 140-worker/719-record evidence; the original failure receipt and code
remain frozen. One fresh NumFreq is running in an append-only continuation,
with zero new optimization dispatches and unchanged physical/coupling gates.
Twenty parser/thirteen completion tests pass, and seventeen prior restart
tests reran. The cumulative budget is six saddles/three NumFreq after dispatch;
no rate or global angular model is admitted. See
[the frequency checkpoint](docs/stability/OH_METHANE_FREQUENCY_COMPLETION_2026-10-13.md).

Native frequency closeout (2026-10-13): the fresh +15° NumFreq completed
normally. Independent raw replay confirms all three constrained pilot points
have one negative/thirteen positive tangent curvatures and pass stationarity,
state, mapped H-transfer and local orthogonal-coupling gates. The +15° full
Hessian has an additional angular-like negative mode: it is a constrained
profile point, not a second independently admitted ordinary TS. The pilot is
`pilot_passed_profile_pending`; periodic coverage, global measure/counting,
angular partition and revised rates remain unadmitted. Exactly one NumFreq/no
optimization was added, totaling six saddles/three NumFreq; eleven saddles/one
selected Hessian remain under the original caps. Fifty focused reader/driver
tests and ten angular-math tests pass. All prior source/artifact bindings
reopen unchanged. This does not promote production defaults or environment support.

Next bounded implementation (planned, 2026-10-13):
[Periodic angular profile and conditional rate assessment](docs/stability/OH_METHANE_PERIODIC_PROFILE_PLAN_2026-10-13.md).
Reuse all three completed pilot points, add six wider candidates and one
independent seam-closure control, and reserve at most four optimization
attempts for refinement/retry plus one selected numerical Hessian. Cumulative
caps remain seventeen optimization attempts/four NumFreq, including failures.
New optimization-only points are not independently saddle-verified. A rate
overlay requires explicit metric, counting, transverse-bath and harmonic-limit
checks; otherwise freeze a diagnostic/conditional assessment. No new chemistry
or production changes are authorized by this planning document.

Periodic-profile implementation (2026-10-15): the bounded controller, durable
whole-campaign journal, pure branch/metric diagnostics and offline rate-refusal
consumer are implemented. All three actual pilot points replay with zero
ORCA dispatch; 241 focused tests pass. Normal native work survives wrapper,
receipt and closeout interruptions without automatic redispatch. No reviewed
quantitative angular Hamiltonian is available, so this increment explicitly
chooses profile-only assessment, not a rate multiplier. Native wider-angle
validation is running in kitty from
`runs/oh-methane-periodic-20261015-ei3ezay9`; the first +60° optimization
reached cycle 2 on four cores at the launch check. The journal is seven
optimization attempts/three Hessians; no new scientific result is admitted.
All 554 source/runtime bindings reopen unchanged. See the
[implementation checkpoint](docs/stability/OH_METHANE_PERIODIC_PROFILE_IMPLEMENTATION_2026-10-15.md).

Candidate-policy correction (2026-10-16): the old +60° native optimization
completed normally in 36 cycles, but the extra 2× displacement acquisition
screen closed that run before later points. New append-only readers retain
native convergence and independent physical checks without that extra cutoff.
They also correct a force-frame assumption: the native property file binds
the last gradient to geometry 36, not final energy-only geometry 37. The actual
evaluated geometry passes unchanged stationarity/angle/state checks and is
reused without a new optimization. The correction carries forward 7/3 under
the original 17/4 caps; 300 focused tests pass. Native +120° is now running in
kitty at cycle 1 in `runs/oh-methane-periodic-candidate-20261016-ha4k1u9s`,
with cumulative spending 8/3. All 572 source/runtime bindings reopen unchanged.
No ordinary
TS, rate or lifetime claim is promoted. See
[the correction checkpoint](docs/stability/OH_METHANE_CANDIDATE_POLICY_CORRECTION_2026-10-16.md).

Runtime-reader correction (2026-10-16): +120° completed normally on four cores
in 43 cycles, then paused solely at the mandatory loader-log audit. An explicit
append-only reader extension makes detailed tracing diagnostic-only while
retaining selected-runtime inventories, owned normal completion and every
physical gate. The actual +120° evaluated geometry (43, not final energy-only
44) passes unchanged checks with zero new jobs and cumulative spending 8/3.
All 330 focused tests pass. The +170° calculation is now running in kitty
on four cores at optimization cycle 2; cumulative spending is 9/3. Frozen campaign sources and budgets
are not rewritten, and no rate/lifetime claim is enabled. See
[the runtime resume checkpoint](docs/stability/OH_METHANE_RUNTIME_RESUME_2026-10-16.md).

## Product objective

Source-admission/reference follow-up (2026-10-06): three actual native solution
candidates reach existing preparation, not TS/rate/propagation. A reproduced
corrupted-copy handoff now refuses changed native screen evidence while intact
sources still prepare, with zero new chemistry and original files unchanged.
Native source/state/energy/geometry checks are evidence admission, not physical
accuracy. The scientific panel gains one audited descriptive JPL limiting-rate
comparator; its theory-informed reference has no measured finite-pressure fate
or STORCA calibration, and the actual prediction remains missing. Ten curation
entries and the provisional P1 comparator remain unfinished. The working gate
tracks these changes separately from full prediction/environment support.

P1 checkpoint (2026-09-30): cyclobutene discovery completed under its frozen
development contract. A native replay revealed unintended Cantera self-heating;
the isothermal solver repair now requires actual T/P evidence at admission,
including reopened runs. Missing legacy evidence requires repropagation.
HPL Arkane bath-gas handling and parent-job timeout forwarding are corrected.
The real four-core parent ORCA check passed in 42.463 seconds and is eligible
for production reuse. The working-pipeline ledger tracks native checks, both
parent attempts and the remaining path/rate/repair stages. These repairs do
not close P1 or expand environment support.
Verification: 749 application tests (21 expected skips) and 32 native tests pass.

P1 checkpoint (2026-10-01): source-bound native RMG recipe traces now reach the
generic verifier, resolving the retained ring-opening atom-map ambiguity in
both directions without changing graph-edit weights or physical gates. The
actual-artifact replay preserves all original files and reuses its proofs.
Process-group cleanup now covers outer bridge/kinetics and RMG-supervisor
interruptions; the orphan reproducer and its correction are retained as tests.
Path/rate/repair and actual chemical-route resume remain the next gate. These
changes do not establish a lifetime, calibrated confidence, completeness or a
new supported environment. Current verification is tracked in the working
pipeline ledger, not inferred from the earlier checkpoint's counts.
Final check: 767 application tests (22 expected skips), 46 native tests, and
`git diff --check` passed. No full path/rate demonstration ran in this checkpoint.

Geometry-audit checkpoint (2026-10-01): a real production-chain retry reached
RMG discovery/mapping but rejected legitimate four-ring diagonals before any
ORCA path job. The exact retained seed/parent boundary now passes bounded
coordinate-graph corroboration without changing cutoffs; all 111 original files
remain unchanged. An expanded audit also closes demonstrated wrong-stereo,
bound-IRC identity, per-endpoint-tolerance and malformed-data gaps under an
explicit v9 verifier policy. The 20-method geometry panel is not an experimental
benchmark or proof of electronic applicability. See
[the brief failure report](docs/stability/GEOMETRY_VALIDATION_AUDIT.md) and
[the acceptance ledger](docs/stability/WORKING_PIPELINE_GATE.md) for current test
results and remaining compute. P1 still requires real path/rate/repair,
controlled route resume, uncertainty and coverage; no supported-lane promotion.
Final software check: 794 application tests (22 expected skips), 68 focused
native tests and `git diff --check` passed.

IRC/reference checkpoint (reconciled 2026-10-05): exact downhill execution and
retry ownership, source/job hashes, integer core-count policies and conflicting
source observations are now checked before branch assembly. A read-only native
replay reconstructs 119 byte-identical frames with 102 files/timestamps and
directories unchanged; its actual source has no execution retry. Retry,
fallback, restoration and mutation tests intercept execution, so they do not
prove physically executed retry or full path/rate/repair equivalence. JPL 2011
limiting-rate source curation corrects reference f(298) to f(300), not STORCA
prediction uncertainty. This checkpoint's completed software logs are application-007:
976 tests in 107.129 seconds, 22 expected skips; native-006: 68 tests in 9.725
seconds. These counts remain historical. No first-ten
scientific closure gate or environment support level is promoted. Evidence and
resource constraints remain in the working-pipeline ledger.

The latest completed software checkpoint includes solvent-bound native screening
and strict cache/phase-handoff repairs: 1018 application tests in 151.461 seconds
(22 expected skips), 115 expanded native/discovery-boundary tests in 9.695
seconds (including the original 68 checks) and 179 affected tests in 1.234
seconds. See [screening evidence](runs/solution-discovery-screen-diagnostic-001/EVIDENCE_2026-10-06.md).
The coefficient producer's 999/68 and pre-repair 1016/113, 1017/114 checkpoints
remain historical. These are software checks, not a positive solution path/rate
benchmark or supported-lane promotion; the first-ten scientific gates remain open.

STORCA should answer a conditional, auditable question:

> Given a chemical system, an environment, specified conditions, a requested
> observable, and a finite compute budget, what chemistry is likely to occur,
> on what timescale, through which mechanisms and products, and with what
> unresolved uncertainty?

The target is an evidence-driven multiscale chemical-fate engine. It must
discover plausible elementary chemistry from more than one source, determine
which reactions influence the requested observable, select the minimum
sufficient calculations needed to resolve those reactions, propagate verified
rates to macroscopic time, detect method inapplicability, and fail closed when
the evidence cannot support a prediction.

The architecture remains:

```text
chemical system + environments + observable + budget
  -> independent reaction-discovery evidence
  -> STORCA-native candidate network
  -> observable sensitivity and uncertainty allocation
  -> path and electronic-structure evidence
  -> applicable rate model
  -> repaired reaction network
  -> lane-specific propagation and convergence
  -> cross-environment report of claims, assumptions, and unresolved states
```

Electronic structure supplies evidence for elementary processes; kinetics
remains the bridge to macroscopic stability. AIMD, enhanced sampling, QM/MM,
and reactive potentials may discover mechanisms or free energies, but do not
replace the elementary-rate/network/propagation architecture.

## Stability campaign, lanes, and rungs

The user-facing stability ladder will evolve into a stability campaign. A
campaign contains several physically distinct lanes, while a lane may contain
an ordered sequence of condition rungs. Homogeneous gas, solution, surface,
and electrochemical results are not interchangeable and must not be collapsed
into a universal lifetime without an explicit exposure or aggregation model.

The existing bounded gas ladder and verified-input solution machinery remain
the first adapters. Future lanes become executable only after implementing and
validating their discovery, verification, rate, propagation, convergence,
uncertainty, provenance, and applicability contracts.

## Status vocabulary

- **planned**: accepted work that has not started.
- **in_progress**: implementation is incomplete.
- **implemented_unvalidated**: code exists, but its scientific acceptance gate
  has not passed.
- **scientifically_validated**: the declared benchmark and calibration gate has
  passed for a restricted scope.
- **observable_supported**: an admitted lane can issue the declared observable
  subject to all retained evidence and applicability gates.
- **blocked**: a recorded dependency prevents responsible progress.
- **deferred**: intentionally outside the current build horizon.

These roadmap states are more detailed than the runtime lane registry's
`declared_only`, `screening_supported`, and `observable_supported` levels. Code
completion must never be reported as scientific validation.

## Invariants that every milestone must preserve

- Fail closed when required evidence, applicability, convergence, or mechanism
  coverage is unresolved.
- Preserve the distinction between not found and shown not to exist.
- Treat missing pathways as first-class uncertainty.
- Keep RMG as a valuable discovery backend rather than the native network type.
- Keep gas rates, solution rates, photolysis rates, and surface rates physically
  distinct.
- Retain bounded search, restart, immutable contracts, calculation hashes,
  deterministic artifacts where possible, and explicit compute ledgers.
- Require endpoint identity, transition-state validation, meaningful imaginary
  modes, and IRC checks where conventional saddle-point kinetics applies.
- Do not force barrierless, spin-changing, multireference, solvent-assisted,
  pressure-dependent, photochemical, surface, or electrochemical chemistry
  through ordinary harmonic TST.
- Optimize calculations for a declared observable and uncertainty reduction,
  not calculation count.
- Never interpret agreement with one experimental observable as proof that the
  mechanism is correct.

## Milestone 0: Roadmap governance and traceability

**Status: implemented_unvalidated.**

Scientific hypothesis: stable requirement identifiers, explicit support gates,
and machine-checked traceability reduce silent scope loss and prevent planned
or partial physics from being advertised as supported.

Deliverables:

- this roadmap;
- the stability requirement catalog and traceability matrix;
- a definition of lane support;
- a machine-readable milestone and benchmark-profile manifest; and
- repository tests that compare support claims with the runtime lane registry.

Exit gate: every registered lane and roadmap requirement is represented in the
manifest, all referenced tests exist, and declared-only lanes cannot carry an
observable or benchmark support claim.

Do not build yet: a release-management service or external project-management
database. Repository artifacts remain the scientific source of truth.

## Milestone 1: Unified stability campaign

**Status: implemented_unvalidated.** Versioned chemical-system, observable, per-lane
environment, evidence-reference, budget, and campaign-request contracts are
available. `storca stability-campaign` writes an immutable preparation plan and
resumes only the same scientific contract. Explicit execution independently
runs ordered gas `t95` lanes through the existing bounded stability ladder and
restricted implicit-solution `t95` lanes through either retained verified
networks or the observable-conditioned solution verification session.

Scientific hypothesis: one immutable campaign can coordinate heterogeneous
lanes without weakening their independent physical assumptions or creating a
false universal lifetime.

Deliverables:

- a versioned `ChemicalSystem` and `StabilityCampaignRequest`;
- preparation-first `storca stability-campaign` orchestration;
- observable, environment, budget, cache, and resume contracts;
- adapters for the existing gas ladder and solution verification session;
- standardized per-lane terminal states and claim summaries;
- runnable preparation and synthetic regression examples; and
- a cross-lane report that preserves incomparable conditional results.

Implemented now: strict RDKit molecular identity and formal-charge validation;
one shared target with lane-specific environments; multiple observable and
same-lane-instance requests; scientific hashes that include methods, conditions,
budgets, assumptions, and evidence content hashes while excluding paths and
display metadata; exact preparation resume; gas, solution, combined, and
declared-only examples; a planner that launches no external programs; doubly
opt-in execution; enforceable gas and restricted solution budgets; shared
parent-minimum reuse plus condition-specific bounded gas discovery; gas stage
checkpoints; append-only solution verification
sessions; exact-contract adapter resume; and condition-scoped verified claims
without cross-environment aggregation. Synthetic mixed-lane regressions prove
the orchestration and fail-closed boundaries, but do not constitute experimental
validation of either physical model.

Milestone 1 hardening now also includes route-level gas verification
continuation. Fresh runs retain a hashed execution contract covering the parent
electronic state, conditions, method profile, ORCA resources, verifier policy,
engine bound, RMG configuration, and content-addressed local libraries.
`storca stability-resume-verification RUN_DIRECTORY` accepts no scientific
overrides, creates a new append-only engine session, reuses the retained
calculation workspace, and permits a scientifically deferred route to be
retried only once under a changed, explicit verifier-policy identifier. Legacy
runs may continue only when their independent artifacts establish the parent,
method, cores, timeout, conditions, resource limits, and absence of unretained
local libraries; any remaining historical defaults are recorded as inference,
not silently reconstructed.

Interrupted route execution is now distinguished from a new scientific retry.
Fresh artifacts and live ORCA process leases block duplicate launch. A quiet
retained `started` route can create a new append-only engine session while
reusing its workspace and policy-attempt identity. Complete ORCA NEB
`*_MEP.allxyz` bands are validated and resumed through a new checkpoint-hashed
restart input. Arkane and repaired-RMG mid-job checkpoint recovery remain
pending; their absence still produces an incomplete result rather than an
invented continuation.

Versioned surface-context upgrades are append-only and hash-linked. Confirmed
dead serial fallbacks can restore the retained parallel resource contract in a
new execution segment; live or indeterminate leases continue to block duplicate
launch. Endpoint-return warnings can reuse one validated NEB checkpoint without
being promoted to false intermediate species.

Exit gate: one request prepares and resumes a gas-plus-solution campaign,
rejects mutated requests, retains separate lane claims, and never promotes an
unadmitted lane or combines incompatible lifetimes.

Do not build yet: a graphical workflow editor or a universal ranking of
environments.

## Milestone 2: Autonomous multi-source candidate assembly

**Status: in_progress.** Native candidates, RMG adapters, bounded graph/xTB
discovery, immutable automatic assembly, source-preserving merging, strict
truncation/failure gates, solution candidate adapters, and an observable-bound
categorical coverage-risk/convergence gate exist. **STAB-008 is
implemented_unvalidated:** every gas stability run now enters the same
source-neutral assembly and coverage policy, so RMG-only success cannot bypass
independent-discovery accounting. The coverage artifact retains symmetric
source disagreement, explicit reaction-class scopes, independent bimolecular
and radical scope gaps, candidate influence states, applicability and
budget-expansion gaps, a separated uncertainty register, immutable content
identity/revisions, and stratified false-reassurance/false-abstention benchmark
metrics. Uncalibrated xTB endpoint-energy exclusions now block observable
convergence when they actually omit a native reaction, instead of becoming
silent negative evidence. A reaction independently retained by RMG remains in
the network. No runtime bypass is available until a hash-bound decision-scope
calibration contract exists.
Blind-spot calibration on curated chemistry, user-path orchestration,
and selective physics-driven escalation remain incomplete.

Retained-run scheduling repair: an exhausted direct queue from partial RMG
enumeration no longer prevents a condition-matched graph/xTB hypothesis using
initial reactants from reaching path verification. Solver and observable gates
remain required after repair. Regression tests and a replay of the NNONC
low-pressure run cover allocation; live chemistry completion remains open.

Live recovery checkpoint (2026-09-08): the retained stage
`runs/stability-ladder-20260831T180238Z/ladder/01-dark-low-pressure`
was continued in append-only `verification-engine-002`, retaining its original
condition contract and four-core configuration. ORCA passed the previous MPI
startup failure and reached converged SCF steps during reactant optimization.
This establishes execution recovery only, not endpoint validation, path
verification, rate credibility, or a supported lifetime. Inspect the descendant
engine state and calculation outputs before resuming again; this checkpoint
does not assert that the continuation has finished.

Recovery audit found an endpoint-specific provenance defect: the legacy
Opt/Freq workflow reran under the original basename, overwriting the earlier
startup-failure logs in this live continuation. The descendant engine record
is append-only, but those original raw endpoint logs are no longer intact;
do not describe this run as an unbroken raw-evidence chain. Reactant and both
product optimizations subsequently converged, and execution reached NEB-TS.
Endpoint Opt/Freq and minimum-recovery jobs now use the path executor through
an internal runner adapter. Startup retries retain the failed output and its
hash, reuse completed retry basenames, and refuse unknown incomplete output.
Regression evidence: 86 path/resume tests passed, including failed-log
preservation, completed-retry reuse, method-mutation rejection, and frequency
input built from the recovered geometry. The already running process predates
this fix; live validation of the new endpoint adapter remains open.

Follow-up contract audit: `_run_or_resume` now rejects changed input/dependency
hashes before replacing a retained contract or invoking ORCA. Endpoint input
generation and seed copying also refuse conflicting retained files. The
path/resume/foundation suite passes 240 tests, including method and referenced
geometry mutations, unchanged raw output/contracts after rejection, and
endpoint seed preservation. This hardens evidence identity; it is not a new
chemical capability or a scientific-validation milestone.

Live continuation outcome: `verification-engine-002` finished with two deferred
routes and no lifetime. The RMG route remained `surface_unresolved`; the
independent graph route was actually attempted and reached a refined TS with
one imaginary mode (-568.53 cm-1), followed by an IRC whose **both directions
exhausted 40 iterations**. Normal program termination did not mean IRC
convergence. The final frames also show an additional hydrogen transfer, so the
declared homolysis products cannot be assumed. A new directional IRC-completion
gate reports `irc_path_incomplete` before connectivity admission; 77 focused
tests pass, and a read-only replay of the retained IRC identifies the two
unconverged directions. The original session summaries remain unchanged.

Next live-failure work: bounded completion/relaxation of retained IRC branches,
then charge/spin/mapping-aware identity resolution of the endpoints actually
reached. An unexpected product must become a new unverified discovery
hypothesis, not be substituted into the original reaction's rate. Preserve the
TS/frequency evidence and do not rerun RMG/NEB merely to extend an incomplete
IRC. Both route repair and network propagation remain open.

Downhill continuation checkpoint: an internal, append-only diagnostic now
extends a capped branch using ORCA's documented `Direction down`, with the
original 40-step bound and four-core model retained for this case. Both branches
were queued in one sequential execution from the saved graph-route IRC under
`transition-state/irc-branch-continuation/{forward,backward}`. The forward
extension has entered the downhill iteration table; inspect actual process and
result evidence before starting anything again. Focused IRC/path/resume tests
pass (93 tests). Automatic trajectory joining, endpoint-identity admission, and
integration into versioned ladder recovery are still open; the continuation
diagnostic itself cannot issue a path, product, rate, or lifetime claim.

Branch joining is now implemented as a non-executing, hash-checked assembly
step. Synthetic bidirectional tests cover correct frame ordering, completed
artifact reuse, rejection of a later unconverged segment, and discontinuous
junctions. A read-only check against the live extension correctly returned
`downhill_branch_not_converged`. Live joining and automatic versioned recovery
wiring remain open; do not report this as completed endpoint or rate validation.

Live IRC completion and v7 wiring: both downhill extensions converged and the
checked assembler produced a 119-frame trajectory in
`irc-branch-continuation/assembled-4d154932259cd749`. Rechecking the final
endpoints confirmed the reactant on one side, but the other side includes
O-N cleavage **and hydrogen transfer**, not the declared radical homolysis
products. The retained `endpoint-validation.json` therefore correctly remains
invalid for that original reaction; no rate or lifetime was admitted.
`generic-orca-path-v7-irc-branch-recovery` now performs the bounded extensions,
assembly, and endpoint recheck automatically in TS verification. v2-v6 runs
enter through explicit hashed amendments, preserving the old execution
contract. Unexpected-channel identity resolution and discovery handoff are
the next unresolved live-failure task, not another attempt to force the
original homolysis route through rate generation.

Unexpected-channel proposal checkpoint: the live endpoint now produces native
candidate `reaction-2cdf1d1a16c24c3a`, with provisional mapped products
`[NH:1]=[NH+:2][O-:3]` and `[NH2:4][CH3:5]`. The proposal retains its parent
route, condition, source hashes, and observed O-N cleavage/H-transfer atom
indices. It does not claim those fragment electronic states or minima are
validated. Canonical hydrogen mapping selected a symmetry-equivalent hydrogen
different from the one in the actual trace, so atom-order alignment remains a
real prerequisite for reusing the retained TS. The original homolysis remains
unresolved. The candidate is stored in the live assembled-IRC directory as
`unexpected-endpoint-discovery-002.json`; earlier evidence is preserved.
The native proposal regression uses the retained final coordinates and passes
with the IRC/path/discovery suite (98 tests). The preceding full suite passed
631 tests before this proposal addition. Next: trace-consistent all-atom mapping,
stationary fragment/state checks, bounded condition-matched queue handoff, and
retained-TS reuse before any new rate or network propagation.

Trace-alignment checkpoint: the completed live IRC now has an append-only
`unexpected-endpoint-discovery-003.json` with validated endpoint correspondence
for candidate `reaction-2cdf1d1a16c24c3a`. Only original hydrogens 5 and 6 are
permuted within their same-parent, same-isotope group; hydrogen 5 is identified
as the observed transfer. Heavy-atom maps and all retained geometries remain
unchanged. The original radical-homolysis product check still fails. Regressions
also reject isotope-changing reassignment and malformed/nonfinite trajectories.
This closes the atom-correspondence prerequisite above, not minimum, state,
stereochemical, path, or rate validation. Next: stationary fragment/state checks,
bounded queue handoff, and reuse of the retained TS under explicit evidence
contracts. No user-supplied pathway import or new command-line flags are needed.
Validation at this checkpoint: full unittest discovery passed 632 tests in
66.221 seconds; the saved 119-frame live trace replay produced the aligned
candidate without launching new electronic-structure calculations. These are
software and retained-evidence checks, not an experimental kinetics benchmark.

The graph/xTB v2 implementation now performs condition-specific bounded
radical association and H-atom abstraction against explicitly declared gas
partners with retained charge/spin contracts. Closed-shell partner chemistry
outside those classes, encounter/capture kinetics, and solvent-assisted paths
remain explicit coverage gaps rather than inferred support.

Schema-v3 discovery slice (2026-10-06; shared with item 9/M4): a named implicit
solution request retains its intended SMD/CPCM path conditions separately from
observed native GFN2-ALPB/`gsolv` endpoint screening. A narrow mapped,
valence-preserving isolated-carbon cyclobutene opening/acyclic diene closure
edit adds neutral candidates; it does not establish orbital/stereochemical
allowedness, a TS or broader pericyclic coverage. Unsupported environments
cannot fall back to gas. Native convergence/connectivity and method/solvent/
invocation/JSON/stdout/source checks govern reuse; endpoint and assembly
contracts bind environment and refuse mismatch rather than rewriting old gas
artifacts. The shared path gate's `phase="solution"` alias repair preserves
required named-solvent, SMD/CPCM and exact method-keyword checks; it corrects
vocabulary without granting new solvent physics. Three actual water/acetonitrile
development cases recovered three
expected products using six one-core xTB endpoint jobs in 0.924418 seconds;
cache checks preserved 24 files per case with zero additional endpoint jobs.
No ORCA path, rate, lifetime, coverage or calibration is demonstrated. See
[screening evidence](runs/solution-discovery-screen-diagnostic-001/EVIDENCE_2026-10-06.md).
M2 and item 9 remain in progress, with held-out discovery benefit and complete
solution prediction still required; the many-environment roadmap is unchanged.

Scientific hypothesis: independent discovery disagreement identifies important
mechanism-coverage risk that cannot be exposed by improving RMG barriers alone.

Deliverables:

- automatic RMG plus graph/xTB execution under one bounded discovery contract;
- user and curated pathway import through the same native representation;
- agreement, disagreement, truncation, and chemical-class coverage evidence;
- candidate-network construction independent of RMG rates;
- environment-specific rules for reusing connectivity hypotheses; and
- selective physics-driven discovery for influential uncovered regions.

Exit gate: preregistered RMG-success and RMG-blind-spot cases demonstrate
improved candidate recall without unacceptable high-level verification cost;
truncated searches remain unresolved rather than complete.

Do not build yet: exhaustive arbitrary bond-edit enumeration or an ML model
presented as a completeness oracle.

### STAB-011 — Bounded reaction-surface resolution and segmented verification

**Status: implemented_unvalidated.** This is a dependency-restoration track
within M2/M5/M6 rather than a new environment lane. It was promoted ahead of
STAB-009 because the retained tert-butoxy calculation exposed a concrete path
resolution failure that blocked trustworthy kinetic verification.

Scientific hypothesis: quantitative diagnosis plus a bounded sequence of
stationary-point recovery, mapped coordinate scans, and elementary
segmentation resolves more chemically meaningful paths than repeating
equivalent full-endpoint NEBs, while decreasing false TS and false-rate claims.

Implemented software scope:

- immutable `generic-orca-path-v6-dissociation-loss-bound` policy, including a
  conservative neutral-gas 1-to-2 detailed-balance loss ceiling, and explicit
  diagnoses/terminal states;
- append-only recovery of normally terminated but unconverged stationary
  optimizations, including Hessian, coordinate-system, and trust-radius
  escalation;
- bounded mapped one/two-coordinate scans whose maxima remain TS seeds only;
- shared OptTS, meaningful-frequency, intended-mode, and bidirectional IRC
  gates;
- geometry-derived intermediate identity with charge/spin/distinctness and
  common atom-map validation;
- two STORCA-native elementary candidates with parent-route provenance;
- independently verified segment paths, recursion prevention, and specialist
  boundaries for separated/barrierless segments;
- elementary Arkane-rate bundles and all-at-once repaired-network
  repropagation, with no effective one-barrier collapse;
- route-wide calculation/core-hour/wall-time accounting from ORCA process
  leases, append-only snapshots, and timeout contraction at the remaining
  budget;
- one common append-only evidence envelope for every standard NEB orientation,
  stationary recovery, guided scan, and segment attempt, linking scientific
  context, exact inputs, leases, budget, results, failures, and claim bounds;
- per-image connectivity and mapped separated-fragment motion diagnostics;
- connectivity-first rejection of fragmented intermediate geometries before
  bond-order assignment;
- pre-ORCA environment applicability gates, including exact CPCM/SMD solvent
  contracts for implicit solution and explicit escalation for unsupported
  phases; and
- precise `segment_path_unresolved`, `segment_rate_model_unavailable`,
  `verified_multistep_path_network_insertion_pending`, barrierless,
  electronic-state, and budget-exhaustion abstentions.

Exit gate: focused contract and negative tests pass; retained profiles cover an
ordinary saddle, multiwell path, barrierless path, shallow mode, and
tert-butoxy beta scission; tert-butoxy stops equivalent full NEB repetition and
reaches either verified elementary rates plus repaired propagation or a precise
evidence-backed specialist abstention. Synthetic fixtures prove software
behavior only and do not close this gate.

Do not build yet: universal capture theory, custom RRKM/master-equation code,
unbounded recursive segmentation, AFIR/dynamics, or an effective multistep
barrier used merely to force a lifetime.

## Milestone 3: Conformer- and temperature-aware rate evidence

**Status: in progress.** Explicit single-temperature production and observed
ORCA free-energy recipe matching and a bounded reactant sampler → native ORCA
refinement → rate-admission diagnostic are implemented and software-tested.
Distinct minima without interconversion/TS-conformer evidence withhold aggregate
rates; conditional scenarios remain diagnostic only. Actual ORCA conformer-rate
benchmarking, TS/product/encounter ensembles, interconversion, unsampled-risk
calibration and the scientific exit gate remain open. Retained installed-RDKit
seeds and read-only gas endpoint cache checks do not close this milestone.

Scientific hypothesis: bounded conformer ensembles and condition-matched
thermochemistry materially improve rate credibility and expose uncertainty that
one endpoint/path conformer at 298.15 K cannot represent.

Deliverables:

- CREST/xTB conformer acquisition for reactants, products, encounter complexes,
  and transition states where applicable;
- ensemble free energies and conformer-weighted rates;
- arbitrary-temperature ORCA thermochemistry contracts;
- conformer-domain coverage and missing-conformer uncertainty;
- bounded TS-conformer and endpoint-orientation allocation; and
- cache identities including environment, temperature, charge, spin, and
  conformer protocol.

Exit gate: conformational benchmark cases improve or retain rate accuracy,
repeatable ensembles produce compatible rates, and inadequate conformer
coverage blocks precision claims.

Do not build yet: unconstrained global conformational exhaustiveness claims.

## Milestone 4: Bimolecular homogeneous-solution kinetics

**Status: in progress.**

Scientific hypothesis: a bounded implicit-solution lane with explicit
concentration, standard-state, reversibility, and diffusion diagnostics can
support useful observables for a restricted class of homogeneous reactions.

Deliverables:

- first- and second-order dimensional rate contracts;
- reversible nonlinear concentration propagation;
- target reformation, branching, product yield, and concentration-at-time
  observables;
- ideal-activity declaration and optional supplied activity corrections;
- diffusion-limit and encounter-complex diagnostics;
- fixed user-supplied speciation and pH contracts; and
- charge, element, nonnegativity, flux, and numerical closure.

Implemented slice: `solution-neutral-implicit-mass-action-v1` is a separate
verified-input lane so the original linear first-order artifact contract is not
weakened. It admits dimensional first- and second-order rates, requires an
evidenced diffusion ceiling for every bimolecular rate, propagates reversible
target consumption/reformation, switches audibly to a linearly implicit method
when an initial-Jacobian stiffness diagnostic crosses its threshold, and retains
concentration, product, extent, uncertainty-scenario, and closure evidence.
Standardized requests now cover t95, retention, concentration, product yield,
dominant product, gross target-consumption branching, and integrated pathway
extent, including several claims from one immutable lane request. Half-life,
selectivity, and induction time remain open. Supplied nonideal activities,
fixed pH/speciation contracts, experimental calibration, and autonomous
solution discovery remain open, so this milestone has not passed its exit gate.

In-progress producer slice (2026-10-05): the existing native solution path
workflow can emit conditional forward/reverse selected-saddle coefficients
from an internally selected neutral one-/two-molecule candidate. A
solvent-bound native adapter needs no preassembled species table, concentration
history, duration or invented bound. Separated-species free-energy sums,
explicit 1 M conversion, units and detailed balance are checked; the existing
CLI adds no new flags. Coefficients have unresolved uncertainty/ensembles/
transport and no network-admission authority or stability bundle. Synthetic
regressions are not real solution validation or autonomous discovery. The
completed 999-test application/68-test native checkpoint covers this initial
software slice only. Applicable
error evidence, diffusion and concentration-dependent experimental production
benchmarks remain exit requirements; the milestone stays `in_progress`.

The 2026-10-06 solvent-bound graph/ALPB screening slice described in M2 now
supplies native neutral connectivity hypotheses for internally selected
solution candidates. Its six actual semiempirical endpoint jobs are discovery
evidence only; they do not complete autonomous path/rate/propagation, quantify
coverage/error or promote the solution lane's scientific applicability.

Exit gate: analytical fixtures and preregistered experimental solution systems
pass declared rate, branching, conservation, and uncertainty-calibration gates.

Do not build yet: automatic universal pKa prediction, implicit conversion of
gas rates, or claims of explicit-solvent mechanism coverage.

## Milestone 5: Auditable method selection and escalation

**Status: planned.**

Scientific hypothesis: explicit diagnostics can prevent influential reactions
from being silently assigned an inapplicable default method and can route them
to a bounded specialist treatment or unresolved state.

Deliverables:

- versioned method-applicability findings and escalation policies;
- diagnostics for barrierless association, tunneling, pressure dependence,
  conformational sensitivity, solvent participation, spin changes,
  multireference risk, excited states, surfaces, and applied potential;
- separation of automatic diagnostics, user declarations, and specialist
  external evidence; and
- complete decision provenance with no LLM authority over physical selection.

Exit gate: a challenge set of ordinary and pathological reactions is routed to
the preregistered treatment or an explicit unresolved state with calibrated
false-escalation and missed-escalation reporting.

Do not build yet: a black-box universal method recommender.

## Milestone 6: Pressure-dependent, barrierless, and tunneling kinetics

**Status: in_progress.** Arkane pressure-dependent job preparation and
barrierless classification safeguards exist, but no complete generally
supported escalation chain is available.

The source-neutral specialist-rate request/result contract, immutable request
export, hash-matched retained-result import, rate-only checkpoint continuation,
and repaired-network handoff are implemented. Imported library files are
checksummed and a specialist continuation reuses retained path evidence instead
of rerunning ORCA. A default calibrated VTST/capture/master-equation producer
and benchmarked uncertainty models remain outstanding.

The M6 development benchmark and provider conformance harness are now
implemented. It measures log-rate error, uncertainty coverage, missing
predictions, and false method applicability against condition- and unit-specific
reference points. It now also represents exact-pressure falloff curves, scores
curve log-rate RMSE and pressure span, verifies cryptographic holdout locks and
prediction-freeze declarations, and retains source-backed cases whose numerical
reference data are still explicitly pending. The retained manifest contains no
locked holdout or completed falloff curve and therefore cannot promote M6 beyond
`in_progress`. A retained JPL 2011 visual/source audit now verifies the four
development reactions' paired limiting-rate recommendations, units, air bath
and evaluated envelope semantics. Applying the published f(T) correction to
300 K changes reference factors only; the original manifest bytes, central
rates, conditions, partitions and acceptance thresholds are retained. These
eight anchors are not independent experimental samples or calibrated STORCA
error distributions. The backend decision record separates a
VaReCoF-class capture producer from a MESS-class master-equation propagator.
The first MESS adapter slice deterministically renders an explicit
one-well/two-fragment model containing source-backed RRHO, collisional, and
external Rotd microcanonical-flux inputs, and parses one exact native T/P rate
table. A second slice now executes that prepared input with a hard timeout,
revalidates immutable evidence, hashes the actual executable binary, checks the
exact requested channel and T/P grid, and reuses a valid retained checkpoint.
The third slice implements a separate scientific admission gate. It requires
binary-bound amount-basis evidence, complete component uncertainty, exact-grid
conditions, full-grid collision ceilings, and independently sourced equilibrium
rate ratios; parses the opposite channel; and checks detailed balance at every
T/P point. An admitted table is fitted to an RMG pressure-dependent model,
reloaded, checked over the full grid, and packaged through the existing
specialist-result and repaired-network continuation. Live scientific falloff
benchmarks, a locked holdout, and calibrated uncertainty remain fail-closed
follow-up gates.

Calculated Arkane central rates are no longer treated as exact during
prediction convergence. Nonzero-flux calculated replacements require explicit
condition-matched rate uncertainty; verified and unresolved rate bounds are
propagated together in one simultaneous prediction envelope. Missing
calibration leaves the observable unresolved rather than assigning a default
factor.

The 1↔2 Arkane path also accepts a route-specific, hash-backed collisional
model for an exact energized well and bath gas. Explicit
`SingleExponentialDown` and Lennard-Jones parameters are written into the
Arkane input; absent evidence remains visibly estimated and unvalidated. This
is a narrow supported input contract, not calibrated parameter prediction.

Scientific hypothesis: specialist rate models for influential reactions change
macroscopic predictions in known regimes and reduce systematic misuse of
ordinary harmonic TST.

Deliverables:

- capture or VTST adapters for loose association;
- RRKM/master-equation execution and validation for pressure dependence;
- tunneling sensitivity diagnostics and bounded corrections;
- explicit bath-gas, energy-transfer, pressure, and temperature evidence; and
- uncertainty propagation from these rate theories.

Exit gate: benchmark rates and falloff curves meet predeclared tolerances and
ordinary TST is rejected in challenge cases where it is known to fail.

Do not build yet: a custom master-equation package when Arkane can be adapted.

## Milestone 7: Microsolvation and solvent-assisted escalation

**Status: planned.** The lane is registered as declared-only.

Scientific hypothesis: selectively adding explicit solvent molecules for
observable-influential reactions can resolve specific solvent-assisted paths
without treating finite clusters as bulk solution.

Deliverables:

- bounded microsolvated cluster generation and solvent stoichiometry;
- cluster/continuum path verification and conformer coverage;
- standard-state and solvent-chemical-potential accounting;
- cluster-size and solvent-placement convergence diagnostics; and
- escalation back to unresolved when bulk free-energy physics is required.

Exit gate: curated proton-shuttle and solvent-assisted benchmarks demonstrate
mechanism and rate improvement over continuum-only calculations without false
bulk claims.

Do not build yet: labeling a few explicit solvent molecules as an equilibrated
liquid.

## Milestone 8: Photochemical stability

**Status: in_progress.** Condition-bound supplied photolysis evidence can enter
the gas ladder, and TD-DFT screening exists; autonomous excited-state reaction
discovery, quantum yields, and nonadiabatic kinetics do not.

Scientific hypothesis: wavelength-resolved absorption, exposure, quantum
yield, and competing thermal kinetics can support conditional photochemical
fate predictions for a restricted domain.

Deliverables:

- immutable light spectra and exposure geometry;
- absorption and quantum-yield evidence contracts;
- photoproduct candidate discovery and branching;
- ground/excited-state applicability diagnostics; and
- propagation of photolysis with thermal networks.

Exit gate: known photolysis systems reproduce rates and dominant products under
matched spectra, with missing quantum-yield evidence failing closed.

Do not build yet: infer quantum yields from oscillator strengths or treat TD-DFT
excitation energies as photochemical rates.

## Milestone 9: Interfaces and heterogeneous exposure

**Status: planned.**

Scientific hypothesis: explicit phase-contact and transport contracts can
separate reaction-controlled from transport-controlled fate without hiding the
interface inside a homogeneous surrogate.

Deliverables:

- phase, geometry, area, volume, partitioning, and mass-transfer contracts;
- gas-liquid and liquid-solid coupling;
- adsorption/desorption evidence and surface-site bookkeeping; and
- reaction/transport sensitivity and convergence diagnostics.

Exit gate: benchmark reactors with known reaction- versus transport-limited
regimes are classified and propagated correctly.

Do not build yet: general CFD or arbitrary reactor engineering.

## Milestone 10: Periodic surface and electrochemical specialist lanes

**Status: deferred.** Both lanes are registered as declared-only.

Scientific hypothesis: specialist periodic and potential-aware evidence can be
integrated through common claims and networks without making molecular ORCA
jobs the architecture.

Deliverables:

- slab, adsorption-site, coverage, surface-charge, and potential contracts;
- periodic path verification adapters;
- microkinetic propagation and site balance;
- constant-potential or explicitly declared approximation evidence; and
- surface/electrochemical benchmark profiles.

Exit gate: one narrowly defined surface system and one electrochemical system
pass mechanism, rate/selectivity, balance, and uncertainty gates against
independent references.

Do not build yet: broad periodic-code support or a universal electrochemistry
workflow.

## Milestone 11: Prediction-conditioned validation and calibration

**Status: planned.** Software regressions and small discovery fixtures exist;
a preregistered cross-environment scientific benchmark does not.

Scientific hypothesis: sensitivity-driven and multi-discovery STORCA improves
observable accuracy, mechanism recall, and uncertainty calibration per unit
compute over RMG-only and current fixed-allocation baselines.

Deliverables:

- locked train, development, and holdout systems;
- gas, solution, competing-path, radical, pressure-dependent, barrierless,
  conformational, and difficult electronic-structure cases;
- mechanism recall/precision, dominant-path accuracy, rates, branching,
  lifetime, calibration, wall time, high-level job count, and information per
  compute metrics; and
- comparisons among RMG only, current STORCA, sensitivity-driven STORCA, and
  multi-discovery STORCA.

Exit gate: preregistered primary endpoints improve without degraded calibration
or hidden exclusion of failed systems. Fail-closed abstention is scored rather
than discarded.

Do not build yet: a large benchmark whose conditions or mechanisms are too
poorly curated to support the claimed comparison.

## Execution order

The dependency order is:

```text
M0 traceability
  -> M1 campaign foundation
  -> M2 autonomous discovery
  -> M3 conformer/temperature credibility
  -> M4 bimolecular solution
  -> M5 method selection
  -> M6/M7/M8 specialist kinetic and environment escalations
  -> M9/M10 heterogeneous specialist lanes

M11 validation begins now and expands alongside every milestone.
```

No milestone may close without updating the requirement catalog, traceability
matrix, machine-readable manifest, tests, benchmark evidence, `METHODS.md`, and
`LIMITATIONS.md` where the scientific claim or applicability boundary changes.

P1/P2 execution checkpoint (2026-09-09): retained IRC fragment seeds now feed
the existing stationary-species preparation/validation helper in declared atom
order; seed reuse no longer overwrites species labels with element labels.
The full software suite passes 632 tests (68.597 s). Automatic unexpected-route
queue admission, validated retained-TS reuse, and end-to-end kinetics remain
open. P2 curation is recorded in `benchmarks/stability/scientific_panel.json`:
eight development entries, including four aqueous saponification measurements,
with explicit missing conditions/uncertainty and no frozen holdout or support
claim. Charged aqueous chemistry is not admitted by the neutral-solution lane.

Unexpected-channel queue checkpoint: engine callback proposals now enter a
persisted, one-generation bounded supplemental queue after identity, reactant,
environment, mapping and geometry-hash checks. Duplicates do not launch new
work and no rate is copied. A retained-coordinate regression selects the parent
then its unexpected channel with rate callbacks prohibited (59 focused tests
pass). A separate replay against the actual NNONC artifacts admits
`reaction-2cdf1d1a16c24c3a` as unverified; the append-only
`unexpected-queue-admission-diagnostic-001.json` records this without running
chemistry or modifying the old engine state. Retained-TS admission, continuation
of that new evidence, and eventual rates remain open.

2026-09-10 rate-execution hardening: Arkane and its RMG syntax/library helpers
now have finite per-command timeouts with POSIX descendant cleanup and retained
failure streams. Low-pressure gas condition spelling is recognized by the
conditional thermochemical screening bound; non-gas models remain excluded.
172 focused tests passed, then 11 execution/thermochemistry tests passed after
adding route timeout rejection coverage. This adds no validated lifetime or
environment support. Live rate generation and repaired-network propagation
remain required; see the comprehensive plan for retained-path checkpoints.
Subsequently all 639 tests passed in 55.456 seconds. The real new-channel input
loaded successfully in `rmg4_env` (RMG 4.0.0), including the pressure-dependence
job. This verifies API/environment compatibility only, not computed kinetics.

Later live checkpoint: an atom-balanced electronic-energy reference correction
allowed the retained new-channel Arkane calculation to finish pressure-dependent
kinetics in `rate-execution-003`. Full suite passed 639 tests after the gauge
change. Export then failed on Arkane's missing native species dictionary; a
bounded graph-serialization fallback is under test. Preserve all three rate
attempts. Completed kinetics are not yet an admitted rate or lifetime; physical
applicability, uncertainty, and repaired-network gates remain.
The subsequent validation-only pass successfully exported and reloaded the
new-channel Chebyshev library, reproducing its condition rate. Eleven focused
Arkane tests passed after the dictionary fallback. Uncalibrated rate uncertainty
and collision parameters remain explicit: no converged lifetime is established.
Full unittest discovery subsequently passed 640 tests in 63.288 seconds after
the export fix. Live rate/library artifacts are retained; network propagation
and scientifically calibrated prediction remain outstanding.

New-channel repair diagnostic (2026-09-10): RMG loaded the verified library into
its edge, but the retained core did not contain that channel and generation
encountered DASPK initial-step singularities. Repair validation rejected it;
there is no independent propagation/flux closure or lifetime. The bounded
attempt is retained under `network-repair-001`. Thirty-seven engine tests pass,
including new discovery -> insertion -> propagation with missing uncertainty
still blocking convergence. The scientific panel now has nine development
entries, including sourced hydrocarbon-solution dimerization rates; it is not
frozen or a claim of supported mixture predictions.

2026-09-11 repaired-core hardening (P1; STAB-PROP-001): validated generated
repair libraries now enter native RMG core seeds; reference libraries remain
on the edge. Name collisions and changed manifests are rejected. The retained
`network-repair-003` confirms initial core insertion of the new channel, but
its first attempt reports twelve core collision-rate violations and the next
attempt failed when RMG cleared its external-library map on restart. A
process-local identity-bound alias adapter now passes a read-only replay with
the real saved libraries and 64 focused tests. Full-suite verification before
that final restart fix passed 646 tests. No repaired-network convergence,
lifetime, or expanded environment support is claimed. Consult the comprehensive
plan for the append-only sequence rather than rerunning completed path/rate work.

Subsequent serialization audit: the retained new-channel Chebyshev rate
reloaded as NaN because RMG's Chemkin writer rounded its minimum pressure to
zero. The launcher now preserves that polynomial block's bounds and coefficients
at binary64 precision, while the bridge rejects nonfinite forward rates and
nonstandard JSON numbers. A separate derived export reproduces 30 source T/P
rates exactly at observed output precision without changing any original file
or running chemistry. Four native RMG roundtrip tests pass, covering gas rate
dimensions and two coefficient shapes. The next handoff gap is source-bound,
resonance-aware reaction matching after RMG species renaming; calibration and
complete propagation remain open. No support claim changes.

2026-09-22 identity-handoff hardening (P1; STAB-PROP-001): source-bound,
resonance-aware replacement correspondence now supplies the same hash-backed
proof to applied-rate verification and verification planning. Reversibility,
charge/spin/isotopes, repeated species, exact library provenance and unique
matching are enforced. An irreversible-library fixture exposed a missing
`reversible=False` declaration also present in the production photolysis writer;
both are corrected, with native loading/rate tests and specific loader-failure
diagnostics. No installed RMG code or retained source artifacts were modified.

The append-only `identity-replay-001/result.json` in the retained NNONC handoff
passes all eleven fail-closed replay checks: the graph candidate is found but
potential stereochemistry remains unresolved, both consumers refuse verified
credit, reopened evidence agrees, and all ten original artifact hashes match.
This is evidence-handling validation, not network convergence or a lifetime.
See the comprehensive plan's dated checkpoint for final test counts and next
work. No lane or milestone is promoted; all environment goals remain intact.

2026-09-29 P1 attempt-evidence checkpoint: failed restarts no longer erase
earlier parsed discoveries. Immutable condition/source-bound attempt records
and a structural-only recovery queue preserve hypotheses without borrowing
earlier kinetic or convergence credit. Raw parser output survives inspection
failure but cannot authorize route verification without endpoint evidence.
The saved native 93-reaction replay passes 14 checks with 45 original files
unchanged. Application suite: 727 tests, 18 expected skips; focused native
suite: 23 passed. These are engineering gates, not a stability prediction.

The focused working-pipeline goal and its remaining acceptance requirements are
tracked in `docs/stability/WORKING_PIPELINE_GATE.md`. A cyclobutene development
case is preregistered before live predictions, with unchanged uncertainty and
coverage gates. A real quantitative result and controlled interruption/resume
are still required. No environment milestone is removed or promoted.

2026-10-07 prediction-first checkpoint: native008 has real completed OptTS and
frequency evidence but timed-out, forward-only IRC and no rate/lifetime. An
explicit evidence-bound derived preparation now reuses those completed stages
without copying leases or treating the partial trace as connected chemistry.
Original conditions, budgets and later matter/exposure milestones remain intact.
Independent-queue pure-target gas alias hardening and an unaudited primary
reference lead also advance; neither expands scientific support. Consult the
working-pipeline and reuse ledgers for measured evidence and remaining gates.

Further 2026-10-07 evidence: a separate primary measured OH + methane reference
slice is source-bound but still has a missing prediction; valid failure-cost
proof now survives malformed prediction metadata. A native partition audit
exposes TS optical-counting disagreement without treating it as a justified
rate correction. These advance reference/admission work, not P1, calibrated
uncertainty or environment support. The shared panel, first-ten acceptance
requirements, later matter/exposure milestones and compute ceilings remain
unchanged; see the dated working-pipeline evidence.

Item-5 allocation hardening additionally requires explicit exact-condition
binding for method calibration and prevents malformed IDs being stringified
into matches. It fixes demonstrated public-consumer admission defects, not the
missing automatic calibrated error producer; no scientific gate is promoted.
