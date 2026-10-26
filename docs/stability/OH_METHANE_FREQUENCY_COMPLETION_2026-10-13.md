# OH/methane frequency-only completion

Status: **fresh +15° NumFreq completed normally; all three constrained pilot
points and the unchanged local coupling checks passed**.
Run: `runs/oh-methane-angular-frequency-20261013-ivb3p6ba`.
[Contract](../../runs/oh-methane-angular-frequency-20261013-ivb3p6ba/contract.json).
[Result](../../runs/oh-methane-angular-frequency-20261013-ivb3p6ba/result.json).

The prior matched-runtime optimization completed normally, satisfying all five
candidate criteria. Its failed MPI-evidence receipt is retained unchanged.
The new [pure loader reader](../../storca/mpi_runtime_evidence.py) stops path
capture before the next intact `dyld` header, including zero-whitespace joins.
It does not reconstruct damaged records or lower the missing/foreign-library
gate. It requires at least four identified gradient workers and a private MPI
library witness at every identified gradient-worker PID. Its scope excludes
unobserved ranks, all-runtime coverage and proof of the original crash cause.
The caller separately binds the actual library files, runtime, native output,
process ownership and completion.

Read-only replay recovers 140 gradient PIDs/719 complete private-MPI records in
the completed optimization. Candidate geometry, electronic state and both
angular targets pass. Center/−15° frequency receipts replay without new jobs.
Their local orthogonal ZPE shift is −60.95 J/mol (limit ±500 J/mol) and maximum
absolute log-partition shift is 0.01229 (limit log(1.20)=0.18232). These are
local consistency checks, not global separability or rate evidence.

The new [driver](../../tests/complete_oh_methane_angular_frequency.py) creates
exclusive source-keyed claims, never overwrites old evidence, and permits
**zero new saddles and one NumFreq**. Cumulative prior counts six saddles/two
NumFreq become six/three after this dispatch; the two earlier execution EnGrad
diagnostics remain separate. Four cores, six hours/job, 1000 MB MaxCore/rank,
the 4 GiB disk reserve, private MPI4.1.1 and no automatic retries are retained.

The fresh unconstrained central Hessian/gradient must pass the original
constrained Lagrangian-Hessian, stationarity, state, mapped H-transfer and
orthogonal-coupling gates. A valid individual point with failed coupling remains
useful evidence, but does not admit a separable angular correction.

Twenty parser tests and thirteen completion tests pass; seventeen preceding
restart tests also reran successfully. Tests include glued loader records,
missing/foreign-library refusal, completed-optimization reuse, frequency-only
dispatch, immutable source checks, exclusive claims, failure accounting,
unchanged physical thresholds and all promotion flags remaining false.
This is not a full-suite or scientific accuracy certification.

Even a passing result is only `pilot_passed_profile_pending`: no periodic
coverage, global metric, mirror/basin count, angular partition or new rate is
admitted. Under the original global caps, eleven saddle attempts and one
selected numerical Hessian remain after this job; a future profile must not
reset those counts or assume eight refinement attempts are still available.

## Native result and independent closeout

The result is `pilot_passed_profile_pending`, not an admitted rate. Exactly one
NumFreq and zero optimization jobs were dispatched by this continuation;
cumulative counts are six saddles/three NumFreq. Raw native completion, exact
mapped input/geometry, Hessian/gradient artifacts, declared electronic state,
runtime witnesses and all old hashes/mtimes were independently reopened.
The frozen physical validator reproduced all three stored receipts, and the
orthogonal comparison replayed exactly. No old artifact or source was changed.

| Offset | Tangent gradient RMS / max, Eh/bohr | Retained transfer overlap |
| --- | --- | --- |
| 0° | 5.0598e−7 / 9.1970e−7 | 0.99999968 |
| −15° | 2.3414e−6 / 5.4230e−6 | 0.99959253 |
| +15° | 3.2715e−6 / 8.6377e−6 | 0.99957846 |

Each constrained 14-dimensional tangent has one negative/thirteen positive
curvatures. The +15° orthogonal ZPE shift is −18.3855 J/mol and its maximum
absolute log-partition shift is 0.00252192, both inside the unchanged
500 J/mol/log(1.20) coupling limits. These are local diagnostics only.

**Do not label +15° an ordinary first-order TS.** Its full numerical Hessian
contains −1239.47 and −38.24 cm⁻¹ modes, while the constrained Lagrangian
tangent has only the intended unstable H-transfer direction. The second full
negative mode has mass-weighted angular-normal overlap 0.817893; the transfer
mode's overlap is 0.0005991. With nonzero angular force multiplier
7.6768873e−6 Eh/rad, this supports an angular-like instability at a constrained
point, not a proved pure rotor or angular barrier. It cannot be inserted as a
new independent harmonic TS basin or used to select a factor of two.

The native gradient runtime proof contains eight identified gradient PIDs and
forty complete private-MPI library records. macOS sleep/dark-wake interruptions
confound elapsed-time performance; the result is not a speed benchmark.

The next scientific stage is a separately frozen bounded periodic profile:
coverage/closure, mapped equivalences, counting, kinetic measure, orthogonal
contributions and harmonic-limit recovery must be explicit before any angular
partition or revised rate. The wider profile and rate reevaluation have not
been launched or implemented by this frequency-only completion.
