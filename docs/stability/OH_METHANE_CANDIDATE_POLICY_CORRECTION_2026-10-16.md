# OH + methane: candidate-policy correction and continuation

Date: 2026-10-16. Status: **implemented; 300 focused tests pass; corrected
continuation running in kitty at +120°**.

Follow-up: +120° completed normally, then paused at the detailed MPI loader-log
reader. The [explicit runtime-reader extension](OH_METHANE_RUNTIME_RESUME_2026-10-16.md)
reuses that calculation, makes loader tracing diagnostic-only, and preserves
this frozen implementation. The running/cycle-1 descriptions below record the
earlier launch checkpoint, not the current process state.

This is a narrow correction to
[the periodic-profile implementation](OH_METHANE_PERIODIC_PROFILE_IMPLEMENTATION_2026-10-15.md),
not a new search, rate model or method sweep. The old run, its closed failure,
source snapshots, native files and attempt journal remain unchanged.

## What stopped, and what is corrected

`runs/oh-methane-periodic-20261015-ei3ezay9` completed its +60° optimization
normally after 36 cycles in approximately 16½ minutes. ORCA accepted its
overachieved-gradient convergence route. STORCA's extra preliminary displacement
screen rejected RMS/MAX steps of 2.23/2.54 times their requested tolerances
because its additional limit was two. This was an acquisition-policy refusal,
not an ORCA crash or proof of an invalid chemical path.

The new `storca/constrained_candidate_eligibility.py` preserves unique native
normal completion, final-cycle table consistency, native convergence and the
native gradient tests. It records displacement ratios and warnings without
another STORCA displacement cutoff. Candidate eligibility still requires the
independent coordinate, state, topology and gradient checks below. It never
establishes saddle character or a rate.

## A necessary correction to the earlier gradient interpretation

Matching the saved `.engrad` coordinate header to final `.xyz` did **not** prove
the gradient was evaluated at that final geometry. In this actual output:

- Property geometry **36** has the last full gradient and the dispersion-inclusive
  energy −116.25936092875381 Eh. The native cycle-36 Cartesian frame, energy and
  printed full gradient independently match these indexed property records.
- Property geometry **37** matches final `.xyz` and has a final energy-only
  evaluation, but no matching full-gradient record. The `.engrad` gradient
  components are those of geometry 36, despite its geometry/energy header
  corresponding to geometry 37.

Consequently the previous status explanation's “final-geometry gradient” label
was too strong. The numerical projection was small, but its frame identity had
not been demonstrated. The new reader fixes the identity, not the tolerance.

`storca/orca_evaluated_gradient.py` joins native property geometry, gradient and
`Single_Point_Data.FinalEnergy` by the **same GeometryIndex**, then checks the
corresponding native output frame, energy and all 21 gradient components. Bare
SCF energy without dispersion is not substituted. The last actually evaluated
frame is selected; a distinct native final XYZ remains separate evidence.
Malformed, duplicate, incomplete or mismatched sources refuse.

The reused +60° candidate therefore explicitly uses **geometry 36**, not an
unverified force assignment to geometry 37. At geometry 36, the unchanged
checks give:

| Diagnostic | Replayed value | Existing limit |
| --- | ---: | ---: |
| Tangent gradient RMS | 1.87318e−6 Eh/bohr | 8e−6 |
| Tangent gradient maximum | 3.78137e−6 Eh/bohr | 3e−5 |
| Angular error | 6.64e−8 rad | 2e−4 |

State/method, spectator topology and observed-worker MPI checks also pass.
This is a **constrained optimizer candidate**, not a Hessian-verified saddle.
No new electronic-structure calculation is required for this reuse. Original
native files, including final XYZ and the mixed-frame EnGrad file, are retained
unchanged; the new evaluated-frame receipt states which source is used.

## Append-only continuation and budget

Active correction directory:
`runs/oh-methane-periodic-candidate-20261016-ha4k1u9s`.
The +120° optimization has an owned running four-core ORCA lease and reached
cycle 1 at the launch check. Its reservation brings cumulative spending to
**8 optimization attempts / 3 NumFreq**. The child has no +60° optimization
directory: the completed source optimization was not repeated. All 572 bound
source/runtime records reopen unchanged after launch. No new saddle, partition
or rate result has yet been admitted.

`tests/continue_oh_methane_angular_profile.py` takes the closed stopped run or
its own resume directory as its sole argument. It verifies the parent result
seal, source claim, immutable contract/source inventory, exact journal history,
the single paid normal +60° execution and its absence of an old analysis
receipt. The new reuse receipt references that old reservation; it does not
create another optimization reservation.

A source-keyed exclusive claim selects one correction directory. Its contract
binds the new readers, the reused evaluated frame and the unchanged parent
execution/scientific policy. The parent journal is locked read-only while the
new controller runs. All sources are checked before the next reservation;
the retained executor checks again before native dispatch. Both journals
remain independently hash-bound. No history is rewritten or budget refunded.

Starting consumption remains **7 optimization attempts / 3 NumFreq jobs** under
the original **17/4 caps**. Thus ten optimization attempts and one Hessian remain,
including failures. The remaining scan starts **+60→+120→+170**, then follows
the negative arm and +170→+190 closure control. Six further optimizations plus
one selected Hessian are planned; unneeded refinement slots are not spent on
an unsupported partition model. Four cores, six hours/job, private MPI 4.1.1,
MaxCore 1000 MB/rank and the 4 GiB reserve stay unchanged.

## Verification and unchanged scientific boundaries

All **300 focused tests pass**: the original regression group, 20 historical
candidate-parser tests, 16 new constrained-candidate tests, 14 evaluated-frame
tests, and nine correction-controller tests. The actual retained +60° output
replays with native chemistry launch forbidden, selects index 36 rather than
37, and spends zero attempts. Synthetic controller tests prove no +60°
redispatch, no counter reset, parent-history preservation, interruption-safe
reuse of later jobs, and refusal before spending on changed evidence.

The fresh selected full Hessian must still establish the appropriate tangent
saddle character and H-transfer mode. Other optimization-only nodes remain
unverified. Closure, coupling, provenance and all no-rate/no-lifetime gates
are unchanged. A reviewed angular Hamiltonian/measure and globally supported
counting/transverse model are still absent; no corrected rate is produced.

To continue in the existing configured environment:

```sh
python -B tests/continue_oh_methane_angular_profile.py runs/oh-methane-periodic-20261015-ei3ezay9
```

The printed correction directory can be supplied directly on resume. It does
not reset counters or automatically retry a failed/ambiguous job.
