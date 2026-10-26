# OH + methane periodic profile: implementation checkpoint

Date: 2026-10-15. Status: **software implemented and offline-tested; bounded
native profile running in kitty**.

Follow-up (2026-10-16): this native run subsequently completed +60° normally
but stopped at its extra displacement screen. The
[append-only correction](OH_METHANE_CANDIDATE_POLICY_CORRECTION_2026-10-16.md)
reuses the evaluated frame with a corrected native-gradient identity and
preserves this historical run. The launch status below is a dated checkpoint,
not a claim that the old process is still running.

Implements the diagnostic branch of
[the bounded profile plan](OH_METHANE_PERIODIC_PROFILE_PLAN_2026-10-13.md).
The full angular thermal/rate model remains unadmitted. This does not promote
production defaults, environment support, network reinjection or a lifetime.

## What is implemented

- `storca/ts_angular_profile.py`: mapped proper-alignment/closure diagnostics,
  spectator permutations without automatic symmetry factors, the existing
  constraint-tangent gradient convention, discrete periodic coverage, explicit
  candidate-versus-Hessian evidence states, deterministic provisional
  final-Hessian ranking, and an aligned relaxed-path metric diagnostic.
- `storca/angular_attempts.py`: a dedicated append-only hash-chained journal,
  exclusive/fsynced reservations, global 17/4 caps including prior 6/3, stable
  controller locking, exact immutable execution/analysis receipts and external
  history bindings. Incomplete records block spending rather than refund it.
- `tests/run_oh_methane_angular_profile.py`: one source/resume argument;
  read-only physical replay of all three completed pilot points; two scan
  arms through ±60°, ±120°, ±170°; independent +170→+190° closure; and one
  selected full NumFreq. Owned normal native results survive wrapper/reader
  interruptions without redispatch. Partial wrapper receipts remain preserved
  behind a separately bound adoption receipt. Final derived reports are
  published atomically without overwriting prior evidence.
- `tests/replay_oh_methane_angular_rates.py`: source-bound offline rate-admission
  assessment. No supplied success flag or arbitrary Q(T) table can turn on a
  nonexistent physical model. Earlier conditional recipes remain source
  evidence, not recalculated or recertified rates.

The ordinary `CompositeStationaryPoint` validator and all historical pilot,
runtime, route and rate artifacts/source versions remain unchanged. The new
owned-job candidate adapter reproduces the frozen candidate reader on all
three actual completed pilot optimizations. It changes execution adoption,
not the native convergence or physical saddle gates.

## Validation completed before new chemistry

**241 focused tests pass**, including 100 new tests across the four new test
files, prior OH/methane drivers, angular mathematics, MPI evidence, composite
rate validation and roadmap traceability. These are software/diagnostic tests,
not a physical accuracy certificate.

The actual `oh-methane-angular-frequency-20261013-ivb3p6ba` predecessor chain
replays with ORCA dispatch forbidden. All 0°/−15°/+15° physical receipts and
local coupling checks pass unchanged. New/frozen candidate adapters agree.
Prior source hashes and mtimes, runtime identities and native ownership are
reopened rather than trusting a saved success flag.

Tests cover concurrent controllers, partial journal writes, exhausted budgets,
the unique remaining Hessian, normal native completion followed by a wrapper
failure, missing/truncated execution receipts, analysis interruption, sealed
report publication interruption, and exact closed-history replay. Ambiguous
jobs cannot be relaunched and no process is killed on PID alone.

## Native campaign contract

Active continuation:
`runs/oh-methane-periodic-20261015-ei3ezay9`.
The +60° optimization launched normally on four witnessed native worker
processes and reached cycle 2 at the launch check. Its owned lease is running;
no new saddle/partition/rate result is yet admitted. The cumulative journal
is seven optimization attempts/three Hessians after dispatch. All 554 bound
source/runtime records reopen unchanged, with approximately 52 GB disk free.
The source claim and `preflight.json` retain the exact tested starting evidence.

Reuse 0°/±15° and their three Hessians. Six wider optimizer candidates plus
one independent closure attempt are planned; then select at most one Hessian.
All attempts, including failures, consume the remaining 11 optimizations and
one Hessian. Four cores, MaxCore 1000 MB/rank, six hours/job, private witnessed
MPI 4.1.1 and the 4 GiB disk reserve remain unchanged. No automatic retries.

The final Hessian checks one additional point only. Unsampled spectra do not
become verified by interpolation. Known chart/physical applicability failures
close `profile_unresolved`; parser/artifact/ambiguous execution problems retain
an attention checkpoint and their original native work.

## Why this does not yet produce a corrected rate

No reviewed angular/reaction/rotation Hamiltonian and thermal measure is
implemented. An inverse-constraint metric is not a relaxed-path inertia, and
the latter diagnostic alone is still not a complete partition measure.
Sparse Hessians do not certify the transverse bath throughout the circle.
Symmetry/channel/mirror counts remain jointly unresolved.

Consequently this campaign explicitly selects **diagnostic-only** output.
After a complete valid coarse profile, the result is
`profile_complete_partition_unsupported`; failed gates give
`profile_unresolved`. Both retain precise missing evidence and the earlier
conditional rate sources. There is no new angular rate multiplier, experimental
fitting, extra imaginary-mode waiver or automatic refinement of an undefined
partition target. Unused refinement budget is preserved rather than spent.

The larger plan's conditional rate overlay, harmonic/free-rotor partition
limits, energy-zero/ZPE replacement and updated forward/reverse tunneling
remain blocked on a defensible reduced model. This checkpoint does not claim
those optional quantitative deliverables were implemented.

## Run and resume

From the repository in the existing `orca-auto-py312` environment, with the
same bound ORCA executable configured as before:

```sh
python -B tests/run_oh_methane_angular_profile.py runs/oh-methane-angular-frequency-20261013-ivb3p6ba
```

The source claim directs repeat invocations to the same continuation. Its
printed work directory can also be supplied directly as the sole argument.
Neither route resets the journal. The controller stops at a scientific failure
or unresolved execution rather than silently starting a replacement run.

On closure, `result.json`, `result-binding.json`, `angular-rate-assessment.json`,
the attempt journal, source snapshots, native logs/leases and checkpoint files
form the retained assessment. The offline consumer needs no RMG environment;
there is no new environment-selection flag. Future actual Arkane work must
continue to use the existing `--rmg-env` mechanism.
