# OH/methane matched-MPI angular continuation

Status: **+15° optimization completed normally at cycle 35 with all five
criteria passing; continuation stopped before NumFreq at a log-reader defect**.
This is a continuation of the bounded three-point pilot, not a production
rate correction or completion of the full periodic angular audit.

Run: `runs/oh-methane-angular-matched-mpi-20261013-kh1a0e5v`.
[Contract](../../runs/oh-methane-angular-matched-mpi-20261013-kh1a0e5v/contract.json).
[Driver](../../tests/continue_oh_methane_matched_mpi.py).

The exact failed cycle-six coordinates passed the preceding serial/MPI
gradient diagnostic. This run restarts from those coordinates, with fresh
SCF and the original verified numerical Hessian as an optimizer estimate.
No interrupted binary optimizer checkpoint, old orbitals or failed derivatives
are admitted. The restart target differs from the original +15° target by
4.3714e−8 radians; final geometry must satisfy both targets within the unchanged
2e−4-radian tolerance. Method, mapped transfer coordinate, candidate policy
and independent physical derivative gates are unchanged.

Center and −15° frequency receipts replayed successfully without new chemistry.
One new OptTS and, only after candidate acceptance, one fresh unconstrained
central NumFreq are permitted. Cumulative upper bounds are six saddle attempts
and three NumFreq (prior five/two); the two preceding EnGrad diagnostics are
separate. Four cores, 1000 MB MaxCore/rank, six hours/job, 4 GiB disk reserve,
owned native completion and no automatic retries remain in force.

The private Open MPI 4.1.1 runtime is selected only for this process. Its
launcher, daemon, libraries, components and configuration are frozen (140
runtime files), alongside old evidence, participating code and ORCA modules.
After each native job, actual gradient-worker loader evidence must identify
the private MPI library, and inventories must reopen unchanged. No shell,
Homebrew or ORCA installation defaults are changed. An idle-sleep assertion
lasts only for this process; closing the lid can still pause the calculation.

Seventeen new focused tests pass; the 18 preceding diagnostic and 11 candidate
continuation tests were also rerun successfully. Complete loader-token
interleaving is tested; zero-whitespace glued records remain a fail-closed
reader limitation, not grounds for a false scientific admission. No frozen
driver or retained native evidence was modified during execution.

## Actual outcome and preserved reader repair

The completed optimization passed owned native completion and, when independently
reopened, its unchanged candidate gate. The original-target angular gap is
4.1933e−8 radians; the saved restart-target gap is 1.7810e−9 radians.
No MPI worker crash recurred. Cumulative counts are six saddle attempts/two
NumFreq. The optimizer Hessian is not an independently verified saddle Hessian.

The frozen reader consumed a second worker's intact loader header as part of
the preceding worker's partial framework path. It then discarded that header,
losing worker 20444's library witness. A corrected, separately versioned parser
reopens the unchanged output and finds all 140 identified gradient-worker PIDs
with 719 complete private-MPI library records. No foreign MPI path is observed.
The failed result, input, output, lease, contract and original parser remain
unchanged; there is no retrospective overwrite of its unresolved receipt.

The [frequency-only continuation](OH_METHANE_FREQUENCY_COMPLETION_2026-10-13.md)
reuses this completed geometry and repaired evidence. It permits one pending
NumFreq and no optimization dispatch, with the original physical/coupling gates.
That continuation has now completed normally and passed the local pilot.
The +15° point is valid under its angular constraint, not an ordinary
unconstrained first-order TS; the frequency checkpoint records that distinction.

Success can establish only `pilot_passed_profile_pending`, after the fresh
constrained Lagrangian-Hessian/transfer/stationarity checks and the original
orthogonal-coupling checks. Periodic coverage, mirror/basin counting, angular
partition, rate accuracy and lifetime remain unadmitted. Failure stops
unresolved without a method sweep or automatic redispatch.
