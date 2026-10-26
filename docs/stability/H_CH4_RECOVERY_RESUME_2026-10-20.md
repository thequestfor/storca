# H + methane recovery resume

Date: 2026-10-20 (client date). Native work resumed; no implementation source
was changed in this increment. This document records a launch checkpoint,
not the final result; verify live ownership and artifacts before another run.

## Completed electronic comparison

`runs/h-ch4-scf-diagnostic-20261019-r3/comparison.json` reopens against both
owned raw outputs. Attempts 10 and 11 completed normally and pass the frozen
numerical energy, spin, electron-count, charge and multiplicity screens.

| Quantity | Original solver | Explicit TRAH |
| --- | --- | --- |
| Energy / Eh | -41.000557502889 | -41.000557502890 |
| Printed S² | 0.759080 | 0.759080 |
| Integrated alpha electrons | 5.999999809553 | 5.999999809553 |
| Integrated beta electrons | 4.999999924969 | 4.999999924969 |
| Owned elapsed time / s | 10.942878 | 7.622303 |

TRAH minus original energy is -9.947598300641403e-13 Eh. Both fresh-guess
solves succeed at the previous failing geometry; this does not reproduce
the lost optimizer DIIS history, establish the original crash cause, or prove
physical rate accuracy. The numerical screens are not an uncertainty estimate.

## Frozen recovery and running work

The explicit request selects `last_evaluated_opt` (cycle four), not the
unevaluated failing cycle-five geometry:
`runs/h-ch4-partial-neb-recovery-20261020/partial-request.json`.
The original context, candidate sources and comparison are hash-bound. An
initial preparation refused a JSON numeric-representation mismatch in the
new request; preserving the original floating-point representations resolved
it before sealing, without source changes, native dispatch or another charge.

The new authorized contract is:
`runs/h-ch4-partial-neb-recovery-20261020/campaign/recovery-contract.json`
with SHA256 `f0ed946002e04082504c20af03c2c58b3747e7b7a91b46a1b952f6840e1a9c2e`.
Preparation and inspection both passed. The original ledger, old source
capsule and exact live implementation inventory remain bound.

The TS-only continuation started attempt 12. The retained native input is
under `campaign/path/route-search/actions/cb5b3175d5513671af34584544fe245a450c07a9faf4533ee5fb46799e1cac80/ts-opt.inp`.
It uses the original WB97X-D3/def2-TZVPP model plus explicit TRAH, four cores,
1000 MB/rank, VeryTightOpt, an initial optimizer Hessian and recalculation every
five cycles. The unchanged independent final Hessian/mode, IRC and reached-
minimum gates remain required. No NEB rerun, threshold relaxation, new ledger,
production-default change or package installation occurred.

At launch the original ledger contains 12 charges: optimization/path 6,
frequency 3, IRC 0, single point 3. All failed historical work remains charged.
The two diagnostic single points leave the reserved six downstream SP slots.
Four actual native workers were observed after startup. Per-job allowance is
still 21,600 seconds, with the original cumulative and resource safeguards.

## Subsequent verified checkpoint

Attempt 12 completed its TS optimization in four cycles. The final maximum
gradient is 7.363e-7 Eh/bohr against the unchanged 3e-5 threshold; all printed
geometry convergence criteria pass. Attempt 13, the independent NumFreq,
terminated normally with exactly one imaginary frequency at -1317.42 cm-1.
The existing mapped reaction-motion gate passed before IRC dispatch.

Attempt 14 is the owned bidirectional IRC, now running with four observed
workers and the original six-hour allowance. The original ledger therefore
has 14 charges: optimization/path 6, frequency 4, IRC 1, single point 3.
There is one active native reservation, not a missing/duplicate calculation.
ORCA IRC parent observed at this checkpoint: PID 94106; controller remains
PID 92843. Final mapped endpoints and reached minima have not yet been verified.
The original crash remains unexplained; successful recovery is not proof of
its cause. No rate or lifetime is available.

Roadmap traceability tests: four passed. `git diff --check` passed. No scientific
implementation edits were made, so the previously recorded 409-test result is
historical implementation evidence, not a new test run in this increment.

## Continuation instructions

Do not start another controller while the current one owns native work.
Controller at launch: PID 92843; ORCA parent: PID 92879. PIDs are historical
observations, not sufficient identity evidence for terminating or resuming work.

After checking that the original owned process has ended, the existing command
is resumable from the repository/application environment:

```sh
python -B tests/run_h_ch4_partial_neb_recovery.py resume runs/h-ch4-partial-neb-recovery-20261020/campaign --rmg-env rmg4_env
```

Do not rerun the old adaptive launcher, edit sealed sources, reset charges,
or repeat the two diagnostic arms. A completed optimizer alone is not a
verified saddle or reaction route. No rate, lifetime or network reinjection is
admitted at this checkpoint.
