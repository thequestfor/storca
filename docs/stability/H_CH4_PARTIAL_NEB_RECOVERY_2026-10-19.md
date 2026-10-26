# H + methane: partial-NEB recovery implementation

Date: 2026-10-19. **Software implemented; native comparison blocked by RAM.**
The user authorized the controlled electronic comparison and narrow paid-stage
recovery. No new ORCA job or charged attempt resulted from this increment.

## Scientific result retained

Attempt nine converged its NEB band in 27 iterations, then aborted in ORCA's
DIIS/Fock diagonalization during internal OptTS cycle five. A native ORCA bug
is possible, but its exact cause remains unproved. The detailed diagnosis is
`runs/h-ch4-adaptive-authorized-20261018/diagnosis-20261019.md`.

The last evaluated geometry is cycle four, not `neb-ts.xyz`'s unevaluated cycle
five. Its remaining maximum gradient, 8.74622e-4 Eh/bohr, does not meet the
3e-5 requirement. Neither this geometry nor the saved climbing image is a
verified TS, rate or lifetime.

## Implemented changes

- `storca/partial_neb_recovery.py`: read-only admission of an owned converged
  NEB followed by failed refinement. It checks the failed native ledger receipt,
  mapped method/state/environment context and coordinates against terminal-owned
  output and MEP, not just the existence of a side XYZ. CI and evaluated OptTS
  candidates are distinguished from the failing diagnostic-only geometry.
- `reaction_path.run_generic_reaction_path` accepts the opt-in
  `retained_partial_neb_request` only with an explicit adaptive policy.
  `route_search_execution.py` refines its admitted candidate through the existing
  TS/frequency/mapped-mode/IRC/reached-minimum pipeline, with no new NEB fallback.
  A selected TRAH override is scoped to the TS verification chain and recorded
  in action identity/evidence; original fragment methods and defaults remain.
- `tests/run_h_ch4_scf_diagnostic.py`: exactly two fixed-geometry, fresh-guess
  single points, original solver versus explicit TRAH. Four cores, original
  method/basis/grid/state/accuracy targets and the original native ledger remain.
  Failed attempts stay charged, and six downstream SP slots remain reserved.
- `tests/run_h_ch4_partial_neb_recovery.py`: source-bound TS-only continuation.
  It refuses without both raw electronic solves, exact source/geometry lineage
  and frozen numerical-equivalence screens. Energy difference <=1e-8 Eh,
  S² difference <=1e-4, doublet deviation <=0.05 and integrated electron counts
  near 6 alpha/5 beta are recovery screens, not calibrated rate error bounds.
- Fixed a settings-scope collision during owned read-only revalidation. An
  identical active profile can be reused; a mismatched profile still refuses.
  Also fixed immutable adaptive publication to compare the same canonical bytes
  it writes, so an unchanged admission can actually reopen on resume.

All final TS, state, gradient, independent Hessian/mode, bidirectional IRC and
actual reached-species requirements remain unchanged. No new chemistry theory,
environment lane, production default, package, or experimental rate input was
added. No full NEB was repeated.

## Tests and preserved lineage

**409 combined regression tests pass**, including new pure admission,
diagnostic, wrapper and adaptive-integration tests. They cover arbitrary/wrong
side XYZ refusal, state/context/source mutation, the diagnostic-only failing
geometry, scoped solver selection, no NEB fallback, paid-work reuse, no partial
success, active-settings restoration and zero-charge resource blocks.
`git diff --check` passes.

Before editing the three previously bound implementation files, all 27 prior
continuation source files and the prior contract were preserved exactly in
`runs/h-ch4-adaptive-source-capsule-20261019`. Its source-capsule SHA256 is
`2de3f3d2cf489b54ae2402e9d9e63cffb32430a45452b09b992999cd7ab1f2dd`.
Old source guards were not waived. New diagnostic/recovery revisions bind the
archived old implementation and the exact new live inventory separately.

The first diagnostic preparation caught a pre-dispatch nested-settings error;
its fixed runner and r2 runner were archived before subsequent edits. The r2
and r3 launch attempts then refused at the unchanged memory guard. These were
all zero-native, zero-charge starts, not chemical failures. Earlier packages
are retained; **r3 is the current executable diagnostic revision**.

The original ledger remains:
`runs/h-ch4-authorized-20261017/campaign/native-control`.
Counts are Opt/path 5, frequency 3, IRC 0, SP 1: **nine charges, none pending**.
No failure was refunded. Remaining maxima are 7/7/2/8. Two successful diagnostic
attempts would leave six SP slots for the five matched energies and one reached
atomic point; no additional SP retry is implied.

## Current handoff and next action

Current authorized diagnostic:
`runs/h-ch4-scf-diagnostic-20261019-r3`.
Contract SHA256:
`a480226d3167a18983abf16ec3eb8dc349888fd6614a6b547669980e69b907e5`.
It reopens with zero native dispatch. `launch-status.json` records the resource
block. The guard requires 5,268,045,824 available bytes at actual dispatch.
A followup snapshot was only 5,290,999,808 bytes; available RAM fluctuates and
controller startup also consumes memory. Free roughly 1 GB of extra headroom;
do not lower the guard or close unrelated applications automatically.

From the repository using the existing application environment:

```sh
python -B tests/run_h_ch4_scf_diagnostic.py run runs/h-ch4-scf-diagnostic-20261019-r3
```

This resumes the same two arms on the same ledger. It will not repeat a charged
arm. Keep four cores, 1000 MB/rank, six hours/job and all cumulative safeguards.

Only after the comparison completes and the frozen equivalence checks pass:
freeze the explicit request selecting the last evaluated native geometry and
reviewed solver, prepare a fresh recovery workspace with the partial-recovery
driver, then run it with `--rmg-env rmg4_env`. Do not run the older adaptive
launcher against edited sources, fabricate the missing comparison report, or
treat passing software tests as evidence that TRAH fixes the native failure.
The recovery wrapper deliberately cannot prepare/run an admitted TS continuation
until those electronic results exist. The original crash remains unexplained
until controlled evidence resolves it; even two successful fresh solves would
not recreate the lost internal DIIS history or prove a particular ORCA defect.
