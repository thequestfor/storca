# Rate-to-network handoff checkpoint — 2026-10-08

The focused handoff implementation is complete. Scientific P1 acceptance and a
physical stability/lifetime conclusion are **not** complete.

## What changed

- Verified path/minimum evidence and conditional forward rates remain visible
  when graph-species network admission is blocked. A conformer-basis refusal is
  not reported as failed ORCA/Arkane execution.
- New Arkane records retain exact-condition forward and reached-conformer
  reverse fit evidence. Old artifacts can be reopened read-only.
- A conservative finite-horizon return-flux assessment distinguishes a bound
  between named fixed models from a physical rate/confidence bound.
- A separately named irreversible approximation library may be admitted only
  with a recomputable, exact-condition, source-hash-bound certificate. The
  original reversible library and unresolved conformer basis are unchanged.
- Repaired-model inventory is checked again. Conditional rates earn no verified
  observable/flux credit, terminal kinetic eligibility, t95 or lifetime.
- Main screen, resume and ladder reports distinguish conditional predictions
  from calculation failures, retain evidence gaps, and continue the ladder.
- Ordinary gas Arrhenius Chemkin exports now retain 17-digit precision. Native
  RMG's previous three-decimal n/Ea export shifted this rate by about 0.022%,
  incorrectly breaking the existing rate-match guard. No guard was weakened,
installed package changed, or historical result rewritten.

Final verification: application suite **1,257 tests / 68.205 s**, 23 expected
skips; selected native RMG suite **70 tests / 3.162 s**, no skips; focused
screen/resume/ladder/handoff integration **102 tests / 1.285 s**. `git diff
--check` passes. Software correctness is not scientific calibration or P1
acceptance.

## Real saved-artifact replay

Conditions remain cyclobutene, 423.15 K, 0.0666611842 bar, pure dark gas,
isothermal batch, 3600 s, 95% retention threshold. These are development model
conditions, not a shelf-life or safety determination.

The saved Arkane forward fit evaluates to 3.548599062138614e-5 s^-1. The two
declared reverse alternatives are 6.202760921694612e-11 s^-1 for the reached
endpoint fit and 1.1152867251873439e-11 s^-1 using native graph-species Kc.
Neither is a bound on unknown physical reverse rates or unsampled conformers.
The omission bound relative to these fixed models is 2.2329939318100603e-7
absolute parent fraction at one hour.

[Full repaired RMG replay](../../runs/p1-cyclobutene-handoff-002/verification-result.json)
is **unresolved**, not a failed path/rate calculation. The rate is exactly
reinjected (relative error zero), but RMG expands the original six-species,
one-reaction core to seven species and three reactions, including a triplet
biradical intermediate. The original narrow certificate does not cover that
expanded mechanism. No physical conclusion is promoted.

[Frozen retained-model probe](../../runs/p1-cyclobutene-frozen-handoff-003/result.json)
passes native graph matching, actual library insertion, precision-preserving
export/reload, structural identity, repaired-model certificate and independent
Cantera propagation. It intentionally retains the original one-reaction model,
excluding the later discovered branches. It is **not** a substitute for the
expanded-network assessment or evidence of completeness.

Its model-only curve has about 88.0% parent remaining at one hour. Native
presented retention is 0.880093014676585; the independently calculated one-step
analytical value is 0.8800737517304962. Their numerical difference is
1.9262946088782762e-5, larger than the reverse-omission bound. The two errors
are recorded separately: the omission bound is not total prediction uncertainty.
Physical uncertainty and high-pressure applicability are unresolved.

All probes preserve original source bytes and modification timestamps. No new
electronic-structure or Arkane jobs were launched. The full replay did run new
RMG mechanism generation; the frozen probe did not.

## Reproduction

After activating the existing application Python environment, run from the
repository root, choosing an output directory that does not already exist:

```sh
python tests/replay_forward_rate_handoff.py --output runs/handoff-replay-new --rmg-env rmg4_env --propagate
python tests/replay_retained_forward_model.py --output runs/frozen-handoff-replay-new --rmg-env rmg4_env
```

The first exercises new RMG growth and may correctly refuse the original
certificate. The second explicitly tests the frozen known-incomplete model.
Neither invokes ORCA or repeats Arkane jobs. These are development replay
commands, not new user-facing stability options.

## Next scientific work

1. Assess the newly discovered initiating/return branches under stable native
   species and electronic-state identities. Do not collapse same-formula singlet
   and triplet species or silently reuse the one-reaction certificate.
2. Determine which additional branches influence the requested retention;
   verify/resolve influential uncertain branches and their method applicability.
3. Establish applicability of the rate model at the requested low pressure and
   supply defensible rate/conformer/counting uncertainty and coverage evidence.

These are distinct from repairing an execution error. Further endpoints, NEB,
IRC or GOAT calculations are not automatically required to redo the completed
ring-opening path.

The subsequent bounded implementation is specified in
[Expanded-network assessment and next-evidence plan](EXPANDED_NETWORK_IMPLEMENTATION_PLAN_2026-10-08.md).
It is now implemented: full native identities, a separate current-model
assessment, declared branch scenarios, pre-path electronic-state audit, and
saved-run/reporting integration. See
[its checkpoint](EXPANDED_NETWORK_ASSESSMENT_2026-10-08.md) for source-bound replay
evidence. It does not allocate chemistry jobs or promote any scientific support
gate; the original unresolved expanded replay remains historical and unchanged.
