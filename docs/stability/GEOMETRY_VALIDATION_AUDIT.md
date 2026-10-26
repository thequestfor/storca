# Why the current stability run stopped

Audit checkpoint: 2026-10-01. This closes a bounded geometry-validation repair,
not the working-prediction milestone.

## The immediate cause

RMG did discover cyclobutene ring opening. The production run then stopped
before launching an ORCA reaction-path calculation: STORCA's endpoint checker
rejected every generated seed. Opposite carbons across the four-membered ring
were about 2.05 Å apart. A generic covalent-radius rule treated those legitimate
nonbonded diagonals as ambiguous extra bonds. The retained ORCA-optimized parent
had the same problem with that rule; this was not evidence that ORCA generated
an invalid ring.

The run therefore returned `verification_incomplete` with no lifetime. It
failed closed, but for the wrong reason. The correction does not relax the
global cutoffs: a closed-shell four-ring diagonal is accepted as nonbonded only
when bounded, independent coordinate-based graph perception exactly matches
the declared atom-ordered bonds, bond orders, formal charges and radical counts.
Failed assignments and other close contacts remain unresolved.

## Your concern exposed a broader problem

The old checks conflated different questions: “are these coordinates close?”,
“is this the declared molecule?”, and “is this a physically valid stationary
point or reaction?” They are not interchangeable. Constructed negative tests
reproduced additional weaknesses:

- A mirror stereoisomer or opposite E/Z geometry could pass connectivity checks.
- Bound IRC endpoints were checked by coordinate similarity without the mapped
  connectivity check used for separated fragments. A broken H–H bond could
  therefore pass that endpoint check.
- Summed endpoint deviations allowed one perfect end to compensate for another
  exceeding the stated per-endpoint tolerance.
- Nonfinite monatomic coordinates and incomplete/corrupted trajectory frames
  were not consistently rejected; an empty trace raised an unclassified error.

These demonstrated false acceptances are repaired. Declared tetrahedral/E/Z
geometry is checked, both IRC ends need mapped chemistry checks, bound ends
also need separate coordinate tolerances in the same direction, and every
frame needs intact element order and finite coordinates. Empty or invalid atom
maps are refused. Geometry evidence explicitly disclaims stationary-point and
electronic-state validity. ORCA convergence, frequencies and IRC completion
remain separate requirements.

A separate reopening defect confused Python tuple pairs with JSON list pairs
in atom-map evidence. Exact integer-pair normalization fixes that without
accepting a changed map. The new verifier policy is versioned; older contracts
remain immutable and require explicit policy amendments rather than acquiring
stronger evidence retroactively.

## What is and is not established

The earlier tests did not challenge these specific identity and malformed-data
cases. Adding software checks did not establish their chemical adequacy; this
audit adds the missing adversarial cases and replays actual retained geometry.

Evidence: the [before audit](../../runs/p1-geometry-audit-before-001/result.json)
retains nine failure entries and two errors across 15 methods. The
[expanded after audit](../../runs/p1-geometry-audit-after-002/result.json) passes
20 methods, including a 16-geometry constructed reference panel. The
[actual boundary replay](../../runs/p1-cyclobutene-ring-003/result.json) accepts
the retained seed and parent while leaving all 111 original files unchanged.
A suspected negative-count parser hang did **not** reproduce: the existing
parser already rejects it; only its regression test was added.
Final full application suite: 794 tests, OK with 22 expected native-only skips.
Focused RMG/Arkane-environment suite: 68 tests, OK. `git diff --check` passed.

This is not a general electronic bond-order validator, a conformer-coverage
study, an experimental benchmark, or proof of spin/multireference applicability.
Unspecified stereochemistry is not resolved by this check. No new quantum path,
rate, repaired propagation, or validated lifetime was produced during this audit.

These defects were tractable; switching goals is not necessary to repair them.
The next scientific test remains the bounded real path → applicable rate →
repaired kinetics chain, followed by interruption/resume. Rate uncertainty and
independent mechanism coverage must still be resolved before claiming a
defensible prediction. See [the acceptance ledger](WORKING_PIPELINE_GATE.md).
