# OH + methane: resume native work without a loader-log gate

Date: 2026-10-16. Status: **implemented; 330 focused tests pass; actual +120°
reused; +170° running in kitty on four cores at optimization cycle 2.**

Follow-up (2026-10-17): the controller completed all remaining scan points and
the selected +170° NumFreq. The sealed outcome is
`profile_complete_partition_unsupported`, cumulative spending 13/4. The
[closeout](OH_METHANE_CLOSEOUT_2026-10-17.md) supersedes the launch state below.

## The correction

The +120° four-core optimization in
`runs/oh-methane-periodic-candidate-20261016-ha4k1u9s` terminated normally after
43 cycles (19 minutes 30 seconds). STORCA subsequently stopped at a mandatory
debugging check of interleaved macOS loader output, not a chemistry or MPI
execution failure. The raw output contains cross-written path fragments.
Determining every worker's loaded library from that stream is not a necessary
condition for accepting an otherwise owned, normally completed calculation.

The new `storca/configured_mpi_runtime.py` validates explicit runtime selection
and records detailed loader evidence as optional diagnostics. Missing or
garbled loader evidence is `unresolved_diagnostic`, not a native failure.
Contradictory runtime configuration still refuses. Unrelated programming errors
are not swallowed. Configuration alone is explicitly **not** proof of the
library loaded by every rank, much less proof of chemical validity.

This continuation retains the bound private Open MPI 4.1.1 installation, ORCA
binary/runtime inventories, exact launch environment, four-core owned lease,
input/output hashes and unique normal termination. The existing environment
still emits tracing; this narrow correction changes its status from mandatory
acceptance evidence to diagnostic evidence, not the environment itself.

## Paid work is reused

The unchanged physical readers admit the actually evaluated +120° geometry
**43**, not final energy-only geometry **44**, as an optimizer candidate:

| Diagnostic | Value | Existing limit |
| --- | ---: | ---: |
| Tangent gradient RMS | 3.09928e−6 Eh/bohr | 8e−6 |
| Tangent gradient maximum | 8.61664e−6 Eh/bohr | 3e−5 |
| Angular error | 7.01019e−8 rad | 2e−4 |

Native method/state, spin continuity and mapped spectator topology also pass.
The zero-job replay preserves cumulative consumption at **8 optimizations /
3 NumFreq**, with original **17/4 caps**. No +60° or +120° job is repeated.
The next new target is +170°, followed by the unchanged negative arm, labelled
seam closure and one selected full Hessian. Six hours/native job, 1000 MB/rank
and the 4 GiB disk reserve remain unchanged.

Launch verified: `point-p170/optimization/opt.process.json` is an owned running
ORCA lease with `ncores: 4`; native output has reached cycle 2. Its reservation
brings cumulative spending to **9 optimizations / 3 NumFreq** (eight unspent
optimization slots and one Hessian remain). The +120° analysis and checkpoint
bind the runtime extension; that completed optimization was not relaunched.
The first desktop-launch approval review timed out before launching anything;
the permitted single retry succeeded. No duplicate native job was started.

## Small, explicit reader extension

`tests/resume_oh_methane_angular_runtime.py` wraps the frozen candidate driver
without editing any source already bound by that campaign. A fixed, exclusive
`runtime-policy-extension-v1.json` and its binding record identify the new
policy, exact reader sources, original contract, attention checkpoint and paid
journal prefix. Future journal entries may append; the bound prefix cannot
change. There is no ledger fork, counter reset or refunded attempt.

Four scoped callbacks cover runtime diagnostics, pre-dispatch validation,
checkpoints and final closeout. Only current journal-owned normal outputs get
the new runtime reader. Historical +60° and pilot receipts retain their exact
original reader and replay. Every new runtime receipt and final result binds
the extension, even when loader tracing happens to parse cleanly. The callbacks
are restored on exit, including exceptions; production defaults are unchanged.

`runtime-policy-plus120-replay.json` retains the actual zero-dispatch replay.
Completed +120° analysis then supplies the upstream seed for +170° through the
unchanged original schedule. Starting at +170° by discarding the +120° node
would lose that seed and is intentionally not done.

## Validation and remaining scientific limits

All **330 focused tests pass**, including 15 new pure-runtime tests and 15 new
adapter tests. They cover historical replay, clean and damaged diagnostic logs,
configuration/ownership/native-completion refusal, callback restoration,
changed source/journal bindings, unchanged counters, zero-job paid replay and
preservation of physical failure gates. The actual +120° physical replay also
passes with native dispatch forbidden.

An initial combined test invocation from standard input could not spawn the
ledger's process-death test. Rerunning through the standard unittest module
entry point resolves that harness limitation; the full 330-test run passes.

An optimizer candidate is not a Hessian-verified saddle. Selected independent
Hessian, H-transfer mode, branch closure, coupling and all rate-admission gates
remain unchanged. No reviewed global angular Hamiltonian/measure has been
admitted, so this increment still produces **no corrected rate or lifetime**.

Resume in the existing configured environment:

```sh
python -B tests/resume_oh_methane_angular_runtime.py runs/oh-methane-periodic-candidate-20261016-ha4k1u9s
```

The visible kitty launcher is
`runs/oh-methane-runtime-resume-20261016.command`. The previous candidate driver
remains available for historical replay; use this explicitly versioned adapter
for the revised runtime policy.
