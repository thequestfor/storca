# Bounded known-reaction execution adapter

Date: 2026-10-17. Status: **implemented, synthetic/native-reader regression
tested; not yet validated by a transfer chemistry campaign**. Zero new ORCA
calculations were launched. The OH + methane closeout, original frozen
candidate and original transfer preparation remain unchanged. This increment
does not admit a physical rate, lifetime, network replacement or new environment.

## What changed

The preparation-only transfer protocol now has a separate opt-in executor:

```text
answer-free known reaction (H + CH4 first; OH + H2 second)
  → separately frozen native settings, assumptions and compute authorization
  → existing generic mapped endpoint / NEB / TS / mode / bidirectional IRC checks
  → independent refinement of the actual reached species
  → five geometry-matched high-level electronic energies
  → source-separated composite recipe and isolated native Arkane replay
  → conditional elementary-rate diagnostic, never automatic network reinjection
```

The known-reaction protocol deliberately supplies both products. It does not
impersonate RMG discovery, flux evidence or autonomous chemical-fate prediction.
Experimental answers and the scoring dossier do not enter its prediction calls.

- `storca/rate_transfer_execution.py`: draft/freeze/reopen/run/resume adapter;
  ordinary route verification remains in `reaction_path.py`. Raw optimization,
  frequency, Hessian, IRC and reached-species sources form a stable scientific
  handoff; mutable resume metadata is not a calculation identity.
- `storca/route_geometry.py`: complete explicit nuclear maps survive RDKit
  parsing, including the specifically transferred H and standalone H/H2.
  Partial, duplicate and isotope-conflicting maps refuse. Legacy heavy-only
  routes retain bounded hydrogen mapping; explicit maps cannot be rescued by
  silently permuting the declared hydrogen identities.
- `reaction_path.py` / `route_thermochemistry.py`: neutral atomic H has typed,
  source-bound native energy, mass, state and ideal-gas translation/spin
  evidence. It has no invented optimization, vibration, Hessian, rotor or ZPE.
  Molecular minima still require their ordinary optimization/frequency evidence.
- `src/orca_job_settings.py`: explicit opt-in VeryTightOpt, SCF and central
  NumFreq settings reach all participating builders and reuse identities.
  Outside this scope existing production inputs/defaults remain unchanged.
- `storca/rate_transfer_control.py`: exclusive controller, immutable attempt
  reservations/terminal receipts, completed-work reuse, conservative failure
  accounting, memory/disk/core and wall/core-time safeguards. Actual owned job
  duration is charged, not idle time between invocations. Stopped controller
  crashes may be reconciled from native process evidence without refund or
  rerun; ambiguous/live ownership blocks. Unknown job categories refuse.
- `composite_rate_recipe.py` / `arkane_rotor.py` / `arkane_launcher.py`: explicit
  atomic-SP variant without a fictitious Hessian; historical molecular-only v1
  artifacts retain their schema and validation. Atomic H uses a declared native
  one-electron HF calculation in the same QZ basis, not an assumed −0.5 Eh or a
  silent failed-CC fallback. Composite D/T remains unsupported, even though
  mapping and low-level atomic thermochemistry can retain isotopes.
- `storca/native_rate_audit.py`: composition-derived native electron/basis/core
  inventories and reference-free TST, detailed-balance, units and E0 checks.
  Both Eckart barriers must be positive. The selected `--rmg-env` is used;
  installed RMG/Arkane is not patched. Failed/ambiguous native replay is retained,
  not silently retried; successful replay is source-checked before reuse.

The numerical procedure is transferred, not OH/methane's numerical symmetry,
optical/channel multipliers, energies or fitted discrepancy. The atomic and
CC native readers remain narrowly scoped to the declared ORCA 6.1 H/C/O basis
inventory; unsupported outputs require evidence, not guessed compatibility.

## Budget correction: proposal only

The original preparation proposed 6 optimization/path and 4 frequency jobs.
The existing generic path requires **at least 8 optimization/path, 7 molecular
frequencies, 1 bidirectional IRC and 7 single points** for either transfer:

- Three molecular minima before the path and three independently refined
  actual IRC-reached molecular minima: 6 optimizations + 6 frequencies.
- NEB and separate TS refinement: 2 optimization/path jobs + 1 TS frequency.
- Two low-level atomic-H energy jobs and five matched high-level energies:
  7 single points. No fictitious atomic frequency is counted.

The new draft proposes **12 optimization/path, 10 frequencies, 2 IRC and 9
single-point attempts per case**, including bounded recovery reserves. Every
dispatch, including failure, consumes an attempt. Internal images/displacements
consume the owning job's time, not uncounted independent jobs. The original 6/4
proposal is not rewritten or treated as sufficient.

Retain four cores, 1000 MB/rank, a 4 GiB disk reserve, 1 GiB additional available
memory reserve and six hours per job. The proposed aggregate electronic-job
ceilings are 48 owned wall-hours / 192 allocated core-hours per case; these are
caps, not estimates or authorization. The single selected-environment native
statmech/rate replay has a separate six-hour ceiling and no automatic retry.
Its duration is recorded separately, not misreported as ORCA time. No production
default changes, old ledger reset or failed-attempt refund follows.

## Regression evidence and limits

The combined focused suite exercises historical OH closeout/freeze, transfer
preparation, new settings/control/mapping/atomic/bridge/audit code and existing
generic-path, reached-IRC, thermochemistry, P1 and physical-gate tests. Synthetic
tests cover an interrupted first-SP handoff and repeated resume with changing
runtime metadata; retained work and final result identities are preserved.
Actual composite staging/reopening is exercised with synthetic mapped OH + H2
native records, without mocking composite validators. Native reader checks
also reopened retained OH/methane correlated outputs; none was recalculated.

The 631-test combined campaign suite and overlapping 458-test broader
compatibility suite pass. The latter includes saved native Arkane loading and
existing mapping/geometry/foundation/verification tests. Whitespace checks pass.
Focused transfer-only reproduction (synthetic chemistry/owned fake engines):

```sh
PYTHONPATH=tests python -B -m unittest test_rate_transfer_execution test_rate_transfer_execution_bridge test_rate_transfer_control test_native_rate_audit test_frequency_transfer test_transfer_mapping_atoms
```

These are software and ingestion tests, **not a real H + CH4 or OH + H2 accuracy
result**, and not proof that a newly searched saddle will pass. The closed
OH/methane candidate reopens `diagnostic_not_admitted`, still 13 optimization
attempts / 4 NumFreq and zero new native jobs.

## Ready artifacts and next decision

Execution drafts (not authorized):

- `runs/h-ch4-execution-draft-20261017/execution-plan.json`
- `runs/oh-h2-execution-draft-20261017/execution-plan.json`

They leave the ORCA runtime binding and case-specific symmetry/optical counts
and physical assumptions unset. No OH/methane counts are copied to fill them.
Electronic-state/SOC and H2 rotational-population assumptions remain explicit
scientific choices. A frozen binary and environment strings do not alone prove
every inherited MPI/library component: bind the selected runtime for the future
campaign, without reinstating noisy loader logs as a physical gate.

Primary experimental tables, exact conditions and the H + CH4 published
correction still require audit. The current temperature grids are diagnostic
selection targets, not manufactured measurements. The executor therefore
labels its model isolated canonical gas and comparison-ready false; it does
not invent 1 bar, a bath mixture, experimental uncertainty or state preparation.

The next task is to resolve/freeze these case-specific declarations, retain an
explicit runtime inventory and obtain separate authorization for the revised
**H + CH4 only** compute campaign. Then execute one frozen campaign and report
all physical failures/costs. OH + H2 follows; do not launch both by default.
Reference audit can proceed independently; it is required before experimental
scoring, not a reason to feed reference answers into prediction.

Zero-calculation draft reproduction, using the ordinary activated STORCA Python:

```sh
python -B tests/run_rate_transfer_benchmark.py prepare runs/rate-transfer-preparation-20261017 h_ch4 runs/h-ch4-execution-draft-20261017
```

After a separately completed and authorized plan exists, freeze it into a fresh
campaign directory with the `freeze PLAN WORKDIR` subcommand. Only then does
`run WORKDIR --rmg-env ENV` or `resume WORKDIR --rmg-env ENV` execute. Those are
future authorized steps, not commands run in this increment. Do not edit the
original closed OH/methane evidence or old preparation to bypass a refusal.
