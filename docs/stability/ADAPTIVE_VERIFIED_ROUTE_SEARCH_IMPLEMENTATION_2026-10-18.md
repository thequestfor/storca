# Adaptive verified-route search: implementation handoff

Date: 2026-10-18. Status: **implemented, opt-in; native scientific performance
unvalidated**. No new chemistry was launched. This closes the software increment
in [the saved plan](ADAPTIVE_VERIFIED_ROUTE_SEARCH_PLAN_2026-10-18.md), not the
H + methane reaction investigation or the full stability roadmap.

Launch follow-up (2026-10-18): the user subsequently authorized execution.
The new `runs/h-ch4-adaptive-authorized-20261018` plan/contract preserves this
handoff's unauthorized r2 package and the original native ledger. Its kitty
launch reached the first native dispatch guard, which refused insufficient
available memory. No new ORCA job or attempt was charged; the eight historical
charges and all physical gates remain. See that workspace's `launch-status.json`
for the last resource snapshot and exact same-campaign resume command. Free
memory before resuming; no implementation or chemical validation failure was
established by this attempt.

Resume follow-up (2026-10-18): a subsequent user-authorized same-campaign
resume passed the unchanged resource guard. Native NEB is running in kitty
with an owned four-core lease and four observed MPI workers, under the original
six-hour job allowance. Original-ledger attempt nine is charged; the original
eight charges remain. The new initial path and endpoint energy/gradient
calculations completed. See `resume-status-20261018.json` in the authorized
workspace for the launch snapshot; later state must be checked live. Execution
progress does not establish a converged saddle, verified route or credible rate.

## What changed

- `storca/reaction_seed_geometry.py`: immutable mapped geometry context and
  paired seed portfolios; local outward H-transfer proposals, general rigid
  fragment translation/rotation, all-contact screening, proper-fit diversity,
  and content-bound reopening. Geometry qualification cannot verify a route.
- `storca/route_search.py`: deterministic persistent action frontier, scientific
  fingerprints, separate execution/reader/physical/applicability outcomes,
  completion-cost forecasts and explicit incomplete/budget/unsupported stops.
  A funded search is not limited to one try plus one retry.
- `storca/route_search_execution.py` and `route_search_recovery.py`: narrow
  adapters to existing NEB, guided TS recovery, frequency/mode, IRC continuation,
  intermediate and reached-minimum machinery. There is no second native runner
  or compute ledger. Stopped owned work is reused; live/ambiguous ownership
  prevents duplicate dispatch. Resource-only blocking does not poison retries.
- `reaction_path.run_generic_reaction_path(..., route_search_policy=...)` opts
  in. Legacy/default callers retain their existing sequence. Independent
  fragment reuse no longer requires a previously verified full path.
- `surface_resolution.EndpointReferenceAssessment` distinguishes a finite
  encounter from a qualified bound/fragment reference. A submerged climbing
  image against an unqualified encounter cannot establish a barrierless or
  dissociation mechanism. Bounded stopping remains available.
- `rate_transfer_continuation.py` binds archived historical sources, current
  sources, original raw evidence and the original native ledger. An intentional
  code revision is explicit lineage, not permission to ignore old hashes.
  Resumed unresolved adaptive reports retain separate immutable summaries and
  closeouts; they neither overwrite history nor remain bound to a mutable report.

Adaptive native inputs disable hidden fragment preparation and inspect the
actual initial path when present. Native movement is evidence to assess, not
proof of an asymptote. Final success still requires the intended mapped mode,
bidirectional mapped IRC and independently validated actual reached species.
A multistep result requires the complete ordered chain and all required segments;
preliminary TS/IRC evidence is never sufficient.

## Actual saved-evidence replay

Artifact: `runs/h-ch4-adaptive-geometry-replay-20261018-final/replay.json`.
The real seven-job fragment evidence was revalidated with all native dispatch
entrypoints forbidden. The charged eighth NEB attempt remains diagnostic.

| Geometry diagnostic | Original seed | Selected new seed |
| --- | ---: | ---: |
| Donor–H / H–attacker alignment cosine | −0.827032 | 1.000000 |
| Incoming H to closest spectator H | 1.030425 Å | 4.372811 Å |
| Incoming H to selected transfer H | 2.800000 Å | 2.800000 Å |

Both mapped endpoint connectivity/stereochemistry screens pass; 18 distinct
geometry proposals qualify. Original campaign file bytes, sizes, mtimes and
journal are unchanged. Native dispatch attempts and new jobs are both zero.
Portfolio SHA256:
`7946aed9a5acba03ec435ea78b07a0c6285c23c8de719d9f8dd855ebf63169f1`.

This establishes repair of the observed inward/crowded starting arrangement.
It does not establish a saddle, successful IRC, lower native cost, lowest
barrier, rate, uncertainty interval or lifetime.

## Tests and defects addressed

The combined campaign/new-feature regression run passed **787 tests**; a
separate broader compatibility run passed **322 tests**. These are overlapping
suite executions, not 1,109 distinct new tests. The final unresolved-report
resume fix passed an additional **63-test** targeted run. The focused reproducible
command below passed **141 tests**, including that regression and roadmap
traceability. Roadmap traceability also passed again after the document updates.

Covered defects include outward mapped placement, treating an encounter as an
asymptote, missing mapped-mode inputs in recovery, accepting partial verification
without reached-species evidence, wrongly skipping resumed paid NEB seeds,
hidden IRC-extension ownership, operational failures triggering geometry sweeps,
resource-blocked receipts preventing later admissible retries, and changed
unresolved summaries conflicting with immutable closeouts. Regression tests
also exercise multistep coherence, crash recovery, source mutation, attempt-nine
charging and unauthorized zero-dispatch refusal. No physical tolerance was
relaxed.

From the repository in the existing application Python environment:

```sh
PYTHONPATH=tests python -B -m unittest test_reaction_seed_geometry test_route_search test_endpoint_reference test_route_search_integration test_route_search_recovery test_rate_transfer_continuation test_h_ch4_route_search_replay test_rate_transfer_execution test_stability_roadmap_traceability
```

The verified interpreter here was
`/opt/homebrew/anaconda3/envs/orca-auto-py312/bin/python`. The two real-artifact
replay tests require this checkout's retained campaign; they skip without it.
Use `-B` to avoid generating bytecode in evidence-bound source directories.

## Preserved spending and current continuation

The source capsule at `runs/h-ch4-source-capsule-20261018/source-capsule.json`
retains 365 source/configuration files, including all 21 originally bound
implementation files at their exact historical hashes. Its SHA256 is
`1d77fca0f927a9eee2bc7834d31111ef142b978a9e25527056efbe4d7197c948`.

The original control root remains
`runs/h-ch4-authorized-20261017/campaign/native-control`. Its 16 immutable events
record eight charged jobs, with no pending reservation:

| Category | Spent, including failures | Original cap | Remaining maximum |
| --- | ---: | ---: | ---: |
| Optimization/path | 4 | 12 | 8 |
| Frequency | 3 | 10 | 7 |
| IRC | 0 | 2 | 2 |
| Single point | 1 | 9 | 8 |

Recorded owned wall time is 3,329.274526 seconds; allocated core time is
13,317.098104 seconds. All failures remain charged. Seven completed native jobs
support reusable H/CH4/H2/CH3 fragments; the stopped eighth NEB is not a verified
route. No historical work was redispatched. Four cores, 1000 MB/rank,
six hours/job, memory/disk safeguards and 48-hour wall/192-core-hour cumulative
ceilings remain unchanged. Counts are ceilings, not guaranteed available work:
the planner must still afford mandatory completion and downstream stages.

Current preparation: **`runs/h-ch4-adaptive-preparation-20261018-r2`**.

- Status: `adaptive_continuation_prepared_not_authorized`.
- Retained fragment package: four fragment rows, 45 raw source bindings.
- Contract SHA256:
  `d115a77e21c39568a10f93e1485d881a58081793bf4c7426d3ef13312801c5c1`.
- Plan SHA256:
  `70b55c85bed0f5bdc7270a2bd8d6150162d7beaecaa40e8277d96ef63b3408ca`.
- Five downstream matched single points remain in completion-cost forecasts;
  search cannot spend their capacity as though they were free.
- The earlier unsuffixed preparation is retained, superseded by this source
  revision. Do not overwrite its seal or bypass its expected source mismatch.

Safe next command, from the repository:

```sh
python -B tests/run_rate_transfer_benchmark.py inspect-continuation runs/h-ch4-adaptive-preparation-20261018-r2/campaign
```

Reproducing preparation in a fresh directory launches no chemistry:

```sh
python -B tests/run_rate_transfer_benchmark.py prepare-continuation runs/h-ch4-authorized-20261017/campaign runs/h-ch4-source-capsule-20261018 runs/h-ch4-adaptive-preparation-replay
```

**Next decision:** authorize a native H + CH4 adaptive continuation using this
implementation and the unchanged original remaining budget. On explicit consent,
create a new authorized plan/workspace, preserve this unauthorized package,
freeze its sources and use `run-continuation ... --rmg-env rmg4_env`. The first
new charged native attempt is number nine. The current sealed package must not
be edited to turn authorization on. The old frozen driver is not a way around
its source checks. OH + H2 and the closed OH/methane angular campaign are not
authorized by this handoff.

## Remaining scientific work

The H + methane route remains **unresolved**. Native continuation must determine
whether the improved approach gives the intended first-order saddle and reached
species. A finite frontier can exhaust without proving that no route exists.
Existing inexpensive preparation does not add general internal conformer
sampling, xTB screening, a new native S-IDPP reader, or specialized barrierless,
spin-crossing, surface or solvent-participation physics. Unsupported motions or
environments must remain explicit; no lane/support level is promoted here.

The separate three-arm efficiency protocol is prepared at
`benchmarks/kinetics/adaptive_route_search_20261018/README.md`: old placement/fixed
sequence, improved placement/fixed sequence, and improved/adaptive. H + methane
is development; OH + H2, H + ethene and cyclobutene are proposed transfer cases.
Before native comparison, finish the fixed-sequence B-arm wrapper and freeze
case maps/states, identical information, failed-inclusive scoring and separate
budgets. No speedup or general route-recovery probability is established yet.
Primary experimental table/condition/correction gaps from the prior transfer
handoff remain; they are not needed to declare geometric route success and
must not be filled with invented values or leaked into prediction inputs.
