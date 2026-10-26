# OH + methane: bounded angular investigation closeout

Date: 2026-10-17. Scientific outcome:
**`profile_complete_partition_unsupported`**.

The native angular campaign is finished. There is no native execution failure
to repair and no additional OH/methane calculation is authorized or needed for
this closeout. The remaining rate-model questions are physical assumptions,
not another reason to tune methods until a measurement is reproduced.

## What completed and what it means chemically

Reaction: OH + CH4 → H2O + CH3 on the declared neutral-doublet gas surface.
The earlier bidirectional IRC and matched plan five-species composite energy/rate
audit remain separate, retained evidence. This investigation tested the soft
relative angular motion near that reaction saddle; it did not discover a new
reaction or establish a methane lifetime.

All nine prescribed coarse nodes are present: 0°, ±15°, ±60°, ±120° and ±170°.
The independently grown +190° control closes against −170° without reflection
or spectator permutation. All optimization candidates pass the unchanged
mapped-angle, evaluated-frame gradient, native-state and spectator guards.
Four nodes—0°, ±15° and +170°—have independent full numerical Hessians; the
other five remain optimizer candidates, not independently verified saddles.

| Diagnostic | Observed result | Predeclared requirement |
| --- | ---: | ---: |
| Sampled DFT electronic-energy span | 0.07576 kJ/mol | Descriptive, not an accuracy bound |
| Closure mapped mass-weighted RMSD | 0.00009238 Å | ≤0.02 Å |
| Closure maximum atomic displacement | 0.0002934 Å | ≤0.05 Å |
| Closure electronic-energy mismatch | 0.001645 kJ/mol | ≤0.10 kJ/mol |
| +170° constrained tangent curvatures | 1 negative, 13 positive | 1 negative, 13 positive |
| +170° H-transfer mode overlap | 0.97388 | ≥0.90 |
| +170° tangent gradient RMS / maximum | 2.9902e−6 / 7.4653e−6 Eh/bohr | ≤8e−6 / ≤3e−5 |
| +170° orthogonal ZPE shift | −194.93 J/mol | Absolute value ≤500 J/mol |
| Maximum absolute orthogonal log-partition shift | 0.014705 | ≤log(1.20) |

The branch is very soft at the sampled geometries, and the wide sampled point
retains the intended carbon–hydrogen breaking / oxygen–hydrogen forming mode.
Small sampled electronic-energy changes and local transverse-mode changes are
useful evidence; they do not prove continuous branch completeness or the
correct kinetic contribution of the soft motion. The largest angular sampling
gap remains 60°, with a 20° seam gap.

The +170° **full** Hessian has imaginary frequencies near −1233.05 and
−6.33 cm⁻¹. The fixed-angle tangent Lagrangian Hessian passes its one-unstable-
direction requirement. This is a verified **constrained** point, not an
ordinary first-order transition state. No negative mode is deleted or hidden.

## What remains unresolved—and is not being fitted away

No reviewed reduced angular Hamiltonian or kinetic measure is implemented.
The local inverse-constraint metric is not demonstrated relaxed-path inertia.
Overall rotation, the reaction normal, transverse modes and angular motion
need a consistent treatment, with harmonic/free-rotor and reparameterization
checks. External symmetry, methane abstraction equivalence, spectator
permutations, mirror configurations and optical counts must be reconciled
together. A full angular integral cannot then be multiplied by the same basin
count again.

Sparse Hessians cannot establish transverse stability everywhere. The earlier
OH two-level treatment also retains assumed TS spin–orbit quenching and
electronic-state separation. Classical TST/Wigner/Eckart alternatives and
one-/two-harmonic-basin conditional branches are retained, not selected by
agreement with the exposed development experiment. An isolated gas model is
not identical to the humid-air, finite-pressure experimental environment.

Consequently there is **no new angular partition, corrected rate, calibrated
uncertainty, production-default change, repaired-network rate or lifetime**.
The existing rate-admission consumer independently records the refusal and
reopens 602 bound sources. This is a completed bounded investigation, not
general validation of chemistry, reaction discovery or any environment lane.

## Software and retained evidence

The source campaign is
`runs/oh-methane-periodic-candidate-20261016-ha4k1u9s`:

- `result.json` and `result-binding.json`: sealed outcome, coverage, closure,
  point evidence, coupling diagnostics and attempt history.
- `angular-rate-assessment.json`: explicit quantitative-model refusal.
- `runtime-policy-extension-v1.json`: explicit diagnostic-only loader policy.
- `attempts/`: the original immutable identity and 21-event journal.
- `point-p170/numerical-Hessian/`: actual final native input, output, Hessian,
  gradient, owned four-core process record and execution receipt.

The final native NumFreq completed in 638.676 seconds; the execution wrapper
recorded 638.823 seconds. All seven optimizations and the one selected Hessian
in the current periodic increment completed normally. Cumulative accounting is
**13 optimization attempts / 4 NumFreq**, under the original **17/4 caps**.
The cumulative totals include historical failures; they are not 17 successful
calculations. Two earlier fixed-frame EnGrad diagnostics are separate costs.
The unused four optimization slots are not spent on a partition model that
has not been admitted. No new Hessian may be started under this cap.

The 17 angular attempts comprise **14 normal native completions and 3 native
execution failures**. All 17 wrapper timings are retained, summing to
13,466.013624 seconds (3.74056 hours) and 14.96224 allocated core-hours at four
requested cores. These are summed job-wall/resource-allocation costs, not
measured CPU utilization or the calendar duration of this investigation. The
two separate EnGrad diagnostics have native run times of 914.167 seconds
(four cores) and 47.345 seconds (one core). Their wrapper timings are absent
and remain null; they are not silently added to the wrapper subtotal. Earlier
geometry/derivative/composite-energy development is outside this angular-cost
scope.

This task observed the existing controller and launched **zero new chemistry**.
A closed-controller replay with native execution forbidden reused the result,
preserved the journal and confirmed the same 13/4 totals. A separate raw
Hessian/gradient/mode replay, orthogonal-coupling replay and labelled closure
replay reproduce the stored receipts. No stored success flag substitutes for
those raw-file calculations. The initial one-off closure replay used an
incorrect Python module attribute; correcting that verification harness gave
identical results without changing sources or launching chemistry.

The final combined **424 focused tests pass** (5.903 seconds): the existing
330-test group, 20 new candidate-freeze tests, 20 new transfer-preparation tests
and 54 existing scientific-panel isolation/reference tests. These are software
and retained-evidence checks, not experimental validation. `git diff --check`
passes. No claim is made that the entire application suite was rerun.

No new native-run implementation defect appeared in this closeout. Previously
corrected STORCA defects remain explicit: an extra displacement acquisition
veto, a mixed-frame OptTS EnGrad interpretation, and a mandatory audit of
interleaved MPI debugging text. Their preserved failures are not relabelled
as ORCA chemistry failures. Frozen source files, methods, physical thresholds
and ordinary-TS validation remain unchanged.

One genuine defect appeared in the **new closeout utility**, not the chemistry
workflow: execution timings referenced inside the explicit journal events were
missing from the first freeze's cost subtotal. The original subtotal was
labelled partial, with unknown timings left null. A small correction opens only
those already-bound event payloads and retains the seven missing execution
bindings. Regression tests cover that case and changed-event rejection. The
first bundle remains intact; the final complete-cost bundle below supersedes
it. Native return code zero alone is not called success: the cost ledger also
requires the normal-termination marker and terminal lease evidence.

## Frozen candidate and reproducibility

Final bundle:
`runs/oh-methane-candidate-freeze-20261017-complete/candidate.json`.
SHA256:
`2cdc2464e0be231775d5f30ede126bd4ac74b18ee64ba0304d6b03c48f09d4fa`.

It binds the sealed profile/refusal, original five-species composite recipe,
electronic-state replay, exact source/runtime evidence and all-outcome costs.
There are 97 copied Python files, covering 98 source-version bindings including
the existing exact historical-generator reconstruction, plus 416 external
artifact bindings. Native data and runtime binaries are hash-bound, not copied
as redundant gigabytes. Copied Python versions can be verified after future
checkout edits; native source mutation still refuses. Byte-integrity reopening
is not a fresh physical certificate or a replacement for exact-version semantic
replay. Tests distinguish those guarantees.

The frozen physical procedure remains canonical CCSD(T)/def2-QZVPP electronic
points on WB97X-D3/def2-TZVPP geometry/vibration sources, with separate energy
and derivative bindings. Native energy inputs explicitly retain `FrozenCore`
and `UseQROs false`; the TS output reports a UHF reference and chemical core
(four frozen electrons). Singlet asymptotes use the retained RHF references.
The QZ TS source used two workers/2500 MB per worker; that historical resource
setting is not relabelled as the angular campaign's four-core allocation.
The angular DFT contract retains RIJONX/def2-J, DefGrid3, ExtremeSCF and central
numerical derivatives (0.001-bohr increment, unscaled frequencies). Native
Arkane's `frequencyScaleFactor/1.014` ZPE convention stays explicit rather
than being replaced by raw harmonic ZPE. Exact inputs and source methods are
retained, not inferred from the method names alone.

All prior electronic-state/count-one/count-two conditional branches remain;
`selected_branch` is null. No new rate is generated by freezing them. The
original first bundle is retained at
`runs/oh-methane-candidate-freeze-20261017`, with its honestly partial timing
coverage. Neither bundle contains new chemistry.

From the repository using the existing application environment:

```sh
python -B tests/freeze_oh_methane_candidate.py reopen runs/oh-methane-candidate-freeze-20261017-complete
```

This reopens the final bundle without dispatch. Its status remains
`diagnostic_not_admitted`.

## Next transfer benchmark: preparation, not dispatch

Order: H + CH4 → H2 + CH3, then OH + H2 → H2O + H. These are known-reaction
rate tests, not autonomous discovery or fate tests. The existing method
procedure is the candidate for transfer; OH/methane's species-specific
symmetries, basin counts and numerical barriers are not transferable values.

Reference access is a real readiness gap. The H/methane study has a
[published 1992 correction](https://pubs.acs.org/doi/10.1021/j100181a600),
listed in the [publisher issue contents](https://pubs.acs.org/toc/jpchax/96/2).
Its effect on the original numbers remains unaudited. The
[OH/H2 institutional abstract](https://esdpubs.nasa.gov/content/Rate_Constant_for_the_Reaction_of_OH_with_H2_between_200_and_480_K)
provides author fits, not an independently audited table of measured rows.
Neither fit-generated values nor secondary database tables are silently
substituted for experiments. Exact pressures, bath/composition, state
preparation, uncertainty semantics and relevant corrections require full
primary documents.

`benchmarks/kinetics/transfer_reference_curation_20261017/curation.json` records
both dossiers, explicit empty measured-row sets, primary/institutional URLs
and precise gaps. For OH/H2, experimental H2 rotational/nuclear-spin
preparation and the OH state model must also be matched. Prediction requests
and scoring evidence remain separate. No reference-fit choice grants a
physical uncertainty interval or a new support claim.

The preparation request is
`benchmarks/kinetics/rate_transfer_20261017/request.json`. Its actual generated
package is `runs/rate-transfer-preparation-20261017`, containing separate
prediction requests, scoring dossier and manifest. Status:
**`prepared_not_authorized`**, `execution_ready=false`. An identical actual
preparation replay preserves file modification times and launches zero jobs.

The proposed budget for each reaction is six optimization/path/recovery jobs,
four full NumFreq, two IRC jobs and five matched electronic single points,
including both asymptotes and the TS. Four cores, six hours/job, 1000 MB/rank,
4 GiB reserve and aggregate 24 native wall-hours / 96 allocated core-hours per
reaction are proposals only, not permission or runtime predictions. Internal
workers consume their owning job's resources; separately invoked screening
and failures consume attempts. A resource preflight must still establish that
the proposal is adequate for the chosen calculation.

The minimal known-reaction preparer has no native dispatcher. Wiring/testing
the existing route verifier and composite/state evaluators under the new
known-reaction contract remains an explicit next software step, not hidden
behind the phrase “ready to run.” Experimental values and their curation
binding stay outside prediction requests; condition/pressure/bath/state and
unit mismatches refuse descriptive scores. The autonomous fate allowlist is
unchanged. Shared candidate evidence is verified once per bundle construction,
not redundantly once per case, with validation again on reopen.

Exact preparation/reopen command (no chemistry):

```sh
python -B tests/prepare_rate_transfer_benchmark.py benchmarks/kinetics/rate_transfer_20261017/request.json runs/rate-transfer-preparation-20261017
```

**Next user decision:** provide lawful full copies/access for Rabinowitz 1991,
its 1992 correction and Orkin 2006; then approve or revise the separate proposed
budget and authorize the bounded known-reaction executor after exact conditions
and case-specific physical/counting assumptions are frozen. No native launch
command is fabricated for an executor that is not yet wired. There is no need
to restart OH/methane or buy more angular Hessians while those prerequisites
are unresolved.

No transfer calculation, production integration, rate reinjection or P1 fate
rerun occurred. This closes the present investigation and leaves a reproducible
next-step package with precise remaining prerequisites; it does not claim the
next quantitative benchmark is already executable or scientifically validated.
