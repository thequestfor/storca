# Bounded rate validation, then return to chemical fate

Latest status (2026-10-17): the bounded angular campaign has closed
`profile_complete_partition_unsupported` at 13 optimization attempts/4 NumFreq.
The selected +170° Hessian and labelled seam closure pass the unchanged gates;
no quantitative angular model is admitted. OH/methane tuning is stopped; the
conditional candidate and preparation-only transfer package are frozen. The
known-reaction executor, full primary-reference/condition audits and separate
compute authorization remain outstanding before transfer execution.
The earlier status below records historical progress. See
[the closeout](OH_METHANE_CLOSEOUT_2026-10-17.md).

Status: **step 1 is a conditional diagnostic; step 2's center/−15° physical
checks passed. A matched-MPI +15° restart converged at cycle 35 (2026-10-13),
but a log-reader defect stopped admission. Its append-only frequency-only
continuation passed all three constrained point/coupling checks. The local
pilot is complete; full step 2 and steps 3–5 remain incomplete**.
See the [electronic-state replay checkpoint](OH_METHANE_ELECTRONIC_STATE_REPLAY_2026-10-11.md).
The [angular pilot checkpoint](OH_METHANE_ANGULAR_PILOT_2026-10-11.md) records
its bounded native campaign and unresolved physical claims.
The remaining step-2 implementation is specified in the
[periodic-profile plan](OH_METHANE_PERIODIC_PROFILE_PLAN_2026-10-13.md): six
wider nodes plus an independent closure control, leaving at most four
refinement/retry attempts and one Hessian under the original cumulative caps.
This is the next focused increment of the
[stability roadmap](../../STABILITY_ROADMAP.md), not a replacement for its
multi-environment objective. No chemistry jobs, production defaults, support
levels or scientific certificates change through this document.

## Decision

Stop treating OH + methane as an unlimited search for settings that reproduce
one experiment. Complete two bounded physical-model checks, freeze a candidate
recipe including its unresolved assumptions, and test transfer to other
reactions. Then return to the complete discovery → verification → rate →
repaired-network → observable workflow.

An unsuccessful physical check is useful evidence. It must identify what is
unresolved, not trigger another unbounded functional/basis/tunneling sweep.

## Starting evidence and the actual remaining problem

The [completed rate audit](OH_METHANE_CONSISTENT_RATE_AUDIT_2026-10-11.md)
establishes the following narrow result:

- The mapped OH + CH4 → H2O + CH3 route is connected by a bidirectional IRC
  on the declared WB97X-D3/def2-TZVPP surface. A scoped numerical-Hessian
  admission is supported by repeats and independent local energy/gradient
  curvature; the analytical weak-mode discrepancy is not explained away.
- All five matched canonical CCSD(T)/def2-QZVPP energies completed. Native
  Arkane statistical mechanics, manual TST arithmetic and detailed balance
  agree. Energy/frequency geometry, state and source bindings reopen correctly.
- At 340.3 K, the composite direct Eckart rate is 3.0623 times the frozen
  measurement. Across the fourteen measured temperatures it is 2.92–3.27
  times the measurement. Those rows are one development study, not fourteen
  independent validations.
- The original software ingestion defects and internally inconsistent source
  combinations have concrete repairs. The remaining discrepancy is not a
  demonstrated RMG/Arkane arithmetic bug. Nor does it prove the electronic
  energy is correct: CCSD(T), a finite basis and a matched calculation still
  have physical limitations.

Two unresolved statistical-mechanical questions now merit investigation:
OH electronic levels/energy reference, and the TS's soft angular motion and
associated basin/counting convention. The native C1/optical-count warning and
76 cm⁻¹ mode cannot be resolved by choosing whichever multiplier matches data.
Classical TST's closer agreement is not a reason to omit tunneling.

The measured humid-air finite-pressure experiment and isolated canonical gas
model are not identical environments. Preserve that applicability warning.
This audit predicts an elementary rate, not a methane lifetime or broad
reaction-network completeness.

## Ordered execution and gates

```text
Retained source/recipe freeze
          ↓
OH electronic-state replay (zero ORCA jobs)
          ↓
Bounded TS angular/basin audit
          ↓
Candidate recipe freeze — unresolved assumptions remain visible
          ↓
Condition-locked transfer tests: H + CH4, then OH + H2
          ↓
Scoped production handoff only if the physical/evidence gates pass
          ↓
Existing cyclobutene P1 end-to-end fate regression
          ↓
Resume wider roadmap with environment-specific validation
```

Reference curation can proceed alongside the first two steps. Execution and
source hashes must be frozen before predictions are compared with references.

### 1. Replay electronic-state thermodynamics consistently

**Question:** How do a declared OH level model and energy reference change the
rate when populations, both barriers and tunneling are recomputed together?

Add a small deterministic electronic-partition evaluator and a new diagnostic
driver. Build on `CompositeRateRecipe` and the retained native modal partition
functions; do not modify raw ORCA logs, bound historical generators or installed
RMG/Arkane. Represent level energies, degeneracies, energy-zero offsets and the
spin-only baseline explicitly. Avoid double-counting the native spin factor.

Evaluate a source-supported OH fine-structure model, with any assumption of
quenched TS spin–orbit coupling labelled as an assumption. A spin-free centroid
and a ground-level reference must yield identical physical predictions when
transformed consistently. Changing an asymptotic energy changes forward/reverse
E0 barriers and Eckart transmission; it is not a posthoc population multiplier.
The [OH/methane theoretical benchmark](https://www2.sci.u-szeged.hu/czako/papers/PCCP_22_14560_2020.pdf)
provides a physical motivation, not a transferable fitted correction.

Replay classical TST, Wigner and Eckart at all retained temperatures. Keep
TS-count branches together with their angular-model assumptions, not as
independent tunable factors. The purpose is consistency and attribution;
agreement with the experiment is not the admission test.

**Gate:** baseline identity, analytic population limits, energy-zero invariance,
dimensional consistency and forward/reverse detailed balance pass; all source
bindings reopen; **zero ORCA dispatches**. Result remains a conditional state
diagnostic until the relevant electronic-state assumptions are reviewed. A
missing TS SOC calculation is not silently declared negligible or calibrated.

### 2. Audit TS angular motion and counting, with a finite campaign

**Question:** Does the soft motion describe local harmonic basins, an accessible
angular domain, or coupled motion for which a one-coordinate correction is
inadequate? Are mirror configurations separately counted without duplication?

Use an experiment-specific audit driver, not a new general conformer engine.
`storca/conformer_rates.py` currently serves bounded neutral-singlet solution
reactant conformers; it is not an OH/radical TS-basin implementation. GOAT or
ordinary minimum optimization cannot establish a constrained reaction saddle.

Before dispatch, bind the exact mapped transferring hydrogen and OH hydrogen,
allowed spectator-H permutations, numerical Hessian, method and coordinate
definition. Inspect the mode and define a nonsingular angular chart, periodic
domain and kinetic metric; do not assume a dihedral remains well-defined near
a collinear H-transfer geometry.

Pilot three constrained saddle points, included in a maximum of nine initial
points. Track the intended unstable reaction mode while fixing the transverse
coordinate. Require native convergence, electronic-state continuity and the
appropriate constraint-projected saddle character. Failed points remain in
the record; silently relaxing to a minimum is not a successful scan.

If the pilot is defensible, complete the periodic profile and permit one
refinement of at most eight added points: **at most seventeen constrained
saddle attempts and four selected numerical-Hessian jobs in this increment**.
Retries consume the cap. Selected Hessians assess orthogonal-mode changes,
coupling and consistency of energy/ZPE conventions, rather than assuming every
other frequency stays fixed. Hash-matched completed jobs are reused.

A reduced angular partition model is permitted only after its metric,
integration domain, orthogonal contribution and state counting are explicit.
Recover the harmonic limit. Replace the corresponding harmonic contribution
when appropriate; do not add a rotor partition on top of that vibration or
multiply again by an optical count for basins already integrated. Profile
geometry alone does not prove a physically accurate rate correction.

**Gate:** domain coverage, mapping/counting and state/derivative checks pass;
the inferred partition contribution changes by less than 10% across the
evaluated temperatures after refinement. This is a numerical convergence
criterion, **not a 10% physical error bound**. Significant coupling, chart
failure, discontinuous states or unconverged coverage yields `unresolved`.
There is no automatic extension to a multidimensional scan.

Preserve the existing IRC as connectivity evidence. Its coarse near-TS spacing
and associated-complex endpoints are not sufficient for quantitative
variational/multidimensional tunneling work. Such work would need a separately
validated denser local path and its own justified campaign.

### 3. Freeze the candidate; stop OH/methane development

Freeze a versioned candidate recipe whether it succeeds, fails or remains
conditional. Its manifest records:

- Geometry/frequency and energy methods, open-shell reference, core treatment,
  grids, derivative admission and the native ZPE convention.
- Electronic levels/energy zero, rotational symmetry, channel/optical counting,
  angular treatment, molecularity, units and standard-state convention.
- Rate theory, temperature/environment domain, applicability exclusions,
  unresolved assumptions, source snapshots and native runtime/module hashes.

Preserve `composite_rate_recipe.py`, the retained helpers and historical
contracts as evidence. Use append-only freeze/diagnostic sidecars. Capture the
participating source versions so future routine edits do not make historical
results unreopenable. Reuse the existing exact parser-revision receipt; do not
rewrite its original hash or admit arbitrary revisions. This needs a small
reproducibility bundle, not a new universal artifact-management platform.

**One optional control, not a new search:** at most three matched TZ electronic
points for OH, CH4 and TS, only if the preceding analysis shows unresolved basis
sensitivity could change the decision. Use the same coordinates/reference/core
model as the five QZ points. This is a forward-barrier sensitivity control,
not a complete alternative Eckart recipe; reverse-barrier consistency requires
products too. Do not infer a calibrated error bar or claim an automatic CBS
limit from two finite-basis results. This control must be declared before
dispatch and does not authorize further basis/reference sweeps.

Give native jobs generous explicitly recorded runtime/memory/scratch budgets;
do not reinstate a ten-minute IRC deadline. Bound this investigation primarily
by its scientific questions and attempt counts. A time/resource stop preserves
partial evidence and returns unresolved, not a physical conclusion.

**Stop rule:** after these checks and, if justified, the one control, no more
OH/methane method tuning in this increment. A frozen conditional recipe may be
evaluated diagnostically; freezing does not promote it to scientifically
admitted or suitable for every reaction.

### 4. Test transfer, with rates separate from autonomous fate

Use two distinct protocols:

1. **Known-reaction rate validation:** supply the declared elementary reaction
   and calculate its verified route/rate without supplying experimental rate,
   barrier or a fitted adjustment. This tests the method, not discovery.
2. **End-to-end fate validation:** supply only system/conditions/observable;
   reference products, mechanisms and kinetics enter the scoring stage only.
   This tests discovery, allocation, verification and propagation together.

Do not weaken `scientific_panel.py`'s reference-free fate request allowlist to
inject expected products. Add a narrowly typed known-reaction protocol if
needed, while retaining the current autonomous protocol. Reuse existing route
verification and the source-bound recipe evaluator, not a parallel production
workflow. The current rate-only overlay is descriptive development scoring;
it is not a locked panel acceptance gate.

| Case | Observable and proposed domain | Evidential role |
| --- | --- | --- |
| H + CH4 → H2 + CH3 | Elementary k(T), exact measured rows near 1000/1200/1500 K | First transfer test; ordinary doublet abstraction without OH fine structure |
| OH + H2 → H2O + H | Elementary k(T) and temperature dependence near 298/375/450 K | Transfer of the OH state and H-transfer model; not a tunneling-free control |
| Cyclobutene → 1,3-butadiene | Existing P1: provisional k comparator; derived first-order parent-loss quantities at 423.15 K/50 Torr | Already exposed singlet engineering control for a future end-to-end regression; not a full-fate reference or a new holdout |

Primary reference leads are [Rabinowitz et al.](https://pubs.acs.org/doi/10.1021/j100155a033)
(direct H/methane measurements over 897–1729 K),
[Orkin et al.](https://pubs.acs.org/doi/10.1021/jp057035b)
(OH/H2 measurements and non-Arrhenius temperature dependence), and the
[existing P1 manifest](../../benchmarks/stability/p1_cyclobutene_development.json).
H/methane full measurement rows, pressure, uncertainty and secondary-chemistry
treatment remain curation-pending; a publisher title is not a numerical
comparator. Audit OH/H2 primary rows, bath conditions and uncertainty semantics
before locking them too. Use an author fit only when explicitly labelled as
such. The proposed temperatures are selections to curate, not invented data.

Cyclobutene is not a physics-trivial control: primary work treats its
[heavy-atom tunneling and anharmonic kinetics](https://pubs.acs.org/doi/10.1021/acs.jctc.1c01143).
Its reference transcription still needs visual audit. Derived first-order t95
is not an independent lifetime measurement. The current scientific panel has
no audited full-fate comparator; running P1 end to end does not manufacture one.
H + H2 exchange can remain a
software/path-counting control; without labelling it is not a useful net-fate
observable. OH/methane isotope tests are optional mechanism diagnostics, not
additional independent reaction families.

**Freeze before predictions:** exact reference slice and units, conditions,
applicability policy, method/diagnostic escalation rules, resource ceilings,
scoring rules and source hashes. This is a condition-locked development/transfer
panel, not a blind prospective holdout: the sources and cyclobutene result have
already been inspected. A future genuine holdout needs a separately locked
protocol and withheld scoring data.

Suggested engineering accuracy targets to preregister for declared-applicable
cases: each reaction's geometric-mean multiplicative error ≤3, with no scored
point exceeding a factor of 10. Report signed and absolute log errors,
temperature-dependence error and the worst point as well. These are usefulness
targets for engineering triage, not experimental uncertainty, a calibrated
confidence interval or sufficient grounds for a supported-lane claim.
If source/applicability curation changes the targets, do so before execution
and record the reason; do not redefine a pass after seeing predictions.

All frozen cases stay in the denominator: absent output, crash, budget stop,
outside-scope and unresolved are different outcomes, not discarded rows.
Passing an applicable subset does not pass the panel. Report declared and
completed coverage separately from accuracy, and count false applicability
decisions. Multiple temperatures from one study are correlated; do not let
them dominate independent-reaction evidence. Prespecify one summary per
reaction for panel aggregation. Add source-independent applicability challenge
fixtures (for example, an explicitly declared spin-changing or condensed-phase
request sent to this gas/single-surface profile): no unsupported rate or false
reassurance is permitted. Withholding is a correct challenge outcome, but is
not an accurate prediction on an in-domain experimental case.

Compare source-bound native RMG estimates, available frozen legacy-STORCA
outputs and the candidate recipe. Missing legacy predictions are missing
baselines, not made-up estimates. Bind RMG database provenance: experimental
information may already be present there, so this is not necessarily a
reference-blind RMG baseline. Record actual wall/core time, failed/reused work
and high-level calculation counts under comparable declared budgets.

**Gate:** transfer accuracy, applicability, counting and reproducibility pass
for the declared narrow domain without tuning to each new answer. A failure
opens a named physical assumption and a new bounded plan; it does not authorize
automatic retuning. Any subsequent method change creates a new frozen profile
and requires reevaluation of the entire cohort, including failures, not just
the bad case; exact-compatible QM artifacts may still be reused. Two transfer
reactions cannot calibrate general uncertainty
or prove the method across chemical classes/environments.

### 5. Integrate only admitted evidence, then close a fate loop

When justified, adapt the frozen recipe to existing
`route_thermochemistry.py`, `arkane_runner.py`, `verification_adapters.py` and
`generated_kinetics.py`. Preserve native statmech receipts, exact state/mapping
and condition identity, route evidence, reversible consistency and replacement
validation. Do not inject the isolated audit's rates merely because execution
completed or a comparator was approached. Diagnostic records stay
`diagnostic_not_admitted`; the existing generated-library gate should refuse
them.

Use `method_applicability.py` for explicit diagnostic facts and scope decisions.
Use `kinetic_uncertainty.py` for source-backed, condition-bound uncertainty when
available. Neither a factor-three engineering target nor a three-case error
spread is a calibrated rate interval. An uncertain method may support a
labelled conditional estimate without supporting a conclusive lifetime.

Resume the existing P1 discovery → mapped TS/IRC → rate → repaired library →
repropagation workflow. Require the requested observable, verification
coverage, applicable rate uncertainty and model assumptions to satisfy the
existing convergence policy. Compare resumed and fresh propagation with the
same model. Do not force a converged stability verdict when coverage or
applicability is unresolved.

**Gate:** one documented end-to-end result whose rate and observable claims
can be traced to source-bound calculations and explicit assumptions, including
an honest unresolved verdict when evidence is insufficient. Passing this
regression is not proof of general pathway completeness.

## First five implementation work items

These are local work packages, not GitHub issues created by this plan.

1. **Freeze retained source versions and add state-aware replay.** Small pure
   evaluator, new diagnostic driver and append-only manifest. Test baseline
   identity, reference invariance, level populations, detailed balance, tamper
   refusal and zero ORCA dispatch. No production integration. Depends on the
   completed rate audit, not additional QM.
2. **Add the bounded angular-saddle audit.** Mapped coordinate/constraint
   contract, native-evidence readers, capped owned execution and conditional
   angular partition report. Test wrong atom/state/constraint refusal, saddle
   inventory, closure/counting, harmonic limit, partial completion and exact
   reuse. Requires item 1's joint counting model. No general rotor/conformer or
   variational-dynamics platform.
3. **Curate and lock the rate-transfer protocol.** Primary-reference audit,
   separate known-reaction request type, recipe/version freeze, explicit
   acceptance and all-outcome/cost accounting. Reuse the scientific-panel
   primitives while retaining autonomous reference exclusion. Test leakage,
   changed references/recipes, missing results and failed-job accounting.
   Curation can start now; predictions wait for item 2's freeze.
4. **Run H/methane, then OH/H2 without answer-driven adjustments.** Retain
   direct rates, route/state/counting evidence and costs; evaluate all selected
   temperatures. Acceptance is the preregistered transfer gate, not one good
   point. No methane lifetime without a specified reactive environment/network.
5. **Admit the narrow recipe and rerun P1 only if justified.** Production
   adapter plus exact-evidence and reversibility regression tests; then one
   end-to-end fate assessment. Requires item 4 and applicable uncertainty/
   coverage evidence. If those are missing, retain conditional outputs and
   next-evidence requests instead of bypassing the gate.

For every item, close out with source identity, tests, actual scientific
evidence, unresolved limitations and the next bounded decision. Software
completion and physical validation are separate checklist entries.

## What this does not authorize or abandon

No new package collection, default-method switch, generic AIMD/VTST/RPMD
integration, blind automatic rotor correction, open-ended basis search or
replacement of RMG/ORCA/Arkane. No confidence claim inferred from experimental
agreement. No LLM authority for physical-model admission.

All environments remain the product objective. These tests establish, at most,
one narrow gas-reaction rate lane. Ordinary solution chemistry, neat liquids,
solids, surfaces, light and other environments remain on the existing roadmap
and require their own evidence. Reuse the validated orchestration/provenance
machinery, not a gas-phase accuracy certificate. The existing solution
rate-only benchmark can guide a separate solvent-specific recipe gate; it is
not evidence that bulk-solution fate already works.

The immediate deliverable is **a transferable, auditable rate recipe—or a
specific diagnosis of why it is not transferable—followed by a real fate-loop
test**, not a perfectly fitted OH/methane number.
