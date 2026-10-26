# Comprehensive, Working STORCA Stability: Implementation Plan

Recorded: 2026-09-09. Status: accepted direction; implementation incomplete.

## Purpose and user intent

Build an autonomous, evidence-driven chemical-fate system that answers, for a
specified chemical system, environment, conditions, observable, and budget:
what chemistry occurs, through which mechanisms, producing what products, on
what timescale, and with what justified uncertainty.

The user wants comprehensive environment coverage AND working predictions,
not exclusively gas-phase chemistry, a collection of package adapters, or a
framework whose scientific inputs must all be supplied manually. Breadth must
grow through complete, benchmarked applications. No finite discovery search
proves universal chemical completeness; report the supported domain and gaps.

This document is the execution plan to reference in a future `/goal` request.
It refines the execution order of `STABILITY_ROADMAP.md`; it does not erase its
requirements or imply that existing planned/deferred capabilities are done.
Read it together with `REQUIREMENTS.md`, `TRACEABILITY.md`,
`DEFINITION_OF_SUPPORTED.md`, `METHODS.md`, `LIMITATIONS.md`, and
`benchmarks/stability/milestones.json`. Paths in this document are relative to
the repository root unless the context clearly identifies a sibling document.
If documents disagree, resolve and record the discrepancy before claiming
completion. Do not silently discard a requirement.

## Prediction-first execution order — 2026-10-01

This section is the current execution order. The P1–P7 requirements and dated
history below remain authoritative; completed repairs must not be repeated.
This is a planning revision, not an implementation or support promotion.

### What we must deliver

The product is a stability ladder/campaign, not a collection of chemistry jobs.
For each requested environment it must report the observable, relevant
products/mechanisms, model assumptions, uncertainty, evidence, and unresolved
requirements. A molecular local minimum is not a lifetime. A verified path is
not a validated rate. A populated lane registry is not environment support.

Retain three different outcomes: conditional model estimates; experimentally
benchmarked predictions within a declared applicability domain; and explicit
unresolved/outside-applicability results. Reuse existing claim/status fields
where possible. An estimated central rate cannot populate a withheld verified
lifetime field or complete the active P1 goal. No output is a laboratory safety
certification or a universal claim that a molecule is stable.

The central development objective is useful, accurate conditional predictions
at bounded cost, with low false reassurance AND useful prediction coverage.
An engine that abstains on every molecule does not meet it. Nor does an engine
that confidently returns numbers for unsupported physics.

### Current baseline and the missing producers

| Boundary | Reusable machinery | Evidence/capability still needed |
|---|---|---|
| Discovery | `reaction_candidates.py`, `reaction_discovery.py`, `discovery_assembly.py`, bounded RMG and graph/xTB; schema-v3 solvent-bound ALPB screens and narrow neutral carbon four-pi proposals | Consequential missed-channel recovery; non-RMG/earlier-only candidates must reach verification even without RMG flux; independent held-out coverage assessment |
| Path | `reaction_path.py`, `route_geometry.py`, `surface_resolution.py`, TS/frequency/IRC and source-bound atom maps | Real selected-case path/rate/repair demonstration; conformer coverage; auditable method applicability beyond a single surface |
| Rates | `arkane_runner.py`, solution Eyring producer with bounded reactant refinement/conditional scenarios, `specialist_rates.py`, `mess_adapter.py` | Applicable, automatically produced rate-error evidence; real conformer-rate benchmarking and specialist physics inputs, not another interface |
| Allocation | `flux_verification.py` already runs two-sided t95 rate perturbations and can rank sensitivity × supplied uncertainty / supplied cost | Credible uncertainty/cost producers; horizon-censored and low-flux initiation handling; candidate influence before a rate exists |
| Solution | `solution_path_verification.py`, verification session, explicit condition-bound thermochemistry-temperature production, bounded reactant sampler/native refinement, linear and nonlinear propagators; in-progress native-candidate conditional forward/reverse coefficient producer | Autonomous discovery; automatic applicable rate-error evidence; real non-298 K/conformer-rate validation, interconversion and TS/encounter ensembles, diffusion and admitted bimolecular production/propagation |
| Validation | Scientific-panel curation, coverage-risk manifest, specialist benchmark harness | Frozen numerical references, locked prospective holdouts, real predictions and calibration results |

The completed 2026-10-07 stage-reuse checkpoint passes 1090 application tests in
110.136 seconds (22 expected native-only skips), 125 relevant native/discovery/
stage tests in 6.542 seconds and 130 affected tests in 2.198 seconds. Actual
preparation and source-metadata negatives launch zero chemistry; they do not
verify IRC, rate, propagation, calibration or coverage. See the stage-reuse
evidence below. The preceding completed 2026-10-06 application checkpoint passes 1018 tests in
151.461 seconds (22 expected native-only skips); the expanded native/discovery
boundary suite passes 115 tests in 9.695 seconds, including the previous 68
focused native checks. The affected suite passes 179 in 1.234 seconds. The
preceding coefficient checkpoint passed 999 in 157.101 seconds and 68 native
tests in 16.729 seconds. The IRC/reference checkpoint passed
976 in 107.129 seconds and 68 native tests in 9.725 seconds; the integrated-conformer
checkpoint passed 972, the recipe checkpoint 948 and the intermediate conformer
checkpoint 970. The latest logs cover solvent-bound screening, environment-bound
cache identity and the strict solution-phase vocabulary repair; see
[screening evidence](../../runs/solution-discovery-screen-diagnostic-001/EVIDENCE_2026-10-06.md).
The preceding [coefficient evidence](../../runs/solution-rate-coefficients-diagnostic-001/EVIDENCE_2026-10-05.md)
remains a separate historical checkpoint.
Exact artifact ownership,
downhill execution-retry/source-policy admission, the descriptive rate-only overlay,
requested-temperature/recipe matching, bounded reactant refinement and stale
session-production-policy refusal are software-tested. Installed-RDKit seeds
are deterministic; the new cache check preserves 804 retained native files.
The read-only ordinary-basename IRC replay preserves 102 native files and the
same 119-frame trace; execution-retry and source-mutation cases are intercepted
regressions, not physical retry/full-chain equivalence. The JPL limiting-rate
source audit corrects f(298) to f(300) in eight evaluated anchors while preserving
central rates, conditions, partitions and scoring thresholds. It supplies no
STORCA error calibration, finite-pressure experiment or locked holdout.
Actual ORCA conformer refinement/rate benchmarking is still unproved. This is
not chemical prediction validation. Cyclobutene
discovery, a real parent minimum and validated mapped endpoints exist. The
corrected native NEB performed 22 completed iterations before its 600-second
job limit, without a verified saddle or rate. The v9 geometry boundary remains
active. Exact retained checkpoint admission is software-tested and replayed;
controlled interruption and the full prediction chain remain open. Dispatch
007 subsequently timed out before a completed iteration, with no own native
band; its result remains null. The current panel has ten pending curation
entries, one source-audited descriptive high-pressure-limiting-rate comparator,
and the provisional P1 comparator. It has zero audited full-fate comparators and
zero prospective holdouts; the full panel is not frozen. The JPL H + O2 reference
is evaluated/theory-informed, not an independent measured finite-pressure rate
or STORCA calibration. Its actual offline outcome remains missing; see
`runs/scientific-panel-jpl-hpl-002/report.json`. A separate descriptive overlay registers eight CD3CN fitted
rates without inventing missing concentration histories; actual predictions are
eight missing/zero scored and it does not promote full-case reference readiness.
The coverage-risk case list is
empty. Arkane explicitly reports `calibration_required` after central-rate
generation. These facts, not the number of implemented modules, set priorities.

2026-10-05 in-progress solution slice (items 9/10): a retained solvent-bound
native candidate can enter the existing path producer and yield conditional
forward/reverse one-/two-molecule coefficients without a species table,
concentrations, duration or fabricated uncertainty bound. Explicit matched
thermochemistry, separated-species sums and 1 M dimensional conversion preserve
shared-TS detailed balance. Packets carry unresolved error/ensembles/transport
and no network-admission authority; no stability bundle is emitted. This is
an internal candidate-to-calculation API and existing CLI extension, not a manual
rate-import replacement product, autonomous discovery or validated solution
prediction. The completed 999-test application and 68-test native checkpoint
covers this initial software slice, not a real positive solution coefficient
benchmark or supported solution fate. Applicable automatic error evidence, diffusion,
coverage, real positive benchmarks and autonomous propagation remain required.

2026-10-06 in-progress item 9/shared-discovery slice: requested SMD/CPCM candidate
conditions are now separate from observed native GFN2-ALPB/`gsolv` endpoint
screens, with no unsupported-solvent gas fallback. Bounded mapped,
valence-preserving isolated-carbon cyclobutene opening and acyclic diene closure
produce neutral hypotheses without supplying an expected product to discovery.
The three-case water/acetonitrile development batch recovered its three expected
products using six actual one-core xTB endpoint jobs in 0.924418 seconds; exact
cache checks preserved 24 files per case and attempted zero new endpoint jobs.
Native convergence/connectivity, method/solvent/JSON/stdout and source hashes
are checked. Schema-v3 endpoint and environment-bound assembly contracts refuse
mismatch in place; old gas artifacts remain historical and are never relabeled.
The shared path gate now recognizes `phase="solution"` only with the same
required named-solvent, SMD/CPCM and exact method-keyword checks; this repairs a
reproduced native handoff mismatch, not physics applicability.
No ORCA path, rate, lifetime, coverage or calibration is established. Full
autonomous solution selection/verification/propagation, applicable error
evidence and the original items 1–10 remain open. See
[screening evidence](../../runs/solution-discovery-screen-diagnostic-001/EVIDENCE_2026-10-06.md).

2026-10-07 evidence update: a separate source-audited **primary measured** OH +
methane rate slice now retains all fourteen original points from one study,
selecting the first table row before prediction. Its actual bounded evaluator
outcome is one missing prediction and zero scored rates. Unknown point-error
semantics, shared systematics and finite-pressure/HPL applicability are explicit;
this is not a calibrated error profile, full-fate comparator or holdout. The
shared eleven-case manifest plus optional P1 remains unchanged. See
[primary-source evidence](../../benchmarks/stability/reference_sources/2026-10-07-berg2024-oh-methane/EVIDENCE_2026-10-07.md).

The evaluator also preserves separately valid measured cost when prediction
metadata is malformed: an actual P1-003-derived corruption replay retains
73.12736991699785 seconds while withholding scientific scores. Missing/corrupt
cost evidence remains unknown, not zero. The combined reference/accounting
suite passes 54 tests in 1.363 seconds, not a new full application checkpoint.
See [accounting evidence](../../runs/scientific-panel-audit-001/EVIDENCE_2026-10-07.md).

An actual installed-Arkane partition audit infers two optical configurations for
the C2 TS candidate, against the producer's forced one; both achiral endpoints
remain one. With identical unscaled modes this gives partition-only forward and
reverse factors of two and unchanged equilibrium. This is a counting-policy
discrepancy, not permission to double a rate or infer a physical uncertainty.
Resolve configurational/channel inventory and degeneracy before changing or
calibrating the recipe. No quantum job, complete IRC or prediction follows;
see [native diagnostic](../../runs/arkane-optical-partition-diagnostic-003/EVIDENCE_2026-10-07.md).

An additional item-5 admission correction requires explicit string condition
binding for `method_calibration` records and prevents specialist/planner ID
stringification from creating false matches. Actual public-consumer synthetic
before/after replays and 126 affected tests (root replay 1.075 s) support this
narrow safety fix, not source/recipe/domain authentication or an automatic error
producer. Old unbound calibration artifacts remain preserved but refuse new
allocation; user/database supplied scenarios are not promoted. See
[admission evidence](../../runs/method-calibration-admission-diagnostic-001/EVIDENCE_2026-10-07.md).

### Track A — close a real, reusable prediction chain

**Next implementation batch:** continue the frozen cyclobutene development case
under the current endpoint policy, reusing `p1-cyclobutene-parent-002` and valid
source-bound evidence. Inspect the actual retained resume contract before
choosing continuation versus a fresh derived dispatch; never edit an old
scientific artifact to make resume work.

1. Check current process ownership, evidence hashes and aggregate budget.
2. Exercise endpoints → NEB/TS → frequencies/IRC → applicable Arkane central
   rate → library insertion → isothermal repaired propagation. Read each real
   boundary, including replacement direction and temperature/pressure evidence.
3. Reproduce and repair only observed engineering failures. If an applicable
   path/rate cannot be found within the recorded bounds, retain the reason;
   do not force a saddle, alter the scoring reference, or restart indefinitely.
4. Demonstrate controlled interruption/resume through rates and repair, not
   merely a restart-file parser. Start with interruption after a completed path
   and before rate/repair; then test an in-flight route interruption. Verify
   immutable evidence and which expensive jobs were actually reused.
5. Bind the resulting model calculation to independently justified rate-error
   and discovery-coverage evidence from Track B before claiming P1 success.

The historical 2026-10-03 reconciliation of terminal production attempts
001–007 consumes 1321.1618874169653 seconds and leaves
1078.8381125830347 seconds of the original ceiling. The source-hashed read-only
reconciliation is `runs/p1-full-validation-2026-10-03/production-budget-001.json`.
Dispatch 007 is terminal, not a planned/live job. Recompute from actual terminal
records before dispatch; subsequent attempts and resume share the remainder,
not another fresh per-dispatch ceiling. Preserve four aggregate cores, finite per-job
ceilings and the original cleanup grace. If insufficient, report measured
remaining work and obtain a new budget; this plan does not grant unbounded
computation. One of the three P1 candidate preflights has been used.

2026-10-07 update: attempt 008 is terminal after a converged OptTS, independent
one-imaginary-mode/mapped-coordinate check and a 600-second IRC timeout. There
is no completed bidirectional connection, Arkane rate, repair or lifetime.
The actual eight-attempt ledger consumes 2037.3441255009966 monotonic seconds
and leaves 362.6558744990034. Source-bound completed-stage reuse now prepares a
new bidirectional IRC from 008's geometry/Hessian without repeating OptTS/frequency
or borrowing partial-IRC credit. Actual preparation preserves 940 original plus
198 retained-stage files. Per-job ceilings and scientific conditions remain
unchanged; no additional compute has been authorized. This is not full-chain
resume equivalence or a P1/first-ten closure. See
[reuse evidence](../../runs/p1-ts-frequency-reuse-diagnostic-001/EVIDENCE_2026-10-07.md).

**Exit evidence:** real path/rate/repaired-network artifacts; unchanged frozen
conditions/reference; the preregistered factor-three development rate criterion
and 1% fresh/resume numerical criterion; applicable uncertainty and coverage;
simple fresh/resume commands. The factor-three scoring threshold is not a
physical confidence interval. A chain demonstration with missing uncertainty
remains an engineering milestone, not completion of P1.

### Track B — establish scientific accuracy while repairing the chain

Do not wait for P1 to finish before curating references and developing rate-
error evidence. Do not expand the active goal into unrelated environment code;
the wider panel defines subsequent work and supplies relevant P1 evidence.

**B1. Freeze and run the chemical panel.** Extend the existing manifests and
`tests/run_p1_*` drivers into a small bounded batch evaluator, not a second
workflow engine. Start with 8–12 development cases spanning ordinary thermal
decomposition, competition, radicals, pressure dependence/barrierless chemistry,
conformer sensitivity, neutral solution and deliberate inapplicability. Retain
the existing difficult cases as failures/challenges, not forced positive cases.

Each scored case needs actual physical conditions, observable definition, rate
units/convention, quantitative reference and provenance, known uncertainty and
unknowns, applicability assessment, budgets and thresholds before predictions.
Observed products are not assumed microscopic-mechanism ground truth. A lifetime
derived analytically from a measured rate is labeled as a derived comparator.
Keep reference-only products/rates inaccessible to discovery and calculation
requests. Report pending references separately; do not quietly remove failed
predictions from an eligible scored cohort.

The existing CPD solution source specifies a hydrocarbon fraction, not pure
hexane; the saponification source requires ionic chemistry outside the neutral
lane. Neither becomes a ready positive solution benchmark by changing labels.
Review or replace their reference readiness while preserving the original entries.

Run RMG-only where applicable, frozen current STORCA, and each proposed change
under matched resource ceilings. Report log-rate and lifetime errors, retention
classification, branching/dominant-product accuracy where evidenced, influential
mechanism recall, justified/false abstention, false reassurance, interval/bound
performance, executed high-level jobs and measured core/wall time. Define false
reassurance separately for an incorrect observable and a consequential omitted
path. Preserve unresolved and crashed cases in outcome accounting.

Lock prospective holdouts before tuning; cases already inspected during
development are not holdouts. A small pilot or zero observed errors cannot
establish universal calibration. Preregister domain-specific sample sizes,
accuracy, prediction-coverage and cost thresholds before formal evaluation;
report sampling uncertainty and per-stratum results.

**B2. Produce credible rates, rather than require users to invent bounds.**
Introduce versioned applicability/calibration profiles consumed by the existing
rate-uncertainty contracts. Profiles bind reaction class, elements/charge/spin,
environment, temperature range, electronic structure, conformer treatment,
statistical mechanics and tunneling/pressure treatment to source-backed error
evidence. Missing or out-of-domain calibration remains unresolved.

Develop these profiles from separate, relevant rate/barrier/thermochemistry
datasets and validate prospectively. Do not fit a bound to cyclobutene's answer
and call the same case validation. Agreement between two DFT methods measures
method sensitivity, not calibrated accuracy. Comparison to higher-level theory
also has applicability/error limits. Keep deterministic sensitivity ranges,
empirical prediction intervals and physical upper bounds distinct.

Complete STAB-009/010 as integrated rate work: bounded endpoint and encounter
conformer searches, identity-preserving deduplication, selective ORCA refinement,
TS searches for relevant conformers, and temperature-matched thermochemistry.
Use population-weighted rates only when conformer interconversion justifies
equilibrium; otherwise represent conformer states and interconversion rates
explicitly. Treat hindered rotors/low modes, symmetry and degeneracy consistently;
avoid counting the same entropy contribution twice. Bound the cost and retain
unsampled-conformer risk. No new JSON-only conformer feature earns this milestone.

Propagate correlated species/TS energy errors and rate-model alternatives through
the existing repaired network. Preserve forward/reverse thermodynamic
consistency. Initially report conditional bounds/scenarios; add probabilistic
confidence only when the error model and reference uncertainty support it.

**B3. Make method selection and discovery improve the observable.**
Add a deterministic, evidence-backed method policy around current surface/rate
gates. Check declarations and available diagnostics for environment mismatch,
spin change, barriers/loose association, collisional stabilization, H transfer/
tunneling and conformer sensitivity. Absence of a warning is not proof of
single-reference or single-surface applicability. Missing diagnostics mean
explicit conditional assumptions or escalation. No LLM decides physical validity.

Use the existing sensitivity/uncertainty/cost ranking before inventing a new
optimizer. Extend beyond local t95 only where benchmarks demand it: concentration
at the requested horizon when t95 is censored, product/branching influence,
low-flux initiation and correlated rate perturbations. Give unresolved candidates
and discovery gaps a separate exploration allocation; zero/absent RMG flux is
not evidence of irrelevance. Apply physical upper bounds only within their
documented scope. Refresh influence after each repaired propagation.

Evaluate one physics-driven discovery backend on bounded missed-path and false-
positive cases. SCINE Chemoton remains the preferred first feasibility study,
not a mandated rewrite or installation: test Apple-Silicon dependency/runtime
feasibility, restart/budget ownership and export of atom/path/electronic-state
evidence before broad integration. It has its own database/scheduler machinery;
STORCA should import discovery evidence, not adopt a competing product architecture.
If feasibility fails, retain the result and evaluate one bounded alternative,
not several parallel integrations. Agreement among engines does not prove
coverage. Screening rejection/absence needs its own false-negative benchmark.

**Exit evidence:** a declared domain with useful prediction coverage and
condition-matched accuracy, honestly validated uncertainty, and recovery or
exposure of consequential missing chemistry at recorded cost. All-refusal and
candidate-count-only improvements fail this gate.

### Track C — grow the ladder through complete environment applications

Once shared path/rate handoff is reliable, autonomous homogeneous solution is
the next implementation track, alongside discovery/calibration work—not after
solving every gas reaction. Preserve the broader environment target throughout.

**C1. Ordinary neutral solution.** Connect native, solvent-bound discovery to
the existing solution candidate/session machinery, so identity plus conditions
does not require a manually assembled reaction network or barrier-error bound.
The narrow solvent-bound graph/ALPB producer supplies an initial discovery
boundary, not the complete autonomous loop or an independently calibrated
mechanism-coverage domain.
Do not transfer gas rates or gas absence claims into solution. Extend the
single-temperature unimolecular producer to scientifically benchmarked condition-matched thermochemistry, then
neutral bimolecular TS/rate production with explicit standard states, molecularity,
dimensional factors, encounter evidence and diffusion diagnostics. Reuse the
existing nonlinear reversible mass-action propagator. Select a source-backed
neutral rearrangement and a neutral dimerization in a specified solvent before
predictions; the current mixed-solvent CPD entry is a research lead, not an
already admitted case. Independently vary concentration/temperature within
the benchmark domain and score decay/product evidence. pH/speciation, ions,
solvent-as-reactant mechanisms and nonideality remain named follow-on milestones.

**C2. A complete specialist-rate producer.** Reuse Arkane's existing capabilities
for conventional pressure-dependent wells where admissible. For a selected
loose/barrierless case, finish a physically sourced capture/microcanonical-flux
producer plus the existing MESS handoff; MESS does not supply missing capture
physics. Validate detailed balance, bath/energy-transfer inputs, falloff on an
exact T/P grid and independent reference cases. Preserve supplied-result support
but do not substitute it for autonomous production. A VaReCoF/ORCA potential
bridge is an explicit feasibility/research dependency, not a presumed capability.

**C3. Additional matter/exposure lanes.** Each row below is a complete scientific
application with its own reference-data and method-evaluation gate. Its numerical
backend is selected after that evaluation, not installed by this plan.

| Lane | Required new physics/evidence | First demonstration must distinguish |
|---|---|---|
| Solvent-assisted solution | Solvent participation, speciation, explicit clusters or free-energy sampling; consistent solvent chemical potentials | Passive continuum versus a participating solvent mechanism |
| Neat liquid | Intermolecular discovery, ensemble sampling, density/composition, activities and rate/free-energy treatment | Dilute isolated-solute approximation versus concentration-dependent liquid fate |
| Crystalline solid | Structure/polymorph, relevant defects, periodic reaction/migration barriers and event kinetics | Molecular decomposition versus packing/defect-controlled transformation |
| Surface/interface | Adsorption, coverage/site balance, partitioning and transport with microkinetics | Homogeneous rate versus contact/catalytic/transport-controlled fate |
| Photochemical exposure | Spectrum/intensity, absorption and excited-state branching/yields, thermal competition | Dark versus illuminated fate in an explicitly supported phase |
| Electrochemical interface | Potential-aware states, interfacial/speciation model and electron-transfer kinetics | Potential-dependent fate versus an ordinary neutral ground-state calculation |

No molecule-only SMILES can specify a unique crystal, interface or exposure.
Ask only for physically indispensable inputs. Prefer internal adaptive defaults
for search effort/method escalation; do not expose each new implementation choice
as a user flag. Keep independent lane predictions separate until an explicit
coupled exposure/transport model exists. Disordered solids, phase transitions
and strongly coupled multiphase dynamics require separate later scope.

### First implementation backlog and closure evidence

These are planned work items, not created GitHub issues or completed features.
They map onto existing requirements; do not build a second engine or another
registration-only lane to satisfy them. Rows can overlap where dependencies
permit. The expanded 2026-10-01 objective authorizes these ten work items; P1
remains the first acceptance gate, and autonomous solution follows reliable
shared prediction boundaries. Later matter/exposure milestones are not erased.

| Order / work item | Existing anchors and work | Closure evidence / dependencies |
|---|---|---|
| 1. Real selected-case path/rate/repair | `run_p1_workflow.py`, `reaction_path.py`, Arkane and repair adapters; correct budget and reuse before actual dispatch | Real artifacts and numerical comparison; final P1 closure also needs items 2/4/5/7 |
| 2. Frozen panel and bounded evaluator | Existing scientific/coverage manifests and P1 drivers; references isolated from predictions, outcome/cost ledger | Frozen scored development cohort, prospective holdout protocol, retained baseline failures; starts now alongside item 1 |
| 3. Full-chain interruption/resume | `verification_resume.py`, engine/ORCA leases, stability report checkpoint; reproduce any initial-interruption gap first | Equivalent rates/propagation within frozen tolerance and audited job reuse; depends on item 1 reaching those boundaries |
| 4. Auditable method applicability | Existing environment/surface/rate gates; versioned diagnostic requirements and decisions | Conventional and deliberately inapplicable cases routed correctly, missing diagnostics explicit; starts with items 1/2 |
| 5. Automatic applicable rate-error evidence | `kinetic_uncertainty.py`, Arkane result and solution barrier contracts; domain-bound sourced profiles, separate calibration/evaluation | No fabricated default error factor; out-of-domain refusal; real propagated bounds and prospective validation; depends on items 2/4 and profile-specific item 6 |
| 6. Conformer/temperature-aware rates | STAB-009/010; bounded sampler → selected ORCA refinement → rates, conditional equilibrium or explicit conformer network | Flexible versus rigid cases, temperature checks, deduplication/detailed balance and measured integrated cost; no schema-only completion |
| 7. Prediction-conditioned native queue/allocation | Native discoveries, earlier-only supplemental candidates, `flux_verification.py` and engine | Consequential absent/low-flux route investigated; irrelevant expensive work reduced without accuracy loss; depends on items 2/5 and reliable item 1 |
| 8. One physics-discovery feasibility/benefit trial | `ReactionDiscoveryBackend`, native mapping/evidence admission; first evaluate Chemoton | Budgeted restartable evidence export and consequential held-out discovery gain, or explicit feasibility rejection; depends on items 2/4/7 |
| 9. Autonomous ordinary solution loop | Existing solution candidate/session/path/kinetics modules; native solvent-bound graph/ALPB discovery slice; automatic applicable error evidence still needed | Identity/solvent/conditions → quantitative neutral unimolecular benchmark without supplied pathways; depends on shared items 1/4/5/6/7 |
| 10. Bimolecular solution production | Extend solution TS/free-energy rate producer and reuse mass-action propagation | Concentration-dependent measured benchmark, standard-state/unit/diffusion/reversibility checks and useful uncertainty; depends on items 2/9 |

For every work item, add a regression that fails before a demonstrated defect
fix, replay the affected real boundary, retain hashes and failures, and record
the user-visible capability separately from tests. A feature item needs both
positive scientific evidence and deliberate negative cases. Do not erase
missing-physics failures or tune references to obtain a passing dashboard.

### Tool and resource decisions

Keep the existing ORCA, RMG/Arkane, xTB and native propagators. CREST is a
candidate bounded conformer sampler, not a replacement rate theory. Prefer
existing RDKit/xTB sampling for the first small integrated slice if CREST setup
is a blocker; assess whether its additional sampling changes the benchmark.
Chemoton is an evaluated discovery adapter, not a mandatory dependency.
MESS/capture software is added only for a complete specialist application.
No mass installation, new service, paid compute or external message is
authorized by this planning revision.

Four local cores can develop and exercise small molecular cases; measured
calculations determine viable budgets. Broad conformer, capture, condensed-phase
and periodic studies may require other compute. Record that resource dependency
instead of quietly truncating the physics or promising completion dates.
Prefer dependency gates to calendar promises: reliable shared loop first;
panel/calibration in parallel; autonomous solution next; specialist/liquid/solid
applications afterward according to benchmark need and resources.

Primary documentation consulted for the backend choices (2026-10-01):

- [SCINE Chemoton](https://scine.ethz.ch/download/chemoton): AFIR/NT exploration,
  conformer generation, filtering and MongoDB/Puffin execution. This does not
  establish local installation readiness or STORCA compatibility.
- [CREST](https://crest-lab.github.io/crest-docs/): bounded conformer/ensemble
  sampling tools. [Entropy limitations](https://crest-lab.github.io/crest-docs/page/examples/entropy.html)
  explicitly caution against general noncovalent-system entropy use.
- [Arkane input documentation](https://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/input.html):
  conventional/statistical-mechanical and pressure-dependent inputs, including
  rotors and collisional networks. Existing program capability is not proof that
  STORCA automatically produces every required physical input.

## Starting evidence, not a completion claim

This is the historical planning snapshot, not the latest implementation status.
For a future goal, read the dated checkpoints below before choosing work; verify
their artifacts and current tests rather than repeating already completed fixes.

At planning time:

- Gas and restricted implicit-solution campaign orchestration exist, but broad
  scientific validation and autonomous solution discovery remain incomplete.
- Native candidates, RMG discovery, bounded graph/xTB enumeration and endpoint
  screening exist. xTB endpoint energies are not barriers or reaction rates.
- ORCA endpoint/path/TS/frequency/IRC machinery, Arkane integration, repaired
  network propagation, uncertainty/allocation slices, and resume contracts exist.
  Reuse them; inspect actual interfaces rather than redesigning from filenames.
- Specialist-rate contracts and a limited MESS adapter exist. They are not a
  complete calibrated capture/master-equation production pipeline.
- The NNONC retained run found an unexpected O-N-cleavage/H-transfer channel.
  Its proposed products and atom correspondence are retained, not admitted as a
  verified reaction or lifetime. Queue handoff and retained-TS reuse are unfinished.
- The last full passing checkpoint had 632 tests. A subsequent retained-fragment
  validation change has a failing test: `KeyError: 'N'` during seed reuse in
  `_stationary_fragment_records`. Inspect the reused `labels` variable before
  fixing it; do not report the current tree as passing from an earlier run.
- Conformer ensembles, common adaptive method selection, and cross-environment
  experimental calibration are not complete. Neat liquids and bulk solids need
  explicit milestones rather than being subsumed under solution/surface labels.

Revalidate this snapshot against the worktree and retained artifacts before work.
Never restart a calculation solely because an observation timed out. Confirm
the actual process/session state and reuse validated checkpoints.

## Scientific and implementation invariants

- Preserve bounded searches, explicit states, route identity, atom/isotope maps,
  TS and IRC validation, repaired-model propagation, and scientific provenance.
- Distinguish not found, unresolved, estimated, verified within a model, and
  outside applicability. Successful program execution is not scientific success.
- Keep discovery, path verification, rate theory, and propagation separable.
- Keep RMG and ORCA where useful; do not rewrite functioning machinery to make
  the architecture look cleaner. Do not assume DFT is ground truth.
- Never copy a rate from a failed proposed reaction to a newly observed channel.
- No universal lifetime across environments without an explicit exposure model.
- No LLM authority over physical-model selection. Use diagnostics, algorithms,
  evidence, and explicit assumptions.
- Keep kinetics as the macroscopic-time bridge. Dynamics and enhanced sampling
  are optional mechanism/free-energy tools, not a brute-force lifetime solution.
- Do not prioritize user-supplied pathway import. The requested workflow is
  autonomous discovery. Existing specialist evidence support is not permission
  to substitute manual reaction assembly for autonomous operation.
- Avoid new user flags unless scientifically necessary information cannot be
  obtained from the declared system, conditions, installed configuration, or
  retained evidence. Preserve the normal environment workflow and `--rmg-env`.
- Preserve unrelated changes and append-only evidence. No destructive cleanup,
  paid service, or expansive installation is implied by this plan.

## Sequence and acceptance gates

### P1 — Close and stabilize the end-to-end prediction loop

Priority: immediate. Roadmap links: M1, M2, STAB-011, applicable M6 rate gates.

Work:

1. Fix and regression-test the retained-fragment reuse bug.
2. Validate unexpected endpoint fragments as minima with identity and explicit
   electronic-state assumptions/diagnostics; geometry alone is insufficient.
3. Queue unexpected native candidates automatically, preserving parent evidence,
   exact conditions, deduplication, bounded depth, and compute accounting.
4. Reuse retained TS/IRC only when atom correspondence, route identity, method,
   environment, convergence, and evidence hashes authorize that reuse.
5. Complete applicable rate generation, repaired-network insertion, propagation,
   uncertainty treatment, and observable reporting. Unsupported rate theory
   remains a named unresolved requirement.
6. Exercise fresh and interrupted/resumed operation through the full chain,
   including rate and repaired-network stages, without manual artifact editing.

Existing anchors: `reaction_path.py`, `route_geometry.py`,
`intermediate_identity.py`, `verification_engine.py`, `verification_resume.py`,
`arkane_runner.py`, and `verification_adapters.py` under `storca/`.

Gate: a suitable fresh case and its resumed equivalent reach equivalent
scientific outcomes, without duplicate expensive calculations or weakened
validation. The NNONC case remains a bounded regression, not the only success
criterion or an obligation to force a particular product/lifetime.

### P2 — Establish a permanent scientific test panel

Priority: begin alongside P1. Roadmap link: M11 begins now, not after features.

Curate approximately 8–12 source-backed systems spanning ordinary gas
decomposition, competing paths, radicals, barrierless/pressure-dependent
chemistry, conformational sensitivity, homogeneous solution, solvent-assisted
chemistry, and deliberately inappropriate electronic structure.

For each case, record reference conditions, rate/lifetime/product/mechanism
evidence, experimental uncertainty, applicability, and provenance. Freeze
acceptance thresholds before inspecting predictions. Separate development from
held-out cases; do not fill missing references with invented numbers.

Compare RMG-only where applicable, current STORCA, and successive improvements.
Where RMG is not applicable, use an appropriate lane baseline and say why.
Measure observable error, dominant-path accuracy, mechanism recall/precision
where reference knowledge allows, uncertainty calibration, false reassurance,
abstention, wall time, and high-level calculation count. Retain failures in the
denominator. Experimental agreement alone does not prove a mechanism.

Gate: demonstrate useful predictions within a declared domain and appropriate
rejection outside it. Publish reproducible commands, resources, artifacts, and
limitations. Synthetic tests and live chemistry benchmarks remain distinct.

### P3 — Add genuinely independent physics-driven discovery

Priority: after candidate handoff can turn additional discoveries into rates;
evaluate early against P2. Roadmap link: M2.

Retain three roles: RMG for known reaction patterns, native graph enumeration
for bounded alternatives, and one physics-driven backend for exploration
without requiring predetermined products.

Preferred first evaluation: SCINE Chemoton, as a discovery adapter rather than
the definition of the STORCA network or a replacement orchestration system.
This is an evaluation recommendation, not a proven platform-compatibility claim
or an unconditional dependency commitment.

Run a bounded feasibility/chemistry comparison before substantial integration:
installation/platform suitability, atom/path evidence export, reproducible
budgeted resume, important missed-channel recovery, false positives, and cost.
Compare its incremental value with the existing graph/xTB backend. If it fails
the gate, record that result and evaluate one alternative; do not accumulate
adapters. Biased searches produce candidates, not unbiased physical rates.

Integrate accepted discoveries through native identities, provenance, screening,
verification, and coverage accounting. Independent engines may share errors;
agreement alone is not calibrated completeness.

Gate: recover consequential missed chemistry on held-out cases and improve a
prediction or correctly expose an earlier unjustified prediction, at recorded
cost. More candidate count by itself is not success.

### P4 — Observable-driven allocation and credible rate evidence

Priority: alongside P3; necessary before substantially increasing search size.
Roadmap links: M2, M3, M6, M11.

Extend existing `flux_verification.py`, `kinetic_uncertainty.py`, and the
verification engine rather than duplicating their existing functionality.
Measure observable response to uncertain rates across the requested time span.
Account for transient intermediates and low-flux initiation routes. Rank work
by influence, reducible uncertainty, expected diagnostic success, and measured
cost; retain auditable score components rather than insisting on one formula.

Add bounded conformer sampling and temperature-consistent thermochemistry for
influential species/routes. Preserve forward/reverse consistency and shared
species-energy correlations. Keep mechanism, state, solvent, electronic-model,
conformer, kinetic-theory, and numerical uncertainties distinguishable. Do not
invent confidence levels for uncalibrated bounds.

Gate: comparable or improved held-out prediction quality at lower expensive
calculation cost than flux-only allocation. Convergence is conditional on the
observable tolerance, coverage assessment, and applicability; budget exhaustion
is reported as budget exhaustion.

### P5 — Make homogeneous solution autonomous

Priority: an active expansion track alongside P3/P4, not after all gas chemistry.
Roadmap links: M4 and M7.

Build on the existing solution candidate, verification-session, mass-action,
observable, and campaign modules. Supply autonomous condition-aware discovery,
consistent solvent/standard-state rates, relevant speciation assumptions,
reversible bimolecular kinetics, and diffusion diagnostics. Start with a
declared tractable chemical domain, then widen it through benchmarks.

Introduce bounded explicit-solvent participation where necessary, including
placement/cluster-size checks and appropriate solvent chemical-potential
accounting. Escalate when cluster/continuum evidence cannot justify a bulk
solution claim. Do not confuse implicit solvent with autonomous solution chemistry.

Gate: starting from chemical identity and solution conditions, produce a
benchmark-supported concentration, yield, or lifetime without a user-supplied
preverified network. Missing speciation or solvent participation is visible.

### P6 — Expand environments through complete applications

Priority: milestone-level expansion after reusable discovery/rate contracts
work; define benchmarks and missing requirements early. Roadmap links: M7–M10
plus new explicit neat-liquid and bulk-solid milestones.

| Domain | Required complete scientific slice |
|---|---|
| Solvent-assisted solution | Explicit solvent mechanisms, appropriate free energies and solution kinetics |
| Neat liquids | Composition/activity treatment, intermolecular mechanisms and adequate sampling |
| Gas–liquid / liquid–solid contact | Partitioning, contact geometry, transport and coupled networks |
| Crystalline solids | Declared structure/polymorph, relevant defects, reaction/migration processes |
| Surfaces | Adsorption, site balance, coverage-dependent reactions and microkinetics |
| Electrochemistry | Potential-aware states/rates and justified interface assumptions |
| Photochemistry | Light spectrum/exposure, excitation, branching/yield evidence and thermal competition |

Register support only after a complete restricted application passes its gate.
Solid structure cannot generally be inferred from a molecular SMILES; request
necessary physical inputs explicitly. Amorphous/disordered solids and strongly
coupled phase changes require separate scope and later validation. No claim of
"all environments" may hide these gaps.

Gate per domain: one independently benchmarked end-to-end application, including
discovery, verification, rates, propagation, uncertainty, failure modes, and
resume. A declaration or an external-code launcher does not satisfy this gate.

### P7 — Auditable method selection before rate admission

Priority: cross-cutting from P1/P2, expanding with observed failures—not a final
step after all lanes. Roadmap links: M5/M6 and specialist lane gates.

Implement explicit policies for ordinary saddle/TST applicability; loose
association/capture or variational treatment; collisional stabilization/master
equations; tunneling; conformer sensitivity; solvent-assisted/free-energy
physics; spin changes; multireference warnings; excited states; periodic and
potential-aware processes. Separate automatic diagnostics, required physical
declarations, and specialized external methods. Lack of a diagnostic is not
evidence of applicability.

Finish a useful benchmarked specialist rate producer before adding more
provider interfaces. Preserve existing Arkane/MESS/evidence contracts where
appropriate. Escalation may end in unresolved, not automatic misuse of TST.

Gate: preregistered ordinary and difficult cases receive the justified treatment
or an intelligible refusal. Measure missed and unnecessary escalations.

## Execution discipline and progress reporting

Every batch must specify its user-visible capability, exact tests/benchmark,
compute budget, stopping condition, and remaining unsupported scope before
expensive work. Determine budgets from available resources and measured jobs;
do not promise an unsupported wall-clock schedule. Use the retained four-core
contract for applicable existing runs; inspect current resources before new work.

A stubborn molecule receives a bounded diagnostic budget. At its limit, retain
the failure and continue the wider panel. Do not spend unbounded time forcing
one case to succeed, drop it from evaluation, or call abstention a prediction.

Maintain two independent records: software verification and scientific evidence.
Update this plan's ledger plus the roadmap, requirement catalog, traceability,
manifest, METHODS, and LIMITATIONS when scope or support claims change. Keep
prior evidence and test history; never replace a failing current result with an
older passing count. Reconcile already implemented work before writing new code.

### Progress ledger

| Work package | Planning-time state | Evidence required to close |
|---|---|---|
| P1 | in progress; unexpected-channel path and central Arkane rate retained; verified repair enters native RMG core; restart mapping and shared source-bound identity proof implemented; saved NNONC identity replay correctly withholds stereo admission; propagation/uncertainty gates remain open | Fresh/resume full-loop results, route/rate/propagation artifacts |
| P2 | nine development entries curated in scientific_panel.json; aqueous and hydrocarbon-solution rate measurements transcribed; panel not frozen | Sourced panel, locked thresholds, complete baseline/failure report |
| P3 | graph/xTB exists; external evaluation not started | Backend feasibility and held-out incremental chemistry/observable results |
| P4 | allocation/uncertainty slices exist | Correlated credible-rate treatment and allocation comparison |
| P5 | restricted solution infrastructure exists | Autonomous solution benchmark without preverified network |
| P6 | mostly planned/declared; explicit bulk milestones missing | Per-domain complete applications and support evidence |
| P7 | individual safeguards exist; common policy incomplete | Auditable policy and method challenge-set results |

## Instructions for a future implementation goal

Execution checkpoint (2026-09-09): repaired seed reuse by keeping XYZ element
labels separate from reaction-species labels in `_stationary_fragment_records`.
The retained-coordinate regression exercises preparation, repeat preparation,
stubbed minimum validation and mutated-seed rejection. This is not new DFT
validation. Full unittest discovery passed 632 tests in 68.597 seconds after
the repair. `benchmarks/stability/scientific_panel.json` starts P2 with source
locations and explicit missing evidence. The aqueous saponification case uses
charged reactants and is deliberately not misrepresented as supported by the
current neutral-solution lane. No baseline prediction or holdout claim exists.

Unexpected-channel handoff checkpoint: the engine now extracts proposals from
its callback evidence, checks identity/reactant/environment/hash contracts,
deduplicates, bounds one discovery generation, and records/persists the new
unverified route. A fixture-based engine run selects the parent and then the
unexpected channel, while refusing any rate generation from either unverified
path. Changed environments, reactants, hashes, unresolved mapping and depth
exhaustion are rejected. Retained-TS admission and its resume validation remain
next; the queue change itself does not close P1 or validate new live chemistry.
The full suite subsequently passed 632 tests in 80.314 seconds. A separate
diagnostic on the actual saved NNONC proposal passed queue admission and saved
`unexpected-queue-admission-diagnostic-001.json` beside the assembled IRC;
no chemistry ran and the old engine state was not changed.

TS-source checkpoint: new unexpected proposals carry a checked, hashed source
bundle for optimization, frequency/Hessian and IRC evidence. Source method,
charge/spin, frequency geometry, TS presence in the path and completed IRC
branches are checked without executing chemistry. Tests reject changed methods,
states, geometries and abnormal frequency termination. The actual NNONC source
passed and is saved append-only as `retained-ts-source-001.json`. This establishes
source evidence only; the new reaction's mode, stationary products, electronic
state applicability and rate admission remain open. Do not promote this to a
verified new reaction or repeat a complete NEB search to work around that gap.
Verification: full suite passed 633 tests in 56.414 seconds; after adding the
nonfinite-frequency guard, the seven IRC/source tests passed again. No new
electronic-structure calculation or end-to-end observable was produced.

Route-specific assessment checkpoint: `assess_unexpected_retained_path` now
checks current source hashes and method/state evidence, exact native identity,
trace/endpoint correspondence, reconstructed bond changes and imaginary-mode
displacement. The real NNONC new-channel assessment passes and is stored under
`new-channel-assessment-001/assessment.json`. Seventeen focused tests pass,
including changed-hash and changed-bond rejection. The diagnostic is not a
verified path or rate. Bounded live isolated-species validation has started in
`runs/unexpected-channel-validation-reaction-2cdf1d1a16c24c3a-001` with four
cores, B3LYP/def2-SVP/NoAutoStart, at most eight calculations, a 1800-second
surface wall budget, two core-hours, and 600 seconds per requested job (clipped
by remaining budget). Inspect the live handle or terminal result before any
restart. Completion and scientific validity are not yet established.

Subsequent live result: session 66044 completed with `stationary-results.json`
status `validated`. All three separated species passed optimization,
frequency-minimum and connectivity checks at the retained method. This is
conditional minimum evidence, not experimental validation or proof of the
assumed singlet/bond-order description. The default ORCA adapter now wires
queued proposals into generic retained-path verification with the ordinary
environment/state/minimum/thermochemistry gates and standard Arkane inputs.
Mocked integration tests prohibit new TS/IRC/NEB execution and reject failed
minima. Live rate generation and network repair remain open; do not rerun the
completed diagnostic or claim the newly wired path was executed end-to-end.
Live source checking exposed an ORCA origin shift between Hessian and TS XYZ.
Translation-only alignment now records that shift and requires unchanged
internal geometry/axes (live residual about 1.3e-8 angstrom). Twelve focused
tests pass after this change; rotated/permuted Hessians are not silently accepted.

Read-only stationary reuse checkpoint (2026-09-10): the generic verifier can
reuse completed fragment folders after checking retained seeds, method/state
metadata, normal completion, minima and connectivity, retaining source hashes.
No output is copied and no electronic job is launched. The actual new-channel
run now has `verified_transition_state_path` and standard Arkane inputs in
`runs/unexpected-channel-rate-handoff-reaction-2cdf1d1a16c24c3a-001/route-verification.json`.
Input generation succeeded in `rate-preparation-001` without rate execution.
The optional thermochemical loss bound reports `unsupported_phase` for the
low-pressure surrogate spelling; it was not used as negative evidence. Inspect
and fix that applicability inconsistency separately. Before launching live
Arkane, implement bounded execution: `src.orca_runner.run_arkane` and its
syntax/validation subprocesses currently lack timeouts. Do not substitute input
preparation for a verified rate, pressure model, uncertainty, or repaired network.

Bounded-rate execution checkpoint (2026-09-10): Arkane execution and all three
RMG/Arkane syntax/library helper subprocesses now use a finite 600-second
per-command ceiling. POSIX timeout cleanup kills the process group (including
conda descendants), retains text streams, and does not retry or admit a library.
This is a safety ceiling, not yet shared campaign-budget accounting. The
low-pressure gas surrogate spelling is admitted by the existing conditional
thermochemical bound with pressure and phase retained; non-gas/unknown models
remain rejected. This does not establish a pressure-dependent rate or validate
the collision-ceiling assumptions. The focused run passed 172 tests; after
adding explicit timeout-to-failed-route coverage, 11 execution/thermochemistry
tests passed. Full-suite verification then passed 639 tests in 55.456 seconds.
Bounded live input/API validation completed in `rmg4_env` (RMG 4.0.0), recognizing
four StatMech jobs, a Kinetics job, and a PressureDependence job for the prepared
new-channel input. Session 70755 exited zero; this was input loading, not job
execution or rate validation. No new electronic-structure or rate calculation
ran in this checkpoint. Next: bounded live Arkane execution using `rmg4_env`,
followed by rate applicability and repaired-network gates; P1 remains incomplete.

Live rate checkpoint (2026-09-10): `rate-execution-001` failed during RMG
equilibrium-weight exponentiation because uncorrected absolute electronic
energies were very large and negative. A process-local log-weight experiment
in `rate-execution-002` removed that overflow but exposed the same reference
problem in density-of-states Boltzmann factors. That workaround was removed
from production rather than extended into multiple upstream patches; RMG's
installed environment was not modified.

The input generator now applies one common per-atom electronic-energy gauge
using the complete reactant channel as reference. Atom balance preserves
channel energy differences and TS barriers for both 1-to-2 and 2-to-1 routes;
the three-part `e_electronic` tuple keeps Arkane responsible for ZPE. Each
species provenance records the gauge. These energies must not be used as
absolute formation thermochemistry or combined across unrelated route gauges.
Generated-expression tests cover both orientations. Ten Arkane tests passed,
then full discovery passed 639 tests in 64.145 seconds.

Live `rate-execution-003` completed Arkane's rate and pressure-dependence jobs
(normal termination at 2026-09-10 14:24:05 local; about 344 MB reported). STORCA
then failed export because this Arkane job wrote no native species dictionary.
A bounded RMG graph-serialization fallback is now under test, with explicit
spin/roundtrip/duplicate-label checks and no rate or thermo estimation. Validation
is being resumed from the completed output; do not rerun the rate calculation.
The export change's focused/full verification and live admission result remain
to be recorded. No lifetime, uncertainty calibration, or P1 completion is claimed.

Subsequent export result: 11 focused Arkane tests passed. Validation-only session
36173 completed without rerunning kinetics, producing `library-validation-001.json`
and a route-unique `verified_RMG_library/.../generated-kinetics.json` beneath
`rate-execution-003`. RMG reloaded the Chebyshev library and reproduced the
condition rate, approximately 3.24e-7 s^-1 at 298.15 K and 1e-6 bar, for the
new channel (not the failed homolysis). This is a central-model diagnostic,
not a lifetime. Its manifest explicitly retains `calibration_required` rate
uncertainty and unvalidated collisional estimates. The original execution
failure record remains untouched; the later validation is separate evidence.
Next: validate the new-candidate rate's repaired-network handoff and uncertainty
gates without allowing unknown uncertainty to imply convergence, and exercise
fresh/resumed runs. The full suite after the export fix is running; record its
result before asserting the current tree passes.
Live numerical cross-check: the high-pressure Arrhenius fit in
`rate-execution-001/output.py` and `rate-execution-003/output.py` has identical
printed A, n, Ea, grid and fit-error values after the gauge correction. This is
agreement at retained output precision, not a pressure-grid convergence study
or a new physical uncertainty estimate.
Final verification for this checkpoint: full unittest discovery passed 640
tests in 63.288 seconds after the dictionary fallback; `git diff --check`
passed. No chemistry process remains running from this checkpoint. Continue
with repaired-network insertion/propagation and the wider scientific panel,
not another calculation of the already completed path or central rate.

Network-handoff checkpoint (2026-09-10): 37 verification-engine tests pass,
including a new case where an initially missing native reaction is verified,
inserted and repropagated, but its unknown Arkane rate uncertainty blocks a
converged observable. The live full-RMG repair is running under
`network-repair-001` alongside the retained rate artifacts, with four requested
cores, one attempt, a three-minute RMG wall budget, 240-second RMG hard timeout,
and a 600-second whole-diagnostic ceiling. Session 28057 is the observation
handle; inspect it before any restart. The original path/rate jobs are not rerun.

P2 adds cyclopentadiene dimerization in a hydrocarbon fraction (Szekeres et al.,
1977), with four prose-reported rate constants and explicit unresolved solvent,
pressure and rate-convention details. The current neutral-solution validation
gap remains: a hydrocarbon fraction is not a pure-solvent input. No experimental
uncertainty was invented, and no predictions have been scored against this case.

Live repair result: session 28057 ended normally as an orchestration process,
but `network-repair-001/repair-result.json` reports incomplete RMG generation
and rejected repair. RMG encountered the recognized DASPK initial-step
singularity (two resurrection attempts). The library did load: RMG.log records
the exact new reaction and products being added to the edge. They are absent
from the retained core mechanism, so applied-rate matching correctly fails.
Independent propagation and flux closure are unavailable; no lifetime is
reported. Do not conflate requested library loading with a core rate actually
used in propagation. Before another live attempt, address how verified native
channels are retained for propagation (seed/core admission versus an explicitly
justified negligible-edge bound), and use the existing bounded DASPK recovery
policy where appropriate. Do not repeat ORCA/Arkane or increase budgets blindly.
Full unittest discovery passed 641 tests in 67.233 seconds after this
checkpoint; panel JSON parses and `git diff --check` passes. The diagnostic
and test sessions have ended. P1 remains incomplete; the wider P2–P7 scope is
unchanged.

Native core-admission checkpoint (2026-09-11): validated generated libraries
marked as kinetic repairs now enter RMG through `seedMechanisms`, not only the
edge. Ordinary reference reaction libraries stay on the edge; duplicate paths
are removed and conflicting registered library names are rejected before
execution. Manifest hashes are checked before admission, and insertion evidence
explicitly retains post-propagation rate matching and uncertainty requirements.
Core placement is not a scientific-validity exemption.

The bounded `network-repair-002` failed at startup because RMG 4.0.0's seed
insertion looked up an absolute external-library path in a basename-keyed
registry. A process-local launcher adapter now resolves only loader-registered
aliases, without basename guessing or modifying the installed RMG environment.
`network-repair-003` then successfully inserted the new reaction and its three
species into the initial core. The first attempt produced a retained model, but
reported twelve core collision-rate violations. Its automatic second attempt
failed when loading restart libraries cleared RMG's external alias map. Neither
attempt establishes repaired-network closure or a lifetime. The full software
suite at that stage passed 646 tests in 63.522 seconds.

Restart resolution now retains loader-established aliases bound to the exact
loaded library objects. Missing libraries and changed object identities fail
closed. A read-only check using the actual `rmg4_env` library loader, generated
rate library, and `network-repair-003` saved restart libraries passed: RMG's
alias map was cleared, yet the adapter resolved the original library object.
The diagnostic is `check-restart-library-resolution-001.py` beside the retained
rate handoff and executed no chemistry. All 64 focused seed/engine/wiring tests
passed in 1.090 seconds. Full-suite verification is pending at this checkpoint.
Before any further NNONC simulation, audit the retained model and collision
violations and choose a bounded diagnostic; do not repeat completed ORCA/Arkane
work or treat successful library loading as successful propagation. P1 remains
incomplete, and broader P2–P7 requirements remain unchanged.

Full-suite result after restart alias retention: 647 tests passed in 64.105
seconds. The saved-library check and all test sessions ended; no chemistry was
running. Subsequent read-only audit established that `network-repair-003` first
stopped because not enough of its three-minute wall budget remained for another
enlargement, not DASPK. It retained 37 core species/93 reactions and 518 edge
species/2187 reactions. Its solver profile reached one year, which is time
coverage only, not mechanism convergence. All twelve collision violations were
RMG family estimates (eleven recombination, one disproportionation), not the new
Arkane channel; several exceeded their collision bounds by roughly 1e13.

Rate-serialization checkpoint (2026-09-11): real bridge inspection found that
RMG exported the verified Chebyshev model with a zero minimum pressure
(`PCHEB/ 0.000 9.869 /`) instead of its source 1e-6 bar bound. The reloaded rate
was NaN even though the library had entered the core. The launcher now
re-renders only the gas Chebyshev polynomial block from its in-memory model,
preserving finite positive T/P bounds and coefficients at binary64 text
precision, dimensional conversion and surrounding RMG reaction metadata.
Unknown output layouts fail closed. Bridge inspection now rejects invalid
forward rates with reaction-specific diagnostics; its protocol cannot emit
NaN/Infinity, and invalid reverse/collision evidence stays unavailable.

The bounded serialization-only replay in
`runs/unexpected-channel-rate-handoff-reaction-2cdf1d1a16c24c3a-001/chebyshev-export-check-001/result.json`
passed. It restored only the checksummed source model for the uniquely
forward-isomorphic reaction in an in-memory copy of the 93-reaction core,
wrote a separate derived export, and compared 30 T/P points after reload:
zero observed relative error. All original file hashes remained unchanged.
No propagation, mechanism growth, ORCA or Arkane execution occurred. This is
not an automatically resumed repaired network or a stability prediction.

Native RMG serialization tests separately cover one-, two-, and three-reactant
synthetic polynomials, two polynomial shapes, and SI/Chemkin rate conversion.
All four tests passed in `rmg4_env` (3.201 seconds); the initial three-reactant
fixture needed the reader's fourth-order unit slot for its temporary (+M)
count. Main-environment full discovery passed 657 tests in 65.595 seconds with
that native-only test skipped; six bridge tests and four existing standalone
numerical tests also passed. After the observed-failure manifest update, 19
focused checks passed with the same native-only skip. The final full-suite
result after all changes remains to be recorded.

Next P1 handoff work, before further NNONC growth:

1. Add source-bound, resonance-aware structural replacement matching. RMG
   renamed the native products and retained `[NH-][NH+]=O` instead of the source
   `N=[NH+][O-]` resonance form; label signatures or exact canonical SMILES alone
   cannot establish identity. Preserve orientation, stoichiometry, spin,
   charge/isotopes, unique source-library provenance and rate agreement.
   Existing manifests need an append-only identity sidecar checked against
   their exact libraries, not edits or new quantum calculations. Reuse the
   proven correspondence in both applied-rate checks and planning replacements.
2. Retain parsed evidence per RMG attempt. The controller currently retains
   earlier artifact paths but returns only the final attempt's parsed evidence;
   a failed restart must not erase useful candidate/inspection evidence or
   promote an earlier incomplete attempt to final success.
3. Exercise the automatic fresh/resume handoff on a suitable bounded panel
   case. NNONC remains a failure regression with uncalibrated rate/collision
   assumptions; do not spend unlimited growth/verification budget forcing it.

The observed failure is registered in `observed_failure_profiles.json` and
linked to gas regression tests. No lane or milestone was promoted. Continue P2
curation/baselines and the full P3–P7 scope; export precision is infrastructure
evidence, not the requested comprehensive scientific end state.

Execution checkpoint (2026-09-22): replacement identity and irreversible-library
hardening. This supersedes the earlier "add structural replacement matching"
next step; do not rebuild that implementation or repeat completed quantum work.

- `rmg_identity_bridge.py` compares resonance-aware graphs in the owning RMG
  environment, preserving spin, charge/isotopes, occurrence stoichiometry,
  direction, exact source-library provenance and unique channel matching.
  `reaction_correspondence.py` retains content-addressed sidecars bound to
  source/model file hashes and exact inspection records. Both rate verification
  and verification planning consume the same proof; source records are unchanged.
- The native irreversible test failed before matching because its `=>` entry
  omitted `reversible=False`. The production photon-source writer had the same
  defect. Both declarations are corrected. A real-writer RMG regression checks
  atom balance, irreversibility and preservation of the supplied forward rate.
  Library loading failures now retain the original exception under the specific
  `replacement_source_library_load_failed` refusal, not a generic identity error.
- A separate native handoff regression exports and reloads a genuine renamed,
  achiral one-reaction Chemkin model in both forward and reverse orientations.
  It checks that graph correspondence, source-bound sidecar, exact directional
  rate, applied-rate validation, planning rebind, and reopened JSON evidence
  agree. Its rates and thermochemistry are synthetic serialization fixtures,
  not scientifically validated chemistry.
- The retained replay at
  `runs/unexpected-channel-rate-handoff-reaction-2cdf1d1a16c24c3a-001/identity-replay-001/result.json`
  passes all 11 checks. Its actual 93-reaction inspection finds the expected
  index-0 graph candidate, but potential stereochemistry is not preserved by
  the RMG artifacts. Both consumers refuse verified credit, reopened evidence
  agrees, the same sidecar is reused, and all ten original source/model/contract
  artifact hashes are unchanged. The earlier precision-preserving export is a
  diagnostic derivative, not a completed repaired network.
- No new quantum calculation, rate-model fit, mechanism growth or propagation
  was performed for that replay. No lifetime, calibrated rate, new supported
  environment or P1 completion is claimed. Potential stereo is an identity-
  evidence limitation, not proof of a distinct long-lived stereoisomer.
- Verification (2026-09-23): application-environment full unittest discovery
  passed 693 tests in 72.656 seconds (18 expected native-only skips). Native
  RMG focused identity, actual photolysis-writer, handoff, and Chebyshev export
  suites passed 23 tests in 3.971 seconds. The full suite first exposed one
  outdated planning test that omitted the now-required condition contract;
  its fixture now supplies exact condition and rate evidence and proves missing
  evidence withholds planning credit. `git diff --check` passed. These software
  tests and synthetic rates do not establish scientific benchmark accuracy.

At that checkpoint, next P1 work was: retain parsed evidence for every bounded RMG attempt, without
promoting an incomplete earlier attempt; then exercise automatic fresh/resumed
full-loop operation on a suitable bounded panel case. NNONC stays a regression
with unresolved stereochemical identity and uncalibrated rate/collision inputs,
not an obligation to force a positive prediction. P2–P7 remain in scope.

Execution checkpoint (2026-09-29): the active goal is the focused P1
fresh/resumed quantitative-prediction gate, not the full-roadmap suggested goal
below. Broader milestones remain intact and unrelated integrations are paused.
`WORKING_PIPELINE_GATE.md` now holds the acceptance ledger and exact commands.

The failed-restart reproducer is fixed: the controller persists source-hashed,
condition-bound parsed records before each retry and records setup failures.
Earlier direct candidates are rebuilt through a structural whitelist and
revalidated by discovery, planning and route selection. Their rates, flux,
verification and lifetime credit are not inherited. Successful Chemkin parsing
also survives a subsequent inspection failure as raw, non-admitted evidence.
Legacy artifacts are not retroactively assigned missing attempt proofs.

The saved 93-reaction native replay in `attempt-evidence-replay-001/result.json`
under the existing NNONC handoff passes all 14 checks, including unchanged
hashes for 45 source files and identical reopened consumers. Its later failure
is simulated and its source mechanism is the earlier diagnostic derivative;
it proves evidence transport, not a completed chemical prediction. A completed
final model retains earlier absent discoveries, but automatic supplemental
selection of those earlier-only routes remains future handoff work.

Verification: 727 application tests ran in 79.482 seconds, OK with 18 expected
native-only skips. The focused native identity, photolysis-writer, handoff and
Chebyshev suites passed all 23 tests in 3.277 seconds. No quantum calculation,
rate generation or new mechanism growth ran during this engineering checkpoint.

`benchmarks/stability/p1_cyclobutene_development.json` freezes the first proposed
live case before prediction: cyclobutene at 423.15 K and 50 Torr, one-hour
isothermal dark gas model. The development accuracy and fresh/resume tolerances
are distinct from calibrated uncertainty; current gates are unchanged. The
1958 experimental table still requires visual auditing. The first bounded
preflight has not yet run at this checkpoint. Rate uncertainty, real path/rate
verification, repaired propagation and interruption/resume remain open; P1 is
not complete.

Execution checkpoint (2026-09-30): the first cyclobutene discovery preflight
completed and independently found ring opening. A saved native replay exposed
energy-on Cantera self-heating in an isothermal request. Explicit energy-off
propagation and retained T/P checks repair this; legacy unproven propagation
must be regenerated before flux/rate/lifetime admission. Original artifacts
remain unchanged. HPL Arkane no longer demands an unused bath gas; pressure-
dependent bath requirements remain. The bounded real parent ORCA check passed
in 42.463 seconds with a validated minimum; the first sandboxed launch failed
at MPI socket startup, and both attempts are retained separately. The completed
minimum is eligible for production reuse. Path/rate/repair is next. See
`WORKING_PIPELINE_GATE.md` for exact compute/test evidence. No path/rate or
quantitative-prediction milestone is promoted.
Final application check: 749 tests, 21 expected native-only skips, 57.079 seconds.
The relevant native suite passed 32 tests in 5.041 seconds.

Suggested goal text:

Latest focused P1 checkpoint (2026-10-01): source-bound RMG family atom tracing
repairs the observed cyclobutene graph-edit ambiguity without forcing expected
chemistry. The unmocked cross-environment replay resolves both orientations,
reuses reopened evidence and leaves 46 original files unchanged. Recipe traces
remain discovery evidence; TS/IRC/rate and observable admission are unchanged.
Bridge and RMG-supervisor interruption ownership is corrected after a real
orphan reproducer. See `WORKING_PIPELINE_GATE.md` for exact commands, limits,
test checkpoints and the still-open production path/rate/repair/resume gate.
Do not repeat the completed parent minimum or the earlier attempt-store fix.
Final verification: 767 application tests in 66.200 seconds (22 expected skips),
46 focused native tests in 5.429 seconds, and `git diff --check` passed. The
retained live proof also passes a read-only reopened admission check. No new
path, calculated rate, propagation or chemical-route resume was completed.

The broader suggested goal below is not the active focused P1 objective.

Geometry checkpoint (2026-10-01): the first real production-chain retry reached
RMG discovery and family mapping but stopped before ORCA paths: the distance
gate misclassified nonbonded cyclobutene ring diagonals. Bounded independent
coordinate-graph corroboration repairs the actual seed/parent rejection without
loosening cutoffs. The retained boundary replay preserves all 111 original files.
The user's broader validation concern also reproduced wrong-stereo, bound-IRC
identity, per-endpoint-tolerance and malformed-data gaps. Current v9 policy
separates mapped chemical geometry, coordinate correspondence, stationarity and
electronic applicability; the expanded 20-method audit passes. These are
engineering/geometry results, not scientific prediction validation. See
`GEOMETRY_VALIDATION_AUDIT.md` and the working-pipeline ledger for evidence,
current tests and the remaining aggregate production budget. Resume the real
path/rate/repair chain next, not another parent minimum or more environment
features. No verified lifetime or actual route resume has been demonstrated.
Final software checkpoint: 794 application tests (22 expected native-only
skips), 68 focused native tests and `git diff --check` passed. These counts do
not establish method accuracy or prediction convergence.

> Implement and harden the comprehensive STORCA stability plan in
> docs/stability/COMPREHENSIVE_IMPLEMENTATION_PLAN.md. Treat that file and the
> linked roadmap/requirements as the full scope. Begin by inspecting current
> state and the latest dated checkpoints, then finish outstanding P1 work and
> execute P2. Do not repeat fixes already evidenced as complete. Develop broader
> independent discovery and autonomous solution
> chemistry as the next expansion tracks, with P4/P7 scientific gates throughout.
> Preserve all unfinished environment milestones. Do not substitute infrastructure
> count, passing software tests, or one successful molecule for scientific
> completion. Maintain the progress ledger and evidence, use bounded computation,
> and report real missing authority/resources when encountered. No new flags or
> manual pathway-import workflow unless genuinely required by a physical input.

A goal is not complete merely because P1/P2 finish. For every package and
explicit environment, inspect authoritative acceptance evidence before closing
the full goal. Missing benchmark references, specialist physics, or permissions
remain recorded dependencies, not permission to invent data or shrink scope.

## Sources behind the discovery recommendation

- SCINE Chemoton official overview: https://scine.ethz.ch/download/chemoton
- Chemoton 2.0 research paper, including exploration-completeness challenges:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11516015/

External backend suitability must be rechecked during its evaluation; these
references do not establish STORCA integration or local installation readiness.
