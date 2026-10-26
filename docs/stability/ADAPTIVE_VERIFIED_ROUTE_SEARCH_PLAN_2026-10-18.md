# Adaptive verified reaction-route search: implementation plan

Date: 2026-10-18. Status: **software implementation and zero-native replay
complete; native route success and efficiency not yet validated**. See the
[implementation handoff](ADAPTIVE_VERIFIED_ROUTE_SEARCH_IMPLEMENTATION_2026-10-18.md)
for the actual APIs, tests, prepared continuation and remaining limits. The
design below is retained as the implementation specification, not evidence
that every proposed scientific benchmark has run. This document authorizes no
native calculation, new physical model, production default change or increase
in compute limits.

## 1. Deliverable and boundary

For a declared mapped reaction, electronic surface and environment, persist
toward a fully verified intended route while an admissible, materially useful
action remains within the authorized compute budget. Select actions using
retained evidence, geometry, explicit diagnostics and measured cost. Minimize
failed-inclusive time to a verified route, rather than merely improving XYZ
files or running a fixed sequence successfully.

The output is one of:

- A verified elementary route with its saddle, meaningful reactive mode,
  bidirectional mapped IRC and independently validated actual reached species.
- A verified multistep route only when every required elementary segment and
  intermediate passes its corresponding existing gates.
- An explicit incomplete search, operational blocker, unsupported physical
  treatment or budget stop, retaining partial evidence and the next action.

Finding a route does not establish the lowest barrier, complete conformer or
mechanism coverage, an accurate rate, or molecular stability. Barrierless
chemistry cannot be forced into a saddle requirement; it needs separately
supported physical evidence. No finite strategy set guarantees every route.
An exhausted search is never evidence that a reaction does not exist.

The user-visible objective is `route_required`. Existing stability allocation
may still decide that an unimportant reaction is not worth verifying. Its
thermochemical/observable pruning must not satisfy an explicit request to find
a route. First-route search and observable-conditioned mechanism convergence
remain different tasks.

## 2. Repository-grounded starting point

| Existing code | Preserve and use |
| --- | --- |
| `storca/route_geometry.py::resolve_route_atom_mapping` | Complete nuclear maps, isotopes, state declarations and ambiguity refusal |
| `route_geometry.py::_mapped_bond_changes` | All forming, breaking and bond-order-changing contacts |
| `route_geometry.py::build_endpoint_complex_seeds` | Public seed adapter and existing orientation manifest shape |
| `reaction_path.py::_stationary_fragment_records` | Validated molecular fragments and typed atomic evidence |
| `reaction_path.py::reuse_stationary_fragment_records` | Read-only native fragment revalidation |
| `reaction_path.py::_run_or_resume` | Owned native execution, valid restart, retained-output reuse and existing input repairs |
| `reaction_path.py::_verify_ts_seed` | TS refinement, frequency/mode and bidirectional IRC requirements |
| `reaction_path.py::_run_bond_change_guided_ts_recovery` | Existing bounded one-/two-coordinate recovery |
| `reaction_path.py::_validate_neb_intermediates` and `_run_segmented_intermediate_nebs` | Real intermediate qualification and segmented verification |
| `reaction_path.py::run_irc_branch_continuation` | Existing evidence-bound downhill continuation |
| `reaction_path.py::_refine_reached_irc_stationary_points` | Mandatory validation of actual reached species |
| `surface_resolution.py::SurfaceResolutionPolicy` | Existing scientific thresholds and bounded strategy configuration |
| `rate_transfer_control.py::NativeExecutionControl` | One exclusive native ledger, failed-attempt charging, resource safeguards and cumulative caps |
| `verification_adapters.py`, `solution_path_verification.py`, `rate_transfer_execution.py` | Existing consumers and scientific handoff contracts |

Concrete deficiencies to repair:

1. Fragment placement chooses a world-axis direction rather than a local
   reaction frame; only the first relevant interfragment contact is used.
2. Connectivity and universal short-distance guards do not qualify approach
   geometry or spectator crowding.
3. Reactant/product fragment rotations are independently chosen; unchanged-core
   alignment alone does not produce an economical paired path.
4. The first acceptable orientations are selected, not a ranked diverse set.
5. Seed files can be regenerated over retained candidate paths on resume.
6. A fixed finite-separation NEB endpoint is treated as an asymptote by the
   submerged-climbing-image diagnostic without qualification.
7. One failed guided recovery can end the orientation loop even if alternatives
   remain. Several action-local failures are treated like campaign-level stops.
8. TS/IRC success is selected before reached-minimum refinement; that later
   validation can fail. Preliminary path success is not the finish line.
9. Independent fragment reuse is unnecessarily coupled to full retained-path
   admission. The transfer adapter also hardcodes one orientation.
10. A frozen transfer contract binds live implementation bytes. Editing those
    files and ignoring the resulting hash failures is not a valid continuation.

The earlier four-member-ring audit remains binding: preserve legitimate strained
geometries, graph-based exceptions and stereochemical checks. Do not replace
those guarantees with stricter universal distance or ideal-angle cutoffs.

## 3. Small architecture, explicit responsibilities

Add two focused modules, not a new electronic-structure framework:

### `storca/reaction_seed_geometry.py`

Pure geometry proposals, paired-endpoint preconditioning, metrics and
deduplication. Consume normalized graph/coordinate data from `route_geometry`;
do not import `reaction_path`, execute chemistry, or assign populations/rates.

Suggested records:

- `ReactionGeometryContext`: exact map, graph edits, fragment sources,
  conformer IDs, state, environment, allowed motions and stereochemistry.
- `SeedCandidate`: paired XYZ bindings, generator and parameter identities,
  reaction-center constraints, geometry metrics, diversity aliases and parent
  evidence. Rejection does not reject the reaction channel.

### `storca/route_search.py`

Deterministic policy, evidence reduction, frontier selection and orchestration
through injected existing stage adapters. Do not import the execution module
at module initialization or create a second process runner/accounting system.

Suggested records:

- `RouteSearchPolicy`: objective, supported strategy set, preparation bounds,
  eligibility/ranking rules, validation policy and source inventory.
- `SearchAction`: kind, immutable inputs, parent evidence, scientific context,
  parameters, preconditions, descendant-cost envelope and selection reason.
- `ActionOutcome`: execution, reader, physical-validation and applicability
  results separately; native reservation IDs and all retained outputs.
- `SearchSnapshot`: reconstructed completed evidence, tried fingerprints,
  pending frontier, unresolved claims and authoritative budget snapshot.

Proposed pure API:

```python
derive_geometry_context(route, mapping, fragments, environment)
propose_seed_pairs(context, policy)
assess_seed_pair(candidate, context, policy)
rank_diverse_candidates(candidates, evidence, policy)

reduce_search_state(policy, decisions, evidence, native_budget)
eligible_actions(state)
choose_next_action(state)
```

`reaction_path.py` remains the execution/validation adapter and public
compatibility boundary. A new explicit search-policy branch of
`run_generic_reaction_path` calls the controller. Legacy policies retain their
existing replay semantics. A route-search policy is an internal frozen recipe,
not a new collection of required user flags.

Extend `surface_resolution.py` with an `EndpointReferenceAssessment`, not a
third general framework. Separate:

- `finite_encounter`: atom-ordered path endpoint, not an energy asymptote.
- `bound_minimum`: genuinely validated bound stationary structure.
- `qualified_fragment_reference`: method/state/environment-matched separated
  fragment evidence, including interaction/separation diagnostics when needed.
- `unresolved`: reference cannot support a mechanistic classification.

Reference qualification permits a particular comparison; it does not itself
prove a rate model or barrierless mechanism.

## 4. Search state machine and action ownership

```text
Reopen immutable sources + acquire existing controller ownership
    -> recover/revalidate paid evidence
    -> reconstruct frontier and remaining budget
    -> choose one eligible evidence-grounded action
    -> execute using existing stage adapter and native controller
    -> seal native outcome, then scientific assessment
    -> checkpoint and reconstruct state
    -> verify complete success OR choose the next action
```

The planner owns search decisions. Existing native control owns dispatch and
costs. Existing physical readers own scientific acceptance.

Do not implement the controller by repeatedly calling the entire current
`run_generic_reaction_path`: it already contains nested endpoint, intermediate,
restart, segmentation, TS/frequency and IRC executions. That would hide work
and duplicate recovery.

Expose or factor the necessary stage boundaries incrementally:

1. Prepare/reopen fragment evidence.
2. Propose/qualify an immutable endpoint pair.
3. Prepare/inspect the initial band.
4. Execute or resume one NEB attempt.
5. Assess a retained saddle/intermediate candidate.
6. Run a bounded guided scan and select a candidate.
7. Run the existing TS/frequency/mode/IRC verifier.
8. Continue a qualifying IRC branch or validate an intermediate/segments.
9. Validate actual reached species and finalize the canonical route result.

An initial adapter may wrap an existing composite executor, such as
`_verify_ts_seed`, provided its possible descendant jobs are explicit and every
native child remains under the same active controller. Do not duplicate its
scientific validation. A geometry failure must not cause a helper to run an
unplanned second search before the outer controller can inspect the outcome.

Preserve canonical `route_verification` fields and the existing
`verified_transition_state_path` / `verified_multistep_path` distinctions.
Construct `arkane_inputs` only after the existing final evidence gates pass.
Search progress never grants rate, lifetime or reinjection permission.

## 5. Generalized, compute-aware geometry preparation

### 5.1 Reuse and context

Reopen applicable optimized fragments, supplied structures, prior bands and
candidate TS evidence first. Native inputs, outputs, Hessians, coordinate
order, method/state, temperature convention and owned process receipts remain
bound. A previously useful TS on another context is only a proposal unless its
full exact-context evidence can legitimately be reused.

For a new molecule, an unvalidated geometry may be used only as a labelled seed
proposal. Do not label it a minimum to avoid paying for validation. Keep the
current fragment-validation sequence in the first implementation; postponing
frequencies is a later measured optimization, not removal of a physical gate.

### 5.2 Proposal portfolio

Use reaction-local proposals where justified, backed by general paired
placement:

- H transfer: place the attacker beyond the selected H along the donor-H axis;
  sample distinct bends/azimuths when the reacting fragments require them.
- Addition/substitution: sample applicable site directions and attack faces;
  do not declare one mechanism from a name or graph pattern alone.
- Multiple contacts: jointly satisfy all relevant contacts, not `relevant[0]`.
- Intramolecular changes: bounded torsion/distance-geometry proposals; rigid
  fragment placement alone cannot solve internal conformational requirements.
- Unknown motifs: joint rigid translation/rotation with contact accessibility,
  spectator repulsion and matched-region alignment. Label geometric ambiguity.

Optimize only geometric placement in this stage. Keep rigid fragment internal
coordinates intact, preserve stereochemistry and use proper rotations. Do not
freeze spectators during the later physical path optimization.

Molecular, linear, atomic and symmetric sites need explicit degenerate-frame
handling. Do not force a unique frame using arbitrary laboratory axes and then
claim rotation invariance. Enumerate unresolved directional freedom in a
bounded, reproducible way.

### 5.3 Cheap screening and ranking

Hard reject malformed identities, nonfinite coordinates, impossible overlaps
and violated declared stereochemistry. Use element-/role-aware spectator
checks; intended forming/breaking contacts have their own narrowly scoped
treatment. Ambiguous close contacts are not automatic chemical absence.

Rank lexicographically at first: integrity/applicability; severe crowding;
reaction-center access/alignment; paired spectator motion/path distortion;
diversity; qualified low-level evidence and expected cost when available.
Do not invent calibrated success probabilities or tune scores to experiment.

Select a diverse frontier rather than only the lowest encounter energy. Cheap
model or force-field energy ranks proposals, not channel existence. Deduplication
must retain map aliases, isotope labels and distinct stereofaces/conformers.
Seed counts never become reaction degeneracy or symmetry multipliers.

### 5.4 Initial path and optional electronic preparation

ORCA already performs IDPP preparation. Expose/reuse its initial-path artifacts
before expensive work where the installed executable supports that separation.
Add explicit IDPP/S-IDPP or retained TS-guess input options only after native
input/output compatibility is tested. Do not pretend IDPP is new to STORCA.

Inspect intermediate images for severe collisions, unnecessary spectator
connectivity changes and avoidable rotations. These are preparation diagnostics,
not universal mechanism exclusions. Record ORCA's actual prepared endpoints:
its `Prepare_Frags` option can move nominal inputs. A qualified separation must
not be silently changed by native preprocessing.

Use restrained xTB proposals/path refinement only when geometry alone leaves
useful ambiguity and the declared model is applicable. A failed low-level path
does not veto the channel. Preserve a diverse alternative that does not depend
on the same approximation. No automatic spin-state ranking from ordinary
GFN1/2-xTB. No unconstrained relaxation that silently changes the declared
reactive approach and is then called a successful seed repair.

The current xTB H-bond cluster sampler and neutral-singlet conformer-rate
producer are not general reaction tools. Reuse their isolation/provenance
patterns or small pure helpers without widening their scientific contracts.

First delivery uses geometry, existing NEB, retained seeds and existing guided
scans. Optional xTB, constrained conformer ensembles and template/DG acceleration
must be integrated only with explicit bounded costs and tests; they are not
required to close the immediate H/methane placement defect.

## 6. Adaptive strategy selection and non-repetition

Initial priority is evidence reuse, not a fixed universal method ordering:

1. Complete/reanalyze an already paid qualifying verification stage.
2. Continue a qualifying IRC or reached-minimum recovery.
3. Refine an eligible retained saddle candidate rather than restart a full path.
4. Resume a useful owned band when evidence supports continuation.
5. Try the best diverse qualified endpoint pair or a supported guided scan.
6. Generate additional distinct placements/conformers or use an applicable
   inexpensive preconditioner when the previous diagnosis warrants it.

Eligibility takes precedence over priority: the crowded old H/methane band is
not automatically the first action just because it exists. Record why paid
evidence is reused, used as a new hypothesis, or deemed unsuitable for restart.

| Evidence/failure | Next admissible action |
| --- | --- |
| Wrong approach or spectator crowding | Change implicated placement/separation; use another distinct pair |
| Paired endpoint misregistration | Joint registration or initial-band repair |
| Internal conformer obstruction | Another relevant conformer within supported scope |
| Credible interior saddle seed | Existing TS verifier, with mapped-mode and IRC gates |
| Useful unconverged band | Existing evidence-bound restart |
| Guided scan produces no candidate | Another eligible distinct action; not unconditional campaign termination |
| TS mode mismatch or wrong mapped IRC | Preserve result, reject intended-route claim, choose another hypothesis |
| Incomplete qualifying IRC | Existing bounded branch continuation |
| Reached molecular endpoint is not a minimum | Existing lawful minimum recovery, then revalidate |
| Genuine validated intermediate | Existing segment verification with explicit recursion limits |
| Unqualified endpoint-energy reference | Geometry/reference repair; no automatic dissociation-model claim |
| Native execution/reader failure | Existing operational repair/offline reanalysis; no invented chemical diagnosis |
| Unsupported surface/environment or solver | Explicit escalation/blocker; no blind geometry loop |

Do not infer that a submerged saddle necessarily requires a particular kinetic
theory merely from a fixed encounter comparison. Keep geometric warnings,
qualified energetic evidence and actual method applicability separate.

An action fingerprint includes action kind, ordered/map-aware geometry or band
hashes, state/environment, method, parameters and validator policy. Directory
names and orientation numbers are not scientific distinctness.

- A previously failed physical action cannot be reissued under a new filename.
- New geometry, conformer, band evidence or another justified input change
  must explain the next attempt.
- Execution-only retry references the stopped native receipt and the existing
  resource/runtime recovery policy.
- Successful reanalysis/reuse must dispatch no native job.
- A new physical model/method is not a routine geometry retry.

Keep cheap proposal spaces finite per recipe, but do not install an arbitrary
global "one attempt plus one retry" rule. Under a sufficiently funded policy,
the controller must demonstrate continuing after multiple failed attempts when
an eligible alternative remains. Unknown multi-coordinate strategies remain
unsupported; the existing two-coordinate scan cap is not silently broadened.

## 7. Budget, progress and stopping

Use the existing native journal as the sole charged execution ledger. Planner
decisions reference it; they do not grant free calculations or become a second
budget authority. All nested work uses the same active control. Policy-local
surface budgets remain additional bounds, never a means to reset campaign cost.

Before exploratory work, compute a conservative **minimum remaining category
capacity** for mandatory completion: remaining TS optimization/frequency, IRC,
actual reached molecular Opt/Freq, atomic energy evidence and applicable
continuations. Preserve explicitly allocated recovery capacity where declared.
Forecasts do not reserve native attempts or charge hypothetical work. Recompute
after each accepted/rejected evidence stage. Wall/core-time forecasts are
estimates, not a guarantee of finishing; never impose another arbitrary short
timeout to manufacture apparent efficiency.

For a route-only task, QZ energies/Arkane are not completion requirements. For
the existing rate-transfer campaign, preserve capacity for its separately
specified downstream rate work or explicitly stop at a route-only checkpoint;
do not silently repurpose its budget. Native preparation or optional xTB cost
must be classified and recorded honestly. The current controller refuses
EnGrad/unknown jobs; do not disguise them as single points or run them outside
accounting. Defer such a backend until its narrow accounting extension is ready.

Stop only for:

- Fully verified intended route, including actual reached species.
- User cancellation.
- Hard authorized budget/strategy limits reached.
- Explicit unsupported physical treatment requiring a new method or authority.
- No eligible action in the supported policy, with every exhausted/inapplicable
  alternative and remaining frontier recorded.

Memory/disk shortages and live/ambiguous ownership are resumable operational
blocks, not chemical exhaustion. No heuristic "low expected gain" stop is
introduced before benchmark evidence supports one. A checkpoint should name
the useful next action and its missing prerequisite/budget when available.

Use states such as `searching`, `awaiting_verification`, `verified_route`,
`verified_multistep_route`, `blocked_resources`, `ownership_unresolved`,
`unsupported_physical_model`, `budget_exhausted`,
`search_space_exhausted_under_policy`, and `cancelled`. Per-action outcomes such
as `ts_mode_mismatch` remain truthful but are not automatically campaign-terminal.

## 8. Persistence, environment and compatibility

Suggested small artifact layout:

```text
route-search/
  contract.json                  immutable search context/policy/sources
  candidates/<content-id>/       immutable seed pair + proposal assessment
  actions/<sequence>-<id>/       immutable decision + existing stage artifacts
  checkpoint.json                reconstructible summary, never sole evidence
route-verification.json          existing canonical result + search provenance
```

Use existing atomic artifact/evidence utilities; introduce no general database,
workflow language or provenance platform. A crash between decision, reservation,
native completion and assessment is recoverable from each separate evidence
layer. Native output without an admitted owner is not adopted automatically.

Cache keys contain exact identity/map/isotope/state, environment, input geometry,
method/settings, restraints, generator/validator version and source hashes.
Changing the search policy does not invalidate raw paid chemistry automatically;
changed scientific assessment needs a new receipt under the new policy. Never
rewrite historical assessments or input files to match the current result.

The sequence carries environment throughout. Gas and already-supported implicit
solution contexts can use applicable molecular adapters. Fixed surfaces/cells,
explicit participating solvent and boundary constraints require their own
allowed-motion model; do not rotate away a substrate or relabel a microsolvated
cluster as bulk solution. No new environment support is claimed in this work.

## 9. Current H + methane: explicit source-bound continuation

Verified planning snapshot, 2026-10-18:

- Campaign: `runs/h-ch4-authorized-20261017/campaign`.
- Execution contract SHA256:
  `d81b7a60973ba20648127614796b2714425922c2bc253e5652a1b817b53e2646`.
- Sixteen native events, eight charged jobs, no pending reservations.
- Spending: 4 optimization/path, 3 frequency, 0 IRC, 1 single point.
- Owned execution: 3329.2745259169023 wall seconds;
  13317.09810366761 allocated core seconds.
- Prefix tip:
  `000016-f85dddf929ec994304ab50957715b88813279a9ed87d59dc367c50892b9d9f98.json`.
- Seven jobs support H, CH4, H2 and CH3 fragment evidence. The eighth is the
  interrupted NEB, still charged; no verified saddle/route/rate follows.
- Original remaining category allowance: 8 optimization/path, 7 frequency,
  2 IRC, 8 single point. Recompute from the live ledger before any execution.

### First implementation action, before editing bound source

Capture the exact bytes of every implementation source bound by the original
execution contract in a sealed source capsule. Each copy must match its
original recorded SHA; retain original path-to-copy mapping. Copy the complete
required local dependency closure for replay as appropriate and distinguish
extra archived dependencies from the original bound inventory. Do not edit
first and then invent a snapshot from changed source or unrelated Git content.

Add a narrow continuation adapter, provisionally
`storca/rate_transfer_continuation.py`, with prepare/freeze/reopen/run entry
points in the existing transfer driver. Preparation is zero-dispatch and
unauthorized. Original `reopen_execution` remains strict.

The continuation seals:

- Original plan/contract bindings and the exact historical source capsule.
- New implementation inventory and search policy.
- Original native control contract/root and immutable journal prefix.
- Original runtime, physical settings, mapped reaction and counting assumptions.
- Read-only revalidated fragment evidence.
- A new geometry/search workspace outside attempted old basenames.
- Separate explicit continuation authorization.

Historical source evidence is reopened through exact-byte archived bindings,
not an "ignore changed hashes" option. Raw calculation identity and completed
native sources must still match. Reinterpret the old endpoint warning only in
a new assessment; preserve its original report as history.

Do not call `run_execution` / `_run_execution` unchanged from the continuation:
they reopen the original live-source-bound inventory, including at closeout.
Factor the shared stage body behind explicitly validated execution context, or
call the stages through the continuation's own strict lineage validator. Never
catch and suppress an integrity failure. Preserve historical
`_orca_job_input_evidence` input-plus-referenced-file digest semantics: native
controller activation rechecks old reservations before allowing cached reuse.
Intentional new live code is separately bound; historical capsule bytes and
original inputs/referenced XYZ files must remain exact.

Proposed driver commands, to be implemented in
`tests/run_rate_transfer_benchmark.py` (not available yet):

```text
archive-sources ORIGINAL_CAMPAIGN SOURCE_CAPSULE
prepare-continuation ORIGINAL_CAMPAIGN SOURCE_CAPSULE DESTINATION
freeze-continuation PLAN WORKDIR
inspect-continuation WORKDIR
run-continuation WORKDIR --rmg-env ENV
resume-continuation WORKDIR --rmg-env ENV
```

The archive subcommand can be added to this currently unbound driver before
editing any bound implementation module. Archive/preparation/inspection never
dispatch native work. Method and search details live in the frozen plan, not
new required user flags. Run/resume use the same existing native controller.

Run under `NativeExecutionControl(original_campaign / 'native-control',
original_control_contract)`. New chemistry appends attempt nine to the same
ledger, with all old charges, failures and original caps intact. Do not create
a zero-count child campaign. Never rerun the old `candidate-01/path/neb-ts`
basename with new coordinates.

Decouple fragment-only reuse from successful full-path reuse. Make retained
endpoint correspondence optional where only fragment evidence is being
imported, while requiring exact route/map/state/method/settings and all owned
native bindings. Missing old evidence refuses; it is not silently regenerated.

Four cores, MaxCore 1000 MB/rank, six hours/native job, memory/disk safeguards
and original cumulative caps remain. Continue using `--rmg-env rmg4_env`.
Do not reopen the closed OH/methane angular campaign or launch OH + H2 here.

Implementation permission is not native continuation authorization. The
implementation run ends with zero-dispatch readiness evidence and the exact
proposed next command/authorization request unless execution is separately
and explicitly authorized.

## 10. Ordered implementation work packages

Each package must leave an independently testable checkpoint; implement the
core sequence rather than stopping after the H-abstraction special case.

1. **Preserve original evidence and qualify endpoint references.** Archive old
   bound sources first; add explicit reference types; prevent unqualified
   encounters from granting a dissociation/capture classification. Preserve
   valid bounded early-stop diagnostics and all final physical requirements.
2. **Implement immutable general seed preparation.** Pure geometry records,
   all-contact paired placement, H-transfer local proposals, generic fallback,
   role-aware screening, diversity and no-overwrite resume. Preserve legacy
   manifests through an adapter.
3. **Expose existing execution/validation stages.** Factor only necessary
   orchestration boundaries. Keep existing readers/executors; make nested
   recovery ownership and possible job costs explicit. Add independent
   fragment-only admission.
4. **Implement deterministic adaptive planning.** Action eligibility/ranking,
   meaningful-difference fingerprints, persistent frontier, failure-specific
   recovery, completion-capacity checks and honest terminal states.
5. **Integrate full-route finalization.** Do not stop at preliminary TS/IRC
   success or one failed guided scan. Require actual reached species. Wire
   opt-in policy through the transfer bridge and existing generic API without
   changing unrelated production defaults or downstream rate contracts.
6. **Implement the narrow historical continuation adapter.** Reopen old source
   capsule/new sources, same native ledger and exact fragments; prepare a new
   unauthorized continuation workspace. No migrated zero-count budget.
7. **Harden by replay and adversarial tests.** Run the required existing and new
   suites; replay the real H/methane evidence with native dispatch forbidden.
   Document validated software behavior separately from scientific performance.
8. **Prepare the scientific evaluation and handoff.** Freeze a separate bounded
   geometry/search benchmark proposal, report missing capabilities and provide
   one exact continuation command once authorized. No automatic benchmark or
   transfer launch.

Optional xTB/conformer/template acceleration is subsequent measured work. Do
not install packages, introduce ML/LLM decisions, switch electronic methods,
add state/environment lanes, or implement new spin-crossing/photochemistry/
pressure-dependent theories during this increment.

## 11. Required tests and acceptance gates

### Geometry and reference tests

- Saved H/methane crowding is rejected before native execution; corrected
  proposals face the selected H from outside its donor bond.
- Proper global rotations/translations and correctly mapped atom/fragment
  reordering preserve approach-quality metrics and equivalent frontier choices.
- Atomic, linear, symmetric, isotopic, radical/ionic, strained-ring,
  intramolecular, stereochemical and multi-contact fixtures remain well defined.
- All relevant contacts participate; internal rigid fragment distances and
  declared stereochemistry remain intact during rigid preparation.
- Distinct conformers, mapped channels and stereofaces are not lost through
  deduplication. Aliases do not change symmetry/degeneracy conventions.
- Malformed proposals never dispatch. A cheap-model failure cannot veto the
  channel. Soft idealized geometry preferences do not reject legitimate rings.
- Elevating only an unqualified encounter energy cannot change the physical
  classification into a proven barrierless/dissociation mechanism.
- Actual native endpoint changes are detected; nominal XYZ qualification
  cannot survive unexamined native fragment repositioning.

### Controller and verification tests

- A bad arrangement leads to a materially different next action.
- A failed guided scan does not end a funded search with another eligible action.
- A funded scenario succeeds after more than two unsuccessful attempts.
- Same state/evidence yields the same next decision. Renaming a workspace does
  not enable an identical failed physical action.
- Retained near-saddle/band/IRC evidence is used only when admitted; inappropriate
  paid work does not get priority merely because it was expensive.
- Native execution failure, reader failure, physical rejection and unsupported
  model produce distinct decisions; live/ambiguous leases prevent duplicates.
- TS/mode/IRC evidence without validated reached minima is not final success.
  Existing allowed endpoint recovery can continue from that partial result.
- A multistep success claim fails if any required segment, intermediate or
  reached-species evidence is missing/invalid, even if another segment passes;
  it must go through the same complete-evidence finalizer.
- Wrong mode or wrong mapped IRC never verifies the intended channel.
- Early loss-pruning cannot satisfy `route_required`.
- All children of composite actions are charged and fit explicit stage
  envelopes. Verification-capacity checks do not create reserved/free attempts.
- Interruption after every decision/reservation/native/assessment boundary
  resumes without duplicate work. Budget exhaustion preserves partial results.
- Supported gas and supported implicit-solution contracts remain distinct;
  unsupported environments do not fall back to gas.

### Continuation tests and actual zero-job replay

- The exact eight-job/16-event fixture reopens without changing the prefix.
  A new successful or failed native attempt is number nine and remains charged.
- Original/current source mutation, missing archive, changed runtime/settings,
  prefix gap/mutation, unowned output and pending work refuse as appropriate.
  This means original artifacts/archived historical bytes or unexpected current
  source changes, not the intentional new implementation explicitly bound by
  the continuation. Historical reservation digests remain identical.
- Reopen seven-job fragment evidence with the native runner patched to raise;
  assert all counts, native hashes and original artifact bytes are unchanged.
- Newly proposed seed files live elsewhere; the original report and crowded
  XYZ remain unchanged. The new assessment names its own policy/source binding.
- Preparation/inspection of an unauthorized continuation dispatches nothing.

Keep existing mapping, geometry/ring/stereo, native input/restart, TS frequency,
mode, IRC completion, reached-minimum, thermochemistry, control/runtime and
composite-consumer tests. Deliberately revise tests whose premise was that
every fixed encounter energy is an asymptote or one failed guided scan must
end all search; do not weaken physical tolerances to obtain passing tests.

Proposed new test files:

```text
tests/test_reaction_seed_geometry.py
tests/test_route_search.py
tests/test_route_search_integration.py
tests/test_rate_transfer_continuation.py
tests/replay_h_ch4_route_search.py       # zero-dispatch real saved evidence
```

Extend the existing `test_generic_reaction_path`, `test_transfer_mapping_atoms`,
`test_rate_transfer_execution`, `test_rate_transfer_control`, relevant geometry
audits and road-map traceability tests rather than reproducing their fixtures
and scientific validators in a separate fake workflow.

## 12. Scientific efficiency benchmark: separate authorization

Compare three arms to separate geometry from adaptive-planning benefit:

A. Existing placement and fixed search sequence.
B. Improved geometry with the fixed sequence.
C. Improved geometry with adaptive search.

Lock methods, final validators, fragment inputs, permitted information and
failed-inclusive compute budgets before predictions. Held-out reference TSs
cannot also be seed proposals in the same test. Sharing a previously paid
fragment is allowed only when equally available to all relevant arms.

Initial proposed panel: H + CH4, OH + H2, radical addition of H to ethene, and
cyclobutene ring opening. Use retained OH + methane only for offline evidence
compatibility, not another angular campaign. Add substitution, multi-contact
cycloaddition and conformer-sensitive cases before claiming broad transfer.

Measure:

- Fraction of fully verified intended routes at fixed compute.
- Time/core-hours/gradient calls to first fully verified route.
- Failed/censored costs for every assigned reaction, not only successful cases.
- Native attempts, reuse/resume, preparation overhead and distinct alternatives.
- False success claims (required zero) and pathways lost through filtering.

Audit a fixed diverse subset of cheap-screening rejections at a higher level
when separately budgeted, to detect false negatives. One verified route is not
proof of the lowest-barrier route. A four-case pilot cannot establish universal
robustness or a calibrated probability of success.

Do not compute QZ energies or Arkane rates just to score geometry/search
efficiency. Freeze a separate panel budget derived from the actual required
native stage counts; this document does not allocate extra native attempts.

## 13. Definition of implementation complete

The implementation run is complete when:

1. The adaptive sequence operates through existing physical validators, persists
   across interruption and continues past recoverable failures without repeats.
2. No seed, warning or preliminary TS/IRC result can grant final route success.
3. Failed-inclusive budgeting, shared-ledger continuation and exact source
   lineage pass adversarial tests.
4. Real saved H/methane fragment reuse/preparation is demonstrated with **zero
   native dispatch**, original evidence unchanged and all eight charges intact.
5. A separately authorized next native action can be prepared reproducibly,
   with an exact command, remaining budget and unresolved assumptions.
6. Roadmap and closeout distinguish implemented software from native route
   success and measured compute savings, neither of which tests alone prove.

No claim that H + methane now has a verified route is permitted without the
actual future native validation.

## Primary methodological references

- Robertson and Habershon, paired position/orientation preconditioning,
  DOI 10.1002/jcc.26495:
  <https://wrap.warwick.ac.uk/id/eprint/150277/7/WRAP-Simple-position-orientation-preconditioning-scheme-minimum-energy-path-calculations-2021.pdf>.
- autodE reaction-local translation/rotation and bond rearrangement:
  <https://duartegroup.github.io/autodE/reference/transition_states/locate_tss.html>.
- ORCA NEB, IDPP/S-IDPP, native fragment preparation, optional low-level paths
  and TS guesses:
  <https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/neb.html>.
  Test the installed 6.1.0 behavior explicitly; a living manual can describe
  later patch-level behavior.
- xTB open-shell/spin-model limitations:
  <https://xtb-docs.readthedocs.io/en/latest/spgfn.html>.

These sources motivate proposals; published geometry statistics are not
STORCA verification rates, and their algorithms do not replace STORCA's final
physical evidence gates.
