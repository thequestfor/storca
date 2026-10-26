# Stability Traceability Matrix

This matrix connects scientific requirements to implementation, tests,
benchmarks, and roadmap milestones. “Implemented” means code exists; it does
not imply experimental validation. GitHub issue numbers should be added when
the milestone issues are created.

| Requirements | Milestone | Current implementation | Tests / benchmark profile | Status |
|---|---|---|---|---|
| STAB-GOV-001–002 | M0 | Roadmap, requirement catalog, lane manifest | `test_stability_roadmap_traceability.py` | implemented_unvalidated |
| STAB-CAM-001–004, STAB-CAM-006, STAB-ENV-001, STAB-PROV-003–005 | M1 | `chemical_system.py`, `stability_campaign.py`, preparation planner/CLI, immutable resume, examples | campaign contract/planning tests | implemented_unvalidated |
| STAB-CAM-005 | M1 | independent gas and restricted implicit-solution campaign adapters, exact resume, conditional claims, no cross-environment aggregation | `test_gas_campaign_adapter.py`, `test_solution_campaign_adapter.py`, `test_stability_campaign_execution.py` | implemented_unvalidated |
| STAB-DISC-001–005, STAB-DISC-008 | M2 | native candidates, RMG and condition-specific graph/xTB adapters, bounded radical-association/H-abstraction proposals with spin compatibility, schema-v3 solvent-bound ALPB screening and narrow neutral carbon four-pi hypotheses, environment-bound immutable assembly, source-preserving merge | `test_reaction_discovery.py`, `test_discovery_assembly.py`, `test_stability_verification_wiring.py`, `test_graph_xtb_solvent_screen.py`; three-case native screening diagnostic, not a rate/coverage benchmark | implemented_unvalidated |
| STAB-DISC-006, STAB-UNC-001 | M2; STAB-008 | symmetric disagreement evidence, empty-search shared-blind-spot gate, independent-scope gaps, uncalibrated endpoint-screen exclusion gate, separated uncertainty register, immutable coverage revisions, observable claim withholding | `test_discovery_coverage.py`, `test_stability_verification_wiring.py`; `coverage_risk_manifest.json` | implemented_unvalidated |
| STAB-DISC-007 | M2 | selective physics-driven discovery escalation pending | reaction-discovery blind-spot benchmark pending | planned |
| STAB-PATH-001–003, STAB-PATH-006 | M2/M5 | Existing endpoint, NEB-TS, frequency, IRC, barrierless, and solution solvent audit; exact downhill execution ownership, retry-policy/job/source hashes and conflicting-source refusal | gas workflow, `test_irc_completion.py`, read-only 119-frame/102-file native replay, and `test_solution_path_verification.py`; actual retry equivalence remains open | implemented_unvalidated |
| STAB-PATH-007–009 | M2/M5/M6; STAB-011 | `surface_resolution.py`; nonconverged stationary recovery; mapped scans; per-image connectivity/fragment-motion diagnosis; connectivity-first intermediate identity/native candidates; independent segment TS/frequency/IRC; common attempt-evidence envelopes; environment applicability gate; multistep Arkane bundle and repaired-network handoff; route-wide lease budget | `test_surface_resolution.py`, `test_intermediate_identity.py`, `test_generic_reaction_path.py`, `test_verification_adapters.py`, `test_verification_engine.py`; `surface_resolution_profiles_v1` | implemented_unvalidated |
| STAB-PATH-004, STAB-RATE-003 | M3 | bounded RDKit reactant seeds → existing native ORCA minimum refinement → matched free energies and rate admission; unresolved interconversion/TS coverage withholds aggregate rates; immutable cache sources and shared budget | `test_conformer_rate_production.py`; installed-RDKit acquisition/gas cache replay; real conformer-rate benchmark and full ensembles pending | in_progress |
| STAB-RATE-004 | M3; STAB-010 | explicit single-temperature ORCA frequency production, immutable route/endpoint/session cache contract and condition-matched solution Eyring conversion | `test_solution_thermochemistry_temperature.py`; real solution temperature/reference benchmark pending | implemented_unvalidated |
| STAB-RATE-001–002, STAB-RATE-008 | M3/M4 | Arkane gas rates, restricted solution Eyring evidence with condition- and observed RRHO/quasi-RRHO recipe matching, and hashed dimensional first/second-order mass-action rate contracts; in-progress selected-saddle forward/reverse coefficient production without network admission | gas regression; linear and nonlinear analytical regression; `test_solution_thermochemistry_recipe.py`, synthetic `test_solution_rate_coefficients.py` and unchanged retained-parent recipe replay; real coefficient benchmark pending | implemented_unvalidated |
| STAB-PROP-001 | M1/M6 | checksum-validated repair libraries admitted as native RMG core seeds, identity-bound restart resolution, precision-preserving Chebyshev export, nonfinite-rate rejection, and source-bound resonance-aware correspondence shared by rate verification and planning; stereo ambiguity fails closed | `test_rmg_repair_seeding.py`, `test_rmg_export_precision.py`, `test_rmg_identity_bridge.py`, `test_reaction_correspondence.py`, `test_replacement_identity_handoff.py`, `test_rmg_bridge_numerics.py`, `test_verification_engine.py`, `test_flux_verification.py`; retained NNONC identity replay refuses unsupported stereo admission; full propagation remains open | implemented_unvalidated |
| STAB-PROP-002 | M1/M4 | exact linear concentration/extent propagation | `test_solution_stability.py` | implemented_unvalidated |
| STAB-PROP-003–004, STAB-ENV-004 | M4 | separate verified-input nonlinear mass-action lane; dimensional first/second-order rates, reversible target reformation, diffusion ceilings, audited stiff-method escalation | analytical regression implemented; experimental solution benchmark pending | in_progress |
| STAB-OBS-001–002 | M1/M4 | standardized t95, retention, mass-action concentration, product-yield, dominant-product, gross-consumption branching, and integrated-pathway-extent claims; half-life, selectivity, induction time, and broader observable sensitivity remain open | uncertainty, solution session, and mass-action analytical/campaign tests | in_progress |
| STAB-UNC-002–003, STAB-ALLOC-001–002 | M2/M3/M6 | rate intervals, analytical sensitivity, score breakdown, budgets, and simultaneous propagation of verified calculated-rate plus unresolved-route uncertainty; no default Arkane error factor | `test_kinetic_uncertainty.py`, `test_flux_verification.py`, `test_verification_engine.py`, `test_solution_verification_session.py` | implemented_unvalidated |
| STAB-CONV-001–002 | M1/M2 | repaired gas convergence and restricted solution convergence/terminal states | verification and solution session tests | implemented_unvalidated |
| STAB-METHOD-001–002, STAB-PATH-005 | M5 | individual fail-closed checks; common selector pending | method challenge set pending | planned |
| STAB-RATE-005–007 | M6; STAB-012; M6-BENCH-001 | Arkane preparation, explicit hash-backed well/bath-gas collisional inputs or visibly unvalidated estimates, barrierless rejection, source-neutral specialist request/result contract, immutable request export, checksummed result import, retained-path continuation, capture/VTST evidence gates, development rate benchmark, calculated-rate uncertainty propagation, and repaired-network handoff; built-in calibrated specialist solvers pending | `test_arkane_verified.py`, `test_specialist_rates.py`, `test_specialist_rate_benchmark.py`, `test_flux_verification.py`, `test_verification_engine.py`, `test_verification_resume.py`; locked falloff/capture holdout pending | in_progress |
| STAB-ENV-005 | M7 | declared-only microsolvation lane | benchmark pending | planned |
| STAB-ENV-006 | M8 | supplied gas photolysis evidence and TD-DFT screen; explicit irreversible photon-source libraries; autonomous fate lane pending | `test_generated_photolysis_library.py` checks native loading and supplied-rate preservation; scientific photolysis benchmark pending | in_progress |
| STAB-ENV-007, STAB-PROP-005 | M9 | no supported interface lane | transport/reaction benchmark pending | planned |
| STAB-ENV-008–009 | M10 | declared-only periodic surface and electrochemical lanes | specialist benchmarks pending | deferred |
| STAB-PROV-001–005 | M1/M2 | retained contracts, evidence-content hashes, condition-scoped claims, exact adapter cache keys, gas checkpoints, append-only solution sessions, ORCA process leases, stale-route continuation, validated NEB `.allxyz` restart evidence | gas and solution restart tests; `test_verification_resume.py`; `test_generic_reaction_path.py` | implemented_unvalidated |
| STAB-BENCH-001–005 | M11 | small discovery fixtures, software regressions and source-audited JPL development limiting-rate references with corrected f(300) factors | cross-environment benchmark, independent falloff curves, locked holdouts and STORCA calibration pending | planned |

## First implementation issue sequence

1. **STAB-001 — Establish roadmap and machine-checked traceability.** Implemented. Addresses
   M0 and requirements STAB-GOV-001–002.
2. **STAB-002 — Define `ChemicalSystem` and stability-campaign request schema.** Implemented.
   Addresses STAB-CAM-001–002 and STAB-ENV-001.
3. **STAB-003 — Add preparation-first campaign planner and immutable resume.** Implemented.
   Addresses STAB-CAM-003–004 and STAB-PROV-003–005.
4. **STAB-004 — Adapt the existing gas ladder into the campaign runner.** Implemented.
   Addresses STAB-CAM-005 and preserves STAB-PROP-001.
5. **STAB-005 — Adapt solution verification sessions into the campaign runner.** Implemented.
   Addresses STAB-CAM-005 and preserves STAB-PROP-002.
6. **STAB-006 — Ship bounded campaign examples and preparation reports.** Implemented with
   preparation-only, direct solution, live gas, and mixed gas/solution contracts.
   Addresses STAB-CAM-006.
7. **STAB-007 — Automatically assemble RMG and graph/xTB candidates.** Implemented.
   Addresses STAB-DISC-002–006.
8. **STAB-008 — Add discovery disagreement and coverage-risk accounting.**
   Implemented and software-tested; curated blind-spot calibration remains open.
   Addresses STAB-DISC-006 and STAB-UNC-001.
9. **STAB-009 — Add CREST-backed conformer contracts and cache identity.**
   Addresses STAB-PATH-004 and STAB-RATE-003.
10. **STAB-010 — Generalize condition-matched thermochemistry beyond 298.15 K.**
    Explicit-temperature producer slice implemented and software-tested;
    real condition-matched thermochemistry/rate benchmark and solvent-domain
    validation remain open. Addresses STAB-RATE-004.
11. **STAB-011 — Bounded reaction-surface resolution and segmented verification.**
    Implemented and software-tested; retained scientific profiles remain open.
    Addresses STAB-PATH-007–009, STAB-ALLOC-002, STAB-RATE-001–002, and
    STAB-PROP-001.
12. **STAB-012 — Retained specialist-rate export/import and rate-only continuation.**
    Implemented and software-tested. Exact condition/path requests, provider-result
    matching, raw-artifact hashes, library checksums, and append-only network repair
    are enforced. Built-in calibrated VTST/capture production remains an M6 research
    dependency. Addresses STAB-RATE-005–007 and STAB-PROP-001.
13. **M6-BENCH-001 — Specialist-rate conformance and development benchmark.**
    Implemented and software-tested. Capture/VTST physical evidence fields,
    rate/uncertainty scoring, method-applicability challenges, and a development
    reference manifest are enforced. A locked holdout, full falloff curves, and
    live VaReCoF/MESS results remain required before scientific validation.
    Addresses STAB-RATE-005–007 and STAB-BENCH-003–005.
14. **M6-MESS-001 — Explicit MESS preparation and native T/P parser.**
    Implemented and software-tested. The adapter renders only hash-backed RRHO,
    collision, energy-transfer, and external Rotd evidence; bundle paths do not
    contaminate scientific cache identity. It remains non-executing and cannot
    claim a rate. Addresses STAB-RATE-005–006 and STAB-PROV-003–004.
15. **M6-MESS-002 — Bounded, checkpointed MESS execution.**
    Implemented and software-tested. Prepared hashes are revalidated; the
    executable binary, stdout/stderr, native outputs, and parsed exact-channel
    grid are retained under a hard timeout. Parsed rates remain unadmitted and
    cannot support a lifetime. Addresses STAB-RATE-005 and STAB-PROV-003–005.
16. **M6-MESS-003 — Evidence admission and repaired-network handoff.**
    Implemented and software-tested. Native bimolecular amount basis is never
    inferred; executable-bound unit evidence, complete uncertainty components,
    collision ceilings, independent equilibrium ratios, reverse-channel
    detailed balance, and exact conditions are required. The fitted library is
    reloaded by RMG and checked over every computed T/P point before a standard
    specialist result can enter network repair and repropagation. Live falloff
    validation and a locked holdout remain open. Addresses STAB-RATE-005–007,
    STAB-PROP-001, and STAB-PROV-003–005.
17. **M6-BENCH-002 — Falloff-curve and holdout-integrity benchmark contract.**
    Infrastructure implemented and software-tested; scientific curation remains
    in progress. Exact finite-pressure grids, curve RMSE and pressure-span
    metrics, content-hashed holdouts, prediction-freeze evidence, and explicit
    reference-data-pending cases are enforced. The first H + O2/N2 case is
    preregistered but has no admitted numerical curve, so it cannot support a
    lane claim. Addresses STAB-RATE-005–007 and STAB-BENCH-003–005.

Each issue must use the scientific acceptance and support-change procedure in
`DEFINITION_OF_SUPPORTED.md`. The **Prediction-first execution order —
2026-10-01** section of `COMPREHENSIVE_IMPLEMENTATION_PLAN.md` is the current
dependency order. STAB-011/P1 real path/rate/repair and controlled resume remain
open; STAB-009's initial integrated reactant producer remains scientifically
unvalidated, not permission to postpone a
working prediction or scientific validation. Curated STAB-008/M11 benchmarks,
rate-error evidence and common method applicability advance alongside the
shared-loop repair. Existing analytical sensitivity and supplied-uncertainty
allocation must be reused, not reimplemented as a new engine. Autonomous
ordinary solution is the next environment implementation once those shared
boundaries work; it need not wait for all gas chemistry. Neat liquids, bulk
solids and specialist exposure lanes require complete applications and keep
their existing unsupported states. No requirement or support level is closed
by this planning revision.

2026-10-01 expanded-goal checkpoint: first ten prediction-first work items are
authorized, with P1 first rather than P1-only scope. Item 2 has a bounded,
reference-isolated evaluator, not a frozen numerical panel. Item 4 has a
hash-bound ordinary Arkane method decision, not universal applicability. Item 7
now admits source-valid earlier-only channels through the existing bounded
queue without historical kinetics; held-out consequential recovery remains open.
Items 1/3 have exact terminal NEB checkpoint admission and 143 passing focused
regressions, plus a 516-file unchanged real admission replay; no saddle, rate,
repaired numerical prediction or controlled resume equivalence is proved.
Item 5 source review supplies calibration prerequisites, not calibrated bounds.
No environment support level or scientific requirement is closed by these checks.

2026-10-03 reconciliation: actual P1 dispatch 007 is terminal with no completed
NEB iteration, own band, rate or observable. The original production-driver
ledger binds all seven attempts and leaves 1078.8381125830347 seconds; no budget
was renewed. Items 1/3 now have exact-execution trajectory/saddle selection and
completed-iteration interruption admission, with failing-before regressions
and an unchanged 45-source no-job timeout replay. These are evidence-integrity
repairs, not full-chain scientific equivalence. Item 2 has eight newly audited
CD3CN rate coefficients but no admitted full-case numerical cohort or holdout.
Item 5 retains a 341-record/5,250-rate source audit with zero admitted calibration
rows. Item 8 has a bounded isolated arm64 dependency/platform audit; runtime,
export, restart and consequential benefit remain unproved, so the backend is
deferred rather than rejected. Retained test logs and detailed limitations are
in `WORKING_PIPELINE_GATE.md`; subsequent source changes require a new full suite.
All first-ten scientific closure gates and later environment milestones remain
explicitly open unless separately demonstrated.

2026-10-03 producer/validation checkpoint: item 6 has explicit single-T frequency
production and condition-matched solution conversion, but not integrated
conformer rates. Items 4/6 now also bind actual ORCA 6.1 RRHO/quasi-RRHO entropy
recipes, reference/cutoff, scaling, pressure, version and source hashes; missing,
contradictory or mismatched recipes withhold rate admission. The installed
parent replay preserves 27 files, launches no chemistry and observes 100 cm^-1
quasi-RRHO reference treatment. Gas Arkane reconstruction is unchanged.
Thirteen recipe methods and 119 affected solution/panel tests pass; final full
application validation passes 948 tests (22 expected skips), and native focused
validation passes 68. The preceding checkpoint passed 935 application tests.

Items 1/3 additionally have exact-execution/hash-bound NEB intermediate cache
selection, seven new methods, 116 combined path/resume tests and an unchanged
83-file actual source replay. Correct seed dispatch is not a new optimized
intermediate or observable. Item 2's descriptive eight-rate overlay passes
32 tests and reports eight missing/zero scored predictions; full-case readiness
and holdouts are unchanged. No first-ten scientific closure gate or environment
support level is promoted. Detailed artifacts and current resource state are
in `WORKING_PIPELINE_GATE.md`.

2026-10-03 integrated-reactant-conformer checkpoint: item 6 now dispatches
bounded identity-preserving seeds through the existing ORCA endpoint executor,
matched temperature/recipe parsing and rate admission. Tests intercept the
electronic launcher, not the whole conformer producer; distinct minima and
known sampling/refinement gaps withhold aggregate rates. Nineteen new conformer
methods and five session-policy methods support 144 affected passing tests.
Exact input/reference/source ownership and complete mode evidence repair
demonstrated cache-admission failures; expired shared budgets do not prevent
completed-source reuse. Old session-produced TST rates cannot gain current
production-policy credit on resume, while exact declared imports remain imports.

Full application validation passes 972 tests (22 expected skips), with 68 native
focused tests. Actual installed-RDKit acquisitions for three systems repeat
selected geometry hashes; two are explicitly refinement-selection-truncated.
Four retained gas opt/freq contracts pass the new cache digest check and 804
native files remain byte/timestamp unchanged, with zero quantum jobs. This
does not establish real solution conformer free energies/rates, interconversion,
TS/encounter/product ensembles, calibrated unsampled risk or an experimental
prediction. No item 1–10 scientific closure gate or environment support level
is promoted. See the conformer checkpoint evidence and working-pipeline ledger.

2026-10-05 IRC/reference reconciliation: items 1/3 now have exact-execution
downhill continuation assembly, native job/reference/source hash checks,
append-only retry ownership and integer core-count policy consistency.
Conflicting repeated source observations cannot overwrite earlier proof.
The 118-test affected suite passes; the reproducible actual NNONC replay
reconstructs 119 byte-identical frames and preserves 102 native files/timestamps
and directories with no chemistry. The actual source has no execution retry;
retry/fallback/restoration/corruption coverage is intercepted workflow evidence.
It does not prove chemical endpoint admission, new rates or full-chain resume.

Items 2/5 and M6 additionally have a retained JPL source/visual audit and
predecessor-byte-hash regression for eight limiting-rate anchors from four
evaluated association recommendations. The f(298) → f(300) reference amendment
does not alter central rates, conditions, partitions or acceptance thresholds.
No exact confidence probability, independent pressure samples, prospective
holdout or STORCA rate-error calibration is supplied. The scientific panel's
full-case readiness and all first-ten closure gates remain unchanged.

At the IRC/reference checkpoint, application-007 passed 976 tests in 107.129 seconds
(22 expected skips); native-006 passed 68 in 9.725 seconds. These counts and
earlier checkpoints remain historical.
Detailed commands, source/log hashes and remaining scope are recorded in the
working-pipeline ledger and both 2026-10-05 evidence documents. No support level
is promoted; later matter/exposure milestones and existing budgets are intact.

2026-10-05 in-progress solution checkpoint (items 9/10; STAB-RATE-001–002/008):
the existing native path producer now accepts a separate candidate coefficient
request without concentrations, duration, species table or invented error
bound. A solvent-bound native adapter and existing CLI reach it without new
flags. Both directions share a verified saddle, matched thermochemistry and
explicit 1 M molecularity/unit conversion using separated two-molecule channel
sums and selected one-molecule endpoints. Shared-thermochemistry detailed
balance and source ownership are retained. Coefficients are conditional estimates
with unresolved error/ensembles/diffusion, not admitted network kinetics or a
stability bundle. Synthetic tests do not establish real paths, physical accuracy,
autonomous discovery or supported solution fate. The completed initial-producer
checkpoint passes 999 application tests in 157.101 seconds (22 expected skips)
and 68 native tests in 16.729 seconds; see
[producer evidence](../../runs/solution-rate-coefficients-diagnostic-001/EVIDENCE_2026-10-05.md).
This software result does not close physical validation gates. Automatic applicable uncertainty, coverage,
transport and actual unimolecular/bimolecular benchmarks remain open.

2026-10-06 item 9/shared M2 discovery checkpoint (STAB-DISC-001–005/008):
requested SMD/CPCM candidate conditions and observed GFN2-ALPB/`gsolv` screen
evidence are separate. Bounded mapped neutral carbon four-pi opening/closure
proposals preserve inventory; no gas fallback or rate inference is allowed.
Native convergence/connectivity, method/version/invocation, stdout/JSON and
hash checks plus environment-bound cache contracts fail closed. Schema-v3
artifacts do not relabel historical gas evidence; mismatches refuse without
overwriting sources. The shared path gate's reproduced `phase="solution"`
mismatch is repaired under unchanged mandatory solvent, SMD/CPCM and exact
keyword checks, not expanded physical applicability. The retained diagnostic used
six actual one-core xTB endpoint jobs (0.924418 seconds batch elapsed), recovered
three expected products and preserved 24 files per case on zero-job cache checks.
No ORCA path, coefficient, lifetime, coverage or calibration was demonstrated.
Final software logs pass 1018 application tests in 151.461 seconds (22 expected
skips), 115 expanded native/discovery-boundary tests in 9.695 seconds and 179
affected tests in 1.234 seconds. Nineteen new regressions include the demonstrated
cache-identity and phase-handoff defects. Earlier 999/68, 1016/113 and 1017/114
checkpoints remain historical, not new scientific acceptance evidence.
The full autonomous solution loop and all first-ten/later-environment gates
remain open; see
[screening evidence](../../runs/solution-discovery-screen-diagnostic-001/EVIDENCE_2026-10-06.md).

2026-10-06 source-admission/reference follow-up (STAB-DISC-001–005/008,
STAB-RATE-001–002/008, STAB-BENCH-001–005): selected solution graph candidates
reopen their native contracts/results, hashes, endpoint states/identities,
observed charge/UHF settings, solvation/invocation and energy sums before native
preparation, then recheck sources afterward. Three actual intact candidates
prepare; three deliberately corrupted derived copies refuse, with all 77
original files and 13 directories unchanged and no chemistry dispatched.
`tests/test_solution_discovery_sources.py` exercises corruption, selection,
source races and immutable request binding; synthetic fixtures are not chemistry
validation. `tests/replay_solution_discovery_sources.py` retains actual sources.

The first audited descriptive gas comparator is frozen in
`benchmarks/stability/reference_sources/2026-10-06-jpl2011-hpl/source-audit.json`;
`scientific_panel.py` binds rates/conditions/scopes to that audit and the batch
driver refuses contradictory method scopes. Actual offline report
`runs/scientific-panel-jpl-hpl-002/report.json` contains one missing prediction,
one audited limiting-rate reference and ten pending entries. The theory-informed
JPL recommendation supplies neither finite-pressure experimental validation nor
STORCA error calibration. No first-ten or later-environment closure is claimed.

2026-10-06 item 7/shared discovery-to-verification queue checkpoint:
`stability.py` and `verification_engine.py` retain independently screened direct
gas hypotheses even without quantitative RMG flux, using the existing bounded
queue. Exact condition IDs, initial-reactant/atom-balance checks, graph-only
native source and all explicit supported gas aliases govern this new lane;
unsupported/conflicting environment declarations do not gain fallback credit.
Six fixture regressions plus the existing affected checks pass **137 tests in
0.904 seconds**. Fixtures are not physical path/rate validation.

The [actual-source 002 replay](../../runs/independent-discovery-flux-fallback-diagnostic-002/result.json)
selects two retained CNONN gas candidates through ordinary engine planning and
actual native prepare-only mapping/geometry. RMG failure is explicitly simulated
and the path callback is injected prepare-only, not the default executor.
Both statuses are `prepared`, not computed/verified paths; all 145 source files
and 17 directories plus the condition source remain hash/mtime unchanged.
There are zero native jobs, rate-builder calls or propagation calls. The
summary is `verification_incomplete`, with empty verified routes and null loss,
coverage/model coverage and lifetime. No rate/model/observable credit follows
from the retained structural handoff. Scientific P1/full-chain, uncertainty,
coverage and first-ten/later-environment acceptance gates and budgets remain
unchanged; details and replay hash are in `WORKING_PIPELINE_GATE.md`.

2026-10-07 items 1/3/shared verification: native008 is terminal with completed
OptTS/frequency and an incomplete timed-out IRC. The new source-bound stage
reuse reopens exact conditions/map/method/context, owned job dependencies,
stationary geometry/Hessian/spectrum/mode and absent terminal processes. It
prepares only a new bidirectional IRC, preserves original job references, and
re-admits sources before kinetics. Synthetic interception proves no OptTS/Freq
regeneration; actual no-job preparation preserves 940+198 source files. Genuine
job-contract-owner and declared-context-pointer/hash negatives are retained,
with actual metadata-substitution replay before/after the narrow context fix.
This does not verify a connection, rate, full-chain equivalence or lifetime.
The original eight-attempt ledger leaves 362.6558744990034 monotonic seconds;
no new allowance or scientific gate is granted. See the reuse evidence ledger.

Item 7's exact existing pure-target gas alias is added to the strict fallback
allowlist; conflict/solvent/surface declarations remain rejected. Actual replay
003 preserves 145 files/17 directories and prepares two candidates with simulated
RMG failure, zero jobs and null scientific claims. The alias itself is tested
with fixtures; the actual source remains its unchanged gas-surrogate model.
Item 2 additionally retains an OH + H2O2 primary-source access/condition audit
lead without changing panel readiness or treating an incomplete visual source
as a reference. All ten and later-environment closure requirements remain open.

2026-10-07 items 2/4/5 evidence update: the separate Berg2024 OH + methane slice
contains a primary measured descriptive comparator with source/hash/T/P/unit
binding, while its actual evaluation stays missing (zero scored). Fourteen
measurements remain one shared-systematics study, not independent calibration
units. Shared panel/holdouts/support are unchanged. A retained P1-derived
malformed-envelope replay demonstrates preserved valid cost without chemical
metrics. Six new source regressions plus accounting/panel tests pass 54/1.363 s;
full application/native validation has not yet been rerun for this update.

The installed-native optical audit records the actual C2 TS candidate's inferred
count two versus forced one, with identical-mode partition-only ratios and
unchanged endpoint ratio. It does not supply an IRC connection, rate multiplier
authority or calibrated error. Configurational/channel inventory and explicit
degeneracy are still required before a recipe correction. Original sources and
compute limits are preserved. See the dated source/accounting/counting evidence
linked from `WORKING_PIPELINE_GATE.md`; no closure or environment promotion.

Item 5 allocation safety now rejects unbound/non-string method-calibration
condition IDs at construction/load, preserves specialist request-ID type, and
does not stringify malformed current planner IDs. Actual public-consumer
synthetic replays retain the original and laundering failures; exact matches
survive and changed/unknown current bindings are excluded. Root's 126-test
affected replay passes in 1.075 s with zero chemistry. Invalid historical records
are preserved but cannot newly authorize allocation; unbound user/database
scenarios stay explicitly supplied, not calibrated. Source/recipe/domain
authentication and automatic applicable uncertainty remain unimplemented; see
the admission evidence and working-pipeline ledger.
