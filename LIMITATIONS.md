# Limitations

- Results depend on the selected electronic-structure method, basis set, and
  input structure; the defaults are not appropriate for every molecule.
- A frequency calculation with no large imaginary modes is evidence for a
  local minimum, not a guarantee of global stability or experimental isolation.
- IR conformer populations currently use optimized electronic energies rather
  than quasi-harmonic conformational Gibbs free energies. Duplicate minima that
  differ only through symmetry-equivalent atom mappings are not yet guaranteed
  to be collapsed.
- The current dimer/trimer environment branch uses bounded, restrained
  finite-temperature trajectories, not a bulk liquid ensemble. Its sufficiency gate prevents an
  undersampled frequency spread from being displayed as a calculated width,
  but passing the minimum gate does not supply unbiased bulk populations, extended
  hydrogen-bond networks, anharmonic lifetimes, or benchmark validation. The
  restrained GFN2-xTB seeds deliberately cover interaction coordinates. The
  reported weights are decorrelated equal-time frame occupancies conditional
  on that finite restrained multi-seed protocol. They are not equilibrium
  liquid occupancies and include no restraint-bias reweighting. Hydrogen mass
  repartitioning is used for integration stability, so these trajectories are
  configurational samplers rather than physical dynamical models. Their ORCA Hessians
  are evaluated at restrained-xTB snapshots rather than unconstrained DFT
  stationary points. Adaptive batch convergence tests numerical stability of
  this bounded representative set; it does not establish bulk thermodynamic
  convergence or turn its conditional occupancies into bulk thermodynamics.
  Strongly nonstationary snapshot Hessians are now marked diagnostic-only. An
  environment-preserving DFT refinement and unrestrained gradient gate are
  implemented. Projected local-mode finite differences are implemented and can
  override matching transfer pairs, but automatic execution for every failed
  representative is not enabled until the additional ORCA invocations fit the
  global run budget and broader molecular validation passes.
- Additive DFT correction transfer is currently limited to local X-H stretches.
  It is exposed only for mode/coordination classes that pass withheld-
  representative error and sampled-domain coverage gates. xTB intensities are
  diagnostic and are not transferred as DFT intensities. The gated hybrid
  diagnostic instead uses the representative-ORCA mean intensity within each
  local mode class, not a learned snapshot-specific intensity correction. A passing harmonic
  transfer corrects xTB toward the chosen ORCA harmonic method, not directly
  toward experiment; it does not add anharmonicity or bulk populations.
- A numerically stable restrained ensemble can still have broad
  bootstrap uncertainty. Ensemble extension is capped and reported as
  statistically imprecise when that occurs; it does not manufacture precision
  or convert conditional restrained-trajectory weights into liquid occupancies.
- General internal-coordinate fingerprints now track stretches, bends,
  torsions, and ring deformations across equal-sized Hessians. Generalized
  cross-size dimer/trimer matching is not complete: mixed cluster sizes still
  use the validated target-local X-H path, and coupled fingerprint-region modes
  require additional atom-mapping and degenerate-subspace validation.
- `spectrum_intrinsic.csv` is separated from the rendered measurement, but the
  current ATR response remains a first-order penetration-depth approximation.
  Exact ATR optics, transmission absorbance, and absolute intensity require
  condition fields such as refractive indices, concentration, and path length;
  recording those fields does not yet mean every transfer model uses them.
- ORCA, xTB, RMG, RDKit, and Open Babel are external or optional dependencies.
  `storca doctor` only checks whether they are discoverable; it does not
  validate an installation or a scientific setup.
- `typerint.py` is archived prototype code. Its broader integrations are not
  supported by STORCA 2.0, and its Multiwfn path is explicitly out of scope.
- Every current RMG stability stage is a homogeneous gas-phase surrogate. The
  low-pressure nitrogen stage is not a literal vacuum or isolated-molecule
  trajectory. Dry/humid gas models do not represent ordinary liquids, solids,
  solutions, droplets, crystal packing, container walls, surfaces, metal or
  impurity catalysis, diffusion, or heat and mass transfer.
- RMG explores chemistry available to its databases, families, species
  representation, tolerances, and finite resource limits. A completed retained
  model with no target-loss flux is not proof that an omitted mechanism does
  not exist. Solver coverage and full repaired propagation prevent false use of
  an incomplete model; they cannot make the reaction search chemically
  exhaustive.
- The current gas propagation contract is isothermal and constant-pressure;
  it does not predict self-heating, thermal runaway, or heat-transfer-limited
  chemistry. Earlier independent Cantera outputs used its energy-on default
  and may have self-heated despite a fixed-temperature request. Their fluxes
  and lifetimes must be repropagated, not relabeled as isothermal. Missing or
  inconsistent reactor-condition evidence now withholds admission. Corrected
  propagation alone does not validate a mechanism, a rate, or its uncertainty.
- The lane registry is not a claim that every registered environment is
  supported. Current gas adapters and two narrowly restricted verified-input
  implicit-solution lanes are `observable_supported`; their evidence gates still
  apply. Microsolvation, periodic-surface, and electrochemical entries are
  `declared_only`; admission fails closed and no lifetime can be retained from
  them.
- `storca stability-campaign` has live adapters for admitted homogeneous-gas
  `t95` lanes and the restricted neutral implicit-solution `t95` lane.
  Preparation launches no external programs. Explicit gas execution reuses the
  existing RMG/ORCA/Arkane/repaired-propagation gates; it does not make RMG
  complete or DFT exact. Solution execution either propagates a retained
  verified first-order network, invokes the condition-bound single-temperature verification
  session, or propagates a separately evidenced ideal-activity first/second-order
  mass-action network. Microsolvation, surface, electrochemical, and other environment
  requests remain blocked or adapter-not-implemented rather than being mapped
  onto gas or implicit-solution physics. Only named enforceable budgets are
  accepted; the solution core-hour ceiling is a declared estimate, not measured
  scheduler accounting, and no universal ORCA-job ledger exists. Evidence
  content hashes are immutable identities, but they do not establish the
  scientific truth of an external claim.
- The implicit-solution lane does not discover solution reactions, run ORCA, or
  derive solution rates from gas calculations. It propagates externally
  verified first-order rates after validating their condition, standard-state,
  path/rate provenance, balance, and restricted network-coverage claims. A
  formally valid external claim can still be scientifically wrong; provenance
  is auditable evidence, not independent validation.
- Solution t95 bounds combine explicit multiplicative bounds on all
  target-consuming first-order rates. They are not calibrated posterior
  intervals and do not include mechanism uncertainty outside the supplied
  restricted coverage claim. A horizon crossing that changes across the bounds
  is reported unresolved.
- The first solution propagator is a constant-volume linear concentration model.
  It excludes bimolecular encounters, diffusion, solvent molecules as reactants,
  activity coefficients, ionic strength, acid/base equilibria, ion pairing,
  phase separation, transport, heat transfer, photochemistry, surfaces, and
  target-reforming networks. Any declared occurrence blocks t95 rather than
  being approximated by a gas or pseudo-first-order rate.
- The separate mass-action solution lane supports verified first- and
  second-order reversible networks and target reformation, but it does not make
  solution mechanism discovery autonomous. It assumes ideal molar activities,
  fixed supplied species and concentrations, constant volume, neutral species,
  and a complete-network claim within that restricted scope. It requires a
  sourced diffusion ceiling for each bimolecular rate but does not calculate
  encounter complexes, viscosity-dependent diffusion, ionic strength, pH,
  activity coefficients, dynamic acid/base speciation, solvent molecules as
  participants, phase separation, or transport. Those omissions can dominate a
  real solution lifetime and therefore remain applicability failures.
- Mass-action uncertainty is a deterministic rate-bound sensitivity set. Every
  bound corner is propagated only for networks with at most ten reactions;
  larger networks use joint and one-at-a-time extremes. The latter is not a
  rigorous interval enclosure, neither mode is a probability distribution, and
  neither covers missing mechanisms. The initial-Jacobian stiffness diagnostic
  can miss stiffness that appears only later in a trajectory. Experimental
  calibration and broader stiff-solver benchmarks are still required.
- CPCM and SMD names in an externally supplied solution bundle identify the
  claimed environment; they remain claims rather than independently validated
  calculations. For one ordinary neutral unimolecular candidate,
  `solution-path-verify` can generate an auditable solvent-consistent ORCA
  endpoint/NEB-TS/frequency/IRC path and conventional Eyring rate with an observed
  condition-matched RRHO/quasi-RRHO free-energy recipe. This producer is
  restricted to one explicitly requested thermochemistry temperature per job,
  one selected TS path, transmission coefficient one, and an
  externally justified composite barrier bound. Named CPCM/SMD parameters stay
  fixed; explicit temperature matching does not prove their temperature-domain,
  isotopic-solvent or phase applicability. The non-298 K producer has software
  tests, not a new real solution accuracy demonstration. It
  does not provide conformer completeness, tunneling, anharmonicity,
  explicit-solvent dynamics, capture/VTST, or proof of mechanism completeness.
  The recipe check is restricted to ORCA 6.1 output evidence, rejects missing
  or inconsistent settings, and does not calibrate a rate. Quasi-RRHO entropy
  is not explicit hindered-rotor/anharmonic treatment. Its ideal-gas
  translational/rotational assumptions are not removed by CPCM/SMD. Existing
  externally supplied-result bundles and historical claims are not silently
  rewritten or retroactively validated; gas Arkane rates retain their separate
  statistical-mechanics recipe.
- The in-progress coefficient-only producer accepts an internally selected
  neutral singlet candidate with one or two separated molecules on each side,
  without concentration history, duration, supplied kinetics or an invented
  uncertainty bound. It reuses the native path producer and reports conditional
  forward/reverse selected-saddle coefficients, with explicit 1 M conversion,
  dimensional units and shared-thermochemistry detailed balance. Summed separated
  species are not replaced by an encounter-complex free energy. These packets
  have unresolved calibration/ensembles/interconversion/path degeneracy,
  tunneling/recrossing, diffusion and activities, and explicitly lack network
  admission authority. They do not bypass the stricter uncertainty/coverage
  requirements for a stability bundle or mass-action prediction. The native
  candidate adapter requires retained solvent-bound discovery evidence; it is
  not an autonomous solution-discovery engine. Synthetic tests and pre-dispatch
  refusal of an actual unsupported retained route are not a positive solution
  path/rate benchmark, a validated solvent domain or a new supported lane.
  The coefficient producer's dated 999-test application/68-test native checkpoint does
  not close those physical gates; see
  [producer evidence](runs/solution-rate-coefficients-diagnostic-001/EVIDENCE_2026-10-05.md).
- The integrated solution reactant-conformer producer is limited to one small
  neutral closed-shell singlet, ordinary isotope treatment and at most 20 heavy
  atoms. It uses bounded RDKit/UFF geometric seeds and selected ORCA refinement,
  not exhaustive or solvent-equilibrated sampling. Rotor/ring triggers can miss
  OH/NH rotation, inversion and tautomerism; the symmetry-match cap can make even
  an otherwise small molecule unresolved. UFF energies never become physical
  populations or rates. Its seed-wall check is cooperative between bounded
  native calls, not hard process preemption. Failed/incomplete/truncated samples
  prevent admission; a single refined cluster still needs an independently
  justified sampling/model bound. Distinct minima require interconversion and
  TS-conformer evidence; a conditional rapid-equilibrium single-path scenario
  is not admitted as an aggregate network rate. Explicit conformer-state
  kinetics, product/encounter/TS ensembles, rotor/anharmonic treatment and
  calibrated unsampled risk remain missing. New refinement admission requires
  complete 3N frequency evidence and immutable exact native job/source contracts;
  this does not retroactively impose those stronger checks on all old endpoints.
  Tests intercept electronic calculations; retained installed-RDKit acquisitions
  and gas cache checks are not a real solution-rate benchmark.
- `solution-verify-session` allocates only the current monotonic first-order
  t95 observable under its explicit single-temperature contract. Its analytical sensitivity correctly assigns zero
  t95 influence to downstream reactions, but it does not yet optimize product
  yields or branching. Approximate candidate rates are scheduling evidence,
  not posterior distributions or admissible solution kinetics. The session
  requires an independent restricted coverage claim and cannot prove that an
  undiscovered pathway is absent. If downstream candidate kinetics remain, a
  t95-only claim may converge but the general concentration/product bundle is
  deliberately withheld.
  Resume now withholds convergence credit for old session-produced TST rates
  with missing/changed temperature/recipe/conformer production policy, preserving
  old artifacts. Explicitly supplied rates remain imports; neither the policy
  check nor imported provenance independently validates their physical accuracy
  or replays all calculation sources.
- The opt-in graph/xTB lane is an independent but deliberately small graph-edit
  search. It covers target-unimolecular edits and, for explicitly declared gas
  partners, pairwise radical association and H abstraction, plus narrow
  valence-preserving carbon four-pi ring-opening/closure hypotheses. It does not
  compute
  encounter/capture rates and omits closed-shell bimolecular substitution,
  solvent participation or catalysis, termolecular chemistry, catalytic/surface
  chemistry, excited states, arbitrary multi-bond edits, and dynamics. A search status of
  `completed_declared_search`, `completed_with_unresolved_candidates`, or
  `discovery_truncated` only describes the declared bounded edit space; none
  establishes reaction-network completeness. Requested backend failure,
  unresolved RMG endpoint conversion, failed xTB candidate screening, or
  truncation now vetoes an otherwise converged lifetime.
- The new solvent-bound graph/xTB screen is not autonomous solution kinetics.
  Requested SMD/CPCM conditions and observed GFN2-ALPB/`gsolv` screening are
  separate evidence; ALPB endpoint totals include solvation terms but are not
  thermal barriers, SMD/CPCM free energies or rate constants. Only parameterized
  named solvents are accepted, with no gas fallback. Three actual development
  cases used six one-core xTB jobs and recovered the expected connectivity
  hypotheses in water/acetonitrile; this is not an experimental mechanism,
  orbital/stereochemical-path, rate, uncertainty or coverage benchmark. The
  isolated-carbon four-pi edit omits specified reacting stereochemistry and
  broader pericyclic chemistry. Native convergence/connectivity and source
  identity checks do not establish xTB accuracy. Graph schema-v3 and
  environment-bound assembly contracts refuse mismatches and mutated sources
  without rewriting old artifacts;
  historical gas screens never become solution evidence. See
  [screening evidence](runs/solution-discovery-screen-diagnostic-001/EVIDENCE_2026-10-06.md).
  The final 1018-test application/115-test expanded native-boundary checkpoint
  checks software, not scientific accuracy. The shared path gate's
  `phase="solution"` alias repair retains mandatory named-solvent/SMD/CPCM and
  exact keyword checks; it is not a solvent-domain or physical-support grant.
- Observable-conditioned coverage is a deterministic fail-closed policy, not a
  calibrated probability that the mechanism is complete. Independent-source
  agreement and candidate-set stability under one bounded expansion can reduce
  identified risk but cannot rule out shared blind spots. The current policy is
  explicitly marked `not_blind_spot_calibrated` until provenance-backed,
  preregistered held-out mechanism cases pass the coverage benchmark.
- RMG-only stability runs now intentionally withhold a lifetime or no-loss
  conclusion under that policy. Even RMG plus graph/xTB cannot support an empty-
  network conclusion without applicable blind-spot calibration or complementary
  evidence. This increases abstention and has not yet been scientifically
  calibrated against the preregistered chemistry strata.
- Every retained graph-only survivor must either receive a verified rate and
  observable-robust repropagation or remain a coverage blocker. Increasing
  `--graph-xtb-max-candidates` changes high-level verification cost, not the
  meaning of unqueued candidates; those candidates are not silently discarded.
- In a condition ladder, graph/xTB discovery is condition-specific so reactive
  O2/water hypotheses are not transferred into inert rungs. This does not make
  mechanisms, barriers, rates, or product yields transferable between gas
  rungs. These gas ladder searches are not liquid-solvent kinetics, surface,
  electrochemical or photochemical discovery. The separate ALPB screen above
  does not add those physical models. A failed, state-unresolved, or truncated
  requested search is retained as an explicit ladder evidence gap.
- GFN2-xTB supplies retained endpoint electronic-energy screening in that lane,
  not reaction prediction, transition-state barriers, rates, conformer
  coverage, or kinetic relevance. A graph-only candidate is evidence to
  investigate, not evidence that it controls a lifetime. Failed xTB screens
  are retained rather than silently removed. Candidates excluded above the
  configured endpoint-energy window also veto an observable conclusion when no
  independent backend retains the same native reaction. An RMG-retained copy
  remains in the candidate network and the xTB result becomes disagreement
  evidence rather than an omitted-pathway blocker. There is currently no
  runtime bypass. Any future exception must require a hash-bound,
  observable-, condition-, environment-, screen-, and chemical-scope benchmark
  calibration; reaction energy is not treated as a barrier or absence proof.
- Automatic ORCA/Arkane verification is default for modelled target-loss
  routes, but it can finish only when RMG retains resolved endpoint structures,
  charges, multiplicities, reaction identity, and quantitative flux. Ambiguous
  atom mappings, multiple unresolved spin surfaces, large mapping searches, or
  missing upstream routes fail closed rather than being guessed.
- Source-bound RMG family atom tracing can resolve maps lost by generic
  graph-edit heuristics, but is only discovery evidence. Library reactions with
  no retained family recipe, unsupported stereochemistry/electronic states,
  inequivalent recipe maps, or hydrogen traces incompatible with the existing
  explicit-H layer cannot receive recipe-mapping credit. This handoff has been
  exercised on retained cyclobutene ring-opening evidence in RMG 4; it is not
  validation of all reaction families, methods, environments or mechanisms.
  Changed source or recipe hashes invalidate reuse, not justify relabeling an
  old calculation. It does not reset the engine's same-policy retry limits.
- Generic NEB-TS searches sample a finite set of generated endpoint complexes
  on one declared adiabatic electronic surface. Separated species are validated
  as individual minima, but their finite-separation encounter geometries are
  path seeds rather than stationary complexes. The search can miss a
  lower-energy conformer, encounter orientation, intermediate, or transition
  state. Spin crossings, electronically nonadiabatic paths, and strongly
  multireference bond breaking are not resolved by this workflow.
- Surface-resolution diagnostics reduce repeated equivalent NEBs but cannot
  prove that an unseen pathway or intermediate is absent. Per-image distance-
  based connectivity changes are diagnostics, not quantum bond orders or proof
  of an elementary mechanism. RDKit connectivity is determined before bond-order
  assignment and is admitted only for a single charge/spin-compatible graph;
  disconnected ion/radical pairs and ambiguous electronic states remain
  unresolved. A graph identical to a declared endpoint is retained as an
  endpoint return and may dismiss only that specific warning image; this does
  not prove that other intermediates are absent.
- One/two-coordinate relaxed scans can rescue poor seeds but can miss collective
  coordinates, solvent participation, conformational gating, or bifurcating
  surfaces. A scan maximum never supports a rate without subsequent OptTS,
  frequency, mode, and IRC evidence.
- Persistent endpoint-highest and submerged-climbing-image diagnostics are
  bounded method-allocation evidence, not proof that a hidden saddle or alternate
  pathway does not exist. They prohibit an inapplicable harmonic-TST rate and
  request dissociation/capture or pressure-dependent treatment; those specialist
  rate models are not yet implemented by the default Arkane adapter.
- Neutral gas-phase 1-to-2 dissociation can now receive a conservative
  condition-specific loss upper bound from validated endpoint thermochemistry,
  detailed balance, and an inflated reverse collision ceiling. This is only a
  screening bound. If it crosses the requested loss threshold, the route
  remains unresolved and it is not relabeled as a VTST, capture, falloff, or
  RRKM rate.
- A validated specialist-rate provider interface and retained-result
  export/import workflow are available, but STORCA does not yet ship a default
  VTST, capture, MESMER/MultiWell, or barrierless Arkane producer. Imported
  evidence must match the exact request, declare applicability and uncertainty,
  hash raw artifacts, and checksum the generated RMG library. These checks do
  not independently establish that an external solver or its calibration is
  scientifically correct. Without such evidence, the default run still returns
  `dissociation_rate_model_required` or `barrierless_rate_model_required`.
- The specialist-rate benchmark contains only development limiting-rate
  anchors and a method-applicability challenge. It has no locked holdout or
  digitized falloff curves, so it tests plumbing and exposes numerical
  disagreement but cannot validate a capture or master-equation backend. A
  provider that passes the current points still does not establish mechanism
  correctness or broad applicability.
  The JPL 2011 audit confirms four pairs of evaluated limiting coefficients and
  corrects their reference envelopes from f(298) to f(300). It does not turn
  those eight anchors into independent experiments, specify an exact confidence
  probability, or supply STORCA rate-error calibration. The air-bath
  recommendations cannot be relabeled N2-only; the OH+OH association coefficient
  is not total OH loss. Evaluated falloff-model points would not become
  independently measured finite-pressure data or a prospective holdout.
- The MESS adapter now has a complete software admission path from a retained
  native table to a checksummed RMG pressure-dependent replacement library.
  It does not supply the missing scientific evidence automatically. In
  particular, the native `cm^3/sec` amount basis must be established from a
  source applicable to the exact executable, equilibrium ratios and collision
  ceilings must be supplied independently, and uncertainty must cover all
  named model components. The requested condition must currently be an exact
  grid point; STORCA does not silently interpolate during admission. No live
  MESS falloff case or locked holdout has yet passed the benchmark, so this is
  an evidence-enforcing integration, not a supported general master-equation
  prediction lane.
- Arkane central rates now require explicit uncertainty propagation before a
  nonzero-flux route can support a terminal observable. STORCA deliberately has
  no universal default multiplier: the necessary calibration depends on the
  reaction class, electronic-structure method, conformer treatment, tunneling,
  and (for falloff) collisional energy-transfer model. A run without applicable
  supplied calibration therefore ends with
  `calculated_rate_uncertainty_missing_for_<route>` rather than false precision.
- Unless route-specific collisional evidence is supplied, the built-in Arkane
  1↔2 lane uses RMG's generic `SingleExponentialDown` energy-transfer estimate
  and estimated transport. The exact bath gas, master-equation algorithm, and
  grain settings are retained, but this is not a calibrated collisional model.
  An explicit supported contract can provide sourced well/bath-gas parameters;
  it currently supports only single-exponential-down energy transfer and
  Lennard-Jones transport. Neither branch supplies a universal uncertainty
  calibration. A central `k(T,P)` may be used for reranking but cannot by itself
  close pressure-dependent uncertainty.
- High-pressure-limit Arkane TST can execute for a pure-target gas without
  inventing an inert bath species. This removes an irrelevant input dependency;
  it does not prove pressure independence or supply self-collision parameters
  for a pure-target master-equation calculation. The latter remains outside
  the current non-target bath-gas adapter.
- The generic surface resolver currently supports homogeneous gas and the
  dedicated single-continuum implicit-solution contract only. Implicit solution
  requires an exact declared CPCM/SMD solvent keyword. Explicit solvent
  participation, liquid microstructure, interfaces, periodic surfaces, solids,
  electrochemical potential, illumination, and transport need their own lanes;
  they return `environment_model_escalation_required` rather than being
  approximated by an undeclared gas or continuum calculation.
- A uniform surface-attempt evidence envelope makes decisions auditable but
  does not improve the physical model by itself. Its environment, mapping, and
  policy links expose assumptions; they do not establish completeness or
  accuracy.
- The endpoint connectivity gate remains geometric, not an electronic bond
  analysis. Its bounded closed-shell four-ring diagonal resolution requires an
  exact independent coordinate-graph match; it is not a general exemption for
  close contacts. Open-shell/unsupported assignments and ambiguous intermolecular
  contacts remain unresolved. Graph agreement does not validate electronic
  structure, a stationary point, or a reaction path. A separate check now
  corroborates declared tetrahedral/E/Z stereochemistry from mapped coordinates,
  conditional on the declared graph and isotope labels. Unspecified or
  unsupported stereochemistry is not independently established. The ordinary
  distance gate still does not assign all bond orders or validate charge, spin,
  multireference character or conformer coverage. The bounded geometry audit is
  a software challenge panel, not an experimental chemical-identity or kinetic
  validation benchmark.
- Verified multistep paths are not collapsed to a single barrier. Each ordinary
  segment needs an Arkane-admissible rate and successful insertion into a full
  repaired network. Channel-specific, barrierless, or rejected network segments
  leave the observable null with a precise pending/specialist state.
- A reproducible downhill ORCA path is not a verified rate. Barrierless
  classification requires two complete endpoint-matched orientations, but the
  default Arkane TST adapter still rejects it until a defensible
  collision/capture and recrossing model is supplied. Such a route therefore
  returns incomplete verification rather than a `t95`.
- Arkane rates inherit ideal-gas, stationary-point, harmonic/rigid-rotor,
  tunnelling, energy-transfer, bath-gas, and electronic-structure assumptions.
  Reinjection proves that the calculated replacement controlled the retained
  RMG model under the exact contract; it does not validate those approximations
  against experiment.
- The gas ORCA-to-Arkane reader corrects demonstrably diatomic external
  rotation and requires a source-bound native mode audit before rate
  reinjection. It does not generally repair near-linear polyatomic modes,
  anharmonicity, internal rotors, electronic-state statistical weights,
  path degeneracy or barrier accuracy. Ambiguous geometry/mode classifications
  fail closed. Historical unaudited Arkane libraries cannot be silently reused
  as verified replacements; a new rate-only workspace is required. A passing
  rotor audit or collision check is not experimental validation or a lifetime
  prediction. See the native rotational-mode admission policy in METHODS.md.
- A verified rate library can be loaded into RMG yet fail to survive exported
  network serialization or species renaming. The launcher now preserves gas
  Chebyshev bounds/coefficients at binary64 text precision, and bridge inspection
  rejects nonfinite forward rates. Old lossy exports are not silently repaired:
  restoration requires the exact checksummed source model and independently
  established reaction correspondence in a separate derived artifact. A
  serialization roundtrip is not kinetic validation. Resonance-aware,
  source-bound replacement matching now checks unique graph correspondence,
  state, direction, stoichiometry, exact provenance and file hashes. It does not
  establish stereochemical identity: potential stereochemistry not preserved
  by RMG artifacts remains unresolved, including the retained NNONC-shaped
  regression. Matching a library name or molecular formula alone is insufficient.
  Historical label-only compatibility is explicitly identified; it does not
  substitute for required structural evidence in new repaired-network runs.
  Fresh/resumed full-network scientific acceptance remains an open gate.
- Generated photon-source libraries now explicitly declare irreversibility for
  native RMG loading. The loading/rate roundtrip validates serialization of a
  supplied rate, not photochemical mechanism, quantum yield, environmental
  applicability, or lifetime. Malformed retained libraries are not silently
  rewritten; corrected outputs require new provenance and checksums.
- Flux-guided verification is bounded by configured iterations and per-job
  timeouts. Surface resolution additionally enforces route-wide calculation,
  elapsed core-hour, wall-time, orientation, intermediate, scan, and recursion
  limits. Process-lease accounting is exact for new streamed ORCA jobs; reused
  legacy calculations are evidence but do not reconstruct historical compute
  cost. Exhaustion returns `surface_resolution_budget_exhausted` and no lifetime.
- Route continuation can resume a validated ORCA NEB `.allxyz` checkpoint, but
  it is not a universal scheduler. New streamed ORCA jobs retain process leases;
  older jobs without a lease rely on a conservative artifact-quiet interval and
  cannot prove process death from files alone. Arkane and repaired-RMG
  interruptions do not yet have equivalent program-level checkpoint recovery.
  Exact terminal NEB time limits can now enter a source-bound checkpoint
  continuation without rewriting failed evidence or renewing the aggregate
  budget. The retained cyclobutene checkpoint-admission replay is not a resumed
  saddle/rate result. Its first SIGINT control did not stop the driver, so real
  controlled interruption/resume equivalence remains unproved.
  The subsequent corrected restart (dispatch 007) timed out before a completed
  NEB iteration and retained no own native band. It is not an additional usable
  checkpoint or evidence that the saddle does not exist. Exact artifact lookup
  no longer displays an older execution's trajectory beside that failure; the
  original incorrect association is retained in immutable historical results.
  Intermediate image selection and validation caches are likewise bound to
  the actual NEB execution. A retained replay changes the dispatched image,
  but deliberately runs no optimizer and proves no changed physical disposition.
  Downhill IRC assembly now uses exact retained execution metadata and source
  hashes, validates retry-chain ownership and core-count policy, and refuses
  conflicting repeated source observations. Retry, serial fallback, restoration
  and corruption handling are software-tested with execution intercepted. The
  real 119-frame/102-file read-only replay uses the ordinary basename, so it does
  not establish physically executed retry equivalence, controlled interruption
  through rate/repair, or a new chemical endpoint/rate/lifetime claim.
- The method-applicability policy before ordinary Arkane execution is an
  auditable declaration/diagnostic gate, not an automatic physical-validity
  oracle. Missing diagnostics remain conditional assumptions; absence of a
  warning does not establish a single-reference surface or calibrated accuracy.
  Source-dataset curation is not a rate-error producer without independent,
  method-matched paired predictions and prospective evaluation.
- Source-consistent solution graph/xTB screening and native preparation are not
  validated paths/rates or an autonomous solution-fate prediction. The selected
  candidate adapter now refuses missing, changed or mismatched native sources
  before path preparation, including observed charge/unpaired-electron settings;
  matching those settings does not establish spin purity or ground-state order.
  Its actual intact/corrupted-copy replay launches no chemistry and proves no
  additional mechanism coverage, lifetime or calibrated uncertainty.
- The scientific panel now contains one audited descriptive high-pressure-limit
  comparator and ten pending entries, with the optional P1 comparator still
  provisional. The new JPL reference is evaluated/theory-informed, not an
  independent measured finite-pressure point. Its prediction is missing. No full
  cohort freeze, prospective holdout, physical support or STORCA calibration
  follows from a correct source transcription or offline scoring infrastructure.
- The scientific-panel evaluator and earlier-only structural candidate queue
  have software tests and saved-source replays. The numerical panel is not yet
  frozen, interval calibration is absent, and consequential held-out missing-path
  recovery has not been demonstrated. Simulated zero-flux queue admission does
  not establish mechanism coverage or a useful stability prediction.
  The eight fitted Rulisek2005 coefficients are separate rate-only reference
  evidence in CD3CN, not executable concentration histories, independent
  Gaussian error bars or validated solution support. Their joint-fit limits do
  not supply STORCA's electronic-structure error model. The RDB7 source audit
  likewise admits no calibration rows; its ring sampling domain excludes the
  nonaromatic four-membered P1 system. A local Chemoton metadata/platform audit
  supports deferral, not runtime admission or scientific backend rejection.
- Recognized ORCA/OpenMPI startup denials receive append-only exact-input
  execution retries and are not interpreted as failed scans or missing saddles.
  This recognition is marker-based and currently covers the retained local
  socket-denial signatures; unrelated cluster, license, filesystem, or launcher
  failures may still require a new explicit execution-failure classifier.
  A stopped serial fallback can restore the immutable requested core count in
  a hash-linked execution segment, but this does not make MPI available: a new
  launcher failure exhausts the bounded retry rather than falling back forever.
- Observable-conditioned allocation uses a local first-order product of retained
  logarithmic `t95` sensitivity and an explicitly supplied multiplicative rate
  bound. It is not a posterior distribution, propagated lifetime interval, or
  expected-information-gain calculation. Sensitivities can change after every
  repaired-network rerun, and externally supplied uncertainty bounds can be
  miscalibrated or outside their stated applicability. Missing coverage disables
  uncertainty ranking rather than making the missing route appear certain.
- Collision-limit rate envelopes evaluate the retained nonlinear network at
  the joint zero-rate and collision-ceiling endpoints while preserving each
  reversible reaction's forward/reverse relationship. They expose whether the
  requested conclusion changes at those physical endpoints, but they are not
  exact rate verification and do not prove reaction-network completeness.
  Multiple interacting uncertain rates or non-monotone interior behavior may
  require additional rate verification or a stronger uncertainty-propagation
  method; mechanism-coverage and applicability gates therefore remain
  independent and fail closed.
- Verification cost estimates influence allocation only when every required
  route uses the same declared units and cost model. Wall time, convergence
  probability, queue behavior, and scientific value are not yet learned or
  predicted. A low-cost score cannot bypass a required upstream route.
- Even `verified_t95_converged` is conditional on the declared gas model,
  reaction-network coverage, electronic-structure method, Arkane model, and
  retention criterion. It is not a laboratory shelf-life, handling, or safety
  guarantee.
- A `conditional_model_prediction` is weaker still: its reverse-return bound
  is an approximation-error certificate between explicitly named fixed models,
  not a physical confidence interval. The bounded handoff supports only closed
  isothermal conservative first-order pure-parent networks and finite-horizon
  parent retention. Unknown conformers/return channels, rate calibration,
  pressure applicability and mechanism completeness remain unresolved. An
  expanded repaired model requires a new model-scoped assessment; the existing
  certificate cannot be carried across new species/channels. No t95 or lifetime
  is admitted from this result. This extension does not add solution, surface,
  photochemical or other environment applicability.
- Expanded retained-network scenario diagnostics now distinguish operative
  forward/reverse rates, electronic states and return/escape competition, but
  they do not supply physical uncertainty bounds. Removing a spin-changing
  branch or one of its directed channels is a declared counterfactual, not an
  uncertainty sample; weak nominal flux or insensitive tested scenarios cannot
  certify an unresolved branch irrelevant. The scope remains conservative,
  closed, isothermal first-order gas models. Numerical agreement is checked
  separately from reverse omission and conversion-model differences, and is
  not total physical error. The state audit recognizes a need for nonadiabatic
  evidence but does not calculate a spin-crossing rate. Pressure adequacy,
  conformer/counting accuracy, calibrated rate uncertainty and mechanism
  coverage remain independent gates.
- Decomposition storyboards are topology diagrams. Prepared animations are
  interpolated endpoint previews, relaxed scans are one-dimensional
  ground-state paths, and only a one-imaginary-mode TS plus endpoint-matched IRC
  is labelled a verified transition-state path. A barrierless animation shows
  a reproducible computed coordinate, not molecular dynamics or irreversibility.
- A computational-light route followed with a ground-state dissociation scan is
  not a simulation of excited-state branching, fluorescence, internal
  conversion, intersystem crossing, conical intersections, or nonadiabatic
  dynamics. Radical-pair dissociation also depends on the declared total-spin
  surface.
- The UV/Vis screen is a single-conformer, vertical full-LR-TDDFT calculation
  with modelled broadening. Solvent shifts, conformer averaging, vibronic
  structure, quantum yields, and nonradiative branching are not calculated.
- The `reactivity` command reports molecule-level proxies derived from ORCA
  Kohn-Sham frontier orbitals. These are not measured IPs, EAs, redox
  potentials, atom-resolved reactivity maps, or selectivity predictions.
- `structure` writes an RDKit/UFF initial geometry. It is not an optimized
  structure and must not be treated as an experimental or quantum geometry.
- `enrich` deliberately accesses an external source only when invoked. Its
  returned identifiers and properties are source-labelled database records, not
  verified identity, safety, or experimental-property determinations.
- Completed TS/frequency reuse after native008's IRC timeout only prepares a
  new bidirectional IRC in a derived workspace. Its source contracts, context,
  modes, geometry and stopped leases are checked; the old partial trace gives
  no connection credit. Real new IRC, rate/repair, full-chain resume equivalence,
  calibration and coverage remain unproved. The original campaign remainder
  is 362.6558744990034 monotonic seconds; larger ceilings require authorization.
  Epoch lease durations and monotonic campaign times differ and are not CPU
  time; clock-domain discrepancy is retained, not silently corrected.
- A subsequently completed IRC can now be reused as tracing evidence in a new
  workspace, with actual reached-endpoint Opt/Freq refinement. Conformer
  mismatch is not automatically wrong chemistry, but its rescue currently
  requires supported neutral closed-shell coordinate-graph corroboration.
  This is not general electronic-state determination, GOAT/global conformer
  coverage, or validation of a thermal conformer ensemble. Distinct/unresolved
  endpoint conformer bases may yield conditional Arkane rates while explicitly
  withholding RMG reinjection and a network lifetime. Matching coordinates do
  not establish correct absolute RMG species thermo or calibrated rate accuracy.
- The 2026-10-07 OH + H2O2 primary-reference lead remains unaudited for exact
  finite-pressure conditions and table/error definitions. The impossible NIST
  pressure field and incomplete reproduced table raster are not execution or
  scoring inputs. No additional case, holdout or rate-error calibration is
  admitted by that curation checkpoint.
- A separate source-audited primary OH + methane rate slice is descriptive only:
  its actual prediction is missing, supplement sigma/covariance/total uncertainty
  are unresolved, and all fourteen points share one study. It supplies no HPL
  applicability, full-fate validation, prospective holdout or STORCA calibration.
  The shared panel is unchanged; a useful numeric reference is not evidence that
  the engine predicts it correctly.
- The producer forces optical-isomer count one. Native inference identifies the
  retained cyclobutene C2 TS candidate as having two optical configurations,
  giving a factor-two identical-mode partition difference with achiral endpoint
  counts unchanged. Mirror-channel/ensemble inventory and degeneracy have not
  been established; neither a blanket inferred count nor doubling the rate is
  justified. The diagnostic is unscaled statistical mechanics, not production
  kinetics or a calibrated uncertainty bound.
- Unbound/invalid `method_calibration` allocation artifacts now refuse admission,
  and malformed current IDs cannot be stringified into matches. Exact opaque
  condition binding alone does not authenticate physical conditions, method
  recipes, source evidence or calibration. Legacy unbound user/database bounds
  are still supplied scenarios and may be inapplicable or miscalibrated; this
  compatibility is not an automatic scientifically justified error producer.
