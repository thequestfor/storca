# Methods

`storca run` performs a geometry optimization followed by a harmonic frequency
calculation with ORCA. The frequency job uses the XYZ geometry written by the
optimization step, not the initial input geometry.

The default input generator uses B3LYP/def2-SVP. Charge and multiplicity are
explicit command options and are passed unchanged to both calculations. The
named `b3lyp-def2-svp` harmonic profile supplies the current 0.970 project
baseline scaling factor; it is recorded in every new spectrum run along with
whether the factor was profile-supplied or explicitly overridden. It is not a
universal literature guarantee. A run is considered a minimum when no parsed
vibrational frequency is below -20 cm^-1; low-magnitude modes are reported
separately because they may reflect rotations or numerical noise.

This project manages calculations and parses their outputs. It does not make
chemical safety, reactivity, or experimental-spectrum guarantees.

## Candidate harmonic profiles

The initial method bakeoff has two additional, intentionally unscaled ORCA
composite profiles: `b3lyp-3c` and `r2scan-3c`. ORCA Manual 6.1 section 3.6.5
describes B3LYP-3c as introduced for efficient gas-phase IR spectra; section
3.6.3 describes r2SCAN-3c as a broad thermochemistry/geometry/noncovalent
method. Those descriptions motivate testing, not a claim that either is more
accurate for this project's heterogeneous condensed-phase references. Each run
records the exact profile and scale-factor provenance in `metadata.json`.

## IR spectra

`storca spectrum` creates RDKit conformers, optimizes each retained conformer
with ORCA, and runs frequencies on the corresponding optimized geometry.
Unit-area Gaussian broadening is applied to positive IR modes, so changing FWHM
does not change integrated band strength. Conformer weights are computed from
successfully completed optimized electronic energies at the selected
temperature, including GOAT-reported degeneracy when available; they are not
Gibbs free-energy populations. Conformers with a
significant imaginary mode are excluded and successful weights are
renormalized.
The default 0.970 harmonic scaling factor is a practical starting point, not a
universal correction; validate it for the chosen method and chemical class.
`spectrum-analyze` performs transparent peak matching: it reports position and
relative-intensity differences, but does not shift, fit, or otherwise modify
either trace. The match window and prominence threshold are stored in its JSON
report.

### Calculation-first practical display

`storca spectrum --spectrum-model practical` first calculates exactly the same
ORCA conformer ensemble as the raw mode. It then applies the versioned
`ambient-organic-v0.3` rule profile to those calculated bands. The initial rules
use visible RDKit descriptors to broaden selected regions. At present, only
the flexible-molecule fingerprint rule is enabled. The initial O-H and N-H
rules are not applied automatically: their limited benchmarks conflicted.
Raw artifacts are always retained.

This is not a model that predicts spectra from SMILES and it is not a
substitute for an explicit environment calculation. It is a small, empirical
calibration layer whose parameters must earn their place through held-out
experimental benchmarks.

Experimental runs serialize a versioned condition contract containing the
resolved phase, temperature, pressure, composition, solvent, concentration,
measurement geometry, resolution, apodization, and available optical geometry.
Unknown fields remain null or absent and cause condition-matched evaluation to
fail closed. The calculated absorption-strength array is retained as
`spectrum_intrinsic.csv` before ATR or instrument transfer; `spectrum.csv` is
the separately rendered measurement.

The experimental FTIR model can additionally evaluate a bounded set of neutral
dimer geometries and, for bounded small protic molecules, linear and cyclic
trimer geometries for eligible neat condensed phases. Dimer and trimer
intensities are divided by two and three respectively. A per-band
environment-sufficiency gate requires at least three independent geometry
clusters, positive frequency variance, effective sample size of at least 2.5,
documented interaction-geometry diversity, and displacement overlap of at
least 0.45 covering at least 80% of the band population. A failed distribution
is retained diagnostically but collapsed to its weighted center before the
small residual linewidth and instrument response are applied. Passing these
minimum-evidence checks does not establish thermodynamic sampling or liquid
convergence. Larger aggregates, unbiased bulk occupancies, and representative
DFT correction transfer remain validation-gated extensions.
Statistics are clustered by independent environment before effective sample
size and geometric diversity are calculated. Local oscillators are separated
into hydrogen-bonded and non-donating spectral classes, with donor-only,
acceptor-only, donor-acceptor, and free coordination retained diagnostically.

The current restrained-xTB stage deterministically proposes 40 (`auto`) or 75
(`balanced`) dimer/trimer poses from the retained monomer atom ordering. It spans
the configured hydrogen-bond distance, angle, axial-rotation, and framework
coordinates, then runs an isolated GFN2-xTB loose optimization with distance,
angle, and, where well-defined, dihedral bias potentials. Each candidate is
content-hashed and resumable. Post-optimization validation checks atom order,
intramolecular connectivity, collisions, association distance, and departure
from its requested sampling stratum. xTB electronic energies are retained only
as sampling features: the stratified seed poses receive no population weights.
Four (`auto`) or six (`balanced`) diverse poses seed restrained NVT GFN2-xTB
trajectories at the spectrum temperature (298.15 K by default), for 2 or 4 ps
respectively. A 0.5 fs step and hydrogen mass 4 are used without SHAKE;
coordinates are dumped every 10 fs. Emergency-exit text or fewer than 95% of
the expected frames fails the trajectory even if xTB returns success. Partial
seed completion is rejected rather than renormalized. The first 25% is discarded,
and an initial-positive-sequence statistical-inefficiency estimate sets the
minimum extraction stride. At most 40/75 decorrelated frames are retained and
assigned equal-time occupancy. These occupancies are conditional on the finite
restrained multi-seed protocol and are not labelled as unbiased bulk-liquid
populations. A failed trajectory stage records the downgrade and falls back to
the static coverage ensemble. Near-duplicates are removed using fixed
interaction-coordinate, local-frequency, and heavy-atom-RMSD thresholds. The
representative acquisition controller identifies prevalent local
bond/coordination classes and first seeks three independent environments for
each class. Remaining jobs maximize local-frequency, interaction-geometry,
topology, and association diversity. Transfer validation after each cumulative
ORCA batch reprioritizes pending representatives that contain failed classes.
The `auto` and `balanced` hard caps remain four and six environments. Each
representative receives the summed occupancy of the decorrelated frames mapped
to its cluster; duplicate removal preserves the removed frames' occupancy
mass. The global ORCA cap reserves their jobs before allocating monomer
jobs. ORCA evaluates a snapshot Hessian at the selected restrained-xTB geometry
without an unconstrained DFT optimization; low-frequency nonstationarity is
retained diagnostically and the result is not labelled an unconstrained DFT
minimum. Representatives are submitted central-first and then by farthest-point
diversity in cumulative batches (`2,3,4` for `auto`; `2,4,6` for `balanced`).
Important band centers and raw distribution widths, integrated intensity, mode
overlap, new significant bands, and a fixed-width whole-spectrum cosine overlap
are compared after every batch. Two consecutive passing comparisons are
required. A sufficient distribution that reaches the hard ORCA budget first is
retained as provisional and explicitly labelled unconverged.
Every retained xTB geometry also receives a resumable unrestrained GFN2-xTB
numerical Hessian. Gaussian-style xTB normal coordinates are parsed and
projected onto the same local X--H coordinates and internal-coordinate
fingerprints as ORCA. At representative geometries, STORCA forms additive
frequency corrections and validates them by leaving one independent
representative out. A mode class fails closed unless it has at least three
representatives, similarity at least 0.50, at least 20% improvement over raw
xTB, high-frequency MAE at most 30 cm^-1, maximum error at most 75 cm^-1, and
applicability coverage. Transfer uncertainty combines withheld error,
correction spread, and environment distance. Cluster intensities are divided
by cluster size, but xTB intensity transfer remains disabled pending its own
validation. Restrained-trajectory occupancies are explicitly not asserted to
be unbiased liquid populations.
Transferred convergence treats one candidate as the independent statistical
unit per mode class, even when a trimer contributes multiple local stretches.
It must pass two consecutive batches under three deterministic stratified
orders; candidate bootstrap intervals are reported for centers and equivalent
FWHM values. A hybrid diagnostic combines corrected local frequencies with
mean representative-ORCA intensity for the same class and retains
representative ORCA nonlocal modes. It becomes the display basis only if
validation, order-robust convergence, and local-intensity coverage all pass.
For equal atom counts, normal modes are projected onto an inferred
internal-coordinate basis containing bond stretches, angle bends, torsions,
and chordless-ring deformations. Normalized coordinate-participation
fingerprints provide the primary assignment, with frequency as a weak guard
and geometry-aligned Cartesian overlap as a fallback. Chemical stretch labels
include C=O, C-O, C-N, C-H, O-H, and N-H where the geometry supports them.
Across differently sized dimer and trimer Hessians, X--H modes still use the
validated target-local bond projection; generalized cross-size fingerprint
matching remains incomplete. Multiple molecules in one cluster do not count
as independent environment configurations.
Snapshot Hessians with modes below -50 cm^-1 are labelled poor-stationarity and
their full Hessian is diagnostic-only; this does not automatically discard a
high-confidence target-local X--H projection. Retained legacy jobs without a
gradient contract explicitly report that gradient reliability is unavailable.
The environment-preserving refinement contract applies interaction restraints
only during optimization. Its separate unrestrained gradient must pass RMS
3e-4 and maximum-component 1e-3 hartree/bohr gates before full-Hessian use is
permitted; otherwise local-mode finite differences are required.
The finite-difference implementation displaces the heavy atom and hydrogen in
opposite mass-weighted directions so their pair center of mass is unchanged.
Analytic ORCA gradients and dipoles are evaluated at plus and minus 0.005 and
0.010 Angstrom displacements. The projected curvature and dipole derivative
are Richardson-extrapolated; the two raw frequencies must agree within 20
cm^-1. Each physical ORCA process attempt is recorded against a stage-local
hard invocation cap. Validated finite-difference frequencies and intensities
override only the matching bond-local snapshot-Hessian pair; the retained
Hessian continues to supply nonlocal modes. Stationary validation compares the
center and summed intensity of the complete localized stretch set with the
corresponding coupled normal-mode subspace, rather than comparing arbitrary
individual symmetric or antisymmetric modes.
The environment fallback planner selects at most one oscillator of a given
mode class from each representative, preventing a cyclic trimer from inflating
independent coverage. It allocates complete three-representative class blocks
only; an allowance too small for a complete block is left unspent. Local jobs
share a namespaced run-level ledger so equal bond indices in different
representatives remain distinct process invocations.
Static xTB ensembles can be extended in independently seeded 20-candidate
acquisition rounds. Each round preserves dimer, linear/cyclic trimer, and
strong/intermediate/weak association coverage. Convergence compares real
cumulative rounds rather than arbitrary prefixes: two consecutive comparisons
must pass the center, width, whole-spectrum-overlap, and new-class gates.
Additionally, every final band must have center and FWHM bootstrap 95% interval
half-widths no larger than 15 and 25 cm^-1 respectively. Ordering permutations
remain sensitivity diagnostics rather than substitutes for new observations.

### Held-out harmonic scaling

`storca calibrate-ir-scale` is a transparent global-scale search over completed
ORCA normal modes. It holds FWHM fixed, chooses a factor only from explicitly
named training molecules by prominent-band position error, and then reports a
separate holdout score. It never invokes ORCA, edits an existing run, or
automatically changes the default profile. This guards against improving a
small, heterogeneous benchmark by overfitting its digitized references.
Schema-v2 reference entries retain a content hash, provenance, condition
contract, and fixed dataset partition. One canonical chemical identity cannot
occur in more than one locked partition. Legacy schema-v1 references and v2
`development_unassigned` entries remain usable for development diagnostics but
do not constitute a held-out validation claim.

## Stability condition models

### Multi-lane campaign contract and gas/solution execution adapters

`storca stability-campaign` schema version 1 normalizes a shared molecular
system, requested observables, multiple lane instances, lane-specific
environments and methods, evidence references, user applicability declarations,
and local/global compute ceilings. RDKit canonicalizes every component,
removes atom-map labels from shared identity, and requires the declared formal
charge to match the molecular graph. The shared system contains identity and
role only; mole fractions, concentrations, solvents, pressure, illumination,
and other thermodynamic assumptions remain within their owning lane.

The scientific campaign hash includes normalized molecular states, observables,
environments, method profiles, budgets, user assumptions, and evidence content
hashes. File locations, display names, and notes are excluded. Reusing a
campaign directory with a changed scientific hash fails before planning;
changed descriptive metadata is retained in a separate history. The existing
gas `ConditionSpec`, implicit-solution environment, solution verification
session, and legacy `CampaignSpec` are not weakened or rehashed.

Without `--execute`, the command is preparation-only. It consults the runtime lane registry,
performs environment preflight, records required inputs and capabilities, and
writes independent lane/observable admissions. It reports declared-only,
unsupported-observable, missing-input, method-mismatch, and unresolved
user-applicability states without launching ORCA, xTB, RMG, Arkane, or Cantera.
It does not estimate a calculation count before a lane adapter resolves a
candidate/evidence inventory, and it emits no observable claim. Gas and
solution plans are conditional records and are never combined into a universal
lifetime.

Live execution is doubly opt-in: the immutable request records
`execution_policy: explicit_execute`, and the CLI must receive `--execute`.
Schema-v1 execution selects only fully admitted homogeneous-gas or restricted
neutral implicit-solution `t95` lane requests. Every adapter retains an
environment-, method-, observable-, budget-, evidence-, and software-bound
execution contract. A changed contract cannot reuse a directory. Gas stages
retain their append-only checkpoints. An exact completed solution-lane rerun
reads the retained terminal result; a solution verification session additionally
retains its existing append-only allocation and repair history.

The adapter translates each campaign lane into the existing gas stage boundary
and calls the unchanged stability screen. RMG remains mandatory; bounded
graph/xTB discovery can be enabled once for the shared target and reused by all
gas conditions. Existing endpoint, NEB/NEB-TS, TS/frequency/IRC, barrierless,
Arkane, verified-library reinjection, independent propagation, and convergence
gates remain authoritative. `worst_case_search` stops only after an explicit
condition-matched verified and converged `t95`; a numerical RMG screening
lifetime is insufficient. `environment_map` executes every selected gas lane.
Only a lane result with `observable_status: verified` creates a condition-scoped
claim. A verified gas result and a verified solution result therefore remain two
conditional claims even when they share the same target and observable. The
report explicitly declines cross-environment aggregation. Failure of one
adapter does not erase a verified claim from another; the overall campaign is
then completed with unresolved lanes rather than advertised as fully resolved.

Accepted gas ceilings are deliberately named by their actual scope:
`rmg_attempt_timeout_seconds`, `rmg_attempts`, `rmg_max_processes`,
`rmg_initial_max_iterations`, `rmg_initial_max_edge_species`,
`rmg_max_objects_per_iteration`, `rmg_maximum_heavy_atoms`,
`rmg_heartbeat_timeout_seconds`, `verification_max_iterations`,
`verification_path_timeout_seconds`, `verification_orientations`, and
`verification_neb_images`. “Initial” RMG bounds may grow only through the
existing declared restart schedule; the per-attempt supervisor timeout remains
fixed. The only schema-v1 global gas ceiling is `maximum_gas_stages`. Requests
using broader but unenforced quantities such as estimated core-hours fail before
chemistry is launched. Restricted solution execution accepts `profile_points`
for direct verified-network propagation, or the exact pair
`maximum_orca_paths` and `maximum_estimated_core_hours` for a verification
session. Global solution ceilings are `maximum_solution_lanes` and
`maximum_solution_estimated_core_hours`; the latter sums the declared session
budgets before any lane launches.

The solution adapter has two mutually exclusive evidence modes. Direct mode
requires exactly one `solution_environment` artifact and one
`verified_solution_network` artifact. Their content hashes, target identity,
environment, method, and coverage-claim identity must match the campaign before
the existing restricted propagator runs; it performs no electronic-structure
calculation. Session mode requires exactly one
`solution_verification_session_request`. Its target, implicit-solvent contract,
method profile, observable, and compute budget must exactly match the lane. The
existing observable-conditioned session then verifies one allocated path at a
time and repropagates the repaired network. This session route remains limited
to 298.15 K. Mixing both evidence modes, supplying extra evidence kinds, using
additional mixture components, or requesting a charged/open-shell target fails
closed under the current solution applicability contract.

Every stability stage is bound to a serialized condition contract containing
the scenario and phase model, temperature, pressure, composition, target mole
fraction, simulation duration, retention fraction, light condition, and
radical-source assumptions. Repaired calculations must return the identical
contract hash; a changed condition cannot inherit the original stage's
verification.

The first ladder stage is `low-pressure-intrinsic-gas-screen`: 1% target in a
nonreactive nitrogen bath at `1e-6` bar, propagated by RMG as a homogeneous gas.
It is a low-collision gas surrogate, not a literal vacuum or isolated-molecule
dynamics calculation. Uniform single-bond homolysis enumeration may be retained
as diagnostic candidate information elsewhere, but it is not this stage's
kinetic model and cannot supply its lifetime. Later dark stages use the same RMG
framework for dry nitrogen, dry air, and humid gas at the declared conditions.
Humid gas contains water vapour; it is not liquid water or a droplet model.
The parent isolated-molecule optimization/frequency dossier is reused between
ladder stages only after exact identity, charge, multiplicity, method, element
inventory, mode coverage, and artifact checks. Reaction-path and kinetic
evidence is never reused across different condition contracts.

The ladder has two explicit campaign policies. `worst_case_search` preserves
the bounded default: after the first condition-bound, verified, converged
`t95`, later rungs are recorded as not run. `environment_map` executes every
requested rung even after a verified loss event, and reports coverage,
condition-specific unresolved evidence, and all verified `t95` records. The
mode changes stopping semantics only; it does not relax any route, rate, or
propagation gate.

When the graph/xTB discovery lane is requested for a ladder, STORCA runs a
bounded search within each gas condition. The target-unimolecular edit space is
retained, while conditions with explicitly reactive species additionally search
pairwise radical association and H-atom abstraction using declared partner
charge and multiplicity. Candidate total-spin surfaces are retained only when
reactant and product spin-coupling sets intersect; ambiguous partner states
remain unresolved. An immutable discovery-assembly contract binds the canonical
target, target state, reactive-partner identities and states, full condition
contract, resolved graph/xTB budget, and RMG bounds/libraries. RMG generation,
ORCA path verification, Arkane rate evidence, repaired-library injection, and
network repropagation remain condition-specific.

### Solvent-bound graph/xTB screening and neutral four-pi proposals

The schema-v3 graph/xTB backend can screen a candidate requested in a named
implicit solvent. The requested SMD/CPCM path environment remains attached to
the candidate; the cheaper calculation separately records observed GFN2-xTB
ALPB screening with an explicit `gsolv` (1 M gas/solution) reference. ALPB totals
include solvation terms, not just electronic energy. They are endpoint screening
quantities, not SMD/CPCM thermal free energies, activation barriers or rates.
Unsupported solvent parameters, inconsistent phase/model declarations, missing
native solvent evidence or gas output cannot fall back to a gas calculation.

An added bounded edit proposes isolated ordinary-carbon cyclobutene ring opening
and the reverse acyclic conjugated-diene closure while preserving mapped atom,
hydrogen, valence, charge and radical inventories. Specified reacting
stereochemistry, fused rings, inappropriate ring sizes, heteroatom/isotopic
reacting sites and incompatible hybridization are outside this narrow edit.
This is a connectivity hypothesis: orbital allowedness, conrotatory branches,
saddles and kinetic relevance still require separate evidence.

Native admission checks optimization convergence, declared optimized
connectivity, GFN2 method/version/invocation, ALPB solvent/reference, and
stdout/JSON total-energy rounding consistency. Exact cache reuse rechecks
source hashes and the observed screen; endpoint and assembly contracts bind
environment as well as target/state/budget. Changed or incomplete contracts are
refused in place, not overwritten. Historical gas artifacts are preserved and
cannot be relabeled as solvent-screen evidence.

An actual native-candidate handoff exposed a vocabulary mismatch in the shared
path gate: `phase="solution"` was not recognized as implicit solution. The
canonical alias is now accepted only under the existing named-solvent,
SMD/CPCM and exact method-keyword checks. This repair does not grant new
solvent physics, waive applicability or prove a path.

The retained three-case native diagnostic recovered the expected neutral
connectivity proposals for cyclobutene opening and diene closure in water, and
cyclobutene opening in acetonitrile. It used six actual semiempirical xTB
endpoint jobs at one core, with 0.924418 seconds recorded batch elapsed time.
Exact cache checks preserved all 24 files/hashes/timestamps per case and
attempted no additional endpoint jobs. No ORCA, rate or propagation calculation
was performed. These are development screening checks, not mechanism recall,
coverage calibration, a positive verified path or a lifetime benchmark. See
[screening evidence](runs/solution-discovery-screen-diagnostic-001/EVIDENCE_2026-10-06.md).
The final 2026-10-06 software checkpoint passes 1018 application tests in
151.461 seconds (22 expected skips), 115 expanded native/discovery-boundary
tests in 9.695 seconds and 179 affected tests in 1.234 seconds. These software
results add no physical prediction acceptance credit.

The assembly retains backend-specific results plus one source-neutral candidate
pool. It explicitly lists source agreement, RMG-only, graph-only, excluded by
the declared xTB screen, and unresolved candidate/assertion identities. RMG
routes without resolved endpoints and graph candidates whose xTB screen failed
are not silently discarded. Backend exceptions, unavailability, unresolved
assertions, and search truncation produce incomplete or fail-closed assembly
states. If the requested assembly is not `completed_declared_search`, the
existing route/rate verifier may retain its model result, but STORCA vetoes a
verified lifetime because discovery coverage under the requested policy remains
unresolved. Completion is scoped to the configured searches and never proves
mechanism completeness.

When partial RMG enumeration exhausts its direct verification queue, the
engine can still attempt a retained graph/xTB survivor whose reactants all have
positive initial composition in the same condition. The candidate must consume
the target, retain passed atom balance and endpoint records, and match any
retained verification condition ID. This fallback supplies path evidence only;
incomplete RMG propagation continues to block a lifetime. A successful path and
rate still require repaired-library insertion and a completed model rerun.
Deferred RMG routes remain unresolved, and the planner distinguishes an empty
direct queue from one exhausted by verified or deferred attempts.

A completed xTB calculation can also place a candidate above the configured
endpoint reaction-energy window. Such a candidate is retained as
`xtb_screened_out`; if no independent backend retains the same native reaction,
the exclusion blocks an observable claim.
GFN2-xTB endpoint energy is neither an activation barrier nor evidence that a
high-energy intermediate cannot feed a consequential downstream path. The
block currently has no runtime bypass. A future relaxation would require a
hash-bound calibration artifact tied to the exact screen, observable,
conditions, environment, and chemical scope; generic discovery calibration is
not sufficient. The current resolution is to expand/remove the screen or verify
the candidate. If RMG independently retained the same reaction, the xTB result
is preserved as source disagreement but does not remove or separately block the
RMG-retained candidate.

Each default stage is also represented by a versioned stability-lane contract.
Lane registration, runtime admission, execution, and observable authorization
are separate decisions. Registration records an environment and intended
method chain. Admission requires an exact environment-kind match, an allowed
observable, and every declared scientific capability. A lane marked
`declared_only` or missing a capability cannot retain a verified observable;
the standardized lane result enforces this independently of stage-specific
assessment text.

The six current gas stages are adapters over the existing implementations:
bounded RMG discovery, optional graph/xTB hypotheses, ORCA path verification,
Arkane thermal rates, verified-library reinjection, independent Cantera gas
propagation, and observable robustness gates. The lane layer therefore adds
scope and auditability without substituting a new kinetic model. The first
implicit-solution lane is executable only for retained verified,
neutral first-order networks; rates may be externally supplied or generated by
the restricted STORCA solution-path producer. Microsolvation, surface, and
electrochemical lanes remain `declared_only` until their own physical models
and validation evidence exist.

## Restricted implicit-solution concentration lane

The `solution-neutral-implicit-v1` executor consumes an immutable solvent,
solvation-model, temperature, duration, retention, volume, dark-condition, and
one-molar standard-state contract. Every species has a canonicalized SMILES;
RDKit-derived elemental composition and formal charge must match its declared
bookkeeping, and the network target must match the bundle target SMILES. It
accepts CPCM or SMD as identities of retained rate evidence. Every species
supplies elemental composition, charge, and initial concentration. Every
reaction has exactly one unit-stoichiometry
reactant, positive integer products, a condition-matched rate in `s^-1`, a
source, and a rate claim ID. Element and charge balance are validated before
propagation. Every target-consuming rate also requires explicit multiplicative
lower and upper bounds and their provenance. The lower bounds and upper bounds
are combined across all parallel target-loss reactions to retain a t95 interval.
If one bound crosses the requested horizon and the other does not, the
observable remains unresolved.

An experimental rate must be identified as a verified experimental solution
rate. A solution Eyring/TST rate additionally requires a verified path claim,
explicit one-molar standard-state correction, and retained calculation hashes.
This validates the evidence interface, not the correctness of an external
calculation. Network coverage is separately required as
`verified_complete_within_declared_unimolecular_scope`; any unresolved candidate
blocks the observable.

For an ordinary neutral unimolecular candidate, `solution-path-verify` can
produce the path and rate claims directly. It normalizes the candidate to the
existing engine-neutral `RouteSpec` and reuses the generic endpoint, NEB-TS,
explicit TS refinement, frequency, and bidirectional IRC workflow. The selected
method profile is augmented with exactly one `SMD(solvent)` or `CPCM(solvent)`
keyword. A post-run audit requires that exact token in every retained ORCA
input and rejects missing, mixed, or `%cpcm`-overridden inputs. Environment
profile, solvent, model, target identity, stoichiometry, charge, and route
identities are immutable preflight gates.

The solution producer requests one explicit thermochemistry temperature through
ORCA 6.1's `%freq Temp` control. That temperature is bound to the route context,
endpoint/workflow metadata, immutable request and session cache identity;
nested endpoint recovery and TS frequency jobs inherit it. Default gas/P1
inputs remain unchanged when this optional contract is absent. Rate conversion
requires the exact frequency input and one normally terminated output section
at the requested temperature. The output-line tolerance is 0.00500001 K for
two-decimal printing, not physical uncertainty. A missing request, conflicting
cache, wrong output temperature or unsupported multi-temperature record is
refused without borrowing another temperature's evidence. After a verified
first-order saddle and matched observed free-energy recipes, STORCA computes

```text
k(T) = (k_B T / h) exp[-(G_TS - G_R)/(R T)]
```

with transmission coefficient one. The 1 M concentration correction is zero
for the unimolecular reactant-to-TS activation step because both sides have
molecularity one; that cancellation is retained explicitly. A user-declared
composite free-energy-barrier bound becomes deterministic lower/upper rate
multipliers. Required uncertainty components force conformer,
electronic-structure, implicit-solvent, and harmonic-TST uncertainty to be
acknowledged, but do not calibrate those terms. Content hashes link the rate
claim to all retained ORCA path artifacts. The emitted bundle is validated and
immediately repropagated by the existing solution lane.

No bundle is emitted for a nonpositive activation free energy, incomplete or
temperature- or recipe-mismatched thermochemistry, barrierless path, unverified IRC,
active method-escalation flag, mixed solvent input, incomplete network
coverage, or unresolved discovery candidate. Path verification never creates
a reaction-network completeness claim.

The solution producer uses `observed_orca_6_1_thermochemistry_v1` to retain the
actual output-selected RRHO or ORCA Grimme quasi-RRHO vibrational entropy,
reference frequency when active, cutoff, applied frequency scaling, pressure
and program version. Exact input/output hashes bind each recipe; explicit input
overrides must agree with observed scalars at printed rounding precision.
Recipe identities include those common settings, not species-specific
rotational symmetry, which is retained separately. Missing, contradictory,
unprojected-mode or out-of-parser-domain evidence blocks rate admission.
Temperature/recipe parsing and final admission reject detected source mutation.
Session candidate cache identity binds the recipe admission contract.

The installed retained ORCA 6.1.0 parent prints quasi-RRHO entropy with a
100 cm^-1 reference, cutoff 1 cm^-1, frequency scale 1.0 and pressure 1 atm.
No quantum input or free-energy total is changed by observing this recipe.
The gas Arkane reconstruction, including its own scale factors/statistical
mechanics, is unchanged. The legacy uncertainty component `harmonic_tst`
remains an input name for supplied conventional-TST model bounds, not proof
that the ORCA entropy is pure RRHO or that its error is calibrated.
Quasi-RRHO is not explicit rotor, conformer or anharmonic sampling; ideal-gas
translational/rotational statistics combined with implicit-solvent energies
remain a model assumption. See the [ORCA thermochemistry documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html).

The observable-conditioned solution verification session places a speculative
`SolutionCandidateNetwork` in front of this path producer. Candidate records
retain discovery evidence, applicability flags, a condition-bearing approximate
rate interval, its permitted use, verification state, and estimated cost. The
strict `VerifiedSolutionNetwork` is not relaxed: a candidate crosses that
boundary only when the path producer returns a `VerifiedSolutionReaction`.
Resolved RMG and STORCA-native/graph candidates are adapted by canonical
molecular identity. Gas rates and graph/xTB screening remain allocation-only
priors; endpoint energies are not converted into barriers.

### Conditional candidate coefficients: in-progress producer

A separate internal `implicit_solution_rate_coefficient_request` can calculate
the forward and reverse coefficients of an already selected neutral singlet
candidate, with one or two separated molecules on either side. It needs the
candidate identity/stoichiometry, named CPCM/SMD solvent, temperature and method
contract, not a supplied verified network, initial concentrations, duration or
fabricated barrier-error bound. `run_solution_candidate_rate_coefficients`
adapts a retained solvent-bound native discovery candidate into that request;
gas/unspecified or differently bound discovery evidence is not silently relabeled.
This is an internal discovery-to-calculation boundary, not a replacement product
based on manual pathway/rate imports and not autonomous solution discovery.

The existing `solution-path-verify` command handles this request without new
flags and reuses the native endpoint → NEB/TS → frequency → IRC workflow. After
the existing path gates pass, the producer reparses condition-/recipe-matched
thermochemistry for the shared TS and both channel sides. Two-molecule sides
use the sum of validated separated-species free energies, not an encounter
complex; repeated stoichiometric copies count repeatedly even if one species
calculation is reused. One-molecule sides retain the selected path endpoint.
The observed thermochemistry reference pressure is not a new operating-pressure
input or evidence that a solution pressure effect was modeled.

For a side with molecularity n, the conversion to the fixed 1 mol/L standard
state is explicit:

```text
delta_species = R*T*ln(1000*C_standard*R*T / (p_thermochemistry_atm*101325))
DeltaG_dagger_C = G_TS - sum(G_separated_i) + (1-n)*delta_species
k_n = (k_B*T/h) * C_standard^(1-n) * exp(-DeltaG_dagger_C/(R*T))
```

Here C_standard is expressed in mol/L. Coefficients have units `s^-1` for n=1
or `L mol^-1 s^-1` for n=2 and use elementary-reaction-event mass action, without
an extra identical-reactant depletion factor. Both directions share one TS and
the same thermochemistry/standard-state contract; their dimensioned ratio is
checked for detailed-balance consistency. Nonpositive corrected barriers,
wrong solvent/method/state/temperature/recipe, source mutation or failed path
evidence withhold coefficient packets.

Successful packets are `estimated_with_unresolved_uncertainty`, with
`calibration_required`, null error multipliers and
`network_rate_admission_allowed=false`. They retain source hashes and separate
path/rate claim IDs, but emit no verified stability bundle and perform no
network propagation. They are selected-saddle conventional-TST scenarios with
transmission coefficient one, not summed conformer/path populations,
diffusion-corrected kinetics or calibrated confidence. Synthetic path-boundary,
dimensional-arithmetic, adapter and CLI tests exercise this implementation;
real positive solution paths, coefficient accuracy, automatic error evidence,
coverage and autonomous concentration/lifetime predictions remain open.

At the coefficient producer's 2026-10-05 checkpoint, 999 application tests passed in
157.101 seconds (22 expected native-only skips) and 68 focused native tests in
16.729 seconds. [Producer evidence](runs/solution-rate-coefficients-diagnostic-001/EVIDENCE_2026-10-05.md)
records commands, hashes and the retained negative boundary. This validates
software contracts, not positive solution paths, coefficient accuracy or lane
applicability; all scientific acceptance requirements above remain open.

### Bounded reactant conformers in solution rate production

`conformer_rates.py` is called by the existing solution path/rate producer,
not a standalone ensemble importer. Preparation and production share an
eligibility check: one neutral closed-shell singlet, ordinary isotope treatment,
at most 20 heavy atoms. A rotor/nonaromatic-ring trigger initiates ETKDGv3
sampling with explicit hydrogens, seeds 42/1729, 12 geometries per seed, one
thread, and finite embedding/UFF iteration caps. UFF energies only rank
unverified geometric seeds. Source-versioned immutable acquisition records
retain failed/incomplete embeddings and UFF relaxations, truncation, measured
sampling time and unresolved coverage. The 30-second seed allowance is checked
between bounded native operations; it is not a hard preemptive C++ timeout.
See the [RDKit embedding API](https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html).

Chirality-preserving graph automorphisms and proper-rotation all-atom RMSD
deduplicate seeds/minima, capped at 1,024 symmetry matches. Reflection is not
an equivalence operation. Seed/refined thresholds are 0.25/0.15 Å; near-duplicate
refined geometries differing by more than 1e-4 Hartree in G are unresolved.
These are protocol clustering settings, not calibrated physical error bounds.
Duplicate seed counts are not molecular degeneracies or populations.

Up to four selected seeds use the existing `_optimize_endpoint` executor,
minimum recovery, explicit temperature contract and inherited route-wide
budget. No per-conformer budget is created. Refined identity/stereochemistry,
composition/order, complete finite 3N frequency list including external modes,
minimum status, observed recipe and exact executed native input/reference
contracts are required. The new conformer producer does not borrow missing
job-contract ownership as legacy evidence. Immutable per-refinement source
records bind geometry, optimization/frequency inputs/outputs and job contracts;
resume rejects mutation and reuses completed evidence even after budget expiry.
The complete solution-input solvent audit is repeated after refinements.

If selected refinements collapse to one cluster without known seed failures,
unrefined selected seeds, selection truncation or a seed-wall stop, the original
single-path estimate may proceed under its externally justified sampling/model
bound. This does not prove unsampled or TS-conformer coverage. If distinct
reactant minima remain, interconversion and TS coverage are unresolved and no
aggregate reaction rate is admitted. A retained diagnostic may compute

```text
k_selected,conditional(T) = (k_B T/h) exp[-(G_TS-G_min)/(RT)]
                            / sum_i exp[-(G_i-G_min)/(RT)]
```

for represented reactant minima and the one selected verified TS. Rapid
equilibration and unit additional state degeneracy are explicit assumptions,
not inferred evidence. This scenario is neither the full reaction rate nor an
uncertainty bound. Non-equilibrated conformer networks, TS/product/encounter
ensembles, hindered rotors, unsampled-risk calibration and real rate benchmarks
remain future work. The rotor/ring trigger can miss OH/NH rotation, inversion,
tautomerism and other conformational risks.

Session state additionally binds its scientific production policy. Native
session-produced TST reactions with missing/changed temperature, recipe or
conformer-policy credit cannot support new prediction convergence. Refusal
preserves existing reports, bundles and state; exact user-declared imports
remain explicitly supplied evidence rather than being retroactively validated.
Policy consistency alone does not replay all underlying calculation evidence.

For the current no-target-reformation first-order model,

```text
k_loss = sum_i k_i
t95 = -ln(retention) / k_loss
d ln(t95) / d ln(k_i) = -k_i / k_loss
```

so the evaluator uses an analytical sensitivity rather than finite differences.
Non-target reactions have exactly zero t95 influence. Allocation applies
auditable priority classes first—unbounded influential rate, horizon-controlling
rate, non-admissible screening interval, then bounded residual uncertainty—and
uses the retained influence, interval width, discovery-source factor, and cost
only to order candidates within a class. Active explicit-solvent, spin,
photochemical, surface, ionic/speciation, diffusion, or multireference findings
are method-escalation outcomes and are not sent to the ordinary verifier.

One candidate is verified per append-only iteration. Successful rate evidence
repairs the candidate network; failed or barrierless results remain unresolved
and cannot disappear from convergence. Exact request hashes prevent a session
directory from being rebound to a changed budget or environment. Convergence
requires the declared t95 interval ratio, robust horizon classification,
verified direct target-loss rates, complete restricted discovery coverage, and
no unresolved influential method finding. Screening central values are never
copied into an unresolved scientific claim. A t95-only claim may converge while
downstream candidates remain, but a general concentration/product bundle is
then withheld because its product kinetics would be incomplete.

For `n` species and `m` first-order reactions, STORCA constructs the linear
concentration matrix `A` and an extent-rate matrix `B`. It propagates the
augmented system

```text
d/dt [c] = [A 0] [c]
     [x]   [B 0] [x]
```

with a deterministic scaling-and-squaring Padé(13) matrix exponential. This
produces species concentrations and integrated reaction extents from the same
solution. A monotonic target-loss contract permits a bounded bisection for the
first retention crossing. Elemental inventory, non-negativity, target-loss flux
closure, finite values, and the Padé linear-solve condition are terminal
numerical gates.

The lane is not admitted scientifically when it encounters ionic species,
bimolecular/diffusion-controlled chemistry, explicit solvent participation,
acid/base or ion-pair speciation, photochemistry, spin crossing, surfaces,
unresolved multireference character, target reformation, mismatched solvent or
temperature evidence, a non-one-molar standard state, or incomplete network
coverage. No gas-phase Arkane or Cantera rate is converted into a solution rate.

Solution bundles may be appended to `storca stability-ladder`. Target SMILES
identity is checked before any gas calculation starts. The additional lane
retains its own environment, evidence, and claim graph and is not represented
as an RMG gas stage. Campaign `environment_map` executes every requested bundle;
`worst_case_search` records later bundles as not run after the first verified
t95.

## Restricted implicit-solution mass-action lane

`solution-neutral-implicit-mass-action-v1` is a sibling of, not a schema
upgrade to, the linear solution lane. It consumes a distinct
`verified_implicit_solution_mass_action_network`, preserving interpretation of
older stored first-order artifacts. Each reaction has explicit reactant and
product stoichiometry, order one or two, condition-matched solvent and
temperature, a dimensional rate (`s^-1` or `L mol^-1 s^-1`), multiplicative
uncertainty bounds, a rate claim, and SHA-256-linked evidence. A second-order
rate additionally requires a sourced diffusion ceiling, and its upper rate
bound may not exceed that ceiling. Gas/RMG kinetics and rates with the wrong
dimensions are rejected rather than converted.

The supported environment is dark, neutral, constant-volume CPCM/SMD implicit
solution with a one-molar standard state and explicitly declared ideal molar
activities. Every species with nonzero initial concentration in a stability
campaign must match a component in the immutable `ChemicalSystem`. Element and
charge balance are checked reaction by reaction; ionic species, nonideal
activities, dynamic speciation, explicit solvent participation, photochemistry,
surfaces, spin crossing, and unresolved multireference flags are outside the
lane's applicability.

The nonlinear propagator evaluates elementary mass action,

```text
r_j = k_j product_i c_i^nu_ij
dc/dt = S r
dx_j/dt = r_j
```

where `S` is the signed stoichiometric matrix and `x` contains integrated
reaction extents. Target-producing directions are allowed; a reversible step is
represented by two separately evidenced rates whose stoichiometry and partner
IDs must be exact reverses. Adaptive RK4 step doubling is used for the ordinary
regime. If the initial concentration Jacobian's spectral rate multiplied by the
horizon exceeds the retained threshold, the run records that diagnostic and
uses adaptive linearly implicit Euler step doubling with a conditioned linear
solve. Step underflow, excessive rejection, ill conditioning, material negative
concentrations, or budget exhaustion is a numerical failure, not a lifetime.

Central propagation retains the first t95 crossing even when the target later
reforms, the full concentration profile, product gains/yields, dominant product,
reaction extents, and each reaction's net target contribution. Element
inventory, charge inventory, nonnegativity, extent reconstruction, and target
flux closure are terminal gates. For ten or fewer reactions every lower/upper
rate-bound corner is propagated. Larger networks receive joint-low/joint-high
and one-at-a-time scenarios. These are deterministic sensitivity scenarios, not
posterior confidence intervals; the reduced large-network set is explicitly not
a rigorous enclosure of all simultaneous rate combinations. A t95 claim is
authorized only when the central model crosses and every retained scenario
crosses within the horizon. Mixed crossing behavior is unresolved.

Lane schema v1 standardizes t95, retention at the declared duration,
component concentration, product yield relative to initial target, dominant
product, gross target-consumption branching fraction, and time-integrated
pathway extent. A single campaign lane request may issue several such explicit
observable claims from the same immutable evidence contract; each has its own
target/parameter semantics and rate-bound scenario evidence. Dominant product
is withheld if its identity changes across scenarios. `pathway_flux` currently
means time-integrated extent in `mol L^-1`, not an instantaneous rate.
Half-life, selectivity, and induction time are not yet standardized. Fixed
pH/speciation, nonideal activities, encounter-complex treatment, experimental
calibration, and autonomous solution reaction discovery remain open. The analytical example at
`examples/stability-campaign/solution-mass-action-executable.json` is software
and equation regression evidence, not experimental chemical validation.

RMG execution is accepted for kinetic interpretation only when retained
artifacts establish normal/recognized completion, requested time coverage, a
parseable mechanism and solver profile, physically usable kinetics, and
quantitative target-loss flux with acceptable numerical closure. Resource- or
coverage-limited models remain incomplete. Rates and `t95` values from the
initial RMG model are screening quantities until the verification procedure
below succeeds.

## Bounded independent graph/xTB candidate discovery

When explicitly enabled, STORCA retains a second candidate-discovery lane
without changing RMG's generation or kinetics machinery. The backend
enumerates bounded target single-bond dissociation, intramolecular formation,
one-bond break/one-bond formation, and atom-mapped hydrogen-transfer
hypotheses. For declared reactive gas partners it combines explicit-H mapped
graphs and attempts radical-radical association plus H abstraction in both
directions. This is a discrete connectivity search, not an encounter complex,
capture calculation, solvent model, or claim of bimolecular completeness.
RDKit structural validation checks molecule and
hydrogen inventories, total charge, sanitization/valence, unintended
fragmentation, identity with the reactant graph, and duplicate product graphs.
Each generated edit—including one rejected before xTB—is retained with a
specific status and reason.

Each structurally valid endpoint is seeded and optimised with a retained
GFN2-xTB calculation contract. The difference of summed endpoint electronic
energies is an optional, bounded screen for grossly endergonic proposals. It
is explicitly neither a transition-state barrier nor a rate coefficient. xTB
errors are candidate-level `xtb_failed` evidence; search/optimization/failure
budgets yield `discovery_truncated` or `search_truncated`, never an assertion
that no reaction exists.

Candidate identity is source-neutral and direction-independent for
deduplication, but retains the directed endpoints and atom maps required for
path verification. The merged record preserves RMG and graph/xTB discovery
assertions. RMG-only candidates continue through flux/dependency ranking.
Graph-only survivors remain `kinetic_relevance: unknown`; their deterministic
pre-queue is based only on retained endpoint-screen energy, not fabricated
kinetics. They can influence a condition observable only after the existing
TS/frequency/IRC, Arkane-rate, verified-library, and full-RMG-repropagation
contracts all succeed.

## Observable-conditioned discovery coverage

For every gas stability run, STORCA constructs a categorical coverage
assessment from the immutable assembly, its explicit backend scopes, the
requested observable, the exact condition contract, and any subsequent
verification summary. The assessment records a reaction-class matrix only
from classes explicitly declared by a backend. It does not infer that an RMG
family was searched from an equation or from the absence of candidates. An
RMG-only run is therefore retained but cannot silently skip the independent-
discovery requirement.

The assessment contains a symmetric disagreement ledger. A candidate is
labelled source agreement, source disagreement, or single-source outside the
other backend's declared scope; disagreement never elects a chemically correct
backend. Two completed empty searches are treated as unresolved shared-blind-
spot risk unless an applicable calibration policy explicitly says otherwise.
Reactive condition partners and open-shell targets add blockers when the
non-RMG search did not declare bimolecular or radical scope.

Graph-only candidates with no verified condition-applicable rate are classified
as `potentially_influential_unknown_rate`. Endpoint reaction energy alone is
not used to make them negligible. A candidate becomes
`verified_and_propagated` only when its route identifier is retained in a
verified rate replacement and the repaired network reaches an observable-
robust terminal state. Deferred candidates and survivors omitted by the
high-level verification budget remain blockers.

Backend failure, unavailability, truncation, unresolved candidate adaptation,
declared method-applicability gaps, changed candidate sets under an optional
bounded search expansion, and nonterminal propagation all remain distinct
blocking findings. `adequate_within_declared_scope` authorizes only the
condition- and observable-specific claim under that policy. Every assessment
sets `mechanism_completeness_claimed` to false, carries a content-derived
identity, and can be recomputed after rate evidence changes without rerunning
unchanged discovery calculations. Each recomputation is written as an immutable
content-addressed revision while `mechanism-coverage.json` remains the latest
pointer. Its uncertainty register isolates mechanism coverage from uncertainty
domains that this assessment does not evaluate.

## Flux-guided ORCA and Arkane verification

Independent gas propagation uses an **isothermal ideal-gas constant-pressure**
Cantera reactor, with energy explicitly disabled. The actual reactor temperature
and pressure are checked at initialization and every retained internal step;
the output records their observed ranges and the physical-model policy ID.
Flux planning, repaired-model admission and the no-loss shortcut require this
evidence to match the immutable condition contract. Legacy outputs without it
require repropagation of the retained mechanism; they cannot acquire an
isothermal stamp from their input temperature or successful process exit.
This repair changes neither elementary rate validation nor uncertainty gates.

Automatic verification is enabled by default for `storca stability` and
`storca stability-ladder`. When the completed RMG model contains target-loss
flux, the verifier ranks retained reactions by integrated target destruction.
Reachability dependencies are part of the ordering: if the controlling route
requires a generated intermediate, a reachable upstream route that forms that
intermediate is verified before the downstream consumer. Only one selected
route is changed per iteration so its effect can be attributed after
repropagation.

When explicit rate-uncertainty evidence is supplied, allocation can additionally
use the retained absolute logarithmic `t95` sensitivity. Each uncertainty record
must identify a route or reaction equation, bracket the retained rate with
`0 <= lower_multiplier <= 1 <= upper_multiplier`, and declare its source,
source kind, and applicability. The allocation proxy is
`|d ln(t95)/d ln(k)|` multiplied by the worst explicit `|ln(k/k0)|` bound.
It is a dimensionless first-order prioritization proxy, not a propagated
confidence interval. A zero lower multiplier remains an unbounded uncertainty.
When that uncertainty coverage is complete but sensitivities are absent,
STORCA estimates them with independent two-sided 0.5x/2x Cantera perturbations
of the retained reaction indices. Both perturbed propagations must complete and
cross `t95`; the complete diagnostic is rerun after each repaired mechanism.
Failure to obtain it is retained and restores the established ordering.

This policy activates only when every selectable unresolved nonzero-flux route
has both sensitivity and rate-uncertainty evidence. Otherwise the planner
retains the previous sensitivity/flux ordering and records the missing route
identities; absent evidence is never scored as zero uncertainty. Optional cost
normalization is applied only if all required routes provide positive estimates
with identical cost units and cost-model identity. Expected information gain is
recorded as `not_estimated` because STORCA does not yet have a calibrated model
for the probability that a calculation will resolve a route. Generated-species
dependencies remain hard constraints: a highly scored downstream route cannot
bypass verification of its reachable upstream producer.

For the selected reaction, the generic path layer retains the original RMG
reaction equation and species identities, resolved SMILES, per-species charges
and multiplicities, and candidate common total-spin surfaces. It requires one
explicitly selected common adiabatic surface. A supplied or uniquely resolved
atom map is extended to explicit hydrogens, including a uniquely defined
hydrogen transfer; ambiguous or search-exhausted maps stop before ORCA.

For selected RMG routes with retained native reaction indices and structural
artifacts, the default adapter first requests a bounded family-recipe atom
trace in the `--rmg-env` environment. It reloads the exact hashed Chemkin and
dictionary, checks endpoint labels/occurrences, charge, spin and graph/isotope
identity, loads only the source family, and preserves atom IDs through its
recipe. Family-canonical direction is explicitly translated back to the
requested reaction direction. Exact endpoint graph isomorphisms translate those
IDs into the existing SMILES atom order; graph-edit weights are not retuned.
Distinct recipe maps remain unresolved. The existing explicit-H layer must
agree with the trace's donor/acceptor/isotope accounting; otherwise no mapping
is admitted. Potential stereochemistry that the RMG artifacts cannot preserve
is not silently assigned. This handoff does not establish rate accuracy,
reaction completeness, or the physical correctness of a recipe mechanism.

Request/result sidecars are append-only and bind the mapping to source,
family and native implementation hashes. Geometry/reuse boundaries reopen and
revalidate them. Legacy routes without recipe provenance retain their original
serialized contract and do not acquire retroactive mapping evidence. Earlier
attempts may provide this structural trace without inheriting historical rates
or flux. Native bridge commands use process-group timeout/interruption cleanup;
the RMG supervisor likewise cleans its owned group before propagating an outer
interruption. Existing TS, frequency, IRC and rate-admission gates still apply.

The v8 ring-connectivity policy keeps the existing covalent-distance thresholds.
An ambiguous heavy-atom contact across a declared four-membered ring can be
resolved as nonbonded only for a closed-shell fragment whose independently
perceived coordinate graph exactly matches every atom-ordered bond and bond
order, formal charge and radical count. The RDKit connect-the-dots backend and
version are retained; bond-order assignment is capped at 1000 iterations.
Definite undeclared bonds, overlaps, missing bonds, interfragment near contacts
and failed/mismatched graph assignments remain fail-closed. Resolved diagonals
and the supporting graph check are explicit evidence, not new inferred bonds
or proof of a minimum, TS, rate, spin state or mechanism. Retained v7/v8 sessions
can enter the current v9 policy only through the append-only scientific-policy
amendment; their original policy contract is not rewritten. IRC branch recovery
remains enabled for those explicitly supported policies.
Serialized recipe maps reopen as exact integer pairs whether represented by
Python tuples or JSON lists; boolean/float indices and changed maps are refused.

The v9 endpoint-identity policy requires finite XYZ coordinates and an explicit
integer-only atom-order bijection. A declared tetrahedral or E/Z configuration
is checked using RDKit's 3D stereo assignment on a copy of the declared graph,
with backend/version and individual comparisons retained. This inference is
conditional on the declared connectivity and isotope labels; it is not an
independent bond-order, isotope, charge or spin assignment. Unspecified stereo
is not resolved. Unsupported declared stereo or failed assignment stays
unresolved. The geometry record explicitly disclaims stationary-point and
electronic-state applicability claims.

IRC trajectory integrity requires at least two frames, unchanged element order
and finite coordinates throughout. Both bound and separated ends must match
mapped declared connectivity and declared stereochemistry in one common
direction. Bound ends additionally need separate aligned RMSDs within the
retained 0.75 Å tolerance; a good end cannot compensate for a bad end. The old
summed RMSD remains diagnostic only. Separated channels retain their existing
exemption from arbitrary encounter-pose RMSD. Parsing/integrity failure withholds
endpoint admission; these checks do not replace ORCA's independent IRC
completion evidence.

The path calculation generates deterministic endpoint-complex orientations in
one atom order. For a one-fragment endpoint, the assembled endpoint must have a
converged ORCA optimization and a frequency-verified local minimum. For an
explicitly separated channel, every physical reactant and product species is
instead optimized and frequency-validated independently. Those validated
fragment geometries are assembled at a retained finite separation as NEB seeds;
the freely translating and rotating encounter geometry is not misrepresented
as a molecular minimum. The activated branch then runs ORCA NEB-TS, refines the
retained geometry with `OptTS`, requires exactly one significant imaginary
frequency, and runs a bidirectional IRC. Bound IRC endpoints are checked against
the optimized geometries and mapped declared chemistry; separated-channel IRC
endpoints are matched to declared fragment chemistry in either direction before the path status
becomes `verified_transition_state_path`. The independently validated species
also provide the thermochemistry Arkane requires.

Each streamed ORCA execution retains a process lease containing its PID,
process group, input, output, and terminal status. Verification continuation
uses a live lease or a fresh calculation-artifact modification time as evidence
that a detached engine may still be active and refuses a duplicate. A retained
`started` route becomes execution-recoverable only after the artifact family is
quiet and no live lease remains. This reuses the existing verifier-policy
attempt rather than manufacturing another scientific retry. If the interrupted
job is NEB/NEB-TS, STORCA validates the complete ORCA `.allxyz` band against the
declared image count and atom order, creates an append-only restart input using
ORCA's `Restart_ALLXYZFile`, hashes the checkpoint in restart evidence, and
continues downstream discovery using the restart job basename.
An ORCA/OpenMPI launcher denial is categorized separately from path chemistry.
The failed output and lease remain immutable; recovery creates an exact-input
`-execution-retry-###` calculation with a hash link to the denied output. This
reuses the same scientific-policy attempt because no electronic-structure
result was obtained. A normally completed calculation can never enter this
recovery branch.
Endpoint optimization/frequency and displaced-minimum workflows use the same
path executor through an internal job-runner adapter. Frequency input consumes
the geometry belonging to the executed retry basename, never a stale original
`opt.xyz`. Legacy endpoint startup failures without job contracts receive
hash-linked retries; a completed retry is reused only while the original input
and failed-output hash still match. Unknown incomplete endpoint output is not
overwritten. This correction cannot recover logs already overwritten by older
endpoint workflows; known incidents must remain explicit provenance gaps.
For contracted path jobs, a changed input or referenced geometry/Hessian now
fails before replacing the job contract or executing ORCA. Endpoint input
generation preserves existing files byte-for-byte when the generated input
matches and rejects conflicting input; an external seed cannot overwrite a
different retained endpoint seed. A different scientific calculation requires
a distinct workspace or the existing explicit amendment/retry mechanism.
IRC admission also requires both printed forward/backward direction tables to
meet their reported maximum and RMS gradient thresholds, normal termination,
and no iteration-limit marker. Missing direction evidence, an unconverged
gradient, or an iteration cap yields `irc_path_incomplete` before endpoint
connectivity can admit a TS path. This separates unfinished tracing from a
completed trace reaching the wrong channel. Geometry matching alone is not
evidence of IRC convergence. Older retained path claims have not automatically
been re-audited by this new check.
A bounded internal downhill-IRC continuation diagnostic can extend either
iteration-limited branch from its retained directional endpoint, following
ORCA's documented `Direction down` mode. It keeps the original method, charge,
spin, requested execution-resource contract, and iteration cap; checks endpoint/trajectory agreement; and
hash-links the source input, output, endpoint, and trajectory. Each original
branch has one fixed continuation workspace, and unknown unfinished output
blocks duplicate execution. Source TS displacement/Hessian directives are
removed because this starts from a nonstationary branch endpoint, not the TS.
This diagnostic does not itself admit endpoint identity, rates, or a lifetime.
A converged downhill segment alone is
not a bidirectional IRC claim. Reference:
https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/irc.html

The companion branch assembler can now join retained completed segments
without executing chemistry. It reparses convergence, checks original
input/output/endpoint/trajectory hashes, requires the exact derived downhill
input and seed, and checks both segment junctions and final coordinates.
Continuation assembly selects the actual executed input/output from the retained
continuation result, including append-only execution retries, rather than assuming
the `irc-downhill` basename. Each selected retry is tied to its exact immutable
job/reference digest and source-linked restart chain; input/output/trajectory
ownership, electronic method and seed are rechecked. Recognized infrastructure
retries may change only execution parallelism under the retained retry policy:
effective/requested or prior core counts must be positive integers matching the
actual inputs, including serial fallback and restoration. Changed, missing or
conflicting source hashes fail closed. A sibling execution's trajectory cannot
fill a missing selected trajectory. Legacy non-retry continuation results without
the newer execution snapshot are checked against their native job contracts and
current sources; they do not acquire invented historical evidence.
The backward extension is reversed, the original full IRC is retained, and the
forward extension is appended without duplicate junction frames. Missing or
unconverged branches stay unresolved. Joined trajectories and their manifests
are content-linked, and their only permitted claim is completed continuous IRC
tracing; endpoint identity and TS/rate admission remain separate gates. This
assembler is wired into `generic-orca-path-v7-irc-branch-recovery`: an eligible
capped IRC receives at most one extension per unfinished direction, still
subject to the surface compute budget. Only a successful assembly proceeds to
the ordinary chemical endpoint checks. Older policies do not silently gain
this behavior; retained v2-v6 sessions require the hashed policy-amendment
mechanism. The assembled trajectory remains a separate artifact rather than
overwriting an earlier trajectory cited by failed verification evidence.
`tests/replay_irc_branch_assembly.py` performs a read-only assembly replay against
an existing source-bound retained trace, intercepting only derived artifact
writers and refusing electronic dispatch or report overwrite. The actual NNONC
replay preserves 102 native files and their timestamps and reconstructs the same
119-frame trace byte for byte. Its real source uses the ordinary basename;
execution-retry, serial-fallback and mutation cases are intercepted workflow
regressions, not new physical calculations or full-chain resume validation.

When a completed trace has exactly one verified reactant endpoint but misses
the declared product channel, STORCA can retain a separate native reaction
candidate from the unexpected geometry. The initial proposal scope is a
neutral singlet surface with neutral, closed-shell fragment graph assignments;
other charge/spin inventories remain unresolved. The source trace hash,
endpoint/frame agreement, parent heavy-atom maps, and observed all-atom bond
changes are retained. Geometry-derived bond orders, fragment states, and
minimum identities are explicitly unvalidated. The parent reaction stays
unresolved, and its rate cannot be copied to the new candidate. Symmetry-related
hydrogen mapping is aligned to the retained endpoints only by permuting
same-isotope hydrogens attached to the same original reactant atom. Heavy-atom
maps remain fixed; isotope swaps, ambiguous contacts, and reassignments across
different donor atoms fail closed. The alignment retains canonical and revised
maps, endpoint hashes, and connectivity checks. This establishes correspondence,
not stereochemical equivalence, electronic-state validity, or path/rate admission.
Malformed, single-frame, inconsistent-atom-order, and nonfinite trajectories
cannot support a proposal. The engine now admits these proposals to a bounded
supplemental verification queue after checking native identity, unchanged
reactants/environment, source geometry hashes, and resolved trace correspondence.
This queue accepts one generation of unexpected chemistry, within its size bound
derived from the engine iteration budget. Duplicates are retained once, and
rejections are recorded alongside the parent path result. The queue is persisted
in engine state; candidates remain kinetically unranked and unverified. Selection
records distinguish IRC discovery from graph/xTB screening. This handoff does
not admit reuse of the retained TS or insert a rate into the network. Those
still require independent path/state/rate admission and repaired propagation.

Unexpected proposals now also retain a checked TS source bundle when available:
hashes of the optimization/frequency/IRC inputs and outputs, TS coordinates,
Hessian, and trajectory; input method and charge/multiplicity checks; normal
optimization/frequency completion; one meaningful imaginary mode at the explicit
-20 cm^-1 source-check threshold; frequency-input/TS atom and coordinate agreement;
TS inclusion in the trajectory; and rechecked original or continued IRC completion.
This bundle identifies reusable source calculations, not an admitted new path.
In particular, reaction-specific mode overlap, product minima and state
applicability, full condition/method-context admission, and rate theory must still
be validated before reusing it for the unexpected reaction. Legacy proposals
without this bundle remain candidates, not automatically reusable TS evidence.

The route-specific retained-path assessment rehashes the source files, reruns
source validation under the requested method/state, and binds the exact native
candidate to that trajectory. It reconstructs endpoint atom correspondence and
bond changes before testing the retained imaginary mode with the existing
bond-coordinate displacement diagnostic. This diagnostic tests whether at
least one changed bond moves appreciably; it does not prove synchronous bond
changes and its RMS-normalized magnitude is not a probability or cosine bounded
by one. A passing assessment is `retained_path_geometry_supported`, explicitly
not `verified_path`; stationary fragments, state applicability, rate theory and
new-route admission remain separate gates. Assessment artifacts live in a new
workspace, and changed retained source files or endpoint seeds are rejected.

The default ORCA adapter now passes queued unexpected-channel proposals to the
ordinary generic verifier's retained-path branch. Failed source assessment does
not silently launch a replacement NEB. Passing assessment supplies aligned
fragment seeds to the existing minimum/frequency checks, while the ordinary
environment, surface-state, and thermochemical loss-bound gates remain in force.
Only after those checks may the branch emit the usual conditional
`verified_transition_state_path` with Arkane stationary-species inputs. This
means a validated connection on the declared electronic surface/method, not
proof of electronic-structure accuracy or exhaustive product/conformer coverage.
Rate applicability, rate uncertainty and repaired-network propagation are still
required. The source check also compares Hessian coordinates (Bohr converted to
angstrom) with the retained TS, not merely matching atom labels.
An origin translation is allowed and recorded; internal coordinates and the
Cartesian axes must agree. Rotations are not accepted by this source policy
because mode vectors would also require a validated rotation. The real saved
Hessian differed by translation only, with approximately 1.3e-8 angstrom maximum
residual after centering.

Completed external stationary folders may be reused internally for a rechecked
unexpected path, without copying jobs or executing ORCA. Reuse requires all
seed/contract/metadata/optimization/frequency artifacts; validates seed content
against the retained endpoints and seed hash; checks method/state/resource
metadata, normal optimization/frequency completion and minimum/identity gates;
and retains source hashes in the new species records. Missing or changed
evidence is an error, not permission to launch replacement calculations.

If the bounded launcher policy reaches a serial fallback, a later continuation
does not silently make that fallback permanent. Once the serial PID and process
group are confirmed absent, STORCA compares its effective core count with the
immutable surface context. A mismatch creates another append-only execution
input that restores the retained parallel resource contract, links the serial
input and context hashes, and records that the scientific model did not change.
A live or indeterminate process still blocks this recovery.

The route, mapping, environment applicability, electronic-structure model, and
versioned surface policy are retained in an immutable scientific context.
Policy upgrades create a hash-linked context record under
`surface-resolution-contexts/`; the predecessor is never rewritten. Execution
core changes may be recorded in that amendment, while route, mapping,
environment, charge/multiplicity, and method changes are rejected. The retired
v4 tuple/list context-comparison failure has one narrowly classified,
pre-calculation same-policy recovery and cannot form an unlimited retry loop.
Verification-engine state files are immutable even when a detached or failed
session is later continued. STORCA therefore writes an append-only session
reconciliation ledger: any historical `running` record with a hash-linked
descendant is reported as `superseded_by_continuation`, while only the leaf
session is assessed as active, interrupted, or terminal. This operational
interpretation never changes scientific evidence or a lifetime claim.

The no-saddle branch is deliberately narrower than a visual downhill curve. A
barrierless classification requires normal ORCA completion, the expected path
coverage, a complete finite energy profile, endpoint agreement, and compatible
profiles from at least two independently oriented endpoint complexes. A partial
or single-orientation curve remains `surface_unresolved`. Even a reproducible
`barrierless_capture_candidate` or `barrierless_dissociation_candidate` is path
evidence rather than a rate: the current default Arkane TST adapter rejects it
until a collision-bounded capture and recrossing model is supplied.

An unresolved NEB now enters the versioned
`generic-orca-path-v6-dissociation-loss-bound` policy rather than an unlimited retry
loop. Retained energies, structural arc length, mapped bond changes, fragment
motion, per-image covalent-connectivity changes, warning history, endpoint
agreement, and independent orientations form a diagnosis record; that record
allocates work but never proves a mechanism. Fragment motion includes relative
translation and orientation where a separated-side atom partition can be
mapped into the reactant atom order.
Normally terminated but unconverged stationary optimizations retain their last
geometry and use append-only redundant/Cartesian, Hessian, and trust-radius
recoveries. Meaningful imaginary modes are displaced in both directions rather
than being dismissed by a shallow-frequency label.

When independent evidence indicates poor path seeding, the resolver derives at
most two relaxed coordinates from the retained broken, formed, or bond-order
changes. Incomplete scans, endpoint maxima, and nonunique maxima are rejected.
An accepted interior maximum is only an OptTS seed and must pass the same
one-meaningful-mode, intended-coordinate, and bidirectional IRC gates as a NEB
seed. A reproducible no-saddle segment is routed to a barrierless specialist
boundary and is never sent to harmonic TST.

A frequency-validated intermediate is converted from XYZ to a native species
only when connectivity is determined first and gives one molecular fragment;
bond orders are assigned only after that gate. The resulting structure must
give one distinct molecular graph
with the retained atom inventory, total charge, and a radical inventory
compatible with the selected multiplicity. Common atom-map numbers link that
graph to the parent route. STORCA then creates reactants-to-intermediate and
intermediate-to-products native candidates and verifies each bound segment
independently. Separated segments retain capture/dissociation requirements.
When both ordinary segments pass, Arkane generates a separate condition-bound
rate/library for each; all libraries and the intermediate are supplied to one
full repaired-model rerun. No effective parent barrier is constructed. Failed
network acceptance is retained as
`verified_multistep_path_network_insertion_pending` and withholds the lifetime.
If an ORCA intermediate warning optimizes back to a declared reactant or product
graph, it is retained as an `endpoint_return`, not promoted to a species and not
reported as an unresolved identity. Only that already validated image warning
may be ignored during one checkpoint continuation; any warning at a new image
stops the band for independent identity validation.

A second bounded allocation diagnostic handles endothermic dissociation paths.
If a NEB highest-energy image remains an endpoint for 24 trailing iterations
and repeatedly crosses below ORCA's climbing-image switch force, STORCA stops
that band as unresolved and reallocates work; this is not a no-barrier claim.
If an interior climbing image forms but remains below the higher fixed endpoint
asymptote for 12 consecutive climbing-image iterations, ordinary forward
harmonic TST is marked inapplicable. The route returns
`dissociation_rate_model_required` with explicit VTST, capture-theory, and
RRKM/master-equation evidence requirements. A submerged saddle may be real; it
simply cannot be the forward energetic bottleneck toward the higher endpoint.

Before spending a separated-channel NEB, the v6 policy also applies a
condition-bound detailed-balance upper limit to neutral gas-phase 1-to-2
dissociation when every fragment has validated ORCA thermochemistry. The
reverse association is bounded by an explicitly inflated hard-sphere
collision ceiling, and the ideal-gas thermochemistry pressure standard state
is retained in the evidence. If even this forward dissociation upper bound
stays below the requested loss threshold, path and rate resolution are not
required for that route. If the bound crosses the threshold, it supplies no
rate estimate and the route still escalates to capture/VTST or
RRKM/master-equation treatment.

Specialist rate solvers now enter through a source-neutral, condition-hashed
provider contract. The request binds route identity, path classification,
allowed rate theories, condition contract, and path-evidence hash. A result is
admitted only when it retains provider name/version, explicit method-
applicability diagnostics, calibrated multiplicative uncertainty, hashed raw
artifacts, and a condition-matched generated kinetics library whose admitted
files have explicit content checksums. RRKM/master-
equation results additionally require a pressure grid containing the requested
condition, bath-gas composition, and an energy-transfer model; VTST results
require a reaction coordinate, sampled dividing surfaces, and a stated
variational criterion. Capture results additionally require a sourced
long-range interaction model, sampled variational surfaces, explicit
recrossing evidence, symmetry and electronic-degeneracy factors, and an
applicable electronic-state assessment. Accepted specialist rates enter
the same library validation, full RMG repair, exact-rate matching,
repropagation, sensitivity, and convergence loop as Arkane rates. No default
specialist executable is silently selected. A deferred route can export its
exact request with `storca stability-export-specialist-request <run>` and can
consume a matching result with `storca stability-resume-verification <run>
--specialist-rate-result <result-or-directory>`. That continuation reuses the
hashed retained path classification and does not rerun ORCA. Relative paths in
a submitted result bundle resolve against the result file, and the submitted
JSON hash is retained as import provenance. This is an evidence-import adapter,
not a claim that STORCA itself performed or validated the numerical solver's
scientific calibration.

The specialist-rate benchmark evaluates provider conformance separately from
numerical accuracy. The repository development set contains evaluated low- and
high-pressure association-rate anchors and an electronic-state challenge, but
no locked holdout. Consequently it cannot authorize a supported capture lane.
The retained JPL 2011 source audit visually checks the relevant Table 2-1 rows
and their bath, unit and uncertainty conventions. Its four reactions supply
paired low-/high-pressure limiting recommendations, not eight independent
experiments. The published multiplicative factors are anchored at 298 K:
`f(T) = f(298) * exp(g * abs(1/T - 1/298))`. The eight 300 K reference factors
therefore use `f(300)`, while central rates, conditions, partitions and scoring
thresholds remain unchanged. The original manifest bytes and source PDF hash
are retained. These evaluator-assigned envelopes specify neither an exact
confidence probability nor independent parameter distributions; they do not
calibrate STORCA's electronic structure, capture or master-equation errors.
See `benchmarks/kinetics/reference_sources/2026-10-03-jpl2011/source-audit.json`.
The backend decision record in
`docs/stability/SPECIALIST_RATE_BACKENDS.md` identifies VaReCoF-class VRC-TST
as the capture producer target and MESS-class master-equation calculation as a
separate pressure-dependent propagation target; neither is silently selected.
The MESS bridge validates and renders a narrow
one-well/two-fragment model with explicit RRHO evidence, collision and energy-
relaxation sources, and an external Rotd microcanonical-flux artifact. Its
native temperature/pressure-table parser requires an exact named channel and
complete finite grid. The bounded execution layer revalidates the prepared
hashes, fingerprints the executable binary, retains timeout/nonzero-exit
evidence, and reuses only an immutable matching checkpoint. Native MESS units
are retained without conversion at execution. A separate admission gate then
requires a binary-bound, source-backed interpretation of the otherwise
ambiguous bimolecular `cm^3/sec` amount basis; no unit basis is inferred. It
also requires an exact-grid condition, a complete rate-uncertainty budget,
full-grid collision ceilings for bimolecular association, and independently
sourced equilibrium rate ratios. The opposite MESS channel is parsed and
forward/reverse detailed balance is checked at every temperature and pressure.
Only an admitted table may be fitted to an RMG Chebyshev or PDepArrhenius model.
The serialized library is reloaded by RMG and reevaluated over the complete
grid under a predeclared log-rate tolerance. The final result uses the existing
specialist-rate contract and therefore enters the ordinary repaired-network,
repropagation, sensitivity, and convergence loop. Execution alone still
explicitly prohibits a verified-rate, falloff-accuracy, or lifetime claim;
admission still does not validate MESS against experiment or permit a lifetime
before network repropagation.

Arkane's verified TST or pressure-dependent output is a verified *central rate
model*, not an exact rate. Every generated Arkane manifest therefore records
`rate_uncertainty.status = calibration_required`. When such a replacement
carries nonzero target-loss flux, the verification engine requires an explicit,
condition-matched uncertainty record and simultaneously propagates its lower
and upper multipliers with any still-unresolved rate bounds. A lifetime or
retention claim is withheld if the uncertainty is missing, the bound
propagations fail, the threshold crossing changes, or the t95 span exceeds the
declared convergence tolerance. STORCA supplies no universal DFT/TST or
master-equation error factor.

For the existing Arkane 1↔2 central-rate lane, the generated manifest exposes
the modified-strong-collision method, bath-gas mixture, grain controls, and
whether collisional inputs were estimated or explicit. In the default branch,
RMG supplies generic `SingleExponentialDown` and transport estimates and their
applicability is explicitly *not* claimed. A route may instead carry a
schema-v1 `storca_explicit_gas_pressure_dependence_model`: it must provide an
applicable temperature/pressure range, a sourced `SingleExponentialDown`
`alpha0`, `T0`, and exponent for the exact energized well, sourced
Lennard-Jones `sigma`, `epsilon`, and shape index for the well and every bath
gas, and hashed evidence artifacts. STORCA writes those values directly into
the Arkane species definitions and retains the normalized model hash. Even an
explicit central model remains uncertain; terminal prediction calibration must
cover either `total_rate_model` or every required component.

One-sided physical ceilings retained for unresolved paths are converted to
rate-uncertainty records only when their units match a positive condition-rate
baseline. The lower multiplier remains zero rather than being assigned a
convenient finite value. After at least one verified network repair, if every
remaining nonzero-flux route has explicit condition-matched uncertainty and a
retained baseline rate, the engine automatically performs combined lower/upper
bound propagations. Observable convergence can then be claimed only if the
full requested classification or t95 remains within the declared tolerance
across those simultaneous bounds. Local one-at-a-time sensitivity remains
allocation evidence and is never substituted for this global robustness test.

Every policy attempt has an immutable budget contract. New ORCA executions are
counted from process leases, elapsed core-hours and wall time are recomputed
before launch, and the next job timeout is reduced to the remaining allowance.
Each decision writes an append-only budget snapshot. A changed policy creates a
new budget contract rather than mutating the previous scientific attempt.
Standard NEB orientations, stationary recoveries, guided scans, and individual
segments also write the same `storca_surface_resolution_attempt` evidence
envelope. It links endpoint/input hashes, mapping, charge/multiplicity, method,
environment applicability, policy, parameters, route-budget consumption,
process leases, detailed artifacts, outcome/failure, and allowed/prohibited
claims while leaving the detailed calculation records intact.

Before any generic path job is launched, the declared phase is checked against
the surface policy. Gas paths reject accidental implicit-solvent keywords.
Implicit-solution paths require a declared solvent, CPCM or SMD model, and one
exact matching ORCA solvent keyword; this is normally supplied by the dedicated
solution-path lane. Unsupported periodic, explicit-solvent, interfacial, solid,
electrochemical, or other environments stop as
`environment_model_escalation_required` rather than silently using the gas or
implicit-solution machinery.

For a verified transition-state path, Arkane reads the validated separated
species and TS frequency artifacts and evaluates a grid containing the exact
stage condition. Pressure-dependent 1↔2 calculations require the declared,
supported non-target bath gas. High-pressure-limit TST does not require or
emit a bath-gas/collisional model, including for a pure-target gas; this does
not establish that the high-pressure approximation is physically applicable.
The
condition rate must have the expected molecularity and units and, for a
bimolecular route, must not exceed the retained collision ceiling. Atom-balanced
routes use explicit artifact validation; missing stationary-point, Hessian,
state, or rate evidence fails closed.

### Native ORCA rotational-mode admission

Verified gas routes launch Arkane through a STORCA-owned, process-local ORCA
reader in the environment selected by `--rmg-env`. No installed RMG/Arkane
package is edited, and no new user flag is needed. The policy
`storca-orca-diatomic-rotor-v1` addresses the exact-zero inertia classification
defect reported in [RMG-Py issue #3016](https://github.com/ReactionMechanismGenerator/RMG-Py/issues/3016).
For two distinct atoms with finite coordinates and positive native atomic
masses, the transverse inertia is computed analytically as
`I = m1*m2/(m1+m2) * |r2-r1|^2`. The native external rotor is replaced with
`LinearRotor`, preserving its symmetry and quantum treatment. Translation,
electronic energy, unscaled frequencies, spin and optical-isomer factors are
unchanged; Arkane still applies its ordinary ZPE, scaling and tunneling model.

Atoms retain no external rotor. Clearly nonlinear molecules retain the native
rotor. Already-correct native linear polyatomics are admitted only when the
geometry is collinear within the policy's numerical boundary and independent
principal moments confirm the linear rotor; they are not corrected. Ambiguous
near-linear polyatomics, invalid mode inventories and inconsistent source
evidence fail closed. This is an external rigid-rotor ingestion policy, not
an anharmonic, hindered-rotor or rovibronic theory.

Each job admits only the hashed stationary-point logs and colocated Hessians
in its contract. Native symmetry uses owned temporary scratch. The older
RMG 3.3 ORCA Hessian reader is called within its admitted artifact directory
to resolve its basename, with hashes rechecked afterward; the process working
directory is restored. The successful receipt binds input, provenance, output,
adapter/native source hashes and every stationary point's before/after audit.
Full native inertia and its actual serialized representation are both retained;
Arkane's printed rounding is checked against that representation rather than
mistaken for a change to the full-precision model.

An Arkane replacement without a valid route-matched receipt cannot enter a
network, even if its rate is below the collision ceiling. Its independently
retained `arkane_rate_source_artifacts` input/output/provenance records must
exactly match that receipt, and supplied conditional forward evidence must
agree with those sources and the route. A receipt from another calculation of
the same route is insufficient. Legacy artifacts
remain readable as diagnostics, but require a fresh audited rate calculation
for replacement admission. Replays reuse validated ORCA/TS/IRC artifacts in
a new workspace; retained calculation evidence is never overwritten. The
receipt establishes correct loading under this policy, not experimental rate
accuracy, physical uncertainty calibration or a lifetime.

## Kinetic reinjection and ladder stopping

Each accepted rate becomes a minimal route-specific RMG kinetics library with
its provenance and original reaction signature. All accumulated verified
libraries are then supplied to a full RMG model-generation and propagation
rerun—not an algebraic rescaling of the previous target profile. The rerun must
retain the exact condition contract, load each replacement at the verified
rate, pass mechanism and collision checks, complete the requested propagation,
and provide new flux attribution. The verifier reranks that repaired network
and repeats upstream/path/rate/reinjection steps as needed.

Flux coverage alone is not a terminal criterion. A condition-specific `t95` is
released only as `verified_t95_converged` when repaired propagation crosses the
retention threshold and the result is robust to the remaining unresolved-route
bounds. A robust repaired model that stays below the loss threshold may report
`verified_below_loss_threshold`, but the ladder continues to the next condition
because no `t95` was found. Likewise, `no_target_loss_in_completed_rmg_model`
describes only that completed retained model and does not stop the ladder as a
universal stability conclusion. Any path, rate, reinjection, coverage, or
robustness failure produces `verification_incomplete` and no final lifetime.

When a retained bimolecular RMG rate exceeds its independently evaluated
collision ceiling, STORCA does not silently clip the mechanism or call that
rate verified. For a completed below-threshold propagation, it may instead run
the exact retained network at two joint rate-bound endpoints: removal of every
indexed violating reversible reaction and common forward/reverse scaling to
each collision ceiling. Numerical flux-closure error is added to each computed
loss before comparison with the immutable threshold. A terminal result is
labelled `collision_bounded_for_observable_not_rate_verified`, retains both
propagations and every multiplier, and remains subordinate to independent
discovery and applicability gates. This endpoint envelope is a conditional
retained-network robustness test; it is not an exact rate, a completeness
claim, or proof of behavior outside the declared bounds.
# Rate execution safety checkpoint (2026-09-10)

Arkane execution and RMG/Arkane syntax/library subprocesses have a finite
600-second per-command safety ceiling. POSIX timeouts terminate their process
group and retain stdout/stderr; they do not trigger automatic retries or rate
admission. This is not aggregate campaign-budget accounting. The neutral-gas
conditional thermochemical screening bound recognizes the declared low-pressure
gas surrogate as well as the ordinary gas surrogate, retaining the actual
pressure. This spelling correction supplies neither a falloff model nor proof
that the assumed collision ceiling is valid for a particular reaction.

## Arkane electronic-energy reference

For generated verified-route inputs, each electronic energy is shifted by its
atom count times the mean electronic energy per atom of the complete reactant
channel. Validated atom balance makes this a common channel reference, preserving
reaction energies and forward/reverse barriers even when molecularity changes.
The values are tagged `e_electronic` so Arkane still adds scaled ZPE. ORCA logs
and Hessians are unchanged; the gauge is recorded per species in provenance.
These are not absolute formation energies and must not be merged across route
gauges as species thermochemistry. The change avoids overflow from raw absolute
energies in multiple RMG pressure-network Boltzmann expressions without changing
the installed RMG code. A live pressure-dependent calculation completed after
this correction; that alone does not validate its physical applicability.

If Arkane omits its native species dictionary, STORCA serializes the declared
reactant/product graphs in the configured RMG environment. Multiplicity,
duplicate-label consistency and adjacency roundtrips are checked, and the
derived dictionary has separate provenance. This performs no thermo estimation
or kinetics calculation. The assembled library must still reload in RMG and
reproduce the retained condition rate before admission. A calculated central
rate can pass that check while remaining uncalibrated for lifetime prediction.

### Isolated source-separated composite-rate recipes

`storca.composite_rate_recipe` provides an explicit gas-phase diagnostic API
for separate electronic-energy and geometry/vibration evidence. Each source
binds native input/output hashes; Hessians, exact nuclear order/coordinates,
charge, multiplicity, stationary-mode inventory and method declarations must
agree. Native Hessian COM-origin translation is ignored, not arbitrary
rotation/permutation. It does not launch calculations, edit ORCA logs, infer
scales/symmetry or grant route/rate/network verification.

Fresh recipes use software-qualified native Arkane method objects, unchanged
frequency/Hessian files, electronic-only scalars and the balanced route gauge
described above. Explicit state/counting/soft-mode assumptions and derivative/
route evidence are retained. The existing rotor contract is reused; these
diagnostic artifacts are not admitted as repaired network libraries.
The OH/methane development study separates raw harmonic ZPE from native
Arkane's `frequencyScaleFactor/1.014` ZPE convention, bare TST/Wigner/Eckart,
direct evaluations and fitted rates. Its five matched CCSD(T)/QZ single points
reduce the direct Eckart discrepancy, but do not certify physical accuracy or
promote production defaults. See
[the completed audit](docs/stability/OH_METHANE_CONSISTENT_RATE_AUDIT_2026-10-11.md).

## Verified native reactions in repaired RMG networks

Checksum-validated generated libraries marked as kinetic repairs are supplied
as native RMG core seeds. This ensures a newly discovered, verified channel is
present in the propagated core rather than merely loaded into the edge.
Ordinary reference reaction libraries retain their existing edge role. Input
preparation deduplicates paths and rejects conflicting registered names.

The STORCA launcher resolves external seed paths through RMG's own loader
registry. Because the RMG 4.0.0 restart loader clears that registry while
retaining loaded libraries, the launcher keeps process-local aliases bound to
the actual library objects. It rejects unloaded or overwritten bindings; it
does not guess a basename or change the installed RMG package. This was checked
with actual saved restart libraries as well as synthetic regressions.

Core admission does not establish valid propagation or a rate's physical
accuracy. Exact post-propagation rate matching, network closure, applicability,
collision diagnostics, and uncertainty gates still apply. The retained NNONC
repair demonstrates initial core insertion, not a converged stability result.

The launcher also preserves the in-memory gas Chebyshev model's temperature
and pressure bounds and polynomial coefficients at binary64 text precision
when RMG writes Chemkin. Its ordinary reaction/comment/duplicate handling is
retained. Bounds must be positive, finite and ordered; coefficient dimensions
must agree. The logarithmic SI-to-cm/mol/s conversion follows molecularity
without mutating the source model. This fixes RMG 4.0.0's three-decimal pressure
format rounding a 1e-6 bar minimum to zero and its low coefficient precision;
it does not refit rates or improve source electronic-structure accuracy.

The bridge rejects nonfinite/negative forward rates with reaction-specific
diagnostics and emits strict JSON. Invalid reverse rates or collision ceilings
remain unavailable, not inferred. A derived export using the retained NNONC
core and the uniquely forward-isomorphic, checksummed source-rate model
reproduced the original source rate at all 30 tested T/P points (zero observed
relative discrepancy). This serialization diagnostic ran no propagation or
new chemistry and did not alter source files. Source-rate precision, network
truncation, other invalid estimates and uncalibrated uncertainty remain open.

Replacement identity is now checked in the owning RMG environment by
`rmg_identity_bridge.py`, independently of species labels. Source and model
species are compared across resonance forms using strict graph matching with
unchanged charge, multiplicity, isotope composition, and occurrence-level
stoichiometry. The source must identify exactly one structural channel in the
inspected model with exact source-library provenance. Reverse correspondence
is permitted only for a reversible source and model; it selects the reverse
rate and molecularity for the subsequent numerical comparison. Parallel
structural channels are not disambiguated by choosing the closest rate.

`reaction_correspondence.py` retains a content-addressed identity sidecar bound
to the exact source library, source manifest, mechanism files, and inspection
records. Both applied-rate verification and verification planning consume this
same evidence; native route identifiers and original source records remain
unchanged. Changed files, ambiguous mappings, failed loading, or unresolved
stereochemistry cannot fall back to label matching. Historical inspections
without required structural evidence retain an explicitly labeled legacy
label-matching path; that compatibility path is not a structural identity proof.

RMG graph artifacts do not establish stereochemical identity. Any supported
resonance form exposing potential stereochemistry withholds admission and
retains only diagnostic structural candidates. This is an evidence limitation,
not a claim that each flagged configuration is physically long-lived. A source
library loading error is separately reported as
`replacement_source_library_load_failed`, with the original exception retained;
it is not a finding that no corresponding reaction exists.

Generated photon-source entries explicitly declare `reversible=False` as well
as an irreversible arrow. Native RMG requires these declarations to agree.
The writer's roundtrip regression verifies loading, atom balance, irreversibility,
and preservation of the supplied forward photolysis rate. A checksum-valid
manifest alone does not prove that RMG can load a library, and successful loading
does not validate the supplied spectrum, quantum yield, or product branching.

## Prediction-first checkpoint and scientific-decision boundaries

The exact NEB time-limit continuation policy admits only a hash-bound native
checkpoint from a terminal owned process. It requires the retained timeout,
finite complete atom-ordered images, unchanged job dependencies and confirmed
absent PID/process group. A source segment can authorize only one continuation;
the retained aggregate budget still applies. Restart evidence hashes the source
input, output, checkpoint, contract and lease without rewriting the failed job.
This is execution recovery, not a saddle, path-completeness or rate claim.

New ordinary Arkane callbacks retain a versioned method-applicability decision
before execution and attach its hashes to returned rate/library evidence.
Explicit environmental, state or path contradictions block admission; missing
diagnostics remain named conditional assumptions or escalation requirements.
This policy does not automatically solve multireference, tunneling, solvent
participation or conformer-sampling questions and does not calibrate rate errors.

The bounded scientific panel exports reference-free prediction requests and
keeps missing, crashed, unresolved and invalid outcomes in its ledger. Invalid
predictions cannot contribute admitted accuracy scores. References and prospective
holdouts need independent curation before formal evaluation. Earlier-only RMG
attempt channels enter the existing supplemental queue as source-revalidated
structural evidence, without inherited historical rate, flux or convergence.

The selected solution graph/xTB candidate adapter now reopens native endpoint
contracts/results, all retained input/output/optimized-geometry/JSON hashes,
stoichiometric endpoint identities, observed charge/unpaired-electron settings,
method/invocation/solvent and energy sums before invoking the existing path
producer. The reopened sources are bound into its immutable request and checked
again after production. This is source-consistent screening provenance, not
spin purity, ground-state ordering, a saddle or rate-accuracy evidence. Actual
three-case preparation and corrupted-derived-source replays preserve the original
77 files/13 directories with zero new chemistry. Independently sourced paths
retain their own discovery contracts; they are not required to pretend to be xTB.

One descriptive JPL H + O2 high-pressure-limiting-rate comparator is now bound
to its exact source table, footnotes, units, conditions and evaluated uncertainty
metadata. Its reference is theory-informed and does not supply a measured
finite-pressure lifetime, DFT error calibration or a prospective holdout. The
offline outcome is explicitly missing. Calculation exports preserve the limiting
pressure regime without leaking the scoring reference or substituting a finite
pressure; comparison scope is checked across methods.

The retained ORCA 6.1 `Restart_ALLXYZFile`/`PreOpt true` parser conflict has a
separate narrow repair. It requires the exact pre-electronic parser failure,
immutable job dependencies and complete ordered source band, and a terminal
absent owned process. Only `PreOpt true` is changed to `false` in an appended,
hash-linked input. Ineligible conflicts cannot fall through to ordinary
same-basename execution. This is syntax recovery, not chemical verification.

NEB saddle/trajectory and IRC trajectory selection is exact-execution scoped:
no prefix or directory-wide fallback may borrow an earlier or descendant job's
artifact. The retained timeout replay verifies this boundary with unchanged
sources. A generic exact-basename XYZ, where accepted, remains only a TS seed;
optimization, frequency and endpoint-matched completed IRC checks are still
required. The controlled-interruption request additionally requires a completed
NEB iteration and retains its output hash; a band file alone is insufficient.

ORCA-flagged intermediate images are also execution-scoped. Validation selects
only `<actual-output-stem>_imN.xyz` and retains output/image hashes and the
requested index. Its cache is keyed by that exact execution and full source
binding, not a shared `intermediates/image-NN` directory. Changed sources during
validation withhold identity and endpoint-return dismissal; legacy caches are
preserved. The real replay demonstrates a corrected dispatch, not a newly
optimized intermediate or a changed lifetime.

Requested-temperature thermochemistry does not adapt named CPCM/SMD solvent
parameters or validate liquid-phase persistence. Those fixed-parameter,
single-selected-TS and harmonic-frequency assumptions remain explicit. Bounded
reactant refinement is now integrated, but ensemble interconversion and
TS-conformer coverage are not established. The temperature
producer is software-tested, including rate arithmetic and propagation with
synthetic verified-path evidence; a real condition-matched solution benchmark
and automatic applicable uncertainty remain required.

### Reached IRC minima and conformer-specific rate handoff

The `generic-orca-path-v10-reached-endpoint-minima` policy separates mapped
chemical identity, reference-conformer RMSD and stationary-point validation.
The reference RMSD tolerance is not loosened. When both mapped connectivity/
stereo checks pass but bound-reference RMSD fails, exact independent coordinate
bond/order/formal-charge/radical assignment can corroborate ordinary neutral,
non-zwitterionic closed-shell endpoints. This bounded RDKit diagnostic does
not establish electronic-state applicability; unsupported or ambiguous
assignment remains `irc_endpoint_identity_unresolved`, not a new-product claim.

Validated IRC ends are oriented, extracted into declared fragment atom order,
locally optimized and frequency-checked at the route method/state/environment.
Both ordinary and historical paths now use these reached stationary records
for Arkane, rather than regenerating a SMILES seed or borrowing the historical
reference minimum. Hash-bound optimized XYZ, optimization output, frequency
output and Hessian artifacts are reopened before and after Arkane handoff.

A different or unresolved refined-reference conformer basis allows retention
of a conditional elementary rate, but forbids graph-species network reinjection
until thermochemical/equilibration evidence resolves the mismatch. No fast
conformer equilibration or global-minimum/ensemble completeness is inferred.
Geometry correspondence alone is not calibration of RMG thermochemistry.

#### Conditional forward evidence and bounded reverse omission

Successful path/minimum and forward-rate evidence is retained independently of
network admission. `conformer_thermo_basis_unresolved` is a scientific handoff
refusal, not an Arkane execution failure. New Arkane manifests retain native
condition-evaluated forward and reached-endpoint reverse fit evidence. Legacy
evidence can be reopened read-only, with exact statmech/geometry/frequency/
Hessian/energy-gauge bindings; missing reverse evidence never erases a valid
conditional forward rate. A rounded printed reverse-rate table is diagnostic,
not a physical upper bound.

The `finite-horizon-reverse-return-v1` policy can bound omission error **relative
to named fixed kinetic models**. It requires a closed, isothermal, conservative
first-order network, pure parent initial inventory, a one-parent/one-product
route and matching reactant-reference conformer evidence. Native species Kc
is evaluated directly from RMG thermochemistry at the requested temperature;
it is not inferred from independently fitted reverse and forward coefficients.
The envelope takes the maximum of the reached-endpoint reverse fit and the
new forward fit divided by the oriented native species Kc. Those central rates
are fixed-model alternatives, not bounds on unknown physical rates.

For product inventory ceiling B=1 relative to initial parent, the absolute
parent-retention perturbation is conservatively bounded by t * k_return_upper
(summed over distinct omitted returns, capped at one). Positivity and
nonexpansiveness of the conservative linear network justify the bound; it does
not apply to bimolecular, amplifying, open, adiabatic or transport-coupled
networks. The default absolute tolerance is 1e-4 and cannot exceed the allowed
parent loss. The exact condition contract, horizon, model inventory, rate
alternatives and reopened source hashes are part of the recomputable evidence.

Only a separately named, checksummed `=>` / `reversible=False` library may use
this approximation. The original reversible library and unresolved conformer
basis remain unchanged. Ordinary single-rate engine handoff tries this path
only when the requisite new evidence exists; legacy/unsupported inputs retain
their scientific blocker. Generated-library and identity admission require
the full condition contract, not temperature/pressure alone. The repaired
model must be rechecked: new species, channels, changed spectator rates or a
reintroduced reverse process invalidate this particular certificate.

An accepted result is `conditional_model_prediction`, with parent profile and
finite-horizon omission error stored separately. It earns no verified flux,
conclusive lifetime, t95, rate calibration, high-pressure applicability,
conformer-completeness or network-completeness credit. A small omission error
does not make the assumed electronic-structure/rate model accurate.
The omission bound also excludes numerical propagation/presentation error and
rate-export roundoff. Repaired inventory matching uses a 1e-12 relative rate
serialization tolerance; the ordinary historical rate-matching gate is not
relaxed. Neither solver tolerances nor a displayed curve establish a bound on
total prediction error.

An expanded closed first-order retained mechanism may now receive a separate
`source_bound_retained_linear_network_v2` assessment, without rewriting or
widening its older library certificate. Full native identifiers and state
evidence are retained; parent/product correspondence requires the existing
modern native proof. Its rates are Cantera operative constants, not RMG's
fitted reverse inspection values. Native conversion differences from the
retained forward fit remain explicit under the existing insertion tolerance.
Only a virtual direct product return is added to the comparison inventory;
actual intermediate recycling stays explicit. All conditions, model exports,
identity/library ancestry and implementation hashes are source-bound.

Declared branch scenarios are diagnostics, not `RateUncertaintyRecord` bounds
or terminal robustness envelopes. Directed return/escape removal does not
preserve the original equilibrium ratio. The scoped native solver uses stiff
BDF with tolerance refinement and Radau, integrated directed flux and
tiny-population resolution checks, plus two tightened exact-horizon Cantera
baseline solves. It rejects nonlinear/collider laws, conflicting environment
declarations, conservation/positivity failures and unresolved numerical
agreement; it never clips, renormalizes or silently applies QSS. Pre-path
endpoint-state audits flag spin-changing channels and withhold ordinary-path
eligibility. Saved responses guide explicit next-evidence requests without
inventing physical uncertainty, cost or information-gain scores. See the
[actual implementation checkpoint](docs/stability/EXPANDED_NETWORK_ASSESSMENT_2026-10-08.md).

Explicit completed-IRC reuse reopens terminal owned jobs, TS/Hessian/mode and
scientific context, reconstructs continued branches and their assembled trace
in memory, and checks pinned hashes. Only tracing evidence is reused. New
endpoint/minimum/rate checks run in a separate workspace; no NEB, OptTS, TS
frequency or IRC is reexecuted. Source contracts are reopened after endpoint
refinement and before rate inputs are exposed. Old artifacts retain their
original policy and conclusions; reuse does not renew a compute budget.

### Completed TS-stage reuse after an incomplete IRC

An explicit source-bound derived historical trial may reopen completed OptTS
and frequency calculations without rerunning them. It checks the declared
scientific-context pointer/hash, exact method/state/map/conditions, job-contract
owners and dependencies, normal converged stationary observations, ordered
XYZ/inline-frequency/Hessian geometry and spectrum, one meaningful imaginary
mode with intended-coordinate overlap, and terminal absent source processes.
Only geometry/Hessian are copied into a new IRC workspace; original outputs
and leases retain their ownership. The stopped forward-only IRC from native008
is not eligible for the existing normally terminated downhill continuation
policy and cannot supply connection evidence. A fresh bidirectional IRC must
pass the unchanged completion and mapped endpoint checks, with source
re-admission before rate handoff. This does not renew compute budgets or prove
full-chain interruption/resume equivalence, experimental accuracy or lifetime.

### Primary-rate comparisons and counting-policy audits

The separate Berg2024 OH + methane reference slice freezes measured supplement
rows and selects its descriptive scalar by original table order, without using
prediction agreement. Total bimolecular OH-loss units, nominal reaction
pressure and exact measured temperature remain source bound. Point-error sigma,
covariance and total-uncertainty unknowns are not assigned confidence levels;
one study's temperature rows are not independent calibration reactions. The
bounded evaluator does not dispatch chemistry, leak rates/products to request
inputs, validate HPL applicability or promote a missing prediction to support.

Independent valid compute evidence is revalidated even when a prediction
envelope is malformed. Failed chemistry remains invalid/unscored; corrupt or
missing cost proof remains unknown rather than zero. Recorded clock/measurement
bases are preserved and finite resource metadata is not scientific evidence.

The installed-Arkane optical diagnostic compares the producer's fixed optical
count with native point-group inference using identical unscaled modes and
spin. It is separate from production-scaled frequencies, Eckart rates and fits.
The actual retained C2 TS candidate gives a factor-two partition-only difference;
the C2v/C2h endpoint counts agree. Native approximate chirality does not resolve
enantiomer populations, explicit mirror conformers/channels or independent path
degeneracy. Production is not changed on this diagnostic alone. A justified
counting correction must version the recipe/cache and recalibrate comparisons;
model disagreement is not itself an uncertainty interval.

Method-calibration allocation records require an explicit nonempty string
`condition_contract_id`; free-text applicability is not an executable domain.
The specialist converter preserves this declaration's type and the allocation
planner does not stringify an invalid current condition into a matching ID.
Exact-condition records match only a valid declared current ID. Legacy unbound
user/database records remain explicit supplied allocation scenarios. This
safety admission does not authenticate sources, verify condition hashes/method
recipes or produce uncertainty; the automatic producer still requires audited
domain-specific evidence and independent paired predictions.
