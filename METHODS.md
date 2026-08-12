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
convergence. Larger aggregates, trajectory occupancies, and representative DFT
correction transfer remain validation-gated extensions.
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
as sampling features: the stratified poses receive no population weights and
are not labelled as a liquid ensemble. Near-duplicates are removed using fixed
interaction-coordinate, local-frequency, and heavy-atom-RMSD thresholds. The
representative acquisition controller identifies prevalent local
bond/coordination classes and first seeks three independent environments for
each class. Remaining jobs maximize local-frequency, interaction-geometry,
topology, and association diversity. Transfer validation after each cumulative
ORCA batch reprioritizes pending representatives that contain failed classes.
The `auto` and `balanced` hard caps remain four and six environments. These representatives
receive equal geometry-stratum weights rather than electronic-energy Boltzmann
weights. The global ORCA cap reserves their jobs before allocating monomer
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
validation. Stratified equal weights are explicitly not liquid populations.
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

RMG execution is accepted for kinetic interpretation only when retained
artifacts establish normal/recognized completion, requested time coverage, a
parseable mechanism and solver profile, physically usable kinetics, and
quantitative target-loss flux with acceptable numerical closure. Resource- or
coverage-limited models remain incomplete. Rates and `t95` values from the
initial RMG model are screening quantities until the verification procedure
below succeeds.

## Flux-guided ORCA and Arkane verification

Automatic verification is enabled by default for `storca stability` and
`storca stability-ladder`. When the completed RMG model contains target-loss
flux, the verifier ranks retained reactions by integrated target destruction.
Reachability dependencies are part of the ordering: if the controlling route
requires a generated intermediate, a reachable upstream route that forms that
intermediate is verified before the downstream consumer. Only one selected
route is changed per iteration so its effect can be attributed after
repropagation.

For the selected reaction, the generic path layer retains the original RMG
reaction equation and species identities, resolved SMILES, per-species charges
and multiplicities, and candidate common total-spin surfaces. It requires one
explicitly selected common adiabatic surface. A supplied or uniquely resolved
atom map is extended to explicit hydrogens, including a uniquely defined
hydrogen transfer; ambiguous or search-exhausted maps stop before ORCA.

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
the optimized geometries; separated-channel IRC endpoints are matched to the
declared fragment connectivity in either direction before the path status
becomes `verified_transition_state_path`. The independently validated species
also provide the thermochemistry Arkane requires.

The no-saddle branch is deliberately narrower than a visual downhill curve. A
barrierless classification requires normal ORCA completion, the expected path
coverage, a complete finite energy profile, endpoint agreement, and compatible
profiles from at least two independently oriented endpoint complexes. A partial
or single-orientation curve remains `surface_unresolved`. Even a reproducible
`barrierless_capture_candidate` or `barrierless_dissociation_candidate` is path
evidence rather than a rate: the current default Arkane TST adapter rejects it
until a collision-bounded capture and recrossing model is supplied.

For a verified transition-state path, Arkane reads the validated separated
species and TS frequency artifacts, evaluates a temperature/pressure grid that
contains the exact stage condition, and retains the declared bath gas. The
condition rate must have the expected molecularity and units and, for a
bimolecular route, must not exceed the retained collision ceiling. Atom-balanced
routes use explicit artifact validation; missing stationary-point, Hessian,
state, or rate evidence fails closed.

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
