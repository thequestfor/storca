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

The current ensemble covers isolated-molecule conformers. Hydrogen-bonded
clusters and explicit-environment snapshots are the next validation-gated
extension; they must be weighted by meaningful conditions rather than added
equally to a monomer spectrum.

### Held-out harmonic scaling

`storca calibrate-ir-scale` is a transparent global-scale search over completed
ORCA normal modes. It holds FWHM fixed, chooses a factor only from explicitly
named training molecules by prominent-band position error, and then reports a
separate holdout score. It never invokes ORCA, edits an existing run, or
automatically changes the default profile. This guards against improving a
small, heterogeneous benchmark by overfitting its digitized references.

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
