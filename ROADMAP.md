# STORCA Experimental Spectrum Roadmap

## Product objective

STORCA should progress from calculating harmonic molecular spectra to predicting
condition-matched experimental FTIR spectra. Given a chemical system, phase,
temperature, composition, and measurement configuration, the target product
should reproduce the experimentally meaningful quantities:

- band positions and assignments;
- band widths, asymmetry, and overlap;
- relative and, when the required sample data exist, absolute intensities;
- conformer, environment, and speciation contributions;
- the sample and instrument response; and
- uncertainty and applicability diagnostics.

The pipeline must keep two products separate:

1. **Intrinsic spectrum**: the predicted absorption of the modeled molecular or
   condensed-phase system before the measurement transfer function.
2. **Rendered spectrum**: the intrinsic prediction transformed through an
   explicit sample and instrument model.

STORCA must not create experimental resemblance using decorative noise,
arbitrary broadening, hidden per-molecule shifts, or fitting to held-out test
spectra. An exact trace is not a well-defined target without the corresponding
phase, temperature, concentration, path length or ATR configuration,
resolution, apodization, and other relevant experimental metadata.

## Target scientific architecture

```text
chemical identity and experimental conditions
  -> conformers, species, clusters, solvent or crystal environments
  -> finite-temperature sampling and defensible populations
  -> mode-character identification across changing geometries
  -> multi-fidelity frequencies and intensities
  -> anharmonicity and resonance treatment
  -> intrinsic line shapes and uncertainty
  -> sample and instrument transfer function
  -> condition-matched validation
```

This architecture generalizes the current hydrogen-bond work. O-H and N-H
stretches remain important validation cases, but the implementation must also
track carbonyls, C-O and C-N modes, C-H stretches and bends, aromatic and ring
modes, torsions, and the coupled fingerprint region.

## Current baseline

The repository already contains important pieces of this architecture:

- reproducible ORCA calculation contracts, caching, and resume behavior;
- GOAT conformer generation and monomer spectra;
- diverse xTB dimer and trimer generation for small protic systems;
- topology-aware environment selection under a hard ORCA-job budget;
- local X-H mode projection and matching across cluster sizes;
- environment-sufficiency and effective-sample-size gates;
- adaptive environment convergence diagnostics;
- separate monomer, dimer, network, raw, ensemble, and final spectrum outputs;
- reliability warnings for poorly stationary snapshot Hessians; and
- an initial condition-aware benchmark and calibration framework.

Important limitations remain:

- environment populations are generally static strata, not trajectory
  occupancies or liquid free-energy populations;
- the sampled clusters do not yet constitute a finite-temperature bulk liquid,
  solution, or crystal ensemble;
- equal-size mode identity now uses a general internal-coordinate description,
  but mixed-size environment matching remains specialized toward local X-H
  stretches;
- full ORCA Hessians still carry too much of the environment workload;
- representative DFT corrections are transferred to the sampled xTB ensemble
  for validated local X--H classes, but remain diagnostic until order-robust
  convergence passes;
- anharmonicity, resonances, and vibrational lifetimes are not modeled at the
  level required for general experimental line shapes; and
- sample and instrument metadata are not yet a complete rendering contract.

## Status vocabulary

Roadmap items use the following states:

- **Available**: implemented and exercised by the current test suite.
- **Partial**: useful implementation exists, but does not yet satisfy the target
  scientific contract.
- **Planned**: not yet implemented as a supported pipeline feature.

No partial component should be presented as a completed physical model. In
particular, vacuum electronic-energy Boltzmann weights must never be labeled as
liquid populations, and a residual display width must never be labeled as an
environment width.

## Milestone 0: Harden the current environment pipeline

**Status: Available, with continuing hardening.**

Purpose: preserve the reliability gains already made and provide a stable base
for broader band classes.

Deliverables:

- versioned environment-sampling, convergence, and band schemas;
- environment sufficiency based on independent clusters, nonzero frequency
  variance, geometric diversity, effective sample size, and matching
  confidence;
- explicit `insufficient_environment_sampling` and unconverged states;
- distinct calculated, residual, and instrument width components;
- deterministic representative selection and resumable calculations; and
- regression fixtures for methanol and a non-hydrogen-bonded control.

Exit gate: identical inputs and cached calculations reproduce the same selected
environments, diagnostics, bands, and spectra; insufficient evidence fails
closed rather than producing a realistic-looking width.

## Milestone 1: Define experimental contracts and general band identity

**Status: Partial. Versioned reference contracts, intrinsic/rendered separation,
and equal-size internal-coordinate fingerprints are available; generalized
mixed-size atom mapping and locked condition-complete partitions remain.**

Purpose: make every comparison condition-matched and replace mode-number
tracking with chemically meaningful mode character.

Deliverables:

- an `ExperimentalCondition` contract containing phase, temperature, pressure,
  composition, concentration, and measurement configuration when known;
- a versioned reference manifest with source provenance, digitization method,
  uncertainty, and immutable train, validation, and holdout membership;
- generalized internal-coordinate participation fingerprints for stretches,
  bends, torsions, ring modes, and collective modes;
- atom-mapped mode-overlap matching across conformers, environments, and
  cluster sizes;
- stable band classes that survive mode reordering and avoided crossings; and
- baseline metrics for centers, FWHM, asymmetry, integrated regional intensity,
  whole-spectrum overlap, missing bands, spurious bands, and accuracy per ORCA
  job.

Exit gate: controlled fixtures correctly track carbonyl, C-O, C-N, C-H,
aromatic, X-H, and coupled fingerprint modes through geometry changes and mode
crossings. Every validation spectrum either has a compatible condition contract
or is explicitly excluded from quantitative scoring.

## Milestone 2: Build large, finite-temperature environment ensembles

**Status: Partial for restrained dimer and trimer pose sampling; planned for
trajectory populations.**

Purpose: use many inexpensive configurations to represent the environment
distribution while keeping quantum-chemistry cost bounded.

Deliverables:

- interaction-site discovery beyond hydrogen bonds, including strong dipoles,
  ionic contacts, carbonyl solvation, and relevant pi interactions;
- restrained xTB scans and short finite-temperature trajectories;
- decorrelated snapshot extraction and clustering by interaction geometry,
  molecular orientation, local field, energy, and estimated mode shifts;
- dimers, trimers, and selected larger central-molecule environments subject to
  atom-count and ORCA-job caps;
- population weights in this order of preference: trajectory occupancy,
  free-energy-informed xTB populations, then explicitly labeled stratified
  equal weights; and
- convergence over the inexpensive ensemble before expensive representatives
  are requested.

The large ensemble is essential. The small ORCA set is not intended to be the
ensemble; it is the calibration set for the ensemble.

Exit gate: repeated seeded runs give statistically compatible centers, widths,
and populations; effective sample size and structural coverage pass configured
thresholds; and no vacuum-energy weighting is described as a bulk population.

## Milestone 3: Transfer DFT accuracy across the xTB ensemble

**Status: Partial: snapshot Hessians and fail-closed local X--H transfer are implemented.**

Purpose: calculate dozens or hundreds of environment contributions without an
ORCA Hessian for every snapshot.

Deliverables:

- environment-preserving DFT micro-relaxation and gradient diagnostics for
  representative snapshots;
- xTB local-mode or full-mode frequency estimates for all retained snapshots;
- full ORCA Hessians for a small number of representatives;
- local-mode finite differences for important modes where a full Hessian is not
  cost-effective;
- mode-character- and environment-class-aware frequency corrections such as

  ```text
  nu_estimated(snapshot) = nu_xTB(snapshot)
                         + nu_DFT(representative)
                         - nu_xTB(representative)
  ```

- an analogous, separately validated treatment for intensity corrections; and
- transfer uncertainty that increases outside the representative set's domain.

Representative selection must cover populated clusters and retain meaningful
finite-temperature higher-energy configurations. Selecting only the
lowest-energy geometries is not acceptable.

Exit gate: cross-validation on withheld DFT representatives shows that transfer
improves upon uncorrected xTB and meets a predeclared error target. Poorly
stationary or out-of-domain representatives are reported, not silently used.

Implemented now: resumable unrestrained GFN2-xTB snapshot Hessians, parsed IR
intensities and normal coordinates, internal-coordinate fingerprints, same-
geometry xTB/ORCA local X--H pairing, additive corrections, applicability
distance and uncertainty, leave-one-representative-out validation, and cheap-
ensemble convergence diagnostics with stratified ordering and candidate
bootstrap intervals. A gated hybrid diagnostic combines transferred local
frequencies, representative-ORCA class intensities, and representative ORCA
nonlocal modes. Environment-preserving constrained DFT refinement and an
unrestrained gradient stationarity gate are also available. Projected local
X--H finite differences now use two mass-weighted displacement sizes, analytic
ORCA gradients and dipoles, Richardson extrapolation, a displacement-stability
gate, and an actual-process invocation ledger. Validated artifacts can replace
the matching local snapshot-Hessian transfer pair. Representative acquisition is now
mode-class-aware: prevalent sampled classes receive coverage before remaining
jobs are assigned by diversity, and failed transfer classes reprioritize the
pending ORCA queue after each batch. Remaining: generalized mixed-size non-X--H
transfer, automatic budget-aware execution of the implemented local-mode
finite differences for structures that fail the gradient gate,
snapshot-specific validated intensity transfer, and trajectory
occupancy weights.

Static ensemble extension is now available as a resumable command. It preserves
the initial ensemble, adds independently seeded balanced 20-candidate rounds,
reuses all prior xTB optimizations and Hessians, evaluates convergence on real
cumulative acquisition rounds, reports bootstrap precision separately, and
requests additional ORCA coverage only for out-of-domain corrected snapshots.
The retained methanol 40→60→80→100 exercise reached numerical round stability
without reaching the declared bootstrap-precision gate, demonstrating the
intended fail-closed distinction.

## Milestone 4: Add adaptive multi-fidelity convergence

**Status: Partial for the current environment representatives.**

Purpose: spend the next ORCA job only when it is likely to change the answer.

The convergence controller will evaluate batches using:

- important band centers;
- frequency standard deviations and equivalent FWHM values;
- integrated intensities;
- whole-spectrum overlap;
- band-class discovery; and
- transfer-model uncertainty and coverage.

The existing target is two consecutive stable batches with center changes below
5 cm-1, width changes below 10 cm-1 or 10%, overlap above 0.98, and no new
significant band class. These thresholds must remain configurable and recorded
in the calculation contract. Reaching the hard budget first produces an
unconverged result, not a false success.

Exit gate: convergence decisions are reproducible, traceable to batch metrics,
and validated against deliberately extended calculations.

## Milestone 5: Model anharmonicity and resonances

**Status: Planned.**

Purpose: correct systematic harmonic errors and recover spectral structure that
cannot arise from harmonic fundamentals alone.

Deliverables:

- a validated scaled-harmonic baseline by method and mode class;
- one-dimensional local potentials for selected strongly anharmonic modes;
- selective VPT2 or an equivalent treatment where stationary structures and
  cost permit it;
- automatic resonance-risk detection;
- Fermi resonance, overtone, and combination-band handling in applicable
  fidelity tiers; and
- temperature-dependent hot-band contributions where experimentally relevant.

This layer must be selective: applying an expensive anharmonic method everywhere
would undermine the ensemble strategy, while ignoring resonance physics would
leave important regions experimentally wrong.

Exit gate: each enabled correction improves held-out condition-matched spectra
for its declared applicability class without increasing spurious-band rates.

## Milestone 6: Predict physical intensities and intrinsic line shapes

**Status: Partial for harmonic intensities and decomposed display widths;
planned for physical lifetime models.**

Purpose: replace generic peak broadening with separately interpretable sources
of line shape.

Deliverables:

- intensity changes across conformers and environments;
- unit-area line profiles whose integrated intensities are conserved;
- separate conformer, environment, homogeneous lifetime, rotational, and
  instrument contributions;
- asymmetric and non-Gaussian environment distributions when supported by the
  sampled frequencies;
- mode-class lifetime or homogeneous-width models with provenance and
  uncertainty; and
- absolute absorbance only when concentration, path length, and the required
  optical data are available.

Exit gate: regional integrated intensity is conserved through broadening,
component widths are auditable, and held-out width and asymmetry metrics improve
over the residual-width baseline.

## Milestone 7: Add phase, sample, and instrument models

**Status: Partial for generic residual and instrument broadening.**

Purpose: turn the intrinsic prediction into the spectrum expected from a stated
measurement rather than from an unspecified instrument.

Deliverables:

- distinct applicability paths for gas, solution, neat liquid, and solid or
  crystal samples;
- solution composition and explicit-solvent contracts;
- crystal or periodic environment support when cluster models are inadequate;
- transmission rendering from concentration and path length;
- ATR rendering from incidence angle, crystal properties, and refractive-index
  information;
- resolution, line-spread function, apodization, sampling grid, and other
  instrument metadata; and
- `spectrum_intrinsic.csv` and `spectrum.csv` as first-class, separate outputs.

Exit gate: one intrinsic spectrum can be rendered through multiple known
measurement configurations with predictable, independently tested changes. A
missing measurement contract produces a standardized simulation, clearly
labeled as such.

## Milestone 8: Run the benchmark and calibration program

**Status: Partial.**

The initial validation set should include methanol, ethanol, water, acetic acid,
methylamine, and a non-hydrogen-bonded control such as hexane. It should then
expand across functional groups, molecular sizes, phases, temperatures, and
measurement modes.

For every case, evaluate:

- band-center error;
- FWHM and envelope-asymmetry error;
- whole-spectrum and regional overlap;
- integrated regional intensity;
- missing and spurious bands;
- uncertainty coverage;
- computational cost; and
- accuracy gained per ORCA job.

Calibration data, validation data, and held-out data must be separated by
molecule or chemical family as appropriate. Calibration versions must record
their training references and content hashes. Thresholds for promotion must be
locked before the held-out evaluation is run.

Exit gate: an ablation report shows the contribution of environment sampling,
DFT transfer, anharmonicity, line shapes, and measurement rendering. Improvements
must persist outside the calibration set.

## Milestone 9: Enable the experimental-spectrum pipeline by default

**Status: Planned.**

The new pipeline becomes the default only when all of the following are true:

- condition-matched held-out center, width, overlap, and intensity metrics
  improve over the current monomer and residual-width baseline;
- non-hydrogen-bonded controls do not acquire artificial broad bands;
- missing- and spurious-band rates do not regress materially;
- reported uncertainty has acceptable empirical coverage;
- convergence and sufficiency statuses remain honest under hard budgets;
- calculations are deterministic or statistically reproducible as applicable;
- cache, resume, schema migration, and provenance tests pass; and
- the improvement per ORCA job justifies the selected default fidelity.

Until those gates pass, the feature remains opt-in and its outputs identify the
experimental capabilities that are incomplete.

## Target fidelity tiers

The exact counts remain capped by the calculation contract and may be reduced by
atom count, charge, method support, or the hard ORCA-job budget.

| Fidelity | Environment sampling target | Expensive reference target | Physical model target |
| --- | --- | --- | --- |
| `fast` | none or monomer conformers | monomers only | harmonic intrinsic spectrum |
| `auto` | roughly 30-50 optimized poses or snapshots | 3-4 representative environments | validated transfer where available |
| `balanced` | roughly 50-100 poses or a short trajectory | 5-7 representative environments | broader transfer and selected corrections |
| `accurate` | longer dynamics with dimers, trimers, and selected larger environments | 8-12 representative environments | selective anharmonicity, resonances, and richer phase treatment |

These are target definitions, not a claim that every row is currently
implemented. Every run must record what was actually executed.

## Target output contract

The eventual output set should include:

- `experimental-condition.json`;
- `environment-sampling.json`;
- `environment-convergence.json`;
- `environment-geometries/`;
- `mode-character.json` or equivalent per-band mode data;
- `spectrum_monomer.csv`;
- `spectrum_dimer.csv`;
- `spectrum_network.csv`;
- `spectrum_intrinsic.csv`;
- `spectrum.csv`;
- `uncertainty.json`; and
- `evaluation.json` when a compatible reference is supplied.

Each reported band should eventually contain at least:

```json
{
  "center_cm-1": 3350.0,
  "center_uncertainty_cm-1": 18.0,
  "conformer_fwhm_cm-1": 8.0,
  "environment_fwhm_cm-1": 210.0,
  "homogeneous_fwhm_cm-1": 12.0,
  "residual_fwhm_cm-1": 6.0,
  "instrument_resolution_cm-1": 4.0,
  "integrated_intensity": 1.0,
  "mode_character": ["O-H stretch", "hydrogen-bonded"],
  "environments": 7,
  "effective_sample_size": 5.4,
  "population_source": "trajectory_occupancy",
  "converged": true,
  "reliability": "validated"
}
```

Fields unsupported by a run must be absent or explicitly unavailable; they must
not be filled with plausible defaults.

## Immediate implementation sequence

The next work should be delivered in small, testable increments:

1. **Experimental and reference contracts**: finish the condition schema,
   reference provenance, fixed data splits, and baseline evaluation report.
2. **General mode-character layer**: implement internal-coordinate fingerprints
   and matching fixtures beyond X-H stretches.
3. **Intrinsic/rendered output split**: make the scientific prediction distinct
   from residual and instrument rendering throughout the schema and reports.
4. **DFT-transfer pilot**: add environment-preserving representative
   micro-relaxation, xTB snapshot frequencies, correction transfer, and
   cross-validation on methanol plus controls.
5. **Short xTB dynamics**: use decorrelated snapshot occupancy for populations
   and extend adaptive convergence to the inexpensive ensemble.
6. **Generalize environment chemistry**: add non-hydrogen-bond interactions and
   validate carbonyl, amine, aromatic, and hydrocarbon cases.
7. **Add selective anharmonic and line-shape layers** only after the ensemble and
   transfer errors are measurable.
8. **Calibrate and run the locked held-out benchmark** before changing defaults.

The order is deliberate: a much larger ensemble without reliable mode identity
cannot be combined correctly; expensive physics without a condition-matched
benchmark cannot be judged; and instrument rendering cannot repair an
incorrect intrinsic spectrum.

## Cross-cutting engineering rules

- Fail closed when sampling, matching, convergence, or transfer evidence is
  insufficient.
- Keep hard atom-count, wall-time, and ORCA-job limits enforceable.
- Preserve deterministic plans, seeds, content hashes, and calculation
  provenance.
- Make cached and resumed runs scientifically equivalent to uninterrupted runs.
- Attach uncertainty and an applicability status to every calibrated result.
- Prevent training-reference leakage into validation and holdout evaluation.
- Require an ablation and held-out improvement before enabling a new physical
  layer by default.
- Keep raw, intrinsic, and rendered outputs available so every transformation is
  auditable.

## What success means

Success is not a curve that merely looks more experimental. It is a spectrum
whose centers, distributions, intensities, line shapes, and measurement effects
are traceable to explicit models, whose limitations are visible, and whose
accuracy is demonstrated on condition-matched molecules that were not used to
calibrate it.
