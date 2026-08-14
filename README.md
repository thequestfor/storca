# STORCA

### Gated direct-local-DFT and bulk-liquid development

The O--H environment workflow now has three fail-closed implementation levels:

- `storca.direct_local_dft` acquires representatives 4--6 with a persistent
  process ledger and assembles `direct-local-dft-ensemble.json`, separate
  intrinsic/transmission/comparison CSV files, and NIST regional metrics.  It
  excludes every harmonic mode in the configured O--H window and replaces that
  window only with validated ORCA projected finite-difference modes; it never
  imports or applies the xTB frequency-transfer model.
- `storca.adaptive_local_modes` selects later representatives from trajectory
  occupancy, coordination, H-bond geometry, topology, geometric diversity, and
  distance from calculated environments.  Batches contain three representatives,
  and uncertainty resamples whole trajectory clusters.  Two consecutive batches
  must pass center, FWHM, envelope-overlap, CI, effective-sample-size, and class
  coverage gates.
- `storca.bulk_embedding` builds a seeded periodic methanol box at 298.15 K
  density, validates NVT observables and autocorrelation, extracts the central
  methanol plus its geometric first H-bond shell, and writes reproducible ORCA
  external-point-charge inputs.  Embedded local modes use four gradients on
  eight cores per central O--H bond.  Promotion requires two independent seeds
  and improvement on held-out trajectory blocks.

`evaluate_spectroscopy_gate_sequence(run_dir)` in
`storca.spectroscopy_gates` reads the retained artifacts and authorizes only the
next justified expensive stage.  Missing references, fewer than three valid
environments in a coordination class, or incomplete finite-difference evidence
block promotion.

The six-representative planner also checks the arithmetic of the calculation
request before launching ORCA.  At four gradients per bond, nine new bonds cost
36 invocations, not 24; a 24-invocation hard cap therefore fails closed unless
the retained representative manifest contains at most six uncached bonds.

STORCA is a small, reproducible workflow for ORCA geometry optimizations and
vibrational frequency calculations. Each calculation writes to its own run
directory rather than cluttering the repository root.

The staged plan for advancing from harmonic calculations to condition-matched
experimental-spectrum predictions is documented in [ROADMAP.md](ROADMAP.md).

## Quick start

Create an isolated Python environment, install the package, then check
prerequisites:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[chemistry]'
storca doctor
```

`doctor` checks whether ORCA, xTB, Open Babel, RMG, RDKit, and the Python Open
Babel binding are discoverable in the active environment. On Debian/Ubuntu,
`/usr/bin/orca` may be the unrelated GNOME screen reader; it is not the ORCA
quantum-chemistry program. Obtain ORCA from FACCTs, then set
`STORCA_ORCA_BIN` to its executable path.

For a machine-specific ORCA location, set `STORCA_ORCA_BIN` to the executable
path instead of changing project code.

Run from an XYZ geometry:

```bash
storca run ethanol.xyz --charge 0 --multiplicity 1 --cores 4
```

Try the included small geometry with `storca run examples/water.xyz`.

Or generate an initial geometry from SMILES (requires the `chemistry` extra):

```bash
python -m pip install -e '.[chemistry]'
storca run --smiles CCO --name ethanol
```

Outputs are placed under `runs/<name>-<timestamp>/` and include the input
geometry, ORCA inputs/outputs, optimized geometry, and `metadata.json`.

## Scope

The core molecular workflow is ORCA optimization followed by a frequency
calculation. The stability suite additionally couples bounded homogeneous-gas
RMG models to default, flux-guided ORCA path verification, Arkane kinetics, and
full repaired-model propagation under an immutable condition contract.
`typerint.py` is archived prototype code and is not a supported interface.
Use `storca --help`; Multiwfn is not used anywhere in the supported path.

See [METHODS.md](METHODS.md) and [LIMITATIONS.md](LIMITATIONS.md) before using
the results for scientific decisions.

## Molecule records

Create a portable, offline starting record before launching a calculation:

```bash
storca structure --smiles CCO --name ethanol
```

This writes `molecule.json`, a 2D PNG, and an RDKit/UFF initial XYZ geometry
in its own run directory. The geometry is only a calculation seed; ORCA must
optimize it before geometry-dependent interpretation.

## Orbital and qualitative reactivity summary

Extract frontier-orbital information from an existing ORCA output without any
post-processing dependency:

```bash
storca orbitals runs/ethanol-.../sp.out --json-output orbitals.json
```

For a molecule-level, explicitly qualitative frontier-orbital summary:

```bash
storca reactivity runs/ethanol-.../freq.out --json-output reactivity.json
```

It reports HOMO/LUMO gap and Koopmans-like chemical-potential,
electronegativity, hardness, softness, and electrophilicity proxies. It does
not provide atom-resolved reactivity, Fukui
functions, redox potentials, or safety predictions.

## IR spectrum prediction

With RDKit installed, generate, optimize, and Boltzmann-weight conformers from
a SMILES string:

```bash
storca spectrum --smiles CCO --name ethanol --cores 4 --conformer-engine goat
```

The command writes `spectrum.csv`, `spectrum.png`, and per-conformer results in its run
directory. It retains failed calculations for diagnosis and only weights
conformers whose optimization converged and whose frequency calculation is a
local minimum. Gaussian bands preserve the integrated ORCA band strength, and
the displayed relative trace is normalized to unit peak. The spectrum is a
harmonic, broadened prediction, not an experimental replacement.

GOAT/xTB is the default engine: it searches a large low-cost ensemble and
passes only conformers covering 95% of the GOAT population (up to 10 by
default) into higher-level ORCA calculations. Run STORCA from the single
`orca-auto-py312` Conda environment so ORCA, xTB, RDKit, and STORCA share one
PATH. Use `--conformer-engine rdkit` only as an explicit fallback.

## Molecule description

Describe a SMILES locally, without a network lookup:

```bash
storca describe --smiles CCO --json-output ethanol-description.json
```

This reports canonical structure identifiers, formula, mass, structural counts,
RDKit descriptors, detected functional-group patterns, and a generated
description. The JSON identifies all values as local RDKit-derived descriptors;
it does not claim they are measured properties. See
[LEGACY_INVENTORY.md](LEGACY_INVENTORY.md) for the migration status of older
utilities.

To intentionally obtain identifiers, names, database properties, and synonyms
from PubChem, use the separate networked command:

```bash
storca enrich --smiles CCO --source pubchem --json-output ethanol-pubchem.json
```

The result records its source URL and retrieval time. It is not an SDS, hazard
assessment, or experimental-property guarantee.

If a spectrum calculation is interrupted after some conformers finish, recover
the available results without rerunning them:

```bash
storca spectrum-finalize runs/spectrum-YYYYMMDDTHHMMSSZ
```

To continue missing ORCA jobs while reusing completed ones:

```bash
storca spectrum-resume runs/spectrum-YYYYMMDDTHHMMSSZ --cores 4
```

Resume and finalize inherit the original method, charge, spin, spectrum model,
temperature, scaling, and display settings. Per-conformer geometry/electronic
signatures prevent completed jobs from being reused under an incompatible
contract. Once `experimental-condition.json` exists, resume/finalize will not
mutate its sample or measurement settings in place; create a new run for a
different measurement contract.

## Spectrum benchmarking

Compare against a numerical two-column reference CSV (wavenumber plus the same
intensity convention as the calculated CSV), recording measurement conditions:

```bash
storca benchmark runs/example/spectrum.csv reference.csv \
  --phase liquid --solvent neat --temperature 298.15 --measurement ATR
```

This creates a JSON report with overlap, RMSE, and correlation. Use it to
calibrate method profiles; it does not replace chemically informed review.

For a transparent peak-by-peak comparison—without shifting or fitting either
spectrum—use:

```bash
storca spectrum-analyze runs/example/spectrum.csv reference.csv \
  --output runs/example/peak-analysis.json \
  --overlay runs/example/peak-analysis.png
```

The report lists matched and unmatched bands, position errors, and relative
intensity differences. The supplied SDBS image digitizer uses its non-linear
printed wavenumber ticks; use `digitize-ir --axis linear` only for an image
whose x axis is truly linear.

New spectrum runs record the harmonic method profile, scale-factor source,
conformer workflow, temperature, and line-shape settings in `metadata.json`.
The default `b3lyp-def2-svp` profile currently supplies the existing 0.970
project baseline. Override it intentionally with `--scale-factor`, for example
`--scale-factor 0.975`; this is recorded as a user override.

Two unscaled ORCA composite-method candidates are also available for the
validation bakeoff: `b3lyp-3c` and `r2scan-3c`. Select one deliberately, for
example:

```bash
storca spectrum --smiles 'c1ccccc1' --method-profile r2scan-3c --cores 4
```

They are not recommended as replacements for the default yet. Their scale
factors and performance must be accepted from a diverse benchmark set before a
project default changes.

### Held-out scale calibration

`calibrate-ir-scale` reuses completed ORCA frequency outputs; it does not run
ORCA or overwrite a run. Provide disjoint training and holdout benchmark IDs:

```bash
storca calibrate-ir-scale benchmarks/manifest.json \
  --training-ids acetic-anhydride-ccl4,hexane-ccl4,ethylamine-liquid-film,ethanol-nist,acetaldehyde-nist \
  --holdout-ids benzene-liquid-film,glycolaldehyde-dimer-kbr,dimethylamine-nist \
  --output benchmarks/scale-calibration.json
```

The selected factor minimizes prominent-band position error on the training
set, with correlation only as a tiebreaker. Accept it only when the untouched
holdout metrics also improve and the reference conditions match the intended
use; STORCA never changes a method profile automatically.
The benchmark manifest is schema v2: numerical references are content-hashed,
carry provenance and condition metadata, and declare a dataset partition. The
current heterogeneous references remain `development_unassigned`; they cannot
support a locked held-out claim until condition-complete partitions are
declared before calibration.

## Practical IR display

The normal `spectrum` command preserves the raw, conformer-weighted harmonic
calculation. Add `--spectrum-model practical` to create a second, explicitly
rule-calibrated display for common organic experimental IR spectra:

```bash
storca spectrum --smiles 'O1CC(O)OCC1O' --spectrum-model practical --cores 4
```

The calculated spectrum remains in `spectrum_raw.csv` and `spectrum_raw.png`.
The displayed `spectrum.csv` records the `ambient-organic-v0.3` rules in
`metadata.json`. At present, this profile has retained only its
flexible-molecule fingerprint rule. The initial broad O-H envelope helped a
KBr-disc benchmark but degraded an ethanol reference, so it is not applied
without an explicit sample-aware model. These rules are an empirical display
calibration—not a solvent, crystal, or
first-principles linewidth calculation—and will be retained only when they
improve held-out benchmarks.

## Calculation-first experimental FTIR

Use the experimental model when the desired output is a simulated measurement
rather than a uniformly broadened harmonic trace:

```bash
storca spectrum --smiles CCO --spectrum-model experimental \
  --phase liquid --measurement atr --instrument-resolution 4 \
  --apodization happ-genzel --pressure 1 \
  --solvent neat --atr-crystal diamond --fidelity auto --cores 4
```

To permit the validated projected-local-mode fallback for nonstationary
environment representatives, supply a separate process allowance, for example
`--local-mode-orca-invocations 36`. Zero is the default and performs no
additional ORCA work. The planner spends four gradient invocations per selected
bond and starts a mode class only when it can cover three distinct
representatives; the persistent run-level ledger is a hard cap.

Extend a completed cheap environment ensemble without repeating its retained
calculations using:

```bash
storca spectrum-extend-environments runs/spectrum-YYYYMMDDTHHMMSSZ \
  --maximum-candidates 100 --batch-size 20 --cores 4
```

Extension rounds use independent deterministic seeds and balanced dimer,
linear-trimer, cyclic-trimer, strong, intermediate, and weak strata. Numerical
batch stability and bootstrap precision are separate gates: stable centers and
widths are not labelled converged while their 95% intervals remain too broad.

The resolved sample and measurement settings are retained in the versioned
`experimental-condition.json` contract. Optional CLI fields record composition,
concentration, transmission path length, ATR incidence angle, and refractive
index when those values are known. Missing values remain explicitly unknown.

This model uses the frequency distribution calculated across the
Boltzmann-weighted conformer ensemble as physical inhomogeneous broadening.
With `--fidelity auto` or `balanced`, a neutral closed-shell neat liquid or
solid containing both donor and acceptor sites receives deterministic,
coverage-oriented dimer sampling. Small protic molecules (at most four heavy
atoms and 36 atoms in the trimer) additionally receive linear donor-acceptor
chains and cyclic trimers where geometrically plausible. The `auto` and
`balanced` tiers propose 40 and 75 total poses respectively across 1.6--2.5 Angstrom H--acceptor distances,
120--180 degree donor--H--acceptor angles, hydrogen-bond-axis rotations, and
alternative framework orientations. Each pose receives an isolated, weakly
restrained GFN2-xTB optimization with a content-hashed resume contract. Four
(`auto`) or six (`balanced`) diverse poses then seed restrained NVT GFN2-xTB
trajectories at the requested spectrum temperature (298.15 K by default). The
bounded defaults are 2 ps (`auto`) and 4 ps (`balanced`), with coordinates
propagated at a conservative 0.5 fs step with hydrogen mass 4 and written every
10 fs. STORCA rejects xTB's nominal-success emergency exits and truncated
trajectories, discards the first 25%, estimates statistical
inefficiency from the restrained-contact time series, and extracts frames at
no less than that autocorrelation stride. The 40/75-frame target ensemble
receives equal-time snapshot occupancies. These are conditional occupancies
within the finite restrained multi-seed protocol, not unbiased bulk-liquid
equilibrium populations. If trajectory sampling fails, the workflow records
the downgrade and falls back to the static coverage ensemble.

Each retained trajectory snapshot receives a separately cached, unrestrained
GFN2-xTB numerical Hessian; the sampling restraints are never included in the
force constants. The resulting local X--H frequencies are included in the
diversity feature space. This inexpensive ensemble consumes no ORCA budget.
Near-duplicates are removed using interaction geometry and aligned heavy-atom
RMSD while preserving their occupancy mass. A coverage-constrained acquisition controller then
selects four (`auto`) or six (`balanced`) representatives. It first targets
three independent representatives for every prevalent local mode/coordination
class, then spends remaining jobs on frequency, geometry, topology, and
association diversity. After each ORCA batch, failed transfer classes are
moved ahead of redundant pending representatives. Representatives receive the
summed occupancy of the decorrelated snapshots assigned to their clusters.

Environment jobs are reserved before GOAT monomer jobs: the default 12-job cap
allocates 8 monomers plus 4 environments for `auto`, or 6 plus 6 for
`balanced`. ORCA evaluates a Hessian directly at each selected restrained-xTB
geometry instead of performing an unconstrained optimization that could
collapse the environments into one minimum. These are recorded as snapshot
Hessians, not unconstrained DFT stationary points. A separate
environment-preserving DFT refinement contract constrains interaction geometry
during optimization and then evaluates an unrestrained gradient. Full-Hessian
use is permitted only when both gradient gates pass.

Normal modes across equal-sized configurations are matched primarily by
internal-coordinate fingerprints containing stretch, bend, torsion, and ring
participation. Geometry-aligned Cartesian displacement overlap remains a
fallback. Mixed dimer/trimer X--H bands are followed through target-local
bond-stretch projections, avoiding an invalid comparison between differently
sized Hessian vectors. The model then
applies a small, recorded residual linewidth, the selected measurement
geometry, and a finite-resolution FTIR instrument line shape. ATR
uses a first-order wavelength-dependent penetration-depth response;
transmission and gas-cell modes do not apply it. The `happ-genzel` choice is a
positive central-lobe approximation, not an emulation of particular instrument
firmware.

A network frequency distribution is exposed as a calculated environment
width only when each band has at least three independent geometry clusters,
nonzero frequency variance, effective sample size of at least 2.5, documented
hydrogen-bond geometry diversity, and normal-mode overlap of at least 0.45
covering at least 80% of its population. These are versioned minimum-evidence
thresholds, not a claim that three environments are a converged liquid model.
If any requirement fails, that band's calculated environment FWHM is reported
as zero with `width_status=insufficient_environment_sampling`. Its sampled
frequencies remain in a diagnostic spectrum, while the displayed band is
collapsed to its weighted center before residual and instrument response are
applied.
Local X--H observations are first grouped by independent geometry for
effective-sample-size and diversity checks, so several oscillators in one
trimer cannot inflate the number of environments. Hydrogen-bonded and
non-donating X--H oscillators are evaluated as separate spectral band classes;
free terminal modes are not averaged into the bonded-network center.
ORCA representatives are evaluated in cumulative diversity-first batches:
`2,3,4` for `auto` and `2,4,6` for `balanced`, capped by the selected and hard
job budgets. After each batch STORCA compares important band centers, raw
frequency-distribution FWHM values, integrated intensity, and a fixed 2 cm^-1
absorption-strength spectrum. Convergence requires center changes below
5 cm^-1, width changes below 10 cm^-1 or 10%, whole-spectrum cosine overlap
above 0.98, adequate mode overlap, and no new significant band class for two
consecutive comparisons. If the hard budget is reached first, a sufficient
distribution remains visible but is labelled
`width_status=environment_width_unconverged`; it is never reported as a
converged environment width.

After representative ORCA Hessians exist, STORCA pairs xTB and ORCA local X--H
modes at the same geometries and tests an additive DFT correction by leaving
each representative out in turn. Transfer requires at least three independent
representatives per mode/coordination class, mode-character similarity of at
least 0.50, at least 20% lower error than raw xTB, a high-frequency MAE no
larger than 30 cm^-1, and no withheld error above 75 cm^-1. Out-of-domain or
failed classes remain diagnostic; the representative ORCA spectrum stays the
display basis. xTB intensities are not treated as DFT intensities, and transfer
does not convert stratified geometry weights into liquid populations.

The run keeps the layers separate:

- `spectrum_raw.csv/png`: the original harmonic display using `--fwhm`.
- `spectrum_ensemble.csv/png`: calculated conformer distribution with the
  phase-default or user-supplied `--residual-fwhm`.
- `spectrum_intrinsic.csv/png`: calculated absorption strength before ATR,
  transmission, or finite-resolution instrument response.
- `spectrum.csv/png`: the simulated measurement after geometry and instrument
  response.
- `experimental-condition.json`: immutable, versioned sample and measurement
  metadata used to render and evaluate the spectrum.
- `mode-character.json`: internal-coordinate character and matching confidence
  for each retained band.
- `ensemble-bands.json`: per-mode calculated mean, standard deviation, and
  equivalent ensemble FWHM, effective sample size, sufficiency decision, and
  mode-matching confidence.
- `environment-sampling.json`: restrained-xTB sampling configuration,
  candidate records, the explicit no-population warning, versioned sufficiency
  thresholds, per-band decisions, failure reasons, and display-width policy.
- `xtb-trajectories.json` and `environment-trajectories/`: restrained NVT
  contracts, raw trajectories, decorrelation diagnostics, extracted snapshots,
  and conditional occupancy weights.
- `environment-convergence.json`: cumulative batch metrics, threshold results,
  consecutive-pass count, manifest hash, budget state, and stop reason.
- `environment-acquisition.json`: sampled mode-class prevalence, required and
  completed DFT coverage, representative decisions, adaptive reorderings,
  transfer feedback after each ORCA batch, and hard-budget stop state.
- Snapshot-Hessian reliability records the material imaginary-mode count and
  labels strongly nonstationary full Hessians as diagnostic-only. Target-local
  X--H projections retain their own matching-confidence assessment.
- `environment-geometries/`: per-candidate initial geometry, xTB restraint
  input, calculation contract, output, optimized geometry, sampling record,
  and an `xtb-frequency/` unrestrained snapshot-Hessian cache.
- `xtb-snapshot-frequencies.json`: parsed xTB modes, local-mode assignments,
  mode fingerprints, intensities, and imaginary-curvature diagnostics.
- `frequency-transfer.json` and `frequency-transfer-validation.json`: paired
  representative corrections, applicability distances, uncertainty estimates,
  and leave-one-representative-out gates.
- `frequency-transfer-convergence.json`: cumulative cheap-ensemble center,
  width, intensity, and whole-spectrum convergence diagnostics under three
  deterministic stratified orderings, with candidate bootstrap intervals.
- `spectrum_xtb_environment.csv` and
  `spectrum_dft_transferred_intrinsic.csv`: explicitly diagnostic raw-xTB and
  validated/partially transferred local-mode ensembles.
- `spectrum_hybrid_multifidelity_intrinsic.csv`: representative ORCA nonlocal
  modes plus corrected xTB-snapshot local frequencies and representative-ORCA
  class intensities. It remains diagnostic unless validation, convergence, and
  complete local-intensity coverage all pass.
- Representative `local-mode-finite-differences/` directories contain
  center-of-mass-preserving bond displacements, ORCA gradient/dipole jobs, a
  versioned calculation contract, `local-modes.json`, and an actual-process
  `orca-invocations.json` hard-cap ledger. Two displacement sizes must agree
  before a projected frequency can override a snapshot-Hessian local mode.
- `spectrum_environment_sampled.csv/png`: the uncollapsed sampled environment
  distribution retained for diagnosis even when it fails the display gate.
- `clusters.json` and `clusters/`: selected xTB medoids, geometry-cluster and
  trajectory-occupancy contracts, and their independent resumable ORCA snapshot
  artifacts.
- `spectrum_monomer*.csv/png`: the monomer-only result retained before a
  successful network environment replaces the displayed spectrum.
- `spectrum_dimer.csv`: the selected, per-molecule dimer ensemble.
- `spectrum_network.csv`: the mixed, per-target-molecule dimer/trimer ensemble.

Inspect the adaptive budget without running ORCA:

```bash
storca spectrum --smiles CCO --spectrum-model experimental \
  --phase liquid --fidelity auto --max-orca-jobs 12 --dry-run
```

Use `--fidelity fast` to disable cluster work. `auto` activates the 40-pose xTB
tier only for the eligible condition/interaction contract; `balanced` requests
75 poses. Charged or open-shell clusters are not generated because their
spin coupling and counterion environment are underspecified.

Phase defaults for the residual linewidth are deliberately small: 1, 4, 6,
and 8 cm^-1 for gas, solution, liquid, and solid respectively. Override them
with `--residual-fwhm` when a condition-matched calibration supports doing so.
The model does not add decorative noise or baselines. Restrained finite-time
neutral dimer/trimer trajectories supply conditional snapshot occupancies, not
exact bulk populations. Larger aggregation, explicit solvent configurations, anharmonic
lifetimes, crystal packing, VPT2, and larger-cluster convergence remain
unimplemented. Their absence is recorded rather than hidden by a large
empirical FWHM.

## RMG stability screen

With RMG configured, screen a SMILES string under stated reactor conditions:

```bash
storca stability --smiles CCO --scenario dry-inert-gas-screen --temperature 298 --pressure 1 --rmg-env rmg_env
```

The dark scenarios are homogeneous gas-phase surrogates. The ladder also uses
a low-pressure nitrogen surrogate and a humid-gas composition; none of these is
a liquid, solid, surface, container, or literal-vacuum model.
Every screen records a condition contract (composition, duration, 95% retention
criterion, dark/no-radical-source assumptions, temperature, and pressure).
An unrepaired RMG target-loss profile or crossing time is screening evidence,
not the reported lifetime. A completed model with no threshold crossing is
reported only as no loss in that retained RMG model, not unconditional
stability.

This runs an ORCA optimization/frequency calculation and an RMG-generated
pathway screen. By default, modeled target-loss routes then enter the generic
ORCA/Arkane verification loop described below. Use `--no-auto-verify-routes`
only when a screening-only run is intentional; such a run cannot issue a
verified `t95`. The result is not a safety assessment or a guarantee of storage
stability.

On the M3 laptop, the default RMG execution is deliberately bounded to one
process, 100 generation iterations, 250 edge species, and ten minutes. One
process avoids the RMG 3.3 macOS worker-database failure observed with multiple
processes. RMG lives in the separate `rmg_env` Conda environment and is invoked explicitly
with `--rmg-env rmg_env`; this does not require adding RMG to the active ORCA
environment. Increase these limits only after reviewing the retained RMG logs.

STORCA automatically consults curated RMG reference libraries by elemental
scope: the bundled nitrogen mechanism for nitrogen-containing targets and the
primary H2/O2 mechanism for H/O-only targets. It extracts only direct reactions
reachable from the declared initial species, preserving a small provenance-rich
local subset instead of injecting the complete library. Disable this with
`--no-auto-reference-libraries`, or add another installed database library with
the repeatable `--rmg-database-library NAME` option. Generated products are
bounded to one target-target or target-environment collision by default. For a
controlled comparison, override that cap with
`--rmg-maximum-heavy-atoms NUMBER`; the selected value and its source are
retained in the report.

### Condition ladder

The ordered ladder runs a low-pressure RMG gas surrogate, dark dry nitrogen,
dark dry air, dark humid air, dry-air sunlight, and humid-air sunlight. The
first stage is a dilute target in nitrogen at `1e-6` bar: it approaches
low-collision gas behavior but is not a vacuum calculation, and single-bond
enumeration is not used as its kinetic model. It stops only after a stage
supplies a condition-bound `t95` from converged repaired propagation whose
controlling kinetics passed ORCA/Arkane verification. An initial RMG `t95`, a
photolysis integral, or a verified route without repaired propagation cannot
stop the ladder. Later stages are recorded as not run, not safe:

Cantera reaction-flux attribution is adaptively refined when its integrated
target loss does not close against the propagated target inventory. A route is
still withheld if the bounded refinement budget cannot satisfy the numerical
closure gate.

```bash
storca stability-ladder --smiles CCO --rmg-env rmg_env
```

The isolated parent optimization/frequency result is validated once and reused
across later rungs when SMILES, charge, spin, and method are unchanged. Route
paths, Arkane rates, and repaired RMG propagation remain condition-specific.

Humid air is water vapour at the declared relative humidity (default 50%), not
an aqueous or surface model. Sunlight is declared as the ASTM G173-23 AM1.5
global reference spectrum. A sunlight lifetime is deliberately withheld unless
`--photolysis-evidence` supplies wavelength-resolved absorption, quantum-yield,
and derived photolysis-rate evidence; ground-state ORCA/RMG alone cannot model
photolysis honestly.

`--photolysis-evidence` can either provide a retained,
`photolysis_rate_constant_s-1`, or point to two CSV files:

```json
{
  "spectrum_csv": "am15g.csv",
  "absorption_csv": "target_absorption.csv",
  "container_transmission": 0.90
}
```

The spectrum CSV requires `wavelength_nm,irradiance_W_m2_nm`; the absorption
CSV requires `wavelength_nm,absorption_cross_section_cm2_molecule,quantum_yield`.
STORCA requires identical wavelength grids, calculates the photon-flux integral,
and retains the wavelength contributions used to obtain the photolysis rate.
Photochemical evidence must also declare the resolved product SMILES, for
example `"photoproducts": [{"label": "oh", "smiles": "[OH]"}]`. STORCA
turns the retained first-order photolysis rate into a validated local RMG
library and only uses the coupled RMG target-loss profile to report a sunlight
`t95`; a standalone absorption integral cannot end the ladder.

### Default flux-guided route verification

When a completed RMG model predicts target loss, STORCA ranks reactions by
quantitative target-destruction flux. If a controlling reaction consumes a
generated intermediate, the reachable upstream formation route is verified
first. One route is handled per iteration; after every repair, flux is
recomputed instead of assuming the original ranking remains valid.

The route verifier preserves the exact RMG species labels, structures,
stoichiometry, charges, fragment multiplicities, and reaction equation. It
resolves a common total-spin surface and a complete atom correspondence,
failing closed when mapping or electronic state is ambiguous. Multiple
deterministic endpoint-complex orientations are used for multi-fragment routes.
Each separated species is optimized and frequency-validated independently; the
assembled encounter geometries are NEB seeds and are not required to be
fictitious bound minima. For an activated route, ORCA runs NEB-TS and an
explicit TS refinement, requires exactly one significant imaginary mode, and
requires a bidirectional IRC to match the declared bound endpoints or separated
fragment channels. The same IRC geometry sequence can be rendered as the
decomposition animation.

The declared ORCA charge and multiplicity constrain the combined system; they
do not by themselves preserve fragment-local radical multiplicities after
fragments are assembled. If a validated intermediate splits a route into an
association or dissociation segment, STORCA therefore retains a
fragment-spin-aware capture/recrossing or constrained-scan requirement instead
of applying an unconstrained fixed-end NEB-TS calculation. A persistent NEB
whose highest-energy image remains an endpoint is stopped early and retained as
an unresolved path-quality result, never as a transition state or rate.

If no saddle point is retained, a barrierless label requires a normally
completed path with complete energies, endpoint agreement, and the same result
from two independent orientations. That classification does not itself provide
a rate. The default Arkane adapter currently accepts only the validated
TS/frequency/IRC branch; a barrierless candidate remains incomplete until a
collision-bounded capture/recrossing rate model is available.

Arkane consumes separately validated reactant, product, and transition-state
ORCA frequency artifacts and calculates condition-containing temperature and
pressure grids. A validated rate is written to a route-specific local RMG
library, checked against the retained collision limit where applicable, and
then reinjected into a **full RMG model rerun** under the exact original
condition contract. STORCA verifies that the replacement was actually loaded
at the intended rate, reranks the repaired flux, and repeats until the `t95` is
robust to remaining routes or verification fails closed. It does not obtain a
final lifetime by merely rescaling the original profile.

### Comparing RMG releases

RMG 4 is installed separately as `rmg4_env` (RMG-Py 4.0.0); the existing
`rmg_env` remains the RMG-Py 3.3.0 baseline. Compare them with matched
conditions and resource limits rather than replacing the baseline blindly:

```bash
storca rmg-compare --smiles OO --rmg3-env rmg_env --rmg4-env rmg4_env \
  --scenario ambient-air-gas-screen --target-duration-hours 0.01
```

The retained `rmg-comparison.json` records exact versions, wall time, model
size, target loss, and the full Chemkin/solver artifacts for each run. It is a
comparison of two bounded models, not an accuracy ranking; inspect any
mechanism difference before adopting one version as the new production default.

### TD-DFT photostability screen

Use ORCA TD-DFT to screen for bright vertical transitions in the default
280–800 nm solar window:

```bash
storca photostability-screen --smiles OOO --nroots 20
```

The command verifies the single ground-state geometry as a local minimum,
requests full linear-response TD-DFT (not the Tamm-Dancoff approximation), and
retains the geometry, ORCA input/output, electric-length-gauge transition
table, and a broadened **relative** absorption CSV. A negative bright-state
result is withheld when the requested roots do not reach the 280 nm edge.
Supply a
declared irradiance CSV (`wavelength_nm,irradiance_W_m2_nm`) with
`--sunlight-spectrum` to retain the corresponding relative solar-overlap
integral. This is an absorption-risk screen only: it cannot report a
photolysis rate or `t95` unless quantum yield and resolved product evidence
are subsequently supplied to the photolysis/RMG workflow.

For a bounded model projection that creates RMG photon reactions, provide a
declared source spectrum:

```bash
storca computational-light --smiles OO --sunlight-spectrum source.csv \
  --rmg-env rmg4_env --scenario ambient-air-gas-screen
```

This uses one full LR-TDDFT calculation and at most three generic homolysis
product-energy calculations, then propagates low (0.1%), central (1%), and
high (10%) reactive-photon branch profiles through RMG. Route quantum yields
are wavelength-resolved and smoothly energy-gated; energetically inaccessible
routes no longer consume the full profile prior. Its `t95` values are a
model sensitivity range—not measured quantum-yield lifetimes. Every computed
candidate route receives a labelled PNG storyboard; the same visual module can
render a future ORCA IRC trajectory as a GIF.

`storca/arkane_runner.py` creates pressure-dependent Arkane jobs from verified
ORCA stationary-point and transition-state artifacts. It evaluates a declared
temperature/pressure grid and bath gas; it does not infer a transition state or
report a rate from fragment topology alone.

### Decomposition explanations and animations

Create a route explanation from a completed ladder, stability screen, or
computational-light result:

```bash
storca stability-explain runs/<run>
```

The command selects the leading executable route (earliest modeled `t95`,
kinetic relevance, and light-route accessibility are considered) or accepts
`--route <id>`. Use `--prepare-only` to create mapped endpoints and a visibly
labelled geometry preview without running an ORCA reaction path.

The current explanation layer uses two evidence levels:

- A single retained bond homolysis receives an ORCA relaxed bond-distance scan.
  The result is labelled `computed_dissociation_path` only when the requested
  scan completes and ORCA terminates normally.
- A resolved unimolecular or multi-fragment route can enter the same generic
  endpoint-complex, NEB-TS, TS-frequency, and IRC verifier used by stability.
  Atom mapping may come from retained map numbers or the bounded graph-edit
  resolver; ambiguous mappings are rejected. A barrierless path is labelled
  only after complete, endpoint-matched profiles agree across two orientations,
  and remains rate-unverified.

Legacy RMG routes containing Chemkin labels but no resolved endpoint structures
return `unsupported_route_physics`. Resolved structures with an ambiguous atom
correspondence fail closed rather than receiving an arbitrary path. Each run
retains `explanation/decomposition-explanation.json`, the ORCA inputs and
outputs, a candidate storyboard, and—when trajectory frames exist—a GIF and
energy profile. `incomplete_dissociation_path` animations are useful diagnostic
paths but are prominently labelled incomplete and do not upgrade stability or
kinetic evidence.

A ground-state dissociation scan is not an excited-state trajectory. In
particular, computational-light route animations show the selected chemical
coordinate on the declared ground-state spin surface, not fluorescence,
internal conversion, intersystem crossing, or nonadiabatic photodynamics.

## Persistence benchmark foundation

`benchmarks/plausibility/` contains a condition-specific persistence benchmark
schema. Its initial entries are explicitly draft candidates: they cannot tune
or validate STORCA until a condition-matched evidence note and source are
reviewed and the entry is marked `accepted`. Evaluate completed dossiers with:

```bash
storca plausibility-benchmark benchmarks/plausibility/manifest.json dossiers/
```

The key metric is false reassurance: reporting an accepted reactive or
transient reference as ordinary-condition persistent.
