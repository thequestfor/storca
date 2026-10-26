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
The stability and chemical-fate program has its own
[scientific roadmap](STABILITY_ROADMAP.md),
[requirement catalog](docs/stability/REQUIREMENTS.md),
[traceability matrix](docs/stability/TRACEABILITY.md), and
[definition of a supported lane](docs/stability/DEFINITION_OF_SUPPORTED.md).
The machine-readable roadmap manifest is retained at
[`benchmarks/stability/milestones.json`](benchmarks/stability/milestones.json)
and is checked against the runtime lane registry by the test suite.

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

The current implementation status and dependency-ordered expansion plan are
tracked in [STABILITY_ROADMAP.md](STABILITY_ROADMAP.md). Registration, software
implementation, scientific validation, and observable support are deliberately
different states; adding a program or lane declaration does not promote a
scientific capability.

Prepare a multi-environment campaign without launching chemistry calculations:

```bash
storca stability-campaign \
  --request examples/stability-campaign/gas-and-solution.json \
  --campaign-dir runs/campaign-preparation
```

Schema version 1 validates a shared molecular system, lane-specific
environments, observables, methods, evidence references, applicability flags,
and compute ceilings. It writes `campaign-request.json`, `campaign-plan.json`,
individual `lane-plans/`, and `campaign-summary.json`. Reusing the directory
with changed scientific inputs fails; display metadata may be retained without
changing the scientific hash. Preparation is the default and emits no
observable. The gas-only, solution-only, combined, and deliberately unsupported
examples under `examples/stability-campaign/` demonstrate ready, missing-input,
and declared-only states.

Live campaign execution supports admitted homogeneous-gas `t95` lanes and two
restricted neutral implicit-solution `t95` lanes. Gas uses the existing bounded
RMG → ORCA route verification → Arkane → repaired-network propagation
machinery. Solution execution propagates a condition-matched verified
first-order network, runs the existing observable-conditioned solution
verification session, or propagates a separately evidenced first/second-order
ideal-solution mass-action network. Execution requires two independent controls:
`execution_policy` must be `explicit_execute` in the immutable request and the
command must include `--execute`.

This synthetic verified-input solution smoke test runs without RMG or ORCA:

```bash
storca stability-campaign \
  --request examples/stability-campaign/solution-verified-executable.json \
  --campaign-dir runs/solution-campaign \
  --execute
```

Its rate and coverage evidence are deliberately synthetic software-regression
fixtures, not a claim about experimental cyclobutane stability. The nonlinear
mass-action smoke test is:

```bash
storca stability-campaign \
  --request examples/stability-campaign/solution-mass-action-executable.json \
  --campaign-dir runs/solution-mass-action-campaign \
  --execute
```

This second synthetic smoke test exercises dimensional bimolecular rates,
diffusion ceilings, nonlinear propagation, and a declared reactive partner. It
does not infer solution rates or establish experimental mechanism coverage.
A live gas example is:

```bash
storca stability-campaign \
  --request examples/stability-campaign/gas-executable.json \
  --campaign-dir runs/gas-campaign \
  --execute --cores 4 --rmg-env rmg_env
```

The gas command is real chemistry work and can be expensive. The mixed
`gas-and-solution-executable.json` example is also live and will launch that gas
work; it is intended for an installed RMG/ORCA/Arkane setup or an injected test
runner, not a cheap smoke test. Only explicitly enforceable per-lane and global
ceilings are accepted. Exact-contract reruns reuse completed gas stages and
completed solution lane results; screening `t95` values cannot terminate the
gas ladder. Other environment lanes remain retained as blocked or
adapter-not-implemented results, and their physics is never replaced with a gas
or implicit-solution calculation. Campaign reports do not combine incomparable
lane values into a universal lifetime.

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

### Optional independent graph/xTB discovery

The default discovery execution remains RMG-only unless xTB is requested, but
it no longer bypasses mechanism-coverage accounting: a completed RMG model is
retained as useful evidence while the final lifetime or no-loss claim is
withheld as `verification_incomplete` because independent discovery is absent.
Add
`--graph-xtb-discovery` to run a second, deliberately bounded discovery lane
alongside the retained RMG model:

```bash
storca stability --smiles CC --scenario dry-inert-gas-screen --graph-xtb-discovery
```

This lane edits the target molecular graph (single-bond breaking, formation,
one-break/one-form rearrangement, and mapped hydrogen transfer). When the gas
condition declares reactive partners with charge and multiplicity, it also
performs bounded pairwise radical association and H-atom abstraction. It rejects
invalid, spin-incompatible, or duplicate products, then uses GFN2-xTB optimised **endpoint
electronic energies** only to screen grossly endergonic proposals. Its retained
artifact records the edit, atom maps, structural pruning decision, calculation
contract, energies, hashes, xTB failures, and any search truncation. It does
not predict a barrier, a rate, or kinetic relevance.

RMG and graph/xTB candidates are merged by net molecular-graph stoichiometry
while preserving both discovery assertions. Agreement is independent supporting
evidence; disagreement is an explicit mechanism-coverage question. A
graph-only survivor is ranked only by its retained xTB reaction-energy screen,
is marked `kinetic_relevance: unknown`, and is sent through the same ORCA
endpoint/NEB-TS/frequency/IRC, Arkane, and full-RMG-repropagation gates. It
cannot create a `t95` by itself. If its path or rate cannot be resolved, the
result remains `verification_incomplete` rather than reverting to an RMG-only
no-loss conclusion.

When graph/xTB is requested, RMG production, graph execution or validated
reuse, source-preserving merge, and the supplemental verification queue are
bound by one immutable discovery-assembly contract. The retained assembly
separates source-agreement, RMG-only, graph-only, screened-out, malformed, and
unresolved assertions. Exact terminal reruns reuse the assembly. A truncated or
unavailable requested backend, a failed xTB endpoint screen, or an RMG route
without resolvable endpoints vetoes an otherwise converged lifetime; backend
completion still describes only the declared bounded search.

Every gas stability run now writes
`reaction-discovery/mechanism-coverage.json`. This condition- and
observable-bound artifact records the explicitly searched reaction classes,
candidate source agreement and disagreement, graph-only candidates, screened exclusions,
backend failures or truncation, applicability gaps, and each candidate's
current influence state. It also retains an uncertainty register that keeps
mechanism coverage separate from electronic-structure, conformer, TS,
environment, rate-theory, tunneling, pressure, electronic-state, propagation,
and numerical uncertainty. Content-addressed revisions preserve both the
pre-verification and post-verification decisions. A graph-only survivor remains
`potentially_influential_unknown_rate` until the ordinary ORCA path, rate,
repaired-network, and observable-robustness gates resolve it. Candidates left
outside a bounded high-level verification queue therefore block an observable
claim instead of disappearing from convergence accounting.

Coverage decisions are categorical. `adequate_within_declared_scope` means
only that the requested observable passed every retained gate within the exact
declared searches, conditions, and methods. It is not a probability of
mechanism completeness. Search-budget stability, when supplied, is likewise
reported only as stability within that declared expansion. The artifact always
retains `mechanism_completeness_claimed: false`.
Two completed backends that both return an empty candidate set still produce
`empty_search_shared_blind_spot_unresolved` until an applicable calibrated
benchmark, bounded expansion, or complementary discovery source addresses the
shared-blind-spot risk. Reactive co-species and open-shell targets likewise
remain blocked when the independent search did not cover bimolecular or radical
chemistry.

Version 2 remains deliberately narrow and gas-workflow-compatible. Its
bimolecular scope is pairwise radical association and H abstraction only; it
does not supply encounter/capture rates or enumerate closed-shell substitution,
solvent-assisted chemistry, termolecular chemistry, surfaces, photochemistry,
arbitrary multi-bond edits, or reaction dynamics. Use
`--graph-xtb-max-candidates` and `--graph-xtb-energy-window` to make its
high-level queue and endpoint-energy screen explicit.

An incomplete retained gas verification can be continued without redeclaring
scientific options:

```bash
storca stability-resume-verification runs/RETAINED-STABILITY-RUN
```

Fresh stability runs retain a hashed verification execution contract. The
continuation command accepts only the run directory, rejects changed
conditions/methods/library content, creates a new engine session, and reuses
the route calculation workspace. A retained v2 terminal deferral is migrated
to the current bounded surface-resolution policy through a separate append-only,
hashed scientific-policy amendment; the immutable execution contract is not
rewritten, and the amendment alone permits no mechanism, rate, or lifetime
claim. A scientifically deferred route is retried only when the explicit
path-verifier policy changed; repeating the same policy fails closed. An
interrupted `started` route is different: STORCA follows the
append-only continuation chain, checks retained calculation activity, and
refuses a duplicate while artifacts are fresh or a persisted ORCA process lease
is live. After a quiet, dead execution it reuses the same policy attempt before
any scientific-policy migration. For an
interrupted NEB/NEB-TS calculation, a complete atom-order-consistent
`*_MEP.allxyz` band is validated and supplied to ORCA through
`Restart_ALLXYZFile` in a new, hashed restart input; the partial output is never
overwritten. Recognized ORCA/OpenMPI startup denials likewise retain the failed
output and launch an exact-input `-execution-retry-###` artifact under the same
scientific-policy attempt; they are never reported as chemical scan failures.
This is route- and NEB-level continuation, not yet a general
restart from every Arkane or repaired-RMG interruption point.

### Condition ladder

The ordered ladder runs a low-pressure RMG gas surrogate, dark dry nitrogen,
dark dry air, dark humid air, dry-air sunlight, and humid-air sunlight. The
first stage is a dilute target in nitrogen at `1e-6` bar: it approaches
low-collision gas behavior but is not a vacuum calculation, and single-bond
enumeration is not used as its kinetic model. By default it stops only after a
stage supplies a condition-bound `t95` from converged repaired propagation whose
controlling kinetics passed ORCA/Arkane verification. An initial RMG `t95`, a
photolysis integral, or a verified route without repaired propagation cannot
stop the `worst_case_search` campaign. Later stages are then recorded as not
run, not safe:

Cantera reaction-flux attribution is adaptively refined when its integrated
target loss does not close against the propagated target inventory. A route is
still withheld if the bounded refinement budget cannot satisfy the numerical
closure gate.

```bash
storca stability-ladder --smiles CCO --rmg-env rmg_env
```

Use `--campaign-mode environment_map` when the deliverable is the requested
condition map rather than the first verified loss condition. It executes every
gas rung and retains the executed/skipped stage lists, every verified `t95`,
and stage-specific evidence gaps in `stability-ladder.json`.

```bash
storca stability-ladder --smiles CCO --rmg-env rmg_env \
  --campaign-mode environment_map --graph-xtb-discovery
```

With `--graph-xtb-discovery`, the ladder performs a bounded graph/xTB search
inside each gas condition. Dry inert rungs search target-unimolecular edits;
air and humid-air rungs additionally search the reactive partners actually
declared in that rung. This prevents target-only evidence from being mistaken
for independent O2/water reaction coverage. Every condition retains its own
assembly, endpoint screens, candidates, NEB/TS/frequency/IRC evidence, Arkane
rates, repaired libraries, and propagation. Missing partner electronic states,
backend failure, or truncation remain explicit evidence gaps.

### Stability lane contracts

Every default ladder rung now carries an immutable lane definition,
environment identity, runtime admission decision, and standardized lane-result
summary. The lane contract identifies the discovery sources, path verifier,
rate model, propagator, convergence policy, supported observables, required
capabilities, applicability boundary, and support level. It wraps the existing
condition and evidence artifacts; it does not replace or relax them.

Inspect the registry with:

```bash
storca stability-lanes
storca stability-lanes --environment-kind homogeneous_solution --json
```

The current gas and illuminated-gas lanes are `observable_supported`, still
conditional on their existing route, rate, propagation, and convergence gates.
Two homogeneous-solution lanes are executable, each in a deliberately narrow
verified-input domain described below. Microsolvation, surface, and
electrochemical entries remain `declared_only`. Moving a lane to
`observable_supported` requires an implemented propagation chain and a
lane-specific validation profile, not merely an installed program or keyword.

### Verified-input implicit-solution lane

`solution-neutral-implicit-v1` supports dark, neutral, homogeneous,
constant-volume solution kinetics only when every retained reaction is an
element/charge-balanced first-order elementary step with an external verified
solution-rate claim or a STORCA-generated solvent-consistent path/rate claim.
The evidence bundle must also carry a restricted
reaction-network coverage claim with no unresolved candidates. The propagator
does not convert gas rates or infer solvation corrections.

STORCA can produce one such rate claim for an ordinary neutral unimolecular
candidate with its existing endpoint/NEB-TS/TS-frequency/IRC engine.
Preparation is the default; add `--execute` only after inspecting the retained
request and endpoint mapping:

```bash
storca solution-path-verify \
  --request ring-opening-smd.json \
  --method-profile b3lyp-def2-svp

storca solution-path-verify \
  --request ring-opening-smd.json \
  --method-profile b3lyp-def2-svp \
  --execute --cores 4
```

The request contains the same environment and network fields shown below,
plus one `candidate_reaction` and one executable `route`. The candidate must
declare a positive composite `uncertainty.barrier_kcal_mol`, its source, and
components covering `electronic_structure_barrier`, `conformer_sampling`,
`implicit_solvent_model`, and `harmonic_tst`. Network coverage remains an
independent claim: verifying one path cannot manufacture mechanism
completeness.

Every generated ORCA input is audited for exactly one matching `SMD(solvent)`
or `CPCM(solvent)` keyword. A verified bundle is emitted only after endpoint
minima, a single-imaginary-mode saddle, bidirectional IRC endpoint matching,
condition-matched thermochemistry, and a positive activation free energy pass.
The producer explicitly requests the environment temperature through ORCA's
single-temperature `%freq Temp` control, including endpoint, recovery and TS
frequency jobs. Input, output and cache temperatures must agree before rate
generation; changing temperature requires a new request workspace. No extra
temperature flag is needed. The producer also observes ORCA 6.1's actual
RRHO/quasi-RRHO entropy recipe, reference/cutoff frequencies, applied frequency
scale and thermochemistry pressure. Reactant and TS recipes must match; missing
or contradictory evidence withholds the rate. No default is inferred from a
method-profile name. This remains conventional TST with transmission
coefficient one and a zero 1 M correction for the unimolecular
reactant-to-TS step. Named CPCM/SMD parameters remain fixed: this is not a
validated arbitrary-temperature solvent model or a new scientific support
claim. Barrierless results are retained for future capture/VTST escalation and
are never forced through this rate model.

Quasi-RRHO entropy is not an explicit hindered-rotor or anharmonic calculation.
The observed recipe is recorded separately from the kinetic theory; gas Arkane
partition functions and legacy supplied-result bundles are not rewritten.

The solution rate producer now automatically samples reactant conformers when
rotatable bonds or nonaromatic rings of at least four atoms trigger the bounded
protocol. Its first production scope is one neutral, closed-shell singlet with
ordinary isotope treatment and at most 20 heavy atoms; known unsupported
reactants are refused before path calculations. RDKit ETKDGv3/UFF generates
geometric seeds, not populations or rate estimates. Up to four selected seeds
are refined with the same solution ORCA method, temperature and shared route
budget. Identity, complete frequency evidence, free-energy recipe and immutable
calculation sources must agree.

Distinct refined minima do **not** automatically become Boltzmann-weighted
network rates. Without interconversion and TS-conformer evidence, the producer
returns `conformer_rate_model_unresolved`, retains a clearly conditional
single-path equilibrium scenario where numerically possible, and emits no
aggregate rate/bundle. Failed, truncated or incomplete sampling also withholds
admission. One refined cluster still requires the separately justified model/
sampling bound; finite sampling does not prove completeness. No extra command
flag is needed. Completed refinements can be revalidated without new jobs after
the shared budget expires. Old session-produced TST rates cannot silently gain
new temperature/recipe/conformer-policy credit; explicitly supplied networks
remain labeled imports. This implementation is software-tested, not yet a
real conformer-rate accuracy demonstration.

### Verified-input solution mass-action lane

`solution-neutral-implicit-mass-action-v1` accepts a distinct evidence schema
for dark, neutral, constant-volume CPCM/SMD solution with ideal molar
activities. It supports first-order (`s^-1`) and second-order
(`L mol^-1 s^-1`) elementary rates, reversible pairs, and target reformation.
Every rate needs explicit uncertainty and SHA-256-linked evidence. Every
second-order rate also needs a sourced diffusion ceiling that bounds its upper
rate interval. Gas and RMG rate contracts are rejected rather than converted.

The nonlinear solver retains and can claim concentrations, reaction extents,
product yields, dominant product, first t95 crossing, retention, gross
target-consumption branching, rate-bound scenarios, and numerical closure. A
single lane request may name several observables. A campaign must declare every
positive-concentration starting species
as a chemical-system component. The lane does not perform solution discovery,
speciation, pH/activity correction, explicit solvent chemistry, or transport;
those cases remain outside applicability.

### Observable-conditioned solution verification session

Several retained native, resolved RMG, graph/xTB, or user candidates can be
allocated through one append-only session:

```bash
storca solution-verify-session \
  --request solution-session.json \
  --session-dir runs/solution-session

storca solution-verify-session \
  --request solution-session.json \
  --session-dir runs/solution-session \
  --execute --cores 4
```

The immutable request declares `environment`, `target_smiles`, `network_context`,
`candidate_reactions`, `observable: t95`, `method_profile`, a convergence
contract, and both ORCA-path and estimated-core-hour budgets. Reusing the exact
session directory resumes only the identical request; a changed request is
rejected instead of mutating prior evidence.

Speculative candidates live in a `SolutionCandidateNetwork`, separately from
the strict `VerifiedSolutionNetwork`. Gas-phase RMG rates and graph/xTB evidence
may rank work but cannot become solution-rate evidence. For the current
monotonic first-order lane, t95 sensitivity is analytical: only reactions that
directly consume the target affect target disappearance. Downstream candidate
reactions therefore receive zero t95 allocation priority.

Each executed iteration retains observable bounds, the complete score
breakdown, one selected candidate, its cache/environment contract, path result,
and the repaired network. Screening t95 values remain allocation evidence and
are withheld from the scientific claim. A final t95 requires verified direct
target-loss rates, robust horizon classification, the requested interval
tolerance, and an independent restricted network-coverage claim. Specialized
physics, unresolved paths, coverage gaps, and exhausted budgets are terminal
explicit states. A general concentration/product bundle is withheld when
unverified downstream candidates remain, even if the t95-only claim has
converged.

The session uses the same explicit requested-temperature contract for neutral
unimolecular chemistry. It does not add conformer completeness, temperature-
dependent solvent parameters, bimolecular rate production, diffusion,
speciation, microsolvation, photochemical, surface, or electrochemical physics.

Run one bundle directly with:

```bash
storca solution-stability \
  --evidence-bundle verified-solution-network.json
```

The bundle has this top-level form (the skeletal fragment below is deliberately
incomplete and cannot produce a t95):

```json
{
  "schema_version": 1,
  "kind": "verified_implicit_solution_stability_bundle",
  "target_smiles": "<required exact target identity>",
  "environment": {
    "solvent": "water",
    "solvent_model": "SMD",
    "temperature_K": 298.15,
    "target_duration_seconds": 86400.0,
    "retention_fraction": 0.95,
    "standard_state_mol_L": 1.0,
    "light_condition": "dark",
    "electronic_structure_profile": "<retained profile identity>"
  },
  "network": {
    "target_species_id": "target",
    "coverage_status": "verified_complete_within_declared_unimolecular_scope",
    "coverage_claim_id": "<evidence claim ID>",
    "discovery_sources": ["<retained discovery or curated-mechanism source>"],
    "unresolved_candidate_ids": [],
    "species": [{
      "species_id": "target",
      "smiles": "<canonicalizable molecular identity>",
      "initial_concentration_mol_L": 0.01,
      "element_counts": {"C": 1, "H": 4},
      "charge": 0
    }],
    "reactions": []
  }
}
```

Every species SMILES is canonicalized and its RDKit-derived elemental formula
and formal charge must match the declared bookkeeping. Each reaction record
uses one `reactant_id`, integer product stoichiometry, a
positive `rate_constant_s`, exact solvent/model/temperature/one-molar standard
state, a source, and a `rate_claim_id`. Every target-consuming rate also needs
explicit lower/upper multipliers bracketing 1.0 and an `uncertainty_source`.
A `verified_solution_tst` rate also
requires `path_claim_id`, `standard_state_correction: one_molar_solution`, and
calculation hashes. Experimental first-order solution rates use
`verified_experimental_solution_rate` and `experimental_solution_rate`.

The native concentration propagator evaluates an augmented linear system with
a scaling-and-squaring matrix exponential. Species concentrations and reaction
extents are solved together, providing t95, product concentrations, branching,
dominant target-loss flux, elemental closure, non-negativity, and flux closure
from one retained calculation. Combined target-rate bounds produce a retained
t95 interval. A within-horizon t95 is withheld when the qualitative threshold
crossing changes across those bounds; the interval is a physical allocation
bound, not a statistical confidence interval.

Attach one or more bundles to the main campaign with:

```bash
storca stability-ladder --smiles CCO --rmg-env rmg_env \
  --campaign-mode environment_map \
  --solution-evidence-bundle ethanol-water-solution.json
```

The bundle `target_smiles` must match the campaign target. Additional lanes are
reported separately from gas stages because their physical models are not
interchangeable. Under `worst_case_search`, solution bundles after the first
verified t95 are retained as not run.

The solution lane fails closed for bimolecular encounter chemistry,
diffusion limitation, explicit solvent participation, acid/base speciation,
ion pairing, ionic species, photochemistry, spin crossing, surfaces,
unresolved multireference behavior, target-reforming networks, condition-rate
mismatches, or incomplete mechanism coverage. Those findings produce no t95.

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

An optional evidence file can make this allocation observable-conditioned:

```bash
storca stability --smiles CCO --scenario dry-air-gas-screen \
  --rate-uncertainty-evidence rate-uncertainty.json
```

```json
{
  "schema_version": 1,
  "kind": "rate_uncertainty_evidence",
  "records": [{
    "reaction_equation": "stability + oxygen <=> P + radical",
    "lower_rate_multiplier": 0.1,
    "upper_rate_multiplier": 10.0,
    "source": "named calibration dataset or user analysis",
    "source_kind": "method_calibration",
    "applicability": "declared reaction class, state, environment, and temperature range",
    "uncertainty_components": ["electronic_structure", "transition_state_search"]
  }]
}
```

STORCA combines an explicit multiplier bound with retained logarithmic `t95`
sensitivity only when every selectable unresolved nonzero-flux route has both
pieces of evidence. Incomplete coverage falls back to the established
sensitivity/flux ordering and lists the missing route IDs. Optional verification
costs affect ranking only when all required records share explicit units and the
same cost model. This feature prioritizes which existing ORCA/Arkane verification
to run next; it does not propagate a lifetime confidence interval and does not
claim expected information gain. Upstream reachability dependencies still take
precedence over the score.

If uncertainty coverage is complete, missing local sensitivities are calculated
automatically from independent 0.5x/2x rate perturbations of the exact retained
Cantera reaction indices. Both runs must complete and retain a `t95`. The
diagnostic is repeated after mechanism repair because network sensitivities can
move; a failed diagnostic is recorded and falls back to the prior ordering.

The same file may be passed to `storca stability-ladder`. Each record is matched
again inside each condition. Set `condition_contract_id` to restrict evidence to
one exact rung; omit it only when the stated applicability defensibly spans the
conditions. Unmatched or condition-mismatched records remain visible evidence
and never force a route match.

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

For an unresolved path, the v3 surface-resolution policy first diagnoses the
retained band—including energy/geometry, per-image connectivity, mapped bond
changes, and separated-fragment motion—and then spends a single route-wide
budget on bounded stationary recovery, mapped one/two-coordinate TS seeding,
or validated-intermediate segmentation. A resolved intermediate becomes two source-neutral native
reaction candidates; each bound segment must independently pass TS, frequency,
and IRC checks. Arkane then builds two elementary rates and STORCA inserts both
with the intermediate before rerunning the full network. If identity, segment
physics, rate theory, budget, or repaired-network insertion remains unresolved,
the report retains the precise state and keeps `t95` null.

Every active surface strategy now emits one common append-only attempt-evidence
record linking the reaction/environment/mapping contract, exact inputs, method
and electronic state, policy and parameters, process leases and budget usage,
outcome/failure, and the claims it does or does not support. Implicit solution
is accepted only with a declared CPCM/SMD solvent and the exact matching ORCA
keyword; unsupported or mismatched environments stop before ORCA as
`environment_model_escalation_required`.

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

When a verified path instead stops at `dissociation_rate_model_required` or
`barrierless_rate_model_required`, export the condition- and path-bound request:

```bash
storca stability-export-specialist-request runs/<run>
```

After an applicable VTST, capture, master-equation, or other allowed specialist
method has produced the schema-v1 result, raw evidence, uncertainty statement,
and checksummed RMG library, continue without repeating ORCA:

```bash
storca stability-resume-verification runs/<run> \
  --specialist-rate-result path/to/result.json
```

The option also accepts a directory containing results for multiple request
hashes. STORCA matches exactly one result per request and fails closed on a
condition mismatch, missing applicability evidence, mutated artifact, or
mutated kinetics-library file. This imports specialist evidence; it is not a
built-in capture/VTST solver.

Provider developers can exercise the separate capture/VTST/master-equation
conformance and numerical gate with:

```bash
storca specialist-rate-benchmark \
  benchmarks/kinetics/specialist_rate_manifest.json \
  path/to/provider-results
```

The bundled cases are a development set and deliberately cannot authorize a
support claim without a locked holdout. See
[`docs/stability/SPECIALIST_RATE_BACKENDS.md`](docs/stability/SPECIALIST_RATE_BACKENDS.md)
for the backend decision record and required retained evidence.

A source-backed one-well MESS model can be rendered without executing it:

```bash
storca stability-prepare-mess-input model.json prepared-mess/
```

The resulting manifest is explicitly `prepared_not_executed` and prohibits
rate and lifetime claims. `storca doctor --rmg-env rmg4_env` checks the same
RMG environment used by stability runs instead of incorrectly relying only on
the active shell's `PATH`.

If MESS is installed on `PATH`, run the prepared input with:

```bash
storca stability-run-mess-input prepared-mess/
```

This command is bounded, checkpointed, and exact-channel parsed. Its result is
`completed_unadmitted_rate_table`, not a verified rate. If the required unit,
uncertainty, collision-limit, and independent thermodynamic evidence has been
assembled, finalize the exact exported specialist request with:

```bash
storca stability-finalize-mess-rate prepared-mess/ \
  mess-admission-evidence.json specialist-rate-request.json \
  --rmg-env rmg4_env
```

Finalization parses the reverse channel, checks detailed balance and collision
ceilings over the complete grid, converts only an explicitly verified amount
basis, fits a pressure-dependent RMG model, reloads the saved library, and
checks every grid point. It emits `specialist-rate-result.json`, which can be
passed to `stability-resume-verification`. See
[`docs/stability/MESS_RATE_ADMISSION.md`](docs/stability/MESS_RATE_ADMISSION.md)
for the evidence contract. This does not make a lifetime claim before the full
network is repaired and repropagated, and the development benchmark still does
not establish general MESS accuracy.

For an ordinary verified 1↔2 transition-state route, the Arkane adapter can
also consume a route-level `storca_explicit_gas_pressure_dependence_model`.
That evidence must declare the applicable T/P range, sourced
single-exponential-down parameters for the energized well, sourced
Lennard-Jones parameters for the well and all bath gases, and immutable
artifact hashes. Without it, Arkane's generic collisional estimates are
retained as an unvalidated central screening model. In both cases, a nonzero-
flux rate needs uncertainty covering `total_rate_model` or all named model
components before the stability observable can converge.

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
