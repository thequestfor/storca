# STORCA

STORCA is a small, reproducible workflow for ORCA geometry optimizations and
vibrational frequency calculations. Each calculation writes to its own run
directory rather than cluttering the repository root.

## Quick start

Install the package in editable mode, then check prerequisites:

```bash
python -m pip install -e .
storca doctor
```

`doctor` checks whether ORCA, xTB, Open Babel, RMG, RDKit, and the Python Open
Babel binding are discoverable in the active environment.

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
contract.

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

## Practical IR display

The normal `spectrum` command preserves the raw, conformer-weighted harmonic
calculation. Add `--spectrum-model practical` to create a second, explicitly
rule-calibrated display for common organic experimental IR spectra:

```bash
storca spectrum --smiles 'O1CC(O)OCC1O' --spectrum-model practical --cores 4
```

The calculated spectrum remains in `spectrum_raw.csv` and `spectrum_raw.png`.
The displayed `spectrum.csv` records the `ambient-organic-v0.2` rules in
`metadata.json`. These rules broaden and, for hydrogen-bond-capable O-H/N-H
regions, shift calculated modes into a practical condensed-organic envelope.
At present, this profile has retained only its flexible-molecule fingerprint
rule. The initial broad O-H envelope helped a KBr-disc benchmark but degraded
an ethanol reference, so it is not applied without an explicit sample-aware
model. These rules are an empirical display calibration—not a solvent, crystal, or
first-principles linewidth calculation—and will be retained only when they
improve held-out benchmarks.

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
