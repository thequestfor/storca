# Limitations

- Results depend on the selected electronic-structure method, basis set, and
  input structure; the defaults are not appropriate for every molecule.
- A frequency calculation with no large imaginary modes is evidence for a
  local minimum, not a guarantee of global stability or experimental isolation.
- IR conformer populations currently use optimized electronic energies rather
  than quasi-harmonic conformational Gibbs free energies. Duplicate minima that
  differ only through symmetry-equivalent atom mappings are not yet guaranteed
  to be collapsed.
- The current dimer/trimer environment branch is a bounded static sample, not a
  finite-temperature liquid ensemble. Its sufficiency gate prevents an
  undersampled frequency spread from being displayed as a calculated width,
  but passing the minimum gate does not supply trajectory populations, extended
  hydrogen-bond networks, anharmonic lifetimes, or benchmark validation. The
  restrained GFN2-xTB poses deliberately cover interaction coordinates and
  therefore do not have thermodynamic populations. Selected medoids receive
  equal geometry-stratum weights, not liquid occupancies. Their ORCA Hessians
  are evaluated at restrained-xTB snapshots rather than unconstrained DFT
  stationary points. Adaptive batch convergence tests numerical stability of
  this bounded representative set; it does not establish bulk thermodynamic
  convergence or repair its equal-stratum population approximation.
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
- A numerically stable restrained static ensemble can still have broad
  bootstrap uncertainty. Static extension is capped and reported as
  statistically imprecise when that occurs; it does not manufacture precision
  or convert stratified equal weights into liquid occupancies.
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
- Automatic ORCA/Arkane verification is default for modelled target-loss
  routes, but it can finish only when RMG retains resolved endpoint structures,
  charges, multiplicities, reaction identity, and quantitative flux. Ambiguous
  atom mappings, multiple unresolved spin surfaces, large mapping searches, or
  missing upstream routes fail closed rather than being guessed.
- Generic NEB-TS searches sample a finite set of generated endpoint complexes
  on one declared adiabatic electronic surface. Separated species are validated
  as individual minima, but their finite-separation encounter geometries are
  path seeds rather than stationary complexes. The search can miss a
  lower-energy conformer, encounter orientation, intermediate, or transition
  state. Spin crossings, electronically nonadiabatic paths, and strongly
  multireference bond breaking are not resolved by this workflow.
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
- Flux-guided verification is bounded by configured iterations and per-job
  timeouts. Complex networks may require more controlling or upstream routes
  than the budget permits and will finish as `verification_incomplete`.
- Even `verified_t95_converged` is conditional on the declared gas model,
  reaction-network coverage, electronic-structure method, Arkane model, and
  retention criterion. It is not a laboratory shelf-life, handling, or safety
  guarantee.
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
