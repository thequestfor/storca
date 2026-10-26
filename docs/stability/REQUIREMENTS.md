# Stability Requirements

This catalog is the durable requirement source for STORCA stability and
chemical-fate development. Requirement identifiers are never reused. A removed
or postponed capability remains recorded as deferred with its rationale.

## Campaign and system contracts

- **STAB-CAM-001 — Multi-lane campaign request.** Represent chemical identity,
  charge, multiplicity, requested environments, conditions, observables, and
  compute budgets in one immutable request.
- **STAB-CAM-002 — Lane/rung separation.** Treat physical environments as
  independent lanes and condition escalation within a lane as rungs.
- **STAB-CAM-003 — Preparation-first execution.** Permit inspection of every
  admitted, blocked, and planned calculation before expensive execution.
- **STAB-CAM-004 — Exact resume contract.** Resume only when the scientific
  request, methods, environment, and budgets match retained identity rules.
- **STAB-CAM-005 — Cross-lane reporting.** Report conditional lane results
  together without manufacturing a universal lifetime.
- **STAB-CAM-006 — Runnable examples.** Ship bounded preparation, synthetic
  regression, and live smoke-test examples with explicit cost warnings.

## Reaction discovery and native networks

- **STAB-DISC-001 — Native candidate identity.** Store reactants, products,
  stoichiometry, mapping, charge, spin, environment, discovery evidence,
  estimates, applicability findings, uncertainty, and provenance independently
  of RMG.
- **STAB-DISC-002 — RMG discovery adapter.** Preserve RMG assertions and bounded
  restart evidence without making RMG the native network representation.
- **STAB-DISC-003 — Independent graph/xTB discovery.** Retain bounded graph-edit
  proposals, endpoint screens, failures, and truncation without calling them
  barriers or rates.
- **STAB-DISC-004 — Multi-source merge.** Deduplicate chemically equivalent
  candidates while retaining all independent assertions and disagreements.
- **STAB-DISC-005 — Autonomous candidate assembly.** Run selected discovery
  backends and construct a candidate network under one explicit budget.
- **STAB-DISC-006 — Coverage uncertainty.** Represent bounded search completion,
  chemical-class gaps, disagreements, and unresolved candidates separately from
  proof of mechanism completeness.
- **STAB-DISC-007 — Physics-driven escalation.** Invoke bounded AFIR-like or
  automated reaction exploration only for influential uncovered regions.
- **STAB-DISC-008 — User pathway import.** Admit mapped user or curated pathways
  through the same validation and provenance contracts.

## Path and electronic-structure evidence

- **STAB-PATH-001 — Endpoint validation.** Validate molecular identity,
  stoichiometry, atom mapping, charge, spin, and environment at both endpoints.
- **STAB-PATH-002 — Conventional saddle verification.** Require an applicable
  NEB/TS search, one meaningful imaginary mode, and bidirectional IRC endpoint
  matching before claiming a conventional elementary path.
- **STAB-PATH-003 — Barrierless separation.** Preserve reproducible barrierless
  classification separately from rate evidence and ordinary TST.
- **STAB-PATH-004 — Conformer coverage.** Bound endpoint, encounter, and TS
  conformer searches and retain unresolved coverage uncertainty.
- **STAB-PATH-005 — Electronic-state diagnostics.** Detect or require declaration
  of spin changes, excited states, and multireference warning signs.
- **STAB-PATH-006 — Environment consistency.** Audit every path calculation
  against the retained environment and method contract.
- **STAB-PATH-007 — Bounded surface-resolution escalation.** Diagnose unresolved
  paths quantitatively and execute deterministic recovery strategies under one
  immutable orientation, calculation, core-hour, wall-time, and recursion
  budget, preserving a precise abstention when the budget is exhausted.
- **STAB-PATH-008 — Validated intermediate segmentation.** Convert a stationary
  intermediate into a native species and mapped elementary candidates only
  after connectivity, distinctness, charge, multiplicity, and environment
  checks; verify every segment independently and propagate it explicitly.
- **STAB-PATH-009 — Bond-change-guided TS recovery.** Build bounded one- or
  two-coordinate scans from retained mapped bond changes, use only justified
  interior maxima as TS seeds, and still require OptTS, a relevant imaginary
  mode, and bidirectional IRC endpoint agreement.

## Rate models

- **STAB-RATE-001 — Condition-bound rates.** Bind every admitted rate to
  temperature, pressure or concentration as applicable, environment, standard
  state, method, rate theory, and evidence hashes.
- **STAB-RATE-002 — Harmonic TST applicability.** Use conventional Eyring or
  Arkane TST only when its saddle, state, and environment assumptions pass.
- **STAB-RATE-003 — Conformer-weighted rates.** Combine bounded conformer
  ensembles without hiding missing conformer uncertainty.
- **STAB-RATE-004 — Arbitrary temperature.** Generate and validate
  temperature-matched thermochemistry rather than reusing 298.15 K evidence.
- **STAB-RATE-005 — Pressure dependence.** Support RRKM/master-equation rates
  with explicit bath-gas and energy-transfer assumptions.
- **STAB-RATE-006 — Barrierless rates.** Route loose association to capture or
  variational treatment instead of forcing a saddle.
- **STAB-RATE-007 — Tunneling.** Diagnose tunneling sensitivity and retain the
  selected correction or unresolved state.
- **STAB-RATE-008 — Solution standard state.** Retain molecularity, dimensions,
  activities or ideality, and standard-state corrections for solution rates.

## Propagation and observables

- **STAB-PROP-001 — Repaired gas propagation.** Reinject verified rates into a
  validated RMG library and propagate the repaired full network independently.
  Preserve rate dimensions and numerical precision through serialization/reload;
  match the propagated reaction's chemical identity, orientation, and rate
  rather than relying on library presence or unchanged species labels.
- **STAB-PROP-002 — Linear solution propagation.** Propagate verified
  first-order solution networks with concentration, extent, element, charge,
  nonnegativity, and flux closure.
- **STAB-PROP-003 — Nonlinear solution propagation.** Support verified first-
  and second-order reversible solution networks with stiff numerical methods.
- **STAB-PROP-004 — Target reformation.** Propagate networks where the target can
  be consumed and regenerated without using the monotonic-loss shortcut.
- **STAB-PROP-005 — Transport coupling.** Separate reaction and transport
  contributions in explicitly declared multiphase or interfacial models.
- **STAB-OBS-001 — Observable contracts.** Define t95, half-life, concentration
  at time, product yield, dominant product, selectivity, branching, induction
  time, and pathway flux as explicit requests.
- **STAB-OBS-002 — Observable influence.** Quantify flux and sensitivity for the
  requested observable, not merely total network activity.

## Allocation, uncertainty, and convergence

- **STAB-UNC-001 — Uncertainty taxonomy.** Separate mechanism coverage,
  electronic structure, conformers, path search, environment, kinetic theory,
  tunneling, pressure dependence, electronic state, propagation, and numerical
  uncertainty.
- **STAB-UNC-002 — Approximate-rate scope.** Keep priors and screening intervals
  usable for allocation but inadmissible as verified kinetics.
- **STAB-UNC-003 — Uncertainty propagation.** Propagate explicit rate and
  mechanism alternatives to observable bounds or distributions without false
  confidence labels.
- **STAB-ALLOC-001 — Observable-conditioned allocation.** Rank eligible work by
  influence, uncertainty, expected information gain, independent discovery
  value, and compute cost using an auditable policy.
- **STAB-ALLOC-002 — Bounded verification.** Enforce calculation, path,
  iteration, wall-time, and estimated core-hour budgets.
- **STAB-CONV-001 — Prediction convergence.** Require robust observable
  classification, requested uncertainty tolerance, verified influential rates,
  adequate numerical closure, and no unresolved influential applicability gap.
- **STAB-CONV-002 — Explicit abstention.** Distinguish convergence, no crossing
  in a bounded model, unresolved coverage, unresolved rate, method escalation,
  path failure, budget exhaustion, and invalid evidence.

## Environment and method applicability

- **STAB-ENV-001 — General environment identity.** Represent phase/model,
  temperature, pressure, composition, concentration, solvent, light, interfaces,
  and other lane-specific parameters without weakening stricter contracts.
- **STAB-ENV-002 — Gas contracts.** Preserve the existing validated homogeneous
  gas surrogates and explicitly reject liquid, surface, and transport claims.
- **STAB-ENV-003 — Implicit solution contracts.** Preserve solvent model,
  concentration, standard state, light, charge/speciation, and volume identity.
- **STAB-ENV-004 — Bimolecular solution applicability.** Diagnose diffusion,
  encounter complexes, ionic strength, activities, and speciation limitations.
- **STAB-ENV-005 — Microsolvation escalation.** Treat explicit solvent clusters
  as bounded cluster/continuum evidence, not bulk solution.
- **STAB-ENV-006 — Photochemical exposure.** Bind wavelength-resolved light,
  absorption, quantum yields, and photoproduct evidence.
- **STAB-ENV-007 — Interface contracts.** Represent phases, geometry, area,
  partitioning, adsorption, and mass transfer.
- **STAB-ENV-008 — Surface contracts.** Represent periodic model, sites,
  coverage, adsorption, and site balance.
- **STAB-ENV-009 — Electrochemical contracts.** Represent electrode, potential,
  charge convention, electrolyte, and constant-potential approximation.
- **STAB-METHOD-001 — Auditable method selection.** Route reactions by explicit
  diagnostics and user declarations; an LLM is never the physical authority.
- **STAB-METHOD-002 — Fail-closed escalation.** Withhold an observable when a
  required specialist method is unavailable or unresolved.

## Evidence, caching, and reproducibility

- **STAB-PROV-001 — Scientific claim graph.** Trace observable claims through
  propagation, rates, paths, calculations, discovery, assumptions, and gaps.
- **STAB-PROV-002 — Evidence over jobs.** Model calculations as evidence for
  claims rather than making ORCA or RMG jobs the domain architecture.
- **STAB-PROV-003 — Immutable identities.** Hash scientific inputs, methods,
  environment, software identity where available, and calculation artifacts.
- **STAB-PROV-004 — Cache admissibility.** Reuse only evidence whose scientific
  identity matches the requesting lane and claim.
- **STAB-PROV-005 — Append-only sessions.** Retain every allocation, repair,
  failure, and convergence decision without rewriting history.

## Scientific validation and governance

- **STAB-BENCH-001 — Locked benchmark partitions.** Maintain condition-complete
  train/development/holdout systems and retain abstentions and failures.
- **STAB-BENCH-002 — Mechanism discovery metrics.** Score recall, precision,
  dominant pathway, discovery disagreement, and high-level false-positive cost.
- **STAB-BENCH-003 — Kinetic observable metrics.** Score rates, branching,
  products, lifetime, pressure dependence, and activation parameters.
- **STAB-BENCH-004 — Calibration and compute metrics.** Score uncertainty
  coverage, abstention, wall time, DFT jobs, core hours, and information gained
  per compute.
- **STAB-BENCH-005 — Baseline comparisons.** Compare RMG only, current STORCA,
  sensitivity-driven STORCA, and multi-discovery STORCA under matched budgets.
- **STAB-GOV-001 — Traceable changes.** Every scientific change names its
  requirements, tests, benchmark, support effect, failure states, and remaining
  limitations.
- **STAB-GOV-002 — Machine-checked claims.** Validate the lane registry,
  milestone manifest, benchmark profiles, and referenced files in repository
  tests.
