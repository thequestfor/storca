# Definition of a Supported Stability Lane

Registration is not support. An installed program, accepted input schema,
successful calculation, or plausible reaction is not by itself support for a
stability observable.

## Runtime support levels

- `declared_only`: recognized for planning and fail-closed routing; no
  scientific execution or observable claim is allowed.
- `screening_supported`: bounded evidence may be generated, but no verified
  observable may be retained.
- `observable_supported`: an admitted execution may retain only the explicitly
  declared observables after all lane-specific evidence gates pass.

`observable_supported` is conditional support within a restricted method
domain, not proof of experimental accuracy, exhaustive reaction coverage, or
universal applicability.

## Promotion checklist

A lane may be promoted to `screening_supported` only when it has:

- a versioned environment and chemical-system identity;
- explicit discovery or evidence inputs;
- bounded compute and restart behavior;
- standardized failure and truncation states;
- immutable provenance and cache rules;
- regression tests; and
- documented scientific limitations.

A lane may be promoted to `observable_supported` only when it additionally has:

- a native candidate/network representation;
- mechanism-coverage evidence or an explicit restricted external coverage
  claim;
- a path-verification policy appropriate to the reaction class;
- an applicable rate model with dimensions and conditions;
- a network propagator with physical and numerical closure checks;
- a declared observable and units;
- first-class uncertainty inputs and propagation;
- an observable convergence policy;
- method-applicability diagnostics and fail-closed escalation;
- a scientific claim/provenance chain;
- exact restart and scientifically valid cache reuse;
- unit, integration, negative, and restart tests;
- a named benchmark profile with retained validation scope;
- a runnable bounded example; and
- synchronized `METHODS.md`, `LIMITATIONS.md`, requirements, traceability, and
  manifest records.

## Required negative behavior

Every supported lane must demonstrate that it withholds an observable for:

- missing or mismatched environment evidence;
- incomplete or invalid molecular identity and balance;
- unresolved influential candidates;
- unsupported reaction physics;
- failed path or rate verification;
- uncertainty that changes the requested conclusion;
- propagation or closure failure;
- exhausted verification budget before convergence; and
- changed inputs presented to an existing append-only session.

## Benchmark interpretation

Three evidence levels must remain distinct:

1. **Software regression**: schemas, algorithms, failure states, and
   reproducibility behave as designed.
2. **Scientific component validation**: discovery, path, rate, or propagation
   components meet a preregistered benchmark in a restricted domain.
3. **End-to-end observable validation**: the complete lane meets declared
   prediction and calibration targets on condition-matched holdout systems.

A workflow regression profile can support safe execution of a conditional lane
without establishing end-to-end experimental validation. Reports and roadmap
status must state which level has actually been reached.

## Support change procedure

Any pull request that changes a lane support level must:

1. name the affected requirement identifiers;
2. update the runtime registry and machine-readable manifest together;
3. add or update the named benchmark profile;
4. show all positive and negative acceptance tests;
5. update scientific methods and limitations;
6. state the new observable claim and its applicability boundary; and
7. receive an explicit scientific review distinct from ordinary code review.

