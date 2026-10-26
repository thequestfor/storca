# Expanded-network assessment and next-evidence plan

Date: 2026-10-08. Status: **implemented; scoped software/native replay validated,
physical prediction unresolved**.

Implementation evidence and the exact retained-model reproduction are in
[the expanded-network checkpoint](EXPANDED_NETWORK_ASSESSMENT_2026-10-08.md).
The scientific follow-on section below remains outstanding; no P1 lifetime,
calibration, coverage or supported-environment gate is closed.

This is the next bounded implementation slice of P1 (working prediction loop),
P4 (observable-directed allocation), and P7 (method applicability) in
[the comprehensive plan](COMPREHENSIVE_IMPLEMENTATION_PLAN.md). It follows
[the rate-to-network handoff checkpoint](RATE_NETWORK_HANDOFF_2026-10-08.md).
It does not supersede the wider environment roadmap, authorize new chemistry
jobs, or close any scientific support gate.

## Outcome and scope

STORCA should be able to reopen the **whole saved expanded model**, preserve its
different electronic states, assess its declared approximation, and explain
which unresolved chemistry could change parent retention. It should not stop
merely because a narrowly scoped, older certificate does not cover new species.
Nor should it silently extend that certificate.

The acceptance replay is the retained cyclobutene model in:

`runs/p1-cyclobutene-handoff-002/engine/verification-engine/iteration-001-rmg-02be421ca3cb/rerun-result.json`

Its exact contract is 423.15 K, 0.0666611842 bar, pure-parent dark homogeneous
gas, isothermal constant-pressure batch, 3600 seconds, and 95% retention.
It contains seven native species and three reactions / four directed channels:

```text
A = stability(1): singlet cyclobutene
B = C4H6(2): singlet butadiene
I = C4H6(3): triplet biradical

A -> B       conditional forward rate from the verified singlet path
A <-> I      estimated electronic-state-changing entry/return
I -> B       estimated electronic-state-changing escape
```

The remaining four species are retained inert species with zero initial
population. They are not silently deleted to simplify the certificate.

The present `_inventory()` in `storca/bounded_handoff.py` strips numeric native
suffixes. It therefore refuses `C4H6(2)` and `C4H6(3)` as colliding labels. That
refusal is protective; native structural correspondence has already succeeded
for the inserted singlet route. The new work needs snapshot-local full native
identifiers, not another geometry matcher or a universal species-ID rewrite.

Initial applicability is deliberately restricted to closed, isothermal,
conservative one-species-to-one-species first-order networks. Bimolecular,
open, illuminated, nonisothermal, or unsupported environment models receive an
explicit refusal. This is a diagnostic capability, not a declaration that all
stability lanes are gas-only. Existing solution and other-lane work is unchanged.

## Implementation sequence

### 1. Preserve and validate the actual retained network

Add a small `storca/retained_linear_network.py`, shared by reassessment and
diagnostics. Do not replace `ReactionCandidate`, the existing propagation
interfaces, or RMG's native mechanism.

The source-bound snapshot contains:

- Full `chemkin_identifier` for each species, native graph/state evidence,
  charge, multiplicity, and elemental/isotopic composition where available.
- Native reaction index plus direction as the channel ID; full endpoint IDs,
  stoichiometry, reversibility, operative rate, units, and exact conditions.
- Initial fractions and the parent ID resolved through the existing strict
  `reaction_correspondence.resolve_replacement_correspondence()` proof.
- Hashes of the retained model, Chemkin, dictionary, applicable transport,
  replacement library, native identity sidecar, and rate snapshot.
- Versioned policy and explicit limitations. Native IDs are local to this
  hashed model, not permanent identifiers across RMG reruns.

Validate the existing replacement with `flux_verification.verify_applied_replacement()`.
Require modern, successful structural/state correspondence. Do not use
normalized equation strings or a legacy-label fallback after a failed proof.
Preserve existing stereo refusals; this increment does not establish stereo
identity that the native artifacts cannot prove.

Add an isolated read-only operation to `storca/rmg_bridge.py`, accessed through
the existing `rmg_bridge_client.run_rmg_bridge()`. Reuse `_cantera_solution()`
to retain Cantera's actual fixed-condition forward/reverse rate constants and
their source attribution. Check one-to-one native-to-Cantera reaction mapping;
refuse conversions that expand into unsupported multiple reactions. Validate
atom/charge balance, first-order units, finite nonnegative rates, initial
composition, and the actual isothermal reactor policy.

Important: `inspect()` currently exposes a **fitted reverse Arrhenius rate**.
That is not necessarily the reverse rate used during Cantera propagation.
The new model snapshot must use operative Cantera constants, not that fit.
Direct native `kf/Kc` is a cross-check, not an unexamined substitute for the
actual propagated thermochemistry. Retain any discrepancy explicitly.

Tests: same-formula singlet/triplet species stay distinct; reordered numbering
requires new correspondence; reaction orientation reversals yield equivalent
directed models; stale hashes, wrong conditions, missing states, ambiguous
mapping, unbalanced/nonlinear/open models and nonfinite rates refuse admission.

### 2. Issue a new expanded-model assessment, not a broader old certificate

Extend `storca/bounded_handoff.py` with separate entry points, provisionally:

```text
assess_expanded_forward_only_model(...)
validate_expanded_forward_only_model_assessment(...)
```

Preserve the current v1 certificate, its validator and its inventory-change
refusal. An expanded assessment is a new append-only sidecar for the current
model, not a rewrite of `generated-kinetics.json` or of historical evidence.

Reopen the original rate/library ancestry, the conditional forward fit and the
reached-conformer reverse evidence. Recheck actual irreversible insertion and
native parent/product correspondence. Construct a **comparison model** by
adding only the omitted `B -> A` return to the complete current inventory.
The existing `I -> A` return remains explicit and must never be counted as an
omitted channel.

Reuse `reverse_relevance.assess_reverse_return_relevance()` for the same
finite-horizon, model-only omission argument. Bind the assessment to this
expanded model and the named reached-conformer/native-thermochemistry return
alternatives. If the return bound exceeds policy tolerance, or source/model
validation fails, retain the assessment as unresolved. No new irreversible
library, equilibrium assumption, or physical rate distribution is invented.

Keep the evidence graph acyclic:

```text
original path/rate artifacts -> original library -> derived library
                                                -> native identity/model
                                                -> expanded assessment
                                                -> conditional prediction
```

In particular, do not put an assessment hashing the native identity report
back into the library manifest that the identity report already hashes.
The new assessment reopens and recomputes its inputs during validation.

Tests: unchanged v1 refusal; complete four-channel inventory; only the intended
virtual fifth channel in the comparison; unchanged spectator rates; rejection
after any species/channel/rate/source/condition/horizon change; refusal for a
large reverse alternative; no recursion or source mutation.

### 3. Measure branch influence without pretending to know its uncertainty

Add `storca/network_influence.py` for scenario definition, source validation,
result interpretation and evidence retention. Execute its small linear solver
through the configured RMG environment, where SciPy is already used by the
native machinery. Do not add a mandatory application dependency or install a
new package. If the selected environment lacks the solver, return unavailable.

For the admitted network, form the generator `Q`:

```text
for a directed channel i -> j with rate k:
    Q[j,i] += k
    Q[i,i] -= k
dx/dt = Q x
```

Use a stiff BDF solver with constant analytic Jacobian, and Radau plus a
tightened native Cantera baseline as independent numerical checks. Matrix
exponentials are suitable analytical/nonstiff test cross-checks, not an
unconditionally reliable solution for this very stiff model.

Retain explicit tolerance refinement. The very small intermediate population
can be around 1e-15 under the central rates: a scalar `atol=1e-15` does not
establish its escape-flux accuracy. Check parent retention, integrated branch
flux and tiny-population resolution separately. Include the exact requested
horizon rather than treating an interpolated presentation point as a solver
endpoint. Do not clip populations, renormalize, or silently eliminate a fast
intermediate. Failed conservation, positivity, convergence or agreement means
`numerical_unresolved`.

Begin with four saved, deterministic scenarios:

1. Full expanded baseline.
2. Additional unresolved triplet branch removed.
3. Triplet return removed while entry/escape are retained.
4. Triplet escape removed while entry/return are retained.

Record the exact changed channels. Branch removal is a model comparison;
one-direction removal is a counterfactual and changes the original equilibrium
ratio. Neither is a calibrated physical perturbation. Optional paired local
barrier perturbations must scale both reversible directions and explicitly
preserve their ratio. Do not label a directional perturbation as a
thermodynamically consistent uncertainty sample.

Output horizon parent fraction, maximum retained-grid change from baseline,
integrated loss/return/escape flux, intermediate population, solver checks and
scenario definitions. A retained internal policy can start with a 1e-4
absolute-retention reporting resolution and require numerical agreement at
least 100 times finer. This is a prospective diagnostic resolution, not a
physical accuracy tolerance or confidence interval.

Do not reuse the terminal `combined_*_sensitivity` envelope type. New records
must carry:

```text
kind: declared_network_scenario_response
physical_rate_bounds_established: false
uncertainty_calibrated: false
terminal_robustness_eligible: false
lifetime_claim_allowed: false
```

An influential tested scenario identifies a concern. An insensitive scenario
does **not** prove physical irrelevance when applicable rate uncertainty is
missing. All-low/all-high joint multipliers are not automatically extrema in
a branching/recycling network.

Tests: analytical two-state decay and reversible models; competing/recycling
three-state models; the saved stiff rate regime; conservation/positivity;
native baseline agreement; paired-ratio preservation; counterfactual labels;
different return/escape responses; source mutation and numerical-failure cases.

### 4. Audit electronic states before requesting another ordinary path

Add a pre-path endpoint-state audit in `storca/method_applicability.py`, sharing
the existing coupled-multiplicity rules. It must not fabricate a path or a
selected surface to satisfy `evaluate_method_applicability()`.

Return explicit evidence-backed states:

- Missing/invalid endpoint states: `electronic_state_identity_unresolved`.
- No common coupled multiplicity: `spin_changing_rate_model_unresolved`.
- Multiple common surfaces without selection: `surface_selection_required`.
- A compatible selected surface: eligible for further ordinary-path review,
  not proof of single-reference validity or an accepted rate.

Both additional reactions in this saved model must be flagged as connecting
singlet/triplet endpoints. Retain their estimated rates only as declared-model
assumptions. Do not automatically schedule ordinary NEB/TST to resolve a known
spin-changing channel, or treat the absence of diagnostics as applicability.

Produce a source-bound next-evidence record, not an automatically executed job:

```text
action: resolve_spin_changing_branch_model_or_supported_relevance_bound
routes: exact native entry/return and escape reaction references
reason: explicit state audit plus observed scenario response
evidence_required:
  supported channel/electronic-state interpretation
  justified return/escape or thermochemistry constraints
  nonadiabatic rate evidence if this branch remains influential
executable_with_current_backend: false
estimated_cost: null
```

Always retain the independent dominant-rate requirements: requested-pressure
adequacy, statistical counting/conformer basis and applicable rate uncertainty.
No information-gain score is manufactured when uncertainty/cost is unavailable.

Tests: singlet-to-triplet isomerization; ordinary singlet path; radical-fragment
spin coupling (not naive equality of fragment multiplicities); malformed state
data; no specialist backend invoked; no ordinary-TST admission for the flagged
triplet routes.

### 5. Wire reporting and prove the saved-run path

In `storca/verification_engine.py`, retain the original v1 refusal. Attempt the
new assessment only after successful insertion/identity/condition validation
and a specific expanded-inventory scope refusal. Never use this path to
bypass a stale source, failed rate match, wrong environment, or corrupt proof.

When the expanded assessment passes, reuse `conditional_model_prediction` and
attach the **new** model assessment, numerical evidence, branch diagnostics
and next-evidence records. Take the omission bound from the accepted current
assessment, never the historical one-path certificate. Known state concerns
remain explicit assumptions/blockers for physical interpretation.

If the model or diagnostics refuse, retain the actual current model and its
specific failure evidence in an incomplete result. Do not discard it in favor
of only the pre-repair model. Do not trigger repeated RMG regeneration, create
a reinjection loop or redo completed quantum stages. A same-input resume
reuses valid sidecars; changed model/policy/conditions require a fresh assessment.

Reuse the existing main/resume/condition-ladder conditional reporting. Extend
those surfaces only where they currently omit diagnostics or next-action
references. No new user-facing flags are required.

Add `tests/replay_expanded_forward_model.py` to reopen the actual 002 artifacts
into a new output directory using `--rmg-env rmg4_env`. It may inspect/repropagate
the retained native model, but cannot call RMG discovery, ORCA or Arkane. Check
all original bytes and mtimes before/after. Persist a replay result plus a brief
human-readable explanation. This development replay is not a new stability CLI.

## Closure criteria and resource boundaries

This increment closes only when all of the following hold:

1. The saved seven-species model is represented without state/identity collapse;
   all three native reactions / four actual directions remain visible.
2. The original certificate still rejects inventory expansion; the new
   assessment passes or refuses by recomputation on its own exact scope.
3. The saved full-model baseline agrees with native exact-horizon propagation
   within the prospective numerical policy, or explicitly returns numerical
   unresolved. A solver exit alone is not a pass.
4. Return and escape scenarios are distinguishable; both triplet reaction
   records identify electronic-state concerns before ordinary-path selection.
5. The report distinguishes verified path evidence, conditional forward rate,
   estimated branches, model-only retention and missing physical evidence.
6. No verified-flux credit, physical t95/lifetime, completeness, calibrated
   confidence, ladder convergence or supported-environment promotion occurs.
7. Focused portable tests, selected native bridge/identity tests, full existing
   application tests and whitespace checks pass. Actual-artifact replay and
   source immutability are recorded separately from fixture tests.

Use internal versioned bounds rather than more CLI options: initially at most
32 species, 128 directed channels, 16 declared scenarios and two native
reference propagations, under the existing finite bridge timeout. A limit hit
is explicit and does not imply scientific irrelevance. No new quantum or
mechanism-generation budget is allocated here. Numerical error, reverse
omission error and unknown physical/model uncertainty stay separate.

## Immediate scientific follow-on, after this increment

This implementation makes the next decision visible; it does not itself make
the lifetime defensible. Follow-on work is driven by retained evidence:

1. **Pressure applicability audit before a new pressure solver.** The existing
   `benchmarks/stability/p1_cyclobutene_development.json` contains an unaudited
   primary-source transcription and a reported approximately pressure-flat
   regime. Visually audit/materialize the actual source pressure series,
   temperature, composition and precision. Distinguish support for local
   pressure flatness from a quantitatively justified HPL approximation error.
   Missing quantitative evidence remains unresolved. Do not fit this case to
   its measured rate or silently turn the reference into uncertainty calibration.
2. **Resolve statistical recipe before a corrected rate.** Reconcile the
   retained optical/path-counting diagnostic and reached-conformer basis in
   `verification_adapters.py` / `arkane_runner.py`. Do not apply a factor-of-two
   correction from the partition-only diagnostic without physical counting
   evidence. Keep immutable old rate recipes and separate new model comparisons.
3. **Resolve consequential state-changing chemistry.** If the diagnostic
   identifies a consequential branch scenario, seek source-backed state/rate
   constraints or a specialized nonadiabatic calculation. Otherwise report
   that its relevance is not certified; do not erase it because nominal net
   flux is small. An absent executable backend is an explicit next-evidence
   dependency, not an ordinary-TST job waiting to be run.
4. **Supply transferable rate and coverage evidence.** Use the existing
   `RateUncertaintyRecord` contracts and scientific panel, not invented bounds
   or agreement with this development case. Only then test terminal observable
   robustness and prospective benchmark performance.

If finite-pressure modeling is genuinely required, reuse Arkane rather than
adding another solver. The current `VerifiedArkaneRouteSpec.use_pressure_dependence`
explicitly excludes one-to-one isomerization: supporting this case would require
a two-well network contract/builder, both wells' transport/energy-transfer
evidence and explicit pure-parent self-collision semantics. Never silently add
nitrogen as a bath gas. That extension is a separately planned, evidence-driven
follow-on, not part of the present increment.

## Non-goals

No new discovery engine, universal network rewrite, GOAT/NEB/IRC redo, MECP
implementation, master-equation backend, probabilistic calibration platform,
new environment lane, software installation, CLI-option expansion or relaxation
of scientific gates. Broader gas/solution/liquid/solid/interface/exposure work
continues under the existing roadmap; this shared-loop repair neither cancels
nor claims it.
