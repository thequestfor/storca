# Descriptive rate-only comparison checkpoint - 2026-10-03

This bounded item-2 slice prepares future condition-matched rate comparison.
It is not a completed scientific cohort, a thresholded acceptance test, an
autonomous solution run or an uncertainty/coverage/support claim. No chemistry
was launched, and no real predicted coefficient or confidence was invented.

## Preserved boundaries

The original source audits, six primary-content reproduction captures, pending
manifest records and earlier curation evidence remain unchanged. The separately
frozen `rate-only-development-preregistration.json` references audit v2 by its
SHA-256 and registers all eight main furan/maleimide observations. A source's
`source_fitted_rate_only_comparison_available` permits this specific descriptive
overlay; `quantitative_scoring_ready=false` still blocks full-case promotion.
The overlay is applied only in memory from a separate pre-prediction,
development-only protocol. Direct manifest injection and holdout reinterpretation
are refused. Dropping an inconvenient rate is refused.

These are source-fitted **concentration-law** coefficients, not proved
activity-based intrinsic rates or independently verified microscopic mechanisms.
Numerical initial concentrations are unnecessary for retaining/comparing a
coefficient, but necessary for reconstructing a specific concentration-time or
lifetime experiment. No inferred concentration, lifetime or product branching
is added. CD3CN and 300/338 K are preserved exactly. The protocol does not admit
ordinary CH3CN as an isotope surrogate; that would require a separately sourced
model-discrepancy/applicability review and a pre-prediction amendment.

## Input isolation

`prediction-requests.json` exports a `storca_scientific_panel_rate_request` with
only initial neutral singlet furan/maleimide identities and physical rate
conditions. Adduct structures, source channels, expected products, coefficient
values and fit limits are not exported. Structural admission remains native,
and `execution_admitted=false` remains fixed. An empty budget grants no chemistry
permission; a separately authorized finite production budget is still necessary.

## Future retained-result contract

A result file named `<case_id>.json` uses schema 1, kind
`storca_scientific_panel_rate_prediction`, matching `case_id` and `request_sha256`,
`disposition`, a finite `execution_budget` with at most four cores, optional
case-level measured `compute`, and a bounded list `rates`. Each rate record has:

- `observation_id`, identifying the **scoring-only** registered observation;
- `disposition` and a reason for withheld/crashed outcomes;
- `production_artifact`: a path/SHA-256 pair for an actual retained producer JSON;
- `assignment_artifact`: a path/SHA-256 pair for a scoring-only structural audit;
- `calculation_artifacts`: the retained path/geometry/rate calculation dependencies.

The evaluator reopens producer JSON and extracts `verified_reaction` rather
than trusting a benchmark-supplied `value`. It reads the existing unimolecular
producer's `rate_constant_s`, or a second-order record's `rate_constant` and
`rate_units`; molecularity follows the record's actual reactant contract.
Rate/path claim identities, temperature, one-molar theoretical standard-state
correction, upstream endpoint/one-imaginary-mode/IRC statuses and calculation
hashes must be consistent. A prepared or speculative path is not sufficient.
This checks retained producer contracts; it does not repeat quantum calculations,
parse Hessians again or establish DFT accuracy.

The structural assignment uses schema 1, kind
`storca_rate_comparator_assignment`, status `audited_structural_assignment`,
the matching observation/channel/adduct/direction, production SHA-256 and
rate/path claim IDs, `rate_convention=concentration_mass_action`, exact
`solvent_isotopic_formula`, and hashed `evidence_artifacts`. Its
`method_applicability` has an artifact-backed explicit `applicable` or
`conditional` assessment, reason and limitations; unresolved/outside assessments
cannot supply a point score. Such records require an explicit algorithm or
scientific structural/model audit, not an LLM's label. Automatic stereoisomer
assignment and the actual bimolecular producer are **not implemented by this
evaluator** and remain real dependencies before their rates can be admitted.

The eight-rate denominator includes missing/invalid/withheld results. Partial
predictions remain partial. Extra or duplicate identities cannot silently replace
a registered coefficient. Temperature, solvent/isotope, units, molecularity,
convention, hash or applicability failures withhold the relevant point score.
Macroscopic result fields are not admitted. Conditional point comparisons stay
separately counted. Wall/core time and high-level jobs retain the existing hashed
measurement-basis requirements, counted once per system rather than per rate.

## Meaning of the report

The `rate_only` section reports signed/absolute log10 point errors and whether
the point lies within the **source-reported joint fit limits**. That boolean is
not a prediction-confidence score or calibration. The probability 0.9545,
non-independent interval semantics, absent covariance and unknown systematic
uncertainty remain explicit. No independent normal distribution is generated.
Accuracy/coverage/cost thresholds and `scientific_acceptance_pass` remain null.
Rate agreement proves neither microscopic mechanism nor network completeness.

## Retained checks and actual progress

Before implementation, three selected contract checks failed because the new
rate-preregistration interfaces did not exist. This is a feature regression,
not an observed chemistry-engine failure. After implementation, **32 focused
scientific-panel tests pass**. Positive producer examples exist only in explicitly
synthetic temporary arithmetic/contract fixtures; they are not retained as real
chemical predictions. Negative checks cover omission/cherry-picking, duplicates,
holdouts, injected overlays, changed sources/production, isotope/T/unit mismatch,
unresolved applicability, invalid path booleans and macroscopic promotion. A cost
fixture verifies once-per-system accounting.

`runs/scientific-panel-rate-only-001/no-job-replay/` contains a real invocation of
the offline driver, with no prediction files supplied. It reports **eight
expected/missing rate-only observations, zero scored**; the original full-case
panel remains **12 entries including P1, eleven pending references, zero audited
full-case comparators and one provisional P1 comparator**. All support,
calibration, holdout and frozen-complete-cohort flags remain false. Evaluator
chemistry jobs are zero; costs and uncertainty are unknown, not zero.

Hashes and before/after checks are retained in `rate-only-evidence.json`. This
implements a bounded comparison interface and locks its descriptive reference
scope. Real condition-matched producer artifacts, structural assignments,
applicability evidence and a scientifically justified acceptance protocol are
still required; tests do not close the expanded goal.
