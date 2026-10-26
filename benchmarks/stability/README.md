# Scientific panel evaluation

`scientific_panel.json` retains eleven development curation entries. They do not
yet have complete condition, uncertainty, convention and applicability review.
`p1_cyclobutene_development.json` registers an additional fixed development
comparator whose primary table still needs visual audit. No entry is a locked
prospective holdout, and the evaluator does not freeze or promote the cohort.

The [2026-10-02 source audit](reference_sources/2026-10-01/source-audit.json)
freezes six myrcene/acrylic-acid rate observations from the primary article and
supporting information, with SHA-256 dependencies and visually checked table
multipliers. It preserves the apparent reciprocal-concentration fit convention,
the main-table/SI rounding difference, high reactant concentrations, unspecified
batch pressure and unknown uncertainty. The data are not yet admitted scoring
records or dilute-solution validation. Their source is [Hornung et al.](https://doi.org/10.3762/bjoc.13.15).

The same audit reclassifies the existing CH3CHOO conformer entry as a qualitative
challenge: the primary paper explicitly cannot determine accurate absolute
rate coefficients because its probe and photolysis volumes differ. Transient
decay ordering must not become a fitted absolute rate or lifetime. See
[Lin et al.](https://doi.org/10.1038/ncomms8012).

The [2026-10-03 neutral-solution source audit](reference_sources/2026-10-03-rulisek2005/source-audit-v2.json)
retains eight fitted reversible rate coefficients for furan/maleimide at 300
and 338 K in **acetonitrile-d3 (CD3CN)**, from [Rulisek et al.](https://doi.org/10.1021/jo050759z).
Original tables, conditions, kinetic model and measurement caveats were visually
checked in a public reproduction of the primary paper; six source captures are
hashed. A complete original publisher/author PDF is not retained, and that source
limitation remains explicit. Forward units are L mol^-1 s^-1; reverse units are
s^-1. The reported 0.9545 joint-fit confidence limits are not independent normal
errors or total experimental uncertainty.

These coefficients are usable as **rate-only experimental evidence** at matched
temperature, exact solvent, stereoisomer/channel and concentration-law convention.
That use does not require inventing initial concentrations. It does not establish
microscopic pathway completeness or validate the current solution lane. Actual
initial concentrations, systematic uncertainty, covariance and reaction pressure
are not fully supplied; trajectory/lifetime/product-distribution comparison and
autonomous execution remain pending. A separate preregistered rate-only scoring
protocol and solvent/model applicability review are still required. The audit
neither substitutes CH3CN nor manufactures a lifetime from a reverse coefficient.

Optional `reference_audit` links bind a case to the frozen observation record
and its primary-source artifacts. The loader reopens their hashes and refuses
changed/missing sources, mismatched case identity, and scoring promotion of an
explicitly unready observation set. This integrity check does not establish
the scientific validity of an otherwise unsupported rate convention or model.
The development cohort, numerical acceptance protocol and prospective holdout
remain unfrozen: eleven pending entries plus the provisional P1 comparator, zero
audited scoring-ready cases.

`tests/run_scientific_panel.py` evaluates at most twelve stored case outcomes.
It launches no chemistry, changes no scoring reference, and reserves a new
output directory. For example, from the repository root, using the installed
Python environment:

```bash
python tests/run_scientific_panel.py \
  --panel benchmarks/stability/scientific_panel.json \
  --preregistration benchmarks/stability/p1_cyclobutene_development.json \
  --results current_storca=runs/current-storca-panel-predictions \
  --output-dir runs/panel-evaluation-001
```

## Separate descriptive rate-only comparison

The [Rulisek rate-only development preregistration](reference_sources/2026-10-03-rulisek2005/rate-only-development-preregistration.json)
locks all eight audited coefficient observation IDs and their source hash. It
attaches an in-memory overlay to the existing pending case; it does not change
the source audit, the full-case readiness gate or the eleven pending full-case
references. Missing initial concentrations do not block this **coefficient-only**
comparison. They still block reconstructing the experimental concentration-time,
lifetime and time-dependent branching experiment.

An optional internal benchmark-driver argument enables the separate slice:

```bash
python tests/run_scientific_panel.py \
  --panel benchmarks/stability/scientific_panel.json \
  --rate-preregistration benchmarks/stability/reference_sources/2026-10-03-rulisek2005/rate-only-development-preregistration.json \
  --results current_storca=runs/current-storca-panel-predictions \
  --output-dir runs/rate-only-evaluation-001
```

This launches no chemistry and grants no calculation budget. Prepared requests
contain only genuine initial furan/maleimide identities, physical conditions and
a generic coefficient request. Rates, source-channel/adduct answers and joint-fit
limits stay in scoring. Native structural and method admission remain required;
no adduct or reaction route is supplied as an execution input.

The separate `rate_only` report retains all eight expected rates, including
missing, invalid, unresolved, outside-applicability and crashed outcomes. It
extracts central values from reopened production rate/path artifacts rather than
trusting a number in a benchmark packet, and requires calculation hashes,
scoring-only audited structural assignment and artifact-backed method assessment.
Exact temperature, CD3CN, direction/molecularity, units, standard state and
concentration-law convention must match. CH3CN substitution is explicitly
unresolved under this protocol. Conditional method comparisons are separately
counted; case-level compute is counted once, not eight times. See the
[record contract and retained checks](reference_sources/2026-10-03-rulisek2005/rate-only-checkpoint.md).

Errors and source-fit-limit comparisons are **descriptive**, not a scientific
acceptance pass. Accuracy, coverage and cost thresholds remain null; there is
no uncertainty calibration or population inference. The source's 0.9545 joint-fit
limits never become independent errors or prediction confidence. No lifetime,
retention, product-distribution, network-coverage or support claim is released.
This slice is not an autonomous solution demonstration or closure of item 2.

An absent prediction directory is allowed: it produces `missing` eligible
outcomes, never benchmark passes. To include a retained production P1 result,
add `--p1-result current_storca=/absolute/path/to/p1-run/result.json`.
The run must retain its matching `execution-contract.json` and
`preregistration.json`. The adapter imports only a released verified lifetime
and the driver's measured elapsed time. It does not manufacture a rate, product,
core time, job count, or verified observable from screening evidence. Failed and
unresolved runs are retained and charged for measured elapsed time.

Repeat `--results METHOD=DIR` for RMG-only, current STORCA and proposed methods.
`comparisons` declares matched resource ceilings only when each retained
execution budget has finite wall time, positive integer cores, and identical
declared ceilings. The evaluator itself is not a baseline execution engine;
the existing bounded P1/lane drivers perform calculations.

The output separates `prediction-requests.json` from `report.json`:

- Request exports contain allowlisted parent identity, physical conditions,
  requested observables, declared budget and request hashes. Reference-only
  rates, products, mechanisms, sources and threshold answers are excluded,
  including extra fields nested in initial species records. These are prepared
  requests; admission to execution still belongs to the existing lane drivers.
- Reports retain every eligible missing, crashed, invalid, unresolved and
  outside-applicability outcome, together with the eleven pending references.
  Prediction coverage uses all cases whose frozen scoring record expects a
  prediction, including unsuccessful cases. Descriptive error statistics name
  their available-score and missing-score counts.
- Rate and lifetime errors use log10 values with exact units and rate convention
  matching. Derived lifetime/retention references retain their origin. Product,
  branching and curated net-reaction scores appear only where a reference is
  explicitly supplied; observed products are not assumed elementary mechanisms.
- False reassurance from incorrect retention and from a known consequential
  omitted net reaction are separate fields. Missing mechanism evidence is
  unassessed, never evidence of zero omission risk. Justified abstention needs
  a preregistered inapplicability expectation, retained reason and valid evidence;
  independently assessed applicable cases that abstain count as false abstentions.
  Unresolved applicability, including the present P1 review, remains unassessed
  for false abstention rather than penalizing justified withholding. Crashes stay separate.
- Intervals retain their physical-bound, sensitivity-range or empirical-interval
  label. Covering one comparator point is a diagnostic, not calibration.
- Wall time, core time and high-level job counts require retained artifact hashes
  and a stated measurement basis. Missing costs are unknown, not zero. A core
  ceiling multiplied by elapsed wall time is not measured core time.

For additional development cases, an explicit `scoring` record is required.
It contains `reference_status` (`audited` or `provisional`),
`expected_disposition` (`prediction` or `outside_applicability`), provenance,
`applicability_assessment` (`applicable`, `outside_applicability` or `unresolved`),
the experimental uncertainty or null, and the actual rate/lifetime/product/
branching/net-reaction comparators. A rate has `value`, `units` and `convention`;
a lifetime also declares `kind` (`experimental_measurement` or
`derived_comparator`). An independently justified consequential-path list is
distinct from ordinary curated net reactions. Acceptance thresholds must be
registered before predictions; the evaluator does not choose them.

Stored predictions use schema version 1, kind
`storca_scientific_panel_prediction`, `case_id`, `request_sha256`,
`disposition` (`prediction`, `unresolved`, `outside_applicability` or `crashed`),
`reason` for withheld outcomes, `execution_budget`, and `evidence_artifacts`
with `path` and plain SHA-256 digest. `observables` contains rate/lifetime
records, a numeric `parent_fraction`, a product identity string, or a product
identity-to-branching-fraction mapping. Optional `candidate_reactions` use the
existing discovery benchmark's `reactants` and `products` representation.
`compute` uses `wall_time_seconds`, `core_time_seconds` and
`executed_high_level_jobs`, plus per-field `measurement_basis` and artifact
records. Relative artifact paths resolve against the prediction directory.

Prospective holdout scoring is refused by this development slice. Source audit,
condition-complete frozen cohorts, baseline calculations, independent rate-error
profiles, locked blinded holdout procedures and formal sample-size/coverage/cost
thresholds remain required scientific work. A small curated panel supplies no
population calibration claim. Every report keeps support and calibration false.
