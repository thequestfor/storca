# Neutral solution reference curation checkpoint - 2026-10-03

This is source curation for first-ten item 2, not autonomous solution production,
a scored development cohort, a prospective holdout or a support/calibration claim.
No native, ORCA, RMG, Arkane, rate-generation or propagation job was launched.

## What was recovered

[Rulisek et al., J. Org. Chem. 2005, 70, 6295-6302](https://doi.org/10.1021/jo050759z)
reports reversible furan/maleimide kinetic coefficients in **CD3CN**, under
argon in screw-cap NMR tubes, at 300 and 338 K. Relevant experimental/model/table
passages of the publicly indexed author proof were read. Original tables,
conditions, model exclusion and fitting limitations were then visually checked
in the public journal-content reproduction. Six captures are retained, reopened
and hashed in `source-audit-v2.json`.

Eight positive point estimates and their reported joint-fit interval half-widths
are frozen. Four are bimolecular forward coefficients in L mol^-1 s^-1; four
are first-order reverse coefficients in s^-1. Table footnotes specify joint
confidence limits at probability 0.9545. These are not independent marginal
Gaussian intervals, total experimental uncertainties or a rate-error profile.
Table 1 maleic-anhydride observations are retained as ancillary source evidence,
not another panel case; its approximately-zero reverse exo entry is not a numeric
hard zero or an upper bound. A narrative/Table 1 reverse-rate unit discrepancy
is explicitly retained, with equations and the table determining the convention.

The original v1 audit remains unchanged. V2 corrects only the wording about
source reading: relevant primary passages were read, not a retained complete
original PDF. The predecessor hash is checked by a regression; all numerical
transcriptions and captures are unchanged.

## Rate-only evidence versus an executable stability test

The source-fitted coefficients are usable as rate-only experimental evidence
under matched temperature, isotopic solvent, channel/stereoisomer and the
authors' concentration-law convention. Quoting/comparing a coefficient does
not require inventing an initial concentration.

This does **not** make the present case scoring- or execution-ready:

- Exact starting concentrations and reactant volumes are not tabulated;
  plots and an ideal-volume preparation assumption are not exact input values.
  Trajectory/lifetime/time-dependent branching scoring needs actual data or a
  separately frozen uncertainty-aware digitization protocol.
- CD3CN is preserved. Ordinary CH3CN parameters require an explicit applicability
  decision and discrepancy assessment; no automatic solvent substitution occurs.
- The fitted reversible coefficients do not prove an elementary microscopic
  mechanism or a complete network. Direct endo-exo conversion was excluded from
  the authors' model by lag/theoretical arguments, not universally disproved.
- Pressure, joint covariance, activities and complete systematic uncertainty
  remain unknown. The paper notes NMR calibration/mass-balance problems in longer
  300 K measurements; those are not covered by assigning independent fit errors.
- A separate pre-prediction scoring protocol and domain-specific acceptance
  criteria remain required. No STORCA prediction was used to choose the numbers
  or thresholds; no threshold has been created in this audit.
- Direct PDF download returned 403; institutional originals were explicitly
  audience-restricted and were not requested. Public HTML retention returned
  503. Retained visual reproductions are not misrepresented as full originals.

For a first rate-only development comparison, the endo reverse coefficient at
338 K is a useful *candidate* (2.05e-4 s^-1, reported joint-fit half-width
2.2e-5 s^-1), not an admitted single-parent lifetime test. Recombination and the
exo channel preclude silently deriving the experimental fate from that one rate.

## Actual checks and counts

Two new data/contract tests fail before the transcription because the case does
not yet exist; this is a missing reference feature, **not an observed chemistry
engine defect**. After curation, all 24 scientific-panel tests pass. Existing
integrity regressions reopen primary artifacts and reject changed sources or
promotion of unready source observations. No broad suite was run.

Retained evidence is in `regression-evidence.json` and
`runs/scientific-panel-rulisek-source-audit-001/`. The latest offline evaluator
report has 12 entries including P1: **11 pending references, zero audited scoring
cases, one provisional P1 comparator**, and one missing eligible prediction.
Support, calibration, frozen-complete-cohort and holdout evaluation remain false.
There are zero chemistry jobs. Missing costs and uncertainty are not zero.

The previous Hornung/Lin audit and all earlier source artifacts are unchanged.
The mixed-hydrocarbon CPD and ionic saponification entries remain pending; neither
was relabeled as a pure-solvent neutral-lane success. Failed Claisen/P1 source
retrieval remains a source-access limitation, not evidence of pathway absence.

The next necessary decision is a **separate rate-only preregistration** with
explicit isotopic-solvent/model applicability, or an actually condition-complete
neutral rearrangement source. This audit alone does not close item 2, item 9 or
item 10, and does not justify an autonomous chemistry launch.
