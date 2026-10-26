# Condition-locked transfer reference curation — 2026-10-17

The preparation order remains **H + CH4, then OH + H2**. This directory is
scoring-side evidence only; no new chemistry was launched. Both cases are
`curation_incomplete`, with empty `measured_rows` and `ready_for_scoring=false`.
This does not mean experiments are absent. It means their individual rows
have not been verified from the full primary documents.

## H + CH4: original article plus correction are required

[Rabinowitz et al., 1991](https://pubs.acs.org/doi/10.1021/j100155a033)
is the intended 897–1729 K flash-photolysis/shock-tube source. Publisher
metadata was accessible, but direct article access returned HTTP 403.

Crucially, an [ACS correction in 1992](https://pubs.acs.org/doi/10.1021/j100181a600)
is listed in the [publisher issue contents](https://pubs.acs.org/toc/jpchax/96/2),
volume 96(2), page 1010. Its body was unavailable. Its effect on numbers is
unknown; the original fit must not be frozen without checking it.

The [NIST record](https://kinetics.nist.gov/kinetics/Detail?id=1991RAB%2FSUT674-681%3A2)
provides retrieval leads for argon and 0.02–0.04 bar, not primary-audited
conditions. Its displayed rate table is calculated from an Arrhenius fit,
not individual experiments. Its upper temperature and parameter-error
transcription also require reconciliation with the primary documents.
No H/methane numerical fit or generated table was admitted here.

## OH + H2: author fits are not measured rows

The [Orkin et al. publisher abstract](https://pubs.acs.org/doi/10.1021/jp057035b)
supports 200–479 K resonance-fluorescence measurements, two temperature-range
fits and the published overall multiplicative uncertainty factor.
Those coefficients are transcribed separately under `author_fits` in
`curation.json`, **not** under `measured_rows`. No fit evaluations were
manufactured as observations.

The [NASA institutional page](https://esdpubs.nasa.gov/content/Rate_Constant_for_the_Reaction_of_OH_with_H2_between_200_and_480_K)
reproduces the abstract; its full-paper link returns to the publisher.
The [NIST entry](https://kinetics.nist.gov/kinetics/Detail?id=2006ORK%2FKOZ6978-6985%3A1)
offers pressure and confidence-level leads, but neither replaces primary
methods/uncertainty audit. Its generated table and rounded fit exponent
must not be substituted for exact observations or full-precision author fits.

## Unlocking the next benchmark

Provide lawful full copies of the three primary documents, or an accessible
author/institutional archive. Then transcribe and visually verify exact
measurement rows nearest the proposed temperatures: H/methane around
1000/1200/1500 K and OH/H2 around 298/375/450 K. These are selection targets,
not claims that measurements occurred at exactly those values.

Record row-specific conditions, units, errors/confidence semantics, shared
systematic uncertainty, secondary-chemistry treatment and the correction's
effects. For OH/H2, check experimental preparation against declared OH
electronic-state and H2 rotational/nuclear-spin populations. Freeze the
reference bytes and scoring policy separately from reference-free prediction
requests. These are elementary gas-reaction transfer tests, not validation
of a stability observable or of every environment in the roadmap.

No full source bytes were retained; this is an access/curation manifest,
not a reproducible primary-table archive. No production defaults, stored
scientific evidence or existing benchmark files were changed.
