# Independent gas reference curation - 2026-10-07

The smallest plausible measured comparator is OH + H2O2 at 296 K, but it is
**not scoring-ready**. Only new evidence files were created. No shared panel,
threshold, partition or reference bytes were changed, and no chemistry ran.

The publisher-indexed abstract resolves the existing malformed fit notation:
`k(T) = 6.8e-13 exp(285/T)` in `cm^3 molecule^-1 s^-1` over 96-296 K, with
source parameter plus/minus values `1.0e-13` and `27 K`. These are empirical-fit
parameters, not a temperature grid of independently measured observations.

Indexed reproduced journal text identifies a direct 296 K lead of
`(1.78 +/- 0.19)e-12 cm^3 molecule^-1 s^-1` and describes the reported error as
twice combined standard uncertainty. This lead remains unaudited: its full
conditions, complete table and error definitions were not visually verified.
No exact confidence probability, parameter independence or method calibration
was inferred. The large erroneous-looking NIST pressure field was rejected,
not silently repaired.

Why not the neutral unimolecular alternatives? Barnard 1974 is physically
relevant to cyclopropane/cyclobutane, but the available primary abstract does
not expose quantitative points or the temperature/falloff corrections. Its
publisher PDF link redirects to the unauthenticated abstract. Lin 2015 already
states that accurate absolute rates are unavailable, so qualitative ordering
does not repair this gap. The existing JPL limiting-rate comparator is not an
independent finite-pressure measurement.

The exact next dependency is DOI 10.1021/jp030424q, complete original pages
10643-10646. Inspect Table 1, measured 296 K conditions, gas density/pressure,
composition, signal normalization and combined-uncertainty definitions. A
legitimately obtained copy or complete journal-page reproduction is sufficient
for that bounded audit; another secondary database is not.

The retained raster asset contains figures and table rules but not the HTML
table text, so it does not satisfy that dependency. The white display render
only makes the original transparency visible; it reconstructs nothing. Browser
discovery found no available browser, and complete public HTML retries failed.

The proposed amendment is metadata-only and has not been applied. Even after a
source audit, a descriptive coefficient comparison would not establish a
measured lifetime, reaction-network completeness, microscopic mechanism,
finite-pressure applicability, calibrated prediction interval or prospective
holdout. Reference values and expected pathways must not enter production.
