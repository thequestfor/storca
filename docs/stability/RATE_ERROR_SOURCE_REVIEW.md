# Rate-error source review, 2026-10-02

Work item 5 / prediction-first B2. This is a source review and proposed narrow
producer slice, not an admitted calibration profile. It authorizes no rate
multiplier and does not close the scientific acceptance gate. Machine-readable
source facts and deliberately pending fields are in
`benchmarks/kinetics/rate_error_source_review.json`.

The inspected local evidence and primary literature do not provide an admissible
full-rate prediction interval for the current plain B3LYP/def2-SVP, 0.970 project
frequency scale, single-conformer RRHO/Eckart gas method. This is a conclusion
about the evidence reviewed here, not proof that no relevant paper exists.

The local scale calibration is a fit to prominent IR band positions in a mixed
phase reference panel. Its selected 0.986 value is neither a ZPE scale calibration
nor a barrier, entropy, tunneling, or rate-error distribution. Arkane defaults its
frequency scale to 1.0, while the ordinary verification adapter passes the named
method profile's 0.970 value. Match the actual retained route settings, not a
method name alone. The generated input disables inferred hindered rotors and bond
corrections. Its `calibration_required` disposition remains correct.

The actual generated input is more specific than “0.970 frequencies.” In both
inspected local Arkane sources (`rmg4_env`, RMG 4.0.0; `rmg_env`, RMG 3.3.0),
the harmonic and TS imaginary frequencies receive that factor, while ZPE receives
`frequencyScaleFactor / 1.014`, giving **0.956607495** for this profile. This is
what the inspected implementation does, not an independent ZPE calibration.
The balanced per-atom electronic gauge preserves barriers and still lets Arkane
add ZPE. Record the owning software version; a profile name does not authorize
transfer between RMG versions.

## What the sources can support

[BH9](https://doi.org/10.1021/acs.jctc.1c00694) supplies electronic barrier
references and structures. Its published B3LYP comparison uses def2-QZVPP;
its def2-SVP comparison uses PBE0-XDM. Neither comparison is the current method.
The [2022 correction](https://doi.org/10.1021/acs.jctc.2c00362) must supersede
the original reference tables and supporting files. Its corrected B3LYP means
are retained only as descriptive context in the manifest. Means are not tail
quantiles or physical bounds. The authors' reference accuracy estimate is also
an estimate, not a certified worst-case bound.

[Updated BHPERI](https://doi.org/10.1002/jcc.23837) and
[GMTKN55](https://doi.org/10.1039/C7CP04913G) are suitable sources for a small
electronic-barrier pilot after molecularity, spin, charge, and identity review.
The original BHPERI includes cyclobutene; exclude the cyclobutene/butadiene pair
in both directions. Older CBS-QB3 reference values must not replace the updated
Weizmann values. GMTKN55 explicitly excludes ZPVE from its reference energies.
No temperature-dependent empirical rate interval follows from that dataset.

[Barone et al.](https://doi.org/10.1021/acs.jctc.1c00406) assess electronic
energies, ZPE, thermal contributions, and rates together, using a different model
chemistry. Their RRHO/Eckart use demonstrates how a complete rate comparison can
be constructed; it does not transfer its accuracy to B3LYP/def2-SVP.
[Tikhonov](https://doi.org/10.1002/cphc.202400547) provides fundamental-frequency
scaling evidence for dispersion-corrected methods. It does not establish a
full-rate bound for plain B3LYP or its ZPE/entropy treatment.

[RDB7](https://doi.org/10.1038/s41597-022-01529-6) is the strongest rate-pilot
candidate identified: 105 rigid reactions have forward high-pressure RRHO/Eckart
rates at 50 temperatures, 300–2000 K. References use
CCSD(T)-F12a/cc-pVDZ-F12//ωB97X-D3/def2-TZVP with a 0.984 ZPE scale, not the
current method. Rigidity filters require no endpoint rotatable bonds, only
aromatic/three-membered rings, no positive TS frequency below 100 cm−1, and
visual TS review. Some source singlet products are not ground-state singlets.
Use the raw rates, not Arrhenius fits: published fit errors reach 233%; fitting
error is not physical uncertainty. The
[v1.0.1 deposit](https://zenodo.org/records/6618262) declares CC BY 4.0; official
file URLs, hashes, sizes, and the verified temperature-column mapping are in the
manifest. The small CSV has 105 rows, but none is yet admitted to our domain.

The [2025 RDB7 revisit](https://doi.org/10.1039/D5CP01181G) identifies unstable
restricted orbitals and reduced reliability of some original restricted
coupled-cluster references. It classifies RHF-stable, κ-OOMP2-stable intermediate,
and strongly correlated cases, refining only selected difficult references.
The paper and supplementary workflow were reviewed; the refined row table has
not been imported or joined to the 105 rates. Its official ZIP link failed on
review. The classifications are separate for forward barriers, reverse barriers,
and reaction energies. An easy forward barrier alone does not establish suitable
product states or the reverse barrier used by Eckart. Do not treat all old rows
as corrected, or transfer its other-functional statistics to B3LYP. Restrict a
first pilot to audited weak-correlation rows.

## Original rate-reference audit, 2026-10-02

The official 7 MB Arkane archive and small temperature table were obtained without
running chemistry. Their published MD5 checksums match; SHA256 hashes are in the
source manifest. The attributed
`benchmarks/kinetics/rdb7_rigid_reference_inventory.json` retains all 105 original
indices, mapped identities, source-member hashes, temperature labels and
unresolved audit states. All 5,250 forward rate entries agree with the retained
Arkane tables at their displayed precision. The CSV has three/four significant
figures; its maximum relative difference from the printed table is about 0.5%.
That is display-rounding consistency, not a physical uncertainty estimate.

The original forward coefficients are `s^-1`, including the published tunneling
correction. The inspected cache records identify RMG 3.0.0. Species/TS inputs
retain symmetry, optical-isomer and multiplicity declarations; cached species
and TS conformers retain harmonic modes. The small archive does not contain the
global Arkane job input or the original electronic logs, so the full frequency-
scale application, path-degeneracy and numerical-method audit remains open.
The [original rate paper](https://doi.org/10.1038/s41597-022-01529-6) establishes
RRHO/Eckart and high-pressure-limit conventions, not applicability at every
laboratory pressure.

Atom-map/element and total charge conservation were checked before canonical
graph comparison. The 105 rows contain only 104 distinct undirected graph pairs:
indices 358 and 4380 are reverse directions of the same triazole tautomerization.
They must not count as independent calibration units. There is no exact
cyclobutene/butadiene pair overlap in either direction. This does not complete
the exclusion of other evaluation identities, conformers or related scaffolds.

Seven product graphs have two RDKit radical electrons: 1240, 2603, 2732, 3416,
3474, 3475 and 4937. This is a structural warning, not a ground-state assignment.
For row 3474, the cached product graph metadata declares multiplicity 3, while
the retained Arkane input and statistical-mechanics conformer use multiplicity
1. The inventory preserves both fields and marks the inconsistency. It does not
change the source or pretend that the rate was computed for the triplet state.
No row is admitted merely because its graph is formally neutral or singlet.

The [publisher landing page](https://pubs.rsc.org/zh-hans/content/articlelanding/2025/cp/d5cp01181g)
still links the unavailable ZIP. Direct publisher reads returned a Cloudflare
challenge, which was not bypassed. The publisher author-version link, the
[author's Berkeley group page](https://mhggroupberkeley.com/members/xiao-liu/),
the MIT bibliography and bounded official Zenodo searches provided no alternative
refined table. The later official RDB7 v1.1.0 metadata adds 2024 individual heats
of formation; it is not the 2025 orbital/reference refinement. These checks are
an access limitation, not proof that no alternative deposit exists.

All 105 forward/reverse/reaction-energy classifications therefore remain null;
there are **zero admitted calibration rows**. Obtaining the publisher-linked
`RDB7_refined.xlsx` (or an author-authenticated equivalent with redistribution
terms) is the specific dependency for the row join. A hash/version/identity join
must precede inspecting residuals or partitioning eligible cases. It must also
identify refined reference changes and require recomputation of affected rates,
not swap an electronic barrier into an old rate table.

The bounded read-only auditor is
`benchmarks/kinetics/audit_rdb7_source_inventory.py`. Give it a directory containing
the three official files under `rdb7-rigid-v1.0.1.csv`,
`rdb7-rigid-v1.0.1.tar.gz` and `rdb7-temperatures-v1.0.1.csv`, followed by a new
output JSON path. It verifies hashes before parsing, does not execute published
Python inputs, writes exclusively, and always leaves calibration admission false.
It is reproducible source curation, not a rate-error producer or a confidence
claim. If the same directory contains `rdb7-authors-code-v1.0.0.tar.gz`, the
auditor additionally checks the pinned author archive and script hashes and
extracts the recipe through AST parsing. Without that optional file, its author-
recipe field is null; it does not infer the missing recipe. The retained v2
inventory includes that optional author archive audit.

The final v2 inventory was regenerated and replayed byte-for-byte on
2026-10-03. Ten bounded auditor tests pass; the two added rotor/model regressions
had first failed against the earlier implementation (missing mode metadata and
silently collapsed multiple external rotors). This verifies source parsing and
fail-closed bookkeeping, not chemical rates or a calibrated uncertainty model.

The authors' original [v1.0.0 code deposit](https://doi.org/10.5281/zenodo.5652085)
was also obtained (4.6 MB, official MD5 verified) and read without execution.
Its `scripts/run_arkane.py` requests `compute_high_p_rate_coefficient`,
`bac_type=None`, and `three_params=False`, with defaults of 50 temperatures from
300 to 2000 K. It obtains the frequency factor through
`assign_frequency_scale_factor(LevelOfTheory(...))`; it does not record a
numeric value. It requests an `update_qchem_parser` RMG branch and ARC without
exact dependency commits. Consequently this code corroborates the HPL and
two-parameter-fit conventions, but does not resolve the historical harmonic,
imaginary-frequency or ZPE scale application. The archive, script and current
producer source hashes are retained in the source-review manifest. The archive
contains an MIT license; its Zenodo record declares CC BY 4.0.

## Match the actual producer before collecting residuals

This inspection follows `make_arkane_rate_builder` →
`VerifiedArkaneRouteSpec` → `create_verified_arkane_input` →
`validate_verified_arkane_library`, not just the method-profile label.

| Convention | Current production | Retained RDB7 reference | Consequence |
|---|---|---|---|
| Electronic calculation | Plain B3LYP/def2-SVP; `NoAutoStart`; `TightSCF` for NEB/OptTS/IRC, but not explicitly in the separate frequency input | Restricted CCSD(T)-F12a energy on source ωB97X-D3 geometry/frequencies | Exact inputs, actual ORCA defaults/version and orbital diagnostics belong to the recipe; a source-geometry single point is not end-to-end calibration |
| State | Source-bound charge/multiplicity checks; one selected surface | Singlet rates, some product-state warnings and one graph/statmech discrepancy | Closed-shell graph labels do not establish ground states or reference reliability |
| Scale | 0.970 harmonic/imaginary; 0.970/1.014 ZPE in the inspected Arkane sources | Paper reports 0.984 ZPE; historical numeric harmonic application unresolved | Do not silently copy the reference scale into production or claim a matched ZPE protocol |
| Entropy counting | `opticalIsomers=1`; external symmetry inferred; no explicit path degeneracy, owning Arkane reaction defaults to one | Declared and cached optical/symmetry values retained | 69 of 105 cached TSs use two optical isomers; 72 rows have at least one nonunit optical value. This is a declared-policy disagreement, not yet a verified physical error |
| Molecularity/pressure | 1↔1 HPL; 1↔2 mandatory pressure dependence; at most two endpoint species | All raw forward rates HPL, including fragmentation | Only 80 one-product rows have the initial matched HPL shape. Exclude the 24 two-product rows from this first comparison; the remaining three-product row exceeds producer cardinality |
| Rate reported | Three-parameter Arkane fit by default, then exported model evaluated in the selected RMG environment | Raw forward RRHO/Eckart rate table, plus a separate two-parameter fit | Compare the actual exported production rate to raw reference data; also record direct-production TST rates to isolate export-fit error |
| Temperature | Eight baseline points 250–1000 K plus the requested temperature | Fifty raw points 300–2000 K | Start with prespecified overlapping raw temperatures; no 300→298.15 K substitution or 2000 K claim from the default producer |

The auditor now retains cached mode-model, external-rotor symmetry, optical,
spin and E0 metadata with member hashes. All scale-application fields remain
null until the source logs/global dependency recipe are joined. Five of the
80 one-product graphs carry the radical warning (2603, 2732, 3474, 3475, 4937),
leaving **75 graph candidates before state, reference and method audits**, not
75 eligible rows. All 105 orbital joins remain null and admission remains zero.
Only 18 of those graph candidates have all source-declared optical values equal
to one; that count is a diagnostic, not permission to select away a known
production-policy disagreement or an adequate calibration sample by itself.

### Smallest paired-prediction experiment

After obtaining the refined spreadsheet and resolving source-state/refinement
changes, freeze eligible 1→1 neutral rigid ground-state-singlet HCNO channels,
reaction/scaffold groups and partitions before computing residuals. Keep P1,
other evaluation identities, reverse duplicates and prespecified related
structures out. The raw reference values and source TS geometries must remain
unavailable to the P1 discovery/path request.

This is a **channel-conditioned rate benchmark**, not a discovery benchmark:
the calibration calculation knows an independently curated mapped candidate
reaction, but not its reference rate/barrier. Generate endpoints and verify the
path using current bounded production machinery, including the retained method
policy, endpoint identities, one meaningful TS mode and bidirectional IRC.
Source-TS-seeded refinement or single points can be separate component pilots;
they cannot be counted as successful autonomous path searches. Preserve crashes,
missing paths, method escalations and budget abstentions in the frozen cohort.

For each predeclared reference temperature within the default producer's overlap
(300–1000 K), reuse the same completed electronic artifacts and invoke the actual
adapter with that requested temperature. Its fit grid must remain the baseline
grid union that temperature; replacing it with all 50 reference points changes
the production fit recipe. Retain input/model/library hashes, requested RMG
environment, inferred symmetry, optical/spin factors, all three stationary-point
E0 values and scaled imaginary frequency. Measure
`ln(k_exported_production(T)/k_raw_reference(T))` as the primary rate residual
and `ln(k_exported_production(T)/k_direct_production_TST(T))` as the fit diagnostic.
Do not use the source two-parameter Arrhenius fit as the reference or call its
fit covariance physical uncertainty. Discrete tested temperatures do not alone
validate a continuous interval or extrapolation.

Audit forward/reverse E0 barriers because Eckart depends on both. Preserve
detailed balance with the same endpoint thermochemistry; do not separately fit
an inconsistent reverse rate or count reverse direction/temperatures as new
independent calibration units. Audit entropy/path counting before accuracy
admission. If that requires a physical producer fix, freeze a new recipe and
regenerate its predictions; do not absorb a known counting defect into generic
“DFT error” or tune the reference to match it.

Even a completed experiment yields theory-relative full-rate disagreement, not
an experimental prediction interval. The smallest immediate blocked prerequisite
is still the hash/version/identity join to author-authenticated
`RDB7_refined.xlsx`. Reference-scale/state audit, independent predictions and
prospective validation then remain required; importing the sheet is not closure.

### Why this pilot cannot automatically close P1

The RDB7 rigidity filter admits only aromatic or three-membered rings. P1
cyclobutene contains a nonaromatic four-membered ring, so it is outside that
structural sampling domain even though its exact graph pair does not leak into
the table. Its frozen 423.15 K condition is also not a raw RDB7 temperature point.
Independent four-membered-ring/conformer and condition-matched full-rate evidence
is needed before transferring any RDB7-derived profile to P1. The separate
corrected pericyclic electronic-barrier pilot may help assess transfer, but
cannot supply entropy/tunneling/physical-model error by itself. Do not promote
RDB7 curation into a calibrated cyclobutene lifetime.

## Admissible producer slice

The first candidate slice is an electronic-barrier error producer for neutral,
closed-shell, ordinary unimolecular pericyclic gas reactions. The final chemical
scope must follow the audited eligible rows; do not call all pericyclic reactions
unimolecular or pool ion/radical transfer, association, and rearrangement errors.
Use corrected source references, retain source hashes and electronic-energy
conventions, and pair them with exact B3LYP/def2-SVP predictions. A single-point
benchmark on source geometries calibrates only that single-point protocol;
production geometry/path errors require separate end-to-end calculations.

Importing references cannot supply missing method predictions. No chemistry was
run during this review. BH9/BHPERI source-file import and redistribution license
review, exact eligible reaction identities, and overlap audit remain pending.
RDB7's CC BY 4.0 declaration does not settle the separate refined-table license.

A second, more direct candidate is a **theory-relative full-rate log-error
producer** for an audited rigid RDB7 subset. It compares exact production-method
rates to the retained higher-level rate model, including frequencies and Eckart
rather than a barrier-only exponential approximation. It can measure the
production method's disagreement with that reference model. It does not validate
errors common to both models, such as RRHO, Eckart, recrossing, or the
high-pressure assumption; an empirical physical prediction interval still needs
independent reference-model accuracy evidence.

Before that producer can run, it needs original raw rate/temperature artifacts,
retained reference inputs and software versions, the row-level orbital-stability
and reference-refinement join, ground-state and mechanism audit, exact
production endpoint/TS/frequency evidence, and condition-matched generated rates.
Bind geometry optimization, numerical grid/convergence, orbital stability,
frequency-scale application (ZPE versus thermal modes), symmetry, path
degeneracy, conformer policy, units and Eckart settings. Changing these choices
creates a different method profile. None of the reviewed sources supplies the
current method's missing rate predictions. The manifest deliberately has empty
eligible IDs, predictions and residuals, not invented pilot cases.

Freeze the eligible domain and partitions before examining residuals. Group
forward/reverse directions, temperature points, conformers, and duplicate or
related source entries by independent reaction/scaffold family. Cyclobutene and
all panel evaluation identities belong outside calibration. Development cases
already inspected are not prospective holdouts. Chemical-family transfer and
exchangeability remain assumptions to evaluate, not consequences of a dataset's
size.

An electronic energy perturbation with other factors fixed can produce the
deterministic sensitivity `k/k0 = exp(-delta_E/(R*T))`. Label that result as an
electronic component sensitivity. Eckart also changes with barrier height/shape;
rerun the full rate expression for its sensitivity. Neither calculation covers
ZPE, entropy, low modes/rotors, symmetry, degeneracy, conformers, recrossing,
tunneling model error, collisional falloff, or missing mechanisms. It must not be
serialized as a full `method_calibration` rate bound that closes observable
convergence. A full-rate profile requires independently justified evidence for
those applicable components, preferably exact-method experimental log-rate
residuals at matching conditions, and prospective validation.

## Sample-size limits and prospective checks

For an absolute residual split-conformal interval, rank
`ceil((n_cal + 1)*(1-alpha))` must not exceed `n_cal`; otherwise use an unbounded
interval. Thus 95% marginal coverage needs at least 19 independent calibration
units even for a finite maximum-residual interval. This is a mathematical minimum,
not evidence of usefulness or conditional coverage. Follow
[Lei et al., Section 2.2](https://arxiv.org/abs/1604.04173) and its exchangeability
assumptions. Prespecified chemical strata need their own counts; correlated
forward/reverse barriers or many temperatures do not multiply the number of
independent reactions.

Separately, with zero observed misses in `n_test` independent prospective trials,
the one-sided 95% binomial upper bound on miss probability is
`1 - 0.05**(1/n_test)`. At least 59 independent trials are needed to put that bound
below 5%. This arithmetic is an illustrative design requirement, not an adopted
acceptance threshold or a measured result. With nonzero misses, use the exact
binomial bound and report the count. A narrow subset of 26 BHPERI barriers before
exclusions cannot alone establish that prospective guarantee. Even zero errors
does not bound chemistry outside the declared sampling domain.

Until method-matched residuals and condition-specific full-rate evidence exist,
the automatic output should name the missing components and preserve
`calibration_required`, with null rate multipliers. A refusal registry or a
component-only pilot is engineering progress, not completion of work item 5.
