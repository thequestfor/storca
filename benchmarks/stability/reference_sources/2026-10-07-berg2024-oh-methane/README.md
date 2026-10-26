# Independent source curation: OH + methane

Frozen reference-only evidence and an unmerged experimental slice. No shared panel
amendment, predicted mechanism, acceptance threshold, calibration claim, or chemistry
execution is authorized by this directory. Its projected request contains only starting
identities and the nominal reported physical conditions, with an empty execution budget.

Primary source: Berg et al. (2024), *Atmospheric Chemistry and Physics* **24**,
13715–13731, https://doi.org/10.5194/acp-24-13715-2024.
The publisher distributes the paper under CC BY 4.0.

The original paper and supplement were independently visually audited on complete
relevant pages (main PDF pages 1-7; supplement pages 1-2). `source-audit.json` retains all
14 actual Table S1 methane observations, their printed errors, canister identities,
conditions, source hashes, and the unresolved uncertainty definitions. Fitted values and
the paper's cited JPL/IUPAC recommendations are not substituted for measurements.

`experimental-rate-slice.json` uses the existing scientific-panel curation schema and
contains one audited descriptive central-rate comparator: the first original Table S1
methane row in the first canister column, selected by source table order before any
prediction inspection. It is 1.40e-14 cm^3 molecule^-1 s^-1 at 340.3 K and nominal
1000 hPa (1.0 bar) in humidified synthetic air. All other points remain frozen in the
source audit; this scalar slice does not validate the full temperature dependence.

The printed S1 plus/minus error has no explicit sigma definition. Methods report
12% methane accuracy at 2 sigma; Table 1's room-temperature 6.1+/-1.0 (times 1e-15)
is explicitly labelled 2-sigma accuracy and differs numerically. All are retained
without selecting a smaller interval. Confidence probability and shared systematic
covariance remain unknown. Exact N2/O2 fractions and per-point pressures are absent;
this is not a complete admitted flow-reactor input.

The shared `scientific_panel.json` is unchanged (SHA-256
`edba373ce60f0c0d7959b968a1a5b687974a85c3a3e2287e66993075195bfbc9`).
It already has 11 entries, or 12 with optional P1; this separate slice does not remove
pending entries or exceed that evaluation bound.

Existing source audit and input projection gates bind the point to the original PDFs.
No product, expected atom mapping, fitted rate law or reference number appears in the
prediction request. The existing gas TS/frequency/IRC-to-Arkane interface is a plausible
producer boundary, not evidence of a valid prediction for this case: total symmetry
degeneracy, tunneling, electronic-state weights, recrossing and finite-pressure
applicability must remain honest conditional assessments. No native calculation was
run during curation.
