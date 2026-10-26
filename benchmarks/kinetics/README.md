# Reaction-discovery benchmark scaffold

The cross-environment development panel is being curated in
`benchmarks/stability/scientific_panel.json`. It links existing specialist
references and adds primary experimental sources; it is not yet frozen or a
completed prediction benchmark. Missing conditions and uncertainty are explicit.

This fixture set compares three candidate pools:

- A — retained RMG candidate routes;
- B — bounded graph-edit proposals screened by GFN2-xTB endpoint energies;
- C — their source-preserving merged pool.

The metrics are candidate recall, unmatched candidates, raw/structurally valid
candidate counts, xTB endpoint calculations, and the number eligible for
high-level verification. The merged comparison also reports source agreement,
RMG-only and graph-only counts, recall gain, and candidate inflation. They do
**not** validate barriers, rates, dominant
mechanisms, or mechanism completeness. Each fixture therefore labels its
limited scientific purpose explicitly.

Expand this scaffold only with curated systems whose expected net
connectivity, charge, spin assumptions, and experimental/mechanistic context
are recorded alongside the fixture.

`coverage_risk_manifest.json` separately preregisters the schema and primary
safety metric for observable-conditioned coverage decisions. Its primary
metric is false reassurance: an observable claim made while a provenance-
backed influential reference pathway is missing. The empty
`curation_required` manifest is intentional; synthetic fixtures may validate
the scoring software but must not be advertised as scientific calibration.
The schema separates missing-path false reassurance from false abstention on
reference-covered cases, reports both by chemistry stratum, and requires
RMG-only, pre-STAB-008, STAB-008, and multi-discovery baselines with compute
accounting. The empty manifest cannot pass its calibration gate.

## Specialist limiting-rate source audit

`specialist_rate_manifest.json` remains a development reference set, with no
locked holdout or completed finite-pressure experimental curve. Its four
association reactions have two evaluated JPL 2011 limiting coefficients each;
the eight numerical anchors are not eight independent measurements. The
original NASA/JPL evaluation, selected-page visual audit, predecessor manifest
and amendment record are retained under
`reference_sources/2026-10-03-jpl2011/`.

The published uncertainty factor is referenced to 298 K, even though the
limiting-rate table is anchored at 300 K. Applying
`f(T) = f(298) * exp(g * abs(1/T - 1/298))` gives:

| Reaction | f(298) | g / K | f(300) |
|---|---:|---:|---:|
| H + O2 association | 1.3 | 50 | 1.30145495228246 |
| OH + OH association | 1.5 | 100 | 1.5033594610833196 |
| CH3 + O2 association | 1.1 | 50 | 1.1012311134697739 |
| C2H5 + O2 association | 1.2 | 50 | 1.2013430328761168 |

Only the reference factors and explicit uncertainty/provenance semantics are
amended. Central rates, physical conditions, case partitions and acceptance
thresholds are unchanged; a regression verifies the original predecessor bytes
and reconstructs the prior manifest after removing the amendment.

These are evaluator-assigned multiplicative envelopes, not a specified
probability distribution or exact confidence probability. Low- and
high-pressure coefficients share a reaction-level evidence group; pressure
points generated from an evaluated falloff expression are not independent
experiments. Air is not an N2-only bath, and the OH association coefficient is
not total OH loss. The audit supplies no STORCA electronic-structure,
capture-potential, energy-transfer or mechanism-coverage calibration, and
promotes no lane or milestone. See `source-audit.json` and
`EVIDENCE_2026-10-05.md` in the retained source directory.
