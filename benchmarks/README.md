# STORCA benchmarks

## Stability and kinetics

`stability/milestones.json` is the machine-readable stability roadmap and
support-traceability manifest. It maps every registered lane to its roadmap
milestones and named regression profile, and maps every stability requirement
to a milestone. Repository tests require it to agree with the runtime lane
registry. The manifest is governance metadata, not experimental validation.
`kinetics/specialist_rate_manifest.json` separates limiting-rate anchors,
finite-pressure falloff curves, method-applicability challenges, locked
holdouts, and reference-data-pending cases. Falloff predictions must match the
exact pressure grid and are scored by point error, log-rate curve RMSE, and
pressure span. A locked holdout requires a content hash, protocol hash,
custodian, and per-provider evidence that predictions were frozen before
reference access. Pending cases categorically prevent support promotion; a
source citation without independently reviewed numerical values is not a
benchmark result.
The executable stability-campaign examples include no-ORCA synthetic linear and
nonlinear mass-action solution smoke tests plus a mixed gas/solution contract.
The nonlinear fixture verifies an analytical equal-concentration second-order
solution, dimensional rate rejection, diffusion ceilings, reversibility,
stiffness routing, uncertainty scenarios, and closure; it is not experimental
validation. The mixed
example launches live gas chemistry when executed and is not a cheap regression
command.

`kinetics/reaction_discovery_manifest.json` is an intentionally small
candidate-recall scaffold for RMG, graph/xTB, and merged reaction discovery.
It does not establish mechanism completeness, kinetic accuracy, or stability
prediction accuracy; see [`kinetics/README.md`](kinetics/README.md).

`kinetics/specialist_rate_manifest.json` is the development acceptance harness
for capture/VTST/master-equation providers. It retains unit-specific limiting
rate anchors, reference provenance, uncertainty factors, and a false-
applicability challenge. It deliberately lacks a locked holdout and full
falloff curves, so no supported-lane claim can be produced from it yet.

`stability/surface_resolution_profiles.json` separates software-regression
profiles from scientific validation. The retained tert-butoxy run currently
ends fail-closed at `intermediate_identity_unresolved`: its recovered stationary
geometry contains disconnected methyl and C3H6O fragments, so no rate or t95 is
released. Reproduce the same fresh four-core benchmark contract from the
repository root with:

```bash
conda run --no-capture-output -n orca-auto-py312 python -m storca stability \
  --smiles 'CC(C)(C)[O]' \
  --scenario low-pressure-intrinsic-gas-screen \
  --output-dir runs \
  --name tert-butoxy-surface-resolution \
  --charge 0 \
  --multiplicity 2 \
  --cores 4 \
  --method-profile b3lyp-def2-svp \
  --rmg-env rmg4_env \
  --temperature 298 \
  --pressure 1e-6 \
  --screen-tier quick-screen \
  --target-duration-hours 24 \
  --retention-fraction 0.95 \
  --verification-max-iterations 12 \
  --verification-timeout-seconds 14400 \
  --verification-orientations 3 \
  --verification-neb-images 8
```

This is a live ORCA/RMG calculation and may take hours. A newly generated
timestamped run must be assessed independently; the command does not promote
the retained abstention into a scientific validation reference.

## Spectrum benchmarks

Store a numerical reference spectrum as a two-column CSV: wavenumber in
cm^-1 followed by absorbance, transmittance, or relative intensity. Match that
same convention in the STORCA plot before comparison. Record sample phase,
solvent, concentration, temperature, and measurement method when running
`storca benchmark`.
