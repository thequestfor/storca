# Specialist rate backend decision record

## Decision status

The integration boundary and development benchmark exist. No backend is yet
scientifically accepted as a default STORCA producer. The current decision is
to evaluate two separable capabilities rather than selecting one program by
name:

1. **VaReCoF-class VRC-TST production** for barrierless fragment association;
2. **MESS-class master-equation propagation** for pressure-dependent wells and
   channels after microcanonical rate evidence is available.

This ordering does not replace Arkane. Arkane remains the working adapter for
ordinary verified transition-state rates and its currently supported
pressure-dependent preparation. A new backend is admitted only where it adds a
physical capability and passes the same evidence and repaired-network gates.

## Why the capabilities are separate

A capture calculation determines flux through a variational dividing surface
for interacting fragments. A master-equation calculation propagates
energy-resolved populations among wells and channels in a collision bath.
They require different evidence. A solver name alone cannot establish either
claim, and a pressure-dependent solver does not create a defensible capture
rate merely because it can consume one.

## Candidate assessment

| Candidate | Officially documented capability | Fit for STORCA | Current decision |
|---|---|---|---|
| [VaReCoF](https://tcg.cse.anl.gov/papr/codes/varecof.html) | Variable-reaction-coordinate TST with fragment pivot surfaces and analytic or direct electronic-structure interaction energies | Best scientific match for loose radical association. Its documented direct interfaces are Gaussian and MOLPRO, not ORCA, so an ORCA-native claim requires a separately validated potential/sampling bridge. | Primary capture-method research target; external retained-result adapter first, native automation only after benchmark evidence. |
| [MESS](https://github.com/Auto-Mech/MESS) | Open-source one-dimensional master-equation solution for complex-forming reactions and pressure-dependent phenomenological rates | Strong fit for pressure-dependent propagation, installable from the Auto-Mech conda channel, and separable from capture-rate generation. | Primary external master-equation integration target after input/output fixtures and balance checks pass. |
| [MultiWell](https://multiwell.engin.umich.edu/about/) | Stochastic master equation plus canonical and microcanonical VTST utilities along supplied paths | Valuable independent comparison and possible fallback, but it is a suite rather than a single automated evidence contract. | Benchmark comparator; do not integrate before MESS/VaReCoF feasibility is resolved. |
| AutoMech/MechDriver | Workflow automation around thermochemistry and kinetics, with archived-data and cluster-oriented testing | Potential source of reusable MESS/VaReCoF preparation logic. Current documented fresh electronic-structure workflow expects Gaussian 16 and development setup is comparatively heavy. | Inspect adapters; do not make it a mandatory STORCA runtime dependency. |
| Arkane | Existing STORCA central TST and selected pressure-dependent preparation | Already integrated and valuable. It does not currently supply STORCA with independently validated barrierless capture evidence. | Preserve; tighten inputs and uncertainty instead of replacing it. |

## Provider acceptance contract

A capture result must retain all of the following before its generated RMG
library can enter network repair:

- exact condition, path-evidence, reaction, charge, and multiplicity hashes;
- an explicit long-range fragment interaction model and parameter source;
- sampled dividing surfaces and the variational criterion;
- a recrossing model and evidence;
- electronic degeneracy, symmetry, and electronic-state applicability;
- the temperature/pressure domain and bath-gas/collision model when relevant;
- calibrated uncertainty components;
- immutable raw artifacts and checksummed generated kinetics library.

The runtime validator now rejects `capture_theory_with_recrossing` results that
omit these physical fields. It also rejects a VTST result without a reaction
coordinate, sampled surfaces, and variational criterion.

## Benchmark gate

`benchmarks/kinetics/specialist_rate_manifest.json` contains development
limiting-rate anchors for H + O2, OH + OH, CH3 + O2, and C2H5 + O2 association,
plus an electronic-state ambiguity challenge. The numerical anchors are from
the cited JPL evaluated kinetic recommendation; H + O2 also points to the
primary multi-bath-gas experiment. These are not ground truth and experimental
agreement alone does not prove a mechanism.

Run a provider result set with:

```bash
storca specialist-rate-benchmark \
  benchmarks/kinetics/specialist_rate_manifest.json \
  path/to/provider-results
```

The current manifest is deliberately `development_reference_set`: it has no
locked holdout and therefore cannot emit `support_claim_allowed: true`, even
if every supplied development prediction agrees. Promotion requires digitized
falloff curves, locked held-out systems, correlated reference uncertainties,
compute accounting, and zero false applicability on the challenge set.

## Next implementation slice

1. ~~Add non-executing MESS request rendering and parsed-output fixtures.~~
   Implemented in `storca/mess_adapter.py`; it accepts only a hash-bound,
   source-backed one-well/two-fragment/Rotd model and makes no rate claim.
2. ~~Add bounded MESS execution and exact native-table parsing.~~ Implemented
   with immutable input revalidation, executable-binary hashing, timeout and
   failure retention, exact channel/grid checks, and checkpoint reuse. The
   resulting table remains explicitly unadmitted.
3. ~~Enforce detailed balance, native-to-RMG unit conversion, uncertainty, and
   generated-library identity before producing a specialist result.~~
   Implemented with full-grid reverse-channel, collision, and serialized-library
   checks. No native bimolecular amount basis is assumed.
4. Reproduce at least one development falloff case without tuning to it. The
   exact-pressure curve metric, holdout-lock protocol, and H + O2/N2 pending
   case are implemented; independent numerical reference extraction and a live
   MESS result remain open.
5. Add and lock independent holdout cases before changing acceptance thresholds.
6. In parallel, produce one VaReCoF retained-result bundle for a simple radical
   association and test the capture-specific admission fields.

Do not install either backend as a silent default and do not translate a
collision ceiling into a capture rate.
