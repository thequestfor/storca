# Expanded-network implementation checkpoint — 2026-10-08

The bounded implementation in
[the accepted plan](EXPANDED_NETWORK_IMPLEMENTATION_PLAN_2026-10-08.md) is complete.
This is software/native retained-model validation, **not** P1 scientific
acceptance or a defensible physical lifetime.

## What now works

- `retained_linear_network.py` preserves full native graph/electronic-state
  identities and the whole source-bound inventory. Same-formula singlet and
  triplet species cannot collide through suffix stripping.
- New native bridge operations retain actual Cantera forward/reverse constants.
  Fitted reverse Arrhenius inspection values are not silently used as the
  propagated rates. Nonlinear/collider, unbalanced, conflicting-environment,
  changed-source and incomplete-state cases refuse the linear diagnostic.
- `bounded_handoff.py` issues and recomputes a separate expanded-model
  assessment. The old certificate still refuses the expanded mechanism; its
  library, hashes and scope are unchanged. The comparison adds only a virtual
  direct product return, retaining the intermediate's actual return explicitly.
- `network_influence.py` retains baseline, branch-removal, return-removal and
  escape-removal responses. Native BDF/refinement/Radau integration checks tiny
  populations, integrated directed flux, conservation and positivity, with two
  tightened exact-horizon Cantera baseline checks. Failed checks stay unresolved.
- A pre-path state audit flags singlet/triplet channels without pretending an
  ordinary single-surface NEB/TST calculation will establish their rates.
- The production engine retains the current assessment, conditional profile,
  diagnostic responses and specific next-evidence requests. Failed assessments
  preserve the current full-model artifacts. Resume/ladder reporting retains
  these records and does not acquire verified observable credit.
- Same-input cache reuse reopens sources, sidecars and implementation hashes,
  recomputes scientific assessments/scenario interpretation, and makes no new
  native call. Changed inputs or implementation require a fresh assessment.

## Actual production-engine replay

[Replay 003](../../runs/p1-cyclobutene-expanded-assessment-003/result.json)
reopens the actual earlier stages and expanded mechanism through the production
verification engine. Path/rate/model callbacks each reopen retained evidence
once; they do not execute those programs. Native operations are exactly one
rate snapshot and one declared-model diagnostic. Cache reopening adds zero
calls. Elapsed time is **13.364 seconds**, with **zero ORCA, Arkane or RMG
discovery jobs**. All **1,321 original files / 114 directories** retain their
bytes and modification timestamps. Earlier replay 001 tests the direct
assessment wrapper; 002 is the preceding production checkpoint. They remain
historical and are not overwritten.

Conditions are unchanged: cyclobutene, 423.15 K, 0.0666611842 bar, pure dark gas,
isothermal batch, 3600 seconds. Seven species and all four directed channels
remain present, including singlet `C4H6(2)` and triplet `C4H6(3)`.

| Declared scenario | Parent fraction at one hour |
| --- | ---: |
| Full expanded baseline | 0.879853826295 |
| State-changing subnetwork removed | 0.880078734912 |
| Triplet return removed | 0.820861106266 |
| Triplet escape removed | 0.880078734910 |

The return-removal scenario changes retention by about **0.059 absolute
fraction**. This identifies dependence on an unresolved assumption, not a
physical uncertainty range or probability that return is absent. Both added
native reactions are flagged as spin-changing and ineligible for ordinary
single-surface path review.

The omitted direct-return bound is **2.2329939318100603e-7 absolute parent
fraction**, relative only to the two named fixed reverse models. The independent
native baseline difference is about **1.62e-12**, a numerical agreement check,
not a certified total numerical/physical error bound.

Operative forward rate is **3.5484417784723735e-5 s^-1**, while the retained
native fit evaluates to **3.548599062138614e-5 s^-1**. The approximately
4.43e-5 relative conversion difference passes the existing 1e-4 insertion
tolerance and is retained separately. No tolerance was weakened, and this
conversion difference is not absorbed into the reverse-omission bound.

[Production summary](../../runs/p1-cyclobutene-expanded-assessment-003/verification-result.json)
is `conditional_model_prediction`, with verified flux credit zero and physical
t95/lifetime null. Full conditions, source hashes, numerical passes, all species
profiles and next-evidence requests are in the replay's engine artifacts.

## Validation

The final application suite passes **1,332 tests in 71.653 seconds**, with
**24 expected skips**. The selected native `rmg4_env` suite passes **135 tests
in 3.114 seconds**, with no skips. The focused engine/reporting suite passes
20 tests, and roadmap traceability passes four tests. These check source and
condition refusals, electronic-state separation, numerical diagnostics, current
assessment handoff, cache reopening and preservation of unresolved claims.
They are software/native-model validation, not experimental accuracy or
uncertainty calibration. No physical acceptance gate is promoted by test counts.

## Reproduction

From the repository root in the existing STORCA application Python environment,
choose an output directory that does not already exist:

```sh
python tests/replay_expanded_forward_model.py --output runs/expanded-network-replay-new --rmg-env rmg4_env
```

This is a development replay, not a new stability CLI option or a new molecular
stability calculation. It does not supply experimental rates, calibration,
spin-crossing kinetics or pressure-adequacy evidence.

## Remaining scientific decisions

1. Resolve electronic-state/channel interpretation and supported return/escape
   constraints, or an applicable relevance bound. Specialized nonadiabatic rate
   evidence may be required; the current backend cannot supply it automatically.
2. Audit requested-pressure applicability before deciding on a two-well Arkane
   extension. A reported local pressure plateau alone is not an HPL error bound.
3. Resolve conformer/statistical/path counting before admitting a corrected rate
   recipe; do not apply the earlier partition-only factor-two diagnostic blindly.
4. Supply independently applicable rate uncertainty and observable-conditioned
   mechanism-coverage evidence, then test prediction robustness and benchmark
   accuracy. No development-case rate fitting or gate relaxation is authorized.

No new environment lane, discovery engine or electronic-structure backend was
added. The wider solution/liquid/solid/interface/exposure roadmap remains intact.
