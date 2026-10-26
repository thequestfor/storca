# OH/methane: consistent electronic-state replay

Status: **step 1 implemented; conditional native diagnostic completed**.
No new ORCA calculations, production defaults, repaired-network rates,
calibrated uncertainty or environment-support claims were introduced.
The next step is the bounded angular/basin audit, not selection of the
best-agreeing branch.

## Implementation

The new [electronic-state module](../../storca/electronic_states.py) provides
immutable, strictly serialized level/state declarations with deterministic
hashes. Level degeneracies include spin already; the evaluator divides out
the native spin factor before replacing it. It supplies:

- Ground-energy offsets and temperature-dependent partitions, with explicit
  energy reference, source pointers and assumptions.
- Electronic changes to G, H, S and Cp relative to the native spin-only model.
- Both barrier shifts, forward/reverse population factors, classical rate
  factors and an equilibrium/detailed-balance consistency relation.
- Log-first arithmetic; unrepresentable linear values are explicitly absent,
  not silently zero/infinite. No native jobs or applicability admission.

The new [experiment replay](../../tests/replay_oh_methane_electronic_states.py)
uses the unchanged [completed five-point audit](OH_METHANE_CONSISTENT_RATE_AUDIT_2026-10-11.md).
It reconstructs only typed native translation/rotation/vibration objects from
the hash-bound statmech snapshots; it never evaluates their repr strings or
executes retained Arkane input. The selected RMG environment supplies native
TST and new Wigner/Eckart objects. Both the matched DFT-only and composite
recipes are evaluated at all 21 retained temperatures, including all fourteen
measurement rows. Experimental rates enter a separate comparison after replay,
not the physical model or a fit.

It reproduces the original spin-only rates, checks equivalent ground/centroid
energy references, and checks joint count-one/count-two branches. The parent
also recomputes rate/partition/barrier/thermo witnesses instead of trusting a
stored success flag. Sources, runtime/module hashes and native constants must
match before comparison. Missing/changed source evidence fails closed.

## Declared physical model

For OH, use two electronic levels each with degeneracy two, separated by
140 cm⁻¹. The explicit level/partition convention is supported by
[Chakraborty–Truhlar, Eq. 29](https://truhlar.chem.umn.edu/sites/truhlar.chem.umn.edu/files/2022-05/final704.pdf)
and the [OH/methane kinetics study](https://arxiv.org/pdf/1503.01257).
The spin-free centroid interpretation and stationary-point SOC quenching are
motivated by [Gruber–Czakó](https://www2.sci.u-szeged.hu/czako/papers/PCCP_22_14560_2020.pdf).
Their computed corrections are not calculations on STORCA's geometry.
Here TS quenching, thermal populations and a separable electronic/rotational
model are **assumed**, not demonstrated.

With splitting Δ and ground lowering δ=Δ/2:

```text
OH electronic partition relative to ground: 2(1 + exp(-Δ/RT))
Native OH Q already includes spin two: replace it by factor (1 + exp(-Δ/RT))
Derived OH ground E0: native E0 − δ
Forward barrier: native barrier + δ
Reverse barrier: unchanged
Reaction E0: native reaction E0 + δ
```

No new TS spin factor or methane ×4 channel factor is added. The imaginary
frequency and all positive harmonic modes stay unchanged. Ground and centroid
declarations give the same physical result when their energy references are
transformed together.

Δ=1674.771919 J/mol uses exact-SI h, c and NA. Native thermodynamic exponents
use the retained RMG R=8.314472 J/mol/K, not a silently substituted modern
value. The readable constants source and actual compiled-runtime values are
both checked. This convention is recorded, not interpreted as physical error.

The conditional composite ground barriers are 22.889064 kJ/mol forward and
80.820220 kJ/mol reverse. A new asymmetric Eckart object uses both asymptotes;
its factor at 340.3 K changes from 3.26184 to 3.29788. Multiplying the old
Eckart rate by a population-only correction would not be this calculation.

Corrected equilibrium is evaluated independently from partitions and from
original native **Species** Gibbs energies plus the electronic ΔG. A native
spin-only Gibbs value after an E0 shift is labelled separately; it does not
contain the extra electronic-level entropy. No inconsistent thermo is exported.

## Results, not branch selection

At the frozen development study's 340.3 K point, the measured coefficient is
8430.997064 m³ mol⁻¹ s⁻¹. Composite direct Eckart results:

| Joint declared model | k, m³ mol⁻¹ s⁻¹ | k / measured |
| --- | ---: | ---: |
| Original spin-only, one harmonic TS basin | 25,818.37 | 3.0623 |
| OH two-level model, one harmonic TS basin | 12,500.34 | 1.4827 |
| OH two-level model, two independent mirror harmonic basins assumed | 25,000.68 | 2.9653 |

Across all fourteen measurement rows, the conditional OH/Eckart ratio is
1.4035–1.5834 for count one and 2.8071–3.1668 for count two. Both are retained.
Geometric chirality does not establish independent thermal harmonic basins.
Count two changes TS Q and **both** directions by exactly two while leaving
barriers, tunneling and equilibrium unchanged. It is not a separately fitted
electronic degeneracy correction.

Classical TST and Wigner are also reported, not silently chosen when they
agree better. Their composite count-one coefficients at 340.3 K are 3790.42
and 8121.05 m³ mol⁻¹ s⁻¹. Neither establishes that the corresponding dynamics
is physically adequate. The maximum native TST/partition/Gibbs identity
residual over both stages/all branches/temperatures is 6.62×10⁻¹².

## Evidence and preservation

The completed fresh run is
[oh-methane-state-replay-20261011-y_sxkygn](../../runs/oh-methane-state-replay-20261011-y_sxkygn/result.json).
Its native command took 4.74 seconds, including the selected environment's
startup. There are **zero electronic-structure jobs**.

- [Frozen contract](../../runs/oh-methane-state-replay-20261011-y_sxkygn/contract.json): 125 declared file bindings,
  33 Python source snapshots and an exact reconstructed original composite
  generator using the existing one-line revision proof. Native binaries are
  hash bound, not copied as redundant large snapshots.
- [Native replay](../../runs/oh-methane-state-replay-20261011-y_sxkygn/native-replay.json): direct rates,
  partitions, ground barriers, corrected thermochemistry, constants and checks.
- [Execution receipt](../../runs/oh-methane-state-replay-20261011-y_sxkygn/execution.json): actual command,
  selected environment, streams and measured time.
- [Conditional comparison](../../runs/oh-methane-state-replay-20261011-y_sxkygn/result.json): every declared
  branch and all fourteen measurement rows; no fitting or promotion.

All original completed scientific files retain their hashes/mtimes. New
output is exclusively created in a fresh directory. Historical contracts,
native logs/Hessians, original helper/generator sources, package installations
and production defaults are unchanged. The saved snapshots preserve source
versions; they are not permission to accept a changed live implementation.

Two development attempts are also preserved: a pre-native constants reader
initially treated a compiled module as text, and a subsequent replay correctly
refused an equilibrium witness because native Conformer Gibbs excludes E0
while Species Gibbs includes it. Both are STORCA diagnostic-adapter mistakes,
not RMG physical failures; neither produced an admitted rate. The final driver
handles compiled/readable constant bindings, uses full Species Gibbs, and
records contract-building failures without dispatch.

## Run and verification

From the repository using the normal STORCA/application Python environment:

```sh
python -B tests/replay_oh_methane_electronic_states.py runs/oh-methane-consistent-rate-20261011-XqYhq3 --rmg-env rmg4_env
```

This creates another fresh diagnostic, never overwrites the completed run,
and launches no ORCA calculations. A different environment is accepted only
if its interpreter/native bindings match the retained recipe; drift is an
explicit refusal, not an automatic package repair.

Focused tests: 22 pass, including analytic populations, energy-reference
invariance, both-direction detailed balance, thermodynamic derivatives,
strict serialization, source/hash/mtime refusal, semantic k/Q/K/barrier
mutation, runtime drift, false admission, exclusive writes, preflight/freeze
failure and forbidden ORCA/Arkane dispatch. Full-suite closeout is recorded
as 1,457 tests in 182.137 seconds, with zero failures/errors and 24 skips.
That run collected the initial nineteen new tests; the three subsequent
runtime/constants/freeze-failure checks also pass. A final 47-test run includes
all twenty-two new tests and the composite/revision/comparison/roadmap
regressions. These are software checks, not physical accuracy certification.

## Remaining gate

The result is `diagnostic_not_admitted` and every accuracy/calibration/network/
lifetime promotion flag remains false. Fourteen temperatures are still one
development study. Humid-air versus isolated-gas applicability, finite-basis
electronic error, TS angular/orthogonal coupling, harmonic validity and
tunneling/recrossing remain unresolved.

Next: the capped transverse angular-saddle/basin audit in
[the bounded plan](BOUNDED_RATE_VALIDATION_PLAN_2026-10-11.md), then recipe
freeze and independent rate transfer tests. This increment is not general
stability prediction and does not narrow the project's multi-environment goal.
