# OH + methane: completed basis/reference audit — 2026-10-10

Status: completed fixed-geometry sensitivity experiment. The missing canonical CCSD(T)/def2-QZVPP transition-state single point completed normally. No rate, lifetime, calibrated uncertainty, production-profile change or support-level promotion follows from this audit.

This resolves the resource blocker in the [dated partial checkpoint](OH_METHANE_BASIS_REFERENCE_AUDIT_2026-10-09.md), which remains unchanged as historical evidence.

## Scientific result

All points use identical saved B3LYP/def2-SVP frequency-input geometries, charge, multiplicity and atom order. These are higher-level single points, not reoptimized higher-level minima or saddles.

| Electronic treatment | Electronic barrier, kJ/mol | Diagnostic barrier with retained low-level ZPE, kJ/mol |
| --- | ---: | ---: |
| Canonical CCSD(T)/def2-TZVPP, UHF/RHF | 30.676408 | 24.137529 |
| Canonical CCSD(T)/def2-TZVPP, QRO-derived open-shell reference / RHF methane | 30.232515 | 23.693635 |
| Canonical CCSD(T)/def2-QZVPP, UHF/RHF | 27.508540 | 20.969660 |

The UHF/RHF TZ→QZ barrier shift is **−3.167869 kJ/mol**. Basis sensitivity is therefore material for this particular calculation; the TZ barrier alone was not sufficiently basis-insensitive. Two basis sizes do not establish the complete-basis limit or a calibrated remaining error bound.

The TZ open-shell-reference shift remains **−0.443893 kJ/mol**. It is smaller than the basis shift here, but does not establish reference independence in general. The preregistered 1 kJ/mol reference-follow-up allocation threshold was not reached, so QRO/QZ was not dispatched. Basis/reference interaction is explicitly **not tested**; QRO here is not optimized ROHF or the published benchmark's ROHF-F12 convention.

The unchanged low-level ZPE contribution is **−6.538880 kJ/mol**. It cancels from the basis/reference shifts. The composite column must not be interpreted as a new high-level vibrational calculation.

The admitted QZ energies are:

| Point | Energy, Eh |
| --- | ---: |
| OH | −75.662498173999 |
| CH4 | −40.450705367065 |
| TS | −116.102726091553 |

The barrier is `(E_TS − E_OH − E_CH4) × 2625.4996394799 kJ/mol per Eh`.

## Comparison, not an accuracy certificate

[Gruber and Czakó's composite theoretical benchmark](https://www2.sci.u-szeged.hu/czako/papers/PCCP_22_14560_2020.pdf) reports classical and adiabatic abstraction barriers of 6.30 and 4.78 kcal/mol, respectively: 26.3592 and 19.99952 kJ/mol. These are theoretical barrier references, not measured activation barriers or experimental rates.

The present QZ electronic and diagnostic composite barriers differ from those respective values by **+1.149340** and **+0.970140 kJ/mol**. These are descriptive cross-model differences, not isolated electronic-structure errors or bounds on rate error. The published treatment has different optimized geometries, explicitly correlated large-basis energies, post-CCSD(T), core, relativistic, spin–orbit and ZPE contributions. Agreement can include cancellation of unlike approximations.

No new Arkane, RMG or kinetic propagation ran. The original approximately 5200-fold rate discrepancy has **not** been recomputed or demonstrated fixed. Correcting a barrier alone does not validate the statistical-mechanical prefactor, soft modes, electronic partition convention, tunneling, reverse barrier or experimental-condition correspondence.

## Execution and admission evidence

The user-authorized scratch retry dispatched exactly one new point: TS/QZ, two MPI workers, MaxCore 2500 MB per worker, six-hour per-job allowance and a 4 GiB runtime free-space reserve. Its input is byte-identical to the preceding disk-interrupted two-worker retry. Relative to the original four-worker attempt, only worker count and MaxCore differ. No RI, local-pair, lower-precision or relaxed-threshold substitution was introduced.

The new native job took **669.06 seconds**, about 11.15 minutes. It passed the previous exchange-transformation memory gate, completed SCF, CC and triples, and terminated normally with return code zero. The final CC residual is **9.44×10⁻⁷**, below the declared 10⁻⁶ threshold. The TS uses the observed UHF reference, four chemically frozen-core electrons and no QRO transformation. Its native T1 diagnostic is 0.021104678, preceding SCF S² is 0.789739, and S diagnostics are 0.219641 / 0.239185 / 0.514713. These are retained diagnostics, not proof that single-reference CCSD(T) is exact or universally applicable.

Native SCF convergence retains ORCA's SOSCF/reasonable-density semantics. Not every individually printed SCF threshold passed; no exhaustive reference-stability test or formal numerical error bound is claimed.

Initial free space was about 59 GiB. Sampled free space during the disk-heavy work was about 35 GiB, and the guard did not trigger. These are whole-volume observations, not a measured scratch peak or a universal resource guarantee. ORCA's normal completion removed its large temporary intermediates: the new attempt occupies about 7.3 MiB, and about 57 GiB was free at the post-run check. No additional user files were deleted.

Four previously completed reference/QZ points were reused with identical input/output hashes, not recalculated. Original benchmark, Arkane replay, earlier energy audit and all prior attempts retain identical content hashes and modification times. Acceptance requires successful driver exit **and** `prior_scientific_artifacts_unchanged=true`, not `status=completed` alone. Both passed. No owned native ORCA/MPI process remained at the final check.

- [Bound contract and reuse provenance](../../runs/oh-methane-basis-reference-20261010-scratch-retry-GZelcb/contract.json)
- [Completed comparison](../../runs/oh-methane-basis-reference-20261010-scratch-retry-GZelcb/comparison.json)
- [Audit result and preservation receipt](../../runs/oh-methane-basis-reference-20261010-scratch-retry-GZelcb/audit-result.json)
- [Native TS output](../../runs/oh-methane-basis-reference-20261010-scratch-retry-GZelcb/uhf-qzvpp/TS/energy.out)
- [TS admission diagnostics and hashes](../../runs/oh-methane-basis-reference-20261010-scratch-retry-GZelcb/uhf-qzvpp/TS/result.json)

## Next bounded scientific task

Basis sensitivity is now measured, not exhausted. The next distinct task remains geometry/Hessian/ZPE and rate-construction assessment, as scoped in the previous checkpoint:

1. Refine reactant/TS geometry under a declared independent treatment; separately assess geometry sensitivity and validate any changed saddle/connectivity.
2. Resolve the retained TS's 21.98 cm⁻¹ real mode with tighter Hessian/refined geometry evidence; only assign a rotor after an explicit constrained scan. No arbitrary frequency floor or mode deletion to fit a rate.
3. Audit OH electronic partition and spin–orbit energy/population conventions without arbitrary degeneracy factors or double-counting.
4. Construct a consistent five-species energy/vibration model, including H2O and CH3, both barrier directions and fresh tunneling/statistical-mechanics evidence, before recalculating rates.

This completed calculation does not authorize those new jobs, change defaults or close the broader stability/chemical-fate validation gates.
