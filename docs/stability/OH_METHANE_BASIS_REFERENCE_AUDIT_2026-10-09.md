# OH + methane: bounded basis/reference audit — 2026-10-09

Status: partial, resource-blocked. The open-shell-reference comparison is complete; the quadruple-zeta barrier is unresolved. No rate, lifetime, calibrated uncertainty, or support-level promotion follows from this experiment.

## Result

All energies use the identical saved B3LYP/def2-SVP frequency-input geometries, charge, multiplicity and atom order. These are single points, not newly optimized higher-level stationary points.

| Electronic treatment | Electronic barrier, kJ/mol | Diagnostic barrier with retained low-level ZPE, kJ/mol |
| --- | ---: | ---: |
| Canonical CCSD(T)/def2-TZVPP, UHF/RHF | 30.676408 | 24.137529 |
| Canonical CCSD(T)/def2-TZVPP, QRO-derived open-shell reference / RHF methane | 30.232515 | 23.693635 |
| Canonical CCSD(T)/def2-QZVPP, UHF/RHF | Unresolved: no TS energy | Unresolved |

The signed reference barrier shift is -0.443893 kJ/mol. This particular reference choice therefore does not explain the approximately 20 kJ/mol fixed-geometry B3LYP-to-CCSD(T) barrier difference identified in the [preceding energy audit](../../runs/oh-methane-energy-audit-20261009-0ksMQL/REPORT.md). It does not rule out other electronic-state/reference problems and is not an accuracy certificate.

The retained ZPE contribution is -6.538880 kJ/mol. It cancels in the reference comparison. Showing it separately prevents the composite column from being mistaken for new high-level vibrational evidence.

The QRO/TZ shift is below the preregistered absolute 1 kJ/mol follow-up allocation threshold. QRO/QZ was therefore not dispatched. The threshold is not a confidence interval; basis/reference interaction remains untested.

Completed QZ species energies are OH −75.662498173999 Eh and CH4 −40.450705367065 Eh. Without a completed TS, neither a QZ barrier nor a TZ→QZ barrier shift can be calculated.

## Scientific scope

[Gruber and Czakó's published composite theoretical benchmark](https://www2.sci.u-szeged.hu/czako/papers/PCCP_22_14560_2020.pdf) reports classical and adiabatic abstraction barriers of 6.30 and 4.78 kcal/mol (26.3592 and 19.99952 kJ/mol). Those are theoretical barrier references, not measured activation barriers or experimental rates. Its high-level optimized geometries, explicitly correlated large-basis energies, post-CCSD(T), core, relativistic, spin–orbit and ZPE contributions differ from this fixed-geometry calculation.

Our QRO/TZ electronic and diagnostic composite differences from those respective references are +3.873315 and +3.694115 kJ/mol. These are descriptive cross-model differences, not isolated remaining error components or defensible bounds on the rate error.

QRO here means a UHF-derived quasi-restricted orbital reference used with native unrestricted canonical coupled-cluster equations; it is not an independently optimized ROHF calculation and not the paper's ROHF-based F12 convention. Actual QRO generation, orbital rereading and the QRO guess flag are present in both hashed native outputs. Preceding UHF SCF S² is not relabelled as the QRO reference's spin contamination. [ORCA MDCI documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/mdci.html)

Native normal termination, SCF and CC convergence, final CCSD(T) totals, residuals, chemical frozen cores and reference formalism are checked for every admitted point. These are calculation-admission checks, not physical validation. The retained native SCF policy uses its SOSCF/reasonable-density termination semantics; the report does not assert that every individually printed SCF threshold passed, reference stability was exhaustively tested, or a formal numerical error bound was established.

## Attempts and resource failure

Four distinct new points completed: QRO/TZ OH and TS, and UHF/RHF QZ OH and CH4. No completed point was dispatched twice.

1. [Initial reference attempt](../../runs/oh-methane-basis-reference-20261009-xtnOC1/audit-result.json): OH completed normally; its checker expected a status string not printed by this ORCA version. The rejection remains recorded. A fresh checker required three actual native QRO-use markers and reused the identical completed input/output.
2. [Four-core continuation](../../runs/oh-methane-basis-reference-20261009-continuation-WFYgpd/audit-result.json): TS/QZ failed in the exchange transformation. Native output reports a minimum 2180.9 MB per rank versus the supplied 1000 MB. No final energy was produced.
3. [Pre-dispatch guard attempt](../../runs/oh-methane-basis-reference-20261009-memory-retry-E5sJIU/audit-result.json): an overescaped comparison expression rejected the resource-only input change before any new native dispatch. The failed record is preserved.
4. [Resource-only continuation](../../runs/oh-methane-basis-reference-20261009-memory-continuation-uPugA9/audit-result.json): two MPI workers with 2500 MB per rank passed SCF and the prior memory gate. Only nprocs and MaxCore changed. The SCF energy agrees with the earlier four-worker attempt to about 3×10⁻¹⁴ Eh. The disk watchdog interrupted the owned process group at the 4 GiB free-space reserve, before CC convergence or a final energy.

All attempts have six-hour per-job limits and a 24-hour overall bound. Neither failure was a scientific time-limit exhaustion. MaxCore is not a hard bound on total resident memory. [ORCA memory guidance](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html#global-memory-use)

Native ORCA reports that KC_MO is unavailable for this UHF calculation and switches to stored AO exchange integrals. The documented lower-disk AO-direct variants likewise fall back for UHF. CDouble retains double-precision matrix-container storage with compression, but its benefit for the failing AOX scratch stream has not been established; it is not claimed as a demonstrated solution. No RI, local-pair, lower-precision or relaxed-threshold substitution was made. [MDCI documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/mdci.html), [storage documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/integralhandling.html)

## Preservation and cleanup

Original benchmark, rotor replay, previous energy audit and prior attempts were checked by content hashes and modification times before/after dispatch and remained unchanged. Core software, installed packages and production defaults were not modified.

After confirming all owned native processes had stopped, only ten explicitly listed large integral scratch files created by the interrupted continuation were removed: 8,347,374,756 bytes, approximately 7.77 GiB. These intermediates can be regenerated from the preserved input, but were not moved to trash and cannot be resumed from the removed files. No unrelated files were targeted. Inputs, outputs, wavefunction, process lease, native MDCI input and audit metadata remained hash-identical. The source manifest records names/sizes, not content hashes of deleted scratch.

- [Cleanup manifest](../../runs/oh-methane-basis-reference-20261009-memory-continuation-uPugA9/scratch-cleanup-manifest.json)
- [Cleanup result](../../runs/oh-methane-basis-reference-20261009-memory-continuation-uPugA9/scratch-cleanup-result.json)
- [Partial scientific comparison](../../runs/oh-methane-basis-reference-20261009-memory-continuation-uPugA9/partial-comparison.json)

## Next bounded work

First unblock the one missing canonical TS/QZ calculation with adequate scratch space: roughly 40–50 GiB free is a practical starting allocation, not a measured peak or guaranteed requirement. Reuse the completed OH/CH4/QRO points and retain the 4 GiB guard. An external SSD is an alternative; no additional mounted writable data volume was found. Do not delete user data autonomously.

After basis sensitivity is resolved, the next scientific stage is geometry/Hessian/ZPE and rate-construction assessment, not feature expansion:

- Refine the reactant and TS geometries under a declared independent treatment; assess the higher-level energies on those geometries separately from the current basis comparison. Verify any changed saddle and its connectivity before treating it as rate evidence.
- Check the native TS's 21.98 cm⁻¹ real mode (scaled to about 21.32 cm⁻¹ in the retained Arkane model) with tighter Hessian numerics/refined geometry. Its printed eigenvector suggests spectator OH orientation; a rotor assignment requires an explicit constrained scan. Do not delete or arbitrarily floor it to match a rate.
- Audit the OH electronic partition convention, including consistency of spin-free versus spin–orbit energy/population treatment. The retained Arkane conformer declares multiplicity two; orbital/fine-structure treatment remains unquantified. Do not add an arbitrary degeneracy factor or double-count a spin–orbit correction.
- Only then construct a consistent all-five-species energy/vibration model (OH, CH4, TS, H2O, CH3), including both barrier directions and a fresh tunneling/statistical-mechanics audit. Separate electronic single-point energy evidence from geometry/frequency evidence explicitly; no ZPE double-counting.
- Compare resulting rates over the declared temperature range and independently benchmark other reactions before changing production profiles or asserting calibrated accuracy.

The original approximately 5200-fold rate discrepancy has not been recomputed or fixed by this audit. Nothing here establishes a conclusive molecular stability prediction.
