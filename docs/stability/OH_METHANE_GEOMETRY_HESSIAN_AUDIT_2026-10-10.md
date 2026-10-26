# OH + methane: geometry and Hessian audit — 2026-10-10

Status: bounded calculations completed; independent first-order TS/route admission remains **unresolved** because analytical and numerical derivatives disagree. All prior artifacts are preserved. No production methods, original artifacts, rates, lifetimes, uncertainty calibration or environment support levels are replaced by this experiment.

Follow-up: the [completed direct energy/gradient curvature audit](OH_METHANE_CURVATURE_AUDIT_2026-10-10.md) supports positive local curvature along the disputed displacement, consistent with both numerical Hessians. That diagnoses one direction; it does not admit the full TS or supersede this stage's preserved unresolved result.

This follows the [completed fixed-geometry basis/reference audit](OH_METHANE_BASIS_REFERENCE_AUDIT_2026-10-10.md). Its question is narrower than overall rate accuracy: are the retained geometry and weak TS vibration reproducible, and can an independently refined saddle pass stationary-character and route checks?

## Retained B3LYP/def2-SVP controls

Three fresh Hessians used exactly the same saved TS geometry, with unscaled frequencies and projected translations/rotations. The changed numerical treatment uses ExtremeSCF, DefGrid3 and RIJONX/def2/J. This also changes exchange treatment relative to the older RIJCOSX source; differences must not be attributed solely to the SCF threshold.

| Same saved geometry | H-transfer imaginary mode, cm⁻¹ | Matched weak mode, cm⁻¹ | Raw positive-mode ZPE, kJ/mol |
| --- | ---: | ---: | ---: |
| Analytical | −1060.3483 | +26.5423 | 131.633735 |
| Central numerical, 0.005 bohr | −1060.2050 | +22.7657 | 131.612942 |
| Central numerical, 0.001 bohr | −1060.2324 | +22.7019 | 131.613107 |

The numerical weak-mode difference is 0.0638 cm⁻¹, versus 3.8403 cm⁻¹ between the analytical and smaller-step numerical results. The numerical modes match each other with overlap greater than 0.9999999. The analytical–numerical raw ZPE difference is only 0.02063 kJ/mol, but this does not establish entropy or rate-model accuracy. Reproducible weak curvature is not proof of a harmonic oscillator or a hindered rotor.

A separately refined B3LYP/SVP candidate passes native optimization and a fresh first-order H-transfer frequency check: −1059.8439 cm⁻¹, with a +23.3816 cm⁻¹ weak real mode. No fresh B3LYP IRC or refined-geometry numerical repeat ran. It remains diagnostic, not a newly verified rate path.

## Independent geometry and derivative checks

The independent treatment is native WB97X-D3/def2-TZVPP, with its built-in D3(0), not WB97X-D3BJ and not an additional dispersion correction. It is a declared alternative geometry/harmonic model, not DFT ground truth or a reproduction of the published explicitly correlated benchmark.

OH, CH4, H2O and CH3 all completed unconstrained optimization and separate analytical frequency calculations. Their declared states and covalent identities are retained, and their physical spectra contain no negative modes. OH has five projected external modes; the nonlinear species have six.

The TS optimization completed normally in 630.4 seconds and passed all five printed VeryTightOpt convergence criteria. Stationarity alone did not settle saddle order:

| Same refined WB97X-D3 geometry | H-transfer mode, cm⁻¹ | Matched angular mode, cm⁻¹ | Physical negatives |
| --- | ---: | ---: | ---: |
| Default analytical response | −1239.3492 | −31.6368 | 2 |
| Central numerical, 0.005 bohr | −1238.5081 | +76.1864 | 1 |
| Central numerical, 0.001 bohr | −1238.5267 | +76.1469 | 1 |
| Analytical, response residual 10⁻⁹ only | −1239.3472 | −31.6135 | 2 |
| Analytical, same residual plus response XC grid 6 / 4.959 | −1246.0294 | −31.2701 | 2 |

The angular mode changes sign despite a displacement overlap of approximately 0.99998. Both numerical Hessians agree: their weak-mode difference is 0.03947 cm⁻¹, overlap is 0.999999971, and raw positive-mode ZPE difference is 0.000412 kJ/mol. Neither analytical control resolves the discrepancy. The controlled analytical–small-step numerical mode gap remains 107.417 cm⁻¹. Therefore neither a physically consistent higher-order saddle nor an agreed first-order saddle is established. Numerical repetition is strong evidence of repeatability, not an automatic declaration that the numerical result is ground truth.

The weak mode is principally relative OH orientation about the noncovalent C–O axis, not ordinary O–H stretching. This characterization does not itself license a one-dimensional rotor treatment.

An independent read-only reconstruction of all five raw Cartesian Hessian matrices reproduces all 15 physical frequencies within 0.000219 cm⁻¹. The matrices are symmetric, have identical masses and agree in centered coordinates within 8.6×10⁻¹³ bohr. Projecting every matrix onto the **same** normalized Cartesian displacement gives opposite curvature signs:

| Hessian | Same-displacement curvature, Eh/bohr² |
| --- | ---: |
| Default analytical | −3.94510×10⁻⁵ |
| Numerical, 0.005 bohr | +2.30247×10⁻⁴ |
| Numerical, 0.001 bohr | +2.30008×10⁻⁴ |
| Tight analytical response residual | −3.93928×10⁻⁵ |
| Larger analytical response grid | −3.85009×10⁻⁵ |

Thus the sign conflict exists in the retained matrices themselves; it is not introduced by STORCA's frequency labels, mode matching or this independent translation/rotation projection. This does not identify the physical cause or rule out issues in native derivative construction. The [read-only matrix replay](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/matrix_replay.py) binds the closed result/Hessian hashes and makes this diagnostic reproducible without launching calculations or writing artifacts.

## Bounded control policy

ORCA's ground-state SCF tolerance and analytical CP-SCF response tolerance are separate controls. The native analytical calculation uses a response residual tolerance of 3×10⁻⁶ despite ExtremeSCF. The [documented CP-SCF controls](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/cpscf.html) permit an isolated response-residual test; [response integration grids](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/numericalintegration.html#changing-td-dft-cp-scf-and-hessian-grids) are controlled separately from the ground-state grid.

After the independent stage closed unresolved, an append-only diagnostic ran exactly two analytical Hessians at identical geometry, charge, multiplicity, functional, basis and ground-state grid:

1. Only tighten the response residual target to 10⁻⁹, retaining its 128-iteration limit.
2. If necessary, retain that target and increase only the response XC quadrature to angular grid 6 / radial accuracy 4.959.

The driver separates derivative consistency from first-order stationary character. Numerical allocation tolerances are mode overlap ≥0.90, no physical sign changes, weak-frequency difference ≤max(2 cm⁻¹, 5% of the reference magnitude), and total raw TS ZPE difference ≤0.1 kJ/mol. These are review thresholds, not calibrated physical error bars. Both controls completed all 21 response perturbations at iteration 18, with final residuals 4.6976×10⁻¹⁰ and 4.6611×10⁻¹⁰, below 10⁻⁹. Both reproduce the same ground-state energy and printed S² exactly. A consistent second negative would stop unresolved; settings are not tightened to force a favorable sign.

These observations exclude insufficient CP-SCF residual convergence as the principal discrepancy in this test. The tested response-grid change is also insufficient; they do **not** establish a grid-independent limit or rule out all integration, native derivative-construction or finite-difference issues. Native ORCA 6.1.0 outputs explicitly include the van der Waals Hessian term, so there is no evidence here for simply omitted dispersion. No RMG or Arkane calculation participates in this derivative disagreement.

Any changed saddle still requires a fresh bidirectional IRC with native endpoint convergence and a unique atom-mapped OH + CH4 / H2O + CH3 assignment. Old IRC evidence cannot verify another electronic-structure surface. XYZ endpoint connectivity does not independently certify electronic states, separated fragments or frequency-validated terminal minima. The prepared conditional IRC driver was exercised in read-only preflight and correctly refused the unresolved derivative evidence **before writes or native dispatch**. No new IRC ran.

The prepared five-point constrained orientation diagnostic is not dispatched before those gates. Its fixed transfer distances would suppress coupling to reaction progress; even a completed local profile could not supply a periodic rotor partition function or a rate correction. [ORCA's constraint documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html#constrained-optimization) explicitly distinguishes constrained optimization from unconstrained stationary-point evidence.

## Execution and remaining claims

All new native jobs are sequential, four MPI workers, MaxCore 1000 MB per worker, with a six-hour per-job allowance and a 4 GiB runtime disk reserve. Numerical Hessians take about ten minutes; they are not subject to the historical ten-minute cutoff. Fresh folders/basenames, input hashes, owned process leases and native normal completion are required.

The B3LYP control, independent-stage and response-control receipts confirm unchanged original benchmark, Arkane replay, earlier energy audit and completed QZ audit artifacts. Bound source hashes and modification times also pass the response-stage check. Nineteen native jobs completed normally: five B3LYP controls, twelve independent optimization/frequency/repeat jobs and two analytical response controls. Execution success is not scientific admission: the independent stage exits unresolved and the response stage records `local_stationary_and_numerical_admission=false`. No original scientific files are overwritten and no additional user data are deleted. About 58 GiB remained free after the last native calculation.

There is still no new canonical energy barrier at the refined geometries, consistent five-species rate construction, OH electronic-partition/spin–orbit audit, rotor/free-energy model, tunneling validation or experimental-condition match. The original approximately 5200-fold rate discrepancy has not been recomputed or demonstrated fixed.

Evidence:

- [B3LYP controls and preservation receipt](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/control-result.json)
- [Independent source/method/budget contract](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/independent-contract.json)
- [Closed independent derivative result](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/independent-result.json)
- [Bound response-control contract](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/response-control-contract.json)
- [Response outcomes, mode comparisons and preservation receipt](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/response-control-result.json)
- [Larger-grid native analytical output](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/response-control/residual-and-XC-grid/freq.out)

## Next smallest discriminating experiment

Before more geometry refinement, rotor assignments or rates, directly test the disputed mode's curvature using independent symmetric energy/gradient points at two small amplitudes along the same fixed-frame displacement. Compare central energy curvature and gradient-projected curvature with Rayleigh projections of both full Hessians. Record linear-gradient drift, SCF solution identity, quadrature and step sensitivity. This tests energy→gradient→Hessian consistency without relying solely on either frequency workflow or using constrained relaxation to change the question.

If the disagreement persists, isolate derivative treatment or native-version effects with a minimal reproducer. The installed executable is 6.1.0. The official [6.1.1 change log](https://www.faccts.de/docs/orca/6.1/manual/contents/appendix/detailedchangelog.html) lists Hessian fixes, but does not identify this exact sign discrepancy; there is no evidence here that an upgrade would fix it. No package upgrade or upstream report was performed.

Only after local derivatives agree should fresh IRC verification and the local orientation/free-energy assessment proceed. A periodic rotor model, consistent five-species high-level energies/vibrations and OH electronic-partition/tunneling conventions remain separately scoped requirements for a defensible revised rate.
