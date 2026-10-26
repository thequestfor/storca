# OH + methane: NoRI control and numerical-Hessian route audit — 2026-10-10

Status: the six-job conventional-integral control is complete. Its stored
analytical Hessian still disagrees with its own direct energy/gradient evidence.
The separately admitted, original-model numerical-Hessian IRC completed in
both directions and verified the mapped OH + CH4 / H2O + CH3 connection.
No revised rate or lifetime follows from route verification alone. Production
method defaults remain unchanged.

## Why this continuation is justified

The [direct-curvature audit](OH_METHANE_CURVATURE_AUDIT_2026-10-10.md) supplies
two-step independent energy and gradient support for the two repeated full
numerical Hessians at the retained RIJONX geometry. The analytical discrepancy
is retained; no imaginary mode is deleted or turned into a rotor.

The saved native optimization passes all five final VeryTightOpt criteria.
Its fresh raw Cartesian gradient has RMS 6.604e-7 and maximum 1.46549e-6
Eh/bohr. Mass-metric translation/rotation projection gives RMS 6.979e-7 and
maximum 1.47492e-6. Both pass the explicit RMS2e-6/MAX3e-6 stationarity review.
This review is an allocation/admission threshold, not a physical accuracy bound.
There is no evidence requiring another geometry optimization in this stage.

The derivative admission additionally replays all fifteen physical-mode
matches between full central NumFreq increments 0.005/0.001 bohr, their signs,
the unique opposing C–H-breaking/O–H-forming imaginary displacement, soft-mode
overlap and ZPE repeat, and the original direct energy/gradient calculations.
Electronic-state/population screens remain diagnostics, not proof of a unique
wavefunction, spin ground state or accuracy against experiment. Every Cartesian
NumFreq child wavefunction has not been independently certified.

## Conventional-integral control

The control removes RIJONX and auxiliary def2/J and uses `NoCOSX NoRI` plus
explicit `%method RI off`. Functional, orbital basis, original Cartesian
geometry, ground-state grid, charge, multiplicity and D3(0) remain fixed.
Response controls match the earlier residual-only analytical test. The
[ORCA RI documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/RI.html)
distinguishes disabling COSX from disabling fitting; requesting only NoCOSX
would not be a conventional-Coulomb control.

Actual native output confirms absent Coulomb and J/K fitting bases, completed
exact two-electron integral Hessian, and exact SHARK two-electron gradients.
All twenty-one response perturbations complete at iteration 18 with maximum
residual 4.7209e-10, below requested 1e-9. The original control reader incorrectly
required a nonexistent `Coulomb Hessian` output phrase; it failed closed after
the frequency job terminated normally. Its driver, contract and refusal remain
unchanged. An append-only reader replay reuses that completed Hessian and runs
only the five already planned EnGrad points. No duplicate frequency job runs.

For the same Cartesian-normalized original Num001 weak vector, units Eh/bohr²:

| NoRI evidence | Energy curvature | Gradient curvature |
| --- | ---: | ---: |
| ±0.01 bohr | +1.602500e-4 | +1.627480e-4 |
| ±0.02 bohr | +1.677125e-4 | +1.776599e-4 |
| Two-size local extrapolation | +1.577625e-4 | +1.577774e-4 |
| Stored analytical matrix, uᵀHu | −1.104510e-4 | — |

The smooth-leading-h² extrapolation is descriptive and uncalibrated, not an
independently established asymptotic limit. Raw estimates pass the original
max(5e-6 Eh/bohr², 10%) review threshold. Analytical/direct agreement fails.
Native analytical frequencies retain −1238.5066 and −52.8384 cm⁻¹. All six
control states print S²=0.758948; displacement charge/spin changes are at most
4e-5. Simple mass-metric T/R projection changes the directional stored-matrix
curvature by less than 3e-16 Eh/bohr², not enough to explain the discrepancy.

Important limitation: this geometry is **not stationary on the changed NoRI
surface**. Its raw gradient RMS is 4.694e-5, maximum 8.993e-5 Eh/bohr.
Do not call its two negative frequencies a validated second-order saddle or
mix its energies into the original rate model. The
[frequency manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/frequencies.html)
also distinguishes acoustic-sum-rule Hessian adjustment from subsequent T/R
projection; the control does not isolate every native intermediate matrix.
The defensible conclusion is that removing fitting did not cure the observed
analytical negative or stored-matrix/direct discrepancy in this control. No
unique defective component, native bug or electronic-structure accuracy claim
is established.

The native frequency job took 120.3 seconds; the five remaining checks took
123.1 seconds. Independent audit verifies 144 bound sources and 991 prior
folder files unchanged, including the failed reader attempt and reused Hessian.

## Original-model IRC stage

This is one fresh bidirectional IRC on the original neutral-doublet
UKS WB97X-D3/def2-TZVPP, D3(0), RIJONX/def2/J, ExtremeSCF, DefGrid3 surface.
No NoRI energies or Hessian are used. It reads the admitted original Num001
Hessian with `%irc InitHess read` and its exact `Hess_Filename`, follows physical
mode 0, and allows 100 iterations per direction. This is the
[documented IRC read-seed interface](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/irc.html).
Native read type, actual filename, mode and completed read are checked, not
merely the requested input. Analytical/frequency reentry is forbidden.

The stage is restricted to one job, four MPI workers, MaxCore1000 MB/worker,
six hours total/job and a 4 GiB disk reserve. It cannot refine geometry,
recalculate Hessians, launch Arkane, replace energies/ZPE or claim a lifetime.
All original scientific folders and bound source bytes/mtimes are protected.
The initial native IRC energy reproduces the original center energy,
−116.259386991363 Eh, and the header confirms the retained Num001 seed read.

Completion requires native convergence in both directions, complete terminal
artifacts and exactly one mapped OH+CH4/H2O+CH3 side assignment. The explicit
map preserves O0–H1, breaks C2–H6 and forms O0–H6; endpoint reference RMSD is
not a gate. Bond graphs alone do not certify fragment electronic states, and
IRC terminal complexes need not be separated frequency-validated minima.

The [closed numerical-route result](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/numerical-continuation-result.json)
reports native forward convergence after 28 steps to OH + CH4 and backward
convergence after 27 steps to H2O + CH3. The mapped hydrogen transfer and
native numerical-seed read pass; prior scientific artifacts and bound source
bytes/mtimes remain unchanged. The analytical-Hessian discrepancy remains
unresolved. This establishes a route on the declared DFT surface, not a
calibrated rate or general STORCA stability claim. The approximately 5200-fold
historical rate discrepancy is not resolved by this IRC.

## Evidence and next boundary

- [Preserved first NoRI result](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/nori-control-result.json)
- [Append-only reader and pure frequency replay](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/nori_replay.py)
- [Closed NoRI replay and direct data](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/nori-control-replay-result.json)
- [IRC-only continuation driver and read-only tests](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/numerical_continue.py)
- [Source-bound IRC contract](../../runs/oh-methane-geometry-frequency-20261010-lX9ABY/numerical-continuation-contract.json)

With the route result closed, the next task is consistent five-species
energy/statistical-mechanical/rate construction at an explicitly declared
recipe. Old fixed-geometry correlated energies cannot silently become energies
at these newly refined structures. Soft-mode treatment, tunneling, electronic
counting and experimental-condition matching remain separate rate questions;
agreement with one reference would not establish the mechanism or calibration.
