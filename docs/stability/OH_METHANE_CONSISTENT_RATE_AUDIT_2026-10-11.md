# OH + methane: consistent five-species rate reconstruction

Status: completed native DFT and matched five-point CCSD(T)/def2-QZVPP rate
reconstruction, with hash-reopening closeout and all 14 frozen measurement
temperatures evaluated. The composite Eckart rate remains 2.92–3.27× the
measured comparator. This is an isolated development comparison, not
production-profile promotion, calibrated uncertainty or a lifetime claim.

## What changed

The numerical-Hessian IRC already verifies the intended mapped OH + CH4 →
H2O + CH3 connection on the declared WB97X-D3/def2-TZVPP surface. This task
reuses those refined structures, the four admitted minimum frequency logs and
the repeated/admitted Num001 TS frequency log. No optimization, NEB, Hessian,
IRC or discovery job is repeated.

The new [composite recipe API](../../storca/composite_rate_recipe.py) separates
hash-bound electronic-energy input/output from geometry/frequency input/output
and Hessian. The two sources must describe the same Cartesian nuclear order,
coordinates, charge and multiplicity. A uniform native Hessian-origin
translation is ignored, not a rotation, permutation or geometry change.
Raw ORCA logs/Hessians are staged byte-for-byte. Electronic energies are
explicit scalar `e_electronic` values under the existing balanced-per-atom
gauge; they are not disguised by editing frequency-log energies.

Arkane's software-qualified `LevelOfTheory`/`CompositeLevelOfTheory` objects
explicitly name ORCA. A first string-only omegaB97X-D3 declaration failed
before any ORCA dispatch because native Arkane requires its software identity.
The [failed attempt](../../runs/oh-methane-composite-rate-20261011-5ivCKk/baseline/execution-result.json)
is unchanged; a fresh corrected recipe executes in `rmg4_env`. No installed
RMG/Arkane files or production method defaults are edited.

The API validates source consistency, not the IRC's scientific meaning,
mechanism coverage, ground-state ordering or rate accuracy. Its artifacts
explicitly deny automatic network reinjection. The existing rotor contract
wire format is reused, not a new route certificate.

## Declared model

- Isolated ideal-gas single-adiabatic-surface RRHO/TST abstraction.
- Frequency/geometry: original admitted UKS/RKS WB97X-D3/def2-TZVPP;
  RIJONX/def2/J, D3(0), DefGrid3/ExtremeSCF as retained in source inputs.
- Four minimum analytical Hessians; one admitted numerical TS Hessian.
- Arkane frequency scale 1.0. Native Arkane nevertheless uses
  `frequencyScaleFactor/1.014` for ZPE; raw harmonic and native E0 barriers
  are distinguished. The final comparison cross-checks this convention.
- Declared rotational symmetries OH1, CH4_12, H2O2, CH3_6, TS1;
  optical counts1; spin multiplicities2/1/1/2/2.
- No added methane times-four path multiplier, inferred hindered rotor,
  low-frequency floor, deleted vibration or arbitrary frequency scaling.
- Bare TST, Wigner and asymmetric Eckart are computed separately; direct
  native rates are distinguished from Arrhenius fits.

The numerical soft mode near76.15cm-1 remains harmonic. Repeatable numerical
curvature and a connected IRC do not establish harmonic validity. Native
symmetry identifies the TS as C1/optical2 while this baseline declares1;
that is retained as an unresolved counting warning, not silently doubled or
selected according to agreement with measurement.

## Completed DFT-only import and arithmetic check

The [native replay](../../runs/oh-methane-consistent-rate-20261011-XqYhq3/baseline/native-statmech-audit.json)
loads OH as a LinearRotor, the four fragments' complete physical vibration
inventories and all14 real TS vibrations. Native E0 barriers are
11.614508kJ/mol forward and62.174127kJ/mol reverse. Direct rates at340.3K:

| Model | k, m³ mol⁻¹ s⁻¹ | Tunneling factor |
| --- | ---: | ---: |
| Bare TST | 316,588.67 | 1 |
| Wigner diagnostic | 678,297.53 | 2.14252 |
| Asymmetric Eckart | 842,460.00 | 2.66106 |

The Eckart fitted value is842,863.86, approximately0.048% above the direct
value at340.3K. The manual partition-function/TST arithmetic and forward/
reverse detailed-balance identities reproduce native rates. These are
internal-consistency checks, not an independent physical validation.

The DFT Eckart result is approximately100× the frozen measured comparator,
not the historical approximately5200×. Old-to-new changes include functional,
basis, geometry, frequencies and scaling, so that reduction is not attributed
to a single isolated correction. The matched correlated-energy comparison
holds the new geometry/vibration/statistical declarations fixed.

## Matched correlated energies and resources

The five fresh single points use the exact respective frequency-input
coordinates, not old B3LYP/SVP geometries or loosely separated IRC terminals.
Canonical CCSD(T)/def2-QZVPP uses RHF singlets and UHF doublets,
VeryTightSCF, chemical frozen core, NoAutoStart and UseQROs false. Native
normal termination, SCF/CC convergence, final CCSD(T) energy, reference,
frozen core, basis/electron inventory and residual≤1e-6 are checked.
T1 and preceding-reference S² are retained diagnostics, not accuracy proofs.

All five jobs terminated normally. The complete sequential stage took
941.74 seconds; the TS took 815.26 seconds. The matched electronic barriers are
28.015965 kJ/mol forward and 81.582651 kJ/mol reverse. All native CC residuals
are ≤9.09×10⁻⁷. The TS T1=.0212036 is larger than the fragments' .00789–.00885;
its preceding UHF S²=.785456 versus the pure-doublet .75. These diagnostics
do not by themselves establish multireference failure or single-reference
accuracy. No reference-stability calculation, CBS extrapolation, computed SOC,
post-CCSD(T) correction or correlated geometry optimization was performed.

The native composite E0 barriers are 22.051678 kJ/mol forward and 80.820220 kJ/mol
reverse. The raw-harmonic forward electronic+ZPE barrier is 21.968175 kJ/mol;
the difference records native Arkane's ZPE scaling, not an altered energy log.
Source/input/Hessian/state checks and direct partition-function/detailed-balance
identities pass. The final report reopens the source graph, including native
module/interpreter hashes, rather than accepting stored success flags alone.

| Model at 340.3 K | Direct k, m³ mol⁻¹ s⁻¹ | k/measured | Tunneling factor |
| --- | ---: | ---: | ---: |
| Bare TST | 7,915.28 | .93883 | 1 |
| Wigner diagnostic | 16,958.64 | 2.01146 | 2.14252 |
| Asymmetric Eckart | 25,818.37 | 3.06232 | 3.26184 |

Changing only the matched electronic energies raises the native forward E0
barrier by 10.437170 kJ/mol; classical k drops about 40-fold. The asymmetric
Eckart factor itself changes because both bound barrier heights change.
The fitted Eckart rate is 26,312.91 at 340.3 K, 1.915% above the direct value;
reference scoring always uses the direct value.

Across all 14 frozen rows (284.4–340.7 K), the composite prediction/measurement
ranges are .55689–1.00545 (classical), 1.46784–2.15150 (Wigner) and
2.91784–3.27025 (Eckart). Mean absolute natural-log errors are .293316,
.567316 and 1.134689 respectively. These are descriptive development errors,
not calibration. Bare TST's better agreement is not a justification to omit
tunneling: error cancellation among electronic, statistical and dynamical
approximations remains possible. This experiment narrows the discrepancy,
but does not isolate its remaining physical cause.

Jobs run sequentially: two MPI workers, MaxCore2500MB/worker, six-hour
per-job allowance and24-hour aggregate bound. A40GiB starting-space guard
and4GiB runtime reserve protect owned processes; these are allocation
policies, not guaranteed scratch or RSS limits. No user files are deleted.

A final one-line CCSD(T) energy-parser fix occurred after the original builder
hash was bound: an f-string needed an escaped regex quantifier
(`\.{2,}`→`\.{{2,}}`). The original contract is unchanged. An append-only
revision receipt proves that reversing exactly this line reproduces that hash,
and that the corrected builder reproduces the entire retained DFT assessment
and recipe ID exactly. The DFT branch never enters the CC-only parser.
A new guarded completion wrapper admits this exact source revision, not
arbitrary changed code. No chemistry was repeated for this software repair.

- [Source-bound contract](../../runs/oh-methane-consistent-rate-20261011-XqYhq3/contract.json)
- [Energy execution contract](../../runs/oh-methane-consistent-rate-20261011-XqYhq3/energy-contract.json)
- [Completed energy results](../../runs/oh-methane-consistent-rate-20261011-XqYhq3/energy-result.json)
- [Source revision admission](../../runs/oh-methane-consistent-rate-20261011-XqYhq3/builder-revision-admission.json)
- [Composite native statmech replay](../../runs/oh-methane-consistent-rate-20261011-XqYhq3/composite/native-statmech-audit.json)
- [Direct-rate comparison and evidence closeout](../../runs/oh-methane-consistent-rate-20261011-XqYhq3/rate-comparison.json)
- [Revision driver](../../tests/run_oh_methane_rate_reconstruction.py)
- [Guarded completion wrapper](../../tests/complete_oh_methane_rate_audit.py)
- [Native statmech replay helper](../../tests/oh_methane_statmech_inspect.py)
- [Hash-reopening comparison/closeout](../../tests/report_oh_methane_rate_audit.py)

## Reference and interpretation

[Berg et al.2024](https://acp.copernicus.org/articles/24/13715/2024/) measures
methane-dependent OH loss in nominal1000hPa humid synthetic air. The frozen
[Supplement S1](https://acp.copernicus.org/articles/24/13715/2024/acp-24-13715-2024-supplement.pdf)
first point is340.3±0.6K and1.40±0.07×10⁻¹⁴cm³ molecule⁻¹ s⁻¹.
Exact-SI conversion gives8430.997064m³ mol⁻¹ s⁻¹. The published S1 error's
sigma semantics are unresolved. It is not replaced by the paper's fitted
Arrhenius value or treated as a calibrated acceptance interval.

All14 frozen measurement temperatures are evaluated directly, without using
measurement rates to fit the model. They belong to one development study,
not14 independent validations or a prospective holdout. Ideal-gas canonical
direct-abstraction versus finite-pressure humid-air applicability remains
conditional. No full flow-tube network or molecular lifetime is predicted.

The report also labels an OH electronic-state sensitivity, following the
[composite theoretical study](https://www2.sci.u-szeged.hu/czako/papers/PCCP_22_14560_2020.pdf):
140cm-1 spin-orbit splitting, equal-degeneracy ground/excited populations,
ground OH lowered by half the splitting from the spin-free centroid, and
quenched TS spin-orbit coupling assumed. Energy and populations change
together. This is an asymptotic bare-TST control, not newly calculated SOC,
a confidence interval or an Eckart correction. Changing the barrier also
changes Eckart transmission; the reporter does not infer that rate.

At 340.3 K this control multiplies bare TST by .478874 (.643804 from populations
alone), giving 3,790.42 m³ mol⁻¹ s⁻¹. It is not applied to the reported central
rate. Likewise, an independent geometric mirror-superposition diagnostic
on the exact TS permits only the three spectator methane H permutations,
holding the mapped transferring H and OH H fixed. Its best proper-rotation
RMSD is .0747174Å (max .139466Å); coordinate rounding changes RMSD only
2.7×10⁻⁷Å. Thus automatic C1 is not merely a printed-coordinate artifact.
This supports geometric chirality, but does not prove independent harmonic
mirror basins for a 76 cm⁻¹ angular mode. A conditional optical-count 2 sensitivity
would double k; it is neither applied nor selected to improve agreement.

## Remaining gate

Software verification: the full application suite ran 1,438 tests in
172.804 seconds, with no failures/errors and 24 skips under `orca-auto-py312`:

```sh
PYTHONPATH=. python -B -m unittest discover -s tests
```

The initial full run exposed a test-only child import-path defect in the
synthetic P1 interruption check: `python -c` did not inherit discovery's
`tests/` entry in the parent `sys.path`. The worker now explicitly locates
its test driver, and startup failures retain bounded diagnostics. Existing
SIGINT/process-lease/descendant cleanup assertions pass; no production
timeout or chemistry logic changed. Recipe/revision/comparison/roadmap-focused
tests also pass. A fresh independent read-only closeout exactly reproduces
the retained report and composite source/geometry/state/mode assessment.
These software checks do not certify physical accuracy.

The authorized five-energy experiment and native reconstruction are complete.
The next physical gate is a jointly consistent electronic-state/TS-basin
counting assessment, followed by justified soft-mode and tunneling/recrossing
treatment where diagnostics require it. These must not be assembled from
posthoc independent multipliers chosen for agreement. Finite-basis/reference/
higher-level geometry errors and coverage remain distinct uncertainties.
After the recipe is scientifically admitted, test a second independent reaction
before production promotion. This does not close the general working
stability gate, replace the broader all-environment ladder roadmap, change
production defaults or justify an automatic numerical-Hessian fallback.
