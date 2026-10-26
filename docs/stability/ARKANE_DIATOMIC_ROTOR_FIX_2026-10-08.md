# Native Arkane diatomic rotor checkpoint — 2026-10-08

Status: **software workaround demonstrated; physical rate accuracy unresolved**.
This checkpoint does not close the working stability gate or expand environment
support. It repairs an observed ingestion defect, not the underlying electronic
structure or rate theory.

## What calculates what

RMG discovers and estimates the mechanism. ORCA supplies electronic energies,
optimized stationary points, frequencies and IRC evidence. Arkane (distributed
with RMG-Py) constructs statistical mechanics and elementary rates. Cantera
propagates the gas mechanism to macroscopic concentrations. Graph enumeration
and xTB provide independent, inexpensive discovery/screening evidence; RDKit
supports molecular graphs and geometry preparation. Multiwfn is not the primary
kinetics or lifetime calculator. STORCA owns evidence admission, selection,
workflow and scientific reporting.

## Bounded implementation

The native ORCA reader's exact-zero principal-inertia test can misclassify an
obliquely oriented diatomic as a nonlinear rotor. The defect is reported in
[RMG-Py #3016](https://github.com/ReactionMechanismGenerator/RMG-Py/issues/3016).

Verified gas rates now use a process-local native-reader adapter under the
existing `--rmg-env` selection. No installed RMG/Arkane files are edited; no
additional user flag or new quantum calculation is required. For exactly two
distinct atoms, the adapter substitutes a native `LinearRotor` with analytic
transverse inertia. Native symmetry and quantum convention are retained, as
are energy, translation, vibrations, spin multiplicity and optical isomers.
Atoms and clearly nonlinear molecules retain their applicable native modes.
Already-correct linear polyatomics are checked and preserved, not repaired;
ambiguous or inconsistent polyatomic classification is refused.

The source/input/Hessian/native-module/adapter hashes, complete stationary-point
inventory and before/after native modes form an append-only contract and
receipt. Admission checks both full-precision analytic inertia and the exact
native serialized output separately. Missing, changed, incomplete or borrowed
route evidence cannot authorize an Arkane rate replacement. A mandatory
`arkane_rate_source_artifacts` record independently binds the rate's own input,
output and provenance to the receipt; matching route labels alone is not enough.
Any supplied conditional forward evidence must agree with these bindings.
Historical
unaudited rates remain diagnostic artifacts and require a fresh rate workspace.
Installed packages and original chemistry artifacts are not rewritten.

The receipt proves mode admission, not an independent reconstruction of every
rate metadata value. The producer still compares native Arkane evaluation with
the loadable RMG library rate. The standalone replacement predicate checks
condition, molecularity, ceiling and evidence, but does not reevaluate arbitrary
numeric rate-field edits. End-to-end rate-value integrity remains a distinct
hardening concern; no general tamper-proof rate claim is made here.

The adapter also owns symmetry scratch directories and resolves the admitted
Hessian directory for the older native ORCA reader. These compatibility
operations do not waive chemistry or physical applicability checks. Portable
geometry reductions use `math.fsum` to avoid a demonstrated Python 3.11/3.12
summation difference while retaining strict evidence comparison.

Implementation: `storca/arkane_rotor.py`, `storca/arkane_launcher.py`, the existing
`src/orca_runner.py` launcher and `storca/arkane_runner.py` producer, and
`storca/generated_kinetics.py` replacement admission. See the rotational-mode
section of `METHODS.md` and the corresponding limitation in `LIMITATIONS.md`.

## Actual saved OH + methane replay

The successful fresh workspace is
[replay 888ajm6u](../../runs/arkane-diatomic-rotor-replay-888ajm6u/replay-result.json).
It reuses the saved endpoint/TS/frequency/IRC evidence from
`runs/oh-methane-benchmark-003/workflow/verification-engine/iteration-001-rmg-53a85fa71324/rate`.

- Selected environment: `rmg4_env`, native RMG/Arkane 4.0.0.
- Route: `rmg-53a85fa71324`; condition: 340.3 K, recorded 1 bar.
- Elapsed replay: **27.604 seconds**; new electronic-structure jobs: **zero**.
- Original artifact inventory unchanged: **true**.
- OH transverse inertia: **0.9126586911957871 amu·Å²**.
- Corrected OH rotational partition function at 340.3 K: **12.805036286251942**.
- Generated central rate: **44,047,847.000151 m³/(mol·s)**, or
  **7.314317076864706e-11 cm³/(molecule·s)**.
- Retained collision ceiling: **190,571,421.80130282 m³/(mol·s)**;
  rate/ceiling: **0.23113563714750998**, passing the existing physical gate.
- Complete native statmech audit: **passed**. Rate uncertainty remains
  **calibration_required**. No full network rerun or physical lifetime is claimed.

The receipt, contract, provenance, native output and generated library are in
that new workspace. Earlier fresh development attempts `dzua9c5t` (native output
serialization mismatch) and `tyevwssc` (cross-interpreter reduction mismatch)
remain retained. They are not successful certified rates and were not overwritten.
Earlier replay `6cj1uyvm` passed native mode validation in 28.594 seconds before
the final independent rate-source binding was added. It remains historical;
its manifest is not rewritten to gain the new admission field.

## Accuracy is not fixed

The frozen, proposed development comparator in
[the Berg 2024 measured-rate slice](../../benchmarks/stability/reference_sources/2026-10-07-berg2024-oh-methane/experimental-rate-slice.json)
is **1.4e-14 cm³/(molecule·s)** at 340.3 K. The descriptive ratio is about
**5,224.5**. The measurement is a nominal finite-pressure loss coefficient;
the generated model is high-pressure-limit gas TST. Their applicability match
is still unresolved, so this ratio is not a calibrated error score.
Neither passing the collision ceiling nor removing the rotor bug establishes
rate accuracy, mechanism completeness or an experimental lifetime.

The next bounded scientific task is to audit the retained barrier/energy and
rate construction: endpoint/TS energy and ZPE conventions, artifact/state
consistency, H-transfer tunneling, symmetry/path degeneracy and applicability.
Use the source-bound discrepancy to distinguish software conversion errors from
physical model errors; do not tune a barrier to this single reference or treat
DFT as ground truth. Escalate calculations only after the audit identifies the
missing evidence. No such new chemistry job is launched by this checkpoint.

## Verification and reproduction

The dedicated rotor suite passes **10 tests in 64.908 seconds**, exercising
actual native factories, statmech loading and CLI execution in both `rmg4_env`
(RMG 4.0.0, Python 3.11) and `rmg_env` (RMG 3.3.0, Python 3.9), including
unchanged installed-module hashes and strict tampering refusals. After the final
independent rate-source binding, **134 targeted integration tests** pass,
including **13 strict admission tests**. The final full application suite runs
**1,417 tests in 135.676 seconds**, with **24 expected skips** and no failures;
`git diff --check` also passes. An intermediate regression exposed a fixture
that precreated `input.py` before the new fresh-workspace guard; input creation
was moved into the mocked creation step, not allowed through the guard.
These are software/native-ingestion checks, not experimental validation.

From the repository root with the existing STORCA application environment active:

```sh
python tests/replay_arkane_diatomic_rotor.py --rmg-env rmg4_env
```

This automatically creates a fresh workspace, reuses completed chemistry and
performs the Arkane-only replay. It is a development replay, not a new stability
CLI option. Application regression command:

```sh
PYTHONPATH=.:tests OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m unittest discover -s tests
```
