# H + methane: exact serial-H2 continuation

Date: 2026-10-23. **H2 and methyl energies passed; final TS energy stopped at
the frozen memory allowance. The attempt cap is exhausted; nothing is running.**

## Scope and concrete repair

The user authorized one core for the H2 matched energy only, and the cumulative
SP attempt cap **10 → 11**. The original contract/ledger is not rewritten.
Failed attempts 24 (PMIx cleanup) and 26 (four ranks exceeding one electron
pair) remain charged, with their original native inputs, outputs and receipts.
The verified route, counting admission and successful atomic-H/methane energies
(23 and 25) are reused. No NEB, TS optimization, frequency or IRC is repeated.

The new H2 input differs from attempt 26 **only** in `%pal nprocs 4 → 1`.
Method, orbital basis, geometry, charge, multiplicity, frozen core, SCF and
diagnostic controls are identical. It uses a fresh append-only basename:
`runs/h-ch4-partial-neb-recovery-20261020/campaign/energies/p1/serial-retry-001/energy.inp`.
Methyl and TS energies retain four workers. The existing PMIx hash workaround
continues for these jobs. Six hours/job, 1000 MB/rank, disk reserve, cumulative
wall/core caps and physical admission gates are unchanged. Even for serial H2,
the conservative original **4.906 GiB** launch-memory requirement is retained;
actual core-time charging uses one core for that job, without refunds.

`storca/rate_transfer_control.py` validates the exact predecessor amendment,
closed journal prefix, failed H2 ownership/error, rank-only input difference
and unchanged other energy bindings. Its one-core permission applies only to
that exact input. Historical readers reopen both amendments and use each
reservation's actual core count. No arbitrary reduction or broad cap increase
is allowed. A further failure still spends an attempt and stops this bounded
campaign; no unlimited retry authority exists.

`storca/rate_transfer_execution.py` now preflights all planned matched-energy
inputs before the first energy dispatch. The known neutral-singlet H2
canonical-CCSD process-count mismatch refuses before native work. This is a
narrow known backend restriction, not a general prediction of optimal MPI
scaling, and does not change the scientific model. Historical four-rank inputs
remain reproducible/readable but are not launchable through this handoff.

## Reproducibility and regression evidence

Before edits, all 32 previously bound source files were archived in
`runs/h-ch4-h2-serial-source-capsule-20261023`. Source-capsule SHA256:
`9f3c013b3a1ddaf427c8b342de89d0a7397cb5d026c35deae29b39aaa29eb028`.
`tests/resume_h_ch4_h2_serial.py` validates that archive and all older evidence,
then binds the new 33-source executing revision. Old contracts and archives
remain intact. No general scheduler, physical model or package was introduced.

Authorized revision:
`runs/h-ch4-h2-serial-continuation-20261023/runtime-contract.json`, SHA256
`506038b2fa94a0720bd6b55517170086fa922167ee784b64021e2d1f58d7952f`.
Amendment SHA256:
`bcefcc7981147a3d9c93a8cb6ea42d26817d0941e797ca5bc70004c65625bd25`.
Original ledger:
`runs/h-ch4-authorized-20261017/campaign/native-control`.

**451 combined regression tests passed** (71.591 s). The final three-test
parallelism subset also passed, including a newly added integrated refusal
before the first energy dispatch. Tests cover actual core-time accounting,
unchanged failed charges, exact input/rank authority, parent/source/model
mutation refusal, the conservative memory guard, original execution behavior,
the retained real campaign and downstream rate/Arkane gates.
`git diff --check` passed. Tests are not scientific rate validation.

## Launch and current checkpoint

Before launch: **9 Opt/path, 7 Freq, 2 IRC, 8 SP**, 26 closed attempts,
zero pending. No native worker was running. At preparation, available memory
was 5.205 GiB against the unchanged 4.906 GiB guard; disk free was 42.03 GiB.
These are snapshots, not guarantees about peak memory. Kitty was opened with
sleep prevention and console capture for:

```sh
python -u -B tests/resume_h_ch4_h2_serial.py run runs/h-ch4-h2-serial-continuation-20261023 --rmg-env rmg4_env
```

Use the existing `orca-auto-py312` environment. Do not start a second controller;
check the owned lease and ledger before resuming. Three energies remain (H2,
methyl, TS), then the existing isolated composite/Arkane audit. Native completion
alone does not establish a defensible rate, lifetime, calibrated uncertainty,
experimental agreement or permission to reinject a rate into a network.

## Actual outcome and bounded stop

The three authorized remaining energy attempts ran:

| Attempt | Species | Cores | Ledger seconds | Outcome |
| --- | --- | --- | --- | --- |
| 27 | H2 | 1 | 1.764 | Normal termination; strict energy reader passed |
| 28 | CH3 | 4 | 24.670 | Normal termination; strict energy reader passed |
| 29 | TS | 4 | 18.709 | Native integral-transformation memory refusal; charged failed |

H2: RHF CCSD(T)/def2-QZVPP, -1.173868289341 Eh, two correlated electrons,
residual 2.17e-7, T1 0.00574261. Output SHA256:
`48a4aed4af46f2a9c3224f9ad9875c50e52979eccd3a062935a44bae65a4f018`.
CH3: UHF CCSD(T)/def2-QZVPP, -39.772675325512 Eh, seven correlated electrons,
residual 6.37e-7, T1 0.008525142. Output SHA256:
`08802b668184edee216f965fcc8d40c3bba1500f0d7005ed7bceccedf8af0012`.
These are source-bound numerical results, not error bars or a demonstration of
single-reference physical accuracy. The existing strict reader admitted both.

The TS log explicitly reports on all four ranks:

```text
ERROR - OUT OF MEMORY !!!
MINIMUM REQUIRED: 1359.7 MB
MAXCORE:          1000.0 MB
```

The failure occurs in the SHARK exchange integral half-transformation, before
a completed CCSD(T) result. This is an explicit ORCA per-rank allocation-limit
refusal, **not evidence that macOS killed the process for exhausting physical
RAM**, and not a TS geometry, Hessian, reaction-mode or IRC failure. The
top-level lease again reports returncode 0, but the native completion check
correctly refuses the abnormal output. No failed energy is admitted.
TS output SHA256:
`37799a3e2089120d592a9ef88b35d732fa13678a70bcdf47c01f188d1887a5b5`;
terminal process-lease SHA256:
`62b3ae108aa6204112ba44ee1fd897ca58655eebce716f8b9b5be435afe843ae`.

Final original-ledger totals: **29 attempts**, with 9 Opt/path, 7 Freq, 2 IRC,
**11 SP**, zero pending, 9039.771 recorded wall seconds and 36153.793 core
seconds. Actual core-time charging reflects the one-core H2 job. This increment
spent three SP attempts (two passed, one failed), reused atomic H and methane,
and launched zero geometry/path/frequency/IRC work. All historical charges and
native files remain. The stopped controller/worker state was checked; no native
process remained. Console and native leases are retained.

Only the matched TS electronic energy remains before the existing composite
and isolated Arkane audit can proceed. Four of five matched energies now exist.
There is still no rate, lifetime, calibrated uncertainty or reinjection claim.
Do not run the above command again: its SP budget is spent and its frozen TS
memory input is inadequate.

The next decision is a **separately authorized TS-only resource continuation**:
proposed 1500 MB/rank on the same four cores, one additional SP slot (11 → 12),
and unchanged method/basis/geometry, six-hour and cumulative time limits. This
exceeds the reported 1359.7 MB minimum, but is not a guarantee every later
allocation will fit. With the same 1 GiB reserve, the guard would require
6.859 GiB available; the post-run snapshot was 7.344 GiB. Re-measure at launch:
MaxCore is not a hard whole-process cap and the margin is limited. No memory
guard may be waived. Preserve failed attempt 29, use a fresh TS-only basename,
reuse all four successful energies, and stop on any new execution/physical
failure. This further memory/budget change has not been authorized or implemented.
