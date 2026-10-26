# H + methane: one TS-only memory retry

Date: 2026-10-23. **Implemented/tested; actual launch refused by the higher RAM
guard. Zero new native attempts; the authorized retry remains unspent.**

The user approved one final matched TS-energy retry at 1500 MB/rank on four
cores and a cumulative SP cap **11 → 12**. This is a resource correction to
attempt 29's explicit 1359.7 MB minimum versus 1000 MB limit, not a change of
electronic method, basis, state, geometry or physical acceptance criteria.
It is not a guarantee that every subsequent native allocation will fit.

## Exact scope

- Fresh TS basename: `runs/h-ch4-partial-neb-recovery-20261020/campaign/energies/TS/memory-retry-001/energy.inp`.
- Only `%maxcore 1000 → 1500` changes from failed attempt 29. Four cores,
  UHF canonical CCSD(T)/def2-QZVPP, frozen core, original geometry, SCF and
  diagnostic settings are retained. No method tuning or extra search occurs.
- Atomic H, methane, H2 and methyl energies are reused, with original source
  bindings/ownership. H2 retains its one-core source; it is not recalculated.
- Failed attempts 24, 26 and 29 remain charged and intact. The original ledger
  and base contract are not rewritten; no attempt or time is refunded.
- Six hours/job, cumulative wall/core limits, 4 GiB disk reserve and 1 GiB
  memory reserve remain. The raised guard is **6.859375 GiB available**
  (1 GiB plus four × 1500 MiB), enforced immediately before dispatch. MaxCore
  is not a hard whole-process ceiling. A failed guard spends zero attempts.
- No path, optimizer, Hessian or IRC is launched. No installation, external
  issue, network reinjection, production default or environment support changes.

`storca/rate_transfer_control.py` validates the exact H2-continuation predecessor,
closed journal prefix/spending, original failed TS ownership/error, memory-only
input difference and unchanged other inputs. Only the exact new TS input may
launch under this amendment; historical readers reopen all three amendments.
The memory allowance and launch guard increase together. Native-job settings
for the retained DFT/path work remain at their original values.

`storca/rate_transfer_execution.py` forwards the explicit TS-only memory map to
the existing single-point generator, then reuses the existing strict energy,
composite and isolated Arkane admission path. All physical gates remain; a
normally completed electronic energy is not a verified rate or accuracy claim.

## Source-bound revision

All 33 previously bound implementation files were archived before edits in
`runs/h-ch4-ts-memory-source-capsule-20261023`, manifest SHA256:
`ffd045b22ab9ed0812e4fc2070645333551b29a7964caa5973ef63f003be5f81`.
The campaign-specific wrapper is `tests/resume_h_ch4_ts_memory.py`; it binds
the new 34-source executing revision alongside the exact historical archive.
No old contract, source capsule or paid native evidence is overwritten.

Authorized revision:
`runs/h-ch4-ts-memory-continuation-20261023/runtime-contract.json`, SHA256
`661871d21963e6fac75fee9cac611a384291f13b5c2bbc7534718acf782f601c`.
Amendment SHA256:
`c50efea63754c088e18b82e788e6f85dcdff124f36828d68241e6810e4de557a`.
Original ledger remains
`runs/h-ch4-authorized-20261017/campaign/native-control`.

Starting checkpoint: 29 closed attempts; **9 Opt/path, 7 Freq, 2 IRC, 11 SP**,
zero pending. The verified route/counting and four completed matched energies
remain valid. Only the TS electronic energy and subsequent composite/rate audit
remain. No rate, lifetime, calibrated uncertainty or physical accuracy is admitted.

## Regression and resource check

**460 combined regression tests passed** (77.671 s), including exact TS-only
dispatch authority, failed-charge retention, unchanged base contract, native
core-time accounting, memory-only input generation, source/model mutation
refusal, higher memory guard, old source reopening, original path/rate/Arkane
gates and the existing H2 exception. `git diff --check` passed.

After the regression process exited, available RAM was **5.740 GiB**, below
the required **6.859 GiB**. No test process or native worker remained. This is
an operational launch constraint, not chemical evidence. The kitty controller
was invoked with the original safeguards; native launch/result must be read
from its retained console/ledger, not assumed from opening a terminal. No user
applications were closed, memory guard lowered, or cores silently reduced.

The actual kitty invocation reached `NativeExecutionControl.reserve` and
refused with `Native memory reserve is unavailable`. No ORCA job launched.
The controller exited; no native worker or pending reservation remained.
Ledger is unchanged at **29 attempts / 11 SP**, including all historical failed
charges. The four completed energies were reused, not recomputed. Retained
console: `runs/h-ch4-ts-memory-continuation-20261023/console.log`.

The final available-memory snapshot was **5.976 GiB**: roughly **0.9 GiB**
short of the guard. Free at least **1.2 GiB** for some additional headroom,
then rerun the same authorized
command below after checking no controller is active. Do not create a new
campaign, raise the cap again, reduce cores or change the memory guard. The
source-bound TS retry and single extra slot remain available because the guard
refused before reservation. This is a resource block, not another chemistry or
implementation failure. No rate can be admitted until the TS energy and the
existing composite/native rate audit complete.

## Reproduce/resume

Use the existing `orca-auto-py312` application environment:

```sh
python -u -B tests/resume_h_ch4_ts_memory.py run runs/h-ch4-ts-memory-continuation-20261023 --rmg-env rmg4_env
```

Do not start a second controller. Check live ownership/ledger and the guard
before resuming. This exact revision permits one TS attempt, not an unlimited
sequence of memory increases. If it fails, retain its evidence/charge and stop;
do not reset the ledger or change the native input in place.
