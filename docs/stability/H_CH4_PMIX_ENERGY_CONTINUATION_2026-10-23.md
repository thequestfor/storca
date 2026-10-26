# H + methane: scoped PMIx energy continuation

Date: 2026-10-23. **Methane retry passed; H2 process-count failure stopped the
continuation. No rate is admitted; no native work remains running/pending.**

The reaction route and reached minima already passed. This continuation does
not repeat discovery, NEB, TS optimization, Hessians, IRC or endpoint refinement.
It reuses the completed matched atomic-H energy and attempts only methane,
H2, CH3 and TS matched energies, followed by the existing isolated rate audit.

## Concrete failure and narrowly scoped workaround

Attempt 24 remains a failed, charged methane calculation. Its macOS crash
report identifies `pmix_ds12_lock_finalize` / `pthread_rwlock_destroy` during
MPI finalization. Printed CCSD(T) energy without normal termination is not
accepted. This is not an established memory-kill or chemical-model failure.

The continuation selects the already installed PMIx `hash` data-store
component using `PMIX_MCA_gds=hash`, explicitly forwarded to MPI workers. The
same macOS cleanup stack and this workaround are documented in
[Open MPI issue 7516](https://github.com/open-mpi/ompi/issues/7516).
This is a scoped workaround to test, not an upstream fix or proof of loaded
worker libraries. No MPI installation, ORCA version or production default changes.

The user authorized **one additional SP slot: 9 → 10**. Original contract,
all 24 historical charges, category caps for other jobs, cumulative wall/core
limits, four cores, 1000 MB/rank, six hours/job, 1 GiB memory reserve and 4 GiB
disk reserve are unchanged. If another attempt fails, stop within this budget;
there is no authority for an unlimited retry sequence.

## Implementation and reproducibility

`storca/rate_transfer_control.py` accepts an explicit, immutable amendment for
these five frozen matched-energy inputs only (including cached H). Reservations
record that amendment; historical readers reopen its original closed journal
prefix and authority. Failed charges cannot be refunded. Other chemistry and
unlisted inputs refuse. The runtime environment is changed per dispatch only.

`storca/rate_transfer_execution.py` can reopen the frozen verified path and
raw handoff with exact native ownership/source checks, then reuse its existing
counting, energy, composite, electronic-state and Arkane admission machinery.
The failed methane basename is preserved; retry input in `pmix-hash-retry-001`
is byte-identical. No normally terminated energy is inferred from a partial log.

The campaign-specific wrapper is `tests/resume_h_ch4_pmix_runtime.py`.
Before source edits, the previous 31 bound implementation files were archived
in `runs/h-ch4-energy-runtime-source-capsule-20261023`; manifest SHA256:
`b0a821b1c586ced1a8c55871860f60274adc612f09f335c63a5dd3ce43e00b64`.
Old source/model/evidence bindings are validated alongside the new 32-source
implementation. Old contracts remain unchanged.

Authorized revision:
`runs/h-ch4-pmix-hash-continuation-20261023-r2/runtime-contract.json`, SHA256
`3fb441629afd2121587c9a4799f3a6f820628ec1af477e06cede97d0203d7468`.
Execution amendment SHA256:
`deadb3e8b1b11c402e678ac5e31e91a8e0eeb06a70b751c557c292abb45fb8c6`.
The original ledger remains
`runs/h-ch4-authorized-20261017/campaign/native-control`.

## Verification and launch

**441 combined regression tests passed** (64.575 s), including bounded
amendment accounting, failed retries, mutation/scope refusal, real retained
route/source reopening without native dispatch, energy-only wiring and existing
path/rate/Arkane gates. The first launch exposed a wrapper output-directory
omission before any native reservation; a captured repeat confirmed the exact
exception, also without a new charge. The wrapper now creates that directory;
**56 focused tests passed**, including the new startup regression. The first
zero-native preparation and its exact old driver bytes remain in the unsuffixed
continuation directory; r2 is separately sealed, not a rewrite of old authority.
`git diff --check` passed.

Before launch: 9 Opt/path, 7 Freq, 2 IRC, 6 SP; no pending native reservation.
Available RAM was 5.734 GiB against the unchanged 4.906 GiB requirement; disk
free was 42.07 GiB. These measurements are launch snapshots, not guarantees of
whole-job peak memory. Kitty was opened with sleep prevention and this command:

```sh
python -u -B tests/resume_h_ch4_pmix_runtime.py run runs/h-ch4-pmix-hash-continuation-20261023-r2 --rmg-env rmg4_env
```

Use the existing `orca-auto-py312` application environment. Never run a second
controller while the first owns the ledger. `inspect` can be used after the
owned attempt closes; the historical lineage deliberately refuses pending work.

R2 launch verified: attempt 25 is the methane retry, input SHA256
`190a821f5f137752a82ba702815aab90a56e932241208d108bb90886ea157b95`,
identical to attempt 24. Its reservation binds the amendment and retains
`timeout_seconds: 21600`, `ncores: 4`. The running process lease identifies ORCA
parent/group 27829; four MDCI workers 27863–27866 were observed under mpirun
27862. These are dated observations, not durable permission to stop a PID.
Ledger at that observation: 9 Opt/path, 7 Freq, 2 IRC, 7 SP (one pending).
The completed atomic energy was reused and no path calculation repeated.
Console output is retained in the r2 directory as `console.log`.

No rate, calibrated uncertainty, physical accuracy, lifetime or network
reinjection permission follows from launching this continuation. A completed
rate audit would still be a conditional model result, not experimental validation.

## Actual outcome and exact next decision

Attempt **25** (methane) terminated normally in 28.865 ledger seconds and passed
the unchanged electronic single-point reader: RHF CCSD(T)/def2-QZVPP,
E = -40.451041475117 Eh, CC residual 4.4e-7, T1 = 0.00788639, two frozen-core
and eight correlated electrons. The workaround succeeded on the previously
failing input; this does not establish a universal MPI repair or physical rate
accuracy. Native output SHA256:
`f8c1a4d8ac97b9af96c284fa04cf16f673cc280e608695e5d344c1464d1053c8`.

Attempt **26** (H2) failed after 3.448 ledger seconds. Its explicit native error:

> Number of processes (4) in parallel calculation exceeds number of pairs (1)

This is an execution/resource-selection defect in the frozen four-rank energy
recipe for a two-electron molecule, not a geometry failure, established memory
kill, or evidence that H2 lacks a valid energy. Its output SHA256:
`1e59f31989e91b59f7b924f8db935140c5b16b31b2777f1923c993d949337f9b`.
ORCA's top-level lease records returncode 0 despite the internal MDCI abort;
the normal-termination check correctly refused it. It remains charged as failed.
No methyl or TS energy was launched.

Both attempts are closed. Final counts: **9 Opt/path, 7 Freq, 2 IRC, 8 SP**;
26 charged attempts, zero pending, 8994.629 recorded wall seconds and
35978.516 core seconds. This increment spent two SP attempts (one passed, one
failed), reused the verified route/counting and atomic matched energy, and
launched zero path/optimization/frequency/IRC work. Attempts 24 and 26 remain
failed and preserved. The stopped native worker/controller state was checked.

Three matched energies remain (H2, methyl, TS), but only two slots remain under
the authorized SP cap of 10. **Do not automatically rerun this command.**
The next decision is authorization for a separately source-bound continuation
that (1) permits one core for this H2 energy only, retaining four for the other
energies, and (2) increases the cumulative SP cap **10 → 11**, without refunds.
Keep the same electronic method, geometries, memory/disk reserves, six-hour/job
and cumulative time limits. No physical threshold or method tuning is needed
to address this particular process-count failure. The next implementation must
preflight the supported rank count and test this exception before dispatch,
reuse attempts 23 and 25, and retain both failed inputs/outputs. That exception
has not been authorized or implemented by this handoff.
