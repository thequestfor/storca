# H + methane: remaining-work resume repair

Date: 2026-10-23. **Resume repair implemented; full route verified; later native
energy work stopped at MPI finalization. No rate is admitted.**

## Concrete defects and changes

1. An interrupted action retained the original full dispatch and completion
   forecast. With both IRC slots spent and only five SP slots left, it was
   incorrectly excluded even though the IRC extension and reached atomic-H
   point had completed. The new forecast reopens exact owned normal outputs,
   validates their journal receipts and subtracts only reusable stages from
   the forecast. Actual charges, failed attempts and cumulative time never
   change. The original five matched-energy reserve remains mandatory.
2. A resumed branch action was dispatched by the outer resume wrapper rather
   than its original operation. Dispatch now retains the original branch
   action. A validated existing assembly is reread without rerunning IRC or
   rewriting the historical continuation result's transient `resumed` flag.
3. The branch-assembly path did not populate the final path-evidence fields
   produced by ordinary completed IRC verification. These fields are now
   derived from the existing frequency count, mapped-mode result and freshly
   revalidated assembled endpoints. All actual reached-minimum gates remain.

Production changes are confined to the adaptive adapter, forwarding its
explicit archived-planner context through `reaction_path`, and accepting an
independently validated archived diagnostic implementation in the narrow
campaign wrapper. No physical model, numerical threshold, counting convention,
production method default, native resource guard or cumulative cap changed.

## Provenance and unchanged budget

Before modifying any bound source, all 30 prior implementation files and the
prior recovery contract/plan were preserved in
`runs/h-ch4-recovery-source-capsule-20261023`. Capsule SHA256:
`53b0f967a95c84cebe724f63cca0d94ada561872bae5cb7b25dcf8435a4ac877`.

`tests/resume_h_ch4_remaining_work.py` is a campaign-specific source-bound
revision, not a new runner or ledger. It validates that exact old archive,
the original authorization/model/ledger and all prior source bindings, then
seals the exact new live implementation separately. Old planner decision
identities are retained only when every scientific/policy field still matches
and its old implementation bytes are verifiably archived. The outer revision
identifies the actual new executing sources. Arbitrary source drift, edited
archives, changed chemistry/resources, missing/duplicate sources and traversal
paths refuse. Old contracts are not edited or hash-waived.

The current authorized revision is
`runs/h-ch4-remaining-work-20261023-r3/resume-contract.json`, SHA256:
`8d5efcd36f5410b82108162e132ef965a0c8ef5839b7310c0c3bdd22a375bb65`.
Earlier same-day r1/r2 preparations were zero-native intermediate revisions;
they are retained, superseded and not the executable handoff.

The original ledger remains
`runs/h-ch4-authorized-20261017/campaign/native-control`.
Before native restart it had 18 charges, no pending work. The saved-run
regression proves the remaining physical work forecast is two optimizations
and two frequency jobs, no IRC or atomic-H SP, plus the unchanged five matched
energy calculations. Four cores, 1000 MB/rank, 21,600 seconds/job, memory/disk
safeguards and cumulative caps all remain unchanged.

## Verification and running checkpoint

**424 combined regression tests passed** against the final implementation.
Coverage includes paid-stage credits without refunds, failed/unconverged-stage
non-credit, accumulated interruptions, wrong/missing/mutated native evidence,
exact external branch ownership, original-action dispatch, no repeated IRC,
physical finalization fields, archived/current source lineage and no weakened
downstream reserve. The real local H + methane fixture was read with native
dispatch forbidden and the full ledger snapshot remained unchanged.
`git diff --check` passed.

Native restart then succeeded. The controller reused the assembled IRC and
reached H/CH4 records. New attempts 19/20 completed H2 Opt/Freq; attempt 21
completed methyl Opt; attempt 22 is methyl NumFreq. At this checkpoint the
ledger counts are 9 optimization/path, 7 frequency, 2 IRC and 4 single points.
Four real MPI workers were observed. Controller PID at observation: 16886;
methyl frequency ORCA parent: 17308. PIDs are historical observations, not
sufficient identity evidence for stopping or restarting a process.

No complete route/rate is claimed at this checkpoint. Methyl frequency and
reached-identity checks, geometric counting, strict composite input admission,
matched energies and Arkane rate auditing still have to pass. The original
ORCA DIIS crash remains unexplained; this repair establishes no new error
bound, lifetime, environmental support or network-reinjection permission.

## Handoff

### Subsequent result: full path success, separate native stop

Both products passed reached-minimum and identity checks. The canonical route
search status is `verified_route`, the path is `verified_transition_state_path`
with `verified_barriered_path`, and all four actual reached species are
validated. Geometric counting admission also passed. The rate-source contract
was published without weakening any scientific admission gate.

Attempt 23 completed the one-electron atomic-H matched energy. Attempt 24,
methane RHF canonical CCSD(T)/def2-QZVPP, converged its CC iterations and printed
the triples-corrected energy, but a worker then failed during finalization.
The output reports rank zero killed with signal 9 and ORCA error termination
in MDCI. No final normal-termination marker exists, so that printed energy is
not a valid completed rate-basis calculation.

The matching macOS report is retained at
`runs/h-ch4-remaining-work-20261023-r3/orca_mdci_mpi-2026-10-23-113522.ips`.
It identifies PID 17699, parent mpirun 17698, with
`BUG IN CLIENT OF LIBPLATFORM: Trying to recursively lock an os_unfair_lock`.
The faulting stack includes `pthread_rwlock_destroy`,
`pmix_ds12_lock_finalize`, `pmix_common_dstor_finalize`, `PMIx_Finalize`,
`ompi_mpi_finalize` and ORCA's `PAL::Finalize`/`TPAL` destructor. This is concrete
MPI/PMIx cleanup-failure evidence, not proof of an ORCA electronic-method bug
or an OOM diagnosis. A narrow kernel-log query provided no matching memory-kill
entry. No runtime workaround or retry has been implemented or tested here.

All 24 attempts are closed, with no pending work and no surviving controller
or energy worker. Counts are 9 optimization/path, 7 frequency, 2 IRC and 6 SP.
Attempt 24 remains charged as failed. Only three SP slots remain; methane,
H2, methyl and TS matched energies still require four completed attempts.
Therefore a retry cannot finish within the existing SP cap. Stop rather than
increase caps, erase the failure, reuse an incompletely terminated result or
alter the runtime silently. Physical accuracy/uncertainty, rate, lifetime and
network reinjection remain unadmitted.

### Next decision

The user must authorize a separately documented native-runtime/budget
continuation before another energy attempt. At least one additional SP slot
is needed; that is not a guarantee a retry will succeed. Diagnose the precise
PMIx finalization failure before choosing a narrowly scoped runtime amendment.
Preserve the verified path, counting result, successful atomic energy and all
failed native work. Do not repeat NEB, TS, IRC or reached-minimum calculations.

Historical command for this exact source revision (not permission to retry a
failed native basename or exceed the original cap):

```sh
python -B tests/resume_h_ch4_remaining_work.py run runs/h-ch4-remaining-work-20261023-r3 --rmg-env rmg4_env
```

Use the existing `orca-auto-py312` application environment. The standalone
revision's `inspect` command validates old candidate evidence and therefore
deliberately refuses while the ledger has pending native work. Do not use the
older strict-live-source launchers after these explicitly archived source
repairs, reset counts, repeat the diagnostic arms or relax a physical gate.
