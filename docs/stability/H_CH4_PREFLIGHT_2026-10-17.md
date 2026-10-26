# H + methane: conditional campaign prepared, not authorized

Date: 2026-10-17. Zero new ORCA, MPI or Arkane jobs. The original OH + methane
campaign, candidate freeze and preparation-only transfer package are preserved.

## What is ready

`runs/h-ch4-preflight-20261017/campaign/execution-contract.json` seals a complete
**unauthorized** H + CH4 execution proposal. Its inputs remain reference-free.
The original incomplete draft is retained under `original-draft/`; the resolved
plan and runtime selection are separate artifacts, not rewritten history.

Counting is H1/CH4_12/H2_2/CH3_6/TS3, optical1 throughout, reaction degeneracy1;
there is no additional site multiplier. This is conditioned on native geometric
symmetry and achirality, not a completed geometry result.
[Primary basis: Fernández-Ramos et al., Eq.6/Table2](https://d-nb.info/1239993528/34).

The new H/CH4-specific geometry check runs after the ordinary verified-path and
reached-minimum checks, **before** five expensive matched energy calculations,
and again at composite construction. It enumerates at most 120 same-isotope
nuclear permutations, fits proper rotations, checks the permutation group and
tests reflected superposability. Fixed numerical tolerances are 0.001 Å RMSD
and 0.002 Å maximum displacement, not physical rate-error bounds. Coordinates
are never repaired, counts are not switched after seeing a failure, and exact
IRC atom mapping remains separate from physical indistinguishability in
partition functions. Incompatible/C1/chiral/rank-deficient structures refuse
this declared model. Older execution-plan versions are also checked before
dispatch for incompatible H/CH4 count declarations.

The frozen conditional model specifies the ground-state doublet surface,
classical high-temperature spin-averaged RRHO, numerical vibrations and
matched electronic energies. It does not claim spin-resolved H2 kinetics,
multireference adequacy, calibrated error, a lifetime or network completeness.
Native statmech must actually use the stated classical rotor controls; missing
or quantum rotor controls are not silently relabelled classical.

The runtime selection reopens 140 retained MPI inventory records and seven
selected ORCA driver/helper bindings from the existing installation. It binds
the curated launcher/library search paths and single-thread-per-worker
environment. Four MPI cores and all prior disk/memory/time safeguards remain.
No packages were installed and no library trace is an acceptance gate. This is
bound selection, **not proof of the libraries loaded by future workers**.

## Budget and authorization

The unchanged proposed ceilings are 12 optimization/path, 10 frequency, 2 IRC
and 9 single-point attempts; four cores, 1000 MB/rank, six hours/job, aggregate
48 electronic-job wall-hours/192 allocated core-hours. Failures consume
attempts. The isolated selected-environment statmech audit has its separately
recorded six-hour ceiling. These caps are not expected runtimes.

No launch occurred and no authorization is inferred from preparation. To run,
the user must separately authorize **H + CH4 only** under this proposal. Then
create a new authorized copy and fresh campaign freeze; do not edit the sealed
unauthorized plan in place or alter an existing calculation ledger. Continue
using the working environment through `--rmg-env`; version-2 runtime selection
supplies the pinned ORCA binary without another user resource flag.

## What is still unavailable

The bounded reference followup acquired no full measured table or correction
body. The 1991 primary article returns HTTP403, its 1992 correction remains
inaccessible through the available tools, and the NASA OH/H2 page provides an
abstract rather than a full measured-row table. The separate
`benchmarks/kinetics/transfer_reference_curation_20261017/followup-access-audit.json`
records this without changing the original frozen curation.

Actual measured rows, conditions, units and uncertainty still need primary
audit before experimental scoring. Fitted/generated tables are not substituted.
This prevents scoring, not a separately authorized isolated-gas diagnostic.
The fixed 1000/1200/1500 K grid remains a declared diagnostic grid, not measured
temperatures. No real H/CH4 accuracy or transfer result exists yet.

## Reproduction, zero calculations

```sh
python -B tests/run_rate_transfer_benchmark.py prepare-h-ch4 runs/rate-transfer-preparation-20261017 runs/oh-methane-periodic-candidate-20261016-ha4k1u9s/contract.json runs/h-ch4-preflight-20261017
python -B tests/run_rate_transfer_benchmark.py inspect runs/h-ch4-preflight-20261017/campaign
```

The combined 662-test campaign/regression suite passes, including 31 added
counting/runtime/preparation checks. Whitespace checks pass. Actual package
reproduction returns the identical contract hash and unchanged file-stat
digest; source bindings reopen successfully without execution.

Identical preparation/reopening preserves artifact bytes and mtimes. The
software tests use ideal/falsifying geometries, synthetic runtime inventories,
the actual preparation/freezing interfaces and a refused-geometry no-SP case;
they do not replace a native benchmark. The original OH/methane frozen candidate
still reopens at 13 optimization attempts / 4 NumFreq and no new native work.
