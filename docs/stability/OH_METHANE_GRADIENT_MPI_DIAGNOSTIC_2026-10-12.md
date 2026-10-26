# OH/methane fixed-gradient MPI diagnostic

Status: **both native jobs completed; numerical comparison passed**.
The original angular pilot remains incomplete. No rate, lifetime, angular
partition or saddle claim is admitted by this execution diagnostic.

## Question and scope

The +15° continuation failed during cycle-six SCF-gradient execution after
printing an energy and gradient. Test whether the saved geometry can be
evaluated with unchanged chemistry in serial and with MPI matched to the
installed ORCA package. This is a fresh SCF comparison, not literal replay of
the failed worker's orbital state or proof of its original crash cause.

Preserved predecessor:
`runs/oh-methane-angular-continuation-20261012-ws7fj_vt`.
New artifacts:
`runs/oh-methane-gradient-mpi-diagnostic-20261012-b6yb4_ss`.
The [result](../../runs/oh-methane-gradient-mpi-diagnostic-20261012-b6yb4_ss/result.json)
binds the contract and both case receipts. The source inventory binds the old
evidence, participating implementations and native ORCA binaries with hashes
and modification times. Those bindings still reopen unchanged.

## Controlled calculations

Exactly two new EnGrad jobs, with no optimization or Hessian dispatch:

- Saved full-precision seven-atom coordinates; identical mapped atom order,
  no alignment, permutation or reflection.
- Neutral doublet, 19 electrons; UKS WB97X-D3/def2-TZVPP, RIJONX/def2-J,
  DefGrid3, NoAutoStart, ExtremeSCF and the retained SCF controls.
- One core versus four MPI workers; 1000 MB MaxCore/rank, six hours/job,
  4 GiB disk reserve. No automatic retries or overwrite.

Each admission required an owned terminal process, unique native normal
completion, exact input hash, mapped EnGrad coordinates, finite derivatives,
native output/EnGrad agreement, declared state/grid/method and unchanged old
sources. The raw native artifacts were independently reparsed before closing
the comparison; receipt JSON was not the only evidence.

Native SCF acceptance is reported truthfully. It does not imply every printed
residual met its individual target under the retained convergence controls.
The small spin diagnostic is a continuity check, not proof of electronic-state
uniqueness or single-reference applicability.

## Private runtime, not a global installation change

Open MPI 4.1.1 was built from the official archive, whose SHA256
`d80b9219e80ea1f8bcfe5ad921bd9014285c4948c5965f4156a3831e60776444`
was verified before extraction. Native Apple Clang and isolated build variables
were used; inherited Conda-toolchain configuration failures were preserved.
No chemical jobs were consumed by these build attempts.

The prefix is `mpi-build/install` inside the new run. Installed MPI binaries
and components are arm64 and link to that private prefix and macOS system
libraries, not Homebrew MPI 5 or Conda. The library compatibility/current
versions are 71.0/71.1, matching ORCA's worker dependency. No ORCA installation,
Homebrew installation or shell startup file was modified.

Job-local `OMPI_MCA_mca_base_env_list` explicitly forwarded library selection
and tracing to workers through any protected shell. The build/result record
preserves those values. Independent token-based parsing of interleaved native
stderr identifies 24 private `libmpi.40.dylib` loading witnesses across 24
helper/worker PIDs, including all four gradient workers (11314–11317). There
is no alternate MPI path in those witnesses. The frozen line-anchored reader
sees 23 distinct PIDs because parallel stderr can interleave. Launcher selection
alone was not accepted as proof.

## Result and limits

| Quantity | Serial | Matched four-core |
| --- | --- | --- |
| Native normal completion | Yes | Yes |
| Energy / Eh | −116.259368933374 | −116.259368933375 |
| Printed S² | 0.758945 | 0.758945 |
| Native SCF acceptance | 53 cycles | 53 cycles |

Absolute energy difference is 9.9476e−13 Eh in binary64 (exactly 1e−12 at the
stored decimal precision). Maximum gradient-component difference is
1.00000001e−12 Eh/bohr. These pass the declared engineering tolerances
1e−8 Eh/1e−7 Eh/bohr, not a physical uncertainty model. Serial also reproduces
the incomplete predecessor's printed derivatives within their rounding.

Runtime performance is not validated: ORCA reported 47.345 seconds serial
versus 914.167 seconds matched MPI. The four-core module breakdown assigns
704.480 seconds to startup, 42.321 seconds to SCF and 165.151 seconds to gradient
evaluation. This is not evidence of a chemistry-setting problem or a reliable
parallel slowdown benchmark. Scoped macOS power logs show clamshell sleep at
21:24:21 EDT for 704 seconds, dark wake at 21:36:05, maintenance sleep at
21:36:50 for 160 seconds and dark wake at 21:39:30. These closely account for
the startup/gradient delays. Subtracting the 864 logged sleeping seconds from
the owned wall interval leaves approximately 50.4 seconds, but that subtraction
is not a controlled active-runtime measurement. Sleep-log resolution, dark-wake
scheduling and tracing preclude a speed claim. The separate
[timing evidence](../../runs/oh-methane-gradient-mpi-diagnostic-20261012-b6yb4_ss/sleep-timing-evidence.json)
records the diagnostic attribution. No power settings were changed.

The original macOS wakeup warning did not establish a process kill, and this
one successful matched-runtime test does not prove MPI 5 caused that failure.
The angular ledger remains five saddle attempts/two NumFreq; the two EnGrad
jobs have a separate explicit diagnostic budget. A separately authorized,
append-only +15° optimization/derivative continuation is still needed.

## Software verification

[Diagnostic driver](../../tests/run_oh_methane_gradient_diagnostic.py) and
[18 focused tests](../../tests/test_oh_methane_gradient_diagnostic.py) cover
fixed-frame inputs, mapped output parsing, mutation/redispatch refusal, actual
runtime evidence, failure retention and declared numerical comparison gates.
All 18 pass. This is not a fresh full-suite or general ORCA-runtime certification.
