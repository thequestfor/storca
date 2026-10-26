# Known-reaction transfer preparation

Order: **H + CH4 → H2 + CH3**, then **OH + H2 → H2O + H**.

This package freezes a proposed protocol, not an executable native campaign or
a scoring-ready experimental benchmark. `request.json` binds the completed
conditional OH/methane candidate and a separate reference-curation dossier.
There are no new electronic-structure, RMG or Arkane calculations.

From the repository with the existing application Python environment:

```sh
python -B tests/prepare_rate_transfer_benchmark.py benchmarks/kinetics/rate_transfer_20261017/request.json runs/rate-transfer-preparation-20261017
```

The command creates/reopens three immutable files:

- `prediction-requests.json`: explicitly mapped elementary reactions, charge,
  spin, proposed temperature targets and conditional method identity. Products
  are intentional for this **known-reaction** protocol; this does not alter
  the autonomous fate protocol or claim discovery success.
- `scoring-dossier.json`: experimental sources and curation gaps, separate
  from prediction inputs, with preregistered descriptive engineering metrics.
- `manifest.json`: source identities, protocol, proposed budgets and readiness
  gaps. Its status is `prepared_not_authorized`, not execution-ready.

Identical preparation is reusable without rewriting the files. Changed source
or prepared bytes refuse. References, fitted adjustments and expected answers
cannot be added to the input allowlist. Conditions, units and all frozen cases
must match during descriptive scoring; missing/crashed/unresolved cases remain
in the denominator. Factor-three mean and factor-ten worst-point targets are
engineering goals, not uncertainty calibration or scientific admission.

## Proposed resources—not authorization

For **each** reaction: at most six optimization/path/recovery jobs (including
NEB), four full numerical Hessians, two IRC jobs and five matched electronic
single points covering both asymptotes and the TS. An isolated H atom needs no
vibrational Hessian. All external native invocations count; failed attempts
consume their caps and valid identical work is reused.

Four cores, six hours/job, 1000 MB/rank, 4 GiB disk reserve, and aggregate
24 native wall-hours / 96 allocated core-hours per reaction are proposals.
Internal image/displacement workers consume their owning job's resources.
These ceilings are not runtime estimates or proof that the proposed memory
is sufficient; a separately authorized executable contract must enforce them.

## What must happen before native execution

1. Audit lawful full sources: Rabinowitz 1991, its 1992 correction and Orkin
   2006. Lock actual measured rows, conditions, units, state preparation and
   uncertainty semantics. Current temperatures are selection targets, not
   verified measurement rows; pressure/bath fields deliberately remain null.
2. Resolve or explicitly retain case-specific diagnostic assumptions. Do not
   copy OH/methane TS symmetry/optical counts into another reaction. For
   OH/H2, declare the OH levels and H2 rotational/nuclear-spin population model.
3. Approve a separate budget and condition/method contract. This preparation
   is not that authorization and is not an angular-campaign budget extension.
4. Wire/test the narrowly scoped known-reaction execution adapter using the
   existing `make_orca_route_verifier` and composite/state rate evaluators.
   There is deliberately no new native dispatcher in this preparation-only
   increment. The existing reference-free OH/methane fate driver is not a
   substitute: injecting expected products into it would change its protocol.

No production rate, library reinjection, calibrated uncertainty, supported lane
or lifetime can be promoted from this package. The candidate is frozen with
its unresolved physical assumptions, not declared universally transferable.
