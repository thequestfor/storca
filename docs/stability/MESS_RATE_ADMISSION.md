# MESS rate admission

MESS execution and scientific rate admission are deliberately separate. A
successful `mess` process produces `completed_unadmitted_rate_table`; it does
not create a verified rate, repair a network, or support a lifetime.

The finalization command is:

```bash
storca stability-finalize-mess-rate \
  prepared-mess/ \
  mess-admission-evidence.json \
  specialist-rate-request.json \
  --rmg-env rmg4_env
```

Only `--rmg-env` is exposed because the owning RMG environment must be
declared. The interpolation model and acceptance tolerance are implementation
policy, not additional chemical conditions. Programmatic callers may select
`Chebyshev` or `PDepArrhenius` and a stricter tolerance.

## Preconditions

`prepared-mess/` must contain an immutable, successful execution made from the
same request:

- `mess-input-manifest.json`;
- `mess-execution-manifest.json`;
- `mess.inp` and the hashed Rotd evidence;
- native `mess.out`, `mess.log`, and `mess.aux` when emitted;
- `mess-rate-table.json`; and
- retained stdout and stderr.

The execution manifest, prepared input, executable hash, and every declared
artifact are revalidated before admission. Existing preparation bundles that
predate retained fragment labels can still be inspected and reused as
execution evidence, but must be regenerated from their original model before
rate admission.

## Required admission evidence

The JSON object has `schema_version: 1` and
`kind: storca_mess_rate_admission_evidence`. It must contain the exact
`request_sha256`, route and reaction identity, endpoint labels and SMILES, and
one requested condition that occurs exactly on the computed grid.

The four physical evidence blocks are mandatory.

### Bimolecular amount basis

Native MESS output labels bimolecular rates as `cm^3/sec`. STORCA does not
infer the missing amount basis. The evidence must explicitly establish:

```json
{
  "bimolecular_unit_interpretation": {
    "status": "verified",
    "native_units": "cm^3/sec",
    "interpreted_units": "cm^3 molecule^-1 s^-1",
    "si_units": "m^3/(mol*s)",
    "conversion_factor_to_si": 6.02214076e17,
    "executable_sha256": "<exact binary SHA-256>",
    "source": "<primary source or executable-specific assessment>",
    "source_kind": "<source category>",
    "applicability": "<why the source applies to this binary>",
    "evidence_artifacts": [
      {"path": "<retained file>", "sha256": "<file SHA-256>"}
    ]
  }
}
```

The conversion factor is accepted only after the evidence asserts a
per-molecule native basis. A per-mole interpretation, missing binary binding,
or unsupported unit spelling is rejected. First-order `1/sec` values convert
directly to `s^-1`, but the bimolecular interpretation remains necessary to
check the reverse channel.

### Collision ceiling

For a bimolecular requested direction, `collision_limit.values` must cover
every temperature in the MESS grid in `m^3/(mol*s)`. Each rate at every
pressure is checked against the corresponding ceiling. The model, source,
applicability, and immutable evidence artifacts must be retained. A collision
ceiling is only a physical upper bound; STORCA never substitutes it for a
capture rate.

### Independent thermodynamic consistency

`thermodynamic_consistency` must provide an independently sourced expected
requested-to-reverse rate ratio for every temperature. Units depend on the
requested direction:

- association divided by dissociation: `m^3/mol`;
- dissociation divided by association: `mol/m^3`.

MESS must have emitted exactly one complete table for the reverse channel.
STORCA converts both directions, compares their ratio with the supplied value
at every temperature and pressure, and requires a declared relative tolerance
no larger than 5%. This checks dimensional and thermodynamic consistency; it
does not prove that the shared thermochemistry is accurate.

### Rate uncertainty

Multiplicative lower and upper rate bounds must satisfy
`0 < lower <= 1 <= upper` and retain source, source kind, and applicability.
The uncertainty record must cover at least:

- `electronic_structure`;
- `conformer_sampling`;
- `statistical_mechanics`;
- `capture_or_transition_state`;
- `energy_transfer`;
- `collision_transport`; and
- `master_equation_numerics`.

No default multiplier is supplied. Listing a component is a claim that the
provided bounds cover it, so that claim must be supported by the named source
and applicability statement.

## Products and claim boundaries

Successful finalization writes:

- `mess-rate-admission.json`, including every check and immutable hash;
- a one-reaction `verified_RMG_library` with `reactions.py`, `dictionary.txt`,
  and `generated-kinetics.json`; and
- `specialist-rate-result.json` for the ordinary retained specialist provider.

RMG fits the complete admitted table, reloads the serialized library, and
reevaluates every original temperature/pressure point. A fit outside the
declared log-rate tolerance is rejected rather than hidden by evaluating only
the requested condition.

The specialist result can continue the original stability run:

```bash
storca stability-resume-verification runs/<run> \
  --specialist-rate-result prepared-mess/specialist-rate-result.json
```

The continuation reuses the RMG environment retained by the original stability
run; it is not a scientific override supplied at resume time.

The existing verification engine then inserts the checksummed library, reruns
the complete RMG model under the identical condition contract, verifies that
the replacement rate was actually applied, repropagates the network, reranks
flux and sensitivity, and propagates the supplied uncertainty.

Admission permits a condition rate and a pressure-dependent model as network
repair inputs. It still prohibits:

- a lifetime before repaired-network propagation;
- use outside the declared temperature, pressure, and bath-gas domain;
- a claim that MESS, Rotd, the electronic-structure method, or the uncertainty
  model has been experimentally validated; and
- promotion of the MESS lane to supported status before preregistered live
  falloff benchmarks and locked holdouts pass.
