# Plausibility benchmarks

This benchmark measures condition-specific persistence assessments, not whether
an isolated geometry optimizer can find a stationary point.

Every accepted entry must declare a charge, multiplicity, conditions, expected
category, evidence note, and evidence source. The initial manifest is a draft
candidate list only. It must not tune thresholds or support scientific claims
until each entry is reviewed and marked `accepted`.

The important safety metric is `false_reassurance_rate`: assigning
`ordinary_condition_persistent` to a reference known to be reactive, transient,
special-condition-only, or unsupported under the same conditions.
