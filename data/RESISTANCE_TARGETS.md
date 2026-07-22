# Resistance target data

`resistance_targets_v1.csv` is the versioned long-form companion to the two legacy wide matrices:

- `resistance_prevalence_values.csv`
- `resistance_average_resistant_values.csv`

Version 1 preserves the cleaned wide-matrix values while making calibration inclusion explicit.
Tests require the long-form file to reproduce every value and missing cell exactly.

The prevalence values are classified as **evidence-informed calibration benchmarks** until cell-level
provenance can be recovered. Their linked records in `resistance_target_sources_v1.csv` preserve the
existing bacterium-level notes without upgrading those notes into cell-level evidence.

Conditional mean `any_r` values are **expert-informed model-scale resistance-severity
placeholders**. They constrain a unitless internal model quantity among resistant-positive active
infection person-days and are not direct MIC or breakpoint-surveillance estimates.

Two subsets have more specific provenance. Forty-eight legacy reserve-drug cells that had copied
their prevalence values now use coarse best-guess placeholders: `0.60` for cefiderocol and `0.70`
for ceftolozane/tazobactam. Five `0.60` values paired with zero prevalence benchmarks are retained
as rare-positive structural priors: they specify conditional severity if a positive phenotype is
simulated, rather than asserting that resistant infections are observed. Both subsets have their
own source and rationale identifiers in the long-form files.

`provenance_class` makes the evidence status machine-readable. The allowed classes are direct
empirical estimates with recovered cell-level sources, evidence-informed benchmarks with
unrecovered cell provenance, expert-informed placeholders, structural priors, and unassigned
cells. Version 1 contains no rows in the direct-empirical class. `source_id` and `rationale`
preserve the identity of the source record and design rationale for every numeric cell; the source
table repeats the provenance class so mismatches can be rejected.

Every bacterium-drug-component cell has a row, including cells represented by `.` in the legacy
matrices. `include_in_score` records static v1 eligibility after the target, organism, rifampicin,
baseline-potency, and phenotype-representability exclusions. Potency and resistance-mechanism
checks are materialized from
`model_potency_matrix.csv`, a deterministic projection of the typed Rust matrix that is checked
against Rust in the test suite, and `model_resistance_reachability_matrix.csv`, a checked projection
of whether any applicable mechanism has a positive phenotypic effect for each pair and the maximum
`any_r` attainable if every such mechanism is present. A retained benchmark is not scored when the
current mechanism architecture cannot represent resistance to it, or when a conditional-severity
benchmark exceeds that structural maximum.
Numeric severity values without a paired prevalence benchmark are retained for provenance but are
also inactive, with status `inactive_unpaired_legacy_benchmark`.
Unavailable simulation denominators can still exclude a row from a particular analysis.
`evidence_weight` remains blank because evidence quality has not yet been assessed;
`score_row_weight` records equal static row weighting only.

`resistance_targets_v1.manifest.json` records SHA-256 hashes and byte sizes for the long-form
targets, source table, schema, legacy wide value matrices, and Rust-derived potency and resistance-
reachability matrices. The production loader verifies this manifest before using the target set.
This binds the score input to the mechanism and potency projections used to determine eligibility.

Regenerate the companion files from the cleaned wide matrices with:

```powershell
cargo run --quiet --bin export_potency_matrix -- data/model_potency_matrix.csv
cargo run --quiet --bin export_resistance_reachability_matrix -- data/model_resistance_reachability_matrix.csv
python amr_simulation_output_analysis/build_resistance_targets_v1.py
```

Any future numerical or semantic target change should create a new target-set version or explicitly
document why v1 was amended. Generated parity tests must pass before the wide matrices can be
retired as production inputs.
