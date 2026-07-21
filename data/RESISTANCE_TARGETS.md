# Resistance target data

`resistance_targets_v1.csv` is the versioned long-form companion to the two legacy wide matrices:

- `resistance_prevalence_values.csv`
- `resistance_average_resistant_values.csv`

Version 1 preserves the cleaned wide-matrix values while making calibration inclusion explicit.
Tests require the long-form file to reproduce every value and missing cell exactly.

The prevalence values are classified as **source-informed/expert anchors** until cell-level
provenance can be recovered. Their linked records in `resistance_target_sources_v1.csv` preserve the
existing bacterium-level notes without upgrading those notes into cell-level evidence.

Conditional mean `any_r` values are **expert-informed model-scale resistance-severity
placeholders**. They constrain a unitless internal model quantity among resistant-positive active
infection person-days and are not direct MIC or breakpoint-surveillance estimates.

Every bacterium-drug-component cell has a row, including cells represented by `.` in the legacy
matrices. `include_in_score` records static v1 eligibility after the target, organism, rifampicin,
and baseline-potency exclusions. Potency exclusions are materialized from
`model_potency_matrix.csv`, a deterministic projection of the typed Rust matrix that is checked
against Rust in the test suite. Unavailable simulation denominators can still exclude a row from a
particular analysis. `evidence_weight` remains blank because evidence quality has not yet been
assessed; `score_row_weight` records equal static row weighting only.

Regenerate the companion files from the cleaned wide matrices with:

```powershell
cargo run --quiet --bin export_potency_matrix -- data/model_potency_matrix.csv
python amr_simulation_output_analysis/build_resistance_targets_v1.py
```

Any future numerical or semantic target change should create a new target-set version or explicitly
document why v1 was amended. Generated parity tests must pass before the wide matrices can be
retired as production inputs.
