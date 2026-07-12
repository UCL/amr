# Hot-Path Optimization Validation

## Purpose

This branch reapplies the independent runtime optimizations from the superseded
`perf/hotpath-cleanup-v1` branch to corrected `UCL/amr` main. It deliberately
does not restore the removed syndrome-driven historical empiric-treatment
behavior.

The optimization scope is limited to implementation efficiency:

- indexed treatment and historical score lookup;
- cached immutable rule parameters and execution context;
- removal of hot-loop scratch allocations;
- allocation-free weighted sampling;
- consolidated population cleanup passes; and
- portable release-profile tuning.

No model configuration, rule order, probability, seed derivation, output
schema, or calibration behavior is intentionally changed.

## Source Authority

- Corrected baseline: `6d9dd50df8c851374c3f0b5023454389f446d4ae`
- Optimized runtime: `60699a9608a40b4db5cfeedac5d7cbbe206339ec`
- Common upstream history: `UCL/amr`

The baseline includes the upstream removal of the incorrect
syndrome-driven historical empiric-treatment behavior. The optimized source
contains neither `empiric_explicit_era_initiation_multiplier` nor
`explicit_era_initiation_multiplier_at_year`.

## Quality Gates

The exact source archives were built and tested independently on
`compute72a` using Rust 1.97.0:

- `cargo fmt --all -- --check`
- `cargo clippy --locked --all-targets -- -D warnings`
- `cargo test --locked`
- fixed-seed short-horizon output equivalence
- fixed-seed full-horizon output equivalence

The runtime-source diff is net negative: 1,313 lines added and 1,316 deleted
across 10 files.

## Full-Horizon Result

Both revisions ran 5,000 people for 35,040 time steps in Full calibration
mode with seed `4242424242`, 72 Rayon threads, CPUs 0-71, memory nodes 0-2,
strict configuration validation, and `RUSTFLAGS=-C target-cpu=native`.

Both produced an exact 112,515,100-byte CSV with SHA-256:

`6bb2f1d7f958824a3a41d7993c1c08fef90ad05d313eef68773814c9fa42b231`

| Revision | Elapsed | Maximum RSS |
| --- | ---: | ---: |
| Corrected main | 1,118.056216 s | 566,324 KiB |
| Optimized | 975.646737 s | 571,448 KiB |

The optimized run was 12.737238% faster with a 0.904782% increase in measured
maximum RSS. Exact output identity is the acceptance condition; timing is
supporting evidence and may vary across hosts.

The complete machine-readable parameters and measurements are in
`hotpath-v2-reference.toml`. The obsolete v1 benchmark is not valid evidence
for corrected main and is superseded by this record.

## Workload-Adapter Compatibility

The optimized runtime source was also submitted to the AMR workload launcher
from job-platform revision `791b63429fbd451e7811192ec6312dd362fa4cb1` on
`compute72a`. The launcher materialized the staged source, built it with its
isolated Rust toolchain, and completed a strict fixed-seed 100-person,
100-step Partial run.

The model-created `simulation_summary_014471.csv` and adapter-stable
`summary.csv` were both 8,806,730 bytes with SHA-256:

`35ef30dd40226c51e77ab498f1f17a3e36ccc533f42a36fa01e0b9089a636675`

The adapter also produced its run config, completed run manifest, report,
progress records, and a 103-member `source_used.tar.gz` containing the
materialized `Cargo.toml` and `src/main.rs` without build output. This is an
adapter build/run/artifact compatibility proof, not an installed VM resource
enforcement proof.
