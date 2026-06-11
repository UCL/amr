# Repository Hygiene Notes

This repository currently contains model source, analysis scripts, empirical calibration inputs, generated paper artifacts, historical archives, and local diagnostic scratch files.

The first hygiene pass is intentionally conservative: it repairs ignore rules and documents file classes without moving or deleting tracked files.

One exception is `panic_log.txt`: it is a runtime crash artifact, so this pass removes it from version control while leaving local copies ignored.

## Tracked Source And Inputs

- Rust model source: `src/`
- Rust tests: `tests/`
- Rust manifests: `Cargo.toml`, `Cargo.lock`
- Python analysis package: `amr_simulation_output_analysis/`
- Curated calibration inputs: selected files under `data/`
- Paper/documentation sources and currently tracked rendered outputs: `MODEL_DESCRIPTION.md`, `MODEL_DESCRIPTION.html`, `paper_tables/`

## Local Or Generated Files

These should not be added in future commits unless a PR explicitly documents why they are source artifacts:

- `_diagnostics/`
- `MODEL_FAILURE_TRACKER.md`
- `STACK_OVERFLOW_CRASH_HUNT_PLAN.md`
- `submit_jump.ps1`
- `panic_log.txt`
- `target/`
- `amr_simulation_output_analysis_outputs/`
- `output_graphs/`
- `infection_journeys/`
- `amr_branch_checkpoints/`
- `simulation_run_log.csv`
- `run_metadata_*.txt`
- root-generated documentation and snapshots such as `MODEL_DESCRIPTION.html`, `appendix_b_generated.md`, `model_overview_slides.html`, and `config_working_changes*.patch`
- large data drops outside the curated `data/` allowlist

## Deferred Cleanup

Future PRs should decide, one class at a time, whether to move or untrack remaining generated material such as:

- tracked paper table HTML/SVG outputs
- historical calibration/config snapshots in `archive/`

Do not combine those moves with model behavior changes.
