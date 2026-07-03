# Provisional Unused File Audit

Date: 2026-07-02

Scope: tracked repository files outside `archive/`. The `archive/` directory is treated as already out of the active runtime path and was not reviewed file-by-file.

No files were removed or moved. This is a conservative starting list for later scrutiny.

## Active Paths Treated As In-Use

- Rust build/model path: `Cargo.toml`, `Cargo.lock`, `src/` active modules, `tests/`, and `.github/workflows/fast-guardrails.yml`.
- Python analysis path: `amr_simulation_output_analysis/amr_analysis.py` and package modules reachable through `amr_simulation_output_analysis/__init__.py`, `data_loader.py`, `calibration_summary.py`, plotting modules, empirical `data_loader.py`, empirical `normalizers.py`, and `polars_loader.py`.
- Paper table path: `make_paper_tables.py` and `parse_calibration.py`.
- Data/input path: tracked curated files under `data/`, plus `potency_audit_matrix.csv` if present locally because `calibration_summary.py` reads it when generating resistance benchmark caveats.
- Output inputs discovered dynamically: `output_graphs/calibration_summary_*.txt` and matching `amr_simulation_output_analysis_outputs/simulation_summary_*.csv` are runtime inputs for `make_paper_tables.py`, but those directories are ignored/generated rather than tracked source.

## Strong Candidates

These appear to be accidental, stale, or scratch files. They are not used by Rust, not imported by `amr_analysis.py`, not imported/read by `make_paper_tables.py`, and are not active tracked data inputs.

- `,cargo/config.toml` - likely typo for `.cargo/config.toml`; Cargo and `submit.ps1` look for `.cargo`, and `.cargo` does not exist.
- `.codex_make_paper_tables_20260701_184107.err.log`
- `.codex_make_paper_tables_20260701_184107.out.log`
- `.codex_make_paper_tables_20260701_184219.err.log`
- `.codex_make_paper_tables_20260701_184219.out.log`
- `.codex_make_paper_tables_venv_20260701_184333.err.log`
- `.codex_make_paper_tables_venv_20260701_184333.out.log`
- `.codex_make_paper_tables_venv_wait.err.log`
- `.codex_make_paper_tables_venv_wait.out.log`
- `.codex_startprocess_test_20260701_184157.err.log`
- `.codex_startprocess_test_20260701_184157.out.log`
- `.codex_startprocess_venv_noargs.err.log`
- `.codex_startprocess_venv_noargs.out.log`
- `src/config.rs.bak`
- `_sim_pre34e39f7.rs`
- `_sim_34e39f7.rs`
- `_sim_5e9c7ac.rs`
- `codex_prompts/drug_share_update.md`

## One-Off Maintenance Or Patch Scripts

These are not on the Rust build path and are not used by `amr_analysis.py` or `make_paper_tables.py`. Keep only if they are still useful as reproducible maintenance tools; otherwise they are good cleanup candidates.

- `_check_target_potency.py`
- `_fix_intrinsic_targets.py`
- `_fix_intrinsic_targets_pass2b.py`
- `_patch_potency_pass2.py`
- `export_potency_audit.py`
- `filter_calibration_sweep.py`
- `make_html.py`
- `_md_to_html.py`
- `txt_to_html_summary.py`
- `update_parameter_appendix.py`

## Standalone Analysis Utilities Not Used By `amr_analysis.py`

These may be useful manual tools, but they are not reached by the main analysis entrypoint.

- `amr_simulation_output_analysis/multi_run_activity_r_plot.py`
- `amr_simulation_output_analysis/empirical/acquire_empirical_data.py`
- `amr_simulation_output_analysis/empirical/enhanced_empirical_loader.py`

## Historical Calibration Snapshots

These are not used by Rust, `amr_analysis.py`, or the usual `make_paper_tables.py output_graphs/calibration_summary_*.txt` path. They look historical and could either move to `archive/` or be removed after review.

- `calibration_configs/calibration_summary_002921.txt`
- `calibration_configs/calibration_summary_476463.txt`
- `calibration_configs/config_002921_calibration_1.txt`
- `calibration_configs/config_002921_calibration_1_full.rs`
- `calibration_configs/config_476463_calibration_2_full.rs`

## Generated Or Review Artifacts

These are not runtime inputs, but may be intentionally tracked for paper/review communication. They should be reviewed as a policy decision rather than treated as obviously stale.

- `paper_tables/` tracked HTML/SVG outputs and `paper_tables/index.html`
- `MODEL_DESCRIPTION.html` - generated from `MODEL_DESCRIPTION.md`
- `appendix_b_generated.md`
- `model_overview_slides.html`
- `CALIBRATION_DECISION_TABLE.md`
- `potency_matrix.csv` - data-like export; no current references found from Rust, `amr_analysis.py`, or `make_paper_tables.py`

## Repo/Editor/Agent Support Files

Not runtime, but possibly useful for collaboration or CI. These are only cleanup candidates if the team wants a minimal repository.

- `.github/agents/github_amr_model_planner.md`
- `.vscode/settings.json`

## Explicitly Not Candidates In This Pass

- `archive/` - already considered out of active use by request.
- `data/` tracked allowlist - used by calibration summary / empirical loaders or intentionally curated as input data.
- `compute_global_antibiotic_activity.py` - imported by `amr_simulation_output_analysis/plotting/detail_plots.py` when `global_antibiotic_activity` is enabled.
- `parse_calibration.py` - imported by `make_paper_tables.py`.
- `requirements.txt`, `.python-version`, `.gitattributes`, `.gitignore`, `README.md`, `docs/repo_hygiene.md`, active Rust and test files.

## Local Generated Folders Not Pushed

These are ignored/local generated folders. They are not the main focus of this GitHub cleanup pass, but they can be reviewed separately for local disk cleanup:

- `.venv/`
- `target/`
- `output_graphs/`
- `amr_simulation_output_analysis_outputs/`
- `infection_journeys/`
- `amr_branch_checkpoints/`
- `__pycache__/`
