Active reference set: **resistance_targets_v2**, reviewed 28 September 2026.

The correction changes 18 prevalence references and the provenance or pairing status of
18 corresponding severity rows. All conditional-severity numbers are preserved. The
[source review and complete before/after table](data/RESISTANCE_TARGETS_V2_REVIEW.md)
explain each organism-drug decision. The
[36-row revision audit](data/resistance_targets_v2.revisions.csv) records old and new
values, statuses, inclusion flags, reasons and sources for both components.

No new positive numerical reference was introduced. T. pallidum tetracycline,
doxycycline and minocycline use structural-prior centres of 0.0, with no invented
uncertainty interval. The minocycline prior is explicitly class-based. E. faecium
tigecycline and linezolid replacements remain unresolved: exact all-region resistant
results and susceptibility criteria could not be recovered from the accessible ATLAS
abstract/highlights. Both old references are now unassigned, as are the 13 individually
audited unsupported C. trachomatis references. Restricted Australian and Turkish
results were not converted into global targets. Source observation years, publication
years, denominators, access levels and limitations are recorded in
[the version-bound review input](data/resistance_targets_v2.review.json).

E. faecium conventional tetracycline and tedizolid remain unchanged. Tedizolid's existing
0.10 prevalence reference remains unresolved for separate review. Ampicillin and
amoxicillin remain missing throughout loading, scoring and comparison outputs.

Static prevalence inclusion is 1,176 to 1,161; severity inclusion is 1,056 to 1,043.
All 15 withdrawn prevalence rows have unpaired severity exclusions. Two corresponding
tigecycline severity rows were already excluded; therefore 13 severity rows newly leave
the score. Explicit zero references receive ordinary finite percentage-point residuals.
Missing/excluded references receive no residual and cannot enter a paired average.

The changed and added implementation files are:

- `amr_simulation_output_analysis/build_resistance_targets_v2.py`
- `amr_simulation_output_analysis/refresh_resistance_targets.py`
- `amr_simulation_output_analysis/calibration_summary.py`
- `amr_simulation_output_analysis/make_paper_tables.py`
- `amr_simulation_output_analysis/parse_calibration.py`
- `amr_simulation_output_analysis/plotting/detail_plots.py`
- `data/calibration_targets.json`, `data/RESISTANCE_TARGETS.md`, and `data/RESISTANCE_TARGETS_V2_REVIEW.md`
- `data/resistance_targets_v2.csv`, `data/resistance_target_sources_v2.csv`, `data/resistance_targets_v2.schema.json`, `data/resistance_targets_v2.manifest.json`, `data/resistance_targets_v2.review.json`, and `data/resistance_targets_v2.revisions.csv`
- `data/resistance_targets_v1_snapshot/`: eight byte-preserved original target/source/schema/manifest and wide/projection inputs
- `tests/test_resistance_targets_v2.py`, `tests/test_resistance_reference_scoring.py`, `tests/test_resistance_paper_comparison.py`, and `tests/test_detail_resistance_references.py`
- `tests/test_resistance_target_snapshot.py`, `tests/test_refresh_resistance_targets.py`, and `tests/test_parse_calibration.py`
- `.gitignore`, `.gitattributes`, this report, and the comparison output directory below

The production loader verifies the selected version's complete manifest. The refresh
command stages generation, validates row/source semantics, publishes under the v2 lock,
and writes the manifest last. Reproduction tests cover both v2 and the independent v1
snapshot. Historical live v1 artifacts and wide matrices are unchanged. New version-bound
artifacts and the live projections have checkout attributes preserving hash-bound bytes.
There is no active fallback to the old wide matrices or embedded summary targets.

The regenerated [comparison bundle](paper_tables_resistance_targets_v2/index.html)
contains Supplementary Table S2 and Figures 2i, 2ii, 2iii, 2A and 2B. Each figure has
HTML, PNG, SVG and a class-comparison CSV. There are 26 artifacts in total, including
the index, row overlay audit, class membership audit, and provenance text/JSON.

The original accepted runs are unchanged: `fff164402`, `fff311750`, `fff591413`,
and `fff928228`, paired with their existing schema-6 simulation CSVs. The original
summary hashes remain unchanged. Original summary metadata and embedded reference
values are retained separately from the v2 comparison version and its target/manifest
hashes. The original summaries did not declare a version hash; that absence is recorded
rather than retrospectively assigning one. No archived `paper_tables/` output is replaced.

The targeted class effects across those four runs are:

| Comparison | Reference mean before (%) | Reference mean after (%) | Paired model mean before (%) | Paired model mean after (%) | Membership effect |
|---|---:|---:|---:|---:|---|
| T. pallidum tetracyclines | 6.333333 | 0 | 16.9625 | 16.9625 | Same three conventional drugs; explicit zero marker |
| E. faecium tetracyclines | 43.75 | 45 | 40.033125 | 41.81 | Tigecycline excluded from both means |
| E. faecium oxazolidinones | 10 | 10 | 29.595 | 29.595 | Tedizolid only; linezolid excluded |
| C. trachomatis tetracyclines, fluoroquinolones, macrolides, lincosamides and chloramphenicol | Previously assigned | Unavailable | Previously paired | Unavailable | Five classes have no usable reference; simulation-only information retained |

Paired means use the same eligible drugs within each run; run weighting and confidence
interval methodology are unchanged. Eight targeted organism-class combinations affect
32 run/class audit rows. The full before/after membership and aggregate values are in
[the class audit](paper_tables_resistance_targets_v2/resistance_class_membership_audit.csv).
Individual drug/run model measurements are unchanged. A changed paired mean reflects
a changed contributing drug subset, not changed model biology.

Validation actually completed:

| Command/check | Result |
|---|---|
| `.venv/Scripts/python.exe -m unittest discover -s tests` | All 306 discovered tests passed in the final full run |
| `cargo test --test target_matrix_invariants --test resistance_reachability_matrix_invariants` | All 3 tests passed |
| `.venv/Scripts/python.exe -m amr_simulation_output_analysis.refresh_resistance_targets` | Successful production staging, validation and publication; projections and v1 remained unchanged |
| V1 and v2 reproducible generation, schema/source consistency and strict manifest verification | Passed within the focused tests |
| Git diff scope check and `git diff --check` | No biological/v1/original-bundle changes or whitespace errors |
| Five-figure and S2 regeneration from the four original runs | Completed; original summary hashes verified |

`cargo check` was not separately run: no compiled reference-loading/scoring code changed,
and the targeted Rust integration tests compiled successfully. No new simulations,
recalibration, commits or pushes were performed. Historical calibration score files were
preserved; this bundle revises their reference comparisons without recomputing archived
overall scores. The inaccessible ATLAS full tables remain the evidence gap, not an
unfinished numerical substitution.

`potency_when_no_r`, resistance mechanisms, simulation rates, prescribing, treatment
effects and original simulation results were not changed. Neither improved agreement
from revised references nor removal of unsupported targets demonstrates improved
simulated biology.
