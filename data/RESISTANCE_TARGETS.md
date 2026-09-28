# Resistance target data

The active reference set is now **`resistance_targets_v2`**. Its reviewed corrections,
row-level before/after table, source access record, and reproduction instructions are in
[RESISTANCE_TARGETS_V2_REVIEW.md](RESISTANCE_TARGETS_V2_REVIEW.md). The default loader and
transactional refresh command use v2. The unchanged legacy wide matrices no longer supply
active references or fallback values for v2.

V2 starts with the byte-preserved inputs in `resistance_targets_v1_snapshot/` and explicit
decisions in `resistance_targets_v2.review.json`. Its manifest binds those versioned inputs,
the v2 schema, target/source/audit outputs, and current Rust eligibility projections.
The historical v1 description below records the earlier workflow; refresh now publishes
v2 under `data/.resistance_targets_v2.update.lock`, validates parsed rows and sources before
publication, and writes the v2 manifest last. Existing live v1 files are preserved, and the
snapshot remains independently reproducible even after later projection refreshes.

## Historical version 1 baseline

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
and baseline-potency exclusions. Potency and resistance-mechanism
checks are materialized from
`model_potency_matrix.csv`, a deterministic projection of the typed Rust matrix that is checked
against Rust in the test suite, and `model_resistance_reachability_matrix.csv`, a checked projection
of whether any applicable mechanism has a positive phenotypic effect for each pair and the maximum
`any_r` attainable if every such mechanism is present. A numeric prevalence benchmark remains in
the score when the current mechanism architecture cannot represent resistance to it: the row is
marked `active_target_model_unrepresentable`, and its zero simulated prevalence contributes to the
reported fit as a structural model gap. A conditional-severity benchmark remains inactive when no
positive phenotype can exist or when it exceeds the structural maximum, because mean `any_r`
conditional on `any_r > 0` is then undefined or unattainable.
Numeric severity values without a paired prevalence benchmark are retained for provenance but are
also inactive, with status `inactive_unpaired_legacy_benchmark`.
Unavailable simulation denominators can still exclude a row from a particular analysis.
`evidence_weight` remains blank because evidence quality has not yet been assessed;
`score_row_weight` records equal static row weighting only.

`resistance_targets_v1.manifest.json` records SHA-256 hashes and byte sizes for the long-form
targets, source table, schema, legacy wide value matrices, and Rust-derived potency and resistance-
reachability matrices. The production loader verifies this manifest before using the target set.
This binds the score input to the mechanism and potency projections used to determine eligibility.

The definitive baseline potencies are in `src/config.rs`. After changing them, refresh both
Rust projections, target eligibility, and its hash manifest together. The cutoff applies
to each bacterium-drug pair: a potency below `0.15` excludes both calibration components, while
a potency of `0.15` or greater removes that exclusion. Missing benchmarks and the other
component-specific exclusions still apply. Crossing the cutoff does not change the benchmark
values in the wide matrices.

Run this command from the repository root (with the project's Python environment active):

```powershell
python -m amr_simulation_output_analysis.refresh_resistance_targets
```

The command compiles and exports the Rust projections into a staging directory, copies the
editable inputs without changing their bytes, rebuilds the target/source CSVs and manifest,
and verifies every staged artifact. It then replaces the five generated files while holding
`data/.resistance_targets_v1.update.lock`, publishing the manifest last. Calibration readers
wait during this short publication step. Do not export individual matrices directly into
`data/` during an analysis: their hashes would temporarily disagree with the old manifest.

A staging failure leaves the live files untouched. A publication failure restores the old
files before releasing the marker. If restoration also fails, the error identifies retained
recovery files and a marker that requires investigation before it can be removed. A second
writer cannot overwrite an existing marker. The lower-level target builder remains available
for development and isolated output directories; use the refresh command for the live target set.

Then rerun the affected calibration analyses and check the Rust matrix/target invariants and
Python target-schema tests. Tests that record the previous configuration's eligibility counts
must be reconciled with the regenerated set; the per-pair checks against Rust remain authoritative.

On 26 September 2026, the projections and derived v1 eligibility were refreshed to match the
definitive Rust potency configuration. This excluded 44 additional prevalence targets and 44
additional conditional-severity targets, leaving 1,184 and 1,064 statically included rows,
respectively. Every benchmark value and missing cell was preserved. Version 1 was retained
because this refresh applies the existing eligibility rules to the updated model configuration.

A subsequent E. faecalis potency update on the same date changed 25 values. Clindamycin moved
below the cutoff, excluding its prevalence and severity benchmarks. Ceftaroline moved above
the cutoff, but remains unscored because both benchmarks are unassigned. After this update,
1,183 prevalence and 1,063 conditional-severity rows are statically included; benchmark values
and missing cells remain unchanged.

The subsequent E. faecium update changed 26 potencies: ceftaroline and fidaxomicin moved above
the cutoff, while gentamicin and trim-sulf moved below it. All four lack assigned prevalence
benchmarks, so the included counts remain 1,183 and 1,063. Their exclusion reasons and resistance
reachability were refreshed; all existing benchmark values and missing cells were preserved.

The subsequent S. aureus update changed 14 potencies. Cefixime moved below the cutoff, excluding
its prevalence and severity benchmarks. Ampicillin, amoxicillin, ceftazidime, nitrofurantoin,
and ceftazidime-avibactam moved above the cutoff but remain unscored because their benchmarks
are unassigned. The included counts are now 1,182 prevalence and 1,062 conditional-severity rows.
All benchmark values and missing cells were preserved.

Seven additional S. aureus effective potencies now follow the definitive matrix entries after
removing later fluoroquinolone/tetracycline assignments that overwrote them with `0.50`. These
seven values remain above the cutoff, so this correction does not change scoring eligibility.

The subsequent S. epidermidis update changed 42 raw potencies and 41 effective potencies;
vancomycin changed from a raw `1.05` to `1.00`, with effective potency remaining `1.00`.
Cefixime moved below the cutoff and its two benchmarks are now excluded. Seven drugs moved
above the cutoff but remain unscored because their benchmarks are unassigned. The included
counts are now 1,181 prevalence and 1,061 conditional-severity rows. Benchmark values and
missing cells were preserved.

The subsequent S. pneumoniae update supplied 53 entries and retained the nine earlier entries
omitted from the pasted block. Ten potencies changed, with no crossings of the `0.15` cutoff.
Resistance reachability and calibration eligibility are unchanged, leaving 1,181 prevalence
and 1,061 conditional-severity rows included. The potency projection and manifest were refreshed.

The subsequent S. pyogenes update supplied all 62 entries and changed eight potencies.
Sulfanilamide and retapamulin moved above the cutoff but remain unscored because both
components lack benchmarks. Reachability and exclusion reasons were refreshed, with the
included counts unchanged at 1,181 prevalence and 1,061 conditional-severity rows.
All benchmark values and missing cells were preserved.

On 27 September 2026, the updated S. agalactiae block already present in `src/config.rs` was
synchronized with the projections and documentation. Eight potencies changed relative to the
previous export. Nitrofurantoin and retapamulin moved above the cutoff but remain unscored
because both components lack benchmarks. The included counts remain 1,181 prevalence and
1,061 conditional-severity rows; all benchmark values and missing cells were preserved.

The subsequent H. influenzae update supplied all 62 entries and changed seven potencies.
No values crossed the cutoff, so resistance reachability and calibration eligibility are
unchanged. The potency projection, hash manifest, and documented potency tables were refreshed;
all benchmark values and missing cells were preserved.

The subsequent C. trachomatis update supplied all 62 entries and changed ten potencies. Six
pairs moved above the cutoff, but all six lack assigned benchmarks; rifampicin also retains
its explicit exclusion. Scoring inclusion therefore remains at 1,181 prevalence and 1,061
conditional-severity rows. The staged refresh command published the projections, targets,
source table, and manifest together. All benchmark values and missing cells were preserved.

The subsequent M. genitalium update supplied all 62 entries and changed 17 potencies.
Fusidic acid and metronidazole moved above the cutoff but remain unscored because both
components lack benchmarks. The staged refresh updated reachability and exclusion reasons;
the included counts remain 1,181 prevalence and 1,061 conditional-severity rows. All benchmark
values and missing cells were preserved.

The subsequent N. meningitidis update supplied all 62 entries and changed 16 potencies.
Seven pairs moved above the cutoff. Nalidixic acid is outside the v1 target roster, and the
other six pairs lack prevalence benchmarks. Existing severity benchmarks for ceftaroline,
tigecycline, and flucloxacillin are preserved but remain inactive without paired prevalence
benchmarks. Included counts remain 1,181 prevalence and 1,061 conditional-severity rows.
The staged refresh updated the projections, exclusion reasons, and manifest without changing
any benchmark values or missing cells.

The subsequent N. gonorrhoeae update supplied all 62 entries and changed ten potencies.
Ceftaroline, tigecycline, and flucloxacillin moved above the cutoff, but their prevalence
benchmarks are unassigned. Their existing severity benchmarks were preserved and remain
inactive without paired prevalence benchmarks. Included counts remain 1,181 prevalence and
1,061 conditional-severity rows. The staged refresh updated reachability and exclusion reasons
without changing any benchmark values or missing cells.

The subsequent L. monocytogenes update supplied all 62 entries and changed 31 potencies.
Reachability, cell statuses, and potency-based exclusion reasons were refreshed. The existing
explicit Listeria exclusion remains in force, so the included counts stay at 1,181 prevalence
and 1,061 conditional-severity rows. All benchmark values and missing cells were preserved.

The subsequent C. difficile update supplied all 62 potencies plus the fidaxomicin initiation
multiplier. Forty-six potencies changed; the initiation multiplier remains `1.05`. The staged
refresh updated reachability and exclusion reasons without changing scoring inclusion, which
remains at 1,181 prevalence and 1,061 conditional-severity rows. All benchmark values and
missing cells were preserved.

The subsequent B. fragilis update supplied all 62 entries and changed 22 potencies. Cefixime
moved below the cutoff, excluding its prevalence and severity benchmarks. Nine pairs moved
above the cutoff but remain unscored without prevalence benchmarks. Included counts are now
1,180 prevalence and 1,060 conditional-severity rows. The staged refresh updated the projections,
exclusion reasons, and manifest while preserving every benchmark value and missing cell.

The subsequent M. catarrhalis update supplied all 62 entries and changed 22 potencies.
Six pairs moved above the cutoff but lack prevalence benchmarks. The existing ceftaroline
severity benchmark was preserved and remains inactive without its paired prevalence benchmark.
The staged refresh updated reachability and exclusion reasons; included counts remain 1,180
prevalence and 1,060 conditional-severity rows. All benchmark values and missing cells were preserved.

The subsequent T. pallidum update supplied all 62 entries and changed 11 potencies. Ciprofloxacin,
levofloxacin, moxifloxacin, and ofloxacin moved below the cutoff, excluding their prevalence and
severity benchmarks. Four other drugs moved above the cutoff but remain unscored without
prevalence benchmarks. Included counts are now 1,176 prevalence and 1,056 conditional-severity
rows. The staged refresh updated the projections, exclusion reasons, and manifest; every
benchmark value and missing cell was preserved.

The subsequent B. pertussis update supplied all 62 entries and changed 20 potencies. Thirteen
pairs moved above the cutoff and aztreonam-avibactam moved below it. The affected pairs lack
prevalence benchmarks, so included counts remain 1,176 prevalence and 1,056 conditional-severity
rows. The staged refresh updated reachability and exclusion reasons while preserving all
benchmark values and missing cells.

The subsequent H. pylori update supplied all 62 entries and changed 25 potencies. Nineteen
pairs moved above the cutoff but lack prevalence benchmarks; rifampicin also retains its
explicit exclusion. Included counts remain 1,176 prevalence and 1,056 conditional-severity
rows. The staged refresh updated the projections, exclusion reasons, and manifest without
changing any benchmark values or missing cells.

The subsequent MDR-TB update supplied all 62 entries and changed 13 potencies. The staged
refresh synchronized three reachability rows, ten target exclusion-reason rows, and the hash
manifest. The existing MDR-TB scoring exclusion remains in force: included counts are still
1,176 prevalence and 1,056 conditional-severity rows. All benchmark values, missing cells,
cell statuses, scoring inclusion flags, and weights were preserved. The existing applicability
rules now yield 5,715 cells, four more than before the potency update; resistance parameters
and applicability overrides were not changed.

The subsequent M. pneumoniae update supplied all 62 entries and changed 16 potencies. The
staged refresh synchronized one reachability row, four target exclusion-reason rows, and
the hash manifest. All benchmark values, missing cells, cell statuses, scoring inclusion
flags, and weights were preserved; included counts remain 1,176 prevalence and 1,056
conditional-severity rows. The existing applicability rules now yield 5,721 cells, six more
than before this potency update. The Markdown and HTML potency tables were synchronized;
Python consumers continue to read the generated tables without source changes.

The subsequent Morganella spp. update supplied all 62 entries and changed 11 potencies.
The staged refresh synchronized one reachability row, two target exclusion-reason rows,
and the hash manifest. All benchmark values, missing cells, cell statuses, scoring inclusion
flags, and weights were preserved; included counts remain 1,176 prevalence and 1,056
conditional-severity rows. The existing applicability rules now yield 5,731 cells, ten more
than before this potency update. The Markdown and HTML potency tables were synchronized;
Python consumers continue to read the generated tables without source changes.

The subsequent Proteus spp. update supplied all 62 entries and changed three potencies.
The staged refresh synchronized the potency projection and its hash manifest. Reachability,
all target rows, and the source table were unchanged; included counts remain 1,176 prevalence
and 1,056 conditional-severity rows. The applicability count remains 5,731. The Markdown and
HTML potency tables were synchronized; Python consumers continue to read the generated
tables without source changes.

Version 1 was amended during model development to count model-unrepresentable numeric prevalence
benchmarks as fit penalties while preserving their explicit structural-gap status. Any future
numerical or semantic target change should create a new target-set version or explicitly document
why v1 was amended. Generated parity tests must pass before the wide matrices can be retired as
production inputs.
