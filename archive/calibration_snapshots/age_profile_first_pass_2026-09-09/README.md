# First-pass age-profile parameter update

Applied on 9 September 2026. This changes existing demographic parameter values
only. The simulation's interaction structure, sampler algorithm and disease
parameters are unchanged.

## What was applied

`proposed_weights.csv` is the applied parameter audit: its `before` and `after`
columns record all 108 demographic weights. Exactly 48 change, eight per region:
initial age bands from `[-36000, -32000)` through `[-8000, -4000)` days.

Each region's total configured weight is preserved exactly, retaining the
regional-weight pass that preceded this update. The other 60 weights are fixed:
the future band `[-40000, -36000)`, restored `[-4000, 0)` cohort and all eight
initially living age bands. Consequently the configured regional shares among
all initial records and among initially living records remain unchanged.

The fitting target is the age distribution **conditional on being below 80**,
using groups 0-5, 6-14, 15-49 and 50-79. It is not a fit of the complete five-group
age distribution. The oldest population is deliberately not independently
calibrated from a baseline that omitted its restored birth cohort.

| Region | Smallest multiplier among the eight edited weights | Largest multiplier |
| --- | ---: | ---: |
| Africa | 0.504 | 1.442 |
| Asia | 0.709 | 1.205 |
| Europe | 0.555 | 1.400 |
| North America | 0.675 | 1.305 |
| South America | 0.640 | 1.251 |
| Oceania | 0.786 | 1.129 |

No fitted multiplier reaches the allowed range limits of 0.25 and 2.0. Decimal
precision is retained for reproducibility and exact preservation of regional
weight sums; it does not express demographic certainty.

## Data and method

- UN source: [World Population Prospects 2024, single-age population, both sexes](https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20%28Standard%29/EXCEL_FILES/2_Population/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx).
- Age targets: 1 July 2023, Estimates worksheet, matching the earlier five-region
  reference year. Full North America now combines Northern America, Central
  America and the Caribbean. The six continents sum exactly to World at each
  single age, with total population 8,091,734,930.
- `age_targets.csv` records the full five-group UN profiles. The first four
  groups are renormalized to sum to one for the fitting target.
- Run 103646 supplies 1,460 daily demographic observations for baseline
  2022-2025. It is treated as a pre-fix reference: the original, unscaled weights
  in the earlier archive and the legacy missing-cohort mask are used to infer
  its response. Current region-adjusted weights are frozen separately in
  `before_weights.csv`; regional multipliers are not applied a second time.
- The existing age-overlap kernel and observed demographic response estimate
  map initial cohorts to the reported age groups. The response coefficients
  combine survival, movement and normalization; they are not independently
  measured cohort-specific survival probabilities.
- The offline optimization minimizes relative changes plus a penalty on abrupt
  differences between neighbouring cohort multipliers. Conditional age-profile
  constraints and each region's free-cohort total are linear constraints.
  This calculation does not add any runtime model rule.

An independent check of the rounded applied values found a maximum conditional
under-80 share residual of 7.16e-11 within this response approximation. This is
a numerical consistency check, not validation of the resulting simulation.

## Verification of the applied update

Source comparison confirms that only the 48 demographic values and explanatory
comments changed in Rust; all other parameters and executable logic are
unchanged. The two demographic sampler regression tests pass, as do the four
Python demographic geometry/fitting tests and `cargo check --bins --offline`.
The wider Rust library suite reports 177 passed and one existing unrelated
parameter-expectation failure in
`supported_vaccine_targets_use_organism_specific_acquisition_parameters`.

## What still needs a new simulation

The actual age profile after changing cohort weights and fixing the sampler
has not yet been measured. Age-dependent survival can change the living
regional population shares despite preserved initial regional weights.
The editable band `[-8000, -4000)` partly contributes to ages 80+, so oldest-age
outcomes can also change even though all pure-80+ overlap columns remain fixed.

Check living population by age and region, including 80+, and then the headline
infection, sepsis, mortality and antibiotic-use outputs. The UN 2023 profile is
a reference for the 2022-2025 window, not an exact daily target for every year.
No new full simulation was run as part of this parameter update.

## Reproduction

From the repository root, with the UN workbook at the stated path:

```powershell
.venv/Scripts/python.exe archive/calibration_snapshots/age_profile_first_pass_2026-09-09/extract_targets.py target/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx
.venv/Scripts/python.exe archive/calibration_snapshots/age_profile_first_pass_2026-09-09/fit_age_profiles.py
```

The second command needs the original simulation CSV identified in the script.
Both scripts write audit data beside themselves; they do not edit `src/config.rs`.
The source workbook is kept out of the archive; its SHA-256, geographic mapping
and all 909 extracted source observations are recorded in the manifest and
`un_source_age_counts.csv`. Existing observed data are also retained in
`observed_demographic_summary.csv`.
