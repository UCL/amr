# Regional calibration review: run 261633

Date: 9 September 2026. Review of the latest unfinished request in
"Continue Rust Python model review": whether infection incidence should be the
next calibration focus after the regional mortality and demographic changes.

## Sampling-key fix: 9 September 2026

Following a separate instruction to fix only the sampling bug, the Rust sampler
now looks up the cohort ending at zero as `age_neg4000_0`. A missing expected
demographic key produces an explicit error instead of silently receiving zero
weight. All 108 configured demographic weights and all other parameter values
remain unchanged; the archived fitted weights have not been applied.

The original review below describes the state before this bug fix. Existing
simulation outputs are unchanged; a new simulation is needed to observe the
corrected sampling distribution.

## Subsequent regional-weight adjustment: 9 September 2026

The user subsequently approved the six regional multipliers from
`DEMOGRAPHIC_PARAMETER_REVIEW_2026-09-09.md`. All 18 demographic weights per
region have now been multiplied by Asia 1.06, Africa 1.05, Europe 0.74,
North America 1.09, South America 1.07 and Oceania 0.22. Within-region initial
cohort ratios are preserved. The archived joint age-profile fit remains
unapplied. See `demographic_region_multiplier_changes_2026-09-09.csv` for all
108 before/after values. The run-261633 findings below predate both changes.

## Subsequent under-80 age-profile fit: 9 September 2026

A parameter-only age-profile pass has now changed eight birth-cohort weights
per region while preserving each region's total configured weight. The full
North American age benchmark has been corrected to include Central America and
the Caribbean. The restored cohort and all initially living cohorts retain
their previous weights; an independent 80+ fit remains deferred until its
response can be measured. See
`calibration_snapshots/age_profile_first_pass_2026-09-09/README.md` for the
applied audit, method and validation limits. The original findings below are
historical observations, not results from these new parameter values.

## Main finding

Infection incidence has deteriorated substantially, but the intended demographic
correction is not present in the current Rust configuration. Establish the
demographic baseline before interpreting another incidence-calibration pass.
This review changes no model parameters or Python configuration.

## Inputs and observation scope

- New snapshot: `output_graphs/calibration_summary_aaa261633.txt`; the
  unprefixed `calibration_summary_261633.txt` is absent from the workspace.
- New source: `amr_simulation_output_analysis_outputs/simulation_summary_261633.csv`,
  schema 6.
- Previous source: `amr_simulation_output_analysis_outputs/simulation_summary_226163.csv`,
  schema 4, with current and `zzz` calibration snapshots.
- Both CSV checks select baseline policy 0 and all 1,460 daily rows from
  timestep 33,580 through 35,039 (2022-2025).
- Annual death counts use window totals / 4 multiplied by
  8.2 billion / mean simulated population.
- Headline incidence and scoped headline deaths can be compared across these
  runs. Historical schema-4 regional/age mortality had broader organism scope,
  and its global resistance observations had different timing. Those quantities
  should not be treated as directly comparable to their schema-6 counterparts.

## Headline comparison

| Metric | Previous 226163 | New 261633 | Existing target |
| --- | ---: | ---: | ---: |
| Annual infection deaths, millions | 6.94 | 13.17 | 7.70 |
| People taking antibiotics on an average day, millions | 100.33 | 123.15 | 100.00 |
| Annual person-level infection acquisition metric, percent | 22.11 | 35.92 | 20.00 |
| Annual sepsis incidence, millions | 74.23 | 165.84 | 70.00 |

Mean simulated population fell from 5,863,090 to 5,726,512, approximately 2.33%.
Raw person-level acquisition counts rose from 5,186,402 to 8,227,139,
approximately 58.63%. Thus the incidence increase predominantly reflects more
recorded events, rather than population scaling alone.

Bacterium-incidence mean absolute error increased from 0.0551 to 0.3283
percentage points. The deterioration varies considerably across organisms;
a uniform adjustment would not address the pattern. Sepsis and mortality also
increased faster than acquisition incidence, so incidence alone does not
explain the entire change in severe outcomes.

## Current regional deaths

Schema-6 regional counters count unique person-level infection deaths using
the headline organism scope. Region means effective location at the start of
the death day, including travel.

| Group | Annual deaths in 261633 | Rounded soft target discussed on 9 September |
| --- | ---: | ---: |
| Africa | 1.512 million | 1.9 million |
| Asia | 10.901 million | 3.9 million |
| Rest of world | 0.753 million | 1.9 million |
| Total | 13.166 million | 7.7 million |

The unrounded raw total is 36,777 deaths, producing 13,165,579.37 annual scaled
deaths. Both age and regional partitions reconcile with it. The soft targets
are discussion-level approximations inferred from GBD super-region bars, not
published continent estimates or newly configured calibration gates.

## Demographic correction remains unapplied

All **108 current demographic weights in `src/config.rs` exactly match** the
saved `baseline_weights.csv`. They differ from the proposed
`fitted_demographic_config.rs` in
`archive/calibration_snapshots/regional_first_pass_2026-09-08/`.

The sampler around `src/config.rs:12787` also retains the previous key-generation
condition. For the initial age interval `(-4000, 0)`, it requests
`demo_<region>_age_neg4000_neg0`. The configured key is
`demo_<region>_age_neg4000_0`. The failed lookup defaults to zero weight, omitting
that future-birth cohort. The archived fit explicitly assumes a corrected
sampler.

The new run's age distribution is still distant from the archived UN-derived
age targets. Examples, using population-weighted mean daily age shares:

| Within-region age share | New 261633 | Archived 2023 age target |
| --- | ---: | ---: |
| Africa, ages 0-5 | 12.04% | 16.99% |
| Africa, ages 50-79 | 20.04% | 11.07% |
| Asia, ages 0-5 | 12.24% | 8.59% |
| Europe, ages 80+ | 2.37% | 5.41% |

This establishes the current source state; the CSV itself does not establish
the exact executable source revision. It cannot serve as validation of the
archived fitted configuration when that configuration is absent from the
checked-out model.

The archived demographic targets also need their geographic convention kept
explicit: their North America age profile uses UN **Northern America**, excluding
Central America and the Caribbean. The six-continent mortality placeholders
discussed subsequently assumed those areas were included in North America.
The archived fit preserves the relative population shares of regions within
"rest of world" rather than fitting all six continental population totals.

## Recommended sequence

1. Establish the intended geographic definitions and complete/verify the
   previously proposed demographic correction, including the missing cohort.
   The archived fit is a first-pass approximation and requires validation in a
   new simulation; it is not an already validated replacement.
2. Assess incidence by organism and region on that demographic baseline,
   distinguishing population composition from the number of acquisition events.
   Global averages alone can conceal regional imbalances.
3. Reassess mortality and sepsis conditional on the resulting infection burden.
   The new run has excessive mortality concentrated in Asia, while mortality in
   the rest-of-world group remains below the discussed soft target.
4. Recheck antibiotic use and resistance after those changes. Keep the
   schema-4 resistance timing limitation separate from changes in model fit.

This is a diagnostic sequence, not a new parameter-setting prescription.
No Rust changes, simulation runs, commits, or pushes were performed for this review.
