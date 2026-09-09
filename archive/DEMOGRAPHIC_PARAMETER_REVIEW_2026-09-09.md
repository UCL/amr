# Demographic parameter review and proposed directions

Date: 9 September 2026. The initial review was a proposal only and changed no
Rust/Python code or parameter values. At that point, source included the
demographic sampling-key correction in commit `2aa772c`. The calibration run
reviewed here, 261633, predates that correction and the regional adjustment below.

## Applied regional adjustment — 9 September 2026

The user approved applying the six region-only multipliers proposed below:
Africa 1.05, Asia 1.06, Europe 0.74, North America 1.09, South America 1.07 and
Oceania 0.22. Each factor has now been applied to all 18 `demo_<region>_age_*`
values in its region, changing 108 demographic weights in `src/config.rs`.
The [parameter audit](demographic_region_multiplier_changes_2026-09-09.csv)
records each original value, multiplier and updated value.

This preserves the relative initial age-cohort weights within each region.
It does not apply the archived fitted demographic configuration or any
age-specific correction proposed in this review. The sampling-key correction
remains in place. Resulting population shares and age distributions have not
yet been validated in a new simulation; none of the observed results below
describe the adjusted weights. Age-profile fitting remains pending, including
reaggregation of the North American age target to cover the full continent.

The original recommendations and comparisons are retained below as the basis
for this first-pass adjustment and the remaining age-profile work.

## Original recommendation

Refit the demographic cohort weights to population and age targets, rather than
using them to fit infection-death targets. The strongest regional changes are
less sampling weight for Oceania and Europe. Africa needs a younger age profile;
Asia, Europe, South America and Oceania need older profiles relative to run
261633. North America's age target must first use the full continent.

All `demo_<region>_age_*` values are relative weights for initial region/age
cohorts at the 1930 start, including future births encoded as negative ages.
They are not percentages of the population alive in 2025. Survival and the
timing of births intervene between these weights and the reported population.

## Regional population targets

Observed shares below are equal-time means of the daily living-population
shares in baseline run 261633, timesteps 33580-35039 (2022-2025).

The reference shares use the sum of UN WPP2024 mid-year populations for 2023
and 2024 divided by the corresponding sum of world populations. These are a
nearby two-year demographic reference, not an exact 2022-2025 daily average.
North America combines Northern America, Central America and the Caribbean.

| Region | Observed population share | UN reference share | Rough region-only weight multiplier |
| --- | ---: | ---: | ---: |
| Africa | 17.57% | 18.43% | 1.05 |
| Asia | 55.59% | 58.97% | 1.06 |
| Europe | 12.41% | 9.17% | 0.74 |
| North America | 6.88% | 7.52% | 1.09 |
| South America | 4.97% | 5.34% | 1.07 |
| Oceania | 2.58% | 0.56% | 0.22 |

The last column is reference share / observed share, rounded. It describes a
possible first-pass proportional change to each region's cohort weights if
their demographic response were fixed. It is not a validated replacement for
the 108 weights. Correcting the missing cohort and changing the age profile
will also affect survival and regional shares, so a combined fit must account
for those effects. Globally scaling all weights alone changes no probabilities.

Sources:

- [UN Population and Vital Statistics Report 2026, Table 1](https://unstats.un.org/unsd/demographic-social/products/vitstats/sets/Series_A_2026.pdf), using WPP2024.
- [UN M49 geographical groupings](https://unstats.un.org/unsd/methodology/m49/).

## Age-profile corrections

Shares are within each region. The reference age shares below come from the
archived WPP2024 single-age extraction for both sexes on 1 July 2023. Observed
values use the same daily-share averaging as the fitting utility; small
differences from person-day-weighted averages are expected.

| Region and age group | Run 261633 | Archived UN age target | Direction |
| --- | ---: | ---: | --- |
| Africa, 0-5 | 12.04% | 16.99% | More young children |
| Africa, 6-14 | 16.88% | 22.44% | More older children |
| Africa, 50-79 | 20.04% | 11.07% | Fewer older adults |
| Africa, 80+ | 0.89% | 0.46% | Reassess after the cohort fix |
| Asia, 0-5 | 12.24% | 8.59% | Fewer young children |
| Asia, 50-79 | 20.33% | 24.05% | More older adults |
| Asia, 80+ | 0.92% | 1.76% | Reassess after the cohort fix |
| Europe, 0-5 | 10.37% | 5.63% | Fewer young children |
| Europe, 50-79 | 26.69% | 35.18% | More older adults |
| Europe, 80+ | 2.36% | 5.41% | Reassess after the cohort fix |
| South America, 0-5 | 12.85% | 8.17% | Fewer young children |
| South America, 50-79 | 19.81% | 24.02% | More older adults |
| Oceania, 0-5 | 11.48% | 9.01% | Fewer young children |
| Oceania, 50-79 | 23.68% | 25.86% | More older adults |

The recent-birth cohort weights, especially those ending in
`age_neg36000_neg32000` and `age_neg32000_neg28000`, influence children in the
calibration window. Earlier birth cohorts influence older adults. Since each
initial band spans roughly eleven years, one weight can affect several current
age groups; changing a single weight by an age-target ratio is not an exact fit.

Age source and extraction provenance:
`archive/calibration_snapshots/regional_first_pass_2026-09-08/demographic_targets_un_wpp2024.csv`
and its accompanying markdown and manifest.

## Why the previous fitted file should not be applied unchanged

1. Its regional targets only explicitly fit Africa and Asia. It preserves the
   relative mix within the rest of the world, leaving Oceania's population share
   at about 2.16%, far above the full-continent demographic reference.
2. Its North American age profile uses UN Northern America only. In 2023 this
   covers 382.903 million people; including Central America and the Caribbean
   raises the total to 608.770 million. The missing areas comprise about 37.1%
   of the full continent, so the age counts need reaggregation.
3. It was derived from run 226163, before the later mortality changes. Its
   demographic-response approximation is no longer a validated description of
   the current model.
4. The restored initial cohort ending at zero affects the 80+ population. Its
   survival response was not observed in the old, buggy-sampler runs. A broad
   80+ response coefficient also cannot identify survival separately for the
   many much older initial cohorts.

The preferable next parameter proposal is a jointly fitted set using consistent
six-continent targets and a simulation with the corrected sampler. Any numerical
fit produced beforehand should remain explicitly provisional. Validate the
resulting population by age and region, and check earlier years if the same
weights are used for historical calibration or policy projections.
