Review of section K, `src/config.rs`, 12 September 2026.

**Assessment: the weights form a coherent, reproducible first-pass demographic calibration, but the older-age and forecast components remain provisional.** All 108 current values match the applied 9 September audit exactly. The strongest evidence for revisiting them is the excess 80+ population in recent outputs, followed by the imposed birth-cohort discontinuity around late 2028. No configuration or model code was changed during this review.

The [complete 108-weight audit](demographic_weight_audit_2026-09-12.csv) records every literal, source line, initial-age interval, implied birth-year interval, age at the start of 2025, normalized sampling share, prior value, fitting status and assessment. The [recent-run comparison](demographic_outcome_review_2026-09-12.csv) covers runs 071999, 545650 and 391322. The [matched-date comparison](demographic_matched_date_review_2026-09-12.csv) checks the two latest runs against the age reference at approximately its own date.

**What the values mean.** The sampler draws from one joint region/initial-age distribution, normalized by the total of all 108 weights. It then draws uniformly within the chosen 4,000-day band. A negative initial age schedules a future birth; a nonnegative initial age represents someone already living at the 1930 start. These are allocated-record probabilities, not the proportions of living people in 2025. See [the sampler](../src/config.rs#L12747) and [initialization](../src/simulation/population.rs#L2083).

All 108 weights are unique, positive and finite; each region has the same 18 contiguous bands. The total is 1.3769, which is valid because sampling normalizes it. Multiplying every weight by the same constant changes no sampling probabilities. The repaired `neg4000_0` key is correctly consumed. Conditional on each region, the initial age distribution is also normalized implicitly by the region's total.

| Region | Raw total weight | Share of allocated records | Mean living-population share, run 071999 | UN population reference |
| --- | ---: | ---: | ---: | ---: |
| Asia | 0.80984 | 58.816% | 60.144% | 58.971% |
| Africa | 0.25725 | 18.683% | 17.629% | 18.432% |
| Europe | 0.12950 | 9.405% | 9.105% | 9.171% |
| North America | 0.09919 | 7.204% | 7.220% | 7.518% |
| South America | 0.07276 | 5.284% | 5.254% | 5.344% |
| Oceania | 0.00836 | 0.607% | 0.648% | 0.564% |

Observed population shares are equal-time means of daily regional population / global population, baseline policy 0, all 1,460 days of 2022–2025. Regional counts use effective location, including travel. Reference shares divide the sum of 2023 and 2024 regional mid-year populations by the corresponding world sum. They are a nearby reference, not an exact four-year average. North America combines Northern America, Central America and the Caribbean. Source: [UN Population and Vital Statistics Report 2026, Table 1](https://unstats.un.org/unsd/demographic-social/products/vitstats/sets/Series_A_2026.pdf#page=7).

**How the choices were made.** The [9 September history](DEMOGRAPHIC_PARAMETER_REVIEW_2026-09-09.md) and both audits agree exactly with current source. First, all 18 weights in each region were multiplied by a rounded target-share / observed-share ratio from run 261633: Asia 1.06, Africa 1.05, Europe 0.74, North America 1.09, South America 1.07, Oceania 0.22. This was a reasonable proportional first correction, conditional on unchanged demographic response. The later age redistribution preserved those regional weight totals, but survival and birth timing can still change living regional shares.

| Choice group | Number of weights | Treatment and assessment |
| --- | ---: | --- |
| `[-40000,-36000)` | 6 | Held fixed. These births are outside the 2022–2025 calibration window; their numerical choices remain inherited priors. They need a forecast-specific justification. |
| `[-36000,-4000)`, eight bands | 48 | Fitted to the age distribution conditional on age below 80, retaining the regional total. The fit is reproducible, but does not independently determine eight birth-cohort sizes. |
| `[-4000,0)` | 6 | Held fixed because the legacy fitting baseline omitted this cohort. These scheduled births occur approximately 1930–1941 and contribute exclusively to 80+ in the calibration window. Their restored contribution requires validation. |
| `[0,32000)`, eight bands | 48 | Held fixed. These are the 1930 living population, not future births. The inspected provenance contains no separate census derivation for each inherited literal. |

The 48 fitted values use a frozen demographic-response approximation from run 103646, accounting explicitly for that run's original weights and legacy sampler. Their targets use WPP2024, both sexes, 1 July 2023, with full-continent North America. The complete geographic components sum to World at every single age in the archived extraction. The accepted [target file](calibration_snapshots/age_profile_first_pass_2026-09-09/age_targets.csv) supersedes the older Northern-America-only extract. [WPP2024](https://population.un.org/wpp/) supplies age-specific population estimates and projections.

For each region, eight fitted weights face only four independent constraints: three conditional age-share equations and the conserved total. Four dimensions are consequently selected by regularization, rather than independently observed age data. The [applied fitting script](calibration_snapshots/age_profile_first_pass_2026-09-09/fit_age_profiles.py#L59) explicitly checks rank four and penalizes squared relative departures plus differences between adjacent multipliers, using smoothness 0.25 and multiplier bounds 0.25–2.0. No archived fitted solution reaches a bound. The many decimal places preserve reproducibility and regional totals; they should not be interpreted as demographic precision.

**The largest current discrepancy is at age 80+.**

| Region | Run 071999, 2022–2025 | Run 545650, 2022–2025 | UN age reference, 2023 |
| --- | ---: | ---: | ---: |
| Asia | 3.148% | 3.132% | 1.765% |
| Africa | 1.521% | 1.525% | 0.458% |
| Europe | 7.011% | 6.966% | 5.405% |
| North America | 4.239% | 4.189% | 3.098% |
| South America | 2.862% | 2.894% | 1.945% |
| Oceania | 4.445% | 4.450% | 3.151% |

These are within-region shares of living people, not regional shares of the world. The source age groups are 0–5, 6–14, 15–49, 50–79 and 80+, exactly matching the model's age boundaries. Most younger-age shares are substantially closer to the reference. The discrepancy persists at model day 34128, approximately mid-2023: for example, Africa is 1.513% versus 0.458%, Asia 3.073% versus 1.765%, and Europe 6.848% versus 5.405% in run 071999. It is therefore not explained simply by comparing a four-year window with a 2023 reference.

The older-age fit should examine the fixed `[-4000,0)` cohort, the overlapping fitted `[-8000,-4000)` cohort, and survival together. At the start of 2025 these represent ages approximately 84–95 and 73–84 respectively. Reducing one old-cohort weight by the observed 80+ ratio would not be an identified correction: several cohorts contribute, survival intervenes, and global normalization affects other populations. No matching per-run parameter/seed metadata was located with these CSVs, so the outcome comparisons do not isolate the causal contribution of the demographic weights from other changes between runs.

**The untouched future cohort imposes an abrupt forecast assumption.** All bands have equal width and uniform sampling, so expected scheduled births per day are proportional to their weights. At model day 36000, around year 2028.63, births enter the unfitted `[-40000,-36000)` cohort. Relative to the preceding fitted cohort, expected daily birth counts change by:

| Region | Imposed change around 2028.63 |
| --- | ---: |
| Asia | +53.93% |
| Africa | −23.74% |
| Europe | +100.11% |
| North America | +61.68% |
| South America | +73.55% |
| Oceania | +27.20% |
| All regions | +32.67% |

This is a mathematical implication of the weights, not a measured forecast or a fertility rate per woman. It is absent from calibration output ending in 2025, but enters the longer policy horizon. It needs a population/birth projection basis before interpreting those later years. More generally, 4,000-day blocks imply piecewise-constant birth counts across roughly eleven-year periods. Smoothing adjacent fitted multipliers does not guarantee a smooth absolute birth schedule, particularly at fixed/fitted boundaries.

**The initially living weights need a historical justification if earlier years matter.** At the 1930 start they encode ages 0 to below about 87.7, with no initial 88+ tail. The repeated flat values in several older-age bands are inherited shapes rather than a demonstrated census reconstruction. By the start of 2025, all initially living cohorts would be at least 95 if surviving; most are much older. Modern broad age shares therefore cannot identify their original sizes. Evaluate these 48 choices against appropriate early population/age data and subsequent historical trajectories. WPP2024 starts in 1950, so a 1930 reconstruction requires another historical source or an explicitly documented initialization/burn-in assumption.

Only about 13.347% of allocated records are initially living under these weights. About 16.174% are still unborn at the end of 2025. These are expected sampling proportions before stochastic variation and deaths, not observed population estimates. Reserving future records is compatible with the model design, but changing an unobserved future weight still alters the global normalizer and therefore all other allocation probabilities. The final cohort runs out around 2039.59; projections beyond that need an explicit continuation of the birth schedule.

**Assessment by region.** Asia's living share remains high despite its normalized allocation share looking close to the reference. Africa's living share remains low and its 80+ excess is the largest relative discrepancy. Europe is close on total regional share, so another broad reduction has little support from this comparison; the stronger issue is the oldest-age mix. North America is modestly underrepresented, with the corrected full-continent target retained. South America is close on regional share but has an older-age excess. Oceania remains about 15% high relative to its small reference share and also has an oldest-age excess. These are directions for a joint demographic reassessment, not six independently justified new multipliers.

**Recommended order of work.** First, validate the full age distribution, including 80+, using current demographics and survival. Then jointly reassess all six living regional shares with a consistent observation window and geographic definition. Give the late-2028 onward cohort a projection-based justification, and inspect the complete birth schedule for discontinuities. Validate the 1930/early-history weights if historical outputs are used. Retain uncertainty/regularization assumptions alongside fitted values. Targets for these weights should be population, birth and age data.

Two smaller implementation observations are separate from the empirical choices: the sampler rebuilds the 108-element distribution for each person, and the current validation does not enforce finite nonnegative demographic weights with a strictly positive total. Neither invalid-weight condition occurs in the present configuration. One comment is stale: Africa's `neg4000_0` value 0.0105 is described as “Future births tapering”, although it rises from the adjacent fitted 0.008904529480.

Validation performed: exact checks of all 108 current literals against the applied age audit; exact regional-multiplier arithmetic and regional-sum preservation; sampler/key/boundary inspection; baseline daily population/age reads from three raw CSVs; a matched-date check from the two fresh subset caches; review of UN region definitions and demographic reference dates. No optimization, simulation run or production-code edit was performed.
