**Suggested demographic trial, 12 September 2026 — not applied.**

For the next calibration run through 2025, reduce the six `demo_<region>_age_neg4000_0` weights as follows. This cohort contributes exclusively to the 80+ population during 2022–2025. Its weight is only one contributor to that age group's eventual population.

| Region | Current weight | Proposed weight | Change |
| --- | ---: | ---: | ---: |
| Asia | 0.028620 | 0.022896 | −20% |
| Africa | 0.010500 | 0.007875 | −25% |
| Europe | 0.006660 | 0.005994 | −10% |
| North America | 0.003270 | 0.0027795 | −15% |
| South America | 0.002140 | 0.001819 | −15% |
| Oceania | 0.000220 | 0.000187 | −15% |

These are judgment-based trial magnitudes, deliberately limited to 10–25%. Africa and Asia receive the larger changes because their oldest-age excess is stronger relative to the reference. They are not target/observed ratios applied to a supposedly identified cohort response, and they are not expected percentage reductions in the whole 80+ population.

Redistribute the removed weight proportionally across the eight younger birth cohorts in the same region, from `age_neg36000_neg32000` through `age_neg8000_neg4000`. For a region with removed mass D and current combined weight S in those eight cohorts, multiply each by `(S + D) / S`.

| Region | Multiplier for each of those eight weights |
| --- | ---: |
| Asia | 1.010000000000 |
| Africa | 1.014534883721 |
| Europe | 1.008256880734 |
| North America | 1.007500000000 |
| South America | 1.006521739130 |
| Oceania | 1.006818181818 |

The [complete proposed values](proposed_weights.csv) retain enough precision to preserve every regional total exactly. After rounding redistributed weights to 12 decimal places, any rounding residual is added to the largest redistributed cohort. The displayed multipliers are explanatory rounded values; use the CSV or patch for the exact proposal.

This changes 54 of the 108 parameters: six older-cohort reductions and 48 small proportional increases. The other 54 parameters are fixed: all eight initially living cohorts per region and the six far-future `neg40000_neg36000` cohorts. All regional allocation shares, the global total of 1.3769 and the initial 1930 living allocation remain unchanged. This preserves the relative weights of the cohorts contributing to under-80 ages in the calibration window. It does not guarantee unchanged simulated age shares, living regional shares or living population size.

I would assess the age response before adding another regional multiplier adjustment. The [previous review](../../DEMOGRAPHIC_WEIGHT_REVIEW_2026-09-12.md) identifies the regional residuals and the separate late-2028 birth discontinuity. This proposal is for the current calibration horizon; it does not resolve that forecast assumption. Changing the untouched future cohort from 2025 observations alone would lack a projection basis and would also change the sampling normalizer unless compensated elsewhere.

The current input file `simulation_summary_638391.csv` contains a complete baseline 2022–2025 window and shows the same qualitative demographic pattern as the two earlier files. Its embedded `run_id` is 479654. The earlier filenames also differ from their embedded IDs, and no matching local run/configuration metadata was found. The observations support the direction of this trial but do not establish a securely matched parameter-response fit. The source weights themselves match all 108 values in the accepted September 9 audit.

Why an exact new fit is not claimed: the existing broad-age response approximation uses one coefficient for all contributors to 80+, mixing very old initially living cohorts with later births. It cannot separately identify their survival or justify exact cohort reductions. An independent arithmetic check found that proportional redistribution does not reverse the intended direction in that approximation and preserves its conditional under-80 mix. Those are mathematical checks, not results from a new simulation.

Files supplied:

- [proposed_weights.csv](proposed_weights.csv): all 108 current/proposed values and reasons.
- [regional_changes.csv](regional_changes.csv): six principal changes, redistribution factors and exact regional totals.
- [proposed_section_k.rs](proposed_section_k.rs): the full suggested section K.
- [section_k.patch](section_k.patch): a patch against the reviewed working configuration; it contains only the demographic changes and a trial note.
- [proposal_manifest.json](proposal_manifest.json): source hash, scope and assumptions.

Validation: all weights remain finite and positive; exactly 54 change; the 54 fixed weights are identical; all six regional totals and the global total are preserved exactly in decimal arithmetic. Source comparison confirms that `src/config.rs` was not modified. No simulation, optimization or disease-parameter change was performed.
