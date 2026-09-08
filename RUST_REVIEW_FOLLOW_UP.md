# Rust review follow-up checklist

Recorded: 7 September 2026. Reviewed commit: `e85012474dbd9e213d0f3e5c82984f9dee950619`.

This checklist records the review and subsequent follow-up work. The initial review made no Rust implementation changes; later completion notes are recorded under the corresponding items. Source line numbers refer to the reviewed revision; recheck them before starting. Existing uncommitted Python configuration and paper-output changes were present during the review.

The priorities below distinguish confirmed correctness defects from confirmed behaviours whose intended modelling semantics still need a decision. Python-only findings from the wider review are outside this checklist.

1. [x] **High priority: collect resistance numerators and infection denominators from the same population snapshot.**

   **Review evidence:** In the reviewed [simulation.rs](src/simulation/simulation.rs), resistance-positive counts were collected before `apply_rules()` (around lines 4964-5008), while infection denominators used post-rule survival and infection state (around lines 5695, 5911 and 5930-6004). The hospital classification in both paths was the pre-rule hospital status; the definite mismatch was the observation time and eligible population.

   **Why change it:** A person can contribute to the resistance numerator and then clear infection or die before contributing to its denominator. The stored schema-3 `simulation_summary_368577.csv` contains 27 resistant S. epidermidis/erythromycin infections against 26 current infections in its first row. Some hospital ratios also exceed 100% across the full window. The exact source revision producing that stored CSV was not established, so its magnitude must not be attributed entirely to the current code without checking provenance.

   **Implementation (schema 5):** Global and hospital/community drug-specific resistance-positive counts and `any_r` sums now use the same post-rule living, born, active-infection population as their denominators. The associated MIC-proxy, infection-mechanism and mechanism-family, infected-on-drug, and identified/AST-ready infection stocks use that snapshot too. Hospital/community fields preserve the pre-rule hospital classification for both numerator and denominator. Schema-4 regional fields retain their post-rule home-region definition; their totals now reconcile with the global resistance fields when both are collected. Resistance-profile sampling, its random draws, microbiome-only reporting and within-day applied-activity observations retain their existing timing, preserving model dynamics.

   **Historical compatibility:** Python accepts schema 5 alongside the existing supported versions. Calibration summaries from schemas 1–4 still run but carry an explicit warning about incompatible global resistance observation times. Schema-5 calibration rejects inconsistent counts and sums rather than silently clipping percentages. Historical CSVs cannot be repaired from their aggregate fields; rebuild and rerun to obtain corrected global results.

   **Verified:** Two new Rust regression tests compare raw statistics with the actual final population across full, standard-calibration and minimal-calibration profiles, and explicitly exercise same-day clearance and death. They check positive-count bounds, hospital/community partitions, regional/global counts and sums, testing stocks, mechanism counts and susceptibility counts. Seventeen new Python tests enforce schema-5 arithmetic and preserve historical handling. A schema-5 Rust export passed the real Python benchmark and regional-table calculations. A separate comparison with the previous compiled binary matched all 64,703 non-version fields in the unchanged short-run fixture. Full Rust and Python suites retained only the three previously recorded parameter/whitespace failures.

2. [ ] **High priority: propagate checkpoint and policy-branch failures to the launcher.**

   **Evidence:** [Simulation::run](src/simulation/simulation.rs), around lines 6835 and 6866-6899, returns `()` and only prints baseline/branch errors. [main.rs](src/main.rs), around lines 416 and 431-435, then exports summaries and can record `completed` if export and hashing succeed.

   **Why change it:** A failed checkpoint write, checksum check or restore can leave truncated output or missing policy branches with successful completion metadata.

   **Proposed work:** Return a `Result` from the run path and propagate failure through the launcher to an unsuccessful process exit and truthful run metadata. Mark completion only after all requested trajectories finish. If partial outputs are retained, identify them explicitly as incomplete.

   **Verification:** Force checkpoint-write and checkpoint-restore failures, including failure after an earlier policy succeeded. Check exit status, metadata, branch coverage and cleanup. Keep successful fixed-seed checkpoint parity tests passing.

3. [x] **High priority: activate the agreed infection/carriage exchange pathway (implemented and regression-verified).**

   **Review evidence:** In the reviewed revision, [rules/mod.rs](src/rules/mod.rs), line 5201, required no infection episode. Inside that block, the transfer guard around line 5526 required a positive infection level. Comments at line 5522 explicitly acknowledged that the pathway was unreachable. [MODEL_DESCRIPTION.md](model_description/MODEL_DESCRIPTION.md), line 1872, nevertheless described it as operative.

   **Accepted contract and implementation:** The existing exchange probability now applies at the start of each organism's daily within-host update, before new acquisition and ordinary infection progression. An already-present positive infection episode, including a fading positive level, and separate carriage of the same organism are both required. When their eligible profiles differ, at most one draw shares the union of their existing eligible mechanisms and raises the corresponding infection/carriage resistance measures where needed. Exchange creates neither compartment, changes no episode dates, and does not itself promote the predominant-strain profile. Equal profiles, an absent compartment, or a zero effective probability make no state change and consume no random draw. Counterfactual scaling is retained, including complete suppression at a zero multiplier. A newly established infection uses its separate inheritance route on the acquisition day and first becomes eligible for exchange the following day.

   **Scope and impact:** This is an intentional change to model dynamics, separate from the schema-5 reporting correction in item 1. The existing exchange parameter is unchanged, and item 4's ordinary carriage-update suspension is retained. Its effect on calibration has not yet been measured; fresh simulations are needed before interpreting new calibration results.

   **Verified (8 September 2026):** All eight focused exchange tests pass. They cover eligible unions, compartment/state preservation, active and fading episodes (including the exact activity threshold), absent/equal/disabled cases, rejected draws, host filtering, resistance values that must not decrease, and execution through `apply_rules()` with carriage still frozen. The once-only call and its position before new acquisition are checked against the source. `cargo test --all-targets --no-fail-fast` reports 215 passes and only the two existing parameter-assertion failures listed in item 6; fixed-seed checkpoint parity and thread determinism pass. `cargo test --doc` has no failures and one ignored test. Python discovery reports 207 passes and the existing whitespace-sensitive `test_current_launcher_population_and_horizons_are_reported` failure. Local logs: [Rust](target/exchange_rust_verification.log), [doctests](target/exchange_rust_doctests.log), [Python](target/exchange_python_verification.log). The calibration impact assessment remains outstanding.

4. [ ] **High priority, modelling decision required: resolve the suspension of ordinary carriage dynamics during infection.**

   **Review evidence:** The same `!has_infection_episode` block in the reviewed [rules/mod.rs](src/rules/mod.rs), starting at line 5201, enclosed ordinary carriage acquisition, clearance, reversion and emergence. Those updates were skipped throughout a same-organism infection, including a fading positive episode. Later HGT and possible carriage clearance at infection resolution were exceptions. The description around [line 511](model_description/MODEL_DESCRIPTION.md) did not describe this suspension.

   **Proposed work:** Explicitly decide which carriage processes continue while that organism also causes an infection, then align the control flow and documentation. Review this together with item 3: repairing the transfer block alone would leave the broader suspension in place.

   **Scope decision for item 3:** Ordinary carriage acquisition, clearance, emergence and reversion intentionally remain suspended throughout a positive same-organism infection episode, including fading positive levels. Same-organism exchange is an explicit exception alongside the existing HGT and drug-assisted resolution pathways. The model description now records this restriction. Whether to change the broader suspension remains an open modelling decision; it is outside the agreed exchange change.

   **Verification:** Compare controlled carriage-only, active-infection-plus-carriage and fading-infection-plus-carriage states. Verify each process against the agreed contract and check for duplicate updates. Measure any resulting changes in carriage duration and calibration.

5. [ ] **Medium priority, modelling decision required: review treatment-failure clock resets on unchanged-drug reselection.**

   **Evidence:** [rules/mod.rs](src/rules/mod.rs), around lines 4697-4707, resets tracking for every positive infection after every successful drug selection, explicitly including reselection of an already-active drug. `mark_new_treatment_course()`, around lines 1462-1473, resets treatment days and the assessment flag and resamples the response multiplier.

   **Why review it:** Continuous unchanged treatment can repeatedly restart the failure-assessment clock. A selection associated with one infection also resets tracking for concurrent infections. This behaviour is explicit in comments, so it is a design concern rather than a proven accidental regression.

   **Proposed work:** Define when a course is genuinely new. The recommended default is to preserve duration and response state on unchanged-drug reselection. Decide separately how additions, switches, restarts and concurrent infections should affect tracking.

   **Verification:** Distinguish repeated unchanged treatment from a genuine switch or restart, and include concurrent infections. Confirm that failure assessment occurs after the intended duration under the chosen definition.

6. [ ] **Medium priority: reconcile the two failing Rust parameter-contract tests.**

   **Observed failures:** `supported_vaccine_targets_use_organism_specific_acquisition_parameters` in [rules/mod.rs](src/rules/mod.rs), around line 9976, expects an H. influenzae acquisition baseline of `-18.47`; the configuration provides `-16.5`. `enteric_incidence_recalibration_preserves_carriage_intercepts` in [config_invariants.rs](tests/config_invariants.rs), around line 341, expects an E. coli acquisition baseline of `-11.353`; the configuration provides `-11.4`.

   **Proposed work:** Establish whether the parameter changes were intentional, then reconcile configuration, documentation and tests. Do not revert calibrated values just to pass a test, or blindly update assertions. Where appropriate, distinguish tests of biological/output contracts from assertions that freeze a particular calibration snapshot.

   **Verification:** Run the complete Rust suite after reconciliation and record any intentionally changed parameter expectations.

7. [ ] **Medium priority: bring the parameter reference back into agreement with Rust.**

   **Evidence:** Generating Appendix B with [dump_parameter_appendix.rs](src/bin/dump_parameter_appendix.rs) and comparing it in memory against [MODEL_DESCRIPTION.md](model_description/MODEL_DESCRIPTION.md) found 189 differing lines. For example, the documented E. coli sepsis baseline is `-9.7`; the current Rust configuration generates `-9.4`.

   **Proposed work:** Once the intended configuration is settled, regenerate the appendix using the existing generator/update script and reconcile explanatory prose. Add a check for unintended drift between generated parameter documentation and the executable configuration. This is a Rust-to-documentation consistency task, not a reason to retune parameters to old prose.

   **Verification:** The generated appendix should match the checked-in reference, and the descriptions of items 3-5 should match the accepted behaviour.

8. [ ] **Medium priority: ensure CI executes the full Rust regression coverage.**

   **Evidence:** [fast-guardrails.yml](.github/workflows/fast-guardrails.yml) checks all targets but executes only six named integration-test suites. It does not run the library unit tests or all other integration suites, so some existing failures and future regressions can escape this workflow.

   **Proposed work:** Add complete Rust test execution, either within this workflow or in an additional job with an appropriate runtime budget. Retain the fast targeted checks where useful.

   **Verification:** Confirm that library tests and every integration target actually run in CI and that failed assertions make the workflow fail.

9. [x] **Requested output addition: export resistance counts and sums by geographical region.**

   **Purpose:** Support a regional version of the calibration summary's "Simulation mean (%)" for infection resistance and "Resistant level (among positives)", and provide inputs for future regional figures. Prior to schema 4, the summary CSV lacked the necessary regional resistance counts and sums.

   **Implemented:** Schema 4 introduced `regional_resistance_collected` and the three field families below. All use home residence and the same end-of-day surviving active-infection snapshot. The existing `regional` content flag enables them in the standard calibration and full-history modes; `FullMinimal` leaves collection disabled. Python adds a six-row "Regional Resistance Summary" with prevalence, conditional severity and contributing-pair counts. Schemas 3, 4 and 5 remain supported for general analysis; older or disabled regional data is explicitly unavailable. This is reporting only and does not alter model dynamics or the calibration score. Item 1 records the schema-5 correction to global resistance timing.

   **Required CSV quantities:** For each retained timestep and policy, export:

   - Active infected counts for each bacterium and region.
   - Infected counts with `any_r > 0` for each bacterium, drug and region.
   - Summed infection-level `any_r` for each bacterium, drug and region, to calculate conditional resistance severity.

   **Implementation locations:** [simulation.rs](src/simulation/simulation.rs) collects and exports the snapshots; [calibration_summary.py](amr_simulation_output_analysis/calibration_summary.py) calculates and writes the table. Regional acquisition counts and circulating resistance-profile caches are different quantities and are not used as substitutes. The exact column patterns and mode coverage are documented in [MODEL_DESCRIPTION.md](model_description/MODEL_DESCRIPTION.md), Appendix C.

   **Downstream calculation:** Over the baseline calibration window, calculate each bacterium-drug-region prevalence as `100 × sum(positive counts) / sum(infected counts)` and conditional severity as `100 × sum(any_r) / sum(positive counts)`. To match the existing "Simulation mean (%)", average the eligible bacterium-drug percentages with equal weight within each region; do not replace that statistic with a pooled fraction across all infections. Handle zero denominators as missing observations.

   **Verification:** Tests cover home versus visited region, living/active eligibility, regional totals against post-rule infection stocks, count/sum consistency, disabled collection, unchanged model state, export layout and invalid array lengths. Python tests cover sum-over-window ratios, equal pair weighting, eligibility, missingness and policy/window filtering. Regional positive counts and `any_r` sums share the global post-rule population from schema 5; schema-4 regional resistance must not be forced to reconcile with its older global pre-rule fields. Existing CSV files cannot recover a missing breakdown; a new simulation run is required for those regional research results.

The review's baseline was 157 passing library tests and one failing library test, plus a separate failing configuration-invariant test. The other exercised Rust targets passed. These are existing failures to reconcile, not evidence of changes introduced by future repairs.

Suggested order for the remaining review work: address failure propagation; resolve the outstanding carriage and treatment-tracking decisions before changing those pathways; then reconcile parameter tests, maintain generated documentation and complete CI coverage. Item 3 records the separately agreed exchange change. Preserve separate reviewable changes for reporting repairs and changes to model dynamics.
