# Parameter lookup hotspots

This note records every place where the simulation pulls numeric parameters from the global `HashMap<String, f64>` that is built in `config.rs`. The goal is to make it clear which groups of lookups should be converted to indexed storage before we begin the refactor.

## Helper APIs and data sources

- `crate::config::PARAMETERS`: single `HashMap<String, f64>` holding all scalar, bacteria-, drug-, age-, and region-specific values.
- `get_global_param(&str) -> Option<f64>`: direct lookup by key string.
- `get_bacteria_param(bacteria: &str, suffix: &str)`: formats `"{bacteria}_{suffix}"` and then calls `PARAMETERS.get`.
- `get_drug_param(drug: &str, suffix: &str)`: formats `"drug_{drug}_{suffix}"` and then calls `PARAMETERS.get`.
- `ParameterKeyCache`: pre-computes string keys (e.g., `drug_bacteria_potency_keys`) but all call sites still fetch the value from `PARAMETERS` on every use.

## High-frequency lookups inside `rules::apply_rules`

`apply_rules` runs once per living individual per simulated day. Most CPU time is spent here, so any `HashMap` fetch in this function is automatically a hotspot.

### Per-individual scalars (executed once per person per day)

These are all read before any loops but still involve 20+ hash lookups per individual:

- Drug initiation knobs: `microbiome_resistance_transfer_probability_per_day`, `drug_base_initiation_rate_per_day`, `drug_infection_present_multiplier`, `already_on_drug_initiation_multiplier`, `drug_test_identified_multiplier`, `double_dose_probability_if_identified_infection`, `random_drug_cessation_probability`.
- Immunodeficiency dynamics: `temporary_immunosuppression_onset_rate_per_day`, `temporary_immunosuppression_recovery_rate_per_day`, `chronic_immunosuppression_onset_rate_per_day`, `chronic_immunosuppression_recovery_rate_per_day`, plus four age-banded probabilities (`chronic_immunodeficiency_probability_age_*`).
- Hospitalization block: `hospitalization_baseline_rate_per_day`, `hospitalization_age_multiplier_per_day`, `hospitalization_recovery_rate_per_day`, `hospitalization_max_days`, `hospitalization_sepsis_admission_multiplier`, `hospitalization_prevent_discharge_with_sepsis`.
- Travel frequency: `travel_probability_per_day`.

Every one of these is a single string lookup that could live in a dedicated struct of cached scalars.

### Hospitalisation and sepsis recovery

- Region-specific sepsis adjustment uses `param_cache.region_sepsis_keys[region]` → `get_global_param(key)` each timestep an individual has sepsis.
- Recovery odds pull the following on every evaluation: `sepsis_base_log_odds_of_recovery_per_day`, `sepsis_log_odds_bacteria_level`, `sepsis_log_odds_in_hospital`, `sepsis_log_odds_age_infant/child/adult/elderly`, and `sepsis_log_odds_immunosuppressed`.

### Region travel block (per individual per day)

- Region multiplier fetched as `param_cache.region_travel_keys[region]` → `get_global_param`.
- Later in the travel block, several regional acquisition multipliers are read through `param_cache.region_bacteria_acquisition_keys[(region,bacteria)]` and fall back to `param_cache.region_bacteria_default_keys[region]`.

### Drug selection & initiation (per individual × drug candidate)

This is the heaviest consumer of parameter lookups:

- For every drug (length of `DRUG_SHORT_NAMES`), we fetch potency and initiation multipliers via `param_cache.drug_bacteria_potency_keys[(drug,bacteria)]` and `param_cache.drug_bacteria_initiation_keys[(drug,bacteria)]`.
- Thresholds and weights for empirical vs targeted therapy: `minimal_potency_threshold_for_drug_selection`, `effective_potency_threshold_for_targeted_therapy`, `effective_potency_threshold_for_empirical_therapy`, `targeted_therapy_narrow_spectrum_bonus`, `targeted_therapy_broad_spectrum_penalty`, `targeted_therapy_ineffective_drug_penalty`, `empiric_therapy_broad_spectrum_bonus`, `empiric_therapy_ineffective_drug_penalty`, `empiric_therapy_low_spectrum_penalty` (fallback via inline `0.6`).
- Regional resistance penalties rely on `regional_resistance_threshold_{very_high,high,moderate}` and `regional_resistance_penalty_{very_high,high,moderate}`—all global lookups inside nested loops over bacteria.
- Stochastic selection uses `drug_selection_temperature` every time a drug is picked.
- Drug properties: `get_drug_param(drug, "initial_level")`, `get_drug_param(drug, "double_dose_multiplier")`, `get_drug_param(drug, "spectrum_breadth")`.
- `get_drug_availability_time_aware` ultimately calls `PARAMETERS.get` for keys like `"{region}_drug_{name}_availability"` each iteration.

### Bacteria acquisition & microbiome update (per individual × bacteria)

For every bacteria in `BACTERIA_LIST` (~30 entries):

- Fetch `get_bacteria_param(bacteria, "acquisition_log_odds_baseline")` with fallback to `acquisition_log_odds_baseline` global.
- Age category adjustments build keys `"{bacteria_clean}_log_odds_{age_category}"` and fall back to `"default_log_odds_{age_category}"`.
- Additional modifiers: `log_odds_vaccinated`, `log_odds_microbiome_present`, `log_odds_hospital_acquired`, `log_odds_microbiome_vs_infection` (mix of bacteria and global keys).
- Age × region interaction keys: `"{bacteria_clean}_{region}_log_odds_{age_category}"` with fallback to region-only key `"{region}_log_odds_{age_category}"`.
- Historical MDR-TB multipliers: `mdr_tb_pre_antibiotic_era_multiplier`, `mdr_tb_early_antibiotic_era_multiplier`, `mdr_tb_modern_era_multiplier`.
- Microbiome acquisition/clearance uses `environmental_majority_r_level_for_new_acquisition`, `max_resistance_level`, `default_microbiome_clearance_probability_per_day`, `microbiome_resistance_emergence_rate_per_day_baseline`, `any_r_emergence_level_on_first_emergence`.

### Resistance transfer & HGT (nested over donor × recipient × drug)

- Horizontal gene transfer probability uses cached keys `param_cache.hgt_keys[(donor,recipient)]` → `PARAMETERS.get` on every pair.
- Mechanism assignment consults `mechanism_assignment_probability_on_any_r_gain` and `param_cache.resistance_mechanism_enhancement_keys[mechanism]` for every mechanism and drug whenever resistance transfers.

### Drug toxicity and pharmacodynamics (per individual × active drug)

- `get_drug_param(drug, "toxicity_per_unit_level_per_day")` with fallback to `default_drug_toxicity_per_unit_level_per_day`.
- Bounding and clearance use `max_toxicity_level` and `toxicity_clearance_rate_per_day`.

### Immune response updates (per individual × bacteria)

- Multiple bacteria-specific keys retrieved each day: `base_bacteria_level_change`, `immunity_effect_on_level_change`, `max_level`, `immunity_base_response`, `immunity_increase_per_infection_day`, `immunity_increase_per_unit_higher_bacteria_level`, `immunity_age_modifier`, `immunity_immunodeficiency_modifier`, `max_immune_response`.
- Global fallback keys: `immune_decay_rate_per_day`, `immune_response_recovery_rate_per_day` (later in the file).

## Other modules

### `simulation::Simulation::new`

- Precomputes potency matrix by looping over all bacteria × drug pairs and performing `PARAMETERS.get(&format!("drug_{drug}_for_bacteria_{bacteria}_potency_when_no_r"))`.
- Also computes MIC thresholds using the same retrieved potency values.

### `simulation::Simulation::run`

- Summary metrics rely on `drug_evaluation_days_post_infection`, and the same value is fetched twice per timestep.
- Additional diagnostic/reporting code pulls a handful of global scalars (e.g., evaluation-day thresholds) within loops over the entire population; these are less performance critical but still candidates for caching.

## Parameter families to migrate to indexed storage

| Parameter family | Key pattern today | Primary call sites | Loop frequency | Natural indexed shape |
| --- | --- | --- | --- | --- |
| Drug potency vs bacteria | `drug_{drug}_for_bacteria_{bacteria}_potency_when_no_r` | `Simulation::new`, `apply_rules` drug scoring | per individual × drugs × infected bacteria | `[bacteria][drug]` matrix of `f32/f64` |
| Drug initiation multipliers | `drug_{drug}_for_bacteria_{bacteria}_initiation_multiplier` | Drug selection | per individual × drugs × infected bacteria | `[bacteria][drug]` |
| HGT probabilities | `hgt_prob_{donor}_to_{recipient}` | HGT loop | per individual × bacteria² | `[donor][recipient]` |
| Region × bacteria acquisition | `{region}_{bacteria}_acquisition_log_odds` | Infection & microbiome acquisition | per individual × bacteria | `[region][bacteria]` |
| Age × bacteria log odds | `{bacteria}_log_odds_{age}` | Infection & microbiome acquisition | per individual × bacteria | `[bacteria][age_category]` |
| Age × region log odds | `{region}_log_odds_{age}` | Infection & microbiome acquisition fallback | per individual × bacteria | `[region][age_category]` |
| Resistance mechanism enhancement | `resistance_mechanism_{mechanism}_enhancement_multiplier` | HGT + mechanism assignment | per individual × bacteria × drug × mechanisms | `[mechanism]` |
| Drug toxicity and pharmacokinetics | `drug_{drug}_toxicity_per_unit_level_per_day`, `drug_{drug}_half_life_days` | Drug decay/toxicity | per individual × active drug | `[drug]` |
| Global scalars | direct keys (see list above) | Everywhere | per individual | Struct of cached scalars |

## Notes for the upcoming refactor

1. The `ParameterKeyCache` already maps (index tuples) → string keys; replacing the backing `HashMap<String, f64>` with structured arrays will let us swap those maps for pure index arithmetic.
2. Almost every hot path currently formats or clones strings before looking up `PARAMETERS`. Eliminating the `HashMap` will remove both allocations and hash work.
3. Many fallback patterns (e.g., bacteria-specific key → default) suggest layered storage: a dense array plus optional overrides. We can model this with per-family default values that are applied during config load time.
4. Once the indexed data structures are in place, helper functions like `get_global_param` can become thin wrappers over the arrays for non-performance-critical code, easing the migration.
