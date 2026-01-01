// Centralized configuration and parameter management for the AMR simulation.
//
// Contains:
//   - Initialization of global, bacteria-specific, and drug-specific parameters
//   - Functions for parameter lookup and cross-resistance group management
//   - Age-specific vaccination, HGT, and other model parameters
//   - Reference for template and override logic
//
// ===================== PARAMETER DEFAULTS QUICK INDEX =====================
//   A) Bacteria baseline defaults & vaccination scaffolding ........... ~4140
//   B) Horizontal gene transfer (HGT) priors ........................... ~4185
//   C) Drug initiation, selection, and pharmacokinetics ................ ~4199
//   D) Drug interaction adjustments & therapy flow ..................... ~4356
//   E) Drug-bacteria potency & emergence settings ...................... ~4381
//   F) Regional acquisition pressure baselines ........................ ~5556
//   G) Microbiome carriage & clearance priors .......................... ~5836
//   H) Resistance emergence, transfer, and mechanism weights ........... ~5922
//   I) Clinical outcome scalars (mortality, sepsis, toxicity) .......... ~6186
//   J) Regional availability & introduction timelines .................. ~6834
//   K) Demographic distribution defaults .............................. ~6950
// ========================================================================
// ===================== READER / LOOKUP STRUCTURE =======================
//   1) Core indices & constants ........................................ ~46
//   2) Parameter store & struct definitions ............................ ~91
//   3) Global scalar readers (from_map accessors) ...................... ~153
//   4) Immunodeficiency / region / syndrome / sex readers .............. ~747
//   5) Vaccination & age category readers .............................. ~880
//   6) Drug parameter blocks (reader structs) .......................... ~1096
//   7) Bacteria-level readers & microbiome logic ....................... ~1210
//   8) Clearance, acquisition, and age readers ......................... ~1466
//   9) Drug-bacteria matrices & cross-resistance readers ............... ~4385
//  10) Regional drug availability & introduction lookups ............... ~6837
//  11) Demographic distribution readers ................................ ~6953
//  12) Helper lookups (drug intro, availability, etc.) ................. ~7125
//  13) HGT matrices & resistance mechanism readers ..................... ~7266
// ========================================================================
// Sections A-K highlight where defaults are inserted. Sections 1-13 show the
// typed readers that consume those defaults (falling back to literals only
// when a configuration key is absent).

// src/config.rs
use crate::simulation::population::{
    AgeCategory,
    AGE_CATEGORY_SEQUENCE,
    Region,
    ResistanceMechanism,
    BACTERIA_LIST,
    DRUG_SHORT_NAMES,
};
use lazy_static::lazy_static;
use std::borrow::Cow;
use std::collections::HashMap; // Import both lists and helper enums

// ---------------- 1) Core indices & constants ----------------

lazy_static! {
    pub static ref BACTERIA_INDEX: HashMap<&'static str, usize> = {
        BACTERIA_LIST
            .iter()
            .enumerate()
            .map(|(idx, &name)| (name, idx))
            .collect()
    };
}

fn is_canonical_bacteria_slug(name: &str) -> bool {
    BACTERIA_INDEX.contains_key(name)
}

pub fn canonicalize_bacteria_slug<'a>(name: &'a str) -> Cow<'a, str> {
    let trimmed = name.trim();
    if trimmed.is_empty() {
        panic!("Bacteria name cannot be empty");
    }

    if is_canonical_bacteria_slug(trimmed) {
        return Cow::Borrowed(trimmed);
    }

    panic!(
        "Unknown bacteria name '{name}'. Use canonical slugs from BACTERIA_LIST (e.g., acinetobacter_baumannii)",
        name = name
    );
}

const REGION_NAMES: [&str; 6] = [
    "north_america",
    "south_america",
    "africa",
    "asia",
    "europe",
    "oceania",
];

const HOME_REGION_NAME: &str = "home";

const REGION_VARIANTS: [Region; 7] = [
    Region::NorthAmerica,
    Region::SouthAmerica,
    Region::Africa,
    Region::Asia,
    Region::Europe,
    Region::Oceania,
    Region::Home,
];

const AGE_BUCKETS: [AgeCategory; AGE_CATEGORY_SEQUENCE.len()] = AGE_CATEGORY_SEQUENCE;
const AGE_BUCKET_COUNT: usize = AGE_BUCKETS.len();

// ---------------- 2) Parameter store & struct definitions ----------------
#[derive(Debug)]
pub struct ParameterStore {
    pub globals: GlobalScalars,
    pub immunodeficiency: ImmunodeficiencyParameters,
    pub region: RegionParameters,
    pub syndrome: SyndromeParameters,
    pub sex: SexParameters,
    pub vaccination: VaccinationParameters,
    pub drug: DrugParameters,
    pub bacteria: BacteriaParameters,
    pub clearance: ClearanceParameters,
    pub drug_bacteria: DrugBacteriaMatrix,
    pub region_bacteria: RegionBacteriaAcquisition,
    pub age_categories: AgeCategoryParameters,
    #[allow(dead_code)]
    pub age_tables: AgeTables,
    pub hgt: HgtMatrix,
    pub resistance_mechanism: ResistanceMechanismParameters,
}

impl ParameterStore {
    fn from_parameter_map(map: &HashMap<String, f64>) -> Self {
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();

        let globals = GlobalScalars::from_map(map);
        let immunodeficiency = ImmunodeficiencyParameters::from_map(map);
        let region = RegionParameters::from_map(map);
        let syndrome = SyndromeParameters::from_map(map);
        let sex = SexParameters::from_map(map);
        let vaccination = VaccinationParameters::from_map(map);
        let drug = DrugParameters::from_map(map, num_drugs);
        let bacteria = BacteriaParameters::from_map(map, num_bacteria);
        let clearance = ClearanceParameters::from_map(map, num_bacteria);
        let drug_bacteria = DrugBacteriaMatrix::from_map(map, num_bacteria, num_drugs);
        let region_bacteria = RegionBacteriaAcquisition::from_map(map, num_bacteria);
        let age_categories = AgeCategoryParameters::from_map(map, num_bacteria);
        let age_tables = AgeTables::from_map(map, num_bacteria);
        let hgt = HgtMatrix::from_map(map, num_bacteria);
        let resistance_mechanism = ResistanceMechanismParameters::from_map(map);

        ParameterStore {
            globals,
            immunodeficiency,
            region,
            syndrome,
            sex,
            vaccination,
            drug,
            bacteria,
            clearance,
            drug_bacteria,
            region_bacteria,
            age_categories,
            age_tables,
            hgt,
            resistance_mechanism,
        }
    }
}

// ---------------- 3) Global scalar defaults & helpers ----------------
#[derive(Debug)]
pub struct GlobalScalars {
    pub drug_base_initiation_rate_per_day: f64,
    pub drug_infection_present_multiplier: f64,
    pub drug_activity_to_bacteria_level_multiplier: f64,
    pub already_on_drug_initiation_multiplier: f64,
    pub drug_test_identified_multiplier: f64,
    pub double_dose_probability_if_identified_infection: f64,
    pub random_drug_cessation_probability: f64,
    pub random_drug_cessation_probability_if_no_active_infection: f64,
    pub immunodeficiency_prophylactic_drug_multiplier: f64,
    pub microbiome_resistance_transfer_probability_per_day: f64,
    pub hospital_baseline_rate_per_day: f64,
    pub hospital_age_multiplier_per_day: f64,
    pub hospital_recovery_rate_per_day: f64,
    pub hospital_max_days: f64,
    pub hospital_sepsis_admission_multiplier: f64,
    pub hospital_prevent_discharge_with_sepsis: f64,
    pub travel_probability_per_day: f64,
    pub antibiotic_infection_prevention_efficacy: f64,
    pub max_resistance_level: f64,
    pub resistance_emergence_bacteria_level_multiplier: f64,
    /// Debug multiplier to keep 2025 resistance prevalence stable when population size changes
    pub resistance_emergence_pop_size_multiplier: f64,
    pub any_r_emergence_level_on_first_emergence: f64,
    pub multi_drug_penalty_threshold_num_drugs: f64,
    pub resistance_development_inhibition_single_drug: f64,
    pub resistance_development_inhibition_partial_cross: f64,
    pub mechanism_assignment_probability_on_any_r_gain: f64,
    pub treatment_failure_enabled: bool,
    pub treatment_failure_assessment_day: i32,
    pub treatment_failure_threshold: f64,
    pub drug_failure_memory_days: i32,
    pub minimal_potency_threshold_for_drug_selection: f64,
    pub drug_selection_temperature: f64,
    pub reserve_drug_score_penalty: f64,
    pub restart_window_enabled: bool,
    pub restart_window_days: i32,
    pub restart_bacteria_level_threshold: f64,
    pub restart_window_probability: f64,
    pub log_odds_sepsis_region_a: f64,
    pub log_odds_sepsis_region_b: f64,
    pub mdr_tb_pre_antibiotic_era_multiplier: f64,
    pub mdr_tb_early_antibiotic_era_multiplier: f64,
    pub mdr_tb_modern_era_multiplier: f64,
    pub microbiome_resistance_emergence_rate_per_day_baseline: f64,
    pub default_toxicity_reservoir_half_life_days: f64,
    pub toxicity_age_multiplier_infant: f64,
    pub toxicity_age_multiplier_child: f64,
    pub toxicity_age_multiplier_adult: f64,
    pub toxicity_age_multiplier_elderly: f64,
    pub toxicity_immunosuppressed_multiplier: f64,
    pub toxicity_hospital_multiplier: f64,
    pub regional_resistance_threshold_very_high: f64,
    pub regional_resistance_threshold_high: f64,
    pub regional_resistance_threshold_moderate: f64,
    pub regional_resistance_penalty_very_high: f64,
    pub regional_resistance_penalty_high: f64,
    pub regional_resistance_penalty_moderate: f64,
    pub targeted_therapy_narrow_spectrum_bonus: f64,
    pub targeted_therapy_broad_spectrum_penalty: f64,
    pub targeted_therapy_ineffective_drug_penalty: f64,
    pub effective_potency_threshold_for_targeted_therapy: f64,
    pub empiric_therapy_broad_spectrum_bonus: f64,
    pub empiric_therapy_ineffective_penalty: f64,
    pub any_r_increase_rate_per_day_when_drug_present: f64,
    pub sepsis_minimum_duration_days: i32,
    pub sepsis_base_log_odds_of_recovery_per_day: f64,
    pub sepsis_log_odds_bacteria_level: f64,
    pub sepsis_log_odds_in_hospital: f64,
    pub sepsis_log_odds_age_infant: f64,
    pub sepsis_log_odds_age_child: f64,
    pub sepsis_log_odds_age_adult: f64,
    pub sepsis_log_odds_age_elderly: f64,
    pub sepsis_log_odds_immunosuppressed: f64,
    pub infection_non_sepsis_base_log_odds: f64,
    pub infection_non_sepsis_log_odds_per_level: f64,
    pub infection_non_sepsis_log_odds_age_infant: f64,
    pub infection_non_sepsis_log_odds_age_child: f64,
    pub infection_non_sepsis_log_odds_age_adult: f64,
    pub infection_non_sepsis_log_odds_age_elderly: f64,
    pub infection_non_sepsis_log_odds_immunosuppressed: f64,
    pub infection_non_sepsis_log_odds_in_hospital: f64,
    pub infection_non_sepsis_minimum_bacteria_level: f64,
    pub background_mortality_baseline_log_odds: f64,
    pub mortality_baseline_1930_multiplier: f64,
    pub mortality_baseline_2035_multiplier: f64,
    pub mortality_improvement_half_life_years: f64,
    pub log_odds_mortality_per_year_of_age: f64,
    pub log_odds_mortality_per_year_of_age_squared: f64,
    pub log_odds_mortality_immunosuppressed: f64,
    pub log_odds_mortality_hospitalized: f64,
    pub base_sepsis_death_risk_per_day: f64,
    pub sepsis_age_mortality_multiplier_infant: f64,
    pub sepsis_age_mortality_multiplier_child: f64,
    pub sepsis_age_mortality_multiplier_adult: f64,
    pub sepsis_age_mortality_multiplier_elderly: f64,
    pub sepsis_immunosuppressed_multiplier: f64,
    // Enhanced microbiome/carriage model parameters
    #[allow(dead_code)]
    pub antibiotic_disruption_decay_half_life_days: f64,
    pub microbiome_resistance_multiplier_on_acquisition: f64,
    pub infection_from_microbiome_dampening: f64,
    pub carriage_duration_log_odds_coefficient: f64,
    pub carriage_duration_max_log_odds_effect: f64,
    pub antibiotic_clearance_log_odds_per_unit_activity: f64,
    pub carrier_resistance_inheritance_probability: f64,
    #[allow(dead_code)]
    pub majority_r_memory_retention_per_day: f64,
    pub microbiome_majority_decay_half_life_days: f64,
    pub microbiome_minority_decay_half_life_days: f64,
    pub microbiome_majority_promotion_rate_per_day: f64,
    pub majority_r_window_days: u32,
    pub majority_r_min_total_samples: u32,
    pub majority_r_freeze_at_last_positive: bool,
}

impl GlobalScalars {
    fn from_map(map: &HashMap<String, f64>) -> Self {
        // Reads configuration values already present in `map`; the fallback literal only applies when no entry exists.
        let majority_r_freeze_at_last_positive =
            get_or_default(map, "majority_r_freeze_at_last_positive", 1.0) > 0.5;
        GlobalScalars {
            drug_base_initiation_rate_per_day: get_or_default(
                map,
                "drug_base_initiation_rate_per_day",
                0.0010,
            ),
            drug_infection_present_multiplier: get_or_default(
                map,
                "drug_infection_present_multiplier",
                260.0,
            ),
            drug_activity_to_bacteria_level_multiplier: get_or_default(
                map,
                "drug_activity_to_bacteria_level_multiplier",
                1.0,
            ),
            already_on_drug_initiation_multiplier: get_or_default(
                map,
                "already_on_drug_initiation_multiplier",
                1.2,
            ),
            drug_test_identified_multiplier: get_or_default(
                map,
                "drug_test_identified_multiplier",
                2.5,
            ),
            double_dose_probability_if_identified_infection: get_or_default(
                map,
                "double_dose_probability_if_identified_infection",
                0.25,
            ),
            random_drug_cessation_probability: get_or_default(
                map,
                "random_drug_cessation_probability",
                0.0045,
            ),
            random_drug_cessation_probability_if_no_active_infection: get_or_default(
                map,
                "random_drug_cessation_probability_if_no_active_infection",
                0.15,
            ),
            immunodeficiency_prophylactic_drug_multiplier: get_or_default(
                map,
                "immunodeficiency_prophylactic_drug_multiplier",
                8.0,
            ),
            microbiome_resistance_transfer_probability_per_day: get_or_default(
                map,
                "microbiome_resistance_transfer_probability_per_day",
                0.0008,
            ),
            hospital_baseline_rate_per_day: get_or_default(
                map,
                "hospitalization_baseline_rate_per_day",
                0.00003,
            ),
            hospital_age_multiplier_per_day: get_or_default(
                map,
                "hospitalization_age_multiplier_per_day",
                0.00000005,
            ),
            hospital_recovery_rate_per_day: get_or_default(
                map,
                "hospitalization_recovery_rate_per_day",
                0.25,
            ),
            hospital_max_days: get_or_default(map, "hospitalization_max_days", 30.0),
            hospital_sepsis_admission_multiplier: get_or_default(
                map,
                "hospitalization_sepsis_admission_multiplier",
                80.0,
            ),
            hospital_prevent_discharge_with_sepsis: get_or_default(
                map,
                "hospitalization_prevent_discharge_with_sepsis",
                0.0,
            ),
            travel_probability_per_day: get_or_default(map, "travel_probability_per_day", 0.0005),
            antibiotic_infection_prevention_efficacy: get_or_default(
                map,
                "antibiotic_infection_prevention_efficacy",
                0.7,
            ),
            max_resistance_level: get_or_default(map, "max_resistance_level", 1.0),
            resistance_emergence_bacteria_level_multiplier: get_or_default(
                map,
                "resistance_emergence_bacteria_level_multiplier",
                0.08,
            ),
            resistance_emergence_pop_size_multiplier: get_or_default(
                map,
                "resistance_emergence_pop_size_multiplier",
                30.0,
            ),
            any_r_emergence_level_on_first_emergence: get_or_default(
                map,
                "any_r_emergence_level_on_first_emergence",
                0.5,
            ),
            multi_drug_penalty_threshold_num_drugs: get_or_default(
                map,
                "multi_drug_penalty_threshold_num_drugs",
                2.0,
            ),
            resistance_development_inhibition_single_drug: get_or_default(
                map,
                "resistance_development_inhibition_single_drug",
                0.05,
            ),
            resistance_development_inhibition_partial_cross: get_or_default(
                map,
                "resistance_development_inhibition_partial_cross",
                0.3,
            ),
            mechanism_assignment_probability_on_any_r_gain: get_or_default(
                map,
                "mechanism_assignment_probability_on_any_r_gain",
                0.8,
            ),
            treatment_failure_enabled: get_or_default(
                map,
                "enable_treatment_failure_assessment",
                1.0,
            ) > 0.5,
            treatment_failure_assessment_day: get_or_default(
                map,
                "treatment_failure_assessment_day",
                4.0,
            ) as i32,
            treatment_failure_threshold: get_or_default(map, "treatment_failure_threshold", 0.5),
            drug_failure_memory_days: get_or_default(map, "drug_failure_memory_days", 14.0) as i32,
            minimal_potency_threshold_for_drug_selection: get_or_default(
                map,
                "minimal_potency_threshold_for_drug_selection",
                0.15,
            ),
            drug_selection_temperature: get_or_default(map, "drug_selection_temperature", 1.1),
            reserve_drug_score_penalty: get_or_default(
                map,
                "reserve_drug_score_penalty",
                0.00005,
            ),
            restart_window_enabled: get_or_default(map, "enable_restart_window", 1.0) > 0.5,
            restart_window_days: get_or_default(map, "restart_window_days", 5.0) as i32,
            restart_bacteria_level_threshold: get_or_default(
                map,
                "restart_bacteria_level_threshold",
                1.5,
            ),
            restart_window_probability: get_or_default(map, "restart_window_probability", 0.3),
            log_odds_sepsis_region_a: get_or_default(map, "log_odds_sepsis_region_a", -0.3),
            log_odds_sepsis_region_b: get_or_default(map, "log_odds_sepsis_region_b", 0.2),
            mdr_tb_pre_antibiotic_era_multiplier: get_or_default(
                map,
                "mdr_tb_pre_antibiotic_era_multiplier",
                0.0001,
            ),
            mdr_tb_early_antibiotic_era_multiplier: get_or_default(
                map,
                "mdr_tb_early_antibiotic_era_multiplier",
                0.01,
            ),
            mdr_tb_modern_era_multiplier: get_or_default(map, "mdr_tb_modern_era_multiplier", 1.0),
            microbiome_resistance_emergence_rate_per_day_baseline: get_or_default(
                map,
                "microbiome_resistance_emergence_rate_per_day_baseline",
                0.0001,
            ),
            default_toxicity_reservoir_half_life_days: get_or_default(
                map,
                "default_toxicity_reservoir_half_life_days",
                1.5,
            )
            .max(0.0),
            toxicity_age_multiplier_infant: get_or_default(
                map,
                "toxicity_age_multiplier_infant",
                1.5,
            ),
            toxicity_age_multiplier_child: get_or_default(
                map,
                "toxicity_age_multiplier_child",
                1.1,
            ),
            toxicity_age_multiplier_adult: get_or_default(
                map,
                "toxicity_age_multiplier_adult",
                1.0,
            ),
            toxicity_age_multiplier_elderly: get_or_default(
                map,
                "toxicity_age_multiplier_elderly",
                2.0,
            ),
            toxicity_immunosuppressed_multiplier: get_or_default(
                map,
                "toxicity_immunosuppressed_multiplier",
                2.5,
            ),
            toxicity_hospital_multiplier: get_or_default(map, "toxicity_hospital_multiplier", 1.2),
            regional_resistance_threshold_very_high: get_or_default(
                map,
                "regional_resistance_threshold_very_high",
                0.5,
            ),
            regional_resistance_threshold_high: get_or_default(
                map,
                "regional_resistance_threshold_high",
                0.3,
            ),
            regional_resistance_threshold_moderate: get_or_default(
                map,
                "regional_resistance_threshold_moderate",
                0.1,
            ),
            regional_resistance_penalty_very_high: get_or_default(
                map,
                "regional_resistance_penalty_very_high",
                0.2,
            ),
            regional_resistance_penalty_high: get_or_default(
                map,
                "regional_resistance_penalty_high",
                0.4,
            ),
            regional_resistance_penalty_moderate: get_or_default(
                map,
                "regional_resistance_penalty_moderate",
                0.7,
            ),
            targeted_therapy_narrow_spectrum_bonus: get_or_default(
                map,
                "targeted_therapy_narrow_spectrum_bonus",
                3.0,
            ),
            targeted_therapy_broad_spectrum_penalty: get_or_default(
                map,
                "targeted_therapy_broad_spectrum_penalty",
                0.10,
            ),
            targeted_therapy_ineffective_drug_penalty: get_or_default(
                map,
                "targeted_therapy_ineffective_drug_penalty",
                0.1,
            ),
            effective_potency_threshold_for_targeted_therapy: get_or_default(
                map,
                "effective_potency_threshold_for_targeted_therapy",
                0.10,
            ),
            empiric_therapy_broad_spectrum_bonus: get_or_default(
                map,
                "empiric_therapy_broad_spectrum_bonus",
                0.85,
            ),
            empiric_therapy_ineffective_penalty: get_or_default(
                map,
                "empiric_therapy_ineffective_drug_penalty",
                0.05,
            ),
            any_r_increase_rate_per_day_when_drug_present: get_or_default(
                map,
                "any_r_increase_rate_per_day_when_drug_present",
                0.045,
            ),
            sepsis_minimum_duration_days: get_or_default(map, "sepsis_minimum_duration_days", 1.0)
                as i32,
            sepsis_base_log_odds_of_recovery_per_day: get_required(
                map,
                "sepsis_base_log_odds_of_recovery_per_day",
            ),
            sepsis_log_odds_bacteria_level: get_required(map, "sepsis_log_odds_bacteria_level"),
            sepsis_log_odds_in_hospital: get_required(map, "sepsis_log_odds_in_hospital"),
            sepsis_log_odds_age_infant: get_required(map, "sepsis_log_odds_age_infant"),
            sepsis_log_odds_age_child: get_required(map, "sepsis_log_odds_age_child"),
            sepsis_log_odds_age_adult: get_required(map, "sepsis_log_odds_age_adult"),
            sepsis_log_odds_age_elderly: get_required(map, "sepsis_log_odds_age_elderly"),
            sepsis_log_odds_immunosuppressed: get_required(map, "sepsis_log_odds_immunosuppressed"),
            infection_non_sepsis_base_log_odds: get_or_default(
                map,
                "infection_non_sepsis_base_log_odds",
                -9.0,
            ),
            infection_non_sepsis_log_odds_per_level: get_or_default(
                map,
                "infection_non_sepsis_log_odds_per_level",
                0.0,
            ),
            infection_non_sepsis_log_odds_age_infant: get_or_default(
                map,
                "infection_non_sepsis_log_odds_age_infant",
                0.0,
            ),
            infection_non_sepsis_log_odds_age_child: get_or_default(
                map,
                "infection_non_sepsis_log_odds_age_child",
                0.0,
            ),
            infection_non_sepsis_log_odds_age_adult: get_or_default(
                map,
                "infection_non_sepsis_log_odds_age_adult",
                0.0,
            ),
            infection_non_sepsis_log_odds_age_elderly: get_or_default(
                map,
                "infection_non_sepsis_log_odds_age_elderly",
                0.0,
            ),
            infection_non_sepsis_log_odds_immunosuppressed: get_or_default(
                map,
                "infection_non_sepsis_log_odds_immunosuppressed",
                0.0,
            ),
            infection_non_sepsis_log_odds_in_hospital: get_or_default(
                map,
                "infection_non_sepsis_log_odds_in_hospital",
                0.0,
            ),
            infection_non_sepsis_minimum_bacteria_level: get_or_default(
                map,
                "infection_non_sepsis_minimum_bacteria_level",
                0.5,
            ),
            background_mortality_baseline_log_odds: get_required(
                map,
                "background_mortality_baseline_log_odds",
            ),
            mortality_baseline_1930_multiplier: get_or_default(
                map,
                "mortality_baseline_1930_multiplier",
                3.0,
            ),
            mortality_baseline_2035_multiplier: get_or_default(
                map,
                "mortality_baseline_2035_multiplier",
                1.0,
            ),
            mortality_improvement_half_life_years: get_or_default(
                map,
                "mortality_improvement_half_life_years",
                35.0,
            ),
            log_odds_mortality_per_year_of_age: get_required(
                map,
                "log_odds_mortality_per_year_of_age",
            ),
            log_odds_mortality_per_year_of_age_squared: get_or_default(
                map,
                "log_odds_mortality_per_year_of_age_squared",
                0.0,
            ),
            log_odds_mortality_immunosuppressed: get_or_default(
                map,
                "log_odds_mortality_immunosuppressed",
                0.0,
            ),
            log_odds_mortality_hospitalized: get_or_default(
                map,
                "log_odds_mortality_hospitalized",
                0.0,
            ),
            base_sepsis_death_risk_per_day: get_required(map, "base_sepsis_death_risk_per_day"),
            sepsis_age_mortality_multiplier_infant: get_or_default(
                map,
                "sepsis_age_mortality_multiplier_infant",
                3.0,
            ),
            sepsis_age_mortality_multiplier_child: get_or_default(
                map,
                "sepsis_age_mortality_multiplier_child",
                0.5,
            ),
            sepsis_age_mortality_multiplier_adult: get_or_default(
                map,
                "sepsis_age_mortality_multiplier_adult",
                1.0,
            ),
            sepsis_age_mortality_multiplier_elderly: get_or_default(
                map,
                "sepsis_age_mortality_multiplier_elderly",
                2.5,
            ),
            sepsis_immunosuppressed_multiplier: get_or_default(
                map,
                "sepsis_immunosuppressed_multiplier",
                3.0,
            ),
            antibiotic_disruption_decay_half_life_days: get_or_default(
                map,
                "antibiotic_disruption_decay_half_life_days",
                30.0,
            ),
            microbiome_resistance_multiplier_on_acquisition: get_or_default(
                map,
                "microbiome_resistance_multiplier_on_acquisition",
                0.18,
            ),
            infection_from_microbiome_dampening: get_or_default(
                map,
                "infection_from_microbiome_dampening",
                0.85,
            ),
            carriage_duration_log_odds_coefficient: get_or_default(
                map,
                "carriage_duration_log_odds_coefficient",
                -0.01,
            ),
            carriage_duration_max_log_odds_effect: get_or_default(
                map,
                "carriage_duration_max_log_odds_effect",
                -2.0,
            ),
            antibiotic_clearance_log_odds_per_unit_activity: get_or_default(
                map,
                "antibiotic_clearance_log_odds_per_unit_activity",
                0.9,
            ),
            carrier_resistance_inheritance_probability: get_or_default(
                map,
                "carrier_resistance_inheritance_probability",
                0.32,
            ),
            majority_r_memory_retention_per_day: get_or_default(
                map,
                "majority_r_memory_retention_per_day",
                0.93,
            ),
            microbiome_majority_decay_half_life_days: get_or_default(
                map,
                "microbiome_majority_decay_half_life_days",
                45.0,
            ),
            microbiome_minority_decay_half_life_days: get_or_default(
                map,
                "microbiome_minority_decay_half_life_days",
                18.0,
            ),
            microbiome_majority_promotion_rate_per_day: get_or_default(
                map,
                "microbiome_majority_promotion_rate_per_day",
                0.02,
            ),
            majority_r_window_days: {
                let fallback = get_or_default(map, "majority_r_tier1_window_days", 1000.0);
                get_or_default(map, "majority_r_window_days", fallback)
                    .max(0.0)
                    .round() as u32
            },
            majority_r_min_total_samples: {
                let fallback = get_or_default(map, "majority_r_tier1_min_total_samples", 10.0);
                get_or_default(map, "majority_r_min_total_samples", fallback)
                    .max(0.0)
                    .round() as u32
            },
            majority_r_freeze_at_last_positive,
        }
    }
}

// ---------------- 4) Immunodeficiency / region / syndrome / sex parameters ----------------
#[derive(Debug)]
pub struct ImmunodeficiencyParameters {
    temporary_onset_rate_per_day: f64,
    temporary_recovery_rate_per_day: f64,
    chronic_onset_rate_per_day: f64,
    chronic_recovery_rate_per_day: f64,
    chronic_probability_age_bands: [f64; 4],
}

#[derive(Debug)]
pub struct RegionParameters {
    travel_multiplier: [f64; RegionParameters::REGION_COUNT],
    cessation_multiplier: [f64; RegionParameters::REGION_COUNT],
    mortality_log_odds: [f64; RegionParameters::REGION_COUNT],
    sepsis_log_odds: [f64; RegionParameters::REGION_COUNT],
    sepsis_mortality_multiplier: [f64; RegionParameters::REGION_COUNT],
    testing_multiplier: [f64; RegionParameters::REGION_COUNT],
}

impl RegionParameters {
    pub const REGION_COUNT: usize = REGION_VARIANTS.len();

    fn from_map(map: &HashMap<String, f64>) -> Self {
        let mut travel_multiplier = [1.0; RegionParameters::REGION_COUNT];
        let mut cessation_multiplier = [1.0; RegionParameters::REGION_COUNT];
        let mut mortality_log_odds = [0.0; RegionParameters::REGION_COUNT];
        let mut sepsis_log_odds = [0.0; RegionParameters::REGION_COUNT];
        let mut sepsis_mortality_multiplier = [1.0; RegionParameters::REGION_COUNT];
        let mut testing_multiplier = [1.0; RegionParameters::REGION_COUNT];

        for (idx, region) in REGION_VARIANTS.iter().enumerate() {
            let key_prefix = region.to_string();
            travel_multiplier[idx] =
                get_or_default(map, &format!("{}_travel_multiplier", key_prefix), 1.0);
            cessation_multiplier[idx] =
                get_or_default(map, &format!("{}_cessation_multiplier", key_prefix), 1.0);
            mortality_log_odds[idx] = get_or_default(
                map,
                &format!("log_odds_mortality_region_{}", key_prefix),
                0.0,
            );
            sepsis_log_odds[idx] =
                get_or_default(map, &format!("sepsis_log_odds_region_{}", key_prefix), 0.0);
            sepsis_mortality_multiplier[idx] = get_or_default(
                map,
                &format!("{}_sepsis_mortality_multiplier", key_prefix),
                1.0,
            );
            testing_multiplier[idx] =
                get_or_default(map, &format!("{}_testing_multiplier", key_prefix), 1.0);
        }

        RegionParameters {
            travel_multiplier,
            cessation_multiplier,
            mortality_log_odds,
            sepsis_log_odds,
            sepsis_mortality_multiplier,
            testing_multiplier,
        }
    }

    #[inline]
    pub fn travel_multiplier(&self, region: Region) -> f64 {
        self.travel_multiplier[Self::region_index(region)]
    }

    #[inline]
    pub fn cessation_multiplier(&self, region: Region) -> f64 {
        self.cessation_multiplier[Self::region_index(region)]
    }

    #[inline]
    pub fn mortality_log_odds(&self, region: Region) -> f64 {
        self.mortality_log_odds[Self::region_index(region)]
    }

    #[inline]
    pub fn sepsis_log_odds(&self, region: Region) -> f64 {
        self.sepsis_log_odds[Self::region_index(region)]
    }

    #[inline]
    pub fn sepsis_mortality_multiplier(&self, region: Region) -> f64 {
        self.sepsis_mortality_multiplier[Self::region_index(region)]
    }

    #[inline]
    pub fn testing_multiplier(&self, region: Region) -> f64 {
        self.testing_multiplier[Self::region_index(region)]
    }

    #[inline]
    pub fn region_index(region: Region) -> usize {
        match region {
            Region::NorthAmerica => 0,
            Region::SouthAmerica => 1,
            Region::Africa => 2,
            Region::Asia => 3,
            Region::Europe => 4,
            Region::Oceania => 5,
            Region::Home => 6,
        }
    }
}

#[derive(Debug)]
pub struct SexParameters {
    male_log_odds: f64,
    female_log_odds: f64,
}

impl SexParameters {
    fn from_map(map: &HashMap<String, f64>) -> Self {
        SexParameters {
            male_log_odds: get_or_default(map, "log_odds_mortality_sex_male", 0.0),
            female_log_odds: get_or_default(map, "log_odds_mortality_sex_female", 0.0),
        }
    }

    #[inline]
    pub fn mortality_log_odds(&self, sex_at_birth: &str) -> f64 {
        if sex_at_birth.eq_ignore_ascii_case("male") {
            self.male_log_odds
        } else if sex_at_birth.eq_ignore_ascii_case("female") {
            self.female_log_odds
        } else {
            0.0
        }
    }
}

// ---------------- 5) Vaccination & age category tables ----------------
#[derive(Debug)]
pub struct VaccinationParameters {
    daily_probabilities: Vec<f64>,
    availability_years: Vec<f64>,
}

impl VaccinationParameters {
    const VACCINES: &'static [&'static str] = &["pneumococcal", "meningococcal", "hib"];

    fn index(vaccine_idx: usize, age_idx: usize) -> usize {
        vaccine_idx * AGE_BUCKETS.len() + age_idx
    }

    pub fn from_map(map: &HashMap<String, f64>) -> Self {
        let mut daily_probabilities = vec![0.0; Self::VACCINES.len() * AGE_BUCKETS.len()];
        let mut availability_years = vec![2100.0; Self::VACCINES.len()];

        for (vaccine_idx, &vaccine) in Self::VACCINES.iter().enumerate() {
            let availability_key = format!("vaccine_{}_availability_year", vaccine);
            availability_years[vaccine_idx] = get_or_default(map, &availability_key, 2100.0);

            for (age_idx, &age_group) in AGE_BUCKETS.iter().enumerate() {
                let key = format!(
                    "vaccine_{}_daily_prob_age_{}",
                    vaccine,
                    age_group.bucket_slug()
                );
                daily_probabilities[Self::index(vaccine_idx, age_idx)] =
                    get_or_default(map, &key, 0.0);
            }
        }

        VaccinationParameters {
            daily_probabilities,
            availability_years,
        }
    }

    #[inline]
    pub fn vaccine_index(vaccine: &str) -> Option<usize> {
        match vaccine {
            "pneumococcal" => Some(0),
            "meningococcal" => Some(1),
            "hib" => Some(2),
            _ => None,
        }
    }

    #[inline]
    pub fn age_group_index(age_years: f64) -> usize {
        if age_years < 1.0 {
            0
        } else if age_years < 5.0 {
            1
        } else if age_years < 18.0 {
            2
        } else if age_years < 50.0 {
            3
        } else if age_years < 70.0 {
            4
        } else {
            5
        }
    }

    #[inline]
    pub fn daily_probability(&self, vaccine_idx: usize, age_idx: usize) -> f64 {
        self.daily_probabilities[Self::index(vaccine_idx, age_idx)]
    }

    #[inline]
    pub fn availability_year(&self, vaccine_idx: usize) -> f64 {
        self.availability_years
            .get(vaccine_idx)
            .copied()
            .unwrap_or(2100.0)
    }
}

#[derive(Debug)]
pub struct SyndromeParameters {
    sepsis_log_odds: Vec<f64>,
    initiation_multiplier: Vec<f64>,
    non_sepsis_mortality_log_odds: Vec<f64>,
    empiric_drug_scores: Vec<Vec<f64>>,
}

impl SyndromeParameters {
    const MAX_SYNDROME_ID: usize = 10;

    fn from_map(map: &HashMap<String, f64>) -> Self {
        let len = Self::MAX_SYNDROME_ID + 1;
        let mut sepsis_log_odds = vec![0.0; len];
        let mut initiation_multiplier = vec![1.0; len];
        let mut non_sepsis_mortality_log_odds = vec![0.0; len];
        let mut empiric_drug_scores =
            vec![vec![1.0; DRUG_SHORT_NAMES.len()]; len];

        for syndrome_id in 1..=Self::MAX_SYNDROME_ID {
            sepsis_log_odds[syndrome_id] = get_or_default(
                map,
                &format!("log_odds_syndrome_{}_sepsis", syndrome_id),
                0.0,
            );
            initiation_multiplier[syndrome_id] = get_or_default(
                map,
                &format!("syndrome_{}_initiation_multiplier", syndrome_id),
                1.0,
            );
            non_sepsis_mortality_log_odds[syndrome_id] = get_or_default(
                map,
                &format!(
                    "syndrome_{}_non_sepsis_infection_death_log_odds",
                    syndrome_id
                ),
                0.0,
            );

            for (drug_idx, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
                let key = format!(
                    "syndrome_{}_empiric_drug_{}_score",
                    syndrome_id,
                    drug
                );
                empiric_drug_scores[syndrome_id][drug_idx] =
                    get_or_default(map, &key, 1.0);
            }
        }

        SyndromeParameters {
            sepsis_log_odds,
            initiation_multiplier,
            non_sepsis_mortality_log_odds,
            empiric_drug_scores,
        }
    }

    #[inline]
    pub fn sepsis_log_odds(&self, syndrome_id: usize) -> f64 {
        self.sepsis_log_odds
            .get(syndrome_id)
            .copied()
            .unwrap_or(0.0)
    }

    #[inline]
    pub fn initiation_multiplier(&self, syndrome_id: usize) -> f64 {
        self.initiation_multiplier
            .get(syndrome_id)
            .copied()
            .unwrap_or(1.0)
    }

    #[inline]
    pub fn non_sepsis_mortality_log_odds(&self, syndrome_id: usize) -> f64 {
        self.non_sepsis_mortality_log_odds
            .get(syndrome_id)
            .copied()
            .unwrap_or(0.0)
    }

    #[inline]
    pub fn empiric_drug_score(&self, syndrome_id: usize, drug_idx: usize) -> f64 {
        self.empiric_drug_scores
            .get(syndrome_id)
            .and_then(|scores| scores.get(drug_idx))
            .copied()
            .unwrap_or(1.0)
    }
}

impl ImmunodeficiencyParameters {
    fn from_map(map: &HashMap<String, f64>) -> Self {
        let temporary_onset_rate_per_day = get_or_default(
            map,
            "temporary_immunosuppression_onset_rate_per_day",
            0.00005,
        );
        let temporary_recovery_rate_per_day = get_or_default(
            map,
            "temporary_immunosuppression_recovery_rate_per_day",
            0.01,
        );
        let chronic_onset_rate_per_day =
            get_or_default(map, "chronic_immunosuppression_onset_rate_per_day", 0.00006);
        let chronic_recovery_rate_per_day = get_or_default(
            map,
            "chronic_immunosuppression_recovery_rate_per_day",
            0.0012,
        );

        let chronic_probability_age_bands = [
            get_or_default(map, "chronic_immunodeficiency_probability_age_0_1", 0.3),
            get_or_default(map, "chronic_immunodeficiency_probability_age_1_18", 0.2),
            get_or_default(map, "chronic_immunodeficiency_probability_age_18_65", 0.4),
            get_or_default(map, "chronic_immunodeficiency_probability_age_65_plus", 0.6),
        ];

        ImmunodeficiencyParameters {
            temporary_onset_rate_per_day,
            temporary_recovery_rate_per_day,
            chronic_onset_rate_per_day,
            chronic_recovery_rate_per_day,
            chronic_probability_age_bands,
        }
    }

    #[inline]
    pub fn temporary_onset_rate(&self) -> f64 {
        self.temporary_onset_rate_per_day
    }

    #[inline]
    pub fn temporary_recovery_rate(&self) -> f64 {
        self.temporary_recovery_rate_per_day
    }

    #[inline]
    pub fn chronic_onset_rate(&self) -> f64 {
        self.chronic_onset_rate_per_day
    }

    #[inline]
    pub fn chronic_recovery_rate(&self) -> f64 {
        self.chronic_recovery_rate_per_day
    }

    #[inline]
    pub fn chronic_probability(&self, age_days: i32) -> f64 {
        let age_days = age_days.max(0);
        if age_days <= 365 {
            self.chronic_probability_age_bands[0]
        } else if age_days <= 6570 {
            self.chronic_probability_age_bands[1]
        } else if age_days <= 23725 {
            self.chronic_probability_age_bands[2]
        } else {
            self.chronic_probability_age_bands[3]
        }
    }
}

// ---------------- 6) Drug parameter blocks (initiation, failure, restart) ----------------
#[derive(Debug)]
pub struct DrugParameters {
    pub initial_level: Vec<f64>,
    pub double_dose_multiplier: Vec<f64>,
    pub spectrum_breadth: Vec<f64>,
    pub half_life_days: Vec<f64>,
    pub toxicity_death_hazard_per_unit_level: Vec<f64>,
    pub toxicity_reservoir_half_life_days: Vec<f64>,
    pub microbiome_disruption_log_odds: Vec<f64>,
}

impl DrugParameters {
    fn from_map(map: &HashMap<String, f64>, num_drugs: usize) -> Self {
        let mut initial_level = Vec::with_capacity(num_drugs);
        let mut double_dose_multiplier = Vec::with_capacity(num_drugs);
        let mut spectrum_breadth = Vec::with_capacity(num_drugs);
        let mut half_life_days = Vec::with_capacity(num_drugs);
        let mut toxicity_death_hazard_per_unit_level = Vec::with_capacity(num_drugs);
        let mut toxicity_reservoir_half_life_days = Vec::with_capacity(num_drugs);
        let mut microbiome_disruption_log_odds = Vec::with_capacity(num_drugs);

        for &drug in DRUG_SHORT_NAMES.iter() {
            let prefix = format!("drug_{}", drug);
            initial_level.push(get_or_default(
                map,
                &format!("{}_initial_level", prefix),
                10.0,
            ));
            double_dose_multiplier.push(get_or_default(
                map,
                &format!("{}_double_dose_multiplier", prefix),
                2.0,
            ));
            spectrum_breadth.push(get_or_default(
                map,
                &format!("{}_spectrum_breadth", prefix),
                3.0,
            ));
            half_life_days.push(get_or_default(
                map,
                &format!("{}_half_life_days", prefix),
                0.25,
            ));
            let hazard = get_or_default(
                map,
                &format!("{}_toxicity_death_hazard_per_unit_level", prefix),
                get_or_default(
                    map,
                    "default_drug_toxicity_death_hazard_per_unit_level",
                    0.0,
                ),
            );
            toxicity_death_hazard_per_unit_level.push(hazard.max(0.0));
            let half_life = get_or_default(
                map,
                &format!("{}_toxicity_reservoir_half_life_days", prefix),
                get_or_default(map, "default_toxicity_reservoir_half_life_days", 1.5),
            );
            toxicity_reservoir_half_life_days.push(half_life.max(0.0));
            microbiome_disruption_log_odds.push(get_or_default(
                map,
                &format!("{}_microbiome_disruption_log_odds", prefix),
                get_or_default(map, "default_microbiome_disruption_log_odds", 0.0),
            ));
        }

        DrugParameters {
            initial_level,
            double_dose_multiplier,
            spectrum_breadth,
            half_life_days,
            toxicity_death_hazard_per_unit_level,
            toxicity_reservoir_half_life_days,
            microbiome_disruption_log_odds,
        }
    }

    #[inline]
    pub fn initial_level(&self, drug_idx: usize) -> f64 {
        self.initial_level[drug_idx]
    }

    #[inline]
    pub fn double_dose_multiplier(&self, drug_idx: usize) -> f64 {
        self.double_dose_multiplier[drug_idx]
    }

    #[inline]
    pub fn spectrum_breadth(&self, drug_idx: usize) -> f64 {
        self.spectrum_breadth[drug_idx]
    }

    #[inline]
    pub fn half_life_days(&self, drug_idx: usize) -> f64 {
        self.half_life_days[drug_idx]
    }

    #[inline]
    pub fn toxicity_death_hazard_per_unit_level(&self, drug_idx: usize) -> f64 {
        self.toxicity_death_hazard_per_unit_level[drug_idx]
    }

    #[inline]
    pub fn toxicity_reservoir_half_life_days(&self, drug_idx: usize) -> f64 {
        self.toxicity_reservoir_half_life_days[drug_idx]
    }

    #[inline]
    pub fn microbiome_disruption_log_odds(&self, drug_idx: usize) -> f64 {
        self.microbiome_disruption_log_odds[drug_idx]
    }
}

// ---------------- 7) Bacteria-level parameters & microbiome logic ----------------
#[derive(Debug)]
#[allow(dead_code)]
pub struct BacteriaParameters {
    pub acquisition_log_odds_baseline: Vec<f64>,
    pub log_odds_vaccinated: Vec<f64>,
    pub log_odds_microbiome_present: Vec<f64>,
    pub log_odds_hospital_acquired: Vec<f64>,
    pub microbiome_clearance_probability_per_day: Vec<f64>,
    pub environmental_acquisition_proportion: Vec<f64>,
    pub initial_infection_level: Vec<f64>,
    pub base_bacteria_level_change: Vec<f64>,
    pub max_level: Vec<f64>,
    pub daily_symptom_onset_probability: Vec<f64>,
    pub symptom_onset_threshold_level: Vec<f64>,
    pub symptom_onset_delay_days: Vec<f64>,
    pub symptom_onset_level_multiplier: Vec<f64>,
    pub microbiome_vs_infection_log_odds: Vec<f64>,
    pub drug_cessation_probability: Vec<f64>,
    pub treatment_recognition_year: Vec<Option<f64>>,
    pub sepsis_baseline_log_odds: Vec<f64>,
    pub sepsis_log_odds_infection_level: Vec<f64>,
    pub sepsis_log_odds_infection_duration: Vec<f64>,
    pub infection_non_sepsis_mortality_log_odds: Vec<f64>,
}

impl BacteriaParameters {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize) -> Self {
        let mut acquisition_log_odds_baseline = Vec::with_capacity(num_bacteria);
        let mut log_odds_vaccinated = Vec::with_capacity(num_bacteria);
        let mut log_odds_microbiome_present = Vec::with_capacity(num_bacteria);
        let mut log_odds_hospital_acquired = Vec::with_capacity(num_bacteria);
        let mut microbiome_clearance_probability_per_day = Vec::with_capacity(num_bacteria);
        let mut environmental_acquisition_proportion = Vec::with_capacity(num_bacteria);
        let mut initial_infection_level = Vec::with_capacity(num_bacteria);
        let mut base_bacteria_level_change = Vec::with_capacity(num_bacteria);
        let mut max_level = Vec::with_capacity(num_bacteria);
        let mut daily_symptom_onset_probability = Vec::with_capacity(num_bacteria);
        let mut symptom_onset_threshold_level = Vec::with_capacity(num_bacteria);
        let mut symptom_onset_delay_days = Vec::with_capacity(num_bacteria);
        let mut symptom_onset_level_multiplier = Vec::with_capacity(num_bacteria);
        let mut microbiome_vs_infection_log_odds = Vec::with_capacity(num_bacteria);
        let mut drug_cessation_probability = Vec::with_capacity(num_bacteria);
        let mut treatment_recognition_year = Vec::with_capacity(num_bacteria);
        let mut sepsis_baseline_log_odds = Vec::with_capacity(num_bacteria);
        let mut sepsis_log_odds_infection_level = Vec::with_capacity(num_bacteria);
        let mut sepsis_log_odds_infection_duration = Vec::with_capacity(num_bacteria);
        let mut infection_non_sepsis_mortality_log_odds = Vec::with_capacity(num_bacteria);

        for &bacteria in BACTERIA_LIST.iter() {
            let prefix = bacteria;
            acquisition_log_odds_baseline.push(get_or_default(
                map,
                &format!("{}_acquisition_log_odds_baseline", prefix),
                get_or_default(map, "acquisition_log_odds_baseline", -30.0),
            ));
            log_odds_vaccinated.push(get_or_default(
                map,
                &format!("{}_log_odds_vaccinated", prefix),
                get_or_default(map, "log_odds_vaccinated", 0.0),
            ));
            log_odds_microbiome_present.push(get_or_default(
                map,
                &format!("{}_log_odds_microbiome_present", prefix),
                get_or_default(map, "log_odds_microbiome_present", 0.0),
            ));
            log_odds_hospital_acquired.push(get_or_default(
                map,
                &format!("{}_log_odds_hospital_acquired", prefix),
                get_or_default(map, "log_odds_hospital_acquired", 0.0),
            ));
            microbiome_clearance_probability_per_day.push(get_or_default(
                map,
                &format!("{}_microbiome_clearance_probability_per_day", prefix),
                get_or_default(
                    map,
                    "default_microbiome_clearance_probability_per_day",
                    0.075,
                ),
            ));
            environmental_acquisition_proportion.push(get_or_default(
                map,
                &format!("{}_environmental_acquisition_proportion", prefix),
                0.1,
            ));
            initial_infection_level.push(get_or_default(
                map,
                &format!("{}_initial_infection_level", prefix),
                0.01,
            ));
            base_bacteria_level_change.push(get_or_default(
                map,
                &format!("{}_base_bacteria_level_change", prefix),
                0.5,
            ));
            max_level.push(get_or_default(map, &format!("{}_max_level", prefix), 5.0));
            daily_symptom_onset_probability.push(get_or_default(
                map,
                &format!("{}_daily_symptom_onset_probability", prefix),
                0.15,
            ));
            symptom_onset_threshold_level.push(get_or_default(
                map,
                &format!("{}_symptom_onset_threshold_level", prefix),
                0.5,
            ));
            symptom_onset_delay_days.push(get_or_default(
                map,
                &format!("{}_symptom_onset_delay_days", prefix),
                1.0,
            ));
            symptom_onset_level_multiplier.push(get_or_default(
                map,
                &format!("{}_symptom_onset_level_multiplier", prefix),
                1.0,
            ));
            microbiome_vs_infection_log_odds.push(get_or_default(
                map,
                &format!(
                    "{}_log_odds_microbiome_vs_infection",
                    prefix
                ),
                get_or_default(map, "log_odds_microbiome_vs_infection", 3.0), // Iter7: was -6.0 (99.75% infection!), changed to 3.0 (~5% infection)
            ));
            let cessation_key = format!(
                "{}_drug_cessation_probability",
                prefix.to_lowercase()
            );
            let default_cessation = get_or_default(map, "random_drug_cessation_probability", 0.001);
            drug_cessation_probability.push(get_or_default(map, &cessation_key, default_cessation));
            let recognition_key =
                format!("{}_treatment_recognition_year", prefix);
            treatment_recognition_year.push(map.get(&recognition_key).copied());
            sepsis_baseline_log_odds.push(get_or_default(
                map,
                &format!("{}_sepsis_baseline_log_odds", prefix),
                get_or_default(map, "sepsis_baseline_log_odds", -5.0),
            ));
            sepsis_log_odds_infection_level.push(get_or_default(
                map,
                &format!("{}_log_odds_sepsis_infection_level", prefix),
                get_or_default(map, "log_odds_sepsis_infection_level", 0.1),
            ));
            sepsis_log_odds_infection_duration.push(get_or_default(
                map,
                &format!("{}_log_odds_sepsis_infection_duration", prefix),
                get_or_default(map, "log_odds_sepsis_infection_duration", 0.01),
            ));
            infection_non_sepsis_mortality_log_odds.push(get_or_default(
                map,
                &format!("{}_non_sepsis_infection_death_log_odds", prefix),
                0.0,
            ));
        }

        BacteriaParameters {
            acquisition_log_odds_baseline,
            log_odds_vaccinated,
            log_odds_microbiome_present,
            log_odds_hospital_acquired,
            microbiome_clearance_probability_per_day,
            environmental_acquisition_proportion,
            initial_infection_level,
            base_bacteria_level_change,
            max_level,
            daily_symptom_onset_probability,
            symptom_onset_threshold_level,
            symptom_onset_delay_days,
            symptom_onset_level_multiplier,
            microbiome_vs_infection_log_odds,
            drug_cessation_probability,
            treatment_recognition_year,
            sepsis_baseline_log_odds,
            sepsis_log_odds_infection_level,
            sepsis_log_odds_infection_duration,
            infection_non_sepsis_mortality_log_odds,
        }
    }

    #[inline]
    #[allow(dead_code)]
    pub fn base_level_change(&self, bacteria_idx: usize) -> f64 {
        self.base_bacteria_level_change[bacteria_idx]
    }

    #[inline]
    pub fn sepsis_baseline_log_odds(&self, bacteria_idx: usize) -> f64 {
        self.sepsis_baseline_log_odds[bacteria_idx]
    }

    #[inline]
    pub fn sepsis_log_odds_infection_level(&self, bacteria_idx: usize) -> f64 {
        self.sepsis_log_odds_infection_level[bacteria_idx]
    }

    #[inline]
    pub fn sepsis_log_odds_infection_duration(&self, bacteria_idx: usize) -> f64 {
        self.sepsis_log_odds_infection_duration[bacteria_idx]
    }

    #[inline]
    pub fn infection_non_sepsis_mortality_log_odds(&self, bacteria_idx: usize) -> f64 {
        self.infection_non_sepsis_mortality_log_odds[bacteria_idx]
    }

    #[inline]
    pub fn microbiome_vs_infection_log_odds(&self, bacteria_idx: usize) -> f64 {
        self.microbiome_vs_infection_log_odds[bacteria_idx]
    }

    #[inline]
    pub fn treatment_recognition_year(&self, bacteria_idx: usize) -> Option<f64> {
        self.treatment_recognition_year[bacteria_idx]
    }

    #[inline]
    pub fn microbiome_clearance_probability_per_day(&self, bacteria_idx: usize) -> f64 {
        self.microbiome_clearance_probability_per_day[bacteria_idx]
    }

    #[inline]
    pub fn environmental_acquisition_proportion(&self, bacteria_idx: usize) -> f64 {
        self.environmental_acquisition_proportion[bacteria_idx]
    }

    #[inline]
    pub fn initial_infection_level(&self, bacteria_idx: usize) -> f64 {
        self.initial_infection_level[bacteria_idx]
    }

    #[inline]
    pub fn max_level(&self, bacteria_idx: usize) -> f64 {
        self.max_level[bacteria_idx]
    }

    #[inline]
    pub fn daily_symptom_onset_probability(&self, bacteria_idx: usize) -> f64 {
        self.daily_symptom_onset_probability[bacteria_idx]
    }

    #[inline]
    pub fn symptom_onset_threshold_level(&self, bacteria_idx: usize) -> f64 {
        self.symptom_onset_threshold_level[bacteria_idx]
    }

    #[inline]
    pub fn symptom_onset_delay_days(&self, bacteria_idx: usize) -> f64 {
        self.symptom_onset_delay_days[bacteria_idx]
    }

    #[inline]
    pub fn symptom_onset_level_multiplier(&self, bacteria_idx: usize) -> f64 {
        self.symptom_onset_level_multiplier[bacteria_idx]
    }
}

// ---------------- 8) Clearance, acquisition, and age tables ----------------
#[derive(Debug)]
pub struct ClearanceParameters {
    base_delay_days: f64,
    base_daily_hazard: f64,
    per_bacteria_delay_days: Vec<Option<f64>>,
    per_bacteria_hazard_multiplier: Vec<f64>,
    age_multipliers: [f64; AGE_BUCKET_COUNT],
    immunodeficient_multiplier: f64,
    level_reference: f64,
    level_exponent: f64,
}

impl ClearanceParameters {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize) -> Self {
        let base_delay_days = get_or_default(map, "default_clearance_delay_days", 3.0);
        let base_daily_hazard = get_or_default(map, "default_clearance_hazard_after_delay", 0.045);

        let mut age_multipliers = [1.0; AGE_BUCKET_COUNT];
        for (idx, category) in AGE_BUCKETS.iter().enumerate() {
            let key = format!("clearance_age_multiplier_{}", category.label());
            age_multipliers[idx] = get_or_default(map, &key, 1.0);
        }

        let mut per_bacteria_delay_days = Vec::with_capacity(num_bacteria);
        let mut per_bacteria_hazard_multiplier = Vec::with_capacity(num_bacteria);
        for &bacteria in BACTERIA_LIST.iter() {
            per_bacteria_delay_days.push(
                map.get(&format!("{}_clearance_delay_days", bacteria))
                    .copied(),
            );
            per_bacteria_hazard_multiplier.push(get_or_default(
                map,
                &format!("{}_clearance_hazard_multiplier", bacteria),
                1.0,
            ));
        }

        let immunodeficient_multiplier =
            get_or_default(map, "clearance_immunodeficient_multiplier", 0.5);
        let level_reference = get_or_default(map, "clearance_level_reference", 1.0);
        let level_exponent = get_or_default(map, "clearance_level_exponent", 0.0);

        ClearanceParameters {
            base_delay_days,
            base_daily_hazard,
            per_bacteria_delay_days,
            per_bacteria_hazard_multiplier,
            age_multipliers,
            immunodeficient_multiplier,
            level_reference,
            level_exponent,
        }
    }

    #[inline]
    pub fn delay_days(&self, bacteria_idx: usize) -> f64 {
        self.per_bacteria_delay_days[bacteria_idx]
            .unwrap_or(self.base_delay_days)
            .max(0.0)
    }

    #[inline]
    pub fn hazard(&self, bacteria_idx: usize) -> f64 {
        (self.base_daily_hazard * self.per_bacteria_hazard_multiplier[bacteria_idx]).clamp(0.0, 1.0)
    }

    #[inline]
    pub fn age_multiplier(&self, age_days: i32) -> f64 {
        let idx =
            AgeCategoryParameters::age_category_index(age_days).min(self.age_multipliers.len() - 1);
        self.age_multipliers[idx]
    }

    #[inline]
    pub fn immunodeficient_multiplier(&self, is_immunodeficient: bool) -> f64 {
        if is_immunodeficient {
            self.immunodeficient_multiplier
        } else {
            1.0
        }
    }

    #[inline]
    pub fn level_modifier(&self, level: f64) -> f64 {
        if self.level_exponent <= 0.0 {
            return 1.0;
        }

        let ratio = self.level_reference / (self.level_reference + level.max(0.0) + f64::EPSILON);
        ratio.powf(self.level_exponent)
    }

    #[inline]
    pub fn hazard_for(
        &self,
        bacteria_idx: usize,
        age_days: i32,
        is_immunodeficient: bool,
        level: f64,
    ) -> f64 {
        let base = self.hazard(bacteria_idx);
        if base <= 0.0 {
            return 0.0;
        }

        let age_factor = self.age_multiplier(age_days).max(0.0);
        let immuno_factor = self.immunodeficient_multiplier(is_immunodeficient).max(0.0);
        let level_factor = self.level_modifier(level).max(0.0);

        (base * age_factor * immuno_factor * level_factor).clamp(0.0, 1.0)
    }
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct DrugBacteriaMatrix {
    pub potency_when_no_r: Vec<f64>,
    pub initiation_multiplier: Vec<f64>,
    pub resistance_emergence_rate: Vec<f64>,
    pub mic_lt2_threshold: Vec<f64>,
    num_bacteria: usize,
    num_drugs: usize,
}

impl DrugBacteriaMatrix {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize, num_drugs: usize) -> Self {
        let mut potency_when_no_r = Vec::with_capacity(num_bacteria * num_drugs);
        let mut initiation_multiplier = Vec::with_capacity(num_bacteria * num_drugs);
        let mut resistance_emergence_rate = Vec::with_capacity(num_bacteria * num_drugs);
        let mut mic_lt2_threshold = Vec::with_capacity(num_bacteria * num_drugs);

        for &bacteria in BACTERIA_LIST.iter() {
            for &drug in DRUG_SHORT_NAMES.iter() {
                // Bacteria names now use underscores consistently
                let key_prefix = format!("drug_{}_for_bacteria_{}", drug, bacteria);
                let potency =
                    get_or_default(map, &format!("{}_potency_when_no_r", key_prefix), 0.1).max(0.0);
                potency_when_no_r.push(potency);
                initiation_multiplier.push(get_or_default(
                    map,
                    &format!("{}_initiation_multiplier", key_prefix),
                    1.0,
                ));
                resistance_emergence_rate.push(get_or_default(
                    map,
                    &format!("{}_resistance_emergence_rate_per_day_baseline", key_prefix),
                    0.0000,
                ));
                let threshold = 2.0 * potency.min(1.0) - 1.0;
                mic_lt2_threshold.push(threshold);
            }
        }

        DrugBacteriaMatrix {
            potency_when_no_r,
            initiation_multiplier,
            resistance_emergence_rate,
            mic_lt2_threshold,
            num_bacteria,
            num_drugs,
        }
    }

    #[inline]
    fn index(&self, bacteria_idx: usize, drug_idx: usize) -> usize {
        bacteria_idx * self.num_drugs + drug_idx
    }

    #[inline]
    pub fn potency(&self, bacteria_idx: usize, drug_idx: usize) -> f64 {
        let idx = self.index(bacteria_idx, drug_idx);
        self.potency_when_no_r[idx]
    }

    #[inline]
    pub fn initiation_multiplier(&self, bacteria_idx: usize, drug_idx: usize) -> f64 {
        self.initiation_multiplier[self.index(bacteria_idx, drug_idx)]
    }

    #[inline]
    pub fn resistance_emergence_rate(&self, bacteria_idx: usize, drug_idx: usize) -> f64 {
        self.resistance_emergence_rate[self.index(bacteria_idx, drug_idx)]
    }

    #[inline]
    #[allow(dead_code)]
    pub fn mic_lt2_threshold(&self, bacteria_idx: usize, drug_idx: usize) -> f64 {
        self.mic_lt2_threshold[self.index(bacteria_idx, drug_idx)]
    }
}

#[derive(Debug)]
pub struct RegionBacteriaAcquisition {
    values: Vec<f64>,
    num_bacteria: usize,
}

impl RegionBacteriaAcquisition {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize) -> Self {
        let mut values = Vec::with_capacity(REGION_NAMES.len() * num_bacteria);

        for &region in REGION_NAMES
            .iter()
            .chain(std::iter::once(&HOME_REGION_NAME))
        {
            for &bacteria in BACTERIA_LIST.iter() {
                let key = format!("{}_{}_acquisition_log_odds", region, bacteria);
                let default_key = format!("{}_acquisition_log_odds_default", region);
                let val = map
                    .get(&key)
                    .copied()
                    .unwrap_or_else(|| get_or_default(map, &default_key, 0.0));
                values.push(val);
            }
        }

        RegionBacteriaAcquisition {
            values,
            num_bacteria,
        }
    }

    fn region_index(region: Region) -> usize {
        match region {
            Region::NorthAmerica => 0,
            Region::SouthAmerica => 1,
            Region::Africa => 2,
            Region::Asia => 3,
            Region::Europe => 4,
            Region::Oceania => 5,
            Region::Home => 6,
        }
    }

    #[inline]
    fn index(region_idx: usize, bacteria_idx: usize, num_bacteria: usize) -> usize {
        region_idx * num_bacteria + bacteria_idx
    }

    #[inline]
    pub fn acquisition_log_odds(&self, region: Region, bacteria_idx: usize) -> f64 {
        let idx = Self::index(Self::region_index(region), bacteria_idx, self.num_bacteria);
        self.values[idx]
    }
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct AgeTables {
    pub bacteria_age_log_odds: Vec<f64>,
    pub region_age_log_odds: Vec<f64>,
    pub default_age_log_odds: Vec<f64>,
    num_bacteria: usize,
}

impl AgeTables {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize) -> Self {
        let mut bacteria_age_log_odds = Vec::with_capacity(num_bacteria * AGE_BUCKETS.len());
        let mut region_age_log_odds =
            Vec::with_capacity((REGION_NAMES.len() + 1) * AGE_BUCKETS.len());
        let mut default_age_log_odds = Vec::with_capacity(AGE_BUCKETS.len());

        for &bacteria in BACTERIA_LIST.iter() {
            for &category in AGE_BUCKETS.iter() {
                let age_slug = category.bucket_slug();
                let key = format!("{}_log_odds_{}", bacteria, age_slug);
                let default_key = format!("default_log_odds_{}", age_slug);
                let value = map
                    .get(&key)
                    .copied()
                    .unwrap_or_else(|| get_or_default(map, &default_key, 0.0));
                bacteria_age_log_odds.push(value);
            }
        }

        for &region in REGION_NAMES
            .iter()
            .chain(std::iter::once(&HOME_REGION_NAME))
        {
            for &category in AGE_BUCKETS.iter() {
                let age_slug = category.bucket_slug();
                let key = format!("{}_log_odds_{}", region, age_slug);
                region_age_log_odds.push(get_or_default(map, &key, 0.0));
            }
        }

        for &category in AGE_BUCKETS.iter() {
            default_age_log_odds.push(get_or_default(
                map,
                &format!("default_log_odds_{}", category.bucket_slug()),
                0.0,
            ));
        }

        AgeTables {
            bacteria_age_log_odds,
            region_age_log_odds,
            default_age_log_odds,
            num_bacteria,
        }
    }

    #[inline]
    #[allow(dead_code)]
    pub fn bacteria_age_log_odds(&self, bacteria_idx: usize, age_bucket: usize) -> f64 {
        self.bacteria_age_log_odds[bacteria_idx * AGE_BUCKETS.len() + age_bucket]
    }

    #[inline]
    #[allow(dead_code)]
    pub fn region_age_log_odds(&self, region_idx: usize, age_bucket: usize) -> f64 {
        self.region_age_log_odds[region_idx * AGE_BUCKETS.len() + age_bucket]
    }

    #[inline]
    #[allow(dead_code)]
    pub fn default_log_odds(&self, age_bucket: usize) -> f64 {
        self.default_age_log_odds[age_bucket]
    }
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct AgeCategoryParameters {
    bacteria_age_log_odds: Vec<f64>,
    region_age_log_odds: Vec<f64>,
    bacteria_region_age_log_odds: Vec<f64>,
    default_age_log_odds: [f64; AGE_BUCKET_COUNT],
    num_bacteria: usize,
}

impl AgeCategoryParameters {
    const AGE_CATEGORY_COUNT: usize = AGE_BUCKET_COUNT;

    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize) -> Self {
        let age_count = Self::AGE_CATEGORY_COUNT;

        let mut default_age_log_odds = [0.0; AGE_BUCKET_COUNT];
        for (idx, category) in AGE_BUCKETS.iter().enumerate() {
            default_age_log_odds[idx] = get_or_default(
                map,
                &format!("default_log_odds_{}", category.label()),
                0.0,
            );
        }

        let mut bacteria_age_log_odds = Vec::with_capacity(num_bacteria * age_count);
        for &bacteria in BACTERIA_LIST.iter() {
            for (age_idx, category) in AGE_BUCKETS.iter().enumerate() {
                let key = format!("{}_log_odds_{}", bacteria, category.label());
                let value = map
                    .get(&key)
                    .copied()
                    .unwrap_or(default_age_log_odds[age_idx]);
                bacteria_age_log_odds.push(value);
            }
        }

        let mut region_age_log_odds =
            Vec::with_capacity(RegionParameters::REGION_COUNT * age_count);
        for region in REGION_VARIANTS.iter() {
            let region_key = region.to_string();
            for category in AGE_BUCKETS.iter() {
                region_age_log_odds.push(get_or_default(
                    map,
                    &format!("{}_log_odds_{}", region_key, category.label()),
                    0.0,
                ));
            }
        }

        let mut bacteria_region_age_log_odds =
            Vec::with_capacity(RegionParameters::REGION_COUNT * num_bacteria * age_count);
        for region in REGION_VARIANTS.iter() {
            let region_key = region.to_string();
            let region_idx = RegionParameters::region_index(*region);
            for (_bacteria_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
                for (age_idx, category) in AGE_BUCKETS.iter().enumerate() {
                    let fallback = region_age_log_odds[region_idx * age_count + age_idx];
                    let key = format!("{}_{}_log_odds_{}", bacteria, region_key, category.label());
                    let value = map.get(&key).copied().unwrap_or(fallback);
                    bacteria_region_age_log_odds.push(value);
                }
            }
        }

        AgeCategoryParameters {
            bacteria_age_log_odds,
            region_age_log_odds,
            bacteria_region_age_log_odds,
            default_age_log_odds,
            num_bacteria,
        }
    }

    #[inline]
    fn age_count() -> usize {
        Self::AGE_CATEGORY_COUNT
    }

    #[inline]
    pub fn age_category_index(age_days: i32) -> usize {
        match AgeCategory::from_age_days(age_days) {
            AgeCategory::Prenatal => 0,
            category => category.order(),
        }
    }

    #[inline]
    #[allow(dead_code)]
    pub fn age_category_name(age_idx: usize) -> &'static str {
        AGE_BUCKETS[age_idx].label()
    }

    #[inline]
    #[allow(dead_code)]
    pub fn default_log_odds(&self, age_idx: usize) -> f64 {
        self.default_age_log_odds[age_idx]
    }

    #[inline]
    pub fn bacteria_age_log_odds(&self, bacteria_idx: usize, age_idx: usize) -> f64 {
        self.bacteria_age_log_odds[bacteria_idx * Self::age_count() + age_idx]
    }

    #[inline]
    #[allow(dead_code)]
    pub fn region_age_log_odds(&self, region: Region, age_idx: usize) -> f64 {
        let region_idx = RegionParameters::region_index(region);
        self.region_age_log_odds[region_idx * Self::age_count() + age_idx]
    }

    #[inline]
    pub fn bacteria_region_age_log_odds(
        &self,
        region: Region,
        bacteria_idx: usize,
        age_idx: usize,
    ) -> f64 {
        let region_idx = RegionParameters::region_index(region);
        let age_count = Self::age_count();
        let idx = (region_idx * self.num_bacteria + bacteria_idx) * age_count + age_idx;
        self.bacteria_region_age_log_odds[idx]
    }
}

#[derive(Debug)]
pub struct HgtMatrix {
    values: Vec<f64>,
    num_bacteria: usize,
}

impl HgtMatrix {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize) -> Self {
        let mut values = Vec::with_capacity(num_bacteria * num_bacteria);
        for &donor in BACTERIA_LIST.iter() {
            for &recipient in BACTERIA_LIST.iter() {
                if donor == recipient {
                    values.push(0.0);
                    continue;
                }
                let key = format!("hgt_prob_{}_to_{}", donor, recipient);
                values.push(get_or_default(map, &key, 0.0));
            }
        }

        HgtMatrix {
            values,
            num_bacteria,
        }
    }

    #[inline]
    pub fn probability(&self, donor_idx: usize, recipient_idx: usize) -> f64 {
        self.values[donor_idx * self.num_bacteria + recipient_idx]
    }
}

#[derive(Debug)]
pub struct ResistanceMechanismParameters {
    pub emergence_rate: Vec<f64>,
    pub enhancement_multiplier: Vec<f64>,
    pub reversion_rate: Vec<f64>,
}

impl ResistanceMechanismParameters {
    fn from_map(map: &HashMap<String, f64>) -> Self {
        let mut emergence_rate = Vec::with_capacity(ResistanceMechanism::all().len());
        let mut enhancement_multiplier = Vec::with_capacity(ResistanceMechanism::all().len());
        let mut reversion_rate = Vec::with_capacity(ResistanceMechanism::all().len());

        for mechanism in ResistanceMechanism::all() {
            let name = mechanism.as_str();
            emergence_rate.push(get_or_default(
                map,
                &format!("resistance_mechanism_{}_emergence_rate", name),
                0.0,
            ));
            enhancement_multiplier.push(get_or_default(
                map,
                &format!("resistance_mechanism_{}_enhancement_multiplier", name),
                0.0,
            ));
            reversion_rate.push(get_or_default(
                map,
                &format!("resistance_mechanism_{}_reversion_rate", name),
                0.0001,
            ));
        }

        ResistanceMechanismParameters {
            emergence_rate,
            enhancement_multiplier,
            reversion_rate,
        }
    }

    #[inline]
    pub fn emergence_rate(&self, mechanism_idx: usize) -> f64 {
        self.emergence_rate[mechanism_idx]
    }

    #[inline]
    pub fn enhancement_multiplier(&self, mechanism_idx: usize) -> f64 {
        self.enhancement_multiplier[mechanism_idx]
    }

    #[inline]
    pub fn reversion_rate(&self, mechanism_idx: usize) -> f64 {
        self.reversion_rate[mechanism_idx]
    }
}

fn get_or_default(map: &HashMap<String, f64>, key: &str, default: f64) -> f64 {
    map.get(key).copied().unwrap_or(default)
}

fn get_required(map: &HashMap<String, f64>, key: &str) -> f64 {
    map.get(key)
        .copied()
        .unwrap_or_else(|| panic!("Missing {} in config", key))
}

fn normalize_label(input: &str) -> String {
    let mut result = String::with_capacity(input.len());
    let mut prev_space = false;
    for ch in input.trim().to_lowercase().chars() {
        match ch {
            'a'..='z' | '0'..='9' | '.' => {
                result.push(ch);
                prev_space = false;
            }
            '_' => {
                if !prev_space {
                    result.push(' ');
                    prev_space = true;
                }
            }
            '-' | '/' => {
                if !prev_space {
                    result.push(' ');
                    prev_space = true;
                }
            }
            ch if ch.is_whitespace() => {
                if !prev_space {
                    result.push(' ');
                    prev_space = true;
                }
            }
            _ => {
                // Ignore other punctuation
            }
        }
    }
    result.trim().to_string()
}

fn normalize_drug_label(input: &str) -> String {
    input.trim().to_lowercase().replace([' ', '-'], "_")
}

fn build_bacteria_lookup() -> HashMap<String, &'static str> {
    let mut lookup = HashMap::new();
    for &name in BACTERIA_LIST.iter() {
        let normalized = normalize_label(name);
        lookup.insert(normalized.clone(), name);
        // Also add space-separated version for compatibility with external data
        lookup.insert(normalized.replace('_', " "), name);
    }
    lookup
}

fn build_drug_lookup() -> HashMap<String, &'static str> {
    let mut lookup = HashMap::new();
    for &drug in DRUG_SHORT_NAMES.iter() {
        let normalized = normalize_drug_label(drug);
        lookup.insert(normalized.clone(), drug);
        lookup.insert(normalized.replace('_', " "), drug);
    }
    lookup
}

// Header count: 52
const POTENCY_EMBEDDED_HEADER: [&str; 52] = [
    "sulfanilamide",
    "penicilling",
    "ampicillin",
    "amoxicillin",
    "piperacillin",
    "ticarcillin",
    "cephalexin",
    "cefazolin",
    "cefuroxime",
    "ceftriaxone",
    "ceftazidime",
    "cefepime",
    "ceftaroline",
    "meropenem",
    "imipenem_c",
    "ertapenem",
    "aztreonam",
    "erythromycin",
    "azithromycin",
    "clarithromycin",
    "clindamycin",
    "gentamicin",
    "tobramycin",
    "amikacin",
    "ciprofloxacin",
    "levofloxacin",
    "moxifloxacin",
    "ofloxacin",
    "tetracycline",
    "doxycycline",
    "minocycline",
    "vancomycin",
    "teicoplanin",
    "dalbavancin",
    "linezolid",
    "tedizolid",
    "quinu_dalfo",
    "trim_sulf",
    "chlorampheni",
    "nitrofurantoin",
    "retapamulin",
    "fusidic_a",
    "metronidazole",
    "furazolidone",
    "rifampicin",
    "amoxicillin_clavulanate",
    "piperacillin_tazobactam",
    "ampicillin_sulbactam",
    "ticarcillin_clavulanate",
    "ceftazidime_avibactam",
    "meropenem_vaborbactam",
    "colistin",
];

const POTENCY_EMBEDDED_DATA: &[(&str, [Option<f64>; 52])] = &[
    (
        "acinetobacter_baumannii",
        [
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.600000),
            Some(0.500000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.100000),
            Some(0.600000),
            Some(0.700000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.750000),
            Some(0.700000),
            Some(0.800000),
            Some(0.700000),
            Some(0.700000),
            Some(0.600000),
            Some(0.600000),
            Some(0.600000),
            Some(0.700000),
            Some(0.800000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.600000),
            Some(0.700000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.600000),
            Some(0.050000),
            Some(0.700000),
            Some(0.700000),
            Some(0.600000),
            Some(0.700000),
            Some(0.800000),
            Some(0.900000),
        ],
    ),
    (
        "citrobacter_spp.",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.800000),
            Some(0.750000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.600000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "enterobacter_spp.",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.500000),
            Some(0.500000),
            Some(0.750000),
            Some(0.700000),
            Some(0.500000),
            Some(0.500000),
            Some(0.600000),
            Some(0.500000),
            Some(0.800000),
            Some(0.850000),
            Some(0.400000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.600000),
            Some(0.700000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.900000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "enterococcus_faecalis",
        [
            Some(0.100000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.750000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.100000),
            Some(0.700000),
            Some(0.800000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.750000),
            Some(0.900000),
            Some(0.700000),
            Some(0.100000),
            Some(0.750000),
            Some(0.050000),
        ],
    ),
    (
        "enterococcus_faecium",
        [
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.750000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.900000),
            Some(0.850000),
            Some(0.850000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
            Some(0.600000),
            Some(0.700000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.100000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
        ],
    ),
    (
        "escherichia_coli",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.900000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.950000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.900000),
            Some(0.970000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "klebsiella_pneumoniae",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.400000),
            Some(0.400000),
            Some(0.800000),
            Some(0.750000),
            Some(0.500000),
            Some(0.500000),
            Some(0.700000),
            Some(0.900000),
            Some(0.850000),
            Some(0.920000),
            Some(0.500000),
            Some(0.940000),
            Some(0.950000),
            Some(0.940000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.600000),
            Some(0.850000),
            Some(0.920000),
            Some(0.750000),
            Some(0.750000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "morganella_spp.",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.500000),
            Some(0.500000),
            Some(0.750000),
            Some(0.700000),
            Some(0.500000),
            Some(0.500000),
            Some(0.600000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.400000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.850000),
            Some(0.700000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.600000),
            Some(0.700000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.900000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "proteus_spp.",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.750000),
            Some(0.850000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "serratia_spp.",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.750000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.600000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.500000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.600000),
            Some(0.700000),
            Some(0.850000),
            Some(0.700000),
            Some(0.750000),
            Some(0.900000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "pseudomonas_aeruginosa",
        [
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.800000),
            Some(0.700000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.900000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.800000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.850000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.800000),
            Some(0.500000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.050000),
            Some(0.900000),
            Some(0.050000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.850000),
        ],
    ),
    (
        "staphylococcus_aureus",
        [
            Some(0.100000),
            Some(0.950000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.600000),
            Some(0.800000),
            Some(0.850000),
            Some(0.700000),
            Some(0.700000),
            Some(0.100000),
            Some(0.600000),
            Some(0.950000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.100000),
            Some(0.800000),
            Some(0.800000),
            Some(0.800000),
            Some(0.800000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.800000),
            Some(0.700000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.600000),
            Some(0.100000),
            Some(0.700000),
            Some(0.050000),
        ],
    ),
    (
        "streptococcus_pneumoniae",
        [
            Some(0.100000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.700000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.050000),
        ],
    ),
    (
        "salmonella_enterica_serovar_typhi",
        [
            Some(0.700000),
            Some(0.100000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "salmonella_enterica_serovar_paratyphi_a",
        [
            Some(0.700000),
            Some(0.100000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "Invasive non-typhoidal Salmonella spp.",
        [
            Some(0.700000),
            Some(0.100000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "shigella_spp.",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.750000),
            Some(0.700000),
            Some(0.600000),
            Some(0.650000),
            Some(0.700000),
            Some(0.900000),
            Some(0.850000),
            Some(0.850000),
            Some(0.600000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.800000),
            Some(0.700000),
            Some(0.850000),
            Some(0.750000),
            Some(0.700000),
            Some(0.800000),
            Some(0.750000),
            Some(0.850000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.900000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.850000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
        ],
    ),
    (
        "neisseria_gonorrhoeae",
        [
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.850000),
            Some(0.800000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.850000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.800000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.850000),
            Some(0.850000),
            Some(0.800000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.050000),
        ],
    ),
    (
        "streptococcus_pyogenes",
        [
            Some(0.100000),
            Some(1.000000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.050000),
        ],
    ),
    (
        "streptococcus_agalactiae",
        [
            Some(0.100000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.050000),
        ],
    ),
    (
        "haemophilus_influenzae",
        [
            Some(0.100000),
            Some(0.700000),
            Some(0.800000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.850000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.700000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.900000),
            Some(0.850000),
            Some(0.900000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.050000),
        ],
    ),
    (
        "chlamydia_trachomatis",
        [
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.950000),
            Some(0.900000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
        ],
    ),
    (
        "vibrio_cholerae",
        [
            Some(0.500000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.900000),
            Some(0.850000),
            Some(0.850000),
            Some(0.700000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.800000),
            Some(0.700000),
            Some(0.800000),
            Some(0.750000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.850000),
            Some(0.900000),
            Some(0.850000),
            Some(0.750000),
            Some(0.850000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.800000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.850000),
            Some(0.900000),
            Some(0.850000),
            Some(0.850000),
            Some(0.900000),
            Some(0.900000),
            Some(0.700000),
        ],
    ),
    (
        "neisseria_meningitidis",
        [
            Some(0.100000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.700000),
            Some(0.800000),
            Some(0.750000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.850000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.850000),
            Some(0.900000),
            Some(0.850000),
            Some(0.900000),
            Some(0.850000),
            Some(0.950000),
            Some(0.950000),
            Some(0.050000),
        ],
    ),
    (
        "listeria_monocytogenes",
        [
            Some(0.100000),
            Some(0.700000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
            Some(0.600000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.100000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.700000),
            Some(0.950000),
            Some(0.600000),
            Some(0.100000),
            Some(0.700000),
            Some(0.050000),
            None,
        ],
    ),
    (
        "clostridioides_difficile",
        [
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.750000),
            Some(0.700000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
        ],
    ),
    (
        "campylobacter_jejuni",
        [
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.800000),
            Some(0.750000),
            Some(0.700000),
            Some(0.750000),
            Some(0.750000),
            Some(0.800000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
        ],
    ),
    (
        "enterobacter_cloacae",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.500000),
            Some(0.500000),
            Some(0.750000),
            Some(0.700000),
            Some(0.500000),
            Some(0.500000),
            Some(0.600000),
            Some(0.400000),
            Some(0.800000),
            Some(0.850000),
            Some(0.400000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.600000),
            Some(0.700000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.900000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "yersinia_enterocolitica",
        [
            Some(0.500000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.750000),
            Some(0.700000),
            Some(0.600000),
            Some(0.650000),
            Some(0.700000),
            Some(0.900000),
            Some(0.850000),
            Some(0.850000),
            Some(0.600000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.950000),
            Some(0.850000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.850000),
            Some(0.850000),
            Some(0.800000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.700000),
        ],
    ),
    (
        "moraxella_catarrhalis",
        [
            Some(0.100000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.800000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.950000),
            Some(0.850000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.950000),
            Some(0.850000),
            Some(0.950000),
            Some(0.850000),
            Some(0.950000),
            Some(0.950000),
            Some(0.050000),
        ],
    ),
    (
        "treponema_pallidum",
        [
            Some(0.100000),
            Some(1.000000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.900000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.750000),
            Some(0.750000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.100000),
            Some(0.950000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.950000),
            Some(0.950000),
            Some(0.050000),
        ],
    ),
    (
        "bordetella_pertussis",
        [
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.700000),
            Some(0.750000),
            Some(0.750000),
            Some(0.700000),
            Some(0.700000),
            Some(0.750000),
            Some(0.750000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.800000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
        ],
    ),
    (
        "helicobacter_pylori",
        [
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.700000),
            Some(0.750000),
            Some(0.700000),
            Some(0.800000),
            Some(0.800000),
            Some(0.850000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.700000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.850000),
            Some(0.100000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
        ],
    ),
    (
        "mdr Mycobacterium tuberculosis",
        [
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.200000),
            Some(0.200000),
            Some(0.200000),
            Some(0.200000),
            Some(0.200000),
            Some(0.250000),
            Some(0.200000),
            Some(0.200000),
            Some(0.250000),
            Some(0.250000),
            Some(0.300000),
            Some(0.400000),
            Some(0.450000),
            Some(0.450000),
            Some(0.400000),
            Some(0.300000),
            Some(0.350000),
            Some(0.350000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.200000),
            Some(0.200000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.200000),
            Some(0.200000),
        ],
    ),
];

fn apply_potency_overrides_from_embedded_table(map: &mut HashMap<String, f64>) {
    let bacteria_lookup = build_bacteria_lookup();
    let drug_lookup = build_drug_lookup();

    let mut column_drugs: Vec<Option<&'static str>> =
        Vec::with_capacity(POTENCY_EMBEDDED_HEADER.len());
    for header in POTENCY_EMBEDDED_HEADER.iter() {
        let normalized = normalize_drug_label(header);
        if let Some(&drug) = drug_lookup.get(&normalized) {
            column_drugs.push(Some(drug));
        } else {
            eprintln!(
                "Unknown drug column '{}' in embedded potency table; skipping",
                header
            );
            column_drugs.push(None);
        }
    }

    let mut override_count = 0usize;

    for (raw_bacteria, potency_values) in POTENCY_EMBEDDED_DATA.iter() {
        let normalized_bacteria = normalize_label(raw_bacteria);
        let Some(&canonical_bacteria) = bacteria_lookup.get(&normalized_bacteria) else {
            eprintln!(
                "Unknown bacteria '{}' in embedded potency table; skipping row",
                raw_bacteria
            );
            continue;
        };

        for (idx, maybe_value) in potency_values.iter().enumerate() {
            let Some(potency) = maybe_value else {
                continue;
            };

            let Some(drug) = column_drugs
                .get(idx)
                .and_then(|entry| entry.as_ref().copied())
            else {
                continue;
            };

            let key = format!(
                "drug_{}_for_bacteria_{}_potency_when_no_r",
                drug, canonical_bacteria
            );
            map.insert(key, *potency);
            override_count += 1;
        }
    }

    println!(
        "Applied {} potency overrides from embedded potency table",
        override_count
    );
}

// --- Global Simulation Parameters ---

// ===========================================================================================================
// Everything inserted into `map` below actively sets defaults that override the fallbacks in from_map.
// ===========================================================================================================

lazy_static! {
    pub static ref PARAMETERS: HashMap<String, f64> = {
        let mut map = HashMap::new();

        // === [A] Bacteria baseline defaults & vaccination scaffolding ===
        // Establishes per-bacteria seed levels, symptom behaviour, and age-aware vaccine priors
        // so scenario templates only need to override deviations instead of rebuilding the grid.
        // --- Default Parameters for ALL Bacteria from BACTERIA_LIST ---
        // These are set first, and can then be overridden by specific entries below.
        for &bacteria in BACTERIA_LIST.iter() {
            map.insert(format!("{}_initial_infection_level", bacteria), 0.01); // 0.01 // bacteria level at initial infection
            map.insert(format!("{}_environmental_acquisition_proportion", bacteria), 0.1); // 0.1  // proportion of new infections from environment
            map.insert(format!("{}_base_bacteria_level_change", bacteria), 0.5); // 0.2 // base change in bacteria level per day
            map.insert(format!("{}_max_level", bacteria), 5.0); // max bacteria level (arbitrary standardized scale)

            // --- Symptom Onset Parameters (Clinical Presentation) ---
            map.insert(format!("{}_daily_symptom_onset_probability", bacteria), 0.15); // Default: 15% chance per day of developing symptoms
            map.insert(format!("{}_symptom_onset_threshold_level", bacteria), 0.5); // Minimum bacteria level needed for symptom onset
            map.insert(format!("{}_symptom_onset_delay_days", bacteria), 1.0); // Minimum days infected before symptoms can start
            map.insert(format!("{}_symptom_onset_level_multiplier", bacteria), 1.0); // How much higher bacteria levels increase symptom probability
            // --- Clearance tuning ---
            // To specialize hazard-based immune clearance, override keys like
            // "{bacteria}_clearance_delay_days" or "{bacteria}_clearance_hazard_multiplier"
            // in template or scenario-specific configurations.

            // Age-related bactera-specific infection risk parameters
            map.insert(format!("{}_age_effect_scaling", bacteria), 1.0); // Scale the template effect (1.0 = full effect)

            // --- Age-specific daily vaccination probability parameters for bacterial vaccines only ---
            // Only vaccines targeting bacteria in our BACTERIA_LIST
            // Age groups: 0-1, 1-5, 5-18, 18-50, 50-70, 70+
            let bacterial_vaccines = vec![
                ("pneumococcal", 1977), // PCV first licensed in 1977 (earlier polysaccharide vaccines)
                ("meningococcal", 1981), // First meningococcal vaccine licensed in 1981
                ("hib", 1985),           // haemophilus_influenzae type b vaccine licensed in 1985
                ("pertussis", 1948),     // DTP vaccine first licensed in 1948
            ];
            for (vaccine, availability_year) in &bacterial_vaccines {
                for category in AGE_BUCKETS.iter() {
                    // Default: 0.0, user should override as needed
                    map.insert(
                        format!(
                            "vaccine_{}_daily_prob_age_{}",
                            vaccine,
                            category.bucket_slug()
                        ),
                        0.0,
                    );
                }
                // Store vaccine availability year for historical modeling
                map.insert(format!("vaccine_{}_availability_year", vaccine), *availability_year as f64);
            }
        }

        // === [B] Horizontal gene transfer (HGT) priors ===
        // Baseline daily probabilities for each donor/recipient pair; tweak here for broad shifts,
        // or override specific pairs in input templates.
        // --- HGT Probabilities for All Donor-Recipient Bacteria Pairs ---
        for &donor in BACTERIA_LIST.iter() {
            for &recipient in BACTERIA_LIST.iter() {
                if donor != recipient {
                    // Default HGT probability (adjust as needed)
                    map.insert(format!("hgt_prob_{}_to_{}", donor, recipient), 0.0001);
                }
            }
        }


        // === [C] Drug initiation, selection, and pharmacokinetics ===
        // Core knobs for therapy behaviour: initiation heuristics, scoring multipliers, and half-lives
        // that drive drug levels. Overwrite these for global experiments; use per-drug keys for specifics.
        // General Drug Parameters
        map.insert("drug_base_initiation_rate_per_day".to_string(), 0.001); // Higher baseline daily initiation to reach usage targets
        map.insert("drug_infection_present_multiplier".to_string(), 260.0); // Encourage more treatment starts when infection detected

        // a non-bacteria-specific parameter that determines how rapidly drugs of a given potency eliminate bacteria level
        // with a value 1 it is nearly always within 1 day 
        map.insert(
            "drug_activity_to_bacteria_level_multiplier".to_string(),
            0.5,
        ); // Global scaling knob for drug-driven bacteria decay
        map.insert("drug_test_identified_multiplier".to_string(), 2.5); // Multiplier when lab diagnostics confirm the pathogen
        map.insert("drug_decay_per_day".to_string(), 1.0); // Legacy parameter - now using drug-specific half-lives

        // Drug Selection Algorithm Parameters
        map.insert("drug_selection_temperature".to_string(), 1.1); // Higher temperature to diversify drug picks and prevent deterministic reserve dominance
        map.insert("reserve_drug_score_penalty".to_string(), 0.001); // Reserve agents get 80% score haircut to keep total usage under stewardship targets

        // Drug-specific half-lives (in days) for realistic pharmacokinetics
        // Beta-lactam/beta-lactamase inhibitor combinations
        map.insert("drug_amoxicillin_clavulanate_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_piperacillin_tazobactam_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ampicillin_sulbactam_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ticarcillin_clavulanate_half_life_days".to_string(), 0.046); // ~1.1 hours
        map.insert("drug_ceftazidime_avibactam_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_meropenem_vaborbactam_half_life_days".to_string(), 0.04); // ~1 hour

        // Polymyxins (Colistin)
        map.insert("drug_colistin_half_life_days".to_string(), 0.08); // ~2 hours

        // Colistin parameters (grouped with other drugs)
        map.insert("drug_colistin_spectrum_breadth".to_string(), 4.0); // Broad spectrum (mainly Gram-negative)
        // Toxicity hazard placeholders (per unit drug level). These represent best-guess daily fatal toxicity odds for active therapy.
        map.insert("drug_colistin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000025); // Colistin-associated nephrotoxicity with high fatal risk
        map.insert("drug_gentamicin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000015); // Aminoglycoside renal failure/ototoxicity
        map.insert("drug_tobramycin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000013);
        map.insert("drug_amikacin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000017);
        map.insert("drug_vancomycin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000006); // Severe nephrotoxicity/red man syndrome rare but serious
        map.insert("drug_chlorampheni_toxicity_death_hazard_per_unit_level".to_string(), 0.00000001); // Aplastic anemia risk
        // Regional availability (assume widely available, adjust as needed)
        map.insert("north_america_drug_colistin_availability".to_string(), 1.0);
        map.insert("europe_drug_colistin_availability".to_string(), 1.0);
        map.insert("asia_drug_colistin_availability".to_string(), 1.0);
        map.insert("oceania_drug_colistin_availability".to_string(), 1.0);
        map.insert("south_america_drug_colistin_availability".to_string(), 1.0);
        map.insert("africa_drug_colistin_availability".to_string(), 1.0);
        map.insert("home_drug_colistin_availability".to_string(), 1.0);

        // Sulfonamides (first antibiotics)
        map.insert("drug_sulfanilamide_half_life_days".to_string(), 0.45); // ~11 hours

        // Beta-lactams (Penicillins)
        map.insert("drug_penicilling_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ampicillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_amoxicillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_piperacillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ticarcillin_half_life_days".to_string(), 0.046); // ~1.1 hours

        // Cephalosporins
        map.insert("drug_cephalexin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_cefazolin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_cefuroxime_half_life_days".to_string(), 0.05); // ~1.3 hours
        map.insert("drug_ceftriaxone_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_ceftazidime_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_cefepime_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_ceftaroline_half_life_days".to_string(), 0.11); // ~2.6 hours

        // Carbapenems
        map.insert("drug_meropenem_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_imipenem_c_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ertapenem_half_life_days".to_string(), 0.17); // ~4 hours

        // Monobactams
        map.insert("drug_aztreonam_half_life_days".to_string(), 0.08); // ~2 hours

        // Macrolides
        map.insert("drug_erythromycin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_azithromycin_half_life_days".to_string(), 2.8); // ~68 hours
        map.insert("drug_clarithromycin_half_life_days".to_string(), 0.25); // ~6 hours

        // Lincosamides
        map.insert("drug_clindamycin_half_life_days".to_string(), 0.125); // ~3 hours

        // Aminoglycosides
        map.insert("drug_gentamicin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_tobramycin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_amikacin_half_life_days".to_string(), 0.08); // ~2 hours

        // Fluoroquinolones
        map.insert("drug_ciprofloxacin_half_life_days".to_string(), 0.17); // ~4 hours
        map.insert("drug_levofloxacin_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_moxifloxacin_half_life_days".to_string(), 0.5); // ~12 hours
        map.insert("drug_ofloxacin_half_life_days".to_string(), 0.25); // ~6 hours

        // Tetracyclines
        map.insert("drug_tetracycline_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_doxycycline_half_life_days".to_string(), 0.75); // ~18 hours
        map.insert("drug_minocycline_half_life_days".to_string(), 0.67); // ~16 hours

        // Glycopeptides
        map.insert("drug_vancomycin_half_life_days".to_string(), 0.25); // ~6 hours
        map.insert("drug_teicoplanin_half_life_days".to_string(), 3.5); // ~83 hours (very long)

        // Oxazolidinones
        map.insert("drug_linezolid_half_life_days".to_string(), 0.21); // ~5 hours
        map.insert("drug_tedizolid_half_life_days".to_string(), 0.5); // ~12 hours

        // Quinolones (older)
        map.insert("drug_quinu_dalfo_half_life_days".to_string(), 0.5); // ~12 hours (quinupristin/dalfopristin)

        // Folate antagonists
        map.insert("drug_trim_sulf_half_life_days".to_string(), 0.5); // ~12 hours (trimethoprim)

        // Other antibiotics
        map.insert("drug_chlorampheni_half_life_days".to_string(), 0.125); // ~3 hours
        map.insert("drug_nitrofurantoin_half_life_days".to_string(), 0.017); // ~20 minutes
        map.insert("drug_retapamulin_half_life_days".to_string(), 0.25); // ~6 hours (topical, limited data)
        map.insert("drug_fusidic_a_half_life_days".to_string(), 0.375); // ~9 hours
        map.insert("drug_metronidazole_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_furazolidone_half_life_days".to_string(), 0.25); // ~6 hours
        map.insert("already_on_drug_initiation_multiplier".to_string(), 1.2); // modest boost for layered therapy when already on treatment
        map.insert("double_dose_probability_if_identified_infection".to_string(), 0.25); // Increased from 0.1 to 0.25 for more aggressive dosing

        // Clinical Decision-Making Potency Thresholds
        map.insert("minimal_potency_threshold_for_drug_selection".to_string(), 0.15); // Higher minimum potency to consider drug (blocks ineffective broad picks)
        map.insert("effective_potency_threshold_for_targeted_therapy".to_string(), 0.10); // Threshold for "good activity" in targeted therapy
        // Drug Evaluation Timing Parameters
        map.insert("drug_evaluation_days_post_infection".to_string(), 7.0); // Number of days after infection to evaluate drug initiation

        // Treatment Failure Assessment Parameters
        map.insert("treatment_failure_assessment_day".to_string(), 4.0); // Days to wait before assessing treatment failure
        map.insert("enable_treatment_failure_assessment".to_string(), 1.0); // Enable/disable treatment failure assessment (1.0=enabled, 0.0=disabled)
        map.insert("drug_failure_memory_days".to_string(), 14.0); // Shorter memory limits rapid escalation to reserve agents
        map.insert("treatment_failure_threshold".to_string(), 0.5); // Threshold for treatment failure (0.5 = failure if bacteria level >= 50% of initial)

        // Restart Window Parameters (for patients who stop drugs early while still infected)
        map.insert("restart_window_days".to_string(), 5.0); // Days after drug cessation to allow "restart" treatment
        map.insert("restart_window_probability".to_string(), 0.3); // Probability that patient returns to care during restart window
        map.insert("restart_bacteria_level_threshold".to_string(), 1.5); // Bacteria level multiplier to trigger restart (current >= cessation * threshold)
        map.insert("enable_restart_window".to_string(), 1.0); // Enable/disable restart window system (1.0=enabled, 0.0=disabled)

        // TB Multi-Drug Synergy Parameters
        // WHY TB NEEDS THESE (and other bacteria don't):
        // 1. TB biology: Intracellular location + thick cell wall + slow metabolism require sustained multi-drug pressure
        // 2. Clinical requirement: WHO guidelines mandate multi-drug therapy; monotherapy = guaranteed treatment failure
        // 3. Pharmacology: TB drugs work synergistically through different mechanisms (cell wall, RNA polymerase, protein synthesis)
        // 4. Resistance prevention: Single drugs lead to rapid resistance; only combination therapy prevents resistance emergence
        // Other bacteria don't have this absolute biological requirement - many can be successfully treated with monotherapy
        map.insert("mdr_mycobacterium_tuberculosis_multi_drug_synergy_threshold".to_string(), 2.0); // Minimum number of active TB drugs for synergy
        map.insert("mdr_mycobacterium_tuberculosis_multi_drug_synergy_multiplier".to_string(), 2.5); // Effectiveness multiplier when ≥2 TB drugs active
        map.insert("mdr_mycobacterium_tuberculosis_background_drug_effectiveness".to_string(), 0.8); // Additional effectiveness from unmodeled TB-specific drugs
        map.insert("mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance".to_string(), 0.90); // Rifampicin resistance level at MDR-TB acquisition

        // Historical MDR TB incidence parameters (time-dependent infection rates)
        map.insert("mdr_tb_pre_antibiotic_era_multiplier".to_string(), 0.0); // Pre-1944: zero MDR TB when antibiotics absent
        map.insert("mdr_tb_early_antibiotic_era_multiplier".to_string(), 0.0);  // 1944-1965: force zero incidence via parameters
        map.insert("mdr_tb_modern_era_multiplier".to_string(), 1.0);            // 1966+: explicitly disabled for diagnostic run
        map.insert("mdr_tb_pre_antibiotic_mortality_multiplier".to_string(), 3.0); // Pre-antibiotic TB had much higher mortality
        map.insert("mdr_tb_ineffective_treatment_mortality_multiplier".to_string(), 2.5); // Ineffective treatment increases mortality

        // === [D] Drug interaction adjustments & therapy flow ===
        // When more than one drug is active, these multipliers capture clinically observed PK interactions
        // and safety-driven dose reductions. Use this section to encode regimen-specific adjustments.
        // --- Drug Level Interaction Parameters ---
        // These are pairwise interactions that affect the effective level of each drug when co-administered
        // Format: "drug_level_multiplier_{drug1}_when_coadministered_with_{drug2}" -> multiplier for drug1's level
        // Important: These are NOT activity bonuses or generic effects - they model actual pharmacokinetic interactions

        // Rifampicin interactions (CYP450 induction reduces levels of co-administered drugs)
        // Rifampicin is a potent CYP3A4 inducer - reduces levels of drugs metabolized by this pathway
        map.insert("drug_level_multiplier_levofloxacin_when_coadministered_with_rifampicin".to_string(), 0.7); // 30% reduction in levofloxacin levels
        map.insert("drug_level_multiplier_moxifloxacin_when_coadministered_with_rifampicin".to_string(), 0.8); // 20% reduction in moxifloxacin levels
        map.insert("drug_level_multiplier_clarithromycin_when_coadministered_with_rifampicin".to_string(), 0.6); // 40% reduction in clarithromycin levels
        map.insert("drug_level_multiplier_azithromycin_when_coadministered_with_rifampicin".to_string(), 0.8); // 20% reduction in azithromycin levels

        // Fluoroquinolone interactions with divalent cations (from combination drugs)
        // Amoxicillin/clavulanate contains potassium clavulanate which can reduce fluoroquinolone absorption
        map.insert("drug_level_multiplier_ciprofloxacin_when_coadministered_with_amoxicillin_clavulanate".to_string(), 0.85); // 15% reduction
        map.insert("drug_level_multiplier_levofloxacin_when_coadministered_with_amoxicillin_clavulanate".to_string(), 0.9); // 10% reduction

        // Macrolide-fluoroquinolone interactions (QT prolongation leads to dose reduction)
        // Clinical practice often reduces doses when combining these drug classes due to cardiac safety
        map.insert("drug_level_multiplier_ciprofloxacin_when_coadministered_with_erythromycin".to_string(), 0.85); // Dose reduction for safety
        map.insert("drug_level_multiplier_levofloxacin_when_coadministered_with_azithromycin".to_string(), 0.9); // Dose reduction for safety

        // === [E] Drug-bacteria potency & emergence settings ===
        // Qualitative potency buckets, initiation multipliers, and baseline resistance emergence
        // rates for every drug/bacteria pair. Override here for wide shifts, or patch specific
        // entries when ingesting empirical potency tables.

        
    // ---------------- 9) Drug-bacteria potency & cross-resistance mappings ----------------
        // --- Drug-Bacteria Potency Matrix: Evidence-Based Approach ---
        // Instead of uniform potency, use clinically relevant potency categories:
        // 1.00+ = Excellent potency (first-line therapy)
        // 0.50-0.99 = Good potency (reliable option)
        // 0.25-0.49 = Moderate potency (situational use)
        // 0.05-0.24 = Poor potency (usually ineffective)
        // 0.05 = Very poor/no activity

        // Potency values stay on this qualitative scale rather than raw 1/MIC inputs.

        // Define drug classes for easier management
        // Polymyxins (currently only Colistin)
        let polymyxins = vec!["colistin"];
        let penicillins = vec!["penicilling", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin",
            // BL/BLI combinations
            "amoxicillin_clavulanate", "piperacillin_tazobactam", "ampicillin_sulbactam", "ticarcillin_clavulanate"
        ];
        let cephalosporins_1_2 = vec!["cephalexin", "cefazolin", "cefuroxime"];
        let _cephalosporins_3_4 = vec!["ceftriaxone", "ceftazidime", "cefepime", "ceftaroline"];
        let cephalosporins_3_4 = vec!["ceftriaxone", "ceftazidime", "cefepime", "ceftaroline",
            // BL/BLI cephalosporin
            "ceftazidime_avibactam"
        ];
        // let _carbapenems = vec!["meropenem", "imipenem_c", "ertapenem"];
        let carbapenems = vec!["meropenem", "imipenem_c", "ertapenem",
            // BL/BLI carbapenem
            "meropenem_vaborbactam"
        ];
        let _monobactams = vec!["aztreonam"];
        let macrolides = vec!["erythromycin", "azithromycin", "clarithromycin"];
        let _lincosamides = vec!["clindamycin"];
        let aminoglycosides = vec!["gentamicin", "tobramycin", "amikacin"];
        let fluoroquinolones = vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"];
        let tetracyclines = vec!["tetracycline", "doxycycline", "minocycline"];
        let glycopeptides = vec!["vancomycin", "teicoplanin"]; // TODO: add dalbavancin potency overrides once parameters are curated
        let oxazolidinones = vec!["linezolid", "tedizolid"];
        let _folate_antagonists = vec!["trim_sulf"];
        let _other_antibiotics = vec!["quinu_dalfo", "chlorampheni", "nitrofurantoin", "retapamulin", "fusidic_a", "metronidazole", "furazolidone"];

        // Define bacterial groups for potency patterns - all names use underscores consistently
        let gram_pos_cocci = vec!["staphylococcus_aureus", "staphylococcus_epidermidis", "streptococcus_pneumoniae", "streptococcus_pyogenes", "streptococcus_agalactiae", "enterococcus_faecalis", "enterococcus_faecium"];
        let gram_neg_enterobacteria = vec!["escherichia_coli", "klebsiella_pneumoniae", "enterobacter_spp.", "citrobacter_spp.", "serratia_spp.", "proteus_spp.", "morganella_spp.", "enterobacter_cloacae"];
        let gram_neg_non_fermenting = vec!["pseudomonas_aeruginosa", "acinetobacter_baumannii", "stenotrophomonas_maltophilia"];
        let fastidious_gram_neg = vec!["haemophilus_influenzae", "moraxella_catarrhalis", "neisseria_gonorrhoeae", "neisseria_meningitidis", "bordetella_pertussis"];
        let enteric_pathogens = vec!["salmonella_enterica_serovar_typhi", "salmonella_enterica_serovar_paratyphi_a", "invasive_non-typhoidal_salmonella_spp.", "shigella_spp.", "vibrio_cholerae", "campylobacter_jejuni", "yersinia_enterocolitica"];
        let atypical_pathogens = vec!["chlamydia_trachomatis"];
        let anaerobes_spore_formers = vec!["clostridioides_difficile"];
        let gram_pos_rods = vec!["listeria_monocytogenes"];
        let gastric_pathogens = vec!["helicobacter_pylori"]; // Unique microaerophilic Gram-negative

        for &drug in DRUG_SHORT_NAMES.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                // Bacteria names now use underscores consistently
                map.insert(format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug, bacteria), 1.0);
                map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.1); // Default low potency 0.1
                // Default resistance emergence rate: calibrated to generate ~20-25% resistance prevalence
                // This is multiplied by bacteria_level_factor, drug_concentration_factor, and multi_drug_penalty
                // Note: Higher than pure mutation rates because represents combined selection + mutation + amplification
                // CALIBRATION: 0.0001 gave 4%, 0.001 gave 37%, 0.0003 gave 15% - trying 0.0005 to aim for ~23%
                // At 0.0005 baseline with modifiers, expect ~0.005-0.025% daily emergence probability when on treatment
                map.insert(format!("drug_{}_for_bacteria_{}_resistance_emergence_rate_per_day_baseline", drug, bacteria), 0.001);
            }
        }

        // === GLOBAL CARBAPENEM AND RESERVE ANTIBIOTIC PENALTIES ===
        // Carbapenems and other reserve agents should have very low base initiation rates
        // to reflect antimicrobial stewardship principles. These are "last resort" drugs
        // that should only be used when first-line agents fail or resistance is documented.
        // Apply across ALL bacteria to ensure consistent stewardship behavior.
        let carbapenem_reserve_drugs = vec![
            "meropenem",
            "meropenem_vaborbactam",
            "imipenem_c",
            "ertapenem",
            "colistin",
            "linezolid",
            "tedizolid",
            "quinu_dalfo",
            "dalbavancin",
            "ceftazidime_avibactam",
        ];
        for &drug in carbapenem_reserve_drugs.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                // CALIBRATION: Reduced from 0.05 to 0.005 (200x less likely than default) to keep reserve usage <5%
                // This will be overridden by specific bacteria-drug combinations where appropriate
                map.insert(
                    format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug, bacteria),
                    0.005,
                );
            }
        }

        // examples of how to override potencies
        // Specific override: sulfanilamide has high potency against haemophilus_influenzae
        // map.insert("drug_sulfanilamide_for_bacteria_haemophilus_influenzae_potency_when_no_r".to_string(), 0.85); // Example high potency

        // Sulfanilamide - historically effective against specific pathogens
        for &drug in DRUG_SHORT_NAMES.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                if drug == "sulfanilamide" {
                    let potency = match bacteria {
                        // Excellent against streptococci (historical primary indication)
                        bacteria if bacteria.contains("streptococcus") => 0.85,
                        // Good against E. coli (UTI treatment)
                        "escherichia_coli" => 0.65,
                        // Moderate against other gram-positives
                        "staphylococcus_aureus" => 0.45,
                        // Limited against enterococci and most gram-negatives
                        _ => 0.20,
                    };
                    map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                }
            }
        }

        // Set specific potencies based on clinical evidence

        // GRAM-POSITIVE COCCI (Staph, Strep, Enterococcus)
        for &bacteria in gram_pos_cocci.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Penicillins - excellent for Strep, strong for E. faecalis, limited for Staph / E. faecium
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("streptococcus") {
                            0.95
                        } else if bacteria == "enterococcus_faecalis" {
                            0.75
                        } else if bacteria == "enterococcus_faecium" {
                            0.20
                        } else if bacteria == "staphylococcus_aureus" {
                            0.15
                        } else if bacteria == "staphylococcus_epidermidis" {
                            0.05
                        } else {
                            0.10
                        };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Cephalosporins - good for most gram-positive (except Enterococcus)
                for &drug in cephalosporins_1_2.iter().chain(cephalosporins_3_4.iter()) {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("enterococcus") { 0.05 } else { 0.75 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Carbapenems - good but reserve for resistant cases
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("enterococcus") { 0.25 } else { 0.80 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Macrolides - good for Strep and atypicals
                for &drug in macrolides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.60);
                    }
                }

                // Glycopeptides - excellent for gram-positive, especially MRSA/VRE
                for &drug in glycopeptides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 1.00);
                    }
                }

                // Oxazolidinones - excellent for resistant gram-positive
                for &drug in oxazolidinones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 1.10);
                    }
                }
            }
        }

        // GRAM-NEGATIVE ENTEROBACTERIA (E. coli, Klebsiella, etc.)
        for &bacteria in gram_neg_enterobacteria.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Penicillins - BL/BLI combinations retain activity, older agents mostly ineffective
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = match drug {
                            "penicilling" => 0.05,
                            "ampicillin" | "amoxicillin" => 0.20,
                            "ampicillin_sulbactam" => 0.55,
                            "amoxicillin_clavulanate" => 0.50,
                            "piperacillin" => 0.75,
                            "piperacillin_tazobactam" => 0.95,
                            "ticarcillin" => 0.15,
                            "ticarcillin_clavulanate" => 0.70,
                            _ => 0.10,
                        };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Cephalosporins - variable by generation
                for &drug in cephalosporins_1_2.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.40);
                    }
                }
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.80);
                    }
                }

                // Carbapenems - excellent broad-spectrum
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 1.05);
                    }
                }
                // Polymyxins (Colistin) - only active against Gram-negatives
                for &drug in polymyxins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = match bacteria {
                            // Gram-positive bacteria - intrinsically resistant (colistin doesn't penetrate Gram-positive cell walls)
                            "enterococcus_faecalis" | "enterococcus_faecium" | "staphylococcus_aureus" |
                            "streptococcus_pneumoniae" | "streptococcus_pyogenes" | "streptococcus_agalactiae" |
                            "listeria_monocytogenes" | "clostridioides_difficile" => 0.0,

                            // Gram-negative intrinsic resistance
                            "morganella_spp." | "proteus_spp." | "serratia_spp." => 0.0,

                            // Gram-negative with variable/reduced susceptibility
                            "salmonella_enterica_serovar_typhi" | "salmonella_enterica_serovar_paratyphi_a" |
                            "invasive_non-typhoidal_salmonella_spp." | "shigella_spp." => 0.5,

                            // Gram-negative normally susceptible
                            "vibrio_cholerae" | "yersinia_enterocolitica" => 1.0,

                            // Most other Gram-negatives (E. coli, Klebsiella, Pseudomonas, Acinetobacter, etc.)
                            _ => 1.0,
                        };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Fluoroquinolones - good broad-spectrum
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.85);
                    }
                }

                // Aminoglycosides - good for serious infections
                for &drug in aminoglycosides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.75);
                    }
                }

                // Trim-sulf - moderate activity
                if DRUG_SHORT_NAMES.contains(&"trim_sulf") {
                    map.insert(format!("drug_trim_sulf_for_bacteria_{}_potency_when_no_r", bacteria), 0.50);
                }
            }
        }

        // klebsiella_pneumoniae - practical inactivity for most early beta-lactams
        let kleb_low_pen = vec![
            "penicilling", "ampicillin", "amoxicillin", "ampicillin_sulbactam",
            "amoxicillin_clavulanate", "ticarcillin", "ticarcillin_clavulanate"
        ];
        for &drug in kleb_low_pen.iter() {
            if DRUG_SHORT_NAMES.contains(&drug) {
                map.insert(format!("drug_{}_for_bacteria_klebsiella_pneumoniae_potency_when_no_r", drug), 0.03);
            }
        }
        let kleb_low_pen_advanced = vec!["piperacillin", "piperacillin_tazobactam"];
        for &drug in kleb_low_pen_advanced.iter() {
            if DRUG_SHORT_NAMES.contains(&drug) {
                map.insert(format!("drug_{}_for_bacteria_klebsiella_pneumoniae_potency_when_no_r", drug), 0.05);
            }
        }
        let kleb_ceph12 = vec!["cephalexin", "cefazolin", "cefuroxime"];
        for &drug in kleb_ceph12.iter() {
            if DRUG_SHORT_NAMES.contains(&drug) {
                map.insert(format!("drug_{}_for_bacteria_klebsiella_pneumoniae_potency_when_no_r", drug), 0.05);
            }
        }

        // PSEUDOMONAS & ACINETOBACTER (Non-fermenting gram-negatives)
        for &bacteria in gram_neg_non_fermenting.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Most beta-lactams poor except specific anti-pseudomonal agents
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = match drug {
                            "piperacillin" => 0.65,
                            "piperacillin_tazobactam" | "ticarcillin" | "ticarcillin_clavulanate" | "ampicillin_sulbactam" => 0.25,
                            _ => 0.025,
                        };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Only specific cephalosporins active
                for &drug in cephalosporins_1_2.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.025);
                    }
                }
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = match drug {
                            "ceftazidime" | "cefepime" => 0.75,
                            "ceftazidime_avibactam" => 0.9,
                            _ => 0.10,
                        };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Carbapenems - good but resistance emerging
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("acinetobacter") { 0.60 } else { 0.80 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                // Polymyxins (Colistin) - high potency for non-fermenters
                for &drug in polymyxins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 1.0);
                    }
                }

                // Fluoroquinolones - good activity
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.75);
                    }
                }

                // Aminoglycosides - good for combination therapy
                for &drug in aminoglycosides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.70);
                    }
                }
            }
        }

        // FASTIDIOUS GRAM-NEGATIVE (H. flu, Moraxella, Neisseria, Bordetella)
        for &bacteria in fastidious_gram_neg.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Penicillins - poor due to beta-lactamase (except amoxicillin-clavulanate)
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if drug == "amoxicillin" { 0.20 } else { 0.05 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Beta-lactam/beta-lactamase inhibitor combinations - excellent
                if DRUG_SHORT_NAMES.contains(&"amoxicillin_clavulanate") {
                    map.insert(format!("drug_amoxicillin_clavulanate_for_bacteria_{}_potency_when_no_r", bacteria), 0.85);
                }

                // Cephalosporins - 2nd/3rd gen good, 1st gen poor
                for &drug in cephalosporins_1_2.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if drug == "cephalexin" || drug == "cefazolin" { 0.15 } else { 0.75 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.80);
                    }
                }

                // Macrolides - excellent for respiratory pathogens
                for &drug in macrolides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.75);
                    }
                }

                // Fluoroquinolones - good activity
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.70);
                    }
                }
            }
        }

        // neisseria_gonorrhoeae potency tuning - cephalosporins remain highly effective, FQs moderate
        let gonorrhea_cephs = vec!["ceftriaxone", "ceftazidime", "cefepime", "ceftazidime_avibactam"];
        for &drug in gonorrhea_cephs.iter() {
            if DRUG_SHORT_NAMES.contains(&drug) {
                map.insert(format!("drug_{}_for_bacteria_neisseria_gonorrhoeae_potency_when_no_r", drug), 0.55);
            }
        }
        let gonorrhea_fqs = vec!["ciprofloxacin", "levofloxacin", "ofloxacin", "moxifloxacin"];
        for &drug in gonorrhea_fqs.iter() {
            if DRUG_SHORT_NAMES.contains(&drug) {
                map.insert(format!("drug_{}_for_bacteria_neisseria_gonorrhoeae_potency_when_no_r", drug), 0.25);
            }
        }

        // ENTERIC PATHOGENS (Salmonella, Shigella, Vibrio, Campylobacter, Yersinia)
        for &bacteria in enteric_pathogens.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Penicillins - poor (intrinsic resistance in many)
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.05);
                    }
                }

                // Cephalosporins - good for most
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.75);
                    }
                }

                // Fluoroquinolones - excellent for enteric pathogens
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.85);
                    }
                }

                // Tetracyclines - good for many enteric pathogens
                for &drug in tetracyclines.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.70);
                    }
                }
            }
        }

        // ATYPICAL PATHOGENS (Chlamydia)
        for &bacteria in atypical_pathogens.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Beta-lactams - no activity (no cell wall)
                for &drug in penicillins.iter().chain(cephalosporins_1_2.iter()).chain(cephalosporins_3_4.iter()).chain(carbapenems.iter()) {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.01);
                    }
                }

                // Macrolides - excellent (first-line)
                for &drug in macrolides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.90);
                    }
                }

                // Tetracyclines - excellent alternative
                for &drug in tetracyclines.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.85);
                    }
                }

                // Fluoroquinolones - good alternative
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.80);
                    }
                }
            }
        }

        // ANAEROBES/SPORE FORMERS (C. difficile)
        for &bacteria in anaerobes_spore_formers.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Most antibiotics - poor (anaerobic environment)
                for &drug in penicillins.iter().chain(cephalosporins_1_2.iter()).chain(cephalosporins_3_4.iter()).chain(aminoglycosides.iter()) {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.05);
                    }
                }

                // Vancomycin - first-line for C. diff
                if DRUG_SHORT_NAMES.contains(&"vancomycin") {
                    map.insert(format!("drug_vancomycin_for_bacteria_{}_potency_when_no_r", bacteria), 0.90);
                }

                // Metronidazole - good for anaerobes
                if DRUG_SHORT_NAMES.contains(&"metronidazole") {
                    map.insert(format!("drug_metronidazole_for_bacteria_{}_potency_when_no_r", bacteria), 0.80);
                }
            }
        }

        // --- neisseria_meningitidis SPECIFIC POTENCY OVERRIDES ---
        // N. meningitidis is typically very penicillin-sensitive, unlike other fastidious gram-negatives
        // Clinical reality: Penicillin G and ampicillin are first-line therapies for sensitive strains
        if BACTERIA_LIST.contains(&"neisseria_meningitidis") {
            // Penicillins - EXCELLENT activity (first-line therapy)
            map.insert("drug_penicilling_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.95); // Excellent first-line
            map.insert("drug_ampicillin_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.90); // Excellent alternative

            // Ensure ceftriaxone remains excellent (current first-line)
            map.insert("drug_ceftriaxone_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.95); // Current first-line

            // Chloramphenicol - historically important alternative
            if DRUG_SHORT_NAMES.contains(&"chlorampheni") {
                map.insert("drug_chlorampheni_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.85); // Good alternative
            }

            // Rifampin - used for prophylaxis (not treatment)
            if DRUG_SHORT_NAMES.contains(&"rifampicin") {
                map.insert("drug_rifampicin_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.80); // Good for prophylaxis
            }

            // Ciprofloxacin - critical for prophylaxis and some treatment
            if DRUG_SHORT_NAMES.contains(&"ciprofloxacin") {
                map.insert("drug_ciprofloxacin_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.90); // Excellent for prophylaxis
            }

            // Levofloxacin - also excellent for meningococcal infections
            if DRUG_SHORT_NAMES.contains(&"levofloxacin") {
                map.insert("drug_levofloxacin_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.85); // Very good alternative
            }

            // Cefotaxime - alternative third-generation cephalosporin
            if DRUG_SHORT_NAMES.contains(&"cefotaxime") {
                map.insert("drug_cefotaxime_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.95); // Excellent alternative to ceftriaxone
            }

            // Sulfanilamide - historically used but resistance common
            if DRUG_SHORT_NAMES.contains(&"sulfanilamide") {
                map.insert("drug_sulfanilamide_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.10); // High resistance rates
            }

            // Reduce inappropriately high inherited macrolide values
            if DRUG_SHORT_NAMES.contains(&"azithromycin") {
                map.insert("drug_azithromycin_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.50); // Lower than inherited 0.75
            }
            if DRUG_SHORT_NAMES.contains(&"clarithromycin") {
                map.insert("drug_clarithromycin_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.50); // Lower than inherited 0.75
            }
        }

        // GRAM-POSITIVE RODS (Listeria)
        for &bacteria in gram_pos_rods.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Penicillins - excellent (first-line)
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.90);
                    }
                }

                // Cephalosporins - poor (intrinsic resistance)
                for &drug in cephalosporins_1_2.iter().chain(cephalosporins_3_4.iter()) {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.05);
                    }
                }

                // Carbapenems - good
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.80);
                    }
                }
            }
        }

        // GASTRIC PATHOGENS (H. pylori)
        for &bacteria in gastric_pathogens.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Most antibiotics - poor in gastric environment
                for &drug in penicillins.iter().chain(cephalosporins_1_2.iter()).chain(cephalosporins_3_4.iter()) {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if drug == "amoxicillin" { 0.70 } else { 0.05 }; // Amoxicillin exception for H. pylori
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                // Specific H. pylori drugs handled in individual overrides below
            }
        }





        // Add specific high-potency combinations for clinical effectiveness
        // These represent particularly effective drug-bacteria pairs

        // Azithromycin for atypicals and some enteric pathogens
        if DRUG_SHORT_NAMES.contains(&"azithromycin") {
            for &bacteria in &["chlamydia_trachomatis", "campylobacter_jejuni"] {
                if BACTERIA_LIST.contains(&bacteria) {
                    map.insert(format!("drug_azithromycin_for_bacteria_{}_potency_when_no_r", bacteria), 1.25);
                }
            }
        }

        if DRUG_SHORT_NAMES.contains(&"clarithromycin") && BACTERIA_LIST.contains(&"campylobacter_jejuni") {
            map.insert(
                "drug_clarithromycin_for_bacteria_campylobacter_jejuni_potency_when_no_r".to_string(),
                1.05,
            );
        }
        if DRUG_SHORT_NAMES.contains(&"erythromycin") && BACTERIA_LIST.contains(&"campylobacter_jejuni") {
            map.insert(
                "drug_erythromycin_for_bacteria_campylobacter_jejuni_potency_when_no_r".to_string(),
                0.95,
            );
        }

        // Nitrofurantoin for urinary E. coli
        if DRUG_SHORT_NAMES.contains(&"nitrofurantoin") && BACTERIA_LIST.contains(&"escherichia_coli") {
            map.insert("drug_nitrofurantoin_for_bacteria_escherichia_coli_potency_when_no_r".to_string(), 0.95);
        }

        // Metronidazole for anaerobes
        if DRUG_SHORT_NAMES.contains(&"metronidazole") && BACTERIA_LIST.contains(&"clostridioides_difficile") {
            map.insert("drug_metronidazole_for_bacteria_clostridioides_difficile_potency_when_no_r".to_string(), 0.90);
        }





        // --- TARGETED CLINICAL FIXES FOR SPECIFIC BACTERIA ---
        // Conservative adjustments for clear clinical issues while preserving regional resistance surveillance

        // vibrio_cholerae - Fix cholera treatment (should be tetracyclines + fluoroquinolones)
        // Boost first-line drugs moderately
        map.insert("drug_doxycycline_for_bacteria_vibrio_cholerae_potency_when_no_r".to_string(), 0.9);  // Excellent activity
        map.insert("drug_tetracycline_for_bacteria_vibrio_cholerae_potency_when_no_r".to_string(), 0.85); // Excellent activity
        map.insert("drug_ciprofloxacin_for_bacteria_vibrio_cholerae_potency_when_no_r".to_string(), 0.8); // Very good activity
        map.insert("drug_levofloxacin_for_bacteria_vibrio_cholerae_potency_when_no_r".to_string(), 0.8);  // Very good activity
        // Reduce inappropriate drugs slightly (regional resistance surveillance will handle local patterns)
        map.insert("drug_penicilling_for_bacteria_vibrio_cholerae_potency_when_no_r".to_string(), 0.05);  // Poor activity
        map.insert("drug_ampicillin_for_bacteria_vibrio_cholerae_potency_when_no_r".to_string(), 0.05);  // Poor activity

        // haemophilus_influenzae - Address beta-lactamase resistance (intrinsic in many strains)
        // Reduce basic penicillins (H. flu commonly produces beta-lactamase)
        map.insert("drug_penicilling_for_bacteria_haemophilus_influenzae_potency_when_no_r".to_string(), 0.03); // Poor due to beta-lactamase
        map.insert("drug_ampicillin_for_bacteria_haemophilus_influenzae_potency_when_no_r".to_string(), 0.15);  // Reduced due to beta-lactamase resistance
        // Boost appropriate alternatives modestly
        map.insert("drug_amoxicillin_clavulanate_for_bacteria_haemophilus_influenzae_potency_when_no_r".to_string(), 0.75); // Good activity with beta-lactamase inhibitor
        map.insert("drug_cefuroxime_for_bacteria_haemophilus_influenzae_potency_when_no_r".to_string(), 0.7);   // Good 2nd gen cephalosporin
        map.insert("drug_azithromycin_for_bacteria_haemophilus_influenzae_potency_when_no_r".to_string(), 0.65); // Good macrolide activity

        // bordetella_pertussis - Macrolide-first treatment (azithromycin, clarithromycin, erythromycin)
        map.insert("bordetella_pertussis_base_bacteria_level_change".to_string(), 0.3); // Catarrhal phase builds over ~1 week
        // Boost macrolides (first-line treatment for pertussis)
        map.insert("drug_azithromycin_for_bacteria_bordetella_pertussis_potency_when_no_r".to_string(), 0.9);   // Excellent activity, first-line
        map.insert("drug_clarithromycin_for_bacteria_bordetella_pertussis_potency_when_no_r".to_string(), 0.85); // Excellent activity, first-line
        map.insert("drug_erythromycin_for_bacteria_bordetella_pertussis_potency_when_no_r".to_string(), 0.8);    // Good activity, traditional first-line
        map.insert("drug_trim_sulf_for_bacteria_bordetella_pertussis_potency_when_no_r".to_string(), 0.7);       // Alternative for macrolide-allergic patients
        // Reduce inappropriate antibiotics
        map.insert("drug_penicilling_for_bacteria_bordetella_pertussis_potency_when_no_r".to_string(), 0.05);     // Poor activity
        map.insert("drug_ampicillin_for_bacteria_bordetella_pertussis_potency_when_no_r".to_string(), 0.05);     // Poor activity

        // helicobacter_pylori - Triple/quadruple therapy drugs (clarithromycin + amoxicillin + metronidazole)
        // Historical note: H. pylori discovered 1982 by Marshall & Warren, triple therapy established ~1990s
        // Before 1982, peptic ulcers attributed to stress/diet, not infection
        // Boost first-line eradication therapy drugs
        map.insert("drug_clarithromycin_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.85);   // Key component of triple therapy
        map.insert("drug_amoxicillin_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.8);       // Key component of triple therapy
        map.insert("drug_metronidazole_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.75);    // Alternative/quadruple therapy
        map.insert("drug_tetracycline_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.7);      // Bismuth quadruple therapy
        map.insert("drug_levofloxacin_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.75);     // Rescue therapy

        // H. pylori drug availability dates - triple therapy not available until 1990s
        map.insert("drug_clarithromycin_for_bacteria_helicobacter_pylori_availability_year".to_string(), 1990.0); // Triple therapy era
        map.insert("drug_amoxicillin_for_bacteria_helicobacter_pylori_availability_year".to_string(), 1990.0);    // Triple therapy era
        map.insert("drug_metronidazole_for_bacteria_helicobacter_pylori_availability_year".to_string(), 1990.0);  // Triple therapy era

        // H. pylori testing parameters - endoscopy/biopsy required for identification
        map.insert("helicobacter_pylori_test_probability_per_day".to_string(), 0.15); // Higher testing rate for chronic symptoms
        map.insert("helicobacter_pylori_test_sensitivity".to_string(), 0.95);        // High sensitivity with proper testing
        // Bacteria-specific test availability years for late-discovered organisms
        map.insert("helicobacter_pylori_test_availability_year".to_string(), 1982.0); // Marshall & Warren discovery
        map.insert("chlamydia_trachomatis_test_availability_year".to_string(), 1966.0); // Chlamydia cell culture methods
        map.insert("campylobacter_jejuni_test_availability_year".to_string(), 1973.0); // Campylobacter isolation techniques

        // Bacteria-specific treatment recognition years (when bacteria was first recognized as needing treatment)
        map.insert("helicobacter_pylori_treatment_recognition_year".to_string(), 1982.0); // H. pylori not treated before Marshall & Warren discovery

        // Bacteria-specific sepsis risk overrides for organisms that don't cause acute sepsis are defined in the log-odds section below
        map.insert("helicobacter_pylori_base_bacteria_level_change".to_string(), 0.2); // Slow-growing chronic colonizer

        // H. pylori-specific drug selection bonuses when bacteria is identified
        map.insert("drug_clarithromycin_for_bacteria_helicobacter_pylori_initiation_multiplier".to_string(), 15.0); // Strong preference for triple therapy
        map.insert("drug_amoxicillin_for_bacteria_helicobacter_pylori_initiation_multiplier".to_string(), 12.0);   // Strong preference for triple therapy
        map.insert("drug_metronidazole_for_bacteria_helicobacter_pylori_initiation_multiplier".to_string(), 8.0);  // Alternative therapy

        // N. meningitidis-specific drug selection bonuses - emergency treatment for meningococcal disease
        map.insert("drug_penicilling_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 25.0); // First-line for sensitive strains
        map.insert("drug_ampicillin_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 22.0);  // Excellent alternative
        map.insert("drug_ceftriaxone_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 30.0); // Current standard of care
        map.insert("drug_cefotaxime_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 28.0);  // Equivalent 3rd generation
        map.insert("drug_chlorampheni_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 18.0); // Important alternative historically
        map.insert("drug_ciprofloxacin_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 15.0); // Prophylaxis and treatment
        map.insert("drug_rifampicin_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 12.0);  // Prophylaxis agent

        // Reduce inappropriate antibiotics
        map.insert("drug_penicilling_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.05);       // Not used for H. pylori
        map.insert("drug_cephalexin_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.05);       // Not effective

        // --- BACTERIA-SPECIFIC SYMPTOM ONSET PARAMETERS ---

        // H. PYLORI - Usually asymptomatic chronic gastritis
        map.insert("helicobacter_pylori_daily_symptom_onset_probability".to_string(), 0.001); // 0.1% per day - very low symptomatic rate
        map.insert("helicobacter_pylori_symptom_onset_threshold_level".to_string(), 2.0);     // High threshold for symptoms
        map.insert("helicobacter_pylori_symptom_onset_delay_days".to_string(), 30.0);         // Long delay before symptoms possible

        // chlamydia_trachomatis - Often asymptomatic
        map.insert("chlamydia_trachomatis_daily_symptom_onset_probability".to_string(), 0.01); // 1% per day - often asymptomatic
        map.insert("chlamydia_trachomatis_symptom_onset_threshold_level".to_string(), 1.5);    // Moderate threshold
        map.insert("chlamydia_trachomatis_base_bacteria_level_change".to_string(), 0.3);       // Slow intracellular replication

        // neisseria_meningitidis - Often asymptomatic carriage
        map.insert("neisseria_meningitidis_daily_symptom_onset_probability".to_string(), 0.25); // 25% per day - increase clinical visibility
        map.insert("neisseria_meningitidis_base_bacteria_level_change".to_string(), 0.65);      // Fulminant meningococcemia progression
        map.insert("neisseria_meningitidis_symptom_onset_threshold_level".to_string(), 3.0);    // High threshold for invasive disease

        // moraxella_catarrhalis - Often just colonization
        map.insert("moraxella_catarrhalis_daily_symptom_onset_probability".to_string(), 0.05);  // 5% per day - often colonizer
        map.insert("moraxella_catarrhalis_symptom_onset_threshold_level".to_string(), 2.0);     // Moderate threshold
        map.insert("moraxella_catarrhalis_base_bacteria_level_change".to_string(), 0.55);       // Rapid otitis/sinusitis onset in children

    // pseudomonas_aeruginosa - Clinically apparent when burden high
        map.insert("pseudomonas_aeruginosa_daily_symptom_onset_probability".to_string(), 0.2);   // 20% per day - improve detection of invasive disease
        map.insert("pseudomonas_aeruginosa_symptom_onset_threshold_level".to_string(), 0.8);     // Higher burden needed before symptoms manifest
        map.insert("pseudomonas_aeruginosa_base_bacteria_level_change".to_string(), 0.55);       // Rapid proliferation in ventilated hosts

        // ACUTE INFECTIONS - High symptomatic rates
        map.insert("haemophilus_influenzae_base_bacteria_level_change".to_string(), 0.55);       // Rapid pediatric respiratory progression
        map.insert("streptococcus_pneumoniae_daily_symptom_onset_probability".to_string(), 0.8); // 80% per day - usually symptomatic
        map.insert("streptococcus_pyogenes_daily_symptom_onset_probability".to_string(), 0.7);   // 70% per day - usually symptomatic
        map.insert("streptococcus_pyogenes_base_bacteria_level_change".to_string(), 0.6);         // Fast doubling in invasive GAS
        map.insert("staphylococcus_aureus_daily_symptom_onset_probability".to_string(), 0.6);    // 60% per day - usually symptomatic
        map.insert("staphylococcus_epidermidis_daily_symptom_onset_probability".to_string(), 0.2); // 20% per day - device-associated pathogen often subacute
        map.insert("staphylococcus_epidermidis_symptom_onset_threshold_level".to_string(), 1.0);   // Needs higher burden for symptoms
        map.insert("staphylococcus_epidermidis_symptom_onset_delay_days".to_string(), 3.0);        // Slight delay before clinical detection
        map.insert("staphylococcus_epidermidis_base_bacteria_level_change".to_string(), 0.35);     // Slower growth kinetics than S. aureus
        map.insert("staphylococcus_epidermidis_max_level".to_string(), 4.0);                        // Lower peak burden due to biofilm focus
        map.insert("staphylococcus_epidermidis_environmental_acquisition_proportion".to_string(), 0.05); // Mostly device/skin origin
        map.insert("staphylococcus_epidermidis_microbiome_clearance_probability_per_day".to_string(), 0.015); // Chronic colonizer of skin/devices
        map.insert("staphylococcus_epidermidis_sepsis_baseline_log_odds".to_string(), -12.5);       // Rarely causes fulminant sepsis; keep orders of magnitude below high-virulence pathogens
        map.insert("staphylococcus_epidermidis_log_odds_sepsis_infection_level".to_string(), 0.04); // Slight level effect on sepsis risk
        map.insert("staphylococcus_epidermidis_log_odds_sepsis_infection_duration".to_string(), 0.005); // Chronic devices slowly accumulate risk
        map.insert("staphylococcus_epidermidis_non_sepsis_infection_death_log_odds".to_string(), -6.0); // Very low direct mortality

        map.insert("stenotrophomonas_maltophilia_daily_symptom_onset_probability".to_string(), 0.35); // 35% per day - clinically apparent in ventilated hosts
        map.insert("stenotrophomonas_maltophilia_symptom_onset_threshold_level".to_string(), 0.9);      // Moderate burden before symptoms
        map.insert("stenotrophomonas_maltophilia_symptom_onset_delay_days".to_string(), 2.5);          // Early signs once established
        map.insert("stenotrophomonas_maltophilia_base_bacteria_level_change".to_string(), 0.45);       // Moderate growth rate
        map.insert("stenotrophomonas_maltophilia_max_level".to_string(), 5.0);                          // Can reach high burdens in lungs
        map.insert("stenotrophomonas_maltophilia_environmental_acquisition_proportion".to_string(), 0.08); // Calibrated lower to curb runaway incidence
        map.insert("stenotrophomonas_maltophilia_microbiome_clearance_probability_per_day".to_string(), 0.06); // Persistent colonizer in ICU settings
        map.insert("stenotrophomonas_maltophilia_sepsis_baseline_log_odds".to_string(), -9.2);          // Opportunistic but still requires high burden or host compromise
        map.insert("stenotrophomonas_maltophilia_log_odds_sepsis_infection_level".to_string(), 0.08);   // Rising burden increases risk notably
        map.insert("stenotrophomonas_maltophilia_log_odds_sepsis_infection_duration".to_string(), 0.012); // Prolonged infection raises odds
        map.insert("stenotrophomonas_maltophilia_non_sepsis_infection_death_log_odds".to_string(), -4.0); // Some mortality via pneumonia progression

        // ENTERIC PATHOGENS - Moderate to high symptomatic rates
        map.insert("salmonella_enterica_serovar_typhi_daily_symptom_onset_probability".to_string(), 0.4);       // 40% per day
        map.insert("salmonella_enterica_serovar_typhi_base_bacteria_level_change".to_string(), 0.45);          // Longer incubation than typical enterics
        map.insert("salmonella_enterica_serovar_paratyphi_a_daily_symptom_onset_probability".to_string(), 0.4); // 40% per day
        map.insert("salmonella_enterica_serovar_paratyphi_a_base_bacteria_level_change".to_string(), 0.45);     // Similar incubation to typhi
        map.insert("shigella_spp._daily_symptom_onset_probability".to_string(), 0.6);                           // 60% per day
        map.insert("shigella_spp._base_bacteria_level_change".to_string(), 0.55);                               // Short incubation dysentery
        map.insert("vibrio_cholerae_daily_symptom_onset_probability".to_string(), 0.5);                         // 50% per day
        map.insert("vibrio_cholerae_base_bacteria_level_change".to_string(), 0.6);                              // Profuse cholera within 1-2 days
        map.insert("campylobacter_jejuni_daily_symptom_onset_probability".to_string(), 0.5);                    // 50% per day
        map.insert("campylobacter_jejuni_base_bacteria_level_change".to_string(), 0.52);                        // Incubation typically 2-4 days

        // CHRONIC/SLOW-ONSET PATHOGENS - Evidence-based symptom presentation rates
        // chlamydia_trachomatis - Most infections asymptomatic (~70-80% in women, ~50% in men)
        map.insert("chlamydia_trachomatis_daily_symptom_onset_probability".to_string(), 0.03);  // Only ~20-30% ever become symptomatic
        map.insert("chlamydia_trachomatis_base_bacteria_level_change".to_string(), 0.25);       // Slow intracellular replication
        map.insert("chlamydia_trachomatis_symptom_onset_threshold_level".to_string(), 0.8);     // Higher threshold before symptoms

        // treponema_pallidum - Syphilis has defined stages with variable presentation
        map.insert("treponema_pallidum_daily_symptom_onset_probability".to_string(), 0.08);     // Primary chancre develops in ~3-4 weeks
        map.insert("treponema_pallidum_base_bacteria_level_change".to_string(), 0.15);          // Very slow spirochete replication (33-hour doubling)
        map.insert("treponema_pallidum_symptom_onset_threshold_level".to_string(), 0.6);        // Moderate threshold

        // bordetella_pertussis - Catarrhal stage followed by paroxysmal cough
        map.insert("bordetella_pertussis_daily_symptom_onset_probability".to_string(), 0.35);   // 1-2 week incubation
        map.insert("bordetella_pertussis_base_bacteria_level_change".to_string(), 0.42);        // Moderate growth during catarrhal phase

        // helicobacter_pylori - Most infections (~80%) are asymptomatic
        map.insert("helicobacter_pylori_daily_symptom_onset_probability".to_string(), 0.005);   // Only ~20% develop symptomatic disease
        map.insert("helicobacter_pylori_symptom_onset_threshold_level".to_string(), 1.5);       // Very high threshold (chronic colonization)

        // MDR-TB - Slow progression; most latent infections never reactivate
        map.insert("mdr_mycobacterium_tuberculosis_daily_symptom_onset_probability".to_string(), 0.001); // ~5-10% lifetime reactivation risk
        map.insert("mdr_mycobacterium_tuberculosis_base_bacteria_level_change".to_string(), 0.08);       // Very slow mycobacterial growth
        map.insert("mdr_mycobacterium_tuberculosis_symptom_onset_threshold_level".to_string(), 2.0);     // High threshold for active disease

        // neisseria_gonorrhoeae - Variable symptoms (~10-20% asymptomatic in men, ~50% in women)
        map.insert("neisseria_gonorrhoeae_daily_symptom_onset_probability".to_string(), 0.25);  // Most symptomatic within 2-7 days
        map.insert("neisseria_gonorrhoeae_base_bacteria_level_change".to_string(), 0.55);       // Rapid mucosal colonization

        // yersinia_enterocolitica - Address intrinsic penicillin resistance
        // Reduce penicillins (intrinsic resistance)
        map.insert("drug_penicilling_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.02); // Intrinsic resistance
        map.insert("drug_ampicillin_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.02);  // Intrinsic resistance
        // Boost appropriate drugs modestly
        map.insert("drug_doxycycline_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.75); // Good activity
        map.insert("drug_ciprofloxacin_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.7); // Good activity
        map.insert("drug_trim_sulf_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.65);    // Good activity

        // streptococcus_pyogenes - Ensure penicillin remains preferred (no resistance ever develops)
        // S. pyogenes has never developed penicillin resistance - boost slightly to counter any drift
        map.insert("drug_penicilling_for_bacteria_streptococcus_pyogenes_potency_when_no_r".to_string(), 0.95); // Excellent and consistent activity

        // ENTERIC PATHOGENS - Modest fluoroquinolone boost for appropriate cases
        // Salmonella - boost fluoroquinolones for invasive disease (conservative increase)
        map.insert("drug_ciprofloxacin_for_bacteria_salmonella_enterica_serovar_typhi_potency_when_no_r".to_string(), 0.8);
        map.insert("drug_levofloxacin_for_bacteria_salmonella_enterica_serovar_typhi_potency_when_no_r".to_string(), 0.8);
        map.insert("drug_ciprofloxacin_for_bacteria_invasive_non-typhoidal_salmonella_spp._potency_when_no_r".to_string(), 0.75);

        // Shigella - boost fluoroquinolones modestly (first-line for severe cases)
        map.insert("drug_ciprofloxacin_for_bacteria_shigella_spp._potency_when_no_r".to_string(), 0.75);
        map.insert("drug_levofloxacin_for_bacteria_shigella_spp._potency_when_no_r".to_string(), 0.75);

        // Campylobacter - ensure fluoroquinolones are recognized as first-line
        map.insert("drug_ciprofloxacin_for_bacteria_campylobacter_jejuni_potency_when_no_r".to_string(), 0.8);
        map.insert("drug_levofloxacin_for_bacteria_campylobacter_jejuni_potency_when_no_r".to_string(), 0.8);

        // --- CLINICALLY APPROPRIATE DRUG-BACTERIA INITIATION MULTIPLIER OVERRIDES ---
        // Based on analysis of drug usage patterns, boost initiation probability for clinically appropriate combinations
        // Higher values = more likely to be selected as first-line therapy

        // ANTI-MRSA AGENTS FOR staphylococcus_aureus
        map.insert("drug_vancomycin_for_bacteria_staphylococcus_aureus_initiation_multiplier".to_string(), 5.0); // First-line for MRSA
        map.insert("drug_linezolid_for_bacteria_staphylococcus_aureus_initiation_multiplier".to_string(), 4.0); // Alternative for MRSA
        map.insert("drug_teicoplanin_for_bacteria_staphylococcus_aureus_initiation_multiplier".to_string(), 4.0); // Alternative for MRSA
        map.insert("drug_vancomycin_for_bacteria_staphylococcus_epidermidis_initiation_multiplier".to_string(), 6.0); // CoNS device infections rely on glycopeptides
        map.insert("drug_linezolid_for_bacteria_staphylococcus_epidermidis_initiation_multiplier".to_string(), 5.0); // Linezolid as IV/oral bridge
        map.insert("drug_teicoplanin_for_bacteria_staphylococcus_epidermidis_initiation_multiplier".to_string(), 5.0); // Teicoplanin effective for CoNS
        map.insert("drug_quinu_dalfo_for_bacteria_staphylococcus_epidermidis_initiation_multiplier".to_string(), 4.0); // Reserved for resistant CoNS
        map.insert("drug_trim_sulf_for_bacteria_staphylococcus_epidermidis_initiation_multiplier".to_string(), 2.5); // Occasionally used oral step-down

        map.insert("drug_penicilling_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.05); // Widespread beta-lactam resistance
        map.insert("drug_ampicillin_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.05);
        map.insert("drug_amoxicillin_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.05);
        map.insert("drug_cefazolin_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.15);
        map.insert("drug_ceftriaxone_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.1);
        map.insert("drug_gentamicin_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.2); // Synergy only
        map.insert("drug_vancomycin_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 1.05);
        map.insert("drug_linezolid_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.95);
        map.insert("drug_tedizolid_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.95);
        map.insert("drug_colistin_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.0);

        // ANTI-PSEUDOMONAL AGENTS FOR pseudomonas_aeruginosa
        map.insert("drug_piperacillin_tazobactam_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 5.0); // First-line anti-pseudomonal
        map.insert("drug_ceftazidime_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 4.5); // Good anti-pseudomonal activity
        map.insert("drug_cefepime_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 4.5); // Good anti-pseudomonal activity
        // CALIBRATION: Reduced carbapenem/reserve multipliers for stewardship - target <10% reserve drug usage
        map.insert("drug_meropenem_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.05); // Carbapenem - reserve agent, strong stewardship
        map.insert("drug_imipenem_c_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.05); // Carbapenem - reserve agent, strong stewardship
        map.insert("drug_colistin_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.05); // Last resort - very restricted

        // MACROLIDES FOR RESPIRATORY PATHOGENS
        map.insert("drug_erythromycin_for_bacteria_campylobacter_jejuni_initiation_multiplier".to_string(), 5.0); // First-line for Campylobacter
        map.insert("drug_azithromycin_for_bacteria_campylobacter_jejuni_initiation_multiplier".to_string(), 4.5); // Alternative for Campylobacter
        map.insert("drug_erythromycin_for_bacteria_haemophilus_influenzae_initiation_multiplier".to_string(), 4.0); // Good for H. influenzae
        map.insert("drug_azithromycin_for_bacteria_haemophilus_influenzae_initiation_multiplier".to_string(), 4.0); // Good for H. influenzae
        map.insert("drug_clarithromycin_for_bacteria_haemophilus_influenzae_initiation_multiplier".to_string(), 4.0); // Good for H. influenzae

        // stenotrophomonas_maltophilia - prefer TMP-SMX or minocycline, avoid carbapenems
        map.insert("drug_trim_sulf_for_bacteria_stenotrophomonas_maltophilia_initiation_multiplier".to_string(), 7.0); // TMP-SMX first line
        map.insert("drug_minocycline_for_bacteria_stenotrophomonas_maltophilia_initiation_multiplier".to_string(), 6.0); // Alternative therapy
        map.insert("drug_doxycycline_for_bacteria_stenotrophomonas_maltophilia_initiation_multiplier".to_string(), 4.5); // Doxy as related option
        map.insert("drug_levofloxacin_for_bacteria_stenotrophomonas_maltophilia_initiation_multiplier".to_string(), 3.5); // Fluoroquinolone rescue
        map.insert("drug_piperacillin_tazobactam_for_bacteria_stenotrophomonas_maltophilia_initiation_multiplier".to_string(), 0.05); // Intrinsic resistance
        map.insert("drug_ceftazidime_for_bacteria_stenotrophomonas_maltophilia_initiation_multiplier".to_string(), 0.1);
        map.insert("drug_meropenem_for_bacteria_stenotrophomonas_maltophilia_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_imipenem_c_for_bacteria_stenotrophomonas_maltophilia_initiation_multiplier".to_string(), 0.01);

        map.insert("drug_trim_sulf_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 1.05);
        map.insert("drug_minocycline_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.95);
        map.insert("drug_doxycycline_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.85);
        map.insert("drug_levofloxacin_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.6);
        map.insert("drug_ciprofloxacin_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.55);
        map.insert("drug_tobramycin_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.05);
        map.insert("drug_gentamicin_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.05);
        map.insert("drug_piperacillin_tazobactam_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.05);
        map.insert("drug_ceftazidime_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.15);
        map.insert("drug_meropenem_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.01);
        map.insert("drug_imipenem_c_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.05);
        map.insert("drug_colistin_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.0);
        map.insert("drug_vancomycin_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r".to_string(), 0.0);

        // ANTI-VRE AGENTS FOR ENTEROCOCCI
        map.insert("drug_vancomycin_for_bacteria_enterococcus_faecalis_initiation_multiplier".to_string(), 4.0); // First-line for Enterococcus
        map.insert("drug_vancomycin_for_bacteria_enterococcus_faecium_initiation_multiplier".to_string(), 3.5); // E. faecium more resistant
        // CALIBRATION: Linezolid is reserve agent - reduced for stewardship
        map.insert("drug_linezolid_for_bacteria_enterococcus_faecalis_initiation_multiplier".to_string(), 0.25); // Reserve - only for VRE
        map.insert("drug_linezolid_for_bacteria_enterococcus_faecium_initiation_multiplier".to_string(), 0.3); // Reserve - VRE E. faecium

        // C. DIFFICILE SPECIFIC AGENTS
        map.insert("drug_metronidazole_for_bacteria_clostridioides_difficile_initiation_multiplier".to_string(), 6.0); // First-line for C. diff
        map.insert("drug_vancomycin_for_bacteria_clostridioides_difficile_initiation_multiplier".to_string(), 5.0); // Oral vancomycin for severe C. diff

        // INTRACELLULAR PATHOGENS (tetracyclines, macrolides)
        map.insert("drug_doxycycline_for_bacteria_chlamydia_trachomatis_initiation_multiplier".to_string(), 5.0); // First-line for Chlamydia
        map.insert("drug_tetracycline_for_bacteria_chlamydia_trachomatis_initiation_multiplier".to_string(), 4.5); // Alternative for Chlamydia
        map.insert("drug_azithromycin_for_bacteria_chlamydia_trachomatis_initiation_multiplier".to_string(), 4.0); // Alternative for Chlamydia

        // CARBAPENEMS FOR ESBL PRODUCERS - reserve agents with strong stewardship
        // CALIBRATION: Reduced to target <10% reserve drug usage; carbapenems reserved for confirmed ESBL
        map.insert("drug_meropenem_for_bacteria_klebsiella_pneumoniae_initiation_multiplier".to_string(), 0.05); // Reserve - ESBL Klebsiella only
        map.insert("drug_imipenem_c_for_bacteria_klebsiella_pneumoniae_initiation_multiplier".to_string(), 0.05); // Reserve - ESBL Klebsiella only
        map.insert("drug_ertapenem_for_bacteria_klebsiella_pneumoniae_initiation_multiplier".to_string(), 0.05); // Reserve - outpatient ESBL option
        map.insert("drug_meropenem_for_bacteria_escherichia_coli_initiation_multiplier".to_string(), 0.05); // Reserve - ESBL E. coli only
        map.insert("drug_ertapenem_for_bacteria_escherichia_coli_initiation_multiplier".to_string(), 0.05); // Reserve - outpatient ESBL option

        // REDUCE INAPPROPRIATE COMBINATIONS
        // Penicillins should not be used for intrinsically resistant gram-negatives
        map.insert("drug_penicilling_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_penicilling_for_bacteria_acinetobacter_baumannii_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_ampicillin_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_ampicillin_for_bacteria_acinetobacter_baumannii_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_amoxicillin_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_amoxicillin_for_bacteria_acinetobacter_baumannii_initiation_multiplier".to_string(), 0.01);

        // Macrolides should not be used for Enterococcus (intrinsically resistant)
        map.insert("drug_erythromycin_for_bacteria_enterococcus_faecalis_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_erythromycin_for_bacteria_enterococcus_faecium_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_azithromycin_for_bacteria_enterococcus_faecalis_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_azithromycin_for_bacteria_enterococcus_faecium_initiation_multiplier".to_string(), 0.01);

        // --- TB-SPECIFIC DRUG POTENCIES ---
        map.insert("mdr_mycobacterium_tuberculosis_base_bacteria_level_change".to_string(), 0.15); // Very slow replication and chronic progression
        // MDR Mycobacterium tuberculosis has unique characteristics: thick cell wall, intracellular location, slow metabolism
        // Most standard antibiotics have poor activity against TB; only specific drugs are effective
        // MDR-TB represents established drug-resistant strains with guaranteed rifampicin resistance

        // FIRST-LINE TB DRUGS (high potency despite resistance)
        map.insert("drug_rifampicin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.6); // Primary first-line TB drug

        // SECOND-LINE TB DRUGS (moderate potency)
        map.insert("drug_levofloxacin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.4);  // Good FQ for TB
        map.insert("drug_moxifloxacin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.5);  // Best FQ for TB
        map.insert("drug_amikacin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.3);     // Injectable second-line
        map.insert("drug_linezolid_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.3);    // Oral second-line

        // OTHER FLUOROQUINOLONES (lower potency)
        map.insert("drug_ciprofloxacin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.3); // Less active than newer FQs
        map.insert("drug_ofloxacin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.35);    // Moderate activity

        // OTHER AMINOGLYCOSIDES (limited activity)
        map.insert("drug_gentamicin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.15);   // Poor TB activity
        map.insert("drug_tobramycin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.15);   // Poor TB activity

        // STANDARD ANTIBIOTICS (poor TB activity - thick cell wall barrier)
        map.insert("drug_penicilling_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.05);  // No TB activity
        map.insert("drug_ampicillin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.05);   // No TB activity
        map.insert("drug_vancomycin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.08);   // Minimal TB activity
        map.insert("drug_ceftriaxone_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.05);  // No TB activity
        map.insert("drug_meropenem_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.05);    // Minimal TB activity

        // RIFAMPICIN POTENCIES FOR OTHER BACTERIA (occasional use for severe staph infections)
        map.insert("drug_rifampicin_for_bacteria_staphylococcus_aureus_potency_when_no_r".to_string(), 0.4);         // Good anti-staph activity
        map.insert("drug_rifampicin_for_bacteria_enterococcus_faecalis_potency_when_no_r".to_string(), 0.2);         // Limited activity
        map.insert("drug_rifampicin_for_bacteria_enterococcus_faecium_potency_when_no_r".to_string(), 0.2);          // Limited activity
        // Most other bacteria: rifampicin has minimal activity (default 0.1 will apply)

        // --- TARGETED RESISTANCE EMERGENCE RATES---

        // NEVER-RESISTANT COMBINATIONS (should never develop resistance)
        // S. pyogenes has never developed penicillin resistance in >80 years of use
        map.insert("drug_penicilling_for_bacteria_streptococcus_pyogenes_resistance_emergence_rate_per_day_baseline".to_string(), 0.0);

        // VERY RARE RESISTANCE (extremely slow emergence)
        // Linezolid resistance in enterococci should remain very rare (~100x lower than baseline)
        map.insert("drug_linezolid_for_bacteria_enterococcus_faecium_resistance_emergence_rate_per_day_baseline".to_string(), 0.000003); 
        map.insert("drug_linezolid_for_bacteria_enterococcus_faecalis_resistance_emergence_rate_per_day_baseline".to_string(), 0.000003); 

        // PROBLEMATIC HIGH-RESISTANCE BACTERIA (higher emergence rates)
        // acinetobacter_baumannii - notorious for rapid resistance development
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_acinetobacter_baumannii_resistance_emergence_rate_per_day_baseline", drug), 0.01); 
        }

        // E. coli - moderate emergence rate (common pathogen with variable resistance)
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_escherichia_coli_resistance_emergence_rate_per_day_baseline", drug), 0.005); 
        }
        let high_pressure_e_coli_drugs = vec!["amoxicillin", "ampicillin", "ampicillin_sulbactam"];
        for &drug in high_pressure_e_coli_drugs.iter() {
            map.insert(format!("drug_{}_for_bacteria_escherichia_coli_resistance_emergence_rate_per_day_baseline", drug), 0.005);
        }

        // klebsiella_pneumoniae - rapid β-lactam resistance with selective retention of novel agents
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_klebsiella_pneumoniae_resistance_emergence_rate_per_day_baseline", drug), 0.0012);
        }
        let kleb_collapse_drugs = vec![
            "amoxicillin", "ampicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate",
            "piperacillin", "piperacillin_tazobactam", "ticarcillin", "ticarcillin_clavulanate",
            "cephalexin", "cephalothin", "cefazolin", "cefaclor", "cefuroxime"
        ];
        for &drug in kleb_collapse_drugs.iter() {
            map.insert(format!("drug_{}_for_bacteria_klebsiella_pneumoniae_resistance_emergence_rate_per_day_baseline", drug), 0.0025);
        }
        let kleb_preserved_drugs = vec!["ceftazidime_avibactam", "meropenem_vaborbactam", "colistin", "cefiderocol"];
        for &drug in kleb_preserved_drugs.iter() {
            map.insert(format!("drug_{}_for_bacteria_klebsiella_pneumoniae_resistance_emergence_rate_per_day_baseline", drug), 0.00018);
        }

        // pseudomonas_aeruginosa 
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_pseudomonas_aeruginosa_resistance_emergence_rate_per_day_baseline", drug), 0.035); 
        }
        let pseudo_reserve_drugs = vec!["colistin", "ceftazidime_avibactam", "meropenem_vaborbactam", "cefiderocol"];
        for &drug in pseudo_reserve_drugs.iter() {
            map.insert(format!("drug_{}_for_bacteria_pseudomonas_aeruginosa_resistance_emergence_rate_per_day_baseline", drug), 0.012);
        }
        let pseudo_problem_drugs = vec!["levofloxacin", "ciprofloxacin", "meropenem", "imipenem", "ceftazidime", "cefepime"];
        for &drug in pseudo_problem_drugs.iter() {
            map.insert(format!("drug_{}_for_bacteria_pseudomonas_aeruginosa_resistance_emergence_rate_per_day_baseline", drug), 0.065);
        }

        // stenotrophomonas_maltophilia gs
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_stenotrophomonas_maltophilia_resistance_emergence_rate_per_day_baseline", drug), 0.0009);         
        }

        // staphylococcus_epidermidis - biofilm former, moderate resistance development
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_staphylococcus_epidermidis_resistance_emergence_rate_per_day_baseline", drug), 0.0006); 
        }

        // proteus_spp.
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_proteus_spp._resistance_emergence_rate_per_day_baseline", drug), 0.2); 
        }

        // pseudomonas_aeruginosa
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_pseudomonas_aeruginosa_resistance_emergence_rate_per_day_baseline", drug), 0.005); 
        }
        
        // salmonella_enterica_serovar_paratyphi_a
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_salmonella_enterica_serovar_paratyphi_a_resistance_emergence_rate_per_day_baseline", drug), 0.1); 
        }
     
        // salmonella_enterica_serovar_typhi
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_salmonella_enterica_serovar_typhi_resistance_emergence_rate_per_day_baseline", drug), 0.005); 
        }
        
       
        // shigella_spp_
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_shigella_spp_resistance_emergence_rate_per_day_baseline", drug), 0.2); 
        }
        
      
        // staphylococcus_aureus
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_staphylococcus_aureus_resistance_emergence_rate_per_day_baseline", drug), 0.8); 
        }
        
        // bordetella pertussis
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_bordetella_pertussis_resistance_emergence_rate_per_day_baseline", drug), 0.2); 
        }

        // campylobacter_jejuni
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_campylobacter_jejuni_resistance_emergence_rate_per_day_baseline", drug), 0.005); 
        }

        // chlamydia_trachomatis
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_chlamydia_trachomatis_resistance_emergence_rate_per_day_baseline", drug), 0.05); 
        }

        // citrobacter_spp._
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_citrobacter_spp._resistance_emergence_rate_per_day_baseline", drug), 0.5); 
        }

        // clostridioides_difficile
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_clostridioides_difficile_resistance_emergence_rate_per_day_baseline", drug), 0.3); 
        }

        // enterobacter_cloacae
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_enterobacter_cloacae_resistance_emergence_rate_per_day_baseline", drug), 0.2); 
        }
        
        // enterobacter_spp.
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_enterobacter_spp._resistance_emergence_rate_per_day_baseline", drug), 0.5); 
        }

        // enterococcus_faecalis
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_enterococcus_faecalis_resistance_emergence_rate_per_day_baseline", drug), 0.3); 
        }

        // enterococcus_faecium
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_enterococcus_faecium_resistance_emergence_rate_per_day_baseline", drug), 0.3); 
        }

        // haemophilus_influenza
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_haemophilus_influenza_resistance_emergence_rate_per_day_baseline", drug), 0.005); 
        }

        // helicobacter_pylori
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_helicobacter_pylori_resistance_emergence_rate_per_day_baseline", drug), 0.1); 
        }

        // listeria_monocytogenes
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_listeria_monocytogenes_resistance_emergence_rate_per_day_baseline", drug), 0.2); 
        }

        // morganella_spp.
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_morganella_spp._resistance_emergence_rate_per_day_baseline", drug), 0.5); 
        }        

        // moraxella_catarrhalis
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_moraxella_catarrhalis_resistance_emergence_rate_per_day_baseline", drug), 0.03); 
        }        
        
        // serratia_spp.
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_serratia_spp._resistance_emergence_rate_per_day_baseline", drug), 0.8); 
        }
       
        // streptococcus_agalactiae
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_streptococcus_agalactiae_resistance_emergence_rate_per_day_baseline", drug), 0.02); 
        }

        // treponema_pallidum
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_treponema_pallidum_resistance_emergence_rate_per_day_baseline", drug), 0.03); 
        }

        // vibrio_cholerae
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_vibrio_cholerae_resistance_emergence_rate_per_day_baseline", drug), 0.3); 
        }

        // yersinia_enterocolitica
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_yersinia_enterocolitica_resistance_emergence_rate_per_day_baseline", drug), 0.8); 
        }

        // neisseria_meningitidis
        for &drug in DRUG_SHORT_NAMES.iter() {
        map.insert(format!("drug_{}_for_bacteria_neisseria_meningitidis_resistance_emergence_rate_per_day_baseline", drug), 0.2); 
        }


/*

        // here here

*/
        

        // SPECIFIC DRUG-BACTERIA COMBINATIONS WITH CLINICAL CONSTRAINTS
        // Colistin resistance should remain rare (last-resort antibiotic)
        let gram_negative_bacteria = vec![
            "acinetobacter_baumannii", "pseudomonas_aeruginosa", "escherichia_coli",
            "klebsiella_pneumoniae", "enterobacter_spp.", "citrobacter_spp.",
            "serratia_spp.", "proteus_spp.", "morganella_spp.", "enterobacter_cloacae"
        ];
        for &bacteria in gram_negative_bacteria.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // CALIBRATION: Colistin resistance should remain rare (~10x lower than baseline)
                map.insert(format!("drug_colistin_for_bacteria_{}_resistance_emergence_rate_per_day_baseline", bacteria), 0.00003); 
            }
        }

        // Nitrofurantoin & amoxicillin/clavulanate resistance in E. coli should remain low (important for UTI treatment)
        map.insert("drug_nitrofurantoin_for_bacteria_escherichia_coli_resistance_emergence_rate_per_day_baseline".to_string(), 0.00003);
        map.insert("drug_amoxicillin_clavulanate_for_bacteria_escherichia_coli_resistance_emergence_rate_per_day_baseline".to_string(), 0.00003);

        // enterococcus_faecium - glycopeptide, oxazolidinone, and lipopeptide resistance emerges slowly but can persist
        map.insert("drug_vancomycin_for_bacteria_enterococcus_faecium_resistance_emergence_rate_per_day_baseline".to_string(), 0.0008);
        map.insert("drug_teicoplanin_for_bacteria_enterococcus_faecium_resistance_emergence_rate_per_day_baseline".to_string(), 0.00065);
        map.insert("drug_daptomycin_for_bacteria_enterococcus_faecium_resistance_emergence_rate_per_day_baseline".to_string(), 0.0005);

        // neisseria_gonorrhoeae - limit emergence to maintain dual-therapy effectiveness
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_for_bacteria_neisseria_gonorrhoeae_resistance_emergence_rate_per_day_baseline", drug), 0.005);
        }
        map.insert("drug_ceftriaxone_for_bacteria_neisseria_gonorrhoeae_resistance_emergence_rate_per_day_baseline".to_string(), 0.005);
        let gonorrhea_fq = vec!["ciprofloxacin", "levofloxacin", "ofloxacin", "moxifloxacin"];
        for &drug in gonorrhea_fq.iter() {
            map.insert(format!("drug_{}_for_bacteria_neisseria_gonorrhoeae_resistance_emergence_rate_per_day_baseline", drug), 0.005);
        }

        // neisseria_meningitidis - penicillin and ceftriaxone resistance remains extremely rare
        let meningo_anchor_drugs = vec!["penicilling", "ampicillin", "ceftriaxone"];
        for &drug in meningo_anchor_drugs.iter() {
            map.insert(format!("drug_{}_for_bacteria_neisseria_meningitidis_resistance_emergence_rate_per_day_baseline", drug), 0.000005);
        }

        // moraxella_catarrhalis - usually susceptible; best keep emergence low for beta-lactams/macrolides
        let moraxella_beta_lactams = vec!["amoxicillin", "amoxicillin_clavulanate", "ampicillin", "ampicillin_sulbactam", "piperacillin_tazobactam", "cefuroxime"];
        for &drug in moraxella_beta_lactams.iter() {
            map.insert(format!("drug_{}_for_bacteria_moraxella_catarrhalis_resistance_emergence_rate_per_day_baseline", drug), 0.005);
        }
        let moraxella_macrolides = vec!["azithromycin", "clarithromycin", "erythromycin"];
        for &drug in moraxella_macrolides.iter() {
            map.insert(format!("drug_{}_for_bacteria_moraxella_catarrhalis_resistance_emergence_rate_per_day_baseline", drug), 0.005);
        }
        let moraxella_reserve = vec!["ceftriaxone", "ceftazidime", "levofloxacin", "ciprofloxacin"];
        for &drug in moraxella_reserve.iter() {
            map.insert(format!("drug_{}_for_bacteria_moraxella_catarrhalis_resistance_emergence_rate_per_day_baseline", drug), 0.005);
        }

        // Vancomycin resistance should be impossible in Gram-negative bacteria (intrinsic resistance handled by potency)
        for &bacteria in gram_negative_bacteria.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                map.insert(format!("drug_vancomycin_for_bacteria_{}_resistance_emergence_rate_per_day_baseline", bacteria), 0.0);
                map.insert(format!("drug_teicoplanin_for_bacteria_{}_resistance_emergence_rate_per_day_baseline", bacteria), 0.0);
                // Set potency to zero - glycopeptides have no activity against gram-negative bacteria
                map.insert(format!("drug_vancomycin_for_bacteria_{}_potency_when_no_r", bacteria), 0.0);
                map.insert(format!("drug_teicoplanin_for_bacteria_{}_potency_when_no_r", bacteria), 0.0);
                // Linezolid and daptomycin also have minimal/no activity against gram-negatives
                map.insert(format!("drug_linezolid_for_bacteria_{}_potency_when_no_r", bacteria), 0.0);
                map.insert(format!("drug_daptomycin_for_bacteria_{}_potency_when_no_r", bacteria), 0.0);
            }
        }

        // Colistin has no activity against gram-positive bacteria (intrinsic resistance)
        let gram_positive_bacteria = vec![
            "staphylococcus_aureus", "streptococcus_pneumoniae", "streptococcus_pyogenes",
            "streptococcus_agalactiae", "enterococcus_faecalis", "enterococcus_faecium",
            "staphylococcus_epidermidis", "streptococcus_viridans"
        ];
        for &bacteria in gram_positive_bacteria.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                map.insert(format!("drug_colistin_for_bacteria_{}_potency_when_no_r", bacteria), 0.0);
            }
        }


        // for each drug-bacteria combination will need a specific multiplier for initiation rate
        // will need changes also in mod.rs

        // --- Evidence-Based Bacteria-Specific Drug Cessation Rates ---
        // Based on typical treatment durations and patient adherence patterns
        // Target: 90% completion for appropriate treatment courses
        //
        // Literature Sources for Regional Variations:
        // - Butler et al. Patient adherence to antibiotic therapy. BMC Infect Dis. 2007;7:72
        // - Kardas et al. A systematic review of adherence with medications for chronic conditions. Patient Prefer Adherence. 2013;7:479-489
        // - WHO Global Health Observatory data on healthcare access and quality
        // - Sabaté E. Adherence to long-term therapies: evidence for action. WHO. 2003
        // - Regional healthcare infrastructure studies (World Bank, OECD Health Statistics)

        // REGIONAL CESSATION MULTIPLIERS (applied to base bacteria-specific rates)
        // High-income regions with strong healthcare systems and universal access
        map.insert("north_america_cessation_multiplier".to_string(), 0.85); // Good healthcare access, medication coverage
        map.insert("europe_cessation_multiplier".to_string(), 0.80); // Universal healthcare, excellent adherence programs
        map.insert("oceania_cessation_multiplier".to_string(), 0.85); // Similar to North America (Australia/NZ dominant)

        // Middle-income regions with variable healthcare access
        map.insert("asia_cessation_multiplier".to_string(), 1.15); // Mixed development levels, variable infrastructure
        map.insert("south_america_cessation_multiplier".to_string(), 1.25); // Economic constraints, healthcare gaps

        // Lower-income regions with significant healthcare challenges
        map.insert("africa_cessation_multiplier".to_string(), 1.40); // Economic barriers, limited infrastructure

        // TB-SPECIFIC REGIONAL ADHERENCE MODIFIERS (applied in addition to general regional multipliers)
        // These capture DOT (Directly Observed Therapy) program effectiveness and TB-specific healthcare infrastructure
        // Applied as: final_tb_cessation_rate = base_tb_rate × regional_multiplier × tb_adherence_modifier
        //
        // WHY TB HAS UNIQUE ADHERENCE MODIFIERS (and other bacteria don't):
        // 1. DOT programs: TB has specific Directly Observed Therapy programs that don't exist for other infections
        // 2. Long treatment duration: 6-24 months vs. days/weeks for other infections creates unique challenges
        // 3. Public health priority: TB control programs are disease-specific with dedicated infrastructure
        // 4. Social support systems: Contact tracing, case management, nutritional support specific to TB
        // 5. Regulatory framework: International TB control standards (WHO, CDC) create region-specific adherence differences
        // Other bacteria use general regional cessation multipliers but don't have disease-specific adherence programs
        map.insert("mdr_mycobacterium_tuberculosis_north_america_adherence_modifier".to_string(), 0.4); // Excellent TB programs → 60% better adherence
        map.insert("mdr_mycobacterium_tuberculosis_europe_adherence_modifier".to_string(), 0.3);        // Best TB programs globally → 70% better adherence
        map.insert("mdr_mycobacterium_tuberculosis_oceania_adherence_modifier".to_string(), 0.4);       // Similar to North America
        map.insert("mdr_mycobacterium_tuberculosis_asia_adherence_modifier".to_string(), 0.6);          // Good but variable DOT programs → 40% better adherence
        map.insert("mdr_mycobacterium_tuberculosis_south_america_adherence_modifier".to_string(), 0.7); // Moderate DOT programs → 30% better adherence
        map.insert("mdr_mycobacterium_tuberculosis_africa_adherence_modifier".to_string(), 0.8);        // Resource constraints limit DOT effectiveness → 20% better adherence

        // DEFAULT CESSATION RATES (for bacteria without specific overrides)
    map.insert("random_drug_cessation_probability".to_string(), 0.0045); // 0.45% daily = ~94% complete 14-day course
    map.insert("random_drug_cessation_probability_if_no_active_infection".to_string(), 0.15); // Higher probability if no active infection

        // SHORT-COURSE INFECTIONS (3-7 days typical treatment)
        // UTI, simple pneumonia, skin infections - target 90% completion for 3-7 day courses
        // Daily cessation rate: 0.035 = 90% complete 3-day course, 0.015 = 90% complete 7-day course
        map.insert("escherichia_coli_drug_cessation_probability".to_string(), 0.025); // UTI: 3-5 days, compromise rate
        map.insert("streptococcus_pneumoniae_drug_cessation_probability".to_string(), 0.015); // Pneumonia: 5-7 days
        map.insert("staphylococcus_aureus_drug_cessation_probability".to_string(), 0.015); // Skin/soft tissue: 7-10 days
        map.insert("streptococcus_pyogenes_drug_cessation_probability".to_string(), 0.015); // Strep throat: 10 days
        map.insert("haemophilus_influenzae_drug_cessation_probability".to_string(), 0.015); // Respiratory: 7-10 days

        // MODERATE-COURSE INFECTIONS (10-14 days typical treatment)
        map.insert("klebsiella_pneumoniae_drug_cessation_probability".to_string(), 0.0075); // Hospital pneumonia: 10-14 days
        map.insert("pseudomonas_aeruginosa_drug_cessation_probability".to_string(), 0.0075); // Complex infections: 14-21 days
        map.insert("acinetobacter_baumannii_drug_cessation_probability".to_string(), 0.0075); // Hospital infections: 14-21 days
        map.insert("enterococcus_faecium_drug_cessation_probability".to_string(), 0.0075); // VRE infections: 10-14 days
        map.insert("enterococcus_faecalis_drug_cessation_probability".to_string(), 0.0075); // Enterococcal infections: 10-14 days

        // EXTENDED-COURSE INFECTIONS (21+ days typical treatment)
        map.insert("clostridioides_difficile_drug_cessation_probability".to_string(), 0.005); // C. diff: 10-14 days + taper
        map.insert("clostridioides_difficile_base_bacteria_level_change".to_string(), 0.55); // Rapid toxin surge after microbiome disruption
        map.insert("helicobacter_pylori_drug_cessation_probability".to_string(), 0.005); // Triple therapy: 14 days + follow-up

        // CHRONIC/PERSISTENT INFECTIONS (weeks to months)
        // Very low cessation rates for infections requiring prolonged treatment
        map.insert("mdr_mycobacterium_tuberculosis_drug_cessation_probability".to_string(), 0.0006); // MDR-TB: 6-24 months (0.06% daily = 90% complete 180 days)
        map.insert("chlamydia_trachomatis_drug_cessation_probability".to_string(), 0.007); // Chlamydia: 7-21 days depending on regimen

        // ENTERIC PATHOGENS (variable courses 3-14 days)
        map.insert("salmonella_enteritidis_drug_cessation_probability".to_string(), 0.0105); // Salmonella: 7-14 days
        map.insert("shigella_spp_drug_cessation_probability".to_string(), 0.0105); // Shigellosis: 3-5 days
        map.insert("campylobacter_jejuni_drug_cessation_probability".to_string(), 0.015); // Campylobacter: 5-7 days
        map.insert("vibrio_cholerae_drug_cessation_probability".to_string(), 0.025); // Cholera: 3 days

        // RESPIRATORY PATHOGENS
        map.insert("bordetella_pertussis_drug_cessation_probability".to_string(), 0.0075); // Pertussis: 14 days
        map.insert("neisseria_meningitidis_drug_cessation_probability".to_string(), 0.01); // Meningitis: 7-10 days
        map.insert("streptococcus_agalactiae_drug_cessation_probability".to_string(), 0.015); // GBS: 7-10 days

        // General Acquisition & Resistance Parameters
        // === [F] Regional acquisition pressure baselines ===
        // Sets consistent fallbacks for infection chance modifiers and then layers region-specific
        // multipliers reflecting surveillance data. Override individual entries when calibrating
        // against new incidence estimates.
        // --- Logistic Model Parameters for Infection and Microbiome Acquisition ---
        // Infection acquisition (site infection)
        map.insert("acquisition_log_odds_baseline".to_string(), -17.0); 
        // This gives ~0.000005% per day per bacteria = ~0.018% per year per bacteria
        // With 34 bacteria: ~0.6% annual baseline, realistic after regional/risk adjustments
        map.insert("neisseria_meningitidis_acquisition_log_odds_baseline".to_string(), -16.5); 
        map.insert("haemophilus_influenzae_acquisition_log_odds_baseline".to_string(), -16.5); 
        map.insert("salmonella_enterica_serovar_typhi_acquisition_log_odds_baseline".to_string(), -17.0); 
        map.insert("bordetella_pertussis_acquisition_log_odds_baseline".to_string(), -13.3);
        map.insert("acinetobacter_baumannii_acquisition_log_odds_baseline".to_string(), -15.5); 
        map.insert("campylobacter_jejuni_acquisition_log_odds_baseline".to_string(), -13.0); 
        map.insert("chlamydia_trachomatis_acquisition_log_odds_baseline".to_string(), -13.2); 
        map.insert("citrobacter_spp._acquisition_log_odds_baseline".to_string(), -16.0);
        map.insert("clostridioides_difficile_acquisition_log_odds_baseline".to_string(), -15.2); 
        map.insert("enterobacter_cloacae_acquisition_log_odds_baseline".to_string(), -16.0); 
        map.insert("enterobacter_spp._acquisition_log_odds_baseline".to_string(), -17.0); 
        map.insert("enterococcus_faecalis_acquisition_log_odds_baseline".to_string(), -16.5); 
        map.insert("enterococcus_faecium_acquisition_log_odds_baseline".to_string(), -17.0); 
        map.insert("escherichia_coli_acquisition_log_odds_baseline".to_string(), -13.5); 
        map.insert("helicobacter_pylori_acquisition_log_odds_baseline".to_string(), -13.5); 
        map.insert("invasive_non-typhoidal_salmonella_spp._acquisition_log_odds_baseline".to_string(), -16.5); 
        map.insert("klebsiella_pneumoniae_acquisition_log_odds_baseline".to_string(), -16.0); 
        map.insert("listeria_monocytogenes_acquisition_log_odds_baseline".to_string(), -15.5); 
        map.insert("mdr_mycobacterium_tuberculosis_acquisition_log_odds_baseline".to_string(), -16.0); 
        map.insert("moraxella_catarrhalis_acquisition_log_odds_baseline".to_string(), -16.5); 
        map.insert("morganella_spp._acquisition_log_odds_baseline".to_string(), -15.5); 
        map.insert("neisseria_gonorrhoeae_acquisition_log_odds_baseline".to_string(), -13.5); 
        map.insert("proteus_spp._acquisition_log_odds_baseline".to_string(), -15.5); 
        map.insert("pseudomonas_aeruginosa_acquisition_log_odds_baseline".to_string(), -15.0);
        map.insert("salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds_baseline".to_string(), -16.0); 
        map.insert("serratia_spp._acquisition_log_odds_baseline".to_string(), -17.0); 
        map.insert("shigella_spp._acquisition_log_odds_baseline".to_string(), -13.5); 
        map.insert("staphylococcus_epidermidis_acquisition_log_odds_baseline".to_string(), -16.0); 
        map.insert("stenotrophomonas_maltophilia_acquisition_log_odds_baseline".to_string(), -17.0); 
        map.insert("staphylococcus_aureus_acquisition_log_odds_baseline".to_string(), -15.0); 
        map.insert("streptococcus_agalactiae_acquisition_log_odds_baseline".to_string(), -17.0); 
        map.insert("streptococcus_pneumoniae_acquisition_log_odds_baseline".to_string(), -13.0); 
        map.insert("streptococcus_pyogenes_acquisition_log_odds_baseline".to_string(), -16.5); 
        map.insert("treponema_pallidum_acquisition_log_odds_baseline".to_string(), -15.0); 
        map.insert("vibrio_cholerae_acquisition_log_odds_baseline".to_string(), -18.0); 
        map.insert("yersinia_enterocolitica_acquisition_log_odds_baseline".to_string(), -17.0); 
        map.insert("log_odds_vaccinated".to_string(), -2.0); // Vaccination reduces log-odds
        map.insert("log_odds_microbiome_present".to_string(), 0.5); // Microbiome presence effect (example)
        map.insert("log_odds_hospital_acquired".to_string(), 2.0); // Hospital-acquired effect (default/fallback)

        // Bacteria-specific hospital acquisition log-odds (healthcare-associated infection risk)
        // These parameters reflect clinical reality where certain bacteria have much higher
        // acquisition risk in hospital settings due to:
        // - Environmental persistence (Acinetobacter, C. diff)
        // - Device-associated transmission (Pseudomonas, Klebsiella)
        // - Healthcare worker transmission (MRSA, VRE)
        // - Antibiotic selection pressure (C. diff, ESBL producers)
        // Values are log-odds multipliers: exp(3.0) = 20x higher risk, exp(2.0) = 7x higher, etc.

        // High HAI risk bacteria (major hospital pathogens)
        map.insert("acinetobacter_baumannii_log_odds_hospital_acquired".to_string(), 3.4); // ~30x higher risk (exp(3.4))
        map.insert("pseudomonas_aeruginosa_log_odds_hospital_acquired".to_string(), 3.0); // 20x higher risk
        map.insert("enterococcus_faecium_log_odds_hospital_acquired".to_string(), 3.3); // 27x higher risk (VRE)
        map.insert("staphylococcus_aureus_log_odds_hospital_acquired".to_string(), 2.3); // 10x higher risk (MRSA)
        map.insert("clostridioides_difficile_log_odds_hospital_acquired".to_string(), 1.4); // Still elevated but scaled back to curb over-incidence
        map.insert("klebsiella_pneumoniae_log_odds_hospital_acquired".to_string(), 2.0); // 7x higher risk
        map.insert("enterobacter_spp._log_odds_hospital_acquired".to_string(), 0.8); // Moderately elevated after calibration
        map.insert("enterobacter_cloacae_log_odds_hospital_acquired".to_string(), 0.9); // Moderately elevated after calibration
        map.insert("serratia_spp._log_odds_hospital_acquired".to_string(), 2.0); // 7x higher risk
        map.insert("citrobacter_spp._log_odds_hospital_acquired".to_string(), 0.4); // Trimmed to reduce excess incidence

        // Moderate HAI risk bacteria
        map.insert("escherichia_coli_log_odds_hospital_acquired".to_string(), 2.0); // 7x higher risk (device-associated)
        map.insert("enterococcus_faecalis_log_odds_hospital_acquired".to_string(), 2.1); // 8x higher risk
        map.insert("streptococcus_pneumoniae_log_odds_hospital_acquired".to_string(), 0.5); // Moderately above community baseline after rollback
        map.insert("proteus_spp._log_odds_hospital_acquired".to_string(), 1.4); // 4x higher risk
        map.insert("morganella_spp._log_odds_hospital_acquired".to_string(), -0.2); // Treat predominately community cases after calibration
        map.insert("listeria_monocytogenes_base_bacteria_level_change".to_string(), 0.25); // Prolonged incubation (up to weeks)
        map.insert("listeria_monocytogenes_log_odds_hospital_acquired".to_string(), -0.4); // Foodborne pathway dominates
        map.insert("neisseria_meningitidis_log_odds_hospital_acquired".to_string(), 0.3); // Lower nosocomial amplification after vaccine scale-up
        map.insert("streptococcus_pyogenes_log_odds_hospital_acquired".to_string(), 1.1); // 3x higher risk
        map.insert("streptococcus_agalactiae_log_odds_hospital_acquired".to_string(), 1.2); // 3.3x higher risk
        map.insert("haemophilus_influenzae_log_odds_hospital_acquired".to_string(), 0.9); // 2.5x higher risk
        map.insert("moraxella_catarrhalis_log_odds_hospital_acquired".to_string(), 0.3); // Slightly elevated nosocomial risk after rollback
        map.insert("yersinia_enterocolitica_log_odds_hospital_acquired".to_string(), 0.5); // 1.6x higher risk

        // Low/No HAI risk bacteria (mostly community pathogens)
        map.insert("chlamydia_trachomatis_log_odds_hospital_acquired".to_string(), -1.0); // 0.37x (lower risk in hospital)
        map.insert("neisseria_gonorrhoeae_log_odds_hospital_acquired".to_string(), -0.8); // 0.45x (lower risk in hospital)
        map.insert("treponema_pallidum_log_odds_hospital_acquired".to_string(), -1.2); // 0.3x (much lower risk)
        map.insert("salmonella_enterica_serovar_typhi_log_odds_hospital_acquired".to_string(), 0.0); // Neutral (travel-related)
        map.insert("salmonella_enterica_serovar_paratyphi_a_log_odds_hospital_acquired".to_string(), 0.0); // Neutral
        map.insert("invasive_non-typhoidal_salmonella_spp._log_odds_hospital_acquired".to_string(), 0.2); // Slight increase
        map.insert("shigella_spp._log_odds_hospital_acquired".to_string(), -0.3); // 0.74x (foodborne, lower in hospital)
        map.insert("vibrio_cholerae_log_odds_hospital_acquired".to_string(), -0.5); // 0.6x (waterborne, lower in hospital)
        map.insert("campylobacter_jejuni_log_odds_hospital_acquired".to_string(), -0.6); // Further reduced nosocomial amplification
        map.insert("mdr_mycobacterium_tuberculosis_log_odds_hospital_acquired".to_string(), 0.5); // Slightly elevated nosocomial risk

        // effect of region on bacteria acquisition risk (vs north america)
        map.insert("south_america_shigella_spp_acquisition_log_odds".to_string(), 0.4);
        map.insert("africa_shigella_spp_acquisition_log_odds".to_string(), 3.0);
        map.insert("europe_shigella_spp_acquisition_log_odds".to_string(), 0.5);
        map.insert("asia_shigella_spp_acquisition_log_odds".to_string(), 2.0);
        map.insert("oceania_shigella_spp_acquisition_log_odds".to_string(), 0.7);

        map.insert("africa_acinetobacter_baumannii_acquisition_log_odds".to_string(), 2.3);
        map.insert("europe_acinetobacter_baumannii_acquisition_log_odds".to_string(), -0.3);
        map.insert("asia_acinetobacter_baumannii_acquisition_log_odds".to_string(), 2.0);
        map.insert("south_america_acinetobacter_baumannii_acquisition_log_odds".to_string(), 1.6);
        map.insert("oceania_acinetobacter_baumannii_acquisition_log_odds".to_string(), -0.1);

        // citrobacter_spp. - Predominantly healthcare-associated, modest regional differences
        map.insert("africa_citrobacter_spp._acquisition_log_odds".to_string(), 0.1);
        map.insert("europe_citrobacter_spp._acquisition_log_odds".to_string(), -0.5);
        map.insert("asia_citrobacter_spp._acquisition_log_odds".to_string(), -0.2);
        map.insert("south_america_citrobacter_spp._acquisition_log_odds".to_string(), -0.3);
        map.insert("oceania_citrobacter_spp._acquisition_log_odds".to_string(), -0.4);

        // enterobacter_spp. - Predominantly healthcare-associated, modest regional differences
        map.insert("africa_enterobacter_spp._acquisition_log_odds".to_string(), 0.2);
        map.insert("europe_enterobacter_spp._acquisition_log_odds".to_string(), -0.5);
        map.insert("asia_enterobacter_spp._acquisition_log_odds".to_string(), 0.0);
        map.insert("south_america_enterobacter_spp._acquisition_log_odds".to_string(), -0.2);
        map.insert("oceania_enterobacter_spp._acquisition_log_odds".to_string(), -0.3);

        // enterococcus_faecalis - Mixed healthcare/community, moderate regional differences
        map.insert("africa_enterococcus_faecalis_acquisition_log_odds".to_string(), 1.6);
        map.insert("europe_enterococcus_faecalis_acquisition_log_odds".to_string(), 0.2);
        map.insert("asia_enterococcus_faecalis_acquisition_log_odds".to_string(), 1.3);
        map.insert("south_america_enterococcus_faecalis_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_enterococcus_faecalis_acquisition_log_odds".to_string(), 0.2);

        // enterococcus_faecium - Predominantly healthcare-associated, high AMR burden
        map.insert("africa_enterococcus_faecium_acquisition_log_odds".to_string(), 1.4);
        map.insert("europe_enterococcus_faecium_acquisition_log_odds".to_string(), 0.2);
        map.insert("asia_enterococcus_faecium_acquisition_log_odds".to_string(), 1.2);
        map.insert("south_america_enterococcus_faecium_acquisition_log_odds".to_string(), 0.9);
        map.insert("oceania_enterococcus_faecium_acquisition_log_odds".to_string(), 0.2);

        // escherichia_coli - Major community and healthcare pathogen, high regional variation
        map.insert("africa_escherichia_coli_acquisition_log_odds".to_string(), 2.4);
        map.insert("europe_escherichia_coli_acquisition_log_odds".to_string(), 0.4);
        map.insert("asia_escherichia_coli_acquisition_log_odds".to_string(), 2.1);
        map.insert("south_america_escherichia_coli_acquisition_log_odds".to_string(), 1.6);
        map.insert("oceania_escherichia_coli_acquisition_log_odds".to_string(), 0.7);

        // klebsiella_pneumoniae - Mixed community/healthcare, major AMR threat
        map.insert("africa_klebsiella_pneumoniae_acquisition_log_odds".to_string(), 1.6);
        map.insert("europe_klebsiella_pneumoniae_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_klebsiella_pneumoniae_acquisition_log_odds".to_string(), 1.3);
        map.insert("south_america_klebsiella_pneumoniae_acquisition_log_odds".to_string(), 0.8);
        map.insert("oceania_klebsiella_pneumoniae_acquisition_log_odds".to_string(), 0.0);

        // morganella_spp. - Predominantly healthcare-associated, urinary tract infections
        map.insert("africa_morganella_spp._acquisition_log_odds".to_string(), -0.3);
        map.insert("europe_morganella_spp._acquisition_log_odds".to_string(), -0.7);
        map.insert("asia_morganella_spp._acquisition_log_odds".to_string(), -0.4);
        map.insert("south_america_morganella_spp._acquisition_log_odds".to_string(), -0.5);
        map.insert("oceania_morganella_spp._acquisition_log_odds".to_string(), -0.6);

        // proteus_spp. - Mixed healthcare/community, urinary tract and wound infections
        map.insert("africa_proteus_spp._acquisition_log_odds".to_string(), 1.1);
        map.insert("europe_proteus_spp._acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_proteus_spp._acquisition_log_odds".to_string(), 0.8);
        map.insert("south_america_proteus_spp._acquisition_log_odds".to_string(), 0.5);
        map.insert("oceania_proteus_spp._acquisition_log_odds".to_string(), 0.0);

        // serratia_spp. - Predominantly healthcare-associated, opportunistic pathogen
        map.insert("africa_serratia_spp._acquisition_log_odds".to_string(), 0.8);
        map.insert("europe_serratia_spp._acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_serratia_spp._acquisition_log_odds".to_string(), 0.6);
        map.insert("south_america_serratia_spp._acquisition_log_odds".to_string(), 0.4);
        map.insert("oceania_serratia_spp._acquisition_log_odds".to_string(), 0.0);

        // pseudomonas_aeruginosa - Predominantly healthcare-associated, major AMR threat
        map.insert("africa_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 1.6);
        map.insert("europe_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 0.3);
        map.insert("asia_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 1.3);
        map.insert("south_america_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 0.4);

        // staphylococcus_aureus - Major community and healthcare pathogen, high regional variation
        map.insert("africa_staphylococcus_aureus_acquisition_log_odds".to_string(), 1.5);
        map.insert("europe_staphylococcus_aureus_acquisition_log_odds".to_string(), -0.2);
        map.insert("asia_staphylococcus_aureus_acquisition_log_odds".to_string(), 1.2);
        map.insert("south_america_staphylococcus_aureus_acquisition_log_odds".to_string(), 0.8);
        map.insert("oceania_staphylococcus_aureus_acquisition_log_odds".to_string(), 0.0);

        // streptococcus_pneumoniae - Predominantly community-acquired, high regional variation
        map.insert("africa_streptococcus_pneumoniae_acquisition_log_odds".to_string(), 2.2);
        map.insert("europe_streptococcus_pneumoniae_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_streptococcus_pneumoniae_acquisition_log_odds".to_string(), 1.8);
        map.insert("south_america_streptococcus_pneumoniae_acquisition_log_odds".to_string(), 1.2);
        map.insert("oceania_streptococcus_pneumoniae_acquisition_log_odds".to_string(), 0.3);

        // salmonella_enterica_serovar_typhi - Typhoid fever, highly endemic in certain regions
        map.insert("africa_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), 3.3);
        map.insert("europe_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), -1.8);
        map.insert("asia_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), 3.0);
        map.insert("south_america_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), 1.6);
        map.insert("oceania_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), -1.1);

        // salmonella_enterica_serovar_paratyphi_a - Paratyphoid fever, similar but less common
        map.insert("africa_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), 3.2);
        map.insert("europe_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), -2.1);
        map.insert("asia_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), 2.9);
        map.insert("south_america_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), 1.5);
        map.insert("oceania_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), -1.3);

        // invasive_non-typhoidal_salmonella_spp. - Bloodstream infections, especially in immunocompromised
        map.insert("africa_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), 4.5); // ~3 million cases/year in Africa; HIV/malnutrition
        map.insert("europe_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), -0.8);
        map.insert("asia_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), 2.3);
        map.insert("south_america_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), 1.7);
        map.insert("oceania_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), -0.3);

        // neisseria_gonorrhoeae - Sexually transmitted infection, moderate regional variation
        map.insert("africa_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 2.1);
        map.insert("europe_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 0.4);
        map.insert("asia_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 1.6);
        map.insert("south_america_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 1.2);
        map.insert("oceania_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 0.5);

        // streptococcus_pyogenes - Group A Strep, community-acquired, moderate regional variation
        map.insert("africa_streptococcus_pyogenes_acquisition_log_odds".to_string(), 1.8);
        map.insert("europe_streptococcus_pyogenes_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_streptococcus_pyogenes_acquisition_log_odds".to_string(), 1.4);
        map.insert("south_america_streptococcus_pyogenes_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_streptococcus_pyogenes_acquisition_log_odds".to_string(), 0.2);

        // streptococcus_agalactiae - Group B Strep, neonatal/maternal infections, moderate variation
        map.insert("africa_streptococcus_agalactiae_acquisition_log_odds".to_string(), 1.4);
        map.insert("europe_streptococcus_agalactiae_acquisition_log_odds".to_string(), 0.0);
        map.insert("asia_streptococcus_agalactiae_acquisition_log_odds".to_string(), 1.1);
        map.insert("south_america_streptococcus_agalactiae_acquisition_log_odds".to_string(), 0.7);
        map.insert("oceania_streptococcus_agalactiae_acquisition_log_odds".to_string(), 0.1);

        // haemophilus_influenzae - Respiratory pathogen, dramatically reduced by Hib vaccine
        map.insert("africa_haemophilus_influenzae_acquisition_log_odds".to_string(), 2.0);
        map.insert("europe_haemophilus_influenzae_acquisition_log_odds".to_string(), -0.3);
        map.insert("asia_haemophilus_influenzae_acquisition_log_odds".to_string(), 1.6);
        map.insert("south_america_haemophilus_influenzae_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_haemophilus_influenzae_acquisition_log_odds".to_string(), -0.2);

        // chlamydia_trachomatis - STI and trachoma, moderate regional variation
        map.insert("africa_chlamydia_trachomatis_acquisition_log_odds".to_string(), 2.1);
        map.insert("europe_chlamydia_trachomatis_acquisition_log_odds".to_string(), 0.2);
        map.insert("asia_chlamydia_trachomatis_acquisition_log_odds".to_string(), 1.7);
        map.insert("south_america_chlamydia_trachomatis_acquisition_log_odds".to_string(), 1.1);
        map.insert("oceania_chlamydia_trachomatis_acquisition_log_odds".to_string(), 0.6);

        // helicobacter_pylori - Gastric colonization with strong regional gradients
        map.insert("africa_helicobacter_pylori_acquisition_log_odds".to_string(), 3.2);
        map.insert("europe_helicobacter_pylori_acquisition_log_odds".to_string(), -0.2);
        map.insert("asia_helicobacter_pylori_acquisition_log_odds".to_string(), 2.8);
        map.insert("south_america_helicobacter_pylori_acquisition_log_odds".to_string(), 1.9);
        map.insert("oceania_helicobacter_pylori_acquisition_log_odds".to_string(), 0.4);

        // vibrio_cholerae - Waterborne disease, extreme regional variation
        map.insert("africa_vibrio_cholerae_acquisition_log_odds".to_string(), 3.5);
        map.insert("europe_vibrio_cholerae_acquisition_log_odds".to_string(), -3.0);
        map.insert("asia_vibrio_cholerae_acquisition_log_odds".to_string(), 2.8);
        map.insert("south_america_vibrio_cholerae_acquisition_log_odds".to_string(), 1.5);
        map.insert("oceania_vibrio_cholerae_acquisition_log_odds".to_string(), -2.2);

        // neisseria_meningitidis - Meningococcal disease, moderate regional variation with vaccine impact
        map.insert("africa_neisseria_meningitidis_acquisition_log_odds".to_string(), 1.5); // Calibrated down to avoid epidemic-scale baseline prevalence
        map.insert("europe_neisseria_meningitidis_acquisition_log_odds".to_string(), -0.8);
        map.insert("asia_neisseria_meningitidis_acquisition_log_odds".to_string(), -0.2);
        map.insert("south_america_neisseria_meningitidis_acquisition_log_odds".to_string(), -0.2);
        map.insert("oceania_neisseria_meningitidis_acquisition_log_odds".to_string(), -0.6);

        // listeria_monocytogenes - Foodborne pathogen, moderate regional variation
        map.insert("africa_listeria_monocytogenes_acquisition_log_odds".to_string(), 0.4);
        map.insert("europe_listeria_monocytogenes_acquisition_log_odds".to_string(), -0.2);
        map.insert("asia_listeria_monocytogenes_acquisition_log_odds".to_string(), 0.3);
        map.insert("south_america_listeria_monocytogenes_acquisition_log_odds".to_string(), 0.1);
        map.insert("oceania_listeria_monocytogenes_acquisition_log_odds".to_string(), -0.3);

        // clostridioides_difficile - Healthcare-associated, antibiotic-driven, modest regional variation
        map.insert("africa_clostridioides_difficile_acquisition_log_odds".to_string(), 0.2);
        map.insert("europe_clostridioides_difficile_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_clostridioides_difficile_acquisition_log_odds".to_string(), 0.1);
        map.insert("south_america_clostridioides_difficile_acquisition_log_odds".to_string(), 0.0);
        map.insert("oceania_clostridioides_difficile_acquisition_log_odds".to_string(), -0.2);

        // campylobacter_jejuni - Foodborne pathogen, moderate regional variation
        map.insert("africa_campylobacter_jejuni_acquisition_log_odds".to_string(), 2.5);
        map.insert("europe_campylobacter_jejuni_acquisition_log_odds".to_string(), 1.0);
        map.insert("asia_campylobacter_jejuni_acquisition_log_odds".to_string(), 1.9);
        map.insert("south_america_campylobacter_jejuni_acquisition_log_odds".to_string(), 1.5);
        map.insert("oceania_campylobacter_jejuni_acquisition_log_odds".to_string(), 0.6);

        // enterobacter_cloacae - Healthcare-associated Enterobacteriaceae, modest regional variation
        map.insert("africa_enterobacter_cloacae_acquisition_log_odds".to_string(), 0.2);
        map.insert("europe_enterobacter_cloacae_acquisition_log_odds".to_string(), -0.5);
        map.insert("asia_enterobacter_cloacae_acquisition_log_odds".to_string(), 0.0);
        map.insert("south_america_enterobacter_cloacae_acquisition_log_odds".to_string(), -0.2);
        map.insert("oceania_enterobacter_cloacae_acquisition_log_odds".to_string(), -0.3);

        // yersinia_enterocolitica - Foodborne/zoonotic pathogen, moderate regional variation
        map.insert("africa_yersinia_enterocolitica_acquisition_log_odds".to_string(), 0.2);
        map.insert("europe_yersinia_enterocolitica_acquisition_log_odds".to_string(), 0.0);
        map.insert("asia_yersinia_enterocolitica_acquisition_log_odds".to_string(), 0.0);
        map.insert("south_america_yersinia_enterocolitica_acquisition_log_odds".to_string(), 0.0);
        map.insert("oceania_yersinia_enterocolitica_acquisition_log_odds".to_string(), -0.1);

        // moraxella_catarrhalis - Respiratory pathogen, moderate regional variation
        map.insert("africa_moraxella_catarrhalis_acquisition_log_odds".to_string(), 0.9);
        map.insert("europe_moraxella_catarrhalis_acquisition_log_odds".to_string(), -0.3);
        map.insert("asia_moraxella_catarrhalis_acquisition_log_odds".to_string(), 0.5);
        map.insert("south_america_moraxella_catarrhalis_acquisition_log_odds".to_string(), 0.2);
        map.insert("oceania_moraxella_catarrhalis_acquisition_log_odds".to_string(), -0.2);

        // treponema_pallidum - Syphilis, moderate-high regional variation
        map.insert("treponema_pallidum_base_bacteria_level_change".to_string(), 0.18); // Very slow spirochete replication (30+ hour doubling)
        map.insert("africa_treponema_pallidum_acquisition_log_odds".to_string(), 0.3);
        map.insert("europe_treponema_pallidum_acquisition_log_odds".to_string(), -0.2);
        map.insert("asia_treponema_pallidum_acquisition_log_odds".to_string(), 0.0);
        map.insert("south_america_treponema_pallidum_acquisition_log_odds".to_string(), 0.1);
        map.insert("oceania_treponema_pallidum_acquisition_log_odds".to_string(), -0.3);

        // === [G] Microbiome carriage & clearance priors ===
        // Links microbiome carriage to infection rates and establishes decay rates so default
        // carriage prevalence matches surveillance. Individual bacteria can override any of these.
    // Bacteria-specific microbiome vs infection acquisition log odds
    // Values chosen so average carriage prevalence aligns with clinical carriage estimates
    // NOTE: Very high values (>7) send almost ALL acquisitions to carriage, blocking infections!
    map.insert("escherichia_coli_log_odds_microbiome_vs_infection".to_string(), 5.5); 
    map.insert("enterococcus_faecalis_log_odds_microbiome_vs_infection".to_string(), 10.0); 
    map.insert("enterococcus_faecium_log_odds_microbiome_vs_infection".to_string(), 11.0); 
    map.insert("klebsiella_pneumoniae_log_odds_microbiome_vs_infection".to_string(), 11.5); 
    map.insert("staphylococcus_aureus_log_odds_microbiome_vs_infection".to_string(), 11.0); 
    map.insert("staphylococcus_epidermidis_log_odds_microbiome_vs_infection".to_string(), 11.0); 
    map.insert("enterobacter_spp._log_odds_microbiome_vs_infection".to_string(), 6.0); 
    map.insert("enterobacter_cloacae_log_odds_microbiome_vs_infection".to_string(), 6.0); 
    map.insert("citrobacter_spp._log_odds_microbiome_vs_infection".to_string(), 5.0); 
    map.insert("proteus_spp._log_odds_microbiome_vs_infection".to_string(), 4.0); 
    map.insert("serratia_spp._log_odds_microbiome_vs_infection".to_string(), 5.0); 
    map.insert("morganella_spp._log_odds_microbiome_vs_infection".to_string(), 7.0); 
    map.insert("streptococcus_pneumoniae_log_odds_microbiome_vs_infection".to_string(), 7.0); 
    map.insert("haemophilus_influenzae_log_odds_microbiome_vs_infection".to_string(), 8.5); 
    map.insert("moraxella_catarrhalis_log_odds_microbiome_vs_infection".to_string(), 14.0); 
    map.insert("streptococcus_pyogenes_log_odds_microbiome_vs_infection".to_string(), 10.0); 
    map.insert("streptococcus_agalactiae_log_odds_microbiome_vs_infection".to_string(), 12.0); 
    map.insert("acinetobacter_baumannii_log_odds_microbiome_vs_infection".to_string(), 3.0); 
    map.insert("pseudomonas_aeruginosa_log_odds_microbiome_vs_infection".to_string(), 3.0); 
    map.insert("clostridioides_difficile_log_odds_microbiome_vs_infection".to_string(), 6.0); 
    map.insert("salmonella_enterica_serovar_typhi_log_odds_microbiome_vs_infection".to_string(), -1.0); 
    map.insert("salmonella_enterica_serovar_paratyphi_a_log_odds_microbiome_vs_infection".to_string(), 1.05); 
    map.insert("invasive_non-typhoidal_salmonella_spp._log_odds_microbiome_vs_infection".to_string(), 2.96); 
    map.insert("shigella_spp._log_odds_microbiome_vs_infection".to_string(), -1.5); 
    map.insert("vibrio_cholerae_log_odds_microbiome_vs_infection".to_string(), 2.96); 
    map.insert("campylobacter_jejuni_log_odds_microbiome_vs_infection".to_string(), 0.5); 
    map.insert("yersinia_enterocolitica_log_odds_microbiome_vs_infection".to_string(), 1.2); 
    map.insert("listeria_monocytogenes_log_odds_microbiome_vs_infection".to_string(), 5.0); 
    map.insert("neisseria_gonorrhoeae_log_odds_microbiome_vs_infection".to_string(), -1.0); 
    map.insert("chlamydia_trachomatis_log_odds_microbiome_vs_infection".to_string(), 2.0); 
    map.insert("treponema_pallidum_log_odds_microbiome_vs_infection".to_string(), 1.5); 
    map.insert("neisseria_meningitidis_log_odds_microbiome_vs_infection".to_string(), 11.5); 
    map.insert("helicobacter_pylori_log_odds_microbiome_vs_infection".to_string(), 6.65); 
    map.insert("mdr_mycobacterium_tuberculosis_log_odds_microbiome_vs_infection".to_string(), 3.0); 

    // Bacteria-specific microbiome clearance probabilities (per day)
    map.insert("escherichia_coli_microbiome_clearance_probability_per_day".to_string(), 0.005); // Persistent gut commensal; years-long colonization
    map.insert("enterococcus_faecalis_microbiome_clearance_probability_per_day".to_string(), 0.008); // Persistent gut flora; rarely cleared
    map.insert("enterococcus_faecium_microbiome_clearance_probability_per_day".to_string(), 0.06);
    map.insert("klebsiella_pneumoniae_microbiome_clearance_probability_per_day".to_string(), 0.03);
    map.insert("staphylococcus_aureus_microbiome_clearance_probability_per_day".to_string(), 0.03); // Nasal carriage persists weeks-months
    map.insert("enterobacter_spp._microbiome_clearance_probability_per_day".to_string(), 0.07);
    map.insert("enterobacter_cloacae_microbiome_clearance_probability_per_day".to_string(), 0.08);
    map.insert("citrobacter_spp._microbiome_clearance_probability_per_day".to_string(), 0.08);
    map.insert("proteus_spp._microbiome_clearance_probability_per_day".to_string(), 0.08);
    map.insert("serratia_spp._microbiome_clearance_probability_per_day".to_string(), 0.1);
    map.insert("morganella_spp._microbiome_clearance_probability_per_day".to_string(), 0.1);
    map.insert("streptococcus_pneumoniae_microbiome_clearance_probability_per_day".to_string(), 0.05);
    map.insert("haemophilus_influenzae_microbiome_clearance_probability_per_day".to_string(), 0.06);
    map.insert("moraxella_catarrhalis_microbiome_clearance_probability_per_day".to_string(), 0.05);
    map.insert("streptococcus_pyogenes_microbiome_clearance_probability_per_day".to_string(), 0.08);
    map.insert("streptococcus_agalactiae_microbiome_clearance_probability_per_day".to_string(), 0.06);
    map.insert("acinetobacter_baumannii_microbiome_clearance_probability_per_day".to_string(), 0.1);
    map.insert("pseudomonas_aeruginosa_microbiome_clearance_probability_per_day".to_string(), 0.12);
    map.insert("clostridioides_difficile_microbiome_clearance_probability_per_day".to_string(), 0.02); // Colonization persists months, especially elderly
    map.insert("salmonella_enterica_serovar_typhi_microbiome_clearance_probability_per_day".to_string(), 0.003); // Chronic carriers persist for years
    map.insert("salmonella_enterica_serovar_paratyphi_a_microbiome_clearance_probability_per_day".to_string(), 0.15);
    map.insert("invasive_non-typhoidal_salmonella_spp._microbiome_clearance_probability_per_day".to_string(), 0.12);
    map.insert("shigella_spp._microbiome_clearance_probability_per_day".to_string(), 0.15);
    map.insert("vibrio_cholerae_microbiome_clearance_probability_per_day".to_string(), 0.15);
    map.insert("campylobacter_jejuni_microbiome_clearance_probability_per_day".to_string(), 0.12);
    map.insert("yersinia_enterocolitica_microbiome_clearance_probability_per_day".to_string(), 0.25);
    map.insert("listeria_monocytogenes_microbiome_clearance_probability_per_day".to_string(), 0.1);
    map.insert("neisseria_gonorrhoeae_microbiome_clearance_probability_per_day".to_string(), 0.2);
    map.insert("chlamydia_trachomatis_microbiome_clearance_probability_per_day".to_string(), 0.2);
    map.insert("treponema_pallidum_microbiome_clearance_probability_per_day".to_string(), 0.35);
    map.insert("neisseria_meningitidis_microbiome_clearance_probability_per_day".to_string(), 0.05);
    map.insert("helicobacter_pylori_microbiome_clearance_probability_per_day".to_string(), 0.001); // Extremely persistent; decades without treatment
    map.insert("mdr_mycobacterium_tuberculosis_microbiome_clearance_probability_per_day".to_string(), 0.0015); // Latent carriage now clears slowly over years
    map.insert("bordetella_pertussis_microbiome_clearance_probability_per_day".to_string(), 0.2);
    map.insert("staphylococcus_epidermidis_clearance_probability_per_day".to_string(), 0.2);


        // Microbiome acquisition now uses infection acquisition parameters plus bacteria-specific offset
        // Fallback parameter for backward compatibility (should rarely be used if bacteria-specific params set)
    map.insert("log_odds_microbiome_vs_infection".to_string(), 2.0); // Fallback: modest carriage boost if no bacteria-specific param

        // Environmental resistance level for new acquisitions


        map.insert("max_resistance_level".to_string(), 1.0);
        map.insert("majority_r_evolution_rate_per_day_when_drug_present".to_string(), 0.18); // faster majority_r emergence under sustained therapy

        // === [H] Resistance emergence, transfer, and mechanism weights ===
        // Tunes how quickly resistance signals appear, decay, and propagate across mechanisms.
        // Use these defaults for broad behaviour; override targeted keys for specific bacteria/drugs.
        // Resistance Emergence and Decay Parameters
        // Resistance reversion parameter: probability per day that resistance reverts to 0 if no drug present
        map.insert("resistance_reversion_rate_per_day".to_string(), 0.0005); // Default: modest, slightly faster decay to curb relapse
        // Microbiome emergence rate: lower than infection emergence because microbiome bacteria
        // experience less intense selection pressure.
        // CALIBRATION: 0.001 gave 26%, 0.003 gave 59% - reverting to 0.001 to target 15-30%
        map.insert("microbiome_resistance_emergence_rate_per_day_baseline".to_string(), 0.001); // Calibrated for microbiome resistance emergence
        map.insert("resistance_emergence_bacteria_level_multiplier".to_string(), 0.08); // Multiplier for bacteria level's effect on emergence

        map.insert("resistance_emergence_pop_size_multiplier".to_string(), 30.0); // Debug knob to keep prevalence steady when population size changes - 3/5 x 500_000 / pop size

        map.insert("any_r_increase_rate_per_day_when_drug_present".to_string(), 0.045); // Growth rate of resistance signal while therapy is active
        map.insert("any_r_emergence_level_on_first_emergence".to_string(), 0.5); // The resistance level 'any_r' starts at upon emergence


        //  Microbiome Resistance Transfer Parameter
    map.insert("microbiome_resistance_transfer_probability_per_day".to_string(), 0.0025); // Probability per day for resistance transfer between infection and microbiome

        // --- Multi-Drug Resistance Emergence Penalty Parameters ---
        // When multiple drugs are active, resistance emergence is reduced because mutations
        // must confer resistance to ALL active drugs to provide survival advantage
        map.insert("multi_drug_penalty_for_single_drug_resistance".to_string(), 0.05); // Penalty when resistance affects only 1 of multiple active drugs (5% survival)
        map.insert("multi_drug_penalty_for_partial_cross_resistance".to_string(), 0.3); // Penalty when resistance affects some but not all active drugs (30% survival)
        map.insert("multi_drug_penalty_threshold_num_drugs".to_string(), 2.0); // Minimum number of active drugs to trigger multi-drug penalty

        // --- Resistance Mechanisms Parameters ---
        // Baseline emergence rates for specific resistance mechanisms (per day when drug present)
        // Empirical basis: rates vary by mechanism complexity and genetic requirements
        // Common mechanisms (single mutations, regulatory changes): ~1e-6
        // Mobile genetic elements: ~1e-7 to 1e-8
        // Complex resistance clusters: ~1e-9
        map.insert("resistance_mechanism_target_site_mutation_emergence_rate".to_string(), 0.003); // Point mutations - most common
        map.insert("resistance_mechanism_efflux_overexpression_emergence_rate".to_string(), 0.003); // Regulatory mutations relatively common
        map.insert("resistance_mechanism_reduced_permeability_emergence_rate".to_string(), 0.003); // Porin loss is common
        map.insert("resistance_mechanism_qnr_emergence_rate".to_string(), 0.003); // Mobile genetic element acquisition
        map.insert("resistance_mechanism_erm_methylation_emergence_rate".to_string(), 0.003); // Common in gram-positives
        map.insert("resistance_mechanism_esbl_emergence_rate".to_string(), 0.001); // Requires specific gene mutations
        map.insert("resistance_mechanism_ampc_emergence_rate".to_string(), 0.003); // Chromosomal or plasmid-mediated
        map.insert("resistance_mechanism_meca_emergence_rate".to_string(), 0.001); // Requires SCCmec element acquisition
        map.insert("resistance_mechanism_carbapenemase_emergence_rate".to_string(), 0.001); // Rare, high-level resistance genes
        map.insert("resistance_mechanism_van_type_emergence_rate".to_string(), 0.001); // Complex vanA/vanB resistance cluster
        map.insert("resistance_mechanism_16s_methyltransferase_emergence_rate".to_string(), 0.001); // Rare, high-level aminoglycoside resistance

        // Resistance enhancement multipliers: how much each mechanism increases resistance level
        map.insert("resistance_mechanism_esbl_enhancement_multiplier".to_string(), 0.4); // Adds 40% resistance
        map.insert("resistance_mechanism_carbapenemase_enhancement_multiplier".to_string(), 0.6); // Adds 60% resistance
        map.insert("resistance_mechanism_ampc_enhancement_multiplier".to_string(), 0.3); // Adds 30% resistance
        map.insert("resistance_mechanism_16s_methyltransferase_enhancement_multiplier".to_string(), 0.5); // Adds 50% resistance
        map.insert("resistance_mechanism_qnr_enhancement_multiplier".to_string(), 0.2); // Adds 20% resistance (low-level)
        map.insert("resistance_mechanism_efflux_overexpression_enhancement_multiplier".to_string(), 0.2); // Adds 20% resistance
        map.insert("resistance_mechanism_erm_methylation_enhancement_multiplier".to_string(), 0.5); // Adds 50% resistance
        map.insert("resistance_mechanism_van_type_enhancement_multiplier".to_string(), 0.8); // Adds 80% resistance (high-level)
        map.insert("resistance_mechanism_meca_enhancement_multiplier".to_string(), 0.7); // Adds 70% resistance
        map.insert("resistance_mechanism_reduced_permeability_enhancement_multiplier".to_string(), 0.2); // Adds 20% resistance
        map.insert("resistance_mechanism_target_site_mutation_enhancement_multiplier".to_string(), 0.4); // Adds 40% resistance

        map.insert("mechanism_assignment_probability_on_any_r_gain".to_string(), 0.8); // Default 80%

        // Mechanism-specific fitness costs (reversion rates per day when drug absent)
        // High-cost mechanisms: Metabolically expensive enzymes
        map.insert("resistance_mechanism_carbapenemase_reversion_rate".to_string(), 0.001); // High cost - large, expensive enzymes
        map.insert("resistance_mechanism_van_type_reversion_rate".to_string(), 0.002); // High cost - complex resistance pathway

        // Medium-cost mechanisms: Moderate metabolic burden
        map.insert("resistance_mechanism_esbl_reversion_rate".to_string(), 0.0006); // Medium cost - beta-lactamase production
        map.insert("resistance_mechanism_meca_reversion_rate".to_string(), 0.0009 ); // Medium cost - altered PBP
        map.insert("resistance_mechanism_erm_methylation_reversion_rate".to_string(), 0.0006); // Medium cost - methyltransferase
        map.insert("resistance_mechanism_16s_methyltransferase_reversion_rate".to_string(), 0.0006); // Medium cost - rRNA modification

        // Low-cost mechanisms: Minimal fitness burden
        map.insert("resistance_mechanism_qnr_reversion_rate".to_string(), 0.0001); // Low cost - point mutation effect
        map.insert("resistance_mechanism_target_site_mutation_reversion_rate".to_string(), 0.0002); // Low cost - single nucleotide changes
        map.insert("resistance_mechanism_reduced_permeability_reversion_rate".to_string(), 0.0003 ); // Low cost - adaptive change
        map.insert("resistance_mechanism_efflux_overexpression_reversion_rate".to_string(), 0.0005); // Low-medium cost - energy for pumping
        map.insert("resistance_mechanism_ampc_reversion_rate".to_string(), 0.00015); // Low cost - chromosomal enzyme

        // Testing Parameters
        map.insert("bacterial_testing_available_from_day".to_string(), 5478.0); // 5478.0  1945 (15 years after 1930) - Bacterial culture/identification becomes available
        map.insert("resistance_testing_available_from_day".to_string(), 9131.0); // 9131.0  1955 (25 years after 1930) - Antibiotic susceptibility testing becomes available
        map.insert("test_delay_days".to_string(), 3.0);
        map.insert("test_rate_per_day".to_string(), 0.2);  // 0.15

        // --- Test result and test_r logic parameters ---
        map.insert("prob_test_r_done".to_string(), 0.95); // Probability test is actually done (per day eligible)
        map.insert("test_r_error_probability".to_string(), 0.02); // Probability of error in test result
        map.insert("test_r_error_value".to_string(), 0.25); // Value to use for error in test_r

        // Syndrome-specific initiation multipliers
        map.insert("syndrome_3_initiation_multiplier".to_string(), 10.0); // Respiratory syndrome
        map.insert("syndrome_7_initiation_multiplier".to_string(), 8.0);  // Gastrointestinal syndrome
        map.insert("syndrome_8_initiation_multiplier".to_string(), 12.0); // Genital syndrome (example ID)

        // Empiric drug scoring tables (clinician-facing heuristics per syndrome ID)
        // These preserve pre-refactor prescribing patterns when organism is unknown.
        let empiric_syndrome_templates: &[(usize, &[(&str, f64)])] = &[
            // 1 = UTI / Genitourinary
            (
                1,
                &[
                    ("nitrofurantoin", 18.0),
                    ("trim_sulf", 14.0),
                    ("ciprofloxacin", 12.0),
                    ("levofloxacin", 10.0),
                    ("amoxicillin_clavulanate", 9.0),
                    ("amoxicillin", 7.0),
                    ("ampicillin", 6.0),
                    ("ceftriaxone", 8.0),
                    ("cefuroxime", 7.0),
                    ("piperacillin_tazobactam", 5.0),
                    ("cefepime", 4.0),
                    ("ceftazidime", 4.0),
                    ("meropenem", 4.0),
                    ("imipenem_c", 4.0),
                    ("ertapenem", 4.0),
                    ("meropenem_vaborbactam", 3.0),
                    ("ceftazidime_avibactam", 3.0),
                    ("colistin", 0.2),
                    ("vancomycin", 0.1),
                    ("linezolid", 0.1),
                ],
            ),
            // 2 = Skin / soft tissue
            (
                2,
                &[
                    ("penicilling", 14.0),
                    ("ampicillin", 11.0),
                    ("amoxicillin", 12.0),
                    ("amoxicillin_clavulanate", 12.0),
                    ("cephalexin", 13.0),
                    ("cefazolin", 12.0),
                    ("clindamycin", 12.0),
                    ("trim_sulf", 9.0),
                    ("doxycycline", 9.0),
                    ("minocycline", 9.0),
                    ("linezolid", 10.0),
                    ("tedizolid", 9.0),
                    ("dalbavancin", 9.0),
                    ("vancomycin", 11.0),
                    ("quinu_dalfo", 8.0),
                    ("rifampicin", 6.0),
                    ("ciprofloxacin", 4.0),
                    ("piperacillin_tazobactam", 3.0),
                ],
            ),
            // 3 = Respiratory
            (
                3,
                &[
                    ("penicilling", 9.0),
                    ("ampicillin", 9.5),
                    ("amoxicillin", 11.0),
                    ("amoxicillin_clavulanate", 12.0),
                    ("cefuroxime", 8.5),
                    ("ceftriaxone", 9.5),
                    ("cefepime", 7.5),
                    ("piperacillin_tazobactam", 7.0),
                    ("meropenem", 6.0),
                    ("imipenem_c", 6.0),
                    ("azithromycin", 11.5),
                    ("clarithromycin", 10.5),
                    ("erythromycin", 8.5),
                    ("doxycycline", 8.0),
                    ("minocycline", 7.0),
                    ("levofloxacin", 11.0),
                    ("moxifloxacin", 11.0),
                    ("ofloxacin", 8.0),
                    ("linezolid", 7.0),
                    ("vancomycin", 6.5),
                ],
            ),
            // 4 = Bloodstream / bacteremia
            (
                4,
                &[
                    ("piperacillin_tazobactam", 14.0),
                    ("meropenem", 13.0),
                    ("imipenem_c", 13.0),
                    ("meropenem_vaborbactam", 13.0),
                    ("ceftazidime_avibactam", 12.5),
                    ("cefepime", 12.0),
                    ("ceftazidime", 11.0),
                    ("ceftriaxone", 10.0),
                    ("vancomycin", 11.0),
                    ("linezolid", 10.0),
                    ("tedizolid", 9.0),
                    ("dalbavancin", 8.0),
                    ("quinu_dalfo", 8.5),
                    ("gentamicin", 7.0),
                    ("tobramycin", 6.5),
                    ("amikacin", 7.0),
                    ("colistin", 6.0),
                    ("ciprofloxacin", 6.0),
                    ("levofloxacin", 5.5),
                    ("rifampicin", 4.0),
                ],
            ),
            // 5 = Intra-abdominal
            (
                5,
                &[
                    ("metronidazole", 15.0),
                    ("piperacillin_tazobactam", 13.0),
                    ("ampicillin_sulbactam", 11.0),
                    ("amoxicillin_clavulanate", 10.0),
                    ("meropenem", 13.0),
                    ("imipenem_c", 12.5),
                    ("ertapenem", 11.0),
                    ("ceftazidime", 9.0),
                    ("cefepime", 9.0),
                    ("ceftriaxone", 9.0),
                    ("ceftazidime_avibactam", 10.0),
                    ("meropenem_vaborbactam", 10.0),
                    ("ciprofloxacin", 7.0),
                    ("levofloxacin", 6.5),
                    ("trim_sulf", 4.0),
                    ("colistin", 3.5),
                ],
            ),
            // 6 = Central nervous system
            (
                6,
                &[
                    ("ceftriaxone", 15.0),
                    ("ceftazidime", 12.0),
                    ("cefepime", 12.0),
                    ("penicilling", 9.0),
                    ("ampicillin", 11.0),
                    ("vancomycin", 13.0),
                    ("linezolid", 10.0),
                    ("meropenem", 11.0),
                    ("imipenem_c", 10.0),
                    ("chlorampheni", 9.0),
                    ("rifampicin", 7.0),
                    ("piperacillin_tazobactam", 6.0),
                ],
            ),
            // 7 = Gastrointestinal (non-invasive)
            (
                7,
                &[
                    ("ciprofloxacin", 12.0),
                    ("levofloxacin", 10.0),
                    ("azithromycin", 10.0),
                    ("doxycycline", 8.5),
                    ("minocycline", 6.5),
                    ("trim_sulf", 8.5),
                    ("furazolidone", 11.0),
                    ("metronidazole", 12.0),
                    ("rifampicin", 5.0),
                    ("ampicillin", 4.0),
                    ("amoxicillin", 4.5),
                    ("amoxicillin_clavulanate", 5.0),
                ],
            ),
            // 8 = Genital / pelvic
            (
                8,
                &[
                    ("azithromycin", 13.0),
                    ("doxycycline", 12.0),
                    ("penicilling", 11.0),
                    ("ceftriaxone", 13.0),
                    ("cefuroxime", 9.0),
                    ("amoxicillin", 8.0),
                    ("amoxicillin_clavulanate", 8.0),
                    ("metronidazole", 12.0),
                    ("clindamycin", 9.0),
                    ("ciprofloxacin", 7.0),
                    ("levofloxacin", 6.5),
                    ("trim_sulf", 5.0),
                    ("rifampicin", 4.0),
                ],
            ),
            // 9 = Bone / joint / hardware-associated
            (
                9,
                &[
                    ("cefazolin", 13.0),
                    ("cephalexin", 11.0),
                    ("ceftriaxone", 11.0),
                    ("vancomycin", 12.0),
                    ("linezolid", 11.0),
                    ("tedizolid", 10.0),
                    ("dalbavancin", 10.0),
                    ("clindamycin", 10.0),
                    ("ciprofloxacin", 9.0),
                    ("levofloxacin", 9.0),
                    ("rifampicin", 9.0),
                    ("trim_sulf", 8.0),
                    ("meropenem", 7.0),
                    ("piperacillin_tazobactam", 6.5),
                ],
            ),
            // 10 = Other severe / device-related catch-all
            (
                10,
                &[
                    ("piperacillin_tazobactam", 8.0),
                    ("cefepime", 8.0),
                    ("ceftriaxone", 8.0),
                    ("meropenem", 8.0),
                    ("imipenem_c", 8.0),
                    ("vancomycin", 8.0),
                    ("linezolid", 7.0),
                    ("ciprofloxacin", 7.0),
                    ("azithromycin", 6.0),
                ],
            ),
        ];

        for (syndrome_id, entries) in empiric_syndrome_templates {
            for &(drug, score) in *entries {
                if DRUG_SHORT_NAMES.contains(&drug) {
                    map.insert(
                        format!("syndrome_{}_empiric_drug_{}_score", syndrome_id, drug),
                        score,
                    );
                }
            }
        }

    // Hospitalization Parameters
    map.insert("hospitalization_baseline_rate_per_day".to_string(), 0.00003); // Baseline daily probability tuned for ~0.3%-0.5% prevalence
    map.insert("hospitalization_age_multiplier_per_day".to_string(), 0.00000005); // Incremental daily hospitalization probability per day of age (~0.07%/day at age 40)
    map.insert("hospitalization_recovery_rate_per_day".to_string(), 0.28); // Slightly shorter stays (~3.6 day avg) to reinforce target occupancy
    map.insert("hospitalization_max_days".to_string(), 30.0); // Max days in hospital before forced discharge (as fallback)
    map.insert("hospitalization_sepsis_admission_multiplier".to_string(), 80.0); // Sepsis substantially increases admission odds
    map.insert("hospitalization_prevent_discharge_with_sepsis".to_string(), 1.0); // 1.0 = block discharge with sepsis, 0.0 = allow discharge

        // Testing Framework Parameters
        // Base testing rates (modern era baseline)
        map.insert("bacterial_testing_base_rate_per_day".to_string(), 0.15); // Modern baseline rate for bacterial identification
        map.insert("resistance_testing_base_rate_per_day".to_string(), 0.95); // Probability of resistance testing given bacterial identification

        // Hospital status multipliers
        map.insert("bacterial_testing_hospital_multiplier".to_string(), 8.0); // Hospitalized patients 8x more likely to get bacterial testing
        map.insert("resistance_testing_hospital_multiplier".to_string(), 5.0); // Hospitalized patients 5x more likely to get resistance testing

        // Regional resource multipliers for testing access
        map.insert("north_america_testing_multiplier".to_string(), 1.1);
        map.insert("europe_testing_multiplier".to_string(), 1.2);
        map.insert("asia_testing_multiplier".to_string(), 0.7);
        map.insert("south_america_testing_multiplier".to_string(), 0.6);
        map.insert("oceania_testing_multiplier".to_string(), 0.8);
        map.insert("africa_testing_multiplier".to_string(), 0.3); // Limited lab infrastructure

        // Clinical status multipliers
        map.insert("testing_immunosuppressed_multiplier".to_string(), 2.5); // Immunosuppressed patients get more testing
        map.insert("testing_sepsis_multiplier".to_string(), 4.0); // Sepsis patients get urgent testing

        // Temporal adoption parameters (testing evolution over time)
        // Using S-curve (sigmoid) model for realistic technology adoption
        // Bacterial testing temporal evolution
        map.insert("bacterial_testing_initial_adoption_rate".to_string(), 0.1); // 1945: 10% of modern rates
        map.insert("bacterial_testing_adoption_rate_per_year".to_string(), 0.025); // DEPRECATED: kept for backward compatibility
        map.insert("bacterial_testing_max_temporal_multiplier".to_string(), 1.0); // Cap at modern rates (100%)

        // Resistance testing temporal evolution (slower adoption)
        map.insert("resistance_testing_initial_adoption_rate".to_string(), 0.05); // 1955: 5% of modern rates
        map.insert("resistance_testing_adoption_rate_per_year".to_string(), 0.015); // DEPRECATED: kept for backward compatibility
        map.insert("resistance_testing_max_temporal_multiplier".to_string(), 1.0); // Cap at modern rates (100%)

        // initiate travel
        map.insert("travel_probability_per_day".to_string(), 0.00005);

        // Region-specific travel multipliers based on income/development level
        // Higher income regions have higher outbound travel rates
        map.insert("north_america_travel_multiplier".to_string(), 3.0);  // High income, high travel
        map.insert("europe_travel_multiplier".to_string(), 3.5);         // High income, highest travel rates
        map.insert("oceania_travel_multiplier".to_string(), 2.5);        // High income, high travel
        map.insert("asia_travel_multiplier".to_string(), 1.5);          // Mixed income levels, moderate travel
        map.insert("south_america_travel_multiplier".to_string(), 0.8);  // Middle income, lower travel
        map.insert("africa_travel_multiplier".to_string(), 0.3);        // Lower income, lowest travel rates



        // Default Initial Drug Levels and Double Dose Multipliers for ALL Drugs
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_initial_level", drug), 10.0); // Default initial level for each drug
            map.insert(format!("drug_{}_double_dose_multiplier", drug), 2.0); // Default double dose multiplier
            map.insert(format!("drug_{}_spectrum_breadth", drug), 3.0); // Default spectrum: 1.0=narrow, 5.0=very broad
        }

        // Bacterial Identification Effect Parameters
    map.insert("empiric_therapy_broad_spectrum_bonus".to_string(), 0.85); // Make broad-spectrum empiric choices noticeably less favored than narrow options
        map.insert("empiric_therapy_ineffective_drug_penalty".to_string(), 0.001); // STRENGTHENED: Heavy penalty for drugs ineffective against actual pathogens (empirical)
        map.insert("targeted_therapy_narrow_spectrum_bonus".to_string(), 5.0); // Further reward narrow agents once pathogen identified
        map.insert("targeted_therapy_broad_spectrum_penalty".to_string(), 0.1); // Stronger penalty for broad-spectrum drugs when bacteria identified
        map.insert("targeted_therapy_ineffective_drug_penalty".to_string(), 0.001); // STRENGTHENED: Strong penalty for drugs ineffective against identified bacteria

        // Regional Resistance Surveillance Parameters for Drug Choice
        // Penalties applied during empirical therapy based on local resistance rates
        map.insert("regional_resistance_penalty_very_high".to_string(), 0.2); // Penalty when >50% regional resistance
        map.insert("regional_resistance_penalty_high".to_string(), 0.4); // Penalty when 30-50% regional resistance
        map.insert("regional_resistance_penalty_moderate".to_string(), 0.7); // Penalty when 10-30% regional resistance
        map.insert("regional_resistance_threshold_very_high".to_string(), 0.6); // Threshold for very high resistance (60%)
        map.insert("regional_resistance_threshold_high".to_string(), 0.45); // Threshold for high resistance (45%)
        map.insert("regional_resistance_threshold_moderate".to_string(), 0.1); // Threshold for moderate resistance (10%)

        // Drug Spectrum Classifications (1.0=narrow, 5.0=very broad)
        map.insert("drug_colistin_spectrum_breadth".to_string(), 4.0); // Broad spectrum (mainly Gram-negative)
        map.insert("drug_penicilling_spectrum_breadth".to_string(), 2.0); // Narrow spectrum
        map.insert("drug_amoxicillin_spectrum_breadth".to_string(), 3.0); // Medium spectrum
        map.insert("drug_azithromycin_spectrum_breadth".to_string(), 4.0); // Broad spectrum
        map.insert("drug_ciprofloxacin_spectrum_breadth".to_string(), 4.5); // Very broad spectrum
        map.insert("drug_trim_sulf_spectrum_breadth".to_string(), 3.5); // Medium-broad spectrum
        map.insert("drug_meropenem_spectrum_breadth".to_string(), 5.0); // Very broad spectrum (carbapenem)
        map.insert("drug_cefepime_spectrum_breadth".to_string(), 4.0); // Broad spectrum (4th gen cephalosporin)
        map.insert("drug_vancomycin_spectrum_breadth".to_string(), 2.5); // Narrow-medium spectrum (gram-positive only)
        map.insert("drug_linezolid_spectrum_breadth".to_string(), 2.0); // Narrow spectrum (gram-positive only)
        map.insert("drug_ceftriaxone_spectrum_breadth".to_string(), 4.0); // Broad spectrum (3rd gen cephalosporin)


        // NEW: Logistic Sepsis Risk Parameters (replacing old linear model)
        map.insert("sepsis_baseline_log_odds".to_string(), -14.0); // Fallback baseline for organisms without explicit intercept


        // sepsis rates     

        // Bacteria-specific sepsis baseline log-odds (best-guess placeholders calibrated by clinical severity)
        let bacteria_sepsis_baseline_overrides: &[(&str, f64)] = &[
            ("acinetobacter_baumannii", -6.0),
            ("citrobacter_spp.", -6.0),
            ("enterobacter_spp.", -6.0),
            ("enterococcus_faecalis", -6.0),
            ("enterococcus_faecium", -6.0),
            ("escherichia_coli", -15.5),
            ("klebsiella_pneumoniae", -6.0),
            ("morganella_spp.", -6.0),
            ("proteus_spp.", -6.0),
            ("serratia_spp.", -6.0),
            ("pseudomonas_aeruginosa", -6.0),
            ("stenotrophomonas_maltophilia", -6.0),
            ("staphylococcus_aureus", -5.0),
            ("staphylococcus_epidermidis", -6.0),
            ("streptococcus_pneumoniae", -14.0),
            ("salmonella_enterica_serovar_typhi", -9.1),
            ("salmonella_enterica_serovar_paratyphi_a", -6.0),
            ("invasive_non-typhoidal_salmonella_spp.", -7.0),
            ("shigella_spp.", -9.0),
            ("neisseria_gonorrhoeae", -18.0),
            ("streptococcus_pyogenes", -6.0),
            ("streptococcus_agalactiae", -6.0),
            ("haemophilus_influenzae", -15.0),
            ("chlamydia_trachomatis", -7.0),
            ("vibrio_cholerae", -6.0),
            ("neisseria_meningitidis", -13.2),
            ("listeria_monocytogenes", -6.0),
            ("clostridioides_difficile", -9.0),
            ("campylobacter_jejuni", -12.0),
            ("enterobacter_cloacae", -6.0),
            ("yersinia_enterocolitica", -6.0),
            ("moraxella_catarrhalis", -6.0),
            ("treponema_pallidum", -10.0),
            ("bordetella_pertussis", -6.0),
            ("helicobacter_pylori", -220.0),
            ("mdr_mycobacterium_tuberculosis", -38.0),
        ];

        for (bacteria, log_odds) in bacteria_sepsis_baseline_overrides {
            map.insert(format!("{}_sepsis_baseline_log_odds", bacteria), *log_odds);
        }
        map.insert("log_odds_sepsis_infection_level".to_string(), 2.0); // Log odds increase per unit bacterial level
        // === [I] Clinical outcome scalars (mortality, sepsis, toxicity) ===
        // Collects mortality/sepsis odds adjustments together so scenario designers can reason about
        // outcome severity in one place. These parameters shape the probability of severe outcomes
        // once infection is established.
        map.insert("log_odds_sepsis_infection_duration".to_string(), 0.001); // Log odds increase per day of infection duration

        // --- AGE-DEPENDENT SEPSIS LOG-ODDS (global baseline + age deltas + bacteria-age deltas) ---
        map.insert("sepsis_age_log_odds_baseline".to_string(), 0.0); // Reference intercept for age adjustments
        map.insert("sepsis_age_log_odds_neonatal".to_string(), 1.10); // ln(3.0): neonates ~3x odds vs reference
        map.insert("sepsis_age_log_odds_pediatric".to_string(), 0.18); // ln(1.2): pediatrics modestly higher odds
        map.insert("sepsis_age_log_odds_young_adult".to_string(), 0.0); // ln(1.0): young adults at reference odds
        map.insert("sepsis_age_log_odds_elderly".to_string(), 0.69); // ln(2.0): elderly ~2x odds vs reference

        // Additional per-bacteria age interactions (values are additive deltas beyond the general age effect)
        // NEONATAL (0-28 days)
        map.insert("streptococcus_agalactiae_neonatal_sepsis_log_odds".to_string(), 0.981); // ln(8) - ln(3)
        map.insert("escherichia_coli_neonatal_sepsis_log_odds".to_string(), 0.511); // ln(5) - ln(3)
        map.insert("listeria_monocytogenes_neonatal_sepsis_log_odds".to_string(), 0.693); // ln(6) - ln(3)
        map.insert("enterococcus_faecalis_neonatal_sepsis_log_odds".to_string(), 0.0); // ln(3) - ln(3)
        map.insert("staphylococcus_aureus_neonatal_sepsis_log_odds".to_string(), 0.288); // ln(4) - ln(3)

        // PEDIATRIC (1 month - 18 years)
        map.insert("streptococcus_pneumoniae_pediatric_sepsis_log_odds".to_string(), 0.916); // ln(3) - ln(1.2)
        map.insert("haemophilus_influenzae_pediatric_sepsis_log_odds".to_string(), 0.734); // ln(2.5) - ln(1.2)
        map.insert("neisseria_meningitidis_pediatric_sepsis_log_odds".to_string(), 0.511); // ln(2) - ln(1.2)
        map.insert("staphylococcus_aureus_pediatric_sepsis_log_odds".to_string(), 0.511); // ln(2) - ln(1.2)

        // ELDERLY (65+ years)
        map.insert("streptococcus_pneumoniae_elderly_sepsis_log_odds".to_string(), 0.693); // ln(4) - ln(2)
        map.insert("escherichia_coli_elderly_sepsis_log_odds".to_string(), 0.223); // ln(2.5) - ln(2)
        map.insert("klebsiella_pneumoniae_elderly_sepsis_log_odds".to_string(), 0.405); // ln(3) - ln(2)
        map.insert("pseudomonas_aeruginosa_elderly_sepsis_log_odds".to_string(), 0.560); // ln(3.5) - ln(2)
        map.insert("acinetobacter_baumannii_elderly_sepsis_log_odds".to_string(), 0.405); // ln(3) - ln(2)
        map.insert("enterococcus_faecium_elderly_sepsis_log_odds".to_string(), 0.336); // ln(2.8) - ln(2)
        map.insert("staphylococcus_aureus_elderly_sepsis_log_odds".to_string(), 0.470); // ln(3.2) - ln(2)

        // YOUNG ADULT (18-65 years)
        map.insert("neisseria_meningitidis_young_adult_sepsis_log_odds".to_string(), 0.336); // ln(1.4) - ln(1.0)
        map.insert("staphylococcus_aureus_young_adult_sepsis_log_odds".to_string(), 0.588); // ln(1.8) - ln(1.0)

        // --- REGIONAL SEPSIS RISK FACTORS ---
        // Account for healthcare infrastructure, population density, socioeconomic factors

        // HEALTHCARE ACCESS AND QUALITY MODIFIERS
        map.insert("log_odds_sepsis_region_a".to_string(), -0.3); // Higher resource region - better sepsis recognition/treatment
        map.insert("log_odds_sepsis_region_b".to_string(), 0.2);  // Lower resource region - delayed recognition/limited resources


        // Syndrome-specific sepsis risk parameters (infectious site effects)
        // CALIBRATED: Reduced sepsis risk for high-volume low-severity syndromes to achieve
        // ~14 million infection deaths/year (vs 18.57M before adjustment)
        // Key changes: UTI -2.0→-3.5, Skin -1.0→-2.5, Respiratory 0.0→-0.5, GI -0.5→-1.5, Genital -1.5→-2.5
        map.insert("log_odds_syndrome_1_sepsis".to_string(), -2.0); // UTI/Genitourinary: Very low sepsis risk (urosepsis is rare)
        map.insert("log_odds_syndrome_2_sepsis".to_string(), -1.0); // Skin/Soft tissue: Low sepsis risk (cellulitis rarely → sepsis)
        map.insert("log_odds_syndrome_3_sepsis".to_string(),  0.0); // Respiratory: Slightly below reference (mixes pneumonia with bronchitis)
        map.insert("log_odds_syndrome_4_sepsis".to_string(), 1.5);  // Bloodstream/Bacteremia: Much higher sepsis risk (UNCHANGED - appropriate)
        map.insert("log_odds_syndrome_5_sepsis".to_string(), 0.8);  // Intra-abdominal: Higher sepsis risk (UNCHANGED - peritonitis is serious)
        map.insert("log_odds_syndrome_6_sepsis".to_string(), 1.2);  // Central nervous system: High sepsis risk (UNCHANGED - meningitis is serious)
        map.insert("log_odds_syndrome_7_sepsis".to_string(), -0.5); // Gastrointestinal: Lower sepsis risk (most self-limiting)
        map.insert("log_odds_syndrome_8_sepsis".to_string(), -1.5); // Genital: Low sepsis risk (localized infections)
        map.insert("log_odds_syndrome_9_sepsis".to_string(), 0.5);  // Bone/Joint: Moderately higher sepsis risk (UNCHANGED)

        // // Background Mortality Parameters (Age, Region, and Sex dependent)

        // These parameters are on the log-odds scale.
        map.insert("background_mortality_baseline_log_odds".to_string(), -14.0); // -15.5
        map.insert("log_odds_mortality_per_year_of_age".to_string(), 0.04); // 0.04  Odds of dying increase by ~4% per year (exp(0.04) ≈ 1.04)
        map.insert("log_odds_mortality_per_year_of_age_squared".to_string(), 0.05); // 0.05  Additional non-linear effect for elderly

        // Time-varying background mortality (1930-2035): reflects dramatic mortality improvements over 105 years
        // Based on historical life expectancy trends showing ~3x mortality decline from 1930 to 2035
        map.insert("mortality_baseline_1930_multiplier".to_string(), 3.0);  // 3x higher mortality in 1930 (life expectancy ~35-45 years)
        map.insert("mortality_baseline_2035_multiplier".to_string(), 1.0);  // Modern baseline by 2035 (life expectancy ~70-80 years)
        map.insert("mortality_improvement_half_life_years".to_string(), 35.0); // Half-life of mortality improvement (exponential decay)

        // Region-specific log-odds adjustments. ln(1.0) = 0.
        map.insert("log_odds_mortality_region_north_america".to_string(), 0.0);      // Reference
        map.insert("log_odds_mortality_region_south_america".to_string(), 0.26);     // ln(1.3) - increased from 1.2
        map.insert("log_odds_mortality_region_africa".to_string(), 0.69);            // ln(2.0) - higher mortality burden
        map.insert("log_odds_mortality_region_asia".to_string(), 0.18);              // ln(1.2) - increased from 1.1
        map.insert("log_odds_mortality_region_europe".to_string(), -0.105);          // ln(0.9) - kept same
        map.insert("log_odds_mortality_region_oceania".to_string(), 0.0);           // Reference

        // Sex-specific log-odds adjustments.
        map.insert("log_odds_mortality_sex_male".to_string(), 0.095);    // ln(1.1)
        map.insert("log_odds_mortality_sex_female".to_string(), -0.105); // ln(0.9)

        // Additional background mortality risk factors (log-odds).
        map.insert("log_odds_mortality_immunosuppressed".to_string(), 0.916); // ln(2.5)
        map.insert("log_odds_mortality_hospitalized".to_string(), 0.262);     // ln(1.3)


        //  Immunosuppression Onset and Recovery Rates
        // Temporary immunodeficiency (e.g., chemotherapy, short-term steroids)
        map.insert("temporary_immunosuppression_onset_rate_per_day".to_string(), 0.00005);  // Calibrated onset for ~5% prevalence
        map.insert("temporary_immunosuppression_recovery_rate_per_day".to_string(), 0.01);   // Faster recovery (10x faster)

        // Chronic immunodeficiency (e.g., HIV, genetic disorders, organ transplant)
        map.insert("chronic_immunosuppression_onset_rate_per_day".to_string(), 0.00006);     // Lower onset rate for chronic
        map.insert("chronic_immunosuppression_recovery_rate_per_day".to_string(), 0.0012);   // Faster chronic recovery to offset prevalence

        // Age effect on immunodeficiency type assignment (probability of chronic vs temporary at onset)
        map.insert("chronic_immunodeficiency_probability_age_0_1".to_string(), 0.3);   // Infants: higher chance of genetic/congenital
        map.insert("chronic_immunodeficiency_probability_age_1_18".to_string(), 0.2);  // Children: moderate chance
        map.insert("chronic_immunodeficiency_probability_age_18_65".to_string(), 0.4); // Adults: higher chance (HIV, transplants)
        map.insert("chronic_immunodeficiency_probability_age_65_plus".to_string(), 0.6); // Elderly: highest chance (multiple conditions)

        // Prophylactic antibiotic use in immunocompromised patients
        map.insert("immunodeficiency_prophylactic_drug_multiplier".to_string(), 8.0);  // 8x higher drug initiation rate for immunocompromised (prophylaxis)
        map.insert("antibiotic_infection_prevention_efficacy".to_string(), 0.7);       // 70% efficacy: allow more breakthrough infections despite prophylaxis


        // Sepsis Mortality Parameters (Age, Region, and Risk Factor dependent)
        map.insert("base_sepsis_death_risk_per_day".to_string(), 0.003); // Base 1% daily death risk for sepsis
        map.insert("sepsis_age_mortality_multiplier_infant".to_string(), 3.0); // 0-1 years: much higher risk
        map.insert("sepsis_age_mortality_multiplier_child".to_string(), 0.5); // 1-18 years: lower risk
        map.insert("sepsis_age_mortality_multiplier_adult".to_string(), 1.0); // 18-65 years: baseline risk
        map.insert("sepsis_age_mortality_multiplier_elderly".to_string(), 2.5); // 65+ years: much higher risk
        // CALIBRATION: Reduced from 30.0 to 15.0 - immunosuppressed multiplier was very high
        map.insert("sepsis_immunosuppressed_multiplier".to_string(), 15.0); // Immunosuppressed: 15x higher risk

        // Region-specific sepsis mortality multipliers (reflecting healthcare quality)
        map.insert("north_america_sepsis_mortality_multiplier".to_string(), 0.8); // Better ICU care
        map.insert("europe_sepsis_mortality_multiplier".to_string(), 0.7); // Excellent healthcare systems
        map.insert("oceania_sepsis_mortality_multiplier".to_string(), 0.8); // Good healthcare
        map.insert("asia_sepsis_mortality_multiplier".to_string(), 1.2); // Variable healthcare quality
        map.insert("south_america_sepsis_mortality_multiplier".to_string(), 1.4); // Limited ICU access
        map.insert("africa_sepsis_mortality_multiplier".to_string(), 2.0); // Limited healthcare infrastructure

        // Sepsis Recovery Parameters (Logistic Model)
        map.insert("sepsis_base_log_odds_of_recovery_per_day".to_string(), -1.5); // Base log odds (low baseline recovery probability ~12%)
        map.insert("sepsis_log_odds_bacteria_level".to_string(), -0.3); // Higher bacteria level decreases recovery (negative coefficient)
        map.insert("sepsis_log_odds_in_hospital".to_string(), 0.8); // Being in hospital increases recovery probability
        map.insert("sepsis_log_odds_age_infant".to_string(), -0.5); // Infants have lower recovery probability
        map.insert("sepsis_log_odds_age_child".to_string(), 0.4); // Children have higher recovery probability
        map.insert("sepsis_log_odds_age_adult".to_string(), 0.0); // Adults baseline (reference category)
        map.insert("sepsis_log_odds_age_elderly".to_string(), -0.7); // Elderly have lower recovery probability
        map.insert("sepsis_log_odds_immunosuppressed".to_string(), -1.0); // Immunosuppressed have much lower recovery probability

        // Region-specific sepsis recovery log odds (reflecting healthcare quality and ICU availability)
        map.insert("sepsis_log_odds_region_north_america".to_string(), 0.4); // Better healthcare systems increase recovery
        map.insert("sepsis_log_odds_region_europe".to_string(), 0.5); // Excellent healthcare systems, best recovery rates
        map.insert("sepsis_log_odds_region_oceania".to_string(), 0.3); // Good healthcare systems
        map.insert("sepsis_log_odds_region_asia".to_string(), 0.0); // Mixed healthcare quality, reference category
        map.insert("sepsis_log_odds_region_south_america".to_string(), -0.3); // Limited ICU access decreases recovery
        map.insert("sepsis_log_odds_region_africa".to_string(), -0.7); // Limited healthcare infrastructure significantly decreases recovery
        map.insert("sepsis_log_odds_region_home".to_string(), 0.0); // Default to reference category

        map.insert("sepsis_minimum_duration_days".to_string(), 1.0); // Minimum sepsis duration (1 day)

        //  Default Toxicity Parameter
        //  Default Microbiome Clearance Parameter (required by simulation logic)
        map.insert("default_microbiome_clearance_probability_per_day".to_string(), 0.01); // E.g., 1% chance to lose carriage per day
        // Probability of clearing microbiome when drug treatment successfully clears infection
        map.insert("microbiome_clearance_probability_on_drug_treatment".to_string(), 0.8); // 80% chance drugs clear microbiome when they clear infection

        // ===========================================================================================
        // --- Enhanced Microbiome/Carriage Model Parameters ---
        // ===========================================================================================
        // These parameters implement a biologically realistic model of bacterial carriage (asymptomatic
        // colonization) and its critical role in antimicrobial resistance dynamics. Carriage matters
        // because: (1) it's 10-100x more common than infection, (2) carriers are the primary reservoir
        // for resistance transmission, and (3) when carriers develop infections, they overwhelmingly
        // inherit their carried strain's resistance profile (carrier amplification effect).

        // --- ANTIBIOTIC DISRUPTION EFFECT ON CARRIAGE ACQUISITION ---
        // Mechanism: Antibiotics kill commensal bacteria, disrupting colonization resistance and creating
        // ecological niches for pathogen colonization. This is why C. difficile infections spike during
        // broad-spectrum antibiotic use, and why MRSA/ESBL colonization increases during antibiotic courses.
        // Empirical basis: 5-15x increased colonization risk during antibiotic therapy, persisting weeks
        // to months after cessation. Studies show antibiotics are the strongest risk factor for MDR carriage.
        map.insert("default_microbiome_disruption_log_odds".to_string(), 0.3);
        map.insert("microbiome_resistance_multiplier_on_acquisition".to_string(), 0.35);
        map.insert("infection_from_microbiome_dampening".to_string(), 0.85);
        // Each active antibiotic adds +0.3 to log-odds of carriage acquisition (multiplicative ~1.35x per drug)
        // Default 0.3 gives ~2x risk with 2 drugs, ~3x with 3 drugs (reasonable based on literature)

        map.insert("antibiotic_disruption_decay_half_life_days".to_string(), 30.0);
        // Half-life for decay of disruption effect after antibiotics stop (reserved for future implementation)
        // Empirical basis: Microbiome recovery takes weeks to months; colonization risk elevated for 1-3 months post-antibiotics

        // --- DURATION EFFECTS ON CARRIAGE CLEARANCE ---
        // Mechanism: Newly acquired bacteria are susceptible to immune clearance and microbial competition.
        // Over time, successful colonizers establish stable niches, form biofilms, and evade immunity,
        // becoming progressively harder to eliminate (established vs. transient colonization).
        // Empirical basis: MRSA decolonization success ~70% for recent carriers vs ~30% for chronic carriers.
        // S. aureus carriage often persists months to years once established.
        map.insert("carriage_duration_log_odds_coefficient".to_string(), -0.01);
        // NEGATIVE coefficient: each day of carriage reduces clearance log-odds by 0.01
        // At 100 days: -1.0 log-odds reduction (clearance ~2.7x less likely)
        // At 200 days: -2.0 log-odds reduction (clearance ~7.4x less likely, hits cap)

        map.insert("carriage_duration_max_log_odds_effect".to_string(), -2.0);
        // Cap duration effect at -2.0 log-odds (~7.4x reduction in clearance probability)
        // Prevents unrealistic complete persistence; even chronic carriers can occasionally clear colonization

        // --- ANTIBIOTIC EFFECT ON CARRIAGE CLEARANCE ---
        // Mechanism: Antibiotics with activity against colonizing bacteria suppress or eliminate them,
        // even at sub-therapeutic concentrations. This is why prophylaxis prevents colonization and
        // why treatment courses often clear carriage as a side effect.
        // Empirical basis: Decolonization protocols use topical antibiotics (mupirocin for MRSA nasal carriage).
        // Systemic treatment often clears S. aureus carriage. However, resistant strains persist.
        map.insert("antibiotic_clearance_log_odds_per_unit_activity".to_string(), 0.5);
        // Each unit of activity_r adds +0.5 to clearance log-odds
        // activity_r already accounts for drug level, potency, and resistance, so this scales appropriately
        // Example: activity_r=2.0 → +1.0 log-odds → ~2.7x higher clearance probability

        // --- CARRIER RESISTANCE INHERITANCE (CARRIER AMPLIFICATION EFFECT) ---
        // Mechanism: When carriers develop infections, the infecting strain is usually their carried strain
        // (endogenous infection), inheriting its resistance profile. This creates a "carrier amplification
        // effect" where resistance rates in infections exceed population prevalence.
        // Empirical basis:
        // - MRSA carriers: 80-90% of S. aureus infections are MRSA (vs ~30% in non-carriers)
        // - ESBL-E. coli carriers: 70-80% of UTIs are ESBL-positive (vs ~10-15% in non-carriers)
        // - VRE carriers: >90% of subsequent bacteremias are VRE
        // Population impact: Carriers maintain resistance without selective pressure (asymptomatic), then
        // amplify resistance rates when they develop infections. This is THE key mechanism for resistance
        // spread in populations, more important than de novo emergence during treatment.
        map.insert("carrier_resistance_inheritance_probability".to_string(), 0.55);
        map.insert(
            "majority_r_memory_retention_per_day".to_string(),
            0.93,
        );
        // Majority_r cache defaults: rolling window horizon and minimum sample threshold.
        map.insert("majority_r_window_days".to_string(), 1000.0);
        map.insert("majority_r_min_total_samples".to_string(), 10.0);
        // Prevent small simulations from catastrophically erasing resistance prevalence once observed;
        // flip to 0 if you want buckets to decay back to zero when no positive samples remain.
        map.insert(
            "majority_r_freeze_at_last_positive".to_string(),
            0.0,
        );
        // 55% probability that carrier's infection inherits microbiome resistance profile
        // Default 0.55 keeps endogenous infections common without locking in microbiome resistance
        // This parameter has MASSIVE impact on population resistance dynamics - most important in the model

        map.insert("default_drug_toxicity_death_hazard_per_unit_level".to_string(), 0.0); // Most drugs have negligible fatal toxicity risk per unit level by default
        map.insert("default_toxicity_reservoir_half_life_days".to_string(), 1.5); // Toxicity hazard persistence when drug stops (days)
        map.insert("toxicity_age_multiplier_infant".to_string(), 1.8); // Neonates more vulnerable to severe toxicity
        map.insert("toxicity_age_multiplier_child".to_string(), 1.2);
        map.insert("toxicity_age_multiplier_adult".to_string(), 1.0);
        map.insert("toxicity_age_multiplier_elderly".to_string(), 2.2);
        map.insert("toxicity_immunosuppressed_multiplier".to_string(), 2.5);
        map.insert("toxicity_hospital_multiplier".to_string(), 1.3); // Hospitalized patients often have greater monitoring but also more severe illness

        // --- Age Category Effects on Infection Acquisition ---
        // Default age category log-odds adjustments (applied to all bacteria unless overridden)
        // Age categories: infant (0-1y), preschool (1-5y), school (5-18y), young_adult (18-50y), middle_age (50-70y), elderly (70+y)
        map.insert("default_log_odds_infant".to_string(), 1.5);      // Infants: higher susceptibility due to immature immune system
        map.insert("default_log_odds_preschool".to_string(), 0.8);   // Preschoolers: moderately higher susceptibility
        map.insert("default_log_odds_school".to_string(), 0.3);      // School age: slightly higher susceptibility
        map.insert("default_log_odds_young_adult".to_string(), 0.0); // Young adults: reference category
        map.insert("default_log_odds_middle_age".to_string(), 0.2);  // Middle age: slightly higher susceptibility
        map.insert("default_log_odds_elderly".to_string(), 0.9);     // Elderly: highest susceptibility

        // --- Age-Region Interaction Effects ---
        // General age-region interactions (applied when bacteria-specific interactions not available)
        // Format: {region}_log_odds_{age_category}

        // Africa: Higher infectious disease burden across all ages, especially in vulnerable populations
        map.insert("africa_log_odds_infant".to_string(), 2.0);       // Very high infant susceptibility (malnutrition, poor healthcare access)
        map.insert("africa_log_odds_preschool".to_string(), 1.2);    // High preschooler susceptibility
        map.insert("africa_log_odds_school".to_string(), 0.6);       // Moderate school age susceptibility
        map.insert("africa_log_odds_young_adult".to_string(), 0.3);  // Slightly higher young adult susceptibility
        map.insert("africa_log_odds_middle_age".to_string(), 0.4);   // Higher middle age susceptibility
        map.insert("africa_log_odds_elderly".to_string(), 1.3);      // High elderly susceptibility

        // Asia: Variable healthcare quality, high population density effects
        map.insert("asia_log_odds_infant".to_string(), 1.0);         // Moderately high infant susceptibility
        map.insert("asia_log_odds_preschool".to_string(), 0.5);      // Moderate preschooler susceptibility
        map.insert("asia_log_odds_school".to_string(), 0.2);         // Slight school age susceptibility increase
        map.insert("asia_log_odds_young_adult".to_string(), 0.1);    // Slight young adult susceptibility increase
        map.insert("asia_log_odds_middle_age".to_string(), 0.2);     // Slight middle age susceptibility increase
        map.insert("asia_log_odds_elderly".to_string(), 0.8);        // High elderly susceptibility

        // Europe: Generally good healthcare, lower infectious disease burden
        map.insert("europe_log_odds_infant".to_string(), -0.2);      // Slightly lower infant susceptibility
        map.insert("europe_log_odds_preschool".to_string(), -0.1);   // Slightly lower preschooler susceptibility
        map.insert("europe_log_odds_school".to_string(), 0.0);       // Neutral school age
        map.insert("europe_log_odds_young_adult".to_string(), 0.0);  // Neutral young adult
        map.insert("europe_log_odds_middle_age".to_string(), 0.0);   // Neutral middle age
        map.insert("europe_log_odds_elderly".to_string(), 0.2);      // Moderate elderly susceptibility

        // North America: Reference region (all zeros)
        map.insert("north_america_log_odds_infant".to_string(), 0.0);       // Reference
        map.insert("north_america_log_odds_preschool".to_string(), 0.0);    // Reference
        map.insert("north_america_log_odds_school".to_string(), 0.0);       // Reference
        map.insert("north_america_log_odds_young_adult".to_string(), 0.0);  // Reference
        map.insert("north_america_log_odds_middle_age".to_string(), 0.0);   // Reference
        map.insert("north_america_log_odds_elderly".to_string(), 0.0);      // Reference

        // South America: Moderate infectious disease burden, variable healthcare access
        map.insert("south_america_log_odds_infant".to_string(), 1.2);       // High infant susceptibility
        map.insert("south_america_log_odds_preschool".to_string(), 0.7);    // Moderate preschooler susceptibility
        map.insert("south_america_log_odds_school".to_string(), 0.3);       // Slight school age susceptibility increase
        map.insert("south_america_log_odds_young_adult".to_string(), 0.2);  // Slight young adult susceptibility increase
        map.insert("south_america_log_odds_middle_age".to_string(), 0.3);   // Moderate middle age susceptibility increase
        map.insert("south_america_log_odds_elderly".to_string(), 0.9);      // High elderly susceptibility

        // Oceania: Generally good healthcare, similar to North America but smaller healthcare systems
        map.insert("oceania_log_odds_infant".to_string(), 0.1);       // Slightly higher infant susceptibility
        map.insert("oceania_log_odds_preschool".to_string(), 0.0);    // Neutral preschooler
        map.insert("oceania_log_odds_school".to_string(), 0.0);       // Neutral school age
        map.insert("oceania_log_odds_young_adult".to_string(), 0.0);  // Neutral young adult
        map.insert("oceania_log_odds_middle_age".to_string(), 0.0);   // Neutral middle age
        map.insert("oceania_log_odds_elderly".to_string(), 0.3);      // Moderate elderly susceptibility

        // --- Bacteria-Specific Age Category Effects ---
        // These override the default age category effects for specific bacteria
        // Format: {bacteria_clean}_log_odds_{age_category}

        // === SEXUALLY TRANSMITTED BACTERIA ===
        // chlamydia_trachomatis - Peak in sexually active young adults
        map.insert("chlamydia_trachomatis_log_odds_infant".to_string(), -3.0);      // Very low risk (vertical transmission only)
        map.insert("chlamydia_trachomatis_log_odds_preschool".to_string(), -3.5);   // Extremely low risk
        map.insert("chlamydia_trachomatis_log_odds_school".to_string(), -1.0);      // Capture adolescent-onset sexual activity
        map.insert("chlamydia_trachomatis_log_odds_young_adult".to_string(), 1.7);  // VERY HIGH risk - peak sexual activity
        map.insert("chlamydia_trachomatis_log_odds_middle_age".to_string(), 0.8);   // Moderate risk - continued sexual activity
        map.insert("chlamydia_trachomatis_log_odds_elderly".to_string(), -1.0);     // Low risk - reduced sexual activity

        // neisseria_gonorrhoeae - Similar pattern to chlamydia
        map.insert("neisseria_gonorrhoeae_log_odds_infant".to_string(), -2.5);      // Very low risk (vertical transmission)
        map.insert("neisseria_gonorrhoeae_log_odds_preschool".to_string(), -3.5);   // Extremely low risk
        map.insert("neisseria_gonorrhoeae_log_odds_school".to_string(), -0.8);      // Capture rising adolescent incidence
        map.insert("neisseria_gonorrhoeae_log_odds_young_adult".to_string(), 2.0);  // VERY HIGH risk - peak sexual activity
        map.insert("neisseria_gonorrhoeae_log_odds_middle_age".to_string(), 0.9);   // Moderate risk
        map.insert("neisseria_gonorrhoeae_log_odds_elderly".to_string(), -0.8);     // Low risk

        // treponema_pallidum (Syphilis) - Similar pattern to other STIs
        map.insert("treponema_pallidum_log_odds_infant".to_string(), -2.2);         // Low risk (congenital syphilis)
        map.insert("treponema_pallidum_log_odds_preschool".to_string(), -4.3);      // Extremely low risk
        map.insert("treponema_pallidum_log_odds_school".to_string(), -2.4);         // Very low risk
        map.insert("treponema_pallidum_log_odds_young_adult".to_string(), 0.4);     // Elevated but moderated risk
        map.insert("treponema_pallidum_log_odds_middle_age".to_string(), -0.6);     // Slightly above baseline risk
        map.insert("treponema_pallidum_log_odds_elderly".to_string(), -2.2);        // Low risk

        // === RESPIRATORY BACTERIA ===
        // streptococcus_pneumoniae - Classic U-shaped age distribution (high in infants and elderly)
        map.insert("streptococcus_pneumoniae_log_odds_infant".to_string(), 1.7);         // Elevated risk - immature immunity
        map.insert("streptococcus_pneumoniae_log_odds_preschool".to_string(), 0.9);      // Elevated risk - daycare exposure
        map.insert("streptococcus_pneumoniae_log_odds_school".to_string(), 0.2);         // Moderate risk - school exposure
        map.insert("streptococcus_pneumoniae_log_odds_young_adult".to_string(), -0.4);   // Low risk - strong immunity
        map.insert("streptococcus_pneumoniae_log_odds_middle_age".to_string(), -0.1);    // Slightly below baseline risk
        map.insert("streptococcus_pneumoniae_log_odds_elderly".to_string(), 1.2);        // HIGH risk - immunosenescence

        // haemophilus_influenzae - Similar to pneumococcus but more pediatric
        map.insert("haemophilus_influenzae_log_odds_infant".to_string(), 2.5);           // VERY HIGH risk - major infant pathogen
        map.insert("haemophilus_influenzae_log_odds_preschool".to_string(), 1.5);        // HIGH risk
        map.insert("haemophilus_influenzae_log_odds_school".to_string(), 0.8);           // Moderate risk
        map.insert("haemophilus_influenzae_log_odds_young_adult".to_string(), -0.5);     // Low risk
        map.insert("haemophilus_influenzae_log_odds_middle_age".to_string(), -0.2);      // Low risk
        map.insert("haemophilus_influenzae_log_odds_elderly".to_string(), 1.0);          // High risk (COPD patients)

        // moraxella_catarrhalis - Primarily pediatric and elderly (COPD)
        map.insert("moraxella_catarrhalis_log_odds_infant".to_string(), 1.8);            // VERY HIGH risk
        map.insert("moraxella_catarrhalis_log_odds_preschool".to_string(), 1.0);         // HIGH risk
        map.insert("moraxella_catarrhalis_log_odds_school".to_string(), 0.3);            // Moderate risk
        map.insert("moraxella_catarrhalis_log_odds_young_adult".to_string(), -0.8);      // Low risk
        map.insert("moraxella_catarrhalis_log_odds_middle_age".to_string(), -0.3);       // Low risk
        map.insert("moraxella_catarrhalis_log_odds_elderly".to_string(), 1.2);           // HIGH risk (COPD)

        // === ENTERIC/FOODBORNE BACTERIA ===
        // Salmonella species - Higher in children and elderly
        map.insert("salmonella_enterica_serovar_typhi_log_odds_infant".to_string(), 1.0);      // High risk - severe disease
        map.insert("salmonella_enterica_serovar_typhi_log_odds_preschool".to_string(), 0.8);   // High risk
        map.insert("salmonella_enterica_serovar_typhi_log_odds_school".to_string(), 0.5);      // Moderate risk
        map.insert("salmonella_enterica_serovar_typhi_log_odds_young_adult".to_string(), 0.0); // Baseline (travel-related)
        map.insert("salmonella_enterica_serovar_typhi_log_odds_middle_age".to_string(), 0.2);  // Moderate risk
        map.insert("salmonella_enterica_serovar_typhi_log_odds_elderly".to_string(), 0.8);     // High risk - severe disease

        // shigella_spp. - Higher in children (daycare transmission)
        map.insert("shigella_spp._log_odds_infant".to_string(), 1.5);            // VERY HIGH risk
        map.insert("shigella_spp._log_odds_preschool".to_string(), 2.0);         // EXTREMELY HIGH risk - daycare outbreaks
        map.insert("shigella_spp._log_odds_school".to_string(), 1.2);            // HIGH risk - school transmission
        map.insert("shigella_spp._log_odds_young_adult".to_string(), 0.3);       // Moderate risk
        map.insert("shigella_spp._log_odds_middle_age".to_string(), 0.0);        // Baseline risk
        map.insert("shigella_spp._log_odds_elderly".to_string(), 0.5);           // Moderate risk

        // campylobacter_jejuni - All ages but higher in young children
        map.insert("campylobacter_jejuni_log_odds_infant".to_string(), 1.6);          // VERY HIGH risk - severe dehydration
        map.insert("campylobacter_jejuni_log_odds_preschool".to_string(), 1.2);       // HIGH risk
        map.insert("campylobacter_jejuni_log_odds_school".to_string(), 0.4);          // Moderate risk
        map.insert("campylobacter_jejuni_log_odds_young_adult".to_string(), 0.2);     // Slightly elevated risk in young travelers
        map.insert("campylobacter_jejuni_log_odds_middle_age".to_string(), 0.0);      // Baseline risk
        map.insert("campylobacter_jejuni_log_odds_elderly".to_string(), 0.5);         // Moderate risk

        // helicobacter_pylori - Accumulates with age due to persistent colonization
        map.insert("helicobacter_pylori_log_odds_infant".to_string(), -2.5);         // Very low risk - rare vertical transmission
        map.insert("helicobacter_pylori_log_odds_preschool".to_string(), -1.5);      // Low risk
        map.insert("helicobacter_pylori_log_odds_school".to_string(), -0.5);         // Slowly rising risk
        map.insert("helicobacter_pylori_log_odds_young_adult".to_string(), 0.8);     // Moderate risk - cumulative exposure
        map.insert("helicobacter_pylori_log_odds_middle_age".to_string(), 1.4);      // High risk - high prevalence
        map.insert("helicobacter_pylori_log_odds_elderly".to_string(), 1.8);         // VERY HIGH risk - chronic colonization

        // === HEALTHCARE-ASSOCIATED BACTERIA ===
        // acinetobacter_baumannii - Higher in critically ill (middle-aged to elderly)
        map.insert("acinetobacter_baumannii_log_odds_infant".to_string(), 0.5);       // Moderate risk (NICU)
        map.insert("acinetobacter_baumannii_log_odds_preschool".to_string(), -0.5);   // Low risk
        map.insert("acinetobacter_baumannii_log_odds_school".to_string(), -0.8);      // Low risk
        map.insert("acinetobacter_baumannii_log_odds_young_adult".to_string(), 0.2);  // Moderate risk (trauma patients)
        map.insert("acinetobacter_baumannii_log_odds_middle_age".to_string(), 0.8);   // HIGH risk (ICU patients)
        map.insert("acinetobacter_baumannii_log_odds_elderly".to_string(), 1.5);      // VERY HIGH risk

        // clostridioides_difficile - Strongly age-associated (elderly)
        map.insert("clostridioides_difficile_log_odds_infant".to_string(), -1.0);     // Low risk (protective microbiome)
        map.insert("clostridioides_difficile_log_odds_preschool".to_string(), -1.5);  // Very low risk
        map.insert("clostridioides_difficile_log_odds_school".to_string(), -2.0);     // Very low risk
        map.insert("clostridioides_difficile_log_odds_young_adult".to_string(), -0.5); // Low risk
        map.insert("clostridioides_difficile_log_odds_middle_age".to_string(), 0.5);  // Moderate risk
        map.insert("clostridioides_difficile_log_odds_elderly".to_string(), 2.0);     // VERY HIGH risk - major pathogen

        // === UROGENITAL BACTERIA ===
        // escherichia_coli (UTI) - Higher in infants, young women, and elderly
        map.insert("escherichia_coli_log_odds_infant".to_string(), 1.2);             // HIGH risk - anatomical factors
        map.insert("escherichia_coli_log_odds_preschool".to_string(), 0.3);          // Moderate risk
        map.insert("escherichia_coli_log_odds_school".to_string(), 0.0);             // Baseline risk
        map.insert("escherichia_coli_log_odds_young_adult".to_string(), 0.8);        // HIGH risk - sexual activity, pregnancy
        map.insert("escherichia_coli_log_odds_middle_age".to_string(), 0.6);         // High risk - continued risk factors
        map.insert("escherichia_coli_log_odds_elderly".to_string(), 1.5);            // VERY HIGH risk - multiple factors

        // === INVASIVE/SYSTEMIC BACTERIA ===
        // neisseria_meningitidis - Bimodal distribution (infants and adolescents/young adults)
        map.insert("neisseria_meningitidis_log_odds_infant".to_string(), 1.8);           // VERY HIGH risk but reduced after vaccine rollout
        map.insert("neisseria_meningitidis_log_odds_preschool".to_string(), 0.4);        // Elevated risk
        map.insert("neisseria_meningitidis_log_odds_school".to_string(), 1.0);           // HIGH risk - school clustering
        map.insert("neisseria_meningitidis_log_odds_young_adult".to_string(), 1.3);      // VERY HIGH risk - dormitories, military
        map.insert("neisseria_meningitidis_log_odds_middle_age".to_string(), -0.2);      // Slightly below baseline
        map.insert("neisseria_meningitidis_log_odds_elderly".to_string(), 0.2);          // Modest rebound risk

        // listeria_monocytogenes - Mainly immunocompromised, pregnant women, elderly
        map.insert("listeria_monocytogenes_log_odds_infant".to_string(), 1.8);          // VERY HIGH risk (neonatal)
        map.insert("listeria_monocytogenes_log_odds_preschool".to_string(), -0.5);      // Low risk
        map.insert("listeria_monocytogenes_log_odds_school".to_string(), -1.0);         // Very low risk
        map.insert("listeria_monocytogenes_log_odds_young_adult".to_string(), 0.5);     // Moderate risk (pregnancy)
        map.insert("listeria_monocytogenes_log_odds_middle_age".to_string(), 0.0);      // Baseline risk
        map.insert("listeria_monocytogenes_log_odds_elderly".to_string(), 1.5);         // VERY HIGH risk

        // --- Bacteria-Specific Age-Region Interaction Overrides ---
        // These override the general age-region interactions for specific bacteria where there's strong evidence
        // Format: {bacteria_clean}_{region}_log_odds_{age_category}

        // === shigella_spp. - STRONG AGE-REGION INTERACTIONS ===
        // Much higher burden in developing regions, especially in children due to poor sanitation

        // Africa - Very high burden, especially in children (poor sanitation, malnutrition)
        map.insert("shigella_spp._africa_log_odds_infant".to_string(), 3.0);        // EXTREMELY HIGH - severe dehydration risk
        map.insert("shigella_spp._africa_log_odds_preschool".to_string(), 3.5);     // HIGHEST RISK - daycare equivalent settings
        map.insert("shigella_spp._africa_log_odds_school".to_string(), 2.8);        // VERY HIGH - school crowding
        map.insert("shigella_spp._africa_log_odds_young_adult".to_string(), 1.0);   // HIGH - caregivers
        map.insert("shigella_spp._africa_log_odds_middle_age".to_string(), 0.8);    // MODERATE-HIGH
        map.insert("shigella_spp._africa_log_odds_elderly".to_string(), 1.2);       // HIGH - vulnerability

        // Asia - High burden in children, variable by subregion
        map.insert("shigella_spp._asia_log_odds_infant".to_string(), 2.5);          // VERY HIGH
        map.insert("shigella_spp._asia_log_odds_preschool".to_string(), 3.0);       // EXTREMELY HIGH
        map.insert("shigella_spp._asia_log_odds_school".to_string(), 2.2);          // VERY HIGH
        map.insert("shigella_spp._asia_log_odds_young_adult".to_string(), 0.8);     // HIGH
        map.insert("shigella_spp._asia_log_odds_middle_age".to_string(), 0.5);      // MODERATE
        map.insert("shigella_spp._asia_log_odds_elderly".to_string(), 0.8);         // HIGH

        // Europe - Much lower burden, mainly travel-related
        map.insert("shigella_spp._europe_log_odds_infant".to_string(), -0.5);       // Low
        map.insert("shigella_spp._europe_log_odds_preschool".to_string(), 0.2);     // Slight increase (daycare)
        map.insert("shigella_spp._europe_log_odds_school".to_string(), -0.2);       // Low
        map.insert("shigella_spp._europe_log_odds_young_adult".to_string(), 0.5);   // Moderate (travel)
        map.insert("shigella_spp._europe_log_odds_middle_age".to_string(), 0.3);    // Moderate (travel)

        // === SALMONELLA TYPHI - STRONG REGIONAL AND AGE PATTERNS ===
        // Endemic in South Asia, Sub-Saharan Africa; severe in children

        // Africa - High burden, especially severe in children
        map.insert("salmonella_enterica_serovar_typhi_africa_log_odds_infant".to_string(), 2.5);      // EXTREMELY HIGH - very severe
        map.insert("salmonella_enterica_serovar_typhi_africa_log_odds_preschool".to_string(), 2.2);   // VERY HIGH
        map.insert("salmonella_enterica_serovar_typhi_africa_log_odds_school".to_string(), 1.8);      // VERY HIGH
        map.insert("salmonella_enterica_serovar_typhi_africa_log_odds_young_adult".to_string(), 1.2); // HIGH
        map.insert("salmonella_enterica_serovar_typhi_africa_log_odds_middle_age".to_string(), 1.0);  // HIGH
        map.insert("salmonella_enterica_serovar_typhi_africa_log_odds_elderly".to_string(), 1.5);     // VERY HIGH

        // Asia - Endemic regions (South Asia), very high burden
        map.insert("salmonella_enterica_serovar_typhi_asia_log_odds_infant".to_string(), 3.0);        // EXTREMELY HIGH
        map.insert("salmonella_enterica_serovar_typhi_asia_log_odds_preschool".to_string(), 2.8);     // EXTREMELY HIGH
        map.insert("salmonella_enterica_serovar_typhi_asia_log_odds_school".to_string(), 2.5);        // EXTREMELY HIGH
        map.insert("salmonella_enterica_serovar_typhi_asia_log_odds_young_adult".to_string(), 1.8);   // VERY HIGH
        map.insert("salmonella_enterica_serovar_typhi_asia_log_odds_middle_age".to_string(), 1.5);    // VERY HIGH
        map.insert("salmonella_enterica_serovar_typhi_asia_log_odds_elderly".to_string(), 2.0);       // EXTREMELY HIGH

        // Europe/North America - Mainly travel-related, much lower endemic burden
        map.insert("salmonella_enterica_serovar_typhi_europe_log_odds_young_adult".to_string(), 0.8); // Moderate (travel)
        map.insert("salmonella_enterica_serovar_typhi_europe_log_odds_middle_age".to_string(), 0.6);  // Moderate (travel)
        map.insert("salmonella_enterica_serovar_typhi_north_america_log_odds_young_adult".to_string(), 0.5); // Moderate (travel)
        map.insert("salmonella_enterica_serovar_typhi_north_america_log_odds_middle_age".to_string(), 0.3);  // Moderate (travel)

        // === vibrio_cholerae - HIGHLY REGION AND AGE SPECIFIC ===
        // Endemic in specific regions, severe in children and elderly

        // Africa - High burden in certain regions, water-related
        map.insert("vibrio_cholerae_africa_log_odds_infant".to_string(), 2.8);       // EXTREMELY HIGH - severe dehydration
        map.insert("vibrio_cholerae_africa_log_odds_preschool".to_string(), 2.5);    // EXTREMELY HIGH
        map.insert("vibrio_cholerae_africa_log_odds_school".to_string(), 2.0);       // VERY HIGH
        map.insert("vibrio_cholerae_africa_log_odds_young_adult".to_string(), 1.2);  // HIGH
        map.insert("vibrio_cholerae_africa_log_odds_middle_age".to_string(), 1.0);   // HIGH
        map.insert("vibrio_cholerae_africa_log_odds_elderly".to_string(), 2.0);      // EXTREMELY HIGH

        // Asia - Endemic regions, high burden
        map.insert("vibrio_cholerae_asia_log_odds_infant".to_string(), 3.2);         // EXTREMELY HIGH
        map.insert("vibrio_cholerae_asia_log_odds_preschool".to_string(), 2.8);      // EXTREMELY HIGH
        map.insert("vibrio_cholerae_asia_log_odds_school".to_string(), 2.2);         // VERY HIGH
        map.insert("vibrio_cholerae_asia_log_odds_young_adult".to_string(), 1.5);    // VERY HIGH
        map.insert("vibrio_cholerae_asia_log_odds_middle_age".to_string(), 1.3);     // HIGH
        map.insert("vibrio_cholerae_asia_log_odds_elderly".to_string(), 2.2);        // EXTREMELY HIGH

        // Europe/North America - Very rare, mainly travel-related
        map.insert("vibrio_cholerae_europe_log_odds_young_adult".to_string(), -1.0); // Low (travel)
        map.insert("vibrio_cholerae_europe_log_odds_middle_age".to_string(), -0.8);  // Low (travel)
        map.insert("vibrio_cholerae_north_america_log_odds_young_adult".to_string(), -1.2); // Very low
        map.insert("vibrio_cholerae_north_america_log_odds_middle_age".to_string(), -1.0);  // Very low

        // === neisseria_meningitidis - SPECIFIC HIGH-RISK COMBINATIONS ===
        // Higher burden in "meningitis belt" of Africa, specific age patterns

        // Africa - "Meningitis belt", extremely high infant and adolescent/young adult risk
        map.insert("neisseria_meningitidis_africa_log_odds_infant".to_string(), 2.5);      // VERY HIGH
        map.insert("neisseria_meningitidis_africa_log_odds_preschool".to_string(), 1.0);   // HIGH
        map.insert("neisseria_meningitidis_africa_log_odds_school".to_string(), 1.5);      // VERY HIGH
        map.insert("neisseria_meningitidis_africa_log_odds_young_adult".to_string(), 1.9); // EXTREMELY HIGH - "meningitis belt"
        map.insert("neisseria_meningitidis_africa_log_odds_middle_age".to_string(), 0.3);  // Moderate
        map.insert("neisseria_meningitidis_africa_log_odds_elderly".to_string(), 0.6);     // HIGH

        // Europe/North America - Much lower overall burden, but maintain age pattern
        map.insert("neisseria_meningitidis_europe_log_odds_infant".to_string(), 1.2);      // Reduced after vaccine impact
        map.insert("neisseria_meningitidis_europe_log_odds_young_adult".to_string(), 0.6); // Moderate (university outbreaks)
        map.insert("neisseria_meningitidis_north_america_log_odds_infant".to_string(), 1.4);     // Lower due to MenACWY uptake
        map.insert("neisseria_meningitidis_north_america_log_odds_young_adult".to_string(), 0.5); // Moderate (dormitories)

        // === haemophilus_influenzae - VACCINATION IMPACT VARIES BY REGION ===
        // Lower burden in regions with good Hib vaccination coverage

        // Africa - Higher burden due to limited vaccination coverage
        map.insert("haemophilus_influenzae_africa_log_odds_infant".to_string(), 3.5);      // EXTREMELY HIGH
        map.insert("haemophilus_influenzae_africa_log_odds_preschool".to_string(), 2.5);   // EXTREMELY HIGH
        map.insert("haemophilus_influenzae_africa_log_odds_school".to_string(), 1.5);      // VERY HIGH

        // Asia - Variable vaccination coverage
        map.insert("haemophilus_influenzae_asia_log_odds_infant".to_string(), 3.0);        // EXTREMELY HIGH
        map.insert("haemophilus_influenzae_asia_log_odds_preschool".to_string(), 2.0);     // VERY HIGH
        map.insert("haemophilus_influenzae_asia_log_odds_school".to_string(), 1.2);        // HIGH

        // Europe - Good vaccination coverage, lower pediatric burden
        map.insert("haemophilus_influenzae_europe_log_odds_infant".to_string(), 1.5);      // Still elevated but lower
        map.insert("haemophilus_influenzae_europe_log_odds_preschool".to_string(), 0.8);   // Moderate
        map.insert("haemophilus_influenzae_europe_log_odds_school".to_string(), 0.3);      // Lower

        // North America - Excellent vaccination coverage
        map.insert("haemophilus_influenzae_north_america_log_odds_infant".to_string(), 1.2); // Lower due to vaccination
        map.insert("haemophilus_influenzae_north_america_log_odds_preschool".to_string(), 0.5); // Much lower
        map.insert("haemophilus_influenzae_north_america_log_odds_school".to_string(), 0.1);    // Very low

        // === [J] Regional availability & introduction timelines ===
    // Defines which drugs are reachable in each region and when they enter the market.
    // Scenario templates can override granular entries to stage roll-outs or supply shocks.
    // ---------------- 10) Regional drug availability & introduction timing ----------------
        // Region-specific drug availability multipliers
        // Format: "{region}_drug_{drug_name}_availability"
        // Values: 1.0 = fully available, 0.5 = limited availability, 0.0 = not available
        // Based on realistic antibiotic access patterns across different healthcare systems

        // North America - Full access to most antibiotics
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("north_america_drug_{}_availability", drug), 1.0);
        }

        // Europe - Full access to most antibiotics
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("europe_drug_{}_availability", drug), 1.0);
        }

        // Asia - Good access to most drugs, some newer drugs may be limited
        for &drug in DRUG_SHORT_NAMES.iter() {
            let availability = match drug {
                // Newer, expensive drugs may have limited availability
                "tedizolid" | "ceftaroline" => 0.3,
                "teicoplanin" => 0.7, // More available in Asia than tedizolid
                // TODO: set dalbavancin availability for Asia once launch data added
                _ => 1.0, // Most other drugs widely available
            };
            map.insert(format!("asia_drug_{}_availability", drug), availability);
        }

        // Oceania - Generally good access, similar to developed regions
        for &drug in DRUG_SHORT_NAMES.iter() {
            let availability = match drug {
                "tedizolid" | "ceftaroline" => 0.5, // Somewhat limited
                // TODO: set dalbavancin availability for Oceania when introduced
                _ => 1.0,
            };
            map.insert(format!("oceania_drug_{}_availability", drug), availability);
        }

        // South America - Variable access, newer/expensive drugs limited
    for &drug in DRUG_SHORT_NAMES.iter() {
            let availability = match drug {
                // Very limited access to newest drugs
        "tedizolid" | "ceftaroline" => 0.1,
                "teicoplanin" => 0.3,
                "linezolid" => 0.5,
        // TODO: configure dalbavancin availability for South America
                // Limited access to some carbapenems
                "ertapenem" => 0.6,
                "meropenem" | "imipenem_c" => 0.7,
                // Moderate access to some newer cephalosporins
                "cefepime" => 0.8,
                // Good access to older, established drugs
                _ => 1.0,
            };
            map.insert(format!("south_america_drug_{}_availability", drug), availability);
        }

        // Africa - Most limited access, mainly basic antibiotics available
        for &drug in DRUG_SHORT_NAMES.iter() {
            let availability = match drug {
                // Basic penicillins - widely available
                "penicilling" | "ampicillin" | "amoxicillin" => 1.0,
                // Basic cephalosporins - good availability
                "cephalexin" | "cefazolin" => 0.9,
                "cefuroxime" => 0.7,
                // Third-generation cephalosporins - limited
                "ceftriaxone" => 0.6,
                "ceftazidime" => 0.4,
                // Basic macrolides and fluoroquinolones - moderate availability
                "erythromycin" | "azithromycin" => 0.8,
                "ciprofloxacin" => 0.7,
                "levofloxacin" => 0.5,
                // Aminoglycosides - basic ones available
                "gentamicin" => 0.8,
                "tobramycin" | "amikacin" => 0.4,
                // Older drugs - generally available
                "tetracycline" | "doxycycline" => 0.9,
                "trim_sulf" => 0.9,
                "chlorampheni" => 0.8,
                "metronidazole" => 0.9,
                // Vancomycin - very limited
                "vancomycin" => 0.3,
                // Newer/expensive drugs - very limited or unavailable
                "meropenem" | "imipenem_c" => 0.2,
                "ertapenem" => 0.1,
                "linezolid" => 0.1,
                "tedizolid" | "ceftaroline" | "teicoplanin" => 0.0,
                // TODO: define dalbavancin availability for Africa once rollout planned
                "aztreonam" => 0.1,
                "cefepime" => 0.3,
                "moxifloxacin" => 0.2,
                "minocycline" => 0.4,
                "quinu_dalfo" => 0.1,
                "nitrofurantoin" => 0.6,
                "retapamulin" | "fusidic_a" => 0.2,
                "furazolidone" => 0.3,
                // Default for any remaining drugs
                _ => 0.1,
            };
            map.insert(format!("africa_drug_{}_availability", drug), availability);
        }

        // Home region - use availability based on region_living
        // This will be handled in the drug initiation logic
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("home_drug_{}_availability", drug), 1.0); // Default fallback
        }

        // Ensure you have multipliers for all variants of your `Region` enum,
        // or add a default handling in the `mod.rs` if a region param isn't found.
        // If `Region::Home` refers to a generic home location not tied to a specific geographical region,
        // you might need to reconsider its role or default it to 1.0 or an average.

        // === [K] Demographic distribution defaults ===
        // Population share assumptions across regions and age bands. Adjust when calibrating
        // against census projections so disease burden and treatment demand scale correctly.
        // ---------------- 11) Demographic distribution defaults ----------------
        // Demographic distribution parameters (108 total: 6 regions × 18 age bands)
        // Each parameter represents probability of being in that region-age combination
        // All 108 parameters should sum to 1.0
        // Age bands span from -40000 to +32000 in 4000-year intervals

        // Asia demographic distribution (18 age bands of 4000 years each)
        // 55% of global population based on 1930-2035 timeframe
        // Massive growth: ~90% at negative ages (future births), only ~10% alive in 1930
        map.insert("demo_asia_age_neg40000_neg36000".to_string(), 0.120); // Heavy weighting on future births
        map.insert("demo_asia_age_neg36000_neg32000".to_string(), 0.110);
        map.insert("demo_asia_age_neg32000_neg28000".to_string(), 0.100);
        map.insert("demo_asia_age_neg28000_neg24000".to_string(), 0.090);
        map.insert("demo_asia_age_neg24000_neg20000".to_string(), 0.080);
        map.insert("demo_asia_age_neg20000_neg16000".to_string(), 0.070);
        map.insert("demo_asia_age_neg16000_neg12000".to_string(), 0.060);
        map.insert("demo_asia_age_neg12000_neg8000".to_string(), 0.050);
        map.insert("demo_asia_age_neg8000_neg4000".to_string(), 0.040);
        map.insert("demo_asia_age_neg4000_0".to_string(), 0.030);          // Future births tapering
        map.insert("demo_asia_age_0_4000".to_string(), 0.008);             // Very small portion alive in 1930
        map.insert("demo_asia_age_4000_8000".to_string(), 0.006);
        map.insert("demo_asia_age_8000_12000".to_string(), 0.005);
        map.insert("demo_asia_age_12000_16000".to_string(), 0.004);
        map.insert("demo_asia_age_16000_20000".to_string(), 0.003);
        map.insert("demo_asia_age_20000_24000".to_string(), 0.002);
        map.insert("demo_asia_age_24000_28000".to_string(), 0.001);
        map.insert("demo_asia_age_28000_32000".to_string(), 0.001);

        // Africa demographic distribution
        // 12% of global population based on 1930-2035 timeframe
        // Explosive growth: ~90% at negative ages (future births), only ~10% alive in 1930
        map.insert("demo_africa_age_neg40000_neg36000".to_string(), 0.026); // Heavy weighting on future births
        map.insert("demo_africa_age_neg36000_neg32000".to_string(), 0.024);
        map.insert("demo_africa_age_neg32000_neg28000".to_string(), 0.022);
        map.insert("demo_africa_age_neg28000_neg24000".to_string(), 0.020);
        map.insert("demo_africa_age_neg24000_neg20000".to_string(), 0.018);
        map.insert("demo_africa_age_neg20000_neg16000".to_string(), 0.016);
        map.insert("demo_africa_age_neg16000_neg12000".to_string(), 0.014);
        map.insert("demo_africa_age_neg12000_neg8000".to_string(), 0.012);
        map.insert("demo_africa_age_neg8000_neg4000".to_string(), 0.010);
        map.insert("demo_africa_age_neg4000_0".to_string(), 0.008);          // Future births tapering
        map.insert("demo_africa_age_0_4000".to_string(), 0.002);             // Very small portion alive in 1930
        map.insert("demo_africa_age_4000_8000".to_string(), 0.002);
        map.insert("demo_africa_age_8000_12000".to_string(), 0.001);
        map.insert("demo_africa_age_12000_16000".to_string(), 0.001);
        map.insert("demo_africa_age_16000_20000".to_string(), 0.001);
        map.insert("demo_africa_age_20000_24000".to_string(), 0.001);
        map.insert("demo_africa_age_24000_28000".to_string(), 0.001);
        map.insert("demo_africa_age_28000_32000".to_string(), 0.001);

        // Europe demographic distribution
        // 15% of global population based on 1930-2035 timeframe
        // Moderate growth: ~70% at negative ages (future births), ~30% alive in 1930
        map.insert("demo_europe_age_neg40000_neg36000".to_string(), 0.020); // Moderate weighting on future births
        map.insert("demo_europe_age_neg36000_neg32000".to_string(), 0.018);
        map.insert("demo_europe_age_neg32000_neg28000".to_string(), 0.016);
        map.insert("demo_europe_age_neg28000_neg24000".to_string(), 0.015);
        map.insert("demo_europe_age_neg24000_neg20000".to_string(), 0.014);
        map.insert("demo_europe_age_neg20000_neg16000".to_string(), 0.013);
        map.insert("demo_europe_age_neg16000_neg12000".to_string(), 0.012);
        map.insert("demo_europe_age_neg12000_neg8000".to_string(), 0.011);
        map.insert("demo_europe_age_neg8000_neg4000".to_string(), 0.010);
        map.insert("demo_europe_age_neg4000_0".to_string(), 0.009);          // Future births tapering
        map.insert("demo_europe_age_0_4000".to_string(), 0.008);             // Larger portion alive in 1930
        map.insert("demo_europe_age_4000_8000".to_string(), 0.007);
        map.insert("demo_europe_age_8000_12000".to_string(), 0.006);
        map.insert("demo_europe_age_12000_16000".to_string(), 0.005);
        map.insert("demo_europe_age_16000_20000".to_string(), 0.004);
        map.insert("demo_europe_age_20000_24000".to_string(), 0.003);
        map.insert("demo_europe_age_24000_28000".to_string(), 0.002);
        map.insert("demo_europe_age_28000_32000".to_string(), 0.002);

        // North America demographic distribution
        // 9% of global population based on 1930-2035 timeframe
        // Significant growth: ~75% at negative ages (future births), ~25% alive in 1930
        map.insert("demo_north_america_age_neg40000_neg36000".to_string(), 0.012);
        map.insert("demo_north_america_age_neg36000_neg32000".to_string(), 0.011);
        map.insert("demo_north_america_age_neg32000_neg28000".to_string(), 0.010);
        map.insert("demo_north_america_age_neg28000_neg24000".to_string(), 0.009);
        map.insert("demo_north_america_age_neg24000_neg20000".to_string(), 0.008);
        map.insert("demo_north_america_age_neg20000_neg16000".to_string(), 0.007);
        map.insert("demo_north_america_age_neg16000_neg12000".to_string(), 0.006);
        map.insert("demo_north_america_age_neg12000_neg8000".to_string(), 0.005);
        map.insert("demo_north_america_age_neg8000_neg4000".to_string(), 0.004);
        map.insert("demo_north_america_age_neg4000_0".to_string(), 0.003);
        map.insert("demo_north_america_age_0_4000".to_string(), 0.003);
        map.insert("demo_north_america_age_4000_8000".to_string(), 0.003);
        map.insert("demo_north_america_age_8000_12000".to_string(), 0.002);
        map.insert("demo_north_america_age_12000_16000".to_string(), 0.002);
        map.insert("demo_north_america_age_16000_20000".to_string(), 0.002);
        map.insert("demo_north_america_age_20000_24000".to_string(), 0.002);
        map.insert("demo_north_america_age_24000_28000".to_string(), 0.001);
        map.insert("demo_north_america_age_28000_32000".to_string(), 0.001);

        // South America demographic distribution
        // 6% of global population based on 1930-2035 timeframe
        // Strong growth: ~80% at negative ages (future births), ~20% alive in 1930
        map.insert("demo_south_america_age_neg40000_neg36000".to_string(), 0.010);
        map.insert("demo_south_america_age_neg36000_neg32000".to_string(), 0.009);
        map.insert("demo_south_america_age_neg32000_neg28000".to_string(), 0.008);
        map.insert("demo_south_america_age_neg28000_neg24000".to_string(), 0.007);
        map.insert("demo_south_america_age_neg24000_neg20000".to_string(), 0.006);
        map.insert("demo_south_america_age_neg20000_neg16000".to_string(), 0.005);
        map.insert("demo_south_america_age_neg16000_neg12000".to_string(), 0.004);
        map.insert("demo_south_america_age_neg12000_neg8000".to_string(), 0.004);
        map.insert("demo_south_america_age_neg8000_neg4000".to_string(), 0.003);
        map.insert("demo_south_america_age_neg4000_0".to_string(), 0.002);
        map.insert("demo_south_america_age_0_4000".to_string(), 0.002);
        map.insert("demo_south_america_age_4000_8000".to_string(), 0.002);
        map.insert("demo_south_america_age_8000_12000".to_string(), 0.001);
        map.insert("demo_south_america_age_12000_16000".to_string(), 0.001);
        map.insert("demo_south_america_age_16000_20000".to_string(), 0.001);
        map.insert("demo_south_america_age_20000_24000".to_string(), 0.001);
        map.insert("demo_south_america_age_24000_28000".to_string(), 0.001);
        map.insert("demo_south_america_age_28000_32000".to_string(), 0.001);

        // Oceania demographic distribution
        // 3% of global population based on 1930-2035 timeframe
        // Moderate growth: ~65% at negative ages (future births), ~35% alive in 1930
        map.insert("demo_oceania_age_neg40000_neg36000".to_string(), 0.004);
        map.insert("demo_oceania_age_neg36000_neg32000".to_string(), 0.004);
        map.insert("demo_oceania_age_neg32000_neg28000".to_string(), 0.003);
        map.insert("demo_oceania_age_neg28000_neg24000".to_string(), 0.003);
        map.insert("demo_oceania_age_neg24000_neg20000".to_string(), 0.003);
        map.insert("demo_oceania_age_neg20000_neg16000".to_string(), 0.003);
        map.insert("demo_oceania_age_neg16000_neg12000".to_string(), 0.002);
        map.insert("demo_oceania_age_neg12000_neg8000".to_string(), 0.002);
        map.insert("demo_oceania_age_neg8000_neg4000".to_string(), 0.002);
        map.insert("demo_oceania_age_neg4000_0".to_string(), 0.001);
        map.insert("demo_oceania_age_0_4000".to_string(), 0.002);
        map.insert("demo_oceania_age_4000_8000".to_string(), 0.002);
        map.insert("demo_oceania_age_8000_12000".to_string(), 0.002);
        map.insert("demo_oceania_age_12000_16000".to_string(), 0.001);
        map.insert("demo_oceania_age_16000_20000".to_string(), 0.001);
        map.insert("demo_oceania_age_20000_24000".to_string(), 0.001);
        map.insert("demo_oceania_age_24000_28000".to_string(), 0.001);
        map.insert("demo_oceania_age_28000_32000".to_string(), 0.001);

    apply_potency_overrides_from_embedded_table(&mut map);

    map
    };

    // --- String Parameters (for template names, etc.) ---
    pub static ref STRING_PARAMETERS: HashMap<String, String> = {
        let mut map = HashMap::new();

        // Default age risk templates for all bacteria
        for &bacteria in BACTERIA_LIST.iter() {
            map.insert(format!("{}_age_risk_template", bacteria), "respiratory".to_string()); // Default template
        }

        // Specific bacteria overrides - assign each bacteria to most appropriate template
        map.insert("strep_pneu_age_risk_template".to_string(), "respiratory".to_string());
        map.insert("haem_infl_age_risk_template".to_string(), "respiratory".to_string());
        map.insert("salm_typhi_age_risk_template".to_string(), "gastrointestinal".to_string());
        map.insert("esch_coli_age_risk_template".to_string(), "urogenital".to_string());
        map.insert("pseud_aerug_age_risk_template".to_string(), "bloodstream".to_string());
        map.insert("staph_aureus_age_risk_template".to_string(), "skin_soft_tissue".to_string());
        map.insert("n_gonorrhoeae_age_risk_template".to_string(), "sexually_transmitted".to_string());
        map.insert("acinetobac_bau_age_risk_template".to_string(), "bloodstream".to_string());
        map.insert("staph_epidermidis_age_risk_template".to_string(), "bloodstream".to_string());
        map.insert("staphylococcus_epidermidis_age_risk_template".to_string(), "bloodstream".to_string());
        map.insert("stenotrophomonas_maltophilia_age_risk_template".to_string(), "respiratory".to_string());

        map
    };

    pub static ref PARAMETER_STORE: ParameterStore = ParameterStore::from_parameter_map(&PARAMETERS);
}

// ---------------- 12) Helper lookups (drug intro, availability, etc.) ----------------
/// Helper accessor for the indexed parameter store.
#[allow(dead_code)]
pub fn parameter_store() -> &'static ParameterStore {
    &PARAMETER_STORE
}

/// Retrieves a global simulation parameter.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_global_param(key: &str) -> Option<f64> {
    PARAMETERS.get(key).copied()
}

/// Retrieves a bacteria-specific simulation parameter.
/// It directly looks up "{bacteria_name}_{param_suffix}".
/// Because all bacteria now have explicit entries, there's no need for a "generic_bacteria_" fallback in this function.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_bacteria_param(bacteria_name: &str, param_suffix: &str) -> Option<f64> {
    let canonical = canonicalize_bacteria_slug(bacteria_name);
    let specific_key = format!("{}_{}", canonical.as_ref(), param_suffix);
    PARAMETERS.get(&specific_key).copied()
}

/// Retrieves multiple bacteria-specific parameters in one pass and applies a callback.
///
/// The callback is invoked with `(bacteria_index, bacteria_name, params)` where `params`
/// is a slice in the same order as `param_suffixes`.
///
/// Returns `true` if the parameter set was available and the callback executed, `false` otherwise.
#[allow(dead_code)]
pub fn with_bacteria_params<F>(
    bacteria_name: &str,
    param_suffixes: &[&str],
    mut callback: F,
) -> bool
where
    F: FnMut(usize, &str, &[f64]),
{
    let canonical = canonicalize_bacteria_slug(bacteria_name);
    let bacteria_name = canonical.as_ref();
    let specific_keys: Vec<String> = param_suffixes
        .iter()
        .map(|suffix| format!("{}_{}", bacteria_name, suffix))
        .collect();

    let mut values = Vec::with_capacity(param_suffixes.len());
    for key in &specific_keys {
        if let Some(value) = PARAMETERS.get(key) {
            values.push(*value);
        } else {
            return false;
        }
    }

    if let Some(index) = BACTERIA_INDEX.get(bacteria_name) {
        callback(*index, bacteria_name, &values);
        true
    } else {
        false
    }
}
/// Retrieves a drug-specific simulation parameter.
/// It looks up "drug_{drug_name}_{param_suffix}".
/// Returns `Some(value)` if found, `None` otherwise.
#[allow(dead_code)]
pub fn get_drug_param(drug_name: &str, param_suffix: &str) -> Option<f64> {
    let specific_key = format!("drug_{}_{}", drug_name, param_suffix);
    PARAMETERS.get(&specific_key).copied()
}

/// Checks if a drug is available in a given region.
/// Returns the availability multiplier (0.0 = not available, 1.0 = fully available).
/// For Home region, uses the individual's region_living.
/// Special handling for colistin's historical discontinuation period.
pub fn get_drug_availability(drug_name: &str, region: &str, region_living: Option<&str>) -> f64 {
    // Handle Home region by using region_living
    let effective_region = if region == "home" {
        region_living.unwrap_or("north_america") // Default fallback if region_living not provided
    } else {
        region
    };

    let availability_key = format!("{}_drug_{}_availability", effective_region, drug_name);
    PARAMETERS.get(&availability_key).copied().unwrap_or(1.0) // Default to available if not specified
}

/// Time-aware drug availability that accounts for historical discontinuation periods.
/// Currently handles colistin's abandonment (1970-1995) and reintroduction.
pub fn get_drug_availability_time_aware(
    drug_name: &str,
    region: &str,
    region_living: Option<&str>,
    time_step: usize,
) -> f64 {
    // Calculate simulation year (assuming time_step 0 = year 1930, one step per day)
    let simulation_year = 1930.0 + (time_step as f64 / 365.0);

    // Special case for colistin's historical discontinuation
    if drug_name == "colistin" {
        // Colistin timeline:
        // 1952-1970: Active use (18 years)
        // 1970-1995: Largely abandoned due to toxicity (25 years)
        // 1995+: Reintroduced as last resort for MDR infections

        if simulation_year < 1952.0 {
            return 0.0; // Not yet introduced
        } else if simulation_year >= 1952.0 && simulation_year < 1970.0 {
            // Active use period - full availability
            return get_drug_availability(drug_name, region, region_living);
        } else if simulation_year >= 1970.0 && simulation_year < 1995.0 {
            // Discontinuation period - very limited availability (research/compassionate use only)
            return get_drug_availability(drug_name, region, region_living) * 0.05;
        // 5% of normal availability
        } else {
            // Reintroduction period - available but as last resort
            return get_drug_availability(drug_name, region, region_living);
        }
    }

    // For all other drugs, use standard availability
    get_drug_availability(drug_name, region, region_living)
}

// --- Age Risk Templates Configuration ---

lazy_static! {
    pub static ref AGE_RISK_TEMPLATES: HashMap<&'static str, Vec<f64>> = {
        let mut m = HashMap::new();

        // Age groups: 0-1, 1-5, 5-18, 18-50, 50-70, 70+
        // Values represent risk multipliers relative to baseline (18-50 age group = 1.0)

        m.insert("respiratory", vec![3.0, 1.8, 0.8, 1.0, 1.3, 2.5]);          // High infant/elderly risk (pneumonia, URI)
        m.insert("gastrointestinal", vec![2.5, 2.0, 1.2, 1.0, 1.1, 1.8]);    // High young child risk (diarrheal diseases)
        m.insert("urogenital", vec![1.2, 0.8, 0.9, 1.0, 1.4, 2.2]);          // Moderate elderly risk (UTIs)
        m.insert("skin_soft_tissue", vec![1.5, 1.3, 1.1, 1.0, 1.2, 1.8]);    // Mild age gradient
        m.insert("bloodstream", vec![4.0, 2.0, 0.7, 1.0, 1.5, 3.0]);         // Very high infant/elderly risk (sepsis)
        m.insert("sexually_transmitted", vec![0.1, 0.2, 0.8, 1.0, 0.8, 0.3]); // Peak in young adults
        m.insert("flat", vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);               // No age effect

        m
    };
}

// ---------------- 13) HGT matrices & resistance mechanism settings ----------------
// --- CROSS-RESISTANCE CONFIGURATION ---
// NOTE: These groups are DIFFERENT from the potency drug classes above!
// Potency classes = therapeutic effectiveness groupings
// Cross-resistance groups = resistance mechanism groupings (bacteria-specific)

lazy_static! {
    static ref CROSS_RESISTANCE_GROUPS: HashMap<&'static str, Vec<Vec<&'static str>>> = {
        let mut m = HashMap::new();

        // E. coli resistance patterns
        m.insert("escherichia_coli", vec![
            // ESBL resistance affects penicillins + some cephalosporins (BL/BLI combinations overcome ESBL)
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Fluoroquinolone resistance (often ciprofloxacin + levofloxacin together)
            vec!["ciprofloxacin", "levofloxacin"],
            // Aminoglycoside resistance (often linked)
            vec!["gentamicin", "tobramycin"],
        ]);

        // acinetobacter_baumannii resistance patterns
        m.insert("acinetobacter_baumannii", vec![
            // β-lactamase affects most β-lactams (BL/BLI combinations included)
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Carbapenemase affects carbapenems (including BL/BLI)
            vec!["meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin"],
            // Aminoglycoside resistance
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // klebsiella_pneumoniae resistance patterns
        m.insert("klebsiella_pneumoniae", vec![
            // ESBL resistance (BL/BLI combinations overcome ESBL)
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "ceftriaxone", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Carbapenemase (KPC, NDM, etc.)
            vec!["meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin"],
        ]);

        // streptococcus_pneumoniae resistance patterns
        m.insert("streptococcus_pneumoniae", vec![
            // Macrolide resistance (erm genes affect all macrolides)
            vec!["erythromycin", "azithromycin", "clarithromycin"],
            // Penicillin resistance (affects β-lactams)
            vec!["penicilling", "ampicillin", "amoxicillin"],
        ]);

        // staphylococcus_aureus resistance patterns
        m.insert("staphylococcus_aureus", vec![
            // β-lactamase affects penicillins
            vec!["penicilling", "ampicillin", "amoxicillin"],
            // MRSA affects most β-lactams
            vec!["cephalexin", "cefazolin", "cefuroxime", "ceftriaxone"],
            // Macrolide-lincosamide resistance
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
        ]);

        m.insert("staphylococcus_epidermidis", vec![
            // mecA-mediated resistance impacts nearly all β-lactams
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "ceftriaxone"],
            // Macrolide/clindamycin cross-resistance common in CoNS
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
            // Multidrug efflux impacting fluoroquinolones when acquired
            vec!["ciprofloxacin", "levofloxacin"],
        ]);

        m.insert("stenotrophomonas_maltophilia", vec![
            // Sulfonamide resistance often linked across TMP-SMX components
            vec!["trim_sulf"],
            // Efflux-mediated tetracycline resistance (minocycline/doxycycline)
            vec!["tetracycline", "doxycycline", "minocycline"],
            // Fluoroquinolone resistance frequently cross-links
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // pseudomonas_aeruginosa resistance patterns
        m.insert("pseudomonas_aeruginosa", vec![
            // β-lactamase affects multiple β-lactams
            vec!["piperacillin", "ceftazidime", "cefepime"],
            // Carbapenemase
            vec!["meropenem", "imipenem_c"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin"],
            // Aminoglycoside resistance
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // Enterobacter species resistance patterns
        m.insert("enterobacter_spp.", vec![
            // AmpC β-lactamase (chromosomal)
            vec!["ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime"],
            // ESBL if acquired
            vec!["ceftriaxone", "ceftazidime", "cefepime"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin"],
        ]);

        // MDR Mycobacterium tuberculosis resistance patterns
        // MDR-TB is defined by resistance to at least rifampicin + isoniazid, with guaranteed rifampicin resistance
        m.insert("mdr_mycobacterium_tuberculosis", vec![
            // Fluoroquinolone resistance (gyrA/gyrB mutations commonly affect all FQs)
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            // Aminoglycoside resistance (16S rRNA mutations can affect multiple AGs)
            vec!["gentamicin", "tobramycin", "amikacin"],
            // Rifampicin resistance (rpoB mutations - single drug class)
            vec!["rifampicin"],
            // Linezolid resistance (rrl gene mutations - single drug currently)
            vec!["linezolid"],
        ]);

        // Add more bacteria as needed...
        // The key insight: resistance groupings are bacteria-specific and mechanism-based,
        // while potency groupings are based on therapeutic similarity.

        m
    };
}

/// Returns the cross-resistance drug groups for each bacterium.
pub fn get_cross_resistance_groups() -> &'static HashMap<&'static str, Vec<Vec<&'static str>>> {
    &CROSS_RESISTANCE_GROUPS
}

/// Retrieves a string parameter (like template names).
/// Returns `Some(value)` if found, `None` otherwise.
#[allow(dead_code)]
pub fn get_string_param(key: &str) -> Option<String> {
    STRING_PARAMETERS.get(key).cloned()
}

/// Calculates the age-based infection risk multiplier for a given bacteria and age.
/// Uses the template system with bacteria-specific scaling.
/// Returns a multiplier (1.0 = baseline risk, >1.0 = increased risk, <1.0 = decreased risk)
#[allow(dead_code)]
pub fn get_age_infection_multiplier(bacteria_name: &str, age_days: i32) -> f64 {
    let canonical = canonicalize_bacteria_slug(bacteria_name);
    let bacteria_name = canonical.as_ref();
    let age_group_idx = AgeCategory::from_age_days(age_days).order();

    // Get the template name for this bacteria
    let template_key = format!("{}_age_risk_template", bacteria_name);
    let template_name =
        get_string_param(&template_key).unwrap_or_else(|| "respiratory".to_string());

    // Get the scaling factor for this bacteria
    let scaling = get_bacteria_param(bacteria_name, "age_effect_scaling").unwrap_or(1.0);

    // Look up the base multiplier from the template
    if let Some(template) = AGE_RISK_TEMPLATES.get(template_name.as_str()) {
        let base_multiplier = template[age_group_idx];
        // Scale the deviation from 1.0 by the scaling factor
        // scaling = 0.0 means no age effect (flat = 1.0)
        // scaling = 1.0 means full template effect
        // scaling > 1.0 means amplified age effect
        1.0 + (base_multiplier - 1.0) * scaling
    } else {
        // Fallback if template not found
        1.0
    }
}

/// Gets age-dependent sepsis risk log-odds for a specific bacteria and age.
/// Accounts for clinically important age-bacteria interactions (e.g., GBS in neonates, pneumococcus in elderly).
pub fn get_age_dependent_bacteria_sepsis_risk_log_odds(
    bacteria_name: &str,
    age_days: u32,
) -> f64 {
    let canonical = canonicalize_bacteria_slug(bacteria_name);
    let bacteria_name = canonical.as_ref();
    // Define age categories
    let age_category = if age_days <= 28 {
        "neonatal"
    } else if age_days <= 365 * 18 {
        "pediatric"
    } else if age_days <= 365 * 65 {
        "young_adult"
    } else {
        "elderly"
    };

    let base = get_global_param("sepsis_age_log_odds_baseline").unwrap_or(0.0);
    let age_key = format!("sepsis_age_log_odds_{}", age_category);
    let age_delta = get_global_param(&age_key).unwrap_or(0.0);
    let bacteria_age_log_key = format!("{}_{}_sepsis_log_odds", bacteria_name, age_category);
    let bacteria_age_delta = get_global_param(&bacteria_age_log_key).unwrap_or(0.0);

    base + age_delta + bacteria_age_delta
}

// Drug introduction dates (as time steps from start of 1930)
// Each time step = 1 day, so multiply years by 365
lazy_static! {
    pub static ref DRUG_INTRODUCTION_DATES: HashMap<&'static str, usize> = {
        let mut map = HashMap::new();

        // Sulfonamides (first antibiotics)
        map.insert("sulfanilamide", 2555);   // 2555 // 1937 (simulation start, 7 years after 1930)

        // Beta-lactams (Penicillins)
        map.insert("penicilling", 3555);     // 3555 // 1942 (12 years after 1930)
        map.insert("ampicillin", 11315);     // 11315 // 1961 (31 years after 1930)
        map.insert("amoxicillin", 13780);    // 1972 (42 years after 1930)
        map.insert("piperacillin", 16065);   // 1981 (51 years after 1930)
        map.insert("ticarcillin", 14600);    // 1977 (47 years after 1930)

        // Beta-lactam/beta-lactamase inhibitor combinations
        map.insert("amoxicillin_clavulanate", 16425); // 1985 (55 years after 1930)
        map.insert("ampicillin_sulbactam", 18250);    // 1990 (60 years after 1930)
        map.insert("piperacillin_tazobactam", 19715); // 1984 (54 years after 1930)
        map.insert("ticarcillin_clavulanate", 18250); // 1990 (60 years after 1930)
        map.insert("meropenem_vaborbactam", 32045);   // 2018 (88 years after 1930)
        map.insert("ceftazidime_avibactam", 27740);   // 2006 (76 years after 1930)

        // Cephalosporins
        map.insert("cephalexin", 14605);     // 1970 (40 years after 1930)
        map.insert("cefazolin", 15700);      // 1973 (43 years after 1930)
        map.insert("cefuroxime", 17525);     // 1978 (48 years after 1930)
        map.insert("ceftriaxone", 19715);    // 1984 (54 years after 1930)
        map.insert("ceftazidime", 20080);    // 1985 (55 years after 1930)
        map.insert("cefepime", 24195);       // 1996 (66 years after 1930)
        map.insert("ceftaroline", 29305);    // 2010 (80 years after 1930)

        // Carbapenems
        map.insert("meropenem", 24195);      // 1996 (66 years after 1930)
        map.insert("imipenem_c", 20080);     // 1985 (55 years after 1930)
        map.insert("ertapenem", 25920);      // 2001 (71 years after 1930)

        // Monobactams
        map.insert("aztreonam", 20445);      // 1986 (56 years after 1930)

        // Macrolides
        map.insert("erythromycin", 8025);    // 1952 (22 years after 1930)
        map.insert("azithromycin", 22260);   // 1991 (61 years after 1930)
        map.insert("clarithromycin", 21895); // 1990 (60 years after 1930)

        // Lincosamides
        map.insert("clindamycin", 13870);    // 1968 (38 years after 1930)

        // Aminoglycosides
        map.insert("gentamicin", 12045);      // 1963 (33 years after 1930)
        map.insert("tobramycin", 16325);     // 1975 (45 years after 1930)
        map.insert("amikacin", 16690);       // 1976 (46 years after 1930)

        // Fluoroquinolones
        map.insert("ciprofloxacin", 20805);  // 1987 (57 years after 1930)
        map.insert("levofloxacin", 24195);   // 1996 (66 years after 1930)
        map.insert("moxifloxacin", 25290);   // 1999 (69 years after 1930)
        map.insert("ofloxacin", 21895);      // 1990 (60 years after 1930)

        // Tetracyclines
        map.insert("tetracycline", 6575);    // 1948 (18 years after 1930)
        map.insert("doxycycline", 13505);   // 1967 (37 years after 1930)
        map.insert("minocycline", 14965);    // 1971 (41 years after 1930)

        // Glycopeptides
        map.insert("vancomycin", 10215);      // 1958 (28 years after 1930)
        map.insert("teicoplanin", 21170);    // 1988 (58 years after 1930)

        // Oxazolidinones
        map.insert("linezolid", 25550);      // 2000 (70 years after 1930)
        map.insert("tedizolid", 30660);      // 2014 (84 years after 1930)

        // Folate antagonists
        map.insert("trim_sulf", 13870);      // 1968 (38 years after 1930) - trimethoprim-sulfamethoxazole

        // Other antibiotics
        map.insert("rifampicin", 13140);     // 1966 (36 years after 1930) - critical for TB treatment
        map.insert("quinu_dalfo", 25290);    // 1999 (69 years after 1930) - quinupristin/dalfopristin
        map.insert("chlorampheni", 6935);    // 1949 (19 years after 1930) - chloramphenicol
        map.insert("nitrofurantoin", 8395);  // 1953 (23 years after 1930)
        map.insert("retapamulin", 28405);    // 2007 (77 years after 1930) - topical antibiotic
        map.insert("fusidic_a", 11680);       // 1962 (32 years after 1930) - fusidic acid
        map.insert("metronidazole", 10965);   // 1960 (30 years after 1930)
        map.insert("furazolidone", 9125);    // 1955 (25 years after 1930)

        // Polymyxins
        map.insert("colistin", 8020);        // 1952 (22 years after 1930) - clinical introduction

        map.insert("dalbavancin", 30660);   // 2014 (84 years after 1930)

        map
    };
}

/// Gets the introduction date for a drug (as time step from 1930)
/// Returns the time step if found, None otherwise
pub fn get_drug_introduction_time_step(drug_name: &str) -> Option<usize> {
    DRUG_INTRODUCTION_DATES.get(drug_name).copied()
}

/// Samples an age and region from the 120-parameter demographic distribution
/// Returns (region, age_in_years)
pub fn sample_age_and_region_from_distribution(
    rng: &mut impl rand::Rng,
) -> (crate::simulation::population::Region, i32) {
    use crate::simulation::population::Region;

    // Build cumulative probability distribution
    let mut cumulative_probs = Vec::new();
    let mut running_total = 0.0;

    // Define regions and age bands (-40000 to +32000 in 4000-year bands)
    let regions = [
        Region::Asia,
        Region::Africa,
        Region::Europe,
        Region::NorthAmerica,
        Region::SouthAmerica,
        Region::Oceania,
    ];

    let age_bands = [
        (-40000, -36000),
        (-36000, -32000),
        (-32000, -28000),
        (-28000, -24000),
        (-24000, -20000),
        (-20000, -16000),
        (-16000, -12000),
        (-12000, -8000),
        (-8000, -4000),
        (-4000, 0),
        (0, 4000),
        (4000, 8000),
        (8000, 12000),
        (12000, 16000),
        (16000, 20000),
        (20000, 24000),
        (24000, 28000),
        (28000, 32000),
    ];

    // Build distribution
    for region in &regions {
        let region_name = match region {
            Region::Asia => "asia",
            Region::Africa => "africa",
            Region::Europe => "europe",
            Region::NorthAmerica => "north_america",
            Region::SouthAmerica => "south_america",
            Region::Oceania => "oceania",
            Region::Home => "north_america", // Default fallback
        };

        for (age_min, age_max) in &age_bands {
            let param_name = if *age_min < 0 && *age_max <= 0 {
                format!(
                    "demo_{}_age_neg{}_neg{}",
                    region_name,
                    (*age_min as i32).abs(),
                    (*age_max as i32).abs()
                )
            } else if *age_min < 0 && *age_max > 0 {
                format!("demo_{}_age_neg{}_0", region_name, (*age_min as i32).abs())
            } else {
                format!("demo_{}_age_{}_{}", region_name, age_min, age_max)
            };
            let prob = get_global_param(&param_name).unwrap_or(0.0);
            running_total += prob;
            cumulative_probs.push((running_total, *region, *age_min, *age_max));
        }
    }

    // Sample from distribution
    let random_value = rng.gen::<f64>() * running_total;

    for (cumulative_prob, region, age_min, age_max) in cumulative_probs {
        if random_value <= cumulative_prob {
            // Sample a random age within the band
            let age = rng.gen_range(age_min..=age_max);
            return (region, age);
        }
    }

    // Fallback (should rarely be reached)
    (Region::Asia, 0)
}

/// Helper function to get drug interaction multiplier between two drugs
/// Returns the multiplier for drug1's level when co-administered with drug2
/// Returns 1.0 (no interaction) if no specific interaction is defined
#[allow(dead_code)]
pub fn get_drug_interaction_multiplier(drug1: &str, drug2: &str) -> f64 {
    let interaction_key = format!(
        "drug_level_multiplier_{}_when_coadministered_with_{}",
        drug1, drug2
    );
    get_global_param(&interaction_key).unwrap_or(1.0)
}

/// Helper function to check if two drugs have a defined interaction
#[allow(dead_code)]
pub fn drugs_have_interaction(drug1: &str, drug2: &str) -> bool {
    let interaction_key = format!(
        "drug_level_multiplier_{}_when_coadministered_with_{}",
        drug1, drug2
    );
    get_global_param(&interaction_key).is_some()
}

/// Helper function to list all active drug interactions for debugging/analysis
#[allow(dead_code)]
pub fn get_all_active_interactions() -> Vec<(String, String, f64)> {
    let mut interactions = Vec::new();

    // This would require iterating through PARAMETERS to find interaction keys
    // For now, return the known interactions from our configuration
    let known_interactions = vec![
        ("levofloxacin", "rifampicin", 0.7),
        ("moxifloxacin", "rifampicin", 0.8),
        ("clarithromycin", "rifampicin", 0.6),
        ("azithromycin", "rifampicin", 0.8),
        ("ciprofloxacin", "amoxicillin_clavulanate", 0.85),
        ("levofloxacin", "amoxicillin_clavulanate", 0.9),
        ("ciprofloxacin", "erythromycin", 0.85),
        ("levofloxacin", "azithromycin", 0.9),
    ];

    for (drug1, drug2, multiplier) in known_interactions {
        interactions.push((drug1.to_string(), drug2.to_string(), multiplier));
    }

    interactions
}
