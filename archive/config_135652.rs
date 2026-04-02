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
    AgeCategory, BacteriaGroup, Region, ResistanceMechanism, AGE_CATEGORY_SEQUENCE,
    BACTERIA_GROUPS, BACTERIA_LIST, DRUG_SHORT_NAMES,
};
use lazy_static::lazy_static;
use rand::Rng;
use std::borrow::Cow;
use std::collections::HashMap; // Import both lists and helper enums
use std::ptr;
use std::sync::atomic::{AtomicPtr, Ordering};

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
    pub bacteria_mechanism_emergence: BacteriaMechanismEmergenceRates,
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
        let bacteria_mechanism_emergence =
            BacteriaMechanismEmergenceRates::from_map(map, BACTERIA_LIST.len());

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
            bacteria_mechanism_emergence,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RunParameterSamplingConfig {
    pub enabled: bool,
}

impl RunParameterSamplingConfig {
    pub const fn disabled() -> Self {
        Self { enabled: false }
    }

    pub const fn enabled() -> Self {
        Self { enabled: true }
    }
}

#[derive(Debug, Clone)]
pub struct SampledParameterRecord {
    pub record_kind: &'static str,
    pub sampled_quantity: &'static str,
    pub applies_to: &'static str,
    pub baseline_value: f64,
    pub sampled_value: f64,
    pub transform: &'static str,
    pub draw_value: f64,
}

#[derive(Debug)]
struct ActiveParameterContext {
    map: HashMap<String, f64>,
    store: ParameterStore,
}

pub const RUN_PATHWAY_INFECTION_DE_NOVO_MULTIPLIER_KEY: &str =
    "run_pathway_infection_de_novo_multiplier";
pub const RUN_PATHWAY_MICROBIOME_DE_NOVO_MULTIPLIER_KEY: &str =
    "run_pathway_microbiome_de_novo_multiplier";
pub const RUN_PATHWAY_HGT_MULTIPLIER_KEY: &str = "run_pathway_hgt_multiplier";
pub const RUN_PATHWAY_REVERSION_RATE_MULTIPLIER_KEY: &str =
    "run_pathway_reversion_rate_multiplier";
pub const RUN_PATHWAY_CARRIER_INHERITANCE_MULTIPLIER_KEY: &str =
    "run_pathway_carrier_inheritance_multiplier";
pub const RUN_PATHWAY_COMMUNITY_DILUTION_MULTIPLIER_KEY: &str =
    "run_pathway_community_dilution_multiplier";
pub const RUN_PATHWAY_MICROBIOME_ACQUISITION_MULTIPLIER_KEY: &str =
    "run_pathway_microbiome_acquisition_multiplier";
pub const RUN_PATHWAY_MICROBIOME_DISRUPTION_MULTIPLIER_KEY: &str =
    "run_pathway_microbiome_disruption_multiplier";
pub const RUN_PATHWAY_DRUG_ACTIVITY_MULTIPLIER_KEY: &str =
    "run_pathway_drug_activity_multiplier";

struct RunLevelSamplingSpec {
    parameter_key: &'static str,
    applies_to: &'static str,
    min_multiplier: f64,
    max_multiplier: f64,
}

const RUN_LEVEL_SAMPLING_SPECS: [RunLevelSamplingSpec; 9] = [
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_INFECTION_DE_NOVO_MULTIPLIER_KEY,
        applies_to: "infection_de_novo_pathway",
        min_multiplier: 0.01,
        max_multiplier: 100.0,
    },
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_MICROBIOME_DE_NOVO_MULTIPLIER_KEY,
        applies_to: "microbiome_de_novo_pathway",
        min_multiplier: 0.01,
        max_multiplier: 100.0,
    },
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_HGT_MULTIPLIER_KEY,
        applies_to: "hgt_pathway",
        min_multiplier: 0.01,
        max_multiplier: 100.0,
    },
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_REVERSION_RATE_MULTIPLIER_KEY,
        applies_to: "mechanism_reversion_rate",
        min_multiplier: 0.1,
        max_multiplier: 10.0,
    },
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_CARRIER_INHERITANCE_MULTIPLIER_KEY,
        applies_to: "carrier_resistance_inheritance",
        min_multiplier: 0.2,
        max_multiplier: 2.0,
    },
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_COMMUNITY_DILUTION_MULTIPLIER_KEY,
        applies_to: "community_resistance_dilution",
        min_multiplier: 0.2,
        max_multiplier: 2.0,
    },
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_MICROBIOME_ACQUISITION_MULTIPLIER_KEY,
        applies_to: "microbiome_acquisition_seeding",
        min_multiplier: 0.1,
        max_multiplier: 2.0,
    },
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_MICROBIOME_DISRUPTION_MULTIPLIER_KEY,
        applies_to: "microbiome_disruption_log_odds",
        min_multiplier: 0.2,
        max_multiplier: 5.0,
    },
    RunLevelSamplingSpec {
        parameter_key: RUN_PATHWAY_DRUG_ACTIVITY_MULTIPLIER_KEY,
        applies_to: "drug_activity_efficacy",
        min_multiplier: 0.5,
        max_multiplier: 2.0,
    },
];

#[derive(Clone, Copy)]
enum SamplingStrategy {
    MultiplicativeLogUniform {
        min_multiplier: f64,
        max_multiplier: f64,
    },
}

fn sample_value<R: Rng + ?Sized>(
    baseline_value: f64,
    strategy: SamplingStrategy,
    rng: &mut R,
) -> Option<(f64, &'static str, f64)> {
    if baseline_value.abs() <= f64::EPSILON {
        return None;
    }

    match strategy {
        SamplingStrategy::MultiplicativeLogUniform {
            min_multiplier,
            max_multiplier,
        } => {
            let min_log = min_multiplier.ln();
            let max_log = max_multiplier.ln();
            let multiplier = rng.gen_range(min_log..=max_log).exp();
            Some((
                baseline_value * multiplier,
                "log_uniform_multiplier",
                multiplier,
            ))
        }
    }
}

fn sample_parameter_records<R: Rng + ?Sized>(
    parameter_map: &mut HashMap<String, f64>,
    rng: &mut R,
) -> Vec<SampledParameterRecord> {
    let mut sampled_records = Vec::with_capacity(RUN_LEVEL_SAMPLING_SPECS.len());

    for spec in RUN_LEVEL_SAMPLING_SPECS.iter() {
        let baseline_value = 1.0;
        let (sampled_value, transform, draw_value) = sample_value(
            baseline_value,
            SamplingStrategy::MultiplicativeLogUniform {
                min_multiplier: spec.min_multiplier,
                max_multiplier: spec.max_multiplier,
            },
            rng,
        )
        .expect("run-level pathway sampling unexpectedly skipped");

        parameter_map.insert(spec.parameter_key.to_string(), sampled_value);
        sampled_records.push(SampledParameterRecord {
            record_kind: "latent_pathway_multiplier",
            sampled_quantity: spec.parameter_key,
            applies_to: spec.applies_to,
            baseline_value,
            sampled_value,
            transform,
            draw_value,
        });
    }

    sampled_records
}

fn get_active_parameter_context() -> Option<&'static ActiveParameterContext> {
    let context_ptr = ACTIVE_PARAMETER_CONTEXT.load(Ordering::Acquire);
    if context_ptr.is_null() {
        None
    } else {
        // Safe because the context is leaked for the lifetime of the process.
        unsafe { context_ptr.as_ref() }
    }
}

fn parameter_map() -> &'static HashMap<String, f64> {
    if let Some(context) = get_active_parameter_context() {
        &context.map
    } else {
        &PARAMETERS
    }
}

fn set_active_parameter_context(parameter_map: HashMap<String, f64>) {
    let store = ParameterStore::from_parameter_map(&parameter_map);
    // Leak the run context so existing read-heavy code can keep borrowing static references.
    let leaked_context = Box::leak(Box::new(ActiveParameterContext {
        map: parameter_map,
        store,
    }));
    ACTIVE_PARAMETER_CONTEXT.store(leaked_context as *mut ActiveParameterContext, Ordering::Release);
}

pub fn clear_active_run_parameters() {
    ACTIVE_PARAMETER_CONTEXT.store(ptr::null_mut(), Ordering::Release);
}

pub fn activate_run_parameter_sampling<R: Rng + ?Sized>(
    rng: &mut R,
    sampling_config: RunParameterSamplingConfig,
) -> Vec<SampledParameterRecord> {
    if !sampling_config.enabled {
        clear_active_run_parameters();
        return Vec::new();
    }

    let mut sampled_parameter_map = PARAMETERS.clone();
    let sampled_records = sample_parameter_records(&mut sampled_parameter_map, rng);
    set_active_parameter_context(sampled_parameter_map);
    sampled_records
}

// ---------------- 3) Global scalar defaults & helpers ----------------
#[derive(Debug)]
pub struct GlobalScalars {
    // Logistic model for antibiotic initiation probability
    // P(initiation) = 1 / (1 + exp(-log_odds)) where log_odds = base + sum of applicable effects
    pub antibiotic_initiation_base_log_odds: f64,
    pub antibiotic_initiation_log_odds_symptomatic_infection: f64,
    pub antibiotic_initiation_log_odds_test_identified: f64,
    pub antibiotic_initiation_log_odds_already_on_drug: f64,
    pub antibiotic_initiation_log_odds_immunodeficiency: f64,
    pub antibiotic_initiation_log_odds_sepsis: f64, // New parameter for sepsis cases
    pub antibiotic_initiation_log_odds_no_indication: f64, // negative effect when prescribing without infection/immunodeficiency
    // Drug activity parameters (still used for bacteria level effects)
    pub drug_activity_to_bacteria_level_multiplier: f64,
    pub drug_activity_slow_clearance_probability: f64,
    pub drug_activity_slow_clearance_multiplier: f64,
    pub double_dose_probability_if_identified_infection: f64,
    pub random_drug_cessation_probability: f64,
    pub random_drug_cessation_probability_if_no_active_infection: f64,
    pub microbiome_resistance_transfer_probability_per_day: f64,
    // Logistic model for hospitalization probability
    // P(hospitalization) = 1 / (1 + exp(-log_odds)) where log_odds = base + age_effect + sepsis_effect + symptomatic_infection_effect
    pub hospitalization_base_log_odds: f64,
    pub hospitalization_log_odds_per_age_year: f64,
    pub hospitalization_log_odds_sepsis: f64,
    pub hospitalization_log_odds_symptomatic_infection: f64,
    pub hospitalization_symptomatic_infection_level_threshold: f64,
    pub hospital_recovery_rate_per_day: f64,
    pub hospital_max_days: f64,
    pub hospital_prevent_discharge_with_sepsis: f64,
    pub travel_probability_per_day: f64,
    pub antibiotic_infection_prevention_efficacy: f64,
    pub max_resistance_level: f64,
    pub resistance_emergence_bacteria_level_multiplier: f64,
    #[allow(dead_code)]
    pub any_r_emergence_level_on_first_emergence: f64,
    pub multi_drug_penalty_threshold_num_drugs: f64,
    pub resistance_development_inhibition_single_drug: f64,
    pub resistance_development_inhibition_partial_cross: f64,
    pub mechanism_assignment_probability_on_any_r_gain: f64,
    pub mechanism_cache_ewma_decay: f64,
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
    // Per-region sepsis onset log-odds for fine-grained regional differentiation
    pub log_odds_sepsis_onset_region_north_america: f64,
    pub log_odds_sepsis_onset_region_europe: f64,
    pub log_odds_sepsis_onset_region_oceania: f64,
    pub log_odds_sepsis_onset_region_asia: f64,
    pub log_odds_sepsis_onset_region_south_america: f64,
    pub log_odds_sepsis_onset_region_africa: f64,
    pub mdr_tb_pre_antibiotic_era_multiplier: f64,
    pub mdr_tb_early_antibiotic_era_multiplier: f64,
    pub mdr_tb_modern_era_multiplier: f64,
    pub default_toxicity_reservoir_half_life_days: f64,
    // Drug toxicity death multiplicative model parameters
    pub toxicity_age_multiplier_infant: f64,
    pub toxicity_age_multiplier_child: f64,
    pub toxicity_age_multiplier_adult: f64,
    pub toxicity_age_multiplier_elderly: f64,
    pub toxicity_immunosuppressed_multiplier: f64,
    pub toxicity_hospital_multiplier: f64,
    pub toxicity_discontinuation_threshold: f64,
    pub toxicity_discontinuation_avoidance_days: i32,
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
    // Sepsis onset additional factors
    pub log_odds_sepsis_onset_immunosuppressed: f64,
    pub log_odds_sepsis_onset_hospitalized: f64,
    pub log_odds_sepsis_onset_not_under_care: f64,
    // Sepsis death logistic model parameters (log-odds scale)
    pub sepsis_death_base_log_odds: f64,
    pub sepsis_death_log_odds_age_infant: f64,
    pub sepsis_death_log_odds_age_child: f64,
    pub sepsis_death_log_odds_age_adult: f64,
    pub sepsis_death_log_odds_age_elderly: f64,
    pub sepsis_death_log_odds_immunosuppressed: f64,
    pub sepsis_death_log_odds_bacteria_level: f64,
    pub sepsis_death_log_odds_duration: f64,
    pub sepsis_death_log_odds_early_phase: f64,
    pub sepsis_death_early_phase_days: f64,
    pub sepsis_death_log_odds_not_under_care: f64,
    // Bacteria growth rate multipliers (affect base_bacteria_level_change)
    pub bacteria_growth_age_multiplier_infant: f64,
    pub bacteria_growth_age_multiplier_child: f64,
    pub bacteria_growth_age_multiplier_adult: f64,
    pub bacteria_growth_age_multiplier_elderly: f64,
    pub bacteria_growth_immunodeficiency_multiplier: f64,
    // Enhanced microbiome/carriage model parameters
    pub antibiotic_disruption_decay_half_life_days: f64,
    pub microbiome_resistance_multiplier_on_acquisition: f64,
    pub infection_from_microbiome_dampening: f64,
    pub carriage_duration_log_odds_coefficient: f64,
    pub carriage_duration_max_log_odds_effect: f64,
    pub antibiotic_clearance_log_odds_per_unit_activity: f64,
    pub carrier_resistance_inheritance_probability: f64,
    pub community_resistance_dilution_factor: f64,
    pub hgt_hospital_multiplier: f64,
    pub hgt_antibiotic_pressure_multiplier: f64,
    pub hgt_coinfection_multiplier: f64,
    pub hgt_microbiome_only_penalty: f64,
    pub hgt_gut_compartment_multiplier: f64,
    pub hgt_minority_donor_multiplier: f64,
    #[allow(dead_code)]
    pub majority_r_memory_retention_per_day: f64,
    pub mechanism_reversion_rate_global_multiplier: f64,
    pub infection_de_novo_multiplier: f64,
    pub microbiome_de_novo_multiplier: f64,
    pub hgt_multiplier: f64,
    #[allow(dead_code)]
    pub microbiome_majority_decay_half_life_days: f64,
    #[allow(dead_code)]
    pub microbiome_minority_decay_half_life_days: f64,
    #[allow(dead_code)]
    pub microbiome_majority_promotion_rate_per_day: f64,
}

impl GlobalScalars {
    fn from_map(map: &HashMap<String, f64>) -> Self {
        // Reads configuration values already present in `map`; the fallback literal only applies when no entry exists.
        GlobalScalars {
            // Logistic antibiotic initiation parameters
            // Defaults calibrated to produce similar rates to old multiplicative model
            // Base: P=0.001 → log(0.001/0.999) ≈ -6.9
            // With symptomatic infection: P=0.26 → need +5.85 log-odds
            antibiotic_initiation_base_log_odds: get_or_default(
                map,
                "antibiotic_initiation_base_log_odds",
                -6.5,
            ),
            antibiotic_initiation_log_odds_symptomatic_infection: get_or_default(
                map,
                "antibiotic_initiation_log_odds_symptomatic_infection",
                6.5,
            ),
             antibiotic_initiation_log_odds_sepsis: get_or_default(
                map,
                "antibiotic_initiation_log_odds_sepsis",
                6.0,
            ),
            antibiotic_initiation_log_odds_test_identified: get_or_default(
                map,
                "antibiotic_initiation_log_odds_test_identified",
                0.92, // log(2.5) ≈ 0.92
            ),
            antibiotic_initiation_log_odds_already_on_drug: get_or_default(
                map,
                "antibiotic_initiation_log_odds_already_on_drug",
                0.18, // log(1.2) ≈ 0.18
            ),
            antibiotic_initiation_log_odds_immunodeficiency: get_or_default(
                map,
                "antibiotic_initiation_log_odds_immunodeficiency",
                2.08, // log(8.0) ≈ 2.08
            ),
            antibiotic_initiation_log_odds_no_indication: get_or_default(
                map,
                "antibiotic_initiation_log_odds_no_indication",
                -1.05, // log(0.35) ≈ -1.05, reduces odds when no clinical indication
            ),
            // Drug activity parameters (still used for bacteria level effects)
            drug_activity_to_bacteria_level_multiplier: get_or_default(
                map,
                "drug_activity_to_bacteria_level_multiplier",
                0.75,
            ),
            drug_activity_slow_clearance_probability: get_or_default(
                map,
                "drug_activity_slow_clearance_probability",
                0.1,
            ),
            drug_activity_slow_clearance_multiplier: get_or_default(
                map,
                "drug_activity_slow_clearance_multiplier",
                0.2,
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
            microbiome_resistance_transfer_probability_per_day: get_or_default(
                map,
                "microbiome_resistance_transfer_probability_per_day",
                0.000_000_000_000_000_000_000_000_000_000_000_000_000_1,  // 0.0008  ***
            ),
            // Logistic hospitalization parameters
            // Calibrated to produce ~3-4% baseline hospitalization rate
            // Base: -5.0 gives ~0.67% daily admission, maintaining steady-state ~3-4% hospitalized
            hospitalization_base_log_odds: get_or_default(
                map,
                "hospitalization_base_log_odds",
                -5.0,
            ),
            hospitalization_log_odds_per_age_year: get_or_default(
                map,
                "hospitalization_log_odds_per_age_year",
                0.02, // ~2% increase in log-odds per year of age
            ),
            hospitalization_log_odds_sepsis: get_or_default(
                map,
                "hospitalization_log_odds_sepsis",
                4.4, // log(80) ≈ 4.4, equivalent to 80x multiplier
            ),
            hospitalization_log_odds_symptomatic_infection: get_or_default(
                map,
                "hospitalization_log_odds_symptomatic_infection",
                2.5, // ~12x multiplier: severe symptomatic infection drives hospitalization
            ),
            hospitalization_symptomatic_infection_level_threshold: get_or_default(
                map,
                "hospitalization_symptomatic_infection_level_threshold",
                3.0, // Only infections above this level contribute to hospitalization
            ),
            hospital_recovery_rate_per_day: get_or_default(
                map,
                "hospitalization_recovery_rate_per_day",
                0.25,
            ),
            hospital_max_days: get_or_default(map, "hospitalization_max_days", 30.0),
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
                0.0,
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
            mechanism_cache_ewma_decay: get_or_default(
                map,
                "mechanism_cache_ewma_decay",
                0.9,
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
            drug_selection_temperature: get_or_default(map, "drug_selection_temperature", 0.70),
            reserve_drug_score_penalty: get_or_default(map, "reserve_drug_score_penalty", 0.00005),
            restart_window_enabled: get_or_default(map, "enable_restart_window", 1.0) > 0.5,
            restart_window_days: get_or_default(map, "restart_window_days", 5.0) as i32,
            restart_bacteria_level_threshold: get_or_default(
                map,
                "restart_bacteria_level_threshold",
                1.5,
            ),
            restart_window_probability: get_or_default(map, "restart_window_probability", 0.3),
            // Per-region sepsis onset log-odds
            log_odds_sepsis_onset_region_north_america: get_or_default(map, "log_odds_sepsis_onset_region_north_america", -0.5),
            log_odds_sepsis_onset_region_europe: get_or_default(map, "log_odds_sepsis_onset_region_europe", -0.6),
            log_odds_sepsis_onset_region_oceania: get_or_default(map, "log_odds_sepsis_onset_region_oceania", -0.5),
            log_odds_sepsis_onset_region_asia: get_or_default(map, "log_odds_sepsis_onset_region_asia", -0.1),
            log_odds_sepsis_onset_region_south_america: get_or_default(map, "log_odds_sepsis_onset_region_south_america", 0.0),
            log_odds_sepsis_onset_region_africa: get_or_default(map, "log_odds_sepsis_onset_region_africa", 0.1),
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
            default_toxicity_reservoir_half_life_days: get_or_default(
                map,
                "default_toxicity_reservoir_half_life_days",
                1.5,
            )
            .max(0.0),
            // Drug toxicity death multiplicative model parameters
            toxicity_age_multiplier_infant: get_or_default(
                map,
                "toxicity_age_multiplier_infant",
                1.8, // Neonates more vulnerable to severe toxicity
            ),
            toxicity_age_multiplier_child: get_or_default(
                map,
                "toxicity_age_multiplier_child",
                1.2,
            ),
            toxicity_age_multiplier_adult: get_or_default(
                map,
                "toxicity_age_multiplier_adult",
                1.0,
            ),
            toxicity_age_multiplier_elderly: get_or_default(
                map,
                "toxicity_age_multiplier_elderly",
                2.2,
            ),
            toxicity_immunosuppressed_multiplier: get_or_default(
                map,
                "toxicity_immunosuppressed_multiplier",
                2.5,
            ),
            toxicity_hospital_multiplier: get_or_default(
                map,
                "toxicity_hospital_multiplier",
                1.3, // Hospitalized patients often sicker but also monitored
            ),
            toxicity_discontinuation_threshold: get_or_default(
                map,
                "toxicity_discontinuation_threshold",
                0.000_01, // Sub-lethal threshold: clinician stops most-toxic drug before death risk is significant
            ),
            toxicity_discontinuation_avoidance_days: get_or_default(
                map,
                "toxicity_discontinuation_avoidance_days",
                30.0, // Days to avoid re-prescribing the toxicity-stopped drug
            ) as i32,
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
            log_odds_sepsis_onset_immunosuppressed: get_or_default(
                map,
                "log_odds_sepsis_onset_immunosuppressed",
                0.7,
            ),
            log_odds_sepsis_onset_hospitalized: get_or_default(
                map,
                "log_odds_sepsis_onset_hospitalized",
                0.5,
            ),
            log_odds_sepsis_onset_not_under_care: get_or_default(
                map,
                "log_odds_sepsis_onset_not_under_care",
                1.0,
            ),
            // Sepsis death logistic model parameters
            sepsis_death_base_log_odds: get_or_default(
                map,
                "sepsis_death_base_log_odds",
                -4.2, // ~1.5% baseline daily mortality
            ),
            sepsis_death_log_odds_age_infant: get_or_default(
                map,
                "sepsis_death_log_odds_age_infant",
                1.1, // ~3x baseline for infants
            ),
            sepsis_death_log_odds_age_child: get_or_default(
                map,
                "sepsis_death_log_odds_age_child",
                -0.7, // ~0.5x baseline for children
            ),
            sepsis_death_log_odds_age_adult: get_or_default(
                map,
                "sepsis_death_log_odds_age_adult",
                0.0, // Reference category
            ),
            sepsis_death_log_odds_age_elderly: get_or_default(
                map,
                "sepsis_death_log_odds_age_elderly",
                0.9, // ~2.5x baseline for elderly
            ),
            sepsis_death_log_odds_immunosuppressed: get_or_default(
                map,
                "sepsis_death_log_odds_immunosuppressed",
                1.5, // ~4.5x baseline for immunosuppressed
            ),
            sepsis_death_log_odds_bacteria_level: get_or_default(
                map,
                "sepsis_death_log_odds_bacteria_level",
                0.35, // Per unit bacteria level (0-5 scale), ~1.4x per unit
            ),
            sepsis_death_log_odds_duration: get_or_default(
                map,
                "sepsis_death_log_odds_duration",
                0.04, // Per day of sepsis (after early phase)
            ),
            sepsis_death_log_odds_early_phase: get_or_default(
                map,
                "sepsis_death_log_odds_early_phase",
                0.8, // Additional risk in first 72 hours (~2.2x)
            ),
            sepsis_death_early_phase_days: get_or_default(
                map,
                "sepsis_death_early_phase_days",
                3.0, // Early phase lasts 3 days
            ),
            sepsis_death_log_odds_not_under_care: get_or_default(
                map,
                "sepsis_death_log_odds_not_under_care",
                1.4, // ~4x mortality without treatment
            ),
            // Bacteria growth rate multipliers by age
            bacteria_growth_age_multiplier_infant: get_or_default(
                map,
                "bacteria_growth_age_multiplier_infant",
                1.3, // Immature immune system → faster bacterial proliferation
            ),
            bacteria_growth_age_multiplier_child: get_or_default(
                map,
                "bacteria_growth_age_multiplier_child",
                1.0, // Baseline
            ),
            bacteria_growth_age_multiplier_adult: get_or_default(
                map,
                "bacteria_growth_age_multiplier_adult",
                1.0, // Baseline
            ),
            bacteria_growth_age_multiplier_elderly: get_or_default(
                map,
                "bacteria_growth_age_multiplier_elderly",
                1.2, // Immunosenescence → reduced containment
            ),
            bacteria_growth_immunodeficiency_multiplier: get_or_default(
                map,
                "bacteria_growth_immunodeficiency_multiplier",
                1.5, // Compromised immunity → faster bacterial proliferation
            ),
            antibiotic_disruption_decay_half_life_days: get_or_default(
                map,
                "antibiotic_disruption_decay_half_life_days",
                30.0,
            ),
            // ── CALIBRATION AXIS 6 (default): microbiome seeding probability ──
            microbiome_resistance_multiplier_on_acquisition: get_or_default(
                map,
                "microbiome_resistance_multiplier_on_acquisition",
                0.18,  // 0.18 ***  (overridden to 0.50 in PARAMETERS block ~line 10744)
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
            // ── CALIBRATION AXIS 4 (default): microbiome→infection bridge ──
            carrier_resistance_inheritance_probability: get_or_default(
                map,
                "carrier_resistance_inheritance_probability",
                0.32,  // 0.32 ***  (overridden to 0.50 in PARAMETERS block ~line 10790)
            ),
            // ── CALIBRATION AXIS 5 (default): fraction from human reservoir ──
            community_resistance_dilution_factor: get_or_default(
                map,
                "community_resistance_dilution_factor",
                0.50,  // 0.50 ***  (overridden to 0.50 in PARAMETERS block ~line 10791)
            ),

            // =====================================================================
            // COUNTERFACTUAL HGT CONTROL
            // To run a "No Horizontal Gene Transfer" counterfactual scenario, set
            // the `hgt_hospital_multiplier` parameter to 0.0 (or set it here below).
            // Since this parameter multiplies every HGT calculation, 0.0 effectively
            // disables all HGT events across the simulation without needing to touch
            // the HgtMatrix pair-probabilities or modify the core Rust simulation loop.
            // =====================================================================
            hgt_hospital_multiplier: get_or_default(map, "hgt_hospital_multiplier", 3.0),  // 3.0  ***
            hgt_antibiotic_pressure_multiplier: get_or_default(
                map,
                "hgt_antibiotic_pressure_multiplier",
                1.5,  // 1.5  ***
            ),
            hgt_coinfection_multiplier: get_or_default(
                map,
                "hgt_coinfection_multiplier",
                1.25,  // 1.25 ***
            ),
            hgt_microbiome_only_penalty: get_or_default(
                map,
                "hgt_microbiome_only_penalty",
                0.65,  // 0.65 ***   ^^^ micro
            ),
            hgt_gut_compartment_multiplier: get_or_default(
                map,
                "hgt_gut_compartment_multiplier",
                2.0,  // Gut has higher bacterial density and more conjugation opportunities
            ),
            hgt_minority_donor_multiplier: get_or_default(
                map,
                "hgt_minority_donor_multiplier",
                0.20,
            ),
            majority_r_memory_retention_per_day: get_or_default(
                map,
                "majority_r_memory_retention_per_day",
                0.93,
            ),
            // ^^^ micro
            // ── CALIBRATION AXIS 1 (default): resistance persistence / reversion speed ──
            mechanism_reversion_rate_global_multiplier: get_or_default(
                map,
                "mechanism_reversion_rate_global_multiplier",
                1.0,
            ),
            // ── CALIBRATION AXIS 2 (default): de-novo emergence in infections ──
            infection_de_novo_multiplier: get_or_default(
                map,
                "infection_de_novo_multiplier",
                1.0,
            ),
            // ── CALIBRATION AXIS 3 (default): de-novo emergence in gut carriage ──
            microbiome_de_novo_multiplier: get_or_default(
                map,
                "microbiome_de_novo_multiplier",
                1.0,
            ),
            // ── CALIBRATION AXIS 7 (default): horizontal gene transfer rate scaling ──
            hgt_multiplier: get_or_default(
                map,
                "hgt_multiplier",
                1.0,
            ),
            microbiome_majority_decay_half_life_days: get_or_default(
                map,
                "microbiome_majority_decay_half_life_days",
                60.0,
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
        }
    }
}

// ---------------- 4) Immunodeficiency / region / syndrome / sex parameters ----------------
#[derive(Debug)]
pub struct ImmunodeficiencyParameters {
    startup_seed_fraction: f64,
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
    antibiotic_initiation_log_odds: [f64; RegionParameters::REGION_COUNT],
    hospitalization_log_odds: [f64; RegionParameters::REGION_COUNT],
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
        let mut antibiotic_initiation_log_odds = [0.0; RegionParameters::REGION_COUNT];
        let mut hospitalization_log_odds = [0.0; RegionParameters::REGION_COUNT];

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
            antibiotic_initiation_log_odds[idx] =
                get_or_default(map, &format!("{}_antibiotic_initiation_log_odds", key_prefix), 0.0);
            hospitalization_log_odds[idx] =
                get_or_default(map, &format!("{}_hospitalization_log_odds", key_prefix), 0.0);
        }

        RegionParameters {
            travel_multiplier,
            cessation_multiplier,
            mortality_log_odds,
            sepsis_log_odds,
            sepsis_mortality_multiplier,
            testing_multiplier,
            antibiotic_initiation_log_odds,
            hospitalization_log_odds,
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
    pub fn antibiotic_initiation_log_odds(&self, region: Region) -> f64 {
        self.antibiotic_initiation_log_odds[Self::region_index(region)]
    }

    #[inline]
    pub fn hospitalization_log_odds(&self, region: Region) -> f64 {
        self.hospitalization_log_odds[Self::region_index(region)]
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
    bacteria_growth_multiplier: Vec<f64>,
    /// Drug penetration multipliers by syndrome: [syndrome_id][drug_idx] -> penetration factor (0.0-1.0)
    /// Accounts for tissue/compartment-specific drug distribution
    drug_penetration: Vec<Vec<f64>>,
}

impl SyndromeParameters {
    const MAX_SYNDROME_ID: usize = 10;

    fn from_map(map: &HashMap<String, f64>) -> Self {
        let len = Self::MAX_SYNDROME_ID + 1;
        let num_drugs = DRUG_SHORT_NAMES.len();
        let mut sepsis_log_odds = vec![0.0; len];
        let mut initiation_multiplier = vec![1.0; len];
        let mut non_sepsis_mortality_log_odds = vec![0.0; len];
        let mut empiric_drug_scores = vec![vec![0.01; num_drugs]; len];
        let mut bacteria_growth_multiplier = vec![1.0; len];
        let mut drug_penetration = vec![vec![1.0; num_drugs]; len];

        // Initialize syndrome-specific defaults for drug penetration
        // Syndromes: 1=UTI, 2=Skin, 3=Respiratory, 4=Bloodstream, 5=Intra-abdominal, 
        //           6=CNS, 7=GI, 8=Genital, 9=Bone/joint, 10=Other
        // Drug class-based penetration factors by syndrome
        Self::initialize_drug_penetration_defaults(&mut drug_penetration);

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
            bacteria_growth_multiplier[syndrome_id] = get_or_default(
                map,
                &format!("syndrome_{}_bacteria_growth_multiplier", syndrome_id),
                1.0,
            );

            for (drug_idx, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
                let key = format!("syndrome_{}_empiric_drug_{}_score", syndrome_id, drug);
                empiric_drug_scores[syndrome_id][drug_idx] = get_or_default(map, &key, 0.01);
                
                // Override drug penetration from config if specified
                let penetration_key = format!("syndrome_{}_drug_{}_penetration", syndrome_id, drug);
                drug_penetration[syndrome_id][drug_idx] = get_or_default(
                    map,
                    &penetration_key,
                    drug_penetration[syndrome_id][drug_idx],
                );
            }
        }

        SyndromeParameters {
            sepsis_log_odds,
            initiation_multiplier,
            non_sepsis_mortality_log_odds,
            empiric_drug_scores,
            bacteria_growth_multiplier,
            drug_penetration,
        }
    }
    
    /// Initialize drug penetration defaults based on pharmacokinetic properties
    /// Values represent fraction of serum concentration achieved at infection site
    fn initialize_drug_penetration_defaults(drug_penetration: &mut Vec<Vec<f64>>) {
        // Drug indices (from DRUG_SHORT_NAMES):
        // 0=sulfanilamide, 1=penicillin_g, 2=ampicillin, 3=amoxicillin, 4=piperacillin, 5=ticarcillin
        // 6=cephalexin, 7=cefazolin, 8=cefuroxime, 9=ceftriaxone, 10=ceftazidime, 11=cefepime
        // 12=ceftaroline, 13=ceftolozane_tazobactam, 14=cefiderocol, 15=meropenem, 16=imipenem_c
        // 17=ertapenem, 18=aztreonam, 19=erythromycin, 20=azithromycin, 21=clarithromycin
        // 22=clindamycin, 23=gentamicin, 24=tobramycin, 25=amikacin
        // 26=ciprofloxacin, 27=levofloxacin, 28=moxifloxacin, 29=ofloxacin
        // 30=tetracycline, 31=doxycycline, 32=minocycline, 33=tigecycline
        // 34=vancomycin, 35=teicoplanin, 36=dalbavancin, 37=linezolid, 38=tedizolid
        // 39=daptomycin, 40=quinu_dalfo, 41=trim_sulf, 42=chloramphenicol, 43=nitrofurantoin
        // 44=fosfomycin, 45=retapamulin, 46=fusidic_a, 47=metronidazole, 48=fidaxomicin
        // 49=furazolidone, 50=rifampicin, 51=amoxicillin_clavulanate, 52=piperacillin_tazobactam
        // 53=ampicillin_sulbactam, 54=ticarcillin_clavulanate, 55=ceftazidime_avibactam
        // 56=meropenem_vaborbactam, 57=colistin, 58=flucloxacillin, 59=aztreonam_avibactam
        // 60=cefixime
        
        // Syndrome indices: 1=UTI, 2=Skin, 3=Resp, 4=BSI, 5=Intra-abd, 6=CNS, 7=GI, 8=Genital, 9=Bone, 10=Other
        
        // Define drug class groupings for easier assignment
        let penicillins = [1, 2, 3, 4, 5, 58]; // penicillin_g through flucloxacillin
        let oral_cephalosporins = [6]; // cephalexin
        let iv_cephalosporins = [7, 8, 9, 10, 11, 12, 13, 14]; // cefazolin through cefiderocol
        let carbapenems = [15, 16, 17]; // meropenem, imipenem, ertapenem
        let aztreonam_idx = 18;
        let ceftriaxone_idx = 9;
        let meropenem_idx = 15;
        let macrolides = [19, 20, 21]; // erythromycin, azithromycin, clarithromycin
        let clindamycin_idx = 22;
        let aminoglycosides = [23, 24, 25]; // gentamicin, tobramycin, amikacin
        let fluoroquinolones = [26, 27, 28, 29]; // cipro, levo, moxi, ofloxacin
        let moxifloxacin_idx = 28;
        let tetracyclines = [30, 31, 32, 33]; // tetracycline, doxycycline, minocycline, tigecycline
        let minocycline_idx = 32;
        let glycopeptides = [34, 35, 36]; // vancomycin, teicoplanin, dalbavancin
        let vancomycin_idx = 34;
        let oxazolidinones = [37, 38]; // linezolid, tedizolid
        let daptomycin_idx = 39;
        let trim_sulf_idx = 41;
        let chloramphenicol_idx = 42;
        let nitrofurantoin_idx = 43;
        let fosfomycin_idx = 44;
        let metronidazole_idx = 47;
        let fidaxomicin_idx = 48;
        let furazolidone_idx = 49;
        let rifampicin_idx = 50;
        let blbli_combinations = [51, 52, 53, 54, 55, 56]; // β-lactam/β-lactamase inhibitor combos
        let colistin_idx = 57;
        let aztreonam_avibactam_idx = 59;
        let cefixime_idx = 60;
        
        // --- CNS (syndrome 6) - Blood-brain barrier severely limits most drugs ---
        for &d in &penicillins { drug_penetration[6][d] = 0.15; } // Poor unless meningeal inflammation
        for &d in &oral_cephalosporins { drug_penetration[6][d] = 0.05; } // Very poor
        for &d in &iv_cephalosporins { drug_penetration[6][d] = 0.20; } // Ceftriaxone better ~30%
        drug_penetration[6][ceftriaxone_idx] = 0.35; // Ceftriaxone - best CSF penetration among cephalosporins
        for &d in &carbapenems { drug_penetration[6][d] = 0.25; } // Meropenem preferred for CNS
        drug_penetration[6][meropenem_idx] = 0.35; // Meropenem - good CNS penetration
        drug_penetration[6][aztreonam_idx] = 0.10; // Poor
        drug_penetration[6][aztreonam_avibactam_idx] = 0.10; // Similar tissue distribution to aztreonam
        drug_penetration[6][cefixime_idx] = 0.10; // Oral 3G cephalosporin - poor CNS penetration
        for &d in &macrolides { drug_penetration[6][d] = 0.15; } // Poor
        drug_penetration[6][clindamycin_idx] = 0.15; // Poor
        for &d in &aminoglycosides { drug_penetration[6][d] = 0.05; } // Very poor - aminoglycosides don't cross BBB
        for &d in &fluoroquinolones { drug_penetration[6][d] = 0.50; } // Good - lipophilic
        drug_penetration[6][moxifloxacin_idx] = 0.60; // Moxifloxacin - excellent CNS penetration
        for &d in &tetracyclines { drug_penetration[6][d] = 0.25; } // Moderate
        drug_penetration[6][minocycline_idx] = 0.40; // Minocycline - good lipophilicity
        for &d in &glycopeptides { drug_penetration[6][d] = 0.15; } // Poor unless inflamed meninges
        for &d in &oxazolidinones { drug_penetration[6][d] = 0.70; } // Linezolid excellent CNS
        drug_penetration[6][trim_sulf_idx] = 0.50; // Good - used for CNS toxoplasmosis
        drug_penetration[6][chloramphenicol_idx] = 0.70; // Excellent CNS penetration
        drug_penetration[6][nitrofurantoin_idx] = 0.05; // No CNS penetration
        drug_penetration[6][metronidazole_idx] = 0.80; // Excellent - used for brain abscess
        drug_penetration[6][rifampicin_idx] = 0.50; // Good
        for &d in &blbli_combinations { drug_penetration[6][d] = 0.15; } // Similar to parent β-lactam
        drug_penetration[6][colistin_idx] = 0.05; // Very poor
        drug_penetration[6][daptomycin_idx] = 0.05;
        drug_penetration[6][fosfomycin_idx] = 0.3;
        drug_penetration[6][fidaxomicin_idx] = 0.0;
        
        // --- Bone/joint (syndrome 9) - Poor vascularity, biofilm ---
        for &d in &penicillins { drug_penetration[9][d] = 0.40; }
        for &d in &oral_cephalosporins { drug_penetration[9][d] = 0.30; }
        for &d in &iv_cephalosporins { drug_penetration[9][d] = 0.45; }
        for &d in &carbapenems { drug_penetration[9][d] = 0.50; }
        drug_penetration[9][aztreonam_idx] = 0.35;
        drug_penetration[9][aztreonam_avibactam_idx] = 0.35;
        drug_penetration[9][cefixime_idx] = 0.40;
        for &d in &macrolides { drug_penetration[9][d] = 0.40; }
        drug_penetration[9][clindamycin_idx] = 0.60; // Clindamycin good bone penetration
        for &d in &aminoglycosides { drug_penetration[9][d] = 0.25; } // Poor bone penetration
        for &d in &fluoroquinolones { drug_penetration[9][d] = 0.70; } // Excellent bone penetration
        for &d in &tetracyclines { drug_penetration[9][d] = 0.50; }
        for &d in &glycopeptides { drug_penetration[9][d] = 0.35; } // Vancomycin moderate
        for &d in &oxazolidinones { drug_penetration[9][d] = 0.75; } // Linezolid excellent
        drug_penetration[9][trim_sulf_idx] = 0.55;
        drug_penetration[9][chloramphenicol_idx] = 0.50;
        drug_penetration[9][nitrofurantoin_idx] = 0.10; // Poor systemic distribution
        drug_penetration[9][metronidazole_idx] = 0.55;
        drug_penetration[9][rifampicin_idx] = 0.80; // Excellent - key for osteomyelitis
        for &d in &blbli_combinations { drug_penetration[9][d] = 0.40; }
        drug_penetration[9][colistin_idx] = 0.20;
        drug_penetration[9][daptomycin_idx] = 0.5;
        drug_penetration[9][fosfomycin_idx] = 0.6;
        drug_penetration[9][fidaxomicin_idx] = 0.0;
        
        // --- Intra-abdominal (syndrome 5) - Abscess cavities, acidic pH ---
        for &d in &penicillins { drug_penetration[5][d] = 0.60; }
        for &d in &oral_cephalosporins { drug_penetration[5][d] = 0.45; }
        for &d in &iv_cephalosporins { drug_penetration[5][d] = 0.65; }
        for &d in &carbapenems { drug_penetration[5][d] = 0.75; }
        drug_penetration[5][aztreonam_idx] = 0.55;
        drug_penetration[5][aztreonam_avibactam_idx] = 0.55;
        drug_penetration[5][cefixime_idx] = 0.55;
        for &d in &macrolides { drug_penetration[5][d] = 0.50; }
        drug_penetration[5][clindamycin_idx] = 0.65;
        for &d in &aminoglycosides { drug_penetration[5][d] = 0.30; } // Inactivated at acidic pH
        for &d in &fluoroquinolones { drug_penetration[5][d] = 0.75; }
        for &d in &tetracyclines { drug_penetration[5][d] = 0.55; }
        for &d in &glycopeptides { drug_penetration[5][d] = 0.45; }
        for &d in &oxazolidinones { drug_penetration[5][d] = 0.70; }
        drug_penetration[5][trim_sulf_idx] = 0.60;
        drug_penetration[5][chloramphenicol_idx] = 0.60;
        drug_penetration[5][nitrofurantoin_idx] = 0.15;
        drug_penetration[5][metronidazole_idx] = 0.90; // Excellent for anaerobic abscesses
        drug_penetration[5][rifampicin_idx] = 0.65;
        for &d in &blbli_combinations { drug_penetration[5][d] = 0.65; }
        drug_penetration[5][colistin_idx] = 0.35;
        drug_penetration[5][daptomycin_idx] = 0.6;
        drug_penetration[5][fosfomycin_idx] = 0.5;
        drug_penetration[5][fidaxomicin_idx] = 0.05;
        
        // --- UTI (syndrome 1) - Renal excretion concentrates many drugs ---
        for &d in &penicillins { drug_penetration[1][d] = 0.80; }
        for &d in &oral_cephalosporins { drug_penetration[1][d] = 0.85; }
        for &d in &iv_cephalosporins { drug_penetration[1][d] = 0.85; }
        for &d in &carbapenems { drug_penetration[1][d] = 0.85; }
        drug_penetration[1][aztreonam_idx] = 0.80;
        drug_penetration[1][aztreonam_avibactam_idx] = 0.80;
        drug_penetration[1][cefixime_idx] = 0.90;
        for &d in &macrolides { drug_penetration[1][d] = 0.40; } // Poor urinary excretion
        drug_penetration[1][clindamycin_idx] = 0.30; // Poor for UTI
        for &d in &aminoglycosides { drug_penetration[1][d] = 0.75; } // Renally excreted
        for &d in &fluoroquinolones { drug_penetration[1][d] = 1.0; } // Excellent urinary concentration
        for &d in &tetracyclines { drug_penetration[1][d] = 0.50; }
        for &d in &glycopeptides { drug_penetration[1][d] = 0.60; }
        for &d in &oxazolidinones { drug_penetration[1][d] = 0.70; }
        drug_penetration[1][trim_sulf_idx] = 1.0; // Excellent - first-line UTI
        drug_penetration[1][chloramphenicol_idx] = 0.40;
        drug_penetration[1][nitrofurantoin_idx] = 1.0; // Concentrated in urine - first-line UTI
        drug_penetration[1][metronidazole_idx] = 0.50;
        drug_penetration[1][rifampicin_idx] = 0.40;
        for &d in &blbli_combinations { drug_penetration[1][d] = 0.80; }
        drug_penetration[1][colistin_idx] = 0.70;
        drug_penetration[1][daptomycin_idx] = 0.1;
        drug_penetration[1][fosfomycin_idx] = 1.0;
        drug_penetration[1][fidaxomicin_idx] = 0.0;
        
        // --- Skin/soft tissue (syndrome 2) - Generally good penetration ---
        for &d in &penicillins { drug_penetration[2][d] = 0.85; }
        for &d in &oral_cephalosporins { drug_penetration[2][d] = 0.80; }
        for &d in &iv_cephalosporins { drug_penetration[2][d] = 0.85; }
        for &d in &carbapenems { drug_penetration[2][d] = 0.85; }
        drug_penetration[2][aztreonam_idx] = 0.75;
        drug_penetration[2][aztreonam_avibactam_idx] = 0.75;
        drug_penetration[2][cefixime_idx] = 0.75;
        for &d in &macrolides { drug_penetration[2][d] = 0.80; }
        drug_penetration[2][clindamycin_idx] = 0.85; // Excellent skin penetration
        for &d in &aminoglycosides { drug_penetration[2][d] = 0.60; }
        for &d in &fluoroquinolones { drug_penetration[2][d] = 0.90; }
        for &d in &tetracyclines { drug_penetration[2][d] = 0.80; }
        for &d in &glycopeptides { drug_penetration[2][d] = 0.75; }
        for &d in &oxazolidinones { drug_penetration[2][d] = 0.90; }
        drug_penetration[2][trim_sulf_idx] = 0.80;
        drug_penetration[2][chloramphenicol_idx] = 0.70;
        drug_penetration[2][nitrofurantoin_idx] = 0.20; // Poor systemic distribution
        drug_penetration[2][metronidazole_idx] = 0.75;
        drug_penetration[2][rifampicin_idx] = 0.80;
        for &d in &blbli_combinations { drug_penetration[2][d] = 0.85; }
        drug_penetration[2][colistin_idx] = 0.50;
        drug_penetration[2][daptomycin_idx] = 0.95;
        drug_penetration[2][fosfomycin_idx] = 0.5;
        drug_penetration[2][fidaxomicin_idx] = 0.0;
        
        // --- Respiratory (syndrome 3) - ELF penetration varies significantly ---
        for &d in &penicillins { drug_penetration[3][d] = 0.65; }
        for &d in &oral_cephalosporins { drug_penetration[3][d] = 0.55; }
        for &d in &iv_cephalosporins { drug_penetration[3][d] = 0.70; }
        for &d in &carbapenems { drug_penetration[3][d] = 0.75; }
        drug_penetration[3][aztreonam_idx] = 0.60;
        drug_penetration[3][aztreonam_avibactam_idx] = 0.60;
        drug_penetration[3][cefixime_idx] = 0.60;
        for &d in &macrolides { drug_penetration[3][d] = 0.95; } // Excellent lung tissue concentration
        drug_penetration[3][clindamycin_idx] = 0.75;
        for &d in &aminoglycosides { drug_penetration[3][d] = 0.40; } // Poor ELF penetration
        for &d in &fluoroquinolones { drug_penetration[3][d] = 0.95; } // Excellent respiratory penetration
        for &d in &tetracyclines { drug_penetration[3][d] = 0.70; }
        for &d in &glycopeptides { drug_penetration[3][d] = 0.50; } // Vancomycin poor ELF
        for &d in &oxazolidinones { drug_penetration[3][d] = 0.90; } // Excellent
        drug_penetration[3][trim_sulf_idx] = 0.80;
        drug_penetration[3][chloramphenicol_idx] = 0.70;
        drug_penetration[3][nitrofurantoin_idx] = 0.15;
        drug_penetration[3][metronidazole_idx] = 0.60;
        drug_penetration[3][rifampicin_idx] = 0.85; // Excellent - TB treatment
        for &d in &blbli_combinations { drug_penetration[3][d] = 0.65; }
        drug_penetration[3][colistin_idx] = 0.30; // Poor systemic, but used inhaled
        drug_penetration[3][daptomycin_idx] = 0.0; // Inactivated by pulmonary surfactant
        drug_penetration[3][fosfomycin_idx] = 0.4;
        drug_penetration[3][fidaxomicin_idx] = 0.0;
        
        // --- Bloodstream (syndrome 4) - Direct access, reference compartment ---
        // All drugs get 1.0 (full serum concentration) by default initialization
        
        // --- GI (syndrome 7) - Luminal vs systemic varies ---
        for &d in &penicillins { drug_penetration[7][d] = 0.55; }
        for &d in &oral_cephalosporins { drug_penetration[7][d] = 0.50; }
        for &d in &iv_cephalosporins { drug_penetration[7][d] = 0.60; }
        for &d in &carbapenems { drug_penetration[7][d] = 0.65; }
        drug_penetration[7][aztreonam_idx] = 0.50;
        drug_penetration[7][aztreonam_avibactam_idx] = 0.50;
        drug_penetration[7][cefixime_idx] = 0.55;
        for &d in &macrolides { drug_penetration[7][d] = 0.70; } // Good GI tissue penetration
        drug_penetration[7][clindamycin_idx] = 0.65;
        for &d in &aminoglycosides { drug_penetration[7][d] = 0.40; }
        for &d in &fluoroquinolones { drug_penetration[7][d] = 0.85; } // Excellent
        for &d in &tetracyclines { drug_penetration[7][d] = 0.60; }
        for &d in &glycopeptides { drug_penetration[7][d] = 0.35; } // Poor oral absorption but good for C.diff
        drug_penetration[7][vancomycin_idx] = 0.90; // Oral vancomycin excellent for C.diff (luminal)
        for &d in &oxazolidinones { drug_penetration[7][d] = 0.75; }
        drug_penetration[7][trim_sulf_idx] = 0.70;
        drug_penetration[7][chloramphenicol_idx] = 0.65;
        drug_penetration[7][nitrofurantoin_idx] = 0.25;
        drug_penetration[7][metronidazole_idx] = 0.95; // Excellent - C.diff, amebiasis
        drug_penetration[7][furazolidone_idx] = 0.90; // Furazolidone - GI specific
        drug_penetration[7][rifampicin_idx] = 0.60;
        for &d in &blbli_combinations { drug_penetration[7][d] = 0.55; }
        drug_penetration[7][colistin_idx] = 0.40;
        drug_penetration[7][daptomycin_idx] = 0.3;
        drug_penetration[7][fosfomycin_idx] = 0.4;
        drug_penetration[7][fidaxomicin_idx] = 1.0; // Excellent - C.diff specific (luminal)
        
        // --- Genital (syndrome 8) - Prostate barrier significant for males ---
        for &d in &penicillins { drug_penetration[8][d] = 0.55; }
        for &d in &oral_cephalosporins { drug_penetration[8][d] = 0.45; }
        for &d in &iv_cephalosporins { drug_penetration[8][d] = 0.55; }
        for &d in &carbapenems { drug_penetration[8][d] = 0.60; }
        drug_penetration[8][aztreonam_idx] = 0.45;
        drug_penetration[8][aztreonam_avibactam_idx] = 0.45;
        drug_penetration[8][cefixime_idx] = 0.50;
        for &d in &macrolides { drug_penetration[8][d] = 0.75; } // Good tissue penetration
        drug_penetration[8][clindamycin_idx] = 0.60;
        for &d in &aminoglycosides { drug_penetration[8][d] = 0.35; } // Poor prostate
        for &d in &fluoroquinolones { drug_penetration[8][d] = 0.90; } // Excellent prostatic penetration
        for &d in &tetracyclines { drug_penetration[8][d] = 0.75; } // Good - used for STIs
        for &d in &glycopeptides { drug_penetration[8][d] = 0.40; }
        for &d in &oxazolidinones { drug_penetration[8][d] = 0.70; }
        drug_penetration[8][trim_sulf_idx] = 0.80; // Good prostatic penetration
        drug_penetration[8][chloramphenicol_idx] = 0.55;
        drug_penetration[8][nitrofurantoin_idx] = 0.30;
        drug_penetration[8][metronidazole_idx] = 0.80; // Good - trichomoniasis, BV
        drug_penetration[8][rifampicin_idx] = 0.60;
        for &d in &blbli_combinations { drug_penetration[8][d] = 0.55; }
        drug_penetration[8][colistin_idx] = 0.30;
        drug_penetration[8][daptomycin_idx] = 0.4;
        drug_penetration[8][fosfomycin_idx] = 0.5;
        drug_penetration[8][fidaxomicin_idx] = 0.0;
        
        // --- Other (syndrome 10) - Use moderate defaults ---
        // Keep at 1.0 (default) or slightly reduced
        for &d in &aminoglycosides { drug_penetration[10][d] = 0.70; }
        drug_penetration[10][nitrofurantoin_idx] = 0.30;
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
            .unwrap_or(0.1)
    }

    #[inline]
    pub fn bacteria_growth_multiplier(&self, syndrome_id: usize) -> f64 {
        self.bacteria_growth_multiplier
            .get(syndrome_id)
            .copied()
            .unwrap_or(1.0)
    }

    /// Get drug penetration factor for a specific syndrome and drug
    /// Returns value between 0.0 and 1.0 representing fraction of serum concentration
    /// achieved at the infection site
    #[inline]
    pub fn drug_penetration(&self, syndrome_id: usize, drug_idx: usize) -> f64 {
        self.drug_penetration
            .get(syndrome_id)
            .and_then(|drugs| drugs.get(drug_idx))
            .copied()
            .unwrap_or(1.0)
    }
}

impl ImmunodeficiencyParameters {
    fn from_map(map: &HashMap<String, f64>) -> Self {
        let startup_seed_fraction = get_or_default(
            map,
            "immunosuppression_startup_seed_fraction",
            0.05,
        );
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
            startup_seed_fraction,
            temporary_onset_rate_per_day,
            temporary_recovery_rate_per_day,
            chronic_onset_rate_per_day,
            chronic_recovery_rate_per_day,
            chronic_probability_age_bands,
        }
    }

    #[inline]
    pub fn startup_seed_fraction(&self) -> f64 {
        self.startup_seed_fraction
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
    pub initial_infection_level: Vec<f64>,
    pub base_bacteria_level_change: Vec<f64>,
    pub max_level: Vec<f64>,
    // Logistic model for symptom onset probability
    // P(symptoms) = 1 / (1 + exp(-log_odds)) where log_odds = base + level_effect
    pub symptom_onset_base_log_odds: Vec<f64>,
    pub symptom_onset_threshold_level: Vec<f64>,
    pub symptom_onset_delay_days: Vec<f64>,
    pub symptom_onset_log_odds_per_level_unit: Vec<f64>,
    pub mechanismless_resistance_reversion_rate: Vec<f64>,
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
        let mut initial_infection_level = Vec::with_capacity(num_bacteria);
        let mut base_bacteria_level_change = Vec::with_capacity(num_bacteria);
        let mut max_level = Vec::with_capacity(num_bacteria);
        let mut symptom_onset_base_log_odds = Vec::with_capacity(num_bacteria);
        let mut symptom_onset_threshold_level = Vec::with_capacity(num_bacteria);
        let mut symptom_onset_delay_days = Vec::with_capacity(num_bacteria);
        let mut symptom_onset_log_odds_per_level_unit = Vec::with_capacity(num_bacteria);
        let mut mechanismless_resistance_reversion_rate = Vec::with_capacity(num_bacteria);
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
            // Symptom onset logistic parameters
            // Default base log-odds: log(0.15/0.85) ≈ -1.73 (equivalent to 15% daily probability)
            symptom_onset_base_log_odds.push(get_or_default(
                map,
                &format!("{}_symptom_onset_base_log_odds", prefix),
                -1.73,
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
            // Log-odds increase per unit of bacteria level above threshold
            symptom_onset_log_odds_per_level_unit.push(get_or_default(
                map,
                &format!("{}_symptom_onset_log_odds_per_level_unit", prefix),
                0.5, // Each unit above threshold adds 0.5 log-odds (~1.6x odds multiplier)
            ));
            mechanismless_resistance_reversion_rate.push(get_or_default(
                map,
                &format!("{}_mechanismless_resistance_reversion_rate", prefix),
                get_or_default(
                    map,
                    "mechanismless_resistance_reversion_rate",
                    0.0004,
                ),
            ));
            microbiome_vs_infection_log_odds.push(get_or_default(
                map,
                &format!("{}_log_odds_microbiome_vs_infection", prefix),
                get_or_default(map, "log_odds_microbiome_vs_infection", 3.0), // Iter7: was -6.0 (99.75% infection!), changed to 3.0 (~5% infection)
            ));
            let cessation_key = format!("{}_drug_cessation_probability", prefix.to_lowercase());
            let default_cessation = get_or_default(map, "random_drug_cessation_probability", 0.001);
            drug_cessation_probability.push(get_or_default(map, &cessation_key, default_cessation));
            let recognition_key = format!("{}_treatment_recognition_year", prefix);
            treatment_recognition_year.push(map.get(&recognition_key).copied());
            sepsis_baseline_log_odds.push(get_or_default(
                map,
                &format!("{}_sepsis_baseline_log_odds", prefix),
                get_or_default(map, "sepsis_baseline_log_odds", -14.0),
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
            initial_infection_level,
            base_bacteria_level_change,
            max_level,
            symptom_onset_base_log_odds,
            symptom_onset_threshold_level,
            symptom_onset_delay_days,
            symptom_onset_log_odds_per_level_unit,
            mechanismless_resistance_reversion_rate,
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
    pub fn initial_infection_level(&self, bacteria_idx: usize) -> f64 {
        self.initial_infection_level[bacteria_idx]
    }

    #[inline]
    pub fn max_level(&self, bacteria_idx: usize) -> f64 {
        self.max_level[bacteria_idx]
    }

    #[inline]
    pub fn symptom_onset_base_log_odds(&self, bacteria_idx: usize) -> f64 {
        self.symptom_onset_base_log_odds[bacteria_idx]
    }

    #[inline]
    pub fn symptom_onset_threshold_level(&self, bacteria_idx: usize) -> f64 {
        self.symptom_onset_threshold_level[bacteria_idx]
    }

    #[inline]
    pub fn mechanismless_resistance_reversion_rate(&self, bacteria_idx: usize) -> f64 {
        self.mechanismless_resistance_reversion_rate[bacteria_idx]
    }

    #[inline]
    pub fn symptom_onset_delay_days(&self, bacteria_idx: usize) -> f64 {
        self.symptom_onset_delay_days[bacteria_idx]
    }

    #[inline]
    pub fn symptom_onset_log_odds_per_level_unit(&self, bacteria_idx: usize) -> f64 {
        self.symptom_onset_log_odds_per_level_unit[bacteria_idx]
    }
}

// ---------------- 8) Clearance, acquisition, and age tables ----------------
// immune_clearance
// ClearanceParameters encodes immune-mediated clearance of active infections.
// Logistic model: P(clearance) = 1 / (1 + exp(-log_odds))
// log_odds = base + bacteria_effect + age_effect + immuno_effect + level_effect
#[derive(Debug)]
pub struct ClearanceParameters {
    _base_delay_days: f64,
    // Logistic model parameters
    base_clearance_log_odds: f64,
    _per_bacteria_delay_days: Vec<Option<f64>>,
    per_bacteria_log_odds_adjustment: Vec<f64>,
    age_log_odds_adjustments: [f64; AGE_BUCKET_COUNT],
    immunodeficient_log_odds_adjustment: f64,
    level_log_odds_per_unit: f64, // negative: higher level → lower clearance
}

impl ClearanceParameters {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize) -> Self {
        // Base post-infection delay before clearance can occur
        let base_delay_days = get_or_default(map, "default_clearance_delay_days", 3.0);
        // Logistic model: base log-odds for clearance probability
        // Default: log(0.015/0.985) ≈ -4.2 (equivalent to 1.5% daily clearance)
        let base_clearance_log_odds = get_or_default(map, "default_clearance_base_log_odds", -4.2);

        // Age-specific log-odds adjustments (additive in log-odds space)
        let mut age_log_odds_adjustments = [0.0; AGE_BUCKET_COUNT];
        for (idx, category) in AGE_BUCKETS.iter().enumerate() {
            let key = format!("clearance_age_log_odds_{}", category.label());
            age_log_odds_adjustments[idx] = get_or_default(map, &key, 0.0);
        }

        let mut per_bacteria_delay_days = Vec::with_capacity(num_bacteria);
        let mut per_bacteria_log_odds_adjustment = Vec::with_capacity(num_bacteria);
        for &bacteria in BACTERIA_LIST.iter() {
            per_bacteria_delay_days.push(
                map.get(&format!("{}_clearance_delay_days", bacteria))
                    .copied(),
            );
            // Bacteria-specific log-odds adjustment (additive)
            per_bacteria_log_odds_adjustment.push(get_or_default(
                map,
                &format!("{}_clearance_log_odds_adjustment", bacteria),
                0.0,
            ));
        }

        // Immunodeficiency effect: negative log-odds adjustment (harder to clear)
        // Default: log(0.5) ≈ -0.69, equivalent to 50% multiplier
        let immunodeficient_log_odds_adjustment =
            get_or_default(map, "clearance_immunodeficient_log_odds", -0.69);
        // Level effect: higher bacteria level → harder to clear (negative coefficient)
        let level_log_odds_per_unit = get_or_default(map, "clearance_level_log_odds_per_unit", -0.3);

        ClearanceParameters {
            _base_delay_days: base_delay_days,
            base_clearance_log_odds,
            _per_bacteria_delay_days: per_bacteria_delay_days,
            per_bacteria_log_odds_adjustment: per_bacteria_log_odds_adjustment,
            age_log_odds_adjustments,
            immunodeficient_log_odds_adjustment,
            level_log_odds_per_unit,
        }
    }

    #[inline]
    #[allow(dead_code)]
    pub fn delay_days(&self, bacteria_idx: usize) -> f64 {
        self._per_bacteria_delay_days[bacteria_idx]
            .unwrap_or(self._base_delay_days)
            .max(0.0)
    }

    #[inline]
    #[allow(dead_code)]
    pub fn base_log_odds(&self) -> f64 {
        self.base_clearance_log_odds
    }

    #[inline]
    pub fn bacteria_log_odds_adjustment(&self, bacteria_idx: usize) -> f64 {
        self.per_bacteria_log_odds_adjustment[bacteria_idx]
    }

    #[inline]
    pub fn age_log_odds_adjustment(&self, age_days: i32) -> f64 {
        let idx = AgeCategoryParameters::age_category_index(age_days)
            .min(self.age_log_odds_adjustments.len() - 1);
        self.age_log_odds_adjustments[idx]
    }

    #[inline]
    #[allow(dead_code)]
    pub fn immunodeficient_log_odds_adjustment(&self) -> f64 {
        self.immunodeficient_log_odds_adjustment
    }

    #[inline]
    pub fn level_log_odds_effect(&self, level: f64) -> f64 {
        // Higher level → more negative log-odds (harder to clear)
        level.max(0.0) * self.level_log_odds_per_unit
    }

    #[inline]
    pub fn hazard_for(
        &self,
        bacteria_idx: usize,
        age_days: i32,
        is_immunodeficient: bool,
        level: f64,
        duration_days: u32,
    ) -> f64 {
        // Logistic model: P(clearance) = 1 / (1 + exp(-log_odds))
        // log_odds = base + bacteria + age + immuno + level
        let mut log_odds = self.base_clearance_log_odds;
        log_odds += self.bacteria_log_odds_adjustment(bacteria_idx);
        log_odds += self.age_log_odds_adjustment(age_days);
        
        if is_immunodeficient {
            log_odds += self.immunodeficient_log_odds_adjustment;
        }
        
        log_odds += self.level_log_odds_effect(level);

        // --- immune_clearance ---
        // Basic innate immunity is the baseline.
        // As infection persists, adaptive immunity recruits (T-cells, antibodies),
        // increasing the clearance log-odds linearly over time, which creates a sigmoidal
        // rise in clearance probability.
        let adaptive_recruit_slope = 0.25; 
        log_odds += (duration_days as f64) * adaptive_recruit_slope;
        
        // Logistic transformation
        1.0 / (1.0 + (-log_odds).exp())
    }
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct DrugBacteriaMatrix {
    pub potency_when_no_r: Vec<f64>,
    pub initiation_multiplier: Vec<f64>,
    pub mic_lt2_threshold: Vec<f64>,
    num_bacteria: usize,
    num_drugs: usize,
}

impl DrugBacteriaMatrix {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize, num_drugs: usize) -> Self {
        let mut potency_when_no_r = Vec::with_capacity(num_bacteria * num_drugs);
        let mut initiation_multiplier = Vec::with_capacity(num_bacteria * num_drugs);
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
                let threshold = 2.0 * potency.min(1.0) - 1.0;
                mic_lt2_threshold.push(threshold);
            }
        }

        DrugBacteriaMatrix {
            potency_when_no_r,
            initiation_multiplier,
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
            default_age_log_odds[idx] =
                get_or_default(map, &format!("default_log_odds_{}", category.label()), 0.0);
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

#[derive(Clone, Copy, PartialEq, Eq)]
enum PlasmidPool {
    None,
    GramPositive,
    EntericGramNegative,
    RespiratoryGramNegative,
    Anaerobe,
}

fn group_has_structural_hgt_exclusion(group: BacteriaGroup) -> bool {
    matches!(
        group,
        BacteriaGroup::Spirochete | BacteriaGroup::Helicobacter | BacteriaGroup::Mycobacteria
    )
}

fn pool_for_group(group: BacteriaGroup) -> PlasmidPool {
    match group {
        BacteriaGroup::GramPositive => PlasmidPool::GramPositive,
        BacteriaGroup::Enterobacterales
        | BacteriaGroup::NonFermenter
        | BacteriaGroup::EntericPathogen => PlasmidPool::EntericGramNegative,
        BacteriaGroup::Fastidious => PlasmidPool::RespiratoryGramNegative,
        BacteriaGroup::Anaerobe => PlasmidPool::Anaerobe,
        BacteriaGroup::Spirochete | BacteriaGroup::Helicobacter | BacteriaGroup::Mycobacteria => {
            PlasmidPool::None
        }
    }
}

// ^^^
fn default_hgt_probability(donor_idx: usize, recipient_idx: usize) -> f64 {
    let donor_group = BACTERIA_GROUPS
        .get(donor_idx)
        .copied()
        .unwrap_or(BacteriaGroup::Enterobacterales);
    let recipient_group = BACTERIA_GROUPS
        .get(recipient_idx)
        .copied()
        .unwrap_or(BacteriaGroup::Enterobacterales);

    if group_has_structural_hgt_exclusion(donor_group)
        || group_has_structural_hgt_exclusion(recipient_group)
    {
        return 0.0;
    }

    let donor_pool = pool_for_group(donor_group);
    let recipient_pool = pool_for_group(recipient_group);

    if donor_pool == PlasmidPool::None || recipient_pool == PlasmidPool::None {
        return 0.0;
    }

    let same_group = donor_group == recipient_group;

    // ***  ^^^ default_hgt_probability
    match (donor_pool, recipient_pool) {
        (PlasmidPool::GramPositive, PlasmidPool::GramPositive) => {
            if same_group { 0.000_000_001 } else { 0.000_000_000_1 }
        }
        (PlasmidPool::EntericGramNegative, PlasmidPool::EntericGramNegative) => {
            if same_group { 0.000_000_001 } else { 0.000_000_000_1 }
        }
        (PlasmidPool::RespiratoryGramNegative, PlasmidPool::RespiratoryGramNegative) => {
            if same_group { 0.000_000_001 } else { 0.000_000_000_1 }
        }
        // Cross-pool exceptions exist, but stay below within-pool mismatches.
        (PlasmidPool::EntericGramNegative, PlasmidPool::RespiratoryGramNegative)
        | (PlasmidPool::RespiratoryGramNegative, PlasmidPool::EntericGramNegative) => 0.000_000_000_03,
        (PlasmidPool::Anaerobe, PlasmidPool::EntericGramNegative)
        | (PlasmidPool::EntericGramNegative, PlasmidPool::Anaerobe) => 0.000_000_000_03,
        // Anaerobes already receive additional ecological opportunity via the gut-compartment multiplier.
        (PlasmidPool::Anaerobe, PlasmidPool::Anaerobe) => 0.000_000_001,
        // Remaining cross-pool transfers are treated as biologically negligible.
        _ => 0.0,
    }
}

#[derive(Debug)]
pub struct ResistanceMechanismParameters {
    /// Per-mechanism per-drug-class enhancement multipliers.
    /// Indexed as [mechanism_idx * DrugClass::NUM_CLASSES + drug_class_idx].
    /// Represents how much resistance a given mechanism confers against each drug class.
    pub enhancement_multiplier: Vec<f64>,
    pub reversion_rate: Vec<f64>,
}

impl ResistanceMechanismParameters {
    fn from_map(map: &HashMap<String, f64>) -> Self {
        use crate::simulation::population::DrugClass;

        let num_mechanisms = ResistanceMechanism::all().len();
        let num_classes = DrugClass::NUM_CLASSES;
        let mut enhancement_multiplier = vec![0.0; num_mechanisms * num_classes];
        let mut reversion_rate = Vec::with_capacity(num_mechanisms);

        for (mech_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
            let name = mechanism.as_str();

            // Load the legacy single-value as fallback default
            let legacy_value = get_or_default(
                map,
                &format!("resistance_mechanism_{}_enhancement_multiplier", name),
                0.0,
            );

            // For each drug class, look for a class-specific key first, fall back to legacy
            for drug_class in DrugClass::all() {
                let class_key = format!(
                    "resistance_mechanism_{}_enhancement_{}", name, drug_class.as_str()
                );
                let value = if let Some(&v) = map.get(&class_key) {
                    v
                } else {
                    legacy_value
                };
                enhancement_multiplier[mech_idx * num_classes + drug_class.index()] = value;
            }

            reversion_rate.push(get_or_default(
                map,
                &format!("resistance_mechanism_{}_reversion_rate", name),
                0.0001,
            ));
        }

        ResistanceMechanismParameters {
            enhancement_multiplier,
            reversion_rate,
        }
    }

    /// Get the enhancement multiplier for a specific mechanism against a specific drug class
    #[inline]
    pub fn enhancement_multiplier(&self, mechanism_idx: usize, drug_class_idx: usize) -> f64 {
        self.enhancement_multiplier[mechanism_idx * crate::simulation::population::DrugClass::NUM_CLASSES + drug_class_idx]
    }

    #[inline]
    pub fn reversion_rate(&self, mechanism_idx: usize) -> f64 {
        self.reversion_rate[mechanism_idx]
    }
}

#[derive(Debug)]
pub struct BacteriaMechanismEmergenceRates {
    values: Vec<f64>,
    num_mechanisms: usize,
}
 
impl BacteriaMechanismEmergenceRates {
    fn from_map(map: &HashMap<String, f64>, num_bacteria: usize) -> Self {
        let mechanisms = ResistanceMechanism::all();
        let num_mechanisms = mechanisms.len();
        let mut values = Vec::with_capacity(num_bacteria * num_mechanisms);

        for &bacteria in BACTERIA_LIST.iter() {
            for mechanism in mechanisms {
                let key = format!(
                    "bacteria_{}_mechanism_{}_emergence_rate",
                    bacteria,
                    mechanism.as_str()
                );
                values.push(get_or_default(map, &key, 0.0));
            }
        }

        Self {
            values,
            num_mechanisms,
        }
    }

    #[inline]
    pub fn rate(&self, bacteria_idx: usize, mechanism_idx: usize) -> f64 {
        let offset = bacteria_idx * self.num_mechanisms + mechanism_idx;
        self.values[offset]
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
// This embedded table only covers the legacy 52-drug matrix. Late-added drugs are
// assigned through the rule-based potency overrides below.
const POTENCY_EMBEDDED_HEADER: [&str; 52] = [
    "sulfanilamide",
    "penicillin_g",
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
    "chloramphenicol",
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.600000),
            Some(0.700000),
            Some(0.100000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.100000),
            Some(0.100000),
            Some(0.800000),
            Some(0.750000),
            Some(0.100000),
            Some(0.100000),
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
            Some(0.100000),
            Some(0.100000),
            Some(0.750000),
            Some(0.700000),
            Some(0.100000),
            Some(0.100000),
            Some(0.600000),
            Some(0.500000),
            Some(0.800000),
            Some(0.850000),
            Some(0.400000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.800000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.850000),
            Some(0.800000),
            Some(0.700000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.100000),
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
            Some(0.000000),
        ],
    ),
    (
        "enterococcus_faecium",
        [
            Some(0.100000),
            Some(0.100000),
            Some(0.300000),
            Some(0.300000),
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
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.900000),
            Some(0.850000),
            Some(0.950000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.100000),
            Some(0.100000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.900000),
            Some(0.850000),
            Some(0.800000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
        "p_stuartii",
        [
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.350000),
            Some(0.300000),
            Some(0.100000),
            Some(0.100000),
            Some(0.200000),
            Some(0.450000),
            Some(0.750000),
            Some(0.850000),
            Some(0.200000),
            Some(0.900000),
            Some(0.900000),
            Some(0.850000),
            Some(0.650000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.450000),
            Some(0.500000),
            Some(0.750000),
            Some(0.550000),
            Some(0.600000),
            Some(0.600000),
            Some(0.550000),
            Some(0.200000),
            Some(0.300000),
            Some(0.350000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.300000),
            Some(0.250000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.200000),
            Some(0.750000),
            Some(0.200000),
            Some(0.450000),
            Some(0.900000),
            Some(0.950000),
            Some(0.050000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.100000),
            Some(0.100000),
            Some(0.050000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
        "stenotrophomonas_maltophilia",
        [
            Some(0.600000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.200000),
            Some(0.250000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.350000),
            Some(0.150000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.100000),
            Some(0.100000),
            Some(0.400000),
            Some(0.750000),
            Some(0.800000),
            Some(0.600000),
            Some(0.350000),
            Some(0.600000),
            Some(0.850000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.950000),
            Some(0.400000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.200000),
            Some(0.050000),
            Some(0.300000),
            Some(0.050000),
            Some(0.700000),
            Some(0.400000),
            Some(0.050000),
            Some(0.050000),
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
            Some(0.000000),
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
            Some(0.000000),
        ],
    ),
    (
        "staphylococcus_epidermidis",
        [
            Some(0.100000),
            Some(0.150000),
            Some(0.150000),
            Some(0.150000),
            Some(0.200000),
            Some(0.200000),
            Some(0.200000),
            Some(0.200000),
            Some(0.200000),
            Some(0.250000),
            Some(0.100000),
            Some(0.150000),
            Some(0.750000),
            Some(0.400000),
            Some(0.500000),
            Some(0.400000),
            Some(0.050000),
            Some(0.450000),
            Some(0.500000),
            Some(0.500000),
            Some(0.600000),
            Some(0.600000),
            Some(0.650000),
            Some(0.700000),
            Some(0.500000),
            Some(0.550000),
            Some(0.600000),
            Some(0.500000),
            Some(0.500000),
            Some(0.750000),
            Some(0.800000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.900000),
            Some(0.750000),
            Some(0.600000),
            Some(0.200000),
            Some(0.800000),
            Some(0.850000),
            Some(0.050000),
            Some(0.050000),
            Some(0.900000),
            Some(0.200000),
            Some(0.400000),
            Some(0.200000),
            Some(0.250000),
            Some(0.100000),
            Some(0.400000),
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
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
        "invasive_non-typhoidal_salmonella_spp.",
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
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.900000),
            Some(0.850000),
            Some(0.100000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.700000),
            Some(0.800000),
            Some(0.100000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
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
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.850000),
            Some(0.800000),
            Some(0.100000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
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
        "mycoplasma_genitalium",
        [
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
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.800000),
            Some(0.900000),
            Some(0.900000),
            Some(0.200000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.300000),
            Some(0.500000),
            Some(0.850000),
            Some(0.450000),
            Some(0.400000),
            Some(0.600000),
            Some(0.700000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.200000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.100000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
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
        "bacteroides_fragilis",
        [
            Some(0.050000),
            Some(0.100000),
            Some(0.200000),
            Some(0.250000),
            Some(0.500000),
            Some(0.400000),
            Some(0.050000),
            Some(0.050000),
            Some(0.200000),
            Some(0.200000),
            Some(0.250000),
            Some(0.250000),
            Some(0.200000),
            Some(0.950000),
            Some(0.950000),
            Some(0.950000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.600000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.250000),
            Some(0.350000),
            Some(0.500000),
            Some(0.250000),
            Some(0.300000),
            Some(0.500000),
            Some(0.500000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.300000),
            Some(0.700000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.950000),
            Some(0.050000),
            Some(0.200000),
            Some(0.750000),
            Some(0.850000),
            Some(0.750000),
            Some(0.800000),
            Some(0.500000),
            Some(0.950000),
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
        "mdr_mycobacterium_tuberculosis",
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
    (
        "mycoplasma_pneumoniae",
        [
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
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.800000),
            Some(0.850000),
            Some(0.800000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
            Some(0.600000),
            Some(0.700000),
            Some(0.750000),
            Some(0.800000),
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
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
        ],
    ),
    (
        "legionella_pneumophila",
        [
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
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.800000),
            Some(0.800000),
            Some(0.900000),
            Some(0.800000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.900000),
            Some(0.950000),
            Some(0.900000),
            Some(0.700000),
            Some(0.800000),
            Some(0.850000),
            Some(0.900000),
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
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
            Some(0.050000),
        ],
    ),
    (
        "burkholderia_cepacia_complex",
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
            Some(0.700000),
            Some(0.750000),
            Some(0.100000),
            Some(0.800000),
            Some(0.800000),
            Some(0.100000),
            Some(0.100000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.700000),
            Some(0.650000),
            Some(0.750000),
            Some(0.600000),
            Some(0.650000),
            Some(0.600000),
            Some(0.600000),
            Some(0.600000),
            Some(0.650000),
            Some(0.700000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.600000),
            Some(0.700000),
            Some(0.100000),
            Some(0.000000),
            Some(0.000000),
            Some(0.000000),
            Some(0.100000),
            Some(0.500000),
            Some(0.050000),
            Some(0.650000),
            Some(0.650000),
            Some(0.600000),
            Some(0.650000),
            Some(0.750000),
            Some(0.800000),
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

    eprintln!(
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
        map.insert(format!("{}_base_bacteria_level_change", bacteria), 0.5); // 0.2 // base change in bacteria level per day
        map.insert(format!("{}_max_level", bacteria), 5.0); // max bacteria level (arbitrary standardized scale)

            // --- Symptom Onset Parameters (Clinical Presentation) ---
        map.insert(format!("{}_daily_symptom_onset_probability", bacteria), 0.15); // Default: 15% chance per day of developing symptoms
        map.insert(format!("{}_symptom_onset_threshold_level", bacteria), 0.5); // Minimum bacteria level needed for symptom onset
        map.insert(format!("{}_symptom_onset_delay_days", bacteria), 1.0); // Minimum days infected before symptoms can start
        map.insert(format!("{}_symptom_onset_level_multiplier", bacteria), 1.0); // How much higher bacteria levels increase symptom probability
        map.insert(
                format!("{}_mechanismless_resistance_reversion_rate", bacteria),
                0.0004,
            ); // Daily probability of losing resistance when no specific mechanism is present
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
        // 
        // =====================================================================
        // COUNTERFACTUAL HGT CONTROL:
        // To run a "No Horizontal Gene Transfer" counterfactual scenario, simply modify the
        // `default_prob` value below to be 0.0 for all pairs, or in your Python launcher script
        // programmatically generate and pass all `hgt_prob_{donor}_to_{recipient}: 0.0` pairs.
        // Doing this will bypass all downstream HGT operations.
        // =====================================================================
        // --- HGT Probabilities for All Donor-Recipient Bacteria Pairs ---
        for (donor_idx, &donor) in BACTERIA_LIST.iter().enumerate() {
            for (recipient_idx, &recipient) in BACTERIA_LIST.iter().enumerate() {
                if donor_idx == recipient_idx {
                    continue;
                }

                let default_prob = default_hgt_probability(donor_idx, recipient_idx);
            map.insert(
                    format!("hgt_prob_{}_to_{}", donor, recipient),
                    default_prob,
                );
            }
        }


        // === [C] Drug initiation, selection, and pharmacokinetics ===
        // Core knobs for therapy behaviour: initiation heuristics, scoring multipliers, and half-lives
        // that drive drug levels. Overwrite these for global experiments; use per-drug keys for specifics.
        
        // *** Logistic model for antibiotic initiation probability ***
        // P(initiation) = 1 / (1 + exp(-log_odds))
        // log_odds = base + sum of applicable effects (additive in log-odds space)
        // This replaces the old multiplicative model and naturally bounds P ∈ (0,1)
        // Defaults calibrated to match historical behavior:
        //   - Base alone: P ≈ 0.1% (background rate)
        //   - With symptomatic infection: P ≈ 26% per day
        //   - With test + immunodeficiency: can approach high certainty
        map.insert("antibiotic_initiation_base_log_odds".to_string(), -5.5); // baseline: log(0.001/0.999) ≈ -6.9
        map.insert("antibiotic_initiation_log_odds_symptomatic_infection".to_string(), 6.0); // +5.85 → brings P from 0.1% to ~26%
        map.insert("antibiotic_initiation_log_odds_sepsis".to_string(), 6.0); // +6.0 -> strong boost for emergency care
        map.insert("antibiotic_initiation_log_odds_test_identified".to_string(), 0.92); // log(2.5) - lab confirmation boost
        map.insert("antibiotic_initiation_log_odds_already_on_drug".to_string(), 0.18); // log(1.2) - modest boost for layered therapy
        map.insert("antibiotic_initiation_log_odds_immunodeficiency".to_string(), 2.08); // log(8.0) - prophylaxis for immunocompromised
        map.insert("antibiotic_initiation_log_odds_no_indication".to_string(), -1.05); // log(0.35) - penalty when no infection/immunodeficiency

        // Drug activity parameters (still used for bacteria level effects)
        // a non-bacteria-specific parameter that determines how rapidly drugs of a given potency eliminate bacteria level
        // with a value 1 it is nearly always within 1 day
        map.insert(
            "drug_activity_to_bacteria_level_multiplier".to_string(),
            0.75,
        ); // Global scaling knob for drug-driven bacteria decay
        map.insert(
            "drug_activity_slow_clearance_probability".to_string(),
            0.25,
        ); // Share of infections with very slow pharmacodynamic response   
        map.insert(
            "drug_activity_slow_clearance_multiplier".to_string(),
            0.2,
        ); // Activity multiplier assigned to difficult-clearance infections
        map.insert("drug_decay_per_day".to_string(), 1.0); // Legacy parameter - now using drug-specific half-lives

        // Drug Selection Algorithm Parameters
        map.insert("drug_selection_temperature".to_string(), 0.55); // Lowered from 0.70 to concentrate prescribing on template drugs, reducing long-tail class leakage
        map.insert("reserve_drug_score_penalty".to_string(), 0.005); // Stronger restriction for reserve agents to achieve ~0.5-2% usage

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

        // Sulfonamides (first antibiotics)
        map.insert("drug_sulfanilamide_half_life_days".to_string(), 0.45); // ~11 hours

        // Beta-lactams (Penicillins)
        map.insert("drug_penicillin_g_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ampicillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_amoxicillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_piperacillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ticarcillin_half_life_days".to_string(), 0.046); // ~1.1 hours

        // Cephalosporins
        map.insert("drug_cephalexin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_cefazolin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_flucloxacillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_cefuroxime_half_life_days".to_string(), 0.05); // ~1.3 hours
        map.insert("drug_ceftriaxone_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_cefixime_half_life_days".to_string(), 0.125); // ~3 hours
        map.insert("drug_ceftazidime_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_cefepime_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_ceftaroline_half_life_days".to_string(), 0.11); // ~2.6 hours

        // Carbapenems
        map.insert("drug_meropenem_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_imipenem_c_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ertapenem_half_life_days".to_string(), 0.17); // ~4 hours

        // Monobactams
        map.insert("drug_aztreonam_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_aztreonam_avibactam_half_life_days".to_string(), 0.08); // ~2 hours

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
        map.insert("drug_chloramphenicol_half_life_days".to_string(), 0.125); // ~3 hours
        map.insert("drug_nitrofurantoin_half_life_days".to_string(), 0.017); // ~20 minutes
        map.insert("drug_retapamulin_half_life_days".to_string(), 0.25); // ~6 hours (topical, limited data)
        map.insert("drug_fusidic_a_half_life_days".to_string(), 0.375); // ~9 hours
        map.insert("drug_metronidazole_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_furazolidone_half_life_days".to_string(), 0.25); // ~6 hours

        // === NEW RESCUE DRUGS ===
        map.insert("drug_ceftolozane_tazobactam_half_life_days".to_string(), 0.125); // ~3 hours
        map.insert("drug_cefiderocol_half_life_days".to_string(), 0.1); // ~2.5 hours
        map.insert("drug_tigecycline_half_life_days".to_string(), 1.75); // ~42 hours
        map.insert("drug_daptomycin_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_fosfomycin_half_life_days".to_string(), 0.15); // ~4 hours
        map.insert("drug_fidaxomicin_half_life_days".to_string(), 0.5); // ~12 hours
        map.insert("drug_dalbavancin_half_life_days".to_string(), 10.0); // ~long active half-life


        // Colistin parameters (grouped with other drugs)
        map.insert("drug_colistin_spectrum_breadth".to_string(), 4.0); // Broad spectrum (mainly Gram-negative)
        // Regional availability (assume widely available, adjust as needed)
        map.insert("north_america_drug_colistin_availability".to_string(), 1.0);
        map.insert("europe_drug_colistin_availability".to_string(), 1.0);
        map.insert("asia_drug_colistin_availability".to_string(), 1.0);
        map.insert("oceania_drug_colistin_availability".to_string(), 1.0);
        map.insert("south_america_drug_colistin_availability".to_string(), 1.0);
        map.insert("africa_drug_colistin_availability".to_string(), 1.0);
        map.insert("home_drug_colistin_availability".to_string(), 1.0);


        // Toxicity hazard placeholders (per unit drug level). These represent best-guess daily fatal toxicity odds for active therapy.
        // Values scaled to produce realistic fatal adverse event rates when combined with age/immunodeficiency multipliers.
        // Reference: ~0.1-1% mortality rate for serious drug toxicity over typical 7-14 day course
        
        // === HIGH-TOXICITY DRUGS (nephrotoxicity, bone marrow suppression) ===
        map.insert("drug_colistin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000025); // Colistin-associated nephrotoxicity with high fatal risk
        map.insert("drug_gentamicin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000015); // Aminoglycoside renal failure/ototoxicity
        map.insert("drug_tobramycin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000013); // Similar aminoglycoside profile
        map.insert("drug_amikacin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000017); // Slightly higher renal toxicity than gentamicin
        map.insert("drug_vancomycin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000006); // Severe nephrotoxicity/red man syndrome rare but serious
        map.insert("drug_chloramphenicol_toxicity_death_hazard_per_unit_level".to_string(), 0.00000001); // Aplastic anemia risk (idiosyncratic, ~1:20,000-40,000)
        
        // === MODERATE-TOXICITY DRUGS (organ-specific toxicity) ===
        // Oxazolidinones - myelosuppression, lactic acidosis, peripheral neuropathy
        map.insert("drug_linezolid_toxicity_death_hazard_per_unit_level".to_string(), 0.000000008); // Thrombocytopenia, lactic acidosis with prolonged use
        map.insert("drug_tedizolid_toxicity_death_hazard_per_unit_level".to_string(), 0.000000004); // Lower toxicity than linezolid
        
        // Fluoroquinolones - tendon rupture, QT prolongation, CNS effects, aortic dissection
        map.insert("drug_ciprofloxacin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000003); // Tendinopathy, QT prolongation
        map.insert("drug_levofloxacin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000003); // Similar to ciprofloxacin
        map.insert("drug_moxifloxacin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000005); // Higher QT prolongation risk
        map.insert("drug_ofloxacin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000003); // Similar to other FQs
        
        // Other moderate-toxicity drugs
        map.insert("drug_rifampicin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000004); // Hepatotoxicity, drug interactions
        map.insert("drug_metronidazole_toxicity_death_hazard_per_unit_level".to_string(), 0.000000002); // Peripheral neuropathy with prolonged use
        map.insert("drug_nitrofurantoin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000003); // Pulmonary fibrosis (chronic use), hepatotoxicity
        map.insert("drug_trim_sulf_toxicity_death_hazard_per_unit_level".to_string(), 0.000000002); // Stevens-Johnson syndrome, bone marrow suppression
        
        // === LOW-TOXICITY DRUGS (rare serious adverse events) ===
        // Tetracyclines
        map.insert("drug_doxycycline_toxicity_death_hazard_per_unit_level".to_string(), 0.000000001); // Esophagitis, photosensitivity (rarely fatal)
        map.insert("drug_tetracycline_toxicity_death_hazard_per_unit_level".to_string(), 0.000000001);
        map.insert("drug_minocycline_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000015); // Slightly higher due to vestibular/autoimmune effects
        
        // Macrolides - QT prolongation, hepatotoxicity (rare)
        map.insert("drug_azithromycin_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000015); // QT prolongation, cardiac arrhythmia
        map.insert("drug_erythromycin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000002); // QT prolongation, hepatotoxicity
        map.insert("drug_clarithromycin_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000015); // Similar to azithromycin
        
        // Clindamycin - C. diff risk captured elsewhere; rare fatal anaphylaxis
        map.insert("drug_clindamycin_toxicity_death_hazard_per_unit_level".to_string(), 0.000000001);
        
        // === VERY LOW TOXICITY DRUGS (beta-lactams generally safe) ===
        // Beta-lactams - mainly anaphylaxis risk (~1:10,000-50,000 fatal reactions)
        map.insert("drug_penicillin_g_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000005); // Anaphylaxis
        map.insert("drug_ampicillin_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000005);
        map.insert("drug_amoxicillin_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000004);
        map.insert("drug_amoxicillin_clavulanate_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000006); // Slightly higher hepatotoxicity
        map.insert("drug_piperacillin_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000005);
        map.insert("drug_piperacillin_tazobactam_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000006);
        
        // Cephalosporins - similar to penicillins, cross-reactivity ~1-2%
        map.insert("drug_cephalexin_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000004);
        map.insert("drug_cefazolin_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000004);
        map.insert("drug_flucloxacillin_toxicity_death_hazard_per_unit_level".to_string(), 0.00000001);
        map.insert("drug_cefuroxime_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000004);
        map.insert("drug_ceftriaxone_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000005); // Biliary sludge in prolonged use
        map.insert("drug_cefixime_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000005);
        map.insert("drug_ceftazidime_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000004);
        map.insert("drug_cefepime_toxicity_death_hazard_per_unit_level".to_string(), 0.000000001); // Neurotoxicity in renal impairment
        map.insert("drug_ceftaroline_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000005);
        
        // Carbapenems - similar to beta-lactams, seizure risk with imipenem
        map.insert("drug_meropenem_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000006);
        map.insert("drug_imipenem_c_toxicity_death_hazard_per_unit_level".to_string(), 0.000000001); // Higher seizure risk
        map.insert("drug_ertapenem_toxicity_death_hazard_per_unit_level".to_string(), 0.0000000006);


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

        // === [D.2] Resistance floor configuration for rare bacteria ===
        // For bacteria with very low infection counts (like S. maltophilia and E. faecium at 100k pop),
        // the cache-based resistance sampling may not sustain observed resistance levels.
        // This feature provides minimum resistance floors that ramp up after drug introduction.
        // 
        // Drug introduction dates are already defined in DRUG_INTRODUCTION_DATES (lazy_static at bottom of file).
        // 
        // Enabled per-bacteria with: bacteria_{name}_resistance_floor_enabled = 1.0 (or 0.0 to disable)
        // Ramp period: bacteria_{name}_resistance_floor_ramp_years = years from drug intro to full floor
        // Per-drug-class floors: bacteria_{name}_{drug_class}_resistance_floor = target floor level (0.0-1.0)
        //
        // The floor is applied as: floor_level * ramp_fraction, where ramp_fraction = 
        // min(1.0, (current_day - drug_intro_day) / (ramp_years * 365)) for the earliest drug in the class.
        // If current_day < drug_intro_day, no floor is applied (resistance can't precede drug).

        // Master enable flag for resistance floors (set to 0.0 to disable globally)
        map.insert("resistance_floor_feature_enabled".to_string(), 1.0);

        // --- Stenotrophomonas maltophilia resistance floors ---
        // S. maltophilia has intrinsic L1/L2 beta-lactamases and multi-drug efflux pumps
        // At 100k population: ~552 infected person-days over 3 years, very sparse data
        // Target floors based on resistance_prevalence_values.csv intrinsic resistance patterns
        map.insert("bacteria_stenotrophomonas_maltophilia_resistance_floor_enabled".to_string(), 1.0);
        map.insert("bacteria_stenotrophomonas_maltophilia_resistance_floor_ramp_years".to_string(), 5.0); 
        // Drug class floors (using drug class names from potency section)
        map.insert("bacteria_stenotrophomonas_maltophilia_penicillins_resistance_floor".to_string(), 0.95); // Intrinsic L1/L2
        map.insert("bacteria_stenotrophomonas_maltophilia_cephalosporins_1_2_resistance_floor".to_string(), 0.95); // Intrinsic L1
        map.insert("bacteria_stenotrophomonas_maltophilia_cephalosporins_3_4_resistance_floor".to_string(), 0.75); // Partial L1/L2 coverage
        map.insert("bacteria_stenotrophomonas_maltophilia_carbapenems_resistance_floor".to_string(), 0.98); // Intrinsic L1 (metalloenzyme)
        map.insert("bacteria_stenotrophomonas_maltophilia_aminoglycosides_resistance_floor".to_string(), 0.80); // Efflux + modifying enzymes
        map.insert("bacteria_stenotrophomonas_maltophilia_fluoroquinolones_resistance_floor".to_string(), 0.45); // Moderate - acquired Smqnr
        map.insert("bacteria_stenotrophomonas_maltophilia_macrolides_resistance_floor".to_string(), 0.95); // Intrinsic efflux
        map.insert("bacteria_stenotrophomonas_maltophilia_tetracyclines_resistance_floor".to_string(), 0.40); // Variable - doxycycline/minocycline active
        map.insert("bacteria_stenotrophomonas_maltophilia_folate_antagonists_resistance_floor".to_string(), 0.15); // TMP-SMX is preferred therapy
        map.insert("bacteria_stenotrophomonas_maltophilia_polymyxins_resistance_floor".to_string(), 0.70); // Moderate colistin resistance

        // --- Enterococcus faecium resistance floors ---
        // E. faecium: intrinsically resistant to cephalosporins, low-level aminoglycosides, clindamycin
        // VRE (vancomycin-resistant) is a major concern globally
        // At 100k population: very low infection counts, resistance not sustained
        map.insert("bacteria_enterococcus_faecium_resistance_floor_enabled".to_string(), 0.0);
        map.insert("bacteria_enterococcus_faecium_resistance_floor_ramp_years".to_string(), 10.0); // Slower ramp - VRE emerged gradually
        // Drug class floors
        map.insert("bacteria_enterococcus_faecium_penicillins_resistance_floor".to_string(), 0.0); // Ampicillin resistance acquired, start at 0
        map.insert("bacteria_enterococcus_faecium_carbapenems_resistance_floor".to_string(), 0.0); // Not intrinsic
        map.insert("bacteria_enterococcus_faecium_aminoglycosides_resistance_floor".to_string(), 0.0); // High-level resistance acquired
        map.insert("bacteria_enterococcus_faecium_fluoroquinolones_resistance_floor".to_string(), 0.1); // High resistance observed
        map.insert("bacteria_enterococcus_faecium_macrolides_resistance_floor".to_string(), 0.1); // Moderate resistance
        map.insert("bacteria_enterococcus_faecium_glycopeptides_resistance_floor".to_string(), 0.05); // VRE - ~35-45% globally
        map.insert("bacteria_enterococcus_faecium_oxazolidinones_resistance_floor".to_string(), 0.01); // Low linezolid resistance
        map.insert("bacteria_enterococcus_faecium_tetracyclines_resistance_floor".to_string(), 0.05); // Moderate tetracycline resistance
        map.insert("bacteria_enterococcus_faecium_folate_antagonists_resistance_floor".to_string(), 0.20); // High TMP-SMX resistance

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
        // Potency values all defined above 

        // Define drug classes 
        // Polymyxins (currently only Colistin)
        let polymyxins = vec!["colistin"];
        let penicillins = vec!["penicillin_g", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin",
            // BL/BLI combinations
            "amoxicillin_clavulanate", "piperacillin_tazobactam", "ampicillin_sulbactam", "ticarcillin_clavulanate"
        ];
        let cephalosporins_1_2 = vec!["cephalexin", "cefazolin", "cefuroxime"];
        let _cephalosporins_3_4 = vec!["ceftriaxone", "cefixime", "ceftazidime", "cefepime", "ceftaroline"];
        let cephalosporins_3_4 = vec!["ceftriaxone", "ceftazidime", "cefixime", "cefepime", "ceftaroline", "ceftolozane_tazobactam", "cefiderocol",
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
        let tetracyclines = vec!["tetracycline", "doxycycline", "minocycline", "tigecycline"];
        let glycopeptides = vec!["vancomycin", "teicoplanin", "dalbavancin"];
        let oxazolidinones = vec!["linezolid", "tedizolid"];
        let _folate_antagonists = vec!["trim_sulf"];
        let _other_antibiotics = vec!["quinu_dalfo", "chloramphenicol", "nitrofurantoin", "retapamulin", "fusidic_a", "metronidazole", "furazolidone", "daptomycin", "fosfomycin", "fidaxomicin"];

        // Define bacterial groups for potency patterns - all names use underscores consistently
        let gram_pos_cocci = vec!["staphylococcus_aureus", "staphylococcus_epidermidis", "streptococcus_pneumoniae", "streptococcus_pyogenes", "streptococcus_agalactiae", "enterococcus_faecalis", "enterococcus_faecium"];
        let gram_neg_enterobacteria = vec!["escherichia_coli", "klebsiella_pneumoniae", "enterobacter_spp.", "citrobacter_spp.", "serratia_spp.", "proteus_spp.", "morganella_spp.", "enterobacter_cloacae", "p_stuartii"];
        let gram_neg_non_fermenting = vec!["pseudomonas_aeruginosa", "acinetobacter_baumannii", "stenotrophomonas_maltophilia", "burkholderia_cepacia_complex"];
        let fastidious_gram_neg = vec!["haemophilus_influenzae", "moraxella_catarrhalis", "neisseria_gonorrhoeae", "neisseria_meningitidis", "bordetella_pertussis", "legionella_pneumophila"];
        let enteric_pathogens = vec!["salmonella_enterica_serovar_typhi", "salmonella_enterica_serovar_paratyphi_a", "invasive_non-typhoidal_salmonella_spp.", "shigella_spp.", "vibrio_cholerae", "campylobacter_jejuni", "yersinia_enterocolitica"];
        let atypical_pathogens = vec!["chlamydia_trachomatis", "mycoplasma_genitalium", "mycoplasma_pneumoniae"];
        let anaerobes_spore_formers = vec!["clostridioides_difficile"];
        let obligate_anaerobes_non_spore = vec!["bacteroides_fragilis"];
        let gram_pos_rods = vec!["listeria_monocytogenes"];
        let gastric_pathogens = vec!["helicobacter_pylori"]; // Unique microaerophilic Gram-negative

        for &drug in DRUG_SHORT_NAMES.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                // Bacteria names now use underscores consistently
            map.insert(format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug, bacteria), 1.0);
            map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.1); // Default low potency 0.1
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
                            "penicillin_g" => 0.05,
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

        // Providencia stuartii is intrinsically resistant to polymyxins
        if BACTERIA_LIST.contains(&"p_stuartii") && DRUG_SHORT_NAMES.contains(&"colistin") {
        map.insert("drug_colistin_for_bacteria_p_stuartii_potency_when_no_r".to_string(), 0.0);
        }

        // klebsiella_pneumoniae - practical inactivity for most early beta-lactams
        let kleb_low_pen = vec![
            "penicillin_g", "ampicillin", "amoxicillin", "ampicillin_sulbactam",
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
        let gonorrhea_cephs = vec!["ceftriaxone", "cefixime", "ceftazidime", "cefepime", "ceftazidime_avibactam"];
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

        // OBLIGATE ANAEROBES (non-spore-forming Bacteroides)
        for &bacteria in obligate_anaerobes_non_spore.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Beta-lactams: BL/BLI and carbapenems are reliable; plain penicillins weak
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = match drug {
                            "piperacillin_tazobactam" | "ampicillin_sulbactam" | "ticarcillin_clavulanate" => 0.9,
                            "amoxicillin_clavulanate" => 0.7,
                            "piperacillin" => 0.75,
                            "ampicillin" | "amoxicillin" | "penicillin_g" => 0.0,
                            _ => 0.05,
                        };
                    map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }

                for &drug in cephalosporins_1_2.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                    map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.1);
                    }
                }
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                    map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.45);
                    }
                }

                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                    map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.95);
                    }
                }

                if DRUG_SHORT_NAMES.contains(&"metronidazole") {
                map.insert(format!("drug_metronidazole_for_bacteria_{}_potency_when_no_r", bacteria), 1.05);
                }
                if DRUG_SHORT_NAMES.contains(&"clindamycin") {
                map.insert(format!("drug_clindamycin_for_bacteria_{}_potency_when_no_r", bacteria), 0.55);
                }

                for &drug in aminoglycosides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                    map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.05);
                    }
                }
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                    map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.35);
                    }
                }
            }
        }

        // --- neisseria_meningitidis SPECIFIC POTENCY OVERRIDES ---
        // N. meningitidis is typically very penicillin-sensitive, unlike other fastidious gram-negatives
        // Clinical reality: Penicillin G and ampicillin are first-line therapies for sensitive strains
        if BACTERIA_LIST.contains(&"neisseria_meningitidis") {
            // Penicillins - EXCELLENT activity (first-line therapy)
        map.insert("drug_penicillin_g_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.95); // Excellent first-line
        map.insert("drug_ampicillin_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.90); // Excellent alternative

            // Ensure ceftriaxone remains excellent (current first-line)
        map.insert("drug_ceftriaxone_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.95); // Current first-line

            // Chloramphenicol - historically important alternative
            if DRUG_SHORT_NAMES.contains(&"chloramphenicol") {
            map.insert("drug_chloramphenicol_for_bacteria_neisseria_meningitidis_potency_when_no_r".to_string(), 0.85); // Good alternative
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

        if DRUG_SHORT_NAMES.contains(&"flucloxacillin") {
            for &bacteria in BACTERIA_LIST.iter() {
                let potency = match bacteria {
                    "staphylococcus_aureus" => 0.95,
                    "staphylococcus_epidermidis" => 0.85,
                    bacteria if bacteria.contains("streptococcus") => 0.80,
                    "enterococcus_faecalis" | "enterococcus_faecium" => 0.05,
                    "listeria_monocytogenes" => 0.05,
                    _ => 0.01,
                };
            map.insert(
                    format!(
                        "drug_flucloxacillin_for_bacteria_{}_potency_when_no_r",
                        bacteria
                    ),
                    potency,
                );
            }
        }

        if DRUG_SHORT_NAMES.contains(&"aztreonam_avibactam") {
            for &bacteria in BACTERIA_LIST.iter() {
                let potency = match bacteria {
                    "escherichia_coli" | "klebsiella_pneumoniae" | "enterobacter_spp."
                    | "citrobacter_spp." | "serratia_spp." | "proteus_spp."
                    | "morganella_spp." | "enterobacter_cloacae" | "p_stuartii" => 1.00,
                    "pseudomonas_aeruginosa" => 0.90,
                    "stenotrophomonas_maltophilia" => 0.75,
                    "burkholderia_cepacia_complex" => 0.60,
                    "acinetobacter_baumannii" => 0.10,
                    "haemophilus_influenzae" | "moraxella_catarrhalis"
                    | "neisseria_gonorrhoeae" | "neisseria_meningitidis"
                    | "bordetella_pertussis" => 0.80,
                    "legionella_pneumophila" | "chlamydia_trachomatis"
                    | "mycoplasma_genitalium" | "mycoplasma_pneumoniae"
                    | "clostridioides_difficile" | "bacteroides_fragilis"
                    | "listeria_monocytogenes" | "helicobacter_pylori"
                    | "staphylococcus_aureus" | "staphylococcus_epidermidis"
                    | "streptococcus_pneumoniae" | "streptococcus_pyogenes"
                    | "streptococcus_agalactiae" | "enterococcus_faecalis"
                    | "enterococcus_faecium" => 0.01,
                    _ => 0.90,
                };
            map.insert(
                    format!(
                        "drug_aztreonam_avibactam_for_bacteria_{}_potency_when_no_r",
                        bacteria
                    ),
                    potency,
                );
            }
        }

        // === [I] Novel Resistance Mechanism Overrides (Examples) ===
        // Use these to define emerging "Other" mechanisms specifically for certain bacteria.
        // This allows defining highly specific resistance threats (e.g., novel pumps or plasmids)
        // without enabling them globally for all bacteria properly.

        // Example 1: Novel Efflux Pump in Pseudomonas (affecting Meropenem)
        // map.insert("mechanism_enzyme_cat_applies_to_meropenem_in_pseudomonas_aeruginosa".to_string(), 1.0);
        // map.insert("mechanism_enzyme_cat_applies_to_imipenem_c_in_pseudomonas_aeruginosa".to_string(), 1.0);

        // Example 2: Novel Plasmid in E. coli (affecting Ciprofloxacin, reusing the SAME flag)
        // map.insert("mechanism_enzyme_cat_applies_to_ciprofloxacin_in_escherichia_coli".to_string(), 1.0);

        // Example 3: Global Override (affects ALL bacteria with this mechanism - use with caution)
        // map.insert("mechanism_enzyme_oxa_48_applies_to_colistin".to_string(), 1.0);





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
        map.insert("drug_penicillin_g_for_bacteria_vibrio_cholerae_potency_when_no_r".to_string(), 0.05);  // Poor activity
        map.insert("drug_ampicillin_for_bacteria_vibrio_cholerae_potency_when_no_r".to_string(), 0.05);  // Poor activity

        // haemophilus_influenzae - Address beta-lactamase resistance (intrinsic in many strains)
        // Reduce basic penicillins (H. flu commonly produces beta-lactamase)
        map.insert("drug_penicillin_g_for_bacteria_haemophilus_influenzae_potency_when_no_r".to_string(), 0.03); // Poor due to beta-lactamase
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
        map.insert("drug_penicillin_g_for_bacteria_bordetella_pertussis_potency_when_no_r".to_string(), 0.05);     // Poor activity
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
        map.insert("mycoplasma_genitalium_test_availability_year".to_string(), 1991.0); // PCR assays became available in the early 1990s

        // Bacteria-specific treatment recognition years (when bacteria was first recognized as needing treatment)
        map.insert("helicobacter_pylori_treatment_recognition_year".to_string(), 1982.0); // H. pylori not treated before Marshall & Warren discovery

        // Bacteria-specific sepsis risk overrides for organisms that don't cause acute sepsis are defined in the log-odds section below
        map.insert("helicobacter_pylori_base_bacteria_level_change".to_string(), 0.2); // Slow-growing chronic colonizer

        // H. pylori-specific drug selection bonuses when bacteria is identified
        map.insert("drug_clarithromycin_for_bacteria_helicobacter_pylori_initiation_multiplier".to_string(), 15.0); // Strong preference for triple therapy
        map.insert("drug_amoxicillin_for_bacteria_helicobacter_pylori_initiation_multiplier".to_string(), 12.0);   // Strong preference for triple therapy
        map.insert("drug_metronidazole_for_bacteria_helicobacter_pylori_initiation_multiplier".to_string(), 8.0);  // Alternative therapy
        map.insert("drug_tetracycline_for_bacteria_helicobacter_pylori_initiation_multiplier".to_string(), 6.0);   // Bismuth quadruple therapy
        map.insert("drug_levofloxacin_for_bacteria_helicobacter_pylori_initiation_multiplier".to_string(), 5.0);  // Rescue therapy

        // N. meningitidis-specific drug selection bonuses - emergency treatment for meningococcal disease
        map.insert("drug_penicillin_g_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 25.0); // First-line for sensitive strains
        map.insert("drug_ampicillin_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 22.0);  // Excellent alternative
        map.insert("drug_ceftriaxone_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 30.0); // Current standard of care
        map.insert("drug_cefotaxime_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 28.0);  // Equivalent 3rd generation
        map.insert("drug_chloramphenicol_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 18.0); // Important alternative historically
        map.insert("drug_ciprofloxacin_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 15.0); // Prophylaxis and treatment
        map.insert("drug_rifampicin_for_bacteria_neisseria_meningitidis_initiation_multiplier".to_string(), 12.0);  // Prophylaxis agent

        map.insert("drug_penicillin_g_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.05);       // Not used for H. pylori
        map.insert("drug_cephalexin_for_bacteria_helicobacter_pylori_potency_when_no_r".to_string(), 0.05);       // Not effective

        // --- BACTERIA-SPECIFIC SYMPTOM ONSET PARAMETERS (Logistic Model) ---
        // P(symptoms) = 1 / (1 + exp(-log_odds)), log_odds = base + level_effect
        // Converted from probability using: log_odds = ln(p / (1-p))

        // H. PYLORI - Usually asymptomatic chronic gastritis
        map.insert("helicobacter_pylori_symptom_onset_base_log_odds".to_string(), -6.9); // ~0.1% per day - very low symptomatic rate
        map.insert("helicobacter_pylori_symptom_onset_threshold_level".to_string(), 2.0);     // High threshold for symptoms
        map.insert("helicobacter_pylori_symptom_onset_delay_days".to_string(), 30.0);         // Long delay before symptoms possible

        // chlamydia_trachomatis - Often asymptomatic
        map.insert("chlamydia_trachomatis_symptom_onset_base_log_odds".to_string(), -4.6); // ~1% per day - often asymptomatic
        map.insert("chlamydia_trachomatis_symptom_onset_threshold_level".to_string(), 1.5);    // Moderate threshold
        map.insert("chlamydia_trachomatis_base_bacteria_level_change".to_string(), 0.3);       // Slow intracellular replication

        // neisseria_meningitidis - Often asymptomatic carriage
        map.insert("neisseria_meningitidis_symptom_onset_base_log_odds".to_string(), -1.1); // ~25% per day - increase clinical visibility
        map.insert("neisseria_meningitidis_base_bacteria_level_change".to_string(), 0.65);      // Fulminant meningococcemia progression
        map.insert("neisseria_meningitidis_symptom_onset_threshold_level".to_string(), 3.0);    // High threshold for invasive disease

        // moraxella_catarrhalis - Often just colonization
        map.insert("moraxella_catarrhalis_symptom_onset_base_log_odds".to_string(), -2.9);  // ~5% per day - often colonizer
        map.insert("moraxella_catarrhalis_symptom_onset_threshold_level".to_string(), 2.0);     // Moderate threshold
        map.insert("moraxella_catarrhalis_base_bacteria_level_change".to_string(), 0.55);       // Rapid otitis/sinusitis onset in children

        // bacteroides_fragilis - intra-abdominal abscesses; symptoms emerge when burden high
        map.insert("bacteroides_fragilis_symptom_onset_base_log_odds".to_string(), -0.2); // ~45% per day
        map.insert("bacteroides_fragilis_symptom_onset_threshold_level".to_string(), 1.2);
        map.insert("bacteroides_fragilis_symptom_onset_delay_days".to_string(), 2.0);
        map.insert("bacteroides_fragilis_base_bacteria_level_change".to_string(), 0.42);

        // p_stuartii - catheter-associated UTI/bacteremia; presents promptly when burdens rise
        map.insert("p_stuartii_symptom_onset_base_log_odds".to_string(), 0.2); // ~55% per day
        map.insert("p_stuartii_symptom_onset_threshold_level".to_string(), 0.75);
        map.insert("p_stuartii_base_bacteria_level_change".to_string(), 0.5);

        // mycoplasma_genitalium - frequently asymptomatic but persistent STI
        map.insert("mycoplasma_genitalium_symptom_onset_base_log_odds".to_string(), -2.0); // ~12% per day
        map.insert("mycoplasma_genitalium_symptom_onset_threshold_level".to_string(), 0.9);
        map.insert("mycoplasma_genitalium_symptom_onset_delay_days".to_string(), 5.0);
        map.insert("mycoplasma_genitalium_base_bacteria_level_change".to_string(), 0.28);

    // pseudomonas_aeruginosa - Clinically apparent when burden high
        map.insert("pseudomonas_aeruginosa_symptom_onset_base_log_odds".to_string(), -1.4);   // ~20% per day - improve detection of invasive disease
        map.insert("pseudomonas_aeruginosa_symptom_onset_threshold_level".to_string(), 0.8);     // Higher burden needed before symptoms manifest
        map.insert("pseudomonas_aeruginosa_base_bacteria_level_change".to_string(), 0.55);       // Rapid proliferation in ventilated hosts

        // ACUTE INFECTIONS - High symptomatic rates
        map.insert("haemophilus_influenzae_base_bacteria_level_change".to_string(), 0.55);       // Rapid pediatric respiratory progression
        map.insert("streptococcus_pneumoniae_symptom_onset_base_log_odds".to_string(), 1.4); // ~80% per day - usually symptomatic
        map.insert("streptococcus_pyogenes_symptom_onset_base_log_odds".to_string(), 0.85);   // ~70% per day - usually symptomatic
        map.insert("streptococcus_pyogenes_base_bacteria_level_change".to_string(), 0.6);         // Fast doubling in invasive GAS
        map.insert("staphylococcus_aureus_symptom_onset_base_log_odds".to_string(), 0.4);    // ~60% per day - usually symptomatic
        map.insert("staphylococcus_epidermidis_symptom_onset_base_log_odds".to_string(), -1.4); // ~20% per day - device-associated pathogen often subacute
        map.insert("staphylococcus_epidermidis_symptom_onset_threshold_level".to_string(), 1.0);   // Needs higher burden for symptoms
        map.insert("staphylococcus_epidermidis_symptom_onset_delay_days".to_string(), 3.0);        // Slight delay before clinical detection
        map.insert("staphylococcus_epidermidis_base_bacteria_level_change".to_string(), 0.35);     // Slower growth kinetics than S. aureus
        map.insert("staphylococcus_epidermidis_max_level".to_string(), 4.0);                        // Lower peak burden due to biofilm focus
        map.insert("staphylococcus_epidermidis_microbiome_clearance_probability_per_day".to_string(), 0.015); // Chronic colonizer of skin/devices
        map.insert("staphylococcus_epidermidis_log_odds_sepsis_infection_level".to_string(), 0.04); // Slight level effect on sepsis risk
        map.insert("staphylococcus_epidermidis_log_odds_sepsis_infection_duration".to_string(), 0.005); // Chronic devices slowly accumulate risk
        map.insert("staphylococcus_epidermidis_non_sepsis_infection_death_log_odds".to_string(), -6.0); // Very low direct mortality

        map.insert("stenotrophomonas_maltophilia_symptom_onset_base_log_odds".to_string(), -0.6); // ~35% per day - clinically apparent in ventilated hosts
        map.insert("stenotrophomonas_maltophilia_symptom_onset_threshold_level".to_string(), 0.9);      // Moderate burden before symptoms
        map.insert("stenotrophomonas_maltophilia_symptom_onset_delay_days".to_string(), 2.5);          // Early signs once established
        map.insert("stenotrophomonas_maltophilia_base_bacteria_level_change".to_string(), 0.45);       // Moderate growth rate
        map.insert("stenotrophomonas_maltophilia_max_level".to_string(), 5.0);                          // Can reach high burdens in lungs
        map.insert("stenotrophomonas_maltophilia_microbiome_clearance_probability_per_day".to_string(), 0.06); // Persistent colonizer in ICU settings
        map.insert("stenotrophomonas_maltophilia_log_odds_sepsis_infection_level".to_string(), 0.08);   // Rising burden increases risk notably
        map.insert("stenotrophomonas_maltophilia_log_odds_sepsis_infection_duration".to_string(), 0.012); // Prolonged infection raises odds
        map.insert("stenotrophomonas_maltophilia_non_sepsis_infection_death_log_odds".to_string(), -4.0); // Some mortality via pneumonia progression

        // Non-sepsis infection death log-odds overrides (bacteria where deaths come from non-sepsis mechanisms)
        // Formula: daily P(death) = logistic(base[-9.0] + bacteria_adj + syndrome_adj + level*coeff + age + ...)
        // Negative values reduce non-sepsis death; positive values increase it.
        // --- Reduce non-sepsis death (over-mortality bacteria) ---
        map.insert("chlamydia_trachomatis_non_sepsis_infection_death_log_odds".to_string(), -5.0);  // Sepsis at -19 is negligible; deaths are non-sepsis artefact (CFR 128× over)
        map.insert("mycoplasma_genitalium_non_sepsis_infection_death_log_odds".to_string(), -4.5);  // STI with essentially no sepsis; deaths are non-sepsis artefact (CFR 66× over)
        map.insert("neisseria_gonorrhoeae_non_sepsis_infection_death_log_odds".to_string(), -2.5);  // Gonorrhea rarely fatal in real life (CFR 11.6× over)
        map.insert("campylobacter_jejuni_non_sepsis_infection_death_log_odds".to_string(), -0.5);   // Sepsis at -20; non-sepsis is main death route (CFR 1.4× over)
        map.insert("mycoplasma_pneumoniae_non_sepsis_infection_death_log_odds".to_string(), -0.7);  // Low-mortality respiratory pathogen (CFR 1.9× over)
        // --- Increase non-sepsis death (under-mortality bacteria where sepsis is not the mechanism) ---
        map.insert("bordetella_pertussis_non_sepsis_infection_death_log_odds".to_string(), 4.0);    // Deaths from respiratory failure in infants, not sepsis (CFR 0.03×)
        map.insert("treponema_pallidum_non_sepsis_infection_death_log_odds".to_string(), 3.5);      // Tertiary/congenital syphilis deaths (CFR 0.06×)
        map.insert("vibrio_cholerae_non_sepsis_infection_death_log_odds".to_string(), 2.5);         // Death from dehydration, not sepsis (CFR 0.07×)
        map.insert("clostridioides_difficile_non_sepsis_infection_death_log_odds".to_string(), 2.0); // Colitis/toxic megacolon deaths (CFR 0.10×)
        map.insert("streptococcus_pyogenes_non_sepsis_infection_death_log_odds".to_string(), 1.5);  // RHD and post-streptococcal complications (CFR 0.20×)
        map.insert("bacteroides_fragilis_non_sepsis_infection_death_log_odds".to_string(), 1.5);    // Abscess mortality; sepsis default -14 with no override (CFR 0.27×)
        map.insert("helicobacter_pylori_non_sepsis_infection_death_log_odds".to_string(), 1.7);     // Gastric cancer deaths; sepsis at -250 (CFR 0.37×)
        map.insert("shigella_spp._non_sepsis_infection_death_log_odds".to_string(), 1.0);           // Dysentery deaths in children; sepsis at -12 contributes little (CFR 0.44×)

        // ENTERIC PATHOGENS - Moderate to high symptomatic rates
        map.insert("salmonella_enterica_serovar_typhi_symptom_onset_base_log_odds".to_string(), -0.4);       // ~40% per day
        map.insert("salmonella_enterica_serovar_typhi_base_bacteria_level_change".to_string(), 0.45);          // Longer incubation than typical enterics
        map.insert("salmonella_enterica_serovar_paratyphi_a_symptom_onset_base_log_odds".to_string(), -0.4); // ~40% per day
        map.insert("salmonella_enterica_serovar_paratyphi_a_base_bacteria_level_change".to_string(), 0.45);     // Similar incubation to typhi
        map.insert("shigella_spp._symptom_onset_base_log_odds".to_string(), 0.4);                           // ~60% per day
        map.insert("shigella_spp._base_bacteria_level_change".to_string(), 0.55);                               // Short incubation dysentery
        map.insert("vibrio_cholerae_symptom_onset_base_log_odds".to_string(), 0.0);                         // ~50% per day
        map.insert("vibrio_cholerae_base_bacteria_level_change".to_string(), 0.6);                              // Profuse cholera within 1-2 days
        map.insert("campylobacter_jejuni_symptom_onset_base_log_odds".to_string(), 0.0);                    // ~50% per day
        map.insert("campylobacter_jejuni_base_bacteria_level_change".to_string(), 0.52);                        // Incubation typically 2-4 days

        // CHRONIC/SLOW-ONSET PATHOGENS - Evidence-based symptom presentation rates
        // chlamydia_trachomatis - Most infections asymptomatic (~70-80% in women, ~50% in men)
        map.insert("chlamydia_trachomatis_symptom_onset_base_log_odds".to_string(), -3.5);  // Only ~3% daily → ~20-30% ever become symptomatic
        map.insert("chlamydia_trachomatis_base_bacteria_level_change".to_string(), 0.25);       // Slow intracellular replication
        map.insert("chlamydia_trachomatis_symptom_onset_threshold_level".to_string(), 0.8);     // Higher threshold before symptoms

        // treponema_pallidum - Syphilis has defined stages with variable presentation
        map.insert("treponema_pallidum_symptom_onset_base_log_odds".to_string(), -2.4);     // ~8% per day - Primary chancre develops in ~3-4 weeks
        map.insert("treponema_pallidum_base_bacteria_level_change".to_string(), 0.15);          // Very slow spirochete replication (33-hour doubling)
        map.insert("treponema_pallidum_symptom_onset_threshold_level".to_string(), 0.6);        // Moderate threshold

        // bordetella_pertussis - Catarrhal stage followed by paroxysmal cough
        map.insert("bordetella_pertussis_symptom_onset_base_log_odds".to_string(), -0.6);   // ~35% per day - 1-2 week incubation
        map.insert("bordetella_pertussis_base_bacteria_level_change".to_string(), 0.42);        // Moderate growth during catarrhal phase

        // helicobacter_pylori - Most infections (~80%) are asymptomatic
        map.insert("helicobacter_pylori_symptom_onset_base_log_odds".to_string(), -5.3);   // Only ~0.5% daily → ~20% develop symptomatic disease
        map.insert("helicobacter_pylori_symptom_onset_threshold_level".to_string(), 1.5);       // Very high threshold (chronic colonization)

        // MDR-TB - Slow progression; most latent infections never reactivate
        map.insert("mdr_mycobacterium_tuberculosis_symptom_onset_base_log_odds".to_string(), -6.9); // ~0.1% daily → ~5-10% lifetime reactivation risk
        map.insert("mdr_mycobacterium_tuberculosis_base_bacteria_level_change".to_string(), 0.08);       // Very slow mycobacterial growth
        map.insert("mdr_mycobacterium_tuberculosis_symptom_onset_threshold_level".to_string(), 2.0);     // High threshold for active disease

        // neisseria_gonorrhoeae - Variable symptoms (~10-20% asymptomatic in men, ~50% in women)
        map.insert("neisseria_gonorrhoeae_symptom_onset_base_log_odds".to_string(), -1.1);  // ~25% daily - Most symptomatic within 2-7 days
        map.insert("neisseria_gonorrhoeae_base_bacteria_level_change".to_string(), 0.55);       // Rapid mucosal colonization

        // --- COMPREHENSIVE BACTERIA GROWTH RATE OVERRIDES ---
        // Base level change per day reflects in vivo growth kinetics, NOT lab doubling times
        // Clinical progression depends on host factors, tissue site, and immune response
        // Default is 0.5/day; values below override specific pathogens based on microbiology
        
        // FULMINANT PATHOGENS (rapidly progressive, often life-threatening)
        map.insert("streptococcus_pyogenes_base_bacteria_level_change".to_string(), 0.7);       // Necrotizing fasciitis progresses in hours; invasive GAS very aggressive
        map.insert("neisseria_meningitidis_base_bacteria_level_change".to_string(), 0.65);      // Fulminant meningococcemia/purpura fulminans
        map.insert("vibrio_cholerae_base_bacteria_level_change".to_string(), 0.7);              // Massive fluid loss within 12-24 hours, extremely rapid toxin production
        map.insert("staphylococcus_aureus_base_bacteria_level_change".to_string(), 0.6);        // Rapid in endocarditis, bacteremia, necrotizing pneumonia
        map.insert("streptococcus_pneumoniae_base_bacteria_level_change".to_string(), 0.6);     // Rapid pneumonia/meningitis progression
        
        // RAPID PROGRESSORS (symptomatic within days)
        map.insert("acinetobacter_baumannii_base_bacteria_level_change".to_string(), 0.55);     // Rapid VAP/bacteremia in ICU patients
        map.insert("pseudomonas_aeruginosa_base_bacteria_level_change".to_string(), 0.55);      // Rapid proliferation in ventilated/immunocompromised hosts
        map.insert("haemophilus_influenzae_base_bacteria_level_change".to_string(), 0.55);      // Rapid otitis/meningitis in children
        map.insert("shigella_spp._base_bacteria_level_change".to_string(), 0.55);               // Dysentery within 1-3 days, low infectious dose
        map.insert("clostridioides_difficile_base_bacteria_level_change".to_string(), 0.55);    // Rapid toxin-mediated colitis after microbiome disruption
        map.insert("moraxella_catarrhalis_base_bacteria_level_change".to_string(), 0.55);       // Rapid otitis/sinusitis onset
        map.insert("neisseria_gonorrhoeae_base_bacteria_level_change".to_string(), 0.55);       // Rapid urethritis/cervicitis (2-7 days)
        map.insert("klebsiella_pneumoniae_base_bacteria_level_change".to_string(), 0.52);       // Rapid progression in pneumonia, can be necrotizing
        map.insert("legionella_pneumophila_base_bacteria_level_change".to_string(), 0.55);      // Acute onset pneumonia
        
        // MODERATE PROGRESSORS (typical acute infections)
        map.insert("escherichia_coli_base_bacteria_level_change".to_string(), 0.5);             // Variable by site; UTI to bacteremia
        map.insert("mycoplasma_pneumoniae_base_bacteria_level_change".to_string(), 0.35);       // Walking pneumonia (slower onset)
        map.insert("burkholderia_cepacia_complex_base_bacteria_level_change".to_string(), 0.45); // Moderate progression (CF exacerbation)
        map.insert("enterobacter_spp._base_bacteria_level_change".to_string(), 0.5);            // Nosocomial infections
        map.insert("enterobacter_cloacae_base_bacteria_level_change".to_string(), 0.5);         // Similar to other Enterobacter
        map.insert("campylobacter_jejuni_base_bacteria_level_change".to_string(), 0.52);        // Gastroenteritis 2-5 day incubation
        map.insert("enterococcus_faecalis_base_bacteria_level_change".to_string(), 0.48);       // Variable; endocarditis slow, UTI faster
        map.insert("enterococcus_faecium_base_bacteria_level_change".to_string(), 0.48);        // Similar to E. faecalis
        map.insert("citrobacter_spp._base_bacteria_level_change".to_string(), 0.5);             // Opportunistic, moderate progression
        map.insert("proteus_spp._base_bacteria_level_change".to_string(), 0.5);                 // UTI with moderate progression
        map.insert("serratia_spp._base_bacteria_level_change".to_string(), 0.48);               // Opportunistic, somewhat slower
        map.insert("morganella_spp._base_bacteria_level_change".to_string(), 0.48);             // Opportunistic, moderate
        map.insert("streptococcus_agalactiae_base_bacteria_level_change".to_string(), 0.52);    // Neonatal sepsis can be rapid
        map.insert("p_stuartii_base_bacteria_level_change".to_string(), 0.5);                   // Catheter-associated UTI
        map.insert("yersinia_enterocolitica_base_bacteria_level_change".to_string(), 0.45);     // 4-7 day incubation, pseudoappendicitis
        map.insert("stenotrophomonas_maltophilia_base_bacteria_level_change".to_string(), 0.45); // Opportunistic, moderate growth
        map.insert("salmonella_enterica_serovar_typhi_base_bacteria_level_change".to_string(), 0.45);      // Longer incubation (1-3 weeks)
        map.insert("salmonella_enterica_serovar_paratyphi_a_base_bacteria_level_change".to_string(), 0.45); // Similar to typhi
        map.insert("invasive_non-typhoidal_salmonella_spp._base_bacteria_level_change".to_string(), 0.5);  // More acute than typhi
        
        // SLOW PROGRESSORS (indolent or chronic infections)
        map.insert("bacteroides_fragilis_base_bacteria_level_change".to_string(), 0.42);        // Abscess formation is gradual
        map.insert("bordetella_pertussis_base_bacteria_level_change".to_string(), 0.42);        // 1-2 week catarrhal phase
        map.insert("staphylococcus_epidermidis_base_bacteria_level_change".to_string(), 0.35);  // Biofilm-associated, indolent
        map.insert("mycoplasma_genitalium_base_bacteria_level_change".to_string(), 0.28);       // Slow-growing, persistent
        map.insert("chlamydia_trachomatis_base_bacteria_level_change".to_string(), 0.25);       // Obligate intracellular, 48-72h cycle
        map.insert("listeria_monocytogenes_base_bacteria_level_change".to_string(), 0.25);      // Long incubation despite fast lab growth (intracellular)
        
        // VERY SLOW PROGRESSORS (chronic infections)
        map.insert("helicobacter_pylori_base_bacteria_level_change".to_string(), 0.2);          // Chronic colonizer, years to decades
        map.insert("treponema_pallidum_base_bacteria_level_change".to_string(), 0.18);          // 30+ hour doubling time, stages over weeks-months
        map.insert("mdr_mycobacterium_tuberculosis_base_bacteria_level_change".to_string(), 0.15); // 18-24 hour doubling, months to years progression

        // yersinia_enterocolitica - Address intrinsic penicillin resistance
        // Reduce penicillins (intrinsic resistance)
        map.insert("drug_penicillin_g_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.02); // Intrinsic resistance
        map.insert("drug_ampicillin_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.02);  // Intrinsic resistance
        // Boost appropriate drugs modestly
        map.insert("drug_doxycycline_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.75); // Good activity
        map.insert("drug_ciprofloxacin_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.7); // Good activity
        map.insert("drug_trim_sulf_for_bacteria_yersinia_enterocolitica_potency_when_no_r".to_string(), 0.65);    // Good activity

        // streptococcus_pyogenes - Ensure penicillin remains preferred (no resistance ever develops)
        // S. pyogenes has never developed penicillin resistance - boost slightly to counter any drift
        map.insert("drug_penicillin_g_for_bacteria_streptococcus_pyogenes_potency_when_no_r".to_string(), 0.95); // Excellent and consistent activity

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

        map.insert("drug_penicillin_g_for_bacteria_staphylococcus_epidermidis_potency_when_no_r".to_string(), 0.05); // Widespread beta-lactam resistance
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

        if BACTERIA_LIST.contains(&"mycoplasma_genitalium") {
        map.insert("drug_azithromycin_for_bacteria_mycoplasma_genitalium_initiation_multiplier".to_string(), 8.0); // Azithro first-line for M. genitalium
        map.insert("drug_moxifloxacin_for_bacteria_mycoplasma_genitalium_initiation_multiplier".to_string(), 4.0); // Fluoroquinolone rescue for macrolide failure
        map.insert("drug_levofloxacin_for_bacteria_mycoplasma_genitalium_initiation_multiplier".to_string(), 2.5); // Secondary FQ option when moxi unavailable
        map.insert("drug_doxycycline_for_bacteria_mycoplasma_genitalium_initiation_multiplier".to_string(), 1.5); // Doxycycline used for debulking per guidelines
        }

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
        map.insert("drug_penicillin_g_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_penicillin_g_for_bacteria_acinetobacter_baumannii_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_ampicillin_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_ampicillin_for_bacteria_acinetobacter_baumannii_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_amoxicillin_for_bacteria_pseudomonas_aeruginosa_initiation_multiplier".to_string(), 0.01);
        map.insert("drug_amoxicillin_for_bacteria_acinetobacter_baumannii_initiation_multiplier".to_string(), 0.01);

        // Beta-lactams irrelevant for Mycoplasma (No cell wall) and Legionella (Intracellular/Intrinsic)
        for &drug in &["penicillin_g", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "ceftriaxone", "meropenem", "ertapenem"] {
        map.insert(format!("drug_{}_for_bacteria_mycoplasma_pneumoniae_initiation_multiplier", drug), 0.001);
        map.insert(format!("drug_{}_for_bacteria_legionella_pneumophila_initiation_multiplier", drug), 0.001);
        }

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
        map.insert("drug_penicillin_g_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.05);  // No TB activity
        map.insert("drug_ampicillin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.05);   // No TB activity
        map.insert("drug_vancomycin_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.08);   // Minimal TB activity
        map.insert("drug_ceftriaxone_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.05);  // No TB activity
        map.insert("drug_meropenem_for_bacteria_mdr_mycobacterium_tuberculosis_potency_when_no_r".to_string(), 0.05);    // Minimal TB activity

        // RIFAMPICIN POTENCIES FOR OTHER BACTERIA (occasional use for severe staph infections)
        map.insert("drug_rifampicin_for_bacteria_staphylococcus_aureus_potency_when_no_r".to_string(), 0.4);         // Good anti-staph activity
        map.insert("drug_rifampicin_for_bacteria_enterococcus_faecalis_potency_when_no_r".to_string(), 0.2);         // Limited activity
        map.insert("drug_rifampicin_for_bacteria_enterococcus_faecium_potency_when_no_r".to_string(), 0.2);          // Limited activity
        // Most other bacteria: rifampicin has minimal activity (default 0.1 will apply)



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
        map.insert("mycoplasma_pneumoniae_drug_cessation_probability".to_string(), 0.015);  // Atypical pneumonia: 7-10 days

        // MODERATE-COURSE INFECTIONS (10-14 days typical treatment)
        map.insert("klebsiella_pneumoniae_drug_cessation_probability".to_string(), 0.0075); // Hospital pneumonia: 10-14 days
        map.insert("pseudomonas_aeruginosa_drug_cessation_probability".to_string(), 0.0075); // Complex infections: 14-21 days
        map.insert("acinetobacter_baumannii_drug_cessation_probability".to_string(), 0.0075); // Hospital infections: 14-21 days
        map.insert("burkholderia_cepacia_complex_drug_cessation_probability".to_string(), 0.0075); // Complex/CF: 14-21 days
        map.insert("legionella_pneumophila_drug_cessation_probability".to_string(), 0.0085); // Legionnaires: 10-14 days
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

        map.insert("neisseria_meningitidis_acquisition_log_odds_baseline".to_string(), -18.3);
        map.insert("haemophilus_influenzae_acquisition_log_odds_baseline".to_string(), -17.3);
        map.insert("salmonella_enterica_serovar_typhi_acquisition_log_odds_baseline".to_string(), -17.0);
        map.insert("bordetella_pertussis_acquisition_log_odds_baseline".to_string(), -13.3);
        map.insert("acinetobacter_baumannii_acquisition_log_odds_baseline".to_string(), -17.5);
        map.insert("campylobacter_jejuni_acquisition_log_odds_baseline".to_string(), -13.0);
        map.insert("chlamydia_trachomatis_acquisition_log_odds_baseline".to_string(), -13.0);
        map.insert("mycoplasma_genitalium_acquisition_log_odds_baseline".to_string(), -12.5);
        map.insert("mycoplasma_pneumoniae_acquisition_log_odds_baseline".to_string(), -12.0); // Periodic epidemics
        map.insert("legionella_pneumophila_acquisition_log_odds_baseline".to_string(), -15.8); // Sporadic environmental
        map.insert("burkholderia_cepacia_complex_acquisition_log_odds_baseline".to_string(), -17.5); // Rare/CF
        map.insert("citrobacter_spp._acquisition_log_odds_baseline".to_string(), -16.0);
        map.insert("clostridioides_difficile_acquisition_log_odds_baseline".to_string(), -15.2);
        map.insert("enterobacter_cloacae_acquisition_log_odds_baseline".to_string(), -15.5);
        map.insert("enterobacter_spp._acquisition_log_odds_baseline".to_string(), -16.0);
        map.insert("enterococcus_faecalis_acquisition_log_odds_baseline".to_string(), -16.5);
        map.insert("enterococcus_faecium_acquisition_log_odds_baseline".to_string(), -16.5);
        map.insert("escherichia_coli_acquisition_log_odds_baseline".to_string(), -11.5);
        map.insert("helicobacter_pylori_acquisition_log_odds_baseline".to_string(), -13.5);
        map.insert("invasive_non-typhoidal_salmonella_spp._acquisition_log_odds_baseline".to_string(), -17.4);
        map.insert("klebsiella_pneumoniae_acquisition_log_odds_baseline".to_string(), -15.3);
        map.insert("listeria_monocytogenes_acquisition_log_odds_baseline".to_string(), -19.6);
        map.insert("mdr_mycobacterium_tuberculosis_acquisition_log_odds_baseline".to_string(), -16.0);
        map.insert("moraxella_catarrhalis_acquisition_log_odds_baseline".to_string(), -15.0);
        map.insert("bacteroides_fragilis_acquisition_log_odds_baseline".to_string(), -15.3);
        map.insert("morganella_spp._acquisition_log_odds_baseline".to_string(), -15.5);
        map.insert("p_stuartii_acquisition_log_odds_baseline".to_string(), -16.1); 
        map.insert("neisseria_gonorrhoeae_acquisition_log_odds_baseline".to_string(), -13.5);
        map.insert("proteus_spp._acquisition_log_odds_baseline".to_string(), -15.7);
        map.insert("pseudomonas_aeruginosa_acquisition_log_odds_baseline".to_string(), -15.5);
        map.insert("salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds_baseline".to_string(), -16.5);
        map.insert("serratia_spp._acquisition_log_odds_baseline".to_string(), -16.5);
        map.insert("shigella_spp._acquisition_log_odds_baseline".to_string(), -12.8);
        map.insert("staphylococcus_epidermidis_acquisition_log_odds_baseline".to_string(), -15.5);
        map.insert("stenotrophomonas_maltophilia_acquisition_log_odds_baseline".to_string(), -16.8);
        map.insert("staphylococcus_aureus_acquisition_log_odds_baseline".to_string(), -12.6);
        map.insert("streptococcus_agalactiae_acquisition_log_odds_baseline".to_string(), -16.0);
        map.insert("streptococcus_pneumoniae_acquisition_log_odds_baseline".to_string(), -12.7);
        map.insert("streptococcus_pyogenes_acquisition_log_odds_baseline".to_string(), -14.8);
        map.insert("treponema_pallidum_acquisition_log_odds_baseline".to_string(), -13.2);
        map.insert("vibrio_cholerae_acquisition_log_odds_baseline".to_string(), -18.0);
        map.insert("yersinia_enterocolitica_acquisition_log_odds_baseline".to_string(), -16.0); // Slight boost in cases to get resistance signal
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
        map.insert("p_stuartii_log_odds_hospital_acquired".to_string(), 2.4); // Catheter-associated Providencia is predominantly nosocomial
        map.insert("bacteroides_fragilis_log_odds_hospital_acquired".to_string(), 1.2); // Post-surgical/abdominal procedures shift risk upward

        // Moderate HAI risk bacteria
        map.insert("escherichia_coli_log_odds_hospital_acquired".to_string(), 1.1); // ~3x higher risk (device-associated)
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
        map.insert("mycoplasma_genitalium_log_odds_hospital_acquired".to_string(), -1.5); // Strongly community-acquired STI; hospital cases rare

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
        map.insert("africa_escherichia_coli_acquisition_log_odds".to_string(), 1.2);
        map.insert("europe_escherichia_coli_acquisition_log_odds".to_string(), -0.3);
        map.insert("asia_escherichia_coli_acquisition_log_odds".to_string(), 1.0);
        map.insert("south_america_escherichia_coli_acquisition_log_odds".to_string(), 0.6);
        map.insert("oceania_escherichia_coli_acquisition_log_odds".to_string(), 0.1);

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
    map.insert("escherichia_coli_log_odds_microbiome_vs_infection".to_string(), 5.7);
    map.insert("enterococcus_faecalis_log_odds_microbiome_vs_infection".to_string(), 10.0);
    map.insert("enterococcus_faecium_log_odds_microbiome_vs_infection".to_string(), 11.0);
    map.insert("klebsiella_pneumoniae_log_odds_microbiome_vs_infection".to_string(), 9.0);
    map.insert("staphylococcus_aureus_log_odds_microbiome_vs_infection".to_string(), 8.5);
    map.insert("staphylococcus_epidermidis_log_odds_microbiome_vs_infection".to_string(), 11.3);
    map.insert("enterobacter_spp._log_odds_microbiome_vs_infection".to_string(), 10.0);
    map.insert("enterobacter_cloacae_log_odds_microbiome_vs_infection".to_string(), 11.5);
    map.insert("citrobacter_spp._log_odds_microbiome_vs_infection".to_string(), 10.0);
    map.insert("proteus_spp._log_odds_microbiome_vs_infection".to_string(), 6.5);
    map.insert("serratia_spp._log_odds_microbiome_vs_infection".to_string(), 8.8);
    map.insert("morganella_spp._log_odds_microbiome_vs_infection".to_string(), 8.7);
    map.insert("streptococcus_pneumoniae_log_odds_microbiome_vs_infection".to_string(), 7.0);
    map.insert("haemophilus_influenzae_log_odds_microbiome_vs_infection".to_string(), 11.2);
    map.insert("moraxella_catarrhalis_log_odds_microbiome_vs_infection".to_string(), 12.3);
    map.insert("bacteroides_fragilis_log_odds_microbiome_vs_infection".to_string(), 11.5);
    map.insert("streptococcus_pyogenes_log_odds_microbiome_vs_infection".to_string(), 8.8);
    map.insert("streptococcus_agalactiae_log_odds_microbiome_vs_infection".to_string(), 10.5);
    map.insert("acinetobacter_baumannii_log_odds_microbiome_vs_infection".to_string(), 7.0);
    map.insert("pseudomonas_aeruginosa_log_odds_microbiome_vs_infection".to_string(), 7.5);
    map.insert("clostridioides_difficile_log_odds_microbiome_vs_infection".to_string(), 7.0);
    map.insert("salmonella_enterica_serovar_typhi_log_odds_microbiome_vs_infection".to_string(), -4.0);
    map.insert("salmonella_enterica_serovar_paratyphi_a_log_odds_microbiome_vs_infection".to_string(), 1.1);
    map.insert("invasive_non-typhoidal_salmonella_spp._log_odds_microbiome_vs_infection".to_string(), 3.2);
    map.insert("shigella_spp._log_odds_microbiome_vs_infection".to_string(), 0.2);
    map.insert("vibrio_cholerae_log_odds_microbiome_vs_infection".to_string(), 2.5);
    map.insert("campylobacter_jejuni_log_odds_microbiome_vs_infection".to_string(), 2.5);
    map.insert("yersinia_enterocolitica_log_odds_microbiome_vs_infection".to_string(), 7.0);
    map.insert("listeria_monocytogenes_log_odds_microbiome_vs_infection".to_string(), 9.2);
    map.insert("p_stuartii_log_odds_microbiome_vs_infection".to_string(), 7.6);
    map.insert("neisseria_gonorrhoeae_log_odds_microbiome_vs_infection".to_string(), 1.0);
    map.insert("chlamydia_trachomatis_log_odds_microbiome_vs_infection".to_string(), 4.5);
    map.insert("mycoplasma_genitalium_log_odds_microbiome_vs_infection".to_string(), 3.5);
    map.insert("treponema_pallidum_log_odds_microbiome_vs_infection".to_string(), 8.0);
    map.insert("neisseria_meningitidis_log_odds_microbiome_vs_infection".to_string(), 9.8);
    map.insert("helicobacter_pylori_log_odds_microbiome_vs_infection".to_string(), 6.65);
    map.insert("mdr_mycobacterium_tuberculosis_log_odds_microbiome_vs_infection".to_string(), 1.0);
    map.insert("bordetella_pertussis_log_odds_microbiome_vs_infection".to_string(), 2.5);
    map.insert("stenotrophomonas_maltophilia_log_odds_microbiome_vs_infection".to_string(), 6.0);

    // Bacteria-specific microbiome clearance probabilities (per day)
    map.insert("escherichia_coli_microbiome_clearance_probability_per_day".to_string(), 0.005); // Persistent gut commensal; years-long colonization
    map.insert("enterococcus_faecalis_microbiome_clearance_probability_per_day".to_string(), 0.008); // Persistent gut flora; rarely cleared
    map.insert("enterococcus_faecium_microbiome_clearance_probability_per_day".to_string(), 0.06);
    map.insert("klebsiella_pneumoniae_microbiome_clearance_probability_per_day".to_string(), 0.03);
    map.insert("staphylococcus_aureus_microbiome_clearance_probability_per_day".to_string(), 0.05); // Nasal carriage persists weeks-months
    map.insert("enterobacter_spp._microbiome_clearance_probability_per_day".to_string(), 0.07);
    map.insert("enterobacter_cloacae_microbiome_clearance_probability_per_day".to_string(), 0.04);
    map.insert("citrobacter_spp._microbiome_clearance_probability_per_day".to_string(), 0.08);
    map.insert("proteus_spp._microbiome_clearance_probability_per_day".to_string(), 0.08);
    map.insert("serratia_spp._microbiome_clearance_probability_per_day".to_string(), 0.1);
    map.insert("morganella_spp._microbiome_clearance_probability_per_day".to_string(), 0.1);
    map.insert("p_stuartii_microbiome_clearance_probability_per_day".to_string(), 0.09); // Catheter biofilms persist but clear with device removal
    map.insert("bacteroides_fragilis_microbiome_clearance_probability_per_day".to_string(), 0.004); // Dominant anaerobe in gut microbiome; rarely displaced
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
    map.insert("mycoplasma_genitalium_microbiome_clearance_probability_per_day".to_string(), 0.18); // Persists for months without therapy
    map.insert("treponema_pallidum_microbiome_clearance_probability_per_day".to_string(), 0.35);
    map.insert("neisseria_meningitidis_microbiome_clearance_probability_per_day".to_string(), 0.05);
    map.insert("helicobacter_pylori_microbiome_clearance_probability_per_day".to_string(), 0.001); // Extremely persistent; decades without treatment
    map.insert("mdr_mycobacterium_tuberculosis_microbiome_clearance_probability_per_day".to_string(), 0.0015); // Latent carriage now clears slowly over years
    map.insert("bordetella_pertussis_microbiome_clearance_probability_per_day".to_string(), 0.2);
    map.insert("staphylococcus_epidermidis_clearance_probability_per_day".to_string(), 0.2);
    map.insert("stenotrophomonas_maltophilia_clearance_probability_per_day".to_string(), 0.2);


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
        // (Resistance reversion is handled by per-mechanism reversion_rate and per-bacteria
        // mechanismless_resistance_reversion_rate — see their respective config sections)
        // Resistance emergence is now entirely mechanism-based — see bacteria-mechanism emergence rates below.

        map.insert("resistance_emergence_bacteria_level_multiplier".to_string(), 9.0); // Multiplier for bacteria level's effect on emergence (ranges 1.0x to 10.0x)

        map.insert("any_r_emergence_level_on_first_emergence".to_string(), 0.5); // The resistance level 'any_r' starts at upon emergence


        //  Microbiome Resistance Transfer Parameter
        // ^^^ micro
        map.insert("microbiome_resistance_transfer_probability_per_day".to_string(), 0.0001); // 0.00001  *** Probability per day for resistance transfer between infection and microbiome

        // --- Multi-Drug Resistance Emergence Penalty Parameters ---
        // When multiple drugs are active, resistance emergence is reduced because mutations
        // must confer resistance to ALL active drugs to provide survival advantage
        map.insert("multi_drug_penalty_for_single_drug_resistance".to_string(), 0.05); // Penalty when resistance affects only 1 of multiple active drugs (5% survival)
        map.insert("multi_drug_penalty_for_partial_cross_resistance".to_string(), 0.3); // Penalty when resistance affects some but not all active drugs (30% survival)
        map.insert("multi_drug_penalty_threshold_num_drugs".to_string(), 2.0); // Minimum number of active drugs to trigger multi-drug penalty



        // --- Resistance Mechanisms Parameters ---
        // Implementation of granular resistance mechanisms (40 types)
        //
        // =====================================================================================
        // BACTERIA-MECHANISM-SPECIFIC EMERGENCE RATES
        // =====================================================================================
        // Direct emergence rates (per day when drug present) for each bacteria-mechanism pair.
        // All 40 resistance mechanisms in standardized biological order for all 42 bacteria.
        //
        // ^^^^
        //
        // ======================================================================
        // Gram-Negative Bacteria — Enterobacterales
        // ======================================================================
        // E. coli — Gram-negative, Enterobacterales
        // Band 6 (x0.6)
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_000_01    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_000_01    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_000_01  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_000_01   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_000_01   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_000_01     ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_000_01     ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_000_01     ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_1            ); // classes: chl
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_000_003    ); // classes: ag
        map.insert("bacteria_escherichia_coli_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_1              ); // classes: fq
        map.insert("bacteria_escherichia_coli_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_000_03  ); // classes: fq
        map.insert("bacteria_escherichia_coli_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_1              ); // classes: fq
        map.insert("bacteria_escherichia_coli_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_000_1    ); // classes: fq, tet, chl
        map.insert("bacteria_escherichia_coli_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_1     ); // classes: fq, mls, tet, chl
        map.insert("bacteria_escherichia_coli_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_000_1     ); // classes: none (currently zeroed)
        map.insert("bacteria_escherichia_coli_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_000_1     ); // classes: poly
        map.insert("bacteria_escherichia_coli_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_001      ); // classes: sulf
        map.insert("bacteria_escherichia_coli_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_000_001      ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_000_01        ); // classes: other (fosfomycin)
        map.insert("bacteria_escherichia_coli_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_1        ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_escherichia_coli_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_000_01    ); // classes: tet
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_01         ); // classes: ag
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; 
        map.insert("bacteria_escherichia_coli_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_000_01   ); // classes: tet
        map.insert("bacteria_escherichia_coli_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_escherichia_coli_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 


         // K. pneumoniae — Gram-negative, Enterobacterales
        // Band 7 (x10)
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_005   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_005    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_005   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_005    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_005    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_005     ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_005        ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_005        ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_cat_emergence_rate".to_string(), 0.002   ); // classes: chl
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.015  ); // classes: ag
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_5    ); // classes: fq
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_5   ); // classes: fq
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_5   ); // classes: fq
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.001   ); // classes: fq, tet, chl
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_3       ); // classes: fq, mls, tet, chl
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.000_001    ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_003 ); // classes: none (currently zeroed)
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.01   ); // classes: poly
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.001_5   ); // classes: sulf
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.002     ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.003    ); // classes: other (fosfomycin)
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_2  ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_2  ); // classes: tet
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.04      ); // classes: ag
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_3  ); // classes: tet
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_klebsiella_pneumoniae_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0


        // Citrobacter spp. — Gram-negative, Enterobacterales
        // Band 8 (x187.5)
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.002    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.002    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.002     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.002   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.002    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_kpc_emergence_rate".to_string(), 0.002    ) ; // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.002    ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.002    ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_cat_emergence_rate".to_string(), 0.05   ); // classes: chl
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.2   ); // classes: ag
        map.insert("bacteria_citrobacter_spp._mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.05   ); // classes: fq
        map.insert("bacteria_citrobacter_spp._mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_3  ); // classes: fq
        map.insert("bacteria_citrobacter_spp._mechanism_protection_qnr_emergence_rate".to_string(), 0.05   ); // classes: fq
        map.insert("bacteria_citrobacter_spp._mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.015  ); // classes: fq, tet, chl
        map.insert("bacteria_citrobacter_spp._mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.002    ); // classes: ceph, carb, fq, ag
        map.insert("bacteria_citrobacter_spp._mechanism_global_efflux_pump_emergence_rate".to_string(), 0.003      ); // classes: fq, mls, tet, chl
        map.insert("bacteria_citrobacter_spp._mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_5 ); // classes: none (currently zeroed)
        map.insert("bacteria_citrobacter_spp._mechanism_modification_mcr_1_emergence_rate".to_string(), 0.05    ); // classes: poly
        map.insert("bacteria_citrobacter_spp._mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.003   ); // classes: sulf
        map.insert("bacteria_citrobacter_spp._mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.015   ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.005    ); // classes: other (fosfomycin)
        map.insert("bacteria_citrobacter_spp._mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.02   ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_citrobacter_spp._mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_protection_tet_m_emergence_rate".to_string(), 0.005    ); // classes: tet
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.15  ); // classes: ag
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.005     ); // classes: tet
        map.insert("bacteria_citrobacter_spp._mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.002  ); // classes: pen, bli, ceph, mono; modest beta-lactam lift without rebuilding broad efflux resistance
        map.insert("bacteria_citrobacter_spp._mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_citrobacter_spp._mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

  // Enterobacter spp. — Gram-negative, Enterobacterales
        // Band 8 (x125)
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_5   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_5   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_5 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_5); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_5   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_5  ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_5    ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_5 ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_cat_emergence_rate".to_string(), 0.005 ); // classes: chl
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.3   ); // classes: ag
        map.insert("bacteria_enterobacter_spp._mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_5  ); // classes: fq
        map.insert("bacteria_enterobacter_spp._mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_5   ); // classes: fq
        map.insert("bacteria_enterobacter_spp._mechanism_protection_qnr_emergence_rate".to_string(), 0.000_5   ); // classes: fq
        map.insert("bacteria_enterobacter_spp._mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_5 ); // classes: fq, tet, chl
        map.insert("bacteria_enterobacter_spp._mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_5     ); // classes: fq, mls, tet, chl
        map.insert("bacteria_enterobacter_spp._mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_5   ); // classes: none (currently zeroed)
        map.insert("bacteria_enterobacter_spp._mechanism_modification_mcr_1_emergence_rate".to_string(), 0.03     ); // classes: poly
        map.insert("bacteria_enterobacter_spp._mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.001  ); // classes: sulf
        map.insert("bacteria_enterobacter_spp._mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.005  ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.02    ); // classes: other (fosfomycin)
        map.insert("bacteria_enterobacter_spp._mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.03    ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_enterobacter_spp._mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_protection_tet_m_emergence_rate".to_string(), 0.01   ); // classes: tet
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.3 ); // classes: ag
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.01  ); // classes: tet
        map.insert("bacteria_enterobacter_spp._mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_5); // classes: pen, bli, ceph, mono; modest beta-lactam lift without rebuilding broad non-beta-lactam resistance
        map.insert("bacteria_enterobacter_spp._mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_spp._mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0


        // E. cloacae — Gram-negative, Enterobacterales
        // Band 8 (x150)
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_01); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_01); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_01); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_01 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_01 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_01     ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_01  ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_01  ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_2   ); // classes: chl
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.03     ); // classes: ag
        map.insert("bacteria_enterobacter_cloacae_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_2 ); // classes: fq
        map.insert("bacteria_enterobacter_cloacae_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_2  ); // classes: fq
        map.insert("bacteria_enterobacter_cloacae_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_2  ); // classes: fq
        map.insert("bacteria_enterobacter_cloacae_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_05 ); // classes: fq, tet, chl
        map.insert("bacteria_enterobacter_cloacae_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.0    ); // classes: fq, mls, tet, chl
        map.insert("bacteria_enterobacter_cloacae_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_001); // classes: none (currently zeroed)
        map.insert("bacteria_enterobacter_cloacae_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.01  ); // classes: poly
        map.insert("bacteria_enterobacter_cloacae_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_05  ); // classes: sulf
        map.insert("bacteria_enterobacter_cloacae_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.001  ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.001    ); // classes: other (fosfomycin)
        map.insert("bacteria_enterobacter_cloacae_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.01   ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_enterobacter_cloacae_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_protection_tet_m_emergence_rate".to_string(), 0.002  ); // classes: tet
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.1     ); // classes: ag
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_1   ); // classes: tet
        map.insert("bacteria_enterobacter_cloacae_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_000_1); // classes: pen, bli, ceph, mono; modest beta-lactam lift without overdriving other classes
        map.insert("bacteria_enterobacter_cloacae_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterobacter_cloacae_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

       // Morganella spp. — Gram-negative, Enterobacterales
        // Band 8 (x300)
        map.insert("bacteria_morganella_spp._mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_3 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_morganella_spp._mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_3 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_morganella_spp._mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_3 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_morganella_spp._mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_1  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_morganella_spp._mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_1  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_morganella_spp._mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_1 ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_morganella_spp._mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_1  ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_morganella_spp._mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_1  ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_morganella_spp._mechanism_enzyme_cat_emergence_rate".to_string(), 0.1    ); // classes: chl
        map.insert("bacteria_morganella_spp._mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 1.0   ); // classes: ag
        map.insert("bacteria_morganella_spp._mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.02 ); // classes: fq
        map.insert("bacteria_morganella_spp._mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.02   ); // classes: fq
        map.insert("bacteria_morganella_spp._mechanism_protection_qnr_emergence_rate".to_string(), 0.02    ); // classes: fq
        map.insert("bacteria_morganella_spp._mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.002   ); // classes: fq, tet, chl
        map.insert("bacteria_morganella_spp._mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.000_5 ); // classes: ceph, carb, fq, ag
        map.insert("bacteria_morganella_spp._mechanism_global_efflux_pump_emergence_rate".to_string(), 0.003  ); // classes: fq, mls, tet, chl
        map.insert("bacteria_morganella_spp._mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_005 ); // classes: none (currently zeroed)
        map.insert("bacteria_morganella_spp._mechanism_modification_mcr_1_emergence_rate".to_string(), 0.02     ); // classes: poly
        map.insert("bacteria_morganella_spp._mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.005 ); // classes: sulf
        map.insert("bacteria_morganella_spp._mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.03    ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_morganella_spp._mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.003   ); // classes: other (fosfomycin)
        map.insert("bacteria_morganella_spp._mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.01   ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_morganella_spp._mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_5  ); // classes: tet
        map.insert("bacteria_morganella_spp._mechanism_enzyme_aac_aph_emergence_rate".to_string(), 1.0    ); // classes: ag
        map.insert("bacteria_morganella_spp._mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.001_5   ); // classes: tet
        map.insert("bacteria_morganella_spp._mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_000_5  ); // classes: pen, bli, ceph, mono; tiny seed; not treated as impossible
        map.insert("bacteria_morganella_spp._mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_morganella_spp._mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

        // Proteus spp. — Gram-negative, Enterobacterales
        // Band 8 (x50)
        map.insert("bacteria_proteus_spp._mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_1    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_proteus_spp._mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_1     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_proteus_spp._mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_03    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_proteus_spp._mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_03    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_proteus_spp._mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_03     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_proteus_spp._mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_02      ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_proteus_spp._mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_02      ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_proteus_spp._mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_02      ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_proteus_spp._mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_35    ); // classes: chl
        map.insert("bacteria_proteus_spp._mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.002         ); // classes: ag
        map.insert("bacteria_proteus_spp._mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.04       ); // classes: fq
        map.insert("bacteria_proteus_spp._mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.01        ); // classes: fq
        map.insert("bacteria_proteus_spp._mechanism_protection_qnr_emergence_rate".to_string(), 0.01        ); // classes: fq
        map.insert("bacteria_proteus_spp._mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.005      ); // classes: fq, tet, chl
        map.insert("bacteria_proteus_spp._mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_global_efflux_pump_emergence_rate".to_string(), 0.05     ); // classes: fq, mls, tet, chl
        map.insert("bacteria_proteus_spp._mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_15     ); // classes: none (currently zeroed)
        map.insert("bacteria_proteus_spp._mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_5      ); // classes: poly
        map.insert("bacteria_proteus_spp._mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.004       ); // classes: sulf
        map.insert("bacteria_proteus_spp._mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_015      ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_proteus_spp._mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0015       ); // classes: other (fosfomycin)
        map.insert("bacteria_proteus_spp._mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_1       ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_proteus_spp._mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_protection_tet_m_emergence_rate".to_string(), 0.003_5    ); // classes: tet
        map.insert("bacteria_proteus_spp._mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_001    ); // classes: ag
        map.insert("bacteria_proteus_spp._mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_001     ); // classes: tet
        map.insert("bacteria_proteus_spp._mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_proteus_spp._mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // Serratia spp. — Gram-negative, Enterobacterales
        // Band 8 (x250)
        map.insert("bacteria_serratia_spp._mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.001  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_serratia_spp._mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.001  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_serratia_spp._mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.001     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_serratia_spp._mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.001     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_serratia_spp._mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.001     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_serratia_spp._mechanism_enzyme_kpc_emergence_rate".to_string(), 0.001      ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_serratia_spp._mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_1    ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_serratia_spp._mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_1  ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_serratia_spp._mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_05 ); // classes: chl
        map.insert("bacteria_serratia_spp._mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.9     ); // classes: ag
        map.insert("bacteria_serratia_spp._mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.02    ); // classes: fq
        map.insert("bacteria_serratia_spp._mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.02    ); // classes: fq
        map.insert("bacteria_serratia_spp._mechanism_protection_qnr_emergence_rate".to_string(), 0.02      ); // classes: fq
        map.insert("bacteria_serratia_spp._mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.015    ); // classes: fq, tet, chl
        map.insert("bacteria_serratia_spp._mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_global_efflux_pump_emergence_rate".to_string(), 0.003   ); // classes: fq, mls, tet, chl
        map.insert("bacteria_serratia_spp._mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_03); // classes: none (currently zeroed)
        map.insert("bacteria_serratia_spp._mechanism_modification_mcr_1_emergence_rate".to_string(), 0.03     ); // classes: poly
        map.insert("bacteria_serratia_spp._mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.001  ); // classes: sulf
        map.insert("bacteria_serratia_spp._mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.002     ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_serratia_spp._mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_2    ); // classes: other (fosfomycin)
        map.insert("bacteria_serratia_spp._mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.03  ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_serratia_spp._mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_protection_tet_m_emergence_rate".to_string(), 0.002  ); // classes: tet
        map.insert("bacteria_serratia_spp._mechanism_enzyme_aac_aph_emergence_rate".to_string(), 1.0 ); // classes: ag
        map.insert("bacteria_serratia_spp._mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_3 ); // classes: tet
        map.insert("bacteria_serratia_spp._mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_serratia_spp._mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // P. stuartii — Gram-negative, Enterobacterales
        // Band 9 (x375)
        map.insert("bacteria_p_stuartii_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_19); // classes: pen, bli, ceph, mono
        map.insert("bacteria_p_stuartii_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_19); // classes: pen, bli, ceph, mono
        map.insert("bacteria_p_stuartii_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_037); // classes: pen, bli, ceph, mono
        map.insert("bacteria_p_stuartii_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_38); // classes: pen, bli, ceph, mono
        map.insert("bacteria_p_stuartii_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_19); // classes: pen, bli, ceph, mono
        map.insert("bacteria_p_stuartii_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_019); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_p_stuartii_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_019); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_p_stuartii_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_019); // classes: pen, bli, ceph, carb
        map.insert("bacteria_p_stuartii_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_19); // classes: chl
        map.insert("bacteria_p_stuartii_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_003_8); // classes: ag
        map.insert("bacteria_p_stuartii_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_38); // classes: fq
        map.insert("bacteria_p_stuartii_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_19); // classes: fq
        map.insert("bacteria_p_stuartii_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_19); // classes: fq
        map.insert("bacteria_p_stuartii_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_19); // classes: fq, tet, chl
        map.insert("bacteria_p_stuartii_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_037); // classes: fq, mls, tet, chl
        map.insert("bacteria_p_stuartii_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_037); // classes: none (currently zeroed)
        map.insert("bacteria_p_stuartii_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_019); // classes: poly
        map.insert("bacteria_p_stuartii_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_38); // classes: sulf
        map.insert("bacteria_p_stuartii_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_037); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_p_stuartii_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_037); // classes: other (fosfomycin)
        map.insert("bacteria_p_stuartii_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_003_8); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_p_stuartii_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_38); // classes: tet
        map.insert("bacteria_p_stuartii_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_000_000_05); // classes: ag
        map.insert("bacteria_p_stuartii_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_000_000_05); // classes: tet
        map.insert("bacteria_p_stuartii_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_p_stuartii_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // S. Typhi — Gram-negative, Enterobacterales
        // Band 7 (x7.9)
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.005    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.15     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.001_5     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_5      ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_5      ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_03       ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_03       ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_03       ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_cat_emergence_rate".to_string(), 0.015      ); // classes: chl
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_1        ); // classes: ag
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.8     ); // classes: fq
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.05       ); // classes: fq
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_protection_qnr_emergence_rate".to_string(), 0.03       ); // classes: fq
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.05     ); // classes: fq, tet, chl
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.005       ); // classes: fq, mls, tet, chl
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_global_porin_loss_emergence_rate".to_string(), 0.001       ); // classes: none (currently zeroed)
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.001_5     ); // classes: poly
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.15    ); // classes: sulf
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.001_5     ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.001_5     ); // classes: other (fosfomycin)
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_15     ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_protection_tet_m_emergence_rate".to_string(), 0.015      ); // classes: tet
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_1     ); // classes: ag
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_002   ); // classes: tet
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.01); // classes: pen, mls, tet, chl; low broad-efflux proxy
        map.insert("bacteria_salmonella_enterica_serovar_typhi_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // S. Paratyphi A — Gram-negative, Enterobacterales
        // Band 7 (x25)
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.002    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.002    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.002    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.001_5  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.001_5  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_05   ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_07 ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_05  ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_cat_emergence_rate".to_string(), 0.02       ); // classes: chl
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.15     ); // classes: ag
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.05     ); // classes: fq
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.04     ); // classes: fq
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_protection_qnr_emergence_rate".to_string(), 0.005    ); // classes: fq
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.03     ); // classes: fq, tet, chl
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.000_8); // classes: ceph, carb, fq, ag
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.004_5  ); // classes: fq, mls, tet, chl
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0     ); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_1  ); // classes: none (currently zeroed)
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.003    ); // classes: poly
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.03     ); // classes: sulf
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_15      ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_15      ); // classes: other (fosfomycin)
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.01     ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_protection_tet_m_emergence_rate".to_string(), 0.02     ); // classes: tet
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.15     ); // classes: ag
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.007  ); // classes: tet
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.001    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_5  ); // classes: pen, mls, tet, chl; low broad-efflux proxy
        map.insert("bacteria_salmonella_enterica_serovar_paratyphi_a_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // iNTS — Gram-negative, Enterobacterales
        // Band 8 (x37.5)
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_25   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_25   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_25   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_25  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_25  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_025 ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_025  ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_025  ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_cat_emergence_rate".to_string(), 0.015  ); // classes: chl
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.8     ); // classes: ag
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.01    ); // classes: fq
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.01     ); // classes: fq
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_protection_qnr_emergence_rate".to_string(), 0.02     ); // classes: fq
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.02   ); // classes: fq, tet, chl
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_global_efflux_pump_emergence_rate".to_string(), 0.02       ); // classes: fq, mls, tet, chl
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_01      ); // classes: none (currently zeroed)
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_modification_mcr_1_emergence_rate".to_string(), 0.03    ); // classes: poly
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.003    ); // classes: sulf
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.003        ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.001       ); // classes: other (fosfomycin)
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.02    ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_protection_tet_m_emergence_rate".to_string(), 0.004  ); // classes: tet
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.8       ); // classes: ag
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.08    ); // classes: tet
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_002_5  ); // classes: pen, mls, tet, chl; low broad-efflux proxy
        map.insert("bacteria_invasive_non-typhoidal_salmonella_spp._mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // Shigella spp. — Gram-negative, Enterobacterales
        // Band 6 (x0.625)
        map.insert("bacteria_shigella_spp._mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.001      ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_shigella_spp._mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.001      ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_shigella_spp._mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.001      ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_shigella_spp._mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.001      ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_shigella_spp._mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.001      ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_shigella_spp._mechanism_enzyme_kpc_emergence_rate".to_string(), 0.001       ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_shigella_spp._mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.001        ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_shigella_spp._mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.001        ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_shigella_spp._mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_25   ); // classes: chl
        map.insert("bacteria_shigella_spp._mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.9    ) ; // classes: ag
        map.insert("bacteria_shigella_spp._mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_target_site_erm_b_emergence_rate".to_string(), 0.8    ); // classes: mls
        map.insert("bacteria_shigella_spp._mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_5     ); // classes: fq
        map.insert("bacteria_shigella_spp._mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_5      ); // classes: fq
        map.insert("bacteria_shigella_spp._mechanism_protection_qnr_emergence_rate".to_string(), 0.000_5      ); // classes: fq
        map.insert("bacteria_shigella_spp._mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_5      ); // classes: fq, tet, chl
        map.insert("bacteria_shigella_spp._mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_global_efflux_pump_emergence_rate".to_string(), 0.9     ); // classes: fq, mls, tet, chl
        map.insert("bacteria_shigella_spp._mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_03    ); // classes: none (currently zeroed)
        map.insert("bacteria_shigella_spp._mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_3    ); // classes: poly; Colistin ~20%
        map.insert("bacteria_shigella_spp._mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.003     ); // classes: sulf
        map.insert("bacteria_shigella_spp._mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_shigella_spp._mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // classes: other (fosfomycin)
        map.insert("bacteria_shigella_spp._mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.04    ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_shigella_spp._mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_protection_tet_m_emergence_rate".to_string(), 0.03    ); // classes: tet
        map.insert("bacteria_shigella_spp._mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.9     ); // classes: ag
        map.insert("bacteria_shigella_spp._mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_shigella_spp._mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.9     ); // classes: mls
        map.insert("bacteria_shigella_spp._mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.03    ); // classes: tet
        map.insert("bacteria_shigella_spp._mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.001      ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_shigella_spp._mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.001      ); // classes: pen, mls, tet, chl; low broad-efflux proxy
        map.insert("bacteria_shigella_spp._mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // Y. enterocolitica — Gram-negative, Enterobacterales
        // Band 8 (x300)
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_000_000_1  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_000_000_03 ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_000_000_03 ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_000_000_03 ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: chl
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_000_000_03 ); // classes: ag
        map.insert("bacteria_yersinia_enterocolitica_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: fq
        map.insert("bacteria_yersinia_enterocolitica_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: fq
        map.insert("bacteria_yersinia_enterocolitica_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: fq
        map.insert("bacteria_yersinia_enterocolitica_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: fq, tet, chl
        map.insert("bacteria_yersinia_enterocolitica_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: fq, mls, tet, chl
        map.insert("bacteria_yersinia_enterocolitica_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_000_1 ); // classes: none (currently zeroed)
        map.insert("bacteria_yersinia_enterocolitica_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_000_000_3  ); // classes: poly
        map.insert("bacteria_yersinia_enterocolitica_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: sulf
        map.insert("bacteria_yersinia_enterocolitica_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_000_001  ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: other (fosfomycin)
        map.insert("bacteria_yersinia_enterocolitica_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_000_03 ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_yersinia_enterocolitica_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_000_000_3 ); // classes: tet
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_000_000_15); // classes: ag
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_000_000_15); // classes: tet
        map.insert("bacteria_yersinia_enterocolitica_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_yersinia_enterocolitica_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_000_000_3); // classes: pen, mls, tet, chl; tiny broad-efflux proxy
        map.insert("bacteria_yersinia_enterocolitica_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // ======================================================================
        // Non-fermenting Gram-Negatives
        // ======================================================================
        // P. aeruginosa — Gram-negative, NonFermenter
        // Band 8 (x37.5)
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_03    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_03    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_03    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_03   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_03   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_03    ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_03  ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_03    ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_4   ); // classes: chl
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_4   ); // classes: ag
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_03    ); // classes: mls
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_03    ); // classes: oxa, mls, chl
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.001      ); // classes: fq
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.001   ); // classes: fq
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_protection_qnr_emergence_rate".to_string(), 0.001    ); // classes: fq
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.000_5     ); // classes: ceph, carb, fq, ag
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_3    ); // classes: fq, mls, tet, chl
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.000_5    ); // classes: carb
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_3  ); // classes: none (currently zeroed)
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.00_3    ); // classes: poly
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.003     ); // classes: sulf
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_03    ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.01      ); // classes: other (fosfomycin)
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.001     ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_protection_tet_m_emergence_rate".to_string(), 0.004     ); // classes: tet
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.001     ); // classes: ag; AMEs extremely common in PA
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_02   ); // classes: tet; TetABC present but MexAB dominates
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, mono; deactivated
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // classes: pen, mls, tet, chl; deactivated
        map.insert("bacteria_pseudomonas_aeruginosa_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; deactivated

       // A. baumannii — Gram-negative, NonFermenter
        // Band 8 (x100)
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, mono
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, mono
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, mono
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, mono
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, mono
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, carb
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_cat_emergence_rate".to_string(), 0.01 ); // classes: chl
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 10.0  ); // classes: ag
        map.insert("bacteria_acinetobacter_baumannii_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.09 ); // classes: fq
        map.insert("bacteria_acinetobacter_baumannii_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.09 ); // classes: fq
        map.insert("bacteria_acinetobacter_baumannii_mechanism_protection_qnr_emergence_rate".to_string(), 0.09 ); // classes: fq
        map.insert("bacteria_acinetobacter_baumannii_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.025  ); // classes: fq, mls, tet, chl
        map.insert("bacteria_acinetobacter_baumannii_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_1); // classes: none (currently zeroed)
        map.insert("bacteria_acinetobacter_baumannii_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.01 ); // classes: poly; Colistin 0%->12%
        map.insert("bacteria_acinetobacter_baumannii_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.025 ); // classes: sulf
        map.insert("bacteria_acinetobacter_baumannii_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.003_5   ); // classes: other (fosfomycin)
        map.insert("bacteria_acinetobacter_baumannii_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.045); // classes: other (rifampicin, fidaxomicin); Rifampicin 0%->55%
        map.insert("bacteria_acinetobacter_baumannii_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_protection_tet_m_emergence_rate".to_string(), 0.05 ); // classes: tet
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 10.0 ); // classes: ag; AMEs frequent in MDR AB
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.007); // classes: carb; OXA-23/40/58 primary carbapenem R in AB
        map.insert("bacteria_acinetobacter_baumannii_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.002_5); // classes: tet
        map.insert("bacteria_acinetobacter_baumannii_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.007); // classes: pen, bli, ceph, mono
        map.insert("bacteria_acinetobacter_baumannii_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_acinetobacter_baumannii_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // S. maltophilia — Gram-negative, NonFermenter
        // Band 9 (x500)
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_02    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_02    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_02    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.05   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_2   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_02    ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.1   ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_02    ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_5); // classes: chl; chlor 68.85% vs target 50%
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.01   ); // classes: ag
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.1); // classes: mls; clindamycin 0% vs target 95%
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_000_05); // classes: fq; FQ 88% vs target 45-55%
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_05); // classes: fq; FQ broad reduction
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_02   ); // classes: fq
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_1); // classes: fq, mls, tet, chl; inflating FQ+tet+chlor
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_000_000_1); // classes: none (currently zeroed)
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.05); // classes: poly; colistin 0% vs target 70%
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.005); // classes: sulf; trim/sulf 34% vs target 15%
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 10.0); // classes: other (metronidazole, nitrofurantoin, furazolidone); nitrofurantoin 0% vs target 90%
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.03      ); // classes: other (fosfomycin)
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_2   ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_000_5); // classes: tet; tet 68% vs target 30-60%
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.05   ); // classes: ag
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_000_002   ); // classes: tet
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_stenotrophomonas_maltophilia_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // B. cepacia complex — Gram-negative, NonFermenter
        // Band 9 (x750)
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_007_5); // classes: pen, bli, ceph, mono
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_007_5); // classes: pen, bli, ceph, mono
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_007_5); // classes: pen, bli, ceph, mono
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_38); // classes: pen, bli, ceph, mono
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_075); // classes: pen, bli, ceph, mono
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_007_5); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_037); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_007_5); // classes: pen, bli, ceph, carb
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_037); // classes: chl
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_007_5); // classes: ag
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_38); // classes: fq
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_075); // classes: fq
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_007_5); // classes: fq
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_75); // classes: fq, mls, tet, chl
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_38); // classes: none (currently zeroed)
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_007_5); // classes: poly
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_075); // classes: sulf
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_037); // classes: other (fosfomycin)
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_007_5); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_037); // classes: tet
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_000_000_05); // classes: ag
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_000_000_05); // classes: tet
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_burkholderia_cepacia_complex_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // ======================================================================
        // Other Gram-Negatives
        // ======================================================================
        // V. cholerae — Gram-negative, EntericPathogen
        // Band 7 (x30)
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_003); // classes: pen, bli, ceph, mono
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_003); // classes: pen, bli, ceph, mono
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_001  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_001  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_001  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_000_1 ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_001  ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_000_1 ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_1); // classes: chl
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_000_03); // classes: ag
        map.insert("bacteria_vibrio_cholerae_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_03); // classes: fq
        map.insert("bacteria_vibrio_cholerae_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_015); // classes: fq
        map.insert("bacteria_vibrio_cholerae_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_003); // classes: fq
        map.insert("bacteria_vibrio_cholerae_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_003); // classes: fq, tet, chl
        map.insert("bacteria_vibrio_cholerae_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_003); // classes: fq, mls, tet, chl
        map.insert("bacteria_vibrio_cholerae_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_001_5); // classes: none (currently zeroed)
        map.insert("bacteria_vibrio_cholerae_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_000_3); // classes: poly
        map.insert("bacteria_vibrio_cholerae_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_15); // classes: sulf
        map.insert("bacteria_vibrio_cholerae_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_000_3); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_000_3); // classes: other (fosfomycin)
        map.insert("bacteria_vibrio_cholerae_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_3); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_vibrio_cholerae_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_05); // classes: tet
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_000_000_05); // classes: ag
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_15); // classes: mls
        map.insert("bacteria_vibrio_cholerae_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_000_000_05); // classes: tet
        map.insert("bacteria_vibrio_cholerae_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_vibrio_cholerae_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_000_001); // classes: pen, mls, tet, chl; low broad-efflux proxy
        map.insert("bacteria_vibrio_cholerae_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // C. jejuni — Gram-negative, Helicobacter group
        // Band 6 (x1.2)
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_6    ); // classes: chl
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_01   ); // classes: mls; low macrolide/lincosamide resistance
        map.insert("bacteria_campylobacter_jejuni_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_01    ); // classes: oxa, mls, chl; supports clindamycin and chloramphenicol tail
        map.insert("bacteria_campylobacter_jejuni_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.005      ); // classes: fq
        map.insert("bacteria_campylobacter_jejuni_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.008      ); // classes: fq
        map.insert("bacteria_campylobacter_jejuni_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_01       ); // classes: fq, mls, tet, chl
        map.insert("bacteria_campylobacter_jejuni_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_02   ); // classes: none (currently zeroed)
        map.insert("bacteria_campylobacter_jejuni_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.6  ); // classes: sulf
        map.insert("bacteria_campylobacter_jejuni_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.6      ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_campylobacter_jejuni_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_3   ); // classes: tet
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.03  ); // classes: ag
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_campylobacter_jejuni_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_003    ); // classes: mls; major macrolide R mechanism
        map.insert("bacteria_campylobacter_jejuni_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_01    ); // classes: tet
        map.insert("bacteria_campylobacter_jejuni_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, mono; PBP mosaic: not a significant mechanism in Campylobacter
        map.insert("bacteria_campylobacter_jejuni_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_2   ); // classes: pen, mls, tet, chl; CmeABC efflux: major efflux pump, functionally analogous to mtrCDE
        map.insert("bacteria_campylobacter_jejuni_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // H. pylori — Gram-negative, Helicobacter group
        // Band 6 (x0.6)
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_cat_emergence_rate".to_string(), 0.0); // classes: chl
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_target_site_erm_b_emergence_rate".to_string(), 100000000.0  ); // classes: mls
        map.insert("bacteria_helicobacter_pylori_mechanism_target_site_cfr_emergence_rate".to_string(), 100000000.0  ); // classes: oxa, mls, chl
        map.insert("bacteria_helicobacter_pylori_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 100000000.0 ); // classes: fq
        map.insert("bacteria_helicobacter_pylori_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 100000000.0 ); // classes: fq
        map.insert("bacteria_helicobacter_pylori_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_global_efflux_pump_emergence_rate".to_string(), 100000000.00  ); // classes: fq, mls, tet, chl
        map.insert("bacteria_helicobacter_pylori_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_global_porin_loss_emergence_rate".to_string(), 100000000.00  ); // classes: none (currently zeroed)
        map.insert("bacteria_helicobacter_pylori_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 100000000.00     ); // classes: sulf
        map.insert("bacteria_helicobacter_pylori_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 100000000.00 ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_mutation_rpo_b_emergence_rate".to_string(), 100000000.00); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_helicobacter_pylori_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_protection_tet_m_emergence_rate".to_string(), 100000000.00 ); // classes: tet
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_bla_z_emergence_rate".to_string(), 100000000.00 ); // classes: pen; PBP1A proxy via blaZ; targets amox/amp only
        map.insert("bacteria_helicobacter_pylori_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 100000000.00 ); // classes: mls; primary clarithromycin R mechanism
        map.insert("bacteria_helicobacter_pylori_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 100000000.00 ); // classes: pen, bli, ceph, mono; low PBP1A-type amoxicillin resistance proxy
        map.insert("bacteria_helicobacter_pylori_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_helicobacter_pylori_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // N. gonorrhoeae — Gram-negative, Fastidious
        // Band 6 (x1.5)
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, mono; GC does not produce ESBLs
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, mono; TEM-1 modeled via blaZ instead
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, mono; GC does not produce SHV
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, mono; GC does not produce AmpC
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, mono; GC does not produce AmpC
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, carb, mono; GC does not produce KPC
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, carb, mono; GC does not produce MBLs
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // classes: pen, bli, ceph, carb; GC does not produce OXA-48
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_3      ); // classes: chl
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 10.0        ); // classes: ag
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.002    ); // classes: mls
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_target_site_cfr_emergence_rate".to_string(), 0.002           ); // classes: oxa, mls, chl
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.9     ); // classes: fq
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.9      ); // classes: fq
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_protection_qnr_emergence_rate".to_string(), 0.9       ); // classes: fq
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0        ); // classes: fq, tet, chl
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.3            ); // classes: fq, mls, tet, chl
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0    ); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_01    ); // classes: none (currently zeroed)
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_3          ); // classes: poly
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.03   ); // classes: sulf
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.03          ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_03        ); // classes: other (fosfomycin)
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.3        ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_protection_tet_m_emergence_rate".to_string(), 0.5          ); // classes: tet
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 10.0   ); // classes: ag
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.000_4   ); // classes: pen; TEM-1 penicillinase; all penicillins
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.002    ); // classes: mls; GC macrolide target-site resistance
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.002           ); // classes: tet
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.004    ); // classes: pen, bli, ceph, mono; PBP mosaic: penA mosaic alleles, major contributor to pen/ceph R
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_5   ); // classes: pen, mls, tet, chl; mtrCDE efflux: common, broad GC resistance driver
        map.insert("bacteria_neisseria_gonorrhoeae_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0 ); // classes: broad placeholder

        // N. meningitidis — Gram-negative, Fastidious
        // Band 8 (x250)
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_000_01 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_000_01); // classes: pen, bli, ceph, mono
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_000_01 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_000_01 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_000_01 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_002); // classes: chl
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_000_1 ); // classes: ag
        map.insert("bacteria_neisseria_meningitidis_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_000_1); // classes: mls
        map.insert("bacteria_neisseria_meningitidis_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_2) ; // classes: oxa, mls, chl
        map.insert("bacteria_neisseria_meningitidis_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_03 ); // classes: fq
        map.insert("bacteria_neisseria_meningitidis_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_01 ); // classes: fq
        map.insert("bacteria_neisseria_meningitidis_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_000_3 ); // classes: fq
        map.insert("bacteria_neisseria_meningitidis_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_000_3); // classes: fq, tet, chl
        map.insert("bacteria_neisseria_meningitidis_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_002); // classes: fq, mls, tet, chl
        map.insert("bacteria_neisseria_meningitidis_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_5); // classes: none (currently zeroed)
        map.insert("bacteria_neisseria_meningitidis_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_000_3 ); // classes: poly
        map.insert("bacteria_neisseria_meningitidis_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_005); // classes: sulf
        map.insert("bacteria_neisseria_meningitidis_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_003  ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_003); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_neisseria_meningitidis_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_005); // classes: tet
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_neisseria_meningitidis_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_000_1); // classes: mls; rare but not impossible
        map.insert("bacteria_neisseria_meningitidis_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.000_005); // classes: tet; low-probability tetracycline efflux
        map.insert("bacteria_neisseria_meningitidis_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_001); // classes: pen, bli, ceph, mono; PBP mosaic: penA alterations, intermediate penicillin R
        map.insert("bacteria_neisseria_meningitidis_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_000_001); // classes: pen, mls, tet, chl; mtrCDE-like efflux: present but less clinically significant
        map.insert("bacteria_neisseria_meningitidis_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

        // M. catarrhalis — Gram-negative, Fastidious
        // Band 7 (x18.8)
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_000_2 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_001  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_000_02 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_000_5 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_000_2 ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_5 ); // classes: chl
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_000_02 ); // classes: ag
        map.insert("bacteria_moraxella_catarrhalis_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_005  ); // classes: mls
        map.insert("bacteria_moraxella_catarrhalis_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_05 ); // classes: oxa, mls, chl
        map.insert("bacteria_moraxella_catarrhalis_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_002  ); // classes: fq
        map.insert("bacteria_moraxella_catarrhalis_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_5 ); // classes: fq
        map.insert("bacteria_moraxella_catarrhalis_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_000_2 ); // classes: fq
        map.insert("bacteria_moraxella_catarrhalis_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_000_5 ); // classes: fq, tet, chl
        map.insert("bacteria_moraxella_catarrhalis_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_005  ); // classes: fq, mls, tet, chl
        map.insert("bacteria_moraxella_catarrhalis_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_2 ); // classes: none (currently zeroed)
        map.insert("bacteria_moraxella_catarrhalis_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_000_02 ); // classes: poly
        map.insert("bacteria_moraxella_catarrhalis_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_01   ); // classes: sulf
        map.insert("bacteria_moraxella_catarrhalis_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_000_2); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_2 ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_moraxella_catarrhalis_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_01   ); // classes: tet
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_moraxella_catarrhalis_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_00_5); // classes: pen, bli, ceph, mono; PBP mosaic: some PBP modifications reported
        map.insert("bacteria_moraxella_catarrhalis_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_000_5); // classes: pen, mls, tet, chl; mtrCDE-like efflux: contributes to macrolide/pen R
        map.insert("bacteria_moraxella_catarrhalis_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

               // H. influenzae — Gram-negative, Fastidious
        // Band 7 (x18.8)
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_001    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_001  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_001   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_001    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_001    ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.000_001 ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.000_001 ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.000_001 ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_1    ); // classes: chl
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.05     ); // classes: ag
        map.insert("bacteria_haemophilus_influenzae_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_03  ); // classes: mls
        map.insert("bacteria_haemophilus_influenzae_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_03 ); // classes: oxa, mls, chl
        map.insert("bacteria_haemophilus_influenzae_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_3    ); // classes: fq
        map.insert("bacteria_haemophilus_influenzae_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_08    ); // classes: fq
        map.insert("bacteria_haemophilus_influenzae_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_1    ); // classes: fq
        map.insert("bacteria_haemophilus_influenzae_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0   ); // classes: fq, tet, chl
        map.insert("bacteria_haemophilus_influenzae_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_3  ); // classes: fq, mls, tet, chl
        map.insert("bacteria_haemophilus_influenzae_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_3    ); // classes: none 
        map.insert("bacteria_haemophilus_influenzae_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_002     ); // classes: poly
        map.insert("bacteria_haemophilus_influenzae_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_8    ); // classes: sulf
        map.insert("bacteria_haemophilus_influenzae_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_02      ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.015     ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_haemophilus_influenzae_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_8     ); // classes: tet
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.05     ); // classes: ag
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.000_001_5   ); // classes: pen
        map.insert("bacteria_haemophilus_influenzae_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_haemophilus_influenzae_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_003   ); // classes: mls
        map.insert("bacteria_haemophilus_influenzae_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // classes: tet; tet signal is already high
        map.insert("bacteria_haemophilus_influenzae_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_002     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_haemophilus_influenzae_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_002     ); // classes: pen, mls, tet, chl
        map.insert("bacteria_haemophilus_influenzae_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0


        // L. pneumophila — Gram-negative, Fastidious
        // Band 8 (x100)
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_000_3); // classes: pen, bli, ceph, mono
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_000_3); // classes: pen, bli, ceph, mono
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_000_3); // classes: pen, bli, ceph, mono
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_000_3); // classes: pen, bli, ceph, mono
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_000_3); // classes: pen, bli, ceph, mono
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_003); // classes: chl
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_000_3); // classes: ag
        map.insert("bacteria_legionella_pneumophila_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_03 ); // classes: mls
        map.insert("bacteria_legionella_pneumophila_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_3); // classes: oxa, mls, chl
        map.insert("bacteria_legionella_pneumophila_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_03); // classes: fq
        map.insert("bacteria_legionella_pneumophila_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_03 ); // classes: fq
        map.insert("bacteria_legionella_pneumophila_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_000_1); // classes: fq
        map.insert("bacteria_legionella_pneumophila_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_000_3); // classes: fq, tet, chl
        map.insert("bacteria_legionella_pneumophila_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_03 ); // classes: fq, mls, tet, chl
        map.insert("bacteria_legionella_pneumophila_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_001); // classes: none (currently zeroed)
        map.insert("bacteria_legionella_pneumophila_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.000_000_3); // classes: poly
        map.insert("bacteria_legionella_pneumophila_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_003); // classes: sulf
        map.insert("bacteria_legionella_pneumophila_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_000_3); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_003); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_legionella_pneumophila_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_03 ); // classes: tet
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // classes: ag; disabled here; aminoglycosides are not a useful clinical lever
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_000_000_05); // classes: mls; rare but possible
        map.insert("bacteria_legionella_pneumophila_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_legionella_pneumophila_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

        // ======================================================================
        // Gram-Positive Bacteria — Staphylococci
        // ======================================================================
        // S. aureus — Gram-positive, Staphylococcus
        // Band 6 (x1.5)
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_01     ); // classes: chl; chlor target 30%
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.3            ); // classes: ag
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.000_05  ); // classes: pen, bli, ceph, carb; driving penicillins ~35-60%
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_van_a_emergence_rate".to_string(), 0.000_003    ); // classes: glyc; targets ~10% for vanco/teico
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_van_b_emergence_rate".to_string(), 0.000_000_3        ); // classes: glyc
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0013      ); // classes: mls; macrolides ~40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_3      ); // classes: oxa, mls, chl; Linezolid 1%
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.05   ); // classes: fq; FQ target ~40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.015  ); // classes: fq; FQ target ~40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_protection_qnr_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_1   ); // classes: fq, mls, tet, chl; bumps FQs, Tet
        map.insert("bacteria_staphylococcus_aureus_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // classes: none (currently zeroed); broad base
        map.insert("bacteria_staphylococcus_aureus_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.005   ); // classes: sulf; trim/sulf target 35%
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.005        ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.005       ); // classes: other (fosfomycin)
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.005        ); // classes: other (daptomycin)
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.005        ); // classes: other (rifampicin, fidaxomicin); Rifampicin ~40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_protection_fus_b_emergence_rate".to_string(), 0.005        ); // classes: other (fusidic acid); Fusidic acid ~20%
        map.insert("bacteria_staphylococcus_aureus_mechanism_protection_tet_m_emergence_rate".to_string(), 0.005    ); // classes: tet; tet target 20-30%
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.3   ); // classes: ag; AG 40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.000_05); // classes: pen; blaZ
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; disabled here; 23S-mediated resistance is possible
        map.insert("bacteria_staphylococcus_aureus_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0


       // S. aureus — Gram-positive, Staphylococcus
        // Band 6 (x1.5)
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0  ); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_01   ); // classes: chl; chlor target 30%
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_01           ); // classes: ag
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.000_000_1); // classes: pen, bli, ceph, carb; driving penicillins ~35-60%
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_van_a_emergence_rate".to_string(), 0.000_000_01    ); // classes: glyc; targets ~10% for vanco/teico
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_van_b_emergence_rate".to_string(), 0.000_000_000_05    ); // classes: glyc
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_004  ); // classes: mls; macrolides ~40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_04  ); // classes: oxa, mls, chl; Linezolid 1%
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_015  ); // classes: fq; FQ target ~40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_015  ); // classes: fq; FQ target ~40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_protection_qnr_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_004  ); // classes: fq, mls, tet, chl; bumps FQs, Tet
        map.insert("bacteria_staphylococcus_aureus_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // classes: none (currently zeroed); broad base
        map.insert("bacteria_staphylococcus_aureus_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_01   ); // classes: sulf; trim/sulf target 35%
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_01        ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.000_01       ); // classes: other (fosfomycin)
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.000_01        ); // classes: other (daptomycin)
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_01        ); // classes: other (rifampicin, fidaxomicin); Rifampicin ~40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_protection_fus_b_emergence_rate".to_string(), 0.000_01        ); // classes: other (fusidic acid); Fusidic acid ~20%
        map.insert("bacteria_staphylococcus_aureus_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_004   ); // classes: tet; tet target 20-30%
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_01    ); // classes: ag; AG 40%
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.000_000_01  ); // classes: pen; blaZ
        map.insert("bacteria_staphylococcus_aureus_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_aureus_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; disabled here; 23S-mediated resistance is possible
        map.insert("bacteria_staphylococcus_aureus_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0

        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_001); // classes: chl
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.000_001); // classes: pen, bli, ceph, carb
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_target_site_van_a_emergence_rate".to_string(), 0.000_000_01); // classes: glyc
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_target_site_van_b_emergence_rate".to_string(), 0.000_000_01); // classes: glyc
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_01); // classes: mls
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_1); // classes: oxa, mls, chl
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_01); // classes: fq
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_001); // classes: fq
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_001); // classes: fq, mls, tet, chl
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_01); // classes: sulf
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.000_001); // classes: other (daptomycin)
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.001); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_protection_fus_b_emergence_rate".to_string(), 0.000_003); // classes: other (fusidic acid)
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_03); // classes: tet
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_000_001); // classes: ag
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.000_000_1); // classes: pen; common in CoNS
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; disabled here; 23S-mediated resistance is possible
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_staphylococcus_epidermidis_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

        // ======================================================================
        // Streptococci
        // ======================================================================
         // S. pneumoniae — Gram-positive, Streptococcus
        // Band 6 (x0.75)
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_12     ); // classes: chl
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.000_000_2 ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_001      ); // classes: mls
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_015    ); // classes: fq
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_015    ); // classes: fq
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_015    ); // classes: fq, mls, tet, chl
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_015  ); // classes: none (currently zeroed)
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.001          ); // classes: sulf
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_5 ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_5   ); // classes: tet
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.000_000_1  ); // classes: pen
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_000_2  ); // classes: mls
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_000_1     ); // classes: pen, bli, ceph, mono; PBP mosaic remains the primary pen-R lever, but was strongly overcalled
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // classes: pen, mls, tet, chl; mtrCDE: not relevant (GramPositive group excluded)
        map.insert("bacteria_streptococcus_pneumoniae_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // S. pyogenes — Gram-positive, Streptococcus
        // Band 7 (x7.5)
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_38); // classes: chl
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_007_5); // classes: mls
        map.insert("bacteria_streptococcus_pyogenes_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_007_5); // classes: oxa, mls, chl
        map.insert("bacteria_streptococcus_pyogenes_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_000_38); // classes: fq
        map.insert("bacteria_streptococcus_pyogenes_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_075); // classes: fq
        map.insert("bacteria_streptococcus_pyogenes_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_75); // classes: fq, mls, tet, chl
        map.insert("bacteria_streptococcus_pyogenes_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_003_8); // classes: sulf
        map.insert("bacteria_streptococcus_pyogenes_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.000_000_007_5); // classes: other (daptomycin)
        map.insert("bacteria_streptococcus_pyogenes_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_007_5); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_streptococcus_pyogenes_mechanism_protection_fus_b_emergence_rate".to_string(), 0.000_000_007_5); // classes: other (fusidic acid)
        map.insert("bacteria_streptococcus_pyogenes_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_007_5); // classes: tet
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_000_000_05); // classes: mls
        map.insert("bacteria_streptococcus_pyogenes_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_pyogenes_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // S. agalactiae — Gram-positive, Streptococcus
        // Band 8 (x50)
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_3); // classes: chl
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.000_000_01); // classes: pen, bli, ceph, carb
        map.insert("bacteria_streptococcus_agalactiae_mechanism_target_site_van_a_emergence_rate".to_string(), 0.000_000_01); // classes: glyc
        map.insert("bacteria_streptococcus_agalactiae_mechanism_target_site_van_b_emergence_rate".to_string(), 0.000_000_01); // classes: glyc
        map.insert("bacteria_streptococcus_agalactiae_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_05); // classes: mls
        map.insert("bacteria_streptococcus_agalactiae_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_01); // classes: oxa, mls, chl
        map.insert("bacteria_streptococcus_agalactiae_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_001); // classes: fq
        map.insert("bacteria_streptococcus_agalactiae_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_001  ); // classes: fq
        map.insert("bacteria_streptococcus_agalactiae_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_001); // classes: fq, mls, tet, chl
        map.insert("bacteria_streptococcus_agalactiae_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_01 ); // classes: sulf
        map.insert("bacteria_streptococcus_agalactiae_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.000_000_1); // classes: other (daptomycin)
        map.insert("bacteria_streptococcus_agalactiae_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_1); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_streptococcus_agalactiae_mechanism_protection_fus_b_emergence_rate".to_string(), 0.000_000_1); // classes: other (fusidic acid)
        map.insert("bacteria_streptococcus_agalactiae_mechanism_protection_tet_m_emergence_rate".to_string(), 0.001   ); // classes: tet
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_000_000_05); // classes: mls
        map.insert("bacteria_streptococcus_agalactiae_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_001); // classes: pen, bli, ceph, mono; tiny seed for rare reduced beta-lactam susceptibility
        map.insert("bacteria_streptococcus_agalactiae_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_streptococcus_agalactiae_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // ======================================================================
        // Enterococci
        // ======================================================================
        // E. faecalis — Gram-positive, Enterococcus
        // Band 8 (x60)
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_2  ); // classes: chl; chlor 35% vs target 25%
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.000_02   ); // classes: pen, bli, ceph, carb; pen 29% vs target 10%
        map.insert("bacteria_enterococcus_faecalis_mechanism_target_site_van_a_emergence_rate".to_string(), 0.02   ); // classes: glyc; vanc 29% vs target 5%
        map.insert("bacteria_enterococcus_faecalis_mechanism_target_site_van_b_emergence_rate".to_string(), 0.02   ); // classes: glyc; vanc 29% vs target 5%
        map.insert("bacteria_enterococcus_faecalis_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.002   ); // classes: mls; macrolide 47% vs target 35%
        map.insert("bacteria_enterococcus_faecalis_mechanism_target_site_cfr_emergence_rate".to_string(), 0.003  ); // classes: oxa, mls, chl
        map.insert("bacteria_enterococcus_faecalis_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.002   ); // classes: fq; FQ 74% vs target 35%
        map.insert("bacteria_enterococcus_faecalis_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_3 ); // classes: fq; FQ broad reduction
        map.insert("bacteria_enterococcus_faecalis_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_3 ); // classes: fq, mls, tet, chl; cuts FQ+tet+chlor
        map.insert("bacteria_enterococcus_faecalis_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_1  ); // classes: none (currently zeroed); low-probability background permeability effect
        map.insert("bacteria_enterococcus_faecalis_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.02   ); // classes: sulf
        map.insert("bacteria_enterococcus_faecalis_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.02   ); // classes: other (metronidazole, nitrofurantoin, furazolidone); rare but not impossible
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.005 ); // classes: other (daptomycin)
        map.insert("bacteria_enterococcus_faecalis_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.005 ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_enterococcus_faecalis_mechanism_protection_fus_b_emergence_rate".to_string(), 0.000_5  ); // classes: other (fusidic acid)
        map.insert("bacteria_enterococcus_faecalis_mechanism_protection_tet_m_emergence_rate".to_string(), 0.002    ); // classes: tet; tet 70% vs target 15%
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.000_02    ); // classes: ag; HLAR in enterococci is well-documented
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_3  ); // classes: mls; rare but not impossible
        map.insert("bacteria_enterococcus_faecalis_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.000_015   ); // classes: pen, bli, ceph, mono; low-level PBP contribution, much weaker than E. faecium
        map.insert("bacteria_enterococcus_faecalis_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecalis_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

        // E. faecium — Gram-positive, Enterococcus
        // Band 8 (x100)
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_cat_emergence_rate".to_string(),  0.002    ); // classes: chl
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.004      ); // classes: ag
        map.insert("bacteria_enterococcus_faecium_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(),  0.0 ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_enterococcus_faecium_mechanism_target_site_van_a_emergence_rate".to_string(),   0.01    ); // classes: glyc
        map.insert("bacteria_enterococcus_faecium_mechanism_target_site_van_b_emergence_rate".to_string(),   0.01      ); // classes: glyc
        map.insert("bacteria_enterococcus_faecium_mechanism_target_site_erm_b_emergence_rate".to_string(),   0.01        ); // classes: mls
        map.insert("bacteria_enterococcus_faecium_mechanism_target_site_cfr_emergence_rate".to_string(),  0.01     ); // classes: oxa, mls, chl
        map.insert("bacteria_enterococcus_faecium_mechanism_mutation_gyra_primary_emergence_rate".to_string(),   0.000_5  ); // classes: fq
        map.insert("bacteria_enterococcus_faecium_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(),   0.002   ) ; // classes: fq
        map.insert("bacteria_enterococcus_faecium_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.005   ); // classes: fq, mls, tet, chl
        map.insert("bacteria_enterococcus_faecium_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.005   ); // classes: sulf
        map.insert("bacteria_enterococcus_faecium_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.4    ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.4    ); // classes: other (fosfomycin)
        map.insert("bacteria_enterococcus_faecium_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.015  ); // classes: other (daptomycin)
        map.insert("bacteria_enterococcus_faecium_mechanism_mutation_rpo_b_emergence_rate".to_string(),   0.15 ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_enterococcus_faecium_mechanism_protection_fus_b_emergence_rate".to_string(),  0.015  ); // classes: other (fusidic acid)
        map.insert("bacteria_enterococcus_faecium_mechanism_protection_tet_m_emergence_rate".to_string(),  0.005   ); // classes: tet
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.03    ); // classes: ag
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.001_5  ); // classes: mls; rare but not impossible
        map.insert("bacteria_enterococcus_faecium_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_enterococcus_faecium_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.015  ); // classes: pen, bli, ceph, mono; PBP mosaic: PBP5 mutations → intrinsic ampicillin R
        map.insert("bacteria_enterococcus_faecium_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.001_5  ); // classes: pen, mls, tet, chl; mtrCDE: not relevant (GramPositive group excluded)
        map.insert("bacteria_enterococcus_faecium_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; deactivated

        // ======================================================================
        // Other Gram-Positives and Anaerobes
        // ======================================================================
        // L. monocytogenes — Gram-positive, Listeria
        // Band 10 (x3750)
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_19); // classes: chl
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_target_site_van_a_emergence_rate".to_string(), 0.000_003_8); // classes: glyc
        map.insert("bacteria_listeria_monocytogenes_mechanism_target_site_van_b_emergence_rate".to_string(), 0.000_003_8); // classes: glyc
        map.insert("bacteria_listeria_monocytogenes_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_38); // classes: mls
        map.insert("bacteria_listeria_monocytogenes_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_003_8); // classes: oxa, mls, chl
        map.insert("bacteria_listeria_monocytogenes_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_19); // classes: fq
        map.insert("bacteria_listeria_monocytogenes_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_038); // classes: fq
        map.insert("bacteria_listeria_monocytogenes_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_19); // classes: fq, mls, tet, chl
        map.insert("bacteria_listeria_monocytogenes_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_38); // classes: sulf
        map.insert("bacteria_listeria_monocytogenes_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.000_038); // classes: other (daptomycin)
        map.insert("bacteria_listeria_monocytogenes_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_038); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_listeria_monocytogenes_mechanism_protection_fus_b_emergence_rate".to_string(), 0.000_038); // classes: other (fusidic acid)
        map.insert("bacteria_listeria_monocytogenes_mechanism_protection_tet_m_emergence_rate".to_string(), 0.001_9); // classes: tet
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; 
        map.insert("bacteria_listeria_monocytogenes_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_listeria_monocytogenes_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

        // C. difficile — Anaerobe
        // Band 8 (x60)
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_006); // classes: chl
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_000_1 ); // classes: ag
        map.insert("bacteria_clostridioides_difficile_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_06); // classes: mls
        map.insert("bacteria_clostridioides_difficile_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_001  ); // classes: oxa, mls, chl
        map.insert("bacteria_clostridioides_difficile_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_06); // classes: fq
        map.insert("bacteria_clostridioides_difficile_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_01 ); // classes: fq
        map.insert("bacteria_clostridioides_difficile_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_01 ); // classes: fq, mls, tet, chl
        map.insert("bacteria_clostridioides_difficile_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_006); // classes: sulf
        map.insert("bacteria_clostridioides_difficile_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_12); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_06); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_clostridioides_difficile_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_06); // classes: tet
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; 
        map.insert("bacteria_clostridioides_difficile_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_clostridioides_difficile_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

         // B. fragilis — Anaerobe
        // Band 8 (x37.5)
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.005  ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.005   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.005     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.001     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.001     ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.001     ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.001   ); // classes: pen, bli, ceph, carb, mono
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.001    ); // classes: pen, bli, ceph, carb
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_02   ); // classes: chl
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 1.0       ); // classes: ag
        map.insert("bacteria_bacteroides_fragilis_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.02   ); // classes: mls
        map.insert("bacteria_bacteroides_fragilis_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_2     ); // classes: oxa, mls, chl
        map.insert("bacteria_bacteroides_fragilis_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_5   ); // classes: fq
        map.insert("bacteria_bacteroides_fragilis_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.01      ); // classes: fq
        map.insert("bacteria_bacteroides_fragilis_mechanism_protection_qnr_emergence_rate".to_string(), 0.000_2     ); // classes: fq
        map.insert("bacteria_bacteroides_fragilis_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_1    ); // classes: fq, tet, chl
        map.insert("bacteria_bacteroides_fragilis_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0   ); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_1   ); // classes: fq, mls, tet, chl
        map.insert("bacteria_bacteroides_fragilis_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0 ); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_01); // classes: none (currently zeroed)
        map.insert("bacteria_bacteroides_fragilis_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.002     ); // classes: poly
        map.insert("bacteria_bacteroides_fragilis_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.005    ); // classes: sulf
        map.insert("bacteria_bacteroides_fragilis_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_5  ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.002      ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_bacteroides_fragilis_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_protection_tet_m_emergence_rate".to_string(), 0.01 ); // classes: tet
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 1.0  ); // classes: ag
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; 
        map.insert("bacteria_bacteroides_fragilis_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.001   ); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bacteroides_fragilis_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bacteroides_fragilis_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

        // ======================================================================
        // Fastidious and Atypical Bacteria
        // ======================================================================
        // B. pertussis — Gram-negative, Fastidious
        // Band 7 (x7.5)
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.000_000_000_1); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.000_000_000_1); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.000_000_000_1); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.000_000_000_1); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.000_000_000_1); // classes: pen, bli, ceph, mono
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_000_003_8); // classes: chl
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.000_000_000_1); // classes: ag
        map.insert("bacteria_bordetella_pertussis_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // classes: mls
        map.insert("bacteria_bordetella_pertussis_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_000_05); // classes: oxa, mls, chl
        map.insert("bacteria_bordetella_pertussis_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_000_000_002); // classes: fq
        map.insert("bacteria_bordetella_pertussis_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_000_001); // classes: fq
        map.insert("bacteria_bordetella_pertussis_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // classes: fq
        map.insert("bacteria_bordetella_pertussis_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.000_000_000_001); // classes: fq, tet, chl
        map.insert("bacteria_bordetella_pertussis_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_000_005); // classes: fq, mls, tet, chl
        map.insert("bacteria_bordetella_pertussis_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_global_porin_loss_emergence_rate".to_string(), 0.000_000_000_001); // classes: none (currently zeroed)
        map.insert("bacteria_bordetella_pertussis_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_000_000_01); // classes: sulf
        map.insert("bacteria_bordetella_pertussis_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_000_000_1); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_000_001); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_bordetella_pertussis_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_000_000_01); // classes: tet
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls
        map.insert("bacteria_bordetella_pertussis_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_bordetella_pertussis_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.000_000_000_002); // classes: pen, mls, tet, chl
        map.insert("bacteria_bordetella_pertussis_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

           // M. genitalium — Atypical (no cell wall), Fastidious
        // Band 7 (x3.8)
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_cat_emergence_rate".to_string(), 0.003    ); // classes: chl
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.005    ); // classes: mls; lower clindamycin/macrolides without collapsing them completely
        map.insert("bacteria_mycoplasma_genitalium_mechanism_target_site_cfr_emergence_rate".to_string(), 0.005    ); // classes: oxa, mls, chl
        map.insert("bacteria_mycoplasma_genitalium_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.05     ); // classes: fq; raise cipro/ofloxacin directly without using broad efflux
        map.insert("bacteria_mycoplasma_genitalium_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.05     ); // classes: fq; primary FQ lever for levo/moxi in this panel
        map.insert("bacteria_mycoplasma_genitalium_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_05    ); // classes: fq, mls, tet, chl
        map.insert("bacteria_mycoplasma_genitalium_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.001   ); // classes: sulf
        map.insert("bacteria_mycoplasma_genitalium_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.001   ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.01   ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_mycoplasma_genitalium_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_protection_tet_m_emergence_rate".to_string(), 0.02     ); // classes: tet  
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.01    ); // classes: mls; keep macrolides elevated, but closer to the target band
        map.insert("bacteria_mycoplasma_genitalium_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_genitalium_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // M. pneumoniae — Atypical (no cell wall), Fastidious
        // Band 6 (x3)
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_003); // classes: chl
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_015); // classes: mls
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_003); // classes: oxa, mls, chl
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_000_3); // classes: fq
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_15); // classes: fq
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_15); // classes: fq, mls, tet, chl
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_000_003); // classes: sulf
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_000_003); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_03); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_000_3); // classes: tet
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_000_005); // classes: mls; increasingly common
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mycoplasma_pneumoniae_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // ======================================================================
        // Obligate Intracellular and Special Cases
        // ======================================================================
        // C. trachomatis — Obligate intracellular, Fastidious
        // Band 6 (x0.9375)
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_001_3 ); // classes: chl
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.000_000_1  ); // classes: mls
        map.insert("bacteria_chlamydia_trachomatis_mechanism_target_site_cfr_emergence_rate".to_string(), 0.000_000_002   ); // classes: oxa, mls, chl
        map.insert("bacteria_chlamydia_trachomatis_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_000_2  ); // classes: fq
        map.insert("bacteria_chlamydia_trachomatis_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_2  ); // classes: fq
        map.insert("bacteria_chlamydia_trachomatis_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_02   ); // classes: fq, mls, tet, chl
        map.insert("bacteria_chlamydia_trachomatis_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_000_002   ); // classes: sulf
        map.insert("bacteria_chlamydia_trachomatis_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.000_000_002   ); // classes: other (metronidazole, nitrofurantoin, furazolidone)
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_02   ); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_chlamydia_trachomatis_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_000_2  ); // classes: tet
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.000_000_000_1 ); // classes: mls; rare
        map.insert("bacteria_chlamydia_trachomatis_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_chlamydia_trachomatis_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // tier 0

        // T. pallidum — Spirochete
        // Band 7 (x15)
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_15); // classes: chl
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_001_5); // classes: fq
        map.insert("bacteria_treponema_pallidum_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_75); // classes: fq
        map.insert("bacteria_treponema_pallidum_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_15); // classes: fq, mls, tet, chl
        map.insert("bacteria_treponema_pallidum_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_000_15); // classes: sulf
        map.insert("bacteria_treponema_pallidum_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_15); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_treponema_pallidum_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_001_5); // classes: tet
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; disabled here; macrolide resistance can be 23S-mediated
        map.insert("bacteria_treponema_pallidum_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_treponema_pallidum_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 

        // ======================================================================
        // Acid-Fast Bacteria
        // ======================================================================
        // MDR M. tuberculosis — Acid-fast, Mycobacteria
        // Band 8 (x300)
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_esbl_ctx_m_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_esbl_tem_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_esbl_shv_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_ampc_cmy_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_ampc_dha_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_kpc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_ndm_vim_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_oxa_48_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_cat_emergence_rate".to_string(), 0.000_000_000_000); // classes: chl
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_16s_rrmt_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_target_site_pbp2a_meca_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_target_site_van_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_target_site_van_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_target_site_erm_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_target_site_cfr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_mutation_gyra_primary_emergence_rate".to_string(), 0.000_000_000_000); // classes: fq
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_mutation_gyra_parc_secondary_emergence_rate".to_string(), 0.000_000_000_000 ); // classes: fq
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_protection_qnr_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_efflux_acrab_tolc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_efflux_mexxy_oprm_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_global_efflux_pump_emergence_rate".to_string(), 0.000_000_000_000); // classes: fq, mls, tet, chl
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_porin_loss_ompk35_36_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_porin_loss_oprd_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_global_porin_loss_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_modification_mcr_1_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_mutation_folate_pathway_emergence_rate".to_string(), 0.000_000_000_000); // classes: sulf
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_mutation_nitroreductase_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_fos_a_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_mutation_mpr_f_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_mutation_rpo_b_emergence_rate".to_string(), 0.000_000_000_000); // classes: other (rifampicin, fidaxomicin)
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_protection_fus_b_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_protection_tet_m_emergence_rate".to_string(), 0.000_000_000_000); // classes: tet
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_aac_aph_emergence_rate".to_string(), 0.0); // tier 0 — TB uses rrs/rpsL mutations
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_bla_z_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_enzyme_oxa_acinetobacter_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_mutation_23s_rrna_emergence_rate".to_string(), 0.0); // classes: mls; disabled here; linezolid resistance can involve 23S rRNA
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_efflux_tet_abc_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_mutation_pbp_mosaic_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_efflux_mtr_cde_emergence_rate".to_string(), 0.0); // tier 0
        map.insert("bacteria_mdr_mycobacterium_tuberculosis_mechanism_as_yet_unknown_emergence_rate".to_string(), 0.0); // classes: broad placeholder; 



        // Resistance enhancement multipliers: how much each mechanism increases resistance level
        
        // Beta-Lactamases
        // map.insert("resistance_mechanism_enzyme_penicillinase_blaz_enhancement_multiplier".to_string(), 0.9); // Staph Penicillins dead (Mech removed)
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_multiplier".to_string(), 0.8); // High for Cephs
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_multiplier".to_string(), 0.6); // TEM
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_multiplier".to_string(), 0.6); // SHV
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_multiplier".to_string(), 0.7); // Good for Cephs
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_multiplier".to_string(), 0.7); 
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_multiplier".to_string(), 0.95); // Resists all
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_multiplier".to_string(), 0.95); // Resists all
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_multiplier".to_string(), 0.6); // Variable

        // Target Mods
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_multiplier".to_string(), 0.99); // MecA is absolute
        // map.insert("resistance_mechanism_mutation_pbp2x_enhancement_multiplier".to_string(), 0.5); // Pneumo Penicillin is gradual (Mech removed)
        map.insert("resistance_mechanism_mutation_gyra_primary_enhancement_multiplier".to_string(), 0.4); // Low level FQ
        map.insert("resistance_mechanism_mutation_gyra_parc_secondary_enhancement_multiplier".to_string(), 0.95); // High level FQ
        // map.insert("resistance_mechanism_mutation_rpob_enhancement_multiplier".to_string(), 0.95); // Rifampicin R is high (Mech removed)
        map.insert("resistance_mechanism_target_site_cfr_enhancement_multiplier".to_string(), 0.95); // Linezolid/Chloramphenicol

        // Enzymatic/Protection
        map.insert("resistance_mechanism_protection_tet_m_enhancement_multiplier".to_string(), 0.9); // TetM/TetO: high-level tetracycline resistance (16-64× MIC increase)
        // map.insert("resistance_mechanism_enzyme_aac_6_enhancement_multiplier".to_string(), 0.85); (Mech removed/merged)
        map.insert("resistance_mechanism_enzyme_16s_rrmt_enhancement_multiplier".to_string(), 0.95); // High level Aminoglycoside
        map.insert("resistance_mechanism_enzyme_cat_enhancement_multiplier".to_string(), 0.9);
        map.insert("resistance_mechanism_target_site_erm_b_enhancement_multiplier".to_string(), 0.9); // MLSb is high
        // map.insert("resistance_mechanism_efflux_mef_enhancement_multiplier".to_string(), 0.4); // Mef is low/moderate (Mech removed)
        map.insert("resistance_mechanism_protection_qnr_enhancement_multiplier".to_string(), 0.2); // Qnr is low

        // Other
        map.insert("resistance_mechanism_porin_loss_oprd_enhancement_multiplier".to_string(), 0.8);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_multiplier".to_string(), 0.8);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_multiplier".to_string(), 0.3); // Efflux
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_multiplier".to_string(), 0.3); // Efflux
        map.insert("resistance_mechanism_modification_mcr_1_enhancement_multiplier".to_string(), 0.85);
        map.insert("resistance_mechanism_target_site_van_a_enhancement_multiplier".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_van_b_enhancement_multiplier".to_string(), 0.99);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_multiplier".to_string(), 0.2);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_multiplier".to_string(), 0.2);

        // --- New mechanisms ---
        map.insert("resistance_mechanism_mutation_folate_pathway_enhancement_multiplier".to_string(), 0.85); // sul/dfr genes: high-level sulfonamide/trimethoprim resistance
        map.insert("resistance_mechanism_mutation_nitroreductase_enhancement_multiplier".to_string(), 0.7);  // nim/nfsAB: moderate-high resistance to nitroimidazoles/nitrofurans
        map.insert("resistance_mechanism_enzyme_fos_a_enhancement_multiplier".to_string(), 0.8);            // FosA: strong fosfomycin inactivation
        map.insert("resistance_mechanism_mutation_mpr_f_enhancement_multiplier".to_string(), 0.6);          // MprF: moderate daptomycin resistance (often heteroresistance)
        map.insert("resistance_mechanism_mutation_rpo_b_enhancement_multiplier".to_string(), 0.95);         // RpoB: high-level rifampicin/fidaxomicin resistance
        map.insert("resistance_mechanism_protection_fus_b_enhancement_multiplier".to_string(), 0.7);        // FusB: moderate fusidic acid resistance
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_multiplier".to_string(), 0.5);      // PBP mosaic: moderate β-lactam resistance via target modification
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_multiplier".to_string(), 0.4);      // mtrCDE efflux: moderate broad efflux (macrolides, penicillins, tetracyclines, chloramphenicol)
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_multiplier".to_string(), 0.5);      // Calibration placeholder 3: partial resistance

        // --- Additional mechanisms (5 new) ---
        map.insert("resistance_mechanism_enzyme_aac_aph_enhancement_multiplier".to_string(), 0.85);       // AAC/APH/ANT: strong aminoglycoside inactivation
        map.insert("resistance_mechanism_enzyme_bla_z_enhancement_multiplier".to_string(), 0.90);         // blaZ: high-level penicillin hydrolysis in Staph
        map.insert("resistance_mechanism_enzyme_oxa_acinetobacter_enhancement_multiplier".to_string(), 0.80); // OXA-23/40/58: variable carbapenem hydrolysis
        map.insert("resistance_mechanism_mutation_23s_rrna_enhancement_multiplier".to_string(), 0.80);    // 23S rRNA: moderate-high macrolide resistance
        map.insert("resistance_mechanism_efflux_tet_abc_enhancement_multiplier".to_string(), 0.70);       // TetA/B/C: moderate tetracycline efflux

        // --- Per-drug-class enhancement multipliers ---
        // Key format: resistance_mechanism_{mechanism}_enhancement_{drug_class}
        // Where drug_class is one of: pen, bli, c1_2g, c3g, c4_5g, cft_avi, mer_vab, azt_avi, carb, mono, fq, ag, mls, glyc, tet, poly, oxa, chl, sulf, other
        // If a class-specific key is absent, the legacy single value above is used as fallback.
        // Only non-zero cross-class values need to be specified; classes where a mechanism has no effect
        // are handled by the binary mechanism_applicable table (0.0 default).

        // ESBL CTX-M: Strong vs penicillins/C3G, weaker vs C4-5G, overcome by BLIs
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_pen".to_string(), 0.90);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_bli".to_string(), 0.25);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_bli_anti_pseudomonal".to_string(), 0.25);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_bli_sulbactam".to_string(), 0.25);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_bli_anti_pseudomonal".to_string(), 0.25);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_bli_sulbactam".to_string(), 0.25);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_c1_2g".to_string(), 0.90);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_c3g".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_c3g_bli".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_c3g_bli".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_c4g".to_string(), 0.35);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_anti_mrsa_ceph".to_string(), 0.35);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_siderophore_ceph".to_string(), 0.35);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_cft_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_mer_vab".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_azt_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_enhancement_mono".to_string(), 0.80);

        // ESBL TEM
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_pen".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_bli".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_bli_anti_pseudomonal".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_bli_sulbactam".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_bli_anti_pseudomonal".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_bli_sulbactam".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_c1_2g".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_c3g".to_string(), 0.65);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_c3g_bli".to_string(), 0.65);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_c3g_bli".to_string(), 0.65);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_c4g".to_string(), 0.25);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_anti_mrsa_ceph".to_string(), 0.25);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_siderophore_ceph".to_string(), 0.25);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_cft_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_mer_vab".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_azt_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_esbl_tem_enhancement_mono".to_string(), 0.60);

        // ESBL SHV
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_pen".to_string(), 0.80);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_bli".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_bli_anti_pseudomonal".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_bli_sulbactam".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_bli_anti_pseudomonal".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_bli_sulbactam".to_string(), 0.20);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_c1_2g".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_c3g".to_string(), 0.65);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_c3g_bli".to_string(), 0.65);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_c3g_bli".to_string(), 0.65);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_c4g".to_string(), 0.30);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_anti_mrsa_ceph".to_string(), 0.30);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_siderophore_ceph".to_string(), 0.30);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_cft_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_mer_vab".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_azt_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_esbl_shv_enhancement_mono".to_string(), 0.55);

        // AmpC CMY: Resists BLIs (clavulanate-stable), spares carbapenems
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_pen".to_string(), 0.70);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_bli".to_string(), 0.60);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_bli_anti_pseudomonal".to_string(), 0.60);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_bli_sulbactam".to_string(), 0.60);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_bli_anti_pseudomonal".to_string(), 0.60);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_bli_sulbactam".to_string(), 0.60);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_c1_2g".to_string(), 0.80);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_c3g".to_string(), 0.80);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_c3g_bli".to_string(), 0.80);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_c3g_bli".to_string(), 0.80);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_c4g".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_anti_mrsa_ceph".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_siderophore_ceph".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_cft_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_mer_vab".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_azt_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_enhancement_mono".to_string(), 0.10);

        // AmpC DHA
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_pen".to_string(), 0.70);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_bli".to_string(), 0.55);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_bli_anti_pseudomonal".to_string(), 0.55);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_bli_sulbactam".to_string(), 0.55);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_bli_anti_pseudomonal".to_string(), 0.55);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_bli_sulbactam".to_string(), 0.55);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_c1_2g".to_string(), 0.75);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_c3g".to_string(), 0.75);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_c3g_bli".to_string(), 0.75);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_c3g_bli".to_string(), 0.75);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_c4g".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_anti_mrsa_ceph".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_siderophore_ceph".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_cft_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_mer_vab".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_azt_avi".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_enzyme_ampc_dha_enhancement_mono".to_string(), 0.10);

        // KPC: Broad hydrolysis, partially inhibited by avibactam/vaborbactam
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_pen".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_bli".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_bli_anti_pseudomonal".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_bli_sulbactam".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_bli_anti_pseudomonal".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_bli_sulbactam".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_c1_2g".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_c3g".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_c3g_bli".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_c3g_bli".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_c4g".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_anti_mrsa_ceph".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_siderophore_ceph".to_string(), 0.85);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_cft_avi".to_string(), 0.30);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_mer_vab".to_string(), 0.30);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_azt_avi".to_string(), 0.30);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_carb_group1".to_string(), 0.90);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_carb_group2".to_string(), 0.90);
        map.insert("resistance_mechanism_enzyme_kpc_enhancement_mono".to_string(), 0.90);

        // NDM/VIM: Metallo-BL, NOT inhibited by standard BLIs, does NOT hydrolyze aztreonam
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_pen".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_bli".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_bli_anti_pseudomonal".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_bli_sulbactam".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_bli_anti_pseudomonal".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_bli_sulbactam".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_c1_2g".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_c3g".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_c3g_bli".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_c3g_bli".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_c4g".to_string(), 0.90);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_anti_mrsa_ceph".to_string(), 0.90);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_siderophore_ceph".to_string(), 0.90);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_cft_avi".to_string(), 0.50);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_mer_vab".to_string(), 0.50);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_azt_avi".to_string(), 0.50);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_carb_group1".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_carb_group2".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_ndm_vim_enhancement_mono".to_string(), 0.10);

        // OXA-48: Spares cephalosporins, hits carbapenems
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_pen".to_string(), 0.80);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_bli".to_string(), 0.50);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_bli_anti_pseudomonal".to_string(), 0.50);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_bli_sulbactam".to_string(), 0.50);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_bli_anti_pseudomonal".to_string(), 0.50);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_bli_sulbactam".to_string(), 0.50);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_c1_2g".to_string(), 0.40);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_c3g".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_c3g_bli".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_c3g_bli".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_c4g".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_anti_mrsa_ceph".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_siderophore_ceph".to_string(), 0.10);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_cft_avi".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_mer_vab".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_azt_avi".to_string(), 0.15);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_carb_group1".to_string(), 0.70);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_carb_group2".to_string(), 0.70);
        map.insert("resistance_mechanism_enzyme_oxa_48_enhancement_mono".to_string(), 0.0);

        // PBP2a/MecA: Complete beta-lactam resistance except ceftaroline/cefiderocol
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_pen".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_bli".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_bli_anti_pseudomonal".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_bli_sulbactam".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_bli_anti_pseudomonal".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_bli_sulbactam".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_c1_2g".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_c3g".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_c3g_bli".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_c3g_bli".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_c4g".to_string(), 0.70);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_anti_mrsa_ceph".to_string(), 0.70);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_siderophore_ceph".to_string(), 0.70);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_cft_avi".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_mer_vab".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_azt_avi".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_carb_group1".to_string(), 0.85);
        map.insert("resistance_mechanism_target_site_pbp2a_meca_enhancement_carb_group2".to_string(), 0.85);

        // VanA: All glycopeptides resistant
        map.insert("resistance_mechanism_target_site_van_a_enhancement_glyc".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_van_a_enhancement_lipoglycopeptides".to_string(), 0.99);
        map.insert("resistance_mechanism_target_site_van_a_enhancement_lipoglycopeptides".to_string(), 0.99);

        // VanB: Vancomycin resistant, teicoplanin active → blended class value
        map.insert("resistance_mechanism_target_site_van_b_enhancement_glyc".to_string(), 0.70);
        map.insert("resistance_mechanism_target_site_van_b_enhancement_lipoglycopeptides".to_string(), 0.70);
        map.insert("resistance_mechanism_target_site_van_b_enhancement_lipoglycopeptides".to_string(), 0.70);

        // GyrA primary: Low-level FQ resistance
        map.insert("resistance_mechanism_mutation_gyra_primary_enhancement_fq".to_string(), 0.40);

        // GyrA+ParC: High-level FQ resistance
        map.insert("resistance_mechanism_mutation_gyra_parc_secondary_enhancement_fq".to_string(), 0.95);

        // ErmB: MLS resistance
        map.insert("resistance_mechanism_target_site_erm_b_enhancement_mls".to_string(), 0.90);
        map.insert("resistance_mechanism_target_site_erm_b_enhancement_lincosamides".to_string(), 0.90);
        map.insert("resistance_mechanism_target_site_erm_b_enhancement_lincosamides".to_string(), 0.90);

        // Cfr: PhLOPSA phenotype — linezolid, MLS cross-resistance, chloramphenicol
        map.insert("resistance_mechanism_target_site_cfr_enhancement_oxa".to_string(), 0.90);
        map.insert("resistance_mechanism_target_site_cfr_enhancement_mls".to_string(), 0.70);
        map.insert("resistance_mechanism_target_site_cfr_enhancement_lincosamides".to_string(), 0.70);
        map.insert("resistance_mechanism_target_site_cfr_enhancement_lincosamides".to_string(), 0.70);
        map.insert("resistance_mechanism_target_site_cfr_enhancement_chl".to_string(), 0.70);

        // 16S rRMT: High-level aminoglycoside resistance
        map.insert("resistance_mechanism_enzyme_16s_rrmt_enhancement_ag_group1".to_string(), 0.95);
        map.insert("resistance_mechanism_enzyme_16s_rrmt_enhancement_ag_group2".to_string(), 0.95);

        // CAT: Chloramphenicol acetyltransferase
        map.insert("resistance_mechanism_enzyme_cat_enhancement_chl".to_string(), 0.90);

        // Qnr: Low-level FQ protection
        map.insert("resistance_mechanism_protection_qnr_enhancement_fq".to_string(), 0.20);

        // MCR-1: Colistin resistance
        map.insert("resistance_mechanism_modification_mcr_1_enhancement_poly".to_string(), 0.85);

        // Efflux AcrAB-TolC: Restricted to primary substrates (FQ, tetracyclines, chloramphenicol)
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_pen".to_string(), 0.0);  // zeroed: marginal cross-class
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_bli_anti_pseudomonal".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_bli_sulbactam".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_bli_anti_pseudomonal".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_bli_sulbactam".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_c3g".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_c3g_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_c3g_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_c4g".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_anti_mrsa_ceph".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_siderophore_ceph".to_string(), 0.0);  // zeroed: marginal cross-class
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_fq".to_string(), 0.25);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_ag_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_ag_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_mls".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_lincosamides".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_lincosamides".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_tet".to_string(), 0.25);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_glycylcyclines".to_string(), 0.25);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_glycylcyclines".to_string(), 0.25);
        map.insert("resistance_mechanism_efflux_acrab_tolc_enhancement_chl".to_string(), 0.20);

        // Efflux MexXY-OprM: Pseudomonas-specific, primarily aminoglycosides and cefepime
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_c4g".to_string(), 0.20);
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_anti_mrsa_ceph".to_string(), 0.20);
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_siderophore_ceph".to_string(), 0.20);
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_carb_group1".to_string(), 0.05);
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_carb_group2".to_string(), 0.05);
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_fq".to_string(), 0.20);
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_ag_group1".to_string(), 0.30);
        map.insert("resistance_mechanism_efflux_mexxy_oprm_enhancement_ag_group2".to_string(), 0.30);

        // Porin loss OmpK35/36 (Klebsiella): Moderate broad-spectrum, strongest for carbapenems
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_pen".to_string(), 0.30);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_bli".to_string(), 0.40);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_bli_anti_pseudomonal".to_string(), 0.40);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_bli_sulbactam".to_string(), 0.40);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_bli_anti_pseudomonal".to_string(), 0.40);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_bli_sulbactam".to_string(), 0.40);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_c3g".to_string(), 0.40);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_c3g_bli".to_string(), 0.40);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_c3g_bli".to_string(), 0.40);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_c4g".to_string(), 0.30);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_anti_mrsa_ceph".to_string(), 0.30);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_siderophore_ceph".to_string(), 0.30);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_cft_avi".to_string(), 0.25);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_mer_vab".to_string(), 0.25);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_azt_avi".to_string(), 0.25);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_fq".to_string(), 0.0);  // zeroed: marginal cross-class
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_ag_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_enhancement_ag_group2".to_string(), 0.0);

        // Porin loss OprD (Pseudomonas): Primarily carbapenem (especially imipenem)
        map.insert("resistance_mechanism_porin_loss_oprd_enhancement_carb_group1".to_string(), 0.80);
        map.insert("resistance_mechanism_porin_loss_oprd_enhancement_carb_group2".to_string(), 0.80);

        // Global Efflux: Restricted to primary substrates (FQ, macrolides, tetracyclines, chloramphenicol)
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_pen".to_string(), 0.0);  // zeroed: marginal cross-class
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_bli_anti_pseudomonal".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_bli_sulbactam".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_bli_anti_pseudomonal".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_bli_sulbactam".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_c3g".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_c3g_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_c3g_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_c4g".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_anti_mrsa_ceph".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_siderophore_ceph".to_string(), 0.0);  // zeroed: marginal cross-class
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_fq".to_string(), 0.15);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_ag_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_ag_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_mls".to_string(), 0.10);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_lincosamides".to_string(), 0.10);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_lincosamides".to_string(), 0.10);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_tet".to_string(), 0.15);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_glycylcyclines".to_string(), 0.15);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_glycylcyclines".to_string(), 0.15);
        map.insert("resistance_mechanism_global_efflux_pump_enhancement_chl".to_string(), 0.15);

        // Global Porin Loss: zeroed — too broad, undermines drug-class differentiation
        map.insert("resistance_mechanism_global_porin_loss_enhancement_pen".to_string(), 0.0);  // zeroed: marginal cross-class
        map.insert("resistance_mechanism_global_porin_loss_enhancement_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_bli_anti_pseudomonal".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_bli_sulbactam".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_bli_anti_pseudomonal".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_bli_sulbactam".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_c3g".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_c3g_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_c3g_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_c4g".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_anti_mrsa_ceph".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_siderophore_ceph".to_string(), 0.0);  // zeroed: marginal cross-class
        map.insert("resistance_mechanism_global_porin_loss_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_fq".to_string(), 0.0);  // zeroed: marginal cross-class
        map.insert("resistance_mechanism_global_porin_loss_enhancement_ag_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_global_porin_loss_enhancement_ag_group2".to_string(), 0.0);

        // --- Per-drug-class enhancement multipliers for NEW mechanisms ---

        // Folate Pathway (sul/dfr): only affects sulfonamide class
        map.insert("resistance_mechanism_mutation_folate_pathway_enhancement_sulf".to_string(), 0.85);

        // Nitroreductase: affects "other" class drugs (metronidazole, nitrofurantoin, furazolidone)
        map.insert("resistance_mechanism_mutation_nitroreductase_enhancement_other".to_string(), 0.7);

        // FosA: affects "other" class (fosfomycin)
        map.insert("resistance_mechanism_enzyme_fos_a_enhancement_other".to_string(), 0.8);

        // MprF: affects "other" class (daptomycin is classified as Other)
        map.insert("resistance_mechanism_mutation_mpr_f_enhancement_other".to_string(), 0.6);

        // RpoB: affects "other" class (rifampicin, fidaxomicin)
        map.insert("resistance_mechanism_mutation_rpo_b_enhancement_other".to_string(), 0.95);

        // FusB: affects "other" class (fusidic acid)
        map.insert("resistance_mechanism_protection_fus_b_enhancement_other".to_string(), 0.7);

        // MutationPbpMosaic (PBP mosaic): target modification → affects penicillins, cephalosporins, aztreonam
        // NOT carbapenems. BL/BLI combos still affected (BLI irrelevant to target modification).
        // Gradient: penicillins > 1-2G ceph > 3G ceph > 4G ceph > combination agents.
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_pen".to_string(), 0.8);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_bli".to_string(), 0.7);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_bli_anti_pseudomonal".to_string(), 0.7);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_bli_sulbactam".to_string(), 0.7);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_bli_anti_pseudomonal".to_string(), 0.7);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_bli_sulbactam".to_string(), 0.7);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_c1_2g".to_string(), 0.6);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_c3g".to_string(), 0.3);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_c3g_bli".to_string(), 0.3);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_c3g_bli".to_string(), 0.3);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_c4g".to_string(), 0.15);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_anti_mrsa_ceph".to_string(), 0.15);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_siderophore_ceph".to_string(), 0.15);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_cft_avi".to_string(), 0.1);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_mer_vab".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_azt_avi".to_string(), 0.5);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_mono".to_string(), 0.5);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_fq".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_ag_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_ag_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_mls".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_lincosamides".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_lincosamides".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_glyc".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_lipoglycopeptides".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_lipoglycopeptides".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_tet".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_glycylcyclines".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_glycylcyclines".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_poly".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_oxa".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_chl".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_sulf".to_string(), 0.0);
        map.insert("resistance_mechanism_mutation_pbp_mosaic_enhancement_other".to_string(), 0.0);

        // EffluxMtrCde (mtrCDE efflux): broad efflux pump → macrolides, penicillins, tetracyclines, chloramphenicol
        // Moderate enhancement levels — efflux provides partial (not complete) resistance.
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_pen".to_string(), 0.3);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_bli_anti_pseudomonal".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_bli_sulbactam".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_bli_anti_pseudomonal".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_bli_sulbactam".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_c1_2g".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_c3g".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_c3g_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_c3g_bli".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_c4g".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_anti_mrsa_ceph".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_siderophore_ceph".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_cft_avi".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_mer_vab".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_azt_avi".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_carb_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_carb_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_mono".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_fq".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_ag_group1".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_ag_group2".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_mls".to_string(), 0.5);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_lincosamides".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_lincosamides".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_glyc".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_lipoglycopeptides".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_lipoglycopeptides".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_tet".to_string(), 0.4);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_glycylcyclines".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_glycylcyclines".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_poly".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_oxa".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_chl".to_string(), 0.4);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_sulf".to_string(), 0.0);
        map.insert("resistance_mechanism_efflux_mtr_cde_enhancement_other".to_string(), 0.0);

            // AsYetUnknown: placeholder — 0.5 across all classes (dormant, no emergence rates set)
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_pen".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_bli".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_bli_anti_pseudomonal".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_bli_sulbactam".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_bli_anti_pseudomonal".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_bli_sulbactam".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_c1_2g".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_c3g".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_c3g_bli".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_c3g_bli".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_c4g".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_anti_mrsa_ceph".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_siderophore_ceph".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_cft_avi".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_mer_vab".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_azt_avi".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_carb_group1".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_carb_group2".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_mono".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_fq".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_ag_group1".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_ag_group2".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_mls".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_lincosamides".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_lincosamides".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_glyc".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_lipoglycopeptides".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_lipoglycopeptides".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_tet".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_glycylcyclines".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_glycylcyclines".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_poly".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_oxa".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_chl".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_sulf".to_string(), 0.5);
        map.insert("resistance_mechanism_as_yet_unknown_enhancement_other".to_string(), 0.5);

        map.insert("mechanism_assignment_probability_on_any_r_gain".to_string(), 0.8); // Default 80%
        map.insert("mechanism_cache_ewma_decay".to_string(), 0.9); // Default EWMA decay for mechanism cache

        // Mechanism-specific fitness costs (reversion rates per day when drug absent)
        
        map.insert("resistance_mechanism_enzyme_kpc_reversion_rate".to_string(), 0.001); 
        map.insert("resistance_mechanism_enzyme_ndm_vim_reversion_rate".to_string(), 0.0015); // Metallo often on large costly plasmids
        map.insert("resistance_mechanism_enzyme_oxa_48_reversion_rate".to_string(), 0.0005);
        
        map.insert("resistance_mechanism_enzyme_esbl_ctx_m_reversion_rate".to_string(), 0.0006);
        map.insert("resistance_mechanism_enzyme_esbl_tem_reversion_rate".to_string(), 0.0006);
        map.insert("resistance_mechanism_enzyme_esbl_shv_reversion_rate".to_string(), 0.0006);
        map.insert("resistance_mechanism_enzyme_ampc_dha_reversion_rate".to_string(), 0.0006);
        map.insert("resistance_mechanism_enzyme_ampc_cmy_reversion_rate".to_string(), 0.0001); // Native gene, low cost to maintain potential

        map.insert("resistance_mechanism_target_site_pbp2a_meca_reversion_rate".to_string(), 0.0009); // SCCmec is large
        // map.insert("resistance_mechanism_mutation_pbp2x_reversion_rate".to_string(), 0.0002); // Point mutations usually
        
        map.insert("resistance_mechanism_mutation_gyra_primary_reversion_rate".to_string(), 0.0001);
        map.insert("resistance_mechanism_mutation_gyra_parc_secondary_reversion_rate".to_string(), 0.0002);

        // map.insert("resistance_mechanism_mutation_rpob_reversion_rate".to_string(), 0.002); // High cost for rpoB
        map.insert("resistance_mechanism_target_site_erm_b_reversion_rate".to_string(), 0.002); 
        
        map.insert("resistance_mechanism_target_site_van_a_reversion_rate".to_string(), 0.002); 
        map.insert("resistance_mechanism_target_site_van_b_reversion_rate".to_string(), 0.002);

        // map.insert("resistance_mechanism_enzyme_erm_reversion_rate".to_string(), 0.0006); // Use limit above
        map.insert("resistance_mechanism_protection_tet_m_reversion_rate".to_string(), 0.0005); // Moderate fitness cost for Tn916/plasmid carriage
        
        map.insert("resistance_mechanism_modification_mcr_1_reversion_rate".to_string(), 0.0015);
        
        // Defaults for others
        // map.insert("resistance_mechanism_enzyme_penicillinase_blaz_reversion_rate".to_string(), 0.0001);
        // map.insert("resistance_mechanism_enzyme_aac_6_reversion_rate".to_string(), 0.0005);
        // map.insert("resistance_mechanism_enzyme_ant_2_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_enzyme_cat_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_enzyme_16s_rrmt_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_target_site_cfr_reversion_rate".to_string(), 0.0005);
        // map.insert("resistance_mechanism_efflux_mef_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_protection_qnr_reversion_rate".to_string(), 0.0001);
        map.insert("resistance_mechanism_porin_loss_oprd_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_porin_loss_ompk35_36_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_global_efflux_pump_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_global_porin_loss_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_efflux_acrab_tolc_reversion_rate".to_string(), 0.0005);
        map.insert("resistance_mechanism_efflux_mexxy_oprm_reversion_rate".to_string(), 0.0005);

        // --- New mechanism reversion rates ---
        map.insert("resistance_mechanism_mutation_folate_pathway_reversion_rate".to_string(), 0.0001);  // Low cost, often on integrons
        map.insert("resistance_mechanism_mutation_nitroreductase_reversion_rate".to_string(), 0.0003);  // Loss-of-function mutations, moderate fitness cost
        map.insert("resistance_mechanism_enzyme_fos_a_reversion_rate".to_string(), 0.0005);            // Plasmid-mediated, moderate cost
        map.insert("resistance_mechanism_mutation_mpr_f_reversion_rate".to_string(), 0.001);            // Membrane modification has fitness cost
        map.insert("resistance_mechanism_mutation_rpo_b_reversion_rate".to_string(), 0.002);            // High fitness cost for rpoB mutations
        map.insert("resistance_mechanism_protection_fus_b_reversion_rate".to_string(), 0.0005);         // Moderate cost
        map.insert("resistance_mechanism_mutation_pbp_mosaic_reversion_rate".to_string(), 0.001);        // Calibration placeholder 1
        map.insert("resistance_mechanism_efflux_mtr_cde_reversion_rate".to_string(), 0.001);        // Calibration placeholder 2
        map.insert("resistance_mechanism_as_yet_unknown_reversion_rate".to_string(), 0.001);        // Calibration placeholder 3

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

        // Syndrome-specific bacteria growth multipliers
        // Syndromes: 1=UTI, 2=Skin/soft tissue, 3=Respiratory, 4=Bloodstream, 5=Intra-abdominal,
        //           6=CNS, 7=GI, 8=Genital, 9=Bone/joint, 10=Other
        map.insert("syndrome_1_bacteria_growth_multiplier".to_string(), 1.0);   // UTI - baseline growth
        map.insert("syndrome_2_bacteria_growth_multiplier".to_string(), 1.1);   // Skin/soft tissue - faster in necrotizing infections
        map.insert("syndrome_3_bacteria_growth_multiplier".to_string(), 1.2);   // Respiratory - rapid progression in pneumonia
        map.insert("syndrome_4_bacteria_growth_multiplier".to_string(), 1.4);   // Bloodstream - fulminant bacteremia/sepsis
        map.insert("syndrome_5_bacteria_growth_multiplier".to_string(), 1.15);  // Intra-abdominal - abscess formation but somewhat contained
        map.insert("syndrome_6_bacteria_growth_multiplier".to_string(), 1.3);   // CNS - meningitis progresses rapidly
        map.insert("syndrome_7_bacteria_growth_multiplier".to_string(), 1.1);   // GI - moderate progression
        map.insert("syndrome_8_bacteria_growth_multiplier".to_string(), 0.9);   // Genital - often chronic/indolent (chlamydia, gonorrhea)
        map.insert("syndrome_9_bacteria_growth_multiplier".to_string(), 0.85);  // Bone/joint - slow progression in osteomyelitis
        map.insert("syndrome_10_bacteria_growth_multiplier".to_string(), 1.0);  // Other - baseline

        // Empiric drug scoring tables (clinician-facing heuristics per syndrome ID)
        // These preserve pre-refactor prescribing patterns when organism is unknown.
        let empiric_syndrome_templates: &[(usize, &[(&str, f64)])] = &[
            // 1 = UTI / Genitourinary
            (
                1,
                &[
                    ("trim_sulf", 11.0),
                    ("amoxicillin_clavulanate", 14.0),
                    ("amoxicillin", 12.0),
                    ("ciprofloxacin", 8.0),
                    ("ampicillin", 10.0),
                    ("levofloxacin", 6.0),
                    ("nitrofurantoin", 5.0),
                    ("cephalexin", 8.0),
                    ("ceftriaxone", 8.0),
                    ("cefazolin", 7.0),
                    ("cefuroxime", 7.0),
                    ("piperacillin_tazobactam", 5.0),
                    ("cefepime", 4.0),
                    ("ceftazidime", 4.0),
                    ("meropenem", 4.0),
                    ("imipenem_c", 4.0),
                    ("ertapenem", 4.0),
                    ("meropenem_vaborbactam", 3.0),
                    ("ceftazidime_avibactam", 3.0),
                    ("aztreonam_avibactam", 3.0),
                    ("cefixime", 7.0),
                    ("colistin", 0.2),
                    ("vancomycin", 0.1),
                    ("linezolid", 0.1),
                ],
            ),
            // 2 = Skin / soft tissue
            (
                2,
                &[
                    ("penicillin_g", 16.0),
                    ("ampicillin", 13.0),
                    ("amoxicillin", 14.0),
                    ("amoxicillin_clavulanate", 14.0),
                    ("cephalexin", 13.0),
                    ("cefazolin", 12.0),
                    ("flucloxacillin", 15.0),
                    ("clindamycin", 12.0),
                    ("trim_sulf", 9.0),
                    ("doxycycline", 9.0),
                    ("minocycline", 9.0),
                    ("linezolid", 10.0),
                    ("tedizolid", 9.0),
                    ("dalbavancin", 9.0),
                    ("vancomycin", 11.0),
                    ("quinu_dalfo", 8.0),
                    ("rifampicin", 0.5),
                    ("ciprofloxacin", 4.0),
                    ("piperacillin_tazobactam", 3.0),
                ],
            ),
            // 3 = Respiratory
            (
                3,
                &[
                    ("amoxicillin_clavulanate", 20.0),
                    ("amoxicillin", 17.0),
                    ("penicillin_g", 16.0),
                    ("ampicillin", 15.0),
                    ("azithromycin", 12.0),
                    ("clarithromycin", 11.0),
                    ("ceftriaxone", 9.5),
                    ("erythromycin", 9.0),
                    ("cefuroxime", 8.5),
                    ("piperacillin_tazobactam", 8.0),
                    ("levofloxacin", 8.0),
                    ("moxifloxacin", 8.0),
                    ("cefixime", 6.5),
                    ("aztreonam_avibactam", 6.0),
                    ("cefepime", 7.5),
                    ("cephalexin", 7.0),
                    ("doxycycline", 6.5),
                    ("vancomycin", 6.5),
                    ("meropenem", 6.0),
                    ("imipenem_c", 6.0),
                    ("ofloxacin", 6.0),
                    ("linezolid", 7.0),
                    ("minocycline", 5.5),
                ],
            ),
            // 4 = Bloodstream / bacteremia
            (
                4,
                &[
                    ("piperacillin_tazobactam", 18.0), // BOOSTED for BL/BLI class share (was 14.0)
                    ("meropenem", 13.0),
                    ("imipenem_c", 13.0),
                    ("meropenem_vaborbactam", 13.0),
                    ("ceftazidime_avibactam", 12.5),
                    ("aztreonam_avibactam", 12.0),
                    ("cefepime", 12.0),
                    ("ceftazidime", 11.0),
                    ("ceftriaxone", 10.0),
                    ("ampicillin_sulbactam", 16.0), // BOOSTED for BL/BLI class share (was 11.5)
                    ("amoxicillin_clavulanate", 16.0), // BOOSTED for BL/BLI class share (was 10.5)
                    ("ampicillin", 10.0),
                    ("amoxicillin", 9.5),
                    ("penicillin_g", 6.5),
                    ("flucloxacillin", 7.5),
                    ("vancomycin", 11.0),
                    ("linezolid", 10.0),
                    ("tedizolid", 9.0),
                    ("dalbavancin", 8.0),
                    ("quinu_dalfo", 8.5),
                    ("gentamicin", 1.0), // REDUCED: aminoglycoside overuse correction (was 7.0)
                    ("tobramycin", 1.0), // REDUCED: aminoglycoside overuse correction (was 6.5)
                    ("amikacin", 1.0), // REDUCED: aminoglycoside overuse correction (was 7.0)
                    ("colistin", 0.1), // REDUCED: polymyxin overuse correction (was 6.0)
                    ("cefazolin", 6.0),
                    ("ciprofloxacin", 6.0),
                    ("levofloxacin", 5.5),
                    ("cephalexin", 4.0),
                    ("rifampicin", 0.5),
                ],
            ),
            // 5 = Intra-abdominal
            (
                5,
                &[
                    ("metronidazole", 2.5),
                    ("piperacillin_tazobactam", 13.0),
                    ("ampicillin_sulbactam", 12.5),
                    ("amoxicillin_clavulanate", 11.5),
                    ("meropenem", 13.0),
                    ("imipenem_c", 12.5),
                    ("ertapenem", 11.0),
                    ("ceftazidime", 9.0),
                    ("cefepime", 9.0),
                    ("ceftriaxone", 9.0),
                    ("ceftazidime_avibactam", 10.0),
                    ("aztreonam_avibactam", 9.5),
                    ("meropenem_vaborbactam", 10.0),
                    ("ciprofloxacin", 7.0),
                    ("levofloxacin", 6.5),
                    ("ampicillin", 8.0),
                    ("amoxicillin", 7.0),
                    ("trim_sulf", 4.0),
                    ("colistin", 0.1), // REDUCED: polymyxin overuse correction (was 3.5)
                ],
            ),
            // 6 = Central nervous system
            (
                6,
                &[
                    ("ceftriaxone", 15.0),
                    ("ceftazidime", 12.0),
                    ("cefepime", 12.0),
                    ("penicillin_g", 11.0),
                    ("ampicillin", 13.0),
                    ("vancomycin", 13.0),
                    ("linezolid", 10.0),
                    ("cefixime", 1.0),
                    ("meropenem", 11.0),
                    ("imipenem_c", 10.0),
                    ("chloramphenicol", 2.0),
                    ("rifampicin", 1.0),
                    ("piperacillin_tazobactam", 6.0),
                ],
            ),
            // 7 = Gastrointestinal (non-invasive)
            (
                7,
                &[
                    ("ciprofloxacin", 8.0),
                    ("azithromycin", 12.0),
                    ("amoxicillin_clavulanate", 11.0),
                    ("amoxicillin", 10.0),
                    ("ampicillin", 10.0),
                    ("levofloxacin", 6.0),
                    ("ampicillin_sulbactam", 9.0),
                    ("trim_sulf", 8.5),
                    ("doxycycline", 8.5),
                    ("minocycline", 6.5),
                    ("cefixime", 4.5),
                    ("penicillin_g", 5.0),
                    ("cephalexin", 5.0),
                    ("cefuroxime", 5.0),
                    ("furazolidone", 0.2),
                    ("metronidazole", 0.2),
                    ("rifampicin", 0.5),
                ],
            ),
            // 8 = Genital / pelvic
            (
                8,
                &[
                    ("penicillin_g", 14.0),
                    ("azithromycin", 13.0),
                    ("ceftriaxone", 13.0),
                    ("cefixime", 10.5),
                    ("doxycycline", 12.0),
                    ("amoxicillin_clavulanate", 12.0),
                    ("amoxicillin", 11.0),
                    ("cefuroxime", 10.0),
                    ("clindamycin", 9.0),
                    ("ampicillin", 9.0),
                    ("ampicillin_sulbactam", 8.0),
                    ("ciprofloxacin", 4.0),
                    ("levofloxacin", 4.0),
                    ("cephalexin", 6.0),
                    ("trim_sulf", 5.0),
                    ("metronidazole", 0.25),
                    ("rifampicin", 0.5),
                ],
            ),
            // 9 = Bone / joint / hardware-associated
            (
                9,
                &[
                    ("penicillin_g", 14.0),
                    ("ampicillin", 12.0),
                    ("cefazolin", 13.0),
                    ("cephalexin", 11.0),
                    ("ceftriaxone", 11.0),
                    ("flucloxacillin", 14.0),
                    ("vancomycin", 12.0),
                    ("linezolid", 11.0),
                    ("tedizolid", 10.0),
                    ("dalbavancin", 10.0),
                    ("clindamycin", 10.0),
                    ("ciprofloxacin", 9.0),
                    ("levofloxacin", 9.0),
                    ("trim_sulf", 8.0),
                    ("meropenem", 7.0),
                    ("piperacillin_tazobactam", 6.5),
                    ("rifampicin", 2.0),
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
                    ("aztreonam_avibactam", 7.5),
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

    // Hospitalization Parameters - Logistic Model
    // P(hospitalization) = 1 / (1 + exp(-log_odds))
    // log_odds = base + (age_years × per_year) + sepsis_effect
    // This naturally bounds P ∈ (0,1) without clamping
    map.insert("hospitalization_base_log_odds".to_string(), -10.4); // baseline: log(0.00003/0.99997) ≈ -10.4
    map.insert("hospitalization_log_odds_per_age_year".to_string(), 0.02); // ~2% increase in log-odds per year
    map.insert("hospitalization_log_odds_sepsis".to_string(), 4.4); // log(80) ≈ 4.4, sepsis strongly increases admission
    map.insert("hospitalization_log_odds_symptomatic_infection".to_string(), 2.5); // ~12x multiplier: severe symptomatic infection (level > threshold) drives pre-antibiotic hospitalization
    map.insert("hospitalization_symptomatic_infection_level_threshold".to_string(), 3.0); // Bacteria level must exceed this AND have symptoms to trigger hospitalization effect
    
    // Regional hospitalization log-odds adjustments - HICs admit sepsis patients more readily
    // Positive values = higher admission rates (better access to hospital care)
    // This improves survival in HICs by getting patients into hospitals faster
    map.insert("north_america_hospitalization_log_odds".to_string(), 0.5); // Good hospital access
    map.insert("europe_hospitalization_log_odds".to_string(), 0.6); // Universal healthcare, excellent access
    map.insert("oceania_hospitalization_log_odds".to_string(), 0.4); // Good access in developed areas
    map.insert("asia_hospitalization_log_odds".to_string(), 0.0); // Reference - mixed access
    map.insert("south_america_hospitalization_log_odds".to_string(), -0.2); // Variable access
    map.insert("africa_hospitalization_log_odds".to_string(), -0.5); // Limited hospital capacity
    map.insert("hospitalization_recovery_rate_per_day".to_string(), 0.28); // Slightly shorter stays (~3.6 day avg) to reinforce target occupancy
    map.insert("hospitalization_max_days".to_string(), 30.0); // Max days in hospital before forced discharge (as fallback)
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

        // Regional antibiotic initiation log-odds adjustments
        // These affect probability of receiving ANY antibiotic when symptomatic
        // Reflects healthcare access disparities across regions
        // Values are additive adjustments to base log-odds; negative = lower access
        map.insert("north_america_antibiotic_initiation_log_odds".to_string(), 0.0); // Reference region
        map.insert("europe_antibiotic_initiation_log_odds".to_string(), 0.0); // Similar access to NA
        map.insert("oceania_antibiotic_initiation_log_odds".to_string(), 0.0); // Similar access
        map.insert("asia_antibiotic_initiation_log_odds".to_string(), -0.5); // ~38% reduction in odds (P: 26% -> ~17%)
        map.insert("south_america_antibiotic_initiation_log_odds".to_string(), -0.8); // ~55% reduction (P: 26% -> ~13%)
        map.insert("africa_antibiotic_initiation_log_odds".to_string(), -1.4); // ~75% reduction (P: 26% -> ~7%)

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
        map.insert("regional_resistance_penalty_very_high".to_string(), 0.3); // Gentler to promote penicillins (was 0.2)
        map.insert("regional_resistance_penalty_high".to_string(), 0.5); // Gentler to promote penicillins (was 0.4)
        map.insert("regional_resistance_penalty_moderate".to_string(), 0.8); // Gentler to promote penicillins (was 0.7)
        map.insert("regional_resistance_threshold_very_high".to_string(), 0.6); // Threshold for very high resistance (60%)
        map.insert("regional_resistance_threshold_high".to_string(), 0.45); // Threshold for high resistance (45%)
        map.insert("regional_resistance_threshold_moderate".to_string(), 0.1); // Threshold for moderate resistance (10%)

        // Drug Spectrum Classifications (1.0=narrow, 5.0=very broad)
        map.insert("drug_colistin_spectrum_breadth".to_string(), 4.0); // Broad spectrum (mainly Gram-negative)
        map.insert("drug_penicillin_g_spectrum_breadth".to_string(), 2.0); // Narrow spectrum
        map.insert("drug_flucloxacillin_spectrum_breadth".to_string(), 1.6); // Very narrow anti-staphylococcal penicillin
        map.insert("drug_amoxicillin_spectrum_breadth".to_string(), 3.0); // Medium spectrum
        map.insert("drug_azithromycin_spectrum_breadth".to_string(), 4.0); // Broad spectrum
        map.insert("drug_ciprofloxacin_spectrum_breadth".to_string(), 4.5); // Very broad spectrum
        map.insert("drug_trim_sulf_spectrum_breadth".to_string(), 3.5); // Medium-broad spectrum
        map.insert("drug_meropenem_spectrum_breadth".to_string(), 5.0); // Very broad spectrum (carbapenem)
        map.insert("drug_cefepime_spectrum_breadth".to_string(), 4.0); // Broad spectrum (4th gen cephalosporin)
        map.insert("drug_vancomycin_spectrum_breadth".to_string(), 2.5); // Narrow-medium spectrum (gram-positive only)
        map.insert("drug_linezolid_spectrum_breadth".to_string(), 2.0); // Narrow spectrum (gram-positive only)
        map.insert("drug_ceftriaxone_spectrum_breadth".to_string(), 4.0); // Broad spectrum (3rd gen cephalosporin)
        map.insert("drug_cefixime_spectrum_breadth".to_string(), 2.8); // 3rd gen oral


        // NEW: Logistic Sepsis Risk Parameters (replacing old linear model)
        map.insert("sepsis_baseline_log_odds".to_string(), -14.0); // Fallback baseline for organisms without explicit intercept


        // sepsis rates
        // ***sepsis_incidence
        // Bacteria-specific sepsis baseline log-odds (best-guess placeholders calibrated by clinical severity)
        let bacteria_sepsis_baseline_overrides: &[(&str, f64)] = &[
            ("acinetobacter_baumannii", -7.7),
            ("citrobacter_spp.", -9.2),
            ("enterobacter_spp.", -7.7),
            ("enterococcus_faecalis", -7.5),       // was -8.0; CFR 0.74× under-death
            ("enterococcus_faecium", -7.0),         // was -8.0; CFR 0.40× under-death
            ("escherichia_coli", -9.5),             // was -8.5; CFR 2.7× over-death
            ("klebsiella_pneumoniae", -7.5),        // was -8.0; CFR 0.83× under-death
            ("morganella_spp.", -7.8),
            ("proteus_spp.", -7.8),
            ("serratia_spp.", -8.0),
            ("pseudomonas_aeruginosa", -6.5),       // was -8.0; CFR 0.38× under-death
            ("stenotrophomonas_maltophilia", -8.0),
            ("staphylococcus_aureus", -7.3),
            ("staphylococcus_epidermidis", -8.0),
            ("streptococcus_pneumoniae", -10.5),
            ("salmonella_enterica_serovar_typhi", -8.0),
            ("salmonella_enterica_serovar_paratyphi_a", -9.0), // was -8.0; CFR 3.4× over-death
            ("invasive_non-typhoidal_salmonella_spp.", -9.2),  // was -8.0; CFR 1.9× over-death
            ("shigella_spp.", -12.0),
            ("neisseria_gonorrhoeae", -23.0),
            ("streptococcus_pyogenes", -6.5),       // was -8.0; CFR 0.20× under-death (invasive GAS)
            ("streptococcus_agalactiae", -7.0),     // was -8.0; CFR 0.55× under-death
            ("haemophilus_influenzae", -9.2),
            ("chlamydia_trachomatis", -19.0),
            ("vibrio_cholerae", -9.0),
            ("neisseria_meningitidis", -8.6),
            ("listeria_monocytogenes", -8.0),
            ("clostridioides_difficile", -11.0),
            ("campylobacter_jejuni", -20.0),
            ("enterobacter_cloacae", -7.8),
            ("yersinia_enterocolitica", -9.5),     // was -8.0; CFR 6.7× over-death
            ("moraxella_catarrhalis", -10.8),
            ("treponema_pallidum", -11.0),
            ("bordetella_pertussis", -11.0),
            ("helicobacter_pylori", -250.0),
            ("mdr_mycobacterium_tuberculosis", -38.0),
        ];

        for (bacteria, log_odds) in bacteria_sepsis_baseline_overrides {
        map.insert(format!("{}_sepsis_baseline_log_odds", bacteria), *log_odds);
        }
        // ***sepsis_incidence
        map.insert("log_odds_sepsis_infection_level".to_string(), 0.93); // 0.8      Log odds increase per unit bacterial level
        // === [I] Clinical outcome scalars (mortality, sepsis, toxicity) ===
        // Collects mortality/sepsis odds adjustments together so scenario designers can reason about
        // outcome severity in one place. These parameters shape the probability of severe outcomes
        // once infection is established.
        map.insert("log_odds_sepsis_infection_duration".to_string(), 0.005); // Log odds increase per day of infection duration (increased from 0.001)

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

        // PER-REGION SEPSIS ONSET LOG-ODDS (replaces 2-tier system for better differentiation)
        // Negative values = lower risk of developing sepsis from infection (better early recognition/treatment)
        // Values calibrated to reduce overall deaths while increasing regional differentiation
        map.insert("log_odds_sepsis_onset_region_north_america".to_string(), -0.5); // Excellent early sepsis recognition
        map.insert("log_odds_sepsis_onset_region_europe".to_string(), -0.6); // Best early warning systems
        map.insert("log_odds_sepsis_onset_region_oceania".to_string(), -0.5); // Good healthcare infrastructure
        map.insert("log_odds_sepsis_onset_region_asia".to_string(), -0.1); // Mixed - improving rapidly
        map.insert("log_odds_sepsis_onset_region_south_america".to_string(), 0.0); // Variable access
        map.insert("log_odds_sepsis_onset_region_africa".to_string(), 0.1); // Limited early detection capacity
        
        // DEPRECATED: Kept for backward compatibility, no longer used in rules/mod.rs
        map.insert("log_odds_sepsis_region_a".to_string(), -0.5); // Higher resource region - better sepsis recognition/treatment
        map.insert("log_odds_sepsis_region_b".to_string(), 0.0);  // Lower resource region - delayed recognition/limited resources


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


        //  Immunosuppression seeding, onset, and recovery rates
        map.insert("immunosuppression_startup_seed_fraction".to_string(), 0.05);      // Initial fraction seeded into the composite higher-risk state

        // Temporary immunodeficiency represents short-lived higher-risk states.
        map.insert("temporary_immunosuppression_onset_rate_per_day".to_string(), 0.00005);  // Calibrated background onset into the temporary state
        map.insert("temporary_immunosuppression_recovery_rate_per_day".to_string(), 0.01);   // Faster recovery than chronic immunosuppression

        // Chronic immunodeficiency represents longer-lived persistent higher-risk states.
        map.insert("chronic_immunosuppression_onset_rate_per_day".to_string(), 0.00006);     // Background onset into the chronic state
        map.insert("chronic_immunosuppression_recovery_rate_per_day".to_string(), 0.0012);   // Slower recovery keeps chronic episodes persistent

        // Age effect on chronic-vs-temporary typing when a new immunodeficiency episode occurs.
        // These are structural assignment weights, not literal prevalence estimates.
        map.insert("chronic_immunodeficiency_probability_age_0_1".to_string(), 0.3);   // Early-life episodes can map to persistent congenital or neonatal high-risk states
        map.insert("chronic_immunodeficiency_probability_age_1_18".to_string(), 0.2);  // Most childhood episodes remain temporary
        map.insert("chronic_immunodeficiency_probability_age_18_65".to_string(), 0.4); // Adult episodes more often persist as chronic higher-risk states
        map.insert("chronic_immunodeficiency_probability_age_65_plus".to_string(), 0.6); // Late-life episodes more often map to persistent frailty/immunosenescence-type states

        // Prophylactic antibiotic use in immunocompromised patients
        map.insert("antibiotic_infection_prevention_efficacy".to_string(), 0.7);       // 70% efficacy: allow more breakthrough infections despite prophylaxis


        // Sepsis onset additional factors (log-odds scale)
        map.insert("log_odds_sepsis_onset_immunosuppressed".to_string(), 0.7); // ~2x higher onset risk for immunocompromised
        map.insert("log_odds_sepsis_onset_hospitalized".to_string(), 0.5); // ~1.6x higher onset risk when hospitalized (sicker patients)
        map.insert("log_odds_sepsis_onset_not_under_care".to_string(), 1.0); // ~2.7x higher onset risk if not receiving treatment

        // Sepsis death logistic model parameters (log-odds scale)
        // The logistic model: P(death) = 1 / (1 + exp(-log_odds))
        // where log_odds = base + age_effect + region_effect + immuno_effect + level_effect + duration_effect + care_effect
        // ***sepsis_death
        map.insert("sepsis_death_base_log_odds".to_string(), -5.0); // REVISED: Lowered from -4.5 to balance population death rates (~0.4% daily base)
        map.insert("sepsis_death_log_odds_age_infant".to_string(), 1.1); // Infants: +1.1 log-odds (~3x baseline)
        map.insert("sepsis_death_log_odds_age_child".to_string(), -0.7); // Children: -0.7 log-odds (~0.5x baseline)
        map.insert("sepsis_death_log_odds_age_adult".to_string(), 0.0); // Adults: reference category
        map.insert("sepsis_death_log_odds_age_elderly".to_string(), 0.9); // Elderly: +0.9 log-odds (~2.5x baseline)
        map.insert("sepsis_death_log_odds_immunosuppressed".to_string(), 1.5); // Immunosuppressed: +1.5 log-odds (~4.5x)
        map.insert("sepsis_death_log_odds_bacteria_level".to_string(), 0.35); // Per unit bacteria level (0-5 scale)
        map.insert("sepsis_death_log_odds_duration".to_string(), 0.04); // Per day of sepsis after early phase
        map.insert("sepsis_death_log_odds_early_phase".to_string(), 0.8); // Additional risk in first 72h (~2.2x)
        map.insert("sepsis_death_early_phase_days".to_string(), 3.0); // Early phase duration (days)
        map.insert("sepsis_death_log_odds_not_under_care".to_string(), 1.4); // Not receiving treatment: +1.4 log-odds (~4x)

        // Region-specific sepsis mortality multipliers (reflecting healthcare quality)
        // Regional sepsis mortality multipliers - CALIBRATED for realistic HIC/LMIC differentiation
        // Wider spread for differentiation but overall lower values to reduce total deaths
        // HIC regions: 0.4-0.5 (better ICU care, early intervention)
        // LMIC regions: 0.8-1.5 (reduced from previous 1.2-2.0)
        map.insert("north_america_sepsis_mortality_multiplier".to_string(), 0.5); // Was 0.8 - better ICU survival
        map.insert("europe_sepsis_mortality_multiplier".to_string(), 0.4); // Was 0.7 - best sepsis outcomes globally
        map.insert("oceania_sepsis_mortality_multiplier".to_string(), 0.5); // Was 0.8 - similar to NA
        map.insert("asia_sepsis_mortality_multiplier".to_string(), 0.9); // Was 1.2 - improving rapidly
        map.insert("south_america_sepsis_mortality_multiplier".to_string(), 1.1); // Was 1.4 - variable but improving
        map.insert("africa_sepsis_mortality_multiplier".to_string(), 1.5); // Was 2.0 - still highest but reduced 

        // Sepsis Recovery Parameters (Logistic Model)
        map.insert("sepsis_base_log_odds_of_recovery_per_day".to_string(), -0.0); // -0.5 Base log odds (low baseline recovery probability ~12%)
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
        // ^^^ micro
        // ══▶▶ CALIBRATION AXIS 6: microbiome seeding probability — CHANGE HERE ◀◀══
        map.insert("microbiome_resistance_multiplier_on_acquisition".to_string(), 0.50);  //  0.50 
        map.insert("infection_from_microbiome_dampening".to_string(), 0.70);  // 0.85  ***
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
        map.insert("antibiotic_clearance_log_odds_per_unit_activity".to_string(), 0.5);  // ***
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
        // ══▶▶ CALIBRATION AXIS 4: microbiome→infection bridge — CHANGE HERE ◀◀══
        map.insert("carrier_resistance_inheritance_probability".to_string(), 0.50);  // 0.50            
        // ══▶▶ CALIBRATION AXIS 5: fraction from human reservoir — CHANGE HERE ◀◀══
        map.insert("community_resistance_dilution_factor".to_string(), 0.50); // 0.50

        // ══▶▶ CALIBRATION AXIS 1: resistance persistence / reversion speed — CHANGE HERE ◀◀══
        map.insert("mechanism_reversion_rate_global_multiplier".to_string(), 1.0);
        // ══▶▶ CALIBRATION AXIS 2: de-novo emergence in infections — CHANGE HERE ◀◀══
        map.insert("infection_de_novo_multiplier".to_string(), 1.0);
        // ══▶▶ CALIBRATION AXIS 3: de-novo emergence in gut carriage — CHANGE HERE ◀◀══
        map.insert("microbiome_de_novo_multiplier".to_string(), 1.0);
        // ══▶▶ CALIBRATION AXIS 7: horizontal gene transfer rate scaling — CHANGE HERE ◀◀══
        map.insert("hgt_multiplier".to_string(), 1.0);

 //     map.insert(
 //         "majority_r_memory_retention_per_day".to_string(),
 //         0.93,
 //     );  // note this is not currently implemented

        // Majority_r cache defaults: rolling window horizon and minimum sample threshold.
        map.insert("majority_r_window_days".to_string(), 100.0);
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
        
        // Drug toxicity death logistic model parameters (log-odds scale)
        // The logistic model: P(death) = 1 / (1 + exp(-log_odds))
        // where log_odds = base + reservoir_effect + age_effect + immuno_effect + hospital_effect
        map.insert("toxicity_death_base_log_odds".to_string(), -8.0); // Base log-odds (~0.03% with no toxicity reservoir)
        map.insert("toxicity_death_log_odds_per_reservoir_unit".to_string(), 2.0); // Per unit toxicity reservoir (~7x per unit)
        map.insert("toxicity_death_log_odds_age_infant".to_string(), 0.6); // Infants: +0.6 log-odds (~1.8x baseline)
        map.insert("toxicity_death_log_odds_age_child".to_string(), 0.2); // Children: +0.2 log-odds (~1.2x baseline)
        map.insert("toxicity_death_log_odds_age_adult".to_string(), 0.0); // Adults: reference category
        map.insert("toxicity_death_log_odds_age_elderly".to_string(), 0.8); // Elderly: +0.8 log-odds (~2.2x baseline)
        map.insert("toxicity_death_log_odds_immunosuppressed".to_string(), 0.9); // Immunosuppressed: +0.9 log-odds (~2.5x)
        map.insert("toxicity_death_log_odds_hospitalized".to_string(), 0.25); // Hospitalized: +0.25 log-odds (~1.3x)

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
        map.insert("escherichia_coli_log_odds_infant".to_string(), 0.6);             // Elevated risk - anatomic factors
        map.insert("escherichia_coli_log_odds_preschool".to_string(), -0.1);         // Slightly above baseline
        map.insert("escherichia_coli_log_odds_school".to_string(), -0.2);            // Near-baseline risk
        map.insert("escherichia_coli_log_odds_young_adult".to_string(), 0.3);        // Moderate risk - sexual activity, pregnancy
        map.insert("escherichia_coli_log_odds_middle_age".to_string(), 0.2);         // Mildly elevated risk
        map.insert("escherichia_coli_log_odds_elderly".to_string(), 0.9);            // Higher risk - multiple factors

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
                "penicillin_g" | "ampicillin" | "amoxicillin" => 1.0,
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
                "chloramphenicol" => 0.8,
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
}

// --- String Parameters (for template names, etc.) ---
lazy_static! {
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
}

lazy_static! {
    pub static ref PARAMETER_STORE: ParameterStore = ParameterStore::from_parameter_map(&PARAMETERS);
}

lazy_static! {
    static ref ACTIVE_PARAMETER_CONTEXT: AtomicPtr<ActiveParameterContext> =
        AtomicPtr::new(ptr::null_mut());
}

// ---------------- 12) Helper lookups (drug intro, availability, etc.) ----------------
/// Helper accessor for the indexed parameter store.
#[allow(dead_code)]
pub fn parameter_store() -> &'static ParameterStore {
    if let Some(context) = get_active_parameter_context() {
        &context.store
    } else {
        &PARAMETER_STORE
    }
}

/// Retrieves a global simulation parameter.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_global_param(key: &str) -> Option<f64> {
    parameter_map().get(key).copied()
}

/// Retrieves a bacteria-specific simulation parameter.
/// It directly looks up "{bacteria_name}_{param_suffix}".
/// Because all bacteria now have explicit entries, there's no need for a "generic_bacteria_" fallback in this function.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_bacteria_param(bacteria_name: &str, param_suffix: &str) -> Option<f64> {
    let canonical = canonicalize_bacteria_slug(bacteria_name);
    let specific_key = format!("{}_{}", canonical.as_ref(), param_suffix);
    parameter_map().get(&specific_key).copied()
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
        if let Some(value) = parameter_map().get(key) {
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
    parameter_map().get(&specific_key).copied()
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
    parameter_map().get(&availability_key).copied().unwrap_or(1.0) // Default to available if not specified
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
            // ESBL resistance affects penicillins + 1st-3rd gen cephalosporins (BL/BLI combinations overcome ESBL, cefepime spared)
            vec!["penicillin_g", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "ceftriaxone", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Fluoroquinolone resistance (gyrA/parC mutations affect all FQs)
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            // Aminoglycoside resistance (AMEs and 16S rRNA methylases)
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // acinetobacter_baumannii resistance patterns
        m.insert("acinetobacter_baumannii", vec![
            // β-lactamase affects most β-lactams (BL/BLI combinations included)
            vec!["penicillin_g", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Carbapenemase affects carbapenems (including BL/BLI)
            vec!["meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            // Aminoglycoside resistance
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // klebsiella_pneumoniae resistance patterns
        m.insert("klebsiella_pneumoniae", vec![
            // ESBL resistance (BL/BLI combinations overcome ESBL)
            vec!["penicillin_g", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "ceftriaxone", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Carbapenemase (KPC, NDM, etc.)
            vec!["meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"],
            // Fluoroquinolone resistance (gyrA/parC mutations)
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // streptococcus_pneumoniae resistance patterns
        m.insert("streptococcus_pneumoniae", vec![
            // Macrolide resistance (erm genes affect all macrolides)
            vec!["erythromycin", "azithromycin", "clarithromycin"],
            // Penicillin resistance (PBP alterations affect all penicillins including BL/BLI)
            vec!["penicillin_g", "ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
        ]);

        // staphylococcus_aureus resistance patterns
        m.insert("staphylococcus_aureus", vec![
            // β-lactamase affects penicillins
            vec!["penicillin_g", "ampicillin", "amoxicillin"],
            // MRSA affects most β-lactams
            vec!["cephalexin", "cefazolin", "cefuroxime", "ceftriaxone"],
            // Macrolide-lincosamide resistance
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
        ]);

        m.insert("staphylococcus_epidermidis", vec![
            // mecA-mediated resistance impacts nearly all β-lactams
            vec!["penicillin_g", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "ceftriaxone"],
            // Macrolide/clindamycin cross-resistance common in CoNS
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
            // Multidrug efflux impacting fluoroquinolones when acquired
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
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
            // β-lactamase affects multiple β-lactams including BL/BLI combinations
            vec!["piperacillin", "piperacillin_tazobactam", "ceftazidime", "ceftazidime_avibactam", "cefepime"],
            // Carbapenemase
            vec!["meropenem", "meropenem_vaborbactam", "imipenem_c"],
            // Fluoroquinolone resistance (gyrA/gyrB mutations)
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            // Aminoglycoside resistance
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // Enterobacter species resistance patterns
        m.insert("enterobacter_spp.", vec![
            // AmpC β-lactamase (chromosomal) - affects BL/BLI combinations as well
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "cephalexin", "cefazolin", "cefuroxime"],
            // ESBL if acquired
            vec!["ceftriaxone", "cefixime", "ceftazidime", "cefepime"],
            // Fluoroquinolone resistance (gyrA/parC mutations)
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
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

        // --- ENTEROCOCCI ---
        // Enterococcus faecalis - typically ampicillin-susceptible, VRE rare
        m.insert("enterococcus_faecalis", vec![
            // Penicillin resistance (PBP alterations) - affects all penicillins
            vec!["penicillin_g", "ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "piperacillin", "piperacillin_tazobactam", "ticarcillin", "ticarcillin_clavulanate"],
            // Macrolide-lincosamide resistance (erm genes) - MLSB phenotype
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
            // Fluoroquinolone resistance (gyrA/parC mutations)
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            // Tetracycline resistance (tet genes)
            vec!["tetracycline", "doxycycline", "minocycline"],
            // Glycopeptide resistance (VanA/VanB) - vancomycin + teicoplanin
            vec!["vancomycin", "teicoplanin", "dalbavancin"],
        ]);

        // Enterococcus faecium - typically ampicillin-resistant, VRE more common
        m.insert("enterococcus_faecium", vec![
            // Penicillin resistance - most E. faecium intrinsically resistant
            vec!["penicillin_g", "ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "piperacillin", "piperacillin_tazobactam"],
            // Macrolide-lincosamide resistance (erm genes)
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            // Tetracycline resistance
            vec!["tetracycline", "doxycycline", "minocycline"],
            // Glycopeptide resistance (VanA/VanB)
            vec!["vancomycin", "teicoplanin", "dalbavancin"],
        ]);

        // --- OTHER ENTEROBACTERALES ---
        // Citrobacter spp. - AmpC producers like Enterobacter
        m.insert("citrobacter_spp.", vec![
            // AmpC β-lactamase (chromosomal, inducible) - affects BL/BLI combinations
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "cephalexin", "cefazolin", "cefuroxime"],
            // Extended-spectrum cephalosporin resistance
            vec!["ceftriaxone", "cefixime", "ceftazidime", "cefepime"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            // Aminoglycoside resistance
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // Enterobacter cloacae - same as Enterobacter spp.
        m.insert("enterobacter_cloacae", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "cephalexin", "cefazolin", "cefuroxime"],
            vec!["ceftriaxone", "cefixime", "ceftazidime", "cefepime"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // Morganella spp. - intrinsic AmpC, chromosomal β-lactamase
        m.insert("morganella_spp.", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "cephalexin", "cefazolin", "cefuroxime"],
            vec!["ceftriaxone", "cefixime", "ceftazidime", "cefepime"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // Proteus spp. - intrinsic ampicillin resistance varies
        m.insert("proteus_spp.", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "cephalexin", "cefazolin", "cefuroxime"],
            vec!["ceftriaxone", "cefixime", "ceftazidime"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // Serratia spp. - chromosomal AmpC
        m.insert("serratia_spp.", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "cephalexin", "cefazolin", "cefuroxime"],
            vec!["ceftriaxone", "cefixime", "ceftazidime", "cefepime"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // Providencia stuartii
        m.insert("p_stuartii", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "cephalexin", "cefazolin"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // --- SALMONELLA & ENTERIC PATHOGENS ---
        // Salmonella Typhi
        m.insert("salmonella_enterica_serovar_typhi", vec![
            // Ampicillin resistance + early-gen cephalosporins (MDR strains) - BL/BLI combinations included
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate", "cephalexin", "cefazolin", "cefuroxime"],
            // Fluoroquinolone resistance (gyrA mutations)
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            // Cephalosporin resistance (ESBL - affects 3rd gen)
            vec!["ceftriaxone", "cefixime", "ceftazidime"],
        ]);

        // Salmonella Paratyphi A - similar to Typhi
        m.insert("salmonella_enterica_serovar_paratyphi_a", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["ceftriaxone", "cefixime", "ceftazidime"],
        ]);

        // Invasive non-typhoidal Salmonella
        m.insert("invasive_non-typhoidal_salmonella_spp.", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["ceftriaxone", "cefixime", "ceftazidime"],
        ]);

        // Shigella spp. - high FQ resistance in many regions
        m.insert("shigella_spp.", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["ceftriaxone", "cefixime", "ceftazidime"],
            vec!["tetracycline", "doxycycline"],
        ]);

        // Vibrio cholerae
        m.insert("vibrio_cholerae", vec![
            vec!["tetracycline", "doxycycline"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["erythromycin", "azithromycin", "clarithromycin"],
        ]);

        // Campylobacter jejuni - macrolide and FQ resistance
        m.insert("campylobacter_jejuni", vec![
            vec!["erythromycin", "azithromycin", "clarithromycin"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["tetracycline", "doxycycline"],
        ]);

        // Yersinia enterocolitica
        m.insert("yersinia_enterocolitica", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["tetracycline", "doxycycline"],
        ]);

        // Helicobacter pylori - clarithromycin and metronidazole resistance key
        m.insert("helicobacter_pylori", vec![
            vec!["amoxicillin", "ampicillin", "amoxicillin_clavulanate", "ampicillin_sulbactam"],
            vec!["clarithromycin", "erythromycin", "azithromycin"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["metronidazole"],
            vec!["tetracycline", "doxycycline"],
        ]);

        // --- STREPTOCOCCI ---
        // Streptococcus pyogenes (Group A Strep) - penicillin still universally susceptible
        m.insert("streptococcus_pyogenes", vec![
            // Macrolide resistance (erm genes)
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
            // Tetracycline resistance
            vec!["tetracycline", "doxycycline"],
            // Fluoroquinolone resistance (rare)
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // Streptococcus agalactiae (Group B Strep)
        m.insert("streptococcus_agalactiae", vec![
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
            vec!["tetracycline", "doxycycline"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // --- FASTIDIOUS ORGANISMS ---
        // Haemophilus influenzae - beta-lactamase common (TEM-1)
        // BL/BLI combinations partially overcome beta-lactamase but not always reliable
        m.insert("haemophilus_influenzae", vec![
            vec!["ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
            vec!["erythromycin", "azithromycin", "clarithromycin"],
            vec!["ciprofloxacin", "levofloxacin"],
            vec!["tetracycline", "doxycycline"],
        ]);

        // Moraxella catarrhalis - nearly all beta-lactamase positive (BRO-1/BRO-2)
        // Beta-lactamase affects both simple penicillins and BL/BLI combinations
        m.insert("moraxella_catarrhalis", vec![
            vec!["ampicillin", "amoxicillin", "penicillin_g", "piperacillin", "ticarcillin",
                 "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            vec!["erythromycin", "azithromycin", "clarithromycin"],
        ]);

        // Neisseria gonorrhoeae - critical resistance concern
        // Beta-lactamase (TEM-1, TEM-135) affects most penicillins including some BL/BLI
        m.insert("neisseria_gonorrhoeae", vec![
            vec!["penicillin_g", "ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["tetracycline", "doxycycline"],
            vec!["erythromycin", "azithromycin"],
            vec!["ceftriaxone", "cefixime"],
        ]);

        // Neisseria meningitidis - reduced penicillin susceptibility spreading
        m.insert("neisseria_meningitidis", vec![
            vec!["penicillin_g", "ampicillin", "amoxicillin", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["rifampicin"],
        ]);

        // --- ANAEROBES & OTHER ---
        // Clostridioides difficile - limited treatment options
        m.insert("clostridioides_difficile", vec![
            vec!["vancomycin"],
            vec!["metronidazole"],
        ]);

        // Bacteroides fragilis
        m.insert("bacteroides_fragilis", vec![
            vec!["metronidazole"],
            vec!["clindamycin"],
            vec!["meropenem", "imipenem_c"],
        ]);

        // Listeria monocytogenes - intrinsic cephalosporin resistance
        // Penicillin resistance affects BL/BLI combinations as well
        m.insert("listeria_monocytogenes", vec![
            vec!["ampicillin", "amoxicillin", "penicillin_g", "ampicillin_sulbactam", "amoxicillin_clavulanate"],
            vec!["tetracycline", "doxycycline"],
        ]);

        // --- ATYPICALS ---
        // Chlamydia trachomatis - macrolides/tetracyclines
        m.insert("chlamydia_trachomatis", vec![
            vec!["azithromycin", "erythromycin", "clarithromycin"],
            vec!["doxycycline", "tetracycline"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        ]);

        // Mycoplasma genitalium - macrolide resistance emerging
        m.insert("mycoplasma_genitalium", vec![
            vec!["azithromycin", "erythromycin", "clarithromycin"],
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
            vec!["doxycycline", "tetracycline"],
        ]);

        // Treponema pallidum - limited options, macrolide resistance emerging
        m.insert("treponema_pallidum", vec![
            vec!["erythromycin", "azithromycin", "clarithromycin"],
            vec!["doxycycline", "tetracycline"],
        ]);

        // Bordetella pertussis
        m.insert("bordetella_pertussis", vec![
            vec!["erythromycin", "azithromycin", "clarithromycin"],
        ]);

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
pub fn get_age_dependent_bacteria_sepsis_risk_log_odds(bacteria_name: &str, age_days: u32) -> f64 {
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
        map.insert("penicillin_g", 3555);     // 3555 // 1942 (12 years after 1930)
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
        map.insert("flucloxacillin", 14600); // 1970
        map.insert("cefuroxime", 17525);     // 1978 (48 years after 1930)
        map.insert("ceftriaxone", 19715);    // 1984 (54 years after 1930)
        map.insert("cefixime", 21535); // ~1989
        map.insert("ceftazidime", 20080);    // 1985 (55 years after 1930)
        map.insert("cefepime", 24195);       // 1996 (66 years after 1930)
        map.insert("ceftaroline", 29305);    // 2010 (80 years after 1930)

        // Carbapenems
        map.insert("meropenem", 24195);      // 1996 (66 years after 1930)
        map.insert("imipenem_c", 20080);     // 1985 (55 years after 1930)
        map.insert("ertapenem", 25920);      // 2001 (71 years after 1930)

        // Monobactams
        map.insert("aztreonam", 20445);      // 1986 (56 years after 1930)
        map.insert("aztreonam_avibactam", 34675); // ~2025

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
        map.insert("chloramphenicol", 6935);    // 1949 (19 years after 1930) - chloramphenicol
        map.insert("nitrofurantoin", 8395);  // 1953 (23 years after 1930)
        map.insert("retapamulin", 28405);    // 2007 (77 years after 1930) - topical antibiotic
        map.insert("fusidic_a", 11680);       // 1962 (32 years after 1930) - fusidic acid
        map.insert("metronidazole", 10965);   // 1960 (30 years after 1930)
        map.insert("furazolidone", 9125);    // 1955 (25 years after 1930)
        
        // Newer generation antibiotics
        map.insert("tigecycline", 28040);        // 2007 (77 years after 1930)
        map.insert("daptomycin", 27375);         // 2005 (75 years after 1930)
        map.insert("fosfomycin", 10590);         // 1959 (29 years after 1930)
        map.insert("fidaxomicin", 29565);        // 2011 (81 years after 1930)
        map.insert("ceftolozane_tazobactam", 30295); // 2014 (84 years after 1930)
        map.insert("cefiderocol", 33510);        // 2019 (89 years after 1930)

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
            // Sample uniformly within the age band
            let age = rng.gen_range(age_min..age_max);
            return (region, age);
        }
    }

    // Fallback (should rarely be reached)
    (Region::Asia, 0)
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
// ============================================================================
// === Resistance Floor Helper Functions ===
// ============================================================================
// These functions support the resistance floor feature for rare bacteria where
// cache-based sampling doesn't sustain observed resistance levels.

/// Get the drug class name for a given drug
/// Returns the class name used in resistance floor parameters
pub fn get_drug_class(drug: &str) -> Option<&'static str> {
    match drug {
        // Penicillins (including BL/BLI)
        "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin" |
        "amoxicillin_clavulanate" | "piperacillin_tazobactam" | "ampicillin_sulbactam" | 
        "ticarcillin_clavulanate" => Some("penicillins"),
        
        // Cephalosporins 1st/2nd gen
        "cephalexin" | "cefazolin" | "cefuroxime" => Some("cephalosporins_1_2"),
        
        // Cephalosporins 3rd/4th gen (including BL/BLI)
        "ceftriaxone" | "ceftazidime" | "cefixime" | "cefepime" | "ceftaroline" | 
        "ceftazidime_avibactam" => Some("cephalosporins_3_4"),
        
        // Carbapenems (including BL/BLI)
        "meropenem" | "imipenem_c" | "ertapenem" | "meropenem_vaborbactam" => Some("carbapenems"),
        
        // Monobactams - no separate floor, treat like cephalosporins 3/4 for coverage
        "aztreonam" => Some("cephalosporins_3_4"),
            "aztreonam_avibactam" => Some("cephalosporins_3_4"),
        
        // Macrolides
        "erythromycin" | "azithromycin" | "clarithromycin" => Some("macrolides"),
        
        // Lincosamides - treat like macrolides (MLSb resistance)
        "clindamycin" => Some("macrolides"),
        
        // Aminoglycosides
        "gentamicin" | "tobramycin" | "amikacin" => Some("aminoglycosides"),
        
        // Fluoroquinolones
        "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin" => Some("fluoroquinolones"),
        
        // Tetracyclines
        "tetracycline" | "doxycycline" | "minocycline" => Some("tetracyclines"),
        
        // Glycopeptides
        "vancomycin" | "teicoplanin" | "dalbavancin" => Some("glycopeptides"),
        
        // Oxazolidinones
        "linezolid" | "tedizolid" => Some("oxazolidinones"),
        
        // Folate antagonists
        "trim_sulf" => Some("folate_antagonists"),
        
        // Polymyxins
        "colistin" => Some("polymyxins"),
        
        // Sulfanilamide - original sulfonamide, treat as folate antagonist
        "sulfanilamide" => Some("folate_antagonists"),
        
        // Others without specific floors
        _ => None,
    }
}

/// Get the introduction day for a drug (uses existing DRUG_INTRODUCTION_DATES)
/// Returns None if drug introduction is not configured
pub fn get_drug_introduction_day(drug: &str) -> Option<i32> {
    // Use the existing DRUG_INTRODUCTION_DATES via get_drug_introduction_time_step
    get_drug_introduction_time_step(drug).map(|ts| ts as i32)
}

/// Get the earliest introduction day for any drug in a class
/// This is used to determine when resistance floors should start ramping
pub fn get_drug_class_introduction_day(drug_class: &str) -> Option<i32> {
    // Map drug class to its constituent drugs and find earliest introduction
    let drugs: &[&str] = match drug_class {
        "penicillins" => &["penicillin_g", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin",
                          "amoxicillin_clavulanate", "piperacillin_tazobactam", "ampicillin_sulbactam", 
                          "ticarcillin_clavulanate"],
        "cephalosporins_1_2" => &["cephalexin", "cefazolin", "cefuroxime"],
        "cephalosporins_3_4" => &["ceftriaxone", "ceftazidime", "cefepime", "ceftaroline", "cefixime", 
                                  "ceftazidime_avibactam", "aztreonam"],
        "carbapenems" => &["meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"],
        "macrolides" => &["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
        "aminoglycosides" => &["gentamicin", "tobramycin", "amikacin"],
        "fluoroquinolones" => &["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
        "tetracyclines" => &["tetracycline", "doxycycline", "minocycline"],
        "glycopeptides" => &["vancomycin", "teicoplanin", "dalbavancin"],
        "oxazolidinones" => &["linezolid", "tedizolid"],
        "folate_antagonists" => &["trim_sulf", "sulfanilamide"],
        "polymyxins" => &["colistin"],
        _ => return None,
    };
    
    drugs.iter()
        .filter_map(|drug| get_drug_introduction_day(drug))
        .min()
}

/// Check if resistance floors are enabled globally
pub fn resistance_floors_enabled() -> bool {
    get_global_param("resistance_floor_feature_enabled").unwrap_or(0.0) > 0.5
}

/// Check if resistance floors are enabled for a specific bacteria
pub fn bacteria_resistance_floor_enabled(bacteria_name: &str) -> bool {
    if !resistance_floors_enabled() {
        return false;
    }
    let canonical = canonicalize_bacteria_slug(bacteria_name);
    let key = format!("bacteria_{}_resistance_floor_enabled", canonical.as_ref());
    get_global_param(&key).unwrap_or(0.0) > 0.5
}

/// Get the resistance floor ramp period for a bacteria (in years)
pub fn get_resistance_floor_ramp_years(bacteria_name: &str) -> f64 {
    let canonical = canonicalize_bacteria_slug(bacteria_name);
    let key = format!("bacteria_{}_resistance_floor_ramp_years", canonical.as_ref());
    get_global_param(&key).unwrap_or(10.0) // Default 10 year ramp
}

/// Get the target resistance floor for a bacteria-drug combination
/// Returns 0.0 if no floor is configured
pub fn get_resistance_floor_target(bacteria_name: &str, drug: &str) -> f64 {
    let drug_class = match get_drug_class(drug) {
        Some(class) => class,
        None => return 0.0,
    };
    
    let canonical = canonicalize_bacteria_slug(bacteria_name);
    let key = format!("bacteria_{}_{}_resistance_floor", canonical.as_ref(), drug_class);
    get_global_param(&key).unwrap_or(0.0)
}

/// Calculate the effective resistance floor for a bacteria-drug pair at a given simulation day
/// 
/// The floor ramps linearly from 0 at drug class introduction to the target floor
/// over the configured ramp period.
/// 
/// Returns 0.0 if:
/// - Resistance floors are disabled
/// - The bacteria doesn't have floors enabled  
/// - The simulation day is before the drug class was introduced
/// - No floor is configured for this bacteria-drug combination
pub fn calculate_resistance_floor(bacteria_name: &str, drug: &str, current_day: i32) -> f64 {
    // Check if floors are enabled for this bacteria
    if !bacteria_resistance_floor_enabled(bacteria_name) {
        return 0.0;
    }
    
    // Get drug class
    let drug_class = match get_drug_class(drug) {
        Some(class) => class,
        None => return 0.0,
    };
    
    // Get drug class introduction day
    let intro_day = match get_drug_class_introduction_day(drug_class) {
        Some(day) => day,
        None => return 0.0,
    };
    
    // If before drug introduction, no floor
    if current_day < intro_day {
        return 0.0;
    }
    
    // Get target floor
    let target_floor = get_resistance_floor_target(bacteria_name, drug);
    if target_floor <= 0.0 {
        return 0.0;
    }
    
    // Calculate ramp fraction
    let ramp_years = get_resistance_floor_ramp_years(bacteria_name);
    let ramp_days = (ramp_years * 365.0) as i32;
    let days_since_intro = current_day - intro_day;
    
    let ramp_fraction = if ramp_days <= 0 {
        1.0
    } else {
        (days_since_intro as f64 / ramp_days as f64).min(1.0)
    };
    
    target_floor * ramp_fraction
}
