/// Generates the Appendix B markdown for MODEL_DESCRIPTION.md.
///
/// Prints structured, thematically organised parameter tables derived from the
/// live Rust configuration as resolved Markdown tables.
use amr_project::config::{
    get_drug_class, get_drug_introduction_time_step, BacteriumMechanismStatus, PARAMETERS,
    PARAMETER_STORE,
};
use amr_project::simulation::population::{
    bacterium_has_separate_microbiome_compartment, DrugClass, Region, ResistanceMechanism,
    AGE_CATEGORY_SEQUENCE, BACTERIA_LIST, DRUG_SHORT_NAMES,
};

const REGION_NAMES: [&str; 6] = [
    "north_america",
    "south_america",
    "africa",
    "asia",
    "europe",
    "oceania",
];

const REGION_VARIANTS: [Region; 7] = [
    Region::NorthAmerica,
    Region::SouthAmerica,
    Region::Africa,
    Region::Asia,
    Region::Europe,
    Region::Oceania,
    Region::Home,
];

const SYNDROME_NAMES: [&str; 11] = [
    "none",
    "uti",
    "skin_soft_tissue",
    "respiratory",
    "bloodstream",
    "intra_abdominal",
    "cns_meningitis",
    "gastrointestinal",
    "genital_sti",
    "bone_joint",
    "other",
];

const VACCINES: [&str; 4] = ["pneumococcal", "meningococcal", "hib", "pertussis"];

fn format_value(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let abs = v.abs();
    if abs >= 0.001 && abs < 1_000_000.0 {
        let s = format!("{:.10}", v);
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        // Scientific notation: trim trailing zeros in the coefficient only
        let s = format!("{:.6e}", v);
        if let Some(e_pos) = s.find('e') {
            let (coeff, exp) = s.split_at(e_pos);
            let coeff = coeff.trim_end_matches('0').trim_end_matches('.');
            format!("{}{}", coeff, exp)
        } else {
            s
        }
    }
}

fn bool_value(value: bool) -> f64 {
    if value {
        1.0
    } else {
        0.0
    }
}

fn mechanism_status_label(status: BacteriumMechanismStatus) -> &'static str {
    match status {
        BacteriumMechanismStatus::ExcludedHost => "excluded host",
        BacteriumMechanismStatus::EligibleNoDeNovo => "eligible; no de novo or HGT",
        BacteriumMechanismStatus::HgtOnly => "eligible; HGT only",
        BacteriumMechanismStatus::DeNovo => "eligible; de novo enabled",
    }
}

fn configured_value(key: &str, fallback: f64) -> f64 {
    PARAMETERS.get(key).copied().unwrap_or(fallback)
}

fn configured_rows_matching<F>(predicate: F) -> Vec<Vec<String>>
where
    F: Fn(&str) -> bool,
{
    let mut entries: Vec<(&String, &f64)> = PARAMETERS
        .iter()
        .filter(|(key, _)| predicate(key))
        .collect();
    entries.sort_by(|a, b| a.0.cmp(b.0));
    entries
        .into_iter()
        .map(|(key, value)| vec![key.clone(), format_value(*value)])
        .collect()
}

/// Print a Markdown table from a header row and data rows.
fn md_table(headers: &[&str], rows: &[Vec<String>]) {
    // Header
    print!("|");
    for h in headers {
        print!(" {} |", h);
    }
    println!();
    // Separator — right-align numeric columns (all except the first)
    print!("|");
    for (i, _) in headers.iter().enumerate() {
        if i == 0 {
            print!(" --- |");
        } else {
            print!(" ---: |");
        }
    }
    println!();
    // Data rows
    for row in rows {
        print!("|");
        for cell in row {
            print!(" {} |", cell);
        }
        println!();
    }
    println!();
}

fn main() {
    let store = &*PARAMETER_STORE;
    let _params = &*PARAMETERS;

    print_heading();
    print_global_scalars(store);
    print_drug_properties(store);
    print_bacteria_properties(store);
    print_drug_bacteria_matrix(store);
    print_regional_parameters(store);
    print_age_dependent_parameters(store);
    print_syndrome_parameters(store);
    print_clearance_parameters(store);
    print_immunodeficiency_sex_vaccination(store);
    print_resistance_mechanisms(store);
    print_hgt_matrix(store);
}

fn print_heading() {
    println!("## Appendix B — Parameter Reference");
    println!();
    println!(
        "This appendix is auto-generated from the live Rust configuration. \
              Parameters are organised thematically into resolved tables \
              derived from the internal data structures. All values shown are \
              the effective defaults before any run-level pathway sensitivity \
              multipliers are applied. Where a family has a uniform fallback, the fallback \
              is stated and only explicit exceptions are listed. Dynamically parsed era \
              overrides and environmental floors are included. Raw compatibility keys that \
              are loaded nowhere in the executable rules are intentionally excluded."
    );
    println!();
}

// ─────────────────────────────────────────────────────────────────────────────
// B.1  Global Scalar Parameters
// ─────────────────────────────────────────────────────────────────────────────

fn print_global_scalars(store: &amr_project::config::ParameterStore) {
    let g = &store.globals;
    println!("### B.1 Global Scalar Parameters");
    println!();
    println!(
        "Scalar parameters that govern cross-cutting model behaviour. \
              Grouped thematically; each row gives the parameter name and its \
              default value."
    );
    println!();
    println!(
        "See: \
              [§6.1 Treatment initiation](#61-treatment-initiation-deciding-to-start-antibiotics), \
              [§6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use), \
              [§6.3 Drug pharmacokinetics](#63-drug-pharmacokinetics), \
              [§6.7 Drug toxicity](#67-drug-toxicity), \
              [§2.4 Hospitalisation](#24-hospitalisation), \
              [§2.5 Travel](#25-travel), \
              [§4.3 Sepsis](#43-sepsis), \
              [§7.3 Resistance emergence](#73-resistance-emergence), \
              [§7.4 Resistance reversion](#74-resistance-reversion-and-fitness-costs), \
              [§8 Microbiome and Carriage](#8-microbiome-and-carriage), \
              [§9 Horizontal Gene Transfer](#9-horizontal-gene-transfer-hgt), \
              [§10 Mortality](#10-mortality)."
    );
    println!();

    print_scalar_group(
        "Treatment Initiation (logistic model)",
        &[
            (
                "antibiotic_initiation_base_log_odds",
                g.antibiotic_initiation_base_log_odds,
            ),
            (
                "antibiotic_initiation_log_odds_symptomatic_infection",
                g.antibiotic_initiation_log_odds_symptomatic_infection,
            ),
            (
                "antibiotic_initiation_log_odds_test_identified",
                g.antibiotic_initiation_log_odds_test_identified,
            ),
            (
                "antibiotic_initiation_log_odds_already_on_drug",
                g.antibiotic_initiation_log_odds_already_on_drug,
            ),
            (
                "antibiotic_initiation_log_odds_immunodeficiency",
                g.antibiotic_initiation_log_odds_immunodeficiency,
            ),
            (
                "antibiotic_initiation_log_odds_sepsis",
                g.antibiotic_initiation_log_odds_sepsis,
            ),
            (
                "antibiotic_initiation_log_odds_hospitalized",
                g.antibiotic_initiation_log_odds_hospitalized,
            ),
            (
                "antibiotic_initiation_log_odds_no_indication",
                g.antibiotic_initiation_log_odds_no_indication,
            ),
        ],
    );

    print_scalar_group(
        "Drug Activity and Cessation",
        &[
            (
                "drug_activity_to_bacteria_level_multiplier",
                g.drug_activity_to_bacteria_level_multiplier,
            ),
            (
                "drug_activity_slow_clearance_probability",
                g.drug_activity_slow_clearance_probability,
            ),
            (
                "drug_activity_slow_clearance_multiplier",
                g.drug_activity_slow_clearance_multiplier,
            ),
            (
                "double_dose_probability_if_identified_infection",
                g.double_dose_probability_if_identified_infection,
            ),
            (
                "random_drug_cessation_probability",
                g.random_drug_cessation_probability,
            ),
            (
                "random_drug_cessation_probability_if_no_active_infection",
                g.random_drug_cessation_probability_if_no_active_infection,
            ),
            (
                "antibiotic_infection_prevention_efficacy",
                g.antibiotic_infection_prevention_efficacy,
            ),
        ],
    );

    print_scalar_group(
        "Drug Selection",
        &[
            (
                "minimal_potency_threshold_for_drug_selection",
                g.minimal_potency_threshold_for_drug_selection,
            ),
            ("drug_selection_temperature", g.drug_selection_temperature),
            ("reserve_drug_score_penalty", g.reserve_drug_score_penalty),
        ],
    );

    print_scalar_group(
        "Treatment Failure and Restart",
        &[
            (
                "treatment_failure_enabled",
                bool_value(g.treatment_failure_enabled),
            ),
            (
                "treatment_failure_assessment_day",
                g.treatment_failure_assessment_day as f64,
            ),
            ("treatment_failure_threshold", g.treatment_failure_threshold),
            (
                "drug_failure_memory_days",
                g.drug_failure_memory_days as f64,
            ),
            (
                "restart_window_enabled",
                bool_value(g.restart_window_enabled),
            ),
            ("restart_window_days", g.restart_window_days as f64),
            (
                "restart_bacteria_level_threshold",
                g.restart_bacteria_level_threshold,
            ),
            ("restart_window_probability", g.restart_window_probability),
            (
                "drug_evaluation_days_post_infection",
                g.drug_evaluation_days_post_infection as f64,
            ),
        ],
    );

    print_scalar_group(
        "Hospitalisation",
        &[
            (
                "hospitalization_base_log_odds",
                g.hospitalization_base_log_odds,
            ),
            (
                "hospitalization_log_odds_per_age_year",
                g.hospitalization_log_odds_per_age_year,
            ),
            (
                "hospitalization_log_odds_sepsis",
                g.hospitalization_log_odds_sepsis,
            ),
            (
                "hospitalization_log_odds_symptomatic_infection",
                g.hospitalization_log_odds_symptomatic_infection,
            ),
            (
                "hospitalization_log_odds_serious_resistance_test_positive",
                g.hospitalization_log_odds_serious_resistance_test_positive,
            ),
            (
                "hospitalization_symptomatic_infection_level_threshold",
                g.hospitalization_symptomatic_infection_level_threshold,
            ),
            (
                "hospital_recovery_rate_per_day",
                g.hospital_recovery_rate_per_day,
            ),
            ("hospital_max_days", g.hospital_max_days),
            (
                "hospital_prevent_discharge_with_sepsis",
                g.hospital_prevent_discharge_with_sepsis,
            ),
        ],
    );

    print_scalar_group(
        "Resistance Emergence and Decay",
        &[
            ("max_resistance_level", g.max_resistance_level),
            (
                "resistance_emergence_bacteria_level_multiplier",
                g.resistance_emergence_bacteria_level_multiplier,
            ),
            (
                "multi_drug_penalty_threshold_num_drugs",
                g.multi_drug_penalty_threshold_num_drugs,
            ),
            (
                "resistance_development_inhibition_single_drug",
                g.resistance_development_inhibition_single_drug,
            ),
            (
                "resistance_development_inhibition_partial_cross",
                g.resistance_development_inhibition_partial_cross,
            ),
            (
                "mechanism_assignment_probability_on_any_r_gain",
                g.mechanism_assignment_probability_on_any_r_gain,
            ),
            (
                "community_profile_cache_retention",
                g.community_profile_cache_retention,
            ),
            (
                "hospital_profile_cache_retention",
                g.hospital_profile_cache_retention,
            ),
            (
                "local_mechanism_persistence_enabled",
                bool_value(g.local_mechanism_persistence_enabled),
            ),
            (
                "local_mechanism_persistence_virtual_profile_mass",
                g.local_mechanism_persistence_virtual_profile_mass,
            ),
            (
                "local_mechanism_persistence_max_sampling_probability",
                g.local_mechanism_persistence_max_sampling_probability,
            ),
            (
                "debug_seed_hospital_cache_resistant_profiles",
                bool_value(g.debug_seed_hospital_cache_resistant_profiles),
            ),
        ],
    );

    print_scalar_group(
        "Microbiome Dynamics",
        &[
            (
                "microbiome_resistance_transfer_probability_per_day",
                g.microbiome_resistance_transfer_probability_per_day,
            ),
            (
                "antibiotic_disruption_decay_half_life_days",
                g.antibiotic_disruption_decay_half_life_days,
            ),
            (
                "infection_from_microbiome_dampening",
                g.infection_from_microbiome_dampening,
            ),
            (
                "carriage_duration_log_odds_coefficient",
                g.carriage_duration_log_odds_coefficient,
            ),
            (
                "carriage_duration_max_log_odds_effect",
                g.carriage_duration_max_log_odds_effect,
            ),
            (
                "antibiotic_clearance_log_odds_per_unit_activity",
                g.antibiotic_clearance_log_odds_per_unit_activity,
            ),
            (
                "carrier_resistance_inheritance_probability",
                g.carrier_resistance_inheritance_probability,
            ),
            (
                "community_resistance_dilution_factor",
                g.community_resistance_dilution_factor,
            ),
        ],
    );

    print_scalar_group(
        "Horizontal Gene Transfer Modifiers",
        &[
            ("hgt_hospital_multiplier", g.hgt_hospital_multiplier),
            (
                "hgt_antibiotic_pressure_multiplier",
                g.hgt_antibiotic_pressure_multiplier,
            ),
            ("hgt_coinfection_multiplier", g.hgt_coinfection_multiplier),
            ("hgt_microbiome_only_penalty", g.hgt_microbiome_only_penalty),
            (
                "hgt_gut_compartment_multiplier",
                g.hgt_gut_compartment_multiplier,
            ),
            (
                "hgt_minority_donor_multiplier",
                g.hgt_minority_donor_multiplier,
            ),
        ],
    );

    print_scalar_group(
        "Travel",
        &[("travel_probability_per_day", g.travel_probability_per_day)],
    );

    print_scalar_group(
        "Bacteria Growth Age Multipliers",
        &[
            (
                "bacteria_growth_age_multiplier_infant",
                g.bacteria_growth_age_multiplier_infant,
            ),
            (
                "bacteria_growth_age_multiplier_child",
                g.bacteria_growth_age_multiplier_child,
            ),
            (
                "bacteria_growth_age_multiplier_adult",
                g.bacteria_growth_age_multiplier_adult,
            ),
            (
                "bacteria_growth_age_multiplier_elderly",
                g.bacteria_growth_age_multiplier_elderly,
            ),
            (
                "bacteria_growth_immunodeficiency_multiplier",
                g.bacteria_growth_immunodeficiency_multiplier,
            ),
        ],
    );

    print_scalar_group(
        "Sepsis Onset",
        &[
            (
                "sepsis_minimum_duration_days",
                g.sepsis_minimum_duration_days as f64,
            ),
            (
                "log_odds_sepsis_onset_immunosuppressed",
                g.log_odds_sepsis_onset_immunosuppressed,
            ),
            (
                "log_odds_sepsis_onset_hospitalized",
                g.log_odds_sepsis_onset_hospitalized,
            ),
            (
                "log_odds_sepsis_onset_not_under_care",
                g.log_odds_sepsis_onset_not_under_care,
            ),
            (
                "log_odds_sepsis_onset_region_north_america",
                g.log_odds_sepsis_onset_region_north_america,
            ),
            (
                "log_odds_sepsis_onset_region_europe",
                g.log_odds_sepsis_onset_region_europe,
            ),
            (
                "log_odds_sepsis_onset_region_oceania",
                g.log_odds_sepsis_onset_region_oceania,
            ),
            (
                "log_odds_sepsis_onset_region_asia",
                g.log_odds_sepsis_onset_region_asia,
            ),
            (
                "log_odds_sepsis_onset_region_south_america",
                g.log_odds_sepsis_onset_region_south_america,
            ),
            (
                "log_odds_sepsis_onset_region_africa",
                g.log_odds_sepsis_onset_region_africa,
            ),
        ],
    );

    print_scalar_group(
        "Sepsis Recovery",
        &[
            (
                "sepsis_recovery_base_log_odds_per_day",
                g.sepsis_recovery_base_log_odds_per_day,
            ),
            (
                "sepsis_recovery_log_odds_bacteria_level",
                g.sepsis_recovery_log_odds_bacteria_level,
            ),
            (
                "sepsis_recovery_log_odds_in_hospital",
                g.sepsis_recovery_log_odds_in_hospital,
            ),
            (
                "sepsis_recovery_log_odds_age_infant",
                g.sepsis_recovery_log_odds_age_infant,
            ),
            (
                "sepsis_recovery_log_odds_age_child",
                g.sepsis_recovery_log_odds_age_child,
            ),
            (
                "sepsis_recovery_log_odds_age_adult",
                g.sepsis_recovery_log_odds_age_adult,
            ),
            (
                "sepsis_recovery_log_odds_age_elderly",
                g.sepsis_recovery_log_odds_age_elderly,
            ),
            (
                "sepsis_recovery_log_odds_immunosuppressed",
                g.sepsis_recovery_log_odds_immunosuppressed,
            ),
        ],
    );

    print_scalar_group(
        "Sepsis Death",
        &[
            ("sepsis_death_base_log_odds", g.sepsis_death_base_log_odds),
            (
                "sepsis_death_log_odds_age_infant",
                g.sepsis_death_log_odds_age_infant,
            ),
            (
                "sepsis_death_log_odds_age_child",
                g.sepsis_death_log_odds_age_child,
            ),
            (
                "sepsis_death_log_odds_age_adult",
                g.sepsis_death_log_odds_age_adult,
            ),
            (
                "sepsis_death_log_odds_age_elderly",
                g.sepsis_death_log_odds_age_elderly,
            ),
            (
                "sepsis_death_log_odds_immunosuppressed",
                g.sepsis_death_log_odds_immunosuppressed,
            ),
            (
                "sepsis_death_log_odds_bacteria_level",
                g.sepsis_death_log_odds_bacteria_level,
            ),
            (
                "sepsis_death_log_odds_duration",
                g.sepsis_death_log_odds_duration,
            ),
            (
                "sepsis_death_log_odds_early_phase",
                g.sepsis_death_log_odds_early_phase,
            ),
            (
                "sepsis_death_early_phase_days",
                g.sepsis_death_early_phase_days,
            ),
            (
                "sepsis_death_log_odds_not_under_care",
                g.sepsis_death_log_odds_not_under_care,
            ),
        ],
    );

    print_scalar_group(
        "Non-Sepsis Infection Mortality",
        &[
            (
                "infection_non_sepsis_base_log_odds",
                g.infection_non_sepsis_base_log_odds,
            ),
            (
                "infection_non_sepsis_log_odds_per_level",
                g.infection_non_sepsis_log_odds_per_level,
            ),
            (
                "infection_non_sepsis_log_odds_age_infant",
                g.infection_non_sepsis_log_odds_age_infant,
            ),
            (
                "infection_non_sepsis_log_odds_age_child",
                g.infection_non_sepsis_log_odds_age_child,
            ),
            (
                "infection_non_sepsis_log_odds_age_adult",
                g.infection_non_sepsis_log_odds_age_adult,
            ),
            (
                "infection_non_sepsis_log_odds_age_elderly",
                g.infection_non_sepsis_log_odds_age_elderly,
            ),
            (
                "infection_non_sepsis_log_odds_immunosuppressed",
                g.infection_non_sepsis_log_odds_immunosuppressed,
            ),
            (
                "infection_non_sepsis_log_odds_in_hospital",
                g.infection_non_sepsis_log_odds_in_hospital,
            ),
            (
                "infection_non_sepsis_minimum_bacteria_level",
                g.infection_non_sepsis_minimum_bacteria_level,
            ),
        ],
    );

    print_scalar_group(
        "Background Mortality",
        &[
            (
                "background_mortality_baseline_log_odds",
                g.background_mortality_baseline_log_odds,
            ),
            (
                "mortality_baseline_1930_multiplier",
                g.mortality_baseline_1930_multiplier,
            ),
            (
                "mortality_baseline_2035_multiplier",
                g.mortality_baseline_2035_multiplier,
            ),
            (
                "mortality_improvement_half_life_years",
                g.mortality_improvement_half_life_years,
            ),
            (
                "log_odds_mortality_per_year_of_age",
                g.log_odds_mortality_per_year_of_age,
            ),
            (
                "log_odds_mortality_per_year_of_age_squared",
                g.log_odds_mortality_per_year_of_age_squared,
            ),
            (
                "log_odds_mortality_immunosuppressed",
                g.log_odds_mortality_immunosuppressed,
            ),
            (
                "log_odds_mortality_hospitalized",
                g.log_odds_mortality_hospitalized,
            ),
        ],
    );

    print_scalar_group(
        "Drug Toxicity",
        &[
            (
                "default_toxicity_reservoir_half_life_days",
                g.default_toxicity_reservoir_half_life_days,
            ),
            (
                "toxicity_age_multiplier_infant",
                g.toxicity_age_multiplier_infant,
            ),
            (
                "toxicity_age_multiplier_child",
                g.toxicity_age_multiplier_child,
            ),
            (
                "toxicity_age_multiplier_adult",
                g.toxicity_age_multiplier_adult,
            ),
            (
                "toxicity_age_multiplier_elderly",
                g.toxicity_age_multiplier_elderly,
            ),
            (
                "toxicity_immunosuppressed_multiplier",
                g.toxicity_immunosuppressed_multiplier,
            ),
            (
                "toxicity_hospital_multiplier",
                g.toxicity_hospital_multiplier,
            ),
            (
                "toxicity_discontinuation_threshold",
                g.toxicity_discontinuation_threshold,
            ),
            (
                "toxicity_discontinuation_avoidance_days",
                g.toxicity_discontinuation_avoidance_days as f64,
            ),
        ],
    );

    print_scalar_group(
        "Regional Resistance Scoring",
        &[
            (
                "regional_resistance_threshold_very_high",
                g.regional_resistance_threshold_very_high,
            ),
            (
                "regional_resistance_threshold_high",
                g.regional_resistance_threshold_high,
            ),
            (
                "regional_resistance_threshold_moderate",
                g.regional_resistance_threshold_moderate,
            ),
            (
                "regional_resistance_penalty_very_high",
                g.regional_resistance_penalty_very_high,
            ),
            (
                "regional_resistance_penalty_high",
                g.regional_resistance_penalty_high,
            ),
            (
                "regional_resistance_penalty_moderate",
                g.regional_resistance_penalty_moderate,
            ),
        ],
    );

    print_scalar_group(
        "Therapy Scoring",
        &[
            (
                "targeted_therapy_narrow_spectrum_bonus",
                g.targeted_therapy_narrow_spectrum_bonus,
            ),
            (
                "targeted_therapy_broad_spectrum_penalty",
                g.targeted_therapy_broad_spectrum_penalty,
            ),
            (
                "targeted_therapy_ineffective_drug_penalty",
                g.targeted_therapy_ineffective_drug_penalty,
            ),
            (
                "effective_potency_threshold_for_targeted_therapy",
                g.effective_potency_threshold_for_targeted_therapy,
            ),
            (
                "empiric_therapy_broad_spectrum_bonus",
                g.empiric_therapy_broad_spectrum_bonus,
            ),
            (
                "empiric_therapy_ineffective_penalty",
                g.empiric_therapy_ineffective_penalty,
            ),
        ],
    );

    print_scalar_group(
        "MDR-TB Era Multipliers",
        &[
            (
                "mdr_tb_pre_antibiotic_era_multiplier",
                g.mdr_tb_pre_antibiotic_era_multiplier,
            ),
            (
                "mdr_tb_early_antibiotic_era_multiplier",
                g.mdr_tb_early_antibiotic_era_multiplier,
            ),
            (
                "mdr_tb_modern_era_multiplier",
                g.mdr_tb_modern_era_multiplier,
            ),
        ],
    );

    print_scalar_group(
        "Gonorrhoea Acquisition Era Multipliers",
        &[
            (
                "neisseria_gonorrhoeae_pre_1980_acquisition_multiplier",
                g.neisseria_gonorrhoeae_pre_1980_acquisition_multiplier,
            ),
            (
                "neisseria_gonorrhoeae_pre_2000_acquisition_multiplier",
                g.neisseria_gonorrhoeae_pre_2000_acquisition_multiplier,
            ),
            (
                "neisseria_gonorrhoeae_modern_acquisition_multiplier",
                g.neisseria_gonorrhoeae_modern_acquisition_multiplier,
            ),
        ],
    );

    print_scalar_group(
        "Diagnostic Testing",
        &[
            (
                "bacterial_testing_available_from_day",
                configured_value("bacterial_testing_available_from_day", 5478.0),
            ),
            (
                "resistance_testing_available_from_day",
                configured_value("resistance_testing_available_from_day", 9131.0),
            ),
            ("test_delay_days", configured_value("test_delay_days", 3.0)),
            (
                "resistance_test_result_delay_days",
                configured_value("resistance_test_result_delay_days", 2.0),
            ),
            (
                "test_r_error_probability",
                configured_value("test_r_error_probability", 0.02),
            ),
            (
                "test_r_error_value",
                configured_value("test_r_error_value", 0.25),
            ),
            (
                "bacterial_testing_base_rate_per_day",
                configured_value("bacterial_testing_base_rate_per_day", 0.15),
            ),
            (
                "bacterial_testing_initial_adoption_rate",
                configured_value("bacterial_testing_initial_adoption_rate", 0.1),
            ),
            (
                "bacterial_testing_max_temporal_multiplier",
                configured_value("bacterial_testing_max_temporal_multiplier", 1.0),
            ),
            (
                "bacterial_testing_hospital_multiplier",
                configured_value("bacterial_testing_hospital_multiplier", 8.0),
            ),
            (
                "resistance_testing_base_rate_per_day",
                configured_value("resistance_testing_base_rate_per_day", 0.95),
            ),
            (
                "resistance_testing_initial_adoption_rate",
                configured_value("resistance_testing_initial_adoption_rate", 0.05),
            ),
            (
                "resistance_testing_max_temporal_multiplier",
                configured_value("resistance_testing_max_temporal_multiplier", 1.0),
            ),
            (
                "resistance_testing_hospital_multiplier",
                configured_value("resistance_testing_hospital_multiplier", 5.0),
            ),
            (
                "testing_immunosuppressed_multiplier",
                configured_value("testing_immunosuppressed_multiplier", 2.5),
            ),
            (
                "testing_sepsis_multiplier",
                configured_value("testing_sepsis_multiplier", 4.0),
            ),
        ],
    );

    print_scalar_group(
        "Additional Resistance and Treatment Controls",
        &[
            (
                "microbiome_majority_threshold",
                configured_value("microbiome_majority_threshold", 0.1),
            ),
            (
                "majority_r_evolution_rate_per_day_when_drug_present",
                configured_value("majority_r_evolution_rate_per_day_when_drug_present", 0.0),
            ),
            (
                "microbiome_clearance_probability_on_drug_treatment",
                configured_value("microbiome_clearance_probability_on_drug_treatment", 0.8),
            ),
            (
                "mdr_mycobacterium_tuberculosis_multi_drug_synergy_threshold",
                configured_value(
                    "mdr_mycobacterium_tuberculosis_multi_drug_synergy_threshold",
                    2.0,
                ),
            ),
            (
                "mdr_mycobacterium_tuberculosis_multi_drug_synergy_multiplier",
                configured_value(
                    "mdr_mycobacterium_tuberculosis_multi_drug_synergy_multiplier",
                    2.5,
                ),
            ),
            (
                "mdr_mycobacterium_tuberculosis_background_drug_effectiveness",
                configured_value(
                    "mdr_mycobacterium_tuberculosis_background_drug_effectiveness",
                    0.8,
                ),
            ),
            (
                "mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance",
                configured_value(
                    "mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance",
                    0.9,
                ),
            ),
        ],
    );

    print_scalar_group(
        "Run-Level Resistance Pathway Controls",
        &[
            (
                "run_pathway_infection_de_novo_multiplier",
                configured_value("run_pathway_infection_de_novo_multiplier", 1.0),
            ),
            (
                "run_pathway_reversion_rate_multiplier",
                configured_value("run_pathway_reversion_rate_multiplier", 1.0),
            ),
            (
                "run_pathway_hgt_multiplier",
                configured_value("run_pathway_hgt_multiplier", 1.0),
            ),
            (
                "run_pathway_microbiome_acquisition_multiplier",
                configured_value("run_pathway_microbiome_acquisition_multiplier", 1.0),
            ),
            (
                "run_pathway_ratchet_enabled",
                configured_value("run_pathway_ratchet_enabled", 1.0),
            ),
        ],
    );
}

fn print_scalar_group(title: &str, items: &[(&str, f64)]) {
    println!("#### {}", title);
    println!();
    let rows: Vec<Vec<String>> = items
        .iter()
        .map(|&(name, value)| vec![name.to_string(), format_value(value)])
        .collect();
    md_table(&["Parameter", "Value"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.2  Drug Properties
// ─────────────────────────────────────────────────────────────────────────────

fn print_drug_properties(store: &amr_project::config::ParameterStore) {
    println!("### B.2 Drug Properties");
    println!();
    println!(
        "Pharmacokinetic and clinical properties for each of the {} modelled \
              antimicrobial agents. The introduction time step is measured in days \
              from 1 January 1930.",
        DRUG_SHORT_NAMES.len()
    );
    println!();
    println!(
        "See: \
              [§6.3 Drug pharmacokinetics](#63-drug-pharmacokinetics), \
              [§6.5 Drug potency matrix](#65-drug-potency-matrix), \
              [§6.6 Drug availability](#66-drug-availability-by-region-and-era), \
              [§6.7 Drug toxicity](#67-drug-toxicity), \
              [§6.8 Antibiotic infection prevention](#68-antibiotic-infection-prevention)."
    );
    println!();

    let headers = &[
        "Drug",
        "Class",
        "Intro (days)",
        "Init level",
        "t½ (days)",
        "2× dose mult",
        "Spectrum",
        "Tox hazard",
        "Tox t½ (days)",
        "Microbiome disrupt",
    ];
    let mut rows = Vec::new();
    for (d_idx, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
        let drug_class = get_drug_class(drug).unwrap_or("unknown");
        let intro = get_drug_introduction_time_step(drug)
            .map(|v| v.to_string())
            .unwrap_or_else(|| "n/a".to_string());
        rows.push(vec![
            drug.to_string(),
            drug_class.to_string(),
            intro,
            format_value(store.drug.initial_level(d_idx)),
            format_value(store.drug.half_life_days(d_idx)),
            format_value(store.drug.double_dose_multiplier(d_idx)),
            format_value(store.drug.spectrum_breadth(d_idx)),
            format_value(store.drug.toxicity_death_hazard_per_unit_level(d_idx)),
            format_value(store.drug.toxicity_reservoir_half_life_days(d_idx)),
            format_value(store.drug.microbiome_disruption_log_odds(d_idx)),
        ]);
    }
    md_table(headers, &rows);

    println!("#### Non-Default Regional Drug Availability");
    println!();
    println!(
        "Regional availability defaults to 1.0. Only configured values that differ from that \
              default are shown. The separate time-aware availability rules described in \
              Section 6.6 are implementation rules rather than entries in this table."
    );
    println!();
    let mut rows = Vec::new();
    for region in REGION_NAMES {
        for &drug in DRUG_SHORT_NAMES.iter() {
            let key = format!("{}_drug_{}_availability", region, drug);
            let value = configured_value(&key, 1.0);
            if (value - 1.0).abs() > 1e-12 {
                rows.push(vec![
                    region.to_string(),
                    drug.to_string(),
                    format_value(value),
                ]);
            }
        }
    }
    md_table(&["Region", "Drug", "Availability multiplier"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.3  Bacteria Properties
// ─────────────────────────────────────────────────────────────────────────────

fn print_bacteria_properties(store: &amr_project::config::ParameterStore) {
    let b = &store.bacteria;
    println!("### B.3 Bacteria Properties");
    println!();
    println!(
        "Per-bacteria parameters governing acquisition, growth, symptom \
              onset, and clinical outcomes for each of the {} bacterial species.",
        BACTERIA_LIST.len()
    );
    println!();
    println!(
        "See: \
              [§3.1 Community acquisition](#31-community-acquisition), \
              [§4.2 Infection dynamics](#42-infection-dynamics), \
              [§4.3 Sepsis](#43-sepsis), \
              [§4.4 Natural clearance and microbiome dynamics](#44-natural-clearance-and-microbiome-dynamics), \
              [§8.1 Carriage compartments](#81-carriage-compartments)."
    );
    println!();

    println!("#### Acquisition, Growth, and Carriage");
    println!();
    let headers = &[
        "Bacteria",
        "Acq log-odds",
        "Vaccinated log-odds",
        "Carriage-present log-odds",
        "Hospital-acquired log-odds",
        "Init level",
        "Delta level/day",
        "Max level",
        "Carriage clearance/day",
        "Carriage vs infection log-odds",
        "Separate carriage state",
    ];
    let mut rows = Vec::new();
    for (idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        rows.push(vec![
            bacteria.to_string(),
            format_value(b.acquisition_log_odds_baseline[idx]),
            format_value(b.log_odds_vaccinated[idx]),
            format_value(b.log_odds_microbiome_present[idx]),
            format_value(b.log_odds_hospital_acquired[idx]),
            format_value(b.initial_infection_level[idx]),
            format_value(b.base_bacteria_level_change[idx]),
            format_value(b.max_level[idx]),
            format_value(b.microbiome_clearance_probability_per_day[idx]),
            format_value(b.microbiome_vs_infection_log_odds[idx]),
            if bacterium_has_separate_microbiome_compartment(idx) {
                "yes".to_string()
            } else {
                "no".to_string()
            },
        ]);
    }
    md_table(headers, &rows);

    println!("#### Symptoms and Treatment Tracking");
    println!();
    let headers = &[
        "Bacteria",
        "Symptom base log-odds",
        "Symptom threshold",
        "Symptom delay (days)",
        "Symptom log-odds/level",
        "Drug cessation probability",
        "Treatment recognition year",
        "Failure: no immediate second line",
    ];
    let mut rows = Vec::new();
    for (idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        rows.push(vec![
            bacteria.to_string(),
            format_value(b.symptom_onset_base_log_odds[idx]),
            format_value(b.symptom_onset_threshold_level[idx]),
            format_value(b.symptom_onset_delay_days[idx]),
            format_value(b.symptom_onset_log_odds_per_level_unit[idx]),
            format_value(b.drug_cessation_probability[idx]),
            b.treatment_recognition_year[idx]
                .map(format_value)
                .unwrap_or_else(|| "none".to_string()),
            format_value(b.treatment_failure_no_second_line_probability[idx]),
        ]);
    }
    md_table(headers, &rows);

    println!("#### Clinical Outcomes and Resistance Ecology");
    println!();
    let headers = &[
        "Bacteria",
        "Sepsis base log-odds",
        "Sepsis log-odds/level",
        "Sepsis log-odds/day",
        "Non-sepsis death log-odds",
        "Sepsis-death override",
        "Mechanismless reversion/day",
        "Community human-profile probability",
        "Hospital susceptible prune %",
        "Community mechanism-reversion multiplier",
    ];
    let mut rows = Vec::new();
    for (idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        rows.push(vec![
            bacteria.to_string(),
            format_value(b.sepsis_baseline_log_odds[idx]),
            format_value(b.sepsis_log_odds_infection_level[idx]),
            format_value(b.sepsis_log_odds_infection_duration[idx]),
            format_value(b.infection_non_sepsis_mortality_log_odds[idx]),
            format_value(b.sepsis_death_log_odds_override[idx]),
            format_value(b.mechanismless_resistance_reversion_rate[idx]),
            format_value(b.community_resistance_dilution_factor[idx]),
            format_value(b.hospital_resistance_prune_susceptible_percent[idx]),
            format_value(b.community_mechanism_reversion_multiplier[idx]),
        ]);
    }
    md_table(headers, &rows);

    println!("#### Bacterium-Specific Testing Availability Years");
    println!();
    println!(
        "Only explicit bacterium-specific overrides are shown; all other organisms use the \
              general bacterial-testing availability date in B.1."
    );
    println!();
    let rows = configured_rows_matching(|key| {
        BACTERIA_LIST
            .iter()
            .any(|bacteria| key == format!("{}_test_availability_year", bacteria))
    });
    md_table(&["Parameter", "Year"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.4  Drug–Bacteria Potency Matrix
// ─────────────────────────────────────────────────────────────────────────────

fn print_drug_bacteria_matrix(store: &amr_project::config::ParameterStore) {
    println!("### B.4 Drug–Bacteria Potency Matrix");
    println!();
    println!(
        "Baseline potency (MIC-derived effectiveness when no resistance is \
              present) and initiation multiplier (stewardship weighting for drug \
              selection) for each drug–bacteria pair. {} bacteria × {} drugs = {} entries.",
        BACTERIA_LIST.len(),
        DRUG_SHORT_NAMES.len(),
        BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len()
    );
    println!();
    println!(
        "See: \
              [§6.5 Drug potency matrix](#65-drug-potency-matrix), \
              [§6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use)."
    );
    println!();

    let headers = &["Bacteria", "Drug", "Potency (no R)", "Init multiplier"];
    let mut rows = Vec::new();
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        for (d_idx, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            let potency = store.drug_bacteria.potency(b_idx, d_idx);
            let init_mult = store.drug_bacteria.initiation_multiplier(b_idx, d_idx);
            rows.push(vec![
                bacteria.to_string(),
                drug.to_string(),
                format_value(potency),
                format_value(init_mult),
            ]);
        }
    }
    md_table(headers, &rows);

    println!("#### Time-Varying Drug-Initiation Overrides");
    println!();
    println!(
        "These values replace the base initiation multiplier before the year encoded in the \
              parameter name. For overlapping cut-offs, the earliest cut-off later than the \
              current simulation year is used."
    );
    println!();
    let rows = configured_rows_matching(|key| key.contains("_initiation_multiplier_before_"));
    md_table(&["Parameter", "Multiplier"], &rows);

    println!("#### Additional Clinical-Preference Multipliers");
    println!();
    println!(
        "These directly read bacterium-drug multipliers default to 1.0. Only explicit \
              overrides are shown."
    );
    println!();
    let rows = configured_rows_matching(|key| key.ends_with("_clinical_preference_multiplier"));
    md_table(&["Parameter", "Multiplier"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.5  Regional Parameters
// ─────────────────────────────────────────────────────────────────────────────

fn print_regional_parameters(store: &amr_project::config::ParameterStore) {
    println!("### B.5 Regional Parameters");
    println!();
    println!(
        "Region-level scalars (applicable to all bacteria) and the per-region \
              per-bacteria acquisition log-odds adjustments."
    );
    println!();
    println!(
        "See: \
              [§2.5 Travel](#25-travel), \
              [§3.1 Community acquisition](#31-community-acquisition)."
    );
    println!();

    // Region scalars
    println!("#### Region Scalars");
    println!();
    let headers = &[
        "Region",
        "Travel mult",
        "Cessation mult",
        "Mortality log-odds",
        "Sepsis log-odds",
        "Sepsis mort mult",
        "Testing mult",
        "Abx init log-odds",
        "Hosp log-odds",
    ];
    let mut rows = Vec::new();
    for (idx, &region) in REGION_VARIANTS.iter().enumerate() {
        let name = if idx < REGION_NAMES.len() {
            REGION_NAMES[idx]
        } else {
            "home"
        };
        rows.push(vec![
            name.to_string(),
            format_value(store.region.travel_multiplier(region)),
            format_value(store.region.cessation_multiplier(region)),
            format_value(store.region.mortality_log_odds(region)),
            format_value(store.region.sepsis_recovery_log_odds(region)),
            format_value(store.region.sepsis_mortality_multiplier(region)),
            format_value(store.region.testing_multiplier(region)),
            format_value(store.region.antibiotic_initiation_log_odds(region)),
            format_value(store.region.hospitalization_log_odds(region)),
        ]);
    }
    md_table(headers, &rows);

    // Region-bacteria acquisition
    println!("#### Region–Bacteria Acquisition Log-Odds");
    println!();
    let headers = &["Region", "Bacteria", "Acquisition log-odds"];
    let mut rows = Vec::new();
    for &region in REGION_VARIANTS.iter() {
        let name = region_name(region);
        for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
            let val = store.region_bacteria.acquisition_log_odds(region, b_idx);
            if val.abs() > 1e-12 {
                rows.push(vec![
                    name.to_string(),
                    bacteria.to_string(),
                    format_value(val),
                ]);
            }
        }
    }
    md_table(headers, &rows);
}

fn region_name(region: Region) -> &'static str {
    match region {
        Region::NorthAmerica => "north_america",
        Region::SouthAmerica => "south_america",
        Region::Africa => "africa",
        Region::Asia => "asia",
        Region::Europe => "europe",
        Region::Oceania => "oceania",
        Region::Home => "home",
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// B.6  Age-Dependent Parameters
// ─────────────────────────────────────────────────────────────────────────────

fn print_age_dependent_parameters(store: &amr_project::config::ParameterStore) {
    println!("### B.6 Age-Dependent Parameters");
    println!();
    println!(
        "Log-odds adjustments by age category for bacteria acquisition and \
              regional effects. Age categories: {}.",
        AGE_CATEGORY_SEQUENCE
            .iter()
            .map(|c| c.label())
            .collect::<Vec<_>>()
            .join(", ")
    );
    println!();
    println!(
        "See: \
              [§2.2 Ageing and age categories](#22-ageing-and-age-categories), \
              [§3.1 Community acquisition](#31-community-acquisition)."
    );
    println!();

    // Default age log-odds
    println!("#### Default Age Log-Odds");
    println!();
    let mut rows = Vec::new();
    for (idx, &cat) in AGE_CATEGORY_SEQUENCE.iter().enumerate() {
        rows.push(vec![
            cat.label().to_string(),
            format_value(store.age_categories.default_log_odds(idx)),
        ]);
    }
    md_table(&["Age category", "Default log-odds"], &rows);

    // Per-bacteria age log-odds
    println!("#### Bacteria–Age Log-Odds");
    println!();
    let mut rows = Vec::new();
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        for (a_idx, &cat) in AGE_CATEGORY_SEQUENCE.iter().enumerate() {
            let val = store.age_categories.bacteria_age_log_odds(b_idx, a_idx);
            if val.abs() > 1e-12 {
                rows.push(vec![
                    bacteria.to_string(),
                    cat.label().to_string(),
                    format_value(val),
                ]);
            }
        }
    }
    md_table(&["Bacteria", "Age category", "Log-odds"], &rows);

    // Region-age log-odds
    println!("#### Region–Age Log-Odds");
    println!();
    let mut rows = Vec::new();
    for &region in REGION_VARIANTS.iter() {
        let name = region_name(region);
        for (a_idx, &cat) in AGE_CATEGORY_SEQUENCE.iter().enumerate() {
            let val = store.age_categories.region_age_log_odds(region, a_idx);
            if val.abs() > 1e-12 {
                rows.push(vec![
                    name.to_string(),
                    cat.label().to_string(),
                    format_value(val),
                ]);
            }
        }
    }
    md_table(&["Region", "Age category", "Log-odds"], &rows);

    println!("#### Explicit Bacterium–Region–Age Overrides");
    println!();
    println!(
        "Only explicitly configured three-way overrides are shown. Every unlisted combination \
              inherits the corresponding region-age value above."
    );
    println!();
    let mut rows = Vec::new();
    for &region in REGION_VARIANTS.iter() {
        let region = region_name(region);
        for &bacteria in BACTERIA_LIST.iter() {
            for category in AGE_CATEGORY_SEQUENCE {
                let key = format!("{}_{}_log_odds_{}", bacteria, region, category.label());
                if let Some(value) = PARAMETERS.get(&key) {
                    rows.push(vec![
                        bacteria.to_string(),
                        region.to_string(),
                        category.label().to_string(),
                        format_value(*value),
                    ]);
                }
            }
        }
    }
    md_table(&["Bacteria", "Region", "Age category", "Log-odds"], &rows);

    println!("#### Sepsis-Onset Age Log-Odds");
    println!();
    let rows = vec![
        vec![
            "sepsis_age_log_odds_baseline".to_string(),
            format_value(configured_value("sepsis_age_log_odds_baseline", 0.0)),
        ],
        vec![
            "sepsis_age_log_odds_neonatal".to_string(),
            format_value(configured_value("sepsis_age_log_odds_neonatal", 0.0)),
        ],
        vec![
            "sepsis_age_log_odds_pediatric".to_string(),
            format_value(configured_value("sepsis_age_log_odds_pediatric", 0.0)),
        ],
        vec![
            "sepsis_age_log_odds_young_adult".to_string(),
            format_value(configured_value("sepsis_age_log_odds_young_adult", 0.0)),
        ],
        vec![
            "sepsis_age_log_odds_elderly".to_string(),
            format_value(configured_value("sepsis_age_log_odds_elderly", 0.0)),
        ],
    ];
    md_table(&["Parameter", "Log-odds"], &rows);

    println!("#### Bacterium-Specific Sepsis-Age Overrides");
    println!();
    println!("Only explicit overrides are shown; all other combinations contribute 0.");
    println!();
    let rows = configured_rows_matching(|key| {
        ["neonatal", "pediatric", "young_adult", "elderly"]
            .iter()
            .any(|age| key.ends_with(&format!("_{}_sepsis_log_odds", age)))
    });
    md_table(&["Parameter", "Log-odds"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.7  Syndrome Parameters
// ─────────────────────────────────────────────────────────────────────────────

fn print_syndrome_parameters(store: &amr_project::config::ParameterStore) {
    println!("### B.7 Syndrome Parameters");
    println!();
    println!(
        "Infection-site (syndrome) specific parameters. Syndromes are: \
              1 = UTI, 2 = skin/soft tissue, 3 = respiratory, 4 = bloodstream, \
              5 = intra-abdominal, 6 = CNS/meningitis, 7 = gastrointestinal, \
              8 = genital/STI, 9 = bone/joint, 10 = other."
    );
    println!();
    println!(
        "See: \
              [§4.1 Syndrome assignment](#41-syndrome-assignment), \
              [§6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use), \
              [§6.4 Drug penetration by syndrome](#64-drug-penetration-by-syndrome)."
    );
    println!();

    println!("#### Syndrome-Level Clinical Scalars");
    println!();
    let mut rows = Vec::new();
    for syndrome_id in 0..=10 {
        rows.push(vec![
            SYNDROME_NAMES[syndrome_id].to_string(),
            format_value(store.syndrome.sepsis_log_odds(syndrome_id)),
            format_value(store.syndrome.initiation_multiplier(syndrome_id)),
            format_value(store.syndrome.non_sepsis_mortality_log_odds(syndrome_id)),
            format_value(store.syndrome.bacteria_growth_multiplier(syndrome_id)),
        ]);
    }
    md_table(
        &[
            "Syndrome",
            "Sepsis log-odds",
            "Initiation multiplier",
            "Non-sepsis death log-odds",
            "Growth multiplier",
        ],
        &rows,
    );

    // Empiric drug scores
    println!("#### Non-Default Syndrome Empiric Drug Scores");
    println!();
    println!("The resolved default for every unlisted syndrome-drug pair is 0.01.");
    println!();
    let mut rows = Vec::new();
    for syndrome_id in 1..=10 {
        for (d_idx, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            let score = store.syndrome.empiric_drug_score(syndrome_id, d_idx);
            if (score - 0.01).abs() > 1e-6 {
                rows.push(vec![
                    SYNDROME_NAMES[syndrome_id].to_string(),
                    drug.to_string(),
                    format_value(score),
                ]);
            }
        }
    }
    md_table(&["Syndrome", "Drug", "Empiric score"], &rows);

    println!("#### Time-Varying Syndrome Empiric-Score Overrides");
    println!();
    println!(
        "These values replace the base syndrome score before the year encoded in the \
              parameter name."
    );
    println!();
    let rows = configured_rows_matching(|key| {
        key.starts_with("syndrome_") && key.contains("_score_before_")
    });
    md_table(&["Parameter", "Empiric score"], &rows);

    // Drug penetration
    println!("#### Non-Default Syndrome Drug Penetration");
    println!();
    println!("The resolved default for every unlisted syndrome-drug pair is 1.0.");
    println!();
    let mut rows = Vec::new();
    for syndrome_id in 1..=10 {
        for (d_idx, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            let pen = store.syndrome.drug_penetration(syndrome_id, d_idx);
            if (pen - 1.0).abs() > 1e-6 {
                rows.push(vec![
                    SYNDROME_NAMES[syndrome_id].to_string(),
                    drug.to_string(),
                    format_value(pen),
                ]);
            }
        }
    }
    md_table(&["Syndrome", "Drug", "Penetration factor"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.8  Clearance Parameters
// ─────────────────────────────────────────────────────────────────────────────

fn print_clearance_parameters(store: &amr_project::config::ParameterStore) {
    println!("### B.8 Clearance Parameters");
    println!();
    println!(
        "Infection clearance model parameters. The clearance hazard is a \
              logistic function of base log-odds, per-bacteria adjustments, \
              age effects, immunodeficiency, bacteria level, and infection duration."
    );
    println!();
    println!(
        "See: [§4.4 Natural clearance and microbiome dynamics](#44-natural-clearance-and-microbiome-dynamics)."
    );
    println!();

    let rows = vec![
        vec![
            "default_clearance_delay_days".to_string(),
            format_value(configured_value("default_clearance_delay_days", 3.0)),
        ],
        vec![
            "base_clearance_log_odds".to_string(),
            format_value(store.clearance.base_log_odds()),
        ],
        vec![
            "immunodeficient_log_odds_adjustment".to_string(),
            format_value(store.clearance.immunodeficient_log_odds_adjustment()),
        ],
        vec![
            "clearance_level_log_odds_per_unit".to_string(),
            format_value(configured_value("clearance_level_log_odds_per_unit", -0.3)),
        ],
        vec![
            "adaptive_recruit_slope_per_infection_day (implementation constant)".to_string(),
            format_value(0.25),
        ],
    ];
    md_table(&["Parameter", "Value"], &rows);

    println!(
        "`default_clearance_delay_days` and any bacterium-specific `*_clearance_delay_days` \
         values are loaded for compatibility but are not consulted by the current clearance \
         hazard. Eligibility is instead immediate after acquisition, and infection duration \
         enters through the fixed +0.25 log-odds/day term."
    );
    println!();

    println!("#### Clearance Age Adjustments");
    println!();
    let mut rows = Vec::new();
    for category in AGE_CATEGORY_SEQUENCE {
        let key = format!("clearance_age_log_odds_{}", category.label());
        rows.push(vec![
            category.label().to_string(),
            format_value(configured_value(&key, 0.0)),
        ]);
    }
    md_table(&["Age category", "Log-odds adjustment"], &rows);

    // Per-bacteria
    println!("#### Per-Bacteria Clearance Adjustments");
    println!();
    let mut rows = Vec::new();
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        let adj = store.clearance.bacteria_log_odds_adjustment(b_idx);
        if adj.abs() > 1e-12 {
            rows.push(vec![bacteria.to_string(), format_value(adj)]);
        }
    }
    md_table(&["Bacteria", "Log-odds adjustment"], &rows);

    println!("#### Configured Per-Bacterium Clearance Delays");
    println!();
    println!(
        "Only explicit overrides are shown. As noted above, these loaded values are currently \
              inactive in the executable hazard."
    );
    println!();
    let rows = configured_rows_matching(|key| {
        key.ends_with("_clearance_delay_days") && key != "default_clearance_delay_days"
    });
    md_table(&["Parameter", "Days"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.9  Immunodeficiency, Sex & Vaccination
// ─────────────────────────────────────────────────────────────────────────────

fn print_immunodeficiency_sex_vaccination(store: &amr_project::config::ParameterStore) {
    println!("### B.9 Immunodeficiency, Sex, and Vaccination Parameters");
    println!();
    println!(
        "See: \
              [§2.3 Immunodeficiency](#23-immunodeficiency), \
              [§10 Mortality](#10-mortality)."
    );
    println!();

    // Immunodeficiency
    println!("#### Immunodeficiency");
    println!();
    let mut rows = vec![
        vec![
            "startup_seed_fraction".to_string(),
            format_value(store.immunodeficiency.startup_seed_fraction()),
        ],
        vec![
            "temporary_onset_rate_per_day".to_string(),
            format_value(store.immunodeficiency.temporary_onset_rate()),
        ],
        vec![
            "temporary_recovery_rate_per_day".to_string(),
            format_value(store.immunodeficiency.temporary_recovery_rate()),
        ],
        vec![
            "chronic_onset_rate_per_day".to_string(),
            format_value(store.immunodeficiency.chronic_onset_rate()),
        ],
        vec![
            "chronic_recovery_rate_per_day".to_string(),
            format_value(store.immunodeficiency.chronic_recovery_rate()),
        ],
    ];
    for &(label, age_days) in &[
        ("age_0_1", 180),
        ("age_1_18", 3650),
        ("age_18_65", 14600),
        ("age_65_plus", 25550),
    ] {
        rows.push(vec![
            format!("chronic_probability_{}", label),
            format_value(store.immunodeficiency.chronic_probability(age_days)),
        ]);
    }
    md_table(&["Parameter", "Value"], &rows);

    // Sex
    println!("#### Sex");
    println!();
    let rows = vec![
        vec![
            "male".to_string(),
            format_value(store.sex.mortality_log_odds("male")),
        ],
        vec![
            "female".to_string(),
            format_value(store.sex.mortality_log_odds("female")),
        ],
    ];
    md_table(&["Sex", "Mortality log-odds"], &rows);

    // Vaccination
    println!("#### Vaccination");
    println!();
    let mut rows = Vec::new();
    for (v_idx, &vaccine) in VACCINES.iter().enumerate() {
        let avail = store.vaccination.availability_year(v_idx);
        let birth_coverage = store.vaccination.birth_coverage_target(v_idx);
        let rollout_years = store.vaccination.rollout_years(v_idx);
        rows.push(vec![
            vaccine.to_string(),
            format_value(avail),
            format_value(birth_coverage),
            format_value(rollout_years),
        ]);
    }
    md_table(
        &[
            "Vaccine",
            "Availability year",
            "Target birth-cohort coverage",
            "Rollout years",
        ],
        &rows,
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// B.10  Resistance Mechanisms
// ─────────────────────────────────────────────────────────────────────────────

fn print_resistance_mechanisms(store: &amr_project::config::ParameterStore) {
    println!("### B.10 Resistance Mechanisms");
    println!();
    println!(
        "Parameters for the {} resistance mechanisms modelled. Each mechanism \
              has a per-day reversion rate, per-drug-class enhancement multipliers, \
              and per-bacteria emergence rates.",
        ResistanceMechanism::all().len()
    );
    println!();
    println!("See: \
              [§7.1 Resistance mechanisms](#71-resistance-mechanisms), \
              [§7.2 Mechanism–drug-class enhancement](#72-mechanismdrug-class-enhancement-multipliers), \
              [§7.3 Resistance emergence](#73-resistance-emergence), \
              [§7.4 Resistance reversion](#74-resistance-reversion-and-fitness-costs).");
    println!();

    let mechanisms = ResistanceMechanism::all();
    let drug_classes = DrugClass::all();

    // Reversion rates
    println!("#### Mechanism Reversion Rates");
    println!();
    let mut rows = Vec::new();
    for (m_idx, mechanism) in mechanisms.iter().enumerate() {
        let rate = store.resistance_mechanism.reversion_rate(m_idx);
        rows.push(vec![mechanism.as_str().to_string(), format_value(rate)]);
    }
    md_table(&["Mechanism", "Reversion rate/day"], &rows);

    // Enhancement multipliers
    println!("#### Mechanism Enhancement Multipliers by Drug Class");
    println!();
    println!(
        "Raw class enhancement values loaded for each mechanism. These values are applied only \
              to bacterium-drug pairs admitted by the executable host and drug-specific \
              applicability gates in `rules::mechanism_applies_to_drug`; non-applicable fallback \
              values shown here are inert. The resolved default for every unlisted \
              mechanism-class pair is 0."
    );
    println!();
    let mut rows = Vec::new();
    for (m_idx, mechanism) in mechanisms.iter().enumerate() {
        for drug_class in drug_classes {
            let enh = store
                .resistance_mechanism
                .enhancement_multiplier(m_idx, drug_class.index());
            if enh.abs() > 1e-12 {
                rows.push(vec![
                    mechanism.as_str().to_string(),
                    drug_class.as_str().to_string(),
                    format_value(enh),
                ]);
            }
        }
    }
    md_table(
        &["Mechanism", "Drug class", "Enhancement multiplier"],
        &rows,
    );

    // Bacteria-mechanism emergence rates
    println!("#### Bacteria–Mechanism Emergence Rates");
    println!();
    println!(
        "Resolved de novo emergence rate and executable pathway status for every \
              bacteria–mechanism pair. A zero rate does not necessarily exclude the host: \
              transferable mechanisms can remain HGT-only, while non-transferable eligible \
              mechanisms can still be inherited in an existing complete profile."
    );
    println!();
    let mut rows = Vec::new();
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        for (m_idx, mechanism) in mechanisms.iter().enumerate() {
            let rate = store.bacteria_mechanism_emergence.rate(b_idx, m_idx);
            let status = store.bacteria_mechanism_status.status(b_idx, m_idx);
            rows.push(vec![
                bacteria.to_string(),
                mechanism.as_str().to_string(),
                format_value(rate),
                mechanism_status_label(status).to_string(),
            ]);
        }
    }
    md_table(
        &["Bacteria", "Mechanism", "Emergence rate/day", "Status"],
        &rows,
    );

    println!("#### Environmental and Exogenous Mechanism Floors");
    println!();
    println!(
        "All unspecified bacteria–mechanism floors resolve to 0. The table lists every \
              explicit base or `_before_YYYY` override, including explicit zeroes that mark \
              the start of an era sequence."
    );
    println!();
    let rows = configured_rows_matching(|key| key.contains("_environmental_floor"));
    md_table(&["Parameter", "Assignment probability"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.11  Horizontal Gene Transfer Matrix
// ─────────────────────────────────────────────────────────────────────────────

fn print_hgt_matrix(store: &amr_project::config::ParameterStore) {
    println!("### B.11 Horizontal Gene Transfer Matrix");
    println!();
    println!(
        "Per-day probability of horizontal gene transfer of resistance \
              between co-colonising bacterial species. Only non-zero entries shown."
    );
    println!();
    println!(
        "See: \
              [§9.1 Transfer compatibility](#91-transfer-compatibility), \
              [§9.2 The HGT process](#92-the-hgt-process)."
    );
    println!();

    let mut rows = Vec::new();
    for (donor_idx, &donor) in BACTERIA_LIST.iter().enumerate() {
        for (recip_idx, &recipient) in BACTERIA_LIST.iter().enumerate() {
            if donor_idx == recip_idx {
                continue;
            }
            let prob = store.hgt.probability(donor_idx, recip_idx);
            if prob.abs() > 1e-20 {
                rows.push(vec![
                    donor.to_string(),
                    recipient.to_string(),
                    format_value(prob),
                ]);
            }
        }
    }
    md_table(&["Donor", "Recipient", "Probability/day"], &rows);
}
