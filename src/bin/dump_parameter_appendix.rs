/// Generates the Appendix B markdown for MODEL_DESCRIPTION.md.
///
/// Prints structured, thematically organized parameter tables derived from the
/// live Rust configuration, replacing the previous monolithic key-value dump
/// with resolved, reader-friendly Markdown tables.
use amr_project::config::{
    get_drug_class, get_drug_introduction_time_step, PARAMETERS, PARAMETER_STORE,
};
use amr_project::simulation::population::{
    DrugClass, Region, ResistanceMechanism, AGE_CATEGORY_SEQUENCE, BACTERIA_LIST,
    DRUG_SHORT_NAMES,
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

const VACCINES: [&str; 3] = ["pneumococcal", "meningococcal", "hib"];

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
    println!("This appendix is auto-generated from the live Rust configuration. \
              Parameters are organized thematically into resolved tables \
              derived from the internal data structures. All values shown are \
              the effective defaults before any run-level sampling multipliers \
              are applied.");
    println!();
}

// ─────────────────────────────────────────────────────────────────────────────
// B.1  Global Scalar Parameters
// ─────────────────────────────────────────────────────────────────────────────

fn print_global_scalars(store: &amr_project::config::ParameterStore) {
    let g = &store.globals;
    println!("### B.1 Global Scalar Parameters");
    println!();
    println!("Scalar parameters that govern cross-cutting model behaviour. \
              Grouped thematically; each row gives the parameter name and its \
              default value.");
    println!();
    println!("See: \
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
              [§10 Mortality](#10-mortality).");
    println!();

    print_scalar_group("Treatment Initiation (logistic model)", &[
        ("antibiotic_initiation_base_log_odds", g.antibiotic_initiation_base_log_odds),
        ("antibiotic_initiation_log_odds_symptomatic_infection", g.antibiotic_initiation_log_odds_symptomatic_infection),
        ("antibiotic_initiation_log_odds_test_identified", g.antibiotic_initiation_log_odds_test_identified),
        ("antibiotic_initiation_log_odds_already_on_drug", g.antibiotic_initiation_log_odds_already_on_drug),
        ("antibiotic_initiation_log_odds_immunodeficiency", g.antibiotic_initiation_log_odds_immunodeficiency),
        ("antibiotic_initiation_log_odds_sepsis", g.antibiotic_initiation_log_odds_sepsis),
        ("antibiotic_initiation_log_odds_hospitalized", g.antibiotic_initiation_log_odds_hospitalized),
        ("antibiotic_initiation_log_odds_no_indication", g.antibiotic_initiation_log_odds_no_indication),
    ]);

    print_scalar_group("Drug Activity and Cessation", &[
        ("drug_activity_to_bacteria_level_multiplier", g.drug_activity_to_bacteria_level_multiplier),
        ("drug_activity_slow_clearance_probability", g.drug_activity_slow_clearance_probability),
        ("drug_activity_slow_clearance_multiplier", g.drug_activity_slow_clearance_multiplier),
        ("double_dose_probability_if_identified_infection", g.double_dose_probability_if_identified_infection),
        ("random_drug_cessation_probability", g.random_drug_cessation_probability),
        ("random_drug_cessation_probability_if_no_active_infection", g.random_drug_cessation_probability_if_no_active_infection),
        ("antibiotic_infection_prevention_efficacy", g.antibiotic_infection_prevention_efficacy),
    ]);

    print_scalar_group("Drug Selection", &[
        ("minimal_potency_threshold_for_drug_selection", g.minimal_potency_threshold_for_drug_selection),
        ("drug_selection_temperature", g.drug_selection_temperature),
        ("reserve_drug_score_penalty", g.reserve_drug_score_penalty),
    ]);

    print_scalar_group("Treatment Failure and Restart", &[
        ("treatment_failure_assessment_day", g.treatment_failure_assessment_day as f64),
        ("treatment_failure_threshold", g.treatment_failure_threshold),
        ("drug_failure_memory_days", g.drug_failure_memory_days as f64),
        ("restart_window_days", g.restart_window_days as f64),
        ("restart_bacteria_level_threshold", g.restart_bacteria_level_threshold),
        ("restart_window_probability", g.restart_window_probability),
    ]);

    print_scalar_group("Hospitalization", &[
        ("hospitalization_base_log_odds", g.hospitalization_base_log_odds),
        ("hospitalization_log_odds_per_age_year", g.hospitalization_log_odds_per_age_year),
        ("hospitalization_log_odds_sepsis", g.hospitalization_log_odds_sepsis),
        ("hospitalization_log_odds_symptomatic_infection", g.hospitalization_log_odds_symptomatic_infection),
        ("hospitalization_log_odds_serious_resistance_test_positive", g.hospitalization_log_odds_serious_resistance_test_positive),
        ("hospitalization_symptomatic_infection_level_threshold", g.hospitalization_symptomatic_infection_level_threshold),
        ("hospital_recovery_rate_per_day", g.hospital_recovery_rate_per_day),
        ("hospital_max_days", g.hospital_max_days),
        ("hospital_prevent_discharge_with_sepsis", g.hospital_prevent_discharge_with_sepsis),
    ]);

    print_scalar_group("Resistance Emergence and Decay", &[
        ("max_resistance_level", g.max_resistance_level),
        ("resistance_emergence_bacteria_level_multiplier", g.resistance_emergence_bacteria_level_multiplier),
        ("any_r_emergence_level_on_first_emergence", g.any_r_emergence_level_on_first_emergence),
        ("multi_drug_penalty_threshold_num_drugs", g.multi_drug_penalty_threshold_num_drugs),
        ("resistance_development_inhibition_single_drug", g.resistance_development_inhibition_single_drug),
        ("resistance_development_inhibition_partial_cross", g.resistance_development_inhibition_partial_cross),
        ("mechanism_assignment_probability_on_any_r_gain", g.mechanism_assignment_probability_on_any_r_gain),
        ("community_profile_cache_retention", g.community_profile_cache_retention),
        ("mechanism_reversion_rate_global_multiplier", g.mechanism_reversion_rate_global_multiplier),
        ("majority_r_memory_retention_per_day", g.majority_r_memory_retention_per_day),
    ]);

    print_scalar_group("Microbiome Dynamics", &[
        ("microbiome_resistance_transfer_probability_per_day", g.microbiome_resistance_transfer_probability_per_day),
        ("antibiotic_disruption_decay_half_life_days", g.antibiotic_disruption_decay_half_life_days),
        ("microbiome_resistance_multiplier_on_acquisition", g.microbiome_resistance_multiplier_on_acquisition),
        ("infection_from_microbiome_dampening", g.infection_from_microbiome_dampening),
        ("carriage_duration_log_odds_coefficient", g.carriage_duration_log_odds_coefficient),
        ("carriage_duration_max_log_odds_effect", g.carriage_duration_max_log_odds_effect),
        ("antibiotic_clearance_log_odds_per_unit_activity", g.antibiotic_clearance_log_odds_per_unit_activity),
        ("carrier_resistance_inheritance_probability", g.carrier_resistance_inheritance_probability),
        ("community_resistance_dilution_factor", g.community_resistance_dilution_factor),
        ("microbiome_majority_decay_half_life_days", g.microbiome_majority_decay_half_life_days),
        ("microbiome_minority_decay_half_life_days", g.microbiome_minority_decay_half_life_days),
        ("microbiome_majority_promotion_rate_per_day", g.microbiome_majority_promotion_rate_per_day),
    ]);

    print_scalar_group("De Novo and HGT Multipliers", &[
        ("infection_de_novo_multiplier", g.infection_de_novo_multiplier),
        ("microbiome_de_novo_multiplier", g.microbiome_de_novo_multiplier),
        ("hgt_multiplier", g.hgt_multiplier),
    ]);

    print_scalar_group("Horizontal Gene Transfer Modifiers", &[
        ("hgt_hospital_multiplier", g.hgt_hospital_multiplier),
        ("hgt_antibiotic_pressure_multiplier", g.hgt_antibiotic_pressure_multiplier),
        ("hgt_coinfection_multiplier", g.hgt_coinfection_multiplier),
        ("hgt_microbiome_only_penalty", g.hgt_microbiome_only_penalty),
        ("hgt_gut_compartment_multiplier", g.hgt_gut_compartment_multiplier),
        ("hgt_minority_donor_multiplier", g.hgt_minority_donor_multiplier),
    ]);

    print_scalar_group("Travel", &[
        ("travel_probability_per_day", g.travel_probability_per_day),
    ]);

    print_scalar_group("Bacteria Growth Age Multipliers", &[
        ("bacteria_growth_age_multiplier_infant", g.bacteria_growth_age_multiplier_infant),
        ("bacteria_growth_age_multiplier_child", g.bacteria_growth_age_multiplier_child),
        ("bacteria_growth_age_multiplier_adult", g.bacteria_growth_age_multiplier_adult),
        ("bacteria_growth_age_multiplier_elderly", g.bacteria_growth_age_multiplier_elderly),
        ("bacteria_growth_immunodeficiency_multiplier", g.bacteria_growth_immunodeficiency_multiplier),
    ]);

    print_scalar_group("Sepsis Onset", &[
        ("sepsis_minimum_duration_days", g.sepsis_minimum_duration_days as f64),
        ("log_odds_sepsis_onset_immunosuppressed", g.log_odds_sepsis_onset_immunosuppressed),
        ("log_odds_sepsis_onset_hospitalized", g.log_odds_sepsis_onset_hospitalized),
        ("log_odds_sepsis_onset_not_under_care", g.log_odds_sepsis_onset_not_under_care),
        ("log_odds_sepsis_onset_region_north_america", g.log_odds_sepsis_onset_region_north_america),
        ("log_odds_sepsis_onset_region_europe", g.log_odds_sepsis_onset_region_europe),
        ("log_odds_sepsis_onset_region_oceania", g.log_odds_sepsis_onset_region_oceania),
        ("log_odds_sepsis_onset_region_asia", g.log_odds_sepsis_onset_region_asia),
        ("log_odds_sepsis_onset_region_south_america", g.log_odds_sepsis_onset_region_south_america),
        ("log_odds_sepsis_onset_region_africa", g.log_odds_sepsis_onset_region_africa),
    ]);

    print_scalar_group("Sepsis Recovery", &[
        ("sepsis_base_log_odds_of_recovery_per_day", g.sepsis_base_log_odds_of_recovery_per_day),
        ("sepsis_log_odds_bacteria_level", g.sepsis_log_odds_bacteria_level),
        ("sepsis_log_odds_in_hospital", g.sepsis_log_odds_in_hospital),
        ("sepsis_log_odds_age_infant", g.sepsis_log_odds_age_infant),
        ("sepsis_log_odds_age_child", g.sepsis_log_odds_age_child),
        ("sepsis_log_odds_age_adult", g.sepsis_log_odds_age_adult),
        ("sepsis_log_odds_age_elderly", g.sepsis_log_odds_age_elderly),
        ("sepsis_log_odds_immunosuppressed", g.sepsis_log_odds_immunosuppressed),
    ]);

    print_scalar_group("Sepsis Death", &[
        ("sepsis_death_base_log_odds", g.sepsis_death_base_log_odds),
        ("sepsis_death_log_odds_age_infant", g.sepsis_death_log_odds_age_infant),
        ("sepsis_death_log_odds_age_child", g.sepsis_death_log_odds_age_child),
        ("sepsis_death_log_odds_age_adult", g.sepsis_death_log_odds_age_adult),
        ("sepsis_death_log_odds_age_elderly", g.sepsis_death_log_odds_age_elderly),
        ("sepsis_death_log_odds_immunosuppressed", g.sepsis_death_log_odds_immunosuppressed),
        ("sepsis_death_log_odds_bacteria_level", g.sepsis_death_log_odds_bacteria_level),
        ("sepsis_death_log_odds_duration", g.sepsis_death_log_odds_duration),
        ("sepsis_death_log_odds_early_phase", g.sepsis_death_log_odds_early_phase),
        ("sepsis_death_early_phase_days", g.sepsis_death_early_phase_days),
        ("sepsis_death_log_odds_not_under_care", g.sepsis_death_log_odds_not_under_care),
    ]);

    print_scalar_group("Non-Sepsis Infection Mortality", &[
        ("infection_non_sepsis_base_log_odds", g.infection_non_sepsis_base_log_odds),
        ("infection_non_sepsis_log_odds_per_level", g.infection_non_sepsis_log_odds_per_level),
        ("infection_non_sepsis_log_odds_age_infant", g.infection_non_sepsis_log_odds_age_infant),
        ("infection_non_sepsis_log_odds_age_child", g.infection_non_sepsis_log_odds_age_child),
        ("infection_non_sepsis_log_odds_age_adult", g.infection_non_sepsis_log_odds_age_adult),
        ("infection_non_sepsis_log_odds_age_elderly", g.infection_non_sepsis_log_odds_age_elderly),
        ("infection_non_sepsis_log_odds_immunosuppressed", g.infection_non_sepsis_log_odds_immunosuppressed),
        ("infection_non_sepsis_log_odds_in_hospital", g.infection_non_sepsis_log_odds_in_hospital),
        ("infection_non_sepsis_minimum_bacteria_level", g.infection_non_sepsis_minimum_bacteria_level),
    ]);

    print_scalar_group("Background Mortality", &[
        ("background_mortality_baseline_log_odds", g.background_mortality_baseline_log_odds),
        ("mortality_baseline_1930_multiplier", g.mortality_baseline_1930_multiplier),
        ("mortality_baseline_2035_multiplier", g.mortality_baseline_2035_multiplier),
        ("mortality_improvement_half_life_years", g.mortality_improvement_half_life_years),
        ("log_odds_mortality_per_year_of_age", g.log_odds_mortality_per_year_of_age),
        ("log_odds_mortality_per_year_of_age_squared", g.log_odds_mortality_per_year_of_age_squared),
        ("log_odds_mortality_immunosuppressed", g.log_odds_mortality_immunosuppressed),
        ("log_odds_mortality_hospitalized", g.log_odds_mortality_hospitalized),
    ]);

    print_scalar_group("Drug Toxicity", &[
        ("default_toxicity_reservoir_half_life_days", g.default_toxicity_reservoir_half_life_days),
        ("toxicity_age_multiplier_infant", g.toxicity_age_multiplier_infant),
        ("toxicity_age_multiplier_child", g.toxicity_age_multiplier_child),
        ("toxicity_age_multiplier_adult", g.toxicity_age_multiplier_adult),
        ("toxicity_age_multiplier_elderly", g.toxicity_age_multiplier_elderly),
        ("toxicity_immunosuppressed_multiplier", g.toxicity_immunosuppressed_multiplier),
        ("toxicity_hospital_multiplier", g.toxicity_hospital_multiplier),
        ("toxicity_discontinuation_threshold", g.toxicity_discontinuation_threshold),
        ("toxicity_discontinuation_avoidance_days", g.toxicity_discontinuation_avoidance_days as f64),
    ]);

    print_scalar_group("Regional Resistance Scoring", &[
        ("regional_resistance_threshold_very_high", g.regional_resistance_threshold_very_high),
        ("regional_resistance_threshold_high", g.regional_resistance_threshold_high),
        ("regional_resistance_threshold_moderate", g.regional_resistance_threshold_moderate),
        ("regional_resistance_penalty_very_high", g.regional_resistance_penalty_very_high),
        ("regional_resistance_penalty_high", g.regional_resistance_penalty_high),
        ("regional_resistance_penalty_moderate", g.regional_resistance_penalty_moderate),
    ]);

    print_scalar_group("Therapy Scoring", &[
        ("targeted_therapy_narrow_spectrum_bonus", g.targeted_therapy_narrow_spectrum_bonus),
        ("targeted_therapy_broad_spectrum_penalty", g.targeted_therapy_broad_spectrum_penalty),
        ("targeted_therapy_ineffective_drug_penalty", g.targeted_therapy_ineffective_drug_penalty),
        ("effective_potency_threshold_for_targeted_therapy", g.effective_potency_threshold_for_targeted_therapy),
        ("empiric_therapy_broad_spectrum_bonus", g.empiric_therapy_broad_spectrum_bonus),
        ("empiric_therapy_ineffective_penalty", g.empiric_therapy_ineffective_penalty),
    ]);

    print_scalar_group("MDR-TB Era Multipliers", &[
        ("mdr_tb_pre_antibiotic_era_multiplier", g.mdr_tb_pre_antibiotic_era_multiplier),
        ("mdr_tb_early_antibiotic_era_multiplier", g.mdr_tb_early_antibiotic_era_multiplier),
        ("mdr_tb_modern_era_multiplier", g.mdr_tb_modern_era_multiplier),
    ]);
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
    println!("Pharmacokinetic and clinical properties for each of the {} modelled \
              antimicrobial agents. The introduction time step is measured in days \
              from 1 January 1930.", DRUG_SHORT_NAMES.len());
    println!();
    println!("See: \
              [§6.3 Drug pharmacokinetics](#63-drug-pharmacokinetics), \
              [§6.5 Drug potency matrix](#65-drug-potency-matrix), \
              [§6.6 Drug availability](#66-drug-availability-by-region-and-era), \
              [§6.7 Drug toxicity](#67-drug-toxicity), \
              [§6.8 Antibiotic infection prevention](#68-antibiotic-infection-prevention).");
    println!();

    let headers = &[
        "Drug", "Class", "Intro (days)", "Init level",
        "t½ (days)", "2× dose mult", "Spectrum",
        "Tox hazard", "Tox t½ (days)", "Microbiome disrupt",
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
}

// ─────────────────────────────────────────────────────────────────────────────
// B.3  Bacteria Properties
// ─────────────────────────────────────────────────────────────────────────────

fn print_bacteria_properties(store: &amr_project::config::ParameterStore) {
    let b = &store.bacteria;
    println!("### B.3 Bacteria Properties");
    println!();
    println!("Per-bacteria parameters governing acquisition, growth, symptom \
              onset, and clinical outcomes for each of the {} bacterial species.", BACTERIA_LIST.len());
    println!();
    println!("See: \
              [§3.1 Community acquisition](#31-community-acquisition), \
              [§4.2 Infection dynamics](#42-infection-dynamics), \
              [§4.3 Sepsis](#43-sepsis), \
              [§4.4 Natural clearance](#44-natural-clearance), \
              [§8.1 Carriage compartments](#81-carriage-compartments).");
    println!();

    let headers = &[
        "Bacteria", "Acq log-odds", "Init level", "Δ level/day",
        "Max level", "Microb clr/day", "Microb vs inf",
        "Drug cess prob", "Sx threshold", "Sx delay (d)",
        "Sepsis log-odds", "Mech-less rev rate", "Comm dilution",
        "Hosp conc", "Hosp prune %",
    ];
    let mut rows = Vec::new();
    for (idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        rows.push(vec![
            bacteria.to_string(),
            format_value(b.acquisition_log_odds_baseline[idx]),
            format_value(b.initial_infection_level[idx]),
            format_value(b.base_bacteria_level_change[idx]),
            format_value(b.max_level[idx]),
            format_value(b.microbiome_clearance_probability_per_day[idx]),
            format_value(b.microbiome_vs_infection_log_odds[idx]),
            format_value(b.drug_cessation_probability[idx]),
            format_value(b.symptom_onset_threshold_level[idx]),
            format_value(b.symptom_onset_delay_days[idx]),
            format_value(b.sepsis_baseline_log_odds[idx]),
            format_value(b.mechanismless_resistance_reversion_rate[idx]),
            format_value(b.community_resistance_dilution_factor[idx]),
            format_value(b.hospital_resistance_concentration_factor[idx]),
            format_value(b.hospital_resistance_prune_susceptible_percent[idx]),
        ]);
    }
    md_table(headers, &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.4  Drug–Bacteria Potency Matrix
// ─────────────────────────────────────────────────────────────────────────────

fn print_drug_bacteria_matrix(store: &amr_project::config::ParameterStore) {
    println!("### B.4 Drug–Bacteria Potency Matrix");
    println!();
    println!("Baseline potency (MIC-derived effectiveness when no resistance is \
              present) and initiation multiplier (stewardship weighting for drug \
              selection) for each drug–bacteria pair. {} bacteria × {} drugs = {} entries.",
             BACTERIA_LIST.len(), DRUG_SHORT_NAMES.len(),
             BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len());
    println!();
    println!("See: \
              [§6.5 Drug potency matrix](#65-drug-potency-matrix), \
              [§6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use).");
    println!();

    let headers = &["Bacteria", "Drug", "Potency (no R)", "Init multiplier"];
    let mut rows = Vec::new();
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        for (d_idx, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            let potency = store.drug_bacteria.potency(b_idx, d_idx);
            let init_mult = store.drug_bacteria.initiation_multiplier(b_idx, d_idx);
            if potency.abs() > 1e-12 || (init_mult - 1.0).abs() > 1e-12 {
                rows.push(vec![
                    bacteria.to_string(),
                    drug.to_string(),
                    format_value(potency),
                    format_value(init_mult),
                ]);
            }
        }
    }
    md_table(headers, &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.5  Regional Parameters
// ─────────────────────────────────────────────────────────────────────────────

fn print_regional_parameters(store: &amr_project::config::ParameterStore) {
    println!("### B.5 Regional Parameters");
    println!();
    println!("Region-level scalars (applicable to all bacteria) and the per-region \
              per-bacteria acquisition log-odds adjustments.");
    println!();
    println!("See: \
              [§2.5 Travel](#25-travel), \
              [§3.1 Community acquisition](#31-community-acquisition).");
    println!();

    // Region scalars
    println!("#### Region Scalars");
    println!();
    let headers = &[
        "Region", "Travel mult", "Cessation mult", "Mortality log-odds",
        "Sepsis log-odds", "Sepsis mort mult", "Testing mult",
        "Abx init log-odds", "Hosp log-odds",
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
            format_value(store.region.sepsis_log_odds(region)),
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
    println!("Log-odds adjustments by age category for bacteria acquisition and \
              regional effects. Age categories: {}.",
             AGE_CATEGORY_SEQUENCE.iter()
                 .map(|c| c.label())
                 .collect::<Vec<_>>()
                 .join(", "));
    println!();
    println!("See: \
              [§2.2 Ageing and age categories](#22-ageing-and-age-categories), \
              [§3.1 Community acquisition](#31-community-acquisition).");
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
}

// ─────────────────────────────────────────────────────────────────────────────
// B.7  Syndrome Parameters
// ─────────────────────────────────────────────────────────────────────────────

fn print_syndrome_parameters(store: &amr_project::config::ParameterStore) {
    println!("### B.7 Syndrome Parameters");
    println!();
    println!("Infection-site (syndrome) specific parameters. Syndromes are: \
              1 = UTI, 2 = skin/soft tissue, 3 = respiratory, 4 = bloodstream, \
              5 = intra-abdominal, 6 = CNS/meningitis, 7 = gastrointestinal, \
              8 = genital/STI, 9 = bone/joint, 10 = other.");
    println!();
    println!("See: \
              [§4.1 Syndrome assignment](#41-syndrome-assignment), \
              [§6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use), \
              [§6.4 Drug penetration by syndrome](#64-drug-penetration-by-syndrome).");
    println!();

    // Empiric drug scores
    println!("#### Syndrome Empiric Drug Scores");
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

    // Drug penetration
    println!("#### Syndrome Drug Penetration");
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
    println!("Infection clearance model parameters. The clearance hazard is a \
              logistic function of base log-odds, per-bacteria adjustments, \
              age effects, immunodeficiency, bacteria level, and treatment duration.");
    println!();
    println!("See: [§4.4 Natural clearance](#44-natural-clearance).");
    println!();

    let rows = vec![
        vec!["base_clearance_log_odds".to_string(), format_value(store.clearance.base_log_odds())],
        vec!["immunodeficient_log_odds_adjustment".to_string(), format_value(store.clearance.immunodeficient_log_odds_adjustment())],
    ];
    md_table(&["Parameter", "Value"], &rows);

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
}

// ─────────────────────────────────────────────────────────────────────────────
// B.9  Immunodeficiency, Sex & Vaccination
// ─────────────────────────────────────────────────────────────────────────────

fn print_immunodeficiency_sex_vaccination(store: &amr_project::config::ParameterStore) {
    println!("### B.9 Immunodeficiency, Sex, and Vaccination Parameters");
    println!();
    println!("See: \
              [§2.3 Immunodeficiency](#23-immunodeficiency), \
              [§10 Mortality](#10-mortality).");
    println!();

    // Immunodeficiency
    println!("#### Immunodeficiency");
    println!();
    let mut rows = vec![
        vec!["startup_seed_fraction".to_string(), format_value(store.immunodeficiency.startup_seed_fraction())],
        vec!["temporary_onset_rate_per_day".to_string(), format_value(store.immunodeficiency.temporary_onset_rate())],
        vec!["temporary_recovery_rate_per_day".to_string(), format_value(store.immunodeficiency.temporary_recovery_rate())],
        vec!["chronic_onset_rate_per_day".to_string(), format_value(store.immunodeficiency.chronic_onset_rate())],
        vec!["chronic_recovery_rate_per_day".to_string(), format_value(store.immunodeficiency.chronic_recovery_rate())],
    ];
    for &(label, age_days) in &[("age_0_1", 180), ("age_1_18", 3650), ("age_18_65", 14600), ("age_65_plus", 25550)] {
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
        vec!["male".to_string(), format_value(store.sex.mortality_log_odds("male"))],
        vec!["female".to_string(), format_value(store.sex.mortality_log_odds("female"))],
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
    println!("Parameters for the {} resistance mechanisms modelled. Each mechanism \
              has a per-day reversion rate, per-drug-class enhancement multipliers, \
              and per-bacteria emergence rates.", ResistanceMechanism::all().len());
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
    println!("How much resistance each mechanism confers against each drug class. \
              Only non-zero entries shown.");
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
    md_table(&["Mechanism", "Drug class", "Enhancement multiplier"], &rows);

    // Bacteria-mechanism emergence rates
    println!("#### Bacteria–Mechanism Emergence Rates");
    println!();
    println!("De novo emergence rate per day for each bacteria–mechanism pair. \
              Only non-zero entries shown.");
    println!();
    let mut rows = Vec::new();
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        for (m_idx, mechanism) in mechanisms.iter().enumerate() {
            let rate = store.bacteria_mechanism_emergence.rate(b_idx, m_idx);
            if rate.abs() > 1e-20 {
                rows.push(vec![
                    bacteria.to_string(),
                    mechanism.as_str().to_string(),
                    format_value(rate),
                ]);
            }
        }
    }
    md_table(&["Bacteria", "Mechanism", "Emergence rate/day"], &rows);
}

// ─────────────────────────────────────────────────────────────────────────────
// B.11  Horizontal Gene Transfer Matrix
// ─────────────────────────────────────────────────────────────────────────────

fn print_hgt_matrix(store: &amr_project::config::ParameterStore) {
    println!("### B.11 Horizontal Gene Transfer Matrix");
    println!();
    println!("Per-day probability of horizontal gene transfer of resistance \
              between co-colonising bacterial species. Only non-zero entries shown.");
    println!();
    println!("See: \
              [§9.1 Transfer compatibility](#91-transfer-compatibility), \
              [§9.2 The HGT process](#92-the-hgt-process).");
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
