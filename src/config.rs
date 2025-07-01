// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES}; // Import both lists

// --- Global Simulation Parameters ---
lazy_static! {
    pub static ref PARAMETERS: HashMap<String, f64> = {
        let mut map = HashMap::new();

        // General Drug Parameters
        map.insert("drug_base_initiation_rate_per_day".to_string(), 0.0000001); // 0.0005
        map.insert("drug_infection_present_multiplier".to_string(), 50.0);
        map.insert("drug_test_identified_multiplier".to_string(), 50.0);
        map.insert("drug_decay_rate_per_day".to_string(), 1.0);
        map.insert("already_on_drug_initiation_multiplier".to_string(), 0.000); // 0.0001
        map.insert("double_dose_probability_if_identified_infection".to_string(), 0.1); // Probability for double dose
      

        for &drug in DRUG_SHORT_NAMES.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                // Default to 1.0 (no multiplier effect) for all combinations
                map.insert(format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug, bacteria), 1.0);
            }
        }


        // todo: for each drug-bacteria combination will need a specific multiplier for initiation rate
        // will need changes also in mod.rs 

        map.insert("random_drug_cessation_probability".to_string(), 0.03); // Probability an individual randomly stops a drug per day

        // General Acquisition & Resistance Parameters
        // this two below will need to change over calendar time - for the hospital acquired may decide to sample from 
        // majority_r of people in hospital with the bacteria  
        map.insert("environmental_majority_r_level_for_new_acquisition".to_string(), 0.0);
        map.insert("hospital_majority_r_level_for_new_acquisition".to_string(), 0.0);

        map.insert("max_resistance_level".to_string(), 1.0);
        map.insert("majority_r_evolution_rate_per_day_when_drug_present".to_string(), 0.001);

        // Resistance Emergence and Decay Parameters
        map.insert("resistance_emergence_rate_per_day_baseline".to_string(), 0.01);  // 0.000001 Baseline probability for de novo resistance emergence
        map.insert("resistance_emergence_bacteria_level_multiplier".to_string(), 0.05); // Multiplier for bacteria level's effect on emergence
        map.insert("any_r_emergence_level_on_first_emergence".to_string(), 0.5); // The resistance level 'any_r' starts at upon emergence

        
        //  Microbiome Resistance Transfer Parameter
        map.insert("microbiome_resistance_transfer_probability_per_day".to_string(), 0.05); // Probability per day for resistance transfer between infection and microbiome
    

        // Testing Parameters
        map.insert("test_delay_days".to_string(), 3.0);
        map.insert("test_rate_per_day".to_string(), 0.20);  // 0.15

        // Syndrome-specific multipliers (example)
        map.insert("syndrome_3_initiation_multiplier".to_string(), 10.0); // Respiratory syndrome
        map.insert("syndrome_7_initiation_multiplier".to_string(), 8.0);  // Gastrointestinal syndrome
        map.insert("syndrome_8_initiation_multiplier".to_string(), 12.0); // Genital syndrome (example ID)        

        // Hospitalization Parameters
        map.insert("hospitalization_baseline_rate_per_day".to_string(), 0.00001); // Baseline daily probability of hospitalization
        map.insert("hospitalization_age_multiplier_per_day".to_string(), 0.000001); // Increase in daily hospitalization probability per year of age
        map.insert("hospitalization_recovery_rate_per_day".to_string(), 0.1); // Daily probability of recovering from hospitalization
        map.insert("hospitalization_max_days".to_string(), 30.0); // Max days in hospital before forced discharge (as fallback)

        // initiate travel
        map.insert("travel_probability_per_day".to_string(), 0.0005);

        // --- Default Parameters for ALL Bacteria from BACTERIA_LIST ---
        // These are inserted first, and can then be overridden by specific entries below.
        for &bacteria in BACTERIA_LIST.iter() {
            map.insert(format!("{}_acquisition_prob_baseline", bacteria), 0.0); // 0.01
            map.insert(format!("{}_initial_infection_level", bacteria), 0.01); // 0.01
            map.insert(format!("{}_environmental_acquisition_proportion", bacteria), 0.1); // 0.1
            map.insert(format!("{}_hospital_acquired_proportion", bacteria), 0.05); // 0.05
            map.insert(format!("{}_decay_rate", bacteria), 0.02);
            map.insert(format!("{}_adult_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_child_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_oral_exposure_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_sexual_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_mosquito_exposure_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_vaccine_efficacy", bacteria), 0.0); // Default to no vaccine effect
            map.insert(format!("{}_level_change_rate_baseline", bacteria), 0.2); // currently if this is non-zero it means that people can have level > 0 and ths appear infected
            map.insert(format!("{}_immunity_effect_on_level_change", bacteria), 0.01);
            map.insert(format!("{}_max_level", bacteria), 5.0);
            map.insert(format!("{}_immunity_increase_rate_baseline", bacteria), 0.1); // 0.001
            map.insert(format!("{}_initial_immunity_on_infection", bacteria), 0.001); 
            map.insert(format!("{}_immunity_increase_rate_per_level", bacteria), 0.05);
            map.insert(format!("{}_immunity_age_modifier", bacteria), 1.0);
            map.insert(format!("{}_baseline_immunity_level", bacteria), 0.00001); // 0.00001
            map.insert(format!("{}_immunity_decay_rate", bacteria), 0.1);
        }

        // Default Initial Drug Levels and Double Dose Multipliers for ALL Drugs
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_initial_level", drug), 10.0); // Default initial level for each drug
            map.insert(format!("drug_{}_double_dose_multiplier", drug), 2.0); // Default double dose multiplier
        }

        // Global defaults, used if a bacteria-specific parameter is not found
        map.insert("default_sepsis_baseline_risk_per_day".to_string(), 0.00001); // Very small baseline daily risk
        map.insert("default_sepsis_level_multiplier".to_string(), 0.005); // Multiplier for bacterial level (e.g., higher level = higher risk)
        map.insert("default_sepsis_duration_multiplier".to_string(), 0.000001); // Multiplier for duration of infection (e.g., longer duration = higher risk)


        // Background Mortality Parameters (Age, Region, and Sex dependent)
        map.insert("base_background_mortality_rate_per_day".to_string(), 0.000005); // Example: 0.0005% chance of death per day, for a baseline individual
        map.insert("age_mortality_multiplier_per_year".to_string(), 0.0000001); // Example: Small increase in daily death risk per year of age

        // Region-specific mortality multipliers. Ensure these match your `Region` enum variants.
        map.insert("northamerica_mortality_multiplier".to_string(), 1.0);
        map.insert("southamerica_mortality_multiplier".to_string(), 1.0);
        map.insert("africa_mortality_multiplier".to_string(), 1.2);    // Example: 20% higher mortality risk
        map.insert("asia_mortality_multiplier".to_string(), 1.1);
        map.insert("europe_mortality_multiplier".to_string(), 0.9);     // Example: 10% lower mortality risk
        map.insert("oceania_mortality_multiplier".to_string(), 1.0);
        map.insert("home_mortality_multiplier".to_string(), 1.0);       // If 'Home' is a specific region in your logic

        // Sex-specific mortality multipliers. Ensure these match your `sex_at_birth` strings.
        map.insert("male_mortality_multiplier".to_string(), 1.1);   // Example: Males have 10% higher mortality risk
        map.insert("female_mortality_multiplier".to_string(), 0.9); // Example: Females have 10% lower mortality risk


        //  Immunosuppression Onset and Recovery Rates
        map.insert("immunosuppression_onset_rate_per_day".to_string(), 0.0001);   // Probability of becoming immunosuppressed daily
        map.insert("immunosuppression_recovery_rate_per_day".to_string(), 0.0005); // Probability of recovering from immunosuppression daily


        // Sepsis Mortality Parameter (Absolute risk, independent of other factors)
        map.insert("sepsis_absolute_death_risk_per_day".to_string(), 0.1); // Example: 10% absolute chance of death per day if septic

        //  Default Toxicity Parameter
        map.insert("default_drug_toxicity_per_unit_level_per_day".to_string(), 0.005); // Adjust this default as needed

        //  Default Microbiome Acquisition Parameter
        // A multiplier for the infection acquisition probability to get microbiome acquisition probability.
        // If > 1.0, microbiome acquisition is more likely than infection for the same factors.
        // If < 1.0, microbiome acquisition is less likely.
        map.insert("default_microbiome_acquisition_multiplier".to_string(), 2.0); // Example: Microbiome acquisition is twice as likely as infection given the same exposure.

        //  Default Microbiome Clearance Parameter (from previous suggestion, ensure it's there)
        map.insert("default_microbiome_clearance_probability_per_day".to_string(), 0.01); // E.g., 1% chance to lose carriage per day

        //  Microbiome Presence Effect on Infection Acquisition
        // A multiplier for infection acquisition probability if the bacteria is already present in the microbiome.
        // Value > 1.0 means microbiome presence increases infection risk.
        // Value < 1.0 means microbiome presence decreases infection risk (e.g., due to local immunity/competition).
        map.insert("default_microbiome_infection_acquisition_multiplier".to_string(), 0.1); // Example: Much harder to get infected if already colonized.

        //  Contact and Exposure Level Parameters
        map.insert("contact_level_daily_fluctuation_range".to_string(), 0.5); // Amount of random daily fluctuation
        map.insert("min_contact_level".to_string(), 0.0); // Minimum possible contact/exposure level
        map.insert("max_contact_level".to_string(), 10.0); // Maximum possible contact/exposure level

        // Sexual Contact Parameters
        map.insert("sexual_contact_baseline".to_string(), 5.0); // Baseline level for a young adult
        map.insert("sexual_contact_age_peak_days".to_string(), 25.0 * 365.0); // Age in days (25 years)
        map.insert("sexual_contact_age_rise_exponent".to_string(), 2.0); // Controls how fast contact rises with age before peak (higher = steeper)
        map.insert("sexual_contact_age_decline_rate".to_string(), 0.00005); // Rate of decline per day after peak age (e.g., 0.00005 means ~1.8% drop per year)
        map.insert("sexual_contact_hospital_multiplier".to_string(), 0.0); // Significantly reduced in hospital (0.0 means no contact)

        // Airborne Contact (Adults) Parameters
        map.insert("airborne_contact_adult_baseline".to_string(), 5.0);
        map.insert("airborne_contact_adult_age_breakpoint_days".to_string(), 18.0 * 365.0); // Age in days (18 years)
        map.insert("airborne_contact_adult_child_multiplier".to_string(), 0.2); // How much less children contact adults (vs. adult-adult baseline)
        map.insert("airborne_contact_in_hospital_multiplier".to_string(), 1.5); // May increase due to healthcare staff contact

        // Airborne Contact (Children) Parameters
        map.insert("airborne_contact_child_baseline".to_string(), 3.0);
        map.insert("airborne_contact_child_age_breakpoint_days".to_string(), 12.0 * 365.0); // Age in days (12 years)
        map.insert("airborne_contact_child_child_multiplier".to_string(), 1.5); // How much more children contact children (vs. child baseline)
        map.insert("airborne_contact_child_adult_multiplier".to_string(), 0.5); // How much less adults contact children (vs. child baseline)

        // Oral Exposure Parameters
        map.insert("oral_exposure_baseline".to_string(), 2.0);
        map.insert("oral_exposure_child_age_breakpoint_days".to_string(), 5.0 * 365.0); // Age in days (5 years)
        map.insert("oral_exposure_child_multiplier".to_string(), 3.0); // Higher for young children
        map.insert("oral_exposure_in_hospital_multiplier".to_string(), 0.8); // Slightly reduced due to hospital hygiene

        // Mosquito Exposure Parameters
        map.insert("mosquito_exposure_baseline".to_string(), 1.0);
        map.insert("mosquito_exposure_in_hospital_multiplier".to_string(), 0.2); // Significantly reduced indoors/hospital
        
        // Region-specific multipliers (example values, adjust as needed based on actual epidemiology)
        map.insert("north_america_mosquito_exposure_multiplier".to_string(), 0.5);
        map.insert("south_america_mosquito_exposure_multiplier".to_string(), 5.0);
        map.insert("africa_mosquito_exposure_multiplier".to_string(), 8.0);
        map.insert("asia_mosquito_exposure_multiplier".to_string(), 6.0);
        map.insert("europe_mosquito_exposure_multiplier".to_string(), 0.2);
        map.insert("oceania_mosquito_exposure_multiplier".to_string(), 3.0);
        // Ensure you have multipliers for all variants of your `Region` enum,
        // or add a default handling in the `mod.rs` if a region param isn't found.
        // If `Region::Home` refers to a generic home location not tied to a specific geographical region,
        // you might need to reconsider its role or default it to 1.0 or an average.

 












        // Drug-specific Adverse Event Death Risk Parameters (Absolute risk, independent and drug-specific)
        // Add an entry for each drug that can cause an adverse event leading to death.
        // The key format is "drug_{drug_short_name}_adverse_event_death_risk"
        map.insert("drug_amoxicillin_adverse_event_death_risk".to_string(), 0.0001); // Example: 0.01% daily absolute risk from amoxicillin
        map.insert("drug_meropenem_adverse_event_death_risk".to_string(), 0.0005);   // Example: 0.05% daily absolute risk from meropenem
        // Add similar lines for other drugs as needed:
        // map.insert("drug_ciprofloxacin_adverse_event_death_risk".to_string(), 0.0002);
        // map.insert("drug_vancomycin_adverse_event_death_risk".to_string(), 0.0003);


        // --- Overrides for Specific Bacteria (Customize these as needed) ---

        // acinetobac_bau Parameters
        map.insert("acinetobac_bau_acquisition_prob_baseline".to_string(), 0.5);
        map.insert("acinetobac_bau_hospital_acquired_proportion".to_string(), 0.15); // Often hospital-acquired
        map.insert("acinetobac_bau_immunity_increase_rate_baseline".to_string(), 0.001); // 0.001
        map.insert("acinetobac_bau_immunity_increase_rate_per_day".to_string(), 0.2);  // 0.2
        map.insert("acinetobac_bau_immunity_increase_rate_per_level".to_string(), 0.2);  // 0.2
        map.insert("acinetobac_bau_immunity_effect_on_level_change".to_string(), 0.005);  
        map.insert("acinetobac_bau_resistance_emergence_rate_per_day_baseline".to_string(), 0.0); // 0.7 Baseline probability for de novo resistance emergence

        map.insert("acinetobac_bau_baseline_risk_per_day".to_string(), 0.00002);
        map.insert("acinetobac_bau_level_multiplier".to_string(), 0.008);
        map.insert("acinetobac_bau_duration_multiplier".to_string(), 0.000002);



        map.insert("drug_cefepime_for_bacteria_acinetobac_bau_initiation_multiplier".to_string(), 10000.0); // High multiplier



        //  Drug-Specific Toxicity Parameters (Examples)
        // Format: "drug_{drug_name}_toxicity_per_unit_level_per_day"
        map.insert("drug_penicilling_toxicity_per_unit_level_per_day".to_string(), 0.002); // Lower toxicity example
        map.insert("drug_cefepime_toxicity_per_unit_level_per_day".to_string(), 0.01);    // Higher toxicity example
        map.insert("drug_meropenem_toxicity_per_unit_level_per_day".to_string(), 0.008);


        //  Bacteria-Specific Microbiome Acquisition Multipliers (Optional Examples)
        // Format: "{bacteria_name}_microbiome_acquisition_multiplier"
        map.insert("strep_pneu_microbiome_acquisition_multiplier".to_string(), 3.0); // Strep Pneu might colonize more easily
        map.insert("salm_typhi_microbiome_acquisition_multiplier".to_string(), 0.5);  // Salmonella might colonize less easily than cause infection

        //  Bacteria-Specific Microbiome Clearance Parameters (Optional Examples)
        // Format: "{bacteria_name}_microbiome_clearance_probability_per_day"
        map.insert("strep_pneu_microbiome_clearance_probability_per_day".to_string(), 0.02);
        map.insert("esch_coli_microbiome_clearance_probability_per_day".to_string(), 0.005);


        
        //  Bacteria-Specific Microbiome Infection Acquisition Multipliers (Optional Examples)
        // Format: "{bacteria_name}_microbiome_infection_acquisition_multiplier"
        map.insert("strep_pneu_microbiome_infection_acquisition_multiplier".to_string(), 0.05); // Maybe Strep Pneu is very protective when colonized
        map.insert("salm_typhi_microbiome_infection_acquisition_multiplier".to_string(), 0.8);  // Salmonella might offer less protection, or even slightly increase risk if certain strains


        // Add specific rates for your drugs (e.g., higher for more toxic drugs)
        map.insert("gentamicin_toxicity_mortality_rate_per_day".to_string(), 0.00001); // Example for a specific drug

        // Add more specific parameters for acinetobac_bau if needed

        // Overrides for Specific Drug Initial Levels & Double Dose Multipliers
        map.insert("drug_penicilling_double_dose_multiplier".to_string(), 2.5); // Example: higher multiplier for penicillin

        map.insert("drug_amoxicillin_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_azithromycin_double_dose_multiplier".to_string(), 1.8);

        map.insert("drug_ciprofloxacin_double_dose_multiplier".to_string(), 2.2);

        map.insert("drug_trim_sulf_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_meropenem_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_cefepime_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_vancomycin_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_linezolid_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_ceftriaxone_double_dose_multiplier".to_string(), 2.0);

        /* ... (rest of your commented out bacteria parameters) ... */

        map
    };

    // --- Bacteria-Drug Specific Parameters ---
    // Key changed from (String, String) to (String, String, String) to include param_suffix
    pub static ref BACTERIA_DRUG_PARAMETERS: HashMap<(String, String, String), f64> = {
        let mut map = HashMap::new();

        // Default antibiotic reduction for ALL bacteria-drug combinations
        for &bacteria in BACTERIA_LIST.iter() {
            for &drug in DRUG_SHORT_NAMES.iter() {
                // Key now includes the parameter suffix
                map.insert((bacteria.to_string(), drug.to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.1); // 0.004 Generic reduction
            }
        }

        // Overrides for specific bacteria-drug combinations (Customize these as needed)
        // All these inserts below must also change to the 3-tuple key!

        // strep_pneu drug effectiveness
        map.insert(("strep_pneu".to_string(), "penicilling".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.1);
        map.insert(("strep_pneu".to_string(), "amoxicillin".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.1);
        map.insert(("strep_pneu".to_string(), "azithromycin".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.1);

        // haem_infl drug effectiveness
        map.insert(("haem_infl".to_string(), "amoxicillin".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.006);
        map.insert(("haem_infl".to_string(), "azithromycin".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.005);

        // salm_typhi drug effectiveness
        map.insert(("salm_typhi".to_string(), "ciprofloxacin".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.009);
        map.insert(("salm_typhi".to_string(), "azithromycin".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.007);

        // esch_coli drug effectiveness
        map.insert(("esch_coli".to_string(), "trim_sulf".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.005);
        map.insert(("esch_coli".to_string(), "ciprofloxacin".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.006);

        // pseud_aerug drug effectiveness (example)
        map.insert(("pseud_aerug".to_string(), "meropenem".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.01);
        map.insert(("pseud_aerug".to_string(), "cefepime".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.008);

        // staph_aureus drug effectiveness (example)
        map.insert(("staph_aureus".to_string(), "vancomycin".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.009);
        map.insert(("staph_aureus".to_string(), "linezolid".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.008);

        // n_gonorrhoeae drug effectiveness (example)
        map.insert(("n_gonorrhoeae".to_string(), "ceftriaxone".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.0095);

        // You can add more specific drug effectiveness for each bacteria as needed
        // For example:
        // map.insert(("acinetobac_bau".to_string(), "imipenem_c".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.009);
        // map.insert(("kleb_pneu".to_string(), "meropenem".to_string(), "bacteria_level_reduction_per_unit_of_drug".to_string()), 0.009);

        map
    };
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
    let specific_key = format!("{}_{}", bacteria_name, param_suffix);
    PARAMETERS.get(&specific_key).copied()
}
/// Retrieves a drug-specific simulation parameter.
/// It looks up "drug_{drug_name}_{param_suffix}".
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_drug_param(drug_name: &str, param_suffix: &str) -> Option<f64> {
    let specific_key = format!("drug_{}_{}", drug_name, param_suffix);
    PARAMETERS.get(&specific_key).copied()
}

/// Retrieves a bacteria-drug-specific simulation parameter.
/// It directly looks up the specific (bacteria_name, drug_name, param_suffix) tuple.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_bacteria_drug_param(bacteria_name: &str, drug_name: &str, param_suffix: &str) -> Option<f64> {
    let specific_key = (bacteria_name.to_string(), drug_name.to_string(), param_suffix.to_string());
    BACTERIA_DRUG_PARAMETERS.get(&specific_key).copied()
}