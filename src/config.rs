// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES}; // Import both lists

// --- Global Simulation Parameters ---
lazy_static! {
    pub static ref PARAMETERS: HashMap<String, f64> = {
        let mut map = HashMap::new();

        // General Drug Parameters
        // REMOVED: map.insert("drug_initial_level".to_string(), 10.0); // This general param is now replaced by drug-specific ones below
        map.insert("drug_base_initiation_rate_per_day".to_string(), 0.0001); // 0.0001
        map.insert("drug_infection_present_multiplier".to_string(), 50.0);
        map.insert("drug_test_identified_multiplier".to_string(), 50.0);
        map.insert("drug_decay_rate_per_day".to_string(), 1.0);
        // ADDED: New drug-related parameters
        map.insert("already_on_drug_initiation_multiplier".to_string(), 0.000); // 0.0001
        map.insert("double_dose_probability_if_identified_infection".to_string(), 0.1); // NEW: Probability for double dose

        map.insert("random_drug_cessation_probability".to_string(), 0.001); // New: Probability an individual randomly stops a drug per day

        // General Acquisition & Resistance Parameters
        // this two below will need to change over calendar time - for the hospital acquired may decide to sample from 
        // majority_r of people in hospital with the bacteria  
        map.insert("environmental_majority_r_level_for_new_acquisition".to_string(), 0.0);
        map.insert("hospital_majority_r_level_for_new_acquisition".to_string(), 0.0);

        map.insert("max_resistance_level".to_string(), 1.0);
        map.insert("majority_r_evolution_rate_per_day_when_drug_present".to_string(), 0.001);

        // Resistance Emergence and Decay Parameters
        map.insert("resistance_emergence_rate_per_day_baseline".to_string(), 0.000001); // Baseline probability for de novo resistance emergence
        map.insert("resistance_emergence_bacteria_level_multiplier".to_string(), 0.05); // Multiplier for bacteria level's effect on emergence
        map.insert("any_r_emergence_level_on_first_emergence".to_string(), 0.5); // The resistance level 'any_r' starts at upon emergence
        map.insert("any_r_decay_rate_per_day".to_string(), 0.5); // Rate at which 'any_r' decays when 'majority_r' is 0

        // Testing Parameters
        map.insert("test_delay_days".to_string(), 3.0);
        map.insert("test_rate_per_day".to_string(), 0.20);  // 0.15

        // Syndrome-specific multipliers (example)
        map.insert("syndrome_3_initiation_multiplier".to_string(), 10.0); // Respiratory syndrome
        map.insert("syndrome_7_initiation_multiplier".to_string(), 8.0);  // Gastrointestinal syndrome
        map.insert("syndrome_8_initiation_multiplier".to_string(), 12.0); // Genital syndrome (example ID)
        
        // --- Default Parameters for ALL Bacteria from BACTERIA_LIST ---
        // These are inserted first, and can then be overridden by specific entries below.
        for &bacteria in BACTERIA_LIST.iter() {
            map.insert(format!("{}_acquisition_prob_baseline", bacteria), 0.0); // 0.01
            map.insert(format!("{}_initial_infection_level", bacteria), 0.01);
            map.insert(format!("{}_environmental_acquisition_proportion", bacteria), 0.1);
            map.insert(format!("{}_hospital_acquired_proportion", bacteria), 0.05);
            map.insert(format!("{}_decay_rate", bacteria), 0.02);
            map.insert(format!("{}_adult_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_child_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_oral_exposure_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_sexual_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_mosquito_exposure_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_vaccine_efficacy", bacteria), 0.0); // Default to no vaccine effect
            map.insert(format!("{}_level_change_rate_baseline", bacteria), 0.2); // Default to no growth/decay
            map.insert(format!("{}_immunity_effect_on_level_change", bacteria), 0.01);
            map.insert(format!("{}_max_level", bacteria), 100.0);
            map.insert(format!("{}_immunity_increase_rate_baseline", bacteria), 0.001);
            map.insert(format!("{}_initial_immunity_on_infection", bacteria), 0.0001);
            map.insert(format!("{}_immunity_increase_rate_per_level", bacteria), 0.05);
            map.insert(format!("{}_immunity_age_modifier", bacteria), 1.0);
            map.insert(format!("{}_baseline_immunity_level", bacteria), 0.00001);
            map.insert(format!("{}_immunity_decay_rate", bacteria), 0.1);
        }

        // ADDED: Default Initial Drug Levels and Double Dose Multipliers for ALL Drugs
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_initial_level", drug), 10.0); // Default initial level for each drug
            map.insert(format!("drug_{}_double_dose_multiplier", drug), 2.0); // Default double dose multiplier
        }

        // --- Overrides for Specific Bacteria (Customize these as needed) ---

        // acinetobac_bau Parameters
        map.insert("acinetobac_bau_acquisition_prob_baseline".to_string(), 0.2);
        map.insert("acinetobac_bau_hospital_acquired_proportion".to_string(), 0.15); // Often hospital-acquired
        map.insert("acinetobac_bau_immunity_increase_rate_baseline".to_string(), 0.0); // 0.001
        map.insert("acinetobac_bau_immunity_increase_rate_per_day".to_string(), 0.0);  // 0.2
        map.insert("acinetobac_bau_immunity_increase_rate_per_level".to_string(), 0.0);  // 0.2
        map.insert("acinetobac_bau_immunity_effect_on_level_change".to_string(), 0.005);  
        map.insert("acinetobac_bau_resistance_emergence_rate_per_day_baseline".to_string(), 0.7); // Baseline probability for de novo resistance emergence
 

        // Add more specific parameters for acinetobac_bau if needed

        // ADDED: Overrides for Specific Drug Initial Levels & Double Dose Multipliers
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
    // MODIFIED: Key changed from (String, String) to (String, String, String) to include param_suffix
    pub static ref BACTERIA_DRUG_PARAMETERS: HashMap<(String, String, String), f64> = {
        let mut map = HashMap::new();

        // Default antibiotic reduction for ALL bacteria-drug combinations
        for &bacteria in BACTERIA_LIST.iter() {
            for &drug in DRUG_SHORT_NAMES.iter() {
                // MODIFIED: Key now includes the parameter suffix
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

// ADDED: This function was missing from your config.rs
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
    // MODIFIED: Key now includes the param_suffix directly, no need for the if check
    let specific_key = (bacteria_name.to_string(), drug_name.to_string(), param_suffix.to_string());
    BACTERIA_DRUG_PARAMETERS.get(&specific_key).copied()
}