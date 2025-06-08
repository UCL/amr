// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;

// --- Global Simulation Parameters ---
lazy_static! {
    pub static ref PARAMETERS: HashMap<String, f64> = {
        let mut map = HashMap::new();

        // General Drug Parameters
        map.insert("drug_initial_level".to_string(), 10.0);
        map.insert("drug_base_initiation_rate_per_day".to_string(), 0.0001);
        map.insert("drug_infection_present_multiplier".to_string(), 50.0);
        map.insert("drug_test_identified_multiplier".to_string(), 20.0);
        map.insert("drug_decay_rate_per_day".to_string(), 0.3);

        // General Acquisition & Resistance Parameters
        map.insert("environmental_majority_r_level_for_new_acquisition".to_string(), 0.0);
        map.insert("hospital_majority_r_level_for_new_acquisition".to_string(), 0.0);
        map.insert("max_resistance_level".to_string(), 10.0);
        map.insert("majority_r_evolution_rate_per_day_when_drug_present".to_string(), 0.001); // Example value, adjust as needed

        // Testing Parameters
        map.insert("test_delay_days".to_string(), 3.0);
        map.insert("test_rate_per_day".to_string(), 0.15);

        // Syndrome-specific multipliers (example)
        // Format: "syndrome_{ID}_initiation_multiplier"
        map.insert("syndrome_3_initiation_multiplier".to_string(), 10.0); // Example: Respiratory syndrome might increase drug initiation

        // --- Bacteria-Specific Parameters ---
        // These are accessed via get_bacteria_param.
        // Naming convention: "{bacteria_name}_{param_suffix}"
        // A "generic_bacteria_" prefix can be used for parameters that apply to all bacteria
        // unless a specific bacteria override is provided.

        // Generic Bacteria Defaults (used if specific bacteria parameter is not found)
        map.insert("generic_bacteria_acquisition_prob_baseline".to_string(), 0.01);
        map.insert("generic_bacteria_initial_infection_level".to_string(), 0.01);
        map.insert("generic_environmental_acquisition_proportion".to_string(), 0.1);
        map.insert("generic_hospital_acquired_proportion".to_string(), 0.05);
        map.insert("generic_bacteria_decay_rate".to_string(), 0.02);
        // Add more generic defaults here as needed

        // Strep Pneumoniae Specific Parameters
        map.insert("strep_pneu_acquisition_prob_baseline".to_string(), 0.005); // Example: different from generic
        map.insert("strep_pneu_adult_contact_acq_rate_ratio_per_unit".to_string(), 1.2);
        map.insert("strep_pneu_child_contact_acq_rate_ratio_per_unit".to_string(), 1.5);
        map.insert("strep_pneu_vaccine_efficacy".to_string(), 0.8);
        map.insert("strep_pneu_initial_infection_level".to_string(), 0.05); // Example: different from generic
        map.insert("strep_pneu_environmental_acquisition_proportion".to_string(), 0.05); // Example: different from generic
        map.insert("strep_pneu_hospital_acquired_proportion".to_string(), 0.1); // Example: different from generic

        map.insert("strep_pneu_baseline_immunity_level".to_string(), 0.1);
        map.insert("strep_pneu_immunity_decay_rate".to_string(), 0.001);
        map.insert("strep_pneu_level_change_rate_baseline".to_string(), 0.005); // Can be positive for growth
        map.insert("strep_pneu_immunity_effect_on_level_change".to_string(), 0.01);
        map.insert("strep_pneu_max_level".to_string(), 100.0);
        map.insert("strep_pneu_immunity_increase_rate_baseline".to_string(), 0.005);
        map.insert("strep_pneu_immunity_increase_rate_per_day".to_string(), 0.0001);
        map.insert("strep_pneu_immunity_increase_rate_per_level".to_string(), 0.002);
        map.insert("strep_pneu_immunity_age_modifier".to_string(), 1.0); // Multiplier applied to immunity increase based on age

        // Add parameters for other bacteria (haem_infl, salm_typhi, esch_coli) here following the same pattern
        // Example for haem_infl (you would fill in actual values):
        map.insert("haem_infl_acquisition_prob_baseline".to_string(), 0.012);
        map.insert("haem_infl_initial_infection_level".to_string(), 0.015);
        map.insert("haem_infl_environmental_acquisition_proportion".to_string(), 0.08);
        map.insert("haem_infl_hospital_acquired_proportion".to_string(), 0.06);
        map.insert("haem_infl_decay_rate".to_string(), 0.025);


        map
    };

    // --- Bacteria-Drug Specific Parameters ---
    // This HashMap is for parameters that depend on both a specific bacteria and a specific drug.
    // The key is a tuple (bacteria_name, drug_name).
    // The param_suffix in get_bacteria_drug_param is currently simplified to just "antibiotic_reduction_per_unit".
    pub static ref BACTERIA_DRUG_PARAMETERS: HashMap<(String, String), f64> = {
        let mut map = HashMap::new();

        // Example: Strep Pneu and "generic_drug" (if you have a catch-all drug effect)
        map.insert(("strep_pneu".to_string(), "generic_drug".to_string()), 0.005); // antibiotic_reduction_per_unit for strep_pneu

        // Add specific bacteria-drug combinations if antibiotic_reduction_per_unit varies by drug
        // map.insert(("strep_pneu".to_string(), "amoxicillin".to_string()), 0.007);
        // map.insert(("strep_pneu".to_string(), "azithromycin".to_string()), 0.003);

        map
    };
}

/// Retrieves a global simulation parameter.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_global_param(key: &str) -> Option<f64> {
    PARAMETERS.get(key).copied()
}

/// Retrieves a bacteria-specific simulation parameter.
/// It first tries to find "{bacteria_name}_{param_suffix}".
/// If not found, it falls back to "generic_bacteria_{param_suffix}" or other generic forms.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_bacteria_param(bacteria_name: &str, param_suffix: &str) -> Option<f64> {
    // 1. Try bacteria-specific key (e.g., "strep_pneu_acquisition_prob_baseline")
    let specific_key = format!("{}_{}", bacteria_name, param_suffix);
    if let Some(&value) = PARAMETERS.get(&specific_key) {
        return Some(value);
    }

    // 2. Fallback to generic_bacteria_ prefix for common parameters if not found specifically
    let generic_key_prefix = "generic_bacteria_";
    let generic_key = format!("{}{}", generic_key_prefix, param_suffix);
    if let Some(&value) = PARAMETERS.get(&generic_key) {
        return Some(value);
    }

    // 3. Special fallback for specific "generic_" parameters that don't follow the "generic_bacteria_" prefix
    // (This handles cases like "generic_environmental_acquisition_proportion" as in your original code)
    match param_suffix {
        "acquisition_prob_baseline" => PARAMETERS.get("generic_bacteria_acquisition_prob_baseline").copied(),
        "initial_infection_level" => PARAMETERS.get("generic_bacteria_initial_infection_level").copied(),
        "environmental_acquisition_proportion" => PARAMETERS.get("generic_environmental_acquisition_proportion").copied(),
        "hospital_acquired_proportion" => PARAMETERS.get("generic_hospital_acquired_proportion").copied(),
        "decay_rate" => PARAMETERS.get("generic_bacteria_decay_rate").copied(),
        _ => None, // No generic fallback for this suffix
    }
}


/// Retrieves a bacteria-drug-specific simulation parameter.
/// Currently focused on "antibiotic_reduction_per_unit".
/// It first tries to find a specific (bacteria, drug) combination.
/// If not found, it tries (bacteria, "generic_drug").
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_bacteria_drug_param(bacteria_name: &str, drug_name: &str, param_suffix: &str) -> Option<f64> {
    // This function assumes a single parameter type (e.g., "antibiotic_reduction_per_unit")
    // If you have multiple bacteria-drug specific parameters, you'll need to expand this.
    if param_suffix == "antibiotic_reduction_per_unit" {
        // Try specific (bacteria, drug) combination first
        if let Some(&value) = BACTERIA_DRUG_PARAMETERS.get(&(bacteria_name.to_string(), drug_name.to_string())) {
            return Some(value);
        }
        // Fallback to "generic_drug" for that bacteria if no specific drug is found
        if let Some(&value) = BACTERIA_DRUG_PARAMETERS.get(&(bacteria_name.to_string(), "generic_drug".to_string())) {
            return Some(value);
        }
    }
    None
}