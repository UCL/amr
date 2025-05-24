// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref PARAMETERS: HashMap<&'static str, f64> = {
        let mut m = HashMap::new();
        // Strep Pneumonia Parameters (focused on level)
        m.insert("strep_pneu_acquisition_prob_baseline", 0.05);
        m.insert("strep_pneu_adult_contact_acq_rate_ratio_per_unit", 1.01);
        m.insert("strep_pneu_child_contact_acq_rate_ratio_per_unit", 1.02);
        m.insert("strep_pneu_initial_infection_level", 0.01);
        m.insert("strep_pneu_vaccine_efficacy", 0.8);
        m.insert("strep_pneu_level_change_rate_baseline", 0.05);
        m.insert("strep_pneu_immunity_effect_on_level_change", 0.01);
        m.insert("strep_pneu_max_level", 100.0);
        m.insert("strep_pneu_immunity_reduction_per_unit", 0.005);
        m.insert("strep_pneu_antibiotic_reduction_per_unit", 0.00); // set to zero for testing
        m.insert("strep_pneu_immunity_increase_rate_baseline", 0.01);
        m.insert("strep_pneu_immunity_increase_rate_per_day", 0.001);
        m.insert("strep_pneu_immunity_increase_rate_per_level", 0.0005);
        m.insert("strep_pneu_immunity_age_modifier", 1.0);
        m.insert("strep_pneu_baseline_immunity_level", 0.1);
        m.insert("strep_pneu_immunity_decay_rate", 0.001);

        // New parameters for Strep Pneumonia acquisition source
        m.insert("strep_pneu_environmental_acquisition_proportion", 0.1); // Proportion of strep_pneu acquisitions that are environmental
        m.insert("strep_pneu_hospital_acquired_proportion", 0.05);      // Proportion of strep_pneu acquisitions that are hospital-acquired

        // Generic Bacteria Parameters (for the '_' match arm in rules/mod.rs)
        m.insert("generic_bacteria_acquisition_prob_baseline", 0.01);  // Baseline daily acquisition probability for other bacteria
        m.insert("generic_bacteria_initial_infection_level", 0.01);     // Initial infection level for other bacteria
        m.insert("generic_environmental_acquisition_proportion", 0.15); // Proportion of generic bacteria acquisitions that are environmental
        m.insert("generic_hospital_acquired_proportion", 0.08);        // Proportion of generic bacteria acquisitions that are hospital-acquired
        m.insert("generic_bacteria_decay_rate", 0.02);                  // Rate at which generic bacterial infection levels decay

        m // Return the initialized HashMap
    };
}