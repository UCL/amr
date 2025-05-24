// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref PARAMETERS: HashMap<String, f64> = { // Changed key to String for consistency
        let mut m = HashMap::new();
        // Strep Pneumonia Parameters (focused on level)
        m.insert("strep_pneu_acquisition_prob_baseline".to_string(), 0.05);
        m.insert("strep_pneu_adult_contact_acq_rate_ratio_per_unit".to_string(), 1.01);
        m.insert("strep_pneu_child_contact_acq_rate_ratio_per_unit".to_string(), 1.02);
        m.insert("strep_pneu_initial_infection_level".to_string(), 0.01);
        m.insert("strep_pneu_vaccine_efficacy".to_string(), 0.8);
        m.insert("strep_pneu_level_change_rate_baseline".to_string(), 0.05);
        m.insert("strep_pneu_immunity_effect_on_level_change".to_string(), 0.01);
        m.insert("strep_pneu_max_level".to_string(), 100.0);
        m.insert("strep_pneu_immunity_reduction_per_unit".to_string(), 0.005);
        m.insert("strep_pneu_antibiotic_reduction_per_unit".to_string(), 0.00); // set to zero for testing
        m.insert("strep_pneu_immunity_increase_rate_baseline".to_string(), 0.01);
        m.insert("strep_pneu_immunity_increase_rate_per_day".to_string(), 0.001);
        m.insert("strep_pneu_immunity_increase_rate_per_level".to_string(), 0.0005);
        m.insert("strep_pneu_immunity_age_modifier".to_string(), 1.0);
        m.insert("strep_pneu_baseline_immunity_level".to_string(), 0.1);
        m.insert("strep_pneu_immunity_decay_rate".to_string(), 0.001);

        // New parameters for Strep Pneumonia acquisition source
        m.insert("strep_pneu_environmental_acquisition_proportion".to_string(), 0.1);
        m.insert("strep_pneu_hospital_acquired_proportion".to_string(), 0.05);

        // Generic Bacteria Parameters (for the '_' match arm in rules/mod.rs)
        m.insert("generic_bacteria_acquisition_prob_baseline".to_string(), 0.01);
        m.insert("generic_bacteria_initial_infection_level".to_string(), 0.01);
        m.insert("generic_environmental_acquisition_proportion".to_string(), 0.15);
        m.insert("generic_hospital_acquired_proportion".to_string(), 0.08);
        m.insert("generic_bacteria_decay_rate".to_string(), 0.02);

        // --- NEW: Drug Initiation Parameters ---
        m.insert("drug_base_initiation_rate_per_day".to_string(), 0.0001); // Very low chance for any drug
        m.insert("drug_infection_present_multiplier".to_string(), 50.0); // Stronger chance if *any* infection exists
        m.insert("drug_test_identified_multiplier".to_string(), 20.0); // Even stronger if identified by test

        // NEW: Syndrome-specific initiation multipliers (keyed by syndrome ID as String)
        // These are multipliers on top of the base rate and infection/identified multipliers.
        // Syndrome IDs are 1-10 as per population.rs comments.
        m.insert("syndrome_1_initiation_multiplier".to_string(), 2.0); // Bloodstream
        m.insert("syndrome_2_initiation_multiplier".to_string(), 3.0); // Meningitis
        m.insert("syndrome_3_initiation_multiplier".to_string(), 1.5); // Lower Respiratory
        m.insert("syndrome_4_initiation_multiplier".to_string(), 2.5); // Endocarditis
        m.insert("syndrome_5_initiation_multiplier".to_string(), 1.8); // Peritoneal/Intra-abdominal
        m.insert("syndrome_6_initiation_multiplier".to_string(), 0.5); // Diarrhoea (might not always lead to drug)
        m.insert("syndrome_7_initiation_multiplier".to_string(), 1.2); // Urinary tract infection
        m.insert("syndrome_8_initiation_multiplier".to_string(), 2.0); // Bones/Joints
        m.insert("syndrome_9_initiation_multiplier".to_string(), 1.0); // Skin/Subcutaneous
        m.insert("syndrome_10_initiation_multiplier".to_string(), 1.5); // Typhoid/Paratyphoid/NTS

        // NEW: Drug Level Parameters (for initial setup in this step)
        m.insert("drug_initial_level".to_string(), 10.0); // Initial concentration/level when a drug is started

        // Testing and Diagnosis Parameters (these were missing from your provided config, adding them back)
        m.insert("test_delay_days".to_string(), 3.0); // Days until a test can identify infection
        m.insert("test_rate_per_day".to_string(), 0.15); // Daily probability a test identifies infection after delay


        m // Return the initialized HashMap
    };
}