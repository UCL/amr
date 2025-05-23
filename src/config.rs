// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref PARAMETERS: HashMap<&'static str, f64> = {
        let mut m = HashMap::new();
        // Strep Pneumonia Parameters (focused on level)
        m.insert("strep_pneu_acquisition_prob_baseline", 0.05);         // Baseline daily acquisition probability 0.0005
        m.insert("strep_pneu_adult_contact_acq_rate_ratio_per_unit", 1.01); // Relative increase in acquisition prob per unit adult contact (>1)
        m.insert("strep_pneu_child_contact_acq_rate_ratio_per_unit", 1.02); // Relative increase in acquisition prob per unit child contact (>1)
        m.insert("strep_pneu_initial_infection_level", 0.01);           // Initial infection level upon acquisition
        m.insert("strep_pneu_vaccine_efficacy", 0.8);
        m.insert("strep_pneu_level_change_rate_baseline", 0.05);         // Baseline daily change in infection level
        m.insert("strep_pneu_immunity_effect_on_level_change", 0.01);      // Effect of immunity on the change in level (higher means more reduction)
        m.insert("strep_pneu_max_level", 100.0);                         // Maximum level strep_pneu can reach
        m.insert("strep_pneu_immunity_reduction_per_unit", 0.005);     // Reduction in level change per unit of immunity
        m.insert("strep_pneu_antibiotic_reduction_per_unit", 0.00);     // set to zero for testing - Reduction in level change per unit of antibiotic activity
        m.insert("strep_pneu_immunity_increase_rate_baseline", 0.01);    // Example baseline daily increase in immune response
        m.insert("strep_pneu_immunity_increase_rate_per_day", 0.001);      // Example increase in immune response per day of infection
        m.insert("strep_pneu_immunity_increase_rate_per_level", 0.0005);    // Example increase in immune response per unit of infection level
        m.insert("strep_pneu_immunity_age_modifier", 1.0);                // Example modifier based on age (can be < 1 to reduce with age)
        m.insert("strep_pneu_baseline_immunity_level", 0.1);            // Baseline level of strep_pneu immunity - at time zero for simulation (currently set to 0.1 for all bacteria in population.rs)
        m.insert("strep_pneu_immunity_decay_rate", 0.001);               // Rate at which immunity decays per time step
        m // Return the initialized HashMap
    };
}