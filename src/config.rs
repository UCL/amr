use std::collections::HashMap;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref PARAMETERS: HashMap<&'static str, f64> = {
        let mut m = HashMap::new();
        // Strep Pneumonia Parameters (focused on level)
        m.insert("strep_pneu_acquisition_prob_baseline", 0.05);         // Baseline daily acquisition probability 0.0005
        m.insert("strep_pneu_adult_contact_acq_rate_ratio_per_unit", 1.01); // Relative increase in acquisition prob per unit adult contact (>1)
        m.insert("strep_pneu_child_contact_acq_rate_ratio_per_unit", 1.02); // Relative increase in acquisition prob per unit child contact (>1)
        m.insert("strep_pneu_initial_infection_level", 0.01);          // Initial infection level upon acquisition
        m.insert("strep_pneu_vaccine_efficacy", 0.8);
        m.insert("strep_pneu_level_change_rate_baseline", 0.05);       // Baseline daily change in infection level
        m.insert("strep_pneu_immunity_effect_on_level_change", 0.01); // Effect of immunity on the change in level (higher means more reduction)
        m.insert("strep_pneu_max_level", 100.0);                     // Maximum level strep_pneu can reach
        m.insert("strep_pneu_immunity_reduction_per_unit", 0.005);   // Reduction in level change per unit of immunity
        m.insert("strep_pneu_antibiotic_reduction_per_unit", 0.00); // set to zero for testing - Reduction in level change per unit of antibiotic activity
        m // Return the initialized HashMap
    };
}