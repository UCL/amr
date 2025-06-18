// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES}; // Import both lists

// --- Global Simulation Parameters ---
lazy_static! {
    pub static ref PARAMETERS: HashMap<String, f64> = {
        let mut map = HashMap::new();

        // General Drug Parameters
        map.insert("drug_initial_level".to_string(), 10.0);
        map.insert("drug_base_initiation_rate_per_day".to_string(), 0.0000); // 0.0001
        map.insert("drug_infection_present_multiplier".to_string(), 50.0);
        map.insert("drug_test_identified_multiplier".to_string(), 20.0);
        map.insert("drug_decay_rate_per_day".to_string(), 0.3);

        // General Acquisition & Resistance Parameters
        map.insert("environmental_majority_r_level_for_new_acquisition".to_string(), 0.0);
        map.insert("hospital_majority_r_level_for_new_acquisition".to_string(), 0.0);
        map.insert("max_resistance_level".to_string(), 10.0);
        map.insert("majority_r_evolution_rate_per_day_when_drug_present".to_string(), 0.001);

        // Testing Parameters
        map.insert("test_delay_days".to_string(), 3.0);
        map.insert("test_rate_per_day".to_string(), 0.15);

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
            map.insert(format!("{}_immunity_increase_rate_baseline", bacteria), 0.0);
            map.insert(format!("{}_immunity_increase_rate_per_day", bacteria), 0.01);
            map.insert(format!("{}_immunity_increase_rate_per_level", bacteria), 0.05);
            map.insert(format!("{}_immunity_age_modifier", bacteria), 1.0);
            map.insert(format!("{}_baseline_immunity_level", bacteria), 0.0);
            map.insert(format!("{}_immunity_decay_rate", bacteria), 0.1);
        }

        // --- Overrides for Specific Bacteria (Customize these as needed) ---

        // acinetobac_bau Parameters
        map.insert("acinetobac_bau_acquisition_prob_baseline".to_string(), 0.2  );
        map.insert("acinetobac_bau_hospital_acquired_proportion".to_string(), 0.15); // Often hospital-acquired

        map.insert("acinetobac_bau_immunity_increase_rate_per_day".to_string(), 2.0 );
        map.insert("acinetobac_bau_immunity_increase_rate_per_level".to_string(), 2.0);
        map.insert("acinetobac_bau_immunity_increase_effect_on_level_change".to_string(), 2.0);



        // Add more specific parameters for acinetobac_bau if needed

/*      // citrobac_spec Parameters
        map.insert("citrobac_spec_acquisition_prob_baseline".to_string(), 0.006);
        // Add more specific parameters for citrobac_spec if needed

        // enterobac_spec Parameters
        map.insert("enterobac_spec_acquisition_prob_baseline".to_string(), 0.009);
        // Add more specific parameters for enterobac_spec if needed

        // enterococ_faeca Parameters
        map.insert("enterococ_faeca_acquisition_prob_baseline".to_string(), 0.004);
        // Add more specific parameters for enterococ_faeca if needed

        // enterococ_faeci Parameters
        map.insert("enterococ_faeci_acquisition_prob_baseline".to_string(), 0.004);
        // Add more specific parameters for enterococ_faeci if needed

        // esch_coli Parameters (already present, re-listing for completeness)
        map.insert("esch_coli_acquisition_prob_baseline".to_string(), 0.008);
        map.insert("esch_coli_initial_infection_level".to_string(), 0.02);
        map.insert("esch_coli_environmental_acquisition_proportion".to_string(), 0.12);
        map.insert("esch_coli_hospital_acquired_proportion".to_string(), 0.07);
        map.insert("esch_coli_decay_rate".to_string(), 0.018);
        map.insert("esch_coli_sexual_contact_acq_rate_ratio_per_unit".to_string(), 1.8);
        map.insert("esch_coli_vaccine_efficacy".to_string(), 0.0);
        map.insert("esch_coli_baseline_immunity_level".to_string(), 0.0);
        map.insert("esch_coli_immunity_decay_rate".to_string(), 0.0);
        map.insert("esch_coli_level_change_rate_baseline".to_string(), 0.05);
        map.insert("esch_coli_immunity_effect_on_level_change".to_string(), 0.01);
        map.insert("esch_coli_max_level".to_string(), 90.0);
        map.insert("esch_coli_immunity_increase_rate_baseline".to_string(), 0.004);
        map.insert("esch_coli_immunity_increase_rate_per_day".to_string(), 0.01);
        map.insert("esch_coli_immunity_increase_rate_per_level".to_string(), 0.05);
        map.insert("esch_coli_immunity_age_modifier".to_string(), 1.0);

        // kleb_pneu Parameters
        map.insert("kleb_pneu_acquisition_prob_baseline".to_string(), 0.011);
        map.insert("kleb_pneu_hospital_acquired_proportion".to_string(), 0.18); // Often hospital-acquired
        // Add more specific parameters for kleb_pneu if needed

        // morg_spec Parameters
        map.insert("morg_spec_acquisition_prob_baseline".to_string(), 0.005);
        // Add more specific parameters for morg_spec if needed

        // prot_spec Parameters
        map.insert("prot_spec_acquisition_prob_baseline".to_string(), 0.007);
        // Add more specific parameters for prot_spec if needed

        // serrat_spec Parameters
        map.insert("serrat_spec_acquisition_prob_baseline".to_string(), 0.006);
        // Add more specific parameters for serrat_spec if needed

        // pseud_aerug Parameters
        map.insert("pseud_aerug_acquisition_prob_baseline".to_string(), 0.015);
        map.insert("pseud_aerug_hospital_acquired_proportion".to_string(), 0.2); // High hospital association
        // Add more specific parameters for pseud_aerug if needed

        // staph_aureus Parameters
        map.insert("staph_aureus_acquisition_prob_baseline".to_string(), 0.01);
        // Add more specific parameters for staph_aureus if needed

        // strep_pneu Parameters (already present, re-listing for completeness)
        map.insert("strep_pneu_acquisition_prob_baseline".to_string(), 0.005);
        map.insert("strep_pneu_adult_contact_acq_rate_ratio_per_unit".to_string(), 1.2);
        map.insert("strep_pneu_child_contact_acq_rate_ratio_per_unit".to_string(), 1.5);
        map.insert("strep_pneu_vaccine_efficacy".to_string(), 0.8);
        map.insert("strep_pneu_initial_infection_level".to_string(), 0.05);
        map.insert("strep_pneu_environmental_acquisition_proportion".to_string(), 0.05);
        map.insert("strep_pneu_hospital_acquired_proportion".to_string(), 0.1);
        map.insert("strep_pneu_baseline_immunity_level".to_string(), 0.1);
        map.insert("strep_pneu_immunity_decay_rate".to_string(), 0.001);
        map.insert("strep_pneu_level_change_rate_baseline".to_string(), 0.05);
        map.insert("strep_pneu_immunity_effect_on_level_change".to_string(), 0.01);
        map.insert("strep_pneu_max_level".to_string(), 100.0);
        map.insert("strep_pneu_immunity_increase_rate_baseline".to_string(), 0.001);
        map.insert("strep_pneu_immunity_increase_rate_per_day".to_string(), 0.01);
        map.insert("strep_pneu_immunity_increase_rate_per_level".to_string(), 0.05);
        map.insert("strep_pneu_immunity_age_modifier".to_string(), 1.0);

        // salm_typhi Parameters (already present, re-listing for completeness)
        map.insert("salm_typhi_acquisition_prob_baseline".to_string(), 0.002);
        map.insert("salm_typhi_initial_infection_level".to_string(), 0.03);
        map.insert("salm_typhi_environmental_acquisition_proportion".to_string(), 0.15);
        map.insert("salm_typhi_hospital_acquired_proportion".to_string(), 0.02);
        map.insert("salm_typhi_decay_rate".to_string(), 0.01);
        map.insert("salm_typhi_oral_exposure_acq_rate_ratio_per_unit".to_string(), 1.3);
        map.insert("salm_typhi_vaccine_efficacy".to_string(), 0.6);
        map.insert("salm_typhi_baseline_immunity_level".to_string(), 0.05);
        map.insert("salm_typhi_immunity_decay_rate".to_string(), 0.002);
        map.insert("salm_typhi_level_change_rate_baseline".to_string(), 0.05);
        map.insert("salm_typhi_immunity_effect_on_level_change".to_string(), 0.01);
        map.insert("salm_typhi_max_level".to_string(), 120.0);
        map.insert("salm_typhi_immunity_increase_rate_baseline".to_string(), 0.003);
        map.insert("salm_typhi_immunity_increase_rate_per_day".to_string(), 0.01);
        map.insert("salm_typhi_immunity_increase_rate_per_level".to_string(), 0.05);
        map.insert("salm_typhi_immunity_age_modifier".to_string(), 1.1);

        // salm_parat_a Parameters
        map.insert("salm_parat_a_acquisition_prob_baseline".to_string(), 0.0025);
        map.insert("salm_parat_a_oral_exposure_acq_rate_ratio_per_unit".to_string(), 1.2);
        // Add more specific parameters for salm_parat_a if needed

        // inv_nt_salm Parameters (Invasive Non-Typhoidal Salmonella)
        map.insert("inv_nt_salm_acquisition_prob_baseline".to_string(), 0.003);
        // Add more specific parameters for inv_nt_salm if needed

        // shig_spec Parameters
        map.insert("shig_spec_acquisition_prob_baseline".to_string(), 0.001);
        map.insert("shig_spec_oral_exposure_acq_rate_ratio_per_unit".to_string(), 1.5);
        // Add more specific parameters for shig_spec if needed

        // n_gonorrhoeae Parameters
        map.insert("n_gonorrhoeae_acquisition_prob_baseline".to_string(), 0.003);
        map.insert("n_gonorrhoeae_sexual_contact_acq_rate_ratio_per_unit".to_string(), 2.5); // Very high for sexual contact
        // Add more specific parameters for n_gonorrhoeae if needed

        // group_a_strep Parameters
        map.insert("group_a_strep_acquisition_prob_baseline".to_string(), 0.008);
        map.insert("group_a_strep_child_contact_acq_rate_ratio_per_unit".to_string(), 1.3);
        // Add more specific parameters for group_a_strep if needed

        // group_b_strep Parameters
        map.insert("group_b_strep_acquisition_prob_baseline".to_string(), 0.007);
        // Add more specific parameters for group_b_strep if needed

        // haem_infl Parameters (already present, re-listing for completeness)
        map.insert("haem_infl_acquisition_prob_baseline".to_string(), 0.012);
        map.insert("haem_infl_initial_infection_level".to_string(), 0.015);
        map.insert("haem_infl_environmental_acquisition_proportion".to_string(), 0.08);
        map.insert("haem_infl_hospital_acquired_proportion".to_string(), 0.06);
        map.insert("haem_infl_decay_rate".to_string(), 0.025);
        map.insert("haem_infl_vaccine_efficacy".to_string(), 0.7);
        map.insert("haem_infl_baseline_immunity_level".to_string(), 0.15);
        map.insert("haem_infl_immunity_decay_rate".to_string(), 0.0005);
        map.insert("haem_infl_level_change_rate_baseline".to_string(), 0.05);
        map.insert("haem_infl_immunity_effect_on_level_change".to_string(), 0.01);
        map.insert("haem_infl_max_level".to_string(), 80.0);
        map.insert("haem_infl_immunity_increase_rate_baseline".to_string(), 0.006);
        map.insert("haem_infl_immunity_increase_rate_per_day".to_string(), 0.01);
        map.insert("haem_infl_immunity_increase_rate_per_level".to_string(), 0.05);
        map.insert("haem_infl_immunity_age_modifier".to_string(), 0.9);
*/

        map
    };

    // --- Bacteria-Drug Specific Parameters ---
    pub static ref BACTERIA_DRUG_PARAMETERS: HashMap<(String, String), f64> = {
        let mut map = HashMap::new();

        // Default antibiotic reduction for ALL bacteria-drug combinations
        for &bacteria in BACTERIA_LIST.iter() {
            for &drug in DRUG_SHORT_NAMES.iter() {
                map.insert((bacteria.to_string(), drug.to_string()), 0.004); // Generic reduction
            }
        }

        // Overrides for specific bacteria-drug combinations (Customize these as needed)

        // strep_pneu drug effectiveness
        map.insert(("strep_pneu".to_string(), "penicilling".to_string()), 0.008);
        map.insert(("strep_pneu".to_string(), "amoxicillin".to_string()), 0.007);
        map.insert(("strep_pneu".to_string(), "azithromycin".to_string()), 0.003);

        // haem_infl drug effectiveness
        map.insert(("haem_infl".to_string(), "amoxicillin".to_string()), 0.006);
        map.insert(("haem_infl".to_string(), "azithromycin".to_string()), 0.005);

        // salm_typhi drug effectiveness
        map.insert(("salm_typhi".to_string(), "ciprofloxacin".to_string()), 0.009);
        map.insert(("salm_typhi".to_string(), "azithromycin".to_string()), 0.007);

        // esch_coli drug effectiveness
        map.insert(("esch_coli".to_string(), "trim_sulf".to_string()), 0.005);
        map.insert(("esch_coli".to_string(), "ciprofloxacin".to_string()), 0.006);

        // pseud_aerug drug effectiveness (example)
        map.insert(("pseud_aerug".to_string(), "meropenem".to_string()), 0.01);
        map.insert(("pseud_aerug".to_string(), "cefepime".to_string()), 0.008);

        // staph_aureus drug effectiveness (example)
        map.insert(("staph_aureus".to_string(), "vancomycin".to_string()), 0.009);
        map.insert(("staph_aureus".to_string(), "linezolid".to_string()), 0.008);

        // n_gonorrhoeae drug effectiveness (example)
        map.insert(("n_gonorrhoeae".to_string(), "ceftriaxone".to_string()), 0.0095);

        // You can add more specific drug effectiveness for each bacteria as needed
        // For example:
        // map.insert(("acinetobac_bau".to_string(), "imipenem_c".to_string()), 0.009);
        // map.insert(("kleb_pneu".to_string(), "meropenem".to_string()), 0.009);

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

/// Retrieves a bacteria-drug-specific simulation parameter.
/// It directly looks up the specific (bacteria_name, drug_name) pair.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_bacteria_drug_param(bacteria_name: &str, drug_name: &str, param_suffix: &str) -> Option<f64> {
    if param_suffix != "antibiotic_reduction_per_unit" {
        return None; // Ensure we only look for the expected parameter type
    }

    let specific_key = (bacteria_name.to_string(), drug_name.to_string());
    BACTERIA_DRUG_PARAMETERS.get(&specific_key).copied()
}