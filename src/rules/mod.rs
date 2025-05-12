// src/rules/mod.rs
use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::config::PARAMETERS;
use rand::Rng;
use std::collections::hash_map::Entry;

/// Applies model rules to an individual for one time step.
pub fn apply_rules(individual: &mut Individual, _time_step: usize) {
    let mut rng = rand::thread_rng();
    let no_infection_probability = 0.95;

    // Update age: increment by 1 day per time step
    individual.age += 1;

    // Update per-bacteria fields
    for &bacteria in BACTERIA_LIST.iter() {
        // Calculate the current total antibiotic activity level for this individual
        let current_antibiotic_activity_level: f64 = individual.cur_use_drug.iter().zip(individual.cur_level_drug.iter())
            .filter(|(&use_drug, _)| use_drug)
            .map(|(_, &level)| level)
            .sum();

        match bacteria {
            "strep_pneu" => {
                // Acquisition rule for strep_pneu
                if individual.age < 0 {
                    // Individuals with negative age cannot acquire strep_pneu
                    continue; // Move to the next bacteria for this individual
                }

                // Only attempt acquisition if the individual's level is 0 or not present
                if individual.level.get(bacteria).map_or(true, |&level| level <= 0.0) {
                    let mut acquisition_probability = *PARAMETERS.get("strep_pneu_acquisition_prob_baseline").unwrap_or(&0.0);

                    // Increased acquisition probability based on contact levels (multiplicative)
                    acquisition_probability *= (*PARAMETERS.get("strep_pneu_adult_contact_acq_rate_ratio_per_unit").unwrap_or(&1.0)).powf(individual.airborne_contact_level_with_adults);
                    acquisition_probability *= (*PARAMETERS.get("strep_pneu_child_contact_acq_rate_ratio_per_unit").unwrap_or(&1.0)).powf(individual.airborne_contact_level_with_children);

                    // Reduce acquisition probability if vaccinated
                    if individual.strep_pneu_vaccination_status {
                        acquisition_probability *= (1.0 - *PARAMETERS.get("strep_pneu_vaccine_efficacy").unwrap_or(&0.0));
                    }

                    if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                        let initial_level = *PARAMETERS.get("strep_pneu_initial_infection_level").unwrap_or(&0.01);
                        individual.level.insert(bacteria, initial_level);
                        individual.date_last_infected.insert(bacteria, _time_step as i32);
                        individual.infectious_syndrome.insert(bacteria, 3);
                        // We are not updating immunity here yet
                    }
                } else {
                    // Rules for when strep_pneu is already present

                    // Update infection level based on previous level, immune response, and antibiotic activity
                    if let Entry::Occupied(mut entry) = individual.level.entry(bacteria) {
                        let current_level = *entry.get();
                        let immunity_level = individual.immune_resp.get(bacteria).unwrap_or(&0.0);
                        let baseline_change = *PARAMETERS.get("strep_pneu_level_change_rate_baseline").unwrap_or(&0.0);
                        let immunity_reduction_per_unit = *PARAMETERS.get("strep_pneu_immunity_reduction_per_unit").unwrap_or(&0.0);
                        let antibiotic_reduction_per_unit = *PARAMETERS.get("strep_pneu_antibiotic_reduction_per_unit").unwrap_or(&0.0);
                        let max_level = *PARAMETERS.get("strep_pneu_max_level").unwrap_or(&100.0);

                        // Calculate the change in level
                        let change = baseline_change - (immunity_level * immunity_reduction_per_unit) - (current_antibiotic_activity_level * antibiotic_reduction_per_unit);
                        let new_level = (current_level + change).max(0.0).min(max_level);
                        *entry.get_mut() = new_level;

                        // Clear infection if level drops below a threshold
                        if *entry.get() < 0.001 {
                            individual.level.remove(bacteria);
                            individual.infectious_syndrome.insert(bacteria, 0);
                        }
                    }

                    // We are not updating immunity here yet
                }
            }
            // Keep the existing random updates for other bacteria for now
            _ => {
                // Only attempt acquisition if the individual's level is 0 or not present
                if individual.level.get(bacteria).map_or(true, |&level| level <= 0.0) {
                    if let Entry::Occupied(mut entry) = individual.date_last_infected.entry(bacteria) {
                        *entry.get_mut() += rng.gen_range(0..=1);
                    }
                    if let Entry::Occupied(mut entry) = individual.infectious_syndrome.entry(bacteria) {
                        let val_ref: &mut i32 = entry.get_mut();
                        if rng.gen::<f64>() > no_infection_probability {
                            *val_ref = rng.gen_range(1..=10);
                        } else {
                            *val_ref = 0;
                        }
                    }
                    if let Entry::Occupied(mut entry) = individual.level.entry(bacteria) {
                        let current_val = *entry.get();
                        let rng_val: f64 = rng.gen_range(0.0..=1.0);
                        *entry.get_mut() = (current_val + rng_val).max(0.0);
                    }
                    if let Entry::Occupied(mut entry) = individual.immune_resp.entry(bacteria) {
                        let current_val = *entry.get();
                        let rng_val: f64 = rng.gen_range(0.0..=1.0);
                        *entry.get_mut() = (current_val + rng_val).max(0.0);
                    }
                    if let Entry::Occupied(mut entry) = individual.sepsis.entry(bacteria) {
                        if rng.gen::<f64>() < 0.1 {
                            *entry.get_mut() = !*entry.get();
                        }
                    }
                    if let Entry::Occupied(mut entry) = individual.level_microbiome.entry(bacteria) {
                        *entry.get_mut() += rng.gen_range(0.0..=1.0);
                    }
                } else {
                    // Rules for when other bacteria are already present
                    // ... (You will add these later) ...
                }
            }
        }
    }

    // Update vaccination statuses
    if rng.gen::<f64>() < 0.1 {
        individual.haem_infl_vaccination_status = !individual.haem_infl_vaccination_status;
    }
    if rng.gen::<f64>() < 0.1 {
        individual.strep_pneu_vaccination_status = !individual.strep_pneu_vaccination_status;
    }
    if rng.gen::<f64>() < 0.1 {
        individual.salm_typhi_vaccination_status = !individual.salm_typhi_vaccination_status;
    }
    if rng.gen::<f64>() < 0.1 {
        individual.esch_coli_vaccination_status = !individual.esch_coli_vaccination_status;
    }

    // The drug use update remains the same for now.
    for i in 0..individual.cur_use_drug.len() {
        individual.cur_use_drug[i] = rng.gen_bool(0.1);
    }
    // Stop the random update of drug levels
    // for i in 0..individual.cur_level_drug.len() {
    //     individual.cur_level_drug[i] += rng.gen_range(0.0..=1.0);
    // }

    // Update other scalar fields
    individual.current_infection_related_death_risk += rng.gen_range(0.0..=1.0);
    individual.background_all_cause_mortality_rate += rng.gen_range(0.0..=1.0);
    individual.sexual_contact_level += rng.gen_range(0.0..=1.0);
    individual.airborne_contact_level_with_adults += rng.gen_range(0.0..=1.0);
    individual.airborne_contact_level_with_children += rng.gen_range(0.0..=1.0);
    individual.oral_exposure_level += rng.gen_range(0.0..=1.0);
    individual.mosquito_exposure_level += rng.gen_range(0.0..=1.0);
    if rng.gen::<f64>() < 0.1 {
        individual.under_care = !individual.under_care;
    }
    if rng.gen::<f64>() < 0.1 {
        individual.infection_hospital_acquired = !individual.infection_hospital_acquired;
    }
    individual.current_toxicity += rng.gen_range(0.0..=1.0);
    individual.mortality_risk_current_toxicity += rng.gen_range(0.0..=1.0);

    // Update resistances
    for i in 0..BACTERIA_LIST.len() {
        for j in 0..DRUG_SHORT_NAMES.len() {
            individual.resistances[i][j].microbiome_r += rng.gen_range(0.0..=1.0);
            individual.resistances[i][j].test_r += rng.gen_range(0.0..=1.0);
            individual.resistances[i][j].activity_r += rng.gen_range(0.0..=1.0);
            individual.resistances[i][j].e_r += rng.gen_range(0.0..=1.0);
            individual.resistances[i][j].c_r += rng.gen_range(0.0..=1.0);
        }
    }
}