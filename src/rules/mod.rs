// src/rules/mod.rs
use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::config::PARAMETERS;
use rand::Rng;
use std::collections::hash_map::Entry;

/// Applies model rules to an individual for one time step.
pub fn apply_rules(individual: &mut Individual, _time_step: usize) {
    let mut rng = rand::thread_rng();
    let no_infection_probability = 0.95;

    // Update non-infection, bacteria or antibiotic-specific variables
    individual.age += 1;
    individual.current_infection_related_death_risk += rng.gen_range(0.0..=1.0);
    individual.background_all_cause_mortality_rate += rng.gen_range(0.0..=1.0);
    individual.sexual_contact_level += rng.gen_range(0.0..=1.0);
    individual.airborne_contact_level_with_adults += rng.gen_range(0.0..=1.0);
    individual.airborne_contact_level_with_children += rng.gen_range(0.0..=1.0);
    individual.oral_exposure_level += rng.gen_range(0.0..=1.0);
    individual.mosquito_exposure_level += rng.gen_range(0.0..=1.0);
    if rng.gen::<f64>() < 0.1 { individual.under_care = !individual.under_care; }
    individual.current_toxicity += rng.gen_range(0.0..=1.0);
    individual.mortality_risk_current_toxicity += rng.gen_range(0.0..=1.0);
    if rng.gen::<f64>() < 0.1 { individual.haem_infl_vaccination_status = !individual.haem_infl_vaccination_status; }
    if rng.gen::<f64>() < 0.1 { individual.strep_pneu_vaccination_status = !individual.strep_pneu_vaccination_status; }
    if rng.gen::<f64>() < 0.1 { individual.salm_typhi_vaccination_status = !individual.salm_typhi_vaccination_status; }
    if rng.gen::<f64>() < 0.1 { individual.esch_coli_vaccination_status = !individual.esch_coli_vaccination_status; }

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
                        acquisition_probability *= 1.0 - *PARAMETERS.get("strep_pneu_vaccine_efficacy").unwrap_or(&0.0); // Fixed parentheses
                    }

                    if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                        let initial_level = *PARAMETERS.get("strep_pneu_initial_infection_level").unwrap_or(&0.01);
                        individual.level.insert(bacteria, initial_level);
                        individual.date_last_infected.insert(bacteria, _time_step as i32);
                        individual.infectious_syndrome.insert(bacteria, 3);

                        if rng.gen::<f64>() < 0.1 {
                            individual.infection_hospital_acquired = true;
                        }
                    }

                    // Gradual decay of immune response towards baseline if above
                    if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
                        let current_immunity = *immune_entry.get();
                        let baseline_immunity = *PARAMETERS.get("strep_pneu_baseline_immunity_level").unwrap_or(&0.1); // New parameter
                        let decay_rate = *PARAMETERS.get("strep_pneu_immunity_decay_rate").unwrap_or(&0.001); // New parameter

                        if current_immunity > baseline_immunity {
                            *immune_entry.get_mut() = (current_immunity - decay_rate).max(baseline_immunity);
                        } else if current_immunity < baseline_immunity {
                            *immune_entry.get_mut() = (current_immunity + decay_rate).min(baseline_immunity);
                        }
                    }
                } else {
                    // Rules for when strep_pneu is already present

                    // Update antibiotic activity_r based on current drug use and resistance
                    if let Some(strep_pneu_index) = BACTERIA_LIST.iter().position(|&b| b == "strep_pneu") {
                        for (drug_index, &use_drug) in individual.cur_use_drug.iter().enumerate() {
                            if use_drug {
                                if let Some(&level) = individual.cur_level_drug.get(drug_index) {
                                    if let Some(resistance) = individual.resistances.get_mut(strep_pneu_index).and_then(|r| r.get_mut(drug_index)) {
                                        resistance.activity_r = level * (1.0 - resistance.c_r);
                                    }
                                }
                            } else {
                                // If drug is not in use, maybe set activity_r to 0?
                                if let Some(resistance) = individual.resistances.get_mut(strep_pneu_index).and_then(|r| r.get_mut(drug_index)) {
                                    resistance.activity_r = 0.0;
                                }
                            }
                        }
                    }

                    // Now, the current total antibiotic activity level for strep_pneu
                    let current_antibiotic_activity_level: f64 = if let Some(strep_pneu_index) = BACTERIA_LIST.iter().position(|&b| b == "strep_pneu") {
                        individual.resistances.get(strep_pneu_index).map_or(0.0, |drug_resistances| {
                            drug_resistances.iter().map(|r| r.activity_r).sum()
                        })
                    } else {
                        0.0
                    };

                    // Update infection level based on previous level, immune response, and antibiotic activity
                    if let Entry::Occupied(mut level_entry) = individual.level.entry(bacteria) {
                        let current_level = *level_entry.get();
                        let immunity_level = individual.immune_resp.get(bacteria).unwrap_or(&0.0);
                        let baseline_change = *PARAMETERS.get("strep_pneu_level_change_rate_baseline").unwrap_or(&0.0);
                        let reduction_due_to_immune_resp = *PARAMETERS.get("strep_pneu_reduction_due_to_immune_resp").unwrap_or(&0.0);
                        let max_level = *PARAMETERS.get("strep_pneu_max_level").unwrap_or(&100.0);

                        // Calculate the change in level
                        let antibiotic_reduction_per_unit = *PARAMETERS.get("strep_pneu_antibiotic_reduction_per_unit").unwrap_or(&0.0);
                        let change = baseline_change - (immunity_level * reduction_due_to_immune_resp) - (current_antibiotic_activity_level * antibiotic_reduction_per_unit);
                        let new_level = (current_level + change).max(0.0).min(max_level);
                        *level_entry.get_mut() = new_level;

                        // Clear infection if level drops below a threshold
                        if *level_entry.get() < 0.001 {
                            individual.level.remove(bacteria);
                            individual.infectious_syndrome.insert(bacteria, 0);
                            individual.infection_hospital_acquired = false;
                        }
                    }

                    // Update immune response for active strep_pneu infection
                    if let (Some(&infection_start_time), Some(&current_level)) = (
                        individual.date_last_infected.get(bacteria),
                        individual.level.get(bacteria),
                    ) {
                        let time_since_infection = (_time_step as i32) - infection_start_time;
                        let age = individual.age;
                        // Calculate the increase in immune response based on time, level, and age
                        let mut immune_increase = *PARAMETERS.get("strep_pneu_immunity_increase_rate_baseline").unwrap_or(&0.0);
                        // Influence of time since infection
                        immune_increase += time_since_infection as f64 * *PARAMETERS.get("strep_pneu_immunity_increase_rate_per_day").unwrap_or(&0.0);
                        // Influence of infection level
                        immune_increase += current_level * *PARAMETERS.get("strep_pneu_immunity_increase_rate_per_level").unwrap_or(&0.0);
                        // Influence of age (example: slower increase with older age)
                        let age_modifier = *PARAMETERS.get("strep_pneu_immunity_age_modifier").unwrap_or(&1.0);
                        immune_increase *= age_modifier.powf((-age as f64 / 365.0) / 50.0);
                        // Update the immune_resp, ensuring it doesn't go below 0.0
                        if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
                            *immune_entry.get_mut() = (*immune_entry.get() + immune_increase).max(0.0);
                        }
                    }
                }
            }
            _ => {
                // Keep the existing random updates for other bacteria for now
                // Only attempt acquisition if the individual's level is 0 or not present
                if individual.level.get(bacteria).map_or(true, |&level| level <= 0.0) {
                    if let Entry::Occupied(mut entry) = individual.date_last_infected.entry(bacteria) {
                        *entry.get_mut() += 0; // Changed from rng.gen_range and 0.0
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
                        let rng_val: f64 = 0.0;
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

    // The drug use update remains the same for now.
    for i in 0..individual.cur_use_drug.len() {
        individual.cur_use_drug[i] = false;
    }
    // Stop the random update of drug levels
    for i in 0..individual.cur_level_drug.len() {
        individual.cur_level_drug[i] += 0.0; // Changed from rng.gen_range
    }

    // Update resistances - STOPPED RANDOM UPDATES
    for i in 0..BACTERIA_LIST.len() {
        for j in 0..DRUG_SHORT_NAMES.len() {
            individual.resistances[i][j].microbiome_r += 0.0;
            individual.resistances[i][j].test_r += 0.0;
            individual.resistances[i][j].activity_r += 0.0;
            individual.resistances[i][j].e_r += 0.0;
            individual.resistances[i][j].c_r += 0.0;
        }
    }

    // Check for death, with separate checks for each cause
    if individual.date_of_death.is_none() {
        // 1. Check for death due to background mortality
        if rng.gen::<f64>() < individual.background_all_cause_mortality_rate.clamp(0.0, 1.0) {
            individual.date_of_death = Some(_time_step);
            individual.cause_of_death = Some("background".to_string());
            return; // Important: Stop processing other causes if death occurs
        }

        // 2. Check for death due to infection-related mortality
        if rng.gen::<f64>() < individual.current_infection_related_death_risk.clamp(0.0, 1.0) {
            individual.date_of_death = Some(_time_step);
            individual.cause_of_death = Some("infection".to_string()); // More specific later
            return; // Stop processing other causes if death occurs
        }

        // 3. Check for death due to toxicity
        if rng.gen::<f64>() < individual.mortality_risk_current_toxicity.clamp(0.0, 1.0) {
            individual.date_of_death = Some(_time_step);
            individual.cause_of_death = Some("toxicity".to_string());
            return; // Stop processing other causes if death occurs
        }
    }
}