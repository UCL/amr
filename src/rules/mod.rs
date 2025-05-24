// src/rules/mod.rs
use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::config::PARAMETERS;
use rand::Rng;
use std::collections::hash_map::Entry;

/// Applies model rules to an individual for one time step.
pub fn apply_rules(individual: &mut Individual, time_step: usize) {
    let mut rng = rand::thread_rng();

    // Update non-infection, bacteria or antibiotic-specific variables
    individual.age += 1;
    // These random updates will likely be replaced by more complex rules later
    individual.current_infection_related_death_risk += rng.gen_range(0.0..=1.0);
    individual.background_all_cause_mortality_rate += rng.gen_range(0.0..=1.0);
    individual.sexual_contact_level += rng.gen_range(0.0..=10.0); // Adjusted range to match population.rs
    individual.airborne_contact_level_with_adults += rng.gen_range(0.0..=10.0); // Adjusted range
    individual.airborne_contact_level_with_children += rng.gen_range(0.0..=10.0); // Adjusted range
    individual.oral_exposure_level += rng.gen_range(0.0..=10.0);
    individual.mosquito_exposure_level += rng.gen_range(0.0..=10.0);
    if rng.gen::<f64>() < 0.1 { individual.under_care = !individual.under_care; }
    individual.current_toxicity += rng.gen_range(0.0..=3.0); // Adjusted range
    individual.mortality_risk_current_toxicity += rng.gen_range(0.0..=1.0);
    if rng.gen::<f64>() < 0.1 { individual.haem_infl_vaccination_status = !individual.haem_infl_vaccination_status; }
    if rng.gen::<f64>() < 0.1 { individual.strep_pneu_vaccination_status = !individual.strep_pneu_vaccination_status; }
    if rng.gen::<f64>() < 0.1 { individual.salm_typhi_vaccination_status = !individual.salm_typhi_vaccination_status; }
    if rng.gen::<f64>() < 0.1 { individual.esch_coli_vaccination_status = !individual.esch_coli_vaccination_status; }

    // --- DRUG LOGIC START (Simplified for cur_use_drug verification) ---

    let drug_initial_level = *PARAMETERS.get("drug_initial_level").unwrap_or(&10.0);
    let drug_base_initiation_rate = *PARAMETERS.get("drug_base_initiation_rate_per_day").unwrap_or(&0.0001);
    let drug_infection_present_multiplier = *PARAMETERS.get("drug_infection_present_multiplier").unwrap_or(&50.0);
    let drug_test_identified_multiplier = *PARAMETERS.get("drug_test_identified_multiplier").unwrap_or(&20.0);
    // Retrieve the drug decay rate parameter
    let drug_decay_rate = *PARAMETERS.get("drug_decay_rate_per_day").unwrap_or(&0.3);

    // Determine general infection/diagnosis status for drug administration probability
    let has_any_infection = individual.level.values().any(|&level| level > 0.0);
    let has_any_identified_infection = individual.test_identified_infection.values().any(|&identified| identified);

    // Determine the highest syndrome-specific multiplier if any syndrome is active
    let mut syndrome_administration_multiplier: f64 = 1.0;
    for (&_bacteria_name, &syndrome_id) in individual.infectious_syndrome.iter() {
        if syndrome_id != 0 {
            let param_name = format!("syndrome_{}_initiation_multiplier", syndrome_id);
            // Get parameter using .get() with the String key
            if let Some(&multiplier) = PARAMETERS.get(&param_name) {
                syndrome_administration_multiplier = syndrome_administration_multiplier.max(multiplier);
            }
        }
    }

    // Loop through each drug to decide administration for THIS day (cur_use_drug) AND apply decay
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        // Step 1: Apply decay to the current drug level
        individual.cur_level_drug[drug_idx] *= (1.0 - drug_decay_rate);
        // Ensure drug level doesn't go below zero due to floating point inaccuracies
        if individual.cur_level_drug[drug_idx] < 0.0001 { // Using a small threshold
            individual.cur_level_drug[drug_idx] = 0.0;
        }

        let mut administration_prob = drug_base_initiation_rate;

        if has_any_infection {
            administration_prob *= drug_infection_present_multiplier;
        }

        if has_any_identified_infection {
            administration_prob *= drug_test_identified_multiplier;
        }

        administration_prob *= syndrome_administration_multiplier;

        // Clamp the probability to be between 0 and 1
        administration_prob = administration_prob.clamp(0.0, 1.0);

        // Step 2: Decide whether to administer the drug today (set cur_use_drug and potentially reset cur_level_drug)
        if rng.gen_bool(administration_prob) {
            individual.cur_use_drug[drug_idx] = true;
            individual.cur_level_drug[drug_idx] = drug_initial_level; // Administer drug: reset level to initial
        } else {
            individual.cur_use_drug[drug_idx] = false; // Not administering today, so only decay applied above matters
        }
    }
    // --- DRUG LOGIC END ---


    // Update per-bacteria fields
    for &bacteria in BACTERIA_LIST.iter() {
        // The antibiotic activity level should now correctly reflect the decayed/new drug levels
        let _current_antibiotic_activity_level_general: f64 = individual.cur_use_drug.iter()
            .zip(individual.cur_level_drug.iter())
            .filter(|(&use_drug, _)| use_drug)
            .map(|(_, &level)| level)
            .sum();

        match bacteria {
            "strep_pneu" => {
                if individual.age < 0 {
                    continue;
                }

                // Check if currently infected with strep_pneu
                let is_infected = individual.level.get(bacteria).map_or(false, |&level| level > 0.0);

                if !is_infected { // Attempt acquisition if not currently infected
                    let mut acquisition_probability = *PARAMETERS.get("strep_pneu_acquisition_prob_baseline").unwrap_or(&0.0);

                    acquisition_probability *= (*PARAMETERS.get("strep_pneu_adult_contact_acq_rate_ratio_per_unit").unwrap_or(&1.0)).powf(individual.airborne_contact_level_with_adults);
                    acquisition_probability *= (*PARAMETERS.get("strep_pneu_child_contact_acq_rate_ratio_per_unit").unwrap_or(&1.0)).powf(individual.airborne_contact_level_with_children);

                    if individual.strep_pneu_vaccination_status {
                        acquisition_probability *= 1.0 - *PARAMETERS.get("strep_pneu_vaccine_efficacy").unwrap_or(&0.0);
                    }

                    if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                        let initial_level = *PARAMETERS.get("strep_pneu_initial_infection_level").unwrap_or(&0.01);
                        individual.level.insert(bacteria, initial_level);
                        individual.date_last_infected.insert(bacteria, time_step as i32);

                        individual.infectious_syndrome.insert(bacteria, 3); // Example: Respiratory syndrome

                        let env_acquisition_chance = *PARAMETERS.get("strep_pneu_environmental_acquisition_proportion").unwrap_or(&0.1);
                        individual.cur_infection_from_environment.insert(bacteria, rng.gen::<f64>() < env_acquisition_chance);

                        let hospital_acquired_chance = *PARAMETERS.get("strep_pneu_hospital_acquired_proportion").unwrap_or(&0.1);
                        individual.infection_hospital_acquired.insert(bacteria, rng.gen::<f64>() < hospital_acquired_chance);

                        individual.test_identified_infection.insert(bacteria, false); // Always false on day of acquisition
                    }

                    if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
                        let current_immunity = *immune_entry.get();
                        let baseline_immunity = *PARAMETERS.get("strep_pneu_baseline_immunity_level").unwrap_or(&0.1);
                        let decay_rate = *PARAMETERS.get("strep_pneu_immunity_decay_rate").unwrap_or(&0.001);

                        if current_immunity > baseline_immunity {
                            *immune_entry.get_mut() = (current_immunity - decay_rate).max(baseline_immunity);
                        } else if current_immunity < baseline_immunity {
                            *immune_entry.get_mut() = (current_immunity + decay_rate).min(baseline_immunity);
                        }
                    }
                } else { // strep_pneu is already present (infection progression)
                    if let Some(strep_pneu_index) = BACTERIA_LIST.iter().position(|&b| b == "strep_pneu") {
                        for (drug_index, &use_drug) in individual.cur_use_drug.iter().enumerate() {
                            if use_drug {
                                if let Some(&level) = individual.cur_level_drug.get(drug_index) {
                                    if let Some(resistance) = individual.resistances.get_mut(strep_pneu_index).and_then(|r| r.get_mut(drug_index)) {
                                        resistance.activity_r = level * (1.0 - resistance.c_r);
                                    }
                                }
                            } else {
                                if let Some(resistance) = individual.resistances.get_mut(strep_pneu_index).and_then(|r| r.get_mut(drug_index)) {
                                    resistance.activity_r = 0.0;
                                }
                            }
                        }
                    }

                    let current_antibiotic_activity_level_for_strep_pneu: f64 = if let Some(strep_pneu_index) = BACTERIA_LIST.iter().position(|&b| b == "strep_pneu") {
                        individual.resistances.get(strep_pneu_index).map_or(0.0, |drug_resistances| {
                            drug_resistances.iter().map(|r| r.activity_r).sum()
                        })
                    } else {
                        0.0
                    };

                    if let (Some(&last_infected_time), Some(current_test_status_entry)) = (
                        individual.date_last_infected.get(bacteria),
                        individual.test_identified_infection.get_mut(bacteria)
                    ) {
                        let test_delay_days = *PARAMETERS.get("test_delay_days").unwrap_or(&3.0) as i32;
                        let test_rate_per_day = *PARAMETERS.get("test_rate_per_day").unwrap_or(&0.15);

                        if !*current_test_status_entry && (time_step as i32) >= (last_infected_time + test_delay_days) {
                            if rng.gen_bool(test_rate_per_day.clamp(0.0, 1.0)) {
                                *current_test_status_entry = true;
                            }
                        }
                    }

                    if let Entry::Occupied(mut level_entry) = individual.level.entry(bacteria) {
                        let current_level = *level_entry.get();
                        let immunity_level = individual.immune_resp.get(bacteria).unwrap_or(&0.0);
                        let baseline_change = *PARAMETERS.get("strep_pneu_level_change_rate_baseline").unwrap_or(&0.0);
                        let reduction_due_to_immune_resp = *PARAMETERS.get("strep_pneu_immunity_effect_on_level_change").unwrap_or(&0.0);
                        let max_level = *PARAMETERS.get("strep_pneu_max_level").unwrap_or(&100.0);

                        let antibiotic_reduction_per_unit = *PARAMETERS.get("strep_pneu_antibiotic_reduction_per_unit").unwrap_or(&0.0);
                        let change = baseline_change - (immunity_level * reduction_due_to_immune_resp) - (current_antibiotic_activity_level_for_strep_pneu * antibiotic_reduction_per_unit);
                        let new_level = (current_level + change).max(0.0).min(max_level);
                        *level_entry.get_mut() = new_level;

                        if *level_entry.get() < 0.001 {
                            individual.level.remove(bacteria);
                            individual.infectious_syndrome.remove(bacteria);
                            individual.date_last_infected.remove(bacteria);
                            individual.immune_resp.remove(bacteria);
                            individual.sepsis.remove(bacteria);
                            individual.level_microbiome.remove(bacteria);

                            individual.infection_hospital_acquired.remove(bacteria);
                            individual.cur_infection_from_environment.remove(bacteria);
                            individual.test_identified_infection.insert(bacteria, false);
                        }
                    }

                    if let (Some(&infection_start_time), Some(&current_level)) = (
                        individual.date_last_infected.get(bacteria),
                        individual.level.get(bacteria),
                    ) {
                        let time_since_infection = (time_step as i32) - infection_start_time;
                        let age = individual.age;
                        let mut immune_increase = *PARAMETERS.get("strep_pneu_immunity_increase_rate_baseline").unwrap_or(&0.0);
                        immune_increase += time_since_infection as f64 * *PARAMETERS.get("strep_pneu_immunity_increase_rate_per_day").unwrap_or(&0.0);
                        immune_increase += current_level * *PARAMETERS.get("strep_pneu_immunity_increase_rate_per_level").unwrap_or(&0.0);
                        let age_modifier = *PARAMETERS.get("strep_pneu_immunity_age_modifier").unwrap_or(&1.0);
                        immune_increase *= age_modifier.powf((-age as f64 / 365.0) / 50.0);
                        if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
                            *immune_entry.get_mut() = (*immune_entry.get() + immune_increase).max(0.0);
                        }
                    }
                }
            }
            _ => {
                let is_infected = individual.level.get(bacteria).map_or(false, |&level| level > 0.0);

                if !is_infected {
                    let generic_acquisition_probability = *PARAMETERS.get("generic_bacteria_acquisition_prob_baseline").unwrap_or(&0.01);

                    if rng.gen_bool(generic_acquisition_probability.clamp(0.0, 1.0)) {
                        let initial_level = *PARAMETERS.get("generic_bacteria_initial_infection_level").unwrap_or(&0.01);
                        individual.level.insert(bacteria, initial_level);
                        individual.date_last_infected.insert(bacteria, time_step as i32);
                        individual.infectious_syndrome.insert(bacteria, rng.gen_range(1..=10));

                        let generic_env_acquisition_chance = *PARAMETERS.get("generic_environmental_acquisition_proportion").unwrap_or(&0.1);
                        individual.cur_infection_from_environment.insert(bacteria, rng.gen::<f64>() < generic_env_acquisition_chance);

                        let generic_hospital_acquired_chance = *PARAMETERS.get("generic_hospital_acquired_proportion").unwrap_or(&0.05);
                        individual.infection_hospital_acquired.insert(bacteria, rng.gen::<f64>() < generic_hospital_acquired_chance);

                        individual.test_identified_infection.insert(bacteria, false);
                    }
                } else {
                    if let (Some(&last_infected_time), Some(current_test_status_entry)) = (
                        individual.date_last_infected.get(bacteria),
                        individual.test_identified_infection.get_mut(bacteria)
                    ) {
                        let test_delay_days = *PARAMETERS.get("test_delay_days").unwrap_or(&3.0) as i32;
                        let test_rate_per_day = *PARAMETERS.get("test_rate_per_day").unwrap_or(&0.15);

                        if !*current_test_status_entry && (time_step as i32) >= (last_infected_time + test_delay_days) {
                            if rng.gen_bool(test_rate_per_day.clamp(0.0, 1.0)) {
                                *current_test_status_entry = true;
                            }
                        }
                    }

                    if let Entry::Occupied(mut level_entry) = individual.level.entry(bacteria) {
                        let current_level = *level_entry.get();
                        let decay_rate_generic = *PARAMETERS.get("generic_bacteria_decay_rate").unwrap_or(&0.02);
                        let new_level = (current_level - decay_rate_generic).max(0.0);
                        *level_entry.get_mut() = new_level;

                        if *level_entry.get() < 0.001 {
                            individual.level.remove(bacteria);
                            individual.infectious_syndrome.remove(bacteria);
                            individual.date_last_infected.remove(bacteria);
                            individual.immune_resp.remove(bacteria);
                            individual.sepsis.remove(bacteria);
                            individual.level_microbiome.remove(bacteria);

                            individual.infection_hospital_acquired.remove(bacteria);
                            individual.cur_infection_from_environment.remove(bacteria);
                            individual.test_identified_infection.insert(bacteria, false);
                        }
                    }
                }
            }
        }
    }

    // Update resistances - STOPPED RANDOM UPDATES
    for i in 0..BACTERIA_LIST.len() {
        for j in 0..DRUG_SHORT_NAMES.len() {
            individual.resistances[i][j].microbiome_r += 0.0;
            individual.resistances[i][j].test_r += 0.0;
            // activity_r for strep_pneu is set within its match arm above.
            // For other bacteria, activity_r will remain 0.0 unless updated there.
            individual.resistances[i][j].e_r += 0.0;
            individual.resistances[i][j].c_r += 0.0;
        }
    }

    // Print all resistances for Individual 0 (if it's individual 0)
    if individual.id == 0 {
        println!("-------------------------------------");
        println!("--- Resistance Status (Individual 0) ---");
        for (bacteria_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
            // ONLY PRINT FOR STREP_PNEU
            if bacteria_name == "strep_pneu" {
                for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                    // Ensure we have valid indices before trying to access
                    if let Some(resistance) = individual.resistances.get(bacteria_idx).and_then(|r_vec| r_vec.get(drug_idx)) {
                        println!("  {} resistance to {}:", bacteria_name, drug_name);
 //                     println!("    microbiome_r: {:.2}", resistance.microbiome_r);
 //                     println!("    test_r: {:.2}", resistance.test_r);
 //                     println!("    activity_r: {:.2}", resistance.activity_r);
 //                     println!("    e_r: {:.2}", resistance.e_r);
                        println!("    c_r: {:.2}", resistance.c_r);
                    }
                }
            }
        }
        println!("-------------------------------------");
    }

    // Print Drug Levels for Individual 0
    if individual.id == 0 {
        println!("--- Drug Levels (Individual 0) ---");
        for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
            println!("  {}: cur_use_drug = {}, cur_level_drug = {:.2}",
                     drug_name,
                     individual.cur_use_drug[drug_idx],
                     individual.cur_level_drug[drug_idx]);
        }
        println!("-------------------------------------");
    }


    // Check for death, with separate checks for each cause
    if individual.date_of_death.is_none() {
        if rng.gen::<f64>() < individual.background_all_cause_mortality_rate.clamp(0.0, 1.0) {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = Some("background".to_string());
            return;
        }

        if rng.gen::<f64>() < individual.current_infection_related_death_risk.clamp(0.0, 1.0) {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = Some("infection".to_string());
            return;
        }

        if rng.gen::<f64>() < individual.mortality_risk_current_toxicity.clamp(0.0, 1.0) {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = Some("toxicity".to_string());
            return;
        }
    }
}