// src/rules/mod.rs
use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES};
// IMPORTANT: Updated this import to bring in the new helper functions
use crate::config::{get_global_param, get_bacteria_param, get_bacteria_drug_param};
use rand::Rng;
use std::collections::hash_map::Entry;
use rand::seq::SliceRandom;
use std::collections::HashMap;

/// Applies model rules to an individual for one time step.
pub fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    global_majority_r_proportions: &HashMap<(usize, usize), f64>,
    majority_r_positive_values_by_combo: &HashMap<(usize, usize), Vec<f64>>,
    bacteria_indices: &HashMap<&'static str, usize>,
    drug_indices: &HashMap<&'static str, usize>,
) {
    let mut rng = rand::thread_rng();

    // Update non-infection, bacteria or antibiotic-specific variables
    individual.age += 1;
    // These lines below are placeholders for complex model logic and will likely be replaced
    // with more sophisticated calculations in a full model.
    // NOTE: For these random ranges, you might want to consider putting min/max values
    // in config.rs as well, if they are meant to be adjustable parameters.
    individual.current_infection_related_death_risk = rng.gen_range(0.0..=0.001);
    individual.background_all_cause_mortality_rate = rng.gen_range(0.0..=0.00001);
    individual.sexual_contact_level = rng.gen_range(0.0..=1.0);
    individual.airborne_contact_level_with_adults = rng.gen_range(0.0..=1.0);
    individual.airborne_contact_level_with_children = rng.gen_range(0.0..=1.0);
    individual.oral_exposure_level = rng.gen_range(0.0..=1.0);
    individual.mosquito_exposure_level = rng.gen_range(0.0..=1.0);
    if rng.gen::<f64>() < 0.01 { individual.under_care = !individual.under_care; }
    individual.current_toxicity = (individual.current_toxicity + rng.gen_range(-0.5..=0.5)).max(0.0);
    individual.mortality_risk_current_toxicity = rng.gen_range(0.0..=0.0001);
    // Vaccination status changes very slowly, if at all, over time.
    if rng.gen::<f64>() < 0.0001 { individual.haem_infl_vaccination_status = !individual.haem_infl_vaccination_status; }
    if rng.gen::<f64>() < 0.0001 { individual.strep_pneu_vaccination_status = !individual.strep_pneu_vaccination_status; }
    if rng.gen::<f64>() < 0.0001 { individual.salm_typhi_vaccination_status = !individual.salm_typhi_vaccination_status; }
    if rng.gen::<f64>() < 0.0001 { individual.esch_coli_vaccination_status = !individual.esch_coli_vaccination_status; }

    // --- DRUG LOGIC START ---
    // Using get_global_param for drug-related parameters
    let drug_initial_level = get_global_param("drug_initial_level").unwrap_or(10.0);
    let drug_base_initiation_rate = get_global_param("drug_base_initiation_rate_per_day").unwrap_or(0.0001);
    let drug_infection_present_multiplier = get_global_param("drug_infection_present_multiplier").unwrap_or(50.0);
    let drug_test_identified_multiplier = get_global_param("drug_test_identified_multiplier").unwrap_or(20.0);
    let drug_decay_rate = get_global_param("drug_decay_rate_per_day").unwrap_or(0.3);

    let has_any_infection = individual.level.values().any(|&level| level > 0.0);
    let has_any_identified_infection = individual.test_identified_infection.values().any(|&identified| identified);

    let mut syndrome_administration_multiplier: f64 = 1.0;
    for (&_bacteria_name, &syndrome_id) in individual.infectious_syndrome.iter() {
        if syndrome_id != 0 {
            let param_name = format!("syndrome_{}_initiation_multiplier", syndrome_id);
            // Use get_global_param for syndrome-specific multipliers
            if let Some(multiplier) = get_global_param(&param_name) {
                syndrome_administration_multiplier = syndrome_administration_multiplier.max(multiplier);
            }
        }
    }

    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        // Drug level decay
        individual.cur_level_drug[drug_idx] *= 1.0 - drug_decay_rate;
        if individual.cur_level_drug[drug_idx] < 0.0001 {
            individual.cur_level_drug[drug_idx] = 0.0;
        }

        // Drug initiation probability
        let mut administration_prob = drug_base_initiation_rate;
        if has_any_infection { administration_prob *= drug_infection_present_multiplier; }
        if has_any_identified_infection { administration_prob *= drug_test_identified_multiplier; }
        administration_prob *= syndrome_administration_multiplier;
        administration_prob = administration_prob.clamp(0.0, 1.0);

        if rng.gen_bool(administration_prob) {
            individual.cur_use_drug[drug_idx] = true;
            individual.cur_level_drug[drug_idx] = drug_initial_level;
        } else {
            // If drug level is already 0, stop using it
            if individual.cur_level_drug[drug_idx] == 0.0 {
                individual.cur_use_drug[drug_idx] = false;
            }
        }
    }
    // --- DRUG LOGIC END ---

    // Update per-bacteria fields
    for &bacteria in BACTERIA_LIST.iter() {
        let is_infected = individual.level.get(bacteria).map_or(false, |&level| level > 0.001);
        let b_idx = *bacteria_indices.get(bacteria).unwrap(); // Get bacteria index

        if !is_infected { // Attempt acquisition if not currently infected
            // Using get_bacteria_param for bacteria-specific acquisition probabilities
            // Using "acquisition_prob_baseline" as the suffix, which get_bacteria_param should handle by checking "bacteria_name_acquisition_prob_baseline"
            // or "generic_bacteria_acquisition_prob_baseline"
            let mut acquisition_probability = get_bacteria_param(bacteria, "acquisition_prob_baseline").unwrap_or(0.01);


            // Apply bacteria-specific acquisition modifiers
            match bacteria {
                "strep_pneu" => {
                    acquisition_probability *= get_bacteria_param(bacteria, "adult_contact_acq_rate_ratio_per_unit").unwrap_or(1.0).powf(individual.airborne_contact_level_with_adults);
                    acquisition_probability *= get_bacteria_param(bacteria, "child_contact_acq_rate_ratio_per_unit").unwrap_or(1.0).powf(individual.airborne_contact_level_with_children);
                    if individual.strep_pneu_vaccination_status {
                        acquisition_probability *= 1.0 - get_bacteria_param(bacteria, "vaccine_efficacy").unwrap_or(0.0);
                    }
                },
                _ => { /* No specific modifiers for generic bacteria yet */ },
            }

            if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                // Using get_bacteria_param for initial infection level
                let initial_level = get_bacteria_param(bacteria, "initial_infection_level").unwrap_or(0.01);
                individual.level.insert(bacteria, initial_level);
                individual.date_last_infected.insert(bacteria, time_step as i32);

                let syndrome_id = match bacteria {
                    "strep_pneu" => 3, // Respiratory syndrome
                    _ => rng.gen_range(1..=10), // Random syndrome for generic
                };
                individual.infectious_syndrome.insert(bacteria, syndrome_id);

                // Using get_bacteria_param for environmental acquisition proportion
                let env_acquisition_chance = get_bacteria_param(bacteria, "environmental_acquisition_proportion").unwrap_or(0.1);
                let is_from_environment = rng.gen::<f64>() < env_acquisition_chance;
                individual.cur_infection_from_environment.insert(bacteria, is_from_environment);

                // Using get_bacteria_param for hospital acquired proportion
                let hospital_acquired_chance = get_bacteria_param(bacteria, "hospital_acquired_proportion").unwrap_or(0.05);
                let is_hospital_acquired = rng.gen::<f64>() < hospital_acquired_chance;
                individual.infection_hospital_acquired.insert(bacteria, is_hospital_acquired);

                // --- any_r AND majority_r SETTING LOGIC ON NEW INFECTION ACQUISITION ---
                let env_majority_r_level = get_global_param("environmental_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let hospital_majority_r_level = get_global_param("hospital_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(10.0);

                for drug_name_static in DRUG_SHORT_NAMES.iter() {
                    let d_idx = *drug_indices.get(drug_name_static).unwrap();
                    let resistance_data = &mut individual.resistances[b_idx][d_idx];

                    if is_from_environment {
                        // If from environment, both any_r and majority_r are set to the environmental level
                        resistance_data.majority_r = env_majority_r_level;
                        resistance_data.any_r = env_majority_r_level;
                    } else if is_hospital_acquired {
                        // If hospital-acquired, both any_r and majority_r are set to the hospital level
                        resistance_data.majority_r = hospital_majority_r_level;
                        resistance_data.any_r = hospital_majority_r_level;
                    } else {
                        // Community acquisition: any_r comes from global majority_r, but individual's majority_r starts at 0
                        if let Some(majority_r_values) = majority_r_positive_values_by_combo.get(&(b_idx, d_idx)) {
                            if let Some(&chosen_majority_r_value) = majority_r_values.choose(&mut rng) {
                                // Ensure sampled any_r is within 0 to max_resistance_level and is an integer.
                                // Round to nearest integer if sampling a float that isn't already an int.
                                resistance_data.any_r = chosen_majority_r_value.round().min(max_resistance_level).max(0.0);
                            } else {
                                resistance_data.any_r = 0.0; // Fallback
                            }
                        } else {
                            // If there's no data in majority_r_positive_values_by_combo for this combo, any_r is 0
                            resistance_data.any_r = 0.0;
                        }
                        resistance_data.majority_r = 0.0; // Individual's majority_r starts at 0 for community acquisition
                    }
                }
                // --- END GENERALIZED any_r AND majority_r SETTING LOGIC ---

                individual.test_identified_infection.insert(bacteria, false);
            }

            // Immunity decay for strep_pneu only
            if bacteria == "strep_pneu" {
                if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
                    let current_immunity = *immune_entry.get();
                    // Using get_bacteria_param for immunity levels and decay rates
                    let baseline_immunity = get_bacteria_param(bacteria, "baseline_immunity_level").unwrap_or(0.1);
                    let decay_rate = get_bacteria_param(bacteria, "immunity_decay_rate").unwrap_or(0.001);

                    if current_immunity > baseline_immunity {
                        *immune_entry.get_mut() = (current_immunity - decay_rate).max(baseline_immunity);
                    } else if current_immunity < baseline_immunity {
                        *immune_entry.get_mut() = (current_immunity + decay_rate).min(baseline_immunity);
                    }
                }
            }
        } else { // Bacteria is already present (infection progression)
            // --- majority_r EVOLUTION LOGIC ---
            let majority_r_evolution_rate = get_global_param("majority_r_evolution_rate_per_day_when_drug_present").unwrap_or(0.0);
            let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(10.0);

            if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                for (drug_index, &use_drug) in individual.cur_use_drug.iter().enumerate() {
                    let resistance_data = &mut individual.resistances[bacteria_full_idx][drug_index];

                    // Rule: If majority_r is 0, but any_r is non-zero, and drug is present, majority_r can "catch up" to any_r
                    if resistance_data.majority_r == 0.0 && resistance_data.any_r > 0.0 && use_drug {
                        if rng.gen_bool(majority_r_evolution_rate) {
                            resistance_data.majority_r = resistance_data.any_r;
                        }
                    }

                    // Rule: If majority_r is non-zero, any_r must be equal to majority_r
                    // This handles cases where majority_r just evolved, or if majority_r was non-zero from acquisition.
                    if resistance_data.majority_r > 0.0 {
                        resistance_data.any_r = resistance_data.majority_r;
                    } else {
                        // If majority_r became 0, any_r also becomes 0 if no drug pressure
                        // This prevents any_r > 0 while majority_r == 0 after evolution if the drug is removed
                        if resistance_data.any_r > 0.0 && !use_drug {
                            resistance_data.any_r = 0.0;
                        }
                    }

                    // Ensure resistance levels are always within 0 to max_resistance_level and are integers
                    resistance_data.majority_r = resistance_data.majority_r.round().min(max_resistance_level).max(0.0);
                    resistance_data.any_r = resistance_data.any_r.round().min(max_resistance_level).max(0.0);

                    // Set activity_r based on current drug level and normalized resistance (any_r)
                    if use_drug {
                        if let Some(&level) = individual.cur_level_drug.get(drug_index) {
                            // Normalize any_r for the activity calculation (e.g., 0-10 becomes 0.0-1.0)
                            let normalized_any_r = resistance_data.any_r / max_resistance_level;
                            resistance_data.activity_r = level * (1.0 - normalized_any_r);
                        }
                    } else {
                        resistance_data.activity_r = 0.0;
                    }
                }
            }


            let current_antibiotic_activity_level_for_bacteria: f64 = if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                individual.resistances.get(bacteria_full_idx).map_or(0.0, |drug_resistances| {
                    drug_resistances.iter().map(|r| r.activity_r).sum()
                })
            } else {
                0.0
            };

            // Testing and Diagnosis logic
            if let (Some(&last_infected_time), Some(current_test_status_entry)) = (
                individual.date_last_infected.get(bacteria),
                individual.test_identified_infection.get_mut(bacteria)
            ) {
                // Using get_global_param for test parameters
                let test_delay_days = get_global_param("test_delay_days").unwrap_or(3.0) as i32;
                let test_rate_per_day = get_global_param("test_rate_per_day").unwrap_or(0.15);

                if !*current_test_status_entry && (time_step as i32) >= (last_infected_time + test_delay_days) {
                    if rng.gen_bool(test_rate_per_day.clamp(0.0, 1.0)) {
                        *current_test_status_entry = true;
                    }
                }
            }

            // Bacteria level change (growth/decay)
            if let Entry::Occupied(mut level_entry) = individual.level.entry(bacteria) {
                let current_level = *level_entry.get();
                let decay_rate = match bacteria {
                    "strep_pneu" => {
                        let immunity_level = individual.immune_resp.get(bacteria).unwrap_or(&0.0);
                        // Using get_bacteria_param for strep_pneu level change parameters
                        let baseline_change = get_bacteria_param(bacteria, "level_change_rate_baseline").unwrap_or(0.0);
                        let reduction_due_to_immune_resp = get_bacteria_param(bacteria, "immunity_effect_on_level_change").unwrap_or(0.0);
                        // For antibiotic reduction, we now use get_bacteria_drug_param
                        // Assuming you want the 'activity_r' from *any* drug applied,
                        // and this parameter scales the total antibiotic activity.
                        // You might need more complex logic here if this parameter is drug-specific.
                        let antibiotic_reduction_per_unit = get_bacteria_drug_param(bacteria, "generic_drug", "antibiotic_reduction_per_unit")
                                                                .unwrap_or(0.0);

                        // For strep_pneu, 'change' can be positive (growth) or negative (decay)
                        baseline_change - (immunity_level * reduction_due_to_immune_resp) - (current_antibiotic_activity_level_for_bacteria * antibiotic_reduction_per_unit)
                    },
                    // For generic bacteria, using 'decay_rate' as the suffix
                    _ => -get_bacteria_param(bacteria, "decay_rate").unwrap_or(0.02), // Generic bacteria only decay for now
                };

                let max_level = match bacteria {
                    "strep_pneu" => get_bacteria_param(bacteria, "max_level").unwrap_or(100.0),
                    _ => 100.0, // Default max for generic
                };

                let new_level = (current_level + decay_rate).max(0.0).min(max_level);
                *level_entry.get_mut() = new_level;

                // If bacteria level drops below a threshold, clear the infection
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
                    // When infection clears, reset resistance for this bacteria to 0
                    if let Some(b_idx_clear) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                        for drug_idx_clear in 0..DRUG_SHORT_NAMES.len() {
                            individual.resistances[b_idx_clear][drug_idx_clear].any_r = 0.0;
                            individual.resistances[b_idx_clear][drug_idx_clear].majority_r = 0.0;
                        }
                    }
                }
            }

            // Immunity increase for strep_pneu only
            if bacteria == "strep_pneu" {
                if let (Some(&infection_start_time), Some(&current_level)) = (
                    individual.date_last_infected.get(bacteria),
                    individual.level.get(bacteria),
                ) {
                    let time_since_infection = (time_step as i32) - infection_start_time;
                    let age = individual.age;
                    // Using get_bacteria_param for strep_pneu immunity increase parameters
                    let mut immune_increase = get_bacteria_param(bacteria, "immunity_increase_rate_baseline").unwrap_or(0.0);
                    immune_increase += time_since_infection as f64 * get_bacteria_param(bacteria, "immunity_increase_rate_per_day").unwrap_or(0.0);
                    immune_increase += current_level * get_bacteria_param(bacteria, "immunity_increase_rate_per_level").unwrap_or(0.0);
                    let age_modifier = get_bacteria_param(bacteria, "immunity_age_modifier").unwrap_or(1.0);
                    immune_increase *= age_modifier.powf((-age as f64 / 365.0) / 50.0); // Example age modifier for immunity
                    if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
                        *immune_entry.get_mut() = (*immune_entry.get() + immune_increase).max(0.0);
                    }
                }
            }
        }
    }

    // microbiome_r, test_r are currently kept at 0 as per previous discussions
    for i in 0..BACTERIA_LIST.len() {
        for j in 0..DRUG_SHORT_NAMES.len() {
            individual.resistances[i][j].microbiome_r = 0.0;
            individual.resistances[i][j].test_r = 0.0;
        }
    }


    // Print all resistances for Individual 0 (for debugging/observation)
    if individual.id == 0 {
        println!("-------------------------------------");
        println!("--- Resistance Status (Individual 0) ---");
        let has_relevant_resistance = individual.resistances.iter().any(|b_res| {
            b_res.iter().any(|res| res.any_r > 0.0 || res.majority_r > 0.0)
        });

        if has_relevant_resistance {
            for (bacteria_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                    if let Some(resistance) = individual.resistances.get(bacteria_idx).and_then(|r_vec| r_vec.get(drug_idx)) {
                        if resistance.any_r > 0.0 || resistance.majority_r > 0.0 { // Only print if any resistance is present
                            println!("    {} resistance to {}:", bacteria_name, drug_name);
                            println!("      any_r: {:.0}", resistance.any_r); // Print as integer
                            println!("      majority_r: {:.0}", resistance.majority_r); // Print as integer
                            println!("      activity_r: {:.4}", resistance.activity_r);
                        }
                    }
                }
            }
        } else {
            println!("    No active resistances (any_r > 0 or majority_r > 0) for Individual 0.");
        }
        println!("-------------------------------------");
    }

    // Print Drug Levels for Individual 0 (for debugging/observation)
    if individual.id == 0 {
        println!("--- Drug Levels (Individual 0) ---");
        let has_active_drugs = individual.cur_use_drug.iter().any(|&in_use| in_use) || individual.cur_level_drug.iter().any(|&level| level > 0.0);
        if has_active_drugs {
            for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                if individual.cur_use_drug[drug_idx] || individual.cur_level_drug[drug_idx] > 0.0 {
                    println!("    {}: cur_use_drug = {}, cur_level_drug = {:.2}",
                             drug_name,
                             individual.cur_use_drug[drug_idx],
                             individual.cur_level_drug[drug_idx]);
                }
            }
        } else {
            println!("    No active drug use for Individual 0.");
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