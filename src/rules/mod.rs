// src/rules/mod.rs
use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::config::{get_global_param, get_bacteria_param, get_bacteria_drug_param};
use rand::Rng;
use std::collections::hash_map::Entry;
use rand::seq::SliceRandom; // REMOVED - not directly used here for sampling
use std::collections::HashMap;

/// Applies model rules to an individual for one time step.
pub fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    _global_majority_r_proportions: &HashMap<(usize, usize), f64>, // Prefixed with _
    majority_r_positive_values_by_combo: &HashMap<(usize, usize), Vec<f64>>,
    bacteria_indices: &HashMap<&'static str, usize>,
    drug_indices: &HashMap<&'static str, usize>,
) {
    let mut rng = rand::thread_rng();

    // Update non-infection, bacteria or antibiotic-specific variables
    individual.age += 1;
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

    // MODIFIED: Loop through all bacteria to update vaccination status dynamically
    // This is already using the HashMap as implemented in population.rs
    for &bacteria in BACTERIA_LIST.iter() {
        if let Entry::Occupied(mut status_entry) = individual.vaccination_status.entry(bacteria) {
            if rng.gen::<f64>() < 0.0001 { // Small chance to change status
                *status_entry.get_mut() = !*status_entry.get();
            }
        }
    }

    // --- DRUG LOGIC START ---
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
            // --- BACTERIA-SPECIFIC ACQUISITION PROBABILITY CALCULATION ---
            // These parameters are now explicitly defined for each bacteria in config.rs
            let mut acquisition_probability = get_bacteria_param(bacteria, "acquisition_prob_baseline").unwrap_or(0.01);

            // Apply contact level modifiers dynamically
            let sexual_contact_multiplier = get_bacteria_param(bacteria, "sexual_contact_acq_rate_ratio_per_unit").unwrap_or(1.0);
            let airborne_adult_contact_multiplier = get_bacteria_param(bacteria, "adult_contact_acq_rate_ratio_per_unit").unwrap_or(1.0);
            let airborne_child_contact_multiplier = get_bacteria_param(bacteria, "child_contact_acq_rate_ratio_per_unit").unwrap_or(1.0);
            let oral_exposure_multiplier = get_bacteria_param(bacteria, "oral_exposure_acq_rate_ratio_per_unit").unwrap_or(1.0);
            let mosquito_exposure_multiplier = get_bacteria_param(bacteria, "mosquito_exposure_acq_rate_ratio_per_unit").unwrap_or(1.0);

            acquisition_probability *= sexual_contact_multiplier.powf(individual.sexual_contact_level);
            acquisition_probability *= airborne_adult_contact_multiplier.powf(individual.airborne_contact_level_with_adults);
            acquisition_probability *= airborne_child_contact_multiplier.powf(individual.airborne_contact_level_with_children);
            acquisition_probability *= oral_exposure_multiplier.powf(individual.oral_exposure_level);
            acquisition_probability *= mosquito_exposure_multiplier.powf(individual.mosquito_exposure_level);

            // Apply vaccination status effect dynamically
            if individual.vaccination_status.get(bacteria).copied().unwrap_or(false) {
                let vaccine_efficacy = get_bacteria_param(bacteria, "vaccine_efficacy").unwrap_or(0.0);
                acquisition_probability *= 1.0 - vaccine_efficacy;
            }
            // --- END BACTERIA-SPECIFIC ACQUISITION PROBABILITY CALCULATION ---


            if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                let initial_level = get_bacteria_param(bacteria, "initial_infection_level").unwrap_or(0.01);
                individual.level.insert(bacteria, initial_level);
                individual.date_last_infected.insert(bacteria, time_step as i32);

                // Determine syndrome ID
                // ADDED: Explicit syndrome IDs for additional bacteria for more realistic modeling
                let syndrome_id = match bacteria {
                    "strep_pneu" => 3, // Respiratory syndrome
                    "haem_infl" => 3, // Can cause respiratory issues
                    "kleb_pneu" => 3, // Can cause pneumonia
                    "salm_typhi" => 7, // Gastrointestinal syndrome
                    "salm_parat_a" => 7,
                    "inv_nt_salm" => 7,
                    "shig_spec" => 7,
                    "esch_coli" => 7, // E.coli can also cause GI issues
                    "n_gonorrhoeae" => 8, // Genital syndrome (example ID)
                    "group_a_strep" => 9, // Skin/throat related (example ID)
                    "group_b_strep" => 10, // Neonatal/sepsis related (example ID)
                    _ => rng.gen_range(1..=10), // Random syndrome for others, adjust as needed
                };
                individual.infectious_syndrome.insert(bacteria, syndrome_id);

                let env_acquisition_chance = get_bacteria_param(bacteria, "environmental_acquisition_proportion").unwrap_or(0.1);
                let is_from_environment = rng.gen::<f64>() < env_acquisition_chance;
                individual.cur_infection_from_environment.insert(bacteria, is_from_environment);

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
                        resistance_data.majority_r = env_majority_r_level;
                        resistance_data.any_r = env_majority_r_level;
                    } else if is_hospital_acquired {
                        resistance_data.majority_r = hospital_majority_r_level;
                        resistance_data.any_r = hospital_majority_r_level;
                    } else {
                        // When acquired from other individuals, sample from global majority_r for any_r
                        if let Some(majority_r_values) = majority_r_positive_values_by_combo.get(&(b_idx, d_idx)) {
                            if let Some(&chosen_majority_r_value) = majority_r_values.choose(&mut rng) {
                                resistance_data.any_r = chosen_majority_r_value.round().min(max_resistance_level).max(0.0);
                            } else {
                                resistance_data.any_r = 0.0;
                            }
                        } else {
                            resistance_data.any_r = 0.0;
                        }
                        resistance_data.majority_r = 0.0; // majority_r starts at 0 for new non-environmental/hospital acquisitions
                    }
                }
                // --- END GENERALIZED any_r AND majority_r SETTING LOGIC ---

                individual.test_identified_infection.insert(bacteria, false);
            }

            // Immunity decay
            // Dynamically apply immunity decay for all bacteria using specific parameters
            if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
                let current_immunity = *immune_entry.get();
                let baseline_immunity = get_bacteria_param(bacteria, "baseline_immunity_level").unwrap_or(0.0);
                let decay_rate = get_bacteria_param(bacteria, "immunity_decay_rate").unwrap_or(0.0);

                if current_immunity > baseline_immunity {
                    *immune_entry.get_mut() = (current_immunity - decay_rate).max(baseline_immunity);
                } else if current_immunity < baseline_immunity {
                    *immune_entry.get_mut() = (current_immunity + decay_rate).min(baseline_immunity);
                }
            }
        } else { // Bacteria is already present (infection progression)
            // --- majority_r EVOLUTION LOGIC ---
            let majority_r_evolution_rate = get_global_param("majority_r_evolution_rate_per_day_when_drug_present").unwrap_or(0.0);
            let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(10.0);

            if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                for (drug_index, &use_drug) in individual.cur_use_drug.iter().enumerate() {
                    let resistance_data = &mut individual.resistances[bacteria_full_idx][drug_index];

                    if resistance_data.majority_r == 0.0 && resistance_data.any_r > 0.0 && use_drug {
                        // If any_r exists and drug is used, majority_r can develop
                        if rng.gen_bool(majority_r_evolution_rate) {
                            resistance_data.majority_r = resistance_data.any_r;
                        }
                    }

                    // If majority_r exists, it sets any_r. If majority_r becomes 0, any_r also reverts to 0.
                    if resistance_data.majority_r > 0.0 {
                        resistance_data.any_r = resistance_data.majority_r;
                    } else {
                        // If no majority_r, and no drug used, any_r should also revert over time if not already 0
                        if resistance_data.any_r > 0.0 && !use_drug {
                            // This logic needs to be carefully considered for how fast any_r dissipates without majority_r
                            // For simplicity, here it reverts immediately. A decay rate might be more appropriate.
                            resistance_data.any_r = 0.0;
                        }
                    }

                    resistance_data.majority_r = resistance_data.majority_r.round().min(max_resistance_level).max(0.0);
                    resistance_data.any_r = resistance_data.any_r.round().min(max_resistance_level).max(0.0);

                    // Calculate activity_r
                    if use_drug {
                        if let Some(&level) = individual.cur_level_drug.get(drug_index) {
                            let normalized_any_r = resistance_data.any_r / max_resistance_level;
                            resistance_data.activity_r = level * (1.0 - normalized_any_r);
                        }
                    } else {
                        resistance_data.activity_r = 0.0;
                    }
                }
            }

            // MODIFIED: Prefixed with _ to silence warning
            let _current_antibiotic_activity_level_for_bacteria: f64 = if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
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
                // Dynamically calculate level change for all bacteria using specific parameters
                let immunity_level = individual.immune_resp.get(bacteria).unwrap_or(&0.0);
                let baseline_change = get_bacteria_param(bacteria, "level_change_rate_baseline").unwrap_or(0.0);
                let reduction_due_to_immune_resp = get_bacteria_param(bacteria, "immunity_effect_on_level_change").unwrap_or(0.0);

                // Use get_bacteria_drug_param for antibiotic reduction, iterating through all drugs
                let mut total_antibiotic_reduction = 0.0;
                for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                    if individual.cur_use_drug[drug_idx] {
                        let antibiotic_reduction_per_unit = get_bacteria_drug_param(bacteria, drug_name, "antibiotic_reduction_per_unit")
                                                                .unwrap_or(0.0);
                        total_antibiotic_reduction += individual.cur_level_drug[drug_idx] * antibiotic_reduction_per_unit;
                    }
                }

                let decay_rate = baseline_change - (immunity_level * reduction_due_to_immune_resp) - total_antibiotic_reduction;

                let max_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);

                let new_level = (current_level + decay_rate).max(0.0).min(max_level);
                *level_entry.get_mut() = new_level;

                // If bacteria level drops below a threshold, clear the infection
                if *level_entry.get() < 0.0001 {
 //                 individual.level.remove(bacteria);
 //                 individual.infectious_syndrome.remove(bacteria);
                    individual.date_last_infected.remove(bacteria);
                    individual.immune_resp.remove(bacteria);
 //                 individual.sepsis.remove(bacteria);
 //                 individual.level_microbiome.remove(bacteria);
 //                 individual.infection_hospital_acquired.remove(bacteria);
 //                 individual.cur_infection_from_environment.remove(bacteria);
                    individual.test_identified_infection.insert(bacteria, false);
                    if let Some(b_idx_clear) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                        for drug_idx_clear in 0..DRUG_SHORT_NAMES.len() {
                            individual.resistances[b_idx_clear][drug_idx_clear].any_r = 0.0;
                            individual.resistances[b_idx_clear][drug_idx_clear].majority_r = 0.0;
                        }
                    }
                }
            }

            // Immunity increase
            // Dynamically apply immunity increase for all bacteria using specific parameters
            if let (Some(&infection_start_time), Some(&current_level)) = (
                individual.date_last_infected.get(bacteria),
                individual.level.get(bacteria),
            ) {
                let time_since_infection = (time_step as i32) - infection_start_time;
                let age = individual.age;
                let mut immune_increase = get_bacteria_param(bacteria, "immunity_increase_rate_baseline").unwrap_or(0.0);
                immune_increase += time_since_infection as f64 * get_bacteria_param(bacteria, "immunity_increase_rate_per_day").unwrap_or(0.0);
                immune_increase += current_level * get_bacteria_param(bacteria, "immunity_increase_rate_per_level").unwrap_or(0.0);
                let age_modifier = get_bacteria_param(bacteria, "immunity_age_modifier").unwrap_or(1.0);
                immune_increase *= age_modifier.powf((-age as f64 / 365.0) / 50.0); // Example age-related decay
                if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
                    *immune_entry.get_mut() = (*immune_entry.get() + immune_increase).max(0.0);
                }
            }
        }
    }

    // microbiome_r, test_r are currently kept at 0
    for i in 0..BACTERIA_LIST.len() {
        for j in 0..DRUG_SHORT_NAMES.len() {
            individual.resistances[i][j].microbiome_r = 0.0;
            individual.resistances[i][j].test_r = 0.0;
        }
    }

    // Print all resistances for Individual 0 (for debugging/observation)
    
/*  if individual.id == 0 {
        println!("-------------------------------------");
        println!("--- Resistance Status (Individual 0) ---");
        let has_relevant_resistance = individual.resistances.iter().any(|b_res| {
            b_res.iter().any(|res| res.any_r > 0.0 || res.majority_r > 0.0)
        });

        if has_relevant_resistance {
            for (bacteria_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                    if let Some(resistance) = individual.resistances.get(bacteria_idx).and_then(|r_vec| r_vec.get(drug_idx)) {
                        if resistance.any_r > 0.0 || resistance.majority_r > 0.0 {
                            println!("    {} resistance to {}:", bacteria_name, drug_name);
                            println!("      any_r: {:.0}", resistance.any_r);
                            println!("      majority_r: {:.0}", resistance.majority_r);
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
*/

/*  // Print Drug Levels for Individual 0 (for debugging/observation)
    if individual.id == 0 {
        println!("--- Drug Levels (Individual 0) ---");
        let has_active_drugs = individual.cur_use_drug.iter().any(|&in_use| in_use) || individual.cur_level_drug.iter().any(|&level| level > 0.0);
        if has_active_drugs {
            for (drug_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() { // MODIFIED: Prefixed drug_name with _ to silence warning
                println!("    {}: cur_use_drug = {}, cur_level_drug = {:.2}",
                         DRUG_SHORT_NAMES[drug_idx], // MODIFIED: Use DRUG_SHORT_NAMES directly
                         individual.cur_use_drug[drug_idx],
                         individual.cur_level_drug[drug_idx]);
            }
        } else {
            println!("    No active drug use for Individual 0.");
        }
        println!("-------------------------------------");
    }
*/

    // Check for death, with separate checks for each cause
    if individual.date_of_death.is_none() {
        if rng.gen::<f64>() < individual.background_all_cause_mortality_rate.clamp(0.0, 1.0) {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = Some("background".to_string());
        } else if rng.gen::<f64>() < individual.current_infection_related_death_risk.clamp(0.0, 1.0) {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = Some("infection".to_string());
        } else if rng.gen::<f64>() < individual.mortality_risk_current_toxicity.clamp(0.0, 1.0) {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = Some("toxicity".to_string());
        }
    }
}