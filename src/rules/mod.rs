// src/rules/mod.rs


use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::config::{get_global_param, get_bacteria_param, get_bacteria_drug_param, get_drug_param};
use rand::Rng;
use std::collections::hash_map::Entry;
use rand::seq::SliceRandom;
use std::collections::HashMap;

/// Applies model rules to an individual for one time step.
pub fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    _global_majority_r_proportions: &HashMap<(usize, usize), f64>,
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

    // Loop through all bacteria to update vaccination status dynamically
    for &bacteria in BACTERIA_LIST.iter() {
        if let Entry::Occupied(mut status_entry) = individual.vaccination_status.entry(bacteria) {
            if rng.gen::<f64>() < 0.0001 { // Small chance to change status
                *status_entry.get_mut() = !*status_entry.get();
            }
        }
    }

    // --- DRUG LOGIC START ---
    let drug_base_initiation_rate = get_global_param("drug_base_initiation_rate_per_day").unwrap_or(0.0001);
    let drug_infection_present_multiplier = get_global_param("drug_infection_present_multiplier").unwrap_or(50.0);
    let already_on_drug_initiation_multiplier = get_global_param("already_on_drug_initiation_multiplier").unwrap_or(0.0001);
    let drug_test_identified_multiplier = get_global_param("drug_test_identified_multiplier").unwrap_or(20.0);
    let drug_decay_rate = get_global_param("drug_decay_rate_per_day").unwrap_or(0.3);
    let double_dose_probability = get_global_param("double_dose_probability_if_identified_infection").unwrap_or(0.1);
    let random_drug_cessation_prob = get_global_param("random_drug_cessation_probability").unwrap_or(0.001);


    let has_any_infection = individual.level.values().any(|&level| level > 0.0);
    let initial_on_any_antibiotic = individual.cur_use_drug.iter().any(|&identified| identified);
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

    let mut drugs_initiated_this_time_step: usize = 0;

    // --- NEW DRUG STOPPING LOGIC (Loop 1) ---
    // Iterate over drugs to decide if they should be stopped.
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];

        if individual.cur_use_drug[drug_idx] { // Only consider stopping drugs that are currently in use
            let mut relevant_infection_active_for_this_drug = false;
            for &bacteria_name in BACTERIA_LIST.iter() {
                // Check if the individual is infected with this bacteria
                if let Some(&bacteria_level) = individual.level.get(bacteria_name) {
                    if bacteria_level > 0.0001 { // Check for active infection
                        // Check if this drug is effective against this bacteria
                        let drug_reduction_efficacy = get_bacteria_drug_param(
                            bacteria_name,
                            drug_name,
                            "bacteria_level_reduction_per_unit_of_drug",
                        ).unwrap_or(0.0);

                        if drug_reduction_efficacy > 0.0 {
                            relevant_infection_active_for_this_drug = true;
                            break; // Found an active, relevant infection for this drug
                        }
                    }
                }
            }

            // Condition to stop the drug:
            // 1. No relevant active infection (primary reason)
            // 2. A small random chance to stop early (e.g., non-adherence, side effects)
            if !relevant_infection_active_for_this_drug || rng.gen_bool(random_drug_cessation_prob) {
                individual.cur_use_drug[drug_idx] = false;
                // println!("DEBUG: Individual {} stopped drug {} at day {}. No relevant infection: {}. Random stop: {}",
                //      individual.id, drug_name, time_step, !relevant_infection_active_for_this_drug, rng.gen_bool(random_drug_cessation_prob));
            }
        }
    }

    // --- DRUG LEVEL UPDATE (Decay or Replenish) (Loop 2) ---
    // Apply decay if stopped, or replenish if still in use.
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];
        let drug_initial_level = get_drug_param(drug_name, "initial_level").unwrap_or(10.0);

        if individual.cur_use_drug[drug_idx] {
            // If actively using, set to initial level (representing daily dose)
            individual.cur_level_drug[drug_idx] = drug_initial_level;
        } else {
            // If not using, decay the level
            individual.cur_level_drug[drug_idx] = (individual.cur_level_drug[drug_idx] - drug_decay_rate).max(0.0);
        }
    }


    // --- NEW DRUG INITIATION LOGIC (Loop 3) ---
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx]; // Get the drug name for parameter lookup

        let mut administration_prob = drug_base_initiation_rate;
        if has_any_infection { administration_prob *= drug_infection_present_multiplier; }
        if has_any_identified_infection { administration_prob *= drug_test_identified_multiplier; }
        
        if initial_on_any_antibiotic || drugs_initiated_this_time_step > 0 {
            administration_prob *= already_on_drug_initiation_multiplier;
        }
        
        administration_prob *= syndrome_administration_multiplier;
        administration_prob = administration_prob.clamp(0.0, 1.0);

        // Only attempt to initiate if the drug is NOT currently in use
        // AND its level is effectively zero (decayed from previous use, or never used)
        if drugs_initiated_this_time_step < 2 && !individual.cur_use_drug[drug_idx] && individual.cur_level_drug[drug_idx] < 0.0001 {
            if rng.gen_bool(administration_prob) {
                individual.cur_use_drug[drug_idx] = true;

                // Determine the initial drug level, potentially with a double dose
                let mut chosen_initial_level = get_drug_param(drug_name, "initial_level").unwrap_or(10.0);
                
                // Apply double dose if infection is identified AND random chance passes
                if has_any_identified_infection && rng.gen_bool(double_dose_probability) {
                    let double_dose_multiplier = get_drug_param(drug_name, "double_dose_multiplier").unwrap_or(2.0);
                    chosen_initial_level *= double_dose_multiplier;
                    // println!("DEBUG: Individual {} initiated double dose of {} at day {}. Level: {:.2}", individual.id, drug_name, time_step, chosen_initial_level);
                } else {
                    // println!("DEBUG: Individual {} initiated drug {} at day {}. Level: {:.2}", individual.id, drug_name, time_step, chosen_initial_level);
                }

                individual.cur_level_drug[drug_idx] = chosen_initial_level;
                drugs_initiated_this_time_step += 1; // Increment the counter
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
                // in a newly infected person we should sample majority_r / any_r from all people in the same region with that 
                // bacteria and assign the newly infected person that level 
                let env_majority_r_level = get_global_param("environmental_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let hospital_majority_r_level = get_global_param("hospital_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0);

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
                        // When acquired from other individuals (neither env nor hospital):
                        // Sample a resistance level from the observed majority_r values in the population.
                        if let Some(majority_r_values_from_population) = majority_r_positive_values_by_combo.get(&(b_idx, d_idx)) {
                            if let Some(&acquired_resistance_level) = majority_r_values_from_population.choose(&mut rng) {
                                // If a resistant strain is acquired from the population,
                                // both any_r and majority_r start at that level.
                                let clamped_level = acquired_resistance_level.min(max_resistance_level).max(0.0); // Removed .round()
                                resistance_data.any_r = clamped_level;
                                resistance_data.majority_r = clamped_level; // Assign sampled level to majority_r too
                            } else {
                                // No majority-resistant strains observed in the population for this combo,
                                // so acquired strain is non-resistant.
                                resistance_data.any_r = 0.0;
                                resistance_data.majority_r = 0.0;
                            }
                        } else {
                            // No data for this bacteria-drug combo in majority_r_positive_values_by_combo,
                            // assume non-resistant acquisition.
                            resistance_data.any_r = 0.0;
                            resistance_data.majority_r = 0.0;
                        }
                    }
                }
                // --- END GENERALIZED any_r AND majority_r SETTING LOGIC ---

                individual.test_identified_infection.insert(bacteria, false);

                // Get the initial immunity level from config, or default to a reasonable value
                let initial_immunity = get_bacteria_param(bacteria, "initial_immunity_on_infection").unwrap_or(1.0); 
                individual.immune_resp.insert(bacteria, initial_immunity.max(0.0001)); // Ensure it starts above the floor

            } // End of if rng.gen_bool(acquisition_probability)
        } else { // Bacteria is already present (infection progression)
            // --- majority_r EVOLUTION LOGIC ---
            let majority_r_evolution_rate = get_global_param("majority_r_evolution_rate_per_day_when_drug_present").unwrap_or(0.0);
            let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0); // Now using 1.0 from your config

            if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                for (drug_index, _use_drug) in individual.cur_use_drug.iter().enumerate() { 
                    let resistance_data = &mut individual.resistances[bacteria_full_idx][drug_index];

                    let drug_current_level = individual.cur_level_drug[drug_index];
                    let drug_currently_present = drug_current_level > 0.0001; // Check if drug is effectively present
                    let current_bacteria_level = *individual.level.get(bacteria).unwrap_or(&0.0);

                    // Existing majority_r evolution based on drug presence
                    if resistance_data.majority_r == 0.0 && resistance_data.any_r > 0.0 && drug_currently_present {
                        if rng.gen_bool(majority_r_evolution_rate) {
                            resistance_data.majority_r = resistance_data.any_r;
                        }
                    }

                    // any_r persistence/decay logic
                    if resistance_data.majority_r > 0.0 {
                        // If majority_r exists, any_r is at least that level.
                        // For simplicity, we set any_r to majority_r here, assuming it encompasses.
                        resistance_data.any_r = resistance_data.majority_r;
                    } else {
                        // If no majority_r, any_r can decay or increase
                        let any_r_decay_rate = get_global_param("any_r_decay_rate_per_day").unwrap_or(0.0);
                        
                        // MODIFICATION 1: Prevent decay if drug is NOT currently present AND majority_r is 0
                        if !drug_currently_present {
                            resistance_data.any_r = (resistance_data.any_r - any_r_decay_rate).max(0.0);
                        }
                        // If drug IS present and majority_r is 0, any_r will not decay here.
                    }

                    // MODIFICATION 2: NEW RULE for any_r to increase towards max_resistance_level
                    // when drug is present and majority_r is still 0
                    if resistance_data.majority_r == 0.0 && // No majority resistance yet
                       resistance_data.any_r > 0.0 && // But some minority resistance exists
                       resistance_data.any_r < max_resistance_level && // And it's not yet full resistance
                       drug_currently_present // And the drug is present, providing selection pressure
                    {
                        let any_r_increase_rate = get_global_param("any_r_increase_rate_per_day_when_drug_present").unwrap_or(0.05); // New parameter
                        resistance_data.any_r = (resistance_data.any_r + any_r_increase_rate).min(max_resistance_level);
                    }


                    // Clamping majority_r and any_r
                    // MODIFICATION 3: Removed .round() for continuous float values
                    resistance_data.majority_r = resistance_data.majority_r.min(max_resistance_level).max(0.0);
                    resistance_data.any_r = resistance_data.any_r.min(max_resistance_level).max(0.0);


                    // --- NEW RESISTANCE EMERGENCE LOGIC (any_r from 0) ---
                    // This section handles the de novo emergence of resistance when it's not already present.
                    // It should come AFTER any_r/majority_r evolution/decay but BEFORE activity_r is fully calculated
                    // for use in bacteria level reduction *this* time step.
                    
                    if resistance_data.any_r < 0.0001 { // Check if any_r is effectively zero
                        // Only consider emergence if there's drug present (either being taken or decaying)
                        // and a positive bacteria level for selection pressure.
                        if drug_current_level > 0.0001 && current_bacteria_level > 0.0001 { 
                            let emergence_rate_baseline = get_global_param("resistance_emergence_rate_per_day_baseline").unwrap_or(0.000001); // Very small baseline
                            let bacteria_level_effect_multiplier = get_global_param("resistance_emergence_bacteria_level_multiplier").unwrap_or(0.05); // How much does bacteria level boost it
                            let any_r_emergence_level_on_first_emergence = get_global_param("any_r_emergence_level_on_first_emergence").unwrap_or(0.5); // User changed to 0.5 (was 1.0)

                            // 1. Bacteria Level Dependency: Higher at higher levels
                            let max_bacteria_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);
                            // Normalize bacteria level to [0,1] and apply multiplier
                            let bacteria_level_factor = (current_bacteria_level / max_bacteria_level).clamp(0.0, 1.0) * bacteria_level_effect_multiplier;
                            
                            // 2. Activity_r Dependency: Bell-shaped curve
                            // Use the drug's initial level for normalization to get a comparable 'activity' scale (0-10)
                            let drug_initial_level_for_normalization = get_drug_param(DRUG_SHORT_NAMES[drug_index], "initial_level").unwrap_or(10.0);
                            
                            // Normalized current drug level as a proxy for 'activity_r' when any_r is 0.
                            let mut norm_drug_level = drug_current_level / drug_initial_level_for_normalization;
                            norm_drug_level = norm_drug_level.clamp(0.0, 10.0); 
                            
                            // Bell curve: 4 * x * (1 - x). Peaks at 0.5, is 0 at 0 and 1.
                            let activity_r_bell_curve_factor = 0.02 * norm_drug_level * (10.0 - norm_drug_level);
                            let final_activity_r_factor = activity_r_bell_curve_factor.clamp(0.0, 1.0);   

                            // Total Emergence Probability
                            // Adding 1.0 to bacteria_level_factor ensures a base contribution even if multiplier is low
                            let total_emergence_prob = emergence_rate_baseline * (1.0 + bacteria_level_factor) * final_activity_r_factor;

                            if rng.gen_bool(total_emergence_prob.clamp(0.0, 1.0)) {
                                resistance_data.any_r = any_r_emergence_level_on_first_emergence;
                                // println!("DEBUG: Resistance to {} for {} emerged in Individual {} at day {}. New any_r: {:.2}",
                                //           DRUG_SHORT_NAMES[drug_index], bacteria, individual.id, time_step, resistance_data.any_r);
                            }
                        }
                    }
                    // --- END NEW RESISTANCE EMERGENCE LOGIC ---


                    // Calculate activity_r (this now uses the potentially newly emerged any_r)
                    // MODIFICATION 4: Simplification of normalized_any_r based on `max_resistance_level`
                    if drug_current_level > 0.0 {
                        // Use max_resistance_level from global params
                        let normalized_any_r = resistance_data.any_r / max_resistance_level; 
                        resistance_data.activity_r = drug_current_level * (1.0 - normalized_any_r);
                    } else {
                        resistance_data.activity_r = 0.0;
                    }
                }
            }

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
                let immunity_level = individual.immune_resp.get(bacteria).unwrap_or(&0.0);
                let baseline_change = get_bacteria_param(bacteria, "level_change_rate_baseline").unwrap_or(0.0);
                let reduction_due_to_immune_resp = get_bacteria_param(bacteria, "immunity_effect_on_level_change").unwrap_or(0.0);

                let mut total_reduction_due_to_antibiotic = 0.0;
                for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                    if individual.cur_level_drug[drug_idx] > 0.0 { // Only apply reduction if drug is actually present in the system
                        let bacteria_level_reduction_per_unit_of_drug = get_bacteria_drug_param(
                            bacteria,
                            drug_name, // Pass the actual drug name string
                            "bacteria_level_reduction_per_unit_of_drug"
                        ).unwrap_or(0.0);
                        total_reduction_due_to_antibiotic += individual.cur_level_drug[drug_idx] * bacteria_level_reduction_per_unit_of_drug;
                    }
                }

                let decay_rate = baseline_change - (immunity_level * reduction_due_to_immune_resp) - total_reduction_due_to_antibiotic;

                let max_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);

                let new_level = (current_level + decay_rate).max(0.0).min(max_level);
                *level_entry.get_mut() = new_level;
                // If bacteria level drops below a threshold, clear the infection
                if *level_entry.get() < 0.0001 {
                    individual.level.remove(bacteria);
                    individual.infectious_syndrome.remove(bacteria);
                    individual.date_last_infected.remove(bacteria);
                    individual.immune_resp.remove(bacteria);
                    individual.sepsis.remove(bacteria);
                    individual.level_microbiome.remove(bacteria);
                    individual.infection_hospital_acquired.remove(bacteria);
                    individual.cur_infection_from_environment.remove(bacteria);
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
                    *immune_entry.get_mut() = (*immune_entry.get() + immune_increase).max(0.0001);
                }
            }
        }
        // Immunity decay (applies regardless of current infection status)
        if let Entry::Occupied(mut immune_entry) = individual.immune_resp.entry(bacteria) {
            let current_immunity = *immune_entry.get();
            let baseline_immunity = get_bacteria_param(bacteria, "baseline_immunity_level").unwrap_or(0.0001);
            let decay_rate = get_bacteria_param(bacteria, "immunity_decay_rate").unwrap_or(0.0);

            if current_immunity > baseline_immunity {
                *immune_entry.get_mut() = (current_immunity - decay_rate).max(baseline_immunity);
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
    
/* if individual.id == 0 {
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

/* // Print Drug Levels for Individual 0 (for debugging/observation)
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