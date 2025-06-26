// src/rules/mod.rs

use crate::simulation::population::{Individual, Population, BACTERIA_LIST, DRUG_SHORT_NAMES, HospitalStatus, Region}; 
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
    // need a variable for vulnerability to serious toxicity ?
    individual.age += 1;
//  individual.current_infection_related_death_risk = rng.gen_range(0.0..=0.001);
//  individual.background_all_cause_mortality_rate = rng.gen_range(0.0..=0.00001);






    // --- NEW: Update Contact and Exposure Levels ---
    // Get general parameters for fluctuations and bounds
    let daily_fluctuation = get_global_param("contact_level_daily_fluctuation_range").unwrap_or(0.5);
    let min_contact_level = get_global_param("min_contact_level").unwrap_or(0.0);
    let max_contact_level = get_global_param("max_contact_level").unwrap_or(10.0);

    // Helper closure for applying fluctuation and clamping
    // This calculates a 'target' or 'base' level, then adds noise and clamps it.
    let mut update_contact_level = |current_level: &mut f64, base_level: f64| {
        *current_level = base_level + rng.gen_range(-daily_fluctuation..=daily_fluctuation);
        *current_level = current_level.clamp(min_contact_level, max_contact_level);
    };

    // 1. Sexual Contact Level
    let sexual_contact_age_peak_days = get_global_param("sexual_contact_age_peak_days").unwrap_or(25.0 * 365.0);
    let sexual_contact_age_decline_rate = get_global_param("sexual_contact_age_decline_rate").unwrap_or(0.00005);
    let sexual_contact_hospital_multiplier = get_global_param("sexual_contact_hospital_multiplier").unwrap_or(0.0); // Typically very low in hospital

    let mut base_sexual_level = get_global_param("sexual_contact_baseline").unwrap_or(5.0);
    if (individual.age as f64) < sexual_contact_age_peak_days {
        // Increase towards peak, but don't exceed baseline before peak
        base_sexual_level *= ((individual.age as f64 / sexual_contact_age_peak_days).min(1.0)).powf(get_global_param("sexual_contact_age_rise_exponent").unwrap_or(2.0));
    } else {
        // Decline after peak
        base_sexual_level *= (1.0 - (individual.age as f64 - sexual_contact_age_peak_days) * sexual_contact_age_decline_rate).max(0.0);
    }

    if individual.hospital_status.is_hospitalized() {
        base_sexual_level *= sexual_contact_hospital_multiplier;
    }
    update_contact_level(&mut individual.sexual_contact_level, base_sexual_level);


    // 2. Airborne Contact Level with Adults
    let airborne_adult_baseline = get_global_param("airborne_contact_adult_baseline").unwrap_or(5.0);
    let airborne_adult_age_breakpoint_days = get_global_param("airborne_contact_adult_age_breakpoint_days").unwrap_or(18.0 * 365.0); // 18 years old
    let airborne_in_hospital_multiplier = get_global_param("airborne_contact_in_hospital_multiplier").unwrap_or(1.5); // Might increase in hospital (e.g., healthcare workers)
    let airborne_adult_child_multiplier = get_global_param("airborne_contact_adult_child_multiplier").unwrap_or(0.2); // How much less children contact adults

    let mut base_airborne_adult_level = airborne_adult_baseline;
    if (individual.age as f64) < airborne_adult_age_breakpoint_days {
        base_airborne_adult_level *= airborne_adult_child_multiplier; // Children have less adult contact
    }
    if individual.hospital_status.is_hospitalized() {
        base_airborne_adult_level *= airborne_in_hospital_multiplier;
    }
    update_contact_level(&mut individual.airborne_contact_level_with_adults, base_airborne_adult_level);


    // 3. Airborne Contact Level with Children
    let airborne_child_baseline = get_global_param("airborne_contact_child_baseline").unwrap_or(3.0);
    let airborne_child_age_breakpoint_days = get_global_param("airborne_contact_child_age_breakpoint_days").unwrap_or(12.0 * 365.0); // 12 years old
    let airborne_child_adult_multiplier = get_global_param("airborne_contact_child_adult_multiplier").unwrap_or(0.5); // How much less adults contact children (than children contact children)

    let mut base_airborne_child_level = airborne_child_baseline;
    if (individual.age as f64) < airborne_child_age_breakpoint_days {
        // Higher for children interacting with children
        base_airborne_child_level *= get_global_param("airborne_contact_child_child_multiplier").unwrap_or(1.5);
    } else {
        // Lower for adults interacting with children (e.g., parents/teachers)
        base_airborne_child_level *= airborne_child_adult_multiplier;
    }
    if individual.hospital_status.is_hospitalized() {
        base_airborne_child_level *= airborne_in_hospital_multiplier; // Same multiplier for simplicity
    }
    update_contact_level(&mut individual.airborne_contact_level_with_children, base_airborne_child_level);


    // 4. Oral Exposure Level
    let oral_exposure_baseline = get_global_param("oral_exposure_baseline").unwrap_or(2.0);
    let oral_exposure_child_age_breakpoint_days = get_global_param("oral_exposure_child_age_breakpoint_days").unwrap_or(5.0 * 365.0); // 5 years old
    let oral_exposure_child_multiplier = get_global_param("oral_exposure_child_multiplier").unwrap_or(3.0);
    let oral_exposure_in_hospital_multiplier = get_global_param("oral_exposure_in_hospital_multiplier").unwrap_or(0.8); // Might decrease in hospital due to hygiene

    let mut base_oral_level = oral_exposure_baseline;
    if (individual.age as f64) < oral_exposure_child_age_breakpoint_days {
        base_oral_level *= oral_exposure_child_multiplier;
    }
    if individual.hospital_status.is_hospitalized() {
        base_oral_level *= oral_exposure_in_hospital_multiplier;
    }
    update_contact_level(&mut individual.oral_exposure_level, base_oral_level);


    // 5. Mosquito Exposure Level
    let mosquito_exposure_baseline = get_global_param("mosquito_exposure_baseline").unwrap_or(1.0);
    let mosquito_exposure_in_hospital_multiplier = get_global_param("mosquito_exposure_in_hospital_multiplier").unwrap_or(0.2); // Usually lower indoors/hospitals

    let mut base_mosquito_level = mosquito_exposure_baseline;

    // Apply region-specific multiplier
    // Convert Region enum to a lowercase, snake_case string for parameter lookup
    let region_name_for_param = individual.region_cur_in.to_string().to_lowercase().replace(" ", "_");
    let region_multiplier_key = format!("{}_mosquito_exposure_multiplier", region_name_for_param);
    let region_multiplier = get_global_param(&region_multiplier_key).unwrap_or(1.0); // Default to 1.0 if region not specified
    base_mosquito_level *= region_multiplier;

    if individual.hospital_status.is_hospitalized() {
        base_mosquito_level *= mosquito_exposure_in_hospital_multiplier;
    }
    update_contact_level(&mut individual.mosquito_exposure_level, base_mosquito_level);

    // --- END Update Contact and Exposure Levels ---





       // NEW: Update 'is_severely_immunosuppressed' status based on onset/recovery rates
    let onset_rate = get_global_param("immunosuppression_onset_rate_per_day").unwrap_or(0.0001);
    let recovery_rate = get_global_param("immunosuppression_recovery_rate_per_day").unwrap_or(0.0005);

    if individual.is_severely_immunosuppressed {
        // If currently immunosuppressed, check for recovery
        if rng.gen_bool(recovery_rate) {
            individual.is_severely_immunosuppressed = false;
        }
    } else { 
        // If not immunosuppressed, check for onset
        if rng.gen_bool(onset_rate) {
            individual.is_severely_immunosuppressed = true;
        }
    }





    if rng.gen::<f64>() < 0.01 { individual.under_care = !individual.under_care; }
    individual.current_toxicity = (individual.current_toxicity + rng.gen_range(-0.5..=0.5)).max(0.0);
//  individual.mortality_risk_current_toxicity = rng.gen_range(0.0..=0.0001);

    // todo: we have this variable under_care - do we need this ? if so we need it to be a pre-requisite of being given drugs  




        // Get parameters from config.rs ONCE per individual for this time step
    let baseline_rate = get_global_param("hospitalization_baseline_rate_per_day")
        .expect("Missing hospitalization_baseline_rate_per_day in config");
    let age_multiplier = get_global_param("hospitalization_age_multiplier_per_day")
        .expect("Missing hospitalization_age_multiplier_per_day in config");
    let recovery_rate = get_global_param("hospitalization_recovery_rate_per_day")
        .expect("Missing hospitalization_recovery_rate_per_day in config");
    let max_days_in_hospital = get_global_param("hospitalization_max_days")
        .expect("Missing hospitalization_max_days in config");

    // Rule 1: Potentially get hospitalized (if not currently hospitalized)
    if !individual.hospital_status.is_hospitalized() { // Use the new helper method
        let prob_hospitalization_today = baseline_rate + (individual.age as f64 * age_multiplier);

        if rng.gen::<f64>() < prob_hospitalization_today {
            individual.hospital_status = HospitalStatus::InHospital; // Assign enum variant
            individual.days_hospitalized = 0; // Initialize days hospitalized
            // println!("DEBUG: Individual {} (Age: {}) hospitalized!", individual.id, individual.age); // Optional: For debugging
        }
    } else { // If already hospitalized, consider recovery or max days limit
        individual.days_hospitalized += 1; // Increment days hospitalized

        // Rule 2: Potentially recover from hospitalization
        if rng.gen::<f64>() < recovery_rate {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
            // println!("DEBUG: Individual {} recovered from hospitalization.", individual.id); // Optional: For debugging
        }
        // Rule 3: Forced discharge after max_days_in_hospital
        else if individual.days_hospitalized >= max_days_in_hospital as u32 {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
            // println!("DEBUG: Individual {} forcibly discharged after {} days.", individual.id, max_days_in_hospital); // Optional: For debugging
        }
    }
    // --- END Hospitalization Rules ---



       // --- NEW: Region Travel Rules ---
    let travel_prob = get_global_param("travel_probability_per_day")
        .expect("Missing travel_probability_per_day in config");
    const VISIT_LENGTH_DAYS: u32 = 30; // Fixed visit length

    // Check if the individual is currently in their home region
    if let Region::Home = individual.region_cur_in {
        // If not hospitalized, consider initiating travel
        if !individual.hospital_status.is_hospitalized() && rng.gen::<f64>() < travel_prob {
            // Initiate travel: select a random new region different from their living region
            let mut new_region: Region;
            loop {
                // rng.gen() for Region will give one of the 6 geographic regions (not Home)
                new_region = rng.gen();
                // Ensure the individual doesn't 'travel' to their own living region
                if new_region != individual.region_living {
                    break; // Found a suitable new region to visit
                }
            }
            individual.region_cur_in = new_region;
            individual.days_visiting = 1; // Start the visit counter at 1
            // Optional: for debugging
            // println!("DEBUG (Day {}): Individual {} (Age: {}) traveled from {:?} to {:?}",
            //     time_step, individual.id, individual.age, individual.region_living, individual.region_cur_in);
        }
    } else {
        // Individual is currently visiting another region
        individual.days_visiting += 1; // Increment the visit duration

        // Check if the visit duration has been reached
        if individual.days_visiting >= VISIT_LENGTH_DAYS {
            // End of visit, return to home region
            individual.region_cur_in = Region::Home; // Set current region back to Home
            individual.days_visiting = 0; // Reset visit counter
            // Optional: for debugging
            // println!("DEBUG (Day {}): Individual {} (Age: {}) returned home from a trip.",
            //     time_step, individual.id, individual.age);
        }
    }
    // --- END Region Travel Rules ---




        // --- NEW: Sepsis Risk Rules ---
    for &bacteria in BACTERIA_LIST.iter() {
        if let Some(&current_level) = individual.level.get(bacteria) {
            // Only consider sepsis if an infection is present (level > 0.0)
            if current_level > 0.0 {
                if let Some(&last_infected_day) = individual.date_last_infected.get(bacteria) {
                    let duration_of_infection = (time_step as i32 - last_infected_day).max(0); // Ensure non-negative duration

                    // Retrieve bacteria-specific parameters, falling back to global defaults
                    let sepsis_baseline_risk = get_bacteria_param(bacteria, "sepsis_baseline_risk_per_day")
                        .unwrap_or_else(|| get_global_param("default_sepsis_baseline_risk_per_day").expect("Missing default_sepsis_baseline_risk_per_day"));
                    let sepsis_level_multiplier = get_bacteria_param(bacteria, "sepsis_level_multiplier")
                        .unwrap_or_else(|| get_global_param("default_sepsis_level_multiplier").expect("Missing default_sepsis_level_multiplier"));
                    let sepsis_duration_multiplier = get_bacteria_param(bacteria, "sepsis_duration_multiplier")
                        .unwrap_or_else(|| get_global_param("default_sepsis_duration_multiplier").expect("Missing default_sepsis_duration_multiplier"));

                    // Calculate daily probability of sepsis
                    let prob_sepsis_today = sepsis_baseline_risk
                                            + (current_level * sepsis_level_multiplier)
                                            + (duration_of_infection as f64 * sepsis_duration_multiplier);

                    // Cap the probability at 1.0
                    let prob_sepsis_today = prob_sepsis_today.min(1.0);

                    // Roll the dice for sepsis
                    if rng.gen::<f64>() < prob_sepsis_today {
                        // Set sepsis status to true for this bacteria
                        individual.sepsis.insert(bacteria, true);
                        // Optional: For debugging
                        // println!("DEBUG (Day {}): Individual {} developed sepsis for {} (Level: {:.2}, Duration: {})",
                        //     time_step, individual.id, bacteria, current_level, duration_of_infection);
                    }
                }
            } else {
                // If infection level is 0 or individual is hospitalized, assume sepsis resolves for this bacteria
                // This assumes sepsis is directly tied to active infection and is managed in hospital.
                // You might need a more complex resolution mechanism if sepsis can persist after infection or hospitalization.
                if individual.sepsis.get(bacteria).copied().unwrap_or(false) {
                    individual.sepsis.insert(bacteria, false);
                }
            }
        }
    }
    // --- END Sepsis Risk Rules ---




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

    // --- Drug Stopping Logic ---
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

            let mut stop_drug = false; // Assume not stopping by default

            // Condition to stop the drug:
            // 1. No relevant active infection (primary reason) OR
            // 2. A small random chance to stop early (e.g., non-adherence, side effects)
            if !relevant_infection_active_for_this_drug || rng.gen_bool(random_drug_cessation_prob) {
                stop_drug = true;
            }

            // NEW RULE: Prevent stopping if the drug was initiated in the *previous* time_step.
            // This ensures it is administered for at least the current day (which is its first full day post-initiation).
            if individual.date_drug_initiated[drug_idx] == (time_step as i32) - 1 {
                stop_drug = false; // Override any decision to stop for this drug
            }


            if stop_drug {
                individual.cur_use_drug[drug_idx] = false;
                // Reset the initiation date to indicate the drug is no longer being taken.
                // This is important to allow re-initiation later and for the stopping rule to re-apply correctly.
                individual.date_drug_initiated[drug_idx] = i32::MIN; 
                // println!("DEBUG: Individual {} stopped drug {} at day {}. No relevant infection: {}. Random stop: {}",
                //      individual.id, drug_name, time_step, !relevant_infection_active_for_this_drug, rng.gen_bool(random_drug_cessation_prob));
            }
        }
    }

    // Apply decay if stopped, or set to initial level if continued/re-initiated.
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


    // --- New Drug Initiation Logic ---
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
        if drugs_initiated_this_time_step < 2 && !individual.cur_use_drug[drug_idx] {
            if rng.gen_bool(administration_prob) {
                individual.cur_use_drug[drug_idx] = true;
                individual.date_drug_initiated[drug_idx] = time_step as i32; // ADDED: Record initiation time

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


        // NEW: Add drug-specific toxicity contributions
    let mut daily_drug_toxicity_increase = 0.0;
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];
        // Only consider toxicity if the drug is currently present in the system
        if individual.cur_level_drug[drug_idx] > 0.0 {
            // Get the drug-specific toxicity parameter.
            // Falls back to a global default if the drug-specific one isn't found.
            let drug_toxicity_per_unit = get_drug_param(drug_name, "toxicity_per_unit_level_per_day")
                .unwrap_or_else(|| {
                    get_global_param("default_drug_toxicity_per_unit_level_per_day")
                        .expect("Missing default_drug_toxicity_per_unit_level_per_day in config")
                });
            
            daily_drug_toxicity_increase += individual.cur_level_drug[drug_idx] * drug_toxicity_per_unit;
        }
    }
    // Add the calculated drug-induced toxicity to the individual's current_toxicity
    individual.current_toxicity = (individual.current_toxicity + daily_drug_toxicity_increase).max(0.0);



    // --- DRUG LOGIC END ---






    // --- DEATH LOGIC START ---
    // Only apply death logic if the individual is not already dead
    if individual.date_of_death.is_none() {
        let mut prob_of_death_today = 0.0;
        let mut cause: Option<String> = None;

        // 1. Background Mortality Risk (Age, Region, and Sex dependent)
        let base_background_rate = get_global_param("base_background_mortality_rate_per_day")
            .expect("Missing base_background_mortality_rate_per_day in config");
        let age_multiplier = get_global_param("age_mortality_multiplier_per_year")
            .expect("Missing age_mortality_multiplier_per_year in config");

        // needs some attention as want age_multiplier to be a rate ratio per additional day of age - and at some point to recognise non-linearity of age effect    
        let mut background_risk = base_background_rate;
        background_risk += (individual.age as f64 / 365.0) * age_multiplier; 

        // Apply region-specific multiplier
        let region_multiplier_key = format!("{}_mortality_multiplier", individual.region_living.to_string().to_lowercase().replace(" ", "_"));
        let region_multiplier = get_global_param(&region_multiplier_key)
            .unwrap_or(1.0); // Default to 1.0 if not found
        background_risk *= region_multiplier;

        // Apply sex-specific multiplier
        let sex_multiplier_key = format!("{}_mortality_multiplier", individual.sex_at_birth.to_lowercase());
        let sex_multiplier = get_global_param(&sex_multiplier_key)
            .unwrap_or(1.0); // Default to 1.0 if not found
        background_risk *= sex_multiplier;

        individual.background_all_cause_mortality_rate = background_risk.min(1.0); // Store for potential logging/debugging

        // Initialize total probability of NOT dying
        let mut prob_not_dying = 1.0 - background_risk;

        // 2. Absolute Risk of Death from Sepsis (Independent)
        let has_sepsis = individual.sepsis.values().any(|&status| status);
        if has_sepsis {
            let sepsis_absolute_death_risk = get_global_param("sepsis_absolute_death_risk_per_day")
                .expect("Missing sepsis_absolute_death_risk_per_day in config");
            prob_not_dying *= 1.0 - sepsis_absolute_death_risk;
            if cause.is_none() { cause = Some("sepsis_related".to_string()); } // Prioritize sepsis as a cause
        }

        // 3. Absolute Risk of Death from Drug Adverse Events (Independent and Drug-specific)
        let mut drug_adverse_event_risk_for_individual = 0.0;
        for drug_idx in 0..DRUG_SHORT_NAMES.len() {
            let drug_name = DRUG_SHORT_NAMES[drug_idx];
            if individual.cur_level_drug[drug_idx] > 0.0 { // If drug is present in system
                let drug_adverse_event_risk = get_drug_param(drug_name, "adverse_event_death_risk")
                    .unwrap_or(0.0); // Default to 0 if not specified for drug
                drug_adverse_event_risk_for_individual = (drug_adverse_event_risk_for_individual + drug_adverse_event_risk).min(1.0);
            }
        }
        individual.mortality_risk_current_toxicity = drug_adverse_event_risk_for_individual; // Store for potential logging/debugging
        if drug_adverse_event_risk_for_individual > 0.0 {
            prob_not_dying *= (1.0 - drug_adverse_event_risk_for_individual);
            if cause.is_none() { cause = Some("drug_toxicity_related".to_string()); } // Prioritize drug toxicity if no sepsis
        }

        // REMOVED: 4. Existing Infection-Related Death Risk (No longer a direct cause, only through sepsis)
        // If you still have a current_infection_related_death_risk field, its value won't be used here for death.
        // if individual.current_infection_related_death_risk > 0.0 {
        //     prob_not_dying *= (1.0 - individual.current_infection_related_death_risk.min(1.0));
        //     if cause.is_none() { cause = Some("infection_related".to_string()); }
        // }

        // Calculate final probability of death
        prob_of_death_today = 1.0 - prob_not_dying;
        prob_of_death_today = prob_of_death_today.clamp(0.0, 1.0); // Ensure it's between 0 and 1

        // Roll the dice for death
        if rng.gen::<f64>() < prob_of_death_today {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = cause.or(Some("background_mortality".to_string())); // Default to background if no specific cause assigned
            // println!("DEBUG: Individual {} died at time step {} due to: {:?}",
            //     individual.id, time_step, individual.cause_of_death); // Optional: For debugging
        }
    }
    // --- DEATH LOGIC END ---






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

        // todo: these will depend on bacteria 
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

        // NEW: Apply microbiome presence effect on infection acquisition risk
        // If the bacteria is already in the individual's microbiome (carriage)
        if individual.presence_microbiome.get(bacteria).copied().unwrap_or(false) {
            let microbiome_infection_multiplier = get_bacteria_param(bacteria, "microbiome_infection_acquisition_multiplier")
                .unwrap_or_else(|| {
                    get_global_param("default_microbiome_infection_acquisition_multiplier")
                        .expect("Missing default_microbiome_infection_acquisition_multiplier in config")
                });
            acquisition_probability *= microbiome_infection_multiplier;
            // Optional: for debugging
            // println!("DEBUG (Day {}): Individual {} has {} in microbiome, infection acquisition_probability modified to {:.4}",
            //     time_step, individual.id, bacteria, acquisition_probability);
        }
        // --- END BACTERIA-SPECIFIC ACQUISITION PROBABILITY CALCULATION ---




 // --- NEW: Update Microbiome Presence (Carriage) ---
// This block should be placed for each 'bacteria' where 'acquisition_probability' is known.

// If the individual does NOT currently have this bacteria in their microbiome
if !individual.presence_microbiome.get(bacteria).copied().unwrap_or(false) {
    // Get a specific multiplier for microbiome acquisition.
    // This allows the absolute risk for microbiome acquisition to differ from infection.
    let microbiome_acquisition_multiplier = get_bacteria_param(bacteria, "microbiome_acquisition_multiplier")
        .unwrap_or_else(|| {
            get_global_param("default_microbiome_acquisition_multiplier")
                .expect("Missing default_microbiome_acquisition_multiplier in config")
        });

    // Calculate the microbiome-specific acquisition probability.
    // It depends on the same underlying factors as infection (via acquisition_probability),
    // but its absolute value is scaled by the new multiplier.
    let microbiome_acquisition_probability = acquisition_probability * microbiome_acquisition_multiplier;

    if rng.gen_bool(microbiome_acquisition_probability.clamp(0.0, 1.0)) {
        individual.presence_microbiome.insert(bacteria, true);
        // Optional: for debugging
        // println!("DEBUG (Day {}): Individual {} acquired {} in microbiome (carriage).", time_step, individual.id, bacteria);
    }
} else {
    // If the individual *does* have this bacteria in their microbiome,
    // there should be a chance to lose it (clearance).
    let microbiome_clearance_prob = get_bacteria_param(bacteria, "microbiome_clearance_probability_per_day")
        .unwrap_or_else(|| {
            get_global_param("default_microbiome_clearance_probability_per_day")
                .expect("Missing default_microbiome_clearance_probability_per_day in config")
        });

    if rng.gen_bool(microbiome_clearance_prob) {
        individual.presence_microbiome.insert(bacteria, false);
        // Optional: for debugging
        // println!("DEBUG (Day {}): Individual {} lost {} from microbiome (carriage).", time_step, individual.id, bacteria);
    }
}
// --- END Microbiome Presence Rules ---




        

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

            // review how hospital acquired different - do we need to treat environment differently for people in hospital ? 
            let hospital_acquired_chance = get_bacteria_param(bacteria, "hospital_acquired_proportion").unwrap_or(0.05);
            let mut is_hospital_acquired = false; // Initialize to false

            // NEW LOGIC: Only consider hospital-acquired if the individual is currently hospitalized
            if individual.hospital_status.is_hospitalized() {
                is_hospital_acquired = rng.gen::<f64>() < hospital_acquired_chance;
            }
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
                            let clamped_level = acquired_resistance_level.min(max_resistance_level).max(0.0);
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

                // this code needs correcting - value for any_r or majority_r for any drug bacteria combination will
                // not decline so long as the bacterial infection is present - even after bacterial infection
                // has gone it may be in microbiome      

                if resistance_data.majority_r > 0.0 {
                // If majority_r exists, any_r is at least that level.
                // For simplicity, we set any_r to majority_r here, assuming it encompasses.
                   resistance_data.any_r = resistance_data.majority_r;
                } else {
                // If no majority_r, any_r should only increase or stay stable while infected.
                // The decay logic has been moved to the infection clearance block.
                // So, nothing to do here in terms of decay.
                // It can only increase via the 'any_r increase towards max_resistance_level' block below.
                }
                // any_r increase towards max_resistance_level
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
                resistance_data.majority_r = resistance_data.majority_r.min(max_resistance_level).max(0.0);
                resistance_data.any_r = resistance_data.any_r.min(max_resistance_level).max(0.0);


                //new resistance emergence ---
                // This section handles the de novo emergence of resistance when it's not already present.
                // It should come before activity_r is fully calculated for use in bacteria level reduction *this* time step.
                
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
                        
                        // Bell curve: 0.02 * x * (10 - x). Peaks at 5.0, is 0 at 0 and 10.
                        let activity_r_bell_curve_factor = 0.02 * norm_drug_level * (10.0 - norm_drug_level);
                        let final_activity_r_factor = activity_r_bell_curve_factor.clamp(0.0, 1.0);  

                        // Total Emergence Probability
                        // Adding 1.0 to bacteria_level_factor ensures a base contribution even if multiplier is low
                        let total_emergence_prob = emergence_rate_baseline * (1.0 + bacteria_level_factor) * final_activity_r_factor;

                        if rng.gen_bool(total_emergence_prob.clamp(0.0, 1.0)) {
                            resistance_data.any_r = any_r_emergence_level_on_first_emergence;
                            // println!("DEBUG: Resistance to {} for {} emerged in Individual {} at day {}. New any_r: {:.2}",
                            //      DRUG_SHORT_NAMES[drug_index], bacteria, individual.id, time_step, resistance_data.any_r);
                        }
                    }
                }
                // --- END NEW RESISTANCE EMERGENCE LOGIC ---


                // Calculate activity_r (this now uses the potentially newly emerged any_r)
                // Simplification of normalized_any_r based on `max_resistance_level`
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
            // We only need to check if individual.cur_level_drug[drug_idx] > 0.0
            // because activity_r will be 0.0 if the drug level is 0.0 according to your formula above.
            if individual.cur_level_drug[drug_idx] > 0.0 {
            // Retrieve the activity_r for this specific bacteria-drug combination for the individual
            // (This activity_r would have been calculated earlier in the rules for this time step)
            let b_idx = bacteria_indices.get(bacteria).expect("Bacteria not found in indices");
            let resistance_data = &individual.resistances[*b_idx][drug_idx];
            let activity_r_for_this_combo = resistance_data.activity_r;

            // Add activity_r directly, as it already encapsulates drug level and resistance
            total_reduction_due_to_antibiotic += activity_r_for_this_combo;
            }
            }

            let decay_rate = baseline_change - (immunity_level * reduction_due_to_immune_resp) - total_reduction_due_to_antibiotic;

            let max_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);

            let new_level = (current_level + decay_rate).max(0.0).min(max_level);
            *level_entry.get_mut() = new_level;

            // If bacteria level drops below a threshold, clear the infection
            if *level_entry.get() < 0.0001 {
                // --- THIS IS THE KEY AREA FOR THE CHANGE ---

                // Immediately set any_r and majority_r to 0.0 for this bacteria across all drugs
                if let Some(b_idx_clear) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                    for drug_idx_clear in 0..DRUG_SHORT_NAMES.len() {
                        // Access the mutable resistance data
                        let resistance_data = &mut individual.resistances[b_idx_clear][drug_idx_clear];
                        resistance_data.any_r = 0.0;
                        resistance_data.majority_r = 0.0;
                        // Also reset activity_r, as it depends on any_r
                        resistance_data.activity_r = 0.0; // Ensure activity_r is also zeroed out
                    }
                }

                // Now proceed with removing other infection-related data
                individual.level.remove(bacteria);
                individual.infectious_syndrome.remove(bacteria);
                individual.date_last_infected.remove(bacteria);
                individual.immune_resp.remove(bacteria);
                individual.sepsis.remove(bacteria);
                individual.presence_microbiome.remove(bacteria); // If you're not considering microbiome persistence for resistance
                individual.infection_hospital_acquired.remove(bacteria);
                individual.cur_infection_from_environment.remove(bacteria);
                individual.test_identified_infection.insert(bacteria, false);

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
// if individual.id == 0 {
//     println!("-------------------------------------");
//     println!("--- Resistance Status (Individual 0) ---");
//     let has_relevant_resistance = individual.resistances.iter().any(|b_res| {
//         b_res.iter().any(|res| res.any_r > 0.0 || res.majority_r > 0.0)
//     });

//     if has_relevant_resistance {
//         for (bacteria_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
//             for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
//                 if let Some(resistance) = individual.resistances.get(bacteria_idx).and_then(|r_vec| r_vec.get(drug_idx)) {
//                     if resistance.any_r > 0.0 || resistance.majority_r > 0.0 {
//                         println!("     {} resistance to {}:", bacteria_name, drug_name);
//                         println!("       any_r: {:.0}", resistance.any_r);
//                         println!("       majority_r: {:.0}", resistance.majority_r);
//                         println!("       activity_r: {:.4}", resistance.activity_r);
//                     }
//                 }
//             }
//         }
//     } else {
//         println!("     No active resistances (any_r > 0 or majority_r > 0) for Individual 0.");
//     }
//     println!("-------------------------------------");
// }

// Print Drug Levels for Individual 0 (for debugging/observation)
// if individual.id == 0 {
//     println!("--- Drug Levels (Individual 0) ---");
//     for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
//         if individual.cur_level_drug[drug_idx] > 0.0 || individual.cur_use_drug[drug_idx] {
//             println!("     Drug {}: cur_level={:.2}, cur_use={}",
//                 drug_name,
//                 individual.cur_level_drug[drug_idx],
//                 individual.cur_use_drug[drug_idx]
//             );
//             // Optionally print date_drug_initiated for debugging
//             // println!("       Date Initiated: {}", individual.date_drug_initiated[drug_idx]);
//         }
//     }
//     println!("-------------------------------------");
// }
} // This closing brace seems to be for the overall function or `impl` block