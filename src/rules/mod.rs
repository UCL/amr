//
// Core model rules and update logic for the AMR simulation.
//
// Contains:
//   - apply_rules: main function to update an individual's state for one time step
//   - Logic for resistance emergence, drug effects, MIC calculation, cross-resistance, HGT, and reversion
//   - Helper functions for parameter lookup and stochastic events
//
// src/rules/mod.rs\
//


// for printing individual 0 per time step replace .id == 1000001 with .id == 1000001 (cntrl h to find and replace)


use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES, HospitalStatus, Region, 
                                   InfectionResolutionType, ImmunodeficiencyType}; 
use crate::config::{
    get_global_param,
    get_bacteria_param,
    get_drug_param,
    get_age_infection_multiplier,
    get_drug_availability_time_aware,
    get_bacteria_sepsis_risk_multiplier,
    get_age_dependent_bacteria_sepsis_risk_multiplier,
    get_drug_introduction_time_step,
};
use rand::Rng;
use rand::seq::SliceRandom;
use std::collections::HashMap;

/// Helper function to update the current number of drugs counter
fn update_drug_counter(individual: &mut Individual) {
    individual.current_number_of_drugs = individual.cur_use_drug.iter().filter(|&&on| on).count() as i32;
}
use rand::distributions::WeightedIndex;
use rand::distributions::Distribution; 

/// Assess treatment failure and switch drugs if necessary
/// Returns true if a drug switch occurred
fn assess_treatment_failure(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    bacteria_indices: &HashMap<&'static str, usize>,
    _drug_indices: &HashMap<&'static str, usize>,
    _cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>,
    param_cache: &ParameterKeyCache,
) -> bool {
    let mut rng = rand::thread_rng();
    
    // Check if treatment failure assessment is enabled
    let assessment_enabled = get_global_param("enable_treatment_failure_assessment").unwrap_or(1.0) > 0.5;
    if !assessment_enabled {
        return false;
    }
    
    // Check if we're on the assessment day
    let assessment_day = get_global_param("treatment_failure_assessment_day").unwrap_or(4.0) as i32;
    if individual.days_on_current_treatment[bacteria_idx] != assessment_day {
        return false;
    }
    
    // Check if we've already assessed this treatment course
    if individual.treatment_failure_assessed[bacteria_idx] {
        return false;
    }
    
    // Check if there's a current infection and bacteria level recorded at drug start
    if individual.level[bacteria_idx] <= 0.0 || individual.bacteria_level_at_drug_start[bacteria_idx].is_none() {
        return false;
    }
    
    let initial_level = individual.bacteria_level_at_drug_start[bacteria_idx].unwrap();
    let current_level = individual.level[bacteria_idx];
    
    // Get failure threshold (default 0.5 = 50% of initial level)
    let failure_threshold = get_global_param("treatment_failure_threshold").unwrap_or(0.5);
    let threshold_level = initial_level * failure_threshold;
    
    // Treatment failure criterion: current bacteria level >= threshold × initial level
    let treatment_failed = current_level >= threshold_level;
    
    // Mark assessment as completed for this treatment course
    individual.treatment_failure_assessed[bacteria_idx] = true;
    
    if !treatment_failed {
        return false; // Treatment is working, no switch needed
    }
    
    // Record drug failure date for this bacteria
    individual.date_last_drug_failure[bacteria_idx] = time_step as i32;
    
    // Find current drugs being used for this bacteria
    let current_drugs: Vec<usize> = individual.cur_use_drug.iter().enumerate()
        .filter(|(_, &is_taking)| is_taking)
        .map(|(drug_idx, _)| drug_idx)
        .collect();
    
    if current_drugs.is_empty() {
        return false; // No current drugs to switch from
    }
    
    // Try to find an alternative drug using the same selection logic as initial prescription
    // but excluding recently failed drugs
    let bacteria_name = BACTERIA_LIST[bacteria_idx];
    let failure_memory_days = get_global_param("drug_failure_memory_days").unwrap_or(30.0) as i32;
    
    // Build list of available alternative drugs
    let mut alternative_scores = Vec::new();
    
    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        // Skip if currently taking this drug
        if current_drugs.contains(&drug_idx) {
            continue;
        }
        
        // Skip if this drug failed recently (within memory period)
        if individual.date_drug_initiated_keep[drug_idx] != i32::MIN {
            let days_since_last_use = (time_step as i32) - individual.date_drug_initiated_keep[drug_idx];
            if days_since_last_use >= 0 && days_since_last_use < failure_memory_days {
                // This is a recently used drug, skip it for now (simple approach)
                continue;
            }
        }
        
        // Check if drug is available (using same logic as original selection)
        let avail = get_drug_availability_time_aware(
            drug_name,
            &individual.region_cur_in.to_string(),
            Some(&individual.region_living.to_string()),
            time_step
        );
        
        if avail < 0.01 { // Drug not sufficiently available
            continue;
        }
        
        // Calculate drug score using same logic as original selection
        let mut score = 0.0;
        
        // Base potency score
        let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
        if let Some(potency_param_key) = param_cache.drug_bacteria_potency_keys.get(&(drug_idx, *bacteria_idx_for_cache)) {
            if let Some(potency) = get_global_param(potency_param_key) {
                if potency >= get_global_param("minimal_potency_threshold_for_drug_selection").unwrap_or(0.10) {
                    score += potency;
                }
            }
        }
        
        // Apply clinical multipliers (same as original logic)
        // Add pathogen-specific preference multipliers
        let bacteria_drug_key = format!("{}_{}_clinical_preference_multiplier", bacteria_name.replace(" ", "_"), drug_name);
        if let Some(preference_multiplier) = get_global_param(&bacteria_drug_key) {
            score *= preference_multiplier;
        }
        
        if score > 0.0 {
            alternative_scores.push((drug_idx, score));
        }
    }
    
    // If we found alternatives, select one and switch
    if !alternative_scores.is_empty() {
        // Use same weighted selection as original logic
        let selection_temperature = get_global_param("drug_selection_temperature").unwrap_or(0.5);
        let weights: Vec<f64> = alternative_scores.iter()
            .map(|(_, score)| (score / selection_temperature).exp())
            .collect();
        
        let total_weight: f64 = weights.iter().sum();
        if total_weight > 0.0 && total_weight.is_finite() {
            let dist = WeightedIndex::new(&weights).unwrap();
            let chosen_idx = dist.sample(&mut rng);
            let new_drug_idx = alternative_scores[chosen_idx].0;
            
            // Stop current drugs
            for &current_drug_idx in &current_drugs {
                individual.cur_use_drug[current_drug_idx] = false;
                individual.date_drug_initiated[current_drug_idx] = i32::MIN;
            }
            
            // Start new drug
            let new_drug_name = DRUG_SHORT_NAMES[new_drug_idx];
            individual.cur_use_drug[new_drug_idx] = true;
            individual.date_drug_initiated[new_drug_idx] = time_step as i32;
            individual.date_drug_initiated_keep[new_drug_idx] = time_step as i32;
            individual.ever_taken_drug[new_drug_idx] = true;
            
            // Update drug counter
            update_drug_counter(individual);
            
            // Set drug level
            let initial_level = get_drug_param(new_drug_name, "initial_level").unwrap_or(10.0);
            individual.cur_level_drug[new_drug_idx] = initial_level;
            
            // Reset treatment failure tracking for this bacteria
            individual.bacteria_level_at_drug_start[bacteria_idx] = Some(current_level);
            individual.days_on_current_treatment[bacteria_idx] = 0;
            individual.treatment_failure_assessed[bacteria_idx] = false;
            
            return true; // Drug switch occurred
        }
    }
    
    false // No switch occurred
}

/// Assess restart window for patients who stopped drugs while still infected
/// Returns true if restart treatment was initiated
fn assess_restart_window(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    bacteria_indices: &HashMap<&'static str, usize>,
    param_cache: &ParameterKeyCache,
) -> bool {
    let mut rng = rand::thread_rng();
    
    // Check if restart window is enabled
    if get_global_param("enable_restart_window").unwrap_or(1.0) < 0.5 {
        return false;
    }
    
    // Check if there's a cessation to assess
    if let Some(cessation_day) = individual.drug_stopped_with_infection_day[bacteria_idx] {
        let restart_window_days = get_global_param("restart_window_days").unwrap_or(5.0) as i32;
        let days_since_cessation = (time_step as i32) - cessation_day;
        
        // Within restart window?
        if days_since_cessation >= 1 && days_since_cessation <= restart_window_days {
            
            // Haven't assessed yet?
            if !individual.restart_window_assessed[bacteria_idx] {
                individual.restart_window_assessed[bacteria_idx] = true;
                
                // Check if bacteria level has worsened enough to trigger restart
                if let Some(cessation_level) = individual.bacteria_level_at_drug_cessation[bacteria_idx] {
                    let current_level = individual.level[bacteria_idx];
                    let threshold_multiplier = get_global_param("restart_bacteria_level_threshold").unwrap_or(1.5);
                    
                    // Restart criteria: bacteria level increased significantly OR still very high
                    let bacteria_worsened = current_level >= (cessation_level * threshold_multiplier);
                    let bacteria_still_high = current_level > 2.0; // Arbitrary high threshold for severe infection
                    
                    if (bacteria_worsened || bacteria_still_high) && individual.level[bacteria_idx] > 0.1 {
                        // Patient decides to return to care?
                        let return_probability = get_global_param("restart_window_probability").unwrap_or(0.3);
                        
                        if rng.gen_bool(return_probability) {
                            // Clear restart tracking
                            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
                            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
                            let stopped_drug_idx = individual.stopped_drug_index[bacteria_idx];
                            individual.stopped_drug_index[bacteria_idx] = None;
                            
                            // Start restart treatment, preferring the previously effective drug
                            return start_restart_treatment(individual, time_step, bacteria_idx, stopped_drug_idx, bacteria_indices, param_cache);
                        }
                    }
                }
            }
        } else if days_since_cessation > restart_window_days {
            // Restart window expired - clear tracking
            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
            individual.stopped_drug_index[bacteria_idx] = None;
            individual.restart_window_assessed[bacteria_idx] = false;
        }
    }
    
    false
}

/// Start restart treatment for a patient who returns to care after stopping drugs early
/// Prefers the previously effective drug that was stopped
fn start_restart_treatment(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    stopped_drug_idx: Option<usize>,
    bacteria_indices: &HashMap<&'static str, usize>,
    param_cache: &ParameterKeyCache,
) -> bool {
    let mut rng = rand::thread_rng();
    
    let bacteria_name = BACTERIA_LIST[bacteria_idx];
    
    // Check if we can restart the previously effective drug
    if let Some(prev_drug_idx) = stopped_drug_idx {
        let prev_drug_name = DRUG_SHORT_NAMES[prev_drug_idx];
        
        // Check if previously effective drug is still available
        let avail = get_drug_availability_time_aware(
            prev_drug_name,
            &individual.region_cur_in.to_string(),
            Some(&individual.region_living.to_string()),
            time_step
        );
        
        if avail >= 0.01 && !individual.cur_use_drug[prev_drug_idx] {
            // Check if drug has adequate potency (basic safety check)
            let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
            if let Some(potency_param_key) = param_cache.drug_bacteria_potency_keys.get(&(prev_drug_idx, *bacteria_idx_for_cache)) {
                if let Some(potency) = get_global_param(potency_param_key) {
                    if potency >= get_global_param("minimal_potency_threshold_for_drug_selection").unwrap_or(0.10) {
                        // Restart the previously effective drug!
                        individual.cur_use_drug[prev_drug_idx] = true;
                        individual.date_drug_initiated[prev_drug_idx] = time_step as i32;
                        individual.date_drug_initiated_keep[prev_drug_idx] = time_step as i32;
                        individual.ever_taken_drug[prev_drug_idx] = true;
                        
                        // Update drug counter
                        update_drug_counter(individual);
                        
                        // Set drug level
                        let initial_level = get_drug_param(prev_drug_name, "initial_level").unwrap_or(10.0);
                        individual.cur_level_drug[prev_drug_idx] = initial_level;
                        
                        // Reset treatment failure tracking for new treatment
                        individual.bacteria_level_at_drug_start[bacteria_idx] = Some(individual.level[bacteria_idx]);
                        individual.days_on_current_treatment[bacteria_idx] = 0;
                        individual.treatment_failure_assessed[bacteria_idx] = false;
                        
                        return true; // Successfully restarted previously effective drug
                    }
                }
            }
        }
    }
    
    // If we can't restart the previous drug, use standard drug selection with preference bonus
    
    // Build list of available drugs for restart treatment
    let mut drug_scores = Vec::new();
    
    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        // Skip if currently taking this drug
        if individual.cur_use_drug[drug_idx] {
            continue;
        }
        
        // Only avoid drugs that actually failed (not drugs that were stopped due to adherence)
        // We'll identify failed drugs by checking against treatment failure history
        // For now, we don't avoid any recently used drugs since stopped ≠ failed
        
        // Check if drug is available
        let avail = get_drug_availability_time_aware(
            drug_name,
            &individual.region_cur_in.to_string(),
            Some(&individual.region_living.to_string()),
            time_step
        );
        
        if avail < 0.01 {
            continue;
        }
        
        // Calculate drug score
        let mut score = 0.0;
        
        // Base potency score
        let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
        if let Some(potency_param_key) = param_cache.drug_bacteria_potency_keys.get(&(drug_idx, *bacteria_idx_for_cache)) {
            if let Some(potency) = get_global_param(potency_param_key) {
                if potency >= get_global_param("minimal_potency_threshold_for_drug_selection").unwrap_or(0.10) {
                    score += potency;
                }
            }
        }
        
        // Apply clinical preference multipliers
        let bacteria_drug_key = format!("{}_{}_clinical_preference_multiplier", bacteria_name.replace(" ", "_"), drug_name);
        if let Some(preference_multiplier) = get_global_param(&bacteria_drug_key) {
            score *= preference_multiplier;
        }
        
        // BONUS: If this was the previously effective drug, give it preference
        if let Some(prev_drug_idx) = stopped_drug_idx {
            if drug_idx == prev_drug_idx {
                let effectiveness_bonus = get_global_param("previously_effective_drug_bonus").unwrap_or(2.0);
                score *= effectiveness_bonus;
            }
        }
        
        if score > 0.0 {
            drug_scores.push((drug_idx, score));
        }
    }
    
    // Select and start restart treatment
    if !drug_scores.is_empty() {
        let selection_temperature = get_global_param("drug_selection_temperature").unwrap_or(0.5);
        let weights: Vec<f64> = drug_scores.iter()
            .map(|(_, score)| (score / selection_temperature).exp())
            .collect();
        
        let total_weight: f64 = weights.iter().sum();
        if total_weight > 0.0 && total_weight.is_finite() {
            let dist = WeightedIndex::new(&weights).unwrap();
            let chosen_idx = dist.sample(&mut rng);
            let new_drug_idx = drug_scores[chosen_idx].0;
            
            // Start restart treatment
            let new_drug_name = DRUG_SHORT_NAMES[new_drug_idx];
            individual.cur_use_drug[new_drug_idx] = true;
            individual.date_drug_initiated[new_drug_idx] = time_step as i32;
            individual.date_drug_initiated_keep[new_drug_idx] = time_step as i32;
            individual.ever_taken_drug[new_drug_idx] = true;
            
            // Update drug counter
            update_drug_counter(individual);
            
            // Set drug level
            let initial_level = get_drug_param(new_drug_name, "initial_level").unwrap_or(10.0);
            individual.cur_level_drug[new_drug_idx] = initial_level;
            
            // Reset treatment failure tracking for new treatment
            individual.bacteria_level_at_drug_start[bacteria_idx] = Some(individual.level[bacteria_idx]);
            individual.days_on_current_treatment[bacteria_idx] = 0;
            individual.treatment_failure_assessed[bacteria_idx] = false;
            
            return true; // Restart treatment started
        }
    }
    
    false // No restart treatment started
}

/// Pre-computed parameter keys to avoid string allocation during simulation
pub struct ParameterKeyCache {
    // Most frequently used keys - drug/bacteria combinations
    drug_bacteria_potency_keys: HashMap<(usize, usize), String>,
    drug_bacteria_initiation_keys: HashMap<(usize, usize), String>,
    
    // Region-based keys
    region_travel_keys: HashMap<String, String>,
    region_mortality_keys: HashMap<String, String>,
    region_sepsis_keys: HashMap<String, String>,
    region_sepsis_multiplier_keys: HashMap<String, String>,
    region_bacteria_acquisition_keys: HashMap<(String, String), String>,
    region_bacteria_default_keys: HashMap<String, String>,
    region_testing_keys: HashMap<String, String>,
    
    // Other frequently used keys
    syndrome_sepsis_keys: HashMap<String, String>,
    vaccine_age_keys: HashMap<(String, String), String>,
    syndrome_initiation_keys: HashMap<String, String>,
    sex_mortality_keys: HashMap<String, String>,
    hgt_keys: HashMap<(String, String), String>,
    resistance_mechanism_emergence_keys: HashMap<String, String>,
    resistance_mechanism_enhancement_keys: HashMap<String, String>,
}

impl ParameterKeyCache {
    pub fn new() -> Self {
        let mut cache = ParameterKeyCache {
            drug_bacteria_potency_keys: HashMap::new(),
            drug_bacteria_initiation_keys: HashMap::new(),
            region_travel_keys: HashMap::new(),
            region_mortality_keys: HashMap::new(),
            region_sepsis_keys: HashMap::new(),
            region_sepsis_multiplier_keys: HashMap::new(),
            region_bacteria_acquisition_keys: HashMap::new(),
            region_bacteria_default_keys: HashMap::new(),
            region_testing_keys: HashMap::new(),
            syndrome_sepsis_keys: HashMap::new(),
            vaccine_age_keys: HashMap::new(),
            syndrome_initiation_keys: HashMap::new(),
            sex_mortality_keys: HashMap::new(),
            hgt_keys: HashMap::new(),
            resistance_mechanism_emergence_keys: HashMap::new(),
            resistance_mechanism_enhancement_keys: HashMap::new(),
        };
        
        // Pre-compute all drug/bacteria combinations
        for (d_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                cache.drug_bacteria_potency_keys.insert(
                    (d_idx, b_idx),
                    format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug_name, bacteria_name)
                );
                cache.drug_bacteria_initiation_keys.insert(
                    (d_idx, b_idx),
                    format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug_name, bacteria_name)
                );
            }
        }
        
        // Pre-compute region-based keys
        let regions = ["north_america", "europe", "asia", "africa", "south_america", "oceania", "home"];
        for region in &regions {
            cache.region_travel_keys.insert(
                region.to_string(),
                format!("{}_travel_multiplier", region)
            );
            cache.region_mortality_keys.insert(
                region.to_string(),
                format!("log_odds_mortality_region_{}", region)
            );
            cache.region_sepsis_keys.insert(
                region.to_string(),
                format!("sepsis_log_odds_region_{}", region)
            );
            cache.region_sepsis_multiplier_keys.insert(
                region.to_string(),
                format!("{}_sepsis_mortality_multiplier", region)
            );
            cache.region_bacteria_default_keys.insert(
                region.to_string(),
                format!("{}_acquisition_log_odds_default", region)
            );
            cache.region_testing_keys.insert(
                region.to_string(),
                format!("{}_testing_multiplier", region)
            );
            
            // Pre-compute region/bacteria combinations
            for &bacteria in BACTERIA_LIST {
                let bacteria_clean = bacteria.replace(" ", "_");
                cache.region_bacteria_acquisition_keys.insert(
                    (region.to_string(), bacteria_clean.clone()),
                    format!("{}_{}_acquisition_log_odds", region, bacteria_clean)
                );
            }
        }
        
        // Pre-compute sex-based keys
        for sex in &["male", "female"] {
            cache.sex_mortality_keys.insert(
                sex.to_string(),
                format!("log_odds_mortality_sex_{}", sex)
            );
        }
        
        // Pre-compute syndrome keys (numeric syndrome IDs)
        for syndrome_id in 1..=10 {
            cache.syndrome_sepsis_keys.insert(
                syndrome_id.to_string(),
                format!("log_odds_syndrome_{}_sepsis", syndrome_id)
            );
            cache.syndrome_initiation_keys.insert(
                syndrome_id.to_string(),
                format!("syndrome_{}_initiation_multiplier", syndrome_id)
            );
        }
        
        // Pre-compute vaccine/age combinations - only bacterial vaccines
        let bacterial_vaccines = vec!["pneumococcal", "meningococcal", "hib"];
        let age_groups = vec!["0_1", "1_5", "5_18", "18_50", "50_70", "70plus"];
        for vaccine in &bacterial_vaccines {
            for age_group in &age_groups {
                cache.vaccine_age_keys.insert(
                    (vaccine.to_string(), age_group.to_string()),
                    format!("vaccine_{}_daily_prob_age_{}", vaccine, age_group)
                );
            }
        }
        
        // Pre-compute HGT keys for bacteria pairs
        for &donor in BACTERIA_LIST {
            for &recipient in BACTERIA_LIST {
                if donor != recipient {
                    cache.hgt_keys.insert(
                        (donor.to_string(), recipient.to_string()),
                        format!("hgt_prob_{}_to_{}", donor, recipient)
                    );
                }
            }
        }
        
        // Pre-compute resistance mechanism keys using actual enum values
        use crate::simulation::population::ResistanceMechanism;
        for mechanism in ResistanceMechanism::all() {
            let mechanism_str = mechanism.as_str();
            cache.resistance_mechanism_emergence_keys.insert(
                mechanism_str.to_string(),
                format!("resistance_mechanism_{}_emergence_rate", mechanism_str)
            );
            cache.resistance_mechanism_enhancement_keys.insert(
                mechanism_str.to_string(),
                format!("resistance_mechanism_{}_enhancement_multiplier", mechanism_str)
            );
        }
        
        cache
    }
}

/// applies model rules to an individual for one time step.
pub fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    majority_r_positive_values_by_combo: &HashMap<(usize, bool, usize, usize), Vec<f64>>, // <-- update type
    bacteria_indices: &HashMap<&'static str, usize>,
    drug_indices: &HashMap<&'static str, usize>,
    cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>, // New parameter
    param_cache: &ParameterKeyCache, // New parameter cache
) {

    if individual.age < 0 {
        individual.age += 1; // Only advance age by 1 day
        return; // Exit the function if unborn
    }

    if individual.date_of_death.is_some() {
        return; // Exit the function if dead
    }

    let mut rng = rand::thread_rng();

    // --- all these parameter lookups at the top so they're in scope everywhere ---
    let transfer_prob = get_global_param("microbiome_resistance_transfer_probability_per_day").unwrap_or(0.05);
    let drug_base_initiation_rate = get_global_param("drug_base_initiation_rate_per_day").unwrap_or(0.0001);
    let drug_infection_present_multiplier = get_global_param("drug_infection_present_multiplier").unwrap_or(50.0);
    let already_on_drug_initiation_multiplier = get_global_param("already_on_drug_initiation_multiplier").unwrap_or(0.0001);
    let drug_test_identified_multiplier = get_global_param("drug_test_identified_multiplier").unwrap_or(20.0);
    let double_dose_probability = get_global_param("double_dose_probability_if_identified_infection").unwrap_or(0.1);
    let random_drug_cessation_prob = get_global_param("random_drug_cessation_probability").unwrap_or(0.001);

    // update non-infection, bacteria or antibiotic-specific variables
    // need a variable for vulnerability to serious toxicity ?
    individual.age += 1;


    // ---  Update Contact and Exposure Levels ---
    //  update immunodeficiency status based on onset/recovery rates and type
    
    // Get rates for both types
    let temp_onset_rate = get_global_param("temporary_immunosuppression_onset_rate_per_day").unwrap_or(0.0002);
    let temp_recovery_rate = get_global_param("temporary_immunosuppression_recovery_rate_per_day").unwrap_or(0.01);
    let chronic_onset_rate = get_global_param("chronic_immunosuppression_onset_rate_per_day").unwrap_or(0.0001);
    let chronic_recovery_rate = get_global_param("chronic_immunosuppression_recovery_rate_per_day").unwrap_or(0.0005);
    
    // Get age-based probability for chronic vs temporary assignment
    let chronic_probability = if individual.age <= 365 { // 0-1 years (365 days)
        get_global_param("chronic_immunodeficiency_probability_age_0_1").unwrap_or(0.3)
    } else if individual.age <= 6570 { // 1-18 years (18*365 days)
        get_global_param("chronic_immunodeficiency_probability_age_1_18").unwrap_or(0.2)
    } else if individual.age <= 23725 { // 18-65 years (65*365 days)
        get_global_param("chronic_immunodeficiency_probability_age_18_65").unwrap_or(0.4)
    } else { // 65+ years
        get_global_param("chronic_immunodeficiency_probability_age_65_plus").unwrap_or(0.6)
    };

    match individual.immunodeficiency_type {
        Some(ImmunodeficiencyType::Temporary) => {
            // Currently has temporary immunodeficiency, check for recovery
            if rng.gen_bool(temp_recovery_rate) {
                individual.immunodeficiency_type = None;
            }
        },
        Some(ImmunodeficiencyType::Chronic) => {
            // Currently has chronic immunodeficiency, check for recovery  
            if rng.gen_bool(chronic_recovery_rate) {
                individual.immunodeficiency_type = None;
            }
        },
        None => {
            // Not currently immunodeficient, check for onset
            let total_onset_rate = temp_onset_rate + chronic_onset_rate;
            if rng.gen_bool(total_onset_rate) {
                // Determine type based on age
                if rng.gen_bool(chronic_probability) {
                    individual.immunodeficiency_type = Some(ImmunodeficiencyType::Chronic);
                } else {
                    individual.immunodeficiency_type = Some(ImmunodeficiencyType::Temporary);
                }
            }
        }
    }


    // current toxicity
    individual.current_toxicity = (individual.current_toxicity + rng.gen_range(-0.5..=0.5)).max(0.0);


    // Get parameters from config.rs once per individual for this time step
    let baseline_rate = get_global_param("hospitalization_baseline_rate_per_day")
        .expect("Missing hospitalization_baseline_rate_per_day in config");
    let age_multiplier_hosp = get_global_param("hospitalization_age_multiplier_per_day")
        .expect("Missing hospitalization_age_multiplier_per_day in config");
    let recovery_rate = get_global_param("hospitalization_recovery_rate_per_day")
        .expect("Missing hospitalization_recovery_rate_per_day in config");
    let max_days_in_hospital = get_global_param("hospitalization_max_days")
        .expect("Missing hospitalization_max_days in config");
    let sepsis_admission_multiplier = get_global_param("hospitalization_sepsis_admission_multiplier")
        .expect("Missing hospitalization_sepsis_admission_multiplier in config");
    let prevent_discharge_with_sepsis = get_global_param("hospitalization_prevent_discharge_with_sepsis")
        .expect("Missing hospitalization_prevent_discharge_with_sepsis in config");

    // Check if individual has any active sepsis
    let has_sepsis = individual.sepsis.iter().any(|&s| s);

    // Potentially get hospitalized (if not currently hospitalized)
    if !individual.hospital_status.is_hospitalized() { 
        let mut prob_hospitalization_today = baseline_rate + (individual.age as f64 * age_multiplier_hosp);
        
        // Strong sepsis admission effect - sepsis patients are very likely to be hospitalized
        if has_sepsis {
            prob_hospitalization_today *= sepsis_admission_multiplier;
        }

        if rng.gen::<f64>() < prob_hospitalization_today {
            individual.hospital_status = HospitalStatus::InHospital; 
            individual.days_hospitalized = 0; // Initialize days hospitalized
        }
    } else { // If already hospitalized, consider recovery or max days limit
        individual.days_hospitalized += 1; // Increment days hospitalized

        // Determine if discharge is allowed
        let can_discharge = if prevent_discharge_with_sepsis > 0.5 {
            !has_sepsis // Cannot discharge if patient has sepsis
        } else {
            true // Can always discharge (old behavior)
        };

        // Potentially recover from hospitalization (only if discharge is allowed)
        if can_discharge && rng.gen::<f64>() < recovery_rate {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
            // println!("individual {} recovered from hospitalization.", individual.id);             
        }
        // discharge after max_days_in_hospital (only if discharge is allowed)
        else if can_discharge && individual.days_hospitalized >= max_days_in_hospital as u32 {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
         }
    }
    // --- end hospitalization Rules ---



    // ---  region travel ---
    let base_travel_prob = get_global_param("travel_probability_per_day")
        .expect("Missing travel_probability_per_day in config");
    
    // Apply region-specific travel multiplier based on individual's home region
    let region_name_for_param = individual.region_living.to_string().to_lowercase().replace(" ", "_");
    let region_travel_multiplier_key = &param_cache.region_travel_keys[&region_name_for_param];
    let region_travel_multiplier = get_global_param(region_travel_multiplier_key).unwrap_or(1.0);
    let travel_prob = base_travel_prob * region_travel_multiplier;
    
    const VISIT_LENGTH_DAYS: u32 = 30; // Fixed visit length

    // Check if the individual is currently in their home region
    if let Region::Home = individual.region_cur_in {
        // If not hospitalized, consider initiating travel
        if !individual.hospital_status.is_hospitalized() && rng.gen::<f64>() < travel_prob {
            // Initiate travel: select a random new region different from their living region
            let mut new_region: Region;
            loop {
                // Select destination based on economic development level (main determinant of travel patterns)
                // Higher-income regions have more global travel; lower-income regions travel more regionally
                let destinations = match individual.region_living {
                    Region::NorthAmerica | Region::Europe | Region::Oceania => {
                        // High-income regions: global travel with preference for other developed regions
                        vec![
                            (Region::Europe, 0.35),       // Strong developed-to-developed flow
                            (Region::Asia, 0.25),         // Major business/tourism destination
                            (Region::NorthAmerica, 0.15), // Cross-Atlantic travel
                            (Region::Oceania, 0.10),      // Tourism/business
                            (Region::SouthAmerica, 0.10), // Tourism/business
                            (Region::Africa, 0.05),       // Lower but still significant
                        ]
                    },
                    Region::Asia => {
                        // Mixed income: regional preference with some global reach
                        vec![
                            (Region::Asia, 0.40),         // Strong regional travel
                            (Region::Europe, 0.20),       // Business/education
                            (Region::NorthAmerica, 0.15), // Business/education
                            (Region::Oceania, 0.10),      // Regional proximity
                            (Region::Africa, 0.08),       // Growing connections
                            (Region::SouthAmerica, 0.07), // Limited
                        ]
                    },
                    Region::SouthAmerica => {
                        // Middle income: regional focus with some international travel
                        vec![
                            (Region::SouthAmerica, 0.40), // Strong regional travel
                            (Region::NorthAmerica, 0.25), // Geographic proximity
                            (Region::Europe, 0.15),       // Historical ties
                            (Region::Asia, 0.10),         // Growing connections
                            (Region::Africa, 0.05),       // Limited
                            (Region::Oceania, 0.05),      // Limited
                        ]
                    },
                    Region::Africa => {
                        // Lower income: primarily regional travel
                        vec![
                            (Region::Africa, 0.50),       // Strong regional travel
                            (Region::Europe, 0.20),       // Historical/economic ties
                            (Region::Asia, 0.15),         // Growing connections
                            (Region::NorthAmerica, 0.08), // Limited
                            (Region::SouthAmerica, 0.04), // Very limited
                            (Region::Oceania, 0.03),      // Very limited
                        ]
                    },
                    Region::Home => {
                        // Should not reach here, but default to global uniform if it does
                        vec![
                            (Region::Asia, 0.167), (Region::Africa, 0.167), (Region::Europe, 0.166),
                            (Region::NorthAmerica, 0.167), (Region::SouthAmerica, 0.166), (Region::Oceania, 0.167),
                        ]
                    },
                };
                
                // Sample from the economic-based destination distribution
                let rand_val = rng.gen::<f64>();
                let mut cumulative_prob = 0.0;
                new_region = Region::Asia; // Default fallback
                
                for (region, prob) in destinations {
                    cumulative_prob += prob;
                    if rand_val < cumulative_prob {
                        new_region = region;
                        break;
                    }
                }
                
                // Ensure the individual doesn't 'travel' to their own living region
                if new_region != individual.region_living {
                    break; // Found a suitable new region to visit
                }
            }
            individual.region_cur_in = new_region;
            individual.days_visiting = 1; // Start the visit counter at 1
        }
    } else {
        // Individual is currently visiting another region
        individual.days_visiting += 1; // Increment the visit duration

        // Check if the visit duration has been reached
        if individual.days_visiting >= VISIT_LENGTH_DAYS {
            // End of visit, rto home region
            individual.region_cur_in = Region::Home; // Set current region back to Home
            individual.days_visiting = 0; // Reset visit counter
            // println!("individual {} (Age: {}) returned home from a trip.",
            //     time_step, individual.id, individual.age);
        }
    }
    // --- end region travel updates ---

    // ---  sepsis risk  ---
    for &bacteria in BACTERIA_LIST.iter() {
        let b_idx = *bacteria_indices.get(bacteria).unwrap();
        let current_level = individual.level[b_idx];
        
        if current_level > 0.0 {
            // Only calculate sepsis onset risk if not already septic from this bacteria
            if !individual.sepsis[b_idx] {
                let last_infected_day = individual.date_last_infected[b_idx];
                let duration_of_infection = (time_step as i32 - last_infected_day).max(0); // Ensure non-negative duration

                // Logistic regression model for sepsis risk
                // Retrieve logistic parameters, falling back to global defaults
                let sepsis_baseline_odds = get_bacteria_param(bacteria, "sepsis_baseline_odds")
                    .unwrap_or_else(|| get_global_param("sepsis_baseline_odds").expect("Missing sepsis_baseline_odds"));
                let log_odds_infection_level = get_bacteria_param(bacteria, "log_odds_sepsis_infection_level")
                    .unwrap_or_else(|| get_global_param("log_odds_sepsis_infection_level").expect("Missing log_odds_sepsis_infection_level"));
                let log_odds_infection_duration = get_bacteria_param(bacteria, "log_odds_sepsis_infection_duration")
                    .unwrap_or_else(|| get_global_param("log_odds_sepsis_infection_duration").expect("Missing log_odds_sepsis_infection_duration"));

                // ENHANCED BACTERIA SEPSIS RISK CALCULATION
                // Combines: 1) Enhanced bacteria-specific risk, 2) Age-dependent interactions, 3) Clinical risk categories
                let bacteria_sepsis_risk = get_bacteria_sepsis_risk_multiplier(bacteria);
                let age_bacteria_sepsis_risk = get_age_dependent_bacteria_sepsis_risk_multiplier(bacteria, individual.age as u32);
                
                // Combined bacteria risk multiplier (bacteria-specific × age-dependent interaction)
                let combined_bacteria_risk = bacteria_sepsis_risk * age_bacteria_sepsis_risk;
                
                // Map combined risk to log odds categories for logistic regression
                let bacteria_log_odds = if combined_bacteria_risk >= 3.0 {
                    // Very high combined risk (e.g., MRSA in elderly, GBS in neonates)
                    get_global_param("log_odds_bacteria_with_high_sepsis_risk").expect("Missing log_odds_bacteria_with_high_sepsis_risk") * 1.5
                } else if combined_bacteria_risk >= 1.8 {
                    // High combined risk 
                    get_global_param("log_odds_bacteria_with_high_sepsis_risk").expect("Missing log_odds_bacteria_with_high_sepsis_risk")
                } else if combined_bacteria_risk >= 0.7 && combined_bacteria_risk <= 1.3 {
                    // Medium combined risk (reference category)
                    get_global_param("log_odds_bacteria_with_medium_sepsis_risk").expect("Missing log_odds_bacteria_with_medium_sepsis_risk")
                } else if combined_bacteria_risk >= 0.3 {
                    // Low combined risk
                    get_global_param("log_odds_bacteria_with_low_sepsis_risk").expect("Missing log_odds_bacteria_with_low_sepsis_risk")
                } else {
                    // Very low combined risk (e.g., Chlamydia, localized infections)
                    get_global_param("log_odds_bacteria_with_low_sepsis_risk").expect("Missing log_odds_bacteria_with_low_sepsis_risk") * 0.5
                };

                // Add syndrome-specific sepsis risk (infection site effect)
                // This allows the same bacteria to have different sepsis risks depending on infection site
                // e.g., E. coli UTI vs E. coli bacteremia have very different sepsis risks
                let syndrome_log_odds = if individual.infectious_syndrome[b_idx] != 0 {
                    let syndrome_id = individual.infectious_syndrome[b_idx].to_string();
                    if let Some(param_name) = param_cache.syndrome_sepsis_keys.get(&syndrome_id) {
                        get_global_param(param_name).unwrap_or(0.0)
                    } else {
                        println!("Warning: Missing syndrome cache key for syndrome ID: {}", syndrome_id);
                        0.0 // Default to no effect if missing
                    }
                } else {
                    0.0 // No syndrome specified, no effect
                };
                
                // Add regional sepsis risk factors (healthcare access, population density, resources)
                let region_log_odds = match individual.region_living {
                    Region::Africa | Region::Asia => get_global_param("log_odds_sepsis_region_b").unwrap_or(0.2),  // Lower resource regions
                    Region::NorthAmerica | Region::Europe | Region::Oceania => get_global_param("log_odds_sepsis_region_a").unwrap_or(-0.3), // Higher resource regions
                    Region::SouthAmerica => get_global_param("log_odds_sepsis_region_b").unwrap_or(0.1),  // Mixed resource region
                    Region::Home => 0.0, // Neutral/no effect for home region
                };
                
                // COMPREHENSIVE SEPSIS RISK CALCULATION
                // Integrates: bacteria risk, age interactions, syndrome site, regional factors
                let log_odds_sepsis = sepsis_baseline_odds
                                    + (current_level * log_odds_infection_level)
                                    + (duration_of_infection as f64 * log_odds_infection_duration)
                                    + bacteria_log_odds
                                    + syndrome_log_odds
                                    + region_log_odds;

                // Convert log odds to probability using logistic function
                let prob_sepsis_today = 1.0 / (1.0 + (-log_odds_sepsis).exp());

                if rng.gen::<f64>() < prob_sepsis_today {
                    // Set sepsis status to true for this bacteria and record onset day
                    individual.sepsis[b_idx] = true;
                    individual.sepsis_onset_day[b_idx] = time_step as i32;
                }
            }
            // Note: Recovery logic will be applied later, after death risk is calculated
        } else {
            // If infection has cleared, sepsis should also clear
            if individual.sepsis[b_idx] {
                individual.sepsis[b_idx] = false;
            }
        }
    }
    // --- end sepsis updates ---




    // Update vaccination status dynamically based on age-appropriate schedules
    // Only bacterial vaccines with historical availability checking
    // Vaccines: pneumococcal (1977+), meningococcal (1981+), hib (1985+)
    // Age groups: 0-1, 1-5, 5-18, 18-50, 50-70, 70+
    let age_years = individual.age as f64 / 365.0;
    let age_group = if age_years < 1.0 {
        "0_1"
    } else if age_years < 5.0 {
        "1_5"
    } else if age_years < 18.0 {
        "5_18"
    } else if age_years < 50.0 {
        "18_50"
    } else if age_years < 70.0 {
        "50_70"
    } else {
        "70plus"
    };

    // Calculate simulation year (assuming time_step 0 = year 1950, one step per day)
    let simulation_year = 1950.0 + (time_step as f64 / 365.0);

    let bacterial_vaccines = vec!["pneumococcal", "meningococcal", "hib"];
    for (b_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
        // For each bacterial vaccine, check if this bacteria is targeted by the vaccine
        for vaccine in &bacterial_vaccines {
            // Check vaccine availability date first
            let availability_year = get_global_param(&format!("vaccine_{}_availability_year", vaccine)).unwrap_or(2100.0); // Default to future if not set
            if simulation_year < availability_year {
                continue; // Vaccine not yet available
            }

            // Correct bacteria name matching (fixing underscore vs space issues)
            let targets_bacteria = match (*vaccine, *bacteria) {
                ("pneumococcal", "streptococcus pneumoniae") => true,
                ("meningococcal", "neisseria_meningitidis") => true,  // Fixed: using underscore version
                ("hib", "haemophilus influenzae") => true,
                ("pertussis", "bordetella pertussis") => true,  // DTaP/Tdap vaccines
                _ => false,
            };
            if targets_bacteria && !individual.vaccination_status[b_idx] {
                let param_key = &param_cache.vaccine_age_keys[&(vaccine.to_string(), age_group.to_string())];
                let daily_prob = get_global_param(param_key).unwrap_or(0.0);
                if rng.gen::<f64>() < daily_prob {
                    individual.vaccination_status[b_idx] = true;
                }
            }
        }
    }

    // --- drug updates---
    let has_any_infection = individual.level.iter().any(|&level| level > 0.0);
    let initial_on_any_antibiotic = individual.cur_use_drug.iter().any(|&identified| identified);
    let has_any_identified_infection = individual.test_identified_infection.iter().any(|&identified| identified);

    // --- count number of drugs currently being used ---
    let num_drugs_currently_used = individual.cur_use_drug.iter().filter(|&&on| on).count();

    let mut syndrome_administration_multiplier: f64 = 1.0;
    for &syndrome_id in individual.infectious_syndrome.iter() {
        if syndrome_id != 0 {
            if let Some(param_name) = param_cache.syndrome_initiation_keys.get(&syndrome_id.to_string()) {
                if let Some(multiplier) = get_global_param(param_name) {
                    syndrome_administration_multiplier = syndrome_administration_multiplier.max(multiplier);
                }
            }
        }
    }

let drugs_initiated_this_time_step: usize = 0;

    // --- drug stopping ---
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        if individual.cur_use_drug[drug_idx] {
            let mut relevant_infection_active_for_this_drug = false;
            let mut primary_bacteria_idx: Option<usize> = None;
            let mut highest_bacteria_level = 0.0;
            
            // Find the most significant bacteria infection relevant to this drug
            for b_idx in 0..BACTERIA_LIST.len() {
                if individual.level[b_idx] > 0.0001 {
                    // Use potency_when_no_r to determine if drug is relevant for this bacteria
                    let potency_param_key = &param_cache.drug_bacteria_potency_keys[&(drug_idx, b_idx)];
                    let drug_potency = get_global_param(potency_param_key).unwrap_or(0.0);
                    if drug_potency > 0.0 {
                        relevant_infection_active_for_this_drug = true;
                        // Track the bacteria with highest level (most significant infection)
                        if individual.level[b_idx] > highest_bacteria_level {
                            highest_bacteria_level = individual.level[b_idx];
                            primary_bacteria_idx = Some(b_idx);
                        }
                    }
                }
            }
            
            let mut stop_drug = false;
            
            if !relevant_infection_active_for_this_drug {
                // No relevant infection - use higher cessation rate
                let random_cessation_if_no_infection = get_global_param("random_drug_cessation_probability_if_no_active_infection").unwrap_or(0.25);
                if rng.gen_bool(random_cessation_if_no_infection) {
                    stop_drug = true;
                }
            } else {
                // Calculate bacteria-specific and region-specific cessation probability
                let base_cessation_prob = if let Some(bacteria_idx) = primary_bacteria_idx {
                    let bacteria_name = BACTERIA_LIST[bacteria_idx];
                    let bacteria_cessation_key = format!("{}_drug_cessation_probability", bacteria_name.to_lowercase().replace(' ', "_"));
                    get_global_param(&bacteria_cessation_key).unwrap_or(random_drug_cessation_prob)
                } else {
                    random_drug_cessation_prob
                };
                
                // Apply regional multiplier based on individual's current region
                let region_multiplier = match individual.region_cur_in {
                    Region::NorthAmerica => get_global_param("north_america_cessation_multiplier").unwrap_or(1.0),
                    Region::SouthAmerica => get_global_param("south_america_cessation_multiplier").unwrap_or(1.0),
                    Region::Africa => get_global_param("africa_cessation_multiplier").unwrap_or(1.0),
                    Region::Asia => get_global_param("asia_cessation_multiplier").unwrap_or(1.0),
                    Region::Europe => get_global_param("europe_cessation_multiplier").unwrap_or(1.0),
                    Region::Oceania => get_global_param("oceania_cessation_multiplier").unwrap_or(1.0),
                    Region::Home => {
                        // Use home region multiplier (region_living is their home region)
                        match individual.region_living {
                            Region::NorthAmerica => get_global_param("north_america_cessation_multiplier").unwrap_or(1.0),
                            Region::SouthAmerica => get_global_param("south_america_cessation_multiplier").unwrap_or(1.0),
                            Region::Africa => get_global_param("africa_cessation_multiplier").unwrap_or(1.0),
                            Region::Asia => get_global_param("asia_cessation_multiplier").unwrap_or(1.0),
                            Region::Europe => get_global_param("europe_cessation_multiplier").unwrap_or(1.0),
                            Region::Oceania => get_global_param("oceania_cessation_multiplier").unwrap_or(1.0),
                            Region::Home => 1.0, // Fallback if nested Home
                        }
                    }
                };
                
                let final_cessation_prob = (base_cessation_prob * region_multiplier).min(0.99); // Cap at 99%
                
                if rng.gen_bool(final_cessation_prob) {
                    stop_drug = true;
                }
            }
            if individual.date_drug_initiated[drug_idx] == (time_step as i32) - 1 {
                stop_drug = false;
            }
            if stop_drug {
                individual.cur_use_drug[drug_idx] = false;
                individual.date_drug_initiated[drug_idx] = i32::MIN;
                
                // Update drug counter
                update_drug_counter(individual);
                
                // Check if stopping while infection persists (restart window logic)
                for bacteria_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[bacteria_idx] > 0.1 && // Still infected (threshold for meaningful infection)
                       individual.bacteria_level_at_drug_start[bacteria_idx].is_some() {
                        
                        // Record cessation for restart window tracking
                        individual.drug_stopped_with_infection_day[bacteria_idx] = Some(time_step as i32);
                        individual.bacteria_level_at_drug_cessation[bacteria_idx] = Some(individual.level[bacteria_idx]);
                        individual.stopped_drug_index[bacteria_idx] = Some(drug_idx); // Track which drug was stopped
                        individual.restart_window_assessed[bacteria_idx] = false;
                    }
                    
                    // Reset treatment failure tracking when drug is stopped naturally
                    if individual.bacteria_level_at_drug_start[bacteria_idx].is_some() {
                        individual.bacteria_level_at_drug_start[bacteria_idx] = None;
                        individual.days_on_current_treatment[bacteria_idx] = -1;
                        individual.treatment_failure_assessed[bacteria_idx] = false;
                    }
                }
            }
        }
    }

    // apply decay if stopped, or set to initial level if continued/re-initiated.
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];
        let drug_initial_level = get_drug_param(drug_name, "initial_level").unwrap_or(10.0);
        if individual.cur_use_drug[drug_idx] {
            individual.cur_level_drug[drug_idx] = drug_initial_level;
        } else {
            // Use exponential decay based on drug-specific half-life
            let half_life_days = get_drug_param(drug_name, "half_life_days").unwrap_or(0.25); // Default ~6 hours
            let decay_constant = (2.0_f64).ln() / half_life_days; // k = ln(2) / t_half
            let decay_factor = (-decay_constant).exp(); // e^(-k*t) where t=1 day
            let new_level = individual.cur_level_drug[drug_idx] * decay_factor;
            // Set levels below 0.001 (0.1% of standard dose) to exactly zero to avoid floating point artifacts
            individual.cur_level_drug[drug_idx] = if new_level < 0.001 { 0.0 } else { new_level };
        }
    }

    // --- drug initiation (two-stage process) ---
    // Stage 1: Decide whether to start any antibiotic
let available_drugs: Vec<usize> = DRUG_SHORT_NAMES.iter().enumerate()
    .filter(|(_, &name)| {
            let avail = get_drug_availability_time_aware(
                name,
                &individual.region_cur_in.to_string(),
                Some(&individual.region_living.to_string()),
                time_step
            );
            let intro_ok = match get_drug_introduction_time_step(name) {
                Some(intro_step) => time_step >= intro_step,
                None => true,
            };
            avail >= 0.01 && intro_ok
        })
        .map(|(idx, _)| idx)
        .collect();
    let available_drugs_count = available_drugs.len();
    let min_available_drugs = 5; // Adjustable threshold
    let scaling_factor = if available_drugs_count < min_available_drugs && available_drugs_count > 0 {
        (min_available_drugs as f64) / (available_drugs_count as f64)
    } else { 1.0 };

    // Restriction: if already using three or more drugs, cannot start another (allow up to 3 drugs for severe infections)
    if num_drugs_currently_used + drugs_initiated_this_time_step < 3 && available_drugs_count > 0 {
        // Stage 1: Calculate probability to start any antibiotic
        let mut start_any_antibiotic_prob = drug_base_initiation_rate * scaling_factor;
        let infection_acquired_this_step = individual.date_last_infected.iter().any(|&d| d == time_step as i32);
        if has_any_infection && !infection_acquired_this_step {
            start_any_antibiotic_prob *= drug_infection_present_multiplier;
        }
        if has_any_identified_infection { start_any_antibiotic_prob *= drug_test_identified_multiplier; }
        if initial_on_any_antibiotic || drugs_initiated_this_time_step > 0 {
            start_any_antibiotic_prob *= already_on_drug_initiation_multiplier;
        }
        // Immunocompromised patients more likely to receive prophylactic antibiotics
        if individual.immunodeficiency_type.is_some() {
            let prophylactic_multiplier = get_global_param("immunodeficiency_prophylactic_drug_multiplier").unwrap_or(8.0);
            start_any_antibiotic_prob *= prophylactic_multiplier;
        }
        start_any_antibiotic_prob *= syndrome_administration_multiplier;
        start_any_antibiotic_prob = start_any_antibiotic_prob.clamp(0.0, 1.0);

        if rng.gen_bool(start_any_antibiotic_prob) {
            // Identify primary bacteria for drug score tracking (highest level among infected bacteria)
            let mut primary_bacteria_idx = -1i32;
            let mut highest_bacteria_level = 0.0;
            for b_idx in 0..BACTERIA_LIST.len() {
                if individual.level[b_idx] > 0.001 && individual.level[b_idx] > highest_bacteria_level {
                    highest_bacteria_level = individual.level[b_idx];
                    primary_bacteria_idx = b_idx as i32;
                }
            }
            
            // Store primary bacteria index for this drug selection event
            individual.bacteria_on_selection_day = primary_bacteria_idx;
            
            // Stage 2: Choose the most appropriate drug using weighted probabilistic selection
            // Score each available drug and collect scores for probabilistic selection
            let mut drug_scores: Vec<(usize, f64)> = Vec::new();
            for &drug_idx in &available_drugs {
                let drug_name = DRUG_SHORT_NAMES[drug_idx];
                // Restriction: do not start drug if resistance test has been performed and resistance detected for any bacteria
                let mut resistance_detected = false;
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.test_for_resistance[b_idx] && individual.resistances[b_idx][drug_idx].test_r > 0.0 {
                        resistance_detected = true;
                        break;
                    }
                }
                if resistance_detected { continue; }

                // Score drug based on spectrum, activity, and clinical scenario
                let mut score = 1.0;
                
                // INTRINSIC ACTIVITY GATE: Block drugs with no meaningful activity against current infections
                let minimal_potency_threshold = get_global_param("minimal_potency_threshold_for_drug_selection").unwrap_or(0.10);
                let mut has_meaningful_activity = false;
                let mut max_potency_against_infections: f64 = 0.0;
                
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[b_idx] > 0.001 {
                        let potency_param_key = &param_cache.drug_bacteria_potency_keys[&(drug_idx, b_idx)];
                        let potency = get_global_param(potency_param_key).unwrap_or(0.0);
                        max_potency_against_infections = max_potency_against_infections.max(potency);
                        if potency >= minimal_potency_threshold {
                            has_meaningful_activity = true;
                        }
                    }
                }
                
                // Block drugs with insufficient activity against any current infection
                if !has_meaningful_activity && has_any_infection {
                    continue; // Skip this drug entirely - no meaningful activity
                }
                
                // PATHOGEN-SPECIFIC CLINICAL GUIDELINES: Boost appropriate drugs, block inappropriate ones
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[b_idx] > 0.001 {
                        let bacteria_name = BACTERIA_LIST[b_idx];
                        match (bacteria_name, drug_name) {
                            // Pseudomonas aeruginosa - strict anti-pseudomonal agents only (MUCH stronger multipliers)
                            ("Pseudomonas aeruginosa", "piperacillin_tazobactam") => score *= 25.0,
                            ("Pseudomonas aeruginosa", "ceftazidime") => score *= 20.0,
                            ("Pseudomonas aeruginosa", "cefepime") => score *= 22.0,
                            ("Pseudomonas aeruginosa", "meropenem") => score *= 25.0,
                            ("Pseudomonas aeruginosa", "imipenem") => score *= 20.0,
                            ("Pseudomonas aeruginosa", "ciprofloxacin") => score *= 18.0,
                            ("Pseudomonas aeruginosa", "tobramycin") => score *= 15.0,
                            ("Pseudomonas aeruginosa", "colistin") => score *= 12.0,
                            ("Pseudomonas aeruginosa", "penicillin" | "ampicillin" | "amoxicillin" | "cephalexin" | "ceftriaxone" | "vancomycin") => {
                                score = 0.0; // Completely block - no intrinsic activity
                                break;
                            },
                            
                            // Staphylococcus aureus - DRAMATICALLY strengthen MSSA vs MRSA logic
                            ("Staphylococcus aureus", "penicillin") => {
                                // Early periods: penicillin should dominate (MSSA era)
                                if time_step < 7300 { // First ~20 years
                                    score *= 50.0; // MASSIVE boost for MSSA
                                } else {
                                    score *= 2.0; // Minimal in MRSA era
                                }
                            },
                            ("Staphylococcus aureus", "methicillin" | "flucloxacillin") => {
                                if time_step < 10950 { // First ~30 years
                                    score *= 40.0; // Major boost before MRSA dominance
                                } else {
                                    score *= 3.0; // Reduced in MRSA era
                                }
                            },
                            ("Staphylococcus aureus", "vancomycin") => {
                                if time_step < 7300 { // Early years
                                    score *= 2.0; // Minimal early use
                                } else { // MRSA era
                                    score *= 35.0; // MASSIVE boost for MRSA
                                }
                            },
                            ("Staphylococcus aureus", "linezolid" | "daptomycin") => {
                                if time_step >= 10950 { // Late period only
                                    score *= 25.0; // Strong alternatives to vancomycin
                                } else {
                                    score *= 0.5; // Minimal early use
                                }
                            },
                            ("Staphylococcus aureus", "clindamycin") => score *= 8.0,
                            
                            // E. coli - MASSIVELY strengthen first-line agents
                            ("Escherichia coli", "ciprofloxacin") => score *= 35.0, // Major UTI drug
                            ("Escherichia coli", "nitrofurantoin") => score *= 30.0, // Cystitis first-line
                            ("Escherichia coli", "trimethoprim_sulfamethoxazole") => score *= 25.0,
                            ("Escherichia coli", "ceftriaxone") => score *= 20.0, // Serious infections
                            ("Escherichia coli", "ampicillin") => {
                                if time_step < 7300 { // Early susceptible era
                                    score *= 25.0;
                                } else {
                                    score *= 3.0; // Resistance emerged
                                }
                            },
                            ("Escherichia coli", "meropenem" | "imipenem") => {
                                // Carbapenems should be rare for E. coli except ESBL era
                                if time_step >= 14600 { // Later periods for ESBL
                                    score *= 8.0;
                                } else {
                                    score *= 0.2; // Minimal early use
                                }
                            },
                            
                            // Klebsiella pneumoniae - strengthen appropriate agents
                            ("Klebsiella pneumoniae", "ceftriaxone") => {
                                if time_step < 10950 { // Before ESBL dominance
                                    score *= 25.0;
                                } else {
                                    score *= 8.0;
                                }
                            },
                            ("Klebsiella pneumoniae", "meropenem" | "imipenem") => {
                                if time_step >= 10950 { // ESBL era
                                    score *= 30.0;
                                } else {
                                    score *= 3.0;
                                }
                            },
                            ("Klebsiella pneumoniae", "ciprofloxacin") => score *= 15.0,
                            ("Klebsiella pneumoniae", "piperacillin_tazobactam") => score *= 18.0,
                            
                            // Enterococcus faecalis - strengthen appropriate agents
                            ("Enterococcus faecalis", "ampicillin") => score *= 40.0, // First-line
                            ("Enterococcus faecalis", "vancomycin") => {
                                if time_step >= 10950 { // VRE era
                                    score *= 30.0;
                                } else {
                                    score *= 8.0;
                                }
                            },
                            ("Enterococcus faecalis", "linezolid") => {
                                if time_step >= 14600 { // Late VRE era
                                    score *= 25.0;
                                } else {
                                    score *= 2.0;
                                }
                            },
                            
                            // Enterococcus faecium - more resistant, different pattern
                            ("Enterococcus faecium", "ampicillin") => score *= 5.0, // Less effective than faecalis
                            ("Enterococcus faecium", "vancomycin") => {
                                if time_step >= 10950 {
                                    score *= 35.0;
                                } else {
                                    score *= 15.0;
                                }
                            },
                            ("Enterococcus faecium", "linezolid") => {
                                if time_step >= 14600 {
                                    score *= 30.0;
                                } else {
                                    score *= 3.0;
                                }
                            },
                            ("Enterococcus faecium", "quinupristin_dalfopristin") => {
                                if time_step >= 16425 { // Very late introduction
                                    score *= 20.0;
                                }
                            },
                            
                            // Acinetobacter baumannii - highly resistant pathogen
                            ("Acinetobacter baumannii", "meropenem" | "imipenem") => {
                                if time_step < 18250 { // Before extensive carbapenem resistance
                                    score *= 40.0;
                                } else {
                                    score *= 15.0;
                                }
                            },
                            ("Acinetobacter baumannii", "colistin") => {
                                if time_step >= 14600 { // Later periods for MDR
                                    score *= 35.0;
                                } else {
                                    score *= 8.0;
                                }
                            },
                            ("Acinetobacter baumannii", "ampicillin_sulbactam") => score *= 25.0, // Intrinsic activity
                            
                            _ => {} // No specific guideline
                        }
                    }
                }
                
                // If drug was blocked by pathogen-specific guidelines, skip it
                if score <= 0.0 {
                    continue;
                }
                
                // CLINICAL CONCENTRATION FORCE: Heavily penalize drugs that aren't first/second-line
                // This creates realistic clinical concentration patterns
                let mut is_first_or_second_line = false;
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[b_idx] > 0.001 {
                        let bacteria_name = BACTERIA_LIST[b_idx];
                        let first_second_line_drugs = match bacteria_name {
                            "Pseudomonas aeruginosa" => vec!["piperacillin_tazobactam", "meropenem", "ceftazidime", "cefepime", "ciprofloxacin", "tobramycin"],
                            "Staphylococcus aureus" => vec!["penicillin", "methicillin", "flucloxacillin", "vancomycin", "linezolid", "daptomycin", "clindamycin"],
                            "Escherichia coli" => vec!["ciprofloxacin", "nitrofurantoin", "trimethoprim_sulfamethoxazole", "ceftriaxone", "ampicillin", "cefuroxime"],
                            "Klebsiella pneumoniae" => vec!["ceftriaxone", "meropenem", "imipenem", "ciprofloxacin", "piperacillin_tazobactam", "ertapenem"],
                            "Enterococcus faecalis" => vec!["ampicillin", "vancomycin", "linezolid", "daptomycin"],
                            "Enterococcus faecium" => vec!["vancomycin", "linezolid", "daptomycin", "quinupristin_dalfopristin"],
                            "Acinetobacter baumannii" => vec!["meropenem", "imipenem", "colistin", "ampicillin_sulbactam", "tigecycline"],
                            _ => vec![], // For other bacteria, no specific restriction
                        };
                        
                        if first_second_line_drugs.contains(&drug_name) {
                            is_first_or_second_line = true;
                            break;
                        }
                    }
                }
                
                // Heavily penalize drugs that aren't first/second-line for current infections
                if has_any_infection && !is_first_or_second_line {
                    score *= 0.05; // 95% penalty for non-standard choices
                }
                
                // POTENCY-BASED POSITIVE REINFORCEMENT: Reward high-potency drugs (MUCH STRONGER)
                if max_potency_against_infections >= 0.5 {
                    score *= 15.0; // Very high potency - MASSIVE boost
                } else if max_potency_against_infections >= 0.3 {
                    score *= 10.0; // High potency - major boost
                } else if max_potency_against_infections >= 0.15 {
                    score *= 6.0; // Moderate potency - significant boost
                } else if max_potency_against_infections >= minimal_potency_threshold {
                    score *= 2.0; // Minimal acceptable potency
                }
                
                let mut max_bacteria_specific_multiplier: f64 = 1.0;
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[b_idx] > 0.001 {
                        let param_key = &param_cache.drug_bacteria_initiation_keys[&(drug_idx, b_idx)];
                        if let Some(specific_multiplier) = get_global_param(param_key) {
                            max_bacteria_specific_multiplier = max_bacteria_specific_multiplier.max(specific_multiplier);
                        }
                    }
                }
                score *= max_bacteria_specific_multiplier;

                // Apply regional resistance surveillance penalty for empirical therapy
                if !has_any_identified_infection {
                    let mut regional_resistance_penalty = 1.0_f64;
                    let region_idx = individual.region_cur_in as usize;
                    let hospital_status = individual.hospital_status.is_hospitalized();
                    
                    // Get configurable resistance penalty thresholds and penalties
                    let very_high_threshold = get_global_param("regional_resistance_threshold_very_high").unwrap_or(0.5);
                    let high_threshold = get_global_param("regional_resistance_threshold_high").unwrap_or(0.3);
                    let moderate_threshold = get_global_param("regional_resistance_threshold_moderate").unwrap_or(0.1);
                    
                    let very_high_penalty = get_global_param("regional_resistance_penalty_very_high").unwrap_or(0.2);
                    let high_penalty = get_global_param("regional_resistance_penalty_high").unwrap_or(0.4);
                    let moderate_penalty = get_global_param("regional_resistance_penalty_moderate").unwrap_or(0.7);
                    
                    // For empirical therapy, clinicians consider local resistance patterns
                    // Check resistance rates for likely bacterial causes
                    for b_idx in 0..BACTERIA_LIST.len() {
                        // Get regional resistance data for this bacteria-drug combination
                        if let Some(resistance_values) = majority_r_positive_values_by_combo.get(&(region_idx, hospital_status, b_idx, drug_idx)) {
                            if !resistance_values.is_empty() {
                                // Calculate resistance prevalence: proportion of cases with resistance > 0
                                let resistance_cases = resistance_values.len() as f64;
                                
                                // Estimate total cases by checking all drugs for this bacteria in this region
                                // (This gives us denominator for prevalence calculation)
                                let mut total_cases_estimate = resistance_cases;
                                for d_idx in 0..DRUG_SHORT_NAMES.len() {
                                    if let Some(other_resistance_values) = majority_r_positive_values_by_combo.get(&(region_idx, hospital_status, b_idx, d_idx)) {
                                        total_cases_estimate = total_cases_estimate.max(other_resistance_values.len() as f64);
                                    }
                                }
                                
                                if total_cases_estimate > 0.0 {
                                    let resistance_prevalence = resistance_cases / total_cases_estimate;
                                    
                                    // Apply graduated penalties based on regional resistance levels
                                    let resistance_penalty = if resistance_prevalence >= very_high_threshold {
                                        very_high_penalty  // Very high resistance - avoid drug
                                    } else if resistance_prevalence >= high_threshold {
                                        high_penalty       // High resistance - large penalty
                                    } else if resistance_prevalence >= moderate_threshold {
                                        moderate_penalty   // Moderate resistance - moderate penalty
                                    } else {
                                        1.0                // Low resistance - no penalty
                                    };
                                    
                                    // Use the most restrictive penalty across all bacteria
                                    regional_resistance_penalty = regional_resistance_penalty.min(resistance_penalty);
                                }
                            }
                        }
                    }
                    score *= regional_resistance_penalty;
                }

                let drug_spectrum = get_drug_param(drug_name, "spectrum_breadth").unwrap_or(3.0);
                if has_any_identified_infection {
                    let targeted_narrow_bonus = get_global_param("targeted_therapy_narrow_spectrum_bonus").unwrap_or(3.0);
                    let targeted_broad_penalty = get_global_param("targeted_therapy_broad_spectrum_penalty").unwrap_or(0.4);
                    let ineffective_drug_penalty = get_global_param("targeted_therapy_ineffective_drug_penalty").unwrap_or(0.1);
                    let effective_potency_threshold = get_global_param("effective_potency_threshold_for_targeted_therapy").unwrap_or(0.10);
                    
                    let mut has_good_activity = false;
                    let mut best_potency: f64 = 0.0;
                    for b_idx in 0..BACTERIA_LIST.len() {
                        if individual.test_identified_infection[b_idx] && individual.level[b_idx] > 0.001 {
                            let potency_param_key = &param_cache.drug_bacteria_potency_keys[&(drug_idx, b_idx)];
                            let potency = get_global_param(potency_param_key).unwrap_or(0.0);
                            best_potency = best_potency.max(potency);
                            if potency > effective_potency_threshold {
                                has_good_activity = true;
                            }
                        }
                    }
                    if has_good_activity {
                        if drug_spectrum <= 2.5 {
                            score *= targeted_narrow_bonus;
                        } else if drug_spectrum >= 4.0 {
                            score *= targeted_broad_penalty;
                        }
                    } else {
                        score *= ineffective_drug_penalty;
                    }
                } else if has_any_infection {
                    // Empirical therapy: check potency against actual infecting bacteria
                    // Even in empirical therapy, shouldn't choose drugs ineffective against the pathogen
                    let empiric_broad_bonus = get_global_param("empiric_therapy_broad_spectrum_bonus").unwrap_or(2.0);
                    let empiric_ineffective_penalty = get_global_param("empiric_therapy_ineffective_drug_penalty").unwrap_or(0.05);
                    let effective_potency_threshold = get_global_param("effective_potency_threshold_for_empirical_therapy").unwrap_or(0.10);
                    
                    let mut has_any_activity = false;
                    for b_idx in 0..BACTERIA_LIST.len() {
                        if individual.level[b_idx] > 0.001 {
                            let potency_param_key = &param_cache.drug_bacteria_potency_keys[&(drug_idx, b_idx)];
                            let potency = get_global_param(potency_param_key).unwrap_or(0.0);
                            if potency > effective_potency_threshold {
                                has_any_activity = true;
                                break;
                            }
                        }
                    }
                    
                    if has_any_activity {
                        if drug_spectrum >= 3.5 {
                            score *= empiric_broad_bonus;
                        } else if drug_spectrum <= 2.0 {
                            score *= 0.6;
                        }
                    } else {
                        // Drug has no activity against any infecting bacteria - heavily penalize
                        score *= empiric_ineffective_penalty;
                    }
                }
                // Apply drug availability multiplier
                let current_region_str = individual.region_cur_in.to_string();
                let living_region_str = individual.region_living.to_string();
                let drug_availability = get_drug_availability_time_aware(
                    drug_name, 
                    &current_region_str, 
                    Some(&living_region_str),
                    time_step
                );
                
                // Check if drug has been introduced yet
                let mut drug_introduced = false;
                if let Some(intro_time) = crate::config::get_drug_introduction_time_step(drug_name) {
                    if time_step >= intro_time {
                        drug_introduced = true;
                    }
                }
                // If no introduction date specified, assume drug is NOT available yet
                // This prevents anachronistic drug use
                
                score *= drug_availability;
                if !drug_introduced {
                    score = 0.0; // Drug not yet introduced, can't be prescribed
                }

                // Store drug score for the primary bacteria
                if primary_bacteria_idx >= 0 {
                    individual.drug_score_on_selection_day[drug_idx] = score;
                }

                // Only include drugs with positive scores for selection
                if score > 0.0 {
                    drug_scores.push((drug_idx, score));
                }
            }

            // Weighted probabilistic selection from scored drugs
            if !drug_scores.is_empty() {
                // Add stochasticity parameter to control randomness vs determinism
                let selection_temperature = get_global_param("drug_selection_temperature").unwrap_or(0.5); // MUCH more deterministic
                
                // Apply temperature scaling: lower temp = more deterministic (clinically realistic)
                // Temperature of 0.5 = strongly favor best drugs, 1.0 = moderate, 2.0+ = random
                let weights: Vec<f64> = drug_scores.iter()
                    .map(|(_, score)| (score / selection_temperature).exp())
                    .collect();
                
                // Handle edge case where all weights are zero or infinite
                let total_weight: f64 = weights.iter().sum();
                if total_weight > 0.0 && total_weight.is_finite() {
                    let dist = WeightedIndex::new(&weights).unwrap();
                    let chosen_idx = dist.sample(&mut rng);
                    let chosen_drug_idx = drug_scores[chosen_idx].0;
                    
                    // Initiate the selected drug
                    let drug_name = DRUG_SHORT_NAMES[chosen_drug_idx];
                    individual.cur_use_drug[chosen_drug_idx] = true;
                    individual.date_drug_initiated[chosen_drug_idx] = time_step as i32;
                    individual.date_drug_initiated_keep[chosen_drug_idx] = time_step as i32; // Persistent record
                    individual.ever_taken_drug[chosen_drug_idx] = true;
                    
                    // Update drug counter
                    update_drug_counter(individual);
                    if individual.id == 1000001  {
                        println!(
                            "mod.rs   started {} - two-stage rate of starting was {:.4} (score: {:.3})",
                            drug_name,
                            start_any_antibiotic_prob,
                            drug_scores[chosen_idx].1
                        );
                    }
                    let mut chosen_initial_level = get_drug_param(drug_name, "initial_level").unwrap_or(10.0);
                    if has_any_identified_infection && rng.gen_bool(double_dose_probability) {
                        let double_dose_multiplier = get_drug_param(drug_name, "double_dose_multiplier").unwrap_or(2.0);
                        chosen_initial_level *= double_dose_multiplier;
                    }
                    individual.cur_level_drug[chosen_drug_idx] = chosen_initial_level;
                    
                    // Update treatment failure tracking for all infected bacteria
                    for bacteria_idx in 0..BACTERIA_LIST.len() {
                        if individual.level[bacteria_idx] > 0.0 {
                            // Record bacteria level at drug start and reset tracking
                            individual.bacteria_level_at_drug_start[bacteria_idx] = Some(individual.level[bacteria_idx]);
                            individual.days_on_current_treatment[bacteria_idx] = 0;
                            individual.treatment_failure_assessed[bacteria_idx] = false;
                        }
                    }
                }
            }
        }
    }

    // drug-specific toxicity
    let mut daily_drug_toxicity_increase = 0.0;
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];
        if individual.cur_level_drug[drug_idx] > 0.0 {
            let drug_toxicity_per_unit = get_drug_param(drug_name, "toxicity_per_unit_level_per_day")
                .unwrap_or_else(|| get_global_param("default_drug_toxicity_per_unit_level_per_day")
                .expect("Missing default_drug_toxicity_per_unit_level_per_day in config"));
            daily_drug_toxicity_increase += individual.cur_level_drug[drug_idx] * drug_toxicity_per_unit;
        }
    }
    individual.current_toxicity = (individual.current_toxicity + daily_drug_toxicity_increase).max(0.0);

    // --- Treatment failure tracking and assessment ---
    // Update treatment days counter and assess treatment failure
    for bacteria_idx in 0..BACTERIA_LIST.len() {
        // Only track treatment days if there's an active infection
        if individual.level[bacteria_idx] > 0.0 {
            // Increment treatment days if we have recorded a drug start
            if individual.bacteria_level_at_drug_start[bacteria_idx].is_some() {
                individual.days_on_current_treatment[bacteria_idx] += 1;
                
                // Assess treatment failure if conditions are met
                assess_treatment_failure(
                    individual,
                    time_step,
                    bacteria_idx,
                    bacteria_indices,
                    drug_indices,
                    cross_resistance_groups,
                    param_cache,
                );
            }
        } else {
            // No active infection - reset all tracking
            individual.bacteria_level_at_drug_start[bacteria_idx] = None;
            individual.days_on_current_treatment[bacteria_idx] = -1;
            individual.treatment_failure_assessed[bacteria_idx] = false;
            
            // Also clear restart window tracking since infection has resolved
            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
            individual.stopped_drug_index[bacteria_idx] = None;
            individual.restart_window_assessed[bacteria_idx] = false;
        }
        
        // Assess restart window (independent of current infection status)
        assess_restart_window(
            individual,
            time_step,
            bacteria_idx,
            bacteria_indices,
            param_cache,
        );
    }

    // --- death     


    if individual.date_of_death.is_none() {
        let mut cause: Option<String> = None;

        // --- New Logistic Background Mortality Model ---
        let baseline_log_odds = get_global_param("background_mortality_baseline_log_odds")
            .expect("Missing background_mortality_baseline_log_odds in config");
       
        let mut total_log_odds = baseline_log_odds;

        // Time-varying mortality component (1930-2035): reflects historical mortality decline
        let years_since_1930 = time_step as f64 / 365.0;
        let start_multiplier = get_global_param("mortality_baseline_1930_multiplier").unwrap_or(3.0);
        let end_multiplier = get_global_param("mortality_baseline_2035_multiplier").unwrap_or(1.0);
        let half_life_years = get_global_param("mortality_improvement_half_life_years").unwrap_or(35.0);
        
        // Exponential decay from start_multiplier to end_multiplier
        let decay_rate = (2.0_f64).ln() / half_life_years; // ln(2) / half_life
        let time_multiplier = end_multiplier + (start_multiplier - end_multiplier) * (-decay_rate * years_since_1930).exp();
        let time_log_odds_adjustment = time_multiplier.ln();
        total_log_odds += time_log_odds_adjustment;

        let age_years = individual.age as f64 / 365.0;

        // Age effects
        let log_odds_per_year = get_global_param("log_odds_mortality_per_year_of_age")
            .expect("Missing log_odds_mortality_per_year_of_age in config");
        total_log_odds += age_years * log_odds_per_year;

        // Non-linear age effect for very elderly
        if age_years > 80.0 {
            let log_odds_age_squared = get_global_param("log_odds_mortality_per_year_of_age_squared")
                .unwrap_or(0.0);
            total_log_odds += (age_years - 80.0).powi(2) * log_odds_age_squared;
        }

        // Regional effects
        let region_name_for_param = individual.region_living.to_string().to_lowercase().replace(" ", "_");
        let region_log_odds_key = &param_cache.region_mortality_keys[&region_name_for_param];
        total_log_odds += get_global_param(region_log_odds_key).unwrap_or(0.0);

        // Sex effects
        let sex_log_odds_key = &param_cache.sex_mortality_keys[&individual.sex_at_birth.to_lowercase()];
        total_log_odds += get_global_param(sex_log_odds_key).unwrap_or(0.0);

        // Immunosuppression effect
        if individual.immunodeficiency_type.is_some() {
            total_log_odds += get_global_param("log_odds_mortality_immunosuppressed").unwrap_or(0.0);
        }

        // Hospital status effect
        if matches!(individual.hospital_status, HospitalStatus::InHospital) {
            total_log_odds += get_global_param("log_odds_mortality_hospitalized").unwrap_or(0.0);
        }

        // Convert total log odds to probability
        let background_risk = 1.0 / (1.0 + (-total_log_odds).exp());
       
        individual.background_all_cause_mortality_rate = background_risk.min(1.0);
        let mut prob_not_dying = 1.0 - background_risk;


        let has_sepsis = individual.sepsis.iter().any(|&status| status);
        if has_sepsis {
            // Calculate age-adjusted sepsis mortality risk
            let base_sepsis_death_risk = get_global_param("base_sepsis_death_risk_per_day")
                .expect("missing base_sepsis_death_risk_per_day in config");
            
            let mut sepsis_death_risk = base_sepsis_death_risk;
            
            // Apply age-based multiplier
            let age_years = individual.age as f64 / 365.0;
            let age_multiplier = if age_years < 1.0 {
                get_global_param("sepsis_age_mortality_multiplier_infant").unwrap_or(3.0)
            } else if age_years < 18.0 {
                get_global_param("sepsis_age_mortality_multiplier_child").unwrap_or(0.5)
            } else if age_years < 65.0 {
                get_global_param("sepsis_age_mortality_multiplier_adult").unwrap_or(1.0)
            } else {
                get_global_param("sepsis_age_mortality_multiplier_elderly").unwrap_or(2.5)
            };
            sepsis_death_risk *= age_multiplier;
            
            // Apply region-based multiplier (healthcare quality)
            let region_name_for_param = individual.region_living.to_string().to_lowercase().replace(" ", "_");
            let region_sepsis_multiplier_key = &param_cache.region_sepsis_multiplier_keys[&region_name_for_param];
            let region_sepsis_multiplier = get_global_param(region_sepsis_multiplier_key).unwrap_or(1.0);
            sepsis_death_risk *= region_sepsis_multiplier;
            
            // Apply immunosuppression multiplier
            if individual.immunodeficiency_type.is_some() {
                let immunosuppressed_multiplier = get_global_param("sepsis_immunosuppressed_multiplier").unwrap_or(3.0);
                sepsis_death_risk *= immunosuppressed_multiplier;
            }
            
            // Cap the risk at 1.0 (100%)
            sepsis_death_risk = sepsis_death_risk.min(1.0);
            
            prob_not_dying *= 1.0 - sepsis_death_risk;
            if cause.is_none() { cause = Some("sepsis_related".to_string()); }
        }
        let mut drug_adverse_event_risk_for_individual = 0.0;
        for drug_idx in 0..DRUG_SHORT_NAMES.len() {
            // Removed unused variable 'drug_name'
            if individual.cur_level_drug[drug_idx] > 0.0 {
                // Use only the global config parameter for drug toxicity death risk
                let drug_toxicity_death_risk = get_global_param("drug_toxicity_death_risk_per_day").unwrap_or(0.0);
                drug_adverse_event_risk_for_individual = (drug_adverse_event_risk_for_individual + drug_toxicity_death_risk).min(1.0);
            }
        }
        individual.mortality_risk_current_toxicity = drug_adverse_event_risk_for_individual;
        if drug_adverse_event_risk_for_individual > 0.0 {
            prob_not_dying *= 1.0 - drug_adverse_event_risk_for_individual;
            if cause.is_none() { cause = Some("drug_toxicity_related".to_string()); }
        }

        let mut prob_of_death_today = 1.0 - prob_not_dying;
        prob_of_death_today = prob_of_death_today.clamp(0.0, 1.0);
        if rng.gen::<f64>() < prob_of_death_today {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = cause.or(Some("background_mortality".to_string()));
            
            // Track death resolution for all current infections
            if let Some(ref death_cause) = individual.cause_of_death {
                let resolution_type = match death_cause.as_str() {
                    "sepsis_related" => InfectionResolutionType::DeathFromSepsis,
                    "drug_toxicity_related" => InfectionResolutionType::DeathFromToxicity,
                    _ => InfectionResolutionType::DeathFromBackground,
                };
                
                // Record resolution for ALL bacteria where person is currently infected
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[b_idx] > 0.001 {
                        let resolution_idx = match resolution_type {
                            InfectionResolutionType::ImmuneClearance => 0,
                            InfectionResolutionType::DrugAssistedClearance => 1,
                            InfectionResolutionType::DeathFromSepsis => 2,
                            InfectionResolutionType::DeathFromBackground => 3,
                            InfectionResolutionType::DeathFromToxicity => 4,
                        };
                        
                        individual.infection_resolution_this_timestep[b_idx][resolution_idx] += 1;
                    }
                }
            }
        }
    }
    // --- death logic end   


    // --- sepsis recovery logic (applied after death risk, only if individual is alive) ---
    if individual.date_of_death.is_none() {
        for &bacteria in BACTERIA_LIST.iter() {
            let b_idx = *bacteria_indices.get(bacteria).unwrap();
            
            // Only consider recovery if individual currently has sepsis from this bacteria
            if individual.sepsis[b_idx] {
                let sepsis_duration = (time_step as i32 - individual.sepsis_onset_day[b_idx]).max(0);
                let minimum_duration = get_global_param("sepsis_minimum_duration_days").unwrap_or(1.0) as i32;
                
                // Only allow recovery after minimum duration
                if sepsis_duration >= minimum_duration {
                    // Logistic regression model for sepsis recovery
                    let base_log_odds = get_global_param("sepsis_base_log_odds_of_recovery_per_day")
                        .expect("Missing sepsis_base_log_odds_of_recovery_per_day");
                    
                    let mut total_log_odds = base_log_odds;
                    
                    // (1) Bacteria level effect - higher bacteria level decreases recovery probability
                    let bacteria_level_coefficient = get_global_param("sepsis_log_odds_bacteria_level")
                        .expect("Missing sepsis_log_odds_bacteria_level");
                    total_log_odds += individual.level[b_idx] * bacteria_level_coefficient;
                    
                    // (2) Hospital status effect - being in hospital increases recovery probability
                    if individual.hospital_status.is_hospitalized() {
                        let hospital_coefficient = get_global_param("sepsis_log_odds_in_hospital")
                            .expect("Missing sepsis_log_odds_in_hospital");
                        total_log_odds += hospital_coefficient;
                    }
                    
                    // (3) Age effects with categories
                    let age_years = individual.age as f64 / 365.0;
                    let age_coefficient = if age_years < 1.0 {
                        get_global_param("sepsis_log_odds_age_infant").expect("Missing sepsis_log_odds_age_infant")
                    } else if age_years < 18.0 {
                        get_global_param("sepsis_log_odds_age_child").expect("Missing sepsis_log_odds_age_child")
                    } else if age_years < 65.0 {
                        get_global_param("sepsis_log_odds_age_adult").expect("Missing sepsis_log_odds_age_adult")
                    } else {
                        get_global_param("sepsis_log_odds_age_elderly").expect("Missing sepsis_log_odds_age_elderly")
                    };
                    total_log_odds += age_coefficient;
                    
                    // (4) Severe immunosuppression effect
                    if individual.immunodeficiency_type.is_some() {
                        let immunosuppressed_coefficient = get_global_param("sepsis_log_odds_immunosuppressed")
                            .expect("Missing sepsis_log_odds_immunosuppressed");
                        total_log_odds += immunosuppressed_coefficient;
                    }
                    
                    // (5) Region-specific effect (healthcare quality and ICU availability)
                    let region_name_for_param = individual.region_living.to_string().to_lowercase().replace(' ', "_");
                    let region_key = &param_cache.region_sepsis_keys[&region_name_for_param];
                    let region_coefficient = get_global_param(region_key).unwrap_or(0.0); // Default to 0.0 if region not found
                    total_log_odds += region_coefficient;
                    
                    // Convert log odds to probability using logistic function
                    let recovery_probability = 1.0 / (1.0 + (-total_log_odds).exp());
                    
                    // Check for recovery
                    if rng.gen::<f64>() < recovery_probability {
                        individual.sepsis[b_idx] = false;
                        // Keep sepsis_onset_day for tracking purposes (don't reset to -1)
                    }
                }
            }
        }
    }
    // --- end sepsis recovery logic ---

    // --- update per-bacteria fields ---
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        let is_infected = individual.level[b_idx] > 0.001;

        if !is_infected {
            // --- Logistic model for bacteria acquisition probability ---
            // All risk factors contribute additively to log-odds, then logistic function is applied.
            let mut log_odds = get_bacteria_param(bacteria, "acquisition_log_odds_baseline").unwrap_or_else(|| get_global_param("acquisition_log_odds_baseline").unwrap_or(-4.0));

            // Get age category and clean bacteria name (used multiple times)
            let age_category_str = crate::simulation::population::get_age_category_str(individual.age);
            let bacteria_clean = bacteria.replace(" ", "_");
            let region_name_for_param = individual.region_cur_in.to_string().to_lowercase().replace(" ", "_");

            // Age category effect (bacteria-specific with fallback to default)
            let age_bacteria_key = format!("{}_log_odds_{}", bacteria_clean, age_category_str);
            let log_odds_age_category = get_global_param(&age_bacteria_key)
                .unwrap_or_else(|| {
                    let default_age_key = format!("default_log_odds_{}", age_category_str);
                    get_global_param(&default_age_key).unwrap_or(0.0)
                });
            log_odds += log_odds_age_category;

            // Vaccination status (binary effect)
            if individual.vaccination_status[b_idx] {
                let log_odds_vaccinated = get_bacteria_param(bacteria, "log_odds_vaccinated").unwrap_or_else(|| get_global_param("log_odds_vaccinated").unwrap_or(0.0));
                log_odds += log_odds_vaccinated;
            }

            // Microbiome presence effect
            if individual.presence_microbiome[b_idx] {
                let log_odds_microbiome = get_bacteria_param(bacteria, "log_odds_microbiome_present").unwrap_or_else(|| get_global_param("log_odds_microbiome_present").unwrap_or(0.0));
                log_odds += log_odds_microbiome;
            }

            // Hospital-acquired effect
            if individual.hospital_status.is_hospitalized() {
                let log_odds_hospital = get_bacteria_param(bacteria, "log_odds_hospital_acquired").unwrap_or_else(|| get_global_param("log_odds_hospital_acquired").unwrap_or(0.0));
                log_odds += log_odds_hospital;
            }

            // Age-based effect (can be a function or coefficient)
            let log_odds_age = get_age_infection_multiplier(bacteria, individual.age); // If this returns log-odds, otherwise replace with appropriate function
            log_odds += log_odds_age;

            // Region-specific effect (bacteria-specific with fallback to general region effect)
            let region_bacteria_log_odds_key = &param_cache.region_bacteria_acquisition_keys[&(region_name_for_param.clone(), bacteria_clean.clone())];
            let region_bacteria_log_odds = get_global_param(region_bacteria_log_odds_key)
                .unwrap_or_else(|| {
                    let default_region_key = &param_cache.region_bacteria_default_keys[&region_name_for_param];
                    get_global_param(default_region_key).unwrap_or(0.0)
                });
            log_odds += region_bacteria_log_odds;

            // Age-Region interaction effect (bacteria-specific with fallback to general age-region)
            let age_region_bacteria_key = format!("{}_{}_log_odds_{}", bacteria_clean, region_name_for_param, age_category_str);
            let log_odds_age_region_bacteria = get_global_param(&age_region_bacteria_key)
                .unwrap_or_else(|| {
                    // Fallback to general age-region interaction (not bacteria-specific)
                    let age_region_key = format!("{}_log_odds_{}", region_name_for_param, age_category_str);
                    get_global_param(&age_region_key).unwrap_or(0.0)
                });
            log_odds += log_odds_age_region_bacteria;

            // Convert log-odds to probability
            let acquisition_probability = 1.0 / (1.0 + (-log_odds).exp());

            // --- microbiome presence (Carriage) ---
            if !individual.presence_microbiome[b_idx] {
                // Use the same log-odds formula as infection acquisition, with an extra microbiome-vs-infection log-odds parameter
                let mut log_odds = get_bacteria_param(bacteria, "acquisition_log_odds_baseline").unwrap_or_else(|| get_global_param("acquisition_log_odds_baseline").unwrap_or(-4.0));

                // Age category effect (reuse variables from above)
                let log_odds_age_category = get_global_param(&age_bacteria_key)
                    .unwrap_or_else(|| {
                        let default_age_key = format!("default_log_odds_{}", age_category_str);
                        get_global_param(&default_age_key).unwrap_or(0.0)
                    });
                log_odds += log_odds_age_category;

                // Vaccination status (binary effect)
                if individual.vaccination_status[b_idx] {
                    let log_odds_vaccinated = get_bacteria_param(bacteria, "log_odds_vaccinated").unwrap_or_else(|| get_global_param("log_odds_vaccinated").unwrap_or(0.0));
                    log_odds += log_odds_vaccinated;
                }

                // Hospital-acquired effect
                if individual.hospital_status.is_hospitalized() {
                    let log_odds_hospital = get_bacteria_param(bacteria, "log_odds_hospital_acquired").unwrap_or_else(|| get_global_param("log_odds_hospital_acquired").unwrap_or(0.0));
                    log_odds += log_odds_hospital;
                }

                // Age-based effect (can be a function or coefficient)
                let log_odds_age = get_age_infection_multiplier(bacteria, individual.age);
                log_odds += log_odds_age;

                // Region-specific effect (reuse variables from above)
                let region_bacteria_log_odds_key = &param_cache.region_bacteria_acquisition_keys[&(region_name_for_param.clone(), bacteria_clean.clone())];
                let region_bacteria_log_odds = get_global_param(region_bacteria_log_odds_key)
                    .unwrap_or_else(|| {
                        let default_region_key = &param_cache.region_bacteria_default_keys[&region_name_for_param];
                        get_global_param(default_region_key).unwrap_or(0.0)
                    });
                log_odds += region_bacteria_log_odds;

                // Age-Region interaction effect (reuse variables from above)
                let log_odds_age_region_bacteria = get_global_param(&age_region_bacteria_key)
                    .unwrap_or_else(|| {
                        // Fallback to general age-region interaction (not bacteria-specific)
                        let age_region_key = format!("{}_log_odds_{}", region_name_for_param, age_category_str);
                        get_global_param(&age_region_key).unwrap_or(0.0)
                    });
                log_odds += log_odds_age_region_bacteria;

                // Add the extra log-odds for microbiome vs infection (bacteria-specific)
                let microbiome_vs_infection_key = format!("{}_log_odds_microbiome_vs_infection", bacteria_clean);
                let log_odds_microbiome_vs_infection = get_global_param(&microbiome_vs_infection_key)
                    .unwrap_or_else(|| get_global_param("log_odds_microbiome_vs_infection").unwrap_or(-6.0)); // Fallback to old global param if bacteria-specific not found
                log_odds += log_odds_microbiome_vs_infection;

                // Convert log-odds to probability
                let microbiome_acquisition_probability = 1.0 / (1.0 + (-log_odds).exp());

                if rng.gen_bool(microbiome_acquisition_probability.clamp(0.0, 1.0)) {
                    individual.presence_microbiome[b_idx] = true;

                    // --- assign microbiome_r on new microbiome acquisition (same logic as infection resistance assignment) ---
                    let env_majority_r_level = get_global_param("environmental_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                    let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0);

                    let is_from_environment = true; // Microbiome acquisition is always from environment in this model
                    let is_hospital_acquired = individual.hospital_status.is_hospitalized();

                    let region_idx = individual.region_cur_in as usize;
                    let hospital_status_bool = individual.hospital_status.is_hospitalized();

                    for drug_name_static in DRUG_SHORT_NAMES.iter() {
                        let d_idx = *drug_indices.get(drug_name_static).unwrap();
                        let resistance_data = &mut individual.resistances[b_idx][d_idx];

                        if is_from_environment {
                            resistance_data.microbiome_r = env_majority_r_level;
                        } else {
                            // --- region/hospital-specific sampling for microbiome (same logic as infections) ---
                            let sampling_hospital_status = if is_hospital_acquired {
                                true // Hospital-acquired microbiome samples from hospitalized population
                            } else {
                                hospital_status_bool // Community-acquired microbiome samples based on current status
                            };

                            if let Some(majority_r_values_from_population) =
                                majority_r_positive_values_by_combo.get(&(region_idx, sampling_hospital_status, b_idx, d_idx))
                            {
                                if let Some(&acquired_resistance_level) = majority_r_values_from_population.choose(&mut rng) {
                                    let clamped_level = acquired_resistance_level.min(max_resistance_level).max(0.0);
                                    resistance_data.microbiome_r = clamped_level;
                                } else {
                                    resistance_data.microbiome_r = 0.0;
                                }
                            } else {
                                resistance_data.microbiome_r = 0.0;
                            }
                        }
                    }
                    // --- end microbiome_r assignment ---
                }
            } else {
                let microbiome_clearance_prob = get_bacteria_param(bacteria, "microbiome_clearance_probability_per_day")
                    .unwrap_or_else(|| get_global_param("default_microbiome_clearance_probability_per_day").expect("Missing default_microbiome_clearance_probability_per_day in config"));
                if rng.gen_bool(microbiome_clearance_prob) {
                    individual.presence_microbiome[b_idx] = false;
                }

                // --- de novo resistance emergence in microbiome when on drug ---
                if individual.presence_microbiome[b_idx] {
                    let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0);
                    for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                        let resistance_data = &mut individual.resistances[b_idx][d_idx];
                        let drug_level = individual.cur_level_drug[d_idx];
                        // Only consider emergence if drug is present and microbiome_r is low
                        if drug_level > 0.0001 && resistance_data.microbiome_r < 0.0001 {
                            // Use a specific parameter for microbiome resistance emergence if present, else fallback to general
                            let emergence_rate_baseline = get_global_param("microbiome_resistance_emergence_rate_per_day_baseline")
                                .or_else(|| None)
                                .unwrap_or(0.000001);
                            let microbiome_r_emergence_level = get_global_param("any_r_emergence_level_on_first_emergence").unwrap_or(0.5);

                            // Optionally, you could scale by drug level or other factors
                            let total_emergence_prob = emergence_rate_baseline; // * (drug_level / 10.0).clamp(0.0, 1.0);

                            if rng.gen_bool(total_emergence_prob.clamp(0.0, 1.0)) {
                                resistance_data.microbiome_r = microbiome_r_emergence_level.min(max_resistance_level);
                            }
                        }
                    }
                }
                // --- end de novo resistance emergence in microbiome ---
            }

            // ...resistance transfer (each way) between infection site and microbiome ...
            for &drug in DRUG_SHORT_NAMES.iter() {
                let d_idx = *drug_indices.get(drug).unwrap();
                if !individual.presence_microbiome[b_idx] {
                    individual.resistances[b_idx][d_idx].microbiome_r = 0.0;
                } else {
                    let infection_present = individual.level[b_idx] > 0.0;
                    if infection_present {
                        let current_any_r = individual.resistances[b_idx][d_idx].any_r;
                        let current_microbiome_r = individual.resistances[b_idx][d_idx].microbiome_r;
                        let possible_transfer_r_microbiome = (current_any_r > 0.0 && current_microbiome_r == 0.0) ||
                                                     (current_microbiome_r > 0.0 && current_any_r == 0.0);
                        if possible_transfer_r_microbiome && rng.gen_bool(transfer_prob) {
                            if current_any_r > 0.0 && current_microbiome_r == 0.0 {
                                individual.resistances[b_idx][d_idx].microbiome_r = current_any_r;
                            } else if current_microbiome_r > 0.0 && current_any_r == 0.0 {
                                individual.resistances[b_idx][d_idx].any_r = current_microbiome_r;
                            }
                        }
                    }
                }
            }

            // --- HORIZONTAL GENE TRANSFER (HGT) BETWEEN DIFFERENT BACTERIA ---
            // For each donor bacteria (with resistance), try to transfer to each other recipient bacteria
            for donor_idx in 0..BACTERIA_LIST.len() {
                // Donor must have resistance (infection or microbiome)
                let donor_has_resistance = individual.level[donor_idx] > 0.001 || individual.presence_microbiome[donor_idx];
                if donor_has_resistance {
                    for recipient_idx in 0..BACTERIA_LIST.len() {
                        if recipient_idx == donor_idx { continue; }
                        let donor_name = BACTERIA_LIST[donor_idx];
                        let recipient_name = BACTERIA_LIST[recipient_idx];
                        let hgt_key = &param_cache.hgt_keys[&(donor_name.to_string(), recipient_name.to_string())];
                        let hgt_prob = crate::config::PARAMETERS.get(hgt_key).copied().unwrap_or(0.0);
                        if hgt_prob > 0.0 && rng.gen::<f64>() < hgt_prob {
                            // Transfer resistance for all drugs
                            for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                                let donor_r = individual.resistances[donor_idx][drug_idx].any_r;
                                if donor_r > 0.0 {
                                    // Transfer to infection
                                    if individual.level[recipient_idx] > 0.001 {
                                        let prev_any_r = individual.resistances[recipient_idx][drug_idx].any_r;
                                        let new_any_r = donor_r.max(prev_any_r);
                                        individual.resistances[recipient_idx][drug_idx].any_r = new_any_r;
                                        if prev_any_r == 0.0 && new_any_r > 0.0 {
                                            // Inline mechanism assignment
                                            use crate::simulation::population::ResistanceMechanism;
                                            let mechanism_prob = get_global_param("mechanism_assignment_probability_on_any_r_gain").unwrap_or(0.8);
                                            for (mech_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                                                let mechanism_str = mechanism.as_str();
                                                let enhancement = get_global_param(&param_cache.resistance_mechanism_enhancement_keys[mechanism_str]).unwrap_or(0.0);
                                                if enhancement <= new_any_r {
                                                    if rng.gen_bool(mechanism_prob) {
                                                        individual.resistance_mechanisms[recipient_idx][mech_idx] = true;
                                                    }
                                                }
                                            }
                                            individual.how_resistance_acquired[recipient_idx][drug_idx] = Some(crate::simulation::population::ResistanceAcquisitionType::Hgt);
                                        }
                                    }
                                    // Transfer to microbiome
                                    if individual.presence_microbiome[recipient_idx] {
                                        let prev_any_r = individual.resistances[recipient_idx][drug_idx].any_r;
                                        let new_any_r = donor_r.max(prev_any_r);
                                        individual.resistances[recipient_idx][drug_idx].any_r = new_any_r;
                                        if prev_any_r == 0.0 && new_any_r > 0.0 {
                                            // Inline mechanism assignment
                                            use crate::simulation::population::ResistanceMechanism;
                                            let mechanism_prob = get_global_param("mechanism_assignment_probability_on_any_r_gain").unwrap_or(0.8);
                                            for (mech_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                                                let mechanism_str = mechanism.as_str();
                                                let enhancement = get_global_param(&param_cache.resistance_mechanism_enhancement_keys[mechanism_str]).unwrap_or(0.0);
                                                if enhancement <= new_any_r {
                                                    if rng.gen_bool(mechanism_prob) {
                                                        individual.resistance_mechanisms[recipient_idx][mech_idx] = true;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                // Check if existing antibiotic therapy prevents this infection
                let mut infection_prevented = false;
                let prevention_efficacy = get_global_param("antibiotic_infection_prevention_efficacy").unwrap_or(0.85);
                
                // Check each drug the person is currently taking
                for (drug_idx, &is_taking_drug) in individual.cur_use_drug.iter().enumerate() {
                    if is_taking_drug {
                        // Calculate effective activity using the same method as activity_r calculation
                        let potency_param_key = &param_cache.drug_bacteria_potency_keys[&(drug_idx, b_idx)];
                        let base_potency = get_global_param(potency_param_key).unwrap_or(0.05);
                        let drug_current_level = individual.cur_level_drug[drug_idx];
                        let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0);
                        let resistance_level = individual.resistances[b_idx][drug_idx].any_r;
                        let normalized_any_r = resistance_level / max_resistance_level;
                        let effective_activity = base_potency * drug_current_level * (1.0 - normalized_any_r);
                        
                        // If drug has effective activity, it can prevent infection
                        if effective_activity > 0.5 { // Threshold for effective prevention
                            if rng.gen_bool(prevention_efficacy) {
                                infection_prevented = true;
                                break; // One effective drug is enough
                            }
                        }
                    }
                }
                
                // Only proceed with infection if not prevented by existing antibiotics
                if !infection_prevented {
                    let initial_level = get_bacteria_param(bacteria, "initial_infection_level").unwrap_or(0.01);
                    individual.level[b_idx] = initial_level;
                    individual.date_last_infected[b_idx] = time_step as i32;
                    individual.date_last_infected_keep[b_idx] = time_step as i32; // Keep persistent record

                    // --- probabilistic syndrome assignment ---
                    let syndrome_id = assign_syndrome_for_bacteria(bacteria, &mut rng);
                    individual.infectious_syndrome[b_idx] = syndrome_id as i32;

                let env_acquisition_chance = get_bacteria_param(bacteria, "environmental_acquisition_proportion").unwrap_or(0.1);
                individual.cur_infection_from_environment[b_idx] = rng.gen::<f64>() < env_acquisition_chance;

                individual.infection_hospital_acquired[b_idx] = individual.hospital_status.is_hospitalized();

                // --- any_r and majority_r setting logic on new infection acquisition ---
                let env_majority_r_level = get_global_param("environmental_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0);


                let is_from_environment = individual.cur_infection_from_environment[b_idx];
                let is_hospital_acquired = individual.infection_hospital_acquired[b_idx];

                let region_idx = individual.region_cur_in as usize;
                let hospital_status_bool = individual.hospital_status.is_hospitalized();

                for drug_name_static in DRUG_SHORT_NAMES.iter() {
                    let d_idx = *drug_indices.get(drug_name_static).unwrap();
                    let resistance_data = &mut individual.resistances[b_idx][d_idx];

                    if is_from_environment {
                        // Check if any drug that selects for resistance to this drug has been introduced
                        let mut any_selecting_drug_introduced = false;
                        
                        // Check if the drug itself has been introduced
                        if let Some(intro_time) = crate::config::get_drug_introduction_time_step(drug_name_static) {
                            if time_step >= intro_time {
                                any_selecting_drug_introduced = true;
                            }
                        }
                        
                        // If not yet introduced by direct drug, check cross-resistance groups
                        if !any_selecting_drug_introduced {
                            if let Some(cross_resistance_drug_groups) = cross_resistance_groups.get(&b_idx) {
                                for group in cross_resistance_drug_groups {
                                    if group.contains(&d_idx) {
                                        // This drug is in a cross-resistance group, check if any other drug in the group has been introduced
                                        for &other_drug_idx in group {
                                            if other_drug_idx != d_idx {
                                                if let Some(other_drug_name) = DRUG_SHORT_NAMES.get(other_drug_idx) {
                                                    if let Some(intro_time) = crate::config::get_drug_introduction_time_step(other_drug_name) {
                                                        if time_step >= intro_time {
                                                            any_selecting_drug_introduced = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if any_selecting_drug_introduced {
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        
                        // Only assign environmental resistance if a selecting drug has been introduced
                        if any_selecting_drug_introduced {
                            resistance_data.majority_r = env_majority_r_level;
                            resistance_data.any_r = env_majority_r_level;
                            // Inline mechanism assignment
                            use crate::simulation::population::ResistanceMechanism;
                            let mechanism_prob = get_global_param("mechanism_assignment_probability_on_any_r_gain").unwrap_or(0.8);
                            for (mech_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                                let mechanism_str = mechanism.as_str();
                                let enhancement = get_global_param(&param_cache.resistance_mechanism_enhancement_keys[mechanism_str]).unwrap_or(0.0);
                                if enhancement <= resistance_data.any_r {
                                    if rng.gen_bool(mechanism_prob) {
                                        individual.resistance_mechanisms[b_idx][mech_idx] = true;
                                    }
                                }
                            }
                            individual.how_resistance_acquired[b_idx][d_idx] = Some(crate::simulation::population::ResistanceAcquisitionType::AtInfectionEnv);
                        } else {
                            resistance_data.majority_r = 0.0;
                            resistance_data.any_r = 0.0;
                        }
                    } else {
                        // --- region/hospital-specific sampling for both hospital-acquired and community-acquired ---
                        // For hospital-acquired infections, we sample from hospitalized people (hospital_status_bool = true)
                        // For community-acquired infections, we sample based on the person's current hospital status
                        let sampling_hospital_status = if is_hospital_acquired {
                            true // Hospital-acquired infections sample from hospitalized population
                        } else {
                            hospital_status_bool // Community-acquired infections sample based on current status
                        };

                        if let Some(majority_r_values_from_population) =
                            majority_r_positive_values_by_combo.get(&(region_idx, sampling_hospital_status, b_idx, d_idx))
                        {
                            if let Some(&acquired_resistance_level) = majority_r_values_from_population.choose(&mut rng) {
                                let clamped_level = acquired_resistance_level.min(max_resistance_level).max(0.0);
                                resistance_data.any_r = clamped_level;
                                resistance_data.majority_r = clamped_level;
                                // Inline mechanism assignment
                                use crate::simulation::population::ResistanceMechanism;
                                let mechanism_prob = get_global_param("mechanism_assignment_probability_on_any_r_gain").unwrap_or(0.8);
                                for (mech_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                                    let mechanism_str = mechanism.as_str();
                                    let enhancement = get_global_param(&param_cache.resistance_mechanism_enhancement_keys[mechanism_str]).unwrap_or(0.0);

                                    if enhancement <= resistance_data.any_r {
                                        if rng.gen_bool(mechanism_prob) {
                                            individual.resistance_mechanisms[b_idx][mech_idx] = true;
                                        }
                                    }
                                }
                                individual.how_resistance_acquired[b_idx][d_idx] = Some(crate::simulation::population::ResistanceAcquisitionType::AtInfectionCommunity);
                            } else {
                                resistance_data.any_r = 0.0;
                                resistance_data.majority_r = 0.0;
                            }
                        } else {
                            resistance_data.any_r = 0.0;
                            resistance_data.majority_r = 0.0;
                        }
                    }
                }
                // --- end generalized any_r and majority_r setting logic ---
                } // End if !infection_prevented block
            } 
        } else { // Bacteria is already present (infection progression)
            // --- majority_r evolution ---
            let majority_r_evolution_rate = get_global_param("majority_r_evolution_rate_per_day_when_drug_present").unwrap_or(0.0);
            let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0); // Now using 1.0 from your config

            if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                for (drug_index, _use_drug) in individual.cur_use_drug.iter().enumerate() { 
                    let resistance_data = &mut individual.resistances[bacteria_full_idx][drug_index];

                    let drug_current_level = individual.cur_level_drug[drug_index];
                    let drug_currently_present = drug_current_level > 0.0; // Check if drug is effectively present
                    let current_bacteria_level = individual.level[b_idx];

                    // existing majority_r evolution based on drug presence
                    if resistance_data.majority_r == 0.0 && resistance_data.any_r > 0.0 && drug_currently_present {
                        if rng.gen_bool(majority_r_evolution_rate) {
                            resistance_data.majority_r = resistance_data.any_r;
                        }
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


                    // majority_r and any_r between 0 and 1
                    resistance_data.majority_r = resistance_data.majority_r.min(max_resistance_level).max(0.0);
                    resistance_data.any_r = resistance_data.any_r.min(max_resistance_level).max(0.0);


                    //new resistance emergence ---
                    // this section handles the de novo emergence of resistance when it's not already present.
                    // it should come before activity_r is fully calculated for use in bacteria level reduction *this* time step.
                    
                    if resistance_data.any_r < 0.0001 { // Check if any_r is effectively zero
                        // only consider emergence if there's drug present (either being taken or decaying)
                        // and a positive bacteria level for selection pressure.
                        if drug_current_level > 0.0 && current_bacteria_level > 0.0001 { 
                            let param_key = format![
                                "drug_{}_for_bacteria_{}_resistance_emergence_rate_per_day_baseline",
                                DRUG_SHORT_NAMES[drug_index],
                                bacteria
                            ];
                            let emergence_rate_baseline = get_global_param(&param_key).unwrap_or(0.000001); // Very small baseline
                            let bacteria_level_effect_multiplier = get_global_param("resistance_emergence_bacteria_level_multiplier").unwrap_or(0.05); // How much does bacteria level boost it
                            let any_r_emergence_level_on_first_emergence = get_global_param("any_r_emergence_level_on_first_emergence").unwrap_or(0.5); // User changed to 0.5 (was 1.0)

                            // bacteria level dependency: Higher at higher levels
                            let max_bacteria_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);
                            // Normalize bacteria level to [0,1] and apply multiplier
                            let bacteria_level_factor = (current_bacteria_level / max_bacteria_level).clamp(0.0, 1.0) * bacteria_level_effect_multiplier;
                            
                            // activity_r dependency: Bell-shaped curve
                            // Use the drug's initial level for normalization to get a comparable drug concentration scale (0-10)
                            let drug_initial_level_for_normalization = get_drug_param(DRUG_SHORT_NAMES[drug_index], "initial_level").unwrap_or(10.0);
                            
                            // Normalize current drug level for bell-shaped emergence probability curve
                            let mut norm_drug_level = drug_current_level / drug_initial_level_for_normalization;
                            norm_drug_level = norm_drug_level.clamp(0.0, 10.0); 
                            
                            // resistance emergence probability
                            // bell-shaped curve: 0.02 * x * (10 - x). Peaks at 5.0, is 0.1 at 0 and 10.
                            let emergence_drug_concentration_factor = 0.1 + 0.02 * norm_drug_level * (10.0 - norm_drug_level);
                            let emergence_drug_factor = emergence_drug_concentration_factor.clamp(0.0, 1.0);  

                            // Calculate multi-drug penalty if multiple drugs are active
                            let active_drug_count = individual.cur_level_drug.iter()
                                .filter(|&&level| level > 0.0)
                                .count();
                            
                            let multi_drug_penalty_threshold = get_global_param("multi_drug_penalty_threshold_num_drugs").unwrap_or(2.0) as usize;
                            let mut multi_drug_penalty_factor = 1.0;
                            let mut drugs_affected_by_this_resistance = 1; // Default: single drug resistance
                            
                            if active_drug_count >= multi_drug_penalty_threshold {
                                // Count how many active drugs this potential resistance would affect
                                drugs_affected_by_this_resistance = if let Some(cross_resistance_groups) = 
                                    crate::config::get_cross_resistance_groups().get(bacteria) {
                                    
                                    let current_drug_name = DRUG_SHORT_NAMES[drug_index];
                                    let mut affected_count = 0;
                                    
                                    // Check if this drug is in any cross-resistance group
                                    for group in cross_resistance_groups {
                                        if group.contains(&current_drug_name) {
                                            // Count how many drugs in this cross-resistance group are currently active
                                            for &group_drug in group {
                                                if let Some(group_drug_idx) = DRUG_SHORT_NAMES.iter().position(|&d| d == group_drug) {
                                                    if individual.cur_level_drug[group_drug_idx] > 0.0 {
                                                        affected_count += 1;
                                                    }
                                                }
                                            }
                                            break; // Found the group, no need to check others
                                        }
                                    }
                                    
                                    // If drug not in any cross-resistance group, resistance only affects this single drug
                                    if affected_count == 0 {
                                        affected_count = 1;
                                    }
                                    
                                    affected_count
                                } else {
                                    // No cross-resistance data for this bacteria, assume single drug resistance
                                    1
                                };
                                
                                // Apply penalty based on how many active drugs the resistance affects
                                if drugs_affected_by_this_resistance < active_drug_count {
                                    // Resistance doesn't affect all active drugs
                                    if drugs_affected_by_this_resistance == 1 {
                                        // Single drug resistance among multiple active drugs
                                        multi_drug_penalty_factor = get_global_param("resistance_development_inhibition_single_drug").unwrap_or(0.05);
                                    } else {
                                        // Partial cross-resistance among multiple active drugs
                                        multi_drug_penalty_factor = get_global_param("resistance_development_inhibition_partial_cross").unwrap_or(0.3);
                                    }
                                }
                                // If drugs_affected_by_this_resistance >= active_drug_count, no penalty (full cross-resistance)
                            }

                            // total emergence probability with multi-drug penalty
                            // adding 1.0 to bacteria_level_factor ensures a base contribution even if multiplier is low
                            let total_emergence_prob = emergence_rate_baseline * (1.0 + bacteria_level_factor) * emergence_drug_factor * multi_drug_penalty_factor;

                            if rng.gen_bool(total_emergence_prob.clamp(0.0, 1.0)) {
                                resistance_data.any_r = any_r_emergence_level_on_first_emergence;
                                // Inline mechanism assignment
                                use crate::simulation::population::ResistanceMechanism;
                                let mechanism_prob = get_global_param("mechanism_assignment_probability_on_any_r_gain").unwrap_or(0.8);
                                for (mech_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                                    let mechanism_str = mechanism.as_str();
                                    let enhancement = get_global_param(&param_cache.resistance_mechanism_enhancement_keys[mechanism_str]).unwrap_or(0.0);
                                    if enhancement <= resistance_data.any_r {
                                        if rng.gen_bool(mechanism_prob) {
                                            individual.resistance_mechanisms[b_idx][mech_idx] = true;
                                        }
                                    }
                                }
                                individual.how_resistance_acquired[b_idx][drug_index] = Some(crate::simulation::population::ResistanceAcquisitionType::FromMicrobiomeR);
                            }
                        }
                    }
                    // --- end new resistance emergence logic ---

                    // --- resistance mechanism emergence logic ---
                    // Check for emergence of specific resistance mechanisms when drug is present
                    if drug_current_level > 0.0 && current_bacteria_level > 0.0001 {
                        use crate::simulation::population::ResistanceMechanism;
                        
                        if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                            for (mechanism_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                                // Skip if mechanism already present
                                if individual.resistance_mechanisms[bacteria_full_idx][mechanism_idx] {
                                    continue;
                                }
                                
                                // Check if this mechanism is relevant for current drug
                                let mechanism_applicable = match (mechanism, DRUG_SHORT_NAMES[drug_index]) {
                                    // ESBL affects beta-lactams (except carbapenems)
                                    (ResistanceMechanism::ESBL, drug) => {
                                        matches!(drug, "penicilling" | "ampicillin" | "amoxicillin" | "piperacillin" | 
                                               "ticarcillin" | "cephalexin" | "cefazolin" | "cefuroxime" | 
                                               "ceftriaxone" | "ceftazidime" | "cefepime" | "ceftaroline" | 
                                               "aztreonam" | "amoxicillin_clavulanate" | "piperacillin_tazobactam" |
                                               "ampicillin_sulbactam" | "ticarcillin_clavulanate")
                                    },
                                    // Carbapenemase affects carbapenems
                                    (ResistanceMechanism::Carbapenemase, drug) => {
                                        matches!(drug, "meropenem" | "imipenem_c" | "ertapenem" | "meropenem_vaborbactam")
                                    },
                                    // 16S methyltransferase affects aminoglycosides
                                    (ResistanceMechanism::SixteenSMethyltransferase, drug) => {
                                        matches!(drug, "gentamicin" | "tobramycin" | "amikacin")
                                    },
                                    // Qnr affects quinolones
                                    (ResistanceMechanism::Qnr, drug) => {
                                        matches!(drug, "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin")
                                    },
                                    // Erm methylation affects macrolides
                                    (ResistanceMechanism::ErmMethylation, drug) => {
                                        matches!(drug, "erythromycin" | "azithromycin" | "clarithromycin")
                                    },
                                    // Van-type affects glycopeptides
                                    (ResistanceMechanism::VanType, drug) => {
                                        matches!(drug, "vancomycin" | "teicoplanin")
                                    },
                                    // mecA affects beta-lactams in Staph aureus
                                    (ResistanceMechanism::MecA, drug) => {
                                        bacteria == "staphylococcus aureus" && 
                                        matches!(drug, "penicilling" | "ampicillin" | "amoxicillin" | "cephalexin" | 
                                               "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | 
                                               "cefepime" | "meropenem" | "imipenem_c" | "ertapenem")
                                    },
                                    // Efflux overexpression can affect multiple drug classes
                                    (ResistanceMechanism::EffluxOverexpression, _) => true,
                                    // Reduced permeability affects many drugs, especially in Gram-negatives
                                    (ResistanceMechanism::ReducedPermeability, _) => {
                                        !matches!(bacteria, "staphylococcus aureus" | "streptococcus pneumoniae" | 
                                                 "streptococcus pyogenes" | "streptococcus agalactiae" | 
                                                 "enterococcus faecalis" | "enterococcus faecium")
                                    },
                                    // Target site mutations can affect various drugs
                                    (ResistanceMechanism::TargetSiteMutation, _) => true,
                                    // AmpC affects beta-lactams
                                    (ResistanceMechanism::AmpC, drug) => {
                                        matches!(drug, "penicilling" | "ampicillin" | "amoxicillin" | "piperacillin" | 
                                               "ticarcillin" | "cephalexin" | "cefazolin" | "cefuroxime" | 
                                               "ceftriaxone" | "amoxicillin_clavulanate" | "piperacillin_tazobactam" |
                                               "ampicillin_sulbactam" | "ticarcillin_clavulanate")
                                    },
                                };
                                
                                if mechanism_applicable {
                                    let mechanism_str = mechanism.as_str();
                                    let mechanism_emergence_rate = get_global_param(
                                        &param_cache.resistance_mechanism_emergence_keys[mechanism_str]
                                    ).unwrap_or(0.001);
                                    
                                    if rng.gen_bool(mechanism_emergence_rate.clamp(0.0, 1.0)) {
                                        individual.resistance_mechanisms[bacteria_full_idx][mechanism_idx] = true;
                                    }
                                }
                            }
                        }
                    }
                    // --- end resistance mechanism emergence logic ---


                    // calculate activity_r (should always be updated)
                    if drug_current_level > 0.0 {
                        // Fetch potency from config, fallback to 0.05 if not found
                        let potency_param_key = &param_cache.drug_bacteria_potency_keys[&(drug_index, bacteria_full_idx)];
                        let base_potency = get_global_param(potency_param_key).unwrap_or(0.05);
                        
                        // Calculate resistance mechanism enhancement
                        let mut mechanism_resistance_boost = 0.0;
                        if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                            use crate::simulation::population::ResistanceMechanism;
                            
                            for (mechanism_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                                if individual.resistance_mechanisms[bacteria_full_idx][mechanism_idx] {
                                    // Check if this mechanism affects the current drug
                                    let mechanism_affects_drug = match (mechanism, DRUG_SHORT_NAMES[drug_index]) {
                                        // ESBL affects beta-lactams (except carbapenems)
                                        (ResistanceMechanism::ESBL, drug) => {
                                            matches!(drug, "penicilling" | "ampicillin" | "amoxicillin" | "piperacillin" | 
                                                   "ticarcillin" | "cephalexin" | "cefazolin" | "cefuroxime" | 
                                                   "ceftriaxone" | "ceftazidime" | "cefepime" | "ceftaroline" | 
                                                   "aztreonam" | "amoxicillin_clavulanate" | "piperacillin_tazobactam" |
                                                   "ampicillin_sulbactam" | "ticarcillin_clavulanate")
                                        },
                                        // Carbapenemase affects carbapenems
                                        (ResistanceMechanism::Carbapenemase, drug) => {
                                            matches!(drug, "meropenem" | "imipenem_c" | "ertapenem" | "meropenem_vaborbactam")
                                        },
                                        // 16S methyltransferase affects aminoglycosides
                                        (ResistanceMechanism::SixteenSMethyltransferase, drug) => {
                                            matches!(drug, "gentamicin" | "tobramycin" | "amikacin")
                                        },
                                        // Qnr affects quinolones
                                        (ResistanceMechanism::Qnr, drug) => {
                                            matches!(drug, "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin")
                                        },
                                        // Erm methylation affects macrolides
                                        (ResistanceMechanism::ErmMethylation, drug) => {
                                            matches!(drug, "erythromycin" | "azithromycin" | "clarithromycin")
                                        },
                                        // Van-type affects glycopeptides
                                        (ResistanceMechanism::VanType, drug) => {
                                            matches!(drug, "vancomycin" | "teicoplanin")
                                        },
                                        // mecA affects beta-lactams in Staph aureus
                                        (ResistanceMechanism::MecA, drug) => {
                                            bacteria == "staphylococcus aureus" && 
                                            matches!(drug, "penicilling" | "ampicillin" | "amoxicillin" | "cephalexin" | 
                                                   "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | 
                                                   "cefepime" | "meropenem" | "imipenem_c" | "ertapenem")
                                        },
                                        // Efflux overexpression can affect multiple drug classes
                                        (ResistanceMechanism::EffluxOverexpression, _) => true,
                                        // Reduced permeability affects many drugs, especially in Gram-negatives
                                        (ResistanceMechanism::ReducedPermeability, _) => {
                                            !matches!(bacteria, "staphylococcus aureus" | "streptococcus pneumoniae" | 
                                                     "streptococcus pyogenes" | "streptococcus agalactiae" | 
                                                     "enterococcus faecalis" | "enterococcus faecium")
                                        },
                                        // Target site mutations can affect various drugs
                                        (ResistanceMechanism::TargetSiteMutation, _) => true,
                                        // AmpC affects beta-lactams
                                        (ResistanceMechanism::AmpC, drug) => {
                                            matches!(drug, "penicilling" | "ampicillin" | "amoxicillin" | "piperacillin" | 
                                                   "ticarcillin" | "cephalexin" | "cefazolin" | "cefuroxime" | 
                                                   "ceftriaxone" | "amoxicillin_clavulanate" | "piperacillin_tazobactam" |
                                                   "ampicillin_sulbactam" | "ticarcillin_clavulanate")
                                        },
                                    };
                                    
                                    if mechanism_affects_drug {
                                        let mechanism_str = mechanism.as_str();
                                        let mechanism_enhancement = get_global_param(
                                            &param_cache.resistance_mechanism_enhancement_keys[mechanism_str]
                                        ).unwrap_or(0.3);
                                        
                                        // Only add enhancement if it would actually increase resistance
                                        // Mechanisms can't decrease resistance, but they also don't add if any_r is already higher
                                        let normalized_any_r = resistance_data.any_r / max_resistance_level;
                                        if mechanism_enhancement > normalized_any_r {
                                            let additional_resistance = mechanism_enhancement - normalized_any_r;
                                            mechanism_resistance_boost += additional_resistance;
                                        }
                                    }
                                }
                            }
                        }
                        
                        // Apply mechanism enhancements to resistance levels if they would increase resistance
                        if mechanism_resistance_boost > 0.0 {
                            let normalized_any_r = resistance_data.any_r / max_resistance_level;
                            let new_resistance_level = (normalized_any_r + mechanism_resistance_boost).min(1.0);
                            let new_any_r = new_resistance_level * max_resistance_level;
                            
                            // Update any_r to the new level
                            resistance_data.any_r = new_any_r;
                            
                            // If majority_r > 0, it must equal any_r (maintain the relationship)
                            if resistance_data.majority_r > 0.0 {
                                resistance_data.majority_r = resistance_data.any_r;
                            }
                        }
                        
                        // Calculate activity_r using the updated resistance levels
                        let normalized_any_r = resistance_data.any_r / max_resistance_level;
                        resistance_data.activity_r = base_potency * drug_current_level * (1.0 - normalized_any_r);
   
                    } else {
                        resistance_data.activity_r = 0.0;
                    }
                }
            }
        }

        // testing and diagnosis - Enhanced testing framework
        let last_infected_time = individual.date_last_infected[b_idx];
        let test_delay_days = get_global_param("test_delay_days").unwrap_or(3.0) as i32;
        
        // Check if bacterial testing is available yet (historically realistic dates)
        let bacterial_testing_available_from_day = get_global_param("bacterial_testing_available_from_day").unwrap_or(5478.0) as i32;
        let bacterial_testing_available = time_step >= bacterial_testing_available_from_day as usize;
        
        if is_infected && !individual.test_identified_infection[b_idx] && last_infected_time > 0 && (time_step as i32) >= (last_infected_time + test_delay_days) && bacterial_testing_available {
            // Calculate comprehensive testing probability
            let testing_probability = calculate_testing_probability(
                individual, 
                time_step, 
                bacterial_testing_available_from_day as usize, 
                param_cache, 
                true // is_bacterial_testing
            );
            
            if rng.gen_bool(testing_probability.clamp(0.0, 1.0)) {
                individual.test_identified_infection[b_idx] = true;
            }
        }

        // --- test_r assignment logic ---
        let test_r_error_prob = get_global_param("test_r_error_probability").unwrap_or(0.02);
        let test_r_error_value = get_global_param("test_r_error_value").unwrap_or(0.25);
        let resistance_test_result_delay_days = get_global_param("resistance_test_result_delay_days").unwrap_or(2.0) as i32;
        
        // Check if resistance testing is available yet (historically realistic dates)
        let resistance_testing_available_from_day = get_global_param("resistance_testing_available_from_day").unwrap_or(9131.0) as i32;
        let resistance_testing_available = time_step >= resistance_testing_available_from_day as usize;

        if individual.test_identified_infection[b_idx] && resistance_testing_available {
            // Check if we should initiate resistance testing (if not already initiated)
            if individual.resistance_test_initiated_day[b_idx] == -1 {
                // Calculate comprehensive resistance testing probability
                let resistance_testing_probability = calculate_testing_probability(
                    individual, 
                    time_step, 
                    resistance_testing_available_from_day as usize, 
                    param_cache, 
                    false // is_bacterial_testing
                );
                
                if rng.gen_bool(resistance_testing_probability.clamp(0.0, 1.0)) {
                    // Set the flag indicating resistance testing was initiated
                    individual.test_for_resistance[b_idx] = true;
                    individual.resistance_test_initiated_day[b_idx] = time_step as i32;
                }
            }
            
            // Check if resistance test results should be available yet
            let test_initiated_day = individual.resistance_test_initiated_day[b_idx];
            if test_initiated_day != -1 && (time_step as i32) >= (test_initiated_day + resistance_test_result_delay_days) {
                let test_r_already_set = individual.resistances[b_idx].iter().any(|r| r.test_r > 0.0);
                if !test_r_already_set {
                    for d_idx in 0..DRUG_SHORT_NAMES.len() {
                        let any_r = individual.resistances[b_idx][d_idx].any_r;
                        let error = rng.gen_bool(test_r_error_prob);
                        let test_r = if error {
                            if any_r < 0.001 { test_r_error_value } else { 0.0 }
                        } else {
                            any_r
                        };
                        individual.resistances[b_idx][d_idx].test_r = test_r;
                    }
                }
            }
        } else {
            // Reset resistance test results if bacterial identification test is negative
            for d_idx in 0..DRUG_SHORT_NAMES.len() {
                individual.resistances[b_idx][d_idx].test_r = 0.0;
            }
        }

        // bacteria level change (growth/decay)
        // This entire block should only execute if the individual is currently infected with this bacteria
        if is_infected { 
            let immunity_level = individual.immune_resp[b_idx];
            let baseline_change = get_bacteria_param(bacteria, "base_bacteria_level_change").unwrap_or(0.0);
            let reduction_due_to_immune_resp = get_bacteria_param(bacteria, "immunity_effect_on_level_change").unwrap_or(0.0);
            let mut total_reduction_due_to_antibiotic = 0.0;

            // --- Resistance reversion logic: revert any_r/majority_r to 0 if not on any drug ---
            let resistance_reversion_rate = get_global_param("resistance_reversion_rate_per_day").unwrap_or(0.0001); // Default: very rare
            let on_any_drug = individual.cur_level_drug.iter().any(|&lvl| lvl > 0.0);
            if !on_any_drug {
                for drug_index in 0..DRUG_SHORT_NAMES.len() {
                    let resistance_data = &mut individual.resistances[b_idx][drug_index];
                    if resistance_data.any_r > 0.0 || resistance_data.majority_r > 0.0 {
                        if rng.gen_bool(resistance_reversion_rate) {
                            resistance_data.any_r = 0.0;
                            resistance_data.majority_r = 0.0;
                        }
                    }
                }
            }


            if individual.id == 1000001 {
                println!(" ");
                println!("mod.rs");  
                println!("bacteria: {}", bacteria);
                println!("immunity level: {:.4}", immunity_level);
                println!("baseline change: {:.4}", baseline_change);
                println!("reduction due to immune response: {:.4}", immunity_level * reduction_due_to_immune_resp);
            }


            for (drug_idx, _drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                if individual.cur_level_drug[drug_idx] > 0.0 {
                    let resistance_data = &individual.resistances[b_idx][drug_idx];
                    total_reduction_due_to_antibiotic += resistance_data.activity_r;


                if individual.id == 1000001 {
                        // Calculate standardized MIC: 1 / ((1 - majority_r) * potency)
                        let potency_param_key = &param_cache.drug_bacteria_potency_keys[&(drug_idx, b_idx)];
                        let potency = get_global_param(potency_param_key).unwrap_or(0.05);
                        let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0);
                        let normalized_majority_r = resistance_data.majority_r / max_resistance_level;
                        let standardized_mic = if (1.0 - normalized_majority_r) * potency > 0.0 {
                            1.0 / ((1.0 - normalized_majority_r) * potency)
                        } else {
                            f64::INFINITY
                        };
                        println!(
                            "mod.rs  {}: current level = {:.4}, activity_r = {:.4}, standardized_mic = {:.4}",
                            DRUG_SHORT_NAMES[drug_idx],
                            individual.cur_level_drug[drug_idx],
                            resistance_data.activity_r,
                            standardized_mic
                        );
                    }
                        
            }
            }
            
            // Apply bacteria-specific treatment response modifier to antibiotic effectiveness
            let treatment_response_modifier = get_bacteria_param(bacteria, "treatment_response_modifier").unwrap_or(1.0);
            let adjusted_antibiotic_effect = total_reduction_due_to_antibiotic * treatment_response_modifier;

             if individual.id == 1000001 {
                println!("mod.rs  total reduction due to antibiotic: {:.4}", total_reduction_due_to_antibiotic);
                println!("mod.rs  treatment response modifier: {:.4}", treatment_response_modifier);
                println!("mod.rs  adjusted antibiotic effect: {:.4}", adjusted_antibiotic_effect);
            }   
            
            let decay = baseline_change - (immunity_level * reduction_due_to_immune_resp) - adjusted_antibiotic_effect;

            let max_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);
            let new_level = (individual.level[b_idx] + decay).max(0.0).min(max_level);

            // Check for infection clearance before updating the level
            let old_level = individual.level[b_idx];
            
            if new_level < 0.0001 {
                // Check if there was an infection before clearance (previous level > 0.001)
                let was_previously_infected = old_level > 0.001;
                
                if was_previously_infected {
                    // Determine resolution type based on drugs and context
                    let has_relevant_drugs = individual.cur_use_drug.iter().enumerate()
                        .any(|(drug_idx, &on_drug)| {
                            if !on_drug { return false; }
                            // Check if this drug has potency against this bacteria
                            let potency_param_key = &param_cache.drug_bacteria_potency_keys[&(drug_idx, b_idx)];
                            let drug_potency = get_global_param(potency_param_key).unwrap_or(0.0);
                            drug_potency > 0.0
                        });
                    
                    let resolution_type = if has_relevant_drugs {
                        InfectionResolutionType::DrugAssistedClearance
                    } else {
                        InfectionResolutionType::ImmuneClearance
                    };
                    
                    let resolution_idx = match resolution_type {
                        InfectionResolutionType::ImmuneClearance => 0,
                        InfectionResolutionType::DrugAssistedClearance => 1,
                        InfectionResolutionType::DeathFromSepsis => 2,
                        InfectionResolutionType::DeathFromBackground => 3,
                        InfectionResolutionType::DeathFromToxicity => 4,
                    };
                    individual.infection_resolution_this_timestep[b_idx][resolution_idx] += 1;
                    
                    // If infection was cleared by drugs and bacteria is present in microbiome, 
                    // consider clearing it from microbiome as well
                    if matches!(resolution_type, InfectionResolutionType::DrugAssistedClearance) && 
                       individual.presence_microbiome[b_idx] {
                        let microbiome_clearance_on_drug_treatment = get_global_param("microbiome_clearance_probability_on_drug_treatment").unwrap_or(0.8);
                        if rng.gen_bool(microbiome_clearance_on_drug_treatment) {
                            individual.presence_microbiome[b_idx] = false;
                        }
                    }
                }
                
                // Clear infection data after tracking resolution
                for drug_idx_clear in 0..DRUG_SHORT_NAMES.len() {
                    let resistance_data = &mut individual.resistances[b_idx][drug_idx_clear];
                    resistance_data.any_r = 0.0;
                    resistance_data.majority_r = 0.0;
                    resistance_data.activity_r = 0.0;
                    individual.how_resistance_acquired[b_idx][drug_idx_clear] = None;
                }
                individual.level[b_idx] = 0.0;
                individual.infectious_syndrome[b_idx] = 0;
                individual.date_last_infected[b_idx] = 0;
                individual.immune_resp[b_idx] = 0.0;
                individual.sepsis[b_idx] = false;
                individual.infection_hospital_acquired[b_idx] = false;
                individual.cur_infection_from_environment[b_idx] = false;
                individual.test_identified_infection[b_idx] = false;
                individual.test_for_resistance[b_idx] = false;
                individual.resistance_test_initiated_day[b_idx] = -1;
            } else {
                // Update level for infections that are continuing
                individual.level[b_idx] = new_level;
            }
        }

        // Safety check: ensure test_identified_infection is false when not infected
        if !is_infected {
            individual.test_identified_infection[b_idx] = false;
        }

        // --- Apply cross-resistance logic ---
        apply_cross_resistance(individual, b_idx, cross_resistance_groups);
        // --- END NEW ---

        // immunity dynamics: increase during infection, decay without infection
        if is_infected {
            // immunity increase with maximum cap (only when infected)
            let infection_start_time = individual.date_last_infected[b_idx];
            let time_since_infection = (time_step as i32) - infection_start_time;
            let age = individual.age;
            let mut immune_increase = get_bacteria_param(bacteria, "immunity_base_response").unwrap_or(0.0);
            immune_increase += time_since_infection as f64 * get_bacteria_param(bacteria, "immunity_increase_per_infection_day").unwrap_or(0.0);
            immune_increase += individual.level[b_idx] * get_bacteria_param(bacteria, "immunity_increase_per_unit_higher_bacteria_level").unwrap_or(0.0);
            let age_modifier = get_bacteria_param(bacteria, "immunity_age_modifier").unwrap_or(1.0);
            immune_increase *= age_modifier.powf((age as f64 / 365.0) / 50.0);
            let immunodeficient_modifier = get_bacteria_param(bacteria, "immunity_immunodeficiency_modifier").unwrap_or(0.1);
            if individual.immunodeficiency_type.is_some() {
                immune_increase *= immunodeficient_modifier;
            }
            let max_immune_response = get_bacteria_param(bacteria, "max_immune_response").unwrap_or(10.0);
            individual.immune_resp[b_idx] = (individual.immune_resp[b_idx] + immune_increase).max(0.0001).min(max_immune_response);
        } else {
            // immunity decay when not infected
            let immunity_decay_rate = get_global_param("immune_decay_rate_per_day").unwrap_or(0.02);
            individual.immune_resp[b_idx] = (individual.immune_resp[b_idx] - immunity_decay_rate).max(0.0);
        }
    }

    // Check for post-infection drug usage evaluation (configurable timing)
    let evaluation_days = get_global_param("drug_evaluation_days_post_infection").unwrap_or(7.0) as i32;
    
    for b_idx in 0..BACTERIA_LIST.len() {
        let infection_start_day = individual.date_last_infected_keep[b_idx];
        
        // Only evaluate if there was an infection and today is exactly the evaluation day after infection start
        if infection_start_day > 0 && (time_step as i32) == (infection_start_day + evaluation_days) {
            // Check if any drug was initiated since the infection started
            let mut drug_used_since_infection = false;
            
            for d_idx in 0..DRUG_SHORT_NAMES.len() {
                let drug_start_day = individual.date_drug_initiated_keep[d_idx];
                
                // Drug was started if it was initiated on or after the infection start day
                if drug_start_day != i32::MIN && drug_start_day >= infection_start_day {
                    drug_used_since_infection = true;
                    break;
                }
            }
            
            // Set the evaluation result for this bacteria (this will be counted once in summary stats)
            individual.day_7_since_last_infection_drug_used[b_idx] = Some(drug_used_since_infection);
        }
    }
    
    // Note: We do NOT reset day_7_since_last_infection_drug_used values here because 
    // the summary statistics need to capture them during this timestep. 
    // They will be reset when a new infection occurs or when the infection clears.
    
    // Update the current number of drugs counter at the end of each timestep
    update_drug_counter(individual);
}

/// New helper function to apply cross-resistance within drug groups for a specific bacteria.
fn apply_cross_resistance(
    individual: &mut Individual,
    b_idx: usize,
    cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>,
) {
    // Check if there are any cross-resistance groups defined for this bacterium
    if let Some(groups) = cross_resistance_groups.get(&b_idx) {
        for group in groups {
            // Find the maximum any_r value in the current group
            let mut max_any_r = 0.0;
            for &d_idx in group {
                if let Some(resistance_data) = individual.resistances.get(b_idx).and_then(|r| r.get(d_idx)) {
                    if resistance_data.any_r > max_any_r {
                        max_any_r = resistance_data.any_r;
                    }
                }
            }

            // If there's any resistance in the group, update all drugs in the group to the max value
            if max_any_r > 0.0 {
                for &d_idx in group {
                    if let Some(resistance_data) = individual.resistances.get_mut(b_idx).and_then(|r| r.get_mut(d_idx)) {
                        resistance_data.any_r = max_any_r;
                    }
                }
            }
        }
    }
}

/// Calculate comprehensive testing probability based on multiple factors
fn calculate_testing_probability(
    individual: &Individual, 
    time_step: usize, 
    testing_available_from_day: usize, 
    param_cache: &ParameterKeyCache,
    is_bacterial_testing: bool
) -> f64 {
    // Get base parameters
    let base_rate = if is_bacterial_testing {
        get_global_param("bacterial_testing_base_rate_per_day").unwrap_or(0.15)
    } else {
        get_global_param("resistance_testing_base_rate_per_day").unwrap_or(0.95)
    };
    
    // Calculate temporal multiplier (testing adoption over time)
    let years_since_availability = (time_step - testing_available_from_day) as f64 / 365.0;
    let (initial_rate, adoption_rate, max_multiplier) = if is_bacterial_testing {
        (
            get_global_param("bacterial_testing_initial_adoption_rate").unwrap_or(0.1),
            get_global_param("bacterial_testing_adoption_rate_per_year").unwrap_or(0.025),
            get_global_param("bacterial_testing_max_temporal_multiplier").unwrap_or(1.0)
        )
    } else {
        (
            get_global_param("resistance_testing_initial_adoption_rate").unwrap_or(0.05),
            get_global_param("resistance_testing_adoption_rate_per_year").unwrap_or(0.015),
            get_global_param("resistance_testing_max_temporal_multiplier").unwrap_or(1.0)
        )
    };
    
    let temporal_multiplier = (initial_rate + years_since_availability * adoption_rate).min(max_multiplier);
    
    // Hospital status multiplier
    let hospital_multiplier = if individual.hospital_status.is_hospitalized() {
        if is_bacterial_testing {
            get_global_param("bacterial_testing_hospital_multiplier").unwrap_or(8.0)
        } else {
            get_global_param("resistance_testing_hospital_multiplier").unwrap_or(5.0)
        }
    } else {
        1.0
    };
    
    // Regional resource multiplier
    let region_name = individual.region_cur_in.to_string().to_lowercase().replace(" ", "_");
    let region_multiplier = if let Some(key) = param_cache.region_testing_keys.get(&region_name) {
        get_global_param(key).unwrap_or(1.0)
    } else {
        1.0 // Default if region not found
    };
    
    // Immunosuppression multiplier
    let immunosuppression_multiplier = if individual.immunodeficiency_type.is_some() {
        get_global_param("testing_immunosuppressed_multiplier").unwrap_or(2.5)
    } else {
        1.0
    };
    
    // Sepsis multiplier
    let sepsis_multiplier = if individual.sepsis.iter().any(|&s| s) {
        get_global_param("testing_sepsis_multiplier").unwrap_or(4.0)
    } else {
        1.0
    };
    
    // Calculate final probability
    let final_probability = base_rate * temporal_multiplier * hospital_multiplier * region_multiplier * immunosuppression_multiplier * sepsis_multiplier;
    
    // Cap at 1.0 (100% probability)
    final_probability.min(1.0)
}

/// Helper function to probabilistically assign a syndrome for a given bacteria.
fn assign_syndrome_for_bacteria<R: Rng>(bacteria: &str, rng: &mut R) -> u32 {
    // Define syndrome probabilities for each bacteria based on clinical epidemiology.
    // Each entry: (syndrome_id, probability)
    // Syndromes: 1=UTI, 2=Skin/soft tissue, 3=Respiratory, 4=Bloodstream, 5=Intra-abdominal, 
    //           6=CNS, 7=GI, 8=Genital, 9=Bone/joint, 10=Other
    let syndrome_probs: &[(u32, f64)] = match bacteria {
        // Gram-positive cocci
        "staphylococcus aureus" => &[(2, 0.35), (4, 0.25), (9, 0.15), (3, 0.10), (5, 0.08), (1, 0.05), (6, 0.02)],
        "streptococcus pneumoniae" => &[(3, 0.70), (6, 0.15), (4, 0.08), (1, 0.04), (2, 0.02), (10, 0.01)],
        "streptococcus pyogenes" => &[(2, 0.50), (3, 0.25), (4, 0.15), (9, 0.05), (5, 0.03), (1, 0.02)],
        "streptococcus agalactiae" => &[(4, 0.40), (6, 0.25), (1, 0.15), (2, 0.10), (3, 0.05), (5, 0.05)],
        "enterococcus faecalis" => &[(1, 0.50), (4, 0.25), (5, 0.15), (2, 0.05), (3, 0.03), (9, 0.02)],
        "enterococcus faecium" => &[(1, 0.45), (4, 0.30), (5, 0.15), (2, 0.05), (3, 0.03), (9, 0.02)],
        
        // Gram-negative Enterobacteriaceae
        "escherichia coli" => &[(1, 0.55), (4, 0.20), (5, 0.12), (7, 0.08), (2, 0.03), (3, 0.02)],
        "klebsiella pneumoniae" => &[(3, 0.40), (1, 0.25), (4, 0.20), (5, 0.10), (2, 0.03), (7, 0.02)],
        "enterobacter spp." => &[(1, 0.35), (3, 0.25), (4, 0.20), (5, 0.10), (7, 0.05), (2, 0.05)],
        "enterobacter_cloacae" => &[(4, 0.30), (3, 0.25), (1, 0.25), (5, 0.12), (7, 0.05), (2, 0.03)],
        "citrobacter spp." => &[(1, 0.30), (3, 0.25), (4, 0.20), (5, 0.15), (7, 0.05), (2, 0.05)],
        "serratia spp." => &[(3, 0.35), (1, 0.25), (4, 0.20), (5, 0.10), (2, 0.05), (7, 0.05)],
        "proteus spp." => &[(1, 0.60), (4, 0.15), (3, 0.10), (5, 0.08), (2, 0.04), (7, 0.03)],
        "morganella spp." => &[(1, 0.50), (4, 0.20), (3, 0.15), (5, 0.08), (2, 0.04), (7, 0.03)],
        
        // Non-fermenting Gram-negatives
        "pseudomonas aeruginosa" => &[(3, 0.45), (4, 0.25), (1, 0.15), (2, 0.08), (5, 0.05), (9, 0.02)],
        "acinetobacter baumannii" => &[(3, 0.40), (4, 0.25), (1, 0.15), (5, 0.10), (2, 0.05), (7, 0.05)],
        
        // Gastrointestinal pathogens
        "salmonella enterica serovar typhi" => &[(7, 0.80), (4, 0.15), (5, 0.03), (3, 0.01), (10, 0.01)],
        "salmonella enterica serovar paratyphi a" => &[(7, 0.85), (4, 0.10), (5, 0.03), (3, 0.01), (10, 0.01)],
        "invasive non-typhoidal salmonella spp." => &[(7, 0.70), (4, 0.20), (5, 0.05), (3, 0.03), (1, 0.02)],
        "shigella spp." => &[(7, 0.95), (4, 0.03), (5, 0.01), (10, 0.01)],
        "vibrio cholerae" => &[(7, 0.98), (5, 0.01), (10, 0.01)],
        "campylobacter_jejuni" => &[(7, 0.70), (9, 0.15), (4, 0.08), (5, 0.05), (3, 0.02)],
        "yersinia_enterocolitica" => &[(7, 0.75), (5, 0.15), (4, 0.05), (9, 0.03), (3, 0.02)],
        "clostridioides_difficile" => &[(7, 0.90), (5, 0.08), (4, 0.01), (10, 0.01)],
        
        // Sexually transmitted pathogens
        "neisseria gonorrhoeae" => &[(8, 0.85), (1, 0.10), (5, 0.03), (4, 0.01), (10, 0.01)],
        "chlamydia trachomatis" => &[(8, 0.70), (1, 0.20), (5, 0.05), (6, 0.03), (10, 0.02)],
        "treponema pallidum" => &[(8, 0.60), (2, 0.20), (6, 0.10), (4, 0.05), (10, 0.05)],
        
        // Respiratory pathogens
        "haemophilus influenzae" => &[(3, 0.70), (6, 0.15), (4, 0.08), (1, 0.04), (2, 0.02), (10, 0.01)],
        "moraxella_catarrhalis" => &[(3, 0.85), (4, 0.08), (1, 0.04), (2, 0.02), (10, 0.01)],
        "neisseria_meningitidis" => &[(6, 0.60), (4, 0.25), (3, 0.10), (2, 0.03), (1, 0.02)],
        "bordetella pertussis" => &[(3, 0.95), (6, 0.03), (4, 0.01), (10, 0.01)], // Primarily respiratory (whooping cough)
        
        // Gastrointestinal pathogens
        "helicobacter pylori" => &[(7, 0.85), (5, 0.10), (4, 0.03), (10, 0.02)], // Primarily GI (peptic ulcer disease)
        
        // Foodborne/systemic pathogens  
        "listeria_monocytogenes" => &[(6, 0.50), (4, 0.30), (7, 0.10), (5, 0.05), (3, 0.03), (1, 0.02)],
        
        // Fallback for any unmatched bacteria (should not occur with complete list above)
        _ => &[(1, 0.1), (2, 0.1), (3, 0.1), (4, 0.1), (5, 0.1), (6, 0.1), (7, 0.1), (8, 0.1), (9, 0.1), (10, 0.1)],
    };

    let weights: Vec<f64> = syndrome_probs.iter().map(|&(_, p)| p).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    syndrome_probs[dist.sample(rng)].0
}