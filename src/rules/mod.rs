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

use crate::config::{
    get_age_dependent_bacteria_sepsis_risk_multiplier, get_bacteria_sepsis_risk_multiplier,
    get_drug_availability_time_aware, get_drug_introduction_time_step, get_global_param,
    parameter_store,
};
use crate::simulation::population::{
    HospitalStatus, ImmunodeficiencyType, Individual, InfectionResolutionType, Region,
    BACTERIA_LIST, DRUG_SHORT_NAMES, MICROBIOME_MAJORITY_THRESHOLD,
};
use rand::seq::SliceRandom;
use rand::Rng;

use crate::simulation::simulation::MajorityRCache;
use std::collections::HashMap;

/// Helper function to update the current number of drugs counter
fn update_drug_counter(individual: &mut Individual) {
    individual.current_number_of_drugs =
        individual.cur_use_drug.iter().filter(|&&on| on).count() as i32;
}

/// Apply pairwise drug level interactions based on pharmacokinetic effects
/// Modifies individual.cur_level_drug in-place to account for drug-drug interactions
fn apply_drug_level_interactions(individual: &mut Individual) {
    let store = parameter_store();
    // Create a copy of current levels to calculate interactions from baseline levels
    let original_levels = individual.cur_level_drug.clone();

    // Identify which drugs have significant levels (>0.001, roughly 0.1% of standard dose)
    let active_drugs: Vec<usize> = original_levels
        .iter()
        .enumerate()
        .filter(|(_, &level)| level > 0.001)
        .map(|(idx, _)| idx)
        .collect();

    // If fewer than 2 drugs active, no interactions possible
    if active_drugs.len() < 2 {
        return;
    }

    // Apply each pairwise interaction
    for &drug1_idx in &active_drugs {
        let drug1_name = DRUG_SHORT_NAMES[drug1_idx];

        for &drug2_idx in &active_drugs {
            if drug1_idx == drug2_idx {
                continue; // Skip self-interactions
            }

            let drug2_name = DRUG_SHORT_NAMES[drug2_idx];

            let interaction_key = format!(
                "drug_level_multiplier_{}_when_coadministered_with_{}",
                drug1_name, drug2_name
            );

            if let Some(multiplier) = get_global_param(&interaction_key) {
                // Apply the interaction multiplier to drug1's level
                // Only apply if it would actually change the level (avoid redundant 1.0 multipliers)
                if (multiplier - 1.0).abs() > 0.001 {
                    individual.cur_level_drug[drug1_idx] *= multiplier;

                    // Ensure levels don't go negative or below detection threshold
                    if individual.cur_level_drug[drug1_idx] < 0.001 {
                        individual.cur_level_drug[drug1_idx] = 0.0;
                    }

                    // Cap levels at reasonable maximum (e.g., 5x standard dose to prevent unrealistic accumulation)
                    let max_level = store.drug.initial_level(drug1_idx) * 5.0;
                    if individual.cur_level_drug[drug1_idx] > max_level {
                        individual.cur_level_drug[drug1_idx] = max_level;
                    }
                }
            }
        }
    }
}
use rand::distributions::Distribution;
use rand::distributions::WeightedIndex;

/// Assess treatment failure and switch drugs if necessary
/// Returns true if a drug switch occurred
fn assess_treatment_failure(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    bacteria_indices: &HashMap<&'static str, usize>,
    _drug_indices: &HashMap<&'static str, usize>,
    _cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>,
    _param_cache: &ParameterKeyCache,
    rng: &mut impl Rng,
) -> bool {
    let store = parameter_store();

    // Check if treatment failure assessment is enabled
    if !store.globals.treatment_failure_enabled {
        return false;
    }

    let bacteria_name = BACTERIA_LIST[bacteria_idx];
    let syndrome_id = individual.infectious_syndrome[bacteria_idx];
    let base_assessment_day = store.globals.treatment_failure_assessment_day;
    let assessment_day =
        treatment_failure_assessment_day_for(bacteria_name, syndrome_id, base_assessment_day);

    // Check if we've reached the assessment window for this organism/syndrome
    if individual.days_on_current_treatment[bacteria_idx] < assessment_day {
        return false;
    }

    // Check if we've already assessed this treatment course
    if individual.treatment_failure_assessed[bacteria_idx] {
        return false;
    }

    // Check if there's a current infection and bacteria level recorded at drug start
    if individual.level[bacteria_idx] <= 0.0
        || individual.bacteria_level_at_drug_start[bacteria_idx].is_none()
    {
        return false;
    }

    let initial_level = individual.bacteria_level_at_drug_start[bacteria_idx].unwrap();
    let current_level = individual.level[bacteria_idx];

    // Get failure threshold (default 0.5 = 50% of initial level)
    let threshold_level = initial_level * store.globals.treatment_failure_threshold;

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
    let current_drugs: Vec<usize> = individual
        .cur_use_drug
        .iter()
        .enumerate()
        .filter(|(_, &is_taking)| is_taking)
        .map(|(drug_idx, _)| drug_idx)
        .collect();

    if current_drugs.is_empty() {
        return false; // No current drugs to switch from
    }

    // Try to find an alternative drug using the same selection logic as initial prescription
    // but excluding recently failed drugs
    let failure_memory_days = store.globals.drug_failure_memory_days;

    // Build list of available alternative drugs
    let mut alternative_scores = Vec::new();

    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        // Skip if currently taking this drug
        if current_drugs.contains(&drug_idx) {
            continue;
        }

        // Skip if this drug failed recently (within memory period)
        if individual.date_drug_initiated_keep[drug_idx] != i32::MIN {
            let days_since_last_use =
                (time_step as i32) - individual.date_drug_initiated_keep[drug_idx];
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
            time_step,
        );

        // Check if drug has been historically introduced (CRITICAL: was missing!)
        let intro_ok = match get_drug_introduction_time_step(drug_name) {
            Some(intro_step) => time_step >= intro_step,
            None => true,
        };

        if avail < 0.01 || !intro_ok {
            // Drug not sufficiently available OR not yet introduced
            continue;
        }

        // Calculate drug score using same logic as original selection
        let mut score = 0.0;

        // Base potency score
        let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
        let potency = store
            .drug_bacteria
            .potency(*bacteria_idx_for_cache, drug_idx);
        if potency >= store.globals.minimal_potency_threshold_for_drug_selection {
            score += potency;
        }

        // Apply clinical multipliers (same as original logic)
        // Add pathogen-specific preference multipliers
        let bacteria_drug_key = format!(
            "{}_{}_clinical_preference_multiplier",
            bacteria_name.replace(" ", "_"),
            drug_name
        );
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
        let selection_temperature = store.globals.drug_selection_temperature;
        let weights: Vec<f64> = alternative_scores
            .iter()
            .map(|(_, score)| (score / selection_temperature).exp())
            .collect();

        let total_weight: f64 = weights.iter().sum();
        if total_weight > 0.0 && total_weight.is_finite() {
            let dist = WeightedIndex::new(&weights).unwrap();
            let chosen_idx = dist.sample(rng);
            let new_drug_idx = alternative_scores[chosen_idx].0;

            // Stop current drugs
            for &current_drug_idx in &current_drugs {
                individual.cur_use_drug[current_drug_idx] = false;
                individual.date_drug_initiated[current_drug_idx] = i32::MIN;
            }

            // Start new drug
            individual.cur_use_drug[new_drug_idx] = true;
            individual.date_drug_initiated[new_drug_idx] = time_step as i32;
            individual.date_drug_initiated_keep[new_drug_idx] = time_step as i32;
            individual.ever_taken_drug[new_drug_idx] = true;

            // Update drug counter
            update_drug_counter(individual);

            // Set drug level
            let initial_level = store.drug.initial_level(new_drug_idx);
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

fn treatment_failure_assessment_day_for(
    bacteria_name: &str,
    syndrome_id: i32,
    default_day: i32,
) -> i32 {
    let mut final_day = default_day.max(1);

    // Rapid infection syndromes: respiratory (3), bloodstream (4), intra-abdominal (5), CNS (6)
    let fast_track_syndromes = [3, 4, 5, 6];
    if fast_track_syndromes.contains(&syndrome_id) {
        final_day = final_day.min(3).max(2);
    }

    // Chronic or slow pathogens: TB and indolent infections get longer assessment windows
    if bacteria_name == "mdr mycobacterium tuberculosis" {
        final_day = final_day.max(10);
    } else if bacteria_name == "helicobacter pylori" || syndrome_id == 9 {
        final_day = final_day.max(6);
    }

    final_day
}

/// Assess restart window for patients who stopped drugs while still infected
/// Returns true if restart treatment was initiated
fn assess_restart_window(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    bacteria_indices: &HashMap<&'static str, usize>,
    param_cache: &ParameterKeyCache,
    rng: &mut impl Rng,
) -> bool {
    let store = parameter_store();

    // Check if restart window is enabled
    if !store.globals.restart_window_enabled {
        return false;
    }

    // Check if there's a cessation to assess
    if let Some(cessation_day) = individual.drug_stopped_with_infection_day[bacteria_idx] {
        let restart_window_days = store.globals.restart_window_days;
        let days_since_cessation = (time_step as i32) - cessation_day;

        // Within restart window?
        if days_since_cessation >= 1 && days_since_cessation <= restart_window_days {
            // Haven't assessed yet?
            if !individual.restart_window_assessed[bacteria_idx] {
                individual.restart_window_assessed[bacteria_idx] = true;

                // Check if bacteria level has worsened enough to trigger restart
                if let Some(cessation_level) =
                    individual.bacteria_level_at_drug_cessation[bacteria_idx]
                {
                    let current_level = individual.level[bacteria_idx];
                    let threshold_multiplier = store.globals.restart_bacteria_level_threshold;

                    // Restart criteria: bacteria level increased significantly OR still very high
                    let bacteria_worsened =
                        current_level >= (cessation_level * threshold_multiplier);
                    let bacteria_still_high = current_level > 2.0; // Arbitrary high threshold for severe infection

                    if (bacteria_worsened || bacteria_still_high)
                        && individual.level[bacteria_idx] > 0.1
                    {
                        // Patient decides to return to care?
                        let return_probability = store.globals.restart_window_probability;

                        if rng.gen_bool(return_probability) {
                            // Clear restart tracking
                            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
                            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
                            let stopped_drug_idx = individual.stopped_drug_index[bacteria_idx];
                            individual.stopped_drug_index[bacteria_idx] = None;

                            // Start restart treatment, preferring the previously effective drug
                            return start_restart_treatment(
                                individual,
                                time_step,
                                bacteria_idx,
                                stopped_drug_idx,
                                bacteria_indices,
                                param_cache,
                                rng,
                            );
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
    _param_cache: &ParameterKeyCache,
    rng: &mut impl Rng,
) -> bool {
    let store = parameter_store();

    let bacteria_name = BACTERIA_LIST[bacteria_idx];
    let minimal_potency_threshold = store.globals.minimal_potency_threshold_for_drug_selection;

    // Check if we can restart the previously effective drug
    if let Some(prev_drug_idx) = stopped_drug_idx {
        let prev_drug_name = DRUG_SHORT_NAMES[prev_drug_idx];

        // Check if previously effective drug is still available
        let avail = get_drug_availability_time_aware(
            prev_drug_name,
            &individual.region_cur_in.to_string(),
            Some(&individual.region_living.to_string()),
            time_step,
        );

        // Check if drug has been historically introduced
        let intro_ok = match get_drug_introduction_time_step(prev_drug_name) {
            Some(intro_step) => time_step >= intro_step,
            None => true,
        };

        if avail >= 0.01 && intro_ok && !individual.cur_use_drug[prev_drug_idx] {
            // Check if drug has adequate potency (basic safety check)
            let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
            let potency = store
                .drug_bacteria
                .potency(*bacteria_idx_for_cache, prev_drug_idx);
            if potency >= minimal_potency_threshold {
                // Restart the previously effective drug!
                individual.cur_use_drug[prev_drug_idx] = true;
                individual.date_drug_initiated[prev_drug_idx] = time_step as i32;
                individual.date_drug_initiated_keep[prev_drug_idx] = time_step as i32;
                individual.ever_taken_drug[prev_drug_idx] = true;

                // Update drug counter
                update_drug_counter(individual);

                // Set drug level
                let initial_level = store.drug.initial_level(prev_drug_idx);
                individual.cur_level_drug[prev_drug_idx] = initial_level;

                // Reset treatment failure tracking for new treatment
                individual.bacteria_level_at_drug_start[bacteria_idx] =
                    Some(individual.level[bacteria_idx]);
                individual.days_on_current_treatment[bacteria_idx] = 0;
                individual.treatment_failure_assessed[bacteria_idx] = false;

                return true; // Successfully restarted previously effective drug
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
            time_step,
        );

        // Check if drug has been historically introduced
        let intro_ok = match get_drug_introduction_time_step(drug_name) {
            Some(intro_step) => time_step >= intro_step,
            None => true,
        };

        if avail < 0.01 || !intro_ok {
            continue;
        }

        // Calculate drug score
        let mut score = 0.0;

        // Base potency score
        let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
        let potency = store
            .drug_bacteria
            .potency(*bacteria_idx_for_cache, drug_idx);
        if potency >= minimal_potency_threshold {
            score += potency;
        }

        // Apply clinical preference multipliers
        let bacteria_drug_key = format!(
            "{}_{}_clinical_preference_multiplier",
            bacteria_name.replace(" ", "_"),
            drug_name
        );
        if let Some(preference_multiplier) = get_global_param(&bacteria_drug_key) {
            score *= preference_multiplier;
        }

        // BONUS: If this was the previously effective drug, give it preference
        if let Some(prev_drug_idx) = stopped_drug_idx {
            if drug_idx == prev_drug_idx {
                let effectiveness_bonus = store.globals.previously_effective_drug_bonus;
                score *= effectiveness_bonus;
            }
        }

        if score > 0.0 {
            drug_scores.push((drug_idx, score));
        }
    }

    // Select and start restart treatment
    if !drug_scores.is_empty() {
        let selection_temperature = store.globals.drug_selection_temperature;
        let weights: Vec<f64> = drug_scores
            .iter()
            .map(|(_, score)| (score / selection_temperature).exp())
            .collect();

        let total_weight: f64 = weights.iter().sum();
        if total_weight > 0.0 && total_weight.is_finite() {
            let dist = WeightedIndex::new(&weights).unwrap();
            let chosen_idx = dist.sample(rng);
            let new_drug_idx = drug_scores[chosen_idx].0;

            // Start restart treatment
            individual.cur_use_drug[new_drug_idx] = true;
            individual.date_drug_initiated[new_drug_idx] = time_step as i32;
            individual.date_drug_initiated_keep[new_drug_idx] = time_step as i32;
            individual.ever_taken_drug[new_drug_idx] = true;

            // Update drug counter
            update_drug_counter(individual);

            // Set drug level
            let initial_level = store.drug.initial_level(new_drug_idx);
            individual.cur_level_drug[new_drug_idx] = initial_level;

            // Reset treatment failure tracking for new treatment
            individual.bacteria_level_at_drug_start[bacteria_idx] =
                Some(individual.level[bacteria_idx]);
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
    // Region-based keys

    // Other frequently used keys
}

impl ParameterKeyCache {
    pub fn new() -> Self {
        let mut cache = ParameterKeyCache {
            drug_bacteria_potency_keys: HashMap::new(),
        };

        // Pre-compute all drug/bacteria combinations
        for (d_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                cache.drug_bacteria_potency_keys.insert(
                    (d_idx, b_idx),
                    format!(
                        "drug_{}_for_bacteria_{}_potency_when_no_r",
                        drug_name, bacteria_name
                    ),
                );
            }
        }

        cache
    }
}

/// applies model rules to an individual for one time step.
pub fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    rng: &mut impl Rng,
    majority_r_cache: &MajorityRCache,
    bacteria_indices: &HashMap<&'static str, usize>,
    drug_indices: &HashMap<&'static str, usize>,
    cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>, // New parameter
    param_cache: &ParameterKeyCache,                           // New parameter cache
) {
    let store = parameter_store();

    if individual.age < 0 {
        individual.age += 1; // Only advance age by 1 day
        return; // Exit the function if unborn
    }

    if individual.date_of_death.is_some() {
        return; // Exit the function if dead
    }

    // Reset microbiome acquisition flags ahead of this timestep's updates
    for flag in &mut individual.microbiome_acquired_today {
        *flag = false;
    }
    for flag in &mut individual.microbiome_acquired_on_drug_today {
        *flag = false;
    }
    for flag in &mut individual.microbiome_cleared_today {
        *flag = false;
    }

    // --- all these parameter lookups at the top so they're in scope everywhere ---
    let transfer_prob = store
        .globals
        .microbiome_resistance_transfer_probability_per_day;
    let drug_base_initiation_rate = store.globals.drug_base_initiation_rate_per_day;
    let drug_infection_present_multiplier = store.globals.drug_infection_present_multiplier;
    let already_on_drug_initiation_multiplier = store.globals.already_on_drug_initiation_multiplier;
    let drug_test_identified_multiplier = store.globals.drug_test_identified_multiplier;
    let double_dose_probability = store
        .globals
        .double_dose_probability_if_identified_infection;
    let random_drug_cessation_prob = store.globals.random_drug_cessation_probability;

    // update non-infection, bacteria or antibiotic-specific variables
    // need a variable for vulnerability to serious toxicity ?
    individual.age += 1;

    // ---  Update Contact and Exposure Levels ---
    //  update immunodeficiency status based on onset/recovery rates and type

    let immunodeficiency_params = &store.immunodeficiency;

    // Get rates for both types
    let temp_onset_rate = immunodeficiency_params.temporary_onset_rate();
    let temp_recovery_rate = immunodeficiency_params.temporary_recovery_rate();
    let chronic_onset_rate = immunodeficiency_params.chronic_onset_rate();
    let chronic_recovery_rate = immunodeficiency_params.chronic_recovery_rate();

    // Get age-based probability for chronic vs temporary assignment
    let chronic_probability = immunodeficiency_params.chronic_probability(individual.age);

    match individual.immunodeficiency_type {
        Some(ImmunodeficiencyType::Temporary) => {
            // Currently has temporary immunodeficiency, check for recovery
            if rng.gen_bool(temp_recovery_rate) {
                individual.immunodeficiency_type = None;
            }
        }
        Some(ImmunodeficiencyType::Chronic) => {
            // Currently has chronic immunodeficiency, check for recovery
            if rng.gen_bool(chronic_recovery_rate) {
                individual.immunodeficiency_type = None;
            }
        }
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

    // Get parameters from config.rs once per individual for this time step
    let baseline_rate = store.globals.hospital_baseline_rate_per_day;
    let age_multiplier_hosp = store.globals.hospital_age_multiplier_per_day;
    let recovery_rate = store.globals.hospital_recovery_rate_per_day;
    let max_days_in_hospital = store.globals.hospital_max_days.max(0.0) as u32;
    let sepsis_admission_multiplier = store.globals.hospital_sepsis_admission_multiplier;
    let prevent_discharge_with_sepsis = store.globals.hospital_prevent_discharge_with_sepsis > 0.5;

    // Check if individual has any active sepsis
    let has_sepsis = individual.sepsis.iter().any(|&s| s);

    // Potentially get hospitalized (if not currently hospitalized)
    if !individual.hospital_status.is_hospitalized() {
        let mut prob_hospitalization_today =
            baseline_rate + (individual.age as f64 * age_multiplier_hosp);

        // Strong sepsis admission effect - sepsis patients are very likely to be hospitalized
        if has_sepsis {
            prob_hospitalization_today *= sepsis_admission_multiplier;
        }

        if rng.gen::<f64>() < prob_hospitalization_today {
            individual.hospital_status = HospitalStatus::InHospital;
            individual.days_hospitalized = 0; // Initialize days hospitalized
        }
    } else {
        // If already hospitalized, consider recovery or max days limit
        individual.days_hospitalized += 1; // Increment days hospitalized

        // Determine if discharge is allowed
        let can_discharge = if prevent_discharge_with_sepsis {
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
        else if can_discharge && individual.days_hospitalized >= max_days_in_hospital {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
        }
    }
    // --- end hospitalization Rules ---

    // ---  region travel ---
    let base_travel_prob = store.globals.travel_probability_per_day;

    // Apply region-specific travel multiplier based on individual's home region
    let travel_prob = base_travel_prob * store.region.travel_multiplier(individual.region_living);

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
                    }
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
                    }
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
                    }
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
                    }
                    Region::Home => {
                        // Should not reach here, but default to global uniform if it does
                        vec![
                            (Region::Asia, 0.167),
                            (Region::Africa, 0.167),
                            (Region::Europe, 0.166),
                            (Region::NorthAmerica, 0.167),
                            (Region::SouthAmerica, 0.166),
                            (Region::Oceania, 0.167),
                        ]
                    }
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
                let sepsis_baseline_log_odds = store.bacteria.sepsis_baseline_log_odds(b_idx);
                let log_odds_infection_level =
                    store.bacteria.sepsis_log_odds_infection_level(b_idx);
                let log_odds_infection_duration =
                    store.bacteria.sepsis_log_odds_infection_duration(b_idx);

                // ENHANCED BACTERIA SEPSIS RISK CALCULATION
                // Combines: 1) Enhanced bacteria-specific risk, 2) Age-dependent interactions, 3) Clinical risk categories
                let bacteria_sepsis_risk = get_bacteria_sepsis_risk_multiplier(bacteria);
                let age_bacteria_sepsis_risk = get_age_dependent_bacteria_sepsis_risk_multiplier(
                    bacteria,
                    individual.age as u32,
                );

                // Combined bacteria risk multiplier (bacteria-specific × age-dependent interaction)
                let combined_bacteria_risk = bacteria_sepsis_risk * age_bacteria_sepsis_risk;

                // Map combined risk to log odds categories for logistic regression
                let bacteria_log_odds = if combined_bacteria_risk >= 3.0 {
                    // Very high combined risk (e.g., MRSA in elderly, GBS in neonates)
                    store.globals.log_odds_bacteria_with_high_sepsis_risk * 1.5
                } else if combined_bacteria_risk >= 1.8 {
                    // High combined risk
                    store.globals.log_odds_bacteria_with_high_sepsis_risk
                } else if (0.7..=1.3).contains(&combined_bacteria_risk) {
                    // Medium combined risk (reference category)
                    store.globals.log_odds_bacteria_with_medium_sepsis_risk
                } else if combined_bacteria_risk >= 0.3 {
                    // Low combined risk
                    store.globals.log_odds_bacteria_with_low_sepsis_risk
                } else {
                    // Very low combined risk (e.g., Chlamydia, localized infections)
                    store.globals.log_odds_bacteria_with_low_sepsis_risk * 0.5
                };

                // Add syndrome-specific sepsis risk (infection site effect)
                // This allows the same bacteria to have different sepsis risks depending on infection site
                // e.g., E. coli UTI vs E. coli bacteremia have very different sepsis risks
                let syndrome_log_odds = if individual.infectious_syndrome[b_idx] > 0 {
                    store
                        .syndrome
                        .sepsis_log_odds(individual.infectious_syndrome[b_idx] as usize)
                } else {
                    0.0 // No syndrome specified, no effect
                };

                // Add regional sepsis risk factors (healthcare access, population density, resources)
                let region_log_odds = match individual.region_living {
                    Region::Africa | Region::Asia => store.globals.log_odds_sepsis_region_b, // Lower resource regions
                    Region::NorthAmerica | Region::Europe | Region::Oceania => {
                        store.globals.log_odds_sepsis_region_a
                    } // Higher resource regions
                    Region::SouthAmerica => store.globals.log_odds_sepsis_region_b, // Mixed resource region
                    Region::Home => 0.0, // Neutral/no effect for home region
                };

                // COMPREHENSIVE SEPSIS RISK CALCULATION
                // Integrates: bacteria risk, age interactions, syndrome site, regional factors
                let log_odds_sepsis = sepsis_baseline_log_odds
                    + (current_level * log_odds_infection_level)
                    + (duration_of_infection as f64 * log_odds_infection_duration)
                    + bacteria_log_odds
                    + syndrome_log_odds
                    + region_log_odds;

                // EXPLICIT H. PYLORI SEPSIS PREVENTION
                // If H. pylori is the only infection, force sepsis risk to zero
                let prob_sepsis_today = if bacteria == "helicobacter pylori" {
                    // Check if this is the only active infection
                    let other_infections_exist = individual
                        .level
                        .iter()
                        .enumerate()
                        .any(|(idx, &level)| idx != b_idx && level > 0.001);

                    if !other_infections_exist {
                        // H. pylori as sole infection = ZERO sepsis risk
                        // Also clear any existing sepsis status from H. pylori
                        if individual.sepsis[b_idx] {
                            individual.sepsis[b_idx] = false;
                        }
                        0.0
                    } else {
                        // H. pylori + other bacteria = use calculated risk
                        1.0 / (1.0 + (-log_odds_sepsis).exp())
                    }
                } else {
                    // Non-H. pylori bacteria = use calculated risk
                    1.0 / (1.0 + (-log_odds_sepsis).exp())
                };

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
    let age_idx = crate::config::VaccinationParameters::age_group_index(age_years);

    // Calculate simulation year (assuming time_step 0 = year 1950, one step per day)
    let simulation_year = 1950.0 + (time_step as f64 / 365.0);

    const BACTERIAL_VACCINES: [&str; 3] = ["pneumococcal", "meningococcal", "hib"];
    for (b_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
        // For each bacterial vaccine, check if this bacteria is targeted by the vaccine
        for &vaccine in &BACTERIAL_VACCINES {
            if let Some(vaccine_idx) = crate::config::VaccinationParameters::vaccine_index(vaccine)
            {
                let availability_year = store.vaccination.availability_year(vaccine_idx);
                if simulation_year < availability_year {
                    continue; // Vaccine not yet available
                }

                // Correct bacteria name matching (fixing underscore vs space issues)
                let targets_bacteria = match (vaccine, *bacteria) {
                    ("pneumococcal", "streptococcus pneumoniae") => true,
                    ("meningococcal", "neisseria_meningitidis") => true, // Fixed: using underscore version
                    ("hib", "haemophilus influenzae") => true,
                    ("pertussis", "bordetella pertussis") => true, // DTaP/Tdap vaccines
                    _ => false,
                };

                if targets_bacteria && !individual.vaccination_status[b_idx] {
                    let daily_prob = store.vaccination.daily_probability(vaccine_idx, age_idx);
                    if rng.gen::<f64>() < daily_prob {
                        individual.vaccination_status[b_idx] = true;
                    }
                }
            }
        }
    }

    // --- drug updates---
    // Only count infections that have caused symptoms for treatment initiation decisions
    let has_any_infection = individual
        .level
        .iter()
        .enumerate()
        .any(|(b_idx, &level)| level > 0.0 && individual.infection_has_caused_symptoms[b_idx]);
    let initial_on_any_antibiotic = individual.cur_use_drug.iter().any(|&identified| identified);
    // Only count identified infections that also have symptoms (can't identify asymptomatic infections clinically)
    let has_any_identified_infection = individual
        .test_identified_infection
        .iter()
        .enumerate()
        .any(|(b_idx, &identified)| identified && individual.infection_has_caused_symptoms[b_idx]);

    // --- count number of drugs currently being used ---
    let num_drugs_currently_used = individual.cur_use_drug.iter().filter(|&&on| on).count();

    let mut syndrome_administration_multiplier: f64 = 1.0;
    for &syndrome_id in individual.infectious_syndrome.iter() {
        if syndrome_id > 0 {
            let multiplier = store.syndrome.initiation_multiplier(syndrome_id as usize);
            syndrome_administration_multiplier = syndrome_administration_multiplier.max(multiplier);
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
                    // Check if bacteria treatment was recognized in current year
                    let current_year = 1930.0 + (time_step as f64 / 365.0);
                    if let Some(recognition_year) = store.bacteria.treatment_recognition_year(b_idx)
                    {
                        if current_year < recognition_year {
                            // Skip this bacteria - treatment not yet recognized, don't continue drugs for it
                            continue;
                        }
                    }

                    // Use potency_when_no_r to determine if drug is relevant for this bacteria
                    let drug_potency = store.drug_bacteria.potency(b_idx, drug_idx);
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
                let random_cessation_if_no_infection = store
                    .globals
                    .random_drug_cessation_probability_if_no_active_infection;
                if rng.gen_bool(random_cessation_if_no_infection) {
                    stop_drug = true;
                }
            } else {
                // Calculate bacteria-specific and region-specific cessation probability
                let base_cessation_prob = primary_bacteria_idx
                    .map(|bacteria_idx| store.bacteria.drug_cessation_probability[bacteria_idx])
                    .unwrap_or(random_drug_cessation_prob);

                // Apply regional multiplier based on individual's current region
                let region_multiplier = if individual.region_cur_in == Region::Home {
                    store.region.cessation_multiplier(individual.region_living)
                } else {
                    store.region.cessation_multiplier(individual.region_cur_in)
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
                       individual.bacteria_level_at_drug_start[bacteria_idx].is_some()
                    {
                        // Record cessation for restart window tracking
                        individual.drug_stopped_with_infection_day[bacteria_idx] =
                            Some(time_step as i32);
                        individual.bacteria_level_at_drug_cessation[bacteria_idx] =
                            Some(individual.level[bacteria_idx]);
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
        let drug_initial_level = store.drug.initial_level(drug_idx);
        if individual.cur_use_drug[drug_idx] {
            individual.cur_level_drug[drug_idx] = drug_initial_level;
        } else {
            // Use exponential decay based on drug-specific half-life
            let half_life_days = store.drug.half_life_days(drug_idx);
            let decay_constant = (2.0_f64).ln() / half_life_days; // k = ln(2) / t_half
            let decay_factor = (-decay_constant).exp(); // e^(-k*t) where t=1 day
            let new_drug_level = individual.cur_level_drug[drug_idx] * decay_factor;
            // Set levels below 0.001 (0.1% of standard dose) to exactly zero to avoid floating point artifacts
            individual.cur_level_drug[drug_idx] = if new_drug_level < 0.001 {
                0.0
            } else {
                new_drug_level
            };
        }
    }

    // --- Apply Drug Level Interactions ---
    // Calculate final drug levels considering pairwise pharmacokinetic interactions
    apply_drug_level_interactions(individual);

    // --- drug initiation (two-stage process) ---
    // Stage 1: Decide whether to start any antibiotic
    let available_drugs: Vec<usize> = DRUG_SHORT_NAMES
        .iter()
        .enumerate()
        .filter(|(_, &name)| {
            let avail = get_drug_availability_time_aware(
                name,
                &individual.region_cur_in.to_string(),
                Some(&individual.region_living.to_string()),
                time_step,
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
    let scaling_factor = if available_drugs_count < min_available_drugs && available_drugs_count > 0
    {
        (min_available_drugs as f64) / (available_drugs_count as f64)
    } else {
        1.0
    };

    // Restriction: if already using three or more drugs, cannot start another (allow up to 3 drugs for severe infections)
    if num_drugs_currently_used + drugs_initiated_this_time_step < 3 && available_drugs_count > 0 {
        // Stage 1: Calculate probability to start any antibiotic
        let mut start_any_antibiotic_prob = drug_base_initiation_rate * scaling_factor;
        let infection_acquired_this_step = individual
            .date_last_infected
            .iter()
            .any(|&d| d == time_step as i32);
        if has_any_infection && !infection_acquired_this_step {
            start_any_antibiotic_prob *= drug_infection_present_multiplier;
        }
        if has_any_identified_infection {
            start_any_antibiotic_prob *= drug_test_identified_multiplier;
        }
        if initial_on_any_antibiotic || drugs_initiated_this_time_step > 0 {
            start_any_antibiotic_prob *= already_on_drug_initiation_multiplier;
        }
        // Immunocompromised patients more likely to receive prophylactic antibiotics
        if individual.immunodeficiency_type.is_some() {
            start_any_antibiotic_prob *=
                store.globals.immunodeficiency_prophylactic_drug_multiplier;
        }
        start_any_antibiotic_prob *= syndrome_administration_multiplier;
        start_any_antibiotic_prob = start_any_antibiotic_prob.clamp(0.0, 1.0);

        if rng.gen_bool(start_any_antibiotic_prob) {
            // Identify primary bacteria for drug score tracking (highest level among infected bacteria)
            let mut primary_bacteria_idx = -1i32;
            let mut highest_bacteria_level = 0.0;
            for b_idx in 0..BACTERIA_LIST.len() {
                if individual.level[b_idx] > 0.001
                    && individual.level[b_idx] > highest_bacteria_level
                {
                    highest_bacteria_level = individual.level[b_idx];
                    primary_bacteria_idx = b_idx as i32;
                }
            }

            // Store primary bacteria index for this drug selection event
            individual.bacteria_on_selection_day = primary_bacteria_idx;

            // Stage 2: Choose the most appropriate drug using weighted probabilistic selection
            // Score each available drug and collect scores for probabilistic selection
            // NOTE: TB-specific multi-drug selection logic is not implemented here because:
            // 1. Current potency-based scoring already favors effective TB drugs (rifampicin=0.6, FQs=0.4-0.5)
            // 2. Multi-drug synergy system activates automatically when ≥2 TB drugs are selected
            // 3. Implementing TB-specific simultaneous multi-drug initiation would require substantial
            //    modification to this single-drug selection framework
            // 4. Clinical TB programs often start with sequential drug addition anyway due to tolerance testing
            let mut drug_scores: Vec<(usize, f64)> = Vec::new();
            for &drug_idx in &available_drugs {
                let drug_name = DRUG_SHORT_NAMES[drug_idx];
                // Restriction: do not start drug if resistance test has been performed and resistance detected for any bacteria
                let mut resistance_detected = false;
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.test_for_resistance[b_idx]
                        && individual.resistances[b_idx][drug_idx].test_r > 0.0
                    {
                        resistance_detected = true;
                        break;
                    }
                }
                if resistance_detected {
                    continue;
                }

                // Score drug based on spectrum, activity, and clinical scenario
                let mut score = 1.0;

                // INTRINSIC ACTIVITY GATE: Block drugs with no meaningful activity against current infections
                let minimal_potency_threshold =
                    store.globals.minimal_potency_threshold_for_drug_selection;
                let mut has_meaningful_activity = false;
                let mut max_potency_against_infections: f64 = 0.0;

                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[b_idx] > 0.001 {
                        // Check if bacteria treatment was recognized in current year
                        let current_year = 1930.0 + (time_step as f64 / 365.0);
                        if let Some(recognition_year) =
                            store.bacteria.treatment_recognition_year(b_idx)
                        {
                            if current_year < recognition_year {
                                // Skip this bacteria - treatment not yet recognized
                                continue;
                            }
                        }

                        let potency = store.drug_bacteria.potency(b_idx, drug_idx);
                        max_potency_against_infections =
                            max_potency_against_infections.max(potency);
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
                            ("Pseudomonas aeruginosa", "imipenem_c") => score *= 20.0,
                            ("Pseudomonas aeruginosa", "ciprofloxacin") => score *= 18.0,
                            ("Pseudomonas aeruginosa", "tobramycin") => score *= 15.0,
                            ("Pseudomonas aeruginosa", "colistin") => score *= 12.0,
                            (
                                "Pseudomonas aeruginosa",
                                "penicilling" | "ampicillin" | "amoxicillin" | "cephalexin"
                                | "ceftriaxone" | "vancomycin",
                            ) => {
                                score = 0.0; // Completely block - no intrinsic activity
                                break;
                            }

                            // Staphylococcus aureus - DRAMATICALLY strengthen MSSA vs MRSA logic
                            ("Staphylococcus aureus", "penicilling") => {
                                // Early periods: penicillin should dominate (MSSA era)
                                if time_step < 7300 {
                                    // First ~20 years
                                    score *= 50.0; // MASSIVE boost for MSSA
                                } else {
                                    score *= 2.0; // Minimal in MRSA era
                                }
                            }
                            (
                                "Staphylococcus aureus",
                                "amoxicillin_clavulanate" | "ampicillin_sulbactam",
                            ) => {
                                if time_step < 10950 {
                                    // First ~30 years
                                    score *= 40.0; // Major boost before MRSA dominance
                                } else {
                                    score *= 3.0; // Reduced in MRSA era
                                }
                            }
                            ("Staphylococcus aureus", "vancomycin") => {
                                if time_step < 7300 {
                                    // Early years
                                    score *= 2.0; // Minimal early use
                                } else {
                                    // MRSA era
                                    score *= 35.0; // MASSIVE boost for MRSA
                                }
                            }
                            ("Staphylococcus aureus", "linezolid" | "tedizolid") => {
                                if time_step >= 10950 {
                                    // Late period only
                                    score *= 25.0; // Strong alternatives to vancomycin
                                } else {
                                    score *= 0.5; // Minimal early use
                                }
                            }
                            ("Staphylococcus aureus", "clindamycin") => score *= 8.0,

                            // E. coli - MASSIVELY strengthen first-line agents
                            ("Escherichia coli", "ciprofloxacin") => score *= 35.0, // Major UTI drug
                            ("Escherichia coli", "nitrofurantoin") => score *= 30.0, // Cystitis first-line
                            ("Escherichia coli", "trim_sulf") => score *= 25.0,
                            ("Escherichia coli", "ceftriaxone") => score *= 20.0, // Serious infections
                            ("Escherichia coli", "ampicillin") => {
                                if time_step < 7300 {
                                    // Early susceptible era
                                    score *= 25.0;
                                } else {
                                    score *= 3.0; // Resistance emerged
                                }
                            }
                            ("Escherichia coli", "meropenem" | "imipenem_c") => {
                                // Carbapenems should be rare for E. coli except ESBL era
                                if time_step >= 14600 {
                                    // Later periods for ESBL
                                    score *= 8.0;
                                } else {
                                    score *= 0.2; // Minimal early use
                                }
                            }

                            // Klebsiella pneumoniae - strengthen appropriate agents
                            ("Klebsiella pneumoniae", "ceftriaxone") => {
                                if time_step < 10950 {
                                    // Before ESBL dominance
                                    score *= 25.0;
                                } else {
                                    score *= 8.0;
                                }
                            }
                            ("Klebsiella pneumoniae", "meropenem" | "imipenem_c") => {
                                if time_step >= 10950 {
                                    // ESBL era
                                    score *= 30.0;
                                } else {
                                    score *= 3.0;
                                }
                            }
                            ("Klebsiella pneumoniae", "ciprofloxacin") => score *= 15.0,
                            ("Klebsiella pneumoniae", "piperacillin_tazobactam") => score *= 18.0,

                            // Enterococcus faecalis - strengthen appropriate agents
                            ("Enterococcus faecalis", "ampicillin") => score *= 40.0, // First-line
                            ("Enterococcus faecalis", "vancomycin") => {
                                if time_step >= 10950 {
                                    // VRE era
                                    score *= 30.0;
                                } else {
                                    score *= 8.0;
                                }
                            }
                            ("Enterococcus faecalis", "linezolid") => {
                                if time_step >= 14600 {
                                    // Late VRE era
                                    score *= 25.0;
                                } else {
                                    score *= 2.0;
                                }
                            }

                            // Enterococcus faecium - more resistant, different pattern
                            ("Enterococcus faecium", "ampicillin") => score *= 5.0, // Less effective than faecalis
                            ("Enterococcus faecium", "vancomycin") => {
                                if time_step >= 10950 {
                                    score *= 35.0;
                                } else {
                                    score *= 15.0;
                                }
                            }
                            ("Enterococcus faecium", "linezolid") => {
                                if time_step >= 14600 {
                                    score *= 30.0;
                                } else {
                                    score *= 3.0;
                                }
                            }
                            ("Enterococcus faecium", "quinu_dalfo") => {
                                if time_step >= 16425 {
                                    // Very late introduction
                                    score *= 20.0;
                                }
                            }

                            // Acinetobacter baumannii - highly resistant pathogen
                            ("Acinetobacter baumannii", "meropenem" | "imipenem_c") => {
                                if time_step < 18250 {
                                    // Before extensive carbapenem resistance
                                    score *= 40.0;
                                } else {
                                    score *= 15.0;
                                }
                            }
                            ("Acinetobacter baumannii", "colistin") => {
                                if time_step >= 14600 {
                                    // Later periods for MDR
                                    score *= 35.0;
                                } else {
                                    score *= 8.0;
                                }
                            }
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
                            "Pseudomonas aeruginosa" => vec![
                                "piperacillin_tazobactam",
                                "meropenem",
                                "imipenem_c",
                                "ceftazidime",
                                "cefepime",
                                "ciprofloxacin",
                                "tobramycin",
                            ],
                            "Staphylococcus aureus" => vec![
                                "penicilling",
                                "amoxicillin_clavulanate",
                                "ampicillin_sulbactam",
                                "vancomycin",
                                "linezolid",
                                "tedizolid",
                                "clindamycin",
                            ],
                            "Escherichia coli" => vec![
                                "ciprofloxacin",
                                "nitrofurantoin",
                                "trim_sulf",
                                "ceftriaxone",
                                "ampicillin",
                                "cefuroxime",
                            ],
                            "Klebsiella pneumoniae" => vec![
                                "ceftriaxone",
                                "meropenem",
                                "imipenem_c",
                                "ciprofloxacin",
                                "piperacillin_tazobactam",
                                "ertapenem",
                            ],
                            "Enterococcus faecalis" => {
                                vec!["ampicillin", "vancomycin", "linezolid", "tedizolid"]
                            }
                            "Enterococcus faecium" => {
                                vec!["vancomycin", "linezolid", "tedizolid", "quinu_dalfo"]
                            }
                            "Acinetobacter baumannii" => vec![
                                "meropenem",
                                "imipenem_c",
                                "colistin",
                                "ampicillin_sulbactam",
                                "minocycline",
                            ],
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
                        // Check if bacteria treatment was recognized in current year
                        let current_year = 1930.0 + (time_step as f64 / 365.0);
                        if let Some(recognition_year) =
                            store.bacteria.treatment_recognition_year(b_idx)
                        {
                            if current_year < recognition_year {
                                // Skip this bacteria - treatment not yet recognized
                                continue;
                            }
                        }

                        let specific_multiplier =
                            store.drug_bacteria.initiation_multiplier(b_idx, drug_idx);
                        max_bacteria_specific_multiplier =
                            max_bacteria_specific_multiplier.max(specific_multiplier);
                    }
                }
                score *= max_bacteria_specific_multiplier;

                // Apply regional resistance surveillance penalty for empirical therapy
                if !has_any_identified_infection {
                    let mut regional_resistance_penalty = 1.0_f64;
                    if !majority_r_cache.is_empty() {
                        let region_idx = individual.region_cur_in as usize;
                        let hospital_status = individual.hospital_status.is_hospitalized();

                        let very_high_threshold =
                            store.globals.regional_resistance_threshold_very_high;
                        let high_threshold = store.globals.regional_resistance_threshold_high;
                        let moderate_threshold =
                            store.globals.regional_resistance_threshold_moderate;

                        let very_high_penalty = store.globals.regional_resistance_penalty_very_high;
                        let high_penalty = store.globals.regional_resistance_penalty_high;
                        let moderate_penalty = store.globals.regional_resistance_penalty_moderate;

                        for b_idx in 0..BACTERIA_LIST.len() {
                            let resistance_values = majority_r_cache.bucket(
                                region_idx,
                                hospital_status,
                                b_idx,
                                drug_idx,
                            );
                            if !resistance_values.is_empty() {
                                let resistance_cases = resistance_values.len() as f64;
                                let mut total_cases_estimate = resistance_cases;
                                for d_idx in 0..DRUG_SHORT_NAMES.len() {
                                    let other_resistance_values = majority_r_cache.bucket(
                                        region_idx,
                                        hospital_status,
                                        b_idx,
                                        d_idx,
                                    );
                                    total_cases_estimate = total_cases_estimate
                                        .max(other_resistance_values.len() as f64);
                                }

                                if total_cases_estimate > 0.0 {
                                    let resistance_prevalence =
                                        resistance_cases / total_cases_estimate;

                                    let resistance_penalty =
                                        if resistance_prevalence >= very_high_threshold {
                                            very_high_penalty
                                        } else if resistance_prevalence >= high_threshold {
                                            high_penalty
                                        } else if resistance_prevalence >= moderate_threshold {
                                            moderate_penalty
                                        } else {
                                            1.0
                                        };

                                    regional_resistance_penalty =
                                        regional_resistance_penalty.min(resistance_penalty);
                                }
                            }
                        }
                    }
                    score *= regional_resistance_penalty;
                }

                let drug_spectrum = store.drug.spectrum_breadth(drug_idx);
                if has_any_identified_infection {
                    let targeted_narrow_bonus =
                        store.globals.targeted_therapy_narrow_spectrum_bonus;
                    let targeted_broad_penalty =
                        store.globals.targeted_therapy_broad_spectrum_penalty;
                    let ineffective_drug_penalty =
                        store.globals.targeted_therapy_ineffective_drug_penalty;
                    let effective_potency_threshold = store
                        .globals
                        .effective_potency_threshold_for_targeted_therapy;

                    let mut has_good_activity = false;
                    let mut best_potency: f64 = 0.0;
                    for b_idx in 0..BACTERIA_LIST.len() {
                        if individual.test_identified_infection[b_idx]
                            && individual.level[b_idx] > 0.001
                        {
                            let potency = store.drug_bacteria.potency(b_idx, drug_idx);
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
                    let empiric_broad_bonus = store.globals.empiric_therapy_broad_spectrum_bonus;
                    let empiric_ineffective_penalty =
                        store.globals.empiric_therapy_ineffective_penalty;
                    let effective_potency_threshold = store
                        .globals
                        .effective_potency_threshold_for_empirical_therapy;

                    let mut has_any_activity = false;
                    for b_idx in 0..BACTERIA_LIST.len() {
                        if individual.level[b_idx] > 0.001 {
                            let potency = store.drug_bacteria.potency(b_idx, drug_idx);
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
                    time_step,
                );

                // Check if drug has been introduced yet
                let mut drug_introduced = false;
                if let Some(intro_time) = crate::config::get_drug_introduction_time_step(drug_name)
                {
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
                let selection_temperature = store.globals.drug_selection_temperature;

                // Apply temperature scaling: lower temp = more deterministic (clinically realistic)
                // Temperature of 0.5 = strongly favor best drugs, 1.0 = moderate, 2.0+ = random
                let weights: Vec<f64> = drug_scores
                    .iter()
                    .map(|(_, score)| (score / selection_temperature).exp())
                    .collect();

                // Handle edge case where all weights are zero or infinite
                let total_weight: f64 = weights.iter().sum();
                if total_weight > 0.0 && total_weight.is_finite() {
                    let dist = WeightedIndex::new(&weights).unwrap();
                    let chosen_idx = dist.sample(rng);
                    let chosen_drug_idx = drug_scores[chosen_idx].0;

                    // Initiate the selected drug
                    let drug_name = DRUG_SHORT_NAMES[chosen_drug_idx];
                    individual.cur_use_drug[chosen_drug_idx] = true;
                    individual.date_drug_initiated[chosen_drug_idx] = time_step as i32;
                    individual.date_drug_initiated_keep[chosen_drug_idx] = time_step as i32; // Persistent record
                    individual.ever_taken_drug[chosen_drug_idx] = true;

                    // Update drug counter
                    update_drug_counter(individual);
                    if individual.id == 1000001 {
                        println!(
                            "mod.rs   started {} - two-stage rate of starting was {:.4} (score: {:.3})",
                            drug_name,
                            start_any_antibiotic_prob,
                            drug_scores[chosen_idx].1
                        );
                    }
                    let mut chosen_initial_level = store.drug.initial_level(chosen_drug_idx);
                    if has_any_identified_infection && rng.gen_bool(double_dose_probability) {
                        let double_dose_multiplier =
                            store.drug.double_dose_multiplier(chosen_drug_idx);
                        chosen_initial_level *= double_dose_multiplier;
                    }
                    individual.cur_level_drug[chosen_drug_idx] = chosen_initial_level;

                    // Update treatment failure tracking for all infected bacteria
                    for bacteria_idx in 0..BACTERIA_LIST.len() {
                        if individual.level[bacteria_idx] > 0.0 {
                            // Record bacteria level at drug start and reset tracking
                            individual.bacteria_level_at_drug_start[bacteria_idx] =
                                Some(individual.level[bacteria_idx]);
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
        if individual.cur_level_drug[drug_idx] > 0.0 {
            let drug_toxicity_per_unit = store.drug.toxicity_per_unit_level_per_day(drug_idx);
            daily_drug_toxicity_increase +=
                individual.cur_level_drug[drug_idx] * drug_toxicity_per_unit;
        }
    }

    // Apply toxicity changes: increase from drugs, natural clearance when no drugs
    if daily_drug_toxicity_increase > 0.0 {
        // Drugs present: accumulate toxicity with maximum cap
        let max_toxicity = store.globals.max_toxicity_level;
        individual.current_toxicity = (individual.current_toxicity + daily_drug_toxicity_increase)
            .max(0.0)
            .min(max_toxicity);
    } else {
        // No drugs present: natural toxicity clearance (liver/kidney function)
        let toxicity_clearance_rate = store.globals.toxicity_clearance_rate_per_day;
        individual.current_toxicity =
            (individual.current_toxicity - toxicity_clearance_rate).max(0.0);
    }

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
                    rng,
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
            rng,
        );
    }

    // --- death

    if individual.date_of_death.is_none() {
        let mut cause: Option<String> = None;

        // --- New Logistic Background Mortality Model ---
        let mut total_log_odds = store.globals.background_mortality_baseline_log_odds;

        // Time-varying mortality component (1930-2035): reflects historical mortality decline
        let years_since_1930 = time_step as f64 / 365.0;
        let start_multiplier = store.globals.mortality_baseline_1930_multiplier;
        let end_multiplier = store.globals.mortality_baseline_2035_multiplier;
        let half_life_years = store.globals.mortality_improvement_half_life_years;

        // Exponential decay from start_multiplier to end_multiplier
        let decay_rate = (2.0_f64).ln() / half_life_years; // ln(2) / half_life
        let time_multiplier = end_multiplier
            + (start_multiplier - end_multiplier) * (-decay_rate * years_since_1930).exp();
        let time_log_odds_adjustment = time_multiplier.ln();
        total_log_odds += time_log_odds_adjustment;

        let age_years = individual.age as f64 / 365.0;

        // Age effects
        let log_odds_per_year = store.globals.log_odds_mortality_per_year_of_age;
        total_log_odds += age_years * log_odds_per_year;

        // Non-linear age effect for very elderly
        if age_years > 80.0 {
            let log_odds_age_squared = store.globals.log_odds_mortality_per_year_of_age_squared;
            total_log_odds += (age_years - 80.0).powi(2) * log_odds_age_squared;
        }

        // Regional effects
        total_log_odds += store.region.mortality_log_odds(individual.region_living);

        // Sex effects
        total_log_odds += store.sex.mortality_log_odds(&individual.sex_at_birth);

        // Immunosuppression effect
        if individual.immunodeficiency_type.is_some() {
            total_log_odds += store.globals.log_odds_mortality_immunosuppressed;
        }

        // Hospital status effect
        if matches!(individual.hospital_status, HospitalStatus::InHospital) {
            total_log_odds += store.globals.log_odds_mortality_hospitalized;
        }

        // Convert total log odds to probability
        let background_risk = 1.0 / (1.0 + (-total_log_odds).exp());

        individual.background_all_cause_mortality_rate = background_risk.min(1.0);
        let mut prob_not_dying = 1.0 - background_risk;

        let mut infection_non_sepsis_prob_not_dying = 1.0;
        let mut has_infection_non_sepsis_risk = false;

        let non_sepsis_level_threshold = store.globals.infection_non_sepsis_minimum_bacteria_level;
        let non_sepsis_level_coefficient = store.globals.infection_non_sepsis_log_odds_per_level;

        for (b_idx, level) in individual.level.iter().enumerate() {
            if *level <= non_sepsis_level_threshold {
                continue;
            }

            // Skip infections already progressing through sepsis pathway
            if individual.sepsis[b_idx] {
                continue;
            }

            has_infection_non_sepsis_risk = true;

            let mut log_odds = store.globals.infection_non_sepsis_base_log_odds;
            log_odds += store
                .bacteria
                .infection_non_sepsis_mortality_log_odds(b_idx);

            let syndrome_id = individual.infectious_syndrome[b_idx].max(0) as usize;
            log_odds += store.syndrome.non_sepsis_mortality_log_odds(syndrome_id);

            log_odds += non_sepsis_level_coefficient * level;

            if matches!(individual.hospital_status, HospitalStatus::InHospital) {
                log_odds += store.globals.infection_non_sepsis_log_odds_in_hospital;
            }

            let age_years = individual.age as f64 / 365.0;
            let age_adjustment = if age_years < 1.0 {
                store.globals.infection_non_sepsis_log_odds_age_infant
            } else if age_years < 18.0 {
                store.globals.infection_non_sepsis_log_odds_age_child
            } else if age_years < 65.0 {
                store.globals.infection_non_sepsis_log_odds_age_adult
            } else {
                store.globals.infection_non_sepsis_log_odds_age_elderly
            };
            log_odds += age_adjustment;

            if individual.immunodeficiency_type.is_some() {
                log_odds += store.globals.infection_non_sepsis_log_odds_immunosuppressed;
            }

            let probability = 1.0 / (1.0 + (-log_odds).exp());
            let probability = probability.clamp(0.0, 1.0);
            infection_non_sepsis_prob_not_dying *= 1.0 - probability;
        }

        let infection_non_sepsis_risk = if has_infection_non_sepsis_risk {
            1.0 - infection_non_sepsis_prob_not_dying.clamp(0.0, 1.0)
        } else {
            0.0
        };

        individual.current_infection_related_death_risk = infection_non_sepsis_risk;

        if infection_non_sepsis_risk > 0.0 {
            prob_not_dying *= 1.0 - infection_non_sepsis_risk;
            if cause.is_none() {
                cause = Some("infection_non_sepsis_related".to_string());
            }
        } else {
            individual.current_infection_related_death_risk = 0.0;
        }

        let has_sepsis = individual.sepsis.iter().any(|&status| status);
        if has_sepsis {
            // Calculate age-adjusted sepsis mortality risk
            let mut sepsis_death_risk = store.globals.base_sepsis_death_risk_per_day;

            // Apply age-based multiplier
            let age_years = individual.age as f64 / 365.0;
            let age_multiplier = if age_years < 1.0 {
                store.globals.sepsis_age_mortality_multiplier_infant
            } else if age_years < 18.0 {
                store.globals.sepsis_age_mortality_multiplier_child
            } else if age_years < 65.0 {
                store.globals.sepsis_age_mortality_multiplier_adult
            } else {
                store.globals.sepsis_age_mortality_multiplier_elderly
            };
            sepsis_death_risk *= age_multiplier;

            // Apply region-based multiplier (healthcare quality)
            let region_sepsis_multiplier = store
                .region
                .sepsis_mortality_multiplier(individual.region_living);
            sepsis_death_risk *= region_sepsis_multiplier;

            // Apply immunosuppression multiplier
            if individual.immunodeficiency_type.is_some() {
                sepsis_death_risk *= store.globals.sepsis_immunosuppressed_multiplier;
            }

            // Cap the risk at 1.0 (100%)
            sepsis_death_risk = sepsis_death_risk.min(1.0);

            prob_not_dying *= 1.0 - sepsis_death_risk;
            if cause.is_none() {
                cause = Some("sepsis_related".to_string());
            }
        }
        let mut drug_adverse_event_risk_for_individual = 0.0;
        let toxicity_death_risk_per_day = store.globals.drug_toxicity_death_risk_per_day;

        for drug_idx in 0..DRUG_SHORT_NAMES.len() {
            // Removed unused variable 'drug_name'
            if individual.cur_level_drug[drug_idx] > 0.0 {
                // Use only the global config parameter for drug toxicity death risk
                drug_adverse_event_risk_for_individual =
                    (drug_adverse_event_risk_for_individual + toxicity_death_risk_per_day).min(1.0);
            }
        }
        individual.mortality_risk_current_toxicity = drug_adverse_event_risk_for_individual;
        if drug_adverse_event_risk_for_individual > 0.0 {
            prob_not_dying *= 1.0 - drug_adverse_event_risk_for_individual;
            if cause.is_none() {
                cause = Some("drug_toxicity_related".to_string());
            }
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
                    "infection_non_sepsis_related" => {
                        InfectionResolutionType::DeathFromInfectionNonSepsis
                    }
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
                            InfectionResolutionType::DeathFromInfectionNonSepsis => 3,
                            InfectionResolutionType::DeathFromBackground => 4,
                            InfectionResolutionType::DeathFromToxicity => 5,
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
                let sepsis_duration =
                    (time_step as i32 - individual.sepsis_onset_day[b_idx]).max(0);
                let minimum_duration = store.globals.sepsis_minimum_duration_days;

                // Only allow recovery after minimum duration
                if sepsis_duration >= minimum_duration {
                    // Logistic regression model for sepsis recovery
                    let base_log_odds = store.globals.sepsis_base_log_odds_of_recovery_per_day;

                    let mut total_log_odds = base_log_odds;

                    // (1) Bacteria level effect - higher bacteria level decreases recovery probability
                    let bacteria_level_coefficient = store.globals.sepsis_log_odds_bacteria_level;
                    total_log_odds += individual.level[b_idx] * bacteria_level_coefficient;

                    // (2) Hospital status effect - being in hospital increases recovery probability
                    if individual.hospital_status.is_hospitalized() {
                        let hospital_coefficient = store.globals.sepsis_log_odds_in_hospital;
                        total_log_odds += hospital_coefficient;
                    }

                    // (3) Age effects with categories
                    let age_years = individual.age as f64 / 365.0;
                    let age_coefficient = if age_years < 1.0 {
                        store.globals.sepsis_log_odds_age_infant
                    } else if age_years < 18.0 {
                        store.globals.sepsis_log_odds_age_child
                    } else if age_years < 65.0 {
                        store.globals.sepsis_log_odds_age_adult
                    } else {
                        store.globals.sepsis_log_odds_age_elderly
                    };
                    total_log_odds += age_coefficient;

                    // (4) Severe immunosuppression effect
                    if individual.immunodeficiency_type.is_some() {
                        let immunosuppressed_coefficient =
                            store.globals.sepsis_log_odds_immunosuppressed;
                        total_log_odds += immunosuppressed_coefficient;
                    }

                    // (5) Region-specific effect (healthcare quality and ICU availability)
                    total_log_odds += store.region.sepsis_log_odds(individual.region_living);

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
            let region = individual.region_cur_in;
            let age_idx = crate::config::AgeCategoryParameters::age_category_index(individual.age);

            let mut log_odds = store.bacteria.acquisition_log_odds_baseline[b_idx]
                + store.age_categories.bacteria_age_log_odds(b_idx, age_idx)
                + store.region_bacteria.acquisition_log_odds(region, b_idx)
                + store
                    .age_categories
                    .bacteria_region_age_log_odds(region, b_idx, age_idx);

            // Vaccination status (binary effect)
            if individual.vaccination_status[b_idx] {
                log_odds += store.bacteria.log_odds_vaccinated[b_idx];
            }

            // Microbiome presence effect
            if individual.presence_microbiome[b_idx] {
                log_odds += store.bacteria.log_odds_microbiome_present[b_idx];
            }

            // Hospital-acquired effect
            if individual.hospital_status.is_hospitalized() {
                log_odds += store.bacteria.log_odds_hospital_acquired[b_idx];
            }

            // Convert log-odds to probability
            let mut acquisition_probability = 1.0 / (1.0 + (-log_odds).exp());

            // Apply historical MDR TB incidence modifier
            if bacteria == "mdr mycobacterium tuberculosis" {
                let simulation_year = 1930.0 + (time_step as f64 / 365.0);
                let mdr_tb_multiplier = if simulation_year < 1944.0 {
                    store.globals.mdr_tb_pre_antibiotic_era_multiplier
                } else if simulation_year < 1966.0 {
                    store.globals.mdr_tb_early_antibiotic_era_multiplier
                } else {
                    store.globals.mdr_tb_modern_era_multiplier
                };
                acquisition_probability *= mdr_tb_multiplier;
            }

            // --- microbiome presence (Carriage) ---
            // Carriage (asymptomatic colonization) is modeled separately from infection because:
            // 1. It's vastly more common than infection (e.g., 20-30% carry S. aureus, only ~1% infected)
            // 2. Carriers are the primary reservoir for resistance transmission in the population
            // 3. Antibiotic use disrupts normal microbiome, creating niches for pathogen colonization
            // 4. When carriers develop infections, they're highly likely to have resistant infections (carrier amplification)
            if !individual.presence_microbiome[b_idx] {
                // Logistic model for carriage acquisition (consistent framework with infection acquisition)
                // Baseline includes same demographic and geographic risk factors as infection, but with different
                // baseline probability (typically higher for carriage than infection)
                let mut log_odds = store.bacteria.acquisition_log_odds_baseline[b_idx]
                    + store.age_categories.bacteria_age_log_odds(b_idx, age_idx)
                    + store.region_bacteria.acquisition_log_odds(region, b_idx)
                    + store
                        .age_categories
                        .bacteria_region_age_log_odds(region, b_idx, age_idx);

                // Vaccination status (binary effect)
                if individual.vaccination_status[b_idx] {
                    log_odds += store.bacteria.log_odds_vaccinated[b_idx];
                }

                // Hospital-acquired effect
                if individual.hospital_status.is_hospitalized() {
                    log_odds += store.bacteria.log_odds_hospital_acquired[b_idx];
                }

                // Add the extra log-odds for microbiome vs infection (bacteria-specific)
                // This parameter shifts the baseline rate between carriage and infection (typically positive for carriage)
                log_odds += store.bacteria.microbiome_vs_infection_log_odds(b_idx);

                // --- Antibiotic disruption effect on carriage acquisition ---
                // MECHANISM: Antibiotics kill commensal bacteria, disrupting colonization resistance and creating
                // ecological niches that pathogenic bacteria can exploit. This is why C. difficile infections
                // spike during/after broad-spectrum antibiotic use, and why antibiotic courses increase MRSA
                // and ESBL-producing bacteria colonization risk.
                // EMPIRICAL BASIS: Studies show 5-15x increased colonization risk during antibiotic therapy,
                // persisting for weeks to months after cessation (we model acute effect here).
                let mut antibiotic_disruption_log_odds = 0.0;
                let mut acquisition_on_drug = false;

                for (d_idx, &drug_level) in individual.cur_level_drug.iter().enumerate() {
                    if drug_level > 0.1 {
                        // Only count drugs with meaningful levels
                        acquisition_on_drug = true;
                        antibiotic_disruption_log_odds +=
                            store.drug.microbiome_disruption_log_odds(d_idx);
                    }
                }
                log_odds += antibiotic_disruption_log_odds;

                // Convert log-odds to probability
                let microbiome_acquisition_probability = 1.0 / (1.0 + (-log_odds).exp());

                if rng.gen_bool(microbiome_acquisition_probability.clamp(0.0, 1.0)) {
                    individual.presence_microbiome[b_idx] = true;
                    // Track acquisition date for duration-dependent clearance modeling
                    // RATIONALE: Recent colonization is more easily cleared by immune response or antibiotics,
                    // while established colonization (months to years) is much more persistent.
                    // This mirrors clinical observations that recent MRSA carriers respond better to
                    // decolonization protocols than chronic carriers.
                    individual.date_microbiome_acquired[b_idx] = time_step as i32;
                    individual.microbiome_acquired_today[b_idx] = true;
                    individual.microbiome_acquired_on_drug_today[b_idx] = acquisition_on_drug;

                    // --- assign microbiome_r on new microbiome acquisition (same logic as infection resistance assignment) ---
                    let max_resistance_level = store.globals.max_resistance_level;

                    let is_hospital_acquired = individual.hospital_status.is_hospitalized();

                    let region_idx = individual.region_cur_in as usize;
                    let hospital_status_bool = individual.hospital_status.is_hospitalized();

                    for drug_name_static in DRUG_SHORT_NAMES.iter() {
                        let d_idx = *drug_indices.get(drug_name_static).unwrap();
                        let resistance_data = &mut individual.resistances[b_idx][d_idx];

                        // --- region/hospital-specific sampling for microbiome (same logic as infections) ---
                        let sampling_hospital_status = if is_hospital_acquired {
                            true // Hospital-acquired microbiome samples from hospitalized population
                        } else {
                            hospital_status_bool // Community-acquired microbiome samples based on current status
                        };

                        let majority_r_values_from_population = majority_r_cache.bucket(
                            region_idx,
                            sampling_hospital_status,
                            b_idx,
                            d_idx,
                        );
                        if let Some(&acquired_resistance_level) =
                            majority_r_values_from_population.choose(rng)
                        {
                            let clamped_level =
                                acquired_resistance_level.min(max_resistance_level).max(0.0);
                            resistance_data.microbiome_r = clamped_level;
                        } else if let Some(fallback_level) = majority_r_cache.fallback_mean(
                            region_idx,
                            sampling_hospital_status,
                            b_idx,
                            d_idx,
                        ) {
                            resistance_data.microbiome_r =
                                fallback_level.min(max_resistance_level).max(0.0);
                        } else {
                            resistance_data.microbiome_r = 0.0;
                        }

                        if resistance_data.microbiome_r > 0.0 {
                            resistance_data.microbiome_r = (resistance_data.microbiome_r
                                * store
                                    .globals
                                    .microbiome_resistance_multiplier_on_acquisition)
                                .min(max_resistance_level)
                                .max(0.0);
                        }
                    }
                    // --- end microbiome_r assignment ---
                }
            } else {
                // --- Enhanced microbiome clearance with logistic model ---
                // RATIONALE FOR LOGISTIC FRAMEWORK: Clearance is influenced by multiple independent factors
                // (duration of carriage, antibiotic pressure, immune response) that combine multiplicatively
                // in probability space, which translates to additive effects in log-odds space.
                // This allows us to model complex interactions while maintaining interpretable parameters.

                // Baseline clearance probability (bacteria-specific or default)
                // Represents spontaneous clearance rate from immune surveillance and microbial competition
                let baseline_clearance_prob = store
                    .bacteria
                    .microbiome_clearance_probability_per_day(b_idx);

                // Convert baseline probability to log-odds for additive modeling
                let baseline_log_odds =
                    (baseline_clearance_prob / (1.0 - baseline_clearance_prob)).ln();
                let mut clearance_log_odds = baseline_log_odds;

                // --- Duration effect: longer carriage = harder to clear (established colonization) ---
                // MECHANISM: Newly acquired bacteria are more susceptible to immune clearance and competition.
                // Over time, successful colonizers establish stable niches, develop biofilms, and evade
                // immune responses, making them progressively harder to eliminate.
                // EMPIRICAL BASIS: MRSA decolonization success: ~70% for recent carriers vs ~30% for chronic carriers.
                // S. aureus carriage often persists for months to years once established.
                // IMPLEMENTATION: Negative coefficient (longer duration → lower clearance probability),
                // with maximum effect cap to prevent unrealistic persistence.
                if individual.date_microbiome_acquired[b_idx] > 0 {
                    let days_carried =
                        (time_step as i32 - individual.date_microbiome_acquired[b_idx]) as f64;
                    let duration_coefficient = store.globals.carriage_duration_log_odds_coefficient; // Negative value
                    let max_duration_effect = store.globals.carriage_duration_max_log_odds_effect; // Negative cap
                    let duration_effect =
                        (days_carried * duration_coefficient).max(max_duration_effect);
                    clearance_log_odds += duration_effect;
                }

                // --- Antibiotic effect: active drugs targeting this bacteria increase clearance ---
                // MECHANISM: Antibiotics with activity against the colonizing pathogen can suppress or eliminate it,
                // even at sub-therapeutic concentrations insufficient to treat infection. This is why antibiotic
                // prophylaxis can prevent colonization, and why treatment of infections often clears carriage.
                // EMPIRICAL BASIS: Decolonization protocols use antibiotics (e.g., mupirocin for MRSA nasal carriage).
                // Treatment courses often clear S. aureus carriage as a side effect.
                // IMPLEMENTATION: activity_r already accounts for drug level, potency, and resistance.
                // Each unit of activity proportionally increases clearance log-odds.
                for (d_idx, &drug_level) in individual.cur_level_drug.iter().enumerate() {
                    if drug_level > 0.1 {
                        let resistance_data = &individual.resistances[b_idx][d_idx];
                        let activity = resistance_data.activity_r;
                        if activity > 0.1 {
                            // Only count drugs with meaningful activity
                            let clearance_boost = activity
                                * store
                                    .globals
                                    .antibiotic_clearance_log_odds_per_unit_activity;
                            clearance_log_odds += clearance_boost;
                        }
                    }
                }

                // Convert log-odds back to probability
                let clearance_probability = 1.0 / (1.0 + (-clearance_log_odds).exp());

                if rng.gen_bool(clearance_probability.clamp(0.0, 1.0)) {
                    individual.presence_microbiome[b_idx] = false;
                    individual.date_microbiome_acquired[b_idx] = 0; // Reset acquisition date for potential re-acquisition
                    individual.microbiome_cleared_today[b_idx] = true;
                }

                // --- de novo resistance emergence in microbiome when on drug ---
                if individual.presence_microbiome[b_idx] {
                    let max_resistance_level = store.globals.max_resistance_level;
                    for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                        let resistance_data = &mut individual.resistances[b_idx][d_idx];
                        let drug_level = individual.cur_level_drug[d_idx];
                        // Only consider emergence if drug is present and microbiome_r is low
                        if drug_level > 0.0001 && resistance_data.microbiome_r < 0.0001 {
                            // Use a specific parameter for microbiome resistance emergence if present, else fallback to general
                            let emergence_rate_baseline = store
                                .globals
                                .microbiome_resistance_emergence_rate_per_day_baseline;
                            let microbiome_r_emergence_level =
                                store.globals.any_r_emergence_level_on_first_emergence;

                            // Optionally, you could scale by drug level or other factors
                            let total_emergence_prob = emergence_rate_baseline; // * (drug_level / 10.0).clamp(0.0, 1.0);

                            if rng.gen_bool(total_emergence_prob.clamp(0.0, 1.0)) {
                                resistance_data.microbiome_r =
                                    microbiome_r_emergence_level.min(max_resistance_level);
                            }
                        }
                    }

                    use crate::simulation::population::ResistanceMechanism;

                    // --- mechanism emergence in microbiome under drug pressure ---
                    for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                        let drug_level = individual.cur_level_drug[d_idx];
                        if drug_level <= 0.0 {
                            continue;
                        }

                        for (mechanism_idx, mechanism) in
                            ResistanceMechanism::all().iter().enumerate()
                        {
                            if individual.resistance_mechanisms[b_idx][mechanism_idx] {
                                continue;
                            }

                            // Check if mechanism is relevant for this drug/bacteria combination
                            let mechanism_applicable = match (mechanism, DRUG_SHORT_NAMES[d_idx]) {
                                (ResistanceMechanism::ESBL, drug) => matches!(
                                    drug,
                                    "penicilling"
                                        | "ampicillin"
                                        | "amoxicillin"
                                        | "piperacillin"
                                        | "ticarcillin"
                                        | "cephalexin"
                                        | "cefazolin"
                                        | "cefuroxime"
                                        | "ceftriaxone"
                                        | "ceftazidime"
                                        | "cefepime"
                                        | "ceftaroline"
                                        | "aztreonam"
                                        | "amoxicillin_clavulanate"
                                        | "piperacillin_tazobactam"
                                        | "ampicillin_sulbactam"
                                        | "ticarcillin_clavulanate"
                                ),
                                (ResistanceMechanism::Carbapenemase, drug) => matches!(
                                    drug,
                                    "meropenem"
                                        | "imipenem_c"
                                        | "ertapenem"
                                        | "meropenem_vaborbactam"
                                ),
                                (ResistanceMechanism::SixteenSMethyltransferase, drug) => {
                                    matches!(drug, "gentamicin" | "tobramycin" | "amikacin")
                                }
                                (ResistanceMechanism::Qnr, drug) => matches!(
                                    drug,
                                    "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin"
                                ),
                                (ResistanceMechanism::ErmMethylation, drug) => matches!(
                                    drug,
                                    "erythromycin" | "azithromycin" | "clarithromycin"
                                ),
                                (ResistanceMechanism::VanType, drug) => {
                                    matches!(drug, "vancomycin" | "teicoplanin")
                                }
                                (ResistanceMechanism::MecA, drug) => {
                                    bacteria == "staphylococcus aureus"
                                        && matches!(
                                            drug,
                                            "penicilling"
                                                | "ampicillin"
                                                | "amoxicillin"
                                                | "cephalexin"
                                                | "cefazolin"
                                                | "cefuroxime"
                                                | "ceftriaxone"
                                                | "ceftazidime"
                                                | "cefepime"
                                                | "meropenem"
                                                | "imipenem_c"
                                                | "ertapenem"
                                        )
                                }
                                (ResistanceMechanism::EffluxOverexpression, _) => true,
                                (ResistanceMechanism::ReducedPermeability, _) => !matches!(
                                    bacteria,
                                    "staphylococcus aureus"
                                        | "streptococcus pneumoniae"
                                        | "streptococcus pyogenes"
                                        | "streptococcus agalactiae"
                                        | "enterococcus faecalis"
                                        | "enterococcus faecium"
                                ),
                                (ResistanceMechanism::TargetSiteMutation, _) => true,
                                (ResistanceMechanism::AmpC, drug) => matches!(
                                    drug,
                                    "penicilling"
                                        | "ampicillin"
                                        | "amoxicillin"
                                        | "piperacillin"
                                        | "ticarcillin"
                                        | "cephalexin"
                                        | "cefazolin"
                                        | "cefuroxime"
                                        | "ceftriaxone"
                                        | "amoxicillin_clavulanate"
                                        | "piperacillin_tazobactam"
                                        | "ampicillin_sulbactam"
                                        | "ticarcillin_clavulanate"
                                ),
                            };

                            if !mechanism_applicable {
                                continue;
                            }

                            let mechanism_emergence_rate =
                                store.resistance_mechanism.emergence_rate(mechanism_idx);

                            if rng.gen_bool(mechanism_emergence_rate.clamp(0.0, 1.0)) {
                                individual.resistance_mechanisms[b_idx][mechanism_idx] = true;
                            }
                        }
                    }
                    // --- end mechanism emergence in microbiome ---
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
                        let current_microbiome_r =
                            individual.resistances[b_idx][d_idx].microbiome_r;
                        let possible_transfer_r_microbiome = (current_any_r > 0.0
                            && current_microbiome_r == 0.0)
                            || (current_microbiome_r > 0.0 && current_any_r == 0.0);
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
                let donor_has_resistance = individual.level[donor_idx] > 0.001
                    || individual.presence_microbiome[donor_idx];
                if donor_has_resistance {
                    for recipient_idx in 0..BACTERIA_LIST.len() {
                        if recipient_idx == donor_idx {
                            continue;
                        }
                        let hgt_prob = store.hgt.probability(donor_idx, recipient_idx);
                        if hgt_prob > 0.0 && rng.gen::<f64>() < hgt_prob {
                            // Transfer resistance for all drugs
                            for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                                let donor_r = individual.resistances[donor_idx][drug_idx].any_r;
                                if donor_r > 0.0 {
                                    // Transfer to infection
                                    if individual.level[recipient_idx] > 0.001 {
                                        let prev_any_r =
                                            individual.resistances[recipient_idx][drug_idx].any_r;
                                        let new_any_r = donor_r.max(prev_any_r);
                                        individual.resistances[recipient_idx][drug_idx].any_r =
                                            new_any_r;
                                        if prev_any_r == 0.0 && new_any_r > 0.0 {
                                            // Inline mechanism assignment
                                            use crate::simulation::population::ResistanceMechanism;
                                            let mechanism_prob = store
                                                .globals
                                                .mechanism_assignment_probability_on_any_r_gain;
                                            for (mech_idx, _mechanism) in
                                                ResistanceMechanism::all().iter().enumerate()
                                            {
                                                let enhancement = store
                                                    .resistance_mechanism
                                                    .enhancement_multiplier(mech_idx);
                                                if enhancement <= new_any_r {
                                                    if rng.gen_bool(mechanism_prob) {
                                                        individual.resistance_mechanisms
                                                            [recipient_idx][mech_idx] = true;
                                                    }
                                                }
                                            }
                                            individual.how_resistance_acquired[recipient_idx][drug_idx] = Some(crate::simulation::population::ResistanceAcquisitionType::Hgt);
                                        }
                                    }
                                    // Transfer to microbiome
                                    if individual.presence_microbiome[recipient_idx] {
                                        let prev_any_r =
                                            individual.resistances[recipient_idx][drug_idx].any_r;
                                        let new_any_r = donor_r.max(prev_any_r);
                                        individual.resistances[recipient_idx][drug_idx].any_r =
                                            new_any_r;
                                        if prev_any_r == 0.0 && new_any_r > 0.0 {
                                            // Inline mechanism assignment
                                            use crate::simulation::population::ResistanceMechanism;
                                            let mechanism_prob = store
                                                .globals
                                                .mechanism_assignment_probability_on_any_r_gain;
                                            for (mech_idx, _mechanism) in
                                                ResistanceMechanism::all().iter().enumerate()
                                            {
                                                let enhancement = store
                                                    .resistance_mechanism
                                                    .enhancement_multiplier(mech_idx);
                                                if enhancement <= new_any_r {
                                                    if rng.gen_bool(mechanism_prob) {
                                                        individual.resistance_mechanisms
                                                            [recipient_idx][mech_idx] = true;
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
                let prevention_efficacy = store.globals.antibiotic_infection_prevention_efficacy;

                // Check each drug the person is currently taking
                for (drug_idx, &is_taking_drug) in individual.cur_use_drug.iter().enumerate() {
                    if is_taking_drug {
                        // Calculate effective activity using the same method as activity_r calculation
                        let base_potency = store.drug_bacteria.potency(b_idx, drug_idx);
                        let drug_current_level = individual.cur_level_drug[drug_idx];
                        let max_resistance_level = store.globals.max_resistance_level;
                        let resistance_level = individual.resistances[b_idx][drug_idx].any_r;
                        let normalized_any_r = resistance_level / max_resistance_level;
                        let effective_activity =
                            base_potency * drug_current_level * (1.0 - normalized_any_r);

                        // If drug has effective activity, it can prevent infection
                        if effective_activity > 0.5 {
                            // Threshold for effective prevention
                            if rng.gen_bool(prevention_efficacy) {
                                infection_prevented = true;
                                individual.infection_prevented_by_drug[b_idx] = true; // Track prevention event
                                break; // One effective drug is enough
                            }
                        }
                    }
                }

                // Only proceed with infection if not prevented by existing antibiotics
                if !infection_prevented {
                    let initial_level = store.bacteria.initial_infection_level(b_idx);
                    individual.level[b_idx] = initial_level;
                    individual.date_last_infected[b_idx] = time_step as i32;
                    individual.date_last_infected_keep[b_idx] = time_step as i32; // Keep persistent record
                    individual.clearance_ready_day[b_idx] = -1;

                    // --- probabilistic syndrome assignment ---
                    let syndrome_id = assign_syndrome_for_bacteria(bacteria, rng);
                    individual.infectious_syndrome[b_idx] = syndrome_id as i32;

                    let env_acquisition_chance =
                        store.bacteria.environmental_acquisition_proportion(b_idx);
                    individual.cur_infection_from_environment[b_idx] =
                        rng.gen::<f64>() < env_acquisition_chance;

                    individual.infection_hospital_acquired[b_idx] =
                        individual.hospital_status.is_hospitalized();

                    // --- any_r and majority_r setting logic on new infection acquisition ---
                    let max_resistance_level = store.globals.max_resistance_level;

                    // --- TB-specific logic: guaranteed rifampicin resistance for MDR-TB ---
                    let is_tb = bacteria == "mdr mycobacterium tuberculosis";

                    // Time-dependent MDR TB incidence (historically accurate)
                    let simulation_year = 1930.0 + (time_step as f64 / 365.0);

                    let guaranteed_rifampicin_resistance = if is_tb && simulation_year >= 1966.0 {
                        // Only apply guaranteed rifampicin resistance after rifampicin is available
                        get_global_param(
                            "mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance",
                        )
                        .unwrap_or(0.90)
                    } else {
                        0.0
                    };

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
                            if let Some(intro_time) =
                                crate::config::get_drug_introduction_time_step(drug_name_static)
                            {
                                if time_step >= intro_time {
                                    any_selecting_drug_introduced = true;
                                }
                            }

                            // If not yet introduced by direct drug, check cross-resistance groups
                            if !any_selecting_drug_introduced {
                                if let Some(cross_resistance_drug_groups) =
                                    cross_resistance_groups.get(&b_idx)
                                {
                                    for group in cross_resistance_drug_groups {
                                        if group.contains(&d_idx) {
                                            // This drug is in a cross-resistance group, check if any other drug in the group has been introduced
                                            for &other_drug_idx in group {
                                                if other_drug_idx != d_idx {
                                                    if let Some(other_drug_name) =
                                                        DRUG_SHORT_NAMES.get(other_drug_idx)
                                                    {
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
                                let sampling_hospital_status = if is_hospital_acquired {
                                    true
                                } else {
                                    hospital_status_bool
                                };

                                let majority_r_values_from_population = majority_r_cache.bucket(
                                    region_idx,
                                    sampling_hospital_status,
                                    b_idx,
                                    d_idx,
                                );

                                let assigned_level = majority_r_values_from_population
                                    .choose(rng)
                                    .copied()
                                    .or_else(|| {
                                        majority_r_cache.fallback_mean(
                                            region_idx,
                                            sampling_hospital_status,
                                            b_idx,
                                            d_idx,
                                        )
                                    });

                                if let Some(level) = assigned_level {
                                    let clamped_level = level.min(max_resistance_level).max(0.0);
                                    resistance_data.any_r = clamped_level;
                                    resistance_data.majority_r = clamped_level;
                                } else {
                                    resistance_data.any_r = 0.0;
                                    resistance_data.majority_r = 0.0;
                                }

                                // Inline mechanism assignment
                                use crate::simulation::population::ResistanceMechanism;
                                let mechanism_prob =
                                    store.globals.mechanism_assignment_probability_on_any_r_gain;
                                for (mech_idx, _mechanism) in
                                    ResistanceMechanism::all().iter().enumerate()
                                {
                                    let enhancement =
                                        store.resistance_mechanism.enhancement_multiplier(mech_idx);
                                    if enhancement <= resistance_data.any_r {
                                        if rng.gen_bool(mechanism_prob) {
                                            individual.resistance_mechanisms[b_idx][mech_idx] =
                                                true;
                                        }
                                    }
                                }
                                individual.how_resistance_acquired[b_idx][d_idx] = Some(
                                    crate::simulation::population::ResistanceAcquisitionType::AtInfectionEnv,
                                );
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

                            let majority_r_values_from_population = majority_r_cache.bucket(
                                region_idx,
                                sampling_hospital_status,
                                b_idx,
                                d_idx,
                            );
                            let assigned_level = majority_r_values_from_population
                                .choose(rng)
                                .copied()
                                .or_else(|| {
                                    majority_r_cache.fallback_mean(
                                        region_idx,
                                        sampling_hospital_status,
                                        b_idx,
                                        d_idx,
                                    )
                                });
                            if let Some(level) = assigned_level {
                                let clamped_level = level.min(max_resistance_level).max(0.0);
                                resistance_data.any_r = clamped_level;
                                resistance_data.majority_r = clamped_level;
                                // Inline mechanism assignment
                                use crate::simulation::population::ResistanceMechanism;
                                let mechanism_prob =
                                    store.globals.mechanism_assignment_probability_on_any_r_gain;
                                for (mech_idx, _mechanism) in
                                    ResistanceMechanism::all().iter().enumerate()
                                {
                                    let enhancement =
                                        store.resistance_mechanism.enhancement_multiplier(mech_idx);

                                    if enhancement <= resistance_data.any_r {
                                        if rng.gen_bool(mechanism_prob) {
                                            individual.resistance_mechanisms[b_idx][mech_idx] =
                                                true;
                                        }
                                    }
                                }
                                individual.how_resistance_acquired[b_idx][d_idx] = Some(
                                crate::simulation::population::ResistanceAcquisitionType::AtInfectionCommunity,
                            );
                            } else {
                                resistance_data.any_r = 0.0;
                                resistance_data.majority_r = 0.0;
                            }
                        }
                    }

                    // --- TB-specific guaranteed rifampicin resistance ---
                    if is_tb && guaranteed_rifampicin_resistance > 0.0 {
                        if let Some(&rifampicin_idx) = drug_indices.get("rifampicin") {
                            let resistance_data =
                                &mut individual.resistances[b_idx][rifampicin_idx];
                            let current_resistance =
                                resistance_data.majority_r.max(resistance_data.any_r);
                            if current_resistance < guaranteed_rifampicin_resistance {
                                resistance_data.majority_r = guaranteed_rifampicin_resistance;
                                resistance_data.any_r = guaranteed_rifampicin_resistance;

                                // Add resistance mechanism for rifampicin resistance
                                use crate::simulation::population::ResistanceMechanism;
                                let mechanism_prob =
                                    store.globals.mechanism_assignment_probability_on_any_r_gain;
                                for (mech_idx, _mechanism) in
                                    ResistanceMechanism::all().iter().enumerate()
                                {
                                    let enhancement =
                                        store.resistance_mechanism.enhancement_multiplier(mech_idx);
                                    if enhancement <= resistance_data.any_r {
                                        if rng.gen_bool(mechanism_prob) {
                                            individual.resistance_mechanisms[b_idx][mech_idx] =
                                                true;
                                        }
                                    }
                                }
                                individual.how_resistance_acquired[b_idx][rifampicin_idx] = Some(crate::simulation::population::ResistanceAcquisitionType::AtInfectionTB);
                            }
                        }
                    }

                    // --- Carrier resistance inheritance (THE KEY MECHANISM FOR RESISTANCE AMPLIFICATION) ---
                    // WHY THIS MATTERS: Carriage is the primary mechanism by which resistance spreads in populations.
                    // Carriers are asymptomatic reservoirs who aren't on antibiotics, so resistant strains face no
                    // selective disadvantage in their microbiome. When carriers develop infections, the infecting
                    // strain is usually the one they carry, inheriting its resistance profile.
                    //
                    // EMPIRICAL BASIS:
                    // - MRSA carriers: 80-90% of their S. aureus infections are MRSA (vs ~30% in non-carriers)
                    // - ESBL-producing E. coli carriers: 70-80% of their UTIs are ESBL-positive
                    // - VRE carriers: >90% of subsequent bacteremias are VRE
                    //
                    // POPULATION-LEVEL IMPACT: This creates a "carrier amplification effect" where:
                    // 1. Antibiotics select for resistance in infections → some become carriers
                    // 2. Carriers maintain resistance without antibiotic pressure (no fitness cost)
                    // 3. When carriers get infected, resistance rates are much higher than population prevalence
                    // 4. This amplifies observed resistance rates beyond what direct selection would predict
                    //
                    // MECHANISM: When infection occurs, the bacteria causing infection typically comes from
                    // the person's own microbiome (endogenous infection) rather than the environment.
                    // We model this with high probability (default 85%) that carriers' infections inherit
                    // their microbiome resistance profile.
                    //
                    // IMPLEMENTATION NOTE: This inheritance occurs AFTER environmental/population-based resistance
                    // assignment, overriding it when carriage is present. This ensures carriers preferentially
                    // develop infections with their carried strain rather than acquiring new strains.
                    if individual.presence_microbiome[b_idx] {
                        let inheritance_prob =
                            store.globals.carrier_resistance_inheritance_probability;
                        if rng.gen_bool(inheritance_prob) {
                            let max_resistance_level = store.globals.max_resistance_level;
                            // Inherit microbiome resistance for all drugs
                            for d_idx in 0..DRUG_SHORT_NAMES.len() {
                                let microbiome_resistance =
                                    individual.resistances[b_idx][d_idx].microbiome_r;
                                if microbiome_resistance > 0.0 {
                                    let infection_resistance_data =
                                        &mut individual.resistances[b_idx][d_idx];
                                    // Inherit the higher of existing infection resistance or dampened microbiome resistance
                                    // (ensures we don't lose resistance already assigned from other sources)
                                    let dampened_microbiome_resistance = (microbiome_resistance
                                        * store.globals.infection_from_microbiome_dampening)
                                        .min(max_resistance_level)
                                        .max(0.0);
                                    let inherited_level = dampened_microbiome_resistance
                                        .max(infection_resistance_data.any_r);
                                    infection_resistance_data.any_r = inherited_level;
                                    infection_resistance_data.majority_r = inherited_level;

                                    // Track that this resistance came from microbiome carriage
                                    individual.how_resistance_acquired[b_idx][d_idx] = Some(
                                    crate::simulation::population::ResistanceAcquisitionType::FromMicrobiomeR,
                                );
                                }
                            }
                        }
                    }

                    // --- end generalized any_r and majority_r setting logic ---
                } // End if !infection_prevented block
            }
        } else {
            // Bacteria is already present (infection progression)
            // --- majority_r evolution ---
            let majority_r_evolution_rate =
                get_global_param("majority_r_evolution_rate_per_day_when_drug_present")
                    .unwrap_or(0.0);
            let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0); // Now using 1.0 from your config

            if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                for (drug_index, _use_drug) in individual.cur_use_drug.iter().enumerate() {
                    let resistance_data =
                        &mut individual.resistances[bacteria_full_idx][drug_index];

                    let drug_current_level = individual.cur_level_drug[drug_index];
                    let drug_currently_present = drug_current_level > 0.0; // Check if drug is effectively present
                    let current_bacteria_level = individual.level[b_idx];

                    // existing majority_r evolution based on drug presence
                    if resistance_data.majority_r == 0.0
                        && resistance_data.any_r > 0.0
                        && drug_currently_present
                    {
                        if rng.gen_bool(majority_r_evolution_rate) {
                            resistance_data.majority_r = resistance_data.any_r;
                        }
                    }

                    // any_r increase towards max_resistance_level
                    // when drug is present and majority_r is still 0
                    if resistance_data.majority_r == 0.0 && // No majority resistance yet
                       resistance_data.any_r > 0.0 && // But some minority resistance exists
                       resistance_data.any_r < max_resistance_level && // And it's not yet full resistance

                       drug_currently_present
                    // And the drug is present, providing selection pressure
                    {
                        let any_r_increase_rate =
                            store.globals.any_r_increase_rate_per_day_when_drug_present;
                        resistance_data.any_r =
                            (resistance_data.any_r + any_r_increase_rate).min(max_resistance_level);
                    }

                    // majority_r and any_r between 0 and 1
                    resistance_data.majority_r = resistance_data
                        .majority_r
                        .min(max_resistance_level)
                        .max(0.0);
                    resistance_data.any_r =
                        resistance_data.any_r.min(max_resistance_level).max(0.0);

                    //new resistance emergence ---
                    // this section handles the de novo emergence of resistance when it's not already present.
                    // it should come before activity_r is fully calculated for use in bacteria level reduction *this* time step.

                    if resistance_data.any_r < 0.0001 {
                        // Check if any_r is effectively zero
                        // only consider emergence if there's drug present (either being taken or decaying)
                        // and a positive bacteria level for selection pressure.
                        if drug_current_level > 0.0 && current_bacteria_level > 0.0001 {
                            let emergence_rate_baseline = store
                                .drug_bacteria
                                .resistance_emergence_rate(bacteria_full_idx, drug_index);
                            let bacteria_level_effect_multiplier =
                                store.globals.resistance_emergence_bacteria_level_multiplier;
                            let any_r_emergence_level_on_first_emergence =
                                store.globals.any_r_emergence_level_on_first_emergence;

                            // bacteria level dependency: Higher at higher levels
                            let max_bacteria_level = store.bacteria.max_level[bacteria_full_idx];
                            // Normalize bacteria level to [0,1] and apply multiplier
                            let bacteria_level_factor =
                                (current_bacteria_level / max_bacteria_level).clamp(0.0, 1.0)
                                    * bacteria_level_effect_multiplier;

                            // activity_r dependency: Bell-shaped curve
                            // Use the drug's initial level for normalization to get a comparable drug concentration scale (0-10)
                            let drug_initial_level_for_normalization =
                                store.drug.initial_level(drug_index);

                            // Normalize current drug level for bell-shaped emergence probability curve
                            let mut norm_drug_level =
                                drug_current_level / drug_initial_level_for_normalization;
                            norm_drug_level = norm_drug_level.clamp(0.0, 10.0);

                            // resistance emergence probability
                            // bell-shaped curve: 0.02 * x * (10 - x). Peaks at 5.0, is 0.1 at 0 and 10.
                            let emergence_drug_concentration_factor =
                                0.1 + 0.02 * norm_drug_level * (10.0 - norm_drug_level);
                            let emergence_drug_factor =
                                emergence_drug_concentration_factor.clamp(0.0, 1.0);

                            // Calculate multi-drug penalty if multiple drugs are active
                            let active_drug_count = individual
                                .cur_level_drug
                                .iter()
                                .filter(|&&level| level > 0.0)
                                .count();

                            let multi_drug_penalty_threshold =
                                store.globals.multi_drug_penalty_threshold_num_drugs as usize;
                            let mut multi_drug_penalty_factor = 1.0;

                            if active_drug_count >= multi_drug_penalty_threshold {
                                // Count how many active drugs this potential resistance would affect
                                let drugs_affected_by_this_resistance =
                                    if let Some(cross_resistance_groups) =
                                        crate::config::get_cross_resistance_groups().get(bacteria)
                                    {
                                        let current_drug_name = DRUG_SHORT_NAMES[drug_index];
                                        let mut affected_count = 0;

                                        // Check if this drug is in any cross-resistance group
                                        for group in cross_resistance_groups {
                                            if group.contains(&current_drug_name) {
                                                // Count how many drugs in this cross-resistance group are currently active
                                                for &group_drug in group {
                                                    if let Some(group_drug_idx) = DRUG_SHORT_NAMES
                                                        .iter()
                                                        .position(|&d| d == group_drug)
                                                    {
                                                        if individual.cur_level_drug[group_drug_idx]
                                                            > 0.0
                                                        {
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
                                        multi_drug_penalty_factor = store
                                            .globals
                                            .resistance_development_inhibition_single_drug;
                                    } else {
                                        // Partial cross-resistance among multiple active drugs
                                        multi_drug_penalty_factor = store
                                            .globals
                                            .resistance_development_inhibition_partial_cross;
                                    }
                                }
                                // If drugs_affected_by_this_resistance >= active_drug_count, no penalty (full cross-resistance)
                            }

                            // total emergence probability with multi-drug penalty
                            // adding 1.0 to bacteria_level_factor ensures a base contribution even if multiplier is low
                            let total_emergence_prob = emergence_rate_baseline
                                * (1.0 + bacteria_level_factor)
                                * emergence_drug_factor
                                * multi_drug_penalty_factor;

                            if rng.gen_bool(total_emergence_prob.clamp(0.0, 1.0)) {
                                resistance_data.any_r = any_r_emergence_level_on_first_emergence;
                                // Inline mechanism assignment
                                use crate::simulation::population::ResistanceMechanism;
                                let mechanism_prob =
                                    store.globals.mechanism_assignment_probability_on_any_r_gain;
                                for (mech_idx, _mechanism) in
                                    ResistanceMechanism::all().iter().enumerate()
                                {
                                    let enhancement =
                                        store.resistance_mechanism.enhancement_multiplier(mech_idx);
                                    if enhancement <= resistance_data.any_r {
                                        if rng.gen_bool(mechanism_prob) {
                                            individual.resistance_mechanisms[b_idx][mech_idx] =
                                                true;
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

                        if let Some(bacteria_full_idx) =
                            BACTERIA_LIST.iter().position(|&b| b == bacteria)
                        {
                            for (mechanism_idx, mechanism) in
                                ResistanceMechanism::all().iter().enumerate()
                            {
                                // Skip if mechanism already present
                                if individual.resistance_mechanisms[bacteria_full_idx]
                                    [mechanism_idx]
                                {
                                    continue;
                                }

                                // Check if this mechanism is relevant for current drug
                                let mechanism_applicable =
                                    match (mechanism, DRUG_SHORT_NAMES[drug_index]) {
                                        // ESBL affects beta-lactams (except carbapenems)
                                        (ResistanceMechanism::ESBL, drug) => {
                                            matches!(
                                                drug,
                                                "penicilling"
                                                    | "ampicillin"
                                                    | "amoxicillin"
                                                    | "piperacillin"
                                                    | "ticarcillin"
                                                    | "cephalexin"
                                                    | "cefazolin"
                                                    | "cefuroxime"
                                                    | "ceftriaxone"
                                                    | "ceftazidime"
                                                    | "cefepime"
                                                    | "ceftaroline"
                                                    | "aztreonam"
                                                    | "amoxicillin_clavulanate"
                                                    | "piperacillin_tazobactam"
                                                    | "ampicillin_sulbactam"
                                                    | "ticarcillin_clavulanate"
                                            )
                                        }
                                        // Carbapenemase affects carbapenems
                                        (ResistanceMechanism::Carbapenemase, drug) => {
                                            matches!(
                                                drug,
                                                "meropenem"
                                                    | "imipenem_c"
                                                    | "ertapenem"
                                                    | "meropenem_vaborbactam"
                                            )
                                        }
                                        // 16S methyltransferase affects aminoglycosides
                                        (ResistanceMechanism::SixteenSMethyltransferase, drug) => {
                                            matches!(drug, "gentamicin" | "tobramycin" | "amikacin")
                                        }
                                        // Qnr affects quinolones
                                        (ResistanceMechanism::Qnr, drug) => {
                                            matches!(
                                                drug,
                                                "ciprofloxacin"
                                                    | "levofloxacin"
                                                    | "moxifloxacin"
                                                    | "ofloxacin"
                                            )
                                        }
                                        // Erm methylation affects macrolides
                                        (ResistanceMechanism::ErmMethylation, drug) => {
                                            matches!(
                                                drug,
                                                "erythromycin" | "azithromycin" | "clarithromycin"
                                            )
                                        }
                                        // Van-type affects glycopeptides
                                        (ResistanceMechanism::VanType, drug) => {
                                            matches!(drug, "vancomycin" | "teicoplanin")
                                        }
                                        // mecA affects beta-lactams in Staph aureus
                                        (ResistanceMechanism::MecA, drug) => {
                                            bacteria == "staphylococcus aureus"
                                                && matches!(
                                                    drug,
                                                    "penicilling"
                                                        | "ampicillin"
                                                        | "amoxicillin"
                                                        | "cephalexin"
                                                        | "cefazolin"
                                                        | "cefuroxime"
                                                        | "ceftriaxone"
                                                        | "ceftazidime"
                                                        | "cefepime"
                                                        | "meropenem"
                                                        | "imipenem_c"
                                                        | "ertapenem"
                                                )
                                        }
                                        // Efflux overexpression can affect multiple drug classes
                                        (ResistanceMechanism::EffluxOverexpression, _) => true,
                                        // Reduced permeability affects many drugs, especially in Gram-negatives
                                        (ResistanceMechanism::ReducedPermeability, _) => !matches!(
                                            bacteria,
                                            "staphylococcus aureus"
                                                | "streptococcus pneumoniae"
                                                | "streptococcus pyogenes"
                                                | "streptococcus agalactiae"
                                                | "enterococcus faecalis"
                                                | "enterococcus faecium"
                                        ),
                                        // Target site mutations can affect various drugs
                                        (ResistanceMechanism::TargetSiteMutation, _) => true,
                                        // AmpC affects beta-lactams
                                        (ResistanceMechanism::AmpC, drug) => {
                                            matches!(
                                                drug,
                                                "penicilling"
                                                    | "ampicillin"
                                                    | "amoxicillin"
                                                    | "piperacillin"
                                                    | "ticarcillin"
                                                    | "cephalexin"
                                                    | "cefazolin"
                                                    | "cefuroxime"
                                                    | "ceftriaxone"
                                                    | "amoxicillin_clavulanate"
                                                    | "piperacillin_tazobactam"
                                                    | "ampicillin_sulbactam"
                                                    | "ticarcillin_clavulanate"
                                            )
                                        }
                                    };

                                if mechanism_applicable {
                                    let mechanism_emergence_rate =
                                        store.resistance_mechanism.emergence_rate(mechanism_idx);

                                    if rng.gen_bool(mechanism_emergence_rate.clamp(0.0, 1.0)) {
                                        individual.resistance_mechanisms[bacteria_full_idx]
                                            [mechanism_idx] = true;
                                    }
                                }
                            }
                        }
                    }
                    // --- end resistance mechanism emergence logic ---

                    // calculate activity_r (should always be updated) - but only when both drug and bacteria are present
                    // First check what the bacteria level will be after this timestep
                    let current_bacteria_level = individual.level[bacteria_full_idx];

                    if drug_current_level > 0.0 && current_bacteria_level > 0.001 {
                        // Fetch potency from indexed parameter store
                        let base_potency =
                            store.drug_bacteria.potency(bacteria_full_idx, drug_index);

                        // Calculate resistance mechanism enhancement
                        let mut mechanism_resistance_boost = 0.0;
                        if let Some(bacteria_full_idx) =
                            BACTERIA_LIST.iter().position(|&b| b == bacteria)
                        {
                            use crate::simulation::population::ResistanceMechanism;

                            for (mechanism_idx, mechanism) in
                                ResistanceMechanism::all().iter().enumerate()
                            {
                                if individual.resistance_mechanisms[bacteria_full_idx]
                                    [mechanism_idx]
                                {
                                    // Check if this mechanism affects the current drug
                                    let mechanism_affects_drug =
                                        match (mechanism, DRUG_SHORT_NAMES[drug_index]) {
                                            // ESBL affects beta-lactams (except carbapenems)
                                            (ResistanceMechanism::ESBL, drug) => {
                                                matches!(
                                                    drug,
                                                    "penicilling"
                                                        | "ampicillin"
                                                        | "amoxicillin"
                                                        | "piperacillin"
                                                        | "ticarcillin"
                                                        | "cephalexin"
                                                        | "cefazolin"
                                                        | "cefuroxime"
                                                        | "ceftriaxone"
                                                        | "ceftazidime"
                                                        | "cefepime"
                                                        | "ceftaroline"
                                                        | "aztreonam"
                                                        | "amoxicillin_clavulanate"
                                                        | "piperacillin_tazobactam"
                                                        | "ampicillin_sulbactam"
                                                        | "ticarcillin_clavulanate"
                                                )
                                            }
                                            // Carbapenemase affects carbapenems
                                            (ResistanceMechanism::Carbapenemase, drug) => {
                                                matches!(
                                                    drug,
                                                    "meropenem"
                                                        | "imipenem_c"
                                                        | "ertapenem"
                                                        | "meropenem_vaborbactam"
                                                )
                                            }
                                            // 16S methyltransferase affects aminoglycosides
                                            (
                                                ResistanceMechanism::SixteenSMethyltransferase,
                                                drug,
                                            ) => {
                                                matches!(
                                                    drug,
                                                    "gentamicin" | "tobramycin" | "amikacin"
                                                )
                                            }
                                            // Qnr affects quinolones
                                            (ResistanceMechanism::Qnr, drug) => {
                                                matches!(
                                                    drug,
                                                    "ciprofloxacin"
                                                        | "levofloxacin"
                                                        | "moxifloxacin"
                                                        | "ofloxacin"
                                                )
                                            }
                                            // Erm methylation affects macrolides
                                            (ResistanceMechanism::ErmMethylation, drug) => {
                                                matches!(
                                                    drug,
                                                    "erythromycin"
                                                        | "azithromycin"
                                                        | "clarithromycin"
                                                )
                                            }
                                            // Van-type affects glycopeptides
                                            (ResistanceMechanism::VanType, drug) => {
                                                matches!(drug, "vancomycin" | "teicoplanin")
                                            }
                                            // mecA affects beta-lactams in Staph aureus
                                            (ResistanceMechanism::MecA, drug) => {
                                                bacteria == "staphylococcus aureus"
                                                    && matches!(
                                                        drug,
                                                        "penicilling"
                                                            | "ampicillin"
                                                            | "amoxicillin"
                                                            | "cephalexin"
                                                            | "cefazolin"
                                                            | "cefuroxime"
                                                            | "ceftriaxone"
                                                            | "ceftazidime"
                                                            | "cefepime"
                                                            | "meropenem"
                                                            | "imipenem_c"
                                                            | "ertapenem"
                                                    )
                                            }
                                            // Efflux overexpression can affect multiple drug classes
                                            (ResistanceMechanism::EffluxOverexpression, _) => true,
                                            // Reduced permeability affects many drugs, especially in Gram-negatives
                                            (ResistanceMechanism::ReducedPermeability, _) => {
                                                !matches!(
                                                    bacteria,
                                                    "staphylococcus aureus"
                                                        | "streptococcus pneumoniae"
                                                        | "streptococcus pyogenes"
                                                        | "streptococcus agalactiae"
                                                        | "enterococcus faecalis"
                                                        | "enterococcus faecium"
                                                )
                                            }
                                            // Target site mutations can affect various drugs
                                            (ResistanceMechanism::TargetSiteMutation, _) => true,
                                            // AmpC affects beta-lactams
                                            (ResistanceMechanism::AmpC, drug) => {
                                                matches!(
                                                    drug,
                                                    "penicilling"
                                                        | "ampicillin"
                                                        | "amoxicillin"
                                                        | "piperacillin"
                                                        | "ticarcillin"
                                                        | "cephalexin"
                                                        | "cefazolin"
                                                        | "cefuroxime"
                                                        | "ceftriaxone"
                                                        | "amoxicillin_clavulanate"
                                                        | "piperacillin_tazobactam"
                                                        | "ampicillin_sulbactam"
                                                        | "ticarcillin_clavulanate"
                                                )
                                            }
                                        };

                                    if mechanism_affects_drug {
                                        let mechanism_enhancement = store
                                            .resistance_mechanism
                                            .enhancement_multiplier(mechanism_idx);

                                        // Only add enhancement if it would actually increase resistance
                                        // Mechanisms can't decrease resistance, but they also don't add if any_r is already higher
                                        let normalized_any_r =
                                            resistance_data.any_r / max_resistance_level;
                                        if mechanism_enhancement > normalized_any_r {
                                            let additional_resistance =
                                                mechanism_enhancement - normalized_any_r;
                                            mechanism_resistance_boost += additional_resistance;
                                        }
                                    }
                                }
                            }
                        }

                        // Apply mechanism enhancements to resistance levels if they would increase resistance
                        if mechanism_resistance_boost > 0.0 {
                            let normalized_any_r = resistance_data.any_r / max_resistance_level;
                            let new_resistance_level =
                                (normalized_any_r + mechanism_resistance_boost).min(1.0);
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
                        resistance_data.activity_r =
                            base_potency * drug_current_level * (1.0 - normalized_any_r);
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
        let bacterial_testing_available_from_day =
            get_global_param("bacterial_testing_available_from_day").unwrap_or(5478.0) as i32;
        let bacterial_testing_available =
            time_step >= bacterial_testing_available_from_day as usize;

        // Check bacteria-specific test availability for late-discovered bacteria (e.g., H. pylori 1982)
        // Most bacteria are available once general bacterial testing is available (~1945)
        // Only specific bacteria have delayed discovery dates
        let bacteria_name = BACTERIA_LIST[b_idx];
        let bacteria_param_name = bacteria_name.to_lowercase().replace(" ", "_");
        let bacteria_test_availability_param =
            format!("{}_test_availability_year", bacteria_param_name);
        let bacteria_specific_available = if let Some(bacteria_discovery_year) =
            get_global_param(&bacteria_test_availability_param)
        {
            let bacteria_discovery_day = ((bacteria_discovery_year - 1930.0) * 365.25) as i32;
            time_step >= bacteria_discovery_day as usize
        } else {
            bacterial_testing_available // For most bacteria, use the general bacterial testing availability
        };

        if is_infected
            && !individual.test_identified_infection[b_idx]
            && last_infected_time > 0
            && (time_step as i32) >= (last_infected_time + test_delay_days)
            && bacterial_testing_available
            && bacteria_specific_available
            && individual.infection_has_caused_symptoms[b_idx]
        {
            // Calculate comprehensive testing probability
            let testing_probability = calculate_testing_probability(
                individual,
                time_step,
                bacterial_testing_available_from_day as usize,
                param_cache,
                true, // is_bacterial_testing
            );

            if rng.gen_bool(testing_probability.clamp(0.0, 1.0)) {
                individual.test_identified_infection[b_idx] = true;
            }
        }

        // --- test_r assignment logic ---
        let test_r_error_prob = get_global_param("test_r_error_probability").unwrap_or(0.02);
        let test_r_error_value = get_global_param("test_r_error_value").unwrap_or(0.25);
        let resistance_test_result_delay_days =
            get_global_param("resistance_test_result_delay_days").unwrap_or(2.0) as i32;

        // Check if resistance testing is available yet (historically realistic dates)
        let resistance_testing_available_from_day =
            get_global_param("resistance_testing_available_from_day").unwrap_or(9131.0) as i32;
        let resistance_testing_available =
            time_step >= resistance_testing_available_from_day as usize;

        if individual.test_identified_infection[b_idx] && resistance_testing_available {
            // Check if we should initiate resistance testing (if not already initiated)
            if individual.resistance_test_initiated_day[b_idx] == -1 {
                // Calculate comprehensive resistance testing probability
                let resistance_testing_probability = calculate_testing_probability(
                    individual,
                    time_step,
                    resistance_testing_available_from_day as usize,
                    param_cache,
                    false, // is_bacterial_testing
                );

                if rng.gen_bool(resistance_testing_probability.clamp(0.0, 1.0)) {
                    // Set the flag indicating resistance testing was initiated
                    individual.test_for_resistance[b_idx] = true;
                    individual.resistance_test_initiated_day[b_idx] = time_step as i32;
                }
            }

            // Check if resistance test results should be available yet
            let test_initiated_day = individual.resistance_test_initiated_day[b_idx];
            if test_initiated_day != -1
                && (time_step as i32) >= (test_initiated_day + resistance_test_result_delay_days)
            {
                let test_r_already_set =
                    individual.resistances[b_idx].iter().any(|r| r.test_r > 0.0);
                if !test_r_already_set {
                    for d_idx in 0..DRUG_SHORT_NAMES.len() {
                        let any_r = individual.resistances[b_idx][d_idx].any_r;
                        let error = rng.gen_bool(test_r_error_prob);
                        let test_r = if error {
                            if any_r < 0.001 {
                                test_r_error_value
                            } else {
                                0.0
                            }
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
            let baseline_change = store.bacteria.base_level_change(b_idx);
            let mut total_reduction_due_to_antibiotic = 0.0;
            let mut immune_hazard = 0.0;
            let mut immune_clearance_triggered = false;

            let clearance_ready_day = individual.clearance_ready_day[b_idx];
            if clearance_ready_day != -1 && (time_step as i32) >= clearance_ready_day {
                immune_hazard = store
                    .clearance
                    .hazard_for(
                        b_idx,
                        individual.age,
                        individual.immunodeficiency_type.is_some(),
                        individual.level[b_idx],
                    )
                    .clamp(0.0, 1.0);

                if immune_hazard > 0.0 && rng.gen_bool(immune_hazard) {
                    immune_clearance_triggered = true;
                }
            }

            individual.clearance_hazard[b_idx] = immune_hazard;

            // --- Mechanism-specific fitness cost reversion logic ---
            let on_any_drug = individual.cur_level_drug.iter().any(|&lvl| lvl > 0.0);
            if !on_any_drug {
                // Check for reversion of specific resistance mechanisms based on their fitness costs
                use crate::simulation::population::ResistanceMechanism;
                let mut mechanisms_reverted = Vec::new();

                for (mechanism_idx, _) in ResistanceMechanism::all().iter().enumerate() {
                    if individual.resistance_mechanisms[b_idx][mechanism_idx] {
                        let mechanism_reversion_rate =
                            store.resistance_mechanism.reversion_rate(mechanism_idx);

                        if rng.gen_bool(mechanism_reversion_rate.clamp(0.0, 1.0)) {
                            individual.resistance_mechanisms[b_idx][mechanism_idx] = false;
                            mechanisms_reverted.push(mechanism_idx);
                        }
                    }
                }

                // If any mechanisms were lost, recalculate resistance levels for all drugs
                if !mechanisms_reverted.is_empty() {
                    for drug_index in 0..DRUG_SHORT_NAMES.len() {
                        let resistance_data = &mut individual.resistances[b_idx][drug_index];

                        // Recalculate mechanism-based resistance enhancement
                        let mut mechanism_resistance_boost = 0.0;
                        let max_resistance_level = store.globals.max_resistance_level;

                        for (mechanism_idx, mechanism) in
                            ResistanceMechanism::all().iter().enumerate()
                        {
                            if individual.resistance_mechanisms[b_idx][mechanism_idx] {
                                // Check if this mechanism affects the current drug (same logic as in calculation)
                                let mechanism_affects_drug =
                                    match (mechanism, DRUG_SHORT_NAMES[drug_index]) {
                                        (ResistanceMechanism::ESBL, drug) => {
                                            matches!(
                                                drug,
                                                "penicilling"
                                                    | "ampicillin"
                                                    | "amoxicillin"
                                                    | "piperacillin"
                                                    | "ticarcillin"
                                                    | "cephalexin"
                                                    | "cefazolin"
                                                    | "cefuroxime"
                                                    | "ceftriaxone"
                                                    | "ceftazidime"
                                                    | "cefepime"
                                                    | "ceftaroline"
                                                    | "aztreonam"
                                                    | "amoxicillin_clavulanate"
                                                    | "piperacillin_tazobactam"
                                                    | "ampicillin_sulbactam"
                                                    | "ticarcillin_clavulanate"
                                            )
                                        }
                                        (ResistanceMechanism::Carbapenemase, drug) => {
                                            matches!(
                                                drug,
                                                "meropenem"
                                                    | "imipenem_c"
                                                    | "ertapenem"
                                                    | "meropenem_vaborbactam"
                                            )
                                        }
                                        (ResistanceMechanism::SixteenSMethyltransferase, drug) => {
                                            matches!(drug, "gentamicin" | "tobramycin" | "amikacin")
                                        }
                                        (ResistanceMechanism::Qnr, drug) => {
                                            matches!(
                                                drug,
                                                "ciprofloxacin"
                                                    | "levofloxacin"
                                                    | "moxifloxacin"
                                                    | "ofloxacin"
                                            )
                                        }
                                        (ResistanceMechanism::ErmMethylation, drug) => {
                                            matches!(
                                                drug,
                                                "erythromycin" | "azithromycin" | "clarithromycin"
                                            )
                                        }
                                        (ResistanceMechanism::VanType, drug) => {
                                            matches!(drug, "vancomycin" | "teicoplanin")
                                        }
                                        (ResistanceMechanism::MecA, drug) => {
                                            bacteria == "staphylococcus aureus"
                                                && matches!(
                                                    drug,
                                                    "penicilling"
                                                        | "ampicillin"
                                                        | "amoxicillin"
                                                        | "cephalexin"
                                                        | "cefazolin"
                                                        | "cefuroxime"
                                                        | "ceftriaxone"
                                                        | "ceftazidime"
                                                        | "cefepime"
                                                        | "meropenem"
                                                        | "imipenem_c"
                                                        | "ertapenem"
                                                )
                                        }
                                        (ResistanceMechanism::EffluxOverexpression, _) => true,
                                        (ResistanceMechanism::ReducedPermeability, _) => !matches!(
                                            bacteria,
                                            "staphylococcus aureus"
                                                | "streptococcus pneumoniae"
                                                | "streptococcus pyogenes"
                                                | "streptococcus agalactiae"
                                                | "enterococcus faecalis"
                                                | "enterococcus faecium"
                                        ),
                                        (ResistanceMechanism::TargetSiteMutation, _) => true,
                                        (ResistanceMechanism::AmpC, drug) => {
                                            matches!(
                                                drug,
                                                "penicilling"
                                                    | "ampicillin"
                                                    | "amoxicillin"
                                                    | "piperacillin"
                                                    | "ticarcillin"
                                                    | "cephalexin"
                                                    | "cefazolin"
                                                    | "cefuroxime"
                                                    | "ceftriaxone"
                                                    | "amoxicillin_clavulanate"
                                                    | "piperacillin_tazobactam"
                                                    | "ampicillin_sulbactam"
                                                    | "ticarcillin_clavulanate"
                                            )
                                        }
                                    };

                                if mechanism_affects_drug {
                                    let mechanism_enhancement = store
                                        .resistance_mechanism
                                        .enhancement_multiplier(mechanism_idx);
                                    mechanism_resistance_boost += mechanism_enhancement;
                                }
                            }
                        }

                        // Update resistance levels based on remaining mechanisms
                        let base_resistance = resistance_data.any_r / max_resistance_level
                            - (resistance_data.any_r / max_resistance_level)
                                .min(mechanism_resistance_boost);
                        let new_resistance_level = (base_resistance + mechanism_resistance_boost)
                            .min(1.0)
                            .max(0.0);
                        let new_any_r = new_resistance_level * max_resistance_level;

                        resistance_data.any_r = new_any_r;
                        if resistance_data.majority_r > 0.0 {
                            resistance_data.majority_r = resistance_data.any_r;
                        }
                    }
                }
            }

            if individual.id == 1000001 {
                println!(" ");
                println!("mod.rs");
                println!("bacteria: {}", bacteria);
                println!("immune clearance hazard: {:.4}", immune_hazard);
                println!("baseline change: {:.4}", baseline_change);
            }

            for (drug_idx, _drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                if individual.cur_level_drug[drug_idx] > 0.0 {
                    let resistance_data = &individual.resistances[b_idx][drug_idx];
                    total_reduction_due_to_antibiotic += resistance_data.activity_r;

                    if individual.id == 1000001 {
                        // Calculate standardized MIC: 1 / ((1 - majority_r) * potency)
                        let potency = store.drug_bacteria.potency(b_idx, drug_idx);
                        let max_resistance_level = store.globals.max_resistance_level;
                        let normalized_majority_r =
                            resistance_data.majority_r / max_resistance_level;
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

            // --- TB-specific multi-drug synergy logic ---
            // WHY TB IS SPECIAL: Unlike other bacteria, TB has an absolute biological requirement for multi-drug therapy.
            // Single-drug TB treatment always fails due to rapid resistance development (~10^-6 mutation rate).
            // Other bacteria can often be treated with monotherapy, but TB biology (intracellular location,
            // thick cell wall, slow metabolism) requires sustained multi-drug pressure through different mechanisms.
            // This synergy bonus captures the mechanistic requirement that TB treatment guidelines mandate
            // ≥4 drugs initially, ≥2 for continuation - reflecting clinical reality, not just preference.
            let mut tb_synergy_bonus = 0.0;
            if bacteria == "mdr mycobacterium tuberculosis" {
                // Count active TB drugs with meaningful potency
                let active_tb_drugs: Vec<_> = DRUG_SHORT_NAMES
                    .iter()
                    .enumerate()
                    .filter(|(drug_idx, _drug_name)| {
                        if individual.cur_level_drug[*drug_idx] <= 0.0 {
                            return false;
                        }

                        let potency = store.drug_bacteria.potency(b_idx, *drug_idx);
                        potency >= 0.1 // Only count drugs with meaningful TB potency
                    })
                    .collect();

                let synergy_threshold =
                    get_global_param("mdr_mycobacterium_tuberculosis_multi_drug_synergy_threshold")
                        .unwrap_or(2.0) as usize;

                if active_tb_drugs.len() >= synergy_threshold {
                    let synergy_multiplier = get_global_param(
                        "mdr_mycobacterium_tuberculosis_multi_drug_synergy_multiplier",
                    )
                    .unwrap_or(2.5);
                    // Background effectiveness represents unmodeled TB-specific drugs (bedaquiline, pretomanid, delamanid,
                    // cycloserine, ethionamide, p-aminosalicylic acid) that are critical for MDR-TB treatment but not
                    // explicitly tracked in this general AMR model. Value reflects their collective contribution when
                    // proper multi-drug TB regimens are used.
                    let mut background_effectiveness = get_global_param(
                        "mdr_mycobacterium_tuberculosis_background_drug_effectiveness",
                    )
                    .unwrap_or(0.8);

                    // Apply historical treatment effectiveness modifier
                    let simulation_year = 1930.0 + (time_step as f64 / 365.0);
                    if simulation_year < 1944.0 {
                        // Pre-antibiotic era: no effective TB treatment available
                        background_effectiveness *= 0.01; // 99% reduction in effectiveness
                    } else if simulation_year < 1966.0 {
                        // Early antibiotic era: limited effectiveness with monotherapy
                        background_effectiveness *= 0.3; // 70% reduction in effectiveness
                    }
                    // Modern era (1966+): full effectiveness (no change needed)

                    // Apply synergy: multiply existing drug effects + add background effectiveness
                    tb_synergy_bonus = (total_reduction_due_to_antibiotic
                        * (synergy_multiplier - 1.0))
                        + background_effectiveness;

                    if individual.id == 1000001 {
                        println!(
                            "mod.rs  TB synergy: {} active drugs, bonus = {:.4}",
                            active_tb_drugs.len(),
                            tb_synergy_bonus
                        );
                    }
                }
            }

            // Antibiotic effectiveness is now determined through bacteria-drug specific potency values
            // rather than a universal treatment response modifier
            // TB synergy bonus is added here because multi-drug synergy is fundamental to TB treatment effectiveness -
            // it's not an optional enhancement but a biological requirement for meaningful bacterial killing
            let adjusted_antibiotic_effect = total_reduction_due_to_antibiotic + tb_synergy_bonus;

            if individual.id == 1000001 {
                println!(
                    "mod.rs  total reduction due to antibiotic: {:.4}",
                    total_reduction_due_to_antibiotic
                );
                println!(
                    "mod.rs  adjusted antibiotic effect: {:.4}",
                    adjusted_antibiotic_effect
                );
            }

            let decay = baseline_change - adjusted_antibiotic_effect;

            let max_level = store.bacteria.max_level(b_idx);
            let new_bacteria_level = (individual.level[b_idx] + decay).max(0.0).min(max_level);

            // Check for infection clearance before updating the level
            let old_level = individual.level[b_idx];

            if new_bacteria_level < 0.0001 || immune_clearance_triggered {
                // Check if there was an infection before clearance (previous level > 0.001)
                let was_previously_infected = old_level > 0.001;

                if was_previously_infected {
                    // Determine resolution type based on actual drug activity accounting for resistance
                    let has_active_drugs =
                        individual
                            .cur_use_drug
                            .iter()
                            .enumerate()
                            .any(|(drug_idx, &on_drug)| {
                                if !on_drug {
                                    return false;
                                }
                                // Check if this drug has current activity_r > 0 (accounts for resistance)
                                let activity_r = individual.resistances[b_idx][drug_idx].activity_r;
                                activity_r > 0.0
                            });

                    let resolution_type = if immune_clearance_triggered {
                        InfectionResolutionType::ImmuneClearance
                    } else if has_active_drugs {
                        InfectionResolutionType::DrugAssistedClearance
                    } else {
                        InfectionResolutionType::ImmuneClearance
                    };

                    let resolution_idx = match resolution_type {
                        InfectionResolutionType::ImmuneClearance => 0,
                        InfectionResolutionType::DrugAssistedClearance => 1,
                        InfectionResolutionType::DeathFromSepsis => 2,
                        InfectionResolutionType::DeathFromInfectionNonSepsis => 3,
                        InfectionResolutionType::DeathFromBackground => 4,
                        InfectionResolutionType::DeathFromToxicity => 5,
                    };
                    individual.infection_resolution_this_timestep[b_idx][resolution_idx] += 1;

                    if individual.resistances[b_idx]
                        .iter()
                        .any(|resistance| resistance.any_r > 0.0)
                    {
                        let category = individual
                            .microbiome_resistance_level(b_idx, MICROBIOME_MAJORITY_THRESHOLD);
                        let category_idx = category.as_index();
                        individual.cleared_any_r_microbiome_categories[b_idx][category_idx] += 1;
                    }

                    // If infection was cleared by drugs and bacteria is present in microbiome,
                    // consider clearing it from microbiome as well
                    if matches!(
                        resolution_type,
                        InfectionResolutionType::DrugAssistedClearance
                    ) && individual.presence_microbiome[b_idx]
                    {
                        let microbiome_clearance_on_drug_treatment =
                            get_global_param("microbiome_clearance_probability_on_drug_treatment")
                                .unwrap_or(0.8);
                        if rng.gen_bool(microbiome_clearance_on_drug_treatment) {
                            individual.presence_microbiome[b_idx] = false;
                            individual.microbiome_cleared_today[b_idx] = true;
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
                individual.clearance_hazard[b_idx] = 0.0;
                individual.clearance_ready_day[b_idx] = -1;
                individual.sepsis[b_idx] = false;
                individual.infection_hospital_acquired[b_idx] = false;
                individual.cur_infection_from_environment[b_idx] = false;
                individual.test_identified_infection[b_idx] = false;
                individual.test_for_resistance[b_idx] = false;
                individual.resistance_test_initiated_day[b_idx] = -1;
                individual.infection_has_caused_symptoms[b_idx] = false; // Reset symptom status when infection clears
            } else {
                // Update level for infections that are continuing
                individual.level[b_idx] = new_bacteria_level;
            }
        }

        // Safety check: ensure test_identified_infection and symptom status are false when not infected
        if !is_infected {
            individual.test_identified_infection[b_idx] = false;
            individual.infection_has_caused_symptoms[b_idx] = false;
        }

        // --- Apply cross-resistance logic ---
        apply_cross_resistance(individual, b_idx, cross_resistance_groups);
        // --- END NEW ---

        // Clearance dynamics: arm hazard once infection persists, reset when cleared
        if is_infected {
            if individual.clearance_ready_day[b_idx] == -1 {
                let delay_days = store.clearance.delay_days(b_idx) as i32;
                individual.clearance_ready_day[b_idx] =
                    individual.date_last_infected[b_idx].saturating_add(delay_days.max(0));
            }

            // --- Symptom onset logic for infected bacteria ---
            if !individual.infection_has_caused_symptoms[b_idx] {
                // Get bacteria-specific symptom parameters
                let daily_symptom_probability =
                    store.bacteria.daily_symptom_onset_probability(b_idx);
                let threshold_level = store.bacteria.symptom_onset_threshold_level(b_idx);
                let delay_days = store.bacteria.symptom_onset_delay_days(b_idx) as i32;
                let level_multiplier = store.bacteria.symptom_onset_level_multiplier(b_idx);

                // Check if minimum delay has passed
                let infection_duration = (time_step as i32) - individual.date_last_infected[b_idx];

                if infection_duration >= delay_days && individual.level[b_idx] >= threshold_level {
                    // Calculate symptom onset probability (base rate × level effect)
                    let level_effect =
                        (individual.level[b_idx] / threshold_level).powf(level_multiplier);
                    let symptom_probability =
                        (daily_symptom_probability * level_effect).clamp(0.0, 1.0);

                    // Roll for symptom onset
                    if rng.gen_bool(symptom_probability) {
                        individual.infection_has_caused_symptoms[b_idx] = true;
                    }
                }
            }
        } else {
            individual.clearance_ready_day[b_idx] = -1;
            individual.clearance_hazard[b_idx] = 0.0;
        }
    }

    // Check for post-infection drug usage evaluation (configurable timing)
    let evaluation_days =
        get_global_param("drug_evaluation_days_post_infection").unwrap_or(7.0) as i32;

    for b_idx in 0..BACTERIA_LIST.len() {
        let infection_start_day = individual.date_last_infected_keep[b_idx];

        // Only evaluate if there was an infection and today is exactly the evaluation day after infection start
        if infection_start_day > 0 && (time_step as i32) == (infection_start_day + evaluation_days)
        {
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
            individual.day_7_since_last_infection_drug_used[b_idx] =
                Some(drug_used_since_infection);
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
                if let Some(resistance_data) =
                    individual.resistances.get(b_idx).and_then(|r| r.get(d_idx))
                {
                    if resistance_data.any_r > max_any_r {
                        max_any_r = resistance_data.any_r;
                    }
                }
            }

            // If there's any resistance in the group, update all drugs in the group to the max value
            if max_any_r > 0.0 {
                for &d_idx in group {
                    if let Some(resistance_data) = individual
                        .resistances
                        .get_mut(b_idx)
                        .and_then(|r| r.get_mut(d_idx))
                    {
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
    _param_cache: &ParameterKeyCache,
    is_bacterial_testing: bool,
) -> f64 {
    let store = parameter_store();
    // Get base parameters
    let base_rate = if is_bacterial_testing {
        get_global_param("bacterial_testing_base_rate_per_day").unwrap_or(0.15)
    } else {
        get_global_param("resistance_testing_base_rate_per_day").unwrap_or(0.95)
    };

    // Calculate temporal multiplier (testing adoption over time)
    let years_since_availability = (time_step - testing_available_from_day) as f64 / 365.0;
    let (initial_rate, max_multiplier) = if is_bacterial_testing {
        (
            get_global_param("bacterial_testing_initial_adoption_rate").unwrap_or(0.1),
            get_global_param("bacterial_testing_max_temporal_multiplier").unwrap_or(1.0),
        )
    } else {
        (
            get_global_param("resistance_testing_initial_adoption_rate").unwrap_or(0.05),
            get_global_param("resistance_testing_max_temporal_multiplier").unwrap_or(1.0),
        )
    };

    // Use sigmoid (S-curve) model for more realistic technology adoption
    // Formula: initial_rate + (max_multiplier - initial_rate) * (1 / (1 + e^(-steepness * (years - midpoint))))
    let adoption_years = if is_bacterial_testing { 40.0 } else { 50.0 }; // Years to reach ~95% adoption
    let midpoint = adoption_years / 2.0; // Inflection point (fastest growth)
    let steepness = 6.0 / adoption_years; // Controls how steep the S-curve is

    let sigmoid_factor = 1.0 / (1.0 + (-steepness * (years_since_availability - midpoint)).exp());
    let temporal_multiplier = initial_rate + (max_multiplier - initial_rate) * sigmoid_factor;

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
    let region_multiplier = store.region.testing_multiplier(individual.region_cur_in);

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
    let final_probability = base_rate
        * temporal_multiplier
        * hospital_multiplier
        * region_multiplier
        * immunosuppression_multiplier
        * sepsis_multiplier;

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
        "staphylococcus aureus" => &[
            (2, 0.35),
            (4, 0.25),
            (9, 0.15),
            (3, 0.10),
            (5, 0.08),
            (1, 0.05),
            (6, 0.02),
        ],
        "streptococcus pneumoniae" => &[
            (3, 0.70),
            (6, 0.15),
            (4, 0.08),
            (1, 0.04),
            (2, 0.02),
            (10, 0.01),
        ],
        "streptococcus pyogenes" => &[
            (2, 0.50),
            (3, 0.25),
            (4, 0.15),
            (9, 0.05),
            (5, 0.03),
            (1, 0.02),
        ],
        "streptococcus agalactiae" => &[
            (4, 0.40),
            (6, 0.25),
            (1, 0.15),
            (2, 0.10),
            (3, 0.05),
            (5, 0.05),
        ],
        "enterococcus faecalis" => &[
            (1, 0.50),
            (4, 0.25),
            (5, 0.15),
            (2, 0.05),
            (3, 0.03),
            (9, 0.02),
        ],
        "enterococcus faecium" => &[
            (1, 0.45),
            (4, 0.30),
            (5, 0.15),
            (2, 0.05),
            (3, 0.03),
            (9, 0.02),
        ],

        // Gram-negative Enterobacteriaceae
        "escherichia coli" => &[
            (1, 0.55),
            (4, 0.20),
            (5, 0.12),
            (7, 0.08),
            (2, 0.03),
            (3, 0.02),
        ],
        "klebsiella pneumoniae" => &[
            (3, 0.40),
            (1, 0.25),
            (4, 0.20),
            (5, 0.10),
            (2, 0.03),
            (7, 0.02),
        ],
        "enterobacter spp." => &[
            (1, 0.35),
            (3, 0.25),
            (4, 0.20),
            (5, 0.10),
            (7, 0.05),
            (2, 0.05),
        ],
        "enterobacter_cloacae" => &[
            (4, 0.30),
            (3, 0.25),
            (1, 0.25),
            (5, 0.12),
            (7, 0.05),
            (2, 0.03),
        ],
        "citrobacter spp." => &[
            (1, 0.30),
            (3, 0.25),
            (4, 0.20),
            (5, 0.15),
            (7, 0.05),
            (2, 0.05),
        ],
        "serratia spp." => &[
            (3, 0.35),
            (1, 0.25),
            (4, 0.20),
            (5, 0.10),
            (2, 0.05),
            (7, 0.05),
        ],
        "proteus spp." => &[
            (1, 0.60),
            (4, 0.15),
            (3, 0.10),
            (5, 0.08),
            (2, 0.04),
            (7, 0.03),
        ],
        "morganella spp." => &[
            (1, 0.50),
            (4, 0.20),
            (3, 0.15),
            (5, 0.08),
            (2, 0.04),
            (7, 0.03),
        ],

        // Non-fermenting Gram-negatives
        "pseudomonas aeruginosa" => &[
            (3, 0.45),
            (4, 0.25),
            (1, 0.15),
            (2, 0.08),
            (5, 0.05),
            (9, 0.02),
        ],
        "acinetobacter baumannii" => &[
            (3, 0.40),
            (4, 0.25),
            (1, 0.15),
            (5, 0.10),
            (2, 0.05),
            (7, 0.05),
        ],

        // Gastrointestinal pathogens
        "salmonella enterica serovar typhi" => {
            &[(7, 0.80), (4, 0.15), (5, 0.03), (3, 0.01), (10, 0.01)]
        }
        "salmonella enterica serovar paratyphi a" => {
            &[(7, 0.85), (4, 0.10), (5, 0.03), (3, 0.01), (10, 0.01)]
        }
        "invasive non-typhoidal salmonella spp." => {
            &[(7, 0.70), (4, 0.20), (5, 0.05), (3, 0.03), (1, 0.02)]
        }
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
        "haemophilus influenzae" => &[
            (3, 0.70),
            (6, 0.15),
            (4, 0.08),
            (1, 0.04),
            (2, 0.02),
            (10, 0.01),
        ],
        "moraxella_catarrhalis" => &[(3, 0.85), (4, 0.08), (1, 0.04), (2, 0.02), (10, 0.01)],
        "neisseria_meningitidis" => &[(6, 0.60), (4, 0.25), (3, 0.10), (2, 0.03), (1, 0.02)],
        "bordetella pertussis" => &[(3, 0.95), (6, 0.03), (4, 0.01), (10, 0.01)], // Primarily respiratory (whooping cough)

        // Gastrointestinal pathogens
        "helicobacter pylori" => &[(7, 0.85), (5, 0.10), (4, 0.03), (10, 0.02)], // Primarily GI (peptic ulcer disease)

        // Foodborne/systemic pathogens
        "listeria_monocytogenes" => &[
            (6, 0.50),
            (4, 0.30),
            (7, 0.10),
            (5, 0.05),
            (3, 0.03),
            (1, 0.02),
        ],

        // Fallback for any unmatched bacteria (should not occur with complete list above)
        _ => &[
            (1, 0.1),
            (2, 0.1),
            (3, 0.1),
            (4, 0.1),
            (5, 0.1),
            (6, 0.1),
            (7, 0.1),
            (8, 0.1),
            (9, 0.1),
            (10, 0.1),
        ],
    };

    let weights: Vec<f64> = syndrome_probs.iter().map(|&(_, p)| p).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    syndrome_probs[dist.sample(rng)].0
}
