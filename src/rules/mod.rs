// src/rules/mod.rs

use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES, HospitalStatus, Region}; 
use crate::config::{get_global_param, get_bacteria_param, get_drug_param, get_age_infection_multiplier, get_drug_availability, get_bacteria_sepsis_risk_multiplier};
use rand::Rng;
use rand::seq::SliceRandom;
use std::collections::HashMap;
use rand::distributions::WeightedIndex;
use rand::distributions::Distribution; 

/// applies model rules to an individual for one time step.
pub fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    _global_majority_r_proportions: &HashMap<(usize, usize), f64>,
    majority_r_positive_values_by_combo: &HashMap<(usize, bool, usize, usize), Vec<f64>>, // <-- update type
    bacteria_indices: &HashMap<&'static str, usize>,
    drug_indices: &HashMap<&'static str, usize>,
    cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>, // New parameter
) {

    if individual.age < 0 {
        individual.age += 1; // Only advance age by 1 day
        return; // Exit the function if unborn
    }

    if individual.date_of_death.is_some() {
        return; // Exit the function if dead
    }

    if individual.id == 0  { 
        println!("   "); println!("mod.rs time step {}", time_step); println!("   "); 
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
    // get general parameters for fluctuations and bounds
    let daily_fluctuation = get_global_param("contact_level_daily_fluctuation_range").unwrap_or(0.5);
    let min_contact_level = get_global_param("min_contact_level").unwrap_or(0.0);
    let max_contact_level = get_global_param("max_contact_level").unwrap_or(10.0);

    // helper closure for applying fluctuation and clamping
    // this calculates a 'target' or 'base' level, then adds noise and clamps it.
    let mut update_contact_level = |current_level: &mut f64, base_level: f64| {
        *current_level = base_level + rng.gen_range(-daily_fluctuation..=daily_fluctuation);
        *current_level = current_level.clamp(min_contact_level, max_contact_level);
    };

    //  sexual contact level
    let sexual_contact_age_peak_days = get_global_param("sexual_contact_age_peak_days").unwrap_or(25.0 * 365.0);
    let sexual_contact_age_decline_rate = get_global_param("sexual_contact_age_decline_rate").unwrap_or(0.00005);
    let sexual_contact_hospital_multiplier = get_global_param("sexual_contact_hospital_multiplier").unwrap_or(0.0); // Typically very low in hospital

    let mut base_sexual_level = get_global_param("sexual_contact_baseline").unwrap_or(5.0);
    if (individual.age as f64) < sexual_contact_age_peak_days {
    // Increase towards peak, but don't exceed baseline before peak
        base_sexual_level *= ((individual.age as f64 / sexual_contact_age_peak_days).min(1.0)).powf(get_global_param("sexual_contact_age_rise_exponent").unwrap_or(2.0));
    } else {
    // decline after peak
        base_sexual_level *= (1.0 - (individual.age as f64 - sexual_contact_age_peak_days) * sexual_contact_age_decline_rate).max(0.0);
    }

    if individual.hospital_status.is_hospitalized() {
        base_sexual_level *= sexual_contact_hospital_multiplier;
    }
    update_contact_level(&mut individual.sexual_contact_level, base_sexual_level);


    // airborne contact level with adults
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


    // Airborne Contact Level with Children
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


    // Oral Exposure Level
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


    // Mosquito Exposure Level
    let mosquito_exposure_baseline = get_global_param("mosquito_exposure_baseline").unwrap_or(1.0);
    let mosquito_exposure_in_hospital_multiplier = get_global_param("mosquito_exposure_in_hospital_multiplier").unwrap_or(0.2); // Usually lower indoors/hospitals

    let mut base_mosquito_level = mosquito_exposure_baseline;

    // Apply region-specific multiplier
    let region_name_for_param = individual.region_cur_in.to_string().to_lowercase().replace(" ", "_");
    let region_multiplier_key = format!("{}_mosquito_exposure_multiplier", region_name_for_param);
    let region_multiplier = get_global_param(&region_multiplier_key).unwrap_or(1.0); // Default to 1.0 if region not specified
    base_mosquito_level *= region_multiplier;

    if individual.hospital_status.is_hospitalized() {
        base_mosquito_level *= mosquito_exposure_in_hospital_multiplier;
    }
    update_contact_level(&mut individual.mosquito_exposure_level, base_mosquito_level);

    // --- end update contact and exposure levels ---



    //  update 'is_severely_immunosuppressed' status based on onset/recovery rates
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


    // current toxicity
    individual.current_toxicity = (individual.current_toxicity + rng.gen_range(-0.5..=0.5)).max(0.0);


    // Get parameters from config.rs once per individual for this time step
    // todo: review this update rule
    let baseline_rate = get_global_param("hospitalization_baseline_rate_per_day")
        .expect("Missing hospitalization_baseline_rate_per_day in config");
    let age_multiplier_hosp = get_global_param("hospitalization_age_multiplier_per_day")
        .expect("Missing hospitalization_age_multiplier_per_day in config");
    let recovery_rate = get_global_param("hospitalization_recovery_rate_per_day")
        .expect("Missing hospitalization_recovery_rate_per_day in config");
    let max_days_in_hospital = get_global_param("hospitalization_max_days")
        .expect("Missing hospitalization_max_days in config");

    // Potentially get hospitalized (if not currently hospitalized)
    if !individual.hospital_status.is_hospitalized() { 
        let prob_hospitalization_today = baseline_rate + (individual.age as f64 * age_multiplier_hosp);

        if rng.gen::<f64>() < prob_hospitalization_today {
            individual.hospital_status = HospitalStatus::InHospital; 
            individual.days_hospitalized = 0; // Initialize days hospitalized
        }
    } else { // If already hospitalized, consider recovery or max days limit
        individual.days_hospitalized += 1; // Increment days hospitalized

        // Potentially recover from hospitalization
        if rng.gen::<f64>() < recovery_rate {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
            // println!("individual {} recovered from hospitalization.", individual.id);             
        }
        // discharge after max_days_in_hospital
        else if individual.days_hospitalized >= max_days_in_hospital as u32 {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
         }
    }
    // --- end hospitalization Rules ---



    // ---  region travel ---
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
            let last_infected_day = individual.date_last_infected[b_idx];
            let duration_of_infection = (time_step as i32 - last_infected_day).max(0); // Ensure non-negative duration

            // Retrieve bacteria-specific parameters, falling back to global defaults
            let sepsis_baseline_risk = get_bacteria_param(bacteria, "sepsis_baseline_risk_per_day")
                .unwrap_or_else(|| get_global_param("default_sepsis_baseline_risk_per_day").expect("Missing default_sepsis_baseline_risk_per_day"));
            let sepsis_level_multiplier = get_bacteria_param(bacteria, "sepsis_level_multiplier")
                .unwrap_or_else(|| get_global_param("default_sepsis_level_multiplier").expect("Missing default_sepsis_level_multiplier"));
            let sepsis_duration_multiplier = get_bacteria_param(bacteria, "sepsis_duration_multiplier")
                .unwrap_or_else(|| get_global_param("default_sepsis_duration_multiplier").expect("Missing default_sepsis_duration_multiplier"));

            // Get bacteria-specific sepsis risk category multiplier
            let bacteria_sepsis_multiplier = get_bacteria_sepsis_risk_multiplier(bacteria);
            
            // Calculate daily probability of sepsis with bacteria-specific risk category
            let prob_sepsis_today = (sepsis_baseline_risk
                                    + (current_level * sepsis_level_multiplier)
                                    + (duration_of_infection as f64 * sepsis_duration_multiplier))
                                    * bacteria_sepsis_multiplier;

            // Cap the probability at 1.0
            let prob_sepsis_today = prob_sepsis_today.min(1.0);

            if rng.gen::<f64>() < prob_sepsis_today {
                // Set sepsis status to true for this bacteria
                individual.sepsis[b_idx] = true;
            }
        } else {
            if individual.sepsis[b_idx] {
                individual.sepsis[b_idx] = false;
            }
        }
    }
    // --- end sepsis updates ---




    // loop through all bacteria to update vaccination status dynamically
    for (b_idx, _bacteria) in BACTERIA_LIST.iter().enumerate() {
        if rng.gen::<f64>() < 0.0001 {
            individual.vaccination_status[b_idx] = !individual.vaccination_status[b_idx];
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
            let param_name = format!("syndrome_{}_initiation_multiplier", syndrome_id);
            if let Some(multiplier) = get_global_param(&param_name) {
                syndrome_administration_multiplier = syndrome_administration_multiplier.max(multiplier);
            }
        }
    }

    let mut drugs_initiated_this_time_step: usize = 0;

    // --- drug stopping ---
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];
        if individual.cur_use_drug[drug_idx] {
            let mut relevant_infection_active_for_this_drug = false;
            for b_idx in 0..BACTERIA_LIST.len() {
                if individual.level[b_idx] > 0.0001 {
                    let bacteria_name = BACTERIA_LIST[b_idx];
                    // Use potency_when_no_r to determine if drug is relevant for this bacteria
                    let potency_param_key = format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug_name, bacteria_name);
                    let drug_potency = get_global_param(&potency_param_key).unwrap_or(0.0);
                    if drug_potency > 0.0 {
                        relevant_infection_active_for_this_drug = true;
                        break;
                    }
                }
            }
            let mut stop_drug = false;
            // todo: maybe re-work this so having no relavant infection is a very strong predictor of stopping but not 
            // something that 100% causes stopping
            if !relevant_infection_active_for_this_drug || rng.gen_bool(random_drug_cessation_prob) {
                stop_drug = true;
            }
            if individual.date_drug_initiated[drug_idx] == (time_step as i32) - 1 {
                stop_drug = false;
            }
            if stop_drug {
                individual.cur_use_drug[drug_idx] = false;
                individual.date_drug_initiated[drug_idx] = i32::MIN;
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

    // --- drug initiation ---
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];

        // --- restriction: if already using two or more drugs, cannot start another ---
        if num_drugs_currently_used + drugs_initiated_this_time_step >= 2 {
            continue;
        }

        // --- restriction: do not start drug if test_r > 0 for any bacteria ---
        let mut test_r_positive = false;
        for b_idx in 0..BACTERIA_LIST.len() {
            if individual.resistances[b_idx][drug_idx].test_r > 0.0 {
                test_r_positive = true;
                break;
            }
        }
        if test_r_positive {
            continue;
        }
        // --- end restriction ---

        // --- NEW: Check drug availability in current region ---
        let current_region_str = individual.region_cur_in.to_string();
        let living_region_str = individual.region_living.to_string();
        let drug_availability = get_drug_availability(
            drug_name, 
            &current_region_str, 
            Some(&living_region_str)
        );
        
        // If drug is not available or has very low availability, skip it
        if drug_availability < 0.01 {
            continue; // Drug not available in this region
        }
        // --- end drug availability check ---

        // start with the base initiation rate for *any* drug
        let mut administration_prob = drug_base_initiation_rate;

        // --- apply bacteria-specific multipliers if the individual has active infections ---
        let mut max_bacteria_specific_multiplier: f64 = 1.0; // Use max to represent the highest relevance
        for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
            // check if individual is infected with this specific bacteria and it's above threshold
            if individual.level[b_idx] > 0.001 {
                // Look up the specific multiplier for this drug-bacteria combination
                let param_key = format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug_name, bacteria_name);
                if let Some(specific_multiplier) = get_global_param(&param_key) {
                    max_bacteria_specific_multiplier = max_bacteria_specific_multiplier.max(specific_multiplier);
                }
            }
        }
        // apply the highest relevant bacteria-specific multiplier
        administration_prob *= max_bacteria_specific_multiplier;


        // apply general multipliers after the base rate and bacteria-specific multiplier
        let infection_acquired_this_step = individual.date_last_infected.iter().any(|&d| d == time_step as i32);
        if has_any_infection && !infection_acquired_this_step {
        administration_prob *= drug_infection_present_multiplier;
        }
       
        if has_any_identified_infection { administration_prob *= drug_test_identified_multiplier; }
        if initial_on_any_antibiotic || drugs_initiated_this_time_step > 0 {
            administration_prob *= already_on_drug_initiation_multiplier;
        }
        administration_prob *= syndrome_administration_multiplier;

        // --- NEW: Apply bacterial identification effects on drug spectrum preference ---
        let drug_spectrum = get_drug_param(drug_name, "spectrum_breadth").unwrap_or(3.0); // 1.0=narrow, 5.0=very broad
        
        if has_any_identified_infection {
            // TARGETED THERAPY: bacteria identified, prefer appropriate narrow-spectrum drugs
            let targeted_narrow_bonus = get_global_param("targeted_therapy_narrow_spectrum_bonus").unwrap_or(3.0);
            let targeted_broad_penalty = get_global_param("targeted_therapy_broad_spectrum_penalty").unwrap_or(0.4);
            let ineffective_drug_penalty = get_global_param("targeted_therapy_ineffective_drug_penalty").unwrap_or(0.1);
            
            // Check if this drug has good activity against any identified bacteria
            let mut has_good_activity = false;
            let mut best_potency: f64 = 0.0;
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                if individual.test_identified_infection[b_idx] && individual.level[b_idx] > 0.001 {
                    let potency_param_key = format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug_name, bacteria_name);
                    let potency = get_global_param(&potency_param_key).unwrap_or(0.0);
                    best_potency = best_potency.max(potency);
                    if potency > 0.02 { // Threshold for "good activity" (above baseline)
                        has_good_activity = true;
                    }
                }
            }
            
            if has_good_activity {
                // Drug is effective against identified bacteria - prefer narrow spectrum
                if drug_spectrum <= 2.5 { // Narrow to medium-narrow spectrum
                    administration_prob *= targeted_narrow_bonus;
                } else if drug_spectrum >= 4.0 { // Broad to very broad spectrum
                    administration_prob *= targeted_broad_penalty;
                }
                // Medium spectrum (2.5-4.0) gets no bonus or penalty
            } else {
                // Drug is not effective against identified bacteria - strong penalty
                administration_prob *= ineffective_drug_penalty;
            }
        } else if has_any_infection {
            // EMPIRIC THERAPY: infection present but bacteria not yet identified, prefer broad-spectrum
            let empiric_broad_bonus = get_global_param("empiric_therapy_broad_spectrum_bonus").unwrap_or(2.0);
            
            if drug_spectrum >= 3.5 { // Broad to very broad spectrum drugs
                administration_prob *= empiric_broad_bonus;
            } else if drug_spectrum <= 2.0 { // Very narrow spectrum drugs
                administration_prob *= 0.6; // Slight penalty for very narrow drugs in empiric therapy
            }
            // Medium spectrum (2.0-3.5) gets no bonus or penalty
        }
        // --- END bacterial identification effects ---

        // --- Apply drug availability multiplier ---
        administration_prob *= drug_availability;
        // --- end drug availability ---

        administration_prob = administration_prob.clamp(0.0, 1.0); // Ensure the probability is between 0 and 1
        if drugs_initiated_this_time_step < 2 && !individual.cur_use_drug[drug_idx] {
            if rng.gen_bool(administration_prob) {
                individual.cur_use_drug[drug_idx] = true;
                individual.date_drug_initiated[drug_idx] = time_step as i32;
                individual.ever_taken_drug[drug_idx] = true;

                // debug print       
                                  
                if individual.id == 0  { 
                    println!(
                        "mod.rs   started {} - rate of starting was {:.4}",
                        drug_name,
                        administration_prob,
                    );
                }
                // --- end debug print

                let mut chosen_initial_level = get_drug_param(drug_name, "initial_level").unwrap_or(10.0);
                if has_any_identified_infection && rng.gen_bool(double_dose_probability) {
                    let double_dose_multiplier = get_drug_param(drug_name, "double_dose_multiplier").unwrap_or(2.0);
                    chosen_initial_level *= double_dose_multiplier;
                }
                individual.cur_level_drug[drug_idx] = chosen_initial_level;
                drugs_initiated_this_time_step += 1;
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

    // --- death     

    // todo: review this update rule
    // e.g. death rate by age will be higher at very young age in africa

    if individual.date_of_death.is_none() {
        let mut cause: Option<String> = None;
        let base_background_rate = get_global_param("base_background_mortality_rate_per_day")
            .expect("Missing base_background_mortality_rate_per_day in config");
        let age_multiplier = get_global_param("age_mortality_multiplier_per_year")
            .expect("Missing age_mortality_multiplier_per_year in config");
        let mut background_risk = base_background_rate;
        background_risk *= (individual.age as f64 / 365.0) * age_multiplier;
        let region_multiplier_key = format!("{}_mortality_multiplier", individual.region_living.to_string().to_lowercase().replace(" ", "_"));
        let region_multiplier = get_global_param(&region_multiplier_key).unwrap_or(1.0);
        background_risk *= region_multiplier;
        let sex_multiplier_key = format!("{}_mortality_multiplier", individual.sex_at_birth.to_lowercase());
        let sex_multiplier = get_global_param(&sex_multiplier_key).unwrap_or(1.0);
        background_risk *= sex_multiplier;
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
            let region_sepsis_multiplier_key = format!("{}_sepsis_mortality_multiplier", 
                individual.region_living.to_string().to_lowercase().replace(" ", "_"));
            let region_sepsis_multiplier = get_global_param(&region_sepsis_multiplier_key).unwrap_or(1.0);
            sepsis_death_risk *= region_sepsis_multiplier;
            
            // Apply immunosuppression multiplier
            if individual.is_severely_immunosuppressed {
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
            let drug_name = DRUG_SHORT_NAMES[drug_idx];
            if individual.cur_level_drug[drug_idx] > 0.0 {
                let drug_adverse_event_risk = get_drug_param(drug_name, "adverse_event_death_risk")
                    .unwrap_or(0.0);
                drug_adverse_event_risk_for_individual = (drug_adverse_event_risk_for_individual + drug_adverse_event_risk).min(1.0);
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
        }
    }
    // --- death logic end   

    // --- update per-bacteria fields ---
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        let is_infected = individual.level[b_idx] > 0.001;

        if !is_infected {
            // --- bacteria-specific acquisition probability ---
            let mut acquisition_probability = get_bacteria_param(bacteria, "acquisition_prob_baseline").unwrap_or(0.01);

            // apply contact level modifiers dynamically
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

            // apply vaccination status effect dynamically
            if individual.vaccination_status[b_idx] {
                let vaccine_efficacy = get_bacteria_param(bacteria, "vaccine_efficacy").unwrap_or(0.0);
                acquisition_probability *= 1.0 - vaccine_efficacy;
            }

            // microbiome presence effect
            if individual.presence_microbiome[b_idx] {
                let microbiome_infection_multiplier = get_bacteria_param(bacteria, "microbiome_infection_acquisition_multiplier")
                    .unwrap_or_else(|| get_global_param("default_microbiome_infection_acquisition_multiplier").expect("Missing default_microbiome_infection_acquisition_multiplier in config"));
                acquisition_probability *= microbiome_infection_multiplier;
            }

            // hospital-acquired multiplier (only if in hospital)
            if individual.hospital_status.is_hospitalized() {
                let hospital_multiplier = get_bacteria_param(bacteria, "hospital_acquired_multiplier").unwrap_or(1.0);
                acquisition_probability *= hospital_multiplier;
            }

            // age-based infection risk multiplier
            let age_multiplier = get_age_infection_multiplier(bacteria, individual.age);
            acquisition_probability *= age_multiplier;

            // region-specific bacterial infection risk multiplier
            let region_name_for_param = individual.region_cur_in.to_string().to_lowercase().replace(" ", "_");
            let region_bacteria_multiplier_key = format!("{}_{}_infection_risk_multiplier", region_name_for_param, bacteria.replace(" ", "_"));
            let region_bacteria_multiplier = get_global_param(&region_bacteria_multiplier_key)
                .unwrap_or_else(|| {
                    // If specific region-bacteria combination not found, try default for this region
                    let default_region_key = format!("{}_infection_risk_multiplier_default", region_name_for_param);
                    get_global_param(&default_region_key).unwrap_or(1.0)  // Default to 1.0 if no region-specific multiplier found
                });
            acquisition_probability *= region_bacteria_multiplier;

            // --- microbiome presence (Carriage) ---
            if !individual.presence_microbiome[b_idx] {
                let microbiome_acquisition_multiplier = get_bacteria_param(bacteria, "microbiome_acquisition_multiplier")
                    .unwrap_or_else(|| get_global_param("default_microbiome_acquisition_multiplier").expect("Missing default_microbiome_acquisition_multiplier in config"));
                let microbiome_acquisition_probability = acquisition_probability * microbiome_acquisition_multiplier;
                if rng.gen_bool(microbiome_acquisition_probability.clamp(0.0, 1.0)) {
                    individual.presence_microbiome[b_idx] = true;

                    // --- assign microbiome_r on new microbiome acquisition (same logic as infection resistance assignment) ---
                    let env_majority_r_level = get_global_param("environmental_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                    let hospital_majority_r_level = get_global_param("hospital_majority_r_level_for_new_acquisition").unwrap_or(0.0);
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
                        } else if is_hospital_acquired {
                            resistance_data.microbiome_r = hospital_majority_r_level;
                        } else {
                            if let Some(majority_r_values_from_population) =
                                majority_r_positive_values_by_combo.get(&(region_idx, hospital_status_bool, b_idx, d_idx))
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
                                .or_else(|| get_global_param("resistance_emergence_rate_per_day_baseline"))
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

            if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                let initial_level = get_bacteria_param(bacteria, "initial_infection_level").unwrap_or(0.01);
                individual.level[b_idx] = initial_level;
                individual.date_last_infected[b_idx] = time_step as i32;

                // --- probabilistic syndrome assignment ---
                let syndrome_id = assign_syndrome_for_bacteria(bacteria, &mut rng);
                individual.infectious_syndrome[b_idx] = syndrome_id as i32;

                let env_acquisition_chance = get_bacteria_param(bacteria, "environmental_acquisition_proportion").unwrap_or(0.1);
                individual.cur_infection_from_environment[b_idx] = rng.gen::<f64>() < env_acquisition_chance;

                individual.infection_hospital_acquired[b_idx] = individual.hospital_status.is_hospitalized();

                // --- any_r and majority_r setting logic on new infection acquisition ---
                // todo: have the posisbility of any_r also for new micribione acquisition of bacteria
                let env_majority_r_level = get_global_param("environmental_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let hospital_majority_r_level = get_global_param("hospital_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0);

                //  todo: drug treatment leads to increase in risk of microbiome_r > 0 (due to allowing more bacteria growth due to killing
                //  other bacteria in microbiome (so can be caused by any drug) or direct selection of resistance to the drug veing taken 
                //  (and those with cross resistance) as occurs for an infection) 


                let is_from_environment = individual.cur_infection_from_environment[b_idx];
                let is_hospital_acquired = individual.infection_hospital_acquired[b_idx];

                let region_idx = individual.region_cur_in as usize;
                let hospital_status_bool = individual.hospital_status.is_hospitalized();

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
                        // --- region/hospital-specific sampling ---
                        if let Some(majority_r_values_from_population) =
                            majority_r_positive_values_by_combo.get(&(region_idx, hospital_status_bool, b_idx, d_idx))
                        {
                            if let Some(&acquired_resistance_level) = majority_r_values_from_population.choose(&mut rng) {
                                let clamped_level = acquired_resistance_level.min(max_resistance_level).max(0.0);
                                resistance_data.any_r = clamped_level;
                                resistance_data.majority_r = clamped_level;
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
            } 
        } else { // Bacteria is already present (infection progression)
            // --- majority_r evolution ---
            let majority_r_evolution_rate = get_global_param("majority_r_evolution_rate_per_day_when_drug_present").unwrap_or(0.0);
            let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0); // Now using 1.0 from your config

            if let Some(bacteria_full_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
                for (drug_index, _use_drug) in individual.cur_use_drug.iter().enumerate() { 
                    let resistance_data = &mut individual.resistances[bacteria_full_idx][drug_index];

                    let drug_current_level = individual.cur_level_drug[drug_index];
                    let drug_currently_present = drug_current_level > 0.0001; // Check if drug is effectively present
                    let current_bacteria_level = individual.level[b_idx];

                    // existing majority_r evolution based on drug presence
                    if resistance_data.majority_r == 0.0 && resistance_data.any_r > 0.0 && drug_currently_present {
                        if rng.gen_bool(majority_r_evolution_rate) {
                            resistance_data.majority_r = resistance_data.any_r;
                        }
                    }

                    // todo: check: value for any_r or majority_r for any drug bacteria combination should 
                    // not decline so long as the bacterial infection is present - even after bacterial infection
                    // has gone it may be in microbiome      

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
                        if drug_current_level > 0.0001 && current_bacteria_level > 0.0001 { 
                            let param_key = format!(
                                "drug_{}_for_bacteria_{}_resistance_emergence_rate_per_day_baseline",
                                DRUG_SHORT_NAMES[drug_index],
                                bacteria
                            );
                            let emergence_rate_baseline = get_global_param(&param_key).unwrap_or(0.000001); // Very small baseline
                            let bacteria_level_effect_multiplier = get_global_param("resistance_emergence_bacteria_level_multiplier").unwrap_or(0.05); // How much does bacteria level boost it
                            let any_r_emergence_level_on_first_emergence = get_global_param("any_r_emergence_level_on_first_emergence").unwrap_or(0.5); // User changed to 0.5 (was 1.0)

                            // bacteria level dependency: Higher at higher levels
                            let max_bacteria_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);
                            // Normalize bacteria level to [0,1] and apply multiplier
                            let bacteria_level_factor = (current_bacteria_level / max_bacteria_level).clamp(0.0, 1.0) * bacteria_level_effect_multiplier;
                            
                            // activity_r dependency: Bell-shaped curve
                            // Use the drug's initial level for normalization to get a comparable 'activity' scale (0-10)
                            let drug_initial_level_for_normalization = get_drug_param(DRUG_SHORT_NAMES[drug_index], "initial_level").unwrap_or(10.0);
                            
                            // normalized current drug level as a proxy for 'activity_r' when any_r is 0.
                            let mut norm_drug_level = drug_current_level / drug_initial_level_for_normalization;
                            norm_drug_level = norm_drug_level.clamp(0.0, 10.0); 
                            
                            // todo: review this code for resistance emergence probability
                            // bell-shaped curve: 0.02 * x * (10 - x). Peaks at 5.0, is 0.1 at 0 and 10.
                            let activity_r_bell_curve_factor = 0.1 + 0.02 * norm_drug_level * (10.0 - norm_drug_level);
                            let final_activity_r_factor = activity_r_bell_curve_factor.clamp(0.0, 1.0);  



                            if individual.id == 0 {
                                println!(" ");
                                println!("mod.rs");  
                                println!("final_activity_r_factor: {:.4}", final_activity_r_factor);
                            }



                            // total emergence probability
                            // adding 1.0 to bacteria_level_factor ensures a base contribution even if multiplier is low
                            let total_emergence_prob = emergence_rate_baseline * (1.0 + bacteria_level_factor) * final_activity_r_factor;

                            if rng.gen_bool(total_emergence_prob.clamp(0.0, 1.0)) {
                                resistance_data.any_r = any_r_emergence_level_on_first_emergence;
                            }
                        }
                    }
                    // --- end new resistance emergence logic ---


                    // calculate activity_r (should always be updated)
                    // todo: may need to specify the parameter 0.05 below in config.rs
                    if drug_current_level > 0.0 {
                        // Fetch potency from config, fallback to 0.05 if not found
                        let potency_param_key = format!(
                            "drug_{}_for_bacteria_{}_potency_when_no_r",
                            DRUG_SHORT_NAMES[drug_index], bacteria
                        );
                        let potency = get_global_param(&potency_param_key).unwrap_or(0.05);
                        let normalized_any_r = resistance_data.any_r / max_resistance_level;
                        resistance_data.activity_r = potency * drug_current_level * (1.0 - normalized_any_r);
   
                    } else {
                        resistance_data.activity_r = 0.0;
                    }
                }
            }
        }

        // testing and diagnosis
        let last_infected_time = individual.date_last_infected[b_idx];
        let current_test_status_entry = &mut individual.test_identified_infection[b_idx];
        let test_delay_days = get_global_param("test_delay_days").unwrap_or(3.0) as i32;
        let test_rate_per_day = get_global_param("test_rate_per_day").unwrap_or(0.15);
        if !*current_test_status_entry && (time_step as i32) >= (last_infected_time + test_delay_days) {
            if rng.gen_bool(test_rate_per_day.clamp(0.0, 1.0)) {
                *current_test_status_entry = true;
            }
        }

        // --- test_r assignment logic ---
        let prob_test_r_done = get_global_param("prob_test_r_done").unwrap_or(0.95);
        let test_r_error_prob = get_global_param("test_r_error_probability").unwrap_or(0.02);
        let test_r_error_value = get_global_param("test_r_error_value").unwrap_or(0.25);

        if *current_test_status_entry {
            let test_r_already_set = individual.resistances[b_idx].iter().any(|r| r.test_r > 0.0);
            if !test_r_already_set {
                if rng.gen_bool(prob_test_r_done) {
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


            if individual.id == 0 {
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


                if individual.id == 0 {
                        println!(
                            "mod.rs  {}: current level = {:.4}, activity_r = {:.4}",
                            DRUG_SHORT_NAMES[drug_idx],
                            individual.cur_level_drug[drug_idx],
                            resistance_data.activity_r
                        );
                    }
                        

             if individual.id == 0 {
                println!("mod.rs  total reduction due to antibiotic: {:.4}", total_reduction_due_to_antibiotic);
            }   
            }
            }
            let decay = baseline_change - (immunity_level * reduction_due_to_immune_resp) - total_reduction_due_to_antibiotic;

            let max_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);
            let new_level = (individual.level[b_idx] + decay).max(0.0).min(max_level);

   
                if individual.id == 0 {

                println!(" "); 
                println!("mod.rs");                    
                println!(" ");  
                println!("bacteria level after previous time step: {:.4}", individual.level[b_idx]);
                println!("bacteria level after this time step: {:.4}", new_level);
                println!("calculated change: {:.4}", decay);

                }


            individual.level[b_idx] = new_level;
        } 

        if individual.level[b_idx] < 0.0001 {
            for drug_idx_clear in 0..DRUG_SHORT_NAMES.len() {
                let resistance_data = &mut individual.resistances[b_idx][drug_idx_clear];
                resistance_data.any_r = 0.0;
                resistance_data.majority_r = 0.0;
                resistance_data.activity_r = 0.0;
            }
            individual.level[b_idx] = 0.0;
            individual.infectious_syndrome[b_idx] = 0;
            individual.date_last_infected[b_idx] = 0;
            individual.immune_resp[b_idx] = 0.0;
            individual.sepsis[b_idx] = false;
            individual.presence_microbiome[b_idx] = false;
            individual.infection_hospital_acquired[b_idx] = false;
            individual.cur_infection_from_environment[b_idx] = false;
            individual.test_identified_infection[b_idx] = false;
        }

        // --- NEW: Apply cross-resistance logic ---
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
            if individual.is_severely_immunosuppressed {
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

/// Helper function to probabilistically assign a syndrome for a given bacteria.
fn assign_syndrome_for_bacteria<R: Rng>(bacteria: &str, rng: &mut R) -> u32 {
    // Define syndrome probabilities for each bacteria.
    // Each entry: (syndrome_id, probability)
    let syndrome_probs: &[(u32, f64)] = match bacteria {
        "strep_pneu" => &[(3, 0.9), (7, 0.1)], // 90% respiratory, 10% GI
        "haem_infl" => &[(3, 1.0)],
        "kleb_pneu" => &[(3, 0.8), (7, 0.2)],
        "salm_typhi" => &[(7, 1.0)],
        "salm_parat_a" => &[(7, 1.0)],
        "inv_nt_salm" => &[(7, 1.0)],
        "shig_spec" => &[(7, 1.0)],
        "esch_coli" => &[(7, 0.7), (8, 0.3)],
        "n_gonorrhoeae" => &[(8, 1.0)],
        "group_a_strep" => &[(9, 1.0)],
        "group_b_strep" => &[(10, 1.0)],
        // Add more bacteria as needed
        _ => &[(1, 0.1), (2, 0.1), (3, 0.1), (4, 0.1), (5, 0.1), (6, 0.1), (7, 0.1), (8, 0.1), (9, 0.1), (10, 0.1)],
    };

    let weights: Vec<f64> = syndrome_probs.iter().map(|&(_, p)| p).collect();
    let dist = WeightedIndex::new(&weights).unwrap();
    syndrome_probs[dist.sample(rng)].0
}