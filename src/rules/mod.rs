// src/rules/mod.rs

use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES, HospitalStatus, Region}; 
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

    // --- Move all these parameter lookups to the top so they're in scope everywhere ---
    let transfer_prob = get_global_param("microbiome_resistance_transfer_probability_per_day").unwrap_or(0.05);
    let drug_base_initiation_rate = get_global_param("drug_base_initiation_rate_per_day").unwrap_or(0.0001);
    let drug_infection_present_multiplier = get_global_param("drug_infection_present_multiplier").unwrap_or(50.0);
    let already_on_drug_initiation_multiplier = get_global_param("already_on_drug_initiation_multiplier").unwrap_or(0.0001);
    let drug_test_identified_multiplier = get_global_param("drug_test_identified_multiplier").unwrap_or(20.0);
    let drug_decay_rate = get_global_param("drug_decay_rate_per_day").unwrap_or(0.3);
    let double_dose_probability = get_global_param("double_dose_probability_if_identified_infection").unwrap_or(0.1);
    let random_drug_cessation_prob = get_global_param("random_drug_cessation_probability").unwrap_or(0.001);

    // Update non-infection, bacteria or antibiotic-specific variables
    // need a variable for vulnerability to serious toxicity ?
    individual.age += 1;
//  individual.current_infection_related_death_risk = rng.gen_range(0.0..=0.001);
//  individual.background_all_cause_mortality_rate = rng.gen_range(0.0..=0.00001);






    // ---  Update Contact and Exposure Levels ---
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





       //  Update 'is_severely_immunosuppressed' status based on onset/recovery rates
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





    individual.current_toxicity = (individual.current_toxicity + rng.gen_range(-0.5..=0.5)).max(0.0);
//  individual.mortality_risk_current_toxicity = rng.gen_range(0.0..=0.0001);





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



       // ---  Region Travel Rules ---
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




        // ---  Sepsis Risk Rules ---
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

            // Calculate daily probability of sepsis
            let prob_sepsis_today = sepsis_baseline_risk
                                    + (current_level * sepsis_level_multiplier)
                                    + (duration_of_infection as f64 * sepsis_duration_multiplier);

            // Cap the probability at 1.0
            let prob_sepsis_today = prob_sepsis_today.min(1.0);

            // Roll the dice for sepsis
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
    // --- END Sepsis Risk Rules ---




    // Loop through all bacteria to update vaccination status dynamically
    for (b_idx, _bacteria) in BACTERIA_LIST.iter().enumerate() {
        if rng.gen::<f64>() < 0.0001 {
            individual.vaccination_status[b_idx] = !individual.vaccination_status[b_idx];
        }
    }

    // --- DRUG LOGIC START ---
    let has_any_infection = individual.level.iter().any(|&level| level > 0.0);
    let initial_on_any_antibiotic = individual.cur_use_drug.iter().any(|&identified| identified);
    let has_any_identified_infection = individual.test_identified_infection.iter().any(|&identified| identified);

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

    // --- Drug Stopping Logic ---
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];
        if individual.cur_use_drug[drug_idx] {
            let mut relevant_infection_active_for_this_drug = false;
            for (b_idx, _bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                if individual.level[b_idx] > 0.0001 {
                    let bacteria_name = BACTERIA_LIST[b_idx];
                    let drug_reduction_efficacy = get_bacteria_drug_param(
                        bacteria_name,
                        drug_name,
                        "bacteria_level_reduction_per_unit_of_drug",
                    ).unwrap_or(0.0);
                    if drug_reduction_efficacy > 0.0 {
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

    // Apply decay if stopped, or set to initial level if continued/re-initiated.
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];
        let drug_initial_level = get_drug_param(drug_name, "initial_level").unwrap_or(10.0);
        if individual.cur_use_drug[drug_idx] {
            individual.cur_level_drug[drug_idx] = drug_initial_level;
        } else {
            individual.cur_level_drug[drug_idx] = (individual.cur_level_drug[drug_idx] - drug_decay_rate).max(0.0);
        }
    }

// --- MODIFIED DRUG INITIATION LOGIC ---
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_name = DRUG_SHORT_NAMES[drug_idx];

        // Start with the base initiation rate for *any* drug
        let mut administration_prob = drug_base_initiation_rate;

        // --- Apply bacteria-specific multipliers if the individual has active infections ---
        let mut max_bacteria_specific_multiplier: f64 = 1.0; // Use max to represent the highest relevance
        for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
            // Check if individual is infected with this specific bacteria and it's above threshold
            if individual.level[b_idx] > 0.001 {
                // Look up the specific multiplier for this drug-bacteria combination
                let param_key = format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug_name, bacteria_name);
                if let Some(specific_multiplier) = get_global_param(&param_key) {
                    // We want the *highest* relevant multiplier to apply, or you could sum them.
                    // Summing them would mean if a drug treats multiple active infections, it's even more likely.
                    // Max: administration_prob *= specific_multiplier; (this applies the max one)
                    // Sum: administration_prob += drug_base_initiation_rate * (specific_multiplier - 1.0); (this adds the 'extra' from the multiplier)
                    // Let's go with summing the 'extra' influence of each relevant bacteria.
                    // Or, simpler for now: find the max, then apply it with other multipliers.
                    max_bacteria_specific_multiplier = max_bacteria_specific_multiplier.max(specific_multiplier);
                }
            }
        }
        // Apply the highest relevant bacteria-specific multiplier
        administration_prob *= max_bacteria_specific_multiplier;


        // Apply general multipliers AFTER the base rate and bacteria-specific multiplier
        if has_any_infection { administration_prob *= drug_infection_present_multiplier; }
        if has_any_identified_infection { administration_prob *= drug_test_identified_multiplier; }
        if initial_on_any_antibiotic || drugs_initiated_this_time_step > 0 {
            administration_prob *= already_on_drug_initiation_multiplier;
        }
        administration_prob *= syndrome_administration_multiplier;

        administration_prob = administration_prob.clamp(0.0, 1.0); // Clamp final probability

        if drugs_initiated_this_time_step < 2 && !individual.cur_use_drug[drug_idx] {
            if rng.gen_bool(administration_prob) {
                individual.cur_use_drug[drug_idx] = true;
                individual.date_drug_initiated[drug_idx] = time_step as i32;

                // Debug print
                println!(
                    "DEBUG: Individual {} started {} on day {}. Initiated due to (admin_prob: {:.4}). Current level: {:.2}",
                    individual.id,
                    drug_name,
                    individual.date_drug_initiated[drug_idx],
                    administration_prob,
                    get_drug_param(drug_name, "initial_level").unwrap_or(10.0)
                );

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



    // Drug-specific toxicity
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

    // --- DEATH LOGIC START ---
    if individual.date_of_death.is_none() {
        let mut prob_of_death_today = 0.0;
        let mut cause: Option<String> = None;
        let base_background_rate = get_global_param("base_background_mortality_rate_per_day")
            .expect("Missing base_background_mortality_rate_per_day in config");
        let age_multiplier = get_global_param("age_mortality_multiplier_per_year")
            .expect("Missing age_mortality_multiplier_per_year in config");
        let mut background_risk = base_background_rate;
        background_risk += (individual.age as f64 / 365.0) * age_multiplier;
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
            let sepsis_absolute_death_risk = get_global_param("sepsis_absolute_death_risk_per_day")
                .expect("Missing sepsis_absolute_death_risk_per_day in config");
            prob_not_dying *= 1.0 - sepsis_absolute_death_risk;
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
        prob_of_death_today = 1.0 - prob_not_dying;
        prob_of_death_today = prob_of_death_today.clamp(0.0, 1.0);
        if rng.gen::<f64>() < prob_of_death_today {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = cause.or(Some("background_mortality".to_string()));
        }
    }
    // --- DEATH LOGIC END ---

    // --- Update per-bacteria fields ---
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        let is_infected = individual.level[b_idx] > 0.001;

        if !is_infected {
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
            if individual.vaccination_status[b_idx] {
                let vaccine_efficacy = get_bacteria_param(bacteria, "vaccine_efficacy").unwrap_or(0.0);
                acquisition_probability *= 1.0 - vaccine_efficacy;
            }

            // Microbiome presence effect
            if individual.presence_microbiome[b_idx] {
                let microbiome_infection_multiplier = get_bacteria_param(bacteria, "microbiome_infection_acquisition_multiplier")
                    .unwrap_or_else(|| get_global_param("default_microbiome_infection_acquisition_multiplier").expect("Missing default_microbiome_infection_acquisition_multiplier in config"));
                acquisition_probability *= microbiome_infection_multiplier;
            }

            // --- Microbiome Presence (Carriage) ---
            if !individual.presence_microbiome[b_idx] {
                let microbiome_acquisition_multiplier = get_bacteria_param(bacteria, "microbiome_acquisition_multiplier")
                    .unwrap_or_else(|| get_global_param("default_microbiome_acquisition_multiplier").expect("Missing default_microbiome_acquisition_multiplier in config"));
                let microbiome_acquisition_probability = acquisition_probability * microbiome_acquisition_multiplier;
                if rng.gen_bool(microbiome_acquisition_probability.clamp(0.0, 1.0)) {
                    individual.presence_microbiome[b_idx] = true;
                }
            } else {
                let microbiome_clearance_prob = get_bacteria_param(bacteria, "microbiome_clearance_probability_per_day")
                    .unwrap_or_else(|| get_global_param("default_microbiome_clearance_probability_per_day").expect("Missing default_microbiome_clearance_probability_per_day in config"));
                if rng.gen_bool(microbiome_clearance_prob) {
                    individual.presence_microbiome[b_idx] = false;
                }
            }

            // ...resistance transfer logic...
            for &drug in DRUG_SHORT_NAMES.iter() {
                let d_idx = *drug_indices.get(drug).unwrap();
                if !individual.presence_microbiome[b_idx] {
                    individual.resistances[b_idx][d_idx].microbiome_r = 0.0;
                } else {
                    let infection_present = individual.level[b_idx] > 0.0;
                    if infection_present {
                        let current_any_r = individual.resistances[b_idx][d_idx].any_r;
                        let current_microbiome_r = individual.resistances[b_idx][d_idx].microbiome_r;
                        let should_attempt_transfer = (current_any_r > 0.0 && current_microbiome_r == 0.0) ||
                                                     (current_microbiome_r > 0.0 && current_any_r == 0.0);
                        if should_attempt_transfer && rng.gen_bool(transfer_prob) {
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

                // FIX: Use 'bacteria' directly, as 'bacteria_name' is not defined here.
    //          println!("                     acquisition_probability {}:", bacteria); // Changed bacteria_name to bacteria
  

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
                individual.infectious_syndrome[b_idx] = syndrome_id;


                let env_acquisition_chance = get_bacteria_param(bacteria, "environmental_acquisition_proportion").unwrap_or(0.1);
                individual.cur_infection_from_environment[b_idx] = rng.gen::<f64>() < env_acquisition_chance;

                let hospital_acquired_chance = get_bacteria_param(bacteria, "hospital_acquired_proportion").unwrap_or(0.05);
                let mut is_hospital_acquired = false;

                // Only consider hospital-acquired if the individual is currently hospitalized
                if individual.hospital_status.is_hospitalized() {
                    is_hospital_acquired = rng.gen::<f64>() < hospital_acquired_chance;
                }
                individual.infection_hospital_acquired[b_idx] = is_hospital_acquired;


                // --- any_r AND majority_r SETTING LOGIC ON NEW INFECTION ACQUISITION ---
                // in a newly infected person we should sample majority_r / any_r from all people in the same region with that 
                // bacteria and assign the newly infected person that level 
                let env_majority_r_level = get_global_param("environmental_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let hospital_majority_r_level = get_global_param("hospital_majority_r_level_for_new_acquisition").unwrap_or(0.0);
                let max_resistance_level = get_global_param("max_resistance_level").unwrap_or(1.0);

                // Use this before the for loop:
                let is_from_environment = individual.cur_infection_from_environment[b_idx];
                let is_hospital_acquired = individual.infection_hospital_acquired[b_idx];

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


                // Get the initial immunity level from config, or default to a reasonable value
                // todo: do we want to over-write the existing level of immunity at start of infection ?  I dont think so 
                // commenting out this code below for this reason but may need to review
      //        let initial_immunity = get_bacteria_param(bacteria, "initial_immunity_on_infection").unwrap_or(1.0); 
      //        individual.immune_resp[b_idx] = initial_immunity.max(0.0001); // Ensure it starts above the floor


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
                    let current_bacteria_level = individual.level[b_idx];

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


                    // Calculate activity_r (should always be updated)
                    // todo: may need to specify the parameter 0.05 below in config.rs
                    if drug_current_level > 0.0 {
                        let normalized_any_r = resistance_data.any_r / max_resistance_level;
                        resistance_data.activity_r = 0.05 * drug_current_level * (1.0 - normalized_any_r);
                    } else {
                        resistance_data.activity_r = 0.0;
                    }
                }
            }
        }

        // Testing and Diagnosis logic
        let last_infected_time = individual.date_last_infected[b_idx];
        let current_test_status_entry = &mut individual.test_identified_infection[b_idx];
        let test_delay_days = get_global_param("test_delay_days").unwrap_or(3.0) as i32;
        let test_rate_per_day = get_global_param("test_rate_per_day").unwrap_or(0.15);
        if !*current_test_status_entry && (time_step as i32) >= (last_infected_time + test_delay_days) {
            if rng.gen_bool(test_rate_per_day.clamp(0.0, 1.0)) {
                *current_test_status_entry = true;
            }
        }

       // Bacteria level change (growth/decay)
        // This entire block should only execute if the individual is currently infected with this bacteria
        if is_infected { // <--- ADD THIS LINE
            let immunity_level = individual.immune_resp[b_idx];
            let baseline_change = get_bacteria_param(bacteria, "level_change_rate_baseline").unwrap_or(0.0);
            let reduction_due_to_immune_resp = get_bacteria_param(bacteria, "immunity_effect_on_level_change").unwrap_or(0.0);
            let mut total_reduction_due_to_antibiotic = 0.0;


                        // --- DEBUGGING START ---
            if individual.id == 0 {
                println!("--- DEBUG (Day {}), Individual {}: Bacteria {} ---", time_step, individual.id, bacteria);
                println!("  Immunity Level: {:.4}", immunity_level);
                println!("  Baseline Change: {:.4}", baseline_change);
                println!("  Reduction due to Immune Response: {:.4}", immunity_level * reduction_due_to_immune_resp);
            }
            // --- DEBUGGING END ---



            for (drug_idx, _drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                if individual.cur_level_drug[drug_idx] > 0.0 {
                    let resistance_data = &individual.resistances[b_idx][drug_idx];
                    total_reduction_due_to_antibiotic += resistance_data.activity_r;



                                            // --- DEBUGGING START ---
                        println!(
                            "  Drug {}: Current Level = {:.4}, Activity_r = {:.4}",
                            DRUG_SHORT_NAMES[drug_idx],
                            individual.cur_level_drug[drug_idx],
                            resistance_data.activity_r
                        );
                        // --- DEBUGGING END ---    



                        
                        
             // --- DEBUGGING START ---
                println!("  Total Reduction Due To Antibiotic: {:.4}", total_reduction_due_to_antibiotic);
                // --- DEBUGGING END ---



                }
            }
            let decay_rate = baseline_change - (immunity_level * reduction_due_to_immune_resp) - total_reduction_due_to_antibiotic;

            let max_level = get_bacteria_param(bacteria, "max_level").unwrap_or(100.0);
            let new_level = (individual.level[b_idx] + decay_rate).max(0.0).min(max_level);



                // --- DEBUGGING START ---
                println!("  Calculated Decay Rate: {:.4}", decay_rate);
                println!("  Old Bacteria Level: {:.4}", individual.level[b_idx]);
                println!("  New Bacteria Level: {:.4}", new_level);
                println!("--------------------------------------------------");
                // --- DEBUGGING END ---





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

        // Immunity increase   todo: need a max for immune response ?
        let infection_start_time = individual.date_last_infected[b_idx];
        let time_since_infection = (time_step as i32) - infection_start_time;
        let age = individual.age;
        let mut immune_increase = get_bacteria_param(bacteria, "immunity_increase_rate_baseline").unwrap_or(0.0);
        immune_increase += time_since_infection as f64 * get_bacteria_param(bacteria, "immunity_increase_rate_per_day").unwrap_or(0.0);
        immune_increase += individual.level[b_idx] * get_bacteria_param(bacteria, "immunity_increase_rate_per_level").unwrap_or(0.0);
        let age_modifier = get_bacteria_param(bacteria, "immunity_age_modifier").unwrap_or(1.0);
        immune_increase *= age_modifier.powf((-age as f64 / 365.0) / 50.0);
        individual.immune_resp[b_idx] = (individual.immune_resp[b_idx] + immune_increase).max(0.0001);
    }


    // todo: decay in immunity when bacteria no longer present - codde below ok but need to condition on infection not present

    // Immunity decay (applies regardless of current infection status)
   
//  for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
//      let current_immunity = individual.immune_resp[b_idx];
//      let baseline_immunity = get_bacteria_param(bacteria, "baseline_immunity_level").unwrap_or(0.0001);
//      let decay_rate = get_bacteria_param(bacteria, "immunity_decay_rate").unwrap_or(0.0);
//      if current_immunity > baseline_immunity {
//          individual.immune_resp[b_idx] = (current_immunity - decay_rate).max(baseline_immunity);
//      }
//  }



    // Set test_r to 0 for all bacteria/drug combos
    for i in 0..BACTERIA_LIST.len() {
        for j in 0..DRUG_SHORT_NAMES.len() {
            individual.resistances[i][j].test_r = 0.0;
        }
    }
}