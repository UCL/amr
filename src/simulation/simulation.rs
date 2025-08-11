// src/simulation/simulation.rs
// Main simulation logic and summary data structures for AMR model.
//
// Contains:
//   - TimeStepSummary: struct for per-timestep summary statistics
//   - Simulation: struct and methods for running the simulation, managing population, and logging
//   - Initialization of lookup tables for bacteria, drugs, and cross-resistance
//   - Debug/print blocks for individual and population state
//

// search below for "printing of variable values for individual 0"
// when want to print variable values for individual 0 for de-bugging


use crate::simulation::population::{Population, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::rules::apply_rules;
use crate::config; // Import the config module
use std::collections::HashMap;
use rayon::prelude::*;
// Removed most atomics by using thread-local aggregation; retain no atomic imports here.
use std::time::Instant;
use std::io::Write;

// Compact structure for time step summary data
#[allow(dead_code)]
#[derive(Clone)]
// Summary statistics for each simulation time step.
//
// Captures population-level and per-bacteria/drug summary metrics for each time step.
pub struct TimeStepSummary {
    pub time_step: usize,
    pub total_population: usize,
    pub total_deaths: usize,
    pub deaths_background: usize,        // Deaths from background mortality
    pub deaths_sepsis: usize,           // Deaths from sepsis
    pub deaths_drug_toxicity: usize,    // Deaths from drug toxicity
    pub deaths_past_year: usize, // all-cause     // Rolling 1-year (365 days) death counts
    pub deaths_background_past_year: usize,     // Rolling 1-year (365 days) death counts
    pub deaths_sepsis_past_year: usize,     // Rolling 1-year (365 days) death counts
    pub deaths_drug_toxicity_past_year: usize,     // Rolling 1-year (365 days) death counts
    pub total_with_resistance: usize,
    pub total_currently_infected: usize, // Number of living people currently infected with any bacteria
    pub currently_taking_drug_count: usize, 
    pub infected_10_days_count: usize,     
    pub infected_30_days_count: usize,     
    pub taking_two_drugs_count: usize,     
    pub number_in_hospital: usize,         
    pub number_severely_immunosuppressed: usize, 
    pub number_with_sepsis: usize,         
    pub infections_by_bacteria: Vec<usize>, // indexed by bacteria
    pub resistance_by_bacteria_drug: Vec<Vec<usize>>, // [bacteria][drug] counts
    pub newly_infected_count: usize, // Number of people newly infected this time step
    pub newly_infected_with_resistance_count: usize, // NEW: Number of newly infected people who acquired resistance
    pub newly_infected_past_year: usize, // Rolling 1-year (365 days) newly infected count
    pub currently_infected_and_on_drug_count: usize, // NEW: intersection of currently infected AND on any drug
    pub num_age_0_5: usize,
    pub num_age_6_14: usize,
    pub num_age_15_49: usize,
    pub num_age_50_79: usize,
    pub num_age_80plus: usize,
    pub num_with_any_bacteria_microbiome: usize, // NEW: number of people with any presence_microbiome=true

    // per-bacteria, per-drug infection and resistance counts (flat, len = bacteria * drugs)
    pub infected_and_standardized_mic_lt2_by_bacteria_drug: Vec<usize>,

    // NEW: per-drug currently on drug counts (indexed by drug)
    pub currently_on_drug_by_drug: Vec<usize>,
} 

// Main simulation struct: holds population, time steps, and lookup tables.
//
// Encapsulates the state and configuration of a simulation run, including population, time steps,
// and lookup tables for bacteria, drugs, and cross-resistance groups.
pub struct Simulation {
    pub population: Population,
    pub time_steps: usize,
    pub log_individuals: bool,
    /// Maps bacteria names to their indices in arrays.
    pub bacteria_indices: HashMap<&'static str, usize>,
    /// Maps drug names to their indices in arrays.
    pub drug_indices: HashMap<&'static str, usize>,
    /// Maps bacteria index to cross-resistance groups (each group is a Vec of drug indices).
    pub cross_resistance_groups: HashMap<usize, Vec<Vec<usize>>>,
    /// Stores majority_r positive values for each (bacteria, is_microbiome, drug, time_step) combo.
    pub current_majority_r_positive_values_by_combo: HashMap<(usize, bool, usize, usize), Vec<f64>>,
    /// Efficient storage for summary data at each time step.
    pub summary_log: Vec<TimeStepSummary>,
    /// Pre-computed parameter keys to avoid string allocation during simulation.
    pub param_cache: crate::rules::ParameterKeyCache,
    /// Precomputed potency values indexed by [bacteria * num_drugs + drug]
    pub potency_matrix: Vec<f64>,
    /// Precomputed majority_r threshold below which standardized MIC < 2 (avoids per-step division)
    pub mic_lt2_majority_r_thresholds: Vec<f64>,
    /// Hint: previous timestep total majority_r entries to reserve capacity
    pub prev_majority_r_entries_len: usize,
}

impl Simulation {
    /// Create a new Simulation instance with initialized population and lookup tables.
    ///
    /// Initializes population, bacteria/drug indices, and cross-resistance groups.
    pub fn new(population_size: usize, time_steps: usize, log_individuals: bool) -> Self {
        let population = Population::new(population_size);
        // public function named new (rust’s conventional constructor pattern).  
        // takes two inputs: population_size: how many individuals to initialize.
        // time_steps: how many time steps the simulation should run.
        // returns Self → shorthand for returning an instance of Simulation.
        // calls a new constructor for the Population struct.  Passes in "population_size", returning a Population instance 
        // and stores it in the local population variable.

        // Initialize bacteria_indices and drug_indices
        let mut bacteria_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &bacteria) in BACTERIA_LIST.iter().enumerate() {
            bacteria_indices.insert(bacteria, i);
        }
        let mut drug_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            drug_indices.insert(drug, i);
        }
        // Load and process cross-resistance groups
        let mut cross_resistance_groups = HashMap::new();
        let raw_groups = config::get_cross_resistance_groups();
        for (bacteria_name, groups) in raw_groups.iter() {
            if let Some(&b_idx) = bacteria_indices.get(bacteria_name) {
                let indexed_groups: Vec<Vec<usize>> = groups.iter().map(|group| {
                    group.iter().filter_map(|drug_name| drug_indices.get(drug_name).copied()).collect()
                }).collect();
                cross_resistance_groups.insert(b_idx, indexed_groups);
            }
        }

        println!(" ");
        println!("--- simulation.rs  initial state of individual 0 ---");
        println!(" ");
        println!("id: {}", population.individuals[0].id);
        println!("age: {} days", population.individuals[0].age);
        println!("sex at birth: {}", population.individuals[0].sex_at_birth);
        println!("region living: {:?}", population.individuals[0].region_living);
        println!("region currently in: {:?}", population.individuals[0].region_cur_in);
        println!("current_infection_related_death_risk: {:.2}", population.individuals[0].current_infection_related_death_risk);
        println!("background_all_cause_mortality_rate: {:.4}", population.individuals[0].background_all_cause_mortality_rate);
        println!("sexual_contact_level: {:.2}", population.individuals[0].sexual_contact_level);
        println!("airborne_contact_level_with_adults: {:.2}", population.individuals[0].airborne_contact_level_with_adults);
        println!("airborne_contact_level_with_children: {:.2}", population.individuals[0].airborne_contact_level_with_children);
        println!("oral_exposure_level: {:.2}", population.individuals[0].oral_exposure_level);
        println!("mosquito_exposure_level: {:.2}", population.individuals[0].mosquito_exposure_level);
        println!("current_toxicity: {:.2}", population.individuals[0].current_toxicity);
        println!("mortality_risk_current_toxicity: {:.2}", population.individuals[0].mortality_risk_current_toxicity);
        println!(" ");

        // Precompute potency matrix to avoid repeated string formatting/hash lookups in hot loop
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();
    let mut potency_matrix = Vec::with_capacity(num_bacteria * num_drugs);
    let mut mic_lt2_majority_r_thresholds = Vec::with_capacity(num_bacteria * num_drugs);
        for b_idx in 0..num_bacteria {
            for d_idx in 0..num_drugs {
                let bacteria_name = BACTERIA_LIST[b_idx];
                let drug_name = DRUG_SHORT_NAMES[d_idx];
                let key = format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug_name, bacteria_name);
                let potency = crate::config::PARAMETERS.get(&key).copied().unwrap_or(0.01);
                potency_matrix.push(potency);
        // standardized_mic = 1 / ((1 - r)*potency) < 2  =>  r < 1 - 0.5 / potency
        // Precompute threshold to avoid division in hot loop; if potency very small threshold will be negative
        let threshold = 1.0 - 0.5 / potency;
        mic_lt2_majority_r_thresholds.push(threshold);
            }
        }

        Simulation { // Constructs and returns a new Simulation instance with the initialized population, time steps, and other data structures.
            population,
            time_steps,
            log_individuals,
            bacteria_indices,
            drug_indices,
            cross_resistance_groups, 
            current_majority_r_positive_values_by_combo: HashMap::new(), // Initialize empty
            summary_log: Vec::new(), // Initialize empty log
            param_cache: crate::rules::ParameterKeyCache::new(),
            potency_matrix,
            mic_lt2_majority_r_thresholds,
            prev_majority_r_entries_len: 0,
        }
    }

    pub fn run(&mut self) {
        // public function named run, which executes the simulation for the specified number of time steps.

        println!(" ");
        println!("--- starting to run over time steps");
        println!(" ");

        for t in 0..self.time_steps {
            let timestep_start = Instant::now();
            
            // --- Setup counters; MIC<2 snapshot will use per-thread local vectors reduced after loop (avoids atomic contention) ---
            let num_bacteria = BACTERIA_LIST.len();
            let num_drugs = DRUG_SHORT_NAMES.len();
            
//             let calculation_time = calculation_start.elapsed();
//             if t % 100 == 0 { // Log every 10th timestep
//                 println!("Time step {}", t);
//             }
//          println!("simulation.rs time step: {}", t);

            // Thread-local aggregation will replace most atomics; keep only minimal atomics if needed (none for now).

            // Use previous time step's resistance data for new acquisitions
            let previous_majority_r_positive_values_by_combo = if t == 0 {
                HashMap::new() // Empty for first time step
            } else {
            // Use the data collected in the previous iteration
                std::mem::take(&mut self.current_majority_r_positive_values_by_combo)
            };

            // LocalTotals structure for thread-local aggregation
            struct LocalTotals {
                mic_lt2_counts: Vec<usize>,
                infections_by_bacteria: Vec<usize>,
                resistance_by_bacteria_drug: Vec<usize>,
                currently_on_drug_by_drug: Vec<usize>,
                majority_r_entries: Vec<((usize, bool, usize, usize), f64)>,
                total_deaths: usize,
                deaths_background: usize,
                deaths_sepsis: usize,
                deaths_drug_toxicity: usize,
                currently_taking_drug_count: usize,
                infected_10_days_count: usize,
                infected_30_days_count: usize,
                taking_two_drugs_count: usize,
                number_in_hospital: usize,
                number_severely_immunosuppressed: usize,
                number_with_sepsis: usize,
                newly_infected_count: usize,
                newly_infected_with_resistance_count: usize,
                total_currently_infected: usize,
                total_with_resistance: usize,
                currently_infected_and_on_drug_count: usize,
                num_with_any_bacteria_microbiome: usize,
                // Integrated previously sequential counts:
                living_population: usize,
                num_age_0_5: usize,
                num_age_6_14: usize,
                num_age_15_49: usize,
                num_age_50_79: usize,
                num_age_80plus: usize,
            }
            impl LocalTotals {
                fn new(num_bacteria: usize, num_drugs: usize, majority_r_capacity: usize) -> Self {
                    Self {
                        mic_lt2_counts: vec![0; num_bacteria * num_drugs],
                        infections_by_bacteria: vec![0; num_bacteria],
                        resistance_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        currently_on_drug_by_drug: vec![0; num_drugs],
                        majority_r_entries: Vec::with_capacity(majority_r_capacity),
                        total_deaths: 0,
                        deaths_background: 0,
                        deaths_sepsis: 0,
                        deaths_drug_toxicity: 0,
                        currently_taking_drug_count: 0,
                        infected_10_days_count: 0,
                        infected_30_days_count: 0,
                        taking_two_drugs_count: 0,
                        number_in_hospital: 0,
                        number_severely_immunosuppressed: 0,
                        number_with_sepsis: 0,
                        newly_infected_count: 0,
                        newly_infected_with_resistance_count: 0,
                        total_currently_infected: 0,
                        total_with_resistance: 0,
                        currently_infected_and_on_drug_count: 0,
                        num_with_any_bacteria_microbiome: 0,
                        living_population: 0,
                        num_age_0_5: 0,
                        num_age_6_14: 0,
                        num_age_15_49: 0,
                        num_age_50_79: 0,
                        num_age_80plus: 0,
                    }
                }
                fn merge(&mut self, other: Self) {
                    for (a,b) in self.mic_lt2_counts.iter_mut().zip(other.mic_lt2_counts) { *a += b; }
                    for (a,b) in self.infections_by_bacteria.iter_mut().zip(other.infections_by_bacteria) { *a += b; }
                    for (a,b) in self.resistance_by_bacteria_drug.iter_mut().zip(other.resistance_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.currently_on_drug_by_drug.iter_mut().zip(other.currently_on_drug_by_drug) { *a += b; }
                    self.majority_r_entries.extend(other.majority_r_entries);
                    self.total_deaths += other.total_deaths;
                    self.deaths_background += other.deaths_background;
                    self.deaths_sepsis += other.deaths_sepsis;
                    self.deaths_drug_toxicity += other.deaths_drug_toxicity;
                    self.currently_taking_drug_count += other.currently_taking_drug_count;
                    self.infected_10_days_count += other.infected_10_days_count;
                    self.infected_30_days_count += other.infected_30_days_count;
                    self.taking_two_drugs_count += other.taking_two_drugs_count;
                    self.number_in_hospital += other.number_in_hospital;
                    self.number_severely_immunosuppressed += other.number_severely_immunosuppressed;
                    self.number_with_sepsis += other.number_with_sepsis;
                    self.newly_infected_count += other.newly_infected_count;
                    self.newly_infected_with_resistance_count += other.newly_infected_with_resistance_count;
                    self.total_currently_infected += other.total_currently_infected;
                    self.total_with_resistance += other.total_with_resistance;
                    self.currently_infected_and_on_drug_count += other.currently_infected_and_on_drug_count;
                    self.num_with_any_bacteria_microbiome += other.num_with_any_bacteria_microbiome;
                    self.living_population += other.living_population;
                    self.num_age_0_5 += other.num_age_0_5;
                    self.num_age_6_14 += other.num_age_6_14;
                    self.num_age_15_49 += other.num_age_15_49;
                    self.num_age_50_79 += other.num_age_50_79;
                    self.num_age_80plus += other.num_age_80plus;
                }
            }

            // Single pass: apply rules and collect all statistics
            let _rules_start = Instant::now();
            
        let mic_lt2_thresholds = &self.mic_lt2_majority_r_thresholds;
        let threads = rayon::current_num_threads().max(1);
        let per_thread_cap = (self.prev_majority_r_entries_len / threads).saturating_add(8);
        let totals = self.population.individuals.par_iter_mut()
            .fold(|| LocalTotals::new(num_bacteria, num_drugs, per_thread_cap), |mut lt, individual| {
                        // Pre-rules MIC snapshot
                        if individual.date_of_death.is_none() {
                            for b_idx in 0..num_bacteria {
                                if individual.level[b_idx] > 0.001 {
                                    let base = b_idx * num_drugs;
                                    for d_idx in 0..num_drugs {
                                        let resistance_data = &individual.resistances[b_idx][d_idx];
                                        let threshold = mic_lt2_thresholds[base + d_idx];
                                        if resistance_data.majority_r < threshold { lt.mic_lt2_counts[base + d_idx] += 1; }
                                    }
                                }
                            }
                        }

                        // Apply rules
                        apply_rules(
                            individual,
                            t,
                            &previous_majority_r_positive_values_by_combo,
                            &self.bacteria_indices,
                            &self.drug_indices,
                            &self.cross_resistance_groups,
                            &self.param_cache,
                        );

                        // Death accounting
                        if let Some(death_time) = individual.date_of_death {
                            if death_time == t {
                                lt.total_deaths += 1;
                                if let Some(ref cause) = individual.cause_of_death {
                                    match cause.as_str() {
                                        "background_mortality" => lt.deaths_background += 1,
                                        "sepsis_related" => lt.deaths_sepsis += 1,
                                        "drug_toxicity_related" => lt.deaths_drug_toxicity += 1,
                                        _ => lt.deaths_background += 1,
                                    }
                                } else { lt.deaths_background += 1; }
                            }
                        }

                        if individual.date_of_death.is_none() {
                            // Integrated living population & age groups
                            lt.living_population += 1;
                            let age_years = individual.age as f64 / 365.0;
                            if (0.0..6.0).contains(&age_years) { lt.num_age_0_5 += 1; }
                            else if (6.0..15.0).contains(&age_years) { lt.num_age_6_14 += 1; }
                            else if (15.0..50.0).contains(&age_years) { lt.num_age_15_49 += 1; }
                            else if (50.0..80.0).contains(&age_years) { lt.num_age_50_79 += 1; }
                            else if age_years >= 80.0 { lt.num_age_80plus += 1; }
                            // Drug usage post-rules
                            let mut on_any_drug = false;
                            for (d_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
                                if is_using { lt.currently_on_drug_by_drug[d_idx] += 1; on_any_drug = true; }
                            }
                            if on_any_drug { lt.currently_taking_drug_count += 1; }

                            if individual.presence_microbiome.iter().any(|&x| x) { lt.num_with_any_bacteria_microbiome += 1; }

                            // Infection & resistance
                            let mut individual_max_infection_duration = 0;
                            let mut individual_has_any_infection = false;
                            let mut individual_has_any_r_positive = false;
                            let mut was_newly_infected = false;
                            let mut was_newly_infected_with_resistance = false;
                            let is_currently_infected_any;
                            {
                                let mut infected_any_tmp = false;
                                for b_idx in 0..num_bacteria {
                                    if individual.level[b_idx] > 0.001 {
                                        infected_any_tmp = true;
                                        individual_has_any_infection = true;
                                        lt.infections_by_bacteria[b_idx] += 1;
                                        let days_since_infection = t as i32 - individual.date_last_infected[b_idx];
                                        if days_since_infection > individual_max_infection_duration { individual_max_infection_duration = days_since_infection; }
                                        if individual.date_last_infected[b_idx] == t as i32 { was_newly_infected = true; }
                                        let base = b_idx * num_drugs;
                                        for d_idx in 0..num_drugs {
                                            let resistance_data = &individual.resistances[b_idx][d_idx];
                                            if resistance_data.majority_r > 0.0 { lt.resistance_by_bacteria_drug[base + d_idx] += 1; lt.majority_r_entries.push(((individual.region_cur_in as usize, individual.hospital_status.is_hospitalized(), b_idx, d_idx), resistance_data.majority_r)); }
                                            if resistance_data.any_r > 0.0 {
                                                individual_has_any_r_positive = true;
                                                if individual.date_last_infected[b_idx] == t as i32 && !was_newly_infected_with_resistance {
                                                    lt.newly_infected_with_resistance_count += 1;
                                                    was_newly_infected_with_resistance = true;
                                                }
                                            }
                                        }
                                    }
                                }
                                is_currently_infected_any = infected_any_tmp;
                            }
                            if is_currently_infected_any && on_any_drug { lt.currently_infected_and_on_drug_count += 1; }
                            if individual_has_any_infection { lt.total_currently_infected += 1; }
                            if individual_has_any_r_positive { lt.total_with_resistance += 1; }
                            if individual_max_infection_duration > 10 { lt.infected_10_days_count += 1; }
                            if individual_max_infection_duration > 30 { lt.infected_30_days_count += 1; }
                            if was_newly_infected { lt.newly_infected_count += 1; }
                            let active_drug_count = individual.cur_use_drug.iter().filter(|&&x| x).count();
                            if active_drug_count >= 2 { lt.taking_two_drugs_count += 1; }
                            if individual.hospital_status.is_hospitalized() { lt.number_in_hospital += 1; }
                            if individual.is_severely_immunosuppressed { lt.number_severely_immunosuppressed += 1; }
                            if individual.sepsis.iter().any(|&s| s) { lt.number_with_sepsis += 1; }
                        }
                        lt
                    })
                    .reduce(|| LocalTotals::new(num_bacteria, num_drugs, per_thread_cap), |mut a, b| { a.merge(b); a });

                // Destructure to move out (avoid cloning large vectors)
                let LocalTotals {
                    mic_lt2_counts: infected_and_standardized_mic_lt2_by_bacteria_drug,
                    infections_by_bacteria: infections_by_bacteria_vec,
                    resistance_by_bacteria_drug: resistance_by_bacteria_drug_flat,
                    currently_on_drug_by_drug,
                    majority_r_entries,
                    total_deaths,
                    deaths_background,
                    deaths_sepsis,
                    deaths_drug_toxicity,
                    currently_taking_drug_count,
                    infected_10_days_count,
                    infected_30_days_count,
                    taking_two_drugs_count,
                    number_in_hospital,
                    number_severely_immunosuppressed,
                    number_with_sepsis,
                    newly_infected_count,
                    newly_infected_with_resistance_count,
                    total_currently_infected,
                    total_with_resistance,
                    currently_infected_and_on_drug_count,
                    num_with_any_bacteria_microbiome,
                    living_population,
                    num_age_0_5,
                    num_age_6_14,
                    num_age_15_49,
                    num_age_50_79,
                    num_age_80plus,
                } = totals;

                // Rebuild 2D resistance structure for summary
                let mut resistance_by_bacteria_drug: Vec<Vec<usize>> = Vec::with_capacity(num_bacteria);
                for b_idx in 0..num_bacteria { resistance_by_bacteria_drug.push(resistance_by_bacteria_drug_flat[b_idx*num_drugs..(b_idx+1)*num_drugs].to_vec()); }

                // Pre-size HashMap: heuristic half number of entries (since entries will group by key)
                let mut new_majority_r_positive_values_by_combo: HashMap<(usize, bool, usize, usize), Vec<f64>> = HashMap::with_capacity(majority_r_entries.len() / 2 + 1);
                for (key, val) in majority_r_entries { new_majority_r_positive_values_by_combo.entry(key).or_insert_with(Vec::new).push(val); }
                // Update capacity hint for next timestep (sum of vector lengths)
                let total_entries: usize = new_majority_r_positive_values_by_combo.values().map(|v| v.len()).sum();
                self.prev_majority_r_entries_len = total_entries;
            
            // let rules_time = rules_start.elapsed();
            // if t % 10 == 0 { // Log every 10th timestep
            //     println!("Time step {}: rules application took {:.3}ms", t, rules_time.as_secs_f64() * 1000.0);
            // }

            // Collect remaining statistics that need sequential access
            // No need for sequential pass for per-bacteria/drug majority_r counts

            // Store for next iteration
            self.current_majority_r_positive_values_by_combo = new_majority_r_positive_values_by_combo;

            // Create summary for this time step
            let infected_10_count = infected_10_days_count;
            let infected_30_count = infected_30_days_count;

            // Optional debug (uncomment if needed)
            // if t % 500 == 0 { println!("Time step {} drug usage counts: {:?}", t, currently_on_drug_by_drug); }

            let summary = TimeStepSummary {
                infected_and_standardized_mic_lt2_by_bacteria_drug,
                currently_on_drug_by_drug,
                num_age_0_5,
                num_age_6_14,
                num_age_15_49,
                num_age_50_79,
                num_age_80plus,
                num_with_any_bacteria_microbiome,
                time_step: t,
                total_population: living_population,
                number_in_hospital,
                number_severely_immunosuppressed,
                number_with_sepsis,
                newly_infected_count,
                newly_infected_with_resistance_count,
                total_currently_infected,
                total_with_resistance,
                infected_10_days_count: infected_10_count,
                infected_30_days_count: infected_30_count,
                currently_taking_drug_count,
                taking_two_drugs_count,
                infections_by_bacteria: infections_by_bacteria_vec,
                resistance_by_bacteria_drug,
                total_deaths,
                deaths_background,
                deaths_sepsis,
                deaths_drug_toxicity,
                // Rolling 1-year (365 days) death counts
                deaths_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.total_deaths)
                        .sum::<usize>()
            + total_deaths
                        - self.summary_log.last().map_or(0, |s| s.total_deaths)
                },
                deaths_background_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.deaths_background)
                        .sum::<usize>()
            + deaths_background
                        - self.summary_log.last().map_or(0, |s| s.deaths_background)
                },
                deaths_sepsis_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.deaths_sepsis)
                        .sum::<usize>()
            + deaths_sepsis
                        - self.summary_log.last().map_or(0, |s| s.deaths_sepsis)
                },
                deaths_drug_toxicity_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.deaths_drug_toxicity)
                        .sum::<usize>()
            + deaths_drug_toxicity
                        - self.summary_log.last().map_or(0, |s| s.deaths_drug_toxicity)
                },
                newly_infected_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.newly_infected_count)
                        .sum::<usize>()
            + newly_infected_count
                        - self.summary_log.last().map_or(0, |s| s.newly_infected_count)
                },
        currently_infected_and_on_drug_count: currently_infected_and_on_drug_count,
            };



            // Comprehensive print block for individual 0
            let _individual_0 = &self.population.individuals[0];
            // println!("--- Individual 0 full state ---");
            // println!("id: {}", individual_0.id);
            // println!("age (days): {}", individual_0.age);
            // println!("sex_at_birth: {}", individual_0.sex_at_birth);
            // println!("region_living: {:?}", individual_0.region_living);
            // println!("region_cur_in: {:?}", individual_0.region_cur_in);
            // println!("current_infection_related_death_risk: {:.4}", individual_0.current_infection_related_death_risk);
            // println!("background_all_cause_mortality_rate: {:.4}", individual_0.background_all_cause_mortality_rate);
            // println!("sexual_contact_level: {:.4}", individual_0.sexual_contact_level);
            // println!("airborne_contact_level_with_adults: {:.4}", individual_0.airborne_contact_level_with_adults);
            // println!("airborne_contact_level_with_children: {:.4}", individual_0.airborne_contact_level_with_children);
            // println!("oral_exposure_level: {:.4}", individual_0.oral_exposure_level);
            // println!("mosquito_exposure_level: {:.4}", individual_0.mosquito_exposure_level);
            // println!("current_toxicity: {:.4}", individual_0.current_toxicity);
            // println!("mortality_risk_current_toxicity: {:.4}", individual_0.mortality_risk_current_toxicity);
            // println!("hospital_status: {:?}", individual_0.hospital_status);
            // println!("is_severely_immunosuppressed: {:?}", individual_0.is_severely_immunosuppressed);
            // println!("date_of_death: {:?}", individual_0.date_of_death);
            // // Arrays
            // println!("level: {:?}", individual_0.level);
            // println!("immune_resp: {:?}", individual_0.immune_resp);
            // println!("presence_microbiome: {:?}", individual_0.presence_microbiome);
            // println!("cur_level_drug: {:?}", individual_0.cur_level_drug);
            // println!("cur_use_drug: {:?}", individual_0.cur_use_drug);
            // println!("ever_taken_drug: {:?}", individual_0.ever_taken_drug);
            // println!("date_last_infected: {:?}", individual_0.date_last_infected);
            // println!("cur_infection_from_environment: {:?}", individual_0.cur_infection_from_environment);
            // println!("infection_hospital_acquired: {:?}", individual_0.infection_hospital_acquired);
            // println!("test_identified_infection: {:?}", individual_0.test_identified_infection);
            // println!("sepsis: {:?}", individual_0.sepsis);
            // // Per-bacteria/drug resistance data
            // for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
            //     for (d_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
            //         let resistance = &individual_0.resistances[b_idx][d_idx];
            //         println!(
            //             "Resistance for bacteria {} and drug {}: any_r = {:.4}, activity_r = {:.4}, majority_r = {:.4}",
            //             bacteria_name, drug_name, resistance.any_r, resistance.activity_r, resistance.majority_r
            //         );
            //     }
            // }


            self.summary_log.push(summary);

            if self.log_individuals {
                use std::fs::OpenOptions;
                let n_log = 10.min(self.population.individuals.len());
                let log_path = "individuals_log.csv";
                let file_exists = std::path::Path::new(log_path).exists();
                let mut file = OpenOptions::new().create(true).append(true).open(log_path).expect("Unable to open individuals_log.csv");
                // Write header if file is new
                if !file_exists {
                    writeln!(file, "time_step,individual_index,id,age,sex_at_birth,region_living,region_cur_in,current_infection_related_death_risk,background_all_cause_mortality_rate,sexual_contact_level,airborne_contact_level_with_adults,airborne_contact_level_with_children,oral_exposure_level,mosquito_exposure_level,current_toxicity,mortality_risk_current_toxicity,hospital_status,is_severely_immunosuppressed,date_of_death,level,immune_resp,presence_microbiome,cur_level_drug,cur_use_drug,ever_taken_drug,date_last_infected,cur_infection_from_environment,infection_hospital_acquired,test_identified_infection,sepsis,resistances_microbiome_r,resistances_test_r,resistances_activity_r,resistances_any_r,resistances_majority_r,resistance_mechanisms").unwrap();
                }
                fn fmt_vec<T: std::fmt::Display>(v: &[T]) -> String {
                    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(";")
                }
                for i in 0..n_log {
                    let ind = &self.population.individuals[i];
                    // Flatten all resistance values for all bacteria/drugs
                    let mut microbiome_r = Vec::new();
                    let mut test_r = Vec::new();
                    let mut activity_r = Vec::new();
                    let mut any_r = Vec::new();
                    let mut majority_r = Vec::new();
                    for bact in &ind.resistances {
                        for res in bact {
                            microbiome_r.push(res.microbiome_r);
                            test_r.push(res.test_r);
                            activity_r.push(res.activity_r);
                            any_r.push(res.any_r);
                            majority_r.push(res.majority_r);
                        }
                    }
                    // Flatten all resistance mechanisms for all bacteria
                    let mut mechanisms = Vec::new();
                    for bact_mechs in &ind.resistance_mechanisms {
                        for &present in bact_mechs {
                            mechanisms.push(if present { "1" } else { "0" });
                        }
                    }
                    writeln!(file, "{},{},{},{},{},{:?},{:?},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:?},{},{:?},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                        t,
                        i,
                        ind.id,
                        ind.age,
                        ind.sex_at_birth,
                        ind.region_living,
                        ind.region_cur_in,
                        ind.current_infection_related_death_risk,
                        ind.background_all_cause_mortality_rate,
                        ind.sexual_contact_level,
                        ind.airborne_contact_level_with_adults,
                        ind.airborne_contact_level_with_children,
                        ind.oral_exposure_level,
                        ind.mosquito_exposure_level,
                        ind.current_toxicity,
                        ind.mortality_risk_current_toxicity,
                        format!("{:?}", ind.hospital_status),
                        ind.is_severely_immunosuppressed,
                        format!("{:?}", ind.date_of_death),
                        fmt_vec(&ind.level),
                        fmt_vec(&ind.immune_resp),
                        fmt_vec(&ind.presence_microbiome),
                        fmt_vec(&ind.cur_level_drug),
                        fmt_vec(&ind.cur_use_drug),
                        fmt_vec(&ind.ever_taken_drug),
                        fmt_vec(&ind.date_last_infected),
                        fmt_vec(&ind.cur_infection_from_environment),
                        fmt_vec(&ind.infection_hospital_acquired),
                        fmt_vec(&ind.test_identified_infection),
                        fmt_vec(&ind.sepsis),
                        fmt_vec(&microbiome_r),
                        fmt_vec(&test_r),
                        fmt_vec(&activity_r),
                        fmt_vec(&any_r),
                        fmt_vec(&majority_r),
                        mechanisms.join(";")
                    ).unwrap();
                }
            }
            
            let _timestep_time = timestep_start.elapsed();
            if t % 100 == 0 { // Log every 100th timestep
                println!("Time step {}", t);
                // println!("Time step {} total time: {:.3}ms", t, timestep_time.as_secs_f64() * 1000.0);
            }

        }

    }
 
    

    pub fn print_summary_statistics(&self) {
        if self.summary_log.is_empty() {
            println!("No summary data logged.");
            return;
        }

  

        // println!("\n--- Simulation Summary Statistics ---");
        // for summary in &self.summary_log {
        //     println!("Time step {}: {} newly infected, {} deaths, {} with resistance", 
        //         summary.time_step, 
        //         summary.newly_infected_count, 
        //         summary.total_deaths, 
        //         summary.total_with_resistance
        //     );
        // }



    }

    pub fn export_summary_to_csv(&self, filename: &str) -> Result<(), std::io::Error> {
        use std::fs::File;
        use std::io::Write;

        let mut file = File::create(filename)?;

        // Write header
        write!(file, "time_step,total_population,number_in_hospital,number_severely_immunosuppressed,number_with_sepsis,total_currently_infected,infected_10_days_count,infected_30_days_count,total_with_resistance,currently_taking_drug_count,currently_infected_and_on_drug_count,taking_two_drugs_count,newly_infected_count,newly_infected_with_resistance_count,newly_infected_past_year,total_deaths,deaths_background,deaths_sepsis,deaths_drug_toxicity,deaths_past_year,deaths_background_past_year,deaths_sepsis_past_year,deaths_drug_toxicity_past_year,num_age_0_5,num_age_6_14,num_age_15_49,num_age_50_79,num_age_80plus,num_with_any_bacteria_microbiome")?;
        // Add per-bacteria infection columns
        for bacteria in BACTERIA_LIST.iter() {
            write!(file, ",{}_currently_infected", bacteria.replace(" ", "_"))?;
        }
        // Add per-drug currently on drug columns
        for drug in DRUG_SHORT_NAMES.iter() {
            write!(file, ",{}_currently_on_drug", drug.replace(" ", "_"))?;
        }
        // Add per-bacteria, per-drug MIC < 2 columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                write!(file, ",{}_infected_and_mic_lt2_{}", bacteria.replace(" ", "_"), drug)?;
            }
        }
        writeln!(file)?;

        // Write data
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();
        for summary in &self.summary_log {
            write!(file, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}", 
                summary.time_step, 
                summary.total_population,
                summary.number_in_hospital,
                summary.number_severely_immunosuppressed,
                summary.number_with_sepsis,
                summary.total_currently_infected,
                summary.infected_10_days_count,
                summary.infected_30_days_count,
                summary.total_with_resistance,
                summary.currently_taking_drug_count,
                summary.currently_infected_and_on_drug_count,
                summary.taking_two_drugs_count,
                summary.newly_infected_count,
                summary.newly_infected_with_resistance_count,
                summary.newly_infected_past_year,
                summary.total_deaths,
                summary.deaths_background,
                summary.deaths_sepsis,
                summary.deaths_drug_toxicity,
                summary.deaths_past_year,
                summary.deaths_background_past_year,
                summary.deaths_sepsis_past_year,
                summary.deaths_drug_toxicity_past_year,
                summary.num_age_0_5,
                summary.num_age_6_14,
                summary.num_age_15_49,
                summary.num_age_50_79,
                summary.num_age_80plus,
                summary.num_with_any_bacteria_microbiome,
            )?;
            // Output per-bacteria infection counts
            for b_idx in 0..num_bacteria {
                write!(file, ",{}", summary.infections_by_bacteria[b_idx])?;
            }
            // Output per-drug currently on drug counts
            for d_idx in 0..num_drugs {
                write!(file, ",{}", summary.currently_on_drug_by_drug[d_idx])?;
            }
            // Output per-bacteria, per-drug infection and resistance counts
            for b_idx in 0..num_bacteria {
                for d_idx in 0..num_drugs {
                    let idx = b_idx * num_drugs + d_idx;
                    write!(file, ",{}", summary.infected_and_standardized_mic_lt2_by_bacteria_drug[idx])?;
                }
            }
            writeln!(file)?;
        }

        println!("Summary data exported to {}", filename);
        Ok(())
    }
}