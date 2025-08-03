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
use std::sync::atomic::{AtomicUsize, Ordering};

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
}

impl Simulation {
    /// Create a new Simulation instance with initialized population and lookup tables.
    ///
    /// Initializes population, bacteria/drug indices, and cross-resistance groups.
    pub fn new(population_size: usize, time_steps: usize) -> Self {
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

        Simulation { // Constructs and returns a new Simulation instance with the initialized population, time steps, and other data structures.
            population,
            time_steps,
            bacteria_indices,
            drug_indices,
            cross_resistance_groups, 
            current_majority_r_positive_values_by_combo: HashMap::new(), // Initialize empty
            summary_log: Vec::new(), // Initialize empty log
        }
    }

    pub fn run(&mut self) {
        // public function named run, which executes the simulation for the specified number of time steps.

        println!(" ");
        println!("--- starting to run over time steps");
        println!(" ");

        for t in 0..self.time_steps {
            // --- Per-bacteria, per-drug infection and resistance counts (parallel-friendly) ---
            let num_bacteria = BACTERIA_LIST.len();
            let num_drugs = DRUG_SHORT_NAMES.len();
            // Each thread will return a (infected, infected_and_standardized_mic_gt5) Vec<Vec<usize>>
            let infected_and_standardized_mic_lt2_by_bacteria_drug: Vec<usize> = self.population.individuals.par_iter()
                .filter(|individual| individual.date_of_death.is_none())
                .map(|individual| {
                    let mut infected_and_standardized_mic_lt2 = vec![0usize; num_bacteria * num_drugs];
                    for b_idx in 0..num_bacteria {
                        if individual.level[b_idx] > 0.001 {
                            for d_idx in 0..num_drugs {
                                let resistance_data = &individual.resistances[b_idx][d_idx];
                                let bacteria_name = BACTERIA_LIST[b_idx];
                                let drug_name = DRUG_SHORT_NAMES[d_idx];
                                let potency_key = format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug_name, bacteria_name);
                                let potency = crate::config::PARAMETERS.get(&potency_key).copied().unwrap_or(0.01);
                                let standardized_mic = 1.0 / ((1.0 - resistance_data.majority_r) * potency);
                                if standardized_mic < 2.0 {
                                    infected_and_standardized_mic_lt2[b_idx * num_drugs + d_idx] = 1;
                                }
                            }
                        }
                    }
                    infected_and_standardized_mic_lt2
                })
                .reduce(
                    || vec![0usize; num_bacteria * num_drugs],
                    |mut acc, x| {
                        for i in 0..acc.len() {
                            acc[i] += x[i];
                        }
                        acc
                    }
                );
//          println!("simulation.rs time step: {}", t);

            // Counter for intersection of currently infected AND on any drug
            let currently_infected_and_on_drug_count = AtomicUsize::new(0);

            // Counter for people with any presence_microbiome=true
            let num_with_any_bacteria_microbiome = AtomicUsize::new(0);

            // Use previous time step's resistance data for new acquisitions
            let previous_majority_r_positive_values_by_combo = if t == 0 {
                HashMap::new() // Empty for first time step
            } else {
            // Use the data collected in the previous iteration
                std::mem::take(&mut self.current_majority_r_positive_values_by_combo)
            };

            // Initialize counters and data structures for this time step
            let new_majority_r_positive_values_by_combo: HashMap<(usize, bool, usize, usize), Vec<f64>> = HashMap::new();
            let log_majority_r_positive_counts: Vec<Vec<AtomicUsize>> = (0..BACTERIA_LIST.len())
                .map(|_| (0..DRUG_SHORT_NAMES.len()).map(|_| AtomicUsize::new(0)).collect())
                .collect();

            // All counters use AtomicUsize for thread-safe parallel processing
            let log_infections_by_bacteria: Vec<AtomicUsize> = (0..BACTERIA_LIST.len()).map(|_| AtomicUsize::new(0)).collect();
            let log_resistance_counts: Vec<Vec<AtomicUsize>> = (0..BACTERIA_LIST.len())
                .map(|_| (0..DRUG_SHORT_NAMES.len()).map(|_| AtomicUsize::new(0)).collect())
                .collect();
            let log_total_deaths = AtomicUsize::new(0);
            let log_deaths_background = AtomicUsize::new(0);
            let log_deaths_sepsis = AtomicUsize::new(0);
            let log_deaths_drug_toxicity = AtomicUsize::new(0);
            let currently_taking_drug_count = AtomicUsize::new(0);
            let infected_10_days_count = AtomicUsize::new(0);
            let infected_30_days_count = AtomicUsize::new(0);
            let taking_two_drugs_count = AtomicUsize::new(0);
            let number_in_hospital = AtomicUsize::new(0);
            let number_severely_immunosuppressed = AtomicUsize::new(0);
            let number_with_sepsis = AtomicUsize::new(0);
            let newly_infected_count = AtomicUsize::new(0);
            let total_currently_infected = AtomicUsize::new(0);
            let total_with_resistance = AtomicUsize::new(0);

            // Per-drug: count number of living people currently taking each drug
            let currently_on_drug_by_drug = (0..num_drugs).map(|_| AtomicUsize::new(0)).collect::<Vec<_>>();

            // Single pass: apply rules and collect all statistics
            self.population.individuals.par_iter_mut().for_each(|individual| {
                // Apply rules first
                apply_rules(
                    individual,
                    t,
                    &previous_majority_r_positive_values_by_combo,
                    &self.bacteria_indices,
                    &self.drug_indices,
                    &self.cross_resistance_groups,
                );

                // Count deaths in this time step only, with cause tracking
                if let Some(death_time) = individual.date_of_death {
                    if death_time == t {
                        log_total_deaths.fetch_add(1, Ordering::Relaxed);
                        
                        // Count by cause of death
                        if let Some(ref cause) = individual.cause_of_death {
                            match cause.as_str() {
                                "background_mortality" => {
                                    log_deaths_background.fetch_add(1, Ordering::Relaxed);
                                },
                                "sepsis_related" => {
                                    log_deaths_sepsis.fetch_add(1, Ordering::Relaxed);
                                },
                                "drug_toxicity_related" => {
                                    log_deaths_drug_toxicity.fetch_add(1, Ordering::Relaxed);
                                },
                                _ => {
                                    // Unknown cause, count as background
                                    log_deaths_background.fetch_add(1, Ordering::Relaxed);
                                }
                            }
                        } else {
                            // No cause specified, count as background
                            log_deaths_background.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }

                // Only collect statistics for alive individuals
                if individual.date_of_death.is_none() {
                    // For each drug, count if this individual is currently taking it
                    for (d_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
                        if is_using {
                            currently_on_drug_by_drug[d_idx].fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    // Check if any presence_microbiome is true
                    if individual.presence_microbiome.iter().any(|&x| x) {
                        num_with_any_bacteria_microbiome.fetch_add(1, Ordering::Relaxed);
                    }
                    // Count individuals currently taking a drug
                    if individual.cur_use_drug.iter().any(|&is_using| is_using) {
                        currently_taking_drug_count.fetch_add(1, Ordering::Relaxed);
                    }

                    // Count individuals who are both currently infected (any bacteria) AND on any drug
                    let is_on_any_drug = individual.cur_use_drug.iter().any(|&is_using| is_using);
                    let is_currently_infected = BACTERIA_LIST.iter().enumerate().any(|(b_idx, _)| individual.level[b_idx] > 0.001);
                    if is_on_any_drug && is_currently_infected {
                        currently_infected_and_on_drug_count.fetch_add(1, Ordering::Relaxed);
                    }

                    // Count individuals infected for >10 days and >30 days
                    let mut individual_max_infection_duration = 0;
                    let mut individual_has_any_infection = false;
                    let mut individual_has_any_r_positive = false;
                    let mut was_newly_infected = false;
                    
                    for (b_idx, _) in BACTERIA_LIST.iter().enumerate() {
                        if individual.level[b_idx] > 0.001 {
                            individual_has_any_infection = true;
                            
                            // Log infections
                            log_infections_by_bacteria[b_idx].fetch_add(1, Ordering::Relaxed);

                            // Check infection duration
                            let days_since_infection = t as i32 - individual.date_last_infected[b_idx];
                            if days_since_infection > individual_max_infection_duration {
                                individual_max_infection_duration = days_since_infection;
                            }

                            // Check if newly infected at this time step
                            if individual.date_last_infected[b_idx] == t as i32 {
                                was_newly_infected = true;
                            }

                // Count resistance
                for (d_idx, _) in DRUG_SHORT_NAMES.iter().enumerate() {
                    let resistance_data = &individual.resistances[b_idx][d_idx];
                    if resistance_data.majority_r > 0.0 {
                        log_resistance_counts[b_idx][d_idx].fetch_add(1, Ordering::Relaxed);
                        log_majority_r_positive_counts[b_idx][d_idx].fetch_add(1, Ordering::Relaxed);
                    }
                    if resistance_data.any_r > 0.0 {
                        individual_has_any_r_positive = true;
                    }
                }
                        }
                    }
                    
                    if individual_has_any_infection {
                        total_currently_infected.fetch_add(1, Ordering::Relaxed);
                    }
                    
                    if individual_has_any_r_positive {
                        total_with_resistance.fetch_add(1, Ordering::Relaxed);
                    }
                    
                    if individual_max_infection_duration > 10 {
                        infected_10_days_count.fetch_add(1, Ordering::Relaxed);
                    }
                    if individual_max_infection_duration > 30 {
                        infected_30_days_count.fetch_add(1, Ordering::Relaxed);
                    }

                    if was_newly_infected {
                        newly_infected_count.fetch_add(1, Ordering::Relaxed);
                    }

                    // Count individuals taking two drugs
                    let active_drug_count = individual.cur_use_drug.iter().filter(|&&is_using| is_using).count();
                    if active_drug_count >= 2 {
                        taking_two_drugs_count.fetch_add(1, Ordering::Relaxed);
                    }
                    
                    // Count individuals in hospital
                    if individual.hospital_status.is_hospitalized() {
                        number_in_hospital.fetch_add(1, Ordering::Relaxed);
                    }
                    
                    // Count individuals severely immunosuppressed
                    if individual.is_severely_immunosuppressed {
                        number_severely_immunosuppressed.fetch_add(1, Ordering::Relaxed);
                    }
                    
                    // Count individuals with sepsis
                    if individual.sepsis.iter().any(|&has_sepsis| has_sepsis) {
                        number_with_sepsis.fetch_add(1, Ordering::Relaxed);
                    }
                }
            });

            // Collect remaining statistics that need sequential access
            // No need for sequential pass for per-bacteria/drug majority_r counts

            // Store for next iteration
            self.current_majority_r_positive_values_by_combo = new_majority_r_positive_values_by_combo;

            // Count living population (age >= 0 and no date_of_death)
            let living_population = self.population.individuals.iter()
                .filter(|individual| individual.age >= 0 && individual.date_of_death.is_none())
                .count();


            // Count living individuals in each age group
            let (mut num_age_0_5, mut num_age_6_14, mut num_age_15_49, mut num_age_50_79, mut num_age_80plus) = (0, 0, 0, 0, 0);
            for individual in self.population.individuals.iter() {
                if individual.date_of_death.is_none() {
            let age_years = individual.age as f64 / 365.0;
            if age_years >= 0.0 && age_years < 6.0 {
                num_age_0_5 += 1;
            } else if age_years >= 6.0 && age_years < 15.0 {
                num_age_6_14 += 1;
            } else if age_years >= 15.0 && age_years < 50.0 {
                num_age_15_49 += 1;
            } else if age_years >= 50.0 && age_years < 80.0 {
                num_age_50_79 += 1;
            } else if age_years >= 80.0 {
                num_age_80plus += 1;
            }
                }
            }

            // Create summary for this time step
            let infected_10_count = infected_10_days_count.load(Ordering::Relaxed);
            let infected_30_count = infected_30_days_count.load(Ordering::Relaxed);

            let summary = TimeStepSummary {
        infected_and_standardized_mic_lt2_by_bacteria_drug: infected_and_standardized_mic_lt2_by_bacteria_drug.clone(),
        currently_on_drug_by_drug: currently_on_drug_by_drug.iter().map(|x| x.load(Ordering::Relaxed)).collect(),
                num_age_0_5,
                num_age_6_14,
                num_age_15_49,
                num_age_50_79,
                num_age_80plus,
                num_with_any_bacteria_microbiome: num_with_any_bacteria_microbiome.load(Ordering::Relaxed),
                time_step: t,
                total_population: living_population,
                number_in_hospital: number_in_hospital.load(Ordering::Relaxed),
                number_severely_immunosuppressed: number_severely_immunosuppressed.load(Ordering::Relaxed),
                number_with_sepsis: number_with_sepsis.load(Ordering::Relaxed),
                newly_infected_count: newly_infected_count.load(Ordering::Relaxed),
                total_currently_infected: total_currently_infected.load(Ordering::Relaxed),
                total_with_resistance: total_with_resistance.load(Ordering::Relaxed),
                infected_10_days_count: infected_10_count,
                infected_30_days_count: infected_30_count,
                currently_taking_drug_count: currently_taking_drug_count.load(Ordering::Relaxed),
                taking_two_drugs_count: taking_two_drugs_count.load(Ordering::Relaxed),
                infections_by_bacteria: log_infections_by_bacteria.iter().map(|x| x.load(Ordering::Relaxed)).collect(),
                resistance_by_bacteria_drug: log_resistance_counts.iter().map(|row| 
                    row.iter().map(|x| x.load(Ordering::Relaxed)).collect()
                ).collect(),
                total_deaths: log_total_deaths.load(Ordering::Relaxed),
                deaths_background: log_deaths_background.load(Ordering::Relaxed),
                deaths_sepsis: log_deaths_sepsis.load(Ordering::Relaxed),
                deaths_drug_toxicity: log_deaths_drug_toxicity.load(Ordering::Relaxed),
                // Rolling 1-year (365 days) death counts
                deaths_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.total_deaths)
                        .sum::<usize>()
                        + log_total_deaths.load(Ordering::Relaxed)
                        - self.summary_log.last().map_or(0, |s| s.total_deaths)
                },
                deaths_background_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.deaths_background)
                        .sum::<usize>()
                        + log_deaths_background.load(Ordering::Relaxed)
                        - self.summary_log.last().map_or(0, |s| s.deaths_background)
                },
                deaths_sepsis_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.deaths_sepsis)
                        .sum::<usize>()
                        + log_deaths_sepsis.load(Ordering::Relaxed)
                        - self.summary_log.last().map_or(0, |s| s.deaths_sepsis)
                },
                deaths_drug_toxicity_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.deaths_drug_toxicity)
                        .sum::<usize>()
                        + log_deaths_drug_toxicity.load(Ordering::Relaxed)
                        - self.summary_log.last().map_or(0, |s| s.deaths_drug_toxicity)
                },
                newly_infected_past_year: {
                    let start = if self.summary_log.len() >= 365 { self.summary_log.len() - 365 } else { 0 };
                    self.summary_log[start..]
                        .iter()
                        .map(|s| s.newly_infected_count)
                        .sum::<usize>()
                        + newly_infected_count.load(Ordering::Relaxed)
                        - self.summary_log.last().map_or(0, |s| s.newly_infected_count)
                },
                currently_infected_and_on_drug_count: currently_infected_and_on_drug_count.load(Ordering::Relaxed),
            };

            // Comprehensive print block for individual 0
            let individual_0 = &self.population.individuals[0];
            println!("--- Individual 0 full state ---");
            println!("id: {}", individual_0.id);
            println!("age (days): {}", individual_0.age);
            println!("sex_at_birth: {}", individual_0.sex_at_birth);
            println!("region_living: {:?}", individual_0.region_living);
            println!("region_cur_in: {:?}", individual_0.region_cur_in);
            println!("current_infection_related_death_risk: {:.4}", individual_0.current_infection_related_death_risk);
            println!("background_all_cause_mortality_rate: {:.4}", individual_0.background_all_cause_mortality_rate);
            println!("sexual_contact_level: {:.4}", individual_0.sexual_contact_level);
            println!("airborne_contact_level_with_adults: {:.4}", individual_0.airborne_contact_level_with_adults);
            println!("airborne_contact_level_with_children: {:.4}", individual_0.airborne_contact_level_with_children);
            println!("oral_exposure_level: {:.4}", individual_0.oral_exposure_level);
            println!("mosquito_exposure_level: {:.4}", individual_0.mosquito_exposure_level);
            println!("current_toxicity: {:.4}", individual_0.current_toxicity);
            println!("mortality_risk_current_toxicity: {:.4}", individual_0.mortality_risk_current_toxicity);
            println!("hospital_status: {:?}", individual_0.hospital_status);
            println!("is_severely_immunosuppressed: {:?}", individual_0.is_severely_immunosuppressed);
            println!("date_of_death: {:?}", individual_0.date_of_death);
            // Arrays
            println!("level: {:?}", individual_0.level);
            println!("immune_resp: {:?}", individual_0.immune_resp);
            println!("presence_microbiome: {:?}", individual_0.presence_microbiome);
            println!("cur_level_drug: {:?}", individual_0.cur_level_drug);
            println!("cur_use_drug: {:?}", individual_0.cur_use_drug);
            println!("ever_taken_drug: {:?}", individual_0.ever_taken_drug);
            println!("date_last_infected: {:?}", individual_0.date_last_infected);
            println!("cur_infection_from_environment: {:?}", individual_0.cur_infection_from_environment);
            println!("infection_hospital_acquired: {:?}", individual_0.infection_hospital_acquired);
            println!("test_identified_infection: {:?}", individual_0.test_identified_infection);
            println!("sepsis: {:?}", individual_0.sepsis);
            // Per-bacteria/drug resistance data
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                for (d_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                    let resistance = &individual_0.resistances[b_idx][d_idx];
                    println!(
                        "Resistance for bacteria {} and drug {}: any_r = {:.4}, activity_r = {:.4}, majority_r = {:.4}",
                        bacteria_name, drug_name, resistance.any_r, resistance.activity_r, resistance.majority_r
                    );
                }
            }

            // Print summary statistics for this time step
            println!("Time step {}: {} newly infected, {} deaths, {} with resistance", 
                t, 
                summary.newly_infected_count, 
                summary.total_deaths, 
                summary.total_with_resistance
            );

            self.summary_log.push(summary);

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
        write!(file, "time_step,total_population,number_in_hospital,number_severely_immunosuppressed,number_with_sepsis,total_currently_infected,infected_10_days_count,infected_30_days_count,total_with_resistance,currently_taking_drug_count,currently_infected_and_on_drug_count,taking_two_drugs_count,newly_infected_count,newly_infected_past_year,total_deaths,deaths_background,deaths_sepsis,deaths_drug_toxicity,deaths_past_year,deaths_background_past_year,deaths_sepsis_past_year,deaths_drug_toxicity_past_year,num_age_0_5,num_age_6_14,num_age_15_49,num_age_50_79,num_age_80plus,num_with_any_bacteria_microbiome")?;
        // Add per-bacteria infection columns
        for bacteria in BACTERIA_LIST.iter() {
            write!(file, ",{}_currently_infected", bacteria.replace(" ", "_"))?;
        }
        // Add per-drug currently on drug columns
        for drug in DRUG_SHORT_NAMES.iter() {
            write!(file, ",{}_currently_on_drug", drug.replace(" ", "_"))?;
        }
        // Add per-bacteria, per-drug MIC < 2 columns
        for (b_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
            for drug in DRUG_SHORT_NAMES.iter() {
                write!(file, ",{}_infected_and_mic_lt2_{}", bacteria.replace(" ", "_"), drug)?;
            }
        }
        writeln!(file)?;

        // Write data
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();
        for summary in &self.summary_log {
            write!(file, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}", 
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