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


use crate::simulation::population::{Population, BACTERIA_LIST, DRUG_SHORT_NAMES, ResistanceMechanism, Region};
use crate::rules::apply_rules;
use crate::config::{self, get_global_param}; // Import the config module and get_global_param function
use std::collections::HashMap;
use rayon::prelude::*;
// Removed most atomics by using thread-local aggregation; retain no atomic imports here.
use std::time::Instant;
use std::io::Write;

// Helper function to convert Region enum to array index
fn region_to_index(region: Region) -> usize {
    match region {
        Region::NorthAmerica => 0,
        Region::SouthAmerica => 1,
        Region::Africa => 2,
        Region::Asia => 3,
        Region::Europe => 4,
        Region::Oceania => 5,
        Region::Home => panic!("Home should be resolved to actual region before calling this function"),
    }
}

// Helper function to get the effective region (resolves Home to actual home region)
fn get_effective_region(individual: &crate::simulation::population::Individual) -> Region {
    match individual.region_cur_in {
        Region::Home => individual.region_living,
        other => other,
    }
}

// Compact structure for time step summary data
#[allow(dead_code)]
#[derive(Clone)]
// Summary statistics for each simulation time step.
//
// Captures population-level and per-bacteria/drug summary metrics for each time step.
pub struct TimeStepSummary {
    // per-bacteria count of people on any drug (infected with each bacteria and on at least one drug)
    pub infected_and_on_any_drug_by_bacteria: Vec<usize>,
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
    pub number_with_sepsis_by_bacteria: Vec<usize>, // per-bacteria counts of people with sepsis
    pub new_sepsis_cases_by_bacteria: Vec<usize>, // per-bacteria counts of people who developed sepsis this timestep
    pub infections_by_bacteria: Vec<usize>, // indexed by bacteria
    pub deaths_by_bacteria: Vec<usize>, // indexed by bacteria
    pub resistance_by_bacteria_drug: Vec<Vec<usize>>, // [bacteria][drug] counts
    /// per-bacteria sum of activity_r values for all individuals (float, indexed by bacteria)
    pub activity_r_sum_by_bacteria: Vec<f64>,
    pub newly_infected_count: usize, // Number of people newly infected this time step
    pub newly_infected_with_resistance_count: usize, // Number of newly infected people who acquired resistance
    pub new_drug_initiations_count: usize, // Number of people who started any new drug this time step
    pub new_drug_initiations_count_infected: usize, // Number of currently infected people who started any new drug this time step
    pub newly_infected_by_bacteria_region: Vec<usize>, // [bacteria * region] = new active infections this timestep by bacteria and home region
    pub deaths_infected_by_bacteria_region: Vec<usize>, // [bacteria * region] = deaths this timestep of people currently infected with bacteria by home region
    pub newly_infected_past_year: usize, // Rolling 1-year (365 days) newly infected count
    pub currently_infected_and_on_drug_count: usize, // intersection of currently infected AND on any drug
    pub num_age_0_5: usize,
    pub num_age_6_14: usize,
    pub num_age_15_49: usize,
    pub num_age_50_79: usize,
    pub num_age_80plus: usize,
    pub num_with_any_bacteria_microbiome: usize, // number of people with any presence_microbiome=true
    pub presence_microbiome_by_bacteria: Vec<usize>, // per-bacteria counts of people with this bacteria in microbiome
    pub presence_microbiome_by_bacteria_by_region: Vec<Vec<usize>>, // [bacteria][region] counts of people with bacteria in microbiome by region
    pub infected_with_test_identified_by_bacteria: Vec<usize>, // per-bacteria counts of infected people with test_identified_infection = true
    pub infected_with_test_for_resistance_by_bacteria: Vec<usize>, // per-bacteria counts of infected people with test_for_resistance = true
    
    // Drug failure tracking: day 5 post-drug-initiation events by bacteria and region
    pub drug_failure_events_by_bacteria_region: Vec<Vec<usize>>, // [bacteria][region] - numerator: day 5, on drug, still infected
    pub drug_treatment_day5_events_by_bacteria_region: Vec<Vec<usize>>, // [bacteria][region] - denominator: day 5 post-drug-initiation

    // per-bacteria, per-drug infection and resistance counts (flat, len = bacteria * drugs)
    pub infected_and_standardized_mic_lt2_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug currently on drug counts (flat, len = bacteria * drugs)
    pub currently_on_drug_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug microbiome_r > 0 counts (flat, len = bacteria * drugs)
    pub microbiome_r_positive_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug any_r sum values for infected individuals (flat, len = bacteria * drugs)
    pub any_r_sum_by_bacteria_drug: Vec<f64>,
    
    // per-bacteria, per-drug any_r sum values for hospital-acquired infected individuals (flat, len = bacteria * drugs)
    pub any_r_sum_by_bacteria_drug_hospital: Vec<f64>,

    // per-bacteria, per-drug counts of infected individuals with any_r > 0 (flat, len = bacteria * drugs)
    pub infected_with_any_r_positive_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug MIC sum values for infected individuals (flat, len = bacteria * drugs)
    pub mic_sum_by_bacteria_drug: Vec<f64>,

    // per-region any_r sum values pooled across all bacteria and drugs (indexed by region)
    pub any_r_sum_by_region: Vec<f64>,
    
    // per-region count of infected individuals (for calculating mean) (indexed by region)
    pub infected_count_by_region: Vec<usize>,

    // per-drug currently on drug counts (indexed by drug)
    pub currently_on_drug_by_drug: Vec<usize>,

    // per-bacteria, per-resistance-mechanism counts (flat, len = bacteria * mechanisms)
    // infected_with_bacteria_and_mechanism[bacteria_idx * num_mechanisms + mechanism_idx] = count
    pub infected_with_bacteria_and_mechanism: Vec<usize>,
    
    // counts of newly acquired resistance by acquisition type this timestep per bacteria-drug combination
    // Each Vec has length = num_bacteria * num_drugs, indexed as [bacteria_idx * num_drugs + drug_idx]
    pub new_resistance_at_infection_community_by_bacteria_drug: Vec<usize>,
    pub new_resistance_at_infection_env_by_bacteria_drug: Vec<usize>,
    pub new_resistance_hgt_by_bacteria_drug: Vec<usize>,
    pub new_resistance_from_microbiome_r_by_bacteria_drug: Vec<usize>,
    
    // infection resolution tracking: counts by bacteria and resolution type
    // Each Vec has length = num_bacteria * num_resolution_types, indexed as [bacteria_idx * num_resolution_types + resolution_type_idx]
    pub infection_resolution_immune_clearance_by_bacteria: Vec<usize>,
    pub infection_resolution_drug_assisted_clearance_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_sepsis_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_background_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_toxicity_by_bacteria: Vec<usize>,
    
    // day-7 drug initiation tracking: counts by bacteria
    pub day_7_evaluations_by_bacteria: Vec<usize>,        // [bacteria_idx] = number of post-infection evaluations (configurable timing)
    pub day_7_drug_used_by_bacteria: Vec<usize>,          // [bacteria_idx] = number where drug was used by day 7
    
    // syndrome tracking: counts by syndrome (1-10)
    pub infected_by_syndrome: Vec<usize>,                 // [syndrome_idx] = number of infected individuals with this syndrome (first infection only)
    
    // bacteria-specific syndrome tracking: counts by bacteria and syndrome (bacteria * 10 syndromes)
    // [bacteria_idx * 10 + syndrome_idx] = number of infected individuals with this bacteria and syndrome
    pub infected_by_syndrome_by_bacteria: Vec<usize>,     // [bacteria][syndrome] = number of infected individuals with this bacteria and syndrome
    
    // regional population tracking: counts by region (6 regions: NorthAmerica, SouthAmerica, Africa, Asia, Europe, Oceania)
    pub living_population_by_region: Vec<usize>,          // [region_idx] = number of living individuals currently in this region
    
    // regional hospital population tracking: counts by region (6 regions)
    pub hospital_population_by_region: Vec<usize>,        // [region_idx] = number of individuals currently in hospital in this region
    
    // hospital-acquired new infection tracking: counts by bacteria and region (bacteria * 6 regions)
    pub newly_infected_hospital_by_bacteria_region: HashMap<(usize, usize), usize>, // (bacteria_idx, region_idx) = count of new hospital infections
    
    // regional age distribution tracking: counts by region and age group (6 regions * 5 age groups = 30 values)
    // [region_idx * 5 + age_group_idx] where age_group_idx: 0=0-5, 1=6-14, 2=15-49, 3=50-79, 4=80+
    pub age_distribution_by_region: Vec<usize>,           // [region][age_group] = number of living individuals in this region and age group
    
    // regional death tracking: counts by region and death type (6 regions * 3 death types = 18 values)
    // [region_idx * 3 + death_type_idx] where death_type_idx: 0=background, 1=sepsis, 2=drug_toxicity
    pub deaths_by_region: Vec<usize>,                     // [region][death_type] = number of deaths in this region by cause
    
    // age-specific death tracking by region: counts by region, age group, and death type (6 regions * 5 age groups * 3 death types = 90 values)
    // [region_idx * 15 + age_group_idx * 3 + death_type_idx] where age_group_idx: 0=0-5, 1=6-14, 2=15-49, 3=50-79, 4=80+
    pub deaths_by_region_age: Vec<usize>,                 // [region][age_group][death_type] = number of deaths
    
    // syndrome population by region: counts by syndrome and region (10 syndromes * 6 regions = 60 values)
    // [syndrome_idx * 6 + region_idx] where syndrome_idx: 0-9 (syndromes 1-10), region_idx: 0-5
    pub syndrome_population_by_region: Vec<usize>,        // [syndrome][region] = number of individuals with this syndrome in this region
    
    // syndrome deaths from sepsis by region: counts by syndrome and region (10 syndromes * 6 regions = 60 values)
    // [syndrome_idx * 6 + region_idx] where syndrome_idx: 0-9 (syndromes 1-10), region_idx: 0-5
    pub syndrome_deaths_sepsis_by_region: Vec<usize>,     // [syndrome][region] = number of sepsis deaths with this syndrome in this region
    
    // regional drug usage tracking: counts by region and drug (6 regions * num_drugs values)
    // [region_idx * num_drugs + drug_idx] = number of people currently taking this drug in this region
    pub currently_on_drug_by_region_drug: Vec<usize>,     // [region][drug] = number of people currently on drug in region
    
    // polypharmacy tracking: counts of people taking 1, 2, or ≥3 drugs simultaneously
    pub people_on_1_drug: usize,                         // number of people taking exactly 1 drug
    pub people_on_2_drugs: usize,                        // number of people taking exactly 2 drugs  
    pub people_on_3plus_drugs: usize,                    // number of people taking 3 or more drugs
    
    // treatment failure tracking: people currently on drug + infected + previously failed treatment
    pub infected_on_drug_with_previous_failure: usize,   // numerator: people currently infected, on drug, with previous treatment failure
    
    // drug score tracking: aggregate statistics for clinical guideline debugging
    pub drug_selection_count_by_bacteria: Vec<usize>,    // [bacteria_idx] = number of drug selections for this bacteria this timestep
    pub drug_score_sums_by_bacteria_drug: Vec<f64>,      // [bacteria_idx * num_drugs + drug_idx] = sum of drug scores for this bacteria-drug combo this timestep
    
    // current number of drugs tracking: histogram of people by number of drugs they're taking
    pub people_by_drug_count: Vec<usize>,                // [0] = people on 0 drugs, [1] = people on 1 drug, etc.
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
                infected_and_on_any_drug_by_bacteria: Vec<usize>,
                mic_lt2_counts: Vec<usize>,
                currently_on_drug_by_bacteria_drug: Vec<usize>,
                microbiome_r_positive_by_bacteria_drug: Vec<usize>,
                infections_by_bacteria: Vec<usize>,
                deaths_by_bacteria: Vec<usize>,
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
                number_with_sepsis_by_bacteria: Vec<usize>,
                new_sepsis_cases_by_bacteria: Vec<usize>,
                newly_infected_count: usize,
                newly_infected_with_resistance_count: usize,
                new_drug_initiations_count: usize,
                new_drug_initiations_count_infected: usize,
                newly_infected_by_bacteria_region: Vec<usize>,
                deaths_infected_by_bacteria_region: Vec<usize>,
                total_currently_infected: usize,
                total_with_resistance: usize,
                currently_infected_and_on_drug_count: usize,
                num_with_any_bacteria_microbiome: usize,
                presence_microbiome_by_bacteria: Vec<usize>,
                presence_microbiome_by_bacteria_by_region: Vec<Vec<usize>>,
                drug_failure_events_by_bacteria_region: Vec<Vec<usize>>,
                drug_treatment_day5_events_by_bacteria_region: Vec<Vec<usize>>,
                infected_with_test_identified_by_bacteria: Vec<usize>,
                infected_with_test_for_resistance_by_bacteria: Vec<usize>,
                // Integrated previously sequential counts:
                living_population: usize,
                num_age_0_5: usize,
                num_age_6_14: usize,
                num_age_15_49: usize,
                num_age_50_79: usize,
                num_age_80plus: usize,
                /// per-bacteria sum of activity_r values for all individuals (float, indexed by bacteria)
                activity_r_sum_by_bacteria: Vec<f64>,
                /// per-bacteria, per-drug sum of any_r values for infected individuals (float, indexed by bacteria * drugs)
                any_r_sum_by_bacteria_drug: Vec<f64>,
                /// per-bacteria, per-drug sum of any_r values for hospital-acquired infected individuals (float, indexed by bacteria * drugs)
                any_r_sum_by_bacteria_drug_hospital: Vec<f64>,
                /// per-bacteria, per-drug counts of infected individuals with any_r > 0 (flat, len = bacteria * drugs)
                infected_with_any_r_positive_by_bacteria_drug: Vec<usize>,
                /// per-bacteria, per-drug sum of MIC values for infected individuals (flat, len = bacteria * drugs)
                mic_sum_by_bacteria_drug: Vec<f64>,
                /// per-region sum of any_r values pooled across all bacteria and drugs (indexed by region)
                any_r_sum_by_region: Vec<f64>,
                /// per-region count of infected individuals (for calculating mean) (indexed by region)
                infected_count_by_region: Vec<usize>,
                /// per-bacteria, per-resistance-mechanism counts (flat, len = bacteria * mechanisms)
                infected_with_bacteria_and_mechanism: Vec<usize>,
                /// counts of newly acquired resistance by acquisition type this timestep per bacteria-drug combination
                new_resistance_at_infection_community_by_bacteria_drug: Vec<usize>,
                new_resistance_at_infection_env_by_bacteria_drug: Vec<usize>,
                new_resistance_hgt_by_bacteria_drug: Vec<usize>,
                new_resistance_from_microbiome_r_by_bacteria_drug: Vec<usize>,
                /// infection resolution tracking: counts by bacteria and resolution type
                infection_resolution_immune_clearance_by_bacteria: Vec<usize>,
                infection_resolution_drug_assisted_clearance_by_bacteria: Vec<usize>,
                infection_resolution_death_from_sepsis_by_bacteria: Vec<usize>,
                infection_resolution_death_from_background_by_bacteria: Vec<usize>,
                infection_resolution_death_from_toxicity_by_bacteria: Vec<usize>,
                /// counts of infected individuals by syndrome (1-10)
                infected_by_syndrome: Vec<usize>,
                /// counts of infected individuals by bacteria and syndrome (bacteria * 10 syndromes)
                infected_by_syndrome_by_bacteria: Vec<usize>,
                /// living population count by region (6 regions)
                living_population_by_region: Vec<usize>,
                /// age distribution by region (6 regions * 5 age groups = 30 values)
                age_distribution_by_region: Vec<usize>,
                /// death tracking by region (6 regions * 3 death types = 18 values)
                deaths_by_region: Vec<usize>,
                /// age-specific death tracking by region (6 regions * 5 age groups * 3 death types = 90 values)
                deaths_by_region_age: Vec<usize>,
                /// drug usage by region (6 regions * num_drugs)
                currently_on_drug_by_region_drug: Vec<usize>,
                /// syndrome deaths from sepsis by region (10 syndromes * 6 regions = 60 values)
                syndrome_deaths_sepsis_by_region: Vec<usize>,
            }
            impl LocalTotals {
                fn new(num_bacteria: usize, num_drugs: usize, majority_r_capacity: usize) -> Self {
                    Self {
                        mic_lt2_counts: vec![0; num_bacteria * num_drugs],
                        currently_on_drug_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        microbiome_r_positive_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        infected_and_on_any_drug_by_bacteria: vec![0; num_bacteria],
                        infections_by_bacteria: vec![0; num_bacteria],
                        deaths_by_bacteria: vec![0; num_bacteria],
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
                        number_with_sepsis_by_bacteria: vec![0; num_bacteria],
                        new_sepsis_cases_by_bacteria: vec![0; num_bacteria],
                        newly_infected_count: 0,
                        newly_infected_with_resistance_count: 0,
                        new_drug_initiations_count: 0,
                        new_drug_initiations_count_infected: 0,
                        newly_infected_by_bacteria_region: vec![0; num_bacteria * 6], // bacteria * regions
                        deaths_infected_by_bacteria_region: vec![0; num_bacteria * 6], // bacteria * regions
                        total_currently_infected: 0,
                        total_with_resistance: 0,
                        currently_infected_and_on_drug_count: 0,
                        num_with_any_bacteria_microbiome: 0,
                        presence_microbiome_by_bacteria: vec![0; num_bacteria],
                        presence_microbiome_by_bacteria_by_region: vec![vec![0; 6]; num_bacteria], // [bacteria][region]
                        drug_failure_events_by_bacteria_region: vec![vec![0; 6]; num_bacteria], // [bacteria][region]
                        drug_treatment_day5_events_by_bacteria_region: vec![vec![0; 6]; num_bacteria], // [bacteria][region]
                        infected_with_test_identified_by_bacteria: vec![0; num_bacteria],
                        infected_with_test_for_resistance_by_bacteria: vec![0; num_bacteria],
                        living_population: 0,
                        num_age_0_5: 0,
                        num_age_6_14: 0,
                        num_age_15_49: 0,
                        num_age_50_79: 0,
                        num_age_80plus: 0,
                        activity_r_sum_by_bacteria: vec![0.0; num_bacteria],
                        any_r_sum_by_bacteria_drug: vec![0.0; num_bacteria * num_drugs],
                        any_r_sum_by_bacteria_drug_hospital: vec![0.0; num_bacteria * num_drugs],
                        infected_with_any_r_positive_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        mic_sum_by_bacteria_drug: vec![0.0; num_bacteria * num_drugs],
                        any_r_sum_by_region: vec![0.0; 6], // 6 regions: NorthAmerica, SouthAmerica, Africa, Asia, Europe, Oceania (excluding Home)
                        infected_count_by_region: vec![0; 6], // 6 regions
                        infected_with_bacteria_and_mechanism: vec![0; num_bacteria * ResistanceMechanism::all().len()],
                        new_resistance_at_infection_community_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        new_resistance_at_infection_env_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        new_resistance_hgt_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        new_resistance_from_microbiome_r_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        infection_resolution_immune_clearance_by_bacteria: vec![0; num_bacteria],
                        infection_resolution_drug_assisted_clearance_by_bacteria: vec![0; num_bacteria],
                        infection_resolution_death_from_sepsis_by_bacteria: vec![0; num_bacteria],
                        infection_resolution_death_from_background_by_bacteria: vec![0; num_bacteria],
                        infection_resolution_death_from_toxicity_by_bacteria: vec![0; num_bacteria],
                        infected_by_syndrome: vec![0; 10], // Syndromes 1-10
                        infected_by_syndrome_by_bacteria: vec![0; num_bacteria * 10], // bacteria * syndromes
                        living_population_by_region: vec![0; 6], // 6 regions: NorthAmerica, SouthAmerica, Africa, Asia, Europe, Oceania
                        age_distribution_by_region: vec![0; 6 * 5], // 6 regions * 5 age groups = 30 values
                        deaths_by_region: vec![0; 6 * 3], // 6 regions * 3 death types = 18 values
                        deaths_by_region_age: vec![0; 6 * 5 * 3], // 6 regions * 5 age groups * 3 death types = 90 values
                        currently_on_drug_by_region_drug: vec![0; 6 * num_drugs], // 6 regions * num_drugs
                        syndrome_deaths_sepsis_by_region: vec![0; 10 * 6], // 10 syndromes * 6 regions = 60 values
                    }
                }
                fn merge(&mut self, other: Self) {
                    for (a,b) in self.mic_lt2_counts.iter_mut().zip(other.mic_lt2_counts) { *a += b; }
                    for (a,b) in self.currently_on_drug_by_bacteria_drug.iter_mut().zip(other.currently_on_drug_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.microbiome_r_positive_by_bacteria_drug.iter_mut().zip(other.microbiome_r_positive_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.infected_and_on_any_drug_by_bacteria.iter_mut().zip(other.infected_and_on_any_drug_by_bacteria) { *a += b; }
                    for (a,b) in self.infections_by_bacteria.iter_mut().zip(other.infections_by_bacteria) { *a += b; }
                    for (a,b) in self.deaths_by_bacteria.iter_mut().zip(other.deaths_by_bacteria) { *a += b; }
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
                    for (a,b) in self.number_with_sepsis_by_bacteria.iter_mut().zip(other.number_with_sepsis_by_bacteria) { *a += b; }
                    for (a,b) in self.new_sepsis_cases_by_bacteria.iter_mut().zip(other.new_sepsis_cases_by_bacteria) { *a += b; }
                    self.newly_infected_count += other.newly_infected_count;
                    self.newly_infected_with_resistance_count += other.newly_infected_with_resistance_count;
                    self.new_drug_initiations_count += other.new_drug_initiations_count;
                    self.new_drug_initiations_count_infected += other.new_drug_initiations_count_infected;
                    for i in 0..self.newly_infected_by_bacteria_region.len() {
                        self.newly_infected_by_bacteria_region[i] += other.newly_infected_by_bacteria_region[i];
                    }
                    for i in 0..self.deaths_infected_by_bacteria_region.len() {
                        self.deaths_infected_by_bacteria_region[i] += other.deaths_infected_by_bacteria_region[i];
                    }
                    self.total_currently_infected += other.total_currently_infected;
                    self.total_with_resistance += other.total_with_resistance;
                    self.currently_infected_and_on_drug_count += other.currently_infected_and_on_drug_count;
                    self.num_with_any_bacteria_microbiome += other.num_with_any_bacteria_microbiome;
                    for (a,b) in self.presence_microbiome_by_bacteria.iter_mut().zip(other.presence_microbiome_by_bacteria) { *a += b; }
                    for (a_bact,b_bact) in self.presence_microbiome_by_bacteria_by_region.iter_mut().zip(other.presence_microbiome_by_bacteria_by_region) { 
                        for (a_reg,b_reg) in a_bact.iter_mut().zip(b_bact) { *a_reg += b_reg; } 
                    }
                    for (a_bact,b_bact) in self.drug_failure_events_by_bacteria_region.iter_mut().zip(other.drug_failure_events_by_bacteria_region) {
                        for (a_reg,b_reg) in a_bact.iter_mut().zip(b_bact) { *a_reg += b_reg; }
                    }
                    for (a_bact,b_bact) in self.drug_treatment_day5_events_by_bacteria_region.iter_mut().zip(other.drug_treatment_day5_events_by_bacteria_region) {
                        for (a_reg,b_reg) in a_bact.iter_mut().zip(b_bact) { *a_reg += b_reg; }
                    }
                    for (a,b) in self.infected_with_test_identified_by_bacteria.iter_mut().zip(other.infected_with_test_identified_by_bacteria) { *a += b; }
                    for (a,b) in self.infected_with_test_for_resistance_by_bacteria.iter_mut().zip(other.infected_with_test_for_resistance_by_bacteria) { *a += b; }
                    self.living_population += other.living_population;
                    self.num_age_0_5 += other.num_age_0_5;
                    self.num_age_6_14 += other.num_age_6_14;
                    self.num_age_15_49 += other.num_age_15_49;
                    self.num_age_50_79 += other.num_age_50_79;
                    self.num_age_80plus += other.num_age_80plus;
                    for (a,b) in self.activity_r_sum_by_bacteria.iter_mut().zip(other.activity_r_sum_by_bacteria) { *a += b; }
                    for (a,b) in self.any_r_sum_by_bacteria_drug.iter_mut().zip(other.any_r_sum_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.any_r_sum_by_bacteria_drug_hospital.iter_mut().zip(other.any_r_sum_by_bacteria_drug_hospital) { *a += b; }
                    for (a,b) in self.infected_with_any_r_positive_by_bacteria_drug.iter_mut().zip(other.infected_with_any_r_positive_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.mic_sum_by_bacteria_drug.iter_mut().zip(other.mic_sum_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.any_r_sum_by_region.iter_mut().zip(other.any_r_sum_by_region) { *a += b; }
                    for (a,b) in self.infected_count_by_region.iter_mut().zip(other.infected_count_by_region) { *a += b; }
                    for (a,b) in self.infected_with_bacteria_and_mechanism.iter_mut().zip(other.infected_with_bacteria_and_mechanism) { *a += b; }
                    for (a,b) in self.new_resistance_at_infection_community_by_bacteria_drug.iter_mut().zip(other.new_resistance_at_infection_community_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.new_resistance_at_infection_env_by_bacteria_drug.iter_mut().zip(other.new_resistance_at_infection_env_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.new_resistance_hgt_by_bacteria_drug.iter_mut().zip(other.new_resistance_hgt_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.new_resistance_from_microbiome_r_by_bacteria_drug.iter_mut().zip(other.new_resistance_from_microbiome_r_by_bacteria_drug) { *a += b; }
                    for (a,b) in self.infection_resolution_immune_clearance_by_bacteria.iter_mut().zip(other.infection_resolution_immune_clearance_by_bacteria) { *a += b; }
                    for (a,b) in self.infection_resolution_drug_assisted_clearance_by_bacteria.iter_mut().zip(other.infection_resolution_drug_assisted_clearance_by_bacteria) { *a += b; }
                    for (a,b) in self.infection_resolution_death_from_sepsis_by_bacteria.iter_mut().zip(other.infection_resolution_death_from_sepsis_by_bacteria) { *a += b; }
                    for (a,b) in self.infection_resolution_death_from_background_by_bacteria.iter_mut().zip(other.infection_resolution_death_from_background_by_bacteria) { *a += b; }
                    for (a,b) in self.infection_resolution_death_from_toxicity_by_bacteria.iter_mut().zip(other.infection_resolution_death_from_toxicity_by_bacteria) { *a += b; }
                    for (a,b) in self.infected_by_syndrome.iter_mut().zip(other.infected_by_syndrome) { *a += b; }
                    for (a,b) in self.infected_by_syndrome_by_bacteria.iter_mut().zip(other.infected_by_syndrome_by_bacteria) { *a += b; }
                    for (a,b) in self.living_population_by_region.iter_mut().zip(other.living_population_by_region) { *a += b; }
                    for (a,b) in self.age_distribution_by_region.iter_mut().zip(other.age_distribution_by_region) { *a += b; }
                    for (a,b) in self.deaths_by_region.iter_mut().zip(other.deaths_by_region) { *a += b; }
                    for (a,b) in self.deaths_by_region_age.iter_mut().zip(other.deaths_by_region_age) { *a += b; }
                    for (a,b) in self.currently_on_drug_by_region_drug.iter_mut().zip(other.currently_on_drug_by_region_drug) { *a += b; }
                    for (a,b) in self.syndrome_deaths_sepsis_by_region.iter_mut().zip(other.syndrome_deaths_sepsis_by_region) { *a += b; }
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
                        if individual.date_of_death.is_none() && individual.age >= 0 {
                                // Check if individual is infected with any bacteria
                                let is_infected = individual.level.iter().any(|&level| level > 0.001);
                                if is_infected {
                                    // Count infected individual by region (only once per individual)
                                    let effective_region = get_effective_region(individual);
                                    let region_idx = region_to_index(effective_region);
                                    lt.infected_count_by_region[region_idx] += 1;
                                }
                                
                                for b_idx in 0..num_bacteria {
                                if individual.level[b_idx] > 0.001 {
                                    let base = b_idx * num_drugs;
                                    // Count if infected with this bacteria and on any drug
                                    if individual.cur_use_drug.iter().any(|&x| x) {
                                        lt.infected_and_on_any_drug_by_bacteria[b_idx] += 1;
                                    }
                                    for d_idx in 0..num_drugs {
                                        let resistance_data = &individual.resistances[b_idx][d_idx];
                                        let threshold = mic_lt2_thresholds[base + d_idx];
                                        if resistance_data.majority_r < threshold { lt.mic_lt2_counts[base + d_idx] += 1; }
                                        // per-bacteria, per-drug currently on drug
                                        if individual.cur_use_drug[d_idx] {
                                            lt.currently_on_drug_by_bacteria_drug[base + d_idx] += 1;
                                        }
                                        // Sum any_r values for infected individuals
                                        lt.any_r_sum_by_bacteria_drug[base + d_idx] += resistance_data.any_r;
                                        // Calculate and sum MIC values for infected individuals
                                        // MIC = 1 / ((1 - majority_r) * potency)
                                        let potency = self.potency_matrix[base + d_idx];
                                        let mic = 1.0 / ((1.0 - resistance_data.majority_r) * potency);
                                        lt.mic_sum_by_bacteria_drug[base + d_idx] += mic;
                                        // Count individuals with any_r > 0 for infected individuals
                                        if resistance_data.any_r > 0.0 {
                                            lt.infected_with_any_r_positive_by_bacteria_drug[base + d_idx] += 1;
                                        }
                                        // Sum any_r values for hospital-acquired infected individuals
                                        if individual.infection_hospital_acquired[b_idx] {
                                            lt.any_r_sum_by_bacteria_drug_hospital[base + d_idx] += resistance_data.any_r;
                                        }
                                        // Sum any_r values by region (pooled across all bacteria and drugs)
                                        let effective_region = get_effective_region(individual);
                                        let region_idx = region_to_index(effective_region);
                                        lt.any_r_sum_by_region[region_idx] += resistance_data.any_r;
                                    }
                                    
                                    // Count resistance mechanisms for this bacteria
                                    let num_mechanisms = ResistanceMechanism::all().len();
                                    for (mech_idx, _mechanism) in ResistanceMechanism::all().iter().enumerate() {
                                        if individual.resistance_mechanisms[b_idx][mech_idx] {
                                            let flat_idx = b_idx * num_mechanisms + mech_idx;
                                            lt.infected_with_bacteria_and_mechanism[flat_idx] += 1;
                                        }
                                    }
                                    
                                    // Count infected people with test status flags
                                    if individual.test_identified_infection[b_idx] {
                                        lt.infected_with_test_identified_by_bacteria[b_idx] += 1;
                                    }
                                    if individual.test_for_resistance[b_idx] {
                                        lt.infected_with_test_for_resistance_by_bacteria[b_idx] += 1;
                                    }
                                }                                // Count microbiome_r > 0 for all bacteria-drug combinations (regardless of infection status)
                                for d_idx in 0..num_drugs {
                                    let resistance_data = &individual.resistances[b_idx][d_idx];
                                    if resistance_data.microbiome_r > 0.0 {
                                        let idx = b_idx * num_drugs + d_idx;
                                        lt.microbiome_r_positive_by_bacteria_drug[idx] += 1;
                                    }
                                }
                            }
                        }

                        // Reset drug scores for this time step (initialize to -1 indicating no drug selection)
                        individual.bacteria_on_selection_day = -1;
                        for d_idx in 0..num_drugs {
                            individual.drug_score_on_selection_day[d_idx] = -1.0;
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
                                
                                // Get region for this death
                                let effective_region = get_effective_region(individual);
                                let region_idx = region_to_index(effective_region);
                                
                                // Get age group for this death (ages in days, convert to years)
                                let age_years = individual.age as f64 / 365.0;
                                let age_group_idx = if (0.0..6.0).contains(&age_years) { 
                                    0 // 0-5 years
                                } else if (6.0..15.0).contains(&age_years) { 
                                    1 // 6-14 years
                                } else if (15.0..50.0).contains(&age_years) { 
                                    2 // 15-49 years
                                } else if (50.0..80.0).contains(&age_years) { 
                                    3 // 50-79 years
                                } else { 
                                    4 // 80+ years
                                };
                                
                                if let Some(ref cause) = individual.cause_of_death {
                                    match cause.as_str() {
                                        "background_mortality" => {
                                            lt.deaths_background += 1;
                                            lt.deaths_by_region[region_idx * 3 + 0] += 1; // background death
                                            lt.deaths_by_region_age[region_idx * 15 + age_group_idx * 3 + 0] += 1; // background death by age
                                        },
                                        "sepsis_related" => {
                                            lt.deaths_sepsis += 1;
                                            lt.deaths_by_region[region_idx * 3 + 1] += 1; // sepsis death
                                            lt.deaths_by_region_age[region_idx * 15 + age_group_idx * 3 + 1] += 1; // sepsis death by age
                                            
                                            // Track sepsis deaths by syndrome and region
                                            for syndrome_idx in 0..10 { // syndromes 1-10 -> indices 0-9
                                                if individual.sepsis[syndrome_idx] {
                                                    let index = syndrome_idx * 6 + region_idx;
                                                    lt.syndrome_deaths_sepsis_by_region[index] += 1;
                                                }
                                            }
                                        },
                                        "drug_toxicity_related" => {
                                            lt.deaths_drug_toxicity += 1;
                                            lt.deaths_by_region[region_idx * 3 + 2] += 1; // toxicity death
                                            lt.deaths_by_region_age[region_idx * 15 + age_group_idx * 3 + 2] += 1; // toxicity death by age
                                        },
                                        _ => {
                                            lt.deaths_background += 1;
                                            lt.deaths_by_region[region_idx * 3 + 0] += 1; // default to background
                                            lt.deaths_by_region_age[region_idx * 15 + age_group_idx * 3 + 0] += 1; // default to background by age
                                        },
                                    }
                                } else { 
                                    lt.deaths_background += 1; 
                                    lt.deaths_by_region[region_idx * 3 + 0] += 1; // default to background
                                    lt.deaths_by_region_age[region_idx * 15 + age_group_idx * 3 + 0] += 1; // default to background by age
                                }
                                // Count deaths by bacteria
                                for b_idx in 0..num_bacteria {
                                    if individual.level[b_idx] > 0.001 {
                                        lt.deaths_by_bacteria[b_idx] += 1;
                                    }
                                }
                                
                                // Count deaths by bacteria and home region for currently infected individuals
                                let home_region_idx = region_to_index(individual.region_living);
                                for b_idx in 0..num_bacteria {
                                    if individual.level[b_idx] > 0.001 {
                                        lt.deaths_infected_by_bacteria_region[b_idx * 6 + home_region_idx] += 1;
                                    }
                                }
                            }
                        }

                        if individual.date_of_death.is_none() && individual.age >= 0 {
                            // Integrated living population & age groups (only count individuals who have been born)
                            lt.living_population += 1;
                            
                            // Count living population by region
                            let effective_region = get_effective_region(individual);
                            let region_idx = region_to_index(effective_region);
                            lt.living_population_by_region[region_idx] += 1;
                            
                            let age_years = individual.age as f64 / 365.0;
                            if (0.0..6.0).contains(&age_years) { 
                                lt.num_age_0_5 += 1; 
                                lt.age_distribution_by_region[region_idx * 5 + 0] += 1;
                            }
                            else if (6.0..15.0).contains(&age_years) { 
                                lt.num_age_6_14 += 1; 
                                lt.age_distribution_by_region[region_idx * 5 + 1] += 1;
                            }
                            else if (15.0..50.0).contains(&age_years) { 
                                lt.num_age_15_49 += 1; 
                                lt.age_distribution_by_region[region_idx * 5 + 2] += 1;
                            }
                            else if (50.0..80.0).contains(&age_years) { 
                                lt.num_age_50_79 += 1; 
                                lt.age_distribution_by_region[region_idx * 5 + 3] += 1;
                            }
                            else if age_years >= 80.0 { 
                                lt.num_age_80plus += 1; 
                                lt.age_distribution_by_region[region_idx * 5 + 4] += 1;
                            }
                            // Drug usage post-rules
                            let mut on_any_drug = false;
                            for (d_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
                                if is_using { 
                                    lt.currently_on_drug_by_drug[d_idx] += 1; 
                                    // Count drug usage by region
                                    let idx = region_idx * DRUG_SHORT_NAMES.len() + d_idx;
                                    lt.currently_on_drug_by_region_drug[idx] += 1;
                                    on_any_drug = true; 
                                }
                            }
                            if on_any_drug { lt.currently_taking_drug_count += 1; }

                            // Check if this person started any drug today
                            let mut started_drug_today = false;
                            for &initiation_date in individual.date_drug_initiated.iter() {
                                if initiation_date == t as i32 {
                                    started_drug_today = true;
                                    break;
                                }
                            }
                            if started_drug_today {
                                lt.new_drug_initiations_count += 1;
                                
                                // Check if this person is currently infected (any level > 0)
                                let is_currently_infected = individual.level.iter().any(|&level| level > 0.0);
                                if is_currently_infected {
                                    lt.new_drug_initiations_count_infected += 1;
                                }
                            }

                            if individual.presence_microbiome.iter().any(|&x| x) { lt.num_with_any_bacteria_microbiome += 1; }
                            
                            // Count presence_microbiome by individual bacteria
                            for (b_idx, &has_bacteria) in individual.presence_microbiome.iter().enumerate() {
                                if has_bacteria {
                                    lt.presence_microbiome_by_bacteria[b_idx] += 1;
                                    // Also count by region
                                    let region_idx = individual.region_living as usize;
                                    lt.presence_microbiome_by_bacteria_by_region[b_idx][region_idx] += 1;
                                }
                            }

                            // Track drug failure events: check for day 5 post-drug-initiation
                            let region_idx = individual.region_living as usize;
                            for (d_idx, &drug_init_day) in individual.date_drug_initiated.iter().enumerate() {
                                if drug_init_day != i32::MIN && t as i32 - drug_init_day == 5 {
                                    // This is day 5 post-drug-initiation for this drug
                                    // Check all bacteria for infection status
                                    for b_idx in 0..individual.level.len() {
                                        lt.drug_treatment_day5_events_by_bacteria_region[b_idx][region_idx] += 1;
                                        
                                        // Check if this is a failure: still on drug AND still infected
                                        if individual.cur_use_drug[d_idx] && individual.level[b_idx] > 0.0 {
                                            lt.drug_failure_events_by_bacteria_region[b_idx][region_idx] += 1;
                                        }
                                    }
                                }
                            }

                            // Infection & resistance
                            let mut individual_max_infection_duration = 0;
                            let mut individual_has_any_infection = false;
                            let mut individual_has_any_r_positive = false;
                            let mut was_newly_infected = false;
                            let mut was_newly_infected_with_resistance = false;
                            let mut individual_has_any_infection_counted_for_syndrome = false;
                            let is_currently_infected_any;
                            {
                                let mut infected_any_tmp = false;
                                for b_idx in 0..num_bacteria {
                                    if individual.level[b_idx] > 0.001 {
                                        infected_any_tmp = true;
                                        individual_has_any_infection = true;
                                        lt.infections_by_bacteria[b_idx] += 1;
                                        
                                        // Count syndrome for this infected individual (take first one if multiple infections)
                                        if !individual_has_any_infection_counted_for_syndrome {
                                            let syndrome_id = individual.infectious_syndrome[b_idx];
                                            if syndrome_id >= 1 && syndrome_id <= 10 {
                                                lt.infected_by_syndrome[(syndrome_id - 1) as usize] += 1;
                                                individual_has_any_infection_counted_for_syndrome = true;
                                            }
                                        }
                                        
                                        // Count syndrome for this bacteria specifically (all infections, not just first)
                                        let syndrome_id = individual.infectious_syndrome[b_idx];
                                        if syndrome_id >= 1 && syndrome_id <= 10 {
                                            let flat_idx = b_idx * 10 + (syndrome_id - 1) as usize;
                                            lt.infected_by_syndrome_by_bacteria[flat_idx] += 1;
                                        }
                                        
                                        // sum activity_r for this bacteria, ONLY for individuals on drug
                                        let mut activity_r_sum = 0.0;
                                        let days_since_infection = t as i32 - individual.date_last_infected[b_idx];
                                        if days_since_infection > individual_max_infection_duration { individual_max_infection_duration = days_since_infection; }
                                        if individual.date_last_infected[b_idx] == t as i32 { 
                                            was_newly_infected = true; 
                                            // Count new active infections by bacteria and home region
                                            let home_region_idx = region_to_index(individual.region_living);
                                            let flat_idx = b_idx * 6 + home_region_idx;
                                            lt.newly_infected_by_bacteria_region[flat_idx] += 1;
                                        }
                                        let base = b_idx * num_drugs;
                                        for d_idx in 0..num_drugs {
                                            let resistance_data = &individual.resistances[b_idx][d_idx];
                                            // Only sum activity_r if individual is currently on this drug
                                            if individual.cur_use_drug[d_idx] {
                                                activity_r_sum += resistance_data.activity_r;
                                            }
                                            if resistance_data.majority_r > 0.0 { lt.resistance_by_bacteria_drug[base + d_idx] += 1; lt.majority_r_entries.push(((individual.region_cur_in as usize, individual.hospital_status.is_hospitalized(), b_idx, d_idx), resistance_data.majority_r)); }
                                            if resistance_data.any_r > 0.0 {
                                                individual_has_any_r_positive = true;
                                                if individual.date_last_infected[b_idx] == t as i32 && !was_newly_infected_with_resistance {
                                                    lt.newly_infected_with_resistance_count += 1;
                                                    was_newly_infected_with_resistance = true;
                                                }
                                                // Count newly acquired resistance by acquisition type per bacteria-drug combination
                                                if let Some(acq_type) = individual.how_resistance_acquired[b_idx][d_idx] {
                                                    use crate::simulation::population::ResistanceAcquisitionType;
                                                    let index = b_idx * num_drugs + d_idx;
                                                    match acq_type {
                                                        ResistanceAcquisitionType::AtInfectionCommunity => lt.new_resistance_at_infection_community_by_bacteria_drug[index] += 1,
                                                        ResistanceAcquisitionType::AtInfectionEnv => lt.new_resistance_at_infection_env_by_bacteria_drug[index] += 1,
                                                        ResistanceAcquisitionType::Hgt => lt.new_resistance_hgt_by_bacteria_drug[index] += 1,
                                                        ResistanceAcquisitionType::FromMicrobiomeR => lt.new_resistance_from_microbiome_r_by_bacteria_drug[index] += 1,
                                                    }
                                                }
                                            }
                                        }
                                        // Only include individuals who are on any drug for this bacteria
                                        if individual.cur_use_drug.iter().any(|&x| x) {
                                            lt.activity_r_sum_by_bacteria[b_idx] += activity_r_sum;
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
                            if individual.immunodeficiency_type.is_some() { lt.number_severely_immunosuppressed += 1; }
                            if individual.sepsis.iter().any(|&s| s) { lt.number_with_sepsis += 1; }
                            
                            // Track sepsis by bacteria and new sepsis cases
                            for b_idx in 0..num_bacteria {
                                if individual.sepsis[b_idx] {
                                    // Current sepsis with this bacteria
                                    lt.number_with_sepsis_by_bacteria[b_idx] += 1;
                                    
                                    // Count as new sepsis case if sepsis started today and person is currently infected
                                    if individual.level[b_idx] > 0.001 && individual.sepsis_onset_day[b_idx] == t as i32 {
                                        lt.new_sepsis_cases_by_bacteria[b_idx] += 1;
                                    }
                                }
                            }
                        }
                        lt
                    })
                    .reduce(|| LocalTotals::new(num_bacteria, num_drugs, per_thread_cap), |mut a, b| { a.merge(b); a });

                // Collect infection resolution data after rules have been applied
                let infection_resolution_totals = self.population.individuals.par_iter()
                    .fold(|| (
                        vec![0usize; num_bacteria], // immune_clearance
                        vec![0usize; num_bacteria], // drug_assisted_clearance
                        vec![0usize; num_bacteria], // death_from_sepsis
                        vec![0usize; num_bacteria], // death_from_background
                        vec![0usize; num_bacteria], // death_from_toxicity
                    ), |mut acc, individual| {
                        for (b_idx, resolution_counts) in individual.infection_resolution_this_timestep.iter().enumerate() {
                            acc.0[b_idx] += resolution_counts[0] as usize;
                            acc.1[b_idx] += resolution_counts[1] as usize;
                            acc.2[b_idx] += resolution_counts[2] as usize;
                            acc.3[b_idx] += resolution_counts[3] as usize;
                            acc.4[b_idx] += resolution_counts[4] as usize;
                        }
                        acc
                    })
                    .reduce(|| (
                        vec![0usize; num_bacteria],
                        vec![0usize; num_bacteria],
                        vec![0usize; num_bacteria],
                        vec![0usize; num_bacteria],
                        vec![0usize; num_bacteria],
                    ), |mut a, b| {
                        for i in 0..num_bacteria {
                            a.0[i] += b.0[i];
                            a.1[i] += b.1[i];
                            a.2[i] += b.2[i];
                            a.3[i] += b.3[i];
                            a.4[i] += b.4[i];
                        }
                        a
                    });

                // Destructure to move out (avoid cloning large vectors)
                let LocalTotals {
                    infected_and_on_any_drug_by_bacteria,
                    mic_lt2_counts: infected_and_standardized_mic_lt2_by_bacteria_drug,
                    currently_on_drug_by_bacteria_drug,
                    microbiome_r_positive_by_bacteria_drug,
                    infections_by_bacteria: infections_by_bacteria_vec,
                    deaths_by_bacteria,
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
                    number_with_sepsis_by_bacteria,
                    new_sepsis_cases_by_bacteria,
                    newly_infected_count,
                    newly_infected_with_resistance_count,
                    new_drug_initiations_count,
                    new_drug_initiations_count_infected,
                    newly_infected_by_bacteria_region,
                    deaths_infected_by_bacteria_region,
                    total_currently_infected,
                    total_with_resistance,
                    currently_infected_and_on_drug_count,
                    num_with_any_bacteria_microbiome,
                    presence_microbiome_by_bacteria,
                    presence_microbiome_by_bacteria_by_region,
                    drug_failure_events_by_bacteria_region,
                    drug_treatment_day5_events_by_bacteria_region,
                    infected_with_test_identified_by_bacteria,
                    infected_with_test_for_resistance_by_bacteria,
                    living_population,
                    num_age_0_5,
                    num_age_6_14,
                    num_age_15_49,
                    num_age_50_79,
                    num_age_80plus,
                    activity_r_sum_by_bacteria,
                    any_r_sum_by_bacteria_drug,
                    any_r_sum_by_bacteria_drug_hospital,
                    infected_with_any_r_positive_by_bacteria_drug,
                    mic_sum_by_bacteria_drug,
                    any_r_sum_by_region,
                    infected_count_by_region,
                    infected_with_bacteria_and_mechanism,
                    new_resistance_at_infection_community_by_bacteria_drug,
                    new_resistance_at_infection_env_by_bacteria_drug,
                    new_resistance_hgt_by_bacteria_drug,
                    new_resistance_from_microbiome_r_by_bacteria_drug,
                    infection_resolution_immune_clearance_by_bacteria: _,
                    infection_resolution_drug_assisted_clearance_by_bacteria: _,
                    infection_resolution_death_from_sepsis_by_bacteria: _,
                    infection_resolution_death_from_background_by_bacteria: _,
                    infection_resolution_death_from_toxicity_by_bacteria: _,
                    infected_by_syndrome,
                    infected_by_syndrome_by_bacteria,
                    living_population_by_region,
                    age_distribution_by_region,
                    deaths_by_region,
                    deaths_by_region_age,
                    currently_on_drug_by_region_drug,
                    syndrome_deaths_sepsis_by_region,
                } = totals;

                // Use the separately collected infection resolution data
                let (
                    infection_resolution_immune_clearance_by_bacteria,
                    infection_resolution_drug_assisted_clearance_by_bacteria,
                    infection_resolution_death_from_sepsis_by_bacteria,
                    infection_resolution_death_from_background_by_bacteria,
                    infection_resolution_death_from_toxicity_by_bacteria,
                ) = infection_resolution_totals;

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
                infected_and_on_any_drug_by_bacteria,
                infected_and_standardized_mic_lt2_by_bacteria_drug,
                currently_on_drug_by_bacteria_drug,
                microbiome_r_positive_by_bacteria_drug,
                any_r_sum_by_bacteria_drug,
                any_r_sum_by_bacteria_drug_hospital,
                infected_with_any_r_positive_by_bacteria_drug,
                mic_sum_by_bacteria_drug,
                any_r_sum_by_region,
                infected_count_by_region,
                currently_on_drug_by_drug,
                num_age_0_5,
                num_age_6_14,
                num_age_15_49,
                num_age_50_79,
                num_age_80plus,
                num_with_any_bacteria_microbiome,
                presence_microbiome_by_bacteria,
                presence_microbiome_by_bacteria_by_region,
                drug_failure_events_by_bacteria_region,
                drug_treatment_day5_events_by_bacteria_region,
                infected_with_test_identified_by_bacteria,
                infected_with_test_for_resistance_by_bacteria,
                time_step: t,
                total_population: living_population,
                number_in_hospital,
                number_severely_immunosuppressed,
                number_with_sepsis,
                number_with_sepsis_by_bacteria,
                new_sepsis_cases_by_bacteria,
                newly_infected_count,
                newly_infected_with_resistance_count,
                new_drug_initiations_count,
                new_drug_initiations_count_infected,
                newly_infected_by_bacteria_region,
                deaths_infected_by_bacteria_region,
                total_currently_infected,
                total_with_resistance,
                infected_10_days_count: infected_10_count,
                infected_30_days_count: infected_30_count,
                currently_taking_drug_count,
                taking_two_drugs_count,
                infections_by_bacteria: infections_by_bacteria_vec,
                deaths_by_bacteria,
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
        activity_r_sum_by_bacteria,
        infected_with_bacteria_and_mechanism,
        new_resistance_at_infection_community_by_bacteria_drug,
        new_resistance_at_infection_env_by_bacteria_drug,
        new_resistance_hgt_by_bacteria_drug,
        new_resistance_from_microbiome_r_by_bacteria_drug,
        infection_resolution_immune_clearance_by_bacteria,
        infection_resolution_drug_assisted_clearance_by_bacteria,
        infection_resolution_death_from_sepsis_by_bacteria,
        infection_resolution_death_from_background_by_bacteria,
        infection_resolution_death_from_toxicity_by_bacteria,
        
        // Calculate day-7 drug initiation statistics
        day_7_evaluations_by_bacteria: {
            let evaluation_days = get_global_param("drug_evaluation_days_post_infection").unwrap_or(7.0) as i32;
            let mut day_7_evals = vec![0; BACTERIA_LIST.len()];
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                for b_idx in 0..BACTERIA_LIST.len() {
                    let infection_start_day = individual.date_last_infected_keep[b_idx];
                    
                    // Only count if today is exactly the evaluation day after infection start (i.e., evaluation happens TODAY)
                    if infection_start_day > 0 && (t as i32) == (infection_start_day + evaluation_days) {
                        day_7_evals[b_idx] += 1;
                    }
                }
            }
            day_7_evals
        },
        day_7_drug_used_by_bacteria: {
            let evaluation_days = get_global_param("drug_evaluation_days_post_infection").unwrap_or(7.0) as i32;
            let mut day_7_used = vec![0; BACTERIA_LIST.len()];
            
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                for b_idx in 0..BACTERIA_LIST.len() {
                    let infection_start_day = individual.date_last_infected_keep[b_idx];
                    
                    // Only count if today is exactly the evaluation day after infection start AND drug was used
                    if infection_start_day > 0 && (t as i32) == (infection_start_day + evaluation_days) {
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
                        
                        if drug_used_since_infection {
                            day_7_used[b_idx] += 1;
                        }
                    }
                }
            }
            
            day_7_used
        },
        infected_by_syndrome,
        infected_by_syndrome_by_bacteria,
        living_population_by_region,
        hospital_population_by_region: {
            let mut hospital_pop_by_region = vec![0; 6]; // 6 regions
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                if individual.hospital_status.is_hospitalized() {
                    let region_idx = get_effective_region(individual) as usize;
                    hospital_pop_by_region[region_idx] += 1;
                }
            }
            hospital_pop_by_region
        },
        newly_infected_hospital_by_bacteria_region: {
            let mut hospital_infections = HashMap::new();
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                if individual.hospital_status.is_hospitalized() {
                    let region_idx = get_effective_region(individual) as usize;
                    
                    for b_idx in 0..BACTERIA_LIST.len() {
                        if individual.date_last_infected_keep[b_idx] == t as i32 {
                            // This is a new infection that occurred today in hospital
                            *hospital_infections.entry((b_idx, region_idx)).or_insert(0) += 1;
                        }
                    }
                }
            }
            hospital_infections
        },
        age_distribution_by_region,
        deaths_by_region,
        deaths_by_region_age,
        syndrome_population_by_region: {
            let mut syndrome_pop_by_region = vec![0; 60]; // 10 syndromes * 6 regions
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                let region_idx = get_effective_region(individual) as usize;
                
                // Count individuals with active infections by syndrome
                for syndrome_idx in 0..10 { // syndromes 1-10 -> indices 0-9
                    if individual.sepsis[syndrome_idx] {
                        let index = syndrome_idx * 6 + region_idx;
                        syndrome_pop_by_region[index] += 1;
                    }
                }
            }
            syndrome_pop_by_region
        },
        syndrome_deaths_sepsis_by_region: {
            syndrome_deaths_sepsis_by_region
        },
        currently_on_drug_by_region_drug,
        
        // Calculate polypharmacy distribution (1, 2, or ≥3 drugs)
        people_on_1_drug: {
            let mut count = 0;
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                if individual.current_number_of_drugs == 1 {
                    count += 1;
                }
            }
            count
        },
        people_on_2_drugs: {
            let mut count = 0;
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                if individual.current_number_of_drugs == 2 {
                    count += 1;
                }
            }
            count
        },
        people_on_3plus_drugs: {
            let mut count = 0;
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                if individual.current_number_of_drugs >= 3 {
                    count += 1;
                }
            }
            count
        },
        
        // Calculate infected people on drug with previous treatment failure
        infected_on_drug_with_previous_failure: {
            let mut count = 0;
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                // Check if person is currently infected
                let currently_infected = individual.level.iter().any(|&level| level > 0.0);
                if !currently_infected { continue; }
                
                // Check if person is currently on any drug
                let on_any_drug = individual.cur_use_drug.iter().any(|&is_on| is_on);
                if !on_any_drug { continue; }
                
                // Check if person has had treatment failure assessed (has previous failure experience)
                let has_previous_failure = individual.treatment_failure_assessed.iter().any(|&assessed| assessed);
                if has_previous_failure {
                    count += 1;
                }
            }
            count
        },
        
        // Drug score tracking for clinical guideline debugging
        drug_selection_count_by_bacteria: {
            let mut counts = vec![0; BACTERIA_LIST.len()];
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                // Count if drug selection occurred for this individual today (bacteria_on_selection_day >= 0)
                if individual.bacteria_on_selection_day >= 0 && (individual.bacteria_on_selection_day as usize) < BACTERIA_LIST.len() {
                    counts[individual.bacteria_on_selection_day as usize] += 1;
                }
            }
            counts
        },
        
        drug_score_sums_by_bacteria_drug: {
            let mut sums = vec![0.0; BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len()];
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                // Add drug scores if drug selection occurred today
                if individual.bacteria_on_selection_day >= 0 && (individual.bacteria_on_selection_day as usize) < BACTERIA_LIST.len() {
                    let bacteria_idx = individual.bacteria_on_selection_day as usize;
                    
                    for (drug_idx, &score) in individual.drug_score_on_selection_day.iter().enumerate() {
                        if drug_idx < DRUG_SHORT_NAMES.len() && score >= 0.0 { // Valid score
                            let flat_idx = bacteria_idx * DRUG_SHORT_NAMES.len() + drug_idx;
                            sums[flat_idx] += score;
                        }
                    }
                }
            }
            sums
        },
        
        people_by_drug_count: {
            let mut drug_count_histogram = vec![0; 4]; // 0, 1, 2, 3+ drugs
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; } // Skip dead individuals
                
                let drug_count = individual.current_number_of_drugs as usize;
                let histogram_index = if drug_count >= 3 { 3 } else { drug_count }; // Cap at 3+ drugs
                drug_count_histogram[histogram_index] += 1;
            }
            drug_count_histogram
        }
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

            // Reset infection resolution counts for next timestep (after data has been aggregated and logged)
            self.population.individuals.par_iter_mut().for_each(|individual| {
                for b_idx in 0..BACTERIA_LIST.len() {
                    for res_idx in 0..crate::simulation::population::InfectionResolutionType::all().len() {
                        individual.infection_resolution_this_timestep[b_idx][res_idx] = 0;
                    }
                }
            });

            if self.log_individuals {
                use std::fs::OpenOptions;
                let n_log = 10.min(self.population.individuals.len());
                let log_path = "individuals_log.csv";
                let is_first_timestep = t == 0;
                let mut file = if is_first_timestep {
                    // Overwrite file on first timestep
                    OpenOptions::new().create(true).write(true).truncate(true).open(log_path).expect("Unable to open individuals_log.csv")
                } else {
                    // Append to file for subsequent timesteps
                    OpenOptions::new().create(true).append(true).open(log_path).expect("Unable to open individuals_log.csv")
                };
                // Write header only on first timestep
                if is_first_timestep {
                    writeln!(file, "time_step,individual_index,id,age,sex_at_birth,region_living,region_cur_in,current_infection_related_death_risk,background_all_cause_mortality_rate,current_toxicity,mortality_risk_current_toxicity,hospital_status,is_severely_immunosuppressed,date_of_death,level,immune_resp,presence_microbiome,cur_level_drug,cur_use_drug,ever_taken_drug,date_last_infected,cur_infection_from_environment,infection_hospital_acquired,test_identified_infection,sepsis,infection_resolution_this_timestep,active_infection_activity_r,day_7_since_last_infection_drug_used,resistances_microbiome_r,resistances_test_r,resistances_activity_r,resistances_any_r,resistances_majority_r,resistance_mechanisms,bacteria_on_selection_day,drug_score_on_selection_day,date_last_drug_failure,current_number_of_drugs").unwrap();
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
                    // Flatten infection resolution data for all bacteria/resolution types
                    let mut infection_resolutions = Vec::new();
                    for bact_resolutions in &ind.infection_resolution_this_timestep {
                        for &count in bact_resolutions {
                            infection_resolutions.push(count);
                        }
                    }
                    
                    // Calculate active infection activity_r value
                    let active_infection_activity_r = {
                        // Find first bacteria where person is infected and on drug
                        let mut result = 0.0;
                        for b_idx in 0..BACTERIA_LIST.len() {
                            if ind.level[b_idx] > 0.0 && ind.cur_use_drug.iter().any(|&on_drug| on_drug) {
                                // Person is infected with this bacteria and on some drug
                                // Use activity_r from first drug for this bacteria
                                for d_idx in 0..DRUG_SHORT_NAMES.len() {
                                    if ind.cur_use_drug[d_idx] {
                                        result = ind.resistances[b_idx][d_idx].activity_r;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        result
                    };

                    // Format day_7_since_last_infection_drug_used field
                    let fmt_day_7_drug_used = ind.day_7_since_last_infection_drug_used.iter()
                        .map(|opt| match opt {
                            Some(true) => "true",
                            Some(false) => "false", 
                            None => "null"
                        })
                        .collect::<Vec<_>>()
                        .join(";");

                    writeln!(file, "{},{},{},{},{},{:?},{:?},{:.4},{:.4},{:.4},{:.4},{:?},{},{:?},{},{},{},{},{},{},{},{},{},{},{},{:.4},{},{},{},{},{},{},{},{},{},{},{},{}",
                        t,
                        i,
                        ind.id,
                        ind.age,
                        ind.sex_at_birth,
                        ind.region_living,
                        ind.region_cur_in,
                        ind.current_infection_related_death_risk,
                        ind.background_all_cause_mortality_rate,
                        ind.current_toxicity,
                        ind.mortality_risk_current_toxicity,
                        format!("{:?}", ind.hospital_status),
                        format!("{:?}", ind.immunodeficiency_type),
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
                        fmt_vec(&infection_resolutions),
                        active_infection_activity_r,
                        fmt_day_7_drug_used,
                        fmt_vec(&microbiome_r),
                        fmt_vec(&test_r),
                        fmt_vec(&activity_r),
                        fmt_vec(&any_r),
                        fmt_vec(&majority_r),
                        mechanisms.join(";"),
                        ind.bacteria_on_selection_day,
                        fmt_vec(&ind.drug_score_on_selection_day),
                        fmt_vec(&ind.date_last_drug_failure),
                        ind.current_number_of_drugs
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
        use std::io::{Write, BufWriter};

        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        // Pre-build header string once
        let mut header = String::with_capacity(50000); // Pre-allocate large capacity
        header.push_str("time_step,total_population,number_in_hospital,number_severely_immunosuppressed,number_with_sepsis,total_currently_infected,infected_10_days_count,infected_30_days_count,total_with_resistance,currently_taking_drug_count,currently_infected_and_on_drug_count,taking_two_drugs_count,newly_infected_count,newly_infected_with_resistance_count,new_drug_initiations_count,new_drug_initiations_count_infected,newly_infected_past_year,total_deaths,deaths_background,deaths_sepsis,deaths_drug_toxicity,deaths_past_year,deaths_background_past_year,deaths_sepsis_past_year,deaths_drug_toxicity_past_year,num_age_0_5,num_age_6_14,num_age_15_49,num_age_50_79,num_age_80plus,num_with_any_bacteria_microbiome,people_on_1_drug,people_on_2_drugs,people_on_3plus_drugs,infected_on_drug_with_previous_failure");
        
        // Add per-bacteria infection columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_currently_infected");
        }
        // Add per-bacteria deaths columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_deaths");
        }
        // Add per-bacteria sepsis prevalence columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_number_with_sepsis");
        }
        // Add per-bacteria sepsis incidence columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_new_sepsis_cases");
        }
        // Add per-bacteria activity_r sum columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_activity_r_sum");
        }
        // Add per-bacteria presence_microbiome columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_presence_microbiome");
        }
        // Add per-bacteria per-region presence_microbiome columns
        for bacteria in BACTERIA_LIST.iter() {
            for region in &["north_america", "south_america", "africa", "asia", "europe", "oceania"] {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_presence_microbiome_");
                header.push_str(region);
            }
        }
        // Add per-bacteria per-region drug failure events columns
        for bacteria in BACTERIA_LIST.iter() {
            for region in &["north_america", "south_america", "africa", "asia", "europe", "oceania"] {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_drug_failure_events_");
                header.push_str(region);
            }
        }
        // Add per-bacteria per-region drug treatment day5 events columns  
        for bacteria in BACTERIA_LIST.iter() {
            for region in &["north_america", "south_america", "africa", "asia", "europe", "oceania"] {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_drug_treatment_day5_events_");
                header.push_str(region);
            }
        }
        // Add per-bacteria infected with test_identified_infection columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infected_with_test_identified");
        }
        // Add per-bacteria infected with test_for_resistance columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infected_with_test_for_resistance");
        }
        
        // Add per-bacteria, per-region newly infected columns
        let regions = ["north_america", "south_america", "africa", "asia", "europe", "oceania"];
        for bacteria in BACTERIA_LIST.iter() {
            for region in &regions {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_newly_infected_");
                header.push_str(region);
            }
        }
        
        // Add per-bacteria, per-region deaths (currently infected) columns
        for bacteria in BACTERIA_LIST.iter() {
            for region in &regions {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_deaths_infected_");
                header.push_str(region);
            }
        }
        // Add per-drug currently on drug columns
        for drug in DRUG_SHORT_NAMES.iter() {
            header.push(',');
            header.push_str(&drug.replace(" ", "_"));
            header.push_str("_currently_on_drug");
        }
        // Add per-bacteria, per-drug MIC < 2 columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infected_and_mic_lt2_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug currently on drug columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_currently_on_drug_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug microbiome_r > 0 columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_microbiome_r_positive_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug any_r sum columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_sum_any_r_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug infected with any_r > 0 count columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infected_with_any_r_positive_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug MIC sum columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_sum_mic_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug any_r sum columns for hospital-acquired infections
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_sum_any_r_hospital_");
                header.push_str(drug);
            }
        }
        // Add per-region any_r sum columns (pooled across all bacteria and drugs)
        let region_names = ["north_america", "south_america", "africa", "asia", "europe", "oceania"];
        for region in region_names.iter() {
            header.push(',');
            header.push_str(region);
            header.push_str("_any_r_sum");
        }
        // Add per-region infected count columns
        for region in region_names.iter() {
            header.push(',');
            header.push_str(region);
            header.push_str("_infected_count");
        }
        // Add per-region, per-drug currently on drug columns
        for region in region_names.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(region);
                header.push('_');
                header.push_str(drug);
                header.push_str("_currently_on_drug");
            }
        }
        // Add per-bacteria infected and on any drug columns to header (after other per-bacteria columns)
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infected_and_on_any_drug");
        }
        // Add per-bacteria, per-resistance-mechanism columns to header
        for bacteria in BACTERIA_LIST.iter() {
            for mechanism in ResistanceMechanism::all() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infected_with_");
                header.push_str(mechanism.as_str());
            }
        }
        // Add per-bacteria, per-drug resistance acquisition columns to header
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push('_');
                header.push_str(drug);
                header.push_str("_new_resistance_at_infection_community");
            }
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push('_');
                header.push_str(drug);
                header.push_str("_new_resistance_at_infection_env");
            }
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push('_');
                header.push_str(drug);
                header.push_str("_new_resistance_hgt");
            }
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push('_');
                header.push_str(drug);
                header.push_str("_new_resistance_from_microbiome_r");
            }
        }
        
        // Add per-bacteria infection resolution columns to header
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_immune_clearance");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_drug_assisted_clearance");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_death_from_sepsis");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_death_from_background");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_death_from_toxicity");
        }
        
        // Add per-bacteria day-7 drug initiation columns to header
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_day_7_evaluations");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_day_7_drug_used");
        }
        
        // Add syndrome columns to header  
        for syndrome_id in 1..=10 {
            header.push(',');
            header.push_str(&format!("syndrome_{}_infected", syndrome_id));
        }
        
        // Add bacteria-specific syndrome columns to header
        for bacteria in BACTERIA_LIST.iter() {
            for syndrome_id in 1..=10 {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str(&format!("_syndrome_{}_infected", syndrome_id));
            }
        }
        
        // Add region population columns to header
        let region_names = ["north_america", "south_america", "africa", "asia", "europe", "oceania"];
        for region_name in &region_names {
            header.push(',');
            header.push_str(&format!("{}_population", region_name));
        }
        
        // Add regional hospital population columns to header
        for region_name in &region_names {
            header.push(',');
            header.push_str(&format!("{}_hospital_population", region_name));
        }
        
        // Add per-bacteria, per-region hospital newly infected columns to header
        for bacteria in BACTERIA_LIST.iter() {
            for region in &region_names {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_newly_infected_hospital_");
                header.push_str(region);
            }
        }
        
        // Add regional age distribution columns to header
        let age_group_names = ["prop_age_0_5", "prop_age_6_14", "prop_age_15_49", "prop_age_50_79", "prop_age_80plus"];
        for region_name in &region_names {
            for age_group_name in &age_group_names {
                header.push(',');
                header.push_str(&format!("{}_{}", region_name, age_group_name));
            }
        }
        
        // Add regional death columns to header
        let death_type_names = ["deaths_background", "deaths_sepsis", "deaths_drug_toxicity"];
        for region_name in &region_names {
            for death_type_name in &death_type_names {
                header.push(',');
                header.push_str(&format!("{}_{}", region_name, death_type_name));
            }
        }
        
        // Add age-specific death columns to header (region x age_group x death_type)
        for region_name in &region_names {
            for age_group_name in &age_group_names {
                for death_type_name in &death_type_names {
                    header.push(',');
                    header.push_str(&format!("{}_{}_{}", region_name, age_group_name, death_type_name));
                }
            }
        }
        
        // Add syndrome population by region columns to header
        for syndrome_id in 1..=10 { // syndromes 1-10
            for region_name in &region_names {
                header.push(',');
                header.push_str(&format!("syndrome_{}_population_{}", syndrome_id, region_name));
            }
        }
        
        // Add syndrome deaths from sepsis by region columns to header
        for syndrome_id in 1..=10 { // syndromes 1-10
            for region_name in &region_names {
                header.push(',');
                header.push_str(&format!("syndrome_{}_deaths_sepsis_{}", syndrome_id, region_name));
            }
        }
        
        // Add regional drug usage columns to header (region x drug)
        for region_name in &region_names {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&format!("{}_currently_on_drug_{}", region_name, drug.replace(" ", "_")));
            }
        }
        
        // Add drug score tracking columns to header
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_drug_selection_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_drug_score_sum_");
                header.push_str(drug);
            }
        }
        
        // Add drug count histogram columns
        header.push_str(",people_on_0_drugs,people_on_1_drugs_new,people_on_2_drugs_new,people_on_3plus_drugs_new");
        
        header.push('\n');
        writer.write_all(header.as_bytes())?;

        // Write data with pre-built strings
        for summary in &self.summary_log {
            let mut row = String::with_capacity(20000); // Pre-allocate for each row
            
            // Write basic summary data
            row.push_str(&format!("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}", 
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
                summary.new_drug_initiations_count,
                summary.new_drug_initiations_count_infected,
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
                summary.people_on_1_drug,
                summary.people_on_2_drugs,
                summary.people_on_3plus_drugs,
                summary.infected_on_drug_with_previous_failure,
            ));
            
            // Remove the duplicate polypharmacy data that was causing mismatch
            // (these values are now included in the main format string above)
            
            // Append all array data efficiently
            for value in &summary.infections_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.deaths_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.number_with_sepsis_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.new_sepsis_cases_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.activity_r_sum_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.presence_microbiome_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            // Add regional presence_microbiome data
            for bacteria_vec in &summary.presence_microbiome_by_bacteria_by_region {
                for value in bacteria_vec { row.push(','); row.push_str(&value.to_string()); }
            }
            // Add regional drug failure events data
            for bacteria_vec in &summary.drug_failure_events_by_bacteria_region {
                for value in bacteria_vec { row.push(','); row.push_str(&value.to_string()); }
            }
            // Add regional drug treatment day5 events data
            for bacteria_vec in &summary.drug_treatment_day5_events_by_bacteria_region {
                for value in bacteria_vec { row.push(','); row.push_str(&value.to_string()); }
            }
            for value in &summary.infected_with_test_identified_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.infected_with_test_for_resistance_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.newly_infected_by_bacteria_region { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.deaths_infected_by_bacteria_region { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.currently_on_drug_by_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.infected_and_standardized_mic_lt2_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.currently_on_drug_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.microbiome_r_positive_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.any_r_sum_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.any_r_sum_by_bacteria_drug_hospital { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.infected_with_any_r_positive_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            

            
            for value in &summary.mic_sum_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.any_r_sum_by_region { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.infected_count_by_region { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.currently_on_drug_by_region_drug { row.push(','); row.push_str(&value.to_string()); }
            

            
            for value in &summary.infected_and_on_any_drug_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.infected_with_bacteria_and_mechanism { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.new_resistance_at_infection_community_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.new_resistance_at_infection_env_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.new_resistance_hgt_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.new_resistance_from_microbiome_r_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            
            for value in &summary.infection_resolution_immune_clearance_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            
            for value in &summary.infection_resolution_drug_assisted_clearance_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.infection_resolution_death_from_sepsis_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.infection_resolution_death_from_background_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.infection_resolution_death_from_toxicity_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            
            // Add day-7 drug initiation data
            for value in &summary.day_7_evaluations_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.day_7_drug_used_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            
            // Add syndrome infection data
            for value in &summary.infected_by_syndrome { row.push(','); row.push_str(&value.to_string()); }
            
            // Add bacteria-specific syndrome infection data
            for value in &summary.infected_by_syndrome_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            
            // Add region population data
            for value in &summary.living_population_by_region { row.push(','); row.push_str(&value.to_string()); }
            
            // Add regional hospital population data
            for value in &summary.hospital_population_by_region { row.push(','); row.push_str(&value.to_string()); }
            
            // Add per-bacteria, per-region hospital newly infected data
            for bacteria_idx in 0..BACTERIA_LIST.len() {
                for region_idx in 0..6 { // 6 regions
                    let count = summary.newly_infected_hospital_by_bacteria_region.get(&(bacteria_idx, region_idx)).unwrap_or(&0);
                    row.push(',');
                    row.push_str(&count.to_string());
                }
            }
            
            // Add regional age distribution data (as proportions)
            for region_idx in 0..6 { // 6 regions
                let region_pop = summary.living_population_by_region[region_idx];
                for age_group_idx in 0..5 { // 5 age groups
                    let age_count = summary.age_distribution_by_region[region_idx * 5 + age_group_idx];
                    let proportion = if region_pop > 0 { age_count as f64 / region_pop as f64 } else { 0.0 };
                    row.push(',');
                    row.push_str(&format!("{:.6}", proportion));
                }
            }
            
            // Add regional death data (as counts)
            for region_idx in 0..6 { // 6 regions
                for death_type_idx in 0..3 { // 3 death types: background, sepsis, drug_toxicity
                    let death_count = summary.deaths_by_region[region_idx * 3 + death_type_idx];
                    row.push(',');
                    row.push_str(&death_count.to_string());
                }
            }
            
            // Add age-specific death data by region (as counts)
            for region_idx in 0..6 { // 6 regions
                for age_group_idx in 0..5 { // 5 age groups
                    for death_type_idx in 0..3 { // 3 death types: background, sepsis, drug_toxicity
                        let death_count = summary.deaths_by_region_age[region_idx * 15 + age_group_idx * 3 + death_type_idx];
                        row.push(',');
                        row.push_str(&death_count.to_string());
                    }
                }
            }
            
            // Add syndrome population by region data
            for syndrome_idx in 0..10 { // syndromes 1-10 -> indices 0-9
                for region_idx in 0..6 { // 6 regions
                    let population_count = summary.syndrome_population_by_region[syndrome_idx * 6 + region_idx];
                    row.push(',');
                    row.push_str(&population_count.to_string());
                }
            }
            
            // Add syndrome deaths from sepsis by region data
            for syndrome_idx in 0..10 { // syndromes 1-10 -> indices 0-9
                for region_idx in 0..6 { // 6 regions
                    let death_count = summary.syndrome_deaths_sepsis_by_region[syndrome_idx * 6 + region_idx];
                    row.push(',');
                    row.push_str(&death_count.to_string());
                }
            }
            
            // Add drug score tracking data
            for value in &summary.drug_selection_count_by_bacteria { row.push(','); row.push_str(&value.to_string()); }
            for value in &summary.drug_score_sums_by_bacteria_drug { row.push(','); row.push_str(&value.to_string()); }
            
            // Add drug count histogram data
            for count in &summary.people_by_drug_count {
                row.push(',');
                row.push_str(&count.to_string());
            }
            
            row.push('\n');
            
            writer.write_all(row.as_bytes())?;
        }

        writer.flush()?;
        println!("Summary data exported to {}", filename);
        Ok(())
    }
}