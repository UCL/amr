// src/simulation/simulation.rs
use crate::simulation::population::{Population, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::rules::apply_rules;
use crate::config; // Import the config module
use std::collections::HashMap;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

// Compact structure for time step summary data
#[allow(dead_code)]
#[derive(Clone)]
pub struct TimeStepSummary {
    pub time_step: usize,
    pub total_population: usize,
    pub total_deaths: usize,
    pub deaths_background: usize,        // Deaths from background mortality
    pub deaths_sepsis: usize,           // Deaths from sepsis
    pub deaths_drug_toxicity: usize,    // Deaths from drug toxicity
    pub total_with_resistance: usize,
    pub total_currently_infected: usize, // Number of living people currently infected with any bacteria
    pub currently_taking_drug_count: usize, // New field
    pub infected_10_days_count: usize,     // New field
    pub infected_30_days_count: usize,     // New field
    pub taking_two_drugs_count: usize,     // New field
    pub number_in_hospital: usize,         // New field
    pub number_severely_immunosuppressed: usize, // New field
    pub number_with_sepsis: usize,         // New field
    pub infections_by_bacteria: Vec<usize>, // indexed by bacteria
    pub resistance_by_bacteria_drug: Vec<Vec<usize>>, // [bacteria][drug] counts
    pub newly_infected_count: usize, // Number of people newly infected this time step
}

pub struct Simulation {  // public rust struct which encapsulates the state and configuration of a simulation run.
    pub population: Population, // specifying the population of individuals in the simulation.
    pub time_steps: usize, // specifying how many discrete time steps the simulation will run.

    // todo: ensure that when we count across individuals that we include only those alive

    pub global_majority_r_proportions: HashMap<(usize, usize), f64>,  // Maps (bacteria_index, drug_index) pairs to a global proportion 
                                                                      // value to track summary statistics over time.
    pub bacteria_indices: HashMap<&'static str, usize>, // A string-to-index map converting bacteria names (&'static str) to integer indices.
    pub drug_indices: HashMap<&'static str, usize>, // as above, but for drugs.
    pub cross_resistance_groups: HashMap<usize, Vec<Vec<usize>>>, // New: (b_idx -> [[d_idx, d_idx], ...])
    pub current_majority_r_positive_values_by_combo: HashMap<(usize, bool, usize, usize), Vec<f64>>, // Store between time steps

    pub summary_log: Vec<TimeStepSummary>, // Efficient storage for summary data
}

impl Simulation {
    pub fn new(population_size: usize, time_steps: usize) -> Self {

        // public function named new (rust’s conventional constructor pattern).  
        // Takes two inputs: population_size: how many individuals to initialize.
        // time_steps: how many time steps the simulation should run.
        // Returns Self → shorthand for returning an instance of Simulation.

        let population = Population::new(population_size); 

        // calls a new constructor for the Population struct.  Passes in "population_size", returning a Population instance 
        // and stores it in the local population variable.

        // Initialize bacteria_indices and drug_indices
        let mut bacteria_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &bacteria) in BACTERIA_LIST.iter().enumerate() { // Iterate over the bacteria list and create a mapping from bacteria names to their indices.
            bacteria_indices.insert(bacteria, i); // Inserts each bacteria name and its index into the HashMap.
        }

        let mut drug_indices: HashMap<&'static str, usize> = HashMap::new(); // Create a HashMap to map drug names to their indices.
        for (i, &drug) in DRUG_SHORT_NAMES.iter().enumerate() { // Iterate over the drug list and create a mapping from drug names to their indices.
            drug_indices.insert(drug, i);
        }

        // New: Load and process cross-resistance groups
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

        let global_majority_r_proportions = HashMap::new(); // Initialize an empty HashMap to store global majority_r proportions for bacteria/drug pairs.

        // --- Initial State Logging for Individual 0

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
            global_majority_r_proportions,
            bacteria_indices,
            drug_indices,
            cross_resistance_groups, // Add new field
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
//          println!("simulation.rs time step: {}", t);

            // Use previous time step's resistance data for new acquisitions
            let previous_majority_r_positive_values_by_combo = if t == 0 {
                HashMap::new() // Empty for first time step
            } else {
                // Use the data collected in the previous iteration
                std::mem::take(&mut self.current_majority_r_positive_values_by_combo)
            };

            // Initialize counters and data structures for this time step
            let mut new_majority_r_positive_values_by_combo: HashMap<(usize, bool, usize, usize), Vec<f64>> = HashMap::new();
            let mut current_infected_counts_with_majority_r: HashMap<(usize, usize), usize> = HashMap::new();
            let mut current_infected_counts_total: HashMap<usize, usize> = HashMap::new();
            let mut _individuals_with_any_bacterial_infection = 0;

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

            // Single pass: apply rules and collect all statistics
            self.population.individuals.par_iter_mut().for_each(|individual| {
                // Apply rules first
                apply_rules(
                    individual,
                    t,
                    &self.global_majority_r_proportions,
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
                    // Count individuals currently taking a drug
                    if individual.cur_use_drug.iter().any(|&is_using| is_using) {
                        currently_taking_drug_count.fetch_add(1, Ordering::Relaxed);
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
            for individual in self.population.individuals.iter() {
                let region_idx = individual.region_cur_in as usize;
                let hospital_status_bool = individual.hospital_status.is_hospitalized();

                // Only count if individual is alive
                let is_alive = individual.date_of_death.is_none();

                if is_alive {
                    for (b_idx, &_bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                        if individual.level[b_idx] > 0.001 {
                            *current_infected_counts_total.entry(b_idx).or_insert(0) += 1;

                            for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                                let resistance_data = &individual.resistances[b_idx][d_idx];
                                if resistance_data.majority_r > 0.0 {
                                    new_majority_r_positive_values_by_combo
                                        .entry((region_idx, hospital_status_bool, b_idx, d_idx))
                                        .or_insert_with(Vec::new)
                                        .push(resistance_data.majority_r);
                                    *current_infected_counts_with_majority_r.entry((b_idx, d_idx)).or_insert(0) += 1;
                                }
                            }
                        }
                    }
                }
            }

            // Store for next iteration
            self.current_majority_r_positive_values_by_combo = new_majority_r_positive_values_by_combo;

            // Count living population (age >= 0 and no date_of_death)
            let living_population = self.population.individuals.iter()
                .filter(|individual| individual.age >= 0 && individual.date_of_death.is_none())
                .count();

            // Create summary for this time step
            let infected_10_count = infected_10_days_count.load(Ordering::Relaxed);
            let infected_30_count = infected_30_days_count.load(Ordering::Relaxed);
            
            // Debug: Always print the values for time steps where we saw issues
            if t >= 15 && t <= 25 {
                println!("DEBUG FINAL: Time step {}: infected_10_count={}, infected_30_count={}", 
                       t, infected_10_count, infected_30_count);
            }
            
            // Validation: 30-day count should never exceed 10-day count
            if infected_30_count > infected_10_count {
                println!("*** CRITICAL ERROR DETECTED ***");
                println!("ERROR at time step {}: infected_30_days_count ({}) > infected_10_days_count ({})", 
                         t, infected_30_count, infected_10_count);
                println!("This is logically impossible and indicates a bug in the counting logic.");
                println!("*** END ERROR REPORT ***");
            }
            
            let summary = TimeStepSummary {
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
            };
            self.summary_log.push(summary);


// /*  per time step printing block

            // --- print activity_r for all infected bacteria/drug pairs for individual 0 after update ---
            let individual_0 = &self.population.individuals[0];
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() { 
                if individual_0.level[b_idx] > 0.0001 {
                    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                        if individual_0.cur_level_drug[drug_idx] > 0.0 {
                            let resistance_data = &individual_0.resistances[b_idx][drug_idx];
                            println!("   "); 
                            println!(
                                "simulation.rs  {} (infected) + {} (present): activity_r = {:.4}, any_r = {:.4}, drug_level = {:.4}",
                                bacteria_name,
                                drug_name,
                                resistance_data.activity_r,
                                resistance_data.any_r,
                                individual_0.cur_level_drug[drug_idx]
                            );
                            println!("   "); 
                        }
                    }
                }
            }


            // Print drug details for individual 0, regardless of infection status
            // Note: Using threshold of 0.001 to avoid showing negligible drug levels
            let mut drugs_present_found_overall = false; // Declare and initialize here
            for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                if individual_0.cur_level_drug[drug_idx] > 0.001 {
                    let status = if individual_0.cur_use_drug[drug_idx] {
                        " simulation.rs (currently being taken)"
                    } else {
                        " simulation.rs (decaying)"
                    };
                    println!("simulation.rs ");
                    println!("{}: level = {:.4}{}", drug_name_static, individual_0.cur_level_drug[drug_idx], status);
                    println!(" ");
                    drugs_present_found_overall = true; // Use the newly declared variable
                }
            }
            if !drugs_present_found_overall {
                println!("simulation.rs  no antibiotics currently in system");
            }


            let mut has_infection = false;
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                let level = individual_0.level[b_idx];
                if level > 0.0001 {
                    has_infection = true;
                    println!(" ");  
                    println!("simulation.rs  ");  
                    println!(" ");  
                    println!("bacteria level = {:.4}", level);
                    println!("bacteria: {}", bacteria_name);
                    println!("infected = true");

                    println!("immune response = {:.4}", individual_0.immune_resp[b_idx]);
                    println!("infection from environment = {}", individual_0.cur_infection_from_environment[b_idx]);
                    println!("hospital acquired infection = {}", individual_0.infection_hospital_acquired[b_idx]);
                    println!("test identified infection = {}", individual_0.test_identified_infection[b_idx]);
                    println!("date_last_infected = {}", individual_0.date_last_infected[b_idx]);
                    let mut drugs_present_found = false;
                    println!("antibiotics present in system (current level > 0):");
                    for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                        if individual_0.cur_level_drug[drug_idx] > 0.0 {
                            let status = if individual_0.cur_use_drug[drug_idx] {
                                " (currently being taken)"
                            } else {
                                " (decaying)"
                            };
                            println!("{}: level = {:.4}{}", drug_name_static, individual_0.cur_level_drug[drug_idx], status);
                            drugs_present_found = true;
                        }
                    }
                    if !drugs_present_found {
                        println!("simulation.rs  no antibiotics currently in system");
                    }
                    let mut effective_antibiotics_found = false;
  
                    for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                        if individual_0.cur_level_drug[drug_idx] > 0.0 {
                            let resistance_data = &individual_0.resistances[b_idx][drug_idx];
                            println!("any_r {}:", bacteria_name);    
                            println!(
                                "simulation.rs  {}: level = {:.4}, any_r = {:.4}, activity_r = {:.4}, majority_r = {:.4}",
                                drug_name_static,
                                individual_0.cur_level_drug[drug_idx],
                                resistance_data.any_r,
                                resistance_data.activity_r,
                                resistance_data.majority_r
                            );
                            if resistance_data.activity_r > 0.0 {
                                effective_antibiotics_found = true;
                            }
                        }
                    }
                    if !effective_antibiotics_found {
                        println!("simulation.rs  no effective antibiotics in system against this bacteria");
                    }
                    println!();
                }
            }
            if !has_infection {
                println!("simulation.rs  no active bacterial infection as of end of the time step");
                println!();
            }


            println!(" ");
            println!("simulation.rs  infection and resistance summary outputs:");
            println!(" ");

            let age_in_years = (self.population.individuals[0].age as f64 / 365.0).round() as i32;
            let ever_taken_drug_vector: Vec<u8> = self.population.individuals[0].ever_taken_drug.iter().map(|&taken| if taken { 1 } else { 0 }).collect();
            println!("                                ");
            println!("age_in_years: {}", age_in_years);
            println!("region_living: {:?}", self.population.individuals[0].region_living);                                      
            println!("region_cur_in: {:?}", self.population.individuals[0].region_cur_in);                                      
            println!("hospital_status: {:?}", self.population.individuals[0].hospital_status);                                      
            println!("is_severely_immunosuppressed: {:?}", self.population.individuals[0].is_severely_immunosuppressed);                                      
            println!("date_last_infected: {:?}", self.population.individuals[0].date_last_infected);                                      
            println!("ever_taken_drug: {:?}", ever_taken_drug_vector);
            println!("date of death: {:?}", self.population.individuals[0].date_of_death);   
            println!("                                ");

            // Print resistance summary for all infected individuals at this time step
//          self.print_resistance_summary(t);


//  */  //  end of per timestep printing block


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
        writeln!(file, "time_step,total_population,number_in_hospital,number_severely_immunosuppressed,number_with_sepsis,total_currently_infected,infected_10_days_count,infected_30_days_count,total_with_resistance,currently_taking_drug_count,taking_two_drugs_count,newly_infected_count,total_deaths,deaths_background,deaths_sepsis,deaths_drug_toxicity")?;

        // Write data
        for summary in &self.summary_log {
            writeln!(file, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}", 
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
                summary.taking_two_drugs_count,
                summary.newly_infected_count,
                summary.total_deaths,
                summary.deaths_background,
                summary.deaths_sepsis,
                summary.deaths_drug_toxicity,
            )?;
        }
        
        println!("Summary data exported to {}", filename);
        Ok(())
    }
}