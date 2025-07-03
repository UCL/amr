// src/simulation/simulation.rs
use crate::simulation::population::{Population, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::rules::apply_rules;
use std::collections::HashMap;
use rayon::prelude::*;

pub struct Simulation {  // public rust struct which encapsulates the state and configuration of a simulation run.
    pub population: Population, // specifying the population of individuals in the simulation.
    pub time_steps: usize, // specifying how many discrete time steps the simulation will run.

    // todo: ensure that when we count across individuals that we include only those alive

    pub global_majority_r_proportions: HashMap<(usize, usize), f64>,  // Maps (bacteria_index, drug_index) pairs to a global proportion 
                                                                      // value to track summary statistics over time.
    pub majority_r_positive_values_by_combo: HashMap<(usize, usize), Vec<f64>>,
    pub bacteria_indices: HashMap<&'static str, usize>, // A string-to-index map converting bacteria names (&'static str) to integer indices.
    pub drug_indices: HashMap<&'static str, usize>, // as above, but for drugs.
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

        let global_majority_r_proportions = HashMap::new(); // Initialize an empty HashMap to store global majority_r proportions for bacteria/drug pairs.
        let majority_r_positive_values_by_combo = HashMap::new(); // Initialize an empty HashMap to store majority_r positive values for each bacteria/drug combination.

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
            majority_r_positive_values_by_combo,
            bacteria_indices,
            drug_indices,
        }
    }

    pub fn run(&mut self) { // public function named run, which executes the simulation for the specified number of time steps.

        println!(" ");
        println!("--- starting to run over time steps");
        println!(" ");

        for t in 0..self.time_steps {

            println!("simulation.rs time step: {}", t);

            // --- parallel application of rules to individuals ---
            // each individual's rules are applied independently.
            self.population.individuals.par_iter_mut().for_each(|individual| { // this uses Rayon to parallelize the application of rules across all individuals.
                apply_rules( // Calls the apply_rules function, passing in the individual and other necessary data structures.
                    individual,
                    t,
                    &self.global_majority_r_proportions,
                    &self.majority_r_positive_values_by_combo,
                    &self.bacteria_indices,
                    &self.drug_indices,
                );
            });

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

            // --- sequential aggregation of global statistics ---
            // This part must be sequential because it collects data from all individuals
            // into shared, mutable HashMaps, which would cause data races if done in parallel
            // without complex synchronization.
            let mut current_majority_r_positive_values_by_combo: HashMap<(usize, usize), Vec<f64>> = HashMap::new();
            let mut current_infected_counts_with_majority_r: HashMap<(usize, usize), usize> = HashMap::new();
            let mut current_infected_counts_total: HashMap<usize, usize> = HashMap::new();

            // --- counters  ---
            let mut individuals_with_any_bacterial_infection = 0;
            let mut individuals_with_any_r_positive_for_any_bacteria = 0;
            // --- ---


            for individual in self.population.individuals.iter() { // Iterate over each individual in the population to collect statistics.
                // --- Flags for this individual's status ---
                let mut individual_has_any_infection = false;
                let mut individual_has_any_r_positive = false;
                // --- END ---

                for (b_idx, &_bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                    // Only count if currently infected
                    if individual.level[b_idx] > 0.001 {
                        // -flag if individual has an infection ---
                        individual_has_any_infection = true;
                        // --- end -

                        *current_infected_counts_total.entry(b_idx).or_insert(0) += 1; // Increment the count of total infections for this bacteria index

                        for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() { // Iterate over each drug for the current bacteria
                            let resistance_data = &individual.resistances[b_idx][d_idx]; // get the resistance data for the current bacteria/drug pair
                            if resistance_data.majority_r > 0.0 {
                                current_majority_r_positive_values_by_combo //
                                    .entry((b_idx, d_idx))
                                    .or_insert_with(Vec::new)
                                    .push(resistance_data.majority_r);
                                *current_infected_counts_with_majority_r.entry((b_idx, d_idx)).or_insert(0) += 1; // Increment the count of infections with majority_r for this bacteria/drug pair
                            }
                            // check for any_r > 0 for ANY bacteria/drug combo for this individual ---
                            if resistance_data.any_r > 0.0 {
                                individual_has_any_r_positive = true;
                            }
                            // --- end ---
                        }
                    }
                }
                // ---Increment overall counters for the individual AFTER checking all their infections ---
                if individual_has_any_infection {
                    individuals_with_any_bacterial_infection += 1;
                }
                if individual_has_any_r_positive {
                    individuals_with_any_r_positive_for_any_bacteria += 1;
                }
                // --- end ---
            }

 
            // Print drug details for individual 0, regardless of infection status

            let mut drugs_present_found_overall = false; // Declare and initialize here
            for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                if individual_0.cur_level_drug[drug_idx] > 0.0 {
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

            // --- Global Infection and Resistance Statistics Output ---
            let total_population_size = self.population.individuals.len();
            let proportion_any_r_positive = if individuals_with_any_bacterial_infection > 0 {
                individuals_with_any_r_positive_for_any_bacteria as f64 / individuals_with_any_bacterial_infection as f64
            } else {
                0.0
            };

            println!(" ");
            println!("simulation.rs  infection and resistance summary outputs:");
            println!(" ");
            println!("total individuals in population: {}", total_population_size);
            println!("number of individuals with any bacterial infection: {}", individuals_with_any_bacterial_infection);
            println!("number of individuals with any bacteria having any_r > 0: {}", individuals_with_any_r_positive_for_any_bacteria);
            println!("proportion of infected individuals with any_r > 0: {:.4}\n", proportion_any_r_positive);
            // --- end  ---

            println!("                                ");
            
        }


    }
}