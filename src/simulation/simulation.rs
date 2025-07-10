// src/simulation/simulation.rs
use crate::simulation::population::{Population, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::rules::apply_rules;
use crate::config; // Import the config module
use std::collections::HashMap;
use rayon::prelude::*;

pub struct Simulation {  // public rust struct which encapsulates the state and configuration of a simulation run.
    pub population: Population, // specifying the population of individuals in the simulation.
    pub time_steps: usize, // specifying how many discrete time steps the simulation will run.

    // todo: ensure that when we count across individuals that we include only those alive

    pub global_majority_r_proportions: HashMap<(usize, usize), f64>,  // Maps (bacteria_index, drug_index) pairs to a global proportion 
                                                                      // value to track summary statistics over time.
    pub bacteria_indices: HashMap<&'static str, usize>, // A string-to-index map converting bacteria names (&'static str) to integer indices.
    pub drug_indices: HashMap<&'static str, usize>, // as above, but for drugs.
    pub cross_resistance_groups: HashMap<usize, Vec<Vec<usize>>>, // New: (b_idx -> [[d_idx, d_idx], ...])
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
        }
    }

    pub fn run(&mut self) {
        // public function named run, which executes the simulation for the specified number of time steps.

        println!(" ");
        println!("--- starting to run over time steps");
        println!(" ");

        for t in 0..self.time_steps {
            println!("simulation.rs time step: {}", t);


            let mut current_majority_r_positive_values_by_combo: HashMap<(usize, bool, usize, usize), Vec<f64>> = HashMap::new();
            let mut current_infected_counts_with_majority_r: HashMap<(usize, usize), usize> = HashMap::new();
            let mut current_infected_counts_total: HashMap<usize, usize> = HashMap::new();

            // --- counters  ---
            let mut _individuals_with_any_bacterial_infection = 0;
            let mut _individuals_with_any_r_positive_for_any_bacteria = 0;
            // --- ---

            for individual in self.population.individuals.iter() {
                let region_idx = individual.region_cur_in as usize;
                let hospital_status_bool = individual.hospital_status.is_hospitalized();
                let mut individual_has_any_infection = false;
                let mut individual_has_any_r_positive = false;
                for (b_idx, &_bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                    if individual.level[b_idx] > 0.001 {
                        individual_has_any_infection = true;
                        *current_infected_counts_total.entry(b_idx).or_insert(0) += 1;
                        for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                            let resistance_data = &individual.resistances[b_idx][d_idx];
                            if resistance_data.majority_r > 0.0 {
                                current_majority_r_positive_values_by_combo
                                    .entry((region_idx, hospital_status_bool, b_idx, d_idx))
                                    .or_insert_with(Vec::new)
                                    .push(resistance_data.majority_r);
                                *current_infected_counts_with_majority_r.entry((b_idx, d_idx)).or_insert(0) += 1;
                            }
                            if resistance_data.any_r > 0.0 {
                                individual_has_any_r_positive = true;
                            }
                        }
                    }
                }
                if individual_has_any_infection {
                    _individuals_with_any_bacterial_infection += 1;
                }
                if individual_has_any_r_positive {
                    _individuals_with_any_r_positive_for_any_bacteria += 1;
                }
            }

            // --- parallel application of rules to individuals ---
            self.population.individuals.par_iter_mut().for_each(|individual| {
                apply_rules(
                    individual,
                    t,
                    &self.global_majority_r_proportions,
                    &current_majority_r_positive_values_by_combo,
                    &self.bacteria_indices,
                    &self.drug_indices,
                    &self.cross_resistance_groups, // Pass new data
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
        }

    }


/*   


    fn print_resistance_summary(&self, time_step: usize) {
        // Calculate bacteria infection counts first
        let mut bacteria_infection_counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
        
        for individual in &self.population.individuals {
            for (bacteria, &b_idx) in self.bacteria_indices.iter() {
                if individual.level[b_idx] > 0.001 {
                    *bacteria_infection_counts.entry(bacteria).or_insert(0) += 1;
                }
            }
        }

        
        // Print bacteria and resistance summary
        println!("\n--- Time step {} - Bacteria infection and resistance summary ---", time_step);
        for (bacteria, &count) in &bacteria_infection_counts {
            println!("{}: {} infected", bacteria, count);   
            for (drug, _) in self.drug_indices.iter() {
                // Collect the full distribution of any_r for this bacteria/drug pair
                let mut any_r_values = Vec::new();
                for individual in &self.population.individuals {
                    if let Some(&b_idx) = self.bacteria_indices.get(bacteria) {
                        if individual.level[b_idx] > 0.001 {
                            if let Some(&d_idx) = self.drug_indices.get(drug) {
                                let any_r = individual.resistances[b_idx][d_idx].any_r;
                                any_r_values.push(any_r);
                            }
                        }
                    }
                }

                // Print summary statistics for the distribution
                if !any_r_values.is_empty() {
                    let n = any_r_values.len() as f64;
                    let mut count_0 = 0;
                    let mut count_001_025 = 0;
                    let mut count_0251_05 = 0;
                    let mut count_0501_075 = 0;
                    let mut count_0751_1 = 0;
                    for &val in &any_r_values {
                        if val == 0.0 {
                            count_0 += 1;
                        } else if val > 0.0 && val <= 0.25 {
                            count_001_025 += 1;
                        } else if val > 0.25 && val <= 0.5 {
                            count_0251_05 += 1;
                        } else if val > 0.5 && val <= 0.75 {
                            count_0501_075 += 1;
                        } else if val > 0.75 && val <= 1.0 {
                            count_0751_1 += 1;
                        }
                    }
                    println!(
                        "    {}: n = {}, prop 0.00 = {:.3}, prop 0.25 = {:.3}, prop 0.5 = {:.3}, prop 0.75 = {:.3}, prop 1.00 = {:.3}",
                        drug,
                        n as usize,
                        count_0 as f64 / n,
                        count_001_025 as f64 / n,
                        count_0251_05 as f64 / n,
                        count_0501_075 as f64 / n,
                        count_0751_1 as f64 / n
                    );
                } else {
                    println!("    {}: n = 0", drug);
                }
            }
        }
    }
*/

}

