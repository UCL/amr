// src/simulation/simulation.rs
use crate::simulation::population::{Population, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::rules::apply_rules;
use std::collections::HashMap;
use rayon::prelude::*;

pub struct Simulation {
    pub population: Population,
    pub time_steps: usize,
    pub global_majority_r_proportions: HashMap<(usize, usize), f64>,
    pub majority_r_positive_values_by_combo: HashMap<(usize, usize), Vec<f64>>,
    pub bacteria_indices: HashMap<&'static str, usize>,
    pub drug_indices: HashMap<&'static str, usize>,
}

impl Simulation {
    pub fn new(population_size: usize, time_steps: usize) -> Self {
        println!("--- AMR SIMULATION ---");
        let population = Population::new(population_size);

        // Initialize bacteria_indices and drug_indices
        let mut bacteria_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &bacteria) in BACTERIA_LIST.iter().enumerate() {
            bacteria_indices.insert(bacteria, i);
        }

        let mut drug_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            drug_indices.insert(drug, i);
        }

        let global_majority_r_proportions = HashMap::new();
        let majority_r_positive_values_by_combo = HashMap::new();

        // --- Initial State Logging for Individual 0
        println!("Initial age of individual 0 AFTER population creation: {} days", population.individuals[0].age);
        println!("--- INITIAL STATE OF INDIVIDUAL 0 (from main.rs) ---");
        println!("  Region Living: {:?}", population.individuals[0].region_living);
        println!("  Region Currently In: {:?}", population.individuals[0].region_cur_in);
        let strep_pneu_idx = bacteria_indices["strep_pneu"];
        let amoxicillin_idx = drug_indices["amoxicillin"];
        let kleb_pneu_idx = bacteria_indices["kleb_pneu"];
        let ceftriaxone_idx = drug_indices["ceftriaxone"];

        println!("  strep_pneu: level = {:.2}", population.individuals[0].level[strep_pneu_idx]);
        println!("  strep_pneu: immune_resp = {:.2}", population.individuals[0].immune_resp[strep_pneu_idx]);
        println!("  strep_pneu: sepsis = {}", population.individuals[0].sepsis[strep_pneu_idx]);
        println!("  strep_pneu: infectious_syndrome = {}", population.individuals[0].infectious_syndrome[strep_pneu_idx]);
        println!("  strep_pneu: date_last_infected = {}", population.individuals[0].date_last_infected[strep_pneu_idx]);
        for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
            println!("  {}_vaccination_status: {}", bacteria_name, population.individuals[0].vaccination_status[b_idx]);
        }
        println!("  cur_use_amoxicillin: {}", population.individuals[0].cur_use_drug[amoxicillin_idx]);
        println!("  cur_level_amoxicillin: {:.2}", population.individuals[0].cur_level_drug[amoxicillin_idx]);
        println!("  current_infection_related_death_risk: {:.2}", population.individuals[0].current_infection_related_death_risk);
        println!("  background_all_cause_mortality_rate: {:.4}", population.individuals[0].background_all_cause_mortality_rate);
        println!("  sexual_contact_level: {:.2}", population.individuals[0].sexual_contact_level);
        println!("  airborne_contact_level_with_adults: {:.2}", population.individuals[0].airborne_contact_level_with_adults);
        println!("  airborne_contact_level_with_children: {:.2}", population.individuals[0].airborne_contact_level_with_children);
        println!("  oral_exposure_level: {:.2}", population.individuals[0].oral_exposure_level);
        println!("  mosquito_exposure_level: {:.2}", population.individuals[0].mosquito_exposure_level);
        println!("  current_toxicity: {:.2}", population.individuals[0].current_toxicity);
        println!("  mortality_risk_current_toxicity: {:.2}", population.individuals[0].mortality_risk_current_toxicity);
        println!("  strep_pneu: infection_hospital_acquired = {}", population.individuals[0].infection_hospital_acquired[strep_pneu_idx]);
        println!("  strep_pneu: cur_infection_from_environment = {}", population.individuals[0].cur_infection_from_environment[strep_pneu_idx]);
        let resistance_data = &population.individuals[0].resistances[strep_pneu_idx][amoxicillin_idx];
        println!("  strep_pneu resistance to amoxicillin:");
        println!("    microbiome_r: {:.2}", resistance_data.microbiome_r);
        println!("    test_r: {:.2}", resistance_data.test_r);
        println!("    activity_r: {:.2}", resistance_data.activity_r);
        println!("    any_r: {:.2}", resistance_data.any_r);
        println!("    majority_r: {:.2}", resistance_data.majority_r);
        println!("  kleb_pneu: infection_hospital_acquired = {}", population.individuals[0].infection_hospital_acquired[kleb_pneu_idx]);
        println!("  kleb_pneu: cur_infection_from_environment = {}", population.individuals[0].cur_infection_from_environment[kleb_pneu_idx]);
        let resistance_data = &population.individuals[0].resistances[kleb_pneu_idx][ceftriaxone_idx];
        println!("  kleb_pneu resistance to ceftriaxone:");
        println!("    microbiome_r: {:.2}", resistance_data.microbiome_r);
        println!("    test_r: {:.2}", resistance_data.test_r);
        println!("    activity_r: {:.2}", resistance_data.activity_r);
        println!("    any_r: {:.2}", resistance_data.any_r);
        println!("    majority_r: {:.2}", resistance_data.majority_r);


        Simulation {
            population,
            time_steps,
            global_majority_r_proportions,
            majority_r_positive_values_by_combo,
            bacteria_indices,
            drug_indices,
        }
    }

    pub fn run(&mut self) {
        println!("--- SIMULATION STARTING ---");
        println!("--- AMR SIMULATION (within run function) ---");
        println!("Initial age of individual 0 BEFORE simulation: {} days", self.population.individuals[0].age);
        println!("--- SIMULATION STARTING (within run function) ---");

        // Temporary data collection for global majority_r proportions
        let mut strep_pneu_amox_majority_r_history: Vec<f64> = Vec::new();
        let mut kleb_pneu_ceftr_majority_r_history: Vec<f64> = Vec::new();

        for t in 0..self.time_steps {
            let step_start = std::time::Instant::now();

            // --- PARALLEL APPLICATION OF RULES TO INDIVIDUALS ---
            // Each individual's rules are applied independently.
            self.population.individuals.par_iter_mut().for_each(|individual| {
                apply_rules(
                    individual,
                    t,
                    &self.global_majority_r_proportions,
                    &self.majority_r_positive_values_by_combo,
                    &self.bacteria_indices,
                    &self.drug_indices,
                );
            });

            // --- Print activity_r for all infected bacteria/drug pairs for Individual 0 after update ---
            let individual_0 = &self.population.individuals[0];
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                if individual_0.level[b_idx] > 0.0001 {
                    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                        if individual_0.cur_level_drug[drug_idx] > 0.0 {
                            let resistance_data = &individual_0.resistances[b_idx][drug_idx];
                            println!(
                                "After time step {}: {} (infected) + {} (present): activity_r = {:.4}, any_r = {:.4}, drug_level = {:.4}",
                                t,
                                bacteria_name,
                                drug_name,
                                resistance_data.activity_r,
                                resistance_data.any_r,
                                individual_0.cur_level_drug[drug_idx]
                            );
                        }
                    }
                }
            }

            // --- SEQUENTIAL AGGREGATION OF GLOBAL STATISTICS ---
            // This part must be sequential because it collects data from all individuals
            // into shared, mutable HashMaps, which would cause data races if done in parallel
            // without complex synchronization.
            let mut current_majority_r_positive_values_by_combo: HashMap<(usize, usize), Vec<f64>> = HashMap::new();
            let mut current_infected_counts_with_majority_r: HashMap<(usize, usize), usize> = HashMap::new();
            let mut current_infected_counts_total: HashMap<usize, usize> = HashMap::new();

            // --- NEW: Counters for requested statistics ---
            let mut individuals_with_any_bacterial_infection = 0;
            let mut individuals_with_any_r_positive_for_any_bacteria = 0;
            // --- END NEW ---


            for individual in self.population.individuals.iter() {
                // --- NEW: Flags for this individual's status ---
                let mut individual_has_any_infection = false;
                let mut individual_has_any_r_positive = false;
                // --- END NEW ---

                for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                    // Only count if currently infected
                    if individual.level[b_idx] > 0.001 {
                        // --- NEW: Set flag if individual has an infection ---
                        individual_has_any_infection = true;
                        // --- END NEW ---

                        *current_infected_counts_total.entry(b_idx).or_insert(0) += 1;

                        for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                            let resistance_data = &individual.resistances[b_idx][d_idx];
                            if resistance_data.majority_r > 0.0 {
                                current_majority_r_positive_values_by_combo
                                    .entry((b_idx, d_idx))
                                    .or_insert_with(Vec::new)
                                    .push(resistance_data.majority_r);
                                *current_infected_counts_with_majority_r.entry((b_idx, d_idx)).or_insert(0) += 1;
                            }
                            // --- NEW: Check for any_r > 0 for ANY bacteria/drug combo for this individual ---
                            if resistance_data.any_r > 0.0 {
                                individual_has_any_r_positive = true;
                            }
                            // --- END NEW ---
                        }
                    }
                }
                // --- NEW: Increment overall counters for the individual AFTER checking all their infections ---
                if individual_has_any_infection {
                    individuals_with_any_bacterial_infection += 1;
                }
                if individual_has_any_r_positive {
                    individuals_with_any_r_positive_for_any_bacteria += 1;
                }
                // --- END NEW ---
            }

            // --- Update global_majority_r_proportions based on the current step's data ---
            let strep_pneu_idx = self.bacteria_indices["strep_pneu"];
            let amoxicillin_idx = self.drug_indices["amoxicillin"];
            let infected_with_strep_pneu = *current_infected_counts_total.get(&strep_pneu_idx).unwrap_or(&0);
            let strep_pneu_amox_majority_r_count = *current_infected_counts_with_majority_r.get(&(strep_pneu_idx, amoxicillin_idx)).unwrap_or(&0);

            let proportion = if infected_with_strep_pneu > 0 {
                strep_pneu_amox_majority_r_count as f64 / infected_with_strep_pneu as f64
            } else {
                0.0
            };
            self.global_majority_r_proportions.insert((strep_pneu_idx, amoxicillin_idx), proportion);
            strep_pneu_amox_majority_r_history.push(proportion);


            let kleb_pneu_idx = self.bacteria_indices["kleb_pneu"];
            let ceftriaxone_idx = self.drug_indices["ceftriaxone"];
            let infected_with_kleb_pneu = *current_infected_counts_total.get(&kleb_pneu_idx).unwrap_or(&0);
            let kleb_pneu_ceftr_majority_r_count = *current_infected_counts_with_majority_r.get(&(kleb_pneu_idx, ceftriaxone_idx)).unwrap_or(&0);

            let proportion = if infected_with_kleb_pneu > 0 {
                kleb_pneu_ceftr_majority_r_count as f64 / infected_with_kleb_pneu as f64
            } else {
                0.0
            };
            self.global_majority_r_proportions.insert((kleb_pneu_idx, ceftriaxone_idx), proportion);
            kleb_pneu_ceftr_majority_r_history.push(proportion);

            // --- Logging for Individual 0
            let individual_0 = &self.population.individuals[0];
            println!("       Time step {}: Individual 0 age = {} days", t, individual_0.age);


            // 1. New block: Print drug details for individual 0, regardless of infection status
            println!("          --- Drug Presence Details (Individual 0) ---");
            let mut drugs_present_found_overall = false; // Declare and initialize here
            for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                if individual_0.cur_level_drug[drug_idx] > 0.0 {
                    let status = if individual_0.cur_use_drug[drug_idx] {
                        " (Currently being taken)"
                    } else {
                        " (Decaying)"
                    };
                    println!("             - {}: Level = {:.4}{}", drug_name_static, individual_0.cur_level_drug[drug_idx], status);
                    drugs_present_found_overall = true; // Use the newly declared variable
                }
            }
            if !drugs_present_found_overall {
                println!("              None (no antibiotics currently in system for Individual 0).");
            }


            println!("          --- Bacteria Infection Details (Individual 0) ---");
            let mut has_infection = false;
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                let level = individual_0.level[b_idx];
                if level > 0.0001 {
                    has_infection = true;
                    println!("             Bacteria: {}", bacteria_name);
                    println!("               Infected = true");
                    println!("               Level = {:.4}", level);
                    println!("               Immune Response = {:.4}", individual_0.immune_resp[b_idx]);
                    println!("               Infection From Environment = {}", individual_0.cur_infection_from_environment[b_idx]);
                    println!("               Hospital Acquired Infection = {}", individual_0.infection_hospital_acquired[b_idx]);
                    println!("               Test Identified Infection = {}", individual_0.test_identified_infection[b_idx]);
                    let mut drugs_present_found = false;
                    println!("                 Antibiotics Present in System (Current Level > 0):");
                    for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                        if individual_0.cur_level_drug[drug_idx] > 0.0 {
                            let status = if individual_0.cur_use_drug[drug_idx] {
                                " (Currently being taken)"
                            } else {
                                " (Decaying)"
                            };
                            println!("                   - {}: Level = {:.4}{}", drug_name_static, individual_0.cur_level_drug[drug_idx], status);
                            drugs_present_found = true;
                        }
                    }
                    if !drugs_present_found {
                        println!("                     None (no antibiotics currently in system).");
                    }
                    let mut effective_antibiotics_found = false;
                    println!("                 Antibiotic Resistance (Any_R) in System against {}:", bacteria_name);
                    for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                        if individual_0.cur_level_drug[drug_idx] > 0.0 {
                            let resistance_data = &individual_0.resistances[b_idx][drug_idx];
                            println!(
                                "                   - {}: Level = {:.4}, Any_R = {:.4}, Activity_R = {:.4}, Majority_R = {:.4}",
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
                        println!("                     None (no effective antibiotics in system against this bacteria).");
                    }
                    println!();
                }
            }
            if !has_infection {
                println!("             Individual 0 has no active bacterial infections.");

            }
            println!("       --------------------------------------------");
            println!("       -------------------------------------");

            // --- NEW SECTION: Global Infection and Resistance Statistics Output ---
            let total_population_size = self.population.individuals.len();
            let proportion_any_r_positive = if individuals_with_any_bacterial_infection > 0 {
                individuals_with_any_r_positive_for_any_bacteria as f64 / individuals_with_any_bacterial_infection as f64
            } else {
                0.0
            };

            println!("\n--- Global Infection and Resistance Statistics (Time step {}) ---", t);
            println!("  Total individuals in population: {}", total_population_size);
            println!("  Number of individuals with any bacterial infection: {}", individuals_with_any_bacterial_infection);
            println!("  Number of individuals with any bacteria having any_r > 0: {}", individuals_with_any_r_positive_for_any_bacteria);
            println!("  Proportion of infected individuals with any_r > 0: {:.4}\n", proportion_any_r_positive);
            // --- END NEW SECTION ---

            let step_duration = step_start.elapsed();
            println!("Time step {} took {:.3?} seconds", t, step_duration);
        }

        println!("--- SIMULATION FINISHED (within run function) ---");
        println!("--- SIMULATION ENDED ---");
    }
}