// src/simulation/simulation.rs
use crate::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES, Population};
use crate::rules::apply_rules;
use rayon::prelude::*;
use std::collections::HashMap;
use rand::seq::SliceRandom; // Still needed for individual rule application

/// Runs the simulation for a given population and number of time steps.
pub fn run(population: &mut Population, num_time_steps: usize, bacteria_to_track: &str) {
    println!("--- AMR SIMULATION (within run function) ---");
    if let Some(ind) = population.individuals.get(0) {
        println!("Initial age of individual 0 BEFORE simulation: {} days", ind.age);
    }
    println!("--- SIMULATION STARTING (within run function) ---");

    // Pre-calculate indices for efficiency
    let bacteria_indices: HashMap<&str, usize> = BACTERIA_LIST.iter().enumerate().map(|(i, &b)| (b, i)).collect();
    let drug_indices: HashMap<&str, usize> = DRUG_SHORT_NAMES.iter().enumerate().map(|(i, &d)| (d, i)).collect();

    // History for plotting/analysis (optional, but good for tracking)
    let mut infection_counts_history: HashMap<&'static str, Vec<usize>> = HashMap::new();
    let mut global_cr_proportion_history: HashMap<(usize, usize), Vec<f64>> = HashMap::new();

    for step in 0..num_time_steps {
        // --- GLOBAL RESISTANCE METRICS CALCULATION START ---
        // These HashMaps will hold the data for the current time step
        let mut total_infected_counts_by_bacteria: HashMap<usize, usize> = HashMap::new();
        let mut cr_positive_infected_counts_by_combo: HashMap<(usize, usize), usize> = HashMap::new();
        let mut cr_positive_values_by_combo: HashMap<(usize, usize), Vec<f64>> = HashMap::new();

        // Iterate through all individuals to collect global resistance data
        for ind in population.individuals.iter() {
            for (&bacteria_name, &level) in ind.level.iter() {
                if level > 0.0 {
                    if let Some(&b_idx) = bacteria_indices.get(bacteria_name) {
                        *total_infected_counts_by_bacteria.entry(b_idx).or_insert(0) += 1;

                        for (&drug_name, &d_idx) in drug_indices.iter() {
                            if let Some(resistance_row) = ind.resistances.get(b_idx) {
                                if let Some(resistance_data) = resistance_row.get(d_idx) {
                                    if resistance_data.majority_r > 0.0 {
                                        *cr_positive_infected_counts_by_combo.entry((b_idx, d_idx)).or_insert(0) += 1;
                                        cr_positive_values_by_combo.entry((b_idx, d_idx)).or_insert_with(Vec::new).push(resistance_data.majority_r);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Calculate global proportions for the current time step
        let mut current_global_cr_proportions: HashMap<(usize, usize), f64> = HashMap::new();
        for b_idx in 0..BACTERIA_LIST.len() {
            let total_infected_count = *total_infected_counts_by_bacteria.get(&b_idx).unwrap_or(&0);
            if total_infected_count > 0 {
                for d_idx in 0..DRUG_SHORT_NAMES.len() {
                    let cr_positive_count = *cr_positive_infected_counts_by_combo.get(&(b_idx, d_idx)).unwrap_or(&0);
                    let proportion = cr_positive_count as f64 / total_infected_count as f64;
                    current_global_cr_proportions.insert((b_idx, d_idx), proportion);

                    // Store history for analysis (optional)
                    global_cr_proportion_history.entry((b_idx, d_idx))
                        .or_insert_with(Vec::new)
                        .push(proportion);
                }
            } else {
                for d_idx in 0..DRUG_SHORT_NAMES.len() {
                    current_global_cr_proportions.insert((b_idx, d_idx), 0.0);
                    global_cr_proportion_history.entry((b_idx, d_idx))
                        .or_insert_with(Vec::new)
                        .push(0.0);
                }
            }
        }
        // --- GLOBAL RESISTANCE METRICS CALCULATION END ---

        // Apply rules to each individual in parallel, passing the global data
        population.individuals.par_iter_mut().for_each(|ind| {
            apply_rules(
                ind,
                step,
                &current_global_cr_proportions,
                &cr_positive_values_by_combo,
                &bacteria_indices,
                &drug_indices,
            );
        });

        // Store the count of infections *after* rules are applied for printing
        for &bacteria_name in BACTERIA_LIST.iter() {
            let current_infected_count = population.individuals.iter().filter(|ind| {
                ind.level.get(bacteria_name).map_or(false, |&level| level > 0.0)
            }).count();
            infection_counts_history.entry(bacteria_name)
                .or_insert_with(Vec::new)
                .push(current_infected_count);

            println!(
                "Time step {}: Total individuals infected with {} = {}",
                step, bacteria_name, current_infected_count
            );
        }

        // Print some sample global proportions for inspection
        println!("Time step {}: Global majority_r Proportions (selected):", step);
        if let Some(&b_idx_strep) = bacteria_indices.get("strep_pneu") {
            if let Some(&d_idx_amox) = drug_indices.get("amoxicillin") {
                if let Some(&prop) = current_global_cr_proportions.get(&(b_idx_strep, d_idx_amox)) {
                    println!("    Strep Pneumonia to Amoxicillin: {:.4}", prop);
                }
            }
        }
        if let Some(&b_idx_generic) = bacteria_indices.get("generic_bacteria") {
            if let Some(&d_idx_amox) = drug_indices.get("amoxicillin") {
                if let Some(&prop) = current_global_cr_proportions.get(&(b_idx_generic, d_idx_amox)) {
                    println!("    Generic Bacteria to Amoxicillin: {:.4}", prop);
                }
            }
        }


        // Print the values for the specified bacteria for individual 0 AFTER applying rules
        if let Some(ind) = population.individuals.get(0) {
            println!("    Time step {}: Individual 0 age = {} days", step, ind.age);
            if let Some(level) = ind.level.get(bacteria_to_track) {
                println!("    {}: level = {:.2}", bacteria_to_track, level);
            }
            if let Some(immune_resp) = ind.immune_resp.get(bacteria_to_track) {
                println!("    {}: immune_resp = {:.2}", bacteria_to_track, immune_resp);
            }
            if let Some(sepsis) = ind.sepsis.get(bacteria_to_track) {
                println!("    {}: sepsis = {}", bacteria_to_track, sepsis);
            }
            if let Some(infectious_syndrome) = ind.infectious_syndrome.get(bacteria_to_track) {
                println!("    {}: infectious_syndrome = {}", bacteria_to_track, infectious_syndrome);
            }
            if let Some(date_last_infected) = ind.date_last_infected.get(bacteria_to_track) {
                println!("    {}: date_last_infected = {}", bacteria_to_track, date_last_infected);
            }

            if let Some(&from_env) = ind.cur_infection_from_environment.get(bacteria_to_track) {
                println!("    {}: cur_infection_from_environment = {}", bacteria_to_track, from_env);
            } else {
                println!("    {}: cur_infection_from_environment = Not applicable (no active infection)", bacteria_to_track);
            }

            if let Some(&hospital_acquired) = ind.infection_hospital_acquired.get(bacteria_to_track) {
                println!("    {}: infection_hospital_acquired = {}", bacteria_to_track, hospital_acquired);
            } else {
                println!("    {}: infection_hospital_acquired = Not applicable (no active infection)", bacteria_to_track);
            }

            if let Some(&test_identified) = ind.test_identified_infection.get(bacteria_to_track) {
                println!("    {}: test_identified_infection = {}", bacteria_to_track, test_identified);
            } else {
                println!("    {}: test_identified_infection = Not applicable (no active infection)", bacteria_to_track);
            }

            println!("    --- Drug Use Status (Individual 0) ---");
            for (i, &use_drug) in ind.cur_use_drug.iter().enumerate() {
                let drug_name = DRUG_SHORT_NAMES[i];
                let drug_level = ind.cur_level_drug[i];
                println!("      {}: cur_use_drug = {}, cur_level_drug = {:.2}", drug_name, use_drug, drug_level);
            }
            println!("    -------------------------------------");
        }
    }

    println!("--- SIMULATION FINISHED (within run function) ---");

    // Example: Print summary of global proportions for Strep Pneumonia to Amoxicillin
    if let (Some(&b_idx_strep), Some(&d_idx_amox)) = (bacteria_indices.get("strep_pneu"), drug_indices.get("amoxicillin")) {
        if let Some(history) = global_cr_proportion_history.get(&(b_idx_strep, d_idx_amox)) {
            println!("\n--- Proportion of Strep Pneumonia Infected with majority_r > 0 for Amoxicillin Over Time ---");
            for (time, proportion) in history.iter().enumerate() {
                println!("Time Step {}: {:.4}", time, proportion);
            }
        }
    }
}