// src/simulation/simulation.rs
use crate::simulation::population::{
    BACTERIA_LIST, DRUG_SHORT_NAMES, Population,
};
use crate::rules::apply_rules;
use rayon::prelude::*; // Bring in the parallel iterator traits

/// Runs the simulation for a given population and number of time steps.
pub fn run(population: &mut Population, num_time_steps: usize, bacteria_to_track: &str) {
    println!("--- AMR SIMULATION (within run function) ---");
    if let Some(ind) = population.individuals.get(0) {
        println!("Initial age of individual 0 BEFORE simulation: {} days", ind.age);
    }
    println!("--- SIMULATION STARTING (within run function) ---");

    // Find the index of gentamicin
    let gentamicin_index = DRUG_SHORT_NAMES.iter().position(|&drug| drug == "gentamicin");
    let gentamicin_present = gentamicin_index.is_some();
    let gentamicin_idx = gentamicin_index.unwrap_or(0);

    // Find the index of acinetobac_bau
    let acinetobac_index = BACTERIA_LIST.iter().position(|&bacteria| bacteria == "acinetobac_bau");
    let acinetobac_present = acinetobac_index.is_some();
    let acinetobac_idx = acinetobac_index.unwrap_or(0);

    for step in 0..num_time_steps {
        // Apply rules to each individual in parallel
        population.individuals.par_iter_mut().for_each(|ind| {
            apply_rules(ind, step);
            // Any other per-individual updates for this time step can go here
        });

        // Print the values for the specified bacteria for individual 0 AFTER applying rules
        if let Some(ind) = population.individuals.get(0) {
            println!("Time step {}: Individual 0 age = {} days", step, ind.age);
            if let Some(level) = ind.level.get(bacteria_to_track) {
                println!("  {}: level = {:.2}", bacteria_to_track, level);
            }
            // ... (rest of your printing code for individual 0) ...
            if acinetobac_present && gentamicin_present {
                let resistance = &ind.resistances[acinetobac_idx][gentamicin_idx];
                println!("  acinetobac_bau resistance to gentamicin:");
                println!("    microbiome_r: {:.2}", resistance.microbiome_r);
                println!("    test_r: {:.2}", resistance.test_r);
                println!("    activity_r: {:.2}", resistance.activity_r);
                println!("    e_r: {:.2}", resistance.e_r);
                println!("    c_r: {:.2}", resistance.c_r);
            } else {
                println!("  Could not find acinetobac_bau or gentamicin in the lists.");
            }
        }
    }

    println!("--- SIMULATION FINISHED (within run function) ---");
}