// src/simulation/simulation.rs
use crate::simulation::population::Population;
use crate::rules::apply_rules; // Corrected import path

/// Runs the simulation for a given population and number of time steps.
pub fn run(population: &mut Population, num_time_steps: usize, bacteria_to_track: &str) {
    println!("--- SIMULATION RUNNING ---");

    for step in 0..num_time_steps {
        // Apply rules to each individual
        for ind in &mut population.individuals {
            apply_rules(ind, step);

            // Print the values for the specified bacteria for individual 0
            if ind.id == 0 {
                println!("Time step {}: Individual 0 age = {} days", step, ind.age);
                println!("  Checking for bacteria: {}", bacteria_to_track);
                if ind.level.contains_key(bacteria_to_track) {
                    if let Some(level) = ind.level.get(bacteria_to_track) {
                        println!("    {}: level = {:.2}", bacteria_to_track, level);
                    }
                } else {
                    println!("    {}: level key NOT found!", bacteria_to_track);
                }
                if ind.immune_resp.contains_key(bacteria_to_track) {
                    if let Some(immune_resp) = ind.immune_resp.get(bacteria_to_track) {
                        println!("    {}: immune_resp = {:.2}", bacteria_to_track, immune_resp);
                    }
                } else {
                    println!("    {}: immune_resp key NOT found!", bacteria_to_track);
                }
                if ind.sepsis.contains_key(bacteria_to_track) {
                    if let Some(sepsis) = ind.sepsis.get(bacteria_to_track) {
                        println!("    {}: sepsis = {}", bacteria_to_track, sepsis);
                    }
                } else {
                    println!("    {}: sepsis key NOT found!", bacteria_to_track);
                }
                if ind.infectious_syndrome.contains_key(bacteria_to_track) {
                    if let Some(infectious_syndrome) = ind.infectious_syndrome.get(bacteria_to_track) {
                        println!("    {}: infectious_syndrome = {:.2}", bacteria_to_track, infectious_syndrome);
                    }
                } else {
                    println!("    {}: infectious_syndrome key NOT found!", bacteria_to_track);
                }
                if ind.date_last_infected.contains_key(bacteria_to_track) {
                    if let Some(date_last_infected) = ind.date_last_infected.get(bacteria_to_track) {
                        println!("    {}: date_last_infected = {}", bacteria_to_track, date_last_infected);
                    }
                } else {
                    println!("    {}: date_last_infected key NOT found!", bacteria_to_track);
                }
            }
        }
    }

    println!("--- SIMULATION FINISHED ---");
}