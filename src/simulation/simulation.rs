
// src/simulation/simulation.rs
use crate::simulation::population::Population;
use crate::rules::apply_rules; // Corrected import path


/// Runs the simulation for a given population and number of time steps.
pub fn run(population: &mut Population, num_time_steps: usize) {
    println!("--- SIMULATION RUNNING ---");

    for step in 0..num_time_steps {
        // Apply rules to each individual
        for ind in &mut population.individuals {
            apply_rules(ind, step);
        }

        // Print age in days of individual with ID 0
        if let Some(ind0) = population.get_individual(0) {
            println!("Time step {}: Individual 0 age = {} days", step, ind0.age);
        }
    }

    println!("--- SIMULATION FINISHED ---");
}
