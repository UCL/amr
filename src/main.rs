// src/main.rs
mod rules;
mod simulation;

use simulation::population::Population;
use simulation::simulation::run; // Import the run function

fn main() {
    println!("--- AMR SIMULATION ---");

    let mut population = Population::new(30_000);
    let bacteria_to_track = "acinetobac_bau";

    println!(
        "Initial age of individual 0 BEFORE simulation: {} days",
        population.individuals[0].age
    );
    println!("--- SIMULATION STARTING ---");

    let num_time_steps = 2;
    run(&mut population, num_time_steps, bacteria_to_track); // Call the run function

    println!("--- SIMULATION ENDED ---");
}
