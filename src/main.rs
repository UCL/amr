mod simulation;
mod rules;

use simulation::population::Population;

fn main() {
    println!("--- AMR SIMULATION ---");

    let mut population = Population::new(30_000);

    if let Some(ind0) = population.get_individual(0) {
        println!("Initial age of individual 0 BEFORE simulation: {} days", ind0.age);
    }

    let num_time_steps = 10;
    simulation::simulation::run(&mut population, num_time_steps);

    println!("--- SIMULATION ENDED ---");
}





/* 

// src/main.rs
mod simulation; // Declares a module named `simulation`. Rust will look for a file named `simulation.rs` or a directory named `simulation` with a `mod.rs` file inside.
mod rules; // Declares a module named `rules`. Rust will look for a file named `rules.rs` or a directory named `rules` with a `mod.rs` file inside.

use simulation::population::Population; // Brings the `Population` struct into the current scope from the `population` module, which is inside the `simulation` module.

fn main() { // The main function, the entry point of the executable program.
    println!("--- AMR SIMULATION ---"); // Prints a banner to the console indicating the start of the simulation.

    let mut population = Population::new(30_000); // Creates a mutable variable named `population` and initializes it with a new instance of the `Population` struct, likely with an initial size of 30,000 individuals.
    let num_time_steps = 10; // Declares a constant variable `num_time_steps` and sets it to 10, representing the number of simulation steps to run.

    simulation::simulation::run(&mut population, num_time_steps); // Calls the `run` function located in the `simulation::simulation` module (likely a file named `simulation.rs` inside the `simulation` directory). It passes a mutable reference to the `population` and the `num_time_steps` as arguments, starting the simulation.

    println!("--- SIMULATION ENDED ---"); // Prints a message to the console indicating the end of the simulation after the `run` function completes.
}

*/
