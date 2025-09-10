
// src/main.rs
//
// Entry point for the AMR simulation.
//
// Responsibilities:
//   - Sets up and runs the simulation using parameters from config.rs
//   - Handles some initial and final reporting
//
// to run on less cores type this in before running: $env:RAYON_NUM_THREADS = "4"
//
//
//
 
mod simulation;
mod rules;
mod config;

// 
// -- additional output graphs ---------------------------------------------------------------------------------
//  
// graphs showing for each bacteria proportion of infections with each drug resistance mechanism present 
//
// there must be more things ...................
//
//
//
//
//
// -- model structure developments to consider ------------------------------------------------------------
//
// fix fig 9 plots of failure and number of drugs - consider including drug failure date (bacteria x drug) in 
// the "individual" struct
//
// consider if need to have non-sepsis death risks with specific bacterial infections
//
// need to better account for how death rate can be high in some infections despite treatment
// (and I think this means despite treatment with a non resistant drug)
//
// need to add fidaxomicin as a drug ?
//
// make sure risk of c diff with use of many antibiotics is accounted for, including its risk of recurrence
//
// consider higher doses and longer courses
//
// do we need a variable for whether under care for infection, as a prerequisite 
// for drug initiation even if no tests done
//
// do we need parameters like multi_drug_penalty_for_partial_cross_resistance given the 
// way we are calculating a total activity r which is supposed to be able to consider 2 drugs
// being taken
//
// investigate this "normalized current drug level as a proxy for 'activity_r' when any_r is 0" in mod.rs
//
// think if model able to capture why prior use of broad-spectrum antibiotics (especially cephalosporins, vancomycin, 
// carbapenems) strongly select for VREfm (enterococcus faecium resistance) 
//
// bear in mind that strep pneu for example has a vaccine against it but this has resulted in
// growth of non-vaccine-covered serotypes
//
// ? need to have a variable to indicate severity of infection - this can determine which drugs to use 
// - or is infectious syndrome sufficient ?
//     
// should syndrome influence sepsis risk independent of bacteria ?
//
// when people on two drugs account for drug drug interactions in determining final drug level of each drug
//
//
//
//
//
//
// calibration data: approx drug usage per 100_000 per calendar year 
//                   incidence of infection with each bacteria by age and calendar year
//                   deaths from each bacteria per 100_000 by region and calendar year
//                   resistance distribution for each used drug for each bacteria by calendar year  
//
// https://ourworldindata.org/antibiotics#:~:text=The%20map%20below%20shows%20the%20data%20collected%20by%20the%20World,(DDDs)%20per%201%2C000%20people.
//
// add age and region-specific all cause death rates from wpp/who and try to subtract bacterial  
// infection rates so they are background death rates
//
//
//
//
// set up automated testing for the simulation (probably not yet though)
//
//
//
//
// 
//
// decide on time zero for mda azithromycin project
//
// work on initial age distribution to reflect start year and end year and population growth - decide on start and end year
// for azithromycin mda project
//
// mda with azithromycin is to reduce community incidence as well as treat existing
// infection
//
// for mda project can base in africa with an "other" region all groued together
//
// to(maybe)do: perhaps introduce an effect whereby drug treatment leads to an increase in risk of microbiome_r > 0 due to   
//              allowing more bacteria growth due to killing other bacteria in microbiome, and so can be caused by any drug 
//              - but not sure yet if this is needed / justified
//
// consider adding tb, consider adding fungi
//
// 

use crate::simulation::simulation::Simulation;
 
fn main() {
    // Create and run the simulation
    let population_size =  5_000; 
    let time_steps =  35_000 ;  // 38_325
    let log_individuals = false ; // Set to false to disable detailed individual logging

    let mut simulation = Simulation::new(population_size, time_steps, log_individuals);

    use std::time::Instant;
    let start = Instant::now();

    simulation.run();

    let duration = start.elapsed();
    
    // Print summary statistics from logged data
    simulation.print_summary_statistics();
    
    // DEVELOPMENT: Use a fixed filename for easier analysis in Python
    // NOTE: The random filename logic below is commented out for now. Restore for large-scale runs.
    // use rand::Rng;
    // let mut rng = rand::thread_rng();
    // let random_id: u32 = rng.gen_range(1_000_000..10_000_000);
    // let csv_filename = format!("simulation_summary_{}.csv", random_id);
    let csv_filename = "simulation_summary.csv".to_string();

    // Export to CSV for analysis
    if let Err(e) = simulation.export_summary_to_csv(&csv_filename) {
        println!("Error exporting CSV: {}", e);
    } else {
        println!("Summary data exported to {}", csv_filename);
    }

    println!("main.rs  final outputs ");


    // --- FINAL SUMMARY/STATISTICS/REPORTING SECTION ---
    // This section previously performed detailed summary, statistics, and plotting at the end of the simulation.
    // It has been commented out for performance reasons, as more information is now collected during the time steps.
    // You can uncomment this block if you need the final summary/statistics/plots again.
    // (See previous git history for the full code.)

// END OF ADDITIONAL FINAL PRINTOUTS */


    // Log the simulation run details
    if let Err(e) = log_simulation_run(population_size, time_steps, duration.as_secs_f64()) {
        eprintln!("Error logging simulation run: {}", e);
    }

    println!("\n--- simulation ended ---");
    println!("--- total simulation time: {:.3} seconds", duration.as_secs_f64());
    println!("                          ");


}

// Function to log simulation run details
fn log_simulation_run(population_size: usize, time_steps: usize, duration_secs: f64) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::OpenOptions;
    use std::io::Write;
    use chrono::Utc;
    
    let timestamp = Utc::now();
    let log_entry = format!(
        "{},{},{},{:.3}\n",
        timestamp.format("%Y-%m-%d %H:%M:%S UTC"),
        population_size,
        time_steps,
        duration_secs
    );
    
    // Check if log file exists, if not create it with headers
    let log_path = "simulation_run_log.csv";
    let file_exists = std::path::Path::new(log_path).exists();
    
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(log_path)?;
    
    // Write header if file is new
    if !file_exists {
        writeln!(file, "timestamp,population_size,time_steps,duration_seconds")?;
    }
    
    // Write the log entry
    file.write_all(log_entry.as_bytes())?;
    
    println!("Simulation run logged to {}", log_path);
    
    Ok(())
}