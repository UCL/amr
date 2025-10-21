
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
// checked total number of death with sepsis (~ 8 million in 2021 (model)  vs 14 million gbd)
// number on drug on any one given day 75 million (model) vs ~ 100 million roughly estimated
// percent with new bacterial infection per year (~ 4% (model) vs ~ 10% rough empiric estimate) 
//
//
// -- additional outputs / graphs ---------------------------------------------------------------------------------
//
// review output plots and consider any changes in format or additions
//
// maybe for each drug for each bacteria only plot those drugs with potency above a certain value  
//
// proportion currently infected with h pylori / mdr tb
//
// calculate total number of sepsis deaths in 2019 (pre covid) scaled up to the world population size to compare with gbd 
// 14.1 million
//
// proportion of sepsis deaths by syndrome (and region and age) to compare with gbd
//
// get out more graphs on microbiome presence and microbiome_r and check on whether
// happy with logic
//
//
//
//
// -- model structure developments to consider ------------------------------------------------------------
//
// seems often immunity level is increasing too rapidly ?
//
// review whether immune response is implausibly strong for some syndromes and bacteria, especially
// for young children or the elderly
//
// should we have possibility of syndrome progressing to bloodstream (with much higher sepsis risk ?)  
//
// ok that pseudomonas aeruginosa can continue for > 90 days with full immune response and no clearance ?
// 
// consider whether infection from the environment should also depend on concurrent majority_r 
//
// introduce a new mechanism of death which is "non sepsis but caused by infection"
// (for h pylori and mdr tb for example, which currently don't increase death risk)
//
// bring back in mic suseptible based potency measures which are commented out in config.rs - this will 
// require re-scaling of activity_r scales
//
// consider realistic / plausible microbiome levels by bacteria
//
// grouped figure 4d y-axis - wrong or expressed as percent ?
//
// need to better account for how death rate can be high in some infections despite treatment ?
// (and I think this means despite treatment with a non resistant drug)
//
// make sure risk of c diff with use of many antibiotics is accounted for, including its risk of recurrence
//
// think if model able to capture why prior use of broad-spectrum antibiotics (especially cephalosporins, vancomycin, 
// carbapenems) strongly select for VREfm (enterococcus faecium resistance) 
//
// 
//
//
//
// consider (but probably only for future iterations):
//
// reduced bacterial growth rates for resistant strains ?
// competition between sensitive and resistant strains in microbiome ?
// mechanism-specific resistance costs: high-cost (carbapenemase) vs. low-cost (point mutations)
// time-dependent costs: higher initially, decreasing with compensatory mutations
// multi-drug cost interactions: costs compound with multiple resistance mechanisms
// differentiate growth rates - fast vs. slow growing bacteria ?
// consider treatment phases - intensive vs. continuation therapy ?
// model dormancy - especially for chronic infections like mdr tb ?
// biofilm resistance - reduced drug effectiveness in chronic infections ?
// use gbd super regions instead of continents ?
// need to add fidaxomicin as a drug ?
// bear in mind that strep pneu for example has a vaccine against it but this has resulted in growth of non-vaccine-covered serotypes
// we may still decide we need to model drug-specific drug levels but not clear how we would get data
// consider whether infection from the environment should also depend on concurrent majority_r - or do we need
// to somehow model bacteria in the environment and the influences on them such as use of antibiotics..... ? 
// (maybe someone can do this in a future interation.....)
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
use crate::simulation::population::BACTERIA_LIST;
use crate::config::get_global_param;

fn main() {
    // Validate bacteria configuration
    validate_bacteria_configuration();
    
    // Create and run the simulation
    let population_size = 125_000 ; 
    let time_steps = 38_325    ;  
    let log_individuals = false  ; // Set to false to disable detailed individual logging
    let log_infection_journeys = false  ; // Set to true to enable infection journey logging
    let infection_journey_sample_rate = 0.001      ; // Log 1% of infections for analysis (0.0-1.0)
    let use_fixed_seed = false ; // Toggle to enable deterministic RNG seeding
    let fixed_seed_value: u64 = 1_234_567_890; // Seed used when use_fixed_seed is true
    let infection_journey_bacteria_filter: Option<&str> = None; // Set to Some("escherichia_coli") to log only specific bacteria
    
    // Examples of bacteria filter values (use lowercase with underscores):
    // Some("escherichia_coli")
    // Some("staphylococcus_aureus") 
    // Some("pseudomonas_aeruginosa")
    // Some("acinetobacter_baumannii")
    // Some("enterococcus_faecium")
    // None - logs all bacteria types

    let mut simulation = Simulation::new(population_size, time_steps, log_individuals);
    if use_fixed_seed {
        simulation.rng_seed = Some(fixed_seed_value);
    }
    
    // Configure infection journey logging
    if log_infection_journeys {
        // Enable journey logging with optional bacteria filter
        match infection_journey_bacteria_filter {
            Some(filter) => simulation.enable_infection_journey_logging_with_filter(
                infection_journey_sample_rate,
                Some(filter.to_string()),
            ),
            None => simulation.enable_infection_journey_logging(infection_journey_sample_rate),
        }
    }

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

/// Validates the current bacteria configuration and provides helpful warnings
fn validate_bacteria_configuration() {
    let num_bacteria = BACTERIA_LIST.len();
    
    println!("=== BACTERIA CONFIGURATION VALIDATION ===");
    println!("Number of bacteria in simulation: {}", num_bacteria);
    
    if num_bacteria == 0 {
        panic!("ERROR: BACTERIA_LIST cannot be empty!");
    }
    
    if num_bacteria == 1 {
        println!("⚠️  SINGLE-BACTERIA MODE: This limits biological realism but is valid for:");
        println!("   • Pathogen-specific resistance studies");
        println!("   • Drug development against specific organisms");
        println!("   • Educational/training scenarios");
        println!("   • Computational efficiency");
        println!("   Note: HGT, microbiome competition, and syndromic treatment are disabled.");
    } else if num_bacteria < 5 {
        println!("⚠️  LIMITED-BACTERIA MODE: Some ecosystem effects may be reduced");
        println!("   Consider if your research question needs more bacterial diversity");
    } else {
        println!("✓ MULTI-BACTERIA MODE: Full ecosystem modeling enabled");
    }
    
    // Check for potential HGT configuration issues
    if num_bacteria == 1 {
        for bacteria in BACTERIA_LIST.iter() {
            for other in BACTERIA_LIST.iter() {
                if bacteria != other {
                    let hgt_param_name = format!("hgt_prob_{}_to_{}", bacteria, other);
                    if let Some(hgt_prob) = get_global_param(&hgt_param_name) {
                        if hgt_prob > 0.0 {
                            println!("⚠️  HGT parameters configured but only 1 bacteria present - HGT will not occur");
                            break;
                        }
                    }
                }
            }
        }
    }
    
    println!("Bacteria included:");
    for (i, bacteria) in BACTERIA_LIST.iter().enumerate() {
        println!("   {}. {}", i + 1, bacteria);
    }
    println!("=====================================\n");
}