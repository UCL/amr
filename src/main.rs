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

mod config;
mod rules;
mod simulation;

// "Please use standard file edits so I can review the diff. Do not use terminal 
// commands or scripts to modify files."
//
//
//
//
//
//
//
//
// -- additional outputs / thoughts on calibration  ---------------------------------------------------------------------------------
//
// look at data on effects of stewardship policies on resistance and see if model can re-produce
//
// consider whether can replicate recent antibiotic use as associated with resistance presence
//
// add mechanisms present by day to infection journeys
//
// need outputs showing distributions of mechansims present in infecting bacteria, including co-presence of 
// multiple mechansims
//
// relative chance of drug start by infection site
//
// need targets for (i) % infections started in hospital by bacteria (ii) % of hospital infections with resistance by bacteria
// add to calibration summary: death rate from infection by bacteria
//
// death within 30 days by bacteria, syndrome, age and region ? - make a formal part of calibration 
// score or just present as an fyi ?
//
// infections with test_r done by x days (by region ?)
// or proportion of drug treatment days which is empiric (by region and hospital status ?)
//
// review choice for antibiotics in people not infected
//
//
//
// -- calibration approach:  
// maybe come up with ~ 10 different configs that lead to a resonable fit in different ways and run the 
// policy comparison several times on each 
//
//
//
//
//
// -- model structure developments to consider ------------------------------------------------------------
//
// footnote in calibration summary about any microbiome resistance vs majority
//
//
//
//
//
//
//
//
//
// consider for future iterations:
//
// model low level "treatment" resulting from antimicrobials in the environment ?
// consider having incidence of infection rising in situations if they occur in future in which infections cannot be treated
// reduced bacterial growth rates for resistant strains ?
// competition between sensitive and resistant strains in microbiome ?
// time-dependent costs: higher initially, decreasing with compensatory mutations
// multi-drug cost interactions: costs compound with multiple resistance mechanisms
// differentiate growth rates - fast vs. slow growing bacteria ?
// consider treatment phases - intensive vs. continuation therapy ?
// model dormancy - especially for chronic infections like mdr tb ?
// biofilm resistance - reduced drug effectiveness in chronic infections ?
// use gbd super regions instead of continents ?
// bear in mind that strep pneu for example has a vaccine against it but this has resulted in growth of non-vaccine-covered serotypes
// we may still decide we need to model drug-specific drug levels but not clear how we would get data
// consider whether infection from the environment should also depend on concurrent majority_r - or do we need
// to somehow model bacteria in the environment and the influences on them such as use of antibiotics..... ?
// (maybe someone can do this in a future interation.....)
// add legionella ?  tick-borne bacteria ? 
// remember that we model use of antibiotics when no modelled bacteria present so in some ways this takes care of bacteria not modelled 
// should population majority_r depend (more) on resistance in microbiome/carriage rather than infections ?
// should we consider some syndromes (from which spread is more likely) more than others for population majority_r ?
// consider more granular breakdown of regions
// consider adding tb, consider adding fungi
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
// set up automated testing for the simulation (?) (probably not yet though)
//
//
// comparisons:
//
// default combination therapy from (i) 1950 (ii) 1975 (iii) 2000 (iv) 2000 (v) 2025
// - consider modelling of lower replication capacity of resistant virus in absence of drug
// - consider if raise this what effect has on emergence rates 
// - consider effect of having the floor for some bacteria
// - consider various approaches to choice of the two drugs 
//
// other counterfactuals to include in first paper:
// - stop all antibiotic use  
// - drug can be started from infection
// - test_r known at infection and drug can be started
// - no further existence of resistance (already in)
//
// decide on time zero for mda azithromycin project
//
// work on initial age distribution to reflect start year and end year and population growth - decide on start and end year
// for azithromycin mda project
//
// ? need poc test for drug level at infection site
//
// consider case for real time infection site drug level monitoring 
//
// mda with azithromycin is to reduce community incidence as well as treat existing
// infection
//
// for mda project can base in africa with an "other" region all groued together
//
//
//
//

use crate::config::get_global_param;
use crate::simulation::population::BACTERIA_LIST;
use crate::simulation::simulation::Simulation;
use crate::simulation::simulation::CalibrationMode;
use std::path::PathBuf;

fn main() {
    let _ = env_logger::builder().is_test(false).try_init();

    // Validate bacteria configuration
    validate_bacteria_configuration();

    // Create and run the simulation
    let population_size = 100_000;
    // CalibrationMode::Full  — sparse CSV (2022-2025 only); fastest calibration runs.
    // CalibrationMode::Partial — all 1930-2025 rows kept; time-series plots still work.
    // CalibrationMode::None  — full run with policy branches to 2035.
    let calibration_mode = CalibrationMode::Full;
    // Calibration runs only need rows through the end of 2025.
    // 35_040 = 96 years * 365 days from 1930 to the start of 2026, so it covers 1930-2025 inclusive.
    // Full run (policy branches to 2035) needs 38_325.
    let time_steps = match calibration_mode {
        CalibrationMode::None => 38_325,
        CalibrationMode::Partial | CalibrationMode::Full => 35_040,
    };
    let log_individuals = false; // Set to false to disable detailed individual logging
    let log_infection_journeys = false  ; // Set to true to enable infection journey logging
    let infection_journey_sample_rate = 1.00; // Log 1% of infections for analysis (0.0-1.0)
    let use_fixed_seed = false; // Toggle to enable deterministic RNG seeding
    let fixed_seed_value: u64 = 1_234_567_890; // Seed used when use_fixed_seed is true
    let infection_journey_bacteria_filter: Option<&str> = None; // Set to Some("escherichia_coli") to log only specific bacteria

    // Examples of bacteria filter values (use lowercase with underscores):
    // Some("escherichia_coli")
    // Some("staphylococcus_aureus")
    // Some("pseudomonas_aeruginosa")
    // Some("acinetobacter_baumannii")
    // Some("enterococcus_faecium")
    // None - logs all bacteria types

    let seed_override = use_fixed_seed.then_some(fixed_seed_value);
    let mut simulation =
        Simulation::new(population_size, time_steps, log_individuals, seed_override, calibration_mode);
    let use_disk_branch_checkpointing = false; // Set to keep the branch checkpoint in memory
    let disk_checkpoint_directory: Option<PathBuf> = None; // Override with Some(path) to specify a custom folder

    if use_disk_branch_checkpointing {
        simulation.enable_disk_branch_checkpointing(disk_checkpoint_directory);
    } else {
        simulation.disable_disk_branch_checkpointing();
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
    let output_dir = std::path::Path::new("amr_simulation_output_analysis_outputs");
    if let Err(err) = std::fs::create_dir_all(output_dir) {
        eprintln!(
            "Warning: unable to create output directory {:?}: {}",
            output_dir, err
        );
    }
    let run_id = simulation.run_id;
    let csv_basename = format!("simulation_summary_{:06}.csv", run_id);
    let csv_path = output_dir.join(&csv_basename);

    // Export to CSV for analysis
    if let Err(e) = simulation.export_summary_to_csv(&csv_path) {
        println!("Error exporting CSV: {}", e);
    } else {
        println!("Summary data exported to {}", csv_path.display());
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
    println!(
        "--- total simulation time: {:.3} seconds",
        duration.as_secs_f64()
    );
    println!("                          ");
}

// Function to log simulation run details
fn log_simulation_run(
    population_size: usize,
    time_steps: usize,
    duration_secs: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    use chrono::Utc;
    use std::fs::OpenOptions;
    use std::io::Write;

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
        writeln!(
            file,
            "timestamp,population_size,time_steps,duration_seconds"
        )?;
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
