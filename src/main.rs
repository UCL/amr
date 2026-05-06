
//
//
//
//  ----  plan for paper -----------------------------------------------------------------------------------------
//
//  get at least three calibration parameter sets - from 2026 three "scenarios" to 2035
//  (i) continue status quo (ii) amr disappears (iii) immediate ast - or as a minimum find 3 distinct parameter sets
//  that calibrate
//
//
//
//
//
//
//
// -- thoughts 1 (needed for paper) ---------------------------------------------------------------------------------
//  
// present a global antibiotic activity for each bacteria over time
// - sum across drugs existing at the time of mean potency x (1 - any_r) amongst those infected with the bacteria
//
//  ? need to include temporal trends in incidence and drug choice for some bacteria ?
//
//
//
//
//
// -- thoughts 2 (consider for paper) ------------------------------------------------------------------------------
//  
// separate age categories of 70-80 and 80+ throughout ?
//
// have increased risk of strep b, e coli and listeria in women of childbearing age ?
//
// additional calibration output - show case fatality by infection if untreated (based )
//
// ? drug start distribution by infection site (or is seeing the drug choice rules enough ?)
//
// death within 30 days by bacteria, syndrome, age and region ? - make a formal part of calibration 
// score or just present as an fyi ?
//
// infections with test_r done by x days (by region ?)
// or proportion of drug treatment days which is empiric (by region and hospital status ?)
//
// https://ourworldindata.org/antibiotics#:~:text=The%20map%20below%20shows%20the%20data%20collected%20by%20the%20World,(DDDs)%20per%201%2C000%20people.
//
//
//
//
//
//
// -- thoughts 3 (am thinking not needed for paper, but re-review before submission) ------------------------------------------------------------------
//  
// ? maybe publish online results from one run with fixed seed and also results using the same fixed seed but changing
// one parameter value only - so all can see influence of each parameter - could refer to this in supplementary material.
//
// look at data on effects of stewardship policies on resistance and see if model can re-produce
//
// ? consider whether can replicate recent antibiotic use as associated with resistance presence
//
// ? add mechanisms present by day to infection journeys
//
// ? need outputs showing distributions of mechansims present in infecting bacteria, including co-presence of 
// multiple mechansims
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
// add age and region-specific all cause death rates from wpp/who and try to subtract bacterial
// infection rates so they are background death rates
//
// decide on time zero for mda azithromycin project
//
// ? need poc test for drug level at infection site
//
// consider case for real time infection site drug level monitoring 
//
// mda with azithromycin is to reduce community incidence as well as treat existing infection
//
// for mda project can base in africa with an "other" region all grouped together
//






// src/main.rs
// Simulation entry point.
//
// This file is intentionally small: it chooses one run configuration, builds a `Simulation`,
// executes it, and writes the main CSV output. Most model behaviour lives in `rules/`,
// `simulation/`, and `config.rs`; edit this file when you want to change how a run is launched.
//
// Useful local workflow note:
//   PowerShell: `$env:RAYON_NUM_THREADS = "4"` before running to cap parallelism.

//
//
//

use amr_project::config::get_global_param;
use amr_project::simulation::population::BACTERIA_LIST;
use amr_project::simulation::simulation::CalibrationMode;
use amr_project::simulation::simulation::Simulation;
use std::path::PathBuf;

fn main() {
    let _ = env_logger::builder().is_test(false).try_init();

    // Fail fast on obviously inconsistent bacteria setups before paying the cost of a full run.
    validate_bacteria_configuration();

    // Main run configuration. This is the quickest place to switch between calibration-sized
    // runs, full policy runs, deterministic debug runs, and journey-logging experiments.
    let population_size = 1_000_000;
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
    let log_individuals = false; // Full individual logging is expensive and mainly useful for narrow debugging.
    let log_infection_journeys = false; // Journey logging captures dense snapshots only for sampled infections.
    let infection_journey_sample_rate = 1.00; // Fraction of eligible infections to log when journey logging is enabled.
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
    let mut simulation = Simulation::new(
        population_size,
        time_steps,
        log_individuals,
        seed_override,
        calibration_mode,
    );
    let use_disk_branch_checkpointing = false; // Set to keep the branch checkpoint in memory
    let disk_checkpoint_directory: Option<PathBuf> = None; // Override with Some(path) to specify a custom folder

    if use_disk_branch_checkpointing {
        simulation.enable_disk_branch_checkpointing(disk_checkpoint_directory);
    } else {
        simulation.disable_disk_branch_checkpointing();
    }

    // Journey logging is optional because it writes a much richer, more expensive trace than the summary CSV.
    if log_infection_journeys {
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

    // Use the simulation's run_id in the filename so Python post-processing can join one
    // summary CSV to one set of sampled parameters and one run-log entry without collisions.
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

    // The summary CSV is the primary handoff to the Python analysis scripts.
    if let Err(e) = simulation.export_summary_to_csv(&csv_path) {
        println!("Error exporting CSV: {}", e);
    } else {
        println!("Summary data exported to {}", csv_path.display());
    }

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

/// Append one row to the lightweight run log so wall-clock cost can be compared across runs.
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

    // Create the file lazily so ad hoc runs do not need any manual setup.
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

/// Validate the bacteria roster before the simulation starts.
///
/// The model can technically run with one or a few organisms, but many ecosystem-level
/// assumptions (HGT, microbiome competition, syndromic empiric therapy) become much less
/// realistic, so this prints warnings rather than silently accepting a misleading setup.
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

