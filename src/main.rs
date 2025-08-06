// src/main.rs
//
// Entry point for the AMR simulation.
//
// Responsibilities:
//   - Sets up and runs the simulation using parameters from config.rs
//   - Handles some initial and final reporting
//


// note:
// when need to follow variable values over time steps for individual 0 
// make the change shown at the top of simulation.rs and rules/mod.rs    
// make the change in population.rs to restrict to small number of bacteria and drugs
// decide which variable values to print out from the list in simulation.rs
// run up to a certain point and should be able to see the previous e.g. 10 time step values
//

mod simulation;
mod rules;
mod config;


//
// // model structure developments to consider //
//
//    see todo in rules/mod.rs 
//
//    ✓ IMPLEMENTED: increased risk of infection with certain bacteria in people currently hospitalized
// 
//    adding resistance mechanisms - steps - add risk of each mechanism appearing 
//    for each bacteria, which will depend partially on drug level as for any_r 
//    appearance - keep all any_r code as is as this will remain the default 
//    mechanism - allow presence of mechanism to over-write the any_r value 
//
// 
// // parameter values (recognising there will be many changes) //
//
//    e coli seems likely to be present in the microbiome of all individuals
//
//    update infection site distribution per bacteria
//
//
//  
//
//
//
// calibration data: approx drug usage per 100_000 per calendar year 
//                   incidence of infection with each bacteria by age and calendar year
//                   deaths from each bacteria per 100_000 by age and calendar year
//                   resistance distribution for each used drug for each bacteria by calendar year  
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
// for mda project can base in africa with an "other" region all groued together
//
// to(maybe)do: perhaps introduce an effect whereby drug treatment leads to an increase in risk of microbiome_r > 0 due to   
//              allowing more bacteria growth due to killing other bacteria in microbiome, and so can be caused by any drug 
//              - but not sure yet if this is needed / justified
//
// consider adding tb, consider adding fungi
//
// ? explicitly model resistance mechanisms and allow those to determine the any_r and majority_r values for each drug for 
// that bacteria - so this will be up to 11 mechanisms - there will be less than 11 variables per bacteria as each bacteria
// is only affected by a subset of the mechanisms - include fitness cost so the possibility that the mechanism is reversed 
// when the bacteria is not replicating in the presence of the drug - would still need the possibility of increases
// in any / majority_r by non-specific mechanisms - have not included this until now due to concern about all these
// mechanisms and others not being sufficiently well understood
//

use crate::simulation::simulation::Simulation;

fn main() {
    // Create and run the simulation
    let population_size =  3_000;
    let time_steps = 300 ; 
 
    let mut simulation = Simulation::new(population_size, time_steps);

    let ind0 = &simulation.population.individuals[0];
    
    // print variable values at time step 0, before starting to go through the time steps

    println!("  ");
    println!("main.rs  variable values at time step 0, before starting to go through the time steps");
    println!("  ");

    for (bacteria, &b_idx) in simulation.bacteria_indices.iter() {
        println!("{}_vaccination_status: {}", bacteria, ind0.vaccination_status[b_idx]);
    }

    println!("background_all_cause_mortality_rate: {:.4}", ind0.background_all_cause_mortality_rate);
    println!("sexual_contact_level: {:.2}", ind0.sexual_contact_level);
    println!("airborne_contact_level_with_adults: {:.2}", ind0.airborne_contact_level_with_adults);
    println!("airborne_contact_level_with_children: {:.2}", ind0.airborne_contact_level_with_children);
    println!("oral_exposure_level: {:.2}", ind0.oral_exposure_level);
    println!("mosquito_exposure_level: {:.2}", ind0.mosquito_exposure_level);
    println!("current_toxicity: {:.2}", ind0.current_toxicity);
    println!("mortality_risk_current_toxicity: {:.2}", ind0.mortality_risk_current_toxicity);

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

    // --- DEATH AND STATUS REPORTING START ---
    let mut total_deaths = 0;
    let mut death_causes_count: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    // New: Track per-bacteria and per-drug resistance counts
    let mut bacteria_infection_counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();

    // New: Track number in hospital and number severely immunosuppressed
    let mut number_in_hospital = 0;
    let mut number_severely_immunosuppressed = 0;

    for individual in &simulation.population.individuals {
        // Death reporting (existing)
        if let Some(_date_of_death) = individual.date_of_death {
            total_deaths += 1;
            if let Some(cause) = &individual.cause_of_death {
                *death_causes_count.entry(cause.clone()).or_insert(0) += 1;
            }
        }

        // Per-bacteria and per-drug resistance reporting
        for (bacteria, &b_idx) in simulation.bacteria_indices.iter() {
            if individual.level[b_idx] > 0.001 {
                // Count as infected with this bacteria
                *bacteria_infection_counts.entry(bacteria).or_insert(0) += 1;
            }
        }

        // Count number in hospital
        if individual.hospital_status.is_hospitalized() {
            number_in_hospital += 1;
        }
        // Count number severely immunosuppressed
        if individual.is_severely_immunosuppressed {
            number_severely_immunosuppressed += 1;
        }
    }

    // Write summary outputs to a file as well as printing
    use std::fs::OpenOptions;
    use std::io::Write;
    let mut report_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("simulation_report.txt")
        .expect("Unable to create/open simulation_report.txt");

    writeln!(report_file, "total deaths during simulation: {}", total_deaths).ok();
    println!("total deaths during simulation: {}", total_deaths);

    writeln!(report_file, "breakdown by cause of death:").ok();
    println!("breakdown by cause of death:");
    for (cause, count) in &death_causes_count {
        writeln!(report_file, "{}: {}", cause, count).ok();
        println!("{}: {}", cause, count);
    }

    writeln!(report_file, "number in hospital: {}", number_in_hospital).ok();
    println!("number in hospital: {}", number_in_hospital);

    writeln!(report_file, "number severely immunosuppressed: {}", number_severely_immunosuppressed).ok();
    println!("number severely immunosuppressed: {}", number_severely_immunosuppressed);



//  /*  ADDITIONAL FINAL PRINTOUTS

  
    // New: Print bacteria and resistance summary
    println!("\n--- Bacteria infection and resistance summary ---");
    for (bacteria, &count) in &bacteria_infection_counts {
        println!("{}: {} infected", bacteria, count);
        writeln!(report_file, "{}: {} infected", bacteria, count).ok();
        for (drug, _) in simulation.drug_indices.iter() {
            // Collect the full distribution of any_r for this bacteria/drug pair
            let mut any_r_values = Vec::new();
            for individual in &simulation.population.individuals {
                if let Some(&b_idx) = simulation.bacteria_indices.get(bacteria) {
                    if individual.level[b_idx] > 0.001 {
                        if let Some(&d_idx) = simulation.drug_indices.get(drug) {
                            let any_r = individual.resistances[b_idx][d_idx].any_r;
                            any_r_values.push(any_r);
                        }
                    }
                }
            }
            // Print and write summary statistics for the distribution
            if !any_r_values.is_empty() {
                let n = any_r_values.len() as f64;
                let mut count_0 = 0;
                let mut count_001_025 = 0;
                let mut count_0251_05 = 0;
                let mut count_0501_075 = 0;
                let mut count_0751_1 = 0;
                for &val in &any_r_values {
                    if val == 0.0 {
                        count_0 += 1;
                    } else if val > 0.0 && val <= 0.25 {
                        count_001_025 += 1;
                    } else if val > 0.25 && val <= 0.5 {
                        count_0251_05 += 1;
                    } else if val > 0.5 && val <= 0.75 {
                        count_0501_075 += 1;
                    } else if val > 0.75 && val <= 1.0 {
                        count_0751_1 += 1;
                    }
                }
                println!(
                    "    {}: n = {}, prop 0.00 = {:.3}, prop 0.25 = {:.3}, prop 0.5 = {:.3}, prop 0.75 = {:.3}, prop 1.00 = {:.3}",
                    drug,
                    n as usize,
                    count_0 as f64 / n,
                    count_001_025 as f64 / n,
                    count_0251_05 as f64 / n,
                    count_0501_075 as f64 / n,
                    count_0751_1 as f64 / n
                );
                writeln!(
                    report_file,
                    "    {}: n = {}, prop 0.00 = {:.3}, prop 0.25 = {:.3}, prop 0.5 = {:.3}, prop 0.75 = {:.3}, prop 1.00 = {:.3}",
                    drug,
                    n as usize,
                    count_0 as f64 / n,
                    count_001_025 as f64 / n,
                    count_0251_05 as f64 / n,
                    count_0501_075 as f64 / n,
                    count_0751_1 as f64 / n
                ).ok();
            } else {
                println!("    {}: n = 0", drug);
                writeln!(report_file, "    {}: n = 0", drug).ok();
            }
        }
    }
    // --- Bacteria/drug pairs: standardized_mic < 2 summary ---
    use crate::config;
    println!("\n--- Bacteria/drug pairs: standardized_mic < 2 summary ---");
    for (bacteria, &count) in &bacteria_infection_counts {
        for (drug, _) in simulation.drug_indices.iter() {
            let mut n_infected = 0;
            let mut n_standardized_mic_lt2 = 0;
            for individual in &simulation.population.individuals {
                if let (Some(&b_idx), Some(&d_idx)) = (simulation.bacteria_indices.get(bacteria), simulation.drug_indices.get(drug)) {
                    if individual.level[b_idx] > 0.001 {
                        n_infected += 1;
                        let resistance_data = &individual.resistances[b_idx][d_idx];
                        let potency_param_key = format!(
                            "drug_{}_for_bacteria_{}_potency_when_no_r",
                            drug, bacteria
                        );
                        let potency = config::get_global_param(&potency_param_key).unwrap_or(0.05);
                        let max_resistance_level = config::get_global_param("max_resistance_level").unwrap_or(1.0);
                        let normalized_majority_r = resistance_data.majority_r / max_resistance_level;
                        let standardized_mic = if (1.0 - normalized_majority_r) * potency > 0.0 {
                            1.0 / ((1.0 - normalized_majority_r) * potency)
                        } else {
                            f64::INFINITY
                        };
                        if standardized_mic < 2.0 {
                            n_standardized_mic_lt2 += 1;
                        }
                    }
                }
            }
            if n_infected > 0 {
                println!(
                    "    {} / {}: n_infected = {}, n_standardized_mic_lt2 = {}, prop = {:.3}",
                    bacteria, drug, n_infected, n_standardized_mic_lt2, n_standardized_mic_lt2 as f64 / n_infected as f64
                );
                writeln!(
                    report_file,
                    "    {} / {}: n_infected = {}, n_standardized_mic_lt2 = {}, prop = {:.3}",
                    bacteria, drug, n_infected, n_standardized_mic_lt2, n_standardized_mic_lt2 as f64 / n_infected as f64
                ).ok();
            }
        }
    }
      // --- end death and resistance reporting ---


    // Example: Plot distribution of any_r for one random bacteria-drug pair using plotters
    // (Requires plotters = "0.3" in Cargo.toml)
    use rand::seq::IteratorRandom;
    use plotters::prelude::*;

    // Pick a random bacteria-drug pair with at least one infected individual and at least one any_r > 0
    let mut rng = rand::thread_rng();
    let mut example_pair: Option<(&str, &str)> = None;
    let mut example_any_r_values = Vec::new();

    let pairs: Vec<(&str, &str)> = simulation.bacteria_indices.keys()
        .flat_map(|&bacteria| simulation.drug_indices.keys().map(move |&drug| (bacteria, drug)))
        .collect();

    // --- DEBUG: Print how many pairs have any_r > 0 ---
    let mut found_pairs = 0;
    for &(bacteria, drug) in &pairs {
        let mut values = Vec::new();
        if let (Some(&b_idx), Some(&d_idx)) = (simulation.bacteria_indices.get(bacteria), simulation.drug_indices.get(drug)) {
            for individual in &simulation.population.individuals {
                if individual.level[b_idx] > 0.001 {
                    let any_r = individual.resistances[b_idx][d_idx].any_r;
                    values.push(any_r);
                }
            }
            if values.iter().any(|&v| v > 0.0) {
                found_pairs += 1;
            }
        }
    }
    println!("Number of bacteria/drug pairs: {}", pairs.len());
    println!("Number of bacteria/drug pairs with any any_r > 0: {}", found_pairs);

    for &(bacteria, drug) in pairs.iter().choose_multiple(&mut rng, pairs.len()) {
        let mut values = Vec::new();
        if let (Some(&b_idx), Some(&d_idx)) = (simulation.bacteria_indices.get(bacteria), simulation.drug_indices.get(drug)) {
            for individual in &simulation.population.individuals {
                if individual.level[b_idx] > 0.001 {
                    let any_r = individual.resistances[b_idx][d_idx].any_r;
                    values.push(any_r);
                }
            }
            // Only use this pair if there is at least one value > 0
            if values.iter().any(|&v| v > 0.0) {
                example_pair = Some((bacteria, drug));
                example_any_r_values = values;
                break;
            }
        }
    }

    if let Some((bacteria, drug)) = example_pair {
        println!("\n--- Example histogram for any_r distribution: {} / {} ---", bacteria, drug);

        // Bin edges: [0, 0.25], (0.25, 0.5], (0.5, 0.75], (0.75, 1.0]
        let mut bins = [0; 5];
        for &val in &example_any_r_values {
            if val == 0.0 {
                bins[0] += 1;
            } else if val > 0.0 && val <= 0.25 {
                bins[1] += 1;
            } else if val > 0.25 && val <= 0.5 {
                bins[2] += 1;
            } else if val > 0.5 && val <= 0.75 {
                bins[3] += 1;
            } else if val > 0.75 && val <= 1.0 {
                bins[4] += 1;
            }
        }
        println!("Bin counts: {:?}", bins);

        // Always plot the histogram, even if only one bin is nonzero
        let root = BitMapBackend::new("any_r_histogram.png", (640, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let max_count = *bins.iter().max().unwrap_or(&1);

        // Set a minimum y-axis height for better visibility of small bars
        let y_axis_max = if max_count < 10 { 10 } else { max_count + 2 };

        // Use f64 for both axes in build_cartesian_2d and Rectangle coordinates
        let mut chart = ChartBuilder::on(&root)
            .caption(format!("any_r distribution for {} / {}", bacteria, drug), ("sans-serif", 18))
            .margin(40)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(
                0f64..5f64,
                0f64..(y_axis_max as f64),
            )
            .unwrap();

        chart
            .configure_mesh()
            .x_desc("any_r bin")
            .y_desc("Count")
            .disable_x_mesh() // Keep y-axis grid lines, but hide x-axis ones
            .disable_x_axis() // Hide the default x-axis line and labels
            .draw()
            .unwrap();

        // Draw all bins with equal width, centered on each bin midpoint
        use plotters::style::RGBColor;
        let bar_color = RGBColor(220, 50, 47); // Solarized red for visibility

        chart.draw_series(
            bins.iter().enumerate().map(|(i, &count)| {
                // Center the bar on the midpoint (i+0.5)
                let x0 = i as f64 + 0.05;
                let x1 = i as f64 + 0.95;
                Rectangle::new(
                    [(x0, 0.0), (x1, count as f64)],
                    if count > 0 { bar_color.filled() } else { WHITE.filled() },
                )
            }),
        ).unwrap();

        // Draw count labels above each bar for clarity, always at the center of the bar
        chart.draw_series(
            bins.iter().enumerate().map(|(i, &count)| {
                let x = i as f64 + 0.5;
                if count > 0 {
                    Text::new(
                        format!("{}", count),
                        (x, count as f64 + 0.5),
                        ("sans-serif", 15).into_font().color(&BLACK),
                    )
                } else {
                    Text::new(
                        String::new(),
                        (x, 0.0),
                        ("sans-serif", 15).into_font().color(&BLACK),
                    )
                }
            }),
        ).unwrap();

        // Manually draw the x-axis labels centered under each bar
        let labels = ["0", "0.25", "0.5", "0.75", "1"];
        chart.draw_series(
            labels.iter().enumerate().map(|(i, &label)| {
                let x = i as f64 + 0.5;
                Text::new(
                    label.to_string(),
                    (x, -0.05 * y_axis_max as f64), // Position below the x-axis
                    ("sans-serif", 15).into_font().color(&BLACK),
                )
            })
        ).unwrap();

        println!("Histogram saved to any_r_histogram.png");
    } else {
        println!("No bacteria/drug pair found with any nonzero any_r values.");
    }


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


