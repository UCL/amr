// src/main.rs

mod simulation;
mod rules;
mod config;

//
//
// decide on time zero for mda azithromycin project
//
// work on initial age distribution to reflect start year and end year and population growth - decide on start and end year
// for azithromycin mda project
//
// for mda project can base in africa with an "other" region all groued together
//
// calibration data: approx drug usage per 100_000 per calendar year 
//                   incidence of infection with each bacteria by age and calendar year
//                   deaths from each bacteria per 100_000 by age and calendar year
//                   resistance distribution for each used drug for each bacteria by calendar year  
//
// to(maybe)do: perhaps introduce an effect whereby drug treatment leads to an increase in risk of microbiome_r > 0 due to   
//              allowing more bacteria growth due to killing other bacteria in microbiome, and so can be caused by any drug 
//              - but not sure yet if this is needed / justified
//
// to consider in future: explicitly model resistance mechanisms and allow that to determine the any_r value for each drug for 
//                        that bacteria
//


use crate::simulation::simulation::Simulation;

fn main() {
    // Create and run the simulation
    let population_size =    10_000 ;
    let time_steps =  20  ;  // Reduced for testing immune response changes

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
    println!("main.rs  final outputs ");

    // --- DEATH REPORTING START ---
    let mut total_deaths = 0;
    let mut death_causes_count: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    // New: Track per-bacteria and per-drug resistance counts
    let mut bacteria_infection_counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();

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
    }

    println!("total deaths during simulation: {}", total_deaths);
    println!("breakdown by cause of death:");
    for (cause, count) in death_causes_count {
        println!("{}: {}", cause, count);
    }

/*


    // New: Print bacteria and resistance summary
    println!("\n--- Bacteria infection and resistance summary ---");
    for (bacteria, &count) in &bacteria_infection_counts {
        println!("{}: {} infected", bacteria, count);
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
            // Print summary statistics for the distribution
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
            } else {
                println!("    {}: n = 0", drug);
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
            .caption(format!("any_r distribution for {} / {}", bacteria, drug), ("sans-serif", 30))
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




*/


    println!("\n--- simulation ended ---");
    println!("\n--- total simulation time: {:.3?} seconds", duration);
    println!("                          ");


}


