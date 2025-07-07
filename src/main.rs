// src/main.rs

mod simulation;
mod rules;
mod config;

//
//
//
// to(maybe)do: perhaps introduce an effect whereby drug treatment leads to increase in risk of microbiome_r > 0 due to allowing more 
//              bacteria growth due to killing other bacteria in microbiome, and so can be caused by any drug - but not sure yet if
//              this is needed / justified
//
// todo: review whether / how test_identified_infection is used
//
// start getting out graphs
//
// 



use crate::simulation::simulation::Simulation;

fn main() {
    // Create and run the simulation
    let population_size =   100_000 ;
    let time_steps = 50;

    let mut simulation = Simulation::new(population_size, time_steps);

    let ind0 = &simulation.population.individuals[0];
    
    // print variable values at time step 0, before starting to go through the time steps

    println!("  ");
    println!("main.rs  variable values at time step 0, before starting to go through the time steps");
    println!("  ");

    for (bacteria, &b_idx) in simulation.bacteria_indices.iter() {
        println!("{}_vaccination_status: {}", bacteria, ind0.vaccination_status[b_idx]);
    }
    let acinetobac_bau_idx = simulation.bacteria_indices["acinetobac_bau"];
    println!("acinetobac_bau: level = {:.2}", ind0.level[acinetobac_bau_idx]);
    println!("acinetobac_bau: immune_resp = {:.2}", ind0.immune_resp[acinetobac_bau_idx]);
    println!("acinetobac_bau: sepsis = {}", ind0.sepsis[acinetobac_bau_idx]);
    println!("acinetobac_bau: infectious_syndrome = {}", ind0.infectious_syndrome[acinetobac_bau_idx]);
    println!("acinetobac_bau: date_last_infected = {}", ind0.date_last_infected[acinetobac_bau_idx]);
    println!("acinetobac_bau: cur_infection_from_environment = {}", ind0.cur_infection_from_environment[acinetobac_bau_idx]);
    println!("acinetobac_bau: infection_hospital_acquired = {}", ind0.infection_hospital_acquired[acinetobac_bau_idx]);
    println!("acinetobac_bau: test_identified_infection = {}", ind0.test_identified_infection[acinetobac_bau_idx]);
    let cefepime_idx = simulation.drug_indices["cefepime"];
    println!("cur_use_cefepime: {}", ind0.cur_use_drug[cefepime_idx]);
    println!("cur_level_cefepime: {:.2}", ind0.cur_level_drug[cefepime_idx]);

    println!("current_infection_related_death_risk: {:.2}", ind0.current_infection_related_death_risk);
    println!("background_all_cause_mortality_rate: {:.4}", ind0.background_all_cause_mortality_rate);
    println!("sexual_contact_level: {:.2}", ind0.sexual_contact_level);
    println!("airborne_contact_level_with_adults: {:.2}", ind0.airborne_contact_level_with_adults);
    println!("airborne_contact_level_with_children: {:.2}", ind0.airborne_contact_level_with_children);
    println!("oral_exposure_level: {:.2}", ind0.oral_exposure_level);
    println!("mosquito_exposure_level: {:.2}", ind0.mosquito_exposure_level);
    println!("current_toxicity: {:.2}", ind0.current_toxicity);
    println!("mortality_risk_current_toxicity: {:.2}", ind0.mortality_risk_current_toxicity);

    let resistance_data = &ind0.resistances[acinetobac_bau_idx][cefepime_idx];
    println!("acinetobac_bau resistance to cefepime:");
    println!("microbiome_r: {:.2}", resistance_data.microbiome_r);
    println!("test_r: {:.2}", resistance_data.test_r);
    println!("activity_r: {:.2}", resistance_data.activity_r);
    println!("any_r: {:.2}", resistance_data.any_r);
    println!("majority_r: {:.2}", resistance_data.majority_r);

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
    let mut bacteria_drug_any_r_counts: std::collections::HashMap<(&str, &str), usize> = std::collections::HashMap::new();

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
                // For each drug, count if any_r > 0
                for (drug, &d_idx) in simulation.drug_indices.iter() {
                    if individual.resistances[b_idx][d_idx].any_r > 0.0 {
                        *bacteria_drug_any_r_counts.entry((bacteria, drug)).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    println!("total deaths during simulation: {}", total_deaths);
    println!("breakdown by cause of death:");
    for (cause, count) in death_causes_count {
        println!("{}: {}", cause, count);
    }

    // New: Print bacteria and resistance summary
    println!("\n--- Bacteria infection and resistance summary ---");
    for (bacteria, &count) in &bacteria_infection_counts {
        println!("{}: {} infected", bacteria, count);
        for (drug, _) in simulation.drug_indices.iter() {
            let any_r_count = *bacteria_drug_any_r_counts.get(&(bacteria, drug)).unwrap_or(&0);
            println!("    {}: {} with any_r > 0", drug, any_r_count);
        }
    }
    // --- end death and resistance reporting ---

    println!("\n--- simulation ended ---");
    println!("\n--- total simulation time: {:.3?} seconds", duration);
    println!("                          ");
}