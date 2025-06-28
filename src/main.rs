// src/main.rs

mod simulation;
mod rules;
mod config;


//
// order: 
// infection with no immune response, multiple concurrent infections
// infection with immune response, leading to eradication of infection
// infection with immune response and drug (but no drug resistance)
// infection with immune_response, drug and drug resistance
// start getting out graphs
//
//






use crate::simulation::simulation::Simulation;

fn main() {
    // Create and run the simulation
    let population_size = 1 ;
    let time_steps = 50;

    let mut simulation = Simulation::new(population_size, time_steps);

    // Initial state check for Individual 0 (from main.rs)
    println!("--- Initial state of Individual 0 (from main.rs) ---");
    let ind0 = &simulation.population.individuals[0];
    println!("  ID: {}", ind0.id);
    println!("  Age: {} days", ind0.age);
    println!("  Sex: {}", ind0.sex_at_birth);
    println!("  Region Living: {:?}", ind0.region_living);
    println!("  Region Currently In: {:?}", ind0.region_cur_in);

    // Print vaccination status for all bacteria
    println!("  --- Vaccination Status ---");
    for (bacteria, &b_idx) in simulation.bacteria_indices.iter() {
        println!("    {}_vaccination_status: {}", bacteria, ind0.vaccination_status[b_idx]);
    }
    let strep_pneu_idx = simulation.bacteria_indices["strep_pneu"];
    println!("  strep_pneu: level = {:.2}", ind0.level[strep_pneu_idx]);
    println!("  strep_pneu: immune_resp = {:.2}", ind0.immune_resp[strep_pneu_idx]);
    println!("  strep_pneu: sepsis = {}", ind0.sepsis[strep_pneu_idx]);
    println!("  strep_pneu: infectious_syndrome = {}", ind0.infectious_syndrome[strep_pneu_idx]);
    println!("  strep_pneu: date_last_infected = {}", ind0.date_last_infected[strep_pneu_idx]);
    println!("  strep_pneu: cur_infection_from_environment = {}", ind0.cur_infection_from_environment[strep_pneu_idx]);
    println!("  strep_pneu: infection_hospital_acquired = {}", ind0.infection_hospital_acquired[strep_pneu_idx]);
    println!("  strep_pneu: test_identified_infection = {}", ind0.test_identified_infection[strep_pneu_idx]);
    let amoxicillin_idx = simulation.drug_indices["amoxicillin"];
    println!("  cur_use_amoxicillin: {}", ind0.cur_use_drug[amoxicillin_idx]);
    println!("  cur_level_amoxicillin: {:.2}", ind0.cur_level_drug[amoxicillin_idx]);

    println!("  current_infection_related_death_risk: {:.2}", ind0.current_infection_related_death_risk);
    println!("  background_all_cause_mortality_rate: {:.4}", ind0.background_all_cause_mortality_rate);
    println!("  sexual_contact_level: {:.2}", ind0.sexual_contact_level);
    println!("  airborne_contact_level_with_adults: {:.2}", ind0.airborne_contact_level_with_adults);
    println!("  airborne_contact_level_with_children: {:.2}", ind0.airborne_contact_level_with_children);
    println!("  oral_exposure_level: {:.2}", ind0.oral_exposure_level);
    println!("  mosquito_exposure_level: {:.2}", ind0.mosquito_exposure_level);
    println!("  current_toxicity: {:.2}", ind0.current_toxicity);
    println!("  mortality_risk_current_toxicity: {:.2}", ind0.mortality_risk_current_toxicity);

    let resistance_data = &ind0.resistances[strep_pneu_idx][amoxicillin_idx];
    println!("  strep_pneu resistance to amoxicillin:");
    println!("    microbiome_r: {:.2}", resistance_data.microbiome_r);
    println!("    test_r: {:.2}", resistance_data.test_r);
    println!("    activity_r: {:.2}", resistance_data.activity_r);
    println!("    any_r: {:.2}", resistance_data.any_r);
    println!("    majority_r: {:.2}", resistance_data.majority_r);
    println!("-------------------------------------------");

    simulation.run();

    println!("\n--- SIMULATION RESULTS ---");

    // --- DEATH REPORTING START ---
    let mut total_deaths = 0;
    let mut death_causes_count: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    for individual in &simulation.population.individuals {
        if let Some(_date_of_death) = individual.date_of_death {
            total_deaths += 1;
            if let Some(cause) = &individual.cause_of_death {
                *death_causes_count.entry(cause.clone()).or_insert(0) += 1;
            }
        }
    }

    println!("\nTotal Deaths during simulation: {}", total_deaths);
    println!("Breakdown by Cause of Death:");
    for (cause, count) in death_causes_count {
        println!("  {}: {}", cause, count);
    }
    // --- DEATH REPORTING END ---

    println!("\n--- SIMULATION FINISHED ---");
}




