// src/main.rs
mod simulation;
mod rules;
mod config;

// Corrected import path for Simulation
use crate::simulation::simulation::Simulation;
// Removed: use crate::config::{init_params}; // No init_params function in your config.rs
use crate::simulation::population::BACTERIA_LIST; // Import BACTERIA_LIST

fn main() {
    // Removed: init_params(); // No init_params function to call in your config.rs

    // Create and run the simulation
    let population_size = 100_000;
    let time_steps = 20;

    let mut simulation = Simulation::new(population_size, time_steps);

    // Initial state check for Individual 0 (from main.rs)
    println!("--- Initial state of Individual 0 (from main.rs) ---");
    let ind0 = &simulation.population.individuals[0];
    println!("  ID: {}", ind0.id);
    println!("  Age: {} days", ind0.age);
    println!("  Sex: {}", ind0.sex_at_birth);
    println!("  Region Living: {:?}", ind0.region_living);
    println!("  Region Visiting: {:?}", ind0.region_visiting);

    // Iterating through the HashMap for vaccination status
//  println!("  --- Vaccination Status ---");
//  for &bacteria in BACTERIA_LIST.iter() {
//      if let Some(&status) = ind0.vaccination_status.get(bacteria) {
//          println!("    {}_vaccination_status: {}", bacteria, status);
//      }
//  }

    if let Some(&level) = ind0.level.get("strep_pneu") {
        println!("  strep_pneu: level = {:.2}", level);
    }
    if let Some(&immune_resp) = ind0.immune_resp.get("strep_pneu") {
        println!("  strep_pneu: immune_resp = {:.2}", immune_resp);
    }
    if let Some(&sepsis) = ind0.sepsis.get("strep_pneu") {
        println!("  strep_pneu: sepsis = {}", ind0.sepsis.get("strep_pneu").unwrap_or(&false));
    }
    if let Some(&infectious_syndrome) = ind0.infectious_syndrome.get("strep_pneu") {
        println!("  strep_pneu: infectious_syndrome = {}", infectious_syndrome);
    }
    if let Some(&date_last_infected) = ind0.date_last_infected.get("strep_pneu") {
        println!("  strep_pneu: date_last_infected = {}", date_last_infected);
    }
    if let Some(&from_env) = ind0.cur_infection_from_environment.get("strep_pneu") {
        println!("  strep_pneu: cur_infection_from_environment = {}", from_env);
    }
    if let Some(&hospital_acquired) = ind0.infection_hospital_acquired.get("strep_pneu") {
        println!("  strep_pneu: infection_hospital_acquired = {}", hospital_acquired);
    }
    if let Some(&test_identified) = ind0.test_identified_infection.get("strep_pneu") {
        println!("  strep_pneu: test_identified_infection = {}", test_identified);
    }


    let drug_indices = &simulation.drug_indices; // Get drug_indices from simulation
    if let Some(&amox_idx) = drug_indices.get("amoxicillin") {
        println!("  cur_use_amoxicillin: {}", ind0.cur_use_drug[amox_idx]);
        println!("  cur_level_amoxicillin: {:.2}", ind0.cur_level_drug[amox_idx]);
    }


    println!("  current_infection_related_death_risk: {:.2}", ind0.current_infection_related_death_risk);
    println!("  background_all_cause_mortality_rate: {:.4}", ind0.background_all_cause_mortality_rate);
    println!("  sexual_contact_level: {:.2}", ind0.sexual_contact_level);
    println!("  airborne_contact_level_with_adults: {:.2}", ind0.airborne_contact_level_with_adults);
    println!("  airborne_contact_level_with_children: {:.2}", ind0.airborne_contact_level_with_children);
    println!("  oral_exposure_level: {:.2}", ind0.oral_exposure_level);
    println!("  mosquito_exposure_level: {:.2}", ind0.mosquito_exposure_level);
    println!("  under_care: {}", ind0.under_care);
    println!("  current_toxicity: {:.2}", ind0.current_toxicity);
    println!("  mortality_risk_current_toxicity: {:.2}", ind0.mortality_risk_current_toxicity);

    let strep_pneu_idx = simulation.bacteria_indices["strep_pneu"];
    let amoxicillin_idx = simulation.drug_indices["amoxicillin"];
    let resistance_data = &ind0.resistances[strep_pneu_idx][amoxicillin_idx];
    println!("  strep_pneu resistance to amoxicillin:");
    println!("    microbiome_r: {:.2}", resistance_data.microbiome_r);
    println!("    test_r: {:.2}", resistance_data.test_r);
    println!("    activity_r: {:.2}", resistance_data.activity_r);
    println!("    any_r: {:.2}", resistance_data.any_r);
    println!("    majority_r: {:.2}", resistance_data.majority_r);
    println!("-------------------------------------------");


    simulation.run();
}