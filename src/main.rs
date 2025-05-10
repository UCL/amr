// src/main.rs
mod rules;
mod simulation;

use simulation::population::{Population, BACTERIA_LIST, DRUG_SHORT_NAMES};
use simulation::simulation::run; // Import the run function

fn main() {
    println!("--- AMR SIMULATION ---");

    let mut population = Population::new(30_000);
    let bacteria_to_track = "acinetobac_bau";

    println!(
        "Initial age of individual 0 AFTER population creation: {} days",
        population.individuals[0].age
    );

    // Print the initial state IMMEDIATELY after creating the population
    println!("--- INITIAL STATE OF INDIVIDUAL 0 (from main.rs) ---");
    let ind0 = &population.individuals[0];
    if let Some(level) = ind0.level.get(bacteria_to_track) {
        println!("  {}: level = {:.2}", bacteria_to_track, level);
    }
    if let Some(immune_resp) = ind0.immune_resp.get(bacteria_to_track) {
        println!("  {}: immune_resp = {:.2}", bacteria_to_track, immune_resp);
    }
    if let Some(sepsis) = ind0.sepsis.get(bacteria_to_track) {
        println!("  {}: sepsis = {}", bacteria_to_track, sepsis);
    }
    if let Some(infectious_syndrome) = ind0.infectious_syndrome.get(bacteria_to_track) {
        println!("  {}: infectious_syndrome = {}", bacteria_to_track, infectious_syndrome);
    }
    if let Some(date_last_infected) = ind0.date_last_infected.get(bacteria_to_track) {
        println!("  {}: date_last_infected = {}", bacteria_to_track, date_last_infected);
    }
    println!("  haem_infl_vaccination_status: {}", ind0.haem_infl_vaccination_status);
    println!("  strep_pneu_vaccination_status: {}", ind0.strep_pneu_vaccination_status);
    println!("  salm_typhi_vaccination_status: {}", ind0.salm_typhi_vaccination_status);
    println!("  esch_coli_vaccination_status: {}", ind0.esch_coli_vaccination_status);

    if let Some(gentamicin_index) = DRUG_SHORT_NAMES.iter().position(|&drug| drug == "gentamicin") {
        println!("  cur_use_gentamicin: {}", ind0.cur_use_drug[gentamicin_index]);
        println!("  cur_level_gentamicin: {:.2}", ind0.cur_level_drug[gentamicin_index]);
    } else {
        println!("  gentamicin not found in DRUG_SHORT_NAMES");
    }

    println!("  current_infection_related_death_risk: {:.2}", ind0.current_infection_related_death_risk);
    println!("  background_all_cause_mortality_rate: {:.4}", ind0.background_all_cause_mortality_rate);
    println!("  sexual_contact_level: {:.2}", ind0.sexual_contact_level);
    println!("  airborne_contact_level_with_adults: {:.2}", ind0.airborne_contact_level_with_adults);
    println!("  airborne_contact_level_with_children: {:.2}", ind0.airborne_contact_level_with_children);
    println!("  oral_exposure_level: {:.2}", ind0.oral_exposure_level);
    println!("  mosquito_exposure_level: {:.2}", ind0.mosquito_exposure_level);
    println!("  under_care: {}", ind0.under_care);
    println!("  infection_hospital_acquired: {}", ind0.infection_hospital_acquired);
    println!("  current_toxicity: {:.2}", ind0.current_toxicity);
    println!("  mortality_risk_current_toxicity: {:.2}", ind0.mortality_risk_current_toxicity);

    if let (Some(acinetobac_index), Some(gentamicin_index)) = (
        BACTERIA_LIST.iter().position(|&bacteria| bacteria == "acinetobac_bau"),
        DRUG_SHORT_NAMES.iter().position(|&drug| drug == "gentamicin"),
    ) {
        let resistance = &ind0.resistances[acinetobac_index][gentamicin_index];
        println!("  acinetobac_bau resistance to gentamicin:");
        println!("    microbiome_r: {:.2}", resistance.microbiome_r);
        println!("    test_r: {:.2}", resistance.test_r);
        println!("    activity_r: {:.2}", resistance.activity_r);
        println!("    e_r: {:.2}", resistance.e_r);
        println!("    c_r: {:.2}", resistance.c_r);
    } else {
        println!("  Could not find acinetobac_bau or gentamicin in the lists.");
    }

    println!("--- SIMULATION STARTING ---");

    let num_time_steps = 1;
    run(&mut population, num_time_steps, bacteria_to_track); // Call the run function

    println!("--- SIMULATION ENDED ---");
}