// src/simulation/simulation.rs
use crate::simulation::population::{
    BACTERIA_LIST, DRUG_SHORT_NAMES, Population,
};
use crate::rules::apply_rules;

/// Runs the simulation for a given population and number of time steps.
pub fn run(population: &mut Population, num_time_steps: usize, bacteria_to_track: &str) {
    println!("--- AMR SIMULATION (within run function) ---");
    if let Some(ind) = population.individuals.get(0) {
        println!("Initial age of individual 0 BEFORE simulation: {} days", ind.age);

        // Removed the initial state printing here to avoid duplicates
    }
    println!("--- SIMULATION STARTING (within run function) ---");

    // Find the index of gentamicin
    let gentamicin_index = DRUG_SHORT_NAMES.iter().position(|&drug| drug == "gentamicin");
    let gentamicin_present = gentamicin_index.is_some();
    let gentamicin_idx = gentamicin_index.unwrap_or(0);

    // Find the index of acinetobac_bau
    let acinetobac_index = BACTERIA_LIST.iter().position(|&bacteria| bacteria == "acinetobac_bau");
    let acinetobac_present = acinetobac_index.is_some();
    let acinetobac_idx = acinetobac_index.unwrap_or(0);

    for step in 0..num_time_steps {
        // Apply rules to each individual
        for ind in &mut population.individuals {
            apply_rules(ind, step);

            // Print the values for the specified bacteria for individual 0 AFTER applying rules
            if ind.id == 0 {
                println!("Time step {}: Individual 0 age = {} days", step, ind.age);
                if let Some(level) = ind.level.get(bacteria_to_track) {
                    println!(
                        "  {}: level = {:.2}",
                        bacteria_to_track, level
                    );
                }
                if let Some(immune_resp) = ind.immune_resp.get(bacteria_to_track) {
                    println!(
                        "  {}: immune_resp = {:.2}",
                        bacteria_to_track, immune_resp
                    );
                }
                if let Some(sepsis) = ind.sepsis.get(bacteria_to_track) {
                    println!(
                        "  {}: sepsis = {}",
                        bacteria_to_track, sepsis
                    );
                }
                if let Some(infectious_syndrome) = ind.infectious_syndrome.get(bacteria_to_track) {
                    println!(
                        "  {}: infectious_syndrome = {}",
                        bacteria_to_track, infectious_syndrome
                    );
                }
                if let Some(date_last_infected) = ind.date_last_infected.get(bacteria_to_track) {
                    println!(
                        "  {}: date_last_infected = {}",
                        bacteria_to_track, date_last_infected
                    );
                }
                println!("  haem_infl_vaccination_status: {}", ind.haem_infl_vaccination_status);
                println!("  strep_pneu_vaccination_status: {}", ind.strep_pneu_vaccination_status);
                println!("  salm_typhi_vaccination_status: {}", ind.salm_typhi_vaccination_status);
                println!("  esch_coli_vaccination_status: {}", ind.esch_coli_vaccination_status);

                if gentamicin_present {
                    println!("  cur_use_gentamicin: {}", ind.cur_use_drug[gentamicin_idx]);
                    println!("  cur_level_gentamicin: {:.2}", ind.cur_level_drug[gentamicin_idx]);
                } else {
                    println!("  gentamicin not found in DRUG_SHORT_NAMES");
                }

                println!("  current_infection_related_death_risk: {:.2}", ind.current_infection_related_death_risk);
                println!("  background_all_cause_mortality_rate: {:.4}", ind.background_all_cause_mortality_rate);
                println!("  sexual_contact_level: {:.2}", ind.sexual_contact_level);
                println!("  airborne_contact_level_with_adults: {:.2}", ind.airborne_contact_level_with_adults);
                println!("  airborne_contact_level_with_children: {:.2}", ind.airborne_contact_level_with_children);
                println!("  oral_exposure_level: {:.2}", ind.oral_exposure_level);
                println!("  mosquito_exposure_level: {:.2}", ind.mosquito_exposure_level);
                println!("  under_care: {}", ind.under_care);
                println!("  infection_hospital_acquired: {}", ind.infection_hospital_acquired);
                println!("  current_toxicity: {:.2}", ind.current_toxicity);
                println!("  mortality_risk_current_toxicity: {:.2}", ind.mortality_risk_current_toxicity);

                if acinetobac_present && gentamicin_present {
                    let resistance = &ind.resistances[acinetobac_idx][gentamicin_idx];
                    println!("  acinetobac_bau resistance to gentamicin:");
                    println!("    microbiome_r: {:.2}", resistance.microbiome_r);
                    println!("    test_r: {:.2}", resistance.test_r);
                    println!("    activity_r: {:.2}", resistance.activity_r);
                    println!("    e_r: {:.2}", resistance.e_r);
                    println!("    c_r: {:.2}", resistance.c_r);
                } else {
                    println!("  Could not find acinetobac_bau or gentamicin in the lists.");
                }
            }
        }
    }

    println!("--- SIMULATION FINISHED (within run function) ---");
}