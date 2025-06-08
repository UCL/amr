// src/main.rs

// general thoughts
// add variable for whether msm ?
// infection risk for a specific bacteria / syndrone will depend on sexual_contact_level, airborne_contact_level_with_adults,
// airborne_contact_level_with_children, oral_exposure_level, mosquito_exposure_level
// intend to differentiate infection from another person (dependent on concurrent population)
// from infection from the environment - assume food, water are mainly environment and sex, airborne from another person ?

/* the any_r value represents the level of resistance of resistance (ie the extent to which drug activity against a single 
bacteria which has the resistance is reduced) majority_r is the same but the difference is that any_r represents whether the person 
has ANY resistant bacteria (which will then replicate better than other bacteria in presence of drug) while majority_r indicates 
that all or the majority of bacteria the person has has this resistance level.  so if majority_r has a non-zero value then any_r 
will always have the same value.  If any_r is zero then majority_r has to be zero.  majority_r can be zero while any_r is non-zero - in that 
case majority_r will have a certain probability of taking the value of any_r each day in which drug is present and hence resistant 
virus has the advantage over non resistant  */

// specific next steps - make any_r and majority_r variables taking value 0 - 10
// initiate the value of any_r at the time of infection 

// initial any_r non zero dependent on current level of bacteria
// also need to consider microbiome_r and how it is set and how it influences any_r
// need a variable for world region (e.g. continent) living and perhaps one for world region currently (ie visiting)




mod config;
mod rules;
mod simulation;

use simulation::population::{Population, BACTERIA_LIST, DRUG_SHORT_NAMES};
use simulation::simulation::run; // Import the run function

fn main() {
    println!("--- AMR SIMULATION ---");

    let num_individuals =   300_000;
    let bacteria_to_track = "strep_pneu"; // Define the bacteria to track for initial printout

    let mut population = Population::new(num_individuals);


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

    if let Some(amoxicillin_index) = DRUG_SHORT_NAMES.iter().position(|&drug| drug == "amoxicillin") {
        println!("  cur_use_amoxicillin: {}", ind0.cur_use_drug[amoxicillin_index]);
        println!("  cur_level_amoxicillin: {:.2}", ind0.cur_level_drug[amoxicillin_index]);
    } else {
        println!("  amoxicillin not found in DRUG_SHORT_NAMES");
    }

    println!("  current_infection_related_death_risk: {:.2}", ind0.current_infection_related_death_risk);
    println!("  background_all_cause_mortality_rate: {:.4}", ind0.background_all_cause_mortality_rate);
    println!("  sexual_contact_level: {:.2}", ind0.sexual_contact_level);
    println!("  airborne_contact_level_with_adults: {:.2}", ind0.airborne_contact_level_with_adults);
    println!("  airborne_contact_level_with_children: {:.2}", ind0.airborne_contact_level_with_children);
    println!("  oral_exposure_level: {:.2}", ind0.oral_exposure_level);
    println!("  mosquito_exposure_level: {:.2}", ind0.mosquito_exposure_level);
    println!("  under_care: {}", ind0.under_care);

    // --- CORRECTED LINES FOR HASHMAP PRINTING ---
    // These lines replace the problematic ones.
    if let Some(&hospital_acquired) = ind0.infection_hospital_acquired.get(bacteria_to_track) {
        println!("  {}: infection_hospital_acquired = {}", bacteria_to_track, hospital_acquired);
    } else {
        println!("  {}: infection_hospital_acquired = Not applicable (no active infection)", bacteria_to_track);
    }

    if let Some(&from_env) = ind0.cur_infection_from_environment.get(bacteria_to_track) {
        println!("  {}: cur_infection_from_environment = {}", bacteria_to_track, from_env);
    } else {
        println!("  {}: cur_infection_from_environment = Not applicable (no active infection)", bacteria_to_track);
    }
    // --- END OF CORRECTED HASHMAP PRINTING ---

    println!("  current_toxicity: {:.2}", ind0.current_toxicity);
    println!("  mortality_risk_current_toxicity: {:.2}", ind0.mortality_risk_current_toxicity);

    if let (Some(strep_pneu_index), Some(amoxicillin_index)) = (
        BACTERIA_LIST.iter().position(|&bacteria| bacteria == "strep_pneu"),
        DRUG_SHORT_NAMES.iter().position(|&drug| drug == "amoxicillin"),
    ) {
        let resistance = &ind0.resistances[strep_pneu_index][amoxicillin_index];
        println!("  strep_pneu resistance to amoxicillin:");
        println!("    microbiome_r: {:.2}", resistance.microbiome_r);
        println!("    test_r: {:.2}", resistance.test_r);
        println!("    activity_r: {:.2}", resistance.activity_r);
        println!("    any_r: {:.2}", resistance.any_r);
        println!("    majority_r: {:.2}", resistance.majority_r);
    } else {
        println!("  Could not find strep_pneu or amoxicillin in the lists.");
    }

    println!("--- SIMULATION STARTING ---");

    let num_time_steps = 5 ; // Increased for demonstration if acquisition happens
    run(&mut population, num_time_steps, bacteria_to_track); // Call the run function

    println!("--- SIMULATION ENDED ---");
}