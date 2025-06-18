// src/simulation/simulation.rs
use crate::simulation::population::{Individual, Population, BACTERIA_LIST, DRUG_SHORT_NAMES};
use crate::rules::apply_rules;
use std::collections::HashMap;
// REMOVED: use rand::seq::SliceRandom; // Not directly used in this file for sampling.

pub struct Simulation {
    pub population: Population,
    pub time_steps: usize,
    pub global_majority_r_proportions: HashMap<(usize, usize), f64>,
    pub majority_r_positive_values_by_combo: HashMap<(usize, usize), Vec<f64>>,
    pub bacteria_indices: HashMap<&'static str, usize>,
    pub drug_indices: HashMap<&'static str, usize>,
}

impl Simulation {
    pub fn new(population_size: usize, time_steps: usize) -> Self {
        println!("--- AMR SIMULATION ---");
        let population = Population::new(population_size);

        // Initialize bacteria_indices and drug_indices
        let mut bacteria_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &bacteria) in BACTERIA_LIST.iter().enumerate() {
            bacteria_indices.insert(bacteria, i);
        }

        let mut drug_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            drug_indices.insert(drug, i);
        }

        // Placeholder initial values for global_majority_r_proportions and majority_r_positive_values_by_combo
        // These would typically be derived from initial population state or external data
        let global_majority_r_proportions = HashMap::new();
        let majority_r_positive_values_by_combo = HashMap::new();

        println!("Initial age of individual 0 AFTER population creation: {} days", population.individuals[0].age);
        println!("--- INITIAL STATE OF INDIVIDUAL 0 (from main.rs) ---");
        println!("  Region Living: {:?}", population.individuals[0].region_living);
        println!("  Region Visiting: {:?}", population.individuals[0].region_visiting); // Should now be Home
        if let Some(&level) = population.individuals[0].level.get("strep_pneu") {
            println!("  strep_pneu: level = {:.2}", level);
        }
        if let Some(&immune_resp) = population.individuals[0].immune_resp.get("strep_pneu") {
            println!("  strep_pneu: immune_resp = {:.2}", immune_resp);
        }
        if let Some(&sepsis) = population.individuals[0].sepsis.get("strep_pneu") {
            println!("  strep_pneu: sepsis = {}", sepsis);
        }
        if let Some(&infectious_syndrome) = population.individuals[0].infectious_syndrome.get("strep_pneu") {
            println!("  strep_pneu: infectious_syndrome = {}", infectious_syndrome);
        }
        if let Some(&date_last_infected) = population.individuals[0].date_last_infected.get("strep_pneu") {
            println!("  strep_pneu: date_last_infected = {}", date_last_infected);
        }
        // MODIFIED: Loop to print all vaccination statuses from the HashMap
        for &bacteria_name in BACTERIA_LIST.iter() {
            if let Some(&status) = population.individuals[0].vaccination_status.get(bacteria_name) {
                println!("  {}_vaccination_status: {}", bacteria_name, status);
            }
        }
        println!("  cur_use_amoxicillin: {}", population.individuals[0].cur_use_drug[drug_indices["amoxicillin"]]);
        println!("  cur_level_amoxicillin: {:.2}", population.individuals[0].cur_level_drug[drug_indices["amoxicillin"]]);
        println!("  current_infection_related_death_risk: {:.2}", population.individuals[0].current_infection_related_death_risk);
        println!("  background_all_cause_mortality_rate: {:.4}", population.individuals[0].background_all_cause_mortality_rate);
        println!("  sexual_contact_level: {:.2}", population.individuals[0].sexual_contact_level);
        println!("  airborne_contact_level_with_adults: {:.2}", population.individuals[0].airborne_contact_level_with_adults);
        println!("  airborne_contact_level_with_children: {:.2}", population.individuals[0].airborne_contact_level_with_children);
        println!("  oral_exposure_level: {:.2}", population.individuals[0].oral_exposure_level);
        println!("  mosquito_exposure_level: {:.2}", population.individuals[0].mosquito_exposure_level);
        println!("  under_care: {}", population.individuals[0].under_care);
        println!("  current_toxicity: {:.2}", population.individuals[0].current_toxicity);
        println!("  mortality_risk_current_toxicity: {:.2}", population.individuals[0].mortality_risk_current_toxicity);

        if let Some(&hospital_acquired) = population.individuals[0].infection_hospital_acquired.get("strep_pneu") {
            println!("  strep_pneu: infection_hospital_acquired = {}", hospital_acquired);
        }
        if let Some(&from_env) = population.individuals[0].cur_infection_from_environment.get("strep_pneu") {
            println!("  strep_pneu: cur_infection_from_environment = {}", from_env);
        }
        let strep_pneu_idx = bacteria_indices["strep_pneu"];
        let amoxicillin_idx = drug_indices["amoxicillin"];
        let resistance_data = &population.individuals[0].resistances[strep_pneu_idx][amoxicillin_idx];
        println!("  strep_pneu resistance to amoxicillin:");
        println!("    microbiome_r: {:.2}", resistance_data.microbiome_r);
        println!("    test_r: {:.2}", resistance_data.test_r);
        println!("    activity_r: {:.2}", resistance_data.activity_r);
        println!("    any_r: {:.2}", resistance_data.any_r);
        println!("    majority_r: {:.2}", resistance_data.majority_r);


       if let Some(&hospital_acquired) = population.individuals[0].infection_hospital_acquired.get("kleb_pneu") {
            println!("  kleb_pneu: infection_hospital_acquired = {}", hospital_acquired);
        }
        if let Some(&from_env) = population.individuals[0].cur_infection_from_environment.get("kleb_pneu") {
            println!("  kleb_pneu: cur_infection_from_environment = {}", from_env);
        }
        let kleb_pneu_idx = bacteria_indices["kleb_pneu"];
        let ceftriaxone_idx = drug_indices["ceftriaxone"];
        let resistance_data = &population.individuals[0].resistances[kleb_pneu_idx][ceftriaxone_idx];
        println!("  kleb_pneu resistance to ceftriaxone:");
        println!("    microbiome_r: {:.2}", resistance_data.microbiome_r);
        println!("    test_r: {:.2}", resistance_data.test_r);
        println!("    activity_r: {:.2}", resistance_data.activity_r);
        println!("    any_r: {:.2}", resistance_data.any_r);
        println!("    majority_r: {:.2}", resistance_data.majority_r);


        Simulation {
            population,
            time_steps,
            global_majority_r_proportions,
            majority_r_positive_values_by_combo,
            bacteria_indices,
            drug_indices,
        }
    }

    pub fn run(&mut self) {
        println!("--- SIMULATION STARTING ---");
        println!("--- AMR SIMULATION (within run function) ---");
        println!("Initial age of individual 0 BEFORE simulation: {} days", self.population.individuals[0].age);
        println!("--- SIMULATION STARTING (within run function) ---");

        // Temporary data collection for global majority_r proportions 
        let mut strep_pneu_amox_majority_r_history: Vec<f64> = Vec::new();
        let mut kleb_pneu_ceftr_majority_r_history: Vec<f64> = Vec::new();

        for t in 0..self.time_steps {
            // Recalculate global majority_r_proportions and positive values at each time step
            // For now, this is a placeholder. In a full model, this would aggregate data
            // from all individuals.
            let mut current_majority_r_positive_values_by_combo: HashMap<(usize, usize), Vec<f64>> = HashMap::new();
            let mut current_infected_counts_with_majority_r: HashMap<(usize, usize), usize> = HashMap::new();
            let mut current_infected_counts_total: HashMap<usize, usize> = HashMap::new();


            for individual in self.population.individuals.iter_mut() {
                apply_rules(
                    individual,
                    t,
                    &self.global_majority_r_proportions, // Pass the global map
                    &self.majority_r_positive_values_by_combo, // Pass the global map
                    &self.bacteria_indices,
                    &self.drug_indices,
                );

                // Update global statistics for the *next* time step
                for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                    // Only count if currently infected
                    if individual.level.get(bacteria_name).map_or(false, |&level| level > 0.001) {
                        *current_infected_counts_total.entry(b_idx).or_insert(0) += 1;

                        for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() { // MODIFIED: Prefixed drug_name with _
                            let resistance_data = &individual.resistances[b_idx][d_idx];
                            if resistance_data.majority_r > 0.0 {
                                current_majority_r_positive_values_by_combo
                                    .entry((b_idx, d_idx))
                                    .or_insert_with(Vec::new)
                                    .push(resistance_data.majority_r);
                                *current_infected_counts_with_majority_r.entry((b_idx, d_idx)).or_insert(0) += 1;
                            }
                        }
                    }
                }
            }
            // Update global_majority_r_proportions based on the current step's data
            // This is a simplified calculation for demonstration.
            let strep_pneu_idx = self.bacteria_indices["strep_pneu"];
            let amoxicillin_idx = self.drug_indices["amoxicillin"];
            let infected_with_strep_pneu = *current_infected_counts_total.get(&strep_pneu_idx).unwrap_or(&0);
            let strep_pneu_amox_majority_r_count = *current_infected_counts_with_majority_r.get(&(strep_pneu_idx, amoxicillin_idx)).unwrap_or(&0);

            let proportion = if infected_with_strep_pneu > 0 {
                strep_pneu_amox_majority_r_count as f64 / infected_with_strep_pneu as f64
            } else {
                0.0
            };
            self.global_majority_r_proportions.insert((strep_pneu_idx, amoxicillin_idx), proportion);
            strep_pneu_amox_majority_r_history.push(proportion);


            // Update global_majority_r_proportions based on the current step's data
            // This is a simplified calculation for demonstration.
            let kleb_pneu_idx = self.bacteria_indices["kleb_pneu"];
            let ceftriaxone_idx = self.drug_indices["ceftriaxone"];
            let infected_with_kleb_pneu = *current_infected_counts_total.get(&kleb_pneu_idx).unwrap_or(&0);
            let kleb_pneu_ceftr_majority_r_count = *current_infected_counts_with_majority_r.get(&(kleb_pneu_idx, ceftriaxone_idx)).unwrap_or(&0);

            let proportion = if infected_with_kleb_pneu > 0 {
                kleb_pneu_ceftr_majority_r_count as f64 / infected_with_kleb_pneu as f64
            } else {
                0.0
            };
            self.global_majority_r_proportions.insert((kleb_pneu_idx, ceftriaxone_idx), proportion);
            kleb_pneu_ceftr_majority_r_history.push(proportion);


            // Log total infected individuals for each bacteria
//          for &bacteria_name in BACTERIA_LIST.iter() {
//              let b_idx = *self.bacteria_indices.get(bacteria_name).unwrap();
//              let infected_count = self.population.individuals.iter()
//                  .filter(|ind| ind.level.get(bacteria_name).map_or(false, |&level| level > 0.001))
//                  .count();
//              println!("Time step {}: Total individuals infected with {} = {}", t, bacteria_name, infected_count);
//          }

//          println!("Time step {}: Global majority_r Proportions (selected):", t);
//          println!("    Strep Pneumonia to Amoxicillin: {:.4}", self.global_majority_r_proportions.get(&(strep_pneu_idx, amoxicillin_idx)).unwrap_or(&0.0));

            // Log individual 0's status
            let individual_0 = &self.population.individuals[0];
            println!("    Time step {}: Individual 0 age = {} days", t, individual_0.age);
/*          if let Some(&level) = individual_0.level.get("strep_pneu") {
                println!("    strep_pneu: level = {:.2}", level);
            }
            if let Some(&immune_resp) = individual_0.immune_resp.get("strep_pneu") {
                println!("    strep_pneu: immune_resp = {:.2}", immune_resp);
            }
            if let Some(&sepsis) = individual_0.sepsis.get("strep_pneu") {
                println!("    strep_pneu: sepsis = {}", sepsis);
            }
            if let Some(&infectious_syndrome) = individual_0.infectious_syndrome.get("strep_pneu") {
                println!("    strep_pneu: infectious_syndrome = {}", infectious_syndrome);
            }
            if let Some(&date_last_infected) = individual_0.date_last_infected.get("strep_pneu") {
                println!("    strep_pneu: date_last_infected = {}", date_last_infected);
            }
            if let Some(&from_env) = individual_0.cur_infection_from_environment.get("strep_pneu") {
                println!("    strep_pneu: cur_infection_from_environment = {}", from_env);
            }
            if let Some(&hospital_acquired) = individual_0.infection_hospital_acquired.get("strep_pneu") {
                println!("    strep_pneu: infection_hospital_acquired = {}", hospital_acquired);
            }
            if let Some(&test_identified) = individual_0.test_identified_infection.get("strep_pneu") {
                println!("    strep_pneu: test_identified_infection = {}", test_identified);
            }
*/


//          println!("    Kleb Pneumonia to ceftriaxone: {:.4}", self.global_majority_r_proportions.get(&(kleb_pneu_idx, ceftriaxone_idx)).unwrap_or(&0.0));

            // Log individual 0's status
            let individual_0 = &self.population.individuals[0];
 /*         if let Some(&level) = individual_0.level.get("kleb_pneu") {
                println!("    kleb_pneu: level = {:.2}", level);
            }
            if let Some(&immune_resp) = individual_0.immune_resp.get("kleb_pneu") {
                println!("    kleb_pneu: immune_resp = {:.2}", immune_resp);
            }
            if let Some(&sepsis) = individual_0.sepsis.get("kleb_pneu") {
                println!("    kleb_pneu: sepsis = {}", sepsis);
            }
            if let Some(&infectious_syndrome) = individual_0.infectious_syndrome.get("kleb_pneu") {
                println!("    kleb_pneu: infectious_syndrome = {}", infectious_syndrome);
            }
            if let Some(&date_last_infected) = individual_0.date_last_infected.get("kleb_pneu") {
                println!("    kleb_pneu: date_last_infected = {}", date_last_infected);
            }
            if let Some(&from_env) = individual_0.cur_infection_from_environment.get("kleb_pneu") {
                println!("    kleb_pneu: cur_infection_from_environment = {}", from_env);
            }
            if let Some(&hospital_acquired) = individual_0.infection_hospital_acquired.get("kleb_pneu") {
                println!("    kleb_pneu: infection_hospital_acquired = {}", hospital_acquired);
            }
            if let Some(&test_identified) = individual_0.test_identified_infection.get("kleb_pneu") {
                println!("    kleb_pneu: test_identified_infection = {}", test_identified);
            }
*/


/*
          // Print Bacteria Infection Status and Level for Individual 0
            println!("        --- Bacteria Infection Status (Individual 0) ---");
            if individual_0.level.is_empty() {
                println!("            Individual 0 has no active bacterial infections.");
            } else {
                for &bacteria_name in BACTERIA_LIST.iter() {
                    if let Some(&level) = individual_0.level.get(bacteria_name) {
                        if level > 0.001 { // Only print if infected (level > a small threshold)
                            println!("            {}: Infected = true, Level = {:.4}", bacteria_name, level);
                        }
                    } else {
                        // Optionally, print if not infected, but to keep output cleaner,
                        // we only show active infections based on the 'if let Some' and 'if level > 0.001'
                        // println!("            {}: Infected = false, Level = 0.0000", bacteria_name);
                    }
                }
            }
*/




// Print Bacteria Infection Status, Level, Immunity, and Active Antibiotics for Individual 0
println!("        --- Bacteria Infection Details (Individual 0) ---");
if individual_0.level.is_empty() {
    println!("            Individual 0 has no active bacterial infections.");
} else {
    for &bacteria_name in BACTERIA_LIST.iter() {
        if let Some(&level) = individual_0.level.get(bacteria_name) {
            // Revert this to your desired threshold, e.g., level > 0.001 or level > 0.0001
            if level > 0.0001 { // Only print if infected (level > a small threshold)
                println!("            Bacteria: {}", bacteria_name);
                println!("              Infected = true");
                println!("              Level = {:.4}", level);
                println!("              Immune Response = {:.4}", individual_0.immune_resp.get(bacteria_name).unwrap_or(&0.0));

                // Print all Antibiotics the person is currently taking
                let mut drugs_being_taken_found = false;
                println!("              Currently Taking Antibiotics (Current Level):");
                for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                    if individual_0.cur_use_drug[drug_idx] {
                        println!("                - {}: Level = {:.4}", drug_name_static, individual_0.cur_level_drug[drug_idx]);
                        drugs_being_taken_found = true;
                    }
                }
                if !drugs_being_taken_found {
                    println!("                None (not currently taking any antibiotics).");
                }

                // Print Antibiotics being used for this bacteria (Existing Logic)
                let b_idx_option = BACTERIA_LIST.iter().position(|&b| b == bacteria_name);
                if let Some(b_idx) = b_idx_option {
                    let mut active_antibiotics_found = false;
                    println!("              Active Antibiotics (Drug Activity against {}):", bacteria_name); // Added clarity here
                    for (drug_idx, &drug_name_static) in DRUG_SHORT_NAMES.iter().enumerate() {
                        if individual_0.cur_use_drug[drug_idx] {
                            // Access activity_r for this specific bacteria-drug combination
                            let activity_r = individual_0.resistances[b_idx][drug_idx].activity_r;
                            if activity_r > 0.0 { // Only show if there's actual activity
                                println!("                - {}: Activity_R = {:.4}", drug_name_static, activity_r);
                                active_antibiotics_found = true;
                            }
                        }
                    }
                    if !active_antibiotics_found {
                        println!("                None (no effective drug activity against this bacteria).");
                    }
                }
                println!(); // Add a blank line for separation between bacteria
            }
        }
    }
}
println!("        --------------------------------------------");






/* 
            println!("    --- Drug Use Status (Individual 0) ---");
            for (drug_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() { // MODIFIED: Prefixed drug_name with _
                println!("      {}: cur_use_drug = {}, cur_level_drug = {:.2}",
                         DRUG_SHORT_NAMES[drug_idx], // MODIFIED: Use DRUG_SHORT_NAMES directly
                         individual_0.cur_use_drug[drug_idx],
                         individual_0.cur_level_drug[drug_idx]);
            }
*/
            println!("    -------------------------------------");
        }

        println!("--- SIMULATION FINISHED (within run function) ---");
/*   
        println!("\n--- Proportion of Strep Pneumonia Infected with majority_r > 0 for Amoxicillin Over Time ---");
        for (t, proportion) in strep_pneu_amox_majority_r_history.iter().enumerate() {
            println!("Time Step {}: {:.4}", t, proportion);
        }
        println!("\n--- Proportion of Kleb Pneumonia Infected with majority_r > 0 for ceftriaxone Over Time ---");
        for (t, proportion) in kleb_pneu_ceftr_majority_r_history.iter().enumerate() {
            println!("Time Step {}: {:.4}", t, proportion);
        }
*/

        println!("--- SIMULATION ENDED ---");
    }
}