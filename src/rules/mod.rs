// src/rules/rules.rs
use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES};
use rand::Rng;

/// Applies model rules to an individual for one time step.
/// Now increments age by 1 day each step.
pub fn apply_rules(individual: &mut Individual, _time_step: usize) {
    let mut rng = rand::thread_rng();
    let flip_probability = 0.1;

    // Update age: increment by 1 day per time step
    individual.age += 1;

    // Update per-bacteria fields
    for &bacteria in BACTERIA_LIST.iter() {
        if let Some(val) = individual.date_last_infected.get_mut(bacteria) {
            *val += rng.gen_range(-1..=1);
        }
        if let Some(val) = individual.infectious_syndrome.get_mut(bacteria) {
            *val += rng.gen_range(-1.0..=1.0);
        }
        if let Some(val) = individual.level.get_mut(bacteria) {
            *val += rng.gen_range(-1.0..=1.0);
        }
        if let Some(val) = individual.immune_resp.get_mut(bacteria) {
            *val += rng.gen_range(-1.0..=1.0);
        }
        if let Some(val) = individual.sepsis.get_mut(bacteria) {
            if rng.gen::<f64>() < flip_probability {
                *val = !*val;
            }
        }
        if let Some(val) = individual.level_microbiome.get_mut(bacteria) {
            *val += rng.gen_range(-1.0..=1.0);
        }
    }

    // Update vaccination statuses
    if rng.gen::<f64>() < flip_probability {
        individual.haem_infl_vaccination_status = !individual.haem_infl_vaccination_status;
    }
    if rng.gen::<f64>() < flip_probability {
        individual.strep_pneu_vaccination_status = !individual.strep_pneu_vaccination_status;
    }
    if rng.gen::<f64>() < flip_probability {
        individual.salm_typhi_vaccination_status = !individual.salm_typhi_vaccination_status;
    }
    if rng.gen::<f64>() < flip_probability {
        individual.esch_coli_vaccination_status = !individual.esch_coli_vaccination_status;
    }

    // Update drug use and levels
    for i in 0..individual.cur_use_drug.len() {
        individual.cur_use_drug[i] = rng.gen_bool(0.1);
    }
    for i in 0..individual.cur_level_drug.len() {
        individual.cur_level_drug[i] += rng.gen_range(-1.0..=1.0);
    }

    // Update other scalar fields
    individual.current_infection_related_death_risk += rng.gen_range(-1.0..=1.0);
    individual.background_all_cause_mortality_rate += rng.gen_range(-1.0..=1.0);
    individual.sexual_contact_level += rng.gen_range(-1.0..=1.0);
    individual.airborne_contact_level_with_adults += rng.gen_range(-1.0..=1.0);
    individual.airborne_contact_level_with_children += rng.gen_range(-1.0..=1.0);
    individual.oral_exposure_level += rng.gen_range(-1.0..=1.0);
    individual.mosquito_exposure_level += rng.gen_range(-1.0..=1.0);
    if rng.gen::<f64>() < flip_probability {
        individual.under_care = !individual.under_care;
    }
    if rng.gen::<f64>() < flip_probability {
        individual.infection_hospital_acquired = !individual.infection_hospital_acquired;
    }
    individual.current_toxicity += rng.gen_range(-1.0..=1.0);
    individual.mortality_risk_current_toxicity += rng.gen_range(-1.0..=1.0);

    // Update resistances
    for i in 0..BACTERIA_LIST.len() {
        for j in 0..DRUG_SHORT_NAMES.len() {
            individual.resistances[i][j].microbiome_r += rng.gen_range(-1.0..=1.0);
            individual.resistances[i][j].test_r += rng.gen_range(-1.0..=1.0);
            individual.resistances[i][j].activity_r += rng.gen_range(-1.0..=1.0);
            individual.resistances[i][j].e_r += rng.gen_range(-1.0..=1.0);
            individual.resistances[i][j].c_r += rng.gen_range(-1.0..=1.0);
        }
    }
}