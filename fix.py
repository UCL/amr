import sys

with open('src/simulation/simulation.rs', 'r', encoding='utf-8') as f:
    text = f.read()

old_str = """                newly_infected_any_r_hospital_by_bacteria: {
                    let mut counts = vec![0; BACTERIA_LIST.len()];
                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() { continue; }
                        for b_idx in 0..BACTERIA_LIST.len() {
                            if individual.date_last_infected_keep[b_idx] == t as i32 && !individual.resistance_mechanisms[b_idx].is_empty() {
                                if individual.infection_hospital_acquired[b_idx] {
                                    counts[b_idx] += 1;
                                }
                            }
                        }
                    }
                    counts
                },
                newly_infected_any_r_community_by_bacteria: {
                    let mut counts = vec![0; BACTERIA_LIST.len()];
                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() { continue; }
                        for b_idx in 0..BACTERIA_LIST.len() {
                            if individual.date_last_infected_keep[b_idx] == t as i32 && !individual.resistance_mechanisms[b_idx].is_empty() {
                                if !individual.infection_hospital_acquired[b_idx] {
                                    counts[b_idx] += 1;
                                }
                            }
                        }
                    }
                    counts
                },"""

new_str = """                newly_infected_any_r_hospital_by_bacteria,
                newly_infected_any_r_community_by_bacteria,"""

old_insert = """            let summary = TimeStepSummary {"""
new_insert = """            let mut newly_infected_any_r_hospital_by_bacteria = vec![0; BACTERIA_LIST.len()];
            let mut newly_infected_any_r_community_by_bacteria = vec![0; BACTERIA_LIST.len()];
            for individual in &self.population.individuals {
                if individual.date_of_death.is_some() { continue; }
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.date_last_infected_keep[b_idx] == t as i32 {
                        if individual.resistance_mechanisms[b_idx].iter().any(|&m| m) {
                            if individual.infection_hospital_acquired[b_idx] {
                                newly_infected_any_r_hospital_by_bacteria[b_idx] += 1;
                            } else {
                                newly_infected_any_r_community_by_bacteria[b_idx] += 1;
                            }
                        }
                    }
                }
            }

            let summary = TimeStepSummary {"""

if old_str in text and old_insert in text:
    text = text.replace(old_str, new_str)
    text = text.replace(old_insert, new_insert)
    with open('src/simulation/simulation.rs', 'w', encoding='utf-8') as f:
        f.write(text)
    print("Success")
else:
    print("Failed to find")
