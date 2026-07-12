use amr_project::simulation::population::{
    drug_class_for_drug, DrugClass, Individual, InfectionResolutionType, Population, Region,
    ResistanceMechanism, BACTERIA_CARRIAGE_COMPARTMENTS, BACTERIA_COUNT, BACTERIA_GROUPS,
    BACTERIA_LIST, DRUG_CLASS_LOOKUP, DRUG_SHORT_NAMES, MICROBIOME_RESISTANCE_LEVEL_COUNT,
    TRACK_RESISTANCE_ACQUISITION_PROVENANCE,
};
use amr_project::simulation::simulation::{
    CalibrationMode, Simulation, DIAGNOSTIC_CASCADE_SETTING_COUNT, DIAGNOSTIC_CASCADE_STAGE_COUNT,
    RESISTANCE_MECHANISM_FAMILY_COUNT,
};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use std::collections::HashSet;

// These invariants are guardrails, not a freeze on model development. If the
// model intentionally adds/removes bacteria, drugs, regions, mechanisms, or
// summary dimensions, update these expectations in the same PR as the model
// change so index-based arrays and exported summaries move coherently.
const MODEL_REGION_COUNT: usize = 6;
const AGE_GROUP_COUNT: usize = 5;
const DEATH_CAUSE_COUNT: usize = 4;
const SYNDROME_COUNT: usize = 10;
const CARRIAGE_DURATION_BIN_COUNT: usize = 5;
const DRUG_COUNT_HISTOGRAM_LEN: usize = 4;

fn assert_len(name: &str, actual: usize, expected: usize) {
    assert_eq!(actual, expected, "{name} length should remain stable");
}

fn assert_unique(name: &str, values: &[&str]) {
    let mut seen = HashSet::new();
    for value in values {
        assert!(
            seen.insert(*value),
            "{name} contains duplicate value {value}"
        );
    }
}

fn assert_individual_dimensions(individual: &Individual) {
    let bacteria = BACTERIA_COUNT;
    let drugs = DRUG_SHORT_NAMES.len();
    let resolution_types = InfectionResolutionType::all().len();

    assert_len(
        "date_last_infected",
        individual.date_last_infected.len(),
        bacteria,
    );
    assert_len(
        "date_last_infected_keep",
        individual.date_last_infected_keep.len(),
        bacteria,
    );
    assert_len(
        "infectious_syndrome",
        individual.infectious_syndrome.len(),
        bacteria,
    );
    assert_len("level", individual.level.len(), bacteria);
    assert_len(
        "predicted_infection_risk",
        individual.predicted_infection_risk.len(),
        bacteria,
    );
    assert_len(
        "clearance_hazard",
        individual.clearance_hazard.len(),
        bacteria,
    );
    assert_len(
        "clearance_ready_day",
        individual.clearance_ready_day.len(),
        bacteria,
    );
    assert_len("sepsis", individual.sepsis.len(), bacteria);
    assert_len(
        "sepsis_onset_day",
        individual.sepsis_onset_day.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_open",
        individual.sepsis_episode_open.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_context_at_onset",
        individual.sepsis_episode_context_at_onset.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_best_activity_at_onset",
        individual.sepsis_episode_best_activity_at_onset.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_effective_at_onset",
        individual.sepsis_episode_effective_at_onset.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_first_effective_day",
        individual.sepsis_episode_first_effective_day.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_delay_bucket_recorded",
        individual.sepsis_episode_delay_bucket_recorded.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_region_at_onset",
        individual.sepsis_episode_region_at_onset.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_hospitalized_at_onset",
        individual.sepsis_episode_hospitalized_at_onset.len(),
        bacteria,
    );
    assert_len(
        "sepsis_episode_age_group_at_onset",
        individual.sepsis_episode_age_group_at_onset.len(),
        bacteria,
    );
    assert_len(
        "diagnostic_cascade_open",
        individual.diagnostic_cascade_open.len(),
        bacteria,
    );
    assert_len(
        "diagnostic_cascade_entry_time_step",
        individual.diagnostic_cascade_entry_time_step.len(),
        bacteria,
    );
    assert_len(
        "diagnostic_cascade_entry_hospitalized",
        individual.diagnostic_cascade_entry_hospitalized.len(),
        bacteria,
    );
    assert_len(
        "diagnostic_cascade_bacterial_identification_recorded",
        individual
            .diagnostic_cascade_bacterial_identification_recorded
            .len(),
        bacteria,
    );
    assert_len(
        "diagnostic_cascade_resistance_testing_recorded",
        individual
            .diagnostic_cascade_resistance_testing_recorded
            .len(),
        bacteria,
    );
    assert_len(
        "diagnostic_cascade_targeted_treatment_recorded",
        individual
            .diagnostic_cascade_targeted_treatment_recorded
            .len(),
        bacteria,
    );
    assert_len(
        "diagnostic_cascade_effective_targeted_treatment_recorded",
        individual
            .diagnostic_cascade_effective_targeted_treatment_recorded
            .len(),
        bacteria,
    );
    assert_len(
        "infection_prevented_by_drug",
        individual.infection_prevented_by_drug.len(),
        bacteria,
    );
    assert_len(
        "presence_microbiome",
        individual.presence_microbiome.len(),
        bacteria,
    );
    assert_len(
        "date_microbiome_acquired",
        individual.date_microbiome_acquired.len(),
        bacteria,
    );
    assert_len(
        "microbiome_acquired_today",
        individual.microbiome_acquired_today.len(),
        bacteria,
    );
    assert_len(
        "microbiome_acquired_on_drug_today",
        individual.microbiome_acquired_on_drug_today.len(),
        bacteria,
    );
    assert_len(
        "microbiome_cleared_today",
        individual.microbiome_cleared_today.len(),
        bacteria,
    );
    assert_len(
        "cleared_any_r_microbiome_categories",
        individual.cleared_any_r_microbiome_categories.len(),
        bacteria,
    );
    assert_len(
        "vaccination_status",
        individual.vaccination_status.len(),
        bacteria,
    );
    assert_len(
        "infection_has_caused_symptoms",
        individual.infection_has_caused_symptoms.len(),
        bacteria,
    );
    assert_len(
        "test_identified_infection",
        individual.test_identified_infection.len(),
        bacteria,
    );
    assert_len(
        "test_for_resistance",
        individual.test_for_resistance.len(),
        bacteria,
    );
    assert_len(
        "resistance_test_initiated_day",
        individual.resistance_test_initiated_day.len(),
        bacteria,
    );
    assert_len(
        "infection_hospital_acquired",
        individual.infection_hospital_acquired.len(),
        bacteria,
    );
    assert_len("mechanism_any", individual.mechanism_any.len(), bacteria);
    assert_len(
        "mechanism_majority",
        individual.mechanism_majority.len(),
        bacteria,
    );
    assert_len(
        "mechanism_microbiome",
        individual.mechanism_microbiome.len(),
        bacteria,
    );
    assert_len(
        "infection_resolution_this_timestep",
        individual.infection_resolution_this_timestep.len(),
        bacteria,
    );
    for (idx, row) in individual
        .infection_resolution_this_timestep
        .iter()
        .enumerate()
    {
        assert_len(
            &format!("infection_resolution_this_timestep[{idx}]"),
            row.len(),
            resolution_types,
        );
    }
    assert_len(
        "day_7_since_last_infection_drug_used",
        individual.day_7_since_last_infection_drug_used.len(),
        bacteria,
    );
    assert_len(
        "bacteria_level_at_drug_start",
        individual.bacteria_level_at_drug_start.len(),
        bacteria,
    );
    assert_len(
        "days_on_current_treatment",
        individual.days_on_current_treatment.len(),
        bacteria,
    );
    assert_len(
        "treatment_failure_assessed",
        individual.treatment_failure_assessed.len(),
        bacteria,
    );
    assert_len(
        "drug_activity_response_multiplier",
        individual.drug_activity_response_multiplier.len(),
        bacteria,
    );
    assert_len(
        "drug_stopped_with_infection_day",
        individual.drug_stopped_with_infection_day.len(),
        bacteria,
    );
    assert_len(
        "bacteria_level_at_drug_cessation",
        individual.bacteria_level_at_drug_cessation.len(),
        bacteria,
    );
    assert_len(
        "stopped_drug_index",
        individual.stopped_drug_index.len(),
        bacteria,
    );
    assert_len(
        "restart_window_assessed",
        individual.restart_window_assessed.len(),
        bacteria,
    );
    assert_len(
        "date_last_drug_failure",
        individual.date_last_drug_failure.len(),
        bacteria,
    );

    assert_len("cur_use_drug", individual.cur_use_drug.len(), drugs);
    assert_len("drug_use_context", individual.drug_use_context.len(), drugs);
    assert_len("cur_level_drug", individual.cur_level_drug.len(), drugs);
    assert_len(
        "date_drug_initiated",
        individual.date_drug_initiated.len(),
        drugs,
    );
    assert_len(
        "date_drug_initiated_keep",
        individual.date_drug_initiated_keep.len(),
        drugs,
    );
    assert_len("ever_taken_drug", individual.ever_taken_drug.len(), drugs);
    assert_len(
        "drug_toxicity_reservoir",
        individual.drug_toxicity_reservoir.len(),
        drugs,
    );
    assert_len(
        "toxicity_stopped_drug_day",
        individual.toxicity_stopped_drug_day.len(),
        drugs,
    );
    assert_len(
        "drug_score_on_selection_day",
        individual.drug_score_on_selection_day.len(),
        drugs,
    );

    assert_len("resistances", individual.resistances.len(), bacteria);
    for (idx, row) in individual.resistances.iter().enumerate() {
        assert_len(&format!("resistances[{idx}]"), row.len(), drugs);
    }

    if TRACK_RESISTANCE_ACQUISITION_PROVENANCE {
        assert_len(
            "how_resistance_acquired",
            individual.how_resistance_acquired.len(),
            bacteria,
        );
        for (idx, row) in individual.how_resistance_acquired.iter().enumerate() {
            assert_len(&format!("how_resistance_acquired[{idx}]"), row.len(), drugs);
        }
    } else {
        assert!(
            individual.how_resistance_acquired.is_empty(),
            "provenance matrix should stay unallocated while tracking is disabled"
        );
    }
}

#[test]
fn static_model_dimensions_are_unique_and_consistent() {
    assert_len("BACTERIA_GROUPS", BACTERIA_GROUPS.len(), BACTERIA_COUNT);
    assert_len(
        "BACTERIA_CARRIAGE_COMPARTMENTS",
        BACTERIA_CARRIAGE_COMPARTMENTS.len(),
        BACTERIA_COUNT,
    );
    assert_unique("BACTERIA_LIST", &BACTERIA_LIST);
    assert_unique("DRUG_SHORT_NAMES", DRUG_SHORT_NAMES);

    assert_len(
        "DrugClass::all",
        DrugClass::all().len(),
        DrugClass::NUM_CLASSES,
    );
    let drug_class_names: Vec<&str> = DrugClass::all().iter().map(DrugClass::as_str).collect();
    assert_unique("DrugClass::all", &drug_class_names);
    for (expected_idx, drug_class) in DrugClass::all().iter().enumerate() {
        assert_eq!(
            drug_class.index(),
            expected_idx,
            "DrugClass::all order should match enum indices"
        );
    }

    assert_len(
        "DRUG_CLASS_LOOKUP",
        DRUG_CLASS_LOOKUP.len(),
        DRUG_SHORT_NAMES.len(),
    );
    for (drug_idx, class_idx) in DRUG_CLASS_LOOKUP.iter().enumerate() {
        assert!(
            *class_idx < DrugClass::NUM_CLASSES,
            "DRUG_CLASS_LOOKUP[{drug_idx}] contains out-of-range drug class index {class_idx}"
        );
        assert_eq!(
            *class_idx,
            drug_class_for_drug(drug_idx).index(),
            "DRUG_CLASS_LOOKUP[{drug_idx}] should match drug_class_for_drug"
        );
    }

    let mechanism_names: Vec<&str> = ResistanceMechanism::all()
        .iter()
        .map(ResistanceMechanism::as_str)
        .collect();
    assert_unique("ResistanceMechanism::all", &mechanism_names);
    assert!(
        ResistanceMechanism::all().len() <= u64::BITS as usize,
        "mechanism bitmasks use u64 storage, so mechanism count must not exceed 64"
    );

    let model_regions = [
        Region::NorthAmerica,
        Region::SouthAmerica,
        Region::Africa,
        Region::Asia,
        Region::Europe,
        Region::Oceania,
    ];
    assert_len(
        "model regions excluding Home",
        model_regions.len(),
        MODEL_REGION_COUNT,
    );
}

#[test]
fn individual_and_population_constructors_preserve_core_dimensions() {
    let mut rng = SmallRng::seed_from_u64(987_654_321);
    let individual = Individual::new(0, 40 * 365, "male".to_string(), &mut rng);
    assert_individual_dimensions(&individual);

    let population = Population::new(8, &mut rng);
    assert_len("population", population.individuals.len(), 8);
    for individual in &population.individuals {
        assert_individual_dimensions(individual);
    }
}

#[test]
fn simulation_constructor_preserves_lookup_and_flat_matrix_dimensions() {
    let simulation = Simulation::new(8, 1, false, Some(135_791_113), CalibrationMode::Partial);
    let bacteria = BACTERIA_COUNT;
    let drugs = DRUG_SHORT_NAMES.len();
    let bacteria_drug = bacteria * drugs;

    assert_len(
        "bacteria_indices",
        simulation.bacteria_indices.len(),
        bacteria,
    );
    for (expected_idx, bacteria_name) in BACTERIA_LIST.iter().enumerate() {
        assert_eq!(
            simulation.bacteria_indices.get(bacteria_name),
            Some(&expected_idx),
            "bacteria index should match BACTERIA_LIST order"
        );
    }

    assert_len("drug_indices", simulation.drug_indices.len(), drugs);
    for (expected_idx, drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        assert_eq!(
            simulation.drug_indices.get(drug_name),
            Some(&expected_idx),
            "drug index should match DRUG_SHORT_NAMES order"
        );
    }

    assert_len(
        "potency_matrix",
        simulation.potency_matrix.len(),
        bacteria_drug,
    );
    assert_len(
        "mic_lt2_majority_r_thresholds",
        simulation.mic_lt2_majority_r_thresholds.len(),
        bacteria_drug,
    );
    assert_len(
        "relevant_drugs_by_bacteria",
        simulation.relevant_drugs_by_bacteria.len(),
        bacteria,
    );
    for (bacteria_idx, relevant_drugs) in simulation.relevant_drugs_by_bacteria.iter().enumerate() {
        for drug_idx in relevant_drugs {
            assert!(
                *drug_idx < drugs,
                "relevant_drugs_by_bacteria[{bacteria_idx}] contains out-of-range drug index {drug_idx}"
            );
        }
    }

    assert_len(
        "cross_resistance_groups",
        simulation.cross_resistance_groups.len(),
        bacteria,
    );
    for (bacteria_idx, groups) in simulation.cross_resistance_groups.iter().enumerate() {
        for group in groups {
            for drug_idx in group {
                assert!(
                    *drug_idx < drugs,
                    "cross_resistance_groups[{bacteria_idx}] contains out-of-range drug index {drug_idx}"
                );
            }
        }
    }
}

#[test]
fn full_summary_rows_preserve_expected_vector_dimensions() {
    let mut simulation = Simulation::new(8, 1, false, Some(112_358_132), CalibrationMode::None);
    simulation.run();

    let summary = simulation
        .summary_log
        .first()
        .expect("one timestep should produce one summary row");
    let bacteria = BACTERIA_COUNT;
    let drugs = DRUG_SHORT_NAMES.len();
    let mechanisms = ResistanceMechanism::all().len();
    let bacteria_drug = bacteria * drugs;
    let bacteria_region = bacteria * MODEL_REGION_COUNT;
    let bacteria_syndrome = bacteria * SYNDROME_COUNT;

    assert_len(
        "infected_and_on_any_drug_by_bacteria",
        summary.infected_and_on_any_drug_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "infections_by_bacteria",
        summary.infections_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "deaths_by_bacteria",
        summary.deaths_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "new_active_infections_by_bacteria",
        summary.new_active_infections_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "active_infection_days_by_bacteria",
        summary.active_infection_days_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "treated_infection_days_by_bacteria",
        summary.treated_infection_days_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "effective_treated_infection_days_by_bacteria",
        summary.effective_treated_infection_days_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "infection_resolution_count_by_bacteria",
        summary.infection_resolution_count_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "sepsis_onset_count_by_bacteria",
        summary.sepsis_onset_count_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "infection_death_count_by_bacteria",
        summary.infection_death_count_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "drug_failure_count_by_bacteria",
        summary.drug_failure_count_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "carrier_at_risk_person_days_by_bacteria",
        summary.carrier_at_risk_person_days_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "non_carrier_at_risk_person_days_by_bacteria",
        summary.non_carrier_at_risk_person_days_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "new_infections_in_carriers_by_bacteria",
        summary.new_infections_in_carriers_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "new_infections_in_non_carriers_by_bacteria",
        summary.new_infections_in_non_carriers_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "new_any_r_infections_in_carriers_by_bacteria",
        summary.new_any_r_infections_in_carriers_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "new_any_r_infections_in_non_carriers_by_bacteria",
        summary
            .new_any_r_infections_in_non_carriers_by_bacteria
            .len(),
        bacteria,
    );
    assert_len(
        "presence_microbiome_by_bacteria",
        summary.presence_microbiome_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "diagnostic_cascade_stage_counts",
        summary.diagnostic_cascade_stage_counts.len(),
        DIAGNOSTIC_CASCADE_STAGE_COUNT,
    );
    assert_len(
        "diagnostic_cascade_stage_counts_by_setting",
        summary.diagnostic_cascade_stage_counts_by_setting.len(),
        DIAGNOSTIC_CASCADE_STAGE_COUNT * DIAGNOSTIC_CASCADE_SETTING_COUNT,
    );
    assert_len(
        "infected_with_test_identified_by_bacteria",
        summary.infected_with_test_identified_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "day_7_evaluations_by_bacteria",
        summary.day_7_evaluations_by_bacteria.len(),
        bacteria,
    );
    assert_len(
        "drug_selection_count_by_bacteria",
        summary.drug_selection_count_by_bacteria.len(),
        bacteria,
    );

    assert_len(
        "resistance_by_bacteria_drug",
        summary.resistance_by_bacteria_drug.len(),
        bacteria_drug,
    );
    assert_len(
        "currently_on_drug_by_bacteria_drug",
        summary.currently_on_drug_by_bacteria_drug.len(),
        bacteria_drug,
    );
    assert_len(
        "any_r_sum_by_bacteria_drug",
        summary.any_r_sum_by_bacteria_drug.len(),
        bacteria_drug,
    );
    assert_len(
        "drug_score_sums_by_bacteria_drug",
        summary.drug_score_sums_by_bacteria_drug.len(),
        bacteria_drug,
    );

    assert_len(
        "newly_infected_by_bacteria_region",
        summary.newly_infected_by_bacteria_region.len(),
        bacteria_region,
    );
    assert_len(
        "deaths_infected_by_bacteria_region",
        summary.deaths_infected_by_bacteria_region.len(),
        bacteria_region,
    );
    assert_len(
        "drug_failure_events_by_bacteria_region",
        summary.drug_failure_events_by_bacteria_region.len(),
        bacteria_region,
    );
    assert_len(
        "presence_microbiome_by_bacteria_by_region",
        summary.presence_microbiome_by_bacteria_by_region.len(),
        bacteria_region,
    );

    assert_len(
        "cleared_any_r_microbiome_categories",
        summary.cleared_any_r_microbiome_categories.len(),
        bacteria * MICROBIOME_RESISTANCE_LEVEL_COUNT,
    );
    assert_len(
        "carriage_duration_bins_by_bacteria",
        summary.carriage_duration_bins_by_bacteria.len(),
        bacteria * CARRIAGE_DURATION_BIN_COUNT,
    );
    assert_len(
        "infected_with_bacteria_and_mechanism",
        summary.infected_with_bacteria_and_mechanism.len(),
        bacteria * mechanisms,
    );
    assert_len(
        "infection_days_with_any_resistance_mechanism_by_bacteria",
        summary
            .infection_days_with_any_resistance_mechanism_by_bacteria
            .len(),
        bacteria,
    );
    assert_len(
        "infection_days_with_resistance_mechanism_family_by_bacteria",
        summary
            .infection_days_with_resistance_mechanism_family_by_bacteria
            .len(),
        bacteria * RESISTANCE_MECHANISM_FAMILY_COUNT,
    );

    assert_len(
        "currently_on_drug_by_drug",
        summary.currently_on_drug_by_drug.len(),
        drugs,
    );
    assert_len(
        "currently_on_drug_by_region_drug",
        summary.currently_on_drug_by_region_drug.len(),
        MODEL_REGION_COUNT * drugs,
    );
    assert_len(
        "any_r_sum_by_region",
        summary.any_r_sum_by_region.len(),
        MODEL_REGION_COUNT,
    );
    assert_len(
        "infected_count_by_region",
        summary.infected_count_by_region.len(),
        MODEL_REGION_COUNT,
    );
    assert_len(
        "living_population_by_region",
        summary.living_population_by_region.len(),
        MODEL_REGION_COUNT,
    );
    assert_len(
        "hospital_population_by_region",
        summary.hospital_population_by_region.len(),
        MODEL_REGION_COUNT,
    );
    assert_len(
        "age_distribution_by_region",
        summary.age_distribution_by_region.len(),
        MODEL_REGION_COUNT * AGE_GROUP_COUNT,
    );
    assert_len(
        "deaths_by_region",
        summary.deaths_by_region.len(),
        MODEL_REGION_COUNT * DEATH_CAUSE_COUNT,
    );
    assert_len(
        "deaths_by_region_age",
        summary.deaths_by_region_age.len(),
        MODEL_REGION_COUNT * AGE_GROUP_COUNT * DEATH_CAUSE_COUNT,
    );

    assert_len(
        "infected_by_syndrome",
        summary.infected_by_syndrome.len(),
        SYNDROME_COUNT,
    );
    assert_len(
        "infected_by_syndrome_by_bacteria",
        summary.infected_by_syndrome_by_bacteria.len(),
        bacteria_syndrome,
    );
    assert_len(
        "syndrome_population_by_region",
        summary.syndrome_population_by_region.len(),
        SYNDROME_COUNT * MODEL_REGION_COUNT,
    );
    assert_len(
        "syndrome_deaths_sepsis_by_region",
        summary.syndrome_deaths_sepsis_by_region.len(),
        SYNDROME_COUNT * MODEL_REGION_COUNT,
    );
    assert_len(
        "people_by_drug_count",
        summary.people_by_drug_count.len(),
        DRUG_COUNT_HISTOGRAM_LEN,
    );
}
