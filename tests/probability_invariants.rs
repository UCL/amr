use amr_project::config::{get_global_param, parameter_store};
use amr_project::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES};

fn assert_probability(label: impl AsRef<str>, value: f64) {
    let label = label.as_ref();
    assert!(value.is_finite(), "{label} must be finite, got {value}");
    assert!(
        (0.0..=1.0).contains(&value),
        "{label} must be in [0, 1], got {value}"
    );
}

#[test]
fn direct_rng_probability_defaults_are_valid() {
    let store = parameter_store();
    let globals = &store.globals;
    let immunodeficiency = &store.immunodeficiency;

    for (label, value) in [
        (
            "drug_activity_slow_clearance_probability",
            globals.drug_activity_slow_clearance_probability,
        ),
        (
            "double_dose_probability_if_identified_infection",
            globals.double_dose_probability_if_identified_infection,
        ),
        (
            "restart_window_probability",
            globals.restart_window_probability,
        ),
        (
            "random_drug_cessation_probability_if_no_active_infection",
            globals.random_drug_cessation_probability_if_no_active_infection,
        ),
        (
            "antibiotic_infection_prevention_efficacy",
            globals.antibiotic_infection_prevention_efficacy,
        ),
        (
            "microbiome_resistance_transfer_probability_per_day",
            globals.microbiome_resistance_transfer_probability_per_day,
        ),
        (
            "carrier_resistance_inheritance_probability",
            globals.carrier_resistance_inheritance_probability,
        ),
        (
            "majority_r_evolution_rate_per_day_when_drug_present",
            get_global_param("majority_r_evolution_rate_per_day_when_drug_present").unwrap_or(0.0),
        ),
        (
            "test_r_error_probability",
            get_global_param("test_r_error_probability").unwrap_or(0.0),
        ),
        (
            "microbiome_clearance_probability_on_drug_treatment",
            get_global_param("microbiome_clearance_probability_on_drug_treatment").unwrap_or(0.0),
        ),
        (
            "immunosuppression_startup_seed_fraction",
            immunodeficiency.startup_seed_fraction(),
        ),
        (
            "temporary_immunosuppression_recovery_rate_per_day",
            immunodeficiency.temporary_recovery_rate(),
        ),
        (
            "chronic_immunosuppression_recovery_rate_per_day",
            immunodeficiency.chronic_recovery_rate(),
        ),
    ] {
        assert_probability(label, value);
    }

    let temporary_onset = immunodeficiency.temporary_onset_rate();
    let chronic_onset = immunodeficiency.chronic_onset_rate();
    assert_probability(
        "temporary_immunosuppression_onset_rate_per_day",
        temporary_onset,
    );
    assert_probability(
        "chronic_immunosuppression_onset_rate_per_day",
        chronic_onset,
    );
    assert_probability(
        "combined_immunosuppression_onset_rate_per_day",
        temporary_onset + chronic_onset,
    );

    for age_days in [0, 365, 366, 6570, 6571, 23725, 23726] {
        assert_probability(
            format!("chronic_immunodeficiency_probability_at_age_{age_days}"),
            immunodeficiency.chronic_probability(age_days),
        );
    }
}

#[test]
fn per_bacteria_probability_defaults_are_valid() {
    let store = parameter_store();

    for (bacteria_idx, bacteria_name) in BACTERIA_LIST.iter().enumerate() {
        assert_probability(
            format!("{bacteria_name}_drug_cessation_probability"),
            store.bacteria.drug_cessation_probability[bacteria_idx],
        );
        assert_probability(
            format!("{bacteria_name}_microbiome_clearance_probability_per_day"),
            store.bacteria.microbiome_clearance_probability_per_day[bacteria_idx],
        );
        assert_probability(
            format!("{bacteria_name}_treatment_failure_no_second_line_probability"),
            store.bacteria.treatment_failure_no_second_line_probability[bacteria_idx],
        );
    }
}

#[test]
fn drug_probability_defaults_are_valid() {
    let store = parameter_store();

    for (drug_idx, drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        assert_probability(
            format!("{drug_name}_toxicity_death_hazard_per_unit_level"),
            store.drug.toxicity_death_hazard_per_unit_level[drug_idx],
        );
    }
}
