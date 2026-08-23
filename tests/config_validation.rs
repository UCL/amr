use amr_project::config::PARAMETERS;
use amr_project::config_validation::{
    validate_parameter_map, ConfigValidationMode, ConfigValidationSeverity,
};

fn assert_has_error(report: &amr_project::config_validation::ConfigValidationReport, key: &str) {
    assert!(
        report.issues().iter().any(|issue| {
            issue.severity() == ConfigValidationSeverity::Error && issue.key() == Some(key)
        }),
        "expected validation error for {key}, got:\n{}",
        report.render(ConfigValidationMode::Strict)
    );
}

#[test]
fn canonical_parameters_pass_runtime_validation() {
    let report = validate_parameter_map(&PARAMETERS);

    assert!(
        !report.has_errors(),
        "canonical PARAMETERS should pass runtime validation:\n{}",
        report.render(ConfigValidationMode::Strict)
    );
}

#[test]
fn missing_literal_fallback_key_is_a_hard_error() {
    let mut parameters = PARAMETERS.clone();
    parameters.remove("toxicity_age_multiplier_elderly");

    let report = validate_parameter_map(&parameters);

    assert_has_error(&report, "toxicity_age_multiplier_elderly");
}

#[test]
fn probability_style_key_outside_unit_interval_is_a_hard_error() {
    let mut parameters = PARAMETERS.clone();
    parameters.insert("random_drug_cessation_probability".to_string(), 1.2);

    let report = validate_parameter_map(&parameters);

    assert_has_error(&report, "random_drug_cessation_probability");
}

#[test]
fn non_finite_value_is_a_hard_error() {
    let mut parameters = PARAMETERS.clone();
    parameters.insert("drug_selection_temperature".to_string(), f64::NAN);

    let report = validate_parameter_map(&parameters);

    assert_has_error(&report, "drug_selection_temperature");
}

#[test]
fn emergence_coefficient_is_unbounded_above_but_must_not_be_negative() {
    let key = "bacteria_escherichia_coli_mechanism_mutation_gyra_primary_emergence_rate";
    let mut high_parameters = PARAMETERS.clone();
    high_parameters.insert(key.to_string(), 30.0);
    let high_report = validate_parameter_map(&high_parameters);
    assert!(
        !high_report.has_errors(),
        "unbounded emergence coefficient should accept values above one:\n{}",
        high_report.render(ConfigValidationMode::Strict)
    );

    let mut negative_parameters = PARAMETERS.clone();
    negative_parameters.insert(key.to_string(), -0.1);
    let negative_report = validate_parameter_map(&negative_parameters);
    assert_has_error(&negative_report, key);
}

#[test]
fn drug_initial_level_must_be_strictly_positive() {
    let key = "drug_amoxicillin_initial_level";
    let mut parameters = PARAMETERS.clone();
    parameters.insert(key.to_string(), 0.0);

    let report = validate_parameter_map(&parameters);

    assert_has_error(&report, key);
}

#[test]
fn penetration_and_emergence_inhibition_modifiers_must_be_in_unit_interval() {
    for (key, value) in [
        ("syndrome_3_drug_daptomycin_penetration", 1.1),
        ("resistance_development_inhibition_single_drug", -0.1),
        ("resistance_development_inhibition_partial_cross", 1.1),
    ] {
        let mut parameters = PARAMETERS.clone();
        parameters.insert(key.to_string(), value);

        let report = validate_parameter_map(&parameters);

        assert_has_error(&report, key);
    }
}

#[test]
fn validation_mode_is_strict_by_default_and_warn_when_explicit() {
    assert_eq!(
        ConfigValidationMode::from_value(None),
        Ok(ConfigValidationMode::Strict)
    );
    assert_eq!(
        ConfigValidationMode::from_value(Some("warn")),
        Ok(ConfigValidationMode::Warn)
    );
    assert!(ConfigValidationMode::from_value(Some("off")).is_err());
}
