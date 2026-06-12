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
