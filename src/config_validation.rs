use std::collections::{HashMap, HashSet};
use std::fmt;

const CONFIG_RS: &str = include_str!("config.rs");

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfigValidationMode {
    Strict,
    Warn,
}

impl ConfigValidationMode {
    pub fn from_env() -> Result<Self, String> {
        Self::from_value(std::env::var("AMR_CONFIG_VALIDATION").ok().as_deref())
    }

    pub fn from_value(value: Option<&str>) -> Result<Self, String> {
        match value.map(str::trim).filter(|value| !value.is_empty()) {
            None => Ok(Self::Strict),
            Some("strict") => Ok(Self::Strict),
            Some("warn") => Ok(Self::Warn),
            Some(value) => Err(format!(
                "AMR_CONFIG_VALIDATION must be 'strict' or 'warn', got '{value}'"
            )),
        }
    }

    pub fn blocks_on_errors(self) -> bool {
        matches!(self, Self::Strict)
    }
}

impl fmt::Display for ConfigValidationMode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Strict => formatter.write_str("strict"),
            Self::Warn => formatter.write_str("warn"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfigValidationSeverity {
    Error,
    Warning,
}

impl fmt::Display for ConfigValidationSeverity {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Error => formatter.write_str("ERROR"),
            Self::Warning => formatter.write_str("WARNING"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConfigValidationIssue {
    severity: ConfigValidationSeverity,
    key: Option<String>,
    message: String,
}

impl ConfigValidationIssue {
    fn error(key: impl Into<String>, message: impl Into<String>) -> Self {
        Self {
            severity: ConfigValidationSeverity::Error,
            key: Some(key.into()),
            message: message.into(),
        }
    }

    fn warning(key: impl Into<String>, message: impl Into<String>) -> Self {
        Self {
            severity: ConfigValidationSeverity::Warning,
            key: Some(key.into()),
            message: message.into(),
        }
    }

    pub fn severity(&self) -> ConfigValidationSeverity {
        self.severity
    }

    pub fn key(&self) -> Option<&str> {
        self.key.as_deref()
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ConfigValidationReport {
    issues: Vec<ConfigValidationIssue>,
}

impl ConfigValidationReport {
    pub fn has_errors(&self) -> bool {
        self.error_count() > 0
    }

    pub fn error_count(&self) -> usize {
        self.issues
            .iter()
            .filter(|issue| issue.severity == ConfigValidationSeverity::Error)
            .count()
    }

    pub fn warning_count(&self) -> usize {
        self.issues
            .iter()
            .filter(|issue| issue.severity == ConfigValidationSeverity::Warning)
            .count()
    }

    pub fn status(&self) -> &'static str {
        if self.has_errors() {
            "failed"
        } else if self.warning_count() > 0 {
            "passed_with_warnings"
        } else {
            "passed"
        }
    }

    pub fn issues(&self) -> &[ConfigValidationIssue] {
        &self.issues
    }

    pub fn render(&self, mode: ConfigValidationMode) -> String {
        let mut output = String::new();
        output.push_str(&format!(
            "[config-validation] {}: {} error(s), {} warning(s), mode={}\n",
            self.status().to_ascii_uppercase(),
            self.error_count(),
            self.warning_count(),
            mode
        ));

        if self.has_errors() && !mode.blocks_on_errors() {
            output.push_str(
                "[config-validation] WARNING: validation errors are not blocking because AMR_CONFIG_VALIDATION=warn\n",
            );
        }

        for issue in &self.issues {
            match issue.key() {
                Some(key) => output.push_str(&format!(
                    "[config-validation] {} key={} {}\n",
                    issue.severity(),
                    key,
                    issue.message()
                )),
                None => output.push_str(&format!(
                    "[config-validation] {} {}\n",
                    issue.severity(),
                    issue.message()
                )),
            }
        }

        output
    }
}

pub fn validate_parameter_map(parameters: &HashMap<String, f64>) -> ConfigValidationReport {
    let mut report = ConfigValidationReport::default();

    validate_literal_config_contracts(parameters, &mut report);
    validate_numeric_values(parameters, &mut report);

    report.issues.sort_by(|left, right| {
        let severity_order = |severity| match severity {
            ConfigValidationSeverity::Error => 0,
            ConfigValidationSeverity::Warning => 1,
        };

        severity_order(left.severity)
            .cmp(&severity_order(right.severity))
            .then_with(|| left.key.cmp(&right.key))
            .then_with(|| left.message.cmp(&right.message))
    });

    report
}

fn validate_literal_config_contracts(
    parameters: &HashMap<String, f64>,
    report: &mut ConfigValidationReport,
) {
    let literal_global_keys = collect_string_after(CONFIG_RS, "get_global_param(\"");
    let literal_required_keys = collect_literal_second_arg(CONFIG_RS, "get_required(", "map");
    let literal_fallback_keys = collect_literal_second_arg(CONFIG_RS, "get_or_default(", "map");

    for key in sorted_missing_keys(parameters, &literal_global_keys) {
        report.issues.push(ConfigValidationIssue::error(
            key,
            "literal get_global_param key is missing from PARAMETERS",
        ));
    }

    for key in sorted_missing_keys(parameters, &literal_required_keys) {
        report.issues.push(ConfigValidationIssue::error(
            key,
            "required parameter used by get_required is missing from PARAMETERS",
        ));
    }

    for key in sorted_missing_keys(parameters, &literal_fallback_keys) {
        report.issues.push(ConfigValidationIssue::error(
            key,
            "literal get_or_default key is missing from PARAMETERS; make the default explicit or use an intentional dynamic fallback",
        ));
    }
}

fn validate_numeric_values(parameters: &HashMap<String, f64>, report: &mut ConfigValidationReport) {
    for (key, value) in parameters {
        if !value.is_finite() {
            report.issues.push(ConfigValidationIssue::error(
                key,
                format!("value must be finite, got {value}"),
            ));
            continue;
        }

        if is_probability_key(key) && !(0.0..=1.0).contains(value) {
            report.issues.push(ConfigValidationIssue::error(
                key,
                format!("probability-style parameter must be in [0, 1], got {value}"),
            ));
        }

        if is_unit_interval_modifier_key(key) && !(0.0..=1.0).contains(value) {
            report.issues.push(ConfigValidationIssue::error(
                key,
                format!("bounded modifier must be in [0, 1], got {value}"),
            ));
        }

        if is_strictly_positive_key(key) && *value <= 0.0 {
            report.issues.push(ConfigValidationIssue::error(
                key,
                format!("strictly positive parameter must be > 0, got {value}"),
            ));
        }

        if is_non_negative_key(key) && *value < 0.0 {
            report.issues.push(ConfigValidationIssue::error(
                key,
                format!("non-negative parameter must be >= 0, got {value}"),
            ));
        }

        if is_boolean_flag_key(key) && !(*value == 0.0 || *value == 1.0) {
            report.issues.push(ConfigValidationIssue::warning(
                key,
                format!("boolean-style flag is expected to be 0.0 or 1.0; runtime treats >0.5 as true, got {value}"),
            ));
        }
    }
}

fn sorted_missing_keys(
    parameters: &HashMap<String, f64>,
    required_keys: &HashSet<String>,
) -> Vec<String> {
    let mut missing = required_keys
        .iter()
        .filter(|key| !parameters.contains_key(*key))
        .cloned()
        .collect::<Vec<_>>();
    missing.sort();
    missing
}

fn is_probability_key(key: &str) -> bool {
    key.starts_with("hgt_prob_")
        || key.ends_with("_probability")
        || key.ends_with("_probability_per_day")
        || key.ends_with("_coverage")
        || key.ends_with("_coverage_target")
        || key.ends_with("_fraction")
        || key.ends_with("_retention")
        || key.ends_with("_efficacy")
        || key.ends_with("_dilution_factor")
}

fn is_boolean_flag_key(key: &str) -> bool {
    key.ends_with("_enabled")
        || key.ends_with("_feature_enabled")
        || key.starts_with("enable_")
        || key.starts_with("debug_")
}

fn is_non_negative_key(key: &str) -> bool {
    if key.contains("log_odds") {
        return false;
    }

    key.ends_with("_days")
        || key.ends_with("_years")
        || key.ends_with("_threshold")
        || key.ends_with("_level")
        || key.ends_with("_half_life_days")
        || key.ends_with("_window_days")
        || key.ends_with("_max_days")
        || key.ends_with("_num_drugs")
        || key.ends_with("_temperature")
        || key.ends_with("_multiplier")
        || key.ends_with("_penalty")
        || key.ends_with("_bonus")
        || key.ends_with("_hazard_per_unit_level")
        || key.ends_with("_emergence_rate")
}

fn is_strictly_positive_key(key: &str) -> bool {
    key.starts_with("drug_") && key.ends_with("_initial_level")
}

fn is_unit_interval_modifier_key(key: &str) -> bool {
    key.ends_with("_penetration")
        || matches!(
            key,
            "resistance_development_inhibition_single_drug"
                | "resistance_development_inhibition_partial_cross"
        )
}

fn skip_ascii_whitespace(source: &str, mut offset: usize) -> usize {
    let bytes = source.as_bytes();
    while offset < bytes.len() && bytes[offset].is_ascii_whitespace() {
        offset += 1;
    }
    offset
}

fn read_quoted_value(source: &str, offset: usize) -> Option<(String, usize)> {
    let bytes = source.as_bytes();
    if bytes.get(offset).copied() != Some(b'"') {
        return None;
    }

    let start = offset + 1;
    let relative_end = source[start..].find('"')?;
    let end = start + relative_end;
    Some((source[start..end].to_string(), end + 1))
}

fn collect_string_after(source: &str, marker: &str) -> HashSet<String> {
    let mut values = HashSet::new();
    let mut offset = 0;

    while let Some(relative_start) = source[offset..].find(marker) {
        let start = offset + relative_start + marker.len();
        let Some(relative_end) = source[start..].find('"') else {
            break;
        };
        values.insert(source[start..start + relative_end].to_string());
        offset = start + relative_end + 1;
    }

    values
}

fn collect_literal_second_arg(
    source: &str,
    marker: &str,
    expected_first_arg: &str,
) -> HashSet<String> {
    let mut values = HashSet::new();
    let mut offset = 0;

    while let Some(relative_start) = source[offset..].find(marker) {
        let mut cursor = offset + relative_start + marker.len();
        cursor = skip_ascii_whitespace(source, cursor);

        if !source[cursor..].starts_with(expected_first_arg) {
            offset = cursor.saturating_add(1).min(source.len());
            continue;
        }
        cursor += expected_first_arg.len();
        cursor = skip_ascii_whitespace(source, cursor);

        if !source[cursor..].starts_with(',') {
            offset = cursor.saturating_add(1).min(source.len());
            continue;
        }
        cursor += 1;
        cursor = skip_ascii_whitespace(source, cursor);

        if let Some((value, end)) = read_quoted_value(source, cursor) {
            values.insert(value);
            offset = end;
        } else {
            offset = cursor.saturating_add(1).min(source.len());
        }
    }

    values
}
