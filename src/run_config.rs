//! Typed runtime configuration for one AMR model execution.
//!
//! The model keeps useful standalone defaults, but a launcher may supply an
//! `amr-run-config/v1` document through [`AMR_RUN_CONFIG_ENV`]. Resource policy
//! such as CPU pinning, NUMA placement, and memory limits is intentionally not
//! part of this model-owned contract.

use crate::config_validation::ConfigValidationMode;
use crate::simulation::simulation::CalibrationMode;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};

pub const AMR_RUN_CONFIG_ENV: &str = "AMR_RUN_CONFIG";

const SCHEMA_VERSION: &str = "amr-run-config/v1";
const CONTRACT: &str = "typed_adapter_runtime_config";
const OVERRIDE_SEMANTICS: &str = "explicit_keys_only_model_defaults_when_omitted";
const DEFAULT_POPULATION_SIZE: usize = 3_000_000;
const DEFAULT_MODEL_OUTPUT_DIR: &str = "amr_simulation_output_analysis_outputs";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunOutputPaths {
    pub run_root: PathBuf,
    pub report_json: PathBuf,
    pub progress_json: PathBuf,
    pub progress_jsonl: PathBuf,
    pub model_output_dir: PathBuf,
}

impl Default for RunOutputPaths {
    fn default() -> Self {
        Self {
            run_root: PathBuf::from("."),
            report_json: PathBuf::from("amr_report.json"),
            progress_json: PathBuf::from("amr_progress.json"),
            progress_jsonl: PathBuf::from("amr_progress.jsonl"),
            model_output_dir: PathBuf::from(DEFAULT_MODEL_OUTPUT_DIR),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunConfig {
    pub population_size: usize,
    pub time_steps: usize,
    pub calibration_mode: CalibrationMode,
    pub log_individuals: bool,
    pub log_infection_journeys: bool,
    pub rng_seed: Option<u64>,
    pub config_validation_mode: Option<ConfigValidationMode>,
    pub outputs: RunOutputPaths,
    sources: RunConfigSources,
    source_path: Option<PathBuf>,
    source_sha256: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfigValueSource {
    ModelDefault,
    RuntimeConfig,
}

impl ConfigValueSource {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ModelDefault => "model_default",
            Self::RuntimeConfig => "runtime_config",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunConfigSources {
    pub population_size: ConfigValueSource,
    pub time_steps: ConfigValueSource,
    pub calibration_mode: ConfigValueSource,
    pub log_individuals: ConfigValueSource,
    pub log_infection_journeys: ConfigValueSource,
    pub rng_seed: ConfigValueSource,
    pub config_validation_mode: ConfigValueSource,
}

impl Default for RunConfigSources {
    fn default() -> Self {
        Self {
            population_size: ConfigValueSource::ModelDefault,
            time_steps: ConfigValueSource::ModelDefault,
            calibration_mode: ConfigValueSource::ModelDefault,
            log_individuals: ConfigValueSource::ModelDefault,
            log_infection_journeys: ConfigValueSource::ModelDefault,
            rng_seed: ConfigValueSource::ModelDefault,
            config_validation_mode: ConfigValueSource::ModelDefault,
        }
    }
}

impl Default for RunConfig {
    fn default() -> Self {
        let calibration_mode = CalibrationMode::Partial;
        Self {
            population_size: DEFAULT_POPULATION_SIZE,
            time_steps: default_time_steps(calibration_mode),
            calibration_mode,
            log_individuals: false,
            log_infection_journeys: false,
            rng_seed: None,
            config_validation_mode: None,
            outputs: RunOutputPaths::default(),
            sources: RunConfigSources::default(),
            source_path: None,
            source_sha256: None,
        }
    }
}

impl RunConfig {
    pub fn load_from_env() -> Result<Self, RunConfigError> {
        let Some(value) = std::env::var_os(AMR_RUN_CONFIG_ENV) else {
            return Ok(Self::default());
        };
        if value.is_empty() {
            return Err(RunConfigError::Invalid(format!(
                "{AMR_RUN_CONFIG_ENV} is set but empty"
            )));
        }

        let path = PathBuf::from(value);
        let bytes = fs::read(&path).map_err(|source| RunConfigError::Read {
            path: path.clone(),
            source,
        })?;
        Self::from_slice(&bytes, Some(path))
    }

    pub fn source_path(&self) -> Option<&Path> {
        self.source_path.as_deref()
    }

    pub fn source_sha256(&self) -> Option<&str> {
        self.source_sha256.as_deref()
    }

    pub const fn sources(&self) -> RunConfigSources {
        self.sources
    }

    fn from_slice(bytes: &[u8], source_path: Option<PathBuf>) -> Result<Self, RunConfigError> {
        let source_sha256 = source_path
            .as_ref()
            .map(|_| format!("{:x}", Sha256::digest(bytes)));
        let document: RunConfigDocument =
            serde_json::from_slice(bytes).map_err(|source| RunConfigError::Parse {
                path: source_path.clone(),
                source,
            })?;

        require_contract_value("schema_version", &document.schema_version, SCHEMA_VERSION)?;
        require_contract_value("contract", &document.contract, CONTRACT)?;
        require_contract_value(
            "override_semantics",
            &document.override_semantics,
            OVERRIDE_SEMANTICS,
        )?;

        let mut config = Self::default();
        if let Some(mode) = document.overrides.calibration_mode {
            config.calibration_mode = parse_calibration_mode(&mode)?;
            config.time_steps = default_time_steps(config.calibration_mode);
            config.sources.calibration_mode = ConfigValueSource::RuntimeConfig;
        }
        if let Some(value) = document.overrides.population_size {
            config.population_size = require_positive("overrides.population_size", value)?;
            config.sources.population_size = ConfigValueSource::RuntimeConfig;
        }
        if let Some(value) = document.overrides.time_steps {
            config.time_steps = require_positive("overrides.time_steps", value)?;
            config.sources.time_steps = ConfigValueSource::RuntimeConfig;
        }
        if let Some(value) = document.overrides.log_individuals {
            config.log_individuals = value;
            config.sources.log_individuals = ConfigValueSource::RuntimeConfig;
        }
        if let Some(value) = document.overrides.log_infection_journeys {
            config.log_infection_journeys = value;
            config.sources.log_infection_journeys = ConfigValueSource::RuntimeConfig;
        }
        config.rng_seed = document.overrides.rng_seed;
        if config.rng_seed.is_some() {
            config.sources.rng_seed = ConfigValueSource::RuntimeConfig;
        }
        config.config_validation_mode = document
            .overrides
            .config_validation_mode
            .map(|value| {
                if value == "strict" {
                    Ok(ConfigValidationMode::Strict)
                } else {
                    Err(RunConfigError::Invalid(format!(
                        "overrides.config_validation_mode must be 'strict', got '{value}'"
                    )))
                }
            })
            .transpose()?;
        if config.config_validation_mode.is_some() {
            config.sources.config_validation_mode = ConfigValueSource::RuntimeConfig;
        }
        config.outputs = document.outputs.try_into()?;
        config.source_path = source_path;
        config.source_sha256 = source_sha256;
        Ok(config)
    }
}

#[derive(Debug)]
pub enum RunConfigError {
    Read {
        path: PathBuf,
        source: std::io::Error,
    },
    Parse {
        path: Option<PathBuf>,
        source: serde_json::Error,
    },
    Invalid(String),
}

impl fmt::Display for RunConfigError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Read { path, source } => {
                write!(formatter, "unable to read {}: {source}", path.display())
            }
            Self::Parse { path, source } => match path {
                Some(path) => write!(
                    formatter,
                    "invalid {SCHEMA_VERSION} document {}: {source}",
                    path.display()
                ),
                None => write!(formatter, "invalid {SCHEMA_VERSION} document: {source}"),
            },
            Self::Invalid(message) => formatter.write_str(message),
        }
    }
}

impl std::error::Error for RunConfigError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Read { source, .. } => Some(source),
            Self::Parse { source, .. } => Some(source),
            Self::Invalid(_) => None,
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RunConfigDocument {
    schema_version: String,
    contract: String,
    override_semantics: String,
    overrides: RunOverrides,
    outputs: RunOutputDocument,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct RunOverrides {
    population_size: Option<usize>,
    time_steps: Option<usize>,
    rng_seed: Option<u64>,
    config_validation_mode: Option<String>,
    calibration_mode: Option<String>,
    log_individuals: Option<bool>,
    log_infection_journeys: Option<bool>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RunOutputDocument {
    run_root: String,
    report_json: String,
    progress_json: String,
    progress_jsonl: String,
    model_output_dir: String,
}

impl TryFrom<RunOutputDocument> for RunOutputPaths {
    type Error = RunConfigError;

    fn try_from(value: RunOutputDocument) -> Result<Self, Self::Error> {
        let paths = Self {
            run_root: required_path("outputs.run_root", value.run_root)?,
            report_json: required_path("outputs.report_json", value.report_json)?,
            progress_json: required_path("outputs.progress_json", value.progress_json)?,
            progress_jsonl: required_path("outputs.progress_jsonl", value.progress_jsonl)?,
            model_output_dir: required_path("outputs.model_output_dir", value.model_output_dir)?,
        };
        paths.validate()?;
        Ok(paths)
    }
}

impl RunOutputPaths {
    fn validate(&self) -> Result<(), RunConfigError> {
        let files = [
            ("outputs.report_json", &self.report_json),
            ("outputs.progress_json", &self.progress_json),
            ("outputs.progress_jsonl", &self.progress_jsonl),
        ];
        for (index, (left_name, left_path)) in files.iter().enumerate() {
            for (right_name, right_path) in files.iter().skip(index + 1) {
                if left_path == right_path {
                    return Err(RunConfigError::Invalid(format!(
                        "{left_name} and {right_name} must use different paths"
                    )));
                }
            }
        }
        Ok(())
    }
}

fn default_time_steps(mode: CalibrationMode) -> usize {
    match mode {
        CalibrationMode::None => 38_325,
        CalibrationMode::Partial | CalibrationMode::FullMinimal | CalibrationMode::Full => 35_040,
    }
}

fn parse_calibration_mode(value: &str) -> Result<CalibrationMode, RunConfigError> {
    match value {
        "none" => Ok(CalibrationMode::None),
        "partial" => Ok(CalibrationMode::Partial),
        "full_minimal" => Ok(CalibrationMode::FullMinimal),
        "full" => Ok(CalibrationMode::Full),
        _ => Err(RunConfigError::Invalid(format!(
            "overrides.calibration_mode must be one of none, partial, full_minimal, or full; got '{value}'"
        ))),
    }
}

fn require_contract_value(field: &str, actual: &str, expected: &str) -> Result<(), RunConfigError> {
    if actual == expected {
        Ok(())
    } else {
        Err(RunConfigError::Invalid(format!(
            "{field} must be '{expected}', got '{actual}'"
        )))
    }
}

fn require_positive(field: &str, value: usize) -> Result<usize, RunConfigError> {
    if value > 0 {
        Ok(value)
    } else {
        Err(RunConfigError::Invalid(format!(
            "{field} must be greater than zero"
        )))
    }
}

fn required_path(field: &str, value: String) -> Result<PathBuf, RunConfigError> {
    if value.trim().is_empty() {
        Err(RunConfigError::Invalid(format!(
            "{field} must not be empty"
        )))
    } else {
        Ok(PathBuf::from(value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{json, Value};

    fn document(overrides: Value) -> Vec<u8> {
        serde_json::to_vec(&json!({
            "schema_version": SCHEMA_VERSION,
            "contract": CONTRACT,
            "override_semantics": OVERRIDE_SEMANTICS,
            "overrides": overrides,
            "outputs": {
                "run_root": "/run",
                "report_json": "/run/amr_report.json",
                "progress_json": "/run/amr_progress.json",
                "progress_jsonl": "/run/amr_progress.jsonl",
                "model_output_dir": "/run/model-output"
            }
        }))
        .expect("test runtime config must serialize")
    }

    #[test]
    fn omitted_overrides_preserve_model_defaults() {
        let config = RunConfig::from_slice(&document(json!({})), None).unwrap();

        assert_eq!(
            config,
            RunConfig {
                outputs: RunOutputPaths {
                    run_root: PathBuf::from("/run"),
                    report_json: PathBuf::from("/run/amr_report.json"),
                    progress_json: PathBuf::from("/run/amr_progress.json"),
                    progress_jsonl: PathBuf::from("/run/amr_progress.jsonl"),
                    model_output_dir: PathBuf::from("/run/model-output"),
                },
                ..RunConfig::default()
            }
        );
    }

    #[test]
    fn explicit_overrides_replace_only_named_model_values() {
        let config = RunConfig::from_slice(
            &document(json!({
                "population_size": 64,
                "time_steps": 4,
                "calibration_mode": "full",
                "rng_seed": 1729,
                "config_validation_mode": "strict",
                "log_individuals": true,
                "log_infection_journeys": true
            })),
            None,
        )
        .unwrap();

        assert_eq!(config.population_size, 64);
        assert_eq!(config.time_steps, 4);
        assert_eq!(config.calibration_mode, CalibrationMode::Full);
        assert_eq!(config.rng_seed, Some(1729));
        assert_eq!(
            config.config_validation_mode,
            Some(ConfigValidationMode::Strict)
        );
        assert!(config.log_individuals);
        assert!(config.log_infection_journeys);
        assert_eq!(
            config.sources(),
            RunConfigSources {
                population_size: ConfigValueSource::RuntimeConfig,
                time_steps: ConfigValueSource::RuntimeConfig,
                calibration_mode: ConfigValueSource::RuntimeConfig,
                log_individuals: ConfigValueSource::RuntimeConfig,
                log_infection_journeys: ConfigValueSource::RuntimeConfig,
                rng_seed: ConfigValueSource::RuntimeConfig,
                config_validation_mode: ConfigValueSource::RuntimeConfig,
            }
        );
    }

    #[test]
    fn mode_override_uses_the_modes_default_horizon_when_time_steps_are_omitted() {
        let config =
            RunConfig::from_slice(&document(json!({"calibration_mode": "none"})), None).unwrap();

        assert_eq!(config.time_steps, 38_325);
    }

    #[test]
    fn malformed_or_extended_contracts_fail_closed() {
        for bytes in [
            document(json!({"population_size": 0})),
            document(json!({"time_steps": 0})),
            document(json!({"calibration_mode": "maximum"})),
            document(json!({"platform_resource_policy": {"cpu": 72}})),
        ] {
            assert!(RunConfig::from_slice(&bytes, None).is_err());
        }

        let mut payload: Value = serde_json::from_slice(&document(json!({}))).unwrap();
        payload["copy_identity"] = json!({"copy_index": 1});
        let bytes = serde_json::to_vec(&payload).unwrap();
        assert!(RunConfig::from_slice(&bytes, None).is_err());

        let mut payload: Value = serde_json::from_slice(&document(json!({}))).unwrap();
        payload["outputs"]["progress_jsonl"] = payload["outputs"]["progress_json"].clone();
        let bytes = serde_json::to_vec(&payload).unwrap();
        assert!(RunConfig::from_slice(&bytes, None).is_err());
    }
}
