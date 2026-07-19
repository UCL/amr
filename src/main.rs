// src/main.rs
// Simulation entry point.
//
// This file is intentionally small: it chooses one run configuration, builds a `Simulation`,
// executes it, and writes the main CSV output. Most model behaviour lives in `rules/`,
// `simulation/`, and `config.rs`; edit this file when you want to change how a run is launched.
//
// Useful local workflow note:
//   PowerShell: `$env:RAYON_NUM_THREADS = "4"` before running to cap parallelism.

//
//
//

use amr_project::config::{get_global_param, PARAMETERS};
use amr_project::config_validation::{validate_parameter_map, ConfigValidationMode};
use amr_project::observability;
use amr_project::run_config::{ConfigValueSource, RunConfig};
use amr_project::simulation::population::BACTERIA_LIST;
use amr_project::simulation::simulation::CalibrationMode;
use amr_project::simulation::simulation::Simulation;
use chrono::Utc;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::backtrace::Backtrace;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::path::PathBuf;

const RAYON_WORKER_STACK_BYTES: usize = 4 * 1024 * 1024;

#[derive(Clone, Copy, Debug)]
struct ResolvedRunSeed {
    value: u64,
    source: &'static str,
}

#[derive(Serialize)]
struct EffectiveValue<T> {
    value: T,
    source: &'static str,
}

#[derive(Serialize)]
struct EffectiveRunConfig {
    population_size: EffectiveValue<usize>,
    time_steps: EffectiveValue<usize>,
    calibration_mode: EffectiveValue<String>,
    log_individuals: EffectiveValue<bool>,
    log_infection_journeys: EffectiveValue<bool>,
    rng_seed: EffectiveValue<u64>,
    config_validation_mode: EffectiveValue<String>,
}

#[derive(Serialize)]
struct AmrFinalReport {
    schema_version: &'static str,
    status: &'static str,
    exit_code: i32,
    finished_at_utc: String,
    source_hash: String,
    runtime_config_consumed: bool,
    runtime_config_sha256: Option<String>,
    run_id: Option<u32>,
    last_timestep: Option<usize>,
    duration_seconds: Option<f64>,
    summary_csv: Option<String>,
    summary_sha256: Option<String>,
    failure_class: String,
    effective_config: EffectiveRunConfig,
    source: &'static str,
}

fn config_source(source: ConfigValueSource) -> &'static str {
    source.as_str()
}

fn write_final_report(path: &std::path::Path, report: &AmrFinalReport) -> std::io::Result<()> {
    observability::write_json_atomically(path, report)
}

fn effective_run_config(
    run_config: &RunConfig,
    seed: ResolvedRunSeed,
    config_validation_mode: ConfigValidationMode,
    config_validation_source: &'static str,
) -> EffectiveRunConfig {
    let sources = run_config.sources();
    EffectiveRunConfig {
        population_size: EffectiveValue {
            value: run_config.population_size,
            source: config_source(sources.population_size),
        },
        time_steps: EffectiveValue {
            value: run_config.time_steps,
            source: config_source(sources.time_steps),
        },
        calibration_mode: EffectiveValue {
            value: run_config.calibration_mode.to_string(),
            source: config_source(sources.calibration_mode),
        },
        log_individuals: EffectiveValue {
            value: run_config.log_individuals,
            source: config_source(sources.log_individuals),
        },
        log_infection_journeys: EffectiveValue {
            value: run_config.log_infection_journeys,
            source: config_source(sources.log_infection_journeys),
        },
        rng_seed: EffectiveValue {
            value: seed.value,
            source: seed.source,
        },
        config_validation_mode: EffectiveValue {
            value: config_validation_mode.to_string(),
            source: config_validation_source,
        },
    }
}

fn configure_rayon_worker_stack() {
    rayon::ThreadPoolBuilder::new()
        .stack_size(RAYON_WORKER_STACK_BYTES)
        .thread_name(|idx| format!("amr-rayon-{}", idx))
        .build_global()
        .expect("failed to configure global Rayon thread pool before first use");

    eprintln!(
        "[startup] configured global Rayon worker stack: {} bytes",
        RAYON_WORKER_STACK_BYTES
    );
}

fn resolve_run_seed(
    runtime_config_seed: Option<u64>,
    environment_seed: Option<&str>,
    use_fixed_seed: bool,
    fixed_seed_value: u64,
) -> Result<ResolvedRunSeed, String> {
    let environment_seed = environment_seed
        .map(|value| {
            value
                .parse::<u64>()
                .map_err(|_| format!("AMR_RNG_SEED must parse as a u64, got '{value}'"))
        })
        .transpose()?;

    if let Some(value) = runtime_config_seed {
        if let Some(environment) = environment_seed {
            if environment != value {
                return Err(format!(
                    "runtime config rng_seed {value} conflicts with AMR_RNG_SEED {environment}"
                ));
            }
        }
        eprintln!("[startup] AMR_RNG_SEED={} source=runtime_config", value);
        return Ok(ResolvedRunSeed {
            value,
            source: "runtime_config",
        });
    }

    if let Some(parsed) = environment_seed {
        eprintln!("[startup] AMR_RNG_SEED={} source=env", parsed);
        return Ok(ResolvedRunSeed {
            value: parsed,
            source: "env",
        });
    }

    if use_fixed_seed {
        eprintln!(
            "[startup] AMR_RNG_SEED={} source=fixed_seed_value",
            fixed_seed_value
        );
        return Ok(ResolvedRunSeed {
            value: fixed_seed_value,
            source: "fixed_seed_value",
        });
    }

    let generated = rand::random::<u64>();
    eprintln!("[startup] AMR_RNG_SEED={} source=generated", generated);
    Ok(ResolvedRunSeed {
        value: generated,
        source: "generated",
    })
}

#[cfg(test)]
mod seed_tests {
    use super::*;

    #[test]
    fn runtime_seed_is_explicit_and_must_match_the_environment() {
        let seed = resolve_run_seed(Some(1729), Some("1729"), false, 99).unwrap();
        assert_eq!(seed.value, 1729);
        assert_eq!(seed.source, "runtime_config");

        let error = resolve_run_seed(Some(1729), Some("1730"), false, 99).unwrap_err();
        assert!(error.contains("conflicts"));
    }

    #[test]
    fn model_generates_a_seed_when_no_override_exists() {
        let seed = resolve_run_seed(None, None, false, 99).unwrap();
        assert_eq!(seed.source, "generated");
    }

    #[test]
    fn legacy_environment_and_fixed_seed_fallbacks_are_preserved() {
        let environment = resolve_run_seed(None, Some("42"), false, 99).unwrap();
        assert_eq!((environment.value, environment.source), (42, "env"));

        let fixed = resolve_run_seed(None, None, true, 99).unwrap();
        assert_eq!((fixed.value, fixed.source), (99, "fixed_seed_value"));
    }
}

fn resolve_source_hash() -> String {
    observability::resolve_source_hash()
}

fn hash_file_sha256(path: &std::path::Path) -> std::io::Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 8192];

    loop {
        let bytes_read = file.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        hasher.update(&buffer[..bytes_read]);
    }

    Ok(format!("{:x}", hasher.finalize()))
}

fn classify_panic(payload: &str, location: &str) -> &'static str {
    let payload_lower = payload.to_ascii_lowercase();
    let location_lower = location.to_ascii_lowercase();

    if payload_lower.contains("stack overflow") {
        "stack_overflow"
    } else if payload_lower.contains("index out of bounds") {
        "index_out_of_bounds"
    } else if payload_lower.contains("attempt to divide by zero") {
        "divide_by_zero"
    } else if payload_lower.contains("assertion failed")
        || location_lower.contains("panic_bounds_check")
    {
        "assertion_or_bounds_check"
    } else {
        "panic"
    }
}

fn write_run_metadata(
    path: &std::path::Path,
    status: &str,
    source_hash: &str,
    seed: ResolvedRunSeed,
    population_size: usize,
    time_steps: usize,
    calibration_mode: CalibrationMode,
    active_policies: &[u8],
    run_id: Option<u32>,
    csv_path: Option<&std::path::Path>,
    duration_secs: Option<f64>,
    summary_hash: Option<&str>,
    failure_class: &str,
    last_timestep: Option<usize>,
    config_validation_mode: &str,
    config_validation_status: &str,
    config_validation_errors: usize,
    config_validation_warnings: usize,
    config_validation_report_path: Option<&std::path::Path>,
) -> std::io::Result<()> {
    let mut file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(path)?;

    writeln!(file, "status={}", status)?;
    writeln!(
        file,
        "updated_utc={}",
        Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
    )?;
    writeln!(file, "source_hash={}", source_hash)?;
    writeln!(file, "rng_seed={}", seed.value)?;
    writeln!(file, "rng_seed_source={}", seed.source)?;
    writeln!(
        file,
        "run_id={}",
        run_id.map_or_else(|| "pending".to_string(), |id| id.to_string())
    )?;
    writeln!(file, "population_size={}", population_size)?;
    writeln!(file, "time_steps={}", time_steps)?;
    writeln!(file, "calibration_mode={}", calibration_mode)?;
    writeln!(file, "active_policies={:?}", active_policies)?;
    writeln!(
        file,
        "last_timestep={}",
        last_timestep.map_or_else(|| "pending".to_string(), |step| step.to_string())
    )?;
    let rayon_threads = rayon::current_num_threads();
    writeln!(file, "rayon_threads={}", rayon_threads)?;
    writeln!(
        file,
        "rayon_worker_stack_bytes={}",
        RAYON_WORKER_STACK_BYTES
    )?;
    writeln!(file, "config_validation_mode={}", config_validation_mode)?;
    writeln!(
        file,
        "config_validation_status={}",
        config_validation_status
    )?;
    writeln!(
        file,
        "config_validation_errors={}",
        config_validation_errors
    )?;
    writeln!(
        file,
        "config_validation_warnings={}",
        config_validation_warnings
    )?;
    writeln!(
        file,
        "config_validation_report={}",
        config_validation_report_path
            .map(|path| path.display().to_string())
            .unwrap_or_else(|| "pending".to_string())
    )?;
    writeln!(
        file,
        "duration_seconds={}",
        duration_secs.map_or_else(|| "pending".to_string(), |secs| format!("{:.3}", secs))
    )?;
    writeln!(
        file,
        "summary_csv={}",
        csv_path
            .map(|path| path.display().to_string())
            .unwrap_or_else(|| "pending".to_string())
    )?;
    writeln!(file, "summary_hash={}", summary_hash.unwrap_or("pending"))?;
    writeln!(file, "failure_class={}", failure_class)?;
    writeln!(file, "replay_env=AMR_RNG_SEED={}", seed.value)?;

    Ok(())
}

fn main() {
    install_panic_log_hook();
    configure_rayon_worker_stack();
    let _ = env_logger::builder().is_test(false).try_init();

    let run_config = match RunConfig::load_from_env() {
        Ok(config) => config,
        Err(err) => {
            eprintln!("[startup] invalid AMR runtime configuration: {err}");
            std::process::exit(2);
        }
    };
    match run_config.source_path() {
        Some(path) => eprintln!(
            "[startup] AMR_RUN_CONFIG={} source=runtime_config",
            path.display()
        ),
        None => eprintln!("[startup] AMR_RUN_CONFIG=unset source=model_defaults"),
    }
    let runtime_config_consumed = run_config.source_path().is_some();
    let runtime_config_sha256 = run_config.source_sha256().map(str::to_owned);

    let population_size = run_config.population_size;
    let calibration_mode = run_config.calibration_mode;
    let time_steps = run_config.time_steps;
    let log_individuals = run_config.log_individuals;
    let log_infection_journeys = run_config.log_infection_journeys;
    let run_config_sources = run_config.sources();
    let infection_journey_sample_rate = 1.00; // Fraction of eligible infections to log when journey logging is enabled.
    let use_fixed_seed = false; // Toggle to enable deterministic RNG seeding
    let fixed_seed_value: u64 = 1_234_567_890; // Seed used when use_fixed_seed is true
    let infection_journey_bacteria_filter: Option<&str> = None; // Set to Some("escherichia_coli") to log only specific bacteria

    // Examples of bacteria filter values (use lowercase with underscores):
    // Some("escherichia_coli")
    // Some("staphylococcus_aureus")
    // Some("pseudomonas_aeruginosa")
    // Some("acinetobacter_baumannii")
    // Some("enterococcus_faecium")
    // None - logs all bacteria types

    let environment_seed = std::env::var("AMR_RNG_SEED").ok();
    let resolved_run_seed = match resolve_run_seed(
        run_config.rng_seed,
        environment_seed.as_deref(),
        use_fixed_seed,
        fixed_seed_value,
    ) {
        Ok(seed) => seed,
        Err(err) => {
            eprintln!("[startup] invalid RNG seed configuration: {err}");
            std::process::exit(2);
        }
    };
    std::env::set_var("AMR_RNG_SEED", resolved_run_seed.value.to_string());
    let seed_override = Some(resolved_run_seed.value);
    let source_hash = resolve_source_hash();
    eprintln!("[startup] source_hash={}", source_hash);

    // ── Policy branch selection ────────────────────────────────────────────────
    // Only active when calibration_mode == CalibrationMode::None (full run).
    // List the policy IDs you want to run; comment out any you don't need.
    //
    //   0 = Baseline continuation      (status quo carried forward to 2035)
    //   1 = Antimicrobial Stewardship  (reduced prescribing, better drug selection)
    //   2 = AMR Counterfactual         (resistance cleared; models a world without AMR)
    //   3 = Perfect Diagnostics        (implausibly complete & immediate testing)
    //   4 = Equal Global Access        (all regions get North America–level access)
    //
    // Policies are independent branches starting from POLICY_BRANCH_YEAR (2027);
    // they do not interact with each other.
    let active_policies: &[u8] = &[
        0, // Baseline continuation
        1, // Stewardship
        2, // AMR counterfactual
        3, // Perfect diagnostics
        4, // Equal global access
    ];

    let output_dir = run_config.outputs.model_output_dir.as_path();
    if let Err(err) = std::fs::create_dir_all(output_dir) {
        eprintln!(
            "[startup] unable to create output directory {}: {err}",
            output_dir.display()
        );
        std::process::exit(2);
    }
    std::env::set_var("AMR_REPORT_JSON", &run_config.outputs.report_json);
    if let Err(err) = observability::configure_progress(observability::ProgressConfig {
        snapshot_path: run_config.outputs.progress_json.clone(),
        events_path: run_config.outputs.progress_jsonl.clone(),
        population_size,
        total_steps: time_steps,
        calibration_mode: calibration_mode.to_string(),
    }) {
        eprintln!("[startup] unable to initialize AMR progress evidence: {err}");
        std::process::exit(2);
    }

    let metadata_stamp = Utc::now().format("%Y%m%dT%H%M%SZ");
    let metadata_path = output_dir.join(format!(
        "run_metadata_{}_seed_{}.txt",
        metadata_stamp, resolved_run_seed.value
    ));

    let config_validation_source = if run_config.config_validation_mode.is_some() {
        config_source(run_config_sources.config_validation_mode)
    } else if std::env::var_os("AMR_CONFIG_VALIDATION").is_some() {
        "environment"
    } else {
        "model_default"
    };
    let config_validation_mode = match run_config.config_validation_mode {
        Some(mode) => mode,
        None => match ConfigValidationMode::from_env() {
            Ok(mode) => mode,
            Err(err) => {
                eprintln!("[config-validation] FAILED: {}", err);
                std::process::exit(2);
            }
        },
    };
    let config_validation_report = validate_parameter_map(&PARAMETERS);
    let rendered_config_validation = config_validation_report.render(config_validation_mode);
    eprint!("{}", rendered_config_validation);
    let config_validation_report_path =
        output_dir.join(format!("config_validation_{}.txt", metadata_stamp));
    let config_validation_report_for_metadata = match File::create(&config_validation_report_path)
        .and_then(|mut file| file.write_all(rendered_config_validation.as_bytes()))
    {
        Ok(()) => Some(config_validation_report_path.as_path()),
        Err(err) => {
            eprintln!(
                "Warning: unable to write config validation report {}: {}",
                config_validation_report_path.display(),
                err
            );
            None
        }
    };

    if config_validation_report.has_errors() && config_validation_mode.blocks_on_errors() {
        if let Err(err) = write_run_metadata(
            &metadata_path,
            "config_validation_failed",
            &source_hash,
            resolved_run_seed,
            population_size,
            time_steps,
            calibration_mode,
            active_policies,
            None,
            None,
            None,
            None,
            "config_validation_failed",
            None,
            &config_validation_mode.to_string(),
            config_validation_report.status(),
            config_validation_report.error_count(),
            config_validation_report.warning_count(),
            config_validation_report_for_metadata,
        ) {
            eprintln!(
                "Warning: unable to write failed run metadata {}: {}",
                metadata_path.display(),
                err
            );
        }

        let report = AmrFinalReport {
            schema_version: "amr-report/v1",
            status: "failed",
            exit_code: 2,
            finished_at_utc: Utc::now().to_rfc3339(),
            source_hash: source_hash.clone(),
            runtime_config_consumed,
            runtime_config_sha256: runtime_config_sha256.clone(),
            run_id: None,
            last_timestep: None,
            duration_seconds: None,
            summary_csv: None,
            summary_sha256: None,
            failure_class: "config_validation_failed".to_string(),
            effective_config: effective_run_config(
                &run_config,
                resolved_run_seed,
                config_validation_mode,
                config_validation_source,
            ),
            source: "amr_model",
        };
        if let Err(err) = write_final_report(&run_config.outputs.report_json, &report) {
            eprintln!("[report] unable to write AMR final report: {err}");
        }
        if let Err(err) =
            observability::publish_progress(observability::ProgressStatus::Failed, None)
        {
            eprintln!("[progress] unable to write terminal AMR progress: {err}");
        }

        println!(
            "[report] status=config_validation_failed run_id=pending source_hash={} rng_seed={} last_timestep=pending summary_hash=pending failure_class=config_validation_failed config_validation_status={} config_validation_mode={} summary_csv=pending",
            source_hash,
            resolved_run_seed.value,
            config_validation_report.status(),
            config_validation_mode
        );
        std::process::exit(2);
    }

    // Fail fast on obviously inconsistent bacteria setups before paying the cost of a full run.
    validate_bacteria_configuration();

    if let Err(err) = write_run_metadata(
        &metadata_path,
        "started",
        &source_hash,
        resolved_run_seed,
        population_size,
        time_steps,
        calibration_mode,
        active_policies,
        None,
        None,
        None,
        None,
        "running",
        None,
        &config_validation_mode.to_string(),
        config_validation_report.status(),
        config_validation_report.error_count(),
        config_validation_report.warning_count(),
        config_validation_report_for_metadata,
    ) {
        eprintln!(
            "Warning: unable to write run metadata {}: {}",
            metadata_path.display(),
            err
        );
    } else {
        eprintln!("[startup] run metadata: {}", metadata_path.display());
    }

    let mut simulation = Simulation::new(
        population_size,
        time_steps,
        log_individuals,
        seed_override,
        calibration_mode,
    );
    let use_disk_branch_checkpointing = calibration_mode == CalibrationMode::None;
    let disk_checkpoint_directory: Option<PathBuf> = None; // Override with Some(path) to specify a custom folder

    if use_disk_branch_checkpointing {
        simulation.enable_disk_branch_checkpointing(disk_checkpoint_directory);
    } else {
        simulation.disable_disk_branch_checkpointing();
    }

    // Journey logging is optional because it writes a much richer, more expensive trace than the summary CSV.
    if log_infection_journeys {
        match infection_journey_bacteria_filter {
            Some(filter) => simulation.enable_infection_journey_logging_with_filter(
                infection_journey_sample_rate,
                Some(filter.to_string()),
            ),
            None => simulation.enable_infection_journey_logging(infection_journey_sample_rate),
        }
    }

    simulation.set_active_policy_branches(active_policies);
    // ──────────────────────────────────────────────────────────────────────────

    use std::time::Instant;
    let start = Instant::now();

    simulation.run();

    let duration = start.elapsed();

    // Print summary statistics from logged data
    simulation.print_summary_statistics();

    // Use the simulation's run_id in the filename so Python post-processing can join one
    // summary CSV to one set of sampled parameters and one run-log entry without collisions.
    let run_id = simulation.run_id;
    let csv_basename = format!("simulation_summary_{:06}.csv", run_id);
    let csv_path = output_dir.join(&csv_basename);

    // The summary CSV is the primary handoff to the Python analysis scripts.
    let (summary_hash, mut failure_class, mut exit_code) =
        match simulation.export_summary_to_csv(&csv_path) {
            Ok(()) => {
                println!("Summary data exported to {}", csv_path.display());
                match hash_file_sha256(&csv_path) {
                    Ok(hash) => (Some(hash), "completed".to_string(), 0),
                    Err(err) => {
                        eprintln!(
                            "Warning: unable to hash summary CSV {}: {}",
                            csv_path.display(),
                            err
                        );
                        (None, "summary_hash_failed".to_string(), 1)
                    }
                }
            }
            Err(err) => {
                println!("Error exporting CSV: {}", err);
                (None, "csv_export_failed".to_string(), 1)
            }
        };
    let last_timestep = observability::current_timestep().or_else(|| time_steps.checked_sub(1));
    let summary_csv_path = csv_path.is_file().then_some(csv_path.as_path());

    let report = AmrFinalReport {
        schema_version: "amr-report/v1",
        status: if exit_code == 0 {
            "completed"
        } else {
            "failed"
        },
        exit_code,
        finished_at_utc: Utc::now().to_rfc3339(),
        source_hash: source_hash.clone(),
        runtime_config_consumed,
        runtime_config_sha256,
        run_id: Some(run_id),
        last_timestep,
        duration_seconds: Some(duration.as_secs_f64()),
        summary_csv: summary_csv_path.map(|path| path.display().to_string()),
        summary_sha256: summary_hash.clone(),
        failure_class: failure_class.clone(),
        effective_config: effective_run_config(
            &run_config,
            resolved_run_seed,
            config_validation_mode,
            config_validation_source,
        ),
        source: "amr_model",
    };
    if let Err(err) = write_final_report(&run_config.outputs.report_json, &report) {
        eprintln!("[report] unable to write AMR final report: {err}");
        failure_class = "report_write_failed".to_string();
        exit_code = 3;
    }
    let final_status = if exit_code == 0 {
        "completed"
    } else {
        "failed"
    };

    if let Err(err) = write_run_metadata(
        &metadata_path,
        final_status,
        &source_hash,
        resolved_run_seed,
        population_size,
        time_steps,
        calibration_mode,
        active_policies,
        Some(run_id),
        summary_csv_path,
        Some(duration.as_secs_f64()),
        summary_hash.as_deref(),
        &failure_class,
        last_timestep,
        &config_validation_mode.to_string(),
        config_validation_report.status(),
        config_validation_report.error_count(),
        config_validation_report.warning_count(),
        config_validation_report_for_metadata,
    ) {
        eprintln!(
            "Warning: unable to update run metadata {}: {}",
            metadata_path.display(),
            err
        );
    }

    let progress_step = if exit_code == 0 {
        Some(time_steps)
    } else {
        last_timestep
    };
    let progress_status = if exit_code == 0 {
        observability::ProgressStatus::Completed
    } else {
        observability::ProgressStatus::Failed
    };
    if let Err(err) = observability::publish_progress(progress_status, progress_step) {
        eprintln!("[progress] unable to write terminal AMR progress: {err}");
    }

    println!(
        "[report] status={} run_id={} source_hash={} rng_seed={} last_timestep={} summary_hash={} failure_class={} config_validation_status={} config_validation_mode={} summary_csv={}",
        final_status,
        run_id,
        source_hash,
        resolved_run_seed.value,
        last_timestep
            .map(|step| step.to_string())
            .unwrap_or_else(|| "pending".to_string()),
        summary_hash.as_deref().unwrap_or("pending"),
        failure_class,
        config_validation_report.status(),
        config_validation_mode,
        csv_path.display()
    );

    // Log the simulation run details
    if let Err(e) = log_simulation_run(population_size, time_steps, duration.as_secs_f64()) {
        eprintln!("Error logging simulation run: {}", e);
    }

    println!("\n--- simulation ended ---");
    println!(
        "--- total simulation time: {:.3} seconds",
        duration.as_secs_f64()
    );
    println!("                          ");
    if exit_code != 0 {
        std::process::exit(exit_code);
    }
}

fn install_panic_log_hook() {
    std::panic::set_hook(Box::new(|panic_info| {
        let timestamp = Utc::now().format("%Y-%m-%d %H:%M:%S UTC");
        let rng_seed = std::env::var("AMR_RNG_SEED").unwrap_or_else(|_| "unset".to_string());
        let source_hash = resolve_source_hash();
        let run_id = observability::current_run_id()
            .map(|id| id.to_string())
            .unwrap_or_else(|| "unset".to_string());
        let current_timestep = observability::current_timestep();
        let last_timestep = current_timestep
            .map(|step| step.to_string())
            .unwrap_or_else(|| "unset".to_string());
        let location = panic_info
            .location()
            .map(|loc| format!("{}:{}:{}", loc.file(), loc.line(), loc.column()))
            .unwrap_or_else(|| "unknown location".to_string());

        let payload = if let Some(message) = panic_info.payload().downcast_ref::<&str>() {
            *message
        } else if let Some(message) = panic_info.payload().downcast_ref::<String>() {
            message.as_str()
        } else {
            "non-string panic payload"
        };
        let failure_class = classify_panic(payload, &location);

        let report = format!(
            "\n===== panic =====\ntimestamp: {}\nsource_hash: {}\nrng_seed: {}\nrun_id: {}\nlast_timestep: {}\nfailure_class: {}\nlocation: {}\npayload: {}\nbacktrace:\n{}\n",
            timestamp,
            source_hash,
            rng_seed,
            run_id,
            last_timestep,
            failure_class,
            location,
            payload,
            Backtrace::force_capture()
        );

        eprintln!("{}", report);
        if let Err(error) =
            observability::publish_progress(observability::ProgressStatus::Failed, current_timestep)
        {
            eprintln!("[progress] unable to write panic progress evidence: {error}");
        }

        if let Ok(mut file) = OpenOptions::new()
            .create(true)
            .append(true)
            .open("panic_log.txt")
        {
            let _ = file.write_all(report.as_bytes());
        }
    }));
}

/// Append one row to the lightweight run log so wall-clock cost can be compared across runs.
fn log_simulation_run(
    population_size: usize,
    time_steps: usize,
    duration_secs: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    use chrono::Utc;
    use std::fs::OpenOptions;
    use std::io::Write;

    let timestamp = Utc::now();
    let log_entry = format!(
        "{},{},{},{:.3}\n",
        timestamp.format("%Y-%m-%d %H:%M:%S UTC"),
        population_size,
        time_steps,
        duration_secs
    );

    // Create the file lazily so ad hoc runs do not need any manual setup.
    let log_path = "simulation_run_log.csv";
    let file_exists = std::path::Path::new(log_path).exists();

    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(log_path)?;

    // Write header if file is new
    if !file_exists {
        writeln!(
            file,
            "timestamp,population_size,time_steps,duration_seconds"
        )?;
    }

    // Write the log entry
    file.write_all(log_entry.as_bytes())?;

    println!("Simulation run logged to {}", log_path);

    Ok(())
}

/// Validate the bacteria roster before the simulation starts.
///
/// The model can technically run with one or a few organisms, but many ecosystem-level
/// assumptions (HGT, microbiome competition, syndromic empiric therapy) become much less
/// realistic, so this prints warnings rather than silently accepting a misleading setup.
fn validate_bacteria_configuration() {
    let num_bacteria = BACTERIA_LIST.len();

    println!("=== BACTERIA CONFIGURATION VALIDATION ===");
    println!("Number of bacteria in simulation: {}", num_bacteria);

    if num_bacteria == 0 {
        panic!("ERROR: BACTERIA_LIST cannot be empty!");
    }

    if num_bacteria == 1 {
        println!("⚠️  SINGLE-BACTERIA MODE: This limits biological realism but is valid for:");
        println!("   • Pathogen-specific resistance studies");
        println!("   • Drug development against specific organisms");
        println!("   • Educational/training scenarios");
        println!("   • Computational efficiency");
        println!("   Note: HGT, microbiome competition, and syndromic treatment are disabled.");
    } else if num_bacteria < 5 {
        println!("⚠️  LIMITED-BACTERIA MODE: Some ecosystem effects may be reduced");
        println!("   Consider if your research question needs more bacterial diversity");
    } else {
        println!("✓ MULTI-BACTERIA MODE: Full ecosystem modeling enabled");
    }

    // Check for potential HGT configuration issues
    if num_bacteria == 1 {
        for bacteria in BACTERIA_LIST.iter() {
            for other in BACTERIA_LIST.iter() {
                if bacteria != other {
                    let hgt_param_name = format!("hgt_prob_{}_to_{}", bacteria, other);
                    if let Some(hgt_prob) = get_global_param(&hgt_param_name) {
                        if hgt_prob > 0.0 {
                            println!("⚠️  HGT parameters configured but only 1 bacteria present - HGT will not occur");
                            break;
                        }
                    }
                }
            }
        }
    }

    println!("Bacteria included:");
    for (i, bacteria) in BACTERIA_LIST.iter().enumerate() {
        println!("   {}. {}", i + 1, bacteria);
    }
    println!("=====================================\n");
}
