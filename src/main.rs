// src/main.rs
// Simulation entry point.
//
// This launcher chooses the run configuration, validates it, configures run-level observability,
// builds and executes a `Simulation`, and writes its outputs. Model behaviour remains in `rules/`,
// `simulation/`, and `config.rs`; edit this file to change how a run is launched.
//
// Useful local workflow note:
//   PowerShell: `$env:RAYON_NUM_THREADS = "4"` before running to cap parallelism.

use amr_project::config::{get_global_param, PARAMETERS};
use amr_project::config_validation::{validate_parameter_map, ConfigValidationMode};
use amr_project::observability;
use amr_project::output_files::publish_summary_no_clobber;
use amr_project::simulation::population::BACTERIA_LIST;
use amr_project::simulation::simulation::CalibrationMode;
use amr_project::simulation::simulation::Simulation;
use chrono::Utc;
use sha2::{Digest, Sha256};
use std::backtrace::Backtrace;
use std::fs::{File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

const RAYON_WORKER_STACK_BYTES: usize = 4 * 1024 * 1024;

#[derive(Clone, Copy)]
struct ResolvedRunSeed {
    value: u64,
    source: &'static str,
}

struct InvocationPaths {
    metadata_path: PathBuf,
    validation_path: PathBuf,
    validation_file: File,
}

/// Reserve both diagnostic paths before allocating the research population.
/// The clock is only a readable label; exclusive creation supplies uniqueness.
fn reserve_invocation_paths(
    output_dir: &Path,
    seed: u64,
    stamp: &str,
) -> io::Result<InvocationPaths> {
    std::fs::create_dir_all(output_dir)?;
    let mut attempt = 0u64;
    loop {
        let token = format!("{stamp}_seed_{seed}_pid_{}_{attempt}", std::process::id());
        let metadata_path = output_dir.join(format!("run_metadata_{token}.txt"));
        let validation_path = output_dir.join(format!("config_validation_{token}.txt"));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&metadata_path)
        {
            Ok(file) => drop(file),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {
                attempt = attempt
                    .checked_add(1)
                    .ok_or_else(|| io::Error::other("invocation filename space exhausted"))?;
                continue;
            }
            Err(error) => return Err(error),
        }
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&validation_path)
        {
            Ok(validation_file) => {
                return Ok(InvocationPaths {
                    metadata_path,
                    validation_path,
                    validation_file,
                })
            }
            Err(error) => {
                // Only remove the empty metadata reservation created by this attempt.
                std::fs::remove_file(&metadata_path)?;
                if error.kind() != io::ErrorKind::AlreadyExists {
                    return Err(error);
                }
                attempt = attempt
                    .checked_add(1)
                    .ok_or_else(|| io::Error::other("invocation filename space exhausted"))?;
            }
        }
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

fn resolve_run_seed(use_fixed_seed: bool, fixed_seed_value: u64) -> ResolvedRunSeed {
    if let Ok(value) = std::env::var("AMR_RNG_SEED") {
        let parsed = value
            .parse::<u64>()
            .expect("AMR_RNG_SEED must parse as a u64");
        eprintln!("[startup] AMR_RNG_SEED={} source=env", parsed);
        return ResolvedRunSeed {
            value: parsed,
            source: "env",
        };
    }

    if use_fixed_seed {
        eprintln!(
            "[startup] AMR_RNG_SEED={} source=fixed_seed_value",
            fixed_seed_value
        );
        return ResolvedRunSeed {
            value: fixed_seed_value,
            source: "fixed_seed_value",
        };
    }

    let generated = rand::random::<u64>();
    eprintln!("[startup] AMR_RNG_SEED={} source=generated", generated);
    ResolvedRunSeed {
        value: generated,
        source: "generated",
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

#[derive(Debug)]
struct RunOutput {
    status: &'static str,
    csv_path: Option<PathBuf>,
    summary_hash: Option<String>,
    failure_detail: Option<String>,
}

impl RunOutput {
    fn fail(&mut self, status: &'static str, error: io::Error) {
        self.status = status;
        self.failure_detail = Some(error.to_string());
    }

    fn exit_code(&self) -> ExitCode {
        if self.status == "completed" {
            ExitCode::SUCCESS
        } else {
            ExitCode::FAILURE
        }
    }
}

fn reserve_incomplete_file(path: &Path) -> io::Result<(PathBuf, File)> {
    let mut attempt = 0u64;
    loop {
        let mut incomplete_name = path.as_os_str().to_os_string();
        incomplete_name.push(format!(".{}.{attempt}.incomplete", std::process::id()));
        let incomplete_path = PathBuf::from(incomplete_name);
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&incomplete_path)
        {
            Ok(file) => return Ok((incomplete_path, file)),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {
                attempt = attempt
                    .checked_add(1)
                    .ok_or_else(|| io::Error::other("staging filename space exhausted"))?;
            }
            Err(error) => return Err(error),
        }
    }
}

fn write_atomic_metadata(
    path: &Path,
    write: impl FnOnce(&mut File) -> io::Result<()>,
) -> io::Result<()> {
    let (incomplete_path, mut file) = reserve_incomplete_file(path)?;
    write(&mut file)?;
    file.sync_all()?;
    drop(file);
    // A failed write must preserve the previous metadata, usually `status=started`.
    std::fs::rename(incomplete_path, path)
}

/// Publish a summary only after every trajectory, CSV write, and checksum succeeds.
/// Failed exports remain explicitly named `.incomplete`, outside the analysis CSV glob.
fn finish_run_output(
    run_result: io::Result<()>,
    csv_path: &Path,
    export: impl FnOnce(&Path) -> io::Result<()>,
) -> RunOutput {
    let mut output = RunOutput {
        status: "simulation_failed",
        csv_path: None,
        summary_hash: None,
        failure_detail: None,
    };
    if let Err(error) = run_result {
        output.fail("simulation_failed", error);
        return output;
    }

    // Reserve a unique file, preserving any stale diagnostic files from prior attempts.
    let incomplete_path = match reserve_incomplete_file(csv_path) {
        Ok((path, file)) => {
            drop(file);
            path
        }
        Err(error) => {
            output.fail("csv_export_failed", error);
            return output;
        }
    };
    output.csv_path = Some(incomplete_path.clone());

    if let Err(error) = export(&incomplete_path) {
        output.fail("csv_export_failed", error);
        return output;
    }
    let hash = match hash_file_sha256(&incomplete_path) {
        Ok(hash) => hash,
        Err(error) => {
            output.fail("summary_hash_failed", error);
            return output;
        }
    };
    let published_path = match publish_summary_no_clobber(&incomplete_path, csv_path) {
        Ok(path) => path,
        Err(error) => {
            output.fail("csv_publish_failed", error);
            return output;
        }
    };
    output.status = "completed";
    output.csv_path = Some(published_path);
    output.summary_hash = Some(hash);
    output
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
    failure_detail: Option<&str>,
    last_timestep: Option<usize>,
    config_validation_mode: &str,
    config_validation_status: &str,
    config_validation_errors: usize,
    config_validation_warnings: usize,
    config_validation_report_path: Option<&std::path::Path>,
) -> std::io::Result<()> {
    write_atomic_metadata(path, |file| {
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
        writeln!(
            file,
            "failure_detail={}",
            failure_detail.unwrap_or("none").replace(['\r', '\n'], " ")
        )?;
        writeln!(file, "replay_env=AMR_RNG_SEED={}", seed.value)?;

        Ok(())
    })
}

fn main() -> ExitCode {
    install_panic_log_hook();
    configure_rayon_worker_stack();
    let _ = env_logger::builder().is_test(false).try_init();

    // Main run configuration. This is the quickest place to switch between calibration-sized
    // runs, full policy runs, deterministic debug runs, and journey-logging experiments.
    let population_size = 10_000_000;
    // CalibrationMode::FullMinimal — sparse 2022-2025 CSV with the lean core per-bacteria profile.
    // CalibrationMode::Full        — sparse 2022-2025 CSV with all fields needed for calibration_summary.txt.
    // CalibrationMode::Partial     — all 1930-2025 rows kept; time-series plots still work.
    // CalibrationMode::Partial25Counterfactual — full 1930-2025 baseline plus no-resistance 2022-2025.
    // CalibrationMode::Full25Counterfactual — sparse baseline and no-resistance rows for 2022-2025.
    // CalibrationMode::None        — full run with policy branches to 2035.
    let calibration_mode = CalibrationMode::Full;
    let time_steps = match calibration_mode {
        CalibrationMode::None => 38_325,
        CalibrationMode::Partial | CalibrationMode::FullMinimal | CalibrationMode::Full => 35_040,
        CalibrationMode::Partial25Counterfactual | CalibrationMode::Full25Counterfactual => 35_040,
    };
    debug_assert_eq!(time_steps, calibration_mode.time_steps());
    let log_individuals = false; // Log rich daily state for the first ten population records.
    let log_infection_journeys = false; // Log dense trajectories for sampled infected people.
    let infection_journey_sample_rate = 1.00; // Daily sampling probability for an untracked infected person; qualifying resistance bypasses it.
    let use_fixed_seed = false; // Use fixed_seed_value only when AMR_RNG_SEED is unset.
    let fixed_seed_value: u64 = 1_234_567_890; // Deterministic fallback selected by use_fixed_seed.
    let infection_journey_bacteria_filter: Option<&str> = None; // Optional filter on the selected primary bacterium.

    // The filter must exactly match a value in BACTERIA_LIST. Examples:
    // Some("escherichia_coli")
    // Some("staphylococcus_aureus")
    // Some("pseudomonas_aeruginosa")
    // Some("acinetobacter_baumannii")
    // Some("enterococcus_faecium")
    // None disables the primary-bacterium filter.

    let resolved_run_seed = resolve_run_seed(use_fixed_seed, fixed_seed_value);
    std::env::set_var("AMR_RNG_SEED", resolved_run_seed.value.to_string());
    let seed_override = Some(resolved_run_seed.value);
    let source_hash = resolve_source_hash();
    eprintln!("[startup] source_hash={}", source_hash);

    // ── Policy branch selection ────────────────────────────────────────────────
    // Mode-specific policy runs. Every branch-enabled mode runs the complete baseline first.
    // The 2025 counterfactual modes then run alternate policy 2; the full policy mode runs
    // alternates 1 through 4. ID 0 in the configured lists is redundant and filtered out.
    //
    //   0 = Baseline continuation      (status quo carried forward after the checkpoint)
    //   1 = Antimicrobial Stewardship  (reduced prescribing, better drug selection)
    //   2 = AMR Counterfactual         (resistance-suppressed comparison branch)
    //   3 = Near-complete Diagnostics  (high testing and more targeted selection)
    //   4 = Equal Global Access        (testing, initiation, and cessation use North America references)
    //
    // Branches are independent and restore the mode-specific checkpoint: 2022 for the
    // 2025 counterfactual modes and 2027 for the full policy mode.
    let active_policies = calibration_mode.active_policy_ids();

    let output_dir = std::path::Path::new("amr_simulation_output_analysis_outputs");
    let metadata_stamp = Utc::now().format("%Y%m%dT%H%M%S%.9fZ").to_string();
    let InvocationPaths {
        metadata_path,
        validation_path: config_validation_report_path,
        validation_file: mut config_validation_file,
    } = match reserve_invocation_paths(output_dir, resolved_run_seed.value, &metadata_stamp) {
        Ok(paths) => paths,
        Err(error) => {
            eprintln!("Unable to reserve output paths before simulation: {error}");
            return ExitCode::FAILURE;
        }
    };

    let config_validation_mode = match ConfigValidationMode::from_env() {
        Ok(mode) => mode,
        Err(err) => {
            eprintln!("[config-validation] FAILED: {}", err);
            std::process::exit(2);
        }
    };
    let config_validation_report = validate_parameter_map(&PARAMETERS);
    let rendered_config_validation = config_validation_report.render(config_validation_mode);
    eprint!("{}", rendered_config_validation);
    let config_validation_report_for_metadata = match config_validation_file
        .write_all(rendered_config_validation.as_bytes())
        .and_then(|()| config_validation_file.sync_all())
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
    drop(config_validation_file);

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
            Some("Configuration validation rejected the run before simulation"),
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

        println!(
            "[report] status=config_validation_failed run_id=pending source_hash={} rng_seed={} last_timestep=pending summary_hash=pending failure_class=config_validation_failed config_validation_status={} config_validation_mode={} summary_csv=pending",
            source_hash,
            resolved_run_seed.value,
            config_validation_report.status(),
            config_validation_mode
        );
        std::process::exit(2);
    }

    // Reject an empty bacterium roster and print the configured roster before the full run.
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
    let use_disk_branch_checkpointing = calibration_mode.uses_policy_branches();
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

    let run_result = simulation.run();

    let duration = start.elapsed();

    // Print journey-logging statistics and alternate-branch coverage.
    simulation.print_summary_statistics();

    // Keep the seed-derived run ID in the CSV. Publication adds a repeat suffix to
    // the filename if that ID already has a result, without changing simulation data.
    let run_id = simulation.run_id;
    let csv_basename = format!("simulation_summary_{:06}.csv", run_id);
    let csv_path = output_dir.join(&csv_basename);

    // A failed trajectory must never be handed to analysis as a completed summary.
    let mut output = finish_run_output(run_result, &csv_path, |path| {
        simulation.export_summary_to_csv(path)
    });
    let last_timestep = observability::current_timestep();

    if let Err(err) = write_run_metadata(
        &metadata_path,
        output.status,
        &source_hash,
        resolved_run_seed,
        population_size,
        time_steps,
        calibration_mode,
        active_policies,
        Some(run_id),
        output.csv_path.as_deref(),
        Some(duration.as_secs_f64()),
        output.summary_hash.as_deref(),
        output.status,
        output.failure_detail.as_deref(),
        last_timestep,
        &config_validation_mode.to_string(),
        config_validation_report.status(),
        config_validation_report.error_count(),
        config_validation_report.warning_count(),
        config_validation_report_for_metadata,
    ) {
        eprintln!(
            "Error: unable to update run metadata {}: {}",
            metadata_path.display(),
            err
        );
        // Preserve the primary failure if the simulation or output already failed.
        if output.status == "completed" {
            output.fail("metadata_write_failed", err);
        }
    }

    if let Some(detail) = &output.failure_detail {
        eprintln!("[run-failure] {}: {}", output.status, detail);
    }
    if output.status == "completed" {
        println!(
            "Summary data exported to {}",
            output.csv_path.as_ref().unwrap().display()
        );
    }

    println!(
        "[report] status={} run_id={} source_hash={} rng_seed={} last_timestep={} summary_hash={} failure_class={} config_validation_status={} config_validation_mode={} summary_csv={}",
        output.status,
        run_id,
        source_hash,
        resolved_run_seed.value,
        last_timestep
            .map(|step| step.to_string())
            .unwrap_or_else(|| "pending".to_string()),
        output.summary_hash.as_deref().unwrap_or("pending"),
        output.status,
        config_validation_report.status(),
        config_validation_mode,
        output.csv_path.as_deref().map(|path| path.display().to_string())
            .unwrap_or_else(|| "pending".to_string())
    );

    if let Err(e) = log_simulation_run(population_size, time_steps, duration.as_secs_f64()) {
        eprintln!("Error logging simulation run: {}", e);
    }

    println!("\n--- simulation ended ---");
    println!(
        "--- total simulation time: {:.3} seconds",
        duration.as_secs_f64()
    );
    println!("                          ");
    output.exit_code()
}

fn install_panic_log_hook() {
    std::panic::set_hook(Box::new(|panic_info| {
        let timestamp = Utc::now().format("%Y-%m-%d %H:%M:%S UTC");
        let rng_seed = std::env::var("AMR_RNG_SEED").unwrap_or_else(|_| "unset".to_string());
        let source_hash = resolve_source_hash();
        let run_id = observability::current_run_id()
            .map(|id| id.to_string())
            .unwrap_or_else(|| "unset".to_string());
        let last_timestep = observability::current_timestep()
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

    if !file_exists {
        writeln!(
            file,
            "timestamp,population_size,time_steps,duration_seconds"
        )?;
    }

    file.write_all(log_entry.as_bytes())?;

    println!("Simulation run logged to {}", log_path);

    Ok(())
}

/// Validate the bacteria roster before the simulation starts.
///
/// The executable currently uses the canonical roster in `BACTERIA_LIST`. Subset runs require
/// coordinated changes to aligned inventories, parameters, output contracts, and tests; this
/// function is only a startup guard and summary, not a complete subset-run validator.
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    struct TestDirectory(PathBuf);

    impl TestDirectory {
        fn new() -> Self {
            static NEXT_ID: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "amr_run_output_{}_{}_{}",
                std::process::id(),
                Utc::now().timestamp_nanos_opt().unwrap(),
                NEXT_ID.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }

        fn csv_path(&self) -> PathBuf {
            self.0.join("simulation_summary_123456.csv")
        }
    }

    impl Drop for TestDirectory {
        fn drop(&mut self) {
            // Only this uniquely created test directory is owned by this fixture.
            assert_eq!(self.0.parent(), Some(std::env::temp_dir().as_path()));
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn write_fixture_metadata(path: &Path, output: &RunOutput) {
        write_run_metadata(
            path,
            output.status,
            "test-source",
            ResolvedRunSeed {
                value: 123,
                source: "test",
            },
            2,
            3,
            CalibrationMode::Partial,
            &[0],
            Some(123456),
            output.csv_path.as_deref(),
            Some(0.1),
            output.summary_hash.as_deref(),
            output.status,
            output.failure_detail.as_deref(),
            None,
            "strict",
            "passed",
            0,
            0,
            None,
        )
        .unwrap();
    }

    #[test]
    fn failed_trajectory_does_not_export_or_replace_existing_summary() {
        let directory = TestDirectory::new();
        let csv_path = directory.csv_path();
        std::fs::write(&csv_path, "previous completed run").unwrap();
        let output = finish_run_output(
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "policy 2 checkpoint checksum mismatch",
            )),
            &csv_path,
            |_| panic!("failed trajectories must not export summaries"),
        );
        assert_eq!(output.status, "simulation_failed");
        assert_eq!(output.exit_code(), ExitCode::FAILURE);
        assert!(output.csv_path.is_none());
        assert!(output.summary_hash.is_none());
        assert_eq!(
            std::fs::read_to_string(&csv_path).unwrap(),
            "previous completed run"
        );
        let metadata_path = directory.0.join("metadata.txt");
        write_fixture_metadata(&metadata_path, &output);
        let metadata = std::fs::read_to_string(metadata_path).unwrap();
        assert!(metadata.contains("status=simulation_failed\n"));
        assert!(metadata.contains("failure_detail=policy 2 checkpoint checksum mismatch\n"));
        assert!(metadata.contains("summary_csv=pending\n"));
        assert!(metadata.contains("last_timestep=pending\n"));
        assert!(!metadata.contains("status=completed\n"));
    }

    #[test]
    fn failed_export_retains_only_an_explicitly_incomplete_file() {
        let directory = TestDirectory::new();
        let csv_path = directory.csv_path();
        let output = finish_run_output(Ok(()), &csv_path, |path| {
            std::fs::write(path, "partial CSV contents")?;
            Err(io::Error::other("injected write failure"))
        });
        assert_eq!(output.status, "csv_export_failed");
        assert_eq!(output.exit_code(), ExitCode::FAILURE);
        assert!(!csv_path.exists());
        let incomplete = output.csv_path.unwrap();
        assert_eq!(incomplete.extension().unwrap(), "incomplete");
        assert_eq!(
            std::fs::read_to_string(incomplete).unwrap(),
            "partial CSV contents"
        );
    }

    #[test]
    fn successful_export_publishes_the_hashed_file() {
        let directory = TestDirectory::new();
        let csv_path = directory.csv_path();
        let output = finish_run_output(Ok(()), &csv_path, |path| {
            std::fs::write(path, "time_step,policy_option\n0,0\n")
        });
        assert_eq!(output.status, "completed");
        assert_eq!(output.exit_code(), ExitCode::SUCCESS);
        assert_eq!(output.csv_path.as_deref(), Some(csv_path.as_path()));
        assert_eq!(
            output.summary_hash,
            Some(hash_file_sha256(&csv_path).unwrap())
        );
        assert!(output.failure_detail.is_none());
        assert_eq!(std::fs::read_dir(&directory.0).unwrap().count(), 1);
    }

    #[test]
    fn stale_incomplete_files_are_preserved_and_do_not_block_new_exports() {
        let directory = TestDirectory::new();
        let requested = directory.csv_path();
        let (stale_path, mut stale_file) = reserve_incomplete_file(&requested).unwrap();
        stale_file
            .write_all(b"previous incomplete attempt")
            .unwrap();
        drop(stale_file);
        let output = finish_run_output(Ok(()), &requested, |path| {
            std::fs::write(path, "new complete data")
        });
        assert_eq!(output.exit_code(), ExitCode::SUCCESS);
        assert_eq!(
            std::fs::read_to_string(stale_path).unwrap(),
            "previous incomplete attempt"
        );
        assert_eq!(
            std::fs::read_to_string(requested).unwrap(),
            "new complete data"
        );
    }

    #[test]
    fn repeat_publication_preserves_previous_csv_and_metadata() {
        let directory = TestDirectory::new();
        let csv_path = directory.csv_path();
        let first = finish_run_output(Ok(()), &csv_path, |path| {
            std::fs::write(path, "time_step,value\n0,10\n")
        });
        let metadata_path = directory.0.join("original_metadata.txt");
        write_fixture_metadata(&metadata_path, &first);
        let original_metadata = std::fs::read(&metadata_path).unwrap();
        let original_csv = std::fs::read(&csv_path).unwrap();

        let second = finish_run_output(Ok(()), &csv_path, |path| {
            std::fs::write(path, "time_step,value\n0,20\n")
        });
        let repeat_path = directory.0.join("simulation_summary_123456_repeat_1.csv");
        assert_eq!(second.exit_code(), ExitCode::SUCCESS);
        assert_eq!(second.csv_path.as_deref(), Some(repeat_path.as_path()));
        assert_eq!(std::fs::read(&csv_path).unwrap(), original_csv);
        assert_eq!(std::fs::read(&metadata_path).unwrap(), original_metadata);
        assert_eq!(
            first.summary_hash,
            Some(hash_file_sha256(&csv_path).unwrap())
        );
        assert_eq!(
            second.summary_hash,
            Some(hash_file_sha256(&repeat_path).unwrap())
        );
        assert_ne!(first.summary_hash, second.summary_hash);
        let repeat_metadata = directory.0.join("repeat_metadata.txt");
        write_fixture_metadata(&repeat_metadata, &second);
        assert!(std::fs::read_to_string(repeat_metadata)
            .unwrap()
            .contains(&format!("summary_csv={}\n", repeat_path.display())));
    }

    #[test]
    fn repeated_fixed_seed_runs_keep_identical_csv_bytes_in_distinct_files() {
        let directory = TestDirectory::new();
        let mut outputs = Vec::new();
        let mut run_ids = Vec::new();
        for _ in 0..2 {
            let mut simulation =
                Simulation::new(2, 2, false, Some(12345), CalibrationMode::Partial);
            simulation.summary_content_flags =
                amr_project::simulation::simulation::SummaryContentFlags::none();
            let run_result = simulation.run();
            run_ids.push(simulation.run_id);
            let requested = directory
                .0
                .join(format!("simulation_summary_{:06}.csv", simulation.run_id));
            let output = finish_run_output(run_result, &requested, |path| {
                simulation.export_summary_to_csv(path)
            });
            assert_eq!(output.exit_code(), ExitCode::SUCCESS);
            outputs.push(output);
        }
        assert_eq!(run_ids[0], run_ids[1]);
        assert_ne!(outputs[0].csv_path, outputs[1].csv_path);
        assert_eq!(outputs[0].summary_hash, outputs[1].summary_hash);
        assert_eq!(
            std::fs::read(outputs[0].csv_path.as_ref().unwrap()).unwrap(),
            std::fs::read(outputs[1].csv_path.as_ref().unwrap()).unwrap()
        );
    }

    #[test]
    fn identical_invocation_timestamps_never_reuse_diagnostic_paths() {
        let directory = TestDirectory::new();
        let mut first = reserve_invocation_paths(&directory.0, 123, "same_timestamp").unwrap();
        std::fs::write(&first.metadata_path, "original metadata").unwrap();
        first
            .validation_file
            .write_all(b"original validation")
            .unwrap();
        let second = reserve_invocation_paths(&directory.0, 123, "same_timestamp").unwrap();
        assert_ne!(first.metadata_path, second.metadata_path);
        assert_ne!(first.validation_path, second.validation_path);
        assert_eq!(
            std::fs::read_to_string(&first.metadata_path).unwrap(),
            "original metadata"
        );
        assert_eq!(
            std::fs::read_to_string(&first.validation_path).unwrap(),
            "original validation"
        );
    }

    #[test]
    fn orphan_validation_report_is_preserved_during_path_reservation() {
        let directory = TestDirectory::new();
        let first = reserve_invocation_paths(&directory.0, 123, "same_timestamp").unwrap();
        drop(first.validation_file);
        std::fs::write(&first.validation_path, "orphan validation").unwrap();
        std::fs::remove_file(&first.metadata_path).unwrap();

        let second = reserve_invocation_paths(&directory.0, 123, "same_timestamp").unwrap();
        assert_ne!(second.validation_path, first.validation_path);
        assert!(!first.metadata_path.exists());
        assert_eq!(
            std::fs::read_to_string(first.validation_path).unwrap(),
            "orphan validation"
        );
    }

    #[test]
    fn concurrent_invocations_reserve_distinct_metadata_and_validation_files() {
        let directory = TestDirectory::new();
        let barrier = std::sync::Arc::new(std::sync::Barrier::new(4));
        let threads: Vec<_> = (0..4)
            .map(|_| {
                let output_dir = directory.0.clone();
                let barrier = barrier.clone();
                std::thread::spawn(move || {
                    barrier.wait();
                    reserve_invocation_paths(&output_dir, 123, "same_timestamp").unwrap()
                })
            })
            .collect();
        let paths: Vec<_> = threads
            .into_iter()
            .map(|thread| thread.join().unwrap())
            .collect();
        let unique: std::collections::HashSet<_> = paths
            .iter()
            .flat_map(|paths| [&paths.metadata_path, &paths.validation_path])
            .collect();
        assert_eq!(unique.len(), 8);
    }

    #[test]
    fn missing_export_cannot_be_hashed_or_reported_as_completed() {
        let directory = TestDirectory::new();
        let csv_path = directory.csv_path();
        let output = finish_run_output(Ok(()), &csv_path, |path| std::fs::remove_file(path));
        assert_eq!(output.status, "summary_hash_failed");
        assert_eq!(output.exit_code(), ExitCode::FAILURE);
        assert!(!csv_path.exists());
        assert!(output.summary_hash.is_none());
    }

    #[test]
    fn failed_publication_retains_the_incomplete_file_and_returns_failure() {
        let directory = TestDirectory::new();
        let csv_path = directory.csv_path();
        std::fs::create_dir(&csv_path).unwrap();
        let output = finish_run_output(Ok(()), &csv_path, |path| {
            std::fs::write(path, "finished CSV")
        });
        assert_eq!(output.status, "csv_publish_failed");
        assert_eq!(output.exit_code(), ExitCode::FAILURE);
        assert!(csv_path.is_dir());
        assert_eq!(
            std::fs::read_to_string(output.csv_path.unwrap()).unwrap(),
            "finished CSV"
        );
    }

    #[test]
    fn metadata_write_failure_preserves_previous_status_and_reports_failure() {
        let directory = TestDirectory::new();
        let metadata_path = directory.0.join("metadata.txt");
        std::fs::write(&metadata_path, "status=started\n").unwrap();
        let error = write_atomic_metadata(&metadata_path, |file| {
            writeln!(file, "status=completed")?;
            Err(io::Error::other("injected metadata write failure"))
        })
        .expect_err("partial metadata must never replace the current record");
        assert_eq!(
            std::fs::read_to_string(&metadata_path).unwrap(),
            "status=started\n"
        );

        let mut output = RunOutput {
            status: "completed",
            csv_path: None,
            summary_hash: None,
            failure_detail: None,
        };
        output.fail("metadata_write_failed", error);
        assert_eq!(output.exit_code(), ExitCode::FAILURE);
        assert!(output
            .failure_detail
            .unwrap()
            .contains("metadata write failure"));

        write_fixture_metadata(
            &metadata_path,
            &RunOutput {
                status: "simulation_failed",
                csv_path: None,
                summary_hash: None,
                failure_detail: Some("checkpoint restore failed".to_string()),
            },
        );
        assert!(std::fs::read_to_string(metadata_path)
            .unwrap()
            .contains("status=simulation_failed\n"));
    }

    // Run in a child test process below to exercise the same ExitCode returned by main.
    #[test]
    fn launcher_exit_probe() -> ExitCode {
        let Ok(case) = std::env::var("AMR_TEST_RUN_OUTPUT_CASE") else {
            return ExitCode::SUCCESS;
        };
        let directory = PathBuf::from(std::env::var_os("AMR_TEST_RUN_OUTPUT_DIR").unwrap());
        let run_result = if case == "simulation_failed" {
            Err(io::Error::other("checkpoint restore failed"))
        } else {
            Ok(())
        };
        let output = finish_run_output(run_result, &directory.join("summary.csv"), |path| {
            std::fs::write(path, "time_step\n0\n")?;
            if case == "csv_export_failed" {
                Err(io::Error::other("summary write failed"))
            } else {
                Ok(())
            }
        });
        write_fixture_metadata(&directory.join("metadata.txt"), &output);
        output.exit_code()
    }

    #[test]
    fn launcher_exit_status_and_metadata_agree_on_success_and_failures() {
        for case in ["completed", "simulation_failed", "csv_export_failed"] {
            let directory = TestDirectory::new();
            let child = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "tests::launcher_exit_probe", "--nocapture"])
                .env("AMR_TEST_RUN_OUTPUT_CASE", case)
                .env("AMR_TEST_RUN_OUTPUT_DIR", &directory.0)
                .output()
                .unwrap();
            assert_eq!(
                child.status.success(),
                case == "completed",
                "{case}: {:?}",
                child
            );
            let metadata = std::fs::read_to_string(directory.0.join("metadata.txt")).unwrap();
            assert!(metadata.contains(&format!("status={case}\n")), "{metadata}");
            assert_eq!(
                directory.0.join("summary.csv").exists(),
                case == "completed"
            );
        }
    }
}
