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
use amr_project::simulation::population::BACTERIA_LIST;
use amr_project::simulation::simulation::CalibrationMode;
use amr_project::simulation::simulation::Simulation;
use chrono::Utc;
use sha2::{Digest, Sha256};
use std::backtrace::Backtrace;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::path::PathBuf;

const RAYON_WORKER_STACK_BYTES: usize = 4 * 1024 * 1024;

#[derive(Clone, Copy)]
struct ResolvedRunSeed {
    value: u64,
    source: &'static str,
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

    // Main run configuration. This is the quickest place to switch between calibration-sized
    // runs, full policy runs, deterministic debug runs, and journey-logging experiments.
    let population_size =  10_000_000;
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
    if let Err(err) = std::fs::create_dir_all(output_dir) {
        eprintln!(
            "Warning: unable to create output directory {:?}: {}",
            output_dir, err
        );
    }

    let metadata_stamp = Utc::now().format("%Y%m%dT%H%M%SZ");
    let metadata_path = output_dir.join(format!(
        "run_metadata_{}_seed_{}.txt",
        metadata_stamp, resolved_run_seed.value
    ));

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

    simulation.run();

    let duration = start.elapsed();

    // Print journey-logging statistics and alternate-branch coverage.
    simulation.print_summary_statistics();

    // Include the pseudo-random run ID in the filename and metadata to associate the summary
    // with its run. The one-million-value ID space does not provide global uniqueness.
    let run_id = simulation.run_id;
    let csv_basename = format!("simulation_summary_{:06}.csv", run_id);
    let csv_path = output_dir.join(&csv_basename);

    // The summary CSV is the primary handoff to the Python analysis scripts.
    let (summary_hash, final_status, failure_class) =
        match simulation.export_summary_to_csv(&csv_path) {
            Ok(()) => {
                println!("Summary data exported to {}", csv_path.display());
                match hash_file_sha256(&csv_path) {
                    Ok(hash) => (Some(hash), "completed", "completed"),
                    Err(err) => {
                        eprintln!(
                            "Warning: unable to hash summary CSV {}: {}",
                            csv_path.display(),
                            err
                        );
                        (
                            None,
                            "completed_with_summary_hash_error",
                            "summary_hash_failed",
                        )
                    }
                }
            }
            Err(err) => {
                println!("Error exporting CSV: {}", err);
                (None, "csv_export_failed", "csv_export_failed")
            }
        };
    let last_timestep = observability::current_timestep().or_else(|| time_steps.checked_sub(1));

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
        Some(&csv_path),
        Some(duration.as_secs_f64()),
        summary_hash.as_deref(),
        failure_class,
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

    if let Err(e) = log_simulation_run(population_size, time_steps, duration.as_secs_f64()) {
        eprintln!("Error logging simulation run: {}", e);
    }

    println!("\n--- simulation ended ---");
    println!(
        "--- total simulation time: {:.3} seconds",
        duration.as_secs_f64()
    );
    println!("                          ");
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
