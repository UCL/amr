use amr_project::simulation::simulation::{CalibrationMode, Simulation, SummaryContentFlags};
use std::env;
use std::time::Instant;

#[derive(Clone, Copy)]
struct Scenario {
    population_size: usize,
    time_steps: usize,
}

fn parse_scenarios(value: &str) -> Vec<Scenario> {
    value
        .split(',')
        .filter_map(|part| {
            let (population, steps) = part.split_once(':')?;
            Some(Scenario {
                population_size: population.trim().parse().ok()?,
                time_steps: steps.trim().parse().ok()?,
            })
        })
        .collect()
}

fn env_usize(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn env_u64(name: &str, default: u64) -> u64 {
    env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn calibration_mode_from_env() -> CalibrationMode {
    match env::var("PERF_MODE")
        .unwrap_or_else(|_| "partial".to_string())
        .to_ascii_lowercase()
        .as_str()
    {
        "full" => CalibrationMode::Full,
        "fullminimal" | "full_minimal" | "minimal" => CalibrationMode::FullMinimal,
        "none" => CalibrationMode::None,
        _ => CalibrationMode::Partial,
    }
}

fn main() {
    let label = env::var("PERF_LABEL").unwrap_or_else(|_| "unknown".to_string());
    let scenarios = parse_scenarios(
        &env::var("PERF_SCENARIOS").unwrap_or_else(|_| "64:4,256:8,1024:16".to_string()),
    );
    let repeats = env_usize("PERF_REPEATS", 3);
    let seed_base = env_u64("PERF_SEED_BASE", 987_654_321);
    let lean_summary = env::var("PERF_SUMMARY")
        .map(|value| value.eq_ignore_ascii_case("none"))
        .unwrap_or(false);
    let calibration_mode = calibration_mode_from_env();

    println!(
        "PERF_HEADER,label,population_size,time_steps,repeat,seed,mode,summary_mode,init_seconds,run_seconds,total_seconds,summary_rows"
    );

    for scenario in scenarios {
        for repeat in 0..repeats {
            let seed = seed_base + repeat as u64;
            let total_start = Instant::now();
            let init_start = Instant::now();
            let mut simulation = Simulation::new(
                scenario.population_size,
                scenario.time_steps,
                false,
                Some(seed),
                calibration_mode,
            );
            if lean_summary {
                simulation.summary_content_flags = SummaryContentFlags::none();
            }
            let init_seconds = init_start.elapsed().as_secs_f64();

            let run_start = Instant::now();
            simulation.run();
            let run_seconds = run_start.elapsed().as_secs_f64();
            let total_seconds = total_start.elapsed().as_secs_f64();
            let summary_rows = simulation.summary_log.len();
            let summary_mode = if lean_summary { "none" } else { "default" };

            println!(
                "PERF_RESULT,{label},{},{},{repeat},{seed},{calibration_mode:?},{summary_mode},{init_seconds:.6},{run_seconds:.6},{total_seconds:.6},{summary_rows}",
                scenario.population_size,
                scenario.time_steps
            );
        }
    }
}
