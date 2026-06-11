use amr_project::simulation::simulation::{CalibrationMode, Simulation, SummaryContentFlags};
use std::fs;
use std::path::PathBuf;

const SIMULATION_RS: &str = include_str!("../src/simulation/simulation.rs");

fn output_path(label: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("determinism_tests");
    let _ = fs::create_dir_all(&path);
    path.push(format!(
        "{}_{}_{}.csv",
        label,
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));
    path
}

fn summary_csv_for_thread_count(threads: usize) -> String {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("test Rayon pool should build");

    pool.install(|| {
        let mut simulation =
            Simulation::new(64, 4, false, Some(987_654_321), CalibrationMode::Partial);
        simulation.summary_content_flags = SummaryContentFlags::none();
        simulation.run();

        let output_path = output_path(&format!("threads_{threads}"));
        simulation
            .export_summary_to_csv(&output_path)
            .expect("summary export should not return an IO error");
        fs::read_to_string(&output_path).expect("exported CSV should be readable")
    })
}

#[test]
fn fixed_seed_rng_does_not_depend_on_rayon_thread_index() {
    assert!(
        !SIMULATION_RS.contains("current_thread_index"),
        "fixed-seed RNG derivation must not depend on Rayon worker assignment"
    );
}

#[test]
fn fixed_seed_summary_is_independent_of_rayon_thread_count() {
    let one_thread = summary_csv_for_thread_count(1);
    let four_threads = summary_csv_for_thread_count(4);

    assert_eq!(
        one_thread, four_threads,
        "fixed-seed summary output should not depend on Rayon worker count"
    );
}

#[test]
fn population_chunk_totals_are_heap_backed() {
    assert!(
        SIMULATION_RS.contains("Vec<Box<LocalTotals>>"),
        "population chunk totals should stay heap-backed to avoid large by-value Rayon accumulators"
    );
    assert!(
        SIMULATION_RS.contains("par_chunks_mut(DETERMINISTIC_POPULATION_CHUNK_SIZE)"),
        "population updates should use fixed-size deterministic chunks"
    );
}
