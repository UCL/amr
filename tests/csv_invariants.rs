use amr_project::simulation::journey_logger::JourneyLogger;
use amr_project::simulation::population::Individual;
use amr_project::simulation::simulation::{CalibrationMode, Simulation};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

static JOURNEY_LOGGER_FILE_LOCK: Mutex<()> = Mutex::new(());

fn output_path(label: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("csv_invariant_tests");
    let _ = fs::create_dir_all(&path);
    path.push(format!("{}_{}.csv", label, std::process::id()));
    path
}

fn assert_csv_rows_match_header_width(path: &Path) {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .expect("CSV should open");
    let mut records = reader.records();
    let header = records
        .next()
        .expect("CSV should include a header row")
        .expect("header row should parse");
    let expected_width = header.len();

    assert!(expected_width > 0, "CSV header should not be empty");

    for (row_idx, record) in records.enumerate() {
        let record = record.expect("CSV data row should parse");
        assert_eq!(
            record.len(),
            expected_width,
            "CSV row {} should have the same width as the header",
            row_idx + 2
        );
    }
}

fn csv_header(path: &Path) -> csv::StringRecord {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .expect("CSV should open");
    reader
        .records()
        .next()
        .expect("CSV should include a header row")
        .expect("header row should parse")
}

fn journey_output_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("infection_journeys");
    path.push("infection_journeys.csv");
    path
}

#[test]
fn summary_csv_rows_match_header_width_for_tiny_run() {
    let mut simulation = Simulation::new(64, 4, false, Some(246_813_579), CalibrationMode::Partial);
    simulation.run();

    let path = output_path("summary_schema");
    simulation
        .export_summary_to_csv(&path)
        .expect("summary export should succeed");

    assert_csv_rows_match_header_width(&path);
}

#[test]
fn journey_csv_rows_match_header_width_for_forced_infection() {
    let _guard = JOURNEY_LOGGER_FILE_LOCK
        .lock()
        .expect("journey logger file lock should not be poisoned");

    let path = journey_output_path();
    let _ = fs::remove_file(&path);

    let mut rng = SmallRng::seed_from_u64(246_813_579);
    let mut individual = Individual::new(42, 45 * 365, "female".to_string(), &mut rng);
    let bacteria_idx = 0;
    let drug_idx = 0;

    individual.level[bacteria_idx] = 2.5;
    individual.date_last_infected[bacteria_idx] = 1;
    individual.date_last_infected_keep[bacteria_idx] = 1;
    individual.infectious_syndrome[bacteria_idx] = 0;
    individual.infection_has_caused_symptoms[bacteria_idx] = true;
    individual.presence_microbiome[bacteria_idx] = true;
    individual.cur_use_drug[drug_idx] = true;
    individual.cur_level_drug[drug_idx] = 5.0;
    individual.date_drug_initiated[drug_idx] = 1;
    individual.days_on_current_treatment[bacteria_idx] = 1;
    individual.resistances[bacteria_idx][drug_idx].any_r = 0.25;
    individual.resistances[bacteria_idx][drug_idx].activity_r = 0.75;
    individual.resistances[bacteria_idx][drug_idx].microbiome_r = 0.50;

    let mut logger = JourneyLogger::new(Some(123_456_789));
    logger
        .enable(1.0)
        .expect("journey logger should create its CSV");
    logger.check_individual(&individual, 1);
    logger.finalize().expect("journey logger should flush");

    let (_, journeys_started, snapshots_logged) = logger.get_stats();
    assert_eq!(
        journeys_started, 1,
        "forced infection should start a journey"
    );
    assert_eq!(
        snapshots_logged, 1,
        "forced infection should write one snapshot"
    );

    logger.close().expect("journey logger should close");

    assert_csv_rows_match_header_width(&path);
    let header = csv_header(&path);
    assert!(
        !header
            .iter()
            .any(|column| column == "resistance_majority_r"),
        "journey CSV should not export the removed majority_r scalar"
    );
    let _ = fs::remove_file(&path);
}
