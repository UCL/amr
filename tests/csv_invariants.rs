use amr_project::simulation::simulation::{CalibrationMode, Simulation};
use std::fs;
use std::path::{Path, PathBuf};

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
