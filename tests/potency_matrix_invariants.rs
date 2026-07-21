use amr_project::config::PARAMETER_STORE;
use amr_project::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES};
use csv::ReaderBuilder;
use std::collections::BTreeSet;
use std::path::PathBuf;

#[test]
fn exported_model_potencies_match_the_typed_rust_matrix() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("model_potency_matrix.csv");
    let mut reader = ReaderBuilder::new()
        .from_path(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let projected_headers: Vec<&str> = reader
        .headers()
        .expect("potency projection header")
        .iter()
        .collect();
    assert_eq!(projected_headers, ["bacteria", "drug", "potency_when_no_r"]);

    let mut seen = BTreeSet::new();
    for record in reader.records() {
        let record = record.expect("valid potency projection row");
        let bacterium = &record[0];
        let drug = &record[1];
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|name| *name == bacterium)
            .unwrap_or_else(|| panic!("unknown projected bacterium {bacterium}"));
        let drug_idx = DRUG_SHORT_NAMES
            .iter()
            .position(|name| *name == drug)
            .unwrap_or_else(|| panic!("unknown projected drug {drug}"));
        assert!(
            seen.insert((bacteria_idx, drug_idx)),
            "duplicate projected potency for {bacterium}/{drug}"
        );

        let projected = record[2]
            .parse::<f64>()
            .unwrap_or_else(|_| panic!("invalid projected potency for {bacterium}/{drug}"));
        let typed = PARAMETER_STORE
            .drug_bacteria
            .potency(bacteria_idx, drug_idx);
        assert_eq!(
            projected, typed,
            "projected potency drift for {bacterium}/{drug}"
        );
    }

    assert_eq!(
        seen.len(),
        BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len(),
        "potency projection must contain every typed bacterium-drug pair"
    );
}
