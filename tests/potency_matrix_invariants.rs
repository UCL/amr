use amr_project::config::PARAMETER_STORE;
use amr_project::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES};
use csv::ReaderBuilder;
use std::collections::BTreeSet;
use std::path::PathBuf;

fn potency(bacterium: &str, drug: &str) -> f64 {
    let bacteria_idx = BACTERIA_LIST
        .iter()
        .position(|name| *name == bacterium)
        .unwrap_or_else(|| panic!("unknown bacterium {bacterium}"));
    let drug_idx = DRUG_SHORT_NAMES
        .iter()
        .position(|name| *name == drug)
        .unwrap_or_else(|| panic!("unknown drug {drug}"));
    PARAMETER_STORE
        .drug_bacteria
        .potency(bacteria_idx, drug_idx)
}

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

#[test]
fn staphylococcus_potencies_preserve_the_authoritative_matrix_entries() {
    let source = include_str!("../src/config.rs");
    for bacterium in ["staphylococcus_aureus", "staphylococcus_epidermidis"] {
        for drug in DRUG_SHORT_NAMES {
            let prefix = format!(
                "map.insert(\"drug_{drug}_for_bacteria_{bacterium}_potency_when_no_r\".to_string(), "
            );
            let entries: Vec<_> = source
                .lines()
                .filter_map(|line| line.trim().strip_prefix(prefix.as_str()))
                .collect();
            assert_eq!(
                entries.len(),
                1,
                "expected one matrix entry for {bacterium}/{drug}"
            );
            let configured: f64 = entries[0]
                .split_once(");")
                .expect("matrix entry should terminate after its numeric value")
                .0
                .trim()
                .parse()
                .expect("matrix potency should be numeric");
            assert_eq!(
                potency(bacterium, drug),
                configured.clamp(0.0, 1.0),
                "a later assignment must not override the {bacterium} matrix entry for {drug}"
            );
        }
    }
}

#[test]
fn ceftolozane_tazobactam_potencies_match_reviewed_spectrum() {
    let reviewed = [
        ("staphylococcus_aureus", 0.10),
        ("staphylococcus_epidermidis", 0.10),
        ("campylobacter_jejuni", 0.10),
        ("legionella_pneumophila", 0.05),
        ("streptococcus_pneumoniae", 0.75),
        ("streptococcus_pyogenes", 0.75),
        ("streptococcus_agalactiae", 0.75),
        ("bacteroides_fragilis", 0.45),
        ("pseudomonas_aeruginosa", 0.95),
    ];

    for (bacterium, expected) in reviewed {
        assert_eq!(
            potency(bacterium, "ceftolozane_tazobactam"),
            expected,
            "reviewed ceftolozane/tazobactam potency drift for {bacterium}"
        );
    }
}

#[test]
fn cefiderocol_potencies_match_reviewed_aerobic_gram_negative_spectrum() {
    let reviewed = [
        ("staphylococcus_aureus", 0.00),
        ("staphylococcus_epidermidis", 0.00),
        ("streptococcus_pneumoniae", 0.00),
        ("streptococcus_pyogenes", 0.00),
        ("streptococcus_agalactiae", 0.00),
        ("bacteroides_fragilis", 0.00),
        ("campylobacter_jejuni", 0.10),
        ("legionella_pneumophila", 0.05),
        ("stenotrophomonas_maltophilia", 0.95),
        ("burkholderia_cepacia_complex", 0.95),
        ("acinetobacter_baumannii", 0.90),
        ("pseudomonas_aeruginosa", 0.90),
    ];

    for (bacterium, expected) in reviewed {
        assert_eq!(
            potency(bacterium, "cefiderocol"),
            expected,
            "reviewed cefiderocol potency drift for {bacterium}"
        );
    }
}
