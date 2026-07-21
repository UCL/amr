use amr_project::config::PARAMETER_STORE;
use amr_project::rules::ParameterKeyCache;
use amr_project::simulation::population::{
    ResistanceMechanism, BACTERIA_LIST, DRUG_CLASS_LOOKUP, DRUG_SHORT_NAMES,
};
use csv::ReaderBuilder;
use std::collections::BTreeSet;
use std::path::PathBuf;

#[test]
fn exported_resistance_reachability_matches_the_typed_rust_model() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("model_resistance_reachability_matrix.csv");
    let mut reader = ReaderBuilder::new()
        .from_path(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let projected_headers: Vec<&str> = reader
        .headers()
        .expect("resistance reachability projection header")
        .iter()
        .collect();
    assert_eq!(
        projected_headers,
        [
            "bacteria",
            "drug",
            "resistance_representable",
            "maximum_any_r",
            "positive_effect_mechanisms",
        ]
    );

    let cache = ParameterKeyCache::new();
    let mut seen = BTreeSet::new();
    for record in reader.records() {
        let record = record.expect("valid resistance reachability projection row");
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
            "duplicate resistance reachability row for {bacterium}/{drug}"
        );

        let expected_mechanisms = ResistanceMechanism::all()
            .iter()
            .enumerate()
            .filter_map(|(mechanism_idx, mechanism)| {
                let effect = PARAMETER_STORE
                    .resistance_mechanism
                    .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_idx]);
                (cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx) && effect > 0.0)
                    .then_some(mechanism.as_str())
            })
            .collect::<Vec<_>>();
        assert_eq!(
            &record[2],
            if expected_mechanisms.is_empty() {
                "false"
            } else {
                "true"
            },
            "projected resistance reachability drift for {bacterium}/{drug}"
        );
        assert_eq!(
            record[3]
                .parse::<f64>()
                .unwrap_or_else(|_| panic!("invalid maximum any_r for {bacterium}/{drug}")),
            ResistanceMechanism::all()
                .iter()
                .enumerate()
                .filter_map(|(mechanism_idx, _)| {
                    let effect = PARAMETER_STORE
                        .resistance_mechanism
                        .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_idx]);
                    (cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                        && effect > 0.0)
                        .then_some(effect)
                })
                .fold(0.0, |combined, effect| {
                    1.0 - (1.0 - combined) * (1.0 - effect)
                }),
            "projected maximum any_r drift for {bacterium}/{drug}"
        );
        assert_eq!(
            &record[4],
            expected_mechanisms.join(";"),
            "projected mechanism list drift for {bacterium}/{drug}"
        );
    }

    assert_eq!(
        seen.len(),
        BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len(),
        "resistance reachability projection must contain every typed bacterium-drug pair"
    );
}
