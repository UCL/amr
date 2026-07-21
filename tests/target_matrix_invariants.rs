use amr_project::config::PARAMETER_STORE;
use amr_project::rules::ParameterKeyCache;
use amr_project::simulation::population::{
    ResistanceMechanism, BACTERIA_LIST, DRUG_CLASS_LOOKUP, DRUG_SHORT_NAMES,
};
use csv::ReaderBuilder;
use std::collections::BTreeSet;
use std::path::PathBuf;

const EXPECTED_TARGET_BACTERIA_COUNT: usize = 42;
const EXPECTED_TARGET_DRUG_COUNT: usize = 61;

// Nalidixic acid is retained in the model as a historical selection proxy, not
// as a current surveillance drug with a 2025 resistance target.
const CURRENT_TARGET_EXCLUSIONS: &[&str] = &["nalidixic_acid"];

struct TargetMatrix {
    bacteria: BTreeSet<String>,
    drugs: BTreeSet<String>,
}

fn read_target_matrix(file_name: &str, expects_notes: bool) -> TargetMatrix {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join(file_name);
    let mut reader = ReaderBuilder::new()
        .flexible(false)
        .from_path(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let headers = reader
        .headers()
        .unwrap_or_else(|error| panic!("failed to read headers from {}: {error}", path.display()))
        .clone();

    assert_eq!(
        headers
            .get(0)
            .map(|header| header.trim_start_matches('\u{feff}')),
        Some("Bacteria"),
        "{} must begin with a Bacteria column",
        path.display()
    );

    let data_end = if expects_notes {
        assert_eq!(
            headers.iter().next_back(),
            Some("notes"),
            "{} must end with a notes column",
            path.display()
        );
        headers.len() - 1
    } else {
        assert!(
            !headers.iter().any(|header| header == "notes"),
            "{} must not contain a notes column",
            path.display()
        );
        headers.len()
    };

    let mut drugs = BTreeSet::new();
    let drug_names: Vec<&str> = headers.iter().take(data_end).skip(1).collect();
    for drug in &drug_names {
        assert!(
            !drug.trim().is_empty(),
            "{} has an empty drug header",
            path.display()
        );
        assert!(
            drugs.insert((*drug).to_string()),
            "{} has duplicate drug header {drug}",
            path.display()
        );
    }
    assert_eq!(
        drugs.len(),
        EXPECTED_TARGET_DRUG_COUNT,
        "{} must contain {EXPECTED_TARGET_DRUG_COUNT} target drugs",
        path.display()
    );

    let mut bacteria = BTreeSet::new();
    for (row_offset, record) in reader.records().enumerate() {
        let row_number = row_offset + 2;
        let record = record.unwrap_or_else(|error| {
            panic!(
                "failed to parse row {row_number} in {}: {error}",
                path.display()
            )
        });
        assert_eq!(
            record.len(),
            headers.len(),
            "row {row_number} in {} has the wrong number of columns",
            path.display()
        );

        let bacterium = record[0].trim();
        assert!(
            !bacterium.is_empty(),
            "row {row_number} in {} has an empty bacterium",
            path.display()
        );
        assert!(
            bacteria.insert(bacterium.to_string()),
            "{} has duplicate bacterium {bacterium}",
            path.display()
        );

        for (drug, cell) in drug_names.iter().zip(record.iter().take(data_end).skip(1)) {
            let cell = cell.trim();
            if cell == "." {
                continue;
            }
            let value = cell.parse::<f64>().unwrap_or_else(|_| {
                panic!(
                    "{} row {row_number}, {bacterium}/{drug}: expected '.' or a number in [0, 1], found {cell:?}",
                    path.display()
                )
            });
            assert!(
                value.is_finite() && (0.0..=1.0).contains(&value),
                "{} row {row_number}, {bacterium}/{drug}: target {value} is outside [0, 1]",
                path.display()
            );
        }

        if expects_notes {
            assert!(
                !record[headers.len() - 1].trim().is_empty(),
                "{} row {row_number}, {bacterium}: provenance notes must not be empty",
                path.display()
            );
        }
    }

    assert_eq!(
        bacteria.len(),
        EXPECTED_TARGET_BACTERIA_COUNT,
        "{} must contain {EXPECTED_TARGET_BACTERIA_COUNT} bacteria",
        path.display()
    );

    TargetMatrix { bacteria, drugs }
}

fn target_name_to_model_slug(name: &str) -> String {
    match name {
        "Providencia stuartii" => "p_stuartii".to_string(),
        _ => name.to_ascii_lowercase().replace(' ', "_"),
    }
}

#[test]
fn resistance_target_matrices_are_valid_and_match_model_dimensions() {
    let prevalence = read_target_matrix("resistance_prevalence_values.csv", true);
    let severity = read_target_matrix("resistance_average_resistant_values.csv", false);

    assert_eq!(
        prevalence.bacteria, severity.bacteria,
        "prevalence and resistance-severity matrices must cover the same bacteria"
    );
    assert_eq!(
        prevalence.drugs, severity.drugs,
        "prevalence and resistance-severity matrices must cover the same drugs"
    );

    let target_bacteria: BTreeSet<String> = prevalence
        .bacteria
        .iter()
        .map(|name| target_name_to_model_slug(name))
        .collect();
    assert_eq!(
        target_bacteria.len(),
        prevalence.bacteria.len(),
        "target bacterium display names must map uniquely to model slugs"
    );
    let model_bacteria: BTreeSet<String> = BACTERIA_LIST
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    assert_eq!(
        target_bacteria, model_bacteria,
        "target matrices must cover every model bacterium exactly once"
    );

    for excluded in CURRENT_TARGET_EXCLUSIONS {
        assert!(
            DRUG_SHORT_NAMES.contains(excluded),
            "target exclusion {excluded} is not a model drug"
        );
    }
    let current_target_drugs: BTreeSet<String> = DRUG_SHORT_NAMES
        .iter()
        .filter(|drug| !CURRENT_TARGET_EXCLUSIONS.contains(drug))
        .map(|drug| (*drug).to_string())
        .collect();
    assert_eq!(
        prevalence.drugs, current_target_drugs,
        "target matrices must cover every current model drug except documented exclusions"
    );
}

#[test]
fn model_informed_target_exclusions_match_the_typed_rust_matrices() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("resistance_targets_v1.csv");
    let mut reader = ReaderBuilder::new()
        .from_path(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let headers = reader.headers().expect("long-form target header").clone();
    let column = |name: &str| {
        headers
            .iter()
            .position(|header| header == name)
            .unwrap_or_else(|| panic!("long-form targets lack {name}"))
    };
    let component_col = column("component");
    let bacteria_col = column("bacteria");
    let drug_col = column("drug");
    let value_col = column("value");
    let included_col = column("include_in_score");
    let reason_col = column("score_exclusion_reason");

    let cache = ParameterKeyCache::new();
    let mut low_potency_numeric_targets = 0;
    let mut unrepresentable_numeric_targets = 0;
    for record in reader.records() {
        let record = record.expect("valid long-form target row");
        if &record[component_col] != "resistance_prevalence_any_r_positive" {
            continue;
        }

        let model_slug = target_name_to_model_slug(&record[bacteria_col]);
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|name| *name == model_slug)
            .unwrap_or_else(|| panic!("unknown target bacterium {model_slug}"));
        let drug = &record[drug_col];
        let drug_idx = DRUG_SHORT_NAMES
            .iter()
            .position(|name| *name == drug)
            .unwrap_or_else(|| panic!("unknown target drug {drug}"));
        let has_value = !record[value_col].is_empty();
        let low_potency = PARAMETER_STORE
            .drug_bacteria
            .potency(bacteria_idx, drug_idx)
            < 0.15;
        let records_potency_exclusion = record[reason_col]
            .split(';')
            .any(|reason| reason == "model_baseline_potency_below_0.15");
        assert_eq!(
            records_potency_exclusion, low_potency,
            "potency exclusion drift for {model_slug}/{drug}"
        );
        let resistance_representable =
            ResistanceMechanism::all()
                .iter()
                .enumerate()
                .any(|(mechanism_idx, _)| {
                    cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                        && PARAMETER_STORE
                            .resistance_mechanism
                            .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_idx])
                            > 0.0
                });
        let records_reachability_exclusion = record[reason_col]
            .split(';')
            .any(|reason| reason == "model_resistance_phenotype_not_representable");
        assert_eq!(
            records_reachability_exclusion,
            has_value && !resistance_representable,
            "resistance representability exclusion drift for {model_slug}/{drug}"
        );

        if has_value && low_potency {
            low_potency_numeric_targets += 1;
        }
        if has_value && !resistance_representable {
            unrepresentable_numeric_targets += 1;
        }
        let expected_included = has_value
            && drug != "rifampicin"
            && !model_slug.contains("tuberculosis")
            && !model_slug.contains("listeria")
            && !low_potency
            && resistance_representable;
        assert_eq!(
            &record[included_col],
            if expected_included { "true" } else { "false" },
            "static inclusion drift for {model_slug}/{drug}"
        );
    }

    assert_eq!(
        low_potency_numeric_targets, 9,
        "review low-potency numeric targets when this count changes"
    );
    assert_eq!(
        unrepresentable_numeric_targets, 64,
        "review unrepresentable numeric targets when this count changes"
    );
}
