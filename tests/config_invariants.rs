use amr_project::config::{parameter_store, BacteriumMechanismStatus, PARAMETERS};
use amr_project::simulation::population::{
    bacterium_mechanism_host_is_eligible, mechanism_is_hgt_transferable, ResistanceMechanism,
    BACTERIA_LIST, DRUG_SHORT_NAMES,
};
use std::collections::{BTreeMap, BTreeSet, HashSet};

const CONFIG_RS: &str = include_str!("../src/config.rs");
const MAIN_RS: &str = include_str!("../src/main.rs");
const RULES_RS: &str = include_str!("../src/rules/mod.rs");
const SIMULATION_RS: &str = include_str!("../src/simulation/simulation.rs");

fn skip_ascii_whitespace(source: &str, mut offset: usize) -> usize {
    let bytes = source.as_bytes();
    while offset < bytes.len() && bytes[offset].is_ascii_whitespace() {
        offset += 1;
    }
    offset
}

fn read_quoted_value(source: &str, offset: usize) -> Option<(String, usize)> {
    let bytes = source.as_bytes();
    if bytes.get(offset).copied() != Some(b'"') {
        return None;
    }

    let start = offset + 1;
    let relative_end = source[start..].find('"')?;
    let end = start + relative_end;
    Some((source[start..end].to_string(), end + 1))
}

fn collect_string_after(source: &str, marker: &str) -> HashSet<String> {
    let mut values = HashSet::new();
    let mut offset = 0;

    while let Some(relative_start) = source[offset..].find(marker) {
        let start = offset + relative_start + marker.len();
        let Some(relative_end) = source[start..].find('"') else {
            break;
        };
        values.insert(source[start..start + relative_end].to_string());
        offset = start + relative_end + 1;
    }

    values
}

fn collect_map_insert_literal_keys(source: &str) -> HashSet<String> {
    let mut values = HashSet::new();
    let mut offset = 0;
    let marker = "map.insert(";

    while let Some(relative_start) = source[offset..].find(marker) {
        let mut cursor = offset + relative_start + marker.len();
        cursor = skip_ascii_whitespace(source, cursor);

        if let Some((value, end)) = read_quoted_value(source, cursor) {
            values.insert(value);
            offset = end;
        } else {
            offset = cursor.saturating_add(1).min(source.len());
        }
    }

    values
}

fn parameter_initializer(source: &str) -> &str {
    let start_marker = "pub static ref PARAMETERS: HashMap<String, f64> = {";
    let end_marker = "// --- String Parameters";
    let (_, after_start) = source
        .split_once(start_marker)
        .expect("PARAMETERS initializer start marker should exist");
    let (initializer, _) = after_start
        .split_once(end_marker)
        .expect("PARAMETERS initializer end marker should exist");
    initializer
}

fn find_insert_closing_parenthesis(source: &str, offset: usize) -> Option<usize> {
    let mut nested_parentheses = 0;

    for (relative_offset, character) in source[offset..].char_indices() {
        match character {
            '(' => nested_parentheses += 1,
            ')' if nested_parentheses == 0 => return Some(offset + relative_offset),
            ')' => nested_parentheses -= 1,
            _ => {}
        }
    }

    None
}

fn collect_literal_parameter_insertions(source: &str) -> BTreeMap<String, Vec<f64>> {
    let source = parameter_initializer(source);
    let mut values = BTreeMap::<String, Vec<f64>>::new();
    let mut offset = 0;
    let marker = "map.insert(";

    while let Some(relative_start) = source[offset..].find(marker) {
        let mut cursor = offset + relative_start + marker.len();
        cursor = skip_ascii_whitespace(source, cursor);

        let Some((key, end)) = read_quoted_value(source, cursor) else {
            offset = cursor.saturating_add(1).min(source.len());
            continue;
        };
        cursor = skip_ascii_whitespace(source, end);

        let to_string = ".to_string()";
        if !source[cursor..].starts_with(to_string) {
            offset = cursor.saturating_add(1).min(source.len());
            continue;
        }
        cursor = skip_ascii_whitespace(source, cursor + to_string.len());

        if source.as_bytes().get(cursor).copied() != Some(b',') {
            offset = cursor.saturating_add(1).min(source.len());
            continue;
        }
        cursor = skip_ascii_whitespace(source, cursor + 1);

        let closing_parenthesis = find_insert_closing_parenthesis(source, cursor)
            .unwrap_or_else(|| panic!("map.insert for {key} should have a closing parenthesis"));
        let expression = source[cursor..closing_parenthesis]
            .trim()
            .trim_end_matches(',')
            .trim();
        let value = expression
            .replace('_', "")
            .parse::<f64>()
            .unwrap_or_else(|_| {
                panic!(
                    "literal parameter {key} should have a numeric literal value, got {expression}"
                )
            });

        values.entry(key).or_default().push(value);
        offset = closing_parenthesis + 1;
    }

    values
}

fn duplicate_literal_parameter_values(source: &str) -> BTreeMap<String, Vec<f64>> {
    collect_literal_parameter_insertions(source)
        .into_iter()
        .filter(|(_, values)| values.len() > 1)
        .collect()
}

fn collect_get_or_default_literal_keys(source: &str) -> HashSet<String> {
    let mut values = HashSet::new();
    let mut offset = 0;
    let marker = "get_or_default(";

    while let Some(relative_start) = source[offset..].find(marker) {
        let mut cursor = offset + relative_start + marker.len();
        cursor = skip_ascii_whitespace(source, cursor);

        if !source[cursor..].starts_with("map") {
            offset = cursor.saturating_add(1).min(source.len());
            continue;
        }
        cursor += "map".len();
        cursor = skip_ascii_whitespace(source, cursor);

        if !source[cursor..].starts_with(',') {
            offset = cursor.saturating_add(1).min(source.len());
            continue;
        }
        cursor += 1;
        cursor = skip_ascii_whitespace(source, cursor);

        if let Some((value, end)) = read_quoted_value(source, cursor) {
            values.insert(value);
            offset = end;
        } else {
            offset = cursor.saturating_add(1).min(source.len());
        }
    }

    values
}

#[test]
fn literal_get_global_param_keys_exist_in_parameters() {
    let parameter_keys = collect_map_insert_literal_keys(CONFIG_RS);
    let lookup_keys = [CONFIG_RS, MAIN_RS, RULES_RS, SIMULATION_RS]
        .into_iter()
        .flat_map(|source| collect_string_after(source, "get_global_param(\""))
        .collect::<HashSet<_>>();

    let mut missing_keys = lookup_keys
        .difference(&parameter_keys)
        .cloned()
        .collect::<Vec<_>>();
    missing_keys.sort();

    assert!(
        missing_keys.is_empty(),
        "literal get_global_param keys missing from PARAMETERS: {missing_keys:?}"
    );
}

#[test]
fn literal_get_or_default_keys_exist_in_parameters() {
    let parameter_keys = collect_map_insert_literal_keys(CONFIG_RS);
    let lookup_keys = collect_get_or_default_literal_keys(CONFIG_RS);

    let mut missing_keys = lookup_keys
        .difference(&parameter_keys)
        .cloned()
        .collect::<Vec<_>>();
    missing_keys.sort();

    assert!(
        missing_keys.is_empty(),
        "literal get_or_default keys missing from PARAMETERS: {missing_keys:?}"
    );
}

#[test]
fn parameters_initializer_has_no_duplicate_literal_keys() {
    let duplicates = duplicate_literal_parameter_values(CONFIG_RS);

    assert!(
        duplicates.is_empty(),
        "literal duplicate keys are forbidden in the PARAMETERS initializer: {duplicates:#?}"
    );
}

#[test]
fn retired_resistance_parameter_names_are_absent() {
    for key in [
        "any_r_emergence_level_on_first_emergence",
        "multi_drug_penalty_for_single_drug_resistance",
        "multi_drug_penalty_for_partial_cross_resistance",
    ] {
        assert!(
            !PARAMETERS.contains_key(key),
            "retired resistance parameter returned to PARAMETERS: {key}"
        );
    }
}

#[test]
fn bacterium_mechanism_emergence_grid_is_complete_and_exact() {
    let expected = BACTERIA_LIST
        .iter()
        .flat_map(|bacteria| {
            ResistanceMechanism::all().iter().map(move |mechanism| {
                format!(
                    "bacteria_{bacteria}_mechanism_{}_emergence_rate",
                    mechanism.as_str()
                )
            })
        })
        .collect::<BTreeSet<_>>();
    let actual = PARAMETERS
        .keys()
        .filter(|key| {
            key.starts_with("bacteria_")
                && key.contains("_mechanism_")
                && key.ends_with("_emergence_rate")
        })
        .cloned()
        .collect::<BTreeSet<_>>();

    let missing = expected.difference(&actual).cloned().collect::<Vec<_>>();
    let unexpected = actual.difference(&expected).cloned().collect::<Vec<_>>();

    assert_eq!(
        actual.len(),
        BACTERIA_LIST.len() * ResistanceMechanism::all().len(),
        "expected one emergence-rate key for every bacterium-mechanism pair"
    );
    assert!(
        missing.is_empty() && unexpected.is_empty(),
        "emergence-rate grid mismatch; missing={missing:#?}, unexpected={unexpected:#?}"
    );
}

#[test]
fn bacterium_mechanism_status_matrix_has_explicit_route_semantics() {
    let store = parameter_store();
    let mut statuses_seen = BTreeSet::new();

    for bacteria_idx in 0..BACTERIA_LIST.len() {
        for (mechanism_idx, &mechanism) in ResistanceMechanism::all().iter().enumerate() {
            let host_is_eligible = bacterium_mechanism_host_is_eligible(bacteria_idx, mechanism);
            let emergence_rate = store
                .bacteria_mechanism_emergence
                .rate(bacteria_idx, mechanism_idx);
            let expected = if !host_is_eligible {
                BacteriumMechanismStatus::ExcludedHost
            } else if emergence_rate > 0.0 {
                BacteriumMechanismStatus::DeNovo
            } else if mechanism_is_hgt_transferable(mechanism) {
                BacteriumMechanismStatus::HgtOnly
            } else {
                BacteriumMechanismStatus::EligibleNoDeNovo
            };
            let actual = store
                .bacteria_mechanism_status
                .status(bacteria_idx, mechanism_idx);
            assert_eq!(
                actual,
                expected,
                "status mismatch for {} / {}",
                BACTERIA_LIST[bacteria_idx],
                mechanism.as_str()
            );

            let bit = 1u64 << mechanism_idx;
            assert_eq!(
                store
                    .bacteria_mechanism_status
                    .host_eligible_mask(bacteria_idx)
                    & bit
                    != 0,
                actual.host_is_eligible()
            );
            assert_eq!(
                store.bacteria_mechanism_status.de_novo_mask(bacteria_idx) & bit != 0,
                actual.allows_de_novo()
            );
            assert_eq!(
                store
                    .bacteria_mechanism_status
                    .hgt_recipient_mask(bacteria_idx)
                    & bit
                    != 0,
                actual.host_is_eligible() && mechanism_is_hgt_transferable(mechanism)
            );
            statuses_seen.insert(format!("{actual:?}"));
        }
    }

    assert_eq!(
        statuses_seen,
        ["DeNovo", "EligibleNoDeNovo", "ExcludedHost", "HgtOnly"]
            .into_iter()
            .map(str::to_string)
            .collect(),
        "all four route statuses should be represented in the current matrix"
    );
}

#[test]
fn reviewed_above_unit_raw_potencies_are_clamped_in_the_typed_matrix() {
    let store = parameter_store();
    let cases = [
        (
            "drug_trim_sulf_for_bacteria_stenotrophomonas_maltophilia_potency_when_no_r",
            "stenotrophomonas_maltophilia",
            "trim_sulf",
        ),
        (
            "drug_vancomycin_for_bacteria_staphylococcus_epidermidis_potency_when_no_r",
            "staphylococcus_epidermidis",
            "vancomycin",
        ),
    ];

    for (key, bacteria, drug) in cases {
        assert_eq!(
            PARAMETERS.get(key).map(|value| value.to_bits()),
            Some(1.05_f64.to_bits()),
            "{key} should preserve the reviewed raw calibration input"
        );

        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|candidate| *candidate == bacteria)
            .unwrap_or_else(|| panic!("unknown bacteria {bacteria}"));
        let drug_idx = DRUG_SHORT_NAMES
            .iter()
            .position(|candidate| *candidate == drug)
            .unwrap_or_else(|| panic!("unknown drug {drug}"));

        assert_eq!(
            store
                .drug_bacteria
                .potency(bacteria_idx, drug_idx)
                .to_bits(),
            1.0_f64.to_bits(),
            "{key} should be clamped before simulation use"
        );
    }
}
