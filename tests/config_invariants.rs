use amr_project::config::PARAMETERS;
use std::collections::{BTreeMap, HashSet};

const CONFIG_RS: &str = include_str!("../src/config.rs");
const MAIN_RS: &str = include_str!("../src/main.rs");
const RULES_RS: &str = include_str!("../src/rules/mod.rs");
const SIMULATION_RS: &str = include_str!("../src/simulation/simulation.rs");
const APPROVED_PARAMETER_DUPLICATES: &str =
    include_str!("fixtures/approved_parameter_duplicates.tsv");

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

fn parse_approved_parameter_duplicates(source: &str) -> BTreeMap<String, Vec<f64>> {
    let mut approved = BTreeMap::new();

    for (line_idx, line) in source.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let (key, serialized_values) = line.split_once('\t').unwrap_or_else(|| {
            panic!(
                "approved duplicate inventory line {} should be tab-separated",
                line_idx + 1
            )
        });
        let values = serialized_values
            .split(',')
            .map(|value| {
                value.parse::<f64>().unwrap_or_else(|_| {
                    panic!(
                        "approved duplicate inventory line {} has invalid value {value}",
                        line_idx + 1
                    )
                })
            })
            .collect::<Vec<_>>();

        assert!(
            values.len() > 1,
            "approved duplicate inventory entry {key} must contain at least two values"
        );
        assert!(
            approved.insert(key.to_string(), values).is_none(),
            "approved duplicate inventory contains {key} more than once"
        );
    }

    approved
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
fn literal_parameter_duplicates_match_approved_inventory() {
    let actual = duplicate_literal_parameter_values(CONFIG_RS);
    let approved = parse_approved_parameter_duplicates(APPROVED_PARAMETER_DUPLICATES);

    let unexpected = actual
        .keys()
        .filter(|key| !approved.contains_key(*key))
        .collect::<Vec<_>>();
    let missing = approved
        .keys()
        .filter(|key| !actual.contains_key(*key))
        .collect::<Vec<_>>();
    let changed = actual
        .iter()
        .filter_map(|(key, values)| {
            approved
                .get(key)
                .filter(|approved_values| *approved_values != values)
                .map(|approved_values| {
                    format!("{key}: approved={approved_values:?}, actual={values:?}")
                })
        })
        .collect::<Vec<_>>();

    assert!(
        unexpected.is_empty() && missing.is_empty() && changed.is_empty(),
        "literal PARAMETERS duplicates changed without explicit review\n\
         unexpected keys: {unexpected:?}\n\
         missing keys: {missing:?}\n\
         changed value sequences: {changed:#?}"
    );
}

#[test]
fn approved_duplicate_final_values_match_effective_parameter_map() {
    let duplicates = duplicate_literal_parameter_values(CONFIG_RS);
    let mismatches = duplicates
        .iter()
        .filter_map(|(key, values)| {
            let final_source_value = values.last().expect("duplicate should have values");
            let effective_value = PARAMETERS.get(key);

            (effective_value.map(|value| value.to_bits()) != Some(final_source_value.to_bits()))
                .then(|| {
                    format!(
                        "{key}: final source value={final_source_value}, effective={effective_value:?}"
                    )
                })
        })
        .collect::<Vec<_>>();

    assert!(
        mismatches.is_empty(),
        "final duplicate insertions should be the effective runtime values: {mismatches:#?}"
    );
}
