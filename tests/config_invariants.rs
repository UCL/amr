use std::collections::HashSet;

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
