use std::collections::HashSet;

const CONFIG_RS: &str = include_str!("../src/config.rs");
const MAIN_RS: &str = include_str!("../src/main.rs");
const RULES_RS: &str = include_str!("../src/rules/mod.rs");
const SIMULATION_RS: &str = include_str!("../src/simulation/simulation.rs");

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

#[test]
fn literal_get_global_param_keys_exist_in_parameters() {
    let parameter_keys = collect_string_after(CONFIG_RS, "map.insert(\"");
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
