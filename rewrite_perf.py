import sys

with open('src/rules/mod.rs', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix 1: tb_guaranteed_rifampicin_resistance field in Cache
content = content.replace("pub tb_guaranteed_rifampicin_resistance: bool,", "pub tb_guaranteed_rifampicin_resistance: f64,")
content = content.replace('crate::config::get_global_param("mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance").unwrap_or(1.0) > 0.5,', 'crate::config::get_global_param("mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance").unwrap_or(0.9),')

# Fix 2: Add resistance_test_result_delay_days to cache
content = content.replace("pub test_delay_days: i32,", "pub test_delay_days: i32,\n    pub resistance_test_result_delay_days: i32,")
content = content.replace('test_delay_days: crate::config::get_global_param("test_delay_days").unwrap_or(3.0) as i32,', 'test_delay_days: crate::config::get_global_param("test_delay_days").unwrap_or(3.0) as i32,\n            resistance_test_result_delay_days: crate::config::get_global_param("resistance_test_result_delay_days").unwrap_or(2.0) as i32,')

# Fix 3: replace resistance_test_result_delay_days inside apply_rules
old_delay = """let resistance_test_result_delay_days =
        get_global_param("resistance_test_result_delay_days").unwrap_or(2.0) as i32;"""
new_delay = "let resistance_test_result_delay_days = param_cache.resistance_test_result_delay_days;"
content = content.replace(old_delay, new_delay)


# Fix 4: TB guaranteed inside apply_rules
old_tb = """let guaranteed_rifampicin_resistance = if is_tb && simulation_year >= 1966.0 {
                        // Only apply guaranteed rifampicin resistance after rifampicin is available
                        get_global_param(
                            "mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance",
                        )
                        .unwrap_or(0.90)
                    } else {
                        0.0
                    };"""
new_tb = """let guaranteed_rifampicin_resistance = if is_tb && simulation_year >= 1966.0 {
                        param_cache.tb_guaranteed_rifampicin_resistance
                    } else {
                        0.0
                    };"""
content = content.replace(old_tb, new_tb)


# Fix 5: rifampicin_idx caching
# We just replace drug_indices.get("rifampicin") with DRUG_SHORT_NAMES.iter().position(|&x| x == "rifampicin") ?
# Actually `drug_indices` in apply_rules is a HashMap passed in OR we can just use drug_indices.get...
# Wait, `drug_indices` is just a standard argument to `apply_rules`: `drug_indices: &HashMap<String, usize>`
# We can just look at `drug_indices.get("rifampicin")`. It's a quick hashmap lookup, not parsing a file.
# But wait, looking up "rifampicin" string constant in a HashMap is still a string hash. We can just do `DRUG_SHORT_NAMES.iter().position(|&n| n == "rifampicin")` earlier and pass it on, OR better: the simulation is already using indices everywhere.
old_rif = """if let Some(&rifampicin_idx) = drug_indices.get("rifampicin") {"""
new_rif = """if let Some(rifampicin_idx) = DRUG_SHORT_NAMES.iter().position(|&n| n == "rifampicin") {"""
content = content.replace(old_rif, new_rif)


# Fix 6: Bacteria specific test availability
old_bac_avail_full = """let bacteria_param_name = bacteria_name.to_lowercase().replace(" ", "_");
        let bacteria_test_availability_param =
            format!("{}_test_availability_year", bacteria_param_name);
        let bacteria_specific_available = if let Some(bacteria_discovery_year) =
            get_global_param(&bacteria_test_availability_param)
        {
            let bacteria_discovery_day = ((bacteria_discovery_year - 1930.0) * 365.25) as i32;
            time_step >= bacteria_discovery_day as usize
        } else {
            bacterial_testing_available // For most bacteria, use the general bacterial testing availability
        };"""

new_bac_avail = """let bacteria_specific_available = if let Some(bacteria_discovery_day) = param_cache.bacteria_test_availability_day[b_idx] {
            time_step >= bacteria_discovery_day
        } else {
            bacterial_testing_available
        };"""
content = content.replace(old_bac_avail_full, new_bac_avail)

# Let me rewrite the file now, we will do the travel loop next.
with open('src/rules/mod.rs', 'w', encoding='utf-8') as f:
    f.write(content)
