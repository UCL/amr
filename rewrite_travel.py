import sys

with open('src/rules/mod.rs', 'r', encoding='utf-8') as f:
    content = f.read()

old_travel_loop = """            let mut new_region: Region;
            loop {
                // Select destination based on economic development level (main determinant of travel patterns)
                // Higher-income regions have more global travel; lower-income regions travel more regionally
                let destinations = match individual.region_living {
                    Region::NorthAmerica | Region::Europe | Region::Oceania => {
                        // High-income regions: global travel with preference for other developed regions
                        vec![
                            (Region::Europe, 0.35),       // Strong developed-to-developed flow
                            (Region::Asia, 0.25),         // Major business/tourism destination
                            (Region::NorthAmerica, 0.15), // Cross-Atlantic travel
                            (Region::Oceania, 0.10),      // Tourism/business
                            (Region::SouthAmerica, 0.10), // Tourism/business
                            (Region::Africa, 0.05),       // Lower but still significant
                        ]
                    }
                    Region::Asia => {
                        // Mixed income: regional preference with some global reach
                        vec![
                            (Region::Asia, 0.40),         // Strong regional travel
                            (Region::Europe, 0.20),       // Business/education
                            (Region::NorthAmerica, 0.15), // Business/education
                            (Region::Oceania, 0.10),      // Regional proximity
                            (Region::Africa, 0.08),       // Growing connections
                            (Region::SouthAmerica, 0.07), // Limited
                        ]
                    }
                    Region::SouthAmerica => {
                        // Middle income: regional focus with some international travel
                        vec![
                            (Region::SouthAmerica, 0.40), // Strong regional travel
                            (Region::NorthAmerica, 0.25), // Geographic proximity
                            (Region::Europe, 0.15),       // Historical ties
                            (Region::Asia, 0.10),         // Growing connections
                            (Region::Africa, 0.05),       // Limited
                            (Region::Oceania, 0.05),      // Limited
                        ]
                    }
                    Region::Africa => {
                        // Lower income: primarily regional travel
                        vec![
                            (Region::Africa, 0.50),       // Strong regional travel
                            (Region::Europe, 0.20),       // Historical/economic ties
                            (Region::Asia, 0.15),         // Growing connections
                            (Region::NorthAmerica, 0.08), // Limited
                            (Region::SouthAmerica, 0.04), // Very limited
                            (Region::Oceania, 0.03),      // Very limited
                        ]
                    }
                    Region::Home => {
                        // Should not reach here, but default to global uniform if it does
                        vec![
                            (Region::Asia, 0.167),
                            (Region::Africa, 0.167),
                            (Region::Europe, 0.166),
                            (Region::NorthAmerica, 0.167),
                            (Region::SouthAmerica, 0.166),
                            (Region::Oceania, 0.167),
                        ]
                    }
                };

                // Sample from the economic-based destination distribution
                let rand_val = rng.gen::<f64>();
                let mut cumulative_prob = 0.0;
                new_region = Region::Asia; // Default fallback

                for (region, prob) in destinations {
                    cumulative_prob += prob;
                    if rand_val < cumulative_prob {
                        new_region = region;
                        break;
                    }
                }

                // Ensure the individual doesn't 'travel' to their own living region
                if new_region != individual.region_living {
                    break; // Found a suitable new region to visit
                }
            }"""

new_travel_loop = """            // We pre-define standard travel matrix probabilities.
            // Notice: it no longer uses dynamic vectors!
            let (raw_destinations, len) = match individual.region_living {
                Region::NorthAmerica | Region::Europe | Region::Oceania => (
                    [
                        (Region::Europe, 0.35),
                        (Region::Asia, 0.25),
                        (Region::NorthAmerica, 0.15),
                        (Region::Oceania, 0.10),
                        (Region::SouthAmerica, 0.10),
                        (Region::Africa, 0.05),
                    ], 6)
                ,
                Region::Asia => (
                    [
                        (Region::Asia, 0.40),
                        (Region::Europe, 0.20),
                        (Region::NorthAmerica, 0.15),
                        (Region::Oceania, 0.10),
                        (Region::Africa, 0.08),
                        (Region::SouthAmerica, 0.07),
                    ], 6)
                ,
                Region::SouthAmerica => (
                    [
                        (Region::SouthAmerica, 0.40),
                        (Region::NorthAmerica, 0.25),
                        (Region::Europe, 0.15),
                        (Region::Asia, 0.10),
                        (Region::Africa, 0.05),
                        (Region::Oceania, 0.05),
                    ], 6)
                ,
                Region::Africa => (
                    [
                        (Region::Africa, 0.50),
                        (Region::Europe, 0.20),
                        (Region::Asia, 0.15),
                        (Region::NorthAmerica, 0.08),
                        (Region::SouthAmerica, 0.04),
                        (Region::Oceania, 0.03),
                    ], 6)
                ,
                Region::Home | _ => (
                    [
                        (Region::Asia, 0.167),
                        (Region::Africa, 0.167),
                        (Region::Europe, 0.166),
                        (Region::NorthAmerica, 0.167),
                        (Region::SouthAmerica, 0.166),
                        (Region::Oceania, 0.167),
                    ], 6)
                ,
            };

            let mut valid_destinations = [(Region::Home, 0.0); 6];
            let mut dest_count = 0;
            let mut total_weight = 0.0;
            for i in 0..len {
                let dest = raw_destinations[i].0;
                let weight = raw_destinations[i].1;
                if dest != individual.region_living {
                    valid_destinations[dest_count] = (dest, weight);
                    total_weight += weight;
                    dest_count += 1;
                }
            }

            let mut rand_val = rng.gen::<f64>() * total_weight;
            let mut new_region = valid_destinations[dest_count - 1].0; // Default to last
            for i in 0..dest_count {
                if rand_val < valid_destinations[i].1 {
                    new_region = valid_destinations[i].0;
                    break;
                }
                rand_val -= valid_destinations[i].1;
            }
"""

if old_travel_loop in content:
    content = content.replace(old_travel_loop, new_travel_loop)
else:
    print("WARNING: EXACT MATCH NOT FOUND!")

# The earlier unused variable bacteria_name (warning from rustc) inside apply_rules 
# since `bacteria_param_name` was removed, we should remove `let bacteria_name = ...`
content = content.replace("let bacteria_name = BACTERIA_LIST[b_idx];\n        let bacteria_specific_available", "let bacteria_specific_available")

print(f"Total Bytes changed: {len(old_travel_loop)} down to {len(new_travel_loop)}")
with open('src/rules/mod.rs', 'w', encoding='utf-8') as f:
    f.write(content)