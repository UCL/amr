# Configuration Parameters

> **Source File**: `src/config.rs` (~13,000+ lines)

This document provides a reference guide to all configurable parameters in the simulation, organized by category.

---

## Table of Contents
1. [Parameter System Overview](#parameter-system-overview)
2. [Section A: Bacteria Parameters](#section-a-bacteria-parameters)
3. [Section B: Infection Parameters](#section-b-infection-parameters)
4. [Section C: Drug Parameters](#section-c-drug-parameters)
5. [Section D: Resistance Parameters](#section-d-resistance-parameters)
6. [Section E: Potency Matrix](#section-e-potency-matrix)
7. [Section F: Mortality Parameters](#section-f-mortality-parameters)
8. [Section G: Regional Parameters](#section-g-regional-parameters)
9. [Section H: Hospital Parameters](#section-h-hospital-parameters)
10. [Section I: Population Parameters](#section-i-population-parameters)
11. [Section J: Drug Availability](#section-j-drug-availability)
12. [Section K: Simulation Control](#section-k-simulation-control)
13. [How to Add/Modify Parameters](#how-to-addmodify-parameters)

---

## Parameter System Overview

### Architecture

```rust
// Parameters stored in a HashMap
pub type Params = HashMap<String, f64>;

// Created at startup
pub fn get_params() -> Params {
    let mut map = HashMap::new();
    
    // Insert all parameters
    map.insert("param_name".to_string(), value);
    // ... thousands of parameters ...
    
    map
}
```

### Naming Conventions

```
{category}_{entity}_{property}

Examples:
bacteria_e_coli_growth_rate
drug_ampicillin_half_life_days
region_north_america_prevalence_multiplier
hospital_icu_acquisition_rate
```

### Accessing Parameters

```rust
// In rules code
let params: &Params = &config::get_params();

// Get single parameter
let growth_rate = params.get("bacteria_e_coli_growth_rate").unwrap_or(&0.0);

// Helper functions exist for common patterns
let potency = get_drug_bacteria_potency(drug_idx, bacteria_idx, params);
```

---

## Section A: Bacteria Parameters

### Basic Bacteria Properties

```rust
// Growth rates (per day when infecting)
bacteria_e_coli_growth_rate = 0.15
bacteria_k_pneumoniae_growth_rate = 0.12
bacteria_s_aureus_growth_rate = 0.18
bacteria_p_aeruginosa_growth_rate = 0.10
// ... for all 39 bacteria

// Natural clearance rates (immune-mediated, per day)
bacteria_e_coli_natural_clearance_rate = 0.05
bacteria_k_pneumoniae_natural_clearance_rate = 0.04
bacteria_s_aureus_natural_clearance_rate = 0.03
// ...
```

### Colonization Properties

```rust
// Baseline colonization acquisition rates
bacteria_e_coli_colonization_rate = 0.001
bacteria_k_pneumoniae_colonization_rate = 0.0005
bacteria_s_aureus_colonization_rate = 0.002

// Colonization equilibrium levels (0-1)
bacteria_e_coli_equilibrium_colonization_level = 0.5
bacteria_k_pneumoniae_equilibrium_colonization_level = 0.3
bacteria_s_aureus_equilibrium_colonization_level = 0.4

// Colonization clearance rates
bacteria_e_coli_colonization_clearance_rate = 0.01
bacteria_s_aureus_colonization_clearance_rate = 0.005
```

### Transmission Properties

```rust
// Contact transmission rates (per contact-day)
bacteria_e_coli_contact_transmission_rate = 0.01
bacteria_s_aureus_contact_transmission_rate = 0.02
bacteria_c_difficile_contact_transmission_rate = 0.05

// Shedding thresholds (minimum level to transmit)
bacteria_e_coli_shedding_threshold = 0.05
```

---

## Section B: Infection Parameters

### Symptom Development

```rust
// Symptom threshold (level needed for symptoms)
bacteria_e_coli_symptom_threshold = 0.1
bacteria_s_aureus_symptom_threshold = 0.15
bacteria_p_aeruginosa_symptom_threshold = 0.08

// Symptom probability (per day above threshold)
bacteria_e_coli_symptom_probability = 0.3
bacteria_s_aureus_symptom_probability = 0.4
```

### Sepsis Parameters

```rust
// Base sepsis probability (per day when symptomatic)
bacteria_e_coli_sepsis_probability = 0.02
bacteria_s_aureus_sepsis_probability = 0.03
bacteria_k_pneumoniae_sepsis_probability = 0.025
bacteria_p_aeruginosa_sepsis_probability = 0.04

// Sepsis risk multipliers
sepsis_age_65_plus_multiplier = 2.0
sepsis_age_80_plus_multiplier = 1.5  // Compounds with 65+
sepsis_untreated_infection_multiplier = 2.0
sepsis_high_level_infection_multiplier = 1.5
```

### Syndrome Distribution

```rust
// Distribution of syndromes by bacteria (sums to 1.0)
bacteria_e_coli_syndrome_uti_probability = 0.6
bacteria_e_coli_syndrome_bacteremia_probability = 0.2
bacteria_e_coli_syndrome_gi_probability = 0.1
bacteria_e_coli_syndrome_other_probability = 0.1

bacteria_s_aureus_syndrome_skin_probability = 0.5
bacteria_s_aureus_syndrome_bacteremia_probability = 0.2
bacteria_s_aureus_syndrome_pneumonia_probability = 0.15
bacteria_s_aureus_syndrome_endocarditis_probability = 0.05
bacteria_s_aureus_syndrome_other_probability = 0.1
```

---

## Section C: Drug Parameters

### Pharmacokinetics

```rust
// Half-lives (days)
drug_ampicillin_half_life_days = 0.5
drug_ceftriaxone_half_life_days = 1.0
drug_ciprofloxacin_half_life_days = 4.0
drug_azithromycin_half_life_days = 7.0
drug_vancomycin_half_life_days = 5.0
// ... for all 52 drugs

// Standard dose level
drug_standard_dose_level = 10.0

// Dose multipliers (for dose adjustments)
drug_ampicillin_dose_multiplier = 1.0
drug_meropenem_dose_multiplier = 1.0
```

### Treatment Duration

```rust
// Standard treatment durations (days) by syndrome
drug_ampicillin_duration_uti = 5
drug_ampicillin_duration_pneumonia = 7
drug_ampicillin_duration_skin = 10
drug_ampicillin_duration_bacteremia = 14

drug_ciprofloxacin_duration_uti = 3
drug_ciprofloxacin_duration_pneumonia = 7
// ...
```

### Toxicity

```rust
// Toxicity accumulation rates (per day per unit drug level)
drug_gentamicin_toxicity_accumulation_rate = 0.02
drug_vancomycin_toxicity_accumulation_rate = 0.015
drug_colistin_toxicity_accumulation_rate = 0.03
drug_ampicillin_toxicity_accumulation_rate = 0.005  // Lower toxicity

// Toxicity decay rates (per day)
drug_gentamicin_toxicity_decay_rate = 0.1
drug_vancomycin_toxicity_decay_rate = 0.1
drug_ampicillin_toxicity_decay_rate = 0.2

// Toxicity mortality parameters
toxicity_mortality_midpoint = 50.0
toxicity_mortality_slope = 0.1
toxicity_mortality_scale = 0.01
```

### Adherence

```rust
// Adherence probabilities by treatment phase
drug_ampicillin_adherence_day_1_to_3 = 0.95
drug_ampicillin_adherence_day_4_to_7 = 0.85
drug_ampicillin_adherence_day_8_plus = 0.75

drug_azithromycin_adherence_day_1_to_3 = 0.98  // Once daily = better
drug_azithromycin_adherence_day_4_to_7 = 0.95
```

---

## Section D: Resistance Parameters

### Resistance Dynamics

```rust
// De novo emergence rates (per day during treatment)
bacteria_e_coli_drug_ciprofloxacin_de_novo_rate = 0.0001
bacteria_p_aeruginosa_drug_meropenem_de_novo_rate = 0.0005
bacteria_s_aureus_drug_vancomycin_de_novo_rate = 0.00001  // Very rare

// Resistance reversion rates (per day without drug)
bacteria_e_coli_drug_ciprofloxacin_reversion_rate = 0.001
bacteria_k_pneumoniae_drug_meropenem_reversion_rate = 0.0005
```

### HGT Parameters

```rust
// Base HGT rates by taxonomic distance
hgt_base_rate_same_genus = 0.001
hgt_base_rate_same_family = 0.0001
hgt_base_rate_different_family = 0.00001
hgt_base_rate_gram_switch = 0.000001
```

### Resistance Mechanisms

```rust
// Mechanism properties
mechanism_esbl_is_mobile = 1.0  // 1.0 = true
mechanism_esbl_transfer_probability = 0.1
mechanism_esbl_fitness_cost = 0.05

mechanism_carbapenemase_is_mobile = 1.0
mechanism_carbapenemase_transfer_probability = 0.05
mechanism_carbapenemase_fitness_cost = 0.1

mechanism_meca_is_mobile = 0.0  // Chromosomal
mechanism_meca_transfer_probability = 0.0
```

### Resistance Floors

```rust
// For rare bacteria with intrinsic resistance
resistance_floor_s_maltophilia_trim_sulf = 0.05
resistance_floor_s_maltophilia_ciprofloxacin = 0.1
resistance_floor_e_faecium_vancomycin = 0.03
// ... etc

// Enable/disable
use_resistance_floors = 1.0  // 1.0 = enabled
```

### MajorityR Cache

```rust
// Population sampling for majority_r
majority_r_sample_size = 100
majority_r_refresh_interval_days = 7
```

---

## Section E: Potency Matrix

### Drug-Bacteria Potency

This is a large matrix (~2000 entries) defining drug effectiveness:

```rust
// Potency when no resistance (0.0-1.0+)
drug_ampicillin_for_bacteria_e_coli_potency_when_no_r = 0.7
drug_ampicillin_for_bacteria_s_aureus_potency_when_no_r = 0.8
drug_ampicillin_for_bacteria_p_aeruginosa_potency_when_no_r = 0.0  // No activity

drug_meropenem_for_bacteria_e_coli_potency_when_no_r = 1.0
drug_meropenem_for_bacteria_p_aeruginosa_potency_when_no_r = 0.9
drug_meropenem_for_bacteria_s_maltophilia_potency_when_no_r = 0.05  // Intrinsic R

drug_vancomycin_for_bacteria_s_aureus_potency_when_no_r = 1.0
drug_vancomycin_for_bacteria_e_coli_potency_when_no_r = 0.0  // Gram-negative
```

### Initiation Multipliers

```rust
// Prescribing preference multipliers
drug_ampicillin_for_bacteria_e_coli_initiation_multiplier = 0.5  // Not first choice
drug_ciprofloxacin_for_bacteria_e_coli_initiation_multiplier = 1.0  // Common
drug_meropenem_for_bacteria_e_coli_initiation_multiplier = 0.3  // Reserved

drug_vancomycin_for_bacteria_s_aureus_initiation_multiplier = 0.5  // MRSA only
drug_cefazolin_for_bacteria_s_aureus_initiation_multiplier = 1.0  // First line
```

---

## Section F: Mortality Parameters

### Infection Mortality

```rust
// Base infection mortality rates (per day, severe infection)
bacteria_e_coli_infection_mortality_rate = 0.01
bacteria_s_aureus_infection_mortality_rate = 0.015
bacteria_p_aeruginosa_infection_mortality_rate = 0.02
bacteria_k_pneumoniae_infection_mortality_rate = 0.012

// Mortality multipliers
sepsis_mortality_multiplier = 5.0
septic_shock_mortality_multiplier = 2.0  // Additional
bacteremia_syndrome_mortality_multiplier = 1.5
meningitis_syndrome_mortality_multiplier = 2.0
endocarditis_syndrome_mortality_multiplier = 2.5
```

### Age-Related Mortality

```rust
age_mortality_multiplier_0_to_1 = 1.5    // Infants
age_mortality_multiplier_1_to_5 = 1.0
age_mortality_multiplier_5_to_18 = 0.8   // Children resilient
age_mortality_multiplier_18_to_65 = 1.0
age_mortality_multiplier_65_to_80 = 2.0
age_mortality_multiplier_80_plus = 3.0
```

### Treatment Effect on Mortality

```rust
effective_treatment_mortality_reduction = 0.8  // 80% reduction
ineffective_treatment_mortality_reduction = 0.3  // 30% reduction
no_treatment_mortality_multiplier = 2.0
```

---

## Section G: Regional Parameters

### Regional Prevalence

```rust
// Bacteria prevalence multipliers by region
north_america_bacteria_e_coli_prevalence = 1.0
north_america_bacteria_k_pneumoniae_prevalence = 1.0
north_america_bacteria_p_aeruginosa_prevalence = 0.8

europe_bacteria_e_coli_prevalence = 1.1
europe_bacteria_k_pneumoniae_prevalence = 1.2

south_asia_bacteria_e_coli_prevalence = 1.3
south_asia_bacteria_k_pneumoniae_prevalence = 1.5

africa_bacteria_s_pneumoniae_prevalence = 1.4
```

### Regional Resistance

```rust
// Baseline resistance levels by region
north_america_bacteria_e_coli_drug_ciprofloxacin_baseline_r = 0.2
europe_bacteria_e_coli_drug_ciprofloxacin_baseline_r = 0.25
south_asia_bacteria_e_coli_drug_ciprofloxacin_baseline_r = 0.5

north_america_bacteria_k_pneumoniae_drug_meropenem_baseline_r = 0.02
south_asia_bacteria_k_pneumoniae_drug_meropenem_baseline_r = 0.15
```

---

## Section H: Hospital Parameters

### Hospital Acquisition

```rust
// Hospital acquisition multipliers
hospital_bacteria_e_coli_acquisition_multiplier = 2.0
hospital_bacteria_k_pneumoniae_acquisition_multiplier = 3.0
hospital_bacteria_p_aeruginosa_acquisition_multiplier = 5.0
hospital_bacteria_a_baumannii_acquisition_multiplier = 10.0
hospital_bacteria_c_difficile_acquisition_multiplier = 8.0

// ICU multipliers (additional)
icu_bacteria_p_aeruginosa_acquisition_multiplier = 2.0
icu_bacteria_a_baumannii_acquisition_multiplier = 3.0
```

### Hospital Resistance

```rust
// Hospital strains are more resistant
hospital_acquired_resistance_multiplier = 1.5
hospital_carbapenem_resistant_probability = 0.2
hospital_esbl_probability = 0.3
```

### Length of Stay

```rust
// LOS effect on acquisition
los_acquisition_multiplier_per_day = 0.05  // 5% increase per day
los_max_multiplier = 3.0
```

---

## Section I: Population Parameters

### Demographics

```rust
// Age distribution (sums to 1.0)
age_group_0_to_1_proportion = 0.012
age_group_1_to_5_proportion = 0.05
age_group_5_to_18_proportion = 0.17
age_group_18_to_65_proportion = 0.62
age_group_65_to_80_proportion = 0.12
age_group_80_plus_proportion = 0.028

// Penicillin allergy
penicillin_allergy_prevalence = 0.10  // 10% of population
```

### Population Size

```rust
population_size = 100000
```

---

## Section J: Drug Availability

### Historical Introduction

```rust
// Drug introduction years (as time steps from simulation start)
// DRUG_INTRODUCTION_DATES HashMap in config.rs

// Example entries:
// "sulfanilamide" => 2555 (1937)
// "penicilling" => 3555 (1942)
// "ampicillin" => 11315 (1961)
// "ciprofloxacin" => 20805 (1987)
// "linezolid" => 25550 (2000)
```

### Regional Availability

```rust
// Drug availability by region (0.0-1.0)
north_america_drug_linezolid_availability = 1.0
north_america_drug_dalbavancin_availability = 0.8

africa_drug_linezolid_availability = 0.1
africa_drug_dalbavancin_availability = 0.05

south_asia_drug_meropenem_availability = 0.5
south_asia_drug_colistin_availability = 0.3
```

---

## Section K: Simulation Control

### Time Parameters

```rust
simulation_start_year = 1930
simulation_end_year = 2030
time_steps_per_day = 1
days_per_year = 365
```

### Random Seed

```rust
random_seed = 12345  // Set for reproducibility, or 0 for random
```

### Output Control

```rust
output_summary_interval_days = 365  // Yearly summaries
output_individual_traces = 0  // 0 = disabled, 1 = enabled
output_detailed_resistance = 1
```

### Calibration Toggles

```rust
// Enable/disable specific mechanisms
enable_hgt = 1.0
enable_de_novo_emergence = 1.0
enable_microbiome_seeding = 1.0
enable_contact_transmission = 1.0
enable_hospital_acquisition = 1.0
```

---

## How to Add/Modify Parameters

### Adding a New Parameter

1. **Add to config.rs**:
```rust
// In get_params() function
map.insert("new_category_new_parameter".to_string(), default_value);
```

2. **Use in rules code**:
```rust
let value = params.get("new_category_new_parameter").unwrap_or(&default);
```

### Creating Helper Functions

For commonly accessed patterns, create helpers:

```rust
pub fn get_bacteria_growth_rate(bacteria_name: &str, params: &Params) -> f64 {
    let key = format!("bacteria_{}_growth_rate", bacteria_name);
    *params.get(&key).unwrap_or(&0.1)  // Default 0.1
}

pub fn get_drug_bacteria_potency(drug_idx: usize, bacteria_idx: usize, params: &Params) -> f64 {
    let drug = DRUG_SHORT_NAMES[drug_idx];
    let bacteria = BACTERIA_LIST[bacteria_idx];
    let key = format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria);
    *params.get(&key).unwrap_or(&0.0)
}
```

### Validation

Add parameter validation on startup:

```rust
pub fn validate_params(params: &Params) -> Result<(), String> {
    // Check required parameters exist
    for bacteria in BACTERIA_LIST {
        let key = format!("bacteria_{}_growth_rate", bacteria);
        if !params.contains_key(&key) {
            return Err(format!("Missing parameter: {}", key));
        }
    }
    
    // Check value ranges
    if let Some(&v) = params.get("population_size") {
        if v < 1.0 {
            return Err("population_size must be >= 1".to_string());
        }
    }
    
    Ok(())
}
```

### Parameter Sensitivity Analysis

To test parameter sensitivity:

```rust
// In main.rs or separate analysis script
fn sensitivity_analysis() {
    let base_params = get_params();
    let test_param = "bacteria_e_coli_growth_rate";
    let base_value = base_params.get(test_param).unwrap();
    
    for multiplier in [0.5, 0.75, 1.0, 1.25, 1.5] {
        let mut params = base_params.clone();
        params.insert(test_param.to_string(), base_value * multiplier);
        
        let result = run_simulation(&params);
        println!("{} = {}: outcome = {}", test_param, base_value * multiplier, result);
    }
}
```

---

## Quick Reference: Parameter Counts

| Section | Approximate Count | Description |
|---------|-------------------|-------------|
| A: Bacteria | ~200 | Per-bacteria properties |
| B: Infection | ~150 | Symptoms, sepsis, syndromes |
| C: Drug | ~300 | PK, duration, toxicity |
| D: Resistance | ~500 | Emergence, reversion, HGT |
| E: Potency | ~2000 | Drug-bacteria matrix |
| F: Mortality | ~100 | Death rates |
| G: Regional | ~500 | Geographic variation |
| H: Hospital | ~100 | Nosocomial factors |
| I: Population | ~20 | Demographics |
| J: Availability | ~200 | Drug access |
| K: Control | ~20 | Simulation settings |
| **Total** | **~4000+** | |
