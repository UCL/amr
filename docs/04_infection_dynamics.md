# Infection Dynamics

> **Source Files**: 
> - `src/rules/mod.rs` → acquisition, progression, clearance
> - `src/config.rs` → infection parameters (Sections A, B)
> - `src/simulation/population.rs` → Individual struct

This document explains how infections are acquired, progress, cause symptoms, and clear.

---

## Table of Contents
1. [Infection Model Overview](#infection-model-overview)
2. [Bacterial Reservoirs](#bacterial-reservoirs)
3. [Infection Acquisition](#infection-acquisition)
4. [Infection Progression](#infection-progression)
5. [Symptoms and Syndromes](#symptoms-and-syndromes)
6. [Sepsis](#sepsis)
7. [Infection Clearance](#infection-clearance)
8. [Mortality](#mortality)
9. [Reinfection and Immunity](#reinfection-and-immunity)
10. [Configuration Parameters](#configuration-parameters)

---

## Infection Model Overview

### Conceptual Model

```
┌────────────────────────────────────────────────────────────────────────┐
│                        Infection Lifecycle                             │
├────────────────────────────────────────────────────────────────────────┤
│                                                                        │
│  No Infection     Colonization      Active Infection      Resolution  │
│                                                                        │
│  ┌────────┐      ┌────────────┐    ┌──────────────────┐   ┌────────┐  │
│  │level=0 │ ──▶  │level=0.01- │ ──▶│ level > threshold│──▶│level=0 │  │
│  │        │      │   0.1      │    │ + symptoms       │   │or death│  │
│  └────────┘      └────────────┘    └──────────────────┘   └────────┘  │
│       │                │                    │                   │      │
│       │   Acquisition  │   Progression      │    Clearance      │      │
│       ▼                ▼                    ▼                   ▼      │
│  Community         Microbiome          Drug Effect         Immune     │
│  Acquisition       Seeding             Sepsis Risk         Recovery   │
│  Nosocomial        HGT                 Mortality           or Death   │
│  Contact                                                              │
│                                                                        │
└────────────────────────────────────────────────────────────────────────┘
```

### Key Variables (per bacteria)

| Variable | Type | Purpose |
|----------|------|---------|
| `level[b]` | f64 | Infection intensity (0.0 = none) |
| `date_last_infected[b]` | Option<usize> | Day infection started |
| `symptoms[b]` | bool | Symptomatic infection |
| `sepsis[b]` | bool | Severe/septic infection |
| `syndrome[b]` | Option<Syndrome> | UTI, pneumonia, etc. |

---

## Bacterial Reservoirs

### The 39 Bacteria

```rust
pub const BACTERIA_LIST: &[&str] = &[
    // Enterobacteriaceae (Gram-negative)
    "e_coli", "k_pneumoniae", "k_oxytoca", "e_cloacae", "e_aerogenes",
    "c_freundii", "c_koseri", "s_marcescens", "p_mirabilis", "p_vulgaris",
    "m_morganii", "p_stuartii",
    
    // Non-fermenting Gram-negatives
    "p_aeruginosa", "a_baumannii", "s_maltophilia",
    
    // Other Gram-negatives
    "h_influenzae", "n_meningitidis", "n_gonorrhoeae", "m_catarrhalis",
    "b_fragilis",
    
    // Staphylococci (Gram-positive)
    "s_aureus", "s_epidermidis", "s_lugdunensis", "s_saprophyticus",
    
    // Streptococci (Gram-positive)
    "s_pneumoniae", "s_pyogenes", "s_agalactiae", "s_anginosus",
    
    // Enterococci (Gram-positive)
    "e_faecalis", "e_faecium",
    
    // Other Gram-positives
    "c_difficile", "l_monocytogenes",
    
    // Atypicals
    "m_tuberculosis", "m_avium", "c_pneumoniae", "m_pneumoniae",
    "l_pneumophila",
    
    // Other
    "h_pylori",
];
```

### Bacteria by Infection Site

| Common Site | Typical Bacteria |
|-------------|------------------|
| UTI | E. coli, K. pneumoniae, P. mirabilis, E. faecalis |
| Pneumonia | S. pneumoniae, H. influenzae, K. pneumoniae, P. aeruginosa |
| Skin/soft tissue | S. aureus, S. pyogenes, P. aeruginosa |
| Bloodstream | S. aureus, E. coli, K. pneumoniae, P. aeruginosa |
| GI | C. difficile, H. pylori, B. fragilis |

---

## Infection Acquisition

### Acquisition Sources

1. **Community Acquisition**: Random based on prevalence
2. **Nosocomial (Hospital)**: Higher rates during hospitalization
3. **Contact Transmission**: From infected individuals
4. **Microbiome Seeding**: Colonization → Infection

### Community Acquisition

```rust
// Daily community acquisition check
fn community_acquisition(individual, bacteria_idx, time_step, rng) {
    // Base rate from community prevalence
    let base_rate = community_infection_rate[bacteria_idx];
    
    // Regional adjustment
    let regional_mult = regional_bacteria_prevalence[individual.region_cur_in][bacteria_idx];
    
    // Age adjustment
    let age_mult = age_infection_susceptibility(individual.age);
    
    // Season adjustment (some bacteria)
    let season_mult = seasonal_factor(time_step, bacteria_idx);
    
    // Immune adjustment
    let immune_mult = individual.immune_susceptibility_multiplier[bacteria_idx];
    
    // Combined rate
    let acquisition_rate = base_rate * regional_mult * age_mult * season_mult * immune_mult;
    
    // Stochastic acquisition
    if rng.gen::<f64>() < acquisition_rate {
        initiate_infection(individual, bacteria_idx, time_step);
    }
}
```

### Hospital-Acquired (Nosocomial)

```rust
// Hospital environment increases acquisition
if individual.hospital_status > 0 {
    // Elevated base rate
    let nosocomial_rate = hospital_acquisition_rate[bacteria_idx];
    
    // Further multiplied by ward type
    let ward_mult = match individual.hospital_ward {
        Ward::ICU => icu_acquisition_multiplier,    // ~3x
        Ward::Surgical => surgical_acquisition_multiplier,  // ~2x
        Ward::Medical => medical_acquisition_multiplier,    // ~1.5x
        Ward::General => 1.0,
    };
    
    // Length of stay increases risk
    let los_mult = 1.0 + (individual.days_in_hospital as f64 * 0.05);
    
    // P. aeruginosa, A. baumannii, K. pneumoniae especially common
    if nosocomial_bacteria.contains(&bacteria_idx) {
        rate *= 2.0;
    }
    
    if rng.gen::<f64>() < nosocomial_rate * ward_mult * los_mult {
        initiate_infection(individual, bacteria_idx, time_step);
    }
}
```

### Contact Transmission

```rust
// Person-to-person spread (household/contact network)
for contact_idx in individual.contact_network.iter() {
    let contact = &population[contact_idx];
    
    // Contact must be infected and shedding
    if contact.level[bacteria_idx] > shedding_threshold {
        let transmission_prob = contact_transmission_rate[bacteria_idx];
        
        // Adjust for contact intensity
        let intensity = contact_intensity[individual.id][contact_idx];
        
        if rng.gen::<f64>() < transmission_prob * intensity {
            initiate_infection(individual, bacteria_idx, time_step);
            break;  // Only one source per day
        }
    }
}
```

### Microbiome Seeding

```rust
// Existing colonization can become infection
if individual.presence_microbiome[bacteria_idx] {
    // Carriage level determines risk
    let carriage_level = individual.microbiome_level[bacteria_idx];
    
    // Seeding rate
    let seed_rate = microbiome_to_infection_rate[bacteria_idx];
    
    // Compromised barriers increase risk (surgery, catheters, etc.)
    let barrier_mult = if individual.has_catheter { 3.0 }
                       else if individual.recent_surgery { 2.0 }
                       else { 1.0 };
    
    if rng.gen::<f64>() < seed_rate * carriage_level * barrier_mult {
        initiate_infection(individual, bacteria_idx, time_step);
    }
}
```

### Infection Initiation

```rust
fn initiate_infection(individual, bacteria_idx, time_step) {
    // Set initial infection level
    let initial_level = infection_initial_level[bacteria_idx];  // e.g., 0.1
    individual.level[bacteria_idx] = initial_level;
    
    // Record start date
    individual.date_last_infected[bacteria_idx] = Some(time_step);
    
    // Initialize symptoms as false (may develop)
    individual.symptoms[bacteria_idx] = false;
    
    // No sepsis initially
    individual.sepsis[bacteria_idx] = false;
    
    // Transfer resistance from microbiome if present
    if individual.presence_microbiome[bacteria_idx] {
        // Infection inherits microbiome resistance
        for drug_idx in 0..N_DRUGS {
            let micro_r = individual.microbiome_r[bacteria_idx][drug_idx];
            if micro_r > individual.resistances[bacteria_idx][drug_idx].any_r {
                individual.resistances[bacteria_idx][drug_idx].any_r = micro_r;
            }
        }
    }
}
```

---

## Infection Progression

### Level Growth

```rust
// Each time step during active infection
fn progress_infection(individual, bacteria_idx, time_step) {
    if individual.level[bacteria_idx] <= 0.0 {
        return;  // No infection
    }
    
    // Bacteria-specific growth rate
    let growth_rate = bacteria_growth_rate[bacteria_idx];  // e.g., 0.1 per day
    
    // Host immune response (slows growth)
    let immune_suppression = individual.immune_response_strength[bacteria_idx];
    
    // Drug effect (if on treatment)
    let drug_effect = calculate_total_drug_effect(individual, bacteria_idx);
    
    // Net growth = growth - immune - drug
    let net_change = growth_rate * (1.0 - immune_suppression) - drug_effect;
    
    // Update level
    individual.level[bacteria_idx] += net_change;
    
    // Clamp to valid range
    individual.level[bacteria_idx] = individual.level[bacteria_idx].max(0.0).min(max_infection_level);
}
```

### Maximum Infection Level

```rust
const MAX_INFECTION_LEVEL: f64 = 10.0;

// Level interpretation
// 0.0 - 0.1: Subclinical
// 0.1 - 1.0: Mild infection
// 1.0 - 3.0: Moderate infection
// 3.0 - 6.0: Severe infection
// 6.0 - 10.0: Critical infection (high mortality)
```

---

## Symptoms and Syndromes

### Symptom Development

```rust
// Daily symptom check
fn check_symptoms(individual, bacteria_idx, rng) {
    if individual.level[bacteria_idx] < symptom_threshold {
        return;  // Level too low for symptoms
    }
    
    if individual.symptoms[bacteria_idx] {
        return;  // Already symptomatic
    }
    
    // Probability increases with level
    let level = individual.level[bacteria_idx];
    let symptom_prob = symptom_probability_base[bacteria_idx] * level;
    
    if rng.gen::<f64>() < symptom_prob {
        individual.symptoms[bacteria_idx] = true;
        assign_syndrome(individual, bacteria_idx, rng);
    }
}
```

### Syndromes

```rust
pub enum Syndrome {
    UTI,                  // Urinary tract infection
    Pneumonia,            // Lower respiratory
    SkinSoftTissue,       // Cellulitis, abscess, wound
    Bacteremia,           // Bloodstream
    GI,                   // Gastrointestinal
    Meningitis,           // CNS infection
    BoneJoint,            // Osteomyelitis, septic arthritis
    Endocarditis,         // Heart valve infection
    Other,
}
```

### Syndrome Assignment

```rust
fn assign_syndrome(individual, bacteria_idx, rng) {
    let bacteria_name = BACTERIA_LIST[bacteria_idx];
    
    // Get syndrome distribution for this bacteria
    let syndrome_probs = bacteria_syndrome_distribution[bacteria_name];
    
    // E.g., E. coli: {UTI: 0.6, Bacteremia: 0.2, GI: 0.1, Other: 0.1}
    // S. aureus: {SkinSoftTissue: 0.5, Bacteremia: 0.2, Pneumonia: 0.15, ...}
    
    // Random selection weighted by probabilities
    let roll = rng.gen::<f64>();
    let mut cumulative = 0.0;
    
    for (syndrome, prob) in syndrome_probs.iter() {
        cumulative += prob;
        if roll < cumulative {
            individual.syndrome[bacteria_idx] = Some(syndrome);
            return;
        }
    }
    
    // Default
    individual.syndrome[bacteria_idx] = Some(Syndrome::Other);
}
```

---

## Sepsis

### Sepsis Development

```rust
// Sepsis check for symptomatic infections
fn check_sepsis(individual, bacteria_idx, time_step, rng) {
    if !individual.symptoms[bacteria_idx] {
        return;  // Must be symptomatic first
    }
    
    if individual.sepsis[bacteria_idx] {
        return;  // Already septic
    }
    
    // Base sepsis probability (bacteria-specific)
    let base_prob = sepsis_probability[bacteria_idx];
    
    // Risk factors
    let mut multiplier = 1.0;
    
    // Age
    if individual.age > 65 { multiplier *= 2.0; }
    if individual.age > 80 { multiplier *= 1.5; }  // Compound
    
    // Infection level
    multiplier *= individual.level[bacteria_idx] / 3.0;  // Higher level = higher risk
    
    // Time since symptom onset
    let symptom_duration = time_step - individual.date_symptoms_started[bacteria_idx]?;
    if symptom_duration > 3 { multiplier *= 1.5; }  // Untreated longer
    
    // Treatment status
    if !individual.on_antibiotics[bacteria_idx] {
        multiplier *= 2.0;  // Untreated
    }
    
    // Hospital status (paradoxically both protective and risky)
    // In hospital = sicker patients but more monitoring
    
    let sepsis_prob = base_prob * multiplier;
    
    if rng.gen::<f64>() < sepsis_prob {
        individual.sepsis[bacteria_idx] = true;
        individual.date_sepsis_started[bacteria_idx] = Some(time_step);
    }
}
```

### Sepsis Effects

```rust
// Sepsis dramatically increases mortality
if individual.sepsis[bacteria_idx] {
    // Base mortality hazard is much higher
    let sepsis_mortality_multiplier = 5.0;  // 5x baseline
    
    // Septic shock further increases risk
    if individual.septic_shock[bacteria_idx] {
        sepsis_mortality_multiplier *= 2.0;
    }
    
    // Drug effect is more critical
    // Under-treatment is especially dangerous
}
```

---

## Infection Clearance

### Clearance Mechanisms

1. **Immune-mediated**: Natural resolution
2. **Drug-assisted**: Antibiotic treatment
3. **Death**: Infection clears because host dies

### Immune-Mediated Clearance

```rust
// Natural clearance without drugs
fn immune_clearance(individual, bacteria_idx, time_step, rng) {
    // Immune strength
    let immune_effect = individual.immune_response_strength[bacteria_idx];
    
    // Clearance probability
    let base_clearance = bacteria_natural_clearance_rate[bacteria_idx];
    
    // Lower levels clear more easily
    let level_factor = 1.0 / (1.0 + individual.level[bacteria_idx]);
    
    let clearance_prob = base_clearance * immune_effect * level_factor;
    
    if rng.gen::<f64>() < clearance_prob {
        clear_infection(individual, bacteria_idx, time_step, "immune");
    }
}
```

### Drug-Assisted Clearance

```rust
// Clearance with antibiotic help
fn drug_clearance(individual, bacteria_idx, time_step) {
    let total_drug_effect = calculate_total_drug_effect(individual, bacteria_idx);
    
    if total_drug_effect > drug_clearance_threshold {
        // Reduce level
        let reduction = total_drug_effect * drug_clearance_coefficient;
        individual.level[bacteria_idx] -= reduction;
        
        // Check if cleared
        if individual.level[bacteria_idx] < INFECTION_EPS {
            clear_infection(individual, bacteria_idx, time_step, "drug");
        }
    }
}

const INFECTION_EPS: f64 = 0.01;  // Below this = cleared
```

### Clearance Actions

```rust
fn clear_infection(individual, bacteria_idx, time_step, mechanism: &str) {
    // Reset infection state
    individual.level[bacteria_idx] = 0.0;
    individual.symptoms[bacteria_idx] = false;
    individual.sepsis[bacteria_idx] = false;
    individual.syndrome[bacteria_idx] = None;
    
    // Record clearance
    individual.date_infection_cleared[bacteria_idx] = Some(time_step);
    individual.clearance_mechanism[bacteria_idx] = Some(mechanism.to_string());
    
    // Immune memory (reduced susceptibility to reinfection)
    let memory_boost = immune_memory_strength[bacteria_idx];
    individual.immune_memory[bacteria_idx] = 
        (individual.immune_memory[bacteria_idx] + memory_boost).min(1.0);
}
```

---

## Mortality

### Mortality Hazard Calculation

```rust
fn calculate_infection_mortality_hazard(individual, bacteria_idx) -> f64 {
    if individual.level[bacteria_idx] <= 0.0 {
        return 0.0;
    }
    
    // Base mortality by bacteria
    let base_hazard = bacteria_infection_mortality[bacteria_idx];
    
    let mut hazard = base_hazard;
    
    // Level effect (higher = more dangerous)
    hazard *= (individual.level[bacteria_idx] / 5.0).powf(1.5);
    
    // Sepsis multiplier
    if individual.sepsis[bacteria_idx] {
        hazard *= sepsis_mortality_multiplier;  // ~5x
    }
    
    // Age effect
    let age_mult = age_mortality_multiplier(individual.age);
    hazard *= age_mult;
    
    // Treatment effect (protection)
    if individual.on_treatment_for[bacteria_idx] {
        let effective_treatment = individual.treatment_is_effective[bacteria_idx];
        if effective_treatment {
            hazard *= 0.2;  // 80% reduction
        } else {
            hazard *= 0.7;  // 30% reduction (some benefit)
        }
    }
    
    // Syndrome effect
    match individual.syndrome[bacteria_idx] {
        Some(Syndrome::Bacteremia) => hazard *= 1.5,
        Some(Syndrome::Meningitis) => hazard *= 2.0,
        Some(Syndrome::Endocarditis) => hazard *= 2.5,
        Some(Syndrome::Pneumonia) => hazard *= 1.2,
        _ => {}
    }
    
    hazard
}
```

### Daily Mortality Roll

```rust
fn mortality_check(individual, time_step, rng) {
    if individual.date_of_death.is_some() {
        return;  // Already dead
    }
    
    // Sum mortality hazards from all infections
    let mut total_hazard = 0.0;
    for bacteria_idx in 0..N_BACTERIA {
        total_hazard += calculate_infection_mortality_hazard(individual, bacteria_idx);
    }
    
    // Add toxicity hazard
    total_hazard += individual.current_toxicity_hazard;
    
    // Add background mortality (age-based)
    total_hazard += background_mortality_rate(individual.age);
    
    // Mortality roll
    if rng.gen::<f64>() < total_hazard {
        individual.date_of_death = Some(time_step);
        individual.cause_of_death = determine_cause_of_death(individual);
    }
}
```

### Cause of Death Attribution

```rust
fn determine_cause_of_death(individual) -> String {
    // Find highest-risk factor
    let mut max_hazard = 0.0;
    let mut cause = "background";
    
    for bacteria_idx in 0..N_BACTERIA {
        let hazard = calculate_infection_mortality_hazard(individual, bacteria_idx);
        if hazard > max_hazard {
            max_hazard = hazard;
            cause = BACTERIA_LIST[bacteria_idx];
        }
    }
    
    if individual.current_toxicity_hazard > max_hazard {
        cause = "toxicity";
    }
    
    cause.to_string()
}
```

---

## Reinfection and Immunity

### Reinfection Window

```rust
// Period of reduced susceptibility after clearance
const REINFECTION_PROTECTION_WINDOW: usize = 30;  // days

fn check_reinfection_protection(individual, bacteria_idx, time_step) -> f64 {
    match individual.date_infection_cleared[bacteria_idx] {
        Some(cleared_day) => {
            let days_since = time_step - cleared_day;
            if days_since < REINFECTION_PROTECTION_WINDOW {
                // Exponential decay of protection
                let protection = (-days_since as f64 / 15.0).exp();
                return protection;  // 0 to 1
            }
            0.0
        }
        None => 0.0,
    }
}
```

### Immune Memory

```rust
// Long-term partial immunity
// Updated when infections clear
individual.immune_memory[bacteria_idx]

// Used in acquisition probability
let acquisition_rate = base_rate * (1.0 - individual.immune_memory[bacteria_idx] * 0.5);
```

### Cross-Immunity

Some bacteria confer partial immunity to related species:
```rust
// If had K. pneumoniae, slightly protected against K. oxytoca
let cross_immunity_groups = [
    ("k_pneumoniae", "k_oxytoca", 0.3),
    ("e_coli", "k_pneumoniae", 0.1),
    ("s_aureus", "s_epidermidis", 0.2),
    // ... etc
];
```

---

## Configuration Parameters

### Section A: Acquisition Rates

```rust
// Community acquisition
bacteria_{bacteria}_community_acquisition_rate
bacteria_{bacteria}_hospital_acquisition_rate
bacteria_{bacteria}_icu_acquisition_multiplier
bacteria_{bacteria}_contact_transmission_rate

// Seasonal factors
bacteria_{bacteria}_winter_multiplier
bacteria_{bacteria}_summer_multiplier
```

### Section B: Progression Parameters

```rust
bacteria_{bacteria}_growth_rate
bacteria_{bacteria}_symptom_threshold
bacteria_{bacteria}_symptom_probability
bacteria_{bacteria}_sepsis_probability
bacteria_{bacteria}_natural_clearance_rate
```

### Section F: Mortality Parameters

```rust
bacteria_{bacteria}_infection_mortality_rate
sepsis_mortality_multiplier
age_mortality_multiplier_65_plus
age_mortality_multiplier_80_plus
effective_treatment_mortality_reduction
```

### Section G: Regional Variation

```rust
{region}_bacteria_{bacteria}_prevalence_multiplier
{region}_hospital_acquisition_multiplier
{region}_community_infection_baseline
```

### Key Thresholds

```rust
// Infection levels
INFECTION_EPS = 0.01                // Below = cleared
SYMPTOM_THRESHOLD = 0.1             // Level for symptoms
SHEDDING_THRESHOLD = 0.05           // Level to transmit
MAX_INFECTION_LEVEL = 10.0          // Cap

// Timing (days)
REINFECTION_PROTECTION_WINDOW = 30
SYMPTOM_TO_SEPSIS_MIN_DAYS = 1
UNTREATED_SEPSIS_MAX_DAYS = 7
```

---

## Debugging Infection Issues

### Common Problems

1. **Infections not starting**: Check acquisition rates, regional multipliers
2. **Infections not clearing**: Check drug effects, immune response
3. **Too much mortality**: Check mortality rates, sepsis multipliers
4. **Unrealistic patterns**: Check seasonal factors, age distributions

### Key Logging Points

```rust
// Acquisition (rules/mod.rs)
// Search for: "initiate_infection" or "acquisition"

// Progression (rules/mod.rs)  
// Search for: "progress_infection" or "level"

// Clearance (rules/mod.rs)
// Search for: "clear_infection" or "clearance"

// Mortality (rules/mod.rs)
// Search for: "mortality" or "death"
```

### Diagnostic Outputs

```rust
// Track infections
println!("Individual {} bacteria {} level: {}", 
    individual.id, BACTERIA_LIST[b], individual.level[b]);

// Track symptoms
println!("Symptoms: {}, Sepsis: {}", 
    individual.symptoms[b], individual.sepsis[b]);

// Track clearance
if individual.level[b] == 0.0 && prev_level > 0.0 {
    println!("Cleared! Mechanism: {:?}", individual.clearance_mechanism[b]);
}
```
