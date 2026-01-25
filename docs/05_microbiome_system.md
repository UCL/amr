# Microbiome System

> **Source Files**: 
> - `src/rules/mod.rs` → colonization, HGT, clearance
> - `src/config.rs` → microbiome parameters (Sections A, D)
> - `src/simulation/population.rs` → microbiome arrays

This document explains how bacterial colonization works, how resistance spreads within the microbiome, and how carriage affects infections.

---

## Table of Contents
1. [Microbiome Model Overview](#microbiome-model-overview)
2. [Colonization (Carriage)](#colonization-carriage)
3. [Colonization Acquisition](#colonization-acquisition)
4. [Colonization Clearance](#colonization-clearance)
5. [Resistance in Microbiome](#resistance-in-microbiome)
6. [Horizontal Gene Transfer (HGT)](#horizontal-gene-transfer-hgt)
7. [Microbiome → Infection Seeding](#microbiome--infection-seeding)
8. [Body Sites](#body-sites)
9. [Drug Effects on Microbiome](#drug-effects-on-microbiome)
10. [Configuration Parameters](#configuration-parameters)

---

## Microbiome Model Overview

### Conceptual Model

```
┌───────────────────────────────────────────────────────────────────────────┐
│                        Microbiome Dynamics                                │
├───────────────────────────────────────────────────────────────────────────┤
│                                                                           │
│   Environment        Colonization            Outcomes                    │
│                                                                           │
│   ┌──────────┐      ┌──────────────┐        ┌───────────────────┐       │
│   │Community │ ──▶  │presence=true │ ───▶   │Seeding to infection│       │
│   │Hospital  │      │level: 0.0-1.0│        │HGT resistance spread│      │
│   │Contact   │      │microbiome_r[]│        │Competitive exclusion│      │
│   └──────────┘      └──────────────┘        └───────────────────┘       │
│        │                  │                         │                    │
│   Acquisition        Natural Decay              Clearance                │
│        │            Drug Perturbation           Drug Effect              │
│        ▼                  ▼                         ▼                    │
│   Resistant         Resistance                 Loss of                   │
│   Strains           Transfer (HGT)             Carriage                  │
│                                                                           │
└───────────────────────────────────────────────────────────────────────────┘
```

### Key Variables (per bacteria)

| Variable | Type | Purpose |
|----------|------|---------|
| `presence_microbiome[b]` | bool | Is bacteria colonizing? |
| `microbiome_level[b]` | f64 | Carriage density (0.0-1.0) |
| `microbiome_r[b][d]` | f64 | Resistance level in carriage |
| `date_microbiome_acquired[b]` | Option<usize> | When colonization started |
| `microbiome_acquisition_source[b]` | Option<String> | How acquired |

---

## Colonization (Carriage)

### What is Colonization?

**Colonization** (also called **carriage**) means bacteria are present in/on the body but NOT causing active infection:
- Present in gut, skin, respiratory tract
- Asymptomatic
- Can persist for months or years
- Can transmit to others
- Can seed future infections

### Examples

| Bacteria | Common Carriage Site | Carriage Rate |
|----------|---------------------|---------------|
| S. aureus | Anterior nares (nose) | ~30% population |
| E. coli | Gut | ~100% population |
| K. pneumoniae | Gut | ~5-30% population |
| S. pneumoniae | Nasopharynx | ~5-20% (higher in children) |
| E. faecalis | Gut | ~80% population |
| P. aeruginosa | Rare in healthy | ~5% (higher in hospital) |

### Carriage State

```rust
// Individual has carriage tracking arrays
presence_microbiome: [bool; N_BACTERIA],
microbiome_level: [f64; N_BACTERIA],  // 0.0 to 1.0

// Example state
individual.presence_microbiome[e_coli_idx] = true;
individual.microbiome_level[e_coli_idx] = 0.5;  // Moderate carriage
individual.presence_microbiome[p_aeruginosa_idx] = false;
```

---

## Colonization Acquisition

### Acquisition Sources

1. **Community**: Normal environmental exposure
2. **Hospital**: Nosocomial acquisition (higher rates)
3. **Contact**: From colonized/infected individuals
4. **Vertical**: Birth (some bacteria)

### Community Acquisition

```rust
fn community_colonization_acquisition(individual, bacteria_idx, time_step, rng) {
    if individual.presence_microbiome[bacteria_idx] {
        return;  // Already colonized
    }
    
    // Base colonization rate
    let base_rate = bacteria_colonization_rate[bacteria_idx];
    
    // Regional prevalence
    let regional_mult = regional_carriage_prevalence[individual.region_cur_in][bacteria_idx];
    
    // Age effect (children often colonized more easily)
    let age_mult = colonization_age_multiplier(individual.age);
    
    // Antibiotic use disrupts normal flora, allows pathogens
    let antibiotic_mult = if individual.on_any_antibiotic {
        antibiotic_colonization_susceptibility_multiplier  // e.g., 2.0
    } else {
        1.0
    };
    
    let acquisition_rate = base_rate * regional_mult * age_mult * antibiotic_mult;
    
    if rng.gen::<f64>() < acquisition_rate {
        initiate_colonization(individual, bacteria_idx, time_step, "community");
    }
}
```

### Hospital Acquisition

```rust
// Much higher acquisition in hospital environment
if individual.hospital_status > 0 {
    let hospital_rate = hospital_colonization_rate[bacteria_idx];
    
    // ICU particularly high risk
    let ward_mult = match individual.hospital_ward {
        Ward::ICU => 3.0,
        Ward::Surgical => 2.0,
        _ => 1.5,
    };
    
    // Certain bacteria are hospital-adapted
    let nosocomial_bacteria = ["p_aeruginosa", "a_baumannii", "c_difficile", "s_maltophilia"];
    let bacteria_mult = if nosocomial_bacteria.contains(&BACTERIA_LIST[bacteria_idx]) {
        3.0
    } else {
        1.0
    };
    
    // Hospital flora is often resistant
    let resistant_acquisition_probability = 0.5;  // 50% chance acquired strain is resistant
    
    if rng.gen::<f64>() < hospital_rate * ward_mult * bacteria_mult {
        initiate_colonization(individual, bacteria_idx, time_step, "hospital");
        
        // Potentially acquire resistant strain
        if rng.gen::<f64>() < resistant_acquisition_probability {
            acquire_resistant_colonization(individual, bacteria_idx);
        }
    }
}
```

### Contact Acquisition

```rust
// From household members or contacts
for contact in individual.contact_network.iter() {
    if contact.presence_microbiome[bacteria_idx] {
        // Contact must have sufficient carriage level
        let contact_level = contact.microbiome_level[bacteria_idx];
        if contact_level > shedding_threshold {
            let transmission_rate = colonization_contact_transmission_rate[bacteria_idx];
            
            if rng.gen::<f64>() < transmission_rate * contact_level {
                initiate_colonization(individual, bacteria_idx, time_step, "contact");
                
                // May acquire contact's resistant strain
                if contact.has_resistant_microbiome(bacteria_idx) {
                    transfer_colonization_resistance(contact, individual, bacteria_idx);
                }
            }
        }
    }
}
```

### Colonization Initiation

```rust
fn initiate_colonization(individual, bacteria_idx, time_step, source: &str) {
    individual.presence_microbiome[bacteria_idx] = true;
    individual.microbiome_level[bacteria_idx] = initial_colonization_level;  // e.g., 0.1
    individual.date_microbiome_acquired[bacteria_idx] = Some(time_step);
    individual.microbiome_acquisition_source[bacteria_idx] = Some(source.to_string());
    
    // Initialize microbiome resistance (usually susceptible at acquisition)
    // Unless acquired from resistant source
    for drug_idx in 0..N_DRUGS {
        if !acquired_from_resistant_source {
            individual.microbiome_r[bacteria_idx][drug_idx] = 0.0;
        }
    }
}
```

---

## Colonization Clearance

### Natural Clearance

```rust
fn natural_colonization_clearance(individual, bacteria_idx, time_step, rng) {
    if !individual.presence_microbiome[bacteria_idx] {
        return;
    }
    
    // Base clearance rate (bacteria-specific, per day)
    let clearance_rate = bacteria_colonization_clearance_rate[bacteria_idx];
    
    // Duration effect (longer carriage = harder to clear)
    let duration = time_step - individual.date_microbiome_acquired[bacteria_idx]?;
    let duration_penalty = 1.0 / (1.0 + duration as f64 / 30.0);
    
    // Level effect (lower levels clear more easily)
    let level_factor = 1.0 - individual.microbiome_level[bacteria_idx];
    
    // Competitive exclusion (other bacteria can displace)
    let competition_factor = calculate_competitive_pressure(individual, bacteria_idx);
    
    let total_clearance_rate = clearance_rate * duration_penalty * level_factor * competition_factor;
    
    if rng.gen::<f64>() < total_clearance_rate {
        clear_colonization(individual, bacteria_idx, time_step);
    }
}
```

### Clearance Action

```rust
fn clear_colonization(individual, bacteria_idx, time_step) {
    individual.presence_microbiome[bacteria_idx] = false;
    individual.microbiome_level[bacteria_idx] = 0.0;
    individual.date_microbiome_cleared[bacteria_idx] = Some(time_step);
    
    // Resistance is lost when bacteria is cleared
    for drug_idx in 0..N_DRUGS {
        individual.microbiome_r[bacteria_idx][drug_idx] = 0.0;
    }
}
```

### Colonization Level Dynamics

```rust
// Level fluctuates even during persistent carriage
fn update_colonization_level(individual, bacteria_idx, rng) {
    if !individual.presence_microbiome[bacteria_idx] {
        return;
    }
    
    let current_level = individual.microbiome_level[bacteria_idx];
    
    // Natural fluctuation
    let fluctuation = (rng.gen::<f64>() - 0.5) * 0.1;  // ±5%
    
    // Tendency toward equilibrium
    let equilibrium_level = bacteria_equilibrium_colonization_level[bacteria_idx];
    let pull_to_equilibrium = (equilibrium_level - current_level) * 0.05;
    
    // Drug effects (see later section)
    let drug_effect = calculate_drug_effect_on_microbiome(individual, bacteria_idx);
    
    // Update
    let new_level = current_level + fluctuation + pull_to_equilibrium - drug_effect;
    individual.microbiome_level[bacteria_idx] = new_level.clamp(0.0, 1.0);
    
    // If level drops very low, may clear
    if individual.microbiome_level[bacteria_idx] < colonization_clearance_threshold {
        if rng.gen::<f64>() < low_level_clearance_probability {
            clear_colonization(individual, bacteria_idx, time_step);
        }
    }
}
```

---

## Resistance in Microbiome

### Microbiome Resistance Tracking

```rust
// Resistance levels in colonizing bacteria (per bacteria, per drug)
microbiome_r: [[f64; N_DRUGS]; N_BACTERIA]

// This is SEPARATE from infection resistance
// Carriage can be resistant without active infection
// When infection occurs from carriage, resistance transfers
```

### Resistance Acquisition in Microbiome

```rust
fn microbiome_resistance_acquisition(individual, bacteria_idx, drug_idx, rng) {
    if !individual.presence_microbiome[bacteria_idx] {
        return;
    }
    
    // Current microbiome resistance
    let current_r = individual.microbiome_r[bacteria_idx][drug_idx];
    if current_r >= 1.0 {
        return;  // Already maximally resistant
    }
    
    // Resistance can increase through:
    // 1. Drug pressure (selection)
    // 2. HGT from other bacteria
    // 3. De novo mutation
    
    let mut resistance_increase = 0.0;
    
    // Drug pressure
    if individual.cur_use_drug[drug_idx] {
        let drug_level = individual.cur_level_drug[drug_idx];
        let selection_pressure = drug_level * microbiome_selection_coefficient[bacteria_idx][drug_idx];
        resistance_increase += selection_pressure * 0.01;  // Small daily increase
    }
    
    // De novo mutation (rare)
    if rng.gen::<f64>() < de_novo_mutation_rate_microbiome {
        resistance_increase += de_novo_mutation_increment;
    }
    
    // Apply increase
    individual.microbiome_r[bacteria_idx][drug_idx] = 
        (current_r + resistance_increase).min(1.0);
}
```

### Resistance Decay in Microbiome

```rust
// Without drug pressure, resistance may revert
fn microbiome_resistance_decay(individual, bacteria_idx, drug_idx, rng) {
    if individual.cur_use_drug[drug_idx] {
        return;  // No decay while under selection
    }
    
    let current_r = individual.microbiome_r[bacteria_idx][drug_idx];
    if current_r <= 0.0 {
        return;
    }
    
    // Fitness cost means resistant bacteria may be outcompeted
    let fitness_cost = resistance_fitness_cost[bacteria_idx][drug_idx];
    let decay_rate = microbiome_resistance_reversion_rate * fitness_cost;
    
    if rng.gen::<f64>() < decay_rate {
        let decay_amount = microbiome_resistance_decay_increment;
        individual.microbiome_r[bacteria_idx][drug_idx] = 
            (current_r - decay_amount).max(0.0);
    }
}
```

---

## Horizontal Gene Transfer (HGT)

### HGT Between Bacteria

Resistance genes can transfer between different bacteria species within the gut:

```rust
fn horizontal_gene_transfer(individual, time_step, rng) {
    // Check all pairs of colonizing bacteria
    for donor_idx in 0..N_BACTERIA {
        if !individual.presence_microbiome[donor_idx] {
            continue;
        }
        
        for recipient_idx in 0..N_BACTERIA {
            if donor_idx == recipient_idx {
                continue;
            }
            if !individual.presence_microbiome[recipient_idx] {
                continue;
            }
            
            // HGT probability depends on taxonomic relatedness
            let hgt_rate = hgt_rate_matrix[donor_idx][recipient_idx];
            
            // Also depends on physical proximity (same body site)
            // Gut-to-gut transfer more likely than gut-to-skin
            
            if rng.gen::<f64>() < hgt_rate {
                attempt_resistance_transfer(individual, donor_idx, recipient_idx, rng);
            }
        }
    }
}
```

### Transfer Mechanics

```rust
fn attempt_resistance_transfer(individual, donor_idx, recipient_idx, rng) {
    // Only transfer mobile resistance elements
    for mechanism_idx in 0..N_MECHANISMS {
        let mechanism = &MECHANISMS[mechanism_idx];
        
        // Check if mechanism is mobile (plasmid-borne, etc.)
        if !mechanism.is_mobile {
            continue;
        }
        
        // Check if donor has this mechanism
        let donor_has = individual.microbiome_mechanisms[donor_idx][mechanism_idx];
        if !donor_has {
            continue;
        }
        
        // Check if recipient already has it
        let recipient_has = individual.microbiome_mechanisms[recipient_idx][mechanism_idx];
        if recipient_has {
            continue;
        }
        
        // Transfer probability
        let transfer_prob = mechanism.transfer_probability;
        if rng.gen::<f64>() < transfer_prob {
            // Transfer mechanism
            individual.microbiome_mechanisms[recipient_idx][mechanism_idx] = true;
            
            // Transfer associated resistance
            for drug_idx in mechanism.conferred_resistance.iter() {
                let donor_r = individual.microbiome_r[donor_idx][*drug_idx];
                individual.microbiome_r[recipient_idx][*drug_idx] = 
                    individual.microbiome_r[recipient_idx][*drug_idx].max(donor_r);
            }
            
            record_hgt_event(donor_idx, recipient_idx, mechanism_idx, time_step);
        }
    }
}
```

### HGT Rate Matrix

Transfer rates vary by bacterial relatedness:

| Transfer | Rate | Notes |
|----------|------|-------|
| Within species | High | Same receptors |
| Within genus | Medium | Similar genetics |
| Within family | Low | Compatible plasmids |
| Between families | Very low | Rare events |
| Gram+ to Gram- | Minimal | Different cell walls |

```rust
// Example rates (per day, when both bacteria present)
hgt_rate_matrix[e_coli][k_pneumoniae] = 0.001;   // Same family (Enterobacteriaceae)
hgt_rate_matrix[e_coli][p_aeruginosa] = 0.0001;  // Different families
hgt_rate_matrix[e_coli][s_aureus] = 0.00001;     // Gram- to Gram+
hgt_rate_matrix[e_coli][e_coli] = 0.0;           // Same species (already have it)
```

---

## Microbiome → Infection Seeding

### Seeding Process

Colonizing bacteria can cause infections:

```rust
fn check_microbiome_seeding(individual, bacteria_idx, time_step, rng) {
    if !individual.presence_microbiome[bacteria_idx] {
        return;
    }
    
    if individual.level[bacteria_idx] > 0.0 {
        return;  // Already infected
    }
    
    // Base seeding rate
    let base_rate = microbiome_to_infection_seeding_rate[bacteria_idx];
    
    // Carriage level effect (higher level = more likely)
    let level_factor = individual.microbiome_level[bacteria_idx];
    
    // Risk factors
    let mut risk_mult = 1.0;
    
    // Hospital admission (procedures, stress)
    if individual.hospital_status > 0 {
        risk_mult *= hospital_seeding_multiplier;  // ~3x
    }
    
    // Surgical procedures
    if individual.recent_surgery {
        risk_mult *= surgery_seeding_multiplier;  // ~5x
    }
    
    // Catheters
    if individual.has_urinary_catheter && is_uti_bacteria(bacteria_idx) {
        risk_mult *= catheter_seeding_multiplier;  // ~10x
    }
    
    // Immunosuppression
    if individual.immunocompromised {
        risk_mult *= immunosuppression_seeding_multiplier;  // ~2x
    }
    
    let seeding_prob = base_rate * level_factor * risk_mult;
    
    if rng.gen::<f64>() < seeding_prob {
        seed_infection_from_microbiome(individual, bacteria_idx, time_step);
    }
}
```

### Seeding with Resistance

```rust
fn seed_infection_from_microbiome(individual, bacteria_idx, time_step) {
    // Initiate infection
    individual.level[bacteria_idx] = initial_infection_level;
    individual.date_last_infected[bacteria_idx] = Some(time_step);
    
    // CRITICAL: Infection inherits microbiome resistance
    for drug_idx in 0..N_DRUGS {
        let microbiome_resistance = individual.microbiome_r[bacteria_idx][drug_idx];
        individual.resistances[bacteria_idx][drug_idx].any_r = microbiome_resistance;
    }
    
    // Transfer mechanisms
    for mechanism_idx in 0..N_MECHANISMS {
        if individual.microbiome_mechanisms[bacteria_idx][mechanism_idx] {
            individual.resistance_mechanisms[bacteria_idx][mechanism_idx] = true;
        }
    }
    
    // Record source
    individual.infection_source[bacteria_idx] = Some("microbiome_seeding".to_string());
}
```

---

## Body Sites

### Site-Specific Colonization

Different bacteria colonize different body sites:

```rust
pub enum BodySite {
    Gut,
    Skin,
    RespiratoryUpper,
    RespiratoryLower,
    Urogenital,
    Oral,
}

// Typical colonization sites
let bacteria_primary_site: HashMap<&str, BodySite> = [
    ("e_coli", BodySite::Gut),
    ("k_pneumoniae", BodySite::Gut),
    ("s_aureus", BodySite::Skin),  // Also nares
    ("s_pneumoniae", BodySite::RespiratoryUpper),
    ("p_aeruginosa", BodySite::RespiratoryLower),  // In hospital
    ("h_influenzae", BodySite::RespiratoryUpper),
    ("c_difficile", BodySite::Gut),
    // ...
].iter().cloned().collect();
```

### HGT by Body Site

HGT is more likely between bacteria at the same body site:

```rust
fn hgt_site_compatibility(donor_site: BodySite, recipient_site: BodySite) -> f64 {
    if donor_site == recipient_site {
        1.0  // Same site = full rate
    } else if adjacent_sites(donor_site, recipient_site) {
        0.1  // Adjacent = reduced rate
    } else {
        0.01  // Distant = rare
    }
}
```

---

## Drug Effects on Microbiome

### Antibiotic Perturbation

Antibiotics affect microbiome composition:

```rust
fn antibiotic_effect_on_microbiome(individual, time_step) {
    for drug_idx in 0..N_DRUGS {
        if !individual.cur_use_drug[drug_idx] {
            continue;
        }
        
        let drug_level = individual.cur_level_drug[drug_idx];
        
        for bacteria_idx in 0..N_BACTERIA {
            if !individual.presence_microbiome[bacteria_idx] {
                continue;
            }
            
            // Drug's spectrum affects colonizing bacteria
            let potency = drug_bacteria_potency[drug_idx][bacteria_idx];
            
            // Resistance provides protection
            let resistance = individual.microbiome_r[bacteria_idx][drug_idx];
            
            // Net effect
            let killing_effect = potency * drug_level * (1.0 - resistance);
            
            // Reduce colonization level
            individual.microbiome_level[bacteria_idx] -= killing_effect * microbiome_drug_coefficient;
            
            // May clear colonization entirely
            if individual.microbiome_level[bacteria_idx] < colonization_clearance_threshold {
                clear_colonization(individual, bacteria_idx, time_step);
            }
        }
    }
}
```

### Selection for Resistance

```rust
// While susceptible bacteria die, resistant ones survive and proliferate
fn microbiome_selection(individual, drug_idx, rng) {
    for bacteria_idx in 0..N_BACTERIA {
        if !individual.presence_microbiome[bacteria_idx] {
            continue;
        }
        
        let current_r = individual.microbiome_r[bacteria_idx][drug_idx];
        
        // Drug kills susceptible population, resistant survive
        // This effectively increases resistance level
        let drug_level = individual.cur_level_drug[drug_idx];
        let selection_strength = drug_level * (1.0 - current_r);
        
        // Small increase in resistance with each day of selection
        let resistance_increase = selection_strength * selection_coefficient_daily;
        individual.microbiome_r[bacteria_idx][drug_idx] = 
            (current_r + resistance_increase).min(1.0);
    }
}
```

### Collateral Damage

Broad-spectrum antibiotics affect the whole microbiome:

```rust
// Broad-spectrum drugs damage protective flora
let broad_spectrum_drugs = ["meropenem", "piperacillin_tazobactam", "ciprofloxacin"];

if broad_spectrum_drugs.contains(&DRUG_SHORT_NAMES[drug_idx]) {
    // Increased susceptibility to C. difficile
    individual.cdiff_susceptibility_multiplier *= 2.0;
    
    // Ecological opportunity for pathogens
    for pathogen in opportunistic_pathogens {
        individual.pathogen_colonization_susceptibility[pathogen] *= 1.5;
    }
}
```

---

## Configuration Parameters

### Section A: Colonization Rates

```rust
// Acquisition
bacteria_{bacteria}_colonization_acquisition_rate
bacteria_{bacteria}_hospital_colonization_multiplier
bacteria_{bacteria}_contact_colonization_transmission_rate

// Clearance
bacteria_{bacteria}_colonization_clearance_rate
bacteria_{bacteria}_equilibrium_colonization_level

// Seeding
bacteria_{bacteria}_microbiome_to_infection_seeding_rate
hospital_seeding_multiplier
surgery_seeding_multiplier
catheter_seeding_multiplier
```

### Section D: HGT Parameters

```rust
// Transfer rates
hgt_base_rate_same_genus
hgt_base_rate_same_family
hgt_base_rate_different_family
hgt_base_rate_gram_switch

// Mechanism mobility
mechanism_{mechanism}_is_mobile
mechanism_{mechanism}_transfer_probability
mechanism_{mechanism}_fitness_cost
```

### Section H: Drug-Microbiome Interaction

```rust
microbiome_drug_killing_coefficient
microbiome_selection_coefficient_daily
microbiome_resistance_decay_rate
colonization_clearance_threshold
```

### Key Thresholds

```rust
COLONIZATION_SHEDDING_THRESHOLD = 0.2      // Minimum level to transmit
COLONIZATION_CLEARANCE_THRESHOLD = 0.05   // Below = may clear
INITIAL_COLONIZATION_LEVEL = 0.1          // Starting level
MAX_COLONIZATION_LEVEL = 1.0              // Cap
```

---

## Debugging Microbiome Issues

### Common Problems

1. **Colonization not persisting**: Check clearance rates
2. **Too much HGT**: Check HGT rate matrix values
3. **Seeding too frequent/rare**: Check seeding rates and risk multipliers
4. **Resistance not building**: Check selection coefficients

### Key Logging Points

```rust
// Colonization acquisition (rules/mod.rs)
// Search for: "colonization" or "microbiome_acquired"

// HGT (rules/mod.rs)
// Search for: "hgt" or "horizontal"

// Seeding (rules/mod.rs)
// Search for: "seeding" or "seed_infection"
```

### Useful Diagnostics

```rust
// Track colonization
println!("Microbiome presence: {:?}", individual.presence_microbiome);

// Track microbiome resistance
for b in 0..N_BACTERIA {
    if individual.presence_microbiome[b] {
        println!("Bacteria {} microbiome_r: {:?}", b, individual.microbiome_r[b]);
    }
}

// Track HGT events
println!("HGT event: {} -> {} (mechanism {})", donor, recipient, mechanism);
```
