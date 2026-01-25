# Resistance System

> **Source Files**: 
> - `src/simulation/population.rs` → `Resistance` struct, `ResistanceMechanism` enum
> - `src/rules/mod.rs` → resistance update logic
> - `src/simulation/simulation.rs` → `MajorityRCache`
> - `src/config.rs` → resistance parameters (Section H)

This document explains how antimicrobial resistance is modeled, tracked, and evolves.

---

## Table of Contents
1. [Resistance Level Model](#resistance-level-model)
2. [Resistance Mechanisms](#resistance-mechanisms)
3. [MajorityRCache System](#majoritycache-system)
4. [Resistance Acquisition Sources](#resistance-acquisition-sources)
5. [De Novo Emergence](#de-novo-emergence)
6. [Horizontal Gene Transfer (HGT)](#horizontal-gene-transfer-hgt)
7. [Resistance Reversion](#resistance-reversion)
8. [Resistance Floor System](#resistance-floor-system)
9. [Cross-Resistance Groups](#cross-resistance-groups)
10. [Configuration Parameters](#configuration-parameters)

---

## Resistance Level Model

### The `Resistance` Struct
```rust
pub struct Resistance {
    pub microbiome_r: f64,  // Resistance in colonizing bacteria (0.0-1.0)
    pub test_r: f64,        // Resistance as detected by testing
    pub activity_r: f64,    // Effective resistance for drug activity
    pub any_r: f64,         // Resistance in any fraction of bacteria (0.0-1.0)
    pub majority_r: f64,    // Resistance in >50% of bacteria (0.0-1.0)
}
```

### Interpretation of Levels

| Level | Clinical Meaning | Drug Effectiveness |
|-------|------------------|-------------------|
| 0.0 | Fully susceptible | Full activity |
| 0.1-0.3 | Minor resistance | Slightly reduced activity |
| 0.4-0.6 | Moderate resistance | Substantially reduced activity |
| 0.7-0.9 | High resistance | Minimal activity |
| 1.0 | Fully resistant | No activity |

### any_r vs majority_r

- **`any_r`**: Any resistant bacteria present (even 1%)
  - Affects treatment success probability
  - More sensitive to resistance emergence
  
- **`majority_r`**: >50% of bacteria resistant
  - Determines if resistance will persist after treatment
  - Feeds into population-level tracking (MajorityRCache)
  - Only set when resistance is dominant

**Update Rule**:
```rust
// When majority_r is non-zero, it equals any_r
if majority_r > 0.0 {
    majority_r == any_r  // Always true
}
```

### How Drug Activity is Calculated
```rust
effective_activity = base_potency * drug_level * (1.0 - any_r)
```

Example:
- Potency = 0.8, Drug level = 10.0, any_r = 0.3
- Activity = 0.8 × 10.0 × 0.7 = 5.6 (reduced from 8.0)

---

## Resistance Mechanisms

### Mechanism Types
```rust
pub enum ResistanceMechanism {
    ESBL,                      // Extended-spectrum beta-lactamase
    Carbapenemase,             // Carbapenem-hydrolyzing enzymes
    AmpC,                      // AmpC beta-lactamase
    SixteenSMethyltransferase, // Aminoglycoside resistance
    Qnr,                       // Quinolone resistance protein
    EffluxOverexpression,      // Multi-drug efflux pumps
    ErmMethylation,            // Macrolide resistance
    VanType,                   // Glycopeptide resistance
    MecA,                      // Methicillin resistance
    ReducedPermeability,       // Outer membrane changes
    TargetSiteMutation,        // Point mutations
}
```

### Mechanism-Bacteria Compatibility

| Mechanism | Compatible Bacteria Groups |
|-----------|---------------------------|
| ESBL, Carbapenemase, AmpC | Enterobacterales, Non-fermenters, Enteric pathogens |
| 16S Methyltransferase, Qnr | Enterobacterales, Non-fermenters, Enteric, Fastidious |
| Efflux, Reduced Permeability, Target Mutation | All bacteria |
| Erm Methylation, VanType, MecA | Gram-positives only |

### Mechanism-Drug Mapping

**ESBL** affects:
- 3rd/4th gen cephalosporins (ceftriaxone, ceftazidime, cefepime)
- Aztreonam
- NOT carbapenems (ESBL are inhibited by carbapenems)

**Carbapenemase** affects:
- All beta-lactams including carbapenems
- Usually not affected by beta-lactamase inhibitors

**VanType** (VanA/VanB) affects:
- Vancomycin
- Teicoplanin (VanA only)

**MecA** affects:
- All beta-lactams (encodes altered PBP2a)
- MRSA = S. aureus with MecA

### Enhancement Multiplier

Mechanisms provide multiplicative boosts to resistance:
```rust
enhancement = mechanism_enhancement_multiplier[mech]  // e.g., 1.5

// If mechanism present and enhancement > current resistance:
if enhancement > any_r {
    any_r = enhancement  // Mechanism guarantees minimum resistance
}
```

---

## MajorityRCache System

### Purpose
The `MajorityRCache` tracks population-level resistance prevalence for sampling when new infections are acquired.

### Structure
```rust
struct MajorityRCache {
    buckets: Vec<Bucket>,         // [region × hospital × bacteria × drug]
    world_buckets: Vec<Bucket>,   // [bacteria × drug] fallback
    bucket_threshold_met: Vec<bool>,
}
```

### Bucket Contents
```rust
struct Bucket {
    positive_samples: u32,   // Count of individuals with majority_r > 0
    total_samples: u32,      // Total individuals sampled
    sum_positive: f64,       // Sum of majority_r values for positives
}
```

### Sampling Process

When a new infection is acquired:
```rust
fn sample(region, hospital, bacteria, drug, rng) -> Option<f64> {
    let idx = index(region, hospital, bacteria, drug);
    
    // Use local bucket if enough data, else world bucket
    let probability = if bucket_threshold_met[idx] {
        buckets[idx].probability()  // positive_samples / total_samples
    } else {
        world_buckets[bacteria][drug].probability()
    };
    
    // Roll for resistance
    if rng.gen::<f64>() < probability {
        // Draw resistance level from distribution
        return Some(buckets[idx].draw_positive(rng));
    }
    
    Some(0.0)  // Susceptible
}
```

### Cache Updates
Updated periodically (every `majority_r_cache_update_interval` days):
```rust
for individual in population {
    for bacteria in infected_bacteria {
        for drug in all_drugs {
            let r_level = individual.resistances[b][d].majority_r;
            cache.add_sample(region, hospital, bacteria, drug, r_level);
        }
    }
}
```

---

## Resistance Acquisition Sources

### Source Types
```rust
pub enum ResistanceAcquisitionType {
    AtInfectionCommunity,  // Acquired resistant strain from community
    AtInfectionEnv,        // Environmental acquisition (food, water)
    AtInfectionTB,         // TB-specific
    Hgt,                   // Horizontal gene transfer from microbiome
    FromMicrobiomeR,       // Direct transfer from colonizing bacteria
    DeNovoInfection,       // Emerged during treatment
}
```

### At-Infection Assignment

When infection is acquired, resistance is sampled:
```rust
// 1. Sample from MajorityRCache
let level = majority_r_cache.sample(region, hospital, bacteria, drug, rng);

// 2. Apply resistance floor for rare bacteria
let floor = calculate_resistance_floor(bacteria, drug, current_day);
let level_with_floor = level.max(floor);

// 3. Clamp to max
let clamped = level_with_floor.min(max_resistance_level).max(0.0);

// 4. Assign
resistances[b][d].any_r = clamped;
resistances[b][d].majority_r = clamped;
```

---

## De Novo Emergence

### When Does It Occur?

Resistance can emerge during treatment when:
1. Individual has active infection
2. Drug is being used
3. Drug has suboptimal activity
4. Bacteria is susceptible (room to gain resistance)

### Emergence Calculation

```rust
// In rules/mod.rs, emergence logic
let emergence_rate = drug_bacteria_emergence_rate[d][b];
let drug_pressure = cur_level_drug[d];
let susceptible_fraction = 1.0 - any_r;

// Probability of emergence this timestep
let p_emerge = emergence_rate * drug_pressure * susceptible_fraction 
               * emergence_multiplier;

if rng.gen::<f64>() < p_emerge {
    // Emergence occurred
    let new_r = any_r + emergence_increment;
    any_r = new_r.min(max_resistance_level);
    how_resistance_acquired[b][d] = Some(DeNovoInfection);
}
```

### Mechanism-Based Emergence

Emergence can also assign specific mechanisms:
```rust
// If ESBL emergence occurs for cephalosporin
if mechanism_can_emerge(ESBL, bacteria, drug) {
    resistance_mechanisms[bacteria][ESBL] = true;
    
    // Cross-resistance to related drugs
    for related_drug in esbl_affected_drugs {
        any_r[related_drug] = any_r[related_drug].max(esbl_enhancement);
    }
}
```

---

## Horizontal Gene Transfer (HGT)

### Overview
Resistance genes can transfer between bacteria via:
- Conjugation (plasmids)
- Transformation (free DNA)
- Transduction (phages)

### HGT Process

```rust
// Daily HGT check for each bacteria pair
for donor_bacteria in BACTERIA_LIST {
    for recipient_bacteria in BACTERIA_LIST {
        if !can_hgt(donor, recipient) { continue; }
        
        let donor_has_r = microbiome_r[donor][drug] > threshold;
        let recipient_susceptible = microbiome_r[recipient][drug] < threshold;
        
        if donor_has_r && recipient_susceptible {
            let hgt_rate = get_hgt_rate(donor, recipient, drug);
            
            if rng.gen::<f64>() < hgt_rate {
                // Transfer resistance
                microbiome_r[recipient][drug] = donor_level * transfer_fraction;
            }
        }
    }
}
```

### HGT Compatibility

HGT is restricted by:
1. **Bacteria group**: Gram-negatives transfer with Gram-negatives
2. **Carriage compartment**: Both must colonize same compartment
3. **Mechanism compatibility**: Some mechanisms transfer more readily

### Configuration
```
// HGT rate parameters (per day)
hgt_base_rate_intra_group = 0.0001     // Within same bacteria group
hgt_base_rate_inter_group = 0.00001   // Between groups (rarer)
```

---

## Resistance Reversion

### Fitness Cost Model

Resistant bacteria often have reduced fitness when drugs absent:
```rust
// Daily reversion check (when not on drug)
if cur_level_drug[d] < DRUG_EFFECT_THRESHOLD {
    let reversion_rate = fitness_cost_reversion_rate[b][d];
    let current_r = any_r;
    
    // Probability of reversion
    let p_revert = reversion_rate * current_r;  // Higher R = faster revert
    
    if rng.gen::<f64>() < p_revert {
        any_r *= (1.0 - reversion_decrement);
        
        if any_r < reversion_floor {
            any_r = 0.0;  // Complete reversion
        }
    }
}
```

### Mechanism Persistence

Some mechanisms are more stable than others:
- **Plasmid-borne** (ESBL): Can be lost, moderate reversion
- **Chromosomal** (AmpC): Very stable, slow reversion  
- **Point mutations**: Variable, depends on fitness cost

---

## Resistance Floor System

### Purpose
For rare bacteria with very few infections (e.g., S. maltophilia, E. faecium), cache-based sampling may not sustain realistic resistance levels. The floor system provides minimum resistance values.

### Floor Calculation

```rust
pub fn calculate_resistance_floor(bacteria: &str, drug: &str, current_day: i32) -> f64 {
    // 1. Check if floors enabled for this bacteria
    if !bacteria_resistance_floor_enabled(bacteria) {
        return 0.0;
    }
    
    // 2. Get drug class (floors are per-class, not per-drug)
    let drug_class = get_drug_class(drug);  // e.g., "carbapenems"
    
    // 3. Get drug class introduction day
    let intro_day = get_drug_class_introduction_day(drug_class);
    
    // 4. No floor before drug existed
    if current_day < intro_day {
        return 0.0;
    }
    
    // 5. Get target floor for this bacteria-class
    let target = get_resistance_floor_target(bacteria, drug);
    
    // 6. Calculate ramp (linear from intro to full floor)
    let ramp_years = get_resistance_floor_ramp_years(bacteria);
    let ramp_days = ramp_years * 365.0;
    let days_since_intro = current_day - intro_day;
    let ramp_fraction = (days_since_intro / ramp_days).min(1.0);
    
    target * ramp_fraction
}
```

### Configured Floors

**S. maltophilia** (intrinsic L1/L2 beta-lactamases):
| Drug Class | Floor |
|------------|-------|
| Carbapenems | 0.98 |
| Penicillins | 0.95 |
| Cephalosporins 1/2 | 0.95 |
| Cephalosporins 3/4 | 0.75 |
| TMP-SMX | 0.15 (preferred therapy) |

**E. faecium** (intrinsic cephalosporin resistance):
| Drug Class | Floor |
|------------|-------|
| All cephalosporins | 0.99 |
| Glycopeptides (VRE) | 0.35 |
| Fluoroquinolones | 0.65 |
| Penicillins | 0.0 (acquired) |

---

## Cross-Resistance Groups

### Purpose
Resistance to one drug often confers resistance to related drugs.

### Group Definitions
```rust
// In config.rs, cross-resistance groups
let esbl_group = ["ceftriaxone", "ceftazidime", "cefepime", "aztreonam"];
let carbapenem_group = ["meropenem", "imipenem_c", "ertapenem"];
let fluoroquinolone_group = ["ciprofloxacin", "levofloxacin", "moxifloxacin"];
let aminoglycoside_group = ["gentamicin", "tobramycin", "amikacin"];
```

### Cross-Resistance Application

When resistance increases for one drug in a group:
```rust
// Acquire resistance to ceftriaxone
resistances[b][ceftriaxone].any_r = 0.8;

// Apply cross-resistance
for related_drug in esbl_group {
    let cross_r = 0.8 * cross_resistance_fraction;  // e.g., 0.8 * 0.9 = 0.72
    resistances[b][related_drug].any_r = 
        resistances[b][related_drug].any_r.max(cross_r);
}
```

---

## Configuration Parameters

### Section H: Resistance Parameters (config.rs)

```rust
// Emergence rates
drug_{drug}_for_bacteria_{bacteria}_resistance_emergence_rate_per_day_baseline

// Mechanism multipliers
mechanism_{mechanism}_for_bacteria_{bacteria}_emergence_multiplier

// Reversion
fitness_cost_reversion_rate_baseline
reversion_decrement_per_event

// HGT
hgt_base_rate_same_compartment
hgt_base_rate_different_compartment
hgt_{mechanism}_transfer_probability

// Cross-resistance
cross_resistance_{group}_fraction

// Floors (Section D.2)
resistance_floor_feature_enabled
bacteria_{bacteria}_resistance_floor_enabled
bacteria_{bacteria}_resistance_floor_ramp_years
bacteria_{bacteria}_{drug_class}_resistance_floor
```

### Key Global Parameters
```rust
max_resistance_level = 1.0           // Ceiling for all resistance
mechanism_assignment_probability_on_any_r_gain = 0.3
emergence_multiplier_under_drug_pressure = 1.5
```

---

## Debugging Resistance Issues

### Common Problems

1. **Resistance too high**: Check emergence rates, floor values
2. **Resistance too low**: Check cache sampling, floor configuration
3. **Resistance not persisting**: Check majority_r assignment, reversion rates
4. **Cross-resistance spreading**: Check group definitions, fractions

### Key Log Points
```rust
// In rules/mod.rs, around line 3480+
// Resistance assignment at infection
// Emergence events
// HGT events
// Reversion events
```

### Useful Diagnostics
```rust
// Print cache state
majority_r_cache.print_summary();

// Track individual
if individual.id == target_id {
    println!("R[b][d] = {}", resistances[b][d].any_r);
}
```

