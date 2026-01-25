# Drug Treatment System

> **Source Files**: 
> - `src/rules/mod.rs` → drug selection, effects, toxicity
> - `src/config.rs` → drug parameters (Sections C, D, E)
> - `src/simulation/population.rs` → `DRUG_SHORT_NAMES`

This document explains how antibiotics are selected, how they affect infections, and how treatment decisions are made.

---

## Table of Contents
1. [Drug List and Classes](#drug-list-and-classes)
2. [Drug Selection Algorithm](#drug-selection-algorithm)
3. [Pharmacokinetics](#pharmacokinetics)
4. [Drug Activity Calculation](#drug-activity-calculation)
5. [Treatment Duration](#treatment-duration)
6. [Treatment Failure](#treatment-failure)
7. [Drug Toxicity](#drug-toxicity)
8. [Drug Interactions](#drug-interactions)
9. [Historical Drug Availability](#historical-drug-availability)
10. [Configuration Parameters](#configuration-parameters)

---

## Drug List and Classes

### Complete Drug List (52 drugs)
```rust
pub const DRUG_SHORT_NAMES: &[&str] = &[
    // Sulfonamides
    "sulfanilamide",
    
    // Penicillins
    "penicilling", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin",
    
    // Cephalosporins
    "cephalexin", "cefazolin", "cefuroxime",           // 1st/2nd gen
    "ceftriaxone", "ceftazidime", "cefepime", "ceftaroline",  // 3rd/4th gen
    
    // Carbapenems
    "meropenem", "imipenem_c", "ertapenem",
    
    // Monobactams
    "aztreonam",
    
    // Macrolides
    "erythromycin", "azithromycin", "clarithromycin",
    
    // Lincosamides
    "clindamycin",
    
    // Aminoglycosides
    "gentamicin", "tobramycin", "amikacin",
    
    // Fluoroquinolones
    "ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin",
    
    // Tetracyclines
    "tetracycline", "doxycycline", "minocycline",
    
    // Glycopeptides
    "vancomycin", "teicoplanin", "dalbavancin",
    
    // Oxazolidinones
    "linezolid", "tedizolid",
    
    // Other
    "quinu_dalfo", "trim_sulf", "chlorampheni", "nitrofurantoin",
    "retapamulin", "fusidic_a", "metronidazole", "furazolidone", "rifampicin",
    
    // BL/BLI Combinations
    "amoxicillin_clavulanate", "piperacillin_tazobactam",
    "ampicillin_sulbactam", "ticarcillin_clavulanate",
    "ceftazidime_avibactam", "meropenem_vaborbactam",
    
    // Polymyxins
    "colistin",
];
```

### Drug Classes
```rust
let polymyxins = ["colistin"];
let penicillins = ["penicilling", "ampicillin", "amoxicillin", "piperacillin", 
                   "ticarcillin", "amoxicillin_clavulanate", "piperacillin_tazobactam",
                   "ampicillin_sulbactam", "ticarcillin_clavulanate"];
let cephalosporins_1_2 = ["cephalexin", "cefazolin", "cefuroxime"];
let cephalosporins_3_4 = ["ceftriaxone", "ceftazidime", "cefepime", "ceftaroline",
                          "ceftazidime_avibactam"];
let carbapenems = ["meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"];
let macrolides = ["erythromycin", "azithromycin", "clarithromycin"];
let aminoglycosides = ["gentamicin", "tobramycin", "amikacin"];
let fluoroquinolones = ["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"];
let tetracyclines = ["tetracycline", "doxycycline", "minocycline"];
let glycopeptides = ["vancomycin", "teicoplanin", "dalbavancin"];
let oxazolidinones = ["linezolid", "tedizolid"];
```

---

## Drug Selection Algorithm

### When Selection Occurs

Drug selection is triggered when:
1. Individual has symptomatic infection
2. No current treatment, OR current treatment failing
3. Drug selection hasn't occurred today

### Selection Process Overview

```rust
fn select_drug_for_bacteria(individual, bacteria_idx, syndrome, rng) -> Option<usize> {
    // 1. Get candidate drugs
    let candidates = get_candidate_drugs(bacteria_idx, syndrome);
    
    // 2. Filter by availability
    let available = filter_by_availability(candidates, region, time_step);
    
    // 3. Filter by allergy
    let safe = filter_by_allergy(available, individual.perceived_penicillin_allergy);
    
    // 4. Score drugs
    let scored = score_drugs(safe, individual, bacteria_idx);
    
    // 5. Probabilistic selection (weighted by score)
    let selected = weighted_random_selection(scored, rng);
    
    selected
}
```

### Drug Scoring Components

```rust
fn score_drug(drug_idx, individual, bacteria_idx) -> f64 {
    let mut score = 1.0;
    
    // Base potency against this bacteria
    let potency = drug_bacteria_potency[drug_idx][bacteria_idx];
    score *= potency;
    
    // Initiation multiplier (prescribing preference)
    let init_mult = drug_bacteria_initiation_multiplier[drug_idx][bacteria_idx];
    score *= init_mult;
    
    // Resistance penalty
    let resistance = individual.resistances[bacteria_idx][drug_idx].any_r;
    score *= (1.0 - resistance);
    
    // Syndrome appropriateness
    let syndrome_mult = syndrome_drug_multiplier[syndrome][drug_idx];
    score *= syndrome_mult;
    
    // Regional availability
    let availability = region_drug_availability[region][drug_idx];
    score *= availability;
    
    // Previous failure penalty
    if had_recent_failure(individual, drug_idx, bacteria_idx) {
        score *= failure_retry_penalty;  // e.g., 0.1
    }
    
    score
}
```

### Penicillin Allergy Handling

```rust
const PENICILLIN_CLASS_DRUGS: &[&str] = &[
    "penicilling", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin",
    "amoxicillin_clavulanate", "ampicillin_sulbactam", 
    "piperacillin_tazobactam", "ticarcillin_clavulanate",
];

if individual.perceived_penicillin_allergy {
    // Skip all penicillin-class drugs
    for drug in PENICILLIN_CLASS_DRUGS {
        scores[drug_idx] = 0.0;
    }
}
```

### Probabilistic Selection

```rust
// Not deterministic "best drug" - weighted random
let total_score: f64 = scores.iter().sum();
let roll = rng.gen::<f64>() * total_score;

let mut cumulative = 0.0;
for (idx, score) in scores.iter().enumerate() {
    cumulative += score;
    if roll < cumulative {
        return Some(idx);
    }
}
```

---

## Pharmacokinetics

### Drug Level Model

```rust
// Standard dose
let standard_dose_level = 10.0;

// On dosing day
if taking_drug_today {
    cur_level_drug[d] = standard_dose_level * dose_multiplier;
}

// Decay between doses
else {
    let half_life = drug_half_life[d];  // in days
    let decay_factor = (-LN_2 / half_life).exp();
    cur_level_drug[d] *= decay_factor;
}
```

### Half-Lives by Drug Class

| Drug Class | Typical Half-Life | Notes |
|------------|-------------------|-------|
| Penicillins | 0.5-1.0 days | Short, frequent dosing |
| Cephalosporins | 1-2 days | Intermediate |
| Carbapenems | 0.5-1.0 days | Short |
| Fluoroquinolones | 4-6 days | Long-acting |
| Aminoglycosides | 2-3 days | Concentration-dependent |
| Glycopeptides | 4-6 days | Long |
| Azithromycin | 7-10 days | Very long tissue half-life |

### Level Interpretation

| Level | Interpretation |
|-------|---------------|
| 0.0-0.1 | Sub-therapeutic |
| 0.1-3.0 | Low therapeutic |
| 3.0-7.0 | Moderate therapeutic |
| 7.0-12.0 | Standard therapeutic |
| >12.0 | High dose |

---

## Drug Activity Calculation

### Basic Activity Formula

```rust
fn calculate_drug_activity(drug_idx, bacteria_idx, individual) -> f64 {
    let drug_level = individual.cur_level_drug[drug_idx];
    let potency = drug_bacteria_potency[drug_idx][bacteria_idx];
    let resistance = individual.resistances[bacteria_idx][drug_idx].any_r;
    let max_r = max_resistance_level;
    
    // Normalize resistance to 0-1 range
    let normalized_r = resistance / max_r;
    
    // Activity = potency × level × (1 - resistance)
    let activity = potency * drug_level * (1.0 - normalized_r);
    
    // Apply individual response variation
    activity * individual.drug_activity_response_multiplier[bacteria_idx]
}
```

### Activity Effects

```rust
// In infection clearance logic
let total_drug_effect = sum(calculate_drug_activity(d, b, individual) 
                            for d in active_drugs);

// Reduce infection level
let reduction = total_drug_effect * drug_clearance_coefficient;
individual.level[bacteria_idx] -= reduction;

// If level drops below threshold
if individual.level[bacteria_idx] < INFECTION_EPS {
    // Infection cleared (drug-assisted)
    record_clearance(DrugAssistedClearance);
}
```

### Potency Matrix

Potency values range 0.0 to 1.0+:
- **1.0+**: First-line therapy (excellent)
- **0.5-0.99**: Good activity (reliable option)
- **0.25-0.49**: Moderate (situational use)
- **0.05-0.24**: Poor (usually ineffective)
- **0.0-0.05**: No activity

Example potencies:
```
E. coli + Meropenem: 1.0 (excellent)
E. coli + Ampicillin: 0.6 (moderate due to resistance)
S. maltophilia + Meropenem: 0.05 (intrinsic resistance)
S. aureus + Vancomycin: 1.0 (excellent)
E. faecium + Ceftriaxone: 0.0 (intrinsic resistance)
```

---

## Treatment Duration

### Standard Durations

```rust
let treatment_duration = drug_treatment_duration[drug_idx][syndrome];

// Typical values (in days)
// UTI: 3-7 days
// Pneumonia: 5-10 days
// Skin/soft tissue: 7-14 days
// Bacteremia: 7-14 days
// Endocarditis: 28-42 days
```

### Early Stopping

Treatment may stop early due to:
1. Infection cleared
2. Patient non-adherence
3. Toxicity
4. Treatment failure (switch drug)

```rust
// Non-adherence model
let adherence_probability = get_adherence_probability(drug_idx, days_on_treatment);
if rng.gen::<f64>() > adherence_probability {
    // Stop drug early
    individual.cur_use_drug[drug_idx] = false;
    individual.drug_stopped_with_infection_day[bacteria_idx] = Some(time_step);
}
```

### Completion

```rust
let days_on_drug = time_step - individual.date_drug_initiated[drug_idx];
if days_on_drug >= treatment_duration {
    // Complete course
    individual.cur_use_drug[drug_idx] = false;
    // Don't reset levels immediately - let them decay naturally
}
```

---

## Treatment Failure

### Failure Definition

Treatment is considered failing when:
1. Day 3+ of treatment
2. Infection level not decreasing adequately
3. Or: New symptoms/worsening

### Failure Assessment

```rust
// On day 3 of treatment
if individual.days_on_current_treatment[bacteria_idx] == 3 {
    if !individual.treatment_failure_assessed[bacteria_idx] {
        let initial_level = individual.bacteria_level_at_drug_start[bacteria_idx];
        let current_level = individual.level[bacteria_idx];
        
        let improvement = (initial_level - current_level) / initial_level;
        
        if improvement < failure_threshold {  // e.g., < 0.3 (30% improvement)
            // Treatment failure
            record_failure(bacteria_idx, drug_idx);
            trigger_drug_switch(individual, bacteria_idx);
        }
        
        individual.treatment_failure_assessed[bacteria_idx] = true;
    }
}
```

### Post-Failure Drug Selection

After failure, the failing drug gets a strong penalty:
```rust
if time_step - individual.date_last_drug_failure[bacteria_idx] < failure_penalty_window {
    // Heavily penalize the failed drug
    scores[failed_drug_idx] *= 0.1;
    
    // Also penalize same-class drugs
    for drug in same_class {
        scores[drug] *= 0.3;
    }
}
```

---

## Drug Toxicity

### Toxicity Accumulation

```rust
// Each day while on drug
let toxicity_rate = drug_toxicity_accumulation_rate[drug_idx];
individual.drug_toxicity_reservoir[drug_idx] += 
    individual.cur_level_drug[drug_idx] * toxicity_rate;

// Natural decay
let decay_rate = drug_toxicity_decay_rate[drug_idx];
individual.drug_toxicity_reservoir[drug_idx] *= (1.0 - decay_rate);
```

### Toxicity Hazard

```rust
// Sum toxicity from all drugs
let total_toxicity: f64 = individual.drug_toxicity_reservoir.iter().sum();

// Convert to hazard (S-curve)
let toxicity_midpoint = toxicity_mortality_midpoint;  // e.g., 50
let toxicity_slope = toxicity_mortality_slope;        // e.g., 0.1

individual.current_toxicity_hazard = 
    1.0 / (1.0 + (-toxicity_slope * (total_toxicity - toxicity_midpoint)).exp());
```

### Mortality from Toxicity

```rust
// Daily mortality roll
if rng.gen::<f64>() < individual.current_toxicity_hazard * toxicity_mortality_scale {
    individual.date_of_death = Some(time_step);
    individual.cause_of_death = Some("toxicity".to_string());
}
```

### High-Toxicity Drugs

| Drug | Toxicity Profile |
|------|-----------------|
| Aminoglycosides | Nephrotoxicity, ototoxicity |
| Vancomycin | Nephrotoxicity, red man syndrome |
| Colistin | Nephrotoxicity, neurotoxicity |
| Linezolid | Myelosuppression, lactic acidosis |
| Chloramphenicol | Aplastic anemia |

---

## Drug Interactions

### Level Interactions

When drugs are co-administered, they can affect each other's levels:

```rust
// CYP450 interactions
if taking("rifampicin") && taking("clarithromycin") {
    // Rifampicin induces CYP3A4, reduces clarithromycin levels
    let multiplier = drug_level_multiplier_clarithromycin_when_coadministered_with_rifampicin;
    // typically 0.6 (40% reduction)
    cur_level_drug[clarithromycin_idx] *= multiplier;
}
```

### Configured Interactions

```rust
// In config.rs
drug_level_multiplier_levofloxacin_when_coadministered_with_rifampicin = 0.7
drug_level_multiplier_moxifloxacin_when_coadministered_with_rifampicin = 0.8
drug_level_multiplier_clarithromycin_when_coadministered_with_rifampicin = 0.6
drug_level_multiplier_ciprofloxacin_when_coadministered_with_erythromycin = 0.85
```

### Synergy/Antagonism

Not currently modeled directly, but could be added:
```rust
// Example synergy (not implemented)
if taking("beta_lactam") && taking("aminoglycoside") {
    activity *= synergy_multiplier;  // e.g., 1.3
}
```

---

## Historical Drug Availability

### Introduction Dates

Drugs are only available after their historical introduction:

```rust
pub static ref DRUG_INTRODUCTION_DATES: HashMap<&'static str, usize> = {
    let mut map = HashMap::new();
    map.insert("sulfanilamide", 2555);   // 1937
    map.insert("penicilling", 3555);     // 1942
    map.insert("ampicillin", 11315);     // 1961
    map.insert("ciprofloxacin", 20805);  // 1987
    map.insert("linezolid", 25550);      // 2000
    map.insert("ceftazidime_avibactam", 27740);  // 2006
    // ... etc
    map
};
```

### Availability Check

```rust
fn drug_available(drug_idx, time_step) -> bool {
    let drug_name = DRUG_SHORT_NAMES[drug_idx];
    match get_drug_introduction_time_step(drug_name) {
        Some(intro_day) => time_step >= intro_day,
        None => true,  // No restriction if not configured
    }
}

// In drug selection
for drug in candidates {
    if !drug_available(drug, time_step) {
        scores[drug] = 0.0;  // Not yet invented
    }
}
```

### Regional Availability

Different regions have different drug access:
```rust
// In config.rs
north_america_drug_linezolid_availability = 1.0
africa_drug_linezolid_availability = 0.1
asia_drug_tedizolid_availability = 0.3
```

---

## Configuration Parameters

### Section C: Drug Parameters

```rust
// Initiation and selection
drug_{drug}_for_bacteria_{bacteria}_initiation_multiplier
drug_{drug}_for_bacteria_{bacteria}_potency_when_no_r
drug_{drug}_half_life_days
drug_{drug}_standard_treatment_duration

// Toxicity
drug_{drug}_toxicity_accumulation_rate
drug_{drug}_toxicity_decay_rate
toxicity_mortality_midpoint
toxicity_mortality_slope

// Adherence
drug_{drug}_adherence_day_1_to_3
drug_{drug}_adherence_day_4_to_7
drug_{drug}_adherence_day_8_plus
```

### Section D: Interactions

```rust
drug_level_multiplier_{drug1}_when_coadministered_with_{drug2}
```

### Section E: Potency Matrix

```rust
drug_{drug}_for_bacteria_{bacteria}_potency_when_no_r
// ~2000 combinations (52 drugs × 39 bacteria)
```

### Section J: Availability

```rust
{region}_drug_{drug}_availability
```

### Key Globals

```rust
max_resistance_level = 1.0
drug_selection_score_min = 0.01
failure_assessment_day = 3
failure_improvement_threshold = 0.3
failure_penalty_window_days = 30
failure_retry_penalty = 0.1
```

---

## Debugging Drug Issues

### Common Problems

1. **Drug not being selected**: Check availability, potency, allergy
2. **Drug not working**: Check resistance, potency, levels
3. **Too much toxicity**: Check accumulation rates, decay rates
4. **Wrong drugs for syndrome**: Check syndrome_drug_multiplier

### Key Log Points

```rust
// Drug selection (rules/mod.rs)
// Search for: "drug selection" or "score"

// Drug effects (rules/mod.rs)
// Search for: "drug_effect" or "activity"

// Toxicity (rules/mod.rs)
// Search for: "toxicity"
```

### Useful Diagnostics

```rust
// Print drug scores
println!("Drug scores for bacteria {}: {:?}", bacteria_idx, scores);

// Track drug levels
println!("Drug levels: {:?}", individual.cur_level_drug);

// Check if drug available
println!("Drug {} available: {}", drug, drug_available(drug_idx, time_step));
```

