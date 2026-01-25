# Simulation Flow

> **Source Files**: 
> - `src/main.rs` → entry point
> - `src/simulation/simulation.rs` → main loop
> - `src/rules/mod.rs` → apply_rules function
> - `src/simulation/population.rs` → data structures

This document explains how the simulation runs from start to finish.

---

## Table of Contents
1. [Execution Overview](#execution-overview)
2. [Initialization Phase](#initialization-phase)
3. [Main Simulation Loop](#main-simulation-loop)
4. [Daily Update Sequence](#daily-update-sequence)
5. [Rule Application Order](#rule-application-order)
6. [MajorityR Cache Updates](#majorityr-cache-updates)
7. [Output Generation](#output-generation)
8. [Performance Considerations](#performance-considerations)
9. [Debugging the Flow](#debugging-the-flow)

---

## Execution Overview

```
┌────────────────────────────────────────────────────────────────────────┐
│                        Simulation Lifecycle                            │
├────────────────────────────────────────────────────────────────────────┤
│                                                                        │
│   ┌──────────────┐    ┌──────────────┐    ┌──────────────────────┐    │
│   │ Initialization│ ─▶ │ Main Loop    │ ─▶ │ Output Generation   │    │
│   │              │    │ (daily)      │    │                      │    │
│   │ • Load config│    │ • For each   │    │ • Write summaries   │    │
│   │ • Create pop │    │   time_step  │    │ • Write traces      │    │
│   │ • Init state │    │ • Apply rules│    │ • Aggregate stats   │    │
│   │ • Seed R     │    │ • Update cache│   │                      │    │
│   │ • Init cache │    │ • Record stats│   │                      │    │
│   └──────────────┘    └──────────────┘    └──────────────────────┘    │
│                                                                        │
│   ~1 second           ~hours to days      ~seconds                     │
│                                                                        │
└────────────────────────────────────────────────────────────────────────┘
```

### Key Timing

| Phase | Typical Duration | Notes |
|-------|------------------|-------|
| Initialization | <1 second | One-time setup |
| Main loop | Hours to days | 100 years × 365 days = 36,500 iterations |
| Per time step | ~10-50ms | Depends on population size |
| Output | Seconds | End of simulation |

---

## Initialization Phase

### Entry Point (main.rs)

```rust
fn main() {
    // 1. Parse command line arguments
    let args = parse_args();
    
    // 2. Load configuration
    let params = config::get_params();
    
    // 3. Initialize random number generator
    let seed = params.get("random_seed").unwrap_or(&0) as u64;
    let mut rng = if seed == 0 {
        StdRng::from_entropy()  // Random seed
    } else {
        StdRng::seed_from_u64(seed)  // Reproducible
    };
    
    // 4. Create simulation
    let mut simulation = Simulation::new(&params, &mut rng);
    
    // 5. Run simulation
    simulation.run(&params, &mut rng);
    
    // 6. Generate output
    simulation.write_output(&params);
}
```

### Population Initialization

```rust
impl Simulation {
    pub fn new(params: &Params, rng: &mut StdRng) -> Self {
        let pop_size = *params.get("population_size").unwrap_or(&100000.0) as usize;
        
        // Create population
        let mut population = Vec::with_capacity(pop_size);
        
        for id in 0..pop_size {
            let individual = Individual::new(id, params, rng);
            population.push(individual);
        }
        
        // Initialize MajorityR cache
        let majority_r_cache = MajorityRCache::new(&population, params);
        
        Simulation {
            population,
            majority_r_cache,
            time_step: 0,
            start_time_step: get_start_time_step(params),
            end_time_step: get_end_time_step(params),
            // ... other fields
        }
    }
}
```

### Individual Initialization

```rust
impl Individual {
    pub fn new(id: usize, params: &Params, rng: &mut StdRng) -> Self {
        // Assign demographics
        let age = sample_initial_age(params, rng);
        let sex_at_birth = sample_sex(rng);
        let region_living = sample_region(params, rng);
        let perceived_penicillin_allergy = rng.gen::<f64>() < 
            *params.get("penicillin_allergy_prevalence").unwrap_or(&0.1);
        
        // Initialize all arrays to default values
        let level = [0.0; N_BACTERIA];
        let symptoms = [false; N_BACTERIA];
        let resistances = [[Resistance::default(); N_DRUGS]; N_BACTERIA];
        // ... etc
        
        // Seed baseline colonization
        let presence_microbiome = initialize_baseline_colonization(region_living, params, rng);
        
        // Seed baseline resistance based on region
        let resistances = initialize_regional_resistance(region_living, params, rng);
        
        Individual {
            id,
            age,
            sex_at_birth,
            region_living,
            region_cur_in: region_living,
            perceived_penicillin_allergy,
            hospital_status: 0,
            level,
            symptoms,
            resistances,
            presence_microbiome,
            // ... all other fields
        }
    }
}
```

### Initial Resistance Seeding

```rust
fn initialize_regional_resistance(region: usize, params: &Params, rng: &mut StdRng) 
    -> [[Resistance; N_DRUGS]; N_BACTERIA] 
{
    let mut resistances = [[Resistance::default(); N_DRUGS]; N_BACTERIA];
    
    for bacteria_idx in 0..N_BACTERIA {
        for drug_idx in 0..N_DRUGS {
            let bacteria = BACTERIA_LIST[bacteria_idx];
            let drug = DRUG_SHORT_NAMES[drug_idx];
            
            // Get regional baseline
            let key = format!("{}_bacteria_{}_drug_{}_baseline_r", 
                             REGIONS[region], bacteria, drug);
            let baseline = *params.get(&key).unwrap_or(&0.0);
            
            // Apply some stochastic variation
            let variation = (rng.gen::<f64>() - 0.5) * 0.1;  // ±5%
            let initial_r = (baseline + variation).clamp(0.0, 1.0);
            
            resistances[bacteria_idx][drug_idx].any_r = initial_r;
        }
    }
    
    resistances
}
```

---

## Main Simulation Loop

### Loop Structure

```rust
impl Simulation {
    pub fn run(&mut self, params: &Params, rng: &mut StdRng) {
        // Main time loop
        while self.time_step < self.end_time_step {
            // Progress indicator
            if self.time_step % 365 == 0 {
                let year = self.time_step / 365 + 1930;
                println!("Simulating year {}", year);
            }
            
            // Process each individual
            for individual in self.population.iter_mut() {
                apply_rules(individual, self.time_step, params, 
                           &self.majority_r_cache, rng);
            }
            
            // Update population-level cache
            if self.time_step % self.cache_refresh_interval == 0 {
                self.majority_r_cache.refresh(&self.population, params);
            }
            
            // Record statistics
            self.record_daily_stats();
            
            // Advance time
            self.time_step += 1;
        }
    }
}
```

### Time Step Interpretation

```rust
// Time step = day number since simulation start
// With 1930 start:
//   time_step 0     = Jan 1, 1930
//   time_step 365   = Jan 1, 1931
//   time_step 36500 = ~Jan 1, 2030 (100 years)

fn time_step_to_year(time_step: usize) -> usize {
    1930 + (time_step / 365)
}

fn time_step_to_date(time_step: usize) -> (usize, usize, usize) {
    let year = 1930 + (time_step / 365);
    let day_of_year = time_step % 365;
    // Approximate month/day (ignoring leap years for simplicity)
    let month = day_of_year / 30 + 1;
    let day = day_of_year % 30 + 1;
    (year, month, day)
}
```

---

## Daily Update Sequence

### The apply_rules Function

```rust
pub fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    params: &Params,
    majority_r_cache: &MajorityRCache,
    rng: &mut StdRng,
) {
    // Skip dead individuals
    if individual.date_of_death.is_some() {
        return;
    }
    
    // ========================================
    // PHASE 1: AGE AND DEMOGRAPHICS
    // ========================================
    update_age(individual, time_step);
    
    // ========================================
    // PHASE 2: LOCATION AND HOSPITALIZATION
    // ========================================
    update_hospitalization(individual, time_step, params, rng);
    update_location(individual, time_step, params, rng);
    
    // ========================================
    // PHASE 3: INFECTION ACQUISITION
    // ========================================
    for bacteria_idx in 0..N_BACTERIA {
        if individual.level[bacteria_idx] == 0.0 {
            // May acquire new infection
            community_acquisition(individual, bacteria_idx, time_step, params, rng);
            hospital_acquisition(individual, bacteria_idx, time_step, params, rng);
            contact_acquisition(individual, bacteria_idx, time_step, params, rng);
            microbiome_seeding(individual, bacteria_idx, time_step, params, rng);
        }
    }
    
    // ========================================
    // PHASE 4: INFECTION PROGRESSION
    // ========================================
    for bacteria_idx in 0..N_BACTERIA {
        if individual.level[bacteria_idx] > 0.0 {
            progress_infection(individual, bacteria_idx, time_step, params, rng);
            check_symptoms(individual, bacteria_idx, time_step, params, rng);
            check_sepsis(individual, bacteria_idx, time_step, params, rng);
        }
    }
    
    // ========================================
    // PHASE 5: DRUG SELECTION AND TREATMENT
    // ========================================
    for bacteria_idx in 0..N_BACTERIA {
        if individual.symptoms[bacteria_idx] && !individual.on_treatment_for[bacteria_idx] {
            select_and_start_drug(individual, bacteria_idx, time_step, params, rng);
        }
    }
    
    // ========================================
    // PHASE 6: DRUG EFFECTS
    // ========================================
    update_drug_levels(individual, time_step, params);
    apply_drug_effects(individual, time_step, params);
    assess_treatment_failure(individual, time_step, params, rng);
    update_toxicity(individual, time_step, params);
    
    // ========================================
    // PHASE 7: INFECTION CLEARANCE
    // ========================================
    for bacteria_idx in 0..N_BACTERIA {
        if individual.level[bacteria_idx] > 0.0 {
            check_clearance(individual, bacteria_idx, time_step, params, rng);
        }
    }
    
    // ========================================
    // PHASE 8: RESISTANCE DYNAMICS
    // ========================================
    for bacteria_idx in 0..N_BACTERIA {
        for drug_idx in 0..N_DRUGS {
            de_novo_emergence(individual, bacteria_idx, drug_idx, time_step, params, rng);
            resistance_reversion(individual, bacteria_idx, drug_idx, time_step, params, rng);
        }
    }
    apply_resistance_floors(individual, time_step, params);
    
    // ========================================
    // PHASE 9: MICROBIOME DYNAMICS
    // ========================================
    for bacteria_idx in 0..N_BACTERIA {
        update_colonization(individual, bacteria_idx, time_step, params, rng);
        colonization_acquisition(individual, bacteria_idx, time_step, params, rng);
        colonization_clearance(individual, bacteria_idx, time_step, params, rng);
    }
    horizontal_gene_transfer(individual, time_step, params, rng);
    
    // ========================================
    // PHASE 10: MORTALITY
    // ========================================
    check_mortality(individual, time_step, params, rng);
}
```

---

## Rule Application Order

### Why Order Matters

The order of operations affects outcomes:

1. **Acquisition before progression**: New infections start at initial level
2. **Progression before symptoms**: Level must increase before symptoms appear
3. **Symptoms before treatment**: Must be symptomatic to get treatment
4. **Treatment before clearance**: Drugs reduce level, then clearance check
5. **Resistance before clearance**: Resistance affects drug efficacy

### Phase Dependencies

```
┌─────────────────────────────────────────────────────────────────────┐
│                     Daily Update Dependencies                        │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│   Age/Demographics ──┐                                              │
│                      ▼                                              │
│   Location ─────────▶ Hospitalization                               │
│                          │                                          │
│                          ▼                                          │
│   Infection Acquisition (community, hospital, contact, seeding)     │
│                          │                                          │
│                          ▼                                          │
│   Infection Progression (level increase, symptoms, sepsis)          │
│                          │                                          │
│                          ▼                                          │
│   Drug Selection (only if symptomatic)                              │
│                          │                                          │
│                          ▼                                          │
│   Drug Effects (levels, activity, toxicity)                         │
│        │                                                            │
│        ├─────────────────┐                                          │
│        ▼                 ▼                                          │
│   Clearance        Resistance Dynamics                              │
│        │                 │                                          │
│        └─────────────────┘                                          │
│                          │                                          │
│                          ▼                                          │
│   Microbiome Dynamics (colonization, HGT)                           │
│                          │                                          │
│                          ▼                                          │
│   Mortality Check                                                   │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

---

## MajorityR Cache Updates

### When Cache Updates

```rust
// Cache refreshed periodically, not every time step
let refresh_interval = *params.get("majority_r_refresh_interval_days").unwrap_or(&7.0) as usize;

if self.time_step % refresh_interval == 0 {
    self.majority_r_cache.refresh(&self.population, params);
}
```

### Refresh Process

```rust
impl MajorityRCache {
    pub fn refresh(&mut self, population: &[Individual], params: &Params) {
        let sample_size = *params.get("majority_r_sample_size").unwrap_or(&100.0) as usize;
        
        // For each bacteria-drug pair
        for bacteria_idx in 0..N_BACTERIA {
            for drug_idx in 0..N_DRUGS {
                // Find individuals with this bacteria
                let with_bacteria: Vec<&Individual> = population.iter()
                    .filter(|ind| ind.level[bacteria_idx] > 0.0)
                    .collect();
                
                if with_bacteria.is_empty() {
                    self.cache[bacteria_idx][drug_idx] = 0.0;
                    continue;
                }
                
                // Sample up to sample_size individuals
                let sample: Vec<&Individual> = with_bacteria.iter()
                    .take(sample_size)
                    .cloned()
                    .collect();
                
                // Calculate median resistance
                let mut r_values: Vec<f64> = sample.iter()
                    .map(|ind| ind.resistances[bacteria_idx][drug_idx].any_r)
                    .collect();
                r_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
                
                let median = if r_values.len() % 2 == 0 {
                    (r_values[r_values.len()/2 - 1] + r_values[r_values.len()/2]) / 2.0
                } else {
                    r_values[r_values.len()/2]
                };
                
                self.cache[bacteria_idx][drug_idx] = median;
            }
        }
    }
}
```

---

## Output Generation

### During Simulation

```rust
fn record_daily_stats(&mut self) {
    // Aggregate counts
    let mut infected_count = [0usize; N_BACTERIA];
    let mut symptomatic_count = [0usize; N_BACTERIA];
    let mut septic_count = [0usize; N_BACTERIA];
    let mut on_treatment_count = [0usize; N_BACTERIA];
    let mut deaths_today = 0;
    
    for individual in &self.population {
        if individual.date_of_death.is_some() {
            if individual.date_of_death == Some(self.time_step) {
                deaths_today += 1;
            }
            continue;
        }
        
        for b in 0..N_BACTERIA {
            if individual.level[b] > 0.0 { infected_count[b] += 1; }
            if individual.symptoms[b] { symptomatic_count[b] += 1; }
            if individual.sepsis[b] { septic_count[b] += 1; }
            if individual.on_treatment_for[b] { on_treatment_count[b] += 1; }
        }
    }
    
    // Store in history
    self.daily_stats.push(DailyStats {
        time_step: self.time_step,
        infected_count,
        symptomatic_count,
        septic_count,
        on_treatment_count,
        deaths_today,
        // ...
    });
}
```

### End of Simulation

```rust
impl Simulation {
    pub fn write_output(&self, params: &Params) {
        // Write summary CSV
        self.write_summary_csv(params);
        
        // Write detailed resistance data
        if *params.get("output_detailed_resistance").unwrap_or(&0.0) > 0.0 {
            self.write_resistance_csv(params);
        }
        
        // Write individual traces if enabled
        if *params.get("output_individual_traces").unwrap_or(&0.0) > 0.0 {
            self.write_individual_traces(params);
        }
    }
    
    fn write_summary_csv(&self, params: &Params) {
        // Aggregate by year
        let mut file = File::create("simulation_summary.csv").unwrap();
        
        writeln!(file, "year,bacteria,total_infections,symptomatic_infections,\
                       deaths,mean_resistance_ampicillin,mean_resistance_ciprofloxacin,...").unwrap();
        
        for year in 1930..2030 {
            let start_step = (year - 1930) * 365;
            let end_step = start_step + 365;
            
            // Aggregate statistics for this year
            // ...
        }
    }
}
```

---

## Performance Considerations

### Bottlenecks

| Operation | Complexity | Notes |
|-----------|------------|-------|
| apply_rules per individual | O(N_BACTERIA × N_DRUGS) | ~2000 operations |
| Population loop | O(population_size) | 100,000 individuals |
| HGT check | O(N_BACTERIA²) | ~1500 pairs |
| MajorityR refresh | O(population_size) | Sampled |

### Optimization Strategies

```rust
// 1. Early exit for dead individuals
if individual.date_of_death.is_some() {
    return;
}

// 2. Skip uninfected in infection-specific loops
if individual.level[bacteria_idx] == 0.0 {
    continue;
}

// 3. Cache frequently accessed parameters
let growth_rates: [f64; N_BACTERIA] = precompute_growth_rates(params);

// 4. Use arrays instead of HashMaps for hot paths
// individual.level[bacteria_idx] instead of individual.level.get(&bacteria_name)

// 5. Parallel processing (if thread-safe)
use rayon::prelude::*;
self.population.par_iter_mut().for_each(|individual| {
    apply_rules(individual, ...);
});
```

### Memory Considerations

```rust
// Individual size estimate
// 60+ arrays × (N_BACTERIA or N_DRUGS) × 8 bytes
// ~60 × 52 × 8 = ~25 KB per individual
// 100,000 individuals = ~2.5 GB

// Consider using smaller types if memory constrained
// f64 → f32 halves memory
// bool array → bitset
```

---

## Debugging the Flow

### Tracing Individual

```rust
// Add debug output for specific individual
fn apply_rules(individual: &mut Individual, ...) {
    let trace_id = 12345;  // ID to trace
    let trace = individual.id == trace_id;
    
    if trace {
        println!("=== Time step {} Individual {} ===", time_step, individual.id);
        println!("  Infections: {:?}", individual.level.iter()
            .enumerate()
            .filter(|(_, &l)| l > 0.0)
            .collect::<Vec<_>>());
    }
    
    // ... rest of function
    
    if trace {
        println!("  After update - Infections: {:?}", ...);
    }
}
```

### Time Step Breakpoints

```rust
// Stop at specific time step for debugging
fn run(&mut self, ...) {
    while self.time_step < self.end_time_step {
        // Breakpoint for specific year
        if self.time_step == 36500 {  // Year 2030
            println!("Reached 2030 - breaking for inspection");
            self.print_summary();
            break;  // Or pause for interactive debugging
        }
        
        // ... normal processing
    }
}
```

### Statistics Validation

```rust
// Sanity checks during simulation
fn validate_state(&self) {
    let alive_count = self.population.iter()
        .filter(|ind| ind.date_of_death.is_none())
        .count();
    
    assert!(alive_count > 0, "All individuals died!");
    
    for individual in &self.population {
        for b in 0..N_BACTERIA {
            assert!(individual.level[b] >= 0.0, "Negative infection level!");
            assert!(individual.level[b] <= MAX_INFECTION_LEVEL, "Excessive infection level!");
            
            for d in 0..N_DRUGS {
                assert!(individual.resistances[b][d].any_r >= 0.0, "Negative resistance!");
                assert!(individual.resistances[b][d].any_r <= 1.0, "Resistance > 1.0!");
            }
        }
    }
}
```

### Logging Levels

```rust
// Control verbosity
enum LogLevel {
    Error,    // Only errors
    Warning,  // Errors + warnings
    Info,     // + yearly summaries
    Debug,    // + daily stats
    Trace,    // + individual events
}

static LOG_LEVEL: LogLevel = LogLevel::Info;

macro_rules! log_debug {
    ($($arg:tt)*) => {
        if LOG_LEVEL >= LogLevel::Debug {
            println!($($arg)*);
        }
    };
}
```
