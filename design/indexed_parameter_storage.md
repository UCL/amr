# Indexed parameter storage design

This document sketches the refactor that will replace the runtime `HashMap<String, f64>` parameter table with cache-friendly indexed arrays. It translates the hotspots captured in `parameter_lookup_hotspots.md` into concrete data structures and initialization steps.

## Goals

1. **Eliminate per-timestep hash lookups** in `rules::apply_rules`, drug scoring, and other hot loops.
2. **Preserve existing configuration flexibility** (CSV/JSON loaders still populate the same logical parameter families).
3. **Introduce a clear, typed API** so that high-level code calls `store.drug_potency(b_idx, d_idx)` instead of formatting keys.
4. **Allow staged migration**: legacy helper functions (`get_global_param`, `get_bacteria_param`, `get_drug_param`) should keep working during the transition, forwarding into the new structures.

## Terminology and indexing

- `BacteriaIdx` – integer index into `BACTERIA_LIST` (`usize`). Number of entries: `NUM_BACTERIA = BACTERIA_LIST.len()`.
- `DrugIdx` – integer index into `DRUG_SHORT_NAMES` (`usize`). Number of entries: `NUM_DRUGS = DRUG_SHORT_NAMES.len()`.
- `RegionIdx` – integer index into the six geographic regions (`0..6`). We can reuse `Region as usize` because the enum is already ordered.
- `AgeBucketIdx` – index into six age buckets (`0-1`, `1-5`, `5-18`, `18-50`, `50-70`, `70+`). A helper already exists in `population.rs` to classify ages.
- `MechanismIdx` – index into `ResistanceMechanism::all()`.

Newtype wrappers (`struct BacteriaIdx(usize);`) would improve type safety but are optional for the first pass.

## Top-level container

Introduce a `ParameterStore` struct inside `config.rs` (or a new module) that owns all indexed data:

```rust
pub struct ParameterStore {
    pub globals: GlobalScalars,
    pub drug: DrugParameters,
    pub bacteria: BacteriaParameters,
    pub drug_bacteria: DrugBacteriaMatrix,
    pub region_bacteria: RegionBacteriaParameters,
    pub age_bacteria: AgeBacteriaParameters,
    pub age_region: AgeRegionParameters,
    pub hgt: HgtMatrix,
    pub resistance_mechanism: ResistanceMechanismParameters,
}
```

Each sub-struct encapsulates a dense allocation plus convenience methods.

### Global scalars

```rust
pub struct GlobalScalars {
    pub drug_base_initiation_rate_per_day: f32,
    pub drug_infection_present_multiplier: f32,
    pub double_dose_probability_if_identified_infection: f32,
    // … all other frequently-read scalars …
}
```

Populate with `f32` to reduce bandwidth unless a value genuinely needs double precision. Construction pulls values from the legacy map once at startup, supplying defaults where needed. Adding a `from_hashmap(&HashMap<String, f64>)` associated function keeps initialization contained.

### Drug-specific parameters

```rust
pub struct DrugParameters {
    pub initial_level: Vec<f32>,
    pub double_dose_multiplier: Vec<f32>,
    pub spectrum_breadth: Vec<f32>,
    pub toxicity_per_unit_level_per_day: Vec<f32>,
    pub half_life_days: Vec<f32>,
    pub availability: RegionDrugAvailability, // see below
}

impl DrugParameters {
    #[inline] pub fn initial_level(&self, drug: DrugIdx) -> f32 { self.initial_level[drug] }
    // … similar accessors …
}
```

`RegionDrugAvailability` becomes a `[NUM_REGIONS][NUM_DRUGS]` matrix of multipliers precomputed per region (including `home`). Time-aware logic (e.g., colistin hiatus) stays in Rust code; the base availability multipliers themselves come from this matrix.

### Bacteria-specific parameters

A mirror of the drug struct:

```rust
pub struct BacteriaParameters {
    pub acquisition_log_odds_baseline: Vec<f32>,
    pub log_odds_vaccinated: Vec<f32>,
    pub log_odds_microbiome_present: Vec<f32>,
    pub log_odds_hospital_acquired: Vec<f32>,
    pub microbiome_clearance_probability_per_day: Vec<f32>,
    pub immunity_effect_on_level_change: Vec<f32>,
    pub base_bacteria_level_change: Vec<f32>,
    pub max_level: Vec<f32>,
    pub immunity_base_response: Vec<f32>,
    pub immunity_increase_per_infection_day: Vec<f32>,
    pub immunity_increase_per_unit_level: Vec<f32>,
    pub immunity_age_modifier: Vec<f32>,
    pub immunity_immunodeficiency_modifier: Vec<f32>,
    pub max_immune_response: Vec<f32>,
    // … add other bacteria-only fields as they’re identified …
}
```

### Drug × bacteria matrices

Store the major cross-products as flat vectors for cache efficiency.

```rust
pub struct DrugBacteriaMatrix {
    pub potency_when_no_r: Vec<f32>,         // len = NUM_BACTERIA * NUM_DRUGS
    pub initiation_multiplier: Vec<f32>,    // same shape
    pub resistance_emergence_rate: Vec<f32>,
    pub mic_lt2_threshold: Vec<f32>,        // precomputed thresholds
}

impl DrugBacteriaMatrix {
    #[inline]
    pub fn potency(&self, bacteria: BacteriaIdx, drug: DrugIdx) -> f32 {
        let idx = bacteria * NUM_DRUGS + drug;
        self.potency_when_no_r[idx]
    }
}
```

`Simulation::new` already builds a potency matrix; the refactor will move that logic into `DrugBacteriaMatrix::from_hashmap` so both the store and the simulation share the same data.

### Region × bacteria and age modifiers

- `RegionBacteriaParameters` → `[NUM_REGIONS * NUM_BACTERIA]` floats for acquisition log-odds.
- `AgeBacteriaParameters` → `[NUM_BACTERIA * NUM_AGE_BUCKETS]`.
- `AgeRegionParameters` → `[NUM_REGIONS * NUM_AGE_BUCKETS]` for fallback multipliers.

We can package the age data in a helper struct:

```rust
pub struct AgeLogOddsTables {
    pub bacteria_specific: Vec<f32>, // bacteria-major order
    pub region_specific: Vec<f32>,   // region-major order
    pub default_by_age: [f32; NUM_AGE_BUCKETS],
}
```

### HGT matrix

`hgt_prob_{donor}_to_{recipient}` becomes a dense `Vec<f32>` sized `NUM_BACTERIA * NUM_BACTERIA` (optionally zero on diagonal). Accessed via `store.hgt.get(donor, recipient)`.

### Resistance mechanism parameters

Two vectors sized to the mechanism count:

```rust
pub struct ResistanceMechanismParameters {
    pub emergence_rate: Vec<f32>,
    pub enhancement_multiplier: Vec<f32>,
}
```

### Parameter cache rewrite

`ParameterKeyCache` loses the string HashMaps and becomes index math wrappers:

```rust
pub struct ParameterCache {
    pub drug_bacteria: DrugBacteriaMatrix,
    // …
}
```

In practice, `ParameterCache` may shrink to just hold precomputed derived values (e.g., WeightedIndex inputs). Most present-day fields become redundant once direct functions exist.

## Initialization flow

1. Load raw overrides as today (templates → JSON → CSV → etc.). They still populate a `HashMap<String, f64>` temporarily.
2. Call `ParameterStore::from_hashmap(&PARAMETERS)` during startup, producing the indexed representation.
3. Replace `lazy_static! { pub static ref PARAMETERS: HashMap<_, _> }` with `lazy_static! { pub static ref PARAMETER_STORE: ParameterStore = ParameterStore::from_sources(); }`. Keep the raw map around if required for backward compatibility (see below).
4. Update helper functions:
   - `get_global_param` reads from `PARAMETER_STORE.globals` via pattern matching on the key OR, preferably, is rewritten to be used only by legacy slow paths.
   - Introduce typed accessors: `pub fn global() -> &'static GlobalScalars`, `pub fn drug_params() -> &'static DrugParameters`, etc.

## Backward-compatibility strategy

To avoid touching every call site at once:

- Retain `PARAMETERS` for now, but mark it deprecated. Construct it from the same `ParameterStore` data so both views stay in sync during migration.
- For functions like `get_bacteria_param`, switch to:
  ```rust
  pub fn get_bacteria_param(name: &str, suffix: &str) -> Option<f64> {
      let idx = bacteria_index(name)?;
      match suffix {
          "acquisition_log_odds_baseline" => Some(store.bacteria.acquisition_log_odds_baseline[idx] as f64),
          // … cases …
          _ => PARAMETERS.get(&format!("{}_{}", name, suffix)).copied(),
      }
  }
  ```
  This lets us migrate callers incrementally while maintaining behavior for rarely-used keys.
- Once high-frequency code paths stop using string-based helpers, we can prune the fallback logic.

## Future-proofing and validation

- Add `debug_assert!` guards in `from_hashmap` constructors to ensure every required key exists (better failure mode than silent defaulting).
- Consider emitting a dump of missing/unused keys to help keep configuration files clean.
- Unit-test the mapping from strings to arrays by feeding a small synthetic parameter map and verifying a handful of retrieved values.

## Migration checklist

1. Implement `ParameterStore` with `from_hashmap` for each sub-struct.
2. Populate `lazy_static` with both the store and, temporarily, the legacy `HashMap` (sourced from the same builder).
3. Update `Simulation::new` to pull potency and MIC thresholds from `store.drug_bacteria`.
4. Update `rules::apply_rules` to use typed accessors. Start with per-person scalars and drug potency loops, then migrate bacteria acquisition, HGT, and toxicity blocks.
5. Remove or gate expensive `format!` usage; the only runtime string operations should be logging or truly dynamic keys.
6. Once usage is fully migrated, delete the legacy `HashMap` and simplify helper APIs.
