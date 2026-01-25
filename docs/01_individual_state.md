# Individual State Variables

> **Source File**: `src/simulation/population.rs` → `struct Individual`

The `Individual` struct represents a single person in the simulation. Each person has demographic attributes, health state, and per-bacteria/per-drug tracking arrays.

---

## Table of Contents
1. [Demographic Variables](#demographic-variables)
2. [Location & Hospitalization](#location--hospitalization)
3. [Infection State Arrays](#infection-state-arrays)
4. [Resistance State Arrays](#resistance-state-arrays)
5. [Drug Treatment Arrays](#drug-treatment-arrays)
6. [Microbiome Arrays](#microbiome-arrays)
7. [Clinical Outcome Variables](#clinical-outcome-variables)
8. [Diagnostic Testing Arrays](#diagnostic-testing-arrays)
9. [Treatment Tracking Arrays](#treatment-tracking-arrays)
10. [Update Processes](#update-processes)

---

## Demographic Variables

### `id: usize`
- **Purpose**: Unique identifier for this individual
- **Initialized**: Sequential integer at population creation
- **Updated**: Never changes

### `age: i32`
- **Purpose**: Age in days since simulation start (1930)
- **Range**: Negative (not yet born) to ~38000+ (very old)
- **Initialized**: From demographic distribution sampling
- **Updated**: `age += 1` at end of each time step
- **Special Values**: 
  - Negative = not yet born (will be born when age reaches 0)
  - `age = 0` = day of birth

### `sex_at_birth: String`
- **Purpose**: Biological sex for sex-specific disease risks
- **Values**: `"male"` or `"female"`
- **Initialized**: 50/50 random at creation
- **Updated**: Never changes

### `perceived_penicillin_allergy: bool`
- **Purpose**: Affects drug selection (avoids penicillin-class drugs)
- **Initialized**: Random based on `penicillin_allergy_probability` (~10%)
- **Updated**: Never changes (persistent attribute)
- **Effect**: If true, drug selection skips PENICILLIN_CLASS_DRUGS

---

## Location & Hospitalization

### `region_living: Region`
- **Purpose**: Permanent home region (affects drug availability, resistance patterns)
- **Values**: `NorthAmerica`, `SouthAmerica`, `Africa`, `Asia`, `Europe`, `Oceania`
- **Initialized**: From demographic distribution
- **Updated**: Never changes

### `region_cur_in: Region`
- **Purpose**: Current location (can differ from home during travel)
- **Initialized**: Same as `region_living`
- **Updated**: Travel logic (not currently active)

### `days_visiting: u32`
- **Purpose**: Days spent in current non-home region
- **Initialized**: 0
- **Updated**: Increment during travel, reset on return home

### `hospital_status: HospitalStatus`
- **Purpose**: Whether currently hospitalized
- **Values**: `InHospital`, `NotInHospital`
- **Initialized**: `NotInHospital`
- **Updated**: 
  - → `InHospital`: On sepsis, severe infection, or random admission
  - → `NotInHospital`: After `hospitalization_duration` days or discharge

### `days_hospitalized: u32`
- **Purpose**: Counter for current hospitalization duration
- **Initialized**: 0
- **Updated**: 
  - `+= 1` each day while hospitalized
  - Reset to 0 on discharge

---

## Infection State Arrays

All arrays are indexed by bacteria index (`b_idx`), matching `BACTERIA_LIST` order.

### `level: Vec<f64>` [bacteria]
- **Purpose**: Current infection intensity (bacterial load proxy)
- **Range**: 0.0 to `max_resistance_level` (typically 1.0)
- **Initialized**: 0.0 for all bacteria
- **Updated**:
  ```
  IF acquiring new infection:
      level[b] = initial_infection_level (typically 0.5)
  ELSE IF infected:
      level[b] -= immune_clearance_effect + drug_effect
      IF level[b] < INFECTION_EPS (0.001):
          level[b] = 0.0  // Cleared
  ```
- **Interpretation**: `level > INFECTION_EPS` means active infection

### `date_last_infected: Vec<i32>` [bacteria]
- **Purpose**: Time step when current infection started
- **Initialized**: -1 (never infected)
- **Updated**: Set to `time_step` when new infection acquired
- **Reset**: Set to -1 when infection clears
- **Usage**: Calculate infection duration, day-7 assessment

### `date_last_infected_keep: Vec<i32>` [bacteria]
- **Purpose**: Persistent record of last infection (not reset on clearance)
- **Initialized**: -1
- **Updated**: Set to `time_step` when new infection acquired
- **Never Reset**: Keeps historical record even after clearance

### `infectious_syndrome: Vec<i32>` [bacteria]
- **Purpose**: Clinical syndrome type (affects drug selection)
- **Values**: Syndrome IDs (0-10+, see syndrome mapping)
- **Initialized**: -1 (no syndrome)
- **Updated**: Set when infection acquired, based on bacteria + random
- **Examples**: UTI=0, Pneumonia=1, Skin/soft tissue=2, etc.

### `predicted_infection_risk: Vec<f64>` [bacteria]
- **Purpose**: Daily infection probability from logistic model
- **Range**: 0.0 to 1.0
- **Updated**: Calculated fresh each day based on:
  - Base acquisition rate for bacteria
  - Age-specific modifiers
  - Region factors
  - Hospital status
  - Immunodeficiency
  - Sanitation era adjustments

### `sepsis: Vec<bool>` [bacteria]
- **Purpose**: Life-threatening complication flag
- **Initialized**: false
- **Updated**:
  - `true`: When sepsis develops (random based on risk factors)
  - `false`: When infection clears or death occurs

### `sepsis_onset_day: Vec<i32>` [bacteria]
- **Purpose**: When sepsis started for mortality timing
- **Initialized**: -1
- **Updated**: Set to `time_step` when sepsis develops

### `cur_infection_from_environment: Vec<bool>` [bacteria]
- **Purpose**: Whether infection was environmentally acquired
- **Initialized**: false
- **Updated**: Set on new infection based on `environmental_acquisition_proportion`
- **Effect**: Affects resistance assignment source

### `infection_hospital_acquired: Vec<bool>` [bacteria]
- **Purpose**: Whether infection was hospital-acquired
- **Initialized**: false
- **Updated**: Set to `hospital_status.is_hospitalized()` on new infection
- **Effect**: Hospital-acquired infections sample from hospitalized population resistance

### `infection_has_caused_symptoms: Vec<bool>` [bacteria]
- **Purpose**: Clinical symptoms manifested (gates treatment)
- **Initialized**: false
- **Updated**: Set true when symptoms develop (probability-based)
- **Effect**: Treatment only initiated if symptoms present

### `clearance_hazard: Vec<f64>` [bacteria]
- **Purpose**: Daily immune clearance probability (for reporting)
- **Range**: 0.0 to 1.0
- **Updated**: Calculated each day based on immune function + drug effects

### `clearance_ready_day: Vec<i32>` [bacteria]
- **Purpose**: Day when hazard-based clearance can activate
- **Initialized**: -1 (not armed)
- **Updated**: Set when clearance conditions met
- **Usage**: Prevents immediate clearance, ensures minimum infection duration

### `infection_prevented_by_drug: Vec<bool>` [bacteria]
- **Purpose**: Tracks if prophylactic drugs prevented infection this step
- **Initialized**: false each timestep
- **Updated**: Set true if active drug prevented infection acquisition
- **Reset**: At start of each timestep

---

## Resistance State Arrays

### `resistances: Vec<Vec<Resistance>>` [bacteria][drug]
- **Purpose**: Per-bacteria, per-drug resistance levels
- **Structure**: 2D array indexed by `[bacteria_idx][drug_idx]`
- **Sub-fields in Resistance struct**:

#### `resistances[b][d].any_r: f64`
- **Purpose**: Resistance in any fraction of infecting bacteria
- **Range**: 0.0 (fully susceptible) to 1.0 (fully resistant)
- **Updated**:
  - Set on new infection from MajorityRCache sampling
  - Increased by de novo emergence under drug pressure
  - Can revert toward 0 when drugs removed (fitness cost)

#### `resistances[b][d].majority_r: f64`
- **Purpose**: Resistance in >50% of bacteria
- **Range**: 0.0 to any_r (never exceeds any_r)
- **Updated**: Set equal to any_r when resistance is majority
- **Usage**: Determines if resistance will persist after clearance

#### `resistances[b][d].microbiome_r: f64`
- **Purpose**: Resistance level in microbiome carriage
- **Range**: 0.0 to 1.0
- **Updated**: From microbiome sampling, HGT events
- **Effect**: Can transfer to infection resistance

#### `resistances[b][d].test_r: f64`
- **Purpose**: Resistance level as detected by testing
- **Updated**: Set when resistance test performed
- **May differ from any_r**: Testing sensitivity/specificity

#### `resistances[b][d].activity_r: f64`
- **Purpose**: Effective resistance for drug activity calculation
- **Typically**: Same as any_r

### `resistance_mechanisms: Vec<Vec<bool>>` [bacteria][mechanism]
- **Purpose**: Which specific resistance mechanisms are present
- **Indexed by**: `ResistanceMechanism::all()` order
- **Mechanisms**:
  - ESBL (extended-spectrum beta-lactamase)
  - Carbapenemase
  - AmpC (chromosomal beta-lactamase)
  - 16S Methyltransferase (aminoglycoside resistance)
  - Qnr (quinolone resistance)
  - Efflux Overexpression
  - Erm Methylation (macrolide resistance)
  - VanType (glycopeptide resistance)
  - MecA (methicillin resistance)
  - Reduced Permeability
  - Target Site Mutation

### `how_resistance_acquired: Vec<Vec<Option<ResistanceAcquisitionType>>>` [bacteria][drug]
- **Purpose**: Tracks source of resistance for each bacteria-drug pair
- **Values**:
  - `None`: Never had resistance
  - `AtInfectionCommunity`: Acquired resistant strain from community
  - `AtInfectionEnv`: Acquired from environmental source
  - `AtInfectionTB`: TB-specific acquisition
  - `Hgt`: Horizontal gene transfer
  - `FromMicrobiomeR`: Transferred from microbiome
  - `DeNovoInfection`: Emerged during treatment

---

## Drug Treatment Arrays

All arrays indexed by drug index (`d_idx`), matching `DRUG_SHORT_NAMES` order.

### `cur_use_drug: Vec<bool>` [drug]
- **Purpose**: Currently taking this drug
- **Initialized**: false
- **Updated**:
  - `true`: When drug initiated
  - `false`: When treatment completed or stopped

### `cur_level_drug: Vec<f64>` [drug]
- **Purpose**: Current drug concentration
- **Range**: 0.0 to 10.0+ (10.0 = standard daily dose)
- **Updated**:
  ```
  IF taking drug today:
      cur_level_drug[d] = 10.0 * dose_multiplier
  ELSE:
      cur_level_drug[d] *= exp(-ln(2) / half_life)  // Decay
  ```

### `date_drug_initiated: Vec<i32>` [drug]
- **Purpose**: When current treatment course started
- **Initialized**: -1
- **Updated**: Set to `time_step` when drug initiated
- **Reset**: Set to -1 when drug stopped

### `date_drug_initiated_keep: Vec<i32>` [drug]
- **Purpose**: Persistent record (not reset on drug stop)
- **Initialized**: -1
- **Updated**: Set to `time_step` when drug initiated
- **Never Reset**: Historical record

### `ever_taken_drug: Vec<bool>` [drug]
- **Purpose**: Lifetime drug exposure history
- **Initialized**: false
- **Updated**: Set true when drug first taken
- **Never Reset**: Permanent flag

### `drug_toxicity_reservoir: Vec<f64>` [drug]
- **Purpose**: Accumulated toxicity from each drug
- **Range**: 0.0 to unbounded
- **Updated**:
  - `+= drug_level * toxicity_accumulation_rate` each day
  - Decays with `toxicity_decay_rate`

### `drug_activity_response_multiplier: Vec<f64>` [bacteria]
- **Purpose**: Individual response variation to drugs
- **Range**: 0.5 to 1.5 (sampled from distribution)
- **Updated**: Sampled when treatment starts
- **Effect**: Multiplies drug activity for this person

---

## Microbiome Arrays

### `presence_microbiome: Vec<bool>` [bacteria]
- **Purpose**: Colonized with this bacteria (not infected)
- **Initialized**: Random based on baseline carriage rates
- **Updated**:
  - `true`: Acquisition from environment/contacts
  - `false`: Clearance (natural or antibiotic-induced)

### `date_microbiome_acquired: Vec<i32>` [bacteria]
- **Purpose**: When carriage was acquired
- **Initialized**: 0
- **Updated**: Set to `time_step` on acquisition

### `microbiome_acquired_today: Vec<bool>` [bacteria]
- **Purpose**: Flag for new acquisitions this timestep
- **Reset**: At start of each timestep
- **Usage**: Aggregation and reporting

### `microbiome_acquired_on_drug_today: Vec<bool>` [bacteria]
- **Purpose**: Whether acquisition occurred during antibiotic use
- **Reset**: At start of each timestep
- **Effect**: May affect resistance level of acquired strain

### `microbiome_cleared_today: Vec<bool>` [bacteria]
- **Purpose**: Flag for clearance events this timestep
- **Reset**: At start of each timestep

---

## Clinical Outcome Variables

### `current_infection_related_death_risk: f64`
- **Purpose**: Daily mortality risk from infections
- **Range**: 0.0 to 1.0
- **Updated**: Calculated based on sepsis status, treatment effectiveness

### `background_all_cause_mortality_rate: f64`
- **Purpose**: Age-adjusted background mortality
- **Updated**: From actuarial tables based on age

### `current_toxicity_hazard: f64`
- **Purpose**: Daily risk from drug toxicity
- **Updated**: Sum of all drug toxicity contributions

### `mortality_risk_current_toxicity: f64`
- **Purpose**: Death probability from current toxicity level
- **Updated**: From toxicity hazard via dose-response model

### `date_of_death: Option<usize>`
- **Purpose**: Time step when individual died
- **Values**: `None` (alive) or `Some(timestep)`
- **Updated**: Set once on death, never changes

### `cause_of_death: Option<String>`
- **Purpose**: Primary cause of death
- **Values**: `"sepsis"`, `"infection"`, `"toxicity"`, `"background"`, etc.

### `immunodeficiency_type: Option<ImmunodeficiencyType>`
- **Purpose**: Immunocompromised status
- **Values**:
  - `None`: Normal immune function
  - `Some(Temporary)`: Chemotherapy, transplant (recoverable)
  - `Some(Chronic)`: Primary immunodeficiency (permanent)
- **Effect**: Increases infection risk and severity

---

## Diagnostic Testing Arrays

### `test_identified_infection: Vec<bool>` [bacteria]
- **Purpose**: Infection confirmed by testing
- **Updated**: Set true when diagnostic test positive
- **Effect**: Enables targeted (rather than empiric) treatment

### `test_for_resistance: Vec<bool>` [bacteria]
- **Purpose**: Resistance testing performed
- **Updated**: Set true when susceptibility testing done

### `resistance_test_initiated_day: Vec<i32>` [bacteria]
- **Purpose**: When resistance testing started
- **Updated**: Set to `time_step` when test ordered
- **Usage**: Delay before results available

---

## Treatment Tracking Arrays

### `bacteria_level_at_drug_start: Vec<Option<f64>>` [bacteria]
- **Purpose**: Infection severity when treatment began
- **Usage**: Assess treatment response

### `days_on_current_treatment: Vec<i32>` [bacteria]
- **Purpose**: Duration of current treatment course
- **Updated**: Increment daily while on treatment
- **Reset**: When treatment changes or stops

### `treatment_failure_assessed: Vec<bool>` [bacteria]
- **Purpose**: Whether day-3 failure assessment done
- **Reset**: When new treatment starts

### `drug_stopped_with_infection_day: Vec<Option<i32>>` [bacteria]
- **Purpose**: When drug stopped while still infected
- **Usage**: Track incomplete treatment courses

### `bacteria_level_at_drug_cessation: Vec<Option<f64>>` [bacteria]
- **Purpose**: Infection level when treatment stopped

### `stopped_drug_index: Vec<Option<usize>>` [bacteria]
- **Purpose**: Which drug was stopped

### `restart_window_assessed: Vec<bool>` [bacteria]
- **Purpose**: Whether restart opportunity was evaluated

### `date_last_drug_failure: Vec<i32>` [bacteria]
- **Purpose**: When treatment last failed
- **Usage**: Affects subsequent drug selection

### `day_7_since_last_infection_drug_used: Vec<Option<bool>>` [bacteria]
- **Purpose**: Was any drug used in first 7 days of infection?
- **Updated**: Set on day 7 post-infection
- **Usage**: Treatment timing analysis

### `infection_resolution_this_timestep: Vec<Vec<u32>>` [bacteria][resolution_type]
- **Purpose**: Counts of resolution events this timestep
- **Resolution Types**: ImmuneClearance, DrugAssistedClearance, DeathFromSepsis, etc.
- **Reset**: At start of each timestep

### `cleared_any_r_microbiome_categories: Vec<[u32; 4]>` [bacteria]
- **Purpose**: Clearance counts by microbiome resistance context
- **Categories**: NoMicrobiome, MicrobiomePresentNoResistance, MicrobiomeMinorityR, MicrobiomeMajorityR
- **Reset**: At start of each timestep

### `bacteria_on_selection_day: i32`
- **Purpose**: Which bacteria triggered today's drug selection
- **Values**: Bacteria index or -1 (no selection)

### `drug_score_on_selection_day: Vec<f64>` [drug]
- **Purpose**: Scores computed during drug selection
- **Usage**: Debugging and analysis

### `current_number_of_drugs: i32`
- **Purpose**: Count of active drug treatments
- **Updated**: Recalculated when drugs start/stop

---

## Update Processes

### Per-Timestep Update Sequence (in `apply_rules`)

1. **Reset daily flags**
   ```rust
   infection_prevented_by_drug[*] = false
   microbiome_acquired_today[*] = false
   microbiome_cleared_today[*] = false
   infection_resolution_this_timestep[*][*] = 0
   ```

2. **Drug level decay** (for all drugs)
   ```rust
   if !cur_use_drug[d] {
       cur_level_drug[d] *= exp(-ln(2) / half_life[d])
   }
   ```

3. **Check birth** (if age < 0)
   ```rust
   if age == 0 {
       // Initialize newborn state
   }
   ```

4. **Vaccination updates** (age-appropriate)

5. **Hospital status updates**
   ```rust
   if hospital_status == InHospital {
       days_hospitalized += 1
       if days_hospitalized >= max_stay || can_discharge {
           hospital_status = NotInHospital
           days_hospitalized = 0
       }
   }
   ```

6. **Microbiome dynamics**
   - Acquisition probability
   - Clearance probability
   - Resistance transfer

7. **Infection acquisition** (for each bacteria)
   - Calculate predicted_infection_risk
   - Roll for new infection
   - If acquired: set level, syndrome, resistance

8. **Infection progression** (for each infected bacteria)
   - Symptom development
   - Sepsis development
   - Drug selection (if symptomatic)
   - Drug effects
   - Clearance check

9. **Resistance dynamics**
   - De novo emergence
   - HGT events
   - Reversion (fitness cost)

10. **Mortality check**
    - Sepsis mortality
    - Toxicity mortality
    - Background mortality

11. **Increment age**
    ```rust
    age += 1
    ```

---

## Memory Layout Notes

- Arrays sized to `BACTERIA_LIST.len()` (39 bacteria)
- Drug arrays sized to `DRUG_SHORT_NAMES.len()` (52 drugs)
- Total per-individual memory: ~10-15 KB
- Population of 100,000: ~1-1.5 GB

