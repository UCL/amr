# AMR Simulation — Model Description



## Contents

1. [Overview](#1-overview)
2. [Population and Demographics](#2-population-and-demographics)
3. [Infection Acquisition](#3-infection-acquisition)
4. [Clinical Progression](#4-clinical-progression)
5. [Diagnostic Testing](#5-diagnostic-testing)
6. [Antibiotic Treatment](#6-antibiotic-treatment)
7. [Resistance Dynamics](#7-resistance-dynamics)
8. [Microbiome and Carriage](#8-microbiome-and-carriage)
9. [Horizontal Gene Transfer](#9-horizontal-gene-transfer)
10. [Mortality](#10-mortality)
11. [Policy Evaluation](#11-policy-evaluation)
12. [Limitations](#12-limitations)
- [Appendix A — Bacteria, Drugs, Mechanisms and Enums](#appendix-a--bacteria-drugs-mechanisms-and-enums)
- [Appendix B — Parameter Reference](#appendix-b--parameter-reference)
- [Appendix C — Output Specification](#appendix-c--output-specification)

---



## 1. Overview

This model simulates the emergence and dynamics of antimicrobial resistance (AMR) across a synthetic human population from 1930 to 2035. It is an individual-based (agent-based) model in which each person can acquire bacterial infections, receive antibiotic treatment, develop resistance through de novo mutation or horizontal gene transfer, and carry resistant organisms in their microbiomes.

The model tracks **42 bacterial species**, **58 antibiotics** (grouped into **36 drug classes**), and **35 resistance mechanisms**. The population is distributed across **6 world regions** (North America, Europe, Asia, Oceania, South America, Africa), each with distinct epidemiological, travel, hospitalization, and healthcare profiles.

Time advances in discrete **daily steps**. On each day, every individual in the population is processed through a sequence of 21 mechanistic rules covering ageing, infection acquisition, clinical progression, treatment, resistance dynamics, and mortality.



### Scope and Purpose

The model is specifically designed for reconstructing the historical emergence and growth of AMR over time by mechanistically linking antibiotic consumption, biological mutability, and transmission. It evaluates the potential impact of antibiotic stewardship policies by recreating empirical observations of resistance incidence and separating resistance acquisition across different care settings (e.g., community-acquired versus hospital-acquired).

---



## 2. Population and Demographics



### 2.1 Initialisation

The population is created at day 0 (representing the calendar year 1930). Each individual is assigned:

- **Age**: Drawn from a continuous demographic distribution that encodes both living individuals and future births. Negative age values at initialisation represent individuals who have not yet been born, entering the simulation exactly when their age reaches zero.
- **Sex**: Male or female, assigned with equal probability.
- **Region**: Sampled from demographic weights reflecting the global population distribution.

The six regions and their approximate population shares determine the starting geographical distribution:

| Region | Population Share |
|--------|------------------|
| Asia | ~55% |
| Europe | ~15% |
| Africa | ~12% |
| North America | ~9% |
| South America | ~6% |
| Oceania | ~3% |



### 2.2 Ageing and Categorisation

Each day, every individual's age increments by one day. Age categories are continually recalculated and dictate risk profiles for infection, testing, and hospitalisation:

| Age Category | Age Range | Variable Suffix |
|--------------|-----------|-----------------|
| Infant | 0–1 year | `infant` |
| Preschool | 1–5 years | `preschool` |
| School Age | 5–18 years | `school` |
| Young Adult | 18–50 years | `young_adult` |
| Middle Age | 50–70 years | `middle_age` |
| Elderly | 70+ years | `elderly` |

An overlapping, more granular classification is used specifically for calculating sepsis onset and infection-associated mortality. This temporary grouping isolates **Neonates (0-28 days)**, who possess vastly different microbiological vulnerabilities, maternal antibody profiles, and case-fatality risks compared to older infants. This prevents neonatal vulnerability from being artificially smoothed over the broader 0-1 year demographic bracket:

| Sepsis/Mortality Category | Age Range | Variable Suffix |
|---------------------------|-----------|-----------------|
| Neonatal | 0–28 days | `neonatal` |
| Pediatric | 28 days–18 years | `pediatric` |
| Young Adult | 18–50 years | `young_adult` |
| Elderly | 50+ years | `elderly` |



### 2.3 Immunodeficiency

Individuals can be immunosuppressed, greatly altering their risk profile across infection acquisition, clinical progression, sepsis, and death. Immunosuppression is divided into temporary and chronic states.

**State Transitions:**
- **Temporary Immunosuppression:** Occurs at a rate of `0.00005` per day (calibrated for ~5% prevalence) and resolves at a rate of `0.01` per day (representing acute illnesses or short-acting treatments).
- **Chronic Immunosuppression:** Occurs at a lower rate of `0.00006` per day but recovers much slower at `0.0012` per day (representing long-term conditions like HIV, autoimmune diseases, or organ transplants).

At simulation initialisation (and for new births), there is a baseline probability of existing chronic immunodeficiency, scaled heavily by age to represent cumulative lifetime risk:
- Infants (0-1): `0.3` (Higher chance representing congenital/genetic immunodeficiencies)
- Children (1-18): `0.2` (Moderate chance)
- Adults (18-65): `0.4` (High chance representing acquired conditions like HIV or immunosuppressive therapies)
- Elderly (65+): `0.6` (Highest chance representing immunosenescence and multiple compounding morbidities)

**Clinical Impacts of Immunodeficiency:**
Immunosuppressed individuals experience multipliers across their entire clinical journey:
- **Antibiotic Initiation:** Moderately more likely to receive empirical antibiotics even without a confirmed indication (`antibiotic_initiation_log_odds_immunodeficiency` = `2.08`, equivalent to ~8.0x odds multiplier; penalty for unindicated prescribing is only `-1.05` instead of higher penalties for immunocompetent individuals).
- **Diagnostic Testing:** Receive 2.5x more diagnostic testing (`testing_immunosuppressed_multiplier` = `2.5`).
- **Sepsis Onset:** ~2x higher onset risk (`log_odds_sepsis_onset_immunosuppressed` = `0.7`).
- **Recovery:** Significantly lower probability of spontaneous recovery from sepsis (`sepsis_log_odds_immunosuppressed` = `-1.0`).
- **Mortality:** Massive compound risks facing death. +1.5 log-odds for sepsis-related death (~4.5x risk), +0.9 log-odds for drug toxicity death (~2.5x risk), and an overall baseline mortality multiplier of ~2.5x (ln(2.5) = `0.916`).



### 2.4 Hospitalisation

Hospital admission is modelled as a daily logistic function evaluating the individual's clinical state, age, and regional healthcare access.

**Probability of Admission:**  
P(admission) = 1 / (1 + e^(-log_odds))

Where `log_odds` is fundamentally driven by:
`log_odds = base + age_effect + sepsis_effect + infection_effect + region_effect`

| Factor | Value | Description |
|--------|-------|-------------|
| **Baseline** | `-10.4` | `hospitalization_base_log_odds` (~0.003% daily admission risk for healthy individuals). |
| **Age** | `+0.02 / yr` | `hospitalization_log_odds_per_age_year` (~2% increase in log-odds per year of age). |
| **Sepsis** | `+4.4` | `hospitalization_log_odds_sepsis` (Massive driver, roughly ~80x admission multiplier). |
| **Infection** | `+2.5` | `hospitalization_log_odds_symptomatic_infection` (~12.2x multiplier for symptomatic cases surpassing the severity threshold of `3.0`). |

Once admitted, patients have an average length of stay of ~3.6 days (`hospitalization_recovery_rate_per_day` = `0.28`), capped at a maximum of 30 days (`hospitalization_max_days` = `30.0`). Discharge is completely blocked if the patient is currently septic (`hospitalization_prevent_discharge_with_sepsis` = `1.0`).

**Regional Capacity Modifiers:**
Healthcare access heavily modulates admission thresholds:
- North America: `0.5` (Good access)
- Europe: `0.6` (Universal healthcare, highest access)
- Oceania: `0.4` (Good access in developed areas)
- Asia: `0.0` (Reference baseline, mixed access)
- South America: `-0.2` (Variable access)
- Africa: `-0.5` (Limited hospital capacity and access)

**Nosocomial (Hospital-Acquired) Risks:**
Being in the hospital dramatically alters infection acquisition risks and clinical outcomes. 
- Mortality and Sepsis: In-hospital patients face higher baseline mortality (`+0.262` log-odds, ~1.3x), higher sepsis onset risk due to sickness acuity (`+0.5` log-odds, ~1.6x), but also counter-intuitively a higher probability of *recovering* from sepsis (`+0.8` log-odds) due to intensive care interventions.
- Pathogen Acquisition: The hospital environment amplifies the risk of specific nosocomial pathogens:
  - *A. baumannii*: `+3.4` log-odds (~30x higher risk)
  - *P. aeruginosa*: `+3.0` log-odds (~20x higher risk)
  - *E. faecium* (VRE): `+3.3` log-odds (~27x higher risk)
  - *S. aureus* (MRSA): `+2.3` log-odds (~10x higher risk)
  - *K. pneumoniae* / *Serratia spp.*: `+2.0` log-odds (~7x higher risk)
  - Extracellular/Community pathgens (*M. genitalium*, *C. trachomatis*, *T. pallidum*, *Campylobacter*): Negative modifiers (`-0.6` to `-1.5`), meaning they are significantly *less* likely to be acquired inside the hospital than in the community.



### 2.5 Travel

Individuals may travel between regions, mixing diverse drug availability phenotypes and importing resistance profiles. The baseline probability of travelling on any given day is `0.00005` (`travel_probability_per_day`).

Travel frequency is strongly tied to the income level and development of the originating region:
- Europe: `3.5x` multiplier
- North America: `3.0x` multiplier
- Oceania: `2.5x` multiplier
- Asia: `1.5x` multiplier
- South America: `0.8x` multiplier
- Africa: `0.3x` multiplier

Age and travel patterns also interact to determine specific pathogen risks. For example, young adults and middle-aged individuals travelling from Europe and North America face specific, elevated log-odds scalar modifiers for acquiring travel-related enteric diseases such as *Salmonella enterica* serovar Typhi (up to `+0.8` log-odds for young European adults) and *Shigella spp.* (up to `+0.5` log-odds). Conversely, risks for *Vibrio cholerae* are suppressed (`-1.0`) for these same travelling demographics unless visiting highly endemic zones.

---



## 3. Infection Acquisition



### 3.1 Community Acquisition

Each day, each non-infected individual has a probability of acquiring each of the 42 bacterial species. The probability is computed via a logistic model combining:

- Base acquisition rate for the bacteria
- Regional modifier
- Age-dependent risk (via templates and overrides)
- Immunodeficiency modifier
- Seasonal variation (sinusoidal for respiratory pathogens)
- Calendar-era effects (temporal multiplier)
- Population-level resistance prevalence (`majority_r`) reducing or altering the probability that the bacteria the individual acquires is definitively resistant

| Variable pattern | Description |
|------------------|-------------|
| `bacteria_{name}_acquisition_log_odds` | Baseline acquisition log-odds for each species |
| `{region}_bacteria_{name}_acquisition_log_odds` | Regional override |
| `bacteria_{name}_log_odds_{age_category}` | Age-specific override |
| `{bacteria}_{region}_log_odds_{age_category}` | Bacteria × region × age interaction |



#### Age Risk Templates

Each bacteria is assigned a base risk template that defines its default age distribution profile before specific overrides are applied:

| Template | Multipliers [0–1y, 1–5y, 5–18y, 18–50y, 50–70y, 70+y] |
|----------|--------------------------------------------------------|
| `respiratory` | [3.0, 1.8, 0.8, 1.0, 1.3, 2.5] |
| `gastrointestinal` | [2.5, 2.0, 1.2, 1.0, 1.1, 1.8] |
| `urogenital` | [1.2, 0.8, 0.9, 1.0, 1.4, 2.2] |
| `skin_soft_tissue` | [1.5, 1.3, 1.1, 1.0, 1.2, 1.8] |
| `bloodstream` | [4.0, 2.0, 0.7, 1.0, 1.5, 3.0] |
| `sexually_transmitted` | [0.1, 0.2, 0.8, 1.0, 0.8, 0.3] |
| `flat` | [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] |

**Template Assignments:** Most bacteria default to `"respiratory"` but are specifically overridden. Ex: `salm_typhi_age_risk_template` (`"gastrointestinal"`), `esch_coli_age_risk_template` (`"urogenital"`), `pseud_aerug_age_risk_template` (`"bloodstream"`), `n_gonorrhoeae_age_risk_template` (`"sexually_transmitted"`).



### 3.2 Hospital Acquisition

Hospitalised individuals bypass standard community acquisition rules and are subject to nosocomial acquisition log-odds rates (`{bacteria}_log_odds_hospital_acquired`), scaling up the likelihood of picking up hospital-associated pathogens (e.g. *A. baumannii*, *P. aeruginosa*, MRSA) while often suppressing typical community pathogens according to explicit parameters documented in Section 2.4.



### 3.3 Carrier-Derived Infection

Individuals carrying bacteria asymptomatically in their microbiomes can develop an active infection driven from their own carriage flora. This endogenous pathway is critical for modeling resistance because it enables a previous resistant strain to cause a future clinical infection.
This pathway is governed by:

| Parameter | Baseline Value | Description |
|----------|---------|-------------|
| `carrier_resistance_inheritance_probability` | 0.50 | Probability the newly triggered infection identically inherits the specific resistance payload stored in the individual's microbiome. |
| `infection_from_microbiome_dampening` | 0.70 | A dampening factor ensuring carriage doesn't overwhelmingly explode into clinical infection cases beyond calibrated biological rates. |



### 3.4 Resistance at Acquisition

When a new infection is acquired from the community (exogenous transmission), its initial resistance profile is generated via a robust transmission mechanism reflecting the macro-state of the population:

1. **Population-level prevalence** (`majority_r`): The model continuously captures a dynamic rolling average of definitive resistance across all infected individuals for every bacteria–drug class combination. The time window is governed by `majority_r_window_days` (defaulting to 1000 days or ~3 years of latency bias).
2. **Community resistance dilution**: A scaling factor (`community_resistance_dilution_factor`, typically 0.50) attenuates strictly localized/clinical resistance signals so that environmental and community reservoirs appropriately dilute the high resistance seen directly in clinical settings.
3. **Correlated Mechanism Profiles**: Instead of assigning individual resistance genes randomly (which would biologically misrepresent plasmid and structural linkage), the engine accesses a `MechanismProfileCache`. It samples completely intact matrices of resistance mechanisms directly from other concurrent individuals, guaranteeing the propagation of realistic multi-drug resistance phenotypes (MDR variants).

---



## 4. Clinical Progression



### 4.1 Syndrome Assignment

Upon progressing to an active infection, an anatomical clinical syndrome is determined. This dramatically controls how the bacteria will organically reproduce, what treatments providers will logically select, and how easily drugs will physically penetrate the tissue (`syndrome_{id}_{drug}_penetration`).

| Syndrome | Index | Biological/Clinical Focus |
|----------|-------|--------------------------|
| UTI | 1 | Urinary tract infection; high penetrance for common oral drugs |
| Skin/soft tissue | 2 | Cellulitis, wounds, burns |
| Respiratory | 3 | Community- and hospital-acquired pneumonia |
| Bloodstream | 4 | Bacteraemia |
| Intra-abdominal | 5 | Peritonitis, post-operative abdomen |
| CNS | 6 | Meningitis; restricted blood-brain barrier drug penetrance |
| Gastrointestinal | 7 | Gastroenteritis |
| Genital/pelvic | 8 | STIs, PID |
| Bone/joint | 9 | Osteomyelitis |
| Other | 10 | Device-related, undifferentiated |



#### Syndrome Progression Multipliers

Each syndrome modulates the probability of initiating empiric treatments, as well as the inherent daily reproductive rate of the bacteria. The simulation departs from flat logic by explicitly speeding up or slowing down pathophysiology depending on the syndrome site:

| Syndrome ID | Name | Initiation Multiplier (`_initiation_multiplier`) | Growth Rate Multiplier (`_bacteria_growth_multiplier`) |
|-------------|------|--------------------------------------------------|-------------------------------------------------------|
| 1 | UTI | 1.0 (Ref) | 1.0 (Baseline) |
| 2 | Skin | 1.0 | 1.1 (Faster in necrotic/cellulitis cases) |
| 3 | Respiratory | 10.0 (Huge boost: immediate treatment-seeking) | 1.2 (Rapid alveolar replication) |
| 4 | Bloodstream | 1.0 | 1.4 (Fulminant reproduction) |
| 5 | Intra-abdominal | 1.0 | 1.15 |
| 6 | CNS | 1.0 | 1.3 |
| 7 | GI | 8.0 | 1.1 |
| 8 | Genital | 12.0 (High distinct care-seeking logic) | 0.9 (Chronic/indolent) |
| 9 | Bone/joint | 1.0 | 0.85 (Slow, buried progression) |



### 4.2 Infection Dynamics

Rather than immediately classifying a host as "sick," the model tracks an internal numerical state of bacterial biomass per infection:
- **`initial_infection_level`**: Generally starts at `0.01` upon contraction.
- **Symptom Threshold**: The pathogen payload actively replicates daily until it reaches the `symptomatic_infection_level_threshold` (typically `3.0`), at which point the person exhibits definable symptoms that compel healthcare interaction. Below this threshold, they act identically to asymptomatic carriers for medical intervention purposes, though they still transmit.



### 4.3 Sepsis

Sepsis denotes a state of profound, life-threatening immunological failure provoked by the infection. Sepsis significantly elevates the risk of admission and mortality. Onset is tested daily as a logistic function summing multiple intersecting risk factors:

$$\text{log\_odds}_{\text{sepsis}} = \text{base}_{\text{bacteria}} + f(\text{bacterial\_level}) + f(\text{duration}) + \text{age} + \text{region} + \text{syndrome} + \text{modifiers}$$

| Sepsis Dynamic Multiplier | Default Val | Description |
|---------------------------|-------------|-------------|
| `log_odds_sepsis_infection_level` | `0.9` | Extremely potent scalable driver; log-odds added *per unit* of bacterial mass over threshold. |
| `log_odds_sepsis_infection_duration`| `0.005` | Tiny compound daily accumulation representing prolonged physiological exhaustion. |
| `sepsis_age_log_odds_neonatal` | `1.10` | ~3x elevated baseline risk factor for newborns. |
| `sepsis_age_log_odds_elderly` | `0.69` | ~2x elevated baseline risk factor for those >70 yrs. |



#### Base Sepsis Log-Odds by Bacteria & Syndrome
Not all bacteria or sites induce systemic failure equally:
- **Baseline Intercept**: A fallback of `-14.0` (virtually 0%) is assigned to rare organisms.
- **Aggressive/Common Organisms**: *E. coli* & *N. gonorrhoeae* are artificially pinned incredibly low (`-21.0`), preventing routine UTIs and STIs from frequently scaling into sepsis. *C. trachomatis* (`-19.0`), *Shigella* (`-12.0`), and *S. pneumoniae* (`-9.0`) reflect empirical likelihoods.
- **Syndrome Modifiers**: Superimposed on the pathogen is the anatomical site: 
  - Bloodstream (`+1.5`), CNS (`+1.2`), and Intra-abdominal (`+0.8`) sharply raise sepsis risk. 
  - Generitourinary (`-2.0`), Genital (`-1.5`) and Skin (`-1.0`) sharply lower it.



#### Regional Mitigation
Healthcare and sanitary environments suppress how fast untreated infections devolve into sepsis organically:
- North America & Oceania: `-0.5`
- Europe: `-0.6`
- Asia: `-0.1` (Baseline)
- Africa: `+0.1` (Poorer infrastructural buffers lead to marginally faster biological deterioration)



### 4.4 Natural Clearance

Infections resolve organically via macrophage/lymphatic destruction when not sufficiently deadly, calculated via daily probabilism:
- `default_microbiome_clearance_probability_per_day` = `0.01` (Roughly 1% per day chance of biological success empty-handed)
- Prolonged establishment solidifies the bug (`carriage_duration_log_odds_coefficient` = `-0.01` with a max `-2.0` log-odds reduction, rendering chronic cases over 7x fundamentally harder to wipe out).
- Active, effective antibiotic therapy drives the `microbiome_clearance_probability_on_drug_treatment` up to `0.80`, yielding extremely rapid bacterial destruction if the correct drug profile matches the resistance phenotype.

---



## 5. Diagnostic Testing



### 5.1 Historical Introduction

Bacterial identification via formal biological culture and specialized antimicrobial susceptibility testing (AST) become available mechanically at specific historical landmarks during the simulation run:

| Technology | Time Step variable | Value | ~ Calendar Year |
|------------|------------------|-------|---------------|
| **Bacterial Culture** | `bacterial_testing_available_from_day` | 5,478 | 1945 (15 yrs post-1930) |
| **AST (Resistance Testing)** | `resistance_testing_available_from_day` | 9,131 | 1955 (25 yrs post-1930) |

*(Note: Parameters governing gradual S-curve historical adoptions over time like `bacterial_testing_adoption_rate_per_year` are deprecated but structurally mapped for backward compatibility. Usage caps at `1.0` via `bacterial_testing_max_temporal_multiplier`.)*



### 5.2 Mechanics and Errors

Once testing occurs, the engine simulates structural lab delays and analytical errors:
- **Delay**: It takes `3.0` days to receive actionable diagnostic feedback (`test_delay_days`).
- **AST Execution**: If a culture is secured, there is a `0.95` probability (`prob_test_r_done`) that antimicrobial resistance testing is physically conducted on the sample.
- **Reporting Error**: When susceptible/resistant results are finalized, they suffer an intrinsic `0.02` (2%) probability of false/erroneous reading (`test_r_error_probability`), producing an artificially skewed error-value coefficient modeled as `0.25` (`test_r_error_value`).



### 5.3 Clinical and Demographical Drivers

The base probability of initiating diagnostic testing per day is governed dynamically in the modern era based on clinical state, setting, and demographics:

| Variable | Default Rate / Multiplier | Description |
|----------|---------|-------------|
| `bacterial_testing_base_rate_per_day` | `0.15` | Default base probability (15% per day) for symptomatic cases to receive culture orders. |
| `resistance_testing_base_rate_per_day`| `0.95` | Secondary rate representing almost guaranteed AST reflex when culture is executed. |
| `testing_sepsis_multiplier` | `4.0x` | Patients in clinical sepsis receive highly accelerated urgent testing. |
| `testing_immunosuppressed_multiplier` | `2.5x` | Immunosuppressed patients are monitored much more aggressively. |
| `bacterial_testing_hospital_multiplier` | `8.0x` | Being hospitalized drastically increases access to formal pathology labs (cultures). |
| `resistance_testing_hospital_multiplier` | `5.0x` | Being hospitalized elevates reflexive AST profiling. |



#### Regional Testing Adjustments
Healthcare capacity controls laboratory bandwidth regionally (`{region}_testing_multiplier`):
- **North America**: `1.1x` (High infrastructure)
- **Europe**: `1.2x` (Highest testing density)
- **Oceania**: `0.8x`
- **Asia**: `0.7x`
- **South America**: `0.6x`
- **Africa**: `0.3x` (Limited lab infrastructure)

---



## 6. Antibiotic Treatment



### 6.1 Treatment Initiation

The model executes the decision to begin prescribing antibiotics as evaluating a daily logistic function:

$$P(\text{initiation}) = \frac{1}{1 + e^{-\text{log\_odds}}}$$

Initiation log-odds accumulate through compound modifiers reflecting a physician's real-world diagnostic logic (and errors):

| Clinical Setup Variable | Default | Net Mechanism |
|----------|---------|-------------|
| `antibiotic_initiation_base_log_odds` | `-5.5` | Baseline likelihood of accidental/non-symptomatic prescribing (~0.4% baseline risk). |
| `antibiotic_initiation_log_odds_symptomatic_infection` | `+6.0` | Very powerful driver escalating probability from ~0.4% to roughly ~62% purely off symptoms breaching threshold. |
| `antibiotic_initiation_log_odds_sepsis` | `+6.0` | Further powerful boost representing emergent, life-saving care. |
| `antibiotic_initiation_log_odds_immunodeficiency` | `+2.08` | (ln 8) Heavy prophylactic bias to initiate regimens even tentatively due to patient vulnerability. |
| `antibiotic_initiation_log_odds_no_indication` | `-1.05` | (ln 0.35) Protective penalty. dampens probability of inappropriate prescribing if patient doesn't functionally have active infection. |
| `antibiotic_initiation_log_odds_test_identified` | `+0.92` | (ln 2.5) Lab confirmations prompt targeted initiation. |
| `antibiotic_initiation_log_odds_already_on_drug` | `+0.18` | (ln 1.2) Modest boost representing combined/layered therapy logic once a patient is already functionally linked to a pharmacy loop. |



#### Care Access by Region
Prescription habits and accessibility vary radically globally, driving profound geographic effects on usage:
- **North America, Europe, Oceania**: `0.0` (Reference models, ~26-62% symptomatic probability).
- **Asia**: `-0.5` (~38% reduction in odds).
- **South America**: `-0.8` (~55% reduction).
- **Africa**: `-1.4` (~75% reduction in probabilistic prescription rates due to access blocks).



### 6.2 Drug Selection and Cessation

Drug selection follows a two-stage scoring system mathematically evaluating clinical logic vs biologically effective combinations:

1. **Empiric therapy** (no test result available): Drugs are scored primarily based on a syndrome-specific historical template indicating usual standards of care.
2. **Targeted therapy** (test result available): Drugs are rigidly scored based on strictly known mechanical susceptibility profiling in the lab, minimizing guesswork.



#### Therapy Scoring Parameters
The simulation intensely penalizes inappropriate choices to mirror clinical governance, and rewards de-escalation:

| Phase | Variable | Value | Description |
|-------|----------|-------|-------------|
| **Empirical** | `empiric_therapy_broad_spectrum_bonus` | `0.85` | In empirical spaces lacking ID, broad-spectrum choices are favored conceptually but heavily penalized implicitly to curb overuse (0.85 multiplier restricts abuse). |
| **Empirical** | `empiric_therapy_ineffective_drug_penalty`| `0.001` | Massive penalty reducing the score if the drug happens to be biologically totally misaligned. |
| **Targeted** | `targeted_therapy_narrow_spectrum_bonus` | `5.0` | Strong positive multiplier to actively reward de-escalating to narrow-spectrum agents (like Penicillin G) once a specific bug is caught. |
| **Targeted** | `targeted_therapy_broad_spectrum_penalty` | `0.1` | Heavy negative weight forcing physicians off "cannon" broad-spectrum drugs when a narrow-spectrum alternative is viable. |
| **Targeted** | `targeted_therapy_ineffective_drug_penalty`| `0.001` | Absolute block against using definitively resistant combinations. |



#### Regional Resistance Surveillance Penalties
If population-level algorithmic resistance (`majority_r`) flags severe failure rates regionally, physicians empirically avoid that drug class entirely:

| Resistance Rate | Penalty Variable | Multiplier |
|-----------------|------------------|------------|
| > 60% Resistant | `regional_resistance_penalty_very_high` | `0.3` |
| > 45% Resistant | `regional_resistance_penalty_high` | `0.5` |
| > 10% Resistant | `regional_resistance_penalty_moderate` | `0.8` |



#### Treatment Cessation
Patients discontinue their antibiotic course probabilistically, reflecting variable course lengths per indication:
- **Baseline Default**: `random_drug_cessation_probability` = `0.0045` (~0.45% daily, representing ~94% compliance with completing a typical 14-day robust course).
- **False Alarm Halt**: `random_drug_cessation_probability_if_no_active_infection` = `0.15` (15% daily stop chance if testing/observation reveals the patient actually has zero active bacterial load).
- **Pathogen Specific Courses** (selected examples matching real guidelines):
  - *Cholera / E. coli*: Fast resolution courses (`0.025` probability; typical 3-5 days).
  - *S. aureus / S. pneumoniae*: Mid-range (`0.015`, ~7 days).
  - *MDR-TB*: Extremely chronic (`0.0006`, strictly modelling an arduous 6-24 month regimen).



#### Syndrome-Specific Empiric Scoring Templates

Each syndrome has a pre-defined scoring table ranking drugs by appropriateness. A higher score means that drug is more likely to be selected empirically for that syndrome.

**Syndrome 1 — UTI**

| Drug | Score |
|------|-------|
| nitrofurantoin | 15.0 |
| trim_sulf | 14.0 |
| ciprofloxacin | 12.0 |
| levofloxacin | 11.0 |
| cephalexin | 10.0 |
| amoxicillin_clavulanate | 9.0 |
| cefuroxime | 8.0 |
| ampicillin | 6.0 |
| amoxicillin | 6.0 |
| ceftriaxone | 5.0 |
| gentamicin | 4.0 |
| meropenem | 3.0 |

**Syndrome 2 — Skin/Soft Tissue**

| Drug | Score |
|------|-------|
| penicillin_g | 13.0 |
| amoxicillin_clavulanate | 12.0 |
| cephalexin | 12.0 |
| clindamycin | 11.0 |
| doxycycline | 10.0 |
| trim_sulf | 9.0 |
| cefazolin | 8.0 |
| vancomycin | 7.0 |
| linezolid | 6.0 |
| ciprofloxacin | 5.0 |
| azithromycin | 5.0 |
| metronidazole | 4.0 |

**Syndrome 3 — Respiratory**

| Drug | Score |
|------|-------|
| amoxicillin_clavulanate | 20.0 |
| amoxicillin | 12.0 |
| azithromycin | 12.0 |
| doxycycline | 11.0 |
| levofloxacin | 10.0 |
| ceftriaxone | 9.0 |
| cefuroxime | 8.0 |
| trim_sulf | 7.0 |
| moxifloxacin | 6.0 |
| penicillin_g | 6.0 |
| clarithromycin | 5.0 |
| ampicillin | 4.0 |
| ciprofloxacin | 3.0 |
| erythromycin | 3.0 |

**Syndrome 4 — Bloodstream**

| Drug | Score |
|------|-------|
| piperacillin_tazobactam | 18.0 |
| meropenem | 14.0 |
| vancomycin | 13.0 |
| ceftriaxone | 12.0 |
| cefepime | 11.0 |
| ampicillin_sulbactam | 10.0 |
| ciprofloxacin | 9.0 |
| levofloxacin | 9.0 |
| amoxicillin_clavulanate | 8.0 |
| linezolid | 7.0 |
| metronidazole | 5.0 |
| trim_sulf | 5.0 |
| azithromycin | 4.0 |
| gentamicin | 1.0 |
| tobramycin | 1.0 |
| amikacin | 1.0 |

**Syndrome 5 — Intra-abdominal**

| Drug | Score |
|------|-------|
| piperacillin_tazobactam | 14.0 |
| meropenem | 13.0 |
| ampicillin_sulbactam | 11.0 |
| metronidazole | 11.0 |
| ceftriaxone | 10.0 |
| ciprofloxacin | 9.0 |
| levofloxacin | 9.0 |
| cefepime | 8.0 |
| amoxicillin_clavulanate | 7.0 |
| clindamycin | 6.0 |
| vancomycin | 4.0 |
| azithromycin | 3.0 |

**Syndrome 6 — CNS**

| Drug | Score |
|------|-------|
| ceftriaxone | 15.0 |
| vancomycin | 14.0 |
| ampicillin | 12.0 |
| meropenem | 11.0 |
| penicillin_g | 10.0 |
| linezolid | 8.0 |
| chloramphenicol | 7.0 |
| rifampicin | 6.0 |
| ciprofloxacin | 5.0 |
| trim_sulf | 5.0 |
| metronidazole | 4.0 |
| doxycycline | 3.0 |

**Syndrome 7 — Gastrointestinal**

| Drug | Score |
|------|-------|
| ciprofloxacin | 12.0 |
| furazolidone | 11.0 |
| azithromycin | 11.0 |
| metronidazole | 10.0 |
| doxycycline | 9.0 |
| trim_sulf | 8.0 |
| ceftriaxone | 7.0 |
| amoxicillin_clavulanate | 7.0 |

**Syndrome 8 — Genital/Pelvic**

| Drug | Score |
|------|-------|
| azithromycin | 13.0 |
| ceftriaxone | 13.0 |
| doxycycline | 12.0 |
| penicillin_g | 12.0 |
| amoxicillin_clavulanate | 9.5 |
| amoxicillin | 9.0 |
| cefuroxime | 9.0 |
| clindamycin | 9.0 |
| ciprofloxacin | 7.0 |
| levofloxacin | 6.5 |
| trim_sulf | 5.0 |
| rifampicin | 4.0 |
| metronidazole | 2.5 |

**Syndrome 9 — Bone/Joint**

| Drug | Score |
|------|-------|
| penicillin_g | 14.0 |
| cefazolin | 13.0 |
| ampicillin | 12.0 |
| vancomycin | 12.0 |
| linezolid | 11.0 |
| cephalexin | 11.0 |
| ceftriaxone | 11.0 |
| tedizolid | 10.0 |
| dalbavancin | 10.0 |
| clindamycin | 10.0 |
| ciprofloxacin | 9.0 |
| levofloxacin | 9.0 |
| rifampicin | 9.0 |
| trim_sulf | 8.0 |
| meropenem | 7.0 |
| piperacillin_tazobactam | 6.5 |

**Syndrome 10 — Other/Device-Related**

| Drug | Score |
|------|-------|
| piperacillin_tazobactam | 8.0 |
| cefepime | 8.0 |
| ceftriaxone | 8.0 |
| meropenem | 8.0 |
| imipenem_c | 8.0 |
| vancomycin | 8.0 |
| linezolid | 7.0 |
| ciprofloxacin | 7.0 |
| azithromycin | 6.0 |



### 6.3 Drug Pharmacokinetics

Each drug has a decay half-life determining how quickly its level falls after administration:

| Parameter pattern | Baseline Value | Description |
|------------------|---------|-------------|
| `drug_{name}_half_life_days` | Drug-specific | PK half-life in days |
| `drug_{name}_initial_level` | 10.0 | Level at administration |
| `drug_{name}_double_dose_multiplier` | 2.0 | Level multiplier for double dose |
| `drug_{name}_spectrum_breadth` | 3.0 | Microbiome disruption potential (higher = broader) |



#### Drug half-lives (selected)

| Drug | Half-life (days) |
|------|-----------------|
| sulfanilamide | 0.29 |
| penicillin_g | 0.042 |
| ampicillin | 0.063 |
| amoxicillin | 0.063 |
| ceftriaxone | 0.33 |
| azithromycin | 2.92 |
| doxycycline | 0.75 |
| vancomycin | 0.25 |
| meropenem | 0.042 |
| ciprofloxacin | 0.17 |
| linezolid | 0.21 |
| colistin | 0.21 |
| dalbavancin | 14.0 |
| cefiderocol | 0.10 |



#### Spectrum breadth overrides

| Drug | Breadth | Classification |
|------|---------|---------------|
| penicillin_g | 2.0 | Narrow |
| vancomycin | 2.5 | Narrow-medium |
| linezolid | 2.0 | Narrow |
| trim_sulf | 3.5 | Medium-broad |
| azithromycin | 4.0 | Broad |
| colistin | 4.0 | Broad |
| ciprofloxacin | 4.5 | Very broad |
| ceftriaxone | 4.0 | Broad |
| cefepime | 4.0 | Broad |
| meropenem | 5.0 | Very broad |



### 6.4 Drug Penetration by Syndrome

Drug efficacy at different infection sites varies by drug class. Penetration values range from 0.0 (no effective concentration) to 1.0 (full systemic levels):

| Syndrome | Best penetration | Worst penetration |
|----------|-----------------|-------------------|
| UTI (1) | FQ, TMP-SMX, nitrofurantoin, fosfomycin (1.0) | Macrolides (0.4), clindamycin (0.3), daptomycin (0.1) |
| Skin (2) | Daptomycin (0.95), FQ (0.9), oxazolidinones (0.9) | Nitrofurantoin (0.2) |
| Respiratory (3) | Macrolides (0.95), FQ (0.95), oxazolidinones (0.9) | Daptomycin (0.0), AG (0.4) |
| Bloodstream (4) | All 1.0 (reference compartment) | — |
| Intra-abdominal (5) | Metronidazole (0.9), FQ (0.75), carbapenems (0.75) | AG (0.3) |
| CNS (6) | Metronidazole (0.80), oxazolidinones (0.70), chloramphenicol (0.70) | AG (0.05), colistin (0.05), daptomycin (0.05) |
| GI (7) | Fidaxomicin (1.0), metronidazole (0.95), oral vancomycin (0.90) | Glycopeptides IV (0.35) |
| Genital (8) | FQ (0.9), metronidazole (0.8), TMP-SMX (0.8) | AG (0.35) |
| Bone/joint (9) | Rifampicin (0.80), oxazolidinones (0.75), FQ (0.70) | AG (0.25), colistin (0.2) |



### 6.5 Drug Potency Matrix

Intrinsic drug activity against each bacterium (when no resistance is present) is defined in a 42×52 potency matrix. Values range from 0.0 (no activity) to 1.0 (maximum activity).

The potency matrix is embedded as `POTENCY_EMBEDDED_DATA` in config.rs and generates parameters with the key pattern:

```
drug_{drug}_for_bacteria_{bacteria}_potency_when_no_r
```

For example, meropenem has potency 0.95 against E. coli but 0.0 against MRSA (when mecA is present). Penicillin G has potency 0.90 against S. pneumoniae but 0.0 against P. aeruginosa (intrinsically resistant).



### 6.6 Drug Availability by Region and Era

Drug availability varies by region and becomes available at historical introduction dates:

| Variable pattern | Range | Description |
|------------------|-------|-------------|
| `{region}_drug_{drug}_availability` | 0.0–1.0 | Regional availability (1.0 = fully available) |

**Regional availability patterns:**

| Region | Pattern |
|--------|---------|
| North America | All drugs 1.0 |
| Europe | All drugs 1.0 |
| Asia | Most 1.0; tedizolid/ceftaroline 0.3, teicoplanin 0.7 |
| Oceania | Most 1.0; tedizolid/ceftaroline 0.5 |
| South America | Limited newer drugs (tedizolid 0.1, linezolid 0.5, carbapenems 0.6–0.7) |
| Africa | Basic antibiotics 0.8–1.0; ceftriaxone 0.6; vancomycin 0.3; carbapenems 0.1–0.2; most newer drugs 0.0–0.1 |



#### Drug introduction dates

Each antibiotic becomes available in a specific year:

|Drug|~Year|
|------|-------|
|sulfanilamide|1937|
|penicillin_g|1942|
|chloramphenicol|1949|
|tetracycline|1948|
|colistin|1952|
|erythromycin|1952|
|nitrofurantoin|1953|
|furazolidone|1955|
|vancomycin|1958|
|fosfomycin|1959|
|metronidazole|1960|
|ampicillin|1961|
|fusidic_a|1962|
|gentamicin|1963|
|rifampicin|1966|
|doxycycline|1967|
|clindamycin|1968|
|trim_sulf|1968|
|amoxicillin|1972|
|cephalexin|1970|
|minocycline|1971|
|cefazolin|1973|
|tobramycin|1975|
|amikacin|1976|
|ticarcillin|1977|
|cefuroxime|1978|
|piperacillin|1981|
|ceftriaxone|1984|
|piperacillin_tazobactam|1984|
|ceftazidime|1985|
|imipenem_c|1985|
|amoxicillin_clavulanate|1985|
|aztreonam|1986|
|ciprofloxacin|1987|
|teicoplanin|1988|
|ampicillin_sulbactam|1990|
|ticarcillin_clavulanate|1990|
|clarithromycin|1990|
|ofloxacin|1990|
|azithromycin|1991|
|cefepime|1996|
|meropenem|1996|
|levofloxacin|1996|
|moxifloxacin|1999|
|quinu_dalfo|1999|
|linezolid|2000|
|ertapenem|2001|
|daptomycin|2005|
|ceftazidime_avibactam|2006|
|tigecycline|2007|
|retapamulin|2007|
|ceftaroline|2010|
|fidaxomicin|2011|
|tedizolid|2014|
|dalbavancin|2014|
|ceftolozane_tazobactam|2014|
|meropenem_vaborbactam|2018|
|cefiderocol|2019|

**Special case — Colistin**: Withdrawn from routine use between ~1970 and ~1995 (availability drops to 5% during that window), reflecting historical concerns about nephrotoxicity before its re-adoption for MDR Gram-negative infections.



### 6.7 Drug Toxicity

Drugs accumulate a toxicity reservoir, which can trigger sub-lethal discontinuation or lethal toxicity death:

| Parameter | Baseline Value | Description |
|----------|---------|-------------|
| `default_drug_toxicity_death_hazard_per_unit_level` | 0.0 | Per-drug hazard contribution (most drugs = 0; high for colistin, aminoglycosides) |
| `default_toxicity_reservoir_half_life_days` | 1.5 | Decay rate of accumulated toxicity |



#### Toxicity discontinuation (sub-lethal)

When accumulated toxicity exceeds a threshold, the treating clinician may discontinue the drug:

| Parameter | Baseline (Log-odds ratio) | Description |
|----------|---------|-------------|
| `toxicity_discontinuation_base_log_odds` | −3.0 | Baseline discontinuation probability |
| `toxicity_discontinuation_log_odds_per_reservoir_unit` | 1.5 | Per unit of toxicity reservoir |
| `toxicity_discontinuation_log_odds_sepsis` | −1.5 | Clinicians tolerate more toxicity during sepsis |
| `toxicity_avoidance_penalty_multiplier` | 0.05 | Recently-stopped drugs are avoided in reselection |
| `toxicity_avoidance_window_days` | 14.0 | Duration of avoidance penalty |



#### Lethal toxicity death

| Parameter | Baseline (Log-odds ratio) | Description |
|----------|---------|-------------|
| `toxicity_death_base_log_odds` | −8.0 | Very low baseline |
| `toxicity_death_log_odds_per_reservoir_unit` | 2.0 | Per unit of toxicity |
| `toxicity_death_log_odds_age_infant` | 0.6 | |
| `toxicity_death_log_odds_age_child` | 0.2 | |
| `toxicity_death_log_odds_age_adult` | 0.0 | Reference |
| `toxicity_death_log_odds_age_elderly` | 0.8 | |
| `toxicity_death_log_odds_immunosuppressed` | 0.9 | |
| `toxicity_death_log_odds_hospitalized` | 0.25 | Slight increase (monitoring enables detection but may indicate frailty) |



### 6.8 Antibiotic infection prevention

When an individual is actively taking an antibiotic, the treatment acts as prophylaxis against new incoming infections. The efficacy of this prevention represents the proportional reduction in infection acquisition probability for incoming sensitive strains.

| Parameter | Baseline Value | Description |
|-----------|----------------|-------------|
| ntibiotic_infection_prevention_efficacy | 0.7 | The relative reduction in new infection establishment probability when the individual is already taking an effective antibiotic. |

---


## 7. Resistance Dynamics



### 7.1 Resistance Mechanisms

The model explicitly tracks exactly **40** distinct resistance mechanisms, each representing a biological pathway for antibiotic resistance. These mechanisms map dynamically to specific therapeutic drugs, intrinsically lowering their efficacy when present. Based on the `mechanism_applies_to_drug` evaluation rules in the source engine, the explicit drugs affected by each mechanism are:

  | Mechanism | Variable name | Description | Explicit Drugs Affected | Bacterial Classes Affected |
  |-----------|--------------|-------------|-------------------------|----------------------------|
  | ESBL CTX-M | `esbl_ctx_m` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | ESBL TEM | `esbl_tem` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | ESBL SHV | `esbl_shv` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AmpC CMY | `ampc_cmy` | Plasmid-mediated AmpC β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AmpC DHA | `ampc_dha` | Plasmid-mediated AmpC β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | KPC | `kpc` | *K. pneumoniae* carbapenemase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | NDM/VIM | `ndm_vim` | Metallo-β-lactamases | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `ceftazidime_avibactam`, `meropenem_vaborbactam`, `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | OXA-48 | `oxa_48` | Oxacillinase-type carbapenemase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `meropenem`, `imipenem_c`, `ertapenem`, `meropenem_vaborbactam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | OXA-Acinetob. | `oxa_acinetobacter` | OXA-23/40/58 carbapenemases (A. baumannii) | `meropenem`, `imipenem_c`, `ertapenem`, `ceftazidime`, `cefepime`, `ceftazidime_avibactam` | Nonfermenters |
  | blaZ | `blaz` | Staphylococcal penicillinase | `penicillin_g`, `ampicillin`, `amoxicillin` | Gram-Positives |
  | PBP2a/MecA | `pbp2a_meca` | PBP alteration (MRSA) | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `cefiderocol`, `ceftazidime_avibactam`, `meropenem_vaborbactam`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem` | Gram-Positives, Helicobacter |
  | VanA | `vana` | High-level vancomycin resistance | `vancomycin`, `teicoplanin`, `dalbavancin` | Gram-Positives, Helicobacter |
  | VanB | `vanb` | Variable-level vancomycin resistance | `vancomycin` | Gram-Positives, Helicobacter |
  | GyrA (pri.) | `gyra_primary` | DNA gyrase mutation (step 1) | `ciprofloxacin`, `ofloxacin` | All |
  | GyrA + ParC | `gyra_parc` | Additional topoisomerase mutation | `ciprofloxacin`, `ofloxacin`, `levofloxacin`, `moxifloxacin` | All |
  | Qnr | `qnr` | Quinolone resistance protein | `ciprofloxacin`, `ofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | 16S rRMT | `16s_rrmt` | 16S rRNA methyltransferase | `gentamicin`, `tobramycin`, `amikacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AAC/APH/ANT | `aac_aph` | Aminoglycoside-modifying enzymes | `gentamicin`, `tobramycin`, `amikacin`, `streptomycin`, `neomycin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Gram-Positives |
  | ErmB | `ermb` | Erythromycin ribosome methylase | `erythromycin`, `azithromycin`, `clarithromycin`, `clindamycin`, `quinu_dalfo` | Gram-Positives, Anaerobes, Fastidious, Helicobacter |
  | 23S rRNA | `23s_rrna` | 23S rRNA point mutation | `erythromycin`, `azithromycin`, `clarithromycin` | Helicobacter, Enteric Pathogens, Fastidious, Gram-Positives |
  | Cfr | `cfr` | 23S rRNA methyltransferase | `linezolid`, `tedizolid`, `chloramphenicol`, `clindamycin`, `retapamulin` | Gram-Positives, Anaerobes, Fastidious, Helicobacter |
  | CAT | `cat` | Chloramphenicol acetyltransferase | `chloramphenicol` | All |
  | MCR-1 | `mcr_1` | Mobilised colistin resistance | `colistin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AcrAB-TolC | `acrab_tolc` | Gram-negative efflux pump | `tetracycline`, `doxycycline`, `minocycline`, `tigecycline`, `chloramphenicol`, `ciprofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | MexXY-OprM | `mexxy_oprm` | Pseudomonas-specific efflux pump | `tetracycline`, `doxycycline`, `minocycline`, `gentamicin`, `tobramycin`, `amikacin`, `chloramphenicol`, `ciprofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | Global eff. | `global_efflux` | Non-specific efflux upregulation | `tetracycline`, `doxycycline`, `minocycline`, `tigecycline`, `chloramphenicol`, `ciprofloxacin` | All |
  | TetA/B/C | `tet_abc` | Gram-negative tetracycline efflux | `tetracycline`, `doxycycline` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious |
  | TetM/TetO | `tetm` | Ribosomal protection | `tetracycline`, `doxycycline`, `minocycline` | All |
  | OmpK35/36 | `ompk35_36` | Outer membrane porin loss (Klebsiella) | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `ceftaroline`, `cefiderocol`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem`, `ciprofloxacin`, `levofloxacin`, `moxifloxacin`, `ofloxacin`, `gentamicin`, `tobramycin`, `amikacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | OprD | `oprd` | Outer membrane porin loss (Pseudomonas) | `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | Global por. | `global_porin_loss` | Non-specific porin downregulation | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `ceftaroline`, `cefiderocol`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem`, `ciprofloxacin`, `levofloxacin`, `moxifloxacin`, `ofloxacin`, `gentamicin`, `tobramycin`, `amikacin` | All |
  | Folate path | `folate_pathway` | Altered dihydrofolate reductase | `sulfanilamide`, `trim_sulf` | All |
  | Nitroreduct | `nitroreductase` | Nitroreductase loss | `metronidazole`, `nitrofurantoin`, `furazolidone` | Enterobacterales, Enteric Pathogens, Anaerobes, Fastidious, Helicobacter |
  | FosA | `fosa` | Fosfomycin-modifying enzyme | `fosfomycin` | Enterobacterales, Nonfermenters, Enteric Pathogens |
  | MprF | `mprf` | Membrane charge modification | `daptomycin` | Gram-Positives |
  | RpoB | `rpob` | RNA polymerase mutation | `fidaxomicin` | All |
  | FusB | `fusb` | Fusidic acid resistance determinant | `fusidic_a` | Gram-Positives |
  | Unknown 1–3 | `as_yet_unknown_{1..3}` | Placeholder mechanisms | Evaluates `true` dynamically for all applied overrides | All (Calibration Placeholders) |



### 7.2 Mechanism–Drug-Class Enhancement Multipliers

Each mechanism reduces the efficacy of specific drug classes. The enhancement multiplier (0.0–1.0) determines how much a mechanism reduces drug activity: 0.0 = no effect, 1.0 = complete resistance.

These are defined per mechanism × drug class (40 × 36 = 1440 values), with variable pattern:

```
mech_{mechanism}_enhancement_{drug_class}
```

**Global (legacy) enhancement multipliers** (used as a default fallback when a specific per-class mathematical override hasn't been explicitly configured yet in the internal Rust parameters):

| Mechanism | Global multiplier |
|-----------|------------------|
| ESBL CTX-M | 0.85 |
| ESBL TEM | 0.80 |
| ESBL SHV | 0.75 |
| AmpC CMY/DHA | 0.70 |
| KPC | 0.90 |
| NDM/VIM | 0.95 |
| OXA-48 | 0.80 |
| PBP2a/MecA | 0.90 |
| VanA | 0.95 |
| VanB | 0.85 |
| GyrA primary | 0.70 |
| GyrA + ParC | 0.85 |
| ErmB | 0.80 |
| Cfr | 0.75 |
| 16S rRMT | 0.85 |
| CAT | 0.70 |
| Qnr | 0.50 |
| MCR-1 | 0.60 |
| AcrAB-TolC | 0.40 |
| MexXY-OprM | 0.45 |
| OmpK35/36 | 0.50 |
| OprD | 0.55 |
| Global efflux | 0.35 |
| Global porin loss | 0.45 |
| Folate pathway | 0.70 |
| Nitroreductase | 0.60 |
| FosA | 0.65 |
| MprF | 0.55 |
| RpoB | 0.80 |
| FusB | 0.70 |
| As-yet-unknown 1–4 | 0.50 each |



### 7.3 Resistance Emergence

De novo resistance emergence occurs strictly when an individual is under antibiotic pressure. The probability is calculated daily as a function of the underlying capability of the pathogen to mutate toward that profile, mediated by current drug concentrations and bacterial population size.

To understand emergence and treatment dynamics, two core terms must be defined:
- **Potency (`base_potency`)**: The intrinsic baseline efficacy of a specific drug against a specific bacterium (e.g., `potency_e_coli_ampicillin`).
- **Activity (`activity_r`)**: The dynamically calculated killing effect of a drug on a given day. It is defined mathematically as:  
  `activity_r = base_potency × effective_drug_level × (1.0 - normalized_any_r)`  
  *(where `normalized_any_r` represents the current resistance level bounded between 0.0 and 1.0, and `effective_drug_level` accounts for pharmacokinetic decay and tissue penetration).*

The daily probability of a specific resistance mechanism emerging during an active infection is evaluated using the following expression:

```rust
emergence_rate = mechanism_rate 
               * (1.0 + bacteria_level_factor) 
               * max_emergence_drug_factor 
               * multi_drug_penalty_factor
```

Where:
- **`mechanism_rate`**: The baseline biological capability of the pathogen to acquire this specific mechanism (`bacteria_{name}_mechanism_{mech}_emergence_rate`).
- **`bacteria_level_factor`**: A logarithmic scaling factor representing population size. Larger symptomatic infections have a proportionally higher absolute number of replicating cells, mathematically increasing the probability of a successful mutation/HGT event.
- **`max_emergence_drug_factor`**: A Gaussian curve peaking at exactly 50% of the standard effective tissue concentration ($\sigma = 0.2$). This explicitly models the **sub-inhibitory concentration danger zone**: emergence probabilities are highest when drug levels are sufficient to provide strong selective pressure, but too low to rapidly eradicate the mutant. At peak/max concentrations, emergence is heavily suppressed (falling back to a `0.01` baseline).
- **`multi_drug_penalty_factor`**: A biological suppression modifier applied when multiple active therapies are used simultaneously, significantly reducing the probability of a single-mechanism bypass if the combined therapy includes a drug unaffected by that mechanism.



#### Incidence Band Multipliers

**Epidemiological Scaling of Mutation Rates**

In population-level dynamic modelling, high-incidence pathogens (e.g., *Escherichia coli*) present a computational challenge. If baseline mutational probabilities are applied uniformly across a simulated population of millions without adjustment for absolute event frequency, ubiquitously carried organisms would accumulate *de novo* resistance mutations at a biologically unrealistic velocity, rapidly trending towards pan-resistance. This computational artifact occurs because real-world ecological constraints—such as clonal competition, mucosal immunity, and fitness costs—limit exponential fixation in ways that raw stochastic simulations might otherwise bypass.

To calibrate the underlying model to reflect real-world epidemiological trajectories, modelled bacteria are stratified into 'Incidence Bands' based on their baseline populational frequency and carriage prevalence. Organisms with extensive endemic incidence receive a mathematical penalty (down-scaling of the emergence rate multiplier) to preserve temporal realism and restrict hyper-mutation. Conversely, rare or highly specialised pathogens receive comparative up-scaling to ensure that emergent resistance remains observable during the simulation timeline. The 42 modelled species are categorised as follows:

| Band | Multiplier | All Assigned Bacteria |
|------|-----------|-----------------|
| **High incidence** (Most common) | ×0.1 | *E. coli*, *S. aureus*, *S. pneumoniae*, *C. trachomatis*, *K. pneumoniae*, *H. pylori*, *H. influenzae*, *M. pneumoniae*, *S. epidermidis* |
| **Moderate incidence** | ×1.0 | *S. pyogenes*, *T. pallidum*, *C. jejuni*, *N. gonorrhoeae*, *M. catarrhalis*, *S. agalactiae*, *Enterobacter spp.*, *Enterobacter cloacae*, *Proteus spp.*, *M. genitalium* |
| **Low incidence** | ×3.0 | *P. aeruginosa*, *C. difficile*, *B. fragilis*, *Citrobacter spp.*, *Serratia spp.*, *Salmonella (non-typhoidal)*, *Shigella spp.*, *E. faecalis*, *B. pertussis*, *L. pneumophila* |
| **Very low incidence** (Rarest cases) | ×10.0 | *A. baumannii*, *E. faecium (VRE)*, *MDR M. tuberculosis*, *N. meningitidis*, *S. Typhi*, *S. Paratyphi A*, *V. cholerae*, *L. monocytogenes*, *Y. enterocolitica*, *S. maltophilia*, *Morganella spp.*, *P. stuartii*, *B. cepacia complex* |

Base emergence rates range from highly rare point mutations to more common events like swapping jumping genes (plasmids/transposons). These specific target-site mutations possess strictly defined rate values in the code: `bacteria_{name}_mechanism_{mech}_emergence_rate`.

It is important to note that most of these rates are set to exactly `0.0` for entirely impossible biological events (for example, there is a zero percent chance that *S. pyogenes* will ever suddenly evolve the NDM/VIM carbapenemase mechanism).

- **Key in code**: `bacteria_{bacteria}_mechanism_{mechanism}_emergence_rate`
- **Mechanism assignment probability**: `mechanism_assignment_probability` = 0.8. When the math successfully triggers a new resistance event, there is an 80% chance the specific structural genome mechanism is formally assigned, otherwise it acts as a generalized fallback.




### 7.4 Resistance Reversion and Biological Fitness Costs

In the prolonged absence of antimicrobial selective pressure—typically within an untreated host’s commensal microbiome—resistant phenotypes often manifest a biological fitness cost compared to wild-type, fully susceptible strains. Over time, resistant sub-populations are out-competed by their susceptible counterparts, leading to a phenotypic reversion. 

To model this dynamic, the simulator calculates a daily stoichiometric `reversion_rate` representing the probabilistic loss or suppression of resistance determinants. Crucially, the model maps this reversion risk to the underlying *resistance mechanism* rather than the drug itself. 

For instance, the energetic burden of maintaining large, plasmid-mediated metallo-$eta$-lactamases characteristically exerts different fitness costs compared to purely chromosomal point mutations (such as single-step *gyrA* mutations in fluoroquinolone resistance). If a pathogen acquires a generic, mechanism-less resistance state, the model defaults to a baseline `mechanismless_resistance_reversion_rate` of `0.0004` per day.

The engine configures precise physiological fitness costs (reversion rates per day) for the following modelled mechanisms:

### Enzymatic Inactivation
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| **KPC** (*bla*~KPC~) | `0.001` | Plasmid-mediated carbapenemase; moderate maintenance cost. |
| **NDM / VIM** | `0.0015` | Metallo-$eta$-lactamases, frequently on large, high-burden mobile genetic elements. |
| **OXA-48** | `0.0005` | Class D carbapenemase; comparatively lower fitness burden. |
| **ESBL CTX-M / TEM / SHV** | `0.0006` | Standard extended-spectrum $eta$-lactamases. |
| **AmpC DHA** | `0.0006` | Plasmid-mediated AmpC; typical cost profile. |
| **AmpC CMY** | `0.0001` | Often native gene upregulation; minimal fitness loss to maintain. |
| **FosA** | `0.0005` | Plasmid-mediated fosfomycin resistance; moderate cost. |
| **CAT** | `0.0005` | Chloramphenicol acetyltransferase. |
| **16S rRMTase** | `0.0005` | Ribosomal RNA methyltransferases conferring high-level aminoglycoside resistance. |

### Target Site Alterations
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| **PBP2a / *mecA*** | `0.0009` | High energetic cost associated with maintaining the staphylococcal cassette chromosome *mec* (SCC*mec*). |
| ***erm(B)*** | `0.002` | High reversion rate; target methylation for macrolide-lincosamide-streptogramin B (MLS~B~) resistance. |
| **VanA / VanB** | `0.002` | Highly complex target reprograming (D-Ala-D-Ala to D-Ala-D-Lac); significant energetic drain in the absence of glycopeptide exposure. |
| **CFR** | `0.0005` | RNA methyltransferase (oxazolidinone/phenicol cross-resistance). |

### Structural Mutations
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| ***gyrA* (Primary)** | `0.0001` | Single-step topoisomerase mutation; essentially stable with negligible fitness penalty. |
| ***parC* (Secondary)** | `0.0002` | Secondary topoisomerase IV mutations; compounding structural cost. |
| **Folate Pathway** | `0.0001` | Low cost; largely integron-associated (e.g., *sul* or *dfrA* elements). |
| **Nitroreductase** | `0.0003` | Loss-of-function mutations affecting nitrofurantoin activation; moderate cost. |
| ***mprF*** | `0.001` | Membrane lipid modification (daptomycin resistance); structural cell wall alterations bear a distinctive fitness cost. |
| ***rpoB*** | `0.002` | High fitness cost due to structurally significant alterations in RNA polymerase (rifampicin resistance). |

### Target Protection & Target Modification
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| ***mcr-1*** | `0.0015` | Plasmid-mediated phosphoethanolamine transferase (colistin resistance); substantial lipid A modification burden. |
| **Tet(M)** | `0.0005` | Ribosomal protection protein; moderate cost to maintain on transposons (e.g., Tn*916*). |
| **Qnr** | `0.0001` | Plasmid-mediated quinolone resistance protein; relatively stable. |
| **FusB** | `0.0005` | Target protection mechanism for fusidic acid resistance. |

### Porin Loss & Efflux Pumps
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| **OprD Loss** | `0.0005` | Loss of outer membrane channel (carbapenem resistance); hinders nutrient acquisition. |
| **OmpK35 / OmpK36 Loss** | `0.0005` | Analogous mechanism in Enterobacterales. |
| **Global / Generic Porin Loss** | `0.0005` | Broad permeability phenotype cost. |
| **AcrAB-TolC** | `0.0005` | Overexpression of major RND-family efflux pump complex. |
| **MexXY-OprM** | `0.0005` | Endogenous efflux system upregulation (common in *Pseudomonas aeruginosa*). |
| **Global / Generic Efflux** | `0.0005` | Broad, non-specific transport energy costs. |

*Note: The system reserves three placeholder variables (`as_yet_unknown_1-3`, baseline rate `0.001`) designated for future, non-prescribed clinical trials or empirical calibration against newly emergent paradigms.*

## 7.5 Resistance Floors

**Addressing Stochastic Fade-out in Rare Pathogens**

In a mathematical simulation of a mid-sized population (e.g., 100,000 individuals), clinically important but low-incidence pathogens (such as *Stenotrophomonas maltophilia* or *Enterococcus faecium*) generate sparse data points. Without intervention, their historically established baseline resistance traits (such as intrinsic genes or highly prevalent endemic plasmids) can mathematically "fade out" because there are simply not enough infection events to sustain the resistance reservoir via selection pressure alone.

To counteract this computational decay and reliably enforce known clinical phenotypes, the model applies **"Resistance Floors"**—hard minimum bounds on the resistance prevalence (`majority_r`) that ramp up chronologically following the historical introduction of respective drug classes.

| Parameter pattern | Baseline Value | Description |
|------------------|---------|-------------|
| `resistance_floor_feature_enabled` | `1.0` | Global engine toggle. |
| `bacteria_{name}_resistance_floor_enabled` | `>0.5` enables | Per-organism toggle. |
| `bacteria_{name}_resistance_floor_ramp_years` | `10.0` | Velocity of resistance fixation (years from drug class introduction to full basal floor). |
| `bacteria_{name}_{drug_class}_resistance_floor` | `0.0` - `1.0`| Target minimum baseline resistance prevalence. |

The floor ramps linearly from zero to the target prevalence over the selected ramp period, strictly gating the starting date to the historical introduction of the earliest therapeutic agent in the affected class.

**Currently Configured Organisms:**

1. ***Stenotrophomonas maltophilia***: **Enabled** (`1.0`). S. maltophilia has a low absolute incidence in the simulated cohorts but carries profound intrinsic resistance. The floor explicitly preserves its intrinsic characteristics, such as the L1/L2 metallo-β-lactamases and multidrug efflux pumps, preventing the model from erroneously rendering the bacteria susceptible to carbapenems or generalized cephalosporins dynamically. 
2. ***Enterococcus faecium***: **Disabled** (`0.0`). The infrastructure exists to enforce resistance floors for *E. faecium*—specifically designed to sustain high-level ampC transitions and vancomycin resistance (VRE) baselines—but the parameters are currently toggled off in the live default configuration.



### 7.6 Cross-Resistance Groups

For each bacteria, drugs that share resistance mechanisms are grouped together. Acquiring resistance to one drug in a cross-resistance group confers resistance to all drugs in that group.

| Bacteria | Group | Drugs sharing resistance |
|----------|-------|------------------------|
| E. coli | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Cephalexin, Cefazolin, Cefuroxime, Ceftriaxone, Amoxicillin Clavulanate, Ampicillin Sulbactam, Piperacillin Tazobactam, Ticarcillin Clavulanate |
| E. coli | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| E. coli | Group 3 | Gentamicin, Tobramycin, Amikacin |
| A. baumannii | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Cephalexin, Cefazolin, Cefuroxime, Amoxicillin Clavulanate, Ampicillin Sulbactam, Piperacillin Tazobactam, Ticarcillin Clavulanate |
| A. baumannii | Group 2 | Meropenem, Imipenem C, Ertapenem, Meropenem Vaborbactam |
| A. baumannii | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin |
| A. baumannii | Group 4 | Gentamicin, Tobramycin, Amikacin |
| K. pneumoniae | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Cephalexin, Cefazolin, Cefuroxime, Ceftriaxone, Amoxicillin Clavulanate, Ampicillin Sulbactam, Piperacillin Tazobactam, Ticarcillin Clavulanate |
| K. pneumoniae | Group 2 | Meropenem, Imipenem C, Ertapenem, Meropenem Vaborbactam |
| K. pneumoniae | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. pneumoniae | Group 1 | Erythromycin, Azithromycin, Clarithromycin |
| S. pneumoniae | Group 2 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| S. aureus | Group 1 | Penicillin G, Ampicillin, Amoxicillin |
| S. aureus | Group 2 | Cephalexin, Cefazolin, Cefuroxime, Ceftriaxone |
| S. aureus | Group 3 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| S. epidermidis | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Cephalexin, Cefazolin, Cefuroxime, Ceftriaxone |
| S. epidermidis | Group 2 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| S. epidermidis | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. maltophilia | Group 1 | Trim Sulf |
| S. maltophilia | Group 2 | Tetracycline, Doxycycline, Minocycline |
| S. maltophilia | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| P. aeruginosa | Group 1 | Piperacillin, Piperacillin Tazobactam, Ceftazidime, Ceftazidime Avibactam, Cefepime |
| P. aeruginosa | Group 2 | Meropenem, Meropenem Vaborbactam, Imipenem C |
| P. aeruginosa | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| P. aeruginosa | Group 4 | Gentamicin, Tobramycin, Amikacin |
| E. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| E. spp. | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| E. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Mdr mycobacterium tuberculosis | Group 1 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Mdr mycobacterium tuberculosis | Group 2 | Gentamicin, Tobramycin, Amikacin |
| Mdr mycobacterium tuberculosis | Group 3 | Rifampicin |
| Mdr mycobacterium tuberculosis | Group 4 | Linezolid |
| E. faecalis | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Piperacillin, Piperacillin Tazobactam, Ticarcillin, Ticarcillin Clavulanate |
| E. faecalis | Group 2 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| E. faecalis | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| E. faecalis | Group 4 | Tetracycline, Doxycycline, Minocycline |
| E. faecalis | Group 5 | Vancomycin, Teicoplanin, Dalbavancin |
| E. faecium | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Piperacillin, Piperacillin Tazobactam |
| E. faecium | Group 2 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| E. faecium | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| E. faecium | Group 4 | Tetracycline, Doxycycline, Minocycline |
| E. faecium | Group 5 | Vancomycin, Teicoplanin, Dalbavancin |
| C. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| C. spp. | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| C. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| C. spp. | Group 4 | Gentamicin, Tobramycin, Amikacin |
| E. cloacae | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| E. cloacae | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| E. cloacae | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| M. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| M. spp. | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| M. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| P. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| P. spp. | Group 2 | Ceftriaxone, Ceftazidime |
| P. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| S. spp. | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| S. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. spp. | Group 4 | Gentamicin, Tobramycin, Amikacin |
| P. stuartii | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin |
| P. stuartii | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Salmonella enterica serovar typhi | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| Salmonella enterica serovar typhi | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Salmonella enterica serovar typhi | Group 3 | Ceftriaxone, Ceftazidime |
| Salmonella enterica serovar paratyphi a | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| Salmonella enterica serovar paratyphi a | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Salmonella enterica serovar paratyphi a | Group 3 | Ceftriaxone, Ceftazidime |
| Invasive non-typhoidal salmonella spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| Invasive non-typhoidal salmonella spp. | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Invasive non-typhoidal salmonella spp. | Group 3 | Ceftriaxone, Ceftazidime |
| S. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| S. spp. | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. spp. | Group 3 | Ceftriaxone, Ceftazidime |
| S. spp. | Group 4 | Tetracycline, Doxycycline |
| V. cholerae | Group 1 | Tetracycline, Doxycycline |
| V. cholerae | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| V. cholerae | Group 3 | Erythromycin, Azithromycin, Clarithromycin |
| C. jejuni | Group 1 | Erythromycin, Azithromycin, Clarithromycin |
| C. jejuni | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| C. jejuni | Group 3 | Tetracycline, Doxycycline |
| Y. enterocolitica | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| Y. enterocolitica | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Y. enterocolitica | Group 3 | Tetracycline, Doxycycline |
| H. pylori | Group 1 | Clarithromycin, Erythromycin, Azithromycin |
| H. pylori | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| H. pylori | Group 3 | Tetracycline, Doxycycline |
| S. pyogenes | Group 1 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| S. pyogenes | Group 2 | Tetracycline, Doxycycline |
| S. pyogenes | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. agalactiae | Group 1 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| S. agalactiae | Group 2 | Tetracycline, Doxycycline |
| S. agalactiae | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| H. influenzae | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| H. influenzae | Group 2 | Erythromycin, Azithromycin, Clarithromycin |
| H. influenzae | Group 3 | Ciprofloxacin, Levofloxacin |
| H. influenzae | Group 4 | Tetracycline, Doxycycline |
| M. catarrhalis | Group 1 | Ampicillin, Amoxicillin, Penicillin G, Piperacillin, Ticarcillin, Amoxicillin Clavulanate, Ampicillin Sulbactam, Piperacillin Tazobactam, Ticarcillin Clavulanate |
| M. catarrhalis | Group 2 | Erythromycin, Azithromycin, Clarithromycin |
| N. gonorrhoeae | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| N. gonorrhoeae | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| N. gonorrhoeae | Group 3 | Tetracycline, Doxycycline |
| N. gonorrhoeae | Group 4 | Erythromycin, Azithromycin |
| N. gonorrhoeae | Group 5 | Ceftriaxone, Cefixime |
| N. meningitidis | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| N. meningitidis | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| N. meningitidis | Group 3 | Rifampicin |
| C. difficile | Group 1 | Vancomycin |
| C. difficile | Group 2 | Metronidazole |
| B. fragilis | Group 1 | Metronidazole |
| B. fragilis | Group 2 | Clindamycin |
| B. fragilis | Group 3 | Meropenem, Imipenem C |
| L. monocytogenes | Group 1 | Ampicillin, Amoxicillin, Penicillin G, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| L. monocytogenes | Group 2 | Tetracycline, Doxycycline |
| C. trachomatis | Group 1 | Azithromycin, Erythromycin, Clarithromycin |
| C. trachomatis | Group 2 | Doxycycline, Tetracycline |
| C. trachomatis | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| M. genitalium | Group 1 | Azithromycin, Erythromycin, Clarithromycin |
| M. genitalium | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| M. genitalium | Group 3 | Doxycycline, Tetracycline |
| T. pallidum | Group 1 | Erythromycin, Azithromycin, Clarithromycin |
| T. pallidum | Group 2 | Doxycycline, Tetracycline |
| B. pertussis | Group 1 | Erythromycin, Azithromycin, Clarithromycin |


---

## 8. Microbiome and Carriage



### 8.1 Overview

Each individual maintains a per-bacteria microbiome resistance state. Bacteria can colonise the microbiome (carriage) without causing active infection. Carriage is the primary reservoir for resistance transmission.



### 8.2 Carriage Compartments

Each bacterial species is assigned to a carriage compartment reflecting its natural ecological niche:

| Compartment | Example bacteria |
|-------------|-----------------|
| Gut | E. coli, K. pneumoniae, Enterococcus spp., Shigella, Salmonella, C. difficile |
| Respiratory | S. pneumoniae, H. influenzae, P. aeruginosa, A. baumannii, M. catarrhalis, M. tuberculosis |
| Skin/Soft Tissue | S. aureus, S. epidermidis |
| Genitourinary | N. gonorrhoeae, C. trachomatis, M. genitalium, T. pallidum, S. agalactiae |
| Systemic | (Reserved; not currently assigned to any bacteria) |



### 8.3 Resistance in the Microbiome

The microbiome serves as the primary reservoir for "hidden" microbial resistance transmission. The engine calculates an asymptomatic `microbiome_r` tracking matrix for every organism in every patient.

**Governing Resistance Dynamics:**
- Under antibiotic eradication regimes, clearing established resistant micro-colonies is significantly harder the longer they have existed (`carriage_duration_log_odds_coefficient` = `-0.01` per day, caps at `-2.0` rendering old colonies ~7.4x harder to eliminate).
- When a carrier structurally contracts a true clinical case from their own background flora, `infection_from_microbiome_dampening` = `0.70` attenuates the sheer volume of infection to accurately mirror biological auto-infection likelihoods. 
- *De novo* evolution inside an asymptomatic microbiome strictly lacking therapy is physiologically zero (`microbiome_resistance_emergence_rate_per_day_baseline` = `1.0e-20`), dictating that macroscopic resistance is driven by infected individuals receiving therapy, not silent environmental mutators.



## 9. Horizontal Gene Transfer

The model simulates Horizontal Gene Transfer (HGT) dynamically as a primary driver of resistance spread, distinct from clonal expansion.



### HGT Matrices
Rather than random probability events, the exchange of plasmids is strictly governed by precomputed transmission compatibility matrices representing the probability of a resistance element transferring between a donor and recipient strain in vivo.
Within the model, the `hgt` matrix establishes a scalar transmission capacity `HGT_{(A, B)}` representing flow from Bacteria A to Bacteria B.

Key structural boundaries defined by `PlasmidPool` mapping include:
- **Intra-species transfer** (e.g., E. coli to E. coli): `0.95` (virtually unobstructed plasmid exchange).
- **Inter-genus Gram-Negative transfer** (e.g., E. coli to Klebsiella): `0.80` to `0.90` (highly fluid exchange of ESBLs and Carbapenemases).
- **Narrow Gram-Positive transfer** (e.g., S. aureus to Enterococcus): `0.10` to `0.20` (heavily restricted, rare transposon integration).
- **Gram-Negative to Gram-Positive transfer**: `0.0` (biologically prohibited boundary).



### The HGT Process
1. **Donor Pool Constraint**: The effective pool of transmissible resistomes requires established presence in either the active clinical infection or the background carrier microbiome.
2. **Transfer Rate**: The base probability of an asymptomatic background transfer event occurring daily is calibrated at `microbiome_resistance_transfer_probability_per_day` = `0.0001`.
3. **Antibiotic Pressure Amplifier**: When the host is actively receiving an antibiotic, the resulting microbiome stress and bacterial lysis triggers the SOS response, mechanically increasing conjugative transfer. This multiplies the daily HGT probability by `hgt_antibiotic_pressure_multiplier` = `1.50` (a 50% relative increase in plasmid exchange during active therapy).



## 10. Mortality

Clinical mortality is evaluated daily for individuals suffering active, systemic infection or progressing into sepsis. 



### Baseline Escalation Model
Background population mortality risk rises predictably with age and historic improvements in public health:
- A non-linear aging penalty applies a roughly `4%` increase per year (`log_odds_mortality_per_year_of_age` = `0.04`, `exp(0.04) ≈ 1.04`).
- An additional quadratic escalation captures systemic elderly frailty (`log_odds_mortality_per_year_of_age_squared` = `0.05`).



### Infection Mortality Risk
Bacterial attributes set the foundation for disease lethality (`base_mortality_intercept = -12.0` if not present):
- **Mild risk**: Non-invasive agents (*C. trachomatis* at `-15.0`).
- **Moderate risk**: *E. coli* (`-9.0`), *S. pneumoniae* (`-6.0`).
- **Critical risk**: *N. meningitidis* (`-1.2`), *M. tuberculosis* (`-0.8`).

When an infection structurally ruptures the `symptomatic_infection_level_threshold`, daily mortality checks apply:
$$	ext{mortality\_prob} = 	ext{base}_{	ext{pop}} 	imes f(	ext{duration, level, path}) 	imes 	ext{syndrome\_multiplier}$$

| Syndrome | Mortality Multiplier |
|----------|-----------------------|
| Genital | `0.05` |
| Skin / Ear | `0.1` |
| UTI | `0.5` |
| Bone/Joint | `0.8` |
| Intrabdominal | `1.5` |
| Respiratory | `1.5` |
| CNS | `3.0` |
| Bloodstream | `4.0` |

If the system progresses to sepsis, these limits are aggressively overridden, driving rapid physiological collapse if correct empirical therapy is not structurally active.



## 11. Policy Evaluation



### 11.1 Branching

At a configurable branch year (default 2027), the simulation saves the full population state and runs three independent scenarios forward to the end of the simulation:

| Branch | Description |
|--------|-------------|
| **Baseline** | No policy changes — business as usual |
| **Stewardship** | Antibiotic stewardship adjustments (e.g., narrower prescribing, enhanced testing) |
| **Counterfactual** | Hypothetical scenario (e.g., a world with no resistance) |



### 11.2 Policy Parameters

Each branch can modify:

| Parameter | Baseline | Stewardship | Counterfactual |
|-----------|----------|-------------|----------------|
| `drug_selection_temperature` | — | ×0.65 (more deterministic) | — |
| `minimal_potency_threshold_for_drug_selection` | — | — | — |
| `bacterial_testing_rate_multiplier` | — | ×1.5 | — |
| `resistance_testing_rate_multiplier` | — | ×1.5 | — |
| `counterfactual_resistance_multiplier` | — | — | 0.0 |
| `clear_all_resistance_on_branch_start` | false | false | true |
| `reserve_drug_penalty_multiplier` | — | ×2.0 | — |
| `drug_initiation_rate_multiplier` | — | ×0.85 | — |
| `drug_cessation_rate_multiplier` | — | ×1.2 | — |



### 11.3 Key Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `SIMULATION_START_YEAR` | 1930.0 | Calendar year at day 0 |
| `POLICY_BRANCH_YEAR` | 2027.0 | Year when policies diverge |
| `INFECTION_EPS` | 0.001 | Minimum meaningful infection level |
| `MICROBIOME_MAJORITY_THRESHOLD` | 0.5 | Threshold for majority resistance |
| `MAX_MECHANISM_PROFILES` | 200 | Reservoir sample size per bacteria for mechanism profile cache |

---



## 12. Limitations

- Drug levels are modelled as abstract units rather than pharmacokinetic concentrations
- Bacteria-bacteria competition within the microbiome is represented implicitly through resistance promotion and decay rather than explicit strain dynamics
- The model does not capture within-host spatial heterogeneity (e.g., biofilm vs planktonic)
- Vaccine effects are explicitly modelled: vaccinated people hold a proportionally lower incidence risk, though the infection risk is not dynamically dependent on background prevalence
- Region definitions are broad continental groupings

---



## Appendix A — Bacteria, Drugs, Mechanisms and Enums



### A.1 Bacteria (42 species)

| Index | Species | Group | Carriage compartment |
|-------|---------|-------|---------------------|
| 0 | Acinetobacter baumannii | NonFermenter | Respiratory |
| 1 | Citrobacter spp. | Enterobacterales | Gut |
| 2 | Enterobacter spp. | Enterobacterales | Gut |
| 3 | Enterococcus faecalis | GramPositive | Gut |
| 4 | Enterococcus faecium | GramPositive | Gut |
| 5 | Escherichia coli | Enterobacterales | Gut |
| 6 | Klebsiella pneumoniae | Enterobacterales | Gut |
| 7 | Morganella spp. | Enterobacterales | Gut |
| 8 | Proteus spp. | Enterobacterales | Gut |
| 9 | Serratia spp. | Enterobacterales | Gut |
| 10 | Providencia stuartii | Enterobacterales | Genitourinary |
| 11 | Pseudomonas aeruginosa | NonFermenter | Respiratory |
| 12 | Stenotrophomonas maltophilia | NonFermenter | Respiratory |
| 13 | Staphylococcus aureus | GramPositive | Skin/Soft Tissue |
| 14 | Staphylococcus epidermidis | GramPositive | Skin/Soft Tissue |
| 15 | Streptococcus pneumoniae | GramPositive | Respiratory |
| 16 | Salmonella enterica serovar Typhi | Enterobacterales | Gut |
| 17 | Salmonella enterica serovar Paratyphi A | Enterobacterales | Gut |
| 18 | Invasive non-typhoidal Salmonella spp. | Enterobacterales | Gut |
| 19 | Shigella spp. | Enterobacterales | Gut |
| 20 | Neisseria gonorrhoeae | Fastidious | Genitourinary |
| 21 | Streptococcus pyogenes | GramPositive | Respiratory |
| 22 | Streptococcus agalactiae | GramPositive | Genitourinary |
| 23 | Haemophilus influenzae | Fastidious | Respiratory |
| 24 | Chlamydia trachomatis | Fastidious | Genitourinary |
| 25 | Mycoplasma genitalium | Fastidious | Genitourinary |
| 26 | Vibrio cholerae | EntericPathogen | Gut |
| 27 | Neisseria meningitidis | Fastidious | Respiratory |
| 28 | Listeria monocytogenes | GramPositive | Gut |
| 29 | Clostridioides difficile | Anaerobe | Gut |
| 30 | Bacteroides fragilis | Anaerobe | Gut |
| 31 | Campylobacter jejuni | Helicobacter | Gut |
| 32 | Enterobacter cloacae | Enterobacterales | Gut |
| 33 | Yersinia enterocolitica | Enterobacterales | Gut |
| 34 | Moraxella catarrhalis | Fastidious | Respiratory |
| 35 | Treponema pallidum | Spirochete | Genitourinary |
| 36 | Bordetella pertussis | Fastidious | Respiratory |
| 37 | Helicobacter pylori | Helicobacter | Gut |
| 38 | MDR Mycobacterium tuberculosis | Mycobacteria | Respiratory |
| 39 | Mycoplasma pneumoniae | Fastidious | Respiratory |
| 40 | Legionella pneumophila | Fastidious | Respiratory |
| 41 | Burkholderia cepacia complex | NonFermenter | Respiratory |



### A.2 Antibiotics (58 drugs in 36 classes)

| Drug | Class |
|------|-------|
| sulfanilamide | Sulfonamides |
| penicillin_g | Penicillins |
| ampicillin | Penicillins |
| amoxicillin | Penicillins |
| piperacillin | Penicillins |
| ticarcillin | Penicillins |
| cephalexin | Cephalosporins 1–2G |
| cefazolin | Cephalosporins 1–2G |
| cefuroxime | Cephalosporins 1–2G |
| ceftriaxone | Cephalosporins 3G |
| ceftazidime | Cephalosporins 3G |
| cefepime | Cephalosporins 4G |
| ceftaroline | Anti-MRSA Cephalosporins (5G) |
| ceftolozane_tazobactam | Cephalosporins 3G/BLI |
| cefiderocol | Siderophore Cephalosporins |
| meropenem | Carbapenems |
| imipenem_c | Carbapenems |
| ertapenem | Carbapenems |
| aztreonam | Monobactams |
| erythromycin | Macrolides |
| azithromycin | Macrolides |
| clarithromycin | Macrolides |
| clindamycin | Macrolides |
| gentamicin | Aminoglycosides |
| tobramycin | Aminoglycosides |
| amikacin | Aminoglycosides |
| ciprofloxacin | Fluoroquinolones |
| levofloxacin | Fluoroquinolones |
| moxifloxacin | Fluoroquinolones |
| ofloxacin | Fluoroquinolones |
| tetracycline | Tetracyclines |
| doxycycline | Tetracyclines |
| minocycline | Tetracyclines |
| tigecycline | Tetracyclines |
| vancomycin | Glycopeptides |
| teicoplanin | Glycopeptides |
| dalbavancin | Glycopeptides |
| linezolid | Oxazolidinones |
| tedizolid | Oxazolidinones |
| daptomycin | Other |
| quinu_dalfo | Other |
| trim_sulf | Sulfonamides |
| chloramphenicol | Chloramphenicol |
| nitrofurantoin | Other |
| fosfomycin | Other |
| retapamulin | Other |
| fusidic_a | Other |
| metronidazole | Other |
| fidaxomicin | Other |
| furazolidone | Other |
| rifampicin | Other |
| amoxicillin_clavulanate | BLI Combinations |
| piperacillin_tazobactam | BLI Combinations |
| ampicillin_sulbactam | BLI Combinations |
| ticarcillin_clavulanate | BLI Combinations |
| ceftazidime_avibactam | Novel BL/BLI |
| meropenem_vaborbactam | Novel BL/BLI |
| colistin | Polymyxins |



### A.3 Drug Classes (36)

| Code | Full name | Drugs |
|------|-----------|-------|
| `pen` | Penicillins | penicillin G, ampicillin, amoxicillin, piperacillin, ticarcillin |
| `bli` | BLI Combinations | amoxicillin-clavulanate, piperacillin-tazobactam, ampicillin-sulbactam, ticarcillin-clavulanate |
| `c1_2g` | Cephalosporins 1–2G | cephalexin, cefazolin, cefuroxime |
| `c3g` | Cephalosporins 3G | ceftriaxone, ceftazidime |
| `c3g_bli` | Cephalosporins 3G/BLI | ceftolozane-tazobactam |
| `c4g` | Cephalosporins 4G | cefepime |
| `c5g` | Anti-MRSA Cephalosporins (5G) | ceftaroline |
| `sid_ceph` | Siderophore Cephalosporins | cefiderocol |
| `bl_ni` | Novel BL/BLI | ceftazidime-avibactam, meropenem-vaborbactam |
| `carb` | Carbapenems | meropenem, imipenem, ertapenem |
| `mono` | Monobactams | aztreonam |
| `fq` | Fluoroquinolones | ciprofloxacin, levofloxacin, moxifloxacin, ofloxacin |
| `ag` | Aminoglycosides | gentamicin, tobramycin, amikacin |
| `mls` | Macrolides/Lincosamides | erythromycin, azithromycin, clarithromycin, clindamycin |
| `glyc` | Glycopeptides | vancomycin, teicoplanin, dalbavancin |
| `tet` | Tetracyclines | tetracycline, doxycycline, minocycline, tigecycline |
| `poly` | Polymyxins | colistin |
| `oxa` | Oxazolidinones | linezolid, tedizolid |
| `chl` | Chloramphenicol | chloramphenicol |
| `sulf` | Sulfonamides | sulfanilamide, trimethoprim-sulfamethoxazole |
| `other` | Other | daptomycin, quinupristin-dalfopristin, nitrofurantoin, fosfomycin, retapamulin, fusidic acid, metronidazole, fidaxomicin, furazolidone, rifampicin |



### A.4 Resistance Mechanisms (35)

See [Section 7.1](#71-resistance-mechanisms) for the full table.



### A.5 Enumerations



#### BacteriaGroup (9 groups)

| Group | Description |
|-------|-------------|
| `Enterobacterales` | Gram-negative enteric rods |
| `NonFermenter` | Non-fermenting Gram-negatives |
| `GramPositive` | Gram-positive cocci and rods |
| `Fastidious` | Fastidious Gram-negatives and atypicals |
| `EntericPathogen` | Specific enteric pathogens (V. cholerae) |
| `Anaerobe` | Obligate anaerobes |
| `Spirochete` | Spirochetes (T. pallidum) |
| `Helicobacter` | Helicobacter/Campylobacter |
| `Mycobacteria` | Mycobacteria (M. tuberculosis) |



#### CarriageCompartment (5)

`Gut`, `Respiratory`, `SkinSoftTissue`, `Genitourinary`, `Systemic`



#### ResistanceAcquisitionType (5)

| Type | Description |
|------|-------------|
| `AtInfection` | Acquired at community infection |
| `AtInfectionHosp` | Acquired at hospital infection |
| `AtInfectionTB` | Acquired at carrier-to-infection (treated-by) conversion |
| `DuringInfection` | De novo emergence during treatment |
| `HGT` | Horizontal gene transfer |



#### InfectionResolutionType (6)

`DrugTreatment`, `NaturalClearance`, `Death`, `SepsisDeath`, `ToxicityDeath`, `BackgroundDeath`



#### ImmunodeficiencyType (2)

`Temporary`, `Chronic`



#### AgeCategory (7)

`Infant`, `Preschool`, `SchoolAge`, `YoungAdult`, `MiddleAge`, `Elderly`, `NotYetBorn`



#### HospitalStatus (2)

`Community`, `Hospital`



#### Region (7)

`NorthAmerica`, `Europe`, `Asia`, `Oceania`, `SouthAmerica`, `Africa`, `Home` (fallback)



#### MicrobiomeResistanceLevel (4)

`None`, `Low`, `Medium`, `High`

---



## Appendix B — Parameter Reference

This appendix lists all parameters defined in the model's configuration. Parameters are stored as key–value pairs in a global HashMap and accessed by string key at runtime.



### B.1 Global Scalars

These are the ~120 top-level parameters stored in the `GlobalScalars` struct:



#### Infection acquisition

| Parameter | Baseline Value | Description |
|----------|---------|-------------|
| `infection_growth_rate_per_day` | 0.1 | Daily bacterial growth increment |
| `infection_initial_level` | 1.0 | Starting bacterial load |
| `infection_clearance_threshold` | 0.5 | Level below which infection resolves |
| `infection_death_threshold` | 50.0 | Level at which death may occur |
| `symptom_onset_threshold` | 3.0 | Level for symptom development |
| `symptom_recheck_interval_days` | 7.0 | Re-evaluation interval |
| `symptom_onset_rate_per_day` | 0.1 | Base symptom development rate |
| `not_under_care_fraction` | 0.05 | Fraction not seeking medical care |



#### Antibiotic treatment initiation

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `antibiotic_initiation_base_log_odds` | −6.5 |
| `antibiotic_initiation_log_odds_symptomatic_infection` | 6.5 |
| `antibiotic_initiation_log_odds_sepsis` | 6.0 |
| `antibiotic_initiation_log_odds_test_identified` | 0.92 |
| `antibiotic_initiation_log_odds_already_on_drug` | 0.18 |
| `antibiotic_initiation_log_odds_immunodeficiency` | 2.08 |
| `antibiotic_initiation_log_odds_no_indication` | −1.05 |



#### Antibiotic treatment cessation

| Parameter | Baseline Value |
|----------|---------|
| `treatment_stop_improvement_threshold` | 2.0 |
| `treatment_stop_rate_per_day` | 0.03 |
| `treatment_duration_base_days` | 7.0 |



#### Drug efficacy

| Parameter | Baseline Value |
|----------|---------|
| `drug_effect_on_bacteria_per_day` | 0.5 |
| `drug_minimum_effective_level` | 0.1 |



#### Drug selection

| Parameter | Baseline Value |
|----------|---------|
| `empiric_therapy_broad_spectrum_bonus` | 0.85 |
| `empiric_therapy_ineffective_drug_penalty` | 0.001 |
| `targeted_therapy_narrow_spectrum_bonus` | 5.0 |
| `targeted_therapy_broad_spectrum_penalty` | 0.1 |
| `targeted_therapy_ineffective_drug_penalty` | 0.001 |
| `regional_resistance_penalty_very_high` | 0.3 |
| `regional_resistance_penalty_high` | 0.5 |
| `regional_resistance_penalty_moderate` | 0.8 |
| `regional_resistance_threshold_very_high` | 0.6 |
| `regional_resistance_threshold_high` | 0.45 |
| `regional_resistance_threshold_moderate` | 0.1 |



#### Testing

| Parameter | Baseline Value |
|----------|---------|
| `bacterial_testing_available_from_day` | 5,478 |
| `resistance_testing_available_from_day` | 9,131 |
| `test_delay_days` | 3.0 |
| `test_rate_per_day` | 0.2 |
| `prob_test_r_done` | 0.95 |
| `test_r_error_probability` | 0.02 |
| `bacterial_testing_base_rate_per_day` | 0.15 |
| `resistance_testing_base_rate_per_day` | 0.95 |
| `bacterial_testing_hospital_multiplier` | 8.0 |
| `resistance_testing_hospital_multiplier` | 5.0 |
| `testing_immunosuppressed_multiplier` | 2.5 |
| `testing_sepsis_multiplier` | 4.0 |
| `bacterial_testing_initial_adoption_rate` | 0.1 |
| `bacterial_testing_max_temporal_multiplier` | 1.0 |
| `resistance_testing_initial_adoption_rate` | 0.05 |
| `resistance_testing_max_temporal_multiplier` | 1.0 |



#### Sepsis onset

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `sepsis_baseline_log_odds` | −14.0 |
| `log_odds_sepsis_infection_level` | 0.8 |
| `log_odds_sepsis_infection_duration` | 0.005 |
| `log_odds_sepsis_onset_immunosuppressed` | 0.7 |
| `log_odds_sepsis_onset_hospitalized` | 0.5 |
| `log_odds_sepsis_onset_not_under_care` | 1.0 |
| `sepsis_age_log_odds_neonatal` | 1.10 |
| `sepsis_age_log_odds_pediatric` | 0.18 |
| `sepsis_age_log_odds_young_adult` | 0.0 |
| `sepsis_age_log_odds_elderly` | 0.69 |



#### Sepsis death

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `sepsis_death_base_log_odds` | −5.0 |
| `sepsis_death_log_odds_age_infant` | 1.1 |
| `sepsis_death_log_odds_age_child` | −0.7 |
| `sepsis_death_log_odds_age_adult` | 0.0 |
| `sepsis_death_log_odds_age_elderly` | 0.9 |
| `sepsis_death_log_odds_immunosuppressed` | 1.5 |
| `sepsis_death_log_odds_bacteria_level` | 0.35 |
| `sepsis_death_log_odds_duration` | 0.04 |
| `sepsis_death_log_odds_early_phase` | 0.8 |
| `sepsis_death_early_phase_days` | 3.0 |
| `sepsis_death_log_odds_not_under_care` | 1.4 |



#### Sepsis recovery

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `sepsis_base_log_odds_of_recovery_per_day` | 0.0 |
| `sepsis_log_odds_bacteria_level` | −0.3 |
| `sepsis_log_odds_in_hospital` | 0.8 |
| `sepsis_log_odds_age_infant` | −0.5 |
| `sepsis_log_odds_age_child` | 0.4 |
| `sepsis_log_odds_age_adult` | 0.0 |
| `sepsis_log_odds_age_elderly` | −0.7 |
| `sepsis_log_odds_immunosuppressed` | −1.0 |
| `sepsis_minimum_duration_days` | 1.0 |



#### Background mortality

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `background_mortality_baseline_log_odds` | −14.0 |
| `log_odds_mortality_per_year_of_age` | 0.04 |
| `log_odds_mortality_per_year_of_age_squared` | 0.05 |
| `mortality_baseline_1930_multiplier` | 3.0 |
| `mortality_baseline_2035_multiplier` | 1.0 |
| `mortality_improvement_half_life_years` | 35.0 |
| `log_odds_mortality_immunosuppressed` | 0.916 |
| `log_odds_mortality_hospitalized` | 0.262 |



#### Hospitalisation

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `hospitalization_base_log_odds` | −10.4 |
| `hospitalization_log_odds_per_age_year` | 0.02 |
| `hospitalization_log_odds_sepsis` | 4.4 |
| `hospitalization_log_odds_symptomatic_infection` | 2.5 |
| `hospitalization_symptomatic_infection_level_threshold` | 3.0 |
| `hospitalization_recovery_rate_per_day` | 0.28 |
| `hospitalization_max_days` | 30.0 |
| `hospitalization_prevent_discharge_with_sepsis` | 1.0 |



#### Immunodeficiency

| Parameter | Baseline Value |
|----------|---------|
| `temporary_immunosuppression_onset_rate_per_day` | 0.00005 |
| `temporary_immunosuppression_recovery_rate_per_day` | 0.01 |
| `chronic_immunosuppression_onset_rate_per_day` | 0.00006 |
| `chronic_immunosuppression_recovery_rate_per_day` | 0.0012 |
| `chronic_immunodeficiency_probability_age_0_1` | 0.3 |
| `chronic_immunodeficiency_probability_age_1_18` | 0.2 |
| `chronic_immunodeficiency_probability_age_18_65` | 0.4 |
| `chronic_immunodeficiency_probability_age_65_plus` | 0.6 |
| `antibiotic_infection_prevention_efficacy` | 0.7 |



#### Resistance

| Parameter | Baseline Value |
|----------|---------|
| `mechanism_assignment_probability` | 0.8 |



#### Microbiome and carriage

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `default_microbiome_clearance_probability_per_day` | 0.01 |
| `microbiome_clearance_probability_on_drug_treatment` | 0.8 |
| `default_microbiome_disruption_log_odds` | 0.3 |
| `microbiome_resistance_multiplier_on_acquisition` | 0.50 |
| `infection_from_microbiome_dampening` | 0.70 |
| `antibiotic_disruption_decay_half_life_days` | 30.0 |
| `carriage_duration_log_odds_coefficient` | −0.01 |
| `carriage_duration_max_log_odds_effect` | −2.0 |
| `antibiotic_clearance_log_odds_per_unit_activity` | 0.5 |
| `carrier_resistance_inheritance_probability` | 0.50 |
| `community_resistance_dilution_factor` | 0.50 |
| `majority_r_window_days` | 100 |
| `majority_r_min_total_samples` | 10 |
| `majority_r_freeze_at_last_positive` | 0.0 |
| `microbiome_majority_promotion_rate_per_day` | 0.02 |



#### Toxicity

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `default_drug_toxicity_death_hazard_per_unit_level` | 0.0 |
| `default_toxicity_reservoir_half_life_days` | 1.5 |
| `toxicity_death_base_log_odds` | −8.0 |
| `toxicity_death_log_odds_per_reservoir_unit` | 2.0 |
| `toxicity_death_log_odds_age_infant` | 0.6 |
| `toxicity_death_log_odds_age_child` | 0.2 |
| `toxicity_death_log_odds_age_adult` | 0.0 |
| `toxicity_death_log_odds_age_elderly` | 0.8 |
| `toxicity_death_log_odds_immunosuppressed` | 0.9 |
| `toxicity_death_log_odds_hospitalized` | 0.25 |
| `toxicity_discontinuation_base_log_odds` | −3.0 |
| `toxicity_discontinuation_log_odds_per_reservoir_unit` | 1.5 |
| `toxicity_discontinuation_log_odds_sepsis` | −1.5 |
| `toxicity_avoidance_penalty_multiplier` | 0.05 |
| `toxicity_avoidance_window_days` | 14.0 |



#### HGT

| Parameter | Baseline Value |
|----------|---------|
| `hgt_base_probability` | 1×10⁻⁵ |
| `hgt_co_infection_multiplier` | 10.0 |
| `hgt_hospital_multiplier` | 5.0 |
| `hgt_microbiome_multiplier` | 2.0 |
| `hgt_gut_compartment_multiplier` | 2.0 |



#### Travel

| Parameter | Baseline Value |
|----------|---------|
| `travel_probability_per_day` | 0.00005 |
| `north_america_travel_multiplier` | 3.0 |
| `europe_travel_multiplier` | 3.5 |
| `oceania_travel_multiplier` | 2.5 |
| `asia_travel_multiplier` | 1.5 |
| `south_america_travel_multiplier` | 0.8 |
| `africa_travel_multiplier` | 0.3 |



### B.2 Per-Bacteria Parameters (42 bacteria × N parameters)

Generated with key pattern `bacteria_{name}_{param}`:

| Parameter suffix | Description |
|------------------|-------------|
| `acquisition_log_odds` | Baseline acquisition rate |
| `hospital_acquisition_log_odds` | Hospital acquisition rate |
| `microbiome_vs_infection_log_odds` | Likelihood of carriage vs infection |
| `clearance_probability_per_day_no_treatment` | Natural clearance rate |
| `growth_rate_multiplier` | Bacteria-specific growth modifier |
| `symptom_onset_override` | Override for symptom onset threshold |
| `hgt_donor_probability` | Probability of being an HGT donor |
| `resistance_floor_enabled` | Whether resistance floors apply |
| `resistance_floor_ramp_years` | Ramp-up period for floors |
| `{drug_class}_resistance_floor` | Per-class floor target |



### B.3 Per-Drug Parameters (58 drugs × N parameters)

Generated with key pattern `drug_{name}_{param}`:

| Parameter suffix | Default | Description |
|------------------|---------|-------------|
| `half_life_days` | Drug-specific | PK half-life |
| `initial_level` | 10.0 | Administration level |
| `double_dose_multiplier` | 2.0 | Double-dose level |
| `spectrum_breadth` | 3.0 | Microbiome disruption breadth |
| `toxicity_death_hazard_per_unit_level` | 0.0 | Per-drug toxicity hazard |
| `toxicity_reservoir_half_life_days` | 1.5 | Per-drug toxicity decay |



### B.4 Per-Region Parameters (6 regions × N parameters)

Generated with key pattern `{region}_{param}`:

| Parameter suffix | Description |
|------------------|-------------|
| `hospitalization_log_odds` | Regional hospitalisation modifier |
| `sepsis_onset_log_odds` | Regional sepsis onset modifier |
| `sepsis_mortality_multiplier` | Regional sepsis death modifier |
| `sepsis_recovery_log_odds` | Regional sepsis recovery modifier |
| `mortality_log_odds` | Regional background mortality modifier |
| `testing_multiplier` | Regional testing rate modifier |
| `antibiotic_initiation_log_odds` | Regional prescribing modifier |
| `drug_{drug}_availability` | Per-drug availability (0.0–1.0) |



### B.5 Per-Mechanism Parameters (35 mechanisms × N parameters)

Generated with key pattern for mechanism-specific parameters:

| Parameter pattern | Description |
|-------------------|-------------|
| `mech_{mechanism}_emergence_rate` | Base emergence rate |
| `mech_{mechanism}_reversion_rate` | Fitness-cost reversion rate |
| `mech_{mechanism}_enhancement_{drug_class}` | Enhancement multiplier (how much resistance the mechanism confers against each drug class) |
| `mech_{mechanism}_global_enhancement_multiplier` | Legacy global multiplier |



### B.6 Cross-Indexed Parameters



#### Potency matrix (42 × 52)

Key: `drug_{drug}_for_bacteria_{bacteria}_potency_when_no_r`

Values: 0.0 (no activity) to 1.0 (maximum activity). See Section 6.5.



#### HGT probability matrix (42 × 42)

Key: `hgt_prob_{source}_to_{target}`

Values: see Section 9.2 for the probability rules.



#### Resistance emergence rates (42 × 35)

Key: `bacteria_{bacteria}_mechanism_{mechanism}_emergence_rate`

Values: base rate × incidence band multiplier. See Section 7.3.



#### Demographic distribution (6 × 18)

Key: `demo_{region}_age_{start}_{end}`

Values: probability weight for each region-age combination.



### B.7 Syndrome Scoring Templates (10 syndromes)

See Section 6.2 for the full empiric scoring tables for all 10 syndromes.



### B.8 Drug Penetration Table (10 syndromes × 36 drug classes)

See Section 6.4 for the penetration values by syndrome and drug class.



### B.9 Cross-Resistance Groups

See Section 7.6 for per-bacteria cross-resistance groupings.



### B.10 Drug Introduction Dates

See Section 6.6 for the full table of 58 drug introduction time steps.



### B.11 Age-Specific Risk Templates

See Section 3.1 for the 7 risk templates and their multiplier vectors.



### B.12 Regional Drug Availability

See Section 6.6 for the per-region availability tables.

---



## Appendix C — Output Specification



### C.1 Output File

Each simulation run produces a single CSV file:

```
amr_simulation_output_analysis_outputs/simulation_summary_NNNNNN.csv
```

where `NNNNNN` is a zero-padded run identifier.



### C.2 Row Structure

Each row represents one simulated day. The number of rows equals the total number of time steps (default 38,325).



### C.3 Column Categories



#### Scalar columns (per-timestep)

| Column | Type | Description |
|--------|------|-------------|
| `day` | int | Simulation day (0-indexed) |
| `year` | float | Calendar year (1930.0 + day/365.25) |
| `total_alive` | int | Living individuals |
| `total_infected` | int | Individuals with ≥1 active infection |
| `total_on_treatment` | int | Individuals receiving ≥1 antibiotic |
| `total_in_hospital` | int | Hospitalised individuals |
| `total_sepsis` | int | Individuals with sepsis |
| `total_died_infection` | int | Cumulative infection deaths |
| `total_died_sepsis` | int | Cumulative sepsis deaths |
| `total_died_background` | int | Cumulative background deaths |
| `total_died_toxicity` | int | Cumulative drug toxicity deaths |
| `total_new_infections` | int | New infections this day |
| `drug_stops_due_to_toxicity` | int | Drug courses stopped for toxicity this day |
| `total_carriers` | int | Individuals carrying ≥1 bacteria |
| `total_immunosuppressed` | int | Immunosuppressed individuals |
| `policy_name` | string | Active policy branch name |



#### Per-bacteria columns (~42 each)

| Pattern | Description |
|---------|-------------|
| `{bacteria}_infected` | Currently infected count |
| `{bacteria}_carriers` | Current carrier count |
| `{bacteria}_deaths` | Cumulative deaths from this species |
| `{bacteria}_new_infections` | New infections this day |
| `{bacteria}_new_infections_community` | From community acquisition |
| `{bacteria}_new_infections_hospital` | From hospital acquisition |
| `{bacteria}_new_infections_carrier` | From carrier-to-infection |
| `{bacteria}_sepsis` | Current sepsis count |



#### Per-drug columns (~58 each)

| Pattern | Description |
|---------|-------------|
| `{drug}_prescribed` | Courses initiated this day |
| `{drug}_active_treatments` | Currently on this drug |



#### Per-bacteria × per-drug columns (~2,436 each)

| Pattern | Description |
|---------|-------------|
| `{bacteria}_{drug}_activity_r` | Mean resistance (activity_r) |
| `{bacteria}_{drug}_majority_r` | Population-level resistance prevalence |



#### Per-region columns (~6 each)

| Pattern | Description |
|---------|-------------|
| `{region}_infected` | Regional infection count |
| `{region}_hospitalized` | Regional hospital count |
| `{region}_deaths` | Regional death count |



### C.4 Total Column Count

With 42 bacteria, 58 drugs, and 6 regions, the CSV contains approximately:

- ~16 scalar columns
- ~336 per-bacteria columns (42 × 8)
- ~116 per-drug columns (58 × 2)
- ~4,872 per-bacteria-per-drug columns (42 × 58 × 2)
- ~18 per-region columns (6 × 3)
- **Total: ~5,358 columns**



### C.5 Infection Journey Logs

When enabled, individual infection journeys are logged to the `infection_journeys/` directory as CSV files, capturing:

- Infection acquisition details
- Resistance profile at acquisition and over time
- Treatment episodes
- Clinical outcome (clearance, death, ongoing)
- Mechanism gains and losses

---

*This document describes the model as implemented in the Rust codebase. All variable names correspond to parameter keys used in `src/config.rs`.*
