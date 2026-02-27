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
- [Appendix A — Bacteria, Drugs, Mechanisms and Enums](#appendix-a--bacteria-drugs-mechanisms-and-enums)
- [Appendix B — Parameter Reference](#appendix-b--parameter-reference)
- [Appendix C — Output Specification](#appendix-c--output-specification)

---

## 1. Overview

This model simulates the emergence, transmission, and dynamics of antimicrobial resistance (AMR) across a synthetic human population from 1930 to 2035. It is an individual-based (agent-based) model in which each person can acquire bacterial infections, receive antibiotic treatment, develop resistance through de novo mutation or horizontal gene transfer, and carry resistant organisms in their microbiome.

The model tracks **42 bacterial species**, **58 antibiotics** (grouped into **18 drug classes**), and **35 resistance mechanisms**. The population is distributed across **6 world regions** (North America, Europe, Asia, Oceania, South America, Africa), each with distinct epidemiological and healthcare profiles.

Time advances in discrete **daily steps**. On each day, every individual in the population is processed through a sequence of 21 rules covering ageing, infection acquisition, clinical progression, treatment, resistance dynamics, and death.

### Scope

The model is designed for:

- Estimating the global burden of AMR over time
- Evaluating the potential impact of antibiotic stewardship policies
- Exploring "what if" counterfactual scenarios (e.g., a world without resistance)
- Understanding how resistance mechanisms spread across bacterial species and regions

### Limitations

- Drug levels are modelled as abstract potency units rather than pharmacokinetic concentrations
- Bacteria-bacteria competition within the microbiome is represented implicitly through resistance promotion and decay rather than explicit strain dynamics
- The model does not capture within-host spatial heterogeneity (e.g., biofilm vs planktonic)
- Vaccine effects are not explicitly modelled (though their population-level impact is partially captured through incidence parameters)
- Region definitions are broad continental groupings

---

## 2. Population and Demographics

### 2.1 Initialisation

The population is created at day 0 (calendar year 1930). Each individual is assigned:

- **Age**: drawn from a demographic distribution that encodes both living individuals and future births (negative age values represent individuals not yet born)
- **Sex**: male or female with equal probability
- **Region**: sampled from demographic weights reflecting the global population distribution

Age bands are specified in 4,000-day (~11-year) intervals. The demographic distribution is the product of per-region weights:

| Variable name | Description |
|---------------|-------------|
| `demo_{region}_age_{start}_{end}` | Probability weight for an individual in a given region and age band |

The six regions and their approximate population shares:

| Region | Share |
|--------|-------|
| Asia | 55% |
| Europe | 15% |
| Africa | 12% |
| North America | 9% |
| South America | 6% |
| Oceania | 3% |

### 2.2 Ageing

Each day, every individual's age increments by one day. Age categories are recalculated:

| Age category | Age range | Variable suffix |
|--------------|-----------|-----------------|
| Infant | 0–1 year | `infant` |
| Preschool | 1–5 years | `preschool` |
| School age | 5–18 years | `school` |
| Young adult | 18–50 years | `young_adult` |
| Middle age | 50–70 years | `middle_age` |
| Elderly | 70+ years | `elderly` |

An additional classification is used for sepsis and mortality:

| Category | Age range | Variable suffix |
|----------|-----------|-----------------|
| Neonatal | 0–28 days | `neonatal` |
| Pediatric | 28 days–18 years | `pediatric` |
| Young adult | 18–50 years | `young_adult` |
| Elderly | 50+ years | `elderly` |

### 2.3 Immunodeficiency

Individuals can be immunosuppressed, increasing their risk of infection acquisition, sepsis, and death. Two types exist:

| Type | Onset rate | Recovery rate |
|------|-----------|---------------|
| Temporary | `temporary_immunosuppression_onset_rate_per_day` = 0.00005 | `temporary_immunosuppression_recovery_rate_per_day` = 0.01 |
| Chronic | `chronic_immunosuppression_onset_rate_per_day` = 0.00006 | `chronic_immunosuppression_recovery_rate_per_day` = 0.0012 |

Chronic immunodeficiency probability at birth varies by age category:

| Variable | Value |
|----------|-------|
| `chronic_immunodeficiency_probability_age_0_1` | 0.3 |
| `chronic_immunodeficiency_probability_age_1_18` | 0.2 |
| `chronic_immunodeficiency_probability_age_18_65` | 0.4 |
| `chronic_immunodeficiency_probability_age_65_plus` | 0.6 |

### 2.4 Hospitalisation

Hospital admission is modelled as a logistic function of clinical state:

$$P(\text{admission}) = \frac{1}{1 + e^{-\text{log\_odds}}}$$

where:

$$\text{log\_odds} = \text{base} + \text{age\_effect} + \text{sepsis\_effect} + \text{infection\_effect} + \text{region\_effect}$$

| Variable | Default | Description |
|----------|---------|-------------|
| `hospitalization_base_log_odds` | −10.4 | Baseline (≈0.003%) |
| `hospitalization_log_odds_per_age_year` | 0.02 | ~2% increase per year of age |
| `hospitalization_log_odds_sepsis` | 4.4 | ln(80), strong driver |
| `hospitalization_log_odds_symptomatic_infection` | 2.5 | ~12× for symptomatic cases |
| `hospitalization_symptomatic_infection_level_threshold` | 3.0 | Minimum infection level |
| `hospitalization_recovery_rate_per_day` | 0.28 | ~3.6-day average stay |
| `hospitalization_max_days` | 30.0 | Maximum stay |
| `hospitalization_prevent_discharge_with_sepsis` | 1.0 | Block discharge during sepsis |

Regional hospital admission log-odds:

| Region | Variable | Value |
|--------|----------|-------|
| North America | `north_america_hospitalization_log_odds` | 0.5 |
| Europe | `europe_hospitalization_log_odds` | 0.6 |
| Oceania | `oceania_hospitalization_log_odds` | 0.4 |
| Asia | `asia_hospitalization_log_odds` | 0.0 |
| South America | `south_america_hospitalization_log_odds` | −0.2 |
| Africa | `africa_hospitalization_log_odds` | −0.5 |

### 2.5 Travel

Individuals may travel between regions, acquiring different drug availability and resistance profiles:

| Variable | Default | Description |
|----------|---------|-------------|
| `travel_probability_per_day` | 0.00005 | Base daily travel probability |
| `north_america_travel_multiplier` | 3.0 | Relative travel frequency |
| `europe_travel_multiplier` | 3.5 | |
| `oceania_travel_multiplier` | 2.5 | |
| `asia_travel_multiplier` | 1.5 | |
| `south_america_travel_multiplier` | 0.8 | |
| `africa_travel_multiplier` | 0.3 | |

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
- Population-level resistance prevalence (`majority_r`) reducing the probability that the bacteria the individual acquires is resistant

| Variable pattern | Description |
|------------------|-------------|
| `bacteria_{name}_acquisition_log_odds` | Baseline acquisition log-odds for each species |
| `{region}_bacteria_{name}_acquisition_log_odds` | Regional override |
| `bacteria_{name}_log_odds_{age_category}` | Age-specific override |
| `{bacteria}_{region}_log_odds_{age_category}` | Bacteria × region × age interaction |

#### Age risk templates

Each bacteria is assigned a risk template that defines its age distribution:

| Template | Multipliers [0–1y, 1–5y, 5–18y, 18–50y, 50–70y, 70+y] |
|----------|--------------------------------------------------------|
| `respiratory` | [3.0, 1.8, 0.8, 1.0, 1.3, 2.5] |
| `gastrointestinal` | [2.5, 2.0, 1.2, 1.0, 1.1, 1.8] |
| `urogenital` | [1.2, 0.8, 0.9, 1.0, 1.4, 2.2] |
| `skin_soft_tissue` | [1.5, 1.3, 1.1, 1.0, 1.2, 1.8] |
| `bloodstream` | [4.0, 2.0, 0.7, 1.0, 1.5, 3.0] |
| `sexually_transmitted` | [0.1, 0.2, 0.8, 1.0, 0.8, 0.3] |
| `flat` | [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] |

Template assignments:

| Variable | Default | Selected overrides |
|----------|---------|-------------------|
| `{bacteria}_age_risk_template` | `"respiratory"` | |
| `salm_typhi_age_risk_template` | `"gastrointestinal"` | |
| `esch_coli_age_risk_template` | `"urogenital"` | |
| `pseud_aerug_age_risk_template` | `"bloodstream"` | |
| `staph_aureus_age_risk_template` | `"skin_soft_tissue"` | |
| `n_gonorrhoeae_age_risk_template` | `"sexually_transmitted"` | |

#### Age-specific overrides (selected examples)

**STIs** — C. trachomatis, N. gonorrhoeae, T. pallidum: peak in young adults (log-odds 0.4–2.0), very low in children (−2.5 to −4.3).

**Respiratory** — S. pneumoniae, H. influenzae, M. catarrhalis: U-shaped, high in infants (1.7–2.5) and elderly (1.0–1.2), low in young adults (−0.4 to −0.8).

**Enteric** — S. Typhi, Shigella, C. jejuni: highest in preschool/infants (1.2–2.0). H. pylori accumulates with age (elderly 1.8).

**Healthcare-associated** — A. baumannii peaks in elderly (1.5), C. difficile extreme elderly peak (2.0).

**Invasive** — N. meningitidis bimodal: infant (1.8) + young adult (1.3). L. monocytogenes: neonatal (1.8) + elderly (1.5).

### 3.2 Hospital Acquisition

Hospitalised individuals have elevated acquisition rates controlled by:

| Variable pattern | Description |
|------------------|-------------|
| `bacteria_{name}_hospital_acquisition_log_odds` | Hospital-specific baseline |

### 3.3 Carrier-Derived Infection

Individuals carrying bacteria in their microbiome can develop active infection from their carriage flora. This pathway is probability-gated:

| Variable | Default | Description |
|----------|---------|-------------|
| `carrier_resistance_inheritance_probability` | 0.50 | Probability the new infection inherits resistance from microbiome |
| `infection_from_microbiome_dampening` | 0.70 | Dampening factor applied to carrier-to-infection conversion |

### 3.4 Resistance at Acquisition

When a new infection is acquired from the community, its initial resistance profile is derived from:

1. **Population-level prevalence** (`majority_r`): the rolling average of resistance observed across all individuals for that bacteria–drug combination
2. **Community resistance dilution**: controlled by `community_resistance_dilution_factor` (default 0.50), which scales down the community resistance signal to avoid over-estimating resistance in newly acquired infections
3. **Mechanism profile sampling**: resistance mechanisms are not assigned independently — instead, a correlated mechanism profile is sampled from a reservoir of up to 200 profiles observed in currently infected individuals (`MechanismProfileCache`), preserving realistic co-resistance patterns

---

## 4. Clinical Progression

### 4.1 Syndrome Assignment

Each new infection is assigned a clinical syndrome (site of infection), which affects drug penetration, treatment selection, and progression:

| Syndrome | Index | Examples |
|----------|-------|---------|
| UTI | 1 | Urinary tract infection |
| Skin/soft tissue | 2 | Cellulitis, wound infection |
| Respiratory | 3 | Pneumonia, bronchitis |
| Bloodstream | 4 | Bacteraemia |
| Intra-abdominal | 5 | Peritonitis, abscess |
| CNS | 6 | Meningitis |
| Gastrointestinal | 7 | Gastroenteritis |
| Genital/pelvic | 8 | PID, STI |
| Bone/joint | 9 | Osteomyelitis |
| Other | 10 | Device-related, other |

Syndrome probabilities are bacteria-specific. For example, E. coli infections are predominantly UTIs and bloodstream, while S. pneumoniae infections are predominantly respiratory.

#### Syndrome initiation multipliers

| Variable | Syndrome | Value | Effect |
|----------|----------|-------|--------|
| `syndrome_3_initiation_multiplier` | Respiratory | 2.5 | Higher treatment-seeking |
| `syndrome_7_initiation_multiplier` | GI | 0.5 | Lower treatment-seeking |
| `syndrome_8_initiation_multiplier` | Genital | 0.25 | Lower treatment-seeking |
| Others | — | 1.0 | Default |

#### Syndrome bacteria growth multipliers

| Variable | Syndrome | Value |
|----------|----------|-------|
| `syndrome_1_bacteria_growth_multiplier` | UTI | 0.7 |
| `syndrome_2_bacteria_growth_multiplier` | Skin | 0.8 |
| `syndrome_3_bacteria_growth_multiplier` | Respiratory | 1.0 |
| `syndrome_4_bacteria_growth_multiplier` | Bloodstream | 1.3 |
| `syndrome_5_bacteria_growth_multiplier` | Intra-abdominal | 1.1 |
| `syndrome_6_bacteria_growth_multiplier` | CNS | 1.2 |
| `syndrome_7_bacteria_growth_multiplier` | GI | 0.9 |
| `syndrome_8_bacteria_growth_multiplier` | Genital | 0.6 |
| `syndrome_9_bacteria_growth_multiplier` | Bone/joint | 0.8 |
| `syndrome_10_bacteria_growth_multiplier` | Other | 1.0 |

### 4.2 Symptom Onset

Infection progression from asymptomatic to symptomatic is driven by bacterial load exceeding a threshold, modified by:

- Bacteria-specific symptom onset rates
- Immune status
- Duration of infection

### 4.3 Sepsis

Sepsis onset is a logistic function of multiple risk factors:

$$\text{log\_odds}_{\text{sepsis}} = \text{baseline} + \text{infection\_level} + \text{duration} + \text{age} + \text{region} + \text{syndrome} + \text{immune} + \text{hospital} + \text{care\_status}$$

| Variable | Default | Description |
|----------|---------|-------------|
| `sepsis_baseline_log_odds` | −14.0 | Very low baseline |
| `log_odds_sepsis_infection_level` | 0.8 | Per unit of bacterial load |
| `log_odds_sepsis_infection_duration` | 0.005 | Per day of infection |
| `log_odds_sepsis_onset_immunosuppressed` | 0.7 | ~2× for immunosuppressed |
| `log_odds_sepsis_onset_hospitalized` | 0.5 | ~1.6× for hospitalised |
| `log_odds_sepsis_onset_not_under_care` | 1.0 | ~2.7× if not under care |

#### Bacteria-specific sepsis baseline log-odds

| Bacteria | Value | Comment |
|----------|-------|---------|
| E. coli | −21.0 | Very common, rarely septic |
| N. gonorrhoeae | −21.0 | Extremely rare sepsis |
| C. trachomatis | −19.0 | |
| C. jejuni | −20.0 | |
| Shigella spp. | −12.0 | |
| C. difficile | −11.0 | |
| M. catarrhalis | −11.0 | |
| T. pallidum | −12.0 | |
| B. pertussis | −11.0 | |
| S. pneumoniae | −9.0 | |
| H. influenzae | −9.0 | |
| V. cholerae | −9.0 | |
| N. meningitidis | −9.0 | |
| H. pylori | −250.0 | Effectively never septic |
| MDR M. tuberculosis | −38.0 | Very rare sepsis |
| Most other bacteria | −8.0 | Default |

#### Age-dependent sepsis log-odds

| Variable | Value | Effect |
|----------|-------|--------|
| `sepsis_age_log_odds_neonatal` | 1.10 | ~3× (ln 3.0) |
| `sepsis_age_log_odds_pediatric` | 0.18 | ~1.2× (ln 1.2) |
| `sepsis_age_log_odds_young_adult` | 0.0 | Reference |
| `sepsis_age_log_odds_elderly` | 0.69 | ~2× (ln 2.0) |

#### Per-bacteria-age sepsis interactions

**Neonatal**: S. agalactiae +0.981, E. coli +0.511, L. monocytogenes +0.693, S. aureus +0.288.

**Pediatric**: S. pneumoniae +0.916, H. influenzae +0.734, N. meningitidis +0.511, S. aureus +0.511.

**Elderly**: S. pneumoniae +0.693, E. coli +0.223, K. pneumoniae +0.405, P. aeruginosa +0.560, A. baumannii +0.405, E. faecium +0.336, S. aureus +0.470.

**Young adult**: N. meningitidis +0.336, S. aureus +0.588.

#### Regional sepsis onset log-odds

| Region | Variable | Value |
|--------|----------|-------|
| North America | `north_america_sepsis_onset_log_odds` | −0.5 |
| Europe | `europe_sepsis_onset_log_odds` | −0.6 |
| Oceania | `oceania_sepsis_onset_log_odds` | −0.5 |
| Asia | `asia_sepsis_onset_log_odds` | −0.1 |
| South America | `south_america_sepsis_onset_log_odds` | 0.0 |
| Africa | `africa_sepsis_onset_log_odds` | 0.1 |

#### Syndrome-specific sepsis log-odds

| Syndrome | Value | Comment |
|----------|-------|---------|
| UTI (1) | −2.0 | Low sepsis risk |
| Skin (2) | −1.0 | Low |
| Respiratory (3) | 0.0 | Moderate |
| Bloodstream (4) | 1.5 | High |
| Intra-abdominal (5) | 0.8 | Moderate-high |
| CNS (6) | 1.2 | High |
| GI (7) | −0.5 | Low-moderate |
| Genital (8) | −1.5 | Low |
| Bone/joint (9) | 0.5 | Moderate |

### 4.4 Natural Clearance

Infections can resolve naturally through immune-mediated clearance, subject to:

- Bacteria-specific clearance probabilities
- Antibiotic activity (enhances clearance)
- Duration of infection (longer infections are harder to clear)
- Immune status

| Variable | Default | Description |
|----------|---------|-------------|
| `default_microbiome_clearance_probability_per_day` | 0.01 | Base clearance rate |
| `microbiome_clearance_probability_on_drug_treatment` | 0.8 | Enhanced clearance on antibiotics |
| `antibiotic_clearance_log_odds_per_unit_activity` | 0.5 | Clearance boost per unit of active drug |
| `carriage_duration_log_odds_coefficient` | −0.01 | Longer carriage → harder to clear |
| `carriage_duration_max_log_odds_effect` | −2.0 | Maximum reduction (~7.4× harder) |

---

## 5. Diagnostic Testing

### 5.1 Testing Availability

Bacterial identification and resistance testing become available at historical time points:

| Variable | Default | Calendar year |
|----------|---------|---------------|
| `bacterial_testing_available_from_day` | 5,478 | ~1945 |
| `resistance_testing_available_from_day` | 9,131 | ~1955 |

### 5.2 Legacy Testing Parameters

| Variable | Default | Description |
|----------|---------|-------------|
| `test_delay_days` | 3.0 | Days between sample and result |
| `test_rate_per_day` | 0.2 | Base probability of ordering a test |
| `prob_test_r_done` | 0.95 | Probability AST is performed given culture |
| `test_r_error_probability` | 0.02 | False result rate |

### 5.3 Modern Testing Model

| Variable | Default | Description |
|----------|---------|-------------|
| `bacterial_testing_base_rate_per_day` | 0.15 | Base culture rate |
| `resistance_testing_base_rate_per_day` | 0.95 | Conditional on positive culture |
| `bacterial_testing_hospital_multiplier` | 8.0 | Hospital increases testing |
| `resistance_testing_hospital_multiplier` | 5.0 | |
| `testing_immunosuppressed_multiplier` | 2.5 | |
| `testing_sepsis_multiplier` | 4.0 | |

#### Regional testing multipliers

| Region | Value |
|--------|-------|
| North America | 1.1 |
| Europe | 1.2 |
| Asia | 0.7 |
| South America | 0.6 |
| Oceania | 0.8 |
| Africa | 0.3 |

#### Temporal adoption (S-curve)

| Variable | Default |
|----------|---------|
| `bacterial_testing_initial_adoption_rate` | 0.1 |
| `bacterial_testing_max_temporal_multiplier` | 1.0 |
| `resistance_testing_initial_adoption_rate` | 0.05 |
| `resistance_testing_max_temporal_multiplier` | 1.0 |

---

## 6. Antibiotic Treatment

### 6.1 Treatment Initiation

The decision to start antibiotics is modelled as a logistic function:

$$P(\text{initiation}) = \frac{1}{1 + e^{-\text{log\_odds}}}$$

| Variable | Default | Description |
|----------|---------|-------------|
| `antibiotic_initiation_base_log_odds` | −6.5 | Baseline (~0.1%) |
| `antibiotic_initiation_log_odds_symptomatic_infection` | 6.5 | Lifts to ~26% |
| `antibiotic_initiation_log_odds_sepsis` | 6.0 | Very strong driver |
| `antibiotic_initiation_log_odds_test_identified` | 0.92 | ln(2.5) |
| `antibiotic_initiation_log_odds_already_on_drug` | 0.18 | ln(1.2) |
| `antibiotic_initiation_log_odds_immunodeficiency` | 2.08 | ln(8.0) |
| `antibiotic_initiation_log_odds_no_indication` | −1.05 | ln(0.35), dampens inappropriate use |

#### Regional antibiotic initiation log-odds

| Region | Value | Approximate effect |
|--------|-------|--------------------|
| North America | 0.0 | Reference |
| Europe | 0.0 | Same |
| Oceania | 0.0 | Same |
| Asia | −0.5 | ~38% reduction |
| South America | −0.8 | ~55% reduction |
| Africa | −1.4 | ~75% reduction |

### 6.2 Drug Selection

Drug selection follows a two-stage scoring system:

1. **Empiric therapy** (no test result available): drugs are scored based on a syndrome-specific template
2. **Targeted therapy** (test result available): drugs are scored based on known susceptibility

#### Empiric scoring parameters

| Variable | Default | Description |
|----------|---------|-------------|
| `empiric_therapy_broad_spectrum_bonus` | 0.85 | Bonus for broader coverage |
| `empiric_therapy_ineffective_drug_penalty` | 0.001 | Penalty for likely-ineffective drugs |

#### Targeted scoring parameters

| Variable | Default | Description |
|----------|---------|-------------|
| `targeted_therapy_narrow_spectrum_bonus` | 5.0 | Strong preference for narrow-spectrum |
| `targeted_therapy_broad_spectrum_penalty` | 0.1 | Penalty for broad-spectrum when narrow available |
| `targeted_therapy_ineffective_drug_penalty` | 0.001 | Penalty for ineffective drugs |

#### Regional resistance surveillance penalties

Drugs with high population-level resistance in a region receive selection penalties:

| Variable | Default | Description |
|----------|---------|-------------|
| `regional_resistance_penalty_very_high` | 0.3 | Score multiplier when resistance >60% |
| `regional_resistance_penalty_high` | 0.5 | When resistance >45% |
| `regional_resistance_penalty_moderate` | 0.8 | When resistance >10% |
| `regional_resistance_threshold_very_high` | 0.6 | |
| `regional_resistance_threshold_high` | 0.45 | |
| `regional_resistance_threshold_moderate` | 0.1 | |

#### Syndrome-specific empiric scoring templates

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

| Variable pattern | Default | Description |
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

Each antibiotic becomes available at a specific time step:

| Drug | Time step | ~Year |
|------|-----------|-------|
| sulfanilamide | 2,555 | 1937 |
| penicillin_g | 3,555 | 1942 |
| chloramphenicol | 6,935 | 1949 |
| tetracycline | 6,575 | 1948 |
| colistin | 8,020 | 1952 |
| erythromycin | 8,025 | 1952 |
| nitrofurantoin | 8,395 | 1953 |
| furazolidone | 9,125 | 1955 |
| vancomycin | 10,215 | 1958 |
| fosfomycin | 10,590 | 1959 |
| metronidazole | 10,965 | 1960 |
| ampicillin | 11,315 | 1961 |
| fusidic_a | 11,680 | 1962 |
| gentamicin | 12,045 | 1963 |
| rifampicin | 13,140 | 1966 |
| doxycycline | 13,505 | 1967 |
| clindamycin | 13,870 | 1968 |
| trim_sulf | 13,870 | 1968 |
| amoxicillin | 13,780 | 1972 |
| cephalexin | 14,605 | 1970 |
| minocycline | 14,965 | 1971 |
| cefazolin | 15,700 | 1973 |
| tobramycin | 16,325 | 1975 |
| amikacin | 16,690 | 1976 |
| ticarcillin | 14,600 | 1977 |
| cefuroxime | 17,525 | 1978 |
| piperacillin | 16,065 | 1981 |
| ceftriaxone | 19,715 | 1984 |
| piperacillin_tazobactam | 19,715 | 1984 |
| ceftazidime | 20,080 | 1985 |
| imipenem_c | 20,080 | 1985 |
| amoxicillin_clavulanate | 16,425 | 1985 |
| aztreonam | 20,445 | 1986 |
| ciprofloxacin | 20,805 | 1987 |
| teicoplanin | 21,170 | 1988 |
| ampicillin_sulbactam | 18,250 | 1990 |
| ticarcillin_clavulanate | 18,250 | 1990 |
| clarithromycin | 21,895 | 1990 |
| ofloxacin | 21,895 | 1990 |
| azithromycin | 22,260 | 1991 |
| cefepime | 24,195 | 1996 |
| meropenem | 24,195 | 1996 |
| levofloxacin | 24,195 | 1996 |
| moxifloxacin | 25,290 | 1999 |
| quinu_dalfo | 25,290 | 1999 |
| linezolid | 25,550 | 2000 |
| ertapenem | 25,920 | 2001 |
| daptomycin | 27,375 | 2005 |
| ceftazidime_avibactam | 27,740 | 2006 |
| tigecycline | 28,040 | 2007 |
| retapamulin | 28,405 | 2007 |
| ceftaroline | 29,305 | 2010 |
| fidaxomicin | 29,565 | 2011 |
| tedizolid | 30,660 | 2014 |
| dalbavancin | 30,660 | 2014 |
| ceftolozane_tazobactam | 30,295 | 2014 |
| meropenem_vaborbactam | 32,045 | 2018 |
| cefiderocol | 33,510 | 2019 |

**Special case — Colistin**: Withdrawn from routine use between ~1970 and ~1995 (availability drops to 5% during that window), reflecting historical concerns about nephrotoxicity before its re-adoption for MDR Gram-negative infections.

### 6.7 Drug Toxicity

Drugs accumulate a toxicity reservoir, which can trigger sub-lethal discontinuation or lethal toxicity death:

| Variable | Default | Description |
|----------|---------|-------------|
| `default_drug_toxicity_death_hazard_per_unit_level` | 0.0 | Per-drug hazard contribution (most drugs = 0; high for colistin, aminoglycosides) |
| `default_toxicity_reservoir_half_life_days` | 1.5 | Decay rate of accumulated toxicity |

#### Toxicity discontinuation (sub-lethal)

When accumulated toxicity exceeds a threshold, the treating clinician may discontinue the drug:

| Variable | Default | Description |
|----------|---------|-------------|
| `toxicity_discontinuation_base_log_odds` | −3.0 | Baseline discontinuation probability |
| `toxicity_discontinuation_log_odds_per_reservoir_unit` | 1.5 | Per unit of toxicity reservoir |
| `toxicity_discontinuation_log_odds_sepsis` | −1.5 | Clinicians tolerate more toxicity during sepsis |
| `toxicity_avoidance_penalty_multiplier` | 0.05 | Recently-stopped drugs are avoided in reselection |
| `toxicity_avoidance_window_days` | 14.0 | Duration of avoidance penalty |

#### Lethal toxicity death

| Variable | Default | Description |
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

| Variable | Default |
|----------|---------|
| `antibiotic_infection_prevention_efficacy` | 0.7 |

---

## 7. Resistance Dynamics

### 7.1 Resistance Mechanisms

The model tracks 35 distinct resistance mechanisms, each representing a biological pathway for antibiotic resistance:

| Mechanism | Variable name | Description |
|-----------|--------------|-------------|
| ESBL CTX-M | `esbl_ctx_m` | Extended-spectrum β-lactamase (most common ESBL globally) |
| ESBL TEM | `esbl_tem` | Extended-spectrum β-lactamase (TEM-type) |
| ESBL SHV | `esbl_shv` | Extended-spectrum β-lactamase (SHV-type) |
| AmpC CMY | `ampc_cmy` | Plasmid-mediated AmpC β-lactamase |
| AmpC DHA | `ampc_dha` | Plasmid-mediated AmpC β-lactamase |
| KPC | `kpc` | Klebsiella pneumoniae carbapenemase |
| NDM/VIM | `ndm_vim` | Metallo-β-lactamases |
| OXA-48 | `oxa_48` | Oxacillinase-type carbapenemase |
| PBP2a/MecA | `pbp2a_meca` | Penicillin-binding protein alteration (MRSA) |
| VanA | `vana` | High-level vancomycin resistance |
| VanB | `vanb` | Variable-level vancomycin resistance |
| GyrA (primary) | `gyra_primary` | DNA gyrase mutation (fluoroquinolone resistance, step 1) |
| GyrA + ParC | `gyra_parc` | Additional topoisomerase mutation (high-level FQ resistance) |
| ErmB | `ermb` | Erythromycin ribosome methylase (macrolide resistance) |
| Cfr | `cfr` | 23S rRNA methyltransferase (multi-class resistance) |
| 16S rRMT | `16s_rrmt` | 16S rRNA methyltransferase (aminoglycoside resistance) |
| CAT | `cat` | Chloramphenicol acetyltransferase |
| Qnr | `qnr` | Quinolone resistance protein |
| MCR-1 | `mcr_1` | Mobilised colistin resistance |
| AcrAB-TolC | `acrab_tolc` | Gram-negative efflux pump |
| MexXY-OprM | `mexxy_oprm` | Pseudomonas-specific efflux pump |
| OmpK35/36 | `ompk35_36` | Outer membrane porin loss (Klebsiella) |
| OprD | `oprd` | Outer membrane porin loss (Pseudomonas) |
| Global efflux | `global_efflux` | Non-specific efflux upregulation |
| Global porin loss | `global_porin_loss` | Non-specific porin downregulation |
| Folate pathway | `folate_pathway` | Altered dihydrofolate reductase (trimethoprim resistance) |
| Nitroreductase | `nitroreductase` | Nitroreductase loss (nitrofurantoin resistance) |
| FosA | `fosa` | Fosfomycin-modifying enzyme |
| MprF | `mprf` | Membrane charge modification (daptomycin resistance) |
| RpoB | `rpob` | RNA polymerase mutation (rifampicin resistance) |
| FusB | `fusb` | Fusidic acid resistance determinant |
| As-yet-unknown 1–4 | `as_yet_unknown_{1..4}` | Placeholder for future/novel mechanisms |

### 7.2 Mechanism–Drug-Class Enhancement Multipliers

Each mechanism reduces the efficacy of specific drug classes. The enhancement multiplier (0.0–1.0) determines how much a mechanism reduces drug activity: 0.0 = no effect, 1.0 = complete resistance.

These are defined per mechanism × drug class (35 × 18 = 630 values), with variable pattern:

```
mech_{mechanism}_enhancement_{drug_class}
```

**Global (legacy) enhancement multipliers** (used when no per-class override exists):

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

De novo resistance emergence occurs when an individual is under antibiotic pressure. The probability is a function of:

- Base emergence rate for the mechanism
- Bacteria-specific incidence band multiplier
- Whether the bacteria is biologically capable of acquiring that mechanism

#### Incidence band multipliers

Bacteria are classified into incidence bands reflecting their frequency in the population. More common bacteria have lower per-capita emergence rates to avoid unrealistically fast resistance evolution:

| Band | Multiplier | Example bacteria |
|------|-----------|-----------------|
| High incidence | ×0.1 | E. coli, Shigella, N. gonorrhoeae, C. jejuni, H. pylori, S. aureus, S. pneumoniae, M. pneumoniae, C. trachomatis |
| Moderate incidence | ×1.0 | K. pneumoniae, S. Typhi, S. pyogenes, B. pertussis, M. genitalium, T. pallidum |
| Low incidence | ×3.0 | Enterobacter, E. cloacae, Proteus, iNTS, P. aeruginosa, V. cholerae, M. catarrhalis, H. influenzae, S. epidermidis, S. agalactiae, E. faecalis, C. difficile, B. fragilis |
| Very low incidence | ×10.0 | Citrobacter, Morganella, Serratia, P. stuartii, A. baumannii, S. maltophilia, B. cepacia, N. meningitidis, L. pneumophila, E. faecium, L. monocytogenes, MDR TB, Y. enterocolitica |

Base emergence rates are mechanism-specific, ranging from ~1×10⁻⁹ (rare chromosomal mutations) to ~1×10⁻⁶ (common plasmid-mediated acquisitions).

**Mechanism assignment probability**: `mechanism_assignment_probability` = 0.8 — when a resistance event occurs, this is the probability the specific mechanism responsible is tracked (versus being attributed to a generic mechanism).

### 7.4 Resistance Reversion

Resistance mechanisms that impose a fitness cost on the bacteria can be lost over time when antibiotic pressure is removed. Reversion is:

- **Drug-gated**: reversion only occurs when the individual is NOT currently exposed to any drug in the relevant drug class
- **Mechanism-specific**: each mechanism has its own reversion rate

| Mechanism | Reversion rate (per day) |
|-----------|-------------------------|
| ESBL CTX-M, TEM, SHV | 0.0005 |
| AmpC CMY, DHA | 0.0003 |
| KPC, NDM/VIM | 0.0002 |
| OXA-48 | 0.0004 |
| PBP2a/MecA | 0.0003 |
| VanA | 0.002 |
| VanB | 0.001 |
| GyrA primary, GyrA + ParC | 0.0001 |
| ErmB | 0.0008 |
| Cfr | 0.001 |
| 16S rRMT | 0.0005 |
| CAT | 0.001 |
| Qnr | 0.001 |
| MCR-1 | 0.001 |
| AcrAB-TolC, MexXY-OprM | 0.0005 |
| OmpK35/36, OprD | 0.0003 |
| Global efflux, Global porin loss | 0.0005 |
| Folate pathway | 0.0008 |
| Nitroreductase | 0.001 |
| FosA | 0.001 |
| MprF | 0.0008 |
| RpoB | 0.002 |
| FusB | 0.001 |
| As-yet-unknown 1–4 | 0.0005 |

### 7.5 Resistance Floors

For clinically important bacteria, the model enforces minimum resistance levels that ramp up after the introduction of each drug class. This prevents unrealistically low resistance for well-established resistance phenotypes:

| Variable pattern | Default | Description |
|------------------|---------|-------------|
| `resistance_floor_feature_enabled` | >0.5 to enable | Global toggle |
| `bacteria_{name}_resistance_floor_enabled` | >0.5 to enable | Per-bacteria toggle |
| `bacteria_{name}_resistance_floor_ramp_years` | 10.0 | Years from drug class introduction to full floor |
| `bacteria_{name}_{drug_class}_resistance_floor` | 0.0 | Target minimum resistance level |

The floor ramps linearly from 0 to the target level over the ramp period, starting from the introduction date of the earliest drug in the relevant class.

### 7.6 Cross-Resistance Groups

For each bacteria, drugs that share resistance mechanisms are grouped together. Acquiring resistance to one drug in a cross-resistance group confers resistance to all drugs in that group. Examples:

| Bacteria | Group | Drugs sharing resistance |
|----------|-------|------------------------|
| E. coli | ESBL | Penicillins + C1-3G + BL/BLI |
| E. coli | FQ | Ciprofloxacin, levofloxacin, moxifloxacin, ofloxacin |
| E. coli | AG | Gentamicin, tobramycin, amikacin |
| S. aureus | β-lactamase | Penicillins |
| S. aureus | MRSA | All cephalosporins |
| S. aureus | MLS | Macrolides + clindamycin |
| A. baumannii | Carbapenemase | All carbapenems + BL/BLI |
| P. aeruginosa | β-lactamase | Piperacillin + cephalosporins + BL/BLI |
| S. pneumoniae | PBP | Penicillin + ampicillin + amoxicillin + BL/BLI |
| E. faecium | Glycopeptide | VanA/VanB (vancomycin + teicoplanin) |

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

The microbiome resistance level (`microbiome_r`) is tracked per bacteria per individual as a value from 0.0 to 1.0.

- **Minority resistance** (<0.5): resistant strains are present but not dominant
- **Majority resistance** (≥0.5): resistant strains dominate the carriage flora

The transition from minority to majority is driven by antibiotic exposure (selective pressure promoting resistant strains):

| Variable | Default | Description |
|----------|---------|-------------|
| `microbiome_majority_promotion_rate_per_day` | 0.02 | Daily probability of minority→majority promotion under drug pressure |

Without antibiotic pressure, resistance in the microbiome decays. This represents fitness cost and competitive displacement by susceptible strains.

### 8.4 Microbiome Disruption

Antibiotic use disrupts the normal microbiome, increasing susceptibility to colonisation by resistant organisms:

| Variable | Default | Description |
|----------|---------|-------------|
| `default_microbiome_disruption_log_odds` | 0.3 | ~1.35× disruption per drug |
| `microbiome_resistance_multiplier_on_acquisition` | 0.50 | Resistance level at acquisition |
| `antibiotic_disruption_decay_half_life_days` | 30.0 | Recovery time after stopping antibiotics |

### 8.5 Majority Resistance Cache

Population-level resistance prevalence (`majority_r`) is computed from a rolling window:

| Variable | Default | Description |
|----------|---------|-------------|
| `majority_r_window_days` | 100 | Window for computing rolling prevalence |
| `majority_r_min_total_samples` | 10 | Minimum samples before prevalence is computed |
| `majority_r_freeze_at_last_positive` | 0.0 | If >0.5, freeze prevalence when data is scarce |

---

## 9. Horizontal Gene Transfer

### 9.1 Overview

Horizontal gene transfer (HGT) allows resistance mechanisms to spread between different bacterial species. The model implements mechanism-driven HGT: only mechanisms that are biologically capable of horizontal transfer (primarily plasmid-borne) can be transferred.

### 9.2 Transfer Probabilities

HGT probability depends on the taxonomic relationship between donor and recipient bacteria. Bacteria are assigned to plasmid pools:

| Pool | Members |
|------|---------|
| `GramPositive` | S. aureus, S. epidermidis, S. pneumoniae, S. pyogenes, S. agalactiae, Enterococcus spp., L. monocytogenes |
| `EntericGramNegative` | E. coli, Klebsiella, Enterobacter, Proteus, Salmonella, Shigella, Citrobacter, Serratia, Morganella, Y. enterocolitica |
| `RespiratoryGramNegative` | P. aeruginosa, A. baumannii, S. maltophilia, B. cepacia, M. catarrhalis, H. influenzae, B. pertussis, Legionella |
| `Anaerobe` | C. difficile, B. fragilis |
| `None` | T. pallidum, H. pylori, C. jejuni, M. tuberculosis, N. gonorrhoeae, N. meningitidis, Chlamydia, Mycoplasma (no HGT) |

| Transfer type | Probability |
|---------------|-------------|
| Same group | 1×10⁻¹⁰ |
| Cross-group (same Gram stain) | 1×10⁻¹¹ |
| **Cross-Gram** | **0.0** (blocked) |
| Excluded organisms (`None` pool) | 0.0 |

### 9.3 HGT Modifiers

| Variable | Default | Description |
|----------|---------|-------------|
| `hgt_base_probability` | 1×10⁻⁵ | Global multiplier |
| `hgt_co_infection_multiplier` | 10.0 | Increased when both species are actively infecting |
| `hgt_hospital_multiplier` | 5.0 | Hospital environment premium |
| `hgt_microbiome_multiplier` | 2.0 | Carriage-to-carriage transfer |
| `hgt_gut_compartment_multiplier` | 2.0 | Additional multiplier for gut-compartment bacteria (dense microbial environment) |

### 9.4 Mechanism Transferability

Only certain mechanisms are eligible for HGT. Transferable mechanisms include plasmid-borne determinants (ESBLs, carbapenemases, MCR-1, etc.). Chromosomal mutations (GyrA, RpoB, porin loss) are not transferable.

---

## 10. Mortality

### 10.1 Background Mortality

Non-infection-related death is modelled to maintain realistic population demographics:

$$\text{log\_odds} = \text{base} + \text{age\_linear} + \text{age\_quadratic} + \text{era} + \text{region} + \text{sex} + \text{immune} + \text{hospital}$$

| Variable | Default | Description |
|----------|---------|-------------|
| `background_mortality_baseline_log_odds` | −14.0 | Very low daily baseline |
| `log_odds_mortality_per_year_of_age` | 0.04 | Linear age effect |
| `log_odds_mortality_per_year_of_age_squared` | 0.05 | Quadratic age effect (accelerating) |
| `mortality_baseline_1930_multiplier` | 3.0 | Higher mortality in 1930 |
| `mortality_baseline_2035_multiplier` | 1.0 | Current baseline |
| `mortality_improvement_half_life_years` | 35.0 | Rate of improvement |
| `log_odds_mortality_immunosuppressed` | 0.916 | ln(2.5), ~2.5× |
| `log_odds_mortality_hospitalized` | 0.262 | ln(1.3), ~1.3× |

#### Regional mortality log-odds

| Region | Value |
|--------|-------|
| North America | 0.0 |
| Europe | −0.105 |
| Oceania | 0.0 |
| Asia | 0.18 |
| South America | 0.26 |
| Africa | 0.69 |

#### Sex effect

| Sex | Value |
|-----|-------|
| Male | 0.095 |
| Female | −0.105 |

### 10.2 Sepsis Mortality

Sepsis death follows a logistic model evaluated daily for individuals in sepsis:

| Variable | Default | Description |
|----------|---------|-------------|
| `sepsis_death_base_log_odds` | −5.0 | Daily baseline during sepsis |
| `sepsis_death_log_odds_age_infant` | 1.1 | |
| `sepsis_death_log_odds_age_child` | −0.7 | Children more resilient |
| `sepsis_death_log_odds_age_adult` | 0.0 | Reference |
| `sepsis_death_log_odds_age_elderly` | 0.9 | |
| `sepsis_death_log_odds_immunosuppressed` | 1.5 | ~4.5× |
| `sepsis_death_log_odds_bacteria_level` | 0.35 | Per unit of bacterial load |
| `sepsis_death_log_odds_duration` | 0.04 | Per day of sepsis |
| `sepsis_death_log_odds_early_phase` | 0.8 | Higher in first 3 days |
| `sepsis_death_early_phase_days` | 3.0 | |
| `sepsis_death_log_odds_not_under_care` | 1.4 | ~4× if not under care |

#### Regional sepsis mortality multipliers

| Region | Value |
|--------|-------|
| North America | 0.5 |
| Europe | 0.4 |
| Oceania | 0.5 |
| Asia | 0.9 |
| South America | 1.1 |
| Africa | 1.5 |

### 10.3 Sepsis Recovery

| Variable | Default | Description |
|----------|---------|-------------|
| `sepsis_base_log_odds_of_recovery_per_day` | 0.0 | Base daily recovery probability |
| `sepsis_log_odds_bacteria_level` | −0.3 | Higher load → harder to recover |
| `sepsis_log_odds_in_hospital` | 0.8 | Hospital care assists recovery |
| `sepsis_log_odds_age_infant` | −0.5 | |
| `sepsis_log_odds_age_child` | 0.4 | Children recover faster |
| `sepsis_log_odds_age_adult` | 0.0 | Reference |
| `sepsis_log_odds_age_elderly` | −0.7 | |
| `sepsis_log_odds_immunosuppressed` | −1.0 | |
| `sepsis_minimum_duration_days` | 1.0 | |

#### Regional sepsis recovery log-odds

| Region | Value |
|--------|-------|
| North America | 0.4 |
| Europe | 0.5 |
| Oceania | 0.3 |
| Asia | 0.0 |
| South America | −0.3 |
| Africa | −0.7 |

---

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

### A.2 Antibiotics (58 drugs in 18 classes)

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
| cefepime | Cephalosporins 4–5G |
| ceftaroline | Cephalosporins 4–5G |
| ceftolozane_tazobactam | Cephalosporins 3G |
| cefiderocol | Cephalosporins 4–5G |
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

### A.3 Drug Classes (18)

| Code | Full name | Drugs |
|------|-----------|-------|
| `pen` | Penicillins | penicillin G, ampicillin, amoxicillin, piperacillin, ticarcillin |
| `bli` | BLI Combinations | amoxicillin-clavulanate, piperacillin-tazobactam, ampicillin-sulbactam, ticarcillin-clavulanate |
| `c1_2g` | Cephalosporins 1–2G | cephalexin, cefazolin, cefuroxime |
| `c3g` | Cephalosporins 3G | ceftriaxone, ceftazidime, ceftolozane-tazobactam |
| `c4_5g` | Cephalosporins 4–5G | cefepime, ceftaroline, cefiderocol |
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

| Variable | Default | Description |
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

| Variable | Default |
|----------|---------|
| `antibiotic_initiation_base_log_odds` | −6.5 |
| `antibiotic_initiation_log_odds_symptomatic_infection` | 6.5 |
| `antibiotic_initiation_log_odds_sepsis` | 6.0 |
| `antibiotic_initiation_log_odds_test_identified` | 0.92 |
| `antibiotic_initiation_log_odds_already_on_drug` | 0.18 |
| `antibiotic_initiation_log_odds_immunodeficiency` | 2.08 |
| `antibiotic_initiation_log_odds_no_indication` | −1.05 |

#### Antibiotic treatment cessation

| Variable | Default |
|----------|---------|
| `treatment_stop_improvement_threshold` | 2.0 |
| `treatment_stop_rate_per_day` | 0.03 |
| `treatment_duration_base_days` | 7.0 |

#### Drug efficacy

| Variable | Default |
|----------|---------|
| `drug_effect_on_bacteria_per_day` | 0.5 |
| `drug_minimum_effective_level` | 0.1 |

#### Drug selection

| Variable | Default |
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

| Variable | Default |
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

| Variable | Default |
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

| Variable | Default |
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

| Variable | Default |
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

| Variable | Default |
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

| Variable | Default |
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

| Variable | Default |
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

| Variable | Default |
|----------|---------|
| `mechanism_assignment_probability` | 0.8 |

#### Microbiome and carriage

| Variable | Default |
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

| Variable | Default |
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

| Variable | Default |
|----------|---------|
| `hgt_base_probability` | 1×10⁻⁵ |
| `hgt_co_infection_multiplier` | 10.0 |
| `hgt_hospital_multiplier` | 5.0 |
| `hgt_microbiome_multiplier` | 2.0 |
| `hgt_gut_compartment_multiplier` | 2.0 |

#### Travel

| Variable | Default |
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

### B.8 Drug Penetration Table (10 syndromes × 18 drug classes)

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
