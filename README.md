
We are developing a stochastic individual-based model for anti-bacterial resistance. The engine tracks every infection, resistance state, and treatment decision for each simulated person on a daily timestep starting in 1942 (the introduction of penicillin), although alternative start dates with pre-existing resistance are also supported.

Age is represented in days; negative ages indicate unborn individuals who remain inert until birth. Each person belongs to a home region but may be temporarily located elsewhere, letting the model differentiate home exposure, travel, and hospital-acquired events. The simulation currently includes 39 bacteria and 52 antibiotics, but both lists can grow without changing core data structures.

### Resistance Calculation

Resistance tracking relies on five per-bacteria/drug metrics defined in [src/simulation/population.rs](src/simulation/population.rs):

- `any_r`: fractional resistance in the infection (0 = completely susceptible). It is the aggregate resistance from all active mechanisms.
- `majority_r`: identical to `any_r` but only recorded when resistant strains are the majority of the infection.
- `activity_r`: current killing power of the drug given its intrinsic potency, resistance level, current drug concentration, **syndrome-specific tissue penetration, and time-based accumulation**.
- `test_r`: resistance level last confirmed by diagnostic testing.
- `microbiome_r`: resistance carried in the microbiome rather than the active infection.

**Resistance Mechanism Stacking (New in v1.1)**
The model now employs **Multiplicative Stacking** for resistance mechanisms. Instead of simply using the strongest mechanism present ("High Water Mark"), multiple mechanisms now combine to increase resistance toward 1.0 (complete immunity).

$$
TotalSusceptibility = (1.0 - Mechanism_A) \times (1.0 - Mechanism_B) \\
TotalResistance = 1.0 - TotalSusceptibility
$$

For example, acquiring two independent mechanisms that each provide 50% resistance results in 75% total resistance ($1.0 - (0.5 \times 0.5) = 0.75$), significantly increasing the survival of MDR strains.

### Drug Activity Calculation

The `activity_r` metric determines the effective antibacterial action **at the infection site**, not just in the bloodstream. It incorporates pharmacokinetic realism through syndrome-specific penetration factors ([src/rules/mod.rs](src/rules/mod.rs)):

```
activity_r = base_potency × effective_drug_level × (1.0 - normalized_resistance)

where:
  effective_drug_level = serum_level × penetration_factor
  
  penetration_factor = syndrome_drug_penetration[syndrome_id][drug_idx]  (0.0–1.0)
  normalized_resistance = any_r / max_resistance_level
```

This captures clinical reality where drug concentration varies wildly by tissue:
- **Site-Specific Efficacy**: A drug with high serum levels may still fail if it cannot penetrate the infection site (e.g., Aminoglycosides in CNS).
- **Aminoglycosides** achieve only 5% CNS penetration due to the blood-brain barrier.
- **Fluoroquinolones** achieve 90% prostatic penetration, making them preferred for prostatitis.
- **Bone/joint infections** have poor vascularity requiring prolonged therapy with good penetrating agents.

### Bacteria Growth and Host Factors

Bacteria level changes incorporate host-specific multipliers ([src/rules/mod.rs#L4048-L4095](src/rules/mod.rs#L4048-L4095)):

| Factor | Categories | Multipliers |
| --- | --- | --- |
| **Age** | Infant (≤1y), Child (≤18y), Adult (≤65y), Elderly (>65y) | 1.3×, 1.0×, 1.0×, 1.2× |
| **Immunodeficiency** | Present/Absent | 1.5× / 1.0× |
| **Syndrome** | Bloodstream, CNS, Bone/joint, UTI, etc. | 1.4×, 1.3×, 0.85×, 1.0×, etc. |

### Bacteria taxonomy & carriage sites
Every bacteria listed in [`BACTERIA_LIST`](src/simulation/population.rs) now carries lightweight metadata that the rules engine can query without scanning strings:

- `BacteriaGroup` enumerations in [src/simulation/population.rs](src/simulation/population.rs) tag each organism as Gram-positive, Enterobacterales, non-fermenter, etc.
- `CarriageCompartment` enumerations in [src/simulation/population.rs](src/simulation/population.rs) capture whether a bacteria’s default reservoir is gut, respiratory, skin/soft tissue, genitourinary, or systemic.
- Helper functions convert those enumerations into bit masks for fast lookup during horizontal gene transfer (HGT) checks.

**Transmission Logic & "Frankenstein" Bug Fix**
The simulation uses a `MechanismPrevalenceCache` to track the prevalence of specific resistance genotypes (mechanism combinations) in the population. When a new infection is transmitted, the recipient samples a *single* coherent genotype from the donor pool for that region/bacteria/hospital-status triplet. This prevents the "Frankenstein Strain" artifact where infections previously sampled a random independent mechanism for *every single drug*, creating impossibly resistant superbugs that did not exist in nature.

**Novel Resistance Mechanisms**
The system supports defining "Other" mechanism slots (e.g., `OtherMechanism1`, `OtherMechanism2`) as bacteria-specific resistance traits in `src/config.rs`. This allows modeling emerging threats like novel efflux pumps or plasmids that are specific to one species without polluting the global mechanism table.

The level of any antibiotic is kept on a standardized 0–10 scale per day of standard dosing (double doses reach 20). Drug levels decay after cessation, but residual levels continue to influence `activity_r`. Testing variables capture whether the causative bacteria has been identified and whether resistance testing was performed. Exposure multipliers (sexual, airborne adult/child, oral, mosquito) combine with age and region to set acquisition risks; when acquisition is person-to-person the new infection inherits `any_r` from a sampled source in the same region.

Daily immune clearance hazards are stored per bacteria, with configurable arming delays and multipliers for age, immunodeficiency, and infection level. Background, sepsis-related, and drug-toxicity mortality components are separated, allowing interventions to target individual pathways.

## Simulation State Reference

All variables below are taken directly from the `Individual` struct in [src/simulation/population.rs#L443-L640](src/simulation/population.rs#L443-L640). Vector fields have length equal to the number of tracked bacteria or drugs.

### Demographics & Mobility
| Variable | Type | Description |
| --- | --- | --- |
| id | usize | Stable identifier assigned when the person is created. |
| age | i32 | Age in days; negative values mean the individual has not yet been born, so no events occur until age reaches 0. |
| sex_at_birth | String | Recorded sex at birth used for mortality modifiers. |
| perceived_penicillin_allergy | bool | Whether clinicians believe the individual has a penicillin allergy, influencing drug choice logic. |
| region_living | Region | Home region (e.g., north_america, asia) used for demographics and acquisition baselines. |
| region_cur_in | Region | Region currently inhabited; differs from home while traveling. |
| days_visiting | u32 | Days spent away from the home region on the current trip. |
| hospital_status | HospitalStatus | Indicates whether the person is currently hospitalized, affecting acquisition risks and mortality. |
| days_hospitalized | u32 | Duration of the current hospital stay for LOS-dependent processes. |

### Infection Timeline & Severity Metrics
| Variable | Type | Description |
| --- | --- | --- |
| date_last_infected | Vec<i32> | Start date of the current/most recent infection per bacteria (0 if never infected). |
| date_last_infected_keep | Vec<i32> | Historical start date that is never reset, enabling cumulative reporting. |
| infectious_syndrome | Vec<i32> | Syndrome/site identifier associated with each infection episode. |
| level | Vec<f64> | Current bacteria burden (0–max) that drives symptoms, testing, and drug efficacy. |
| predicted_infection_risk | Vec<f64> | Logistic model prediction for acquisition probability on the current day. |
| cur_infection_from_environment | Vec<bool> | Flags infections currently attributed to environmental sources. |
| infection_hospital_acquired | Vec<bool> | Marks infections caught within hospital settings. |
| infection_has_caused_symptoms | Vec<bool> | Indicates whether an active infection has already produced clinical symptoms; stays true until clearance. |

### Clearance, Sepsis, and Mortality Inputs
| Variable | Type | Description |
| --- | --- | --- |
| clearance_hazard | Vec<f64> | Daily immune-mediated clearance hazard stored for reporting. |
| clearance_ready_day | Vec<i32> | Simulation day when clearance hazards activate (−1 = not yet armed). |
| sepsis | Vec<bool> | Whether each infection has progressed to sepsis or an equivalent life-threatening condition. |
| sepsis_onset_day | Vec<i32> | Day sepsis started per bacteria (−1 if never observed). |
| infection_prevented_by_drug | Vec<bool> | True when an infection attempt was blocked by ongoing therapy during the current timestep. |
| current_infection_related_death_risk | f64 | Aggregated daily probability that infection-driven causes will kill the individual. |
| background_all_cause_mortality_rate | f64 | Baseline mortality hazard from non-infection causes, derived from age, region, and calendar year. |

### Microbiome & Vaccination Traces
| Variable | Type | Description |
| --- | --- | --- |
| presence_microbiome | Vec<bool> | Records whether the person carries each bacteria in their microbiome. |
| date_microbiome_acquired | Vec<i32> | Day microbiome carriage was gained (0 if never). |
| microbiome_acquired_today | Vec<bool> | Flags new microbiome acquisition events during the current day; cleared after aggregation. |
| microbiome_acquired_on_drug_today | Vec<bool> | Notes whether acquisition happened while any antibiotic was active. |
| microbiome_cleared_today | Vec<bool> | Flags microbiome clearance events during the current day. |
| cleared_any_r_microbiome_categories | Vec<[u32;4]> | Counts infection clearances broken down by microbiome resistance category for reporting. |
| vaccination_status | Vec<bool> | Per-bacteria vaccination flag, updated dynamically based on eligibility windows. |

### Diagnostics & Testing
| Variable | Type | Description |
| --- | --- | --- |
| test_identified_infection | Vec<bool> | True once laboratory testing has identified the causative bacteria for the current infection. |
| test_for_resistance | Vec<bool> | Tracks whether resistance testing has been performed for each bacteria. |
| resistance_test_initiated_day | Vec<i32> | Day resistance testing started (−1 if never initiated). |

### Drug Exposure & Toxicity
| Variable | Type | Description |
| --- | --- | --- |
| cur_use_drug | Vec<bool> | Indicates whether each drug is currently being administered. |
| cur_level_drug | Vec<f64> | Standardized drug level (typically 0–10) after accounting for dosing and decay. |
| date_drug_initiated | Vec<i32> | Simulation day that treatment with each drug last began. |
| date_drug_initiated_keep | Vec<i32> | Historical record of drug starts that is never reset. |
| ever_taken_drug | Vec<bool> | True if the person has ever received each drug. |
| drug_toxicity_reservoir | Vec<f64> | Sub-acute toxicity burden that decays per drug-specific half-life. |
| current_toxicity_hazard | f64 | Daily probability of drug toxicity death given the summed reservoir. |
| mortality_risk_current_toxicity | f64 | Instantaneous mortality risk attributable to toxicity in the current timestep. |

### Resistance, Outcomes, and Surveillance
| Variable | Type | Description |
| --- | --- | --- |
| resistances | Vec<Vec<Resistance>> | Full resistance metrics (`any_r`, `majority_r`, `activity_r`, `test_r`, `microbiome_r`) for each bacteria/drug combination. |
| resistance_mechanisms | Vec<Vec<bool>> | Whether each enumerated mechanism is active for the given bacteria. |
| how_resistance_acquired | Vec<Vec<Option<ResistanceAcquisitionType>>> | Records whether resistance emerged endogenously, via transmission, or through microbiome transfer. |
| infection_resolution_this_timestep | Vec<Vec<u32>> | Counts of resolution outcomes (e.g., immune clearance, drug clearance) accrued during the current day. |
| day_7_since_last_infection_drug_used | Vec<Option<bool>> | On day 7 after infection onset, records whether any drug was initiated (Some(true/false)); None otherwise. |
| date_last_drug_failure | Vec<i32> | Most recent day a treatment failure was logged for each bacteria. |

### Treatment Progression & Cessation Tracking
| Variable | Type | Description |
| --- | --- | --- |
| bacteria_level_at_drug_start | Vec<Option<f64>> | Infection level when the current treatment began (None if not on therapy). |
| days_on_current_treatment | Vec<i32> | Days elapsed since the active drug course started (−1 if no active treatment). |
| treatment_failure_assessed | Vec<bool> | Whether treatment failure has already been evaluated for the ongoing course. |
| drug_activity_response_multiplier | Vec<f64> | Per-bacteria scaling of how drug activity translates to bacterial kill rate, sampled at initiation. |
| drug_stopped_with_infection_day | Vec<Option<i32>> | Day a drug was stopped while infection was still present (None if never happened). |
| bacteria_level_at_drug_cessation | Vec<Option<f64>> | Infection level recorded when a drug was stopped due to non-adherence or adverse events. |
| stopped_drug_index | Vec<Option<usize>> | Identifies which drug was stopped while infection remained. |
| restart_window_assessed | Vec<bool> | True once the restart eligibility window has been evaluated after a cessation. |

### Drug Selection, Scoring, and Load
| Variable | Type | Description |
| --- | --- | --- |
| bacteria_on_selection_day | i32 | Bacteria index that triggered drug selection during the current day (−1 if no selection). |
| drug_score_on_selection_day | Vec<f64> | Cached empiric scores for all drugs with respect to the triggering bacteria (−1.0 if inactive). |
| current_number_of_drugs | i32 | Count of simultaneously administered antibiotics. |

### Mortality & Terminal Status
| Variable | Type | Description |
| --- | --- | --- |
| date_of_death | Option<usize> | Simulation day of death if the individual has died. |
| cause_of_death | Option<String> | Textual cause recorded at death (infection, toxicity, baseline, etc.). |
| immunodeficiency_type | Option<ImmunodeficiencyType> | Current severe immunodeficiency state, driving prophylaxis and clearance modifiers. |

### Resistance Metric Fields
| Field | Range | Meaning |
| --- | --- | --- |
| any_r | 0–1 | Fraction of the infection exhibiting resistance to the drug  *** ; copied at transmission or increased through emergence. |
| majority_r | 0–1 | Equals `any_r` when resistant strains are the majority; otherwise remains 0 to distinguish minority carriage. |
| activity_r | 0–1 | Effective potency of the drug against the infection after accounting for current concentration and resistance. |
| test_r | 0–1 | Resistant fraction last measured via diagnostic testing, used to drive targeted therapy choices. |
| microbiome_r | 0–1 | Resistance level carried in microbiome bacteria rather than the active infection. |

## Parameter Reference

All parameter sets originate from `ParameterStore` in [src/config.rs#L108-L173](src/config.rs#L108-L173). Each block is instantiated by parsing a keyed map, so every field shown in code can be tuned in external configuration files. The table below summarizes how each block influences model behavior and points to its definition.

| Parameter block | Definition | Primary influence on the simulation | Example knobs |
| --- | --- | --- | --- |
| ParameterStore | [src/config.rs#L108-L173](src/config.rs#L108-L173) | Top-level container that wires every parameter subset into the simulation. Populated once at startup and shared globally. | Controls which parameter files are loaded and ensures bacteria/drug counts stay consistent with `BACTERIA_LIST` and `DRUG_SHORT_NAMES`. |
| GlobalScalars | [src/config.rs#L174-L739](src/config.rs#L174-L739) | Whole-population multipliers controlling drug initiation/stopping, hospital flows, mortality, resistance emergence, testing bonuses, and sepsis logic. | `drug_base_initiation_rate_per_day` sets the baseline hazard for starting antibiotics; `hospital_*` values shape admission/discharge; `resistance_emergence_*` tune how quickly `any_r` grows; `majority_r_*` configure the reporting window for majority resistance. |
| ImmunodeficiencyParameters | [src/config.rs#L872-L931](src/config.rs#L872-L931) | Rates of temporary/chronic immunosuppression onset/recovery along with age-stratified prevalence, feeding both prophylaxis and clearance modifiers. | `temporary_immunosuppression_onset_rate_per_day` controls how often people enter short-term suppression; `chronic_probability_age_*` sets baseline prevalence per age band. |
| RegionParameters | [src/config.rs#L740-L834](src/config.rs#L740-L834) | Region-specific multipliers for travel, treatment cessation, mortality, sepsis lethality, and testing intensity. | `asia_travel_multiplier` influences how often travellers leave Asia; `europe_testing_multiplier` boosts identification rates for infections in Europe. |
| SexParameters | [src/config.rs#L836-L860](src/config.rs#L836-L860) | Sex-at-birth mortality adjustments that stack with age/region baselines. | `log_odds_mortality_sex_male` and `log_odds_mortality_sex_female` shift mortality hazards for each sex. |
| VaccinationParameters | [src/config.rs#L862-L940](src/config.rs#L862-L940) | Age- and vaccine-specific daily immunization probabilities plus availability start years for pneumococcal, meningococcal, and Hib vaccines. | `vaccine_pneumococcal_daily_prob_age_child` sets routine uptake; `vaccine_hib_availability_year` defines rollout timing. |
| SyndromeParameters | [src/config.rs#L1071-L1460](src/config.rs#L1071-L1460) | Syndrome-level sepsis odds, initiation multipliers, non-sepsis mortality, empiric drug scoring matrices, bacteria growth multipliers, and drug penetration factors by compartment. | `syndrome_3_initiation_multiplier` accelerates treatment for syndrome 3; `syndrome_6_drug_gentamicin_penetration` (default 0.05) reflects poor BBB penetration. |
| DrugParameters | [src/config.rs#L1100-L1194](src/config.rs#L1100-L1194) | Drug-specific pharmacokinetics and toxicity, including starting level, double-dose multiplier, half-life, microbiome disruption, and toxicity reservoir decay. | `drug_ciprofloxacin_half_life_days` drives decay of drug levels; `drug_vancomycin_microbiome_disruption_log_odds` captures dysbiosis risk. |
| BacteriaParameters | [src/config.rs#L1196-L1480](src/config.rs#L1196-L1480) | Per-bacteria acquisition odds, vaccination effects, symptom dynamics, sepsis risks, microbiome clearances, and infection level kinetics. | `klebsiella_pneumoniae_acquisition_log_odds_baseline` controls baseline exposure; `staphylococcus_aureus_max_level` caps infection burden. |
| ClearanceParameters | [src/config.rs#L1506-L1597](src/config.rs#L1506-L1597) | Immune clearance timing and hazard, including age multipliers, immunodeficiency scaling, and bacteria-specific overrides. | `default_clearance_delay_days` postpones immune action after acquisition; `pseudomonas_aeruginosa_clearance_hazard_multiplier` tailors persistence. |
| DrugBacteriaMatrix | [src/config.rs#L1623-L1698](src/config.rs#L1623-L1698) | Potency lookup, initiation multipliers, resistance emergence rates, and MIC-derived thresholds for every bacteria–drug pair. | `drug_meropenem_for_bacteria_acinetobacter_baumannii_potency_when_no_r` encodes kill strength when susceptible; `_resistance_emergence_rate_per_day_baseline` tunes mutation pressure. |
| RegionBacteriaAcquisition | [src/config.rs#L1700-L1755](src/config.rs#L1700-L1755) | Region-specific acquisition log-odds for each bacteria, allowing local epidemiology or healthcare settings to differ. | `asia_klebsiella_pneumoniae_acquisition_log_odds` raises or lowers Asian exposure independent of other regions. |
| AgeTables | [src/config.rs#L1756-L1829](src/config.rs#L1756-L1829) | Pre-computed age log-odds tables used when detailed multipliers are required for reporting or defaults. | `default_log_odds_adult` supplies fallback risk when bacteria-specific data is missing. |
| AgeCategoryParameters | [src/config.rs#L1831-L1933](src/config.rs#L1831-L1933) | Combined bacteria, region, and age log-odds that drive acquisition, enabling fine-grained age-structured risks. | `streptococcus_pneumoniae_log_odds_child` modifies odds for children; `asia_log_odds_elderly` adjusts regional age effects. |
| HgtMatrix | [src/config.rs#L1972-L2074](src/config.rs#L1972-L2074) | Horizontal gene transfer probabilities from donor to recipient bacteria, with defaults derived from group-based compatibility and plasmid pools. | Override `hgt_prob_*` keys (seeded by `default_hgt_probability`) to deviate from the Gram-positive/Enteric/Respiratory baseline. |
| ResistanceMechanismParameters | [src/config.rs#L2077-L2129](src/config.rs#L2077-L2129) | Global emergence, enhancement, and reversion rates for each enumerated mechanism (e.g., ESBL, carbapenemase). | `resistance_mechanism_esbl_emergence_rate` tunes how often the mechanism activates; `_reversion_rate` governs loss when pressure lifts. |
| BacteriaMechanismEmergenceMultipliers | [src/config.rs#L2132-L2166](src/config.rs#L2132-L2166) | Bacteria-specific multipliers applied on top of global mechanism rates, capturing organism-level propensities. | `bacteria_klebsiella_pneumoniae_mechanism_esbl_emergence_multiplier` increases ESBL emergence relative to baseline. |

The parameter system is extensible: new bacteria, drugs, or mechanisms become available simply by extending the relevant lists in [src/simulation/population.rs](src/simulation/population.rs) and adding keyed entries to the configuration map.

## Drug Selection Algorithms

Antibiotic initiation happens in two explicit stages implemented in [src/rules/mod.rs](src/rules/mod.rs#L1410-L2540):

1. **Start decision** – once per person per day we compute a probability of starting any antibiotic based on the globals discussed below.
2. **Drug choice** – if stage 1 succeeds, we score every available drug and draw one using a temperature-controlled softmax so higher scores dominate but clinical variability remains realistic.

### Stage 1: Deciding to treat

The model filters out drugs that have not yet been introduced or are unavailable in either the traveller’s current region or home region ([src/rules/mod.rs#L1410-L1436](src/rules/mod.rs#L1410-L1436)). With that candidate list, the start probability `start_any_antibiotic_prob` is assembled as follows ([src/rules/mod.rs#L1437-L1505](src/rules/mod.rs#L1437-L1505)):

- Begin with `drug_base_initiation_rate_per_day` and scale it up if only a handful of drugs are on formulary or if the individual has symptomatic infection, a confirmed pathogen, or is already on therapy.
- Multiply by `syndrome_administration_multiplier` for the active syndrome, and by the prophylaxis multiplier when the person is immunodeficient. A separate misdiagnosis multiplier allows low-level empiric starts even without symptoms.
- Clamp to $[0,1]$ and draw a Bernoulli. We hard-cap concurrent regimens at three to keep polypharmacy realistic.

### Stage 2A: Empiric (syndrome-driven) selection

When no infection has been lab-confirmed, the selection logic treats the case as empiric therapy ([src/rules/mod.rs#L1516-L2288](src/rules/mod.rs#L1516-L2288)). Important steps:

- **Syndrome multipliers** – all active syndromes contribute multiplicative weights from `SyndromeParameters.empiric_drug_score()` so, for example, meningitis syndromes boost ceftriaxone and penicillin ([src/rules/mod.rs#L1516-L1605](src/rules/mod.rs#L1516-L1605)).
- **Reserve and allergy gates** – penicillin allergies block the entire penicillin class; reserve agents (carbapenems, colistin, linezolid, etc.) cannot be chosen unless there was a documented drug failure in the recent `drug_failure_memory_days` window *and* the regional majority-resistance cache shows high resistance pressure ([src/rules/mod.rs#L1526-L1558](src/rules/mod.rs#L1526-L1558), [src/rules/mod.rs#L2143-L2235](src/rules/mod.rs#L2143-L2235)).
- **Regional resistance penalty** – before species are known we down-weight drugs that the surveillance cache marks as ineffective in the person’s region/hospital context using the configurable moderate/high/very-high thresholds ([src/rules/mod.rs#L2056-L2139](src/rules/mod.rs#L2056-L2139)).
- **Spectrum preferences** – empiric broad-spectrum bonuses and narrow-spectrum penalties are applied using `spectrum_breadth`, while drugs with no syndrome signal get the empiric-ineffective penalty ([src/rules/mod.rs#L2235-L2275](src/rules/mod.rs#L2235-L2275)).

Only drugs with positive final scores survive to the softmax draw, and the selected agent is started at its initial level (with future double-dose options reserved for targeted therapy). Score traces are written back to `drug_score_on_selection_day` for downstream analytics.

### Stage 2B: Targeted therapy after pathogen identification

If any infection has been both symptomatic and lab-confirmed, the algorithm switches to targeted mode ([src/rules/mod.rs#L1516-L2288](src/rules/mod.rs#L1516-L2288)). Key differences from empiric care include:

- **Resistance-test filter** – drugs with a positive `test_r` result for any active infection are removed once the configured lab delay has elapsed, preventing re-use of known failures ([src/rules/mod.rs#L1526-L1554](src/rules/mod.rs#L1526-L1554)).
- **Species-specific guidelines** – every identified bacteria applies aggressive multipliers or hard blocks that reflect clinical guidelines (e.g., anti-pseudomonal agents, MRSA-era vancomycin/linezolid boosts, TMP-SMX for *Stenotrophomonas*) as coded in [src/rules/mod.rs#L1640-L2040](src/rules/mod.rs#L1640-L2040).
- **First/second-line enforcement** – after guideline boosts we penalize off-guideline drugs unless they are in the curated first/second-line set for that pathogen ([src/rules/mod.rs#L1960-L2036](src/rules/mod.rs#L1960-L2036)).
- **Potency reinforcement** – a drug’s maximum potency against the confirmed infections yields large multiplicative rewards (6–15×) and must exceed `minimal_potency_threshold` to stay eligible ([src/rules/mod.rs#L1897-L2036](src/rules/mod.rs#L1897-L2036)).
- **Reserve stewardship** – reserve agents remain heavily penalized unless there is a recorded failure within `drug_failure_memory_days`, even during targeted therapy ([src/rules/mod.rs#L2143-L2192](src/rules/mod.rs#L2143-L2192)).
- **Spectrum balancing** – effective narrow-spectrum drugs get `targeted_therapy_narrow_spectrum_bonus`, while broad agents get penalized to keep stewardship-intentional prescribing ([src/rules/mod.rs#L2192-L2234](src/rules/mod.rs#L2192-L2234)). Drugs lacking potency trigger the `targeted_therapy_ineffective_drug_penalty`.
- **Final guards** – we also multiply in drug availability/launch timing and apply the global reserve penalty before running the softmax and initiating the chosen regimen ([src/rules/mod.rs#L2245-L2340](src/rules/mod.rs#L2245-L2340)).

Both empiric and targeted paths share the same stochastic selection temperature (`drug_selection_temperature`). After a drug is initiated the model records start dates, optionally applies double-dosing for targeted therapy (`double_dose_probability_if_identified_infection`), resets treatment-failure tracking, and immediately updates toxicity reservoirs for downstream mortality logic ([src/rules/mod.rs#L2334-L2540](src/rules/mod.rs#L2334-L2540)).

### GlobalScalars parameters

Derived from [src/config.rs#L174-L739](src/config.rs#L174-L739). Values affect every individual each timestep.

#### Drug initiation, cessation, and travel
| Parameter | Description |
| --- | --- |
| `drug_base_initiation_rate_per_day` | Baseline daily hazard of starting any antibiotic before multiplying by infection-specific factors. |
| `drug_infection_present_multiplier` | Multiplier applied when a pathogen is currently present, sharply increasing initiation odds. |
| `already_on_drug_initiation_multiplier` | Additional multiplier if the person already has another drug active. |
| `drug_test_identified_multiplier` | Boost applied once laboratory testing identifies the pathogen. |
| `double_dose_probability_if_identified_infection` | Probability that initiation uses a double dose when the infection is confirmed. |
| `random_drug_cessation_probability` | Daily hazard for halting therapy regardless of infection status. |
| `random_drug_cessation_probability_if_no_active_infection` | Higher cessation hazard once the infection has resolved. |
| `misdiagnosis_antibiotic_initiation_multiplier` | Multiplier applied when the infection is misclassified, enabling empiric use. |
| `immunodeficiency_prophylactic_drug_multiplier` | Raises initiation propensity for immunosuppressed individuals. |
| `travel_probability_per_day` | Baseline probability that the person begins travel on any day. |
| `hospital_baseline_rate_per_day` | Baseline per-day hospitalization probability. |
| `hospital_age_multiplier_per_day` | Age-based multiplier applied on top of the baseline admission hazard. |
| `hospital_recovery_rate_per_day` | Per-day discharge probability once admitted. |
| `hospital_max_days` | Cap on continuous days in hospital regardless of hazards. |
| `hospital_sepsis_admission_multiplier` | Admission hazard multiplier when the person is septic. |
| `hospital_prevent_discharge_with_sepsis` | Fraction of discharges that are suppressed while sepsis is ongoing. |

#### Drug activity, restart, and scoring
| Parameter | Description |
| --- | --- |
| `drug_activity_to_bacteria_level_multiplier` | Converts summed drug activity into bacteria-level reductions. |
| `drug_activity_slow_clearance_probability` | Chance of sampling the "slow clearance" path when drug activity is low. |
| `drug_activity_slow_clearance_multiplier` | Activity multiplier applied when slow-clearance is triggered. |
| `antibiotic_infection_prevention_efficacy` | Fraction of exposure events prevented by ongoing therapy. |
| `minimal_potency_threshold_for_drug_selection` | Drops drugs whose potency falls below this level during selection. |
| `drug_selection_temperature` | Softmax temperature controlling exploration vs exploitation when ranking drugs. |
| `reserve_drug_score_penalty` | Penalty applied to reserve-class drugs to keep them in reserve. |
| `restart_window_enabled` | Enables the logic that restarts a recently stopped drug if infection persists. |
| `restart_window_days` | Length of the restart eligibility window. |
| `restart_bacteria_level_threshold` | Minimum bacteria level required to consider restarting. |
| `restart_window_probability` | Daily chance of restarting within the eligible window. |

#### Resistance emergence and mechanisms
| Parameter | Description |
| --- | --- |
| `max_resistance_level` | Upper bound for `any_r`/`majority_r` values. |
| `resistance_emergence_bacteria_level_multiplier` | Scales emergence probability by infection burden. |
| `resistance_emergence_pop_size_multiplier` | Debug multiplier to keep prevalence stable when population size changes. |
| `any_r_emergence_level_on_first_emergence` | Initial resistance level assigned the first time a strain turns resistant. |
| `multi_drug_penalty_threshold_num_drugs` | Number of concurrent drugs that triggers additional penalties. |
| `resistance_development_inhibition_single_drug` | Fractional suppression of resistance growth when only one drug is in use. |
| `resistance_development_inhibition_partial_cross` | Suppression term modeling partial cross-protection between drugs. |
| `mechanism_assignment_probability_on_any_r_gain` | Chance of tagging a resistance mechanism each time `any_r` increases. |
| `microbiome_resistance_transfer_probability_per_day` | Baseline probability that microbiome resistance transfers into infection. |
| `microbiome_resistance_emergence_rate_per_day_baseline` | Base rate at which microbiome resistance escalates even without infection. |
| `mdr_tb_pre_antibiotic_era_multiplier` | Historical adjustment to MDR TB levels before antibiotics. |
| `mdr_tb_early_antibiotic_era_multiplier` | Adjustment during the early antibiotic era. |
| `mdr_tb_modern_era_multiplier` | Adjustment during the modern era. |

#### Treatment failure and restart recording
| Parameter | Description |
| --- | --- |
| `treatment_failure_enabled` | Master switch for treatment-failure tracking. |
| `treatment_failure_assessment_day` | Number of days after initiation when failure is evaluated. |
| `treatment_failure_threshold` | Bacteria level above which failure is declared. |
| `drug_failure_memory_days` | Lookback window when counting past failures for selection penalties. |

#### Toxicity and regional resistance penalties
| Parameter | Description |
| --- | --- |
| `default_toxicity_reservoir_half_life_days` | Default decay half-life for the toxicity reservoir per drug. |
| `toxicity_age_multiplier_infant` | Infant-specific multiplier on toxicity hazards. |
| `toxicity_age_multiplier_child` | Child multiplier. |
| `toxicity_age_multiplier_adult` | Adult multiplier. |
| `toxicity_age_multiplier_elderly` | Elderly multiplier. |
| `toxicity_immunosuppressed_multiplier` | Hazard boost for immunosuppressed patients. |
| `toxicity_hospital_multiplier` | Hazard boost while hospitalized. |
| `regional_resistance_threshold_very_high` | Upper threshold that triggers the strongest empiric penalty. |
| `regional_resistance_threshold_high` | Threshold for high-resistance penalty. |
| `regional_resistance_threshold_moderate` | Threshold for moderate penalty. |
| `regional_resistance_penalty_very_high` | Score penalty applied when the region exceeds the very high threshold. |
| `regional_resistance_penalty_high` | Penalty when resistance is high. |
| `regional_resistance_penalty_moderate` | Penalty when resistance is moderate. |

#### Therapy scoring bonuses/penalties
| Parameter | Description |
| --- | --- |
| `targeted_therapy_narrow_spectrum_bonus` | Score bonus when targeted therapy uses narrow-spectrum drugs. |
| `targeted_therapy_broad_spectrum_penalty` | Penalty for using broad-spectrum drugs in targeted mode. |
| `targeted_therapy_ineffective_drug_penalty` | Strong penalty when the chosen targeted drug lacks potency. |
| `effective_potency_threshold_for_targeted_therapy` | Potency cutoff that determines whether a drug counts as effective when targeted. |
| `empiric_therapy_broad_spectrum_bonus` | Bonus for empiric use of broad-spectrum drugs. |
| `empiric_therapy_ineffective_penalty` | Penalty when empiric therapy chooses a low-potency drug. |

#### Sepsis, infection mortality, and background deaths
| Parameter | Description |
| --- | --- |
| `any_r_increase_rate_per_day_when_drug_present` | Daily growth in resistance while treatment pressure is present. |
| `sepsis_minimum_duration_days` | Minimum number of days before sepsis resolution can occur. |
| `sepsis_base_log_odds_of_recovery_per_day` | Baseline log-odds of sepsis recovery per day. |
| `sepsis_log_odds_bacteria_level` | Incremental log-odds change per unit bacteria level. |
| `sepsis_log_odds_in_hospital` | Modifier applied when hospitalized. |
| `sepsis_log_odds_age_infant` | Infant-specific modifier. |
| `sepsis_log_odds_age_child` | Child modifier. |
| `sepsis_log_odds_age_adult` | Adult modifier. |
| `sepsis_log_odds_age_elderly` | Elderly modifier. |
| `sepsis_log_odds_immunosuppressed` | Modifier for immunosuppressed patients. |
| `infection_non_sepsis_base_log_odds` | Baseline non-sepsis infection death log-odds. |
| `infection_non_sepsis_log_odds_per_level` | Increment per bacteria level for non-sepsis mortality. |
| `infection_non_sepsis_log_odds_age_infant` | Infant modifier. |
| `infection_non_sepsis_log_odds_age_child` | Child modifier. |
| `infection_non_sepsis_log_odds_age_adult` | Adult modifier. |
| `infection_non_sepsis_log_odds_age_elderly` | Elderly modifier. |
| `infection_non_sepsis_log_odds_immunosuppressed` | Immunosuppression modifier. |
| `infection_non_sepsis_log_odds_in_hospital` | Hospital modifier. |
| `infection_non_sepsis_minimum_bacteria_level` | Bacteria level required before non-sepsis mortality is considered. |
| `background_mortality_baseline_log_odds` | Baseline all-cause mortality intercept. |
| `mortality_baseline_1930_multiplier` | Early-era multiplier applied near 1930. |
| `mortality_baseline_2035_multiplier` | Future-era multiplier applied near 2035. |
| `mortality_improvement_half_life_years` | Half-life controlling interpolation between baseline multipliers. |
| `log_odds_mortality_per_year_of_age` | Linear age term. |
| `log_odds_mortality_per_year_of_age_squared` | Quadratic age term. |
| `log_odds_mortality_immunosuppressed` | Mortality penalty for immunosuppression. |
| `log_odds_mortality_hospitalized` | Mortality penalty while hospitalized. |
| `base_sepsis_death_risk_per_day` | Base daily risk once in sepsis. |
| `sepsis_age_mortality_multiplier_infant` | Age-specific multipliers on top of the base risk. |
| `sepsis_age_mortality_multiplier_child` | Age-specific multiplier. |
| `sepsis_age_mortality_multiplier_adult` | Age-specific multiplier. |
| `sepsis_age_mortality_multiplier_elderly` | Age-specific multiplier. |
| `sepsis_immunosuppressed_multiplier` | Mortality multiplier when immunosuppressed. |
| `log_odds_sepsis_region_a` | Region-level coefficient (A) used when evaluating sepsis odds. |
| `log_odds_sepsis_region_b` | Region-level coefficient (B). |
| `log_odds_sepsis_onset_immunosuppressed` | Log-odds increase for sepsis onset when immunocompromised (~2× risk at 0.7). |
| `log_odds_sepsis_onset_hospitalized` | Log-odds increase for sepsis onset when hospitalized (~1.6× risk at 0.5). |
| `log_odds_sepsis_onset_not_under_care` | Log-odds increase for sepsis onset when not receiving any antibiotic therapy (~2.7× risk at 1.0). |
| `log_odds_sepsis_infection_duration` | Log-odds contribution per day of infection duration before sepsis onset. |
| `sepsis_death_bacteria_level_coefficient` | Fractional increase in sepsis death risk per unit bacteria level (10% at 0.1). |
| `sepsis_death_duration_coefficient` | Fractional increase in sepsis death risk per day of sepsis duration (2% at 0.02). |
| `sepsis_death_not_under_care_multiplier` | Multiplier on sepsis death risk when not receiving any antibiotic therapy (2× at 2.0). |

#### Microbiome and majority-resistance tracking
| Parameter | Description |
| --- | --- |
| `antibiotic_disruption_decay_half_life_days` | Half-life for microbiome disruption effects from antibiotics. |
| `microbiome_resistance_multiplier_on_acquisition` | Multiplier applied to resistance levels when microbiome carriage is first acquired. |
| `infection_from_microbiome_dampening` | Factor reducing infection risk attributable to microbiome sources. |
| `carriage_duration_log_odds_coefficient` | Slope of the log-odds model used to determine carriage persistence. |
| `carriage_duration_max_log_odds_effect` | Cap for the carriage duration effect. |
| `antibiotic_clearance_log_odds_per_unit_activity` | Incremental effect of antibiotic activity on microbiome clearance log-odds. |
| `carrier_resistance_inheritance_probability` | Probability that a microbiome acquisition inherits the donor's resistance. |
| `majority_r_memory_retention_per_day` | Daily retention factor for the rolling majority-resistance metric. |
| `microbiome_majority_decay_half_life_days` | Half-life for majority microbiome resistance. |
| `microbiome_minority_decay_half_life_days` | Half-life for minority microbiome resistance. |
| `microbiome_majority_promotion_rate_per_day` | Rate at which minority resistance graduates to majority status. |
| `majority_r_window_days` | Observation window for reporting majority resistance. |
| `majority_r_min_total_samples` | Minimum number of samples required before majority statistics are stable. |
| `majority_r_freeze_at_last_positive` | If true, freezes the majority metric once resistance is last observed during the window. |

#### Horizontal gene transfer context multipliers
| Parameter | Description |
| --- | --- |
| `hgt_hospital_multiplier` | Multiplier applied inside [`hgt_context_multiplier`](src/rules/mod.rs#L316-L341) when either participant is hospitalized, reflecting higher plasmid exchange risk in acute-care environments. |
| `hgt_antibiotic_pressure_multiplier` | Boost applied when any antibiotic is active in the host, modeling selective pressure that favors plasmid uptake. |
| `hgt_coinfection_multiplier` | Additional multiplier when both donor and recipient have active infections, mirroring higher DNA release and uptake during symptomatic disease. |
| `hgt_microbiome_only_penalty` | Penalty applied when the transfer attempt is purely microbiome-to-microbiome (no infections on either side), preventing unrealistic amplification outside infections. |

### ImmunodeficiencyParameters

Defined in [src/config.rs#L770-L833](src/config.rs#L770-L833). All values are per-day hazards or age-stratified probabilities.

| Parameter | Description |
| --- | --- |
| `temporary_onset_rate_per_day` | Baseline probability of entering temporary immunosuppression. |
| `temporary_recovery_rate_per_day` | Daily probability of exiting temporary immunosuppression. |
| `chronic_onset_rate_per_day` | Rate of developing chronic immunodeficiency. |
| `chronic_recovery_rate_per_day` | Rate of recovering from chronic immunodeficiency. |
| `chronic_probability_age_bands` | Four-element array giving chronic prevalence at age bands (≤1, ≤18, ≤65, 65+ years). |

### RegionParameters

Defined in [src/config.rs#L779-L834](src/config.rs#L779-L834). Arrays index by `Region`.

| Parameter | Description |
| --- | --- |
| `travel_multiplier` | Multiplier applied to `travel_probability_per_day` for each region. |
| `cessation_multiplier` | Regional multiplier on drug cessation hazards. |
| `mortality_log_odds` | Region-specific adjustments to baseline mortality log-odds. |
| `sepsis_log_odds` | Adjustment to sepsis incidence odds per region. |
| `sepsis_mortality_multiplier` | Region-specific multiplier on sepsis death risk. |
| `testing_multiplier` | Factor applied to diagnostic testing probabilities in each region. |

### SexParameters

Defined in [src/config.rs#L836-L860](src/config.rs#L836-L860).

| Parameter | Description |
| --- | --- |
| `male_log_odds` | Additive log-odds shift applied to baseline mortality for individuals recorded male at birth. |
| `female_log_odds` | Equivalent shift for individuals recorded female at birth. |

### VaccinationParameters

Defined in [src/config.rs#L862-L940](src/config.rs#L862-L940). Arrays are sized `3 × AGE_BUCKETS`.

| Parameter | Description |
| --- | --- |
| `daily_probabilities` | Flattened array of per-vaccine, per-age-group daily vaccination probabilities for pneumococcal, meningococcal, and Hib vaccines. |
| `availability_years` | Year when each vaccine becomes available; prior to this date, uptake remains zero. |

### SyndromeParameters

Defined in [src/config.rs#L1071-L1460](src/config.rs#L1071-L1460) with indices keyed by syndrome id.

Syndromes represent the clinical presentation/infection site:
- **1** = UTI (urinary tract infection)
- **2** = Skin/soft tissue
- **3** = Respiratory
- **4** = Bloodstream (reference compartment)
- **5** = Intra-abdominal
- **6** = CNS (central nervous system)
- **7** = GI (gastrointestinal)
- **8** = Genital
- **9** = Bone/joint
- **10** = Other

| Parameter | Description |
| --- | --- |
| `sepsis_log_odds` | Baseline sepsis log-odds per syndrome. |
| `initiation_multiplier` | Syndrome-specific multiplier on drug initiation hazards. |
| `non_sepsis_mortality_log_odds` | Non-sepsis mortality adjustment per syndrome. |
| `empiric_drug_scores` | Matrix of empiric preference scores per syndrome and drug, used during empiric selection. |
| `bacteria_growth_multiplier` | Syndrome-specific multiplier on bacteria growth rate (e.g., CNS 1.3×, bone 0.85×). |
| `drug_penetration` | Matrix of drug penetration factors `[syndrome_id][drug_idx]` representing fraction of serum concentration achieved at infection site (0.0–1.0). |

#### Syndrome-Specific Drug Penetration

Drug penetration factors account for pharmacokinetic differences at different infection sites. Key defaults:

| Syndrome | Notable Drug Penetration |
| --- | --- |
| **UTI (1)** | Fluoroquinolones/TMP-SMX/nitrofurantoin = 1.0 (urinary concentration) |
| **Skin (2)** | Most drugs 0.80–0.90 (good perfusion) |
| **Respiratory (3)** | Macrolides/FQ = 0.95, aminoglycosides = 0.40 |
| **Bloodstream (4)** | All drugs = 1.0 (reference compartment) |
| **Intra-abdominal (5)** | Metronidazole = 0.90, aminoglycosides = 0.30 (acidic pH) |
| **CNS (6)** | Linezolid/metronidazole = 0.70–0.80, aminoglycosides = 0.05 (blood-brain barrier) |
| **GI (7)** | Oral vancomycin = 0.90 (luminal), metronidazole = 0.95 |
| **Genital (8)** | Fluoroquinolones = 0.90 (prostatic), aminoglycosides = 0.35 |
| **Bone/joint (9)** | Rifampicin = 0.80, linezolid = 0.75, FQ = 0.70, aminoglycosides = 0.25 |
| **Other (10)** | Moderate reduction for aminoglycosides/nitrofurantoin |

The effective drug level at the infection site is calculated as:
```
effective_drug_level = serum_level × penetration_factor
```
where `penetration_factor` represents the fraction of serum concentration achieving therapeutic effect at the protected site.

#### Syndrome Assignment by Bacteria

Each bacteria has clinically-appropriate syndrome probability distributions defined in [src/rules/mod.rs#L4830-L5060](src/rules/mod.rs#L4830-L5060). All 39 bacteria now have explicit entries including:

| Bacteria | Primary Syndrome | Distribution Highlights |
| --- | --- | --- |
| **p_stuartii** | UTI (70%) | Catheter-associated UTI specialist |
| **stenotrophomonas_maltophilia** | Respiratory (50%) | VAP/pneumonia + BSI (32%) |
| **staphylococcus_epidermidis** | Bloodstream (55%) | CLABSI, prosthetic infections |
| **mycoplasma_genitalium** | Genital (85%) | STI causing urethritis/cervicitis |
| **bacteroides_fragilis** | Intra-abdominal (65%) | Peritonitis, abscesses |
| **mdr_mycobacterium_tuberculosis** | Respiratory (82%) | Pulmonary TB + CNS (5%) |

Enteric fever organisms (S. Typhi, S. Paratyphi, iNTS) have corrected distributions emphasizing bloodstream involvement consistent with their systemic nature.

### DrugParameters

Defined in [src/config.rs#L1100-L1194](src/config.rs#L1100-L1194) and indexed by drug.

| Parameter | Description |
| --- | --- |
| `initial_level` | Standardized level assigned on days a standard dose is administered. |
| `double_dose_multiplier` | Factor applied when a double dose is prescribed. |
| `spectrum_breadth` | Relative breadth used when scoring empiric vs targeted therapy. |
| `half_life_days` | Daily decay half-life for active drug levels. |
| `toxicity_death_hazard_per_unit_level` | Hazard contribution per unit of current drug level. |
| `toxicity_reservoir_half_life_days` | Half-life for the toxicity reservoir (cumulative exposure). |
| `microbiome_disruption_log_odds` | Log-odds increment for disrupting microbiome carriage per exposure day. |

### BacteriaParameters

Defined in [src/config.rs#L1196-L1505](src/config.rs#L1196-L1505) for each bacteria entry.

| Parameter | Description |
| --- | --- |
| `acquisition_log_odds_baseline` | Baseline log-odds of acquisition. |
| `log_odds_vaccinated` | Effect of vaccination status on acquisition log-odds. |
| `log_odds_microbiome_present` | Effect of already carrying the microbiome form. |
| `log_odds_hospital_acquired` | Adjustment when infection source is hospital. |
| `microbiome_clearance_probability_per_day` | Baseline probability of clearing microbiome carriage. |
| `environmental_acquisition_proportion` | Share of cases attributed to environmental sources. |
| `initial_infection_level` | Starting bacteria level upon acquisition. |
| `base_bacteria_level_change` | Default daily change in bacteria level without treatment. |
| `max_level` | Maximum allowed bacteria load. |
| `daily_symptom_onset_probability` | Baseline chance of symptoms appearing each day. |
| `symptom_onset_threshold_level` | Bacteria level required before symptoms can occur. |
| `symptom_onset_delay_days` | Delay after infection before symptoms may appear. |
| `symptom_onset_level_multiplier` | Multiplier to the bacteria level when symptoms start. |
| `mechanismless_resistance_reversion_rate` | Daily probability that residual resistance without an active mechanism collapses once drug pressure stops. |
| `microbiome_vs_infection_log_odds` | Log-odds of an exposure resulting in microbiome carriage instead of infection. |
| `drug_cessation_probability` | Override for cessation hazard specific to this pathogen. |
| `treatment_recognition_year` | Optional year when targeted therapy becomes routinely recognized. |
| `sepsis_baseline_log_odds` | Baseline log-odds of developing sepsis once infected. |
| `sepsis_log_odds_infection_level` | Added log-odds per bacteria level unit. |
| `sepsis_log_odds_infection_duration` | Added log-odds per day since infection start. |
| `infection_non_sepsis_mortality_log_odds` | Bacteria-specific non-sepsis mortality modifier. |

### ClearanceParameters

Defined in [src/config.rs#L1506-L1597](src/config.rs#L1506-L1597).

| Parameter | Description |
| --- | --- |
| `base_delay_days` | Default immune arming delay after infection. |
| `base_daily_hazard` | Baseline daily clearance hazard once armed. |
| `per_bacteria_delay_days` | Optional overrides per bacteria for the arming delay. |
| `per_bacteria_hazard_multiplier` | Per-bacteria multipliers on the base hazard. |
| `age_multipliers` | Per-age-bucket multipliers applied to the hazard. |
| `immunodeficient_multiplier` | Factor applied when the host is immunodeficient. |
| `level_reference` | Reference bacteria level used when computing the level modifier. |
| `level_exponent` | Exponent governing how strongly infection level modulates hazard. |

### DrugBacteriaMatrix

Defined in [src/config.rs#L1623-L1698](src/config.rs#L1623-L1698). Indexed by bacteria × drug.

| Parameter | Description |
| --- | --- |
| `potency_when_no_r` | Potency score when the infection is fully susceptible. |
| `initiation_multiplier` | Adjustment to initiation probability for this specific bacteria–drug pair. |
| `resistance_emergence_rate` | Baseline daily emergence rate for resistance under this pairing. |
| `mic_lt2_threshold` | Derived threshold used for MIC-style reporting when potency falls below 2. |

### RegionBacteriaAcquisition

Defined in [src/config.rs#L1700-L1755](src/config.rs#L1700-L1755).

| Parameter | Description |
| --- | --- |
| `values` | Flattened table of acquisition log-odds per region and bacteria, merged with defaults when entries are missing. |

### AgeTables

Defined in [src/config.rs#L1756-L1829](src/config.rs#L1756-L1829).

| Parameter | Description |
| --- | --- |
| `bacteria_age_log_odds` | Per-bacteria, per-age-bucket log-odds adjustments. |
| `region_age_log_odds` | Per-region, per-age-bucket adjustments applied on top of bacteria effects. |
| `default_age_log_odds` | Age-bucket defaults used when bacteria- or region-specific data is unavailable. |

### AgeCategoryParameters

Defined in [src/config.rs#L1831-L1933](src/config.rs#L1831-L1933).

| Parameter | Description |
| --- | --- |
| `bacteria_age_log_odds` | Same as in `AgeTables`, stored for runtime lookups. |
| `region_age_log_odds` | Region-specific log-odds arrays. |
| `bacteria_region_age_log_odds` | Combined bacteria × region × age entries for the most granular tuning. |
| `default_age_log_odds` | Compact array of defaults (prenatal through elderly). |
| `num_bacteria` | Cached bacteria count for index math. |

### HgtMatrix

Defined in [src/config.rs#L1972-L2074](src/config.rs#L1972-L2074).

| Parameter | Description |
| --- | --- |
| `values` | Flattened donor × recipient probability matrix for horizontal gene transfer. |
| `num_bacteria` | Dimension used for indexing into `values`. |

`ParameterStore::builder()` seeds every `hgt_prob_*` key using [`default_hgt_probability`](src/config.rs#L2033-L2074), which groups bacteria into plasmid-compatibility pools and applies realistic structural exclusions ([src/config.rs#L4353-L4369](src/config.rs#L4353-L4369)). Gram-positive, enteric gram-negative, respiratory gram-negative, and anaerobe pools receive different default magnitudes, while spirochetes, *Helicobacter*, and *Mycobacterium* never donate or receive.

At runtime, HGT only triggers when the donor and recipient occupy at least one shared carriage compartment. That check leverages [`bacteria_presence_compartment_mask`](src/rules/mod.rs#L286-L343), which combines the static compartment metadata from [src/simulation/population.rs#L398-L455](src/simulation/population.rs#L398-L455) with per-person infection and microbiome state. The main HGT loop in [src/rules/mod.rs#L4335-L4415](src/rules/mod.rs#L4335-L4415) enforces this overlap before sampling transfers.

Environmental context further scales the transfer probability through [`hgt_context_multiplier`](src/rules/mod.rs#L316-L341): hospitalization, concurrent antibiotic pressure, and donor/recipient co-infections boost the chance via the new global scalars documented below, whereas microbiome-only interactions incur a penalty. This keeps plasmid exchange concentrated in plausible clinical settings while still allowing overrides through the configuration surface.

### ResistanceMechanismParameters

Defined in [src/config.rs#L1984-L2037](src/config.rs#L1984-L2037) and indexed by mechanism enumeration order.

| Parameter | Description |
| --- | --- |
| `emergence_rate` | Baseline probability that a mechanism activates on any given day. |
| `enhancement_multiplier` | Multiplier applied to the potency penalty when the mechanism is present. |
| `reversion_rate` | Probability that the mechanism deactivates when pressure is removed. |

### BacteriaMechanismEmergenceMultipliers

Defined in [src/config.rs#L2039-L2085](src/config.rs#L2039-L2085).

| Parameter | Description |
| --- | --- |
| `values` | Flattened bacteria × mechanism multiplier array applied on top of the global mechanism emergence rate. |
| `num_mechanisms` | Mechanism count used for indexing the flattened vector. |

## Policy Scenarios

The simulation supports branching policy scenarios that fork from the baseline trajectory at a configurable time point, allowing counterfactual comparisons. Policy adjustments are defined in [src/simulation/simulation.rs#L56-L120](src/simulation/simulation.rs#L56-L120) via the `PolicyAdjustments` struct.

### Policy 0: Baseline (Status Quo)

The reference scenario with no interventions. All policy multipliers default to 1.0 or `None`, meaning baseline parameters from `GlobalScalars` and other stores apply without modification. This branch always runs and serves as the comparator for alternate policies.

### Policy 1: Antimicrobial Stewardship Intervention

Models a comprehensive stewardship program with six coordinated interventions:

| Lever | Value | Mechanism | Clinical Rationale |
| --- | --- | --- | --- |
| `drug_selection_temperature` | 0.65× baseline | Lower softmax temperature → more deterministic prescribing | Clinicians follow guidelines more closely; less random prescribing |
| `bacterial_testing_rate_multiplier` | 1.5× | 50% increase in bacterial culture/identification | Better pathogen identification enables targeted therapy |
| `resistance_testing_rate_multiplier` | 1.5× | 50% increase in antimicrobial susceptibility testing | Matched to bacterial testing; AST enables de-escalation |
| `reserve_drug_penalty_multiplier` | 2.0× | Reserve penalty squared (penalty^2.0) | Carbapenems, colistin, linezolid restricted unless failure documented |
| `drug_initiation_rate_multiplier` | 0.85× | 15% reduction in initiation log-odds | Reduces unnecessary prescribing ("antibiotic time-out") |
| `drug_cessation_rate_multiplier` | 1.2× | 20% increase in cessation probability | Shorter courses (evidence supports 5-7 days for most infections) |

**Implementation details:**
- Reserve drug penalty uses exponential scaling: `penalty^multiplier`, so 2.0× squares the penalty making reserve drugs much harder to select empirically
- Initiation reduction applies in log-odds space: `ln(0.85) ≈ -0.16` reduces odds by ~15%
- Cessation multiplier applies to both "active infection" and "no infection" pathways
- Testing rate increases apply multiplicatively to baseline regional testing rates

**Expected effects:**
- Reduced broad-spectrum and reserve antibiotic use
- Faster de-escalation to narrow-spectrum agents
- Shorter treatment courses reducing selection pressure
- Better diagnostic-driven prescribing
- Slower resistance emergence at population level

### Policy 2: AMR Counterfactual (World Without Resistance)

A hypothetical scenario estimating the burden attributable to antimicrobial resistance by simulating a world where resistance never emerged:

| Lever | Value | Effect |
| --- | --- | --- |
| `resistance_emergence_multiplier` | 0.0 | Resistance never emerges de novo |
| `clear_all_resistance_on_branch_start` | true | All existing resistance cleared when branch forks |

**Use case:** Comparing deaths, treatment failures, and healthcare costs between Policy 0 (with resistance) and Policy 2 (without resistance) quantifies the mortality and morbidity burden directly attributable to AMR.

### Adding Custom Policies

New policies can be added by implementing additional constructor functions in `PolicyAdjustments`:

```rust
fn custom_policy() -> Self {
    Self {
        policy_option: 3,
        drug_selection_temperature: Some(0.3),  // Very deterministic
        bacterial_testing_rate_multiplier: Some(2.0),  // Universal testing
        resistance_testing_rate_multiplier: Some(2.0),
        resistance_emergence_multiplier: None,
        clear_all_resistance_on_branch_start: false,
        reserve_drug_penalty_multiplier: Some(3.0),  // Strict reserve restrictions
        drug_initiation_rate_multiplier: Some(0.7),  // 30% reduction
        drug_cessation_rate_multiplier: Some(1.5),   // Much shorter courses
    }
}
```

Then register the policy in `Simulation::new()`:
```rust
let branch_policies = vec![
    PolicyAdjustments::alternate_example(globals),
    PolicyAdjustments::amr_counterfactual(),
    PolicyAdjustments::custom_policy(),  // Add here
];
```

### Policy Output

Each policy branch writes its own summary CSV with `policy_option` column distinguishing branches. The baseline (policy_option=0) summary is always written; alternate branches (policy_option=1, 2, ...) are written to separate files or appended with their option identifier. Downstream analysis can compare trajectories by joining on `time_step` and filtering by `policy_option`.


# AMR Simulation and Analysis Project

This project simulates the spread and treatment of antimicrobial resistance (AMR) in a population using an agent-based model written in Rust, with comprehensive Python-based analysis and visualization tools.

## Project Structure

### Rust Simulation Model (`src/`)
The core simulation engine written in Rust for high-performance modeling of large populations.

**Key Components:**
- `simulation/`: Orchestrates the simulation, manages populations, and aggregates statistics
- `population.rs`: Defines Individual and Population structs with comprehensive attributes
- `rules/`: Core logic for updating individuals each time step (infection, resistance, drug initiation)
- `config/`: Model parameters and lookup functions

### Python Analysis System (`amr_simulation_output_analysis/`)
Modular analysis and visualization system for processing simulation outputs.

**Architecture:**
```
amr_simulation_output_analysis/
├── __init__.py
├── config.py              # Configuration management with granular plot controls
├── data_loader.py         # Centralized data loading and caching
├── utils.py              # Common utilities and helper functions
├── plotting/
│   ├── __init__.py
│   ├── grouped_plots.py   # Main grouped figures (Figures 1-9)
│   ├── detail_plots.py    # Individual detailed visualizations
│   └── utils.py          # Plotting utilities and styling
└── empirical/
    ├── __init__.py
    ├── data_config.py     # Empirical data source configuration
    └── parsers.py         # Data parsing and integration functions
```

## Quick Start

### Running the Simulation
```bash
cargo run --release
```

Configure population size and time steps in `main.rs`:
```rust
let population_size = 100_000;
let time_steps = 20;
```

### Analyzing Results
```python
from amr_simulation_output_analysis import create_all_plots, PlotConfig

# Generate all plots with default configuration
create_all_plots()

# Generate specific plot categories
config = PlotConfig(
    grouped_plots=True,
    detail_plots=False,
    age_specific_plots=True,
    drug_failure_plots=True
)
create_all_plots(config)
```

For complete usage examples, see `amr_analysis.py`.

## Analysis Capabilities

### Grouped Visualizations (Figures 1-9)
- **Figure 1**: Infection resolution patterns by bacteria
- **Figure 2**: Mean activity R analysis by bacteria
- **Figure 3**: Day-7 drug initiation patterns
- **Figure 4**: Syndrome distribution patterns  
- **Figure 5**: Population dynamics tracking
- **Figure 6**: Drug initiation timing analysis
- **Figure 7**: Regional syndrome patterns
- **Figure 8**: Comprehensive resistance tracking
- **Figure 9**: Treatment outcome analysis

### Detailed Analysis Categories
- **Age-specific death rates** (7 visualizations)
- **Regional analysis** (age distribution, death rates, drug failure)
- **Bacteria-specific analysis** (death rates, drug failure, infection resolution)
- **Drug analysis** (usage patterns, resistance emergence, MIC analysis)
- **Resistance mechanisms** (source tracking, emergence patterns)
- **Population health metrics** (mortality, syndrome distribution)

### Configuration Options
The system provides granular control over plot generation:

```python
config = PlotConfig(
    # Main categories
    grouped_plots=True,
    detail_plots=True,
    
    # Specific plot types
    drug_failure_rate_by_bacteria_region=True,
    mean_mic_by_drug_for_each_bacteria=True,
    incidence_of_infection_hospital=True,
    proportion_of_people_taking_each_drug=True,
    # ... 20+ additional granular controls
)
```

## Model Features

### Biological Processes
- **Infection Acquisition**: Contact-based transmission with regional, age, and exposure factors
- **Resistance Emergence**: purely mechanistic ("bottom-up") model using physics-based Gaussian dose response curves.
- **Microbiology**: 39 tracked bacteria with species-specific fitness costs, emergence probabilities, and HGT rates (including "superbug" profiles for *Pseudomonas*, *N. gonorrhoeae*, etc.).
- **Drug Pharmacodynamics**: Multi-drug interactions with realistic decay and site-specific efficacy
- **Testing and Diagnosis**: Laboratory delays, identification accuracy, resistance testing
- **Clearance Dynamics**: Bacteria-specific clearance hazards with age and immunodeficiency modifiers

### Population Dynamics
- **Demographics**: Age, sex, region with realistic population structures
- **Healthcare**: Hospitalization, healthcare access, treatment protocols
- **Mobility**: Inter-regional travel affecting exposure and transmission
- **Mortality**: Background, sepsis-related, and drug toxicity mortality

### Drug and Resistance Modeling
- **Multi-drug Support**: 52 drugs across major antibiotic classes
- **Resistance Mechanisms**: Explicit genetic mechanisms (efflux, target mod, etc.) drive phenotypic resistance (`any_r`).
- **Treatment Protocols**: Empirical and targeted therapy based on testing
- **Microbiome Effects**: Commensal bacteria resistance affecting treatment
- **Pharmacokinetic Modeling**: Syndrome-specific drug penetration and time-to-therapeutic accumulation

## Data Integration

The analysis system integrates multiple empirical data sources:
- **CDC AR Threats 2019**: US antimicrobial resistance surveillance
- **ECDC Data 2023**: European resistance and consumption data
- **IQVIA Sales 2023**: Global antibiotic usage patterns
- **WHO Global Data**: International resistance surveillance
- **GBD Study Data**: Global burden of disease metrics

## Output and Visualization

The system generates comprehensive visualizations:
- **2,000+ individual plots** across all analysis categories
- **Multi-format support**: PNG, PDF, SVG output options
- **Publication-ready figures** with consistent styling
- **Interactive dashboards** (planned enhancement)
- **Statistical summaries** with empirical data overlays

## Development and Customization

### Adding New Analysis
1. Implement new functions in `detail_plots.py`
2. Add configuration flags in `config.py`
3. Update plot dispatcher in `detail_plots.py`
4. Add empirical data sources as needed

### Extending the Simulation
- **Parameters**: Edit `src/config.rs` for model parameters
- **Bacteria/Drugs**: Modify lists in `population.rs` and update `config.rs`
- **Rules**: Update `src/rules/mod.rs` for biological or treatment logic
- **Output**: Extend CSV output in simulation for new metrics

## Archive and Legacy
Previous monolithic analysis scripts have been archived in `archive/legacy_scripts/`:
- `analyze_simulation.py` (6,623-line original script)
- `empirical_enhancement.py` (legacy enhancement tools)
- Development and migration scripts in `archive/development_scripts/`

## Requirements

**Rust**: Latest stable version
**Python**: 3.8+ with packages listed in `requirements.txt`

This project is for research and educational use in antimicrobial resistance modeling. 
