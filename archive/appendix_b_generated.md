## Appendix B ΓÇö Parameter Reference

This appendix is auto-generated from the live Rust configuration. Parameters are organized thematically into resolved tables derived from the internal data structures. All values shown are the effective defaults before any run-level sampling multipliers are applied.

### B.1 Global Scalar Parameters

Scalar parameters that govern cross-cutting model behaviour. Grouped thematically; each row gives the parameter name and its default value.

See: [┬º6.1 Treatment initiation](#61-treatment-initiation-deciding-to-start-antibiotics), [┬º6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use), [┬º6.3 Drug pharmacokinetics](#63-drug-pharmacokinetics), [┬º6.7 Drug toxicity](#67-drug-toxicity), [┬º2.4 Hospitalisation](#24-hospitalisation), [┬º2.5 Travel](#25-travel), [┬º4.3 Sepsis](#43-sepsis), [┬º7.3 Resistance emergence](#73-resistance-emergence), [┬º7.4 Resistance reversion](#74-resistance-reversion-and-fitness-costs), [┬º8 Microbiome and Carriage](#8-microbiome-and-carriage), [┬º9 Horizontal Gene Transfer](#9-horizontal-gene-transfer-hgt), [┬º10 Mortality](#10-mortality).

#### Treatment Initiation (logistic model)

| Parameter | Value |
| --- | ---: |
| antibiotic_initiation_base_log_odds | -5.5 |
| antibiotic_initiation_log_odds_symptomatic_infection | 6 |
| antibiotic_initiation_log_odds_test_identified | 0.92 |
| antibiotic_initiation_log_odds_already_on_drug | 0.18 |
| antibiotic_initiation_log_odds_immunodeficiency | -0.75 |
| antibiotic_initiation_log_odds_sepsis | 6 |
| antibiotic_initiation_log_odds_no_indication | -1.05 |

#### Drug Activity and Cessation

| Parameter | Value |
| --- | ---: |
| drug_activity_to_bacteria_level_multiplier | 0.75 |
| drug_activity_slow_clearance_probability | 0.25 |
| drug_activity_slow_clearance_multiplier | 0.2 |
| double_dose_probability_if_identified_infection | 0.25 |
| random_drug_cessation_probability | 0.0045 |
| random_drug_cessation_probability_if_no_active_infection | 0.15 |
| antibiotic_infection_prevention_efficacy | 0.7 |

#### Drug Selection

| Parameter | Value |
| --- | ---: |
| minimal_potency_threshold_for_drug_selection | 0.15 |
| drug_selection_temperature | 0.55 |
| reserve_drug_score_penalty | 0.005 |

#### Treatment Failure and Restart

| Parameter | Value |
| --- | ---: |
| treatment_failure_assessment_day | 4 |
| treatment_failure_threshold | 0.5 |
| drug_failure_memory_days | 14 |
| restart_window_days | 5 |
| restart_bacteria_level_threshold | 1.5 |
| restart_window_probability | 0.3 |

#### Hospitalization

| Parameter | Value |
| --- | ---: |
| hospitalization_base_log_odds | -10.4 |
| hospitalization_log_odds_per_age_year | 0.02 |
| hospitalization_log_odds_sepsis | 13 |
| hospitalization_log_odds_symptomatic_infection | 9.5 |
| hospitalization_symptomatic_infection_level_threshold | 3 |
| hospital_recovery_rate_per_day | 0.28 |
| hospital_max_days | 30 |
| hospital_prevent_discharge_with_sepsis | 1 |

#### Resistance Emergence and Decay

| Parameter | Value |
| --- | ---: |
| max_resistance_level | 1 |
| resistance_emergence_bacteria_level_multiplier | 9 |
| any_r_emergence_level_on_first_emergence | 0.5 |
| multi_drug_penalty_threshold_num_drugs | 2 |
| resistance_development_inhibition_single_drug | 0.05 |
| resistance_development_inhibition_partial_cross | 0.3 |
| mechanism_assignment_probability_on_any_r_gain | 0.8 |
| community_profile_cache_retention | 0.99 |
| mechanism_reversion_rate_global_multiplier | 1 |
| majority_r_memory_retention_per_day | 0.93 |

#### Microbiome Dynamics

| Parameter | Value |
| --- | ---: |
| microbiome_resistance_transfer_probability_per_day | 1e-4 |
| antibiotic_disruption_decay_half_life_days | 30 |
| microbiome_resistance_multiplier_on_acquisition | 0.5 |
| infection_from_microbiome_dampening | 0.7 |
| carriage_duration_log_odds_coefficient | -0.01 |
| carriage_duration_max_log_odds_effect | -2 |
| antibiotic_clearance_log_odds_per_unit_activity | 0.5 |
| carrier_resistance_inheritance_probability | 0.5 |
| community_resistance_dilution_factor | 0.3 |
| microbiome_majority_decay_half_life_days | 60 |
| microbiome_minority_decay_half_life_days | 18 |
| microbiome_majority_promotion_rate_per_day | 0.02 |

#### De Novo and HGT Multipliers

| Parameter | Value |
| --- | ---: |
| infection_de_novo_multiplier | 3 |
| microbiome_de_novo_multiplier | 1 |
| hgt_multiplier | 1 |

#### Horizontal Gene Transfer Modifiers

| Parameter | Value |
| --- | ---: |
| hgt_hospital_multiplier | 3 |
| hgt_antibiotic_pressure_multiplier | 1.5 |
| hgt_coinfection_multiplier | 1.25 |
| hgt_microbiome_only_penalty | 0.65 |
| hgt_gut_compartment_multiplier | 2 |
| hgt_minority_donor_multiplier | 0.2 |

#### Travel

| Parameter | Value |
| --- | ---: |
| travel_probability_per_day | 5e-5 |

#### Bacteria Growth Age Multipliers

| Parameter | Value |
| --- | ---: |
| bacteria_growth_age_multiplier_infant | 1.3 |
| bacteria_growth_age_multiplier_child | 1 |
| bacteria_growth_age_multiplier_adult | 1 |
| bacteria_growth_age_multiplier_elderly | 1.2 |
| bacteria_growth_immunodeficiency_multiplier | 1.5 |

#### Sepsis Onset

| Parameter | Value |
| --- | ---: |
| sepsis_minimum_duration_days | 1 |
| log_odds_sepsis_onset_immunosuppressed | 0.7 |
| log_odds_sepsis_onset_hospitalized | 0.5 |
| log_odds_sepsis_onset_not_under_care | 1 |
| log_odds_sepsis_onset_region_north_america | -0.5 |
| log_odds_sepsis_onset_region_europe | -0.6 |
| log_odds_sepsis_onset_region_oceania | -0.5 |
| log_odds_sepsis_onset_region_asia | -0.1 |
| log_odds_sepsis_onset_region_south_america | 0 |
| log_odds_sepsis_onset_region_africa | 0.1 |

#### Sepsis Recovery

| Parameter | Value |
| --- | ---: |
| sepsis_base_log_odds_of_recovery_per_day | 0 |
| sepsis_log_odds_bacteria_level | -0.3 |
| sepsis_log_odds_in_hospital | 0.8 |
| sepsis_log_odds_age_infant | -0.5 |
| sepsis_log_odds_age_child | 0.4 |
| sepsis_log_odds_age_adult | 0 |
| sepsis_log_odds_age_elderly | -0.7 |
| sepsis_log_odds_immunosuppressed | -1 |

#### Sepsis Death

| Parameter | Value |
| --- | ---: |
| sepsis_death_base_log_odds | -5 |
| sepsis_death_log_odds_age_infant | 1.1 |
| sepsis_death_log_odds_age_child | -0.7 |
| sepsis_death_log_odds_age_adult | 0 |
| sepsis_death_log_odds_age_elderly | 0.9 |
| sepsis_death_log_odds_immunosuppressed | 1.5 |
| sepsis_death_log_odds_bacteria_level | 0.35 |
| sepsis_death_log_odds_duration | 0.04 |
| sepsis_death_log_odds_early_phase | 0.8 |
| sepsis_death_early_phase_days | 3 |
| sepsis_death_log_odds_not_under_care | 1.4 |

#### Non-Sepsis Infection Mortality

| Parameter | Value |
| --- | ---: |
| infection_non_sepsis_base_log_odds | -9 |
| infection_non_sepsis_log_odds_per_level | 0 |
| infection_non_sepsis_log_odds_age_infant | 0 |
| infection_non_sepsis_log_odds_age_child | 0 |
| infection_non_sepsis_log_odds_age_adult | 0 |
| infection_non_sepsis_log_odds_age_elderly | 0 |
| infection_non_sepsis_log_odds_immunosuppressed | 0 |
| infection_non_sepsis_log_odds_in_hospital | 0 |
| infection_non_sepsis_minimum_bacteria_level | 0.5 |

#### Background Mortality

| Parameter | Value |
| --- | ---: |
| background_mortality_baseline_log_odds | -14.3 |
| mortality_baseline_1930_multiplier | 3 |
| mortality_baseline_2035_multiplier | 1 |
| mortality_improvement_half_life_years | 35 |
| log_odds_mortality_per_year_of_age | 0.055 |
| log_odds_mortality_per_year_of_age_squared | 0.008 |
| log_odds_mortality_immunosuppressed | 0.916 |
| log_odds_mortality_hospitalized | 0.262 |

#### Drug Toxicity

| Parameter | Value |
| --- | ---: |
| default_toxicity_reservoir_half_life_days | 1.5 |
| toxicity_age_multiplier_infant | 1.8 |
| toxicity_age_multiplier_child | 1.2 |
| toxicity_age_multiplier_adult | 1 |
| toxicity_age_multiplier_elderly | 2.2 |
| toxicity_immunosuppressed_multiplier | 2.5 |
| toxicity_hospital_multiplier | 1.3 |
| toxicity_discontinuation_threshold | 1e-5 |
| toxicity_discontinuation_avoidance_days | 30 |

#### Regional Resistance Scoring

| Parameter | Value |
| --- | ---: |
| regional_resistance_threshold_very_high | 0.6 |
| regional_resistance_threshold_high | 0.45 |
| regional_resistance_threshold_moderate | 0.1 |
| regional_resistance_penalty_very_high | 0.3 |
| regional_resistance_penalty_high | 0.5 |
| regional_resistance_penalty_moderate | 0.8 |

#### Therapy Scoring

| Parameter | Value |
| --- | ---: |
| targeted_therapy_narrow_spectrum_bonus | 5 |
| targeted_therapy_broad_spectrum_penalty | 0.1 |
| targeted_therapy_ineffective_drug_penalty | 0.001 |
| effective_potency_threshold_for_targeted_therapy | 0.1 |
| empiric_therapy_broad_spectrum_bonus | 0.85 |
| empiric_therapy_ineffective_penalty | 0.001 |

#### MDR-TB Era Multipliers

| Parameter | Value |
| --- | ---: |
| mdr_tb_pre_antibiotic_era_multiplier | 0 |
| mdr_tb_early_antibiotic_era_multiplier | 0 |
| mdr_tb_modern_era_multiplier | 1 |

### B.2 Drug Properties

Pharmacokinetic and clinical properties for each of the 61 modelled antimicrobial agents. The introduction time step is measured in days from 1 January 1930.

See: [┬º6.3 Drug pharmacokinetics](#63-drug-pharmacokinetics), [┬º6.5 Drug potency matrix](#65-drug-potency-matrix), [┬º6.6 Drug availability](#66-drug-availability-by-region-and-era), [┬º6.7 Drug toxicity](#67-drug-toxicity), [┬º6.8 Antibiotic infection prevention](#68-antibiotic-infection-prevention).

| Drug | Class | Intro (days) | Init level | t┬╜ (days) | 2├ù dose mult | Spectrum | Tox hazard | Tox t┬╜ (days) | Microbiome disrupt |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| sulfanilamide | folate_antagonists | 2555 | 10 | 0.45 | 2 | 3 | 0 | 1.5 | 0.3 |
| penicillin_g | penicillins | 3555 | 10 | 0.04 | 2 | 2 | 5e-10 | 1.5 | 0.3 |
| ampicillin | penicillins | 11315 | 10 | 0.04 | 2 | 3 | 5e-10 | 1.5 | 0.3 |
| amoxicillin | penicillins | 13780 | 10 | 0.04 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| piperacillin | penicillins | 16065 | 10 | 0.04 | 2 | 3 | 5e-10 | 1.5 | 0.3 |
| ticarcillin | penicillins | 14600 | 10 | 0.046 | 2 | 3 | 0 | 1.5 | 0.3 |
| cephalexin | cephalosporins_1_2 | 14605 | 10 | 0.04 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| cefazolin | cephalosporins_1_2 | 15700 | 10 | 0.08 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| cefuroxime | cephalosporins_1_2 | 17525 | 10 | 0.05 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| ceftriaxone | cephalosporins_3_4 | 19715 | 10 | 0.33 | 2 | 4 | 5e-10 | 1.5 | 0.3 |
| ceftazidime | cephalosporins_3_4 | 20080 | 10 | 0.08 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| cefepime | cephalosporins_3_4 | 24195 | 10 | 0.08 | 2 | 4 | 1e-9 | 1.5 | 0.3 |
| ceftaroline | cephalosporins_3_4 | 29305 | 10 | 0.11 | 2 | 3 | 5e-10 | 1.5 | 0.3 |
| ceftolozane_tazobactam | cephalosporins_3_4 | 30295 | 10 | 0.125 | 2 | 3 | 0 | 1.5 | 0.3 |
| cefiderocol | unknown | 33510 | 10 | 0.1 | 2 | 3 | 0 | 1.5 | 0.3 |
| meropenem | carbapenems | 24195 | 10 | 0.04 | 2 | 5 | 6e-10 | 1.5 | 0.3 |
| imipenem_c | carbapenems | 20080 | 10 | 0.04 | 2 | 3 | 1e-9 | 1.5 | 0.3 |
| ertapenem | carbapenems | 25920 | 10 | 0.17 | 2 | 3 | 6e-10 | 1.5 | 0.3 |
| aztreonam | cephalosporins_3_4 | 20445 | 10 | 0.08 | 2 | 3 | 0 | 1.5 | 0.3 |
| erythromycin | macrolides | 8025 | 10 | 0.08 | 2 | 3 | 2e-9 | 1.5 | 0.3 |
| azithromycin | macrolides | 22260 | 10 | 2.8 | 2 | 4 | 1.5e-9 | 1.5 | 0.3 |
| clarithromycin | macrolides | 21895 | 10 | 0.25 | 2 | 3 | 1.5e-9 | 1.5 | 0.3 |
| clindamycin | macrolides | 13870 | 10 | 0.125 | 2 | 3 | 1e-9 | 1.5 | 0.3 |
| gentamicin | aminoglycosides | 12045 | 10 | 0.08 | 2 | 3 | 1.5e-8 | 1.5 | 0.3 |
| tobramycin | aminoglycosides | 16325 | 10 | 0.08 | 2 | 3 | 1.3e-8 | 1.5 | 0.3 |
| amikacin | aminoglycosides | 16690 | 10 | 0.08 | 2 | 3 | 1.7e-8 | 1.5 | 0.3 |
| ciprofloxacin | fluoroquinolones | 20805 | 10 | 0.17 | 2 | 4.5 | 3e-9 | 1.5 | 0.3 |
| levofloxacin | fluoroquinolones | 24195 | 10 | 0.33 | 2 | 3 | 3e-9 | 1.5 | 0.3 |
| moxifloxacin | fluoroquinolones | 25290 | 10 | 0.5 | 2 | 3 | 5e-9 | 1.5 | 0.3 |
| ofloxacin | fluoroquinolones | 21895 | 10 | 0.25 | 2 | 3 | 3e-9 | 1.5 | 0.3 |
| tetracycline | tetracyclines | 6575 | 10 | 0.33 | 2 | 3 | 1e-9 | 1.5 | 0.3 |
| doxycycline | tetracyclines | 13505 | 10 | 0.75 | 2 | 3 | 1e-9 | 1.5 | 0.3 |
| minocycline | tetracyclines | 14965 | 10 | 0.67 | 2 | 3 | 1.5e-9 | 1.5 | 0.3 |
| tigecycline | unknown | 28040 | 10 | 1.75 | 2 | 3 | 0 | 1.5 | 0.3 |
| vancomycin | glycopeptides | 10215 | 10 | 0.25 | 2 | 2.5 | 6e-9 | 1.5 | 0.3 |
| teicoplanin | lipoglycopeptides | 21170 | 10 | 3.5 | 2 | 3 | 0 | 1.5 | 0.3 |
| dalbavancin | lipoglycopeptides | 30660 | 10 | 10 | 2 | 3 | 0 | 1.5 | 0.3 |
| linezolid | oxazolidinones | 25550 | 10 | 0.21 | 2 | 2 | 8e-9 | 1.5 | 0.3 |
| tedizolid | oxazolidinones | 30660 | 10 | 0.5 | 2 | 3 | 4e-9 | 1.5 | 0.3 |
| daptomycin | unknown | 27375 | 10 | 0.33 | 2 | 3 | 0 | 1.5 | 0.3 |
| quinu_dalfo | unknown | 25290 | 10 | 0.5 | 2 | 3 | 0 | 1.5 | 0.3 |
| trim_sulf | folate_antagonists | 13870 | 10 | 0.5 | 2 | 3.5 | 2e-9 | 1.5 | 0.3 |
| chloramphenicol | unknown | 6935 | 10 | 0.125 | 2 | 3 | 1e-8 | 1.5 | 0.3 |
| nitrofurantoin | unknown | 8395 | 10 | 0.017 | 2 | 3 | 3e-9 | 1.5 | 0.3 |
| fosfomycin | unknown | 10590 | 10 | 0.15 | 2 | 3 | 0 | 1.5 | 0.3 |
| retapamulin | unknown | 28405 | 10 | 0.25 | 2 | 3 | 0 | 1.5 | 0.3 |
| fusidic_a | unknown | 11680 | 10 | 0.375 | 2 | 3 | 0 | 1.5 | 0.3 |
| metronidazole | nitroimidazoles | 10965 | 10 | 0.33 | 2 | 3 | 2e-9 | 1.5 | 0.3 |
| fidaxomicin | unknown | 29565 | 10 | 0.5 | 2 | 3 | 0 | 1.5 | 0.3 |
| furazolidone | unknown | 9125 | 10 | 0.25 | 2 | 3 | 0 | 1.5 | 0.3 |
| rifampicin | unknown | 13140 | 10 | 0.25 | 2 | 3 | 4e-9 | 1.5 | 0.3 |
| amoxicillin_clavulanate | penicillins | 16425 | 10 | 0.04 | 2 | 3 | 6e-10 | 1.5 | 0.3 |
| piperacillin_tazobactam | penicillins | 19715 | 10 | 0.04 | 2 | 3 | 6e-10 | 1.5 | 0.3 |
| ampicillin_sulbactam | penicillins | 18250 | 10 | 0.04 | 2 | 3 | 0 | 1.5 | 0.3 |
| ticarcillin_clavulanate | penicillins | 18250 | 10 | 0.046 | 2 | 3 | 0 | 1.5 | 0.3 |
| ceftazidime_avibactam | cephalosporins_3_4 | 27740 | 10 | 0.08 | 2 | 3 | 0 | 1.5 | 0.3 |
| meropenem_vaborbactam | carbapenems | 32045 | 10 | 0.04 | 2 | 3 | 0 | 1.5 | 0.3 |
| colistin | polymyxins | 8020 | 10 | 0.08 | 2 | 4 | 2.5e-8 | 1.5 | 0.3 |
| flucloxacillin | penicillins | 14600 | 10 | 0.04 | 2 | 1.6 | 1e-8 | 1.5 | 0.3 |
| aztreonam_avibactam | cephalosporins_3_4 | 34675 | 10 | 0.08 | 2 | 3 | 0 | 1.5 | 0.3 |
| cefixime | cephalosporins_3_4 | 21535 | 10 | 0.125 | 2 | 2.8 | 5e-10 | 1.5 | 0.3 |

### B.3 Bacteria Properties

Per-bacteria parameters governing acquisition, growth, symptom onset, and clinical outcomes for each of the 42 bacterial species.

See: [┬º3.1 Community acquisition](#31-community-acquisition), [┬º4.2 Infection dynamics](#42-infection-dynamics), [┬º4.3 Sepsis](#43-sepsis), [┬º4.4 Natural clearance](#44-natural-clearance), [┬º8.1 Carriage compartments](#81-carriage-compartments).

| Bacteria | Acq log-odds | Init level | ╬ö level/day | Max level | Microb clr/day | Microb vs inf | Drug cess prob | Sx threshold | Sx delay (d) | Sepsis log-odds | Mech-less rev rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| acinetobacter_baumannii | -18.7 | 0.01 | 0.55 | 5 | 0.1 | 8 | 0.0075 | 0.5 | 1 | -7.7 | 4e-4 |
| citrobacter_spp. | -16.3 | 0.01 | 0.5 | 5 | 0.08 | 11.3 | 0.0045 | 0.5 | 1 | -9.2 | 4e-4 |
| enterobacter_spp. | -16.8 | 0.01 | 0.5 | 5 | 0.07 | 11.5 | 0.0045 | 0.5 | 1 | -7.7 | 4e-4 |
| enterococcus_faecalis | -17.1 | 0.01 | 0.48 | 5 | 0.008 | 11.5 | 0.0075 | 0.5 | 1 | -7.5 | 4e-4 |
| enterococcus_faecium | -17.6 | 0.01 | 0.48 | 5 | 0.06 | 13 | 0.0075 | 0.5 | 1 | -7 | 4e-4 |
| escherichia_coli | -11.6 | 0.01 | 0.5 | 5 | 0.005 | 6 | 0.025 | 0.5 | 1 | -9.5 | 4e-4 |
| klebsiella_pneumoniae | -14.7 | 0.01 | 0.52 | 5 | 0.03 | 9.2 | 0.0075 | 0.5 | 1 | -7.5 | 4e-4 |
| morganella_spp. | -17.2 | 0.01 | 0.48 | 5 | 0.1 | 10 | 0.0045 | 0.5 | 1 | -7.8 | 4e-4 |
| proteus_spp. | -16.1 | 0.01 | 0.5 | 5 | 0.08 | 8.5 | 0.0045 | 0.5 | 1 | -7.8 | 4e-4 |
| serratia_spp. | -17.3 | 0.01 | 0.48 | 5 | 0.1 | 10 | 0.0045 | 0.5 | 1 | -8 | 4e-4 |
| p_stuartii | -17.5 | 0.01 | 0.5 | 5 | 0.09 | 8.7 | 0.0045 | 0.75 | 1 | -14 | 4e-4 |
| pseudomonas_aeruginosa | -16 | 0.01 | 0.55 | 5 | 0.12 | 9 | 0.0075 | 0.8 | 1 | -6.5 | 4e-4 |
| stenotrophomonas_maltophilia | -19 | 0.01 | 0.45 | 5 | 0.06 | 7 | 0.0045 | 0.9 | 2.5 | -8 | 4e-4 |
| staphylococcus_aureus | -12.9 | 0.01 | 0.6 | 5 | 0.05 | 7.5 | 0.015 | 0.5 | 1 | -7.3 | 4e-4 |
| staphylococcus_epidermidis | -16 | 0.01 | 0.35 | 4 | 0.015 | 13.5 | 0.0045 | 1 | 3 | -8 | 4e-4 |
| streptococcus_pneumoniae | -12.6 | 0.01 | 0.6 | 5 | 0.05 | 7 | 0.015 | 0.5 | 1 | -10.5 | 4e-4 |
| salmonella_enterica_serovar_typhi | -17.3 | 0.01 | 0.45 | 5 | 0.003 | -7 | 0.0045 | 0.5 | 1 | -8 | 4e-4 |
| salmonella_enterica_serovar_paratyphi_a | -16.8 | 0.01 | 0.45 | 5 | 0.15 | -0.5 | 0.0045 | 0.5 | 1 | -9 | 4e-4 |
| invasive_non-typhoidal_salmonella_spp. | -17.8 | 0.01 | 0.5 | 5 | 0.12 | 3.2 | 0.0045 | 0.5 | 1 | -9.2 | 4e-4 |
| shigella_spp. | -12.6 | 0.01 | 0.55 | 5 | 0.15 | -0.5 | 0.0045 | 0.5 | 1 | -12 | 4e-4 |
| neisseria_gonorrhoeae | -13.5 | 0.01 | 0.55 | 5 | 0.2 | 3 | 0.0045 | 0.5 | 1 | -23 | 4e-4 |
| streptococcus_pyogenes | -14.4 | 0.01 | 0.7 | 5 | 0.08 | 8 | 0.015 | 0.5 | 1 | -6.5 | 4e-4 |
| streptococcus_agalactiae | -15.9 | 0.01 | 0.52 | 5 | 0.06 | 11 | 0.015 | 0.5 | 1 | -7 | 4e-4 |
| haemophilus_influenzae | -17.4 | 0.01 | 0.55 | 5 | 0.06 | 12.5 | 0.015 | 0.5 | 1 | -9.2 | 4e-4 |
| chlamydia_trachomatis | -12.8 | 0.01 | 0.25 | 5 | 0.2 | 4.5 | 0.007 | 0.8 | 1 | -19 | 4e-4 |
| mycoplasma_genitalium | -12.1 | 0.01 | 0.28 | 5 | 0.18 | 4.7 | 0.0045 | 0.9 | 5 | -14 | 4e-4 |
| vibrio_cholerae | -18.65 | 0.01 | 0.7 | 5 | 0.15 | 1 | 0.025 | 0.5 | 1 | -9 | 4e-4 |
| neisseria_meningitidis | -18.5 | 0.01 | 0.65 | 5 | 0.05 | 10.5 | 0.01 | 3 | 1 | -7.9 | 4e-4 |
| listeria_monocytogenes | -19 | 0.01 | 0.25 | 5 | 0.1 | 12.5 | 0.0045 | 0.5 | 1 | -8 | 4e-4 |
| clostridioides_difficile | -15.15 | 0.01 | 0.55 | 5 | 0.02 | 7 | 0.005 | 0.5 | 1 | -11 | 4e-4 |
| bacteroides_fragilis | -15.1 | 0.01 | 0.42 | 5 | 0.004 | 11 | 0.0045 | 1.2 | 2 | -14 | 4e-4 |
| campylobacter_jejuni | -13 | 0.01 | 0.52 | 5 | 0.12 | 2.5 | 0.015 | 0.5 | 1 | -20 | 4e-4 |
| enterobacter_cloacae | -16.8 | 0.01 | 0.5 | 5 | 0.04 | 13 | 0.0045 | 0.5 | 1 | -7.8 | 4e-4 |
| yersinia_enterocolitica | -16.6 | 0.01 | 0.45 | 5 | 0.25 | 5.5 | 0.0045 | 0.5 | 1 | -9.5 | 4e-4 |
| moraxella_catarrhalis | -14.6 | 0.01 | 0.55 | 5 | 0.05 | 11 | 0.0045 | 2 | 1 | -10.8 | 4e-4 |
| treponema_pallidum | -12.7 | 0.01 | 0.18 | 5 | 0.35 | 5.5 | 0.0045 | 0.6 | 1 | -11 | 4e-4 |
| bordetella_pertussis | -12.85 | 0.01 | 0.42 | 5 | 0.2 | 2.5 | 0.0075 | 0.5 | 1 | -11 | 4e-4 |
| helicobacter_pylori | -13.5 | 0.01 | 0.2 | 5 | 0.001 | 6.65 | 0.005 | 1.5 | 30 | -250 | 4e-4 |
| mdr_mycobacterium_tuberculosis | -16.9 | 0.01 | 0.15 | 5 | 0.0015 | -1 | 6e-4 | 2 | 1 | -38 | 4e-4 |
| mycoplasma_pneumoniae | -12 | 0.01 | 0.35 | 5 | 0.01 | 0.5 | 0.015 | 0.5 | 1 | -14 | 4e-4 |
| legionella_pneumophila | -15.5 | 0.01 | 0.55 | 5 | 0.01 | -2 | 0.0085 | 0.5 | 1 | -14 | 4e-4 |
| burkholderia_cepacia_complex | -20 | 0.01 | 0.45 | 5 | 0.01 | 0.5 | 0.0075 | 0.5 | 1 | -14 | 4e-4 |

### B.4 DrugΓÇôBacteria Potency Matrix

Baseline potency (MIC-derived effectiveness when no resistance is present) and initiation multiplier (stewardship weighting for drug selection) for each drugΓÇôbacteria pair. 42 bacteria ├ù 61 drugs = 2562 entries.

See: [┬º6.5 Drug potency matrix](#65-drug-potency-matrix), [┬º6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use).

| Bacteria | Drug | Potency (no R) | Init multiplier |
| --- | ---: | ---: | ---: |
| acinetobacter_baumannii | sulfanilamide | 0.1 | 0.02 |
| acinetobacter_baumannii | penicillin_g | 0.05 | 0.01 |
| acinetobacter_baumannii | ampicillin | 0.05 | 0.01 |
| acinetobacter_baumannii | amoxicillin | 0.05 | 0.01 |
| acinetobacter_baumannii | piperacillin | 0.6 | 1 |
| acinetobacter_baumannii | ticarcillin | 0.5 | 1 |
| acinetobacter_baumannii | cephalexin | 0.05 | 0.3 |
| acinetobacter_baumannii | cefazolin | 0.05 | 0.3 |
| acinetobacter_baumannii | cefuroxime | 0.1 | 0.3 |
| acinetobacter_baumannii | ceftriaxone | 0.1 | 0.2 |
| acinetobacter_baumannii | ceftazidime | 0.6 | 2 |
| acinetobacter_baumannii | cefepime | 0.7 | 2 |
| acinetobacter_baumannii | ceftaroline | 0.1 | 0.002 |
| acinetobacter_baumannii | ceftolozane_tazobactam | 0.1 | 1 |
| acinetobacter_baumannii | cefiderocol | 0.1 | 1 |
| acinetobacter_baumannii | meropenem | 0.85 | 6 |
| acinetobacter_baumannii | imipenem_c | 0.8 | 4 |
| acinetobacter_baumannii | ertapenem | 0.1 | 0.005 |
| acinetobacter_baumannii | aztreonam | 0.1 | 0.003 |
| acinetobacter_baumannii | gentamicin | 0.75 | 10 |
| acinetobacter_baumannii | tobramycin | 0.7 | 10 |
| acinetobacter_baumannii | amikacin | 0.8 | 15 |
| acinetobacter_baumannii | ciprofloxacin | 0.7 | 1 |
| acinetobacter_baumannii | levofloxacin | 0.7 | 1 |
| acinetobacter_baumannii | moxifloxacin | 0.6 | 1 |
| acinetobacter_baumannii | ofloxacin | 0.6 | 1 |
| acinetobacter_baumannii | tetracycline | 0.6 | 0.25 |
| acinetobacter_baumannii | doxycycline | 0.7 | 0.25 |
| acinetobacter_baumannii | minocycline | 0.8 | 0.25 |
| acinetobacter_baumannii | tigecycline | 0.1 | 1 |
| acinetobacter_baumannii | dalbavancin | 0 | 0.005 |
| acinetobacter_baumannii | linezolid | 0 | 0.005 |
| acinetobacter_baumannii | tedizolid | 0 | 0.005 |
| acinetobacter_baumannii | daptomycin | 0.1 | 1 |
| acinetobacter_baumannii | quinu_dalfo | 0 | 0.005 |
| acinetobacter_baumannii | trim_sulf | 0.6 | 0.04 |
| acinetobacter_baumannii | chloramphenicol | 0.7 | 1 |
| acinetobacter_baumannii | nitrofurantoin | 0.1 | 1 |
| acinetobacter_baumannii | fosfomycin | 0.1 | 1 |
| acinetobacter_baumannii | fidaxomicin | 0.1 | 1 |
| acinetobacter_baumannii | furazolidone | 0.1 | 1 |
| acinetobacter_baumannii | rifampicin | 0.6 | 1 |
| acinetobacter_baumannii | amoxicillin_clavulanate | 0.05 | 1 |
| acinetobacter_baumannii | piperacillin_tazobactam | 0.7 | 8 |
| acinetobacter_baumannii | ampicillin_sulbactam | 0.7 | 1 |
| acinetobacter_baumannii | ticarcillin_clavulanate | 0.6 | 3 |
| acinetobacter_baumannii | ceftazidime_avibactam | 0.7 | 0.005 |
| acinetobacter_baumannii | meropenem_vaborbactam | 0.8 | 0.005 |
| acinetobacter_baumannii | colistin | 0.9 | 0.005 |
| acinetobacter_baumannii | flucloxacillin | 0.01 | 1 |
| acinetobacter_baumannii | aztreonam_avibactam | 0.1 | 0.04 |
| acinetobacter_baumannii | cefixime | 0.1 | 0.2 |
| citrobacter_spp. | sulfanilamide | 0.5 | 0.02 |
| citrobacter_spp. | penicillin_g | 0.1 | 1 |
| citrobacter_spp. | ampicillin | 0.1 | 1 |
| citrobacter_spp. | amoxicillin | 0.1 | 1 |
| citrobacter_spp. | piperacillin | 0.8 | 1 |
| citrobacter_spp. | ticarcillin | 0.75 | 1 |
| citrobacter_spp. | cephalexin | 0.1 | 0.3 |
| citrobacter_spp. | cefazolin | 0.1 | 0.3 |
| citrobacter_spp. | cefuroxime | 0.8 | 0.3 |
| citrobacter_spp. | ceftriaxone | 0.85 | 0.2 |
| citrobacter_spp. | ceftazidime | 0.8 | 0.2 |
| citrobacter_spp. | cefepime | 0.9 | 0.35 |
| citrobacter_spp. | ceftaroline | 0.6 | 0.002 |
| citrobacter_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| citrobacter_spp. | cefiderocol | 0.8 | 1 |
| citrobacter_spp. | meropenem | 0.95 | 4 |
| citrobacter_spp. | imipenem_c | 0.95 | 0.005 |
| citrobacter_spp. | ertapenem | 0.9 | 3 |
| citrobacter_spp. | aztreonam | 0.85 | 0.003 |
| citrobacter_spp. | erythromycin | 0.1 | 1 |
| citrobacter_spp. | azithromycin | 0.1 | 1 |
| citrobacter_spp. | clarithromycin | 0.1 | 1 |
| citrobacter_spp. | clindamycin | 0.1 | 1 |
| citrobacter_spp. | gentamicin | 0.85 | 10 |
| citrobacter_spp. | tobramycin | 0.8 | 1 |
| citrobacter_spp. | amikacin | 0.9 | 1 |
| citrobacter_spp. | ciprofloxacin | 0.9 | 1 |
| citrobacter_spp. | levofloxacin | 0.85 | 1 |
| citrobacter_spp. | moxifloxacin | 0.7 | 1 |
| citrobacter_spp. | ofloxacin | 0.8 | 1 |
| citrobacter_spp. | tetracycline | 0.8 | 0.25 |
| citrobacter_spp. | doxycycline | 0.85 | 0.25 |
| citrobacter_spp. | minocycline | 0.85 | 0.25 |
| citrobacter_spp. | tigecycline | 0.1 | 1 |
| citrobacter_spp. | vancomycin | 0.1 | 1 |
| citrobacter_spp. | teicoplanin | 0.1 | 1 |
| citrobacter_spp. | dalbavancin | 0.1 | 0.005 |
| citrobacter_spp. | linezolid | 0.1 | 0.005 |
| citrobacter_spp. | tedizolid | 0.1 | 0.005 |
| citrobacter_spp. | daptomycin | 0.1 | 1 |
| citrobacter_spp. | quinu_dalfo | 0.1 | 0.005 |
| citrobacter_spp. | trim_sulf | 0.9 | 0.04 |
| citrobacter_spp. | chloramphenicol | 0.85 | 1 |
| citrobacter_spp. | nitrofurantoin | 0.8 | 1 |
| citrobacter_spp. | fosfomycin | 0.1 | 1 |
| citrobacter_spp. | retapamulin | 0.05 | 1 |
| citrobacter_spp. | fusidic_a | 0.05 | 1 |
| citrobacter_spp. | metronidazole | 0.05 | 1 |
| citrobacter_spp. | fidaxomicin | 0.1 | 1 |
| citrobacter_spp. | furazolidone | 0.1 | 1 |
| citrobacter_spp. | rifampicin | 0.7 | 1 |
| citrobacter_spp. | amoxicillin_clavulanate | 0.9 | 6 |
| citrobacter_spp. | piperacillin_tazobactam | 0.9 | 8 |
| citrobacter_spp. | ampicillin_sulbactam | 0.85 | 1 |
| citrobacter_spp. | ticarcillin_clavulanate | 0.8 | 1 |
| citrobacter_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| citrobacter_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| citrobacter_spp. | colistin | 0.7 | 0.005 |
| citrobacter_spp. | flucloxacillin | 0.01 | 1 |
| citrobacter_spp. | aztreonam_avibactam | 1 | 0.003 |
| citrobacter_spp. | cefixime | 0.8 | 0.2 |
| enterobacter_spp. | sulfanilamide | 0.5 | 0.02 |
| enterobacter_spp. | penicillin_g | 0.1 | 1 |
| enterobacter_spp. | ampicillin | 0.1 | 1 |
| enterobacter_spp. | amoxicillin | 0.1 | 1 |
| enterobacter_spp. | piperacillin | 0.75 | 1 |
| enterobacter_spp. | ticarcillin | 0.7 | 1 |
| enterobacter_spp. | cephalexin | 0.1 | 0.3 |
| enterobacter_spp. | cefazolin | 0.1 | 0.3 |
| enterobacter_spp. | cefuroxime | 0.6 | 0.3 |
| enterobacter_spp. | ceftriaxone | 0.5 | 0.2 |
| enterobacter_spp. | ceftazidime | 0.8 | 0.2 |
| enterobacter_spp. | cefepime | 0.85 | 2.5 |
| enterobacter_spp. | ceftaroline | 0.4 | 0.002 |
| enterobacter_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| enterobacter_spp. | cefiderocol | 0.8 | 1 |
| enterobacter_spp. | meropenem | 0.95 | 5 |
| enterobacter_spp. | imipenem_c | 0.95 | 3 |
| enterobacter_spp. | ertapenem | 0.9 | 3 |
| enterobacter_spp. | aztreonam | 0.8 | 0.003 |
| enterobacter_spp. | gentamicin | 0.85 | 10 |
| enterobacter_spp. | tobramycin | 0.8 | 1 |
| enterobacter_spp. | amikacin | 0.9 | 8 |
| enterobacter_spp. | ciprofloxacin | 0.9 | 1 |
| enterobacter_spp. | levofloxacin | 0.85 | 1 |
| enterobacter_spp. | moxifloxacin | 0.7 | 1 |
| enterobacter_spp. | ofloxacin | 0.8 | 1 |
| enterobacter_spp. | tetracycline | 0.8 | 0.25 |
| enterobacter_spp. | doxycycline | 0.85 | 0.25 |
| enterobacter_spp. | minocycline | 0.85 | 0.25 |
| enterobacter_spp. | tigecycline | 0.1 | 1 |
| enterobacter_spp. | dalbavancin | 0 | 0.005 |
| enterobacter_spp. | linezolid | 0 | 0.005 |
| enterobacter_spp. | tedizolid | 0 | 0.005 |
| enterobacter_spp. | daptomycin | 0.1 | 1 |
| enterobacter_spp. | quinu_dalfo | 0 | 0.005 |
| enterobacter_spp. | trim_sulf | 0.85 | 0.04 |
| enterobacter_spp. | chloramphenicol | 0.8 | 1 |
| enterobacter_spp. | nitrofurantoin | 0.7 | 1 |
| enterobacter_spp. | fosfomycin | 0.1 | 1 |
| enterobacter_spp. | fidaxomicin | 0.1 | 1 |
| enterobacter_spp. | furazolidone | 0.1 | 1 |
| enterobacter_spp. | rifampicin | 0.6 | 1 |
| enterobacter_spp. | amoxicillin_clavulanate | 0.7 | 6 |
| enterobacter_spp. | piperacillin_tazobactam | 0.85 | 8 |
| enterobacter_spp. | ampicillin_sulbactam | 0.7 | 1 |
| enterobacter_spp. | ticarcillin_clavulanate | 0.8 | 1 |
| enterobacter_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| enterobacter_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| enterobacter_spp. | colistin | 0.7 | 0.005 |
| enterobacter_spp. | flucloxacillin | 0.01 | 1 |
| enterobacter_spp. | aztreonam_avibactam | 1 | 0.003 |
| enterobacter_spp. | cefixime | 0.8 | 0.2 |
| enterococcus_faecalis | sulfanilamide | 0.1 | 0.02 |
| enterococcus_faecalis | penicillin_g | 0.8 | 1 |
| enterococcus_faecalis | ampicillin | 0.9 | 6 |
| enterococcus_faecalis | amoxicillin | 0.9 | 1 |
| enterococcus_faecalis | piperacillin | 0.75 | 1 |
| enterococcus_faecalis | ticarcillin | 0.7 | 1 |
| enterococcus_faecalis | cephalexin | 0.1 | 0.3 |
| enterococcus_faecalis | cefazolin | 0.1 | 0.3 |
| enterococcus_faecalis | cefuroxime | 0.1 | 0.3 |
| enterococcus_faecalis | ceftriaxone | 0.1 | 0.2 |
| enterococcus_faecalis | ceftazidime | 0.1 | 0.2 |
| enterococcus_faecalis | cefepime | 0.1 | 0.35 |
| enterococcus_faecalis | ceftaroline | 0.1 | 0.002 |
| enterococcus_faecalis | ceftolozane_tazobactam | 0.05 | 1 |
| enterococcus_faecalis | cefiderocol | 0.05 | 1 |
| enterococcus_faecalis | meropenem | 0.7 | 0.005 |
| enterococcus_faecalis | imipenem_c | 0.7 | 0.005 |
| enterococcus_faecalis | ertapenem | 0.1 | 0.005 |
| enterococcus_faecalis | aztreonam | 0 | 0.003 |
| enterococcus_faecalis | erythromycin | 0.7 | 0.01 |
| enterococcus_faecalis | azithromycin | 0.7 | 0.01 |
| enterococcus_faecalis | clarithromycin | 0.7 | 1 |
| enterococcus_faecalis | clindamycin | 0.7 | 1 |
| enterococcus_faecalis | gentamicin | 0.1 | 20 |
| enterococcus_faecalis | tobramycin | 0.1 | 1 |
| enterococcus_faecalis | amikacin | 0.1 | 1 |
| enterococcus_faecalis | ciprofloxacin | 0.7 | 1 |
| enterococcus_faecalis | levofloxacin | 0.7 | 1 |
| enterococcus_faecalis | moxifloxacin | 0.7 | 1 |
| enterococcus_faecalis | ofloxacin | 0.7 | 1 |
| enterococcus_faecalis | tetracycline | 0.8 | 0.25 |
| enterococcus_faecalis | doxycycline | 0.8 | 0.25 |
| enterococcus_faecalis | minocycline | 0.85 | 0.25 |
| enterococcus_faecalis | tigecycline | 0.1 | 1 |
| enterococcus_faecalis | vancomycin | 0.95 | 4 |
| enterococcus_faecalis | teicoplanin | 0.9 | 5 |
| enterococcus_faecalis | dalbavancin | 0.9 | 0.005 |
| enterococcus_faecalis | linezolid | 0.9 | 2 |
| enterococcus_faecalis | tedizolid | 0.9 | 0.005 |
| enterococcus_faecalis | daptomycin | 0.1 | 3 |
| enterococcus_faecalis | quinu_dalfo | 0.1 | 0.005 |
| enterococcus_faecalis | trim_sulf | 0.1 | 0.04 |
| enterococcus_faecalis | chloramphenicol | 0.8 | 1 |
| enterococcus_faecalis | nitrofurantoin | 0.9 | 20 |
| enterococcus_faecalis | fosfomycin | 0.1 | 15 |
| enterococcus_faecalis | retapamulin | 0.1 | 1 |
| enterococcus_faecalis | fusidic_a | 0.1 | 1 |
| enterococcus_faecalis | metronidazole | 0.1 | 1 |
| enterococcus_faecalis | fidaxomicin | 0.1 | 1 |
| enterococcus_faecalis | furazolidone | 0.1 | 1 |
| enterococcus_faecalis | rifampicin | 0.1 | 1 |
| enterococcus_faecalis | amoxicillin_clavulanate | 0.9 | 12 |
| enterococcus_faecalis | piperacillin_tazobactam | 0.75 | 1 |
| enterococcus_faecalis | ampicillin_sulbactam | 0.9 | 1 |
| enterococcus_faecalis | ticarcillin_clavulanate | 0.7 | 1 |
| enterococcus_faecalis | ceftazidime_avibactam | 0.1 | 0.005 |
| enterococcus_faecalis | meropenem_vaborbactam | 0.75 | 0.005 |
| enterococcus_faecalis | colistin | 0 | 0.005 |
| enterococcus_faecalis | flucloxacillin | 0.05 | 1 |
| enterococcus_faecalis | aztreonam_avibactam | 0.01 | 0.003 |
| enterococcus_faecalis | cefixime | 0.05 | 0.2 |
| enterococcus_faecium | sulfanilamide | 0.1 | 0.02 |
| enterococcus_faecium | penicillin_g | 0.1 | 1 |
| enterococcus_faecium | ampicillin | 0.3 | 1 |
| enterococcus_faecium | amoxicillin | 0.3 | 1 |
| enterococcus_faecium | piperacillin | 0.1 | 1 |
| enterococcus_faecium | ticarcillin | 0.1 | 1 |
| enterococcus_faecium | cephalexin | 0.1 | 0.3 |
| enterococcus_faecium | cefazolin | 0.1 | 0.3 |
| enterococcus_faecium | cefuroxime | 0.1 | 0.3 |
| enterococcus_faecium | ceftriaxone | 0.1 | 0.2 |
| enterococcus_faecium | ceftazidime | 0.1 | 0.2 |
| enterococcus_faecium | cefepime | 0.1 | 0.35 |
| enterococcus_faecium | ceftaroline | 0.1 | 0.002 |
| enterococcus_faecium | ceftolozane_tazobactam | 0.05 | 1 |
| enterococcus_faecium | cefiderocol | 0.05 | 1 |
| enterococcus_faecium | meropenem | 0.1 | 0.005 |
| enterococcus_faecium | imipenem_c | 0.1 | 0.005 |
| enterococcus_faecium | ertapenem | 0.1 | 0.005 |
| enterococcus_faecium | aztreonam | 0 | 0.003 |
| enterococcus_faecium | erythromycin | 0.7 | 0.01 |
| enterococcus_faecium | azithromycin | 0.7 | 0.01 |
| enterococcus_faecium | clarithromycin | 0.7 | 1 |
| enterococcus_faecium | clindamycin | 0.7 | 1 |
| enterococcus_faecium | gentamicin | 0.75 | 1 |
| enterococcus_faecium | tobramycin | 0.1 | 1 |
| enterococcus_faecium | amikacin | 0.1 | 1 |
| enterococcus_faecium | ciprofloxacin | 0.7 | 1 |
| enterococcus_faecium | levofloxacin | 0.7 | 1 |
| enterococcus_faecium | moxifloxacin | 0.7 | 1 |
| enterococcus_faecium | ofloxacin | 0.7 | 1 |
| enterococcus_faecium | tetracycline | 0.8 | 0.25 |
| enterococcus_faecium | doxycycline | 0.8 | 0.25 |
| enterococcus_faecium | minocycline | 0.85 | 0.25 |
| enterococcus_faecium | tigecycline | 0.1 | 1 |
| enterococcus_faecium | vancomycin | 0.9 | 5 |
| enterococcus_faecium | teicoplanin | 0.85 | 5 |
| enterococcus_faecium | dalbavancin | 0.85 | 0.005 |
| enterococcus_faecium | linezolid | 0.9 | 3 |
| enterococcus_faecium | tedizolid | 0.9 | 0.005 |
| enterococcus_faecium | daptomycin | 0.1 | 4 |
| enterococcus_faecium | quinu_dalfo | 0.7 | 0.005 |
| enterococcus_faecium | trim_sulf | 0.6 | 0.04 |
| enterococcus_faecium | chloramphenicol | 0.7 | 1 |
| enterococcus_faecium | nitrofurantoin | 0.7 | 15 |
| enterococcus_faecium | fosfomycin | 0.1 | 1 |
| enterococcus_faecium | retapamulin | 0.1 | 1 |
| enterococcus_faecium | fusidic_a | 0.1 | 1 |
| enterococcus_faecium | metronidazole | 0.1 | 1 |
| enterococcus_faecium | fidaxomicin | 0.1 | 1 |
| enterococcus_faecium | furazolidone | 0.1 | 1 |
| enterococcus_faecium | rifampicin | 0.1 | 1 |
| enterococcus_faecium | amoxicillin_clavulanate | 0.8 | 1 |
| enterococcus_faecium | piperacillin_tazobactam | 0.1 | 1 |
| enterococcus_faecium | ampicillin_sulbactam | 0.85 | 1 |
| enterococcus_faecium | ticarcillin_clavulanate | 0.1 | 1 |
| enterococcus_faecium | ceftazidime_avibactam | 0.1 | 0.005 |
| enterococcus_faecium | meropenem_vaborbactam | 0.1 | 0.005 |
| enterococcus_faecium | colistin | 0 | 0.005 |
| enterococcus_faecium | flucloxacillin | 0.05 | 1 |
| enterococcus_faecium | aztreonam_avibactam | 0.01 | 0.003 |
| enterococcus_faecium | cefixime | 0.05 | 0.2 |
| escherichia_coli | sulfanilamide | 0.5 | 0.02 |
| escherichia_coli | penicillin_g | 0.1 | 1 |
| escherichia_coli | ampicillin | 0.8 | 1 |
| escherichia_coli | amoxicillin | 0.8 | 6 |
| escherichia_coli | piperacillin | 0.85 | 1 |
| escherichia_coli | ticarcillin | 0.8 | 1 |
| escherichia_coli | cephalexin | 0.7 | 0.3 |
| escherichia_coli | cefazolin | 0.75 | 1.5 |
| escherichia_coli | cefuroxime | 0.8 | 0.3 |
| escherichia_coli | ceftriaxone | 0.9 | 1.5 |
| escherichia_coli | ceftazidime | 0.9 | 0.2 |
| escherichia_coli | cefepime | 0.9 | 2 |
| escherichia_coli | ceftaroline | 0.7 | 0.002 |
| escherichia_coli | ceftolozane_tazobactam | 0.8 | 1 |
| escherichia_coli | cefiderocol | 0.8 | 1 |
| escherichia_coli | meropenem | 0.95 | 8 |
| escherichia_coli | imipenem_c | 0.95 | 4 |
| escherichia_coli | ertapenem | 0.95 | 5 |
| escherichia_coli | aztreonam | 0.9 | 0.003 |
| escherichia_coli | gentamicin | 0.9 | 20 |
| escherichia_coli | tobramycin | 0.85 | 8 |
| escherichia_coli | amikacin | 0.9 | 12 |
| escherichia_coli | ciprofloxacin | 0.95 | 5 |
| escherichia_coli | levofloxacin | 0.9 | 1 |
| escherichia_coli | moxifloxacin | 0.8 | 1 |
| escherichia_coli | ofloxacin | 0.9 | 3 |
| escherichia_coli | tetracycline | 0.8 | 0.25 |
| escherichia_coli | doxycycline | 0.8 | 0.25 |
| escherichia_coli | minocycline | 0.85 | 0.25 |
| escherichia_coli | tigecycline | 0.1 | 1 |
| escherichia_coli | dalbavancin | 0 | 0.005 |
| escherichia_coli | linezolid | 0 | 0.005 |
| escherichia_coli | tedizolid | 0 | 0.005 |
| escherichia_coli | daptomycin | 0.1 | 1 |
| escherichia_coli | quinu_dalfo | 0 | 0.005 |
| escherichia_coli | trim_sulf | 0.9 | 0.06 |
| escherichia_coli | chloramphenicol | 0.85 | 1 |
| escherichia_coli | nitrofurantoin | 0.95 | 40 |
| escherichia_coli | fosfomycin | 0.1 | 20 |
| escherichia_coli | fidaxomicin | 0.1 | 1 |
| escherichia_coli | furazolidone | 0.1 | 1 |
| escherichia_coli | rifampicin | 0.7 | 1 |
| escherichia_coli | amoxicillin_clavulanate | 0.9 | 12 |
| escherichia_coli | piperacillin_tazobactam | 0.97 | 8 |
| escherichia_coli | ampicillin_sulbactam | 0.9 | 1 |
| escherichia_coli | ticarcillin_clavulanate | 0.9 | 3 |
| escherichia_coli | ceftazidime_avibactam | 0.95 | 0.005 |
| escherichia_coli | meropenem_vaborbactam | 0.95 | 0.005 |
| escherichia_coli | colistin | 0.7 | 0.005 |
| escherichia_coli | flucloxacillin | 0.01 | 1 |
| escherichia_coli | aztreonam_avibactam | 1 | 0.003 |
| escherichia_coli | cefixime | 0.8 | 0.2 |
| klebsiella_pneumoniae | sulfanilamide | 0.5 | 0.02 |
| klebsiella_pneumoniae | penicillin_g | 0.1 | 1 |
| klebsiella_pneumoniae | ampicillin | 0.1 | 1 |
| klebsiella_pneumoniae | amoxicillin | 0.1 | 1 |
| klebsiella_pneumoniae | piperacillin | 0.8 | 1 |
| klebsiella_pneumoniae | ticarcillin | 0.75 | 1 |
| klebsiella_pneumoniae | cephalexin | 0.5 | 0.3 |
| klebsiella_pneumoniae | cefazolin | 0.5 | 1.5 |
| klebsiella_pneumoniae | cefuroxime | 0.7 | 0.3 |
| klebsiella_pneumoniae | ceftriaxone | 0.9 | 1.5 |
| klebsiella_pneumoniae | ceftazidime | 0.85 | 0.2 |
| klebsiella_pneumoniae | cefepime | 0.92 | 2 |
| klebsiella_pneumoniae | ceftaroline | 0.5 | 0.002 |
| klebsiella_pneumoniae | ceftolozane_tazobactam | 0.8 | 1 |
| klebsiella_pneumoniae | cefiderocol | 0.8 | 1 |
| klebsiella_pneumoniae | meropenem | 0.94 | 8 |
| klebsiella_pneumoniae | imipenem_c | 0.95 | 4 |
| klebsiella_pneumoniae | ertapenem | 0.94 | 5 |
| klebsiella_pneumoniae | aztreonam | 0.85 | 0.005 |
| klebsiella_pneumoniae | gentamicin | 0.9 | 15 |
| klebsiella_pneumoniae | tobramycin | 0.85 | 8 |
| klebsiella_pneumoniae | amikacin | 0.9 | 12 |
| klebsiella_pneumoniae | ciprofloxacin | 0.9 | 4 |
| klebsiella_pneumoniae | levofloxacin | 0.85 | 1 |
| klebsiella_pneumoniae | moxifloxacin | 0.7 | 1 |
| klebsiella_pneumoniae | ofloxacin | 0.8 | 1 |
| klebsiella_pneumoniae | tetracycline | 0.8 | 0.25 |
| klebsiella_pneumoniae | doxycycline | 0.8 | 0.25 |
| klebsiella_pneumoniae | minocycline | 0.85 | 0.25 |
| klebsiella_pneumoniae | tigecycline | 0.1 | 1 |
| klebsiella_pneumoniae | dalbavancin | 0 | 0.005 |
| klebsiella_pneumoniae | linezolid | 0 | 0.005 |
| klebsiella_pneumoniae | tedizolid | 0 | 0.005 |
| klebsiella_pneumoniae | daptomycin | 0.1 | 1 |
| klebsiella_pneumoniae | quinu_dalfo | 0 | 0.005 |
| klebsiella_pneumoniae | trim_sulf | 0.9 | 0.04 |
| klebsiella_pneumoniae | chloramphenicol | 0.85 | 1 |
| klebsiella_pneumoniae | nitrofurantoin | 0.8 | 25 |
| klebsiella_pneumoniae | fosfomycin | 0.1 | 10 |
| klebsiella_pneumoniae | fidaxomicin | 0.1 | 1 |
| klebsiella_pneumoniae | furazolidone | 0.1 | 1 |
| klebsiella_pneumoniae | rifampicin | 0.6 | 1 |
| klebsiella_pneumoniae | amoxicillin_clavulanate | 0.85 | 12 |
| klebsiella_pneumoniae | piperacillin_tazobactam | 0.92 | 8 |
| klebsiella_pneumoniae | ampicillin_sulbactam | 0.75 | 1 |
| klebsiella_pneumoniae | ticarcillin_clavulanate | 0.75 | 3 |
| klebsiella_pneumoniae | ceftazidime_avibactam | 0.95 | 0.005 |
| klebsiella_pneumoniae | meropenem_vaborbactam | 0.95 | 0.005 |
| klebsiella_pneumoniae | colistin | 0.7 | 0.005 |
| klebsiella_pneumoniae | flucloxacillin | 0.01 | 1 |
| klebsiella_pneumoniae | aztreonam_avibactam | 1 | 0.005 |
| klebsiella_pneumoniae | cefixime | 0.8 | 0.2 |
| morganella_spp. | sulfanilamide | 0.5 | 0.02 |
| morganella_spp. | penicillin_g | 0.1 | 1 |
| morganella_spp. | ampicillin | 0.5 | 1 |
| morganella_spp. | amoxicillin | 0.5 | 1 |
| morganella_spp. | piperacillin | 0.75 | 1 |
| morganella_spp. | ticarcillin | 0.7 | 1 |
| morganella_spp. | cephalexin | 0.5 | 0.3 |
| morganella_spp. | cefazolin | 0.5 | 0.3 |
| morganella_spp. | cefuroxime | 0.6 | 0.3 |
| morganella_spp. | ceftriaxone | 0.8 | 0.2 |
| morganella_spp. | ceftazidime | 0.8 | 0.2 |
| morganella_spp. | cefepime | 0.85 | 0.35 |
| morganella_spp. | ceftaroline | 0.4 | 0.002 |
| morganella_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| morganella_spp. | cefiderocol | 0.8 | 1 |
| morganella_spp. | meropenem | 0.95 | 4 |
| morganella_spp. | imipenem_c | 0.95 | 0.005 |
| morganella_spp. | ertapenem | 0.9 | 0.005 |
| morganella_spp. | aztreonam | 0.8 | 0.003 |
| morganella_spp. | erythromycin | 0.1 | 1 |
| morganella_spp. | azithromycin | 0.1 | 1 |
| morganella_spp. | clarithromycin | 0.1 | 1 |
| morganella_spp. | clindamycin | 0.1 | 1 |
| morganella_spp. | gentamicin | 0.85 | 10 |
| morganella_spp. | tobramycin | 0.8 | 1 |
| morganella_spp. | amikacin | 0.9 | 1 |
| morganella_spp. | ciprofloxacin | 0.9 | 1 |
| morganella_spp. | levofloxacin | 0.85 | 1 |
| morganella_spp. | moxifloxacin | 0.7 | 1 |
| morganella_spp. | ofloxacin | 0.8 | 1 |
| morganella_spp. | tetracycline | 0.8 | 0.25 |
| morganella_spp. | doxycycline | 0.85 | 0.25 |
| morganella_spp. | minocycline | 0.85 | 0.25 |
| morganella_spp. | tigecycline | 0.1 | 1 |
| morganella_spp. | vancomycin | 0.1 | 1 |
| morganella_spp. | teicoplanin | 0.1 | 1 |
| morganella_spp. | dalbavancin | 0.1 | 0.005 |
| morganella_spp. | linezolid | 0.1 | 0.005 |
| morganella_spp. | tedizolid | 0.1 | 0.005 |
| morganella_spp. | daptomycin | 0.1 | 1 |
| morganella_spp. | quinu_dalfo | 0.1 | 0.005 |
| morganella_spp. | trim_sulf | 0.8 | 0.04 |
| morganella_spp. | chloramphenicol | 0.85 | 1 |
| morganella_spp. | nitrofurantoin | 0.7 | 1 |
| morganella_spp. | fosfomycin | 0.1 | 1 |
| morganella_spp. | retapamulin | 0.05 | 1 |
| morganella_spp. | fusidic_a | 0.05 | 1 |
| morganella_spp. | metronidazole | 0.05 | 1 |
| morganella_spp. | fidaxomicin | 0.1 | 1 |
| morganella_spp. | furazolidone | 0.1 | 1 |
| morganella_spp. | rifampicin | 0.6 | 1 |
| morganella_spp. | amoxicillin_clavulanate | 0.7 | 6 |
| morganella_spp. | piperacillin_tazobactam | 0.85 | 8 |
| morganella_spp. | ampicillin_sulbactam | 0.7 | 1 |
| morganella_spp. | ticarcillin_clavulanate | 0.8 | 1 |
| morganella_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| morganella_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| morganella_spp. | colistin | 0.7 | 0.005 |
| morganella_spp. | flucloxacillin | 0.01 | 1 |
| morganella_spp. | aztreonam_avibactam | 1 | 0.003 |
| morganella_spp. | cefixime | 0.8 | 0.2 |
| proteus_spp. | sulfanilamide | 0.5 | 0.02 |
| proteus_spp. | penicillin_g | 0.1 | 1 |
| proteus_spp. | ampicillin | 0.8 | 1 |
| proteus_spp. | amoxicillin | 0.8 | 1 |
| proteus_spp. | piperacillin | 0.85 | 1 |
| proteus_spp. | ticarcillin | 0.8 | 1 |
| proteus_spp. | cephalexin | 0.7 | 0.3 |
| proteus_spp. | cefazolin | 0.75 | 0.3 |
| proteus_spp. | cefuroxime | 0.8 | 0.3 |
| proteus_spp. | ceftriaxone | 0.95 | 0.2 |
| proteus_spp. | ceftazidime | 0.9 | 0.2 |
| proteus_spp. | cefepime | 0.9 | 0.35 |
| proteus_spp. | ceftaroline | 0.7 | 0.002 |
| proteus_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| proteus_spp. | cefiderocol | 0.8 | 1 |
| proteus_spp. | meropenem | 0.95 | 4 |
| proteus_spp. | imipenem_c | 0.95 | 0.005 |
| proteus_spp. | ertapenem | 0.95 | 3 |
| proteus_spp. | aztreonam | 0.9 | 0.003 |
| proteus_spp. | gentamicin | 0.8 | 10 |
| proteus_spp. | tobramycin | 0.75 | 1 |
| proteus_spp. | amikacin | 0.85 | 1 |
| proteus_spp. | ciprofloxacin | 0.9 | 4 |
| proteus_spp. | levofloxacin | 0.85 | 1 |
| proteus_spp. | moxifloxacin | 0.7 | 1 |
| proteus_spp. | ofloxacin | 0.8 | 1 |
| proteus_spp. | tetracycline | 0.8 | 0.25 |
| proteus_spp. | doxycycline | 0.8 | 0.25 |
| proteus_spp. | minocycline | 0.85 | 0.25 |
| proteus_spp. | tigecycline | 0.1 | 1 |
| proteus_spp. | dalbavancin | 0 | 0.005 |
| proteus_spp. | linezolid | 0 | 0.005 |
| proteus_spp. | tedizolid | 0 | 0.005 |
| proteus_spp. | daptomycin | 0.1 | 1 |
| proteus_spp. | quinu_dalfo | 0 | 0.005 |
| proteus_spp. | trim_sulf | 0.9 | 0.04 |
| proteus_spp. | chloramphenicol | 0.85 | 1 |
| proteus_spp. | nitrofurantoin | 0.8 | 0.05 |
| proteus_spp. | fosfomycin | 0.1 | 1 |
| proteus_spp. | fidaxomicin | 0.1 | 1 |
| proteus_spp. | furazolidone | 0.1 | 1 |
| proteus_spp. | rifampicin | 0.7 | 1 |
| proteus_spp. | amoxicillin_clavulanate | 0.9 | 12 |
| proteus_spp. | piperacillin_tazobactam | 0.95 | 8 |
| proteus_spp. | ampicillin_sulbactam | 0.9 | 1 |
| proteus_spp. | ticarcillin_clavulanate | 0.9 | 3 |
| proteus_spp. | ceftazidime_avibactam | 0.95 | 0.005 |
| proteus_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| proteus_spp. | colistin | 0.7 | 0.005 |
| proteus_spp. | flucloxacillin | 0.01 | 1 |
| proteus_spp. | aztreonam_avibactam | 1 | 0.003 |
| proteus_spp. | cefixime | 0.8 | 0.2 |
| serratia_spp. | sulfanilamide | 0.5 | 0.02 |
| serratia_spp. | penicillin_g | 0.1 | 1 |
| serratia_spp. | ampicillin | 0.1 | 1 |
| serratia_spp. | amoxicillin | 0.1 | 1 |
| serratia_spp. | piperacillin | 0.75 | 1 |
| serratia_spp. | ticarcillin | 0.7 | 1 |
| serratia_spp. | cephalexin | 0.1 | 0.3 |
| serratia_spp. | cefazolin | 0.1 | 0.3 |
| serratia_spp. | cefuroxime | 0.6 | 0.3 |
| serratia_spp. | ceftriaxone | 0.8 | 0.2 |
| serratia_spp. | ceftazidime | 0.85 | 0.2 |
| serratia_spp. | cefepime | 0.85 | 0.35 |
| serratia_spp. | ceftaroline | 0.5 | 0.002 |
| serratia_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| serratia_spp. | cefiderocol | 0.8 | 1 |
| serratia_spp. | meropenem | 0.95 | 4 |
| serratia_spp. | imipenem_c | 0.95 | 0.005 |
| serratia_spp. | ertapenem | 0.9 | 3 |
| serratia_spp. | aztreonam | 0.85 | 0.003 |
| serratia_spp. | erythromycin | 0.1 | 1 |
| serratia_spp. | azithromycin | 0.1 | 1 |
| serratia_spp. | clarithromycin | 0.1 | 1 |
| serratia_spp. | clindamycin | 0.1 | 1 |
| serratia_spp. | gentamicin | 0.85 | 10 |
| serratia_spp. | tobramycin | 0.8 | 1 |
| serratia_spp. | amikacin | 0.9 | 8 |
| serratia_spp. | ciprofloxacin | 0.85 | 1 |
| serratia_spp. | levofloxacin | 0.8 | 1 |
| serratia_spp. | moxifloxacin | 0.7 | 1 |
| serratia_spp. | ofloxacin | 0.75 | 1 |
| serratia_spp. | tetracycline | 0.8 | 0.25 |
| serratia_spp. | doxycycline | 0.8 | 0.25 |
| serratia_spp. | minocycline | 0.85 | 0.25 |
| serratia_spp. | tigecycline | 0.1 | 1 |
| serratia_spp. | vancomycin | 0.1 | 1 |
| serratia_spp. | teicoplanin | 0.1 | 1 |
| serratia_spp. | dalbavancin | 0.1 | 0.005 |
| serratia_spp. | linezolid | 0.1 | 0.005 |
| serratia_spp. | tedizolid | 0.1 | 0.005 |
| serratia_spp. | daptomycin | 0.1 | 1 |
| serratia_spp. | quinu_dalfo | 0.1 | 0.005 |
| serratia_spp. | trim_sulf | 0.85 | 0.04 |
| serratia_spp. | chloramphenicol | 0.8 | 1 |
| serratia_spp. | nitrofurantoin | 0.7 | 1 |
| serratia_spp. | fosfomycin | 0.1 | 1 |
| serratia_spp. | retapamulin | 0.05 | 1 |
| serratia_spp. | fusidic_a | 0.05 | 1 |
| serratia_spp. | metronidazole | 0.05 | 1 |
| serratia_spp. | fidaxomicin | 0.1 | 1 |
| serratia_spp. | furazolidone | 0.1 | 1 |
| serratia_spp. | rifampicin | 0.6 | 1 |
| serratia_spp. | amoxicillin_clavulanate | 0.7 | 6 |
| serratia_spp. | piperacillin_tazobactam | 0.85 | 8 |
| serratia_spp. | ampicillin_sulbactam | 0.7 | 1 |
| serratia_spp. | ticarcillin_clavulanate | 0.75 | 1 |
| serratia_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| serratia_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| serratia_spp. | colistin | 0.7 | 0.005 |
| serratia_spp. | flucloxacillin | 0.01 | 1 |
| serratia_spp. | aztreonam_avibactam | 1 | 0.003 |
| serratia_spp. | cefixime | 0.8 | 0.2 |
| p_stuartii | sulfanilamide | 0.1 | 0.02 |
| p_stuartii | penicillin_g | 0.05 | 1 |
| p_stuartii | ampicillin | 0.05 | 1 |
| p_stuartii | amoxicillin | 0.05 | 1 |
| p_stuartii | piperacillin | 0.35 | 1 |
| p_stuartii | ticarcillin | 0.3 | 1 |
| p_stuartii | cephalexin | 0.1 | 0.3 |
| p_stuartii | cefazolin | 0.1 | 0.3 |
| p_stuartii | cefuroxime | 0.2 | 0.3 |
| p_stuartii | ceftriaxone | 0.45 | 0.2 |
| p_stuartii | ceftazidime | 0.75 | 0.2 |
| p_stuartii | cefepime | 0.85 | 0.35 |
| p_stuartii | ceftaroline | 0.2 | 0.002 |
| p_stuartii | ceftolozane_tazobactam | 0.8 | 1 |
| p_stuartii | cefiderocol | 0.8 | 1 |
| p_stuartii | meropenem | 0.9 | 0.005 |
| p_stuartii | imipenem_c | 0.9 | 0.005 |
| p_stuartii | ertapenem | 0.85 | 0.005 |
| p_stuartii | aztreonam | 0.65 | 0.003 |
| p_stuartii | erythromycin | 0.05 | 1 |
| p_stuartii | azithromycin | 0.05 | 1 |
| p_stuartii | clarithromycin | 0.05 | 1 |
| p_stuartii | clindamycin | 0.05 | 1 |
| p_stuartii | gentamicin | 0.45 | 1 |
| p_stuartii | tobramycin | 0.5 | 1 |
| p_stuartii | amikacin | 0.75 | 1 |
| p_stuartii | ciprofloxacin | 0.55 | 1 |
| p_stuartii | levofloxacin | 0.6 | 1 |
| p_stuartii | moxifloxacin | 0.6 | 1 |
| p_stuartii | ofloxacin | 0.55 | 1 |
| p_stuartii | tetracycline | 0.2 | 0.25 |
| p_stuartii | doxycycline | 0.3 | 0.25 |
| p_stuartii | minocycline | 0.35 | 0.25 |
| p_stuartii | tigecycline | 0.1 | 1 |
| p_stuartii | vancomycin | 0.05 | 1 |
| p_stuartii | teicoplanin | 0.05 | 1 |
| p_stuartii | dalbavancin | 0.05 | 0.005 |
| p_stuartii | linezolid | 0.05 | 0.005 |
| p_stuartii | tedizolid | 0.05 | 0.005 |
| p_stuartii | daptomycin | 0.1 | 1 |
| p_stuartii | quinu_dalfo | 0.05 | 0.005 |
| p_stuartii | trim_sulf | 0.3 | 0.04 |
| p_stuartii | chloramphenicol | 0.25 | 1 |
| p_stuartii | nitrofurantoin | 0.05 | 1 |
| p_stuartii | fosfomycin | 0.1 | 1 |
| p_stuartii | retapamulin | 0.05 | 1 |
| p_stuartii | fusidic_a | 0.05 | 1 |
| p_stuartii | metronidazole | 0.05 | 1 |
| p_stuartii | fidaxomicin | 0.1 | 1 |
| p_stuartii | furazolidone | 0.05 | 1 |
| p_stuartii | rifampicin | 0.1 | 1 |
| p_stuartii | amoxicillin_clavulanate | 0.2 | 1 |
| p_stuartii | piperacillin_tazobactam | 0.75 | 1 |
| p_stuartii | ampicillin_sulbactam | 0.2 | 1 |
| p_stuartii | ticarcillin_clavulanate | 0.45 | 1 |
| p_stuartii | ceftazidime_avibactam | 0.9 | 0.005 |
| p_stuartii | meropenem_vaborbactam | 0.95 | 0.005 |
| p_stuartii | colistin | 0.05 | 0.005 |
| p_stuartii | flucloxacillin | 0.01 | 1 |
| p_stuartii | aztreonam_avibactam | 1 | 0.003 |
| p_stuartii | cefixime | 0.8 | 0.2 |
| pseudomonas_aeruginosa | sulfanilamide | 0.1 | 0.02 |
| pseudomonas_aeruginosa | penicillin_g | 0.05 | 0.01 |
| pseudomonas_aeruginosa | ampicillin | 0.05 | 0.01 |
| pseudomonas_aeruginosa | amoxicillin | 0.05 | 0.01 |
| pseudomonas_aeruginosa | piperacillin | 0.8 | 1 |
| pseudomonas_aeruginosa | ticarcillin | 0.7 | 1 |
| pseudomonas_aeruginosa | cephalexin | 0.05 | 0.3 |
| pseudomonas_aeruginosa | cefazolin | 0.05 | 0.3 |
| pseudomonas_aeruginosa | cefuroxime | 0.1 | 0.3 |
| pseudomonas_aeruginosa | ceftriaxone | 0.1 | 0.2 |
| pseudomonas_aeruginosa | ceftazidime | 0.85 | 3 |
| pseudomonas_aeruginosa | cefepime | 0.9 | 3 |
| pseudomonas_aeruginosa | ceftaroline | 0.1 | 0.002 |
| pseudomonas_aeruginosa | ceftolozane_tazobactam | 0.1 | 1 |
| pseudomonas_aeruginosa | cefiderocol | 0.1 | 1 |
| pseudomonas_aeruginosa | meropenem | 0.9 | 5 |
| pseudomonas_aeruginosa | imipenem_c | 0.85 | 3 |
| pseudomonas_aeruginosa | ertapenem | 0.1 | 0.005 |
| pseudomonas_aeruginosa | aztreonam | 0.8 | 0.05 |
| pseudomonas_aeruginosa | gentamicin | 0.85 | 12 |
| pseudomonas_aeruginosa | tobramycin | 0.9 | 15 |
| pseudomonas_aeruginosa | amikacin | 0.9 | 15 |
| pseudomonas_aeruginosa | ciprofloxacin | 0.9 | 5 |
| pseudomonas_aeruginosa | levofloxacin | 0.8 | 1 |
| pseudomonas_aeruginosa | moxifloxacin | 0.5 | 1 |
| pseudomonas_aeruginosa | ofloxacin | 0.7 | 1 |
| pseudomonas_aeruginosa | tetracycline | 0.1 | 0.25 |
| pseudomonas_aeruginosa | doxycycline | 0.1 | 0.25 |
| pseudomonas_aeruginosa | minocycline | 0.1 | 0.25 |
| pseudomonas_aeruginosa | tigecycline | 0.1 | 1 |
| pseudomonas_aeruginosa | dalbavancin | 0 | 0.005 |
| pseudomonas_aeruginosa | linezolid | 0 | 0.005 |
| pseudomonas_aeruginosa | tedizolid | 0 | 0.005 |
| pseudomonas_aeruginosa | daptomycin | 0.1 | 1 |
| pseudomonas_aeruginosa | quinu_dalfo | 0 | 0.005 |
| pseudomonas_aeruginosa | trim_sulf | 0.1 | 0.04 |
| pseudomonas_aeruginosa | chloramphenicol | 0.1 | 1 |
| pseudomonas_aeruginosa | nitrofurantoin | 0.05 | 0.01 |
| pseudomonas_aeruginosa | fosfomycin | 0.1 | 0.05 |
| pseudomonas_aeruginosa | fidaxomicin | 0.1 | 1 |
| pseudomonas_aeruginosa | furazolidone | 0.05 | 1 |
| pseudomonas_aeruginosa | rifampicin | 0.1 | 1 |
| pseudomonas_aeruginosa | amoxicillin_clavulanate | 0.05 | 1 |
| pseudomonas_aeruginosa | piperacillin_tazobactam | 0.9 | 8 |
| pseudomonas_aeruginosa | ampicillin_sulbactam | 0.05 | 1 |
| pseudomonas_aeruginosa | ticarcillin_clavulanate | 0.8 | 3 |
| pseudomonas_aeruginosa | ceftazidime_avibactam | 0.95 | 0.005 |
| pseudomonas_aeruginosa | meropenem_vaborbactam | 0.9 | 0.005 |
| pseudomonas_aeruginosa | colistin | 0.85 | 0.05 |
| pseudomonas_aeruginosa | flucloxacillin | 0.01 | 1 |
| pseudomonas_aeruginosa | aztreonam_avibactam | 0.9 | 0.04 |
| pseudomonas_aeruginosa | cefixime | 0.1 | 0.2 |
| stenotrophomonas_maltophilia | sulfanilamide | 0.6 | 0.02 |
| stenotrophomonas_maltophilia | penicillin_g | 0.05 | 1 |
| stenotrophomonas_maltophilia | ampicillin | 0.05 | 1 |
| stenotrophomonas_maltophilia | amoxicillin | 0.05 | 1 |
| stenotrophomonas_maltophilia | piperacillin | 0.2 | 1 |
| stenotrophomonas_maltophilia | ticarcillin | 0.25 | 1 |
| stenotrophomonas_maltophilia | cephalexin | 0.05 | 0.3 |
| stenotrophomonas_maltophilia | cefazolin | 0.05 | 0.3 |
| stenotrophomonas_maltophilia | cefuroxime | 0.05 | 0.3 |
| stenotrophomonas_maltophilia | ceftriaxone | 0.05 | 0.2 |
| stenotrophomonas_maltophilia | ceftazidime | 0.35 | 0.2 |
| stenotrophomonas_maltophilia | cefepime | 0.15 | 0.35 |
| stenotrophomonas_maltophilia | ceftaroline | 0.05 | 0.002 |
| stenotrophomonas_maltophilia | ceftolozane_tazobactam | 0.1 | 1 |
| stenotrophomonas_maltophilia | cefiderocol | 0.1 | 1 |
| stenotrophomonas_maltophilia | meropenem | 0.05 | 0.01 |
| stenotrophomonas_maltophilia | imipenem_c | 0.05 | 0.01 |
| stenotrophomonas_maltophilia | ertapenem | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | aztreonam | 0.1 | 0.003 |
| stenotrophomonas_maltophilia | erythromycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | azithromycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | clarithromycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | clindamycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | gentamicin | 0.1 | 1 |
| stenotrophomonas_maltophilia | tobramycin | 0.1 | 1 |
| stenotrophomonas_maltophilia | amikacin | 0.1 | 1 |
| stenotrophomonas_maltophilia | ciprofloxacin | 0.4 | 1 |
| stenotrophomonas_maltophilia | levofloxacin | 0.75 | 3.5 |
| stenotrophomonas_maltophilia | moxifloxacin | 0.8 | 1 |
| stenotrophomonas_maltophilia | ofloxacin | 0.6 | 1 |
| stenotrophomonas_maltophilia | tetracycline | 0.35 | 0.25 |
| stenotrophomonas_maltophilia | doxycycline | 0.6 | 4.5 |
| stenotrophomonas_maltophilia | minocycline | 0.85 | 6 |
| stenotrophomonas_maltophilia | tigecycline | 0.1 | 1 |
| stenotrophomonas_maltophilia | vancomycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | teicoplanin | 0.05 | 1 |
| stenotrophomonas_maltophilia | dalbavancin | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | linezolid | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | tedizolid | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | daptomycin | 0.1 | 1 |
| stenotrophomonas_maltophilia | quinu_dalfo | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | trim_sulf | 0.95 | 5 |
| stenotrophomonas_maltophilia | chloramphenicol | 0.4 | 1 |
| stenotrophomonas_maltophilia | nitrofurantoin | 0.05 | 1 |
| stenotrophomonas_maltophilia | fosfomycin | 0.1 | 1 |
| stenotrophomonas_maltophilia | retapamulin | 0.05 | 1 |
| stenotrophomonas_maltophilia | fusidic_a | 0.05 | 1 |
| stenotrophomonas_maltophilia | metronidazole | 0.05 | 1 |
| stenotrophomonas_maltophilia | fidaxomicin | 0.1 | 1 |
| stenotrophomonas_maltophilia | furazolidone | 0.05 | 1 |
| stenotrophomonas_maltophilia | rifampicin | 0.2 | 1 |
| stenotrophomonas_maltophilia | amoxicillin_clavulanate | 0.05 | 1 |
| stenotrophomonas_maltophilia | piperacillin_tazobactam | 0.3 | 0.05 |
| stenotrophomonas_maltophilia | ampicillin_sulbactam | 0.05 | 1 |
| stenotrophomonas_maltophilia | ticarcillin_clavulanate | 0.7 | 1 |
| stenotrophomonas_maltophilia | ceftazidime_avibactam | 0.4 | 0.005 |
| stenotrophomonas_maltophilia | meropenem_vaborbactam | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | colistin | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | flucloxacillin | 0.01 | 1 |
| stenotrophomonas_maltophilia | aztreonam_avibactam | 0.75 | 0.003 |
| stenotrophomonas_maltophilia | cefixime | 0.1 | 0.2 |
| staphylococcus_aureus | sulfanilamide | 0.1 | 0.02 |
| staphylococcus_aureus | penicillin_g | 0.95 | 1 |
| staphylococcus_aureus | ampicillin | 0.1 | 1 |
| staphylococcus_aureus | amoxicillin | 0.1 | 6 |
| staphylococcus_aureus | piperacillin | 0.7 | 1 |
| staphylococcus_aureus | ticarcillin | 0.6 | 1 |
| staphylococcus_aureus | cephalexin | 0.8 | 5 |
| staphylococcus_aureus | cefazolin | 0.85 | 5 |
| staphylococcus_aureus | cefuroxime | 0.7 | 3 |
| staphylococcus_aureus | ceftriaxone | 0.7 | 0.2 |
| staphylococcus_aureus | ceftazidime | 0.1 | 0.2 |
| staphylococcus_aureus | cefepime | 0.6 | 0.35 |
| staphylococcus_aureus | ceftaroline | 0.95 | 0.15 |
| staphylococcus_aureus | ceftolozane_tazobactam | 0.75 | 1 |
| staphylococcus_aureus | cefiderocol | 0.75 | 1 |
| staphylococcus_aureus | meropenem | 0.7 | 0.005 |
| staphylococcus_aureus | imipenem_c | 0.7 | 0.005 |
| staphylococcus_aureus | ertapenem | 0.7 | 0.005 |
| staphylococcus_aureus | aztreonam | 0 | 0.003 |
| staphylococcus_aureus | erythromycin | 0.8 | 1 |
| staphylococcus_aureus | azithromycin | 0.8 | 1 |
| staphylococcus_aureus | clarithromycin | 0.8 | 1 |
| staphylococcus_aureus | clindamycin | 0.8 | 1 |
| staphylococcus_aureus | gentamicin | 0.7 | 15 |
| staphylococcus_aureus | tobramycin | 0.7 | 1 |
| staphylococcus_aureus | amikacin | 0.7 | 1 |
| staphylococcus_aureus | ciprofloxacin | 0.7 | 1 |
| staphylococcus_aureus | levofloxacin | 0.7 | 1 |
| staphylococcus_aureus | moxifloxacin | 0.8 | 1 |
| staphylococcus_aureus | ofloxacin | 0.7 | 1 |
| staphylococcus_aureus | tetracycline | 0.8 | 0.25 |
| staphylococcus_aureus | doxycycline | 0.85 | 0.25 |
| staphylococcus_aureus | minocycline | 0.85 | 0.25 |
| staphylococcus_aureus | tigecycline | 0.1 | 1 |
| staphylococcus_aureus | vancomycin | 0.95 | 5 |
| staphylococcus_aureus | teicoplanin | 0.9 | 4 |
| staphylococcus_aureus | dalbavancin | 0.9 | 0.005 |
| staphylococcus_aureus | linezolid | 0.9 | 5 |
| staphylococcus_aureus | tedizolid | 0.9 | 0.005 |
| staphylococcus_aureus | daptomycin | 0.1 | 4 |
| staphylococcus_aureus | quinu_dalfo | 0.85 | 0.005 |
| staphylococcus_aureus | trim_sulf | 0.7 | 0.04 |
| staphylococcus_aureus | chloramphenicol | 0.8 | 1 |
| staphylococcus_aureus | nitrofurantoin | 0.1 | 8 |
| staphylococcus_aureus | fosfomycin | 0.1 | 6 |
| staphylococcus_aureus | retapamulin | 0.9 | 1 |
| staphylococcus_aureus | fusidic_a | 0.85 | 1 |
| staphylococcus_aureus | metronidazole | 0.1 | 1 |
| staphylococcus_aureus | fidaxomicin | 0.1 | 1 |
| staphylococcus_aureus | furazolidone | 0.1 | 1 |
| staphylococcus_aureus | rifampicin | 0.8 | 1 |
| staphylococcus_aureus | amoxicillin_clavulanate | 0.85 | 12 |
| staphylococcus_aureus | piperacillin_tazobactam | 0.7 | 1 |
| staphylococcus_aureus | ampicillin_sulbactam | 0.8 | 1 |
| staphylococcus_aureus | ticarcillin_clavulanate | 0.6 | 1 |
| staphylococcus_aureus | ceftazidime_avibactam | 0.1 | 0.005 |
| staphylococcus_aureus | meropenem_vaborbactam | 0.7 | 0.005 |
| staphylococcus_aureus | colistin | 0 | 0.005 |
| staphylococcus_aureus | flucloxacillin | 0.95 | 4 |
| staphylococcus_aureus | aztreonam_avibactam | 0.01 | 0.003 |
| staphylococcus_aureus | cefixime | 0.75 | 0.2 |
| staphylococcus_epidermidis | sulfanilamide | 0.1 | 0.02 |
| staphylococcus_epidermidis | penicillin_g | 0.15 | 1 |
| staphylococcus_epidermidis | ampicillin | 0.15 | 1 |
| staphylococcus_epidermidis | amoxicillin | 0.15 | 6 |
| staphylococcus_epidermidis | piperacillin | 0.2 | 1 |
| staphylococcus_epidermidis | ticarcillin | 0.2 | 1 |
| staphylococcus_epidermidis | cephalexin | 0.2 | 0.3 |
| staphylococcus_epidermidis | cefazolin | 0.2 | 0.3 |
| staphylococcus_epidermidis | cefuroxime | 0.2 | 0.3 |
| staphylococcus_epidermidis | ceftriaxone | 0.25 | 0.2 |
| staphylococcus_epidermidis | ceftazidime | 0.1 | 0.2 |
| staphylococcus_epidermidis | cefepime | 0.15 | 0.35 |
| staphylococcus_epidermidis | ceftaroline | 0.75 | 0.04 |
| staphylococcus_epidermidis | ceftolozane_tazobactam | 0.75 | 1 |
| staphylococcus_epidermidis | cefiderocol | 0.75 | 1 |
| staphylococcus_epidermidis | meropenem | 0.4 | 0.005 |
| staphylococcus_epidermidis | imipenem_c | 0.5 | 0.005 |
| staphylococcus_epidermidis | ertapenem | 0.4 | 0.005 |
| staphylococcus_epidermidis | aztreonam | 0.05 | 0.003 |
| staphylococcus_epidermidis | erythromycin | 0.45 | 1 |
| staphylococcus_epidermidis | azithromycin | 0.5 | 1 |
| staphylococcus_epidermidis | clarithromycin | 0.5 | 1 |
| staphylococcus_epidermidis | clindamycin | 0.6 | 1 |
| staphylococcus_epidermidis | gentamicin | 0.6 | 1 |
| staphylococcus_epidermidis | tobramycin | 0.65 | 1 |
| staphylococcus_epidermidis | amikacin | 0.7 | 1 |
| staphylococcus_epidermidis | ciprofloxacin | 0.5 | 1 |
| staphylococcus_epidermidis | levofloxacin | 0.55 | 1 |
| staphylococcus_epidermidis | moxifloxacin | 0.6 | 1 |
| staphylococcus_epidermidis | ofloxacin | 0.5 | 1 |
| staphylococcus_epidermidis | tetracycline | 0.5 | 0.25 |
| staphylococcus_epidermidis | doxycycline | 0.75 | 0.25 |
| staphylococcus_epidermidis | minocycline | 0.8 | 0.25 |
| staphylococcus_epidermidis | tigecycline | 0.1 | 1 |
| staphylococcus_epidermidis | vancomycin | 0.95 | 4 |
| staphylococcus_epidermidis | teicoplanin | 0.95 | 4 |
| staphylococcus_epidermidis | dalbavancin | 0.95 | 0.005 |
| staphylococcus_epidermidis | linezolid | 0.95 | 5 |
| staphylococcus_epidermidis | tedizolid | 0.95 | 0.005 |
| staphylococcus_epidermidis | daptomycin | 0.1 | 3 |
| staphylococcus_epidermidis | quinu_dalfo | 0.9 | 4 |
| staphylococcus_epidermidis | trim_sulf | 0.75 | 1.2 |
| staphylococcus_epidermidis | chloramphenicol | 0.6 | 1 |
| staphylococcus_epidermidis | nitrofurantoin | 0.2 | 8 |
| staphylococcus_epidermidis | fosfomycin | 0.1 | 6 |
| staphylococcus_epidermidis | retapamulin | 0.8 | 1 |
| staphylococcus_epidermidis | fusidic_a | 0.85 | 1 |
| staphylococcus_epidermidis | metronidazole | 0.05 | 1 |
| staphylococcus_epidermidis | fidaxomicin | 0.1 | 1 |
| staphylococcus_epidermidis | furazolidone | 0.05 | 1 |
| staphylococcus_epidermidis | rifampicin | 0.9 | 1 |
| staphylococcus_epidermidis | amoxicillin_clavulanate | 0.2 | 12 |
| staphylococcus_epidermidis | piperacillin_tazobactam | 0.4 | 1 |
| staphylococcus_epidermidis | ampicillin_sulbactam | 0.2 | 1 |
| staphylococcus_epidermidis | ticarcillin_clavulanate | 0.25 | 1 |
| staphylococcus_epidermidis | ceftazidime_avibactam | 0.1 | 0.005 |
| staphylococcus_epidermidis | meropenem_vaborbactam | 0.4 | 0.005 |
| staphylococcus_epidermidis | colistin | 0.05 | 0.005 |
| staphylococcus_epidermidis | flucloxacillin | 0.85 | 3 |
| staphylococcus_epidermidis | aztreonam_avibactam | 0.01 | 0.003 |
| staphylococcus_epidermidis | cefixime | 0.75 | 0.2 |
| streptococcus_pneumoniae | sulfanilamide | 0.1 | 0.02 |
| streptococcus_pneumoniae | penicillin_g | 0.95 | 6 |
| streptococcus_pneumoniae | ampicillin | 0.95 | 6 |
| streptococcus_pneumoniae | amoxicillin | 0.95 | 6 |
| streptococcus_pneumoniae | piperacillin | 0.9 | 1 |
| streptococcus_pneumoniae | ticarcillin | 0.9 | 1 |
| streptococcus_pneumoniae | cephalexin | 0.85 | 0.3 |
| streptococcus_pneumoniae | cefazolin | 0.9 | 0.3 |
| streptococcus_pneumoniae | cefuroxime | 0.9 | 3 |
| streptococcus_pneumoniae | ceftriaxone | 0.95 | 3 |
| streptococcus_pneumoniae | ceftazidime | 0.7 | 0.2 |
| streptococcus_pneumoniae | cefepime | 0.8 | 0.35 |
| streptococcus_pneumoniae | ceftaroline | 0.95 | 0.05 |
| streptococcus_pneumoniae | ceftolozane_tazobactam | 0.75 | 1 |
| streptococcus_pneumoniae | cefiderocol | 0.75 | 1 |
| streptococcus_pneumoniae | meropenem | 0.95 | 0.005 |
| streptococcus_pneumoniae | imipenem_c | 0.95 | 0.005 |
| streptococcus_pneumoniae | ertapenem | 0.95 | 0.005 |
| streptococcus_pneumoniae | aztreonam | 0 | 0.003 |
| streptococcus_pneumoniae | erythromycin | 0.8 | 5 |
| streptococcus_pneumoniae | azithromycin | 0.85 | 7 |
| streptococcus_pneumoniae | clarithromycin | 0.85 | 7 |
| streptococcus_pneumoniae | clindamycin | 0.8 | 1 |
| streptococcus_pneumoniae | gentamicin | 0.1 | 1 |
| streptococcus_pneumoniae | tobramycin | 0.1 | 1 |
| streptococcus_pneumoniae | amikacin | 0.1 | 1 |
| streptococcus_pneumoniae | ciprofloxacin | 0.9 | 1 |
| streptococcus_pneumoniae | levofloxacin | 0.95 | 5 |
| streptococcus_pneumoniae | moxifloxacin | 0.95 | 5 |
| streptococcus_pneumoniae | ofloxacin | 0.9 | 1 |
| streptococcus_pneumoniae | tetracycline | 0.8 | 0.25 |
| streptococcus_pneumoniae | doxycycline | 0.85 | 0.25 |
| streptococcus_pneumoniae | minocycline | 0.85 | 0.25 |
| streptococcus_pneumoniae | tigecycline | 0.1 | 1 |
| streptococcus_pneumoniae | vancomycin | 0.95 | 3 |
| streptococcus_pneumoniae | teicoplanin | 0.9 | 1 |
| streptococcus_pneumoniae | dalbavancin | 0.9 | 0.005 |
| streptococcus_pneumoniae | linezolid | 0.9 | 0.005 |
| streptococcus_pneumoniae | tedizolid | 0.9 | 0.005 |
| streptococcus_pneumoniae | daptomycin | 0.1 | 1 |
| streptococcus_pneumoniae | quinu_dalfo | 0.85 | 0.005 |
| streptococcus_pneumoniae | trim_sulf | 0.7 | 0.04 |
| streptococcus_pneumoniae | chloramphenicol | 0.8 | 1 |
| streptococcus_pneumoniae | nitrofurantoin | 0.1 | 1 |
| streptococcus_pneumoniae | fosfomycin | 0.1 | 1 |
| streptococcus_pneumoniae | retapamulin | 0.1 | 1 |
| streptococcus_pneumoniae | fusidic_a | 0.1 | 1 |
| streptococcus_pneumoniae | metronidazole | 0.1 | 1 |
| streptococcus_pneumoniae | fidaxomicin | 0.1 | 1 |
| streptococcus_pneumoniae | furazolidone | 0.1 | 1 |
| streptococcus_pneumoniae | rifampicin | 0.8 | 1 |
| streptococcus_pneumoniae | amoxicillin_clavulanate | 0.95 | 12 |
| streptococcus_pneumoniae | piperacillin_tazobactam | 0.9 | 1 |
| streptococcus_pneumoniae | ampicillin_sulbactam | 0.95 | 1 |
| streptococcus_pneumoniae | ticarcillin_clavulanate | 0.9 | 1 |
| streptococcus_pneumoniae | ceftazidime_avibactam | 0.95 | 0.005 |
| streptococcus_pneumoniae | meropenem_vaborbactam | 0.95 | 0.005 |
| streptococcus_pneumoniae | colistin | 0 | 0.005 |
| streptococcus_pneumoniae | flucloxacillin | 0.8 | 1 |
| streptococcus_pneumoniae | aztreonam_avibactam | 0.01 | 0.003 |
| streptococcus_pneumoniae | cefixime | 0.75 | 0.2 |
| salmonella_enterica_serovar_typhi | sulfanilamide | 0.7 | 0.02 |
| salmonella_enterica_serovar_typhi | penicillin_g | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | ampicillin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | amoxicillin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | piperacillin | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | ticarcillin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | cephalexin | 0.7 | 0.3 |
| salmonella_enterica_serovar_typhi | cefazolin | 0.75 | 0.3 |
| salmonella_enterica_serovar_typhi | cefuroxime | 0.8 | 0.3 |
| salmonella_enterica_serovar_typhi | ceftriaxone | 0.95 | 4 |
| salmonella_enterica_serovar_typhi | ceftazidime | 0.9 | 0.2 |
| salmonella_enterica_serovar_typhi | cefepime | 0.9 | 0.35 |
| salmonella_enterica_serovar_typhi | ceftaroline | 0.7 | 0.002 |
| salmonella_enterica_serovar_typhi | ceftolozane_tazobactam | 0.75 | 1 |
| salmonella_enterica_serovar_typhi | cefiderocol | 0.75 | 1 |
| salmonella_enterica_serovar_typhi | meropenem | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | imipenem_c | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | ertapenem | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | aztreonam | 0.9 | 0.003 |
| salmonella_enterica_serovar_typhi | erythromycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | azithromycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | clarithromycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | gentamicin | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | tobramycin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | amikacin | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | ciprofloxacin | 0.9 | 4 |
| salmonella_enterica_serovar_typhi | levofloxacin | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | moxifloxacin | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | ofloxacin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | tetracycline | 0.8 | 0.25 |
| salmonella_enterica_serovar_typhi | doxycycline | 0.85 | 0.25 |
| salmonella_enterica_serovar_typhi | minocycline | 0.85 | 0.25 |
| salmonella_enterica_serovar_typhi | tigecycline | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | dalbavancin | 0 | 0.005 |
| salmonella_enterica_serovar_typhi | linezolid | 0 | 0.005 |
| salmonella_enterica_serovar_typhi | tedizolid | 0 | 0.005 |
| salmonella_enterica_serovar_typhi | daptomycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | quinu_dalfo | 0 | 0.005 |
| salmonella_enterica_serovar_typhi | trim_sulf | 0.9 | 0.04 |
| salmonella_enterica_serovar_typhi | chloramphenicol | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | nitrofurantoin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | fosfomycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | fidaxomicin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | furazolidone | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | rifampicin | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | amoxicillin_clavulanate | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | piperacillin_tazobactam | 0.95 | 1 |
| salmonella_enterica_serovar_typhi | ampicillin_sulbactam | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | ticarcillin_clavulanate | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | ceftazidime_avibactam | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | meropenem_vaborbactam | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | colistin | 0.7 | 0.005 |
| salmonella_enterica_serovar_typhi | flucloxacillin | 0.01 | 1 |
| salmonella_enterica_serovar_typhi | aztreonam_avibactam | 0.9 | 0.003 |
| salmonella_enterica_serovar_typhi | cefixime | 0.75 | 3 |
| salmonella_enterica_serovar_paratyphi_a | sulfanilamide | 0.7 | 0.02 |
| salmonella_enterica_serovar_paratyphi_a | penicillin_g | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ampicillin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | amoxicillin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | piperacillin | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ticarcillin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | cephalexin | 0.7 | 0.3 |
| salmonella_enterica_serovar_paratyphi_a | cefazolin | 0.75 | 0.3 |
| salmonella_enterica_serovar_paratyphi_a | cefuroxime | 0.8 | 0.3 |
| salmonella_enterica_serovar_paratyphi_a | ceftriaxone | 0.95 | 4 |
| salmonella_enterica_serovar_paratyphi_a | ceftazidime | 0.9 | 0.2 |
| salmonella_enterica_serovar_paratyphi_a | cefepime | 0.9 | 0.35 |
| salmonella_enterica_serovar_paratyphi_a | ceftaroline | 0.7 | 0.002 |
| salmonella_enterica_serovar_paratyphi_a | ceftolozane_tazobactam | 0.75 | 1 |
| salmonella_enterica_serovar_paratyphi_a | cefiderocol | 0.75 | 1 |
| salmonella_enterica_serovar_paratyphi_a | meropenem | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | imipenem_c | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | ertapenem | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | aztreonam | 0.9 | 0.003 |
| salmonella_enterica_serovar_paratyphi_a | erythromycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | azithromycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | clarithromycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | gentamicin | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | tobramycin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | amikacin | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ciprofloxacin | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | levofloxacin | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | moxifloxacin | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ofloxacin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | tetracycline | 0.8 | 0.25 |
| salmonella_enterica_serovar_paratyphi_a | doxycycline | 0.85 | 0.25 |
| salmonella_enterica_serovar_paratyphi_a | minocycline | 0.85 | 0.25 |
| salmonella_enterica_serovar_paratyphi_a | tigecycline | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | dalbavancin | 0 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | linezolid | 0 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | tedizolid | 0 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | daptomycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | quinu_dalfo | 0 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | trim_sulf | 0.9 | 0.04 |
| salmonella_enterica_serovar_paratyphi_a | chloramphenicol | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | nitrofurantoin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | fosfomycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | fidaxomicin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | furazolidone | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | rifampicin | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | amoxicillin_clavulanate | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | piperacillin_tazobactam | 0.95 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ampicillin_sulbactam | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ticarcillin_clavulanate | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ceftazidime_avibactam | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | meropenem_vaborbactam | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | colistin | 0.7 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | flucloxacillin | 0.01 | 1 |
| salmonella_enterica_serovar_paratyphi_a | aztreonam_avibactam | 0.9 | 0.003 |
| salmonella_enterica_serovar_paratyphi_a | cefixime | 0.75 | 0.2 |
| invasive_non-typhoidal_salmonella_spp. | sulfanilamide | 0.7 | 0.02 |
| invasive_non-typhoidal_salmonella_spp. | penicillin_g | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ampicillin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | amoxicillin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | piperacillin | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ticarcillin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | cephalexin | 0.7 | 0.3 |
| invasive_non-typhoidal_salmonella_spp. | cefazolin | 0.75 | 0.3 |
| invasive_non-typhoidal_salmonella_spp. | cefuroxime | 0.8 | 0.3 |
| invasive_non-typhoidal_salmonella_spp. | ceftriaxone | 0.95 | 3 |
| invasive_non-typhoidal_salmonella_spp. | ceftazidime | 0.9 | 0.2 |
| invasive_non-typhoidal_salmonella_spp. | cefepime | 0.9 | 0.35 |
| invasive_non-typhoidal_salmonella_spp. | ceftaroline | 0.7 | 0.002 |
| invasive_non-typhoidal_salmonella_spp. | ceftolozane_tazobactam | 0.75 | 1 |
| invasive_non-typhoidal_salmonella_spp. | cefiderocol | 0.75 | 1 |
| invasive_non-typhoidal_salmonella_spp. | meropenem | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | imipenem_c | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | ertapenem | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | aztreonam | 0.9 | 0.003 |
| invasive_non-typhoidal_salmonella_spp. | erythromycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | azithromycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | clarithromycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | clindamycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | gentamicin | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | tobramycin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | amikacin | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ciprofloxacin | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | levofloxacin | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | moxifloxacin | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ofloxacin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | tetracycline | 0.8 | 0.25 |
| invasive_non-typhoidal_salmonella_spp. | doxycycline | 0.85 | 0.25 |
| invasive_non-typhoidal_salmonella_spp. | minocycline | 0.85 | 0.25 |
| invasive_non-typhoidal_salmonella_spp. | tigecycline | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | vancomycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | teicoplanin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | dalbavancin | 0.1 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | linezolid | 0.1 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | tedizolid | 0.1 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | daptomycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | quinu_dalfo | 0.1 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | trim_sulf | 0.9 | 0.04 |
| invasive_non-typhoidal_salmonella_spp. | chloramphenicol | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | nitrofurantoin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | fosfomycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | retapamulin | 0.05 | 1 |
| invasive_non-typhoidal_salmonella_spp. | fusidic_a | 0.05 | 1 |
| invasive_non-typhoidal_salmonella_spp. | fidaxomicin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | furazolidone | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | rifampicin | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | amoxicillin_clavulanate | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | piperacillin_tazobactam | 0.95 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ampicillin_sulbactam | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ticarcillin_clavulanate | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ceftazidime_avibactam | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | colistin | 0.7 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | flucloxacillin | 0.01 | 1 |
| invasive_non-typhoidal_salmonella_spp. | aztreonam_avibactam | 0.9 | 0.003 |
| invasive_non-typhoidal_salmonella_spp. | cefixime | 0.75 | 0.2 |
| shigella_spp. | sulfanilamide | 0.5 | 0.02 |
| shigella_spp. | penicillin_g | 0.1 | 1 |
| shigella_spp. | ampicillin | 0.7 | 1 |
| shigella_spp. | amoxicillin | 0.7 | 1 |
| shigella_spp. | piperacillin | 0.75 | 1 |
| shigella_spp. | ticarcillin | 0.7 | 1 |
| shigella_spp. | cephalexin | 0.6 | 0.3 |
| shigella_spp. | cefazolin | 0.65 | 0.3 |
| shigella_spp. | cefuroxime | 0.7 | 0.3 |
| shigella_spp. | ceftriaxone | 0.9 | 0.2 |
| shigella_spp. | ceftazidime | 0.85 | 0.2 |
| shigella_spp. | cefepime | 0.85 | 0.35 |
| shigella_spp. | ceftaroline | 0.6 | 0.002 |
| shigella_spp. | ceftolozane_tazobactam | 0.75 | 1 |
| shigella_spp. | cefiderocol | 0.75 | 1 |
| shigella_spp. | meropenem | 0.9 | 0.005 |
| shigella_spp. | imipenem_c | 0.9 | 0.005 |
| shigella_spp. | ertapenem | 0.9 | 0.005 |
| shigella_spp. | aztreonam | 0.8 | 0.003 |
| shigella_spp. | erythromycin | 0.7 | 1 |
| shigella_spp. | azithromycin | 0.85 | 1 |
| shigella_spp. | clarithromycin | 0.75 | 1 |
| shigella_spp. | gentamicin | 0.8 | 1 |
| shigella_spp. | tobramycin | 0.75 | 1 |
| shigella_spp. | amikacin | 0.85 | 1 |
| shigella_spp. | ciprofloxacin | 0.95 | 4 |
| shigella_spp. | levofloxacin | 0.9 | 1 |
| shigella_spp. | moxifloxacin | 0.8 | 1 |
| shigella_spp. | ofloxacin | 0.9 | 1 |
| shigella_spp. | tetracycline | 0.8 | 0.25 |
| shigella_spp. | doxycycline | 0.85 | 0.25 |
| shigella_spp. | minocycline | 0.85 | 0.25 |
| shigella_spp. | tigecycline | 0.7 | 1 |
| shigella_spp. | dalbavancin | 0 | 0.005 |
| shigella_spp. | linezolid | 0 | 0.005 |
| shigella_spp. | tedizolid | 0 | 0.005 |
| shigella_spp. | daptomycin | 0.1 | 1 |
| shigella_spp. | quinu_dalfo | 0 | 0.005 |
| shigella_spp. | trim_sulf | 0.9 | 0.04 |
| shigella_spp. | chloramphenicol | 0.85 | 1 |
| shigella_spp. | nitrofurantoin | 0.1 | 1 |
| shigella_spp. | fosfomycin | 0.1 | 1 |
| shigella_spp. | fidaxomicin | 0.1 | 1 |
| shigella_spp. | furazolidone | 0.1 | 1 |
| shigella_spp. | rifampicin | 0.7 | 1 |
| shigella_spp. | amoxicillin_clavulanate | 0.8 | 1 |
| shigella_spp. | piperacillin_tazobactam | 0.85 | 1 |
| shigella_spp. | ampicillin_sulbactam | 0.8 | 1 |
| shigella_spp. | ticarcillin_clavulanate | 0.85 | 1 |
| shigella_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| shigella_spp. | meropenem_vaborbactam | 0.9 | 0.005 |
| shigella_spp. | colistin | 0.7 | 0.005 |
| shigella_spp. | flucloxacillin | 0.01 | 1 |
| shigella_spp. | aztreonam_avibactam | 0.9 | 0.003 |
| shigella_spp. | cefixime | 0.75 | 0.2 |
| neisseria_gonorrhoeae | sulfanilamide | 0.1 | 0.02 |
| neisseria_gonorrhoeae | penicillin_g | 0.9 | 4 |
| neisseria_gonorrhoeae | ampicillin | 0.85 | 1 |
| neisseria_gonorrhoeae | amoxicillin | 0.85 | 2.5 |
| neisseria_gonorrhoeae | piperacillin | 0.8 | 1 |
| neisseria_gonorrhoeae | ticarcillin | 0.8 | 1 |
| neisseria_gonorrhoeae | cephalexin | 0.7 | 0.3 |
| neisseria_gonorrhoeae | cefazolin | 0.75 | 0.3 |
| neisseria_gonorrhoeae | cefuroxime | 0.85 | 0.3 |
| neisseria_gonorrhoeae | ceftriaxone | 0.95 | 6 |
| neisseria_gonorrhoeae | ceftazidime | 0.9 | 0.2 |
| neisseria_gonorrhoeae | cefepime | 0.9 | 0.35 |
| neisseria_gonorrhoeae | ceftaroline | 0.8 | 0.002 |
| neisseria_gonorrhoeae | ceftolozane_tazobactam | 0.8 | 1 |
| neisseria_gonorrhoeae | cefiderocol | 0.8 | 1 |
| neisseria_gonorrhoeae | meropenem | 0.9 | 0.005 |
| neisseria_gonorrhoeae | imipenem_c | 0.9 | 0.005 |
| neisseria_gonorrhoeae | ertapenem | 0.9 | 0.005 |
| neisseria_gonorrhoeae | aztreonam | 0.9 | 0.003 |
| neisseria_gonorrhoeae | erythromycin | 0.7 | 1 |
| neisseria_gonorrhoeae | azithromycin | 0.7 | 5 |
| neisseria_gonorrhoeae | clarithromycin | 0.7 | 1 |
| neisseria_gonorrhoeae | gentamicin | 0.7 | 2 |
| neisseria_gonorrhoeae | tobramycin | 0.7 | 1 |
| neisseria_gonorrhoeae | amikacin | 0.7 | 1 |
| neisseria_gonorrhoeae | ciprofloxacin | 0.9 | 3 |
| neisseria_gonorrhoeae | levofloxacin | 0.85 | 1 |
| neisseria_gonorrhoeae | moxifloxacin | 0.8 | 1 |
| neisseria_gonorrhoeae | ofloxacin | 0.85 | 2 |
| neisseria_gonorrhoeae | tetracycline | 0.8 | 0.25 |
| neisseria_gonorrhoeae | doxycycline | 0.9 | 0.25 |
| neisseria_gonorrhoeae | minocycline | 0.85 | 0.25 |
| neisseria_gonorrhoeae | tigecycline | 0.1 | 1 |
| neisseria_gonorrhoeae | dalbavancin | 0 | 0.005 |
| neisseria_gonorrhoeae | linezolid | 0 | 0.005 |
| neisseria_gonorrhoeae | tedizolid | 0 | 0.005 |
| neisseria_gonorrhoeae | daptomycin | 0.1 | 1 |
| neisseria_gonorrhoeae | quinu_dalfo | 0 | 0.005 |
| neisseria_gonorrhoeae | trim_sulf | 0.7 | 0.04 |
| neisseria_gonorrhoeae | chloramphenicol | 0.8 | 1 |
| neisseria_gonorrhoeae | nitrofurantoin | 0.1 | 1 |
| neisseria_gonorrhoeae | fosfomycin | 0.1 | 1 |
| neisseria_gonorrhoeae | fidaxomicin | 0.1 | 1 |
| neisseria_gonorrhoeae | furazolidone | 0.1 | 1 |
| neisseria_gonorrhoeae | rifampicin | 0.7 | 1 |
| neisseria_gonorrhoeae | amoxicillin_clavulanate | 0.85 | 1 |
| neisseria_gonorrhoeae | piperacillin_tazobactam | 0.85 | 1 |
| neisseria_gonorrhoeae | ampicillin_sulbactam | 0.8 | 1 |
| neisseria_gonorrhoeae | ticarcillin_clavulanate | 0.8 | 1 |
| neisseria_gonorrhoeae | ceftazidime_avibactam | 0.9 | 0.005 |
| neisseria_gonorrhoeae | meropenem_vaborbactam | 0.9 | 0.005 |
| neisseria_gonorrhoeae | colistin | 0.05 | 0.005 |
| neisseria_gonorrhoeae | flucloxacillin | 0.01 | 1 |
| neisseria_gonorrhoeae | aztreonam_avibactam | 0.8 | 0.003 |
| neisseria_gonorrhoeae | cefixime | 0.55 | 5 |
| streptococcus_pyogenes | sulfanilamide | 0.1 | 0.02 |
| streptococcus_pyogenes | penicillin_g | 1 | 6 |
| streptococcus_pyogenes | ampicillin | 0.95 | 6 |
| streptococcus_pyogenes | amoxicillin | 0.95 | 6 |
| streptococcus_pyogenes | piperacillin | 0.9 | 1 |
| streptococcus_pyogenes | ticarcillin | 0.9 | 1 |
| streptococcus_pyogenes | cephalexin | 0.9 | 4 |
| streptococcus_pyogenes | cefazolin | 0.9 | 0.3 |
| streptococcus_pyogenes | cefuroxime | 0.95 | 0.3 |
| streptococcus_pyogenes | ceftriaxone | 0.95 | 0.2 |
| streptococcus_pyogenes | ceftazidime | 0.7 | 0.2 |
| streptococcus_pyogenes | cefepime | 0.8 | 0.35 |
| streptococcus_pyogenes | ceftaroline | 0.95 | 0.002 |
| streptococcus_pyogenes | ceftolozane_tazobactam | 0.75 | 1 |
| streptococcus_pyogenes | cefiderocol | 0.75 | 1 |
| streptococcus_pyogenes | meropenem | 0.95 | 0.005 |
| streptococcus_pyogenes | imipenem_c | 0.95 | 0.005 |
| streptococcus_pyogenes | ertapenem | 0.95 | 0.005 |
| streptococcus_pyogenes | aztreonam | 0 | 0.003 |
| streptococcus_pyogenes | erythromycin | 0.9 | 1 |
| streptococcus_pyogenes | azithromycin | 0.9 | 5 |
| streptococcus_pyogenes | clarithromycin | 0.9 | 4.5 |
| streptococcus_pyogenes | clindamycin | 0.85 | 1 |
| streptococcus_pyogenes | gentamicin | 0.1 | 1 |
| streptococcus_pyogenes | tobramycin | 0.1 | 1 |
| streptococcus_pyogenes | amikacin | 0.1 | 1 |
| streptococcus_pyogenes | ciprofloxacin | 0.8 | 1 |
| streptococcus_pyogenes | levofloxacin | 0.9 | 1 |
| streptococcus_pyogenes | moxifloxacin | 0.9 | 1 |
| streptococcus_pyogenes | ofloxacin | 0.85 | 1 |
| streptococcus_pyogenes | tetracycline | 0.8 | 0.25 |
| streptococcus_pyogenes | doxycycline | 0.85 | 0.25 |
| streptococcus_pyogenes | minocycline | 0.85 | 0.25 |
| streptococcus_pyogenes | tigecycline | 0.1 | 1 |
| streptococcus_pyogenes | vancomycin | 0.95 | 2 |
| streptococcus_pyogenes | teicoplanin | 0.9 | 1 |
| streptococcus_pyogenes | dalbavancin | 0.9 | 0.005 |
| streptococcus_pyogenes | linezolid | 0.9 | 0.005 |
| streptococcus_pyogenes | tedizolid | 0.9 | 0.005 |
| streptococcus_pyogenes | daptomycin | 0.1 | 1 |
| streptococcus_pyogenes | quinu_dalfo | 0.85 | 0.005 |
| streptococcus_pyogenes | trim_sulf | 0.7 | 0.04 |
| streptococcus_pyogenes | chloramphenicol | 0.8 | 1 |
| streptococcus_pyogenes | nitrofurantoin | 0.1 | 1 |
| streptococcus_pyogenes | fosfomycin | 0.1 | 1 |
| streptococcus_pyogenes | retapamulin | 0.1 | 1 |
| streptococcus_pyogenes | fusidic_a | 0.1 | 1 |
| streptococcus_pyogenes | metronidazole | 0.1 | 1 |
| streptococcus_pyogenes | fidaxomicin | 0.1 | 1 |
| streptococcus_pyogenes | furazolidone | 0.1 | 1 |
| streptococcus_pyogenes | rifampicin | 0.8 | 1 |
| streptococcus_pyogenes | amoxicillin_clavulanate | 0.95 | 12 |
| streptococcus_pyogenes | piperacillin_tazobactam | 0.9 | 1 |
| streptococcus_pyogenes | ampicillin_sulbactam | 0.95 | 1 |
| streptococcus_pyogenes | ticarcillin_clavulanate | 0.9 | 1 |
| streptococcus_pyogenes | ceftazidime_avibactam | 0.95 | 0.005 |
| streptococcus_pyogenes | meropenem_vaborbactam | 0.95 | 0.005 |
| streptococcus_pyogenes | colistin | 0 | 0.005 |
| streptococcus_pyogenes | flucloxacillin | 0.8 | 1 |
| streptococcus_pyogenes | aztreonam_avibactam | 0.01 | 0.003 |
| streptococcus_pyogenes | cefixime | 0.75 | 0.2 |
| streptococcus_agalactiae | sulfanilamide | 0.1 | 0.02 |
| streptococcus_agalactiae | penicillin_g | 0.95 | 6 |
| streptococcus_agalactiae | ampicillin | 0.95 | 6 |
| streptococcus_agalactiae | amoxicillin | 0.95 | 6 |
| streptococcus_agalactiae | piperacillin | 0.9 | 1 |
| streptococcus_agalactiae | ticarcillin | 0.9 | 1 |
| streptococcus_agalactiae | cephalexin | 0.9 | 4 |
| streptococcus_agalactiae | cefazolin | 0.9 | 0.3 |
| streptococcus_agalactiae | cefuroxime | 0.95 | 0.3 |
| streptococcus_agalactiae | ceftriaxone | 0.95 | 0.2 |
| streptococcus_agalactiae | ceftazidime | 0.7 | 0.2 |
| streptococcus_agalactiae | cefepime | 0.8 | 0.35 |
| streptococcus_agalactiae | ceftaroline | 0.95 | 0.002 |
| streptococcus_agalactiae | ceftolozane_tazobactam | 0.75 | 1 |
| streptococcus_agalactiae | cefiderocol | 0.75 | 1 |
| streptococcus_agalactiae | meropenem | 0.95 | 0.005 |
| streptococcus_agalactiae | imipenem_c | 0.95 | 0.005 |
| streptococcus_agalactiae | ertapenem | 0.95 | 0.005 |
| streptococcus_agalactiae | aztreonam | 0 | 0.003 |
| streptococcus_agalactiae | erythromycin | 0.8 | 1 |
| streptococcus_agalactiae | azithromycin | 0.85 | 1 |
| streptococcus_agalactiae | clarithromycin | 0.85 | 1 |
| streptococcus_agalactiae | clindamycin | 0.8 | 1 |
| streptococcus_agalactiae | gentamicin | 0.1 | 1 |
| streptococcus_agalactiae | tobramycin | 0.1 | 1 |
| streptococcus_agalactiae | amikacin | 0.1 | 1 |
| streptococcus_agalactiae | ciprofloxacin | 0.8 | 1 |
| streptococcus_agalactiae | levofloxacin | 0.9 | 1 |
| streptococcus_agalactiae | moxifloxacin | 0.9 | 1 |
| streptococcus_agalactiae | ofloxacin | 0.85 | 1 |
| streptococcus_agalactiae | tetracycline | 0.8 | 0.25 |
| streptococcus_agalactiae | doxycycline | 0.85 | 0.25 |
| streptococcus_agalactiae | minocycline | 0.85 | 0.25 |
| streptococcus_agalactiae | tigecycline | 0.1 | 1 |
| streptococcus_agalactiae | vancomycin | 0.95 | 2 |
| streptococcus_agalactiae | teicoplanin | 0.9 | 1 |
| streptococcus_agalactiae | dalbavancin | 0.9 | 0.005 |
| streptococcus_agalactiae | linezolid | 0.9 | 0.005 |
| streptococcus_agalactiae | tedizolid | 0.9 | 0.005 |
| streptococcus_agalactiae | daptomycin | 0.1 | 1 |
| streptococcus_agalactiae | quinu_dalfo | 0.85 | 0.005 |
| streptococcus_agalactiae | trim_sulf | 0.7 | 0.04 |
| streptococcus_agalactiae | chloramphenicol | 0.8 | 1 |
| streptococcus_agalactiae | nitrofurantoin | 0.1 | 1 |
| streptococcus_agalactiae | fosfomycin | 0.1 | 1 |
| streptococcus_agalactiae | retapamulin | 0.1 | 1 |
| streptococcus_agalactiae | fusidic_a | 0.1 | 1 |
| streptococcus_agalactiae | metronidazole | 0.1 | 1 |
| streptococcus_agalactiae | fidaxomicin | 0.1 | 1 |
| streptococcus_agalactiae | furazolidone | 0.1 | 1 |
| streptococcus_agalactiae | rifampicin | 0.8 | 1 |
| streptococcus_agalactiae | amoxicillin_clavulanate | 0.95 | 12 |
| streptococcus_agalactiae | piperacillin_tazobactam | 0.9 | 1 |
| streptococcus_agalactiae | ampicillin_sulbactam | 0.95 | 1 |
| streptococcus_agalactiae | ticarcillin_clavulanate | 0.9 | 1 |
| streptococcus_agalactiae | ceftazidime_avibactam | 0.95 | 0.005 |
| streptococcus_agalactiae | meropenem_vaborbactam | 0.95 | 0.005 |
| streptococcus_agalactiae | colistin | 0 | 0.005 |
| streptococcus_agalactiae | flucloxacillin | 0.8 | 1 |
| streptococcus_agalactiae | aztreonam_avibactam | 0.01 | 0.003 |
| streptococcus_agalactiae | cefixime | 0.75 | 0.2 |
| haemophilus_influenzae | sulfanilamide | 0.1 | 0.02 |
| haemophilus_influenzae | penicillin_g | 0.7 | 1 |
| haemophilus_influenzae | ampicillin | 0.8 | 6 |
| haemophilus_influenzae | amoxicillin | 0.9 | 6 |
| haemophilus_influenzae | piperacillin | 0.85 | 1 |
| haemophilus_influenzae | ticarcillin | 0.8 | 1 |
| haemophilus_influenzae | cephalexin | 0.7 | 0.3 |
| haemophilus_influenzae | cefazolin | 0.75 | 0.3 |
| haemophilus_influenzae | cefuroxime | 0.85 | 3 |
| haemophilus_influenzae | ceftriaxone | 0.95 | 3 |
| haemophilus_influenzae | ceftazidime | 0.9 | 0.2 |
| haemophilus_influenzae | cefepime | 0.9 | 0.35 |
| haemophilus_influenzae | ceftaroline | 0.8 | 0.002 |
| haemophilus_influenzae | ceftolozane_tazobactam | 0.8 | 1 |
| haemophilus_influenzae | cefiderocol | 0.8 | 1 |
| haemophilus_influenzae | meropenem | 0.95 | 0.005 |
| haemophilus_influenzae | imipenem_c | 0.95 | 0.005 |
| haemophilus_influenzae | ertapenem | 0.95 | 0.005 |
| haemophilus_influenzae | aztreonam | 0.9 | 0.003 |
| haemophilus_influenzae | erythromycin | 0.7 | 6 |
| haemophilus_influenzae | azithromycin | 0.9 | 7 |
| haemophilus_influenzae | clarithromycin | 0.85 | 7 |
| haemophilus_influenzae | gentamicin | 0.7 | 1 |
| haemophilus_influenzae | tobramycin | 0.7 | 1 |
| haemophilus_influenzae | amikacin | 0.7 | 1 |
| haemophilus_influenzae | ciprofloxacin | 0.9 | 1 |
| haemophilus_influenzae | levofloxacin | 0.85 | 4 |
| haemophilus_influenzae | moxifloxacin | 0.8 | 4 |
| haemophilus_influenzae | ofloxacin | 0.85 | 1 |
| haemophilus_influenzae | tetracycline | 0.85 | 0.25 |
| haemophilus_influenzae | doxycycline | 0.85 | 0.25 |
| haemophilus_influenzae | minocycline | 0.85 | 0.25 |
| haemophilus_influenzae | tigecycline | 0.1 | 1 |
| haemophilus_influenzae | dalbavancin | 0 | 0.005 |
| haemophilus_influenzae | linezolid | 0 | 0.005 |
| haemophilus_influenzae | tedizolid | 0 | 0.005 |
| haemophilus_influenzae | daptomycin | 0.1 | 1 |
| haemophilus_influenzae | quinu_dalfo | 0 | 0.005 |
| haemophilus_influenzae | trim_sulf | 0.85 | 0.04 |
| haemophilus_influenzae | chloramphenicol | 0.8 | 1 |
| haemophilus_influenzae | nitrofurantoin | 0.1 | 1 |
| haemophilus_influenzae | fosfomycin | 0.1 | 1 |
| haemophilus_influenzae | fidaxomicin | 0.1 | 1 |
| haemophilus_influenzae | furazolidone | 0.1 | 1 |
| haemophilus_influenzae | rifampicin | 0.7 | 1 |
| haemophilus_influenzae | amoxicillin_clavulanate | 0.9 | 12 |
| haemophilus_influenzae | piperacillin_tazobactam | 0.85 | 8 |
| haemophilus_influenzae | ampicillin_sulbactam | 0.9 | 1 |
| haemophilus_influenzae | ticarcillin_clavulanate | 0.8 | 1 |
| haemophilus_influenzae | ceftazidime_avibactam | 0.95 | 0.005 |
| haemophilus_influenzae | meropenem_vaborbactam | 0.95 | 0.005 |
| haemophilus_influenzae | colistin | 0.05 | 0.005 |
| haemophilus_influenzae | flucloxacillin | 0.01 | 1 |
| haemophilus_influenzae | aztreonam_avibactam | 0.8 | 0.003 |
| haemophilus_influenzae | cefixime | 0.8 | 0.2 |
| chlamydia_trachomatis | sulfanilamide | 0.1 | 0.02 |
| chlamydia_trachomatis | penicillin_g | 0.1 | 1 |
| chlamydia_trachomatis | ampicillin | 0.1 | 1 |
| chlamydia_trachomatis | amoxicillin | 0.1 | 1 |
| chlamydia_trachomatis | piperacillin | 0.1 | 1 |
| chlamydia_trachomatis | ticarcillin | 0.1 | 1 |
| chlamydia_trachomatis | cephalexin | 0.1 | 0.3 |
| chlamydia_trachomatis | cefazolin | 0.1 | 0.3 |
| chlamydia_trachomatis | cefuroxime | 0.1 | 0.3 |
| chlamydia_trachomatis | ceftriaxone | 0.1 | 0.2 |
| chlamydia_trachomatis | ceftazidime | 0.1 | 0.2 |
| chlamydia_trachomatis | cefepime | 0.1 | 0.35 |
| chlamydia_trachomatis | ceftaroline | 0.1 | 0.002 |
| chlamydia_trachomatis | ceftolozane_tazobactam | 0.01 | 1 |
| chlamydia_trachomatis | cefiderocol | 0.01 | 1 |
| chlamydia_trachomatis | meropenem | 0.1 | 0.005 |
| chlamydia_trachomatis | imipenem_c | 0.1 | 0.005 |
| chlamydia_trachomatis | ertapenem | 0.1 | 0.005 |
| chlamydia_trachomatis | aztreonam | 0.1 | 0.003 |
| chlamydia_trachomatis | erythromycin | 0.8 | 1 |
| chlamydia_trachomatis | azithromycin | 0.95 | 8 |
| chlamydia_trachomatis | clarithromycin | 0.9 | 1 |
| chlamydia_trachomatis | clindamycin | 0.7 | 1 |
| chlamydia_trachomatis | gentamicin | 0.1 | 1 |
| chlamydia_trachomatis | tobramycin | 0.1 | 1 |
| chlamydia_trachomatis | amikacin | 0.1 | 1 |
| chlamydia_trachomatis | ciprofloxacin | 0.8 | 1 |
| chlamydia_trachomatis | levofloxacin | 0.85 | 1 |
| chlamydia_trachomatis | moxifloxacin | 0.85 | 1 |
| chlamydia_trachomatis | ofloxacin | 0.8 | 3 |
| chlamydia_trachomatis | tetracycline | 0.95 | 2 |
| chlamydia_trachomatis | doxycycline | 0.95 | 1.5 |
| chlamydia_trachomatis | minocycline | 0.9 | 0.25 |
| chlamydia_trachomatis | tigecycline | 0.85 | 1 |
| chlamydia_trachomatis | vancomycin | 0.1 | 1 |
| chlamydia_trachomatis | teicoplanin | 0.1 | 1 |
| chlamydia_trachomatis | dalbavancin | 0.1 | 0.005 |
| chlamydia_trachomatis | linezolid | 0.1 | 0.005 |
| chlamydia_trachomatis | tedizolid | 0.1 | 0.005 |
| chlamydia_trachomatis | daptomycin | 0.1 | 1 |
| chlamydia_trachomatis | quinu_dalfo | 0.1 | 0.005 |
| chlamydia_trachomatis | trim_sulf | 0.1 | 0.04 |
| chlamydia_trachomatis | chloramphenicol | 0.8 | 1 |
| chlamydia_trachomatis | nitrofurantoin | 0.1 | 1 |
| chlamydia_trachomatis | fosfomycin | 0.1 | 1 |
| chlamydia_trachomatis | retapamulin | 0.1 | 1 |
| chlamydia_trachomatis | fusidic_a | 0.1 | 1 |
| chlamydia_trachomatis | metronidazole | 0.1 | 1 |
| chlamydia_trachomatis | fidaxomicin | 0.1 | 1 |
| chlamydia_trachomatis | furazolidone | 0.1 | 1 |
| chlamydia_trachomatis | rifampicin | 0.1 | 1 |
| chlamydia_trachomatis | amoxicillin_clavulanate | 0.1 | 1 |
| chlamydia_trachomatis | piperacillin_tazobactam | 0.1 | 1 |
| chlamydia_trachomatis | ampicillin_sulbactam | 0.1 | 1 |
| chlamydia_trachomatis | ticarcillin_clavulanate | 0.1 | 1 |
| chlamydia_trachomatis | ceftazidime_avibactam | 0.1 | 0.005 |
| chlamydia_trachomatis | meropenem_vaborbactam | 0.1 | 0.005 |
| chlamydia_trachomatis | colistin | 0.1 | 0.005 |
| chlamydia_trachomatis | flucloxacillin | 0.01 | 1 |
| chlamydia_trachomatis | aztreonam_avibactam | 0.01 | 0.003 |
| chlamydia_trachomatis | cefixime | 0.01 | 0.2 |
| mycoplasma_genitalium | sulfanilamide | 0.05 | 0.02 |
| mycoplasma_genitalium | penicillin_g | 0.05 | 1 |
| mycoplasma_genitalium | ampicillin | 0.05 | 1 |
| mycoplasma_genitalium | amoxicillin | 0.05 | 1 |
| mycoplasma_genitalium | piperacillin | 0.05 | 1 |
| mycoplasma_genitalium | ticarcillin | 0.05 | 1 |
| mycoplasma_genitalium | cephalexin | 0.05 | 0.3 |
| mycoplasma_genitalium | cefazolin | 0.05 | 0.3 |
| mycoplasma_genitalium | cefuroxime | 0.05 | 0.3 |
| mycoplasma_genitalium | ceftriaxone | 0.05 | 0.2 |
| mycoplasma_genitalium | ceftazidime | 0.05 | 0.2 |
| mycoplasma_genitalium | cefepime | 0.05 | 0.35 |
| mycoplasma_genitalium | ceftaroline | 0.05 | 0.002 |
| mycoplasma_genitalium | ceftolozane_tazobactam | 0.01 | 1 |
| mycoplasma_genitalium | cefiderocol | 0.01 | 1 |
| mycoplasma_genitalium | meropenem | 0.05 | 0.005 |
| mycoplasma_genitalium | imipenem_c | 0.05 | 0.005 |
| mycoplasma_genitalium | ertapenem | 0.05 | 0.005 |
| mycoplasma_genitalium | aztreonam | 0.05 | 0.003 |
| mycoplasma_genitalium | erythromycin | 0.8 | 1 |
| mycoplasma_genitalium | azithromycin | 0.9 | 5 |
| mycoplasma_genitalium | clarithromycin | 0.9 | 1 |
| mycoplasma_genitalium | clindamycin | 0.2 | 1 |
| mycoplasma_genitalium | gentamicin | 0.05 | 1 |
| mycoplasma_genitalium | tobramycin | 0.05 | 1 |
| mycoplasma_genitalium | amikacin | 0.05 | 1 |
| mycoplasma_genitalium | ciprofloxacin | 0.3 | 1 |
| mycoplasma_genitalium | levofloxacin | 0.5 | 2.5 |
| mycoplasma_genitalium | moxifloxacin | 0.85 | 4 |
| mycoplasma_genitalium | ofloxacin | 0.45 | 1 |
| mycoplasma_genitalium | tetracycline | 0.4 | 0.25 |
| mycoplasma_genitalium | doxycycline | 0.6 | 1.5 |
| mycoplasma_genitalium | minocycline | 0.7 | 0.25 |
| mycoplasma_genitalium | tigecycline | 0.85 | 1 |
| mycoplasma_genitalium | vancomycin | 0.05 | 1 |
| mycoplasma_genitalium | teicoplanin | 0.05 | 1 |
| mycoplasma_genitalium | dalbavancin | 0.05 | 0.005 |
| mycoplasma_genitalium | linezolid | 0.05 | 0.005 |
| mycoplasma_genitalium | tedizolid | 0.05 | 0.005 |
| mycoplasma_genitalium | daptomycin | 0.1 | 1 |
| mycoplasma_genitalium | quinu_dalfo | 0.05 | 0.005 |
| mycoplasma_genitalium | trim_sulf | 0.05 | 0.04 |
| mycoplasma_genitalium | chloramphenicol | 0.2 | 1 |
| mycoplasma_genitalium | nitrofurantoin | 0.05 | 1 |
| mycoplasma_genitalium | fosfomycin | 0.1 | 1 |
| mycoplasma_genitalium | retapamulin | 0.05 | 1 |
| mycoplasma_genitalium | fusidic_a | 0.05 | 1 |
| mycoplasma_genitalium | metronidazole | 0.05 | 1 |
| mycoplasma_genitalium | fidaxomicin | 0.1 | 1 |
| mycoplasma_genitalium | furazolidone | 0.05 | 1 |
| mycoplasma_genitalium | rifampicin | 0.1 | 1 |
| mycoplasma_genitalium | amoxicillin_clavulanate | 0.05 | 1 |
| mycoplasma_genitalium | piperacillin_tazobactam | 0.05 | 1 |
| mycoplasma_genitalium | ampicillin_sulbactam | 0.05 | 1 |
| mycoplasma_genitalium | ticarcillin_clavulanate | 0.05 | 1 |
| mycoplasma_genitalium | ceftazidime_avibactam | 0.05 | 0.005 |
| mycoplasma_genitalium | meropenem_vaborbactam | 0.05 | 0.005 |
| mycoplasma_genitalium | colistin | 0.05 | 0.005 |
| mycoplasma_genitalium | flucloxacillin | 0.01 | 1 |
| mycoplasma_genitalium | aztreonam_avibactam | 0.01 | 0.003 |
| mycoplasma_genitalium | cefixime | 0.01 | 0.2 |
| vibrio_cholerae | sulfanilamide | 0.5 | 0.02 |
| vibrio_cholerae | penicillin_g | 0.7 | 1 |
| vibrio_cholerae | ampicillin | 0.8 | 1 |
| vibrio_cholerae | amoxicillin | 0.8 | 1 |
| vibrio_cholerae | piperacillin | 0.85 | 1 |
| vibrio_cholerae | ticarcillin | 0.8 | 1 |
| vibrio_cholerae | cephalexin | 0.7 | 0.3 |
| vibrio_cholerae | cefazolin | 0.75 | 0.3 |
| vibrio_cholerae | cefuroxime | 0.8 | 0.3 |
| vibrio_cholerae | ceftriaxone | 0.9 | 0.2 |
| vibrio_cholerae | ceftazidime | 0.85 | 0.2 |
| vibrio_cholerae | cefepime | 0.85 | 0.35 |
| vibrio_cholerae | ceftaroline | 0.7 | 0.002 |
| vibrio_cholerae | ceftolozane_tazobactam | 0.75 | 1 |
| vibrio_cholerae | cefiderocol | 0.75 | 1 |
| vibrio_cholerae | meropenem | 0.9 | 0.005 |
| vibrio_cholerae | imipenem_c | 0.9 | 0.005 |
| vibrio_cholerae | ertapenem | 0.9 | 0.005 |
| vibrio_cholerae | aztreonam | 0.8 | 0.003 |
| vibrio_cholerae | erythromycin | 0.7 | 1 |
| vibrio_cholerae | azithromycin | 0.8 | 1 |
| vibrio_cholerae | clarithromycin | 0.75 | 1 |
| vibrio_cholerae | clindamycin | 0.1 | 1 |
| vibrio_cholerae | gentamicin | 0.85 | 1 |
| vibrio_cholerae | tobramycin | 0.8 | 1 |
| vibrio_cholerae | amikacin | 0.85 | 1 |
| vibrio_cholerae | ciprofloxacin | 0.9 | 1 |
| vibrio_cholerae | levofloxacin | 0.85 | 1 |
| vibrio_cholerae | moxifloxacin | 0.75 | 1 |
| vibrio_cholerae | ofloxacin | 0.85 | 1 |
| vibrio_cholerae | tetracycline | 0.95 | 0.25 |
| vibrio_cholerae | doxycycline | 0.95 | 0.25 |
| vibrio_cholerae | minocycline | 0.9 | 0.25 |
| vibrio_cholerae | tigecycline | 0.7 | 1 |
| vibrio_cholerae | vancomycin | 0.1 | 1 |
| vibrio_cholerae | teicoplanin | 0.1 | 1 |
| vibrio_cholerae | dalbavancin | 0.1 | 0.005 |
| vibrio_cholerae | linezolid | 0.1 | 0.005 |
| vibrio_cholerae | tedizolid | 0.1 | 0.005 |
| vibrio_cholerae | daptomycin | 0.1 | 1 |
| vibrio_cholerae | quinu_dalfo | 0.1 | 0.005 |
| vibrio_cholerae | trim_sulf | 0.8 | 0.04 |
| vibrio_cholerae | chloramphenicol | 0.8 | 1 |
| vibrio_cholerae | nitrofurantoin | 0.1 | 1 |
| vibrio_cholerae | fosfomycin | 0.1 | 1 |
| vibrio_cholerae | retapamulin | 0.05 | 1 |
| vibrio_cholerae | fusidic_a | 0.05 | 1 |
| vibrio_cholerae | metronidazole | 0.05 | 1 |
| vibrio_cholerae | fidaxomicin | 0.1 | 1 |
| vibrio_cholerae | furazolidone | 0.1 | 1 |
| vibrio_cholerae | rifampicin | 0.7 | 1 |
| vibrio_cholerae | amoxicillin_clavulanate | 0.85 | 1 |
| vibrio_cholerae | piperacillin_tazobactam | 0.9 | 1 |
| vibrio_cholerae | ampicillin_sulbactam | 0.85 | 1 |
| vibrio_cholerae | ticarcillin_clavulanate | 0.85 | 1 |
| vibrio_cholerae | ceftazidime_avibactam | 0.9 | 0.005 |
| vibrio_cholerae | meropenem_vaborbactam | 0.9 | 0.005 |
| vibrio_cholerae | colistin | 0.7 | 0.005 |
| vibrio_cholerae | flucloxacillin | 0.01 | 1 |
| vibrio_cholerae | aztreonam_avibactam | 0.9 | 0.003 |
| vibrio_cholerae | cefixime | 0.75 | 0.2 |
| neisseria_meningitidis | sulfanilamide | 0.1 | 0.02 |
| neisseria_meningitidis | penicillin_g | 0.95 | 6 |
| neisseria_meningitidis | ampicillin | 0.9 | 6 |
| neisseria_meningitidis | amoxicillin | 0.9 | 1 |
| neisseria_meningitidis | piperacillin | 0.85 | 1 |
| neisseria_meningitidis | ticarcillin | 0.8 | 1 |
| neisseria_meningitidis | cephalexin | 0.8 | 0.3 |
| neisseria_meningitidis | cefazolin | 0.85 | 0.3 |
| neisseria_meningitidis | cefuroxime | 0.9 | 0.3 |
| neisseria_meningitidis | ceftriaxone | 0.95 | 5 |
| neisseria_meningitidis | ceftazidime | 0.9 | 0.2 |
| neisseria_meningitidis | cefepime | 0.9 | 0.35 |
| neisseria_meningitidis | ceftaroline | 0.8 | 0.002 |
| neisseria_meningitidis | ceftolozane_tazobactam | 0.8 | 1 |
| neisseria_meningitidis | cefiderocol | 0.8 | 1 |
| neisseria_meningitidis | meropenem | 0.95 | 0.005 |
| neisseria_meningitidis | imipenem_c | 0.95 | 0.005 |
| neisseria_meningitidis | ertapenem | 0.95 | 0.005 |
| neisseria_meningitidis | aztreonam | 0.9 | 0.003 |
| neisseria_meningitidis | erythromycin | 0.7 | 1 |
| neisseria_meningitidis | azithromycin | 0.8 | 1 |
| neisseria_meningitidis | clarithromycin | 0.75 | 1 |
| neisseria_meningitidis | clindamycin | 0.1 | 1 |
| neisseria_meningitidis | gentamicin | 0.7 | 1 |
| neisseria_meningitidis | tobramycin | 0.7 | 1 |
| neisseria_meningitidis | amikacin | 0.7 | 1 |
| neisseria_meningitidis | ciprofloxacin | 0.9 | 15 |
| neisseria_meningitidis | levofloxacin | 0.85 | 1 |
| neisseria_meningitidis | moxifloxacin | 0.8 | 1 |
| neisseria_meningitidis | ofloxacin | 0.85 | 1 |
| neisseria_meningitidis | tetracycline | 0.8 | 0.25 |
| neisseria_meningitidis | doxycycline | 0.8 | 0.25 |
| neisseria_meningitidis | minocycline | 0.85 | 0.25 |
| neisseria_meningitidis | tigecycline | 0.1 | 1 |
| neisseria_meningitidis | vancomycin | 0.1 | 1 |
| neisseria_meningitidis | teicoplanin | 0.1 | 1 |
| neisseria_meningitidis | dalbavancin | 0.1 | 0.005 |
| neisseria_meningitidis | linezolid | 0.1 | 0.005 |
| neisseria_meningitidis | tedizolid | 0.1 | 0.005 |
| neisseria_meningitidis | daptomycin | 0.1 | 1 |
| neisseria_meningitidis | quinu_dalfo | 0.1 | 0.005 |
| neisseria_meningitidis | trim_sulf | 0.7 | 0.04 |
| neisseria_meningitidis | chloramphenicol | 0.85 | 18 |
| neisseria_meningitidis | nitrofurantoin | 0.1 | 1 |
| neisseria_meningitidis | fosfomycin | 0.1 | 1 |
| neisseria_meningitidis | retapamulin | 0.05 | 1 |
| neisseria_meningitidis | fusidic_a | 0.05 | 1 |
| neisseria_meningitidis | metronidazole | 0.05 | 1 |
| neisseria_meningitidis | fidaxomicin | 0.1 | 1 |
| neisseria_meningitidis | furazolidone | 0.1 | 1 |
| neisseria_meningitidis | rifampicin | 0.85 | 12 |
| neisseria_meningitidis | amoxicillin_clavulanate | 0.9 | 1 |
| neisseria_meningitidis | piperacillin_tazobactam | 0.85 | 1 |
| neisseria_meningitidis | ampicillin_sulbactam | 0.9 | 1 |
| neisseria_meningitidis | ticarcillin_clavulanate | 0.85 | 1 |
| neisseria_meningitidis | ceftazidime_avibactam | 0.95 | 0.005 |
| neisseria_meningitidis | meropenem_vaborbactam | 0.95 | 0.005 |
| neisseria_meningitidis | colistin | 0.05 | 0.005 |
| neisseria_meningitidis | flucloxacillin | 0.01 | 1 |
| neisseria_meningitidis | aztreonam_avibactam | 0.8 | 0.003 |
| neisseria_meningitidis | cefixime | 0.8 | 0.2 |
| listeria_monocytogenes | sulfanilamide | 0.1 | 0.02 |
| listeria_monocytogenes | penicillin_g | 0.7 | 1 |
| listeria_monocytogenes | ampicillin | 0.95 | 6 |
| listeria_monocytogenes | amoxicillin | 0.95 | 1 |
| listeria_monocytogenes | piperacillin | 0.7 | 1 |
| listeria_monocytogenes | ticarcillin | 0.6 | 1 |
| listeria_monocytogenes | cephalexin | 0.1 | 0.3 |
| listeria_monocytogenes | cefazolin | 0.1 | 0.3 |
| listeria_monocytogenes | cefuroxime | 0.1 | 0.3 |
| listeria_monocytogenes | ceftriaxone | 0.1 | 0.2 |
| listeria_monocytogenes | ceftazidime | 0.1 | 0.2 |
| listeria_monocytogenes | cefepime | 0.1 | 0.35 |
| listeria_monocytogenes | ceftaroline | 0.1 | 0.002 |
| listeria_monocytogenes | ceftolozane_tazobactam | 0.05 | 1 |
| listeria_monocytogenes | cefiderocol | 0.05 | 1 |
| listeria_monocytogenes | meropenem | 0.7 | 0.005 |
| listeria_monocytogenes | imipenem_c | 0.7 | 0.005 |
| listeria_monocytogenes | ertapenem | 0.7 | 0.005 |
| listeria_monocytogenes | aztreonam | 0.1 | 0.003 |
| listeria_monocytogenes | erythromycin | 0.8 | 1 |
| listeria_monocytogenes | azithromycin | 0.85 | 1 |
| listeria_monocytogenes | clarithromycin | 0.8 | 1 |
| listeria_monocytogenes | clindamycin | 0.1 | 1 |
| listeria_monocytogenes | gentamicin | 0.1 | 1 |
| listeria_monocytogenes | tobramycin | 0.1 | 1 |
| listeria_monocytogenes | amikacin | 0.1 | 1 |
| listeria_monocytogenes | ciprofloxacin | 0.8 | 1 |
| listeria_monocytogenes | levofloxacin | 0.85 | 1 |
| listeria_monocytogenes | moxifloxacin | 0.8 | 1 |
| listeria_monocytogenes | ofloxacin | 0.8 | 1 |
| listeria_monocytogenes | tetracycline | 0.8 | 0.25 |
| listeria_monocytogenes | doxycycline | 0.85 | 0.25 |
| listeria_monocytogenes | minocycline | 0.85 | 0.25 |
| listeria_monocytogenes | tigecycline | 0.1 | 1 |
| listeria_monocytogenes | vancomycin | 0.1 | 1 |
| listeria_monocytogenes | teicoplanin | 0.1 | 1 |
| listeria_monocytogenes | dalbavancin | 0.1 | 0.005 |
| listeria_monocytogenes | linezolid | 0.1 | 0.005 |
| listeria_monocytogenes | tedizolid | 0.1 | 0.005 |
| listeria_monocytogenes | daptomycin | 0.1 | 1 |
| listeria_monocytogenes | quinu_dalfo | 0.1 | 0.005 |
| listeria_monocytogenes | trim_sulf | 0.9 | 1.5 |
| listeria_monocytogenes | chloramphenicol | 0.85 | 1 |
| listeria_monocytogenes | nitrofurantoin | 0.1 | 1 |
| listeria_monocytogenes | fosfomycin | 0.1 | 1 |
| listeria_monocytogenes | retapamulin | 0.1 | 1 |
| listeria_monocytogenes | fusidic_a | 0.1 | 1 |
| listeria_monocytogenes | metronidazole | 0.1 | 1 |
| listeria_monocytogenes | fidaxomicin | 0.1 | 1 |
| listeria_monocytogenes | furazolidone | 0.1 | 1 |
| listeria_monocytogenes | rifampicin | 0.8 | 1 |
| listeria_monocytogenes | amoxicillin_clavulanate | 0.7 | 6 |
| listeria_monocytogenes | piperacillin_tazobactam | 0.95 | 1 |
| listeria_monocytogenes | ampicillin_sulbactam | 0.6 | 1 |
| listeria_monocytogenes | ticarcillin_clavulanate | 0.1 | 1 |
| listeria_monocytogenes | ceftazidime_avibactam | 0.7 | 0.005 |
| listeria_monocytogenes | meropenem_vaborbactam | 0.05 | 0.005 |
| listeria_monocytogenes | colistin | 0.1 | 0.005 |
| listeria_monocytogenes | flucloxacillin | 0.05 | 1 |
| listeria_monocytogenes | aztreonam_avibactam | 0.01 | 0.003 |
| listeria_monocytogenes | cefixime | 0.05 | 0.2 |
| clostridioides_difficile | sulfanilamide | 0.1 | 0.02 |
| clostridioides_difficile | penicillin_g | 0.1 | 1 |
| clostridioides_difficile | ampicillin | 0.1 | 1 |
| clostridioides_difficile | amoxicillin | 0.1 | 1 |
| clostridioides_difficile | piperacillin | 0.1 | 1 |
| clostridioides_difficile | ticarcillin | 0.1 | 1 |
| clostridioides_difficile | cephalexin | 0.1 | 0.3 |
| clostridioides_difficile | cefazolin | 0.1 | 0.3 |
| clostridioides_difficile | cefuroxime | 0.1 | 0.3 |
| clostridioides_difficile | ceftriaxone | 0.1 | 0.2 |
| clostridioides_difficile | ceftazidime | 0.1 | 0.2 |
| clostridioides_difficile | cefepime | 0.1 | 0.35 |
| clostridioides_difficile | ceftaroline | 0.1 | 0.002 |
| clostridioides_difficile | ceftolozane_tazobactam | 0.05 | 1 |
| clostridioides_difficile | cefiderocol | 0.05 | 1 |
| clostridioides_difficile | meropenem | 0.1 | 0.005 |
| clostridioides_difficile | imipenem_c | 0.1 | 0.005 |
| clostridioides_difficile | ertapenem | 0.1 | 0.005 |
| clostridioides_difficile | aztreonam | 0.1 | 0.003 |
| clostridioides_difficile | erythromycin | 0.7 | 1 |
| clostridioides_difficile | azithromycin | 0.75 | 1 |
| clostridioides_difficile | clarithromycin | 0.7 | 1 |
| clostridioides_difficile | clindamycin | 0.1 | 1 |
| clostridioides_difficile | gentamicin | 0.1 | 1 |
| clostridioides_difficile | tobramycin | 0.1 | 1 |
| clostridioides_difficile | amikacin | 0.1 | 1 |
| clostridioides_difficile | ciprofloxacin | 0.1 | 1 |
| clostridioides_difficile | levofloxacin | 0.1 | 1 |
| clostridioides_difficile | moxifloxacin | 0.1 | 1 |
| clostridioides_difficile | ofloxacin | 0.1 | 1 |
| clostridioides_difficile | tetracycline | 0.7 | 0.25 |
| clostridioides_difficile | doxycycline | 0.7 | 0.25 |
| clostridioides_difficile | minocycline | 0.7 | 0.25 |
| clostridioides_difficile | tigecycline | 0.1 | 1 |
| clostridioides_difficile | vancomycin | 0.95 | 5 |
| clostridioides_difficile | teicoplanin | 0.9 | 1 |
| clostridioides_difficile | dalbavancin | 0.9 | 0.005 |
| clostridioides_difficile | linezolid | 0.85 | 0.005 |
| clostridioides_difficile | tedizolid | 0.85 | 0.005 |
| clostridioides_difficile | daptomycin | 0.1 | 1 |
| clostridioides_difficile | quinu_dalfo | 0.1 | 0.005 |
| clostridioides_difficile | trim_sulf | 0.1 | 0.04 |
| clostridioides_difficile | chloramphenicol | 0.1 | 1 |
| clostridioides_difficile | nitrofurantoin | 0.1 | 1 |
| clostridioides_difficile | fosfomycin | 0.1 | 1 |
| clostridioides_difficile | retapamulin | 0.1 | 1 |
| clostridioides_difficile | fusidic_a | 0.1 | 1 |
| clostridioides_difficile | metronidazole | 0.9 | 5 |
| clostridioides_difficile | fidaxomicin | 0.1 | 1 |
| clostridioides_difficile | furazolidone | 0.1 | 1 |
| clostridioides_difficile | rifampicin | 0.1 | 1 |
| clostridioides_difficile | amoxicillin_clavulanate | 0.1 | 1 |
| clostridioides_difficile | piperacillin_tazobactam | 0.1 | 1 |
| clostridioides_difficile | ampicillin_sulbactam | 0.1 | 1 |
| clostridioides_difficile | ticarcillin_clavulanate | 0.1 | 1 |
| clostridioides_difficile | ceftazidime_avibactam | 0.1 | 0.005 |
| clostridioides_difficile | meropenem_vaborbactam | 0.1 | 0.005 |
| clostridioides_difficile | colistin | 0.05 | 0.005 |
| clostridioides_difficile | flucloxacillin | 0.01 | 1 |
| clostridioides_difficile | aztreonam_avibactam | 0.01 | 0.003 |
| clostridioides_difficile | cefixime | 0.05 | 0.2 |
| bacteroides_fragilis | sulfanilamide | 0.05 | 0.02 |
| bacteroides_fragilis | penicillin_g | 0.1 | 1 |
| bacteroides_fragilis | ampicillin | 0.2 | 1 |
| bacteroides_fragilis | amoxicillin | 0.25 | 1 |
| bacteroides_fragilis | piperacillin | 0.5 | 1 |
| bacteroides_fragilis | ticarcillin | 0.4 | 1 |
| bacteroides_fragilis | cephalexin | 0.05 | 0.3 |
| bacteroides_fragilis | cefazolin | 0.05 | 0.3 |
| bacteroides_fragilis | cefuroxime | 0.2 | 0.3 |
| bacteroides_fragilis | ceftriaxone | 0.2 | 0.2 |
| bacteroides_fragilis | ceftazidime | 0.25 | 0.2 |
| bacteroides_fragilis | cefepime | 0.25 | 0.35 |
| bacteroides_fragilis | ceftaroline | 0.2 | 0.002 |
| bacteroides_fragilis | ceftolozane_tazobactam | 0.45 | 1 |
| bacteroides_fragilis | cefiderocol | 0.45 | 1 |
| bacteroides_fragilis | meropenem | 0.95 | 0.005 |
| bacteroides_fragilis | imipenem_c | 0.95 | 0.005 |
| bacteroides_fragilis | ertapenem | 0.95 | 0.005 |
| bacteroides_fragilis | aztreonam | 0.05 | 0.003 |
| bacteroides_fragilis | erythromycin | 0.05 | 1 |
| bacteroides_fragilis | azithromycin | 0.05 | 1 |
| bacteroides_fragilis | clarithromycin | 0.05 | 1 |
| bacteroides_fragilis | clindamycin | 0.6 | 1 |
| bacteroides_fragilis | gentamicin | 0.05 | 1 |
| bacteroides_fragilis | tobramycin | 0.05 | 1 |
| bacteroides_fragilis | amikacin | 0.05 | 1 |
| bacteroides_fragilis | ciprofloxacin | 0.25 | 1 |
| bacteroides_fragilis | levofloxacin | 0.35 | 1 |
| bacteroides_fragilis | moxifloxacin | 0.5 | 1 |
| bacteroides_fragilis | ofloxacin | 0.25 | 1 |
| bacteroides_fragilis | tetracycline | 0.3 | 0.25 |
| bacteroides_fragilis | doxycycline | 0.5 | 0.25 |
| bacteroides_fragilis | minocycline | 0.5 | 0.25 |
| bacteroides_fragilis | tigecycline | 0.1 | 1 |
| bacteroides_fragilis | vancomycin | 0.05 | 1 |
| bacteroides_fragilis | teicoplanin | 0.05 | 1 |
| bacteroides_fragilis | dalbavancin | 0.05 | 0.005 |
| bacteroides_fragilis | linezolid | 0.05 | 0.005 |
| bacteroides_fragilis | tedizolid | 0.05 | 0.005 |
| bacteroides_fragilis | daptomycin | 0.1 | 1 |
| bacteroides_fragilis | quinu_dalfo | 0.05 | 0.005 |
| bacteroides_fragilis | trim_sulf | 0.3 | 0.04 |
| bacteroides_fragilis | chloramphenicol | 0.7 | 1 |
| bacteroides_fragilis | nitrofurantoin | 0.05 | 1 |
| bacteroides_fragilis | fosfomycin | 0.1 | 1 |
| bacteroides_fragilis | retapamulin | 0.05 | 1 |
| bacteroides_fragilis | fusidic_a | 0.05 | 1 |
| bacteroides_fragilis | metronidazole | 0.95 | 15 |
| bacteroides_fragilis | fidaxomicin | 0.1 | 1 |
| bacteroides_fragilis | furazolidone | 0.05 | 1 |
| bacteroides_fragilis | rifampicin | 0.2 | 1 |
| bacteroides_fragilis | amoxicillin_clavulanate | 0.75 | 6 |
| bacteroides_fragilis | piperacillin_tazobactam | 0.85 | 8 |
| bacteroides_fragilis | ampicillin_sulbactam | 0.75 | 1 |
| bacteroides_fragilis | ticarcillin_clavulanate | 0.8 | 1 |
| bacteroides_fragilis | ceftazidime_avibactam | 0.5 | 0.005 |
| bacteroides_fragilis | meropenem_vaborbactam | 0.95 | 0.005 |
| bacteroides_fragilis | colistin | 0.05 | 0.005 |
| bacteroides_fragilis | flucloxacillin | 0.01 | 1 |
| bacteroides_fragilis | aztreonam_avibactam | 0.01 | 0.003 |
| bacteroides_fragilis | cefixime | 0.45 | 0.2 |
| campylobacter_jejuni | sulfanilamide | 0.1 | 0.02 |
| campylobacter_jejuni | penicillin_g | 0.1 | 1 |
| campylobacter_jejuni | ampicillin | 0.1 | 1 |
| campylobacter_jejuni | amoxicillin | 0.1 | 1 |
| campylobacter_jejuni | piperacillin | 0.1 | 1 |
| campylobacter_jejuni | ticarcillin | 0.1 | 1 |
| campylobacter_jejuni | cephalexin | 0.1 | 0.3 |
| campylobacter_jejuni | cefazolin | 0.1 | 0.3 |
| campylobacter_jejuni | cefuroxime | 0.1 | 0.3 |
| campylobacter_jejuni | ceftriaxone | 0.1 | 0.2 |
| campylobacter_jejuni | ceftazidime | 0.1 | 0.2 |
| campylobacter_jejuni | cefepime | 0.1 | 0.35 |
| campylobacter_jejuni | ceftaroline | 0.1 | 0.002 |
| campylobacter_jejuni | ceftolozane_tazobactam | 0.75 | 1 |
| campylobacter_jejuni | cefiderocol | 0.75 | 1 |
| campylobacter_jejuni | meropenem | 0.1 | 0.005 |
| campylobacter_jejuni | imipenem_c | 0.1 | 0.005 |
| campylobacter_jejuni | ertapenem | 0.1 | 0.005 |
| campylobacter_jejuni | aztreonam | 0.1 | 0.003 |
| campylobacter_jejuni | erythromycin | 0.85 | 5 |
| campylobacter_jejuni | azithromycin | 0.9 | 5 |
| campylobacter_jejuni | clarithromycin | 0.85 | 1 |
| campylobacter_jejuni | clindamycin | 0.7 | 1 |
| campylobacter_jejuni | gentamicin | 0.7 | 1 |
| campylobacter_jejuni | tobramycin | 0.7 | 1 |
| campylobacter_jejuni | amikacin | 0.7 | 1 |
| campylobacter_jejuni | ciprofloxacin | 0.8 | 4 |
| campylobacter_jejuni | levofloxacin | 0.75 | 1 |
| campylobacter_jejuni | moxifloxacin | 0.7 | 1 |
| campylobacter_jejuni | ofloxacin | 0.75 | 1 |
| campylobacter_jejuni | tetracycline | 0.75 | 0.25 |
| campylobacter_jejuni | doxycycline | 0.8 | 0.25 |
| campylobacter_jejuni | minocycline | 0.8 | 0.25 |
| campylobacter_jejuni | tigecycline | 0.7 | 1 |
| campylobacter_jejuni | vancomycin | 0.1 | 1 |
| campylobacter_jejuni | teicoplanin | 0.1 | 1 |
| campylobacter_jejuni | dalbavancin | 0.1 | 0.005 |
| campylobacter_jejuni | linezolid | 0.1 | 0.005 |
| campylobacter_jejuni | tedizolid | 0.1 | 0.005 |
| campylobacter_jejuni | daptomycin | 0.1 | 1 |
| campylobacter_jejuni | quinu_dalfo | 0.1 | 0.005 |
| campylobacter_jejuni | trim_sulf | 0.1 | 0.04 |
| campylobacter_jejuni | chloramphenicol | 0.7 | 1 |
| campylobacter_jejuni | nitrofurantoin | 0.1 | 1 |
| campylobacter_jejuni | fosfomycin | 0.1 | 1 |
| campylobacter_jejuni | retapamulin | 0.05 | 1 |
| campylobacter_jejuni | fusidic_a | 0.05 | 1 |
| campylobacter_jejuni | metronidazole | 0.05 | 1 |
| campylobacter_jejuni | fidaxomicin | 0.1 | 1 |
| campylobacter_jejuni | furazolidone | 0.1 | 1 |
| campylobacter_jejuni | rifampicin | 0.1 | 1 |
| campylobacter_jejuni | amoxicillin_clavulanate | 0.1 | 1 |
| campylobacter_jejuni | piperacillin_tazobactam | 0.1 | 1 |
| campylobacter_jejuni | ampicillin_sulbactam | 0.1 | 1 |
| campylobacter_jejuni | ticarcillin_clavulanate | 0.1 | 1 |
| campylobacter_jejuni | ceftazidime_avibactam | 0.1 | 0.005 |
| campylobacter_jejuni | meropenem_vaborbactam | 0.1 | 0.005 |
| campylobacter_jejuni | colistin | 0.05 | 0.005 |
| campylobacter_jejuni | flucloxacillin | 0.01 | 1 |
| campylobacter_jejuni | aztreonam_avibactam | 0.9 | 0.003 |
| campylobacter_jejuni | cefixime | 0.75 | 0.2 |
| enterobacter_cloacae | sulfanilamide | 0.5 | 0.02 |
| enterobacter_cloacae | penicillin_g | 0.1 | 1 |
| enterobacter_cloacae | ampicillin | 0.5 | 1 |
| enterobacter_cloacae | amoxicillin | 0.5 | 1 |
| enterobacter_cloacae | piperacillin | 0.75 | 1 |
| enterobacter_cloacae | ticarcillin | 0.7 | 1 |
| enterobacter_cloacae | cephalexin | 0.5 | 0.3 |
| enterobacter_cloacae | cefazolin | 0.5 | 0.3 |
| enterobacter_cloacae | cefuroxime | 0.6 | 0.3 |
| enterobacter_cloacae | ceftriaxone | 0.4 | 0.2 |
| enterobacter_cloacae | ceftazidime | 0.8 | 0.2 |
| enterobacter_cloacae | cefepime | 0.85 | 2.5 |
| enterobacter_cloacae | ceftaroline | 0.4 | 0.002 |
| enterobacter_cloacae | ceftolozane_tazobactam | 0.8 | 1 |
| enterobacter_cloacae | cefiderocol | 0.8 | 1 |
| enterobacter_cloacae | meropenem | 0.95 | 5 |
| enterobacter_cloacae | imipenem_c | 0.95 | 3 |
| enterobacter_cloacae | ertapenem | 0.9 | 3 |
| enterobacter_cloacae | aztreonam | 0.8 | 0.003 |
| enterobacter_cloacae | erythromycin | 0.1 | 1 |
| enterobacter_cloacae | azithromycin | 0.1 | 1 |
| enterobacter_cloacae | clarithromycin | 0.1 | 1 |
| enterobacter_cloacae | clindamycin | 0.1 | 1 |
| enterobacter_cloacae | gentamicin | 0.85 | 10 |
| enterobacter_cloacae | tobramycin | 0.8 | 1 |
| enterobacter_cloacae | amikacin | 0.9 | 1 |
| enterobacter_cloacae | ciprofloxacin | 0.9 | 1 |
| enterobacter_cloacae | levofloxacin | 0.85 | 1 |
| enterobacter_cloacae | moxifloxacin | 0.7 | 1 |
| enterobacter_cloacae | ofloxacin | 0.8 | 1 |
| enterobacter_cloacae | tetracycline | 0.8 | 0.25 |
| enterobacter_cloacae | doxycycline | 0.85 | 0.25 |
| enterobacter_cloacae | minocycline | 0.85 | 0.25 |
| enterobacter_cloacae | tigecycline | 0.1 | 1 |
| enterobacter_cloacae | vancomycin | 0.1 | 1 |
| enterobacter_cloacae | teicoplanin | 0.1 | 1 |
| enterobacter_cloacae | dalbavancin | 0.1 | 0.005 |
| enterobacter_cloacae | linezolid | 0.1 | 0.005 |
| enterobacter_cloacae | tedizolid | 0.1 | 0.005 |
| enterobacter_cloacae | daptomycin | 0.1 | 1 |
| enterobacter_cloacae | quinu_dalfo | 0.1 | 0.005 |
| enterobacter_cloacae | trim_sulf | 0.85 | 0.04 |
| enterobacter_cloacae | chloramphenicol | 0.8 | 1 |
| enterobacter_cloacae | nitrofurantoin | 0.7 | 1 |
| enterobacter_cloacae | fosfomycin | 0.1 | 1 |
| enterobacter_cloacae | retapamulin | 0.05 | 1 |
| enterobacter_cloacae | fusidic_a | 0.05 | 1 |
| enterobacter_cloacae | metronidazole | 0.05 | 1 |
| enterobacter_cloacae | fidaxomicin | 0.1 | 1 |
| enterobacter_cloacae | furazolidone | 0.1 | 1 |
| enterobacter_cloacae | rifampicin | 0.6 | 1 |
| enterobacter_cloacae | amoxicillin_clavulanate | 0.7 | 6 |
| enterobacter_cloacae | piperacillin_tazobactam | 0.85 | 8 |
| enterobacter_cloacae | ampicillin_sulbactam | 0.7 | 1 |
| enterobacter_cloacae | ticarcillin_clavulanate | 0.8 | 1 |
| enterobacter_cloacae | ceftazidime_avibactam | 0.9 | 0.005 |
| enterobacter_cloacae | meropenem_vaborbactam | 0.95 | 0.005 |
| enterobacter_cloacae | colistin | 0.7 | 0.005 |
| enterobacter_cloacae | flucloxacillin | 0.01 | 1 |
| enterobacter_cloacae | aztreonam_avibactam | 1 | 0.003 |
| enterobacter_cloacae | cefixime | 0.8 | 0.2 |
| yersinia_enterocolitica | sulfanilamide | 0.5 | 0.02 |
| yersinia_enterocolitica | penicillin_g | 0.1 | 1 |
| yersinia_enterocolitica | ampicillin | 0.7 | 1 |
| yersinia_enterocolitica | amoxicillin | 0.7 | 1 |
| yersinia_enterocolitica | piperacillin | 0.75 | 1 |
| yersinia_enterocolitica | ticarcillin | 0.7 | 1 |
| yersinia_enterocolitica | cephalexin | 0.6 | 0.3 |
| yersinia_enterocolitica | cefazolin | 0.65 | 0.3 |
| yersinia_enterocolitica | cefuroxime | 0.7 | 0.3 |
| yersinia_enterocolitica | ceftriaxone | 0.9 | 0.2 |
| yersinia_enterocolitica | ceftazidime | 0.85 | 0.2 |
| yersinia_enterocolitica | cefepime | 0.85 | 0.35 |
| yersinia_enterocolitica | ceftaroline | 0.6 | 0.002 |
| yersinia_enterocolitica | ceftolozane_tazobactam | 0.75 | 1 |
| yersinia_enterocolitica | cefiderocol | 0.75 | 1 |
| yersinia_enterocolitica | meropenem | 0.95 | 0.005 |
| yersinia_enterocolitica | imipenem_c | 0.95 | 0.005 |
| yersinia_enterocolitica | ertapenem | 0.95 | 0.005 |
| yersinia_enterocolitica | aztreonam | 0.85 | 0.003 |
| yersinia_enterocolitica | erythromycin | 0.1 | 1 |
| yersinia_enterocolitica | azithromycin | 0.1 | 1 |
| yersinia_enterocolitica | clarithromycin | 0.1 | 1 |
| yersinia_enterocolitica | clindamycin | 0.1 | 1 |
| yersinia_enterocolitica | gentamicin | 0.85 | 1 |
| yersinia_enterocolitica | tobramycin | 0.8 | 1 |
| yersinia_enterocolitica | amikacin | 0.9 | 1 |
| yersinia_enterocolitica | ciprofloxacin | 0.9 | 1 |
| yersinia_enterocolitica | levofloxacin | 0.85 | 1 |
| yersinia_enterocolitica | moxifloxacin | 0.7 | 1 |
| yersinia_enterocolitica | ofloxacin | 0.8 | 1 |
| yersinia_enterocolitica | tetracycline | 0.8 | 0.25 |
| yersinia_enterocolitica | doxycycline | 0.85 | 2 |
| yersinia_enterocolitica | minocycline | 0.85 | 0.25 |
| yersinia_enterocolitica | tigecycline | 0.7 | 1 |
| yersinia_enterocolitica | vancomycin | 0.1 | 1 |
| yersinia_enterocolitica | teicoplanin | 0.1 | 1 |
| yersinia_enterocolitica | dalbavancin | 0.1 | 0.005 |
| yersinia_enterocolitica | linezolid | 0.1 | 0.005 |
| yersinia_enterocolitica | tedizolid | 0.1 | 0.005 |
| yersinia_enterocolitica | daptomycin | 0.1 | 1 |
| yersinia_enterocolitica | quinu_dalfo | 0.1 | 0.005 |
| yersinia_enterocolitica | trim_sulf | 0.95 | 0.04 |
| yersinia_enterocolitica | chloramphenicol | 0.85 | 1 |
| yersinia_enterocolitica | nitrofurantoin | 0.1 | 1 |
| yersinia_enterocolitica | fosfomycin | 0.1 | 1 |
| yersinia_enterocolitica | retapamulin | 0.05 | 1 |
| yersinia_enterocolitica | fusidic_a | 0.05 | 1 |
| yersinia_enterocolitica | metronidazole | 0.05 | 1 |
| yersinia_enterocolitica | fidaxomicin | 0.1 | 1 |
| yersinia_enterocolitica | furazolidone | 0.1 | 1 |
| yersinia_enterocolitica | rifampicin | 0.7 | 1 |
| yersinia_enterocolitica | amoxicillin_clavulanate | 0.85 | 1 |
| yersinia_enterocolitica | piperacillin_tazobactam | 0.85 | 1 |
| yersinia_enterocolitica | ampicillin_sulbactam | 0.8 | 1 |
| yersinia_enterocolitica | ticarcillin_clavulanate | 0.8 | 1 |
| yersinia_enterocolitica | ceftazidime_avibactam | 0.95 | 0.005 |
| yersinia_enterocolitica | meropenem_vaborbactam | 0.95 | 0.005 |
| yersinia_enterocolitica | colistin | 0.7 | 0.005 |
| yersinia_enterocolitica | flucloxacillin | 0.01 | 1 |
| yersinia_enterocolitica | aztreonam_avibactam | 0.9 | 0.003 |
| yersinia_enterocolitica | cefixime | 0.75 | 0.2 |
| moraxella_catarrhalis | sulfanilamide | 0.1 | 0.02 |
| moraxella_catarrhalis | penicillin_g | 0.9 | 1 |
| moraxella_catarrhalis | ampicillin | 0.9 | 1 |
| moraxella_catarrhalis | amoxicillin | 0.9 | 6 |
| moraxella_catarrhalis | piperacillin | 0.8 | 1 |
| moraxella_catarrhalis | ticarcillin | 0.8 | 1 |
| moraxella_catarrhalis | cephalexin | 0.8 | 0.3 |
| moraxella_catarrhalis | cefazolin | 0.85 | 0.3 |
| moraxella_catarrhalis | cefuroxime | 0.9 | 0.3 |
| moraxella_catarrhalis | ceftriaxone | 0.95 | 0.2 |
| moraxella_catarrhalis | ceftazidime | 0.9 | 0.2 |
| moraxella_catarrhalis | cefepime | 0.9 | 0.35 |
| moraxella_catarrhalis | ceftaroline | 0.8 | 0.002 |
| moraxella_catarrhalis | ceftolozane_tazobactam | 0.8 | 1 |
| moraxella_catarrhalis | cefiderocol | 0.8 | 1 |
| moraxella_catarrhalis | meropenem | 0.95 | 0.005 |
| moraxella_catarrhalis | imipenem_c | 0.95 | 0.005 |
| moraxella_catarrhalis | ertapenem | 0.95 | 0.005 |
| moraxella_catarrhalis | aztreonam | 0.9 | 0.003 |
| moraxella_catarrhalis | erythromycin | 0.8 | 1 |
| moraxella_catarrhalis | azithromycin | 0.85 | 5 |
| moraxella_catarrhalis | clarithromycin | 0.8 | 5 |
| moraxella_catarrhalis | clindamycin | 0.1 | 1 |
| moraxella_catarrhalis | gentamicin | 0.1 | 1 |
| moraxella_catarrhalis | tobramycin | 0.1 | 1 |
| moraxella_catarrhalis | amikacin | 0.1 | 1 |
| moraxella_catarrhalis | ciprofloxacin | 0.9 | 1 |
| moraxella_catarrhalis | levofloxacin | 0.85 | 4 |
| moraxella_catarrhalis | moxifloxacin | 0.8 | 1 |
| moraxella_catarrhalis | ofloxacin | 0.85 | 1 |
| moraxella_catarrhalis | tetracycline | 0.8 | 0.25 |
| moraxella_catarrhalis | doxycycline | 0.8 | 0.25 |
| moraxella_catarrhalis | minocycline | 0.85 | 0.25 |
| moraxella_catarrhalis | tigecycline | 0.1 | 1 |
| moraxella_catarrhalis | vancomycin | 0.1 | 1 |
| moraxella_catarrhalis | teicoplanin | 0.1 | 1 |
| moraxella_catarrhalis | dalbavancin | 0.1 | 0.005 |
| moraxella_catarrhalis | linezolid | 0.1 | 0.005 |
| moraxella_catarrhalis | tedizolid | 0.1 | 0.005 |
| moraxella_catarrhalis | daptomycin | 0.1 | 1 |
| moraxella_catarrhalis | quinu_dalfo | 0.1 | 0.005 |
| moraxella_catarrhalis | trim_sulf | 0.95 | 0.04 |
| moraxella_catarrhalis | chloramphenicol | 0.85 | 1 |
| moraxella_catarrhalis | nitrofurantoin | 0.1 | 1 |
| moraxella_catarrhalis | fosfomycin | 0.1 | 1 |
| moraxella_catarrhalis | retapamulin | 0.05 | 1 |
| moraxella_catarrhalis | fusidic_a | 0.05 | 1 |
| moraxella_catarrhalis | metronidazole | 0.05 | 1 |
| moraxella_catarrhalis | fidaxomicin | 0.1 | 1 |
| moraxella_catarrhalis | furazolidone | 0.1 | 1 |
| moraxella_catarrhalis | rifampicin | 0.7 | 1 |
| moraxella_catarrhalis | amoxicillin_clavulanate | 0.95 | 12 |
| moraxella_catarrhalis | piperacillin_tazobactam | 0.85 | 1 |
| moraxella_catarrhalis | ampicillin_sulbactam | 0.95 | 1 |
| moraxella_catarrhalis | ticarcillin_clavulanate | 0.85 | 1 |
| moraxella_catarrhalis | ceftazidime_avibactam | 0.95 | 0.005 |
| moraxella_catarrhalis | meropenem_vaborbactam | 0.95 | 0.005 |
| moraxella_catarrhalis | colistin | 0.05 | 0.005 |
| moraxella_catarrhalis | flucloxacillin | 0.01 | 1 |
| moraxella_catarrhalis | aztreonam_avibactam | 0.8 | 0.003 |
| moraxella_catarrhalis | cefixime | 0.8 | 0.2 |
| treponema_pallidum | sulfanilamide | 0.1 | 0.02 |
| treponema_pallidum | penicillin_g | 1 | 6 |
| treponema_pallidum | ampicillin | 0.95 | 1 |
| treponema_pallidum | amoxicillin | 0.95 | 1 |
| treponema_pallidum | piperacillin | 0.9 | 1 |
| treponema_pallidum | ticarcillin | 0.9 | 1 |
| treponema_pallidum | cephalexin | 0.9 | 0.3 |
| treponema_pallidum | cefazolin | 0.9 | 0.3 |
| treponema_pallidum | cefuroxime | 0.95 | 0.3 |
| treponema_pallidum | ceftriaxone | 0.95 | 0.2 |
| treponema_pallidum | ceftazidime | 0.9 | 0.2 |
| treponema_pallidum | cefepime | 0.9 | 0.35 |
| treponema_pallidum | ceftaroline | 0.9 | 0.002 |
| treponema_pallidum | ceftolozane_tazobactam | 0.1 | 1 |
| treponema_pallidum | cefiderocol | 0.1 | 1 |
| treponema_pallidum | meropenem | 0.95 | 0.005 |
| treponema_pallidum | imipenem_c | 0.95 | 0.005 |
| treponema_pallidum | ertapenem | 0.95 | 0.005 |
| treponema_pallidum | aztreonam | 0.9 | 0.003 |
| treponema_pallidum | erythromycin | 0.8 | 1 |
| treponema_pallidum | azithromycin | 0.85 | 1 |
| treponema_pallidum | clarithromycin | 0.8 | 1 |
| treponema_pallidum | clindamycin | 0.1 | 1 |
| treponema_pallidum | gentamicin | 0.1 | 1 |
| treponema_pallidum | tobramycin | 0.1 | 1 |
| treponema_pallidum | amikacin | 0.1 | 1 |
| treponema_pallidum | ciprofloxacin | 0.7 | 1 |
| treponema_pallidum | levofloxacin | 0.75 | 1 |
| treponema_pallidum | moxifloxacin | 0.75 | 1 |
| treponema_pallidum | ofloxacin | 0.7 | 1 |
| treponema_pallidum | tetracycline | 0.8 | 0.25 |
| treponema_pallidum | doxycycline | 0.8 | 2 |
| treponema_pallidum | minocycline | 0.85 | 0.25 |
| treponema_pallidum | tigecycline | 0.1 | 1 |
| treponema_pallidum | vancomycin | 0.1 | 1 |
| treponema_pallidum | teicoplanin | 0.1 | 1 |
| treponema_pallidum | dalbavancin | 0.1 | 0.005 |
| treponema_pallidum | linezolid | 0.1 | 0.005 |
| treponema_pallidum | tedizolid | 0.1 | 0.005 |
| treponema_pallidum | daptomycin | 0.1 | 1 |
| treponema_pallidum | quinu_dalfo | 0.1 | 0.005 |
| treponema_pallidum | trim_sulf | 0.1 | 0.04 |
| treponema_pallidum | chloramphenicol | 0.8 | 1 |
| treponema_pallidum | nitrofurantoin | 0.1 | 1 |
| treponema_pallidum | fosfomycin | 0.1 | 1 |
| treponema_pallidum | retapamulin | 0.05 | 1 |
| treponema_pallidum | fusidic_a | 0.05 | 1 |
| treponema_pallidum | metronidazole | 0.05 | 1 |
| treponema_pallidum | fidaxomicin | 0.1 | 1 |
| treponema_pallidum | furazolidone | 0.1 | 1 |
| treponema_pallidum | rifampicin | 0.1 | 1 |
| treponema_pallidum | amoxicillin_clavulanate | 0.95 | 1 |
| treponema_pallidum | piperacillin_tazobactam | 0.9 | 1 |
| treponema_pallidum | ampicillin_sulbactam | 0.95 | 1 |
| treponema_pallidum | ticarcillin_clavulanate | 0.9 | 1 |
| treponema_pallidum | ceftazidime_avibactam | 0.95 | 0.005 |
| treponema_pallidum | meropenem_vaborbactam | 0.95 | 0.005 |
| treponema_pallidum | colistin | 0.05 | 0.005 |
| treponema_pallidum | flucloxacillin | 0.01 | 1 |
| treponema_pallidum | aztreonam_avibactam | 0.9 | 0.003 |
| treponema_pallidum | cefixime | 0.1 | 0.2 |
| bordetella_pertussis | sulfanilamide | 0.1 | 0.02 |
| bordetella_pertussis | penicillin_g | 0.1 | 1 |
| bordetella_pertussis | ampicillin | 0.1 | 1 |
| bordetella_pertussis | amoxicillin | 0.1 | 1 |
| bordetella_pertussis | piperacillin | 0.1 | 1 |
| bordetella_pertussis | ticarcillin | 0.1 | 1 |
| bordetella_pertussis | cephalexin | 0.1 | 0.3 |
| bordetella_pertussis | cefazolin | 0.1 | 0.3 |
| bordetella_pertussis | cefuroxime | 0.1 | 0.3 |
| bordetella_pertussis | ceftriaxone | 0.1 | 0.2 |
| bordetella_pertussis | ceftazidime | 0.1 | 0.2 |
| bordetella_pertussis | cefepime | 0.1 | 0.35 |
| bordetella_pertussis | ceftaroline | 0.1 | 0.002 |
| bordetella_pertussis | ceftolozane_tazobactam | 0.8 | 1 |
| bordetella_pertussis | cefiderocol | 0.8 | 1 |
| bordetella_pertussis | meropenem | 0.1 | 0.005 |
| bordetella_pertussis | imipenem_c | 0.1 | 0.005 |
| bordetella_pertussis | ertapenem | 0.1 | 0.005 |
| bordetella_pertussis | aztreonam | 0.1 | 0.003 |
| bordetella_pertussis | erythromycin | 0.9 | 7 |
| bordetella_pertussis | azithromycin | 0.95 | 8 |
| bordetella_pertussis | clarithromycin | 0.9 | 7 |
| bordetella_pertussis | clindamycin | 0.1 | 1 |
| bordetella_pertussis | gentamicin | 0.7 | 1 |
| bordetella_pertussis | tobramycin | 0.7 | 1 |
| bordetella_pertussis | amikacin | 0.7 | 1 |
| bordetella_pertussis | ciprofloxacin | 0.7 | 1 |
| bordetella_pertussis | levofloxacin | 0.75 | 1 |
| bordetella_pertussis | moxifloxacin | 0.75 | 1 |
| bordetella_pertussis | ofloxacin | 0.7 | 1 |
| bordetella_pertussis | tetracycline | 0.7 | 0.25 |
| bordetella_pertussis | doxycycline | 0.75 | 1.5 |
| bordetella_pertussis | minocycline | 0.75 | 0.25 |
| bordetella_pertussis | tigecycline | 0.1 | 1 |
| bordetella_pertussis | vancomycin | 0.1 | 1 |
| bordetella_pertussis | teicoplanin | 0.1 | 1 |
| bordetella_pertussis | dalbavancin | 0.1 | 0.005 |
| bordetella_pertussis | linezolid | 0.1 | 0.005 |
| bordetella_pertussis | tedizolid | 0.1 | 0.005 |
| bordetella_pertussis | daptomycin | 0.1 | 1 |
| bordetella_pertussis | quinu_dalfo | 0.1 | 0.005 |
| bordetella_pertussis | trim_sulf | 0.7 | 0.03 |
| bordetella_pertussis | chloramphenicol | 0.8 | 1 |
| bordetella_pertussis | nitrofurantoin | 0.1 | 1 |
| bordetella_pertussis | fosfomycin | 0.1 | 1 |
| bordetella_pertussis | retapamulin | 0.05 | 1 |
| bordetella_pertussis | fusidic_a | 0.05 | 1 |
| bordetella_pertussis | metronidazole | 0.05 | 1 |
| bordetella_pertussis | fidaxomicin | 0.1 | 1 |
| bordetella_pertussis | furazolidone | 0.1 | 1 |
| bordetella_pertussis | rifampicin | 0.7 | 1 |
| bordetella_pertussis | amoxicillin_clavulanate | 0.1 | 1 |
| bordetella_pertussis | piperacillin_tazobactam | 0.1 | 1 |
| bordetella_pertussis | ampicillin_sulbactam | 0.1 | 1 |
| bordetella_pertussis | ticarcillin_clavulanate | 0.1 | 1 |
| bordetella_pertussis | ceftazidime_avibactam | 0.1 | 0.005 |
| bordetella_pertussis | meropenem_vaborbactam | 0.1 | 0.005 |
| bordetella_pertussis | colistin | 0.05 | 0.005 |
| bordetella_pertussis | flucloxacillin | 0.01 | 1 |
| bordetella_pertussis | aztreonam_avibactam | 0.8 | 0.003 |
| bordetella_pertussis | cefixime | 0.8 | 0.2 |
| helicobacter_pylori | sulfanilamide | 0.1 | 0.02 |
| helicobacter_pylori | penicillin_g | 0.1 | 1 |
| helicobacter_pylori | ampicillin | 0.7 | 1 |
| helicobacter_pylori | amoxicillin | 0.85 | 12 |
| helicobacter_pylori | piperacillin | 0.1 | 1 |
| helicobacter_pylori | ticarcillin | 0.1 | 1 |
| helicobacter_pylori | cephalexin | 0.1 | 0.3 |
| helicobacter_pylori | cefazolin | 0.1 | 0.3 |
| helicobacter_pylori | cefuroxime | 0.1 | 0.3 |
| helicobacter_pylori | ceftriaxone | 0.1 | 0.2 |
| helicobacter_pylori | ceftazidime | 0.1 | 0.2 |
| helicobacter_pylori | cefepime | 0.1 | 0.35 |
| helicobacter_pylori | ceftaroline | 0.1 | 0.002 |
| helicobacter_pylori | ceftolozane_tazobactam | 0.05 | 1 |
| helicobacter_pylori | cefiderocol | 0.05 | 1 |
| helicobacter_pylori | meropenem | 0.1 | 0.005 |
| helicobacter_pylori | imipenem_c | 0.1 | 0.005 |
| helicobacter_pylori | ertapenem | 0.1 | 0.005 |
| helicobacter_pylori | aztreonam | 0.1 | 0.003 |
| helicobacter_pylori | erythromycin | 0.8 | 1 |
| helicobacter_pylori | azithromycin | 0.85 | 1 |
| helicobacter_pylori | clarithromycin | 0.8 | 5 |
| helicobacter_pylori | clindamycin | 0.1 | 1 |
| helicobacter_pylori | gentamicin | 0.1 | 1 |
| helicobacter_pylori | tobramycin | 0.1 | 1 |
| helicobacter_pylori | amikacin | 0.1 | 1 |
| helicobacter_pylori | ciprofloxacin | 0.7 | 1 |
| helicobacter_pylori | levofloxacin | 0.7 | 5 |
| helicobacter_pylori | moxifloxacin | 0.75 | 1 |
| helicobacter_pylori | ofloxacin | 0.7 | 1 |
| helicobacter_pylori | tetracycline | 0.8 | 2 |
| helicobacter_pylori | doxycycline | 0.8 | 0.25 |
| helicobacter_pylori | minocycline | 0.85 | 0.25 |
| helicobacter_pylori | tigecycline | 0.1 | 1 |
| helicobacter_pylori | vancomycin | 0.1 | 1 |
| helicobacter_pylori | teicoplanin | 0.1 | 1 |
| helicobacter_pylori | dalbavancin | 0.1 | 0.005 |
| helicobacter_pylori | linezolid | 0.1 | 0.005 |
| helicobacter_pylori | tedizolid | 0.1 | 0.005 |
| helicobacter_pylori | daptomycin | 0.1 | 1 |
| helicobacter_pylori | quinu_dalfo | 0.1 | 0.005 |
| helicobacter_pylori | trim_sulf | 0.1 | 0.04 |
| helicobacter_pylori | chloramphenicol | 0.7 | 1 |
| helicobacter_pylori | nitrofurantoin | 0.1 | 1 |
| helicobacter_pylori | fosfomycin | 0.1 | 1 |
| helicobacter_pylori | retapamulin | 0.05 | 1 |
| helicobacter_pylori | fusidic_a | 0.05 | 1 |
| helicobacter_pylori | metronidazole | 0.8 | 10 |
| helicobacter_pylori | fidaxomicin | 0.1 | 1 |
| helicobacter_pylori | furazolidone | 0.1 | 6 |
| helicobacter_pylori | rifampicin | 0.1 | 1 |
| helicobacter_pylori | amoxicillin_clavulanate | 0.85 | 1 |
| helicobacter_pylori | piperacillin_tazobactam | 0.1 | 1 |
| helicobacter_pylori | ampicillin_sulbactam | 0.7 | 1 |
| helicobacter_pylori | ticarcillin_clavulanate | 0.1 | 1 |
| helicobacter_pylori | ceftazidime_avibactam | 0.1 | 0.005 |
| helicobacter_pylori | meropenem_vaborbactam | 0.1 | 0.005 |
| helicobacter_pylori | colistin | 0.05 | 0.005 |
| helicobacter_pylori | flucloxacillin | 0.01 | 1 |
| helicobacter_pylori | aztreonam_avibactam | 0.01 | 0.003 |
| helicobacter_pylori | cefixime | 0.05 | 0.2 |
| mdr_mycobacterium_tuberculosis | sulfanilamide | 0.1 | 0.02 |
| mdr_mycobacterium_tuberculosis | penicillin_g | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ampicillin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | amoxicillin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | piperacillin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ticarcillin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | cephalexin | 0.05 | 0.3 |
| mdr_mycobacterium_tuberculosis | cefazolin | 0.05 | 0.3 |
| mdr_mycobacterium_tuberculosis | cefuroxime | 0.05 | 0.3 |
| mdr_mycobacterium_tuberculosis | ceftriaxone | 0.05 | 0.2 |
| mdr_mycobacterium_tuberculosis | ceftazidime | 0.05 | 0.2 |
| mdr_mycobacterium_tuberculosis | cefepime | 0.05 | 0.35 |
| mdr_mycobacterium_tuberculosis | ceftaroline | 0.05 | 0.002 |
| mdr_mycobacterium_tuberculosis | ceftolozane_tazobactam | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | cefiderocol | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | meropenem | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | imipenem_c | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | ertapenem | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | aztreonam | 0.2 | 0.003 |
| mdr_mycobacterium_tuberculosis | erythromycin | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | azithromycin | 0.25 | 1 |
| mdr_mycobacterium_tuberculosis | clarithromycin | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | clindamycin | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | gentamicin | 0.25 | 1 |
| mdr_mycobacterium_tuberculosis | tobramycin | 0.25 | 1 |
| mdr_mycobacterium_tuberculosis | amikacin | 0.3 | 1 |
| mdr_mycobacterium_tuberculosis | ciprofloxacin | 0.4 | 1 |
| mdr_mycobacterium_tuberculosis | levofloxacin | 0.45 | 1 |
| mdr_mycobacterium_tuberculosis | moxifloxacin | 0.45 | 1 |
| mdr_mycobacterium_tuberculosis | ofloxacin | 0.4 | 1 |
| mdr_mycobacterium_tuberculosis | tetracycline | 0.3 | 0.25 |
| mdr_mycobacterium_tuberculosis | doxycycline | 0.35 | 0.25 |
| mdr_mycobacterium_tuberculosis | minocycline | 0.35 | 0.25 |
| mdr_mycobacterium_tuberculosis | tigecycline | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | vancomycin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | teicoplanin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | dalbavancin | 0.1 | 0.005 |
| mdr_mycobacterium_tuberculosis | linezolid | 0.1 | 0.005 |
| mdr_mycobacterium_tuberculosis | tedizolid | 0.1 | 0.005 |
| mdr_mycobacterium_tuberculosis | daptomycin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | quinu_dalfo | 0.1 | 0.005 |
| mdr_mycobacterium_tuberculosis | trim_sulf | 0.2 | 0.04 |
| mdr_mycobacterium_tuberculosis | chloramphenicol | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | nitrofurantoin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | fosfomycin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | retapamulin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | fusidic_a | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | metronidazole | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | fidaxomicin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | furazolidone | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | rifampicin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | amoxicillin_clavulanate | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | piperacillin_tazobactam | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ampicillin_sulbactam | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ticarcillin_clavulanate | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ceftazidime_avibactam | 0.05 | 0.005 |
| mdr_mycobacterium_tuberculosis | meropenem_vaborbactam | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | colistin | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | flucloxacillin | 0.01 | 1 |
| mdr_mycobacterium_tuberculosis | aztreonam_avibactam | 0.9 | 0.003 |
| mdr_mycobacterium_tuberculosis | cefixime | 0.1 | 0.2 |
| mycoplasma_pneumoniae | sulfanilamide | 0.05 | 0.02 |
| mycoplasma_pneumoniae | penicillin_g | 0.05 | 0.001 |
| mycoplasma_pneumoniae | ampicillin | 0.05 | 0.001 |
| mycoplasma_pneumoniae | amoxicillin | 0.05 | 0.001 |
| mycoplasma_pneumoniae | piperacillin | 0.05 | 1 |
| mycoplasma_pneumoniae | ticarcillin | 0.05 | 1 |
| mycoplasma_pneumoniae | cephalexin | 0.05 | 0.3 |
| mycoplasma_pneumoniae | cefazolin | 0.05 | 0.3 |
| mycoplasma_pneumoniae | cefuroxime | 0.05 | 0.3 |
| mycoplasma_pneumoniae | ceftriaxone | 0.05 | 0.2 |
| mycoplasma_pneumoniae | ceftazidime | 0.05 | 0.2 |
| mycoplasma_pneumoniae | cefepime | 0.05 | 0.35 |
| mycoplasma_pneumoniae | ceftaroline | 0.05 | 0.002 |
| mycoplasma_pneumoniae | ceftolozane_tazobactam | 0.01 | 1 |
| mycoplasma_pneumoniae | cefiderocol | 0.01 | 1 |
| mycoplasma_pneumoniae | meropenem | 0.05 | 0.001 |
| mycoplasma_pneumoniae | imipenem_c | 0.05 | 0.005 |
| mycoplasma_pneumoniae | ertapenem | 0.05 | 0.001 |
| mycoplasma_pneumoniae | aztreonam | 0.05 | 0.003 |
| mycoplasma_pneumoniae | erythromycin | 0.8 | 1 |
| mycoplasma_pneumoniae | azithromycin | 0.85 | 8 |
| mycoplasma_pneumoniae | clarithromycin | 0.8 | 7 |
| mycoplasma_pneumoniae | clindamycin | 0.05 | 1 |
| mycoplasma_pneumoniae | gentamicin | 0.05 | 1 |
| mycoplasma_pneumoniae | tobramycin | 0.05 | 1 |
| mycoplasma_pneumoniae | amikacin | 0.05 | 1 |
| mycoplasma_pneumoniae | ciprofloxacin | 0.7 | 1 |
| mycoplasma_pneumoniae | levofloxacin | 0.75 | 4 |
| mycoplasma_pneumoniae | moxifloxacin | 0.8 | 4 |
| mycoplasma_pneumoniae | ofloxacin | 0.6 | 1 |
| mycoplasma_pneumoniae | tetracycline | 0.7 | 0.25 |
| mycoplasma_pneumoniae | doxycycline | 0.75 | 1.5 |
| mycoplasma_pneumoniae | minocycline | 0.8 | 0.25 |
| mycoplasma_pneumoniae | tigecycline | 0.85 | 1 |
| mycoplasma_pneumoniae | vancomycin | 0.05 | 1 |
| mycoplasma_pneumoniae | teicoplanin | 0.05 | 1 |
| mycoplasma_pneumoniae | dalbavancin | 0.05 | 0.005 |
| mycoplasma_pneumoniae | linezolid | 0.05 | 0.005 |
| mycoplasma_pneumoniae | tedizolid | 0.05 | 0.005 |
| mycoplasma_pneumoniae | daptomycin | 0.1 | 1 |
| mycoplasma_pneumoniae | quinu_dalfo | 0.05 | 0.005 |
| mycoplasma_pneumoniae | trim_sulf | 0.05 | 0.04 |
| mycoplasma_pneumoniae | chloramphenicol | 0.05 | 1 |
| mycoplasma_pneumoniae | nitrofurantoin | 0.05 | 1 |
| mycoplasma_pneumoniae | fosfomycin | 0.1 | 1 |
| mycoplasma_pneumoniae | retapamulin | 0.05 | 1 |
| mycoplasma_pneumoniae | fusidic_a | 0.05 | 1 |
| mycoplasma_pneumoniae | metronidazole | 0.05 | 1 |
| mycoplasma_pneumoniae | fidaxomicin | 0.1 | 1 |
| mycoplasma_pneumoniae | furazolidone | 0.05 | 1 |
| mycoplasma_pneumoniae | rifampicin | 0.05 | 1 |
| mycoplasma_pneumoniae | amoxicillin_clavulanate | 0.05 | 1 |
| mycoplasma_pneumoniae | piperacillin_tazobactam | 0.05 | 1 |
| mycoplasma_pneumoniae | ampicillin_sulbactam | 0.05 | 1 |
| mycoplasma_pneumoniae | ticarcillin_clavulanate | 0.05 | 1 |
| mycoplasma_pneumoniae | ceftazidime_avibactam | 0.05 | 0.005 |
| mycoplasma_pneumoniae | meropenem_vaborbactam | 0.05 | 0.005 |
| mycoplasma_pneumoniae | colistin | 0.05 | 0.005 |
| mycoplasma_pneumoniae | flucloxacillin | 0.01 | 1 |
| mycoplasma_pneumoniae | aztreonam_avibactam | 0.01 | 0.003 |
| mycoplasma_pneumoniae | cefixime | 0.01 | 0.2 |
| legionella_pneumophila | sulfanilamide | 0.05 | 0.02 |
| legionella_pneumophila | penicillin_g | 0.05 | 0.001 |
| legionella_pneumophila | ampicillin | 0.05 | 0.001 |
| legionella_pneumophila | amoxicillin | 0.05 | 0.001 |
| legionella_pneumophila | piperacillin | 0.05 | 1 |
| legionella_pneumophila | ticarcillin | 0.05 | 1 |
| legionella_pneumophila | cephalexin | 0.05 | 0.3 |
| legionella_pneumophila | cefazolin | 0.05 | 0.3 |
| legionella_pneumophila | cefuroxime | 0.05 | 0.3 |
| legionella_pneumophila | ceftriaxone | 0.05 | 0.2 |
| legionella_pneumophila | ceftazidime | 0.05 | 0.2 |
| legionella_pneumophila | cefepime | 0.05 | 0.35 |
| legionella_pneumophila | ceftaroline | 0.05 | 0.002 |
| legionella_pneumophila | ceftolozane_tazobactam | 0.8 | 1 |
| legionella_pneumophila | cefiderocol | 0.8 | 1 |
| legionella_pneumophila | meropenem | 0.05 | 0.001 |
| legionella_pneumophila | imipenem_c | 0.05 | 0.005 |
| legionella_pneumophila | ertapenem | 0.05 | 0.001 |
| legionella_pneumophila | aztreonam | 0.8 | 0.003 |
| legionella_pneumophila | erythromycin | 0.8 | 1 |
| legionella_pneumophila | azithromycin | 0.9 | 1 |
| legionella_pneumophila | clarithromycin | 0.8 | 1 |
| legionella_pneumophila | clindamycin | 0.05 | 1 |
| legionella_pneumophila | gentamicin | 0.05 | 1 |
| legionella_pneumophila | tobramycin | 0.05 | 1 |
| legionella_pneumophila | amikacin | 0.05 | 1 |
| legionella_pneumophila | ciprofloxacin | 0.9 | 1 |
| legionella_pneumophila | levofloxacin | 0.95 | 6 |
| legionella_pneumophila | moxifloxacin | 0.9 | 1 |
| legionella_pneumophila | ofloxacin | 0.7 | 1 |
| legionella_pneumophila | tetracycline | 0.8 | 0.25 |
| legionella_pneumophila | doxycycline | 0.85 | 2 |
| legionella_pneumophila | minocycline | 0.9 | 0.25 |
| legionella_pneumophila | tigecycline | 0.1 | 1 |
| legionella_pneumophila | vancomycin | 0.05 | 1 |
| legionella_pneumophila | teicoplanin | 0.05 | 1 |
| legionella_pneumophila | dalbavancin | 0.05 | 0.005 |
| legionella_pneumophila | linezolid | 0.05 | 0.005 |
| legionella_pneumophila | tedizolid | 0.05 | 0.005 |
| legionella_pneumophila | daptomycin | 0.1 | 1 |
| legionella_pneumophila | quinu_dalfo | 0.05 | 0.005 |
| legionella_pneumophila | trim_sulf | 0.05 | 0.04 |
| legionella_pneumophila | chloramphenicol | 0.05 | 1 |
| legionella_pneumophila | nitrofurantoin | 0.05 | 1 |
| legionella_pneumophila | fosfomycin | 0.1 | 1 |
| legionella_pneumophila | retapamulin | 0.05 | 1 |
| legionella_pneumophila | fusidic_a | 0.05 | 1 |
| legionella_pneumophila | metronidazole | 0.05 | 1 |
| legionella_pneumophila | fidaxomicin | 0.1 | 1 |
| legionella_pneumophila | furazolidone | 0.05 | 1 |
| legionella_pneumophila | rifampicin | 0.05 | 1 |
| legionella_pneumophila | amoxicillin_clavulanate | 0.05 | 1 |
| legionella_pneumophila | piperacillin_tazobactam | 0.05 | 1 |
| legionella_pneumophila | ampicillin_sulbactam | 0.05 | 1 |
| legionella_pneumophila | ticarcillin_clavulanate | 0.05 | 1 |
| legionella_pneumophila | ceftazidime_avibactam | 0.05 | 0.005 |
| legionella_pneumophila | meropenem_vaborbactam | 0.05 | 0.005 |
| legionella_pneumophila | colistin | 0.05 | 0.005 |
| legionella_pneumophila | flucloxacillin | 0.01 | 1 |
| legionella_pneumophila | aztreonam_avibactam | 0.01 | 0.003 |
| legionella_pneumophila | cefixime | 0.8 | 0.2 |
| burkholderia_cepacia_complex | sulfanilamide | 0.1 | 0.02 |
| burkholderia_cepacia_complex | penicillin_g | 0.05 | 1 |
| burkholderia_cepacia_complex | ampicillin | 0.05 | 1 |
| burkholderia_cepacia_complex | amoxicillin | 0.05 | 1 |
| burkholderia_cepacia_complex | piperacillin | 0.6 | 1 |
| burkholderia_cepacia_complex | ticarcillin | 0.5 | 1 |
| burkholderia_cepacia_complex | cephalexin | 0.05 | 0.3 |
| burkholderia_cepacia_complex | cefazolin | 0.05 | 0.3 |
| burkholderia_cepacia_complex | cefuroxime | 0.1 | 0.3 |
| burkholderia_cepacia_complex | ceftriaxone | 0.1 | 0.2 |
| burkholderia_cepacia_complex | ceftazidime | 0.7 | 0.2 |
| burkholderia_cepacia_complex | cefepime | 0.75 | 0.35 |
| burkholderia_cepacia_complex | ceftaroline | 0.1 | 0.002 |
| burkholderia_cepacia_complex | ceftolozane_tazobactam | 0.1 | 1 |
| burkholderia_cepacia_complex | cefiderocol | 0.1 | 1 |
| burkholderia_cepacia_complex | meropenem | 0.8 | 0.005 |
| burkholderia_cepacia_complex | imipenem_c | 0.8 | 0.005 |
| burkholderia_cepacia_complex | ertapenem | 0.1 | 0.005 |
| burkholderia_cepacia_complex | aztreonam | 0.1 | 0.003 |
| burkholderia_cepacia_complex | gentamicin | 0.7 | 1 |
| burkholderia_cepacia_complex | tobramycin | 0.65 | 1 |
| burkholderia_cepacia_complex | amikacin | 0.75 | 1 |
| burkholderia_cepacia_complex | ciprofloxacin | 0.6 | 1 |
| burkholderia_cepacia_complex | levofloxacin | 0.65 | 1 |
| burkholderia_cepacia_complex | moxifloxacin | 0.6 | 1 |
| burkholderia_cepacia_complex | ofloxacin | 0.6 | 1 |
| burkholderia_cepacia_complex | tetracycline | 0.6 | 0.25 |
| burkholderia_cepacia_complex | doxycycline | 0.65 | 0.25 |
| burkholderia_cepacia_complex | minocycline | 0.7 | 0.25 |
| burkholderia_cepacia_complex | tigecycline | 0.1 | 1 |
| burkholderia_cepacia_complex | dalbavancin | 0 | 0.005 |
| burkholderia_cepacia_complex | linezolid | 0 | 0.005 |
| burkholderia_cepacia_complex | tedizolid | 0 | 0.005 |
| burkholderia_cepacia_complex | daptomycin | 0.1 | 1 |
| burkholderia_cepacia_complex | quinu_dalfo | 0 | 0.005 |
| burkholderia_cepacia_complex | trim_sulf | 0.6 | 0.04 |
| burkholderia_cepacia_complex | chloramphenicol | 0.7 | 1 |
| burkholderia_cepacia_complex | nitrofurantoin | 0.1 | 1 |
| burkholderia_cepacia_complex | fosfomycin | 0.1 | 1 |
| burkholderia_cepacia_complex | fidaxomicin | 0.1 | 1 |
| burkholderia_cepacia_complex | furazolidone | 0.1 | 1 |
| burkholderia_cepacia_complex | rifampicin | 0.5 | 1 |
| burkholderia_cepacia_complex | amoxicillin_clavulanate | 0.05 | 1 |
| burkholderia_cepacia_complex | piperacillin_tazobactam | 0.65 | 1 |
| burkholderia_cepacia_complex | ampicillin_sulbactam | 0.65 | 1 |
| burkholderia_cepacia_complex | ticarcillin_clavulanate | 0.6 | 1 |
| burkholderia_cepacia_complex | ceftazidime_avibactam | 0.65 | 0.005 |
| burkholderia_cepacia_complex | meropenem_vaborbactam | 0.75 | 0.005 |
| burkholderia_cepacia_complex | colistin | 0.8 | 0.005 |
| burkholderia_cepacia_complex | flucloxacillin | 0.01 | 1 |
| burkholderia_cepacia_complex | aztreonam_avibactam | 0.6 | 0.003 |
| burkholderia_cepacia_complex | cefixime | 0.1 | 0.2 |

### B.5 Regional Parameters

Region-level scalars (applicable to all bacteria) and the per-region per-bacteria acquisition log-odds adjustments.

See: [┬º2.5 Travel](#25-travel), [┬º3.1 Community acquisition](#31-community-acquisition).

#### Region Scalars

| Region | Travel mult | Cessation mult | Mortality log-odds | Sepsis log-odds | Sepsis mort mult | Testing mult | Abx init log-odds | Hosp log-odds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| north_america | 3 | 0.85 | 0 | 0.4 | 0.5 | 1.1 | 0 | 0.5 |
| south_america | 0.8 | 1.25 | 0.26 | -0.3 | 1.1 | 0.6 | -0.8 | -0.2 |
| africa | 0.3 | 1.4 | 0.69 | -0.7 | 1.5 | 0.3 | -1.4 | -0.5 |
| asia | 1.5 | 1.15 | 0.18 | 0 | 0.9 | 0.7 | -0.5 | 0 |
| europe | 3.5 | 0.8 | -0.105 | 0.5 | 0.4 | 1.2 | 0 | 0.6 |
| oceania | 2.5 | 0.85 | 0 | 0.3 | 0.5 | 0.8 | 0 | 0.4 |
| home | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 0 |

#### RegionΓÇôBacteria Acquisition Log-Odds

| Region | Bacteria | Acquisition log-odds |
| --- | ---: | ---: |
| south_america | acinetobacter_baumannii | 1.6 |
| south_america | citrobacter_spp. | -0.3 |
| south_america | enterobacter_spp. | -0.2 |
| south_america | enterococcus_faecalis | 1 |
| south_america | enterococcus_faecium | 0.9 |
| south_america | escherichia_coli | 0.6 |
| south_america | klebsiella_pneumoniae | 0.8 |
| south_america | morganella_spp. | -0.5 |
| south_america | proteus_spp. | 0.5 |
| south_america | serratia_spp. | 0.4 |
| south_america | pseudomonas_aeruginosa | 1 |
| south_america | staphylococcus_aureus | 0.8 |
| south_america | streptococcus_pneumoniae | 1.2 |
| south_america | salmonella_enterica_serovar_typhi | 1.6 |
| south_america | salmonella_enterica_serovar_paratyphi_a | 1.5 |
| south_america | invasive_non-typhoidal_salmonella_spp. | 1.7 |
| south_america | neisseria_gonorrhoeae | 1.2 |
| south_america | streptococcus_pyogenes | 1 |
| south_america | streptococcus_agalactiae | 0.7 |
| south_america | haemophilus_influenzae | 1 |
| south_america | chlamydia_trachomatis | 1.1 |
| south_america | vibrio_cholerae | 1.5 |
| south_america | neisseria_meningitidis | -0.2 |
| south_america | listeria_monocytogenes | 0.1 |
| south_america | campylobacter_jejuni | 1.5 |
| south_america | enterobacter_cloacae | -0.2 |
| south_america | moraxella_catarrhalis | 0.2 |
| south_america | treponema_pallidum | 0.1 |
| south_america | helicobacter_pylori | 1.9 |
| africa | acinetobacter_baumannii | 2.3 |
| africa | citrobacter_spp. | 0.1 |
| africa | enterobacter_spp. | 0.2 |
| africa | enterococcus_faecalis | 1.6 |
| africa | enterococcus_faecium | 1.4 |
| africa | escherichia_coli | 1.2 |
| africa | klebsiella_pneumoniae | 1.6 |
| africa | morganella_spp. | -0.3 |
| africa | proteus_spp. | 1.1 |
| africa | serratia_spp. | 0.8 |
| africa | pseudomonas_aeruginosa | 1.6 |
| africa | staphylococcus_aureus | 1.5 |
| africa | streptococcus_pneumoniae | 2.2 |
| africa | salmonella_enterica_serovar_typhi | 3.3 |
| africa | salmonella_enterica_serovar_paratyphi_a | 3.2 |
| africa | invasive_non-typhoidal_salmonella_spp. | 4.5 |
| africa | neisseria_gonorrhoeae | 2.1 |
| africa | streptococcus_pyogenes | 1.8 |
| africa | streptococcus_agalactiae | 1.4 |
| africa | haemophilus_influenzae | 2 |
| africa | chlamydia_trachomatis | 2.1 |
| africa | vibrio_cholerae | 3.5 |
| africa | neisseria_meningitidis | 1.5 |
| africa | listeria_monocytogenes | 0.4 |
| africa | clostridioides_difficile | 0.2 |
| africa | campylobacter_jejuni | 2.5 |
| africa | enterobacter_cloacae | 0.2 |
| africa | yersinia_enterocolitica | 0.2 |
| africa | moraxella_catarrhalis | 0.9 |
| africa | treponema_pallidum | 0.3 |
| africa | helicobacter_pylori | 3.2 |
| asia | acinetobacter_baumannii | 2 |
| asia | citrobacter_spp. | -0.2 |
| asia | enterococcus_faecalis | 1.3 |
| asia | enterococcus_faecium | 1.2 |
| asia | escherichia_coli | 1 |
| asia | klebsiella_pneumoniae | 1.3 |
| asia | morganella_spp. | -0.4 |
| asia | proteus_spp. | 0.8 |
| asia | serratia_spp. | 0.6 |
| asia | pseudomonas_aeruginosa | 1.3 |
| asia | staphylococcus_aureus | 1.2 |
| asia | streptococcus_pneumoniae | 1.8 |
| asia | salmonella_enterica_serovar_typhi | 3 |
| asia | salmonella_enterica_serovar_paratyphi_a | 2.9 |
| asia | invasive_non-typhoidal_salmonella_spp. | 2.3 |
| asia | neisseria_gonorrhoeae | 1.6 |
| asia | streptococcus_pyogenes | 1.4 |
| asia | streptococcus_agalactiae | 1.1 |
| asia | haemophilus_influenzae | 1.6 |
| asia | chlamydia_trachomatis | 1.7 |
| asia | vibrio_cholerae | 2.8 |
| asia | neisseria_meningitidis | -0.2 |
| asia | listeria_monocytogenes | 0.3 |
| asia | clostridioides_difficile | 0.1 |
| asia | campylobacter_jejuni | 1.9 |
| asia | moraxella_catarrhalis | 0.5 |
| asia | helicobacter_pylori | 2.8 |
| europe | acinetobacter_baumannii | -0.3 |
| europe | citrobacter_spp. | -0.5 |
| europe | enterobacter_spp. | -0.5 |
| europe | enterococcus_faecalis | 0.2 |
| europe | enterococcus_faecium | 0.2 |
| europe | escherichia_coli | -0.3 |
| europe | klebsiella_pneumoniae | -0.1 |
| europe | morganella_spp. | -0.7 |
| europe | proteus_spp. | -0.1 |
| europe | serratia_spp. | -0.1 |
| europe | pseudomonas_aeruginosa | 0.3 |
| europe | staphylococcus_aureus | -0.2 |
| europe | streptococcus_pneumoniae | -0.1 |
| europe | salmonella_enterica_serovar_typhi | -1.8 |
| europe | salmonella_enterica_serovar_paratyphi_a | -2.1 |
| europe | invasive_non-typhoidal_salmonella_spp. | -0.8 |
| europe | neisseria_gonorrhoeae | 0.4 |
| europe | streptococcus_pyogenes | -0.1 |
| europe | haemophilus_influenzae | -0.3 |
| europe | chlamydia_trachomatis | 0.2 |
| europe | vibrio_cholerae | -3 |
| europe | neisseria_meningitidis | -0.8 |
| europe | listeria_monocytogenes | -0.2 |
| europe | clostridioides_difficile | -0.1 |
| europe | campylobacter_jejuni | 1 |
| europe | enterobacter_cloacae | -0.5 |
| europe | moraxella_catarrhalis | -0.3 |
| europe | treponema_pallidum | -0.2 |
| europe | helicobacter_pylori | -0.2 |
| oceania | acinetobacter_baumannii | -0.1 |
| oceania | citrobacter_spp. | -0.4 |
| oceania | enterobacter_spp. | -0.3 |
| oceania | enterococcus_faecalis | 0.2 |
| oceania | enterococcus_faecium | 0.2 |
| oceania | escherichia_coli | 0.1 |
| oceania | morganella_spp. | -0.6 |
| oceania | pseudomonas_aeruginosa | 0.4 |
| oceania | streptococcus_pneumoniae | 0.3 |
| oceania | salmonella_enterica_serovar_typhi | -1.1 |
| oceania | salmonella_enterica_serovar_paratyphi_a | -1.3 |
| oceania | invasive_non-typhoidal_salmonella_spp. | -0.3 |
| oceania | neisseria_gonorrhoeae | 0.5 |
| oceania | streptococcus_pyogenes | 0.2 |
| oceania | streptococcus_agalactiae | 0.1 |
| oceania | haemophilus_influenzae | -0.2 |
| oceania | chlamydia_trachomatis | 0.6 |
| oceania | vibrio_cholerae | -2.2 |
| oceania | neisseria_meningitidis | -0.6 |
| oceania | listeria_monocytogenes | -0.3 |
| oceania | clostridioides_difficile | -0.2 |
| oceania | campylobacter_jejuni | 0.6 |
| oceania | enterobacter_cloacae | -0.3 |
| oceania | yersinia_enterocolitica | -0.1 |
| oceania | moraxella_catarrhalis | -0.2 |
| oceania | treponema_pallidum | -0.3 |
| oceania | helicobacter_pylori | 0.4 |

### B.6 Age-Dependent Parameters

Log-odds adjustments by age category for bacteria acquisition and regional effects. Age categories: infant, preschool, school, young_adult, middle_age, elderly.

See: [┬º2.2 Ageing and age categories](#22-ageing-and-age-categories), [┬º3.1 Community acquisition](#31-community-acquisition).

#### Default Age Log-Odds

| Age category | Default log-odds |
| --- | ---: |
| infant | 1.5 |
| preschool | 0.8 |
| school | 0.3 |
| young_adult | 0 |
| middle_age | 0.2 |
| elderly | 0.9 |

#### BacteriaΓÇôAge Log-Odds

| Bacteria | Age category | Log-odds |
| --- | ---: | ---: |
| acinetobacter_baumannii | infant | 0.5 |
| acinetobacter_baumannii | preschool | -0.5 |
| acinetobacter_baumannii | school | -0.8 |
| acinetobacter_baumannii | young_adult | 0.2 |
| acinetobacter_baumannii | middle_age | 0.8 |
| acinetobacter_baumannii | elderly | 1.5 |
| citrobacter_spp. | infant | 1.5 |
| citrobacter_spp. | preschool | 0.8 |
| citrobacter_spp. | school | 0.3 |
| citrobacter_spp. | middle_age | 0.2 |
| citrobacter_spp. | elderly | 0.9 |
| enterobacter_spp. | infant | 1.5 |
| enterobacter_spp. | preschool | 0.8 |
| enterobacter_spp. | school | 0.3 |
| enterobacter_spp. | middle_age | 0.2 |
| enterobacter_spp. | elderly | 0.9 |
| enterococcus_faecalis | infant | 1.5 |
| enterococcus_faecalis | preschool | 0.8 |
| enterococcus_faecalis | school | 0.3 |
| enterococcus_faecalis | middle_age | 0.2 |
| enterococcus_faecalis | elderly | 0.9 |
| enterococcus_faecium | infant | 1.5 |
| enterococcus_faecium | preschool | 0.8 |
| enterococcus_faecium | school | 0.3 |
| enterococcus_faecium | middle_age | 0.2 |
| enterococcus_faecium | elderly | 0.9 |
| escherichia_coli | infant | 0.6 |
| escherichia_coli | preschool | -0.1 |
| escherichia_coli | school | -0.2 |
| escherichia_coli | young_adult | 0.3 |
| escherichia_coli | middle_age | 0.2 |
| escherichia_coli | elderly | 0.9 |
| klebsiella_pneumoniae | infant | 1.5 |
| klebsiella_pneumoniae | preschool | 0.8 |
| klebsiella_pneumoniae | school | 0.3 |
| klebsiella_pneumoniae | middle_age | 0.2 |
| klebsiella_pneumoniae | elderly | 0.9 |
| morganella_spp. | infant | 1.5 |
| morganella_spp. | preschool | 0.8 |
| morganella_spp. | school | 0.3 |
| morganella_spp. | middle_age | 0.2 |
| morganella_spp. | elderly | 0.9 |
| proteus_spp. | infant | 1.5 |
| proteus_spp. | preschool | 0.8 |
| proteus_spp. | school | 0.3 |
| proteus_spp. | middle_age | 0.2 |
| proteus_spp. | elderly | 0.9 |
| serratia_spp. | infant | 1.5 |
| serratia_spp. | preschool | 0.8 |
| serratia_spp. | school | 0.3 |
| serratia_spp. | middle_age | 0.2 |
| serratia_spp. | elderly | 0.9 |
| p_stuartii | infant | 1.5 |
| p_stuartii | preschool | 0.8 |
| p_stuartii | school | 0.3 |
| p_stuartii | middle_age | 0.2 |
| p_stuartii | elderly | 0.9 |
| pseudomonas_aeruginosa | infant | 1.5 |
| pseudomonas_aeruginosa | preschool | 0.8 |
| pseudomonas_aeruginosa | school | 0.3 |
| pseudomonas_aeruginosa | middle_age | 0.2 |
| pseudomonas_aeruginosa | elderly | 0.9 |
| stenotrophomonas_maltophilia | infant | 1.5 |
| stenotrophomonas_maltophilia | preschool | 0.8 |
| stenotrophomonas_maltophilia | school | 0.3 |
| stenotrophomonas_maltophilia | middle_age | 0.2 |
| stenotrophomonas_maltophilia | elderly | 0.9 |
| staphylococcus_aureus | infant | 1.5 |
| staphylococcus_aureus | preschool | 0.8 |
| staphylococcus_aureus | school | 0.3 |
| staphylococcus_aureus | middle_age | 0.2 |
| staphylococcus_aureus | elderly | 0.9 |
| staphylococcus_epidermidis | infant | 1.5 |
| staphylococcus_epidermidis | preschool | 0.8 |
| staphylococcus_epidermidis | school | 0.3 |
| staphylococcus_epidermidis | middle_age | 0.2 |
| staphylococcus_epidermidis | elderly | 0.9 |
| streptococcus_pneumoniae | infant | 1.7 |
| streptococcus_pneumoniae | preschool | 0.9 |
| streptococcus_pneumoniae | school | 0.2 |
| streptococcus_pneumoniae | young_adult | -0.4 |
| streptococcus_pneumoniae | middle_age | -0.1 |
| streptococcus_pneumoniae | elderly | 1.2 |
| salmonella_enterica_serovar_typhi | infant | 1 |
| salmonella_enterica_serovar_typhi | preschool | 0.8 |
| salmonella_enterica_serovar_typhi | school | 0.5 |
| salmonella_enterica_serovar_typhi | middle_age | 0.2 |
| salmonella_enterica_serovar_typhi | elderly | 0.8 |
| salmonella_enterica_serovar_paratyphi_a | infant | 1.5 |
| salmonella_enterica_serovar_paratyphi_a | preschool | 0.8 |
| salmonella_enterica_serovar_paratyphi_a | school | 0.3 |
| salmonella_enterica_serovar_paratyphi_a | middle_age | 0.2 |
| salmonella_enterica_serovar_paratyphi_a | elderly | 0.9 |
| invasive_non-typhoidal_salmonella_spp. | infant | 1.5 |
| invasive_non-typhoidal_salmonella_spp. | preschool | 0.8 |
| invasive_non-typhoidal_salmonella_spp. | school | 0.3 |
| invasive_non-typhoidal_salmonella_spp. | middle_age | 0.2 |
| invasive_non-typhoidal_salmonella_spp. | elderly | 0.9 |
| shigella_spp. | infant | 1.5 |
| shigella_spp. | preschool | 2 |
| shigella_spp. | school | 1.2 |
| shigella_spp. | young_adult | 0.3 |
| shigella_spp. | elderly | 0.5 |
| neisseria_gonorrhoeae | infant | -2.5 |
| neisseria_gonorrhoeae | preschool | -3.5 |
| neisseria_gonorrhoeae | school | -0.8 |
| neisseria_gonorrhoeae | young_adult | 2 |
| neisseria_gonorrhoeae | middle_age | 0.9 |
| neisseria_gonorrhoeae | elderly | -0.8 |
| streptococcus_pyogenes | infant | 1.5 |
| streptococcus_pyogenes | preschool | 0.8 |
| streptococcus_pyogenes | school | 0.3 |
| streptococcus_pyogenes | middle_age | 0.2 |
| streptococcus_pyogenes | elderly | 0.9 |
| streptococcus_agalactiae | infant | 1.5 |
| streptococcus_agalactiae | preschool | 0.8 |
| streptococcus_agalactiae | school | 0.3 |
| streptococcus_agalactiae | middle_age | 0.2 |
| streptococcus_agalactiae | elderly | 0.9 |
| haemophilus_influenzae | infant | 2.5 |
| haemophilus_influenzae | preschool | 1.5 |
| haemophilus_influenzae | school | 0.8 |
| haemophilus_influenzae | young_adult | -0.5 |
| haemophilus_influenzae | middle_age | -0.2 |
| haemophilus_influenzae | elderly | 1 |
| chlamydia_trachomatis | infant | -3 |
| chlamydia_trachomatis | preschool | -3.5 |
| chlamydia_trachomatis | school | -1 |
| chlamydia_trachomatis | young_adult | 1.7 |
| chlamydia_trachomatis | middle_age | 0.8 |
| chlamydia_trachomatis | elderly | -1 |
| mycoplasma_genitalium | infant | 1.5 |
| mycoplasma_genitalium | preschool | 0.8 |
| mycoplasma_genitalium | school | 0.3 |
| mycoplasma_genitalium | middle_age | 0.2 |
| mycoplasma_genitalium | elderly | 0.9 |
| vibrio_cholerae | infant | 1.5 |
| vibrio_cholerae | preschool | 0.8 |
| vibrio_cholerae | school | 0.3 |
| vibrio_cholerae | middle_age | 0.2 |
| vibrio_cholerae | elderly | 0.9 |
| neisseria_meningitidis | infant | 1.8 |
| neisseria_meningitidis | preschool | 0.4 |
| neisseria_meningitidis | school | 1 |
| neisseria_meningitidis | young_adult | 1.3 |
| neisseria_meningitidis | middle_age | -0.2 |
| neisseria_meningitidis | elderly | 0.2 |
| listeria_monocytogenes | infant | 1.8 |
| listeria_monocytogenes | preschool | -0.5 |
| listeria_monocytogenes | school | -1 |
| listeria_monocytogenes | young_adult | 0.5 |
| listeria_monocytogenes | elderly | 1.5 |
| clostridioides_difficile | infant | -1 |
| clostridioides_difficile | preschool | -1.5 |
| clostridioides_difficile | school | -2 |
| clostridioides_difficile | young_adult | -0.5 |
| clostridioides_difficile | middle_age | 0.5 |
| clostridioides_difficile | elderly | 2 |
| bacteroides_fragilis | infant | 1.5 |
| bacteroides_fragilis | preschool | 0.8 |
| bacteroides_fragilis | school | 0.3 |
| bacteroides_fragilis | middle_age | 0.2 |
| bacteroides_fragilis | elderly | 0.9 |
| campylobacter_jejuni | infant | 1.6 |
| campylobacter_jejuni | preschool | 1.2 |
| campylobacter_jejuni | school | 0.4 |
| campylobacter_jejuni | young_adult | 0.2 |
| campylobacter_jejuni | elderly | 0.5 |
| enterobacter_cloacae | infant | 1.5 |
| enterobacter_cloacae | preschool | 0.8 |
| enterobacter_cloacae | school | 0.3 |
| enterobacter_cloacae | middle_age | 0.2 |
| enterobacter_cloacae | elderly | 0.9 |
| yersinia_enterocolitica | infant | 1.5 |
| yersinia_enterocolitica | preschool | 0.8 |
| yersinia_enterocolitica | school | 0.3 |
| yersinia_enterocolitica | middle_age | 0.2 |
| yersinia_enterocolitica | elderly | 0.9 |
| moraxella_catarrhalis | infant | 1.8 |
| moraxella_catarrhalis | preschool | 1 |
| moraxella_catarrhalis | school | 0.3 |
| moraxella_catarrhalis | young_adult | -0.8 |
| moraxella_catarrhalis | middle_age | -0.3 |
| moraxella_catarrhalis | elderly | 1.2 |
| treponema_pallidum | infant | -2.2 |
| treponema_pallidum | preschool | -4.3 |
| treponema_pallidum | school | -2.4 |
| treponema_pallidum | young_adult | 0.4 |
| treponema_pallidum | middle_age | -0.6 |
| treponema_pallidum | elderly | -2.2 |
| bordetella_pertussis | infant | 1.5 |
| bordetella_pertussis | preschool | 0.8 |
| bordetella_pertussis | school | 0.3 |
| bordetella_pertussis | middle_age | 0.2 |
| bordetella_pertussis | elderly | 0.9 |
| helicobacter_pylori | infant | -2.5 |
| helicobacter_pylori | preschool | -1.5 |
| helicobacter_pylori | school | -0.5 |
| helicobacter_pylori | young_adult | 0.8 |
| helicobacter_pylori | middle_age | 1.4 |
| helicobacter_pylori | elderly | 1.8 |
| mdr_mycobacterium_tuberculosis | infant | 1.5 |
| mdr_mycobacterium_tuberculosis | preschool | 0.8 |
| mdr_mycobacterium_tuberculosis | school | 0.3 |
| mdr_mycobacterium_tuberculosis | middle_age | 0.2 |
| mdr_mycobacterium_tuberculosis | elderly | 0.9 |
| mycoplasma_pneumoniae | infant | 1.5 |
| mycoplasma_pneumoniae | preschool | 0.8 |
| mycoplasma_pneumoniae | school | 0.3 |
| mycoplasma_pneumoniae | middle_age | 0.2 |
| mycoplasma_pneumoniae | elderly | 0.9 |
| legionella_pneumophila | infant | 1.5 |
| legionella_pneumophila | preschool | 0.8 |
| legionella_pneumophila | school | 0.3 |
| legionella_pneumophila | middle_age | 0.2 |
| legionella_pneumophila | elderly | 0.9 |
| burkholderia_cepacia_complex | infant | 1.5 |
| burkholderia_cepacia_complex | preschool | 0.8 |
| burkholderia_cepacia_complex | school | 0.3 |
| burkholderia_cepacia_complex | middle_age | 0.2 |
| burkholderia_cepacia_complex | elderly | 0.9 |

#### RegionΓÇôAge Log-Odds

| Region | Age category | Log-odds |
| --- | ---: | ---: |
| south_america | infant | 1.2 |
| south_america | preschool | 0.7 |
| south_america | school | 0.3 |
| south_america | young_adult | 0.2 |
| south_america | middle_age | 0.3 |
| south_america | elderly | 0.9 |
| africa | infant | 2 |
| africa | preschool | 1.2 |
| africa | school | 0.6 |
| africa | young_adult | 0.3 |
| africa | middle_age | 0.4 |
| africa | elderly | 1.3 |
| asia | infant | 1 |
| asia | preschool | 0.5 |
| asia | school | 0.2 |
| asia | young_adult | 0.1 |
| asia | middle_age | 0.2 |
| asia | elderly | 0.8 |
| europe | infant | -0.2 |
| europe | preschool | -0.1 |
| europe | elderly | 0.2 |
| oceania | infant | 0.1 |
| oceania | elderly | 0.3 |

### B.7 Syndrome Parameters

Infection-site (syndrome) specific parameters. Syndromes are: 1 = UTI, 2 = skin/soft tissue, 3 = respiratory, 4 = bloodstream, 5 = intra-abdominal, 6 = CNS/meningitis, 7 = gastrointestinal, 8 = genital/STI, 9 = bone/joint, 10 = other.

See: [┬º4.1 Syndrome assignment](#41-syndrome-assignment), [┬º6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use), [┬º6.4 Drug penetration by syndrome](#64-drug-penetration-by-syndrome).

#### Syndrome Empiric Drug Scores

| Syndrome | Drug | Empiric score |
| --- | ---: | ---: |
| uti | ampicillin | 10 |
| uti | amoxicillin | 12 |
| uti | cephalexin | 8 |
| uti | cefazolin | 7 |
| uti | cefuroxime | 7 |
| uti | ceftriaxone | 8 |
| uti | ceftazidime | 4 |
| uti | cefepime | 4 |
| uti | meropenem | 4 |
| uti | imipenem_c | 4 |
| uti | ertapenem | 4 |
| uti | ciprofloxacin | 8 |
| uti | levofloxacin | 6 |
| uti | vancomycin | 0.1 |
| uti | linezolid | 0.1 |
| uti | trim_sulf | 1 |
| uti | nitrofurantoin | 8 |
| uti | fosfomycin | 10 |
| uti | amoxicillin_clavulanate | 14 |
| uti | piperacillin_tazobactam | 5 |
| uti | ceftazidime_avibactam | 3 |
| uti | meropenem_vaborbactam | 3 |
| uti | colistin | 0.2 |
| uti | aztreonam_avibactam | 3 |
| uti | cefixime | 7 |
| skin_soft_tissue | penicillin_g | 16 |
| skin_soft_tissue | ampicillin | 13 |
| skin_soft_tissue | amoxicillin | 14 |
| skin_soft_tissue | cephalexin | 13 |
| skin_soft_tissue | cefazolin | 12 |
| skin_soft_tissue | clindamycin | 12 |
| skin_soft_tissue | ciprofloxacin | 4 |
| skin_soft_tissue | doxycycline | 3.5 |
| skin_soft_tissue | minocycline | 3 |
| skin_soft_tissue | vancomycin | 11 |
| skin_soft_tissue | dalbavancin | 9 |
| skin_soft_tissue | linezolid | 10 |
| skin_soft_tissue | tedizolid | 9 |
| skin_soft_tissue | quinu_dalfo | 8 |
| skin_soft_tissue | trim_sulf | 0.5 |
| skin_soft_tissue | rifampicin | 0.5 |
| skin_soft_tissue | amoxicillin_clavulanate | 14 |
| skin_soft_tissue | piperacillin_tazobactam | 3 |
| skin_soft_tissue | flucloxacillin | 15 |
| respiratory | penicillin_g | 16 |
| respiratory | ampicillin | 15 |
| respiratory | amoxicillin | 17 |
| respiratory | cephalexin | 7 |
| respiratory | cefuroxime | 8.5 |
| respiratory | ceftriaxone | 9.5 |
| respiratory | cefepime | 7.5 |
| respiratory | meropenem | 6 |
| respiratory | imipenem_c | 6 |
| respiratory | erythromycin | 9 |
| respiratory | azithromycin | 12 |
| respiratory | clarithromycin | 11 |
| respiratory | levofloxacin | 8 |
| respiratory | moxifloxacin | 8 |
| respiratory | ofloxacin | 6 |
| respiratory | doxycycline | 3 |
| respiratory | minocycline | 2.5 |
| respiratory | vancomycin | 6.5 |
| respiratory | linezolid | 7 |
| respiratory | amoxicillin_clavulanate | 20 |
| respiratory | piperacillin_tazobactam | 8 |
| respiratory | cefixime | 6.5 |
| bloodstream | penicillin_g | 6.5 |
| bloodstream | ampicillin | 10 |
| bloodstream | amoxicillin | 9.5 |
| bloodstream | cephalexin | 4 |
| bloodstream | cefazolin | 6 |
| bloodstream | ceftriaxone | 10 |
| bloodstream | ceftazidime | 11 |
| bloodstream | cefepime | 12 |
| bloodstream | meropenem | 17 |
| bloodstream | imipenem_c | 15 |
| bloodstream | gentamicin | 10 |
| bloodstream | tobramycin | 9 |
| bloodstream | amikacin | 10 |
| bloodstream | ciprofloxacin | 6 |
| bloodstream | levofloxacin | 5.5 |
| bloodstream | vancomycin | 11 |
| bloodstream | dalbavancin | 8 |
| bloodstream | linezolid | 10 |
| bloodstream | tedizolid | 9 |
| bloodstream | quinu_dalfo | 8.5 |
| bloodstream | rifampicin | 0.5 |
| bloodstream | amoxicillin_clavulanate | 16 |
| bloodstream | piperacillin_tazobactam | 13 |
| bloodstream | ampicillin_sulbactam | 16 |
| bloodstream | ceftazidime_avibactam | 12.5 |
| bloodstream | meropenem_vaborbactam | 13 |
| bloodstream | colistin | 0.1 |
| bloodstream | flucloxacillin | 7.5 |
| bloodstream | aztreonam_avibactam | 2 |
| intra_abdominal | ampicillin | 8 |
| intra_abdominal | amoxicillin | 7 |
| intra_abdominal | ceftriaxone | 9 |
| intra_abdominal | ceftazidime | 9 |
| intra_abdominal | cefepime | 9 |
| intra_abdominal | meropenem | 15 |
| intra_abdominal | imipenem_c | 14 |
| intra_abdominal | ertapenem | 11 |
| intra_abdominal | gentamicin | 7 |
| intra_abdominal | amikacin | 7 |
| intra_abdominal | ciprofloxacin | 7 |
| intra_abdominal | levofloxacin | 6.5 |
| intra_abdominal | trim_sulf | 0.1 |
| intra_abdominal | metronidazole | 2.5 |
| intra_abdominal | amoxicillin_clavulanate | 11.5 |
| intra_abdominal | piperacillin_tazobactam | 13 |
| intra_abdominal | ampicillin_sulbactam | 12.5 |
| intra_abdominal | ceftazidime_avibactam | 10 |
| intra_abdominal | meropenem_vaborbactam | 10 |
| intra_abdominal | colistin | 0.1 |
| intra_abdominal | aztreonam_avibactam | 2 |
| cns_meningitis | penicillin_g | 11 |
| cns_meningitis | ampicillin | 13 |
| cns_meningitis | ceftriaxone | 15 |
| cns_meningitis | ceftazidime | 12 |
| cns_meningitis | cefepime | 12 |
| cns_meningitis | meropenem | 11 |
| cns_meningitis | imipenem_c | 10 |
| cns_meningitis | vancomycin | 13 |
| cns_meningitis | linezolid | 10 |
| cns_meningitis | chloramphenicol | 2 |
| cns_meningitis | rifampicin | 1 |
| cns_meningitis | piperacillin_tazobactam | 6 |
| cns_meningitis | cefixime | 1 |
| gastrointestinal | penicillin_g | 5 |
| gastrointestinal | ampicillin | 10 |
| gastrointestinal | amoxicillin | 10 |
| gastrointestinal | cephalexin | 5 |
| gastrointestinal | cefuroxime | 5 |
| gastrointestinal | azithromycin | 12 |
| gastrointestinal | ciprofloxacin | 8 |
| gastrointestinal | levofloxacin | 6 |
| gastrointestinal | doxycycline | 4 |
| gastrointestinal | minocycline | 2.5 |
| gastrointestinal | trim_sulf | 0.5 |
| gastrointestinal | metronidazole | 0.2 |
| gastrointestinal | furazolidone | 0.2 |
| gastrointestinal | rifampicin | 0.5 |
| gastrointestinal | amoxicillin_clavulanate | 11 |
| gastrointestinal | ampicillin_sulbactam | 9 |
| gastrointestinal | cefixime | 4.5 |
| genital_sti | penicillin_g | 14 |
| genital_sti | ampicillin | 9 |
| genital_sti | amoxicillin | 11 |
| genital_sti | cephalexin | 6 |
| genital_sti | cefuroxime | 10 |
| genital_sti | ceftriaxone | 13 |
| genital_sti | azithromycin | 13 |
| genital_sti | clindamycin | 9 |
| genital_sti | ciprofloxacin | 4 |
| genital_sti | levofloxacin | 4 |
| genital_sti | doxycycline | 5.5 |
| genital_sti | trim_sulf | 0.4 |
| genital_sti | metronidazole | 0.25 |
| genital_sti | rifampicin | 0.5 |
| genital_sti | amoxicillin_clavulanate | 12 |
| genital_sti | ampicillin_sulbactam | 8 |
| genital_sti | cefixime | 10.5 |
| bone_joint | penicillin_g | 14 |
| bone_joint | ampicillin | 12 |
| bone_joint | cephalexin | 11 |
| bone_joint | cefazolin | 13 |
| bone_joint | ceftriaxone | 11 |
| bone_joint | meropenem | 7 |
| bone_joint | clindamycin | 10 |
| bone_joint | ciprofloxacin | 9 |
| bone_joint | levofloxacin | 9 |
| bone_joint | vancomycin | 12 |
| bone_joint | dalbavancin | 10 |
| bone_joint | linezolid | 11 |
| bone_joint | tedizolid | 10 |
| bone_joint | trim_sulf | 0.5 |
| bone_joint | rifampicin | 2 |
| bone_joint | piperacillin_tazobactam | 6.5 |
| bone_joint | flucloxacillin | 14 |
| other | ceftriaxone | 8 |
| other | cefepime | 8 |
| other | meropenem | 8 |
| other | imipenem_c | 8 |
| other | azithromycin | 6 |
| other | ciprofloxacin | 7 |
| other | vancomycin | 8 |
| other | linezolid | 7 |
| other | piperacillin_tazobactam | 8 |
| other | aztreonam_avibactam | 2 |

#### Syndrome Drug Penetration

| Syndrome | Drug | Penetration factor |
| --- | ---: | ---: |
| uti | penicillin_g | 0.8 |
| uti | ampicillin | 0.8 |
| uti | amoxicillin | 0.8 |
| uti | piperacillin | 0.8 |
| uti | ticarcillin | 0.8 |
| uti | cephalexin | 0.85 |
| uti | cefazolin | 0.85 |
| uti | cefuroxime | 0.85 |
| uti | ceftriaxone | 0.85 |
| uti | ceftazidime | 0.85 |
| uti | cefepime | 0.85 |
| uti | ceftaroline | 0.85 |
| uti | ceftolozane_tazobactam | 0.85 |
| uti | cefiderocol | 0.85 |
| uti | meropenem | 0.85 |
| uti | imipenem_c | 0.85 |
| uti | ertapenem | 0.85 |
| uti | aztreonam | 0.8 |
| uti | erythromycin | 0.4 |
| uti | azithromycin | 0.4 |
| uti | clarithromycin | 0.4 |
| uti | clindamycin | 0.3 |
| uti | gentamicin | 0.75 |
| uti | tobramycin | 0.75 |
| uti | amikacin | 0.75 |
| uti | tetracycline | 0.5 |
| uti | doxycycline | 0.5 |
| uti | minocycline | 0.5 |
| uti | tigecycline | 0.5 |
| uti | vancomycin | 0.6 |
| uti | teicoplanin | 0.6 |
| uti | dalbavancin | 0.6 |
| uti | linezolid | 0.7 |
| uti | tedizolid | 0.7 |
| uti | daptomycin | 0.1 |
| uti | chloramphenicol | 0.4 |
| uti | metronidazole | 0.5 |
| uti | fidaxomicin | 0 |
| uti | rifampicin | 0.4 |
| uti | amoxicillin_clavulanate | 0.8 |
| uti | piperacillin_tazobactam | 0.8 |
| uti | ampicillin_sulbactam | 0.8 |
| uti | ticarcillin_clavulanate | 0.8 |
| uti | ceftazidime_avibactam | 0.8 |
| uti | meropenem_vaborbactam | 0.8 |
| uti | colistin | 0.7 |
| uti | flucloxacillin | 0.8 |
| uti | aztreonam_avibactam | 0.8 |
| uti | cefixime | 0.9 |
| skin_soft_tissue | penicillin_g | 0.85 |
| skin_soft_tissue | ampicillin | 0.85 |
| skin_soft_tissue | amoxicillin | 0.85 |
| skin_soft_tissue | piperacillin | 0.85 |
| skin_soft_tissue | ticarcillin | 0.85 |
| skin_soft_tissue | cephalexin | 0.8 |
| skin_soft_tissue | cefazolin | 0.85 |
| skin_soft_tissue | cefuroxime | 0.85 |
| skin_soft_tissue | ceftriaxone | 0.85 |
| skin_soft_tissue | ceftazidime | 0.85 |
| skin_soft_tissue | cefepime | 0.85 |
| skin_soft_tissue | ceftaroline | 0.85 |
| skin_soft_tissue | ceftolozane_tazobactam | 0.85 |
| skin_soft_tissue | cefiderocol | 0.85 |
| skin_soft_tissue | meropenem | 0.85 |
| skin_soft_tissue | imipenem_c | 0.85 |
| skin_soft_tissue | ertapenem | 0.85 |
| skin_soft_tissue | aztreonam | 0.75 |
| skin_soft_tissue | erythromycin | 0.8 |
| skin_soft_tissue | azithromycin | 0.8 |
| skin_soft_tissue | clarithromycin | 0.8 |
| skin_soft_tissue | clindamycin | 0.85 |
| skin_soft_tissue | gentamicin | 0.6 |
| skin_soft_tissue | tobramycin | 0.6 |
| skin_soft_tissue | amikacin | 0.6 |
| skin_soft_tissue | ciprofloxacin | 0.9 |
| skin_soft_tissue | levofloxacin | 0.9 |
| skin_soft_tissue | moxifloxacin | 0.9 |
| skin_soft_tissue | ofloxacin | 0.9 |
| skin_soft_tissue | tetracycline | 0.8 |
| skin_soft_tissue | doxycycline | 0.8 |
| skin_soft_tissue | minocycline | 0.8 |
| skin_soft_tissue | tigecycline | 0.8 |
| skin_soft_tissue | vancomycin | 0.75 |
| skin_soft_tissue | teicoplanin | 0.75 |
| skin_soft_tissue | dalbavancin | 0.75 |
| skin_soft_tissue | linezolid | 0.9 |
| skin_soft_tissue | tedizolid | 0.9 |
| skin_soft_tissue | daptomycin | 0.95 |
| skin_soft_tissue | trim_sulf | 0.8 |
| skin_soft_tissue | chloramphenicol | 0.7 |
| skin_soft_tissue | nitrofurantoin | 0.2 |
| skin_soft_tissue | fosfomycin | 0.5 |
| skin_soft_tissue | metronidazole | 0.75 |
| skin_soft_tissue | fidaxomicin | 0 |
| skin_soft_tissue | rifampicin | 0.8 |
| skin_soft_tissue | amoxicillin_clavulanate | 0.85 |
| skin_soft_tissue | piperacillin_tazobactam | 0.85 |
| skin_soft_tissue | ampicillin_sulbactam | 0.85 |
| skin_soft_tissue | ticarcillin_clavulanate | 0.85 |
| skin_soft_tissue | ceftazidime_avibactam | 0.85 |
| skin_soft_tissue | meropenem_vaborbactam | 0.85 |
| skin_soft_tissue | colistin | 0.5 |
| skin_soft_tissue | flucloxacillin | 0.85 |
| skin_soft_tissue | aztreonam_avibactam | 0.75 |
| skin_soft_tissue | cefixime | 0.75 |
| respiratory | penicillin_g | 0.65 |
| respiratory | ampicillin | 0.65 |
| respiratory | amoxicillin | 0.65 |
| respiratory | piperacillin | 0.65 |
| respiratory | ticarcillin | 0.65 |
| respiratory | cephalexin | 0.55 |
| respiratory | cefazolin | 0.7 |
| respiratory | cefuroxime | 0.7 |
| respiratory | ceftriaxone | 0.7 |
| respiratory | ceftazidime | 0.7 |
| respiratory | cefepime | 0.7 |
| respiratory | ceftaroline | 0.7 |
| respiratory | ceftolozane_tazobactam | 0.7 |
| respiratory | cefiderocol | 0.7 |
| respiratory | meropenem | 0.75 |
| respiratory | imipenem_c | 0.75 |
| respiratory | ertapenem | 0.75 |
| respiratory | aztreonam | 0.6 |
| respiratory | erythromycin | 0.95 |
| respiratory | azithromycin | 0.95 |
| respiratory | clarithromycin | 0.95 |
| respiratory | clindamycin | 0.75 |
| respiratory | gentamicin | 0.4 |
| respiratory | tobramycin | 0.4 |
| respiratory | amikacin | 0.4 |
| respiratory | ciprofloxacin | 0.95 |
| respiratory | levofloxacin | 0.95 |
| respiratory | moxifloxacin | 0.95 |
| respiratory | ofloxacin | 0.95 |
| respiratory | tetracycline | 0.7 |
| respiratory | doxycycline | 0.7 |
| respiratory | minocycline | 0.7 |
| respiratory | tigecycline | 0.7 |
| respiratory | vancomycin | 0.5 |
| respiratory | teicoplanin | 0.5 |
| respiratory | dalbavancin | 0.5 |
| respiratory | linezolid | 0.9 |
| respiratory | tedizolid | 0.9 |
| respiratory | daptomycin | 0 |
| respiratory | trim_sulf | 0.8 |
| respiratory | chloramphenicol | 0.7 |
| respiratory | nitrofurantoin | 0.15 |
| respiratory | fosfomycin | 0.4 |
| respiratory | metronidazole | 0.6 |
| respiratory | fidaxomicin | 0 |
| respiratory | rifampicin | 0.85 |
| respiratory | amoxicillin_clavulanate | 0.65 |
| respiratory | piperacillin_tazobactam | 0.65 |
| respiratory | ampicillin_sulbactam | 0.65 |
| respiratory | ticarcillin_clavulanate | 0.65 |
| respiratory | ceftazidime_avibactam | 0.65 |
| respiratory | meropenem_vaborbactam | 0.65 |
| respiratory | colistin | 0.3 |
| respiratory | flucloxacillin | 0.65 |
| respiratory | aztreonam_avibactam | 0.6 |
| respiratory | cefixime | 0.6 |
| intra_abdominal | penicillin_g | 0.6 |
| intra_abdominal | ampicillin | 0.6 |
| intra_abdominal | amoxicillin | 0.6 |
| intra_abdominal | piperacillin | 0.6 |
| intra_abdominal | ticarcillin | 0.6 |
| intra_abdominal | cephalexin | 0.45 |
| intra_abdominal | cefazolin | 0.65 |
| intra_abdominal | cefuroxime | 0.65 |
| intra_abdominal | ceftriaxone | 0.65 |
| intra_abdominal | ceftazidime | 0.65 |
| intra_abdominal | cefepime | 0.65 |
| intra_abdominal | ceftaroline | 0.65 |
| intra_abdominal | ceftolozane_tazobactam | 0.65 |
| intra_abdominal | cefiderocol | 0.65 |
| intra_abdominal | meropenem | 0.75 |
| intra_abdominal | imipenem_c | 0.75 |
| intra_abdominal | ertapenem | 0.75 |
| intra_abdominal | aztreonam | 0.55 |
| intra_abdominal | erythromycin | 0.5 |
| intra_abdominal | azithromycin | 0.5 |
| intra_abdominal | clarithromycin | 0.5 |
| intra_abdominal | clindamycin | 0.65 |
| intra_abdominal | gentamicin | 0.3 |
| intra_abdominal | tobramycin | 0.3 |
| intra_abdominal | amikacin | 0.3 |
| intra_abdominal | ciprofloxacin | 0.75 |
| intra_abdominal | levofloxacin | 0.75 |
| intra_abdominal | moxifloxacin | 0.75 |
| intra_abdominal | ofloxacin | 0.75 |
| intra_abdominal | tetracycline | 0.55 |
| intra_abdominal | doxycycline | 0.55 |
| intra_abdominal | minocycline | 0.55 |
| intra_abdominal | tigecycline | 0.55 |
| intra_abdominal | vancomycin | 0.45 |
| intra_abdominal | teicoplanin | 0.45 |
| intra_abdominal | dalbavancin | 0.45 |
| intra_abdominal | linezolid | 0.7 |
| intra_abdominal | tedizolid | 0.7 |
| intra_abdominal | daptomycin | 0.6 |
| intra_abdominal | trim_sulf | 0.6 |
| intra_abdominal | chloramphenicol | 0.6 |
| intra_abdominal | nitrofurantoin | 0.15 |
| intra_abdominal | fosfomycin | 0.5 |
| intra_abdominal | metronidazole | 0.9 |
| intra_abdominal | fidaxomicin | 0.05 |
| intra_abdominal | rifampicin | 0.65 |
| intra_abdominal | amoxicillin_clavulanate | 0.65 |
| intra_abdominal | piperacillin_tazobactam | 0.65 |
| intra_abdominal | ampicillin_sulbactam | 0.65 |
| intra_abdominal | ticarcillin_clavulanate | 0.65 |
| intra_abdominal | ceftazidime_avibactam | 0.65 |
| intra_abdominal | meropenem_vaborbactam | 0.65 |
| intra_abdominal | colistin | 0.35 |
| intra_abdominal | flucloxacillin | 0.6 |
| intra_abdominal | aztreonam_avibactam | 0.55 |
| intra_abdominal | cefixime | 0.55 |
| cns_meningitis | penicillin_g | 0.15 |
| cns_meningitis | ampicillin | 0.15 |
| cns_meningitis | amoxicillin | 0.15 |
| cns_meningitis | piperacillin | 0.15 |
| cns_meningitis | ticarcillin | 0.15 |
| cns_meningitis | cephalexin | 0.05 |
| cns_meningitis | cefazolin | 0.2 |
| cns_meningitis | cefuroxime | 0.2 |
| cns_meningitis | ceftriaxone | 0.35 |
| cns_meningitis | ceftazidime | 0.2 |
| cns_meningitis | cefepime | 0.2 |
| cns_meningitis | ceftaroline | 0.2 |
| cns_meningitis | ceftolozane_tazobactam | 0.2 |
| cns_meningitis | cefiderocol | 0.2 |
| cns_meningitis | meropenem | 0.35 |
| cns_meningitis | imipenem_c | 0.25 |
| cns_meningitis | ertapenem | 0.25 |
| cns_meningitis | aztreonam | 0.1 |
| cns_meningitis | erythromycin | 0.15 |
| cns_meningitis | azithromycin | 0.15 |
| cns_meningitis | clarithromycin | 0.15 |
| cns_meningitis | clindamycin | 0.15 |
| cns_meningitis | gentamicin | 0.05 |
| cns_meningitis | tobramycin | 0.05 |
| cns_meningitis | amikacin | 0.05 |
| cns_meningitis | ciprofloxacin | 0.5 |
| cns_meningitis | levofloxacin | 0.5 |
| cns_meningitis | moxifloxacin | 0.6 |
| cns_meningitis | ofloxacin | 0.5 |
| cns_meningitis | tetracycline | 0.25 |
| cns_meningitis | doxycycline | 0.25 |
| cns_meningitis | minocycline | 0.4 |
| cns_meningitis | tigecycline | 0.25 |
| cns_meningitis | vancomycin | 0.15 |
| cns_meningitis | teicoplanin | 0.15 |
| cns_meningitis | dalbavancin | 0.15 |
| cns_meningitis | linezolid | 0.7 |
| cns_meningitis | tedizolid | 0.7 |
| cns_meningitis | daptomycin | 0.05 |
| cns_meningitis | trim_sulf | 0.5 |
| cns_meningitis | chloramphenicol | 0.7 |
| cns_meningitis | nitrofurantoin | 0.05 |
| cns_meningitis | fosfomycin | 0.3 |
| cns_meningitis | metronidazole | 0.8 |
| cns_meningitis | fidaxomicin | 0 |
| cns_meningitis | rifampicin | 0.5 |
| cns_meningitis | amoxicillin_clavulanate | 0.15 |
| cns_meningitis | piperacillin_tazobactam | 0.15 |
| cns_meningitis | ampicillin_sulbactam | 0.15 |
| cns_meningitis | ticarcillin_clavulanate | 0.15 |
| cns_meningitis | ceftazidime_avibactam | 0.15 |
| cns_meningitis | meropenem_vaborbactam | 0.15 |
| cns_meningitis | colistin | 0.05 |
| cns_meningitis | flucloxacillin | 0.15 |
| cns_meningitis | aztreonam_avibactam | 0.1 |
| cns_meningitis | cefixime | 0.1 |
| gastrointestinal | penicillin_g | 0.55 |
| gastrointestinal | ampicillin | 0.55 |
| gastrointestinal | amoxicillin | 0.55 |
| gastrointestinal | piperacillin | 0.55 |
| gastrointestinal | ticarcillin | 0.55 |
| gastrointestinal | cephalexin | 0.5 |
| gastrointestinal | cefazolin | 0.6 |
| gastrointestinal | cefuroxime | 0.6 |
| gastrointestinal | ceftriaxone | 0.6 |
| gastrointestinal | ceftazidime | 0.6 |
| gastrointestinal | cefepime | 0.6 |
| gastrointestinal | ceftaroline | 0.6 |
| gastrointestinal | ceftolozane_tazobactam | 0.6 |
| gastrointestinal | cefiderocol | 0.6 |
| gastrointestinal | meropenem | 0.65 |
| gastrointestinal | imipenem_c | 0.65 |
| gastrointestinal | ertapenem | 0.65 |
| gastrointestinal | aztreonam | 0.5 |
| gastrointestinal | erythromycin | 0.7 |
| gastrointestinal | azithromycin | 0.7 |
| gastrointestinal | clarithromycin | 0.7 |
| gastrointestinal | clindamycin | 0.65 |
| gastrointestinal | gentamicin | 0.4 |
| gastrointestinal | tobramycin | 0.4 |
| gastrointestinal | amikacin | 0.4 |
| gastrointestinal | ciprofloxacin | 0.85 |
| gastrointestinal | levofloxacin | 0.85 |
| gastrointestinal | moxifloxacin | 0.85 |
| gastrointestinal | ofloxacin | 0.85 |
| gastrointestinal | tetracycline | 0.6 |
| gastrointestinal | doxycycline | 0.6 |
| gastrointestinal | minocycline | 0.6 |
| gastrointestinal | tigecycline | 0.6 |
| gastrointestinal | vancomycin | 0.9 |
| gastrointestinal | teicoplanin | 0.35 |
| gastrointestinal | dalbavancin | 0.35 |
| gastrointestinal | linezolid | 0.75 |
| gastrointestinal | tedizolid | 0.75 |
| gastrointestinal | daptomycin | 0.3 |
| gastrointestinal | trim_sulf | 0.7 |
| gastrointestinal | chloramphenicol | 0.65 |
| gastrointestinal | nitrofurantoin | 0.25 |
| gastrointestinal | fosfomycin | 0.4 |
| gastrointestinal | metronidazole | 0.95 |
| gastrointestinal | furazolidone | 0.9 |
| gastrointestinal | rifampicin | 0.6 |
| gastrointestinal | amoxicillin_clavulanate | 0.55 |
| gastrointestinal | piperacillin_tazobactam | 0.55 |
| gastrointestinal | ampicillin_sulbactam | 0.55 |
| gastrointestinal | ticarcillin_clavulanate | 0.55 |
| gastrointestinal | ceftazidime_avibactam | 0.55 |
| gastrointestinal | meropenem_vaborbactam | 0.55 |
| gastrointestinal | colistin | 0.4 |
| gastrointestinal | flucloxacillin | 0.55 |
| gastrointestinal | aztreonam_avibactam | 0.5 |
| gastrointestinal | cefixime | 0.55 |
| genital_sti | penicillin_g | 0.55 |
| genital_sti | ampicillin | 0.55 |
| genital_sti | amoxicillin | 0.55 |
| genital_sti | piperacillin | 0.55 |
| genital_sti | ticarcillin | 0.55 |
| genital_sti | cephalexin | 0.45 |
| genital_sti | cefazolin | 0.55 |
| genital_sti | cefuroxime | 0.55 |
| genital_sti | ceftriaxone | 0.55 |
| genital_sti | ceftazidime | 0.55 |
| genital_sti | cefepime | 0.55 |
| genital_sti | ceftaroline | 0.55 |
| genital_sti | ceftolozane_tazobactam | 0.55 |
| genital_sti | cefiderocol | 0.55 |
| genital_sti | meropenem | 0.6 |
| genital_sti | imipenem_c | 0.6 |
| genital_sti | ertapenem | 0.6 |
| genital_sti | aztreonam | 0.45 |
| genital_sti | erythromycin | 0.75 |
| genital_sti | azithromycin | 0.75 |
| genital_sti | clarithromycin | 0.75 |
| genital_sti | clindamycin | 0.6 |
| genital_sti | gentamicin | 0.35 |
| genital_sti | tobramycin | 0.35 |
| genital_sti | amikacin | 0.35 |
| genital_sti | ciprofloxacin | 0.9 |
| genital_sti | levofloxacin | 0.9 |
| genital_sti | moxifloxacin | 0.9 |
| genital_sti | ofloxacin | 0.9 |
| genital_sti | tetracycline | 0.75 |
| genital_sti | doxycycline | 0.75 |
| genital_sti | minocycline | 0.75 |
| genital_sti | tigecycline | 0.75 |
| genital_sti | vancomycin | 0.4 |
| genital_sti | teicoplanin | 0.4 |
| genital_sti | dalbavancin | 0.4 |
| genital_sti | linezolid | 0.7 |
| genital_sti | tedizolid | 0.7 |
| genital_sti | daptomycin | 0.4 |
| genital_sti | trim_sulf | 0.8 |
| genital_sti | chloramphenicol | 0.55 |
| genital_sti | nitrofurantoin | 0.3 |
| genital_sti | fosfomycin | 0.5 |
| genital_sti | metronidazole | 0.8 |
| genital_sti | fidaxomicin | 0 |
| genital_sti | rifampicin | 0.6 |
| genital_sti | amoxicillin_clavulanate | 0.55 |
| genital_sti | piperacillin_tazobactam | 0.55 |
| genital_sti | ampicillin_sulbactam | 0.55 |
| genital_sti | ticarcillin_clavulanate | 0.55 |
| genital_sti | ceftazidime_avibactam | 0.55 |
| genital_sti | meropenem_vaborbactam | 0.55 |
| genital_sti | colistin | 0.3 |
| genital_sti | flucloxacillin | 0.55 |
| genital_sti | aztreonam_avibactam | 0.45 |
| genital_sti | cefixime | 0.5 |
| bone_joint | penicillin_g | 0.4 |
| bone_joint | ampicillin | 0.4 |
| bone_joint | amoxicillin | 0.4 |
| bone_joint | piperacillin | 0.4 |
| bone_joint | ticarcillin | 0.4 |
| bone_joint | cephalexin | 0.3 |
| bone_joint | cefazolin | 0.45 |
| bone_joint | cefuroxime | 0.45 |
| bone_joint | ceftriaxone | 0.45 |
| bone_joint | ceftazidime | 0.45 |
| bone_joint | cefepime | 0.45 |
| bone_joint | ceftaroline | 0.45 |
| bone_joint | ceftolozane_tazobactam | 0.45 |
| bone_joint | cefiderocol | 0.45 |
| bone_joint | meropenem | 0.5 |
| bone_joint | imipenem_c | 0.5 |
| bone_joint | ertapenem | 0.5 |
| bone_joint | aztreonam | 0.35 |
| bone_joint | erythromycin | 0.4 |
| bone_joint | azithromycin | 0.4 |
| bone_joint | clarithromycin | 0.4 |
| bone_joint | clindamycin | 0.6 |
| bone_joint | gentamicin | 0.25 |
| bone_joint | tobramycin | 0.25 |
| bone_joint | amikacin | 0.25 |
| bone_joint | ciprofloxacin | 0.7 |
| bone_joint | levofloxacin | 0.7 |
| bone_joint | moxifloxacin | 0.7 |
| bone_joint | ofloxacin | 0.7 |
| bone_joint | tetracycline | 0.5 |
| bone_joint | doxycycline | 0.5 |
| bone_joint | minocycline | 0.5 |
| bone_joint | tigecycline | 0.5 |
| bone_joint | vancomycin | 0.35 |
| bone_joint | teicoplanin | 0.35 |
| bone_joint | dalbavancin | 0.35 |
| bone_joint | linezolid | 0.75 |
| bone_joint | tedizolid | 0.75 |
| bone_joint | daptomycin | 0.5 |
| bone_joint | trim_sulf | 0.55 |
| bone_joint | chloramphenicol | 0.5 |
| bone_joint | nitrofurantoin | 0.1 |
| bone_joint | fosfomycin | 0.6 |
| bone_joint | metronidazole | 0.55 |
| bone_joint | fidaxomicin | 0 |
| bone_joint | rifampicin | 0.8 |
| bone_joint | amoxicillin_clavulanate | 0.4 |
| bone_joint | piperacillin_tazobactam | 0.4 |
| bone_joint | ampicillin_sulbactam | 0.4 |
| bone_joint | ticarcillin_clavulanate | 0.4 |
| bone_joint | ceftazidime_avibactam | 0.4 |
| bone_joint | meropenem_vaborbactam | 0.4 |
| bone_joint | colistin | 0.2 |
| bone_joint | flucloxacillin | 0.4 |
| bone_joint | aztreonam_avibactam | 0.35 |
| bone_joint | cefixime | 0.4 |
| other | gentamicin | 0.7 |
| other | tobramycin | 0.7 |
| other | amikacin | 0.7 |
| other | nitrofurantoin | 0.3 |

### B.8 Clearance Parameters

Infection clearance model parameters. The clearance hazard is a logistic function of base log-odds, per-bacteria adjustments, age effects, immunodeficiency, bacteria level, and treatment duration.

See: [┬º4.4 Natural clearance](#44-natural-clearance).

| Parameter | Value |
| --- | ---: |
| base_clearance_log_odds | -4.2 |
| immunodeficient_log_odds_adjustment | -0.69 |

#### Per-Bacteria Clearance Adjustments

| Bacteria | Log-odds adjustment |
| --- | ---: |

### B.9 Immunodeficiency, Sex, and Vaccination Parameters

See: [┬º2.3 Immunodeficiency](#23-immunodeficiency), [┬º10 Mortality](#10-mortality).

#### Immunodeficiency

| Parameter | Value |
| --- | ---: |
| startup_seed_fraction | 0.05 |
| temporary_onset_rate_per_day | 5e-5 |
| temporary_recovery_rate_per_day | 0.01 |
| chronic_onset_rate_per_day | 6e-5 |
| chronic_recovery_rate_per_day | 0.0012 |
| chronic_probability_age_0_1 | 0.3 |
| chronic_probability_age_1_18 | 0.2 |
| chronic_probability_age_18_65 | 0.4 |
| chronic_probability_age_65_plus | 0.6 |

#### Sex

| Sex | Mortality log-odds |
| --- | ---: |
| male | 0.095 |
| female | -0.105 |

#### Vaccination

| Vaccine | Availability year | Target birth-cohort coverage | Rollout years |
| --- | ---: | ---: | ---: |
| pneumococcal | 1977 | 0.75 | 20 |
| meningococcal | 1981 | 0.55 | 20 |
| hib | 1985 | 0.85 | 15 |

### B.10 Resistance Mechanisms

Parameters for the 40 resistance mechanisms modelled. Each mechanism has a per-day reversion rate, per-drug-class enhancement multipliers, and per-bacteria emergence rates.

See: [┬º7.1 Resistance mechanisms](#71-resistance-mechanisms), [┬º7.2 MechanismΓÇôdrug-class enhancement](#72-mechanismdrug-class-enhancement-multipliers), [┬º7.3 Resistance emergence](#73-resistance-emergence), [┬º7.4 Resistance reversion](#74-resistance-reversion-and-fitness-costs).

#### Mechanism Reversion Rates

| Mechanism | Reversion rate/day |
| --- | ---: |
| enzyme_esbl_ctx_m | 6e-4 |
| enzyme_esbl_tem | 6e-4 |
| enzyme_esbl_shv | 6e-4 |
| enzyme_kpc | 0.001 |
| enzyme_ndm_vim | 0.0015 |
| enzyme_oxa_48 | 5e-4 |
| enzyme_ampc_cmy | 1e-4 |
| enzyme_ampc_dha | 6e-4 |
| target_site_pbp2a_meca | 9e-4 |
| target_site_van_a | 0.002 |
| target_site_van_b | 0.002 |
| mutation_gyra_primary | 1e-4 |
| mutation_gyra_parc_secondary | 2e-4 |
| protection_qnr | 1e-4 |
| enzyme_16s_rrmt | 5e-4 |
| target_site_erm_b | 0.002 |
| target_site_cfr | 5e-4 |
| enzyme_cat | 5e-4 |
| efflux_acrab_tolc | 5e-4 |
| efflux_mexxy_oprm | 5e-4 |
| porin_loss_ompk35_36 | 5e-4 |
| porin_loss_oprd | 5e-4 |
| modification_mcr_1 | 0.0015 |
| global_efflux_pump | 5e-4 |
| global_porin_loss | 5e-4 |
| mutation_folate_pathway | 1e-4 |
| mutation_nitroreductase | 3e-4 |
| enzyme_fos_a | 5e-4 |
| mutation_mpr_f | 0.001 |
| mutation_rpo_b | 0.002 |
| protection_fus_b | 5e-4 |
| protection_tet_m | 5e-4 |
| enzyme_aac_aph | 1e-4 |
| enzyme_bla_z | 1e-4 |
| enzyme_oxa_acinetobacter | 1e-4 |
| mutation_23s_rrna | 1e-4 |
| efflux_tet_abc | 1e-4 |
| mutation_pbp_mosaic | 0.001 |
| efflux_mtr_cde | 0.001 |
| as_yet_unknown | 0.001 |

#### Mechanism Enhancement Multipliers by Drug Class

How much resistance each mechanism confers against each drug class. Only non-zero entries shown.

| Mechanism | Drug class | Enhancement multiplier |
| --- | ---: | ---: |
| enzyme_esbl_ctx_m | pen | 0.9 |
| enzyme_esbl_ctx_m | bli | 0.25 |
| enzyme_esbl_ctx_m | bli_anti_pseudomonal | 0.25 |
| enzyme_esbl_ctx_m | bli_sulbactam | 0.25 |
| enzyme_esbl_ctx_m | c1_2g | 0.9 |
| enzyme_esbl_ctx_m | c3g | 0.85 |
| enzyme_esbl_ctx_m | c3g_bli | 0.85 |
| enzyme_esbl_ctx_m | c4g | 0.35 |
| enzyme_esbl_ctx_m | anti_mrsa_ceph | 0.35 |
| enzyme_esbl_ctx_m | siderophore_ceph | 0.35 |
| enzyme_esbl_ctx_m | cft_avi | 0.1 |
| enzyme_esbl_ctx_m | mer_vab | 0.1 |
| enzyme_esbl_ctx_m | azt_avi | 0.1 |
| enzyme_esbl_ctx_m | mono | 0.8 |
| enzyme_esbl_ctx_m | fq | 0.8 |
| enzyme_esbl_ctx_m | ag_group1 | 0.8 |
| enzyme_esbl_ctx_m | ag_group2 | 0.8 |
| enzyme_esbl_ctx_m | mls | 0.8 |
| enzyme_esbl_ctx_m | lincosamides | 0.8 |
| enzyme_esbl_ctx_m | glyc | 0.8 |
| enzyme_esbl_ctx_m | lipoglycopeptides | 0.8 |
| enzyme_esbl_ctx_m | tet | 0.8 |
| enzyme_esbl_ctx_m | glycylcyclines | 0.8 |
| enzyme_esbl_ctx_m | poly | 0.8 |
| enzyme_esbl_ctx_m | oxa | 0.8 |
| enzyme_esbl_ctx_m | chl | 0.8 |
| enzyme_esbl_ctx_m | sulf | 0.8 |
| enzyme_esbl_ctx_m | lipopeptides | 0.8 |
| enzyme_esbl_ctx_m | streptogramins | 0.8 |
| enzyme_esbl_ctx_m | nitrofurans | 0.8 |
| enzyme_esbl_ctx_m | phosphonic_acids | 0.8 |
| enzyme_esbl_ctx_m | nitroimidazoles | 0.8 |
| enzyme_esbl_ctx_m | rifamycins | 0.8 |
| enzyme_esbl_ctx_m | macrocycles | 0.8 |
| enzyme_esbl_ctx_m | steroid_antibacterials | 0.8 |
| enzyme_esbl_ctx_m | pleuromutilins | 0.8 |
| enzyme_esbl_ctx_m | other | 0.8 |
| enzyme_esbl_tem | pen | 0.85 |
| enzyme_esbl_tem | bli | 0.2 |
| enzyme_esbl_tem | bli_anti_pseudomonal | 0.2 |
| enzyme_esbl_tem | bli_sulbactam | 0.2 |
| enzyme_esbl_tem | c1_2g | 0.85 |
| enzyme_esbl_tem | c3g | 0.65 |
| enzyme_esbl_tem | c3g_bli | 0.65 |
| enzyme_esbl_tem | c4g | 0.25 |
| enzyme_esbl_tem | anti_mrsa_ceph | 0.25 |
| enzyme_esbl_tem | siderophore_ceph | 0.25 |
| enzyme_esbl_tem | cft_avi | 0.1 |
| enzyme_esbl_tem | mer_vab | 0.1 |
| enzyme_esbl_tem | azt_avi | 0.1 |
| enzyme_esbl_tem | mono | 0.6 |
| enzyme_esbl_tem | fq | 0.6 |
| enzyme_esbl_tem | ag_group1 | 0.6 |
| enzyme_esbl_tem | ag_group2 | 0.6 |
| enzyme_esbl_tem | mls | 0.6 |
| enzyme_esbl_tem | lincosamides | 0.6 |
| enzyme_esbl_tem | glyc | 0.6 |
| enzyme_esbl_tem | lipoglycopeptides | 0.6 |
| enzyme_esbl_tem | tet | 0.6 |
| enzyme_esbl_tem | glycylcyclines | 0.6 |
| enzyme_esbl_tem | poly | 0.6 |
| enzyme_esbl_tem | oxa | 0.6 |
| enzyme_esbl_tem | chl | 0.6 |
| enzyme_esbl_tem | sulf | 0.6 |
| enzyme_esbl_tem | lipopeptides | 0.6 |
| enzyme_esbl_tem | streptogramins | 0.6 |
| enzyme_esbl_tem | nitrofurans | 0.6 |
| enzyme_esbl_tem | phosphonic_acids | 0.6 |
| enzyme_esbl_tem | nitroimidazoles | 0.6 |
| enzyme_esbl_tem | rifamycins | 0.6 |
| enzyme_esbl_tem | macrocycles | 0.6 |
| enzyme_esbl_tem | steroid_antibacterials | 0.6 |
| enzyme_esbl_tem | pleuromutilins | 0.6 |
| enzyme_esbl_tem | other | 0.6 |
| enzyme_esbl_shv | pen | 0.8 |
| enzyme_esbl_shv | bli | 0.2 |
| enzyme_esbl_shv | bli_anti_pseudomonal | 0.2 |
| enzyme_esbl_shv | bli_sulbactam | 0.2 |
| enzyme_esbl_shv | c1_2g | 0.85 |
| enzyme_esbl_shv | c3g | 0.65 |
| enzyme_esbl_shv | c3g_bli | 0.65 |
| enzyme_esbl_shv | c4g | 0.3 |
| enzyme_esbl_shv | anti_mrsa_ceph | 0.3 |
| enzyme_esbl_shv | siderophore_ceph | 0.3 |
| enzyme_esbl_shv | cft_avi | 0.1 |
| enzyme_esbl_shv | mer_vab | 0.1 |
| enzyme_esbl_shv | azt_avi | 0.1 |
| enzyme_esbl_shv | mono | 0.55 |
| enzyme_esbl_shv | fq | 0.6 |
| enzyme_esbl_shv | ag_group1 | 0.6 |
| enzyme_esbl_shv | ag_group2 | 0.6 |
| enzyme_esbl_shv | mls | 0.6 |
| enzyme_esbl_shv | lincosamides | 0.6 |
| enzyme_esbl_shv | glyc | 0.6 |
| enzyme_esbl_shv | lipoglycopeptides | 0.6 |
| enzyme_esbl_shv | tet | 0.6 |
| enzyme_esbl_shv | glycylcyclines | 0.6 |
| enzyme_esbl_shv | poly | 0.6 |
| enzyme_esbl_shv | oxa | 0.6 |
| enzyme_esbl_shv | chl | 0.6 |
| enzyme_esbl_shv | sulf | 0.6 |
| enzyme_esbl_shv | lipopeptides | 0.6 |
| enzyme_esbl_shv | streptogramins | 0.6 |
| enzyme_esbl_shv | nitrofurans | 0.6 |
| enzyme_esbl_shv | phosphonic_acids | 0.6 |
| enzyme_esbl_shv | nitroimidazoles | 0.6 |
| enzyme_esbl_shv | rifamycins | 0.6 |
| enzyme_esbl_shv | macrocycles | 0.6 |
| enzyme_esbl_shv | steroid_antibacterials | 0.6 |
| enzyme_esbl_shv | pleuromutilins | 0.6 |
| enzyme_esbl_shv | other | 0.6 |
| enzyme_kpc | pen | 0.95 |
| enzyme_kpc | bli | 0.85 |
| enzyme_kpc | bli_anti_pseudomonal | 0.85 |
| enzyme_kpc | bli_sulbactam | 0.85 |
| enzyme_kpc | c1_2g | 0.95 |
| enzyme_kpc | c3g | 0.95 |
| enzyme_kpc | c3g_bli | 0.95 |
| enzyme_kpc | c4g | 0.85 |
| enzyme_kpc | anti_mrsa_ceph | 0.85 |
| enzyme_kpc | siderophore_ceph | 0.85 |
| enzyme_kpc | cft_avi | 0.3 |
| enzyme_kpc | mer_vab | 0.3 |
| enzyme_kpc | azt_avi | 0.3 |
| enzyme_kpc | carb_group1 | 0.9 |
| enzyme_kpc | carb_group2 | 0.9 |
| enzyme_kpc | mono | 0.9 |
| enzyme_kpc | fq | 0.95 |
| enzyme_kpc | ag_group1 | 0.95 |
| enzyme_kpc | ag_group2 | 0.95 |
| enzyme_kpc | mls | 0.95 |
| enzyme_kpc | lincosamides | 0.95 |
| enzyme_kpc | glyc | 0.95 |
| enzyme_kpc | lipoglycopeptides | 0.95 |
| enzyme_kpc | tet | 0.95 |
| enzyme_kpc | glycylcyclines | 0.95 |
| enzyme_kpc | poly | 0.95 |
| enzyme_kpc | oxa | 0.95 |
| enzyme_kpc | chl | 0.95 |
| enzyme_kpc | sulf | 0.95 |
| enzyme_kpc | lipopeptides | 0.95 |
| enzyme_kpc | streptogramins | 0.95 |
| enzyme_kpc | nitrofurans | 0.95 |
| enzyme_kpc | phosphonic_acids | 0.95 |
| enzyme_kpc | nitroimidazoles | 0.95 |
| enzyme_kpc | rifamycins | 0.95 |
| enzyme_kpc | macrocycles | 0.95 |
| enzyme_kpc | steroid_antibacterials | 0.95 |
| enzyme_kpc | pleuromutilins | 0.95 |
| enzyme_kpc | other | 0.95 |
| enzyme_ndm_vim | pen | 0.95 |
| enzyme_ndm_vim | bli | 0.95 |
| enzyme_ndm_vim | bli_anti_pseudomonal | 0.95 |
| enzyme_ndm_vim | bli_sulbactam | 0.95 |
| enzyme_ndm_vim | c1_2g | 0.95 |
| enzyme_ndm_vim | c3g | 0.95 |
| enzyme_ndm_vim | c3g_bli | 0.95 |
| enzyme_ndm_vim | c4g | 0.9 |
| enzyme_ndm_vim | anti_mrsa_ceph | 0.9 |
| enzyme_ndm_vim | siderophore_ceph | 0.9 |
| enzyme_ndm_vim | cft_avi | 0.5 |
| enzyme_ndm_vim | mer_vab | 0.5 |
| enzyme_ndm_vim | azt_avi | 0.5 |
| enzyme_ndm_vim | carb_group1 | 0.95 |
| enzyme_ndm_vim | carb_group2 | 0.95 |
| enzyme_ndm_vim | mono | 0.1 |
| enzyme_ndm_vim | fq | 0.95 |
| enzyme_ndm_vim | ag_group1 | 0.95 |
| enzyme_ndm_vim | ag_group2 | 0.95 |
| enzyme_ndm_vim | mls | 0.95 |
| enzyme_ndm_vim | lincosamides | 0.95 |
| enzyme_ndm_vim | glyc | 0.95 |
| enzyme_ndm_vim | lipoglycopeptides | 0.95 |
| enzyme_ndm_vim | tet | 0.95 |
| enzyme_ndm_vim | glycylcyclines | 0.95 |
| enzyme_ndm_vim | poly | 0.95 |
| enzyme_ndm_vim | oxa | 0.95 |
| enzyme_ndm_vim | chl | 0.95 |
| enzyme_ndm_vim | sulf | 0.95 |
| enzyme_ndm_vim | lipopeptides | 0.95 |
| enzyme_ndm_vim | streptogramins | 0.95 |
| enzyme_ndm_vim | nitrofurans | 0.95 |
| enzyme_ndm_vim | phosphonic_acids | 0.95 |
| enzyme_ndm_vim | nitroimidazoles | 0.95 |
| enzyme_ndm_vim | rifamycins | 0.95 |
| enzyme_ndm_vim | macrocycles | 0.95 |
| enzyme_ndm_vim | steroid_antibacterials | 0.95 |
| enzyme_ndm_vim | pleuromutilins | 0.95 |
| enzyme_ndm_vim | other | 0.95 |
| enzyme_oxa_48 | pen | 0.8 |
| enzyme_oxa_48 | bli | 0.5 |
| enzyme_oxa_48 | bli_anti_pseudomonal | 0.5 |
| enzyme_oxa_48 | bli_sulbactam | 0.5 |
| enzyme_oxa_48 | c1_2g | 0.4 |
| enzyme_oxa_48 | c3g | 0.15 |
| enzyme_oxa_48 | c3g_bli | 0.15 |
| enzyme_oxa_48 | c4g | 0.1 |
| enzyme_oxa_48 | anti_mrsa_ceph | 0.1 |
| enzyme_oxa_48 | siderophore_ceph | 0.1 |
| enzyme_oxa_48 | cft_avi | 0.15 |
| enzyme_oxa_48 | mer_vab | 0.15 |
| enzyme_oxa_48 | azt_avi | 0.15 |
| enzyme_oxa_48 | carb_group1 | 0.7 |
| enzyme_oxa_48 | carb_group2 | 0.7 |
| enzyme_oxa_48 | fq | 0.6 |
| enzyme_oxa_48 | ag_group1 | 0.6 |
| enzyme_oxa_48 | ag_group2 | 0.6 |
| enzyme_oxa_48 | mls | 0.6 |
| enzyme_oxa_48 | lincosamides | 0.6 |
| enzyme_oxa_48 | glyc | 0.6 |
| enzyme_oxa_48 | lipoglycopeptides | 0.6 |
| enzyme_oxa_48 | tet | 0.6 |
| enzyme_oxa_48 | glycylcyclines | 0.6 |
| enzyme_oxa_48 | poly | 0.6 |
| enzyme_oxa_48 | oxa | 0.6 |
| enzyme_oxa_48 | chl | 0.6 |
| enzyme_oxa_48 | sulf | 0.6 |
| enzyme_oxa_48 | lipopeptides | 0.6 |
| enzyme_oxa_48 | streptogramins | 0.6 |
| enzyme_oxa_48 | nitrofurans | 0.6 |
| enzyme_oxa_48 | phosphonic_acids | 0.6 |
| enzyme_oxa_48 | nitroimidazoles | 0.6 |
| enzyme_oxa_48 | rifamycins | 0.6 |
| enzyme_oxa_48 | macrocycles | 0.6 |
| enzyme_oxa_48 | steroid_antibacterials | 0.6 |
| enzyme_oxa_48 | pleuromutilins | 0.6 |
| enzyme_oxa_48 | other | 0.6 |
| enzyme_ampc_cmy | pen | 0.7 |
| enzyme_ampc_cmy | bli | 0.6 |
| enzyme_ampc_cmy | bli_anti_pseudomonal | 0.6 |
| enzyme_ampc_cmy | bli_sulbactam | 0.6 |
| enzyme_ampc_cmy | c1_2g | 0.8 |
| enzyme_ampc_cmy | c3g | 0.8 |
| enzyme_ampc_cmy | c3g_bli | 0.8 |
| enzyme_ampc_cmy | c4g | 0.15 |
| enzyme_ampc_cmy | anti_mrsa_ceph | 0.15 |
| enzyme_ampc_cmy | siderophore_ceph | 0.15 |
| enzyme_ampc_cmy | cft_avi | 0.1 |
| enzyme_ampc_cmy | mer_vab | 0.1 |
| enzyme_ampc_cmy | azt_avi | 0.1 |
| enzyme_ampc_cmy | mono | 0.1 |
| enzyme_ampc_cmy | fq | 0.7 |
| enzyme_ampc_cmy | ag_group1 | 0.7 |
| enzyme_ampc_cmy | ag_group2 | 0.7 |
| enzyme_ampc_cmy | mls | 0.7 |
| enzyme_ampc_cmy | lincosamides | 0.7 |
| enzyme_ampc_cmy | glyc | 0.7 |
| enzyme_ampc_cmy | lipoglycopeptides | 0.7 |
| enzyme_ampc_cmy | tet | 0.7 |
| enzyme_ampc_cmy | glycylcyclines | 0.7 |
| enzyme_ampc_cmy | poly | 0.7 |
| enzyme_ampc_cmy | oxa | 0.7 |
| enzyme_ampc_cmy | chl | 0.7 |
| enzyme_ampc_cmy | sulf | 0.7 |
| enzyme_ampc_cmy | lipopeptides | 0.7 |
| enzyme_ampc_cmy | streptogramins | 0.7 |
| enzyme_ampc_cmy | nitrofurans | 0.7 |
| enzyme_ampc_cmy | phosphonic_acids | 0.7 |
| enzyme_ampc_cmy | nitroimidazoles | 0.7 |
| enzyme_ampc_cmy | rifamycins | 0.7 |
| enzyme_ampc_cmy | macrocycles | 0.7 |
| enzyme_ampc_cmy | steroid_antibacterials | 0.7 |
| enzyme_ampc_cmy | pleuromutilins | 0.7 |
| enzyme_ampc_cmy | other | 0.7 |
| enzyme_ampc_dha | pen | 0.7 |
| enzyme_ampc_dha | bli | 0.55 |
| enzyme_ampc_dha | bli_anti_pseudomonal | 0.55 |
| enzyme_ampc_dha | bli_sulbactam | 0.55 |
| enzyme_ampc_dha | c1_2g | 0.75 |
| enzyme_ampc_dha | c3g | 0.75 |
| enzyme_ampc_dha | c3g_bli | 0.75 |
| enzyme_ampc_dha | c4g | 0.15 |
| enzyme_ampc_dha | anti_mrsa_ceph | 0.15 |
| enzyme_ampc_dha | siderophore_ceph | 0.15 |
| enzyme_ampc_dha | cft_avi | 0.1 |
| enzyme_ampc_dha | mer_vab | 0.1 |
| enzyme_ampc_dha | azt_avi | 0.1 |
| enzyme_ampc_dha | mono | 0.1 |
| enzyme_ampc_dha | fq | 0.7 |
| enzyme_ampc_dha | ag_group1 | 0.7 |
| enzyme_ampc_dha | ag_group2 | 0.7 |
| enzyme_ampc_dha | mls | 0.7 |
| enzyme_ampc_dha | lincosamides | 0.7 |
| enzyme_ampc_dha | glyc | 0.7 |
| enzyme_ampc_dha | lipoglycopeptides | 0.7 |
| enzyme_ampc_dha | tet | 0.7 |
| enzyme_ampc_dha | glycylcyclines | 0.7 |
| enzyme_ampc_dha | poly | 0.7 |
| enzyme_ampc_dha | oxa | 0.7 |
| enzyme_ampc_dha | chl | 0.7 |
| enzyme_ampc_dha | sulf | 0.7 |
| enzyme_ampc_dha | lipopeptides | 0.7 |
| enzyme_ampc_dha | streptogramins | 0.7 |
| enzyme_ampc_dha | nitrofurans | 0.7 |
| enzyme_ampc_dha | phosphonic_acids | 0.7 |
| enzyme_ampc_dha | nitroimidazoles | 0.7 |
| enzyme_ampc_dha | rifamycins | 0.7 |
| enzyme_ampc_dha | macrocycles | 0.7 |
| enzyme_ampc_dha | steroid_antibacterials | 0.7 |
| enzyme_ampc_dha | pleuromutilins | 0.7 |
| enzyme_ampc_dha | other | 0.7 |
| target_site_pbp2a_meca | pen | 0.99 |
| target_site_pbp2a_meca | bli | 0.99 |
| target_site_pbp2a_meca | bli_anti_pseudomonal | 0.99 |
| target_site_pbp2a_meca | bli_sulbactam | 0.99 |
| target_site_pbp2a_meca | c1_2g | 0.99 |
| target_site_pbp2a_meca | c3g | 0.99 |
| target_site_pbp2a_meca | c3g_bli | 0.99 |
| target_site_pbp2a_meca | c4g | 0.7 |
| target_site_pbp2a_meca | anti_mrsa_ceph | 0.7 |
| target_site_pbp2a_meca | siderophore_ceph | 0.7 |
| target_site_pbp2a_meca | cft_avi | 0.99 |
| target_site_pbp2a_meca | mer_vab | 0.99 |
| target_site_pbp2a_meca | azt_avi | 0.99 |
| target_site_pbp2a_meca | carb_group1 | 0.85 |
| target_site_pbp2a_meca | carb_group2 | 0.85 |
| target_site_pbp2a_meca | mono | 0.99 |
| target_site_pbp2a_meca | fq | 0.99 |
| target_site_pbp2a_meca | ag_group1 | 0.99 |
| target_site_pbp2a_meca | ag_group2 | 0.99 |
| target_site_pbp2a_meca | mls | 0.99 |
| target_site_pbp2a_meca | lincosamides | 0.99 |
| target_site_pbp2a_meca | glyc | 0.99 |
| target_site_pbp2a_meca | lipoglycopeptides | 0.99 |
| target_site_pbp2a_meca | tet | 0.99 |
| target_site_pbp2a_meca | glycylcyclines | 0.99 |
| target_site_pbp2a_meca | poly | 0.99 |
| target_site_pbp2a_meca | oxa | 0.99 |
| target_site_pbp2a_meca | chl | 0.99 |
| target_site_pbp2a_meca | sulf | 0.99 |
| target_site_pbp2a_meca | lipopeptides | 0.99 |
| target_site_pbp2a_meca | streptogramins | 0.99 |
| target_site_pbp2a_meca | nitrofurans | 0.99 |
| target_site_pbp2a_meca | phosphonic_acids | 0.99 |
| target_site_pbp2a_meca | nitroimidazoles | 0.99 |
| target_site_pbp2a_meca | rifamycins | 0.99 |
| target_site_pbp2a_meca | macrocycles | 0.99 |
| target_site_pbp2a_meca | steroid_antibacterials | 0.99 |
| target_site_pbp2a_meca | pleuromutilins | 0.99 |
| target_site_pbp2a_meca | other | 0.99 |
| target_site_van_a | pen | 0.99 |
| target_site_van_a | bli | 0.99 |
| target_site_van_a | bli_anti_pseudomonal | 0.99 |
| target_site_van_a | bli_sulbactam | 0.99 |
| target_site_van_a | c1_2g | 0.99 |
| target_site_van_a | c3g | 0.99 |
| target_site_van_a | c3g_bli | 0.99 |
| target_site_van_a | c4g | 0.99 |
| target_site_van_a | anti_mrsa_ceph | 0.99 |
| target_site_van_a | siderophore_ceph | 0.99 |
| target_site_van_a | cft_avi | 0.99 |
| target_site_van_a | mer_vab | 0.99 |
| target_site_van_a | azt_avi | 0.99 |
| target_site_van_a | carb_group1 | 0.99 |
| target_site_van_a | carb_group2 | 0.99 |
| target_site_van_a | mono | 0.99 |
| target_site_van_a | fq | 0.99 |
| target_site_van_a | ag_group1 | 0.99 |
| target_site_van_a | ag_group2 | 0.99 |
| target_site_van_a | mls | 0.99 |
| target_site_van_a | lincosamides | 0.99 |
| target_site_van_a | glyc | 0.99 |
| target_site_van_a | lipoglycopeptides | 0.99 |
| target_site_van_a | tet | 0.99 |
| target_site_van_a | glycylcyclines | 0.99 |
| target_site_van_a | poly | 0.99 |
| target_site_van_a | oxa | 0.99 |
| target_site_van_a | chl | 0.99 |
| target_site_van_a | sulf | 0.99 |
| target_site_van_a | lipopeptides | 0.99 |
| target_site_van_a | streptogramins | 0.99 |
| target_site_van_a | nitrofurans | 0.99 |
| target_site_van_a | phosphonic_acids | 0.99 |
| target_site_van_a | nitroimidazoles | 0.99 |
| target_site_van_a | rifamycins | 0.99 |
| target_site_van_a | macrocycles | 0.99 |
| target_site_van_a | steroid_antibacterials | 0.99 |
| target_site_van_a | pleuromutilins | 0.99 |
| target_site_van_a | other | 0.99 |
| target_site_van_b | pen | 0.99 |
| target_site_van_b | bli | 0.99 |
| target_site_van_b | bli_anti_pseudomonal | 0.99 |
| target_site_van_b | bli_sulbactam | 0.99 |
| target_site_van_b | c1_2g | 0.99 |
| target_site_van_b | c3g | 0.99 |
| target_site_van_b | c3g_bli | 0.99 |
| target_site_van_b | c4g | 0.99 |
| target_site_van_b | anti_mrsa_ceph | 0.99 |
| target_site_van_b | siderophore_ceph | 0.99 |
| target_site_van_b | cft_avi | 0.99 |
| target_site_van_b | mer_vab | 0.99 |
| target_site_van_b | azt_avi | 0.99 |
| target_site_van_b | carb_group1 | 0.99 |
| target_site_van_b | carb_group2 | 0.99 |
| target_site_van_b | mono | 0.99 |
| target_site_van_b | fq | 0.99 |
| target_site_van_b | ag_group1 | 0.99 |
| target_site_van_b | ag_group2 | 0.99 |
| target_site_van_b | mls | 0.99 |
| target_site_van_b | lincosamides | 0.99 |
| target_site_van_b | glyc | 0.7 |
| target_site_van_b | lipoglycopeptides | 0.7 |
| target_site_van_b | tet | 0.99 |
| target_site_van_b | glycylcyclines | 0.99 |
| target_site_van_b | poly | 0.99 |
| target_site_van_b | oxa | 0.99 |
| target_site_van_b | chl | 0.99 |
| target_site_van_b | sulf | 0.99 |
| target_site_van_b | lipopeptides | 0.99 |
| target_site_van_b | streptogramins | 0.99 |
| target_site_van_b | nitrofurans | 0.99 |
| target_site_van_b | phosphonic_acids | 0.99 |
| target_site_van_b | nitroimidazoles | 0.99 |
| target_site_van_b | rifamycins | 0.99 |
| target_site_van_b | macrocycles | 0.99 |
| target_site_van_b | steroid_antibacterials | 0.99 |
| target_site_van_b | pleuromutilins | 0.99 |
| target_site_van_b | other | 0.99 |
| mutation_gyra_primary | pen | 0.4 |
| mutation_gyra_primary | bli | 0.4 |
| mutation_gyra_primary | bli_anti_pseudomonal | 0.4 |
| mutation_gyra_primary | bli_sulbactam | 0.4 |
| mutation_gyra_primary | c1_2g | 0.4 |
| mutation_gyra_primary | c3g | 0.4 |
| mutation_gyra_primary | c3g_bli | 0.4 |
| mutation_gyra_primary | c4g | 0.4 |
| mutation_gyra_primary | anti_mrsa_ceph | 0.4 |
| mutation_gyra_primary | siderophore_ceph | 0.4 |
| mutation_gyra_primary | cft_avi | 0.4 |
| mutation_gyra_primary | mer_vab | 0.4 |
| mutation_gyra_primary | azt_avi | 0.4 |
| mutation_gyra_primary | carb_group1 | 0.4 |
| mutation_gyra_primary | carb_group2 | 0.4 |
| mutation_gyra_primary | mono | 0.4 |
| mutation_gyra_primary | fq | 0.4 |
| mutation_gyra_primary | ag_group1 | 0.4 |
| mutation_gyra_primary | ag_group2 | 0.4 |
| mutation_gyra_primary | mls | 0.4 |
| mutation_gyra_primary | lincosamides | 0.4 |
| mutation_gyra_primary | glyc | 0.4 |
| mutation_gyra_primary | lipoglycopeptides | 0.4 |
| mutation_gyra_primary | tet | 0.4 |
| mutation_gyra_primary | glycylcyclines | 0.4 |
| mutation_gyra_primary | poly | 0.4 |
| mutation_gyra_primary | oxa | 0.4 |
| mutation_gyra_primary | chl | 0.4 |
| mutation_gyra_primary | sulf | 0.4 |
| mutation_gyra_primary | lipopeptides | 0.4 |
| mutation_gyra_primary | streptogramins | 0.4 |
| mutation_gyra_primary | nitrofurans | 0.4 |
| mutation_gyra_primary | phosphonic_acids | 0.4 |
| mutation_gyra_primary | nitroimidazoles | 0.4 |
| mutation_gyra_primary | rifamycins | 0.4 |
| mutation_gyra_primary | macrocycles | 0.4 |
| mutation_gyra_primary | steroid_antibacterials | 0.4 |
| mutation_gyra_primary | pleuromutilins | 0.4 |
| mutation_gyra_primary | other | 0.4 |
| mutation_gyra_parc_secondary | pen | 0.95 |
| mutation_gyra_parc_secondary | bli | 0.95 |
| mutation_gyra_parc_secondary | bli_anti_pseudomonal | 0.95 |
| mutation_gyra_parc_secondary | bli_sulbactam | 0.95 |
| mutation_gyra_parc_secondary | c1_2g | 0.95 |
| mutation_gyra_parc_secondary | c3g | 0.95 |
| mutation_gyra_parc_secondary | c3g_bli | 0.95 |
| mutation_gyra_parc_secondary | c4g | 0.95 |
| mutation_gyra_parc_secondary | anti_mrsa_ceph | 0.95 |
| mutation_gyra_parc_secondary | siderophore_ceph | 0.95 |
| mutation_gyra_parc_secondary | cft_avi | 0.95 |
| mutation_gyra_parc_secondary | mer_vab | 0.95 |
| mutation_gyra_parc_secondary | azt_avi | 0.95 |
| mutation_gyra_parc_secondary | carb_group1 | 0.95 |
| mutation_gyra_parc_secondary | carb_group2 | 0.95 |
| mutation_gyra_parc_secondary | mono | 0.95 |
| mutation_gyra_parc_secondary | fq | 0.95 |
| mutation_gyra_parc_secondary | ag_group1 | 0.95 |
| mutation_gyra_parc_secondary | ag_group2 | 0.95 |
| mutation_gyra_parc_secondary | mls | 0.95 |
| mutation_gyra_parc_secondary | lincosamides | 0.95 |
| mutation_gyra_parc_secondary | glyc | 0.95 |
| mutation_gyra_parc_secondary | lipoglycopeptides | 0.95 |
| mutation_gyra_parc_secondary | tet | 0.95 |
| mutation_gyra_parc_secondary | glycylcyclines | 0.95 |
| mutation_gyra_parc_secondary | poly | 0.95 |
| mutation_gyra_parc_secondary | oxa | 0.95 |
| mutation_gyra_parc_secondary | chl | 0.95 |
| mutation_gyra_parc_secondary | sulf | 0.95 |
| mutation_gyra_parc_secondary | lipopeptides | 0.95 |
| mutation_gyra_parc_secondary | streptogramins | 0.95 |
| mutation_gyra_parc_secondary | nitrofurans | 0.95 |
| mutation_gyra_parc_secondary | phosphonic_acids | 0.95 |
| mutation_gyra_parc_secondary | nitroimidazoles | 0.95 |
| mutation_gyra_parc_secondary | rifamycins | 0.95 |
| mutation_gyra_parc_secondary | macrocycles | 0.95 |
| mutation_gyra_parc_secondary | steroid_antibacterials | 0.95 |
| mutation_gyra_parc_secondary | pleuromutilins | 0.95 |
| mutation_gyra_parc_secondary | other | 0.95 |
| protection_qnr | pen | 0.2 |
| protection_qnr | bli | 0.2 |
| protection_qnr | bli_anti_pseudomonal | 0.2 |
| protection_qnr | bli_sulbactam | 0.2 |
| protection_qnr | c1_2g | 0.2 |
| protection_qnr | c3g | 0.2 |
| protection_qnr | c3g_bli | 0.2 |
| protection_qnr | c4g | 0.2 |
| protection_qnr | anti_mrsa_ceph | 0.2 |
| protection_qnr | siderophore_ceph | 0.2 |
| protection_qnr | cft_avi | 0.2 |
| protection_qnr | mer_vab | 0.2 |
| protection_qnr | azt_avi | 0.2 |
| protection_qnr | carb_group1 | 0.2 |
| protection_qnr | carb_group2 | 0.2 |
| protection_qnr | mono | 0.2 |
| protection_qnr | fq | 0.2 |
| protection_qnr | ag_group1 | 0.2 |
| protection_qnr | ag_group2 | 0.2 |
| protection_qnr | mls | 0.2 |
| protection_qnr | lincosamides | 0.2 |
| protection_qnr | glyc | 0.2 |
| protection_qnr | lipoglycopeptides | 0.2 |
| protection_qnr | tet | 0.2 |
| protection_qnr | glycylcyclines | 0.2 |
| protection_qnr | poly | 0.2 |
| protection_qnr | oxa | 0.2 |
| protection_qnr | chl | 0.2 |
| protection_qnr | sulf | 0.2 |
| protection_qnr | lipopeptides | 0.2 |
| protection_qnr | streptogramins | 0.2 |
| protection_qnr | nitrofurans | 0.2 |
| protection_qnr | phosphonic_acids | 0.2 |
| protection_qnr | nitroimidazoles | 0.2 |
| protection_qnr | rifamycins | 0.2 |
| protection_qnr | macrocycles | 0.2 |
| protection_qnr | steroid_antibacterials | 0.2 |
| protection_qnr | pleuromutilins | 0.2 |
| protection_qnr | other | 0.2 |
| enzyme_16s_rrmt | pen | 0.95 |
| enzyme_16s_rrmt | bli | 0.95 |
| enzyme_16s_rrmt | bli_anti_pseudomonal | 0.95 |
| enzyme_16s_rrmt | bli_sulbactam | 0.95 |
| enzyme_16s_rrmt | c1_2g | 0.95 |
| enzyme_16s_rrmt | c3g | 0.95 |
| enzyme_16s_rrmt | c3g_bli | 0.95 |
| enzyme_16s_rrmt | c4g | 0.95 |
| enzyme_16s_rrmt | anti_mrsa_ceph | 0.95 |
| enzyme_16s_rrmt | siderophore_ceph | 0.95 |
| enzyme_16s_rrmt | cft_avi | 0.95 |
| enzyme_16s_rrmt | mer_vab | 0.95 |
| enzyme_16s_rrmt | azt_avi | 0.95 |
| enzyme_16s_rrmt | carb_group1 | 0.95 |
| enzyme_16s_rrmt | carb_group2 | 0.95 |
| enzyme_16s_rrmt | mono | 0.95 |
| enzyme_16s_rrmt | fq | 0.95 |
| enzyme_16s_rrmt | ag_group1 | 0.95 |
| enzyme_16s_rrmt | ag_group2 | 0.95 |
| enzyme_16s_rrmt | mls | 0.95 |
| enzyme_16s_rrmt | lincosamides | 0.95 |
| enzyme_16s_rrmt | glyc | 0.95 |
| enzyme_16s_rrmt | lipoglycopeptides | 0.95 |
| enzyme_16s_rrmt | tet | 0.95 |
| enzyme_16s_rrmt | glycylcyclines | 0.95 |
| enzyme_16s_rrmt | poly | 0.95 |
| enzyme_16s_rrmt | oxa | 0.95 |
| enzyme_16s_rrmt | chl | 0.95 |
| enzyme_16s_rrmt | sulf | 0.95 |
| enzyme_16s_rrmt | lipopeptides | 0.95 |
| enzyme_16s_rrmt | streptogramins | 0.95 |
| enzyme_16s_rrmt | nitrofurans | 0.95 |
| enzyme_16s_rrmt | phosphonic_acids | 0.95 |
| enzyme_16s_rrmt | nitroimidazoles | 0.95 |
| enzyme_16s_rrmt | rifamycins | 0.95 |
| enzyme_16s_rrmt | macrocycles | 0.95 |
| enzyme_16s_rrmt | steroid_antibacterials | 0.95 |
| enzyme_16s_rrmt | pleuromutilins | 0.95 |
| enzyme_16s_rrmt | other | 0.95 |
| target_site_erm_b | pen | 0.9 |
| target_site_erm_b | bli | 0.9 |
| target_site_erm_b | bli_anti_pseudomonal | 0.9 |
| target_site_erm_b | bli_sulbactam | 0.9 |
| target_site_erm_b | c1_2g | 0.9 |
| target_site_erm_b | c3g | 0.9 |
| target_site_erm_b | c3g_bli | 0.9 |
| target_site_erm_b | c4g | 0.9 |
| target_site_erm_b | anti_mrsa_ceph | 0.9 |
| target_site_erm_b | siderophore_ceph | 0.9 |
| target_site_erm_b | cft_avi | 0.9 |
| target_site_erm_b | mer_vab | 0.9 |
| target_site_erm_b | azt_avi | 0.9 |
| target_site_erm_b | carb_group1 | 0.9 |
| target_site_erm_b | carb_group2 | 0.9 |
| target_site_erm_b | mono | 0.9 |
| target_site_erm_b | fq | 0.9 |
| target_site_erm_b | ag_group1 | 0.9 |
| target_site_erm_b | ag_group2 | 0.9 |
| target_site_erm_b | mls | 0.9 |
| target_site_erm_b | lincosamides | 0.9 |
| target_site_erm_b | glyc | 0.9 |
| target_site_erm_b | lipoglycopeptides | 0.9 |
| target_site_erm_b | tet | 0.9 |
| target_site_erm_b | glycylcyclines | 0.9 |
| target_site_erm_b | poly | 0.9 |
| target_site_erm_b | oxa | 0.9 |
| target_site_erm_b | chl | 0.9 |
| target_site_erm_b | sulf | 0.9 |
| target_site_erm_b | lipopeptides | 0.9 |
| target_site_erm_b | streptogramins | 0.9 |
| target_site_erm_b | nitrofurans | 0.9 |
| target_site_erm_b | phosphonic_acids | 0.9 |
| target_site_erm_b | nitroimidazoles | 0.9 |
| target_site_erm_b | rifamycins | 0.9 |
| target_site_erm_b | macrocycles | 0.9 |
| target_site_erm_b | steroid_antibacterials | 0.9 |
| target_site_erm_b | pleuromutilins | 0.9 |
| target_site_erm_b | other | 0.9 |
| target_site_cfr | pen | 0.95 |
| target_site_cfr | bli | 0.95 |
| target_site_cfr | bli_anti_pseudomonal | 0.95 |
| target_site_cfr | bli_sulbactam | 0.95 |
| target_site_cfr | c1_2g | 0.95 |
| target_site_cfr | c3g | 0.95 |
| target_site_cfr | c3g_bli | 0.95 |
| target_site_cfr | c4g | 0.95 |
| target_site_cfr | anti_mrsa_ceph | 0.95 |
| target_site_cfr | siderophore_ceph | 0.95 |
| target_site_cfr | cft_avi | 0.95 |
| target_site_cfr | mer_vab | 0.95 |
| target_site_cfr | azt_avi | 0.95 |
| target_site_cfr | carb_group1 | 0.95 |
| target_site_cfr | carb_group2 | 0.95 |
| target_site_cfr | mono | 0.95 |
| target_site_cfr | fq | 0.95 |
| target_site_cfr | ag_group1 | 0.95 |
| target_site_cfr | ag_group2 | 0.95 |
| target_site_cfr | mls | 0.7 |
| target_site_cfr | lincosamides | 0.7 |
| target_site_cfr | glyc | 0.95 |
| target_site_cfr | lipoglycopeptides | 0.95 |
| target_site_cfr | tet | 0.95 |
| target_site_cfr | glycylcyclines | 0.95 |
| target_site_cfr | poly | 0.95 |
| target_site_cfr | oxa | 0.9 |
| target_site_cfr | chl | 0.7 |
| target_site_cfr | sulf | 0.95 |
| target_site_cfr | lipopeptides | 0.95 |
| target_site_cfr | streptogramins | 0.95 |
| target_site_cfr | nitrofurans | 0.95 |
| target_site_cfr | phosphonic_acids | 0.95 |
| target_site_cfr | nitroimidazoles | 0.95 |
| target_site_cfr | rifamycins | 0.95 |
| target_site_cfr | macrocycles | 0.95 |
| target_site_cfr | steroid_antibacterials | 0.95 |
| target_site_cfr | pleuromutilins | 0.95 |
| target_site_cfr | other | 0.95 |
| enzyme_cat | pen | 0.9 |
| enzyme_cat | bli | 0.9 |
| enzyme_cat | bli_anti_pseudomonal | 0.9 |
| enzyme_cat | bli_sulbactam | 0.9 |
| enzyme_cat | c1_2g | 0.9 |
| enzyme_cat | c3g | 0.9 |
| enzyme_cat | c3g_bli | 0.9 |
| enzyme_cat | c4g | 0.9 |
| enzyme_cat | anti_mrsa_ceph | 0.9 |
| enzyme_cat | siderophore_ceph | 0.9 |
| enzyme_cat | cft_avi | 0.9 |
| enzyme_cat | mer_vab | 0.9 |
| enzyme_cat | azt_avi | 0.9 |
| enzyme_cat | carb_group1 | 0.9 |
| enzyme_cat | carb_group2 | 0.9 |
| enzyme_cat | mono | 0.9 |
| enzyme_cat | fq | 0.9 |
| enzyme_cat | ag_group1 | 0.9 |
| enzyme_cat | ag_group2 | 0.9 |
| enzyme_cat | mls | 0.9 |
| enzyme_cat | lincosamides | 0.9 |
| enzyme_cat | glyc | 0.9 |
| enzyme_cat | lipoglycopeptides | 0.9 |
| enzyme_cat | tet | 0.9 |
| enzyme_cat | glycylcyclines | 0.9 |
| enzyme_cat | poly | 0.9 |
| enzyme_cat | oxa | 0.9 |
| enzyme_cat | chl | 0.9 |
| enzyme_cat | sulf | 0.9 |
| enzyme_cat | lipopeptides | 0.9 |
| enzyme_cat | streptogramins | 0.9 |
| enzyme_cat | nitrofurans | 0.9 |
| enzyme_cat | phosphonic_acids | 0.9 |
| enzyme_cat | nitroimidazoles | 0.9 |
| enzyme_cat | rifamycins | 0.9 |
| enzyme_cat | macrocycles | 0.9 |
| enzyme_cat | steroid_antibacterials | 0.9 |
| enzyme_cat | pleuromutilins | 0.9 |
| enzyme_cat | other | 0.9 |
| efflux_acrab_tolc | c1_2g | 0.3 |
| efflux_acrab_tolc | cft_avi | 0.3 |
| efflux_acrab_tolc | mer_vab | 0.3 |
| efflux_acrab_tolc | azt_avi | 0.3 |
| efflux_acrab_tolc | carb_group1 | 0.3 |
| efflux_acrab_tolc | carb_group2 | 0.3 |
| efflux_acrab_tolc | mono | 0.3 |
| efflux_acrab_tolc | fq | 0.25 |
| efflux_acrab_tolc | glyc | 0.3 |
| efflux_acrab_tolc | lipoglycopeptides | 0.3 |
| efflux_acrab_tolc | tet | 0.25 |
| efflux_acrab_tolc | glycylcyclines | 0.25 |
| efflux_acrab_tolc | poly | 0.3 |
| efflux_acrab_tolc | oxa | 0.3 |
| efflux_acrab_tolc | chl | 0.2 |
| efflux_acrab_tolc | sulf | 0.3 |
| efflux_acrab_tolc | lipopeptides | 0.3 |
| efflux_acrab_tolc | streptogramins | 0.3 |
| efflux_acrab_tolc | nitrofurans | 0.3 |
| efflux_acrab_tolc | phosphonic_acids | 0.3 |
| efflux_acrab_tolc | nitroimidazoles | 0.3 |
| efflux_acrab_tolc | rifamycins | 0.3 |
| efflux_acrab_tolc | macrocycles | 0.3 |
| efflux_acrab_tolc | steroid_antibacterials | 0.3 |
| efflux_acrab_tolc | pleuromutilins | 0.3 |
| efflux_acrab_tolc | other | 0.3 |
| efflux_mexxy_oprm | pen | 0.3 |
| efflux_mexxy_oprm | bli | 0.3 |
| efflux_mexxy_oprm | bli_anti_pseudomonal | 0.3 |
| efflux_mexxy_oprm | bli_sulbactam | 0.3 |
| efflux_mexxy_oprm | c1_2g | 0.3 |
| efflux_mexxy_oprm | c3g | 0.3 |
| efflux_mexxy_oprm | c3g_bli | 0.3 |
| efflux_mexxy_oprm | c4g | 0.2 |
| efflux_mexxy_oprm | anti_mrsa_ceph | 0.2 |
| efflux_mexxy_oprm | siderophore_ceph | 0.2 |
| efflux_mexxy_oprm | cft_avi | 0.3 |
| efflux_mexxy_oprm | mer_vab | 0.3 |
| efflux_mexxy_oprm | azt_avi | 0.3 |
| efflux_mexxy_oprm | carb_group1 | 0.05 |
| efflux_mexxy_oprm | carb_group2 | 0.05 |
| efflux_mexxy_oprm | mono | 0.3 |
| efflux_mexxy_oprm | fq | 0.2 |
| efflux_mexxy_oprm | ag_group1 | 0.3 |
| efflux_mexxy_oprm | ag_group2 | 0.3 |
| efflux_mexxy_oprm | mls | 0.3 |
| efflux_mexxy_oprm | lincosamides | 0.3 |
| efflux_mexxy_oprm | glyc | 0.3 |
| efflux_mexxy_oprm | lipoglycopeptides | 0.3 |
| efflux_mexxy_oprm | tet | 0.3 |
| efflux_mexxy_oprm | glycylcyclines | 0.3 |
| efflux_mexxy_oprm | poly | 0.3 |
| efflux_mexxy_oprm | oxa | 0.3 |
| efflux_mexxy_oprm | chl | 0.3 |
| efflux_mexxy_oprm | sulf | 0.3 |
| efflux_mexxy_oprm | lipopeptides | 0.3 |
| efflux_mexxy_oprm | streptogramins | 0.3 |
| efflux_mexxy_oprm | nitrofurans | 0.3 |
| efflux_mexxy_oprm | phosphonic_acids | 0.3 |
| efflux_mexxy_oprm | nitroimidazoles | 0.3 |
| efflux_mexxy_oprm | rifamycins | 0.3 |
| efflux_mexxy_oprm | macrocycles | 0.3 |
| efflux_mexxy_oprm | steroid_antibacterials | 0.3 |
| efflux_mexxy_oprm | pleuromutilins | 0.3 |
| efflux_mexxy_oprm | other | 0.3 |
| porin_loss_ompk35_36 | pen | 0.3 |
| porin_loss_ompk35_36 | bli | 0.4 |
| porin_loss_ompk35_36 | bli_anti_pseudomonal | 0.4 |
| porin_loss_ompk35_36 | bli_sulbactam | 0.4 |
| porin_loss_ompk35_36 | c1_2g | 0.8 |
| porin_loss_ompk35_36 | c3g | 0.4 |
| porin_loss_ompk35_36 | c3g_bli | 0.4 |
| porin_loss_ompk35_36 | c4g | 0.3 |
| porin_loss_ompk35_36 | anti_mrsa_ceph | 0.3 |
| porin_loss_ompk35_36 | siderophore_ceph | 0.3 |
| porin_loss_ompk35_36 | cft_avi | 0.25 |
| porin_loss_ompk35_36 | mer_vab | 0.25 |
| porin_loss_ompk35_36 | azt_avi | 0.25 |
| porin_loss_ompk35_36 | mono | 0.8 |
| porin_loss_ompk35_36 | mls | 0.8 |
| porin_loss_ompk35_36 | lincosamides | 0.8 |
| porin_loss_ompk35_36 | glyc | 0.8 |
| porin_loss_ompk35_36 | lipoglycopeptides | 0.8 |
| porin_loss_ompk35_36 | tet | 0.8 |
| porin_loss_ompk35_36 | glycylcyclines | 0.8 |
| porin_loss_ompk35_36 | poly | 0.8 |
| porin_loss_ompk35_36 | oxa | 0.8 |
| porin_loss_ompk35_36 | chl | 0.8 |
| porin_loss_ompk35_36 | sulf | 0.8 |
| porin_loss_ompk35_36 | lipopeptides | 0.8 |
| porin_loss_ompk35_36 | streptogramins | 0.8 |
| porin_loss_ompk35_36 | nitrofurans | 0.8 |
| porin_loss_ompk35_36 | phosphonic_acids | 0.8 |
| porin_loss_ompk35_36 | nitroimidazoles | 0.8 |
| porin_loss_ompk35_36 | rifamycins | 0.8 |
| porin_loss_ompk35_36 | macrocycles | 0.8 |
| porin_loss_ompk35_36 | steroid_antibacterials | 0.8 |
| porin_loss_ompk35_36 | pleuromutilins | 0.8 |
| porin_loss_ompk35_36 | other | 0.8 |
| porin_loss_oprd | pen | 0.8 |
| porin_loss_oprd | bli | 0.8 |
| porin_loss_oprd | bli_anti_pseudomonal | 0.8 |
| porin_loss_oprd | bli_sulbactam | 0.8 |
| porin_loss_oprd | c1_2g | 0.8 |
| porin_loss_oprd | c3g | 0.8 |
| porin_loss_oprd | c3g_bli | 0.8 |
| porin_loss_oprd | c4g | 0.8 |
| porin_loss_oprd | anti_mrsa_ceph | 0.8 |
| porin_loss_oprd | siderophore_ceph | 0.8 |
| porin_loss_oprd | cft_avi | 0.8 |
| porin_loss_oprd | mer_vab | 0.8 |
| porin_loss_oprd | azt_avi | 0.8 |
| porin_loss_oprd | carb_group1 | 0.8 |
| porin_loss_oprd | carb_group2 | 0.8 |
| porin_loss_oprd | mono | 0.8 |
| porin_loss_oprd | fq | 0.8 |
| porin_loss_oprd | ag_group1 | 0.8 |
| porin_loss_oprd | ag_group2 | 0.8 |
| porin_loss_oprd | mls | 0.8 |
| porin_loss_oprd | lincosamides | 0.8 |
| porin_loss_oprd | glyc | 0.8 |
| porin_loss_oprd | lipoglycopeptides | 0.8 |
| porin_loss_oprd | tet | 0.8 |
| porin_loss_oprd | glycylcyclines | 0.8 |
| porin_loss_oprd | poly | 0.8 |
| porin_loss_oprd | oxa | 0.8 |
| porin_loss_oprd | chl | 0.8 |
| porin_loss_oprd | sulf | 0.8 |
| porin_loss_oprd | lipopeptides | 0.8 |
| porin_loss_oprd | streptogramins | 0.8 |
| porin_loss_oprd | nitrofurans | 0.8 |
| porin_loss_oprd | phosphonic_acids | 0.8 |
| porin_loss_oprd | nitroimidazoles | 0.8 |
| porin_loss_oprd | rifamycins | 0.8 |
| porin_loss_oprd | macrocycles | 0.8 |
| porin_loss_oprd | steroid_antibacterials | 0.8 |
| porin_loss_oprd | pleuromutilins | 0.8 |
| porin_loss_oprd | other | 0.8 |
| modification_mcr_1 | pen | 0.85 |
| modification_mcr_1 | bli | 0.85 |
| modification_mcr_1 | bli_anti_pseudomonal | 0.85 |
| modification_mcr_1 | bli_sulbactam | 0.85 |
| modification_mcr_1 | c1_2g | 0.85 |
| modification_mcr_1 | c3g | 0.85 |
| modification_mcr_1 | c3g_bli | 0.85 |
| modification_mcr_1 | c4g | 0.85 |
| modification_mcr_1 | anti_mrsa_ceph | 0.85 |
| modification_mcr_1 | siderophore_ceph | 0.85 |
| modification_mcr_1 | cft_avi | 0.85 |
| modification_mcr_1 | mer_vab | 0.85 |
| modification_mcr_1 | azt_avi | 0.85 |
| modification_mcr_1 | carb_group1 | 0.85 |
| modification_mcr_1 | carb_group2 | 0.85 |
| modification_mcr_1 | mono | 0.85 |
| modification_mcr_1 | fq | 0.85 |
| modification_mcr_1 | ag_group1 | 0.85 |
| modification_mcr_1 | ag_group2 | 0.85 |
| modification_mcr_1 | mls | 0.85 |
| modification_mcr_1 | lincosamides | 0.85 |
| modification_mcr_1 | glyc | 0.85 |
| modification_mcr_1 | lipoglycopeptides | 0.85 |
| modification_mcr_1 | tet | 0.85 |
| modification_mcr_1 | glycylcyclines | 0.85 |
| modification_mcr_1 | poly | 0.85 |
| modification_mcr_1 | oxa | 0.85 |
| modification_mcr_1 | chl | 0.85 |
| modification_mcr_1 | sulf | 0.85 |
| modification_mcr_1 | lipopeptides | 0.85 |
| modification_mcr_1 | streptogramins | 0.85 |
| modification_mcr_1 | nitrofurans | 0.85 |
| modification_mcr_1 | phosphonic_acids | 0.85 |
| modification_mcr_1 | nitroimidazoles | 0.85 |
| modification_mcr_1 | rifamycins | 0.85 |
| modification_mcr_1 | macrocycles | 0.85 |
| modification_mcr_1 | steroid_antibacterials | 0.85 |
| modification_mcr_1 | pleuromutilins | 0.85 |
| modification_mcr_1 | other | 0.85 |
| global_efflux_pump | c1_2g | 0.2 |
| global_efflux_pump | cft_avi | 0.2 |
| global_efflux_pump | mer_vab | 0.2 |
| global_efflux_pump | azt_avi | 0.2 |
| global_efflux_pump | carb_group1 | 0.2 |
| global_efflux_pump | carb_group2 | 0.2 |
| global_efflux_pump | mono | 0.2 |
| global_efflux_pump | fq | 0.15 |
| global_efflux_pump | mls | 0.1 |
| global_efflux_pump | lincosamides | 0.1 |
| global_efflux_pump | glyc | 0.2 |
| global_efflux_pump | lipoglycopeptides | 0.2 |
| global_efflux_pump | tet | 0.15 |
| global_efflux_pump | glycylcyclines | 0.15 |
| global_efflux_pump | poly | 0.2 |
| global_efflux_pump | oxa | 0.2 |
| global_efflux_pump | chl | 0.15 |
| global_efflux_pump | sulf | 0.2 |
| global_efflux_pump | lipopeptides | 0.2 |
| global_efflux_pump | streptogramins | 0.2 |
| global_efflux_pump | nitrofurans | 0.2 |
| global_efflux_pump | phosphonic_acids | 0.2 |
| global_efflux_pump | nitroimidazoles | 0.2 |
| global_efflux_pump | rifamycins | 0.2 |
| global_efflux_pump | macrocycles | 0.2 |
| global_efflux_pump | steroid_antibacterials | 0.2 |
| global_efflux_pump | pleuromutilins | 0.2 |
| global_efflux_pump | other | 0.2 |
| global_porin_loss | c1_2g | 0.2 |
| global_porin_loss | cft_avi | 0.2 |
| global_porin_loss | mer_vab | 0.2 |
| global_porin_loss | azt_avi | 0.2 |
| global_porin_loss | mono | 0.2 |
| global_porin_loss | mls | 0.2 |
| global_porin_loss | lincosamides | 0.2 |
| global_porin_loss | glyc | 0.2 |
| global_porin_loss | lipoglycopeptides | 0.2 |
| global_porin_loss | tet | 0.2 |
| global_porin_loss | glycylcyclines | 0.2 |
| global_porin_loss | poly | 0.2 |
| global_porin_loss | oxa | 0.2 |
| global_porin_loss | chl | 0.2 |
| global_porin_loss | sulf | 0.2 |
| global_porin_loss | lipopeptides | 0.2 |
| global_porin_loss | streptogramins | 0.2 |
| global_porin_loss | nitrofurans | 0.2 |
| global_porin_loss | phosphonic_acids | 0.2 |
| global_porin_loss | nitroimidazoles | 0.2 |
| global_porin_loss | rifamycins | 0.2 |
| global_porin_loss | macrocycles | 0.2 |
| global_porin_loss | steroid_antibacterials | 0.2 |
| global_porin_loss | pleuromutilins | 0.2 |
| global_porin_loss | other | 0.2 |
| mutation_folate_pathway | pen | 0.85 |
| mutation_folate_pathway | bli | 0.85 |
| mutation_folate_pathway | bli_anti_pseudomonal | 0.85 |
| mutation_folate_pathway | bli_sulbactam | 0.85 |
| mutation_folate_pathway | c1_2g | 0.85 |
| mutation_folate_pathway | c3g | 0.85 |
| mutation_folate_pathway | c3g_bli | 0.85 |
| mutation_folate_pathway | c4g | 0.85 |
| mutation_folate_pathway | anti_mrsa_ceph | 0.85 |
| mutation_folate_pathway | siderophore_ceph | 0.85 |
| mutation_folate_pathway | cft_avi | 0.85 |
| mutation_folate_pathway | mer_vab | 0.85 |
| mutation_folate_pathway | azt_avi | 0.85 |
| mutation_folate_pathway | carb_group1 | 0.85 |
| mutation_folate_pathway | carb_group2 | 0.85 |
| mutation_folate_pathway | mono | 0.85 |
| mutation_folate_pathway | fq | 0.85 |
| mutation_folate_pathway | ag_group1 | 0.85 |
| mutation_folate_pathway | ag_group2 | 0.85 |
| mutation_folate_pathway | mls | 0.85 |
| mutation_folate_pathway | lincosamides | 0.85 |
| mutation_folate_pathway | glyc | 0.85 |
| mutation_folate_pathway | lipoglycopeptides | 0.85 |
| mutation_folate_pathway | tet | 0.85 |
| mutation_folate_pathway | glycylcyclines | 0.85 |
| mutation_folate_pathway | poly | 0.85 |
| mutation_folate_pathway | oxa | 0.85 |
| mutation_folate_pathway | chl | 0.85 |
| mutation_folate_pathway | sulf | 0.85 |
| mutation_folate_pathway | lipopeptides | 0.85 |
| mutation_folate_pathway | streptogramins | 0.85 |
| mutation_folate_pathway | nitrofurans | 0.85 |
| mutation_folate_pathway | phosphonic_acids | 0.85 |
| mutation_folate_pathway | nitroimidazoles | 0.85 |
| mutation_folate_pathway | rifamycins | 0.85 |
| mutation_folate_pathway | macrocycles | 0.85 |
| mutation_folate_pathway | steroid_antibacterials | 0.85 |
| mutation_folate_pathway | pleuromutilins | 0.85 |
| mutation_folate_pathway | other | 0.85 |
| mutation_nitroreductase | pen | 0.7 |
| mutation_nitroreductase | bli | 0.7 |
| mutation_nitroreductase | bli_anti_pseudomonal | 0.7 |
| mutation_nitroreductase | bli_sulbactam | 0.7 |
| mutation_nitroreductase | c1_2g | 0.7 |
| mutation_nitroreductase | c3g | 0.7 |
| mutation_nitroreductase | c3g_bli | 0.7 |
| mutation_nitroreductase | c4g | 0.7 |
| mutation_nitroreductase | anti_mrsa_ceph | 0.7 |
| mutation_nitroreductase | siderophore_ceph | 0.7 |
| mutation_nitroreductase | cft_avi | 0.7 |
| mutation_nitroreductase | mer_vab | 0.7 |
| mutation_nitroreductase | azt_avi | 0.7 |
| mutation_nitroreductase | carb_group1 | 0.7 |
| mutation_nitroreductase | carb_group2 | 0.7 |
| mutation_nitroreductase | mono | 0.7 |
| mutation_nitroreductase | fq | 0.7 |
| mutation_nitroreductase | ag_group1 | 0.7 |
| mutation_nitroreductase | ag_group2 | 0.7 |
| mutation_nitroreductase | mls | 0.7 |
| mutation_nitroreductase | lincosamides | 0.7 |
| mutation_nitroreductase | glyc | 0.7 |
| mutation_nitroreductase | lipoglycopeptides | 0.7 |
| mutation_nitroreductase | tet | 0.7 |
| mutation_nitroreductase | glycylcyclines | 0.7 |
| mutation_nitroreductase | poly | 0.7 |
| mutation_nitroreductase | oxa | 0.7 |
| mutation_nitroreductase | chl | 0.7 |
| mutation_nitroreductase | sulf | 0.7 |
| mutation_nitroreductase | lipopeptides | 0.7 |
| mutation_nitroreductase | streptogramins | 0.7 |
| mutation_nitroreductase | nitrofurans | 0.7 |
| mutation_nitroreductase | phosphonic_acids | 0.7 |
| mutation_nitroreductase | nitroimidazoles | 0.7 |
| mutation_nitroreductase | rifamycins | 0.7 |
| mutation_nitroreductase | macrocycles | 0.7 |
| mutation_nitroreductase | steroid_antibacterials | 0.7 |
| mutation_nitroreductase | pleuromutilins | 0.7 |
| mutation_nitroreductase | other | 0.7 |
| enzyme_fos_a | pen | 0.8 |
| enzyme_fos_a | bli | 0.8 |
| enzyme_fos_a | bli_anti_pseudomonal | 0.8 |
| enzyme_fos_a | bli_sulbactam | 0.8 |
| enzyme_fos_a | c1_2g | 0.8 |
| enzyme_fos_a | c3g | 0.8 |
| enzyme_fos_a | c3g_bli | 0.8 |
| enzyme_fos_a | c4g | 0.8 |
| enzyme_fos_a | anti_mrsa_ceph | 0.8 |
| enzyme_fos_a | siderophore_ceph | 0.8 |
| enzyme_fos_a | cft_avi | 0.8 |
| enzyme_fos_a | mer_vab | 0.8 |
| enzyme_fos_a | azt_avi | 0.8 |
| enzyme_fos_a | carb_group1 | 0.8 |
| enzyme_fos_a | carb_group2 | 0.8 |
| enzyme_fos_a | mono | 0.8 |
| enzyme_fos_a | fq | 0.8 |
| enzyme_fos_a | ag_group1 | 0.8 |
| enzyme_fos_a | ag_group2 | 0.8 |
| enzyme_fos_a | mls | 0.8 |
| enzyme_fos_a | lincosamides | 0.8 |
| enzyme_fos_a | glyc | 0.8 |
| enzyme_fos_a | lipoglycopeptides | 0.8 |
| enzyme_fos_a | tet | 0.8 |
| enzyme_fos_a | glycylcyclines | 0.8 |
| enzyme_fos_a | poly | 0.8 |
| enzyme_fos_a | oxa | 0.8 |
| enzyme_fos_a | chl | 0.8 |
| enzyme_fos_a | sulf | 0.8 |
| enzyme_fos_a | lipopeptides | 0.8 |
| enzyme_fos_a | streptogramins | 0.8 |
| enzyme_fos_a | nitrofurans | 0.8 |
| enzyme_fos_a | phosphonic_acids | 0.8 |
| enzyme_fos_a | nitroimidazoles | 0.8 |
| enzyme_fos_a | rifamycins | 0.8 |
| enzyme_fos_a | macrocycles | 0.8 |
| enzyme_fos_a | steroid_antibacterials | 0.8 |
| enzyme_fos_a | pleuromutilins | 0.8 |
| enzyme_fos_a | other | 0.8 |
| mutation_mpr_f | pen | 0.6 |
| mutation_mpr_f | bli | 0.6 |
| mutation_mpr_f | bli_anti_pseudomonal | 0.6 |
| mutation_mpr_f | bli_sulbactam | 0.6 |
| mutation_mpr_f | c1_2g | 0.6 |
| mutation_mpr_f | c3g | 0.6 |
| mutation_mpr_f | c3g_bli | 0.6 |
| mutation_mpr_f | c4g | 0.6 |
| mutation_mpr_f | anti_mrsa_ceph | 0.6 |
| mutation_mpr_f | siderophore_ceph | 0.6 |
| mutation_mpr_f | cft_avi | 0.6 |
| mutation_mpr_f | mer_vab | 0.6 |
| mutation_mpr_f | azt_avi | 0.6 |
| mutation_mpr_f | carb_group1 | 0.6 |
| mutation_mpr_f | carb_group2 | 0.6 |
| mutation_mpr_f | mono | 0.6 |
| mutation_mpr_f | fq | 0.6 |
| mutation_mpr_f | ag_group1 | 0.6 |
| mutation_mpr_f | ag_group2 | 0.6 |
| mutation_mpr_f | mls | 0.6 |
| mutation_mpr_f | lincosamides | 0.6 |
| mutation_mpr_f | glyc | 0.6 |
| mutation_mpr_f | lipoglycopeptides | 0.6 |
| mutation_mpr_f | tet | 0.6 |
| mutation_mpr_f | glycylcyclines | 0.6 |
| mutation_mpr_f | poly | 0.6 |
| mutation_mpr_f | oxa | 0.6 |
| mutation_mpr_f | chl | 0.6 |
| mutation_mpr_f | sulf | 0.6 |
| mutation_mpr_f | lipopeptides | 0.6 |
| mutation_mpr_f | streptogramins | 0.6 |
| mutation_mpr_f | nitrofurans | 0.6 |
| mutation_mpr_f | phosphonic_acids | 0.6 |
| mutation_mpr_f | nitroimidazoles | 0.6 |
| mutation_mpr_f | rifamycins | 0.6 |
| mutation_mpr_f | macrocycles | 0.6 |
| mutation_mpr_f | steroid_antibacterials | 0.6 |
| mutation_mpr_f | pleuromutilins | 0.6 |
| mutation_mpr_f | other | 0.6 |
| mutation_rpo_b | pen | 0.95 |
| mutation_rpo_b | bli | 0.95 |
| mutation_rpo_b | bli_anti_pseudomonal | 0.95 |
| mutation_rpo_b | bli_sulbactam | 0.95 |
| mutation_rpo_b | c1_2g | 0.95 |
| mutation_rpo_b | c3g | 0.95 |
| mutation_rpo_b | c3g_bli | 0.95 |
| mutation_rpo_b | c4g | 0.95 |
| mutation_rpo_b | anti_mrsa_ceph | 0.95 |
| mutation_rpo_b | siderophore_ceph | 0.95 |
| mutation_rpo_b | cft_avi | 0.95 |
| mutation_rpo_b | mer_vab | 0.95 |
| mutation_rpo_b | azt_avi | 0.95 |
| mutation_rpo_b | carb_group1 | 0.95 |
| mutation_rpo_b | carb_group2 | 0.95 |
| mutation_rpo_b | mono | 0.95 |
| mutation_rpo_b | fq | 0.95 |
| mutation_rpo_b | ag_group1 | 0.95 |
| mutation_rpo_b | ag_group2 | 0.95 |
| mutation_rpo_b | mls | 0.95 |
| mutation_rpo_b | lincosamides | 0.95 |
| mutation_rpo_b | glyc | 0.95 |
| mutation_rpo_b | lipoglycopeptides | 0.95 |
| mutation_rpo_b | tet | 0.95 |
| mutation_rpo_b | glycylcyclines | 0.95 |
| mutation_rpo_b | poly | 0.95 |
| mutation_rpo_b | oxa | 0.95 |
| mutation_rpo_b | chl | 0.95 |
| mutation_rpo_b | sulf | 0.95 |
| mutation_rpo_b | lipopeptides | 0.95 |
| mutation_rpo_b | streptogramins | 0.95 |
| mutation_rpo_b | nitrofurans | 0.95 |
| mutation_rpo_b | phosphonic_acids | 0.95 |
| mutation_rpo_b | nitroimidazoles | 0.95 |
| mutation_rpo_b | rifamycins | 0.95 |
| mutation_rpo_b | macrocycles | 0.95 |
| mutation_rpo_b | steroid_antibacterials | 0.95 |
| mutation_rpo_b | pleuromutilins | 0.95 |
| mutation_rpo_b | other | 0.95 |
| protection_fus_b | pen | 0.7 |
| protection_fus_b | bli | 0.7 |
| protection_fus_b | bli_anti_pseudomonal | 0.7 |
| protection_fus_b | bli_sulbactam | 0.7 |
| protection_fus_b | c1_2g | 0.7 |
| protection_fus_b | c3g | 0.7 |
| protection_fus_b | c3g_bli | 0.7 |
| protection_fus_b | c4g | 0.7 |
| protection_fus_b | anti_mrsa_ceph | 0.7 |
| protection_fus_b | siderophore_ceph | 0.7 |
| protection_fus_b | cft_avi | 0.7 |
| protection_fus_b | mer_vab | 0.7 |
| protection_fus_b | azt_avi | 0.7 |
| protection_fus_b | carb_group1 | 0.7 |
| protection_fus_b | carb_group2 | 0.7 |
| protection_fus_b | mono | 0.7 |
| protection_fus_b | fq | 0.7 |
| protection_fus_b | ag_group1 | 0.7 |
| protection_fus_b | ag_group2 | 0.7 |
| protection_fus_b | mls | 0.7 |
| protection_fus_b | lincosamides | 0.7 |
| protection_fus_b | glyc | 0.7 |
| protection_fus_b | lipoglycopeptides | 0.7 |
| protection_fus_b | tet | 0.7 |
| protection_fus_b | glycylcyclines | 0.7 |
| protection_fus_b | poly | 0.7 |
| protection_fus_b | oxa | 0.7 |
| protection_fus_b | chl | 0.7 |
| protection_fus_b | sulf | 0.7 |
| protection_fus_b | lipopeptides | 0.7 |
| protection_fus_b | streptogramins | 0.7 |
| protection_fus_b | nitrofurans | 0.7 |
| protection_fus_b | phosphonic_acids | 0.7 |
| protection_fus_b | nitroimidazoles | 0.7 |
| protection_fus_b | rifamycins | 0.7 |
| protection_fus_b | macrocycles | 0.7 |
| protection_fus_b | steroid_antibacterials | 0.7 |
| protection_fus_b | pleuromutilins | 0.7 |
| protection_fus_b | other | 0.7 |
| protection_tet_m | pen | 0.9 |
| protection_tet_m | bli | 0.9 |
| protection_tet_m | bli_anti_pseudomonal | 0.9 |
| protection_tet_m | bli_sulbactam | 0.9 |
| protection_tet_m | c1_2g | 0.9 |
| protection_tet_m | c3g | 0.9 |
| protection_tet_m | c3g_bli | 0.9 |
| protection_tet_m | c4g | 0.9 |
| protection_tet_m | anti_mrsa_ceph | 0.9 |
| protection_tet_m | siderophore_ceph | 0.9 |
| protection_tet_m | cft_avi | 0.9 |
| protection_tet_m | mer_vab | 0.9 |
| protection_tet_m | azt_avi | 0.9 |
| protection_tet_m | carb_group1 | 0.9 |
| protection_tet_m | carb_group2 | 0.9 |
| protection_tet_m | mono | 0.9 |
| protection_tet_m | fq | 0.9 |
| protection_tet_m | ag_group1 | 0.9 |
| protection_tet_m | ag_group2 | 0.9 |
| protection_tet_m | mls | 0.9 |
| protection_tet_m | lincosamides | 0.9 |
| protection_tet_m | glyc | 0.9 |
| protection_tet_m | lipoglycopeptides | 0.9 |
| protection_tet_m | tet | 0.9 |
| protection_tet_m | glycylcyclines | 0.9 |
| protection_tet_m | poly | 0.9 |
| protection_tet_m | oxa | 0.9 |
| protection_tet_m | chl | 0.9 |
| protection_tet_m | sulf | 0.9 |
| protection_tet_m | lipopeptides | 0.9 |
| protection_tet_m | streptogramins | 0.9 |
| protection_tet_m | nitrofurans | 0.9 |
| protection_tet_m | phosphonic_acids | 0.9 |
| protection_tet_m | nitroimidazoles | 0.9 |
| protection_tet_m | rifamycins | 0.9 |
| protection_tet_m | macrocycles | 0.9 |
| protection_tet_m | steroid_antibacterials | 0.9 |
| protection_tet_m | pleuromutilins | 0.9 |
| protection_tet_m | other | 0.9 |
| enzyme_aac_aph | pen | 0.85 |
| enzyme_aac_aph | bli | 0.85 |
| enzyme_aac_aph | bli_anti_pseudomonal | 0.85 |
| enzyme_aac_aph | bli_sulbactam | 0.85 |
| enzyme_aac_aph | c1_2g | 0.85 |
| enzyme_aac_aph | c3g | 0.85 |
| enzyme_aac_aph | c3g_bli | 0.85 |
| enzyme_aac_aph | c4g | 0.85 |
| enzyme_aac_aph | anti_mrsa_ceph | 0.85 |
| enzyme_aac_aph | siderophore_ceph | 0.85 |
| enzyme_aac_aph | cft_avi | 0.85 |
| enzyme_aac_aph | mer_vab | 0.85 |
| enzyme_aac_aph | azt_avi | 0.85 |
| enzyme_aac_aph | carb_group1 | 0.85 |
| enzyme_aac_aph | carb_group2 | 0.85 |
| enzyme_aac_aph | mono | 0.85 |
| enzyme_aac_aph | fq | 0.85 |
| enzyme_aac_aph | ag_group1 | 0.85 |
| enzyme_aac_aph | ag_group2 | 0.85 |
| enzyme_aac_aph | mls | 0.85 |
| enzyme_aac_aph | lincosamides | 0.85 |
| enzyme_aac_aph | glyc | 0.85 |
| enzyme_aac_aph | lipoglycopeptides | 0.85 |
| enzyme_aac_aph | tet | 0.85 |
| enzyme_aac_aph | glycylcyclines | 0.85 |
| enzyme_aac_aph | poly | 0.85 |
| enzyme_aac_aph | oxa | 0.85 |
| enzyme_aac_aph | chl | 0.85 |
| enzyme_aac_aph | sulf | 0.85 |
| enzyme_aac_aph | lipopeptides | 0.85 |
| enzyme_aac_aph | streptogramins | 0.85 |
| enzyme_aac_aph | nitrofurans | 0.85 |
| enzyme_aac_aph | phosphonic_acids | 0.85 |
| enzyme_aac_aph | nitroimidazoles | 0.85 |
| enzyme_aac_aph | rifamycins | 0.85 |
| enzyme_aac_aph | macrocycles | 0.85 |
| enzyme_aac_aph | steroid_antibacterials | 0.85 |
| enzyme_aac_aph | pleuromutilins | 0.85 |
| enzyme_aac_aph | other | 0.85 |
| enzyme_bla_z | pen | 0.9 |
| enzyme_bla_z | bli | 0.9 |
| enzyme_bla_z | bli_anti_pseudomonal | 0.9 |
| enzyme_bla_z | bli_sulbactam | 0.9 |
| enzyme_bla_z | c1_2g | 0.9 |
| enzyme_bla_z | c3g | 0.9 |
| enzyme_bla_z | c3g_bli | 0.9 |
| enzyme_bla_z | c4g | 0.9 |
| enzyme_bla_z | anti_mrsa_ceph | 0.9 |
| enzyme_bla_z | siderophore_ceph | 0.9 |
| enzyme_bla_z | cft_avi | 0.9 |
| enzyme_bla_z | mer_vab | 0.9 |
| enzyme_bla_z | azt_avi | 0.9 |
| enzyme_bla_z | carb_group1 | 0.9 |
| enzyme_bla_z | carb_group2 | 0.9 |
| enzyme_bla_z | mono | 0.9 |
| enzyme_bla_z | fq | 0.9 |
| enzyme_bla_z | ag_group1 | 0.9 |
| enzyme_bla_z | ag_group2 | 0.9 |
| enzyme_bla_z | mls | 0.9 |
| enzyme_bla_z | lincosamides | 0.9 |
| enzyme_bla_z | glyc | 0.9 |
| enzyme_bla_z | lipoglycopeptides | 0.9 |
| enzyme_bla_z | tet | 0.9 |
| enzyme_bla_z | glycylcyclines | 0.9 |
| enzyme_bla_z | poly | 0.9 |
| enzyme_bla_z | oxa | 0.9 |
| enzyme_bla_z | chl | 0.9 |
| enzyme_bla_z | sulf | 0.9 |
| enzyme_bla_z | lipopeptides | 0.9 |
| enzyme_bla_z | streptogramins | 0.9 |
| enzyme_bla_z | nitrofurans | 0.9 |
| enzyme_bla_z | phosphonic_acids | 0.9 |
| enzyme_bla_z | nitroimidazoles | 0.9 |
| enzyme_bla_z | rifamycins | 0.9 |
| enzyme_bla_z | macrocycles | 0.9 |
| enzyme_bla_z | steroid_antibacterials | 0.9 |
| enzyme_bla_z | pleuromutilins | 0.9 |
| enzyme_bla_z | other | 0.9 |
| enzyme_oxa_acinetobacter | pen | 0.8 |
| enzyme_oxa_acinetobacter | bli | 0.8 |
| enzyme_oxa_acinetobacter | bli_anti_pseudomonal | 0.8 |
| enzyme_oxa_acinetobacter | bli_sulbactam | 0.8 |
| enzyme_oxa_acinetobacter | c1_2g | 0.8 |
| enzyme_oxa_acinetobacter | c3g | 0.8 |
| enzyme_oxa_acinetobacter | c3g_bli | 0.8 |
| enzyme_oxa_acinetobacter | c4g | 0.8 |
| enzyme_oxa_acinetobacter | anti_mrsa_ceph | 0.8 |
| enzyme_oxa_acinetobacter | siderophore_ceph | 0.8 |
| enzyme_oxa_acinetobacter | cft_avi | 0.8 |
| enzyme_oxa_acinetobacter | mer_vab | 0.8 |
| enzyme_oxa_acinetobacter | azt_avi | 0.8 |
| enzyme_oxa_acinetobacter | carb_group1 | 0.8 |
| enzyme_oxa_acinetobacter | carb_group2 | 0.8 |
| enzyme_oxa_acinetobacter | mono | 0.8 |
| enzyme_oxa_acinetobacter | fq | 0.8 |
| enzyme_oxa_acinetobacter | ag_group1 | 0.8 |
| enzyme_oxa_acinetobacter | ag_group2 | 0.8 |
| enzyme_oxa_acinetobacter | mls | 0.8 |
| enzyme_oxa_acinetobacter | lincosamides | 0.8 |
| enzyme_oxa_acinetobacter | glyc | 0.8 |
| enzyme_oxa_acinetobacter | lipoglycopeptides | 0.8 |
| enzyme_oxa_acinetobacter | tet | 0.8 |
| enzyme_oxa_acinetobacter | glycylcyclines | 0.8 |
| enzyme_oxa_acinetobacter | poly | 0.8 |
| enzyme_oxa_acinetobacter | oxa | 0.8 |
| enzyme_oxa_acinetobacter | chl | 0.8 |
| enzyme_oxa_acinetobacter | sulf | 0.8 |
| enzyme_oxa_acinetobacter | lipopeptides | 0.8 |
| enzyme_oxa_acinetobacter | streptogramins | 0.8 |
| enzyme_oxa_acinetobacter | nitrofurans | 0.8 |
| enzyme_oxa_acinetobacter | phosphonic_acids | 0.8 |
| enzyme_oxa_acinetobacter | nitroimidazoles | 0.8 |
| enzyme_oxa_acinetobacter | rifamycins | 0.8 |
| enzyme_oxa_acinetobacter | macrocycles | 0.8 |
| enzyme_oxa_acinetobacter | steroid_antibacterials | 0.8 |
| enzyme_oxa_acinetobacter | pleuromutilins | 0.8 |
| enzyme_oxa_acinetobacter | other | 0.8 |
| mutation_23s_rrna | pen | 0.8 |
| mutation_23s_rrna | bli | 0.8 |
| mutation_23s_rrna | bli_anti_pseudomonal | 0.8 |
| mutation_23s_rrna | bli_sulbactam | 0.8 |
| mutation_23s_rrna | c1_2g | 0.8 |
| mutation_23s_rrna | c3g | 0.8 |
| mutation_23s_rrna | c3g_bli | 0.8 |
| mutation_23s_rrna | c4g | 0.8 |
| mutation_23s_rrna | anti_mrsa_ceph | 0.8 |
| mutation_23s_rrna | siderophore_ceph | 0.8 |
| mutation_23s_rrna | cft_avi | 0.8 |
| mutation_23s_rrna | mer_vab | 0.8 |
| mutation_23s_rrna | azt_avi | 0.8 |
| mutation_23s_rrna | carb_group1 | 0.8 |
| mutation_23s_rrna | carb_group2 | 0.8 |
| mutation_23s_rrna | mono | 0.8 |
| mutation_23s_rrna | fq | 0.8 |
| mutation_23s_rrna | ag_group1 | 0.8 |
| mutation_23s_rrna | ag_group2 | 0.8 |
| mutation_23s_rrna | mls | 0.8 |
| mutation_23s_rrna | lincosamides | 0.8 |
| mutation_23s_rrna | glyc | 0.8 |
| mutation_23s_rrna | lipoglycopeptides | 0.8 |
| mutation_23s_rrna | tet | 0.8 |
| mutation_23s_rrna | glycylcyclines | 0.8 |
| mutation_23s_rrna | poly | 0.8 |
| mutation_23s_rrna | oxa | 0.8 |
| mutation_23s_rrna | chl | 0.8 |
| mutation_23s_rrna | sulf | 0.8 |
| mutation_23s_rrna | lipopeptides | 0.8 |
| mutation_23s_rrna | streptogramins | 0.8 |
| mutation_23s_rrna | nitrofurans | 0.8 |
| mutation_23s_rrna | phosphonic_acids | 0.8 |
| mutation_23s_rrna | nitroimidazoles | 0.8 |
| mutation_23s_rrna | rifamycins | 0.8 |
| mutation_23s_rrna | macrocycles | 0.8 |
| mutation_23s_rrna | steroid_antibacterials | 0.8 |
| mutation_23s_rrna | pleuromutilins | 0.8 |
| mutation_23s_rrna | other | 0.8 |
| efflux_tet_abc | pen | 0.7 |
| efflux_tet_abc | bli | 0.7 |
| efflux_tet_abc | bli_anti_pseudomonal | 0.7 |
| efflux_tet_abc | bli_sulbactam | 0.7 |
| efflux_tet_abc | c1_2g | 0.7 |
| efflux_tet_abc | c3g | 0.7 |
| efflux_tet_abc | c3g_bli | 0.7 |
| efflux_tet_abc | c4g | 0.7 |
| efflux_tet_abc | anti_mrsa_ceph | 0.7 |
| efflux_tet_abc | siderophore_ceph | 0.7 |
| efflux_tet_abc | cft_avi | 0.7 |
| efflux_tet_abc | mer_vab | 0.7 |
| efflux_tet_abc | azt_avi | 0.7 |
| efflux_tet_abc | carb_group1 | 0.7 |
| efflux_tet_abc | carb_group2 | 0.7 |
| efflux_tet_abc | mono | 0.7 |
| efflux_tet_abc | fq | 0.7 |
| efflux_tet_abc | ag_group1 | 0.7 |
| efflux_tet_abc | ag_group2 | 0.7 |
| efflux_tet_abc | mls | 0.7 |
| efflux_tet_abc | lincosamides | 0.7 |
| efflux_tet_abc | glyc | 0.7 |
| efflux_tet_abc | lipoglycopeptides | 0.7 |
| efflux_tet_abc | tet | 0.7 |
| efflux_tet_abc | glycylcyclines | 0.7 |
| efflux_tet_abc | poly | 0.7 |
| efflux_tet_abc | oxa | 0.7 |
| efflux_tet_abc | chl | 0.7 |
| efflux_tet_abc | sulf | 0.7 |
| efflux_tet_abc | lipopeptides | 0.7 |
| efflux_tet_abc | streptogramins | 0.7 |
| efflux_tet_abc | nitrofurans | 0.7 |
| efflux_tet_abc | phosphonic_acids | 0.7 |
| efflux_tet_abc | nitroimidazoles | 0.7 |
| efflux_tet_abc | rifamycins | 0.7 |
| efflux_tet_abc | macrocycles | 0.7 |
| efflux_tet_abc | steroid_antibacterials | 0.7 |
| efflux_tet_abc | pleuromutilins | 0.7 |
| efflux_tet_abc | other | 0.7 |
| mutation_pbp_mosaic | pen | 0.8 |
| mutation_pbp_mosaic | bli | 0.7 |
| mutation_pbp_mosaic | bli_anti_pseudomonal | 0.7 |
| mutation_pbp_mosaic | bli_sulbactam | 0.7 |
| mutation_pbp_mosaic | c1_2g | 0.6 |
| mutation_pbp_mosaic | c3g | 0.3 |
| mutation_pbp_mosaic | c3g_bli | 0.3 |
| mutation_pbp_mosaic | c4g | 0.15 |
| mutation_pbp_mosaic | anti_mrsa_ceph | 0.15 |
| mutation_pbp_mosaic | siderophore_ceph | 0.15 |
| mutation_pbp_mosaic | cft_avi | 0.1 |
| mutation_pbp_mosaic | azt_avi | 0.5 |
| mutation_pbp_mosaic | mono | 0.5 |
| mutation_pbp_mosaic | lipopeptides | 0.5 |
| mutation_pbp_mosaic | streptogramins | 0.5 |
| mutation_pbp_mosaic | nitrofurans | 0.5 |
| mutation_pbp_mosaic | phosphonic_acids | 0.5 |
| mutation_pbp_mosaic | nitroimidazoles | 0.5 |
| mutation_pbp_mosaic | rifamycins | 0.5 |
| mutation_pbp_mosaic | macrocycles | 0.5 |
| mutation_pbp_mosaic | steroid_antibacterials | 0.5 |
| mutation_pbp_mosaic | pleuromutilins | 0.5 |
| efflux_mtr_cde | pen | 0.3 |
| efflux_mtr_cde | mls | 0.5 |
| efflux_mtr_cde | tet | 0.4 |
| efflux_mtr_cde | chl | 0.4 |
| efflux_mtr_cde | lipopeptides | 0.4 |
| efflux_mtr_cde | streptogramins | 0.4 |
| efflux_mtr_cde | nitrofurans | 0.4 |
| efflux_mtr_cde | phosphonic_acids | 0.4 |
| efflux_mtr_cde | nitroimidazoles | 0.4 |
| efflux_mtr_cde | rifamycins | 0.4 |
| efflux_mtr_cde | macrocycles | 0.4 |
| efflux_mtr_cde | steroid_antibacterials | 0.4 |
| efflux_mtr_cde | pleuromutilins | 0.4 |
| as_yet_unknown | pen | 0.5 |
| as_yet_unknown | bli | 0.5 |
| as_yet_unknown | bli_anti_pseudomonal | 0.5 |
| as_yet_unknown | bli_sulbactam | 0.5 |
| as_yet_unknown | c1_2g | 0.5 |
| as_yet_unknown | c3g | 0.5 |
| as_yet_unknown | c3g_bli | 0.5 |
| as_yet_unknown | c4g | 0.5 |
| as_yet_unknown | anti_mrsa_ceph | 0.5 |
| as_yet_unknown | siderophore_ceph | 0.5 |
| as_yet_unknown | cft_avi | 0.5 |
| as_yet_unknown | mer_vab | 0.5 |
| as_yet_unknown | azt_avi | 0.5 |
| as_yet_unknown | carb_group1 | 0.5 |
| as_yet_unknown | carb_group2 | 0.5 |
| as_yet_unknown | mono | 0.5 |
| as_yet_unknown | fq | 0.5 |
| as_yet_unknown | ag_group1 | 0.5 |
| as_yet_unknown | ag_group2 | 0.5 |
| as_yet_unknown | mls | 0.5 |
| as_yet_unknown | lincosamides | 0.5 |
| as_yet_unknown | glyc | 0.5 |
| as_yet_unknown | lipoglycopeptides | 0.5 |
| as_yet_unknown | tet | 0.5 |
| as_yet_unknown | glycylcyclines | 0.5 |
| as_yet_unknown | poly | 0.5 |
| as_yet_unknown | oxa | 0.5 |
| as_yet_unknown | chl | 0.5 |
| as_yet_unknown | sulf | 0.5 |
| as_yet_unknown | lipopeptides | 0.5 |
| as_yet_unknown | streptogramins | 0.5 |
| as_yet_unknown | nitrofurans | 0.5 |
| as_yet_unknown | phosphonic_acids | 0.5 |
| as_yet_unknown | nitroimidazoles | 0.5 |
| as_yet_unknown | rifamycins | 0.5 |
| as_yet_unknown | macrocycles | 0.5 |
| as_yet_unknown | steroid_antibacterials | 0.5 |
| as_yet_unknown | pleuromutilins | 0.5 |
| as_yet_unknown | other | 0.5 |

#### BacteriaΓÇôMechanism Emergence Rates

De novo emergence rate per day for each bacteriaΓÇômechanism pair. Only non-zero entries shown.

| Bacteria | Mechanism | Emergence rate/day |
| --- | ---: | ---: |
| acinetobacter_baumannii | enzyme_esbl_ctx_m | 0.003 |
| acinetobacter_baumannii | enzyme_esbl_tem | 0.003 |
| acinetobacter_baumannii | enzyme_esbl_shv | 0.003 |
| acinetobacter_baumannii | enzyme_kpc | 0.002 |
| acinetobacter_baumannii | enzyme_ndm_vim | 0.002 |
| acinetobacter_baumannii | enzyme_oxa_48 | 0.002 |
| acinetobacter_baumannii | enzyme_ampc_cmy | 0.003 |
| acinetobacter_baumannii | enzyme_ampc_dha | 0.003 |
| acinetobacter_baumannii | mutation_gyra_primary | 0.09 |
| acinetobacter_baumannii | mutation_gyra_parc_secondary | 0.09 |
| acinetobacter_baumannii | protection_qnr | 0.09 |
| acinetobacter_baumannii | enzyme_16s_rrmt | 2 |
| acinetobacter_baumannii | enzyme_cat | 0.03 |
| acinetobacter_baumannii | modification_mcr_1 | 0.1 |
| acinetobacter_baumannii | global_efflux_pump | 0.025 |
| acinetobacter_baumannii | global_porin_loss | 1e-4 |
| acinetobacter_baumannii | mutation_folate_pathway | 0.1 |
| acinetobacter_baumannii | enzyme_fos_a | 0.0035 |
| acinetobacter_baumannii | mutation_rpo_b | 1 |
| acinetobacter_baumannii | protection_tet_m | 0.03 |
| acinetobacter_baumannii | enzyme_aac_aph | 2 |
| acinetobacter_baumannii | enzyme_oxa_acinetobacter | 0.002 |
| acinetobacter_baumannii | efflux_tet_abc | 0.03 |
| acinetobacter_baumannii | mutation_pbp_mosaic | 0.003 |
| citrobacter_spp. | enzyme_esbl_ctx_m | 1e-4 |
| citrobacter_spp. | enzyme_esbl_tem | 1e-4 |
| citrobacter_spp. | enzyme_esbl_shv | 1e-4 |
| citrobacter_spp. | enzyme_kpc | 1e-4 |
| citrobacter_spp. | enzyme_ndm_vim | 1e-4 |
| citrobacter_spp. | enzyme_oxa_48 | 1e-4 |
| citrobacter_spp. | enzyme_ampc_cmy | 1e-4 |
| citrobacter_spp. | enzyme_ampc_dha | 1e-4 |
| citrobacter_spp. | mutation_gyra_primary | 0.05 |
| citrobacter_spp. | mutation_gyra_parc_secondary | 0.05 |
| citrobacter_spp. | protection_qnr | 0.05 |
| citrobacter_spp. | enzyme_16s_rrmt | 0.1 |
| citrobacter_spp. | enzyme_cat | 0.015 |
| citrobacter_spp. | efflux_acrab_tolc | 0.01 |
| citrobacter_spp. | efflux_mexxy_oprm | 1e-4 |
| citrobacter_spp. | modification_mcr_1 | 0.05 |
| citrobacter_spp. | global_efflux_pump | 0.005 |
| citrobacter_spp. | global_porin_loss | 3e-4 |
| citrobacter_spp. | mutation_folate_pathway | 0.003 |
| citrobacter_spp. | mutation_nitroreductase | 0.01 |
| citrobacter_spp. | enzyme_fos_a | 0.003 |
| citrobacter_spp. | mutation_rpo_b | 0.015 |
| citrobacter_spp. | protection_tet_m | 0.005 |
| citrobacter_spp. | enzyme_aac_aph | 0.1 |
| citrobacter_spp. | efflux_tet_abc | 0.005 |
| citrobacter_spp. | mutation_pbp_mosaic | 1e-4 |
| enterobacter_spp. | enzyme_esbl_ctx_m | 0.001 |
| enterobacter_spp. | enzyme_esbl_tem | 0.001 |
| enterobacter_spp. | enzyme_esbl_shv | 0.001 |
| enterobacter_spp. | enzyme_kpc | 3e-4 |
| enterobacter_spp. | enzyme_ndm_vim | 3e-4 |
| enterobacter_spp. | enzyme_oxa_48 | 3e-4 |
| enterobacter_spp. | enzyme_ampc_cmy | 0.001 |
| enterobacter_spp. | enzyme_ampc_dha | 0.001 |
| enterobacter_spp. | mutation_gyra_primary | 0.01 |
| enterobacter_spp. | mutation_gyra_parc_secondary | 0.001 |
| enterobacter_spp. | protection_qnr | 0.01 |
| enterobacter_spp. | enzyme_16s_rrmt | 0.2 |
| enterobacter_spp. | enzyme_cat | 0.003 |
| enterobacter_spp. | efflux_acrab_tolc | 0.01 |
| enterobacter_spp. | modification_mcr_1 | 0.03 |
| enterobacter_spp. | global_efflux_pump | 0.01 |
| enterobacter_spp. | global_porin_loss | 5e-4 |
| enterobacter_spp. | mutation_folate_pathway | 0.001 |
| enterobacter_spp. | mutation_nitroreductase | 0.005 |
| enterobacter_spp. | enzyme_fos_a | 0.02 |
| enterobacter_spp. | mutation_rpo_b | 0.03 |
| enterobacter_spp. | protection_tet_m | 0.008 |
| enterobacter_spp. | enzyme_aac_aph | 0.2 |
| enterobacter_spp. | efflux_tet_abc | 0.008 |
| enterobacter_spp. | mutation_pbp_mosaic | 3e-4 |
| enterococcus_faecalis | target_site_pbp2a_meca | 5e-6 |
| enterococcus_faecalis | target_site_van_a | 0.02 |
| enterococcus_faecalis | target_site_van_b | 0.02 |
| enterococcus_faecalis | mutation_gyra_primary | 0.01 |
| enterococcus_faecalis | mutation_gyra_parc_secondary | 0.01 |
| enterococcus_faecalis | target_site_erm_b | 0.002 |
| enterococcus_faecalis | target_site_cfr | 0.003 |
| enterococcus_faecalis | enzyme_cat | 2e-4 |
| enterococcus_faecalis | global_efflux_pump | 0.001 |
| enterococcus_faecalis | global_porin_loss | 1e-4 |
| enterococcus_faecalis | mutation_folate_pathway | 0.02 |
| enterococcus_faecalis | mutation_nitroreductase | 0.02 |
| enterococcus_faecalis | mutation_mpr_f | 0.005 |
| enterococcus_faecalis | mutation_rpo_b | 0.005 |
| enterococcus_faecalis | protection_fus_b | 5e-4 |
| enterococcus_faecalis | protection_tet_m | 0.004 |
| enterococcus_faecalis | enzyme_aac_aph | 2e-5 |
| enterococcus_faecalis | mutation_23s_rrna | 3e-4 |
| enterococcus_faecalis | mutation_pbp_mosaic | 2e-6 |
| enterococcus_faecium | target_site_van_a | 0.001 |
| enterococcus_faecium | target_site_van_b | 0.001 |
| enterococcus_faecium | mutation_gyra_primary | 0.015 |
| enterococcus_faecium | mutation_gyra_parc_secondary | 0.015 |
| enterococcus_faecium | enzyme_16s_rrmt | 0.01 |
| enterococcus_faecium | target_site_erm_b | 0.01 |
| enterococcus_faecium | target_site_cfr | 0.005 |
| enterococcus_faecium | enzyme_cat | 0.01 |
| enterococcus_faecium | global_efflux_pump | 0.015 |
| enterococcus_faecium | mutation_folate_pathway | 0.01 |
| enterococcus_faecium | mutation_nitroreductase | 0.3 |
| enterococcus_faecium | enzyme_fos_a | 0.5 |
| enterococcus_faecium | mutation_mpr_f | 0.03 |
| enterococcus_faecium | mutation_rpo_b | 0.3 |
| enterococcus_faecium | protection_fus_b | 0.012 |
| enterococcus_faecium | protection_tet_m | 0.01 |
| enterococcus_faecium | enzyme_aac_aph | 0.1 |
| enterococcus_faecium | mutation_23s_rrna | 0.002 |
| enterococcus_faecium | mutation_pbp_mosaic | 0.03 |
| enterococcus_faecium | efflux_mtr_cde | 0.0012 |
| escherichia_coli | enzyme_esbl_ctx_m | 5e-7 |
| escherichia_coli | enzyme_esbl_tem | 5e-7 |
| escherichia_coli | enzyme_esbl_shv | 5e-7 |
| escherichia_coli | enzyme_kpc | 5e-8 |
| escherichia_coli | enzyme_ndm_vim | 5e-8 |
| escherichia_coli | enzyme_oxa_48 | 5e-8 |
| escherichia_coli | enzyme_ampc_cmy | 5e-7 |
| escherichia_coli | enzyme_ampc_dha | 5e-7 |
| escherichia_coli | mutation_gyra_primary | 0.02 |
| escherichia_coli | mutation_gyra_parc_secondary | 0.003 |
| escherichia_coli | protection_qnr | 0.03 |
| escherichia_coli | enzyme_16s_rrmt | 0.01 |
| escherichia_coli | enzyme_cat | 1e-6 |
| escherichia_coli | efflux_acrab_tolc | 0.003 |
| escherichia_coli | modification_mcr_1 | 1e-5 |
| escherichia_coli | global_efflux_pump | 0.003 |
| escherichia_coli | global_porin_loss | 1e-7 |
| escherichia_coli | mutation_folate_pathway | 2e-4 |
| escherichia_coli | mutation_nitroreductase | 1e-6 |
| escherichia_coli | enzyme_fos_a | 3e-6 |
| escherichia_coli | mutation_rpo_b | 1e-4 |
| escherichia_coli | protection_tet_m | 0.01 |
| escherichia_coli | enzyme_aac_aph | 0.01 |
| escherichia_coli | efflux_tet_abc | 0.003 |
| klebsiella_pneumoniae | enzyme_esbl_ctx_m | 1e-5 |
| klebsiella_pneumoniae | enzyme_esbl_tem | 1e-5 |
| klebsiella_pneumoniae | enzyme_esbl_shv | 1e-5 |
| klebsiella_pneumoniae | enzyme_kpc | 3e-7 |
| klebsiella_pneumoniae | enzyme_ndm_vim | 3e-7 |
| klebsiella_pneumoniae | enzyme_oxa_48 | 3e-7 |
| klebsiella_pneumoniae | enzyme_ampc_cmy | 1e-5 |
| klebsiella_pneumoniae | enzyme_ampc_dha | 1e-5 |
| klebsiella_pneumoniae | mutation_gyra_primary | 0.001 |
| klebsiella_pneumoniae | mutation_gyra_parc_secondary | 0.0015 |
| klebsiella_pneumoniae | protection_qnr | 0.0015 |
| klebsiella_pneumoniae | enzyme_16s_rrmt | 5e-4 |
| klebsiella_pneumoniae | enzyme_cat | 0.002 |
| klebsiella_pneumoniae | efflux_acrab_tolc | 0.0015 |
| klebsiella_pneumoniae | porin_loss_ompk35_36 | 3e-7 |
| klebsiella_pneumoniae | modification_mcr_1 | 0.01 |
| klebsiella_pneumoniae | global_efflux_pump | 0.0015 |
| klebsiella_pneumoniae | global_porin_loss | 3e-9 |
| klebsiella_pneumoniae | mutation_folate_pathway | 0.001 |
| klebsiella_pneumoniae | mutation_nitroreductase | 0.002 |
| klebsiella_pneumoniae | enzyme_fos_a | 0.003 |
| klebsiella_pneumoniae | mutation_rpo_b | 2e-4 |
| klebsiella_pneumoniae | protection_tet_m | 0.001 |
| klebsiella_pneumoniae | enzyme_aac_aph | 5e-4 |
| klebsiella_pneumoniae | efflux_tet_abc | 0.001 |
| morganella_spp. | enzyme_esbl_ctx_m | 3e-5 |
| morganella_spp. | enzyme_esbl_tem | 3e-5 |
| morganella_spp. | enzyme_esbl_shv | 3e-5 |
| morganella_spp. | enzyme_kpc | 1e-5 |
| morganella_spp. | enzyme_ndm_vim | 1e-5 |
| morganella_spp. | enzyme_oxa_48 | 1e-5 |
| morganella_spp. | enzyme_ampc_cmy | 3e-5 |
| morganella_spp. | enzyme_ampc_dha | 3e-5 |
| morganella_spp. | mutation_gyra_primary | 0.15 |
| morganella_spp. | mutation_gyra_parc_secondary | 0.15 |
| morganella_spp. | protection_qnr | 0.15 |
| morganella_spp. | enzyme_16s_rrmt | 0.8 |
| morganella_spp. | enzyme_cat | 0.1 |
| morganella_spp. | efflux_acrab_tolc | 0.003 |
| morganella_spp. | efflux_mexxy_oprm | 1e-5 |
| morganella_spp. | modification_mcr_1 | 0.02 |
| morganella_spp. | global_efflux_pump | 0.01 |
| morganella_spp. | global_porin_loss | 1e-5 |
| morganella_spp. | mutation_folate_pathway | 0.003 |
| morganella_spp. | mutation_nitroreductase | 0.03 |
| morganella_spp. | enzyme_fos_a | 0.003 |
| morganella_spp. | mutation_rpo_b | 0.01 |
| morganella_spp. | protection_tet_m | 0.003 |
| morganella_spp. | enzyme_aac_aph | 0.8 |
| morganella_spp. | efflux_tet_abc | 0.003 |
| morganella_spp. | mutation_pbp_mosaic | 1e-5 |
| proteus_spp. | enzyme_esbl_ctx_m | 3e-5 |
| proteus_spp. | enzyme_esbl_tem | 3e-5 |
| proteus_spp. | enzyme_esbl_shv | 3e-5 |
| proteus_spp. | enzyme_kpc | 1.5e-5 |
| proteus_spp. | enzyme_ndm_vim | 1.5e-5 |
| proteus_spp. | enzyme_oxa_48 | 1.5e-5 |
| proteus_spp. | enzyme_ampc_cmy | 3e-5 |
| proteus_spp. | enzyme_ampc_dha | 3e-5 |
| proteus_spp. | mutation_gyra_primary | 0.04 |
| proteus_spp. | mutation_gyra_parc_secondary | 0.01 |
| proteus_spp. | protection_qnr | 0.01 |
| proteus_spp. | enzyme_16s_rrmt | 0.8 |
| proteus_spp. | enzyme_cat | 3.5e-4 |
| proteus_spp. | efflux_acrab_tolc | 0.003 |
| proteus_spp. | modification_mcr_1 | 5e-4 |
| proteus_spp. | global_efflux_pump | 0.03 |
| proteus_spp. | global_porin_loss | 1.5e-4 |
| proteus_spp. | mutation_folate_pathway | 0.002 |
| proteus_spp. | mutation_nitroreductase | 1.5e-5 |
| proteus_spp. | enzyme_fos_a | 0.0015 |
| proteus_spp. | mutation_rpo_b | 1e-4 |
| proteus_spp. | protection_tet_m | 0.002 |
| proteus_spp. | enzyme_aac_aph | 0.08 |
| proteus_spp. | efflux_tet_abc | 1e-6 |
| serratia_spp. | enzyme_esbl_ctx_m | 5e-4 |
| serratia_spp. | enzyme_esbl_tem | 5e-4 |
| serratia_spp. | enzyme_esbl_shv | 5e-4 |
| serratia_spp. | enzyme_kpc | 1e-4 |
| serratia_spp. | enzyme_ndm_vim | 1e-4 |
| serratia_spp. | enzyme_oxa_48 | 1e-4 |
| serratia_spp. | enzyme_ampc_cmy | 5e-4 |
| serratia_spp. | enzyme_ampc_dha | 5e-4 |
| serratia_spp. | mutation_gyra_primary | 0.025 |
| serratia_spp. | mutation_gyra_parc_secondary | 0.01 |
| serratia_spp. | protection_qnr | 0.01 |
| serratia_spp. | enzyme_16s_rrmt | 0.5 |
| serratia_spp. | enzyme_cat | 3e-5 |
| serratia_spp. | efflux_acrab_tolc | 0.005 |
| serratia_spp. | modification_mcr_1 | 0.03 |
| serratia_spp. | global_efflux_pump | 0.005 |
| serratia_spp. | global_porin_loss | 3e-5 |
| serratia_spp. | mutation_folate_pathway | 0.001 |
| serratia_spp. | mutation_nitroreductase | 0.002 |
| serratia_spp. | enzyme_fos_a | 2e-4 |
| serratia_spp. | mutation_rpo_b | 0.03 |
| serratia_spp. | protection_tet_m | 0.003 |
| serratia_spp. | enzyme_aac_aph | 0.5 |
| serratia_spp. | efflux_tet_abc | 0.003 |
| p_stuartii | enzyme_esbl_ctx_m | 1.9e-4 |
| p_stuartii | enzyme_esbl_tem | 1.9e-4 |
| p_stuartii | enzyme_esbl_shv | 3.7e-5 |
| p_stuartii | enzyme_kpc | 1.9e-5 |
| p_stuartii | enzyme_ndm_vim | 1.9e-5 |
| p_stuartii | enzyme_oxa_48 | 1.9e-5 |
| p_stuartii | enzyme_ampc_cmy | 3.8e-4 |
| p_stuartii | enzyme_ampc_dha | 1.9e-4 |
| p_stuartii | mutation_gyra_primary | 3.8e-4 |
| p_stuartii | mutation_gyra_parc_secondary | 1.9e-4 |
| p_stuartii | protection_qnr | 1.9e-4 |
| p_stuartii | enzyme_16s_rrmt | 3.8e-6 |
| p_stuartii | enzyme_cat | 1.9e-4 |
| p_stuartii | efflux_acrab_tolc | 1.9e-4 |
| p_stuartii | modification_mcr_1 | 1.9e-5 |
| p_stuartii | global_efflux_pump | 3.7e-5 |
| p_stuartii | global_porin_loss | 3.7e-5 |
| p_stuartii | mutation_folate_pathway | 3.8e-4 |
| p_stuartii | mutation_nitroreductase | 3.7e-5 |
| p_stuartii | enzyme_fos_a | 3.7e-5 |
| p_stuartii | mutation_rpo_b | 3.8e-6 |
| p_stuartii | protection_tet_m | 3.8e-4 |
| p_stuartii | enzyme_aac_aph | 5e-11 |
| p_stuartii | efflux_tet_abc | 5e-11 |
| pseudomonas_aeruginosa | enzyme_esbl_ctx_m | 1.5e-4 |
| pseudomonas_aeruginosa | enzyme_esbl_tem | 1.5e-4 |
| pseudomonas_aeruginosa | enzyme_esbl_shv | 1.5e-4 |
| pseudomonas_aeruginosa | enzyme_kpc | 1.5e-4 |
| pseudomonas_aeruginosa | enzyme_ndm_vim | 1.5e-4 |
| pseudomonas_aeruginosa | enzyme_oxa_48 | 1.5e-4 |
| pseudomonas_aeruginosa | enzyme_ampc_cmy | 1.5e-4 |
| pseudomonas_aeruginosa | enzyme_ampc_dha | 1.5e-4 |
| pseudomonas_aeruginosa | mutation_gyra_primary | 0.002 |
| pseudomonas_aeruginosa | mutation_gyra_parc_secondary | 0.002 |
| pseudomonas_aeruginosa | protection_qnr | 0.002 |
| pseudomonas_aeruginosa | enzyme_16s_rrmt | 1e-4 |
| pseudomonas_aeruginosa | target_site_erm_b | 2e-5 |
| pseudomonas_aeruginosa | target_site_cfr | 2e-5 |
| pseudomonas_aeruginosa | enzyme_cat | 2e-4 |
| pseudomonas_aeruginosa | efflux_mexxy_oprm | 1.5e-4 |
| pseudomonas_aeruginosa | porin_loss_oprd | 1.5e-4 |
| pseudomonas_aeruginosa | modification_mcr_1 | 5e-4 |
| pseudomonas_aeruginosa | global_efflux_pump | 5e-4 |
| pseudomonas_aeruginosa | global_porin_loss | 2e-4 |
| pseudomonas_aeruginosa | mutation_folate_pathway | 0.002 |
| pseudomonas_aeruginosa | mutation_nitroreductase | 2e-5 |
| pseudomonas_aeruginosa | enzyme_fos_a | 0.01 |
| pseudomonas_aeruginosa | mutation_rpo_b | 5e-4 |
| pseudomonas_aeruginosa | protection_tet_m | 0.002 |
| pseudomonas_aeruginosa | enzyme_aac_aph | 1e-4 |
| pseudomonas_aeruginosa | efflux_tet_abc | 1e-5 |
| stenotrophomonas_maltophilia | enzyme_esbl_ctx_m | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_esbl_tem | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_esbl_shv | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_kpc | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_ndm_vim | 0.1 |
| stenotrophomonas_maltophilia | enzyme_oxa_48 | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_ampc_cmy | 0.05 |
| stenotrophomonas_maltophilia | enzyme_ampc_dha | 2e-4 |
| stenotrophomonas_maltophilia | mutation_gyra_primary | 5e-8 |
| stenotrophomonas_maltophilia | mutation_gyra_parc_secondary | 5e-8 |
| stenotrophomonas_maltophilia | protection_qnr | 1.5e-5 |
| stenotrophomonas_maltophilia | enzyme_16s_rrmt | 0.01 |
| stenotrophomonas_maltophilia | target_site_erm_b | 0.1 |
| stenotrophomonas_maltophilia | enzyme_cat | 5e-7 |
| stenotrophomonas_maltophilia | modification_mcr_1 | 0.05 |
| stenotrophomonas_maltophilia | global_efflux_pump | 1e-7 |
| stenotrophomonas_maltophilia | global_porin_loss | 1e-13 |
| stenotrophomonas_maltophilia | mutation_folate_pathway | 0.003 |
| stenotrophomonas_maltophilia | mutation_nitroreductase | 10 |
| stenotrophomonas_maltophilia | enzyme_fos_a | 0.03 |
| stenotrophomonas_maltophilia | mutation_rpo_b | 2e-4 |
| stenotrophomonas_maltophilia | protection_tet_m | 2e-7 |
| stenotrophomonas_maltophilia | enzyme_aac_aph | 0.05 |
| stenotrophomonas_maltophilia | efflux_tet_abc | 2e-9 |
| staphylococcus_aureus | target_site_pbp2a_meca | 2e-7 |
| staphylococcus_aureus | target_site_van_a | 1e-6 |
| staphylococcus_aureus | target_site_van_b | 2e-7 |
| staphylococcus_aureus | mutation_gyra_primary | 0.02 |
| staphylococcus_aureus | mutation_gyra_parc_secondary | 0.01 |
| staphylococcus_aureus | enzyme_16s_rrmt | 0.2 |
| staphylococcus_aureus | target_site_erm_b | 1.5e-4 |
| staphylococcus_aureus | target_site_cfr | 1e-4 |
| staphylococcus_aureus | enzyme_cat | 4e-5 |
| staphylococcus_aureus | global_efflux_pump | 0.003 |
| staphylococcus_aureus | mutation_folate_pathway | 0.003 |
| staphylococcus_aureus | mutation_nitroreductase | 0.003 |
| staphylococcus_aureus | enzyme_fos_a | 0.003 |
| staphylococcus_aureus | mutation_mpr_f | 0.003 |
| staphylococcus_aureus | mutation_rpo_b | 0.003 |
| staphylococcus_aureus | protection_fus_b | 0.003 |
| staphylococcus_aureus | protection_tet_m | 0.002 |
| staphylococcus_aureus | enzyme_aac_aph | 0.2 |
| staphylococcus_aureus | enzyme_bla_z | 5e-7 |
| staphylococcus_epidermidis | target_site_pbp2a_meca | 1e-6 |
| staphylococcus_epidermidis | target_site_van_a | 1e-8 |
| staphylococcus_epidermidis | target_site_van_b | 1e-8 |
| staphylococcus_epidermidis | mutation_gyra_primary | 1e-5 |
| staphylococcus_epidermidis | mutation_gyra_parc_secondary | 1e-6 |
| staphylococcus_epidermidis | target_site_erm_b | 1e-5 |
| staphylococcus_epidermidis | target_site_cfr | 1e-7 |
| staphylococcus_epidermidis | enzyme_cat | 1e-6 |
| staphylococcus_epidermidis | global_efflux_pump | 1e-6 |
| staphylococcus_epidermidis | mutation_folate_pathway | 1e-5 |
| staphylococcus_epidermidis | mutation_mpr_f | 1e-6 |
| staphylococcus_epidermidis | mutation_rpo_b | 0.001 |
| staphylococcus_epidermidis | protection_fus_b | 3e-6 |
| staphylococcus_epidermidis | protection_tet_m | 3e-5 |
| staphylococcus_epidermidis | enzyme_aac_aph | 1e-9 |
| staphylococcus_epidermidis | enzyme_bla_z | 1e-7 |
| streptococcus_pneumoniae | target_site_pbp2a_meca | 3e-8 |
| streptococcus_pneumoniae | mutation_gyra_primary | 3e-5 |
| streptococcus_pneumoniae | mutation_gyra_parc_secondary | 3e-5 |
| streptococcus_pneumoniae | target_site_erm_b | 1e-6 |
| streptococcus_pneumoniae | enzyme_cat | 3e-4 |
| streptococcus_pneumoniae | global_efflux_pump | 1e-5 |
| streptococcus_pneumoniae | global_porin_loss | 3e-8 |
| streptococcus_pneumoniae | mutation_folate_pathway | 0.003 |
| streptococcus_pneumoniae | mutation_rpo_b | 0.001 |
| streptococcus_pneumoniae | protection_tet_m | 0.001 |
| streptococcus_pneumoniae | enzyme_bla_z | 3e-8 |
| streptococcus_pneumoniae | mutation_23s_rrna | 5e-7 |
| streptococcus_pneumoniae | mutation_pbp_mosaic | 3e-8 |
| salmonella_enterica_serovar_typhi | enzyme_esbl_ctx_m | 0.005 |
| salmonella_enterica_serovar_typhi | enzyme_esbl_tem | 0.15 |
| salmonella_enterica_serovar_typhi | enzyme_esbl_shv | 0.0015 |
| salmonella_enterica_serovar_typhi | enzyme_kpc | 3e-5 |
| salmonella_enterica_serovar_typhi | enzyme_ndm_vim | 3e-5 |
| salmonella_enterica_serovar_typhi | enzyme_oxa_48 | 3e-5 |
| salmonella_enterica_serovar_typhi | enzyme_ampc_cmy | 5e-4 |
| salmonella_enterica_serovar_typhi | enzyme_ampc_dha | 5e-4 |
| salmonella_enterica_serovar_typhi | mutation_gyra_primary | 0.8 |
| salmonella_enterica_serovar_typhi | mutation_gyra_parc_secondary | 0.05 |
| salmonella_enterica_serovar_typhi | protection_qnr | 0.03 |
| salmonella_enterica_serovar_typhi | enzyme_16s_rrmt | 1e-4 |
| salmonella_enterica_serovar_typhi | enzyme_cat | 0.015 |
| salmonella_enterica_serovar_typhi | efflux_acrab_tolc | 0.05 |
| salmonella_enterica_serovar_typhi | modification_mcr_1 | 0.0015 |
| salmonella_enterica_serovar_typhi | global_efflux_pump | 0.005 |
| salmonella_enterica_serovar_typhi | global_porin_loss | 0.001 |
| salmonella_enterica_serovar_typhi | mutation_folate_pathway | 0.15 |
| salmonella_enterica_serovar_typhi | mutation_nitroreductase | 0.0015 |
| salmonella_enterica_serovar_typhi | enzyme_fos_a | 0.0015 |
| salmonella_enterica_serovar_typhi | mutation_rpo_b | 1.5e-4 |
| salmonella_enterica_serovar_typhi | protection_tet_m | 0.015 |
| salmonella_enterica_serovar_typhi | enzyme_aac_aph | 1e-4 |
| salmonella_enterica_serovar_typhi | efflux_tet_abc | 2e-6 |
| salmonella_enterica_serovar_typhi | efflux_mtr_cde | 0.01 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_esbl_ctx_m | 0.002 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_esbl_tem | 0.002 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_esbl_shv | 0.002 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_kpc | 5e-5 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_ndm_vim | 7e-5 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_oxa_48 | 5e-5 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_ampc_cmy | 0.0015 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_ampc_dha | 0.0015 |
| salmonella_enterica_serovar_paratyphi_a | mutation_gyra_primary | 0.05 |
| salmonella_enterica_serovar_paratyphi_a | mutation_gyra_parc_secondary | 0.04 |
| salmonella_enterica_serovar_paratyphi_a | protection_qnr | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_16s_rrmt | 0.15 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_cat | 0.02 |
| salmonella_enterica_serovar_paratyphi_a | efflux_acrab_tolc | 0.03 |
| salmonella_enterica_serovar_paratyphi_a | efflux_mexxy_oprm | 8e-4 |
| salmonella_enterica_serovar_paratyphi_a | modification_mcr_1 | 0.003 |
| salmonella_enterica_serovar_paratyphi_a | global_efflux_pump | 0.0045 |
| salmonella_enterica_serovar_paratyphi_a | global_porin_loss | 1e-4 |
| salmonella_enterica_serovar_paratyphi_a | mutation_folate_pathway | 0.03 |
| salmonella_enterica_serovar_paratyphi_a | mutation_nitroreductase | 1.5e-4 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_fos_a | 1.5e-4 |
| salmonella_enterica_serovar_paratyphi_a | mutation_rpo_b | 0.01 |
| salmonella_enterica_serovar_paratyphi_a | protection_tet_m | 0.02 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_aac_aph | 0.15 |
| salmonella_enterica_serovar_paratyphi_a | efflux_tet_abc | 0.007 |
| salmonella_enterica_serovar_paratyphi_a | mutation_pbp_mosaic | 0.001 |
| salmonella_enterica_serovar_paratyphi_a | efflux_mtr_cde | 5e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_esbl_ctx_m | 1.5e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_esbl_tem | 1.5e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_esbl_shv | 1.5e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_kpc | 1.5e-5 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_ndm_vim | 1.5e-5 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_oxa_48 | 1.5e-5 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_ampc_cmy | 1.5e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_ampc_dha | 1.5e-4 |
| invasive_non-typhoidal_salmonella_spp. | mutation_gyra_primary | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | mutation_gyra_parc_secondary | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | protection_qnr | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_16s_rrmt | 0.8 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_cat | 0.015 |
| invasive_non-typhoidal_salmonella_spp. | efflux_acrab_tolc | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | modification_mcr_1 | 0.03 |
| invasive_non-typhoidal_salmonella_spp. | global_efflux_pump | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | global_porin_loss | 1e-5 |
| invasive_non-typhoidal_salmonella_spp. | mutation_folate_pathway | 0.003 |
| invasive_non-typhoidal_salmonella_spp. | mutation_nitroreductase | 0.003 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_fos_a | 0.001 |
| invasive_non-typhoidal_salmonella_spp. | mutation_rpo_b | 0.02 |
| invasive_non-typhoidal_salmonella_spp. | protection_tet_m | 0.004 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_aac_aph | 0.8 |
| invasive_non-typhoidal_salmonella_spp. | efflux_tet_abc | 0.08 |
| invasive_non-typhoidal_salmonella_spp. | efflux_mtr_cde | 5e-5 |
| shigella_spp. | enzyme_esbl_ctx_m | 0.001 |
| shigella_spp. | enzyme_esbl_tem | 0.001 |
| shigella_spp. | enzyme_esbl_shv | 0.001 |
| shigella_spp. | enzyme_kpc | 2e-4 |
| shigella_spp. | enzyme_ndm_vim | 2e-4 |
| shigella_spp. | enzyme_oxa_48 | 2e-4 |
| shigella_spp. | enzyme_ampc_cmy | 0.001 |
| shigella_spp. | enzyme_ampc_dha | 0.001 |
| shigella_spp. | mutation_gyra_primary | 0.8 |
| shigella_spp. | mutation_gyra_parc_secondary | 0.8 |
| shigella_spp. | protection_qnr | 0.8 |
| shigella_spp. | enzyme_16s_rrmt | 0.9 |
| shigella_spp. | target_site_erm_b | 3 |
| shigella_spp. | enzyme_cat | 2.5e-4 |
| shigella_spp. | efflux_acrab_tolc | 0.003 |
| shigella_spp. | modification_mcr_1 | 3e-4 |
| shigella_spp. | global_efflux_pump | 3 |
| shigella_spp. | global_porin_loss | 3e-5 |
| shigella_spp. | mutation_folate_pathway | 0.003 |
| shigella_spp. | mutation_rpo_b | 0.04 |
| shigella_spp. | protection_tet_m | 0.03 |
| shigella_spp. | enzyme_aac_aph | 0.9 |
| shigella_spp. | mutation_23s_rrna | 0.9 |
| shigella_spp. | efflux_tet_abc | 0.03 |
| shigella_spp. | mutation_pbp_mosaic | 0.001 |
| shigella_spp. | efflux_mtr_cde | 0.001 |
| neisseria_gonorrhoeae | mutation_gyra_primary | 100 |
| neisseria_gonorrhoeae | mutation_gyra_parc_secondary | 100 |
| neisseria_gonorrhoeae | protection_qnr | 100 |
| neisseria_gonorrhoeae | enzyme_16s_rrmt | 100 |
| neisseria_gonorrhoeae | target_site_erm_b | 0.01 |
| neisseria_gonorrhoeae | target_site_cfr | 0.01 |
| neisseria_gonorrhoeae | enzyme_cat | 0.001 |
| neisseria_gonorrhoeae | modification_mcr_1 | 0.005 |
| neisseria_gonorrhoeae | global_efflux_pump | 100 |
| neisseria_gonorrhoeae | global_porin_loss | 0.001 |
| neisseria_gonorrhoeae | mutation_folate_pathway | 0.1 |
| neisseria_gonorrhoeae | mutation_nitroreductase | 0.1 |
| neisseria_gonorrhoeae | enzyme_fos_a | 3e-4 |
| neisseria_gonorrhoeae | mutation_rpo_b | 100 |
| neisseria_gonorrhoeae | protection_tet_m | 100 |
| neisseria_gonorrhoeae | enzyme_aac_aph | 100 |
| neisseria_gonorrhoeae | enzyme_bla_z | 0.02 |
| neisseria_gonorrhoeae | mutation_23s_rrna | 0.03 |
| neisseria_gonorrhoeae | efflux_tet_abc | 0.003 |
| neisseria_gonorrhoeae | mutation_pbp_mosaic | 0.02 |
| neisseria_gonorrhoeae | efflux_mtr_cde | 0.02 |
| streptococcus_pyogenes | mutation_gyra_primary | 3.8e-7 |
| streptococcus_pyogenes | mutation_gyra_parc_secondary | 7.5e-8 |
| streptococcus_pyogenes | target_site_erm_b | 7.5e-6 |
| streptococcus_pyogenes | target_site_cfr | 7.5e-9 |
| streptococcus_pyogenes | enzyme_cat | 3.8e-7 |
| streptococcus_pyogenes | global_efflux_pump | 7.5e-7 |
| streptococcus_pyogenes | mutation_folate_pathway | 3.8e-6 |
| streptococcus_pyogenes | mutation_mpr_f | 7.5e-9 |
| streptococcus_pyogenes | mutation_rpo_b | 7.5e-9 |
| streptococcus_pyogenes | protection_fus_b | 7.5e-9 |
| streptococcus_pyogenes | protection_tet_m | 7.5e-6 |
| streptococcus_pyogenes | mutation_23s_rrna | 5e-11 |
| streptococcus_agalactiae | target_site_pbp2a_meca | 1e-8 |
| streptococcus_agalactiae | target_site_van_a | 1e-8 |
| streptococcus_agalactiae | target_site_van_b | 1e-8 |
| streptococcus_agalactiae | mutation_gyra_primary | 1e-6 |
| streptococcus_agalactiae | mutation_gyra_parc_secondary | 1e-6 |
| streptococcus_agalactiae | target_site_erm_b | 5e-5 |
| streptococcus_agalactiae | target_site_cfr | 1e-8 |
| streptococcus_agalactiae | enzyme_cat | 3e-7 |
| streptococcus_agalactiae | global_efflux_pump | 1e-6 |
| streptococcus_agalactiae | mutation_folate_pathway | 1e-5 |
| streptococcus_agalactiae | mutation_mpr_f | 1e-7 |
| streptococcus_agalactiae | mutation_rpo_b | 1e-7 |
| streptococcus_agalactiae | protection_fus_b | 1e-7 |
| streptococcus_agalactiae | protection_tet_m | 0.001 |
| streptococcus_agalactiae | mutation_23s_rrna | 5e-11 |
| streptococcus_agalactiae | mutation_pbp_mosaic | 1e-6 |
| haemophilus_influenzae | enzyme_esbl_ctx_m | 1e-7 |
| haemophilus_influenzae | enzyme_esbl_tem | 1e-7 |
| haemophilus_influenzae | enzyme_esbl_shv | 1e-7 |
| haemophilus_influenzae | enzyme_kpc | 5e-8 |
| haemophilus_influenzae | enzyme_ndm_vim | 5e-8 |
| haemophilus_influenzae | enzyme_oxa_48 | 5e-8 |
| haemophilus_influenzae | enzyme_ampc_cmy | 1e-7 |
| haemophilus_influenzae | enzyme_ampc_dha | 1e-7 |
| haemophilus_influenzae | mutation_gyra_primary | 8e-4 |
| haemophilus_influenzae | mutation_gyra_parc_secondary | 8e-5 |
| haemophilus_influenzae | protection_qnr | 2e-4 |
| haemophilus_influenzae | enzyme_16s_rrmt | 0.1 |
| haemophilus_influenzae | target_site_erm_b | 7e-5 |
| haemophilus_influenzae | target_site_cfr | 5e-5 |
| haemophilus_influenzae | enzyme_cat | 1e-4 |
| haemophilus_influenzae | modification_mcr_1 | 2e-6 |
| haemophilus_influenzae | global_efflux_pump | 5e-4 |
| haemophilus_influenzae | global_porin_loss | 3e-7 |
| haemophilus_influenzae | mutation_folate_pathway | 0.002 |
| haemophilus_influenzae | mutation_nitroreductase | 2e-5 |
| haemophilus_influenzae | mutation_rpo_b | 0.015 |
| haemophilus_influenzae | protection_tet_m | 0.001 |
| haemophilus_influenzae | enzyme_aac_aph | 0.1 |
| haemophilus_influenzae | enzyme_bla_z | 1e-7 |
| haemophilus_influenzae | mutation_23s_rrna | 3e-6 |
| haemophilus_influenzae | mutation_pbp_mosaic | 1e-7 |
| haemophilus_influenzae | efflux_mtr_cde | 1e-7 |
| chlamydia_trachomatis | mutation_gyra_primary | 2e-7 |
| chlamydia_trachomatis | mutation_gyra_parc_secondary | 2e-7 |
| chlamydia_trachomatis | target_site_erm_b | 1e-7 |
| chlamydia_trachomatis | target_site_cfr | 2e-9 |
| chlamydia_trachomatis | enzyme_cat | 1.3e-9 |
| chlamydia_trachomatis | global_efflux_pump | 2e-8 |
| chlamydia_trachomatis | mutation_folate_pathway | 2e-9 |
| chlamydia_trachomatis | mutation_nitroreductase | 2e-9 |
| chlamydia_trachomatis | mutation_rpo_b | 2e-8 |
| chlamydia_trachomatis | protection_tet_m | 2e-7 |
| chlamydia_trachomatis | mutation_23s_rrna | 1e-10 |
| mycoplasma_genitalium | mutation_gyra_primary | 0.7 |
| mycoplasma_genitalium | mutation_gyra_parc_secondary | 0.7 |
| mycoplasma_genitalium | target_site_erm_b | 0.05 |
| mycoplasma_genitalium | target_site_cfr | 0.05 |
| mycoplasma_genitalium | enzyme_cat | 0.03 |
| mycoplasma_genitalium | global_efflux_pump | 0.2 |
| mycoplasma_genitalium | mutation_folate_pathway | 0.0035 |
| mycoplasma_genitalium | mutation_nitroreductase | 0.0035 |
| mycoplasma_genitalium | mutation_rpo_b | 0.035 |
| mycoplasma_genitalium | protection_tet_m | 0.2 |
| mycoplasma_genitalium | mutation_23s_rrna | 0.05 |
| vibrio_cholerae | enzyme_esbl_ctx_m | 3e-6 |
| vibrio_cholerae | enzyme_esbl_tem | 3e-6 |
| vibrio_cholerae | enzyme_esbl_shv | 1e-6 |
| vibrio_cholerae | enzyme_kpc | 1e-7 |
| vibrio_cholerae | enzyme_ndm_vim | 1e-6 |
| vibrio_cholerae | enzyme_oxa_48 | 1e-7 |
| vibrio_cholerae | enzyme_ampc_cmy | 1e-6 |
| vibrio_cholerae | enzyme_ampc_dha | 1e-6 |
| vibrio_cholerae | mutation_gyra_primary | 3e-5 |
| vibrio_cholerae | mutation_gyra_parc_secondary | 1.5e-5 |
| vibrio_cholerae | protection_qnr | 3e-6 |
| vibrio_cholerae | enzyme_16s_rrmt | 3e-8 |
| vibrio_cholerae | enzyme_cat | 1e-4 |
| vibrio_cholerae | efflux_acrab_tolc | 3e-6 |
| vibrio_cholerae | modification_mcr_1 | 3e-7 |
| vibrio_cholerae | global_efflux_pump | 3e-6 |
| vibrio_cholerae | global_porin_loss | 1.5e-6 |
| vibrio_cholerae | mutation_folate_pathway | 1.5e-4 |
| vibrio_cholerae | mutation_nitroreductase | 3e-7 |
| vibrio_cholerae | enzyme_fos_a | 3e-7 |
| vibrio_cholerae | mutation_rpo_b | 3e-7 |
| vibrio_cholerae | protection_tet_m | 5e-5 |
| vibrio_cholerae | enzyme_aac_aph | 5e-11 |
| vibrio_cholerae | mutation_23s_rrna | 1.5e-4 |
| vibrio_cholerae | efflux_tet_abc | 5e-11 |
| vibrio_cholerae | efflux_mtr_cde | 1e-9 |
| neisseria_meningitidis | enzyme_esbl_ctx_m | 1e-8 |
| neisseria_meningitidis | enzyme_esbl_tem | 1e-8 |
| neisseria_meningitidis | enzyme_esbl_shv | 1e-8 |
| neisseria_meningitidis | enzyme_ampc_cmy | 1e-8 |
| neisseria_meningitidis | enzyme_ampc_dha | 1e-8 |
| neisseria_meningitidis | mutation_gyra_primary | 3e-5 |
| neisseria_meningitidis | mutation_gyra_parc_secondary | 1e-5 |
| neisseria_meningitidis | protection_qnr | 3e-7 |
| neisseria_meningitidis | enzyme_16s_rrmt | 1e-7 |
| neisseria_meningitidis | target_site_erm_b | 1e-7 |
| neisseria_meningitidis | target_site_cfr | 2e-7 |
| neisseria_meningitidis | enzyme_cat | 2e-6 |
| neisseria_meningitidis | efflux_acrab_tolc | 3e-7 |
| neisseria_meningitidis | modification_mcr_1 | 3e-7 |
| neisseria_meningitidis | global_efflux_pump | 2e-6 |
| neisseria_meningitidis | global_porin_loss | 5e-7 |
| neisseria_meningitidis | mutation_folate_pathway | 5e-6 |
| neisseria_meningitidis | mutation_nitroreductase | 3e-6 |
| neisseria_meningitidis | mutation_rpo_b | 3e-6 |
| neisseria_meningitidis | protection_tet_m | 5e-6 |
| neisseria_meningitidis | mutation_23s_rrna | 1e-7 |
| neisseria_meningitidis | efflux_tet_abc | 5e-6 |
| neisseria_meningitidis | mutation_pbp_mosaic | 1e-6 |
| neisseria_meningitidis | efflux_mtr_cde | 1e-9 |
| listeria_monocytogenes | target_site_van_a | 3.8e-6 |
| listeria_monocytogenes | target_site_van_b | 3.8e-6 |
| listeria_monocytogenes | mutation_gyra_primary | 1.9e-4 |
| listeria_monocytogenes | mutation_gyra_parc_secondary | 3.8e-5 |
| listeria_monocytogenes | target_site_erm_b | 3.8e-4 |
| listeria_monocytogenes | target_site_cfr | 3.8e-6 |
| listeria_monocytogenes | enzyme_cat | 1.9e-4 |
| listeria_monocytogenes | global_efflux_pump | 1.9e-4 |
| listeria_monocytogenes | mutation_folate_pathway | 3.8e-4 |
| listeria_monocytogenes | mutation_mpr_f | 3.8e-5 |
| listeria_monocytogenes | mutation_rpo_b | 3.8e-5 |
| listeria_monocytogenes | protection_fus_b | 3.8e-5 |
| listeria_monocytogenes | protection_tet_m | 0.0019 |
| clostridioides_difficile | mutation_gyra_primary | 6e-5 |
| clostridioides_difficile | mutation_gyra_parc_secondary | 1e-5 |
| clostridioides_difficile | enzyme_16s_rrmt | 1e-7 |
| clostridioides_difficile | target_site_erm_b | 6e-5 |
| clostridioides_difficile | target_site_cfr | 1e-6 |
| clostridioides_difficile | enzyme_cat | 6e-6 |
| clostridioides_difficile | global_efflux_pump | 1e-5 |
| clostridioides_difficile | mutation_folate_pathway | 6e-6 |
| clostridioides_difficile | mutation_nitroreductase | 1.2e-4 |
| clostridioides_difficile | mutation_rpo_b | 6e-5 |
| clostridioides_difficile | protection_tet_m | 6e-5 |
| bacteroides_fragilis | enzyme_esbl_ctx_m | 3e-4 |
| bacteroides_fragilis | enzyme_esbl_tem | 3e-4 |
| bacteroides_fragilis | enzyme_esbl_shv | 3e-4 |
| bacteroides_fragilis | enzyme_kpc | 1e-5 |
| bacteroides_fragilis | enzyme_ndm_vim | 1e-5 |
| bacteroides_fragilis | enzyme_oxa_48 | 1e-5 |
| bacteroides_fragilis | enzyme_ampc_cmy | 3e-4 |
| bacteroides_fragilis | enzyme_ampc_dha | 3e-4 |
| bacteroides_fragilis | mutation_gyra_primary | 0.03 |
| bacteroides_fragilis | mutation_gyra_parc_secondary | 0.03 |
| bacteroides_fragilis | protection_qnr | 0.03 |
| bacteroides_fragilis | enzyme_16s_rrmt | 3 |
| bacteroides_fragilis | target_site_erm_b | 0.01 |
| bacteroides_fragilis | target_site_cfr | 1e-4 |
| bacteroides_fragilis | enzyme_cat | 1e-5 |
| bacteroides_fragilis | efflux_acrab_tolc | 0.003 |
| bacteroides_fragilis | modification_mcr_1 | 5e-4 |
| bacteroides_fragilis | global_efflux_pump | 0.003 |
| bacteroides_fragilis | global_porin_loss | 5e-6 |
| bacteroides_fragilis | mutation_folate_pathway | 0.03 |
| bacteroides_fragilis | mutation_nitroreductase | 1e-4 |
| bacteroides_fragilis | mutation_rpo_b | 0.001 |
| bacteroides_fragilis | protection_tet_m | 0.01 |
| bacteroides_fragilis | enzyme_aac_aph | 3 |
| bacteroides_fragilis | mutation_pbp_mosaic | 2e-4 |
| campylobacter_jejuni | mutation_gyra_primary | 0.03 |
| campylobacter_jejuni | mutation_gyra_parc_secondary | 0.03 |
| campylobacter_jejuni | target_site_erm_b | 1e-5 |
| campylobacter_jejuni | target_site_cfr | 1e-5 |
| campylobacter_jejuni | enzyme_cat | 6e-4 |
| campylobacter_jejuni | global_efflux_pump | 0.03 |
| campylobacter_jejuni | global_porin_loss | 2e-5 |
| campylobacter_jejuni | mutation_folate_pathway | 0.6 |
| campylobacter_jejuni | mutation_rpo_b | 0.6 |
| campylobacter_jejuni | protection_tet_m | 0.01 |
| campylobacter_jejuni | enzyme_aac_aph | 0.1 |
| campylobacter_jejuni | mutation_23s_rrna | 3e-4 |
| campylobacter_jejuni | efflux_tet_abc | 0.001 |
| campylobacter_jejuni | efflux_mtr_cde | 2e-4 |
| enterobacter_cloacae | enzyme_esbl_ctx_m | 5e-5 |
| enterobacter_cloacae | enzyme_esbl_tem | 5e-5 |
| enterobacter_cloacae | enzyme_esbl_shv | 5e-5 |
| enterobacter_cloacae | enzyme_kpc | 5e-6 |
| enterobacter_cloacae | enzyme_ndm_vim | 5e-6 |
| enterobacter_cloacae | enzyme_oxa_48 | 5e-6 |
| enterobacter_cloacae | enzyme_ampc_cmy | 5e-5 |
| enterobacter_cloacae | enzyme_ampc_dha | 5e-5 |
| enterobacter_cloacae | mutation_gyra_primary | 0.003 |
| enterobacter_cloacae | mutation_gyra_parc_secondary | 0.003 |
| enterobacter_cloacae | protection_qnr | 0.003 |
| enterobacter_cloacae | enzyme_16s_rrmt | 0.005 |
| enterobacter_cloacae | enzyme_cat | 3e-4 |
| enterobacter_cloacae | efflux_acrab_tolc | 0.003 |
| enterobacter_cloacae | modification_mcr_1 | 0.01 |
| enterobacter_cloacae | global_porin_loss | 1e-4 |
| enterobacter_cloacae | mutation_folate_pathway | 0.003 |
| enterobacter_cloacae | mutation_nitroreductase | 0.002 |
| enterobacter_cloacae | enzyme_fos_a | 0.002 |
| enterobacter_cloacae | mutation_rpo_b | 0.02 |
| enterobacter_cloacae | protection_tet_m | 0.01 |
| enterobacter_cloacae | enzyme_aac_aph | 0.005 |
| enterobacter_cloacae | efflux_tet_abc | 0.01 |
| enterobacter_cloacae | mutation_pbp_mosaic | 5e-5 |
| yersinia_enterocolitica | enzyme_esbl_ctx_m | 3e-10 |
| yersinia_enterocolitica | enzyme_esbl_tem | 3e-10 |
| yersinia_enterocolitica | enzyme_esbl_shv | 1e-10 |
| yersinia_enterocolitica | enzyme_kpc | 3e-11 |
| yersinia_enterocolitica | enzyme_ndm_vim | 3e-11 |
| yersinia_enterocolitica | enzyme_oxa_48 | 3e-11 |
| yersinia_enterocolitica | enzyme_ampc_cmy | 3e-10 |
| yersinia_enterocolitica | enzyme_ampc_dha | 3e-10 |
| yersinia_enterocolitica | mutation_gyra_primary | 3e-10 |
| yersinia_enterocolitica | mutation_gyra_parc_secondary | 3e-10 |
| yersinia_enterocolitica | protection_qnr | 3e-10 |
| yersinia_enterocolitica | enzyme_16s_rrmt | 3e-11 |
| yersinia_enterocolitica | enzyme_cat | 3e-10 |
| yersinia_enterocolitica | efflux_acrab_tolc | 3e-10 |
| yersinia_enterocolitica | modification_mcr_1 | 3e-10 |
| yersinia_enterocolitica | global_efflux_pump | 3e-10 |
| yersinia_enterocolitica | global_porin_loss | 1e-10 |
| yersinia_enterocolitica | mutation_folate_pathway | 3e-10 |
| yersinia_enterocolitica | mutation_nitroreductase | 1e-9 |
| yersinia_enterocolitica | enzyme_fos_a | 3e-10 |
| yersinia_enterocolitica | mutation_rpo_b | 3e-11 |
| yersinia_enterocolitica | protection_tet_m | 3e-10 |
| yersinia_enterocolitica | enzyme_aac_aph | 1.5e-10 |
| yersinia_enterocolitica | efflux_tet_abc | 1.5e-10 |
| yersinia_enterocolitica | efflux_mtr_cde | 3e-10 |
| moraxella_catarrhalis | enzyme_esbl_ctx_m | 2e-7 |
| moraxella_catarrhalis | enzyme_esbl_tem | 5e-7 |
| moraxella_catarrhalis | enzyme_esbl_shv | 2e-8 |
| moraxella_catarrhalis | enzyme_ampc_cmy | 5e-7 |
| moraxella_catarrhalis | enzyme_ampc_dha | 2e-7 |
| moraxella_catarrhalis | mutation_gyra_primary | 1e-5 |
| moraxella_catarrhalis | mutation_gyra_parc_secondary | 1e-6 |
| moraxella_catarrhalis | protection_qnr | 1e-6 |
| moraxella_catarrhalis | enzyme_16s_rrmt | 1e-7 |
| moraxella_catarrhalis | target_site_erm_b | 1e-5 |
| moraxella_catarrhalis | target_site_cfr | 1e-7 |
| moraxella_catarrhalis | enzyme_cat | 1e-6 |
| moraxella_catarrhalis | efflux_acrab_tolc | 1e-6 |
| moraxella_catarrhalis | modification_mcr_1 | 2e-7 |
| moraxella_catarrhalis | global_efflux_pump | 1e-5 |
| moraxella_catarrhalis | global_porin_loss | 2e-6 |
| moraxella_catarrhalis | mutation_folate_pathway | 1e-5 |
| moraxella_catarrhalis | mutation_nitroreductase | 2e-6 |
| moraxella_catarrhalis | mutation_rpo_b | 2e-6 |
| moraxella_catarrhalis | protection_tet_m | 1e-4 |
| moraxella_catarrhalis | mutation_pbp_mosaic | 5e-7 |
| moraxella_catarrhalis | efflux_mtr_cde | 5e-7 |
| treponema_pallidum | mutation_gyra_primary | 1.5e-6 |
| treponema_pallidum | mutation_gyra_parc_secondary | 7.5e-7 |
| treponema_pallidum | enzyme_cat | 1.5e-7 |
| treponema_pallidum | global_efflux_pump | 1.5e-7 |
| treponema_pallidum | mutation_folate_pathway | 1.5e-7 |
| treponema_pallidum | mutation_rpo_b | 1.5e-7 |
| treponema_pallidum | protection_tet_m | 1.5e-6 |
| bordetella_pertussis | enzyme_esbl_ctx_m | 1e-10 |
| bordetella_pertussis | enzyme_esbl_tem | 1e-10 |
| bordetella_pertussis | enzyme_esbl_shv | 1e-10 |
| bordetella_pertussis | enzyme_ampc_cmy | 1e-10 |
| bordetella_pertussis | enzyme_ampc_dha | 1e-10 |
| bordetella_pertussis | mutation_gyra_primary | 2e-12 |
| bordetella_pertussis | mutation_gyra_parc_secondary | 1e-12 |
| bordetella_pertussis | enzyme_16s_rrmt | 1e-10 |
| bordetella_pertussis | target_site_cfr | 5e-11 |
| bordetella_pertussis | enzyme_cat | 3.8e-12 |
| bordetella_pertussis | efflux_acrab_tolc | 1e-12 |
| bordetella_pertussis | global_efflux_pump | 5e-12 |
| bordetella_pertussis | global_porin_loss | 1e-12 |
| bordetella_pertussis | mutation_folate_pathway | 1e-11 |
| bordetella_pertussis | mutation_nitroreductase | 1e-10 |
| bordetella_pertussis | mutation_rpo_b | 1e-12 |
| bordetella_pertussis | protection_tet_m | 1e-11 |
| bordetella_pertussis | efflux_mtr_cde | 2e-12 |
| helicobacter_pylori | mutation_gyra_primary | 2e-4 |
| helicobacter_pylori | mutation_gyra_parc_secondary | 2e-4 |
| helicobacter_pylori | target_site_erm_b | 0.001 |
| helicobacter_pylori | target_site_cfr | 3e-4 |
| helicobacter_pylori | global_efflux_pump | 2e-4 |
| helicobacter_pylori | mutation_folate_pathway | 0.01 |
| helicobacter_pylori | mutation_nitroreductase | 0.03 |
| helicobacter_pylori | mutation_rpo_b | 0.01 |
| helicobacter_pylori | protection_tet_m | 2e-4 |
| helicobacter_pylori | enzyme_bla_z | 1e-5 |
| helicobacter_pylori | mutation_23s_rrna | 2e-4 |
| helicobacter_pylori | mutation_pbp_mosaic | 1e-5 |
| mycoplasma_pneumoniae | mutation_gyra_primary | 3e-7 |
| mycoplasma_pneumoniae | mutation_gyra_parc_secondary | 1.5e-7 |
| mycoplasma_pneumoniae | target_site_erm_b | 1.5e-5 |
| mycoplasma_pneumoniae | target_site_cfr | 3e-9 |
| mycoplasma_pneumoniae | enzyme_cat | 3e-9 |
| mycoplasma_pneumoniae | global_efflux_pump | 1.5e-7 |
| mycoplasma_pneumoniae | mutation_folate_pathway | 3e-9 |
| mycoplasma_pneumoniae | mutation_nitroreductase | 3e-9 |
| mycoplasma_pneumoniae | mutation_rpo_b | 3e-8 |
| mycoplasma_pneumoniae | protection_tet_m | 3e-7 |
| mycoplasma_pneumoniae | mutation_23s_rrna | 5e-9 |
| legionella_pneumophila | enzyme_esbl_ctx_m | 3e-7 |
| legionella_pneumophila | enzyme_esbl_tem | 3e-7 |
| legionella_pneumophila | enzyme_esbl_shv | 3e-7 |
| legionella_pneumophila | enzyme_ampc_cmy | 3e-7 |
| legionella_pneumophila | enzyme_ampc_dha | 3e-7 |
| legionella_pneumophila | mutation_gyra_primary | 3e-5 |
| legionella_pneumophila | mutation_gyra_parc_secondary | 3e-5 |
| legionella_pneumophila | protection_qnr | 1e-7 |
| legionella_pneumophila | enzyme_16s_rrmt | 3e-7 |
| legionella_pneumophila | target_site_erm_b | 3e-5 |
| legionella_pneumophila | target_site_cfr | 3e-7 |
| legionella_pneumophila | enzyme_cat | 3e-6 |
| legionella_pneumophila | efflux_acrab_tolc | 3e-7 |
| legionella_pneumophila | modification_mcr_1 | 3e-7 |
| legionella_pneumophila | global_efflux_pump | 3e-5 |
| legionella_pneumophila | global_porin_loss | 1e-6 |
| legionella_pneumophila | mutation_folate_pathway | 3e-6 |
| legionella_pneumophila | mutation_nitroreductase | 3e-7 |
| legionella_pneumophila | mutation_rpo_b | 3e-6 |
| legionella_pneumophila | protection_tet_m | 3e-5 |
| legionella_pneumophila | mutation_23s_rrna | 5e-11 |
| burkholderia_cepacia_complex | enzyme_esbl_ctx_m | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_esbl_tem | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_esbl_shv | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_kpc | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_ndm_vim | 3.7e-5 |
| burkholderia_cepacia_complex | enzyme_oxa_48 | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_ampc_cmy | 3.8e-4 |
| burkholderia_cepacia_complex | enzyme_ampc_dha | 7.5e-5 |
| burkholderia_cepacia_complex | mutation_gyra_primary | 3.8e-4 |
| burkholderia_cepacia_complex | mutation_gyra_parc_secondary | 7.5e-5 |
| burkholderia_cepacia_complex | protection_qnr | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_16s_rrmt | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_cat | 3.7e-5 |
| burkholderia_cepacia_complex | modification_mcr_1 | 7.5e-6 |
| burkholderia_cepacia_complex | global_efflux_pump | 7.5e-4 |
| burkholderia_cepacia_complex | global_porin_loss | 3.8e-4 |
| burkholderia_cepacia_complex | mutation_folate_pathway | 7.5e-5 |
| burkholderia_cepacia_complex | enzyme_fos_a | 3.7e-5 |
| burkholderia_cepacia_complex | mutation_rpo_b | 7.5e-6 |
| burkholderia_cepacia_complex | protection_tet_m | 3.7e-5 |
| burkholderia_cepacia_complex | enzyme_aac_aph | 5e-11 |
| burkholderia_cepacia_complex | efflux_tet_abc | 5e-11 |

### B.11 Horizontal Gene Transfer Matrix

Per-day probability of horizontal gene transfer of resistance between co-colonising bacterial species. Only non-zero entries shown.

See: [┬º9.1 Transfer compatibility](#91-transfer-compatibility), [┬º9.2 The HGT process](#92-the-hgt-process).

| Donor | Recipient | Probability/day |
| --- | ---: | ---: |
| acinetobacter_baumannii | citrobacter_spp. | 1e-10 |
| acinetobacter_baumannii | enterobacter_spp. | 1e-10 |
| acinetobacter_baumannii | escherichia_coli | 1e-10 |
| acinetobacter_baumannii | klebsiella_pneumoniae | 1e-10 |
| acinetobacter_baumannii | morganella_spp. | 1e-10 |
| acinetobacter_baumannii | proteus_spp. | 1e-10 |
| acinetobacter_baumannii | serratia_spp. | 1e-10 |
| acinetobacter_baumannii | p_stuartii | 1e-10 |
| acinetobacter_baumannii | pseudomonas_aeruginosa | 1e-9 |
| acinetobacter_baumannii | stenotrophomonas_maltophilia | 1e-9 |
| acinetobacter_baumannii | salmonella_enterica_serovar_typhi | 1e-10 |
| acinetobacter_baumannii | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| acinetobacter_baumannii | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| acinetobacter_baumannii | shigella_spp. | 1e-10 |
| acinetobacter_baumannii | neisseria_gonorrhoeae | 3e-11 |
| acinetobacter_baumannii | haemophilus_influenzae | 3e-11 |
| acinetobacter_baumannii | chlamydia_trachomatis | 3e-11 |
| acinetobacter_baumannii | mycoplasma_genitalium | 3e-11 |
| acinetobacter_baumannii | vibrio_cholerae | 1e-10 |
| acinetobacter_baumannii | neisseria_meningitidis | 3e-11 |
| acinetobacter_baumannii | clostridioides_difficile | 3e-11 |
| acinetobacter_baumannii | bacteroides_fragilis | 3e-11 |
| acinetobacter_baumannii | enterobacter_cloacae | 1e-10 |
| acinetobacter_baumannii | yersinia_enterocolitica | 1e-10 |
| acinetobacter_baumannii | moraxella_catarrhalis | 3e-11 |
| acinetobacter_baumannii | bordetella_pertussis | 3e-11 |
| acinetobacter_baumannii | mycoplasma_pneumoniae | 3e-11 |
| acinetobacter_baumannii | legionella_pneumophila | 3e-11 |
| acinetobacter_baumannii | burkholderia_cepacia_complex | 1e-9 |
| citrobacter_spp. | acinetobacter_baumannii | 1e-10 |
| citrobacter_spp. | enterobacter_spp. | 1e-9 |
| citrobacter_spp. | escherichia_coli | 1e-9 |
| citrobacter_spp. | klebsiella_pneumoniae | 1e-9 |
| citrobacter_spp. | morganella_spp. | 1e-9 |
| citrobacter_spp. | proteus_spp. | 1e-9 |
| citrobacter_spp. | serratia_spp. | 1e-9 |
| citrobacter_spp. | p_stuartii | 1e-9 |
| citrobacter_spp. | pseudomonas_aeruginosa | 1e-10 |
| citrobacter_spp. | stenotrophomonas_maltophilia | 1e-10 |
| citrobacter_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| citrobacter_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| citrobacter_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| citrobacter_spp. | shigella_spp. | 1e-9 |
| citrobacter_spp. | neisseria_gonorrhoeae | 3e-11 |
| citrobacter_spp. | haemophilus_influenzae | 3e-11 |
| citrobacter_spp. | chlamydia_trachomatis | 3e-11 |
| citrobacter_spp. | mycoplasma_genitalium | 3e-11 |
| citrobacter_spp. | vibrio_cholerae | 1e-10 |
| citrobacter_spp. | neisseria_meningitidis | 3e-11 |
| citrobacter_spp. | clostridioides_difficile | 3e-11 |
| citrobacter_spp. | bacteroides_fragilis | 3e-11 |
| citrobacter_spp. | enterobacter_cloacae | 1e-9 |
| citrobacter_spp. | yersinia_enterocolitica | 1e-9 |
| citrobacter_spp. | moraxella_catarrhalis | 3e-11 |
| citrobacter_spp. | bordetella_pertussis | 3e-11 |
| citrobacter_spp. | mycoplasma_pneumoniae | 3e-11 |
| citrobacter_spp. | legionella_pneumophila | 3e-11 |
| citrobacter_spp. | burkholderia_cepacia_complex | 1e-10 |
| enterobacter_spp. | acinetobacter_baumannii | 1e-10 |
| enterobacter_spp. | citrobacter_spp. | 1e-9 |
| enterobacter_spp. | escherichia_coli | 1e-9 |
| enterobacter_spp. | klebsiella_pneumoniae | 1e-9 |
| enterobacter_spp. | morganella_spp. | 1e-9 |
| enterobacter_spp. | proteus_spp. | 1e-9 |
| enterobacter_spp. | serratia_spp. | 1e-9 |
| enterobacter_spp. | p_stuartii | 1e-9 |
| enterobacter_spp. | pseudomonas_aeruginosa | 1e-10 |
| enterobacter_spp. | stenotrophomonas_maltophilia | 1e-10 |
| enterobacter_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| enterobacter_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| enterobacter_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| enterobacter_spp. | shigella_spp. | 1e-9 |
| enterobacter_spp. | neisseria_gonorrhoeae | 3e-11 |
| enterobacter_spp. | haemophilus_influenzae | 3e-11 |
| enterobacter_spp. | chlamydia_trachomatis | 3e-11 |
| enterobacter_spp. | mycoplasma_genitalium | 3e-11 |
| enterobacter_spp. | vibrio_cholerae | 1e-10 |
| enterobacter_spp. | neisseria_meningitidis | 3e-11 |
| enterobacter_spp. | clostridioides_difficile | 3e-11 |
| enterobacter_spp. | bacteroides_fragilis | 3e-11 |
| enterobacter_spp. | enterobacter_cloacae | 1e-9 |
| enterobacter_spp. | yersinia_enterocolitica | 1e-9 |
| enterobacter_spp. | moraxella_catarrhalis | 3e-11 |
| enterobacter_spp. | bordetella_pertussis | 3e-11 |
| enterobacter_spp. | mycoplasma_pneumoniae | 3e-11 |
| enterobacter_spp. | legionella_pneumophila | 3e-11 |
| enterobacter_spp. | burkholderia_cepacia_complex | 1e-10 |
| enterococcus_faecalis | enterococcus_faecium | 1e-9 |
| enterococcus_faecalis | staphylococcus_aureus | 1e-9 |
| enterococcus_faecalis | staphylococcus_epidermidis | 1e-9 |
| enterococcus_faecalis | streptococcus_pneumoniae | 1e-9 |
| enterococcus_faecalis | streptococcus_pyogenes | 1e-9 |
| enterococcus_faecalis | streptococcus_agalactiae | 1e-9 |
| enterococcus_faecalis | listeria_monocytogenes | 1e-9 |
| enterococcus_faecium | enterococcus_faecalis | 1e-9 |
| enterococcus_faecium | staphylococcus_aureus | 1e-9 |
| enterococcus_faecium | staphylococcus_epidermidis | 1e-9 |
| enterococcus_faecium | streptococcus_pneumoniae | 1e-9 |
| enterococcus_faecium | streptococcus_pyogenes | 1e-9 |
| enterococcus_faecium | streptococcus_agalactiae | 1e-9 |
| enterococcus_faecium | listeria_monocytogenes | 1e-9 |
| escherichia_coli | acinetobacter_baumannii | 1e-10 |
| escherichia_coli | citrobacter_spp. | 1e-9 |
| escherichia_coli | enterobacter_spp. | 1e-9 |
| escherichia_coli | klebsiella_pneumoniae | 1e-9 |
| escherichia_coli | morganella_spp. | 1e-9 |
| escherichia_coli | proteus_spp. | 1e-9 |
| escherichia_coli | serratia_spp. | 1e-9 |
| escherichia_coli | p_stuartii | 1e-9 |
| escherichia_coli | pseudomonas_aeruginosa | 1e-10 |
| escherichia_coli | stenotrophomonas_maltophilia | 1e-10 |
| escherichia_coli | salmonella_enterica_serovar_typhi | 1e-9 |
| escherichia_coli | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| escherichia_coli | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| escherichia_coli | shigella_spp. | 1e-9 |
| escherichia_coli | neisseria_gonorrhoeae | 3e-11 |
| escherichia_coli | haemophilus_influenzae | 3e-11 |
| escherichia_coli | chlamydia_trachomatis | 3e-11 |
| escherichia_coli | mycoplasma_genitalium | 3e-11 |
| escherichia_coli | vibrio_cholerae | 1e-10 |
| escherichia_coli | neisseria_meningitidis | 3e-11 |
| escherichia_coli | clostridioides_difficile | 3e-11 |
| escherichia_coli | bacteroides_fragilis | 3e-11 |
| escherichia_coli | enterobacter_cloacae | 1e-9 |
| escherichia_coli | yersinia_enterocolitica | 1e-9 |
| escherichia_coli | moraxella_catarrhalis | 3e-11 |
| escherichia_coli | bordetella_pertussis | 3e-11 |
| escherichia_coli | mycoplasma_pneumoniae | 3e-11 |
| escherichia_coli | legionella_pneumophila | 3e-11 |
| escherichia_coli | burkholderia_cepacia_complex | 1e-10 |
| klebsiella_pneumoniae | acinetobacter_baumannii | 1e-10 |
| klebsiella_pneumoniae | citrobacter_spp. | 1e-9 |
| klebsiella_pneumoniae | enterobacter_spp. | 1e-9 |
| klebsiella_pneumoniae | escherichia_coli | 1e-9 |
| klebsiella_pneumoniae | morganella_spp. | 1e-9 |
| klebsiella_pneumoniae | proteus_spp. | 1e-9 |
| klebsiella_pneumoniae | serratia_spp. | 1e-9 |
| klebsiella_pneumoniae | p_stuartii | 1e-9 |
| klebsiella_pneumoniae | pseudomonas_aeruginosa | 1e-10 |
| klebsiella_pneumoniae | stenotrophomonas_maltophilia | 1e-10 |
| klebsiella_pneumoniae | salmonella_enterica_serovar_typhi | 1e-9 |
| klebsiella_pneumoniae | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| klebsiella_pneumoniae | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| klebsiella_pneumoniae | shigella_spp. | 1e-9 |
| klebsiella_pneumoniae | neisseria_gonorrhoeae | 3e-11 |
| klebsiella_pneumoniae | haemophilus_influenzae | 3e-11 |
| klebsiella_pneumoniae | chlamydia_trachomatis | 3e-11 |
| klebsiella_pneumoniae | mycoplasma_genitalium | 3e-11 |
| klebsiella_pneumoniae | vibrio_cholerae | 1e-10 |
| klebsiella_pneumoniae | neisseria_meningitidis | 3e-11 |
| klebsiella_pneumoniae | clostridioides_difficile | 3e-11 |
| klebsiella_pneumoniae | bacteroides_fragilis | 3e-11 |
| klebsiella_pneumoniae | enterobacter_cloacae | 1e-9 |
| klebsiella_pneumoniae | yersinia_enterocolitica | 1e-9 |
| klebsiella_pneumoniae | moraxella_catarrhalis | 3e-11 |
| klebsiella_pneumoniae | bordetella_pertussis | 3e-11 |
| klebsiella_pneumoniae | mycoplasma_pneumoniae | 3e-11 |
| klebsiella_pneumoniae | legionella_pneumophila | 3e-11 |
| klebsiella_pneumoniae | burkholderia_cepacia_complex | 1e-10 |
| morganella_spp. | acinetobacter_baumannii | 1e-10 |
| morganella_spp. | citrobacter_spp. | 1e-9 |
| morganella_spp. | enterobacter_spp. | 1e-9 |
| morganella_spp. | escherichia_coli | 1e-9 |
| morganella_spp. | klebsiella_pneumoniae | 1e-9 |
| morganella_spp. | proteus_spp. | 1e-9 |
| morganella_spp. | serratia_spp. | 1e-9 |
| morganella_spp. | p_stuartii | 1e-9 |
| morganella_spp. | pseudomonas_aeruginosa | 1e-10 |
| morganella_spp. | stenotrophomonas_maltophilia | 1e-10 |
| morganella_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| morganella_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| morganella_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| morganella_spp. | shigella_spp. | 1e-9 |
| morganella_spp. | neisseria_gonorrhoeae | 3e-11 |
| morganella_spp. | haemophilus_influenzae | 3e-11 |
| morganella_spp. | chlamydia_trachomatis | 3e-11 |
| morganella_spp. | mycoplasma_genitalium | 3e-11 |
| morganella_spp. | vibrio_cholerae | 1e-10 |
| morganella_spp. | neisseria_meningitidis | 3e-11 |
| morganella_spp. | clostridioides_difficile | 3e-11 |
| morganella_spp. | bacteroides_fragilis | 3e-11 |
| morganella_spp. | enterobacter_cloacae | 1e-9 |
| morganella_spp. | yersinia_enterocolitica | 1e-9 |
| morganella_spp. | moraxella_catarrhalis | 3e-11 |
| morganella_spp. | bordetella_pertussis | 3e-11 |
| morganella_spp. | mycoplasma_pneumoniae | 3e-11 |
| morganella_spp. | legionella_pneumophila | 3e-11 |
| morganella_spp. | burkholderia_cepacia_complex | 1e-10 |
| proteus_spp. | acinetobacter_baumannii | 1e-10 |
| proteus_spp. | citrobacter_spp. | 1e-9 |
| proteus_spp. | enterobacter_spp. | 1e-9 |
| proteus_spp. | escherichia_coli | 1e-9 |
| proteus_spp. | klebsiella_pneumoniae | 1e-9 |
| proteus_spp. | morganella_spp. | 1e-9 |
| proteus_spp. | serratia_spp. | 1e-9 |
| proteus_spp. | p_stuartii | 1e-9 |
| proteus_spp. | pseudomonas_aeruginosa | 1e-10 |
| proteus_spp. | stenotrophomonas_maltophilia | 1e-10 |
| proteus_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| proteus_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| proteus_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| proteus_spp. | shigella_spp. | 1e-9 |
| proteus_spp. | neisseria_gonorrhoeae | 3e-11 |
| proteus_spp. | haemophilus_influenzae | 3e-11 |
| proteus_spp. | chlamydia_trachomatis | 3e-11 |
| proteus_spp. | mycoplasma_genitalium | 3e-11 |
| proteus_spp. | vibrio_cholerae | 1e-10 |
| proteus_spp. | neisseria_meningitidis | 3e-11 |
| proteus_spp. | clostridioides_difficile | 3e-11 |
| proteus_spp. | bacteroides_fragilis | 3e-11 |
| proteus_spp. | enterobacter_cloacae | 1e-9 |
| proteus_spp. | yersinia_enterocolitica | 1e-9 |
| proteus_spp. | moraxella_catarrhalis | 3e-11 |
| proteus_spp. | bordetella_pertussis | 3e-11 |
| proteus_spp. | mycoplasma_pneumoniae | 3e-11 |
| proteus_spp. | legionella_pneumophila | 3e-11 |
| proteus_spp. | burkholderia_cepacia_complex | 1e-10 |
| serratia_spp. | acinetobacter_baumannii | 1e-10 |
| serratia_spp. | citrobacter_spp. | 1e-9 |
| serratia_spp. | enterobacter_spp. | 1e-9 |
| serratia_spp. | escherichia_coli | 1e-9 |
| serratia_spp. | klebsiella_pneumoniae | 1e-9 |
| serratia_spp. | morganella_spp. | 1e-9 |
| serratia_spp. | proteus_spp. | 1e-9 |
| serratia_spp. | p_stuartii | 1e-9 |
| serratia_spp. | pseudomonas_aeruginosa | 1e-10 |
| serratia_spp. | stenotrophomonas_maltophilia | 1e-10 |
| serratia_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| serratia_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| serratia_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| serratia_spp. | shigella_spp. | 1e-9 |
| serratia_spp. | neisseria_gonorrhoeae | 3e-11 |
| serratia_spp. | haemophilus_influenzae | 3e-11 |
| serratia_spp. | chlamydia_trachomatis | 3e-11 |
| serratia_spp. | mycoplasma_genitalium | 3e-11 |
| serratia_spp. | vibrio_cholerae | 1e-10 |
| serratia_spp. | neisseria_meningitidis | 3e-11 |
| serratia_spp. | clostridioides_difficile | 3e-11 |
| serratia_spp. | bacteroides_fragilis | 3e-11 |
| serratia_spp. | enterobacter_cloacae | 1e-9 |
| serratia_spp. | yersinia_enterocolitica | 1e-9 |
| serratia_spp. | moraxella_catarrhalis | 3e-11 |
| serratia_spp. | bordetella_pertussis | 3e-11 |
| serratia_spp. | mycoplasma_pneumoniae | 3e-11 |
| serratia_spp. | legionella_pneumophila | 3e-11 |
| serratia_spp. | burkholderia_cepacia_complex | 1e-10 |
| p_stuartii | acinetobacter_baumannii | 1e-10 |
| p_stuartii | citrobacter_spp. | 1e-9 |
| p_stuartii | enterobacter_spp. | 1e-9 |
| p_stuartii | escherichia_coli | 1e-9 |
| p_stuartii | klebsiella_pneumoniae | 1e-9 |
| p_stuartii | morganella_spp. | 1e-9 |
| p_stuartii | proteus_spp. | 1e-9 |
| p_stuartii | serratia_spp. | 1e-9 |
| p_stuartii | pseudomonas_aeruginosa | 1e-10 |
| p_stuartii | stenotrophomonas_maltophilia | 1e-10 |
| p_stuartii | salmonella_enterica_serovar_typhi | 1e-9 |
| p_stuartii | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| p_stuartii | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| p_stuartii | shigella_spp. | 1e-9 |
| p_stuartii | neisseria_gonorrhoeae | 3e-11 |
| p_stuartii | haemophilus_influenzae | 3e-11 |
| p_stuartii | chlamydia_trachomatis | 3e-11 |
| p_stuartii | mycoplasma_genitalium | 3e-11 |
| p_stuartii | vibrio_cholerae | 1e-10 |
| p_stuartii | neisseria_meningitidis | 3e-11 |
| p_stuartii | clostridioides_difficile | 3e-11 |
| p_stuartii | bacteroides_fragilis | 3e-11 |
| p_stuartii | enterobacter_cloacae | 1e-9 |
| p_stuartii | yersinia_enterocolitica | 1e-9 |
| p_stuartii | moraxella_catarrhalis | 3e-11 |
| p_stuartii | bordetella_pertussis | 3e-11 |
| p_stuartii | mycoplasma_pneumoniae | 3e-11 |
| p_stuartii | legionella_pneumophila | 3e-11 |
| p_stuartii | burkholderia_cepacia_complex | 1e-10 |
| pseudomonas_aeruginosa | acinetobacter_baumannii | 1e-9 |
| pseudomonas_aeruginosa | citrobacter_spp. | 1e-10 |
| pseudomonas_aeruginosa | enterobacter_spp. | 1e-10 |
| pseudomonas_aeruginosa | escherichia_coli | 1e-10 |
| pseudomonas_aeruginosa | klebsiella_pneumoniae | 1e-10 |
| pseudomonas_aeruginosa | morganella_spp. | 1e-10 |
| pseudomonas_aeruginosa | proteus_spp. | 1e-10 |
| pseudomonas_aeruginosa | serratia_spp. | 1e-10 |
| pseudomonas_aeruginosa | p_stuartii | 1e-10 |
| pseudomonas_aeruginosa | stenotrophomonas_maltophilia | 1e-9 |
| pseudomonas_aeruginosa | salmonella_enterica_serovar_typhi | 1e-10 |
| pseudomonas_aeruginosa | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| pseudomonas_aeruginosa | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| pseudomonas_aeruginosa | shigella_spp. | 1e-10 |
| pseudomonas_aeruginosa | neisseria_gonorrhoeae | 3e-11 |
| pseudomonas_aeruginosa | haemophilus_influenzae | 3e-11 |
| pseudomonas_aeruginosa | chlamydia_trachomatis | 3e-11 |
| pseudomonas_aeruginosa | mycoplasma_genitalium | 3e-11 |
| pseudomonas_aeruginosa | vibrio_cholerae | 1e-10 |
| pseudomonas_aeruginosa | neisseria_meningitidis | 3e-11 |
| pseudomonas_aeruginosa | clostridioides_difficile | 3e-11 |
| pseudomonas_aeruginosa | bacteroides_fragilis | 3e-11 |
| pseudomonas_aeruginosa | enterobacter_cloacae | 1e-10 |
| pseudomonas_aeruginosa | yersinia_enterocolitica | 1e-10 |
| pseudomonas_aeruginosa | moraxella_catarrhalis | 3e-11 |
| pseudomonas_aeruginosa | bordetella_pertussis | 3e-11 |
| pseudomonas_aeruginosa | mycoplasma_pneumoniae | 3e-11 |
| pseudomonas_aeruginosa | legionella_pneumophila | 3e-11 |
| pseudomonas_aeruginosa | burkholderia_cepacia_complex | 1e-9 |
| stenotrophomonas_maltophilia | acinetobacter_baumannii | 1e-9 |
| stenotrophomonas_maltophilia | citrobacter_spp. | 1e-10 |
| stenotrophomonas_maltophilia | enterobacter_spp. | 1e-10 |
| stenotrophomonas_maltophilia | escherichia_coli | 1e-10 |
| stenotrophomonas_maltophilia | klebsiella_pneumoniae | 1e-10 |
| stenotrophomonas_maltophilia | morganella_spp. | 1e-10 |
| stenotrophomonas_maltophilia | proteus_spp. | 1e-10 |
| stenotrophomonas_maltophilia | serratia_spp. | 1e-10 |
| stenotrophomonas_maltophilia | p_stuartii | 1e-10 |
| stenotrophomonas_maltophilia | pseudomonas_aeruginosa | 1e-9 |
| stenotrophomonas_maltophilia | salmonella_enterica_serovar_typhi | 1e-10 |
| stenotrophomonas_maltophilia | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| stenotrophomonas_maltophilia | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| stenotrophomonas_maltophilia | shigella_spp. | 1e-10 |
| stenotrophomonas_maltophilia | neisseria_gonorrhoeae | 3e-11 |
| stenotrophomonas_maltophilia | haemophilus_influenzae | 3e-11 |
| stenotrophomonas_maltophilia | chlamydia_trachomatis | 3e-11 |
| stenotrophomonas_maltophilia | mycoplasma_genitalium | 3e-11 |
| stenotrophomonas_maltophilia | vibrio_cholerae | 1e-10 |
| stenotrophomonas_maltophilia | neisseria_meningitidis | 3e-11 |
| stenotrophomonas_maltophilia | clostridioides_difficile | 3e-11 |
| stenotrophomonas_maltophilia | bacteroides_fragilis | 3e-11 |
| stenotrophomonas_maltophilia | enterobacter_cloacae | 1e-10 |
| stenotrophomonas_maltophilia | yersinia_enterocolitica | 1e-10 |
| stenotrophomonas_maltophilia | moraxella_catarrhalis | 3e-11 |
| stenotrophomonas_maltophilia | bordetella_pertussis | 3e-11 |
| stenotrophomonas_maltophilia | mycoplasma_pneumoniae | 3e-11 |
| stenotrophomonas_maltophilia | legionella_pneumophila | 3e-11 |
| stenotrophomonas_maltophilia | burkholderia_cepacia_complex | 1e-9 |
| staphylococcus_aureus | enterococcus_faecalis | 1e-9 |
| staphylococcus_aureus | enterococcus_faecium | 1e-9 |
| staphylococcus_aureus | staphylococcus_epidermidis | 1e-9 |
| staphylococcus_aureus | streptococcus_pneumoniae | 1e-9 |
| staphylococcus_aureus | streptococcus_pyogenes | 1e-9 |
| staphylococcus_aureus | streptococcus_agalactiae | 1e-9 |
| staphylococcus_aureus | listeria_monocytogenes | 1e-9 |
| staphylococcus_epidermidis | enterococcus_faecalis | 1e-9 |
| staphylococcus_epidermidis | enterococcus_faecium | 1e-9 |
| staphylococcus_epidermidis | staphylococcus_aureus | 1e-9 |
| staphylococcus_epidermidis | streptococcus_pneumoniae | 1e-9 |
| staphylococcus_epidermidis | streptococcus_pyogenes | 1e-9 |
| staphylococcus_epidermidis | streptococcus_agalactiae | 1e-9 |
| staphylococcus_epidermidis | listeria_monocytogenes | 1e-9 |
| streptococcus_pneumoniae | enterococcus_faecalis | 1e-9 |
| streptococcus_pneumoniae | enterococcus_faecium | 1e-9 |
| streptococcus_pneumoniae | staphylococcus_aureus | 1e-9 |
| streptococcus_pneumoniae | staphylococcus_epidermidis | 1e-9 |
| streptococcus_pneumoniae | streptococcus_pyogenes | 1e-9 |
| streptococcus_pneumoniae | streptococcus_agalactiae | 1e-9 |
| streptococcus_pneumoniae | listeria_monocytogenes | 1e-9 |
| salmonella_enterica_serovar_typhi | acinetobacter_baumannii | 1e-10 |
| salmonella_enterica_serovar_typhi | citrobacter_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | enterobacter_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | escherichia_coli | 1e-9 |
| salmonella_enterica_serovar_typhi | klebsiella_pneumoniae | 1e-9 |
| salmonella_enterica_serovar_typhi | morganella_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | proteus_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | serratia_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | p_stuartii | 1e-9 |
| salmonella_enterica_serovar_typhi | pseudomonas_aeruginosa | 1e-10 |
| salmonella_enterica_serovar_typhi | stenotrophomonas_maltophilia | 1e-10 |
| salmonella_enterica_serovar_typhi | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| salmonella_enterica_serovar_typhi | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | shigella_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | neisseria_gonorrhoeae | 3e-11 |
| salmonella_enterica_serovar_typhi | haemophilus_influenzae | 3e-11 |
| salmonella_enterica_serovar_typhi | chlamydia_trachomatis | 3e-11 |
| salmonella_enterica_serovar_typhi | mycoplasma_genitalium | 3e-11 |
| salmonella_enterica_serovar_typhi | vibrio_cholerae | 1e-10 |
| salmonella_enterica_serovar_typhi | neisseria_meningitidis | 3e-11 |
| salmonella_enterica_serovar_typhi | clostridioides_difficile | 3e-11 |
| salmonella_enterica_serovar_typhi | bacteroides_fragilis | 3e-11 |
| salmonella_enterica_serovar_typhi | enterobacter_cloacae | 1e-9 |
| salmonella_enterica_serovar_typhi | yersinia_enterocolitica | 1e-9 |
| salmonella_enterica_serovar_typhi | moraxella_catarrhalis | 3e-11 |
| salmonella_enterica_serovar_typhi | bordetella_pertussis | 3e-11 |
| salmonella_enterica_serovar_typhi | mycoplasma_pneumoniae | 3e-11 |
| salmonella_enterica_serovar_typhi | legionella_pneumophila | 3e-11 |
| salmonella_enterica_serovar_typhi | burkholderia_cepacia_complex | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | acinetobacter_baumannii | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | citrobacter_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | enterobacter_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | escherichia_coli | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | klebsiella_pneumoniae | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | morganella_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | proteus_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | serratia_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | p_stuartii | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | pseudomonas_aeruginosa | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | stenotrophomonas_maltophilia | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | salmonella_enterica_serovar_typhi | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | shigella_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | neisseria_gonorrhoeae | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | haemophilus_influenzae | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | chlamydia_trachomatis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | mycoplasma_genitalium | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | vibrio_cholerae | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | neisseria_meningitidis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | clostridioides_difficile | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | bacteroides_fragilis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | enterobacter_cloacae | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | yersinia_enterocolitica | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | moraxella_catarrhalis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | bordetella_pertussis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | mycoplasma_pneumoniae | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | legionella_pneumophila | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | burkholderia_cepacia_complex | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | acinetobacter_baumannii | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | citrobacter_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | enterobacter_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | escherichia_coli | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | klebsiella_pneumoniae | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | morganella_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | proteus_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | serratia_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | p_stuartii | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | pseudomonas_aeruginosa | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | stenotrophomonas_maltophilia | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | shigella_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | neisseria_gonorrhoeae | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | haemophilus_influenzae | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | chlamydia_trachomatis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | mycoplasma_genitalium | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | vibrio_cholerae | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | neisseria_meningitidis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | clostridioides_difficile | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | bacteroides_fragilis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | enterobacter_cloacae | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | yersinia_enterocolitica | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | moraxella_catarrhalis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | bordetella_pertussis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | mycoplasma_pneumoniae | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | legionella_pneumophila | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | burkholderia_cepacia_complex | 1e-10 |
| shigella_spp. | acinetobacter_baumannii | 1e-10 |
| shigella_spp. | citrobacter_spp. | 1e-9 |
| shigella_spp. | enterobacter_spp. | 1e-9 |
| shigella_spp. | escherichia_coli | 1e-9 |
| shigella_spp. | klebsiella_pneumoniae | 1e-9 |
| shigella_spp. | morganella_spp. | 1e-9 |
| shigella_spp. | proteus_spp. | 1e-9 |
| shigella_spp. | serratia_spp. | 1e-9 |
| shigella_spp. | p_stuartii | 1e-9 |
| shigella_spp. | pseudomonas_aeruginosa | 1e-10 |
| shigella_spp. | stenotrophomonas_maltophilia | 1e-10 |
| shigella_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| shigella_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| shigella_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| shigella_spp. | neisseria_gonorrhoeae | 3e-11 |
| shigella_spp. | haemophilus_influenzae | 3e-11 |
| shigella_spp. | chlamydia_trachomatis | 3e-11 |
| shigella_spp. | mycoplasma_genitalium | 3e-11 |
| shigella_spp. | vibrio_cholerae | 1e-10 |
| shigella_spp. | neisseria_meningitidis | 3e-11 |
| shigella_spp. | clostridioides_difficile | 3e-11 |
| shigella_spp. | bacteroides_fragilis | 3e-11 |
| shigella_spp. | enterobacter_cloacae | 1e-9 |
| shigella_spp. | yersinia_enterocolitica | 1e-9 |
| shigella_spp. | moraxella_catarrhalis | 3e-11 |
| shigella_spp. | bordetella_pertussis | 3e-11 |
| shigella_spp. | mycoplasma_pneumoniae | 3e-11 |
| shigella_spp. | legionella_pneumophila | 3e-11 |
| shigella_spp. | burkholderia_cepacia_complex | 1e-10 |
| neisseria_gonorrhoeae | acinetobacter_baumannii | 3e-11 |
| neisseria_gonorrhoeae | citrobacter_spp. | 3e-11 |
| neisseria_gonorrhoeae | enterobacter_spp. | 3e-11 |
| neisseria_gonorrhoeae | escherichia_coli | 3e-11 |
| neisseria_gonorrhoeae | klebsiella_pneumoniae | 3e-11 |
| neisseria_gonorrhoeae | morganella_spp. | 3e-11 |
| neisseria_gonorrhoeae | proteus_spp. | 3e-11 |
| neisseria_gonorrhoeae | serratia_spp. | 3e-11 |
| neisseria_gonorrhoeae | p_stuartii | 3e-11 |
| neisseria_gonorrhoeae | pseudomonas_aeruginosa | 3e-11 |
| neisseria_gonorrhoeae | stenotrophomonas_maltophilia | 3e-11 |
| neisseria_gonorrhoeae | salmonella_enterica_serovar_typhi | 3e-11 |
| neisseria_gonorrhoeae | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| neisseria_gonorrhoeae | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| neisseria_gonorrhoeae | shigella_spp. | 3e-11 |
| neisseria_gonorrhoeae | haemophilus_influenzae | 1e-9 |
| neisseria_gonorrhoeae | chlamydia_trachomatis | 1e-9 |
| neisseria_gonorrhoeae | mycoplasma_genitalium | 1e-9 |
| neisseria_gonorrhoeae | vibrio_cholerae | 3e-11 |
| neisseria_gonorrhoeae | neisseria_meningitidis | 1e-9 |
| neisseria_gonorrhoeae | enterobacter_cloacae | 3e-11 |
| neisseria_gonorrhoeae | yersinia_enterocolitica | 3e-11 |
| neisseria_gonorrhoeae | moraxella_catarrhalis | 1e-9 |
| neisseria_gonorrhoeae | bordetella_pertussis | 1e-9 |
| neisseria_gonorrhoeae | mycoplasma_pneumoniae | 1e-9 |
| neisseria_gonorrhoeae | legionella_pneumophila | 1e-9 |
| neisseria_gonorrhoeae | burkholderia_cepacia_complex | 3e-11 |
| streptococcus_pyogenes | enterococcus_faecalis | 1e-9 |
| streptococcus_pyogenes | enterococcus_faecium | 1e-9 |
| streptococcus_pyogenes | staphylococcus_aureus | 1e-9 |
| streptococcus_pyogenes | staphylococcus_epidermidis | 1e-9 |
| streptococcus_pyogenes | streptococcus_pneumoniae | 1e-9 |
| streptococcus_pyogenes | streptococcus_agalactiae | 1e-9 |
| streptococcus_pyogenes | listeria_monocytogenes | 1e-9 |
| streptococcus_agalactiae | enterococcus_faecalis | 1e-9 |
| streptococcus_agalactiae | enterococcus_faecium | 1e-9 |
| streptococcus_agalactiae | staphylococcus_aureus | 1e-9 |
| streptococcus_agalactiae | staphylococcus_epidermidis | 1e-9 |
| streptococcus_agalactiae | streptococcus_pneumoniae | 1e-9 |
| streptococcus_agalactiae | streptococcus_pyogenes | 1e-9 |
| streptococcus_agalactiae | listeria_monocytogenes | 1e-9 |
| haemophilus_influenzae | acinetobacter_baumannii | 3e-11 |
| haemophilus_influenzae | citrobacter_spp. | 3e-11 |
| haemophilus_influenzae | enterobacter_spp. | 3e-11 |
| haemophilus_influenzae | escherichia_coli | 3e-11 |
| haemophilus_influenzae | klebsiella_pneumoniae | 3e-11 |
| haemophilus_influenzae | morganella_spp. | 3e-11 |
| haemophilus_influenzae | proteus_spp. | 3e-11 |
| haemophilus_influenzae | serratia_spp. | 3e-11 |
| haemophilus_influenzae | p_stuartii | 3e-11 |
| haemophilus_influenzae | pseudomonas_aeruginosa | 3e-11 |
| haemophilus_influenzae | stenotrophomonas_maltophilia | 3e-11 |
| haemophilus_influenzae | salmonella_enterica_serovar_typhi | 3e-11 |
| haemophilus_influenzae | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| haemophilus_influenzae | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| haemophilus_influenzae | shigella_spp. | 3e-11 |
| haemophilus_influenzae | neisseria_gonorrhoeae | 1e-9 |
| haemophilus_influenzae | chlamydia_trachomatis | 1e-9 |
| haemophilus_influenzae | mycoplasma_genitalium | 1e-9 |
| haemophilus_influenzae | vibrio_cholerae | 3e-11 |
| haemophilus_influenzae | neisseria_meningitidis | 1e-9 |
| haemophilus_influenzae | enterobacter_cloacae | 3e-11 |
| haemophilus_influenzae | yersinia_enterocolitica | 3e-11 |
| haemophilus_influenzae | moraxella_catarrhalis | 1e-9 |
| haemophilus_influenzae | bordetella_pertussis | 1e-9 |
| haemophilus_influenzae | mycoplasma_pneumoniae | 1e-9 |
| haemophilus_influenzae | legionella_pneumophila | 1e-9 |
| haemophilus_influenzae | burkholderia_cepacia_complex | 3e-11 |
| chlamydia_trachomatis | acinetobacter_baumannii | 3e-11 |
| chlamydia_trachomatis | citrobacter_spp. | 3e-11 |
| chlamydia_trachomatis | enterobacter_spp. | 3e-11 |
| chlamydia_trachomatis | escherichia_coli | 3e-11 |
| chlamydia_trachomatis | klebsiella_pneumoniae | 3e-11 |
| chlamydia_trachomatis | morganella_spp. | 3e-11 |
| chlamydia_trachomatis | proteus_spp. | 3e-11 |
| chlamydia_trachomatis | serratia_spp. | 3e-11 |
| chlamydia_trachomatis | p_stuartii | 3e-11 |
| chlamydia_trachomatis | pseudomonas_aeruginosa | 3e-11 |
| chlamydia_trachomatis | stenotrophomonas_maltophilia | 3e-11 |
| chlamydia_trachomatis | salmonella_enterica_serovar_typhi | 3e-11 |
| chlamydia_trachomatis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| chlamydia_trachomatis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| chlamydia_trachomatis | shigella_spp. | 3e-11 |
| chlamydia_trachomatis | neisseria_gonorrhoeae | 1e-9 |
| chlamydia_trachomatis | haemophilus_influenzae | 1e-9 |
| chlamydia_trachomatis | mycoplasma_genitalium | 1e-9 |
| chlamydia_trachomatis | vibrio_cholerae | 3e-11 |
| chlamydia_trachomatis | neisseria_meningitidis | 1e-9 |
| chlamydia_trachomatis | enterobacter_cloacae | 3e-11 |
| chlamydia_trachomatis | yersinia_enterocolitica | 3e-11 |
| chlamydia_trachomatis | moraxella_catarrhalis | 1e-9 |
| chlamydia_trachomatis | bordetella_pertussis | 1e-9 |
| chlamydia_trachomatis | mycoplasma_pneumoniae | 1e-9 |
| chlamydia_trachomatis | legionella_pneumophila | 1e-9 |
| chlamydia_trachomatis | burkholderia_cepacia_complex | 3e-11 |
| mycoplasma_genitalium | acinetobacter_baumannii | 3e-11 |
| mycoplasma_genitalium | citrobacter_spp. | 3e-11 |
| mycoplasma_genitalium | enterobacter_spp. | 3e-11 |
| mycoplasma_genitalium | escherichia_coli | 3e-11 |
| mycoplasma_genitalium | klebsiella_pneumoniae | 3e-11 |
| mycoplasma_genitalium | morganella_spp. | 3e-11 |
| mycoplasma_genitalium | proteus_spp. | 3e-11 |
| mycoplasma_genitalium | serratia_spp. | 3e-11 |
| mycoplasma_genitalium | p_stuartii | 3e-11 |
| mycoplasma_genitalium | pseudomonas_aeruginosa | 3e-11 |
| mycoplasma_genitalium | stenotrophomonas_maltophilia | 3e-11 |
| mycoplasma_genitalium | salmonella_enterica_serovar_typhi | 3e-11 |
| mycoplasma_genitalium | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| mycoplasma_genitalium | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| mycoplasma_genitalium | shigella_spp. | 3e-11 |
| mycoplasma_genitalium | neisseria_gonorrhoeae | 1e-9 |
| mycoplasma_genitalium | haemophilus_influenzae | 1e-9 |
| mycoplasma_genitalium | chlamydia_trachomatis | 1e-9 |
| mycoplasma_genitalium | vibrio_cholerae | 3e-11 |
| mycoplasma_genitalium | neisseria_meningitidis | 1e-9 |
| mycoplasma_genitalium | enterobacter_cloacae | 3e-11 |
| mycoplasma_genitalium | yersinia_enterocolitica | 3e-11 |
| mycoplasma_genitalium | moraxella_catarrhalis | 1e-9 |
| mycoplasma_genitalium | bordetella_pertussis | 1e-9 |
| mycoplasma_genitalium | mycoplasma_pneumoniae | 1e-9 |
| mycoplasma_genitalium | legionella_pneumophila | 1e-9 |
| mycoplasma_genitalium | burkholderia_cepacia_complex | 3e-11 |
| vibrio_cholerae | acinetobacter_baumannii | 1e-10 |
| vibrio_cholerae | citrobacter_spp. | 1e-10 |
| vibrio_cholerae | enterobacter_spp. | 1e-10 |
| vibrio_cholerae | escherichia_coli | 1e-10 |
| vibrio_cholerae | klebsiella_pneumoniae | 1e-10 |
| vibrio_cholerae | morganella_spp. | 1e-10 |
| vibrio_cholerae | proteus_spp. | 1e-10 |
| vibrio_cholerae | serratia_spp. | 1e-10 |
| vibrio_cholerae | p_stuartii | 1e-10 |
| vibrio_cholerae | pseudomonas_aeruginosa | 1e-10 |
| vibrio_cholerae | stenotrophomonas_maltophilia | 1e-10 |
| vibrio_cholerae | salmonella_enterica_serovar_typhi | 1e-10 |
| vibrio_cholerae | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| vibrio_cholerae | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| vibrio_cholerae | shigella_spp. | 1e-10 |
| vibrio_cholerae | neisseria_gonorrhoeae | 3e-11 |
| vibrio_cholerae | haemophilus_influenzae | 3e-11 |
| vibrio_cholerae | chlamydia_trachomatis | 3e-11 |
| vibrio_cholerae | mycoplasma_genitalium | 3e-11 |
| vibrio_cholerae | neisseria_meningitidis | 3e-11 |
| vibrio_cholerae | clostridioides_difficile | 3e-11 |
| vibrio_cholerae | bacteroides_fragilis | 3e-11 |
| vibrio_cholerae | enterobacter_cloacae | 1e-10 |
| vibrio_cholerae | yersinia_enterocolitica | 1e-10 |
| vibrio_cholerae | moraxella_catarrhalis | 3e-11 |
| vibrio_cholerae | bordetella_pertussis | 3e-11 |
| vibrio_cholerae | mycoplasma_pneumoniae | 3e-11 |
| vibrio_cholerae | legionella_pneumophila | 3e-11 |
| vibrio_cholerae | burkholderia_cepacia_complex | 1e-10 |
| neisseria_meningitidis | acinetobacter_baumannii | 3e-11 |
| neisseria_meningitidis | citrobacter_spp. | 3e-11 |
| neisseria_meningitidis | enterobacter_spp. | 3e-11 |
| neisseria_meningitidis | escherichia_coli | 3e-11 |
| neisseria_meningitidis | klebsiella_pneumoniae | 3e-11 |
| neisseria_meningitidis | morganella_spp. | 3e-11 |
| neisseria_meningitidis | proteus_spp. | 3e-11 |
| neisseria_meningitidis | serratia_spp. | 3e-11 |
| neisseria_meningitidis | p_stuartii | 3e-11 |
| neisseria_meningitidis | pseudomonas_aeruginosa | 3e-11 |
| neisseria_meningitidis | stenotrophomonas_maltophilia | 3e-11 |
| neisseria_meningitidis | salmonella_enterica_serovar_typhi | 3e-11 |
| neisseria_meningitidis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| neisseria_meningitidis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| neisseria_meningitidis | shigella_spp. | 3e-11 |
| neisseria_meningitidis | neisseria_gonorrhoeae | 1e-9 |
| neisseria_meningitidis | haemophilus_influenzae | 1e-9 |
| neisseria_meningitidis | chlamydia_trachomatis | 1e-9 |
| neisseria_meningitidis | mycoplasma_genitalium | 1e-9 |
| neisseria_meningitidis | vibrio_cholerae | 3e-11 |
| neisseria_meningitidis | enterobacter_cloacae | 3e-11 |
| neisseria_meningitidis | yersinia_enterocolitica | 3e-11 |
| neisseria_meningitidis | moraxella_catarrhalis | 1e-9 |
| neisseria_meningitidis | bordetella_pertussis | 1e-9 |
| neisseria_meningitidis | mycoplasma_pneumoniae | 1e-9 |
| neisseria_meningitidis | legionella_pneumophila | 1e-9 |
| neisseria_meningitidis | burkholderia_cepacia_complex | 3e-11 |
| listeria_monocytogenes | enterococcus_faecalis | 1e-9 |
| listeria_monocytogenes | enterococcus_faecium | 1e-9 |
| listeria_monocytogenes | staphylococcus_aureus | 1e-9 |
| listeria_monocytogenes | staphylococcus_epidermidis | 1e-9 |
| listeria_monocytogenes | streptococcus_pneumoniae | 1e-9 |
| listeria_monocytogenes | streptococcus_pyogenes | 1e-9 |
| listeria_monocytogenes | streptococcus_agalactiae | 1e-9 |
| clostridioides_difficile | acinetobacter_baumannii | 3e-11 |
| clostridioides_difficile | citrobacter_spp. | 3e-11 |
| clostridioides_difficile | enterobacter_spp. | 3e-11 |
| clostridioides_difficile | escherichia_coli | 3e-11 |
| clostridioides_difficile | klebsiella_pneumoniae | 3e-11 |
| clostridioides_difficile | morganella_spp. | 3e-11 |
| clostridioides_difficile | proteus_spp. | 3e-11 |
| clostridioides_difficile | serratia_spp. | 3e-11 |
| clostridioides_difficile | p_stuartii | 3e-11 |
| clostridioides_difficile | pseudomonas_aeruginosa | 3e-11 |
| clostridioides_difficile | stenotrophomonas_maltophilia | 3e-11 |
| clostridioides_difficile | salmonella_enterica_serovar_typhi | 3e-11 |
| clostridioides_difficile | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| clostridioides_difficile | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| clostridioides_difficile | shigella_spp. | 3e-11 |
| clostridioides_difficile | vibrio_cholerae | 3e-11 |
| clostridioides_difficile | bacteroides_fragilis | 1e-9 |
| clostridioides_difficile | enterobacter_cloacae | 3e-11 |
| clostridioides_difficile | yersinia_enterocolitica | 3e-11 |
| clostridioides_difficile | burkholderia_cepacia_complex | 3e-11 |
| bacteroides_fragilis | acinetobacter_baumannii | 3e-11 |
| bacteroides_fragilis | citrobacter_spp. | 3e-11 |
| bacteroides_fragilis | enterobacter_spp. | 3e-11 |
| bacteroides_fragilis | escherichia_coli | 3e-11 |
| bacteroides_fragilis | klebsiella_pneumoniae | 3e-11 |
| bacteroides_fragilis | morganella_spp. | 3e-11 |
| bacteroides_fragilis | proteus_spp. | 3e-11 |
| bacteroides_fragilis | serratia_spp. | 3e-11 |
| bacteroides_fragilis | p_stuartii | 3e-11 |
| bacteroides_fragilis | pseudomonas_aeruginosa | 3e-11 |
| bacteroides_fragilis | stenotrophomonas_maltophilia | 3e-11 |
| bacteroides_fragilis | salmonella_enterica_serovar_typhi | 3e-11 |
| bacteroides_fragilis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| bacteroides_fragilis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| bacteroides_fragilis | shigella_spp. | 3e-11 |
| bacteroides_fragilis | vibrio_cholerae | 3e-11 |
| bacteroides_fragilis | clostridioides_difficile | 1e-9 |
| bacteroides_fragilis | enterobacter_cloacae | 3e-11 |
| bacteroides_fragilis | yersinia_enterocolitica | 3e-11 |
| bacteroides_fragilis | burkholderia_cepacia_complex | 3e-11 |
| enterobacter_cloacae | acinetobacter_baumannii | 1e-10 |
| enterobacter_cloacae | citrobacter_spp. | 1e-9 |
| enterobacter_cloacae | enterobacter_spp. | 1e-9 |
| enterobacter_cloacae | escherichia_coli | 1e-9 |
| enterobacter_cloacae | klebsiella_pneumoniae | 1e-9 |
| enterobacter_cloacae | morganella_spp. | 1e-9 |
| enterobacter_cloacae | proteus_spp. | 1e-9 |
| enterobacter_cloacae | serratia_spp. | 1e-9 |
| enterobacter_cloacae | p_stuartii | 1e-9 |
| enterobacter_cloacae | pseudomonas_aeruginosa | 1e-10 |
| enterobacter_cloacae | stenotrophomonas_maltophilia | 1e-10 |
| enterobacter_cloacae | salmonella_enterica_serovar_typhi | 1e-9 |
| enterobacter_cloacae | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| enterobacter_cloacae | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| enterobacter_cloacae | shigella_spp. | 1e-9 |
| enterobacter_cloacae | neisseria_gonorrhoeae | 3e-11 |
| enterobacter_cloacae | haemophilus_influenzae | 3e-11 |
| enterobacter_cloacae | chlamydia_trachomatis | 3e-11 |
| enterobacter_cloacae | mycoplasma_genitalium | 3e-11 |
| enterobacter_cloacae | vibrio_cholerae | 1e-10 |
| enterobacter_cloacae | neisseria_meningitidis | 3e-11 |
| enterobacter_cloacae | clostridioides_difficile | 3e-11 |
| enterobacter_cloacae | bacteroides_fragilis | 3e-11 |
| enterobacter_cloacae | yersinia_enterocolitica | 1e-9 |
| enterobacter_cloacae | moraxella_catarrhalis | 3e-11 |
| enterobacter_cloacae | bordetella_pertussis | 3e-11 |
| enterobacter_cloacae | mycoplasma_pneumoniae | 3e-11 |
| enterobacter_cloacae | legionella_pneumophila | 3e-11 |
| enterobacter_cloacae | burkholderia_cepacia_complex | 1e-10 |
| yersinia_enterocolitica | acinetobacter_baumannii | 1e-10 |
| yersinia_enterocolitica | citrobacter_spp. | 1e-9 |
| yersinia_enterocolitica | enterobacter_spp. | 1e-9 |
| yersinia_enterocolitica | escherichia_coli | 1e-9 |
| yersinia_enterocolitica | klebsiella_pneumoniae | 1e-9 |
| yersinia_enterocolitica | morganella_spp. | 1e-9 |
| yersinia_enterocolitica | proteus_spp. | 1e-9 |
| yersinia_enterocolitica | serratia_spp. | 1e-9 |
| yersinia_enterocolitica | p_stuartii | 1e-9 |
| yersinia_enterocolitica | pseudomonas_aeruginosa | 1e-10 |
| yersinia_enterocolitica | stenotrophomonas_maltophilia | 1e-10 |
| yersinia_enterocolitica | salmonella_enterica_serovar_typhi | 1e-9 |
| yersinia_enterocolitica | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| yersinia_enterocolitica | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| yersinia_enterocolitica | shigella_spp. | 1e-9 |
| yersinia_enterocolitica | neisseria_gonorrhoeae | 3e-11 |
| yersinia_enterocolitica | haemophilus_influenzae | 3e-11 |
| yersinia_enterocolitica | chlamydia_trachomatis | 3e-11 |
| yersinia_enterocolitica | mycoplasma_genitalium | 3e-11 |
| yersinia_enterocolitica | vibrio_cholerae | 1e-10 |
| yersinia_enterocolitica | neisseria_meningitidis | 3e-11 |
| yersinia_enterocolitica | clostridioides_difficile | 3e-11 |
| yersinia_enterocolitica | bacteroides_fragilis | 3e-11 |
| yersinia_enterocolitica | enterobacter_cloacae | 1e-9 |
| yersinia_enterocolitica | moraxella_catarrhalis | 3e-11 |
| yersinia_enterocolitica | bordetella_pertussis | 3e-11 |
| yersinia_enterocolitica | mycoplasma_pneumoniae | 3e-11 |
| yersinia_enterocolitica | legionella_pneumophila | 3e-11 |
| yersinia_enterocolitica | burkholderia_cepacia_complex | 1e-10 |
| moraxella_catarrhalis | acinetobacter_baumannii | 3e-11 |
| moraxella_catarrhalis | citrobacter_spp. | 3e-11 |
| moraxella_catarrhalis | enterobacter_spp. | 3e-11 |
| moraxella_catarrhalis | escherichia_coli | 3e-11 |
| moraxella_catarrhalis | klebsiella_pneumoniae | 3e-11 |
| moraxella_catarrhalis | morganella_spp. | 3e-11 |
| moraxella_catarrhalis | proteus_spp. | 3e-11 |
| moraxella_catarrhalis | serratia_spp. | 3e-11 |
| moraxella_catarrhalis | p_stuartii | 3e-11 |
| moraxella_catarrhalis | pseudomonas_aeruginosa | 3e-11 |
| moraxella_catarrhalis | stenotrophomonas_maltophilia | 3e-11 |
| moraxella_catarrhalis | salmonella_enterica_serovar_typhi | 3e-11 |
| moraxella_catarrhalis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| moraxella_catarrhalis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| moraxella_catarrhalis | shigella_spp. | 3e-11 |
| moraxella_catarrhalis | neisseria_gonorrhoeae | 1e-9 |
| moraxella_catarrhalis | haemophilus_influenzae | 1e-9 |
| moraxella_catarrhalis | chlamydia_trachomatis | 1e-9 |
| moraxella_catarrhalis | mycoplasma_genitalium | 1e-9 |
| moraxella_catarrhalis | vibrio_cholerae | 3e-11 |
| moraxella_catarrhalis | neisseria_meningitidis | 1e-9 |
| moraxella_catarrhalis | enterobacter_cloacae | 3e-11 |
| moraxella_catarrhalis | yersinia_enterocolitica | 3e-11 |
| moraxella_catarrhalis | bordetella_pertussis | 1e-9 |
| moraxella_catarrhalis | mycoplasma_pneumoniae | 1e-9 |
| moraxella_catarrhalis | legionella_pneumophila | 1e-9 |
| moraxella_catarrhalis | burkholderia_cepacia_complex | 3e-11 |
| bordetella_pertussis | acinetobacter_baumannii | 3e-11 |
| bordetella_pertussis | citrobacter_spp. | 3e-11 |
| bordetella_pertussis | enterobacter_spp. | 3e-11 |
| bordetella_pertussis | escherichia_coli | 3e-11 |
| bordetella_pertussis | klebsiella_pneumoniae | 3e-11 |
| bordetella_pertussis | morganella_spp. | 3e-11 |
| bordetella_pertussis | proteus_spp. | 3e-11 |
| bordetella_pertussis | serratia_spp. | 3e-11 |
| bordetella_pertussis | p_stuartii | 3e-11 |
| bordetella_pertussis | pseudomonas_aeruginosa | 3e-11 |
| bordetella_pertussis | stenotrophomonas_maltophilia | 3e-11 |
| bordetella_pertussis | salmonella_enterica_serovar_typhi | 3e-11 |
| bordetella_pertussis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| bordetella_pertussis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| bordetella_pertussis | shigella_spp. | 3e-11 |
| bordetella_pertussis | neisseria_gonorrhoeae | 1e-9 |
| bordetella_pertussis | haemophilus_influenzae | 1e-9 |
| bordetella_pertussis | chlamydia_trachomatis | 1e-9 |
| bordetella_pertussis | mycoplasma_genitalium | 1e-9 |
| bordetella_pertussis | vibrio_cholerae | 3e-11 |
| bordetella_pertussis | neisseria_meningitidis | 1e-9 |
| bordetella_pertussis | enterobacter_cloacae | 3e-11 |
| bordetella_pertussis | yersinia_enterocolitica | 3e-11 |
| bordetella_pertussis | moraxella_catarrhalis | 1e-9 |
| bordetella_pertussis | mycoplasma_pneumoniae | 1e-9 |
| bordetella_pertussis | legionella_pneumophila | 1e-9 |
| bordetella_pertussis | burkholderia_cepacia_complex | 3e-11 |
| mycoplasma_pneumoniae | acinetobacter_baumannii | 3e-11 |
| mycoplasma_pneumoniae | citrobacter_spp. | 3e-11 |
| mycoplasma_pneumoniae | enterobacter_spp. | 3e-11 |
| mycoplasma_pneumoniae | escherichia_coli | 3e-11 |
| mycoplasma_pneumoniae | klebsiella_pneumoniae | 3e-11 |
| mycoplasma_pneumoniae | morganella_spp. | 3e-11 |
| mycoplasma_pneumoniae | proteus_spp. | 3e-11 |
| mycoplasma_pneumoniae | serratia_spp. | 3e-11 |
| mycoplasma_pneumoniae | p_stuartii | 3e-11 |
| mycoplasma_pneumoniae | pseudomonas_aeruginosa | 3e-11 |
| mycoplasma_pneumoniae | stenotrophomonas_maltophilia | 3e-11 |
| mycoplasma_pneumoniae | salmonella_enterica_serovar_typhi | 3e-11 |
| mycoplasma_pneumoniae | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| mycoplasma_pneumoniae | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| mycoplasma_pneumoniae | shigella_spp. | 3e-11 |
| mycoplasma_pneumoniae | neisseria_gonorrhoeae | 1e-9 |
| mycoplasma_pneumoniae | haemophilus_influenzae | 1e-9 |
| mycoplasma_pneumoniae | chlamydia_trachomatis | 1e-9 |
| mycoplasma_pneumoniae | mycoplasma_genitalium | 1e-9 |
| mycoplasma_pneumoniae | vibrio_cholerae | 3e-11 |
| mycoplasma_pneumoniae | neisseria_meningitidis | 1e-9 |
| mycoplasma_pneumoniae | enterobacter_cloacae | 3e-11 |
| mycoplasma_pneumoniae | yersinia_enterocolitica | 3e-11 |
| mycoplasma_pneumoniae | moraxella_catarrhalis | 1e-9 |
| mycoplasma_pneumoniae | bordetella_pertussis | 1e-9 |
| mycoplasma_pneumoniae | legionella_pneumophila | 1e-9 |
| mycoplasma_pneumoniae | burkholderia_cepacia_complex | 3e-11 |
| legionella_pneumophila | acinetobacter_baumannii | 3e-11 |
| legionella_pneumophila | citrobacter_spp. | 3e-11 |
| legionella_pneumophila | enterobacter_spp. | 3e-11 |
| legionella_pneumophila | escherichia_coli | 3e-11 |
| legionella_pneumophila | klebsiella_pneumoniae | 3e-11 |
| legionella_pneumophila | morganella_spp. | 3e-11 |
| legionella_pneumophila | proteus_spp. | 3e-11 |
| legionella_pneumophila | serratia_spp. | 3e-11 |
| legionella_pneumophila | p_stuartii | 3e-11 |
| legionella_pneumophila | pseudomonas_aeruginosa | 3e-11 |
| legionella_pneumophila | stenotrophomonas_maltophilia | 3e-11 |
| legionella_pneumophila | salmonella_enterica_serovar_typhi | 3e-11 |
| legionella_pneumophila | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| legionella_pneumophila | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| legionella_pneumophila | shigella_spp. | 3e-11 |
| legionella_pneumophila | neisseria_gonorrhoeae | 1e-9 |
| legionella_pneumophila | haemophilus_influenzae | 1e-9 |
| legionella_pneumophila | chlamydia_trachomatis | 1e-9 |
| legionella_pneumophila | mycoplasma_genitalium | 1e-9 |
| legionella_pneumophila | vibrio_cholerae | 3e-11 |
| legionella_pneumophila | neisseria_meningitidis | 1e-9 |
| legionella_pneumophila | enterobacter_cloacae | 3e-11 |
| legionella_pneumophila | yersinia_enterocolitica | 3e-11 |
| legionella_pneumophila | moraxella_catarrhalis | 1e-9 |
| legionella_pneumophila | bordetella_pertussis | 1e-9 |
| legionella_pneumophila | mycoplasma_pneumoniae | 1e-9 |
| legionella_pneumophila | burkholderia_cepacia_complex | 3e-11 |
| burkholderia_cepacia_complex | acinetobacter_baumannii | 1e-9 |
| burkholderia_cepacia_complex | citrobacter_spp. | 1e-10 |
| burkholderia_cepacia_complex | enterobacter_spp. | 1e-10 |
| burkholderia_cepacia_complex | escherichia_coli | 1e-10 |
| burkholderia_cepacia_complex | klebsiella_pneumoniae | 1e-10 |
| burkholderia_cepacia_complex | morganella_spp. | 1e-10 |
| burkholderia_cepacia_complex | proteus_spp. | 1e-10 |
| burkholderia_cepacia_complex | serratia_spp. | 1e-10 |
| burkholderia_cepacia_complex | p_stuartii | 1e-10 |
| burkholderia_cepacia_complex | pseudomonas_aeruginosa | 1e-9 |
| burkholderia_cepacia_complex | stenotrophomonas_maltophilia | 1e-9 |
| burkholderia_cepacia_complex | salmonella_enterica_serovar_typhi | 1e-10 |
| burkholderia_cepacia_complex | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| burkholderia_cepacia_complex | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| burkholderia_cepacia_complex | shigella_spp. | 1e-10 |
| burkholderia_cepacia_complex | neisseria_gonorrhoeae | 3e-11 |
| burkholderia_cepacia_complex | haemophilus_influenzae | 3e-11 |
| burkholderia_cepacia_complex | chlamydia_trachomatis | 3e-11 |
| burkholderia_cepacia_complex | mycoplasma_genitalium | 3e-11 |
| burkholderia_cepacia_complex | vibrio_cholerae | 1e-10 |
| burkholderia_cepacia_complex | neisseria_meningitidis | 3e-11 |
| burkholderia_cepacia_complex | clostridioides_difficile | 3e-11 |
| burkholderia_cepacia_complex | bacteroides_fragilis | 3e-11 |
| burkholderia_cepacia_complex | enterobacter_cloacae | 1e-10 |
| burkholderia_cepacia_complex | yersinia_enterocolitica | 1e-10 |
| burkholderia_cepacia_complex | moraxella_catarrhalis | 3e-11 |
| burkholderia_cepacia_complex | bordetella_pertussis | 3e-11 |
| burkholderia_cepacia_complex | mycoplasma_pneumoniae | 3e-11 |
| burkholderia_cepacia_complex | legionella_pneumophila | 3e-11 |
