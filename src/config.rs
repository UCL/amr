// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES}; // Import both lists

// --- Global Simulation Parameters ---
lazy_static! {
    pub static ref PARAMETERS: HashMap<String, f64> = {
        let mut map = HashMap::new();


        // --- Default Parameters for ALL Bacteria from BACTERIA_LIST ---
        // These are inserted first, and can then be overridden by specific entries below.
        for &bacteria in BACTERIA_LIST.iter() {
            map.insert(format!("{}_acquisition_prob_baseline", bacteria), 0.001); // 0.0001
            map.insert(format!("{}_initial_infection_level", bacteria), 0.01); // 0.01
            map.insert(format!("{}_environmental_acquisition_proportion", bacteria), 0.8); // 0.1
            map.insert(format!("{}_hospital_acquired_multiplier", bacteria), 10.0); // multiplier for hospital-acquired risk
            map.insert(format!("{}_adult_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_child_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_oral_exposure_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_sexual_contact_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_mosquito_exposure_acq_rate_ratio_per_unit", bacteria), 1.0);
            map.insert(format!("{}_vaccine_efficacy", bacteria), 0.0); // Default to no vaccine effect
            map.insert(format!("{}_base_bacteria_level_change", bacteria), 0.5); // 0.2 
            map.insert(format!("{}_max_level", bacteria), 5.0);
            map.insert(format!("{}_immunity_effect_on_level_change", bacteria), 0.005); // 0.05 is strong effect
            map.insert(format!("{}_immunity_base_response", bacteria), 0.1); // 0.001
            map.insert(format!("{}_immunity_increase_per_unit_higher_bacteria_level", bacteria), 0.05);
            map.insert(format!("{}_immunity_increase_per_infection_day", bacteria), 0.05);
            map.insert(format!("{}_immunity_age_modifier", bacteria), 1.0);
            map.insert(format!("{}_immunity_immunodeficiency_modifier", bacteria), 0.1);
            map.insert(format!("{}_max_immune_response", bacteria), 10.0); // Maximum immune response level
            
            // Age-related infection risk parameters
            map.insert(format!("{}_age_effect_scaling", bacteria), 1.0); // Scale the template effect (1.0 = full effect)
        }



        // General Drug Parameters
        map.insert("drug_base_initiation_rate_per_day".to_string(), 0.1); // 0.0001
        map.insert("drug_infection_present_multiplier".to_string(), 50.0);
        map.insert("drug_test_identified_multiplier".to_string(), 50.0);
        map.insert("drug_decay_per_day".to_string(), 1.0); // Legacy parameter - now using drug-specific half-lives
        
        // Drug-specific half-lives (in days) for realistic pharmacokinetics
        
        // Beta-lactams (Penicillins)
        map.insert("drug_penicilling_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ampicillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_amoxicillin_half_life_days".to_string(), 0.04); // ~1 hour  
        map.insert("drug_piperacillin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ticarcillin_half_life_days".to_string(), 0.046); // ~1.1 hours
        
        // Cephalosporins
        map.insert("drug_cephalexin_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_cefazolin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_cefuroxime_half_life_days".to_string(), 0.05); // ~1.3 hours
        map.insert("drug_ceftriaxone_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_ceftazidime_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_cefepime_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_ceftaroline_half_life_days".to_string(), 0.11); // ~2.6 hours
        
        // Carbapenems
        map.insert("drug_meropenem_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_imipenem_c_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ertapenem_half_life_days".to_string(), 0.17); // ~4 hours
        
        // Monobactams
        map.insert("drug_aztreonam_half_life_days".to_string(), 0.08); // ~2 hours
        
        // Macrolides
        map.insert("drug_erythromycin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_azithromycin_half_life_days".to_string(), 2.8); // ~68 hours
        map.insert("drug_clarithromycin_half_life_days".to_string(), 0.25); // ~6 hours
        
        // Lincosamides
        map.insert("drug_clindamycin_half_life_days".to_string(), 0.125); // ~3 hours
        
        // Aminoglycosides
        map.insert("drug_gentamicin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_tobramycin_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_amikacin_half_life_days".to_string(), 0.08); // ~2 hours
        
        // Fluoroquinolones
        map.insert("drug_ciprofloxacin_half_life_days".to_string(), 0.17); // ~4 hours
        map.insert("drug_levofloxacin_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_moxifloxacin_half_life_days".to_string(), 0.5); // ~12 hours
        map.insert("drug_ofloxacin_half_life_days".to_string(), 0.25); // ~6 hours
        
        // Tetracyclines
        map.insert("drug_tetracycline_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_doxyclycline_half_life_days".to_string(), 0.75); // ~18 hours
        map.insert("drug_minocycline_half_life_days".to_string(), 0.67); // ~16 hours
        
        // Glycopeptides
        map.insert("drug_vancomycin_half_life_days".to_string(), 0.25); // ~6 hours
        map.insert("drug_teicoplanin_half_life_days".to_string(), 3.5); // ~83 hours (very long)
        
        // Oxazolidinones
        map.insert("drug_linezolid_half_life_days".to_string(), 0.21); // ~5 hours
        map.insert("drug_tedizolid_half_life_days".to_string(), 0.5); // ~12 hours
        
        // Quinolones (older)
        map.insert("drug_quinu_dalfo_half_life_days".to_string(), 0.5); // ~12 hours (quinupristin/dalfopristin)
        
        // Folate antagonists
        map.insert("drug_trim_sulf_half_life_days".to_string(), 0.5); // ~12 hours (trimethoprim)
        
        // Other antibiotics
        map.insert("drug_chlorampheni_half_life_days".to_string(), 0.125); // ~3 hours
        map.insert("drug_nitrofurantoin_half_life_days".to_string(), 0.017); // ~20 minutes
        map.insert("drug_retapamulin_half_life_days".to_string(), 0.25); // ~6 hours (topical, limited data)
        map.insert("drug_fusidic_a_half_life_days".to_string(), 0.375); // ~9 hours
        map.insert("drug_metronidazole_half_life_days".to_string(), 0.33); // ~8 hours
        map.insert("drug_furazolidone_half_life_days".to_string(), 0.25); // ~6 hours
        map.insert("already_on_drug_initiation_multiplier".to_string(), 1.000); // 0.0001
        map.insert("double_dose_probability_if_identified_infection".to_string(), 0.1); // Probability for double dose
        
        // Global Immune System Parameters
        map.insert("immune_decay_rate_per_day".to_string(), 0.02); // Rate at which immunity decays when not actively fighting infection
      

        for &drug in DRUG_SHORT_NAMES.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                map.insert(format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug, bacteria), 1.0); // 0.0
                map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.01); // 0.05
                map.insert(format!("drug_{}_for_bacteria_{}_resistance_emergence_rate_per_day_baseline", drug, bacteria), 0.8); // 0.0001
            }
        }


        // todo: for each drug-bacteria combination will need a specific multiplier for initiation rate
        // will need changes also in mod.rs 

        map.insert("random_drug_cessation_probability".to_string(), 0.03); // Probability an individual randomly stops a drug per day

        // General Acquisition & Resistance Parameters
        // this two below will need to change over calendar time - for the hospital acquired may decide to sample from 
        // majority_r of people in hospital with the bacteria  
        map.insert("environmental_majority_r_level_for_new_acquisition".to_string(), 0.0);
        map.insert("hospital_majority_r_level_for_new_acquisition".to_string(), 0.0);

        map.insert("max_resistance_level".to_string(), 1.0);
        map.insert("majority_r_evolution_rate_per_day_when_drug_present".to_string(), 0.001);

        // Resistance Emergence and Decay Parameters
        map.insert("resistance_emergence_rate_per_day_baseline".to_string(), 0.01);  // 0.000001 Baseline probability for de novo resistance emergence
        map.insert("microbiome_resistance_emergence_rate_per_day_baseline".to_string(), 0.005); // Separate baseline for microbiome resistance emergence
        map.insert("resistance_emergence_bacteria_level_multiplier".to_string(), 0.05); // Multiplier for bacteria level's effect on emergence
        map.insert("any_r_emergence_level_on_first_emergence".to_string(), 0.5); // The resistance level 'any_r' starts at upon emergence

        
        //  Microbiome Resistance Transfer Parameter
        map.insert("microbiome_resistance_transfer_probability_per_day".to_string(), 0.05); // Probability per day for resistance transfer between infection and microbiome
    

        // Testing Parameters
        map.insert("test_delay_days".to_string(), 3.0);
        map.insert("test_rate_per_day".to_string(), 0.20);  // 0.15

        // --- Test result and test_r logic parameters ---
        map.insert("prob_test_r_done".to_string(), 0.95); // Probability test is actually done (per day eligible)
        map.insert("test_r_error_probability".to_string(), 0.02); // Probability of error in test result
        map.insert("test_r_error_value".to_string(), 0.25); // Value to use for error in test_r

        // Syndrome-specific multipliers (example)
        map.insert("syndrome_3_initiation_multiplier".to_string(), 10.0); // Respiratory syndrome
        map.insert("syndrome_7_initiation_multiplier".to_string(), 8.0);  // Gastrointestinal syndrome
        map.insert("syndrome_8_initiation_multiplier".to_string(), 12.0); // Genital syndrome (example ID)        

        // Hospitalization Parameters
        map.insert("hospitalization_baseline_rate_per_day".to_string(), 0.00001); // 0.00001  Baseline daily probability of hospitalization
        map.insert("hospitalization_age_multiplier_per_day".to_string(), 0.000001); // Increase in daily hospitalization probability per year of age
        map.insert("hospitalization_recovery_rate_per_day".to_string(), 0.1); // Daily probability of recovering from hospitalization
        map.insert("hospitalization_max_days".to_string(), 30.0); // Max days in hospital before forced discharge (as fallback)

        // initiate travel
        map.insert("travel_probability_per_day".to_string(), 0.00005);



        // Default Initial Drug Levels and Double Dose Multipliers for ALL Drugs
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("drug_{}_initial_level", drug), 10.0); // Default initial level for each drug
            map.insert(format!("drug_{}_double_dose_multiplier", drug), 2.0); // Default double dose multiplier
            map.insert(format!("drug_{}_spectrum_breadth", drug), 3.0); // Default spectrum: 1.0=narrow, 5.0=very broad
        }

        // Bacterial Identification Effect Parameters
        map.insert("empiric_therapy_broad_spectrum_bonus".to_string(), 2.0); // Multiplier for broad-spectrum drugs when no bacteria identified
        map.insert("targeted_therapy_narrow_spectrum_bonus".to_string(), 3.0); // Multiplier for narrow-spectrum drugs when bacteria identified  
        map.insert("targeted_therapy_broad_spectrum_penalty".to_string(), 0.4); // Penalty for broad-spectrum drugs when bacteria identified
        map.insert("targeted_therapy_ineffective_drug_penalty".to_string(), 0.1); // Strong penalty for drugs ineffective against identified bacteria

        // Drug Spectrum Classifications (1.0=narrow, 5.0=very broad)
        map.insert("drug_penicilling_spectrum_breadth".to_string(), 2.0); // Narrow spectrum
        map.insert("drug_amoxicillin_spectrum_breadth".to_string(), 3.0); // Medium spectrum  
        map.insert("drug_azithromycin_spectrum_breadth".to_string(), 4.0); // Broad spectrum
        map.insert("drug_ciprofloxacin_spectrum_breadth".to_string(), 4.5); // Very broad spectrum
        map.insert("drug_trim_sulf_spectrum_breadth".to_string(), 3.5); // Medium-broad spectrum
        map.insert("drug_meropenem_spectrum_breadth".to_string(), 5.0); // Very broad spectrum (carbapenem)
        map.insert("drug_cefepime_spectrum_breadth".to_string(), 4.0); // Broad spectrum (4th gen cephalosporin)
        map.insert("drug_vancomycin_spectrum_breadth".to_string(), 2.5); // Narrow-medium spectrum (gram-positive only)
        map.insert("drug_linezolid_spectrum_breadth".to_string(), 2.0); // Narrow spectrum (gram-positive only)
        map.insert("drug_ceftriaxone_spectrum_breadth".to_string(), 4.0); // Broad spectrum (3rd gen cephalosporin)

        // Global defaults, used if a bacteria-specific parameter is not found
        map.insert("default_sepsis_baseline_risk_per_day".to_string(), 0.00001); // Very small baseline daily risk
        map.insert("default_sepsis_level_multiplier".to_string(), 0.005); // Multiplier for bacterial level (e.g., higher level = higher risk)
        map.insert("default_sepsis_duration_multiplier".to_string(), 0.000001); // Multiplier for duration of infection (e.g., longer duration = higher risk)


        // Background Mortality Parameters (Age, Region, and Sex dependent)
        map.insert("base_background_mortality_rate_per_day".to_string(), 0.00001); // 0.000005  Example: 0.0005% chance of death per day, for a baseline individual
        map.insert("age_mortality_multiplier_per_year".to_string(), 1.01); // 0.0000001 Example: Small increase in daily death risk per year of age

        // Region-specific mortality multipliers. Ensure these match your `Region` enum variants.
        map.insert("northamerica_mortality_multiplier".to_string(), 1.0);
        map.insert("southamerica_mortality_multiplier".to_string(), 1.0);
        map.insert("africa_mortality_multiplier".to_string(), 1.2);   
        map.insert("asia_mortality_multiplier".to_string(), 1.1);
        map.insert("europe_mortality_multiplier".to_string(), 0.9);     
        map.insert("oceania_mortality_multiplier".to_string(), 1.0);    

        // Sex-specific mortality multipliers. Ensure these match your `sex_at_birth` strings.
        map.insert("male_mortality_multiplier".to_string(), 1.1);   // Example: Males have 10% higher mortality risk
        map.insert("female_mortality_multiplier".to_string(), 0.9); // Example: Females have 10% lower mortality risk


        //  Immunosuppression Onset and Recovery Rates
        map.insert("immunosuppression_onset_rate_per_day".to_string(), 0.0001);   // Probability of becoming immunosuppressed daily
        map.insert("immunosuppression_recovery_rate_per_day".to_string(), 0.0005); // Probability of recovering from immunosuppression daily


        // Sepsis Mortality Parameter (Absolute risk, independent of other factors)
        map.insert("sepsis_absolute_death_risk_per_day".to_string(), 0.1); // Example: 10% absolute chance of death per day if septic

        //  Default Toxicity Parameter
        map.insert("default_drug_toxicity_per_unit_level_per_day".to_string(), 0.005); // Adjust this default as needed

        //  Default Microbiome Acquisition Parameter
        // A multiplier for the infection acquisition probability to get microbiome acquisition probability.
        // If > 1.0, microbiome acquisition is more likely than infection for the same factors.
        // If < 1.0, microbiome acquisition is less likely.
        map.insert("default_microbiome_acquisition_multiplier".to_string(), 2.0); // Example: Microbiome acquisition is twice as likely as infection given the same exposure.

        //  Default Microbiome Clearance Parameter (from previous suggestion, ensure it's there)
        map.insert("default_microbiome_clearance_probability_per_day".to_string(), 0.01); // E.g., 1% chance to lose carriage per day

        //  Microbiome Presence Effect on Infection Acquisition
        // A multiplier for infection acquisition probability if the bacteria is already present in the microbiome.
        // Value > 1.0 means microbiome presence increases infection risk.
        // Value < 1.0 means microbiome presence decreases infection risk (e.g., due to local immunity/competition).
        map.insert("default_microbiome_infection_acquisition_multiplier".to_string(), 0.1); // Example: Much harder to get infected if already colonized.

        //  Contact and Exposure Level Parameters
        map.insert("contact_level_daily_fluctuation_range".to_string(), 0.5); // Amount of random daily fluctuation
        map.insert("min_contact_level".to_string(), 0.0); // Minimum possible contact/exposure level
        map.insert("max_contact_level".to_string(), 10.0); // Maximum possible contact/exposure level

        // Sexual Contact Parameters
        map.insert("sexual_contact_baseline".to_string(), 5.0); // Baseline level for a young adult
        map.insert("sexual_contact_age_peak_days".to_string(), 25.0 * 365.0); // Age in days (25 years)
        map.insert("sexual_contact_age_rise_exponent".to_string(), 2.0); // Controls how fast contact rises with age before peak (higher = steeper)
        map.insert("sexual_contact_age_decline_rate".to_string(), 0.00005); // Rate of decline per day after peak age (e.g., 0.00005 means ~1.8% drop per year)
        map.insert("sexual_contact_hospital_multiplier".to_string(), 0.0); 

        // Airborne Contact (Adults) Parameters
        map.insert("airborne_contact_adult_baseline".to_string(), 5.0);
        map.insert("airborne_contact_adult_age_breakpoint_days".to_string(), 18.0 * 365.0); // Age in days (18 years)
        map.insert("airborne_contact_adult_child_multiplier".to_string(), 0.2); // How much less children contact adults (vs. adult-adult baseline)
        map.insert("airborne_contact_in_hospital_multiplier".to_string(), 1.5); // May increase due to healthcare staff contact

        // Airborne Contact (Children) Parameters
        map.insert("airborne_contact_child_baseline".to_string(), 3.0);
        map.insert("airborne_contact_child_age_breakpoint_days".to_string(), 12.0 * 365.0); // Age in days (12 years)
        map.insert("airborne_contact_child_child_multiplier".to_string(), 1.5); // How much more children contact children (vs. child baseline)
        map.insert("airborne_contact_child_adult_multiplier".to_string(), 0.5); // How much less adults contact children (vs. child baseline)

        // Oral Exposure Parameters
        map.insert("oral_exposure_baseline".to_string(), 2.0);
        map.insert("oral_exposure_child_age_breakpoint_days".to_string(), 5.0 * 365.0); // Age in days (5 years)
        map.insert("oral_exposure_child_multiplier".to_string(), 3.0); // Higher for young children
        map.insert("oral_exposure_in_hospital_multiplier".to_string(), 0.8); // Slightly reduced due to hospital hygiene

        // Mosquito Exposure Parameters
        map.insert("mosquito_exposure_baseline".to_string(), 1.0);
        map.insert("mosquito_exposure_in_hospital_multiplier".to_string(), 0.2); // Significantly reduced indoors/hospital
        
        // Region-specific multipliers (example values, adjust as needed based on actual epidemiology)
        map.insert("north_america_mosquito_exposure_multiplier".to_string(), 0.5);
        map.insert("south_america_mosquito_exposure_multiplier".to_string(), 5.0);
        map.insert("africa_mosquito_exposure_multiplier".to_string(), 8.0);
        map.insert("asia_mosquito_exposure_multiplier".to_string(), 6.0);
        map.insert("europe_mosquito_exposure_multiplier".to_string(), 0.2);
        map.insert("oceania_mosquito_exposure_multiplier".to_string(), 3.0);
        // Ensure you have multipliers for all variants of your `Region` enum,
        // or add a default handling in the `mod.rs` if a region param isn't found.
        // If `Region::Home` refers to a generic home location not tied to a specific geographical region,
        // you might need to reconsider its role or default it to 1.0 or an average.

 












        // Drug-specific Adverse Event Death Risk Parameters (Absolute risk, independent and drug-specific)
        // Add an entry for each drug that can cause an adverse event leading to death.
        // The key format is "drug_{drug_short_name}_adverse_event_death_risk"
        map.insert("drug_amoxicillin_adverse_event_death_risk".to_string(), 0.0001); // Example: 0.01% daily absolute risk from amoxicillin
        map.insert("drug_meropenem_adverse_event_death_risk".to_string(), 0.0005);   // Example: 0.05% daily absolute risk from meropenem
        // Add similar lines for other drugs as needed:
        // map.insert("drug_ciprofloxacin_adverse_event_death_risk".to_string(), 0.0002);
        // map.insert("drug_vancomycin_adverse_event_death_risk".to_string(), 0.0003);


        // --- Overrides for Specific Bacteria (Customize these as needed) ---

        // acinetobac_bau Parameters
        map.insert("acinetobac_bau_acquisition_prob_baseline".to_string(), 0.001 ); // 0.2
        map.insert("acinetobac_bau_hospital_acquired_multiplier".to_string(), 5.0); // higher risk in hospital
        map.insert("acinetobac_bau_immunity_base_response".to_string(), 0.001); // 0.001
        map.insert("acinetobac_bau_immunity_increase_per_infection_day".to_string(), 0.2  );  // 0.2
        map.insert("acinetobac_bau_immunity_increase_per_unit_higher_bacteria_level".to_string(), 0.2  );  // 0.2
        map.insert("acinetobac_bau_immunity_effect_on_level_change".to_string(), 0.005  );  // 0.005
        map.insert("drug_cefepime_for_bacteria_acinetobac_bau_resistance_emergence_rate_per_day_baseline".to_string(), 0.01); // 0.2 0.0001 Baseline probability for de novo resistance emergence
        map.insert("drug_cefepime_for_bacteria_acinetobac_bau_initiation_multiplier".to_string(), 1.0); // 1000000.0



        // Drug-Specific Toxicity Parameters (Examples)
        // Format: "drug_{drug_name}_toxicity_per_unit_level_per_day"
        map.insert("drug_penicilling_toxicity_per_unit_level_per_day".to_string(), 0.002); // Lower toxicity example
        map.insert("drug_cefepime_toxicity_per_unit_level_per_day".to_string(), 0.01);    // Higher toxicity example
        map.insert("drug_meropenem_toxicity_per_unit_level_per_day".to_string(), 0.008);


        //  Bacteria-Specific Microbiome Acquisition Multipliers (Optional Examples)
        // Format: "{bacteria_name}_microbiome_acquisition_multiplier"
        map.insert("strep_pneu_microbiome_acquisition_multiplier".to_string(), 3.0); // Strep Pneu might colonize more easily
        map.insert("salm_typhi_microbiome_acquisition_multiplier".to_string(), 0.5);  // Salmonella might colonize less easily than cause infection

        //  Bacteria-Specific Microbiome Clearance Parameters (Optional Examples)
        // Format: "{bacteria_name}_microbiome_clearance_probability_per_day"
        map.insert("strep_pneu_microbiome_clearance_probability_per_day".to_string(), 0.02);
        map.insert("esch_coli_microbiome_clearance_probability_per_day".to_string(), 0.005);


        
        //  Bacteria-Specific Microbiome Infection Acquisition Multipliers (Optional Examples)
        // Format: "{bacteria_name}_microbiome_infection_acquisition_multiplier"
        map.insert("strep_pneu_microbiome_infection_acquisition_multiplier".to_string(), 0.05); // Maybe Strep Pneu is very protective when colonized
        map.insert("salm_typhi_microbiome_infection_acquisition_multiplier".to_string(), 0.8);  // Salmonella might offer less protection, or even slightly increase risk if certain strains


        // Add specific rates for your drugs (e.g., higher for more toxic drugs)
        map.insert("gentamicin_toxicity_mortality_rate_per_day".to_string(), 0.00001); // Example for a specific drug

        // Add more specific parameters for acinetobac_bau if needed

        // --- Age Effect Scaling Overrides for Specific Bacteria ---
        map.insert("strep_pneu_age_effect_scaling".to_string(), 1.2); // Stronger age effect for pneumonia
        map.insert("salm_typhi_age_effect_scaling".to_string(), 0.8); // Weaker age effect 
        map.insert("n_gonorrhoeae_age_effect_scaling".to_string(), 1.5); // Strong age effect for STI
        map.insert("acinetobac_bau_age_effect_scaling".to_string(), 1.3); // Strong age effect for nosocomial pathogen

        // Overrides for Specific Drug Initial Levels & Double Dose Multipliers
        map.insert("drug_penicilling_double_dose_multiplier".to_string(), 2.5); // Example: higher multiplier for penicillin

        map.insert("drug_amoxicillin_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_azithromycin_double_dose_multiplier".to_string(), 1.8);

        map.insert("drug_ciprofloxacin_double_dose_multiplier".to_string(), 2.2);

        map.insert("drug_trim_sulf_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_meropenem_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_cefepime_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_vancomycin_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_linezolid_double_dose_multiplier".to_string(), 2.0);

        map.insert("drug_ceftriaxone_double_dose_multiplier".to_string(), 2.0);

        /* ... (rest of your commented out bacteria parameters) ... */

        map
    };

    // --- String Parameters (for template names, etc.) ---
    pub static ref STRING_PARAMETERS: HashMap<String, String> = {
        let mut map = HashMap::new();
        
        // Default age risk templates for all bacteria
        for &bacteria in BACTERIA_LIST.iter() {
            map.insert(format!("{}_age_risk_template", bacteria), "respiratory".to_string()); // Default template
        }

        // Specific bacteria overrides - assign each bacteria to most appropriate template
        map.insert("strep_pneu_age_risk_template".to_string(), "respiratory".to_string());
        map.insert("haem_infl_age_risk_template".to_string(), "respiratory".to_string());
        map.insert("salm_typhi_age_risk_template".to_string(), "gastrointestinal".to_string());
        map.insert("esch_coli_age_risk_template".to_string(), "urogenital".to_string());
        map.insert("pseud_aerug_age_risk_template".to_string(), "bloodstream".to_string());
        map.insert("staph_aureus_age_risk_template".to_string(), "skin_soft_tissue".to_string());
        map.insert("n_gonorrhoeae_age_risk_template".to_string(), "sexually_transmitted".to_string());
        map.insert("acinetobac_bau_age_risk_template".to_string(), "bloodstream".to_string());

        map
    };
}

/// Retrieves a global simulation parameter.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_global_param(key: &str) -> Option<f64> {
    PARAMETERS.get(key).copied()
}

/// Retrieves a bacteria-specific simulation parameter.
/// It directly looks up "{bacteria_name}_{param_suffix}".
/// Because all bacteria now have explicit entries, there's no need for a "generic_bacteria_" fallback in this function.
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_bacteria_param(bacteria_name: &str, param_suffix: &str) -> Option<f64> {
    let specific_key = format!("{}_{}", bacteria_name, param_suffix);
    PARAMETERS.get(&specific_key).copied()
}
/// Retrieves a drug-specific simulation parameter.
/// It looks up "drug_{drug_name}_{param_suffix}".
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_drug_param(drug_name: &str, param_suffix: &str) -> Option<f64> {
    let specific_key = format!("drug_{}_{}", drug_name, param_suffix);
    PARAMETERS.get(&specific_key).copied()
}

// --- Age Risk Templates Configuration ---

lazy_static! {
    pub static ref AGE_RISK_TEMPLATES: HashMap<&'static str, Vec<f64>> = {
        let mut m = HashMap::new();
        
        // Age groups: 0-1, 1-5, 5-18, 18-50, 50-70, 70+
        // Values represent risk multipliers relative to baseline (18-50 age group = 1.0)
        
        m.insert("respiratory", vec![3.0, 1.8, 0.8, 1.0, 1.3, 2.5]);          // High infant/elderly risk (pneumonia, URI)
        m.insert("gastrointestinal", vec![2.5, 2.0, 1.2, 1.0, 1.1, 1.8]);    // High young child risk (diarrheal diseases)
        m.insert("urogenital", vec![1.2, 0.8, 0.9, 1.0, 1.4, 2.2]);          // Moderate elderly risk (UTIs)
        m.insert("skin_soft_tissue", vec![1.5, 1.3, 1.1, 1.0, 1.2, 1.8]);    // Mild age gradient
        m.insert("bloodstream", vec![4.0, 2.0, 0.7, 1.0, 1.5, 3.0]);         // Very high infant/elderly risk (sepsis)
        m.insert("vector_borne", vec![1.8, 1.5, 1.0, 1.0, 1.1, 1.4]);        // Moderate child/elderly risk (mosquito-borne)
        m.insert("sexually_transmitted", vec![0.1, 0.2, 0.8, 1.0, 0.8, 0.3]); // Peak in young adults
        m.insert("flat", vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);               // No age effect
        
        m
    };
}

// --- NEW: Cross-Resistance Configuration ---

lazy_static! {
    static ref CROSS_RESISTANCE_GROUPS: HashMap<&'static str, Vec<Vec<&'static str>>> = {
        let mut m = HashMap::new();

        // Example: For E. coli, group penicillins and fluoroquinolones
        m.insert("esch_coli", vec![
            vec!["penicillin", "amoxicillin", "piperacillin_tazobactam"],
            vec!["ciprofloxacin", "levofloxacin"],
        ]);

        // Example: For Strep. pneumoniae, group macrolides
        m.insert("strep_pneu", vec![
            vec!["azithromycin", "clarithromycin"],
        ]);


        m.insert("acinetobacter baumannii", vec![
            vec!["penicilling", "ampicillin", "amoxicillin"],
        ]);
        
        

        // Add other bacteria-specific groups here
        // If a bacteria is not listed, it has no cross-resistance groups.

        m
    };
}

/// Returns the cross-resistance drug groups for each bacterium.
pub fn get_cross_resistance_groups() -> &'static HashMap<&'static str, Vec<Vec<&'static str>>> {
    &CROSS_RESISTANCE_GROUPS
}

/// Retrieves a string parameter (like template names).
/// Returns `Some(value)` if found, `None` otherwise.
pub fn get_string_param(key: &str) -> Option<String> {
    STRING_PARAMETERS.get(key).cloned()
}

/// Calculates the age-based infection risk multiplier for a given bacteria and age.
/// Uses the template system with bacteria-specific scaling.
/// Returns a multiplier (1.0 = baseline risk, >1.0 = increased risk, <1.0 = decreased risk)
pub fn get_age_infection_multiplier(bacteria_name: &str, age_days: i32) -> f64 {
    let age_years = age_days as f64 / 365.0;
    
    // Determine age group index (0-5 for the six age groups)
    let age_group_idx = match age_years {
        x if x < 1.0 => 0,   // 0-1 years
        x if x < 5.0 => 1,   // 1-5 years  
        x if x < 18.0 => 2,  // 5-18 years
        x if x < 50.0 => 3,  // 18-50 years (reference group)
        x if x < 70.0 => 4,  // 50-70 years
        _ => 5,              // 70+ years
    };
    
    // Get the template name for this bacteria
    let template_key = format!("{}_age_risk_template", bacteria_name);
    let template_name = get_string_param(&template_key).unwrap_or_else(|| "respiratory".to_string());
    
    // Get the scaling factor for this bacteria
    let scaling = get_bacteria_param(bacteria_name, "age_effect_scaling").unwrap_or(1.0);
    
    // Look up the base multiplier from the template
    if let Some(template) = AGE_RISK_TEMPLATES.get(template_name.as_str()) {
        let base_multiplier = template[age_group_idx];
        // Scale the deviation from 1.0 by the scaling factor
        // scaling = 0.0 means no age effect (flat = 1.0)
        // scaling = 1.0 means full template effect
        // scaling > 1.0 means amplified age effect
        1.0 + (base_multiplier - 1.0) * scaling
    } else {
        // Fallback if template not found
        1.0
    }
}

