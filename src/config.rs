// src/config.rs
use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES}; // Import both lists

// --- Global Simulation Parameters ---
lazy_static! {
    pub static ref PARAMETERS: HashMap<String, f64> = {
        let mut map = HashMap::new();


        // --- Default Parameters for ALL Bacteria from BACTERIA_LIST ---
        // These are set first, and can then be overridden by specific entries below.
        for &bacteria in BACTERIA_LIST.iter() {
            map.insert(format!("{}_acquisition_prob_baseline", bacteria), 0.001); // 0.001
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
        map.insert("drug_base_initiation_rate_per_day".to_string(), 0.0001); // 0.0001
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
      

        // --- Drug-Bacteria Potency Matrix: Evidence-Based Approach ---
        // Instead of uniform potency, use clinically relevant potency categories:
        // 0.20+ = Excellent potency (first-line therapy)
        // 0.10-0.19 = Good potency (reliable option)
        // 0.05-0.09 = Moderate potency (situational use)
        // 0.01-0.04 = Poor potency (usually ineffective)
        // 0.005 = Very poor/no activity
        
        // Define drug classes for easier management
        let penicillins = vec!["penicilling", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin"];
        let cephalosporins_1_2 = vec!["cephalexin", "cefazolin", "cefuroxime"];
        let cephalosporins_3_4 = vec!["ceftriaxone", "ceftazidime", "cefepime", "ceftaroline"];
        let carbapenems = vec!["meropenem", "imipenem_c", "ertapenem"];
        let _monobactams = vec!["aztreonam"];
        let macrolides = vec!["erythromycin", "azithromycin", "clarithromycin"];
        let _lincosamides = vec!["clindamycin"];
        let aminoglycosides = vec!["gentamicin", "tobramycin", "amikacin"];
        let fluoroquinolones = vec!["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"];
        let _tetracyclines = vec!["tetracycline", "doxyclycline", "minocycline"];
        let glycopeptides = vec!["vancomycin", "teicoplanin"];
        let oxazolidinones = vec!["linezolid", "tedizolid"];
        let _folate_antagonists = vec!["trim_sulf"];
        let _other_antibiotics = vec!["quinu_dalfo", "chlorampheni", "nitrofurantoin", "retapamulin", "fusidic_a", "metronidazole", "furazolidone"];

        // Define bacterial groups for potency patterns
        let gram_pos_cocci = vec!["staphylococcus aureus", "streptococcus pneumoniae", "streptococcus pyogenes", "streptococcus agalactiae", "enterococcus faecalis", "enterococcus faecium"];
        let gram_neg_enterobacteria = vec!["escherichia coli", "klebsiella pneumoniae", "enterobacter spp.", "citrobacter spp.", "serratia spp.", "proteus spp.", "morganella spp.", "enterobacter_cloacae"];
        let gram_neg_non_fermenting = vec!["pseudomonas aeruginosa", "acinetobacter baumannii"];
        let _fastidious_gram_neg = vec!["haemophilus influenzae", "moraxella_catarrhalis", "neisseria gonorrhoeae", "neisseria_meningitidis"];
        let _enteric_pathogens = vec!["salmonella enterica serovar typhi", "salmonella enterica serovar paratyphi a", "invasive non-typhoidal salmonella spp.", "shigella spp.", "vibrio cholerae", "campylobacter_jejuni", "yersinia_enterocolitica"];
        let _atypical_pathogens = vec!["chlamydia trachomatis"];
        let _anaerobes_spore_formers = vec!["clostridioides_difficile"];
        let _gram_pos_rods = vec!["listeria_monocytogenes"];

        // Set default low potency for all combinations first
        for &drug in DRUG_SHORT_NAMES.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                map.insert(format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug, bacteria), 1.0);
                map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.01); // Default low potency
                map.insert(format!("drug_{}_for_bacteria_{}_resistance_emergence_rate_per_day_baseline", drug, bacteria), 0.8);
            }
        }

        // Now set specific potencies based on clinical evidence
        
        // GRAM-POSITIVE COCCI (Staph, Strep, Enterococcus)
        for &bacteria in gram_pos_cocci.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Penicillins - excellent for Strep (if sensitive), poor for Staph due to beta-lactamase
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("streptococcus") { 0.18 } else { 0.02 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Cephalosporins - good for most gram-positive (except Enterococcus)
                for &drug in cephalosporins_1_2.iter().chain(cephalosporins_3_4.iter()) {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("enterococcus") { 0.01 } else { 0.15 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Carbapenems - good but reserve for resistant cases
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("enterococcus") { 0.05 } else { 0.16 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Macrolides - good for Strep and atypicals
                for &drug in macrolides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.12);
                    }
                }
                
                // Glycopeptides - excellent for gram-positive, especially MRSA/VRE
                for &drug in glycopeptides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.20);
                    }
                }
                
                // Oxazolidinones - excellent for resistant gram-positive
                for &drug in oxazolidinones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.22);
                    }
                }
            }
        }

        // GRAM-NEGATIVE ENTEROBACTERIA (E. coli, Klebsiella, etc.)
        for &bacteria in gram_neg_enterobacteria.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Penicillins - poor except piperacillin
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if drug == "piperacillin" { 0.14 } else { 0.02 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Cephalosporins - variable by generation
                for &drug in cephalosporins_1_2.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.08);
                    }
                }
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.16);
                    }
                }
                
                // Carbapenems - excellent broad-spectrum
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.21);
                    }
                }
                
                // Fluoroquinolones - good broad-spectrum
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.17);
                    }
                }
                
                // Aminoglycosides - good for serious infections
                for &drug in aminoglycosides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.15);
                    }
                }
                
                // Trim-sulf - moderate activity
                if DRUG_SHORT_NAMES.contains(&"trim_sulf") {
                    map.insert(format!("drug_trim_sulf_for_bacteria_{}_potency_when_no_r", bacteria), 0.10);
                }
            }
        }

        // PSEUDOMONAS & ACINETOBACTER (Non-fermenting gram-negatives)
        for &bacteria in gram_neg_non_fermenting.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Most beta-lactams poor except specific anti-pseudomonal agents
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if drug == "piperacillin" { 0.13 } else { 0.005 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Only specific cephalosporins active
                for &drug in cephalosporins_1_2.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.005);
                    }
                }
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if drug == "ceftazidime" || drug == "cefepime" { 0.14 } else { 0.02 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Carbapenems - good but resistance emerging
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("acinetobacter") { 0.12 } else { 0.16 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Fluoroquinolones - good activity
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.15);
                    }
                }
                
                // Aminoglycosides - good for combination therapy
                for &drug in aminoglycosides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.14);
                    }
                }
            }
        }

        // Add specific high-potency combinations for clinical effectiveness
        // These represent particularly effective drug-bacteria pairs
        
        // Azithromycin for atypicals and some enteric pathogens
        if DRUG_SHORT_NAMES.contains(&"azithromycin") {
            for &bacteria in &["chlamydia trachomatis", "campylobacter_jejuni"] {
                if BACTERIA_LIST.contains(&bacteria) {
                    map.insert(format!("drug_azithromycin_for_bacteria_{}_potency_when_no_r", bacteria), 0.25);
                }
            }
        }
        
        // Nitrofurantoin for urinary E. coli
        if DRUG_SHORT_NAMES.contains(&"nitrofurantoin") && BACTERIA_LIST.contains(&"escherichia coli") {
            map.insert("drug_nitrofurantoin_for_bacteria_escherichia coli_potency_when_no_r".to_string(), 0.19);
        }
        
        // Metronidazole for anaerobes
        if DRUG_SHORT_NAMES.contains(&"metronidazole") && BACTERIA_LIST.contains(&"clostridioides_difficile") {
            map.insert("drug_metronidazole_for_bacteria_clostridioides_difficile_potency_when_no_r".to_string(), 0.18);
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
        
        // Region-specific travel multipliers based on income/development level
        // Higher income regions have higher outbound travel rates
        map.insert("north_america_travel_multiplier".to_string(), 3.0);  // High income, high travel
        map.insert("europe_travel_multiplier".to_string(), 3.5);         // High income, highest travel rates
        map.insert("oceania_travel_multiplier".to_string(), 2.5);        // High income, high travel
        map.insert("asia_travel_multiplier".to_string(), 1.5);          // Mixed income levels, moderate travel
        map.insert("south_america_travel_multiplier".to_string(), 0.8);  // Middle income, lower travel
        map.insert("africa_travel_multiplier".to_string(), 0.3);        // Lower income, lowest travel rates



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

        // Sepsis Risk Category Multipliers (for bacteria-specific sepsis risk)
        map.insert("high_sepsis_risk_multiplier".to_string(), 2.0);     // High-virulence pathogens (e.g., Staph aureus, Pseudomonas)
        map.insert("moderate_sepsis_risk_multiplier".to_string(), 1.0); // Moderate-virulence pathogens (default)
        map.insert("low_sepsis_risk_multiplier".to_string(), 0.3);      // Low-virulence pathogens (e.g., Chlamydia, Gonorrhea)


        // Background Mortality Parameters (Age, Region, and Sex dependent)
        map.insert("base_background_mortality_rate_per_day".to_string(), 0.00001); // 0.000005  Example: 0.0005% chance of death per day, for a baseline individual
        map.insert("age_mortality_multiplier_per_year".to_string(), 1.01); // 0.0000001 Example: Small increase in daily death risk per year of age

        // Region-specific mortality multipliers. Ensure these match your `Region` enum variants.
        map.insert("north_america_mortality_multiplier".to_string(), 1.0);
        map.insert("south_america_mortality_multiplier".to_string(), 1.0);
        map.insert("africa_mortality_multiplier".to_string(), 1.2);   
        map.insert("asia_mortality_multiplier".to_string(), 1.1);
        map.insert("europe_mortality_multiplier".to_string(), 0.9);     
        map.insert("oceania_mortality_multiplier".to_string(), 1.0);    

        // Sex-specific mortality multipliers. Ensure these match your `sex_at_birth` strings.
        map.insert("male_mortality_multiplier".to_string(), 1.1);   // Example: Males have 10% higher mortality risk
        map.insert("female_mortality_multiplier".to_string(), 0.9); // Example: Females have 10% lower mortality risk

        // Additional background mortality risk factors
        map.insert("immunosuppressed_mortality_multiplier".to_string(), 2.5); // Severely immunosuppressed individuals have higher background mortality
        map.insert("hospital_mortality_multiplier".to_string(), 1.3); // Hospitalized individuals have higher baseline mortality (proxy for comorbidities)
        map.insert("age_squared_mortality_multiplier".to_string(), 0.000001); // Additional non-linear age effect for very elderly


        //  Immunosuppression Onset and Recovery Rates
        map.insert("immunosuppression_onset_rate_per_day".to_string(), 0.0001);   // Probability of becoming immunosuppressed daily
        map.insert("immunosuppression_recovery_rate_per_day".to_string(), 0.0005); // Probability of recovering from immunosuppression daily


        // Sepsis Mortality Parameters (Age, Region, and Risk Factor dependent)
        map.insert("base_sepsis_death_risk_per_day".to_string(), 0.02); // Base 2% daily death risk for sepsis (much more realistic than 10%)
        map.insert("sepsis_age_mortality_multiplier_infant".to_string(), 3.0); // 0-1 years: much higher risk
        map.insert("sepsis_age_mortality_multiplier_child".to_string(), 0.5); // 1-18 years: lower risk  
        map.insert("sepsis_age_mortality_multiplier_adult".to_string(), 1.0); // 18-65 years: baseline risk
        map.insert("sepsis_age_mortality_multiplier_elderly".to_string(), 2.5); // 65+ years: much higher risk
        map.insert("sepsis_immunosuppressed_multiplier".to_string(), 3.0); // Immunosuppressed: 3x higher risk
        
        // Region-specific sepsis mortality multipliers (reflecting healthcare quality)
        map.insert("north_america_sepsis_mortality_multiplier".to_string(), 0.8); // Better ICU care
        map.insert("europe_sepsis_mortality_multiplier".to_string(), 0.7); // Excellent healthcare systems
        map.insert("oceania_sepsis_mortality_multiplier".to_string(), 0.8); // Good healthcare
        map.insert("asia_sepsis_mortality_multiplier".to_string(), 1.2); // Variable healthcare quality
        map.insert("south_america_sepsis_mortality_multiplier".to_string(), 1.4); // Limited ICU access
        map.insert("africa_sepsis_mortality_multiplier".to_string(), 2.0); // Limited healthcare infrastructure

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
        
        // Region-specific multipliers for mosquito exposure (example values, adjust as needed based on actual epidemiology)
        map.insert("north_america_mosquito_exposure_multiplier".to_string(), 0.5);
        map.insert("south_america_mosquito_exposure_multiplier".to_string(), 5.0);
        map.insert("africa_mosquito_exposure_multiplier".to_string(), 8.0);
        map.insert("asia_mosquito_exposure_multiplier".to_string(), 6.0);
        map.insert("europe_mosquito_exposure_multiplier".to_string(), 0.2);
        map.insert("oceania_mosquito_exposure_multiplier".to_string(), 3.0);
        
        // Region-specific bacterial infection risk multipliers
        // Based on real-world epidemiological patterns and regional prevalence
        // Format: "{region}_{bacteria_name}_infection_risk_multiplier"
        // Note: Region names use underscore format (e.g., "north_america", "south_america")
        // and bacteria names have spaces replaced with underscores
        
        // Acinetobacter baumannii - higher in tropical/subtropical regions, hospitals
        map.insert("north_america_acinetobacter_baumannii_infection_risk_multiplier".to_string(), 0.8);
        map.insert("south_america_acinetobacter_baumannii_infection_risk_multiplier".to_string(), 1.5);
        map.insert("africa_acinetobacter_baumannii_infection_risk_multiplier".to_string(), 2.0);
        map.insert("asia_acinetobacter_baumannii_infection_risk_multiplier".to_string(), 1.8);
        map.insert("europe_acinetobacter_baumannii_infection_risk_multiplier".to_string(), 0.9);
        map.insert("oceania_acinetobacter_baumannii_infection_risk_multiplier".to_string(), 1.0);
        
        // Citrobacter spp. - more common in tropical regions
        map.insert("north_america_citrobacter_spp._infection_risk_multiplier".to_string(), 0.9);
        map.insert("south_america_citrobacter_spp._infection_risk_multiplier".to_string(), 1.4);
        map.insert("africa_citrobacter_spp._infection_risk_multiplier".to_string(), 1.8);
        map.insert("asia_citrobacter_spp._infection_risk_multiplier".to_string(), 1.6);
        map.insert("europe_citrobacter_spp._infection_risk_multiplier".to_string(), 0.8);
        map.insert("oceania_citrobacter_spp._infection_risk_multiplier".to_string(), 1.1);
        
        // Enterobacter spp. - globally distributed but higher in developing regions
        map.insert("north_america_enterobacter_spp._infection_risk_multiplier".to_string(), 1.0);
        map.insert("south_america_enterobacter_spp._infection_risk_multiplier".to_string(), 1.3);
        map.insert("africa_enterobacter_spp._infection_risk_multiplier".to_string(), 1.7);
        map.insert("asia_enterobacter_spp._infection_risk_multiplier".to_string(), 1.5);
        map.insert("europe_enterobacter_spp._infection_risk_multiplier".to_string(), 0.9);
        map.insert("oceania_enterobacter_spp._infection_risk_multiplier".to_string(), 1.0);
        
        // Enterococcus faecalis - globally distributed, slightly higher in temperate regions
        map.insert("north_america_enterococcus_faecalis_infection_risk_multiplier".to_string(), 1.1);
        map.insert("south_america_enterococcus_faecalis_infection_risk_multiplier".to_string(), 1.0);
        map.insert("africa_enterococcus_faecalis_infection_risk_multiplier".to_string(), 0.9);
        map.insert("asia_enterococcus_faecalis_infection_risk_multiplier".to_string(), 1.0);
        map.insert("europe_enterococcus_faecalis_infection_risk_multiplier".to_string(), 1.2);
        map.insert("oceania_enterococcus_faecalis_infection_risk_multiplier".to_string(), 1.1);
        
        // Enterococcus faecium - higher in developed regions with heavy antibiotic use
        map.insert("north_america_enterococcus_faecium_infection_risk_multiplier".to_string(), 1.3);
        map.insert("south_america_enterococcus_faecium_infection_risk_multiplier".to_string(), 1.0);
        map.insert("africa_enterococcus_faecium_infection_risk_multiplier".to_string(), 0.7);
        map.insert("asia_enterococcus_faecium_infection_risk_multiplier".to_string(), 1.1);
        map.insert("europe_enterococcus_faecium_infection_risk_multiplier".to_string(), 1.4);
        map.insert("oceania_enterococcus_faecium_infection_risk_multiplier".to_string(), 1.2);
        
        // Escherichia coli - globally distributed, slightly higher in developing regions
        map.insert("north_america_escherichia_coli_infection_risk_multiplier".to_string(), 0.9);
        map.insert("south_america_escherichia_coli_infection_risk_multiplier".to_string(), 1.3);
        map.insert("africa_escherichia_coli_infection_risk_multiplier".to_string(), 1.6);
        map.insert("asia_escherichia_coli_infection_risk_multiplier".to_string(), 1.4);
        map.insert("europe_escherichia_coli_infection_risk_multiplier".to_string(), 0.8);
        map.insert("oceania_escherichia_coli_infection_risk_multiplier".to_string(), 1.0);
        
        // Klebsiella pneumoniae - higher in tropical/subtropical regions
        map.insert("north_america_klebsiella_pneumoniae_infection_risk_multiplier".to_string(), 0.9);
        map.insert("south_america_klebsiella_pneumoniae_infection_risk_multiplier".to_string(), 1.4);
        map.insert("africa_klebsiella_pneumoniae_infection_risk_multiplier".to_string(), 1.8);
        map.insert("asia_klebsiella_pneumoniae_infection_risk_multiplier".to_string(), 1.6);
        map.insert("europe_klebsiella_pneumoniae_infection_risk_multiplier".to_string(), 0.8);
        map.insert("oceania_klebsiella_pneumoniae_infection_risk_multiplier".to_string(), 1.1);
        
        // Pseudomonas aeruginosa - higher in humid/warm climates and developed healthcare systems
        map.insert("north_america_pseudomonas_aeruginosa_infection_risk_multiplier".to_string(), 1.1);
        map.insert("south_america_pseudomonas_aeruginosa_infection_risk_multiplier".to_string(), 1.3);
        map.insert("africa_pseudomonas_aeruginosa_infection_risk_multiplier".to_string(), 1.0);
        map.insert("asia_pseudomonas_aeruginosa_infection_risk_multiplier".to_string(), 1.2);
        map.insert("europe_pseudomonas_aeruginosa_infection_risk_multiplier".to_string(), 1.0);
        map.insert("oceania_pseudomonas_aeruginosa_infection_risk_multiplier".to_string(), 1.2);
        
        // Staphylococcus aureus - globally distributed, slightly higher in crowded/poor sanitation areas
        map.insert("north_america_staphylococcus_aureus_infection_risk_multiplier".to_string(), 0.9);
        map.insert("south_america_staphylococcus_aureus_infection_risk_multiplier".to_string(), 1.2);
        map.insert("africa_staphylococcus_aureus_infection_risk_multiplier".to_string(), 1.5);
        map.insert("asia_staphylococcus_aureus_infection_risk_multiplier".to_string(), 1.3);
        map.insert("europe_staphylococcus_aureus_infection_risk_multiplier".to_string(), 0.8);
        map.insert("oceania_staphylococcus_aureus_infection_risk_multiplier".to_string(), 1.0);
        
        // Streptococcus pneumoniae - slightly higher in cold/dry climates and crowded conditions
        map.insert("north_america_streptococcus_pneumoniae_infection_risk_multiplier".to_string(), 1.1);
        map.insert("south_america_streptococcus_pneumoniae_infection_risk_multiplier".to_string(), 1.0);
        map.insert("africa_streptococcus_pneumoniae_infection_risk_multiplier".to_string(), 1.4);
        map.insert("asia_streptococcus_pneumoniae_infection_risk_multiplier".to_string(), 1.2);
        map.insert("europe_streptococcus_pneumoniae_infection_risk_multiplier".to_string(), 1.2);
        map.insert("oceania_streptococcus_pneumoniae_infection_risk_multiplier".to_string(), 1.0);
        
        // Salmonella enterica serovar typhi - much higher in developing regions with poor sanitation
        map.insert("north_america_salmonella_enterica_serovar_typhi_infection_risk_multiplier".to_string(), 0.2);
        map.insert("south_america_salmonella_enterica_serovar_typhi_infection_risk_multiplier".to_string(), 2.0);
        map.insert("africa_salmonella_enterica_serovar_typhi_infection_risk_multiplier".to_string(), 5.0);
        map.insert("asia_salmonella_enterica_serovar_typhi_infection_risk_multiplier".to_string(), 4.0);
        map.insert("europe_salmonella_enterica_serovar_typhi_infection_risk_multiplier".to_string(), 0.1);
        map.insert("oceania_salmonella_enterica_serovar_typhi_infection_risk_multiplier".to_string(), 0.8);
        
        // Salmonella enterica serovar paratyphi a - similar pattern to typhi
        map.insert("north_america_salmonella_enterica_serovar_paratyphi_a_infection_risk_multiplier".to_string(), 0.3);
        map.insert("south_america_salmonella_enterica_serovar_paratyphi_a_infection_risk_multiplier".to_string(), 1.8);
        map.insert("africa_salmonella_enterica_serovar_paratyphi_a_infection_risk_multiplier".to_string(), 3.5);
        map.insert("asia_salmonella_enterica_serovar_paratyphi_a_infection_risk_multiplier".to_string(), 4.5);
        map.insert("europe_salmonella_enterica_serovar_paratyphi_a_infection_risk_multiplier".to_string(), 0.2);
        map.insert("oceania_salmonella_enterica_serovar_paratyphi_a_infection_risk_multiplier".to_string(), 1.0);
        
        // Invasive non-typhoidal salmonella - highest in sub-Saharan Africa
        map.insert("north_america_invasive_non-typhoidal_salmonella_spp._infection_risk_multiplier".to_string(), 0.5);
        map.insert("south_america_invasive_non-typhoidal_salmonella_spp._infection_risk_multiplier".to_string(), 1.2);
        map.insert("africa_invasive_non-typhoidal_salmonella_spp._infection_risk_multiplier".to_string(), 8.0);
        map.insert("asia_invasive_non-typhoidal_salmonella_spp._infection_risk_multiplier".to_string(), 1.5);
        map.insert("europe_invasive_non-typhoidal_salmonella_spp._infection_risk_multiplier".to_string(), 0.3);
        map.insert("oceania_invasive_non-typhoidal_salmonella_spp._infection_risk_multiplier".to_string(), 1.0);
        
        // Shigella spp. - higher in regions with poor sanitation
        map.insert("north_america_shigella_spp._infection_risk_multiplier".to_string(), 0.6);
        map.insert("south_america_shigella_spp._infection_risk_multiplier".to_string(), 1.8);
        map.insert("africa_shigella_spp._infection_risk_multiplier".to_string(), 3.0);
        map.insert("asia_shigella_spp._infection_risk_multiplier".to_string(), 2.5);
        map.insert("europe_shigella_spp._infection_risk_multiplier".to_string(), 0.4);
        map.insert("oceania_shigella_spp._infection_risk_multiplier".to_string(), 1.0);
        
        // Neisseria gonorrhoeae - varies by region with different sexual health practices
        map.insert("north_america_neisseria_gonorrhoeae_infection_risk_multiplier".to_string(), 1.2);
        map.insert("south_america_neisseria_gonorrhoeae_infection_risk_multiplier".to_string(), 1.1);
        map.insert("africa_neisseria_gonorrhoeae_infection_risk_multiplier".to_string(), 2.0);
        map.insert("asia_neisseria_gonorrhoeae_infection_risk_multiplier".to_string(), 0.8);
        map.insert("europe_neisseria_gonorrhoeae_infection_risk_multiplier".to_string(), 0.9);
        map.insert("oceania_neisseria_gonorrhoeae_infection_risk_multiplier".to_string(), 1.3);
        
        // Vibrio cholerae - much higher in regions with poor water/sanitation
        map.insert("north_america_vibrio_cholerae_infection_risk_multiplier".to_string(), 0.1);
        map.insert("south_america_vibrio_cholerae_infection_risk_multiplier".to_string(), 2.5);
        map.insert("africa_vibrio_cholerae_infection_risk_multiplier".to_string(), 6.0);
        map.insert("asia_vibrio_cholerae_infection_risk_multiplier".to_string(), 4.0);
        map.insert("europe_vibrio_cholerae_infection_risk_multiplier".to_string(), 0.05);
        map.insert("oceania_vibrio_cholerae_infection_risk_multiplier".to_string(), 1.5);
        
        // Chlamydia trachomatis - sexually transmitted, varies by region
        map.insert("north_america_chlamydia_trachomatis_infection_risk_multiplier".to_string(), 1.3);
        map.insert("south_america_chlamydia_trachomatis_infection_risk_multiplier".to_string(), 1.0);
        map.insert("africa_chlamydia_trachomatis_infection_risk_multiplier".to_string(), 1.8);
        map.insert("asia_chlamydia_trachomatis_infection_risk_multiplier".to_string(), 0.7);
        map.insert("europe_chlamydia_trachomatis_infection_risk_multiplier".to_string(), 1.1);
        map.insert("oceania_chlamydia_trachomatis_infection_risk_multiplier".to_string(), 1.4);
        
        // Campylobacter jejuni - higher in regions with poor food safety
        map.insert("north_america_campylobacter_jejuni_infection_risk_multiplier".to_string(), 0.8);
        map.insert("south_america_campylobacter_jejuni_infection_risk_multiplier".to_string(), 1.5);
        map.insert("africa_campylobacter_jejuni_infection_risk_multiplier".to_string(), 2.2);
        map.insert("asia_campylobacter_jejuni_infection_risk_multiplier".to_string(), 1.8);
        map.insert("europe_campylobacter_jejuni_infection_risk_multiplier".to_string(), 0.9);
        map.insert("oceania_campylobacter_jejuni_infection_risk_multiplier".to_string(), 1.1);
        
        // Add region-specific multipliers for remaining bacteria types
        // Using more conservative variations for less well-studied regional patterns
        
        // Morganella spp.
        map.insert("north_america_morganella_spp._infection_risk_multiplier".to_string(), 1.0);
        map.insert("south_america_morganella_spp._infection_risk_multiplier".to_string(), 1.2);
        map.insert("africa_morganella_spp._infection_risk_multiplier".to_string(), 1.4);
        map.insert("asia_morganella_spp._infection_risk_multiplier".to_string(), 1.3);
        map.insert("europe_morganella_spp._infection_risk_multiplier".to_string(), 0.9);
        map.insert("oceania_morganella_spp._infection_risk_multiplier".to_string(), 1.0);
        
        // Proteus spp.
        map.insert("north_america_proteus_spp._infection_risk_multiplier".to_string(), 0.9);
        map.insert("south_america_proteus_spp._infection_risk_multiplier".to_string(), 1.3);
        map.insert("africa_proteus_spp._infection_risk_multiplier".to_string(), 1.6);
        map.insert("asia_proteus_spp._infection_risk_multiplier".to_string(), 1.4);
        map.insert("europe_proteus_spp._infection_risk_multiplier".to_string(), 0.8);
        map.insert("oceania_proteus_spp._infection_risk_multiplier".to_string(), 1.0);
        
        // Serratia spp.
        map.insert("north_america_serratia_spp._infection_risk_multiplier".to_string(), 1.0);
        map.insert("south_america_serratia_spp._infection_risk_multiplier".to_string(), 1.3);
        map.insert("africa_serratia_spp._infection_risk_multiplier".to_string(), 1.5);
        map.insert("asia_serratia_spp._infection_risk_multiplier".to_string(), 1.4);
        map.insert("europe_serratia_spp._infection_risk_multiplier".to_string(), 0.9);
        map.insert("oceania_serratia_spp._infection_risk_multiplier".to_string(), 1.0);
        
        // Default multiplier for Home region and any missing region-bacteria combinations
        map.insert("home_infection_risk_multiplier_default".to_string(), 1.0);
        
        // Region-specific drug availability multipliers
        // Format: "{region}_drug_{drug_name}_availability"
        // Values: 1.0 = fully available, 0.5 = limited availability, 0.0 = not available
        // Based on realistic antibiotic access patterns across different healthcare systems
        
        // North America - Full access to most antibiotics
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("north_america_drug_{}_availability", drug), 1.0);
        }
        
        // Europe - Full access to most antibiotics
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("europe_drug_{}_availability", drug), 1.0);
        }
        
        // Asia - Good access to most drugs, some newer drugs may be limited
        for &drug in DRUG_SHORT_NAMES.iter() {
            let availability = match drug {
                // Newer, expensive drugs may have limited availability
                "tedizolid" | "ceftaroline" => 0.3,
                "teicoplanin" => 0.7, // More available in Asia than tedizolid
                _ => 1.0, // Most other drugs widely available
            };
            map.insert(format!("asia_drug_{}_availability", drug), availability);
        }
        
        // Oceania - Generally good access, similar to developed regions
        for &drug in DRUG_SHORT_NAMES.iter() {
            let availability = match drug {
                "tedizolid" | "ceftaroline" => 0.5, // Somewhat limited
                _ => 1.0,
            };
            map.insert(format!("oceania_drug_{}_availability", drug), availability);
        }
        
        // South America - Variable access, newer/expensive drugs limited
        for &drug in DRUG_SHORT_NAMES.iter() {
            let availability = match drug {
                // Very limited access to newest drugs
                "tedizolid" | "ceftaroline" => 0.1,
                "teicoplanin" => 0.3,
                "linezolid" => 0.5,
                // Limited access to some carbapenems
                "ertapenem" => 0.6,
                "meropenem" | "imipenem_c" => 0.7,
                // Moderate access to some newer cephalosporins
                "cefepime" => 0.8,
                // Good access to older, established drugs
                _ => 1.0,
            };
            map.insert(format!("south_america_drug_{}_availability", drug), availability);
        }
        
        // Africa - Most limited access, mainly basic antibiotics available
        for &drug in DRUG_SHORT_NAMES.iter() {
            let availability = match drug {
                // Basic penicillins - widely available
                "penicilling" | "ampicillin" | "amoxicillin" => 1.0,
                // Basic cephalosporins - good availability
                "cephalexin" | "cefazolin" => 0.9,
                "cefuroxime" => 0.7,
                // Third-generation cephalosporins - limited
                "ceftriaxone" => 0.6,
                "ceftazidime" => 0.4,
                // Basic macrolides and fluoroquinolones - moderate availability
                "erythromycin" | "azithromycin" => 0.8,
                "ciprofloxacin" => 0.7,
                "levofloxacin" => 0.5,
                // Aminoglycosides - basic ones available
                "gentamicin" => 0.8,
                "tobramycin" | "amikacin" => 0.4,
                // Older drugs - generally available
                "tetracycline" | "doxyclycline" => 0.9,
                "trim_sulf" => 0.9,
                "chlorampheni" => 0.8,
                "metronidazole" => 0.9,
                // Vancomycin - very limited
                "vancomycin" => 0.3,
                // Newer/expensive drugs - very limited or unavailable
                "meropenem" | "imipenem_c" => 0.2,
                "ertapenem" => 0.1,
                "linezolid" => 0.1,
                "tedizolid" | "ceftaroline" | "teicoplanin" => 0.0,
                "aztreonam" => 0.1,
                "cefepime" => 0.3,
                "moxifloxacin" => 0.2,
                "minocycline" => 0.4,
                "quinu_dalfo" => 0.1,
                "nitrofurantoin" => 0.6,
                "retapamulin" | "fusidic_a" => 0.2,
                "furazolidone" => 0.3,
                // Default for any remaining drugs
                _ => 0.1,
            };
            map.insert(format!("africa_drug_{}_availability", drug), availability);
        }
        
        // Home region - use availability based on region_living
        // This will be handled in the drug initiation logic
        for &drug in DRUG_SHORT_NAMES.iter() {
            map.insert(format!("home_drug_{}_availability", drug), 1.0); // Default fallback
        }
        
        // Ensure you have multipliers for all variants of your `Region` enum,
        // or add a default handling in the `mod.rs` if a region param isn't found.
        // If `Region::Home` refers to a generic home location not tied to a specific geographical region,
        // you might need to reconsider its role or default it to 1.0 or an average.
        
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

/// Checks if a drug is available in a given region.
/// Returns the availability multiplier (0.0 = not available, 1.0 = fully available).
/// For Home region, uses the individual's region_living.
pub fn get_drug_availability(drug_name: &str, region: &str, region_living: Option<&str>) -> f64 {
    // Handle Home region by using region_living
    let effective_region = if region == "home" {
        region_living.unwrap_or("north_america") // Default fallback if region_living not provided
    } else {
        region
    };
    
    let availability_key = format!("{}_drug_{}_availability", effective_region, drug_name);
    PARAMETERS.get(&availability_key).copied().unwrap_or(1.0) // Default to available if not specified
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

// --- CROSS-RESISTANCE CONFIGURATION ---
// NOTE: These groups are DIFFERENT from the potency drug classes above!
// Potency classes = therapeutic effectiveness groupings
// Cross-resistance groups = resistance mechanism groupings (bacteria-specific)

lazy_static! {
    static ref CROSS_RESISTANCE_GROUPS: HashMap<&'static str, Vec<Vec<&'static str>>> = {
        let mut m = HashMap::new();

        // E. coli resistance patterns
        m.insert("escherichia coli", vec![
            // ESBL resistance affects penicillins + some cephalosporins
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin"],
            // Fluoroquinolone resistance (often ciprofloxacin + levofloxacin together)
            vec!["ciprofloxacin", "levofloxacin"],
            // Aminoglycoside resistance (often linked)
            vec!["gentamicin", "tobramycin"],
        ]);

        // Acinetobacter baumannii resistance patterns
        m.insert("acinetobacter baumannii", vec![
            // β-lactamase affects most β-lactams
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime"],
            // Carbapenemase affects carbapenems
            vec!["meropenem", "imipenem_c", "ertapenem"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin"],
            // Aminoglycoside resistance
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // Klebsiella pneumoniae resistance patterns  
        m.insert("klebsiella pneumoniae", vec![
            // ESBL resistance
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "ceftriaxone"],
            // Carbapenemase (KPC, NDM, etc.)
            vec!["meropenem", "imipenem_c", "ertapenem"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin"],
        ]);

        // Streptococcus pneumoniae resistance patterns
        m.insert("streptococcus pneumoniae", vec![
            // Macrolide resistance (erm genes affect all macrolides)
            vec!["erythromycin", "azithromycin", "clarithromycin"],
            // Penicillin resistance (affects β-lactams)
            vec!["penicilling", "ampicillin", "amoxicillin"],
        ]);

        // Staphylococcus aureus resistance patterns
        m.insert("staphylococcus aureus", vec![
            // β-lactamase affects penicillins
            vec!["penicilling", "ampicillin", "amoxicillin"],
            // MRSA affects most β-lactams
            vec!["cephalexin", "cefazolin", "cefuroxime", "ceftriaxone"],
            // Macrolide-lincosamide resistance
            vec!["erythromycin", "azithromycin", "clarithromycin", "clindamycin"],
        ]);

        // Pseudomonas aeruginosa resistance patterns
        m.insert("pseudomonas aeruginosa", vec![
            // β-lactamase affects multiple β-lactams
            vec!["piperacillin", "ceftazidime", "cefepime"],
            // Carbapenemase
            vec!["meropenem", "imipenem_c"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin"],
            // Aminoglycoside resistance
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // Enterobacter species resistance patterns
        m.insert("enterobacter spp.", vec![
            // AmpC β-lactamase (chromosomal)
            vec!["ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime"],
            // ESBL if acquired
            vec!["ceftriaxone", "ceftazidime", "cefepime"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin"],
        ]);

        // Add more bacteria as needed...
        // The key insight: resistance groupings are bacteria-specific and mechanism-based,
        // while potency groupings are based on therapeutic similarity.

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

/// Gets the sepsis risk category multiplier for a bacteria.
/// Categorizes bacteria into high/moderate/low sepsis risk groups.
/// Returns the appropriate risk multiplier.
pub fn get_bacteria_sepsis_risk_multiplier(bacteria_name: &str) -> f64 {
    // High sepsis risk: bloodstream pathogens, highly virulent
    let high_risk_bacteria = [
        "staphylococcus aureus",
        "pseudomonas aeruginosa", 
        "acinetobacter baumannii",
        "enterococcus faecium",
        "streptococcus pneumoniae",
        "enterobacter spp.",
        "klebsiella pneumoniae"
    ];
    
    // Low sepsis risk: less invasive, more localized infections
    let low_risk_bacteria = [
        "chlamydia trachomatis",
        "neisseria gonorrhoeae",
        "campylobacter_jejuni",
        "shigella spp.",
        "moraxella_catarrhalis",
        "haemophilus influenzae"
    ];
    
    if high_risk_bacteria.contains(&bacteria_name) {
        get_global_param("high_sepsis_risk_multiplier").unwrap_or(2.0)
    } else if low_risk_bacteria.contains(&bacteria_name) {
        get_global_param("low_sepsis_risk_multiplier").unwrap_or(0.3)
    } else {
        // Default to moderate risk for all other bacteria
        get_global_param("moderate_sepsis_risk_multiplier").unwrap_or(1.0)
    }
}

