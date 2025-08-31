/* 



*/


//
// Centralized configuration and parameter management for the AMR simulation.
//
// Contains:
//   - Initialization of global, bacteria-specific, and drug-specific parameters
//   - Functions for parameter lookup and cross-resistance group management
//   - Age-specific vaccination, HGT, and other model parameters
//   - Reference for template and override logic
//


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
            map.insert(format!("{}_initial_infection_level", bacteria), 0.01); // 0.01 // bacteria level at initial infection
            map.insert(format!("{}_environmental_acquisition_proportion", bacteria), 0.1); // 0.1  // proportion of new infections from environment
            map.insert(format!("{}_base_bacteria_level_change", bacteria), 0.5); // 0.2 // base change in bacteria level per day
            map.insert(format!("{}_max_level", bacteria), 5.0); // max bacteria level (arbitrary standardized scale)
            map.insert(format!("{}_immunity_effect_on_level_change", bacteria), 0.08); // 0.1  0.005 / 0.05 is strong effect // effect of the immune response on bacteria level
            map.insert(format!("{}_immunity_base_response", bacteria), 0.03); // 0.01  0.001 // base immune response
            map.insert(format!("{}_immunity_increase_per_unit_higher_bacteria_level", bacteria), 0.05); // effect of bacteria level on immune response 
            map.insert(format!("{}_immunity_increase_per_infection_day", bacteria), 0.1); // effect of time infected on immune response
            map.insert(format!("{}_immunity_age_modifier", bacteria), 1.0); // effect of age on immune response
            map.insert(format!("{}_immunity_immunodeficiency_modifier", bacteria), 0.1); // effect of being immunodeficient on immune response 
            map.insert(format!("{}_max_immune_response", bacteria), 10.0); // Maximum immune response level (arbitrary scale)
             
            // Age-related bactera-specific infection risk parameters
            map.insert(format!("{}_age_effect_scaling", bacteria), 1.0); // Scale the template effect (1.0 = full effect)

            // --- Age-specific daily vaccination probability parameters for bacterial vaccines only ---
            // Only vaccines targeting bacteria in our BACTERIA_LIST
            // Age groups: 0-1, 1-5, 5-18, 18-50, 50-70, 70+
            let bacterial_vaccines = vec![
                ("pneumococcal", 1977), // PCV first licensed in 1977 (earlier polysaccharide vaccines)
                ("meningococcal", 1981), // First meningococcal vaccine licensed in 1981
                ("hib", 1985),           // Haemophilus influenzae type b vaccine licensed in 1985
            ];
            let age_groups = vec!["0_1", "1_5", "5_18", "18_50", "50_70", "70plus"];
            
            for (vaccine, availability_year) in &bacterial_vaccines {
                for age in &age_groups {
                    // Default: 0.0, user should override as needed
                    map.insert(format!("vaccine_{}_daily_prob_age_{}", vaccine, age), 0.0);
                }
                // Store vaccine availability year for historical modeling
                map.insert(format!("vaccine_{}_availability_year", vaccine), *availability_year as f64);
            }
        }

        // --- HGT Probabilities for All Donor-Recipient Bacteria Pairs ---
        for &donor in BACTERIA_LIST.iter() {
            for &recipient in BACTERIA_LIST.iter() {
                if donor != recipient {
                    // Default HGT probability (adjust as needed)
                    map.insert(format!("hgt_prob_{}_to_{}", donor, recipient), 0.0001);
                }
            }
        }



        // General Drug Parameters
        map.insert("drug_base_initiation_rate_per_day".to_string(), 0.0005 ); // 0.0001
        map.insert("drug_infection_present_multiplier".to_string(), 300.0); // 1000
        map.insert("drug_test_identified_multiplier".to_string(), 2.0);
        map.insert("drug_decay_per_day".to_string(), 1.0); // Legacy parameter - now using drug-specific half-lives
        
        // Drug-specific half-lives (in days) for realistic pharmacokinetics
        // Beta-lactam/beta-lactamase inhibitor combinations
        map.insert("drug_amoxicillin_clavulanate_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_piperacillin_tazobactam_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ampicillin_sulbactam_half_life_days".to_string(), 0.04); // ~1 hour
        map.insert("drug_ticarcillin_clavulanate_half_life_days".to_string(), 0.046); // ~1.1 hours
        map.insert("drug_ceftazidime_avibactam_half_life_days".to_string(), 0.08); // ~2 hours
        map.insert("drug_meropenem_vaborbactam_half_life_days".to_string(), 0.04); // ~1 hour

        // Polymyxins (Colistin)
        map.insert("drug_colistin_half_life_days".to_string(), 0.08); // ~2 hours

        // Colistin parameters (grouped with other drugs)
        map.insert("drug_colistin_spectrum_breadth".to_string(), 4.0); // Broad spectrum (mainly Gram-negative)
        map.insert("drug_colistin_toxicity_per_unit_level_per_day".to_string(), 0.02); // Higher toxicity
        // Colistin: higher risk of death from toxicity
        map.insert("drug_colistin_toxicity_death_risk_per_day".to_string(), 0.002); // 0.2% daily risk (higher than default)
        // Regional availability (assume widely available, adjust as needed)
        map.insert("north_america_drug_colistin_availability".to_string(), 1.0);
        map.insert("europe_drug_colistin_availability".to_string(), 1.0);
        map.insert("asia_drug_colistin_availability".to_string(), 1.0);
        map.insert("oceania_drug_colistin_availability".to_string(), 1.0);
        map.insert("south_america_drug_colistin_availability".to_string(), 1.0);
        map.insert("africa_drug_colistin_availability".to_string(), 1.0);
        map.insert("home_drug_colistin_availability".to_string(), 1.0);
        
        // Sulfonamides (first antibiotics)
        map.insert("drug_sulfanilamide_half_life_days".to_string(), 0.45); // ~11 hours
        
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
        
        // Drug Evaluation Timing Parameters
        map.insert("drug_evaluation_days_post_infection".to_string(), 7.0); // Number of days after infection to evaluate drug initiation
      

        // --- Drug-Bacteria Potency Matrix: Evidence-Based Approach ---
        // Instead of uniform potency, use clinically relevant potency categories:
        // 1.00+ = Excellent potency (first-line therapy)
        // 0.50-0.99 = Good potency (reliable option)
        // 0.25-0.49 = Moderate potency (situational use)
        // 0.05-0.24 = Poor potency (usually ineffective)
        // 0.05 = Very poor/no activity
        
        // Define drug classes for easier management
        // Polymyxins (currently only Colistin)
        let polymyxins = vec!["colistin"];
        let penicillins = vec!["penicilling", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin",
            // BL/BLI combinations
            "amoxicillin_clavulanate", "piperacillin_tazobactam", "ampicillin_sulbactam", "ticarcillin_clavulanate"
        ];
        let cephalosporins_1_2 = vec!["cephalexin", "cefazolin", "cefuroxime"];
        let _cephalosporins_3_4 = vec!["ceftriaxone", "ceftazidime", "cefepime", "ceftaroline"];
        let cephalosporins_3_4 = vec!["ceftriaxone", "ceftazidime", "cefepime", "ceftaroline",
            // BL/BLI cephalosporin
            "ceftazidime_avibactam"
        ];
        // let _carbapenems = vec!["meropenem", "imipenem_c", "ertapenem"];
        let carbapenems = vec!["meropenem", "imipenem_c", "ertapenem",
            // BL/BLI carbapenem
            "meropenem_vaborbactam"
        ];
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

        for &drug in DRUG_SHORT_NAMES.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                map.insert(format!("drug_{}_for_bacteria_{}_initiation_multiplier", drug, bacteria), 1.0);
                map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.1); // Default low potency 0.1 
                map.insert(format!("drug_{}_for_bacteria_{}_resistance_emergence_rate_per_day_baseline", drug, bacteria), 0.005);  // 0.005
            } 
        }

        // examples of how to override potencies
        // Specific override: sulfanilamide has high potency against haemophilus influenzae
        // map.insert("drug_sulfanilamide_for_bacteria_haemophilus influenzae_potency_when_no_r".to_string(), 0.85); // Example high potency
        
        // Sulfanilamide - historically effective against specific pathogens
        for &drug in DRUG_SHORT_NAMES.iter() {
            for &bacteria in BACTERIA_LIST.iter() {
                if drug == "sulfanilamide" {
                    let potency = match bacteria {
                        // Excellent against streptococci (historical primary indication)
                        bacteria if bacteria.contains("streptococcus") => 0.85,
                        // Good against E. coli (UTI treatment)
                        "escherichia coli" => 0.65,
                        // Moderate against other gram-positives
                        "staphylococcus aureus" => 0.45,
                        // Limited against enterococci and most gram-negatives
                        _ => 0.20,
                    };
                    map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                }
            }
        }

        // Set specific potencies based on clinical evidence
        
        // GRAM-POSITIVE COCCI (Staph, Strep, Enterococcus)
        for &bacteria in gram_pos_cocci.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Penicillins - excellent for Strep (if sensitive), poor for Staph due to beta-lactamase
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("streptococcus") { 0.90 } else { 0.10 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Cephalosporins - good for most gram-positive (except Enterococcus)
                for &drug in cephalosporins_1_2.iter().chain(cephalosporins_3_4.iter()) {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("enterococcus") { 0.05 } else { 0.75 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Carbapenems - good but reserve for resistant cases
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("enterococcus") { 0.25 } else { 0.80 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Macrolides - good for Strep and atypicals
                for &drug in macrolides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.60);
                    }
                }
                
                // Glycopeptides - excellent for gram-positive, especially MRSA/VRE
                for &drug in glycopeptides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 1.00);
                    }
                }
                
                // Oxazolidinones - excellent for resistant gram-positive
                for &drug in oxazolidinones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 1.10);
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
                        let potency = if drug == "piperacillin" { 0.70 } else { 0.10 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Cephalosporins - variable by generation
                for &drug in cephalosporins_1_2.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.40);
                    }
                }
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.80);
                    }
                }
                
                // Carbapenems - excellent broad-spectrum
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 1.05);
                    }
                }
                // Polymyxins (Colistin) - only active against Gram-negatives
                for &drug in polymyxins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = match bacteria {
                            // Gram-positive bacteria - intrinsically resistant (colistin doesn't penetrate Gram-positive cell walls)
                            "enterococcus faecalis" | "enterococcus faecium" | "staphylococcus aureus" | 
                            "streptococcus pneumoniae" | "streptococcus pyogenes" | "streptococcus agalactiae" |
                            "listeria_monocytogenes" | "clostridioides_difficile" => 0.0,
                            
                            // Gram-negative intrinsic resistance
                            "morganella spp." | "proteus spp." | "serratia spp." => 0.0,
                            
                            // Gram-negative with variable/reduced susceptibility
                            "salmonella enterica serovar typhi" | "salmonella enterica serovar paratyphi a" | 
                            "invasive non-typhoidal salmonella spp." | "shigella spp." => 0.5,
                            
                            // Gram-negative normally susceptible
                            "vibrio cholerae" | "yersinia_enterocolitica" => 1.0,
                            
                            // Most other Gram-negatives (E. coli, Klebsiella, Pseudomonas, Acinetobacter, etc.)
                            _ => 1.0,
                        };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Fluoroquinolones - good broad-spectrum
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.85);
                    }
                }
                
                // Aminoglycosides - good for serious infections
                for &drug in aminoglycosides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.75);
                    }
                }
                
                // Trim-sulf - moderate activity
                if DRUG_SHORT_NAMES.contains(&"trim_sulf") {
                    map.insert(format!("drug_trim_sulf_for_bacteria_{}_potency_when_no_r", bacteria), 0.50);
                }
            }
        }

        // PSEUDOMONAS & ACINETOBACTER (Non-fermenting gram-negatives)
        for &bacteria in gram_neg_non_fermenting.iter() {
            if BACTERIA_LIST.contains(&bacteria) {
                // Most beta-lactams poor except specific anti-pseudomonal agents
                for &drug in penicillins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if drug == "piperacillin" { 0.65 } else { 0.025 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Only specific cephalosporins active
                for &drug in cephalosporins_1_2.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.025);
                    }
                }
                for &drug in cephalosporins_3_4.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if drug == "ceftazidime" || drug == "cefepime" { 0.70 } else { 0.10 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                
                // Carbapenems - good but resistance emerging
                for &drug in carbapenems.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        let potency = if bacteria.contains("acinetobacter") { 0.60 } else { 0.80 };
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), potency);
                    }
                }
                // Polymyxins (Colistin) - high potency for non-fermenters
                for &drug in polymyxins.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 1.0);
                    }
                }
                
                // Fluoroquinolones - good activity
                for &drug in fluoroquinolones.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.75);
                    }
                }
                
                // Aminoglycosides - good for combination therapy
                for &drug in aminoglycosides.iter() {
                    if DRUG_SHORT_NAMES.contains(&drug) {
                        map.insert(format!("drug_{}_for_bacteria_{}_potency_when_no_r", drug, bacteria), 0.70);
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
                    map.insert(format!("drug_azithromycin_for_bacteria_{}_potency_when_no_r", bacteria), 1.25);
                }
            }
        }
        
        // Nitrofurantoin for urinary E. coli
        if DRUG_SHORT_NAMES.contains(&"nitrofurantoin") && BACTERIA_LIST.contains(&"escherichia coli") {
            map.insert("drug_nitrofurantoin_for_bacteria_escherichia coli_potency_when_no_r".to_string(), 0.95);
        }
        
        // Metronidazole for anaerobes
        if DRUG_SHORT_NAMES.contains(&"metronidazole") && BACTERIA_LIST.contains(&"clostridioides_difficile") {
            map.insert("drug_metronidazole_for_bacteria_clostridioides_difficile_potency_when_no_r".to_string(), 0.90);
        }





        // for each drug-bacteria combination will need a specific multiplier for initiation rate
        // will need changes also in mod.rs 

        map.insert("random_drug_cessation_probability".to_string(), 0.03); // Probability an individual randomly stops a drug per day
        map.insert("random_drug_cessation_probability_if_no_active_infection".to_string(), 0.2); // Higher probability if no active infection

        // General Acquisition & Resistance Parameters
        // --- Logistic Model Parameters for Infection and Microbiome Acquisition ---
        // Infection acquisition (site infection)
        map.insert("acquisition_log_odds_baseline".to_string(), -13.5); // -13.5 Default baseline log-odds for infection acquisition
        map.insert("log_odds_sexual_contact_per_unit".to_string(), 0.10); // Per unit sexual contact
        map.insert("log_odds_airborne_adult_contact_per_unit".to_string(), 0.08); // Per unit airborne adult contact
        map.insert("log_odds_airborne_child_contact_per_unit".to_string(), 0.08); // Per unit airborne child contact
        map.insert("log_odds_oral_exposure_per_unit".to_string(), 0.06); // Per unit oral exposure
        map.insert("log_odds_vaccinated".to_string(), -2.0); // Vaccination reduces log-odds
        map.insert("log_odds_microbiome_present".to_string(), 0.5); // Microbiome presence effect (example)
        map.insert("log_odds_hospital_acquired".to_string(), 2.0); // Hospital-acquired effect (default/fallback)
        
        // Bacteria-specific hospital acquisition log-odds (healthcare-associated infection risk)
        // These parameters reflect clinical reality where certain bacteria have much higher
        // acquisition risk in hospital settings due to:
        // - Environmental persistence (Acinetobacter, C. diff)
        // - Device-associated transmission (Pseudomonas, Klebsiella)
        // - Healthcare worker transmission (MRSA, VRE)
        // - Antibiotic selection pressure (C. diff, ESBL producers)
        // Values are log-odds multipliers: exp(3.0) = 20x higher risk, exp(2.0) = 7x higher, etc.
        
        // High HAI risk bacteria (major hospital pathogens)
        map.insert("acinetobacter_baumannii_log_odds_hospital_acquired".to_string(), 3.0); // 20x higher risk (exp(3.0) ≈ 20)
        map.insert("pseudomonas_aeruginosa_log_odds_hospital_acquired".to_string(), 2.5); // 12x higher risk
        map.insert("enterococcus_faecium_log_odds_hospital_acquired".to_string(), 2.8); // 16x higher risk (VRE)
        map.insert("staphylococcus_aureus_log_odds_hospital_acquired".to_string(), 2.3); // 10x higher risk (MRSA)
        map.insert("clostridioides_difficile_log_odds_hospital_acquired".to_string(), 3.2); // 25x higher risk (C. diff)
        map.insert("klebsiella_pneumoniae_log_odds_hospital_acquired".to_string(), 2.0); // 7x higher risk
        map.insert("enterobacter_spp._log_odds_hospital_acquired".to_string(), 2.2); // 9x higher risk
        map.insert("enterobacter_cloacae_log_odds_hospital_acquired".to_string(), 2.2); // 9x higher risk
        map.insert("serratia_spp._log_odds_hospital_acquired".to_string(), 2.0); // 7x higher risk
        map.insert("citrobacter_spp._log_odds_hospital_acquired".to_string(), 1.8); // 6x higher risk
        
        // Moderate HAI risk bacteria
        map.insert("escherichia_coli_log_odds_hospital_acquired".to_string(), 1.5); // 4.5x higher risk (device-associated)
        map.insert("enterococcus_faecalis_log_odds_hospital_acquired".to_string(), 1.6); // 5x higher risk
        map.insert("streptococcus_pneumoniae_log_odds_hospital_acquired".to_string(), 1.2); // 3.3x higher risk
        map.insert("proteus_spp._log_odds_hospital_acquired".to_string(), 1.4); // 4x higher risk
        map.insert("morganella_spp._log_odds_hospital_acquired".to_string(), 1.3); // 3.7x higher risk
        map.insert("listeria_monocytogenes_log_odds_hospital_acquired".to_string(), 1.0); // 2.7x higher risk
        map.insert("neisseria_meningitidis_log_odds_hospital_acquired".to_string(), 0.8); // 2.2x higher risk
        map.insert("streptococcus_pyogenes_log_odds_hospital_acquired".to_string(), 1.1); // 3x higher risk
        map.insert("streptococcus_agalactiae_log_odds_hospital_acquired".to_string(), 1.2); // 3.3x higher risk
        map.insert("haemophilus_influenzae_log_odds_hospital_acquired".to_string(), 0.9); // 2.5x higher risk
        map.insert("moraxella_catarrhalis_log_odds_hospital_acquired".to_string(), 0.7); // 2x higher risk
        map.insert("yersinia_enterocolitica_log_odds_hospital_acquired".to_string(), 0.5); // 1.6x higher risk
        
        // Low/No HAI risk bacteria (mostly community pathogens)
        map.insert("chlamydia_trachomatis_log_odds_hospital_acquired".to_string(), -1.0); // 0.37x (lower risk in hospital)
        map.insert("neisseria_gonorrhoeae_log_odds_hospital_acquired".to_string(), -0.8); // 0.45x (lower risk in hospital)
        map.insert("treponema_pallidum_log_odds_hospital_acquired".to_string(), -1.2); // 0.3x (much lower risk)
        map.insert("salmonella_enterica_serovar_typhi_log_odds_hospital_acquired".to_string(), 0.0); // Neutral (travel-related)
        map.insert("salmonella_enterica_serovar_paratyphi_a_log_odds_hospital_acquired".to_string(), 0.0); // Neutral
        map.insert("invasive_non-typhoidal_salmonella_spp._log_odds_hospital_acquired".to_string(), 0.2); // Slight increase
        map.insert("shigella_spp._log_odds_hospital_acquired".to_string(), -0.3); // 0.74x (foodborne, lower in hospital)
        map.insert("vibrio_cholerae_log_odds_hospital_acquired".to_string(), -0.5); // 0.6x (waterborne, lower in hospital)
        map.insert("campylobacter_jejuni_log_odds_hospital_acquired".to_string(), -0.4); // 0.67x (foodborne, lower in hospital)
        
        // effect of region on bacteria acquisition risk (vs north america)
        map.insert("south_america_bacteria_acquisition_log_odds_default".to_string(), 0.4);
        map.insert("africa_shigella_spp_acquisition_log_odds".to_string(), 3.0);
        map.insert("europe_shigella_spp_acquisition_log_odds".to_string(), 0.5); 
        map.insert("asia_shigella_spp_acquisition_log_odds".to_string(), 2.0);
        map.insert("oceania_shigella_spp_acquisition_log_odds".to_string(), 0.7);
 
        map.insert("africa_acinetobacter_baumannii_acquisition_log_odds".to_string(), 1.8);
        map.insert("europe_acinetobacter_baumannii_acquisition_log_odds".to_string(), -0.3);
        map.insert("asia_acinetobacter_baumannii_acquisition_log_odds".to_string(), 1.5);
        map.insert("south_america_acinetobacter_baumannii_acquisition_log_odds".to_string(), 1.2);
        map.insert("oceania_acinetobacter_baumannii_acquisition_log_odds".to_string(), -0.1);
        
        // Citrobacter spp. - Predominantly healthcare-associated, modest regional differences
        map.insert("africa_citrobacter_spp._acquisition_log_odds".to_string(), 0.8);
        map.insert("europe_citrobacter_spp._acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_citrobacter_spp._acquisition_log_odds".to_string(), 0.6);
        map.insert("south_america_citrobacter_spp._acquisition_log_odds".to_string(), 0.4);
        map.insert("oceania_citrobacter_spp._acquisition_log_odds".to_string(), 0.0);
        
        // Enterobacter spp. - Predominantly healthcare-associated, modest regional differences
        map.insert("africa_enterobacter_spp._acquisition_log_odds".to_string(), 0.9);
        map.insert("europe_enterobacter_spp._acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_enterobacter_spp._acquisition_log_odds".to_string(), 0.7);
        map.insert("south_america_enterobacter_spp._acquisition_log_odds".to_string(), 0.5);
        map.insert("oceania_enterobacter_spp._acquisition_log_odds".to_string(), 0.0);
        
        // Enterococcus faecalis - Mixed healthcare/community, moderate regional differences
        map.insert("africa_enterococcus_faecalis_acquisition_log_odds".to_string(), 1.2);
        map.insert("europe_enterococcus_faecalis_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_enterococcus_faecalis_acquisition_log_odds".to_string(), 0.9);
        map.insert("south_america_enterococcus_faecalis_acquisition_log_odds".to_string(), 0.6);
        map.insert("oceania_enterococcus_faecalis_acquisition_log_odds".to_string(), 0.0);
        
        // Enterococcus faecium - Predominantly healthcare-associated, high AMR burden
        map.insert("africa_enterococcus_faecium_acquisition_log_odds".to_string(), 1.0);
        map.insert("europe_enterococcus_faecium_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_enterococcus_faecium_acquisition_log_odds".to_string(), 0.8);
        map.insert("south_america_enterococcus_faecium_acquisition_log_odds".to_string(), 0.5);
        map.insert("oceania_enterococcus_faecium_acquisition_log_odds".to_string(), 0.0);
        
        // Escherichia coli - Major community and healthcare pathogen, high regional variation
        map.insert("africa_escherichia_coli_acquisition_log_odds".to_string(), 1.8);
        map.insert("europe_escherichia_coli_acquisition_log_odds".to_string(), -0.2);
        map.insert("asia_escherichia_coli_acquisition_log_odds".to_string(), 1.5);
        map.insert("south_america_escherichia_coli_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_escherichia_coli_acquisition_log_odds".to_string(), 0.1);
        
        // Klebsiella pneumoniae - Mixed community/healthcare, major AMR threat
        map.insert("africa_klebsiella_pneumoniae_acquisition_log_odds".to_string(), 1.6);
        map.insert("europe_klebsiella_pneumoniae_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_klebsiella_pneumoniae_acquisition_log_odds".to_string(), 1.3);
        map.insert("south_america_klebsiella_pneumoniae_acquisition_log_odds".to_string(), 0.8);
        map.insert("oceania_klebsiella_pneumoniae_acquisition_log_odds".to_string(), 0.0);
        
        // Morganella spp. - Predominantly healthcare-associated, urinary tract infections
        map.insert("africa_morganella_spp._acquisition_log_odds".to_string(), 0.7);
        map.insert("europe_morganella_spp._acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_morganella_spp._acquisition_log_odds".to_string(), 0.5);
        map.insert("south_america_morganella_spp._acquisition_log_odds".to_string(), 0.3);
        map.insert("oceania_morganella_spp._acquisition_log_odds".to_string(), 0.0);
        
        // Proteus spp. - Mixed healthcare/community, urinary tract and wound infections
        map.insert("africa_proteus_spp._acquisition_log_odds".to_string(), 1.1);
        map.insert("europe_proteus_spp._acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_proteus_spp._acquisition_log_odds".to_string(), 0.8);
        map.insert("south_america_proteus_spp._acquisition_log_odds".to_string(), 0.5);
        map.insert("oceania_proteus_spp._acquisition_log_odds".to_string(), 0.0);
        
        // Serratia spp. - Predominantly healthcare-associated, opportunistic pathogen
        map.insert("africa_serratia_spp._acquisition_log_odds".to_string(), 0.8);
        map.insert("europe_serratia_spp._acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_serratia_spp._acquisition_log_odds".to_string(), 0.6);
        map.insert("south_america_serratia_spp._acquisition_log_odds".to_string(), 0.4);
        map.insert("oceania_serratia_spp._acquisition_log_odds".to_string(), 0.0);
        
        // Pseudomonas aeruginosa - Predominantly healthcare-associated, major AMR threat
        map.insert("africa_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 1.2);
        map.insert("europe_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 0.9);
        map.insert("south_america_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 0.6);
        map.insert("oceania_pseudomonas_aeruginosa_acquisition_log_odds".to_string(), 0.0);
        
        // Staphylococcus aureus - Major community and healthcare pathogen, high regional variation
        map.insert("africa_staphylococcus_aureus_acquisition_log_odds".to_string(), 1.5);
        map.insert("europe_staphylococcus_aureus_acquisition_log_odds".to_string(), -0.2);
        map.insert("asia_staphylococcus_aureus_acquisition_log_odds".to_string(), 1.2);
        map.insert("south_america_staphylococcus_aureus_acquisition_log_odds".to_string(), 0.8);
        map.insert("oceania_staphylococcus_aureus_acquisition_log_odds".to_string(), 0.0);
        
        // Streptococcus pneumoniae - Predominantly community-acquired, high regional variation
        map.insert("africa_streptococcus_pneumoniae_acquisition_log_odds".to_string(), 2.2);
        map.insert("europe_streptococcus_pneumoniae_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_streptococcus_pneumoniae_acquisition_log_odds".to_string(), 1.8);
        map.insert("south_america_streptococcus_pneumoniae_acquisition_log_odds".to_string(), 1.2);
        map.insert("oceania_streptococcus_pneumoniae_acquisition_log_odds".to_string(), 0.3);
        
        // Salmonella enterica serovar typhi - Typhoid fever, highly endemic in certain regions
        map.insert("africa_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), 2.8);
        map.insert("europe_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), -2.5);
        map.insert("asia_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), 2.5);
        map.insert("south_america_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_salmonella_enterica_serovar_typhi_acquisition_log_odds".to_string(), -1.8);
        
        // Salmonella enterica serovar paratyphi a - Paratyphoid fever, similar but less common
        map.insert("africa_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), 2.5);
        map.insert("europe_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), -2.8);
        map.insert("asia_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), 2.2);
        map.insert("south_america_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), 0.8);
        map.insert("oceania_salmonella_enterica_serovar_paratyphi_a_acquisition_log_odds".to_string(), -2.0);
        
        // Invasive non-typhoidal salmonella spp. - Bloodstream infections, especially in immunocompromised
        map.insert("africa_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), 3.2);
        map.insert("europe_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), -1.0);
        map.insert("asia_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), 1.8);
        map.insert("south_america_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), 1.2);
        map.insert("oceania_invasive_non-typhoidal_salmonella_spp._acquisition_log_odds".to_string(), -0.5);
        
        // Neisseria gonorrhoeae - Sexually transmitted infection, moderate regional variation
        map.insert("africa_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 1.6);
        map.insert("europe_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 0.1);
        map.insert("asia_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 1.2);
        map.insert("south_america_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 0.8);
        map.insert("oceania_neisseria_gonorrhoeae_acquisition_log_odds".to_string(), 0.2);
        
        // Streptococcus pyogenes - Group A Strep, community-acquired, moderate regional variation
        map.insert("africa_streptococcus_pyogenes_acquisition_log_odds".to_string(), 1.8);
        map.insert("europe_streptococcus_pyogenes_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_streptococcus_pyogenes_acquisition_log_odds".to_string(), 1.4);
        map.insert("south_america_streptococcus_pyogenes_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_streptococcus_pyogenes_acquisition_log_odds".to_string(), 0.2);
        
        // Streptococcus agalactiae - Group B Strep, neonatal/maternal infections, moderate variation
        map.insert("africa_streptococcus_agalactiae_acquisition_log_odds".to_string(), 1.4);
        map.insert("europe_streptococcus_agalactiae_acquisition_log_odds".to_string(), 0.0);
        map.insert("asia_streptococcus_agalactiae_acquisition_log_odds".to_string(), 1.1);
        map.insert("south_america_streptococcus_agalactiae_acquisition_log_odds".to_string(), 0.7);
        map.insert("oceania_streptococcus_agalactiae_acquisition_log_odds".to_string(), 0.1);
        
        // Haemophilus influenzae - Respiratory pathogen, dramatically reduced by Hib vaccine
        map.insert("africa_haemophilus_influenzae_acquisition_log_odds".to_string(), 2.0);
        map.insert("europe_haemophilus_influenzae_acquisition_log_odds".to_string(), -0.3);
        map.insert("asia_haemophilus_influenzae_acquisition_log_odds".to_string(), 1.6);
        map.insert("south_america_haemophilus_influenzae_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_haemophilus_influenzae_acquisition_log_odds".to_string(), -0.2);
        
        // Chlamydia trachomatis - STI and trachoma, moderate regional variation
        map.insert("africa_chlamydia_trachomatis_acquisition_log_odds".to_string(), 1.4);
        map.insert("europe_chlamydia_trachomatis_acquisition_log_odds".to_string(), 0.0);
        map.insert("asia_chlamydia_trachomatis_acquisition_log_odds".to_string(), 1.0);
        map.insert("south_america_chlamydia_trachomatis_acquisition_log_odds".to_string(), 0.6);
        map.insert("oceania_chlamydia_trachomatis_acquisition_log_odds".to_string(), 0.1);
        
        // Vibrio cholerae - Waterborne disease, extreme regional variation
        map.insert("africa_vibrio_cholerae_acquisition_log_odds".to_string(), 3.5);
        map.insert("europe_vibrio_cholerae_acquisition_log_odds".to_string(), -3.0);
        map.insert("asia_vibrio_cholerae_acquisition_log_odds".to_string(), 2.8);
        map.insert("south_america_vibrio_cholerae_acquisition_log_odds".to_string(), 1.5);
        map.insert("oceania_vibrio_cholerae_acquisition_log_odds".to_string(), -2.2);
        
        // Neisseria meningitidis - Meningococcal disease, moderate regional variation with vaccine impact
        map.insert("africa_neisseria_meningitidis_acquisition_log_odds".to_string(), 2.5);
        map.insert("europe_neisseria_meningitidis_acquisition_log_odds".to_string(), -0.2);
        map.insert("asia_neisseria_meningitidis_acquisition_log_odds".to_string(), 1.2);
        map.insert("south_america_neisseria_meningitidis_acquisition_log_odds".to_string(), 0.8);
        map.insert("oceania_neisseria_meningitidis_acquisition_log_odds".to_string(), -0.1);
        
        // Listeria monocytogenes - Foodborne pathogen, moderate regional variation
        map.insert("africa_listeria_monocytogenes_acquisition_log_odds".to_string(), 1.0);
        map.insert("europe_listeria_monocytogenes_acquisition_log_odds".to_string(), 0.1);
        map.insert("asia_listeria_monocytogenes_acquisition_log_odds".to_string(), 0.8);
        map.insert("south_america_listeria_monocytogenes_acquisition_log_odds".to_string(), 0.5);
        map.insert("oceania_listeria_monocytogenes_acquisition_log_odds".to_string(), 0.0);
        
        // Clostridioides difficile - Healthcare-associated, antibiotic-driven, modest regional variation
        map.insert("africa_clostridioides_difficile_acquisition_log_odds".to_string(), 0.6);
        map.insert("europe_clostridioides_difficile_acquisition_log_odds".to_string(), 0.1);
        map.insert("asia_clostridioides_difficile_acquisition_log_odds".to_string(), 0.4);
        map.insert("south_america_clostridioides_difficile_acquisition_log_odds".to_string(), 0.3);
        map.insert("oceania_clostridioides_difficile_acquisition_log_odds".to_string(), 0.0);
        
        // Campylobacter jejuni - Foodborne pathogen, moderate regional variation
        map.insert("africa_campylobacter_jejuni_acquisition_log_odds".to_string(), 1.3);
        map.insert("europe_campylobacter_jejuni_acquisition_log_odds".to_string(), 0.2);
        map.insert("asia_campylobacter_jejuni_acquisition_log_odds".to_string(), 1.0);
        map.insert("south_america_campylobacter_jejuni_acquisition_log_odds".to_string(), 0.7);
        map.insert("oceania_campylobacter_jejuni_acquisition_log_odds".to_string(), 0.1);
        
        // Enterobacter cloacae - Healthcare-associated Enterobacteriaceae, modest regional variation
        map.insert("africa_enterobacter_cloacae_acquisition_log_odds".to_string(), 0.9);
        map.insert("europe_enterobacter_cloacae_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_enterobacter_cloacae_acquisition_log_odds".to_string(), 0.7);
        map.insert("south_america_enterobacter_cloacae_acquisition_log_odds".to_string(), 0.5);
        map.insert("oceania_enterobacter_cloacae_acquisition_log_odds".to_string(), 0.0);
        
        // Yersinia enterocolitica - Foodborne/zoonotic pathogen, moderate regional variation
        map.insert("africa_yersinia_enterocolitica_acquisition_log_odds".to_string(), 1.1);
        map.insert("europe_yersinia_enterocolitica_acquisition_log_odds".to_string(), 0.3);
        map.insert("asia_yersinia_enterocolitica_acquisition_log_odds".to_string(), 0.8);
        map.insert("south_america_yersinia_enterocolitica_acquisition_log_odds".to_string(), 0.5);
        map.insert("oceania_yersinia_enterocolitica_acquisition_log_odds".to_string(), 0.1);
        
        // Moraxella catarrhalis - Respiratory pathogen, moderate regional variation
        map.insert("africa_moraxella_catarrhalis_acquisition_log_odds".to_string(), 1.6);
        map.insert("europe_moraxella_catarrhalis_acquisition_log_odds".to_string(), -0.1);
        map.insert("asia_moraxella_catarrhalis_acquisition_log_odds".to_string(), 1.2);
        map.insert("south_america_moraxella_catarrhalis_acquisition_log_odds".to_string(), 0.8);
        map.insert("oceania_moraxella_catarrhalis_acquisition_log_odds".to_string(), 0.0);
        
        // Treponema pallidum - Syphilis, moderate-high regional variation
        map.insert("africa_treponema_pallidum_acquisition_log_odds".to_string(), 2.0);
        map.insert("europe_treponema_pallidum_acquisition_log_odds".to_string(), 0.0);
        map.insert("asia_treponema_pallidum_acquisition_log_odds".to_string(), 1.5);
        map.insert("south_america_treponema_pallidum_acquisition_log_odds".to_string(), 1.0);
        map.insert("oceania_treponema_pallidum_acquisition_log_odds".to_string(), 0.2);
        
        // Bacteria-specific microbiome vs infection acquisition log odds
        // High carriage bacteria (common gut/skin commensals)
        map.insert("escherichia_coli_log_odds_microbiome_vs_infection".to_string(), 2.5); // Very high carriage rate
        map.insert("enterococcus_faecalis_log_odds_microbiome_vs_infection".to_string(), 2.2); // High gut carriage
        map.insert("enterococcus_faecium_log_odds_microbiome_vs_infection".to_string(), 2.0); // High gut carriage
        map.insert("klebsiella_pneumoniae_log_odds_microbiome_vs_infection".to_string(), 1.8); // Moderate-high gut carriage
        map.insert("staphylococcus_aureus_log_odds_microbiome_vs_infection".to_string(), 1.5); // ~30% nasal carriage
        
        // Moderate carriage bacteria (opportunistic commensals)
        map.insert("enterobacter_spp._log_odds_microbiome_vs_infection".to_string(), 1.3);
        map.insert("enterobacter_cloacae_log_odds_microbiome_vs_infection".to_string(), 1.2);
        map.insert("citrobacter_spp._log_odds_microbiome_vs_infection".to_string(), 1.1);
        map.insert("proteus_spp._log_odds_microbiome_vs_infection".to_string(), 1.0);
        map.insert("serratia_spp._log_odds_microbiome_vs_infection".to_string(), 0.8);
        map.insert("morganella_spp._log_odds_microbiome_vs_infection".to_string(), 0.7);
        
        // Respiratory tract commensals (episodic carriage)
        map.insert("streptococcus_pneumoniae_log_odds_microbiome_vs_infection".to_string(), 1.2); // Nasopharyngeal carriage
        map.insert("haemophilus_influenzae_log_odds_microbiome_vs_infection".to_string(), 1.0); // Respiratory carriage
        map.insert("moraxella_catarrhalis_log_odds_microbiome_vs_infection".to_string(), 0.9); // Upper respiratory carriage
        map.insert("streptococcus_pyogenes_log_odds_microbiome_vs_infection".to_string(), 0.5); // Transient throat carriage
        map.insert("streptococcus_agalactiae_log_odds_microbiome_vs_infection".to_string(), 0.8); // GI/genital carriage
        
        // Healthcare-associated, low community carriage
        map.insert("acinetobacter_baumannii_log_odds_microbiome_vs_infection".to_string(), 0.3); // Mainly hospital environment
        map.insert("pseudomonas_aeruginosa_log_odds_microbiome_vs_infection".to_string(), 0.2); // Low carriage, environmental
        map.insert("clostridioides_difficile_log_odds_microbiome_vs_infection".to_string(), 0.8); // Spore-forming, gut carriage
        
        // Foodborne/environmental, minimal carriage
        map.insert("salmonella_enterica_serovar_typhi_log_odds_microbiome_vs_infection".to_string(), -1.0); // Chronic carriage rare
        map.insert("salmonella_enterica_serovar_paratyphi_a_log_odds_microbiome_vs_infection".to_string(), -1.2); // Minimal carriage
        map.insert("invasive_non-typhoidal_salmonella_spp._log_odds_microbiome_vs_infection".to_string(), -0.5); // Some gut carriage
        map.insert("shigella_spp._log_odds_microbiome_vs_infection".to_string(), -0.8); // Minimal carriage
        map.insert("vibrio_cholerae_log_odds_microbiome_vs_infection".to_string(), -2.0); // Almost no carriage
        map.insert("campylobacter_jejuni_log_odds_microbiome_vs_infection".to_string(), -1.0); // Minimal human carriage
        map.insert("yersinia_enterocolitica_log_odds_microbiome_vs_infection".to_string(), -0.7); // Low carriage
        map.insert("listeria_monocytogenes_log_odds_microbiome_vs_infection".to_string(), -1.5); // Rare carriage
        
        // Sexually transmitted, no meaningful carriage
        map.insert("neisseria_gonorrhoeae_log_odds_microbiome_vs_infection".to_string(), -2.5); // No carriage
        map.insert("chlamydia_trachomatis_log_odds_microbiome_vs_infection".to_string(), -2.0); // Intracellular, no carriage
        map.insert("treponema_pallidum_log_odds_microbiome_vs_infection".to_string(), -3.0); // No carriage
        
        // Other specialized pathogens
        map.insert("neisseria_meningitidis_log_odds_microbiome_vs_infection".to_string(), 0.5); // Nasopharyngeal carriage
        
        // Microbiome acquisition now uses infection acquisition parameters plus bacteria-specific offset
        // Fallback parameter for backward compatibility
        map.insert("log_odds_microbiome_vs_infection".to_string(), 1.0); // Fallback if bacteria-specific parameter not found

        // Environmental resistance level for new acquisitions
        map.insert("environmental_majority_r_level_for_new_acquisition".to_string(), 0.001); // 0.01 

  
        map.insert("max_resistance_level".to_string(), 1.0);
        map.insert("majority_r_evolution_rate_per_day_when_drug_present".to_string(), 0.01); // 0.01

        // Resistance Emergence and Decay Parameters
        // Resistance reversion parameter: probability per day that resistance reverts to 0 if no drug present
        map.insert("resistance_reversion_rate_per_day".to_string(), 0.0001); // Default: very rare, increase for more rapid reversion
        map.insert("microbiome_resistance_emergence_rate_per_day_baseline".to_string(), 0.005); // Separate baseline for microbiome resistance emergence
        map.insert("resistance_emergence_bacteria_level_multiplier".to_string(), 0.05); // Multiplier for bacteria level's effect on emergence
        map.insert("any_r_emergence_level_on_first_emergence".to_string(), 0.5); // The resistance level 'any_r' starts at upon emergence

        
        //  Microbiome Resistance Transfer Parameter
        map.insert("microbiome_resistance_transfer_probability_per_day".to_string(), 0.05); // Probability per day for resistance transfer between infection and microbiome
    
        // --- Multi-Drug Resistance Emergence Penalty Parameters ---
        // When multiple drugs are active, resistance emergence is reduced because mutations
        // must confer resistance to ALL active drugs to provide survival advantage
        map.insert("multi_drug_penalty_for_single_drug_resistance".to_string(), 0.05); // Penalty when resistance affects only 1 of multiple active drugs (5% survival)
        map.insert("multi_drug_penalty_for_partial_cross_resistance".to_string(), 0.3); // Penalty when resistance affects some but not all active drugs (30% survival)
        map.insert("multi_drug_penalty_threshold_num_drugs".to_string(), 2.0); // Minimum number of active drugs to trigger multi-drug penalty

        // --- Resistance Mechanisms Parameters ---
        // Baseline emergence rates for specific resistance mechanisms (per day when drug present)
        map.insert("resistance_mechanism_esbl_emergence_rate".to_string(), 0.001); // ESBL emergence with beta-lactam pressure
        map.insert("resistance_mechanism_carbapenemase_emergence_rate".to_string(), 0.0005); // Carbapenemase emergence (rarer)
        map.insert("resistance_mechanism_ampc_emergence_rate".to_string(), 0.002); // AmpC emergence with beta-lactam pressure
        map.insert("resistance_mechanism_16s_methyltransferase_emergence_rate".to_string(), 0.001); // Aminoglycoside resistance
        map.insert("resistance_mechanism_qnr_emergence_rate".to_string(), 0.001); // Quinolone resistance
        map.insert("resistance_mechanism_efflux_overexpression_emergence_rate".to_string(), 0.003); // More common mechanism
        map.insert("resistance_mechanism_erm_methylation_emergence_rate".to_string(), 0.001); // Macrolide resistance
        map.insert("resistance_mechanism_van_type_emergence_rate".to_string(), 0.0002); // Vancomycin resistance (rare)
        map.insert("resistance_mechanism_meca_emergence_rate".to_string(), 0.0008); // MRSA emergence
        map.insert("resistance_mechanism_reduced_permeability_emergence_rate".to_string(), 0.002); // Common adaptive mechanism
        map.insert("resistance_mechanism_target_site_mutation_emergence_rate".to_string(), 0.0015); // Point mutations

        // Resistance enhancement multipliers: how much each mechanism increases resistance level
        map.insert("resistance_mechanism_esbl_enhancement_multiplier".to_string(), 0.4); // Adds 40% resistance
        map.insert("resistance_mechanism_carbapenemase_enhancement_multiplier".to_string(), 0.6); // Adds 60% resistance  
        map.insert("resistance_mechanism_ampc_enhancement_multiplier".to_string(), 0.3); // Adds 30% resistance
        map.insert("resistance_mechanism_16s_methyltransferase_enhancement_multiplier".to_string(), 0.5); // Adds 50% resistance
        map.insert("resistance_mechanism_qnr_enhancement_multiplier".to_string(), 0.2); // Adds 20% resistance (low-level)
        map.insert("resistance_mechanism_efflux_overexpression_enhancement_multiplier".to_string(), 0.3); // Adds 30% resistance
        map.insert("resistance_mechanism_erm_methylation_enhancement_multiplier".to_string(), 0.5); // Adds 50% resistance
        map.insert("resistance_mechanism_van_type_enhancement_multiplier".to_string(), 0.8); // Adds 80% resistance (high-level)
        map.insert("resistance_mechanism_meca_enhancement_multiplier".to_string(), 0.7); // Adds 70% resistance
        map.insert("resistance_mechanism_reduced_permeability_enhancement_multiplier".to_string(), 0.2); // Adds 20% resistance
        map.insert("resistance_mechanism_target_site_mutation_enhancement_multiplier".to_string(), 0.4); // Adds 40% resistance
    
        map.insert("mechanism_assignment_probability_on_any_r_gain".to_string(), 0.8); // Default 80%

        // Testing Parameters
        map.insert("bacterial_testing_available_from_day".to_string(), 5478.0); // 5478.0  1945 (15 years after 1930) - Bacterial culture/identification becomes available
        map.insert("resistance_testing_available_from_day".to_string(), 9131.0); // 9131.0  1955 (25 years after 1930) - Antibiotic susceptibility testing becomes available
        map.insert("test_delay_days".to_string(), 3.0);
        map.insert("test_rate_per_day".to_string(), 0.2);  // 0.15

        // --- Test result and test_r logic parameters ---
        map.insert("prob_test_r_done".to_string(), 0.95); // Probability test is actually done (per day eligible)
        map.insert("test_r_error_probability".to_string(), 0.02); // Probability of error in test result
        map.insert("test_r_error_value".to_string(), 0.25); // Value to use for error in test_r

        // Syndrome-specific multipliers (example)
        map.insert("syndrome_3_initiation_multiplier".to_string(), 10.0); // Respiratory syndrome
        map.insert("syndrome_7_initiation_multiplier".to_string(), 8.0);  // Gastrointestinal syndrome
        map.insert("syndrome_8_initiation_multiplier".to_string(), 12.0); // Genital syndrome (example ID)        

        // Hospitalization Parameters
        map.insert("hospitalization_baseline_rate_per_day".to_string(), 0.0000005); // 0.0000001  Baseline daily probability of hospitalization
        map.insert("hospitalization_age_multiplier_per_day".to_string(), 0.0000005); // 0.0000001  Increase in daily hospitalization probability per year of age
        map.insert("hospitalization_recovery_rate_per_day".to_string(), 0.2); // 0.2  Daily probability of recovering from hospitalization
        map.insert("hospitalization_max_days".to_string(), 30.0); // 30.0  Max days in hospital before forced discharge (as fallback)
        map.insert("hospitalization_sepsis_admission_multiplier".to_string(), 1000.0); // Strong predictor: sepsis patients very likely to be hospitalized
        map.insert("hospitalization_prevent_discharge_with_sepsis".to_string(), 1.0); // 1.0 = true, 0.0 = false: prevent discharge of patients with active sepsis

        // Testing Framework Parameters
        // Base testing rates (modern era baseline)
        map.insert("bacterial_testing_base_rate_per_day".to_string(), 0.15); // Modern baseline rate for bacterial identification
        map.insert("resistance_testing_base_rate_per_day".to_string(), 0.95); // Probability of resistance testing given bacterial identification
        
        // Hospital status multipliers
        map.insert("bacterial_testing_hospital_multiplier".to_string(), 8.0); // Hospitalized patients 8x more likely to get bacterial testing
        map.insert("resistance_testing_hospital_multiplier".to_string(), 5.0); // Hospitalized patients 5x more likely to get resistance testing
        
        // Regional resource multipliers for testing access
        map.insert("north_america_testing_multiplier".to_string(), 1.1);
        map.insert("europe_testing_multiplier".to_string(), 1.2);
        map.insert("asia_testing_multiplier".to_string(), 0.7);
        map.insert("south_america_testing_multiplier".to_string(), 0.6);
        map.insert("oceania_testing_multiplier".to_string(), 0.8);
        map.insert("africa_testing_multiplier".to_string(), 0.3); // Limited lab infrastructure
        
        // Clinical status multipliers
        map.insert("testing_immunosuppressed_multiplier".to_string(), 2.5); // Immunosuppressed patients get more testing
        map.insert("testing_sepsis_multiplier".to_string(), 4.0); // Sepsis patients get urgent testing
        
        // Temporal adoption parameters (testing evolution over time)
        // Bacterial testing temporal evolution
        map.insert("bacterial_testing_initial_adoption_rate".to_string(), 0.1); // 1945: 10% of modern rates
        map.insert("bacterial_testing_adoption_rate_per_year".to_string(), 0.025); // 2.5% improvement per year
        map.insert("bacterial_testing_max_temporal_multiplier".to_string(), 1.0); // Cap at modern rates
        
        // Resistance testing temporal evolution (slower adoption)
        map.insert("resistance_testing_initial_adoption_rate".to_string(), 0.05); // 1955: 5% of modern rates
        map.insert("resistance_testing_adoption_rate_per_year".to_string(), 0.015); // 1.5% improvement per year
        map.insert("resistance_testing_max_temporal_multiplier".to_string(), 1.0); // Cap at modern rates

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
        map.insert("empiric_therapy_ineffective_drug_penalty".to_string(), 0.05); // Heavy penalty for drugs ineffective against actual pathogens (empirical)
        map.insert("targeted_therapy_narrow_spectrum_bonus".to_string(), 3.0); // Multiplier for narrow-spectrum drugs when bacteria identified  
        map.insert("targeted_therapy_broad_spectrum_penalty".to_string(), 0.4); // Penalty for broad-spectrum drugs when bacteria identified
        map.insert("targeted_therapy_ineffective_drug_penalty".to_string(), 0.1); // Strong penalty for drugs ineffective against identified bacteria

        // Drug Spectrum Classifications (1.0=narrow, 5.0=very broad)
        map.insert("drug_colistin_spectrum_breadth".to_string(), 4.0); // Broad spectrum (mainly Gram-negative)
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


        // NEW: Logistic Sepsis Risk Parameters (replacing old linear model)
        map.insert("sepsis_baseline_odds".to_string(), -10.0); // -10.0 Baseline log odds (very low baseline probability)
        map.insert("log_odds_sepsis_infection_level".to_string(), 2.0); // Log odds increase per unit bacterial level
        map.insert("log_odds_sepsis_infection_duration".to_string(), 0.001); // Log odds increase per day of infection duration
        map.insert("log_odds_bacteria_with_high_sepsis_risk".to_string(), 1.0); // Log odds for high-risk bacteria (e.g., exp(1.0) = 2.7x odds ratio)
        map.insert("log_odds_bacteria_with_medium_sepsis_risk".to_string(), 0.0); // Log odds for medium-risk bacteria (reference category)
        map.insert("log_odds_bacteria_with_low_sepsis_risk".to_string(), -1.2); // Log odds for low-risk bacteria (e.g., exp(-1.2) = 0.3x odds ratio)

        // Syndrome-specific sepsis risk parameters (infectious site effects)
        map.insert("log_odds_syndrome_1_sepsis".to_string(), -2.0); // UTI/Genitourinary: Much lower sepsis risk
        map.insert("log_odds_syndrome_2_sepsis".to_string(), -1.0); // Skin/Soft tissue: Lower sepsis risk
        map.insert("log_odds_syndrome_3_sepsis".to_string(), 0.0);  // Respiratory: Reference category
        map.insert("log_odds_syndrome_4_sepsis".to_string(), 1.5);  // Bloodstream/Bacteremia: Much higher sepsis risk
        map.insert("log_odds_syndrome_5_sepsis".to_string(), 0.8);  // Intra-abdominal: Higher sepsis risk
        map.insert("log_odds_syndrome_6_sepsis".to_string(), 1.2);  // Central nervous system: High sepsis risk
        map.insert("log_odds_syndrome_7_sepsis".to_string(), -0.5); // Gastrointestinal: Somewhat lower sepsis risk
        map.insert("log_odds_syndrome_8_sepsis".to_string(), -1.5); // Genital: Lower sepsis risk
        map.insert("log_odds_syndrome_9_sepsis".to_string(), 0.5);  // Bone/Joint: Moderately higher sepsis risk

        // // Background Mortality Parameters (Age, Region, and Sex dependent)

        // These parameters are on the log-odds scale.
        map.insert("background_mortality_baseline_log_odds".to_string(), -14.0); // -15.5
        map.insert("log_odds_mortality_per_year_of_age".to_string(), 0.04); // 0.04  Odds of dying increase by ~4% per year (exp(0.04) ≈ 1.04)
        map.insert("log_odds_mortality_per_year_of_age_squared".to_string(), 0.05); // 0.05  Additional non-linear effect for elderly

        // Time-varying background mortality (1930-2035): reflects dramatic mortality improvements over 105 years
        // Based on historical life expectancy trends showing ~3x mortality decline from 1930 to 2035
        map.insert("mortality_baseline_1930_multiplier".to_string(), 3.0);  // 3x higher mortality in 1930 (life expectancy ~35-45 years)
        map.insert("mortality_baseline_2035_multiplier".to_string(), 1.0);  // Modern baseline by 2035 (life expectancy ~70-80 years)
        map.insert("mortality_improvement_half_life_years".to_string(), 35.0); // Half-life of mortality improvement (exponential decay)

        // Region-specific log-odds adjustments. ln(1.0) = 0.
        map.insert("log_odds_mortality_region_north_america".to_string(), 0.0);      // Reference
        map.insert("log_odds_mortality_region_south_america".to_string(), 0.26);     // ln(1.3) - increased from 1.2
        map.insert("log_odds_mortality_region_africa".to_string(), 0.41);            // ln(1.5) - increased from 1.2 to reflect TB, malaria, parasitic diseases
        map.insert("log_odds_mortality_region_asia".to_string(), 0.18);              // ln(1.2) - increased from 1.1 
        map.insert("log_odds_mortality_region_europe".to_string(), -0.105);          // ln(0.9) - kept same
        map.insert("log_odds_mortality_region_oceania".to_string(), 0.0);           // Reference

        // Sex-specific log-odds adjustments.
        map.insert("log_odds_mortality_sex_male".to_string(), 0.095);    // ln(1.1)
        map.insert("log_odds_mortality_sex_female".to_string(), -0.105); // ln(0.9)

        // Additional background mortality risk factors (log-odds).
        map.insert("log_odds_mortality_immunosuppressed".to_string(), 0.916); // ln(2.5)
        map.insert("log_odds_mortality_hospitalized".to_string(), 0.262);     // ln(1.3)


        //  Immunosuppression Onset and Recovery Rates
        // Temporary immunodeficiency (e.g., chemotherapy, short-term steroids)
        map.insert("temporary_immunosuppression_onset_rate_per_day".to_string(), 0.0002);   // Slightly lower onset rate for temporary
        map.insert("temporary_immunosuppression_recovery_rate_per_day".to_string(), 0.01);   // Faster recovery (10x faster)
        
        // Chronic immunodeficiency (e.g., HIV, genetic disorders, organ transplant)
        map.insert("chronic_immunosuppression_onset_rate_per_day".to_string(), 0.0001);     // Lower onset rate for chronic
        map.insert("chronic_immunosuppression_recovery_rate_per_day".to_string(), 0.0005);  // Much slower recovery (10x slower)
        
        // Age effect on immunodeficiency type assignment (probability of chronic vs temporary at onset)
        map.insert("chronic_immunodeficiency_probability_age_0_1".to_string(), 0.3);   // Infants: higher chance of genetic/congenital
        map.insert("chronic_immunodeficiency_probability_age_1_18".to_string(), 0.2);  // Children: moderate chance
        map.insert("chronic_immunodeficiency_probability_age_18_65".to_string(), 0.4); // Adults: higher chance (HIV, transplants)
        map.insert("chronic_immunodeficiency_probability_age_65_plus".to_string(), 0.6); // Elderly: highest chance (multiple conditions)

        // Prophylactic antibiotic use in immunocompromised patients
        map.insert("immunodeficiency_prophylactic_drug_multiplier".to_string(), 8.0);  // 8x higher drug initiation rate for immunocompromised (prophylaxis)
        map.insert("antibiotic_infection_prevention_efficacy".to_string(), 0.85);      // 85% efficacy: existing antibiotic prevents new susceptible infections


        // Sepsis Mortality Parameters (Age, Region, and Risk Factor dependent)
        map.insert("base_sepsis_death_risk_per_day".to_string(), 0.02); // 0.02 Base 2% daily death risk for sepsis 
        map.insert("sepsis_age_mortality_multiplier_infant".to_string(), 3.0); // 0-1 years: much higher risk
        map.insert("sepsis_age_mortality_multiplier_child".to_string(), 0.5); // 1-18 years: lower risk  
        map.insert("sepsis_age_mortality_multiplier_adult".to_string(), 1.0); // 18-65 years: baseline risk
        map.insert("sepsis_age_mortality_multiplier_elderly".to_string(), 2.5); // 65+ years: much higher risk
        map.insert("sepsis_immunosuppressed_multiplier".to_string(), 30.0); // Immunosuppressed: 30x higher risk
        
        // Region-specific sepsis mortality multipliers (reflecting healthcare quality)
        map.insert("north_america_sepsis_mortality_multiplier".to_string(), 0.8); // Better ICU care
        map.insert("europe_sepsis_mortality_multiplier".to_string(), 0.7); // Excellent healthcare systems
        map.insert("oceania_sepsis_mortality_multiplier".to_string(), 0.8); // Good healthcare
        map.insert("asia_sepsis_mortality_multiplier".to_string(), 1.2); // Variable healthcare quality
        map.insert("south_america_sepsis_mortality_multiplier".to_string(), 1.4); // Limited ICU access
        map.insert("africa_sepsis_mortality_multiplier".to_string(), 2.0); // Limited healthcare infrastructure

        // Sepsis Recovery Parameters (Logistic Model)
        map.insert("sepsis_base_log_odds_of_recovery_per_day".to_string(), -1.5); // Base log odds (low baseline recovery probability ~12%)
        map.insert("sepsis_log_odds_bacteria_level".to_string(), -0.3); // Higher bacteria level decreases recovery (negative coefficient)
        map.insert("sepsis_log_odds_in_hospital".to_string(), 0.8); // Being in hospital increases recovery probability
        map.insert("sepsis_log_odds_age_infant".to_string(), -0.5); // Infants have lower recovery probability
        map.insert("sepsis_log_odds_age_child".to_string(), 0.4); // Children have higher recovery probability  
        map.insert("sepsis_log_odds_age_adult".to_string(), 0.0); // Adults baseline (reference category)
        map.insert("sepsis_log_odds_age_elderly".to_string(), -0.7); // Elderly have lower recovery probability
        map.insert("sepsis_log_odds_immunosuppressed".to_string(), -1.0); // Immunosuppressed have much lower recovery probability
        
        // Region-specific sepsis recovery log odds (reflecting healthcare quality and ICU availability)
        map.insert("sepsis_log_odds_region_north_america".to_string(), 0.4); // Better healthcare systems increase recovery
        map.insert("sepsis_log_odds_region_europe".to_string(), 0.5); // Excellent healthcare systems, best recovery rates
        map.insert("sepsis_log_odds_region_oceania".to_string(), 0.3); // Good healthcare systems
        map.insert("sepsis_log_odds_region_asia".to_string(), 0.0); // Mixed healthcare quality, reference category
        map.insert("sepsis_log_odds_region_south_america".to_string(), -0.3); // Limited ICU access decreases recovery
        map.insert("sepsis_log_odds_region_africa".to_string(), -0.7); // Limited healthcare infrastructure significantly decreases recovery
        map.insert("sepsis_log_odds_region_home".to_string(), 0.0); // Default to reference category
        
        map.insert("sepsis_minimum_duration_days".to_string(), 1.0); // Minimum sepsis duration (1 day)

        //  Default Toxicity Parameter
        //  Default Microbiome Clearance Parameter (required by simulation logic)
        map.insert("default_microbiome_clearance_probability_per_day".to_string(), 0.01); // E.g., 1% chance to lose carriage per day
        // Probability of clearing microbiome when drug treatment successfully clears infection  
        map.insert("microbiome_clearance_probability_on_drug_treatment".to_string(), 0.8); // 80% chance drugs clear microbiome when they clear infection
        map.insert("default_drug_toxicity_per_unit_level_per_day".to_string(), 0.005); // Adjust this default as needed

        //  Probability per day of death due to adverse drug effect (toxicity)
        //  This is the baseline daily risk of death for an individual experiencing drug toxicity.
        //  You can tune this value to make drug toxicity more or less lethal.
        map.insert("drug_toxicity_death_risk_per_day".to_string(), 0.0003); 



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
        
        // Demographic distribution parameters (108 total: 6 regions × 18 age bands)
        // Each parameter represents probability of being in that region-age combination
        // All 108 parameters should sum to 1.0
        // Age bands span from -40000 to +32000 in 4000-year intervals
        
        // Asia demographic distribution (18 age bands of 4000 years each)
        // 55% of global population based on 1930-2035 timeframe
        // Massive growth: ~90% at negative ages (future births), only ~10% alive in 1930
        map.insert("demo_asia_age_neg40000_neg36000".to_string(), 0.120); // Heavy weighting on future births
        map.insert("demo_asia_age_neg36000_neg32000".to_string(), 0.110);
        map.insert("demo_asia_age_neg32000_neg28000".to_string(), 0.100);
        map.insert("demo_asia_age_neg28000_neg24000".to_string(), 0.090);
        map.insert("demo_asia_age_neg24000_neg20000".to_string(), 0.080);
        map.insert("demo_asia_age_neg20000_neg16000".to_string(), 0.070);
        map.insert("demo_asia_age_neg16000_neg12000".to_string(), 0.060);
        map.insert("demo_asia_age_neg12000_neg8000".to_string(), 0.050);
        map.insert("demo_asia_age_neg8000_neg4000".to_string(), 0.040);
        map.insert("demo_asia_age_neg4000_0".to_string(), 0.030);          // Future births tapering
        map.insert("demo_asia_age_0_4000".to_string(), 0.008);             // Very small portion alive in 1930
        map.insert("demo_asia_age_4000_8000".to_string(), 0.006);
        map.insert("demo_asia_age_8000_12000".to_string(), 0.005);
        map.insert("demo_asia_age_12000_16000".to_string(), 0.004);
        map.insert("demo_asia_age_16000_20000".to_string(), 0.003);
        map.insert("demo_asia_age_20000_24000".to_string(), 0.002);
        map.insert("demo_asia_age_24000_28000".to_string(), 0.001);
        map.insert("demo_asia_age_28000_32000".to_string(), 0.001);        
        
        // Africa demographic distribution
        // 12% of global population based on 1930-2035 timeframe
        // Explosive growth: ~90% at negative ages (future births), only ~10% alive in 1930
        map.insert("demo_africa_age_neg40000_neg36000".to_string(), 0.026); // Heavy weighting on future births
        map.insert("demo_africa_age_neg36000_neg32000".to_string(), 0.024);
        map.insert("demo_africa_age_neg32000_neg28000".to_string(), 0.022);
        map.insert("demo_africa_age_neg28000_neg24000".to_string(), 0.020);
        map.insert("demo_africa_age_neg24000_neg20000".to_string(), 0.018);
        map.insert("demo_africa_age_neg20000_neg16000".to_string(), 0.016);
        map.insert("demo_africa_age_neg16000_neg12000".to_string(), 0.014);
        map.insert("demo_africa_age_neg12000_neg8000".to_string(), 0.012);
        map.insert("demo_africa_age_neg8000_neg4000".to_string(), 0.010);
        map.insert("demo_africa_age_neg4000_0".to_string(), 0.008);          // Future births tapering
        map.insert("demo_africa_age_0_4000".to_string(), 0.002);             // Very small portion alive in 1930
        map.insert("demo_africa_age_4000_8000".to_string(), 0.002);
        map.insert("demo_africa_age_8000_12000".to_string(), 0.001);
        map.insert("demo_africa_age_12000_16000".to_string(), 0.001);
        map.insert("demo_africa_age_16000_20000".to_string(), 0.001);
        map.insert("demo_africa_age_20000_24000".to_string(), 0.001);
        map.insert("demo_africa_age_24000_28000".to_string(), 0.001);
        map.insert("demo_africa_age_28000_32000".to_string(), 0.001);
        
        // Europe demographic distribution
        // 15% of global population based on 1930-2035 timeframe
        // Moderate growth: ~70% at negative ages (future births), ~30% alive in 1930
        map.insert("demo_europe_age_neg40000_neg36000".to_string(), 0.020); // Moderate weighting on future births
        map.insert("demo_europe_age_neg36000_neg32000".to_string(), 0.018);
        map.insert("demo_europe_age_neg32000_neg28000".to_string(), 0.016);
        map.insert("demo_europe_age_neg28000_neg24000".to_string(), 0.015);
        map.insert("demo_europe_age_neg24000_neg20000".to_string(), 0.014);
        map.insert("demo_europe_age_neg20000_neg16000".to_string(), 0.013);
        map.insert("demo_europe_age_neg16000_neg12000".to_string(), 0.012);
        map.insert("demo_europe_age_neg12000_neg8000".to_string(), 0.011);
        map.insert("demo_europe_age_neg8000_neg4000".to_string(), 0.010);
        map.insert("demo_europe_age_neg4000_0".to_string(), 0.009);          // Future births tapering
        map.insert("demo_europe_age_0_4000".to_string(), 0.008);             // Larger portion alive in 1930
        map.insert("demo_europe_age_4000_8000".to_string(), 0.007);
        map.insert("demo_europe_age_8000_12000".to_string(), 0.006);
        map.insert("demo_europe_age_12000_16000".to_string(), 0.005);
        map.insert("demo_europe_age_16000_20000".to_string(), 0.004);
        map.insert("demo_europe_age_20000_24000".to_string(), 0.003);
        map.insert("demo_europe_age_24000_28000".to_string(), 0.002);
        map.insert("demo_europe_age_28000_32000".to_string(), 0.002);
        
        // North America demographic distribution
        // 9% of global population based on 1930-2035 timeframe
        // Significant growth: ~75% at negative ages (future births), ~25% alive in 1930
        map.insert("demo_north_america_age_neg40000_neg36000".to_string(), 0.012);
        map.insert("demo_north_america_age_neg36000_neg32000".to_string(), 0.011);
        map.insert("demo_north_america_age_neg32000_neg28000".to_string(), 0.010);
        map.insert("demo_north_america_age_neg28000_neg24000".to_string(), 0.009);
        map.insert("demo_north_america_age_neg24000_neg20000".to_string(), 0.008);
        map.insert("demo_north_america_age_neg20000_neg16000".to_string(), 0.007);
        map.insert("demo_north_america_age_neg16000_neg12000".to_string(), 0.006);
        map.insert("demo_north_america_age_neg12000_neg8000".to_string(), 0.005);
        map.insert("demo_north_america_age_neg8000_neg4000".to_string(), 0.004);
        map.insert("demo_north_america_age_neg4000_0".to_string(), 0.003);
        map.insert("demo_north_america_age_0_4000".to_string(), 0.003);
        map.insert("demo_north_america_age_4000_8000".to_string(), 0.003);
        map.insert("demo_north_america_age_8000_12000".to_string(), 0.002);
        map.insert("demo_north_america_age_12000_16000".to_string(), 0.002);
        map.insert("demo_north_america_age_16000_20000".to_string(), 0.002);
        map.insert("demo_north_america_age_20000_24000".to_string(), 0.002);
        map.insert("demo_north_america_age_24000_28000".to_string(), 0.001);
        map.insert("demo_north_america_age_28000_32000".to_string(), 0.001);
        
        // South America demographic distribution
        // 6% of global population based on 1930-2035 timeframe
        // Strong growth: ~80% at negative ages (future births), ~20% alive in 1930
        map.insert("demo_south_america_age_neg40000_neg36000".to_string(), 0.010);
        map.insert("demo_south_america_age_neg36000_neg32000".to_string(), 0.009);
        map.insert("demo_south_america_age_neg32000_neg28000".to_string(), 0.008);
        map.insert("demo_south_america_age_neg28000_neg24000".to_string(), 0.007);
        map.insert("demo_south_america_age_neg24000_neg20000".to_string(), 0.006);
        map.insert("demo_south_america_age_neg20000_neg16000".to_string(), 0.005);
        map.insert("demo_south_america_age_neg16000_neg12000".to_string(), 0.004);
        map.insert("demo_south_america_age_neg12000_neg8000".to_string(), 0.004);
        map.insert("demo_south_america_age_neg8000_neg4000".to_string(), 0.003);
        map.insert("demo_south_america_age_neg4000_0".to_string(), 0.002);
        map.insert("demo_south_america_age_0_4000".to_string(), 0.002);
        map.insert("demo_south_america_age_4000_8000".to_string(), 0.002);
        map.insert("demo_south_america_age_8000_12000".to_string(), 0.001);
        map.insert("demo_south_america_age_12000_16000".to_string(), 0.001);
        map.insert("demo_south_america_age_16000_20000".to_string(), 0.001);
        map.insert("demo_south_america_age_20000_24000".to_string(), 0.001);
        map.insert("demo_south_america_age_24000_28000".to_string(), 0.001);
        map.insert("demo_south_america_age_28000_32000".to_string(), 0.001);
        
        // Oceania demographic distribution
        // 3% of global population based on 1930-2035 timeframe
        // Moderate growth: ~65% at negative ages (future births), ~35% alive in 1930
        map.insert("demo_oceania_age_neg40000_neg36000".to_string(), 0.004);
        map.insert("demo_oceania_age_neg36000_neg32000".to_string(), 0.004);
        map.insert("demo_oceania_age_neg32000_neg28000".to_string(), 0.003);
        map.insert("demo_oceania_age_neg28000_neg24000".to_string(), 0.003);
        map.insert("demo_oceania_age_neg24000_neg20000".to_string(), 0.003);
        map.insert("demo_oceania_age_neg20000_neg16000".to_string(), 0.003);
        map.insert("demo_oceania_age_neg16000_neg12000".to_string(), 0.002);
        map.insert("demo_oceania_age_neg12000_neg8000".to_string(), 0.002);
        map.insert("demo_oceania_age_neg8000_neg4000".to_string(), 0.002);
        map.insert("demo_oceania_age_neg4000_0".to_string(), 0.001);
        map.insert("demo_oceania_age_0_4000".to_string(), 0.002);
        map.insert("demo_oceania_age_4000_8000".to_string(), 0.002);
        map.insert("demo_oceania_age_8000_12000".to_string(), 0.002);
        map.insert("demo_oceania_age_12000_16000".to_string(), 0.001);
        map.insert("demo_oceania_age_16000_20000".to_string(), 0.001);
        map.insert("demo_oceania_age_20000_24000".to_string(), 0.001);
        map.insert("demo_oceania_age_24000_28000".to_string(), 0.001);
        map.insert("demo_oceania_age_28000_32000".to_string(), 0.001);
        
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
/// Special handling for colistin's historical discontinuation period.
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

/// Time-aware drug availability that accounts for historical discontinuation periods.
/// Currently handles colistin's abandonment (1970-1995) and reintroduction.
pub fn get_drug_availability_time_aware(drug_name: &str, region: &str, region_living: Option<&str>, time_step: usize) -> f64 {
    // Calculate simulation year (assuming time_step 0 = year 1930, one step per day)
    let simulation_year = 1930.0 + (time_step as f64 / 365.0);
    
    // Special case for colistin's historical discontinuation
    if drug_name == "colistin" {
        // Colistin timeline:
        // 1952-1970: Active use (18 years)
        // 1970-1995: Largely abandoned due to toxicity (25 years)  
        // 1995+: Reintroduced as last resort for MDR infections
        
        if simulation_year < 1952.0 {
            return 0.0; // Not yet introduced
        } else if simulation_year >= 1952.0 && simulation_year < 1970.0 {
            // Active use period - full availability
            return get_drug_availability(drug_name, region, region_living);
        } else if simulation_year >= 1970.0 && simulation_year < 1995.0 {
            // Discontinuation period - very limited availability (research/compassionate use only)
            return get_drug_availability(drug_name, region, region_living) * 0.05; // 5% of normal availability
        } else {
            // Reintroduction period - available but as last resort
            return get_drug_availability(drug_name, region, region_living);
        }
    }
    
    // For all other drugs, use standard availability
    get_drug_availability(drug_name, region, region_living)
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
            // ESBL resistance affects penicillins + some cephalosporins (BL/BLI combinations overcome ESBL)
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Fluoroquinolone resistance (often ciprofloxacin + levofloxacin together)
            vec!["ciprofloxacin", "levofloxacin"],
            // Aminoglycoside resistance (often linked)
            vec!["gentamicin", "tobramycin"],
        ]);

        // Acinetobacter baumannii resistance patterns
        m.insert("acinetobacter baumannii", vec![
            // β-lactamase affects most β-lactams (BL/BLI combinations included)
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Carbapenemase affects carbapenems (including BL/BLI)
            vec!["meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"],
            // Fluoroquinolone resistance
            vec!["ciprofloxacin", "levofloxacin", "moxifloxacin"],
            // Aminoglycoside resistance
            vec!["gentamicin", "tobramycin", "amikacin"],
        ]);

        // Klebsiella pneumoniae resistance patterns  
        m.insert("klebsiella pneumoniae", vec![
            // ESBL resistance (BL/BLI combinations overcome ESBL)
            vec!["penicilling", "ampicillin", "amoxicillin", "cephalexin", "cefazolin", "cefuroxime", "ceftriaxone", "amoxicillin_clavulanate", "ampicillin_sulbactam", "piperacillin_tazobactam", "ticarcillin_clavulanate"],
            // Carbapenemase (KPC, NDM, etc.)
            vec!["meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"],
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





/*

// Drug introduction dates - 
// FOR CODE DEVELOPMENT DIVIDING TIMES TO INTRODUCTION BY 10

lazy_static! {
    pub static ref DRUG_INTRODUCTION_DATES: HashMap<&'static str, usize> = {
        let mut map = HashMap::new();
        

        // Sulfonamides (first antibiotics)
        map.insert("sulfanilamide", 255);   // 2555 // 1937 (simulation start, 7 years after 1930)

        // Beta-lactams (Penicillins)
        map.insert("penicilling", 355);     // 3555 // 1942 (12 years after 1930)
        map.insert("ampicillin", 1131);     // 11315 // 1961 (31 years after 1930)
        map.insert("amoxicillin", 1378);    // 1972 (42 years after 1930)
        map.insert("piperacillin", 1606);   // 1981 (51 years after 1930)
        map.insert("ticarcillin", 1460);    // 1977 (47 years after 1930)
        // Beta-lactam/beta-lactamase inhibitor combinations
        map.insert("amoxicillin_clavulanate", 1642); // 1985 (55 years after 1930)
        map.insert("ampicillin_sulbactam", 1825);    // 1990 (60 years after 1930)
        map.insert("piperacillin_tazobactam", 1971); // 1984 (54 years after 1930)
        map.insert("ticarcillin_clavulanate", 1825); // 1990 (60 years after 1930)
        map.insert("meropenem_vaborbactam", 3204);   // 2018 (88 years after 1930)
        map.insert("ceftazidime_avibactam", 2774);   // 2006 (76 years after 1930)

        // Cephalosporins
        map.insert("cephalexin", 1460);     // 1970 (40 years after 1930)
        map.insert("cefazolin", 1570);      // 1973 (43 years after 1930)
        map.insert("cefuroxime", 1752);     // 1978 (48 years after 1930)
        map.insert("ceftriaxone", 1971);    // 1984 (54 years after 1930)
        map.insert("ceftazidime", 2008);    // 1985 (55 years after 1930)
        map.insert("cefepime", 2419);       // 1996 (66 years after 1930)
        map.insert("ceftaroline", 2930);    // 2010 (80 years after 1930)

        // Carbapenems
        map.insert("meropenem", 2419);      // 1996 (66 years after 1930)
        map.insert("imipenem_c", 2008);     // 1985 (55 years after 1930)
        map.insert("ertapenem", 2592);      // 2001 (71 years after 1930)

        // Monobactams
        map.insert("aztreonam", 2044);      // 1986 (56 years after 1930)

        // Macrolides
        map.insert("erythromycin", 802);    // 1952 (22 years after 1930)
        map.insert("azithromycin", 2226);   // 1991 (61 years after 1930)
        map.insert("clarithromycin", 2189); // 1990 (60 years after 1930)

        // Lincosamides
        map.insert("clindamycin", 1387);    // 1968 (38 years after 1930)

        // Aminoglycosides
        map.insert("gentamicin", 1204);      // 1963 (33 years after 1930)
        map.insert("tobramycin", 1632);     // 1975 (45 years after 1930)
        map.insert("amikacin", 1669);       // 1976 (46 years after 1930)

        // Fluoroquinolones
        map.insert("ciprofloxacin", 2080);  // 1987 (57 years after 1930)
        map.insert("levofloxacin", 2419);   // 1996 (66 years after 1930)
        map.insert("moxifloxacin", 2529);   // 1999 (69 years after 1930)
        map.insert("ofloxacin", 2189);      // 1990 (60 years after 1930)

        // Tetracyclines
        map.insert("tetracycline", 657);    // 1948 (18 years after 1930)
        map.insert("doxyclycline", 1350);   // 1967 (37 years after 1930)
        map.insert("minocycline", 1496);    // 1971 (41 years after 1930)

        // Glycopeptides
        map.insert("vancomycin", 1021);      // 1958 (28 years after 1930)
        map.insert("teicoplanin", 2117);    // 1988 (58 years after 1930)

        // Oxazolidinones
        map.insert("linezolid", 2555);      // 2000 (70 years after 1930)
        map.insert("tedizolid", 3066);      // 2014 (84 years after 1930)

        // Folate antagonists
        map.insert("trim_sulf", 1387);      // 1968 (38 years after 1930) - trimethoprim-sulfamethoxazole

        // Other antibiotics
        map.insert("quinu_dalfo", 2529);    // 1999 (69 years after 1930) - quinupristin/dalfopristin
        map.insert("chlorampheni", 693);    // 1949 (19 years after 1930) - chloramphenicol
        map.insert("nitrofurantoin", 839);  // 1953 (23 years after 1930)
        map.insert("retapamulin", 2840);    // 2007 (77 years after 1930) - topical antibiotic
        map.insert("fusidic_a", 1168);       // 1962 (32 years after 1930) - fusidic acid
        map.insert("metronidazole", 1096);   // 1960 (30 years after 1930)
        map.insert("furazolidone", 912);    // 1955 (25 years after 1930)
        
        // Polymyxins  
        map.insert("colistin", 802);        // 1952 (22 years after 1930) - clinical introduction
        
        map
    };
}

*/




// Drug introduction dates (as time steps from start of 1930)
// Each time step = 1 day, so multiply years by 365
lazy_static! {
    pub static ref DRUG_INTRODUCTION_DATES: HashMap<&'static str, usize> = {
        let mut map = HashMap::new();
        
        // Sulfonamides (first antibiotics)
        map.insert("sulfanilamide", 2555);   // 2555 // 1937 (simulation start, 7 years after 1930)

        // Beta-lactams (Penicillins)
        map.insert("penicilling", 3555);     // 3555 // 1942 (12 years after 1930)
        map.insert("ampicillin", 11315);     // 11315 // 1961 (31 years after 1930)
        map.insert("amoxicillin", 13780);    // 1972 (42 years after 1930)
        map.insert("piperacillin", 16065);   // 1981 (51 years after 1930)
        map.insert("ticarcillin", 14600);    // 1977 (47 years after 1930)

        // Beta-lactam/beta-lactamase inhibitor combinations
        map.insert("amoxicillin_clavulanate", 16425); // 1985 (55 years after 1930)
        map.insert("ampicillin_sulbactam", 18250);    // 1990 (60 years after 1930)
        map.insert("piperacillin_tazobactam", 19715); // 1984 (54 years after 1930)
        map.insert("ticarcillin_clavulanate", 18250); // 1990 (60 years after 1930)
        map.insert("meropenem_vaborbactam", 32045);   // 2018 (88 years after 1930)
        map.insert("ceftazidime_avibactam", 27740);   // 2006 (76 years after 1930)

        // Cephalosporins
        map.insert("cephalexin", 14605);     // 1970 (40 years after 1930)
        map.insert("cefazolin", 15700);      // 1973 (43 years after 1930)
        map.insert("cefuroxime", 17525);     // 1978 (48 years after 1930)
        map.insert("ceftriaxone", 19715);    // 1984 (54 years after 1930)
        map.insert("ceftazidime", 20080);    // 1985 (55 years after 1930)
        map.insert("cefepime", 24195);       // 1996 (66 years after 1930)
        map.insert("ceftaroline", 29305);    // 2010 (80 years after 1930)

        // Carbapenems
        map.insert("meropenem", 24195);      // 1996 (66 years after 1930)
        map.insert("imipenem_c", 20080);     // 1985 (55 years after 1930)
        map.insert("ertapenem", 25920);      // 2001 (71 years after 1930)

        // Monobactams
        map.insert("aztreonam", 20445);      // 1986 (56 years after 1930)

        // Macrolides
        map.insert("erythromycin", 8025);    // 1952 (22 years after 1930)
        map.insert("azithromycin", 22260);   // 1991 (61 years after 1930)
        map.insert("clarithromycin", 21895); // 1990 (60 years after 1930)

        // Lincosamides
        map.insert("clindamycin", 13870);    // 1968 (38 years after 1930)

        // Aminoglycosides
        map.insert("gentamicin", 12045);      // 1963 (33 years after 1930)
        map.insert("tobramycin", 16325);     // 1975 (45 years after 1930)
        map.insert("amikacin", 16690);       // 1976 (46 years after 1930)

        // Fluoroquinolones
        map.insert("ciprofloxacin", 20805);  // 1987 (57 years after 1930)
        map.insert("levofloxacin", 24195);   // 1996 (66 years after 1930)
        map.insert("moxifloxacin", 25290);   // 1999 (69 years after 1930)
        map.insert("ofloxacin", 21895);      // 1990 (60 years after 1930)

        // Tetracyclines
        map.insert("tetracycline", 6575);    // 1948 (18 years after 1930)
        map.insert("doxyclycline", 13505);   // 1967 (37 years after 1930)
        map.insert("minocycline", 14965);    // 1971 (41 years after 1930)

        // Glycopeptides
        map.insert("vancomycin", 10215);      // 1958 (28 years after 1930)
        map.insert("teicoplanin", 21170);    // 1988 (58 years after 1930)

        // Oxazolidinones
        map.insert("linezolid", 25550);      // 2000 (70 years after 1930)
        map.insert("tedizolid", 30660);      // 2014 (84 years after 1930)

        // Folate antagonists
        map.insert("trim_sulf", 13870);      // 1968 (38 years after 1930) - trimethoprim-sulfamethoxazole

        // Other antibiotics
        map.insert("quinu_dalfo", 25290);    // 1999 (69 years after 1930) - quinupristin/dalfopristin
        map.insert("chlorampheni", 6935);    // 1949 (19 years after 1930) - chloramphenicol
        map.insert("nitrofurantoin", 8395);  // 1953 (23 years after 1930)
        map.insert("retapamulin", 28405);    // 2007 (77 years after 1930) - topical antibiotic
        map.insert("fusidic_a", 11680);       // 1962 (32 years after 1930) - fusidic acid
        map.insert("metronidazole", 10965);   // 1960 (30 years after 1930)
        map.insert("furazolidone", 9125);    // 1955 (25 years after 1930)
        
        // Polymyxins  
        map.insert("colistin", 8020);        // 1952 (22 years after 1930) - clinical introduction
        
        map
    };
}








/// Gets the introduction date for a drug (as time step from 1930)
/// Returns the time step if found, None otherwise
pub fn get_drug_introduction_time_step(drug_name: &str) -> Option<usize> {
    DRUG_INTRODUCTION_DATES.get(drug_name).copied()
}

/// Samples an age and region from the 120-parameter demographic distribution
/// Returns (region, age_in_years)
pub fn sample_age_and_region_from_distribution(rng: &mut impl rand::Rng) -> (crate::simulation::population::Region, i32) {
    use crate::simulation::population::Region;
    
    // Build cumulative probability distribution
    let mut cumulative_probs = Vec::new();
    let mut running_total = 0.0;
    
    // Define regions and age bands (-40000 to +32000 in 4000-year bands)
    let regions = [
        Region::Asia,
        Region::Africa, 
        Region::Europe,
        Region::NorthAmerica,
        Region::SouthAmerica,
        Region::Oceania,
    ];
    
    let age_bands = [
        (-40000, -36000), (-36000, -32000), (-32000, -28000), (-28000, -24000), (-24000, -20000),
        (-20000, -16000), (-16000, -12000), (-12000, -8000), (-8000, -4000), (-4000, 0),
        (0, 4000), (4000, 8000), (8000, 12000), (12000, 16000), (16000, 20000),
        (20000, 24000), (24000, 28000), (28000, 32000)
    ];
    
    // Build distribution
    for region in &regions {
        let region_name = match region {
            Region::Asia => "asia",
            Region::Africa => "africa",
            Region::Europe => "europe", 
            Region::NorthAmerica => "north_america",
            Region::SouthAmerica => "south_america",
            Region::Oceania => "oceania",
            Region::Home => "north_america", // Default fallback
        };
        
        for (age_min, age_max) in &age_bands {
            let param_name = if *age_min < 0 && *age_max <= 0 {
                format!("demo_{}_age_neg{}_neg{}", region_name, (*age_min as i32).abs(), (*age_max as i32).abs())
            } else if *age_min < 0 && *age_max > 0 {
                format!("demo_{}_age_neg{}_0", region_name, (*age_min as i32).abs())
            } else {
                format!("demo_{}_age_{}_{}", region_name, age_min, age_max)
            };
            let prob = get_global_param(&param_name).unwrap_or(0.0);
            running_total += prob;
            cumulative_probs.push((running_total, *region, *age_min, *age_max));
        }
    }
    
    // Sample from distribution
    let random_value = rng.gen::<f64>() * running_total;
    
    for (cumulative_prob, region, age_min, age_max) in cumulative_probs {
        if random_value <= cumulative_prob {
            // Sample a random age within the band
            let age = rng.gen_range(age_min..=age_max);
            return (region, age);
        }
    }
    
    // Fallback (should rarely be reached)
    (Region::Asia, 0)
}

