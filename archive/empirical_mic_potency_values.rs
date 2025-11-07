// Empirical MIC-based potency values for drug-bacteria combinations
// Format: potency = 1.0 / empirical_mic_susceptible (mg/L)
// This anchors our 0-10 drug level scale to real-world MIC breakpoints

use std::collections::HashMap;

pub fn add_empirical_mic_based_potencies(map: &mut HashMap<String, f64>) {
    
    // ============================================================================
    // ESCHERICHIA COLI - Most clinically important gram-negative
    // ============================================================================
    
    // Beta-lactams (penicillins)
    map.insert("drug_amoxicillin_for_bacteria_escherichia_coli_potency".to_string(), 0.125); // 1.0/8.0 mg/L
    map.insert("drug_ampicillin_for_bacteria_escherichia_coli_potency".to_string(), 0.125);  // 1.0/8.0 mg/L
    map.insert("drug_piperacillin_for_bacteria_escherichia_coli_potency".to_string(), 0.0625); // 1.0/16.0 mg/L
    
    // Beta-lactam combinations
    map.insert("drug_amoxicillin_clavulanate_for_bacteria_escherichia_coli_potency".to_string(), 0.25); // 1.0/4.0 mg/L
    map.insert("drug_piperacillin_tazobactam_for_bacteria_escherichia_coli_potency".to_string(), 0.25); // 1.0/4.0 mg/L
    
    // Cephalosporins
    map.insert("drug_cephalexin_for_bacteria_escherichia_coli_potency".to_string(), 0.125); // 1.0/8.0 mg/L
    map.insert("drug_cefazolin_for_bacteria_escherichia_coli_potency".to_string(), 0.5);    // 1.0/2.0 mg/L
    map.insert("drug_cefuroxime_for_bacteria_escherichia_coli_potency".to_string(), 0.25);  // 1.0/4.0 mg/L
    map.insert("drug_ceftriaxone_for_bacteria_escherichia_coli_potency".to_string(), 1.0);  // 1.0/1.0 mg/L
    map.insert("drug_ceftazidime_for_bacteria_escherichia_coli_potency".to_string(), 1.0);  // 1.0/1.0 mg/L
    map.insert("drug_cefepime_for_bacteria_escherichia_coli_potency".to_string(), 1.0);    // 1.0/1.0 mg/L
    
    // Carbapenems
    map.insert("drug_meropenem_for_bacteria_escherichia_coli_potency".to_string(), 4.0);    // 1.0/0.25 mg/L
    map.insert("drug_imipenem_c_for_bacteria_escherichia_coli_potency".to_string(), 2.0);   // 1.0/0.5 mg/L
    map.insert("drug_ertapenem_for_bacteria_escherichia_coli_potency".to_string(), 2.0);    // 1.0/0.5 mg/L
    
    // Fluoroquinolones
    map.insert("drug_ciprofloxacin_for_bacteria_escherichia_coli_potency".to_string(), 4.0);    // 1.0/0.25 mg/L
    map.insert("drug_levofloxacin_for_bacteria_escherichia_coli_potency".to_string(), 4.0);     // 1.0/0.25 mg/L
    map.insert("drug_moxifloxacin_for_bacteria_escherichia_coli_potency".to_string(), 2.0);     // 1.0/0.5 mg/L
    
    // Aminoglycosides
    map.insert("drug_gentamicin_for_bacteria_escherichia_coli_potency".to_string(), 0.5);   // 1.0/2.0 mg/L
    map.insert("drug_tobramycin_for_bacteria_escherichia_coli_potency".to_string(), 0.5);   // 1.0/2.0 mg/L
    map.insert("drug_amikacin_for_bacteria_escherichia_coli_potency".to_string(), 0.0625); // 1.0/16.0 mg/L
    
    // Other agents
    map.insert("drug_trim_sulf_for_bacteria_escherichia_coli_potency".to_string(), 2.0);    // 1.0/0.5 mg/L
    map.insert("drug_nitrofurantoin_for_bacteria_escherichia_coli_potency".to_string(), 0.03125); // 1.0/32.0 mg/L
    map.insert("drug_colistin_for_bacteria_escherichia_coli_potency".to_string(), 0.5);     // 1.0/2.0 mg/L

    // ============================================================================
    // STAPHYLOCOCCUS AUREUS - Most important gram-positive
    // ============================================================================
    
    // Beta-lactams (MSSA)
    map.insert("drug_penicilling_for_bacteria_staphylococcus_aureus_potency".to_string(), 8.0);     // 1.0/0.125 mg/L
    map.insert("drug_ampicillin_for_bacteria_staphylococcus_aureus_potency".to_string(), 1.0);      // 1.0/1.0 mg/L
    map.insert("drug_amoxicillin_for_bacteria_staphylococcus_aureus_potency".to_string(), 1.0);     // 1.0/1.0 mg/L
    
    // Cephalosporins
    map.insert("drug_cephalexin_for_bacteria_staphylococcus_aureus_potency".to_string(), 0.25);     // 1.0/4.0 mg/L
    map.insert("drug_cefazolin_for_bacteria_staphylococcus_aureus_potency".to_string(), 0.5);       // 1.0/2.0 mg/L
    map.insert("drug_ceftriaxone_for_bacteria_staphylococcus_aureus_potency".to_string(), 0.125);   // 1.0/8.0 mg/L
    
    // Anti-MRSA agents
    map.insert("drug_vancomycin_for_bacteria_staphylococcus_aureus_potency".to_string(), 0.5);      // 1.0/2.0 mg/L
    map.insert("drug_teicoplanin_for_bacteria_staphylococcus_aureus_potency".to_string(), 0.5);     // 1.0/2.0 mg/L
    map.insert("drug_linezolid_for_bacteria_staphylococcus_aureus_potency".to_string(), 0.25);      // 1.0/4.0 mg/L
    map.insert("drug_tedizolid_for_bacteria_staphylococcus_aureus_potency".to_string(), 2.0);       // 1.0/0.5 mg/L
    
    // Fluoroquinolones
    map.insert("drug_ciprofloxacin_for_bacteria_staphylococcus_aureus_potency".to_string(), 1.0);   // 1.0/1.0 mg/L
    map.insert("drug_levofloxacin_for_bacteria_staphylococcus_aureus_potency".to_string(), 2.0);    // 1.0/0.5 mg/L
    map.insert("drug_moxifloxacin_for_bacteria_staphylococcus_aureus_potency".to_string(), 4.0);    // 1.0/0.25 mg/L
    
    // Macrolides
    map.insert("drug_erythromycin_for_bacteria_staphylococcus_aureus_potency".to_string(), 2.0);    // 1.0/0.5 mg/L
    map.insert("drug_azithromycin_for_bacteria_staphylococcus_aureus_potency".to_string(), 0.5);    // 1.0/2.0 mg/L
    map.insert("drug_clarithromycin_for_bacteria_staphylococcus_aureus_potency".to_string(), 4.0);  // 1.0/0.25 mg/L
    
    // Other agents
    map.insert("drug_clindamycin_for_bacteria_staphylococcus_aureus_potency".to_string(), 4.0);     // 1.0/0.25 mg/L
    map.insert("drug_fusidic_a_for_bacteria_staphylococcus_aureus_potency".to_string(), 2.0);       // 1.0/0.5 mg/L
    map.insert("drug_retapamulin_for_bacteria_staphylococcus_aureus_potency".to_string(), 10.0);    // 1.0/0.1 mg/L

    // ============================================================================
    // KLEBSIELLA PNEUMONIAE - Important ESKAPE pathogen
    // ============================================================================
    
    // Beta-lactams
    map.insert("drug_ampicillin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.001);   // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_amoxicillin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.001);  // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_piperacillin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.0625); // 1.0/16.0 mg/L
    
    // Beta-lactam combinations
    map.insert("drug_amoxicillin_clavulanate_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.125); // 1.0/8.0 mg/L
    map.insert("drug_piperacillin_tazobactam_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.25);  // 1.0/4.0 mg/L
    
    // Cephalosporins
    map.insert("drug_cefazolin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 1.0);       // 1.0/1.0 mg/L
    map.insert("drug_cefuroxime_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.125);    // 1.0/8.0 mg/L
    map.insert("drug_ceftriaxone_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 1.0);     // 1.0/1.0 mg/L
    map.insert("drug_ceftazidime_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 1.0);     // 1.0/1.0 mg/L
    map.insert("drug_cefepime_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 1.0);       // 1.0/1.0 mg/L
    
    // Carbapenems
    map.insert("drug_meropenem_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 4.0);       // 1.0/0.25 mg/L
    map.insert("drug_imipenem_c_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 2.0);      // 1.0/0.5 mg/L
    map.insert("drug_ertapenem_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 2.0);       // 1.0/0.5 mg/L
    
    // Fluoroquinolones
    map.insert("drug_ciprofloxacin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 4.0);   // 1.0/0.25 mg/L
    map.insert("drug_levofloxacin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 4.0);    // 1.0/0.25 mg/L
    
    // Aminoglycosides
    map.insert("drug_gentamicin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.5);      // 1.0/2.0 mg/L
    map.insert("drug_tobramycin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.5);      // 1.0/2.0 mg/L
    map.insert("drug_amikacin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.0625);    // 1.0/16.0 mg/L
    
    // Other agents
    map.insert("drug_trim_sulf_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 2.0);       // 1.0/0.5 mg/L
    map.insert("drug_colistin_for_bacteria_klebsiella_pneumoniae_potency".to_string(), 0.5);        // 1.0/2.0 mg/L

    // ============================================================================
    // PSEUDOMONAS AERUGINOSA - Intrinsically resistant gram-negative
    // ============================================================================
    
    // Beta-lactams (limited spectrum)
    map.insert("drug_penicilling_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.001);  // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_ampicillin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.001);   // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_amoxicillin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.001);  // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_piperacillin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.0625); // 1.0/16.0 mg/L
    
    // Beta-lactam combinations
    map.insert("drug_piperacillin_tazobactam_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.0625); // 1.0/16.0 mg/L
    
    // Cephalosporins (limited spectrum)
    map.insert("drug_cephalexin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.001);   // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_cefazolin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.001);    // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_cefuroxime_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.001);   // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_ceftriaxone_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.001);  // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_ceftazidime_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.125);  // 1.0/8.0 mg/L
    map.insert("drug_cefepime_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.125);    // 1.0/8.0 mg/L
    
    // Carbapenems
    map.insert("drug_meropenem_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.5);      // 1.0/2.0 mg/L
    map.insert("drug_imipenem_c_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.25);    // 1.0/4.0 mg/L
    map.insert("drug_ertapenem_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.001);    // 1.0/1000 mg/L (intrinsic resistance)
    
    // Fluoroquinolones
    map.insert("drug_ciprofloxacin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 4.0);  // 1.0/0.25 mg/L
    map.insert("drug_levofloxacin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.5);   // 1.0/2.0 mg/L
    
    // Aminoglycosides
    map.insert("drug_gentamicin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.25);    // 1.0/4.0 mg/L
    map.insert("drug_tobramycin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.25);    // 1.0/4.0 mg/L
    map.insert("drug_amikacin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.0625);   // 1.0/16.0 mg/L
    
    // Other agents
    map.insert("drug_colistin_for_bacteria_pseudomonas_aeruginosa_potency".to_string(), 0.5);       // 1.0/2.0 mg/L

    // ============================================================================
    // STREPTOCOCCUS PNEUMONIAE - Important respiratory pathogen
    // ============================================================================
    
    // Beta-lactams
    map.insert("drug_penicilling_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 16.0); // 1.0/0.0625 mg/L
    map.insert("drug_ampicillin_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 4.0);   // 1.0/0.25 mg/L
    map.insert("drug_amoxicillin_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 4.0);  // 1.0/0.25 mg/L
    
    // Cephalosporins
    map.insert("drug_ceftriaxone_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 1.0);  // 1.0/1.0 mg/L
    map.insert("drug_cefepime_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 0.5);    // 1.0/2.0 mg/L
    
    // Macrolides
    map.insert("drug_erythromycin_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 4.0); // 1.0/0.25 mg/L
    map.insert("drug_azithromycin_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 2.0); // 1.0/0.5 mg/L
    map.insert("drug_clarithromycin_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 4.0); // 1.0/0.25 mg/L
    
    // Fluoroquinolones
    map.insert("drug_levofloxacin_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 0.5); // 1.0/2.0 mg/L
    map.insert("drug_moxifloxacin_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 8.0); // 1.0/0.125 mg/L
    
    // Other agents
    map.insert("drug_vancomycin_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 1.0);   // 1.0/1.0 mg/L
    map.insert("drug_linezolid_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 0.5);    // 1.0/2.0 mg/L
    map.insert("drug_trim_sulf_for_bacteria_streptococcus_pneumoniae_potency".to_string(), 2.0);    // 1.0/0.5 mg/L

    // ============================================================================
    // ACINETOBACTER BAUMANNII - Highly resistant nosocomial pathogen
    // ============================================================================
    
    // Beta-lactams (limited activity)
    map.insert("drug_ampicillin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.001);  // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_amoxicillin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.001); // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_piperacillin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.0625); // 1.0/16.0 mg/L
    
    // Beta-lactam combinations
    map.insert("drug_piperacillin_tazobactam_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.0625); // 1.0/16.0 mg/L
    
    // Carbapenems (when susceptible)
    map.insert("drug_meropenem_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.125);   // 1.0/8.0 mg/L
    map.insert("drug_imipenem_c_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.25);   // 1.0/4.0 mg/L
    
    // Fluoroquinolones
    map.insert("drug_ciprofloxacin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 1.0); // 1.0/1.0 mg/L
    map.insert("drug_levofloxacin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.5);  // 1.0/2.0 mg/L
    
    // Aminoglycosides
    map.insert("drug_gentamicin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.25);   // 1.0/4.0 mg/L
    map.insert("drug_tobramycin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.25);   // 1.0/4.0 mg/L
    map.insert("drug_amikacin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.0625);  // 1.0/16.0 mg/L
    
    // Last resort agents
    map.insert("drug_colistin_for_bacteria_acinetobacter_baumannii_potency".to_string(), 0.5);      // 1.0/2.0 mg/L

    // ============================================================================
    // ENTEROCOCCUS FAECALIS - Gram-positive with intrinsic resistance
    // ============================================================================
    
    // Beta-lactams
    map.insert("drug_penicilling_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.125);   // 1.0/8.0 mg/L
    map.insert("drug_ampicillin_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.5);      // 1.0/2.0 mg/L
    map.insert("drug_amoxicillin_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.5);     // 1.0/2.0 mg/L
    
    // Glycopeptides
    map.insert("drug_vancomycin_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.25);     // 1.0/4.0 mg/L
    map.insert("drug_teicoplanin_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.125);   // 1.0/8.0 mg/L
    
    // Oxazolidinones
    map.insert("drug_linezolid_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.5);       // 1.0/2.0 mg/L
    map.insert("drug_tedizolid_for_bacteria_enterococcus_faecalis_potency".to_string(), 2.0);       // 1.0/0.5 mg/L
    
    // Quinupristin/dalfopristin (limited activity against E. faecalis)
    map.insert("drug_quinu_dalfo_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.03125); // 1.0/32.0 mg/L
    
    // Intrinsic resistance to many agents
    map.insert("drug_cephalexin_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.001);    // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_cefazolin_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.001);     // 1.0/1000 mg/L (intrinsic resistance)
    map.insert("drug_trim_sulf_for_bacteria_enterococcus_faecalis_potency".to_string(), 0.001);     // 1.0/1000 mg/L (intrinsic resistance)

    // ============================================================================
    // ENTEROCOCCUS FAECIUM - More resistant than E. faecalis
    // ============================================================================
    
    // Beta-lactams (high-level resistance)
    map.insert("drug_penicilling_for_bacteria_enterococcus_faecium_potency".to_string(), 0.0625);   // 1.0/16.0 mg/L
    map.insert("drug_ampicillin_for_bacteria_enterococcus_faecium_potency".to_string(), 0.0625);    // 1.0/16.0 mg/L
    map.insert("drug_amoxicillin_for_bacteria_enterococcus_faecium_potency".to_string(), 0.0625);   // 1.0/16.0 mg/L
    
    // Glycopeptides (VRE concern)
    map.insert("drug_vancomycin_for_bacteria_enterococcus_faecium_potency".to_string(), 0.25);      // 1.0/4.0 mg/L (when susceptible)
    map.insert("drug_teicoplanin_for_bacteria_enterococcus_faecium_potency".to_string(), 0.125);    // 1.0/8.0 mg/L (when susceptible)
    
    // Oxazolidinones
    map.insert("drug_linezolid_for_bacteria_enterococcus_faecium_potency".to_string(), 0.5);        // 1.0/2.0 mg/L
    map.insert("drug_tedizolid_for_bacteria_enterococcus_faecium_potency".to_string(), 2.0);        // 1.0/0.5 mg/L
    
    // Quinupristin/dalfopristin (good activity)
    map.insert("drug_quinu_dalfo_for_bacteria_enterococcus_faecium_potency".to_string(), 1.0);      // 1.0/1.0 mg/L

    // ============================================================================
    // MDR MYCOBACTERIUM TUBERCULOSIS - Special case
    // ============================================================================
    
    // First-line anti-TB drugs
    map.insert("drug_rifampicin_for_bacteria_mdr_mycobacterium_tuberculosis_potency".to_string(), 1.0);  // 1.0/1.0 mg/L
    
    // Second-line and backup agents (limited activity by definition in MDR-TB)
    map.insert("drug_moxifloxacin_for_bacteria_mdr_mycobacterium_tuberculosis_potency".to_string(), 2.0); // 1.0/0.5 mg/L
    map.insert("drug_levofloxacin_for_bacteria_mdr_mycobacterium_tuberculosis_potency".to_string(), 1.0); // 1.0/1.0 mg/L
    map.insert("drug_amikacin_for_bacteria_mdr_mycobacterium_tuberculosis_potency".to_string(), 0.25);    // 1.0/4.0 mg/L
    
    // Most other drugs have minimal activity against TB
    map.insert("drug_vancomycin_for_bacteria_mdr_mycobacterium_tuberculosis_potency".to_string(), 0.001); // 1.0/1000 mg/L
    map.insert("drug_ampicillin_for_bacteria_mdr_mycobacterium_tuberculosis_potency".to_string(), 0.001); // 1.0/1000 mg/L

    // ============================================================================
    // PLACEHOLDER VALUES for remaining combinations
    // These should be updated with actual empirical data as it becomes available
    // ============================================================================
    
    // For combinations not yet specified, use conservative default based on:
    // - Intrinsic resistance: 0.001 (1/1000 mg/L)
    // - Moderate activity: 0.1 (1/10 mg/L) 
    // - Good activity: 1.0 (1/1 mg/L)
    // - Excellent activity: 10.0 (1/0.1 mg/L)
    
    // TODO: Complete remaining ~1400 combinations systematically
    // Priority order: Clinical importance × data availability
    
    println!("Loaded {} empirical MIC-based potency values", map.len());
}